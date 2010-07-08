/*
 *  read_parameters.c
 *  Spike
 *
 *  Created by Ben Evans on 28/11/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "read_parameters.h"

int read_parameters(PARAMS * params, char * paramfile)
{
	int count = 0;
	int l = 0;
	char * string, buff[BUFFER];
	FILE * fp;
	
	fp = myfopen(paramfile, "r");
	while ((string = fgets(buff, sizeof(buff), fp)) != NULL) 	/* Read next line */
		count += parse_string(params, buff);
	fclose(fp);
	
	/* Additional calculations */
	float transP = (params->transP_Train >= params->transP_Test ? params->transP_Train : params->transP_Test);
	params->TotalTime = params->nStimuli * params->nTransPS * transP;
	params->TotalMS = ceil(params->TotalTime * 1000);
	params->TotalTS = ceil(params->TotalTime/DT);
	params->TSperMS = ceil(1/(DT*1000));
	params->spkBuffer = ceil(params->TotalTime / params->refract);
	params->SigmaE = params->noiseScale * (params->ThreshE - params->VhyperE);
	params->SigmaI = params->noiseScale * (params->ThreshI - params->VhyperI);
	//params->inpSpkBuff = ceil(transP / params->refract);
	
	if (params->useFilteredImages)
	{	
		fp = myfopen(IPFILE, "r");	/* Read imgParams file */
		while ((string = fgets(buff, sizeof(buff), fp)) != NULL) 	/* Read next line */
			count += parse_string(params, buff);
		fclose(fp); 
		params->sInputs = params->nPhases*params->nScales*params->nOrients*params->nRows*params->nCols;
	}
	else
	{
		params->sInputs = params->nExcit; // Calculate input size automatically?
		if (floor(params->a * params->nExcit) < 1)
			params->nFiringNeurons = floor(params->nExcit/(params->nStimuli * params->nTransPS));
		else
			params->nFiringNeurons = floor(params->nExcit * params->a);		
	}
	
	// Change all myalloc to myrealloc?
	
	/* Parse vectors of neurons in each layer */
	if (params->LvExcit == 0) // Initialize to 0 in defaults.m
	{
		params->vExcit = myalloc(params->nLayers*sizeof(int));
		params->vExcit[0] = params->sInputs;
		for (l=1; l<params->nLayers; l++)
			params->vExcit[l] = params->nExcit;
		params->LvExcit = params->nLayers;
	}
	else if (params->LvExcit == (params->nLayers - 1))
	{
		params->vExcit = myrealloc(params->vExcit, params->nLayers*sizeof(int));
		memmove((params->vExcit)+1, params->vExcit, params->LvExcit*sizeof(int));
		params->vExcit[0] = params->sInputs;
		params->LvExcit++;
	}
	else if (params->LvExcit == params->nLayers) // Implies vExcit != NULL
	{
		if (params->sInputs != params->vExcit[0]) // Implement reduced sampling here
			params->vExcit[0] = params->sInputs;
	}
	else
		exit_error("read_parameters", "vExcit size mismatch!");
	
	// *** Calculations for arrays of Inhibitory neurons
	/*if (params->rInhib > EPS) // Pass ratio for each layer?
	{
		params->vInhib = myalloc(params->nLayers*sizeof(int));
		for (l=0; l<params->nLayers; l++)
			params->vInhib[l] = (l==0 && !params->inputInhib) ? 0 : round(params->rInhib*params->vExcit[l]);
		params->LvInhib = params->nLayers;
	}*/
	
	
	if (params->LvInhib == 0)
	{
		params->vInhib = myalloc(params->nLayers*sizeof(int));
		for (l=0; l<params->nLayers; l++)
			params->vInhib[l] = params->nInhib;
	}
	else if (params->LvInhib && params->LvInhib != params->nLayers)
		exit_error("read_parameters", "vInhib size mismatch!");
	
	if (params->rInhib > EPS)
	{
		for (l=0; l<params->nLayers; l++)
			params->vInhib[l] = round(params->rInhib*params->vExcit[l]);
		params->LvInhib = params->nLayers;
	}
	
	if (!params->inputInhib)
		params->vInhib[0] = 0;
	
	
	
	/*if (params->LvInhib == 0)
	{
		params->vInhib = myalloc(params->nLayers*sizeof(int));
		for (l=0; l<params->nLayers; l++)
			params->vInhib[l] = (l==0 && !params->inputInhib) ? 0 : params->nInhib;
		params->LvInhib = params->nLayers;
	}
	else if (params->LvInhib == params->nLayers)
	{
		if (!params->inputInhib)
			params->vInhib[0] = 0;
	}
	else
		exit_error("read_parameters", "vInhib size mismatch!");*/
	
	if (params->LpEfE>0 || params->LpElE>0 || params->LpEI>0 || params->LpIE>0 || params->LpII>0)
		params->probConnect = true;
	else
		params->probConnect = false;
	
	if (params->LpEfE == 0)
	{
		assert(params->nSynEfE <= params->nExcit);
		params->pCnxEfE = myalloc(params->nLayers*sizeof(float));
		params->pCnxEfE[0] = 0.0;
		for (l=1; l<params->nLayers; l++)
		{
			params->pCnxEfE[l] = ((float) params->nSynEfE)/params->vExcit[l-1];
			assert((params->pCnxEfE[l] >= 0.0) && (params->pCnxEfE[l] <= 1.0));
		}
		params->LpEfE = l;
	}
	
	if (params->LpElE == 0)
	{
		assert(params->nSynElE <= params->nExcit);
		params->pCnxElE = myalloc(params->nLayers*sizeof(float));
		for (l=0; l<params->nLayers; l++)
		{
			params->pCnxElE[l] = ((float) params->nSynElE)/params->vExcit[l];
			assert((params->pCnxElE[l] >= 0.0) && (params->pCnxElE[l] <= 1.0));
		}
		params->LpElE = l;
	}
	
	if (params->LpEI == 0)
	{
		assert(params->nSynEI <= params->nExcit);
		params->pCnxEI = myalloc(params->nLayers*sizeof(float));
		for (l=0; l<params->nLayers; l++)
		{
			params->pCnxEI[l] = (params->vInhib[l]>0) ? ((float) params->nSynEI)/params->vExcit[l] : 0.0;
			assert((params->pCnxEI[l] >= 0.0) && (params->pCnxEI[l] <= 1.0));
		}
		params->LpEI = l;
	}
	
	if (params->LpIE == 0)
	{
		assert(params->nSynIE <= params->nInhib);
		params->pCnxIE = myalloc(params->nLayers*sizeof(float));
		for (l=0; l<params->nLayers; l++)
		{
			params->pCnxIE[l] = (params->vInhib[l]>0) ? ((float) params->nSynIE)/params->vInhib[l] : 0.0;
			assert((params->pCnxIE[l] >= 0.0) && (params->pCnxIE[l] <= 1.0));
		}
		params->LpIE = l;
	}
	
	if (params->LpII == 0)
	{
		assert(params->nSynII <= params->nInhib);
		params->pCnxII = myalloc(params->nLayers*sizeof(float));
		for (l=0; l<params->nLayers; l++)
		{
			params->pCnxII[l] = (params->vInhib[l]>0) ? ((float) params->nSynII)/params->vInhib[l] : 0.0;
			assert((params->pCnxII[l] >= 0.0) && (params->pCnxII[l] <= 1.0));
		}
		params->LpII = l;
	}
	
	assert(params->LpEfE == params->nLayers); // nWLayers
	assert(params->LpElE == params->nLayers);
	assert(params->LpEI == params->nLayers);
	assert(params->LpIE == params->nLayers);
	assert(params->LpII == params->nLayers);
	return count;
}

char * trim(char * string)
{ // Remove leading and trailing whitespace
	/* Initialize start, end pointers */
	if (string == NULL)
		return string;
	
	char *s1 = string, *s2 = &string[strlen(string) - 1];
	
	/* Trim and delimit right side */
	while ( (isspace(*s2) || *s2 == ';' || *s2 == ',' || *s2 == '\'') && (s2 >= s1) )
		s2--;
	*(s2+1) = '\0';
	
	/* Trim left side */
	while ( (isspace(*s1) || *s1 == ',' || *s1 == '`' || *s1 == '\'') && (s1 < s2) )
		s1++;
	
	/* Copy finished string */
	strcpy(string, s1);
	return string;
}

int parseIntVector(char * string, int ** array)
{ // Function to parse a row vector of integers, store it in an array and return the number of elements
	char * tstr = NULL;
	const char * delims = "[, ]";
	int numElements = 0;
	int block = 1;
	
	myfree(*array); // Call free if *array!=NULL in case the array already exists
	*array = myalloc(VECBUFF*sizeof(int));
	tstr = strtok(string, delims); // Get first token starting with the first nondelimiter char
	
	while ( (tstr != NULL) && (*trim(tstr) != ']') && (*tstr != '\0') ) // Trim and check not end
	{
		if (numElements == block*VECBUFF) // Check array is not full up
			*array = myrealloc(*array, ++block*VECBUFF*sizeof(int)); // Reallocate memory
		
		(*array)[numElements++] = atoi(tstr); // convert to int for assignment
		tstr = strtok(NULL, delims); // Get the next token
	}
	
	if (numElements == 0)
		fprintf(stderr, "Error! Vector: \"%s\" is empty",tstr);
	else
		*array = myrealloc(*array, numElements*sizeof(int));
	
	return numElements;
}

int parseFloatVector(char * string, float ** array)
{ // Function to parse a row vector of integers, store it in an array and return the number of elements
	char * tstr = NULL;
	const char * delims = "[, ]";
	int numElements = 0;
	int block = 1;
	
	myfree(*array); // Call free if *array!=NULL in case the array already exists
	*array = myalloc(VECBUFF*sizeof(float));
	tstr = strtok(string, delims); // Get first token starting with the first nondelimiter char
	
	while ( (tstr != NULL) && (*trim(tstr) != ']') && (*tstr != '\0') ) // Trim and check not end
	{
		if (numElements == block*VECBUFF) // Check array is not full up
			*array = myrealloc(*array, ++block*VECBUFF*sizeof(float)); // Reallocate memory
			
		(*array)[numElements++] = atof(tstr); // convert to int for assignment
		tstr = strtok(NULL, delims); // Get the next token
	}
	
	if (numElements == 0)
		fprintf(stderr, "Error! Vector: \"%s\" is empty",tstr);
	else
		*array = myrealloc(*array, numElements*sizeof(float));

	return numElements;
}
/*
struct data {
	int count;
	int * ints;
	float * floats;
}

typedef enum dataType {
	INTEGERS,
	FLOATS
}

int parseVector(char * string, data * vectors, dataType type)
{ // Function to parse a row vector of integers, store it in an array and return the number of elements
	char * tstr = NULL;
	const char * delims = "[, ]";
	int numElements = 0;
	int block = 1;
	
	free(array); // In case the array already exists
	*array = myalloc(VECBUFF*sizeof(int *));
	tstr = strtok(string, delims); // Get first token starting with the first nondelimiter char
	
	while ( (*trim(tstr) != ']') && (*tstr != '\0') && (tstr != NULL) ) // Trim and check not end
	{
		if (numElements == block*VECBUFF) // Check array is not full up
			*array = myrealloc(*array, ++block*VECBUFF*sizeof(int *)) // Reallocate memory
			
			*array[numElements++] = atoi(tstr); // convert to int for assignment
		tstr = strtok(NULL, delims); // Get the next token
	}
	
	if (numElements == 0)
		fprintf(stderr, "Error! Vector: \"%s\" is empty",tstr);
	else
		if ( (*array = realloc(*array, numElements*sizeof(int *))) == NULL ) // Set array to correct size
			exit_error("REALLOC", "NULL pointer from realloc");
	
	return numElements;
}*/

int parse_string(PARAMS * params, char * string)
{
	static char name[MAXLEN], value[MAXLEN];
	char * sptr;
	
	if (string[0] == '\n' || string[0] == '#' || string[0] == '%') 	/* Skip blank lines and comments */
		return 0; //continue;
	/* Parse name/value pair from line */
	sptr = strtok(string, "=");
	if (sptr==NULL)
	{
		fprintf(stderr, "Error! String: \"%s\" has no name!\n", sptr);
		return 0; // No parameter name on this line
	}
	else
		strncpy(name, sptr, MAXLEN);
	sptr = strtok(NULL, "=");
	if (sptr==NULL)
	{
		fprintf(stderr, "Error! Name: \"%s\" has no associated value!\n", name);
		return 0;
	}
	else
		strncpy(value, sptr, MAXLEN);
	trim(name); // Necessary if using indented variables
	trim(value);
	
	/* Copy into correct entry in parameters struct */ // Use strcmpi() for case insensitive comparisons
	/* Simulation */		
	if (strcmp(name, "loops")==0)
		params->loops = atoi(value);
	else if (strcmp(name, "train")==0)
		params->train = atoi(value);
	else if (strcmp(name, "pretrain")==0)
		params->pretrain = atoi(value);
	else if (strcmp(name, "trainPause")==0)
		params->trainPause = atoi(value);
	else if (strcmp(name, "noise")==0)
		params->noise = atoi(value);
	else if (strcmp(name, "noiseScale")==0)
		params->noiseScale = atof(value);
	else if (strcmp(name, "nRecordsPL")==0)
		params->nRecordsPL = atoi(value);
	else if (strcmp(name, "printConnections")==0)
		params->printConnections = atoi(value);
	else if (strcmp(name, "probConnect")==0)
		params->probConnect = atoi(value);
	
	/* Stimuli */
	else if (strcmp(name, "imageList")==0)
	{
		params->imageList = myalloc((strlen(trim(value))+1)*sizeof(char));
		strcpy(params->imageList, value);
		params->imageList[strlen(value)] = '\0';
	}
	else if (strcmp(name, "randStimOrder")==0)
		params->randStimOrder = atoi(value);
	else if (strcmp(name, "randTransOrder")==0)
		params->randTransOrder = atoi(value);
	else if (strcmp(name, "interleaveTrans")==0)
		params->interleaveTrans = atoi(value);
	else if (strcmp(name, "localRep")==0)
		params->localRep = atoi(value);
	else if (strcmp(name, "current")==0)
		params->current = atof(value);
	else if (strcmp(name, "currentSpread")==0)
		params->currentSpread = atof(value);
	else if (strcmp(name, "nStimuli")==0)
		params->nStimuli = atoi(value);
	else if (strcmp(name, "nTransPS")==0)
		params->nTransPS = atoi(value);
	else if (strcmp(name, "newTestSet")==0)
		params->newTestSet = atoi(value);
	else if (strcmp(name, "nTestStimuli")==0)
		params->nTestStimuli = atoi(value);
	else if (strcmp(name, "nTestTransPS")==0)
		params->nTestTransPS = atoi(value);
	else if (strcmp(name, "transP_Train")==0)
		params->transP_Train = atof(value);
	else if (strcmp(name, "transP_Test")==0)
		params->transP_Test = atof(value);
	else if (strcmp(name, "shift")==0)
		params->shift = atoi(value);
	else if (strcmp(name, "a")==0)
		params->a = atof(value);
	else if (strcmp(name, "useFilteredImages")==0)
		params->useFilteredImages = atoi(value);
	else if (strcmp(name, "gabor")==0)
		params->gabor = atoi(value);
	else if (strcmp(name, "vScales")==0)
	{
		/*if (params->nScales)
		{
			//params->vScales = myfree(params->vScales);
			params->nScales = 0;
		}*/
		params->nScales = parseIntVector(value, &params->vScales);
	}
	else if (strcmp(name, "vOrients")==0)
	{
		/*if (params->nOrients)
		{
			params->vOrients = myfree(params->vOrients);
			params->nOrients = 0;
		}*/
		params->nOrients = parseIntVector(value, &params->vOrients);
	}
	else if (strcmp(name, "vPhases")==0)
	{
		/*if (params->nPhases)
		{
			params->vPhases = myfree(params->vPhases);
			params->nPhases = 0;
		}*/
		params->nPhases = parseIntVector(value, &params->vPhases);
	}
	else if (strcmp(name, "nRows")==0)
		params->nRows = atoi(value);
	else if (strcmp(name, "nCols")==0)
		params->nCols = atoi(value);
	
	/* Network */
	else if (strcmp(name, "nLayers")==0)
	{
		params->nLayers = atoi(value);
		params->nWLayers = params->nLayers - 1;
	}
	else if (strcmp(name, "inputInhib")==0)
		params->inputInhib = atoi(value);
	else if (strcmp(name, "nExcit")==0)
		params->nExcit = atoi(value);
	else if (strcmp(name, "vExcit")==0)
	{
		/*if (params->LvExcit)
		{
			params->vExcit = myfree(params->vExcit);
			params->LvExcit = 0;
		}*/
		params->LvExcit = parseIntVector(value, &params->vExcit);
	}
	else if (strcmp(name, "nSynEfE")==0)
		params->nSynEfE = atoi(value);
	else if (strcmp(name, "pCnxEfE")==0)
	{
		/*if (params->LpEfE)
		{
			params->pCnxEfE = myfree(params->pCnxEfE);
			params->LpEfE = 0;
		}*/
		params->LpEfE = parseFloatVector(value, &params->pCnxEfE);
	}
	else if (strcmp(name, "nSynElE")==0)
		params->nSynElE = atoi(value);
	else if (strcmp(name, "pCnxElE")==0)
	{
		/*if (params->LpElE)
		{
			params->pCnxElE = myfree(params->pCnxElE);
			params->LpElE = 0;
		}*/
		params->LpElE = parseFloatVector(value, &params->pCnxElE);
	}
	else if (strcmp(name, "nSynIE")==0)
		params->nSynIE = atoi(value);
	else if (strcmp(name, "pCnxIE")==0)
	{
		/*if (params->LpIE)
		{
			params->pCnxIE = myfree(params->pCnxIE);
			params->LpIE = 0;
		}*/
		params->LpIE = parseFloatVector(value, &params->pCnxIE);
	}
	else if (strcmp(name, "nInhib")==0)
		params->nInhib = atoi(value);
	else if (strcmp(name, "rInhib")==0)
		params->rInhib = atof(value);
	else if (strcmp(name, "vInhib")==0)
	{
		/*if (params->LvInhib)
		{
			params->vInhib = myfree(params->vInhib);
			params->LvInhib = 0;
		}*/
		params->LvInhib = parseIntVector(value, &params->vInhib);
	}
	else if (strcmp(name, "nSynEI")==0)
		params->nSynEI = atoi(value);
	else if (strcmp(name, "pCnxEI")==0)
	{
		/*if (params->LpEI)
		{
			params->pCnxEI = myfree(params->pCnxEI);
			params->LpEI = 0;
		}*/
		params->LpEI = parseFloatVector(value, &params->pCnxEI);
	}
	else if (strcmp(name, "nSynII")==0)
		params->nSynII = atoi(value);
	else if (strcmp(name, "pCnxII")==0)
	{
		/*if (params->LpII)
		{
			params->pCnxII = myfree(params->pCnxII);
			params->LpII = 0;
		}*/
		params->LpII = parseFloatVector(value, &params->pCnxII);
	}
	else if (strcmp(name, "noise")==0)
		params->noise = atoi(value);
	else if (strcmp(name, "axonDelay")==0)
		params->axonDelay = atoi(value);
	else if (strcmp(name, "d_const")==0)
		params->d_const = atof(value);
	else if (strcmp(name, "d_min")==0)
		params->d_min = atof(value);
	else if (strcmp(name, "d_max")==0)
		params->d_max = atof(value);
	else if (strcmp(name, "d_mean")==0)
		params->d_mean = atof(value);
	else if (strcmp(name, "d_sd")==0)
		params->d_sd = atof(value);
	
	/* Cell bodies */
	else if (strcmp(name, "capE")==0)
		params->capE = atof(value);
	else if (strcmp(name, "capI")==0)
		params->capI = atof(value);
	else if (strcmp(name, "gLeakE")==0)
		params->gLeakE = atof(value);
	else if (strcmp(name, "gLeakI")==0)
		params->gLeakI = atof(value);
	else if (strcmp(name, "VrestE")==0)
		params->VrestE = atof(value);
	else if (strcmp(name, "VrestI")==0)
		params->VrestI = atof(value);
	else if (strcmp(name, "VhyperE")==0)
		params->VhyperE = atof(value);
	else if (strcmp(name, "VhyperI")==0)
		params->VhyperI = atof(value);
	else if (strcmp(name, "ThreshE")==0)
		params->ThreshE = atof(value);
	else if (strcmp(name, "ThreshI")==0)
		params->ThreshI = atof(value);
	else if (strcmp(name, "VrevE")==0)
		params->VrevE = atof(value);
	else if (strcmp(name, "VrevI")==0)
		params->VrevI = atof(value);
	else if (strcmp(name, "refract")==0)
		params->refract = atof(value);
	
	/* Synapses */
	else if (strcmp(name, "alphaC")==0)
		params->alphaC = atof(value);
	else if (strcmp(name, "tauC")==0)
		params->tauC = atof(value);
	else if (strcmp(name, "alphaD")==0)
		params->alphaD = atof(value);
	else if (strcmp(name, "tauD")==0)
		params->tauD = atof(value);
	else if (strcmp(name, "learnR")==0)
		params->learnR = atof(value);
	else if (strcmp(name, "tauEE")==0)
		params->tauEE = atof(value);
	else if (strcmp(name, "Dg_IE")==0)
		params->Dg_IE = atof(value);
	else if (strcmp(name, "tauIE")==0)
		params->tauIE = atof(value);
	else if (strcmp(name, "Dg_EI")==0)
		params->Dg_EI = atof(value);
	else if (strcmp(name, "tauEI")==0)
		params->tauEI = atof(value);
	else if (strcmp(name, "Dg_II")==0)
		params->Dg_II = atof(value);
	else if (strcmp(name, "tauII")==0)
		params->tauII = atof(value);
	else if (strcmp(name, "gMax")==0)
		params->gMax = atof(value);
	
	else
	{
		fprintf(stderr, "WARNING: %s/%s: Unknown name/value pair!\n", name, value);
		exit (-1); //return -1;
	}
	
	return 1; // Parameter succesfully loaded
}

int printParameters(PARAMS * mp, char * paramfile)
{
	char * rfpath = NULL;
	FILE * pFile = NULL;
	FILE * parameters_ptr, * params_FP;
	int c=0; // count
	int l=0;
	
	if (paramfile != NULL)
		sprintf(rfpath, "%s/%s/", RESDIR, paramfile);// strtok(pfile, "."));
	else
		sprintf(rfpath, "%s/", RESDIR);
	
	parameters_ptr = myfopen(MPFILE, "w"); // Variables to read into Matlab
	fprintf(parameters_ptr, "DT=%f;\n",DT); c++;
	fprintf(parameters_ptr, "TotalMS=%d;\n",mp->TotalMS); c++;
	fprintf(parameters_ptr, "rfpath='%s';\n\n",rfpath); c++;
	if (paramfile != NULL)
	{
		params_FP = myfopen(paramfile, "r");
		filecopy(params_FP, parameters_ptr);
		fclose(params_FP);
	}
	
	/*if (p_flag)
	{
		cli_FP = myfopen(cpfile, "r");
		fprintf(parameters_ptr, "\n");
		filecopy(cli_FP, parameters_ptr);
		fclose(cli_FP);
	}*/
	fclose(parameters_ptr);
	
	//fprintf(pFile, "MP.loops = %d;\n", mp->loops);
	FPRINT_INT(pFile, mp->loops); c++;
	FPRINT_INT(pFile, mp->train); c++;
	FPRINT_INT(pFile, mp->pretrain); c++;
	FPRINT_INT(pFile, mp->trainPause); c++;
	FPRINT_INT(pFile, mp->noise); c++;
	FPRINT_INT(pFile, mp->nRecordsPL); c++;
	FPRINT_INT(pFile, mp->probConnect); c++;
	FPRINT_STRING(pFile, mp->imageList); c++;
	FPRINT_INT(pFile, mp->randStimOrder); c++;
	FPRINT_INT(pFile, mp->randTransOrder); c++;
	FPRINT_INT(pFile, mp->interleaveTrans); c++;
	FPRINT_INT(pFile, mp->localRep); c++;
	FPRINT_FLOAT(pFile, mp->current); c++;
	FPRINT_FLOAT(pFile, mp->currentSpread); c++;
	//FPRINT_INT(pFile, mp->nStimuli);
	//FPRINT_INT(pFile, mp->nTransPS);
	FPRINT_INT(pFile, mp->newTestSet); c++;
	if (mp->newTestSet)
	{
		FPRINT_INT(pFile, mp->nTestStimuli); c++;
		FPRINT_INT(pFile, mp->nTestTransPS); c++;
	}
	FPRINT_FLOAT(pFile, mp->transP_Train); c++;
	FPRINT_FLOAT(pFile, mp->transP_Test); c++;
	FPRINT_INT(pFile, mp->useFilteredImages); c++;
	if (mp->useFilteredImages)
	{
		FPRINT_INT(pFile, mp->gabor); c++;
		// load other parameters from imageParams.m
	}
	else
		FPRINT_INT(pFile, mp->shift); c++;
	
	FPRINT_INT(pFile, mp->nLayers); c++;
	FPRINT_INT(pFile, mp->nWLayers); c++;
	FPRINT_INT(pFile, mp->inputInhib); c++;
	fprintf(pFile, "vExcit = [%d",mp->vExcit[0]); c++;
	for (l=1; l<mp->nLayers; l++)
		fprintf(pFile, ",%d", mp->vExcit[l]); c++;
	fprintf(pFile, "];\n");
	fprintf(pFile, "vInhib = [%d",mp->vInhib[0]); c++;
	for (l=1; l<mp->nLayers; l++)
		fprintf(pFile, ",%d", mp->vInhib[l]); c++;
	fprintf(pFile, "];\n");
	FPRINT_FLOAT(pFile, mp->rInhib); c++;
	fprintf(pFile, "pCnxEfE = [%G",mp->pCnxEfE[0]); c++;
	for (l=1; l<mp->nLayers; l++)
		fprintf(pFile, ",%G", mp->pCnxEfE[l]); c++;
	fprintf(pFile, "];\n");
	fprintf(pFile, "pCnxElE = [%G",mp->pCnxElE[0]); c++;
	for (l=1; l<mp->nLayers; l++)
		fprintf(pFile, ",%G", mp->pCnxElE[l]); c++;
	fprintf(pFile, "];\n");
	printFloatArray(pFile, "pCnxIE", mp->pCnxIE, mp->nLayers); c++; // nSyn?
	printFloatArray(pFile, "pCnxEI", mp->pCnxEI, mp->nLayers); c++;
	printFloatArray(pFile, "pCnxII", mp->pCnxII, mp->nLayers); c++;
	switch (mp->axonDelay)
	{
		case MinD:
			fprintf(pFile, "axonDelay = MinD;\n"); c++;
			FPRINT_FLOAT(pFile, mp->d_const); c++;
			break;
		case ConstD:
			fprintf(pFile, "axonDelay = ConstD;\n"); c++;
			FPRINT_FLOAT(pFile, mp->d_const); c++;
			break;
		case UniformD:
			fprintf(pFile, "axonDelay = UniformD;\n"); c++;
			FPRINT_FLOAT(pFile, mp->d_min); c++;
			FPRINT_FLOAT(pFile, mp->d_max); c++;
			break;
		case GaussD:
			fprintf(pFile, "axonDelay = GaussD;\n"); c++;
			FPRINT_FLOAT(pFile, mp->d_mean); c++;
			FPRINT_FLOAT(pFile, mp->d_sd); c++;
			break;
		default:
			break;
	}
	
	/* Cell bodies */
	FPRINT_FLOAT(pFile, mp->capE); c++;
	FPRINT_FLOAT(pFile, mp->capI); c++;
	FPRINT_FLOAT(pFile, mp->gLeakE); c++;
	FPRINT_FLOAT(pFile, mp->gLeakI); c++;
	FPRINT_FLOAT(pFile, mp->VrestE); c++;
	FPRINT_FLOAT(pFile, mp->VrestI); c++;
	FPRINT_FLOAT(pFile, mp->VhyperE); c++;
	FPRINT_FLOAT(pFile, mp->VhyperI); c++;
	FPRINT_FLOAT(pFile, mp->ThreshE); c++;
	FPRINT_FLOAT(pFile, mp->ThreshI); c++;
	FPRINT_FLOAT(pFile, mp->VrevE); c++;
	FPRINT_FLOAT(pFile, mp->VrevI); c++;
	FPRINT_FLOAT(pFile, mp->refract); c++;
	
	/* Synapses */
	FPRINT_FLOAT(pFile, mp->alphaC); c++;
	FPRINT_FLOAT(pFile, mp->tauC); c++;
	FPRINT_FLOAT(pFile, mp->alphaD); c++;
	FPRINT_FLOAT(pFile, mp->tauD); c++;
	FPRINT_FLOAT(pFile, mp->learnR); c++;
	FPRINT_FLOAT(pFile, mp->tauEE); c++;
	FPRINT_FLOAT(pFile, mp->Dg_IE); c++;
	FPRINT_FLOAT(pFile, mp->tauIE); c++;
	FPRINT_FLOAT(pFile, mp->Dg_EI); c++;
	FPRINT_FLOAT(pFile, mp->tauEI); c++;
	FPRINT_FLOAT(pFile, mp->Dg_II); c++;
	FPRINT_FLOAT(pFile, mp->tauII); c++;
	FPRINT_FLOAT(pFile, mp->gMax); c++;

	return c;
}

void printIntArray(FILE * fp, char * name, int * array, int len)
{
	int l=0;
	fprintf(fp, "%s = [%d",name,array[0]);
	for (l=1; l<len; l++)
		fprintf(fp, ",%d", array[l]);
	fprintf(fp, "];\n");
	fflush(fp);	
}

void printFloatArray(FILE * fp, char * name, float * array, int len)
{
	int l=0;
	fprintf(fp, "%s = [%G",name,array[0]);
	for (l=1; l<len; l++)
		fprintf(fp, ",%G", array[l]);
	fprintf(fp, "];\n");
	fflush(fp);	
	
}
