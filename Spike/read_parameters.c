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
	//int square = 0;
	char * string, buff[BUFFER];
	FILE * fp;
	
	fp = myfopen(paramfile, "r");
	while ((string = fgets(buff, sizeof(buff), fp)) != NULL) 	/* Read next line */
		count += parse_string(params, buff);
	fclose(fp);
	
	/* Additional calculations */
	
	if (mp->K > 1)		// Assumption: when training with multiples, test with individual stimuli
	{
		mp->newTestSet = true;
		mp->nTestStimuli = mp->nStimuli;
		int nCombs = gsl_sf_choose(mp->nTestStimuli-1, mp->K-1);
		if (mp->M > nCombs)
			mp->M = nCombs;
		mp->M = (mp->nStimuli == 1) ? 1 : mp->M;
		mp->nStimuli = mp->nTestStimuli * mp->M;
		/*nCombs = gsl_sf_choose(mp->nTestStimuli, mp->K);
		 mp->nStimuli = (nCombs == 1) ? nCombs : mp->nStimuli;*/
		
		//mp->nStimuli = mp->M * gsl_sf_choose(mp->nTestStimuli, mp->K) / mp->nTestStimuli; // Check ////////
		
		mp->nTestTransPS = mp->nTransPS;
	}
	else
	{
		mp->newTestSet = false;
		mp->nTestStimuli = mp->nStimuli;
		mp->nTestTransPS = mp->nTransPS;
	}
	
	float transP = (params->transP_Train >= params->transP_Test ? params->transP_Train : params->transP_Test);
	params->TotalTime = params->nStimuli * params->nTransPS * transP;
	params->TotalMS = round(params->TotalTime * 1000); // Was ceil, used for recording structures
	params->TotalTS = round(params->TotalTime/mp->DT); // Was ceil
	params->TSperMS = round(1/(mp->DT*1000)); // Was ceil, used for recording structures
	params->spkBuffer = ceil(params->TotalTime / params->refract);
	params->SigmaE = params->noiseScale * (params->ThreshE - params->VhyperE);
	params->SigmaI = params->noiseScale * (params->ThreshI - params->VhyperI);
	//params->inpSpkBuff = ceil(transP / params->refract);
	
	if (params->useFilteredImages && !(strcmp(paramfile, MPFILE)==0)) // Do not parse IPFILE if rerunning with parameters.m
	{
		int slen = strlen(mp->imgDir)+1+strlen(IPFILE)+1; // '/' & '\0'
		char * imgParams = myalloc(slen);
		if (snprintf(imgParams, slen, "%s/%s", mp->imgDir, IPFILE) >= slen)
			fprintf(stderr, "Warning! Undersized buffer: %s", imgParams);
		/*strncpy(imgParams, mp->imgDir, strlen(mp->imgDir)+1);
		strncat(imgParams, "/", 1);
		strncat(imgParams, IPFILE, strlen(IPFILE)+1);
		imgParams[slen-1]='\0';*/
		fp = myfopen(imgParams, "r");	/* Read imgParams file */
		while ((string = fgets(buff, sizeof(buff), fp)) != NULL) 	/* Read next line */
			count += parse_string(params, buff);
		fclose(fp); 
		params->sInputs = params->nPhases*params->nScales*params->nOrients*params->nRows*params->nCols;
	}
	else
	{
		params->sInputs = (params->LvExcit) ? params->vExcit[0] : params->nExcit;
		if (floor(params->a * params->sInputs) < 1) // Change to use and output nFiringNeurons instead
		{
			params->nFiringNeurons = floor(params->sInputs/(params->nStimuli * params->nTransPS));
			params->a = params->nFiringNeurons / params->sInputs;
		}
		else
			params->nFiringNeurons = floor(params->sInputs * params->a);	
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
		params->vExcit = myrealloc(params->vExcit, params->nLayers*sizeof(params->vExcit[0]));
		memmove((params->vExcit)+1, params->vExcit, params->LvExcit*sizeof(params->vExcit[0]));
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
	
	//if (mp->SOM) // && !mp->useFilteredImages)
	if (!mp->layDim)
		mp->layDim = myalloc(mp->nLayers * sizeof(DIM)); //get_2D_iarray(mp->nLayers, 2, 0);
	else
	{
		mp->layDim = myrealloc(mp->layDim, mp->nLayers * sizeof(DIM));
		mp->vSquare = myrealloc(mp->vSquare, mp->nLayers * sizeof(mp->vSquare[0])); // Should just change to bool sqInput
		for (l=2; l<mp->nLayers; l++)
			mp->vSquare[l] = mp->vSquare[1]; // Assume that if l==1 is square, then all l>1 will be
	}
		
	
	for (l=0; l<mp->nLayers; l++)
	{
		if (mp->useFilteredImages && l==0) // 3D (input) layer
		{
			mp->layDim[l].nRows = mp->nRows;
			mp->layDim[l].nCols = mp->nCols;
			mp->layDim[l].nFilt = mp->nScales * mp->nOrients * mp->nPhases;
		}
		else
		{
			if (mp->vSquare[l]) // 2D (square) layer
			{
				mp->layDim[l].nRows = mp->layDim[l].nCols = round(sqrt(mp->vExcit[l]));
				mp->vExcit[l] = mp->layDim[l].nRows * mp->layDim[l].nCols;
				mp->layDim[l].nFilt = 1;
			}
			else // 1D (row) layer
			{
				mp->layDim[l].nRows = 1;
				mp->layDim[l].nCols = mp->vExcit[l];
				mp->layDim[l].nFilt = 1;
			}
		}
	}

	mp->spatialScale = MAX(mp->layDim[0].nRows, mp->layDim[0].nCols);
	mp->condSpeed = mp->spatialScale / (2.0 * mp->maxDelay); // maxDelay is across 1/2 largest dimension

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
		params->probConnect = true;
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
		params->probConnect = true;
	}
	
	if (mp->SOM && !mp->SOMinput)
		params->pCnxElE[0] = 0.0;
	
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
		params->probConnect = true;
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
		params->probConnect = true;
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
		params->probConnect = true;
	}
	
	assert(params->LpEfE == params->nLayers); // nWLayers
	assert(params->LpElE == params->nLayers);
	assert(params->LpEI == params->nLayers);
	assert(params->LpIE == params->nLayers);
	assert(params->LpII == params->nLayers);
	
	if (params->nRecordsPL && params->vRecords==NULL)
	{
		params->nRecords = 0;
		params->vRecords = myalloc(params->nLayers * params->nRecordsPL * sizeof(int));
		for (l=0; l<params->nLayers; l++)
		{
			params->vRecords[l] = params->nRecordsPL;
			assert((params->vRecords[l] >= 0) && (params->vRecords[l] <= params->vExcit[l]));
			params->nRecords += params->vRecords[l];
		}
	}
	
	if (params->vRecords)
	{
		params->nRecords = 0;
		for (l=0; l<params->nLayers; l++)
		{
			//assert((params->vRecords[l] >= 0) && (params->vRecords[l] <= params->vExcit[l]));
			if (params->vRecords[l] < 0)
				params->vRecords[l] = 0;
			if (params->vRecords[l] > params->vExcit[l])
				params->vRecords[l] = params->vExcit[l];
			params->nRecords += params->vRecords[l];
		}
	}
	
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
	strncpy(string, s1, 1+strlen(s1)); // +1 to copy terminator
	//strcpy(string, s1);
	return string;
}


int parseIntVector(char * string, int ** array)
{ // Function to parse a row vector of integers, store it in an array and return the number of elements
	char * tstr = NULL;
	const char * delims = "{[, ]}";
	int numElements = 0;
	int block = 1;
	
	myfree(*array); // Call free if *array!=NULL in case the array already exists
	*array = myalloc(VECBUFF*sizeof(int));
	tstr = strtok(string, delims); // Get first token starting with the first nondelimiter char
	
	while ( (tstr != NULL) && (*trim(tstr) != (']' || '}')) && (*tstr != '\0') ) // Trim and check not end
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
	const char * delims = "{[, ]}";
	int numElements = 0;
	int block = 1;
	
	myfree(*array); // Call free if *array!=NULL in case the array already exists
	*array = myalloc(VECBUFF*sizeof(float));
	tstr = strtok(string, delims); // Get first token starting with the first nondelimiter char
	
	while ( (tstr != NULL) && (*trim(tstr) != (']' || '}')) && (*tstr != '\0') ) // Trim and check not end
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

int parseVector(char * string, data * vectors, dataType type, const char * delims)
{ // Function to parse a row vector of integers, store it in an array and return the number of elements
	char * tstr = NULL;
	//const char * delims = "[, ]";
	int numElements = 0;
	int block = 1;
	if (!delims)
		char * delims = "[, ]";
	
	switch (dataType) {
		case 'int':
			dsize = sizeof(int);
			dRead = atoi();
			break;
		case 'float':
			dsize = sizeof(float);
			dRead = atof();
			break
		default:
			break;
	}
	
	free(array); // In case the array already exists
	*array = myalloc(VECBUFF*dsize); //sizeof(int *)
	tstr = strtok(string, delims); // Get first token starting with the first nondelimiter char
	
	while ( (*trim(tstr) != ']') && (*tstr != '\0') && (tstr != NULL) ) // Trim and check not end
	{
		if (numElements == block*VECBUFF) // Check array is not full up
			*array = myrealloc(*array, ++block*VECBUFF*dsize) // Reallocate memory
			
			*array[numElements++] = dRead(tstr); //atoi(tstr); // convert to int for assignment
		tstr = strtok(NULL, delims); // Get the next token
	}
	
	if (numElements == 0)
		fprintf(stderr, "Error! Vector: \"%s\" is empty",tstr);
	else
		if ( (*array = realloc(*array, numElements*dsize)) == NULL ) // Set array to correct size
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
	if (name[0] == 'M' && name[1] == 'P' && name[2] == '.')
		strncpy(name, &name[3], 1+strlen(name)-3); // +1 to copy terminator//name = name[3];
		
	trim(value);
	
	/* Copy into correct entry in parameters struct */ // Use strcmpi() for case insensitive comparisons
//#define Match(arg)	(strcmp(argv[cur_arg], (arg)) == 0)
//retval = (input) ? assign(param, value) : writeout(param, fp) ; 
	/* Simulation */		
	if (strcmp(name, "DT")==0)
		params->DT = atof(value);
	else if (strcmp(name, "loops")==0)
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
	else if (strcmp(name, "normalise")==0)
		params->normalise = atoi(value);
	else if (strcmp(name, "nRecordsPL")==0)
		params->nRecordsPL = atoi(value);
	else if (strcmp(name, "vRecords")==0)
		parseIntVector(value, &params->vRecords);
	else if (strcmp(name, "printConnections")==0)
		params->printConnections = atoi(value);
	else if (strcmp(name, "probConnect")==0)
		params->probConnect = atoi(value);
	
	/* Stimuli */
	else if ((strcmp(name, "imgList")==0)||(strcmp(name, "imageList")==0))
	{
		params->imgList = myalloc(strlen(trim(value))+1);
		strncpy(params->imgList, value, strlen(value));
		params->imgList[strlen(value)] = '\0'; // Needed for strncpy
		
		//myalloc(strlen(IMGDIR)+strlen(trim(value))+2); //*sizeof(char)
		/*strncpy(params->imgList, IMGDIR, strlen(IMGDIR));
		strncat(params->imgList, "/", 1);
		strncat(params->imgList, value, strlen(value));*/
		//strcpy(params->imageList, value);
		//params->imgList[strlen(IMGDIR)+1+strlen(value)] = '\0'; // Needed for strncpy
	}
	else if (strcmp(name, "randStimOrder")==0)
		params->randStimOrder = atoi(value);
	else if (strcmp(name, "randTransOrder")==0)
		params->randTransOrder = atoi(value);
	else if (strcmp(name, "randTransDirection")==0)
		params->randTransDirection = atoi(value);
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
	else if (strcmp(name, "M")==0)
		params->M = atoi(value);
	else if (strcmp(name, "K")==0)
		params->K = atoi(value);
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
		params->nScales = parseIntVector(value, &params->vScales);
	else if (strcmp(name, "vOrients")==0)
		params->nOrients = parseIntVector(value, &params->vOrients);
	else if (strcmp(name, "vPhases")==0)
		params->nPhases = parseIntVector(value, &params->vPhases);
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
		params->LvExcit = parseIntVector(value, &params->vExcit);
	else if (strcmp(name, "nSynEfE")==0)
		params->nSynEfE = atoi(value);
	else if (strcmp(name, "pCnxEfE")==0)
		params->LpEfE = parseFloatVector(value, &params->pCnxEfE);
	else if (strcmp(name, "nSynElE")==0)
		params->nSynElE = atoi(value);
	else if (strcmp(name, "pCnxElE")==0)
		params->LpElE = parseFloatVector(value, &params->pCnxElE);
	else if (strcmp(name, "nSynIE")==0)
		params->nSynIE = atoi(value);
	else if (strcmp(name, "pCnxIE")==0)
		params->LpIE = parseFloatVector(value, &params->pCnxIE);
	else if (strcmp(name, "nInhib")==0)
		params->nInhib = atoi(value);
	else if (strcmp(name, "rInhib")==0)
		params->rInhib = atof(value);
	else if (strcmp(name, "vInhib")==0)
		params->LvInhib = parseIntVector(value, &params->vInhib);
	else if (strcmp(name, "nSynEI")==0)
		params->nSynEI = atoi(value);
	else if (strcmp(name, "pCnxEI")==0)
		params->LpEI = parseFloatVector(value, &params->pCnxEI);
	else if (strcmp(name, "nSynII")==0)
		params->nSynII = atoi(value);
	else if (strcmp(name, "pCnxII")==0)
		params->LpII = parseFloatVector(value, &params->pCnxII);
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
	else if (strcmp(name, "spatialScale")==0)
		params->spatialScale = atof(value);
	else if (strcmp(name, "condSpeed")==0)
		params->condSpeed = atof(value);
	else if (strcmp(name, "maxDelay")==0)
		params->maxDelay = atof(value);
	else if (strcmp(name, "SOM")==0)
		params->SOM = atoi(value);
	else if (strcmp(name, "SOMinput")==0)
		params->SOMinput = atoi(value);
	else if (strcmp(name, "SOMsigE")==0)
		params->SOMsigE = atof(value);
	else if (strcmp(name, "SOMsigI")==0)
		params->SOMsigI = atof(value);
	else if (strcmp(name, "SOMclip")==0)
		params->SOMclip = atof(value);
	else if (strcmp(name, "vSquare")==0)
		parseIntVector(value, &params->vSquare);
	
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
	else if (strcmp(name, "adaptation")==0)
		params->adaptation = atoi(value);
	else if (strcmp(name, "alphaCa")==0)
		params->alphaCa = atof(value);
	else if (strcmp(name, "tauCa")==0)
		params->tauCa = atof(value);
	else if (strcmp(name, "gAHP")==0)
		params->gAHP = atof(value);
	else if (strcmp(name, "VK")==0)
		params->VK = atof(value);	
	
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
	else if (strcmp(name, "modEf")==0)
		params->modEf = atof(value);
	else if (strcmp(name, "tauEE")==0)
		params->tauEE = atof(value);
	else if (strcmp(name, "Dg_ElE")==0)
		params->Dg_ElE = atof(value);
	////////////////////////////////////// tauElE ?
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

int printParameters(PARAMS * mp, char * paramfile) // Update list of parameters
{
	FILE * pFile = NULL;
	int c=0; // count
	int l=0;
	PARAMS MP = *mp;
	
	pFile = myfopen(paramfile, "w"); // Variables to read into Matlab
	//fprintf(pFile, "DT=%f;\n",mp->DT); c++;
	FPRINT_FLOAT(pFile, MP.DT); c++;
	fprintf(pFile, "TotalMS=%d;\n",mp->TotalMS); c++;
	//char * rfpath = NULL;
	/*if (paramfile != NULL)
		sprintf(rfpath, "%s/%s/", RESDIR, paramfile);// strtok(pfile, "."));
	else
		sprintf(rfpath, "%s/", RESDIR);
	fprintf(pFile, "rfpath='%s';\n\n",rfpath); c++;*/

	fprintf(pFile, "\n%%%% Simulation Parameters %%%%\n");
	FPRINT_INT(pFile, MP.loops); c++; //fprintf(pFile, "MP.loops = %d;\n", mp->loops);
	FPRINT_INT(pFile, MP.train); c++;
	FPRINT_INT(pFile, MP.pretrain); c++;
	FPRINT_INT(pFile, MP.trainPause); c++;
	FPRINT_INT(pFile, MP.noise); c++;
	FPRINT_FLOAT(pFile, MP.noiseScale); c++;
	FPRINT_INT(pFile, MP.nRecordsPL); c++;
	FPRINT_INT(pFile, MP.nRecords); c++;
	if (mp->nRecords)
		printIntArray(pFile, "MP.vRecords", mp->vRecords, mp->nLayers);
	FPRINT_INT(pFile, MP.printConnections); c++;
	FPRINT_INT(pFile, MP.probConnect); c++;
	
	fprintf(pFile, "\n%%%% Stimulus Parameters %%%%\n");
	FPRINT_INT(pFile, MP.randStimOrder); c++;
	FPRINT_INT(pFile, MP.randTransOrder); c++;
	FPRINT_INT(pFile, MP.randTransDirection); c++;
	FPRINT_INT(pFile, MP.interleaveTrans); c++;
	FPRINT_INT(pFile, MP.localRep); c++;
	FPRINT_FLOAT(pFile, MP.current); c++;
	FPRINT_FLOAT(pFile, MP.currentSpread); c++;
	FPRINT_INT(pFile, MP.nStimuli);
	FPRINT_INT(pFile, MP.nTransPS);
	FPRINT_INT(pFile, MP.nRows);
	FPRINT_INT(pFile, MP.nCols);
	FPRINT_INT(pFile, MP.newTestSet); c++;
	if (mp->newTestSet)
	{
		FPRINT_INT(pFile, MP.nTestStimuli); c++;
		FPRINT_INT(pFile, MP.nTestTransPS); c++;
	}
	if (mp->M)
	{
		FPRINT_INT(pFile, MP.M); c++;
		FPRINT_INT(pFile, MP.K); c++;
	}
	FPRINT_FLOAT(pFile, MP.transP_Train); c++;
	FPRINT_FLOAT(pFile, MP.transP_Test); c++;
	FPRINT_INT(pFile, MP.useFilteredImages); c++;
	if (mp->useFilteredImages) // load other parameters from imageParams.m
	{
		FPRINT_STRING(pFile, MP.imgList); c++;
		FPRINT_INT(pFile, MP.gabor); c++;
		printIntArray(pFile, "MP.vPhases", mp->vPhases, mp->nPhases); c++;
		printIntArray(pFile, "MP.vScales", mp->vScales, mp->nScales); c++;
		printIntArray(pFile, "MP.vOrients", mp->vOrients, mp->nOrients); c++;
	}
	else
	{
		FPRINT_INT(pFile, MP.shift); c++;
		FPRINT_INT(pFile, MP.nFiringNeurons); c++;
		FPRINT_FLOAT(pFile, MP.a); c++;
	}
		
	
	fprintf(pFile, "\n%%%% Network Parameters %%%%\n");
	FPRINT_INT(pFile, MP.nLayers); c++;
	FPRINT_INT(pFile, MP.nWLayers); c++;
	FPRINT_INT(pFile, MP.inputInhib); c++;
	printIntArray(pFile, "MP.vExcit", mp->vExcit, mp->nLayers); c++;
	printIntArray(pFile, "MP.vInhib", mp->vInhib, mp->nLayers); c++;
	FPRINT_FLOAT(pFile, MP.rInhib); c++;
	printFloatArray(pFile, "MP.pCnxEfE", mp->pCnxEfE, mp->nLayers); c++;
	printFloatArray(pFile, "MP.pCnxElE", mp->pCnxElE, mp->nLayers); c++;
	printFloatArray(pFile, "MP.pCnxIE", mp->pCnxIE, mp->nLayers); c++; // nSyn?
	printFloatArray(pFile, "MP.pCnxEI", mp->pCnxEI, mp->nLayers); c++;
	printFloatArray(pFile, "MP.pCnxII", mp->pCnxII, mp->nLayers); c++;
	switch (mp->axonDelay)
	{
		case MinD:
			fprintf(pFile, "MP.axonDelay = 'MinD';\n"); c++;
			//FPRINT_FLOAT(pFile, MP.d_const); c++;
			break;
		case ConstD:
			fprintf(pFile, "MP.axonDelay = 'ConstD';\n"); c++;
			FPRINT_FLOAT(pFile, MP.d_const); c++;
			break;
		case UniformD:
			fprintf(pFile, "MP.axonDelay = 'UniformD';\n"); c++;
			FPRINT_FLOAT(pFile, MP.d_min); c++;
			FPRINT_FLOAT(pFile, MP.d_max); c++;
			break;
		case GaussD:
			fprintf(pFile, "MP.axonDelay = 'GaussD';\n"); c++;
			FPRINT_FLOAT(pFile, MP.d_mean); c++;
			FPRINT_FLOAT(pFile, MP.d_sd); c++;
			break;
		case SOM:
			fprintf(pFile, "MP.axonDelay = 'SOM';\n"); c++;
			FPRINT_FLOAT(pFile, MP.condSpeed); c++;
			FPRINT_FLOAT(pFile, MP.maxDelay); c++;
		default:
			break;
	}
	
	FPRINT_FLOAT(pFile, MP.spatialScale); c++; // Needed for plotting connectivity
	
	FPRINT_INT(pFile, MP.SOM); c++;
	if (MP.SOM)
	{
		FPRINT_INT(pFile, MP.SOMinput); c++;
		FPRINT_FLOAT(pFile, MP.SOMsigE); c++;
		FPRINT_FLOAT(pFile, MP.SOMsigI); c++;
		FPRINT_FLOAT(pFile, MP.SOMclip); c++;
	}
	printIntArray(pFile, "MP.vSquare", mp->vSquare, mp->nLayers); c++;
	for (l=0; l<MP.nLayers; l++)
	{
		fprintf(pFile, "%%%% Layer %d Dimensions %%%%\n", l);
		fprintf(pFile, "%% MP.layDim[%d].nRows = %d;\n",l,MP.layDim[l].nRows); c++;
		fprintf(pFile, "%% MP.layDim[%d].nCols = %d;\n",l,MP.layDim[l].nCols); c++;
		fprintf(pFile, "%% MP.layDim[%d].nFilt = %d;\n",l,MP.layDim[l].nFilt); c++;
	}
	
	/* Cell bodies */
	fprintf(pFile, "\n%%%% Cell Body Parameters %%%%\n");
	FPRINT_FLOAT(pFile, MP.capE); c++;
	FPRINT_FLOAT(pFile, MP.capI); c++;
	FPRINT_FLOAT(pFile, MP.gLeakE); c++;
	FPRINT_FLOAT(pFile, MP.gLeakI); c++;
	FPRINT_FLOAT(pFile, MP.VrestE); c++;
	FPRINT_FLOAT(pFile, MP.VrestI); c++;
	FPRINT_FLOAT(pFile, MP.VhyperE); c++;
	FPRINT_FLOAT(pFile, MP.VhyperI); c++;
	FPRINT_FLOAT(pFile, MP.ThreshE); c++;
	FPRINT_FLOAT(pFile, MP.ThreshI); c++;
	FPRINT_FLOAT(pFile, MP.VrevE); c++;
	FPRINT_FLOAT(pFile, MP.VrevI); c++;
	FPRINT_FLOAT(pFile, MP.refract); c++;
	FPRINT_INT(pFile, MP.adaptation); c++;
	if (MP.adaptation)
	{
		FPRINT_FLOAT(pFile, MP.alphaCa); c++;
		FPRINT_FLOAT(pFile, MP.tauCa); c++;
		FPRINT_FLOAT(pFile, MP.gAHP); c++;
		FPRINT_FLOAT(pFile, MP.VK); c++;
	}
	
	/* Synapses */
	fprintf(pFile, "\n%%%% Synapse Parameters %%%%\n");
	FPRINT_FLOAT(pFile, MP.alphaC); c++;
	FPRINT_FLOAT(pFile, MP.tauC); c++;
	FPRINT_FLOAT(pFile, MP.alphaD); c++;
	FPRINT_FLOAT(pFile, MP.tauD); c++;
	FPRINT_FLOAT(pFile, MP.learnR); c++;
	FPRINT_FLOAT(pFile, MP.modEf); c++;
	FPRINT_FLOAT(pFile, MP.tauEE); c++;
	FPRINT_FLOAT(pFile, MP.Dg_ElE); c++;
	FPRINT_FLOAT(pFile, MP.Dg_IE); c++;
	FPRINT_FLOAT(pFile, MP.tauIE); c++;
	FPRINT_FLOAT(pFile, MP.Dg_EI); c++;
	FPRINT_FLOAT(pFile, MP.tauEI); c++;
	FPRINT_FLOAT(pFile, MP.Dg_II); c++;
	FPRINT_FLOAT(pFile, MP.tauII); c++;
	FPRINT_FLOAT(pFile, MP.gMax); c++;

	fclose(pFile);
	return c;
}

void printIntArray(FILE * fp, char * name, int * array, int len)
{
	if (!len || !array)
		return;
	int l=0;
	fprintf(fp, "%s = [%d",name,array[0]);
	for (l=1; l<len; l++)
		fprintf(fp, ",%d", array[l]);
	fprintf(fp, "];\n");
	fflush(fp);	
}

void printFloatArray(FILE * fp, char * name, float * array, int len)
{
	if (!len || !array)
		return;
	int l=0;
	fprintf(fp, "%s = [%G",name,array[0]);
	for (l=1; l<len; l++)
		fprintf(fp, ",%G", array[l]);
	fprintf(fp, "];\n");
	fflush(fp);	
	
}
