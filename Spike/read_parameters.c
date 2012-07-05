/*
 *  read_parameters.c
 *  Spike
 *
 *  Created by Ben Evans on 28/11/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "read_parameters.h"
extern unsigned long int seed;
extern int nThreads;
extern char * pFile;

int read_parameters(PARAMS * params, char * paramfile)
{
	int count = 0;
	int l = 0;
	//int square = 0;
	char * string, buff[BUFFER];
	
	int prevNLayers = (mp->initialised) ? mp->nLayers : 0;
	
	FILE * fp = myfopen(paramfile, "r");
	while ((string = fgets(buff, sizeof(buff), fp)) != NULL) // Read each line
		count += parse_string(params, buff);
	fclose(fp);
	
	if (params->useFilteredImages && !(strcmp(paramfile, MPFILE)==0)) // Do not parse IPFILE if rerunning with parameters.m
	{
		int slen = strlen(mp->imgDir)+1+strlen(IPFILE)+1; // '/' & '\0'
		char * imgParams = myalloc(slen);
		if (snprintf(imgParams, slen, "%s/%s", mp->imgDir, IPFILE) >= slen)
			fprintf(stderr, "*** Warning! Undersized buffer: %s ***", imgParams);
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
	else // *** RETHINK this. Replace all with transSize?
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
	
	assert(mp->nLayers > 0); // nLayers must be set in (default) paramfile
	int deltaLayers = (mp->initialised) ? (mp->nLayers - prevNLayers) : 0;
	
	/* Additional calculations */
	
	/*if (!mp->stimGroups)
	{*/
		if (mp->K > 1)		// Assumption: when training with multiples, test with individual stimuli
		{
			mp->newTestSet = true;
			mp->nTestStimuli = mp->nStimuli;
			int nCombs = gsl_sf_choose(mp->nTestStimuli-1, mp->K-1);
			if (mp->M > nCombs)
				mp->M = nCombs;
			mp->M = (mp->nStimuli == 1) ? 1 : mp->M; // Just present the same stimuli again in a loop instead if M>nStimuli?
			mp->nStimuli = mp->nTestStimuli * mp->M;
			
			/*nCombs = gsl_sf_choose(mp->nTestStimuli, mp->K);
			 mp->nStimuli = (nCombs == 1) ? nCombs : mp->nStimuli;*/
			//mp->nStimuli = mp->M * gsl_sf_choose(mp->nTestStimuli, mp->K) / mp->nTestStimuli; // Check ////////
			
			mp->nTestTransPS = mp->nTransPS;
		}
		/*else
		{
			mp->newTestSet = false;
			mp->nTestStimuli = mp->nStimuli;
			mp->nTestTransPS = mp->nTransPS;
		}
	}*/
	
	if (!mp->newTestSet)
	{
		mp->nTestStimuli = mp->nStimuli;
		mp->nTestTransPS = mp->nTransPS;
	}
		
	mp->EpochTime = mp->nStimuli * mp->nTransPS * mp->transP_Train; //Train period per loop
	mp->EpochMS = round(mp->EpochTime * 1000);
	mp->TestTime = mp->nTestStimuli * mp->nTestTransPS * mp->transP_Test;
	mp->TestMS = round(mp->TestTime * 1000);
	mp->MaxTime = MAX(mp->TestTime, mp->EpochTime);
	//float transP = (params->transP_Train >= params->transP_Test ? params->transP_Train : params->transP_Test);
	//params->TotalTime = params->nStimuli * params->nTransPS * transP;
	mp->TotalMS = round(mp->MaxTime * 1000); // Was ceil, used for recording structures
	mp->RecordMS = mp->TotalMS; //Relabel these
	//params->TotalTS = round(params->TotalTime/mp->DT); // Was ceil
	params->TSperMS = round(1/(mp->DT*1000)); // Was ceil, used for recording structures
	mp->spkBuffer = ceil(mp->MaxTime / mp->refract);
	params->SigmaE = params->noiseScale * (params->ThreshE - params->VhyperE);
	params->SigmaI = params->noiseScale * (params->ThreshI - params->VhyperI);
	//params->inpSpkBuff = ceil(transP / params->refract);
	
    if (params->interleaveTrans) // Stop resetting between stimuli (redundant?)
        params->trainPause = false;
    

	
	// Change all myalloc to myrealloc?
    
    // mp->initialised may make Lv<>'s redundant... collapse into if (mp->initialised) {} else {} and move up
    // Hard code nLayers = 2 to start to avoid changes when initialising with defaults.m?

    if (deltaLayers != 0) //(deltaLayers < 0 || deltaLayers > 0) only if mp->initialised
    {
        assert(mp->nLayers > 0);
        mp->vRecords = myrealloc(mp->vRecords, mp->nLayers * sizeof(mp->vRecords[0]));
        mp->vExcit = myrealloc(mp->vExcit, mp->nLayers * sizeof(mp->vExcit[0]));
        mp->vInhib = myrealloc(mp->vInhib, mp->nLayers * sizeof(mp->vInhib[0]));
        mp->pCnxEfE = myrealloc(mp->pCnxEfE, mp->nLayers * sizeof(mp->pCnxEfE[0]));
        mp->pCnxElE = myrealloc(mp->pCnxElE, mp->nLayers * sizeof(mp->pCnxElE[0]));
        mp->pCnxIE = myrealloc(mp->pCnxIE, mp->nLayers * sizeof(mp->pCnxIE[0]));
        mp->pCnxEI = myrealloc(mp->pCnxEI, mp->nLayers * sizeof(mp->pCnxEI[0]));
        mp->pCnxII = myrealloc(mp->pCnxII, mp->nLayers * sizeof(mp->pCnxII[0]));
        mp->layDim = myrealloc(mp->layDim, mp->nLayers * sizeof(mp->layDim[0]));
        mp->vSquare = myrealloc(mp->vSquare, mp->nLayers * sizeof(mp->vSquare[0]));
		mp->LvExcit = mp->LvInhib = mp->LpEfE = mp->LpElE = mp->LpIE = mp->LpEI = mp->LpII = mp->nLayers;
        
        if (deltaLayers > 0) // Layers added 
        {
            int prevN = prevNLayers; //mp->nLayers - deltaLayers; // Or pass this instead...
            for (l=prevN; l<mp->nLayers; l++)
            {
                mp->vRecords[l] = mp->vRecords[prevN-1];
                mp->vExcit[l] = mp->vExcit[prevN-1];
                mp->vInhib[l] = mp->vInhib[prevN-1];
                mp->pCnxEfE[l] = mp->pCnxEfE[prevN-1];
                mp->pCnxElE[l] = mp->pCnxElE[prevN-1];
                mp->pCnxIE[l] = mp->pCnxIE[prevN-1];
                mp->pCnxEI[l] = mp->pCnxEI[prevN-1];
                mp->pCnxII[l] = mp->pCnxII[prevN-1];
                mp->layDim[l].nCols = mp->layDim[prevN-1].nCols;
                mp->layDim[l].nRows = mp->layDim[prevN-1].nRows;
                mp->layDim[l].nFilt = (prevN>1) ? mp->layDim[prevN-1].nFilt : 1;
                mp->vSquare[l] = mp->vSquare[prevN-1];
            }
        }
        deltaLayers = 0;//prevN = mp->nLayers;
    }
    else // No change - perform consistency checks
    {
        //mp->nWLayers = mp->nLayers;
    }
	
	/* Parse vectors of neurons in each layer */
	if (params->LvExcit == 0) // Initialize to 0 in defaults.m // Need to consider old parameter files where nExcit is passed after initialization
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
			assert((0.0 <= params->pCnxEfE[l]) && (params->pCnxEfE[l] <= 1.0));
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
			assert((0.0 <= params->pCnxElE[l]) && (params->pCnxElE[l] <= 1.0));
			if (params->pCnxElE[l] > EPS) // 
				params->SOM = true;
		}
		params->LpElE = l;
		params->probConnect = true;
	}
	
	float latprob=0.0;
	for (l=0; l<params->nLayers; l++)
		latprob += params->pCnxElE[l];
	mp->SOM = (latprob>EPS) ? true : false;
	
	if (mp->initElE || mp->axonDelay) // Include delayEfE, delayElE, delayEI
	{
		mp->SOM = true;
		if (!mp->SOMinput)
			params->pCnxElE[0] = 0.0;
	}
		
	
	/*if (mp->SOM && !mp->SOMinput) //(mp->pCnxElE[0] > EPS && !mp->SOMinput)
		params->pCnxElE[0] = 0.0;*/
	
	/*if (mp->SOM)
		for (l=(mp->SOMinput)?0:1; l<mp->nLayers; l++)
			params->pCnxElE[l] = 1.0;*/
	
	if (params->LpEI == 0)
	{
		assert(params->nSynEI <= params->nExcit);
		params->pCnxEI = myalloc(params->nLayers*sizeof(float));
		for (l=0; l<params->nLayers; l++)
		{
			params->pCnxEI[l] = (params->vInhib[l]>0) ? ((float) params->nSynEI)/params->vExcit[l] : 0.0;
			assert((0.0 <= params->pCnxEI[l]) && (params->pCnxEI[l] <= 1.0));
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
			assert((0.0 <= params->pCnxIE[l]) && (params->pCnxIE[l] <= 1.0));
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
			assert((0.0 <= params->pCnxII[l]) && (params->pCnxII[l] <= 1.0));
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
		params->vRecords = myalloc(params->nLayers * sizeof(int)); // * params->nRecordsPL
		for (l=0; l<params->nLayers; l++)
		{
			params->vRecords[l] = params->nRecordsPL;
			assert((0 <= params->vRecords[l]) && (params->vRecords[l] <= params->vExcit[l]));
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
	
	if (mp->stimGroups && strcmp(STFILE, "")==0) // No STFILE name passed
	{
		//slen = strlen(value);
		//strncpy(STFILE, value, slen*sizeof(char)); // copies <= slen bytes (bytes following null byte are not copied)
		//STFILE[slen] = '\0';
		FILEPARTS * fp = myalloc(sizeof(FILEPARTS));
		getFileParts(pFile, fp);
		strncpy(STFILE, fp->fname, BUFSIZ);
		strncat(STFILE, ".stm", BUFSIZ-strlen(STFILE));
		assert(strlen(STFILE)<BUFSIZ);
		if (mp->priorPhases && strcmp(PPSTFILE, "")==0) // If PPstimFile is not passed, assume it is PP<STFILE>
		{
			strncpy(PPSTFILE, "PP", BUFSIZ);
			strncat(PPSTFILE, fp->fname, BUFSIZ-strlen(PPSTFILE));
			strncat(PPSTFILE, ".stm", BUFSIZ-strlen(PPSTFILE));
			assert(strlen(PPSTFILE)<BUFSIZ);
		}
		myfree(fp);
	}
	
	// Check membrane capacitances and leakages for membrane time constants
	if (mp->gLeakE && (mp->capE/mp->gLeakE < SIM.minTau))
		SIM.minTau = mp->capE/mp->gLeakE;
	
	if (mp->gLeakI && (mp->capI/mp->gLeakI < SIM.minTau))
		SIM.minTau = mp->capI/mp->gLeakI;
	
	if (mp->train)
		mp->trainEfE = true; // By default
	
	if (count > 0)
		mp->initialised = true;
	
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
	//int slen=0;
	
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
	if (strcmp(name, "mSeed")==0)
        printf("\n<Original mSeed: %ld>\n",atol(value));
    else if (strcmp(name, "nThreads")==0)
        printf("<Original nThreads: %d>\n",atoi(value));
    else if (strcmp(name, "RecordMS")==0 || strcmp(name, "EpochMS")==0 || strcmp(name, "TestMS")==0)
        ; // Skip - Matlab variables
	else if (strcmp(name, "DT")==0)
		params->DT = atof(value);
	else if (strcmp(name, "loops")==0 || strcmp(name, "nEpochs")==0)
		params->loops = atoi(value);
	else if (strcmp(name, "pretrain")==0)
		params->pretrain = atoi(value);
	else if (strcmp(name, "train")==0)
		params->train = atoi(value);
	else if (strcmp(name, "priorPhases")==0)
		params->priorPhases = atoi(value);
	else if (strcmp(name, "isolateEfE")==0)
		mp->isolateEfE = atoi(value);
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
    else if (strcmp(name, "nRecords")==0)
        ; // Skip - Matlab variable
	else if (strcmp(name, "vRecords")==0)
		parseIntVector(value, &params->vRecords);
	else if (strcmp(name, "printConnections")==0)
		params->printConnections = atoi(value);
	else if (strcmp(name, "probConnect")==0)
		params->probConnect = atoi(value);
	else if (strcmp(name, "loadWeights")==0)
		params->loadWeights = atoi(value);
	
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
	else if (strcmp(name, "stimFile")==0)
    {
		/*slen = strlen(value);
		STFILE = myalloc((slen+1)*sizeof(char));
		//fprintf(stdout, "BEFORE: dest=%s; src=%s; len=%d",STFILE,value,slen); fflush(stdout);
		strncpy(STFILE, value, slen*sizeof(char)); // copies <= slen bytes (bytes following null byte are not copied)
		STFILE[slen] = '\0';
		//fprintf(stdout, "AFTER: dest=%s; src=%s; len=%d",STFILE,value,slen); fflush(stdout);*/	
		
		strncpy(STFILE, value, BUFSIZ); // Check NULL byte is added
		assert(strlen(STFILE)<BUFSIZ);
		// If PPstimFile is not passed, assume it is PP<STFILE>
		strncpy(PPSTFILE, "PP", BUFSIZ);
		strncat(PPSTFILE, value, BUFSIZ-strlen(PPSTFILE));
		assert(strlen(PPSTFILE)<BUFSIZ);
    }
	
	else if (strcmp(name, "PPstimFile")==0)
	{
		/*slen = strlen(value);
		myfree(ELESTFILE); // In case the default filename has already been allocated
		ELESTFILE = myalloc((slen+1)*sizeof(char));
		strncpy(ELESTFILE, value, slen*sizeof(char));
		ELESTFILE[slen] = '\0';*/
		
		strncpy(PPSTFILE, value, BUFSIZ);
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
	else if (strcmp(name, "stimGroups")==0)
		params->stimGroups = atoi(value);
	else if (strcmp(name, "nBG")==0)
		params->nBG = atoi(value);
	else if (strcmp(name, "nWG")==0)
		params->nWG = atoi(value);
	else if (strcmp(name, "nGroups")==0) //Internal
		params->nGroups = (atoi(value) < 1) ? 1 : atoi(value);
	else if (strcmp(name, "nStimuli")==0)
		params->nStimuli = (atoi(value) < 1) ? 1 : atoi(value);
	else if (strcmp(name, "nTransPS")==0)
		params->nTransPS = (atoi(value) < 1) ? 1 : atoi(value);
	else if (strcmp(name, "newTestSet")==0)
		params->newTestSet = atoi(value);
	else if (strcmp(name, "M")==0)
		params->M = atoi(value);
	else if (strcmp(name, "K")==0)
		params->K = atoi(value);
	/*else if (strcmp(name, "nTestGroups")==0)
		params->nTestGroups = (atoi(value) < 1) ? 1 : atoi(value);*/
	else if (strcmp(name, "nTestStimuli")==0)
		params->nTestStimuli = (atoi(value) < 1) ? 1 : atoi(value);
	else if (strcmp(name, "nTestTransPS")==0)
		params->nTestTransPS = (atoi(value) < 1) ? 1 : atoi(value);
	else if (strcmp(name, "transP_Train")==0)
		params->transP_Train = atof(value);
	else if (strcmp(name, "transP_Test")==0)
		params->transP_Test = atof(value);
	else if (strcmp(name, "shift")==0)
		params->shift = atoi(value);
    else if (strcmp(name, "nFiringNeurons")==0)
    {
        params->nFiringNeurons = atoi(value);
        params->a = (float) params->nFiringNeurons / params->sInputs;
    }
	else if (strcmp(name, "a")==0)
	{
		params->a = atof(value);
		assert(0 <= mp->a && mp->a <= 1);
	}
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
		// change = new - old
		// realloc all vectors
		params->nLayers = atoi(value);
		params->nWLayers = params->nLayers - 1;
	}
    else if (strcmp(name, "nWLayers")==0)
        assert(params->nWLayers == params->nLayers-1);
	else if (strcmp(name, "inputInhib")==0)
		params->inputInhib = atoi(value);
	else if (strcmp(name, "nExcit")==0)
		params->nExcit = atoi(value);
	else if (strcmp(name, "sInputs")==0)
		params->vExcit[0] = atoi(value);
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
	else if (strcmp(name, "initEfE")==0)
		params->initEfE = atoi(value);
	else if (strcmp(name, "iEfE")==0)
		params->iEfE = atof(value);
	else if (strcmp(name, "axonDelay")==0)
		params->axonDelay = atoi(value);
    else if (strcmp(name, "delayEfE")==0)
    {
        params->delayEfE = atoi(value);
        params->axonDelay = true;
    }
    else if (strcmp(name, "delayElE")==0)
    {
        params->delayElE = atoi(value);
        params->axonDelay = true;
    }
    else if (strcmp(name, "delayEI")==0)
    {
        params->delayEI = atoi(value);
        params->axonDelay = true;
    }
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
	else if (strcmp(name, "SOM")==0) // Internal variable
		printf("SOM is now an internal variable!\n");//params->SOM = atoi(value);
	else if (strcmp(name, "SOMinput")==0)
		params->SOMinput = atoi(value);
	else if (strcmp(name, "SOMsigE")==0)
		params->SOMsigE = atof(value);
	else if (strcmp(name, "SOMsigI")==0)
		params->SOMsigI = atof(value);
	else if (strcmp(name, "SOMclip")==0)
		params->SOMclip = atof(value);
	else if (strcmp(name, "trainElE")==0)
		params->trainElE = atoi(value);
	else if (strcmp(name, "initElE")==0)
		params->initElE = atoi(value);
	else if (strcmp(name, "iElE")==0)
		params->iElE = atof(value);
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
	{
		params->tauCa = atof(value);
		if (mp->tauCa < SIM.minTau)
			SIM.minTau = mp->tauCa;
	}
		
	else if (strcmp(name, "gAHP")==0)
		params->gAHP = atof(value);
	else if (strcmp(name, "VK")==0)
		params->VK = atof(value);	
	
	/* Synapses */
	else if (strcmp(name, "alphaC")==0)
		params->alphaC = atof(value);
	else if (strcmp(name, "tauC")==0)
	{
		params->tauC = atof(value);
		if (mp->tauC < SIM.minTau)
			SIM.minTau = mp->tauC;
	}
	else if (strcmp(name, "alphaD")==0)
		params->alphaD = atof(value);
	else if (strcmp(name, "tauD")==0)
	{
		params->tauD = atof(value);
		if (mp->tauD < SIM.minTau)
			SIM.minTau = mp->tauD;
	}
	else if (strcmp(name, "learnR")==0)
	{
		params->learnR = atof(value);
		assert(0 <= mp->learnR && mp->learnR <= 1);
	}
	else if (strcmp(name, "DgEfE")==0 || strcmp(name, "modEf")==0)
		params->DgEfE = atof(value);
	else if (strcmp(name, "tauEfE")==0 || strcmp(name, "tauEE")==0) // Set both to tauEE if used
	{
		params->tauEE = atof(value); // Implement tauEfE		
		if (mp->tauEE < SIM.minTau)
			SIM.minTau = mp->tauEE;
	}
	else if (strcmp(name, "DgElE")==0 || strcmp(name, "Dg_ElE")==0)
		params->DgElE = atof(value);
	else if (strcmp(name, "tauElE")==0 || strcmp(name, "tauEE")==0) // Set both to tauEE if used
	{
		params->tauEE = atof(value); // Implement tauElE
		if (mp->tauEE < SIM.minTau)
			SIM.minTau = mp->tauEE;
	}
	else if (strcmp(name, "DgIE")==0 || strcmp(name, "Dg_IE")==0)
		params->DgIE = atof(value);
	else if (strcmp(name, "tauIE")==0)
	{
		params->tauIE = atof(value);
		if (mp->tauIE < SIM.minTau)
			SIM.minTau = mp->tauIE;
	}
	else if (strcmp(name, "DgEI")==0 || strcmp(name, "Dg_EI")==0)
		params->DgEI = atof(value);
	else if (strcmp(name, "tauEI")==0)
	{
		params->tauEI = atof(value);
		if (mp->tauEI < SIM.minTau)
			SIM.minTau = mp->tauEI;
	}
	else if (strcmp(name, "DgII")==0 || strcmp(name, "Dg_II")==0)
		params->DgII = atof(value);
	else if (strcmp(name, "tauII")==0)
	{
		params->tauII = atof(value);
		if (mp->tauII < SIM.minTau)
			SIM.minTau = mp->tauII;
	}
	else if (strcmp(name, "gMax")==0)
		params->gMax = atof(value);
	
	else
	{
		fprintf(stderr, "*** WARNING: %s/%s: Unknown name/value pair! ***\n", name, value);
		exit (-1); //return -1;
	}
	
	return 1; // Parameter succesfully loaded
}

int printParameters(PARAMS * mp, char * paramfile) // Update list of parameters
{
	FILE * pFP = NULL;
	int c=0; // count
	int l=0;
	PARAMS MP = *mp; // Needed to match Matlab syntax for reading in variables
	
	pFP = myfopen(paramfile, "w"); // Variables to read into Matlab
	//fprintf(pFile, "DT=%f;\n",mp->DT); c++;
	fprintf(pFP, "mSeed = %lu;\n", seed);
	fprintf(pFP, "nThreads = %d;\n", nThreads);
	FPRINT_FLOAT(pFP, MP.DT); c++; fprintf(pFP, "%% %f ms\n",MP.DT*1000);
	fprintf(pFP, "RecordMS = %d;\n",mp->RecordMS); c++;
	fprintf(pFP, "EpochMS = %d;\n",mp->EpochMS); c++;
	fprintf(pFP, "TestMS = %d;\n",mp->TestMS); c++;
	//char * rfpath = NULL;
	/*if (paramfile != NULL)
		sprintf(rfpath, "%s/%s/", RESDIR, paramfile);// strtok(pfile, "."));
	else
		sprintf(rfpath, "%s/", RESDIR);
	fprintf(pFile, "rfpath='%s';\n\n",rfpath); c++;*/

	fprintf(pFP, "\n%%%% Simulation Parameters %%%%\n");
	FPRINT_INT(pFP, MP.loops); c++; //fprintf(pFile, "MP.loops = %d;\n", mp->loops);
	FPRINT_INT(pFP, MP.pretrain); c++;
	FPRINT_INT(pFP, MP.train); c++;
	FPRINT_INT(pFP, MP.priorPhases); c++;
	FPRINT_INT(pFP, MP.isolateEfE); c++;
	FPRINT_INT(pFP, MP.trainPause); c++;
	FPRINT_INT(pFP, MP.noise); c++;
	FPRINT_FLOAT(pFP, MP.noiseScale); c++;
	FPRINT_INT(pFP, MP.nRecordsPL); c++;
	FPRINT_INT(pFP, MP.nRecords); c++;
	if (mp->nRecords)
		printIntArray(pFP, "MP.vRecords", mp->vRecords, mp->nLayers);
	FPRINT_INT(pFP, MP.printConnections); c++;
	FPRINT_INT(pFP, MP.probConnect); c++;
	FPRINT_INT(pFP, MP.loadWeights); c++;
	
	fprintf(pFP, "\n%%%% Stimulus Parameters %%%%\n");
	FPRINT_INT(pFP, MP.randStimOrder); c++;
	FPRINT_INT(pFP, MP.randTransOrder); c++;
	FPRINT_INT(pFP, MP.randTransDirection); c++;
	FPRINT_INT(pFP, MP.interleaveTrans); c++;
	FPRINT_INT(pFP, MP.localRep); c++;
	FPRINT_FLOAT(pFP, MP.current); c++;
	FPRINT_FLOAT(pFP, MP.currentSpread); c++;
	FPRINT_INT(pFP, MP.stimGroups); c++;
	if (MP.stimGroups)
	{
		FPRINT_INT(pFP, MP.nBG); c++;
		FPRINT_INT(pFP, MP.nWG); c++;
		FPRINT_INT(pFP, MP.nGroups); c++;
	}
	FPRINT_INT(pFP, MP.nStimuli); c++;
	FPRINT_INT(pFP, MP.nTransPS); c++;
	FPRINT_INT(pFP, MP.nRows); c++;
	FPRINT_INT(pFP, MP.nCols); c++;
	FPRINT_INT(pFP, MP.newTestSet); c++;
	//if (mp->newTestSet)
	//{
		FPRINT_INT(pFP, MP.nTestStimuli); c++;
		FPRINT_INT(pFP, MP.nTestTransPS); c++;
	//}
	if (mp->M)
	{
		FPRINT_INT(pFP, MP.M); c++;
		FPRINT_INT(pFP, MP.K); c++;
	}
	//FPRINT_INT(pFile, MP.nTestGroups); c++;
	FPRINT_FLOAT(pFP, MP.transP_Train); c++;
	FPRINT_FLOAT(pFP, MP.transP_Test); c++;
	FPRINT_INT(pFP, MP.useFilteredImages); c++;
	if (mp->useFilteredImages) // load other parameters from imageParams.m
	{
		FPRINT_STRING(pFP, MP.imgList); c++;
		FPRINT_INT(pFP, MP.gabor); c++;
		printIntArray(pFP, "MP.vPhases", mp->vPhases, mp->nPhases); c++;
		printIntArray(pFP, "MP.vScales", mp->vScales, mp->nScales); c++;
		printIntArray(pFP, "MP.vOrients", mp->vOrients, mp->nOrients); c++;
	}
	else
	{
		FPRINT_INT(pFP, MP.shift); c++;
		FPRINT_INT(pFP, MP.nFiringNeurons); c++;
		FPRINT_FLOAT(pFP, MP.a); c++;
	}
		
	
	fprintf(pFP, "\n%%%% Network Parameters %%%%\n");
	FPRINT_INT(pFP, MP.nLayers); c++;
	FPRINT_INT(pFP, MP.nWLayers); c++;
	FPRINT_INT(pFP, MP.inputInhib); c++;
	FPRINT_INT(pFP, MP.sInputs); c++;
	printIntArray(pFP, "MP.vExcit", mp->vExcit, mp->nLayers); c++;
	printIntArray(pFP, "MP.vInhib", mp->vInhib, mp->nLayers); c++;
	FPRINT_FLOAT(pFP, MP.rInhib); c++;
	printFloatArray(pFP, "MP.pCnxEfE", mp->pCnxEfE, mp->nLayers); c++;
	printFloatArray(pFP, "MP.pCnxElE", mp->pCnxElE, mp->nLayers); c++;
	printFloatArray(pFP, "MP.pCnxIE", mp->pCnxIE, mp->nLayers); c++; // nSyn?
	printFloatArray(pFP, "MP.pCnxEI", mp->pCnxEI, mp->nLayers); c++;
	printFloatArray(pFP, "MP.pCnxII", mp->pCnxII, mp->nLayers); c++;

    FPRINT_INT(pFP, MP.axonDelay); c++;
    
    if (mp->axonDelay)
    {
        switch (mp->delayEfE) //(mp->axonDelay)
        {
            case MinD:
                fprintf(pFP, "MP.delayEfE = 'MinD';\n"); c++;
                //FPRINT_FLOAT(pFile, MP.d_const); c++;
                break;
            case ConstD:
                fprintf(pFP, "MP.delayEfE = 'ConstD';\n"); c++;
                FPRINT_FLOAT(pFP, MP.d_const); c++;
                break;
            case UniformD:
                fprintf(pFP, "MP.delayEfE = 'UniformD';\n"); c++;
                FPRINT_FLOAT(pFP, MP.d_min); c++;
                FPRINT_FLOAT(pFP, MP.d_max); c++;
                break;
            case GaussD:
                fprintf(pFP, "MP.delayEfE = 'GaussD';\n"); c++;
                FPRINT_FLOAT(pFP, MP.d_mean); c++;
                FPRINT_FLOAT(pFP, MP.d_sd); c++;
                break;
            case SOMD:
                fprintf(pFP, "MP.delayEfE = 'SOM';\n"); c++; // SOMD would break Matlab code
                FPRINT_FLOAT(pFP, MP.condSpeed); c++;
                FPRINT_FLOAT(pFP, MP.maxDelay); c++;
            default:
                break;
        }
        
        switch (mp->delayElE) //(mp->axonDelay)
        {
            case MinD:
                fprintf(pFP, "MP.delayElE = 'MinD';\n"); c++;
                //FPRINT_FLOAT(pFile, MP.d_const); c++;
                break;
            case ConstD:
                fprintf(pFP, "MP.delayElE = 'ConstD';\n"); c++;
                FPRINT_FLOAT(pFP, MP.d_const); c++;
                break;
            case UniformD:
                fprintf(pFP, "MP.delayElE = 'UniformD';\n"); c++;
                FPRINT_FLOAT(pFP, MP.d_min); c++;
                FPRINT_FLOAT(pFP, MP.d_max); c++;
                break;
            case GaussD:
                fprintf(pFP, "MP.delayElE = 'GaussD';\n"); c++;
                FPRINT_FLOAT(pFP, MP.d_mean); c++;
                FPRINT_FLOAT(pFP, MP.d_sd); c++;
                break;
            case SOMD:
                fprintf(pFP, "MP.delayElE = 'SOM';\n"); c++; // SOMD would break Matlab code
                FPRINT_FLOAT(pFP, MP.condSpeed); c++;
                FPRINT_FLOAT(pFP, MP.maxDelay); c++;
            default:
                break;
        }
        
        switch (mp->delayEI) //(mp->axonDelay)
        {
            case MinD:
                fprintf(pFP, "MP.delayEI = 'MinD';\n"); c++;
                //FPRINT_FLOAT(pFile, MP.d_const); c++;
                break;
            case ConstD:
                fprintf(pFP, "MP.delayEI = 'ConstD';\n"); c++;
                FPRINT_FLOAT(pFP, MP.d_const); c++;
                break;
            case UniformD:
                fprintf(pFP, "MP.delayEI = 'UniformD';\n"); c++;
                FPRINT_FLOAT(pFP, MP.d_min); c++;
                FPRINT_FLOAT(pFP, MP.d_max); c++;
                break;
            case GaussD:
                fprintf(pFP, "MP.delayEI = 'GaussD';\n"); c++;
                FPRINT_FLOAT(pFP, MP.d_mean); c++;
                FPRINT_FLOAT(pFP, MP.d_sd); c++;
                break;
            case SOMD:
                fprintf(pFP, "MP.delayEI = 'SOM';\n"); c++; // SOMD would break Matlab code
                FPRINT_FLOAT(pFP, MP.condSpeed); c++;
                FPRINT_FLOAT(pFP, MP.maxDelay); c++;
            default:
                break;
        }
    }
    else
        fprintf(pFP, "%% Axon delays set to minimum (1 tstep)\n");
        
	FPRINT_FLOAT(pFP, MP.spatialScale); c++; // Needed for plotting connectivity
	
	FPRINT_INT(pFP, MP.SOM); c++;
	if (MP.SOM)
	{
		FPRINT_INT(pFP, MP.SOMinput); c++;
		FPRINT_FLOAT(pFP, MP.SOMsigE); c++;
		FPRINT_FLOAT(pFP, MP.SOMsigI); c++;
		FPRINT_FLOAT(pFP, MP.SOMclip); c++;
	}
	FPRINT_INT(pFP, MP.trainElE); c++;
	//FPRINT_INT(pFP, MP.initElE); c++;
	switch (mp->initElE)
	{
		case Constant:
			fprintf(pFP, "MP.initElE = 'Constant';\n"); c++;
			FPRINT_FLOAT(pFP, MP.iElE); c++;
			break;
		case Uniform:
			fprintf(pFP, "MP.initElE = 'Uniform';\n"); c++;
			break;
		case Gaussian:
			fprintf(pFP, "MP.initElE = 'Gaussian';\n"); c++;
			FPRINT_FLOAT(pFP, MP.SOMsigE); c++;
			break;
		case SOM:
			fprintf(pFP, "MP.initElE = 'SOM';\n"); c++;
			FPRINT_FLOAT(pFP, MP.SOMsigE); c++;
			break;
		default:
			break;
	}
	
	//FPRINT_INT(pFP, MP.initEfE); c++;
	//FPRINT_FLOAT(pFP, MP.iEfE); c++;
	switch (mp->initEfE)
	{
		case Constant:
			fprintf(pFP, "MP.initEfE = 'Constant';\n"); c++;
			FPRINT_FLOAT(pFP, MP.iEfE); c++;
			break;
		case Uniform:
			fprintf(pFP, "MP.initEfE = 'Uniform';\n"); c++;
			break;
		case Gaussian:
			fprintf(pFP, "MP.initEfE = 'Gaussian';\n"); c++;
			FPRINT_FLOAT(pFP, MP.SOMsigE); c++;
			break;
		case SOM:
			fprintf(pFP, "MP.initEfE = 'SOM';\n"); c++;
			FPRINT_FLOAT(pFP, MP.SOMsigE); c++;
			break;
		default:
			break;
	}
	
	printIntArray(pFP, "MP.vSquare", mp->vSquare, mp->nLayers); c++;
	for (l=0; l<MP.nLayers; l++)
	{
		fprintf(pFP, "%%%% Layer %d Dimensions %%%%\n", l);
		fprintf(pFP, "%% MP.layDim[%d].nRows = %d;\n",l,MP.layDim[l].nRows); c++;
		fprintf(pFP, "%% MP.layDim[%d].nCols = %d;\n",l,MP.layDim[l].nCols); c++;
		fprintf(pFP, "%% MP.layDim[%d].nFilt = %d;\n",l,MP.layDim[l].nFilt); c++;
	}
	
	/* Cell bodies */
	fprintf(pFP, "\n%%%% Cell Body Parameters %%%%\n");
	FPRINT_FLOAT(pFP, MP.capE); c++;
	FPRINT_FLOAT(pFP, MP.capI); c++;
	FPRINT_FLOAT(pFP, MP.gLeakE); c++;
	FPRINT_FLOAT(pFP, MP.gLeakI); c++;
	FPRINT_FLOAT(pFP, MP.VrestE); c++;
	FPRINT_FLOAT(pFP, MP.VrestI); c++;
	FPRINT_FLOAT(pFP, MP.VhyperE); c++;
	FPRINT_FLOAT(pFP, MP.VhyperI); c++;
	FPRINT_FLOAT(pFP, MP.ThreshE); c++;
	FPRINT_FLOAT(pFP, MP.ThreshI); c++;
	FPRINT_FLOAT(pFP, MP.VrevE); c++;
	FPRINT_FLOAT(pFP, MP.VrevI); c++;
	FPRINT_FLOAT(pFP, MP.refract); c++;
	FPRINT_INT(pFP, MP.adaptation); c++;
	if (MP.adaptation)
	{
		FPRINT_FLOAT(pFP, MP.alphaCa); c++;
		FPRINT_FLOAT(pFP, MP.tauCa); c++;
		FPRINT_FLOAT(pFP, MP.gAHP); c++;
		FPRINT_FLOAT(pFP, MP.VK); c++;
	}
	
	/* Synapses */
	fprintf(pFP, "\n%%%% Synapse Parameters %%%%\n");
	FPRINT_FLOAT(pFP, MP.alphaC); c++;
	FPRINT_FLOAT(pFP, MP.tauC); c++;
	FPRINT_FLOAT(pFP, MP.alphaD); c++;
	FPRINT_FLOAT(pFP, MP.tauD); c++;
	FPRINT_FLOAT(pFP, MP.learnR); c++;
	FPRINT_FLOAT(pFP, MP.DgEfE); c++;
	//FPRINT_FLOAT(pFile, MP.tauEfE); c++;
	FPRINT_FLOAT(pFP, MP.DgElE); c++;
	//FPRINT_FLOAT(pFile, MP.tauElE); c++;
	FPRINT_FLOAT(pFP, MP.tauEE); c++;
	FPRINT_FLOAT(pFP, MP.DgIE); c++;
	FPRINT_FLOAT(pFP, MP.tauIE); c++;
	FPRINT_FLOAT(pFP, MP.DgEI); c++;
	FPRINT_FLOAT(pFP, MP.tauEI); c++;
	FPRINT_FLOAT(pFP, MP.DgII); c++;
	FPRINT_FLOAT(pFP, MP.tauII); c++;
	FPRINT_FLOAT(pFP, MP.gMax); c++;

	fclose(pFP);
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
