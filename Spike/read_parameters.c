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
	char * s, buff[BUFFER];
	FILE * fp;
	
	fp = myfopen(paramfile, "r");
	while ((s = fgets(buff, sizeof(buff), fp)) != NULL) 	/* Read next line */
		count += parse_string(params, buff);
	fclose(fp);
	
	/* Additional calculations */
	float transP = (params->transP_Train >= params->transP_Test ? params->transP_Train : params->transP_Test);
	params->TotalTime = params->nStimuli * params->nTransPS * transP;
	params->TotalMS = ceil(params->TotalTime * 1000);
	params->TotalTS = ceil(params->TotalTime/DT);
	params->TSperMS = ceil(1/(DT*1000));
	params->spkBuffer = ceil(params->TotalTime / params->refract);
	if (floor(params->a * params->nExcit) < 1)
		params->nFiringNeurons = floor(params->nExcit/(params->nStimuli * params->nTransPS));
	else
		params->nFiringNeurons = floor(params->nExcit * params->a);
	//params->nFiringNeurons = (params->a == 0.0) ? floor(params->nExcit/(params->nStimuli * params->nTransPS)) : floor(params->nExcit * params->a);

	/*if (params->interleaveTrans)
		params->randStimOrder = false;*/
	
	return count;
}

char * trim(char * s)
{ // Remove leading and trailing whitespace
	/* Initialize start, end pointers */
	char *s1 = s, *s2 = &s[strlen (s) - 1];
	
	/* Trim and delimit right side */
	while ( (isspace (*s2) || *s2 == ';') && (s2 >= s1) )
		s2--;
	*(s2+1) = '\0';
	
	/* Trim left side */
	while ( (isspace(*s1)) && (s1 < s2) )
		s1++;
	
	/* Copy finished string */
	strcpy(s, s1);
	return s;
}

int parse_string(PARAMS * params, char * string)
{
	static char name[MAXLEN], value[MAXLEN];
	char * s;
	
	if (string[0] == '\n' || string[0] == '#' || string[0] == '%') 	/* Skip blank lines and comments */
		return 0; //continue;
	/* Parse name/value pair from line */
	s = strtok(string, "=");
	if (s==NULL)
	{
		fprintf(stderr, "Error! String: \"%s\" has no name!\n", s);
		return 0; // No parameter name on this line
	}
	else
		strncpy(name, s, MAXLEN);
	s = strtok(NULL, "=");
	if (s==NULL)
	{
		fprintf(stderr, "Error! Name: \"%s\" has no associated value!\n", name);
		return 0;
	}
	else
		strncpy(value, s, MAXLEN);
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
	else if (strcmp(name, "nRecordsPL")==0)
		params->nRecordsPL = atoi(value);
	
	/* Stimuli */
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
	else if (strcmp(name, "a")==0)
		params->a = atof(value);
	else if (strcmp(name, "nStimuli")==0)
		params->nStimuli = atoi(value);
	else if (strcmp(name, "nTransPS")==0)
		params->nTransPS = atoi(value);
	else if (strcmp(name, "transP_Train")==0)
		params->transP_Train = atof(value);
	else if (strcmp(name, "transP_Test")==0)
		params->transP_Test = atof(value);
	else if (strcmp(name, "shift")==0)
		params->shift = atoi(value);
	
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
	else if (strcmp(name, "nSynEfE")==0)
		params->nSynEfE = atoi(value);
	else if (strcmp(name, "nSynElE")==0)
		params->nSynElE = atoi(value);
	else if (strcmp(name, "nSynIE")==0)
		params->nSynIE = atoi(value);
	else if (strcmp(name, "nInhib")==0)
		params->nInhib = atoi(value);
	else if (strcmp(name, "nSynEI")==0)
		params->nSynEI = atoi(value);
	else if (strcmp(name, "nSynII")==0)
		params->nSynII = atoi(value);
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
