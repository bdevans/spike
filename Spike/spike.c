/*
 *  spike.c
 *  Spike
 *
 *  Created by Ben Evans on 6/19/08.
 *  Copyright 2008 University of Oxford. All rights reserved.
 *
 */

#include "spike.h"

int spike(PARAMS * mp)
{
	/*** Declare variables ***/
	int error = 0;
	STIMULI * stim = NULL;
	STIMULI * gStim = NULL;
	
	//RECORD *RECSP = RECS;
	//RECORD *ptr = &RECS[0][0];
	//RECORD *r_ptr = RECS;
	
	/*** Declare file pointers ***/
	FILE * stimuli_FP = NULL;
	
	
	/*************** Calculate RAM requirements ****************/
	
	calcMemory(mp);
	
	printf("Ventral Visual Stream Spiking Neural Network Simulation starting...\n");
	
	
	/*************** Build & Initialize Network ****************/
	
	printf("\tBuilding the network...");
	
	n_E = allocn(mp->nLayers, mp->vExcit, EXCIT); // Create 2D array of Excitatory neuron structures
	n_I = allocn(mp->nLayers, mp->vInhib, INHIB); // Create 2D array of Inhibatory neuron structures
	if (mp->SOM) //(mp->initElE == SOM || mp->initEfE == SOM || mp->axonDelay == SOMD)
		distE = getLowTriF(mp->nLayers, mp->vExcit, 0.0);
	calcConnectivity(mp->probConnect);	// Calculate network connectivity

	printf("\tBuilding complete!\n");
	
	printf("\tInitialising the network...");
	
	if (mp->nRecords)
		setRecords(mp, n_E, mSeed);
	setWeights(mp, n_E, n_I, ""); //regime = Learning;	// 0: Testing (No STDP); 1: Training (STDP);
	initNetwork(Hard); //(regime);		// Initialise parameters
	
	
	/*************** Load stimuli ****************/
	
	printf("\tNetwork initialised!\n");
	
	printf("\tCreating stimuli structures...");
	if (mp->priorPhases) // Load stimuli for ElE training
	{
		gStim = myalloc(sizeof(*gStim));
		gStim->trn_stimuli = gStim->tst_stimuli = NULL;
		gStim->stimShuffle = NULL; 
		gStim->transShuffle = NULL;
		gStim->groups = NULL;
		gStim->trnImages = gStim->tstImages = NULL;
	}
	
	stim = myalloc(sizeof(*stim));
	stim->trn_stimuli = stim->tst_stimuli = NULL;
	stim->stimShuffle = NULL; 
	stim->transShuffle = NULL;
	stim->groups = NULL;
	//stim->trnGrpStim = stim->tstGrpStim = NULL;
	stim->trnImages = stim->tstImages = NULL;
	
	/*stim->nStim = mp->nStimuli;
	stim->nTrans = mp->nTransPS; 
	stim->nTestStim = mp->nTestStimuli; // CHECK
	stim->nTestTrans = mp->nTestTransPS; // CHECK*/
	
	printf("\tCreated!\n");
	
	if (mp->useFilteredImages)
	{
		printf("\tLoading the images...");
		stim->trnImages = get_7D_farray(mp->nStimuli, mp->nTransPS, \
										mp->nScales, mp->nOrients, mp->nPhases, mp->nRows, mp->nCols, 0.0);
		
		error = loadImages(stim, mp); // if (!error)...
		printf("\t{S%d,T%d}", mp->nStimuli, mp->nTransPS);
		if (stim->newTestSet)
			printf(" Test: {S%d,T%d}", mp->nTestStimuli, mp->nTestTransPS);
		else //if (!mp->newTestSet) // Move inside loadImages?
		{
			stim->tstImages = stim->trnImages;
			mp->nTestStimuli = mp->nStimuli;
			mp->nTestTransPS = mp->nTransPS;
		}
		
		stim->nStim = mp->nStimuli;
		stim->nTrans = mp->nTransPS;
		stim->nTestStim = mp->nTestStimuli;
		stim->nTestTrans = mp->nTestTransPS;
		
		if (!error)
			printf("\tImages loaded!\n");
		else
			exit_error("spike", "Error loading Images");
	}
	else if (mp->stimGroups)
	{
		if (mp->loadStimuli)
		{
            if (mp->priorPhases) // Load training & testing stimuli for prior phases
            {
                printf("\tLoading PP grouped stimuli...");
                
                
                /*if (!PPSTFILE) // Use STIMULIFILE as the default file name if one has not been passed
                 {
                 printf("\t\"%s\" (default)",STIMULIFILE);
                 int slen = strlen(STIMULIFILE);
                 PPSTFILE = myalloc((slen+1)*sizeof(char));
                 strncpy(PPSTFILE, STIMULIFILE, slen);
                 PPSTFILE[slen] = '\0';
                 }
                 else
                 printf("\t\"%s\"",ELESTFILE);*/
                
                printf("\t\"%s\"",PPSTFILE);
                gStim->nStim = gStim->nTestStim = 0;
                gStim->nTrans = gStim->nTestTrans = 0;
                loadGroups(gStim, mp, PPSTFILE); // Pass filename here *****************************
                printf("\tPP stimuli loaded!\n");
            }
            
            printf("\tLoading grouped stimuli..."); 
            /*if (!STFILE) // Use STIMULIFILE as the default file name if one has not been passed
             {
             printf("\t\"%s\" (default)",STIMULIFILE);
             int slen = strlen(STIMULIFILE);
             STFILE = myalloc((slen+1)*sizeof(char));
             strncpy(STFILE, STIMULIFILE, slen);
             STFILE[slen] = '\0';
             }
             else
             printf("\t\"%s\"",STFILE);*/
            
            printf("\t\"%s\"",STFILE); 
            loadGroups(stim, mp, STFILE);
            assert(stim->nStim == mp->nStimuli);
            assert(stim->nTrans == mp->nTransPS);
            printf("\tGroups loaded!\n"); 
            
            // Print prototypes for neuron labelling
            stimuli_FP = myfopen("prototypes.stm", "w");
            print_iarray(stimuli_FP, stim->groups, mp->nGroups, mp->sInputs);
            fclose(stimuli_FP);
            
            if (stim->newTestSet)
            {
                assert(stim->nTestStim == mp->nTestStimuli);
                assert(stim->nTestTrans == mp->nTestTransPS);
            }
            else
            {
                mp->nTestStimuli = mp->nStimuli;
                mp->nTestTransPS = mp->nTransPS;
                stim->tst_stimuli = stim->trn_stimuli;
                stim->nTestStim = mp->nTestStimuli;
                stim->nTestTrans = mp->nTestTransPS;
            }
		}
		else // Generate groups of stimuli
		{
			if (mp->priorPhases)
			{
				printf("\tGenerating PP grouped stimuli...");
				printf("\t\"%s\"",PPSTFILE);
				gStim->nStim = gStim->nTestStim = 0;
				gStim->nTrans = gStim->nTestTrans = 0;
				genGroups(gStim, mp); // Pass filename here *****************************
				printGroups(gStim, mp, PPSTFILE);
				printf("\tPP stimuli saved!\n");
			}
			
			printf("\tGenerating grouped stimuli..."); 
			printf("\t\"%s\"",STFILE); 
			genGroups(stim, mp);
			printGroups(stim, mp, STFILE);
			assert(stim->nStim == mp->nStimuli);
			assert(stim->nTrans == mp->nTransPS);
			printf("\tGroups saved!\n"); 
			
			// Print prototypes for neuron labelling
			stimuli_FP = myfopen("prototypes.stm", "w");
			print_iarray(stimuli_FP, stim->groups, mp->nGroups, mp->sInputs);
			fclose(stimuli_FP);
			genGroups(stim, mp);
			printGroups(stim, mp, STFILE);
			
			if (stim->newTestSet)
			{
				assert(stim->nTestStim == mp->nTestStimuli);
				assert(stim->nTestTrans == mp->nTestTransPS);
			}
			else
			{
				mp->nTestStimuli = mp->nStimuli;
				mp->nTestTransPS = mp->nTransPS;
				stim->tst_stimuli = stim->trn_stimuli;
				stim->nTestStim = mp->nTestStimuli;
				stim->nTestTrans = mp->nTestTransPS;
			}
		}

	}
	else
	{
		printf("\tCreating the stimuli...");
		gen_stimuli(mp->localRep, stim, mp);			// Generate Patterns
		printStimuli(stim, mp);
		printf("\tStimuli saved!\n");
	}
	
	
	// Generate shuffles 
	genShuffles(stim, mp);
	if (mp->priorPhases)
		genShuffles(gStim, mp);

	// Set up variables for estimating percentage completion
	SIM.tally = 0;
	SIM.ptTS = (mp->pretrain) ? mp->nTestStimuli * mp->nTestTransPS * mp->transP_Test * ceil(1/mp->DT): 0;
	SIM.trainTS = (mp->train) ? mp->loops * mp->nStimuli * mp->nTransPS * mp->transP_Train * ceil(1/mp->DT): 0;
	SIM.testTS = mp->nTestStimuli * mp->nTestTransPS * mp->transP_Test * ceil(1/mp->DT);
	if (mp->priorPhases) // Recalculate (as EfE phase may have differnt numbers of stimuli and trans...
	{
		SIM.trainTS += mp->loops * gStim->nStim * gStim->nTrans * mp->transP_Train * ceil(1/mp->DT); // REVISE
		SIM.testTS += 2 * gStim->nTestStim * gStim->nTestTrans * mp->transP_Test * ceil(1/mp->DT);
	}
	SIM.totTS = SIM.ptTS + SIM.trainTS + SIM.testTS;
	assert(SIM.totTS <= pow(2, 8*sizeof(tstep)-1)-1); // assumes tstep will be unsigned otherwise pow(2, 8*sizeof(tstep))-1
	
	
	/*************** Simulation Phases ****************/
	
	if (mp->priorPhases && !mp->loadWeights) // Lateral weight training
	{
		mp->trainElE = true;
		if (mp->isolateEfE)
			mp->trainEfE = false;
		
		mp->nRecords = 0;	// Disable records (or recalculate buffers)
		
		if (mp->pretrain)
		{
			printf("\tNow beginning PP PreTraining phase...\n"); 
			simulatePhase(Testing, "pt_PP_", gStim);
			printf("\tPP PreTraining complete!\n");
		}
		
		// Incorporate additional time steps into SIM.totTS
		printf("\tNow beginning PP Training phase...\n"); 
		simulatePhase(Training, "_PP_", gStim); // Train ElE <only> with exemplars individually
		printf("\tPP Training complete!\n");
		
		if (mp->pretrain) // Test to confirm desynchronised representations for novel exemplars
		{
			printf("\tNow beginning PP Testing phase...\n"); 
			simulatePhase(Testing, "_PP_", gStim); // Test with novel exemplars combined
			printf("\tPP Testing complete!\n");
		}
		
		printf("\tFixing PP weights...\t");
		mp->trainElE = false; // May be better to let ElE adjust during EfE training...
		printf("PP weights fixed!\n");
		
		if (mp->isolateEfE)
			mp->trainEfE = true;
		
		if (mp->vRecords)
		{
			int l=0;
			for (l=0; l<mp->nLayers; l++) // Reset Records
				mp->nRecords += mp->vRecords[l];
		}
	}
	
	if (mp->pretrain)
	{
		printf("\tNow beginning PreTraining phase...\n"); 
		simulatePhase(Testing, "pt", stim);
		printf("\tPreTraining complete!\n");
	}
		
	if (mp->train)
	{
		printf("\tNow beginning Training phase...\n"); 
		simulatePhase(Training, "", stim);
		printf("\tTraining complete!\n");
	}
		
	printf("\tNow beginning Testing phase...\n"); 
	simulatePhase(Testing, "", stim);
	printf("\tTesting complete!\n");
	
	
	/*************** Deallocate Memory ****************/
	printf("\tDeallocating memory...");
	
	unallocn(n_E, mp->nLayers, mp->vExcit);
	unallocn(n_I, mp->nLayers, mp->vInhib);
	if (mp->SOM) // (mp->initElE == SOM || mp->initEfE == SOM || mp->axonDelay == SOMD)
		freeTriF(distE, mp->nLayers);
	
	if (mp->randStimOrder)
		free_2D_iarray(stim->stimShuffle);//, mp->loops);
	if (mp->randTransOrder) // || mp->randTransDirection)
		free_3D_iarray(stim->transShuffle, mp->loops);//, mp->nStimuli);
	if (mp->useFilteredImages)
	{
		free_7D_farray(stim->trnImages, mp->nStimuli, mp->nTransPS, mp->nScales, mp->nOrients, mp->nPhases);
		if (stim->newTestSet)
			free_7D_farray(stim->tstImages, mp->nTestStimuli, mp->nTestTransPS, mp->nScales, mp->nOrients, mp->nPhases);
	}
	else
	{
		if (mp->stimGroups)
			myfree(stim->groups);
		if (stim->newTestSet)
			free_3D_farray(stim->trn_stimuli, stim->nStim);//, mp->nTransPS);
		free_3D_farray(stim->tst_stimuli, stim->nTestStim);//, mp->nTransPS);
	}
	myfree(stim); // *** Free new image arrays too

	if (mp->priorPhases) // Free second set of stimuli
	{
		if (mp->randStimOrder)
			free_2D_iarray(gStim->stimShuffle);
		if (mp->randTransOrder)
			free_3D_iarray(gStim->transShuffle, mp->loops);
		free_3D_farray(gStim->trn_stimuli, gStim->nStim);
		if (gStim->newTestSet)
			free_3D_farray(gStim->tst_stimuli, gStim->nTestStim);
		if (mp->stimGroups) // Always true for PP?
			myfree(gStim->groups);
		myfree(gStim);
	}
	
	printf("\tMemory Deallocated!\n");

	return 0;
}


void calcMemory(PARAMS * mp)
{
    fprintf(stdout, "--------------------------------------------------------------------------------\n");
#if DEBUG > 1
	fprintf(stdout, "Variable type:\tNEURON\tAXON  \t*     \tfloat \ttstep \tint\n");
	fprintf(stdout, "Size (bytes): \t%-6lu\t%-6lu\t%-6lu\t%-6lu\t%-6lu\t%-6lu\t\n",\
			sizeof(NEURON),sizeof(AXON),sizeof(int*),sizeof(float),sizeof(tstep),sizeof(int));
#endif
	//size_t size = sizeof();
	float EsynE = 0;
    float EsynEfE = 0;
  	float EsynElE = 0;
    float EsynIE = 0;
	float EsynI = 0;
    float EsynEI = 0;
    float EsynII = 0;
	float memE = 0.0;
	float memI = 0.0;
	float memMisc = 0.0;
	float memTrain = 0.0;
	float memTest = 0.0;
	float Tmem = 0.0;
    float avqEfE = 1;
    float avqElE = 1;
    float avqEI = 1;
    int mult = 0;
    float base = 0.0;
	
	int MB = 1024*1024;
	int l=0;
#if DEBUG > 1
	fprintf(stdout, "\nLayer\tExcit \tInhib \tMem (MB)\n");
#endif
	for (l=0; l<mp->nLayers; l++)
	{
		Tmem += (mp->vExcit[l]+mp->vInhib[l])*(sizeof(NEURON)+(mp->spkBuffer*(float)sizeof(tstep)))/MB;
#if DEBUG > 1
		fprintf(stdout,"%-6d\t%-6d\t%-6d\t%-6.2f\n",l,mp->vExcit[l],mp->vInhib[l],\
				(mp->vExcit[l]+mp->vInhib[l])*(sizeof(NEURON)+(mp->spkBuffer*(float)sizeof(tstep)))/MB);
#endif
	}
	
    if (mp->axonDelay)
    {
        avqEfE = (mp->delayEfE) ? meanQueue(mp, mp->delayEfE) : 1;
        avqElE = (mp->delayElE) ? meanQueue(mp, mp->delayElE) : 1;
        avqEI = (mp->delayEI) ? meanQueue(mp, mp->delayEI) : 1;
    }

	
#if DEBUG > 1
	//fprintf(stderr,"\nPresynaptic connection probabilites for Excitatory postsynaptic cells\n");
	fprintf(stdout,"\n\tp(EfE)\tp(ElE)\tp(IE) \tE[syn] \tp(EI) \tp(II) \tE[syn]\tMem (MB)\n");
#endif
	// Assumes minimum delay model i.e. 1 spike bin per axon for all except EfE, ElE & EI
	for (l=0; l<mp->nLayers; l++)
	{
        EsynEfE = ((l>0) ? (mp->vExcit[l-1]*mp->pCnxEfE[l]) : 0) * mp->vExcit[l];
        EsynElE = pow(mp->vExcit[l],2) * mp->pCnxElE[l];
        EsynIE = mp->vInhib[l] * mp->pCnxIE[l] * mp->vExcit[l];
		EsynE = EsynEfE + EsynElE + EsynIE; 
        memE = EsynE * (sizeof(AXON) + (sizeof(NEURON*) * 3))/MB;
        memE += ((avqEfE*EsynEfE) + (avqElE*EsynElE) + EsynIE) * sizeof(tstep)/MB;

		EsynEI = mp->vExcit[l] * mp->pCnxEI[l] * mp->vInhib[l];
        EsynII = pow(mp->vInhib[l], 2) * mp->pCnxII[l];
        EsynI = EsynEI + EsynII;
        memI = EsynI * (sizeof(AXON) + (sizeof(NEURON*) * 3))/MB;
        memI += (avqEI*EsynEI + EsynII) * sizeof(tstep)/MB;
#if DEBUG > 1
		fprintf(stdout,"L%d:\t%-6.3f\t%-6.3f\t%-6.3f\t%-6.2G\t%-6.3f\t%-6.3f\t%-6.2G\t%-6.2f\n", \
				l,mp->pCnxEfE[l],mp->pCnxElE[l],mp->pCnxIE[l],(float)EsynE, \
				mp->pCnxEI[l],mp->pCnxII[l],(float)EsynI,memE+memI);
#endif
		Tmem += (memE + memI);
	}
	
#if DEBUG > 1
	fprintf(stdout, "\nStimuli:\tStruct.\tTrain \tTest  \tRecords\n");
#endif
	if (mp->SOM) // Triangle of Distances // (mp->initElE == SOM || mp->initEfE == SOM || mp->axonDelay == SOMD)
	{
		for (l=0; l<mp->nLayers; l++)
			memMisc += mp->vExcit[l]*(mp->vExcit[l]+1)/2; // Gauss' method 
		memMisc *= sizeof(float)/MB;
	}
	
	memE = 0.0;
	if (mp->nRecords)
	{
        mult = 1; // V
        base = (mp->RecordMS+1) * sizeof(float);
        mult += ((mp->adaptation) ? 1 : 0) + ((mp->train) ? 1 : 0); // cCa & D
        for (l=0; l<mp->nLayers; l++)
            memE += mp->vRecords[l] * base;
        mult = (mp->train && mp->trainElE) ? 3 : 1; // Lateral g (& Dg, C)
        for (l=0; l<mp->nLayers; l++)
            memE += mp->vRecords[l] * mp->vExcit[l] * mp->pCnxElE[l] * mult * base;
        for (l=(mp->inputInhib ? 0 : 1); l<mp->nLayers; l++) // Sigma g_I
            memE += mp->vRecords[l] * base;
        for (l=1; l<mp->nLayers; l++)
			memE += mp->vRecords[l]*((mp->train)?3:1)*(mp->vExcit[l-1]*mp->pCnxEfE[l])*base; // g (& Dg, C)
		for (l=0; l<mp->nLayers; l++)
            memE += mp->vRecords[l] * sizeof(RECORD); // Record structures
        memE /= MB;
	}
	
	memMisc += (float)(sizeof(PARAMS) + \
					   ((mp->randStimOrder)?mp->loops*mp->nStimuli*sizeof(int):0) + \
					   ((mp->randTransOrder)?mp->loops*mp->nStimuli*mp->nTransPS*sizeof(int):0))/MB;
	
	memTrain = (float)(mp->sInputs*mp->nStimuli*mp->nTransPS*sizeof(float))/MB;
	memTest = (float)(mp->sInputs*mp->nTestStimuli*mp->nTestTransPS*sizeof(float))/MB;
#if DEBUG > 1
	fprintf(stdout, "Size (MB)\t%-6.2f\t%-6.2f\t%-6.2f\t%-6.2f\n\n",memMisc,memTrain,memTest,memE);
#endif
	Tmem += (memMisc + memTrain + memTest + memE);	
	fprintf(stdout, "Total memory requirements (approx.):\t%-6.3G MB\n",Tmem);
	fprintf(stdout, "--------------------------------------------------------------------------------\n"); // Replace with horizontal line cmd?
}


float meanQueue(PARAMS * mp, DELAY synClass)
{
    float avq = 1;
    if (mp->axonDelay)
    {
        switch (synClass) 
        {
            case MinD:
                avq = 1;
                break;
            case ConstD:
                avq = round(mp->d_const/mp->refract);
                break;
            case UniformD:
                avq = (mp->d_min + ((mp->d_max - mp->d_min) / 2))/mp->refract;
                break;
            case GaussD:
                avq = mp->d_mean / mp->refract;
                break;
            case SOMD:
                avq = mp->maxDelay / (2.0 * mp->refract); // Reasonable for 1D layer
                break;
            default:
                exit_error("create_axons", "Unknown axonal delay model!");
                break;
        }
        //avq = (avq < 1) ? 1 : avq;
    }

    return avq; //(avq < 1) ? 1 : ceil(avq);
}

NEURON ** allocn (int nLays, int * vNeurons, NTYPE type)
{
	int l, n;
	int rowLen = 0;
	int colLen = 0;
	float sp_x = 0.0;
	float sp_y = 0.0;
	NEURON *space;
	NEURON **narray;
	// sizeof(*array) = sizeof(int **)
	// sizeof(**array) = sizeof(int *)
	// sizeof(***array) = sizeof(int)
	
	int totNeurons = 0;
	for (l=0; l<nLays; l++)
		totNeurons += vNeurons[l];
	
	space = myalloc(totNeurons * sizeof(*space)); 	/*** Ensures array is contiguous in memory ***/
	narray = myalloc(nLays * sizeof(*narray));

	totNeurons = 0;
	for (l=0; l<nLays; l++)
	{
		narray[l] = space + totNeurons; //(l * nneurons); 
		totNeurons += vNeurons[l];
	}
		
	/* Allocate spike time bins and initialise connectors */
#pragma omp parallel default(shared) private(l,n,rowLen,colLen,sp_x,sp_y)
	{
	for (l=0; l<nLays; l++)
	{
		rowLen = mp->layDim[l].nCols;	//rowSize = (mp->SOM && !(l==0 && mp->SOMinput)) ? mp->layDim[l].nCols : vNeurons[l]; 
		colLen = mp->layDim[l].nRows;
		sp_x = 1.0/mp->layDim[l].nCols; // x spacing
		sp_y = 1.0/mp->layDim[l].nRows; // y spacing
#if DEBUG > 3
#pragma omp master
		printf("\nLayer %d: nRows = %d; nCols = %d; sp_x = %f; sp_y = %f\n",l,mp->layDim[l].nRows,mp->layDim[l].nCols,sp_x,sp_y);
#pragma omp barrier
#endif
#pragma omp for
		for (n=0; n<vNeurons[l]; n++)
		{
			narray[l][n].spkbin = 0;
			//if (mp->useFilteredImages && type==EXCIT && l==0)
			//	narray[l][n].spikeTimes = myalloc(mp->inpSpkBuff * sizeof(narray[l][n].spikeTimes[0]));
			//else
			narray[l][n].spikeTimes = myalloc(mp->spkBuffer * sizeof(narray[l][n].spikeTimes[0]));
			narray[l][n].type = type; //(type==EXCIT) ? EXCIT : INHIB; 
			narray[l][n].nFAff_E = 0;
			narray[l][n].FAffs_E = NULL; // **
			narray[l][n].lm1presyn_E = NULL; // **
			narray[l][n].nLAff_E = 0;
			narray[l][n].LAffs_E = NULL; // **
			narray[l][n].lm0presyn_E = NULL; // **
			narray[l][n].nLAff_I = 0;
			narray[l][n].LAffs_I = NULL; // **
			narray[l][n].lm0presyn_I = NULL; // **
			narray[l][n].n = n; // Linear index
			narray[l][n].l = l;
			// The 2D indexes apply to simple inputs (no hypercolumns) // Consider 1D/2D/3D layers and square/rectangular
			narray[l][n].row = n / rowLen; // implied floor()
			narray[l][n].col = n % rowLen; // n % b := b - (n * floor(b/n))
			narray[l][n].x = (sp_x/2 + (narray[l][n].col * sp_x)) * mp->spatialScale; // col: 0,1,...,rowLen-1
			narray[l][n].y = (sp_y/2 + (((colLen - 1) - narray[l][n].row) * sp_y)) * mp->spatialScale; // row indices and y coordinates run opposite
			// Multiply x & y coords by a spatial scaling factor so condSpeed can be biologically accurate?
#if DEBUG > 3
			printf("L%dN%d [R%d,C%d]:(%f,%f); ",l,n,narray[l][n].row,narray[l][n].col,narray[l][n].x,narray[l][n].y);
#endif
			narray[l][n].nFEff_E = 0;
			narray[l][n].FEffs_E = NULL;
			narray[l][n].lp1postsyn_E = NULL; // **
			narray[l][n].nLEff_E = 0;
			narray[l][n].LEffs_E = NULL;
			narray[l][n].lp0postsyn_E = NULL; // **
			narray[l][n].nLEff_I = 0;
			narray[l][n].LEffs_I = NULL;
			narray[l][n].lp0postsyn_I = NULL; // **
			narray[l][n].rec_flag = false;
			narray[l][n].rec = NULL;
		}
	}
	}
	// Include connectivity arrays?
    // Create vRecords before allocn, then here calcConnectivity & allocRecords?
	
	return narray;
}


int unallocn (NEURON ** narray, int nLays, int * vNeurons)
{
	int l, n, s;
	for (l=0; l<nLays; l++)
	{
		for (n=0; n<vNeurons[l]; n++)
		{
			myfree(narray[l][n].spikeTimes);
			for (s=0; s<narray[l][n].nFEff_E; s++)
				myfree(narray[l][n].FEffs_E[s].queue);
			myfree(narray[l][n].FEffs_E);
			myfree(narray[l][n].lp1postsyn_E);
			myfree(narray[l][n].FAffs_E);
			myfree(narray[l][n].lm1presyn_E);
			for (s=0; s<narray[l][n].nLEff_E; s++)
				myfree(narray[l][n].LEffs_E[s].queue);
			myfree(narray[l][n].LEffs_E);
			myfree(narray[l][n].lp0postsyn_E);
			myfree(narray[l][n].LAffs_E);
			myfree(narray[l][n].lm0presyn_E);
			for (s=0; s<narray[l][n].nLEff_I; s++)
				myfree(narray[l][n].LEffs_I[s].queue);
			myfree(narray[l][n].LEffs_I);
			myfree(narray[l][n].lp0postsyn_I);
			myfree(narray[l][n].LAffs_I);
			myfree(narray[l][n].lm0presyn_I);
			if (narray[l][n].rec_flag)
			{
				myfree(narray[l][n].rec->cellV);//, mp->loops);
				if (mp->adaptation)
					myfree(narray[l][n].rec->cellcCa);
                if (narray[l][n].nLAff_I)
                    myfree(narray[l][n].rec->LsigGI);
                if (narray[l][n].nFAff_E) //(l>0)
                    free_2D_farray(narray[l][n].rec->FSynG);
                if (narray[l][n].nLAff_E) // Synaptic variables are stored with the post-synaptic neuron
                    free_2D_farray(narray[l][n].rec->LSynG);
                if (mp->train)
                {
                    myfree(narray[l][n].rec->cellD);//, mp->loops);
                    if (narray[l][n].nFAff_E) //(l>0)  // Here the synaptic variables are stored with the post-synaptic neuron
                    {
                        free_2D_farray(narray[l][n].rec->FSynDG);
                        free_2D_farray(narray[l][n].rec->FSynC);
                    }
                    
                    if (mp->trainElE && narray[l][n].nLAff_E) //mp->train)
                    {
                        free_2D_farray(narray[l][n].rec->LSynDG);
                        free_2D_farray(narray[l][n].rec->LSynC);
                    }
                }

					
				/*if (l<mp->nWLayers)
				{
					free_3D_farray(narray[l][n].rec->SynC, mp->loops);//, narray[l][n].nFAff_E);
					free_3D_farray(narray[l][n].rec->SynG, mp->loops);//, narray[l][n].nFAff_E);
					free_3D_farray(narray[l][n].rec->SynDG, mp->loops);//, narray[l][n].nFAff_E);
				}*/
				myfree(narray[l][n].rec);
			}
		}
	}
	myfree(*narray);
	myfree(narray);
	return 0;
}


void calcConnectivity(bool probConnect)
{
	int l = 0; 
	int n = 0;
	int s = 0;
	int slen = 0;
	char filename[FNAMEBUFF];
	FILE * connections_FP;
	int i, j, nRows, nCols; //, d_h, d_w, h_min, w_min;
	i = j = nRows = nCols = 0; // d_h = d_w = h_min = w_min = 0;
	
	if (mp->SOM) //(mp->initElE == SOM || mp->initEfE == SOM || mp->axonDelay == SOMD)
	{
		/*if (mp->SOM == 2)
			phiScale = 1/(mp->SOMsigE*sqrt(2*M_PI)); */
		
		// Build distance matrix - could build a much smaller one based upon distE[l][abs(r1-r2)][abs(c1-c2)]
		for (l=0; l<mp->nLayers; l++)
		{
#if DEBUG > 3
			printf("\nLayer %d: nRows = %d; nCols = %d;\n",l,mp->layDim[l].nRows,mp->layDim[l].nCols);
#endif
			for (i=0; i<mp->vExcit[l]; i++)
			{
#if DEBUG > 3
				printf("L%dN%d [R%d,C%d]:(%f,%f); ",l,i,n_E[l][i].row,n_E[l][i].col,n_E[l][i].x,n_E[l][i].y);
#endif
				for (j=0; j<=i; j++)
				{
					distE[l][i][j] = calcDistance(&n_E[l][i], &n_E[l][j], mp->spatialScale);//mp->layDim);
#if DEBUG > 3
					printf("Dist: N%d[%d,%d] --> N%d[%d,%d] = %f\n", i, n_E[l][i].row, n_E[l][i].col, \
						   j, n_E[l][j].row, n_E[l][j].col, readLowTriF(distE, l, i, j));
#endif
					/*if (mp->SOM == 2) // Probabilistically connect
						probE[l][i][j] = probE[l][j][i] = phiScale * exp(-pow(distance,2)/(2*pow(mp->SOMsigE,2))); // mu = 0*/
				}
			}
		}
	}
	
	if (!probConnect)
		mp->probConnect = true;
	/*if (!probConnect) // Deprecated - remove mp->probConnect
	{
		calc_connectivity();
		return;
	}*/
	
	/* Wire afferents */
	if (mp->loadWeights)
	{
#if DEBUG > 1
		printf("\nWarning: Loading weights requires the same size network as the loaded simulation!");
#endif
		loadAfferents(""); // Optionally pass suffix string
	}
	else
	{
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vExcit[l]; n++)
			wireAfferents(&n_E[l][n], mp->pCnxEfE[l], mp->pCnxElE[l], mp->pCnxIE[l]);
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vInhib[l]; n++)
			wireAfferents(&n_I[l][n], 0.0, mp->pCnxEI[l], mp->pCnxII[l]);
	}
		
	/* Allocate efferents */
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vExcit[l]; n++)
			alloc_efferents(&n_E[l][n]);
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vInhib[l]; n++)
			alloc_efferents(&n_I[l][n]);
	
	/* Wire efferents */
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vExcit[l]; n++)
			wire_efferents(&n_E[l][n]);
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vInhib[l]; n++)
			wire_efferents(&n_I[l][n]);
	
	/*** Axonal delays ***/
/*#if DEBUG >	1
	switch (mp->axonDelay) 
	{
		case MinD:
			printf("\nSetting axonal delays to minimum\n");
			break;
			
		case ConstD:
			printf("\nSetting axonal delays to %f seconds\n",mp->d_const);
			break;
			
		case UniformD:
			printf("\nDrawing axonal delays from [%f, %f]\n",mp->d_min, mp->d_max);
			break;
			
		case GaussD:
			printf("\nDrawing axonal delays from N(%f,%f)\n",mp->d_mean, mp->d_sd);
			break;
			
		case SOMD:
			printf("\nSetting axonal delays (ElE) proportional to Euclidean distances\n");
			break;
			
		default:
			printf("\nError: Unknown axonal delay model!\n");
			break;
	}
#endif*/
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vExcit[l]; n++)
			create_axons(&n_E[l][n], mp);
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vInhib[l]; n++)
			create_axons(&n_I[l][n], mp);
	
	/* Given a particular postsynaptic neuron <E|I>, n, and its synapse, s, the presynaptic 
	 cell forming the synapse <E|I> is given by presyncnx_<E|I><E|I>[l][n][s]. */
	/* NB If connections from the previous and the same layer are required, can label each neuron from 0 
	 to (NEXCIT*NLAYERS)-1 (the neuron ID) with l=floor(NID/NEXCIT) and n=((NID+1)%NEXCIT)-1.	*/
	
	/*	For affNeurons_EfE[0][n][s] : layer 0 cells --> layer 1 cells (NWLAYERS)
	 For affNeurons_ElE[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS)
	 For affNeurons_IE[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS)
	 For affNeurons_EI[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS)
	 For affNeurons_II[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS) */
	/*
	int s=0;
	int tot = 0;
	for (l=0; l<mp->nWLayers; l++)
		for (n=0; n<mp->vExcit[l+1]; n++)
			tot += n_E[l+1][n].nFAff_E;
	int *space = myalloc(tot*sizeof(*space));
	tot = 0;
	affNeurons_EfE = myalloc(mp->nWLayers*sizeof(*affNeurons_EfE));
	for (l=0; l<mp->nWLayers; l++)
	{
		affNeurons_EfE[l] = myalloc(mp->vExcit[l+1]*sizeof(**affNeurons_EfE));
		for (n=0; n<mp->vExcit[l+1]; n++)
		{
			affNeurons_EfE[l][n] = space + tot;
			tot += n_E[l+1][n].nFAff_E;
		}
	}
	
	for (l=0; l<mp->nWLayers; l++)
		for (n=0; n<mp->vExcit[l+1]; n++)
			for (s=0; s<n_E[l+1][n].nFAff_E; s++)
				affNeurons_EfE[l][n][s] = (n_E[l+1][n].lm1presyn_E[s])->n;
	
	
	tot = 0;
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vExcit[l]; n++)
			tot = tot + n_E[l][n].nLAff_E;
	space = myalloc(tot*sizeof(*space));
	tot = 0;
	affNeurons_ElE = myalloc(mp->nLayers*sizeof(*affNeurons_ElE));
	for (l=0; l<mp->nLayers; l++)
	{
		affNeurons_ElE[l] = myalloc(mp->vExcit[l]*sizeof(**affNeurons_ElE));
		for (n=0; n<mp->vExcit[l]; n++)
		{
			affNeurons_ElE[l][n] = space + tot;
			tot += n_E[l][n].nLAff_E;
		}
	}
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vExcit[l]; n++)
			for (s=0; s<n_E[l][n].nLAff_E; s++)
				affNeurons_ElE[l][n][s] = (n_E[l][n].lm0presyn_E[s])->n;
	*/
	
	/************** FILE OUTPUT **************/
	if (mp->printConnections && !mp->loadWeights)
	{
		for (l=0; l<mp->nWLayers; l++) // Loop up to nWLayers
		{
			slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsEfE.dat", l+1);
			assert(slen < FNAMEBUFF);
			connections_FP = myfopen(filename, "w");
			for (n=0; n<mp->vExcit[l+1]; n++)
			{
				fprintf(connections_FP, "%d\t", n_E[l+1][n].nFAff_E); // Print number of EfE synapses first
				for (s=0; s<n_E[l+1][n].nFAff_E; s++)
					fprintf(connections_FP, "%d ", n_E[l+1][n].lm1presyn_E[s]->n);
				fprintf(connections_FP, "\n");
			}
			fclose(connections_FP);
		}
		
        if (mp->delayEfE) // Print EfE delays (anything other than minimum)
        {
            for (l=0; l<mp->nWLayers; l++) // Loop up to nWLayers
            {
                slen = snprintf(filename, FNAMEBUFF, "L%daffDelaysEfE.dat", l+1);
                assert(slen < FNAMEBUFF);
                connections_FP = myfopen(filename, "w");
                for (n=0; n<mp->vExcit[l+1]; n++)
                {
                    fprintf(connections_FP, "%d\t", n_E[l+1][n].nFAff_E); // Print number of EfE synapses first
                    for (s=0; s<n_E[l+1][n].nFAff_E; s++)
                        fprintf(connections_FP, "%d ", n_E[l+1][n].FAffs_E[s]->delay);
                    fprintf(connections_FP, "\n");
                }
                fclose(connections_FP);
            }
        }
        
		for (l=0; l<mp->nLayers; l++)
		{
			slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsElE.dat", l);
			assert(slen < FNAMEBUFF);
			connections_FP = myfopen(filename, "w");
			for (n=0; n<mp->vExcit[l]; n++)
			{
				fprintf(connections_FP, "%d\t", n_E[l][n].nLAff_E); // Print number of ElE synapses first
				for (s=0; s<n_E[l][n].nLAff_E; s++) // if (mp->pCnxElE[l] > EPS)
					fprintf(connections_FP, "%d ", n_E[l][n].lm0presyn_E[s]->n);
				fprintf(connections_FP, "\n");
			}
			fclose(connections_FP);
		}
		
		if (mp->SOM) //(mp->initElE == SOM || mp->initEfE == SOM || mp->axonDelay == SOMD) // Print ElE distance
		{
			for (l=0; l<mp->nLayers; l++)
			{
				slen = snprintf(filename, FNAMEBUFF, "L%ddistElE.dat", l);
				assert(slen < FNAMEBUFF);
				connections_FP = myfopen(filename, "w");
				for (i=0; i<mp->vExcit[l]; i++)
				{
					for (j=0; j<=i; j++)
						fprintf(connections_FP, "%f ", readLowTriF(distE, l, i, j));//distE[l][i][j]); 
					fprintf(connections_FP, "\n");
				}
				fclose(connections_FP);
			}
						
            if (mp->delayElE) // Print ElE delays
            {
                for (l=0; l<mp->nLayers; l++)
                {
                    slen = snprintf(filename, FNAMEBUFF, "L%daffDelaysElE.dat", l); 
                    assert(slen < FNAMEBUFF);
                    connections_FP = myfopen(filename, "w");
                    for (n=0; n<mp->vExcit[l]; n++)
                    {
                        fprintf(connections_FP, "%d\t", n_E[l][n].nLAff_E); // Print number of ElE synapses first (in case any have 0)
                        for (s=0; s<n_E[l][n].nLAff_E; s++)
                            fprintf(connections_FP, "%d ", n_E[l][n].LAffs_E[s]->delay);
                        fprintf(connections_FP, "\n");
                    }
                    fclose(connections_FP);
                }
            }
		}
		
		for (l=0; l<mp->nLayers; l++)
		{
			slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsEI.dat", l);
			assert(slen < FNAMEBUFF);
			connections_FP = myfopen(filename, "w");
			for (n=0; n<mp->vInhib[l]; n++)
			{
				fprintf(connections_FP, "%d\t", n_I[l][n].nLAff_E); // Print number of EI synapses first
				for (s=0; s<n_I[l][n].nLAff_E; s++) // if (mp->pCnxEI[l] > EPS)
					fprintf(connections_FP, "%d ", n_I[l][n].lm0presyn_E[s]->n);
				fprintf(connections_FP, "\n");
			}
			fclose(connections_FP);
		}
        
        if (mp->delayEI) // Print EI delays
        {
            for (l=0; l<mp->nLayers; l++)
			{
				slen = snprintf(filename, FNAMEBUFF, "L%daffDelaysEI.dat", l); 
				assert(slen < FNAMEBUFF);
				connections_FP = myfopen(filename, "w");
				for (n=0; n<mp->vInhib[l]; n++)
				{
					fprintf(connections_FP, "%d\t", n_I[l][n].nLAff_E); // Print number of EI synapses first (in case any have 0)
					for (s=0; s<n_I[l][n].nLAff_E; s++)
						fprintf(connections_FP, "%d ", n_I[l][n].LAffs_E[s]->delay);
					fprintf(connections_FP, "\n");
				}
				fclose(connections_FP);
			}
        }
		
		for (l=0; l<mp->nLayers; l++)
		{
			slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsIE.dat", l);
			assert(slen < FNAMEBUFF);
			connections_FP = myfopen(filename, "w");
			for (n=0; n<mp->vExcit[l]; n++)
			{
				fprintf(connections_FP, "%d\t", n_E[l][n].nLAff_I); // Print number of IE synapses first
				for (s=0; s<n_E[l][n].nLAff_I; s++) // if (mp->pCnxIE[l] > EPS)
					fprintf(connections_FP, "%d ", n_E[l][n].lm0presyn_I[s]->n);
				fprintf(connections_FP, "\n");
			}
			fclose(connections_FP);
		}
		
		for (l=0; l<mp->nLayers; l++)
		{
			slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsII.dat", l);
			assert(slen < FNAMEBUFF);
			connections_FP = fopen(filename, "w");
			for (n=0; n<mp->vInhib[l]; n++)
			{
				fprintf(connections_FP, "%d\t", n_I[l][n].nLAff_I); // Print number of II synapses first
				for (s=0; s<n_I[l][n].nLAff_I; s++) // if (mp->pCnxII[l] > EPS)
					fprintf(connections_FP, "%d ", n_I[l][n].lm0presyn_I[s]->n);
				fprintf(connections_FP, "\n");
			}
			fclose(connections_FP);
		}
	}
	
	/************** END OF FILE OUTPUT **************/
	
	return; //void?
}


float calcDistance(NEURON * n1, NEURON * n2, float scale) //DIM * D)
{
	float deltaW = fabs(n1->x - n2->x);
	float deltaH = fabs(n1->y - n2->y);
	float deltaD = (n1->l - n2->l); //* scale;
	//int l = MIN(n1->l, n2->l); // Use lower layer since convergence is FF
	/*fprintf(stderr, "dW = %f; dH = %f; scale = %f; mp->spatialScale = %f;",deltaW,deltaH,scale,mp->spatialScale);
	float minW = (fabs(deltaW) < fabs(deltaW + mp->spatialScale)) ? (deltaW) : (deltaW + mp->spatialScale);
	float minH = (fabs(deltaH) < fabs(deltaH + mp->spatialScale)) ? (deltaH) : (deltaH + mp->spatialScale);
	fprintf(stderr, "minW = %f; minH = %f;\n",minW, minH);*/
	deltaW = MIN(fabs(deltaW), fabs(deltaW - scale)); //D[l].nCols)); // Periodic boundary conditions		
	deltaH = MIN(fabs(deltaH), fabs(deltaH - scale)); //D[l].nRows)); // Periodic boundary conditions
	return sqrt(pow(deltaH,2) + pow(deltaW,2) + pow(deltaD,2));
}


void loadAfferents(const char * suffix)
{
	char * str, fname[FNAMEBUFF], buff[131072]; //128KB
	char * delims = " \t";
	int sLen=0, l=0, n=0, s=0;
	FILE * affFP;
	NEURON * N;
	
	for (l=0; l<mp->nLayers; l++) // Same order of processing to preserve efferent synapses order
	{
		// Check for passed filename or archive file???

		if (l>0) // Load (afferent) EfE connections
		{
			sLen = snprintf(fname, FNAMEBUFF, "L%daffNeuronsEfE%s.dat", l, suffix);
			assert(sLen < FNAMEBUFF);
			affFP = myfopen(fname, "r"); //Feed-forward weights
			for (n=0; n<mp->vExcit[l]; n++)
			{
				if (!(str = fgets(buff, sizeof(buff), affFP)))
					exit_error("loadConnections", "NULL string reading EfE connections");
				N = &n_E[l][n];
				N->nFAff_E = atoi(strtok(buff, delims)); // First number is #presyn connections
				assert(N->nFAff_E <= mp->vExcit[l-1]);
				N->lm1presyn_E = myalloc(N->nFAff_E * sizeof(NEURON *));
				N->FAffs_E = myalloc(N->nFAff_E * sizeof(AXON *));
				for (s=0; s<N->nFAff_E && (str = strtok(NULL, delims))!=NULL; s++)
				{
					n_E[l][n].lm1presyn_E[s] = &n_E[l-1][atoi(str)];
					n_E[l-1][atoi(str)].nFEff_E++;
				}
			}
			fclose(affFP);
		}
		
		// Load (afferent) ElE connections
		sLen = snprintf(fname, FNAMEBUFF, "L%daffNeuronsElE%s.dat", l, suffix);
		assert(sLen < FNAMEBUFF);
		affFP = myfopen(fname, "r"); //Lateral weights
		for (n=0; n<mp->vExcit[l]; n++)
		{
			if (!(str = fgets(buff, sizeof(buff), affFP)))
				exit_error("loadConnections", "NULL string reading ElE connections");
			N = &n_E[l][n];
			N->nLAff_E = atoi(strtok(buff, delims)); // First number is #presyn connections
			assert(N->nLAff_E <= mp->vExcit[l]);
			N->lm0presyn_E = myalloc(N->nLAff_E * sizeof(NEURON *));
			N->LAffs_E = myalloc(N->nLAff_E * sizeof(AXON *));
			for (s=0; s<N->nLAff_E && (str = strtok(NULL, delims))!=NULL; s++)
			{
				N->lm0presyn_E[s] = &n_E[l][atoi(str)];
				n_E[l][atoi(str)].nLEff_E++;
			}
		}
		fclose(affFP);
		
		// Load (afferent) IlE connections
		sLen = snprintf(fname, FNAMEBUFF, "L%daffNeuronsIE%s.dat", l, suffix);
		assert(sLen < FNAMEBUFF);
		affFP = myfopen(fname, "r"); //Lateral weights
		for (n=0; n<mp->vExcit[l]; n++)
		{
			if (!(str = fgets(buff, sizeof(buff), affFP)))
				exit_error("loadConnections", "NULL string reading IE connections");
			N = &n_E[l][n];
			N->nLAff_I = atoi(strtok(buff, delims)); // First number is #presyn connections
			assert(N->nLAff_I <= mp->vInhib[l]);
			N->lm0presyn_I = myalloc(N->nLAff_I * sizeof(NEURON *));
			N->LAffs_I = myalloc(N->nLAff_I * sizeof(AXON *));
			for (s=0; s<N->nLAff_I && (str = strtok(NULL, delims))!=NULL; s++)
			{
				N->lm0presyn_I[s] = &n_I[l][atoi(str)];
				n_I[l][atoi(str)].nLEff_E++;
			}
		}
		fclose(affFP);
		
		// Load (afferent) ElI connections
		sLen = snprintf(fname, FNAMEBUFF, "L%daffNeuronsEI%s.dat", l, suffix);
		assert(sLen < FNAMEBUFF);
		affFP = myfopen(fname, "r"); //Lateral weights
		for (n=0; n<mp->vInhib[l]; n++)
		{
			if (!(str = fgets(buff, sizeof(buff), affFP)))
				exit_error("loadConnections", "NULL string reading EI connections");
			N = &n_I[l][n];
			N->nLAff_E = atoi(strtok(buff, delims)); // First number is #presyn connections
			assert(N->nLAff_E <= mp->vExcit[l]);
			N->lm0presyn_E = myalloc(N->nLAff_E * sizeof(NEURON *));
			N->LAffs_E = myalloc(N->nLAff_E * sizeof(AXON *));
			for (s=0; s<N->nLAff_E && (str = strtok(NULL, delims))!=NULL; s++)
			{
				N->lm0presyn_E[s] = &n_E[l][atoi(str)];
				n_E[l][atoi(str)].nLEff_I++;
			}
		}
		fclose(affFP);
		
		// Load (afferent) IlI connections
		sLen = snprintf(fname, FNAMEBUFF, "L%daffNeuronsII%s.dat", l, suffix);
		assert(sLen < FNAMEBUFF);
		affFP = myfopen(fname, "r"); //Lateral weights
		for (n=0; n<mp->vInhib[l]; n++)
		{
			if (!(str = fgets(buff, sizeof(buff), affFP)))
				exit_error("loadConnections", "NULL string reading II connections");
			N = &n_I[l][n];
			N->nLAff_I = atoi(strtok(buff, delims)); // First number is #presyn connections
			assert(N->nLAff_I <= mp->vInhib[l]);
			N->lm0presyn_I = myalloc(N->nLAff_I * sizeof(NEURON *));
			N->LAffs_I = myalloc(N->nLAff_I * sizeof(AXON *));
			for (s=0; s<N->nLAff_I && (str = strtok(NULL, delims))!=NULL; s++)
			{
				N->lm0presyn_I[s] = &n_I[l][atoi(str)];
				n_I[l][atoi(str)].nLEff_I++;
			}
		}
		fclose(affFP);
	}
	
}

// Try using this macro trick to condense wiring routines by replacing lm0presynE etc.
// #define BUILD_FIELD(field) my_struct.inner_struct.union_a.##field
// Now, when used with a particular field name, it will expand to something like
// my_struct.inner_struct.union_a.field1

void wireAfferents(NEURON * n, float pEfn, float pEln, float pIln)
{
	int s;//,l;
	int buff;
	
	/*** Ef synapses ***/
	if (pEfn > EPS) // (n->type == EXCIT && n->l>0)
	{
		buff = ceil(pEfn * mp->vExcit[n->l-1]);
		n->lm1presyn_E = myalloc(buff * sizeof(NEURON *));
		n->nFAff_E = 0;
		for (s=0; s<mp->vExcit[n->l-1]; s++) //if (n->nFAff_E > 0)
		{
			if (gsl_rng_uniform(mSeed) < pEfn) // Make presynaptic connection //ran3(&idum)
			{
				n->lm1presyn_E[n->nFAff_E++] = &n_E[n->l-1][s];
				n_E[n->l-1][s].nFEff_E++; // Critical section if this function is parallelised
				if (n->nFAff_E == buff)
				{
					buff = (buff < mp->vExcit[n->l-1]/2) ? ceil(2 * buff) : mp->vExcit[n->l-1];
					n->lm1presyn_E = myrealloc(n->lm1presyn_E, buff * sizeof(NEURON *));
				}
			}
		}
		n->lm1presyn_E = myrealloc(n->lm1presyn_E, n->nFAff_E * sizeof(NEURON *));
		n->FAffs_E = myalloc(n->nFAff_E * sizeof(AXON *)); // Create array of AXON pointers
	}
	
	/*if (pEln < 0.0) // SOM architecture
	{
		
	}*/
	
	/*** Lateral El synapses ***/
	if (pEln > EPS)
	{
		buff = ceil(pEln * mp->vExcit[n->l]);
		n->lm0presyn_E = myalloc(buff * sizeof(NEURON *)); // Create array of NEURON pointers
		n->nLAff_E = 0;
		float cutoff = mp->SOMclip * mp->SOMsigE; //int?
		for (s=0; s<mp->vExcit[n->l]; s++)
		{
			if ((n->n == n_E[n->l][s].n) && n->type==EXCIT) // Do not self synapse (i != j)
				continue;
			
			/* Collapse this branch!!! */
			if (n->type==EXCIT && (readLowTriF(distE, n->l, n->n, n_E[n->l][s].n) < cutoff)) // ElE of SOM within range //mp->SOM && 
			{
				if(gsl_rng_uniform(mSeed) < pEln) // Connect //ran3(&idum)
				{
					n->lm0presyn_E[n->nLAff_E++] = &n_E[n->l][s];
					n_E[n->l][s].nLEff_E++;
					if (n->nLAff_E == buff)
					{
						buff = (buff < mp->vExcit[n->l]/2) ? ceil(2 * buff) : mp->vExcit[n->l];
						n->lm0presyn_E = myrealloc(n->lm0presyn_E, buff * sizeof(NEURON *));
					}
				}
			}
			else // ElI || !SOM && ElE
			{
				if(gsl_rng_uniform(mSeed) < pEln) // Connect //ran3(&idum)
				{//if (mp->SOM && n->type==EXCIT && (readLowTriF(distE, n->l, n->n, n_E[n->l][s].n) > mp->SOMclip*mp->SOMsigE)); continue; // Out of range
					n->lm0presyn_E[n->nLAff_E++] = &n_E[n->l][s];
					if (n->type == EXCIT)
						n_E[n->l][s].nLEff_E++;
					else if (n->type == INHIB)
						n_E[n->l][s].nLEff_I++;
					if (n->nLAff_E == buff)
					{
						buff = (buff < mp->vExcit[n->l]/2) ? ceil(2 * buff) : mp->vExcit[n->l];
						n->lm0presyn_E = myrealloc(n->lm0presyn_E, buff * sizeof(NEURON *));
					}
				}
			}
		}
		n->lm0presyn_E = myrealloc(n->lm0presyn_E, n->nLAff_E * sizeof(NEURON *));
		n->LAffs_E = myalloc(n->nLAff_E * sizeof(AXON *)); // Create array of AXON pointers
	}
	
	/*** Lateral I Synapses ***/
	if (pIln > EPS)
	{
		buff = ceil(pIln * mp->vInhib[n->l]);
		n->lm0presyn_I = myalloc(buff * sizeof(NEURON *));
		n->nLAff_I = 0;
		for (s=0; s<mp->vInhib[n->l]; s++)
		{
			if ((n->n == n_I[n->l][s].n) && n->type==INHIB) // Do not self synapse
				continue;
			
			if (gsl_rng_uniform(mSeed) < pIln) //ran3(&idum)
			{
				n->lm0presyn_I[n->nLAff_I++] = &n_I[n->l][s]; // Point to presynaptic neuron
				if (n->type == EXCIT)
					n_I[n->l][s].nLEff_E++;
				else if (n->type == INHIB)
					n_I[n->l][s].nLEff_I++;
				if (n->nLAff_I == buff)
				{
					buff = (buff < mp->vInhib[n->l]/2) ? ceil(2 * buff) : mp->vInhib[n->l];
					n->lm0presyn_I = myrealloc(n->lm0presyn_I, buff * sizeof(NEURON *));
				}
			}
		}
		n->lm0presyn_I = myrealloc(n->lm0presyn_I, n->nLAff_I * sizeof(NEURON *));
		n->LAffs_I = myalloc(n->nLAff_I * sizeof(AXON *)); // Create array of AXON pointers
	}
		
	return;
}


void alloc_efferents(NEURON * n)
{
	//if (n->FEffs_E > 0)
	n->FEffs_E = myalloc(n->nFEff_E * sizeof(AXON));
	n->lp1postsyn_E = myalloc(n->nFEff_E * sizeof(NEURON *));
	n->nFEff_E = 0; // Reset to 0 as a synapse counter in wire_afferents
	
	//if (n->LEffs_E > 0)
	n->LEffs_E = myalloc(n->nLEff_E * sizeof(AXON));
	n->lp0postsyn_E = myalloc(n->nLEff_E * sizeof(NEURON *));
	n->nLEff_E = 0; // Reset to 0 as a synapse counter in wire_afferents
	
	//if (n->LEffs_I > 0)
	n->LEffs_I = myalloc(n->nLEff_I * sizeof(AXON));
	n->lp0postsyn_I = myalloc(n->nLEff_I * sizeof(NEURON *));
	n->nLEff_I = 0; // Reset to 0 as a synapse counter in wire_afferents
	return;
}

// Parallelizable? // Split into two functions?
void wire_efferents(NEURON * n)
{
	int s, effCount;
	NEURON * pre_n;
	/* This function takes a NEURON * and wires the efferents of all its pre-synaptic neurons to it */
	if (n->type == EXCIT)		/* Excitatory Post-synaptic neurons */
	{
		// Pre-synaptic : EXCIT (f) | Post-synaptic : EXCIT
		for (s=0; s<n->nFAff_E; s++)
		{
			pre_n = (NEURON *) n->lm1presyn_E[s];
			effCount = pre_n->nFEff_E++;
			n->FAffs_E[s] = &pre_n->FEffs_E[effCount];	// Point to presynaptic efferent axon
			pre_n->lp1postsyn_E[effCount] = n;			// Point to postsynaptic neuron
		}
		
		// Pre-synaptic : EXCIT (l) | Post-synaptic : EXCIT
		for (s=0; s<n->nLAff_E; s++)
		{
			pre_n = (NEURON *) n->lm0presyn_E[s];
			effCount = pre_n->nLEff_E++; // records how many post-syn connections have been wired up
			n->LAffs_E[s] = &pre_n->LEffs_E[effCount];	// Point to presynaptic efferent axon
			pre_n->lp0postsyn_E[effCount] = n;			// Point to postsynaptic neuron
		}
		
		// Pre-synaptic : INHIB | Post-synaptic: EXCIT
		for (s=0; s<n->nLAff_I; s++)
		{
			pre_n = (NEURON *) n->lm0presyn_I[s];
			effCount = pre_n->nLEff_E++; // records how many post-syn connections have been wired up
			n->LAffs_I[s] = &pre_n->LEffs_E[effCount];	// Point to presynaptic efferent axon
			pre_n->lp0postsyn_E[effCount] = n;			// Point to postsynaptic neuron
		}
	}
	
	else if (n->type == INHIB)		/* Inhibitory Post-synaptic neurons */
	{
		// Pre-synaptic : EXCIT | Post-synaptic : INHIB
		for (s=0; s<n->nLAff_E; s++)
		{
			pre_n = (NEURON *) n->lm0presyn_E[s];
			assert(n->lm0presyn_E[s]->type == EXCIT && n->type == INHIB);
			effCount = pre_n->nLEff_I++; // records how many post-syn connections have been wired up
			n->LAffs_E[s] = &pre_n->LEffs_I[effCount];	// Point to presynaptic efferent axon
			pre_n->lp0postsyn_I[effCount] = n;			// Point to postsynaptic neuron
		}
		
		// Pre-synaptic : INHIB | Post-synaptic : INHIB
		for (s=0; s<n->nLAff_I; s++)
		{
			pre_n = (NEURON *) n->lm0presyn_I[s];
			assert(n->lm0presyn_I[s]->type == INHIB && n->type == INHIB);
			effCount = pre_n->nLEff_I++; // records how many post-syn connections have been wired up
			n->LAffs_I[s] = &pre_n->LEffs_I[effCount];	// Point to presynaptic efferent axon
			pre_n->lp0postsyn_I[effCount] = n;			// Point to postsynaptic neuron
		}
	}
}

void create_axons(NEURON * n, PARAMS * mp)
{
	/*** Axonal delays specified in seconds must be converted to timesteps ***/
	tstep delay = 0;
	//int nBins = 0;
	int span = 0;
	int s = 0;

	
    if (n->type == EXCIT) // Excitatory neuron -> {E,I} delays
    {
        // EfE Delays
        switch (mp->delayEfE) //(mp->axonDelay) 
        {
            case MinD:
                delay = 1;
                break;
            case ConstD:
                delay = round(mp->d_const/mp->DT);
                delay = (!delay) ? 1 : delay;
                break;
            case UniformD:
                span = mp->d_max - mp->d_min;
                break;
            case GaussD:
                break;
            case SOMD:
                break;
            default:
                exit_error("create_axons", "Unknown axonal delay model!");
                break;
        }
        for (s=0; s<n->nFEff_E; s++) // nFEff_E = 0 for l = nLayers and n->type == INHIB
        {
            switch (mp->delayEfE)//(mp->axonDelay) 
            {
                case UniformD:
                    delay = round(((gsl_rng_uniform(mSeed)*span)+mp->d_min)/mp->DT);
                    break;
                case GaussD:
                    delay = round((mp->d_mean + gsl_ran_gaussian(mSeed,mp->d_sd))/mp->DT);
                    break;
                case SOMD:
                    delay = 1;
                    break;
                default:
                    break;
            }
            n->FEffs_E[s].delay = (delay<1) ? 1 : delay; //(!delay) ? 1 : delay;
            n->FEffs_E[s].nBins = ceil((n->FEffs_E[s].delay*mp->DT)/mp->refract);
            n->FEffs_E[s].queue = myalloc(n->FEffs_E[s].nBins * sizeof(tstep));
            init_queue(&(n->FEffs_E[s]));
        }
        
        // ElE Delays
        switch (mp->delayElE) //(mp->axonDelay) 
        {
            case MinD:
                delay = 1;
                break;
            case ConstD:
                delay = round(mp->d_const/mp->DT);
                delay = (!delay) ? 1 : delay;
                break;
            case UniformD:
                span = mp->d_max - mp->d_min;
                break;
            case GaussD:
                break;
            case SOMD:
                break;
            default:
                exit_error("create_axons", "Unknown axonal delay model!");
                break;
        }
        for (s=0; s<n->nLEff_E; s++)
        {
            switch (mp->delayElE) //(mp->axonDelay) 
            {
                case UniformD:
                    delay = round(((gsl_rng_uniform(mSeed)*span)+mp->d_min)/mp->DT);
                    break;
                case GaussD:
                    delay = round((mp->d_mean + gsl_ran_gaussian(mSeed,mp->d_sd))/mp->DT);
                    break;
                case SOMD: 	// Wire ElE connections with delays proportional to their Euclidean distances
                    if (n->type == EXCIT)
                        delay = round(readLowTriF(distE, n->l, n->n, n->lp0postsyn_E[s]->n) / ((float) mp->condSpeed * mp->DT));
                    else
                        delay = 1; // *** Need to add topology to Inhibitory neurons ***
					//delay = round(distE[n->l][n->n][n->lp0postsyn_E[s]->n] / mp->condSpeed);
                    break;
                default:
                    break;
            }
            n->LEffs_E[s].delay = (delay<1) ? 1 : delay;
            n->LEffs_E[s].nBins = ceil((n->LEffs_E[s].delay*mp->DT)/mp->refract);
            n->LEffs_E[s].queue = myalloc(n->LEffs_E[s].nBins * sizeof(tstep));
            init_queue(&(n->LEffs_E[s]));
        }
        
        // EI Delays
        switch (mp->delayEI) //(mp->axonDelay) 
        {
            case MinD:
                delay = 1;
                break;
            case ConstD:
                delay = round(mp->d_const/mp->DT);
                delay = (!delay) ? 1 : delay;
                break;
            case UniformD:
                span = mp->d_max - mp->d_min;
                break;
            case GaussD:
                break;
            case SOMD:
                break;
            default:
                exit_error("create_axons", "Unknown axonal delay model!");
                break;
        }
        for (s=0; s<n->nLEff_I; s++)
        {
            switch (mp->delayEI) //(mp->axonDelay) 
            {
                case UniformD:
                    delay = round(((gsl_rng_uniform(mSeed)*span)+mp->d_min)/mp->DT);
                    break;
                case GaussD:
                    delay = round((mp->d_mean + gsl_ran_gaussian(mSeed,mp->d_sd))/mp->DT);
                    break;
                case SOMD:
                    delay = 1;
                    break;
                default:
                    break;
            }
            n->LEffs_I[s].delay = (delay<1) ? 1 : delay;
            n->LEffs_I[s].nBins = ceil((n->LEffs_I[s].delay*mp->DT)/mp->refract);
            n->LEffs_I[s].queue = myalloc(n->LEffs_I[s].nBins * sizeof(tstep));
            init_queue(&(n->LEffs_I[s]));
        }
    }
    else // Inhibitory Neuron -> {E,I} delays
    {
        for (s=0; s<n->nLEff_E; s++)
        {
            n->LEffs_E[s].delay = 1;
            n->LEffs_E[s].nBins = 1;
            n->LEffs_E[s].queue = myalloc(n->LEffs_E[s].nBins * sizeof(tstep));
            init_queue(&(n->LEffs_E[s]));
        }
        
        for (s=0; s<n->nLEff_I; s++)
        {
            n->LEffs_I[s].delay = 1;
            n->LEffs_I[s].nBins = 1;
            n->LEffs_I[s].queue = myalloc(n->LEffs_I[s].nBins * sizeof(tstep));
            init_queue(&(n->LEffs_I[s]));
        }
    }
	return;
}


void initNetwork(SETSV resetBuffers)
{
	int s, n, l;
	
	/********** Training and Testing **********/
	// Executed for all types of initialization (except Settle where shown)
#pragma omp parallel default(shared) private(l,n,s)//,seed)
	{

	/* Membrane potentials, conductances, learning parameters and spike time arrays */
	for (l=0; l<mp->nLayers; l++)
	{
#pragma omp for nowait
		for (n=0; n<mp->vExcit[l]; n++)
		{
			if (mp->noise) // mp->VrestE + gsl_ran_gaussian(mSeed, mp->SigmaE) old noise model
			{
#ifdef _OPENMP
				n_E[l][n].V = n_E[l][n].V_tm1 = (gsl_rng_uniform(states[omp_get_thread_num()]) * fabs(mp->VhyperE - mp->ThreshE)) + mp->VhyperE;
				//seed = states[omp_get_thread_num()];
#else
				n_E[l][n].V = n_E[l][n].V_tm1 = (gsl_rng_uniform(mSeed) * fabs(mp->VhyperE - mp->ThreshE)) + mp->VhyperE;
				//seed = mSeed;
#endif
				//n_E[l][n].V = n_E[l][n].V_tm1 = (gsl_rng_uniform(seed) * fabs(mp->VhyperE - mp->ThreshE)) + mp->VhyperE;
			}
			else
				n_E[l][n].V = n_E[l][n].V_tm1 = mp->VrestE;

			if (mp->adaptation)
				n_E[l][n].cCa = n_E[l][n].cCa_tm1 = 0.0;
			n_E[l][n].D = n_E[l][n].D_tm1 = 0.0;
			for (s=0; s<n_E[l][n].nFEff_E; s++) // nFEff_E should be 0 for the last layer
			{
				n_E[l][n].FEffs_E[s].C = n_E[l][n].FEffs_E[s].C_tm1 = 0.0;
				n_E[l][n].FEffs_E[s].g = n_E[l][n].FEffs_E[s].g_tm1 = 0.0; // Randomly initialise conductances?
				init_queue(&(n_E[l][n].FEffs_E[s]));
			}
			for (s=0; s<n_E[l][n].nLEff_E; s++)
			{
				n_E[l][n].LEffs_E[s].C = n_E[l][n].LEffs_E[s].C_tm1 = 0.0;
				n_E[l][n].LEffs_E[s].g = n_E[l][n].LEffs_E[s].g_tm1 = 0.0;
				init_queue(&(n_E[l][n].LEffs_E[s]));
			}
			for (s=0; s<n_E[l][n].nLEff_I; s++)
			{
				n_E[l][n].LEffs_I[s].g = n_E[l][n].LEffs_I[s].g_tm1 = 0.0;
				init_queue(&(n_E[l][n].LEffs_I[s]));
			}
			
			n_E[l][n].lastSpike = -BIG;
			n_E[l][n].nextUpdate = -BIG;
			
			if (resetBuffers == Hard) // Reinitialize spike buffers between epochs and phases
			{
				n_E[l][n].spkbin = 0;
				//bins = (mp->useFilteredImages && l==0) ? mp->inpSpkBuff : mp->spkBuffer;
#ifndef __llvm__ // The new LLVM-GCC compiler has a problem with this memset
				memset(n_E[l][n].spikeTimes, 0, mp->spkBuffer*sizeof(n_E[l][n].spikeTimes[0])); //[0]?
#endif
			}
		}
	}

	for (l=0; l<mp->nLayers; l++)
	{
#pragma omp for nowait
		for (n=0; n<mp->vInhib[l]; n++)
		{
			if (mp->noise)
			{
#ifdef _OPENMP
				n_I[l][n].V = n_I[l][n].V_tm1 = (gsl_rng_uniform(states[omp_get_thread_num()]) * fabs(mp->VhyperI - mp->ThreshI)) + mp->VhyperI;
#else
				n_I[l][n].V = n_I[l][n].V_tm1 = (gsl_rng_uniform(mSeed) * fabs(mp->VhyperI - mp->ThreshI)) + mp->VhyperI;
#endif
			}
			else
				n_I[l][n].V = n_I[l][n].V_tm1 = mp->VrestI;	

			for (s=0; s<n_I[l][n].nLEff_E; s++)
			{
				n_I[l][n].LEffs_E[s].g = n_I[l][n].LEffs_E[s].g_tm1 = 0.0;
				init_queue(&(n_I[l][n].LEffs_E[s]));
			}
			for (s=0; s<n_I[l][n].nLEff_I; s++)
			{
				n_I[l][n].LEffs_I[s].g = n_I[l][n].LEffs_I[s].g_tm1 = 0.0;
				init_queue(&(n_I[l][n].LEffs_I[s]));
			}
			
			n_I[l][n].lastSpike = -BIG;
			n_I[l][n].nextUpdate = -BIG;
			
			if (resetBuffers == Hard) //(regime != Settle)
			{
				n_I[l][n].spkbin = 0;
#ifndef __llvm__ //&& #ifdef __GNUC__
				memset(n_I[l][n].spikeTimes, 0, mp->spkBuffer*sizeof(n_I[l][n].spikeTimes[0]));
#endif
			}
		}
	}
	
	} // End of parallel region

	return;
}


void setRecords(PARAMS * mp, NEURON ** n_E, gsl_rng * mSeed) // if (mp->nRecords)
{
	int l, n, r;
	int * choices = NULL;
	int * chosen = NULL;
	
	// Print Records file
	char rString[BUFSIZ];
	int slen = 0;
	FILE * rFile = myfopen("records.m", "w");
	fprintf(rFile, "MP.Records = cell(1,MP.nLayers);\n"); 
		
	/* Set up recording bins */
	printf("\n");
	for (l=0; l<mp->nLayers; l++)
	{
		assert((0 <= mp->vRecords[l]) && (mp->vRecords[l] <= mp->vExcit[l]));
		chosen = myalloc(mp->vRecords[l] * sizeof(*chosen));
		choices = myalloc(mp->vExcit[l] * sizeof(*choices));
		for (n=0; n<mp->vExcit[l]; n++)
			choices[n] = n;
	
		gsl_ran_choose(mSeed, chosen, mp->vRecords[l], choices, mp->vExcit[l], sizeof(*choices));
		for (r=0; r<mp->vRecords[l]; r++) // Randomly set NRECORDS flags
		{
			n_E[l][chosen[r]].rec_flag = true;
			printf("\tLayer #%d, Record %d Assigned nID: #%d\n",l,r+1,chosen[r]);
		}
		
		slen = snprintf(rString, BUFSIZ, "MP.Records{%d}", l+1);
		assert(slen < BUFSIZ);
		//for (r=0; r<mp->vRecords[l]; r++)
		//	chosen[r] += 1; // Make matlab friendly to match layer indexing?
		printIntArray(rFile, rString, chosen, mp->vRecords[l]); 
		myfree(chosen);
		myfree(choices);
	}
	fclose(rFile); // Close Records file
	
	/* Create recording structures - (mp->RecordMS+1) so that initial conditions are recorded */
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vExcit[l]; n++)
			if (n_E[l][n].rec_flag) // Could move inside flag loop above
			{
				n_E[l][n].rec = myalloc(sizeof(RECORD));
                n_E[l][n].rec->cellcCa = NULL;
                n_E[l][n].rec->cellD = NULL;
                n_E[l][n].rec->LsigGI = NULL;
                n_E[l][n].rec->FSynDG = NULL;
                n_E[l][n].rec->FSynC = NULL;
                n_E[l][n].rec->LSynDG = NULL;
                n_E[l][n].rec->LSynC = NULL;
				n_E[l][n].rec->bin = 0;
				n_E[l][n].rec->cellV = myalloc((mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->cellV)));
				memset(n_E[l][n].rec->cellV, 0, (mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->cellV)));
				if (mp->adaptation)
				{
					n_E[l][n].rec->cellcCa = myalloc((mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->cellcCa)));
					memset(n_E[l][n].rec->cellcCa, 0, (mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->cellcCa)));
				}
                if (n_E[l][n].nLAff_I)
                {
                    n_E[l][n].rec->LsigGI = myalloc((mp->RecordMS+1)*sizeof(*(n_E[l][n].rec->LsigGI)));
                    memset(n_E[l][n].rec->LsigGI, 0, (mp->RecordMS+1)*sizeof(*(n_E[l][n].rec->LsigGI)));
                }
                if (n_E[l][n].nFAff_E) //(l>0) // Here the synaptic variables are stored with the post-synaptic neuron
                    n_E[l][n].rec->FSynG = get_2D_farray(n_E[l][n].nFAff_E, (mp->RecordMS+1), 0.0);
                if (n_E[l][n].nLAff_E) //(mp->pCnxElE[l]>EPS)
                    n_E[l][n].rec->LSynG = get_2D_farray(n_E[l][n].nLAff_E, (mp->RecordMS+1), 0.0);
                if (mp->train)
                {
                    n_E[l][n].rec->cellD = myalloc((mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->cellD)));
                    memset(n_E[l][n].rec->cellD, 0, (mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->cellD)));
                    if (n_E[l][n].nFAff_E) //(l>0)
                    {
                        n_E[l][n].rec->FSynDG = get_2D_farray(n_E[l][n].nFAff_E, (mp->RecordMS+1), 0.0);
                        n_E[l][n].rec->FSynC = get_2D_farray(n_E[l][n].nFAff_E, (mp->RecordMS+1), 0.0);
                    }
                    if (mp->trainElE && n_E[l][n].nLAff_E) // Create ElE records //(mp->pCnxElE[l]>EPS))//mp->train) 
                    {
                        n_E[l][n].rec->LSynDG = get_2D_farray(n_E[l][n].nLAff_E, (mp->RecordMS+1), 0.0);
                        n_E[l][n].rec->LSynC = get_2D_farray(n_E[l][n].nLAff_E, (mp->RecordMS+1), 0.0);
                    }
                }
			}
	return;
}

/*int resetRecords() // Reinitialise Records
{
	if (mp->nRecords) // Should be training //sPhase==Testing && 
	{
		if (sPhase == Training)
		{
			slen = snprintf(prefix, FNAMEBUFF, "RE%d", loop);
			assert(slen && slen < FNAMEBUFF); // Check non-negative
		}
		for (l=0; l<mp->nLayers; l++)
		{
			for (n=0; n<mp->vExcit[l]; n++)
			{
				if (n_E[l][n].rec_flag)
				{
					slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dV.dat", prefix, l, n); // Change
					assert(slen && slen < FNAMEBUFF);
					recFile = myfopen(filename, "w"); //r_cellV_ptr = fopen(filename, "wb");
					print_frow(recFile, n_E[l][n].rec->cellV, n_E[l][n].rec->bin); // bin already points to the next free slot
					// print_farray(rCellVout, n_E[l][n].rec->cellV, mp->loops, mp->RecordMS); //fwrite(n_E[l][n].rec->cellV, sizeof(float), mp->loops*mp->RecordMS, r_cellV_ptr); // See nifty trick #1
					fclose(recFile);
					memset(n_E[l][n].rec->cellV, 0, (mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->cellV)));
					
					if (mp->adaptation)
					{
						slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dcCa.dat", prefix, l, n);
						assert(slen && slen < FNAMEBUFF);
						recFile = myfopen(filename, "w");
						print_frow(recFile, n_E[l][n].rec->cellcCa, n_E[l][n].rec->bin);
						fclose(recFile);
						memset(n_E[l][n].rec->cellcCa, 0, (mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->cellcCa)));
					}
					
					slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dD.dat", prefix, l, n);
					assert(slen && slen < FNAMEBUFF);
					recFile = myfopen(filename, "w"); //r_D_ptr = fopen(filename, "wb");
					print_frow(recFile, n_E[l][n].rec->cellD, n_E[l][n].rec->bin); // (rDout, n_E[l][n].rec->cellD, mp->loops, mp->RecordMS);
					fclose(recFile);
					memset(n_E[l][n].rec->cellD, 0, (mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->cellD)));
					
					if (l>0) // The presynaptic cell's values are attached to each record
					{
						slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dFAffC.dat", prefix, l, n);
						assert(slen && slen < FNAMEBUFF);
						recFile = myfopen(filename, "w"); //r_C_ptr = fopen(filename, "wb");
						//for (loop=0; loop<mp->loops; loop++)
						print_farray(recFile, n_E[l][n].rec->FSynC, n_E[l][n].nFAff_E, n_E[l][n].rec->bin); //print_farray(recFile, n_E[l][n].rec->SynC[loop], n_E[l][n].nFAff_E, mp->RecordMS);
						fclose(recFile);
						memset(n_E[l][n].rec->FSynC[0], 0, n_E[l][n].nFAff_E*(mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->FSynC)));
						
						slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dFAffg.dat", prefix, l, n);
						assert(slen && slen < FNAMEBUFF);
						recFile = myfopen(filename, "w");
						//for (loop=0; loop<mp->loops; loop++)
						print_farray(recFile, n_E[l][n].rec->FSynG, n_E[l][n].nFAff_E, n_E[l][n].rec->bin);
						fclose(recFile);
						memset(n_E[l][n].rec->FSynG[0], 0, n_E[l][n].nFAff_E*(mp->RecordMS+1) * sizeof(**(n_E[l][n].rec->FSynG)));
						
						slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dFAffdg.dat", prefix, l, n);
						assert(slen && slen < FNAMEBUFF);
						recFile = myfopen(filename, "w");
						//for (loop=0; loop<mp->loops; loop++)
						print_farray(recFile, n_E[l][n].rec->FSynDG, n_E[l][n].nFAff_E, n_E[l][n].rec->bin);
						fclose(recFile);
						memset(n_E[l][n].rec->FSynDG[0], 0, n_E[l][n].nFAff_E*(mp->RecordMS+1) * sizeof(**(n_E[l][n].rec->FSynDG)));
					}
					
					if (mp->trainElE && mp->train && mp->pCnxElE[l] > EPS)
					{
						slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dLAffC.dat", prefix, l, n);
						assert(slen && slen < FNAMEBUFF);
						recFile = myfopen(filename, "w"); 
						print_farray(recFile, n_E[l][n].rec->LSynC, n_E[l][n].nLAff_E, n_E[l][n].rec->bin); 
						fclose(recFile);
						memset(n_E[l][n].rec->LSynC[0], 0, n_E[l][n].nLAff_E*(mp->RecordMS+1) * sizeof(**(n_E[l][n].rec->LSynC)));
						
						slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dLAffg.dat", prefix, l, n);
						assert(slen && slen < FNAMEBUFF);
						recFile = myfopen(filename, "w");
						print_farray(recFile, n_E[l][n].rec->LSynG, n_E[l][n].nLAff_E, n_E[l][n].rec->bin);
						fclose(recFile);
						memset(n_E[l][n].rec->LSynG[0], 0, n_E[l][n].nLAff_E*(mp->RecordMS+1) * sizeof(**(n_E[l][n].rec->LSynG)));
						
						slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dLAffdg.dat", prefix, l, n);
						assert(slen && slen < FNAMEBUFF);
						recFile = myfopen(filename, "w");
						print_farray(recFile, n_E[l][n].rec->LSynDG, n_E[l][n].nLAff_E, n_E[l][n].rec->bin);
						fclose(recFile);
						memset(n_E[l][n].rec->LSynDG[0], 0, n_E[l][n].nLAff_E*(mp->RecordMS+1) * sizeof(**(n_E[l][n].rec->LSynDG)));
					}
					n_E[l][n].rec->bin = 0; // Reset record counter ready for next phase/epoch
				}
			}
		}
	} // End of state variable records
}*/


void setWeights(PARAMS * mp, NEURON ** n_E, NEURON ** n_I, const char * suffix)
{
	int l, n, s;
	
	if (mp->loadWeights)
	{
		char * str, fname[FNAMEBUFF], buff[131072]; //128KB
		char * delims = " \t";
		int sLen=0, nSyn=0;
		FILE * weightsFP = NULL;
		//const char * suffix = ""; // Check for passed filename or archive file???
		for (l=0; l<mp->nLayers; l++)
		{
			if (l>0) // Load Feed-forward weights
			{
				sLen = snprintf(fname, FNAMEBUFF, "L%dweightsEfE%s.dat", l, suffix);
				assert(sLen < FNAMEBUFF);
				weightsFP = myfopen(fname, "r"); //Feed-forward weights
				for (n=0; n<mp->vExcit[l]; n++)
				{
					if (!(str = fgets(buff, sizeof(buff), weightsFP)))
						exit_error("loadWeights", "NULL string reading EfE weights");
					nSyn = atoi(strtok(buff, delims)); // First number in each row is #presyn connections
					assert(n_E[l][n].nFAff_E == nSyn);
					for (s=0; s<nSyn && (str = strtok(NULL, delims))!=NULL; s++)
						n_E[l][n].FAffs_E[s]->delta_g = n_E[l][n].FAffs_E[s]->delta_g_tm1 = atof(str);
				}
				fclose(weightsFP);
			}
			
			// Load Lateral weights
			sLen = snprintf(fname, FNAMEBUFF, "L%dweightsElE%s.dat", l, suffix);
			assert(sLen < FNAMEBUFF);
			weightsFP = myfopen(fname, "r"); //Lateral weights
			for (n=0; n<mp->vExcit[l]; n++)
			{
				if (!(str = fgets(buff, sizeof(buff), weightsFP)))
					exit_error("loadWeights", "NULL string reading ElE weights");
				nSyn = atoi(strtok(buff, delims)); // First number in each row is #presyn connections
				assert(n_E[l][n].nLAff_E == nSyn);
				for (s=0; s<nSyn && (str = strtok(NULL, delims))!=NULL; s++)
					n_E[l][n].LAffs_E[s]->delta_g = n_E[l][n].LAffs_E[s]->delta_g_tm1 = atof(str);
			}
			fclose(weightsFP);
		}
	}
	else // Set new weights
	{
		/* E_ Synaptic weights */
		
		switch (mp->initEfE) {
				/*case Zero:
				 for (l=0; l<mp->nWLayers; l++)
				 for (n=0; n<mp->vExcit[l]; n++)
				 for (s=0; s<n_E[l][n].nFEff_E; s++)
				 n_E[l][n].FEffs_E[s].delta_g = n_E[l][n].FEffs_E[s].delta_g_tm1 = 0.0;
				 break;*/
			case Constant:
				for (l=0; l<mp->nWLayers; l++)
					for (n=0; n<mp->vExcit[l]; n++)
						for (s=0; s<n_E[l][n].nFEff_E; s++)
							n_E[l][n].FEffs_E[s].delta_g = n_E[l][n].FEffs_E[s].delta_g_tm1 = mp->iEfE; //mp->DgEfE;
				break;
			case Uniform:
				for (l=0; l<mp->nWLayers; l++)
					for (n=0; n<mp->vExcit[l]; n++)
						for (s=0; s<n_E[l][n].nFEff_E; s++)
							n_E[l][n].FEffs_E[s].delta_g = n_E[l][n].FEffs_E[s].delta_g_tm1 = gsl_rng_uniform(mSeed);
				break;
			case Gaussian:
				for (l=0; l<mp->nWLayers; l++)
					for (n=0; n<mp->vExcit[l]; n++)
						for (s=0; s<n_E[l][n].nFEff_E; s++)
							n_E[l][n].FEffs_E[s].delta_g = n_E[l][n].FEffs_E[s].delta_g_tm1 = gsl_ran_gaussian(mSeed, mp->SOMsigE); // Create seperate sigma?
				break;
			case SOM: // Convergent feed-forward connections
				exit_error("setWeights", "Feed-forward convergent connections not yet implemented!");
				//fprintf(stderr, "Warning: Feed-forward convergent connections not yet implemented!\n");
				break;
			default:
				exit_error("setWeights", "Illegal EfE weight initialisation arguement!");
				break;
		}

		float distance = 0.0;
		float phiScale = 0.0; 
		switch (mp->initElE) {
				/*case Zero:
				 for (l=0; l<mp->nLayers; l++)
				 for (n=0; n<mp->vExcit[l]; n++)
				 for (s=0; s<n_E[l][n].nLEff_E; s++)
				 n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = 0.0;
				 break;*/
			case Constant: // mp->DgElE
				for (l=0; l<mp->nLayers; l++)
					for (n=0; n<mp->vExcit[l]; n++)
						for (s=0; s<n_E[l][n].nLEff_E; s++)
							n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = mp->iElE; //mp->DgElE;
				break;
			case Uniform: // mSeed
				for (l=0; l<mp->nLayers; l++)
					for (n=0; n<mp->vExcit[l]; n++)
						for (s=0; s<n_E[l][n].nLEff_E; s++)
							n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = gsl_rng_uniform(mSeed);
				break;
			case Gaussian: // mSeed, mp->SOMsigE
				for (l=0; l<mp->nLayers; l++)
					for (n=0; n<mp->vExcit[l]; n++)
						for (s=0; s<n_E[l][n].nLEff_E; s++)
							n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = gsl_ran_gaussian(mSeed, mp->SOMsigE); //states[th]
				break;
			case SOM: // distE, mp->DgElE, mp->SOMsigE
				phiScale = (mp->trainElE) ? mp->DgElE/(mp->SOMsigE*sqrt(2*M_PI)) : 1/(mp->SOMsigE*sqrt(2*M_PI));
				for (l=0; l<mp->nLayers; l++)
					for (n=0; n<mp->vExcit[l]; n++)
						for (s=0; s<n_E[l][n].nLEff_E; s++)
						{
							distance = readLowTriF(distE, l, n, n_E[l][n].lp0postsyn_E[s]->n); // distE[l][n][n_E[l][n].lp0postsyn_E[s]->n];
							n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = phiScale * exp(-pow(distance,2)/(2*pow(mp->SOMsigE,2)));
						}
				break;					
			default:
				exit_error("setWeights", "Illegal ElE weight initialisation arguement!");
				break;
		}
	}
	
	// Set remaining (non-plastic) weights
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vExcit[l]; n++)
			for (s=0; s<n_E[l][n].nLEff_I; s++)
				n_E[l][n].LEffs_I[s].delta_g = n_E[l][n].LEffs_I[s].delta_g_tm1 = mp->DgEI;
	
	/* I_ Synaptic weights */
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vInhib[l]; n++)
			for (s=0; s<n_I[l][n].nLEff_E; s++)
				n_I[l][n].LEffs_E[s].delta_g = n_I[l][n].LEffs_E[s].delta_g_tm1 = mp->DgIE;
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vInhib[l]; n++)
			for (s=0; s<n_I[l][n].nLEff_I; s++)
				n_I[l][n].LEffs_I[s].delta_g = n_I[l][n].LEffs_I[s].delta_g_tm1 = mp->DgII;
}


int loadImages(STIMULI * stim, PARAMS * mp)
{
	FILE * iList;
	FILEPARTS * fp;
	char * str, buff[BUFSIZ], filename[FNAMEBUFF]; //dir[DIRBUFF],
	int sLen = 0;
	int sc, or, ph;
	int count = 0;
	int nStimuli = 0;
	int * sCount = &nStimuli;
	int nTestStimuli = 0;
	int nTransPS = 0;
	int * tCount = &nTransPS;
	int nTestTransPS = 0;
	mp->newTestSet = stim->newTestSet = false;
	float ******* array = stim->trnImages;
	
	sLen = snprintf(filename, FNAMEBUFF, "%s/%s",mp->imgDir,mp->imgList);
	assert(sLen < FNAMEBUFF);
	iList = myfopen(filename, "r");
	fp = myalloc(sizeof(FILEPARTS));
	while ((str = fgets(buff, sizeof(buff), iList)) != NULL) 	/* Read next line */
	{
		if (str[0] == '\n' || str[0] == '#' || str[0] == '%') 	/* Skip blank lines and comments */
			continue;
		else if (str[0] == '*') // New set of object transforms
		{
			(*sCount)++;
			*tCount = 0;
			continue;
		}
		else if (str[0] == '+') // Seperate testing stimuli
		{
			mp->newTestSet = true;
			stim->newTestSet = true;
			sCount = &nTestStimuli;
			tCount = &nTestTransPS;
			*sCount = 0;
			*tCount = 0;
			stim->tstImages = get_7D_farray(mp->nTestStimuli, mp->nTestTransPS, mp->nScales, mp->nOrients, mp->nPhases, mp->nRows, mp->nCols, 0.0);
			array = stim->tstImages; //[st][tr][sc][or][ph][0], mp->nRows*mp->nCols
			continue;
		}
		else
		{
			trim(strtok(str,"%")); // Chop off comments and whitespace
			getFileParts(str, fp); //getFileParts((char*)trim(strtok(str,"%")), fp);
		}
		
		for (sc=0; sc<mp->nScales; sc++)
			for (or=0; or<mp->nOrients; or++)
				for (ph=0; ph<mp->nPhases; ph++)
					if (mp->gabor) // Use Gabor filter outputs
					{
						sLen = snprintf(filename, FNAMEBUFF, "%s/%s.flt/%s.%d.%d.%d.gbo", fp->path, fp->fname, fp->fname, \
										  mp->vScales[sc], mp->vOrients[or], mp->vPhases[ph]);
						assert(sLen < FNAMEBUFF);
						count += loadGaborOutput(filename, array[*sCount-1][*tCount][sc][or][ph][0], mp->nRows*mp->nCols);
						//print_farray(stderr, array[*sCount-1][*tCount][sc][or][ph], mp->nRows, mp->nCols);
					}
					else // Use output from Roger Watt's filtering software
					{
						sLen = snprintf(filename, FNAMEBUFF, "%s/%s.%s.filtered/%s.%s.%d.%d.%c", fp->path, fp->fname, fp->fext, fp->fname, fp->fext, \
										  mp->vScales[sc], mp->vOrients[or], (mp->vPhases[ph])?'p':'n');
						assert(sLen < FNAMEBUFF);
						count += loadDoGoutput(filename, array[*sCount-1][*tCount][sc][or][ph][0], mp->nRows*mp->nCols);
					}
		
		if (count == mp->nScales*mp->nOrients*mp->nPhases)
		{
			(*tCount)++;
			count = 0;
		}
		else
			exit_error("loadImages", "Wrong number of filters");
	}
	
	assert(mp->nStimuli == nStimuli);
	assert(mp->nTransPS == nTransPS);
	if (stim->newTestSet)
	{ // Test stimuli parameters must be read in from image parameter file
		assert(mp->nTestStimuli == nTestStimuli);
		assert(mp->nTestTransPS == nTestTransPS);
	}
	
	myfree(fp);
	fclose(iList);
	return 0;
}

int loadDoGoutput(const char * filename, float * array, int size)
{
	FILE * dgo;
	int e = 0;
	int err = 0;
	uchar * buff = myalloc(sizeof(uchar) * size);
	dgo = myfopen(filename, "r");
	fread(array, sizeof(uchar), size, dgo);
	assert(feof(dgo));
	err = fclose(dgo);
	for (e=0; e<size; e++)
		array[e] = mp->current + ((2.0 * buff[e] / 255.0) - 1.0) * mp->currentSpread; // Convert uchars [0,255] to floats
	// Scale for current injected too
	return (err) ? 0 : 1;
	/*if (!err)
		return 1;
	else
		return 0;*/
}

int loadGaborOutput(const char * filename, float * array, int size)
{
	FILE * gbo;
	int e = 0;
	int err = 0;
	gbo = myfopen(filename, "r");
	fread(array, sizeof(float), size, gbo);
	//assert(feof(gbo));
	err = fclose(gbo);
	for (e=0; e<size; e++) 	// Scale for current injected too
	{
		//assert(-1 <= array[e] && array[e] <= 1);
		array[e] = (array[e] * mp->currentSpread) + mp->current;
		if (array[e] < 0) // Check that negative currents are not injected
			array[e] = 0;
	}

	return (err) ? 0 : 1;
}

int loadGroups(STIMULI * stim, PARAMS * mp, char * filename)
{
	char * str, buff[32768]; //32KB//BUFSIZ, filename[FNAMEBUFF]; 
	char dummy[BUFSIZ]; // Can the character be disgarded without a dummy variable? See *
	//char * filename = "groupStimuli.dat";
	int n=0, g=0, count=0;
	//int nGroups = 0;
	//int * gCount = &nGroups;
	//int nTestGroups = 0;
	int nStimuli = 0, nTransPS = 0;
	int * sCount = &nStimuli;
	int * tCount = &nTransPS;
	int nTestStimuli = 0, nTestTransPS = 0;
	float *** array = NULL;

	char * delims = "^ \t";
	mp->nGroups = 0;
	stim->newTestSet = false; //bool newTestSet = false;
	
	FILE * sList = myfopen(filename, "r");
	
	// Load nGroups/nStimuli/nTransPS/sInputs from first line as a check?
	while ((str = fgets(buff, sizeof(buff), sList)) != NULL) 	/* Read next line */
	{
		if (str[0] == '\n' || str[0] == '#' || str[0] == '%') 	/* Skip blank lines and comments */
			continue;
		
		else if (str[0] == '@') // Metadata
		{
			if ((count = sscanf(str, "%c %d:%d:%d", dummy,&mp->nGroups,&nStimuli,&nTransPS)) != 4)
				exit_error("loadGroups", "Invalid file format");
			
#if DEBUG > 2
			printf("\n%s Metadata: nGroups=%d; nStimuli=%d; nTransPS=%d;\n",filename,mp->nGroups,nStimuli,nTransPS); fflush(stdout);
#endif
			stim->nStim = nStimuli; //nStimuliPG * mp->nGroups; // Reconsider ************************
			stim->nTrans = nTransPS;
			nStimuli = 0;
			nTransPS = 0;
			/*if (stim->nStim == 0)
			{
				stim->nStim = nStimuliPG * mp->nGroups;
				stim->nTrans = nTransPS;
			}
			else
			{
				assert(nStimuliPG * mp->nGroups == mp->nStimuli); //mp->nStimuli *= mp->nGroups;
				assert(nTransPS == mp->nTransPS);
				//nStimuliPG = 0;
				//nTransPS = 0;
			}*/
			
			stim->trn_stimuli = get_3D_farray(stim->nStim, stim->nTrans, mp->sInputs, 0.0);
			array = stim->trn_stimuli;
			stim->groups = (int **) array2d(mp->nGroups, mp->sInputs, sizeof(**(stim->groups))); // For labelling neurons //(bool **) //(mp->nTransPS>1)?mp->nCols:mp->sInputs
		}
		
		else if (str[0] == '^') // Prototype
		{
			str = strtok(buff, delims);
			for (n=0; n<mp->sInputs; n++) //((mp->nTransPS>1)?mp->nCols:mp->sInputs)
			{
				stim->groups[g][n] = atoi(str);
				str = strtok(NULL, delims);
			}
			assert(g < mp->nGroups);

#if DEBUG > 3
			printf("P%d:\t",g);
			print_irow(stdout, stim->groups[g], mp->sInputs);
#endif
			g++; //mp->nGroups++;
			continue;
		}
		/*else if (str[0] == '^') // New group
		{
			(*gCount)++;
			//if (*gCount > 0 && sCount != mp->nStimuli) // error
			*sCount = 0;
			*tCount = 0;
			continue;
		}*/
		
		else if (str[0] == '*') // New stimulus (set of object transforms)
		{
			if (*sCount > 0 && *tCount != ((stim->newTestSet) ? stim->nTestTrans : stim->nTrans))
				exit_error("loadGroups", "Wrong number of transforms in header");
			
			(*sCount)++;
			*tCount = 0;
			continue;
		}
		
		else if (str[0] == '+') // Seperate testing stimuli
		{
			if (*sCount != stim->nStim) // Check # Training stimuli
				exit_error("loadGroups", "Wrong number of stimuli in header");
			
			stim->newTestSet = true;
			//mp->newTestSet = true; // Present a novel stimulus from each simultaneously
			
			// Scan for nTestStimuli and nTestTrans
			if ((count = sscanf(str, "%c %d:%d", dummy,&nTestStimuli,&nTestTransPS)) != 3)
				exit_error("loadGroups", "Invalid file format (testing)");

			/*assert(*sCount == mp->nTestStimuli);
			assert(*tCount == mp->nTestTransPS);
			stim->tst_stimuli = get_3D_farray(mp->nTestStimuli, mp->nTestTransPS, mp->sInputs, 0.0);*/
			stim->nTestStim = nTestStimuli;
			stim->nTestTrans = nTestTransPS;
			stim->tst_stimuli = get_3D_farray(stim->nTestStim, stim->nTestTrans, mp->sInputs, 0.0);

			// Reset counters and populate Testing stimuli array
			array = stim->tst_stimuli;
			//gCount = &nTestGroups;
			sCount = &nTestStimuli;
			tCount = &nTestTransPS;
			*sCount = 0;
			*tCount = 0;
			continue;
		}
		
		else // Load the transform into stimulus array
		{
			(*tCount)++;
			//count=0;
			// ASCII
			str = strtok(buff, delims);
			for (n=0; n<mp->sInputs; n++)
			{
				/*//if (!str)
				//	break;
				//fprintf(stderr,"n=%d; %s;\n", n,str);
				//fprintf(stderr,"n=%d; str=%s;\tgCount=%d; sCount=%d; tCount=%d;\n",n,str,*gCount,*sCount,*tCount); fflush(stderr);
				//stim->trnGrpStim[*gCount-1][*sCount-1][*tCount-1][n] *= atof(str) * mp->current;
				//count += sscanf(str, "%f", &stim->trnGrpStim[*gCount-1][*sCount-1][*tCount-1][n]);
				//stim->trnGrpStim[*gCount-1][*sCount-1][*tCount-1][n] *= mp->current;*/
				array[*sCount-1][*tCount-1][n] = atof(str) * mp->current;
				str = strtok(NULL, delims); // Tokenise the string to get the next value
			}
			
			/*if (count != mp->sInputs)
				exit_error("loadGroups", "Wrong size inputs");*/
			
			/* Binary
			//loadStimuli(filename, stim->trnGrpStim[*gCount-1][*sCount][*tCount], mp->sInputs);
			//fread(stim->trnGrpStim[g][s][t][0], sizeof(float), mp->sInputs, sList);*/
			
#if DEBUG > 3
			printf("%sS%dT%d:\t",(stim->newTestSet)?"Test\t":"",*sCount-1,*tCount-1);
			print_frow(stdout, array[*sCount-1][*tCount-1], mp->sInputs);
#endif
			
			continue;
		}


	}
	fclose(sList);
	
	if (g != mp->nGroups)
		exit_error("loadGroups", "Inconsistent numbers of groups");
	
	if (!stim->newTestSet)
	{
		stim->nTestStim = stim->nStim;
		stim->nTestTrans = stim->nTrans;
	}
	
	/*if (*sCount * mp->nGroups != mp->nStimuli)
		exit_error("loadGroups", "Wrong number of stimuli in header");
	*sCount = 0;
	
	if (*tCount != mp->nTransPS)
		exit_error("loadGroups", "Wrong number of transforms in header");
	*tCount = 0;*/
	
	return 0;
}

/*int loadStimuli(const char * filename, float * array, int size)
{
	FILE * stim;
	int e=0;
	bool err;
	if (file_exists("stimuli.tbz"))
		system("tar -xvf stimuli.tbz");
	stim = myfopen(filename, "rb");
	fread(array, sizeof(float), size, stim);
	for (e=0; e<size; e++)
		array[e] = array[e] * mp->current;
	err = fclose(stim);
	return (err) ? 0 : 1;
}*/


void genGroups(STIMULI * stim, PARAMS * mp)
{
	// Generate stimuli from prototypes rather than load from files (generate new stimuli with each seed)
	int g=0, s=0, n=0;
	
	// Reconsider sInputs (nCols) and nFiringNeurons (nSP)
	// Consider translating stimuli and PP stimuli
	
	// Error checking to see that constraints are satisfiable
	assert(mp->sInputs >= mp->nBG + (mp->nGroups * mp->nWG));
	
	// Allocate stimulus arrays and set parameters in stim structure
	stim->trn_stimuli = get_3D_farray(mp->nStimuli, mp->nTransPS, mp->sInputs, 0.0);
	stim->tst_stimuli = get_3D_farray(mp->nTestStimuli, mp->nTestTransPS, mp->sInputs, 0.0);
	stim->nStim = mp->nStimuli;
	stim->nTrans = mp->nTransPS;
	stim->nTestStim = mp->nTestStimuli;
	stim->nTestTrans = mp->nTestTransPS;
 
	stim->newTestSet = true;
	
	// Create prototypes
	// mp->layDim[0].nCols or mp->sInputs if not translating
	stim->groups = (int **) array2d(mp->nGroups, mp->sInputs, sizeof(**(stim->groups))); // change to bool
	
	int * bag = myalloc(mp->sInputs * sizeof(*bag)); // Zeroed in function
	for (n=0; n<mp->sInputs; n++)
		bag[n] = n;
	gsl_ran_shuffle(mSeed, bag, mp->sInputs, sizeof(*bag));
	
	int nGPool = mp->nBG + mp->nWG; // Group pool size
	assert(mp->nFiringNeurons <= nGPool);
	int ** gPools = (int **) array2d(mp->nGroups, nGPool, sizeof(**gPools));
	
	// Set common elements for all prototypes (nBG)
	for (n=0; n<mp->nBG; n++) // Skipped if nBG==0
		for (g=0; g<mp->nGroups; g++)
		{
			stim->groups[g][bag[n]] = 1;
			gPools[g][n] = bag[n];
		}
			
	int bStart = mp->nBG; //==n // Bag offset to prevent reuse of assigned neurons
	

	// Set the remaining group specific elements for each group (nWG)
	for (g=0; g<mp->nGroups; g++)
	{
		for (n=0; n<mp->nWG; n++)
		{
			stim->groups[g][bag[n+bStart]] = 1;
			gPools[g][n+bStart] = bag[n+bStart];
		}
		bStart += mp->nWG; // Update the bag offset		
	}
	myfree(bag);

	
	// Generate individual training exemplars
	int nStimPG = mp->nStimuli / mp->nGroups;
	assert(mp->nStimuli % mp->nGroups == 0); //assert(nStimPG * mp->nGroups == mp->nStimuli);
	//int * exemplar = myalloc(mp->nFiringNeurons * sizeof(*exemplar));
	for (g=0; g<mp->nGroups; g++)
		for (s=0; s<nStimPG; s++)
		{
			gsl_ran_shuffle(mSeed, gPools[g], nGPool, sizeof(**gPools));
			for (n=0; n<mp->nFiringNeurons; n++)
				stim->trn_stimuli[s][0][gPools[g][n]] = mp->current;
			
			//gsl_ran_choose(mSeed, exemplar, mp->nFiringNeurons, gPools[g], nGPool, sizeof(*exemplar));
			//for (n=0; n<mp->nFiringNeurons; n++)
			//	stim->trn_stimuli[s][0][exemplar[n]] = mp->current;
		}
	
	
	// Present all groups simultaneously in novel test stimuli
	//int * shuffle = myalloc(nGPool * sizeof(*shuffle));
	//for (n=0; n<nGPool; n++)
	//	shuffle[n] = n;
	for (s=0; s<mp->nTestStimuli; s++)
	{
		for (g=0; g<mp->nGroups; g++)
		{
			gsl_ran_shuffle(mSeed, gPools[g], nGPool, sizeof(**gPools));
			for (n=0; n<mp->nFiringNeurons; n++)
				stim->tst_stimuli[s][0][gPools[g][n]] = mp->current;
			
			//gsl_ran_shuffle(mSeed, shuffle, nGPool, sizeof(*shuffle));
			//stim->tst_stimuli[s][0][gPools[g][shuffle[n]]] = mp->current;
		}
	}
	
	myfree(gPools);
	//myfree(shuffle);
	//myfree(exemplar);
	
}


void printGroups(STIMULI * stim, PARAMS * mp, const char * filename) // Merge with printStimuli
{
	FILE * stimFP = NULL;
	int g=0, s=0, n=0, spg=0;
	
	
	if (mp->priorPhases) // Move outside?
	{
		
	}
	
	stimFP = myfopen(filename, "w");
	fprintf(stimFP, "@ %d:%d:%d\n", mp->nGroups, stim->nStim, stim->nTrans); // Print Metadata
	for (g=0; g<mp->nGroups; g++)
	{
		fprintf(stimFP, "^\t"); // Print Prototypes
		for (n=0; n<mp->sInputs; n++) // Change to printIrow()
			fprintf(stimFP, "%d ", stim->groups[g][n]);
		fprintf(stimFP, "\n");
	}
	
	int nStimPG = stim->nStim / mp->nGroups;
	for (s=0, g=0; g<mp->nGroups; g++)
	{
		for (spg=0; spg<nStimPG; spg++)
		{
			fprintf(stimFP, "* Stimulus %d (g#%d,s#%d)\n", ++s, g, spg);
			if (stim->nTrans > 1)
			{
				
			}
			else
			{
				fprintf(stimFP, "\t");
				for (n=0; n<mp->sInputs; n++) // Change to printIrow()
					fprintf(stimFP, "%d ", stim->trn_stimuli[s][0][n]?1:0); //%1.0f with ceil()
				fprintf(stimFP, "\n");
			}
		}
		
	}
	if (stim->newTestSet)
	{
	fprintf(stimFP, "+ %d:%d Testing Set\n", stim->nTestStim, stim->nTestTrans);
    for (s=0; s<stim->nTestStim; s++)
	{
		fprintf(stimFP,"* Test Stimulus %d\n", s+1);
		if (stim->nTestTrans > 1)
		{
			
		}
		else
		{
			fprintf(stimFP, "\t");
			for (n=0; n<mp->sInputs; n++) // Change to printIrow()
				fprintf(stimFP, "%d ", stim->tst_stimuli[s][0][n]?1:0);
			fprintf(stimFP, "\n");
		}
	}
    }
 
	fclose(stimFP);
	
	
	// Print prototypes for neuron labelling
	stimFP = myfopen("prototypes.stm", "w");
	print_iarray(stimFP, stim->groups, mp->nGroups, mp->sInputs);
	fclose(stimFP);
}


void gen_stimuli(bool rep, STIMULI * stim, PARAMS * mp)
{
	// Assumes 1D inputs
	// Place array arguements in a patterns structure
	int trans, n, p, block=0, slen=0;
	int * choices = NULL;
	int * chosen = NULL;
	char stimStr[BUFSIZ];

	/* Generate the training stimuli to present to the network */
	// Patterns could be generated in Matlab and loaded from dat files... or read in list of pairs from another array
	
	stim->tst_stimuli = get_3D_farray(mp->nTestStimuli, mp->nTestTransPS, mp->sInputs, 0.0);
	
	// Generate test stimuli
	switch (rep) 
	{
		case 0: /* Distributed Patterns */ // *** Use prototypes?
			trans=0;
			choices = myalloc(mp->sInputs * sizeof(*choices));
			for (n=0; n<mp->sInputs; n++)
				choices[n] = n;
			chosen = myalloc(mp->nFiringNeurons * sizeof(*chosen));
			for (p=0; p<mp->nTestStimuli; p++)
			{
				gsl_ran_choose(mSeed, chosen, mp->nFiringNeurons, choices, mp->sInputs, sizeof(*choices));
				for (n=0; n<mp->nFiringNeurons; n++)
					stim->tst_stimuli[p][trans][chosen[n]] = mp->current;
			}
			myfree(choices);
			myfree(chosen);
			break;
			
		case 1: /* Local Patterns */
			block = floor(mp->nFiringNeurons + ((mp->nTransPS - 1) * mp->shift));
			// This assumes stimuli are a contiguous block of 1's
			for (p=0; p<mp->nTestStimuli; p++)
				for (trans=0; trans<mp->nTransPS; trans++)
					for (n=(p*block)+(trans*mp->shift); n<(mp->nFiringNeurons+(p*block)+(trans*mp->shift)); n++)
						stim->tst_stimuli[p][trans][n] = mp->current;
			break;
			/*	***111111111222222222222333333333333
				111***111111222222222222333333333333
				...
				------------***222222222------------
				------------222***222222------------
				...
				------------------------333333333***
				*/ // Extend this to 2D - cf MSc code
	}
	
	// Generate training stimuli
	if (mp->K > 1) // Present multiple stimuli simultaneously
	{
		assert(mp->K <= mp->nTestStimuli);
		if (mp->K == mp->nTestStimuli) // mp->nStimuli = 1
			assert(mp->M == 1);
		
		/*if (mp->K == mp->nTestStimuli) // Try this!
		{
			assert(mp->M == 1);
			mp->nStimuli = 1;
			stim->trn_stimuli = get_3D_farray(mp->nStimuli, mp->nTransPS, mp->sInputs, 0.0);
			memcpy(stim->trn_stimuli[0][0], stim->tst_stimuli[0][0], mp->nTestTransPS*mp->sInputs*sizeof(stim->tst_stimuli[0][0][0]));
			
			for (c=0; c<mp->K-1; c++)
				for (trans=0; trans<mp->nTestTransPS; trans++) // memcpy pth test stimulus ? - does not work for dist. stim
					for (n=0; n<mp->sInputs; n++)
						stim->trn_stimuli[0][trans][n] = (stim->tst_stimuli[c][trans][n]) ? mp->current : stim->trn_stimuli[0][trans][n];			
		}
		else ...*/
		
		/*if (abs(mp->nStimuli - mp->K)==1) // || mp->K == 1
		{
			mp->nStimuli = mp->nTestStimuli
		}
		else
			// Calculate nCombs and select a subset according to mp->M
		 */
		stim->trn_stimuli = get_3D_farray(mp->nStimuli, mp->nTransPS, mp->sInputs, 0.0); // rethink nStimuli for limited training...
		choices = myalloc((mp->nTestStimuli-1) * sizeof(*choices));
		chosen = myalloc((mp->K-1) * sizeof(*chosen));
		int ind=0, c=0, m=0, q=0;
		fprintf(stdout, "\n");
		for (p=0; p<mp->nTestStimuli; p++)
		{
			for (c=0, ind=0; c<mp->nTestStimuli-1; c++, ind++)
				choices[c] = (c==p) ? ++ind : ind;

			for (m=0; m<mp->M; m++) // This method could have duplicate combinations
			{
				gsl_ran_choose(mSeed, chosen, mp->K-1, choices, mp->nTestStimuli-1, sizeof(*choices));
				slen = snprintf(stimStr, BUFSIZ, "\tTraining stimulus %d (inc. %d)", p, p);
				assert(slen < BUFSIZ);
				printIntArray(stdout, stimStr, chosen, mp->K-1);
				
				q = (p * mp->M) + m;
				memcpy(stim->trn_stimuli[q][0], stim->tst_stimuli[p][0], mp->nTestTransPS*mp->sInputs*sizeof(stim->tst_stimuli[p][0][0]));
				
				for (c=0; c<mp->K-1; c++)
					for (trans=0; trans<mp->nTestTransPS; trans++) // memcpy pth test stimulus ? - does not work for dist. stim
						for (n=0; n<mp->sInputs; n++)
							stim->trn_stimuli[q][trans][n] = (stim->tst_stimuli[chosen[c]][trans][n]) ? mp->current : stim->trn_stimuli[q][trans][n];
			}
		}
		
		myfree(choices);
		myfree(chosen);
		
		stim->nStim = mp->nStimuli;
		stim->nTrans = mp->nTransPS;
		stim->nTestStim = mp->nTestStimuli;
		stim->nTestTrans = mp->nTestTransPS;
	}
	else	// Training stimuli presented individually
	{
		stim->trn_stimuli = stim->tst_stimuli;
		stim->nStim = stim->nTestStim = mp->nTestStimuli;
		stim->nTrans = stim->nTestTrans = mp->nTestTransPS;
	}
	
	stim->newTestSet = mp->newTestSet;
	
	// Generate all ways of choosing K from N: 
	// http://phoxis.org/2009/10/13/allcombgen/ 
	// http://compprog.wordpress.com/2007/10/17/generating-combinations-1/ 
	// http://www.cs.utexas.edu/users/djimenez/utsa/cs3343/lecture25.html
	// Select M of those ways
	
	/*if (mp->M) // Generalise to pair with K partner stimuli
	{
		assert(mp->K <= mp->nTestStimuli);
		assert(mp->M < mp->nTestStimuli);
		stim->trn_stimuli = get_3D_farray(mp->nStimuli, mp->nTransPS, mp->sInputs, 0.0);

		// Keep track of combinations otherwise there will be duplicate training patterns (violating calculated nTestStimuli above)
		
		choices = myalloc(mp->nTestStimuli-1 * sizeof(*choices));
		chosen = myalloc(mp->M * sizeof(*chosen));
		int c = 0;
		int ind = 0;
		for (p=0; p<mp->nTestStimuli; p++)
		{
			ind = 0;
			for (c=0; c<mp->nTestStimuli-1; c++)
			{
				if (c == p)
					ind++;
				choices[c] = ind++;
			}
			gsl_ran_choose(mSeed, chosen, (mp->M>mp->nTestStimuli-1)?mp->nTestStimuli-1:mp->M, choices, mp->nTestStimuli-1, sizeof(*choices)); // Choose M (up to N-1) from N
			gsl_ran_shuffle(mSeed, chosen, (mp->M>mp->nTestStimuli-1)?mp->nTestStimuli-1:mp->M, sizeof(*chosen)); // Permute partners
			for (c=0; c<mp->M; c++)
				for (trans=0; trans<mp->nTestTransPS; trans++)
					for (n=0; n<mp->sInputs; n++) // Redo for K stimuli below
						stim->trn_stimuli[p][trans][n] = (stim->tst_stimuli[chosen[c]][trans][n] || stim->tst_stimuli[p][trans][n]) ? mp->current : 0.0;
		}
		myfree(choices);
		myfree(chosen);
	}
	else
	{
		stim->trn_stimuli = stim->tst_stimuli;
	}
	/////*/
	
#if DEBUG>3 // Level 4
	fprintf(stderr, "\nPrinting generated training stimuli...\n");
	for (p=0; p<mp->nStimuli; p++)
		for (trans=0; trans<mp->nTransPS; trans++)
		{
			fprintf(stderr, "S%dT%d: ",p+1,trans+1);
			print_frow(stderr, stim->trn_stimuli[p][trans], mp->sInputs);
		}
#endif
	
	/*** Testing stimuli ***/
	// Modify to generate single stimuli test patterns when there are multi-stimulus training patterns
	//memcpy(*tst_stimuli, *trn_stimuli, sizeof(*tst_stimuli));
	
	return; // void;
}


void printStimuli(STIMULI * stim, PARAMS * mp)
{
	FILE * stFP = NULL;
	int p=0;
	
	
	stFP = myfopen("trn_stimuli.dat", "w");
	for (p=0; p<stim->nStim; p++)
	{
		fprintf(stFP, "*** Stimulus #%d (%d/%d) ***\n", p, p+1, stim->nStim);
		print_farray(stFP, (stim->trn_stimuli)[p], stim->nTrans, mp->sInputs);
		fprintf(stFP, "\n");
	}
	fclose(stFP);
	
	if (stim->newTestSet)
	{
		stFP = myfopen("tst_stimuli.dat", "w");
		for (p=0; p<stim->nTestStim; p++)
		{
			fprintf(stFP, "*** Stimulus #%d (%d/%d) ***\n", p, p+1, stim->nTestStim);
			print_farray(stFP, (stim->tst_stimuli)[p], stim->nTestTrans, mp->sInputs);
			fprintf(stFP, "\n");
		}
		fclose(stFP);
	}
	
	// Matlab friendly stimuli output
	
	// Change according to input layer dimensions i.e. print out 1D, 2D or 7D(?) patterns
	char label[BUFSIZ];
	int t=0;
	stFP = myfopen("stimuli.m", "w");
	fprintf(stFP, "%% *** Training Stimuli ***\n");
	fprintf(stFP, "STIM.train = cell(%d,%d);\n",stim->nStim,stim->nTrans);
	for (p=0; p<stim->nStim; p++)
	{
		for (t=0; t<stim->nTrans; t++)
		{
			//fprintf(stFP, "STIM.train{%d,%d} =", p, p+1, stim->nStim);
			snprintf(label, BUFSIZ, "STIM.train{%d,%d}",p+1,t+1);
			printFloatArray(stFP, label, (stim->trn_stimuli)[p][t], mp->sInputs);
		}
	}
	
	if (stim->newTestSet)
	{
		fprintf(stFP, "\n%% *** Testing Stimuli ***\n");
		fprintf(stFP, "STIM.test = cell(%d,%d);\n",stim->nTestStim,stim->nTestTrans);
		for (p=0; p<stim->nTestStim; p++)
		{
			for (t=0; t<stim->nTestTrans; t++)
			{
				snprintf(label, BUFSIZ, "STIM.test{%d,%d}",p+1,t+1);
				printFloatArray(stFP, label, (stim->tst_stimuli)[p][t], mp->sInputs);
			}
		}
	}
	
	fprintf(stFP, "\n%% *** Training schedule ***\n");
	fprintf(stFP, "%% Transforms are presented sequentially during (pre)testing.\n");
	fprintf(stFP, "STIM.schedule = cell(%d,1);\n", mp->loops);
	
	fclose(stFP);
	return;
}


int genShuffles(STIMULI * stim, PARAMS * mp)
{
	int loop=0, p=0, trans=0, count=0;
	int * choices = NULL;
	
	if (mp->randStimOrder)
	{
		choices = myalloc(stim->nStim * sizeof(int));
		for (p=0; p<stim->nStim; p++)
			choices[p] = p;
		stim->stimShuffle = get_2D_iarray(mp->loops, stim->nStim, 0);
		for (loop=0; loop<mp->loops; loop++)
		{
			gsl_ran_shuffle(mSeed, choices, stim->nStim, sizeof(int));
			memcpy(stim->stimShuffle[loop], choices, stim->nStim * sizeof(int));
		}
		count++;
		myfree(choices);
	}
	
	// Currently assumes the same number of transforms for all stimuli
	if (mp->randTransOrder)
	{
		choices = myalloc(stim->nTrans * sizeof(int));
		for (trans=0; trans<stim->nTrans; trans++)
			choices[trans] = trans;
		stim->transShuffle = get_3D_iarray(mp->loops, stim->nStim, stim->nTrans, 0);
		for (loop=0; loop<mp->loops; loop++)
		{
			for (p=0; p<stim->nStim; p++)
			{
				gsl_ran_shuffle(mSeed, choices, stim->nTrans, sizeof(int));
				memcpy(stim->transShuffle[loop][p], choices, stim->nTrans * sizeof(int));
			}
		}
		count++;
		myfree(choices);
	}
	
	/*if (mp->randTransDirection) // Mutually exclusive with randTransOrder
	{
		stim->transShuffle = get_3D_iarray(mp->loops, stim->nStim, stim->nTrans, 0);
		for (loop=0; loop<mp->loops; loop++)
		{
			for (p=0; p<stim->nStim; p++)
			{
				reverse = (gsl_rng_uniform_int(mSeed, 2)) ? true : false;
				for (trans=0; trans<stim->nTrans; trans++)
					stim->transShuffle[loop][p][trans] = reverse ? (stim->nTrans-1)-trans : trans;
			}
		}
	}*/
	
	return count;
}


/*void genSchedule(STIMULI * stim, PARAMS * mp) // if (mp->train)
{
	//size_t index_size = (sizeof(***stim->sched) * mp->loops) + (sizeof(**stim->sched) * mp->loops * stim->nStim);
	size_t index_size = (sizeof(***stim->sched) * mp->loops) * (1 + stim->nStim);
    size_t store_size = sizeof(SCHEDULE) * mp->loops * stim->nStim * stim->nTrans;
	
    stim->sched = myalloc(index_size + store_size);
    //if(!a) return NULL;
	
    //memset(stim->sched + index_size, 0, store_size); // Be careful with memsets rezeroing the array
	size_t l=0, s=0;
	for (l=0; l<mp->loops; l++)
		for (s=0; s<stim->nStim; s++)
			stim->sched[l][s] = index_size + (l*stim->nStim*stim->nTrans) + (s*stim->nTrans);
			//((void **)a)[i] = a + index_size + i * cols * value_size;
	
    //return (void **)a;
	
	size_t t=0;
	sChoices = myalloc(stim->nStim * sizeof(int));
	for (s=0; s<stim->nStim; s++)
		sChoices[s] = s;
	
	tChoices = myalloc(stim->nTrans * sizeof(int));
	for (t=0; t<stim->nTrans; t++)
		choices[t] = t;
	
	for (l=0; l<mp->loops; l++)
	{
		if (mp->randStimOrder)
		{
			gsl_ran_shuffle(mSeed, sChoices, stim->nStim, sizeof(int));
		}
		else if (mp->interleaveTrans)
		{
			for (t=0; t<stim->nTrans; t++)
			{
				for (s=0; s<stim->nStim; s++)
				{
					stim->sched[l][s][t]
				}
			}
		}

		for (s=0; s<stim->nStim; s++)
		{
			if (mp->randTransOrder)
			{
				gsl_ran_shuffle(mSeed, tChoices, stim->nTrans, sizeof(int));
			}
			else if (mp->randTransDirection)
			{
				reverse = (gsl_rng_uniform_int(mSeed, 2)) ? true : false;
			}
			
			for (t=0; t<stim->nTrans; t++)
			{
				stim->sched[l][s][t].st = sChoices[s];
				stim->sched[l][s][t].tr = (reverse) ? (stim->nTrans-1)-t : tChoices[t];
			}
		}
	}
	
	return;	
}*/


void printSchedule(STIMULI * stim, const char * filename)
{
	int l=0, p=0, slen=0;
	char label[BUFSIZ];
	FILE * fp = myfopen(filename, "w");
	fprintf(fp, "STIM.stimShuffle = cell(MP.loops,1);\n"); 
	fprintf(fp, "STIM.transShuffle = cell(MP.loops,%d);\n", stim->nStim);
	for (l=0; l<mp->loops; l++)
	{
		slen = snprintf(label, BUFSIZ, "stimShuffle{%d}", l+1);
		assert(slen < BUFSIZ);
		printIntArray(fp, label, stim->stimShuffle[l], stim->nStim);
		for (p=0; p<stim->nStim; p++)
		{
			slen = snprintf(label, BUFSIZ, "transShuffle{%d,%d}", l+1, p+1);
			assert(slen < BUFSIZ);
			printIntArray(fp, label, stim->transShuffle[l][stim->stimShuffle[l][p]], stim->nTrans);
		}
	}
	fclose(fp);
}


void calcInput(PARAMS * mp, int loop, int pat, int trans, STIMULI * stim, float ** input, int regime) // return int * input?
{	// Assumes all stimuli translate

	switch (regime)
	{
		case Testing: // Testing stimuli
			if (mp->useFilteredImages)
				*input = ****(stim->tstImages[pat][trans]);
			else
				*input = stim->tst_stimuli[pat][trans];
			//*input = (mp->useFilteredImages) ? ****(stim->tstImages[pat][trans]) : stim->tst_stimuli[pat][trans];
			break;
			
		case Training: // Training stimuli
			pat = (mp->randStimOrder) ? stim->stimShuffle[loop][pat] : pat;
			trans = (mp->randTransOrder) ? stim->transShuffle[loop][pat][trans] : trans;

			*input = (mp->useFilteredImages) ? ****(stim->trnImages[pat][trans]) : stim->trn_stimuli[pat][trans];
			
			FILE * fp = myfopen("stimuli.m", "a+"); //schedule.m
			//fprintf(fp, "%d %d %d\n", loop, pat, trans);
			fprintf(fp, "%d,%d; ",pat, trans);
			fclose(fp);
			
			break;
			
		default:
			exit_error("calcInput", "Unknown regime");
			break;
	}
	
	
#if DEBUG > 3
	print_frow(stderr, *input, mp->sInputs);
#endif
	
	return; // void;
}

void simulatePhase(LEARNREGIME regime, const char * prefix, STIMULI * stim)
{
	int * o = NULL;
	int * i = NULL;
	int oCount=0, iCount=0;
	int p=0, tr=0; //g=0
	int nStimuli=0, nTrans=0; //nGroups,
	float transP = 0.0;
	tstep t_start = 0, t_end = 0;
	int loop, l, wl, n, syn;
	int nLoops = 0;
	int result = 0;
	//LEARNREGIME regime = Continuous; //NoLearning;
	int slen = 0;
	//char phaseString[FNAMEBUFF];
	char stimFile[FNAMEBUFF];
    char filename[FNAMEBUFF];
	char buff[BUFSIZ];
	FILE * excitOutput;
	FILE * inhibOutput;
	FILE * weightOutput;
	FILE * recFile;
	//char fullprefix[FNAMEBUFF];
	bool reverse = false;
	double percentage = 0.0;
	float * input = NULL; //myalloc(mp->sInputs * sizeof(*input));
#if DEBUG > 1 // Level 2
#ifdef _OPENMP
	double secs = 0.0;
	char timeStr[BUFSIZ];
	char remainStr[BUFSIZ];
#endif
#endif
	
	
	switch (regime) //(sPhase) 
	{
		case Training:
			nStimuli = stim->nStim;
			nTrans = stim->nTrans;
			transP = mp->transP_Train;
			nLoops = mp->loops;
			break;
			
		case Testing:
			nStimuli = stim->nTestStim;
			nTrans = stim->nTestTrans;
			transP = mp->transP_Test;
			nLoops = 1;
			break;
			
		default:
			exit_error("simulatePhase", "Unknown phase!");
			break;
	}
	
	if (mp->interleaveTrans && regime == Training)
	{
		o = &tr;
		i = &p;
		oCount = nTrans;
		iCount = nStimuli;
	}
	else
	{
		o = &p;
		i = &tr;
		oCount = nStimuli;
		iCount = nTrans;
	}
	
	//printf("\tNow beginning %s phase...\n",phaseString); // *** Move outside
	if (regime == Training) //(sPhase == Training)
		printf("\tSimulating %2.3f s per loop for %d loop%s.\n",mp->EpochTime,nLoops,(nLoops>1)?"s":"");
	for (loop=0; loop<nLoops; loop++)
	{
		initNetwork(Hard);//(NoLearning); // Reset V's, g's, C's, D's and spike buffers - NoLearning even for Training! 
		if (regime == Training) //(sPhase == Training)
		{
			printf("\tLoop #%d (%d/%d)\n", loop, loop+1, mp->loops);
			
			slen = snprintf(stimFile, FNAMEBUFF, "%sstimuli.m", prefix);
			assert(slen < FNAMEBUFF);
			slen = snprintf(buff, BUFSIZ, "STIM.schedule{%d} = [", loop+1);
			assert(slen < FNAMEBUFF);
			append(stimFile, buff);
		}
			
		/*for (g=0; g<((mp->stimGroups)?nGroups:1); g++)
		{
			if (g>nGroups)
				getchar();
			if (mp->stimGroups)
				printf("\tGroup %d/%d\n", g+1, nGroups);*/
			for (*o=0; *o<oCount; (*o)++)
			{
				if (regime==Training)
				{
					if (mp->trainPause && !mp->interleaveTrans)
						initNetwork(Soft);
					if (mp->randTransDirection) // Generate [0,n-1] with equal probability
						reverse = (gsl_rng_uniform_int(mSeed, 2)) ? true : false;
				}
				
				for (*i=0; *i<iCount; (*i)++)
				{
					if (regime == Testing) //if (sPhase != Training) //Testing only
						initNetwork(Soft);
					calcInput(mp, loop, p, ((reverse) ? (nTrans-1)-tr : tr), stim, &input, regime);
					if ( tr==0 || (mp->interleaveTrans && regime) )
						printf("\t\tPresenting stimulus %d/%d...\n", (p+1), nStimuli);
					if (nTrans > 1)
						printf("\t\t\tTransform %d/%d...\n", (reverse ? nTrans-tr : tr+1), nTrans); //\r
					t_start = round((*i + (*o * iCount)) * transP * ceil(1/mp->DT));
					t_end = t_start + round(transP * ceil(1/mp->DT));


					
#if DEBUG > 1 // Level 2
					fprintf(stderr, "Updating network from timestep %d to %d.\n", t_start, t_end-1); //%lld
#ifdef _OPENMP
					SIM.elapsed = omp_get_wtime() - SIM.start;
					getTimeString(timeStr, BUFSIZ, SIM.elapsed, "ms");
					printf("[%s]",timeStr); // Time stamp
#endif
#endif

					updateNetwork(t_start, t_end, input, regime); //loop, 
					
                    // Update to normalise ElE weights and skip if not Training
					if (mp->normalise) // Normalise plastic weights
						result = normalise(n_E, mp);
					
					/*if (mp->saveInputSpikes)
					 {
					 // Print spikes to file (see below) ...
					 slen = snprintf(filename, FNAMEBUFF, "E%dL0ExcitSpikes.dat", loop);
					 
					 for (n=0; n<mp->sInputs; n++)
					 {
					 n_E[0][n].spkbin = 0;
					 memset(n_E[0][n].spikeTimes, 0, mp->inpSpkBuff*sizeof(n_E[0][n].spikeTimes[0])); //[0]?
					 }
					 }*/
					
					SIM.tally += round(transP/mp->DT); //transP * ceil(1/mp->DT);
					percentage = (100.0*SIM.tally)/SIM.totTS;
#if DEBUG > 1 // Level 2
#ifdef _OPENMP
					secs = omp_get_wtime() - SIM.start - SIM.elapsed;
					SIM.realSecPerSimSec = (secs) / transP;
					getTimeString(remainStr, BUFSIZ, SIM.realSecPerSimSec * (SIM.totTS-SIM.tally) * mp->DT, "s");
					// Format percentage to 3 s.f.
					int width, precision;
					double logv = log10(percentage);
					if (logv < 1)
					{
						if (logv < 0) // 3 d.p.
						{
							precision = 3;
							width = 5; // 0.xxx
						}
						else //if (logv >= 0 && logv < 1) // 2 d.p. x.xx
						{
							precision = 2;
							width = 4;
						}
					}
					else
					{
						if (logv < 2) //if (logv >= 1 && logv < 2) // 1 d.p. xx.x
						{
							precision = 1;
							width = 4;
						}
						else // 0 d.p. 
						{
							precision = 0;
							width = floor(logv) + 1;
						}
					}

					printf("\t%G/s\t%*.*lf%%\tEst. time remaining: %s\n",round(SIM.realSecPerSimSec),width,precision,percentage,remainStr);
#endif
#endif
					if (SIM.Xgrid)
					{
						//printf("\n");
						printf("<xgrid>{control = statusUpdate; percentDone = %.0lf; }</xgrid>\n", percentage);
						fflush(stdout);
					}
				}
			}
		//}
		//printf("\t%s complete!\n",phaseString); // *** Move outside
		
		if (regime == Training)
		{
			//slen = snprintf(filename, FNAMEBUFF, "%sstimuli.m", prefix); // filename not changed yet
			//assert(slen < FNAMEBUFF);
			//slen = snprintf(buff, BUFSIZ, "schedule{%d} = [", loop+1);
			//assert(slen < FNAMEBUFF);
			append(stimFile, "]';\n"); // schedule(1,:) := stimulus; schedule(2,:) := transform;
		}
		
		printf("\tSaving results..."); // Output results to dat files
		
		for (l=0; l<mp->nLayers; l++) // Save Excitatory spikes 
		{
			if (regime == Training)
				slen = snprintf(filename, FNAMEBUFF, "%sE%dL%dExcitSpikes.dat", prefix, loop, l);
			else // <pre>Testing
				slen = snprintf(filename, FNAMEBUFF, "%sL%dExcitSpikes.dat", prefix, l);
			assert(slen < FNAMEBUFF);
			excitOutput = myfopen(filename, "w");
			for (n=0; n<mp->vExcit[l]; n++)
			{
				fprintf(excitOutput, "%d\t", n_E[l][n].spkbin); // Print number of spikes first (in case any at t=0)
				print_irow(excitOutput, (int*) n_E[l][n].spikeTimes, n_E[l][n].spkbin); // spikeTimes[0] = -BIG?
			}
			fclose(excitOutput);
		}
		
		if (regime == Testing) // Save Inhibitory spikes
		{
			for (l=0; l<mp->nLayers; l++)
			{
				slen = snprintf(filename, FNAMEBUFF, "%sL%dInhibSpikes.dat", prefix, l);
				assert(slen < FNAMEBUFF);
				inhibOutput = myfopen(filename, "w"); // Make either HR "w" or Binary "wb"
				for (n=0; n<mp->vInhib[l]; n++)
				{
					fprintf(inhibOutput, "%d\t", n_I[l][n].spkbin); // Print number of spikes first (in case any at t=0)
					print_irow(inhibOutput, (int*) n_I[l][n].spikeTimes, n_I[l][n].spkbin);
				}
				fclose(inhibOutput);
			}
		}
		
		if (regime != Training && !mp->loadWeights) // Save weights for EfE synapses - Pretraining and Testing
		{
			for (wl=1; wl<mp->nLayers; wl++)
			{
				slen = snprintf(filename, FNAMEBUFF, "%sL%dweightsEfE.dat", prefix, wl);
				assert(slen < FNAMEBUFF);
				weightOutput = myfopen(filename, "w");
				for (n=0; n<mp->vExcit[wl]; n++)
				{
					fprintf(weightOutput, "%d\t", n_E[wl][n].nFAff_E); // Print number of EfE synapses first
					for (syn=0; syn<n_E[wl][n].nFAff_E; syn++)
						fprintf(weightOutput, "%f ", n_E[wl][n].FAffs_E[syn]->delta_g);
					fprintf(weightOutput, "\n");
				}
				fclose(weightOutput);
			}
			
			//(mp->trainElE) // *** Print out anyway?
			//{
            for (l=0; l<mp->nLayers; l++)
            {
                if (mp->pCnxElE[l] > EPS)
                {
					slen = snprintf(filename, FNAMEBUFF, "%sL%dweightsElE.dat", prefix, l);
					assert(slen < FNAMEBUFF);
					weightOutput = myfopen(filename, "w");
					for (n=0; n<mp->vExcit[l]; n++)
					{
						fprintf(weightOutput, "%d\t", n_E[l][n].nLAff_E); // Print number of ElE synapses first
						for (syn=0; syn<n_E[l][n].nLAff_E; syn++)
							fprintf(weightOutput, "%f ", n_E[l][n].LAffs_E[syn]->delta_g);
						fprintf(weightOutput, "\n");
					}
					fclose(weightOutput);
                }
            }
			//}
		}
		
		if (mp->nRecords) // Should be training //sPhase==Testing && 
		{
			char pStr[BUFSIZ];
			if (regime == Training)
				slen = snprintf(pStr, FNAMEBUFF, "RE%d", loop);
			else // <pre>Testing
				slen = snprintf(pStr, FNAMEBUFF, "R%s", prefix);

			assert(slen && slen < FNAMEBUFF); // Check non-negative
            
            // Could collapse these loops (and if statement) into one loop over array of N* ?
			for (l=0; l<mp->nLayers; l++)
			{
				for (n=0; n<mp->vExcit[l]; n++)
				{
					if (n_E[l][n].rec_flag)
					{
						slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dV.dat", pStr, l, n); // Change
						assert(slen && slen < FNAMEBUFF);
						recFile = myfopen(filename, "w"); //r_cellV_ptr = fopen(filename, "wb");
						print_frow(recFile, n_E[l][n].rec->cellV, n_E[l][n].rec->bin); // bin already points to the next free slot
						// print_farray(rCellVout, n_E[l][n].rec->cellV, mp->loops, mp->RecordMS); //fwrite(n_E[l][n].rec->cellV, sizeof(float), mp->loops*mp->RecordMS, r_cellV_ptr); // See nifty trick #1
						fclose(recFile);
						memset(n_E[l][n].rec->cellV, 0, (mp->RecordMS+1)*sizeof(*(n_E[l][n].rec->cellV)));
						
						if (mp->adaptation)
						{
							slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dcCa.dat", pStr, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w");
							print_frow(recFile, n_E[l][n].rec->cellcCa, n_E[l][n].rec->bin);
							fclose(recFile);
							memset(n_E[l][n].rec->cellcCa, 0, (mp->RecordMS+1)*sizeof(*(n_E[l][n].rec->cellcCa)));
						}
						
                        if (regime == Training) // *** && mp->train for preTraining?
                        {
                            slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dD.dat", pStr, l, n);
                            assert(slen && slen < FNAMEBUFF);
                            recFile = myfopen(filename, "w"); //r_D_ptr = fopen(filename, "wb");
                            print_frow(recFile, n_E[l][n].rec->cellD, n_E[l][n].rec->bin); // (rDout, n_E[l][n].rec->cellD, mp->loops, mp->RecordMS);
                            fclose(recFile);
                            memset(n_E[l][n].rec->cellD, 0, (mp->RecordMS+1)*sizeof(*(n_E[l][n].rec->cellD)));
                        }
                        
                        if (n_E[l][n].nLAff_I)
                        {
                            slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dLAffsigGI.dat", pStr, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w");
                            print_frow(recFile, n_E[l][n].rec->LsigGI, n_E[l][n].rec->bin);
							fclose(recFile);
							memset(n_E[l][n].rec->LsigGI, 0, (mp->RecordMS+1) * sizeof(*(n_E[l][n].rec->LsigGI)));
                        }
						
						if (n_E[l][n].nFAff_E) //(l>0) // The presynaptic cell's values are attached to each record
						{
                            slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dFAffg.dat", pStr, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w");
							print_farray(recFile, n_E[l][n].rec->FSynG, n_E[l][n].nFAff_E, n_E[l][n].rec->bin);
							fclose(recFile);
							memset(n_E[l][n].rec->FSynG[0], 0, n_E[l][n].nFAff_E*(mp->RecordMS+1)*sizeof(**(n_E[l][n].rec->FSynG)));
                            
                            if (regime == Training) // *** && mp->train for preTraining?
                            {
                                slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dFAffC.dat", pStr, l, n);
                                assert(slen && slen < FNAMEBUFF);
                                recFile = myfopen(filename, "w"); //r_C_ptr = fopen(filename, "wb");
                                print_farray(recFile, n_E[l][n].rec->FSynC, n_E[l][n].nFAff_E, n_E[l][n].rec->bin); //print_farray(recFile, n_E[l][n].rec->SynC[loop], n_E[l][n].nFAff_E, mp->RecordMS);
                                fclose(recFile);
                                memset(n_E[l][n].rec->FSynC[0], 0, n_E[l][n].nFAff_E*(mp->RecordMS+1)*sizeof(**(n_E[l][n].rec->FSynC)));
                                
                                slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dFAffdg.dat", pStr, l, n);
                                assert(slen && slen < FNAMEBUFF);
                                recFile = myfopen(filename, "w");
                                print_farray(recFile, n_E[l][n].rec->FSynDG, n_E[l][n].nFAff_E, n_E[l][n].rec->bin);
                                fclose(recFile);
                                memset(n_E[l][n].rec->FSynDG[0], 0, n_E[l][n].nFAff_E*(mp->RecordMS+1)*sizeof(**(n_E[l][n].rec->FSynDG)));
                            }
						}
						
						if (n_E[l][n].nLAff_E) //(mp->pCnxElE[l] > EPS)
						{
                            slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dLAffg.dat", pStr, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w");
							print_farray(recFile, n_E[l][n].rec->LSynG, n_E[l][n].nLAff_E, n_E[l][n].rec->bin);
							fclose(recFile);
							memset(n_E[l][n].rec->LSynG[0], 0, n_E[l][n].nLAff_E*(mp->RecordMS+1)*sizeof(**(n_E[l][n].rec->LSynG)));
                            
                            if (mp->trainElE && regime == Training) //mp->train
                            {
                                slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dLAffC.dat", pStr, l, n);
                                assert(slen && slen < FNAMEBUFF);
                                recFile = myfopen(filename, "w"); 
                                print_farray(recFile, n_E[l][n].rec->LSynC, n_E[l][n].nLAff_E, n_E[l][n].rec->bin); 
                                fclose(recFile);
                                memset(n_E[l][n].rec->LSynC[0], 0, n_E[l][n].nLAff_E*(mp->RecordMS+1)*sizeof(**(n_E[l][n].rec->LSynC)));
                                
                                slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dLAffdg.dat", pStr, l, n);
                                assert(slen && slen < FNAMEBUFF);
                                recFile = myfopen(filename, "w");
                                print_farray(recFile, n_E[l][n].rec->LSynDG, n_E[l][n].nLAff_E, n_E[l][n].rec->bin);
                                fclose(recFile);
                                memset(n_E[l][n].rec->LSynDG[0], 0, n_E[l][n].nLAff_E*(mp->RecordMS+1)*sizeof(**(n_E[l][n].rec->LSynDG)));
                            }
						}
						n_E[l][n].rec->bin = 0; // Reset record counter ready for next phase/epoch
					}
				}
			}
		} // End of state variable records
		printf("\tResults saved!\n");
	} // End of loops
	//myfree(input);
	return;
}

void updateNetwork(tstep t_start, tstep t_end, float input[], int regime) //int loop, 
{
	int l, n, syn; //, l; //, wl;
	int bin = 0;
	tstep t=0; tstep lstart=0;
	float decayRate, decay_E, decay_I;//, gLeak, Vrest, Thresh, Vhyper;
	
	
#pragma omp parallel default(shared) private(t,l,n,syn,bin,decayRate,decay_E,decay_I,lstart)//,gLeak,Vrest,Thresh,Vhyper)
	{
		
	if (t_start==0 && mp->nRecords) // Save intial neuron states (then every ms - see below) // Records printed to file after each loop
	{
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait 
			for (n=0; n<mp->vExcit[l]; n++) //for (r=0; r<NRECORDS_PL; r++)
				if (n_E[l][n].rec_flag) // Can use tm1 variables since they are initialised equal to tm0 variables
				{
					bin = n_E[l][n].rec->bin++;
					n_E[l][n].rec->cellV[bin] = n_E[l][n].V_tm1;
					if (mp->adaptation)
						n_E[l][n].rec->cellcCa[bin] = n_E[l][n].cCa_tm1;
                    for (syn=0; syn<n_E[l][n].nLAff_I; syn++)
                        n_E[l][n].rec->LsigGI[bin] += n_E[l][n].LAffs_I[syn]->g_tm1; // Sum
                    for (syn=0; syn<n_E[l][n].nFAff_E; syn++) // nFAff_E == 0 for l=0
                        n_E[l][n].rec->FSynG[syn][bin] = n_E[l][n].FAffs_E[syn]->g_tm1;
                    for (syn=0; syn<n_E[l][n].nLAff_E; syn++) // Implicit: if (mp->pCnxElE[l] > EPS) //mp->trainElE && mp->train && 
                        n_E[l][n].rec->LSynG[syn][bin] = n_E[l][n].LAffs_E[syn]->g_tm1;
                    if (regime == Training)
                    {
                        n_E[l][n].rec->cellD[bin] = n_E[l][n].D_tm1;
                        for (syn=0; syn<n_E[l][n].nFAff_E; syn++)
                        {
                            n_E[l][n].rec->FSynDG[syn][bin] = n_E[l][n].FAffs_E[syn]->delta_g_tm1;
                            n_E[l][n].rec->FSynC[syn][bin] = n_E[l][n].FAffs_E[syn]->C_tm1;
                        }
                        if (mp->trainElE) //&& mp->train
                            for (syn=0; syn<n_E[l][n].nLAff_E; syn++)
                            {
                                n_E[l][n].rec->LSynDG[syn][bin] = n_E[l][n].LAffs_E[syn]->delta_g_tm1;
                                n_E[l][n].rec->LSynC[syn][bin] = n_E[l][n].LAffs_E[syn]->C_tm1;
                            }
                    }
				}
		}
//#pragma omp barrier // Only need a barrier before solution variables are copied
	}


	for (t=t_start; t<t_end; t++)
	{
	
#if DEBUG>3
#pragma omp single nowait
		{
			fprintf(stderr, "\rTimestep: %d",t);
			fflush(stderr);
		}
#endif

		/* Update Excitatory cell potentials */
		decayRate = mp->DT/mp->capE;
		/*gLeak = mp->gLeakE;
		Vrest = mp->VrestE;
		Thresh = mp->ThreshE;
		Vhyper = mp->VhyperE;*/
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait schedule(static) // Needed for thread-safe noise
			for (n=0; n<mp->vExcit[l]; n++)
				update_V(&n_E[l][n], t, decayRate, mp->gLeakE, mp->VrestE, mp->ThreshE, mp->VhyperE, (l==0)?input[n]:0.0);
				//update_V(&n_E[l][n], t, decayRate, gLeak, Vrest, Thresh, Vhyper, (l==0)?input[n]:0.0);
		}
		
		/* Update Inhibitory cell potentials */
		decayRate = mp->DT/mp->capI;
		/*gLeak = mp->gLeakI;
		Vrest = mp->VrestI;
		Thresh = mp->ThreshI;
		Vhyper = mp->VhyperI;*/
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait schedule(static)// It is ok to have a nowait here with a minimum conduction delay of 1 tstep
			for (n=0; n<mp->vInhib[l]; n++) // Larger input layer separation as for Excit cells?
				update_V(&n_I[l][n], t, decayRate, mp->gLeakI, mp->VrestI, mp->ThreshI, mp->VhyperI, 0.0);
				//update_V(&n_I[l][n], t, decayRate, gLeak, Vrest, Thresh, Vhyper, 0.0);
		}

		if (mp->adaptation)
		{
			decayRate = mp->DT/mp->tauCa;
			for (l=0; l<mp->nLayers; l++)
			{
#pragma omp for nowait //schedule(runtime)
				for (n=0; n<mp->vExcit[l]; n++)
					update_cCa(&n_E[l][n], t, decayRate);
			}
		}
		
		/* Update presynaptic conductances */
		decay_E = (mp->DT/mp->tauEE);
		decay_I = (mp->DT/mp->tauIE);
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait //schedule(guided)//schedule(runtime) //- experiment!!
			for (n=0; n<mp->vExcit[l]; n++)
				update_g(&n_E[l][n], t, decay_E, decay_I); // Uses same decay_E for EfE and ElE synapses
		}
		
		decay_E = (mp->DT/mp->tauEI);
		decay_I = (mp->DT/mp->tauII);
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait //schedule(guided)//schedule(runtime)
			for (n=0; n<mp->vInhib[l]; n++)
				update_g(&n_I[l][n], t, decay_E, decay_I);
		}
#pragma omp barrier
		
		if (regime==Training) // Learning
		{
			lstart = (mp->trainElE) ? 0 : 1;
			/* Update synaptic weights */ // Move after C & D and use instantaeous values? N.B. Redo nowait clauses
			for (l=lstart; l<mp->nLayers; l++) // Skip 1st layer of afferent processing only for FF plasticity
			{
#pragma omp for nowait //schedule(runtime)
				for (n=0; n<mp->vExcit[l]; n++) // Make parallel and n private
					update_weights(&n_E[l][n], t);
			}
			
			// -->|| Update C for current neuron's outgoing synapses (skip last layer)
			decayRate = mp->DT/mp->tauC;
			for (l=lstart; l<mp->nLayers; l++) // Only nLayers-1 of synapses
			{
#pragma omp for nowait //schedule(runtime) //ok
				for (n=0; n<mp->vExcit[l]; n++)
					update_C(&n_E[l][n], t, decayRate);
			}
			
			// ||--> Update D for current neuron's incoming synapses (skip first layer)
			decayRate = mp->DT/mp->tauD; 
			for (l=lstart; l<mp->nLayers; l++) // Only nLayers-1 of weight layers
			{
#pragma omp for nowait //schedule(runtime) //ok
				for (n=0; n<mp->vExcit[l]; n++)
					update_D(&n_E[l][n], t, decayRate);
			}
#pragma omp barrier
		}
//		#pragma omp barrier
		
		/* Copy solution variables to _tm1 counterparts and reset spike flags / axons */
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait //schedule(runtime) //ok
			for (n=0; n<mp->vExcit[l]; n++)
			{
				//n_E[l][n].fired = 0;
				n_E[l][n].V_tm1 = n_E[l][n].V;
				if (mp->adaptation)
					n_E[l][n].cCa_tm1 = n_E[l][n].cCa;
				for (syn=0; syn<n_E[l][n].nFEff_E; syn++)
				{
					if (t == next_spike(&n_E[l][n].FEffs_E[syn]))
						dequeue(&n_E[l][n].FEffs_E[syn]);
					n_E[l][n].FEffs_E[syn].g_tm1 = n_E[l][n].FEffs_E[syn].g;
				}
				for (syn=0; syn<n_E[l][n].nLEff_E; syn++)
				{
					if (t == next_spike(&n_E[l][n].LEffs_E[syn]))
						dequeue(&n_E[l][n].LEffs_E[syn]);
					n_E[l][n].LEffs_E[syn].g_tm1 = n_E[l][n].LEffs_E[syn].g;
				}
				for (syn=0; syn<n_E[l][n].nLEff_I; syn++)
				{
					if (t == next_spike(&n_E[l][n].LEffs_I[syn]))
						dequeue(&n_E[l][n].LEffs_I[syn]);
					n_E[l][n].LEffs_I[syn].g_tm1 = n_E[l][n].LEffs_I[syn].g;
				}
                if (regime == Training) // Only copy C, D & Dg during Training
                {
                    n_E[l][n].D_tm1 = n_E[l][n].D;
                    for (syn=0; syn<n_E[l][n].nFEff_E; syn++)
                    {
                        n_E[l][n].FEffs_E[syn].delta_g_tm1 = n_E[l][n].FEffs_E[syn].delta_g;
                        n_E[l][n].FEffs_E[syn].C_tm1 = n_E[l][n].FEffs_E[syn].C;
                    }
                    if (mp->trainElE)
                        for (syn=0; syn<n_E[l][n].nLEff_E; syn++)
                        {
                            n_E[l][n].LEffs_E[syn].delta_g_tm1 = n_E[l][n].LEffs_E[syn].delta_g;
                            n_E[l][n].LEffs_E[syn].C_tm1 = n_E[l][n].LEffs_E[syn].C;
                        }
                }
			}
		}
		
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait //schedule(runtime)
			for (n=0; n<mp->vInhib[l]; n++)
			{
				n_I[l][n].V_tm1 = n_I[l][n].V;
				for (syn=0; syn<n_I[l][n].nLEff_E; syn++)
				{
					if (t == next_spike(&n_I[l][n].LEffs_E[syn]))
						dequeue(&n_I[l][n].LEffs_E[syn]);
					n_I[l][n].LEffs_E[syn].g_tm1 = n_I[l][n].LEffs_E[syn].g;
				}
				for (syn=0; syn<n_I[l][n].nLEff_I; syn++)
				{
					if (t == next_spike(&n_I[l][n].LEffs_I[syn]))
						dequeue(&n_I[l][n].LEffs_I[syn]);
					n_I[l][n].LEffs_I[syn].g_tm1 = n_I[l][n].LEffs_I[syn].g;
				}
			}
		}
//#pragma omp barrier
		
		if (mp->nRecords && ((t+1) % mp->TSperMS == 0)) // && regime // Save neuron states (dump to file after each loop)
		{
			for (l=0; l<mp->nLayers; l++)
			{
#pragma omp for nowait //schedule(runtime)//private(bin, syn, n) //l, // Already private!
				for (n=0; n<mp->vExcit[l]; n++) //for (r=0; r<NRECORDS_PL; r++)
					if (n_E[l][n].rec_flag)
					{
						bin = n_E[l][n].rec->bin++;
						n_E[l][n].rec->cellV[bin] = n_E[l][n].V;
						if (mp->adaptation)
							n_E[l][n].rec->cellcCa[bin] = n_E[l][n].cCa;
                        for (syn=0; syn<n_E[l][n].nLAff_I; syn++)
                            n_E[l][n].rec->LsigGI[bin] += n_E[l][n].LAffs_I[syn]->g; // Sum
                        for (syn=0; syn<n_E[l][n].nFAff_E; syn++)  // nFAff_E == 0 for l=0
                            n_E[l][n].rec->FSynG[syn][bin] = n_E[l][n].FAffs_E[syn]->g;
                        for (syn=0; syn<n_E[l][n].nLAff_E; syn++) // Implicit: if (mp->pCnxElE[l] > EPS) //mp->trainElE && mp->train && 
                            n_E[l][n].rec->LSynG[syn][bin] = n_E[l][n].LAffs_E[syn]->g;
                        if (regime == Training)
                        {
                            n_E[l][n].rec->cellD[bin] = n_E[l][n].D;
                            for (syn=0; syn<n_E[l][n].nFAff_E; syn++) // nFAff_E == 0 for l=0
                            {
                                n_E[l][n].rec->FSynDG[syn][bin] = n_E[l][n].FAffs_E[syn]->delta_g;
                                n_E[l][n].rec->FSynC[syn][bin] = n_E[l][n].FAffs_E[syn]->C;
                            }
                            if (mp->trainElE) // && mp->train && mp->pCnxElE[l] > EPS
                                for (syn=0; syn<n_E[l][n].nLAff_E; syn++)
                                {
                                    n_E[l][n].rec->LSynDG[syn][bin] = n_E[l][n].LAffs_E[syn]->delta_g;
                                    n_E[l][n].rec->LSynC[syn][bin] = n_E[l][n].LAffs_E[syn]->C;
                                }
                        }
					}
			}
			//#pragma omp barrier
		}
		
#pragma omp barrier
	} // End of timesteps loop
	
	} // End of parallel section
	return; // void;
}


void update_V(NEURON * n, tstep t, float decay_rate, float gLeak, float Vrest, float Thresh, float Vhyper, float inj)
{	/* This function updates a neuron's cell potential and applies to all layers */

	if (t > n->nextUpdate) //if ((((t - n->lastSpike) * DT) > mp->refract) || n->spkbin==0) // Calculated at spike time
	{
		int syn;
		float tot_g_E = 0.0;
		float tot_g_I = 0.0;
		
		/* Feed-forward connections (n->type==EXCIT && l>0) */
		for (syn=0; syn<n->nFAff_E; syn++)
			tot_g_E += n->FAffs_E[syn]->g_tm1;
		
		/* Lateral connections */
		for (syn=0; syn<n->nLAff_E; syn++)
			tot_g_E += n->LAffs_E[syn]->g_tm1;
		
		for (syn=0; syn<n->nLAff_I; syn++)
			tot_g_I += n->LAffs_I[syn]->g_tm1;
		
		/*float adapt = 0.0;
		if (mp->adaptation && n->type==EXCIT)
			adapt = mp->gAHP * (mp->VK - n->V_tm1);*/
		
		n->V += decay_rate * ( (gLeak * (Vrest - n->V_tm1)) \
							  + (tot_g_E * (mp->VrevE - n->V_tm1)) \
							  + (tot_g_I * (mp->VrevI - n->V_tm1)) \
							  + ((mp->adaptation && n->type==EXCIT) ? (mp->gAHP * n->cCa_tm1 * (mp->VK - n->V_tm1)) : 0.0) \
							  + inj );
		
		if (mp->noise)
		{
			float sigma = 0.0;
			int th = 0;
#ifdef _OPENMP
			th = omp_get_thread_num();
#endif
			sigma = (n->type == EXCIT) ? mp->SigmaE : mp->SigmaI;
			n->V += (gsl_ran_gaussian(states[th], sqrt(gLeak*decay_rate)) * sigma); // decay_rate = DT/mp->cap{E,I}
		}
		
		n->V = ((n->V < mp->VrevI) ? mp->VrevI : n->V); // Neurons can not be more -ve than Inhib reversal potential
		
		if (n->V >= Thresh) // n->V_tm1 ??? - gives one timestep to depolarise (peak)
		{
#if DEBUG>3
			fprintf(stderr," L%dN%d%c: t=%d",n->l,n->n,(n->type==EXCIT)?'E':'I',t);
#endif
			n->V = Vhyper;
			n->lastSpike = n->spikeTimes[n->spkbin++] = t; //++(n->spkbin)
			n->nextUpdate = t + ceil(mp->refract/mp->DT);
			// Enqueue spike for all post-synaptic neurons
			for (syn=0; syn<n->nFEff_E; syn++)
				enqueue(&n->FEffs_E[syn], t);
						
			for (syn=0; syn<n->nLEff_E; syn++)
				enqueue(&n->LEffs_E[syn], t);
						
			for (syn=0; syn<n->nLEff_I; syn++)
				enqueue(&n->LEffs_I[syn], t);
			
			/*if (mp->adaptation && n->type==EXCIT)
				n->cCa = n->cCa_tm1 + mp->alphaCa;*/
		} // End of spike
	} // End of REFRACT check 
	return;
}

inline void update_cCa(NEURON * n, tstep t, float decayRate)
{
	//float impulse = ((t == n->lastSpike) ? mp->alphaCa : 0.0); // private
	n->cCa += (((t == n->lastSpike) ? mp->alphaCa : 0.0) - (n->cCa_tm1 * decayRate)); //mp->DT/mp->tauCa)); // Unbounded!
	return;
}

inline void update_g(NEURON * n, tstep t, float decay_E, float decay_I)
{
	/* This function updates a neuron's afferent (pre-synaptic) conductances and applies to l>0 */
	float impulse;
	int syn;
	float scale = 0.0;
	
	/* Update EfE synapses */
	scale = mp->DgEfE * mp->gMax;
	for (syn=0; syn<n->nFAff_E; syn++) // Make private
	{
		impulse = (next_spike(n->FAffs_E[syn]) == t) ? n->FAffs_E[syn]->delta_g_tm1 * scale : 0.0;
		n->FAffs_E[syn]->g += (impulse - (decay_E * n->FAffs_E[syn]->g_tm1));
	}
	
	/* Update El_ synapses */
	//scale = (mp->trainElE && n->type==EXCIT) ? mp->DgElE * mp->gMax : mp->gMax;
	scale = (n->type==EXCIT) ? mp->DgElE * mp->gMax : mp->gMax; // THINK!!!
	for (syn=0; syn<n->nLAff_E; syn++)
	{
		impulse = (next_spike(n->LAffs_E[syn]) == t) ? n->LAffs_E[syn]->delta_g_tm1 * scale : 0.0;
		n->LAffs_E[syn]->g += (impulse - (decay_E * n->LAffs_E[syn]->g_tm1));
	}
	
	/* Update I synapses */
	for (syn=0; syn<n->nLAff_I; syn++) // Already private
	{
		impulse = (next_spike(n->LAffs_I[syn]) == t) ? n->LAffs_I[syn]->delta_g_tm1 * mp->gMax : 0.0;
		n->LAffs_I[syn]->g += (impulse - (decay_I * n->LAffs_I[syn]->g_tm1));
	}
	
	// Check to make sure conductances have not become -ve through a coarse decay rate???
	return;
}

inline void update_weights(NEURON * n, tstep t)
{
	/*** Learning at Excitatory to Excitatory synapses ***/
	/* Adapted from Perrinet, Delorme, Samuelides & Thorpe 2001 */
	float LTD; //contrib_D
	float LTP; //contrib_C;
	int syn;
	float Dg_tm1 = 0.0;
	//float Dg = 0.0;
	if (mp->trainEfE)
	{
		for (syn=0; syn<n->nFAff_E; syn++)
		{
			Dg_tm1 = n->FAffs_E[syn]->delta_g_tm1;
#ifndef __llvm__
			assert(0 <= Dg_tm1 && Dg_tm1 <= 1);
#endif
			LTD = (t == next_spike(n->FAffs_E[syn])) ? n->FAffs_E[syn]->delta_g_tm1 * n->D_tm1 : 0.0; //tm0?
			LTP = (t == n->lastSpike) ? (1 - n->FAffs_E[syn]->delta_g_tm1) * n->FAffs_E[syn]->C_tm1 : 0.0; //tm0?
			n->FAffs_E[syn]->delta_g += (LTP - LTD) * mp->learnR; //*DT/TAU_DG;
			//Dg = n->FAffs_E[syn]->delta_g;
		}
	}
	
	if (mp->trainElE) // 
	{
		for (syn=0; syn<n->nLAff_E; syn++)
		{
			Dg_tm1 = n->LAffs_E[syn]->delta_g_tm1;
#ifndef __llvm__
			assert(0 <= Dg_tm1 && Dg_tm1 <= 1);
#endif
			LTD = (t == next_spike(n->LAffs_E[syn])) ? n->LAffs_E[syn]->delta_g_tm1 * n->D_tm1 : 0.0; //tm0?
			LTP = (t == n->lastSpike) ? (1 - n->LAffs_E[syn]->delta_g_tm1) * n->LAffs_E[syn]->C_tm1 : 0.0; //tm0?
			n->LAffs_E[syn]->delta_g += (LTP - LTD) * mp->learnR; //*DT/TAU_DG;
			//assert(0 <= n->LAffs_E[syn]->delta_g && n->LAffs_E[syn]->delta_g <= 1);
		}
	}
	return;
}
	
inline void update_C(NEURON * n, tstep t, float decayRate)
{
	// -->|| Update C for current neuron's outgoing synapses (skip last layer)
	int syn;
	float impulse;
	float C_tm1 = 0.0;
	//float decayRate = mp->DT/mp->tauC;

	if (mp->trainEfE) // Loop over afferent synapses (skip first layer)
	{
		for (syn=0; syn<n->nFAff_E; syn++)
		{
			C_tm1 = n->FAffs_E[syn]->C_tm1;
			impulse = (t == next_spike(n->FAffs_E[syn])) ? mp->alphaC * (1 - C_tm1) : 0.0;
			n->FAffs_E[syn]->C += (impulse - (C_tm1 * decayRate));
		}
	}
	
	if (mp->trainElE) // Update Excitatory lateral Afferent synapses
	{
		for (syn=0; syn<n->nLAff_E; syn++)
		{
			C_tm1 = n->LAffs_E[syn]->C_tm1;
			//n->LAffs_E[syn]->C += ((t==next_spike(n->LAffs_E[syn])) ? (mp->alphaC*(1-C_tm1)) : 0.0) - (C_tm1*decayRate);
			impulse = (t == next_spike(n->LAffs_E[syn])) ? (mp->alphaC * (1 - C_tm1)) : 0.0;
			n->LAffs_E[syn]->C += (impulse - (C_tm1 * decayRate));
			//assert(n->LAffs_E[syn]->C <= 1.0 && n->LAffs_E[syn]->C >= 0.0);
		}
	}
	return;
}
	
inline void update_D(NEURON * n, tstep t, float decayRate)
{
	// ||--> Update D for current neuron's incoming synapses (skip first layer)
	/*float impulse = ((t == n->lastSpike) ? mp->alphaD : 0.0); // private
	n->D += (impulse * (1 - n->D_tm1) - (n->D_tm1 * mp->DT/mp->tauD));*/
	//n->D += ((((t == n->lastSpike) ? mp->alphaD : 0.0) * (1 - n->D_tm1)) - (n->D_tm1 * decayRate)); //mp->DT/mp->tauD));
	n->D += ((t == n->lastSpike) ? (mp->alphaD * (1 - n->D_tm1)) : 0.0) - (n->D_tm1 * decayRate);
	return;
}


int normalise(NEURON ** narray, PARAMS * mp)
{
	int l = 0;
	int n = 0;
	int s = 0;
	double sf = 0.0;
	double sum = 0.0;
	
	switch (mp->normalise) 
	{
		case None:
		{
			printf("No normalisation.\n");
			return 0;
		}
		case MaintainLength: // Standard normalisation: set sum of squares to 1
		{
			// EfE weights
			for (l=1; l<mp->nLayers; l++)
			{
				for (n=0; n<mp->vExcit[l]; n++) // Parallelise
				{
					sf = 0.0;
					for (s=0; s<narray[l][n].nFAff_E; s++)
						sf += (narray[l][n].FAffs_E[s]->delta_g_tm1 * n_E[l][n].FAffs_E[s]->delta_g_tm1); //tm1?
					if (sf) // Do not divide by 0!
					{
						sf = 1/sqrt(sf);
						for (s=0; s<narray[l][n].nFAff_E; s++)
						{
							narray[l][n].FAffs_E[s]->delta_g_tm1 *= sf;
							if (narray[l][n].FAffs_E[s]->delta_g_tm1 < 0.0)
								narray[l][n].FAffs_E[s]->delta_g_tm1 = 0.0;
							if (narray[l][n].FAffs_E[s]->delta_g_tm1 > 1.0)
								narray[l][n].FAffs_E[s]->delta_g_tm1 = 1.0;
						}
					}
				}
			}
			break;
		}
			
		case MaintainSum: // Maintain sum of weights
		{
			for (l=1; l<mp->nLayers; l++)
			{
				for (n=0; n<mp->vExcit[l]; n++)
				{
					sum = 0.0;
					for (s=0; s<narray[l][n].nFAff_E; s++)
					{
						sum += narray[l][n].FAffs_E[s]->delta_g_tm1; // tm1 since normalisation comes after solution vars are copied
					}
					if (sum) // Do not divide by 0!
					{
						sf = narray[l][n].nFAff_E / (2 * sum);
						for (s=0; s<narray[l][n].nFAff_E; s++)
						{
							narray[l][n].FAffs_E[s]->delta_g_tm1 *= sf;
							if (narray[l][n].FAffs_E[s]->delta_g_tm1 < 0.0)
								narray[l][n].FAffs_E[s]->delta_g_tm1 = 0.0;
							if (narray[l][n].FAffs_E[s]->delta_g_tm1 > 1.0)
								narray[l][n].FAffs_E[s]->delta_g_tm1 = 1.0;
						}
					}
				}
			}
			break;
		}
			
		default:
		{
			printf("Unknown normalisation mode.\n");
			return 1;
			break;
		}
	}
	return 0;
}


inline void init_queue(AXON *a)
{
	int bin;
	a->next = 0;
	a->last = a->nBins-1;
	a->count = 0;
	for (bin=0; bin<a->nBins; bin++)
		a->queue[bin] = -BIG;
	return; // Necessary?
}

inline void enqueue(AXON *a, tstep t)
{
#ifndef __llvm__
	assert(a->count < a->nBins);
#else
    #ifndef NDEBUG
    if (a->count >= a->nBins)
        exit_error("enqueue", "Queue full");
    #endif
#endif
	a->last = (a->last+1) % a->nBins;
	a->queue[ a->last ] = t + a->delay;
	a->count++;
	return; //a->count++;
}

inline int dequeue(AXON *a)
{
	int t;
#ifndef __llvm__
	assert(a->count > 0);
#else
    #ifndef NDEBUG
    if (a->count <= 0) // isempty(a)
        exit_error("dequeue", "Queue empty");
    #endif
#endif
	t = a->queue[ a->next ];
	a->next = (a->next+1) % a->nBins;
	a->count--;
	return(t);
}

inline int next_spike(AXON * a)
{
	return a->queue[a->next];
}

inline bool isempty(AXON *a)
{
	return (a->count == 0) ? true : false; // Originally <=
}

inline int print_queue(AXON *a)
{
	int i = a->next; 
	while (i != a->last) 
	{
		printf("%d ",a->queue[i]); // Was %c
		i = (i+1) % a->nBins;
	}
	printf("%d ",a->queue[i]);
	printf("\n");
	return a->count;
}
