/*
 *  spike.c
 *  Spike
 *
 *  Created by Ben Evans on 6/19/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "spike.h"

int spike(PARAMS * mp)
{
	/*** Declare variables ***/
	int p;//, l;
	int error = 0;
	REGIMETYPE regime;
	STIMULI * stim;

	//RECORD *RECSP = RECS;
	//RECORD *ptr = &RECS[0][0];
	//RECORD *r_ptr = RECS;
	
	/*** Declare file pointers ***/
	FILE * stimuli_FP;
	
	
	/*************** Calculate RAM requirements ****************/
	
#if DEBUG > 1
	fprintf(stdout, "Variable type:\tNEURON\tAXON  \t*     \tfloat \ttstep \tint\n");
	fprintf(stdout, "Size (bytes): \t%-6lu\t%-6lu\t%-6lu\t%-6lu\t%-6lu\t%-6lu\t\n",\
			sizeof(NEURON),sizeof(AXON),sizeof(int*),sizeof(float),sizeof(tstep),sizeof(int));
#endif
	//size_t size = sizeof();
	int EsynE = 0;
	int EsynI = 0;
	float memE = 0.0;
	float memI = 0.0;
	float memMisc = 0.0;
	float memTrain = 0.0;
	float memTest = 0.0;
	float Tmem = 0.0;
	float avq = 0;
	int EsynElE = 0;
	
	int MB = 1024*1024;
	int l=0;
#if DEBUG > 1
	fprintf(stdout, "\nLayer\tExcit \tInhib \tSize (MB)\n");
#endif
	for (l=0; l<mp->nLayers; l++)
	{
		Tmem += (mp->vExcit[l]+mp->vInhib[l])*(sizeof(NEURON)+(mp->spkBuffer*(float)sizeof(tstep)))/MB;
#if DEBUG > 1
		fprintf(stdout,"%-6d\t%-6d\t%-6d\t%-6.2f\n",l,mp->vExcit[l],mp->vInhib[l],\
				(mp->vExcit[l]+mp->vInhib[l])*(sizeof(NEURON)+(mp->spkBuffer*(float)sizeof(tstep)))/MB);
#endif
	}
	
	switch (mp->axonDelay) 
	{
		case MinD:
			avq = 1;
			break;
		case ConstD:
			avq = round(mp->d_const/mp->DT);
			break;
		case UniformD:
			avq = (mp->d_min + ((mp->d_max - mp->d_min) / 2))/mp->DT;
			break;
		case GaussD:
			avq = mp->d_mean / mp->DT;
			break;
		case SOMD:
			avq = mp->maxDelay / (2.0 * mp->DT); // Reasonable for 1D layer
			break;
		default:
			exit_error("create_axons", "Unknown axonal delay model!");
			break;
	}
	avq = (!avq) ? 1 : avq;
	
#if DEBUG > 1
	//fprintf(stderr,"\nPresynaptic connection probabilites for Excitatory postsynaptic cells\n");
	fprintf(stdout,"\n\tp(EfE)\tp(ElE)\tp(IE) \tE[syn] \tp(EI) \tp(II) \tE[syn]\tSize (MB)\n");
#endif
	// Assumes minimum delay model i.e. 1 spike bin per axon for all except ElE
	for (l=0; l<mp->nLayers; l++)
	{
		EsynE = (((l>0)?(mp->vExcit[l-1]*mp->pCnxEfE[l]):0)+(mp->vInhib[l]*mp->pCnxIE[l]))*mp->vExcit[l];
		EsynElE = (pow(mp->vExcit[l],2)*mp->pCnxElE[l]);
		memE = ((EsynE + EsynElE)*(sizeof(AXON) + (sizeof(NEURON*)*3)) + (sizeof(tstep)*(EsynE + avq*EsynElE)))/MB;
		EsynI = ((mp->vExcit[l]*mp->pCnxEI[l])+(mp->vInhib[l]*mp->pCnxII[l]))*mp->vInhib[l];
		memI = (sizeof(AXON) + (sizeof(NEURON*)*3) + sizeof(tstep))*(float)EsynI/MB;
#if DEBUG > 1
		fprintf(stdout,"L%d:\t%-6.3f\t%-6.3f\t%-6.3f\t%-6.2G\t%-6.3f\t%-6.3f\t%-6.2G\t%-6.2f\n", \
				l,mp->pCnxEfE[l],mp->pCnxElE[l],mp->pCnxIE[l],(float)EsynE+EsynElE, \
				mp->pCnxEI[l],mp->pCnxII[l],(float)EsynI,memE+memI);
#endif
		Tmem += (memE + memI);
	}

#if DEBUG > 1
	fprintf(stdout, "\nStimuli:\tStruct.\tTrain \tTest  \tRecords\n");
#endif
	if (mp->SOM) // (mp->initElE == SOM || mp->initEfE == SOM || mp->axonDelay == SOMD) //mp->SOM // Triangle of Distances
	{
		for (l=0; l<mp->nLayers; l++)
			memMisc += mp->vExcit[l]*(mp->vExcit[l]+1)/2; // Gauss' method 
		memMisc *= sizeof(float)/MB;
	}
	
	memE = 0.0;
	if (mp->nRecords)
	{
		memE = (mp->vRecords[0]*(sizeof(RECORD)+((mp->adaptation)?3:2)*(mp->TotalMS+1)*sizeof(float)))/MB;
		for (l=1; l<mp->nLayers; l++)
			memE += mp->vRecords[l]*3*(mp->TotalMS+1)*(mp->vExcit[l-1]*mp->pCnxEfE[l])*sizeof(float)/MB;
		for (l=0; l<mp->nLayers; l++)
			memE += mp->vRecords[l]*3*(mp->TotalMS+1)*(mp->vExcit[l]*mp->pCnxElE[l])*sizeof(float)/MB;
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
	
	
	printf("Ventral Visual Stream Spiking Neural Network Simulation starting...\n");
	
	/*************** Build & Initialize Network ****************/
	
	printf("\tNow building the network...");
	
	n_E = allocn(mp->nLayers, mp->vExcit, EXCIT); // Create 2D array of Excitatory neuron structures
	n_I = allocn(mp->nLayers, mp->vInhib, INHIB); // Create 2D array of Inhibatory neuron structures
	if (mp->SOM) //(mp->initElE == SOM || mp->initEfE == SOM || mp->axonDelay == SOMD)
		distE = getLowTriF(mp->nLayers, mp->vExcit, 0.0);
	calcConnectivity(mp->probConnect);	// Calculate network connectivity

	printf("\tBuilding complete!\n");
	
	printf("\tNow initialising the network...");
	
	regime = Learning;			// 0: Testing (No STDP); 1: Training (STDP);
	init_network(regime);		// Initialise parameters
	
	/*************** Load stimuli ****************/
	
	printf("\tNetwork initialisation complete!\n");
	
	stim = myalloc(sizeof(*stim));
	stim->trn_stimuli = stim->tst_stimuli = NULL;
	stim->stimShuffle = NULL; 
	stim->transShuffle = NULL;
	stim->trnImages = stim->tstImages = NULL;

	if (mp->useFilteredImages)
	{
		printf("\tNow loading the images...");
		stim->trnImages = get_7D_farray(mp->nStimuli, mp->nTransPS, \
										mp->nScales, mp->nOrients, mp->nPhases, mp->nRows, mp->nCols, 0.0);
		
		error = loadImages(stim, mp); // if (!error)...
		printf("\t{S%d,T%d}", mp->nStimuli, mp->nTransPS);
		if (mp->newTestSet)
			printf(" Test: {S%d,T%d}", mp->nTestStimuli, mp->nTestTransPS);
		else //if (!mp->newTestSet)
		{
			stim->tstImages = stim->trnImages;
			mp->nTestStimuli = mp->nStimuli;
			mp->nTestTransPS = mp->nTransPS;
		}
		
		if (!error)
			printf("\tImages now loaded!\n");
		else
			exit_error("spike", "Error loading Images");
	}
	else
	{
		printf("\tNow creating the stimuli...");
		gen_stimuli(mp->localRep, stim, mp);			// Generate Patterns
		
		stimuli_FP = fopen("trn_stimuli.dat", "w");
		for (p=0; p<mp->nStimuli; p++)
		{
			fprintf(stimuli_FP, "*** Stimulus %d ***\n", p+1);
			print_farray(stimuli_FP, (stim->trn_stimuli)[p], mp->nTransPS, mp->sInputs);
			fprintf(stimuli_FP, "\n");
		}
		fclose(stimuli_FP);
		
		if (mp->newTestSet)
		{
			stimuli_FP = fopen("tst_stimuli.dat", "w");
			for (p=0; p<mp->nTestStimuli; p++)
			{
				fprintf(stimuli_FP, "*** Stimulus %d ***\n", p+1);
				print_farray(stimuli_FP, (stim->tst_stimuli)[p], mp->nTestTransPS, mp->sInputs);
				fprintf(stimuli_FP, "\n");
			}
			fclose(stimuli_FP);
		}
		
		printf("\tStimuli now saved!\n");
	}
	
	genShuffles(stim, mp);

	if (SIM.Xgrid) // Set up variables for estimating percentage completion
	{
		SIM.tally = 0;
		SIM.ptTS = (mp->pretrain) ? mp->nTestStimuli * mp->nTestTransPS * mp->transP_Test * ceil(1/mp->DT): 0;
		SIM.trainTS = (mp->train) ? mp->loops * mp->nStimuli * mp->nTransPS * mp->transP_Train * ceil(1/mp->DT): 0;
		SIM.testTS = mp->nTestStimuli * mp->nTestTransPS * mp->transP_Test * ceil(1/mp->DT);
		SIM.totTS = SIM.ptTS + SIM.trainTS + SIM.testTS;
		assert(SIM.totTS <= pow(2, 8*sizeof(tstep)-1)-1); // assumes tstep will be unsigned otherwise pow(2, 8*sizeof(tstep))-1
	}
	
	/*************** Simulation Phases ****************/
	
	if (mp->pretrain)
		simulatePhase(PreTraining, stim);
	
	if (mp->train)
		simulatePhase(Training, stim);
	
	simulatePhase(Testing, stim);

	
	/*************** Deallocate Memory ****************/
	
	unallocn(n_E, mp->nLayers, mp->vExcit);
	unallocn(n_I, mp->nLayers, mp->vInhib);
	if (mp->SOM) // (mp->initElE == SOM || mp->initEfE == SOM || mp->axonDelay == SOMD)
		freeTriF(distE, mp->nLayers);
	
	if (mp->randStimOrder)
		free_2D_iarray(stim->stimShuffle);//, mp->loops);
	if (mp->randTransOrder)
		free_3D_iarray(stim->transShuffle, mp->loops);//, mp->nStimuli);
	if (mp->useFilteredImages)
	{
		free_7D_farray(stim->trnImages, mp->nStimuli, mp->nTransPS, mp->nScales, mp->nOrients, mp->nPhases);
		if (mp->newTestSet)
			free_7D_farray(stim->tstImages, mp->nTestStimuli, mp->nTestTransPS, mp->nScales, mp->nOrients, mp->nPhases);
	}
	else
	{
		if (mp->newTestSet)
			free_3D_farray(stim->trn_stimuli, mp->nStimuli);//, mp->nTransPS);
		free_3D_farray(stim->tst_stimuli, mp->nTestStimuli);//, mp->nTransPS);
	}
	myfree(stim); // *** Free new image arrays too


	printf("\tMemory Deallocated!\n");

	return 0;
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
				myfree(narray[l][n].rec->cellD);//, mp->loops);
				if (l>0)  // Here the synaptic variables are stored with the post-synaptic neuron
				{
					free_2D_farray(narray[l][n].rec->FSynC);
					free_2D_farray(narray[l][n].rec->FSynG);
					free_2D_farray(narray[l][n].rec->FSynDG);
				}
				
				if (mp->trainElE && mp->train)
				{
					free_2D_farray(narray[l][n].rec->LSynC);
					free_2D_farray(narray[l][n].rec->LSynG);
					free_2D_farray(narray[l][n].rec->LSynDG);
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
					/*d_h = n_E[l][i].row - n_E[l][j].row;
					d_w = n_E[l][i].col - n_E[l][j].col;
					h_min = MIN(abs(d_h), abs(d_h + mp->layDim[l].nRows)); // Periodic boundary conditions
					w_min = MIN(abs(d_w), abs(d_w + mp->layDim[l].nCols)); // Periodic boundary conditions
					distE[l][i][j] = sqrt(pow(h_min,2) + pow(w_min,2)); // j is always > i*/
					
					distE[l][i][j] = calcDistance(&n_E[l][i], &n_E[l][j], mp->spatialScale);//mp->layDim);
#if DEBUG > 3
					printf("Dist: N%d[%d,%d] --> N%d[%d,%d] = %f\n", i, n_E[l][i].row, n_E[l][i].col, \
						   j, n_E[l][j].row, n_E[l][j].col, readLowTriF(distE, l, i, j));//distE[l][i][j]);
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
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vExcit[l]; n++)
			wireAfferents(&n_E[l][n], mp->pCnxEfE[l], mp->pCnxElE[l], mp->pCnxIE[l]);
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vInhib[l]; n++)
			wireAfferents(&n_I[l][n], 0.0, mp->pCnxEI[l], mp->pCnxII[l]);
	
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
#if DEBUG >	1
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
#endif
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vExcit[l]; n++)
			create_axons(&n_E[l][n]);
	
	for (l=0; l<mp->nLayers; l++)
		for (n=0; n<mp->vInhib[l]; n++)
			create_axons(&n_I[l][n]);
	
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
			tot = tot + n_E[l+1][n].nFAff_E;
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
	if (mp->printConnections)
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
					fprintf(connections_FP, "%d\t", n_E[l+1][n].lm1presyn_E[s]->n);
				fprintf(connections_FP, "\n");
			}
			fclose(connections_FP);
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
					fprintf(connections_FP, "%d\t", n_E[l][n].lm0presyn_E[s]->n);
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
						fprintf(connections_FP, "%f\t", readLowTriF(distE, l, i, j));//distE[l][i][j]); 
					fprintf(connections_FP, "\n");
				}
				fclose(connections_FP);
			}
			
			// Print delays
			for (l=0; l<mp->nLayers; l++)
			{
				slen = snprintf(filename, FNAMEBUFF, "L%daffDelaysElE.dat", l); 
				assert(slen < FNAMEBUFF);
				connections_FP = myfopen(filename, "w");
				for (n=0; n<mp->vExcit[l]; n++)
				{
					fprintf(connections_FP, "%d\t", n_E[l][n].nLAff_E); // Print number of ElE synapses first (in case any have 0)
					for (s=0; s<n_E[l][n].nLAff_E; s++)
						fprintf(connections_FP, "%d\t", n_E[l][n].LAffs_E[s]->delay);
					fprintf(connections_FP, "\n");
				}
				fclose(connections_FP);
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
					fprintf(connections_FP, "%d\t", n_I[l][n].lm0presyn_E[s]->n);
				fprintf(connections_FP, "\n");
			}
			fclose(connections_FP);
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
					fprintf(connections_FP, "%d\t", n_E[l][n].lm0presyn_I[s]->n);
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
					fprintf(connections_FP, "%d\t", n_I[l][n].lm0presyn_I[s]->n);
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

/*void calc_connectivity() //* This function randomly connects the neurons together *
{
	int n, l;
	int slen = 0;	
	char filename[FNAMEBUFF];
	FILE * connections_FP;
	
	//* Given a particular postsynaptic neuron <E|I>, n, and its synapse, s, the presynaptic cell forming the synapse <E|I> is given by presyncnx_<E|I><E|I>[l][n][s]. *
	//* NB If connections from the previous and the same layer are required, can label each neuron from 0 to (NEXCIT*NLAYERS)-1 (the neuron ID) with l=floor(NID/NEXCIT) and n=((NID+1)%NEXCIT)-1.	*
	
	//*	For affNeurons_EfE[0][n][s] : layer 0 cells --> layer 1 cells (NWLAYERS)
	//	For affNeurons_ElE[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS)
	//	For affNeurons_IE[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS)
	//	For affNeurons_EI[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS)
	//	For affNeurons_II[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS) *
	
	affNeurons_EfE = get_3D_iarray(mp->nWLayers, mp->nExcit, mp->nSynEfE, 0);
	affNeurons_ElE = get_3D_iarray(mp->nLayers, mp->nExcit, mp->nSynElE, 0);
	affNeurons_IE = get_3D_iarray(mp->nLayers, mp->nExcit, mp->nSynIE, 0);
	affNeurons_EI = get_3D_iarray(mp->nLayers, mp->nInhib, mp->nSynEI, 0);
	affNeurons_II = get_3D_iarray(mp->nLayers, mp->nInhib, mp->nSynII, 0);
	
#pragma omp for nowait private(n, l)
	for (n=0; n<mp->nExcit; n++)
	{
		for (l=0; l<mp->nLayers; l++)
		{
			n_E[l][n].nFEff_E = 0;
			n_E[l][n].nLEff_E = 0;
			n_E[l][n].nLEff_I = 0;
			n_E[l][n].n = n;
			n_E[l][n].l = l;
			n_E[l][n].nFAff_E = (l!=0) ? mp->nSynEfE : 0;
			n_E[l][n].nLAff_E = mp->nSynElE;
			n_E[l][n].nLAff_I = (l==0 && !mp->inputInhib) ? 0 : mp->nSynIE;
		}
	}
	
#pragma omp for private(n, l)
	for (n=0; n<mp->nInhib; n++)
	{
		for (l=0; l<mp->nLayers; l++)
		{
			n_I[l][n].nFEff_E = 0;
			n_I[l][n].nLEff_E = 0;
			n_I[l][n].nLEff_I = 0;
			n_I[l][n].n = n;
			n_I[l][n].l = l;
			n_I[l][n].nFAff_E = 0;
			n_I[l][n].nLAff_E = (l==0 && !mp->inputInhib) ? 0 : mp->nSynEI;
			n_I[l][n].nLAff_I = (l==0 && !mp->inputInhib) ? 0 : mp->nSynII;
		}
	}
	
	//***** Can rands1_new be used in parallel? *****
	
	for (n=0; n<mp->nExcit; n++)
		for (l=0; l<mp->nLayers; l++) // int * EfE_ptr = (l>0) ? affNeurons_EfE[l-1][n] : NULL;
			wire_afferents(&n_E[l][n], l, (l>0)?affNeurons_EfE[l-1][n]:NULL, (affNeurons_ElE)?affNeurons_ElE[l][n]:NULL, affNeurons_IE[l][n]);
		//((l>0)?affNeurons_EfE[l-1][n]:NULL)
			
	for (n=0; n<mp->nInhib; n++)
		for (l=0; l<mp->nLayers; l++)
			wire_afferents(&n_I[l][n], l, NULL, affNeurons_EI[l][n], affNeurons_II[l][n]);
	
#pragma omp for nowait private(n, l)
	for (n=0; n<mp->nExcit; n++)
		for (l=0; l<mp->nLayers; l++)
			alloc_efferents(&n_E[l][n]);
	
#pragma omp for private(n, l)
	for (n=0; n<mp->nInhib; n++)
		for (l=0; l<mp->nLayers; l++)
			alloc_efferents(&n_I[l][n]);
	
	for (n=0; n<mp->nExcit; n++)
		for (l=0; l<mp->nLayers; l++)
			wire_efferents(&n_E[l][n]);
	
	for (n=0; n<mp->nInhib; n++)
		for (l=0; l<mp->nLayers; l++)
			wire_efferents(&n_I[l][n]);
	
	//*** Axonal delays ***
#if DEBUG >	1
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
			
		default:
			printf("\nError: Unknown axonal delay model!\n");
			break;
	}
#endif
	
	for (n=0; n<mp->nExcit; n++)
		for (l=0; l<mp->nLayers; l++)
			create_axons(&n_E[l][n]);
	
	for (n=0; n<mp->nInhib; n++)
		for (l=0; l<mp->nLayers; l++)
			create_axons(&n_I[l][n]);
		
	//************** FILE OUTPUT **************
	
	for (l=0; l<mp->nWLayers; l++) // Loop up to nWLayers
	{
		slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsEfE.dat", l+1);
		assert(slen < FNAMEBUFF);
		connections_FP = fopen(filename, "w");
		print_iarray(connections_FP, affNeurons_EfE[l], mp->nExcit, mp->nSynEfE);
		fclose(connections_FP);
	}
	
	for (l=0; l<mp->nLayers; l++)
	{
		slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsElE.dat", l);
		assert(slen < FNAMEBUFF);
		connections_FP = fopen(filename, "w");
		print_iarray(connections_FP, (affNeurons_ElE)?affNeurons_ElE[l]:NULL, mp->nExcit, mp->nSynElE);
		fclose(connections_FP);
	}
	
	for (l=0; l<mp->nLayers; l++)
	{
		slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsEI.dat", l);
		assert(slen < FNAMEBUFF);
		connections_FP = fopen(filename, "w");
		print_iarray(connections_FP, affNeurons_EI[l], mp->nInhib, mp->nSynEI);
		fclose(connections_FP);
	}
	
	for (l=0; l<mp->nLayers; l++)
	{
		slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsIE.dat", l);
		assert(slen < FNAMEBUFF);
		connections_FP = fopen(filename, "w");
		print_iarray(connections_FP, affNeurons_IE[l], mp->nExcit, mp->nSynIE);
		fclose(connections_FP);
	}
	
	for (l=0; l<mp->nLayers; l++)
	{
		slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsII.dat", l);
		assert(slen < FNAMEBUFF);
		connections_FP = fopen(filename, "w");
		print_iarray(connections_FP, affNeurons_II[l], mp->nInhib, mp->nSynII);
		fclose(connections_FP);
	}
	
	//************** END OF FILE OUTPUT **************
	
	return; //void?
}*/


/*void wire_afferents(NEURON * n, int l, int * affNcnx_fE, int * affNcnx_lE, int * affNcnx_I)
{ // n_E, n_I, 
	int s, choice;
	
	//*** Ef synapses ***
	rands1_new(0,0,&iv,1);	// Initialise RNG (clear internal array of numbers drawn)
	n->lm1presyn_E = myalloc(n->nFAff_E * sizeof(NEURON *));
	n->FAffs_E = myalloc(n->nFAff_E * sizeof(AXON *));
	for (s=0; s<n->nFAff_E; s++) //if (n->nFAff_E > 0)
	{
		choice=rands1_new(0,mp->nExcit-1,&iv,1); // Sample from presyn neurons without replacement
		affNcnx_fE[s] = choice;
		n->lm1presyn_E[s] = &n_E[l-1][choice];
		n_E[l-1][choice].nFEff_E++;
	}
	
	//*** Lateral El synapses ***
	rands1_new(0,0,&iv,1);	// Initialise RNG (clear internal array of numbers drawn)
	n->lm0presyn_E = myalloc(n->nLAff_E * sizeof(NEURON *)); // Create array of NEURON pointers
	n->LAffs_E = myalloc(n->nLAff_E * sizeof(AXON *)); // Create array of AXON pointers
	for (s=0; s<n->nLAff_E; s++)
	{
		choice=rands1_new(0,mp->nExcit-1,&iv,1); // Sample from pre-synaptic neurons without replacement
		affNcnx_lE[s] = choice; //*** NEW ARRAY!
		n->lm0presyn_E[s] = &n_E[l][choice]; // Point to presynaptic neuron
		if (n->type == EXCIT)
			n_E[l][choice].nLEff_E++;
		else if (n->type == INHIB)
			n_E[l][choice].nLEff_I++;
	}
	
	//*** Lateral I Synapses ***
	rands1_new(0,0,&iv,1);	// Initialise RNG (clear internal array of numbers drawn)
	n->lm0presyn_I = myalloc(n->nLAff_I * sizeof(NEURON *)); // Create array of NEURON pointers
	n->LAffs_I = myalloc(n->nLAff_I * sizeof(AXON *)); // Create array of AXON pointers
	for (s=0; s<n->nLAff_I; s++)
	{
		choice=rands1_new(0,mp->nInhib-1,&iv,1); // Sample from pre-synaptic neurons without replacement
		affNcnx_I[s] = choice;
		n->lm0presyn_I[s] = &n_I[l][choice]; // Point to presynaptic neuron
		if (n->type == EXCIT)
			n_I[l][choice].nLEff_E++;
		else if (n->type == INHIB)
			n_I[l][choice].nLEff_I++;
	}
	return;
}*/


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

void create_axons(NEURON * n)
{
	/*** Axonal delays specified in seconds must be converted to timesteps ***/
	tstep delay = 0;
	//int nBins = 0;
	int span = 0;
	int s = 0;
	switch (mp->axonDelay) 
	{
		case MinD:
			delay = 1;
			//nBins = 1;
			break;
		case ConstD:
			delay = round(mp->d_const/mp->DT);
			delay = (!delay) ? 1 : delay;
			//nBins = ceil(mp->d_const/mp->refract);
			//nBins = (!nBins) ? 1 : nBins;
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
		switch (mp->axonDelay) 
		{
			case UniformD:
				//delay = round(((ran3(&idum)*span)+mp->d_min)/DT);
				delay = round(((gsl_rng_uniform(mSeed)*span)+mp->d_min)/mp->DT);
				//delay = (!delay) ? 1 : delay;
				//nbins = ceil((delay*DT)/mp->refract);
				break;
			case GaussD:
				//delay = round((mp->d_mean + mp->d_sd*gasdev(&idum))/DT);
				delay = round((mp->d_mean + gsl_ran_gaussian(mSeed,mp->d_sd))/mp->DT);
				//delay = (!delay) ? 1 : delay;
				//nbins = ceil((delay*DT)/mp->refract);
				break;
			case SOMD:
				delay = 1;
				break;
			default:
				break;
		}
		n->FEffs_E[s].delay = (!delay) ? 1 : delay;
		n->FEffs_E[s].nBins = ceil((n->FEffs_E[s].delay*mp->DT)/mp->refract);
		n->FEffs_E[s].queue = myalloc(n->FEffs_E[s].nBins * sizeof(tstep));
		init_queue(&(n->FEffs_E[s]));
	}
	for (s=0; s<n->nLEff_E; s++)
	{
		switch (mp->axonDelay) 
		{
			case UniformD:
				//delay = round(((ran3(&idum)*span)+mp->d_min)/DT);
				delay = round(((gsl_rng_uniform(mSeed)*span)+mp->d_min)/mp->DT);
				break;
			case GaussD:
				//delay = round((mp->d_mean + mp->d_sd*gasdev(&idum))/DT);
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
		n->LEffs_E[s].delay = (!delay) ? 1 : delay;
		n->LEffs_E[s].nBins = ceil((n->LEffs_E[s].delay*mp->DT)/mp->refract);
		n->LEffs_E[s].queue = myalloc(n->LEffs_E[s].nBins * sizeof(tstep));
		init_queue(&(n->LEffs_E[s]));
	}
	for (s=0; s<n->nLEff_I; s++)
	{
		switch (mp->axonDelay) 
		{
			case UniformD:
				//delay = round(((ran3(&idum)*span)+mp->d_min)/DT);
				delay = round(((gsl_rng_uniform(mSeed)*span)+mp->d_min)/mp->DT);
				break;
			case GaussD:
				//delay = round((mp->d_mean + mp->d_sd*gasdev(&idum))/DT);
				delay = round((mp->d_mean + gsl_ran_gaussian(mSeed,mp->d_sd))/mp->DT);
				break;
			case SOMD:
				delay = 1;
				break;
			default:
				break;
		}
		n->LEffs_I[s].delay = (!delay) ? 1 : delay;
		n->LEffs_I[s].nBins = ceil((n->LEffs_I[s].delay*mp->DT)/mp->refract);
		n->LEffs_I[s].queue = myalloc(n->LEffs_I[s].nBins * sizeof(tstep));
		init_queue(&(n->LEffs_I[s]));
	}
	return;
}

void init_network(int regime)
{
	int s, n, l;
	int r; //, syn, t;
	//int choice;
	
	if (regime == Learning) // Only initialise before training, not testing
	{
		if (mp->nRecords)
		{
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
				assert((mp->vRecords[l] >= 0) && (mp->vRecords[l] <= mp->vExcit[l]));
				chosen = myalloc(mp->vRecords[l] * sizeof(*chosen));
				choices = myalloc(mp->vExcit[l] * sizeof(*choices));
				for (n=0; n<mp->vExcit[l]; n++)
					choices[n] = n;
				//choice=rands1_new(0,0,&iv,1);	// Initialise RNG
				gsl_ran_choose(mSeed, chosen, mp->vRecords[l], choices, mp->vExcit[l], sizeof(*choices));
				for (r=0; r<mp->vRecords[l]; r++)
				{
					// Randomly set NRECORDS flags
					//choice=rands1_new(0,mp->vExcit[l]-1,&iv,1);
					//n_E[l][choice].rec_flag = true;
					n_E[l][chosen[r]].rec_flag = true;
					printf("\tLayer %d, Record %d Assigned nID: %d\n",l,r+1,chosen[r]);
				}
				slen = snprintf(rString, BUFSIZ, "MP.Records{%d}", l+1);
				assert(slen < BUFSIZ);
				printIntArray(rFile, rString, chosen, mp->vRecords[l]); 
				myfree(chosen);
				myfree(choices);
			}
			//myfree(chosen);
			fclose(rFile); // Close Records file
			/* Create recording structures - (mp->TotalMS+1) so that initial conditions are recorded */
			for (l=0; l<mp->nLayers; l++)
				for (n=0; n<mp->vExcit[l]; n++)
					if (n_E[l][n].rec_flag)
					{
						n_E[l][n].rec = myalloc(sizeof(RECORD));
						n_E[l][n].rec->bin = 0;
						n_E[l][n].rec->cellV = myalloc((mp->TotalMS+1) * sizeof(*(n_E[l][n].rec->cellV))); // get_2D_farray(mp->loops, mp->TotalMS, 0.0);
						memset(n_E[l][n].rec->cellV, 0, (mp->TotalMS+1) * sizeof(*(n_E[l][n].rec->cellV)));
						if (mp->adaptation)
						{
							n_E[l][n].rec->cellcCa = myalloc((mp->TotalMS+1) * sizeof(*(n_E[l][n].rec->cellcCa)));
							memset(n_E[l][n].rec->cellcCa, 0, (mp->TotalMS+1) * sizeof(*(n_E[l][n].rec->cellcCa)));
						}
						else
							n_E[l][n].rec->cellcCa = NULL;
						n_E[l][n].rec->cellD = myalloc((mp->TotalMS+1) * sizeof(*(n_E[l][n].rec->cellD))); // get_2D_farray(mp->loops, mp->TotalMS, 0.0);
						memset(n_E[l][n].rec->cellD, 0, (mp->TotalMS+1) * sizeof(*(n_E[l][n].rec->cellD)));
						if (l>0) // Here the synaptic variables are stored with the post-synaptic neuron
						{
							n_E[l][n].rec->FSynC = get_2D_farray(n_E[l][n].nFAff_E, (mp->TotalMS+1), 0.0); // get_3D_farray(mp->loops, n_E[l][n].nFAff_E, mp->TotalMS, 0.0);
							n_E[l][n].rec->FSynG = get_2D_farray(n_E[l][n].nFAff_E, (mp->TotalMS+1), 0.0);
							n_E[l][n].rec->FSynDG = get_2D_farray(n_E[l][n].nFAff_E, (mp->TotalMS+1), 0.0);
						}
						if (mp->trainElE && mp->train) // Create recording structures for ElE variables
						{
							n_E[l][n].rec->LSynC = get_2D_farray(n_E[l][n].nLAff_E, (mp->TotalMS+1), 0.0);
							n_E[l][n].rec->LSynG = get_2D_farray(n_E[l][n].nLAff_E, (mp->TotalMS+1), 0.0);
							n_E[l][n].rec->LSynDG = get_2D_farray(n_E[l][n].nLAff_E, (mp->TotalMS+1), 0.0);
						}
					}
		} // End of Record initialisation
		
		/* E_ Synaptic weights */
		/*for (l=0; l<mp->nWLayers; l++)
			for (n=0; n<mp->vExcit[l]; n++)
				for (s=0; s<n_E[l][n].nFEff_E; s++)
					n_E[l][n].FEffs_E[s].delta_g = n_E[l][n].FEffs_E[s].delta_g_tm1 = gsl_rng_uniform(mSeed); //ran3(&idum);*/

		switch (mp->initEfE) {
			case Constant:
				for (l=0; l<mp->nWLayers; l++)
					for (n=0; n<mp->vExcit[l]; n++)
						for (s=0; s<n_E[l][n].nFEff_E; s++)
							n_E[l][n].FEffs_E[s].delta_g = n_E[l][n].FEffs_E[s].delta_g_tm1 = mp->DgEfE;
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
				exit_error("init_network", "Feed-forward convergent connections not yet implemented!");
				//fprintf(stderr, "Warning: Feed-forward convergent connections not yet implemented!\n");
				break;
			default:
				exit_error("init_network", "Illegal EfE weight initialisation arguement!");
				break;
		}
		/*if (mp->SOM)
		{
			float distance = 0.0;
			float phiScale = mp->Dg_ElE/(mp->SOMsigE*sqrt(2*M_PI)); 
			for (l=0; l<mp->nLayers; l++)
				for (n=0; n<mp->vExcit[l]; n++)
					for (s=0; s<n_E[l][n].nLEff_E; s++)
					{
						distance = readLowTriF(distE, l, n, n_E[l][n].lp0postsyn_E[s]->n); // distE[l][n][n_E[l][n].lp0postsyn_E[s]->n];
						n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = phiScale * exp(-pow(distance,2)/(2*pow(mp->SOMsigE,2)));
					}
			
//			if (!mp->SOMinput)
//				for (n=0; n<mp->sInputs; n++)
//					for (s=0; s<n_E[0][n].nLEff_E; s++)
//						n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = mp->Dg_ElE; 
		}
		else
		{
			for (l=0; l<mp->nLayers; l++)
				for (n=0; n<mp->vExcit[l]; n++)
					for (s=0; s<n_E[l][n].nLEff_E; s++)
						n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = mp->Dg_ElE; // * gsl_rng_uniform(mSeed); //ran3(&idum);
		}*/
		
		/*if (mp->SOM) // pCnxElE ??? nLEff_E == 0 if pCnxElE == 0.0
		{*/
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
							n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = mp->DgElE;
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
				phiScale = mp->DgElE/(mp->SOMsigE*sqrt(2*M_PI));
				for (l=0; l<mp->nLayers; l++)
					for (n=0; n<mp->vExcit[l]; n++)
						for (s=0; s<n_E[l][n].nLEff_E; s++)
						{
							distance = readLowTriF(distE, l, n, n_E[l][n].lp0postsyn_E[s]->n); // distE[l][n][n_E[l][n].lp0postsyn_E[s]->n];
							n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = phiScale * exp(-pow(distance,2)/(2*pow(mp->SOMsigE,2)));
						}
				break;					
			default:
				exit_error("init_network", "Illegal ElE weight initialisation arguement!");
				break;
		}
		/*}*/
		
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
		
	} // End of if (regime==1) clause
	

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
			
			if (regime != Settle) // Reinitialize spike buffers between epochs and phases
			{
				n_E[l][n].spkbin = 0;
				//bins = (mp->useFilteredImages && l==0) ? mp->inpSpkBuff : mp->spkBuffer;
				memset(n_E[l][n].spikeTimes, 0, mp->spkBuffer*sizeof(n_E[l][n].spikeTimes[0])); //[0]?
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
			
			if (regime != Settle)
			{
				n_I[l][n].spkbin = 0;
				memset(n_I[l][n].spikeTimes, 0, mp->spkBuffer*sizeof(n_I[l][n].spikeTimes[0]));
			}
		}
	}
	
	/* Recording structures */ // Moved to simulatePhase
/*	if (mp->nRecordsPL && regime != Settle)
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for
			for (n=0; n<mp->vExcit[l]; n++)
				if (n_E[l][n].rec_flag)
					n_E[l][n].rec->bin = 0;
		}*/
		
	} // End of parallel region

	return;
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
	mp->newTestSet = false;
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
	if (mp->newTestSet)
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

void gen_stimuli(bool rep, STIMULI * stim, PARAMS * mp)
{
	// Assumes 1D inputs
	// Place array arguements in a patterns structure
	int trans;
	int n, p;
	int * choices = NULL;
	int * chosen = NULL;
	int block = 0;
	char stimStr[BUFSIZ];
	int slen = 0;

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
		int ind = 0;
		int c = 0;
		int m = 0;
		int q = 0;
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
	}
	else	// Training stimuli presented individually
		stim->trn_stimuli = stim->tst_stimuli;

	
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

int genShuffles(STIMULI * stim, PARAMS * mp)
{
	int loop = 0;
	int p = 0;
	int trans = 0;
	//int i = 0;
	int count = 0;
	int * choices = NULL;
	//int * chosen = NULL;
	
	if (mp->randStimOrder)
	{
		choices = myalloc(mp->nStimuli * sizeof(int));
		for (p=0; p<mp->nStimuli; p++)
			choices[p] = p;
		stim->stimShuffle = get_2D_iarray(mp->loops, mp->nStimuli, 0);
		for (loop=0; loop<mp->loops; loop++)
		{
			/*i=rands1_new(0,0,&iv,1);	// Initialise RNG (clear internal array of numbers drawn)
			for (p=0; p<mp->nStimuli; p++)
				stim->stimShuffle[loop][p] = rands1_new(0,mp->nStimuli-1,&iv,1);*/
			gsl_ran_shuffle(mSeed, choices, mp->nStimuli, sizeof(int));
			memcpy(stim->stimShuffle[loop], choices, mp->nStimuli * sizeof(int));
		}
		count++;
		myfree(choices);
	}
	
	if (mp->randTransOrder)
	{
		choices = myalloc(mp->nTransPS * sizeof(int));
		for (trans=0; trans<mp->nTransPS; trans++)
			choices[trans] = trans;
		stim->transShuffle = get_3D_iarray(mp->loops, mp->nStimuli, mp->nTransPS, 0);
		for (loop=0; loop<mp->loops; loop++)
		{
			for (p=0; p<mp->nStimuli; p++)
			{
				/*i=rands1_new(0,0,&iv,1);	// Initialise RNG (clear internal array of numbers drawn)
				 for (trans=0; trans<mp->nTransPS; trans++)
				 stim->transShuffle[loop][p][trans] = rands1_new(0,mp->nTransPS-1,&iv,1);*/
				gsl_ran_shuffle(mSeed, choices, mp->nTransPS, sizeof(int));
				memcpy(stim->transShuffle[loop][p], choices, mp->nTransPS * sizeof(int));
			}
		}
		count++;
		myfree(choices);
	}	
	return count;
}

void calc_input(int loop, int pat, int trans, STIMULI * stim, float ** input, int regime) // return int * input?
{	// Assumes all stimuli translate
	float *****patArray = NULL;
	if (mp->useFilteredImages) // Select stimuli from array of filtered images
	{
		switch (regime) 
		{
			case 0:
				*input = ****(stim->tstImages[pat][trans]);	//memcpy(input, ****(stim->tstImages[pat][trans]), mp->sInputs*sizeof(input[0])); // Use pointer not memcpy
				break;
				
			case 1:
				pat = (mp->randStimOrder) ? stim->stimShuffle[loop][pat] : pat;
				trans = (mp->randTransOrder) ? stim->transShuffle[loop][pat][trans] : trans;
				patArray = (mp->newTestSet) ? stim->tstImages[pat][trans] : stim->trnImages[pat][trans];
				*input = ****patArray; //memcpy(input, ****patArray, mp->sInputs*sizeof(input[0]));
				break;
		}
	}
	else
	{
		switch (regime)
		{
			case 0: // Testing stimuli
				*input = stim->tst_stimuli[pat][trans]; //Test me!	//memcpy(input, stim->tst_stimuli[pat][trans], mp->nExcit*sizeof(input[0])); // Copy row of tst_stimuli to input
				break;
				
			case 1: // Training stimuli
				pat = (mp->randStimOrder) ? stim->stimShuffle[loop][pat] : pat;
				trans = (mp->randTransOrder) ? stim->transShuffle[loop][pat][trans] : trans;
				*input = stim->trn_stimuli[pat][trans];	//memcpy(input, stim->trn_stimuli[pat][trans], mp->nExcit*sizeof(input[0]));
				break;
		}
	}

	/*// Move print statements outside function?
	if ( trans==0 || (mp->interleaveTrans && regime) )
		printf("\t\tPresenting stimulus %d/%d...\n", (pat+1), mp->nStimuli);
	printf("\t\t\tTransform %d/%d...\n", (trans+1), mp->nTransPS);*/

	return; // void;
}

void simulatePhase(PHASE sPhase, STIMULI * stim)
{
	int * o = NULL;
	int * i = NULL;
	int oCount = 0;
	int iCount = 0;
	int p=0;
	int tr=0;
	int nStimuli = 0;
	int nTrans = 0;
	float transP = 0.0;
	tstep t_start = 0;
	tstep t_end = 0;
	int loop, l, wl, n, syn;
	int nLoops = 0;
	int result = 0;
	REGIMETYPE regime = NoLearning;
	int slen = 0;
	char phaseString[FNAMEBUFF];
	char filename[FNAMEBUFF];
	FILE * excitOutput;
	FILE * inhibOutput;
	FILE * weightOutput;

	FILE * recFile;
	char prefix[FNAMEBUFF];
	/*FILE * rCellVout;
	FILE * rDout;
	FILE * rCout;
	FILE * rGout;
	FILE * rWeightOut;*/
	
	bool reverse = false;
	
	double percentage = 0.0;
	
	float * input = NULL; //myalloc(mp->sInputs * sizeof(*input));
	
	switch (sPhase) 
	{
		case PreTraining:
			strncpy(phaseString, "Pretraining", FNAMEBUFF);
			strncpy(prefix, "Rpt", FNAMEBUFF);
			nStimuli = mp->nTestStimuli;
			nTrans = mp->nTestTransPS;
			transP = mp->transP_Test;
			regime = NoLearning;
			nLoops = 1;
			break;
		case Training:
			strncpy(phaseString, "Training", FNAMEBUFF);
			nStimuli = mp->nStimuli;
			nTrans = mp->nTransPS;
			transP = mp->transP_Train;
			regime = Learning; // Weights are initialised after building the network
			nLoops = mp->loops;
			break;
		case Testing:
			strncpy(phaseString, "Testing", FNAMEBUFF);
			strncpy(prefix, "R", FNAMEBUFF);
			nStimuli = mp->nTestStimuli;
			nTrans = mp->nTestTransPS;
			transP = mp->transP_Test;
			regime = NoLearning;
			nLoops = 1;
			break;
		default:
			exit_error("simulatePhase", "Unknown phase!");
			break;
	}
	
	if (mp->interleaveTrans && sPhase == Training) // Could move to switch
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
	
	printf("\tNow beginning %s phase...\n",phaseString);
	if (sPhase == Training)
		printf("\tSimulating %2.2f s per loop for %d loop%s.\n",mp->TotalTime,nLoops,(nLoops==1)?"":"s"); // replace TotalTime?
	for (loop=0; loop<nLoops; loop++)
	{
		init_network(NoLearning); // Reset V's, g's, C's, D's and spike buffers - NoLearning even for Training! 
		if (sPhase == Training)
			printf("\tLoop %d/%d\n", loop+1, mp->loops);
		for (*o=0; *o<oCount; (*o)++)
		{
			if (sPhase==Training && !mp->interleaveTrans)
			{
				if (mp->trainPause)
					init_network(Settle);
				if (mp->randTransDirection) // Generate [0,n-1] with equal probability
					reverse = (gsl_rng_uniform_int(mSeed, 2)) ? true : false;
			}

			for (*i=0; *i<iCount; (*i)++)
			{
				if (sPhase != Training) //Testing only
					init_network(Settle);
				calc_input(loop, p, ((reverse) ? (nTrans-1)-tr : tr), stim, &input, regime);
				if ( tr==0 || (mp->interleaveTrans && regime) )
					printf("\t\tPresenting stimulus %d/%d...\n", (p+1), mp->nStimuli);
				printf("\r\t\t\tTransform %d/%d...\n", (reverse ? nTrans-tr : tr+1), mp->nTransPS);
				t_start = round((*i + (*o * iCount)) * transP * ceil(1/mp->DT));
				t_end = t_start + round(transP * ceil(1/mp->DT));
#if DEBUG > 1 // Level 2
				fprintf(stderr, "Updating network from timestep %d to %d.\n", t_start, t_end-1); //%lld
#endif
				
#if DEBUG > 3
				print_frow(stderr, input, mp->sInputs);
#endif
				updateNetwork(t_start, t_end, input, regime); //loop, 
				
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
				
				if (SIM.Xgrid)
				{
					SIM.tally += round(transP/mp->DT); //transP * ceil(1/mp->DT);
					percentage = (100.0*SIM.tally)/SIM.totTS;
					//printf("\n");
					printf("<xgrid>{control = statusUpdate; percentDone = %.0lf; }</xgrid>\n", percentage);
					fflush(stdout);
				}
			}
		}
		printf("\t%s complete!\n",phaseString);
		
		printf("\tSaving results..."); // Output results to dat files
		
		for (l=0; l<mp->nLayers; l++) // Rethink l init //l=(sPhase==Training)?0:1;
		{
			if (sPhase == PreTraining) // Only print out Layer 1 onwards
				slen = snprintf(filename, FNAMEBUFF, "ptL%dExcitSpikes.dat", l);
			else if (sPhase == Training) /* Save input layer spikes after each epoch */
				slen = snprintf(filename, FNAMEBUFF, "E%dL%dExcitSpikes.dat", loop, l);
			else if (sPhase == Testing) // Only print out Layer 1 onwards
				slen = snprintf(filename, FNAMEBUFF, "L%dExcitSpikes.dat", l);

			assert(slen < FNAMEBUFF);
			excitOutput = myfopen(filename, "w");
			for (n=0; n<mp->vExcit[l]; n++)
			{
				fprintf(excitOutput, "%d\t", n_E[l][n].spkbin); // Print number of spikes first (in case any at t=0)
				print_irow(excitOutput, (int*) n_E[l][n].spikeTimes, n_E[l][n].spkbin); // spikeTimes[0] = -BIG?
			}

			fclose(excitOutput);
		}
		
		if (sPhase == Testing) // Save Inhibitory spikes for all layers
		{
			for (l=0; l<mp->nLayers; l++)
			{
				slen = snprintf(filename, FNAMEBUFF, "L%dInhibSpikes.dat", l);
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
		
		if (sPhase != Training) // Save weights for EFE synapses - Pretraining and Testing
		{
			for (wl=1; wl<mp->nLayers; wl++)
			{
				if (sPhase == PreTraining)
					slen = snprintf(filename, FNAMEBUFF, "ptL%dweightsEfE.dat", wl);
				else if (sPhase == Testing)
					slen = snprintf(filename, FNAMEBUFF, "L%dweightsEfE.dat", wl);
				
				assert(slen < FNAMEBUFF);
				weightOutput = myfopen(filename, "w");
				for (n=0; n<mp->vExcit[wl]; n++)
				{
					for (syn=0; syn<n_E[wl][n].nFAff_E; syn++)
						fprintf(weightOutput, "%f\t", n_E[wl][n].FAffs_E[syn]->delta_g);
					fprintf(weightOutput, "\n");
				}
				fclose(weightOutput);
			}
			
			if (mp->trainElE)
			{
				for (l=0; l<mp->nLayers; l++)
				{
					if (sPhase == PreTraining)
						slen = snprintf(filename, FNAMEBUFF, "ptL%dweightsElE.dat", wl);
					else if (sPhase == Testing)
						slen = snprintf(filename, FNAMEBUFF, "L%dweightsElE.dat", wl);
					assert(slen < FNAMEBUFF);
					weightOutput = myfopen(filename, "w");
					for (n=0; n<mp->vExcit[l]; n++)
					{
						for (syn=0; syn<n_E[l][n].nLAff_E; syn++)
							fprintf(weightOutput, "%f\t", n_E[l][n].LAffs_E[syn]->delta_g);
						fprintf(weightOutput, "\n");
					}
					fclose(weightOutput);
				}
			}
		}
		
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
						// print_farray(rCellVout, n_E[l][n].rec->cellV, mp->loops, mp->TotalMS); //fwrite(n_E[l][n].rec->cellV, sizeof(float), mp->loops*mp->TotalMS, r_cellV_ptr); // See nifty trick #1
						fclose(recFile);
						memset(n_E[l][n].rec->cellV, (mp->TotalMS+1), sizeof(*(n_E[l][n].rec->cellV)));
						
						if (mp->adaptation)
						{
							slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dcCa.dat", prefix, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w");
							print_frow(recFile, n_E[l][n].rec->cellcCa, n_E[l][n].rec->bin);
							fclose(recFile);
							memset(n_E[l][n].rec->cellcCa, (mp->TotalMS+1), sizeof(*(n_E[l][n].rec->cellcCa)));
						}
						
						slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dD.dat", prefix, l, n);
						assert(slen && slen < FNAMEBUFF);
						recFile = myfopen(filename, "w"); //r_D_ptr = fopen(filename, "wb");
						print_frow(recFile, n_E[l][n].rec->cellD, n_E[l][n].rec->bin); // (rDout, n_E[l][n].rec->cellD, mp->loops, mp->TotalMS);
						fclose(recFile);
						memset(n_E[l][n].rec->cellD,(mp->TotalMS+1), sizeof(*(n_E[l][n].rec->cellD)));
						
						if (l>0) // The presynaptic cell's values are attached to each record
						{
							slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dFAffC.dat", prefix, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w"); //r_C_ptr = fopen(filename, "wb");
							//for (loop=0; loop<mp->loops; loop++)
							print_farray(recFile, n_E[l][n].rec->FSynC, n_E[l][n].nFAff_E, n_E[l][n].rec->bin); //print_farray(recFile, n_E[l][n].rec->SynC[loop], n_E[l][n].nFAff_E, mp->TotalMS);
							fclose(recFile);
							memset(n_E[l][n].rec->FSynC[0], n_E[l][n].nFAff_E*(mp->TotalMS+1), sizeof(*(n_E[l][n].rec->FSynC)));
							
							slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dFAffg.dat", prefix, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w");
							//for (loop=0; loop<mp->loops; loop++)
							print_farray(recFile, n_E[l][n].rec->FSynG, n_E[l][n].nFAff_E, n_E[l][n].rec->bin);
							fclose(recFile);
							memset(n_E[l][n].rec->FSynG[0], n_E[l][n].nFAff_E*(mp->TotalMS+1), sizeof(*(n_E[l][n].rec->FSynG)));
							
							slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dFAffdg.dat", prefix, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w");
							//for (loop=0; loop<mp->loops; loop++)
							print_farray(recFile, n_E[l][n].rec->FSynDG, n_E[l][n].nFAff_E, n_E[l][n].rec->bin);
							fclose(recFile);
							memset(n_E[l][n].rec->FSynDG[0], n_E[l][n].nFAff_E*(mp->TotalMS+1), sizeof(*(n_E[l][n].rec->FSynDG)));
						}
						
						if (mp->trainElE && mp->train)
						{
							slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dLAffC.dat", prefix, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w"); 
							print_farray(recFile, n_E[l][n].rec->LSynC, n_E[l][n].nLAff_E, n_E[l][n].rec->bin); 
							fclose(recFile);
							memset(n_E[l][n].rec->LSynC[0], n_E[l][n].nLAff_E*(mp->TotalMS+1), sizeof(*(n_E[l][n].rec->LSynC)));
							
							slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dLAffg.dat", prefix, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w");
							print_farray(recFile, n_E[l][n].rec->LSynG, n_E[l][n].nLAff_E, n_E[l][n].rec->bin);
							fclose(recFile);
							memset(n_E[l][n].rec->LSynG[0], n_E[l][n].nLAff_E*(mp->TotalMS+1), sizeof(*(n_E[l][n].rec->LSynG)));
							
							slen = snprintf(filename, FNAMEBUFF, "%sL%dN%dLAffdg.dat", prefix, l, n);
							assert(slen && slen < FNAMEBUFF);
							recFile = myfopen(filename, "w");
							print_farray(recFile, n_E[l][n].rec->LSynDG, n_E[l][n].nLAff_E, n_E[l][n].rec->bin);
							fclose(recFile);
							memset(n_E[l][n].rec->LSynDG[0], n_E[l][n].nLAff_E*(mp->TotalMS+1), sizeof(*(n_E[l][n].rec->LSynDG)));
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
		
	if (!t_start && mp->nRecords) //(t_start==0) Save intials neuron states (then every ms after - see below)
	{ // Dump records to file after each loop
//#pragma omp parallel private(bin, syn, n, l)
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait 
			for (n=0; n<mp->vExcit[l]; n++) //for (r=0; r<NRECORDS_PL; r++)
				if (n_E[l][n].rec_flag) // Can use tm1 variables since they are initialised equal to tm0 variables
				{
					bin = n_E[l][n].rec->bin++;
					n_E[l][n].rec->cellV[bin] = n_E[l][n].V_tm1;
					n_E[l][n].rec->cellD[bin] = n_E[l][n].D_tm1;
					if (mp->adaptation)
						n_E[l][n].rec->cellcCa[bin] = n_E[l][n].cCa_tm1;
					for (syn=0; syn<n_E[l][n].nFAff_E; syn++) // nFAff_E == 0 for l=0
					{
						n_E[l][n].rec->FSynC[syn][bin] = n_E[l][n].FAffs_E[syn]->C_tm1;
						n_E[l][n].rec->FSynG[syn][bin] = n_E[l][n].FAffs_E[syn]->g_tm1;
						n_E[l][n].rec->FSynDG[syn][bin] = n_E[l][n].FAffs_E[syn]->delta_g_tm1;
					}
					if (mp->trainElE && mp->train)
					{
						for (syn=0; syn<n_E[l][n].nLAff_E; syn++) // nFAff_E == 0 for l=0
						{
							n_E[l][n].rec->LSynC[syn][bin] = n_E[l][n].LAffs_E[syn]->C_tm1;
							n_E[l][n].rec->LSynG[syn][bin] = n_E[l][n].LAffs_E[syn]->g_tm1;
							n_E[l][n].rec->LSynDG[syn][bin] = n_E[l][n].LAffs_E[syn]->delta_g_tm1;
						}
					}
				}
		}
//#pragma omp barrier // Only need a barrier before solution variables are copied
	}


	for (t=t_start; t<t_end; t++)
	{
	
#if DEBUG>2
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
		
		if (regime==Learning) // Learning
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
		
		
		/* Copy solution variables to _tm1 counterparts and reset spike flags / axons */
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait //schedule(runtime)// private(l, n, s) //ok
			for (n=0; n<mp->vExcit[l]; n++)
			{
				//n_E[l][n].fired = 0;
				n_E[l][n].V_tm1 = n_E[l][n].V;
				n_E[l][n].D_tm1 = n_E[l][n].D;
				if (mp->adaptation)
					n_E[l][n].cCa_tm1 = n_E[l][n].cCa;
				for (syn=0; syn<n_E[l][n].nFEff_E; syn++)
				{
					if (t == next_spike(&n_E[l][n].FEffs_E[syn]))
						dequeue(&n_E[l][n].FEffs_E[syn]);
					n_E[l][n].FEffs_E[syn].g_tm1 = n_E[l][n].FEffs_E[syn].g;
					n_E[l][n].FEffs_E[syn].delta_g_tm1 = n_E[l][n].FEffs_E[syn].delta_g;
					n_E[l][n].FEffs_E[syn].C_tm1 = n_E[l][n].FEffs_E[syn].C;
				}
				for (syn=0; syn<n_E[l][n].nLEff_E; syn++)
				{
					if (t == next_spike(&n_E[l][n].LEffs_E[syn]))
						dequeue(&n_E[l][n].LEffs_E[syn]);
					n_E[l][n].LEffs_E[syn].g_tm1 = n_E[l][n].LEffs_E[syn].g;
					if (mp->trainElE)
					{
						n_E[l][n].LEffs_E[syn].delta_g_tm1 = n_E[l][n].LEffs_E[syn].delta_g;
						n_E[l][n].LEffs_E[syn].C_tm1 = n_E[l][n].LEffs_E[syn].C;
					}
				}
				for (syn=0; syn<n_E[l][n].nLEff_I; syn++)
				{
					if (t == next_spike(&n_E[l][n].LEffs_I[syn]))
						dequeue(&n_E[l][n].LEffs_I[syn]);
					n_E[l][n].LEffs_I[syn].g_tm1 = n_E[l][n].LEffs_I[syn].g;
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
						n_E[l][n].rec->cellD[bin] = n_E[l][n].D;
						if (mp->adaptation)
							n_E[l][n].rec->cellcCa[bin] = n_E[l][n].cCa;
						for (syn=0; syn<n_E[l][n].nFAff_E; syn++) // nFAff_E == 0 for l=0
						{
							n_E[l][n].rec->FSynC[syn][bin] = n_E[l][n].FAffs_E[syn]->C;
							n_E[l][n].rec->FSynG[syn][bin] = n_E[l][n].FAffs_E[syn]->g;
							n_E[l][n].rec->FSynDG[syn][bin] = n_E[l][n].FAffs_E[syn]->delta_g;
						}
						if (mp->trainElE && mp->train)
						{
							for (syn=0; syn<n_E[l][n].nLAff_E; syn++) // nFAff_E == 0 for l=0
							{
								n_E[l][n].rec->LSynC[syn][bin] = n_E[l][n].LAffs_E[syn]->C;
								n_E[l][n].rec->LSynG[syn][bin] = n_E[l][n].LAffs_E[syn]->g;
								n_E[l][n].rec->LSynDG[syn][bin] = n_E[l][n].LAffs_E[syn]->delta_g;
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
	/* Update E synapses */
	scale = mp->DgEfE * mp->gMax;
	for (syn=0; syn<n->nFAff_E; syn++) // Make private
	{
		impulse = (next_spike(n->FAffs_E[syn]) == t) ? n->FAffs_E[syn]->delta_g_tm1 * scale : 0.0;
		n->FAffs_E[syn]->g += (impulse - (decay_E * n->FAffs_E[syn]->g_tm1));
	}
	
	scale = (mp->trainElE) ? mp->DgElE * mp->gMax : mp->gMax;
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
		// Shunting inhibition ?
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
	for (syn=0; syn<n->nFAff_E; syn++)
	{
		LTD = (t == next_spike(n->FAffs_E[syn])) ? n->FAffs_E[syn]->delta_g_tm1 * n->D_tm1 : 0.0; //tm0?
		LTP = (t == n->lastSpike) ? (1 - n->FAffs_E[syn]->delta_g_tm1) * n->FAffs_E[syn]->C_tm1 : 0.0; //tm0?
		n->FAffs_E[syn]->delta_g += (LTP - LTD) * mp->learnR; //*DT/TAU_DG;
	}
	
	if (mp->trainElE) // 
	{
		for (syn=0; syn<n->nLAff_E; syn++)
		{
			LTD = (t == next_spike(n->LAffs_E[syn])) ? n->LAffs_E[syn]->delta_g_tm1 * n->D_tm1 : 0.0; //tm0?
			LTP = (t == n->lastSpike) ? (1 - n->LAffs_E[syn]->delta_g_tm1) * n->LAffs_E[syn]->C_tm1 : 0.0; //tm0?
			n->LAffs_E[syn]->delta_g += (LTP - LTD) * mp->learnR; //*DT/TAU_DG;
		}
	}
	return;
}
	
inline void update_C(NEURON * n, tstep t, float decayRate)
{
	// -->|| Update C for current neuron's outgoing synapses (skip last layer)
	int syn;
	float impulse;
	//float C_tm1 = 0.0;
	//float decayRate = mp->DT/mp->tauC;
	// Loop over efferent synapses
	/*for (syn=0; syn<n->nFEff_E; syn++)
	{
		impulse = (t == next_spike(&n->FEffs_E[syn])) ? mp->alphaC : 0.0;
		n->FEffs_E[syn]->C += (impulse * (1 - n->FEffs_E[syn]->C_tm1) - (n->FEffs_E[syn]->C_tm1 * DT/mp->tauC));
	}*/
	// Loop over afferent synapses (skip first layer)
	for (syn=0; syn<n->nFAff_E; syn++)
	{
		impulse = (t == next_spike(n->FAffs_E[syn])) ? mp->alphaC * (1 - n->FAffs_E[syn]->C_tm1) : 0.0;
		n->FAffs_E[syn]->C += (impulse - (n->FAffs_E[syn]->C_tm1 * decayRate));
	}
	
	if (mp->trainElE) // Update Excitatory lateral Afferent synapses
	{
		for (syn=0; syn<n->nLAff_E; syn++)
		{
			//C_tm1 = n->LAffs_E[syn]->C_tm1;
			//n->LAffs_E[syn]->C += ((t==next_spike(n->LAffs_E[syn])) ? (mp->alphaC*(1-C_tm1)) : 0.0) - (C_tm1*decayRate);
			impulse = (t == next_spike(n->LAffs_E[syn])) ? (mp->alphaC * (1 - n->LAffs_E[syn]->C_tm1)) : 0.0;
			n->LAffs_E[syn]->C += (impulse - (n->LAffs_E[syn]->C_tm1 * decayRate));
			/*if (n->LAffs_E[syn]->C>1.0)
			{
				printf("n=%d; syn=%d; C_tm1=%f; C=%f; ",n->n,syn,n->LAffs_E[syn]->C_tm1,n->LAffs_E[syn]->C);
				printf("t=%d",t);
				fflush(stdout);
				printf("alpha=%f; ",mp->alphaC);
			}*/
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
	assert(a->count < a->nBins);
	a->last = (a->last+1) % a->nBins;
	a->queue[ a->last ] = t + a->delay;
	a->count++;
	return;
}

inline int dequeue(AXON *a)
{
	int t;
	assert(a->count > 0);
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
