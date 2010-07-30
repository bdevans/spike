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
	
#if DEBUG > 1
	fprintf(stderr, "*** DT = %f ms ***\n", DT*1000);
	fprintf(stderr, "NEURON:\t%ld bytes.\n", sizeof(NEURON));
	fprintf(stderr, "AXON:\t%ld bytes.\n", sizeof(AXON));
	fprintf(stderr, "*:\t%ld bytes.\n", sizeof(NEURON*));
	fprintf(stderr, "float:\t%ld bytes.\n", sizeof(float));
	fprintf(stderr, "tstep:\t%ld bytes.\n", sizeof(tstep));
	fprintf(stderr, "int:\t%ld bytes.\n", sizeof(int));
	
	int MB = 1024*1024;
	int l=0;
	fprintf(stderr, "\nLayer\tExcit\tInhib\tSize (MB)\n");
	for (l=0; l<mp->nLayers; l++)
		fprintf(stderr,"%d\t%d\t%d\t%.2f\n",l,mp->vExcit[l],mp->vInhib[l],\
				(mp->vExcit[l]+mp->vInhib[l])*(sizeof(NEURON)+(mp->spkBuffer*(float)sizeof(tstep)))/MB);
	//fprintf(stderr,"%d\t%d\t%d\t%.2f\n",l,mp->vExcit[l],mp->vInhib[l],(mp->vExcit[l]+mp->vInhib[l])*(sizeof(NEURON)+(((l>0)?mp->spkBuffer:mp->inpSpkBuff)*(float)sizeof(tstep)))/MB);
	
	int Esyn = 0;
	float mem = 0.0;
	
	fprintf(stderr,"\nPresynaptic connection probabilites for Excitatory postsynaptic cells\n");
	fprintf(stderr,"Synapse:\tEfE \tElE \tIE \tE(syn) \tSize (MB)\n");
	// Assumes minimum delay model i.e. 1 spike bin per axon
	for (l=0; l<mp->nLayers; l++)
	{
		Esyn = (((l>0)?(mp->vExcit[l-1]*mp->pCnxEfE[l]):0)+(mp->vExcit[l]*mp->pCnxElE[l])+(mp->vInhib[l]*mp->pCnxIE[l]))*mp->vExcit[l];
		mem = (sizeof(AXON) + (sizeof(NEURON*)*3) + sizeof(tstep))*(float)Esyn/MB;
		fprintf(stderr,"Layer %d:\t%.3f\t%.3f\t%.3f\t%.2G\t%.2f\n",l,mp->pCnxEfE[l],mp->pCnxElE[l],mp->pCnxIE[l],(float)Esyn,mem);
	}
	
	fprintf(stderr,"\nPresynaptic connection probabilites for Inhibitory postsynaptic cells\n");
	fprintf(stderr,"Synapse:\tEI \tII \tE(syn) \tSize (MB)\n");
	for (l=0; l<mp->nLayers; l++)
	{
		Esyn = ((mp->vExcit[l]*mp->pCnxEI[l])+(mp->vInhib[l]*mp->pCnxII[l]))*mp->vInhib[l];
		mem = (sizeof(AXON) + (sizeof(NEURON*)*3) + sizeof(tstep))*(float)Esyn/MB;
		fprintf(stderr,"Layer %d:\t%.3f\t%.3f\t%.2G\t%.2f\n",l,mp->pCnxEI[l],mp->pCnxII[l],(float)Esyn,mem);
	}
	
	fprintf(stderr,"\nStimulus structures:\t%.2f", \
			(float)(sizeof(PARAMS) + ((mp->randStimOrder)?mp->loops*mp->nStimuli*sizeof(int):0) + \
			((mp->randTransOrder)?mp->loops*mp->nStimuli*mp->nTransPS*sizeof(int):0))/MB);
	fprintf(stderr,"\nTraining Stimuli:\t%.2f", (float)(mp->sInputs*mp->nStimuli*mp->nTransPS*sizeof(float))/MB);
	fprintf(stderr,"\nTesting Stimuli:\t%.2f\n\n", (float)(mp->sInputs*mp->nTestStimuli*mp->nTestTransPS*sizeof(float))/MB);
#endif
	
	printf("Ventral Visual Stream Spiking Neural Network Simulation starting...\n");
	
	/*************** Build & Initialize Network ****************/
	
	printf("\tNow building the network...");
	
	n_E = allocn(mp->nLayers, mp->vExcit, EXCIT); // Create 2D array of Excitatory neuron structures
	n_I = allocn(mp->nLayers, mp->vInhib, INHIB); // Create 2D array of Inhibatory neuron structures
	calcConnectivity(mp->probConnect);	// Calculate network connectivity

	printf("\tBuilding complete!\n");
	
	printf("\tNow initialising the network...");
	
	regime = Learning;			// 0: Testing (No STDP); 1: Training (STDP);
	init_network(regime);		// Initialise parameters
	
	/*************** Load stimuli ****************/
	
	printf("\tNetwork initialisation complete!\n");
	if (mp->useFilteredImages)
	{
		printf("\tNow loading the images...");
		
		stim = myalloc(sizeof(*stim));
		stim->trnImages = get_7D_farray(mp->nStimuli, mp->nTransPS, \
										mp->nScales, mp->nOrients, mp->nPhases, mp->nRows, mp->nCols, 0.0);
		
		error = loadImages(stim, mp); // if (!error)...
		printf("\t{S%d,T%d}", mp->nStimuli, mp->nTransPS);
		if (mp->newTestSet)
			printf(" Test: {S%d,T%d}", mp->nTestStimuli, mp->nTestTransPS);
		
		printf("\tImages now loaded!\n");
	}
	else
	{
		printf("\tNow creating the stimuli...");
		
		stim = myalloc(sizeof(*stim));
		gen_stimuli(mp->localRep, stim, mp);			// Generate Patterns
		
		stimuli_FP = fopen("trn_stimuli.dat", "w");
		for (p=0; p<mp->nStimuli; p++)
		{
			fprintf(stimuli_FP, "*** Stimulus %d ***\n", p+1);
			print_farray(stimuli_FP, (stim->trn_stimuli)[p], mp->nTransPS, mp->nExcit);
			fprintf(stimuli_FP, "\n");
		}
		fclose(stimuli_FP);
		
		stimuli_FP = fopen("tst_stimuli.dat", "w");
		for (p=0; p<mp->nStimuli; p++)
		{
			fprintf(stimuli_FP, "*** Stimulus %d ***\n", p+1);
			print_farray(stimuli_FP, stim->tst_stimuli[p], mp->nTransPS, mp->nExcit);
			fprintf(stimuli_FP, "\n");
		}
		fclose(stimuli_FP);
		
		printf("\tStimuli now saved!\n");
	}
	
	genShuffles(stim, mp);

	if (!mp->newTestSet)
	{
		mp->nTestStimuli = mp->nStimuli;
		mp->nTestTransPS = mp->nTransPS;
		stim->tstImages = stim->trnImages;
	}
	
	if (SIM.Xgrid) // Set up variables for estimating percentage completion
	{
		SIM.tally = 0;
		SIM.ptTS = (mp->pretrain) ? mp->nTestStimuli * mp->nTestTransPS * mp->transP_Test * ceil(1/DT): 0;
		SIM.trainTS = (mp->train) ? mp->loops * mp->nStimuli * mp->nTransPS * mp->transP_Train * ceil(1/DT): 0;
		SIM.testTS = mp->nTestStimuli * mp->nTestTransPS * mp->transP_Test * ceil(1/DT);
		SIM.totTS = SIM.ptTS + SIM.trainTS + SIM.testTS;
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
		free_3D_farray(stim->trn_stimuli, mp->nStimuli);//, mp->nTransPS);
		free_3D_farray(stim->tst_stimuli, mp->nStimuli);//, mp->nTransPS);
	}
	myfree(stim); // *** Free new image arrays too
	
	printf("\tMemory Deallocated!\n");

	return 0;
}


NEURON ** allocn (int nLays, int * vNeurons, NTYPE type)
{
	int l, n;
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
	
	
	/* Allocate spike time bins and initialise connectors*/
#pragma omp parallel default(shared) private(l,n)
	{
	for (l=0; l<nLays; l++)
	{
#pragma omp for
		for (n=0; n<vNeurons[l]; n++)
		{
			narray[l][n].spkbin = 0;
			//if (mp->useFilteredImages && type==EXCIT && l==0)
			//	narray[l][n].spikeTimes = myalloc(mp->inpSpkBuff * sizeof(narray[l][n].spikeTimes[0]));
			//else
			narray[l][n].spikeTimes = myalloc(mp->spkBuffer * sizeof(narray[l][n].spikeTimes[0]));
			narray[l][n].type = (type==EXCIT) ? EXCIT : INHIB;
			narray[l][n].nFAff_E = 0;
			narray[l][n].FAffs_E = NULL; // **
			narray[l][n].lm1presyn_E = NULL; // **
			narray[l][n].nLAff_E = 0;
			narray[l][n].LAffs_E = NULL; // **
			narray[l][n].lm0presyn_E = NULL; // **
			narray[l][n].nLAff_I = 0;
			narray[l][n].LAffs_I = NULL; // **
			narray[l][n].lm0presyn_I = NULL; // **
			narray[l][n].n = n;
			narray[l][n].l = l;
			//narray[l][n].row = 
			//narray[l][n].col = 
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
				free_2D_farray(narray[l][n].rec->cellV);//, mp->loops);
				if (l>0)
					free_2D_farray(narray[l][n].rec->cellD);//, mp->loops);
				if (l<mp->nWLayers)
				{
					free_3D_farray(narray[l][n].rec->SynC, mp->loops);//, narray[l][n].nFAff_E);
					free_3D_farray(narray[l][n].rec->SynG, mp->loops);//, narray[l][n].nFAff_E);
					free_3D_farray(narray[l][n].rec->SynDG, mp->loops);//, narray[l][n].nFAff_E);
				}
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
	
	if (!probConnect)
	{
		calc_connectivity();
		return;
	}
	
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
				for (s=0; s<n_E[l][n].nLAff_E; s++) // if (mp->pCnxElE[l] > EPS)
					fprintf(connections_FP, "%d\t", n_E[l][n].lm0presyn_E[s]->n);
				fprintf(connections_FP, "\n");
			}
			fclose(connections_FP);
		}
		
		for (l=0; l<mp->nLayers; l++)
		{
			slen = snprintf(filename, FNAMEBUFF, "L%daffNeuronsEI.dat", l);
			assert(slen < FNAMEBUFF);
			connections_FP = myfopen(filename, "w");
			for (n=0; n<mp->vInhib[l]; n++)
			{
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
			if (ran3(&idum) < pEfn) // Make presynaptic connection
			{
				n->lm1presyn_E[n->nFAff_E++] = &n_E[n->l-1][s];
				n_E[n->l-1][s].nFEff_E++;
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
	
	/*** Lateral El synapses ***/
	if (pEln > EPS)
	{
		buff = ceil(pEln * mp->vExcit[n->l]);
		n->lm0presyn_E = myalloc(buff * sizeof(NEURON *)); // Create array of NEURON pointers
		n->nLAff_E = 0;
		for (s=0; s<mp->vExcit[n->l]; s++)
		{
			if (ran3(&idum) < pEln) // Make presynaptic connection
			{
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
			if (ran3(&idum) < pIln)
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

void calc_connectivity() /* This function randomly connects the neurons together */
{
	int n, l;
	int slen = 0;	
	char filename[FNAMEBUFF];
	FILE * connections_FP;
	
	/* Given a particular postsynaptic neuron <E|I>, n, and its synapse, s, the presynaptic 
	 cell forming the synapse <E|I> is given by presyncnx_<E|I><E|I>[l][n][s]. */
	/* NB If connections from the previous and the same layer are required, can label each neuron from 0 
	 to (NEXCIT*NLAYERS)-1 (the neuron ID) with l=floor(NID/NEXCIT) and n=((NID+1)%NEXCIT)-1.	*/
	
	/*	For affNeurons_EfE[0][n][s] : layer 0 cells --> layer 1 cells (NWLAYERS)
		For affNeurons_ElE[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS)
		For affNeurons_IE[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS)
		For affNeurons_EI[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS)
		For affNeurons_II[0][n][s] : layer 0 cells --> layer 0 cells (NLAYERS) */
	
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
	
	/***** Can rands1_new be used in parallel? *****/
	
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
		
	/************** FILE OUTPUT **************/
	
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
	
	/************** END OF FILE OUTPUT **************/
	
	return; //void?
}


void wire_afferents(NEURON * n, int l, int * affNcnx_fE, int * affNcnx_lE, int * affNcnx_I)
{ // n_E, n_I, 
	int s, choice;
	
	/*** Ef synapses ***/
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
	
	/*** Lateral El synapses ***/
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
	
	/*** Lateral I Synapses ***/
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

void create_axons(NEURON * n)
{
	/*** Axonal delays ***/
	// Delays specified in seconds must be converted to timesteps
	tstep delay = 0;
	int nbins = 0;
	int span = 0;
	int s;
	switch (mp->axonDelay) 
	{
		case MinD:
			delay = 1;
			nbins = 1;
			break;
		case ConstD:
			delay = ceil(mp->d_const/DT);
			delay = (!delay) ? 1 : delay;
			nbins = ceil(mp->d_const/mp->refract);
			nbins = (!nbins) ? 1 : nbins;
			break;
		case UniformD:
			span = mp->d_max - mp->d_min;
			break;
		case GaussD:
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
				delay = round(((ran3(&idum)*span)+mp->d_min)/DT);
				delay = (!delay) ? 1 : delay;
				nbins = ceil((delay*DT)/mp->refract);
				break;
			case GaussD:
				delay = round((mp->d_mean + mp->d_sd*gasdev(&idum))/DT);
				delay = (!delay) ? 1 : delay;
				nbins = ceil((delay*DT)/mp->refract);
				break;
			default:
				break;
		}
		n->FEffs_E[s].delay = delay;
		n->FEffs_E[s].size = nbins;
		n->FEffs_E[s].queue = myalloc(nbins * sizeof(tstep));
		init_queue(&(n->FEffs_E[s]));
	}
	for (s=0; s<n->nLEff_E; s++)
	{
		switch (mp->axonDelay) 
		{
			case UniformD:
				delay = round(((ran3(&idum)*span)+mp->d_min)/DT);
				delay = (!delay) ? 1 : delay;
				nbins = ceil((delay*DT)/mp->refract);
				break;
			case GaussD:
				delay = round((mp->d_mean + mp->d_sd*gasdev(&idum))/DT);
				delay = (!delay) ? 1 : delay;
				nbins = ceil((delay*DT)/mp->refract);
				break;
			default:
				break;
		}
		n->LEffs_E[s].delay = delay;
		n->LEffs_E[s].size = nbins;
		n->LEffs_E[s].queue = myalloc(nbins * sizeof(tstep));
		init_queue(&(n->LEffs_E[s]));
	}
	for (s=0; s<n->nLEff_I; s++)
	{
		switch (mp->axonDelay) 
		{
			case UniformD:
				delay = round(((ran3(&idum)*span)+mp->d_min)/DT);
				delay = (!delay) ? 1 : delay;
				nbins = ceil((delay*DT)/mp->refract);
				break;
			case GaussD:
				delay = round((mp->d_mean + mp->d_sd*gasdev(&idum))/DT);
				delay = (!delay) ? 1 : delay;
				nbins = ceil((delay*DT)/mp->refract);
				break;
			default:
				break;
		}
		n->LEffs_I[s].delay = delay;
		n->LEffs_I[s].size = nbins;
		n->LEffs_I[s].queue = myalloc(nbins * sizeof(tstep));
		init_queue(&(n->LEffs_I[s]));
	}
	return;
}

void init_network(int regime)
{
	int s, n, l;
	int r; //, syn, t;
	int choice;
	
	if (regime == Learning) // Only initialise before training, not testing
	{
		if (mp->nRecordsPL)
		{
			/* Set up recording bins */
			printf("\n");
			for (l=0; l<mp->nLayers; l++)
			{
				choice=rands1_new(0,0,&iv,1);	// Initialise RNG
				for (r=0; r<mp->nRecordsPL; r++)
				{
					// Randomly set NRECORDS flags
					choice=rands1_new(0,mp->vExcit[l]-1,&iv,1);
					n_E[l][choice].rec_flag = true;
					printf("\tLayer %d, Record %d Assigned nID: %d\n",l,r,choice);
				}
			}
			/* Create recording structures */
			for (l=0; l<mp->nLayers; l++)
				for (n=0; n<mp->vExcit[l]; n++)
					if (n_E[l][n].rec_flag)
					{
						n_E[l][n].rec = myalloc(sizeof(RECORD));
						n_E[l][n].rec->bin = 0;
						n_E[l][n].rec->cellV = get_2D_farray(mp->loops, mp->TotalMS, 0.0);
						n_E[l][n].rec->cellD = get_2D_farray(mp->loops, mp->TotalMS, 0.0);
						if (l>0)
						{
							n_E[l][n].rec->SynC = get_3D_farray(mp->loops, n_E[l][n].nFAff_E, mp->TotalMS, 0.0);
							n_E[l][n].rec->SynG = get_3D_farray(mp->loops, n_E[l][n].nFAff_E, mp->TotalMS, 0.0);
							n_E[l][n].rec->SynDG = get_3D_farray(mp->loops, n_E[l][n].nFAff_E, mp->TotalMS, 0.0);
						}
					}
		}
		
		/* E_ Synaptic weights */
		for (l=0; l<mp->nWLayers; l++)
			for (n=0; n<mp->vExcit[l]; n++)
				for (s=0; s<n_E[l][n].nFEff_E; s++)
					n_E[l][n].FEffs_E[s].delta_g = n_E[l][n].FEffs_E[s].delta_g_tm1 = ran3(&idum);
		
		for (l=0; l<mp->nLayers; l++)
			for (n=0; n<mp->vExcit[l]; n++)
				for (s=0; s<n_E[l][n].nLEff_E; s++)
					n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = ran3(&idum);
		
		for (l=0; l<mp->nLayers; l++)
			for (n=0; n<mp->vExcit[l]; n++)
				for (s=0; s<n_E[l][n].nLEff_I; s++)
					n_E[l][n].LEffs_I[s].delta_g = n_E[l][n].LEffs_I[s].delta_g_tm1 = mp->Dg_EI;
		
		/* I_ Synaptic weights */
		for (l=0; l<mp->nLayers; l++)
			for (n=0; n<mp->vInhib[l]; n++)
				for (s=0; s<n_I[l][n].nLEff_E; s++)
					n_I[l][n].LEffs_E[s].delta_g = n_I[l][n].LEffs_E[s].delta_g_tm1 = mp->Dg_IE;
		
		for (l=0; l<mp->nLayers; l++)
			for (n=0; n<mp->vInhib[l]; n++)
				for (s=0; s<n_I[l][n].nLEff_I; s++)
					n_I[l][n].LEffs_I[s].delta_g = n_I[l][n].LEffs_I[s].delta_g_tm1 = mp->Dg_II;
		
		/*switch (mp->noise)
		{
			case 0: // No noise // Already initialised to 0 when declared
				break;
			case 1: // Uniform interval
				break;
			case 2: // Gaussian noise
				break;
		} // End of NOISE switch	*/	
	} // End of if (regime==1) clause
	
	
	/********** Training and Testing **********/
	// Executed for all types of initialization (except Settle where shown)
	
	/* Membrane potentials, conductances, learning parameters and spike time arrays */
	for (l=0; l<mp->nLayers; l++)
	{
		for (n=0; n<mp->vExcit[l]; n++)
		{
			n_E[l][n].V = n_E[l][n].V_tm1 = mp->VrestE;
			n_E[l][n].D = n_E[l][n].D_tm1 = 0.0;
			for (s=0; s<n_E[l][n].nFEff_E; s++) // nFEff_E should be 0 for the last layer
			{
				n_E[l][n].FEffs_E[s].C = n_E[l][n].FEffs_E[s].C_tm1 = 0.0;
				n_E[l][n].FEffs_E[s].g = n_E[l][n].FEffs_E[s].g_tm1 = 0.0;
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
		for (n=0; n<mp->vInhib[l]; n++)
		{
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

	/* Recording structures */
	if (mp->nRecordsPL && regime != Settle)
		for (l=0; l<mp->nLayers; l++)
			for (n=0; n<mp->vExcit[l]; n++)
				if (n_E[l][n].rec_flag)
					n_E[l][n].rec->bin = 0;
	
	switch (mp->noise) // Re-initialise
	{
		case 0: // No noise - no need to reinit in this case
			break;
		case 1: // Uniform interval
			break;
		case 2: // Gaussian noise
			break;
	}
	
	// somLeak_E[][]

	return;
}

int loadImages(STIMULI * stim, PARAMS * mp)
{
	FILE * iList;
	FILEPARTS * fp;
	char * str, buff[BUFFER], filename[FNAMEBUFF]; //dir[DIRBUFF],
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
	
	iList = myfopen(mp->imageList, "r");
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
	if (!err)
		return 1;
	else
		return 0;
}

int loadGaborOutput(const char * filename, float * array, int size)
{
	FILE * gbo;
	int e = 0; 
	int err = 0;
	gbo = myfopen(filename, "r");
	fread(array, sizeof(float), size, gbo);
	//assert(feof(gbo));
	fclose(gbo);
	for (e=0; e<size; e++) 	// Scale for current injected too
	{
		assert(-1 <= array[e] && array[e] <= 1);
		array[e] = (array[e] * mp->currentSpread) + mp->current;
		if (array[e] < 0) // Check that negative currents are not injected
			array[e] = 0;
	}

	if (!err)
		return 1;
	else
		return 0;
}

void gen_stimuli(bool rep, STIMULI * stim, PARAMS * mp)
{
	// Place array arguements in a patterns structure
	int trans;
	int n, p;
	int choice;
	int block = 0;

	/* Generate the training stimuli to present to the network */
	// Distributed or orthogonal will depend on where the balls are replaced i.e. rng reinitialised...
	// Patterns could be generated in Matlab and loaded from dat files...
	// Or read in list of pairs from another array
	
	stim->trn_stimuli = get_3D_farray(mp->nStimuli, mp->nTransPS, mp->nExcit, 0.0);
	stim->tst_stimuli = get_3D_farray(mp->nStimuli, mp->nTransPS, mp->nExcit, 0.0);
	
	//int nFiringNeurons = floor(mp->a * mp->nExcit);
	
	switch (rep) 
	{
		case 0: /* Distributed Patterns */ // *** FINISH THIS!!! ***
			trans=0;
			for (p=0; p<mp->nStimuli; p++)
			{
				i=rands1_new(0,0,&iv,1);	// Initialise random number generator - see p11
				for (n=0; n<mp->nFiringNeurons; n++)
				{
					// Sample from possible pre-synaptic neurons without replacement
					choice=rands1_new(0,mp->nExcit-1,&iv,1); // see p11
					stim->trn_stimuli[p][trans][choice] = stim->tst_stimuli[p][trans][choice] = 1.00 * mp->current;
				}
			}			
			break;
			
		case 1: /* Local Patterns */
			block = floor(mp->nFiringNeurons + ((mp->nTransPS - 1) * mp->shift));
			// This assumes stimuli are a contiguous block of 1's
			for (p=0; p<mp->nStimuli; p++)
				for (trans=0; trans<mp->nTransPS; trans++)
					for (n=(p*block)+(trans*mp->shift); n<(mp->nFiringNeurons+(p*block)+(trans*mp->shift)); n++)
						stim->trn_stimuli[p][trans][n] = stim->tst_stimuli[p][trans][n] = 1.00 * mp->current;
			break;
	}
	
	/*	***111111111222222222222333333333333
		111***111111222222222222333333333333
		...
		------------***222222222------------
		------------222***222222------------
		...
		------------------------333333333***
	 */ // Extend this to 2D - cf MSc code
	
#if DEBUG>1 // Level 2
	fprintf(stderr, "\nPrinting generated training stimuli...\n");
	for (p=0; p<mp->nStimuli; p++)
		for (trans=0; trans<mp->nTransPS; trans++)
		{
			fprintf(stderr, "S%dT%d: ",p+1,trans+1);
			print_frow(stderr, stim->trn_stimuli[p][trans], mp->nExcit);
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
	int i = 0;
	int count = 0;
	
	if (mp->randStimOrder)
	{
		stim->stimShuffle = get_2D_iarray(mp->loops, mp->nStimuli, 0);
		for (loop=0; loop<mp->loops; loop++)
		{
			i=rands1_new(0,0,&iv,1);	// Initialise RNG (clear internal array of numbers drawn)
			for (p=0; p<mp->nStimuli; p++)
				stim->stimShuffle[loop][p] = rands1_new(0,mp->nStimuli-1,&iv,1);
		}
		count++;
	}
	
	if (mp->randTransOrder)
	{
		stim->transShuffle = get_3D_iarray(mp->loops, mp->nStimuli, mp->nTransPS, 0);
		for (loop=0; loop<mp->loops; loop++)
		{
			for (p=0; p<mp->nStimuli; p++)
			{
				i=rands1_new(0,0,&iv,1);	// Initialise RNG (clear internal array of numbers drawn)
				for (trans=0; trans<mp->nTransPS; trans++)
					stim->transShuffle[loop][p][trans] = rands1_new(0,mp->nTransPS-1,&iv,1);
			}
		}
		count++;
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
				*input = ****(stim->tstImages[pat][trans]);
				//memcpy(input, ****(stim->tstImages[pat][trans]), mp->sInputs*sizeof(input[0])); // Use pointer not memcpy
				break;
				
			case 1:
				pat = (mp->randStimOrder) ? stim->stimShuffle[loop][pat] : pat;
				trans = (mp->randTransOrder) ? stim->transShuffle[loop][pat][trans] : trans;
				patArray = (mp->newTestSet) ? stim->tstImages[pat][trans] : stim->trnImages[pat][trans];
				*input = ****patArray;
				//memcpy(input, ****patArray, mp->sInputs*sizeof(input[0]));
				break;
		}
	}
	else
	{
		switch (regime)
		{
			case 0: // Testing stimuli
				*input = stim->tst_stimuli[pat][trans]; //Test me!
				//memcpy(input, stim->tst_stimuli[pat][trans], mp->nExcit*sizeof(input[0])); // Copy row of tst_stimuli to input
				break;
				
			case 1: // Training stimuli
				pat = (mp->randStimOrder) ? stim->stimShuffle[loop][pat] : pat;
				trans = (mp->randTransOrder) ? stim->transShuffle[loop][pat][trans] : trans;
				*input = stim->trn_stimuli[pat][trans];
				//memcpy(input, stim->trn_stimuli[pat][trans], mp->nExcit*sizeof(input[0]));
				break;
		}
	}

	// Move print statements outside function?
	if ( trans==0 || (mp->interleaveTrans && regime) )
		printf("\t\tPresenting stimulus %d/%d...\n", (pat+1), mp->nStimuli);
	printf("\t\t\tTransform %d/%d...\n", (trans+1), mp->nTransPS);

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
	REGIMETYPE regime = NoLearning;
	int slen = 0;
	char phaseString[FNAMEBUFF];
	char filename[FNAMEBUFF];
	FILE * excitOutput;
	FILE * inhibOutput;
	FILE * weightOutput;

	FILE * rCellVout;
	FILE * rDout;
	FILE * rCout;
	FILE * rGout;
	FILE * rWeightOut;
	
	double percentage = 0.0;
	
	float * input = NULL; //myalloc(mp->sInputs * sizeof(*input));
	
	switch (sPhase) 
	{
		case PreTraining:
			strncpy(phaseString, "Pretraining", FNAMEBUFF);
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
			if (mp->trainPause && sPhase==Training && !mp->interleaveTrans)
				init_network(Settle);
			for (*i=0; *i<iCount; (*i)++)
			{
				if (sPhase != Training) //Testing only
					init_network(Settle);
				calc_input(loop, p, tr, stim, &input, regime);
				t_start = round((*i + (*o * iCount)) * transP * ceil(1/DT));
				t_end = t_start + round(transP * ceil(1/DT));
#if DEBUG > 1 // Level 2
				fprintf(stderr, "\nUpdating network from timestep %d to %d.\n", t_start, t_end-1); //%lld
#endif
				
#if DEBUG > 3
				print_frow(stderr, input, mp->sInputs);
#endif
				updateNetwork(t_start, t_end, loop, input, regime);
				/*if (mp->saveInputSpikes)
				{
					// Print spikes to file (see below) ...
					slen = snprintf(filename, FNAMEBUFF, "E%dL0ExcitSpikes.dat", loop);
				}
				for (n=0; n<mp->sInputs; n++)
				{
					n_E[0][n].spkbin = 0;
					memset(n_E[0][n].spikeTimes, 0, mp->inpSpkBuff*sizeof(n_E[0][n].spikeTimes[0])); //[0]?
				}*/
				if (SIM.Xgrid)
				{
					SIM.tally += transP * ceil(1/DT);
					percentage = (100.0*SIM.tally)/SIM.totTS;
					/*if (sPhase == PreTraining)
						percentage = (double) 100*(t_end + 1)/SIM.totTS;
					else if (sPhase == Training)
						percentage = (double) 100*((mp->pretrain)?SIM.ptTS:0.0 + (SIM.trainTS*loop/(double)mp->loops) + t_end + 1)/SIM.totTS;
					else if (sPhase == Testing)
						percentage = (double) 100*((mp->pretrain)?SIM.ptTS:0.0 + (mp->train)?SIM.trainTS:0.0 + t_end + 1)/SIM.totTS;*/
					
					printf("<xgrid>{control = statusUpdate; percentDone = %.0lf; }</xgrid>\n", percentage);
					fflush(stdout);
				}
			}
		}
		printf("\t%s complete!\n",phaseString);
		printf("\tSaving results...");
		// Output results to dat files
		
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
				print_irow(excitOutput, (int*) n_E[l][n].spikeTimes, n_E[l][n].spkbin+1); // spikeTimes[0] = -BIG?
			fclose(excitOutput);
		}
		// fprintf(pt_output_spikes_ptr, "L%dN%d\t", NLAYERS-1, n);
		
		if (sPhase == Testing) // Save Inhibitory spikes for all layers
		{
			for (l=0; l<mp->nLayers; l++)
			{
				slen = snprintf(filename, FNAMEBUFF, "L%dInhibSpikes.dat", l);
				assert(slen < FNAMEBUFF);
				inhibOutput = myfopen(filename, "w"); // Make either HR "w" or Binary "wb"
				for (n=0; n<mp->vInhib[l]; n++)
					print_irow(inhibOutput, (int*) n_I[l][n].spikeTimes, n_I[l][n].spkbin+1);
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
		}
		
		if (sPhase==Testing && mp->nRecordsPL) // Move outside loops loop?
		{
			for (l=0; l<mp->nLayers; l++)
			{
				for (n=0; n<mp->vExcit[l]; n++)
				{
					if (n_E[l][n].rec_flag)
					{
						slen = snprintf(filename, FNAMEBUFF, "rL%dN%dcellV.dat", l, n);
						assert(slen < FNAMEBUFF);
						rCellVout = myfopen(filename, "w"); //r_cellV_ptr = fopen(filename, "wb");
						print_farray(rCellVout, n_E[l][n].rec->cellV, mp->loops, mp->TotalMS); //fwrite(n_E[l][n].rec->cellV, sizeof(float), mp->loops*mp->TotalMS, r_cellV_ptr); // See nifty trick #1
						fclose(rCellVout);
						
						if (l>0)
						{
							slen = snprintf(filename, FNAMEBUFF, "rL%dN%dD.dat", l, n);
							assert(slen < FNAMEBUFF);
							rDout = myfopen(filename, "w"); //r_D_ptr = fopen(filename, "wb");
							print_farray(rDout, n_E[l][n].rec->cellD, mp->loops, mp->TotalMS); //fwrite(n_E[l][n].rec->cellD, sizeof(float), mp->loops*mp->TotalMS, r_D_ptr); // See nifty trick #1 //sizeof(n_E[l][n].rec->cellD)/sizeof(float)
							fclose(rDout);
							
							// The presynaptic cell's values are attached to each record
							slen = snprintf(filename, FNAMEBUFF, "rL%dN%dC.dat", l, n);
							assert(slen < FNAMEBUFF);
							rCout = myfopen(filename, "w"); //r_C_ptr = fopen(filename, "wb");
							for (loop=0; loop<mp->loops; loop++)
								print_farray(rCout, n_E[l][n].rec->SynC[loop], n_E[l][n].nFAff_E, mp->TotalMS);	//fwrite(n_E[l][n].rec->SynC, sizeof(float), mp->loops*n_E[l][n].nFAff_E*mp->TotalMS, r_C_ptr); //fwrite(n_E[l][n].rec->SynC, sizeof(float), sizeof(n_E[l][n].rec->SynC)/sizeof(float), r_C_ptr); // See nifty trick #1
							fclose(rCout);
							
							slen = snprintf(filename, FNAMEBUFF, "rL%dN%dg.dat", l, n);
							assert(slen < FNAMEBUFF);
							rGout = myfopen(filename, "w");
							for (loop=0; loop<mp->loops; loop++)
								print_farray(rGout, n_E[l][n].rec->SynG[loop], n_E[l][n].nFAff_E, mp->TotalMS);
							fclose(rGout);
							
							slen = snprintf(filename, FNAMEBUFF, "rL%dN%ddg.dat", l, n);
							assert(slen < FNAMEBUFF);
							rWeightOut = myfopen(filename, "w");
							for (loop=0; loop<mp->loops; loop++)
								print_farray(rWeightOut, n_E[l][n].rec->SynDG[loop], n_E[l][n].nFAff_E, mp->TotalMS);
							fclose(rWeightOut);
						}
					}
				}
			}
		} // End of state variable records
		printf("\tResults saved!\n");
	} // End of loops
	//myfree(input);
	return;
}

void updateNetwork(tstep t_start, tstep t_end, int loop, float input[], int regime)
{
	int n, syn, l; //, wl;
	float decay_rate, gLeak, Vrest, Thresh, Vhyper;
	float decay_E, decay_I;
	tstep t=0;
	int bin = 0;

#pragma omp parallel default(shared) private(t,l,n,syn,bin,decay_rate,gLeak,Vrest,Thresh,Vhyper,decay_E,decay_I)
	{
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
		decay_rate = DT/mp->capE;
		gLeak = mp->gLeakE;
		Vrest = mp->VrestE;
		Thresh = mp->ThreshE;
		Vhyper = mp->VhyperE;
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait // schedule(static) // Needed for thread-safe noise
			for (n=0; n<mp->vExcit[l]; n++)
				update_V(t, decay_rate, gLeak, Vrest, Thresh, Vhyper, (l==0)?input[n]:0.0, &n_E[l][n]);
		}
		
		/* Update Inhibitory cell potentials */
		decay_rate = DT/mp->capI;
		gLeak = mp->gLeakI;
		Vrest = mp->VrestI;
		Thresh = mp->ThreshI;
		Vhyper = mp->VhyperI;
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait // It is ok to have a nowait here with a minimum conduction delay of 1 tstep
			for (n=0; n<mp->vInhib[l]; n++) // Larger input layer separation as for Excit cells?
				update_V(t, decay_rate, gLeak, Vrest, Thresh, Vhyper, 0.0, &n_I[l][n]);
		}
//#pragma omp barrier
		
		/* Update presynaptic conductances */
		decay_E = (DT/mp->tauEE);
		decay_I = (DT/mp->tauIE);
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait //schedule(runtime) - experiment!!
			for (n=0; n<mp->vExcit[l]; n++)
				update_g(&n_E[l][n], decay_E, decay_I, t);
		}
		
		decay_E = (DT/mp->tauEI);
		decay_I = (DT/mp->tauII);
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for
			for (n=0; n<mp->vInhib[l]; n++)
				update_g(&n_I[l][n], decay_E, decay_I, t);
		}
//#pragma omp barrier
		
		if (regime==Learning) // Learning
		{
			/* Update synaptic weights */ // Move after C & D and use instantaeous values? N.B. Redo nowait clauses
			for (l=1; l<mp->nLayers; l++) // Only nLayers-1 of weight layers
			{
#pragma omp for nowait
				for (n=0; n<mp->vExcit[l]; n++) // Make parallel and n private
					update_weights(&n_E[l][n], t);
			}
			
			// -->|| Update C for current neuron's outgoing synapses (skip last layer)
			for (l=1; l<mp->nLayers; l++) // Only nLayers-1 of synapses
			{
#pragma omp for nowait
				for (n=0; n<mp->vExcit[l]; n++)
					update_C(&n_E[l][n], t);
			}
			
			// ||--> Update D for current neuron's incoming synapses (skip first layer)
			for (l=1; l<mp->nLayers; l++) // Only nLayers-1 of weight layers
			{
#pragma omp for
				for (n=0; n<mp->vExcit[l]; n++)
					update_D(&n_E[l][n], t);
				//#pragma omp barrier
			}
			//#pragma omp barrier
		}
		
		/* Copy solution variables to _tm1 counterparts and reset spike flags / axons */
		for (l=0; l<mp->nLayers; l++)
		{
#pragma omp for nowait// private(l, n, s)
			for (n=0; n<mp->vExcit[l]; n++)
			{
				//n_E[l][n].fired = 0;
				n_E[l][n].V_tm1 = n_E[l][n].V;
				n_E[l][n].D_tm1 = n_E[l][n].D;
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
					//n_E[l][n].LEffs_E[s].delta_g_tm1 = n_E[l][n].LEffs_E[s].delta_g;
					//n_E[l][n].LEffs_E[s].C_tm1 = n_E[l][n].LEffs_E[s].C;
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
#pragma omp for
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
//#pragma omp barrier
		}
		
		if (mp->nRecordsPL && regime && (t % mp->TSperMS == 0)) // && loop==RLOOP) 		// Save neuron states
		{
			//int bin = 0;
			// Dump records to file after each loop?
			for (l=0; l<mp->nLayers; l++)
			{
#pragma omp for// private(bin, syn, l, n)
				for (n=0; n<mp->vExcit[l]; n++) //for (r=0; r<NRECORDS_PL; r++)
					if (n_E[l][n].rec_flag)
					{
						bin = n_E[l][n].rec->bin++;
						n_E[l][n].rec->cellV[loop][bin] = n_E[l][n].V;
						n_E[l][n].rec->cellD[loop][bin] = n_E[l][n].D;
						for (syn=0; syn<n_E[l][n].nFAff_E; syn++) //for (syn=0; syn<NSYN_PN; syn++)
						{
							n_E[l][n].rec->SynC[loop][syn][bin] = n_E[l][n].FAffs_E[syn]->C;
							n_E[l][n].rec->SynG[loop][syn][bin] = n_E[l][n].FAffs_E[syn]->g;
							n_E[l][n].rec->SynDG[loop][syn][bin] = n_E[l][n].FAffs_E[syn]->delta_g;
						}
					}
			}
//#pragma omp barrier
		}		
	}
	}
	return; // void;
}


void update_network(tstep t, int loop, float input[], int regime)
{
	int n, syn, s, l; //, wl;
	int nNeurons;
	float decay_rate, gLeak, Vrest, Thresh, Vhyper;
	float decay_E, decay_I;

#if DEBUG>2
#pragma omp single nowait
	{
	fprintf(stderr,"\rTimestep: %d",t);
	fflush(stderr);
	}
#endif
	
	// *** Make the update functions inline? ***
/*	if (mp->useFilteredImages)
	{
		/* Update Input neurons cell potentials */
/*		decay_rate = DT/mp->capE;
		gLeak = mp->gLeakE;
		Vrest = mp->VrestE;
		Thresh = mp->ThreshE;
		Vhyper = mp->VhyperE;
#pragma omp for private(n,decay_rate,gLeak,Vrest,Thresh,Vhyper)
		for (n=0; n<mp->sInputs; n++) // nExcit or sInputs
			update_V(t, decay_rate, gLeak, Vrest, Thresh, Vhyper, input[n], &n_E[0][n]);
//#pragma omp barrier
	}*/
	
	/* Update Excitatory cell potentials */
	decay_rate = DT/mp->capE;
	gLeak = mp->gLeakE;
	Vrest = mp->VrestE;
	Thresh = mp->ThreshE;
	Vhyper = mp->VhyperE;
	//for (l=(mp->useFilteredImages)?1:0; l<mp->nLayers; l++)
#pragma omp parallel default(shared) private(l,n)//,nNeurons)
	for (l=0; l<mp->nLayers; l++)
	{
		//nNeurons = mp->vExcit[l];
//#pragma omp for nowait private(n,l,nNeurons,decay_rate,gLeak,Vrest,Thresh,Vhyper)
#pragma omp for
		for (n=0; n<mp->vExcit[l]; n++)
			update_V(t, decay_rate, gLeak, Vrest, Thresh, Vhyper, (l==0)?input[n]:0.0, &n_E[l][n]);
	}
	
	/* Update Inhibitory cell potentials */
	decay_rate = DT/mp->capI;
	gLeak = mp->gLeakI;
	Vrest = mp->VrestI;
	Thresh = mp->ThreshI;
	Vhyper = mp->VhyperI;
	for (l=0; l<mp->nLayers; l++)
	{
		nNeurons = mp->vInhib[l];
//#pragma omp for private(n,l,nNeurons,decay_rate,gLeak,Vrest,Thresh,Vhyper)
		for (n=0; n<nNeurons; n++) // Larger input layer separation as for Excit cells?
			update_V(t, decay_rate, gLeak, Vrest, Thresh, Vhyper, 0.0, &n_I[l][n]);
	}
//#pragma omp barrier

	/* Update presynaptic conductances */
	decay_E = (DT/mp->tauEE);
	decay_I = (DT/mp->tauIE);
	for (l=0; l<mp->nLayers; l++)
	{
		nNeurons = mp->vExcit[l];
//#pragma omp for nowait private(n,l,nNeurons,decay_E,decay_I)
		for (n=0; n<nNeurons; n++)
			update_g(&n_E[l][n], decay_E, decay_I, t);
	}
		
	decay_E = (DT/mp->tauEI);
	decay_I = (DT/mp->tauII);
	for (l=0; l<mp->nLayers; l++)
	{
		nNeurons = mp->vInhib[l];
//#pragma omp for private(n,l,nNeurons,decay_E,decay_I)
		for (n=0; n<nNeurons; n++)
			update_g(&n_I[l][n], decay_E, decay_I, t);
	}
//#pragma omp barrier
	
	if (regime==Learning) // Learning
	{
		/* Update synaptic weights */ // Move after C & D and use instantaeous values? N.B. Redo nowait clauses
		for (l=1; l<mp->nLayers; l++) // Only nLayers-1 of weight layers
		{
			nNeurons = mp->vExcit[l];
//#pragma omp for nowait private(n,l,nNeurons) // shared(n_E[][])
			for (n=0; n<nNeurons; n++) // Make parallel and n private
				update_weights(&n_E[l][n], t);
		}
		
		// -->|| Update C for current neuron's outgoing synapses (skip last layer)
		for (l=1; l<mp->nLayers; l++) // Only nLayers-1 of synapses
		{
			nNeurons = mp->vExcit[l];
//#pragma omp for nowait private(n,l,nNeurons) // shared(n_E[][])
			for (n=0; n<nNeurons; n++)
				update_C(&n_E[l][n], t);
		}
			
		// ||--> Update D for current neuron's incoming synapses (skip first layer)
		for (l=1; l<mp->nLayers; l++) // Only nLayers-1 of weight layers
		{
			nNeurons = mp->vExcit[l];
//#pragma omp for private(n,l,nNeurons) // shared(n_E[][]) // Last loop must be synchronised i.e. !(nowait)
			for (n=0; n<nNeurons; n++)
				update_D(&n_E[l][n], t);
//#pragma omp barrier
		}
	}
		
	/* Copy solution variables to _tm1 counterparts and reset spike flags / axons */
	for (l=0; l<mp->nLayers; l++)
	{
//#pragma omp for nowait private(l, n, s)
		for (n=0; n<mp->vExcit[l]; n++)
		{
			//n_E[l][n].fired = 0;
			n_E[l][n].V_tm1 = n_E[l][n].V;
			n_E[l][n].D_tm1 = n_E[l][n].D;
			for (s=0; s<n_E[l][n].nFEff_E; s++)
			{
				if (t == next_spike(&n_E[l][n].FEffs_E[s]))
					dequeue(&n_E[l][n].FEffs_E[s]);
				n_E[l][n].FEffs_E[s].g_tm1 = n_E[l][n].FEffs_E[s].g;
				n_E[l][n].FEffs_E[s].delta_g_tm1 = n_E[l][n].FEffs_E[s].delta_g;
				n_E[l][n].FEffs_E[s].C_tm1 = n_E[l][n].FEffs_E[s].C;
			}
			for (s=0; s<n_E[l][n].nLEff_E; s++)
			{
				if (t == next_spike(&n_E[l][n].LEffs_E[s]))
					dequeue(&n_E[l][n].LEffs_E[s]);
				n_E[l][n].LEffs_E[s].g_tm1 = n_E[l][n].LEffs_E[s].g;
				//n_E[l][n].LEffs_E[s].delta_g_tm1 = n_E[l][n].LEffs_E[s].delta_g;
				//n_E[l][n].LEffs_E[s].C_tm1 = n_E[l][n].LEffs_E[s].C;
			}
			for (s=0; s<n_E[l][n].nLEff_I; s++)
			{
				if (t == next_spike(&n_E[l][n].LEffs_I[s]))
					dequeue(&n_E[l][n].LEffs_I[s]);
				n_E[l][n].LEffs_I[s].g_tm1 = n_E[l][n].LEffs_I[s].g;
			}
		}
	}
	
	for (l=0; l<mp->nLayers; l++)
	{
//#pragma omp for private(l, n, s)
		for (n=0; n<mp->vInhib[l]; n++)
		{
			//n_I[l][n].fired = 0;
			n_I[l][n].V_tm1 = n_I[l][n].V;
			for (s=0; s<n_I[l][n].nLEff_E; s++)
			{
				if (t == next_spike(&n_I[l][n].LEffs_E[s]))
					dequeue(&n_I[l][n].LEffs_E[s]);
				n_I[l][n].LEffs_E[s].g_tm1 = n_I[l][n].LEffs_E[s].g;
			}
			for (s=0; s<n_I[l][n].nLEff_I; s++)
			{
				if (t == next_spike(&n_I[l][n].LEffs_I[s]))
					dequeue(&n_I[l][n].LEffs_I[s]);
				n_I[l][n].LEffs_I[s].g_tm1 = n_I[l][n].LEffs_I[s].g;
			}
		}
	}
//#pragma omp barrier
	
	if (mp->nRecordsPL && regime && (t % mp->TSperMS == 0)) // && loop==RLOOP) 		// Save neuron states
	{
		int bin = 0;
		// Dump records to file after each loop?
		for (l=0; l<mp->nLayers; l++)
		{
//#pragma omp for private(bin, syn, l, n)
			for (n=0; n<mp->vExcit[l]; n++) //for (r=0; r<NRECORDS_PL; r++)
				if (n_E[l][n].rec_flag)
				{
					bin = n_E[l][n].rec->bin++;
					n_E[l][n].rec->cellV[loop][bin] = n_E[l][n].V;
					n_E[l][n].rec->cellD[loop][bin] = n_E[l][n].D;
					for (syn=0; syn<n_E[l][n].nFAff_E; syn++) //for (syn=0; syn<NSYN_PN; syn++)
					{
						n_E[l][n].rec->SynC[loop][syn][bin] = n_E[l][n].FAffs_E[syn]->C;
						n_E[l][n].rec->SynG[loop][syn][bin] = n_E[l][n].FAffs_E[syn]->g;
						n_E[l][n].rec->SynDG[loop][syn][bin] = n_E[l][n].FAffs_E[syn]->delta_g;
					}
				}
		}
//#pragma omp barrier
	}

	return; // void;
}

void update_V(tstep t, float decay_rate, float gLeak, float Vrest, float Thresh, float Vhyper, float inj, NEURON * n)
{	/* This function updates a neuron's cell potential and applies to all layers */
	int syn;
	float tot_g_E = 0.0;
	float tot_g_I = 0.0;
	if (t > n->nextUpdate) //if ((((t - n->lastSpike) * DT) > mp->refract) || n->spkbin==0) // Calculated at spike time
	{
		/* Feed-forward connections (n->type==EXCIT && l>0) */
		for (syn=0; syn<n->nFAff_E; syn++)
			tot_g_E += n->FAffs_E[syn]->g_tm1;
		
		/* Lateral connections */
		for (syn=0; syn<n->nLAff_E; syn++)
			tot_g_E += n->LAffs_E[syn]->g_tm1;
		
		for (syn=0; syn<n->nLAff_I; syn++)
			tot_g_I += n->LAffs_I[syn]->g_tm1;
		
		n->V += decay_rate * (gLeak * (Vrest - n->V_tm1) \
							  + (tot_g_E * (mp->VrevE - n->V_tm1)) \
							  + (tot_g_I * (mp->VrevI - n->V_tm1)) \
							  + inj);
		
		if (mp->noise)
		{
			int th = 0;
			float sigma = 0.0;
#ifdef _OPENMP
			th = omp_get_thread_num();
#endif
			sigma = (n->type == EXCIT) ? mp->SigmaE : mp->SigmaI;
			n->V += (gsl_ran_gaussian(states[th], sqrt(DT)) * sigma);
		}
		
		n->V = ((n->V < mp->VrevI) ? mp->VrevI : n->V); // Neurons can not be more -ve than Inhib reversal potential
		
		if (n->V >= Thresh) // n->V_tm1 ??? - gives one timestep to depolarise (peak)
		{
#if DEBUG>3
			fprintf(stderr," L%dN%d%c: t=%d",n->l,n->n,(n->type==EXCIT)?'E':'I',t);
#endif
			n->V = Vhyper;
			n->lastSpike = n->spikeTimes[++(n->spkbin)] = t;
			n->nextUpdate = t + ceil(mp->refract/DT);
			// Enqueue spike for all post-synaptic neurons
			for (syn=0; syn<n->nFEff_E; syn++)
				enqueue(&n->FEffs_E[syn], t);
						
			for (syn=0; syn<n->nLEff_E; syn++)
				enqueue(&n->LEffs_E[syn], t);
						
			for (syn=0; syn<n->nLEff_I; syn++)
				enqueue(&n->LEffs_I[syn], t);
		} // End of spike
	} // End of REFRACT check 
	return;
}


void update_g(NEURON * n, float decay_E, float decay_I, tstep t)
{
	/* This function updates a neuron's afferent (pre-synaptic) conductances and applies to l>0 */
	float impulse;
	int syn;
	/* Update E synapses */
	for (syn=0; syn<n->nFAff_E; syn++) // Make private
	{
		impulse = (next_spike(n->FAffs_E[syn]) == t) ? n->FAffs_E[syn]->delta_g_tm1 * mp->gMax : 0.0;
		n->FAffs_E[syn]->g += (impulse - (decay_E * n->FAffs_E[syn]->g_tm1));
	}
	
	for (syn=0; syn<n->nLAff_E; syn++)
	{
		impulse = (next_spike(n->LAffs_E[syn]) == t) ? n->LAffs_E[syn]->delta_g_tm1 * mp->gMax : 0.0;
		n->LAffs_E[syn]->g += (impulse - (decay_E * n->LAffs_E[syn]->g_tm1));
	}
	
	/* Update I synapses */
	for (syn=0; syn<n->nLAff_I; syn++) // Already private
	{
		impulse = (next_spike(n->LAffs_I[syn]) == t) ? n->LAffs_I[syn]->delta_g_tm1 * mp->gMax : 0.0;
		n->LAffs_I[syn]->g += (impulse - (decay_I * n->LAffs_I[syn]->g_tm1));
		// Shunting inhibition ?
	}
	return;
}

void update_weights(NEURON * n, tstep t)
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
	return;
}
	
void update_C(NEURON * n, tstep t)
{
	// -->|| Update C for current neuron's outgoing synapses (skip last layer)
	int syn;
	float impulse;
	float decayRate = DT/mp->tauC;
	// Loop over efferent synapses
	/*for (syn=0; syn<n->nFEff_E; syn++)
	{
		impulse = (t == next_spike(&n->FEffs_E[syn])) ? mp->alphaC : 0.0;
		n->FEffs_E[syn]->C += (impulse * (1 - n->FEffs_E[syn]->C_tm1) - (n->FEffs_E[syn]->C_tm1 * DT/mp->tauC));
	}*/
	// Loop over afferent synapses (skip first layer)
	for (syn=0; syn<n->nFAff_E; syn++)
	{
		impulse = (t == next_spike(n->FAffs_E[syn])) ? mp->alphaC : 0.0;
		n->FAffs_E[syn]->C += (impulse * (1 - n->FAffs_E[syn]->C_tm1) - (n->FAffs_E[syn]->C_tm1 * decayRate));
	}
	return;
}
	
void update_D(NEURON * n, tstep t)
{
	// ||--> Update D for current neuron's incoming synapses (skip first layer)
	float impulse = ((t == n->lastSpike) ? mp->alphaD : 0.0); // private
	n->D += (impulse * (1 - n->D_tm1) - (n->D_tm1 * DT/mp->tauD));
	return;
}


void init_queue(AXON *a)
{
	int bin;
	a->next = 0;
	a->last = a->size-1;
	a->count = 0;
	for (bin=0; bin<a->size; bin++)
		a->queue[bin] = -BIG;
	return; // Necessary?
}

inline void enqueue(AXON *a, tstep t)
{
	assert(a->count < a->size);
	a->last = (a->last+1) % a->size;
	a->queue[ a->last ] = t + a->delay;
	a->count++;
	return;
}

int dequeue(AXON *a)
{
	int t;
	assert(a->count > 0);
	t = a->queue[ a->next ];
	a->next = (a->next+1) % a->size;
	a->count--;
	return(t);
}

inline int next_spike(AXON * a)
{
	return a->queue[a->next];
}

bool isempty(AXON *a)
{
	return (a->count <= 0) ? true : false;
}

void print_queue(AXON *a)
{
	int i;
	i=a->next; 
	
	while (i != a->last) {
		printf("%d ",a->queue[i]); // Was %c
		i = (i+1) % a->size;
	}
	printf("%d ",a->queue[i]);
	printf("\n");
	return;
}
