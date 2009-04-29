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
	tstep t;	//unsigned long long t = 0; // Update spikeTimes and set spikeTimes[0] = -BIG and update all functions of t to tstep t not int t
	tstep t_start, t_end;
	int loop = 0;
	int trans = 0;
	int n, p, syn, l, wl;
	REGIMETYPE regime;
	char *plural = "";
	
	char filename[FNAMEBUFF];

	int *** trn_stimuli, *** tst_stimuli; // Place this group in a structure
	int ** shuffle;
	int * input;
	
	//PARAMS * mp; // const?
	
	//RECORD *RECSP = RECS;
	//RECORD *ptr = &RECS[0][0];
	//RECORD *r_ptr = RECS;
	
	/*** Declare file pointers ***/
	//FILE * parameters_ptr;
	FILE * output_spikes_ptr;
	FILE * weights_ptr;
	FILE * stimuli_FP;
	FILE * connections_FP;
	
	FILE * pt_weights_ptr;
	FILE * pt_output_spikes_ptr;
	
	//FILE *r_flags_ptr;
	FILE * r_g_ptr;
	FILE * r_cellV_ptr;
	FILE * r_weights_ptr;
	FILE * r_C_ptr;
	FILE * r_D_ptr;
	
	printf("Ventral Visual Stream Spiking Neural Network Simulation starting...\n");
	
	/*************** Build Network ****************/
	
	printf("\tNow building the network...");
	
	n_E = allocn(mp->nLayers, mp->nExcit, EXCIT); // Create 2D array of Excitatory neuron structures
	n_I = allocn(mp->nLayers, mp->nInhib, INHIB); // Create 2D array of Inhibatory neuron structures
	calc_connectivity();		// Calculate network connectivity
	
	/************** FILE OUTPUT **************/
	
	for (l=0; l<mp->nWLayers; l++) // Loop up to nWLayers
	{
		sprintf(filename, "L%daffNeuronsEfE.dat", l+1);
		connections_FP = fopen(filename, "w");
		print_iarray(connections_FP, affNeurons_EfE[l], mp->nExcit, mp->nSynEfE);
		fclose(connections_FP);
	}
	
	for (l=0; l<mp->nLayers; l++)
	{
		sprintf(filename, "L%daffNeuronsElE.dat", l);
		connections_FP = fopen(filename, "w");
		print_iarray(connections_FP, affNeurons_ElE[l], mp->nExcit, mp->nSynElE);
		fclose(connections_FP);
	}
	
	for (l=0; l<mp->nLayers; l++)
	{
		sprintf(filename, "L%daffNeuronsEI.dat", l);
		connections_FP = fopen(filename, "w");
		print_iarray(connections_FP, affNeurons_EI[l], mp->nInhib, mp->nSynEI);
		fclose(connections_FP);
	}
	
	for (l=0; l<mp->nLayers; l++)
	{
		sprintf(filename, "L%daffNeuronsIE.dat", l);
		connections_FP = fopen(filename, "w");
		print_iarray(connections_FP, affNeurons_IE[l], mp->nExcit, mp->nSynIE);
		fclose(connections_FP);
	}
	
	for (l=0; l<mp->nLayers; l++)
	{
		sprintf(filename, "L%daffNeuronsII.dat", l);
		connections_FP = fopen(filename, "w");
		print_iarray(connections_FP, affNeurons_II[l], mp->nInhib, mp->nSynII);
		fclose(connections_FP);
	}

	/************** END OF FILE OUTPUT **************/
	
	printf("\tBuilding complete!\n");
#if DEBUG > 0
	printf("NEURON structures require %ld bytes.\n",sizeof(NEURON));
	printf("An int requires %ld bytes.\t A tstep requires %ld bytes.\n",sizeof(int), sizeof(tstep));
#endif
	printf("\tNow initialising the network...");
	
	regime = Learning;			// 0: Testing (No STDP); 1: Training (STDP);
	init_network(regime);		// Initialise parameters
	
	printf("\tNetwork initialisation complete!\n");
	
	printf("\tNow creating the stimuli...");
	
	gen_stimuli(mp->localRep, &trn_stimuli, &tst_stimuli, &input, &shuffle);			// Generate Patterns
	
	stimuli_FP = fopen("trn_stimuli.dat", "w");
	for (p=0; p<mp->nStimuli; p++)
	{
		fprintf(stimuli_FP, "*** Stimulus %d ***\n", p+1);
		print_iarray(stimuli_FP, trn_stimuli[p], mp->nTransPS, mp->nExcit);
		fprintf(stimuli_FP, "\n");
	}
	fclose(stimuli_FP);
	
	stimuli_FP = fopen("tst_stimuli.dat", "w");
	for (p=0; p<mp->nStimuli; p++)
	{
		fprintf(stimuli_FP, "*** Stimulus %d ***\n", p+1);
		print_iarray(stimuli_FP, tst_stimuli[p], mp->nTransPS, mp->nExcit);
		fprintf(stimuli_FP, "\n");
	}
	fclose(stimuli_FP);
	
	printf("\tStimuli now saved!\n");

	if (mp->pretrain) // == 1
	{
		/**************** Pretraining ****************/
		
		regime = NoLearning;	// 0: Testing (No STDP); 1: Training (STDP);
		init_network(regime);
		
		printf("\tNow beginning pretraining phase...\n"); 		// Burn in period of random activity ?
		for (p=0; p<mp->nStimuli; p++) 		// Present each transform of each stimulus once
		{
			for (trans=0; trans<mp->nTransPS; trans++)
			{
				calc_input(loop, p, trans, tst_stimuli, input, shuffle, regime); // Move print statements outside
				t_start = round(((mp->transP_Test * trans) + (mp->transP_Test * mp->nTransPS * p)) * ceil(1/DT));
				t_end = t_start + round(mp->transP_Test * ceil(1/DT));
#if DEBUG>0 // Level 1
				fprintf(stderr, "\nUpdating network from timestep %lld to %lld.\n", t_start, t_end-1);
				for (n=0; n<mp->nExcit; n++)
					fprintf(stderr, "%d ",input[n]);
				fprintf(stderr, "\n");
#endif
//#pragma omp barrier
#pragma omp parallel default(shared) private (t) // Parallelize p loop?
				for (t=t_start; t<t_end; t++)
					update_network(t, loop, input, regime);
			}
		}
/*#pragma omp parallel default(shared) private (t)
		{
		for (t=0; t<TEST_PERIOD_TS; t++) // Think about the presentation time
		{
#pragma omp single
			{
			if ((t % cueTS) == 0) // Impose stimulus
				calc_input(t, loop, tst_stimuli, input, regime); // Calculate input pattern
			}
			update_network(t, loop, input, regime);
		}
		}*/
		printf("\tTesting complete!\n");
		
		printf("\tSaving results...");
		// Output results to dat files
	
		for (l=0; l<mp->nLayers; l++)
		{
			sprintf(filename, "ptL%dExcitSpikes.dat", l);
			pt_output_spikes_ptr = fopen(filename, "w");
			for (n=0; n<mp->nExcit; n++)
				print_irow(pt_output_spikes_ptr, n_E[l][n].spikeTimes, n_E[l][n].spkbin+1); // spikeTimes[0] = -BIG?
			fclose(pt_output_spikes_ptr);
		}
		// fprintf(pt_output_spikes_ptr, "L%dN%d\t", NLAYERS-1, n);
		
		for (wl=1; wl<mp->nLayers; wl++)
		{
			sprintf(filename, "ptL%dweightsEfE.dat", wl);
			pt_weights_ptr = fopen(filename, "w");
			for (n=0; n<mp->nExcit; n++)
			{
				for (syn=0; syn<n_E[wl][n].nFAff_E; syn++)
					fprintf(pt_weights_ptr, "%f\t", n_E[wl][n].FAffs_E[syn]->delta_g);
				fprintf(pt_weights_ptr, "\n");
			}
			fclose(pt_weights_ptr);
		}
		
		printf("\tResults saved!\n"); // Print path
		
	}
		/**************** End of pretraining ****************/
	
	if (mp->train)
	{
		/*************** Learning phase ***************/
		regime = Learning;
		printf("\tProgram now beginning learning phase...\n");
		plural = ((mp->loops==1) ? "" : "s");
		printf("\tThe program will now simulate %2.2f s each loop for %d loop%s.\n", mp->TotalTime, mp->loops, plural);
		
		for (loop=0; loop < mp->loops; loop++)
		{
			init_network(0);	// Reinitialise parameters here to reset V's, g's, C's, D's and spike buffers
			printf("\tLoop %d/%d\n", loop+1, mp->loops);
			
			int nPatterns = (!mp->trainPause) ? mp->nStimuli : (2*mp->nStimuli)-1;
			for (p=0; p<nPatterns; p++)
			{
				for (trans=0; trans<mp->nTransPS; trans++)
				{
					if (mp->trainPause && (p % 2 == 1))
					{
						memset(input, 0, mp->nExcit * sizeof(*input));
						printf("\t\tLetting the network settle for %.2fs after stimulus %d...\n",mp->transP_Train,(p+1)/2);
					}
					else
						calc_input(loop, (mp->trainPause)?(p/2)+1:p, trans, trn_stimuli, input, shuffle, regime);
					t_start = round(((mp->transP_Train * trans) + (mp->transP_Train * mp->nTransPS * p)) * ceil(1/DT));
					t_end = t_start + round((mp->transP_Train * ceil(1/DT)));
#if DEBUG>0 // Level 1
					fprintf(stderr, "\nUpdating network from timestep %lld to %lld.\n", t_start, t_end-1);
					for (n=0; n<mp->nExcit; n++)
						fprintf(stderr, "%d ",input[n]);
					fprintf(stderr, "\n");
#endif
#pragma omp parallel default(shared) private (t) // Parallelize p loop?
					for (t=t_start; t<t_end; t++)
						update_network(t, loop, input, regime);
				}
			}
		}
		printf("\tThe network has been trained!\n");
	}
	
	/*************** Testing phase ***************/
	
	regime = NoLearning;		// 0: Testing (No STDP); 1: Training (STDP);
	loop = 0;
	init_network(regime);
	
	printf("\tNow beginning testing phase...\n");
	for (p=0; p<mp->nStimuli; p++)
	{
		for (trans=0; trans<mp->nTransPS; trans++)
		{
			calc_input(loop, p, trans, tst_stimuli, input, shuffle, regime);
			t_start = round(((mp->transP_Test * trans) + (mp->transP_Test * mp->nTransPS * p)) * ceil(1/DT));
			t_end = t_start + round((mp->transP_Test * ceil(1/DT)));
#pragma omp parallel default(shared) private (t) // Parallelize p loop?
			for (t=t_start; t<t_end; t++)
				update_network(t, loop, input, regime);
		}
	}
	printf("\tTesting complete!\n");

	printf("\tSaving results...");
	// Output results to dat files
	
	// Save Excitatory spikes for all layers
	for (l=0; l<mp->nLayers; l++)
	{
		snprintf(filename, FNAMEBUFF, "L%dExcitSpikes.dat", l);
		output_spikes_ptr = fopen(filename, "w"); // Make either HR "w" or Binary "wb"
		for (n=0; n<mp->nExcit; n++)
			print_irow(output_spikes_ptr, n_E[l][n].spikeTimes, n_E[l][n].spkbin+1);
		fclose(output_spikes_ptr);
	}
	
	// Save Inhibitory spikes for all layers
	for (l=0; l<mp->nLayers; l++)
	{
		snprintf(filename, FNAMEBUFF, "L%dInhibSpikes.dat", l); 		//asprintf(&filename, "spikes_I_L%d.dat", l);
		output_spikes_ptr = fopen(filename, "w"); // Make either HR "w" or Binary "wb"
		for (n=0; n<mp->nInhib; n++)
			print_irow(output_spikes_ptr, n_I[l][n].spikeTimes, n_I[l][n].spkbin+1);
		fclose(output_spikes_ptr);
	}
	
	// Save weights for EFE synapses
	for (wl=1; wl<mp->nLayers; wl++) // Start at 1
	{
		sprintf(filename, "L%dweightsEfE.dat", wl);
		weights_ptr = fopen(filename, "w");
		for (n=0; n<mp->nExcit; n++)
		{
			for (syn=0; syn<n_E[wl][n].nFAff_E; syn++)
				fprintf(weights_ptr, "%f\t", n_E[wl][n].FAffs_E[syn]->delta_g);
			fprintf(weights_ptr, "\n");
		}
		fclose(weights_ptr);
	}
	
	/* Save weights in binary format */
	/*bweights_ptr = fopen("bweights.dat", "wb");
	fwrite(delta_g_EE, sizeof(double), sizeof(delta_g_EE)/sizeof(double), bweights_ptr); // See nifty trick #1
	fclose(bweights_ptr);*/
	
	/* Save recorded values */
	if (mp->nRecordsPL)
	{
		int loop=0;
		for (l=0; l<mp->nLayers; l++)
		{
			for (n=0; n<mp->nExcit; n++)
			{
				if (n_E[l][n].rec_flag)
				{
					snprintf(filename, FNAMEBUFF, "rL%dN%dcellV.dat", l, n); 	//asprintf(&filename, "r_cellV_L%dR%d.dat", l, n);
					r_cellV_ptr = fopen(filename, "w"); //r_cellV_ptr = fopen(filename, "wb");
					//fwrite(n_E[l][n].rec->cellV, sizeof(float), mp->loops*mp->TotalMS, r_cellV_ptr); // See nifty trick #1
					print_farray(r_cellV_ptr, n_E[l][n].rec->cellV, mp->loops, mp->TotalMS);
					fclose(r_cellV_ptr);
					//free(filename);
					
					//if (l<mp->nWLayers)
					
					if (l>0)
					{
						snprintf(filename, FNAMEBUFF, "rL%dN%dD.dat", l, n);
						r_D_ptr = fopen(filename, "w"); //r_D_ptr = fopen(filename, "wb");
						print_farray(r_D_ptr, n_E[l][n].rec->cellD, mp->loops, mp->TotalMS);
						//fwrite(n_E[l][n].rec->cellD, sizeof(float), mp->loops*mp->TotalMS, r_D_ptr); // See nifty trick #1 //sizeof(n_E[l][n].rec->cellD)/sizeof(float)
						fclose(r_D_ptr);
						
						// The presynaptic cell's values are attached to each record
						snprintf(filename, FNAMEBUFF, "rL%dN%dC.dat", l, n);
						r_C_ptr = fopen(filename, "w"); //r_C_ptr = fopen(filename, "wb");
						for (loop=0; loop<mp->loops; loop++)
							print_farray(r_C_ptr, n_E[l][n].rec->SynC[loop], n_E[l][n].nFAff_E, mp->TotalMS);
						//fwrite(n_E[l][n].rec->SynC, sizeof(float), mp->loops*n_E[l][n].nFAff_E*mp->TotalMS, r_C_ptr); //fwrite(n_E[l][n].rec->SynC, sizeof(float), sizeof(n_E[l][n].rec->SynC)/sizeof(float), r_C_ptr); // See nifty trick #1
						fclose(r_C_ptr);
						
						snprintf(filename, FNAMEBUFF, "rL%dN%dg.dat", l, n);
						r_g_ptr = fopen(filename, "w");
						for (loop=0; loop<mp->loops; loop++)
							print_farray(r_g_ptr, n_E[l][n].rec->SynG[loop], n_E[l][n].nFAff_E, mp->TotalMS);
						fclose(r_g_ptr);
						
						snprintf(filename, FNAMEBUFF, "rL%dN%ddg.dat", l, n);
						r_weights_ptr = fopen(filename, "w");
						for (loop=0; loop<mp->loops; loop++)
							print_farray(r_weights_ptr, n_E[l][n].rec->SynDG[loop], n_E[l][n].nFAff_E, mp->TotalMS);
						fclose(r_weights_ptr);
					}
				}
			}
			/*
			for (n=0; n<mp->nRecordsPL; n++)
			{
				r_flags_ptr = fopen("r_flags.dat", "wb");
				fwrite(REC_FLAGS, sizeof(int), sizeof(REC_FLAGS)/sizeof(int), r_flags_ptr); // See nifty trick #1
				fclose(r_flags_ptr);
				
				snprintf(filename, FNAMEBUFF, "r_cellV_L%dR%d.dat", l, n);
				//asprintf(&filename, "r_cellV_L%dR%d.dat", l, n);
				r_cellV_ptr = fopen(filename, "wb");
				fwrite(RECSP[l][n]->cellV, sizeof(double), LOOPS*TOT_MS, r_cellV_ptr); // See nifty trick #1
				fclose(r_cellV_ptr);
				//free(filename);
				
				if (l<NWLAYERS)
				{
					snprintf(filename, FNAMEBUFF, "r_C_L%dR%d.dat", l, n);
					//asprintf(&filename, "r_C_L%dR%d.dat", l, n);
					r_weights_ptr = fopen(filename, "wb");
					r_C_ptr = fopen(filename, "wb");
					fwrite(RECSP[l][n]->C, sizeof(double), sizeof(RECSP[l][n]->C)/sizeof(double), r_C_ptr); // See nifty trick #1
					fclose(r_C_ptr);
					//free(filename);
				}
				
				if (l>0)
				{
					//wl = l-1;
					snprintf(filename, FNAMEBUFF, "r_D_L%dR%d.dat", l, n);
					//asprintf(&filename, "r_D_L%dR%d.dat", l, n);
					r_weights_ptr = fopen(filename, "wb");
					r_D_ptr = fopen(filename, "wb");
					fwrite(RECSP[l][n]->D, sizeof(double), sizeof(RECSP[l][n]->D)/sizeof(double), r_D_ptr); // See nifty trick #1
					fclose(r_D_ptr);
					//free(filename);
					
					snprintf(filename, FNAMEBUFF, "r_weights_L%dR%d.dat", l, n);
					//asprintf(&filename, "r_weights_L%dR%d.dat", l, n); // l or wl?
					r_weights_ptr = fopen(filename, "wb");
					//r_weights_ptr = fopen("r_weights.dat", "wb");
					//fwrite(RECSP[l][n]->delta_g, sizeof(double), sizeof(RECSP[l][n]->delta_g)/sizeof(double), r_weights_ptr);
					fwrite(RECSP[l][n]->delta_g, sizeof(double), LOOPS*NSYN_PN*TOT_MS, r_weights_ptr); // See nifty trick #1
					fclose(r_weights_ptr);
					//free(filename);
				}
			}*/
		}
	}
		
	printf("\tResults saved!\n"); // Print path
	
	//unallocn(n_E, mp->nLayers, mp->nExcit);
	//unallocn(n_I, mp->nLayers, mp->nInhib);
	// Free other arrays...

	return 0;
}


NEURON ** allocn (int nlays, int nneurons, NTYPE type)
{
	int l, n;
	NEURON *space;
	NEURON **narray;
	
	// sizeof(*array) = sizeof(int **)
	// sizeof(**array) = sizeof(int *)
	// sizeof(***array) = sizeof(int)
	
	/*** Ensures array is contiguous in memory ***/
	space = myalloc(nlays * nneurons * sizeof(*space));
	/*********************************************/
	narray = myalloc(nlays * sizeof(*narray));

	for (l=0; l<nlays; l++)
		narray[l] = space + (l * nneurons); 
	
	/* Allocate spike time bins */
	for (l=0; l<nlays; l++)
		for (n=0; n<nneurons; n++)
		{
			narray[l][n].spikeTimes = myalloc(mp->spkBuffer * sizeof(int)); //sizeof(int)//sizeof(narray[l][n].spikeTimes[0])
			narray[l][n].type = (type==EXCIT) ? EXCIT : INHIB;
		}
	
	return narray;
}


int unallocn (NEURON ** narray, int nlays, int nneurons)
{
	int l, n, s;
	for (l=0; l<nlays; l++)
	{
		for (n=0; n<nneurons; n++)
		{
			free(narray[l][n].spikeTimes);
			for (s=0; s<narray[l][n].nFEff_E; s++)
				free(narray[l][n].FEffs_E[s].queue);
			if (narray[l][n].FEffs_E != NULL)
				free(narray[l][n].FEffs_E);
			for (s=0; s<narray[l][n].nLEff_E; s++)
				free(narray[l][n].LEffs_E[s].queue);
			if (narray[l][n].LEffs_E != NULL)
				free(narray[l][n].LEffs_E);
			for (s=0; s<narray[l][n].nLEff_I; s++)
				free(narray[l][n].LEffs_I[s].queue);
			if (narray[l][n].LEffs_I != NULL)
				free(narray[l][n].LEffs_I);
			if (narray[l][n].rec_flag)
			{
				free_2D_farray(narray[l][n].rec->cellV, mp->loops);
				if (l>0)
					free_2D_farray(narray[l][n].rec->cellD, mp->loops);
				if (l<mp->nWLayers)
				{
					free_3D_farray(narray[l][n].rec->SynC, mp->loops, narray[l][n].nFAff_E);
					free_3D_farray(narray[l][n].rec->SynG, mp->loops, narray[l][n].nFAff_E);
					free_3D_farray(narray[l][n].rec->SynDG, mp->loops, narray[l][n].nFAff_E);
				}
				free(narray[l][n].rec);
			}
		}
		free(narray[l]);
	}
	free(narray);
	return 0;
}


void calc_connectivity() /* This function randomly connects the neurons together */
{
	//int n, s, l, wl;
	//int choice;
	int n, l;
	
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
			wire_afferents(&n_E[l][n], l, (l>0)?affNeurons_EfE[l-1][n]:NULL, affNeurons_ElE[l][n], affNeurons_IE[l][n]);
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
	for (n=0; n<mp->nExcit; n++)
		for (l=0; l<mp->nLayers; l++)
			create_axons(&n_E[l][n]);
	
	for (n=0; n<mp->nInhib; n++)
		for (l=0; l<mp->nLayers; l++)
			create_axons(&n_I[l][n]);
		
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
	n->nFEff_E = 0; // Reset to as a synapse counter in wire_afferents
	
	//if (n->LEffs_E > 0)
	n->LEffs_E = myalloc(n->nLEff_E * sizeof(AXON));
	n->lp0postsyn_E = myalloc(n->nLEff_E * sizeof(NEURON *));
	n->nLEff_E = 0; // Reset to as a synapse counter in wire_afferents
	
	//if (n->LEffs_I > 0)
	n->LEffs_I = myalloc(n->nLEff_I * sizeof(AXON));
	n->lp0postsyn_I = myalloc(n->nLEff_I * sizeof(NEURON *));
	n->nLEff_I = 0; // Reset to as a synapse counter in wire_afferents
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
	int delay, nbins, span;
	int s;
	switch (mp->axonDelay) 
	{
		case MinD:
#if DEBUG > 1
			printf("\nSetting axonal delays to minimum\n");
#endif
			delay = 1;
			nbins = 1;
			break;
		case ConstD:
#if DEBUG > 1			
			printf("\nSetting axonal delays to %f seconds\n",mp->d_const);
#endif
			delay = ceil(mp->d_const/DT);
			nbins = ceil(mp->d_const/mp->refract);
			break;
		case UniformD:
#if DEBUG > 1
			printf("\nDrawing axonal delays from [%f, %f]\n",mp->d_min, mp->d_max);
#endif
			span = mp->d_max - mp->d_min;
			break;
		case GaussD:
#if DEBUG > 1
			printf("\nDrawing axonal delays from N(%f,%f)\n",mp->d_mean, mp->d_sd);
#endif
			break;
		default:
			printf("\nError: Unknown axonal delay model!\n");
			break;
	}
	
	for (s=0; s<n->nFEff_E; s++) // nFEff_E = 0 for l = nLayers and n->type == INHIB
	{
		switch (mp->axonDelay) 
		{
			case UniformD:
				delay = round(((ran3(&idum)*span)+mp->d_min)/DT);
				nbins = ceil((delay*DT)/mp->refract);
				break;
			case GaussD:
				delay = round((mp->d_mean + mp->d_sd*gasdev(&idum))/DT);
				nbins = ceil((delay*DT)/mp->refract);
				break;
			default:
				break;
		}
		n->FEffs_E[s].delay = delay;
		n->FEffs_E[s].size = nbins;
		n->FEffs_E[s].queue = myalloc(nbins * sizeof(int));
		init_queue(&(n->FEffs_E[s]));
	}
	for (s=0; s<n->nLEff_E; s++)
	{
		switch (mp->axonDelay) 
		{
			case UniformD:
				delay = round(((ran3(&idum)*span)+mp->d_min)/DT);
				nbins = ceil((delay*DT)/mp->refract);
				break;
			case GaussD:
				delay = round((mp->d_mean + mp->d_sd*gasdev(&idum))/DT);
				nbins = ceil((delay*DT)/mp->refract);
				break;
			default:
				break;
		}
		n->LEffs_E[s].delay = delay;
		n->LEffs_E[s].size = nbins;
		n->LEffs_E[s].queue = myalloc(nbins * sizeof(int));
		init_queue(&(n->LEffs_E[s]));
	}
	for (s=0; s<n->nLEff_I; s++)
	{
		switch (mp->axonDelay) 
		{
			case UniformD:
				delay = round(((ran3(&idum)*span)+mp->d_min)/DT);
				nbins = ceil((delay*DT)/mp->refract);
				break;
			case GaussD:
				delay = round((mp->d_mean + mp->d_sd*gasdev(&idum))/DT);
				nbins = ceil((delay*DT)/mp->refract);
				break;
			default:
				break;
		}
		n->LEffs_I[s].delay = delay;
		n->LEffs_I[s].size = nbins;
		n->LEffs_I[s].queue = myalloc(nbins * sizeof(int));
		init_queue(&(n->LEffs_I[s]));
	}
	return;
}

void init_network(int regime)
{

	int s, n, l;
	int r; //, syn, t;
	//int loop;
	int choice;
	//int nbins, delay;
	//float span;

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
					choice=rands1_new(0,mp->nExcit-1,&iv,1);
					n_E[l][choice].rec_flag = 1;
					printf("\tLayer %d, Record %d Assigned nID: %d\n",l,r,choice);
				}
			}
			/* Create recording structures */
			for (n=0; n<mp->nExcit; n++)
				for (l=0; l<mp->nLayers; l++)
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
			// Record the first NSYN_PN synapses (since they are randomly connected)
			/* Copy first array onto all others */
			/*for (l=0; l<NLAYERS; l++)
			{
				n = ((l==0) ? 1 : 0);
				while (n < NRECORDS_PL)
					*RECSP[l][n++] = *RECSP[0][0];
			}*/
		}
		
		/* E_ Synaptic weights */
		for (n=0; n<mp->nExcit; n++)
			for (l=0; l<mp->nWLayers; l++)
				for (s=0; s<n_E[l][n].nFEff_E; s++)
					n_E[l][n].FEffs_E[s].delta_g = n_E[l][n].FEffs_E[s].delta_g_tm1 = ran3(&idum);
		
		for (n=0; n<mp->nExcit; n++)
			for (l=0; l<mp->nLayers; l++)
				for (s=0; s<n_E[l][n].nLEff_E; s++)
					n_E[l][n].LEffs_E[s].delta_g = n_E[l][n].LEffs_E[s].delta_g_tm1 = ran3(&idum);
		
		for (n=0; n<mp->nExcit; n++)
			for (l=0; l<mp->nLayers; l++)
				for (s=0; s<n_E[l][n].nLEff_I; s++)
					n_E[l][n].LEffs_I[s].delta_g = n_E[l][n].LEffs_I[s].delta_g_tm1 = mp->Dg_EI;
		
		/* I_ Synaptic weights */
		for (n=0; n<mp->nInhib; n++)
			for (l=0; l<mp->nLayers; l++)
				for (s=0; s<n_I[l][n].nLEff_E; s++)
					n_I[l][n].LEffs_E[s].delta_g = n_I[l][n].LEffs_E[s].delta_g_tm1 = mp->Dg_IE;
		
		for (n=0; n<mp->nInhib; n++)
			for (l=0; l<mp->nLayers; l++)
				for (s=0; s<n_I[l][n].nLEff_I; s++)
					n_I[l][n].LEffs_I[s].delta_g = n_I[l][n].LEffs_I[s].delta_g_tm1 = mp->Dg_II;
		
		
		switch (mp->noise)
		{
			case 0: // No noise // Already initialised to 0 when declared
				break;
			case 1: // Uniform interval
				break;
			case 2: // Gaussian noise
				break;
		} // End of NOISE switch		
	} // End of if (regime==1) clause
	

	
	/********** Training and Testing **********/
	
	/* Membrane potentials, conductances, learning parameters and spike time arrays */
	for (n=0; n<mp->nExcit; n++)
	{
		for (l=0; l<mp->nLayers; l++)
		{
			n_E[l][n].V = n_E[l][n].V_tm1 = mp->VrestE;
			n_E[l][n].D = n_E[l][n].D_tm1 = 0.0;
			for (s=0; s<n_E[l][n].nFEff_E; s++)
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
			n_E[l][n].lastSpike = 0; //-BIG;
			n_E[l][n].spkbin = 0;
			memset(n_E[l][n].spikeTimes, 0, mp->spkBuffer*sizeof(n_E[l][n].spikeTimes[0])); //[0]?
			/*for (s=0; s<mp->spkBuffer; s++)
				n_E[l][n].spikeTimes[s] = 0;*/
		}
	}
	
	for (n=0; n<mp->nInhib; n++)
	{
		for (l=0; l<mp->nLayers; l++)
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
			n_I[l][n].lastSpike = 0; //-BIG;
			n_I[l][n].spkbin = 0;
			memset(n_I[l][n].spikeTimes, 0, mp->spkBuffer*sizeof(n_I[l][n].spikeTimes[0]));
		}
	}

	/* Recording structures */
	if (mp->nRecordsPL)
		for (n=0; n<mp->nExcit; n++)
			for (l=0; l<mp->nLayers; l++)
				if (n_E[l][n].rec_flag)
					n_E[l][n].rec->bin = 0;
	
	/*	for (l=0; l<NLAYERS; l++)
			for (n=0; n<NRECORDS_PL; n++)
				RECSP[l][n]->bin = 0; */
	
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

void gen_stimuli(bool rep, int ****trn_stimuli, int ****tst_stimuli, int ** input, int *** shuffle)
{
	// Place array arguements in a patterns structure
	int loop;
	int trans;
	int n, p;
	int choice;
	
	int block = 0;

	/* Generate the training stimuli to present to the network */
	// Distributed or orthogonal will depend on where the balls are replaced i.e. rng reinitialised...
	// Patterns could be generated in Matlab and loaded from dat files...
	// Or read in list of pairs from another array
	
	*input = myalloc(mp->nExcit * sizeof(**input));
	*trn_stimuli = get_3D_iarray(mp->nStimuli, mp->nTransPS, mp->nExcit, 0);
	*tst_stimuli = get_3D_iarray(mp->nStimuli, mp->nTransPS, mp->nExcit, 0);
	
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
					(*trn_stimuli)[p][trans][choice] = (*tst_stimuli)[p][trans][choice] = 1;
				}
			}			
			break;
			
		case 1: /* Local Patterns */
			block = floor(mp->nFiringNeurons + ((mp->nTransPS - 1) * mp->shift));
			// This assumes stimuli are a contiguous block of 1's
			for (p=0; p<mp->nStimuli; p++)
				for (trans=0; trans<mp->nTransPS; trans++)
					for (n=(p*block)+(trans*mp->shift); n<(mp->nFiringNeurons+(p*block)+(trans*mp->shift)); n++)
						(*trn_stimuli)[p][trans][n] = (*tst_stimuli)[p][trans][n] = 1;
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
	
	if (mp->random_order)
	{
		*shuffle = get_2D_iarray(mp->loops, mp->nStimuli, 0);
		for (loop=0; loop<mp->loops; loop++)
		{
			i=rands1_new(0,0,&iv,1);	// Initialise RNG (clear internal array of numbers drawn)
			for (p=0; p<mp->nStimuli; p++)
				(*shuffle)[loop][p] = rands1_new(0,mp->nStimuli-1,&iv,1);
		}
	}
	
	
	/*** Testing stimuli ***/
	// Modify to generate single stimuli test patterns when there are multi-stimulus training patterns
	
	//memcpy(*tst_stimuli, *trn_stimuli, sizeof(*tst_stimuli));
	
	
	return; // void;
}

void calc_input(int loop, int pat, int trans, int ***stimuli, int * input, int ** shuffle, int regime) // return int * input?
{	// Assumes all stimuli translate
	// No need to pass stimuli since tst_stimuli and trn_stimuli are global
	/*
	int p, pat;
	int trans;
	float transP;
	
	transP = ((regime==Learning) ? mp->transP_Train : mp->transP_Test);
	p = floor(t*DT/(transP*mp->nTransPS));
	
	trans = floor((t*DT)/transP);
	trans %= mp->nTransPS;
	*/
	switch (regime)
	{
		case 0: // Testing stimuli
			//pat=p;
			memcpy(input, stimuli[pat][trans], mp->nExcit*sizeof(input[0])); // Copy row of tst_stimuli to input array
			break;
			
		case 1: // Training stimuli
			pat=((mp->random_order) ? (shuffle[loop][pat]) : pat );
			memcpy(input, stimuli[pat][trans], mp->nExcit*sizeof(input[0]));
			break;
	}

	if (trans==0)
		printf("\t\tPresenting stimulus %d/%d...\n", (pat+1), mp->nStimuli);
	printf("\t\t\tTransform %d/%d...\n", (trans+1), mp->nTransPS);
	
	return; // void;
}

void update_network(int t, int loop, int input[], int regime)
{
	int n, syn, s, l; //, wl;
	//int bin = 0;
	float decay_rate, gLeak, Vrest, Thresh, Vhyper;
	float decay_E, decay_I;
	//int totSpikes_E, totSpikes_I; // Place totSpikes in bins

	//memcpy(DestArray, SourceArray, sizeof(DestArray));
/*
#pragma omp sections
	{
#pragma omp section
	memcpy(tm1cellV_E, cellV_E, sizeof(cellV_E));

	}
*/
	
	/* Update Excitatory cell potentials */
	decay_rate = DT/mp->capE;
	gLeak = mp->gLeakE;
	Vrest = mp->VrestE;
	Thresh = mp->ThreshE;
	Vhyper = mp->VhyperE;
#pragma omp for nowait private(n, l)
	for (n=0; n<mp->nExcit; n++)
		for (l=0; l<mp->nLayers; l++)
			update_V(t, decay_rate, gLeak, Vrest, Thresh, Vhyper, (l==0 && input[n])?1:0, &n_E[l][n]);
	
	/* Update Inhibitory cell potentials */
	decay_rate = DT/mp->capI;
	gLeak = mp->gLeakI;
	Vrest = mp->VrestI;
	Thresh = mp->ThreshI;
	Vhyper = mp->VhyperI;
#pragma omp for private(n, l)
	for (n=0; n<mp->nInhib; n++)
		for (l=0; l<mp->nLayers; l++)
			update_V(t, decay_rate, gLeak, Vrest, Thresh, Vhyper, 0, &n_I[l][n]);
//#pragma omp barrier

	/* Update presynaptic conductances */
	decay_E = (DT/mp->tauEE);
	decay_I = (DT/mp->tauIE);
#pragma omp for nowait private(n, l)
	for (n=0; n<mp->nExcit; n++)
		for (l=0; l<mp->nLayers; l++)
			update_g(&n_E[l][n], decay_E, decay_I, t);
	
	decay_E = (DT/mp->tauEI);
	decay_I = (DT/mp->tauII);
#pragma omp for private(n, l)
	for (n=0; n<mp->nInhib; n++)
		for (l=0; l<mp->nLayers; l++)
			update_g(&n_I[l][n], decay_E, decay_I, t);
//#pragma omp barrier
	
	if (regime==Learning) // Learning
	{
		/* Update synaptic weights */ // Move after C & D and use instantaeous values? N.B. Redo nowait clauses
#pragma omp for nowait private(n, l) // shared(n_E[][])
		for (n=0; n<mp->nExcit; n++) // Make parallel and n private
			for (l=1; l<mp->nLayers; l++) // Only nLayers-1 of weight layers
				update_weights(&n_E[l][n], t);
		
		// -->|| Update C for current neuron's outgoing synapses (skip last layer)
#pragma omp for nowait private(n, l) // shared(n_E[][])
		for (n=0; n<mp->nExcit; n++)
			for (l=1; l<mp->nLayers; l++) // Only nLayers-1 of synapses
				update_C(&n_E[l][n], t);
		
		// ||--> Update D for current neuron's incoming synapses (skip first layer)
#pragma omp for private(n, l) // shared(n_E[][]) // Last loop must be synchronised i.e. !(nowait)
		for (n=0; n<mp->nExcit; n++)
			for (l=1; l<mp->nLayers; l++) // Only nLayers-1 of weight layers
				update_D(&n_E[l][n], t);
//#pragma omp barrier
	}
		
	/* Copy solution variables to _tm1 counterparts and reset spike flags / axons */
#pragma omp for nowait private(l, n, s)
	for (n=0; n<mp->nExcit; n++)
	{
		for (l=0; l<mp->nLayers; l++)
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
		
#pragma omp for private(l, n, s)
	for (n=0; n<mp->nInhib; n++)
	{
		for (l=0; l<mp->nLayers; l++)
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
# pragma omp for private(bin, syn, l, n)
		for (n=0; n<mp->nExcit; n++) //for (r=0; r<NRECORDS_PL; r++)
			for (l=0; l<mp->nLayers; l++)
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
//#pragma omp barrier
	}

	return; // void;
}

void update_V(int t, float decay_rate, float gLeak, float Vrest, float Thresh, float Vhyper, int inj, NEURON * n)
{
	/* This function updates a neuron's cell potential and applies to all layers */
	int syn;
	float tot_g_E = 0.0;
	float tot_g_I = 0.0;
	float stim;
	if ((((t - n->lastSpike) * DT) > mp->refract) || n->spkbin==0)
	{
		/***** Add driving current or clamp voltage when the stimulus is presented to layer 0 *****/
		stim = (inj) ? mp->current : 0.0;
		
		/* Feed-forward connections */
		//if (n->type==EXCIT && l>0)
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
							  + stim);
							// + noise_E[n][l]);
		
		n->V = ((n->V < mp->VrevI) ? mp->VrevI : n->V); // Neurons can not be more -ve than Inhib reversal potential
		
		if (n->V >= Thresh) // tm1cellV_E? - gives one timestep to depolarise (peak)
		{
			n->V = Vhyper;
			n->lastSpike = n->spikeTimes[++(n->spkbin)] = t;
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


void update_g(NEURON * n, float decay_E, float decay_I, int t)
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

void update_weights(NEURON * n, int t)
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
	
void update_C(NEURON * n, int t)
{
	// -->|| Update C for current neuron's outgoing synapses (skip last layer)
	int syn;
	float impulse;
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
		n->FAffs_E[syn]->C += (impulse * (1 - n->FAffs_E[syn]->C_tm1) - (n->FAffs_E[syn]->C_tm1 * DT/mp->tauC));
	}
	return;
}
	
void update_D(NEURON * n, int t)
{
	// ||--> Update D for current neuron's incoming synapses (skip first layer)
	float impulse = ((t == n->lastSpike) ? mp->alphaD : 0.0); // private
	n->D += (impulse * (1 - n->D_tm1) - (n->D_tm1 * DT/mp->tauD));
	return;
}


//#include "queue.h"
//#include "bool.h"

/* Move to fifo.c */

void init_queue(AXON *a)
{
	int bin;
	//a->size = 0;
	//a->delay = 0;
	a->next = 0;
	a->last = a->size-1;
	a->count = 0;
	for (bin=0; bin<a->size; bin++)
		a->queue[bin] = -BIG;
	return; // Necessary?
}

void enqueue(AXON *a, int x)
{
	
	if (a->count >= a->size)
		printf("Warning: queue overflow enqueue x=%d\n",x);
	//else {
	 
	assert(a->count < a->size);
	a->last = (a->last+1) % a->size;
	a->queue[ a->last ] = x + a->delay;
	a->count++;
	//}
	return;
}

int dequeue(AXON *a)
{
	int x;
	/*
	if (q->count <= 0) printf("Warning: empty queue dequeue.\n");
	else {
	*/
	assert(a->count > 0);
	x = a->queue[ a->next ];
	a->next = (a->next+1) % a->size;
	a->count--; // = a->count - 1;
	//}
	
	return(x);
}

int next_spike(AXON * a)
{
	return a->queue[a->next];
}

int isempty(AXON *a)
{
	/*if(q->count <= 0) return (TRUE);
	else return (FALSE);*/
	return (a->count <= 0) ? 1 : 0;
}

void print_queue(AXON *a)
{
	int i; //,j;
	
	i=a->next; 
	
	while (i != a->last) {
		printf("%c ",a->queue[i]);
		i = (i+1) % a->size;
	}
	
	printf("%2d ",a->queue[i]);
	printf("\n");
	return;
}
