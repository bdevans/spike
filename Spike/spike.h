/*
 *  spike.h
 *  Spike
 *
 *  Created by Ben Evans on 6/19/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _SPIKE_H
#define _SPIKE_H

/*** Include files ***/

#include <stdbool.h> // Necessary for bool variable type
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//#undef _OPENMP
#ifdef _OPENMP
#include <omp.h>
#endif

//#include "main.h"
//#include "utils.h"
#include "array_utils.h"
#include "parameters.h"
#include "globals.h"
#include "read_parameters.h" // One use of printIntArray()

/*********************/

/* Must use these multiline comments */
//extern struct SIMULATION SIM;
//SIMULATION SIM;

typedef enum {
	INHIB,
	EXCIT
} NTYPE;

typedef enum {
	Testing,
	Training
	//NoLearning, 	//Testing,
	//Learning		//Training
	//Continuous,
	//Settle
} LEARNREGIME;

typedef enum {
	Soft,	// Reset solution variables only
	Hard	// Reset Spike buffers too
} SETSV;

/*typedef enum {
	PreTraining,
	Training,
	Testing,
	EfETraining,
	EfETesting
} PHASE; */

//int output_spikes[TBINS][NEXCIT]={0};

//#if RECF
//int REC_FLAGS[NLAYERS][NEXCIT]={0}; // Make this ragged-edged with NTYPES to save Inhibtory cells

typedef struct {
	int bin;
	float * cellV;
	float * cellcCa;
	float * cellD;
	float ** FSynC; // Feed-forward Afferents
	float ** FSynG;
	float ** FSynDG;
	float ** LSynC; // Lateral Afferents
	float ** LSynG;
	float ** LSynDG;
} RECORD;

//RECORD *RECSP[NLAYERS][NRECORDS_PL];
//#endif


// Place in affNeurons structure?
int *** affNeurons_EfE;
int *** affNeurons_ElE;
int *** affNeurons_IE;
int *** affNeurons_EI;
int *** affNeurons_II;

typedef struct AXON {
	tstep delay; // Actual axonal delay in timesteps
	tstep * queue;
	unsigned short int nBins;
	unsigned short int count;
	unsigned short int next;
	unsigned short int last;
	
	/* Synapse */
	float C, C_tm1;
	float g, g_tm1;
	float delta_g, delta_g_tm1;		/* Constrained to lie in the interval [0,1] */
} AXON;

/*// Non plastic Axons
typedef struct AXON_NP {
	tstep delay; // Actual axonal delay in timesteps
	unsigned short int size;
	tstep * queue;
	unsigned short int count;
	unsigned short int next;
	unsigned short int last;
	
	// Synapse
	float g, g_tm1;
	float delta_g;		// Constrained to lie in the interval [0,1]
} AXON_NP;*/

//typedef struct NEURON; 
typedef struct NEURON {
	NTYPE type;
	int l; // Layer index
	int n; // Index ?? l_ind, n_ind
	int row; // Row index (topological)
	int col; // Column index (topological)
	//int scale;
	//int orient;
	//int phase;
	float x; // Neuron centre
	float y; // Neuron centre
	float V, V_tm1;
	float D, D_tm1; // afferent
	float cCa, cCa_tm1;
	
	int spkbin;
	tstep * spikeTimes;	// Array of spike times
	tstep lastSpike;
	tstep nextUpdate;
	
	// Afferent (presynaptic) connections
	struct NEURON ** lm1presyn_E;	// Array of pointers to presynaptic excitatory neurons
	struct NEURON ** lm0presyn_E;
	struct NEURON ** lm0presyn_I;	// Array of pointers to presynaptic inhibitory neurons
	
	// Efferent (post-synaptic) connections
	struct NEURON ** lp1postsyn_E;
	struct NEURON ** lp0postsyn_E;
	struct NEURON ** lp0postsyn_I;
	
	int nFAff_E; // Counters for the number of presynaptic neurons
	//int nAffs_Ef; F&R!!
	int nLAff_E;
	int nLAff_I;

	//AXON ** aff_Ef;
	/* Afferent axons */
	AXON ** FAffs_E;	// Array of pointers to feed-forward excitatory afferent connections
	AXON ** LAffs_E;	// Array of lateral excitatory afferent connections
	AXON ** LAffs_I;	// Array of lateral inhibitory afferent connections
	
	//int nEffs_fE
	int nFEff_E;
	int nLEff_E;
	int nLEff_I;
	
	//AXON * Eff_Ef
	/* Efferent axons */
	AXON * FEffs_E;
	AXON * LEffs_E;
	AXON * LEffs_I;
	
	bool rec_flag;
	RECORD * rec; 	// Recording structure
} NEURON;

NEURON ** n_E;
NEURON ** n_I;

/*typedef struct {
//	int gr;
	int st;
	int tr;
} SCHEDULE;*/

typedef struct { //STIMULI {
	int nStim;
	int nTrans;
	int nTestStim;
	int nTestTrans;
	bool newTestSet;
	int ** groups; //bool
	float *** trn_stimuli;
	float *** tst_stimuli;
	//float **** trnGrpStim;
	//float **** tstGrpStim;
	int ** stimShuffle;
	int *** transShuffle;
	//SCHEDULE *** sched;
	float ******* trnImages; // Arrays for filtered images
	float ******* tstImages; // Arrays for filtered images
} STIMULI;

float *** distE;

// Array of arrays of structures?
// Expand to RECORDS[NLAYERS][NTYPES][NRECORDS] where NTYPES:= {0,1} (Excitatory, Inhibitory)

extern char * trim(char * string);

/* Prototypes from spike.c */

extern int spike(PARAMS * mp);
extern void calcMemory(PARAMS * mp);
extern NEURON ** allocn(int nlays, int * vNeurons, NTYPE type);
extern int unallocn(NEURON ** narray, int nlays, int * vNeurons);
extern void calcConnectivity(bool probConnect);
extern void wireAfferents(NEURON * n, float pEfn, float pEln, float pIln);
//extern void calc_connectivity(void);
extern float calcDistance(NEURON * n1, NEURON * n2, float scale); // DIM * D);
//extern void wire_afferents(NEURON * n, int l, int * affNcnx_fE, int * affNcnx_lE, int * affNcnx_I);
extern void alloc_efferents(NEURON * n);
extern void wire_efferents(NEURON * n);
extern void create_axons(NEURON * n);
extern void initNetwork(SETSV resetBuffers);
extern void setRecords(PARAMS * mp, NEURON ** n_E, gsl_rng * mSeed);
extern void setWeights(PARAMS * mp, NEURON ** n_E, NEURON ** n_I, const char * suffix);
extern int loadImages(STIMULI * stim, PARAMS * mp);
extern int loadDoGoutput(const char * filename, float * array, int size);
extern int loadGaborOutput(const char * filename, float * array, int size);
extern int loadGroups(STIMULI * stim, PARAMS * mp, char * filename);
extern void genGroups(STIMULI * stim, PARAMS * mp);
extern void printGroups(STIMULI * stim, PARAMS * mp, const char * filename);
extern void gen_stimuli(bool rep, STIMULI * stim, PARAMS * mp);
extern void printStimuli(STIMULI * stim, PARAMS * mp);
extern int genShuffles(STIMULI * stim, PARAMS * mp);
extern void printSchedule(STIMULI * stim, const char * filename);
extern void calcInput(PARAMS * mp, int loop, int pat, int trans, STIMULI * stim, float ** input, int regime);
extern void loadAfferents(const char * suffix);
extern void simulatePhase(LEARNREGIME regime, const char * prefix, STIMULI * stim);
extern void updateNetwork(tstep t_start, tstep t_end, float input[], int regime); //int loop, 
//extern void update_network(tstep t, int loop, float input[], int regime);
extern void update_V(NEURON * n, tstep t, float decay_rate, float gLeak, float Vrest, float Thresh, float Vhyper, float inj);
extern inline void update_cCa(NEURON * n, tstep t, float decayRate);
extern inline void update_g(NEURON * n, tstep t, float decay_E, float decay_I);
extern inline void update_weights(NEURON * n, tstep t);
extern inline void update_C(NEURON * n, tstep t, float decayRate);
extern inline void update_D(NEURON * n, tstep t, float decayRate);
extern int normalise(NEURON ** narray, PARAMS * mp);
extern inline void init_queue(AXON * a);
extern inline void enqueue(AXON * a, tstep t);
extern inline int dequeue(AXON * a);
extern inline int next_spike(AXON * a);
extern inline bool isempty(AXON * a);
extern inline int print_queue(AXON * a);

#endif
