/*
 *  spike.h
 *  Spike
 *
 *  Created by Ben Evans on 6/19/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

/* Must use these multiline comments */

//int output_spikes[TBINS][NEXCIT]={0};

//#if RECF
//int REC_FLAGS[NLAYERS][NEXCIT]={0}; // Make this ragged-edged with NTYPES to save Inhibtory cells

typedef struct {
	int bin;
	float ** cellV; 	//double cellV[LOOPS][TOT_MS];
	float ** cellD;		//double D[LOOPS][TOT_MS];
	float *** SynC;		//double C[LOOPS][TOT_MS];
	float *** SynG;		//double synG[LOOPS][NSYN_EE][(int)(1000*TOT_TIME)]; // Generalise
	float *** SynDG; 	//double delta_g[LOOPS][NSYN_PN][TOT_MS]; // Generalise
} RECORD;

//RECORD *RECSP[NLAYERS][NRECORDS_PL];
//#endif

typedef long long tstep; // Signed to allow spikeTimes[0] = -BIG

// Place in affNeurons structure?
int *** affNeurons_EfE;
int *** affNeurons_ElE;
int *** affNeurons_IE;
int *** affNeurons_EI;
int *** affNeurons_II;

typedef struct {
	int delay; // Actual axonal delay in timesteps
	int size;
	int * queue;
	int count;
	int next;
	int last;
	
	/* Synapse */
	float C, C_tm1;
	float g, g_tm1;
	float delta_g, delta_g_tm1;		/* Constrained to lie in the interval [0,1] */
} AXON;

//typedef struct NEURON; 
typedef struct NEURON{
	//enum NTYPE type;
	NTYPE type;
	int n; // Index ?? l_ind, n_ind
	float V, V_tm1;
	float D, D_tm1; // afferent
	//float C, C_tm1; // efferent
	int lastSpike;
	int spkbin;
	int * spikeTimes;	// Array of spike times
	
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
	
	int rec_flag;
	//bool rec_flag;
	RECORD * rec; 	// Recording structure
} NEURON;

NEURON ** n_E;
NEURON ** n_I;

// Array of arrays of structures?
// Expand to RECORDS[NLAYERS][NTYPES][NRECORDS] where NTYPES:= {0,1} (Excitatory, Inhibitory)

/* Prototypes from spike.c */

extern int spike(char * pfile);
extern NEURON ** allocn(int nlays, int nneurons, NTYPE type);
extern int unallocn(NEURON ** narray, int nlays, int nneurons);
extern void calc_connectivity();
extern void wire_afferents(NEURON * n, int l, int * affNcnx_fE, int * affNcnx_lE, int * affNcnx_I);
extern void alloc_efferents(NEURON * n);
extern void wire_efferents(NEURON * n);
extern void create_axons(NEURON * n);
extern void init_network(int regime);
extern void gen_stimuli(bool rep, int **** trn_stimuli, int **** tst_stimuli, int ** input, int *** shuffle);
extern void calc_input(int loop, int pat, int trans, int *** stimuli, int * input, int ** shuffle, int regime);
extern void update_network(int t, int loop, int input[], int regime);
extern void update_V(int t, float decay_rate, float gLeak, float Vrest, float Thresh, float Vhyper, int inj, NEURON * n);
extern void update_g(NEURON * n, float decay_E, float decay_I, int t);
extern void update_weights(NEURON * n, int t);
extern void update_C(NEURON * n, int t);
extern void update_D(NEURON * n, int t);
extern void init_queue(AXON * a);
extern void enqueue(AXON * a, int x);
extern int dequeue(AXON * a);
extern int next_spike(AXON * a);
extern int isempty(AXON * a);
extern void print_queue(AXON * a);