/*
 *  parameters.h
 *  Spike
 *
 *  Created by Ben Evans on 18/08/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include <stdbool.h>

typedef enum {
	INHIB,
	EXCIT
} NTYPE;

typedef enum {
	NoLearning,
	Learning
} REGIMETYPE;

typedef enum {
	Off,
	Uniform,
	Gaussian
} NOISE;

typedef enum {
	MinD,
	ConstD,
	UniformD,
	GaussD
} DELAY;

/********************** Main Network Parmeters *************************/
// Update read_parameters and defaults.m when a new parameter is added
typedef struct {
	// Simulation
	int loops;
	bool train;// = true;
	bool pretrain;// = true;
	bool trainPause;
	int noise;
	int nRecordsPL;// = 0;
	float TotalTime;
	int TotalMS;
	int TotalTS;
	int TSperMS;
	int spkBuffer;
	
	// Stimuli
	bool random_order;// = true;
	bool localRep;// = true;
	float current;// = 1.25e-9;
	int nStimuli;
	int nTransPS;
	float transP_Train;
	float transP_Test;
	int shift;
	int nFiringNeurons;
	float a; // Sparseness of stimuli
	
	// Network
	int nLayers;// = 2;
	int nWLayers;
	bool inputInhib;
	int nExcit;// = 120;
	int nSynEfE;
	int nSynElE;
	int nSynIE;
	int nInhib;// = 40;
	int nSynEI;
	int nSynII;
	//bool noise;// = false;
	DELAY axonDelay;
	float d_const;
	float d_min;
	float d_max;
	float d_mean;
	float d_sd;
	
	// Cell bodies
	float capE;// = 2.0e-10;
	float capI;// = 10.0e-12;
	float gLeakE;// = 10.0e-9;
	float gLeakI;// = 5.0e-9;
	float VrestE;
	float VrestI;
	float VhyperE;
	float VhyperI;
	float ThreshE;
	float ThreshI;
	float VrevE;
	float VrevI;
	float refract;// = 0.002
	
	// Synapses (afferent axons)
	float alphaC;// = 0.5;
	float tauC;// = 0.005;
	float alphaD;// = 0.5;
	float tauD;// = 0.009;
	float learnR;// = 0.1;
	float tauEE;// = 0.01;
	float Dg_IE;// = 0.5;
	float tauIE;// = 0.001
	float Dg_EI;// = 0.5;
	float tauEI;// = 0.001;
	float Dg_II;// = 0.5;
	float tauII;// = 0.001;
	float gMax;// = 48.0e-10;
} PARAMS;

//PARAMS * mp;
extern PARAMS * mp;

#endif
