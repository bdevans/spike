/*
 *  rng.h
 *  Spike
 *
 *  Created by Ben Evans on 18/08/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _RNG_H
#define _RNG_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

FILE *random_seeds_ptr; 

/************ RNG macros ************/
#ifndef RSFILE
#define RSFILE	"random_seeds.rsd"	/* Default random seed file name */
#endif

#ifndef BUFFER
#define BUFFER	512
#endif

#define MBIG	1000000000
#define MSEED	161803398
#define MZ		0
#define FAC		(1.0/MBIG)

#define IA		16807
#define IM		2147483647
#define AM		(1.0/IM)
#define IQ		127773
#define IR		2836
#define NTAB	32
#define NDIV	(1+(IM-1)/NTAB)
#ifndef EPS
#define EPS		1.2e-7
#endif
#define RNMX	(1.0 - EPS)

#define SEED_RANDS1_NEW		0		/* Seed for rands1_new [+ve int] (0: randomly generated seed.)	*/
#define SEED_RAN3			0		/* Seed for ran3 [-ve int] (0: randomly generated seed.)		*/
#define SEED_GASDEV			0		/* Seed for gasdev [-ve int] (0: randomly generated seed.)		*/

#define CLASS_SIZE			1024	/* From rands1_new. Maximum number of neurons in network|layer.	*/
/* Should be max(NEXCIT,NINHIB) since only one synapse type is sampled from at once -see rands1_new	*/

/* RNG Variables */

long idum;
long idum_gasdev;
int i, iv;

long ans[3]; // array to store values read from file
extern int RERUN;

/* Function Prototypes */

// Functions for uniform and Gaussian random variables & to set the random seeds -copy functions in...
extern float ran3(long *idum);				/* RNG: Uniform disribution								*/
extern float gasdev(long *idum);			/* RNG: Gaussian distribution							*/
extern float ran1(long *idum);				/* Required for gasdev									*/
extern int set_random_seeds(int rerun);				/*														*/
extern int rands1_new(int mn,int mx,int *iv,int mode);/* Draw ints from [mn,mx] without replacement */
extern int rnd_new(int mn, int mx);			/* Call at start of program to initialise RNG			*/
#endif
