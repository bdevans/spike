/*
 *  globals.h
 *  Spike
 *
 *  Created by Ben Evans on 09/12/2008.
 *  Copyright 2008 University of Oxford. All rights reserved.
 *
 */

/*** Prevent recursive includes ***/
#ifndef _GLOBALS_H
#define _GLOBALS_H
/**********************************/

#define EPS		1.2e-7
#define BIG		999999

/* Debugging output */
#ifndef DEBUG
#define DEBUG	1	/* Debug output level (set to 3 in the Xcode Debug TARGETS build settings) */
#endif
#ifndef __OPTIMIZE__
#define __OPTIMIZE__ -1
#endif
#if DEBUG == 0
#define NDEBUG		/* Define to disable assert statements - must preceed #include <assert.h> */
#endif

/* Default parameter file names */
#define RESULTSDIRECTORY	"Results"	/* Output file directory name */
extern char * RESDIR;
#define DEFAULTPFILE	"defaults.m" 	/* Default input parameter file name */
extern char * DPFILE;
#define OUTPUTPFILE		"parameters.m"	/* Default output parameter file name */
extern char * MPFILE;
#define IMGDIRECTORY	"Images"
extern char * IMGDIR;
#define IMGPARAMFILE	"imageParams.m"
extern char * IPFILE;
//#define STIMULIFILE		"groups.stm"
extern char STFILE[BUFSIZ];
//#define PPSTIMULIFILE	"ElE_groups.stm"
extern char PPSTFILE[BUFSIZ];
#ifndef RSFILE
#define RSFILE	"random_seeds.rnd"	/* Default random seed file name */
#endif
extern char * rsfile; // = RSFILE;
//extern char * pFile;
/* Define a diagnostic parameter file to run for debugging */

#define DIRBUFF		128
#define FNAMEBUFF	256
//#define ARRBUFF		4096

#endif // _FILE_H

/* Notes */
/*
 %d:	Decimal Integer
 %f:	Float			e.g. %8.2f: 8 spaces wide with 2 d.p. Use '-' before e.g. '8' to right justify
 %c:	Character
 %s:	String

 arr[0] == *arr
 arr[n] == *(arr + n)

 Zero spike bins: Try calloc or double Array[size]={0} for the first time, memset(Array, 0, size), bzero.
 May be quicker just to loop.

 memcpy(DestArray, SourceArray, sizeof(DestArray)); // Quick way of copying contents of an array

 Nifty tricks
 ============

 #1 // Allocate memory or write files without keeping track of variable type
 fwrite(array, sizeof(array), arraydim, ArrayFile);

 #2
 if(!(ptr = malloc(...)))
	fprintf(stderr, "Error: Out of memory!");

 #3 // Allocate memory and test allocation (or fopen)
 if ( ( a = (int *)malloc(i*sizeof(int)) ) == NULL )
 {
	printf("\nError, memory not allocated.\n");
	exit(1);
 }
 */
