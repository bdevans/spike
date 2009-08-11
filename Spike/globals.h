/*
 *  globals.h
 *  Spike
 *
 *  Created by Ben Evans on 09/12/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

/*** Prevent recursive includes ***/
#ifndef _GLOBALS_H
#define _GLOBALS_H
/**********************************/

#define DT		0.00001
#define EPS		1.2e-7
#define BIG		999999

/* Debugging output */
#ifndef DEBUG
#define DEBUG	0	/* Debug output level (defined as 5 in the Xcode Debug build)						*/
#endif
#if DEBUG == 0
#define NDEBUG		/* Define to turn off assert statements	- must come before #include <assert.h>		*/
#endif

/* Default parameter file names */
#define RESULTSDIRECTORY	"Results"	/* Output file directory name */
extern char * RESDIR;
#define DEFAULTPFILE	"defaults.m" //"DEFAULTS.prm"	/* Default input parameter file name */
extern char * DPFILE;
#define OUTPUTPFILE		"parameters.m"	/* Default output parameter file name */
extern char * MPFILE;
/* Define a diagnostic parameter file to run for debugging */

#define FNAMEBUFF	32

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
