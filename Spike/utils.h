/*
 *  utils.h
 *  Spike
 *
 *  Created by Ben Evans on 19/02/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


/*** Prevent recursive includes ***/
#ifndef _UTILS_H
#define _UTILS_H
/**********************************/


#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#define PATHBUFF	1024

extern void exit_error(const char * process, const char * statement);
extern void * myalloc(int n);
extern void * myrealloc(void * mem, int n);
extern void * myfree(void * mem);
extern FILE * myfopen(const char * filename, const char * args);
extern void filecopy(FILE * ifp, FILE * ofp);
extern bool file_exists(const char * filename);
extern void dprint(const char * debug_str);

typedef struct {
	char path[PATHBUFF];
	char fname[PATHBUFF];
	char fext[PATHBUFF];
}FILEPARTS;

extern void getFileParts(char * fullname, FILEPARTS * fp);

#define PRINT_INT(token) printf(#token " = %d;\n", token);
#define PRINT_FLOAT(token) printf(#token " = %G;\n", token);
#define FPRINT_INT(fptr, token) fprintf(fptr, #token " = %d;\n", token);
#define FPRINT_FLOAT(fptr, token) fprintf(fptr, #token " = %G;\n", token);
#define FPRINT_STRING(fptr, token) fprintf(fptr, #token " = '%s';\n", token);

// http://www.xgc.com/manuals/gcc-erc32-ug/p2node29.html
/* #define max(a,b) \
 ({typedef _ta = (a), _tb = (b);  \
 _ta _a = (a); _tb _b = (b);     \
 _a > _b ? _a : _b; })
 */
#if !defined(MAX) 
#define MAX(A, B) ((A) > (B) ? (A) : (B)) 
#endif 

#if !defined(MIN) 
#define MIN(A, B) ((A) < (B) ? (A) : (B)) 
#endif 

#endif // _UTILS_H
