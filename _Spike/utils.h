/*
 *  utils.h
 *  Spike
 *
 *  Created by Ben Evans on 19/02/2009.
 *  Copyright 2009 University of Oxford. All rights reserved.
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
//#include <stdarg.h> // For passing variable number or arguments

#define PATHBUFF	1024

extern void exit_error(const char * process, const char * statement);
extern void * myalloc(int n);
extern void * myrealloc(void * mem, int n);
extern void * myfree(void * mem);
extern FILE * myfopen(const char * filename, const char * args);
extern void filecopy(FILE * ifp, FILE * ofp);
extern void append(const char * fname, const char * str);
extern bool file_exists(const char * filename);
extern void dprint(const char * debug_str, const char * args);

typedef struct {
	char path[PATHBUFF];
	char fname[PATHBUFF];
	char fext[PATHBUFF];
}FILEPARTS;

extern void getFileParts(char * fullname, FILEPARTS * fp);
extern void getTimeString(char * timeStr, size_t buffer, double secs, const char * format);

extern char errloc[PATHBUFF]; 
#define EE(errmsg)   do { \
                    snprintf(errloc, PATHBUFF, "\n[%s (%u) : %s] --> ", \
                        strrchr(__FILE__, '/')+1, __LINE__, __FUNCTION__); \
                    exit_error(errloc,errmsg); } while(0)
// Replace ;'s with ,'s so it can be used as EE();? c.f. assert.h
//const char *buildString = "Compiled at "__DATE__", "__TIME__".";

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

#if !defined(SWAP) // XOR swap - ints only. 
#define SWAP(a, b) (((a) == (b)) || (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b))))
#endif

#endif // _UTILS_H

// The snprintf() and vsnprintf() functions will write at most n-1 of the characters printed into the output string (the n'th character then gets the terminating `\0'); if the return value is greater than or equal to the n argument, the string was too short and some of the printed characters were discarded.  The output is always null-terminated.
// The strncpy() function copies at most n characters from s2 into s1. If s2 is less than n characters long, the remainder of s1 is filled with `\0' characters. Otherwise, s1 is not terminated.
// If src contains n or more characters, strncat() writes  n+1  characters to  dest  (n  from src plus the terminating null byte).  Therefore, the size of dest must be at least strlen(dest)+n+1
