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
//#include <string.h>

extern void exit_error(const char * process, const char * statement);
extern void * myalloc(int n);
extern FILE * myfopen(const char * filename, const char * args);
extern void filecopy(FILE * ifp, FILE * ofp);
extern bool file_exists(const char * filename);

#endif // _UTILS_H
