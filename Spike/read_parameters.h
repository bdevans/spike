/*
 *  read_parameters.h
 *  Spike
 *
 *  Created by Ben Evans on 28/11/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

/*
 * parse: parse simple name/value pairs
 *
 * SAMPLE BUILD:
 * cc -g -Wall -o parse parse.c
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>

/*
 #if DEBUG == 0
 #define NDEBUG	// Define to turn off debugging with assert statements
 #endif
 */
#include "parameters.h"
//#include "spike.h"
#include "globals.h"

#define MAXLEN 256
#define BUFFER 512 // Make larger for comments or break up comments

//extern DPFILE;
//#define CONFIG_FILE "params.ini"

//#include <stdbool.h>

char * trim(char * s);
int parse_config(PARAMS * params, char * paramfile);
int read_parameters(PARAMS * params, char * pfile);