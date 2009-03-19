/*
 *  utils.c
 *  Spike
 *
 *  Created by Ben Evans on 19/02/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "utils.h"

// -> #include <errno.h> for errno and perror(), and <string.h> for strerror(). 
void exit_error(const char * process, const char * statement)
{
	fprintf(stderr, "%s: %s. Exiting!\n", process, statement);
	exit(EXIT_FAILURE);
}

void * myalloc(int n)
{
	void * mem;
	if (n > 0)
	{
		if ((mem = malloc((unsigned) n)) == NULL)
			exit_error("MYALLOC", "NULL pointer from malloc");
	}
	else
		mem = NULL;
	
	return mem;
}

/*FILE * myfopen(const char * filename, const char * args)
{
	FILE * fp;
	fp = fopen(filename, args);
	if (fp == NULL) // exit_error("MYFOPEN", "NULL pointer from fopen");
	{
		fprintf(stderr, "MYFOPEN ERROR: Could not open file %s.\n", filename);
		exit (EXIT_FAILURE);
	}
	else
		return fp;
}*/

void filecopy(FILE * ifp, FILE * ofp)
{
	int c;
	while ((c = getc(ifp)) != EOF)
		putc(c, ofp);
	fflush(ofp);
	return;
}

/*bool file_exists(const char * filename)
{
	FILE * file;
	if (file = fopen(filename, "r"))
	{
		fclose(file);
		return true;
	}
	return false;
}*/