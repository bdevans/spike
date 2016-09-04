/*
 *  utils.c
 *  Spike
 *
 *  Created by Ben Evans on 19/02/2009.
 *  Copyright 2009 University of Oxford. All rights reserved.
 *
 */

#include "utils.h"

char errloc[PATHBUFF];
// -> #include <errno.h> for errno and perror(), and <string.h> for strerror(). 
// TSF: int strerror_r(int errnum, char *strerrbuf, size_t buflen);
void exit_error(const char * process, const char * statement)
{
	fprintf(stderr, "%s: %s. Exiting!\n%s (%u)\n", process, statement, strerror(-errno), -errno); // System generated error codes are usually negative??
	perror(statement);
	fflush(stderr);
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


void * myrealloc(void * mem, int n)
{ // N.B. the returned pointer must be assigned to the passed pointer as memory could be moved to a new location
	void * new;
	/*if (!mem)	mem = myalloc(n); // Implicit in realloc */
	
	if (!n) // if size==0 free memory
		mem = myfree(mem);
	else
		if ( (new = realloc(mem, n)) == NULL ) // realloc and check for NULL pointer
		{
			mem = myfree(mem);
			exit_error("MYREALLOC", "NULL pointer from realloc");
		}
		else
		{
			mem = new;
		}
			
	return mem;
}


void * myfree(void * mem)
{
	if (mem) // != NULL
	{
		free(mem);
		mem = NULL;
	}
	return mem;
}


FILE * myfopen(const char * filename, const char * args)
{
	FILE * fp = NULL;
	if ((fp = fopen(filename, args)) == NULL)
		exit_error("MYFOPEN", filename); //return NULL; else
	return fp;
}


void filecopy(FILE * ifp, FILE * ofp)
{
	int c;
	while ((c = getc(ifp)) != EOF)
		putc(c, ofp);
	fflush(ofp);
	return;
}


void append(const char * fname, const char * str)
{
	FILE * fp = myfopen(fname, "a+");
	fprintf(fp, "%s", str);
	fclose(fp);
	return;
}


bool file_exists(const char * filename)
{
	FILE * file;
	file = fopen(filename, "r");
	//if (file = fopen(filename, "r"))
	if(file)
	{
		fclose(file);
		return true;
	}
	return false;
}

/*sprint(char * filename, const char * format);
{
	int slen = 0;
	if (snprintf(filename, FNAMEBUFF, format) >= FNAMEBUFF)
		exit_error("sprintf", "String too long for FNAMEBUFF");
	
	assert(slen < FNAMEBUFF);
	
	if (snprintf(syscmd, BUFSIZ, "mkdir %s", mp->imgDir) >= BUFSIZ)
		fprintf(stderr, "Warning! Undersized buffer: %s", syscmd);
 if (snprintf(dst, sizeof(dst) - 1, "%s", src) > sizeof(dst) - 1) { 
 // Overflow
}
}*/

void dprint(const char * debug_str, const char * args)
{
	fprintf(stderr, debug_str, args);
	fflush(stderr);
	return;
}

/*
void getFilenameNoPath(char * fullname, char * shortname)
{
	char * cptr;
	cptr = strrchr(fullname, '/');
	if (cptr)
		strcpy(shortname, (cptr+1));
	else
		strcpy(shortname, fullname);
}

void getFileExtension(char * fullname, char * ext)
{
	char * cptr;
	cptr = strrchr(fullname, '.');
	if (cptr) // !=NULL
		strcpy(ext, (cptr+1));
	else
		ext[0] = '\0';
}*/

void getFileParts(char * fullname, FILEPARTS * fp)
{
	char * cptr;
	char * homepath;
	unsigned int len = 0, homelen = 0, offset = 0;
	fp->path[0] = fp->fname[0] = fp->fext[0] = '\0';
	assert(fullname);
	
	cptr = strrchr(fullname, '/');
	if (cptr) // != NULL
	{
		strncpy(fp->fname, (cptr+1), strlen(cptr+1)+1);
		if (fullname[0] == '~')
		{
			homelen = strlen(getenv("HOME"));
			homepath = getenv("HOME");
			strncpy(fp->path, homepath, homelen);
			fp->path[homelen] = '\0';
			offset = 1;  // Skip leading '~'
		}
		len = cptr - (fullname+offset); // Skip leading '~'
		strncat(fp->path, fullname+offset, len); // Adds \0 at end
	}
	else
		strncpy(fp->fname, fullname, strlen(fullname)+1);
	
	cptr = strrchr(fp->fname, '.');
	if (cptr)
	{
		strncpy(fp->fext, (cptr+1), strlen(cptr+1)+1);
		fp->fext[strlen(cptr+1)] = '\0'; // Null terminate extension
		*cptr = '\0'; // Null terminate filename
	}
	return;
}


void getTimeString(char * timeStr, size_t buffer, double secs, const char * format)
{
	int hours = (int)(secs/3600); //floor
	secs -= hours*3600;
	int mins = (int)(secs/60);
	secs -= mins*60; //secs = elapsed - (mins*60) - (hours*3600);
	if (!format || strcmp(format, "ms")==0) // Millisecond precision by default
		snprintf(timeStr, buffer, "%d:%02d:%06.3lf",hours,mins,secs);
	else if (strcmp(format, "s")==0) // Seconds
		snprintf(timeStr, buffer, "%d:%02d:%02.0lf",hours,mins,secs);
	else if (strcmp(format, "m")==0) // Minutes
		snprintf(timeStr, buffer, "%d:%02d",hours,mins);
	else
		exit_error("getTimeString", "Unknown format argument");
	
	return;
}
