#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>

#include "utils.h"
#include "globals.h"
#include "rng.h"
#include "parameters.h"
#include "read_parameters.h"

//#include "spike.h"
extern int spike(PARAMS * mp);
//extern void set_random_seeds();

char * RESDIR = RESULTSDIRECTORY;
char * DPFILE = DEFAULTPFILE;
char * MPFILE = OUTPUTPFILE;

int RERUN = 0; // rng.c requires this to be global

PARAMS * mp;

int main (int argc, const char * argv[])
{	
	int set_error = 0;
	int hours = 0;
	int mins = 0;
	int secs = 0;
	int nThreads = 0;
	char *plural = "";
	//char *scheduling = "FALSE";
	char *dthreads = "FALSE";
	char *nested = "FALSE";
	char *pfile = NULL;
	char *sfile = NULL;
	char *cpfile = "CLIparams.m";
	//char *tmpfile = "tmp";
	FILE * cli_FP;
	FILE * params_FP;
	char rfpath[FNAMEBUFF]; // Results file path
	FILE * parameters_ptr;
	//int NRECORDS;
	bool dynamic = false;

	bool p_flag = false;	// Parameter (CLI) flag
	bool pf_flag = false;	// Parameter file flag
	int pINcount = 0;
	int dcount = 0; //Default parameter count
	int fcount = 0; // File parameter count
	int pcount = 0; // Parameter count
	
	time_t start;
	time_t finish;
	double begin, end;
    int c, result;
	//int th_id;

	int num_thd = 0;
	int skip_arg = 0;
	
	while (--argc > 0 && (*++argv)[0] == '-')
	{
		//skip_arg = 0;
		while (!skip_arg && (c = *++argv[0]))
		{
			switch (c)
			{
				case 'f':	// Parameter file name
					pf_flag = true;
					pfile = myalloc(strlen(*++argv)+1);
					strcpy(pfile, *argv);
					skip_arg = 1;
					argc--;
					break;
				
				case 'r':	// Rerun with last seed
					RERUN = 1;
					break;
					
				case 's':	// Random seed file name - Finish this
					sfile = myalloc(strlen(*++argv)+1);
					strcpy(sfile, *argv);
					skip_arg = 1;
					argc--;
					break;
					
				case 'k':	// Read in list of neurons
					// Dynamically alloc memory for array of structures
					//NRECORDS = *++argv[0];
					//argc--;
					break;
					
				case 'd':	// May reduce num_thd depending on system load
					dynamic = true;
					omp_set_dynamic(dynamic);
					break;
				
				case 't':	// Set the number of threads from the CLI
					nThreads = atoi(*++argv);
					if (nThreads >= omp_get_num_procs())
						omp_set_dynamic(true);
					else
						omp_set_num_threads(nThreads);
					skip_arg = 1;
					argc--;
					break;
					
				case 'p':	// Code to pass a parameter string e.g. "param=0"
					if (!p_flag)
					{
						cli_FP = myfopen(cpfile, "w");
						p_flag = true;
					}
					fprintf(cli_FP, "%s;\n", *++argv);
					skip_arg = 1;
					argc--;
					pINcount++;
					break;
				
				default:
					printf("Illegal arguement: %c\n", c);
					argc = 0;
					break;
			}
			if (skip_arg)
			{
				skip_arg = 0;
				break;
			}
		}
	}
	
#pragma omp parallel //private (th_id)
	{
		//th_id = omp_get_thread_num();
		num_thd = omp_get_num_threads();
#pragma omp single
		{
			plural = ((num_thd==1) ? "" : "s");
			printf("\n*** Program executing with %d/%d thread%s. ***\n", \
				   num_thd, omp_get_num_procs(), plural);
			printf("Dynamic threads (OMP_DYNAMIC): %s\n", \
				   dthreads = (omp_get_dynamic() ? "TRUE" : "FALSE"));
			//printf("Scheduling: %s"); Scheduling is set per loop
			printf("Nested parallelism (OMP_NESTED): %s\n", \
				   nested = (omp_get_nested() ? "TRUE" : "FALSE"));
		}
	}
	
#if DEBUG > 0
	printf("*** Executing with Debug level %d ***\n", DEBUG);
#endif
#ifdef NDEBUG
	printf("Warning: Executing without error checking!\n");
#endif
	
	// Read in parameters from .prm file
	printf("Reading parameters file: \"%s\"", !pfile ? DPFILE : pfile);
	if (p_flag)
		fclose(cli_FP);
	mp = myalloc(sizeof(*mp));
	dcount = read_parameters(mp, DPFILE); // Count defaults?
	fcount = (pfile != NULL) ? read_parameters(mp, pfile) : 0;
	pcount = (p_flag) ? read_parameters(mp, cpfile) : 0;
	assert(pcount == pINcount);
	printf(" {%d},{%d},{%d}\tParsing complete!\n",dcount, fcount, pcount);
	
	if (pfile != NULL)
		sprintf(rfpath, "%s/%s/", RESDIR, pfile); // strtok(pfile, "."));
	else
		sprintf(rfpath, "%s/", RESDIR);
	
	parameters_ptr = myfopen(MPFILE, "w"); // Variables to read into Matlab
	fprintf(parameters_ptr, "DT=%f;\n",DT);
	fprintf(parameters_ptr, "TotalMS=%d;\n",mp->TotalMS);
	fprintf(parameters_ptr, "rfpath='%s';\n\n",rfpath);
	if (pfile != NULL)
	{
		params_FP = myfopen(pfile, "r");
		filecopy(params_FP, parameters_ptr);
		fclose(params_FP);
	}
	
	if (p_flag)
	{
		cli_FP = myfopen(cpfile, "r");
		fprintf(parameters_ptr, "\n");
		filecopy(cli_FP, parameters_ptr);
		fclose(cli_FP);
	}
	fclose(parameters_ptr);

	printf("Now setting random seeds\n");
	set_error = set_random_seeds(RERUN);
	if (set_error)
	{
		RERUN = 0;
		fprintf(stderr,  "WARNING: Creating new seeds...\n");
		set_random_seeds(RERUN);
	}
	start = time(NULL); // Use omp function omp_get_wtime
	begin = omp_get_wtime();
	
	result = spike(mp);
	
	free(mp);
	
	end = omp_get_wtime();
	printf("Total wall time = %lf\n", end - begin);
	finish = round(time(NULL) - start);
	hours = floor(finish/3600);
	finish -= hours*3600;
	mins = floor(finish/60);
	finish -= mins*60;
	secs = finish;
		
	if (result==0)
		printf("Simulation completed in %d:%02d:%02d!\n",hours,mins,secs);
	else
	{
		printf("Simulation aborted at %d:%02d:%02d!\n",hours,mins,secs);
		return 1;
	}
	
    return 0;
}
