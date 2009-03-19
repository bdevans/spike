#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>

#include "rng.h"
//#include "parameters.h"
//#include "spike.h"
#include "utils.h"
#include "globals.h"

char * RESDIR = RESULTSDIRECTORY;
char * DPFILE = DEFAULTPFILE;
char * MPFILE = OUTPUTPFILE;
//extern void set_random_seeds();
extern int spike(char * pfile);

int hours = 0;
int mins = 0;
int secs = 0;
char *plural = "";
//char *scheduling = "FALSE";
char *dthreads = "FALSE";
char *nested = "FALSE";
char *pfile = NULL;
char *sfile = NULL;
int RERUN;
int NRECORDS;
bool dynamic = false;

int main (int argc, const char * argv[])
{
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
				case 'f': // Parameter file name
					pfile = myalloc(strlen(*++argv)+1);
					strcpy(pfile, *argv);
					skip_arg = 1;
					argc--;
					break;
				
				case 'r':
					RERUN = 1;
					// Code to test next argument does not have a '-'

					break;
					
				case 's':
					// Finish this
					sfile = myalloc(strlen(*++argv)+1);
					strcpy(sfile, *argv);
					skip_arg = 1;
					argc--;
					// Read in list of neurons
					// Dynamically alloc memory for array of structures
					//NRECORDS = *++argv[0];
					//argc--;
					break;
					
				case 'd': // May reduce num_thd depending on system load
					dynamic = true;
					omp_set_dynamic(dynamic);
					break;
				/*	
				case 'ALPHA_C':
					ALPHA_C = *++argv[0];
					argc--;
					break;
				*/	
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
	
    printf("Now setting random seeds\n");

	set_random_seeds(RERUN);
	start = time(NULL); // Use omp function omp_get_wtime
	begin = omp_get_wtime();
	
	/* Load filename if passed, otherwise load the default file */
	result = spike(pfile);
	
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