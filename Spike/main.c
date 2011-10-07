#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#ifdef SERIAL
#undef _OPENMP
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <pwd.h>

#include "utils.h"
#include "globals.h"
//#include "rng.h"
#include "parameters.h"
#include "read_parameters.h"

extern int spike(PARAMS * mp);

char * RESDIR = RESULTSDIRECTORY;
char * DPFILE = DEFAULTPFILE;
char * MPFILE = OUTPUTPFILE;
char * IMGDIR = IMGDIRECTORY;
char * IPFILE = IMGPARAMFILE;

//int RERUN = 0; // rng.c requires this to be global

// otool -L Spike // to test which dynamic libraries are being used
// http://discussions.apple.com/thread.jspa?threadID=1741520

PARAMS * mp;
unsigned long int seed = 0;
gsl_rng * mSeed = NULL;
gsl_rng ** states = NULL;
int nThreads = 1;

int main (int argc, const char * argv[])
{	
	//int set_error = 0;
	bool rerun = false;
	bool compress = true;
	bool seedFlag = false;
	int hours = 0;
	int mins = 0;
	char *pfile = NULL;
	char *imageArchive = NULL;
	//char *sfile = NULL;
	char *cpfile = "CLIparams.m";
	FILE * cli_FP = NULL;
	char syscmd[BUFSIZ]; // stdio.h : 1024
	int th = 0;
	float proporption = 1.0;
#ifdef _OPENMP
	bool dynamic = false;
#endif
	bool genNewSeed = false;
	bool p_flag = false;	// Parameter (CLI) flag
	bool pf_flag = false;	// Parameter file flag
	bool ia_flag = false;	// Image archive flag
	int pINcount = 0;
	int dcount = 0; //Default parameter count
	int fcount = 0; // File parameter count
	int pcount = 0; // Parameter count
	int err = 0;

	bool recFlag = false;
	int slen = 0;
	char * recList = NULL;
	
	bool skip_arg = false;
	//double wtime, begin;
	char timeStr[FNAMEBUFF];
    int c, result;
	struct tm *ts;
	char hostname[FNAMEBUFF];
	
	char schedule[BUFSIZ];
	strncpy(schedule, "<default>", BUFSIZ);
	
	time_t now = time(NULL);
	ts = localtime(&now);
	strftime(timeStr, FNAMEBUFF, "%a %d/%b/%Y %H:%M:%S", ts);
	err = gethostname(hostname, FNAMEBUFF);
	//char *user = getenv("USER");
	struct passwd *userinfo = getpwuid(geteuid());
	char *user = userinfo->pw_name;
	if (user && err != -1)
		printf("Program started by %s@%s : %s\n", user, hostname, timeStr);
	if (strcmp(user, "nobody")==0)
		SIM.Xgrid = true;
	else
		if (getenv("OMP_SCHEDULE")) // not NULL string
			strncpy(schedule, getenv("OMP_SCHEDULE"), BUFSIZ);
	
	char * rsfile = RSFILE;
	
	//printf("Optimization level: %d", );	
	
	if (argc==1)
	{
		printf("\nSpike usage:\n");
		printf("Checking for \"%s\" in current directory... \t\t\t    [%s]\n",DPFILE,\
			   (file_exists(DPFILE))?"OK":"NO");
		/*if (file_exists(DPFILE))
			printf("[OK]\n");
		else
		{
			printf("[MISSING]\n");
			printf("This may be symlinked e.g. \"ln -s ~/model/code/%s .\"\n",DPFILE);
		}*/
		printf("Checking for \"%s\" in current directory... \t\t    [%s]\n",rsfile,\
			   (file_exists(rsfile))?"OK":"NO");
		printf("-f <filename>\t: Pass parameter filename\n");
		printf("-r[erun]\t: Rerun simulation with the random seed in %s\n",rsfile);
		printf("-g[enerate]\t: Generate new random seed in %s and exit\n",rsfile);
		printf("-s <seed>\t: Explicitly pass random seed [0, (2^32)-1]\n");
		printf("-k <record list>: Pass list of neurons to be recorded\n");
		printf("-d[ynamic]\t: Set number of threads to be dynamic\n");
		printf("-m <proportion>\t: Set number of threads to a proportion of cores [0.0, 1.0]\n");
		printf("-t <threads>\t: Explicitly set the number of threads to use\n");
		printf("-p <parameter>\t: Pass a parameter string <name>=<value>\n");
		//printf("-i[mage] <directory>\t: Pass directory of filtered images [***incomplete***]\n");
		printf("-j <images>.tbz\t: Pass compressed image archive\n");
		printf("-u[ncompressed]\t: Prevent data compression\n");
		printf("-x[grid]\t: Set as an Xgrid simulation i.e. print progress information\n");
		printf("--------------------------------------------------------------------------------\n");
		return 0;
	}
	
	while (--argc > 0 && (*++argv)[0] == '-')
	{
		//skip_arg = 0;
		while (!skip_arg && (c = *++argv[0]))
		{
			switch (c)
			{
				case 'f':	// Parameter file name
					pf_flag = true;
					pfile = myalloc(strlen(*++argv)+1); //sizeof(char)==1 guaranteed
					strcpy(pfile, *argv);
					skip_arg = true;
					argc--;
					break;
				
				case 'r':	// Rerun with last seed
					rerun = true; //RERUN = 1;
					break;
					
				case 'g':	// Generate a new random seed and exit
					genNewSeed = true;
					break;
					
				case 's':	// Explicitly pass random seed (takes precedence over -r)
					seed = atol(*++argv);
					seedFlag = true;
					skip_arg = true;
					argc--;
					break;
					
				case 'k':	// Read in list of neurons
					//int * recordSet = NULL;
					recFlag = true;
					slen = strlen(*++argv);
					recList = myalloc(slen+1);
					//char * tstr = NULL;
					//int count = 0;
					strncpy(recList, *argv, slen);
					recList[slen] = '\0'; // NULL terminate last byte
					//strtok(list, ";");
					//while (tstr && (*trim(tstr) != ('[' || ']')) && (tstr != '\0'))
					//	count = parseIntVector(list, &recordSet);
					skip_arg = true;
					argc--;
					break;
					
				case 'd':	// May reduce num_thd depending on system load
#ifdef _OPENMP
					dynamic = true;
					omp_set_dynamic(dynamic); // Can not be used with RNG
					// if nCores > 4
					// if nThreads > nCores-1 -> set nThreads = nCores - 2...
#else
					fprintf(stderr, "OpenMP disabled\n");
#endif
					break;
					
				case 'm':	// Set the proportion of threads from the CLI [0.0, 1.0]
#ifdef _OPENMP
					proporption = atof(*++argv);
					assert(proporption > 0.0 && proporption <= 1.0);
					omp_set_num_threads(round(omp_get_num_procs()*proporption));
#else
					fprintf(stderr, "OpenMP disabled\n");
#endif
					skip_arg = true;
					argc--;
					break;
				
				case 't':	// Set the number of threads from the CLI
#ifdef _OPENMP
					nThreads = atoi(*++argv);
					/*if (nThreads >= omp_get_num_procs())
						omp_set_dynamic(true);
					else
						omp_set_num_threads(nThreads);*/
					omp_set_num_threads(nThreads);
					if (nThreads >= omp_get_num_procs())
						printf("Warning: nThreads (%d) >= nProcessors (%d)!\n",\
							   nThreads, omp_get_num_procs());
#else
					fprintf(stderr, "OpenMP disabled\n");
#endif				
					skip_arg = true;
					argc--;
					break;
					
				case 'p':	// Code to pass a parameter string e.g. "param=0"
					if (!p_flag)
					{
						cli_FP = myfopen(cpfile, "w");
						p_flag = true;
					}
					fprintf(cli_FP, "%s;\n", *++argv);
					skip_arg = true;
					argc--;
					pINcount++;
					break;
					
				case 'i':	// Pass Image directory
					
					break;
					
				case 'j':	// Alternatively pass compressed tar (cjvf) of images
					
					ia_flag = true;
					slen = strlen(*++argv);
					imageArchive = myalloc(slen+1);
					strncpy(imageArchive, *argv, slen);
					imageArchive[slen] = '\0'; // NULL terminate last byte
					skip_arg = true;
					argc--;
					break;
					
				case 'u':	// Keep data uncompressed
					compress = false;
					break;
					
				case 'x':	// Xgrid simulation
					SIM.Xgrid = true;
					break;
				
				default:
					printf("Illegal arguement: %c\n", c);
					argc = 0;
					break;
			}
			if (skip_arg)
			{
				skip_arg = false;
				break;
			}
		}
	}
	
#ifdef _OPENMP
#pragma omp parallel //private (th_id)
	{
		//th_id = omp_get_thread_num();
		nThreads = omp_get_num_threads(); //num_thd
#pragma omp single
		{
			printf("Threads: %d/%d\t{OMP_DYNAMIC=%s, OMP_NESTED=%s, OMP_SCHEDULE=%s}\n", \
				   nThreads, omp_get_num_procs(), \
				   (omp_get_dynamic() ? "TRUE" : "FALSE"), \
				   (omp_get_nested() ? "TRUE" : "FALSE"), \
				   (schedule));
		}
	}
#else
	nThreads = 1;
	printf("Executing in serial.\n");
#endif
	
	FILE * randSeedFP;
	char * sString, buffer[BUFSIZ];
	unsigned long long seedMod = pow(2, 32); // 32-bit unsigned seeds (max value = pow(2, 32)-1)	
	
	if (genNewSeed) // Generate a new random seed file and exit
	{
		//set_random_seeds(0);
		// Move rng routine up here?
		seed = (unsigned long) time((time_t *) NULL);
		seed %= seedMod;
		printf("Warning: Creating new seed in %s: %ld\n", rsfile, seed);
		randSeedFP = myfopen(rsfile, "w");
		fprintf(randSeedFP, "mSeed: \t%ld\n", seed);
		fclose(randSeedFP);
		// Error checking... if (file_exists())
		printf("New seeds generated!\n"); // print location and filename
		return 0;
	}
	
	// Allocate and initialise model parameters structure
	mp = myalloc(sizeof(*mp)); // Place in getParameters with default init?
	mp->imgList = NULL;
	mp->LvExcit = mp->LvInhib = mp->LpEfE = mp->LpElE = mp->LpEI = mp->LpIE = mp->LpII = 0;
	mp->vExcit = mp->vInhib = mp->vScales = mp->vOrients = mp->vPhases = NULL;
	mp->pCnxEfE = mp->pCnxElE = mp->pCnxIE = mp->pCnxEI = mp->pCnxII = NULL;
	mp->layDim = NULL; 
	mp->vSquare = mp->vRecords = NULL;
	mp->rInhib = 0.0;
	
	int syserr = 0;
	
	if (ia_flag) //assert(mp->useFilteredImages);
	{
		FILEPARTS * fp = myalloc(sizeof(*fp));
		getFileParts(imageArchive, fp);
		slen = strlen(fp->fname);
		mp->imgDir = myalloc(slen+1);
		strncpy(mp->imgDir, fp->fname, slen);
		mp->imgDir[slen] = '\0'; // NULL terminate last byte
		assert(file_exists(imageArchive));
		if (snprintf(syscmd, BUFSIZ, "mkdir %s", mp->imgDir) >= BUFSIZ)
			fprintf(stderr, "Warning! Undersized buffer: %s", syscmd);
		if((syserr = system(syscmd))==0)
		{
			printf("Now extracting %s...\t", imageArchive);
			if(snprintf(syscmd, BUFSIZ, "tar -xf %s -C %s/",imageArchive, mp->imgDir) >= BUFSIZ)
				fprintf(stderr, "Warning! Undersized buffer: %s", syscmd);
			syserr = system(syscmd);
			if (syserr)
				exit_error("main.c", "Error extracting image archive");
			else
				printf("Images successfully extracted to %s\n", mp->imgDir);
		}
	}
	
	// Read in parameters from .m file
	printf("Reading parameters file: \"%s\"", !pfile ? DPFILE : pfile);
	if (p_flag)
		fclose(cli_FP);
	dcount = read_parameters(mp, DPFILE);
	fcount = (pfile != NULL) ? read_parameters(mp, pfile) : 0;
	pcount = (p_flag) ? read_parameters(mp, cpfile) : 0;
	
	if (!mp->useFilteredImages)
		assert(pcount == pINcount);
	printf(" {%d,%d,%d}\tParsing complete!\n", dcount, fcount, pcount);
	
	pcount = printParameters(mp, MPFILE); // Variables to read into Matlab
	printf("%d parameters written to %s\n", pcount, MPFILE);
	
	const gsl_rng_type * T = gsl_rng_default; // Set RNG type
	gsl_rng_env_setup(); // http://www.gnu.org/software/gsl/manual/html_node/Random-number-environment-variables.html
	mSeed = gsl_rng_alloc(T); // Used for serial sections with randomness
	
	if (!seedFlag)
	{
		if (rerun) // Also rerun with parameters.m?
		{
			randSeedFP = myfopen(rsfile, "r");
			if ((sString = fgets(buffer, sizeof(buffer), randSeedFP)) != NULL) //while
				seed = atol(strrchr(sString,':')+1); //ans[count++]
			fclose(randSeedFP); 
			printf("Rerunning simulation with %s: %ld (%s)\n", rsfile, seed, gsl_rng_name(mSeed));
		}
		else
		{
			seed = (unsigned long) time((time_t *) NULL);
			seed %= seedMod;
			printf("Warning: Creating new seed in %s: %ld (%s)\n", rsfile, seed, gsl_rng_name(mSeed));
			randSeedFP = myfopen(rsfile, "w");
			fprintf(randSeedFP, "mSeed: \t%ld\n", seed);
			fclose(randSeedFP);
		}
	}
	gsl_rng_set(mSeed, seed); //gsl_rng_set(mSeed, -idum);
	//printf("Random number master seed is %ld\n", seed);
	//printf("GSL generator type: %s\n", gsl_rng_name(mSeed));
	
	if (mp->noise) // Could add a state to every neuron to achieve same results with different threads
	{
#ifdef _OPENMP
#pragma omp parallel
		{
#pragma omp single
			{
				omp_set_dynamic(false); // Do not adjust number of threads according to system load
			}
		}
#endif
		states = myalloc(nThreads * sizeof(**states));
		for (th=0; th<nThreads; th++)
		{
			states[th] = gsl_rng_alloc(T);
			gsl_rng_set(states[th], seed+th+1);
		}
	}
	
	/*if (recFlag)
	{
		int ** recordSet = myalloc(mp->nLayers * sizeof(*recordSet));
		char * tstr = NULL;
		int count = 0;
		strtok(list, ";");
		while (tstr && (*trim(tstr) != ']') && (tstr != '\0'))
		for (l=0; l<mp->nLayers; l++)
		{
			mp->vRecords[l] = parseIntVector(list, &recordSet[l]);
			tstr = strtok(NULL, ";");
		}
			
		
		myfree(recList);
	}*/
	
#if DEBUG > 0
	printf("*** Executing with Debug level %d ***\t", DEBUG);
#endif
#ifdef NDEBUG
	printf("Warning: Executing without error checking!\t");
#endif
	fprintf(stdout, "DT = %f ms\n", mp->DT*1000);
	// Display dynamic libraries: otool -L ~/bin/SpikeNet/Debug/Spike
	
#ifdef _OPENMP // Use omp function omp_get_wtime
	double begin = omp_get_wtime();
#else
	time_t start = time(NULL);
#endif
	
	if (!SIM.Xgrid)	// Remove *.dat and *.tbz // Add this to arg switch (-c)?
		system("rm *.dat"); //system("rm *.dat *.tbz");
	
	result = spike(mp);
	
	// Compress data files for crash-free xgrid! '-j' Uses bzip (*.tbz equivalent to *.tar.bz2)
	// Append files to fileList and call system(syscmd); once and keep fileList
	//snprintf(syscmd, BUFSIZ, "tar -cjvf %s.tbz %s > fileList","connectivity","*affNeurons.dat");
	

	if (!SIM.Xgrid) // /sbin/md5
		system("md5 *.dat > datHashs");
	/*snprintf(syscmd, BUFSIZ, "xargs rm < fileList");*/
	//--remove-files (remove files after adding them to the archive) : only 10.5
	// Check that system() returned 0 (no errors) Bash: echo $?

#pragma omp parallel sections private(syserr) // Experimental!
	{
#pragma omp section
	{
	if (compress)
	{
		printf("\tCompressing data to .tbz archives...\t");
		fflush(stdout);
		
		if (!mp->useFilteredImages)
			if ((syserr = system("tar -cjf stimuli.tbz *stimuli.dat")) == 0)
				system("tar -tf stimuli.tbz | xargs rm");
		
		//#pragma omp section
		if(mp->nRecordsPL)
			if ((syserr = system("tar -cjf records.tbz R*.dat")) == 0)
				system("tar -tf records.tbz | xargs rm");
		
		//#pragma omp section		
		if (mp->printConnections)
		{
			if (mp->SOM)
				syserr = system("tar -cjf connectivity.tbz *affNeurons*.dat *affDelays*.dat *dist*.dat");
			else
				syserr = system("tar -cjf connectivity.tbz *affNeurons*.dat"); //system("tar --remove-files -cjvf connectivity.tbz *affNeurons*.dat > fileList");
			if (!syserr)
				system("tar -tf connectivity.tbz | xargs rm");
		}
		
		//#pragma omp section
		if (mp->pretrain)
			if ((syserr = system("tar -cjf preTraining.tbz pt*.dat")) == 0)
				system("tar -tf preTraining.tbz | xargs rm");
		
		//#pragma omp section
		if (mp->train)
			if ((syserr = system("tar -cjf training.tbz E*.dat")) == 0) // 2> tar_err
				system("tar -tf training.tbz | xargs rm");
		
		//#pragma omp section
		if ((syserr = system("tar -cjf postTraining.tbz L*Spikes.dat L*weights*.dat")) == 0)
			system("tar -tf postTraining.tbz | xargs rm");
		//system(syscmd);
		//system("rm fileList");
		
		printf("Data Compressed!\n");
	}
	//#pragma omp section
	if (!SIM.Xgrid) // Print md5 #'s // /sbin/md5
	{
		//system("md5 Spike");
		system("md5 parameters.m");
		system("md5 datHashs");
		//system("md5 *.tbz"); // Contains metadata (e.g. timestamps) which will give different #s
	}
	} // End of section
		
	// Clean up
#pragma omp section
	if (pf_flag)
		myfree(pfile);
	
#pragma omp section
	if (ia_flag)
	{
		myfree(imageArchive);
		if(snprintf(syscmd, BUFSIZ, "rm -R %s/", mp->imgDir) >= BUFSIZ)
			fprintf(stderr, "Warning! Undersized buffer: %s", syscmd);
		system(syscmd);	// Delete expand image files
	}
	
#pragma omp section
	{
	gsl_rng_free(mSeed);
	if (mp->noise)
	{
		for (th=0; th<nThreads; th++)
			gsl_rng_free(states[th]); // Free all memory associated with generator
		myfree(states);
	}
	}
	
	//if (recFlag)	// Free list of records
} // End of parallel sections
	
	if (mp->useFilteredImages)
	{
		myfree(mp->imgDir);
		myfree(mp->imgList);
		myfree(mp->vScales);
		myfree(mp->vOrients);
		myfree(mp->vPhases);
	}
	
	myfree(mp->vRecords);
	myfree(mp->vExcit);
	myfree(mp->vInhib);
	myfree(mp->pCnxEfE);
	myfree(mp->pCnxElE);
	myfree(mp->pCnxIE);
	myfree(mp->pCnxEI);
	myfree(mp->pCnxII);
	myfree(mp->layDim);
	myfree(mp->vSquare);
		
	myfree(mp);
	
#ifdef _OPENMP
	//getTimeString(timeStr, FNAMEBUFF, omp_get_wtime()-begin);
	double wtime = omp_get_wtime() - begin;
	//double integral;
	//double fraction = modf(wtime, &integral);
	//duration = (time_t) round(integral);
	hours = floor(wtime/3600);
	wtime -= hours*3600;
	mins = floor(wtime/60);
	wtime -= mins*60; //secs = wtime - (mins*60) - (hours*3600);
	snprintf(timeStr, FNAMEBUFF, "%d:%02d:%05.2lf (%d Threads)",\
			 hours,mins,wtime,nThreads);
#else
	time_t duration = time(NULL) - start; //	finish = round(time(NULL) - start);
	hours = floor(duration/3600);
	duration -= hours*3600;
	mins = floor(duration/60);
	duration -= mins*60;
	int secs = duration;
	snprintf(timeStr, FNAMEBUFF, "%d:%02d:%02d (Serial)",hours,mins,secs);
#endif
	
	if (result==0)
		printf("Simulation completed in %s!\n",timeStr);
	else
	{
		printf("Simulation aborted after %s!\n",timeStr);
		return 1;
	}
	printf("--------------------------------------------------------------------------------\n");
	
    return 0;
}
