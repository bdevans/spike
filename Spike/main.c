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
#include <dlfcn.h>

#include "utils.h"
#include "globals.h"
#include "parameters.h"
#include "read_parameters.h"

extern int spike(PARAMS * mp);

char * RESDIR = RESULTSDIRECTORY;
char * DPFILE = DEFAULTPFILE;
char * MPFILE = OUTPUTPFILE;
char * IMGDIR = IMGDIRECTORY;
char * IPFILE = IMGPARAMFILE;
char STFILE[BUFSIZ] = ""; // = STIMULIFILE; //char * STFILE = NULL;
char PPSTFILE[BUFSIZ] = "";

#define GSLDIR  "/opt/local/lib/" // Unnecessary? See checks below
// otool -L Spike // to test which dynamic libraries are being used
// http://discussions.apple.com/thread.jspa?threadID=1741520
// http://www.cprogramming.com/tutorial/shared-libraries-linux-gcc.html

PARAMS * mp;
unsigned long int seed = 0;
gsl_rng * mSeed = NULL;
gsl_rng ** states = NULL;
int nThreads = 1;
char *pFile = NULL;

int main (int argc, const char * argv[])
{	
	//int set_error = 0;
	bool rerun = false;
	bool compress = true;
	bool seedFlag = false;
	int hours = 0;
	int mins = 0;
	
	char *imageArchive = NULL;
	//char *sfile = NULL;
	char *cpfile = "CLIparams.m";
	FILE * cli_FP = NULL;
	FILE * pipeFP = NULL;
	char syscmd[BUFSIZ]; // stdio.h : 1024
	char dlver[BUFSIZ];
	char *bptr = NULL;
	int th = 0;
	float proportion = 1.0;
#ifdef _OPENMP
	bool dynamic = false;
#endif
	bool genNewSeed = false;
	bool p_flag = false;	// Parameter (CLI) flag
	bool pf_flag = false;	// Parameter file flag
	bool ia_flag = false;	// Image archive flag
	int pINcount = 0;
	int dcount = 0; // Default parameter count
	int fcount = 0; // File parameter count
	int pcount = 0; // Parameter count
    int ocount = 0; // Parameters out count
	int err = 0;
	int slen = 0;
    char * error = NULL;
    void * dylib = NULL;
    
	bool recFlag = false;
	char * recList = NULL;
	
	bool skip_arg = false;
	//double wtime, begin;
	char timeStr[FNAMEBUFF];
    int c, result;
	struct tm *ts;
	char hostname[FNAMEBUFF];
	
	char schedule[BUFSIZ];
	strncpy(schedule, "<default>", BUFSIZ);
	
    printf("--------------------------------------------------------------------------------\n");
	
	time_t now = time(NULL);
	ts = localtime(&now);
	strftime(timeStr, FNAMEBUFF, "%a %d/%b/%Y %H:%M:%S", ts);
	err = gethostname(hostname, FNAMEBUFF);
	//char *user = getenv("USER");
	struct passwd *userinfo = getpwuid(geteuid());
	char *user = userinfo->pw_name;
	if (user && err != -1)
		printf("[%s] : Program started by %s@%s\n", timeStr, user, hostname);
	char cwd[BUFSIZ];
	if (getcwd(cwd, sizeof(cwd)) != NULL)
		fprintf(stdout, "DIR: %s\n", cwd); // Print directory
	else
		perror("getcwd() error");
    
// Move this section to a seperate header e.g. compiler.h 
#ifdef _OPENMP
#define OMP "T"
#else
#define OMP "F"    
#endif
    
#ifdef __GNUC__ // N.B. __GNUC__ is for any compiler implementing GNU compiler extensions, hence is defined for clang and llvm-gcc
#ifndef __has_feature 
#define __has_feature(x) 0 
#endif
    
#ifdef __llvm__ // Using LLVM backend
    // http://clang.llvm.org/docs/LanguageExtensions.html
    //printf("%d\n",__COUNTER__);
#ifdef __clang__ // Using Clang-LLVM
    // For a list of builtin defines type: clang -x c /dev/null -dM -E
    printf("Compiler: Clang-LLVM %s\n", __clang_version__);

#else // Using GCC-LLVM
    printf("Compiler: GCC-LLVM %s\n", __VERSION__);
#endif

    // Time of last modification of current source file...
    printf("Compiled on: %s | Optimization: %d | Debug: %d | OpenMP: %s\n", \
           __TIMESTAMP__, __OPTIMIZE__, DEBUG, OMP); 
    
#if __has_feature(c_static_assert) // Working? Relevent?
    printf("Includes support for compile-time assertions\n");
#else
    fprintf(stderr, "*** Warning: assert() disabled in parallel regions! ***\n");
#endif
    
#else // Using GCC
    printf("Compiler: %s | Optimization: %d | Debug: %d | OpenMP: %s\n", \
		   __VERSION__, __OPTIMIZE__, DEBUG, OMP);
    printf("Source modified on: %s\n",__TIMESTAMP__);
    printf("Compiled on: %s at %s\n", __DATE__, __TIME__);
#endif
#endif
    
#ifdef NDEBUG
	fprintf(stderr, "*** Warning: Executing without error checking! ***\n");
#endif
    
	if (strcmp(user, "nobody")==0)
		SIM.Xgrid = true;
	else
		if (getenv("OMP_SCHEDULE")) // not NULL string
			strncpy(schedule, getenv("OMP_SCHEDULE"), BUFSIZ);
	
	char * rsfile = RSFILE;
	
	
	printf("--------------------------------------------------------------------------------\n");
    printf("Checking for \"%s\" in current directory... \t\t\t    [%s]\n",DPFILE,\
           (file_exists(DPFILE))?"OK":"NO");
    printf("Checking for \"%s\" in current directory... \t\t    [%s]\n",rsfile,\
           (file_exists(rsfile))?"OK":"NO");
    // Check for GSL
    dylib = dlopen(GSLDIR"libgsl.dylib",RTLD_NOW);
    printf("Checking %s for GSL dyamic libraries... \t\t\t    [%s]\n",GSLDIR,(dylib)?"OK":"NO");
    if ((error = dlerror()) != NULL || !dylib)
        exit_error("main: libgsl.dylib check", error);
    else // dylib != NULL
        dlclose(dylib);
    // Check for System libraries
    dylib = dlopen("libSystem.dylib", RTLD_NOW);
    printf("Checking for System dynamic libraries... \t\t\t\t    [%s]\n",(dylib)?"OK":"NO");
    if ((error = dlerror()) != NULL || !dylib)
        exit_error("main: libSystem.dylib check", error);
    else // dylib != NULL
        dlclose(dylib);
    // Runtime OpenMP check using int omp_in_parallel(void);
#ifdef _OPENMP
#pragma omp parallel
	{
#pragma omp master//single
		{
            printf("Checking for OpenMP runtime parallelism... \t\t\t\t    [%s]\n",\
                   omp_in_parallel()?"OK":"NO");
        }
    }
#endif
    printf("--------------------------------------------------------------------------------\n");
    
	char exec[BUFSIZ];
	strncpy(exec, argv[0], sizeof(exec)-1);
	
	if (argc==1)
	{
		printf("%s usage:\n",argv[0]);
		printf("-c[lean]\t: Clean all dat and tbz files (including image archives!)\n");
		printf("-f <filename>\t: Pass parameter filename\n");
		printf("-r[erun]\t: Rerun simulation with the random seed in %s\n",rsfile);
		printf("-g[enerate]\t: Generate new random seed in %s and exit\n",rsfile);
		printf("-s <seed>\t: Explicitly pass random seed [0, (2^32)-1]\n");
		printf("-k <record list>: Pass list of neurons to be recorded\n");
		printf("-d[ynamic]\t: Set number of threads to be dynamic\n");
		printf("-m <proportion>\t: Set number of threads to a proportion of cores [0.0, 1.0]\n");
		printf("-t <threads>\t: Explicitly set the number of threads to use\n");
		printf("-p <parameter>\t: Pass a parameter string <name>=<value>\n");
		printf("--<parameter>\t: Pass a parameter string <name>=<value>\n");
		//printf("-i[mage] <directory>\t: Pass directory of filtered images [***incomplete***]\n");
		printf("-j <images>.tbz\t: Pass compressed image archive\n");
		printf("-u[ncompressed]\t: Prevent data compression\n");
		printf("-x[grid]\t: Set as an Xgrid simulation i.e. print progress information\n");
		printf("================================================================================\n");
		return 0;
	}
	
	while (--argc > 0 && (*++argv)[0] == '-')
	{
		//skip_arg = 0;
		while (!skip_arg && (c = *++argv[0]))
		{
			switch (c)
			{
				case 'c':	// Clean directory of .dat and .tbz files
					system("rm *.dat *.tbz");
					break;
					
				case 'f':	// Parameter file name
					pf_flag = true;
					pFile = myalloc(strlen(*++argv)+1); //sizeof(char)==1 guaranteed
					strcpy(pFile, *argv);
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
					fprintf(stderr, "*** -k: Specifying neurons for recording is not yet implemented! ***\n");
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
					fprintf(stderr, "*** -d: OpenMP disabled! ***\n");
#endif
					break;
					
				case 'm':	// Set the proportion of threads from the CLI [0.0, 1.0]
#ifdef _OPENMP
					proportion = atof(*++argv);
//#ifndef __llvm__
					assert(0.0 < proportion && proportion <= 1.0);
//#endif
					omp_set_num_threads(round(omp_get_num_procs()*proportion));
#else
					fprintf(stderr, "*** -m %f: OpenMP disabled! ***\n",proportion);
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
					fprintf(stderr, "-t %d: OpenMP disabled\n",nThreads);
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
					
				case '-':	// Equivalent to '-p ' but combines the arguments
					if (!p_flag)
					{
						cli_FP = myfopen(cpfile, "w");
						p_flag = true;
					}
					fprintf(cli_FP, "%s;\n", ++argv[0]); // Advance to next char address
					skip_arg = true;
					//argc--;
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
					//printf("Warning: tbz archives should be removed to prevent analysis of them!\n");
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
			printf("OMP: (%d/%d)\t{OMP_DYNAMIC=%s, OMP_NESTED=%s, OMP_SCHEDULE=%s}\n", \
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
	
//	printf("This program is compiled with GSL version %s.\n", GSL_VERSION);
//#if DEBUG > 1	// Print GSL verison and location
//	system("gsl-config --prefix --version");
//#endif

    //if (!SIM.Xgrid)
    //{
	// Print GSL verison and location
	if ((pipeFP = popen("/opt/local/bin/gsl-config --prefix --version", "r")))
	{
		fgets(syscmd, sizeof(syscmd)-1, pipeFP);
		fgets(dlver, sizeof(dlver)-1, pipeFP);
	}
	pclose(pipeFP);
	if ((bptr = strpbrk(syscmd, "\r\n"))) //strstr(syscmd, '\n')
		*bptr = '\0';
	printf("GSL: Compiled with v%s, found dynamic libraries v%4.2f at: %s\n", \
           GSL_VERSION,atof(dlver),syscmd);
	// if atof(dlver) < GSL_MIN
    //}
    
    // Initialise random seed
	const gsl_rng_type * T = gsl_rng_default; // Set RNG type
	gsl_rng_env_setup(); // http://www.gnu.org/software/gsl/manual/html_node/Random-number-environment-variables.html
	mSeed = gsl_rng_alloc(T); // Used for serial sections with randomness
    
	if (genNewSeed) // Generate a new random seed file and exit
	{
		seed = (unsigned long) time((time_t *) NULL);
		seed %= seedMod;
		randSeedFP = myfopen(rsfile, "w");
		fprintf(randSeedFP, "mSeed: \t%ld\n", seed);
		fclose(randSeedFP);
        printf("New seed generated in %s: %ld (%s) <GSL v%s>\n", rsfile, seed, gsl_rng_name(mSeed), GSL_VERSION); // Generator type not strictly necessary here
		return 0;
	}
	
    if (!seedFlag)
	{
		if (rerun) // Also rerun with parameters.m?
		{
			randSeedFP = myfopen(rsfile, "r");
			if ((sString = fgets(buffer, sizeof(buffer), randSeedFP)) != NULL) //while
				seed = atol(strrchr(sString,':')+1); //ans[count++]
			fclose(randSeedFP); 
			printf("Rerunning simulation with %s: %ld (%s) <GSL v%s>\n", rsfile, seed, gsl_rng_name(mSeed), GSL_VERSION);
		}
		else
		{
			seed = (unsigned long) time((time_t *) NULL);
			seed %= seedMod;
			fprintf(stderr, "*** Warning: Creating new seed in %s: %ld (%s) <GSL v%s> ***\n", rsfile, seed, gsl_rng_name(mSeed), GSL_VERSION);
			randSeedFP = myfopen(rsfile, "w");
			fprintf(randSeedFP, "mSeed: \t%ld\n", seed);
			fclose(randSeedFP);
		}
	}
	gsl_rng_set(mSeed, seed); //gsl_rng_set(mSeed, -idum);
    
	// Allocate and initialise model parameters structure
	mp = myalloc(sizeof(*mp)); // Place in getParameters with default init?
    mp->initialised = false;
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
				fprintf(stderr, "*** Warning! Undersized buffer: %s ***", syscmd);
			syserr = system(syscmd);
			if (syserr)
				EE("Error extracting image archive");
			else
				printf("Images successfully extracted to %s\n", mp->imgDir);
		}
	}
    
	SIM.minTau = BIG;
	
	// Read in parameters from .m file
	printf("I/O: Processing parameters: \"%s\"", !pFile ? DPFILE : pFile);
	if (p_flag)
		fclose(cli_FP);
	dcount = read_parameters(mp, DPFILE);
	fcount = (pFile != NULL) ? read_parameters(mp, pFile) : 0;
	pcount = (p_flag) ? read_parameters(mp, cpfile) : 0;
	
	if (!mp->useFilteredImages)
		assert(pcount == pINcount);
	//printf(" {%d,%d,%d}\tParsing complete!\n", dcount, fcount, pcount);
		
	// Print parameters to MPFILE (parameters.m)
	ocount = printParameters(mp, MPFILE); // Variables to read into Matlab
	//printf("%d parameters written to %s\n", pcount, MPFILE);
	
    printf(" {%d,%d,%d} --> \"%s\" {%d} Done!\n", dcount, fcount, pcount, MPFILE, ocount);
	
	// Create a random seed for each thread to ensure thread safety
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
	
	if (recFlag) // Make recordSet global, prevent random choice in init_network and free at the end
	{
		/*int ** recordSet = myalloc(mp->nLayers * sizeof(*recordSet));
		char * tstr = NULL;
		int count = 0;
		strtok(list, ";");
		while (tstr && (*trim(tstr) != ']') && (tstr != '\0'))
		for (l=0; l<mp->nLayers; l++)
		{
			mp->vRecords[l] = parseIntVector(list, &recordSet[l]);
			tstr = strtok(NULL, ";");
		}*/
		
		myfree(recList);
	}
	
	// Print minimum tau and DT to nearest microsecond
	printf("TAU: Smallest time constant = %.3f ms | DT = %.3f ms\n", SIM.minTau*1000, mp->DT*1000); 
	if (mp->DT >= 2*SIM.minTau) // CHECK THIS
		fprintf(stderr, "*** Warning: Forward Euler stability condition violated! ***\n");
	// Display dynamic libraries: otool -L ~/bin/SpikeNet/Debug/Spike
	
#ifdef _OPENMP // Use omp function omp_get_wtime
	//double begin = omp_get_wtime();
	SIM.start = omp_get_wtime();
	SIM.elapsed = 0.0;
#else
	time_t start = time(NULL);
#endif
	
	if (!SIM.Xgrid && !mp->loadWeights) // Remove *.dat and *.tbz
	{
		system("rm *.dat"); //system("rm *.dat *.tbz");
		if (ia_flag)
		{
			if(snprintf(syscmd, BUFSIZ, "find *.tbz ! -name %s -delete",imageArchive) >= BUFSIZ)
				fprintf(stderr, "*** Warning! Undersized buffer: %s ***", syscmd);
			if (system(syscmd)) // Delete *.tbz except image archive
				printf("Archive files successfully cleaned!\n");
			else
				EE("Error cleaning archive files!"); //exit_error("main.c", "Error cleaning archive files!\n");		
		}
		else
			system("rm *.tbz"); // Delete *.tbz
	}

	if (mp->loadWeights)
	{
		// Pass an archive with all relevant dat files with CLI flag e.g. network.tbz
		const char * suffix = "";
		char fname[FNAMEBUFF];
		slen = snprintf(fname, FNAMEBUFF, "L0affNeuronsElE%s.dat", suffix);
		assert(slen < FNAMEBUFF);
		if (!file_exists(fname))
		{
			if (file_exists("connectivity.tbz"))
				system("tar -xvf connectivity.tbz");
			else
				EE("No connectivity files to load"); //exit_error("main", "No connectivity files to load");
		}
		if (mp->nLayers > 1)
		{
			slen = snprintf(fname, FNAMEBUFF, "L1affNeuronsEfE%s.dat", suffix);
			assert(slen < FNAMEBUFF);
			if (!file_exists(fname))
			{
				if (file_exists("postTraining.tbz"))
					system("tar -xvf postTraining.tbz");
				else
					EE("No weights files to load"); //exit_error("main", "No weights files to load");
			}
		}
	}
	
	
	/***** RUN SIMULATION *****/
	result = spike(mp);
	/**************************/
	
	
	// Compress data files for crash-free xgrid! '-j' Uses bzip (*.tbz equivalent to *.tar.bz2)
	// Append files to fileList and call system(syscmd); once and keep fileList
	//snprintf(syscmd, BUFSIZ, "tar -cjvf %s.tbz %s > fileList","connectivity","*affNeurons.dat");
	
	/*if (!SIM.Xgrid) // /sbin/md5
		system("md5 *.dat > datHashs.txt");*/
	system("shasum *.dat > datHashs.txt"); // /usr/bin/shasum
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
		
		if (!(mp->useFilteredImages || mp->stimGroups))
			if ((syserr = system("tar -cjf stimuli.tbz *stimuli.dat stimuli.m")) == 0)
				system("tar -tf stimuli.tbz | xargs rm");
		
		if(mp->nRecordsPL)
		{
			if (mp->priorPhases)
				if ((syserr = system("tar -cjf PPrecords.tbz R*_PP_*.dat")) == 0)
					system("tar -tf PPrecords.tbz | xargs rm");

			if ((syserr = system("tar -cjf records.tbz R*.dat")) == 0)
				system("tar -tf records.tbz | xargs rm");
		}
			
		if (mp->printConnections)
		{
			if (mp->SOM)
				syserr = system("tar -cjf connectivity.tbz *affNeurons*.dat *affDelays*.dat *dist*.dat");
			else
            {
                if (mp->axonDelay)
                    syserr = system("tar -cjf connectivity.tbz *affNeurons*.dat *affDelays*.dat"); //system("tar --remove-files -cjvf connectivity.tbz *affNeurons*.dat > fileList");
                else
                    syserr = system("tar -cjf connectivity.tbz *affNeurons*.dat");
            }
			if (!syserr)
				system("tar -tf connectivity.tbz | xargs rm");
		}
		
		if (mp->pretrain)
		{
			if (mp->priorPhases)
				if ((syserr = system("tar -cjf PPpreTraining.tbz pt_PP_*.dat")) == 0)
					system("tar -tf PPpreTraining.tbz | xargs rm");
		
			if ((syserr = system("tar -cjf preTraining.tbz pt*.dat")) == 0)
				system("tar -tf preTraining.tbz | xargs rm");
		}
		
		if (mp->train)
		{
			if (mp->priorPhases)
				if ((syserr = system("tar -cjf PPtraining.tbz _PP_E*.dat")) == 0) // 2> tar_err
					system("tar -tf PPtraining.tbz | xargs rm");
			
			if ((syserr = system("tar -cjf training.tbz E*.dat")) == 0) // 2> tar_err
				system("tar -tf training.tbz | xargs rm");
		}
		
		if (mp->priorPhases)
			if ((syserr = system("tar -cjf PPpostTraining.tbz _PP_*.dat")) == 0) // 2> tar_err
				system("tar -tf PPpostTraining.tbz | xargs rm");
				
		if ((syserr = system("tar -cjf postTraining.tbz L*Spikes.dat L*weights*.dat")) == 0)
			system("tar -tf postTraining.tbz | xargs rm");
		//system(syscmd);
		//system("rm fileList");
		
		printf("Data Compressed!\n");
		fflush(stdout);
	}
	//#pragma omp section
	/*if (!SIM.Xgrid) // Print md5 #'s // /sbin/md5
	{
		//system("md5 Spike");
		system("md5 parameters.m"); // shasum
		system("md5 datHashs.txt");
		//system("md5 *.tbz"); // Contains metadata (e.g. timestamps) which will give different #s
	 }*/
		
	slen = snprintf(syscmd, sizeof(syscmd)-1, "shasum %s", exec);
#ifndef __llvm__
	assert(slen < (signed) sizeof(syscmd));
#endif
	system(syscmd);
	system("shasum parameters.m");
	system("shasum datHashs.txt");
	} // End of section
		
	// Clean up
#pragma omp section
	if (pf_flag)
		myfree(pFile);
	
#pragma omp section
	if (ia_flag)
	{
		myfree(imageArchive);
		if(snprintf(syscmd, BUFSIZ, "rm -R %s/", mp->imgDir) >= BUFSIZ)
			fprintf(stderr, "*** Warning! Undersized buffer: %s ***", syscmd);
		system(syscmd);	// Delete expand image files
	}
        // Print out input/output file list? array of structs with a bool and filename string...
	
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
	double wtime = omp_get_wtime() - SIM.start; //begin;
	//double integral;
	//double fraction = modf(wtime, &integral);
	//duration = (time_t) round(integral);
	hours = floor(wtime/3600);
	wtime -= hours*3600;
	mins = floor(wtime/60);
	wtime -= mins*60; //secs = wtime - (mins*60) - (hours*3600);
	snprintf(timeStr, FNAMEBUFF, "%d:%02d:%06.3lf (%d Threads)",\
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
		fprintf(stderr, "*** Simulation aborted after %s! ***\n",timeStr);
		return 1;
	}
	//printf("--------------------------------------------------------------------------------\n");
    printf("================================================================================\n");
	
    return 0;
}
