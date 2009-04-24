/*
 *  rng.c
 *  Spike
 *
 *  Created by Ben Evans on 18/08/2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "rng.h"
//#include "parameters.h"

/* Code for Random Number Generators */
int set_random_seeds(int RERUN)
{
	/* This routine sets the random seeds */
	
	char * dummy;
	char * string, buff[BUFFER];
	int count = 0;
	char * rsfile = RSFILE;
	// FILE * random_seeds_ptr;
	//long file_location; // Variable to hold the current location in the file
	
	if (RERUN == 1)
	{	// Generalise the code to optionally accept a filename e.g. random_seeds5b.dat
		/* File Format:
		 rands1_new
		 ran3
		 gasdev */
		printf("Rerunning simulation with %s\n",rsfile);
		if ((random_seeds_ptr = fopen(rsfile, "r")) == NULL)
		{	
			fprintf(stderr, "Error: %s could not be opened!\n", rsfile);
			fclose(random_seeds_ptr);
			return 1;
		}
		while ((string = fgets(buff, sizeof(buff), random_seeds_ptr)) != NULL) 	/* Read next line */
		{
			//sscanf(s, "%s %ld", dummy, &ans[count]);
			dummy = (char *) strtok(string, ":");
			ans[count] = atoi(strtok(NULL, ":"));
			count++;
		}
		/*while(!feof(random_seeds_ptr))
		{
			fscanf(random_seeds_ptr,"%s %ld", dummy, &ans[count]);
			count++;	
		}*/
		fclose(random_seeds_ptr); 
	}
	else
	{
		random_seeds_ptr = fopen(rsfile, "w");
	}	
	
	(void) rnd_new(0, 0);	/* Reseed random number generators */
		
	if (SEED_RAN3 == 0)
	{
		if (RERUN == 0)
		{
			idum = 0;
			while (idum == 0)
			{
				idum = (long) -time((time_t *) NULL);
				idum %= MSEED;
			}
		}
		else if (RERUN == 1)
		{
			idum = ans[1]; // Second line
		}
	}
	else
	{
		idum = SEED_RAN3;
	}
	printf("Random number seed for ran3 is %ld\n", idum);
	if (RERUN != 1)
	{
		fprintf(random_seeds_ptr, "ran3: \t\t%ld\n", idum);
	}
	(void) ran3(&idum);
	
	if (SEED_GASDEV == 0)
	{
		if (RERUN == 0)
		{
			idum_gasdev = 0;
			while (idum_gasdev == 0)
			{
				idum_gasdev = (long) -time((time_t *) NULL);
				idum_gasdev %= MSEED;
			}
		}
		else if (RERUN == 1)
		{
			idum_gasdev = ans[2]; // Third line
		}
	}
	else
	{
		idum_gasdev = SEED_GASDEV;
	}
	printf("Random number seed for gasdev is %ld\n", idum_gasdev);
	if (RERUN != 1)
	{
		fprintf(random_seeds_ptr, "gasdev: \t%ld\n", idum_gasdev);
	}
	printf("\n");
	(void) gasdev(&idum_gasdev);
	
	fclose(random_seeds_ptr);
	
	return 0;
}


int rands1_new(int mn, int mx, int *iv, int mode)
{
	static int ar[CLASS_SIZE];	/* Record of previous stimuli */
	static int nnov;		/* Number of novel stimuli */
	static int ntr;
	int     ret;
	const int nmax = CLASS_SIZE;
	int     i, n;
	int     ok_flag = 0;
	int     end_flag = 0;
	
	n = mx - mn + 1;
	if (n >= nmax) // if (n<1) //corrected BDE 30/7/08
		return mn;
	
	/* Prepare for new random sequence by clearing array if mx <= mn */
	if (n <= 1) 
	{				/* Clear array */
		for (i = 0; i < nmax; i++)
			ar[i] = 0;
		nnov = 1;
		ntr = 0;
		return 0;
	}
	
	ntr++;
	do
	{				/* Get a random number */
		ret = rnd_new(mn, mx);
		i = ret - mn;
		if (ar[i] != 0)
		{				/* This number has been used before */
			if (mode == 1 || ar[i] == -1)
			{			/* Check for all numbers used up */
				i = 0;
				while (i < n)
				{			/* We've found an unused one */
					if (ar[i] != -1)
						break;	/* must NOT use the ++ here !!! */
					i++;
				}
				if (i == n)
				{			/* They've all been used, restart */
					for (i = 0; i < nmax; i++)
					{
						ar[i] = 0;
					}
					end_flag = 1;
				}
			}
			else
			{			/* Calculate number of intervening stimuli */
				*iv = ntr - ar[i] - 1;
				ar[i] = -1;
				nnov = 0;		/* Reset novel count */
				ok_flag = 1;
			}
		}
		else
		{
			if (mode == 1)
			{
				ar[i] = -1;
			}
			else
			{
				nnov++;
				if (nnov > 3)
					continue;	/* Too many novels */
				ar[i] = ntr;
			}
			if (end_flag)
				ret = (-ret);	/* Show end of sequence */
			ok_flag = 1;
			*iv = -1;
		}
	}
	while (!ok_flag);
	
	/* printf(" rands1_new %4d returned \n", ret);      */
	return ret;
}


/* Return a random number in the range [mn, mx]
 Reseed generator if mx <= mn */
int rnd_new(int mn, int mx)
{
	int d, seed;
	
	d = mx - mn + 1;
	
	/* seed random number generator if mx <= mn */
	if (mx <= mn)
	{
		if (SEED_RANDS1_NEW == 0)
		{
			if (RERUN == 0)
			{
				srand(seed = (int) time((time_t *) NULL));
			}
			else if (RERUN == 1)
			{
				srand(seed = ans[0]);
			}
		}
		else
		{
			seed = SEED_RANDS1_NEW;
			srand(SEED_RANDS1_NEW);
		}
		printf("Random number seed for rands1_new is %ld\n", seed);
		if (RERUN != 1)
		{
			fprintf(random_seeds_ptr, "rands1_new: \t%ld\n", seed);
		}
		return 0;
	}
	return (rand() % d) + mn;
}


float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;
	
	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;
	
	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


float ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;
	
	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}