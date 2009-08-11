/*
 *  array_utils.c
 *  Spike
 *
 *  Created by Ben Evans on 04/02/2009.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "array_utils.h"


void normalize(double *vec, int size)
{
	int i;
	double sf = 0.0; // Scaling factor
	
	for (i=0; i<size; i++)
		sf += vec[i] * vec[i];
	
	sf = sqrt(sf);
	
	if (sf>=EPS)
		for (i=0; i<size; i++)
			vec[i] /= sf;
	
	return;
}

double ndp(double *vec_x, double *vec_y, int size)
{
	int i;
	double sxsq = 0.0;
	double sysq = 0.0;
	double sxy = 0.0;
	double root = 0.0;
	// double eps = ...	
	
	if (size <= 0)
	{
		return 9.999999;
		printf("NDP error: Vector lengths must be >= 1");
	}
	
	for (i=0; i<size; i++)
	{
		sxsq += vec_x[i] * vec_x[i];
		sysq += vec_y[i] * vec_y[i];
		sxy += vec_x[i] * vec_y[i];
	}
	
	root = sqrt(sxsq * sysq);
	
	// double cos_theta = ((fabs(root)<=eps) ? 0.0 : (sxy/root));
	// return cos_theta;
	
	
	return ((fabs(root)<=EPS) ? 0.0 : (sxy/root));
	/*
	 if (fabs(root) <= EPS)
	 {
	 return 0.0;
	 }
	 else
	 return (sxy / root); // returns cos @
	 }
	 */
}

void print_irow(FILE * file_ptr, int * data, int ncols)
{
	int j;
	for (j=0; j<ncols; j++)
		fprintf(file_ptr, "%d\t", data[j]);
	fprintf(file_ptr, "\n");
	fflush(file_ptr);
	return;
}

void print_iarray(FILE *file_ptr, int **data, int nrows, int ncols)
{ // Change to void **data in order to take ints and doubles
	int i = 0;
	int j = 0;
	//FILE * file_ptr;
	//file_ptr = fopen(file_name, "w"); // See nifty trick #2
	assert(file_ptr != NULL);
	for (i=0; i<nrows; i++)
	{
		for (j=0; j<ncols; j++)
		{
			fprintf(file_ptr, "%d\t", data[i][j]);
		}
		fprintf(file_ptr, "\n");
		fflush(file_ptr); // Is this necessary?
	}
	//fclose(file_ptr);
	return;
}

void print_farray(FILE *file_ptr, float **data, int nrows, int ncols)
{ // Change to void **data in order to take ints and doubles
	int i = 0;
	int j = 0;
	//FILE * file_ptr;
	//file_ptr = fopen(file_name, "w"); // See nifty trick #2
	assert(file_ptr != NULL);
	for (i=0; i<nrows; i++)
	{
		for (j=0; j<ncols; j++)
		{
			fprintf(file_ptr, "%f\t", data[i][j]);
		}
		fprintf(file_ptr, "\n");
		fflush(file_ptr); // Is this necessary?
	}
	//fclose(file_ptr);
	return;
}


int ** get_2D_iarray(int nrows, int ncols, int init)
{
	int i, j;
	// sizeof(*array) = sizeof(int *)
	// sizeof(**array) = sizeof(int)
	
	int **array2D = myalloc(nrows * sizeof(int *));
	array2D[0] = myalloc(nrows * ncols * sizeof(int));
	for (i=1; i<nrows; i++) 
		array2D[i] = array2D[0] + (i * ncols);
	
	if (init == 0)
		memset(array2D[0], 0, nrows*ncols*sizeof(int)); //sizeof(**array2D)); 
	else
	{
		for (i=0; i<nrows; i++)
			for (j=0; j<ncols; j++)
				array2D[i][j] = init; //(i * ncols) + j;
	}
	
	return array2D;	
}

void free_2D_iarray(int ** array2D, int nrows)
{
	int i;
	
/*	for (i=0; i<nrows; i++)
		free(array2D[i]);*/
	free(*array2D); // Since the memory is contiguous
	
	free(array2D);
	return;
}

int *** get_3D_iarray(int nlays, int nrows, int ncols, int init)
{
	int x, y, z;

	// sizeof(*array) = sizeof(int **)
	// sizeof(**array) = sizeof(int *)
	// sizeof(***array) = sizeof(int)
	
	int *space = myalloc(nlays * nrows * ncols * sizeof(*space)); // sizeof(int)
	int ***array3D = myalloc(nlays * sizeof(*array3D)); // sizeof(int **)
	
	for (z=0; z<nlays; z++)
	{
		array3D[z] = myalloc(nrows * sizeof(**array3D)); // sizeof(int *)
		for (y=0; y<nrows; y++)
			array3D[z][y] = space + (z*ncols*nrows) + (y*ncols);
	}
	
	/* Assign values to 3D array */
	if (init == 0)
		memset(array3D[0][0], 0, nlays*nrows*ncols*sizeof(*space));
	else
	{
		for (z=0; z<nlays; z++)
			for (y=0; y<nrows; y++)
				for (x=0; x<ncols; x++)
					array3D[z][y][x] = init; //(z * nrows * ncols) + (y * ncols) + x;
	}
	
	return array3D;
}


void free_3D_iarray(int *** array3D, int nlays, int nrows)
{
	int y,z;
	free(**array3D);
	for (z=0; z<nlays; z++)
	{
		/*for (y=0; y<nrows; y++)
			free(array3D[z][y]);*/
		free(array3D[z]);
	}
	free(array3D);
	return;
}


float ** get_2D_farray(int nrows, int ncols, float init)
{
	int i, j;
	// sizeof(*array) = sizeof(int *)
	// sizeof(**array) = sizeof(int)
	
	float **array2D = myalloc(nrows * sizeof(float *));
	array2D[0] = myalloc(nrows * ncols * sizeof(float));
	for (i=1; i<nrows; i++) 
		array2D[i] = array2D[0] + (i * ncols);
	
	if (fabs(init) <= EPS) // (int) init == 0
		memset(array2D, 0, nrows*ncols*sizeof(**array2D));
	else
	{
		for (i=0; i<nrows; i++)
			for (j=0; j<ncols; j++)
				array2D[i][j] = init; //(i * ncols) + j;
	}
	
	return array2D;	
}

void free_2D_farray(float ** array2D, int nrows)
{
	int i;
	
	/*for (i=0; i<nrows; i++)
		free(array2D[i]);*/
	free(array2D[0]); // Since memory is contiguous
	free(array2D);
	return;
}

float *** get_3D_farray(int nlays, int nrows, int ncols, float init)
{
	int x, y, z;
	
	// sizeof(*array) = sizeof(int **)
	// sizeof(**array) = sizeof(int *)
	// sizeof(***array) = sizeof(int)
	
	float *space = myalloc(nlays * nrows * ncols * sizeof(*space)); // sizeof(int)
	float ***array3D = myalloc(nlays * sizeof(*array3D)); // sizeof(int **)
	
	for (z=0; z<nlays; z++)
	{
		array3D[z] = myalloc(nrows * sizeof(**array3D)); // sizeof(int *)
		for (y=0; y<nrows; y++)
			array3D[z][y] = space + (z*ncols*nrows) + (y*ncols);
	}
	
	/* Assign values to 3D array */
	if (fabs(init) <= EPS)
		memset(array3D, 0, nlays*nrows*ncols*sizeof(***array3D));
	else
	{
		for (z=0; z<nlays; z++)
			for (y=0; y<nrows; y++)
				for (x=0; x<ncols; x++)
					array3D[z][y][x] = init; //(z * nrows * ncols) + (y * ncols) + x;
	}
	
	return array3D;
}


void free_3D_farray(float *** array3D, int nlays, int nrows)
{
	int y,z;
	free(**array3D);
	for (z=0; z<nlays; z++)
	{
		/*for (y=0; y<nrows; y++)
			free(array3D[z][y]);*/
		free(array3D[z]);
	}
	free(array3D);
	return;
}


