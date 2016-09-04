/*
 *  array_utils.c
 *  Spike
 *
 *  Created by Ben Evans on 04/02/2009.
 *  Copyright 2009 University of Oxford. All rights reserved.
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
    if (sf>EPS)
    {
        sf = 1/sf;
        for (i=0; i<size; i++)
			vec[i] *= sf;
    }
	
	return;
}

double ndp(double *vec_x, double *vec_y, int size)
{
    // Return the normalised dot product (cos \theta) between two vectors
	int i;
	double sxsq = 0.0;
	double sysq = 0.0;
	double sxy = 0.0;
	double root = 0.0;

    assert(size >= 1);
	
	for (i=0; i<size; i++)
	{
		sxsq += vec_x[i] * vec_x[i];
		sysq += vec_y[i] * vec_y[i];
		sxy += vec_x[i] * vec_y[i];
	}
	
	root = sqrt(sxsq * sysq);
	
	return (fabs(root)>EPS) ? (sxy/root) : 0.0;
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
	assert(file_ptr != NULL);
	for (i=0; i<nrows; i++)
	{
		for (j=0; j<ncols; j++)
			fprintf(file_ptr, "%d\t", data[i][j]);
		fprintf(file_ptr, "\n");
		fflush(file_ptr);
	}
	return;
}

void print_frow(FILE * file_ptr, float * data, int ncols)
{
	int j;
	for (j=0; j<ncols; j++)
		fprintf(file_ptr, "%G\t", data[j]);
	fprintf(file_ptr, "\n");
	fflush(file_ptr);
	return;
}

void print_farray(FILE *file_ptr, float **data, int nrows, int ncols)
{ // Change to void **data in order to take ints and doubles
	int i = 0;
	int j = 0;
	assert(file_ptr != NULL);
	for (i=0; i<nrows; i++)
	{
		for (j=0; j<ncols; j++)
			fprintf(file_ptr, "%G\t", data[i][j]);
		fprintf(file_ptr, "\n");
		fflush(file_ptr);
	}
	return;
}


void print_darray(FILE *file_ptr, double **data, int nrows, int ncols)
{ // Change to void **data in order to take ints and doubles
	int i = 0;
	int j = 0;
	assert(file_ptr != NULL);
	for (i=0; i<nrows; i++)
	{
		for (j=0; j<ncols; j++)
			fprintf(file_ptr, "%lf\t", data[i][j]);
		fprintf(file_ptr, "\n");
		fflush(file_ptr);
	}
	return;
}


void ** array2d(size_t rows, size_t cols, size_t value_size)
{
	if (rows < 1 || cols < 1) // Catch errors where one or more dimensions is 0 or -ve
		return NULL;
	
    size_t index_size = sizeof(void *) * rows;
    size_t store_size = value_size * rows * cols;
	
    char * a = myalloc(index_size + store_size);
	
    memset(a + index_size, 0, store_size); // Be careful with memsets rezeroing the array
	size_t i=0;
    for(i = 0; i < rows; ++i)
        ((void **)a)[i] = a + index_size + i * cols * value_size;
	
    return (void **)a;
}


int ** get_2D_iarray(int nrows, int ncols, int init)
{
	int i, j;
	// sizeof(*array) = sizeof(int *)
	// sizeof(**array) = sizeof(int)
	
	int **array2D = myalloc(nrows * sizeof(int *));
	array2D[0] = myalloc(nrows * ncols * sizeof(int));
	if (!array2D[0]) // if (array2D[0] == NULL)
		return NULL;
	for (i=1; i<nrows; i++) 
		array2D[i] = array2D[0] + (i * ncols);
	
	if (init == 0)
		memset(array2D[0], 0, nrows*ncols*sizeof(int)); //sizeof(**array2D)); 
	else
		for (i=0; i<nrows; i++)
			for (j=0; j<ncols; j++)
				array2D[i][j] = init; //(i * ncols) + j;
	
	return array2D;	
}

void free_2D_iarray(int ** array2D)//, int nrows)
{
	if (!array2D)
		return;
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
	if (!space)
		return NULL;
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
		for (z=0; z<nlays; z++)
			for (y=0; y<nrows; y++)
				for (x=0; x<ncols; x++)
					array3D[z][y][x] = init; //(z * nrows * ncols) + (y * ncols) + x;
	
	return array3D;
}


void free_3D_iarray(int *** array3D, int nlays)//, int nrows)
{
	int z;
	if (!array3D)
		return;
	free(**array3D);
	for (z=0; z<nlays; z++)
		free(array3D[z]);
	free(array3D);
	return;
}


float ** get_2D_farray(int nrows, int ncols, float init)
{
	int i, j;
	// sizeof(*array) = sizeof(int *)
	// sizeof(**array) = sizeof(int)
	
	if (nrows == 0 || ncols == 0)
		return NULL;
	
	float **array2D = myalloc(nrows * sizeof(float *));
	array2D[0] = myalloc(nrows * ncols * sizeof(float));
	if (!array2D[0])
		return NULL;
	for (i=1; i<nrows; i++) 
		array2D[i] = array2D[0] + (i * ncols);
	
	if (fabs(init) <= EPS) // (int) init == 0
		memset(array2D[0], 0, nrows*ncols*sizeof(**array2D));
	else
		for (i=0; i<nrows; i++)
			for (j=0; j<ncols; j++)
				array2D[i][j] = init; //(i * ncols) + j;
	
	return array2D;	
}

void free_2D_farray(float ** array2D)//, int nrows)
{
	if (!array2D)
		return;
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
	
	float *space = myalloc(nlays * nrows * ncols * sizeof(*space)); // sizeof(float)
	if (!space)
		return NULL;
	float ***array3D = myalloc(nlays * sizeof(*array3D)); // sizeof(float **)
	
	for (z=0; z<nlays; z++)
	{
		array3D[z] = myalloc(nrows * sizeof(**array3D)); // sizeof(float *)
		for (y=0; y<nrows; y++)
			array3D[z][y] = space + (z*ncols*nrows) + (y*ncols);
	}
	
	/* Assign values to 3D array */
	if (fabs(init) <= EPS)
		memset(array3D[0][0], 0, nlays*nrows*ncols*sizeof(***array3D));
	else
		for (z=0; z<nlays; z++)
			for (y=0; y<nrows; y++)
				for (x=0; x<ncols; x++)
					array3D[z][y][x] = init; //(z * nrows * ncols) + (y * ncols) + x;
	
	return array3D;
}


void free_3D_farray(float *** array3D, int nlays)//, int nrows)
{
	int z;
	if (!array3D)
		return;
	free(**array3D);
	for (z=0; z<nlays; z++)
		free(array3D[z]);
	free(array3D);
	return;
}


float **** get4Dfarray(int D1, int D2, int D3, int D4, float init)
{
	int w, x, y, z;
	float *space = myalloc(D1 * D2 * D3 * D4 * sizeof(*space));
	float ****array4D = myalloc(D1 * sizeof(*array4D));
	for (w=0; w<D1; w++)
	{
		array4D[w] = myalloc(D2 * sizeof(**array4D));
		for (x=0; x<D2; x++)
		{
			array4D[w][x] = myalloc(D3 * sizeof(***array4D));
			for (y=0; y<D3; y++)
				array4D[w][x][y] = space + (w*D2*D3*D4) + (x*D3*D4) + (y*D4);
		}
	}
	
	//Initialise
	if (fabs(init) < EPS)
		memset(array4D[0][0][0], 0, D1*D2*D3*sizeof(****array4D));
	else
		for (w=0; w<D1; w++)
			for (x=0; x<D2; x++)
				for (y=0; y<D3; y++)
					for (z=0; z<D4; z++)
						array4D[w][x][y][z] = init;
	return array4D;
}


void free4Dfarray(float **** array4D, int D1, int D2) // Can the free functions be type independent? i.e. (void **** 4Darray, int D1, int D2)
{
	int w, x;
	if (!array4D)
		return;
	myfree(***array4D); // Free storage
	for (w=0; w<D1; w++)
	{
		for (x=0; x<D2; x++)
			myfree(array4D[w][x]);
		myfree(array4D[w]);
	}
	myfree(array4D);
	return;
}


float ******* get_7D_farray(int D1, int D2, int D3, int D4, int D5, int D6, int D7, float init)
{
	int t, u, v, w, x, y, z;
	
	// sizeof(*array) = sizeof(int **)
	// sizeof(**array) = sizeof(int *)
	// sizeof(***array) = sizeof(int)
	
	float *space = myalloc(D1 * D2 * D3 * D4 * D5 * D6 * D7 * sizeof(*space)); // sizeof(int)
	if (!space)
		return NULL;
	float *******array7D = myalloc(D1 * sizeof(*array7D)); // sizeof(int ******)
	
	for (t=0; t<D1; t++)
	{
		array7D[t] = myalloc(D2 * sizeof(**array7D)); // sizeof(int *****)
		for (u=0; u<D2; u++)
		{
			array7D[t][u] = myalloc(D3 * sizeof(***array7D));
			for (v=0; v<D3; v++)
			{
				array7D[t][u][v] = myalloc(D4 * sizeof(****array7D));
				for (w=0; w<D4; w++)
				{
					array7D[t][u][v][w] = myalloc(D5 * sizeof(*****array7D));
					for (x=0; x<D5; x++)
					{
						array7D[t][u][v][w][x] = myalloc(D6 * sizeof(******array7D));
						for (y=0; y<D6; y++)
						{
							array7D[t][u][v][w][x][y] = space + (t*D2*D3*D4*D5*D6*D7) \
							+ (u*D3*D4*D5*D6*D7) + (v*D4*D5*D6*D7) + (w*D5*D6*D7) + (x*D6*D7) + (y*D7);
						}
					}
				}
			}
		}
	}
	
	/* Assign values to 3D array */
	if (fabs(init) <= EPS)
		memset(array7D[0][0][0][0][0][0], 0, D1*D2*D3*D4*D5*D6*D7*sizeof(*******array7D));
	else
		for (t=0; t<D1; t++)
			for (u=0; u<D2; u++)
				for (v=0; v<D3; v++)
					for (w=0; w<D4; w++)
						for (x=0; x<D5; x++)
							for (y=0; y<D6; y++)
								for (z=0; z<D7; z++)
									array7D[t][u][v][w][x][y][z] = init; 
	
	return array7D;
}


void free_7D_farray(float ******* array7D, int D1, int D2, int D3, int D4, int D5)
{
	int t, u, v, w, x;//, y;
	if (!array7D)
		return;
	
	free(******array7D); //free space
	for (t=0; t<D1; t++)
	{
		for (u=0; u<D2; u++)
		{
			for (v=0; v<D3; v++)
			{
				for (w=0; w<D4; w++)
				{
					for (x=0; x<D5; x++)
						free(array7D[t][u][v][w][x]);
					free(array7D[t][u][v][w]);
				}
				free(array7D[t][u][v]);
			}
			free(array7D[t][u]);
		}
		free(array7D[t]);
	}
	free(array7D);
	return;
}


float *** getLowTriF(int nlays, int * lDims, float init)
{
	int tot = 0;
	int z, y, x, laySum, rowSum, n;
	z = y = x = laySum = rowSum = n = 0;
	for (z=0; z<nlays; z++)
		tot += lDims[z]*(lDims[z]+1)/2; // Gauss' method
	float *space = myalloc(tot * sizeof(*space));
	float ***lowTri = myalloc(nlays * sizeof(*lowTri));
	for (z=0; z<nlays; z++)
	{
		laySum += (z>0) ? lDims[z-1]*(lDims[z-1]+1)/2 : 0; // layer size
		lowTri[z] = myalloc(lDims[z] * sizeof(**lowTri));
		rowSum = 0;
		n = 0; // lDims[z];
		for (y=0; y<lDims[z]; y++)
		{
			rowSum += n++; //(y>0) ? n-- : 0;
			lowTri[z][y] = space + laySum + rowSum;
		}
	}
	
	/* Assign values to 3D array */
	if (fabs(init) <= EPS)
		memset(lowTri[0][0], 0, tot*sizeof(***lowTri));
	else
		for (z=0; z<nlays; z++)
			for (y=0; y<lDims[z]; y++)
				for (x=0; x<y+1; x++)
					lowTri[z][y][x] = init; //(z * nrows * ncols) + (y * ncols) + x;
	
	return lowTri;
}

inline float readLowTriF(float *** lowTri, int l, int x, int y)
{
	return (x < y) ? lowTri[l][y][x] : lowTri[l][x][y];
}


/*int printLowTriF()
{
	for (l=0; l<mp->nLayers; l++)
	{
		slen = snprintf(filename, FNAMEBUFF, "L%ddistElE.dat", l);
		assert(slen < FNAMEBUFF);
		connections_FP = myfopen(filename, "w");
		for (i=0; i<mp->vExcit[l]; i++)
		{
			for (j=0; j<=i; j++) // 0 -> i
				fprintf(connections_FP, "%f\t", distE[l][i][j]);
			fprintf(connections_FP, "\n");
		}
		fclose(connections_FP);
	}
	
}*/


float *** getUppTriF(int nlays, int * lDims, float init)
{ // Fix me!!!
	int tot = 0;
	int z, y, x, laySum, rowSum, n;
	z = y = x = laySum = rowSum = n = 0;
	for (z=0; z<nlays; z++)
		tot += lDims[z]*(lDims[z]+1)/2; // Gauss' method
	float *space = myalloc(tot * sizeof(*space));
	float ***uppTri = myalloc(nlays * sizeof(*uppTri));
	for (z=0; z<nlays; z++)
	{
		laySum += (z>0) ? lDims[z-1]*(lDims[z-1]+1)/2 : 0; // layer size
		uppTri[z] = myalloc(lDims[z] * sizeof(**uppTri));
		rowSum = 0;
		n = lDims[z];
		for (y=0; y<lDims[z]; y++)
		{
			rowSum += (y>0) ? n-- : 0;
			uppTri[z][y] = space + laySum + rowSum;
		}
	}
	
	/* Assign values to 3D array */
	if (fabs(init) <= EPS)
		memset(uppTri[0][0], 0, tot*sizeof(***uppTri));
	else
		for (z=0; z<nlays; z++)
			for (y=0; y<lDims[z]; y++)
				for (x=y; x<lDims[z]; x++) // Needs to start at 0 or memory is skipped
					uppTri[z][y][x] = init; //(z * nrows * ncols) + (y * ncols) + x;
	
	return uppTri;
}

inline float readUppTriF(float *** uppTri, int l, int x, int y)
{
	return (y < x) ? uppTri[l][y][x] : uppTri[l][x][y];
}

/*int printUppTriF(FILE * outfile, float *** uppTri, int nlays, int * lDims)
{
	int l=0;
	int i=0;
	int j=0;
	for (l=0; l<mp->nLayers; l++)
	{
		slen = snprintf(filename, FNAMEBUFF, "L%ddistElE.dat", l);
		assert(slen < FNAMEBUFF);
		connections_FP = myfopen(filename, "w");
		for (i=0; i<mp->vExcit[l]; i++)
		{
			for (j=i; j<mp->vExcit[l]; j++) // i -> Dim (j_max)
				fprintf(connections_FP, "%f\t", distE[l][i][j]);
			fprintf(connections_FP, "\n");
		}
		fclose(connections_FP);
	}
}*/

void freeTriF(float ***tri, int nlays)
{
	if (!tri)
		return;
	int z = 0;
	free(**tri);
	for (z=0; z<nlays; z++)
		free(tri[z]);
	free(tri);
	return;
}

double ** get_2D_darray(int nrows, int ncols, double init)
{
	int i, j;
	// sizeof(*array) = sizeof(int *)
	// sizeof(**array) = sizeof(int)
	
	double **array2D = myalloc(nrows * sizeof(double *));
	array2D[0] = myalloc(nrows * ncols * sizeof(double));
	if (!array2D[0])
		return NULL;
	for (i=1; i<nrows; i++) 
		array2D[i] = array2D[0] + (i * ncols);
	
	if (fabs(init) <= EPS) // (int) init == 0
		memset(array2D[0], 0, nrows*ncols*sizeof(**array2D));
	else
		for (i=0; i<nrows; i++)
			for (j=0; j<ncols; j++)
				array2D[i][j] = init; //(i * ncols) + j;
	
	return array2D;	
}

void free_2D_darray(double ** array2D)
{
	if (!array2D)
		return;
	free(array2D[0]); // Since memory is contiguous
	free(array2D);
	return;
}

double *** get_3D_darray(int nlays, int nrows, int ncols, double init)
{
	int x, y, z;
	
	// sizeof(*array) = sizeof(int **)
	// sizeof(**array) = sizeof(int *)
	// sizeof(***array) = sizeof(int)
	
	double *space = myalloc(nlays * nrows * ncols * sizeof(*space)); // sizeof(int)
	if (!space)
		return NULL;
	double ***array3D = myalloc(nlays * sizeof(*array3D)); // sizeof(int **)
	
	for (z=0; z<nlays; z++)
	{
		array3D[z] = myalloc(nrows * sizeof(**array3D)); // sizeof(int *)
		for (y=0; y<nrows; y++)
			array3D[z][y] = space + (z*ncols*nrows) + (y*ncols);
	}
	
	/* Assign values to 3D array */
	if (fabs(init) <= EPS)
		memset(array3D[0][0], 0, nlays*nrows*ncols*sizeof(***array3D));
	else
		for (z=0; z<nlays; z++)
			for (y=0; y<nrows; y++)
				for (x=0; x<ncols; x++)
					array3D[z][y][x] = init; //(z * nrows * ncols) + (y * ncols) + x;
	
	return array3D;
}


void free_3D_darray(double *** array3D, int nlays)//, int nrows)
{
	//int y,z;
	int z;
	if (!array3D)
		return;
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
