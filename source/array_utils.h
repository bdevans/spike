/*
 *  array_utils.h
 *  Spike
 *
 *  Created by Ben Evans on 04/02/2009.
 *  Copyright 2009 University of Oxford. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
//#include <errno.h>

#include "utils.h"
//#include "globals.h"

#ifndef EPS
#define EPS		1.2e-7
#endif
//const double eps = 1.0e-8;

extern void normalize(double *vec, int size);
extern double ndp(double *vec_x, double *vec_y, int size);
extern void print_irow(FILE * file_ptr, int * data, int ncols);
extern void print_iarray(FILE *file_ptr, int **data, int nrows, int ncols);
extern void print_frow(FILE * file_ptr, float * data, int ncols);
extern void print_farray(FILE *file_ptr, float **data, int nrows, int ncols);
extern void print_darray(FILE *file_ptr, double **data, int nrows, int ncols);
extern void ** array2d(size_t rows, size_t cols, size_t value_size);
extern int ** get_2D_iarray(int nrows, int ncols, int init);
extern void free_2D_iarray(int ** array2D);
extern int *** get_3D_iarray(int nlays, int nrows, int ncols, int init);
extern void free_3D_iarray(int *** array3D, int nlays);
extern float ** get_2D_farray(int nrows, int ncols, float init);
extern void free_2D_farray(float ** array2D);
extern float *** get_3D_farray(int nlays, int nrows, int ncols, float init);
extern void free_3D_farray(float *** array3D, int nlays);
extern float **** get4Dfarray(int D1, int D2, int D3, int D4, float init);
extern void free4Dfarray(float **** array4D, int D1, int D2);
extern float ******* get_7D_farray(int D1, int D2, int D3, int D4, int D5, int D6, int D7, float init);
extern void free_7D_farray(float ******* array7D, int D1, int D2, int D3, int D4, int D5);
extern float *** getLowTriF(int nlays, int * lDims, float init);
extern inline float readLowTriF(float *** lowTri, int l, int x, int y);
extern float *** getUppTriF(int nlays, int * lDims, float init);
extern inline float readUppTriF(float *** uppTri, int l, int x, int y);
extern void freeTriF(float ***tri, int nlays);
extern double ** get_2D_darray(int nrows, int ncols, double init);
extern void free_2D_darray(double ** array2D);
extern double *** get_3D_darray(int nlays, int nrows, int ncols, double init);
extern void free_3D_darray(double *** array3D, int nlays);
