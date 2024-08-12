#ifndef HELPER_H_
#define HELPER_H_

#include <stdio.h>
#include <math.h>

#define PI M_PI
#define BohrToAng 0.529177210544
#define DIE(message) { printf(message); printf("\nAt line %d in file %s\n", __LINE__, __FILE__); fflush(stdout); exit(1); }
#define PARALLEL_FOR for
#define BOX_LENGTH 4.0

void dgemv(const char TransA, const int m, const int n, const double alpha, const double* a, const int lda, const double* x, const int incx, const double beta, double* y, const int incy);
void dgemm(const char TransA, const char TransB, const int M, const int N, const int K, const  double alpha, const  double* A, const int lda, const  double* B, const int ldb, const double beta, double* C, const int ldc);

#include <cuda_runtime.h>
#define CUERR { cudaError_t err = cudaGetLastError(); if (err != cudaSuccess) { printf("CUDA Error: %s\n", cudaGetErrorString(err)); fflush(stdout); exit(1); } }

#endif