
#include <cblas.h>

void dgemv(const char TransA, const int m, const int n, const double alpha, const double* a, const int lda, const double* x, const int incx, const double beta, double* y, const int incy)
{
    const CBLAS_TRANSPOSE TransA_converted = (TransA == 't') ? CblasTrans : CblasNoTrans;
    cblas_dgemv(CblasColMajor, TransA_converted, m, n, alpha, a, lda, x, incx, beta, y, incy);
}

void dgemm(const char TransA, const char TransB, const int M, const int N, const int K, const  double alpha, const  double* A, const int lda, const  double* B, const int ldb, const double beta, double* C, const int ldc)
{
    const CBLAS_TRANSPOSE TransA_converted = (TransA == 't') ? CblasTrans : CblasNoTrans;
    const CBLAS_TRANSPOSE TransB_converted = (TransB == 't') ? CblasTrans : CblasNoTrans;
    cblas_dgemm(CblasColMajor, TransA_converted, TransB_converted, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}
