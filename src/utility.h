#ifndef utility_h
#define utility_h

#include <stdio.h>

/* DGELS prototype */
extern "C"
{
void dgels_( char* trans, int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* work, int* lwork, int* info );
}
void print2Darray(double *arr, int s1, int s2);
void print2Darray2(double *arr, int rowNum, int colmNum, int r1, int r2, int c1, int c2);
void partialPrint_3Darray(double *arr, int *Lth, int fixedInd);
double rowMaxNorm(double *Mat, int row, int numRows, int numColms);
void scaleRows(double *Mat, double *scale, int numScaledRows, int numRows, int numColms);
void extractBlocks(double *Mat, double *Block, int MatRowNum, int MatColmNum, int extractRows[2], int extractColms[2], int scale);
void LSsolve(double *A, double *b, double *x, int m, int n, int nrhs); 
void getIDmatrix(double *ID, int N);
void set2DArrayToZero(double **A, int s1, int s2);
void set1DArrayToZero(double *A, int s);
void getVarAxis(int varAxis[2], int fixedAxis);

double dotProduct(double *v1, double *v2, int lth);
#endif /* utility_hpp */
