#ifndef utility_h
#define utility_h

#include <stdio.h>
#include <math.h>
#include "LCBC_data.h"
#include "LCBCmacros.h"

/* DGELS prototype */
extern "C"
{
void dgels_( char* trans, int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* work, int* lwork, int* info );
}

/* NEW Functions */
//void transposeMatrix(double *&A_transponse, double *A, int R, int C);
void LSsolve(double *A, double *b, double *&x, int m, int n, int nrhs, char trans);
void multiplyMatrices(double *&C, double *A, double *B, int R1, int C1, int R2, int C2, double s);

void print2Darray(double *arr, int s1, int s2);
void print1Darray(double *arr, int s);
void stopHere();

void print2Darray2(double *arr, int rowNum, int colmNum, int r1, int r2, int c1, int c2);
void partialPrint_3Darray(double *arr, int *Lth, int fixedInd);
double rowMaxNorm(double *Mat, int row, int numRows, int numColms);
void scaleRows(double *&Mat, double scale[], int numScaledRows, int numRows, int numColms);
void extractBlocks(double *Mat, double *&Block, int MatRowNum, int MatColmNum, int extractRows[2], int extractColms[2], int scale);
void getIDmatrix(double *&ID, int N);
void set2DArrayToZero(double **&A, int s1, int s2);
void getVarAxis(int varAxis[2], int fixedAxis);

inline void set1DArrayToZero(double *&A, int s){
    for(int i = 0; i<s; i++){
        A[i] = 0;
    }
}

inline void transposeMultiplyMatrices(double *&C, double *A, double *B, int R1, int C1, int R2, int C2, double s) {
        for (int i = 0; i < C1; i++) {
            for (int j = 0; j < C2; j++) {
                C[ind2(i,j,C1,C2)] = 0;

                for (int k = 0; k < R2; k++) {
                    C[ind2(i,j,C1,C2)] += s*A[ind2(k,i,R1,C1)]*B[ind2(k,j,R2,C2)];
                }// end of k loop
            }// end of j loop
        }// end of i loop
        return;
}// end of function

inline double dotProduct(double *v1, double *v2, int v1_param, int paramLth, int lth){
    double dotProduct = 0;
    for(int i = 0; i<lth; i++){
        dotProduct = dotProduct + v1[ind2(v1_param,i,paramLth,lth)]*v2[i];
    }// end of i loop
    return dotProduct;
}// end of dotProduct

inline double dotProduct(double *v1, double *v2, int lth){
    double dotProduct = 0;
    for(int i = 0; i<lth; i++){
        dotProduct = dotProduct + v1[i]*v2[i];
    }// end of i loop
    return dotProduct;
}// end of dotProduct

#endif /* utility_hpp */
