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


void getOtherAxes(int otherAxis[2], int fixedAxis); 
int intVectorSum(int *vec, int size);
void printIndexRange(int indexRange[3][2]);
void LSsolve(double *A, double *b, double *&x, int m, int n, int nrhs, char trans);
void multiplyMatrices(double *&C, double *A, double *B, int R1, int C1, int R2, int C2, double s);
void print2Darray(double *arr, int s1, int s2);
void print2Darray(int *arr, int s1, int s2);
void print1Darray(double *arr, int s);
void stopHere();
double rowMaxNorm(double *Mat, int row, int numRows, int numColms);
void extractBlocks(double *Mat, double *&Block, int MatRowNum, int MatColmNum, int extractRows[2], int extractColms[2], int scale);
void getIDmatrix(double *&ID, int N);
void set2DArrayToZero(double **&A, int s1, int s2);

/// \brief Set a 1D array to zero
/// \param A (output): the matrix we want to set to zero
/// \param s (input): the size of the 1D matrix
inline void set1DArrayToZero(double *&A, int s){
    for(int i = 0; i<s; i++){
        A[i] = 0;
    }
}

/// \brief Compute C = s x (A^T) x B, where s is some scalar
/// \param C (output)
/// \param A (input)
/// \param B (input)
/// \param R1 (input): the number of rows of A
/// \param C1 (input): the number of columns of A
/// \param R2 (input): the number of rows of B
/// \param C2 (input): the number of columns of B
/// \param s (input)
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

/// \brief Compute the dot product of vectors v1 and v2 where v1 may have parameters (components)
/// \param v1 (input)
/// \param v2 (input)
/// \param v1_param (input): the parameter number of v1
/// \param paramLth (input): the number of parameters v1 can take
/// \param lth (input): the length of v2 and each component of the vector v1
inline double dotProduct(double *v1, double *v2, int v1_param, int paramLth, int lth){
    double dotProduct = 0;
    for(int i = 0; i<lth; i++){
        dotProduct = dotProduct + v1[ind2(v1_param,i,paramLth,lth)]*v2[i];
    }// end of i loop
    return dotProduct;
}// end of dotProduct

/// \brief Compute the dot product of vectors v1 and v2 where v1 may have parameters (components)
/// \param v1 (input)
/// \param v2 (input)
/// \param lth (input): the length of v2 and v1
inline double dotProduct(double *v1, double *v2, int lth){
    double dotProduct = 0;
    for(int i = 0; i<lth; i++){
        dotProduct = dotProduct + v1[i]*v2[i];
    }// end of i loop
    return dotProduct;
}// end of dotProduct

#endif /* utility_hpp */
