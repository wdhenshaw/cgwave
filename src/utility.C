#include "utility.h"
#include <string.h>
#include <assert.h>
//#include <malloc.h>  // *wdh* 2022/03/16
//#include <stdlib.h>  // *wdh* for exit

/// \brief Given a fixed axis. Get a vector containing the other axes.
/// \param otherAxis (output): a vector containing integers corresponding to the variable axes
/// \param fixedAxis (input)
void getOtherAxes(int otherAxis[2], int fixedAxis){
    int cnt = 0;
    for(int d  = 0; d<3; d++){
        if(d!=fixedAxis){
            otherAxis[cnt] = d;
            cnt++;
        }// end of if d statement
    }// end of d for loop
}// end of getVarAxis

/// \brief Sum the element of a vector of a given size
/// \param vec (input): the vector for which we want to sum the elements
/// \param size (input): the size of the input vector 'vec'
int intVectorSum(int *vec, int size){
    int sum = 0;
    for(int i = 0; i<size; i++){
        sum = sum + vec[i];
    }
    return sum;
}

/// \brief A function to print the index range. Used for debugging
/// \param indexRange (input): the index range object we want to print
void printIndexRange(int indexRange[3][2]){
    printf("[%d,%d]x[%d,%d]x[%d,%d]\n",indexRange[0][0],indexRange[0][1],indexRange[1][0],indexRange[1][1],indexRange[2][0],indexRange[2][1]); 
}

/// \brief Get the identity matrix of size NxN
/// \param ID (output): the identity matrix of size NxN
/// \param N (input): size of the identity matrix
void getIDmatrix(double *&ID, int N){
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            if(i == j){
                ID[ind2(i,j,N,N)] = 1.0;
            }
            else{
                ID[ind2(i,j,N,N)] = 0;
            }
        }
    }
    return;
}

/// \brief Set a 2D array to zero
/// \param A (input/output): the array that we want to set to zero
/// \param s1 (input): size of the array in the first dimension
/// \param s2 (input): size of the array in the second dimension
void set2DArrayToZero(double **&A, int s1, int s2){
    for(int i = 0; i<s1; i++){
        for(int j = 0; j<s2; j++){
            A[i][j] = 0;
        }
    }
}

/// \brief Solve Ax = b using least-squares where A is mxn and the number of right-hand-sides is nrhs
/// \param A (input)
/// \param b (input)
/// \param x (input)
/// \param m (input)
/// \param n (input)
/// \param nrhs (input): number of right-hand-side vectors
/// \param trans (input): set to 'T' to solve the transpose problem A^T x = b instead. Set to 'N' to solve the original problem
void LSsolve(double *A, double *b, double *&x, int m, int n, int nrhs, char trans){
    int M = m;
    int N = n;
    int NRHS = nrhs;
    int b_rows = (trans == 'N')?(m):(n);
    
    int LDA = m; // length of the first dimension of A
    int maxDim = MAX(m,n);
    int LDB = MAX(1,maxDim);
    int info, lwork;
    
    double *Q = new double[(LDA*N)];    // vector to store A values
    double *B = new double[(LDB*NRHS)]; // vector to store b
    
    /* copy A and b into new locations */
    for(int i = 0; i<(LDA*N); i++){
        Q[i] = A[i];
    }
    
    for(int j = 0; j<nrhs; j++){
        for(int i = 0; i<b_rows; i++){
            B[ind2(i,j,LDB,nrhs)] = b[ind2(i,j,b_rows,nrhs)];
        }
        for(int i = b_rows; i<LDB; i++){
            B[ind2(i,j,LDB,nrhs)] = 0.0;
        }
    }

    double wkopt;
    double* work;
    lwork = -1;
    dgels_(&trans, &M, &N, &NRHS, Q, &LDA, B, &LDB, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work  = (double*)malloc( lwork*sizeof(double) );

    /* Solve the equations A*X = B */
    dgels_(&trans, &M, &N, &NRHS, Q, &LDA, B, &LDB, work, &lwork, &info );
    /* the function above returns A as its QR factorization and b as the solution of Ax = b */
    
    /* double check that the matrix is not rank deficient */
    if( info > 0 ) {
        printf( "The diagonal element %i of the triangular factor ", info );
        printf( "of A is zero, so that A does not have full rank;\n" );
        printf( "the least squares solution could not be computed.\n" );
        exit( 1 );
    }
    
    /* copy the values of bc into x to avoid confusion */
    
    int rowsOfx = (trans == 'N')?(n):(m);
    for(int i = 0; i<rowsOfx ; i++){
        for(int j = 0; j<nrhs; j++){
            x[ind2(i,j,rowsOfx ,nrhs)] = B[ind2(i,j,LDB,nrhs)];
        }
    }
    
    free( (void*)work );
    delete [] Q;
    delete [] B;
}

/// \brief A function to extract block matrices from a bigger matrix Mat
/// \param Mat (input): the original matrix that we'd like to extract blocks from
/// \param Block (output): the block matrix that we want to extract
/// \param MatRowNum (input): the number of rows of Mat
/// \param MatColmNum (input): the number of columns of Mat
/// \param extractRows (input): a vector containing the range of rows that we want to extract
/// \param extractColms (input): a vector containing the range of columns that we want to extract
/// \param scale (input): multiply the extracted block by a constant 'scale'
void extractBlocks(double *Mat, double *&Block, int MatRowNum, int MatColmNum, int extractRows[2], int extractColms[2], int scale){
    int BlockRowNum  = extractRows[1] - extractRows[0];
    int BlockColmNum = extractColms[1] - extractColms[0];
    
    int row = 0; int colm = 0;
    for(int MatRow=extractRows[0]; MatRow<extractRows[1]; MatRow++ ){
        for(int MatColm = extractColms[0]; MatColm<extractColms[1]; MatColm++){
            Block[ind2(row,colm,BlockRowNum,BlockColmNum)] = scale*Mat[ind2(MatRow,MatColm,MatRowNum,MatColmNum)];
            colm++;
        }
        colm = 0;
        row++;
    }
}

/// \brief Find the maximum norm of each row of a matrix Mat
/// \param Mat (input)
/// \param row (input): the row of the matrix Mat we want to find the maximum norm of
/// \param numRows (input): the number of rows of the matrix Mat
/// \param numColms (input): the number of columns of the matrix Mat
double rowMaxNorm(double *Mat, int row, int numRows, int numColms){
    /* this functions finds the maximum of the absolute values of rows of a matrix Mat */
    double max = Mat[ind2(row,0,numRows,numColms)];

    for(int colm = 1; colm<numColms; colm++){
        max = MAX(max,fabs(Mat[ind2(row,colm,numRows,numColms)]));
    } // end of j loop

    return max;
}// end of function

/// \brief Print a 2D array of doubles
/// \param arr (input): the array we want to print
/// \param s1 (input): the size of the array arr in the first dimension
/// \param s2 (input): the size of the array arr in the second dimension
void print2Darray(double *arr, int s1, int s2){
    for(int i = 0; i<s1; i++){
        for(int j = 0; j<s2; j++){
            printf("%1.4f|",arr[ind2(i,j,s1,s2)]);
        }
        printf("\n");
    }
    printf("\n");
}

/// \brief Print a 2D array of integers
/// \param arr (input): the array we want to print
/// \param s1 (input): the size of the array arr in the first dimension
/// \param s2 (input): the size of the array arr in the second dimension
void print2Darray(int *arr, int s1, int s2){
    for(int i = 0; i<s1; i++){
        for(int j = 0; j<s2; j++){
            printf("%d|",arr[ind2(i,j,s1,s2)]);
        }
        printf("\n");
    }
    printf("\n");
}

/// \brief Print a 1D array of doubles
/// \param arr (input): the array we want to print
/// \param s (input): the size of the array arr
void print1Darray(double *arr, int s){
    for(int i = 0; i<s; i++){
        printf("%1.4f\n",arr[i]);
    }
    printf("\n");
}

/// \brief A debugging function used to print the message "STOP HERE"
void stopHere(){
    printf("\n STOP HERE \n"); 
}

