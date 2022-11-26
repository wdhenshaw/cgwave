#include "utility.h"
#include <string.h>
#include <assert.h>
#include <malloc.h>  // *wdh* 2022/03/16
#include <stdlib.h>  // *wdh* for exit

//void multiplyMatrices(double *&C, double *A, double *B, int R1, int C1, int R2, int C2, double s) {
////    if (C1 != R2) {
////        printf("Matrices are incompatible for multiplication\n");
////        exit(EXIT_FAILURE);
////    }
////    else{
//        for (int i = 0; i < R1; i++) {
//            for (int j = 0; j < C2; j++) {
//                C[ind2(i,j,R1,C2)] = 0;
//
//                for (int k = 0; k < R2; k++) {
//                    C[ind2(i,j,R1,C2)] += s*A[ind2(i,k,R1,C1)]*B[ind2(k,j,R2,C2)];
//                }// end of k loop
//            }// end of j loop
//        }// end of i loop
//        return;
////    }
//}// end of function

void getVarAxis(int varAxis[2], int fixedAxis){
    int cnt = 0;
    for(int d  = 0; d<3; d++){
        if(d!=fixedAxis){
            varAxis[cnt] = d;
            cnt++;
        }// end of if d statement
    }// end of d for loop
}// end of getVarAxis

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

void set2DArrayToZero(double **&A, int s1, int s2){
    for(int i = 0; i<s1; i++){
        for(int j = 0; j<s2; j++){
            A[i][j] = 0;
        }
    }
}

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
    
//    print2Darray(B, LDB, NRHS);

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

void scaleRows(double *&Mat, double scale[], int numScaledRows, int numRows, int numColms){
    /* this function scales the rows of the matrix Mat by their max norm and saves the scaling factor in scale */
    
    int cnt = 0;
    for(int row = 0; row<numScaledRows; row++){
        scale[cnt] = rowMaxNorm(Mat, row, numRows, numColms); // this finds the max norm of each row of Mat
        if(scale[cnt] == 0){
            scale[cnt] = 1;
        }
        for(int colm = 0; colm<numColms; colm++){
            Mat[ind2(row,colm,numRows,numColms)] = Mat[ind2(row,colm,numRows,numColms)]/scale[cnt];
        }// end of c loop
        cnt = cnt + 1;
    }// end of r loop
}// end of the function

double rowMaxNorm(double *Mat, int row, int numRows, int numColms){
    /* this functions finds the maximum of the absolute values of rows of a matrix Mat */
    double max = Mat[ind2(row,0,numRows,numColms)];

    for(int colm = 1; colm<numColms; colm++){
        max = MAX(max,fabs(Mat[ind2(row,colm,numRows,numColms)]));
    } // end of j loop

    return max;
}// end of function

void print2Darray(double *arr, int s1, int s2){
    for(int i = 0; i<s1; i++){
        for(int j = 0; j<s2; j++){
            printf("%1.4f|",arr[ind2(i,j,s1,s2)]);
        }
        printf("\n");
    }
    printf("\n");
}

void print1Darray(double *arr, int s){
    for(int i = 0; i<s; i++){
        printf("%1.4f\n",arr[i]);
    }
    printf("\n");
}

void stopHere(){
    printf("\n STOP HERE \n"); 
}

void partialPrint_3Darray(double *arr, int *Lth, int fixedInd){
    int i[3];
    i[2] = fixedInd;
    for(i[0] = 0; i[0]<Lth[0]; i[0]++){
        for(i[1] = 0; i[1]<Lth[1]; i[1]++){
            printf("%1.2f|",arr[ind(i, Lth)]);
        }// end of i[1]
        printf("\n");
    }// end of i[0]
    
    printf("\n");
}// end of fn

void print2Darray2(double *arr, int rowNum, int colmNum, int r1, int r2, int c1, int c2){
    for(int i = r1; i<r2; i++){
        for(int j = c1; j<c2; j++){
            printf("%1.2f|",arr[ind2(i,j,rowNum,colmNum)]);
        }// end of j loop
        printf("\n");
    }// end of i loop
    printf("\n");
}


