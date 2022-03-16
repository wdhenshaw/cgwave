#include "utility.h"
#include <math.h>
#include <malloc.h>  // *wdh* 2022/03/16
#include <stdlib.h>  // *wdh* for exit

#define ind2(i,j,n1,n2) (((j)*(n1))+(i))
#define ind(I,N) (I[0] + N[0]*(I[1]+ I[2]*N[1]))
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y));
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y));

double dotProduct(double *v1, double *v2, int lth){
    double dotProduct = 0;
    for(int i = 0; i<lth; i++){
        dotProduct = dotProduct + v1[i]*v2[i];
    }// end of i loop
    return dotProduct;
}// end of dotProduct

void getVarAxis(int varAxis[2], int fixedAxis){
    int cnt = 0;
    for(int d  = 0; d<3; d++){
        if(d!=fixedAxis){
            varAxis[cnt] = d;
            cnt++;
        }// end of if d statement
    }// end of d for loop
}// end of getVarAxis

void getIDmatrix(double *ID, int N){
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

void set2DArrayToZero(double **A, int s1, int s2){
    for(int i = 0; i<s1; i++){
        for(int j = 0; j<s2; j++){
            A[i][j] = 0;
        }
    }
}

void set1DArrayToZero(double *A, int s){
    for(int i = 0; i<s; i++){
            A[i] = 0;
    }
}

void LSsolve(double *A, double *b, double *x, int m, int n, int nrhs){
    
    /* copy A and b into new locations */
    double *Q = new double[(m*n)];
    for(int i = 0; i<(m*n); i++){
        Q[i] = A[i];
    }
    
    /* allocate work space */
    int lda = m, ldb = m, info, lwork;
    double wkopt;
    double* work;
    lwork = -1;
    char prop[] = "No transpose";
    dgels_(prop, &m, &n, &nrhs, Q, &lda, b, &ldb, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work  = (double*)malloc( lwork*sizeof(double) );
    
    /* Solve the equations A*X = B */
    dgels_(prop, &m, &n, &nrhs, Q, &lda, b, &ldb, work, &lwork, &info );
    /* the function above returns A as its QR factorization and b as the solution of Ax = b */
    
    /* double check that the matrix is not rank deficient */
    if( info > 0 ) {
        printf( "The diagonal element %i of the triangular factor ", info );
        printf( "of A is zero, so that A does not have full rank;\n" );
        printf( "the least squares solution could not be computed.\n" );
        exit( 1 );
         
    }
    
    /* copy the values of bc into x to avoid confusion */
    for(int i = 0; i<n; i++){
        for(int j = 0; j<nrhs; j++){
            x[ind2(i,j,n,nrhs)] = b[ind2(i,j,m,nrhs)];
        }
    }
    
    free( (void*)work );
    delete [] Q; 
}

void extractBlocks(double *Mat, double *Block, int MatRowNum, int MatColmNum, int extractRows[2], int extractColms[2], int scale){
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

void scaleRows(double *Mat, double *scale, int numScaledRows, int numRows, int numColms){
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
            printf("%1.5f|",arr[ind2(i,j,s1,s2)]);
        }
        printf("\n");
    }
    printf("\n"); 
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
