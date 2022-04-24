#include "utility.h"
#include <math.h>
#include <malloc.h>  // *wdh* 2022/03/16
#include <stdlib.h>  // *wdh* for exit
#include <string.h>
#include "LCBCmacros.h"

double extrapolate(double *gridFn, int index[3], int lth[3], int axis, int sign, int order){
    
    double E = 0;
    
    if(order == 2){
        int i1[3]; memcpy(i1, index, sizeof(i1)); i1[axis] = i1[axis] + sign*1;
        int i2[3]; memcpy(i2, index, sizeof(i1)); i2[axis] = i2[axis] + sign*2;
        int i3[3]; memcpy(i3, index, sizeof(i1)); i3[axis] = i3[axis] + sign*3;
        
        E = 3.0*gridFn[ind(i1,lth)]
        -3.0*gridFn[ind(i2,lth)]
        +1.0*gridFn[ind(i3,lth)];
    }
    
    else if(order == 4){
        
        int i1[3]; memcpy(i1, index, sizeof(i1)); i1[axis] = i1[axis] + sign*1;
        int i2[3]; memcpy(i2, index, sizeof(i1)); i2[axis] = i2[axis] + sign*2;
        int i3[3]; memcpy(i3, index, sizeof(i1)); i3[axis] = i3[axis] + sign*3;
        int i4[3]; memcpy(i4, index, sizeof(i1)); i4[axis] = i4[axis] + sign*4;
        int i5[3]; memcpy(i5, index, sizeof(i1)); i5[axis] = i5[axis] + sign*5;
        
        E =  5.0*gridFn[ind(i1,lth)]
        -10.0*gridFn[ind(i2,lth)]
        +10.0*gridFn[ind(i3,lth)]
        -5.0*gridFn[ind(i4,lth)]
        +1.0*gridFn[ind(i5,lth)];
    }
    
    else if(order == 6){
        int i1[3]; memcpy(i1, index, sizeof(i1)); i1[axis] = i1[axis] + sign*1;
        int i2[3]; memcpy(i2, index, sizeof(i1)); i2[axis] = i2[axis] + sign*2;
        int i3[3]; memcpy(i3, index, sizeof(i1)); i3[axis] = i3[axis] + sign*3;
        int i4[3]; memcpy(i4, index, sizeof(i1)); i4[axis] = i4[axis] + sign*4;
        int i5[3]; memcpy(i5, index, sizeof(i1)); i5[axis] = i5[axis] + sign*5;
        int i6[3]; memcpy(i6, index, sizeof(i1)); i6[axis] = i6[axis] + sign*6;
        int i7[3]; memcpy(i7, index, sizeof(i1)); i7[axis] = i7[axis] + sign*7;
        
        E = 7.0*gridFn[ind(i1,lth)]
        -21.0*gridFn[ind(i2,lth)]
        +35.0*gridFn[ind(i3,lth)]
        -35.0*gridFn[ind(i4,lth)]
        +21.0*gridFn[ind(i5,lth)]
        -7.0*gridFn[ind(i6,lth)]
        +1.0*gridFn[ind(i7,lth)];
        
    }else{
        printf("Extrapolation order %d is not supported (extrapolation/utility.C)\n",order);
        exit(-1);
    }
    
    return E;
}

void getForcingGridLth(int indexRange[3][2], int lth[3], int wth[3], int bdryRange[3][2], int fixedAxis, int fixedSide, int dim, int p){
    for(int axis = 0; axis<3; axis++){
        wth[axis] = ((axis==fixedAxis)?(p):(2*p));
        if(axis==dim){
            wth[axis] = 0;
        }
        for(int side = 0; side<2; side++){
            if(axis == fixedAxis){
                bdryRange[axis][side] = indexRange[fixedAxis][fixedSide];
            }
            else{
                bdryRange[axis][side] = indexRange[axis][side];
            }// end of if axis
        }// end of side
        int bdryRangeLth = bdryRange[axis][1] - bdryRange[axis][0] + 1;
        lth[axis] = (2*wth[axis] + bdryRangeLth);
    }// end of axis
}

void fixCoefGridFunctions(double **&coef, int indexRange[3][2], int dim, int p, int *faceEval){
    
    int faceNum = 2*dim;
    int maxCoefNum = ((dim+1)*(dim + 2))/2;
    
    for(int axis = 0; axis<dim; axis++){
        int varAxis[2]; getVarAxis(varAxis, axis);
        int v0 = varAxis[0], v1 = varAxis[1];
        
        for(int side = 0; side<2; side++){
            int face = (side + 2*axis);
            if(faceEval[face]>0){
                int bdryRange[3][2],lth[3], wth[3]; getForcingGridLth(indexRange, lth, wth, bdryRange, axis, side, dim, p);
                
                int i[3];
                for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
                    
                    /* LOWER y axis over partial z axis */
                    for(i[v0] = (p-1); i[v0]>=0; i[v0]--){
                        for(i[v1] = dimBasedValue(dim, 0, p); i[v1]<dimBasedValue(dim, 1, (lth[v1]-p)); i[v1]++){
                            for(i[axis] = 0; i[axis]<lth[axis]; i[axis]++){
                                
                                coef[ind2(face,coefNum,faceNum,maxCoefNum)][ind(i,lth)] = extrapolate(coef[ind2(face,coefNum,faceNum,maxCoefNum)], i, lth, v0, 1, (2*p));
                                
                            }// end of i[0] loop
                        }// end of i[1] loop
                    }// end of i[2] loop
                    
                    /* Upper y axis over partial z axis */
                    for(i[v0] = (lth[v0]-p); i[v0]<lth[v0]; i[v0]++){
                        for(i[v1] = dimBasedValue(dim, 0, p); i[v1]<dimBasedValue(dim, 1, (lth[v1]-p)); i[v1]++){
                            for(i[axis] = 0; i[axis]<lth[axis]; i[axis]++){
                                
                                coef[ind2(face,coefNum,faceNum,maxCoefNum)][ind(i,lth)] = extrapolate(coef[ind2(face,coefNum,faceNum,maxCoefNum)], i, lth, v0, (-1), (2*p));
                                
                            }// end of i[0] loop
                        }// end of i[1] loop
                    }// end of i[2] loop
                    
                    if(dim == 3){
                        /* Lower z */
                        for(i[v1] = (p-1); i[v1]>=0; i[v1]--){
                            for(i[v0] = 0; i[v0]<lth[v0]; i[v0]++){
                                for(i[axis] = 0; i[axis]<lth[axis]; i[axis]++){
                                    
                                    coef[ind2(face,coefNum,faceNum,maxCoefNum)][ind(i,lth)] = extrapolate(coef[ind2(face,coefNum,faceNum,maxCoefNum)], i, lth, v1, 1, (2*p));
                                    
                                }// end of i[axis]
                            }// end of i[v0]
                        }// end of i[v1]
                        
                        /* Upper z */
                        for(i[v1] = (lth[v1]-p); i[v1]< lth[v1]; i[v1]++){
                            for(i[v0] = 0; i[v0]<lth[v0]; i[v0]++){
                                for(i[axis] = 0; i[axis]<lth[axis]; i[axis]++){
                                    
                                    coef[ind2(face,coefNum,faceNum,maxCoefNum)][ind(i,lth)] = extrapolate(coef[ind2(face,coefNum,faceNum,maxCoefNum)], i, lth, v1, (-1), (2*p));
                                    
                                }// end of i[axis]
                            }// end of i[v0]
                        }// end of i[v1]
                    }// end of if dim statment
                    
                }// end of nu loop
            }// end of if faceEval
        }// end of side loop
    }// end of axis loop
    
    
}// getBdryGridFunctions

void fixForcingGridFunctions(double **&fn, int indexRange[3][2], int dim, int p, int *faceEval){
    
    int faceNum = 2*dim;
    int NU= p;
    
    for(int axis = 0; axis<dim; axis++){
        int varAxis[2]; getVarAxis(varAxis, axis);
        int v0 = varAxis[0], v1 = varAxis[1];
        
        for(int side = 0; side<2; side++){
            int face = (side + 2*axis);
            if(faceEval[face]>0){
                int bdryRange[3][2],lth[3], wth[3]; getForcingGridLth(indexRange, lth, wth, bdryRange, axis, side, dim, p);
                
                int i[3];
                for(int nu = 0; nu<NU; nu++){
                    
                    /* LOWER y axis over partial z axis */
                    for(i[v0] = (p-1); i[v0]>=0; i[v0]--){
                        for(i[v1] = dimBasedValue(dim, 0, p); i[v1]<dimBasedValue(dim, 1, (lth[v1]-p)); i[v1]++){
                            for(i[axis] = 0; i[axis]<lth[axis]; i[axis]++){
                                
                                fn[ind2(face,nu,faceNum,NU)][ind(i,lth)] = extrapolate(fn[ind2(face,nu,faceNum,NU)], i, lth, v0, 1, (2*p));
                                
                            }// end of i[0] loop
                        }// end of i[1] loop
                    }// end of i[2] loop
                    
                    /* Upper y axis over partial z axis */
                    for(i[v0] = (lth[v0]-p); i[v0]<lth[v0]; i[v0]++){
                        for(i[v1] = dimBasedValue(dim, 0, p); i[v1]<dimBasedValue(dim, 1, (lth[v1]-p)); i[v1]++){
                            for(i[axis] = 0; i[axis]<lth[axis]; i[axis]++){
                                
                                fn[ind2(face,nu,faceNum,NU)][ind(i,lth)] = extrapolate(fn[ind2(face,nu,faceNum,NU)], i, lth, v0, (-1), (2*p));
                                
                            }// end of i[0] loop
                        }// end of i[1] loop
                    }// end of i[2] loop
                    
                    if(dim == 3){
                        /* Lower z */
                        for(i[v1] = (p-1); i[v1]>=0; i[v1]--){
                            for(i[v0] = 0; i[v0]<lth[v0]; i[v0]++){
                                for(i[axis] = 0; i[axis]<lth[axis]; i[axis]++){
                                    
                                    fn[ind2(face,nu,faceNum,NU)][ind(i,lth)] = extrapolate(fn[ind2(face,nu,faceNum,NU)], i, lth, v1, 1, (2*p));
                                    
                                }// end of i[axis]
                            }// end of i[v0]
                        }// end of i[v1]
                        
                        /* Upper z */
                        for(i[v1] = (lth[v1]-p); i[v1]< lth[v1]; i[v1]++){
                            for(i[v0] = 0; i[v0]<lth[v0]; i[v0]++){
                                for(i[axis] = 0; i[axis]<lth[axis]; i[axis]++){
                                    
                                    fn[ind2(face,nu,faceNum,NU)][ind(i,lth)] = extrapolate(fn[ind2(face,nu,faceNum,NU)], i, lth, v1, (-1), (2*p));
                                    
                                }// end of i[axis]
                            }// end of i[v0]
                        }// end of i[v1]
                    }// end of if dim statment
                    
                }// end of nu loop
            }// end of if faceEval
        }// end of side loop
    }// end of axis loop
    
}// getBdryGridFunctions

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

void set1DArrayToZero(double *&A, int s){
    for(int i = 0; i<s; i++){
        A[i] = 0;
    }
}

void LSsolve(double *A, double *b, double *&x, int m, int n, int nrhs){
    
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
