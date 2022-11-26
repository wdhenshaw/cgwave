#include "LCBC.h"
#include "utility.h"
#include "derivatives.h"
#include <string.h>
#include "LagrangeDerivFunctions.h"

void scaleRows_new(double *&A11, double *&A12, double *&D_scaled, double *D, int numRows, int numColms1, int numColms2){
    
    for(int row = 0; row<numRows; row++){
        double scale = MAX(rowMaxNorm(A11, row, numRows, numColms1),rowMaxNorm(A12, row, numRows, numColms2)); // this finds the max norm of each row of Mat
        if(scale == 0){
            scale = 1;
        }
        /* Scale A11 row */
        for(int colm = 0; colm<numColms1; colm++){
            A11[ind2(row,colm,numRows,numColms1)] = A11[ind2(row,colm,numRows,numColms1)]/scale;
        }
        
        /* Scale A12 row */
        for(int colm = 0; colm<numColms2; colm++){
            A12[ind2(row,colm,numRows,numColms2)] = A12[ind2(row,colm,numRows,numColms2)]/scale;
        }
        
        /* Scale D row */
        for(int colm = 0; colm<numRows; colm++){
            D_scaled[ind2(row,colm,numRows,numRows)] = D[ind2(row,colm,numRows,numRows)]/scale;
        }
    }// end of r loop
}// end of the function


void Lcbc::getSideMatrix(double **&CaVec, double **&CbVec, int *eqNum, int bdryRange[3][2], int bdryNgx[3], int NU, int axis, int side){
    
    LagrangeDerivFun LagrangeDeriv = pickLagrangeDerivFun(dim, p, axis);
    
    /* parameters */
    int n = (2*p+1);
    int face = side + 2*axis;
    int compCondNum = faceParam[face].compCondNum;
    int auxiliaryEqNum = param.auxiliaryEqNum;
    int interiorEqNum = param.interiorEqNum;
    int approxEqNum = faceParam[face].approxEqNum; // number of equations treated via least squares
    int unknownVarNum = param.unknownVarNum;
    
    /* variables */
    double *A11 = new double[(approxEqNum*unknownVarNum)];
    double *A12 = new double[(approxEqNum*interiorEqNum)];
    
    /* Prepare the D matrix */
    double *D = new double[(approxEqNum*approxEqNum)];
    double *D_scaled = new double[(approxEqNum*approxEqNum)];
    getD(D, axis, side);
    
    int ghostIndRange[2] = {sideBasedValue(side, 0, (1+p)), sideBasedValue(side, (p-1), (2*p))}; // The part of the stencil where there are ghost points
    int eqNumLth[] = {n,n,dimBasedValue(dim,1,n)};   // The stencil size in each axis
    int eqNumInd[] = {p,p,dimBasedValue(dim, 0, p)}; // The stencil center: (p,p,p) used to pick out the correct row number from Ca and Cb
    
    double *Et = new double[(unknownVarNum*p)]; set1DArrayToZero(Et, (unknownVarNum*p));
    double *w = new double[(approxEqNum*p)];
    
    int ghostRow = 0;
    for(int ghostInd = ghostIndRange[0]; ghostInd<=ghostIndRange[1]; ghostInd++){
        eqNumInd[axis] = ghostInd;
        int row = eqNum[ind(eqNumInd,eqNumLth)];
        Et[ind2(row,ghostRow,unknownVarNum,p)] = 1.0;
        ghostRow++;
    }// end of ghostVal
    
    int Ind[3];
    for(Ind[2] = bdryRange[2][0]; Ind[2]<=bdryRange[2][1]; Ind[2]++){
        for(Ind[1] = bdryRange[1][0]; Ind[1]<=bdryRange[1][1]; Ind[1]++){
            for(Ind[0] = bdryRange[0][0]; Ind[0]<=bdryRange[0][1]; Ind[0]++){
                
                getBlockMatrices(A11, A12, approxEqNum, interiorEqNum, compCondNum, auxiliaryEqNum, Ind, axis, side, eqNum, NU, LagrangeDeriv);
                
                /* scale the matrices A11, A12 and D */
                scaleRows_new(A11, A12, D_scaled, D, approxEqNum, unknownVarNum, interiorEqNum);
                
                /* do a LS solve to obtain Ca and Cb */
                int copy = Ind[axis]; Ind[axis] = 0;
                int cVecInd = (cstCoef)?(0):(ind(Ind,bdryNgx)); // index of CaVec and CbVec

                LSsolve(A11, Et, w, approxEqNum, unknownVarNum, p, 'T');

                transposeMultiplyMatrices(CbVec[cVecInd], w, A12, approxEqNum, p, approxEqNum, interiorEqNum, (-1));
                transposeMultiplyMatrices(CaVec[cVecInd], w, D_scaled, approxEqNum, p, approxEqNum, approxEqNum, (1));

                Ind[axis] = copy;
                
                /* If cstCoef, break out of the loops */
                if(cstCoef){
                    for(int ax = 0; ax<3; ax++){
                        Ind[ax] = bdryRange[ax][1];
                    }// end ax loop (short for axis)
                }//end if cstCoef
                /*------------------------------------*/
                
            }// Ind[0] loop
        }// Ind[1] loop
    }// Ind[2] loop
    
    delete [] A11;
    delete [] A12;
    delete [] D_scaled;
    delete [] D;
    
    delete [] Et;
    delete [] w;
}// end of getSideMatrix

void Lcbc::getBlockMatrices(double *&A11, double *&A12, int approxEqNum, int interiorEqNum, int compCondNum, int auxiliaryEqNum, int *Ind, int axis, int side, int *eqNum, int NU, LagrangeDerivFun LagrangeDeriv){
    
    int face = side + 2*axis;
    int n = (2*p+1), MU = n;
    int totalVarNum = param.totalVarNum;
    int unknownVarNum = param.unknownVarNum;
    
    /* Prepare Variables that will hold the Lagrange derivative information */
    double Z[totalVarNum][compCondNum];
    
    /* Fill Z and Y with appropriate values */
    int LagInd[3], LagIndLth[3] = {n,n,dimBasedValue(dim, 1, n)};
    
    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim, 1, n); LagInd[2]++){
        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                
                LagrangeDeriv(Z[ind(LagInd,LagIndLth)],Ind, LagInd, LagrangeData, LagrangeData_center, coef[face].Fn, coef[face].wth, coef[face].lth, cstCoef, G.dx, axis, side, NU, memory);
                
            }// LagInd[0]
        }// LagInd[1]<
    }// LagInd[2]
    
    int row = 0;
    for(int nu = 0; nu<NU; nu++){
        for(int mu1 = 0; mu1<dimBasedValue(dim,1,MU); mu1++){
            for(int mu0 = 0; mu0<MU; mu0++){
                
                for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim,1,n); LagInd[2]++){
                    for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                        for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                            int colm = eqNum[ind(LagInd,LagIndLth)];
                            if(colm<unknownVarNum){
                                A11[ind2(row,colm,approxEqNum,unknownVarNum)] = Z[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU,MU,dimBasedValue(dim, 1, MU))];
                            }
                            else{
                                int A12colm = colm - unknownVarNum;
                                A12[ind2(row,A12colm,approxEqNum,interiorEqNum)] = Z[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU,MU,dimBasedValue(dim, 1, MU))];
                            }
                        }// LagInd[0]
                    }// LagInd[1]
                }// LagInd[2]
                
                row++;
            }// mu0 loop
        }// mu1 loop
    }// end of nu loop
    
}// end of fillMatrix_LagrangeDeriv

void Lcbc::getCornerMatrix(double **&CaVec, double **&CbVec, int *eqNum, int bdryRange[3][2], int NU1, int NU2, int varAxis, int fixedAxis[2], int side1, int side2){
    
    LagrangeDerivFun LagrangeDeriv1 = pickLagrangeDerivFun(dim, p, fixedAxis[0]);
    LagrangeDerivFun LagrangeDeriv2 = pickLagrangeDerivFun(dim, p, fixedAxis[1]);
    
    /* parameters */
    int n = (2*p+1), MU = n;
    int totalVarNum = (n*n*dimBasedValue(dim, 1, n));
    int compCondNum1 = NU1*MU*dimBasedValue(dim, 1, MU);
    int compCondNum2 = NU2*MU*dimBasedValue(dim, 1, MU);
    int auxiliaryEqNum = param.auxiliaryEqNum;
    int interiorEqNum = (p+1)*(p+1)*dimBasedValue(dim, 1,n);
    int approxEqNum = (compCondNum1 + compCondNum2 + auxiliaryEqNum);
    int unknownVarNum = totalVarNum - interiorEqNum;
    int ghostPointsNum = unknownVarNum - 2*p;
    
    /* variables */
    double *A11 = new double[(approxEqNum*unknownVarNum)];
    double *A12 = new double[(approxEqNum*interiorEqNum)];
    
    /* Prepare the D matrix */
    double *D = new double[(approxEqNum*approxEqNum)];
    double *D_scaled = new double[(approxEqNum*approxEqNum)];
    getD_corner(D, varAxis, fixedAxis, NU1, NU2);
    
    int ghostIndRange1[2] = {sideBasedValue(side1, 0, 1), sideBasedValue(side1, (2*p-1), (2*p))}; // The part of the stencil where there are ghost points
    int ghostIndRange2[2] = {sideBasedValue(side2, 0, 1), sideBasedValue(side2, (2*p-1), (2*p))}; // The part of the stencil where there are ghost points
    int eqNumLth[] = {n,n,dimBasedValue(dim,1,n)};   // The stencil size in each axis
    int eqNumInd[] = {p,p,dimBasedValue(dim, 0, p)}; // The stencil center: (p,p,p) used to pick out the correct row number from Ca and Cb
    
    double *Et = new double[(unknownVarNum*ghostPointsNum)];
    set1DArrayToZero(Et, (unknownVarNum*ghostPointsNum));
    double *w = new double[(approxEqNum*ghostPointsNum)];
    
    int ghostRow = 0; // number of relevant rows
    for(int g2 = ghostIndRange2[0]; g2<=ghostIndRange2[1]; g2++){
        for(int g1 = ghostIndRange1[0]; g1<=ghostIndRange1[1]; g1++){
            if((COND(g1, side1)||COND(g2, side2))){
                eqNumInd[fixedAxis[0]] = g1;
                eqNumInd[fixedAxis[1]] = g2;
                int k = eqNum[ind(eqNumInd,eqNumLth)];
                Et[ind2(k,ghostRow,unknownVarNum,ghostPointNum)] = 1.0;
                ghostRow++;
            }// end of if statement
        }// end of g1 loop
    }// end of g2 loop
    
    int Ind[3];
    for(Ind[2] = bdryRange[2][0]; Ind[2]<=bdryRange[2][1]; Ind[2]++){
        for(Ind[1] = bdryRange[1][0]; Ind[1]<=bdryRange[1][1]; Ind[1]++){
            for(Ind[0] = bdryRange[0][0]; Ind[0]<=bdryRange[0][1]; Ind[0]++){
                
                getBlockMatrices_edge(A11, A12, approxEqNum, interiorEqNum, unknownVarNum, compCondNum1, compCondNum2, auxiliaryEqNum, Ind, fixedAxis, side1, side2, eqNum, NU1, NU2, LagrangeDeriv1, LagrangeDeriv2);
                
                /* scale the matrices A11, A12 and D */
                scaleRows_new(A11, A12, D_scaled, D, approxEqNum, unknownVarNum, interiorEqNum);
                
                /* do a LS solve to obtain Ca and Cb */
                int cVecInd = (cstCoef)?(0):(Ind[varAxis]); // index of CaVec and CbVec // index of CaVec and CbVec

                LSsolve(A11, Et, w, approxEqNum, unknownVarNum, ghostPointsNum, 'T');
                
                transposeMultiplyMatrices(CbVec[cVecInd], w, A12, approxEqNum, ghostPointsNum, approxEqNum, interiorEqNum, (-1));
                transposeMultiplyMatrices(CaVec[cVecInd], w, D_scaled, approxEqNum, ghostPointsNum, approxEqNum, approxEqNum, (1));
                
                /* If cstCoef, break out of the loops */
                if(cstCoef){
                    for(int ax = 0; ax<3; ax++){
                        Ind[ax] = bdryRange[ax][1];
                    }// end ax loop (short for axis)
                }//end if cstCoef
                /*------------------------------------*/
                
            }// Ind[0] loop
        }// Ind[1] loop
    }// Ind[2] loop
    
    delete [] A11;
    delete [] A12;
    delete [] D_scaled;
    delete [] D;
    
    delete [] Et;
    delete [] w;
}// end of getSideMatrix

void Lcbc::getBlockMatrices_edge(double *&A11, double *&A12, int approxEqNum, int interiorEqNum, int unknownVarNum, int compCondNum1, int compCondNum2, int auxiliaryEqNum, int *Ind, int fixedAxis[2], int side1, int side2, int *eqNum, int NU1, int NU2, LagrangeDerivFun LagrangeDeriv1, LagrangeDerivFun LagrangeDeriv2){
    
    int n = (2*p+1), NU = (p+1), MU = n;
    int totalVarNum = param.totalVarNum;
    int face1 = side1 + 2*fixedAxis[0];
    int face2 = side2 + 2*fixedAxis[1];
    
    /* Prepare Variables that will hold the Lagrange derivative information */
    double Z1[totalVarNum][compCondNum1];
    double Z2[totalVarNum][compCondNum2];
    
    /* Fill Z and Y with appropriate values */
    int LagInd[3], LagIndLth[3] = {n,n,dimBasedValue(dim, 1, n)};
    
    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim, 1, n); LagInd[2]++){
        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){

                LagrangeDeriv1(Z1[ind(LagInd,LagIndLth)],Ind, LagInd, LagrangeData, LagrangeData_center, coef[face1].Fn, coef[face1].wth, coef[face1].lth, cstCoef, G.dx, fixedAxis[0], side1, NU1, memory);

                LagrangeDeriv2(Z2[ind(LagInd,LagIndLth)],Ind, LagInd, LagrangeData, LagrangeData_center, coef[face2].Fn, coef[face2].wth, coef[face2].lth, cstCoef, G.dx, fixedAxis[1], side2, NU2, memory);
                
            }// LagInd[0]
        }// LagInd[1]<
    }// LagInd[2]
    
    int row = 0;
    for(int nu = 0; nu<NU; nu++){
        for(int mu1 = 0; mu1<dimBasedValue(dim,1,MU); mu1++){
            for(int mu0 = 0; mu0<MU; mu0++){
                
                
                if(nu<NU1){
                    
                    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim,1,n); LagInd[2]++){
                        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                int colm = eqNum[ind(LagInd,LagIndLth)];
                                if(colm<unknownVarNum){
                                    A11[ind2(row,colm,approxEqNum,unknownVarNum)] = Z1[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU1,MU,dimBasedValue(dim, 1, MU))];
                                }
                                else{
                                    int A12colm = colm - unknownVarNum;
                                    A12[ind2(row,A12colm,approxEqNum,interiorEqNum)] = Z1[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU1,MU,dimBasedValue(dim, 1, MU))];
                                }
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    
                    row++;
                } // end if nu 1
                
                
                if(nu<NU2){
                    
                    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim,1,n); LagInd[2]++){
                        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                int colm = eqNum[ind(LagInd,LagIndLth)];
                                if(colm<unknownVarNum){
                                    A11[ind2(row,colm,approxEqNum,unknownVarNum)] = Z2[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU2,MU,dimBasedValue(dim, 1, MU))];
                                }
                                else{
                                    int A12colm = colm - unknownVarNum;
                                    A12[ind2(row,A12colm,approxEqNum,interiorEqNum)] = Z2[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU2,MU,dimBasedValue(dim, 1, MU))];
                                }
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    
                    row++;
                } // end if nu 1
                
            }// mu0 loop
        }// mu1 loop
    }// end of nu loop
    
}// end of fillMatrix_LagrangeDeriv

void Lcbc::getVertexMatrix(double *&CaVec, double *&CbVec, int *eqNum, int side0, int side1, int side2, int NU[3], int compCondNum[3]){
    
    LagrangeDerivFun LagrangeDeriv0 = pickLagrangeDerivFun(dim, p, 0);
    LagrangeDerivFun LagrangeDeriv1 = pickLagrangeDerivFun(dim, p, 1);
    LagrangeDerivFun LagrangeDeriv2 = pickLagrangeDerivFun(dim, p, 2);
    
    int n = (2*p+1);
    
    /* parameters */
    int totalVarNum    = (n*n*n);
    int auxiliaryEqNum = param.auxiliaryEqNum;
    int interiorEqNum  = (p+1)*(p+1)*(p+1);
    int approxEqNum    = (compCondNum[0] + compCondNum[1] + compCondNum[2] + auxiliaryEqNum);
    int totalEqNum     = approxEqNum + interiorEqNum;
    int unknownVarNum  = totalVarNum - interiorEqNum;
    int ghostPointsNum = unknownVarNum - 3*p;
    
    /* variables */
    double *A11 = new double[(approxEqNum*unknownVarNum)];
    double *A12 = new double[(approxEqNum*interiorEqNum)];
    
    /* prepare the D matrix here */
    double *D = new double[(approxEqNum*approxEqNum)];
    double *D_scaled = new double[(approxEqNum*approxEqNum)];
    getD_vertex(D, NU, compCondNum, approxEqNum);
        
    /* prepare index needed and ghost point limits */
    int eqNumLth[] = {n,n,n};   // The stencil size in each axis
    int Ind[3] = {G.indexRange[0][side0], G.indexRange[1][side1], G.indexRange[2][side2]};
    int ghostIndRange0[2] = {sideBasedValue(side0, 0, 1), sideBasedValue(side0, (2*p-1), (2*p))}; // The part of the stencil where there are ghost points
    int ghostIndRange1[2] = {sideBasedValue(side1, 0, 1), sideBasedValue(side1, (2*p-1), (2*p))};
    int ghostIndRange2[2] = {sideBasedValue(side2, 0, 1), sideBasedValue(side2, (2*p-1), (2*p))};
    
    double *Et = new double[(unknownVarNum*ghostPointsNum)];
    set1DArrayToZero(Et, (unknownVarNum*ghostPointsNum));
    double *w = new double[(approxEqNum*ghostPointsNum)];
    
    int g[3], vecNum = 0;
    for(g[2] = ghostIndRange2[0]; g[2]<=ghostIndRange2[1]; g[2]++){
        for(g[1] = ghostIndRange1[0]; g[1]<=ghostIndRange1[1]; g[1]++){
            for(g[0] = ghostIndRange0[0]; g[0]<=ghostIndRange0[1]; g[0]++){
                if((COND(g[0], side0)||COND(g[1], side1)||COND(g[2], side2))){
                    int k = eqNum[ind(g,eqNumLth)];
                    Et[ind2(k,vecNum,unknownVarNum,ghostPointNum)] = 1.0;
                    vecNum++;
                }// end of if statement
            }// end of g0 loop
        }// end of g1 loop
    }// end of g2 loop
    
    getBlockMatrices_vertex(A11, A12, Ind, side0, side1, side2, eqNum, NU, compCondNum, auxiliaryEqNum, totalEqNum, approxEqNum, unknownVarNum, LagrangeDeriv0, LagrangeDeriv1, LagrangeDeriv2);
    
    /* scale the matrices A11, A12 and D */
    scaleRows_new(A11, A12, D_scaled, D, approxEqNum, unknownVarNum, interiorEqNum);
    
    /* Work out Ca and Cb via LS */
    LSsolve(A11, Et, w, approxEqNum, unknownVarNum, ghostPointsNum, 'T');
    
    transposeMultiplyMatrices(CbVec, w, A12, approxEqNum, ghostPointsNum, approxEqNum, interiorEqNum, (-1));
    transposeMultiplyMatrices(CaVec, w, D_scaled, approxEqNum, ghostPointsNum, approxEqNum, approxEqNum, (1));
    

    delete [] A11;
    delete [] A12;
    delete [] D_scaled;
    delete [] D;
    delete [] Et;
    delete [] w;
    
}// end of getVertexMatrix

void Lcbc::getBlockMatrices_vertex(double *&A11, double *&A12, int *Ind, int side0, int side1, int side2, int *eqNum, int NU[3], int compCondNum[3], int auxiliaryEqNum, int totalEqNum, int approxEqNum, int unknownVarNum, LagrangeDerivFun LagrangeDeriv0, LagrangeDerivFun LagrangeDeriv1, LagrangeDerivFun LagrangeDeriv2){
    
    int n = (2*p+1), MU = n;
    int totalVarNum = (n*n*n);
    int NUmax = (p+1);
    int face[3] = {(side0), (side1 + 2), (side2 + 4)};
    
    /* Prepare Variables that will hold the Lagrange derivative information */
    double Z0[totalVarNum][compCondNum[0]];
    double Z1[totalVarNum][compCondNum[1]];
    double Z2[totalVarNum][compCondNum[2]];
    
    /* Fill Z and Y with appropriate values */
    int LagInd[3], LagIndLth[3] = {n,n,n};
    
    for(LagInd[2] = 0; LagInd[2]<n; LagInd[2]++){
        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){ 
                LagrangeDeriv0(Z0[ind(LagInd,LagIndLth)],Ind, LagInd, LagrangeData, LagrangeData_center, coef[face[0]].Fn, coef[face[0]].wth, coef[face[0]].lth, cstCoef, G.dx, 0, side0, NU[0], memory);
                LagrangeDeriv1(Z1[ind(LagInd,LagIndLth)],Ind, LagInd, LagrangeData, LagrangeData_center, coef[face[1]].Fn, coef[face[1]].wth, coef[face[1]].lth, cstCoef, G.dx, 1, side1, NU[1], memory);
                LagrangeDeriv2(Z2[ind(LagInd,LagIndLth)],Ind, LagInd, LagrangeData, LagrangeData_center, coef[face[2]].Fn, coef[face[2]].wth, coef[face[2]].lth, cstCoef, G.dx, 2, side2, NU[2], memory);
                
            }// LagInd[0]
        }// LagInd[1]
    }// LagInd[2]
    
    int row = 0;
    for(int nu = 0; nu<NUmax; nu++){
        for(int mu1 = 0; mu1<MU; mu1++){
            for(int mu0 = 0; mu0<MU; mu0++){
                
                if(nu<NU[0]){
                    
                    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim,1,n); LagInd[2]++){
                        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                int colm = eqNum[ind(LagInd,LagIndLth)];
                                if(colm<unknownVarNum){
                                    A11[ind2(row,colm,approxEqNum,unknownVarNum)] = Z0[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU[0],MU,MU)];
                                }
                                else{
                                    int A12colm = colm - unknownVarNum;
                                    A12[ind2(row,A12colm,approxEqNum,interiorEqNum)] = Z0[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU[0],MU,MU)];
                                }
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    
                    row++;
                    
                }// end of if nu<NU[0]
                
                if(nu<NU[1]){
                    
                    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim,1,n); LagInd[2]++){
                        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                int colm = eqNum[ind(LagInd,LagIndLth)];
                                if(colm<unknownVarNum){
                                    A11[ind2(row,colm,approxEqNum,unknownVarNum)] = Z1[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU[1],MU,MU)];
                                }
                                else{
                                    int A12colm = colm - unknownVarNum;
                                    A12[ind2(row,A12colm,approxEqNum,interiorEqNum)] = Z1[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU[1],MU,MU)];
                                }
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    
                    row++;
                    
                }// end of if nu<NU[1]
                
                if(nu<NU[2]){
                    
                    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim,1,n); LagInd[2]++){
                        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                int colm = eqNum[ind(LagInd,LagIndLth)];
                                if(colm<unknownVarNum){
                                    A11[ind2(row,colm,approxEqNum,unknownVarNum)] = Z2[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU[2],MU,MU)];
                                }
                                else{
                                    int A12colm = colm - unknownVarNum;
                                    A12[ind2(row,A12colm,approxEqNum,interiorEqNum)] = Z2[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU[2],MU,MU)];
                                }
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    
                    row++;

                }// end of if nu<NU[2]
            }// mu0 loop
        }// mu1 loop
    }// end of nu loop
}// end of fillVertexMatrix_LagrangeDeriv

