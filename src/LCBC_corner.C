#include "LCBC.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utility.h"
#include "LCBCmacros.h"

void getFixedAxis(int fixedAxis[2], int varAxis);

void Lcbc::updateEdgeGhostPeriodic(double *&unp1, int side1, int side2, int varAxis, int fixedAxis[2]){
    int fixedSides[2] = {side1, side2};
    int periodicAxis = -1, perAxisInd = -1;
    for(int i = 0; i<2; i++){
        int edgeFace = fixedSides[i] + 2*fixedAxis[i];
        if(faceEval[edgeFace] == -1){
            periodicAxis = fixedAxis[i];
            perAxisInd = i;
        }// end of faceEval
    }// end of i
    
    assert(periodicAxis>=0 && perAxisInd>=0);
    
    int bdryRange[3][2]; getBdryRange_corner(bdryRange, fixedAxis, fixedSides, varAxis, 0, 0);
    int perIndexRange[3][2];
    perIndexRange[fixedAxis[0]][0] = bdryRange[fixedAxis[0]][0] + sideBasedValue(side1, (-p), 1);
    perIndexRange[fixedAxis[0]][1] = bdryRange[fixedAxis[0]][1] + sideBasedValue(side1, (-1), p);
    perIndexRange[fixedAxis[1]][0] = bdryRange[fixedAxis[1]][0] + sideBasedValue(side2, (-p), 1);
    perIndexRange[fixedAxis[1]][1] = bdryRange[fixedAxis[1]][0] + sideBasedValue(side2, (-1), p);
    perIndexRange[varAxis][0]      = bdryRange[varAxis][0];
    perIndexRange[varAxis][1]      = bdryRange[varAxis][1];

    
    int i[3];
    for(i[2] = perIndexRange[2][0]; i[2]<=perIndexRange[2][1]; i[2]++){
        for(i[1] = perIndexRange[1][0]; i[1]<=perIndexRange[1][1]; i[1]++){
            for(i[0] = perIndexRange[0][0]; i[0]<=perIndexRange[0][1]; i[0]++){
                int I[3]; memcpy(I, i, sizeof(I));
          
                I[periodicAxis] = I[periodicAxis] + sideBasedValue(fixedSides[perAxisInd], G.Nx[periodicAxis], (-G.Nx[periodicAxis]));
                
                unp1[solInd(i, G.Ngx)] = unp1[solInd(I, G.Ngx)];
                
            }// end of i[0]
        }// end of i[1]
    }// end of i[2]
}// end of updateEdgeGhostPeriodic

void Lcbc::updateEdgeGhost(double **R, double *&unp1, double t, double dt, int side1, int side2, int varAxis){

    int fixedAxis[2]; getFixedAxis(fixedAxis, varAxis);
    int face1 = side1 + 2*fixedAxis[0];
    int face2 = side2 + 2*fixedAxis[1];

    int NU1 = faceParam[face1].NU;
    int NU2 = faceParam[face2].NU;
    int auxiliaryEqNum = param.auxiliaryEqNum;
    int compCondNum1   = faceParam[face1].compCondNum;
    int compCondNum2   = faceParam[face2].compCondNum;
    int approxEqNum    = compCondNum1 + compCondNum2 + auxiliaryEqNum;
    
    int bdryPointsNum = dimBasedValue(dim, 1, G.Ngx[varAxis]);
    double **Rv = new double*[bdryPointsNum];
    
    for(int bdryPt = 0; bdryPt<bdryPointsNum; bdryPt++){
        Rv[bdryPt] = new double[approxEqNum];
    }
    
    int corner = ind3(side1,side2,dimBasedValue(dim, 0, varAxis),2,2,3);
    
    /* prepare the matrices for the corner */
    if(CornerMat[corner].flag == false){
        prepCornerMatrix(side1, side2, varAxis, approxEqNum, NU1, NU2);
    }
    
    /* prepare the RHS vector for the corner */
    prepDataVec_corner(R, Rv, t, dt, side1, side2, varAxis, approxEqNum, NU1, NU2);
    
    /* evaluate the ghost points */
    getCornerGhost(unp1, Rv, CornerMat[corner].CaVec, CornerMat[corner].CbVec, CornerMat[corner].eqNum, side1, side2, varAxis, approxEqNum);
    
    /* Free any allocated variables */
    for(int bdryPoint = 0; bdryPoint<bdryPointsNum; bdryPoint++)
        delete [] Rv[bdryPoint];
    
    delete [] Rv;
    
}// end of updateGhost function

void Lcbc::prepCornerMatrix(int side1, int side2, int varAxis, int approxEqNum, int NU1, int NU2){
    int corner = ind3(side1,side2,dimBasedValue(dim, 0, varAxis),2,2,3);
    int bdryNg = (cstCoef)?(1):(dimBasedValue(dim, 1, G.Ngx[varAxis]));
    CornerMat[corner].flag = true;
    
    int fixedSides[] = {side1, side2};
    int fixedAxis[2]; getFixedAxis(fixedAxis, varAxis);
    int bdryRange[3][2]; getBdryRange_corner(bdryRange, fixedAxis, fixedSides, varAxis, (p), (-p));
    
    int n = (2*p+1);
    int totalVarNum = n*n*dimBasedValue(dim, 1, n);
    int knownVarNum   = (p+1)*(p+1)*dimBasedValue(dim, 1, n);
    int unknownVarNum = (totalVarNum - knownVarNum);
    int ghostPointsNum = unknownVarNum - 2*p;
    
    CornerMat[corner].eqNum = new int[totalVarNum];
    CornerMat[corner].CaVec = new double**[bdryNg];
    CornerMat[corner].CbVec = new double**[bdryNg];
    
    for(int bdryPoint = 0; bdryPoint<bdryNg; bdryPoint++){
        CornerMat[corner].CaVec[bdryPoint] = new double*[ghostPointsNum];
        CornerMat[corner].CbVec[bdryPoint] = new double*[ghostPointsNum];
        
        for(int gp = 0; gp<ghostPointsNum; gp++){
            CornerMat[corner].CaVec[bdryPoint][gp] = new double[approxEqNum];
            CornerMat[corner].CbVec[bdryPoint][gp] = new double[knownVarNum];
        }// end of gp loop
    }// end of bdryPoint loop
    
    /* Assign the number of equations in order in eqNum vectors */
    
    int k[3], lth[3] = {n,n,dimBasedValue(dim, 1, n)}, unknownVar = 0, knownVar = 0;
    
    for(k[2] = 0; k[2]<dimBasedValue(dim, 1, n); k[2]++){
        for(k[1] = 0; k[1]<n; k[1]++){
        for(k[0] = 0; k[0]<n; k[0]++){
                if((COND(k[fixedAxis[0]],side1)||COND(k[fixedAxis[1]],side2))){
                    CornerMat[corner].eqNum[ind(k,lth)] = unknownVar;
                    unknownVar++;
                }
                else{
                    CornerMat[corner].eqNum[ind(k,lth)] = unknownVarNum + knownVar;
                    knownVar++;
                }// end of if statement
            }// end of k[0] loop
        }// end of k[1] loop
    }// end of k[2] loop
    
    getCornerMatrix(CornerMat[corner].CaVec, CornerMat[corner].CbVec, CornerMat[corner].eqNum, bdryRange, bdryNg, side1, side2, varAxis, fixedAxis, NU1, NU2);
}// end of prepCornerMatrix

void Lcbc::getCornerMatrix(double ***&CaVec, double ***&CbVec, int *eqNum, int bdryRange[3][2], int bdryNg, int side1, int side2, int varAxis, int fixedAxis[2], int NU1, int NU2){
    
    int n = (2*p+1), MU = n;
    
    /* Some needed numbers */
    int totalVarNum = (n*n*dimBasedValue(dim, 1, n));
    int compCondNum1 = NU1*MU*dimBasedValue(dim, 1, MU);
    int compCondNum2 = NU2*MU*dimBasedValue(dim, 1, MU);
    int auxiliaryEqNum = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
    int interiorEqNum = (p+1)*(p+1)*dimBasedValue(dim, 1,n);
    int approxEqNum = (compCondNum1 + compCondNum2 + auxiliaryEqNum);
    int totalEqNum  = approxEqNum + interiorEqNum;
    int unknownVarNum = totalVarNum - interiorEqNum;
    
    /* define the LCBC matrix and its blocks */
    double *A = new double[(totalEqNum*totalVarNum)];
    double *A_scale = new double[approxEqNum];
    double *A11 = new double[(approxEqNum*unknownVarNum)];
    double *A12 = new double[(approxEqNum*interiorEqNum)];
    
    /* Prepare the D matrix */
    double *D = new double[(approxEqNum*approxEqNum)];
    set1DArrayToZero(D, (approxEqNum*approxEqNum));
    getD_corner(D, varAxis, fixedAxis, NU1, NU2);
    
    /* prepare some values for the interior */
    int interiorRange[3][2] = {{0,(2*p)},{0,(2*p)},{0,dimBasedValue(dim,0,(2*p))}};
    interiorRange[fixedAxis[0]][0] = sideBasedValue(side1, p, 0);
    interiorRange[fixedAxis[0]][1] = sideBasedValue(side1, (2*p), p);
    interiorRange[fixedAxis[1]][0] = sideBasedValue(side2, p, 0);
    interiorRange[fixedAxis[1]][1] = sideBasedValue(side2, (2*p), p);
    
    int eqNumLth[] = {n,n,dimBasedValue(dim,1,n)}, intInd[3];
    int eqNumInd[] = {p,p,dimBasedValue(dim, 0, p)};
    
    int Ind[3];
    for(Ind[2] = bdryRange[2][0]; Ind[2]<=bdryRange[2][1]; Ind[2]++){
        for(Ind[1] = bdryRange[1][0]; Ind[1]<=bdryRange[1][1]; Ind[1]++){
            for(Ind[0] = bdryRange[0][0]; Ind[0]<=bdryRange[0][1]; Ind[0]++){
                
                set1DArrayToZero(A,(totalEqNum*totalVarNum)); 
                
                /* Assign the matrix part corresponding to approximated equations */
                fillCornerMatrix_LagrangeDeriv(A, Ind, side1, side2, varAxis, fixedAxis, eqNum, NU1, NU2, compCondNum1, compCondNum2, auxiliaryEqNum,totalEqNum);
                
                int row = approxEqNum;
                
                /* Assign the matrix part corresponding to interior values */
                for(intInd[2] = interiorRange[2][0]; intInd[2]<=interiorRange[2][1]; intInd[2]++){
                    for(intInd[1] = interiorRange[1][0]; intInd[1]<=interiorRange[1][1]; intInd[1]++){
                        for(intInd[0] = interiorRange[0][0]; intInd[0]<=interiorRange[0][1]; intInd[0]++){
                            
                            int colm = eqNum[ind(intInd,eqNumLth)];
                            
                            A[ind2(row,colm,totalEqNum,totalVarNum)] = 1.0;
                            row++;
                        }// intInd[0]
                    }// intInd[1]
                }// intInd[2]
                
                /* scale the matrix A */
                scaleRows(A, A_scale, approxEqNum, totalEqNum, totalVarNum);
                
                /* extract the matrix A11 */
                int A11_extractRows[]  = {0,approxEqNum};
                int A11_extractColms[] = {0,unknownVarNum};
                extractBlocks(A, A11, totalEqNum, totalVarNum, A11_extractRows, A11_extractColms, 1);
                
                /* extract the matrix A12 */
                int A12_extractRows[]  = {0,approxEqNum};
                int A12_extractColms[] = {unknownVarNum,totalVarNum};
                extractBlocks(A, A12, totalEqNum, totalVarNum, A12_extractRows, A12_extractColms, (-1));
                
                /* Scale the matrix D */
                double *D_scaled = new double[(approxEqNum*approxEqNum)];
                scaleD(D_scaled, D, A_scale, approxEqNum);
                
                /* Work out Ca and Cb via LS */
                int cVecInd = (cstCoef)?(0):(Ind[varAxis]); // index of CaVec and CbVec
                double *Ca = new double[(unknownVarNum*approxEqNum)];
                double *Cb = new double[(unknownVarNum*interiorEqNum)];
                LSsolve(A11, D_scaled, Ca, approxEqNum, unknownVarNum, approxEqNum);
                LSsolve(A11,A12,Cb,approxEqNum,unknownVarNum,interiorEqNum);
                
                delete [] D_scaled;
                
                /* Extract the row vectors out of Ca and Cb needed in ghost evaluation */
                int vecNum = 0;
                
                for(int g2 = sideBasedValue(side2, 0, 1); g2<=sideBasedValue(side2, (2*p-1), (2*p)); g2++){
                    for(int g1 = sideBasedValue(side1, 0, 1); g1<=sideBasedValue(side1, (2*p-1), (2*p)); g1++){
                        if((COND(g1, side1)||COND(g2, side2))){
                            eqNumInd[fixedAxis[0]] = g1;
                            eqNumInd[fixedAxis[1]] = g2;
                            int k = eqNum[ind(eqNumInd,eqNumLth)];
                            
                            /* extract the kth row from C_alpha */
                            int Ca_extractedRows[] = {k,(k+1)};
                            int Ca_extractedColms[] = {0,approxEqNum};
                            extractBlocks(Ca, CaVec[cVecInd][vecNum], unknownVarNum, approxEqNum, Ca_extractedRows, Ca_extractedColms, 1);
                            
                            /* extract the kth row from Cbeta */
                            int Cb_extractedRows[] = {k,(k+1)};
                            int Cb_extractedColms[] = {0,interiorEqNum};
                            extractBlocks(Cb, CbVec[cVecInd][vecNum], unknownVarNum, interiorEqNum, Cb_extractedRows, Cb_extractedColms, 1);
                            
                            vecNum++;
                        }// end of if statement
                    }// end of g1 loop
                }// end of g2 loop
                
                delete [] Ca;
                delete [] Cb;
                
                /* If cstCoef, break out of the loops */
                if(cstCoef){
                    for(int ax = 0; ax<3; ax++){
                        Ind[ax] = bdryRange[ax][1];
                    }// end ax loop (short for axis)
                }//end if cstCoef
                /*------------------------------------*/
                
            }// end of Ind[0] loop
        }// end of Ind[1] loop
    }// end of Ind[2] loop
    
    delete [] A;
    delete [] A_scale;
    delete [] A11;
    delete [] A12;
    delete [] D;
    
}// end of getCornerMatrix Function

void Lcbc::getD_corner(double *&D, int varAxis, int fixedAxis[2], int NU1, int NU2){
    int n = (2*p+1);
    int compCondNum1   = NU1*n*dimBasedValue(dim,1,n);
    int compCondNum2   = NU2*n*dimBasedValue(dim,1,n);
    int auxiliaryEqNum = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
    int approxEqNum    = compCondNum1 + compCondNum2 + auxiliaryEqNum;
    int dof            = (n*dimBasedValue(dim, 1, n));
    int NU[] = {NU1, NU2};
    
    /* use the delta approach to determine D */
    double *Id = new double[(dof*dof)];
    getIDmatrix(Id, dof);
    
    /* R corresponds to derivatives of data functions on stencil. b is the RHS vector */
    double b[(compCondNum1 + compCondNum2)], **R = new double*[(NU1 + NU2)];
    
    for(int nu = 0; nu<(NU1 + NU2); nu++)
        R[nu] = new double[dof];
    
    int colm = 0;
    for(int fixedAxisInd = 0; fixedAxisInd<2; fixedAxisInd++){
        int Rlth[] = {n,n,dimBasedValue(dim, 1, n)}; // this is a vector that carries the length of R in each axis
        Rlth[fixedAxis[fixedAxisInd]] = 1;
        for(int nu = 0; nu<NU[fixedAxisInd]; nu++){
            set2DArrayToZero(R, (NU1 + NU2), dof);
            
            for(int l = 0; l<dof; l++){
                int row = 0;
                int stenInd[] = {0,0,0};
                for(stenInd[varAxis] = 0; stenInd[varAxis]<dimBasedValue(dim, 1, n); stenInd[varAxis]++){
                    for(stenInd[fixedAxis[(1-fixedAxisInd)]] = 0; stenInd[fixedAxis[(1-fixedAxisInd)]]<n; stenInd[fixedAxis[(1-fixedAxisInd)]]++){
        
                        R[ind2(nu,fixedAxisInd,NU[0],2)][ind(stenInd,Rlth)] = Id[ind2(row,l,dof,dof)];
                        
                        row++;
                    }// end of i loop
                }// end of j loop
                
                getbVec_corner(b, R, varAxis, fixedAxis, NU1, NU2);
                
                for(int k = 0; k<(compCondNum1 + compCondNum2); k++){
                    D[ind2(k,colm,approxEqNum,approxEqNum)] = b[k];
                }// end of k loop
                colm++;
                
            }// end of l loop
        }// end of nu loop
    }// end of fixedAxisInd loop
    
    for(int nu = 0; nu<(NU1 + NU2); nu++){
        delete [] R[nu];
        R[nu] = NULL;
    }
    delete [] R;
    delete [] Id;
    
}// end of getD function

void Lcbc::getbVec_corner(double b[], double **R, int varAxis, int fixedAxis[2], int NU1, int NU2){
    int NU = (p+1), n = (2*p+1), order1, order2;
    int MU0 = n, MU1 = dimBasedValue(dim, 1, n);
    int compCondNum1 = NU1*MU0*MU1;
    int compCondNum2 = NU2*MU0*MU1;
    int compCondNum  = compCondNum1 + compCondNum2;
    memset(b, 0, (compCondNum*sizeof(double)));

    int Ind[] = {p,dimBasedValue(dim,0,p),0};
    int sizeRnu[] = {n,dimBasedValue(dim, 1, n),1};
    
    int var1 = MIN(fixedAxis[1], varAxis);
    int var2 = MAX(fixedAxis[1], varAxis);
    double hx0[] = {G.dx[var1], G.dx[var2],0};
    
    var1 = MIN(fixedAxis[0], varAxis);
    var2 = MAX(fixedAxis[0], varAxis);
    double hx1[] = {G.dx[var1], G.dx[var2],0};
    
    int eqnNum = 0, mu[3] = {0,0,0};
    for(int nu = 0; nu<NU; nu++){
        order1 = (nu == 0)?(2*p):(2*(NU1 - nu));
        order2 = (nu == 0)?(2*p):(2*(NU2 - nu));
        
        int ord1[] = {order1,dimBasedValue(dim,0,order1),0};
        int ord2[] = {order2,dimBasedValue(dim,0,order2),0};
        
        for(mu[1] = 0; mu[1]<MU1; mu[1]++){
            for(mu[0] = 0; mu[0]<MU0; mu[0]++){
                
                if(nu<NU1){
                    b[eqnNum] = mixedDeriv(R[ind2(nu,0,NU1,2)], Ind, mu, hx0, ord1, sizeRnu);
                    eqnNum++;
                }
                
                if(nu<NU2){
                    b[eqnNum] =mixedDeriv(R[ind2(nu,1,NU1,2)], Ind, mu, hx1, ord2, sizeRnu);
                    eqnNum++;
                }
                
                if( (mu[0]!=0) && (mu[0]%2)==0){
                    if(ord1[0]>2){ ord1[0] = ord1[0] - 2;}
                    if(ord2[0]>2){ ord2[0] = ord2[0] - 2;}
                }// end of if mu0
            }// end of mu0 loop
            if( (mu[1]!=0) && (mu[1]%2)==0){
                if(ord1[1]>2){ ord1[1] = ord1[1] - 2;}
                if(ord2[1]>2){ ord2[1] = ord2[1] - 2;}
            }// end of if mu1
        }// end of mu1 loop
    }// end of nu loop
}// end of getbVec

void Lcbc::fillCornerMatrix_LagrangeDeriv(double *&Matrix, int *Ind, int side1, int side2, int varAxis, int fixedAxis[2], int *eqNum, int NU1, int NU2, int compCondNum1, int compCondNum2, int auxiliaryEqNum, int totalEqNum){
    
    int n = (2*p+1), NU = (p+1), MU = n;
    int totalVarNum = (n*n*dimBasedValue(dim, 1, n));
    
    /* Prepare Variables that will hold the Lagrange derivative information */
    double Z1[totalVarNum][compCondNum1];
    double Z2[totalVarNum][compCondNum2];
    double Y[totalVarNum][auxiliaryEqNum];
    
    /* Fill Z and Y with appropriate values */
    int LagInd[3], LagIndLth[3] = {n,n,dimBasedValue(dim, 1, n)};
    
    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim, 1, n); LagInd[2]++){
        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                getLagrangeDeriv(Y[ind(LagInd,LagIndLth)], Z1[ind(LagInd,LagIndLth)], Ind, LagInd, fixedAxis[0], side1);
                
                getLagrangeDeriv(Y[ind(LagInd,LagIndLth)], Z2[ind(LagInd,LagIndLth)], Ind, LagInd, fixedAxis[1], side2,false);
            }// LagInd[0]
        }// LagInd[1]
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
                                Matrix[ind2(row,colm,totalEqNum,totalVarNum)] = Z1[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU1,MU,dimBasedValue(dim, 1, MU))];
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    row++;
                }
                
                if(nu<NU2){
                    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim,1,n); LagInd[2]++){
                        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                        for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                int colm = eqNum[ind(LagInd,LagIndLth)];
                                    Matrix[ind2(row,colm,totalEqNum,totalVarNum)] =     Z2[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU2,MU,dimBasedValue(dim, 1, MU))];
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    row++;
                }
            }// mu0 loop
        }// mu1 loop
    }// end of nu loop
    
    int auxEqNum = 0;
    for(int m2 = 0; m2<dimBasedValue(dim,1,(2*p+1)); m2++){
        for(int m1 = 0; m1<(2*p+1); m1++){
            for(int m0 = 0; m0<(2*p+1); m0++){
                if((m0 + m1 + m2) > (2*p)){
                    
                    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim,1,n); LagInd[2]++){
                        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                        for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                int colm = eqNum[ind(LagInd,LagIndLth)];
                                Matrix[ind2(row,colm,totalEqNum,totalVarNum)] = Y[ind(LagInd,LagIndLth)][auxEqNum];
                                
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    row++;
                    auxEqNum++;
                }// end of if statement
            }// end of m0 loop
        }// end of m1 loop
    }// end of m2 loop
}// end of fillMatrix_LagrangeDeriv


void Lcbc::prepDataVec_corner(double **R, double **&Rv, double t, double dt, int side1, int side2, int varAxis, int approxEqNum, int NU1, int NU2){
    
    /* prepare needed parameters */
    int NU[] = {NU1, NU2};
    int fixedAxis[2]; getFixedAxis(fixedAxis, varAxis);
    int fixedSide[] = {side1, side2};
    int fixedFace[] = {ind2(side1,fixedAxis[0],2,3), ind2(side2,fixedAxis[1],2,3)};
    
    /* prepare the boundary range at a corner or edge */
    int bdryRange[3][2]; getBdryRange_corner(bdryRange, fixedAxis, fixedSide, varAxis, (+p), (-p));
    int bdryPoints = dimBasedValue(dim, 1, G.Ngx[varAxis]);
    
    /* prepare values needed to retrieve the data derivative information */
    int faceBdryPts[2][3];
    for(int face = 0; face<2; face++){
        memcpy(faceBdryPts[face], G.Ngx, sizeof(faceBdryPts[face]));
        faceBdryPts[face][fixedAxis[face]] = 1;
    }// end of face loop
    
    /* Assign values of R to Rv */
    int Ind[3];
    for(Ind[2] = bdryRange[2][0]; Ind[2]<=bdryRange[2][1]; Ind[2]++){
        for(Ind[1] = bdryRange[1][0]; Ind[1]<=bdryRange[1][1]; Ind[1]++){
            for(Ind[0] = bdryRange[0][0]; Ind[0]<=bdryRange[0][1]; Ind[0]++){
                
                int r = 0;
                
                for(int fixedAxisInd = 0; fixedAxisInd<2; fixedAxisInd++){
                    for(int nu = 0; nu<NU[fixedAxisInd]; nu++){
                        int stenInd[3] = {0,0,0};
                        for(stenInd[varAxis] = dimBasedValue(dim, 0, (-p)); stenInd[varAxis]<= dimBasedValue(dim, 0, p); stenInd[varAxis]++){
                            
                            for(stenInd[fixedAxis[(1-fixedAxisInd)]] = (-p); stenInd[fixedAxis[(1-fixedAxisInd)]]<= p; stenInd[fixedAxis[(1-fixedAxisInd)]]++){
  
                                int Rind[] = sumVectors(Ind, stenInd);
                                Rind[fixedAxis[fixedAxisInd]] = 0;
                                Rv[Ind[varAxis]][r] =  R[ind2(nu,fixedFace[fixedAxisInd],(p+1),6)][ind(Rind,faceBdryPts[fixedAxisInd])];
                                r++;
                            }// end of stenInd loop
                        }// end of stenInd loop
                        
                    }// end of nu loop
                }// end of fixedAxisInd loop
                
                for(int row = r; row<approxEqNum; row++){
                    Rv[Ind[varAxis]][row] = 0;
                }
                
            }// end of Ind[0] loop
        }// end of Ind[1] loop
    }// end of Ind[2] loop
    
}// end of prepDataVec_corner Function

void Lcbc::getCornerGhost(double *&un, double **Rv, double ***CaVec, double ***CbVec, int *eqNum, int side1, int side2, int varAxis, int approxEqNum){
    
    /* define necessary variables */
    int interiorEqNum = (p+1)*(p+1)*dimBasedValue(dim, 1, (2*p+1));
    int fixedAxis[2]; getFixedAxis(fixedAxis, varAxis);
    int fixedSide[2] = {side1, side2};
    int bdryRange[3][2]; getBdryRange_corner(bdryRange, fixedAxis, fixedSide, varAxis, p, (-p));
    
    double b[interiorEqNum];
    
    /* Evaluate the solution at the ghost points */
    int Ind[3];
    
    for(Ind[2] = bdryRange[2][0]; Ind[2]<=bdryRange[2][1]; Ind[2]++){
        for(Ind[1] = bdryRange[1][0]; Ind[1]<=bdryRange[1][1]; Ind[1]++){
            for(Ind[0] = bdryRange[0][0]; Ind[0]<=bdryRange[0][1]; Ind[0]++){
                
                /* get the b interior here */
                getbInt_corner(b, un, eqNum, Ind, side1, side2, fixedAxis);
                
                int vecNum = 0;
                
                int cVecInd = (cstCoef)?(0):(Ind[varAxis]);
                
                for(int g2 = sideBasedValue(side2, (-p), (-p+1)); g2<=sideBasedValue(side2, (p-1), p); g2++){
                    for(int g1 = sideBasedValue(side1, (-p), (-p+1)); g1<=sideBasedValue(side1, (p-1), p); g1++){
                        if((COND((g1+p), side1)||COND((g2+p), side2))){
                            
                            double Ca_vecRv = dotProduct(CaVec[cVecInd][vecNum], Rv[Ind[varAxis]], approxEqNum);
                            double Cb_vecb  = dotProduct(CbVec[cVecInd][vecNum], b, interiorEqNum);
                            
                            int uInd[3]; uInd[varAxis] = Ind[varAxis];
                            uInd[fixedAxis[0]] = Ind[fixedAxis[0]] + g1;
                            uInd[fixedAxis[1]] = Ind[fixedAxis[1]] + g2;
                            
                            un[solInd(uInd,G.Ngx)] = Ca_vecRv + Cb_vecb;
                            vecNum++;
                        }// end of if statement
                    }// end of g1 loop
                }// end of g2 loop
                
            }// end of Ind[0] loop
        }// end of Ind[1] loop
    }// end of Ind[2] loop
    
}// end of getCornerGhost

void Lcbc::getbInt_corner(double b[], double *un, int *eqNum, int *Ind, int side1, int side2, int fixedAxis[2]){
    int n = (2*p+1);
    int interiorEqNum = (p+1)*(p+1)*dimBasedValue(dim, 1, n);
    int totalVarNum   = n*n*dimBasedValue(dim, 1, n);
    int unknownVarNum = totalVarNum - interiorEqNum;
    
    int interiorRange[3][2] = {{0,(2*p)},{0,(2*p)},{0,dimBasedValue(dim,0,(2*p))}};
    interiorRange[fixedAxis[0]][0] = sideBasedValue(side1, p, 0);
    interiorRange[fixedAxis[0]][1] = sideBasedValue(side1, (2*p), p);
    interiorRange[fixedAxis[1]][0] = sideBasedValue(side2, p, 0);
    interiorRange[fixedAxis[1]][1] = sideBasedValue(side2, (2*p), p);
    
    int eqNumLth[] = {n,n,dimBasedValue(dim,1,n)}, intInd[3];
    
    for(intInd[2] = interiorRange[2][0]; intInd[2]<=interiorRange[2][1]; intInd[2]++){
        for(intInd[1] = interiorRange[1][0]; intInd[1]<=interiorRange[1][1]; intInd[1]++){
            for(intInd[0] = interiorRange[0][0]; intInd[0]<=interiorRange[0][1]; intInd[0]++){
                int r  = eqNum[ind(intInd,eqNumLth)] - unknownVarNum;
                
                int u_ind[] = {(Ind[0] - p + intInd[0]),(Ind[1] - p + intInd[1]),dimBasedValue(dim, 0, (Ind[2] - p + intInd[2]))};
                
                b[r]  =  un[solInd(u_ind,G.Ngx)];
            }// intInd[0]
        }// intInd[1]
    }// intInd[2]
    
}// end of getbInt_corner



void Lcbc::getBdryRange_corner(int bdryRange[3][2], int fixedAxis[2], int fixedSide[2], int varAxis, int addOnSide0, int addOnSide1){
    
    for(int side = 0; side<2; side++){
        bdryRange[fixedAxis[0]][side] = G.indexRange[fixedAxis[0]][fixedSide[0]];
        bdryRange[fixedAxis[1]][side] = G.indexRange[fixedAxis[1]][fixedSide[1]];
        if(dim == 3){
            if(faceEval[(2*varAxis)]>=0){
                bdryRange[varAxis][side] = G.indexRange[varAxis][side] + (1-side)*addOnSide0 + addOnSide1*side;
            }
            else{
                bdryRange[varAxis][side] = G.indexRange[varAxis][side];
            }
        }else{
            bdryRange[2][side] = G.indexRange[2][side];
        }
    }
}// end of getBdryRange_corner

void Lcbc::getFixedAxis(int fixedAxis[2], int varAxis){
    int cnt = 0;
    for(int axis = 0; axis<3; axis++){
        if(axis != varAxis){
            fixedAxis[cnt] = axis;
            cnt++;
        }
    }// end of axis loop
}// end of getFixedAxis function


