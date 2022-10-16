#include "LCBC.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utility.h"
#include "LCBCmacros.h"

void Lcbc::updateVertexGhostPeriodic(double *&unp1, int side0, int side1, int side2){
    
    int fixedSides[3] = {side0, side1, side2};
    int perIndexRange[3][2];
    for(int axis = 0; axis<3; axis++){
        perIndexRange[axis][0] = G.indexRange[axis][fixedSides[axis]] + sideBasedValue(fixedSides[axis], (-p), 1);
        perIndexRange[axis][1] = G.indexRange[axis][fixedSides[axis]] + sideBasedValue(fixedSides[axis], (-1), p);
    }
    
    int vertices[3] = {side0, (side1 + 2), (side2 + 4)};

    int i[3];
    for(i[2] = perIndexRange[2][0]; i[2]<=perIndexRange[2][1]; i[2]++){
        for(i[1] = perIndexRange[1][0]; i[1]<=perIndexRange[1][1]; i[1]++){
            for(i[0] = perIndexRange[0][0]; i[0]<=perIndexRange[0][1]; i[0]++){
                int I[3]; memcpy(I, i, sizeof(I));
                
                for(int axis = 0; axis<3; axis ++){
                    if(faceEval[vertices[axis]] == -1){ // periodic boundary
                        I[axis] = I[axis] + sideBasedValue(fixedSides[axis], G.Nx[axis], (-G.Nx[axis]));
                    }
                }// end of axis loop

                unp1[solInd(i, G.Ngx)] = unp1[solInd(I, G.Ngx)];

            }// end of i[0]
        }// end of i[1]
    }// end of i[2]
}// end of updateVertexGhostPeriodic


void Lcbc::updateVertexGhost(double **R, double *&unp1, double t, double dt, int side0, int side1, int side2){
    
    int face[] = {(side0), (side1+2), (side2+4)};
    int NU[3], n = (2*p+1), compCondNum[3];
    for(int axis = 0; axis<3; axis++){
        NU[axis] = (faceEval[face[axis]] == 1)?(p+1):(p);
        compCondNum[axis] = NU[axis]*n*n;
    }// end of axis loop
    
    int auxiliaryEqNum = (p*(2*p+1)*(10*p+7))/3;
    int approxEqNum    = compCondNum[0] + compCondNum[1] + compCondNum[2] + auxiliaryEqNum;
    
    double *Rv = new double[approxEqNum];
    int vertex = ind3(side0,side1,side2,2,2,2);
    
    /* prepare the matrix for a vertex */
    if(VertexMat[vertex].flag == false){
        prepVertexMatrix(side0,side1,side2, NU, compCondNum, approxEqNum);
    }
    
    /* prepare the RHS vector for a vertex*/
    prepDataVec_vertex(R, Rv, t, dt, approxEqNum, side0, side1, side2, NU);
    
    /* evaluate the ghost points */
    getVertexGhost(unp1, Rv, VertexMat[vertex].CaVec[0], VertexMat[vertex].CbVec[0], VertexMat[vertex].eqNum, approxEqNum, side0, side1, side2);
    
    /* delete allocated space */
    delete [] Rv;
}// end of updateVertexGhost

void Lcbc::prepVertexMatrix(int side0, int side1, int side2, int NU[3], int compCondNum[3], int approxEqNum){
    
    int vertex = ind3(side0,side1,side2,2,2,2);
    VertexMat[vertex].flag = true;
    
    int n = (2*p+1);
    int totalVarNum   = n*n*n;
    int knownVarNum   = (p+1)*(p+1)*(p+1);
    int unknownVarNum = (totalVarNum - knownVarNum);
    int ghostPointsNum = unknownVarNum - 3*p;
    
    VertexMat[vertex].eqNum = new int[totalVarNum];
    VertexMat[vertex].CaVec = new double**[1];
    VertexMat[vertex].CbVec = new double**[1];
    
    VertexMat[vertex].CaVec[0] = new double*[ghostPointsNum];
    VertexMat[vertex].CbVec[0] = new double*[ghostPointsNum];
    
    for(int gp = 0; gp<ghostPointsNum; gp++){
        VertexMat[vertex].CaVec[0][gp] = new double[approxEqNum];
        VertexMat[vertex].CbVec[0][gp] = new double[knownVarNum];
    }// end of gp loop
    
    /* Assign the number of equations in order in eqNum vectors */
    
    int k[3], lth[3] = {n,n,n}, unknownVar = 0, knownVar = 0;

    for(k[2] = 0; k[2]<n; k[2]++){
        for(k[1] = 0; k[1]<n; k[1]++){
            for(k[0] = 0; k[0]<n; k[0]++){
                
                if((COND(k[0],side0)||COND(k[1],side1)||COND(k[2],side2))){
                    VertexMat[vertex].eqNum[ind(k,lth)] = unknownVar;
                    unknownVar++;
                }
                else{
                    VertexMat[vertex].eqNum[ind(k,lth)] = unknownVarNum + knownVar;
                    knownVar++;
                }// end of if statement
            }// end of k[0] loop
        }// end of k[1] loop
    }// end of k[2] loop
    
    /* get the vertex Matrix */
    getVertexMatrix(VertexMat[vertex].CaVec[0], VertexMat[vertex].CbVec[0], VertexMat[vertex].eqNum, side0, side1, side2, NU, compCondNum);
    /* -------------- */
}// end of prepVertexMatrix

void Lcbc::getVertexMatrix(double **&CaVec, double **&CbVec, int *eqNum, int side0, int side1, int side2, int NU[3], int compCondNum[3]){
    
    int n = (2*p+1);
    
    /* prepare needed parameters */
    int totalVarNum    = (n*n*n);
    int auxiliaryEqNum = (p*(2*p+1)*(10*p+7))/3;
    int interiorEqNum  = (p+1)*(p+1)*(p+1);
    int approxEqNum    = (compCondNum[0] + compCondNum[1] + compCondNum[2] + auxiliaryEqNum);
    int totalEqNum     = approxEqNum + interiorEqNum;
    int unknownVarNum  = totalVarNum - interiorEqNum;
    
    /* define the LCBC matrix and its blocks */
    double *A = new double[(totalEqNum*totalVarNum)];
    double *A_scale = new double[approxEqNum];
    double *A11 = new double[(approxEqNum*unknownVarNum)];
    double *A12 = new double[(approxEqNum*interiorEqNum)];
    
    /* prepare the D matrix here */
    double *D = new double[(approxEqNum*approxEqNum)];
    getD_vertex(D, NU, compCondNum, approxEqNum);
        
    /* prepare the interior indices 2D array */
    int interiorRange[3][2], side[] = {side0,side1,side2};
    for(int axis = 0; axis<3; axis++){
        interiorRange[axis][0] = sideBasedValue(side[axis], p, 0);
        interiorRange[axis][1] = sideBasedValue(side[axis], (2*p), p);
    }

    int Ind[3] = {G.indexRange[0][side0], G.indexRange[1][side1], G.indexRange[2][side2]};
    
    /* fill the Vertex Matrix by the derivatives of Lagrange polynomials */
    set1DArrayToZero(A,(totalEqNum*totalVarNum));
    set1DArrayToZero(A_scale, approxEqNum);
    fillVertexMatrix_LagrangeDeriv(A, Ind, side0, side1, side2, eqNum, NU, compCondNum, auxiliaryEqNum, totalEqNum);
    
    int row = approxEqNum;
    
    /* Assign the matrix part corresponding to interior values */
    int intInd[3], eqNumLth[] = {n,n,n};
    
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
    double *Ca = new double[(unknownVarNum*approxEqNum)];
    double *Cb = new double[(unknownVarNum*interiorEqNum)];
    LSsolve(A11, D_scaled, Ca, approxEqNum, unknownVarNum, approxEqNum);
    LSsolve(A11,A12,Cb,approxEqNum,unknownVarNum,interiorEqNum);
    
    delete [] D_scaled;
    
    int g[3], vecNum = 0;
    for(g[2] = sideBasedValue(side2, 0, 1); g[2]<=sideBasedValue(side2, (2*p-1), (2*p)); g[2]++){
        for(g[1] = sideBasedValue(side1, 0, 1); g[1]<=sideBasedValue(side1, (2*p-1), (2*p)); g[1]++){
            for(g[0] = sideBasedValue(side0, 0, 1); g[0]<=sideBasedValue(side0, (2*p-1), (2*p)); g[0]++){
                if((COND(g[0], side0)||COND(g[1], side1)||COND(g[2], side2))){
                    int k = eqNum[ind(g,eqNumLth)];
                    
                    /* extract the kth row from C_alpha */
                    int Ca_extractedRows[] = {k,(k+1)};
                    int Ca_extractedColms[] = {0,approxEqNum};
                    extractBlocks(Ca, CaVec[vecNum], unknownVarNum, approxEqNum, Ca_extractedRows, Ca_extractedColms, 1);
                    
                    /* extract the kth row from Cbeta */
                    int Cb_extractedRows[] = {k,(k+1)};
                    int Cb_extractedColms[] = {0,interiorEqNum};
                    extractBlocks(Cb, CbVec[vecNum], unknownVarNum, interiorEqNum, Cb_extractedRows, Cb_extractedColms, 1);
                    
                    vecNum++;
                }// end of if statement
            }// end of g0 loop
        }// end of g1 loop
    }// end of g2 loop
    
    delete [] Ca;
    delete [] Cb;
    
    delete [] A;
    delete [] A_scale;
    delete [] A11;
    delete [] A12;
    delete [] D;
    
}// end of getVertexMatrix

void Lcbc::getD_vertex(double *&D, int NU[3], int compCondNum[3], int approxEqNum){

    int n = (2*p+1);
    int dof = (n*n);
    int totalCompCondNum = compCondNum[0] + compCondNum[1] + compCondNum[2];
    int totalNuNum = NU[0] + NU[1] + NU[2];
    
    set1DArrayToZero(D, (approxEqNum*approxEqNum));
    
    /* use the delta approach to determine D */
    double *Id = new double[(dof*dof)];
    getIDmatrix(Id, dof);
    
    /* R corresponds to derivatives of data functions on stencil. b is the RHS vector */
    double b[totalCompCondNum], **R = new double*[totalNuNum];
    
    for(int nu = 0; nu<(totalNuNum); nu++)
        R[nu] = new double[dof];
    
    int colm = 0;
    for(int nu = 0; nu<(totalNuNum); nu++){
        set2DArrayToZero(R, (totalNuNum), dof);
        for(int l = 0; l<dof; l++){
            int row = 0;
            
            for(int j = 0; j<n; j++){
                for(int i = 0; i<n; i++){
                    R[nu][ind2(i,j,n,n)] = Id[ind2(row,l,dof,dof)];
                    
                    row++;
                }// end of i loop
            }// end of j loop
            
            getbVec_vertex(b, R, NU, totalCompCondNum);
            
            for(int k = 0; k<(totalCompCondNum); k++){
                D[ind2(k,colm,approxEqNum,approxEqNum)] = b[k];
            }// end of k loop
            colm++;
            
        }// end of l loop
    }// end of nu loop

    /* delete variables */
    
    for(int nu = 0; nu<totalNuNum; nu++)
        delete [] R[nu];
    
    delete [] R;
    delete [] Id; 
}// end of getD_vertex

void Lcbc::getbVec_vertex(double b[], double **R, int NU[3], int totalCompCondNum){
    int n = (2*p+1), order[3], ordVec[3][3];
    int MU0 = n, MU1 = n;
    memset(b, 0, (totalCompCondNum*sizeof(double)));

    int Ind[] = {p,p,0};
    int sizeRnu[] = {n,n,1};

    double hx0[] = {G.dx[1], G.dx[2],0}; // axis 0 fixed
    double hx1[] = {G.dx[0], G.dx[2],0}; // axis 1 fixed
    double hx2[] = {G.dx[0], G.dx[1],0}; // axis 2 fixed
    
    int eqnNum = 0, mu[3] = {0,0,0};
    for(int nu = 0; nu<(p+1); nu++){
        
        for(int axis = 0; axis<3; axis++){
            order[axis] = (nu==0)?(2*p):(2*(NU[axis] - nu));
            
            ordVec[axis][0] = order[axis];
            ordVec[axis][1] = order[axis];
            ordVec[axis][2] = 0;
        }
        
        for(mu[1] = 0; mu[1]<MU1; mu[1]++){
            for(mu[0] = 0; mu[0]<MU0; mu[0]++){
                if(nu<NU[0]){
                    b[eqnNum] = mixedDeriv(R[nu], Ind, mu, hx0, ordVec[0], sizeRnu);
                    eqnNum++;
                }
                
                if(nu<NU[1]){
                    b[eqnNum] = mixedDeriv(R[(nu + NU[0])], Ind, mu, hx1, ordVec[1], sizeRnu);
                    eqnNum++;
                }
                
                if(nu<NU[2]){
                    b[eqnNum] = mixedDeriv(R[(nu + NU[0] + NU[1])], Ind, mu, hx2, ordVec[2], sizeRnu);
                    eqnNum++;
                }
                
                if( (mu[0]!=0) && (mu[0]%2)==0 ){
                    for(int axis = 0; axis<3; axis++){
                        if(ordVec[axis][0]>2){ordVec[axis][0] = ordVec[axis][0] - 2;}
                    }
                }// end of if mu0
            }// end of mu0 loop
            if( (mu[1]!=0) && (mu[1]%2)==0 ){
                for(int axis = 0; axis<3; axis++){
                    if(ordVec[axis][1]>2){ordVec[axis][1] = ordVec[axis][1] - 2;}
                }
            }// end of if mu1
        }// end of mu1 loop
    }// end of nu loop
}// end of getbVec_vertex

void Lcbc::fillVertexMatrix_LagrangeDeriv(double *&Matrix, int *Ind, int side0, int side1, int side2, int *eqNum, int NU[3], int compCondNum[3], int auxiliaryEqNum, int totalEqNum){
    
    int n = (2*p+1), MU = n;
    int totalVarNum = (n*n*n);
    int NUmax = (p+1);
    
    /* Prepare Variables that will hold the Lagrange derivative information */
    double Z0[totalVarNum][compCondNum[0]];
    double Z1[totalVarNum][compCondNum[1]];
    double Z2[totalVarNum][compCondNum[2]];
    double Y[totalVarNum][auxiliaryEqNum];
    
    /* Fill Z and Y with appropriate values */
    int LagInd[3], LagIndLth[3] = {n,n,n};
    
    for(LagInd[2] = 0; LagInd[2]<n; LagInd[2]++){
        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                getLagrangeDeriv(Y[ind(LagInd,LagIndLth)], Z0[ind(LagInd,LagIndLth)], Ind, LagInd, 0, side0);
                
                getLagrangeDeriv(Y[ind(LagInd,LagIndLth)], Z1[ind(LagInd,LagIndLth)], Ind, LagInd, 1, side1, false);
                
                getLagrangeDeriv(Y[ind(LagInd,LagIndLth)], Z2[ind(LagInd,LagIndLth)], Ind, LagInd, 2, side2, false);
            }// LagInd[0]
        }// LagInd[1]
    }// LagInd[2]
    
    int row = 0;
    for(int nu = 0; nu<NUmax; nu++){
        for(int mu1 = 0; mu1<MU; mu1++){
            for(int mu0 = 0; mu0<MU; mu0++){
                
                if(nu<NU[0]){
                    for(LagInd[2] = 0; LagInd[2]<n; LagInd[2]++){
                        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                int colm = eqNum[ind(LagInd,LagIndLth)];
                                    Matrix[ind2(row,colm,totalEqNum,totalVarNum)] = Z0[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU[0],MU,MU)];
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    row++;
                }// end of if nu<NU[0]
                
                if(nu<NU[1]){
                    for(LagInd[2] = 0; LagInd[2]<n; LagInd[2]++){
                        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                int colm = eqNum[ind(LagInd,LagIndLth)];
                                Matrix[ind2(row,colm,totalEqNum,totalVarNum)] = Z1[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU[1],MU,MU)];
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    row++;
                }// end of if nu<NU[1]
                
                if(nu<NU[2]){
                    for(LagInd[2] = 0; LagInd[2]<n; LagInd[2]++){
                        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                int colm = eqNum[ind(LagInd,LagIndLth)];
                                Matrix[ind2(row,colm,totalEqNum,totalVarNum)] = Z2[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU[2],MU,MU)];
                            }// LagInd[0]
                        }// LagInd[1]
                    }// LagInd[2]
                    row++;
                }// end of if nu<NU[2]
            }// mu0 loop
        }// mu1 loop
    }// end of nu loop
    
    int auxEqNum = 0;
    for(int m2 = 0; m2<n; m2++){
        for(int m1 = 0; m1<n; m1++){
            for(int m0 = 0; m0<n; m0++){
                if((m0 + m1 + m2) > (2*p)){
                    
                    for(LagInd[2] = 0; LagInd[2]<n; LagInd[2]++){
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
}// end of fillVertexMatrix_LagrangeDeriv

void Lcbc::prepDataVec_vertex(double **R, double *&Rv, double t, double dt, int approxEqNum, int side0, int side1, int side2, int NU[3]){
    
    /* prepare needed parameters */
    int face[] = {ind2(side0,0,2,3), ind2(side1,1,2,3), ind2(side2,2,2,3)};
    
    int Ind[3] = {G.indexRange[0][side0], G.indexRange[1][side1], G.indexRange[2][side2]};
    
    int r = 0;
    
    for(int axis = 0; axis<3; axis++){
        
        int varAxis[2]; getVarAxis(varAxis, axis);
        int stenInd[3] = {0,0,0};
        int facePoints[3] = {G.Ngx[0],G.Ngx[1],G.Ngx[2]};
        facePoints[axis] = 1;
        
        for(int nu = 0; nu<NU[axis]; nu++){
            for(stenInd[varAxis[1]] = (-p); stenInd[varAxis[1]]<=p; stenInd[varAxis[1]]++){
                for(stenInd[varAxis[0]] = (-p); stenInd[varAxis[0]]<=p; stenInd[varAxis[0]]++){
                    
                    int Rind[] = sumVectors(Ind, stenInd);
                    Rind[axis] = 0;
                    Rv[r] =  R[ind2(nu,face[axis],(p+1),6)][ind(Rind, facePoints)];
                    r++;
                }// end of stenInd0 loop
            }// end of stenInd1 loop
            
        }// end of nu loop
    }// end of axis loop
    
    for(int row = r; row<approxEqNum; row++){
        Rv[row] = 0;
    }
}// end of prepDataVec_vertex

void Lcbc::getVertexGhost(double *&un, double *Rv, double **CaVec, double **CbVec, int *eqNum, int approxEqNum, int side0, int side1, int side2){
    
    /* define parameters */
    int interiorEqNum = (p+1)*(p+1)*(p+1);
    
    double b[interiorEqNum];
    
    int Ind[3] = {G.indexRange[0][side0], G.indexRange[1][side1], G.indexRange[2][side2]};
    
    /* get the b interior here */
    getbInt_vertex(b, un, eqNum, Ind, side0, side1, side2);

    int vecNum = 0, g[3];
    for(g[2] = sideBasedValue(side2, (-p), (-p+1)); g[2]<=sideBasedValue(side2, (p-1), p); g[2]++){
        for(g[1] = sideBasedValue(side1, (-p), (-p+1)); g[1]<=sideBasedValue(side1, (p-1), p); g[1]++){
            for(g[0] = sideBasedValue(side0, (-p), (-p+1)); g[0]<=sideBasedValue(side0, (p-1), p); g[0]++){
                if((COND((g[0]+p), side0)||COND((g[1]+p), side1)||COND((g[2]+p), side2))){
                    double Ca_vecRv = dotProduct(CaVec[vecNum], Rv, approxEqNum);
                    double Cb_vecb  = dotProduct(CbVec[vecNum], b, interiorEqNum);
                    
                    int uInd[3] = sumVectors(Ind, g);
                    un[solInd(uInd,G.Ngx)] = Ca_vecRv + Cb_vecb;
                    vecNum++;
                } // end of if conditions
            }// end of g0 loop
        }// end of g1 loop
    }// end of g2 loop
    
}// end of getVertexGhost

void Lcbc::getbInt_vertex(double b[], double *un, int *eqNum, int *Ind, int side0, int side1, int side2){
    
    int n = (2*p+1);
    int interiorEqNum = (p+1)*(p+1)*(p+1);
    int totalVarNum = n*n*n;
    int unknownVarNum = totalVarNum - interiorEqNum;
    
    int interiorRange[3][2], side[] = {side0,side1,side2};
    for(int axis = 0; axis<3; axis++){
        interiorRange[axis][0] = sideBasedValue(side[axis], p, 0);
        interiorRange[axis][1] = sideBasedValue(side[axis], (2*p), p);
    }
    
    int eqNumLth[] = {n,n,n}, intInd[3];
    
    for(intInd[2] = interiorRange[2][0]; intInd[2]<=interiorRange[2][1]; intInd[2]++){
        for(intInd[1] = interiorRange[1][0]; intInd[1]<=interiorRange[1][1]; intInd[1]++){
            for(intInd[0] = interiorRange[0][0]; intInd[0]<=interiorRange[0][1]; intInd[0]++){
                int r  = eqNum[ind(intInd,eqNumLth)] - unknownVarNum;
                
                int u_ind[] = {(Ind[0] - p + intInd[0]),(Ind[1] - p + intInd[1]), (Ind[2] - p + intInd[2])};
                
                b[r]  =  un[solInd(u_ind,G.Ngx)];
            }// intInd[0]
        }// intInd[1]
    }// intInd[2]

}// end of getbInt_vertex
