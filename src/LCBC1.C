#include "LCBC.h"
#include <string.h>
#include <math.h>
#include "utility.h"
#include "LCBCmacros.h"

void getWidth(int wth[3], int lth[3], int fixedAxis, int nu, int k, int p, int dim);

void Lcbc::updateFaceGhost(double *&unp1, double **R, double t, double dt, int approxEqNum, int axis, int side){

    int bdryPointsNum = getBdryPointsNum(axis);
    double **Rv = new double*[bdryPointsNum];
    
    for(int bdryPt = 0; bdryPt<bdryPointsNum; bdryPt++){
        Rv[bdryPt] = new double[approxEqNum];
    }
    
    int face = ind2(side, axis, 2, 3);
    
    /* prepare the LCBC matrix*/
    if(FaceMat[face].flag == false){
        prepSideMatrix(axis, side);
    }
    
    /* prepare the RHS data vector */
    prepDataVec(R, Rv, t, dt, axis, side);

    /* Evaluate the solution at the ghost points */
    getFaceGhost(unp1, Rv, FaceMat[face].CaVec, FaceMat[face].CbVec, FaceMat[face].eqNum, axis, side);
    
    /* Free any allocated variables */
    for(int bdryPoint = 0; bdryPoint<bdryPointsNum; bdryPoint++)
        delete [] Rv[bdryPoint];
    
    delete [] Rv;
}// end of update FaceGhost

void Lcbc::prepSideMatrix(int axis, int side){
    
    int face = ind2(side,axis,2,3);
    FaceMat[face].flag = true;

    int bdryRange[3][2], bdryNgx[3]; getBdryRange(bdryRange, bdryNgx, axis, side, p, (-p));
    int bdryNg = bdryNgx[0]*bdryNgx[1]*bdryNgx[2];
    
    int n = (2*p+1);
    int totalVarNum = n*n*dimBasedValue(dim, 1, n);
    int knownVarNum   = (p+1)*n*dimBasedValue(dim, 1, n);
    int unknownVarNum = (totalVarNum - knownVarNum);
    int compCondNum = (p+1)*n*dimBasedValue(dim, 1, n);
    int auxiliaryEqNum = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
    int approxEqNum = compCondNum + auxiliaryEqNum; // number of equations treated via least squares
    
    FaceMat[face].eqNum = new int[totalVarNum];
    FaceMat[face].CaVec = new double**[bdryNg];
    FaceMat[face].CbVec = new double**[bdryNg];
    
    for(int bdryPoint = 0; bdryPoint<bdryNg; bdryPoint++){
        FaceMat[face].CaVec[bdryPoint] = new double*[p];
        FaceMat[face].CbVec[bdryPoint] = new double*[p];
        
        for(int k = 0; k<p; k++){
            FaceMat[face].CaVec[bdryPoint][k] = new double[approxEqNum];
            FaceMat[face].CbVec[bdryPoint][k] = new double[knownVarNum];
        }
    }
    
    /* Assign the number of equations in order in eqNum vectors */
    int k[3], lth[3] = {n,n,dimBasedValue(dim, 1, n)}, unknownVar = 0, knownVar = 0;

    for(k[0] = 0; k[0]<n; k[0]++){
        for(k[1] = 0; k[1]<n; k[1]++){
            for(k[2] = 0; k[2]<dimBasedValue(dim, 1, n); k[2]++){
                if(COND(k[axis],side)){
                    FaceMat[face].eqNum[ind(k,lth)] = unknownVar;
                    unknownVar++;
                }
                else{
                    FaceMat[face].eqNum[ind(k,lth)] = unknownVarNum + knownVar;
                    knownVar++;
                }// end of if statement
            }// end of k[0] loop
        }// end of k[1] loop
    }// end of k[2] loop

    getSideMatrix(FaceMat[face].CaVec,FaceMat[face].CbVec,FaceMat[face].eqNum,bdryRange,bdryNgx,axis,side);
}// end of prepSideMatrix

void Lcbc::getSideMatrix(double ***&CaVec, double ***&CbVec, int *eqNum, int bdryRange[3][2], int bdryNgx[3], int axis, int side){

    int n = (2*p+1), NU = (p+1), MU = n;
    /* prepare some needed numbers */
    int totalVarNum = (n*n*dimBasedValue(dim, 1, n));
    int compCondNum = NU*MU*dimBasedValue(dim, 1, MU);
    int auxiliaryEqNum = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
    int interiorEqNum = (p+1)*n*dimBasedValue(dim, 1,n);
    int approxEqNum = compCondNum + auxiliaryEqNum; // number of equations treated via least squares
    int totalEqNum  = approxEqNum + interiorEqNum;
    int unknownVarNum = totalVarNum - interiorEqNum;
    
    /* define needed variables */
    double A_scale[approxEqNum];
    double *A   = new double[(totalEqNum*totalVarNum)];
    double *A11 = new double[(approxEqNum*unknownVarNum)];
    double *A12 = new double[(approxEqNum*interiorEqNum)];
    
    /* Prepare the D matrix */
    double *D = new double[(approxEqNum*approxEqNum)];
    set1DArrayToZero(D, (approxEqNum*approxEqNum));
    getD(D, axis);
    
    int interiorRange[3][2] = {{0,(2*p)},{0,(2*p)},{0,dimBasedValue(dim,0,(2*p))}};
    interiorRange[axis][0] = sideBasedValue(side, p, 0); // if side = 0, choose p else choose 0
    interiorRange[axis][1] = sideBasedValue(side, (2*p), p);
    int eqNumLth[] = {n,n,dimBasedValue(dim,1,n)}, intInd[3];
    int eqNumInd[] = {p,p,dimBasedValue(dim, 0, p)};
 
    int Ind[3];
    for(Ind[2] = bdryRange[2][0]; Ind[2]<=bdryRange[2][1]; Ind[2]++){
        for(Ind[1] = bdryRange[1][0]; Ind[1]<=bdryRange[1][1]; Ind[1]++){
            for(Ind[0] = bdryRange[0][0]; Ind[0]<=bdryRange[0][1]; Ind[0]++){

                fillMatrix_LagrangeDeriv(A, totalEqNum, compCondNum, auxiliaryEqNum, Ind, axis, side, eqNum);

                int row = approxEqNum;

                for(intInd[0] = interiorRange[0][0]; intInd[0]<=interiorRange[0][1]; intInd[0]++){
                    for(intInd[1] = interiorRange[1][0]; intInd[1]<=interiorRange[1][1]; intInd[1]++){
                        for(intInd[2] = interiorRange[2][0]; intInd[2]<=interiorRange[2][1]; intInd[2]++){
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
            
                /* Scale the D matrix */
                double *D_scaled = new double[(approxEqNum*approxEqNum)];
                scaleD(D_scaled, D, A_scale, approxEqNum);

                /* do a LS solve to obtain Ca and Cb */
                int copy = Ind[axis]; Ind[axis] = 0;
                double *Ca = new double[(unknownVarNum*approxEqNum)];
                double *Cb = new double[(unknownVarNum*interiorEqNum)];
                LSsolve(A11, D_scaled, Ca, approxEqNum, unknownVarNum, approxEqNum);
                LSsolve(A11,A12,Cb,approxEqNum,unknownVarNum,interiorEqNum);

                delete [] D_scaled;

                int vecNum = 0;
                for(int ghostInd = sideBasedValue(side, 0, (1+p)); ghostInd<=sideBasedValue(side, (p-1), (2*p)); ghostInd++){
                    eqNumInd[axis] = ghostInd;
                    int k = eqNum[ind(eqNumInd,eqNumLth)];

                    /* extract the kth row from C_alpha */
                    int Ca_extractedRows[] = {k,(k+1)};
                    int Ca_extractedColms[] = {0,approxEqNum};
                    extractBlocks(Ca, CaVec[ind(Ind,bdryNgx)][vecNum], unknownVarNum, approxEqNum, Ca_extractedRows, Ca_extractedColms, 1);

                    /* extract the kth row from Cbeta */
                    int Cb_extractedRows[] = {k,(k+1)};
                    int Cb_extractedColms[] = {0,interiorEqNum};
                    extractBlocks(Cb, CbVec[ind(Ind,bdryNgx)][vecNum], unknownVarNum, interiorEqNum, Cb_extractedRows, Cb_extractedColms, 1);
                    vecNum++;
                }// end of ghostVal
                
                delete [] Ca;
                delete [] Cb; 

                Ind[axis] = copy;
            }// Ind[0] loop
        }// Ind[1] loop
    }// Ind[2] loop

    delete [] A;
    delete [] A11;
    delete [] A12;
    delete [] D;
    
}// end of getSideMatrix

void Lcbc::getD(double *&D, int axis){

    int n = (2*p+1);
    int compCondNum    = (p+1)*n*dimBasedValue(dim,1,n);
    int auxiliaryEqNum = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
    int approxEqNum    = compCondNum + auxiliaryEqNum;
    int dof            = (n*dimBasedValue(dim, 1, n));
    int NU             = (p+1);

    /* use the delta approach to determine D */
    double *Id = new double[(dof*dof)];
    getIDmatrix(Id, dof);

    double *b = new double[compCondNum], **R = new double*[NU];

    for(int nu = 0; nu<NU; nu++)
        R[nu] = new double[dof];

    int colm = 0;
    for(int nu = 0; nu<NU; nu++){
        set2DArrayToZero(R, NU, dof);
        for(int l = 0; l<dof; l++){
            int row = 0;
            for(int j = 0; j<dimBasedValue(dim, 1, n); j++){
                for(int i = 0; i<n; i++){
                    R[nu][ind2(i,j,n,n)] = Id[ind2(row,l,dof,dof)];
                    row++;
                }// end of i loop
            }// end of j loop
            getbVec(b, R, axis);
            for(int k = 0; k<compCondNum; k++){
                D[ind2(k,colm,approxEqNum,approxEqNum)] = b[k];
            }// end of k loop
            colm++;
        }// end of l loop
    }// end of nu loop
    
    for(int nu = 0; nu<NU; nu++){
        delete [] R[nu];
        R[nu] = NULL;
    }// end of nu loop
    delete [] R;
    delete [] Id;
    delete [] b;
}// end of getD function

void Lcbc::getbVec(double *&b, double **R, int axis){
    int NU = (p+1), n = (2*p+1), order;
    int MU0 = n, MU1 = dimBasedValue(dim, 1, n);
    int compCondNum = NU*MU0*MU1;
    memset(b, 0, (compCondNum*sizeof(double)));

    int varAxis[2]; getVarAxis(varAxis, axis);
    int Ind[] = {p,dimBasedValue(dim,0,p),0};
    int sizeRnu[] = {n,dimBasedValue(dim, 1, n),1};
    double hx[] = {G.dx[varAxis[0]], G.dx[varAxis[1]],0};

    int eqnNum = 0, mu[3] = {0,0,0};
    for(int nu = 0; nu<NU; nu++){
        (nu == 0) ? (order = 2*p) : (order = (-2*nu + 2*p + 2));

        int ord[] = {order,dimBasedValue(dim,0,order),0};

        for(mu[1] = 0; mu[1]<MU1; mu[1]++){
            for(mu[0] = 0; mu[0]<MU0; mu[0]++){

                b[eqnNum] = der.mixedNumDeriv(R[nu], 2, Ind, mu, ord, hx, sizeRnu);
                eqnNum++;

                if( (mu[0]!=0) && (mu[0]%2)==0 && ord[0]>2 ){
                    ord[0] = ord[0] - 2;
                }// end of if mu0
            }// end of mu0 loop
            if( (mu[1]!=0) && (mu[1]%2)==0 && ord[1]>2 ){
                ord[1] = ord[1] - 2;
            }// end of if mu1
        }// end of mu1 loop
    }// end of nu loop
}// end of getbVec

void Lcbc::scaleD(double *&D_scaled, double *D, double *S, int D_size){
    for(int j = 0; j<D_size; j++){
        for(int i = 0; i<D_size; i++){
            D_scaled[ind2(i,j,D_size,D_size)] = D[ind2(i,j,D_size,D_size)]/S[i];
        }
    }
    return;
}

void Lcbc::fillMatrix_LagrangeDeriv(double *&Matrix, int totalEqNum, int compCondNum, int auxiliaryEqNum, int *Ind, int axis, int side, int *eqNum){
    
    int n = (2*p+1), NU = (p+1), MU = n;
    int totalVarNum = (n*n*dimBasedValue(dim, 1, n));
    
    /* Prepare Variables that will hold the Lagrange derivative information */
    double Z[totalVarNum][compCondNum];
    double Y[totalVarNum][auxiliaryEqNum];
    
    /* Fill Z and Y with appropriate values */
    int LagInd[3], LagIndLth[3] = {n,n,dimBasedValue(dim, 1, n)};
    
    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim, 1, n); LagInd[2]++){
        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                getLagrangeDeriv_Dirichlet(Y[ind(LagInd,LagIndLth)], Z[ind(LagInd,LagIndLth)], Ind, LagInd, axis, side);
            }// LagInd[0]
        }// LagInd[1]
    }// LagInd[2]
    
    int row = 0;
    for(int nu = 0; nu<NU; nu++){
        for(int mu1 = 0; mu1<dimBasedValue(dim,1,MU); mu1++){
            for(int mu0 = 0; mu0<MU; mu0++){
                
                for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim,1,n); LagInd[2]++){
                    for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                        for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                            int colm = eqNum[ind(LagInd,LagIndLth)];
                            Matrix[ind2(row,colm,totalEqNum,totalVarNum)] =     Z[ind(LagInd,LagIndLth)][ind3(nu,mu0,mu1,NU,MU,dimBasedValue(dim, 1, MU))];
                        }// LagInd[0]
                    }// LagInd[1]
                }// LagInd[2]
                row++;
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

void Lcbc::getLagrangeDeriv_Dirichlet(double Y[], double Z[], int *Ind, int *LInd, int axis, int side){
    /* This function finds and saves derivatives of the Lagrange polynomial centered at a boundary point [Ind] */
    
    int NU = (p+1), K = (p+1), MU = 2*p+1;
    int varAxis[2]; getVarAxis(varAxis, axis);
    
    double *V[(NU*K)]; // V[nu,k][Grid_index]
    double *W[(K*p*p*p)];

    int wth[3], lth[3];
    int nu = 0; int ll = LagrangeData_center, Lw = (2*ll + 1);
    
    for(int k = 1; k<=p; k++){
        getWidth(wth, lth, axis, nu, k, p, dim);
        
        V[ind2(nu,k,NU,K)] = new double[(lth[0]*lth[1]*lth[2])];
        int i[3];
        for(i[2] = -wth[2]; i[2]<=wth[2]; i[2]++){
            for(i[1] = -wth[1]; i[1]<=wth[1]; i[1]++){
                for(i[0] = -wth[0]; i[0]<=wth[0]; i[0]++){
                    double Ln = 1;
                    for(int d = 0; d<dim; d++){
                        Ln = Ln*LagrangeData[ind2(LInd[d],(ll+i[d]),(2*p+1),Lw)];
                    }
                    int V_ind[] = sumVectors(wth, i);
                    V[ind2(nu,k,NU,K)][ind(V_ind,lth)] = Ln;
                }// end of i0
            }// end of i1
        }// end of i2
    }// end of k
    
    int ord[] = {2,2,dimBasedValue(dim, 0, 2)};

    int V_wth[3], V_lth[3];
    for(int nu = 0; nu<=(p-1); nu++){
        for(int k = 1; k<=(p-nu); k++){
            if(k != 1){
                getWidth(wth, lth, axis, nu, k, p, dim);
                for(int l = 1; l<=(k-1); l++){
                    getWidth(V_wth, V_lth, axis, nu, l, p, dim);

                    int m = (k-l);
                    for(int Dx = 0; Dx<=m; Dx++){
                        for(int Dy = dimBasedValue(dim, (m-Dx), 0); Dy<=(m-Dx); Dy++){
                            int Dz = dimBasedValue(dim, 0, (m - Dx - Dy));
                            W[ind4(l,Dx,Dy,Dz,K,p,p,p)] = new double[(lth[0]*lth[1]*lth[2])];

                                    int Deriv[] = {(2*Dx),(2*Dy),(2*Dz)};
                                    int i[3];
                                        for(i[2] = (-wth[2]); i[2]<=(wth[2]); i[2]++){
                                            for(i[1] = (-wth[1]); i[1]<=(wth[1]); i[1]++){
                                                for(i[0] = (-wth[0]); i[0]<=(wth[0]); i[0]++){
                                                    int V_ind[] = sumVectors(V_wth, i);
                                                    int W_ind[] = sumVectors(wth, i);

                                                    W[ind4(l,Dx,Dy,Dz,K,p,p,p)][ind(W_ind,lth)] = der.mixedNumDeriv(V[ind2(nu,l,NU,K)], 3, V_ind, Deriv, ord, G.dx, V_lth);
                                                }// end of i0 loop
                                            }// end of i1 loop
                                        }// end of i2 loop

                        }// end of Dy
                    }// end of Dx
                }// end of l loop
            }// end of if k statement
            getWidth(wth, lth, axis, (nu+1), k, p, dim);
            V[ind2((nu+1),k,NU,K)] = new double[(lth[0]*lth[1]*lth[2])];

            applyQh(V[ind2((nu+1),k,NU,K)], V[ind2(nu,k,NU,K)], W, Ind, nu, k, axis, side);

            if(k!= 1){
                for(int l = 1; l<=(k-1); l++){
                    int m = (k-l);
                    for(int Dx = 0; Dx<=m; Dx++){
                        for(int Dy = dimBasedValue(dim, (m-Dx), 0); Dy<=(m-Dx); Dy++){
                            int Dz = dimBasedValue(dim, 0, (m - Dx - Dy));
                            delete [] W[ind4(l,Dx,Dy,Dz,K,p,p,p)];
                            W[ind4(l,Dx,Dy,Dz,K,p,p,p)]= NULL;
                        }// end of Dx[1]
                    }// end of Dx[0]
                }// end of l loop
            }// end of if k statement

        }// end of k loop
    }// end of nu loop

    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    int MU0 = MU, MU1 = dimBasedValue(dim, 1, MU), k, mu[2]; double *Vcopy;

    int Deriv[] = {0,0,0}, order[3] = {0,0,0};
    for(int nu = 0; nu<NU; nu++){
        (nu == 0) ? (k = p):(k = (p+1-nu));
        getWidth(wth, lth, axis, nu, k, p, dim);

        Vcopy = V[ind2(nu,k,NU,K)];

        int kappa[2] = {k,k};

        for(mu[1] = 0; mu[1]<MU1; mu[1]++ ){
            for(mu[0] = 0; mu[0]<MU0; mu[0]++){

                for(int l = 0; l<2; l++){
                    Deriv[varAxis[l]] = mu[l];
                    order[varAxis[l]] = 2*kappa[l];
                }

                Z[ind3(nu,mu[0],mu[1],NU,MU0,MU1)] = der.mixedNumDeriv(Vcopy, 3, wth, Deriv, order, G.dx, lth);

                if( (mu[0]!=0) && (mu[0]%2)==0 && kappa[0]>1 ){
                    kappa[0] = kappa[0] - 1;
                }// end of if mu0
            }// end of mu1 loop
            if( (mu[1]!=0) && (mu[1]%2)==0 && kappa[1]>1 ){
                kappa[1] = kappa[1] - 1;
            }// end of if mu1
        }// end of mu1
    }// end of nu loop

    /* Calculate \partial_x^{m1}\partial_y^{m2}L_iL_j evaluated at boundary point when (m1+m2)>(2p) */
    int cnt = 0, m[3], M_lth[3] ={MU,MU,dimBasedValue(dim, 1, MU)};
    getWidth(wth, lth, axis, 0, 1, p, dim);
    for(m[2] = 0; m[2]<M_lth[2]; m[2]++){
        for(m[1] = 0; m[1]<M_lth[1]; m[1]++){
            for(m[0] = 0; m[0]<M_lth[0]; m[0]++){
                if((m[0] + m[1] + m[2])>(2*p)){
                    int order[] = {2*getOrder(m[0],p),2*getOrder(m[1],p),0};
                    (dim == 3) ? (order[2] = (2*getOrder(m[2], p))) : (order[2] = 0);
                    Y[cnt] = der.mixedNumDeriv(V[ind2(0,1,NU,K)], 3, wth, m, order, G.dx, lth);
                    cnt = cnt + 1;
                }// end if statement
            }// end m[0] loop
        }// end m[1] loop
    }// end m[2] loop
    
    /* free all the pointers allocated in this function */
    for(int k = 1; k<=p; k++){
        delete [] V[ind2(0,k,NU,K)];
        V[ind2(0,k,NU,K)] = NULL;
    }
    for(int nu = 0; nu<=(p-1); nu++){
        for(int k = 1; k<=(p-nu); k++){
            delete [] V[ind2((nu+1),k,NU,K)];
            V[ind2((nu+1),k,NU,K)] = NULL;
        }// end of k loop
    }// end of nu loop
    /* end of pointer freeing */
}// end of Dirichlet_LagrangeDer3 function

void Lcbc::getLagrangeDeriv_Dirichlet(double Z[], int *Ind, int *LInd, int axis, int side){
    /* This function finds and saves derivatives of the Lagrange polynomial centered at a boundary point [Ind] */
    
    int NU = (p+1), K = (p+1), MU = 2*p+1;
    int varAxis[2]; getVarAxis(varAxis, axis);
    
    double *V[(NU*K)]; // V[nu,k][Grid_index]
    double *W[(K*p*p*p)];

    int wth[3], lth[3];
    int nu = 0; int ll = LagrangeData_center, Lw = (2*ll + 1);
    
    for(int k = 1; k<=p; k++){
        getWidth(wth, lth, axis, nu, k, p, dim);
        
        V[ind2(nu,k,NU,K)] = new double[(lth[0]*lth[1]*lth[2])];
        int i[3];
        for(i[2] = -wth[2]; i[2]<=wth[2]; i[2]++){
            for(i[1] = -wth[1]; i[1]<=wth[1]; i[1]++){
                for(i[0] = -wth[0]; i[0]<=wth[0]; i[0]++){
                    double Ln = 1;
                    for(int d = 0; d<dim; d++){
                        Ln = Ln*LagrangeData[ind2(LInd[d],(ll+i[d]),(2*p+1),Lw)];
                    }
                    int V_ind[] = sumVectors(wth, i);
                    V[ind2(nu,k,NU,K)][ind(V_ind,lth)] = Ln;
                }// end of i0
            }// end of i1
        }// end of i2
    }// end of k
    
    int ord[] = {2,2,dimBasedValue(dim, 0, 2)};
    
    int V_wth[3], V_lth[3];
    for(int nu = 0; nu<=(p-1); nu++){
        for(int k = 1; k<=(p-nu); k++){
            if(k != 1){
                getWidth(wth, lth, axis, nu, k, p, dim);
                for(int l = 1; l<=(k-1); l++){
                    getWidth(V_wth, V_lth, axis, nu, l, p, dim);

                    int m = (k-l);
                    for(int Dx = 0; Dx<=m; Dx++){
                        for(int Dy = dimBasedValue(dim, (m-Dx), 0); Dy<=(m-Dx); Dy++){
                            int Dz = dimBasedValue(dim, 0, (m - Dx - Dy));
                                    W[ind4(l,Dx,Dy,Dz,K,p,p,p)] = new double[(lth[0]*lth[1]*lth[2])];
                            
                                    int Deriv[] = {(2*Dx),(2*Dy),(2*Dz)};
                                    int i[3];
                                        for(i[2] = (-wth[2]); i[2]<=(wth[2]); i[2]++){
                                            for(i[1] = (-wth[1]); i[1]<=(wth[1]); i[1]++){
                                                for(i[0] = (-wth[0]); i[0]<=(wth[0]); i[0]++){
                                                    int V_ind[] = sumVectors(V_wth, i);
                                                    int W_ind[] = sumVectors(wth, i);

                                                    W[ind4(l,Dx,Dy,Dz,K,p,p,p)][ind(W_ind,lth)] = der.mixedNumDeriv(V[ind2(nu,l,NU,K)], 3, V_ind, Deriv, ord, G.dx, V_lth);
                                                }// end of i0 loop
                                            }// end of i1 loop
                                        }// end of i2 loop

                        }// end of Dy
                    }// end of Dx
                }// end of l loop
            }// end of if k statement
            getWidth(wth, lth, axis, (nu+1), k, p, dim);
            V[ind2((nu+1),k,NU,K)] = new double[(lth[0]*lth[1]*lth[2])];
            
            applyQh(V[ind2((nu+1),k,NU,K)], V[ind2(nu,k,NU,K)], W, Ind, nu, k, axis, side);

            if(k!= 1){
                for(int l = 1; l<=(k-1); l++){
                    int m = (k-l);
                    for(int Dx = 0; Dx<=m; Dx++){
                        for(int Dy = dimBasedValue(dim, (m-Dx), 0); Dy<=(m-Dx); Dy++){
                            int Dz = dimBasedValue(dim, 0, (m - Dx - Dy));
                            delete [] W[ind4(l,Dx,Dy,Dz,K,p,p,p)];
                            W[ind4(l,Dx,Dy,Dz,K,p,p,p)]= NULL;
                        }// end of Dx[1]
                    }// end of Dx[0]
                }// end of l loop
            }// end of if k statement
  
        }// end of k loop
    }// end of nu loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    int MU0 = MU, MU1 = dimBasedValue(dim, 1, MU), k, mu[2]; double *Vcopy;
    
    int Deriv[] = {0,0,0}, order[3] = {0,0,0};
    for(int nu = 0; nu<NU; nu++){
        (nu == 0) ? (k = p):(k = (p+1-nu));
        getWidth(wth, lth, axis, nu, k, p, dim);
        
        Vcopy = V[ind2(nu,k,NU,K)];
        
        int kappa[2] = {k,k};
        
        for(mu[1] = 0; mu[1]<MU1; mu[1]++ ){
            for(mu[0] = 0; mu[0]<MU0; mu[0]++){
                
                for(int l = 0; l<2; l++){
                    Deriv[varAxis[l]] = mu[l];
                    order[varAxis[l]] = 2*kappa[l];
                }
                
                Z[ind3(nu,mu[0],mu[1],NU,MU0,MU1)] = der.mixedNumDeriv(Vcopy, 3, wth, Deriv, order, G.dx, lth);
                
                if( (mu[0]!=0) && (mu[0]%2)==0 && kappa[0]>1 ){
                    kappa[0] = kappa[0] - 1;
                }// end of if mu0
            }// end of mu1 loop
            if( (mu[1]!=0) && (mu[1]%2)==0 && kappa[1]>1 ){
                kappa[1] = kappa[1] - 1;
            }// end of if mu1
        }// end of mu1
    }// end of nu loop

    /* free all the pointers allocated in this function */
    for(int k = 1; k<=p; k++){
        delete [] V[ind2(0,k,NU,K)];
        V[ind2(0,k,NU,K)] = NULL;
    }
    for(int nu = 0; nu<=(p-1); nu++){
        for(int k = 1; k<=(p-nu); k++){
            delete [] V[ind2((nu+1),k,NU,K)];
            V[ind2((nu+1),k,NU,K)] = NULL;
        }// end of k loop
    }// end of nu loop
    /* end of pointer freeing */
}// end of Dirichlet_LagrangeDer3 function

void Lcbc::applyQh(double *&QV, double *V, double **W, int *Ind, int nu, int k, int axis, int side){
    
    int face = side + 2*axis;
    int faceNum = 2*dim;
    
    int newWth[3], newLth[3], oldWth[3], oldLth[3];
    getWidth(newWth, newLth, axis, (nu+1), k, p, dim);
    getWidth(oldWth, oldLth, axis, nu, k, p, dim);
    
    /* Get the correction terms and save them in VS */
    double **VS;
    VS = new double*[(maxCoefNum-1)];
    for(int coefNum = 0; coefNum<(maxCoefNum-1); coefNum++){
        VS[coefNum] = new double[(oldLth[0]*oldLth[1]*oldLth[2])];
    }
    getCorrectionTerms(VS,V, W, oldLth, k,axis);
    
    int coefLth[3]; getCoefGridLth(G.indexRange, coefLth, axis, side, dim, p);
    int varAxis[2]; getVarAxis(varAxis, axis);
    int v0 = varAxis[0];
    int v1 = varAxis[1];

    /* Find Qh of V */
    int i[3], order[] = {2,2,2}; int cInd[3];
    for(i[2] = (-newWth[2]); i[2]<=newWth[2]; i[2]++){
        for(i[1] = (-newWth[1]); i[1]<=newWth[1]; i[1]++){
            for(i[0] = (-newWth[0]); i[0]<=newWth[0]; i[0]++){
                cInd[axis] = p + i[axis];
                cInd[v0] = Ind[v0] + p + i[v0];
                cInd[v1] = dimBasedValue(dim, 0, (Ind[v1] + p + i[v1]));
                int vInd[] = sumVectors(oldWth, i);
                int qInd[] = sumVectors(newWth,i);
                
                int cI = ind(cInd,coefLth);
                int qI = ind(qInd,newLth);

                QV[qI] = 0; int coefNum = 0;
                for(int degree = 2; degree>0; degree = degree - 1){
                    for(int d = 0; d<dim; d++){
                        QV[qI] = QV[qI] + coef[ind2(face,coefNum,faceNum,maxCoefNum)][cI]*der.numDeriv(VS[coefNum], dim, vInd, d, degree,2, G.dx[d], oldLth);
                        coefNum++;
                    }
                }

                for(int d = (2*(dim-2)); d>=0; d--){
                    int deriv[] = {1,1,1};
                    (dim>2) ? (deriv[d] = 0) : (deriv[2] = 0);
                    QV[qI] = QV[qI] + 2*coef[ind2(face,coefNum,faceNum,maxCoefNum)][cI]*der.mixedNumDeriv(VS[coefNum], dim, vInd, deriv, order, G.dx, oldLth);
                    coefNum++;
                }

                QV[qI] = QV[qI] + coef[ind2(face,coefNum,faceNum,maxCoefNum)][cI]*V[ind(vInd,oldLth)];

            }// end of i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Free the pointers */
    for(int i = 0; i<(maxCoefNum-1); i++){
        delete [] VS[i];
        VS[i] = NULL;
    }
    delete [] VS;
}


void Lcbc::getCorrectionTerms(double **&VS, double *V, double **W, int *lth, int k, int axis){
    /* Note that 0,1,2,3,4 correspond to c11, c22, c1, c2, c12 correction terms */
#define coefVec(num) ((num < dim) ? (a) : (b))
    
    double a[] = {0,(-1.0/12.0),(1.0/90.0)};
    double b[] = {(1.0),(-1.0/6.0),(1.0/30.0)};
    
    int i[3];
    
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<
            lth[1]; i[1]++){
            for(i[0] = 0; i[0]<lth[0]; i[0]++){
                
                int coefNum = 0;
                while(coefNum<(2*dim)){
                    int alt = (coefNum % dim);
                    int m[] = {0,0,0};
                    VS[coefNum][ind(i,lth)] = 0;
                    for(int n = 1; n<=(k-1); n++){
                        m[alt] = n;
                        VS[coefNum][ind(i,lth)] = VS[coefNum][ind(i,lth)] + coefVec(coefNum)[n]*pow(G.dx[alt],((double) (2.0*n)))*W[ind4((k-n),m[0],m[1],m[2],(p+1),p,p,p)][ind(i,lth)];
                    }// end of n loop
                    VS[coefNum][ind(i,lth)] = VS[coefNum][ind(i,lth)] + V[ind(i,lth)];
                    coefNum++;
                }// end of while coefNum loop
                
                for(int d1 = 0; d1<(dim - 1); d1++){
                    for(int d2 = (d1+1); d2<dim; d2++){
                        VS[coefNum][ind(i,lth)] = 0;
                        int m[] = {0,0,0};

                        for(int n = 1; n<=(k-1); n++){
                            for(int l=0; l<=n; l++){
                                m[d1] = l;
                                m[d2] = (n-l);

                                VS[coefNum][ind(i,lth)] = VS[coefNum][ind(i,lth)] + (b[l]*b[n-l]*(pow(G.dx[d1],((double)(2.0*l))))*(pow(G.dx[d2],((double) 2.0*(n-l)))))*W[ind4((k-n),m[0],m[1],m[2],(p+1),p,p,p)][ind(i,lth)];
                            }// end of l loop
                        }// end of n loop
                        VS[coefNum][ind(i,lth)] = VS[coefNum][ind(i,lth)] + V[ind(i,lth)];
                        coefNum++;
                    }// end of d2 loop
                }// end of d1 loop
            }// end of i[0] loop
        }// end of i[1] loop
    } // end of i[2] loop
    
    
}// end of getCorrectionTerms

void getWidth(int wth[3], int lth[3], int fixedAxis, int nu, int k, int p, int dim){
    for(int axis = 0; axis<3; axis ++){
        if(axis<dim){
            if(axis == fixedAxis){
                wth[axis] = (p+1 - (nu+k));
                
            }else{
                wth[axis] = (2*p+1 - (nu+k));
            }
        }else{
            wth[axis] = 0;
        }
        lth[axis] = (2*wth[axis] + 1);
    }
}

void Lcbc::getDataBdryRange(int bdryRange[3][2], int bdryNgx[3], int fixedAxis, int fixedSide, int addOnSide0, int addOnSide1){
    for(int axis = 0; axis<3; axis++){
        for(int side = 0; side<2; side++){
            if(axis == fixedAxis){
                bdryRange[axis][side] = G.indexRange[fixedAxis][fixedSide];
                bdryNgx[axis] = 1;
            }
            else{
                if(axis<dim){
                    bdryRange[axis][side] = G.indexRange[axis][side] + (1-side)*addOnSide0 + addOnSide1*side;
                }else{
                    bdryRange[axis][side] = G.indexRange[axis][side];
                }
                bdryNgx[axis] = G.Ngx[axis];
            }// end of if axis
        }// end of side
    }// end of axis
}// end of getDataBdryRange

void Lcbc::getBdryRange(int bdryRange[3][2], int bdryNgx[3], int fixedAxis, int fixedSide, int addOnSide0, int addOnSide1){
    for(int axis = 0; axis<3; axis++){
        for(int side = 0; side<2; side++){
            if(axis == fixedAxis){
                bdryRange[axis][side] = G.indexRange[fixedAxis][fixedSide];
                bdryNgx[axis] = 1;
            }
            else{
                if(axis<dim && (faceEval[(2*axis)]>=0)){
                    bdryRange[axis][side] = G.indexRange[axis][side] + (1-side)*addOnSide0 + addOnSide1*side;
                }else{
                    bdryRange[axis][side] = G.indexRange[axis][side];
                }
                bdryNgx[axis] = G.Ngx[axis];
            }// end of if axis
        }// end of side
    }// end of axis
}// end of getBdryRange

int Lcbc::getBdryPointsNum(int fixedAxis){
    int bdryNg = 1;
    for(int axis = 0; axis<3; axis++){
            if(axis != fixedAxis){
                bdryNg = bdryNg*G.Ngx[axis];
            }// end of if axis
    }// end of axis
    return bdryNg;
}// end of getBdryPointsNum
