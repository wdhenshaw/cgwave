#include "LCBC.h"
#include <string.h>
#include <math.h>
#include "utility.h"
#include "LCBCmacros.h"
#include <assert.h>

void getWidth(int (&wth)[3], int (&lth)[3], int fixedAxis, int nu, int k, int p, int dim, int bdryRangeLth[3]);

void Lcbc::prepDataVec(double **R, double **&Rv, double t, double dt, int axis, int side){
    int NU = (p+1);
    int n  = (2*p+1);
    int face = ind2(side,axis,2,3);
    int faceNum = 2*dim;
    
    int bdryRange[3][2], bdryNgx[3]; getBdryRange(bdryRange, bdryNgx, axis, side, (+p), (-p));
    
    int auxiliaryEqNum = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
    int compCondNum     = (p+1)*n*dimBasedValue(dim,1,n);
    int approxEqNum     = compCondNum + auxiliaryEqNum;
    int varAxis[2]; getVarAxis(varAxis, axis);
    
    int i[3], j[3] = {0,0,0};
    for(i[2] = bdryRange[2][0]; i[2]<=bdryRange[2][1]; i[2]++){
        for(i[1] = bdryRange[1][0]; i[1]<=bdryRange[1][1]; i[1]++){
            for(i[0] = bdryRange[0][0]; i[0]<=bdryRange[0][1]; i[0]++){
                int copy = i[axis]; i[axis] = 0;
                int r = 0;
                for(int nu = 0; nu<NU; nu++){
                    for(j[varAxis[1]] = dimBasedValue(dim, 0, (-p)); j[varAxis[1]]<= dimBasedValue(dim, 0, p); j[varAxis[1]]++){
                        for(j[varAxis[0]] = -p; j[varAxis[0]]<=(p); j[varAxis[0]]++){
                            int Rind[] = sumVectors(i, j);
                            Rv[ind(i,bdryNgx)][r] =  R[ind2(nu,face,NU,faceNum)][ind(Rind,bdryNgx)];
                            r++;
                        }// end of k loop
                    }// end of j loop
                }// end of nu loop
            
                for(int row = r; row<approxEqNum; row++){
                    Rv[ind(i,bdryNgx)][row] = 0;
                }
            
                i[axis] = copy;
            }// end of i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
}// end of prepDataVec function

void Lcbc::getFaceGhost(double *&un, double *Rv[], double ***CaVec, double ***CbVec, int *eqNum, int axis, int side){
    /* prepare the parameters */
    int n = 2*p+1;
    int compCondNum = (p+1)*n*dimBasedValue(dim, 1, n);
    int auxiliaryEqNum = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
    int approxEqNum = compCondNum + auxiliaryEqNum;
    int interiorEqNum = (p+1)*n*dimBasedValue(dim, 1,n);
    int bdryRange[3][2], bdryNgx[3];
    getBdryRange(bdryRange, bdryNgx, axis, side, p, (-p));
    
    double b[interiorEqNum];

    /* Evaluate the solution at the ghost points */
    int i[3];
    for(i[2] = bdryRange[2][0]; i[2]<=bdryRange[2][1]; i[2]++){
        for(i[1] = bdryRange[1][0]; i[1]<= bdryRange[1][1]; i[1]++){
            for(i[0] = bdryRange[0][0]; i[0]<=bdryRange[0][1]; i[0]++){
                getbInt(b, un, eqNum, i, axis, side);
                                
                int vecNum = 0;
                for(int ghostInd = sideBasedValue(side, (-p), 1); ghostInd<=sideBasedValue(side, (-1), p); ghostInd++){
                    int copy = i[axis]; i[axis] = 0;
                    
                    double Ca_vecRv = dotProduct(CaVec[ind(i,bdryNgx)][vecNum], Rv[ind(i,bdryNgx)], approxEqNum);
                    double Cb_vecb  = dotProduct(CbVec[ind(i,bdryNgx)][vecNum], b, interiorEqNum);

                    i[axis] = (copy + ghostInd);
                    un[solInd(i,G.Ngx)] = Ca_vecRv + Cb_vecb;
                    i[axis] = copy;
                    
                    vecNum++; 
                }// end of ghostVal
            }// end of i[0]
        }// end of i[1]
    }// end of i[2]
}// end of getFaceGhost

void Lcbc::getbInt(double b[], double *un, int *EqNum, int *Ind, int axis, int side){
    int n = (2*p+1);
    int interiorEqNum = (p+1)*n*dimBasedValue(dim, 1, n);
    int totalVarNum   = (n*n*dimBasedValue(dim, 1, n));
    int unknownVarNum = totalVarNum - interiorEqNum;
    
    int interiorRange[3][2] = {{0,(2*p)},{0,(2*p)},{0,dimBasedValue(dim,0,(2*p))}};
    interiorRange[axis][0] = sideBasedValue(side, p, 0); // if side = 0, choose p else choose 0
    interiorRange[axis][1] = sideBasedValue(side, (2*p), p);
    int eqNumLth[] = {n,n,dimBasedValue(dim,1,n)}, intInd[3];

    for(intInd[0] = interiorRange[0][0]; intInd[0]<=interiorRange[0][1]; intInd[0]++){
        for(intInd[1] = interiorRange[1][0]; intInd[1]<=interiorRange[1][1]; intInd[1]++){
            for(intInd[2] = interiorRange[2][0]; intInd[2]<=interiorRange[2][1]; intInd[2]++){
                int r  = EqNum[ind(intInd,eqNumLth)] - unknownVarNum;
               
                int u_ind[] = {(Ind[0] - p + intInd[0]),(Ind[1] - p + intInd[1]),dimBasedValue(dim, 0, (Ind[2] - p + intInd[2]))};
                
                b[r]  =  un[solInd(u_ind,G.Ngx)];
            }// intInd[0]
        }// intInd[1]
    }// intInd[2]
} // end of getbInt function

void Lcbc::getDataDeriv(double **&R, double t, double dt, int axis, int side, double **gn, double **fn){
    int face = ind2(side,axis,2,3);
    int faceNum = 2*dim;
    int bdryRange[3][2], bdryNgx[3]; getDataBdryRange(bdryRange, bdryNgx, axis, side, (-p), (+p));
    int NU = (p+1), K = (p+1);
    
    for(int nu = 0; nu<NU; nu++){
        R[ind2(nu,face,NU,faceNum)] = new double[(bdryNgx[0]*bdryNgx[1]*bdryNgx[2])];
    }
    
    double *F[NU*K];
    double *W[(K*p*p*p)];
    
    int i[3];
    for(int nu = 0; nu<NU; nu++){
        for(i[2] = bdryRange[2][0]; i[2]<=bdryRange[2][1]; i[2]++){
            for(i[1] = bdryRange[1][0]; i[1]<=bdryRange[1][1]; i[1]++){
                for(i[0] = bdryRange[0][0]; i[0]<=bdryRange[0][1]; i[0]++){

                    int copy = i[axis];
                    if(gn == NULL){
                        double arg[] = {((double)(side + 2*axis)), G.x[0][i[0]], G.x[1][i[1]], G.x[2][i[2]], t};
                        arg[(dim+1)] = t;
                        
                        i[axis] = 0;
                        R[ind2(nu,face,NU,faceNum)][ind(i, bdryNgx)] = der.numDerivFn(Gn, arg, (dim+1), (q*nu), orderInTime, dt);
                    }else{
                        i[axis] = 0;
                        R[ind2(nu,face,NU,faceNum)][ind(i, bdryNgx)] = gn[ind2(face,nu,faceNum,NU)][bdryInd(i,bdryNgx,axis)];
                    }
                    i[axis] = copy;
                }// end of i[0] loop s
            }// end of i[1] loop
        }// end of i[2] loop
    }// end of nu loop
    
    int bdryRangeLth[3];
    getBdryRange(bdryRange, bdryRangeLth, axis, side);
    
    
    int fnLth[3];
    
    if(!(fn==NULL)){
    fnLth[0] = (bdryRangeLth[0]  + 4*p);
    fnLth[1] = (bdryRangeLth[1]  + 4*p);
    fnLth[2] = dimBasedValue(dim, 1, (bdryRangeLth[2]  + 4*p));
    fnLth[axis] = (bdryRangeLth[axis]  + 2*p);
    }

    int wth[3], lth[3];
    for(int n = 0; n<=(p-1); n++){
        for(int k = 1; k<=p; k++){
            getWidth(wth, lth, axis, 0, k, p, dim, bdryRangeLth);
            assert((lth[0]*lth[1]*lth[2])>0);
            F[ind2(0,k,NU,K)] = new double[(lth[0]*lth[1]*lth[2])];
            
            int i[3];
            for(i[2] = 0; i[2]<lth[2]; i[2]++){
                for(i[1] = 0; i[1]<lth[1]; i[1]++){
                    for(i[0] = 0; i[0]<lth[0]; i[0]++){
                        
                        if(fn == NULL){
                            double arg[] = {0,G.x[0][bdryRange[0][0]] + ((- wth[0] + i[0])*G.dx[0]),
                                              G.x[1][bdryRange[1][0]] + ((- wth[1] + i[1])*G.dx[1]),
                                              G.x[2][bdryRange[2][0]] + ((- wth[2] + i[2])*G.dx[2]),t};
                            arg[(dim+1)] = t;
                            
                            F[ind2(0,k,NU,K)][ind(i,lth)] = der.numDerivFn(Fn, arg, (dim+1), (q*n), orderInTime, dt);
                        }else{
                            int Ind[3] = {(2*p-wth[0] + i[0]),(2*p-wth[1] + i[1]),dimBasedValue(dim, 0, (2*p-wth[2] + i[2]))};
                            Ind[axis] = (p - wth[axis] + i[axis]);
                            
                            F[ind2(0,k,NU,K)][ind(i,lth)] = fn[ind2(face,n,faceNum,p)][ind(Ind,fnLth)];
                        }
                    }// end of i[0]
                }// end of i[1]
            }// end of i[2]
        }// end of k loop
        
        int ord[] = {2,2,dimBasedValue(dim, 0, 2)};
        int W_wth[3], W_lth[3];
        for(int nub = 0; nub<=(p-n-2); nub++){
            for(int k = 1; k<=(p-nub); k++){
                if(k!= 1){
                    getWidth(W_wth, W_lth, axis, nub, k, p, dim, bdryRangeLth);
                    for(int l = 1; l<=(k-1); l++){
                        getWidth(wth, lth, axis, nub, l, p, dim, bdryRangeLth);
                        int m = (k - l);
                        for(int Dx = 0; Dx<=m; Dx++){
                            for(int Dy = dimBasedValue(dim, (m - Dx), 0); Dy<=(m - Dx); Dy++){
                                int Dz = dimBasedValue(dim, 0, (m - Dx - Dy));
                                W[ind4(l,Dx,Dy,Dz,K,p,p,p)] = new double[(W_lth[0]*W_lth[1]*W_lth[2])];

                                int Deriv[] = {(2*Dx),(2*Dy),(2*Dz)};
                                int i[3];
                                for(i[2] = 0; i[2]<W_lth[2]; i[2]++){
                                    for(i[1] = 0; i[1]<W_lth[1]; i[1]++){
                                        for(i[0] = 0; i[0]<W_lth[0]; i[0]++){
                                            int F_ind[] = {(wth[0] - W_wth[0]+ i[0]),
                                                           (wth[1] - W_wth[1]+ i[1]),
                                                           (wth[2] - W_wth[2]+ i[2])};
                
                                            W[ind4(l,Dx,Dy,Dz,K,p,p,p)][ind(i,W_lth)] = der.mixedNumDeriv(F[ind2(nub,l,NU,K)], 3, F_ind, Deriv, ord, G.dx, lth);
                                        }// end of i0 loop
                                    }// end of i1 loop
                                }// end of i2 loop
                            }// end of Dy loop
                        }// end of Dx loop
                    }// end of l loop
                }// end of if k statement
                
                
                /* Here's where we apply Qf */
                getWidth(wth, lth, axis, (nub+1), k, p, dim, bdryRangeLth);
                assert((lth[0]*lth[1]*lth[2])>0);
                F[ind2((nub+1),k,NU,K)] = new double[(lth[0]*lth[1]*lth[2])];
                
                applyQhF(F[ind2((nub+1),k,NU,K)], F[ind2(nub,k,NU,K)], W, bdryRange, bdryRangeLth, nub, k, axis, side);

                /* Free the W variable */
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
        }// end of nub
        
        /* Here's where we finalize the construction of R */
        int R_wth[3], R_lth[3];
        for(int nu = (n+1); nu<=p; nu++){
            int nub = nu - n - 1;
            int k = MIN((p+1-nub),p);
            getWidth(wth, lth, axis, nub, k, p, dim, bdryRangeLth);
            getWidth(R_wth, R_lth, axis, p, 1, p, dim, bdryRangeLth);
            
            int i[3];
            for(i[2] = 0; i[2]<R_lth[2]; i[2]++){
                for(i[1] = 0; i[1]<R_lth[1]; i[1]++){
                    for(i[0] = 0; i[0]<R_lth[0]; i[0]++){
                        int Rind[] = {(bdryRange[0][0] - R_wth[0] + i[0]),
                                      (bdryRange[1][0] - R_wth[1] + i[1]),
                                      (bdryRange[2][0] - R_wth[2] + i[2])};
                        Rind[axis] = 0;
                        int RI = ind(Rind,bdryNgx);
                        int Find[] = {(wth[0] - R_wth[0] + i[0]),
                                      (wth[1] - R_wth[1] + i[1]),
                                      (wth[2] - R_wth[2] + i[2])};
                        int FI = ind(Find,lth);
                        
                        R[ind2(nu,face,NU,faceNum)][RI] = R[ind2(nu,face,NU,6)][RI] - F[ind2(nub,k,NU,K)][FI];
                    }// end of i[0]
                }// end of i[1]
            }// end of i[2]
        }// end of nu loop
        
        /* free the space from F */
        for(int k = 1; k<=p; k++){
            delete [] F[ind2(0,k,NU,K)];
            F[ind2(0,k,NU,K)] = NULL;
        }// end of k loop
        for(int nub = 0; nub<=(p-n-2); nub++){
            for(int k = 1; k<=(p-nub); k++){
                delete [] F[ind2((nub+1),k,NU,K)];
                F[ind2((nub+1),k,NU,K)] = NULL;
            }// end of k loop
        }// end of nub loop
        
    }// end of n loop
}// end of getDataDeriv function

void Lcbc::applyQhF(double *&QF, double *F, double **W, int bdryRange[3][2], int bdryRangeLth[3], int nu, int k, int axis, int side){
    
    int face = side + 2*axis;
    int faceNum = 2*dim;
    
    int newWth[3], newLth[3], oldWth[3], oldLth[3];
    getWidth(newWth, newLth, axis, (nu+1), k, p, dim, bdryRangeLth);
    getWidth(oldWth, oldLth, axis, nu, k, p, dim, bdryRangeLth);
    
    /* Get the correction terms and save them in FS */
    double **FS;
    FS = new double*[(maxCoefNum - 1)];
    for(int coefNum = 0; coefNum<(maxCoefNum - 1); coefNum++){
        FS[coefNum] = new double[(oldLth[0]*oldLth[1]*oldLth[2])];
    }// end of coefNum loop
    getCorrectionTerms(FS,F, W, oldLth, k,axis);
    
    int coefLth[3]; getCoefGridLth(G.indexRange, coefLth, axis, side, dim, p);
    int varAxis[2]; getVarAxis(varAxis, axis);
    int v0 = varAxis[0];
    int v1 = varAxis[1];
    
    int i[3], order[] = {2,2,2}; int cInd[3];
    for(i[2] = 0; i[2]<newLth[2]; i[2]++){
        for(i[1] = 0; i[1]<newLth[1]; i[1]++){
            for(i[0] = 0; i[0]<newLth[0]; i[0]++){
                
                cInd[axis] = p - newWth[axis] + i[axis];
                cInd[v0]   = bdryRange[v0][0] + p - newWth[v0] + i[v0];
                cInd[v1]   = dimBasedValue(dim, 0, (bdryRange[v1][0] + p - newWth[v1] + i[v1]));
 
                int FInd[] = {(oldWth[0] - newWth[0] + i[0]),
                              (oldWth[1] - newWth[1] + i[1]),
                              (oldWth[2] - newWth[2] + i[2])};
                
                int cI  = ind(cInd,coefLth);
                int FI  = ind(FInd, oldLth);
                int I = ind(i,newLth);
                
                QF[I] = 0; int coefNum = 0;
                for(int degree = 2; degree>0; degree = degree - 1){
                    for(int d = 0; d<dim; d++){
                        QF[I] = QF[I] + coef[ind2(face,coefNum,faceNum,maxCoefNum)][cI]*der.numDeriv(FS[coefNum], dim, FInd, d, degree,2, G.dx[d], oldLth);
                        coefNum++;
                    }// end of d loop
                }// end of degree loop
                
                for(int d = (2*(dim-2)); d>=0; d--){
                    int deriv[] = {1,1,1};
                    (dim>2) ? (deriv[d] = 0) : (deriv[2] = 0);
                    QF[I] = QF[I] + 2*coef[ind2(face,coefNum,faceNum,maxCoefNum)][cI]*der.mixedNumDeriv(FS[coefNum], dim, FInd, deriv, order, G.dx, oldLth);
                    coefNum++;
                }// end of d loop
                
                QF[I] = QF[I] + coef[ind2(face,coefNum,faceNum,maxCoefNum)][cI]*F[FI];
                
            }// end of i[0]
        }// end of i[1]
    }// end of i[2]
    
    for(int coefNum = 0; coefNum<(maxCoefNum - 1); coefNum++){
        delete [] FS[coefNum]; FS[coefNum] = NULL;
    }// end of coefNum loop
    delete [] FS;
}// end of applyQhF function

void Lcbc::getBdryRange(int (&bdryRange)[3][2], int (&bdryRangeLth)[3], int fixedAxis, int fixedSide){
    for(int axis = 0; axis<3; axis++){
        for(int side = 0; side<2; side++){
            if(axis == fixedAxis){
                bdryRange[axis][side] = G.indexRange[fixedAxis][fixedSide];
            }
            else{
                if(axis<dim){
                    bdryRange[axis][side] = G.indexRange[axis][side];
                }else{
                    bdryRange[axis][side] = G.indexRange[axis][side];
                }
            }// end of if axis
        }// end of side
        bdryRangeLth[axis] = bdryRange[axis][1] - bdryRange[axis][0] + 1;
    }// end of axis
}// end of getBdryRange

void getWidth(int (&wth)[3], int (&lth)[3], int fixedAxis, int nu, int k, int p, int dim, int bdryRangeLth[3]){
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
        lth[axis] = (2*wth[axis] + bdryRangeLth[axis]);
    }
}
