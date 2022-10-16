#include "LCBC.h"
#include <string.h>
#include <math.h>
#include "utility.h"
#include "LCBCmacros.h"
#include <assert.h>

void getWidth(int (&wth)[3], int (&lth)[3], int fixedAxis, int nu, int k, int p, int dim, int bdryRangeLth[3], int faceType);
void getDataFnLth(int fnLth[3], int fnWth[3], int fixedAxis, int bdryRangeLth[3], int dim, int p, int faceType, int extraDataGhost);

void Lcbc::prepDataVec(double **R, double **&Rv, double t, double dt, int LcbcBdryRange[3][2], int bdryNgx[3], int axis, int side){
    
    int face = side + 2*axis;
    int NU = faceParam[face].NU;
    
    int approxEqNum = faceParam[face].approxEqNum;
    int varAxis[2]  = {faceParam[face].otherAxis[0], faceParam[face].otherAxis[1]};
    
    int i[3], j[3] = {0,0,0};
    for(i[2] = LcbcBdryRange[2][0]; i[2]<=LcbcBdryRange[2][1]; i[2]++){
        for(i[1] = LcbcBdryRange[1][0]; i[1]<=LcbcBdryRange[1][1]; i[1]++){
            for(i[0] = LcbcBdryRange[0][0]; i[0]<=LcbcBdryRange[0][1]; i[0]++){
                int copy = i[axis]; i[axis] = 0;
                int r = 0;
                for(int nu = 0; nu<NU; nu++){
                    for(j[varAxis[1]] = dimBasedValue(dim, 0, (-p)); j[varAxis[1]]<= dimBasedValue(dim, 0, p); j[varAxis[1]]++){
                        for(j[varAxis[0]] = -p; j[varAxis[0]]<=(p); j[varAxis[0]]++){
                            int Rind[] = sumVectors(i, j);
                            Rv[ind(i,bdryNgx)][r] =  R[ind2(nu,face,(p+1),faceNum)][ind(Rind,bdryNgx)];
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

void Lcbc::getFaceGhost(double *&un, double *Rv[], double ***CaVec, double ***CbVec, int *eqNum, int LcbcBdryRange[3][2], int bdryNgx[3], int axis, int side){
    
    /* prepare the parameters */
    int face = side + 2*axis;
    int approxEqNum = faceParam[face].approxEqNum;
    int interiorEqNum = param.interiorEqNum;
    
    double b[interiorEqNum];
    
    /* Evaluate the solution at the ghost points */
    int i[3];
    for(i[2] = LcbcBdryRange[2][0]; i[2]<=LcbcBdryRange[2][1]; i[2]++){
        for(i[1] = LcbcBdryRange[1][0]; i[1]<= LcbcBdryRange[1][1]; i[1]++){
            for(i[0] = LcbcBdryRange[0][0]; i[0]<=LcbcBdryRange[0][1]; i[0]++){
                
                getbInt(b, un, eqNum, i, faceParam[face].interiorRange, axis, side);
                
                int copy = i[axis]; i[axis] = 0;
                int bdInd = ind(i,bdryNgx);
                int cVecInd = (cstCoef)?(0):(bdInd);
                
                int vecNum = 0;
                for(int ghostInd = sideBasedValue(side, (-p), 1); ghostInd<=sideBasedValue(side, (-1), p); ghostInd++){

                    double Ca_vecRv = dotProduct(CaVec[cVecInd][vecNum], Rv[bdInd], approxEqNum);
                    double Cb_vecb  = dotProduct(CbVec[cVecInd][vecNum], b, interiorEqNum);

                    i[axis] = (copy + ghostInd);
                    un[solInd(i,G.Ngx)] = Ca_vecRv + Cb_vecb;

                    vecNum++; 
                }// end of ghostVal
                i[axis] = copy;
            }// end of i[0]
        }// end of i[1]
    }// end of i[2]
}// end of getFaceGhost

void Lcbc::getDataDeriv(double **&R, double t, double dt, int axis, int side, LcbcData &gn, LcbcData &fn){
    int face = ind2(side,axis,2,3);
    
    int faceNum = param.faceNum;
    int NU = faceParam[face].NU;
    
    int K = (p+1);
    int R_NU = (p+1);
    
    int faceType = faceEval[face];
    
    for(int nu = 0; nu<NU; nu++){
        R[ind2(nu,face,R_NU,faceNum)] = new double[(faceParam[face].bdryNgx[0]*faceParam[face].bdryNgx[1]*faceParam[face].bdryNgx[2])];
    }
    if(faceType == 2){ // Neumann BC case
        R[ind2(p,face,R_NU,faceNum)] = NULL;
    }
    
    int i[3];
    for(int nu = 0; nu<NU; nu++){
        int derivOrder = (nu == 0)?(2*p):(2*(NU - nu));
        for(i[2] = faceParam[face].bdryRangeExt[2][0]; i[2]<=faceParam[face].bdryRangeExt[2][1]; i[2]++){
            for(i[1] = faceParam[face].bdryRangeExt[1][0]; i[1]<=faceParam[face].bdryRangeExt[1][1]; i[1]++){
                for(i[0] = faceParam[face].bdryRangeExt[0][0]; i[0]<=faceParam[face].bdryRangeExt[0][1]; i[0]++){

                    int copy = i[axis];
                    
                    if(zeroBC[face]){
                        i[axis] = 0;
                        R[ind2(nu,face,R_NU,faceNum)][ind(i, faceParam[face].bdryNgx)] = 0;
                    }else{
                        if(!(gn.initialized)){
                            double arg[] = {((double)(side + 2*axis)), G.x[0][i[0]], G.x[1][i[1]], G.x[2][i[2]], t};
                            arg[(dim+1)] = t;
                            
                            i[axis] = 0;
                            R[ind2(nu,face,R_NU,faceNum)][ind(i, faceParam[face].bdryNgx)] = Deriv(Gn, arg, dim, (q*nu), dt, derivOrder);
                        }else{
                            i[axis] = 0;
                            int Ind[3] = {(gn.wth[0]-p + i[0]),(gn.wth[1]-p + i[1]),dimBasedValue(dim, 0, (gn.wth[2]-p + i[2]))}; Ind[axis] = 0;
                            R[ind2(nu,face,R_NU,faceNum)][ind(i, faceParam[face].bdryNgx)] = gn.Fn[nu][ind(Ind,gn.lth)];
                        }
                    }
                    
                    i[axis] = copy;
                }// end of i[0] loop s
            }// end of i[1] loop
        }// end of i[2] loop
    }// end of nu loop
    
    if(!noForcing){
    
    double *F[NU*K];
    double *W[(K*p*p*p)];
    
    /*--------------------------------------------------*/

    int wth[3], lth[3];
    for(int n = 0; n<=(p-1); n++){
        int derivOrder = (n == 0)?(2*p):(2*(NU - n));
        for(int k = 1; k<=p; k++){
            getWidth(wth, lth, axis, 0, k, p, dim, faceParam[face].bdryRangeLth, faceType);
            F[ind2(0,k,NU,K)] = new double[(lth[0]*lth[1]*lth[2])];
            
            int i[3];
            for(i[2] = 0; i[2]<lth[2]; i[2]++){
                for(i[1] = 0; i[1]<lth[1]; i[1]++){
                    for(i[0] = 0; i[0]<lth[0]; i[0]++){
                        
                        if(!(fn.initialized)){
                            double arg[] = {0,G.x[0][faceParam[face].bdryRange[0][0]] + ((- wth[0] + i[0])*G.dx[0]),
                                              G.x[1][faceParam[face].bdryRange[1][0]] + ((- wth[1] + i[1])*G.dx[1]),
                                              G.x[2][faceParam[face].bdryRange[2][0]] + ((- wth[2] + i[2])*G.dx[2]),t};
                            arg[(dim+1)] = t;
                            F[ind2(0,k,NU,K)][ind(i,lth)] = Deriv(Fn, arg, dim, (q*n), dt, derivOrder);
                        }
                        else{
                            
                            int Ind[3] = {(fn.wth[0]-wth[0] + i[0]),(fn.wth[1]-wth[1] + i[1]),dimBasedValue(dim, 0, (fn.wth[2]-wth[2] + i[2]))};

                            F[ind2(0,k,NU,K)][ind(i,lth)] = fn.Fn[n][ind(Ind,fn.lth)];
                        }
                    }// end of i[0]
                }// end of i[1]
            }// end of i[2]
        }// end of k loop
        
        int ord[] = {2,2,dimBasedValue(dim, 0, 2)};
        int W_wth[3], W_lth[3];
        for(int nub = 0; nub<=(p-n-2); nub++){
            for(int k = 1; k<=(p-nub-1); k++){
                if(k!= 1){
                    getWidth(W_wth, W_lth, axis, nub, k, p, dim, faceParam[face].bdryRangeLth, faceType);
                    for(int l = 1; l<=(k-1); l++){
                        getWidth(wth, lth, axis, nub, l, p, dim, faceParam[face].bdryRangeLth, faceType);
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
                
                                            W[ind4(l,Dx,Dy,Dz,K,p,p,p)][ind(i,W_lth)] = mixedDeriv(F[ind2(nub,l,NU,K)], F_ind, Deriv, G.dx, ord, lth);
                                            
                                        }// end of i0 loop
                                    }// end of i1 loop
                                }// end of i2 loop
                            }// end of Dy loop
                        }// end of Dx loop
                    }// end of l loop
                }// end of if k statement
                
                /* Here's where we apply Qf */
                int oldWth[3], oldLth[3];
                getWidth(oldWth, oldLth, axis, nub, k, p, dim, faceParam[face].bdryRangeLth, faceType);
                
                getWidth(wth, lth, axis, (nub+1), k, p, dim, faceParam[face].bdryRangeLth, faceType);
                F[ind2((nub+1),k,NU,K)] = new double[(lth[0]*lth[1]*lth[2])];
                
                if(k<(p - nub - 1) || faceType == 1){
                    applyQhF(F[ind2((nub+1),k,NU,K)], F[ind2(nub,k,NU,K)], W, faceParam[face].bdryRange, faceParam[face].bdryRangeLth, nub, k, wth, lth, oldWth, oldLth, axis, side);
                }else{
                    applyDhF(F[ind2((nub+1),k,NU,K)], F[ind2(nub,k,NU,K)], W, faceParam[face].bdryRange, faceParam[face].bdryRangeLth, nub, k, wth, lth, oldWth, oldLth, axis, side);
                }
                
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
        for(int nu = (n+1); nu<=(NU-1); nu++){
            
            int nub = nu - n - 1, k;
            if(faceEval[face] == 1){
                k = MIN((p-nub),p);
            }else{
                k = MIN((p-1-nub),p);
            }
            if(faceType == 1){
                getWidth(wth, lth, axis, nub, k, p, dim, faceParam[face].bdryRangeLth, faceType); // length for F
            }else if(faceType == 2){
                getWidth(wth, lth, axis, (nub+1), k, p, dim, faceParam[face].bdryRangeLth, faceType); // length for F
            }
            getWidth(R_wth, R_lth, axis, (p-1), 1, p, dim, faceParam[face].bdryRangeLth, 1); // the final needed lth and wth
            
            int i[3];
            for(i[2] = 0; i[2]<R_lth[2]; i[2]++){
                for(i[1] = 0; i[1]<R_lth[1]; i[1]++){
                    for(i[0] = 0; i[0]<R_lth[0]; i[0]++){
                        int Rind[] = {(faceParam[face].bdryRange[0][0] - R_wth[0] + i[0]),
                                      (faceParam[face].bdryRange[1][0] - R_wth[1] + i[1]),
                                      (faceParam[face].bdryRange[2][0] - R_wth[2] + i[2])};
                        Rind[axis] = 0;
                        int RI = ind(Rind,faceParam[face].bdryNgx);
                        int Find[] = {(wth[0] - R_wth[0] + i[0]),
                                      (wth[1] - R_wth[1] + i[1]),
                                      (wth[2] - R_wth[2] + i[2])};
                        int FI = ind(Find,lth);
                        if(faceEval[face] == 1){
                            R[ind2(nu,face,R_NU,faceNum)][RI] = R[ind2(nu,face,R_NU,6)][RI] - F[ind2(nub,k,NU,K)][FI];
                        }
                        else if(faceEval[face] == 2){
                            R[ind2(nu,face,R_NU,faceNum)][RI] = R[ind2(nu,face,R_NU,6)][RI] - F[ind2((nub+1),k,NU,K)][FI];
                        }
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
            for(int k = 1; k<=(p-nub-1); k++){
                delete [] F[ind2((nub+1),k,NU,K)];
                F[ind2((nub+1),k,NU,K)] = NULL;
            }// end of k loop
        }// end of nub loop
        
    }// end of n loop
    }// end of if(!noForcing)
}// end of getDataDeriv function

void Lcbc::applyQhF(double *&QF, double *F, double **W, int bdryRange[3][2], int bdryRangeLth[3], int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side){
    
    int face = side + 2*axis;
    
    /* Get the correction terms and save them in FS */
    double **FS;
    FS = new double*[(maxCoefNum - 1)];
    for(int coefNum = 0; coefNum<(maxCoefNum - 1); coefNum++){
        FS[coefNum] = new double[(oldLth[0]*oldLth[1]*oldLth[2])];
    }// end of coefNum loop
    getCorrectionTerms(FS,F, W, oldLth, k,axis);
    
    int varAxis[2]; getVarAxis(varAxis, axis);
    int v0 = varAxis[0];
    int v1 = varAxis[1];
    int cI = 0;
    
    int i[3], order[] = {2,2,2}; int cInd[3];
    for(i[2] = 0; i[2]<newLth[2]; i[2]++){
        for(i[1] = 0; i[1]<newLth[1]; i[1]++){
            for(i[0] = 0; i[0]<newLth[0]; i[0]++){
            
                if(!cstCoef){
                    cInd[axis] = coef[face].wth[axis] - newWth[axis] + i[axis];
                    cInd[v0]   = (bdryRange[v0][0] - p + coef[face].wth[v0]) - newWth[v0] + i[v0];
                    cInd[v1]   = dimBasedValue(dim, 0, (bdryRange[v1][0] - p + coef[face].wth[v1] - newWth[v1] + i[v1]));
                    cI = ind(cInd, coef[face].lth);
                }
 
                int FInd[] = {(oldWth[0] - newWth[0] + i[0]),
                              (oldWth[1] - newWth[1] + i[1]),
                              (oldWth[2] - newWth[2] + i[2])};
                

                int FI  = ind(FInd, oldLth);
                int I = ind(i,newLth);
                
                QF[I] = 0; int coefNum = 0;
                for(int degree = 2; degree>0; degree = degree - 1){
                    for(int d = 0; d<dim; d++){
                        
                        QF[I] = QF[I] + coef[face].Fn[coefNum][cI]*Deriv(FS[coefNum], FInd, d, degree, G.dx[d], 2, oldLth);

                        coefNum++;
                    }// end of d loop
                }// end of degree loop
                
                for(int d = (2*(dim-2)); d>=0; d--){
                    int deriv[] = {1,1,1};
                    (dim>2) ? (deriv[d] = 0) : (deriv[2] = 0);
                    
                    QF[I] = QF[I] + 2*coef[face].Fn[coefNum][cI]*mixedDeriv(FS[coefNum], FInd, deriv, G.dx, order, oldLth);
                    
                    coefNum++;
                }// end of d loop
                
                QF[I] = QF[I] + coef[face].Fn[coefNum][cI]*F[FI];
                
            }// end of i[0]
        }// end of i[1]
    }// end of i[2]
    
    for(int coefNum = 0; coefNum<(maxCoefNum - 1); coefNum++){
        delete [] FS[coefNum]; FS[coefNum] = NULL;
    }// end of coefNum loop
    delete [] FS;
}// end of applyQhF function

void Lcbc::applyDhF(double *&QF, double *F, double **W, int bdryRange[3][2], int bdryRangeLth[3], int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side){
    
    int face = side + 2*axis;

    /* Get the correction terms and save them in FS */
    double **FS;
    FS = new double*[dim];
    for(int axisNum = 0; axisNum<dim; axisNum++){
        FS[axisNum] = new double[(oldLth[0]*oldLth[1]*oldLth[2])];
    }// end of coefNum loop
    getCorrectionTermsForDh(FS,F, W, oldLth, k,axis);
    
    int i[3];
    for(i[2] = 0; i[2]<newLth[2]; i[2]++){
        for(i[1] = 0; i[1]<newLth[1]; i[1]++){
            for(i[0] = 0; i[0]<newLth[0]; i[0]++){
        
                int FInd[] = {(oldWth[0] - newWth[0] + i[0]),
                              (oldWth[1] - newWth[1] + i[1]),
                              (oldWth[2] - newWth[2] + i[2])};
                
                int I = ind(i,newLth);

                double arg[] = {(double)face, (G.x[0][bdryRange[0][0]] + (i[0]-newWth[0])*G.dx[0]), (G.x[1][bdryRange[1][0]] + (i[1]-newWth[1])*G.dx[1]), dimBasedValue(dim,0,(G.x[2][bdryRange[2][0]] + (i[2]-newWth[2])*G.dx[2]))};
                double bn_vec[3];
                bn(bn_vec, arg);
            
                QF[I] =  bn_vec[0]*Deriv(FS[0], FInd, 0, 1, G.dx[0], 2, oldLth)
                       + bn_vec[1]*Deriv(FS[1], FInd, 1, 1, G.dx[1], 2, oldLth);

                if(dim == 3){
                    QF[I] = QF[I] + bn_vec[2]*Deriv(FS[2], FInd, 2, 1, G.dx[2], 2, oldLth);
                }// end of if dim = 3
                
            }// end of i[0]
        }// end of i[1]
    }// end of i[2]
    
    for(int axisNum = 0; axisNum<dim; axisNum++){
        delete [] FS[axisNum]; FS[axisNum] = NULL;
    }// end of coefNum loop
    delete [] FS;
}// end of applyQhF function

void getWidth(int (&wth)[3], int (&lth)[3], int fixedAxis, int nu, int k, int p, int dim, int bdryRangeLth[3], int faceType){
    
    for(int axis = 0; axis<3; axis ++){
        if(axis<dim){
            if(axis == fixedAxis){
                wth[axis] = (p - (nu+k));
            }else{
                wth[axis] = (2*p - (nu+k));
            }
        }else{
            wth[axis] = 0;
        }
        lth[axis] = (2*wth[axis] + bdryRangeLth[axis]);
    }
}

void getDataFnLth(int fnLth[3], int fnWth[3], int fixedAxis, int bdryRangeLth[3], int dim, int p, int faceType, int extraDataGhost){
    
    for(int axis = 0; axis<3; axis ++){
        if(axis<dim){
            if(axis == fixedAxis){
                fnWth[axis] = (p);
            }else{
                fnWth[axis] = (p + extraDataGhost);
            }
        }else{
            fnWth[axis] = 0;
        }
        fnLth[axis] = (2*fnWth[axis] + bdryRangeLth[axis]);
    }
}// end of getDataFnLth
