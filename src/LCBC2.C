#include "LCBC.h"
#include <string.h>
#include <math.h>
#include "utility.h"
#include "LCBCmacros.h"
#include <assert.h>

/* Note: repeated functions parameters may be documented once when they first appear */

void getWidth(int (&wth)[3], int (&lth)[3], int fixedAxis, int nu, int k, int p, int dim, int bdryRangeLth[3], int faceType);

/// \brief prepare the data vector R(t). See the stencil approach in the lcbc.pdf file page 20
/// \param R (input): a vector carrying  boundary and forcing values from the right-hand-side of each BC and primary CBC on each boundary point
/// \param Rv (output): the vector described in the lcbc.pdf document page 20. It contains the values in the vector R arranged in the needed order
/// \param t (input): time
/// \param dt (input): time-step
/// \param bdryRange (input): the boundary range over which the data vector needs to be prepared
/// \param bdryNgx (input): the number of grid points at the specific boundary face in each dimension
/// \param bdryParam (input): the type of boundary needed for the LCBC evaluation (face, edge, corner,... see LCBC.h for more information)
/// \param fixedAxis (input): a vector containing the fixed axes
/// \param fixedSide (input): a vector containing the fixed sides
/// \param approxEqNum (input): the number of equations that are applied approximately in the LCBC procedure via least-squares
/// \param NU (input): a vector containing the number of primary CBCs at the faces corresponding to each of the fixed axes
void Lcbc::prepDataVec(double **R, double **&Rv, double t, double dt, int bdryRange[3][2], int bdryNgx[3], bdryTypeParam &bdryParam, int fixedAxis[], int fixedSide[], int approxEqNum, int NU[]){
    
    /* prepare values needed to retrieve the data derivative information */
    int knownAxesNum = bdryParam.knownAxesNum;
    int faceBdryPts[knownAxesNum][3];
    for(int face = 0; face<knownAxesNum; face++){
        memcpy(faceBdryPts[face], G.Ngx, sizeof(faceBdryPts[face]));
        faceBdryPts[face][fixedAxis[face]] = 1;
    }// end of face loop
    
    /* Assign values of R to Rv */
    int Ind[3];
    
    for(Ind[2] = bdryRange[2][0]; Ind[2]<=bdryRange[2][1]; Ind[2]++){
        for(Ind[1] = bdryRange[1][0]; Ind[1]<=bdryRange[1][1]; Ind[1]++){
            for(Ind[0] = bdryRange[0][0]; Ind[0]<=bdryRange[0][1]; Ind[0]++){
                int r = 0;
        
                int RvInd[] = {Ind[0], Ind[1], Ind[2]};
                for(int axis = 0; axis<bdryParam.fixedAxesNum; axis++){
                    RvInd[fixedAxis[axis]] = 0;
                }
                
                for(int fixedAxisInd = 0; fixedAxisInd<knownAxesNum; fixedAxisInd++){
                    int otherAxis[2]; getOtherAxes(otherAxis, fixedAxis[fixedAxisInd]);
                    for(int nu = 0; nu<NU[fixedAxisInd]; nu++){
                        int stenInd[3] = {0,0,0};
                        for(stenInd[otherAxis[1]] = dimBasedValue(dim, 0, (-p)); stenInd[otherAxis[1]]<= dimBasedValue(dim, 0, p); stenInd[otherAxis[1]]++){
                            for(stenInd[otherAxis[0]] = (-p); stenInd[otherAxis[0]]<= p; stenInd[otherAxis[0]]++){
                                
                                int Rind[] = sumVectors(Ind, stenInd);
                                Rind[fixedAxis[fixedAxisInd]] = 0;
   
                                int fixedFace = fixedSide[fixedAxisInd] + 2*fixedAxis[fixedAxisInd];

                                Rv[ind(RvInd,bdryNgx)][r] = R[ind2(nu,fixedFace,(p+1),6)][ind(Rind,faceBdryPts[fixedAxisInd])];
                                r++;
                            }// end of stenInd loop
                        }// end of stenInd loop
                        
                    }// end of nu loop
                }// end of fixedAxisInd loop
                
                for(int row = r; row<approxEqNum; row++){
                    Rv[ind(RvInd,bdryNgx)][row] = 0;
                }
            }// end of Ind[0] loop
        }// end of Ind[1] loop
    }// end of Ind[2] loop
    
}// end of prepDataVec_corner Function

/// \brief Update the solution values at the ghost points.
/// This is based on the formula (53) in the lcbc.pdf file page 19. Refer to the file to understand the input parameters
/// \param un (input/output): solution values needed at the current time
/// \param mask (input): the mask to determine the type of each grid point
/// \param Rv (input): the vector containing the values related to the data needed from boundary and forcing
/// \param Mat (input): the LcbcMat object containing the LCBC stencil coefficients and other information defined in LCBC.h
void Lcbc::getGhost(double *&un, int *mask, double *Rv[], LcbcMat &Mat, bdryTypeParam &bdryParam, int fixedAxis[], int fixedSide[], int approxEqNum){
    
    /* prepare the parameters */
    int interiorEqNum = bdryParam.interiorEqNum; // number of equations from the interior and boundary
    double b[interiorEqNum]; // a vector that will contain solution values on the interior

    /* Evaluate the solution at the ghost points */
    int i[3];
    for(i[2] = Mat.bdryRange[2][0]; i[2]<=Mat.bdryRange[2][1]; i[2]++){
        for(i[1] = Mat.bdryRange[1][0]; i[1]<= Mat.bdryRange[1][1]; i[1]++){
            for(i[0] = Mat.bdryRange[0][0]; i[0]<=Mat.bdryRange[0][1]; i[0]++){

                /* Fill the solution values on the interior and boundary in the vector b */
                getbInt(b, un, Mat.eqNum, i, Mat.interiorRange, bdryParam.unknownVarNum);
                    
                    int copy[bdryParam.fixedAxesNum];
                    for(int axis = 0; axis<bdryParam.fixedAxesNum; axis++){
                        copy[axis] = i[fixedAxis[axis]];
                        i[fixedAxis[axis]] = 0;
                    }
                
                    int bdInd = ind(i,Mat.bdryNgx); // index used in Rv
                
                    int cVecInd = (cstCoef)?(0):(bdInd); // index used in CaVec and CbVec
                    
                    int vecNum = 0; int g[3];
                    for(g[2] = Mat.ghostRange[2][0]; g[2]<=Mat.ghostRange[2][1]; g[2]++){
                        for(g[1] = Mat.ghostRange[1][0]; g[1]<=Mat.ghostRange[1][1]; g[1]++){
                            for(g[0] = Mat.ghostRange[0][0]; g[0]<=Mat.ghostRange[0][1]; g[0]++){
                                if(CONDITION(g, fixedAxis, fixedSide, bdryParam.knownAxesNum)){
                                    double Ca_vecRv = dotProduct(Mat.CaVec[cVecInd], Rv[bdInd], vecNum, bdryParam.ghostPointNum, approxEqNum);
                                    double Cb_vecb  = dotProduct(Mat.CbVec[cVecInd], b, vecNum, bdryParam.ghostPointNum, interiorEqNum);
                                    
                                    for(int axis = 0; axis<bdryParam.fixedAxesNum; axis++){
                                        i[fixedAxis[axis]] = (copy[axis]+g[fixedAxis[axis]]-p);
                                    }
                                    
                                    if(mask[solInd(i,G.Ngx)]>0)
                                        un[solInd(i,G.Ngx)] = Ca_vecRv + Cb_vecb;  // the definition of solInd is in LCBC_macros.h
                                    
                                    vecNum++;
                                }// end of if statement
                            }// end of g0 loop
                        }// end of g1 loop
                    }// end of g2 loop
                                    
                    for(int axis = 0; axis<bdryParam.fixedAxesNum; axis++){
                        i[fixedAxis[axis]] = copy[axis];
                    }
            }// end of i[0]
        }// end of i[1]
    }// end of i[2]
    
}

/// \brief Find derivatives of data functions (boundary and forcing) needed in the evaluation of CBC equations and their tangential derivatives
/// The function below follows Algorithm 4 in lcbc.pdf (page: 19).
/// \param R (output): a vector carrying  boundary and forcing values from the right-hand-side of each BC and primary CBC on each boundary point
/// \param t (input): current time
/// \param dt (input): time-step
/// \param axis (input)
/// \param side (input)
/// \param gn (input): the LcbcData object carrying boundary data on the grid
/// \param fn (input): the LcbcData object carrying forcing data on the grid
void Lcbc::getDataDeriv(double **&R, double t, double dt, int axis, int side, LcbcData &gn, LcbcData &fn){
    int face = ind2(side,axis,2,3);

    int NU = faceParam[face].NU; // number of primary CBCs (depends on whether the boundary is Dirichlet or Neumann)
    
    int K = (p+1); // number of orders of accuracy needed in the discretization of CBCs
    int R_NU = (p+1); // number of components of the vector R (to correspond to nu CBCs)
    
    int faceType = faceEval[face]; // determine type of BC
    
    /* Allocate space for each component of the vector R */
    for(int nu = 0; nu<NU; nu++){
        R[ind2(nu,face,R_NU,faceNum)] = new double[(faceParam[face].bdryNgx[0]*faceParam[face].bdryNgx[1]*faceParam[face].bdryNgx[2])];
    }
    if(faceType == 2){ // Neumann BC case
        R[ind2(p,face,R_NU,faceNum)] = NULL;
    }
    
    /* Fill the vector R up with values from boundary data */
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
        /* Fill the vector R up with values from forcing data */
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
        
        /* Finalize the construction of R by adding boundary data to forcing data based on right-hand-side of BCs and CBCs */
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

/// \brief Apply Q to a vector F to order 2k accuracy using prepared correction terms
/// \param QF (output): the operator Q applied to the grid function F
/// \param F (input): the grid function F to which we want to apply Q
/// \param W (input): order 2 accurate numerical differences needed as correction terms in the high-order accurate approximation
/// \param bdryRange (input): the range of indices at the boundary
/// \param bdryRangeLth (input): the length in each of the axes at the boundary
/// \param nu (input): represents the power of the Q operator
/// \param k (input): half the order of accuracy needed to discretiza a specific CBC
/// \param newWth (input): center of the new QF grid function
/// \param newLth (input): length of the new QF grid function
/// \param oldWth (input): center of the old grid function F
/// \param oldLth (input): length of the old grid function F
/// \param axis (input)
/// \param side (input)
void Lcbc::applyQhF(double *&QF, double *F, double **W, int bdryRange[3][2], int bdryRangeLth[3], int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side){
    
    int face = side + 2*axis;
    
    /* Get the correction terms and save them in FS */
    double **FS;
    FS = new double*[(maxCoefNum - 1)];
    for(int coefNum = 0; coefNum<(maxCoefNum - 1); coefNum++){
        FS[coefNum] = new double[(oldLth[0]*oldLth[1]*oldLth[2])];
    }// end of coefNum loop
    getCorrectionTerms(FS,F, W, oldLth, k,axis);
    
    int varAxis[2]; getOtherAxes(varAxis, axis);
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
                    
                    QF[I] = QF[I] + coef[face].Fn[coefNum][cI]*mixedDeriv(FS[coefNum], FInd, deriv, G.dx, order, oldLth);
                    
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

/// \brief Apply mapped Neumann operator to a vector F to order 2k accuracy using prepared correction terms
/// \param QF (output): the mapped Neumann operator applied to the grid function F
/// \param F (input): the grid function F to which we want to apply the operator
/// \param W (input): order 2 accurate numerical differences needed as correction terms in the high-order accurate approximation
/// \param bdryRange (input): the range of indices at the boundary
/// \param bdryRangeLth (input): the length in each of the axes at the boundary
/// \param nu (input): represents the power of the Q operator
/// \param k (input): half the order of accuracy needed to discretiza a specific CBC
/// \param newWth (input): center of the new QF grid function
/// \param newLth (input): length of the new QF grid function
/// \param oldWth (input): center of the old grid function F
/// \param oldLth (input): length of the old grid function F
/// \param axis (input)
/// \param side (input)
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

/// \brief find the center and length of the grid functions needed in Algorithm 4 in the lcbc.pdf file based on nu, k and boundary
/// \param wth (output): returns the number of ghost points in each axis of the grid function corresponding to (nu,k) pair
/// \param lth (output): returns the size of the grid function in each direction
/// \param fixedAxis (input): the fixed axis on the given boundary face
/// \param nu (input): parameter invovled in Algorithm 3 (represents powers of the operator Q)
/// \param k (input): parameter invovled in Algorithm 3 (represents the order of accuracy divided by 2)
///  Note that k varies depending on the order of derivatives in the CBCs. Rule of thumb, pick k such that the approximation of the CBCs fits in the stencil
/// \param p (input): orderInSpace/2 (order in space of the full scheme)
/// \param dim (input): dimension
/// \param bdryRangeLth (input): number of boundary grid points in each dimension
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

/// \brief Get the known solution values on the stencil (interior and boundary) and put them in the vector b
/// \param b (output): vector holding known interior values of the solution
/// \param un (input): the given grid solution
/// \param EqNum (input): a vector that holds the numbering of equations used in the LCBC procedure
/// \param Ind (input): the indices of the boundary grid point over which the LCBC procedure is centered
/// \param interiorRange (input): the range of the interior indices
/// \param unknownVarNum (input): the number of unknown variables
void Lcbc::getbInt(double b[], double *un, int *EqNum, int *Ind, int interiorRange[3][2], int unknownVarNum){
    
    int n = (2*p+1);
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
