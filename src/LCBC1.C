#include "LCBC.h"
#include <string.h>
#include <math.h>
#include "utility.h"
#include "LCBCmacros.h"
#include "derivatives.h"
#include "LagrangeDerivFunctions.h"
//#include "timer.h"

/* Note: repeated functions parameters may be documented once when they first appear */

void getWidth(int wth[3], int lth[3], int fixedAxis, int nu, int k, int p, int dim, int faceType);

/// \brief A function to prepare the stencil coefficients used in the LCBC procedure
/// \param Mat (input/output): the LcbcMat object defined in LCBC.h. It carries the LCBC stencil coefficients and other parameters needed for the LCBC evaluation at each face
/// \param bdryParam (input): the type of boundary needed for the LCBC evaluation (face, edge, corner,... see LCBC.h for more information)
/// \param fixedAxis (input): a vector containing the numbers of axes that are fixed
/// \param fixedSide (input): a vector containing the labels of sides that are fixed
/// \param NU (input): a vector containing the numbers of primary CBCs at each fixed axis
/// \param compCondNum (input): a vector containing the number of CBCs and their tangential derivatives for each fixed axis
void Lcbc::prepLcbcStencil(LcbcMat &Mat, bdryTypeParam &bdryParam, int fixedAxis[], int fixedSide[], int NU[], int compCondNum[]){
    
    int knownAxesNum = bdryParam.knownAxesNum;
    int bdryNg = Mat.bdryNgx[0]*Mat.bdryNgx[1]*Mat.bdryNgx[2];
    
    Mat.flag = true;
    Mat.eqNum = new int[totalVarNum];
    Mat.CaVec = new double*[bdryNg];
    Mat.CbVec = new double*[bdryNg];
    
    int approxEqNum = intVectorSum(compCondNum, knownAxesNum) + auxiliaryEqNum;
    
    /* Initialize the objects to carry the LCBC stencil coefficients */
    for(int bdryPoint = 0; bdryPoint<bdryNg; bdryPoint++){
        Mat.CaVec[bdryPoint] = new double[bdryParam.ghostPointNum*approxEqNum];
        Mat.CbVec[bdryPoint] = new double[bdryParam.ghostPointNum*bdryParam.interiorEqNum];
    }// end of bdryPoint loop
    
    /* Assign the number of equations in order in eqNum vectors */
    int n = (2*p + 1);
    int k[3], lth[3] = {n,n,dimBasedValue(dim, 1, n)}, unknownVar = 0, knownVar = 0;
    
            for(k[2] = 0; k[2]<dimBasedValue(dim, 1, n); k[2]++){
                for(k[1] = 0; k[1]<n; k[1]++){
                    for(k[0] = 0; k[0]<n; k[0]++){
                if(CONDITION(k, fixedAxis, fixedSide, knownAxesNum)){
                    Mat.eqNum[ind(k,lth)] = unknownVar;
                    unknownVar++;
                }
                else{
                    Mat.eqNum[ind(k,lth)] = bdryParam.unknownVarNum + knownVar;
                    knownVar++;
                }// end of if statement
            }// end of k[0] loop
        }// end of k[1] loop
    }// end of k[2] loop
    
    getLcbcStencilCoef(Mat.CaVec, Mat.CbVec, Mat.eqNum, Mat.bdryRange, Mat.bdryNgx, Mat.ghostRange, knownAxesNum, bdryParam, fixedAxis, fixedSide, NU, compCondNum, approxEqNum);
}

/// \brief Scale the rows of matrices needed in the LCBC implementation
/// \param A11 (input/output): the upper left block matrix in the LCBC matrix A
/// \param A12 (input/output): the upper right block matrix in the LCBC matrix A
/// \param D_scaled (output): the scaled matrix operator D
/// \param D (input): the matrix operator D
/// \param numRows (input): number of rows of A11 and A12 (the number of equations treated approximately)
/// \param numColms1 (input): number of columns of A11
/// \param numColms2 (input): number of columns of A12
void scaleRows(double *&A11, double *&A12, double *&D_scaled, double *D, int numRows, int numColms1, int numColms2){
    
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

/// \brief Compute the stencil coefficients alpha and beta (page 19 in lcbc.pdf)
/// \param CaVec (output): alpha coefficients (stencil coefficients multiplying values coming from data)
/// \param CbVec (output): beta coefficients (stencil coefficients multiplying values of solution on the interior and at the boundary)
/// \param eqNum (input): a vector containing a numbering of coefficients of the interpolating polynomial
/// \param bdryRange (input): range of indices of grid points on the face
/// \param bdryNgx (input): number of boundary grid points in each dimension
/// \param ghostIndRange (input): range of indices of ghost grid points
/// \param knownAxesNum (input): the number of fixed axes for which the boundary data is known
/// \param bdryParam (input): contains parameters related to the specific boundary type
/// \param fixedAxis (input): a vector containing the axes that are fixed
/// \param fixedSide (input): a vector containing the side labels that are fixed
/// \param NU (input): number of primary CBCs at the faces where the axes are fixed
/// \param compCondNum (input): number of CBCs and their tangential derivatives at the faces where the axes are fixed
/// \param approxEqNum (input): number of equations that are approximated via least-squares in the LCBC procedure
void Lcbc::getLcbcStencilCoef(double **&CaVec, double **&CbVec, int *eqNum, int bdryRange[3][2], int bdryNgx[3], int ghostIndRange[3][2], int knownAxesNum, bdryTypeParam &bdryParam, int fixedAxis[], int fixedSide[], int NU[], int compCondNum[], int approxEqNum){
    
    /* Get the function handle to obtain derivatives of Lagrange polynomials based on the number of fixed axes */
    LagrangeDerivFun LagrangeDeriv[knownAxesNum];
    for(int axisInd = 0; axisInd<knownAxesNum; axisInd++){
        LagrangeDeriv[axisInd] = pickLagrangeDerivFun(dim, p, fixedAxis[axisInd]);
    }
    
    /* Define the block matrices used in the LCBC procedure (see page 15 in lcbc.pdf) */
    double *A11 = new double[(approxEqNum*bdryParam.unknownVarNum)];
    double *A12 = new double[(approxEqNum*bdryParam.interiorEqNum)];
    
    /* Prepare the D matrix operator described in page 20 in lcbc.pdf */
    /* D is the operator multiplying the data vector R(t) */
    double *D = new double[(approxEqNum*approxEqNum)];
    double *D_scaled = new double[(approxEqNum*approxEqNum)];
    
    getD(D, fixedAxis, fixedSide, NU, compCondNum, approxEqNum, knownAxesNum); // get the D operator matrix
    
    int n = (2*p+1);
    int eqNumLth[] = {n,n,dimBasedValue(dim,1,n)};   // The stencil size in each axis
    
    /// Note: Instead of solving the full system of equations in Equation (55) in lcbc.pdf, we pick the specific rows corresponding to the coefficients values needed to evaluate the solution at the ghost points. Extracting the needed rows can be done by forming a matrix E which holds elementary row vectors corresponding to the row number of the needed ghost value.
    /// See LCBCsup.pdf to best understand the procedure. The notation used in this code follows that used in the file LCBCsup.pdf
    double *Et = new double[(bdryParam.unknownVarNum*bdryParam.ghostPointNum)];
    set1DArrayToZero(Et, (bdryParam.unknownVarNum*bdryParam.ghostPointNum));
    double *w = new double[(approxEqNum*bdryParam.ghostPointNum)];
    
    int g[3], vecNum = 0;
    for(g[2] = ghostIndRange[2][0]; g[2]<=ghostIndRange[2][1]; g[2]++){
        for(g[1] = ghostIndRange[1][0]; g[1]<=ghostIndRange[1][1]; g[1]++){
            for(g[0] = ghostIndRange[0][0]; g[0]<=ghostIndRange[0][1]; g[0]++){
                if(CONDITION(g, fixedAxis, fixedSide, knownAxesNum)){
                    int k = eqNum[ind(g,eqNumLth)];
                    Et[ind2(k,vecNum,bdryParam.unknownVarNum,bdryParam.ghostPointNum)] = 1.0;
                    vecNum++;
                }// end of if statement
            }// end of g0 loop
        }// end of g1 loop
    }// end of g2 loop
    
    int Ind[3];
    for(Ind[2] = bdryRange[2][0]; Ind[2]<=bdryRange[2][1]; Ind[2]++){
        for(Ind[1] = bdryRange[1][0]; Ind[1]<=bdryRange[1][1]; Ind[1]++){
            for(Ind[0] = bdryRange[0][0]; Ind[0]<=bdryRange[0][1]; Ind[0]++){
                
                getBlockMatrices(A11, A12, Ind, eqNum, bdryParam, NU, compCondNum, approxEqNum, LagrangeDeriv, fixedSide, fixedAxis); // get the block matrices A11 and A12 needed in the LCBC implementation
       
                /* scale the matrices A11, A12 and D */
                scaleRows(A11, A12, D_scaled, D, approxEqNum, bdryParam.unknownVarNum, bdryParam.interiorEqNum);

                /* get the index for the vectors that contain the stencil coefficients */
                int copy[bdryParam.fixedAxesNum];
                for(int axis = 0; axis<bdryParam.fixedAxesNum; axis++){
                    copy[axis] = Ind[fixedAxis[axis]];
                    Ind[fixedAxis[axis]] = 0;
                }
                
                int cVecInd = (cstCoef)?(0):(ind(Ind,bdryNgx)); // index of CaVec and CbVec

                /* do a LS solve to obtain Ca and Cb */
                LSsolve(A11, Et, w, approxEqNum, bdryParam.unknownVarNum, bdryParam.ghostPointNum, 'T');
                // solve the transpose system (see LCBCsup.pdf for more information)
                
                transposeMultiplyMatrices(CbVec[cVecInd], w, A12, approxEqNum, bdryParam.ghostPointNum, approxEqNum, bdryParam.interiorEqNum, (-1));
                transposeMultiplyMatrices(CaVec[cVecInd], w, D_scaled, approxEqNum, bdryParam.ghostPointNum, approxEqNum, approxEqNum, (1));
                
                /* return to the original index */
                for(int axis = 0; axis<bdryParam.fixedAxesNum; axis++){
                    Ind[fixedAxis[axis]] = copy[axis];
                }
                
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
}

/// \brief Get the block matrices needed in the LCBC procedure
/// In this function derivatives of Lagrange polynomials need to be computed. See the corresponding functions used in LagrangeDerivFunctions.C/h
/// \param A11 (output): the upper left block matrix defined in page 20 in the lcbc.pdf file
/// \param A12 (output): the upper right block matrix defined in page 20 in the lcbc.pdf file
/// \param Ind (input): the index of the boundary grid point over which the LCBC procedure is centered
/// \param eqNum (input): a vector containing a numbering of coefficients of the interpolating polynomial
/// \param bdryType (input): the type of boundary needed for the LCBC evaluation (face, edge, corner,... see LCBC.h and LCBC_param.h for more information)
/// \param approxEqNum (input): number of equations approximated through least-squares
/// \param compCondNum (input): number of equations coming in from compatiblity boundary conditions and their tangentail derivatives
/// \param fixedAxis (input): vector containing the fixed axes at the specific boundary object (face, edge or corner)
/// \param fixedSide (input): vector containing the fixed sides at the specific boundary object (face, edge or corner)
/// \param NU (input): the number of primary CBCs on the specific face
/// \param LagrangeDeriv (input): the function we will use to find the derivative of Lagrange polynomial (see list in LagrangeDerivFunctions.C/h)
/// See the defined LagrangeDerivFun in LCBCmacros.h
void Lcbc::getBlockMatrices(double *&A11, double *&A12, int *Ind, int *eqNum, bdryTypeParam &bdryType, int NU[3], int compCondNum[3], int approxEqNum, LagrangeDerivFun LagrangeDeriv[], int fixedSide[], int fixedAxis[]){
    
    int unknownVarNum = bdryType.unknownVarNum;
    int knownAxesNum = bdryType.knownAxesNum;

    int n = (2*p+1), MU = n;
    int totalVarNum = n*n*dimBasedValue(dim, 1, n);
    int NUmax = (p+1);
    
    int fixedFace[3]; // faces corresponding to each of the fixed axes
    
    /* Prepare Variables that will hold the Lagrange derivative information */
    double *Z[knownAxesNum*totalVarNum];
    for(int axis = 0; axis<knownAxesNum; axis++){
        fixedFace[axis] = fixedSide[axis] + 2*fixedAxis[axis];
    }

    int LagInd[3], LagIndLth[3] = {n,n,dimBasedValue(dim, 1, n)};
    
    for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim, 1, n); LagInd[2]++){
        for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
            for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                
                for(int axis = 0; axis<knownAxesNum; axis++){
                    
                    int zInd = axis + knownAxesNum*ind(LagInd,LagIndLth);
                    Z[zInd] = new double[compCondNum[axis]];
                
                    LagrangeDeriv[axis](Z[zInd],Ind, LagInd, LagrangeData, LagrangeData_center, coef[fixedFace[axis]].Fn, coef[fixedFace[axis]].wth, coef[fixedFace[axis]].lth, cstCoef, G.dx, fixedAxis[axis], fixedSide[axis], NU[axis], memory);
                }
            }// LagInd[0]
        }// LagInd[1]
    }// LagInd[2]
    
    int row = 0;
    for(int nu = 0; nu<NUmax; nu++){
        for(int mu1 = 0; mu1<dimBasedValue(dim,1,MU); mu1++){
            for(int mu0 = 0; mu0<MU; mu0++){
                
                for(int axis = 0; axis<knownAxesNum; axis++){
                    if(nu<NU[axis]){
                        
                        for(LagInd[2] = 0; LagInd[2]<dimBasedValue(dim,1,n); LagInd[2]++){
                            for(LagInd[1] = 0; LagInd[1]<n; LagInd[1]++){
                                for(LagInd[0] = 0; LagInd[0]<n; LagInd[0]++){
                                    int colm = eqNum[ind(LagInd,LagIndLth)];
                                    
                                    int zInd = axis + knownAxesNum*ind(LagInd,LagIndLth);
                                    
                                    if(colm<unknownVarNum){
                                        A11[ind2(row,colm,approxEqNum,unknownVarNum)] = Z[zInd][ind3(nu,mu0,mu1,NU[axis],MU,MU)];
                                    }
                                    else{
                                        int A12colm = colm - unknownVarNum;
                                        A12[ind2(row,A12colm,approxEqNum,interiorEqNum)] = Z[zInd][ind3(nu,mu0,mu1,NU[axis],MU,MU)];
                                    }
                                }// LagInd[0]
                            }// LagInd[1]
                        }// LagInd[2]
                        
                        row++;
                        
                    }// end of if nu<NU[0]
                }

            }// mu0 loop
        }// mu1 loop
    }// end of nu loop
    
    for(int k = 0; k<(knownAxesNum*totalVarNum); k++){
        delete [] Z[k];
    }
}// end of fillVertexMatrix_LagrangeDeriv

/// \brief Get the matrix operator D representing discrete derivatives of the data vector R(t).
/// See page 20 of the lcbc.pdf document for a definition of D
/// \param D (output): the matrix operator D
/// \param fixedAxis (input)
/// \param fixedSide (input)
/// \param NU (input)
/// \param compCondNum (input)
/// \param approxEqNum (input)
/// \param knownAxesNum (input): number of fixed axes for which the boundary data is given
void Lcbc::getD(double *&D, int fixedAxis[], int fixedSide[], int NU[], int compCondNum[], int approxEqNum, int knownAxesNum){
    
    int n = (2*p+1);
    int dof = (n*dimBasedValue(dim, 1, n));
    int totalCompCondNum = intVectorSum(compCondNum, knownAxesNum);
    int totalNuNum = intVectorSum(NU, knownAxesNum);;
    
    /* use the delta approach to determine D */
    double *Id = new double[(dof*dof)];
    getIDmatrix(Id, dof);
    
    /* allocate space for two vectors */
    // Vector R: carries derivatives of data functions from the RHS of primary CBCs. See equation (44) on page 15 of the lcbc.pdf ArXiv document
    // Vector b: the right-hand-side vector of the LCBC system in equation (46) in lcbc.pdf corresponding to the right-hand-side coming from the CBCs and their tangential derivatives only (not interior)
    double b[totalCompCondNum], **R = new double*[totalNuNum];
    
    /* start the delta approach procedure */
    for(int nu = 0; nu<(totalNuNum); nu++)
        R[nu] = new double[dof];
    
    int colm = 0;
    for(int fixedAxisInd = 0; fixedAxisInd<knownAxesNum; fixedAxisInd++){
        
        int Rlth[] = {n,n,dimBasedValue(dim, 1, n)}; // this is a vector that carries the length of R in each axis
        Rlth[fixedAxis[fixedAxisInd]] = 1;
        
        int otherAxis[2]; getOtherAxes(otherAxis, fixedAxis[fixedAxisInd]);
        
        for(int nu = 0; nu<NU[fixedAxisInd]; nu++){
            set2DArrayToZero(R, (totalNuNum), dof);
            for(int l = 0; l<dof; l++){
                int row = 0;
                
                int stenInd[] = {0,0,0};
                for(stenInd[otherAxis[1]] = 0; stenInd[otherAxis[1]]<dimBasedValue(dim, 1, n); stenInd[otherAxis[1]]++){
                    for(stenInd[otherAxis[0]] = 0; stenInd[otherAxis[0]]<n; stenInd[otherAxis[0]]++){
        
                        R[ind2(nu,fixedAxisInd,NU[0],knownAxesNum)][ind(stenInd,Rlth)] = Id[ind2(row,l,dof,dof)];
                        
                        row++;
                    }// end of i loop
                }// end of j loop
                
                getbVec(b, R, NU, totalCompCondNum, knownAxesNum, fixedAxis); // get the b from the values of R
                
                for(int k = 0; k<(totalCompCondNum); k++){
                    D[ind2(k,colm,approxEqNum,approxEqNum)] = b[k];
                }// end of k loop
                for(int k = totalCompCondNum; k<approxEqNum; k++){
                    D[ind2(k,colm,approxEqNum,approxEqNum)] = 0.0;
                }
                colm++;
                
            }// end of l loop
        }// end of nu loop
        
        for(int c = colm; c<approxEqNum; c++){
            for(int r = 0; r<approxEqNum; r++){
                D[ind2(r,c,approxEqNum,approxEqNum)] = 0.0;
            }
        }
    }// end of fixedAxisInd loop
    
    /* delete variables */
    for(int nu = 0; nu<totalNuNum; nu++){
        delete [] R[nu];
        R[nu] = NULL;
    }// end of nu loop
    
    delete [] R;
    delete [] Id;
}// end of getD_vertex

/// \brief Get the right-hand-side of the LCBC matrix system corresponding to equations involving forcing and boundary functions.
/// This function is used in the delta approach procedure carried out in getD function above and written for that purpose.
/// \param b (output): The right-hand-side vector described above (see lcbc.pdf for more info)
/// \param R (input): a vector carrying data infomation from the right-hand-side of primary CBCs
void Lcbc::getbVec(double b[], double **R, int NU[3], int totalCompCondNum, int knownAxesNum, int fixedAxis[]){
    int n = (2*p+1), order[knownAxesNum], ordVec[knownAxesNum][3];
    int MU0 = n, MU1 = dimBasedValue(dim, 1, n);
    memset(b, 0, (totalCompCondNum*sizeof(double)));

    int Ind[] = {p,dimBasedValue(dim,0,p),0};
    int sizeRnu[] = {n,dimBasedValue(dim, 1, n),1};
    
    double hx[knownAxesNum][3];
    for(int axis = 0; axis<knownAxesNum; axis++){
        int otherAxis[2]; getOtherAxes(otherAxis, fixedAxis[axis]);
        hx[axis][0] = G.dx[otherAxis[0]];
        hx[axis][1] = G.dx[otherAxis[1]];
        hx[axis][2] = 0;
    }
    
    int eqnNum = 0, mu[3] = {0,0,0};
    for(int nu = 0; nu<(p+1); nu++){
        
        for(int axis = 0; axis<knownAxesNum; axis++){
            order[axis] = (nu==0)?(2*p):(2*(NU[axis] - nu));
            
            ordVec[axis][0] = order[axis];
            ordVec[axis][1] = dimBasedValue(dim,0,order[axis]);
            ordVec[axis][2] = 0;
        }
        
        for(mu[1] = 0; mu[1]<MU1; mu[1]++){
            for(mu[0] = 0; mu[0]<MU0; mu[0]++){
                for(int axis = 0; axis<knownAxesNum; axis++){
                    if(nu<NU[axis]){
                        
                        int addon = intVectorSum(NU, axis);
                        
                        b[eqnNum] = mixedDeriv(R[(nu + addon)], Ind, mu, hx[axis], ordVec[axis], sizeRnu);
                        eqnNum++;
                    }
                }
                
                if( (mu[0]!=0) && (mu[0]%2)==0 ){
                    for(int axis = 0; axis<knownAxesNum; axis++){
                        if(ordVec[axis][0]>2){ordVec[axis][0] = ordVec[axis][0] - 2;}
                    }
                }// end of if mu0
            }// end of mu0 loop
            
            for(int axis = 0; axis<knownAxesNum; axis++){
                ordVec[axis][0] = order[axis];
            }
            
            if( (mu[1]!=0) && (mu[1]%2)==0 ){
                for(int axis = 0; axis<knownAxesNum; axis++){
                    if(ordVec[axis][1]>2){ordVec[axis][1] = ordVec[axis][1] - 2;}
                }
            }// end of if mu1
        }// end of mu1 loop
    }// end of nu loop
}// end of getbVec_vertex

/// \brief Get the Lagrange derivatives corresponding to the CBC equation in the LCBC procedure
/// The function follows Algorithm 3 in the LCBC paper lcbc.pdf and the notation used in it
/// \param Y (output): vector carrying the values corresponding to the auxiliary equations
/// \param Z (output): vector carrying the values corresponding to the CBC equations
/// \param Ind (input): the index of the boundary point over which the LCBC procedure is centered
/// \param LInd (input): the index of the Lagrange polynomial functions
/// \param axis (input)
/// \param side (input)
/// \param evaluateY (input): a boolean to tell whether to fill the values of Y or leave it out
void Lcbc::getLagrangeDeriv(double Y[], double Z[], int *Ind, int *LInd, int axis, int side, bool evaluateY){
    /* This function finds and saves derivatives of the Lagrange polynomial centered at a boundary point [Ind] */
    
    int face = side + 2*axis;
    int NU = faceParam[face].NU; // number of primary CBCs at a given face
    int NUm = p+1; // a fixed value equal to p + 1 (does not depend on type of BC)
    int faceType = faceEval[face]; // type of boundary condition on the face
    
    int K = (p+1), MU = 2*p+1;
    /* prepare a vector containing the variable axes */
    int varAxis[2] = {faceParam[face].otherAxis[0], faceParam[face].otherAxis[1]};
    
    /* prepare a vector V and W (see Algorithm 3 in lcbc.pdf) */
    // V: contains Q^{nu}_{h,2k}(Lagrange grid function)
    // W: contains differences of V used in correction terms of high-order accurate approximations
    double *V[(NUm*K)]; // V[nu,k][Grid_index]
    double *W[(K*p*p*p)];
    
    /* prepare length (lth) and width (wth) vectors */
    // lth: the length of the grid functions in V
    // wth: the center of the grid functions in V
    int wth[3], lth[3];
    int nu = 0; int ll = LagrangeData_center, Lw = (2*ll + 1); // center and length of Lagrange data given
    
    /* Fill values of V when nu = 0 with Lagrange polynomial data */
    for(int k = 1; k<=p; k++){
        getWidth(wth, lth, axis, nu, k, p, dim, faceType); // get the lth and width for a given nu and k
        
        V[ind2(nu,k,NUm,K)] = new double[(lth[0]*lth[1]*lth[2])];
        int i[3];
        for(i[2] = -wth[2]; i[2]<=wth[2]; i[2]++){
            for(i[1] = -wth[1]; i[1]<=wth[1]; i[1]++){
                for(i[0] = -wth[0]; i[0]<=wth[0]; i[0]++){
                    double Ln = 1;
                    for(int d = 0; d<dim; d++){
                        Ln = Ln*LagrangeData[ind2(LInd[d],(ll+i[d]),(2*p+1),Lw)];
                    }
                    int V_ind[] = sumVectors(wth, i);
                    V[ind2(nu,k,NUm,K)][ind(V_ind,lth)] = Ln;
                }// end of i0
            }// end of i1
        }// end of i2
    }// end of k
    
    int ord[] = {2,2,dimBasedValue(dim, 0, 2)};
    
    int V_wth[3], V_lth[3];
    for(int nu = 0; nu<=(p-1); nu++){
        for(int k = 1; k<=(p-nu); k++){
            if(k != 1){
                getWidth(wth, lth, axis, nu, k, p, dim, faceType);
                for(int l = 1; l<=(k-1); l++){
                    getWidth(V_wth, V_lth, axis, nu, l, p, dim, faceType);
                    
                    int m = (k-l);
                    for(int Dx = 0; Dx<=m; Dx++){
                        for(int Dy = dimBasedValue(dim, (m-Dx), 0); Dy<=(m-Dx); Dy++){
                            int Dz = dimBasedValue(dim, 0, (m - Dx - Dy));
                            W[ind4(l,Dx,Dy,Dz,K,p,p,p)] = new double[(lth[0]*lth[1]*lth[2])];
                            
                            int deriv[] = {(2*Dx),(2*Dy),(2*Dz)};
                            int i[3];
                            for(i[2] = (-wth[2]); i[2]<=(wth[2]); i[2]++){
                                for(i[1] = (-wth[1]); i[1]<=(wth[1]); i[1]++){
                                    for(i[0] = (-wth[0]); i[0]<=(wth[0]); i[0]++){
                                        int V_ind[] = sumVectors(V_wth, i);
                                        int W_ind[] = sumVectors(wth, i);
                                        
                                        W[ind4(l,Dx,Dy,Dz,K,p,p,p)][ind(W_ind,lth)] = mixedDeriv(V[ind2(nu,l,NUm,K)], V_ind, deriv, G.dx, ord, V_lth);
                                        
                                    }// end of i0 loop
                                }// end of i1 loop
                            }// end of i2 loop
                            
                        }// end of Dy
                    }// end of Dx
                }// end of l loop
            }// end of if k statement
            
            int oldWth[3], oldLth[3];
            getWidth(oldWth, oldLth, axis, nu, k, p, dim, faceType);
            
            getWidth(wth, lth, axis, (nu+1), k, p, dim, faceType);
            V[ind2((nu+1),k,NUm,K)] = new double[(lth[0]*lth[1]*lth[2])];
            
            if(k<(p-nu) || faceType == 1){
                applyQh(V[ind2((nu+1),k,NUm,K)], V[ind2(nu,k,NUm,K)], W, Ind, nu, k, wth, lth, oldWth, oldLth, axis, side);
            }else{
                applyDh(V[ind2((nu+1),k,NUm,K)], V[ind2(nu,k,NUm,K)], W, Ind, nu, k, wth, lth, oldWth, oldLth, axis, side);
            }
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
    int MU0 = MU, MU1 = dimBasedValue(dim, 1, MU), k, mu[2];
    
    /* Create an array that will hold the values of du/dn or u */
    int V1_lth[] = {(2*p+1),(2*p+1),dimBasedValue(dim, 1, (2*p+1))};
    int V1_wth[] = {p,p,dimBasedValue(dim, 0, p)};
    V1_lth[axis] = 1; V1_wth[axis] = 0;
    double *V1 = new double[(V1_lth[0]*V1_lth[1]*V1_lth[2])];
    
    int deriv[] = {0,0,0}, order[3] = {0,0,0};
    for(int nu = 0; nu<NU; nu++){
        
        k = (nu == 0)?(p):(NU - nu);
        
        if(faceType == 1){
            getWidth(wth, lth, axis, nu, k, p, dim, faceType);
        }else{
            getWidth(wth, lth, axis, (nu+1), k, p, dim, faceType);
        }
        
        int i[3];
        for(i[2] = 0; i[2]<V1_lth[2]; i[2]++){
            for(i[1] = 0; i[1]<V1_lth[1]; i[1]++){
                for(i[0] = 0; i[0]<V1_lth[0]; i[0]++){
                    int vInd[3];
                    vInd[axis] = wth[axis];
                    vInd[varAxis[0]] = wth[varAxis[0]] - p + i[varAxis[0]];
                    vInd[varAxis[1]] = dimBasedValue(dim, 0, (wth[varAxis[1]] - p + i[varAxis[1]]));
                    if(faceEval[face] == 2){
                        V1[ind(i,V1_lth)] = V[ind2((nu+1),k,NUm,K)][ind(vInd,lth)];
                    }// end of Neumann condition
                    else if(faceEval[face] == 1){
                        V1[ind(i,V1_lth)] = V[ind2(nu,k,NUm,K)][ind(vInd,lth)];
                    }// end of Dirichlet condition
                }// end of i[0]
            }// end of i[1]
        }// end of i[2]
        
        /* prepare the tangential derivatives here */
        int kappa[2] = {k,k};
        for(mu[1] = 0; mu[1]<MU1; mu[1]++ ){
            for(mu[0] = 0; mu[0]<MU0; mu[0]++){
                
                for(int l = 0; l<2; l++){
                    deriv[varAxis[l]] = mu[l];
                    order[varAxis[l]] = 2*kappa[l];
                }// end of l loop

                Z[ind3(nu,mu[0],mu[1],NU,MU0,MU1)] = mixedDeriv(V1, V1_wth, deriv, G.dx, order, V1_lth);
     
                if( (mu[0]!=0) && (mu[0]%2)==0 && kappa[0]>1 ){
                    kappa[0] = kappa[0] - 1;
                }// end of if mu0
            }// end of mu0 loop
            kappa[0] = k;
            if( (mu[1]!=0) && (mu[1]%2)==0 && kappa[1]>1 ){
                kappa[1] = kappa[1] - 1;
            }// end of if mu1
        }// end of mu1
    }// end of nu loop

    delete [] V1;
    
    /* prepare the values related to auxiliary equations here (if used) */
    /* Calculate \partial_x^{m1}\partial_y^{m2}L_iL_j evaluated at boundary point when (m1+m2)>(2p) */
    if(evaluateY){
        int cnt = 0, m[3], M_lth[3] ={MU,MU,dimBasedValue(dim, 1, MU)};
        int order[3];
        getWidth(wth, lth, axis, 0, 1, p, dim, faceType);
        for(m[2] = 0; m[2]<M_lth[2]; m[2]++){
            for(m[1] = 0; m[1]<M_lth[1]; m[1]++){
                for(m[0] = 0; m[0]<M_lth[0]; m[0]++){
                    if((m[0] + m[1] + m[2])>(2*p)){
                        order[0] = (m[0]!=0)?(2*getOrder(m[0],p)):(0);
                        order[1] = (m[1]!=0)?(2*getOrder(m[1],p)):(0);
                        order[2] = (m[2]!=0)?(2*getOrder(m[2],p)):(0);
                        Y[cnt] =mixedDeriv(V[ind2(0,1,NUm,K)], wth, m, G.dx, order, lth);
                        cnt = cnt + 1;
                    }// end if statement
                }// end m[0] loop
            }// end m[1] loop
        }// end m[2] loop
    }// end if evaluateY
    
    /* free all the pointers allocated in this function */
    for(int k = 1; k<=p; k++){
        delete [] V[ind2(0,k,NUm,K)];
        V[ind2(0,k,NUm,K)] = NULL;
    }
    for(int nu = 0; nu<=(p-1); nu++){
        for(int k = 1; k<=(p - nu); k++){
            delete [] V[ind2((nu+1),k,NUm,K)];
            V[ind2((nu+1),k,NUm,K)] = NULL;
        }// end of k loop
    }// end of nu loop
    /* end of pointer freeing */
}// end of Dirichlet_LagrangeDer3 function

/// \brief function to apply Q to a grid function V and to order 2k
/// \param QV (output): Q applied to the grid function V
/// \param V (input): given grid function V
/// \param W (input): correction terms needed in the evaluation of QV to the desired order of accuracy
/// \param Ind (input): index of the boundary point at which the LCBC procedure is centered
/// \param nu (input): the parameter utilized in the Algorithm 3 in lcbc.pdf. It represents the power of the operator Q
/// \param k (input): order of accracy needed / 2
/// \param newWth (input): center of the new QV grid function
/// \param newLth (input): length of the new QV grid function
/// \param oldWth (input): center of the old grid function V
/// \param oldLth (input): length of the old grid function V
/// \param axis (input)
/// \param side (input)
void Lcbc::applyQh(double *&QV, double *V, double **W, int *Ind, int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side){
    
    int face = side + 2*axis;
    
    /* Get the correction terms and save them in VS */
    double **VS;
    VS = new double*[(maxCoefNum-1)];
    for(int coefNum = 0; coefNum<(maxCoefNum-1); coefNum++){
        VS[coefNum] = new double[(oldLth[0]*oldLth[1]*oldLth[2])];
    }
    getCorrectionTerms(VS,V, W, oldLth, k,axis);
    
    int v0 = faceParam[face].otherAxis[0];
    int v1 = faceParam[face].otherAxis[1];
    int cI = 0;
    
    /* Find Qh of V */
    int i[3], order[] = {2,2,2}; int cInd[3];
    for(i[2] = (-newWth[2]); i[2]<=newWth[2]; i[2]++){
        for(i[1] = (-newWth[1]); i[1]<=newWth[1]; i[1]++){
            for(i[0] = (-newWth[0]); i[0]<=newWth[0]; i[0]++){
                
                if(!cstCoef){
                    cInd[axis] = coef[face].wth[axis] + i[axis];
                    cInd[v0] = Ind[v0] -p + coef[face].wth[v0] + i[v0];
                    cInd[v1] = dimBasedValue(dim, 0, (Ind[v1] - p + coef[face].wth[v1] + i[v1]));
                    cI = ind(cInd, coef[face].lth);
                }
                
                int vInd[] = sumVectors(oldWth, i);
                int qInd[] = sumVectors(newWth,i);

                int qI = ind(qInd,newLth);
                
                /* Parts of Q corresponding to derivatives with respect to one variable */
                QV[qI] = 0; int coefNum = 0;
                for(int degree = 2; degree>0; degree = degree - 1){
                    for(int d = 0; d<dim; d++){
                        
                        QV[qI] = QV[qI] + coef[face].Fn[coefNum][cI]*Deriv(VS[coefNum], vInd, d, degree, G.dx[d], 2, oldLth);
                        
                        coefNum++;
                    }
                }
                
                /* parts of Q corresponding to the mixed derivative terms */
                for(int d = (2*(dim-2)); d>=0; d--){
                    int deriv[] = {1,1,1};
                    (dim>2) ? (deriv[d] = 0) : (deriv[2] = 0);
                
                    QV[qI] = QV[qI] + coef[face].Fn[coefNum][cI]*mixedDeriv(VS[coefNum], vInd, deriv, G.dx, order, oldLth);

                    coefNum++;
                }
                QV[qI] = QV[qI] + coef[face].Fn[coefNum][cI]*V[ind(vInd,oldLth)];
                
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

/// \brief Apply normal derivative to high-order accuracy for Neumann BC type faces
/// \param QV (output): the result of the normal derivative applied to the grid function V
/// \param V (input): the input grid function V for which we want to apply the normal derivative
/// \param W (input): correction terms needed in the evaluation of QV to the desired order of accuracy
/// \param Ind (input): index of the boundary point at which the LCBC procedure is centered
/// \param nu (input): the parameter utilized in the Algorithm 3 in lcbc.pdf. It represents the power of the operator Q
/// \param k (input): order of accracy needed / 2
/// \param newWth (input): center of the new QV grid function
/// \param newLth (input): length of the new QV grid function
/// \param oldWth (input): center of the old grid function V
/// \param oldLth (input): length of the old grid function V
/// \param axis (input)
/// \param side (input)
void Lcbc::applyDh(double *&QV, double *V, double **W, int *Ind, int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side){
    
    int face = side + 2*axis;
    /* Get the correction terms and save them in VS */
    double **VS;
    VS = new double*[dim];
    for(int axisNum = 0; axisNum<dim; axisNum++){
        VS[axisNum] = new double[(oldLth[0]*oldLth[1]*oldLth[2])];
    }
    getCorrectionTermsForDh(VS,V, W, oldLth, k,axis);

    /* Find Qh of V */
    int i[3];
    for(i[2] = (-newWth[2]); i[2]<=newWth[2]; i[2]++){
        for(i[1] = (-newWth[1]); i[1]<=newWth[1]; i[1]++){
            for(i[0] = (-newWth[0]); i[0]<=newWth[0]; i[0]++){
                
                int vInd[] = sumVectors(oldWth,i);
                int qInd[] = sumVectors(newWth,i);
                int qI = ind(qInd,newLth);
                
                double arg[] = {(double)face, (G.x[0][Ind[0]] + i[0]*G.dx[0]), (G.x[1][Ind[1]] + i[1]*G.dx[1]), (G.x[2][Ind[2]] + i[2]*G.dx[2])};
                
                double bn_vec[3];
                bn(bn_vec, arg); // get the mapped Neumann boundary operator
                
                QV[qI] =  bn_vec[0]*Deriv(VS[0], vInd, 0, 1, G.dx[0], 2, oldLth)
                        + bn_vec[1]*Deriv(VS[1], vInd, 1, 1, G.dx[1], 2, oldLth);
                
                if(dim == 3){
                    QV[qI] = QV[qI] + bn_vec[2]*Deriv(VS[2], vInd, 2, 1, G.dx[2], 2, oldLth);
                }// end of if dim=3
                
            }// end of i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Free the pointers */
    for(int i = 0; i<dim; i++){
        delete [] VS[i];
        VS[i] = NULL;
    }
    delete [] VS;
}

/// \brief get the correction terms needed for the 2kth-order accurate derivative approximations in Q^nu
/// \param VS (output): the sum of correction terms needed in the approximation
/// \param V (input): the grid function to which we are applying Q to
/// \param W (input): objects holding order 2 differences needed in the correction terms
/// \param lth (input): the length of the grid functions VS, V and W
/// \param k (input): the order of accuracy divided by 2
/// \param axis (input): the fixed axis of the boundary
void Lcbc::getCorrectionTerms(double **&VS, double *V, double **W, int *lth, int k, int axis){
    /* Note that 0,1,2,3,4 correspond to c11, c22, c1, c2, c12 correction terms */
#define coefVec(num) ((num < dim) ? (a) : (b))
    
    /* Prepare coefficients needed in the correction terms */
    double a[] = {0,(-1.0/12.0),(1.0/90.0), (-1.0/560.0), (1.0/3150.0)};
    double b[] = {(1.0),(-1.0/6.0),(1.0/30.0), (-1.0/140.0), (1.0/630.0)};
    
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

/// \brief get the correction terms needed for the 2kth-order accurate derivative approximations in the normal derivative of V
/// \param VS (output): the sum of correction terms needed in the approximation
/// \param V (input): the grid function to which we are applying Q to
/// \param W (input): objects holding order 2 differences needed in the correction terms
/// \param lth (input): the length of the grid functions VS, V and W
/// \param k (input): the order of accuracy divided by 2
/// \param axis (input): the fixed axis of the boundary
void Lcbc::getCorrectionTermsForDh(double **&VS, double *V, double **W, int *lth, int k, int axis){
    /* Note that 0,1,2,3,4 correspond to c11, c22, c1, c2, c12 correction terms */
    
    /* prepare the coefficients of the correction terms in the vector b */
    double b[] = {(1.0),(-1.0/6.0),(1.0/30.0)};
    
    int i[3];
    
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<lth[0]; i[0]++){
                
                for(int var = 0; var<dim; var++){
                    int m[] = {0,0,0};
                    VS[var][ind(i,lth)] = 0;
                    for(int n = 1; n<=(k-1); n++){
                        m[var] = n;
                        VS[var][ind(i,lth)] = VS[var][ind(i,lth)] + b[n]*pow(G.dx[var],((double) (2.0*n)))*W[ind4((k-n),m[0],m[1],m[2],(p+1),p,p,p)][ind(i,lth)];
                    }// end of n loop
                    VS[var][ind(i,lth)] = VS[var][ind(i,lth)] + V[ind(i,lth)];
                }// end of var loop
                
            }// end of i[0] loop
        }// end of i[1] loop
    } // end of i[2] loop
}// end of getCorrectionTermsForDh

/// \brief find the center and length of the grid functions needed in Algorithm 3 in the lcbc.pdf file based on nu, k and boundary
/// \param wth (output): returns the indices of the center point of the grid function corresponding to (nu,k) pair
/// \param lth (output): returns the size of the grid function in each direction
/// \param fixedAxis (input): the fixed axis on the given boundary face
/// \param nu (input): parameter invovled in Algorithm 3 (represents powers of the operator Q)
/// \param k (input): parameter invovled in Algorithm 3 (represents the order of accuracy divided by 2)
///  Note that k varies depending on the order of derivatives in the CBCs. Rule of thumb, pick k such that the approximation of the CBCs fits in the stencil
/// \param p (input): orderInSpace/2 (order in space of the full scheme)
/// \param dim (input): dimension
void getWidth(int wth[3], int lth[3], int fixedAxis, int nu, int k, int p, int dim, int faceType){
    
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
}// end of getWidth
