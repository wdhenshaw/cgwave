#ifndef LCBC_h
#define LCBC_h

#include <stdio.h>
#include <stdlib.h>
#include "numericalDeriv.h"
#include "grid.h"
#include "LCBCmacros.h"

typedef double (*F1)(double *);

struct LcbcMat{
    bool flag = false;
    int *eqNum; 
    double ***CaVec;
    double ***CbVec;
};

class Lcbc{
public:
    /* Main Function to update ghost */
    void updateGhost(double *unp1, double t, double dt);
    
    /* LCBC 3D Vertex */
    void updateVertexGhost(double *R[], double *unp1, double t, double dt, int approxEqNum, int side0, int side1, int side2);
    void prepVertexMatrix(int side0, int side1, int side2);
    void getVertexMatrix(double **CaVec, double **CbVec, int *eqNum, int side0, int side1, int side2);
    void fillVertexMatrix_LagrangeDeriv(double *Matrix, int *Ind, int side0, int side1, int side2, int *eqNum, int compCondNum, int auxiliaryEqNum, int totalEqNum);
    void getD_vertex(double *D);
    void getbVec_vertex(double *b, double **R);
    void prepDataVec_vertex(double *R[], double *Rv, double t, double dt, int approxEqNum, int side0, int side1, int side2);
    void getVertexGhost(double *un, double *Rv, double **CaVec, double **CbVec, int *eqNum, int approxEqNum, int side0, int side1, int side2);
    void getbInt_vertex(double *b, double *un, int *eqNum, int *Ind, int side0, int side1, int side2);
    void freeVertexVariables();
    /* --------------------------------------*/

    /* LCBC 2D Corner or 3D edge */
    void updateEdgeGhost(double *R[], double *unp1, double t, double dt, int approxEqNum, int side1, int side2, int varAxis);
    void prepCornerMatrix(int side1, int side2, int varAxis);
    void getCornerMatrix(double ***CaVec, double ***CbVec, int *eqNum, int bdryRange[3][2], int bdryNg, int side1, int side2, int varAxis, int fixedAxis[2]);
    void getD_corner(double *D, int varAxis, int fixedAxis[2]);
    void getbVec_corner(double *b, double **R, int varAxis, int fixedAxis[2]);
    void fillCornerMatrix_LagrangeDeriv(double Matrix[], int *Ind, int side1, int side2, int varAxis, int fixedAxis[2], int *eqNum, int compCondNum, int auxiliaryEqNum, int totalEqNum);
    void prepDataVec_corner(double *R[], double *Rv[], double t, double dt, int side1, int side2, int varAxis, int approxEqNum);
    void getCornerGhost(double *un, double *Rv[], double ***CaVec, double ***CbVec, int *eqNum, int side1, int side2, int varAxis, int approxEqNum);
    void getbInt_corner(double *b, double *un, int *eqNum, int *Ind, int side1, int side2, int fixedAxis[2]);
    void getBdryRange_corner(int bdryRange[3][2], int fixedAxis[2], int fixedSide[2], int varAxis, int addOnSide0, int addOnSide1);
    void freeCornerVariables();
    /* ------------------------------------- */

    /* LCBC 2D or 3D face  */
    /* part 1: matrix preparation */
    void updateFaceGhost(double *R[], double *unp1, double t, double dt, int approxEqNum, int axis, int side);
    void prepSideMatrix(int axis, int side);
    void getSideMatrix(double ***CaVec, double ***CbVec, int *eqNum, int bdryRange[3][2], int bdryNgx[3], int axis, int side);
    void getD(double *D, int axis);
    void getbVec(double *bd, double **R, int axis);
    void scaleD(double *D_scaled, double *D, double *S, int D_size);
    void fillMatrix_LagrangeDeriv(double Matrix[], int totalEqNum, int compCondNum, int auxiliaryEqNum, int *Ind, int axis, int *eqNum);
    void getLagrangeDeriv_Dirichlet(double *Y, double *Z, int *Ind, int *LInd, int axis);
    void getLagrangeDeriv_Dirichlet(double *Z, int *Ind, int *LInd, int axis);
    void applyQh(double *QV, double *V, double **W, int *Ind, int nu, int k, int axis);
    void getCorrectionTerms(double **VS, double *V, double **W, int *lth, int k, int axis);
    void getBdryRange(int bdryRange[3][2], int bdryNgx[3], int fixedAxis, int fixedSide, int addOnSide0, int addOnSide1);
    int getBdryPointsNum(int fixedAxis);

    /* part2: rhs preparation */
    void prepDataVec(double *R[], double *Rv[], double t, double dt, int axis, int side);
    void getFaceGhost(double *un, double *Rv[], double ***CaVec, double ***CbVec, int *eqNum, int axis, int side);
    void getbInt(double *b, double *un, int *EqNum, int *Ind, int axis, int side);
    void prepDataDeriv(double *R[], double t, double dt, int axis, int side);
    void getDataDeriv(double **R, double t, double dt, int axis, int side);
    void applyQhF(double *QF, double *F, double **W, int bdryRange[3][2], int bdryRangeLth[3], int nu, int k, int axis);
    void getBdryRange(int bdryRange[3][2], int bdryRangeLth[3], int fixedAxis, int fixedSide);
    void freeFaceVariables();
    /* ------------------------------------- */
    
    /* LCBC Part 0 functions */
    static double Default_Cn(double *arg);
    void getCoef(Grid& G);
    double LagrangePoly(int z,int index);
    void getLagrangeData();
    void freeVariables();

    /* variables */
    int p, dim, maxCoefNum, LagrangeData_center, orderInTime;
    const static int q = 1;
    
    LcbcMat FaceMat[6];
    LcbcMat CornerMat[12];
    LcbcMat VertexMat[8];
    
    F1 Cn, Gn, Fn;
    double *LagrangeData, **coef;
    bool prepMatrixFlag[6];
    NumericalDeriv der;
    Grid G;
    
    /* Constructor */
    Lcbc(int pIn, int orderInTimeIn, int *numGridPoints, int numGhost, int dimIn = 2, F1 CnIn = NULL, F1 GnIn = NULL, F1 FnIn = NULL){
        p = pIn;
        orderInTime = orderInTimeIn; 
        dim = dimIn;
        getLagrangeData();
        G.prepareGrid(numGridPoints,numGhost,dim);

        Cn = CnIn;
        Gn = GnIn;
        Fn = FnIn;
        
        maxCoefNum = ((dim+1)*(dim+2)/2);
        getCoef(G);
    }
    
    ~Lcbc(){
        freeVariables();
        freeFaceVariables();
        freeCornerVariables();
        if(dim == 3){
            freeVertexVariables(); 
        }
    }// end of destructor
};
#endif /* LCBC_hpp */
