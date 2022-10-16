#ifndef LCBC_h
#define LCBC_h

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "grid.h"
#include "LCBC_param.h"
#include "LCBCmacros.h"
#include "LCBC_data.h"

typedef double (*F1)(double *);
typedef void (*F2)(double[3], double *);

struct LcbcMat{
    bool flag = false;
    int *eqNum;
    double ***CaVec;
    double ***CbVec;
};

class Lcbc{
public:
    
    /* TESTING */
    
    double mixedDeriv(double *u, int index[3], int deriv[3], double dx[3], int order[3], int size[3]);
    double Deriv(double *u, int index[3], int pos, int deriv, double dx, int order, int size[3]);
    double Deriv(double (*F) (double *), double *arg, int pos, int deriv, double dx, int order);
    double dotProduct(double *v1, double *v2, int lth);
    
    /* Main Function to update ghost */
    void updateGhost(double *&unp1, double t, double dt, LcbcData *&gn, LcbcData *&fn);
    
    /* LCBC 3D Vertex */
    void updateVertexGhost(double **R, double *&unp1, double t, double dt, int side0, int side1, int side2);
    void updateVertexGhostPeriodic(double *&unp1, int side0, int side1, int side2);
    void prepVertexMatrix(int side0, int side1, int side2, int NU[3], int compCondNum[3], int approxEqNum);
    void getVertexMatrix(double **&CaVec, double **&CbVec, int *eqNum, int side0, int side1, int side2, int NU[3], int compCondNum[3]);
    void fillVertexMatrix_LagrangeDeriv(double *&Matrix, int *Ind, int side0, int side1, int side2, int *eqNum, int NU[3], int compCondNum[3], int auxiliaryEqNum, int totalEqNum);
    void getD_vertex(double *&D, int NU[3], int compCondNum[3], int approxEqNum);
    void getbVec_vertex(double b[], double **R, int NU[3], int totalCompCondNum);
    void prepDataVec_vertex(double **R, double *&Rv, double t, double dt, int approxEqNum, int side0, int side1, int side2, int NU[3]);
    void getVertexGhost(double *&un, double *Rv, double **CaVec, double **CbVec, int *eqNum, int approxEqNum, int side0, int side1, int side2);
    void getbInt_vertex(double b[], double *un, int *eqNum, int *Ind, int side0, int side1, int side2);
    void freeVertexVariables();
    /* --------------------------------------*/

    /* LCBC 2D Corner or 3D edge */
    void updateEdgeGhost(double **R, double *&unp1, double t, double dt, int side1, int side2, int varAxis);
    void updateEdgeGhostPeriodic(double *&unp1, int side1, int side2, int varAxis, int fixedAxis[2]);
    void prepCornerMatrix(int side1, int side2, int varAxis, int approxEqNum, int NU1, int NU2);
    void getCornerMatrix(double ***&CaVec, double ***&CbVec, int *eqNum, int bdryRange[3][2], int bdryNg, int side1, int side2, int varAxis, int fixedAxis[2], int NU1, int NU2);
    void getD_corner(double *&D, int varAxis, int fixedAxis[2], int NU1, int NU2);
    void getbVec_corner(double b[], double **R, int varAxis, int fixedAxis[2], int NU1, int NU2);
    void fillCornerMatrix_LagrangeDeriv(double *&Matrix, int *Ind, int side1, int side2, int varAxis, int fixedAxis[2], int *eqNum, int NU1, int NU2, int compCondNum1, int compCondNum2, int auxiliaryEqNum, int totalEqNum);
    void prepDataVec_corner(double **R, double **&Rv, double t, double dt, int side1, int side2, int varAxis, int approxEqNum, int NU1, int NU2);
    void getCornerGhost(double *&un, double **Rv, double ***CaVec, double ***CbVec, int *eqNum, int side1, int side2, int varAxis, int approxEqNum);
    void getbInt_corner(double b[], double *un, int *eqNum, int *Ind, int side1, int side2, int fixedAxis[2]);
    void getBdryRange_corner(int bdryRange[3][2], int fixedAxis[2], int fixedSide[2], int varAxis, int addOnSide0, int addOnSide1);
    void getFixedAxis(int fixedAxis[2], int varAxis);
    void freeCornerVariables();
    /* ------------------------------------- */

    /* LCBC 2D or 3D face  */
    /* part 1: matrix preparation */
    void updateFaceGhost(double *&unp1, double **R, double t, double dt, int axis, int side);
    void prepSideMatrix(int axis, int side);
    void getSideMatrix(double ***&CaVec, double ***&CbVec, int *eqNum, int bdryRange[3][2], int interiorRange[3][2], int bdryNgx[3], int NU, int axis, int side);
    void getD(double *&D, int axis, int side);
    void getbVec(double *&b, double **R, int axis, int side);
    void scaleD(double *&D_scaled, double *D, double *S, int D_size);
    void fillMatrix_LagrangeDeriv(double *&Matrix, int totalEqNum, int compCondNum, int auxiliaryEqNum, int *Ind, int axis, int side, int *eqNum, int NU);
    void getLagrangeDeriv(double Y[], double Z[], int *Ind, int *LInd, int axis, int side, bool evaluateY = true);
    void applyQh(double *&QV, double *V, double **W, int *Ind, int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side);
    void applyDh(double *&QV, double *V, double **W, int *Ind, int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side);
    void getCorrectionTerms(double **&VS, double *V, double **W, int *lth, int k, int axis);
    void getCorrectionTermsForDh(double **&VS, double *V, double **W, int *lth, int k, int axis);

    /* part2: rhs preparation */
    void prepDataVec(double **R, double **&Rv, double t, double dt, int bdryRange[3][2], int bdryNgx[3], int axis, int side);
    void getFaceGhost(double *&un, double *Rv[], double ***CaVec, double ***CbVec, int *eqNum, int LcbcBdryRange[3][2], int bdryNgx[3], int axis, int side);
    void getbInt(double b[], double *un, int *EqNum, int *Ind, int interiorRange[3][2], int axis, int side);
    void getDataDeriv(double **&R, double t, double dt, int axis, int side, LcbcData &gn, LcbcData &fn);
    void applyQhF(double *&QF, double *F, double **W, int bdryRange[3][2], int bdryRangeLth[3], int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side);
    void applyDhF(double *&QF, double *F, double **W, int bdryRange[3][2], int bdryRangeLth[3], int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side);
    void freeFaceVariables();
    /* ------------------------------------- */
    
    /* LCBC Part 0 functions */
    void initialize(int dimIn, int orderInSpace, int orderInTimeIn, int *numGridPoints, int numGhostIn, int *faceEvalIn, LcbcData *coefIn = NULL, F1 CnIn = NULL, F1 GnIn = NULL, F1 FnIn = NULL, F2 bnIn = NULL, bool cstCoefIn = false, bool *zeroBCIn = NULL, bool noForcingIn = false);
    void analyzeUserInput(LcbcData *gn, LcbcData *fn);
    
    static double Default_Cn(double *arg);
    static double Default_Gn(double *arg); bool defaultedGn = false;
    static double Default_Fn(double *arg); bool defaultedFn = false;
    void getCoef();
    void getCoefGridLth(int indexRange[3][2], int lth[3], int wth[3], int fixedAxis, int fixedSide, int dim, int p);
    void getCoefGridLth(int indexRange[3][2], int lth[3], int wth[3], int bdryRange[3][2], int fixedAxis, int fixedSide, int dim, int p);
    double LagrangePoly(int z,int index);
    void getLagrangeData();
    void freeVariables();
    
    /* variables */
    int p, dim, numGhost, userNumGhost, maxCoefNum, LagrangeData_center, orderInTime;
    const static int q = 1;
    int extraDataGhost;
    
    Param param; FaceParam faceParam[6];
    LcbcMat FaceMat[6];
    LcbcMat CornerMat[12];
    LcbcMat VertexMat[8];

    
    F1 Cn, Gn, Fn; // Cn: coef, Gn: bdry, Fn: forcing, bn: bdryOperator
    F2 bn;
    double *LagrangeData;
    bool noForcing, *zeroBC, zeroBCnewed = false, cstCoef, isInitialized = false, faceEvalNewed = false, initializedBdryData = false, initializedForcingData = false, analyzedUserData = false;
    int *faceEval;
    Grid G;
    LcbcData *coef, *bdryData = NULL, *forcingData = NULL;
    enum LcbcDataType{boundary = 1, forcing = 2, coefficient = 3};

    bool debug = false;
    
    /* Constructor */
    
    Lcbc(){}
        
    Lcbc(int dimIn, int orderInSpace, int orderInTimeIn, int *numGridPoints, int numGhostIn, int *faceEvalIn = NULL, LcbcData *coefIn = NULL, F1 CnIn = NULL, F1 GnIn = NULL, F1 FnIn = NULL, F2 bnIn = NULL, bool cstCoefIn = false, bool *zeroBCIn = NULL, bool noForcingIn = false){
        
        initialize(dimIn, orderInSpace, orderInTimeIn, numGridPoints, numGhostIn, faceEvalIn, coefIn, CnIn, GnIn, FnIn, bnIn, cstCoefIn, zeroBCIn, noForcingIn);
    }
    
    ~Lcbc(){
        if(isInitialized){
            freeVariables();
            freeFaceVariables();
            freeCornerVariables();
            if(dim == 3){
                freeVertexVariables();
            }// end of if dim
        }// end of if initialized
    }// end of destructor
};

inline double Lcbc::mixedDeriv(double *u, int index[3], int deriv[3], double dx[3], int order[3], int size[3]){

    double D = 0;

    int r[3] = {param.derivRange[pInd((order[0]/2),deriv[0])],
                param.derivRange[pInd((order[1]/2),deriv[1])],
                param.derivRange[pInd((order[2]/2),deriv[2])]};
    
    int cnt[3] = {0,0,0};
    int i[3];
    for(i[2] = (-r[2]); i[2]<=r[2]; i[2]++){
        for(i[1] = (-r[1]); i[1]<=r[1]; i[1]++){
            for(i[0] = (-r[0]); i[0]<=r[0]; i[0]++ ){
                int Ind[3] = {(index[0] + i[0]),(index[1] + i[1]),(index[2] + i[2])};
                double c[3] = {param.derivCoef[pInd((order[0]/2),deriv[0])][cnt[0]],
                               param.derivCoef[pInd((order[1]/2),deriv[1])][cnt[1]],
                               param.derivCoef[pInd((order[2]/2),deriv[2])][cnt[2]]};
                
                D = D + c[0]*c[1]*c[2]*u[ind(Ind,size)];

                cnt[0]++;
            }// end of i2
            cnt[0] = 0;
            cnt[1]++;
        }// end of i1
        cnt[0] = 0;
        cnt[1] = 0;
        cnt[2]++;
    }// end of i0
    
    D = D*(1.0/pow(dx[0],(double) deriv[0]))*(1.0/pow(dx[1],(double) deriv[1]))*(1.0/pow(dx[2],(double) deriv[2]));

    return D;
}// end of mixedDeriv

inline double Lcbc::Deriv(double (*F) (double *), double *arg, int pos, int deriv, double dx, int order){
    double D = 0;

    int r = param.derivRange[pInd((order/2),deriv)];

    double argVar = arg[(pos+1)];
    int cnt = 0;
    for(int i = (-r); i<=r; i++ ){
        arg[(pos + 1)] = argVar + dx*i;
        D = D + param.derivCoef[pInd((order/2),deriv)][cnt]*F(arg);
        cnt++;
    }// end of i2
    
    arg[(pos + 1)] = argVar;
    D = D*(1.0/pow(dx,(double) deriv));

    return D;
}// end of Deriv

inline double Lcbc::Deriv(double *u, int index[3], int pos, int deriv, double dx, int order, int size[3]){
    double D = 0;

    int r = param.derivRange[pInd((order/2),deriv)];

    int varIndex = index[pos];
    int cnt = 0;
    for(int i = (-r); i<=r; i++ ){
        index[pos] = varIndex + i;
        D = D + param.derivCoef[pInd((order/2),deriv)][cnt]*u[ind(index,size)];
        cnt++;
    }// end of i2
    index[pos] = varIndex;
    
    D = D*(1.0/pow(dx,(double) deriv));

    return D;
}// end of Deriv

inline double Lcbc::dotProduct(double *v1, double *v2, int lth){
    double dotProduct = 0;
    for(int i = 0; i<lth; i++){
        dotProduct = dotProduct + v1[i]*v2[i];
    }// end of i loop
    return dotProduct;
}// end of dotProduct

inline void Lcbc::getbInt(double b[], double *un, int *EqNum, int *Ind, int interiorRange[3][2], int axis, int side){
    
    int n = (2*p+1);
    int unknownVarNum = param.unknownVarNum;

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


#endif /* LCBC_h */
