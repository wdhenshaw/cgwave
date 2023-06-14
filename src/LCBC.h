#ifndef LCBC_h
#define LCBC_h

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "grid.h"
#include "LCBC_param.h"
#include "LCBC_data.h"
//#include "debuggingMacros.h"

/* Function handles used to identify input data functions (forcing, PDE coefficients, BCs)*/
typedef double (*F1)(double *);
typedef void (*F2)(double[3], double *);

/// \brief LcbcMat is a structure that holds the LCBC stencil coefficients and other needed information for the LCBC procedure. The LCBC class creates an LcbcMat object at each face, edge or vertex where the LCBC procedure will be implemented
struct LcbcMat{
    bool flag = false; // to tell if the LcbcMat is created and space is allocated
    int *eqNum;     // the ordering of coefficients of the LCBC interpolating polynomial such that the unknown coefficients come first
    double **CaVec; // the rows of C_alpha needed for the stencil formula to evaluate the solution at ghost points (see lcbc.pdf stencil formula) at each boundary point
    double **CbVec; // the rows of C_beta needed for the stencil formula to evaluate the solution at ghost points (see lcbc.pdf stencil formula) at each boundary point
    
    int bdryRange[3][2]; // the starting and ending indices at the boundary
    int bdryNgx[3] = {1,1,1}; // the number of grid points on the boundary in each axis
    int interiorRange[3][2];  // the range of indices corresponding to the part of the stencil where the solution is known on the interior or boundary
    int ghostRange[3][2]; // the range of indices at ghost points where the solution needs to be determined
};

class Lcbc{
private:
    
    /* ----- variables ----- */
    int p, dim, numGhost, userNumGhost, maxCoefNum, LagrangeData_center, orderInTime, auxiliaryEqNum, faceNum, totalVarNum;
    // LagrangeData_center: the width of the stencil for the Lagrange data functions prepared
    // auxiliaryEqNum: the number of auxiliary equations added in the LCBC implementation if any
    // totalVarNum: the total number of coefficients of the interpolating polynomial
    const static int q = 2;
    // q: the order of the temporal operator of the PDE. Here, q = 2 because we're solving the wave equation
    int extraDataGhost; // extra ghost values needed for data grid functions
    
    FaceParam faceParam[6]; // An object defined in LCBC_param.h which holds parameters needed throughout the LCBC process specific to each face
    
    LcbcMat FaceMat[6]; // LcbcMat object (defined at the top of LCBC.h) at faces
    LcbcMat EdgeMat[12]; // LcbcMat object at 3D edges or 2D corners
    LcbcMat VertexMat[8]; // LcbcMat object at 3D vertices
    
    /* bdryTypeParam: an object that carries information related to the boundary type (defined in LCBC_param.h) */
    bdryTypeParam oneAxis; // one fixed axis, one known axis
    bdryTypeParam twoAxis; // two fixed axes, two known axes
    bdryTypeParam threeAxis; // three fixed axes, three known axes
    bdryTypeParam nearEdge; // two fixed axes, one known axes
    bdryTypeParam nearVertex[2]; // three fixed axes, [0] one known axes, [1] two known axes
    
    F1 Cn, Gn, Fn; // Cn: coef, Gn: bdry, Fn: forcing, bn: bdryOperator
    F2 bn; // mapped Neumann boundary operator
    double *LagrangeData; // holds data from Lagrange basis polynomial evaluated at integer values
    bool noForcing, *zeroBC, zeroBCnewed = false, cstCoef, isInitialized = false, faceEvalNewed = false, initializedBdryData = false, initializedForcingData = false, analyzedUserData = false, addAuxEqns = false;
    bool defaultedGn = false;
    bool defaultedFn = false;
    
    // isInitialized is boolean to identify if LCBC object is initialized or not.
    // addAuxEqns is a boolean to tell if we need to set high-order derivatives of the LCBC interpolating polynomial to zero
    
    int *faceEval; // holds the type of BCs on each face
    // faceEval[face] = -1: Periodic BC
    // faceEval[face] =  0: No BC
    // faceEval[face] =  1: Dirichlet BC
    // faceEval[face] =  2: Neumann BC
    
    Grid G; // holds infomation about the grid
    LcbcData *coef, *bdryData = NULL, *forcingData = NULL; // LcbcData object holds PDE coefficient data, boundary and forcing data (see LCBC_data.h where the LcbcData class is defined)
    enum LcbcDataType{boundary = 1, forcing = 2, coefficient = 3};
    char **memory; int memoryComponents; bool preallocatedMemory = false; // memory objects needed if we want to preallocated memory for certain functions
    
    /* the following parameters are needed for numerical derivatives */
    static const int orderNum = 5; // orders = 0, 2, 4, 6
    static const int derivNum = 9; // derivatives = 0, 1, ..., 6
    static const int maxStencilWth = 9;
    double derivCoef[(orderNum*derivNum)][maxStencilWth]; int derivRange[(orderNum*derivNum)];
    
    /* ----- LCBC_deriv.C ----- */
    void getDerivCoef(int p);
    double mixedDeriv(double *u, int index[3], int deriv[3], double dx[3], int order[3], int size[3]);
    double Deriv(double *u, int index[3], int pos, int deriv, double dx, int order, int size[3]);
    double Deriv(double (*F) (double *), double *arg, int pos, int deriv, double dx, int order);
    
    /* ----- LCBC1.C ----- */
    void prepLcbcStencil(LcbcMat &Mat, bdryTypeParam &bdryParam, int fixedAxis[], int fixedSide[], int NU[], int compCondNum[]);
    void getD(double *&D, int fixedAxis[], int fixedSide[], int NU[], int compCondNum[], int approxEqNum, int knownAxesNum);
    void getbVec(double b[], double **R, int NU[3], int totalCompCondNum, int knownAxesNum, int fixedAxis[]);
    void getLagrangeDeriv(double Y[], double Z[], int *Ind, int *LInd, int axis, int side, bool evaluateY = true);
    void applyQh(double *&QV, double *V, double **W, int *Ind, int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side);
    void applyDh(double *&QV, double *V, double **W, int *Ind, int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side);
    void getCorrectionTerms(double **&VS, double *V, double **W, int *lth, int k, int axis);
    void getCorrectionTermsForDh(double **&VS, double *V, double **W, int *lth, int k, int axis);
    void getLcbcStencilCoef(double **&CaVec, double **&CbVec, int *eqNum, int bdryRange[3][2], int bdryNgx[3], int ghostIndRange[3][2], int knownAxesNum, bdryTypeParam &bdryParam, int fixedAxis[], int fixedSide[], int NU[], int compCondNum[], int approxEqNum);
    void getBlockMatrices(double *&A11, double *&A12, int *Ind, int *eqNum, bdryTypeParam &bdryType, int NU[3], int compCondNum[3], int approxEqNum, LagrangeDerivFun LagrangeDeriv[], int fixedSide[], int fixedAxis[]);
    
    /* ----- LCBC2.C ----- */
    void prepDataVec(double **R, double **&Rv, double t, double dt, int bdryRange[3][2], int bdryNgx[3], bdryTypeParam &bdryParam, int fixedAxis[], int fixedSide[], int approxEqNum, int NU[]);
    void getGhost(double *&un, int *mask, double *Rv[], LcbcMat &Mat, bdryTypeParam &bdryParam, int fixedAxis[], int fixedSide[], int approxEqNum);
    void getbInt(double b[], double *un, int *EqNum, int *Ind, int interiorRange[3][2], int unknownVarNum);
    void getDataDeriv(double **&R, double t, double dt, int axis, int side, LcbcData &gn, LcbcData &fn);
    void applyQhF(double *&QF, double *F, double **W, int bdryRange[3][2], int bdryRangeLth[3], int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side);
    void applyDhF(double *&QF, double *F, double **W, int bdryRange[3][2], int bdryRangeLth[3], int nu, int k, int newWth[3], int newLth[3], int oldWth[3], int oldLth[3], int axis, int side);
    
    /* ----- LCBC.C ----- */
    void updateNeededParameters();
    void preallocateMemory(int components, int componentSize);
    void deleteMemory(int components);
    void analyzeUserInput(LcbcData *gn, LcbcData *fn);
    void getCoef();
    void getCoefGridLth(int indexRange[3][2], int lth[3], int wth[3], int bdryRange[3][2], int fixedAxis, int fixedSide, int dim, int p);
    double LagrangePoly(int z,int index);
    void getLagrangeData();
    void freeVariables();
    void freeFaceVariables();
    void freeEdgeVariables();
    void freeVertexVariables();
    
    /* ----- LCBC_updateGhost.C ----- */
    void updateBdryGhost(double *&unp1, int *mask, double **R, bdryTypeParam &bdryType, LcbcMat &Mat, double t, double dt, int fixedAxis[], int fixedSide[]);
    void updateVertexGhostPeriodic(double *&unp1, int *mask, int side0, int side1, int side2);
    void updateEdgeGhostPeriodic(double *&unp1, int *mask, int side1, int side2, int varAxis, int fixedAxis[2]);
    void getBdryRange_corner(int bdryRange[3][2], int fixedAxis[2], int fixedSide[2], int varAxis, int addOnSide0, int addOnSide1);
    void updateGhostNearZeroBdry(double *&unp1, int *mask, double **R, bdryTypeParam &bdryType, LcbcMat &Mat, double t, double dt, int fixedAxis[], int fixedSide[]);
    
public:
    
    /* ----- Constructor ----- */
    Lcbc(){}
    
    /// \brief Constructor: see the documentation of initialize in LCBC.C
    Lcbc(int dimIn, int orderInSpace, int orderInTimeIn, int *numGridPoints, int numGhostIn, int *faceEvalIn = NULL, LcbcData *coefIn = NULL, F1 CnIn = NULL, F1 GnIn = NULL, F1 FnIn = NULL, F2 bnIn = NULL, bool cstCoefIn = false, bool *zeroBCIn = NULL, bool noForcingIn = false){
        
        initialize(dimIn, orderInSpace, orderInTimeIn, numGridPoints, numGhostIn, faceEvalIn, coefIn, CnIn, GnIn, FnIn, bnIn, cstCoefIn, zeroBCIn, noForcingIn);
    }
    
    /* ----- LCBC.C ----- */
    void initialize(int dimIn, int orderInSpace, int orderInTimeIn, int *numGridPoints, int numGhostIn, int *faceEvalIn, LcbcData *coefIn = NULL, F1 CnIn = NULL, F1 GnIn = NULL, F1 FnIn = NULL, F2 bnIn = NULL, bool cstCoefIn = false, bool *zeroBCIn = NULL, bool noForcingIn = false);
    void updateGhost(double *&unp1, int *mask, double t, double dt, LcbcData *&gn, LcbcData *&fn);
    void setExactGhost(double *&un, double *ue, int *mask);
    
    /* ----- Destructor ----- */
    /// \brief Destructor: free allocated memory of variables
    ~Lcbc(){

        if(isInitialized){
            freeVariables();
            freeFaceVariables();
            freeEdgeVariables();
            if(dim == 3){
                freeVertexVariables();
            }// end of if dim
        }// end of if initialized
    }// end of destructor
};

#endif /* LCBC_h */
