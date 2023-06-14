#include <stdio.h>
#include <string.h>
#include <math.h>
#include "LCBC.h"
#include "utility.h"
#include "LCBCmacros.h"

void getGhostIndRange(int ghostIndRange[3][2], int fixedAxis[], int fixedSide[], int fixedAxesNum, int knownAxesNum, int p, int dim);
void getInteriorIndRange(int interiorRange[3][2], int fixedAxis[], int fixedSide[], int knownAxesNum, int p, int dim);
void getBdryRange(int bdryRange[3][2], int LcbcBdryRange[3][2], int faceBdryRange[3][2], int bdryNgx[3], int faceBdryNgx[3], int fixedAxis[], int fixedSide[], int fixedAxesNum);

/// \brief update the ghost values on an edge or a vertex formed with at least one interpolation boundary
/// \param unp1 (input/output): solution grid function
/// \param mask (input): mask to determine grid point type
/// \param R (input): vector holding values of data functions involved in the LCBC procedure. Used to compute R(t) in the stencil approach in the lcbc.pdf paper
/// \param bdryType (input): carries parameters depending on the number of fixed axes and axes where boundary data is known
/// \param Mat (input): An object defined in LCBC.h that carries the LCBC stencil coefficients and other information needed in the LCBC procedure
/// \param t (input): the time at which the solution will be updated
/// \param dt (input): the time-step
/// \param fixedAxis (input): a vector containing the fixed axes at the specific edge or vertex
/// \param fixedSide (input): a vector containing the fixed sides at the specific edge or vertex
void Lcbc::updateGhostNearZeroBdry(double *&unp1, int *mask, double **R, bdryTypeParam &bdryType, LcbcMat &Mat, double t, double dt, int fixedAxis[], int fixedSide[]){
    
    int knownFaceNum = bdryType.knownAxesNum; // the number of axes where boundary data is provided
    
    /* known Axis: fixed axis at the boundary where boundary data is known */
    int sortedAxis[3] = {(-1),(-1),(-1)}; // sorted by known axis first then unknown
    int sortedSide[3] = {(-1),(-1),(-1)}; // sorted by known side first then unknown
    
    int NU[2] = {0,0}, compCondNum[2] = {0,0};
    // NU: the number of primary CBCs at each known axis
    // compCondNum: the number of CBCs and their tangential derivatives at each known axis
    int knownCnt = 0, unknownCnt = 0; // numbers of known and unknown axes

    /* sort the fixed axes starting from known axes */
        for(int axis = 0; axis<bdryType.fixedAxesNum; axis++){
            int face = fixedSide[axis] + 2*fixedAxis[axis];
            if(faceEval[face]>0){
                NU[knownCnt] = faceParam[face].NU;
                compCondNum[knownCnt] =  faceParam[face].compCondNum;
                sortedAxis[knownCnt] = fixedAxis[axis];
                sortedSide[knownCnt] = fixedSide[axis];
                knownCnt++;
            }
            else if(faceEval[face] == 0){
                int axisInd = knownFaceNum + unknownCnt;
                sortedAxis[axisInd] = fixedAxis[axis];
                sortedSide[axisInd] = fixedSide[axis];
                unknownCnt++;
            }
        }

    /* The number of equations that will be approximated using the least-squares procedure */
    int approxEqNum  = intVectorSum(compCondNum, knownFaceNum) + auxiliaryEqNum;
    
    /* prepare the matrices for the corner */
    if(Mat.flag == false){

        int donorFace = sortedSide[0] + 2*sortedAxis[0]; // one of the known axes that will provide the boundary data information for the LCBC procedure
        
        /* get the boundary range */
        getBdryRange(Mat.bdryRange, faceParam[donorFace].LcbcBdryRange, faceParam[donorFace].maskBdryRange, Mat.bdryNgx, faceParam[donorFace].bdryNgx, sortedAxis, sortedSide, bdryType.fixedAxesNum);
        
#if PRINT_PROCESS == 1
        printf("bdryRange: ");
        printIndexRange(Mat.bdryRange);
#endif
        /* get the range of indices where the solution is known on the interior and the boundary */
        getInteriorIndRange(Mat.interiorRange, sortedAxis, sortedSide, bdryType.knownAxesNum, p, dim);
        
        /* get the range of indices where the solution is to be determined on the ghost points */
        getGhostIndRange(Mat.ghostRange, sortedAxis, sortedSide, bdryType.fixedAxesNum, bdryType.knownAxesNum, p, dim);
        
        /* Prepare the LCBC stencil coefficients */
        prepLcbcStencil(Mat, bdryType, sortedAxis, sortedSide, NU, compCondNum);
    }
    
    int bdryPointsNum = Mat.bdryNgx[0]*Mat.bdryNgx[1]*Mat.bdryNgx[2];
    double **Rv = new double*[bdryPointsNum];
    
    for(int bdryPt = 0; bdryPt<bdryPointsNum; bdryPt++){
        Rv[bdryPt] = new double[approxEqNum];
    }
    
    /* prepare the RHS vector for the near corner */
    prepDataVec(R, Rv, t, dt, Mat.bdryRange, Mat.bdryNgx, bdryType, sortedAxis, sortedSide, approxEqNum, NU);
    
    /* evaluate the ghost points */
    getGhost(unp1, mask, Rv, Mat, bdryType, sortedAxis, sortedSide, approxEqNum);
    
    /* Free any allocated variables */
    for(int bdryPoint = 0; bdryPoint<bdryPointsNum; bdryPoint++)
        delete [] Rv[bdryPoint];
    
    delete [] Rv;
}

/// \brief update the ghost values on an edge or a vertex where all the fixed axis correspond to physical boundaries where the boundary data is known
/// \param unp1 (input/output): solution grid function
/// \param mask (input): mask to determine grid point type
/// \param R (input): vector holding values of data functions involved in the LCBC procedure. Used to compute R(t) in the stencil approach in the lcbc.pdf paper
/// \param bdryType (input): carries parameters depending on the number of fixed axes and axes where boundary data is known
/// \param Mat (input): An object defined in LCBC.h that carries the LCBC stencil coefficients and other information needed in the LCBC procedure
/// \param t (input): the time at which the solution will be updated
/// \param dt (input): the time-step
/// \param fixedAxis (input): a vector containing the fixed axes at the specific edge or vertex
/// \param fixedSide (input): a vector containing the fixed sides at the specific edge or vertex
void Lcbc::updateBdryGhost(double *&unp1, int *mask, double **R, bdryTypeParam &bdryType, LcbcMat &Mat, double t, double dt, int fixedAxis[], int fixedSide[]){
    
    /* find the number of axes where the data is known at the boundary */
    int knownAxesNum = bdryType.knownAxesNum;

    /* Prepare parameters needed in the Lcbc process */
    int fixedFace[knownAxesNum];
    int NU[knownAxesNum]; // the number of primary CBCs
    int compCondNum[knownAxesNum]; // the number of CBCs and their tangential derivatives
    for(int face = 0; face<knownAxesNum; face++){
        fixedFace[face] = fixedSide[face] + 2*fixedAxis[face];
        NU[face] = faceParam[fixedFace[face]].NU;
        compCondNum[face] = faceParam[fixedFace[face]].compCondNum;
    }
    int approxEqNum = intVectorSum(compCondNum, knownAxesNum) + auxiliaryEqNum; // the number of equations approximated using least-squares

    /* prepare the LCBC matrix*/
    if(Mat.flag == false){
        getBdryRange(Mat.bdryRange, faceParam[fixedFace[0]].LcbcBdryRange, faceParam[fixedFace[0]].bdryRange, Mat.bdryNgx, faceParam[fixedFace[0]].bdryNgx, fixedAxis, fixedSide, bdryType.fixedAxesNum);
        
#if PRINT_PROCESS == 1
        printf("bdryRange: ");
        printIndexRange(Mat.bdryRange);
#endif
        
        getInteriorIndRange(Mat.interiorRange, fixedAxis, fixedSide, bdryType.knownAxesNum, p, dim);
        getGhostIndRange(Mat.ghostRange, fixedAxis, fixedSide, bdryType.fixedAxesNum, bdryType.knownAxesNum, p, dim);
        
        /* prepare the stencil coefficients used in the LCBC procedure (at the first time-step only) */
        prepLcbcStencil(Mat, bdryType, fixedAxis, fixedSide, NU, compCondNum);
    }
    
    int bdryPoints = Mat.bdryNgx[0]*Mat.bdryNgx[1]*Mat.bdryNgx[2];
    double **Rv = new double*[bdryPoints];

    for(int bdryPt = 0; bdryPt<bdryPoints; bdryPt++){
        Rv[bdryPt] = new double[approxEqNum];
    }

    /* prepare the RHS data vector */
    prepDataVec(R, Rv, t, dt, Mat.bdryRange, Mat.bdryNgx, bdryType, fixedAxis, fixedSide, approxEqNum, NU);
    
    /* Get the solution values at the ghost points */
    getGhost(unp1, mask, Rv, Mat, bdryType, fixedAxis, fixedSide, approxEqNum);

    /* Free any allocated variables */
    for(int bdryPoint = 0; bdryPoint<bdryPoints; bdryPoint++)
        delete [] Rv[bdryPoint];

    delete [] Rv;
}// end of update FaceGhost

/// \brief Update solution values near a 2D corner or a 3D edge with at least one periodic BC
/// \param unp1 (input/output): the solution at the current time-step
/// \param side1 (input): the side corresponding to the first face
/// \param side2 (input): the side corresponding to the second face
/// \param varAxis (input): the axis corresponding to the variable that varies on the edge if applicable
/// \param fixedAxis (input): a vector containing the fixed axes
void Lcbc::updateEdgeGhostPeriodic(double *&unp1, int *mask, int side1, int side2, int varAxis, int fixedAxis[2]){
    
    /* At the meeting of two face, determine which edge is periodic */
    int fixedSides[2] = {side1, side2};
    int periodicAxis = -1, perAxisInd = -1;
    for(int i = 0; i<2; i++){
        int edgeFace = fixedSides[i] + 2*fixedAxis[i];
        if(faceEval[edgeFace] == -1){
            periodicAxis = fixedAxis[i];
            perAxisInd = i;
        }// end of faceEval
    }// end of i
    
    /* Ensure that indeed one of the edges is periodic */
    assert(periodicAxis>=0 && perAxisInd>=0);
    
    /* Find the indices, in each axes, needed at this edge and put them in bdryRange. */
    int bdryRange[3][2]; getBdryRange_corner(bdryRange, fixedAxis, fixedSides, varAxis, 0, 0);
    
    /* For the periodic implementation, adjust the bdryRange above to support extra evaluations */
    int perIndexRange[3][2];
    perIndexRange[fixedAxis[0]][0] = bdryRange[fixedAxis[0]][0] + sideBasedValue(side1, (-p), 1);
    perIndexRange[fixedAxis[0]][1] = bdryRange[fixedAxis[0]][1] + sideBasedValue(side1, (-1), p);
    perIndexRange[fixedAxis[1]][0] = bdryRange[fixedAxis[1]][0] + sideBasedValue(side2, (-p), 1);
    perIndexRange[fixedAxis[1]][1] = bdryRange[fixedAxis[1]][1] + sideBasedValue(side2, (-1), p);
    perIndexRange[varAxis][0]      = bdryRange[varAxis][0];
    perIndexRange[varAxis][1]      = bdryRange[varAxis][1];
    
    int i[3];
    for(i[2] = perIndexRange[2][0]; i[2]<=perIndexRange[2][1]; i[2]++){
        for(i[1] = perIndexRange[1][0]; i[1]<=perIndexRange[1][1]; i[1]++){
            for(i[0] = perIndexRange[0][0]; i[0]<=perIndexRange[0][1]; i[0]++){
                
                if(mask[solInd(i, G.Ngx)]>0){
                    int I[3]; memcpy(I, i, sizeof(I));
                    
                    I[periodicAxis] = I[periodicAxis] + sideBasedValue(fixedSides[perAxisInd], G.Nx[periodicAxis], (-G.Nx[periodicAxis]));
                    
                    unp1[solInd(i, G.Ngx)] = unp1[solInd(I, G.Ngx)];
                }
                
            }// end of i[0]
        }// end of i[1]
    }// end of i[2]
}// end of updateEdgeGhostPeriodic

/// \brief Update solution values near a  3D vertex with at least one periodic BC
/// \param unp1 (input/output): the solution at the current time-step
/// \param side0 (input): the side corresponding to the first face
/// \param side1 (input): the side corresponding to the second face
/// \param side2 (input): the side corresponding to the third face
void Lcbc::updateVertexGhostPeriodic(double *&unp1, int *mask, int side0, int side1, int side2){
    
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
                
                
                if( mask[solInd(i, G.Ngx)]>0){
                    
                    for(int axis = 0; axis<3; axis ++){
                        if(faceEval[vertices[axis]] == -1){ // periodic boundary
                            I[axis] = I[axis] + sideBasedValue(fixedSides[axis], G.Nx[axis], (-G.Nx[axis]));
                        }
                    }// end of axis loop
                    
                    unp1[solInd(i, G.Ngx)] = unp1[solInd(I, G.Ngx)];
                }

            }// end of i[0]
        }// end of i[1]
    }// end of i[2]
}// end of updateVertexGhostPeriodic

/// \brief Find the range of indices at the corner or edge bdryRange
/// \param bdryRange (output)
/// \param addOnSide0 (input): number of extra grid points added at side 0 if needed
/// \param addOnSide1 (input): number of extra grid points added at side 1 if needed
void Lcbc::getBdryRange_corner(int bdryRange[3][2], int fixedAxis[2], int fixedSide[2], int varAxis, int addOnSide0, int addOnSide1){
    for(int side = 0; side<2; side++){
        bdryRange[fixedAxis[0]][side] = G.indexRange[fixedAxis[0]][fixedSide[0]];
        bdryRange[fixedAxis[1]][side] = G.indexRange[fixedAxis[1]][fixedSide[1]];
        if(dim == 3){
            if(faceEval[(side + 2*varAxis)]>0){
                bdryRange[varAxis][side] = G.indexRange[varAxis][side] + (1-side)*addOnSide0 + addOnSide1*side;
            }
//            else if(faceEval[(side + 2*varAxis)]==0){
//                bdryRange[varAxis][side] = maskRange[varAxis][side];
//            }
            else{
                bdryRange[varAxis][side] = G.indexRange[varAxis][side];
            }
        }else{
            bdryRange[2][side] = G.indexRange[2][side];
        }
    }
}// end of getBdryRange_corner

/// \brief get the range  of indices at the ghost points where the solution is needed
/// \param ghostIndRange (input/output): a vector to carry the range of ghost indices
/// \param fixedAxis (input): the axes that are fixed on the specific face, edge or corner
/// \param fixedSide (input): the fixed sides at the specific face, edge or corner
/// \param fixedAxesNum (input): the number of fixed axes
/// \param knownAxesNum (input): the number of fixed axes where the boundary data is known
/// \param p (input): orderInSpace/2
/// \param dim (input): the dimension in space
void getGhostIndRange(int ghostIndRange[3][2], int fixedAxis[], int fixedSide[], int fixedAxesNum, int knownAxesNum, int p, int dim){
    if(fixedAxesNum == 1 && knownAxesNum == 1){
        ghostIndRange[0][0] = p;
        ghostIndRange[0][1] = p;
        ghostIndRange[1][0] = p;
        ghostIndRange[1][1] = p;
        ghostIndRange[2][0] = dimBasedValue(dim, 0, p);
        ghostIndRange[2][1] = dimBasedValue(dim, 0, p);
        ghostIndRange[fixedAxis[0]][0] = sideBasedValue(fixedSide[0], 0, (1+p));
        ghostIndRange[fixedAxis[0]][1] = sideBasedValue(fixedSide[0], (p-1), (2*p));
    }
    else if(fixedAxesNum == 2 && knownAxesNum == 1){
        ghostIndRange[0][0] = p;
        ghostIndRange[0][1] = p;
        ghostIndRange[1][0] = p;
        ghostIndRange[1][1] = p;
        ghostIndRange[2][0] = dimBasedValue(dim, 0, p);
        ghostIndRange[2][1] = dimBasedValue(dim, 0, p);

        ghostIndRange[fixedAxis[0]][0] = sideBasedValue(fixedSide[0], 0, (1+p));
        ghostIndRange[fixedAxis[0]][1] = sideBasedValue(fixedSide[0], (p-1), (2*p));
        ghostIndRange[fixedAxis[1]][0] = sideBasedValue(fixedSide[1], 0, (1+p));
        ghostIndRange[fixedAxis[1]][1] = sideBasedValue(fixedSide[1], (p-1), (2*p));
    }
    else if(fixedAxesNum == 2 && knownAxesNum == 2){
        ghostIndRange[0][0] = p;
        ghostIndRange[0][1] = p;
        ghostIndRange[1][0] = p;
        ghostIndRange[1][1] = p;
        ghostIndRange[2][0] = dimBasedValue(dim, 0, p);
        ghostIndRange[2][1] = dimBasedValue(dim, 0, p);
        
        ghostIndRange[fixedAxis[0]][0] = sideBasedValue(fixedSide[0], 0, 1);
        ghostIndRange[fixedAxis[0]][1] = sideBasedValue(fixedSide[0], (2*p-1), (2*p));
        ghostIndRange[fixedAxis[1]][0] = sideBasedValue(fixedSide[1], 0, 1);
        ghostIndRange[fixedAxis[1]][1] = sideBasedValue(fixedSide[1], (2*p-1), (2*p));
    }
    else if(fixedAxesNum == 3 && knownAxesNum == 3){
        ghostIndRange[0][0] = sideBasedValue(fixedSide[0], 0, 1);
        ghostIndRange[0][1] = sideBasedValue(fixedSide[0], (2*p-1), (2*p));
        ghostIndRange[1][0] = sideBasedValue(fixedSide[1], 0, 1);
        ghostIndRange[1][1] = sideBasedValue(fixedSide[1], (2*p-1), (2*p));
        ghostIndRange[2][0] = sideBasedValue(fixedSide[2], 0, 1);
        ghostIndRange[2][1] = sideBasedValue(fixedSide[2], (2*p-1), (2*p));
    }
    else if(((fixedAxesNum == 3) && (knownAxesNum == 1))){
        ghostIndRange[fixedAxis[0]][0] = sideBasedValue(fixedSide[0], 0, (1+p));
        ghostIndRange[fixedAxis[0]][1] = sideBasedValue(fixedSide[0], (p-1), (2*p));
        ghostIndRange[fixedAxis[1]][0] = sideBasedValue(fixedSide[1], 0, (1+p));
        ghostIndRange[fixedAxis[1]][1] = sideBasedValue(fixedSide[1], (p-1), (2*p));
        ghostIndRange[fixedAxis[2]][0] = sideBasedValue(fixedSide[2], 0, (1+p));
        ghostIndRange[fixedAxis[2]][1] = sideBasedValue(fixedSide[2], (p-1), (2*p));
    }
    else if(((fixedAxesNum == 3) && (knownAxesNum == 2))){
        ghostIndRange[fixedAxis[0]][0] = sideBasedValue(fixedSide[0], 0, 1);
        ghostIndRange[fixedAxis[0]][1] = sideBasedValue(fixedSide[0], (2*p-1), (2*p));
        ghostIndRange[fixedAxis[1]][0] = sideBasedValue(fixedSide[1], 0, 1);
        ghostIndRange[fixedAxis[1]][1] = sideBasedValue(fixedSide[1], (2*p-1), (2*p));
        ghostIndRange[fixedAxis[2]][0] = sideBasedValue(fixedSide[2], 0, (1+p));
        ghostIndRange[fixedAxis[2]][1] = sideBasedValue(fixedSide[2], (p-1), (2*p));
    }
}

/// \brief Get the range  of indices on the stencil where the solution is known
/// \param interiorRange (input/output): a vector to carry the range of indices where the solution is known on the interior
/// \param fixedAxis (input): the axes that are fixed on the specific face, edge or corner
/// \param fixedSide (input): the fixed sides at the specific face, edge or corner
/// \param knownAxesNum (input): the number of fixed axes where the boundary data is known
/// \param p (input): orderInSpace/2
/// \param dim (input): the dimension in space
void getInteriorIndRange(int interiorRange[3][2], int fixedAxis[], int fixedSide[], int knownAxesNum, int p, int dim){
    if(knownAxesNum == 1){
        for(int axis = 0; axis<3; axis++){
            if(axis<dim){
                interiorRange[axis][0] = 0;
                interiorRange[axis][1] = (2*p);
            }else{
                interiorRange[axis][0] = 0;
                interiorRange[axis][1] = 0;
            }
        }
        interiorRange[fixedAxis[0]][0] = sideBasedValue(fixedSide[0], p, 0); // if side = 0, choose p else choose 0
        interiorRange[fixedAxis[0]][1] = sideBasedValue(fixedSide[0], (2*p), p);
    }
    else if(knownAxesNum == 2){
        
        interiorRange[0][0] = 0;
        interiorRange[0][1] = (2*p);
        interiorRange[1][0] = 0;
        interiorRange[1][1] = (2*p);
        interiorRange[2][0] = 0;
        interiorRange[2][1] = dimBasedValue(dim,0,(2*p));
        
        interiorRange[fixedAxis[0]][0] = sideBasedValue(fixedSide[0], p, 0);
        interiorRange[fixedAxis[0]][1] = sideBasedValue(fixedSide[0], (2*p), p);
        interiorRange[fixedAxis[1]][0] = sideBasedValue(fixedSide[1], p, 0);
        interiorRange[fixedAxis[1]][1] = sideBasedValue(fixedSide[1], (2*p), p);
    }
    else{
        for(int axis = 0; axis<3; axis++){
            interiorRange[axis][0] = sideBasedValue(fixedSide[axis], p, 0);
            interiorRange[axis][1] = sideBasedValue(fixedSide[axis], (2*p), p);
        }
    }
}

/// \brief Given a two fixed axes, this function finds and returns the number of the varying axis
/// \param fixedAxis (input): a vector holding the numbers of the fixed axes
int getVarAxis(int fixedAxis[2]){
    int varAxis = -1;
    for(int axis = 0; axis<3; axis++){
        if((axis != fixedAxis[0]) && (axis != fixedAxis[1]))
        {
            varAxis = axis;
            break;
        }
    }
    return varAxis;
}

/// \brief A function that finds the boundary range of indices at a given boundary (face, edge or vertex)
/// \param bdryRange (input/output): carries the range of indices at the boundary
/// \param LcbcBdryRange (input): the boundary range used in the LCBC procedure (may be smaller that the actual boundary range at a face)
/// \param faceBdryRange (input): the actual boundary range at a face
/// \param bdryNgx (input/output): the number of grid points at a given boundary (face, edge or vertex)
/// \param faceBdryNgx (input/output): the number of grid points at a given face (where boundary data is known)
void getBdryRange(int bdryRange[3][2], int LcbcBdryRange[3][2], int faceBdryRange[3][2], int bdryNgx[3], int faceBdryNgx[3], int fixedAxis[], int fixedSide[], int fixedAxesNum){
    if(fixedAxesNum==1){
        for(int ax = 0; ax<3; ax++){
            for(int si = 0; si<2; si++){
                bdryRange[ax][si] = LcbcBdryRange[ax][si];
            }// end of si loop
            bdryNgx[ax] = faceBdryNgx[ax];
        }// end of ax loop
    }
    else if(fixedAxesNum == 2){
        for(int side = 0; side<2; side++){
            for(int fixedAxInd = 0; fixedAxInd < 2; fixedAxInd++){
                bdryRange[fixedAxis[fixedAxInd]][side] = faceBdryRange[fixedAxis[fixedAxInd]][fixedSide[fixedAxInd]];
            }
            int varAxis = getVarAxis(fixedAxis);
            bdryRange[varAxis][side] = LcbcBdryRange[varAxis][side];
            bdryNgx[varAxis] = faceBdryNgx[varAxis];
        }
    }
    else if(fixedAxesNum == 3){
        for(int side = 0; side<2; side++){
            for(int fixedAxInd = 0; fixedAxInd < 3; fixedAxInd++){
                bdryRange[fixedAxis[fixedAxInd]][side] = faceBdryRange[fixedAxis[fixedAxInd]][fixedSide[fixedAxInd]];
            }
        }
    }
}
