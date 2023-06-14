#ifndef LCBC_param_h
#define LCBC_param_h
#include "LCBCmacros.h"

/// \brief In this file, we define bdryTypeParam and FaceParam.
/// bdryTypeParam carries parameters to identify the type of boundary for which the LCBC procedure will be applied
/// FaceParam class carries parameters specific to a face

/// \brief bdryTypeParam is an object which carries parameters to identify the type of boundary for which the LCBC procedure will be applied.
/// For example, an edge formed from two physical faces is one type of boundary evaluation while an edge formed from one physical face and one interpolation face is another.
/// We identify the boundary type used the number of fixed axes and known axes (fixed axes for which boundary data is known)
struct bdryTypeParam{
    int interiorEqNum;
    int unknownVarNum;
    int ghostPointNum;
    int fixedAxesNum;
    int knownAxesNum;
};

class FaceParam{
public:
    /// bdryRangeExt: the range of the boundary indices extended by p points
    /// LcbcBdryRange: the range of the boundary indices that are p points away from corners and edges
    /// bdryRange: the range of boundary indices
    /// maskBdryRange: the range of the boundary over which the mask is positive
    /// solnBdryRange: this is the boundary range for which the solution will be set to exact added for debugging purposes
    /// otherAxis: a vector carrying numbers corresponding to axes other than the fixed axis on the specific face
    /// bdryRangeLth: the length of the bdryRange object
    /// interiorRange: the range of interior indices
    /// bdryNg: the total number of grid points on the boundary
    /// NU: the number of primary CBCs specific to a face (depends on BC type, Dirichlet or Neumann)
    /// compCondNum: the number of equations in the LCBC procedure coming from CBCs and their tangential derivatives
    /// approxEqNum: the number of equations to be approximated in the LCBC procedure
    /// auxiliaryEqNum: the number of auxiliary equation in the LCBC procedure if needed. These are equations that set the high-order derivatives of the interpolating polynomial to zero - the derivatives that are not needed for accuracy
    int bdryRangeExt[3][2], LcbcBdryRange[3][2], bdryRange[3][2], bdryNgx[3], otherAxis[2], bdryRangeLth[3], interiorRange[3][2], maskBdryRange[3][2], bdryNg;
    int solnBdryRange[3][2];
    int NU, compCondNum, approxEqNum, auxiliaryEqNum;
    bool exists = false;
    
    FaceParam(){}
    
    /// \brief Initialize a FaceParam object that carries parameters needed for the LCBC procedure at a specific face
    /// \param axis (input)
    /// \param side (input)
    /// \param indexRange (input): range of indices on the grid
    /// \param Ngx (input): the total number of grid points in each dimension (including ghost points)
    /// \param p (input): order in space/2
    /// \param dim (input): dimension
    /// \param faceEval (input): a vector carrying numbers corresponding to the BC type at each boundary face (-1: periodic, 0: none, 1: Dirichlet, 2: Neumann)
    /// \param addAuxEqns (input): a boolean to tell whether auxiliary equations need to be added or not
    void initialize(int axis, int side, int *mask, int numGhost, int userNumGhost, int indexRange[3][2], int Ngx[3], int p, int dim, int faceEval[], bool addAuxEqns){
        if(!exists){
            int n = (2*p+1);
            int face = side + 2*axis;
            NU = (faceEval[face] == 1)?(p+1):(p);
            
            if(addAuxEqns)
                auxiliaryEqNum = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
            else
                auxiliaryEqNum = 0;
            
            compCondNum = NU*n*dimBasedValue(dim, 1, n);
            approxEqNum = compCondNum + auxiliaryEqNum;
            getVarAxis(otherAxis, axis);
            getBdryRange(indexRange, Ngx, axis, side, faceEval, dim, p);
            getMaskBdryRange(mask, Ngx, numGhost, userNumGhost, axis, side, dim);
            getLcbcBdryRange(axis, side, faceEval, dim, p);
            bdryNg = bdryNgx[0]*bdryNgx[1]*bdryNgx[2];
            getInteriorRange(axis, side, p, dim);
            exists = true;
        }else{
            printf("LCBC FaceParam object already exists\n");
        }
    }
    
    /// \brief Get the range of indices in each dimension at a boundary face: bdryRange, bdryRangeExt and bdryNgx (all defined within the FaceParam class above)
    /// \param indexRange (input): range of indices
    /// \param Ngx (input): the total number of grid points (including ghost) in each dimension
    /// \param fixedAxis (input): fixed axis at the given boundary face
    /// \param fixedSide (input): fixed side at the given boundary face
    /// \param faceEval (input): a vector carrying numbers corresponding to the BC type at each boundary face (-1: periodic, 0: none, 1: Dirichlet, 2: Neumann)
    /// \param dim (input): dimension
    /// \param p (input): order in space/2
    void getBdryRange(int indexRange[3][2], int Ngx[3], int fixedAxis, int fixedSide, int *faceEval, int dim, int p){
        for(int axis = 0; axis<3; axis++){
            for(int side = 0; side<2; side++){
                
                if(axis == fixedAxis){
                    bdryRange[axis][side]     = indexRange[fixedAxis][fixedSide];
                    bdryRangeExt[axis][side]  = indexRange[fixedAxis][fixedSide];
                    bdryNgx[axis] = 1;
                }
                else{
                    if(axis<dim){
                        bdryRange[axis][side]     = indexRange[axis][side];
                        bdryRangeExt[axis][side]  = indexRange[axis][side] + (1-side)*(-p) + (p)*side;
        
                    }else{
                        bdryRange[axis][side]     = indexRange[axis][side];
                        bdryRangeExt[axis][side]  = indexRange[axis][side];
                    }
                    bdryNgx[axis] = Ngx[axis];
                }// end of if axis
            }// end of side
            bdryRangeLth[axis] = bdryRange[axis][1] - bdryRange[axis][0] + 1;
        }// end of axis
    }// end of getBdryRange
    
    /// \brief Get the range of indices in each dimension at a boundary face: LcbcBdryRange and solnBdryRange (defined above)
    /// \param fixedAxis (input): fixed axis at the given boundary face
    /// \param fixedSide (input): fixed side at the given boundary face
    /// \param faceEval (input): a vector carrying numbers corresponding to the BC type at each boundary face (-1: periodic, 0: none, 1: Dirichlet, 2: Neumann)
    /// \param dim (input): dimension
    /// \param p (input): order in space/2
    void getLcbcBdryRange(int fixedAxis, int fixedSide, int *faceEval, int dim, int p){
        
        for(int a = 0; a<3; a++){
            for(int s = 0; s<2; s++){
                LcbcBdryRange[a][s] = bdryRange[a][s];
                solnBdryRange[a][s] = bdryRange[a][s];
            }
        }
        
        for(int axis = 0; axis<3; axis++){
            for(int side = 0; side<2; side++){
                int face = side + 2*axis;
                
                if(axis != fixedAxis){
                    if(axis<dim){
                        if(faceEval[face]>0){ // adjacent faces are Dirichlet or Neumann
                                LcbcBdryRange[axis][side] = LcbcBdryRange[axis][side] + (1-side)*(p) + (-p)*side;
                        }
                        else if (faceEval[face]==0){ // adjacent faces are zero-type boundary
                            LcbcBdryRange[axis][side] = maskBdryRange[axis][side];
                            solnBdryRange[axis][side] = solnBdryRange[axis][side] + (1-side)*(-p) + (p)*side;
                        }
                    }
                }// end of if axis
            }// end of side
        }// end of axis
    }// end of getBdryRange
    
    /// \brief Get the boundary range overwhich the mask is positive
    /// \param mask (input): a label to identify the type of grid point
    /// \param Ngx (input): a vector containing the number of grid points in each dimension
    /// \param numGhost (input): the number of ghost points used which is orderInSpace/2
    /// \param userNumGhost (input): the number of ghost points utilized by the user
    /// \param axis (input): the fixed axis at the specific face
    /// \param side (input): the fixed side at the specific face
    /// \param dim (input): dimension
    void getMaskBdryRange(int *mask, int Ngx[3], int numGhost, int userNumGhost, int axis, int side, int dim){
        int i[3];
        
        for(int a = 0; a<3; a++){
            for(int s = 0; s<2; s++){
                maskBdryRange[a][s] = bdryRange[a][s];
            }
        }
        
        i[axis] = bdryRange[axis][side];
        
        for(int vaxis = 0; vaxis<2; vaxis++){
            
            int posCount = 0;
            int negCount = 0;
            
            int v0 = otherAxis[vaxis];
            int v1 = otherAxis[(1-vaxis)];
            
            i[v1]   = bdryRange[v1][0];
            for(i[v0] = bdryRange[v0][0]; i[v0] <= bdryRange[v0][1]; i[v0]++){
                if(mask[solInd(i,Ngx)]>0 && posCount == 0){
                    maskBdryRange[v0][0] = i[v0];
                    posCount++;
                }
                if(mask[solInd(i,Ngx)]<0 && negCount == 0 && posCount>0){
                    maskBdryRange[v0][1] = i[v0]-1;
                    negCount++;
                }
                if((negCount>0) && (posCount>0))
                    break;
            }
        }
    }

    /// \brief Find the other axes given a fixed axis
    /// \param otherAxis (output)
    /// \param fixedAxis (input)
    void getVarAxis(int otherAxis[2], int fixedAxis){
        int cnt = 0;
        for(int d  = 0; d<3; d++){
            if(d!=fixedAxis){
                otherAxis[cnt] = d;
                cnt++;
            }// end of if d statement
        }// end of d for loop
    }// end of getVarAxis
    
    /// \brief Get the range of indices corresponding to grid points on the interior of the grid
    /// \param fixedAxis (input)
    /// \param fixedSide (input)
    /// \param p (input): order in space /2
    /// \param dim (input): dimension
    void getInteriorRange(int fixedAxis, int fixedSide, int p, int dim){
        for(int axis = 0; axis<3; axis++){
            if(axis<dim){
                interiorRange[axis][0] = 0;
                interiorRange[axis][1] = (2*p);
            }else{
                interiorRange[axis][0] = 0;
                interiorRange[axis][1] = 0;
            }
        }
        interiorRange[fixedAxis][0] = sideBasedValue(fixedSide, p, 0); // if side = 0, choose p else choose 0
        interiorRange[fixedAxis][1] = sideBasedValue(fixedSide, (2*p), p);
    }
    
    ~FaceParam(){}
    
};

#endif /* LCBC_param_h */
