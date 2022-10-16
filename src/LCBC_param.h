#ifndef LCBC_param_h
#define LCBC_param_h
#include "LCBCmacros.h"

#define pInd(k,deriv) (k + deriv*Param::orderNum)

class Param{
public:
    int faceNum, NU, totalVarNum, interiorEqNum, unknownVarNum, auxiliaryEqNum;
    bool exists = false;
    static const int orderNum = 5; // orders = 0, 2, 4, 6
    static const int derivNum = 9; // derivatives = 0, 1, ..., 6
    static const int maxStencilWth = 9;
    double derivCoef[(orderNum*derivNum)][maxStencilWth]; int derivRange[(orderNum*derivNum)];
    
    Param(){
    }
    
    void initialize(int dim, int p){
        if(!exists){
            faceNum = 2*dim;
            NU = p+1;
            
            int n = (2*p+1);
            totalVarNum    = n*n*dimBasedValue(dim, 1, n);
            interiorEqNum  = (p+1)*n*dimBasedValue(dim, 1, n);
            unknownVarNum  = (totalVarNum - interiorEqNum);
            auxiliaryEqNum = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
            getDerivCoef(p);
            exists = true;
            
        }else{
            printf("LCBC Param object already exists\n");
        }
    }
    void getDerivCoef(int p){
        int k = 0, deriv;
        for(deriv = 0; deriv<=6; deriv++){
            derivRange[pInd(k,deriv)] = 0;
            derivCoef[pInd(k,deriv)][0] = 1.0;
        }// end of deriv loop
        
        k = 1;
        deriv = 0; derivRange[pInd(k,deriv)] = 0;
        derivCoef[pInd(k,deriv)][0] = 1.0;
        deriv = 1; derivRange[pInd(k,deriv)] = 1;
        derivCoef[pInd(k,deriv)][0] = -0.5;
        derivCoef[pInd(k,deriv)][1] =  0.0;
        derivCoef[pInd(k,deriv)][2] =  0.5;
        deriv = 2; derivRange[pInd(k,deriv)] = 1;
        derivCoef[pInd(k,deriv)][0] =  1.0;
        derivCoef[pInd(k,deriv)][1] = -2.0;
        derivCoef[pInd(k,deriv)][2] =  1.0;
        if(p > 1){
            k = 1;
            deriv = 3; derivRange[pInd(k,deriv)] = 2;
            derivCoef[pInd(k,deriv)][0] = -0.5;
            derivCoef[pInd(k,deriv)][1] =  1.0;
            derivCoef[pInd(k,deriv)][2] =  0.0;
            derivCoef[pInd(k,deriv)][3] = -1.0;
            derivCoef[pInd(k,deriv)][4] =  0.5;
            deriv = 4; derivRange[pInd(k,deriv)] = 2;
            derivCoef[pInd(k,deriv)][0] =  1.0;
            derivCoef[pInd(k,deriv)][1] = -4.0;
            derivCoef[pInd(k,deriv)][2] =  6.0;
            derivCoef[pInd(k,deriv)][3] = -4.0;
            derivCoef[pInd(k,deriv)][4] =  1.0;
            
            k = 2;
            deriv = 0; derivRange[pInd(k,deriv)] = 0;
            derivCoef[pInd(k,deriv)][0] =  1.0;
            deriv = 1; derivRange[pInd(k,deriv)] = 2;
            derivCoef[pInd(k,deriv)][0] =  1.0/12.0;
            derivCoef[pInd(k,deriv)][1] = -2.0/3.0;
            derivCoef[pInd(k,deriv)][2] =  0.0;
            derivCoef[pInd(k,deriv)][3] =  2.0/3.0;
            derivCoef[pInd(k,deriv)][4] = -1.0/12.0;
            deriv = 2; derivRange[pInd(k,deriv)] = 2;
            derivCoef[pInd(k,deriv)][0] = -1.0/12.0;
            derivCoef[pInd(k,deriv)][1] =  4.0/3.0;
            derivCoef[pInd(k,deriv)][2] = -5.0/2.0;
            derivCoef[pInd(k,deriv)][3] =  4.0/3.0;
            derivCoef[pInd(k,deriv)][4] = -1.0/12.0;
        }
        if(p > 2){
            k = 1;
            deriv = 5; derivRange[pInd(k,deriv)] = 3;
            derivCoef[pInd(k,deriv)][0] = -0.5;
            derivCoef[pInd(k,deriv)][1] =  2.0;
            derivCoef[pInd(k,deriv)][2] = -2.5;
            derivCoef[pInd(k,deriv)][3] =  0.0;
            derivCoef[pInd(k,deriv)][4] =  2.5;
            derivCoef[pInd(k,deriv)][5] = -2.0;
            derivCoef[pInd(k,deriv)][6] =  0.5;
            deriv = 6; derivRange[pInd(k,deriv)] = 3;
            derivCoef[pInd(k,deriv)][0] =   1.0;
            derivCoef[pInd(k,deriv)][1] =  -6.0;
            derivCoef[pInd(k,deriv)][2] =  15.0;
            derivCoef[pInd(k,deriv)][3] = -20.0;
            derivCoef[pInd(k,deriv)][4] =  15.0;
            derivCoef[pInd(k,deriv)][5] =  -6.0;
            derivCoef[pInd(k,deriv)][6] =   1.0;
            
            k = 2;
            deriv = 3; derivRange[pInd(k,deriv)] = 3;
            derivCoef[pInd(k,deriv)][0] =   1.0/8.0;
            derivCoef[pInd(k,deriv)][1] =  -1.0;
            derivCoef[pInd(k,deriv)][2] =  13.0/8.0;
            derivCoef[pInd(k,deriv)][3] =   0.0;
            derivCoef[pInd(k,deriv)][4] = -13.0/8.0;
            derivCoef[pInd(k,deriv)][5] =   1.0;
            derivCoef[pInd(k,deriv)][6] =  -1.0/8.0;
            deriv = 4; derivRange[pInd(k,deriv)] = 3;
            derivCoef[pInd(k,deriv)][0] =  -1.0/6.0;
            derivCoef[pInd(k,deriv)][1] =   2.0;
            derivCoef[pInd(k,deriv)][2] = -13.0/2.0;
            derivCoef[pInd(k,deriv)][3] =  28.0/3.0;
            derivCoef[pInd(k,deriv)][4] = -13.0/2.0;
            derivCoef[pInd(k,deriv)][5] =   2.0;
            derivCoef[pInd(k,deriv)][6] =  -1.0/6.0;
            
            k = 3;
            deriv = 0; derivRange[pInd(k,deriv)] = 0;
            derivCoef[pInd(k,deriv)][0] =  1.0;
            deriv = 1; derivRange[pInd(k,deriv)] = 3;
                derivCoef[pInd(k,deriv)][0] = -1.0/60.0;
                derivCoef[pInd(k,deriv)][1] =  3.0/20.0;
                derivCoef[pInd(k,deriv)][2] = -3.0/4.0;
                derivCoef[pInd(k,deriv)][3] =  0.0;
                derivCoef[pInd(k,deriv)][4] =  3.0/4.0;
                derivCoef[pInd(k,deriv)][5] = -3.0/20.0;
                derivCoef[pInd(k,deriv)][6] =  1.0/60.0;
            deriv = 2; derivRange[pInd(k,deriv)] = 3;
                derivCoef[pInd(k,deriv)][0] =   1.0/90.0;
                derivCoef[pInd(k,deriv)][1] =  -3.0/20.0;
                derivCoef[pInd(k,deriv)][2] =   3.0/2.0;
                derivCoef[pInd(k,deriv)][3] = -49.0/18.0;
                derivCoef[pInd(k,deriv)][4] =   3.0/2.0;
                derivCoef[pInd(k,deriv)][5] =  -3.0/20.0;
                derivCoef[pInd(k,deriv)][6] =   1.0/90.0;
        }
        if(p > 3){
            k = 1;
            deriv = 7; derivRange[pInd(k,deriv)] = 4;
            derivCoef[pInd(k,deriv)][0] = -0.5;
            derivCoef[pInd(k,deriv)][1] =  3.0;
            derivCoef[pInd(k,deriv)][2] = -7.0;
            derivCoef[pInd(k,deriv)][3] =  7.0;
            derivCoef[pInd(k,deriv)][4] =  0.0;
            derivCoef[pInd(k,deriv)][5] = -7.0;
            derivCoef[pInd(k,deriv)][6] =  7.0;
            derivCoef[pInd(k,deriv)][7] = -3.0;
            derivCoef[pInd(k,deriv)][8] =  0.5;
            deriv = 8; derivRange[pInd(k,deriv)] = 4;
            derivCoef[pInd(k,deriv)][0] =   1.0;
            derivCoef[pInd(k,deriv)][1] =  -8.0;
            derivCoef[pInd(k,deriv)][2] =  28.0;
            derivCoef[pInd(k,deriv)][3] = -56.0;
            derivCoef[pInd(k,deriv)][4] =  70.0;
            derivCoef[pInd(k,deriv)][5] = -56.0;
            derivCoef[pInd(k,deriv)][6] =  28.0;
            derivCoef[pInd(k,deriv)][7] =  -8.0;
            derivCoef[pInd(k,deriv)][8] =   1.0;

            k = 2;
            deriv = 5; derivRange[pInd(k,deriv)] = 4;
            derivCoef[pInd(k,deriv)][0] =   1.0/6.0;
            derivCoef[pInd(k,deriv)][1] =  -3.0/2.0;
            derivCoef[pInd(k,deriv)][2] =  13.0/3.0;
            derivCoef[pInd(k,deriv)][3] = -29.0/6.0;
            derivCoef[pInd(k,deriv)][4] =   0.0;
            derivCoef[pInd(k,deriv)][5] =  29.0/6.0;
            derivCoef[pInd(k,deriv)][6] = -13.0/3.0;
            derivCoef[pInd(k,deriv)][7] =   3.0/2.0;
            derivCoef[pInd(k,deriv)][8] =  -1.0/6.0;
            deriv = 6; derivRange[pInd(k,deriv)] = 4;
            derivCoef[pInd(k,deriv)][0] = -1.0/4.0;
            derivCoef[pInd(k,deriv)][1] =  3.0;
            derivCoef[pInd(k,deriv)][2] = -13.0;
            derivCoef[pInd(k,deriv)][3] =  29.0;
            derivCoef[pInd(k,deriv)][4] = -75/2.0;
            derivCoef[pInd(k,deriv)][5] =  29.0;
            derivCoef[pInd(k,deriv)][6] = -13.0;
            derivCoef[pInd(k,deriv)][7] =   3.0;
            derivCoef[pInd(k,deriv)][8] = -1.0/4.0;

            k = 3;
            deriv = 3; derivRange[pInd(k,deriv)] = 4;
            derivCoef[pInd(k,deriv)][0] = -7.0/240.0;
            derivCoef[pInd(k,deriv)][1] =  3.0/10.0;
            derivCoef[pInd(k,deriv)][2] = -169.0/120.0;
            derivCoef[pInd(k,deriv)][3] =  61.0/30.0;
            derivCoef[pInd(k,deriv)][4] =  0.0;
            derivCoef[pInd(k,deriv)][5] = -61.0/30.0;
            derivCoef[pInd(k,deriv)][6] =  169.0/120.0;
            derivCoef[pInd(k,deriv)][7] = -3.0/10.0;
            derivCoef[pInd(k,deriv)][8] =  7.0/240.0;
            deriv = 4; derivRange[pInd(k,deriv)] = 4;
            derivCoef[pInd(k,deriv)][0] =  7.0/240.0;
            derivCoef[pInd(k,deriv)][1] = -2.0/5.0;
            derivCoef[pInd(k,deriv)][2] =  169.0/60.0;
            derivCoef[pInd(k,deriv)][3] = -122.0/15.0;
            derivCoef[pInd(k,deriv)][4] =  91.0/8.0;
            derivCoef[pInd(k,deriv)][5] = -122.0/15.0;
            derivCoef[pInd(k,deriv)][6] =  169.0/60.0;
            derivCoef[pInd(k,deriv)][7] = -2.0/5.0;
            derivCoef[pInd(k,deriv)][8] =  7.0/240.0;

            k = 4;
            deriv = 0; derivRange[pInd(k,deriv)] = 0;
            derivCoef[pInd(k,deriv)][0] =  1.0;
            deriv = 1; derivRange[pInd(k,deriv)] = 4;
            derivCoef[pInd(k,deriv)][0] =  1.0/280.0;
            derivCoef[pInd(k,deriv)][1] = -4.0/105.0;
            derivCoef[pInd(k,deriv)][2] =  1.0/5.0;
            derivCoef[pInd(k,deriv)][3] = -4.0/5.0;
            derivCoef[pInd(k,deriv)][4] =  0.0;
            derivCoef[pInd(k,deriv)][5] =  4.0/5.0;
            derivCoef[pInd(k,deriv)][6] = -1.0/5.0;
            derivCoef[pInd(k,deriv)][7] =  4.0/105.0;
            derivCoef[pInd(k,deriv)][8] = -1.0/280.0;
            deriv = 2; derivRange[pInd(k,deriv)] = 4;
            derivCoef[pInd(k,deriv)][0] = -1.0/560.0;
            derivCoef[pInd(k,deriv)][1] =  8.0/315.0;
            derivCoef[pInd(k,deriv)][2] = -1.0/5.0;
            derivCoef[pInd(k,deriv)][3] =  8.0/5.0;
            derivCoef[pInd(k,deriv)][4] = -205.0/72.0;
            derivCoef[pInd(k,deriv)][5] =  8.0/5.0;
            derivCoef[pInd(k,deriv)][6] = -1.0/5.0;
            derivCoef[pInd(k,deriv)][7] =  8.0/315.0;
            derivCoef[pInd(k,deriv)][8] = -1.0/560.0;
        }
    }
    
    ~Param(){}
};

class FaceParam{
public:
    int bdryRangeExt[3][2], LcbcBdryRange[3][2], bdryRange[3][2], bdryNgx[3], otherAxis[2], bdryRangeLth[3], interiorRange[3][2], bdryNg;
    int NU, compCondNum, approxEqNum, auxiliaryEqNum;
    bool exists = false;
    
    FaceParam(){}
    
    void initialize(int axis, int side, int indexRange[3][2], int Ngx[3], int p, int dim, int faceEval[]){
        if(!exists){
            int n = (2*p+1);
            int face = side + 2*axis;
            NU = (faceEval[face] == 1)?(p+1):(p);
            auxiliaryEqNum = (p*(2*p+1)*dimBasedValue(dim, 1, (10*p+7)))/dimBasedValue(dim, 1,3);
            compCondNum = NU*n*dimBasedValue(dim, 1, n);
            approxEqNum = compCondNum + auxiliaryEqNum;
            getBdryRange(indexRange, Ngx, axis, side, faceEval, dim, p);
            bdryNg = bdryNgx[0]*bdryNgx[1]*bdryNgx[2];
            getVarAxis(otherAxis, axis);
            getInteriorRange(axis, side, p, dim);
            exists = true;
        }else{
            printf("LCBC FaceParam object already exists\n");
        }
    }
    
    void getBdryRange(int indexRange[3][2], int Ngx[3], int fixedAxis, int fixedSide, int *faceEval, int dim, int p){
        for(int axis = 0; axis<3; axis++){
            for(int side = 0; side<2; side++){
                if(axis == fixedAxis){
                    bdryRange[axis][side]     = indexRange[fixedAxis][fixedSide];
                    LcbcBdryRange[axis][side] = indexRange[fixedAxis][fixedSide];
                    bdryRangeExt[axis][side]  = indexRange[fixedAxis][fixedSide];
                    bdryNgx[axis] = 1;
                }
                else{
                    if(axis<dim){
                        bdryRange[axis][side]     = indexRange[axis][side];
                        bdryRangeExt[axis][side]  = indexRange[axis][side] + (1-side)*(-p) + (p)*side;
                        
                        if(faceEval[(2*axis)]>=0){
                            LcbcBdryRange[axis][side] = indexRange[axis][side] + (1-side)*(p) + (-p)*side;
                        }else{
                            LcbcBdryRange[axis][side] = indexRange[axis][side];
                        }
                      
                    }else{
                        bdryRange[axis][side]     = indexRange[axis][side];
                        LcbcBdryRange[axis][side] = indexRange[axis][side];
                        bdryRangeExt[axis][side]  = indexRange[axis][side];
                    }
                    bdryNgx[axis] = Ngx[axis];
                }// end of if axis
            }// end of side
            bdryRangeLth[axis] = bdryRange[axis][1] - bdryRange[axis][0] + 1;
        }// end of axis
    }// end of getBdryRange
    
    void getVarAxis(int otherAxis[2], int fixedAxis){
        int cnt = 0;
        for(int d  = 0; d<3; d++){
            if(d!=fixedAxis){
                otherAxis[cnt] = d;
                cnt++;
            }// end of if d statement
        }// end of d for loop
    }// end of getVarAxis
    
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
