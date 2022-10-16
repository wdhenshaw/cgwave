#include "LCBC_data.h"
#include "LCBCmacros.h"
#include <string.h>

void equateVectors(int v1[3], int v2[3]){
    for(int i = 0; i<3; i++){
        v1[i] = v2[i];
    }
}
/* Start here */

void LcbcData::initialize( int lthIn[3], int wthIn[3], int indexLthIn, double **FnIn, bool fixIn){
    
    if(!initialized){
        equateVectors(lth, lthIn);
        equateVectors(wth, wthIn);
        indexLth = indexLthIn;
        
        if(FnIn == NULL){
            allocateSpace(lth, indexLth);
        }else{
            Fn = FnIn;
        }
        
        fix = fixIn;
        initialized = true;
    }
}// end of initialize

void LcbcData::allocateSpace(int lth[3], int indexLth){
    if(!allocatedSpace){ 
        Fn = new double*[indexLth];
        for(int i = 0; i<indexLth; i++){
            Fn[i] = new double[(lth[0]*lth[1]*lth[2])];
        }
        allocatedSpace = true;
    }
}

void LcbcData::deleteLcbcData(){
    if(allocatedSpace){
        for(int i = 0; i<indexLth; i++){
            delete [] Fn[i];
        }
        delete [] Fn;
    }
}// end of deleteLcbcData

/* Check Data functions */

enum lcbcFnType{boundary=1, forcing=2, coefficient=3};
enum dataInput{less=-1, more=1, sufficient=0, insufficient=999};

bool debug = false;

void checkUserData(LcbcData *&classData, LcbcData *userData, int Nx[3], int Ngx[3], int faceEval[], int dim, int p, int extraDataGhost, int dataType, bool fix, bool overRide){
    
    int axis = -1;
    for(int face = 0; face<(2*dim); face++){
        if((face%2) == 0){
            axis++;
        }
        if(faceEval[face]>0){
            checkData(classData[face], userData[face], Nx, Ngx, axis, p, dim, extraDataGhost, dataType, overRide);
            if(fix){
                if(classData[face].fill){
                    fillData(classData[face], userData[face]);
                }
                if(classData[face].fix){
                    fixData(classData[face], p, dim, axis);
                }
            }
        }
    }// end of face loop
}// end of checkUserData

void checkData(LcbcData &classData, LcbcData &userData, int Nx[3], int Ngx[3], int axis, int p, int dim, int extraDataGhost, int FnType, bool overRide){
    
    int lcbcLth[3], lcbcWth[3], minLth[3];
    
    if(overRide){
        classData = userData;
    }
    else{
        if(FnType == boundary){
            getBdryLcbcLthWth(lcbcLth, lcbcWth, minLth, Ngx, axis, p, dim);
        }else if(FnType == forcing || FnType == coefficient){
            getForcingLthWth(Nx, lcbcLth, lcbcWth, minLth, axis, dim, p, extraDataGhost);
        }else{
            printf("LCBC_data: ERROR data function type is incorrect\nUse 1:boudary, 2:forcing, 3:coefficient\n");
            exit(-1);
        }
        
        int userInput = compareVectors(userData.lth, lcbcLth, minLth);
        if(userInput == less){
            if(debug)
                printf("userInput less\n");
            
            classData.initialize(lcbcLth, lcbcWth, userData.indexLth);
            classData.fill = true;
            classData.fix = true;
        }
        else if((userInput == more) || (userInput == sufficient)){
            if(debug)
                printf("userInput is sufficient\n");
            
            classData.initialize(userData.lth, userData.wth, userData.indexLth, userData.Fn, userData.fix);
        }
        else if(userInput == insufficient){
            if(debug)
                printf("userInput insufficient\n");
            if(FnType == boundary){
                printf("LCBC_data: ERROR in given boundary data.\nneed at least orderInSpace/2 ghost points in tangential directions\n");
                exit(-1);
            }
            else if(FnType == forcing || FnType == coefficient){
                printf("LCBC_data: ERROR in given forcing and coefficient data. Need at least orderInSpace/2 ghost points in each axis\n");
                exit(-1);
            }
        }
        else{
            printf("LCBC_data: ERROR - unidentified user input type\n");
            exit(-1);
        }
    }// end of else to overRide
}// end of checkData

void fixData(LcbcData &data, int p, int dim, int axis){
    
    if(debug)
        printf("fixing Data\n");
    
        int varAxis[2]; getVarAxis_data(varAxis, axis);
        int v0 = varAxis[0], v1 = varAxis[1];
        int i[3];
        
        for(int nu = 0; nu<data.indexLth; nu++){
            /* LOWER y axis over partial z axis */
            for(i[v0] = (data.wth[v0] - p-1); i[v0]>=0; i[v0]--){
                for(i[v1] = dimBasedValue(dim, 0, (data.wth[v1] - p)); i[v1]<dimBasedValue(dim, 1, (data.lth[v1] - data.wth[v1]+p)); i[v1]++){
                    for(i[axis] = 0; i[axis]<data.lth[axis]; i[axis]++){
                        
                        data.Fn[nu][ind(i,data.lth)] = extrapolate_data(data.Fn[nu], i, data.lth, v0, 1, (2*p));
                        
                    }// end of i[0] loop
                }// end of i[1] loop
            }// end of i[2] loop
            
            /* Upper y axis over partial z axis */
            for(i[v0] = (data.lth[v0]-data.wth[v0] + p); i[v0]<data.lth[v0]; i[v0]++){
                for(i[v1] = dimBasedValue(dim, 0, (data.wth[v1] - p)); i[v1]<dimBasedValue(dim, 1, (data.lth[v1] - data.wth[v1] + p)); i[v1]++){
                    for(i[axis] = 0; i[axis]<data.lth[axis]; i[axis]++){
                        
                        data.Fn[nu][ind(i,data.lth)] = extrapolate_data(data.Fn[nu], i, data.lth, v0, (-1), (2*p));
                        
                    }// end of i[0] loop
                }// end of i[1] loop
            }// end of i[2] loop
            
            if(dim == 3){
                /* Lower z */
                for(i[v1] = (data.wth[v1] - p-1); i[v1]>=0; i[v1]--){
                    for(i[v0] = 0; i[v0]<data.lth[v0]; i[v0]++){
                        for(i[axis] = 0; i[axis]<data.lth[axis]; i[axis]++){
                            
                            data.Fn[nu][ind(i,data.lth)] = extrapolate_data(data.Fn[nu], i, data.lth, v1, 1, (2*p));
                            
                        }// end of i[axis]
                    }// end of i[v0]
                }// end of i[v1]
                
                /* Upper z */
                for(i[v1] = (data.lth[v1]-data.wth[v1]+p); i[v1]< data.lth[v1]; i[v1]++){
                    for(i[v0] = 0; i[v0]<data.lth[v0]; i[v0]++){
                        for(i[axis] = 0; i[axis]<data.lth[axis]; i[axis]++){
                            
                            data.Fn[nu][ind(i,data.lth)] = extrapolate_data(data.Fn[nu], i, data.lth, v1, (-1), (2*p));
                            
                        }// end of i[axis]
                    }// end of i[v0]
                }// end of i[v1]
            }// end of if dim statment
            
        }// end of nu loop
}// end of fixData

void fillData(LcbcData &classData, LcbcData &userData){
    
    if(debug)
        printf("filling Data\n");
    
    int i[3], j[3];
    for(int nu = 0; nu<classData.indexLth; nu++){
        for(i[2]=0; i[2]<userData.lth[2]; i[2]++){
            for(i[1]=0; i[1]<userData.lth[1]; i[1]++){
                for(i[0]=0; i[0]<userData.lth[0]; i[0]++){
                    j[0] = i[0] + (classData.wth[0] - userData.wth[0]);
                    j[1] = i[1] + (classData.wth[1] - userData.wth[1]);
                    j[2] = i[2] + (classData.wth[2] - userData.wth[2]);
                    classData.Fn[nu][ind(j,classData.lth)] = userData.Fn[nu][ind(i,userData.lth)];
                }// end of i[0]
            }// end of i[1]
        }// end of i[2]
    }// end of nu
}// end of fill

void getBdryLcbcLthWth(int (&lcbcLth)[3], int (&lcbcWth)[3], int (&minLth)[3], int Ngx[3], int fixedAxis, int p, int dim){
    for(int axis = 0; axis<3; axis++){
        lcbcLth[axis] = (axis == fixedAxis)?(1):Ngx[axis];
        minLth[axis] = lcbcLth[axis];
        if(axis<dim){
            lcbcWth[axis] = (axis== fixedAxis)?(0):(p);
        }else{
            lcbcWth[axis] = 0;
        }
    }// end of axis
}

void getForcingLthWth(int Nx[3], int (&lcbcLth)[3], int (&lcbcWth)[3], int (&minLth)[3], int fixedAxis, int dim, int p, int extraDataGhost){
    
    for(int axis = 0; axis<3; axis++){
        if(axis<dim){
            lcbcWth[axis] = ((axis==fixedAxis)?(p):(p + extraDataGhost));
            int bdryLth = ((axis==fixedAxis)?(1):(Nx[axis]+1));
            lcbcLth[axis] = (2*lcbcWth[axis] + bdryLth);
            minLth[axis] = (2*p + bdryLth);
        }else{
            lcbcWth[axis] = 0;
            lcbcLth[axis] = 1;
            minLth[axis]  = 1;
        }
    }// end of axis
}

void getVarAxis_data(int varAxis[2], int fixedAxis){
    int cnt = 0;
    for(int d  = 0; d<3; d++){
        if(d!=fixedAxis){
            varAxis[cnt] = d;
            cnt++;
        }// end of if d statement
    }// end of d for loop
}// end of getVarAxis_data

int compareVectors(int vec[3], int cstVec[3], int minVec[3]){
    for(int i = 0; i<3; i++){
        if(vec[i] != cstVec[i]){
            if(vec[i]>cstVec[i]){
                return 1;
            }
            else if(vec[i]< minVec[i]){
                return 999;
            }
            else if((minVec[i]<=vec[i]) && (vec[i]<cstVec[i])){
                return (-1);
            }
        }// end of if not equal
    }// end of for loop
    return 0;
}// end of compareVectors

double extrapolate_data(double *gridFn, int index[3], int lth[3], int axis, int sign, int order){
    
    double E = 0;
    
    if(order == 2){
        int i1[3]; memcpy(i1, index, sizeof(i1)); i1[axis] = i1[axis] + sign*1;
        int i2[3]; memcpy(i2, index, sizeof(i1)); i2[axis] = i2[axis] + sign*2;
        int i3[3]; memcpy(i3, index, sizeof(i1)); i3[axis] = i3[axis] + sign*3;
        
        E = 3.0*gridFn[ind(i1,lth)]
           -3.0*gridFn[ind(i2,lth)]
           +1.0*gridFn[ind(i3,lth)];
    }
    
    else if(order == 4){
        
        int i1[3]; memcpy(i1, index, sizeof(i1)); i1[axis] = i1[axis] + sign*1;
        int i2[3]; memcpy(i2, index, sizeof(i1)); i2[axis] = i2[axis] + sign*2;
        int i3[3]; memcpy(i3, index, sizeof(i1)); i3[axis] = i3[axis] + sign*3;
        int i4[3]; memcpy(i4, index, sizeof(i1)); i4[axis] = i4[axis] + sign*4;
        int i5[3]; memcpy(i5, index, sizeof(i1)); i5[axis] = i5[axis] + sign*5;
        
        E =  5.0*gridFn[ind(i1,lth)]
           -10.0*gridFn[ind(i2,lth)]
           +10.0*gridFn[ind(i3,lth)]
            -5.0*gridFn[ind(i4,lth)]
            +1.0*gridFn[ind(i5,lth)];
    }
    
    else if(order == 6){
        int i1[3]; memcpy(i1, index, sizeof(i1)); i1[axis] = i1[axis] + sign*1;
        int i2[3]; memcpy(i2, index, sizeof(i1)); i2[axis] = i2[axis] + sign*2;
        int i3[3]; memcpy(i3, index, sizeof(i1)); i3[axis] = i3[axis] + sign*3;
        int i4[3]; memcpy(i4, index, sizeof(i1)); i4[axis] = i4[axis] + sign*4;
        int i5[3]; memcpy(i5, index, sizeof(i1)); i5[axis] = i5[axis] + sign*5;
        int i6[3]; memcpy(i6, index, sizeof(i1)); i6[axis] = i6[axis] + sign*6;
        int i7[3]; memcpy(i7, index, sizeof(i1)); i7[axis] = i7[axis] + sign*7;
        
        E = 7.0*gridFn[ind(i1,lth)]
        -21.0*gridFn[ind(i2,lth)]
        +35.0*gridFn[ind(i3,lth)]
        -35.0*gridFn[ind(i4,lth)]
        +21.0*gridFn[ind(i5,lth)]
        -7.0*gridFn[ind(i6,lth)]
        +1.0*gridFn[ind(i7,lth)];
        
    }else if(order == 8){
            int i1[3]; memcpy(i1, index, sizeof(i1)); i1[axis] = i1[axis] + sign*1;
            int i2[3]; memcpy(i2, index, sizeof(i1)); i2[axis] = i2[axis] + sign*2;
            int i3[3]; memcpy(i3, index, sizeof(i1)); i3[axis] = i3[axis] + sign*3;
            int i4[3]; memcpy(i4, index, sizeof(i1)); i4[axis] = i4[axis] + sign*4;
            int i5[3]; memcpy(i5, index, sizeof(i1)); i5[axis] = i5[axis] + sign*5;
            int i6[3]; memcpy(i6, index, sizeof(i1)); i6[axis] = i6[axis] + sign*6;
            int i7[3]; memcpy(i7, index, sizeof(i1)); i7[axis] = i7[axis] + sign*7;
            int i8[3]; memcpy(i8, index, sizeof(i1)); i8[axis] = i8[axis] + sign*8;
            int i9[3]; memcpy(i9, index, sizeof(i1)); i9[axis] = i9[axis] + sign*9;
            
            E = 9.0*gridFn[ind(i1,lth)]
              -36.0*gridFn[ind(i2,lth)]
              +84.0*gridFn[ind(i3,lth)]
             -126.0*gridFn[ind(i4,lth)]
             +126.0*gridFn[ind(i5,lth)]
              -84.0*gridFn[ind(i6,lth)]
              +36.0*gridFn[ind(i7,lth)]
               -9.0*gridFn[ind(i8,lth)]
               +1.0*gridFn[ind(i9,lth)];
        
    }else{
        printf("Error in LCBC_data: extrapolation of data functions (coefficients and forcing) to order %d is not supported (extrapolation/utility.C)\n",order);
        exit(-1);
    }
    
    return E;
}


