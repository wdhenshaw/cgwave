#ifndef numericalDeriv_h
#define numericalDeriv_h

#include <stdio.h>
#include <cstdlib>

class NumericalDeriv{
private:
    bool allocateSpaceFlag;
    int **coefVecLocator, coefVecCount; // coefVecLocator(range,derivative)
    double **coefVec; // coefVec(range,derivative)
    
    void centeredDifferenceWeights(int r, int d, double *c);
    double * getDerivCoef(int r,int d);
    void allocateSpace();
    void reallocateSpace(int newMaxRange);
    
public:
    int maxRange;
    void expandSpace(int maxRangeIn);
    double numDeriv(double *u, int narg, int index[3], int pos, int d, int order, double h, int *size);
    double numDerivFn(double (*F) (double *), double *arg, int pos, int d, int order, double h);
    double mixedNumDeriv(double *u, int narg, int index[3], int Deriv[3], int order[3], double h[3], int *size);
    
    
    NumericalDeriv(int maxRangeIn = 3){
        allocateSpaceFlag = false;
        coefVecCount = 1;
        maxRange = maxRangeIn;
    }
    
    ~NumericalDeriv(){
        if(allocateSpaceFlag){
        for(int r = 1; r<=maxRange; r++){
            for(int d = 1; d<=(2*r); d++ ){
                if(!(coefVecLocator[r][d] == 0)){
                    free(coefVec[(coefVecLocator[r][d]-1)]); coefVec[(coefVecLocator[r][d]-1)] = NULL;
                }
            }
            free(coefVecLocator[r]); coefVecLocator[r] = NULL;
        }
        free(coefVecLocator); coefVecLocator = NULL;
        free(coefVec); coefVec = NULL;
        }
    }
};

#endif /* numericalDeriv_hpp */
