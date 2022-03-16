#include "numericalDeriv.h"
#include <math.h>
#include <string.h>
#include <cstring>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define RANGE(d,ord) ((d==0 || ord == 0) ? (0) : (floor(((double)(d+1))/2.0)-1+(((double) ord)/2.0)))
#define ind(I,N) (I[0] + N[0]*(I[1]+ I[2]*N[1]))
#define sumVectors(v1,v2) {(v1[0]+v2[0]),(v1[1]+v2[1]),(v1[2]+v2[2])}
#define setVecTo(v) {v[0],v[1],v[2]}

void NumericalDeriv::expandSpace(int maxRangeIn){
    if(allocateSpaceFlag == false){
        maxRange = maxRangeIn;
        allocateSpace();
    }
    else{
        reallocateSpace(maxRangeIn);
        maxRange = maxRangeIn;
    }
}

void NumericalDeriv::allocateSpace(){
    coefVec     = (double**) malloc(0*sizeof(double*));
    coefVecLocator   = (int **) malloc((maxRange+1)*sizeof(int *));
    for(int r = 1; r<=maxRange; r++){
        coefVecLocator[r] = (int *) malloc((2*r+1)*sizeof(int));
        for(int d = 1; d<=(2*r); d++ ){
            coefVecLocator[r][d] = 0;
        }
    }
}
void NumericalDeriv::reallocateSpace(int newMaxRange){
    if(newMaxRange>maxRange){
        coefVecLocator   = (int **) realloc(coefVecLocator,(newMaxRange+1)*sizeof(int *));
        for(int r = (maxRange+1); r<=newMaxRange; r++){
            coefVecLocator[r] = (int *) malloc((2*r+1)*sizeof(int));
            for(int d = 1; d<=(2*r); d++ ){
                coefVecLocator[r][d] = 0;
            }
        }
    }
    else{
        printf("Error the new stencil must be larger than the old stencil to reallocate space in NumericalDerivative class\n");
        exit(-1);
    }
}

double NumericalDeriv::mixedNumDeriv(double *u, int narg, int index[3], int Deriv[3], int order[3], double h[3], int *size){
        
    if(narg<3){
        for(int d = narg; d<3; d++){
            index[d] = 0; Deriv[d] = 0; order[d] = 0; h[d] = 0;
        }
    }
    double D;
    
    int r[3], s[3];
    for(int d = 0; d<3; d++){
        if((order[d]%2)!=0)
            order[d] = order[d] + 1;
        r[d] = RANGE(Deriv[d], order[d]);
        s[d] = (2*r[d] + 1);
    }
    double *V[3];
    int so[] = setVecTo(size);
    int ro[] = setVecTo(index);
    for(int d = 0; d<3; d++){
        s[d] = 1; r[d] = 0;
        V[d] = (double *) malloc((s[0]*s[1]*s[2])*sizeof(double));
            int i[3];
            for(i[2] = (-r[2]); i[2]<=(r[2]); i[2]++){
                for(i[1] = (-r[1]); i[1]<=(r[1]); i[1]++){
                    for(i[0] = (-r[0]); i[0]<=(r[0]); i[0]++){
                        int Inew[] = sumVectors(r,i);
                        int Iold[] = sumVectors(ro,i);
                        if(d == 0){
                            V[d][ind(Inew,s)] = numDeriv(u,3,Iold,d, Deriv[d], order[d], h[d],so);
                        }else{
                            V[d][ind(Inew,s)] = numDeriv(V[(d-1)],3,Iold,d, Deriv[d], order[d], h[d],so);
                        }
                    }// end of i[0] loop
                }// end of i[1] loop
            }// end of i[2] loop
        if(d<2){
            memcpy(so, s, sizeof(s));
            memcpy(ro, r, sizeof(r));
        }
    }
    D = V[2][0];
    
    for(int d = 0; d<3; d++)
    {free(V[d]); V[d] = NULL;}
    
    return D;
}

double NumericalDeriv::numDeriv(double *u, int narg, int index[3], int pos, int d, int order, double h, int *size){
    
    /* This function computes the [d]th centered numerical derivative of a grid function [u] of [size]. The derivative is computed with respect to the [pos]th argument evaluated at [index] to order [order] with step size [h] */
    int ord = order;
    
    if(narg<3){
        for(int d = narg; d<3; d++)
            index[d] = 0;
    }
    
    double D; // the variable that will hold the result
    
    double *c;
    if(d!=0){
        if ((ord%2) != 0) {
            ord = ord + 1;
        }
        int r = RANGE(d,ord);
        if(r!=0){
            c = getDerivCoef(r,d);}
        else{
            D = u[ind(index,size)];
            return D;
        }

        int varIndex = index[pos];
        D = 0;
        int count = 0;
        for(int k=(-r); k<=r; k++){
            
            index[pos] =(varIndex + k);
            
            if(index[pos]<0){
                printf("DANGER:negative index in numDeriv!\n");
                exit(-1);
            }
            
            D = D + c[count]*u[ind(index,size)];
            count = count + 1;
        }
        index[pos] = varIndex;
        D = D*(1.0/pow(h,(double) d));
        
        c = NULL;
    }
    else{ D = u[ind(index,size)];}
    
    return D;
}

double NumericalDeriv::numDerivFn(double (*F) (double *), double *arg, int pos, int d, int order, double h){
    
    /* This function computes the [d]th centered numerical derivative of a function pointer [F]. The derivative is computed with respect to the [pos]th argument evaluated at [arg] to order [order] with step size [h] */
    
    int r = 0; // the range needed for the stencil
    double D; // the variable that will hold the result
    
    double *c;
    if(d!=0){
        if ((order%2) != 0) {
            order = order + 1;
        }
        r = RANGE(d,order);
        if(r!=0){
            c = getDerivCoef(r,d);}
        else{c = 0;}
    }else{c = 0;}
    
    if(d!=0){
        double varArg = arg[pos];
        D = 0;
        int count = 0;
        for(int k=(-r); k<=r; k++){
            
            arg[pos] =(varArg + k*h);
            
            D = D + c[count]*F(arg);
            count = count + 1;
        }
        arg[pos] = varArg;
        D = D*(1.0/pow(h,(double) d));
    }
    else{ D = F(arg);}
    c = NULL;
    return D;
}

double * NumericalDeriv::getDerivCoef(int r,int d){
    if(r>maxRange){
        printf("Error in getDerivCoef: the stencil of derivative must be within -%d:%d stencil. To expand the stencil use NumericalDerivative member function expandSpace(%d) to get a %d:%d stencil.\n",maxRange,maxRange,r,r,r);
        exit(-1);
    }
    double *c = nullptr;
    if(allocateSpaceFlag == false){
        allocateSpace();
        allocateSpaceFlag = true;
    }
    if(coefVecLocator[r][d] == 0){
        coefVecLocator[r][d] = coefVecCount;
        coefVec = (double**) realloc(coefVec,coefVecCount*sizeof(double*));
        coefVec[(coefVecCount-1)] = (double *) malloc((2*r+1)*sizeof(double));
        centeredDifferenceWeights(r,d,(coefVec[(coefVecCount-1)]));
        c = coefVec[(coefVecCount-1)];
        coefVecCount++;
    }
    else{
        c =coefVec[(coefVecLocator[r][d]-1)];
    }
    return c;
}

void NumericalDeriv::centeredDifferenceWeights(int r, int d, double *c)
{
    // r is the stencil for the numerical derivative
    // d is the number of derivatives
    // c is the vector that will store the coefficients
    
    double c1,c2,c3,c4,c5;
    int mn;
    
    int n = 2*r;
    int x[(2*r+1)];
    int coefVecCount = 0;
    for(int i = -r; i<=r; i++){
        x[coefVecCount] = i;
        coefVecCount = coefVecCount + 1;
    }
    
    double cl[(n+1)][(d+1)];
    memset(cl, 0, sizeof(cl));
    
    c1 = 1;
    c4 = x[0];
    
    cl[0][0] = 1;
    for (int i=1; i<=n; i++) {
        mn = MIN(i,d);
        c2 = 1;
        c5 = c4;
        c4 = x[i];
        
        for (int j = 0; j<=(i-1); j++) {
            c3 = x[i]-x[j];
            c2 = c2*c3;
            if (j == (i-1)) {
                for (int k = mn; k>=1; k=k-1) {
                    cl[i][k] = c1*(k*cl[(i-1)][(k-1)]-c5*cl[(i-1)][k])/c2;
                }
                cl[i][0] = -c1*c5*cl[(i-1)][0]/c2;
            }
            for (int k = mn; k>=1; k=k-1) {
                cl[j][k] = (c4*cl[j][k]-k*cl[j][(k-1)])/c3;
            }
            cl[j][0] = c4*cl[j][0]/c3;
        }
        c1 = c2;
    }
    
    for(int i = 0; i<(n+1); i++)
    {
        c[i] = cl[i][d];
    }
    return;
}


