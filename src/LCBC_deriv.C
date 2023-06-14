#include "LCBC.h"

#define pInd(k,deriv) (k + deriv*orderNum)

double Lcbc::mixedDeriv(double *u, int index[3], int deriv[3], double dx[3], int order[3], int size[3]){

    double D = 0;

    int r[3] = {derivRange[pInd((order[0]/2),deriv[0])],
                derivRange[pInd((order[1]/2),deriv[1])],
                derivRange[pInd((order[2]/2),deriv[2])]};
    
    int cnt[3] = {0,0,0};
    int i[3];
    for(i[2] = (-r[2]); i[2]<=r[2]; i[2]++){
        for(i[1] = (-r[1]); i[1]<=r[1]; i[1]++){
            for(i[0] = (-r[0]); i[0]<=r[0]; i[0]++ ){
                int Ind[3] = {(index[0] + i[0]),(index[1] + i[1]),(index[2] + i[2])};
                double c[3] = {derivCoef[pInd((order[0]/2),deriv[0])][cnt[0]],
                               derivCoef[pInd((order[1]/2),deriv[1])][cnt[1]],
                               derivCoef[pInd((order[2]/2),deriv[2])][cnt[2]]};
                
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

double Lcbc::Deriv(double (*F) (double *), double *arg, int pos, int deriv, double dx, int order){
    double D = 0;

    int r = derivRange[pInd((order/2),deriv)];

    double argVar = arg[(pos+1)];
    int cnt = 0;
    for(int i = (-r); i<=r; i++ ){
        arg[(pos + 1)] = argVar + dx*i;
        D = D + derivCoef[pInd((order/2),deriv)][cnt]*F(arg);
        cnt++;
    }// end of i2
    
    arg[(pos + 1)] = argVar;
    D = D*(1.0/pow(dx,(double) deriv));

    return D;
}// end of Deriv

double Lcbc::Deriv(double *u, int index[3], int pos, int deriv, double dx, int order, int size[3]){
    double D = 0;

    int r = derivRange[pInd((order/2),deriv)];

    int varIndex = index[pos];
    int cnt = 0;
    for(int i = (-r); i<=r; i++ ){
        index[pos] = varIndex + i;
        D = D + derivCoef[pInd((order/2),deriv)][cnt]*u[ind(index,size)];
        cnt++;
    }// end of i2
    index[pos] = varIndex;
    
    D = D*(1.0/pow(dx,(double) deriv));

    return D;
}// end of Deriv


void Lcbc::getDerivCoef(int p){
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
