#include "LagrangeDerivFunctions.h"

void LagrangeDeriv_2D_8_0(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    
    int p = 4;
    int wth[3] = {2*p,2*p,0}; wth[axis] = p;
    int lth[3] = {(2*wth[0]+1), (2*wth[1]+1), 1};
    
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 1; // other axis 1
    int v1 = 2; // other axis 2
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];

    double d20[Lth], d02[Lth], d22[Lth];
    double d40[Lth], d04[Lth];
    double Q11_20[Lth], Q11_02[Lth];
    double d42[Lth], d24[Lth], d60[Lth], d06[Lth];
    double Q11_22[Lth], Q11_40[Lth], Q11_04[Lth];
    double Q21_20[Lth], Q21_02[Lth], Q12_20[Lth], Q12_02[Lth];

    int i[3]; i[2] = 0;
    for(i[1] = 0; i[1]<lth[1]; i[1]++){
        for(i[0] = 0; i[0]<=lth[0]; i[0]++){
            Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)];
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
        for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
            int I = ind(i,lth);
            d20[I] = DERIV20(Q[ind2(0,1,NU,K)],i,lth,dx);
            d02[I] = DERIV02(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            double d10 = DERIV10(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d01 = DERIV01(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d11 = DERIV11(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d20[I])+coef[4][cI]*(d11)+coef[1][cI]*(d02[I])+coef[2][cI]*(d10)+coef[3][cI]*(d01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
        for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
            int I = ind(i,lth);
            d40[I] = DERIV20(d20, i, lth, dx);
            d04[I] = DERIV02(d02, i, lth, dx);
            d22[I] = DERIV20(d02, i, lth, dx);
            
            double d13 = DERIV11(d02, i, lth, dx);
            double d31 = DERIV11(d20, i, lth, dx);
            double d03 = DERIV01(d02, i, lth, dx);
            double d30 = DERIV10(d20, i, lth, dx);
            
            Q11_20[I] = DERIV20(Q[ind2(1,1,NU,K)],i, lth, dx);
            Q11_02[I] = DERIV02(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_11 = DERIV11(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_10 = DERIV10(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_01 = DERIV01(Q[ind2(1,1,NU,K)],i, lth, dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d40[I])+coef[4][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d31)+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d04[I])+coef[2][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d30)+coef[3][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d03);
            Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_20[I])+coef[4][cI]*(Q11_11)+coef[1][cI]*(Q11_02[I])+coef[2][cI]*(Q11_10)+coef[3][cI]*(Q11_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 3; i[1]<(lth[1]-3); i[1]++){
        for(i[0] = 3; i[0]<(lth[0]-3); i[0]++){
            int I = ind(i,lth);
            
            d42[I] = DERIV20(d22,i,lth,dx);
            d24[I] = DERIV02(d22,i,lth,dx);
            d60[I] = DERIV20(d40,i,lth,dx);
            d06[I] = DERIV02(d04,i,lth,dx);

            double d33 = DERIV11(d22,i,lth,dx);
            double d50 = DERIV10(d40,i,lth,dx);
            double d05 = DERIV01(d04,i,lth,dx);
            double d51 = DERIV11(d40,i,lth,dx);
            double d15 = DERIV11(d04,i,lth,dx);
            
            Q11_22[I] = DERIV20(Q11_02,i,lth,dx);
            Q11_40[I] = DERIV20(Q11_20,i, lth, dx);
            Q11_04[I] = DERIV02(Q11_02,i, lth, dx);
            Q21_20[I] = DERIV20(Q[ind2(2,1,NU,K)],i, lth, dx);
            Q21_02[I] = DERIV02(Q[ind2(2,1,NU,K)],i, lth, dx);
            Q12_20[I] = DERIV20(Q[ind2(1,2,NU,K)],i, lth, dx);
            Q12_02[I] = DERIV02(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q11_30 = DERIV10(Q11_20,i, lth, dx);
            double Q11_03 = DERIV01(Q11_02,i, lth, dx);
            double Q11_31 = DERIV11(Q11_20,i, lth, dx);
            double Q11_13 = DERIV11(Q11_02,i, lth, dx);
            
            double Q21_10 = DERIV10(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q21_01 = DERIV01(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q21_11 = DERIV11(Q[ind2(2,1,NU,K)],i, lth, dx);
            
            double Q12_10 = DERIV10(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q12_01 = DERIV01(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q12_11 = DERIV11(Q[ind2(1,2,NU,K)],i, lth, dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,3,NU,K)][I] = Q[ind2(1,2,NU,K)][I]+coef[0][cI]*((1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d60[I])+coef[4][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d15+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*d33+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d51)+coef[1][cI]*((1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d06[I])+coef[2][cI]*((1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d50)+coef[3][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d05);
            Q[ind2(2,2,NU,K)][I] = coef[0][cI]*(Q12_20[I]+((-1.0/12.0)*(dx[0]*dx[0])*Q11_40[I]))+coef[4][cI]*(Q12_11+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q11_13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_31))+coef[1][cI]*(Q12_02[I]+((-1.0/12.0)*(dx[1]*dx[1])*Q11_04[I]))+coef[2][cI]*(Q12_10+((-1.0/6.0)*(dx[0]*dx[0])*Q11_30))+coef[3][cI]*(Q12_01+((-1.0/6.0)*(dx[1]*dx[1])*Q11_03));
            Q[ind2(3,1,NU,K)][I] = coef[0][cI]*(Q21_20[I])+coef[4][cI]*(Q21_11)+coef[1][cI]*(Q21_02[I])+coef[2][cI]*(Q21_10)+coef[3][cI]*(Q21_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 4; i[1]<(lth[1]-4); i[1]++){
        for(i[0] = 4; i[0]<(lth[0]-4); i[0]++){
            int I = ind(i,lth);
            
            double d80 = DERIV20(d60,i,lth,dx);
            double d08 = DERIV02(d06,i,lth,dx);

            double d53 = DERIV11(d42,i,lth,dx);
            double d35 = DERIV11(d24,i,lth,dx);
            double d70 = DERIV10(d60,i,lth,dx);
            double d07 = DERIV01(d06,i,lth,dx);
            double d71 = DERIV11(d60,i,lth,dx);
            double d17 = DERIV11(d06,i,lth,dx);

            double Q11_60 = DERIV20(Q11_40,i,lth,dx);
            double Q11_06 = DERIV02(Q11_04,i,lth,dx);
            double Q21_40 = DERIV20(Q21_20,i,lth,dx);
            double Q21_04 = DERIV02(Q21_02,i,lth,dx);
            double Q12_40 = DERIV20(Q12_20,i,lth,dx);
            double Q12_04 = DERIV02(Q12_02,i,lth,dx);
            double Q22_20 = DERIV20(Q[ind2(2,2,NU,K)],i,lth,dx);
            double Q22_02 = DERIV02(Q[ind2(2,2,NU,K)],i,lth,dx);
            double Q31_20 = DERIV20(Q[ind2(3,1,NU,K)],i,lth,dx);
            double Q31_02 = DERIV02(Q[ind2(3,1,NU,K)],i,lth,dx);
            double Q13_20 = DERIV20(Q[ind2(1,3,NU,K)],i,lth,dx);
            double Q13_02 = DERIV02(Q[ind2(1,3,NU,K)],i,lth,dx);
            double Q11_33 = DERIV11(Q11_22,i,lth,dx);
            double Q11_50 = DERIV10(Q11_40,i,lth,dx);
            double Q11_05 = DERIV01(Q11_04,i,lth,dx);
            double Q11_51 = DERIV11(Q11_40,i,lth,dx);
            double Q11_15 = DERIV11(Q11_04,i,lth,dx);

            double Q21_30 = DERIV10(Q21_20,i,lth,dx);
            double Q21_03 = DERIV01(Q21_02,i,lth,dx);
            double Q21_31 = DERIV11(Q21_20,i,lth,dx);
            double Q21_13 = DERIV11(Q21_02,i,lth,dx);

            double Q12_30 = DERIV10(Q12_20,i,lth,dx);
            double Q12_03 = DERIV01(Q12_02,i,lth,dx);
            double Q12_31 = DERIV11(Q12_20,i,lth,dx);
            double Q12_13 = DERIV11(Q12_02,i,lth,dx);

            double Q22_10 = DERIV10(Q[ind2(2,2,NU,K)],i,lth,dx);
            double Q22_01 = DERIV01(Q[ind2(2,2,NU,K)],i,lth,dx);
            double Q22_11 = DERIV11(Q[ind2(2,2,NU,K)],i,lth,dx);

            double Q31_10 = DERIV10(Q[ind2(3,1,NU,K)],i,lth,dx);
            double Q31_01 = DERIV01(Q[ind2(3,1,NU,K)],i,lth,dx);
            double Q31_11 = DERIV11(Q[ind2(3,1,NU,K)],i,lth,dx);

            double Q13_10 = DERIV10(Q[ind2(1,3,NU,K)],i,lth,dx);
            double Q13_01 = DERIV01(Q[ind2(1,3,NU,K)],i,lth,dx);
            double Q13_11 = DERIV11(Q[ind2(1,3,NU,K)],i,lth,dx);

            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,4,NU,K)][I] = Q[ind2(1,3,NU,K)][I]+coef[0][cI]*((-1.0/560.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d80)+coef[4][cI]*((-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d17+(-1.0/6.0)*(1.0/30.0)*(dx[0]*dx[0])*(dx[1]*dx[1]*dx[1]*dx[1])*d35+(1.0/30.0)*(-1.0/6.0)*(dx[0]*dx[0]*dx[0]*dx[0])*(dx[1]*dx[1])*d53+(-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d71)+coef[1][cI]*((-1.0/560.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d08)+coef[2][cI]*((-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d70)+coef[3][cI]*((-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d07);
            Q[ind2(2,3,NU,K)][I] = coef[0][cI]*(Q13_20+((-1.0/12.0)*(dx[0]*dx[0])*Q12_40+(1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_60))+coef[4][cI]*(Q13_11+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q12_13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q12_31+1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_15+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*Q11_33+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*Q11_51))+coef[1][cI]*(Q13_02+((-1.0/12.0)*(dx[1]*dx[1])*Q12_04+(1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_06))+coef[2][cI]*(Q13_10+((-1.0/6.0)*(dx[0]*dx[0])*Q12_30+(1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_50))+coef[3][cI]*(Q13_01+((-1.0/6.0)*(dx[1]*dx[1])*Q12_03+(1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_05));
            Q[ind2(3,2,NU,K)][I] = coef[0][cI]*(Q22_20+((-1.0/12.0)*(dx[0]*dx[0])*Q21_40))+coef[4][cI]*(Q22_11+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q21_13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q21_31))+coef[1][cI]*(Q22_02+((-1.0/12.0)*(dx[1]*dx[1])*Q21_04))+coef[2][cI]*(Q22_10+((-1.0/6.0)*(dx[0]*dx[0])*Q21_30))+coef[3][cI]*(Q22_01+((-1.0/6.0)*(dx[1]*dx[1])*Q21_03));
            Q[ind2(4,1,NU,K)][I] = coef[0][cI]*(Q31_20)+coef[4][cI]*(Q31_11)+coef[1][cI]*(Q31_02)+coef[2][cI]*(Q31_10)+coef[3][cI]*(Q31_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind2(0,0,NU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind2(1,0,NU,MU)] = Q[ind2(1,4,NU,K)][I];
    Z[ind2(2,0,NU,MU)] = Q[ind2(2,3,NU,K)][I];
    Z[ind2(3,0,NU,MU)] = Q[ind2(3,2,NU,K)][I];
    Z[ind2(4,0,NU,MU)] = Q[ind2(4,1,NU,K)][I];
    
    double d2[Lth];
    double d4[Lth];
    double d6[Lth];
    int j[3] = {wth[0],wth[1],wth[2]};

    int nu = 0;
    double d01 = DERIV01(Q[ind2(nu,1,NU,K)], wth, lth, dx);
    double d03 = DERIV01(d02, wth, lth, dx);
    double d05 = DERIV01(d04, wth, lth, dx);
    double d07 = DERIV01(d06, wth, lth, dx);
    double d08 = DERIV02(d06, wth, lth, dx);
    Z[ind2(nu,1,NU,MU)] = d01 - ((dx[1]*dx[1])/6.0)*d03 + ((dx[1]*dx[1]*dx[1]*dx[1])/30.0)*d05 - ((dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])/140.0)*d07 ;
    Z[ind2(nu,2,NU,MU)] = d02[I] - ((dx[1]*dx[1])/12.0)*d04[I] + ((dx[1]*dx[1]*dx[1]*dx[1])/90.0)*d06[I] - ((dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])/560.0)*d08;
    Z[ind2(nu,3,NU,MU)] = d03 - ((dx[1]*dx[1])/4.0)*d05 + ((7.0*dx[1]*dx[1]*dx[1]*dx[1])/120.0)*d07;
    Z[ind2(nu,4,NU,MU)] = d04[I] - ((dx[1]*dx[1])/6.0)*d06[I] + ((7.0*dx[1]*dx[1]*dx[1]*dx[1])/240.0)*d08;
    Z[ind2(nu,5,NU,MU)] = d05 - ((dx[1]*dx[1])/3.0)*d07; // to order 4
    Z[ind2(nu,6,NU,MU)] = d06[I] - ((dx[1]*dx[1])/4.0)*d08; // to order 4
    Z[ind2(nu,7,NU,MU)] = d07; // to order 2
    Z[ind2(nu,8,NU,MU)] = d08; // to order 2

    int kappa1[] = {1, 1, 1, 0};
    int kappa2[] = {1, 1, 0, 0};
    int kappa3[] = {1, 0, 0, 0};

    for(int nu = 1; nu<=p; nu++){
            int k = (p+1-nu);
            j[0] = wth[0]; j[2] = wth[2];
            for(j[1] = (wth[1]-3); j[1]<=(wth[1]+3); j[1]++){
                int I = ind(j,lth);
                d2[I] = DERIV02(Q[ind2(nu,k,NU,K)], j, lth, dx);
            }
            for(j[1] = (wth[1]-2); j[1]<=(wth[1]+2); j[1]++){
                int I = ind(j,lth);
                d4[I] = DERIV02(d2, j, lth, dx);
            }
            for(j[1] = (wth[1]-1); j[1]<=(wth[1]+1); j[1]++){
                int I = ind(j,lth);
                d6[I] = DERIV02(d4, j, lth, dx);
            }
            int I = ind(wth,lth);
            double d1 = DERIV01(Q[ind2(nu,k,NU,K)], wth, lth, dx);
            double d3 = DERIV01(d2, wth, lth, dx);
            double d5 = DERIV01(d4, wth, lth, dx);
            double d7 = DERIV01(d6, wth, lth, dx);
            double d8 = DERIV02(d6, wth, lth, dx);

            Z[ind2(nu,1,NU,MU)] = d1 - ((dx[1]*dx[1])/6.0)*d3*kappa1[(nu-1)] +((dx[1]*dx[1]*dx[1]*dx[1])/30.0)*d5*kappa2[(nu-1)]- ((dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])/140.0)*d7*kappa3[(nu-1)];
            Z[ind2(nu,2,NU,MU)] = d2[I] - ((dx[1]*dx[1])/12.0)*d4[I]*kappa1[(nu-1)] +((dx[1]*dx[1]*dx[1]*dx[1])/90.0)*d6[I]*kappa2[(nu-1)] - ((dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])/560.0)*d8*kappa3[(nu-1)];
            Z[ind2(nu,3,NU,MU)] = d3 - ((dx[1]*dx[1])/4.0)*d5*kappa2[(nu-1)] + ((7.0*dx[1]*dx[1]*dx[1]*dx[1])/120.0)*d7*kappa3[(nu-1)];
            Z[ind2(nu,4,NU,MU)] = d4[I] - ((dx[1]*dx[1])/6.0)*d6[I]*kappa2[(nu-1)]+ ((7.0*dx[1]*dx[1]*dx[1]*dx[1])/240.0)*d8*kappa3[(nu-1)];
            Z[ind2(nu,5,NU,MU)] = d5 - ((dx[1]*dx[1])/3.0)*d7*kappa3[(nu-1)];
            Z[ind2(nu,6,NU,MU)] = d6[I] - ((dx[1]*dx[1])/4.0)*d8*kappa3[(nu-1)];
            Z[ind2(nu,7,NU,MU)] = d7;
            Z[ind2(nu,8,NU,MU)] = d8;
    }
}

void LagrangeDeriv_2D_8_1(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    
    int p = 4;
    int wth[3] = {2*p,2*p,0}; wth[axis] = p;
    int lth[3] = {(2*wth[0]+1), (2*wth[1]+1), 1};
    
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 0; // other axis 1
    int v1 = 2; // other axis 2
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d20[Lth], d02[Lth], d22[Lth];
    double d40[Lth], d04[Lth];
    double Q11_20[Lth], Q11_02[Lth];
    double d42[Lth], d24[Lth], d60[Lth], d06[Lth];
    double Q11_22[Lth], Q11_40[Lth], Q11_04[Lth];
    double Q21_20[Lth], Q21_02[Lth], Q12_20[Lth], Q12_02[Lth];

    int i[3]; i[2] = 0;
    for(i[1] = 0; i[1]<lth[1]; i[1]++){
        for(i[0] = 0; i[0]<=lth[0]; i[0]++){
            Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)];
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
        for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
            int I = ind(i,lth);
            d20[I] = DERIV20(Q[ind2(0,1,NU,K)],i,lth,dx);
            d02[I] = DERIV02(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            double d10 = DERIV10(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d01 = DERIV01(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d11 = DERIV11(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d20[I])+coef[4][cI]*(d11)+coef[1][cI]*(d02[I])+coef[2][cI]*(d10)+coef[3][cI]*(d01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
        for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
            int I = ind(i,lth);
            d40[I] = DERIV20(d20, i, lth, dx);
            d04[I] = DERIV02(d02, i, lth, dx);
            d22[I] = DERIV20(d02, i, lth, dx);
            
            double d13 = DERIV11(d02, i, lth, dx);
            double d31 = DERIV11(d20, i, lth, dx);
            double d03 = DERIV01(d02, i, lth, dx);
            double d30 = DERIV10(d20, i, lth, dx);
            
            Q11_20[I] = DERIV20(Q[ind2(1,1,NU,K)],i, lth, dx);
            Q11_02[I] = DERIV02(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_11 = DERIV11(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_10 = DERIV10(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_01 = DERIV01(Q[ind2(1,1,NU,K)],i, lth, dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d40[I])+coef[4][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d31)+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d04[I])+coef[2][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d30)+coef[3][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d03);
            Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_20[I])+coef[4][cI]*(Q11_11)+coef[1][cI]*(Q11_02[I])+coef[2][cI]*(Q11_10)+coef[3][cI]*(Q11_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 3; i[1]<(lth[1]-3); i[1]++){
        for(i[0] = 3; i[0]<(lth[0]-3); i[0]++){
            int I = ind(i,lth);
            
            d42[I] = DERIV20(d22,i,lth,dx);
            d24[I] = DERIV02(d22,i,lth,dx);
            d60[I] = DERIV20(d40,i,lth,dx);
            d06[I] = DERIV02(d04,i,lth,dx);

            double d33 = DERIV11(d22,i,lth,dx);
            double d50 = DERIV10(d40,i,lth,dx);
            double d05 = DERIV01(d04,i,lth,dx);
            double d51 = DERIV11(d40,i,lth,dx);
            double d15 = DERIV11(d04,i,lth,dx);
            
            Q11_22[I] = DERIV20(Q11_02,i,lth,dx);
            Q11_40[I] = DERIV20(Q11_20,i, lth, dx);
            Q11_04[I] = DERIV02(Q11_02,i, lth, dx);
            Q21_20[I] = DERIV20(Q[ind2(2,1,NU,K)],i, lth, dx);
            Q21_02[I] = DERIV02(Q[ind2(2,1,NU,K)],i, lth, dx);
            Q12_20[I] = DERIV20(Q[ind2(1,2,NU,K)],i, lth, dx);
            Q12_02[I] = DERIV02(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q11_30 = DERIV10(Q11_20,i, lth, dx);
            double Q11_03 = DERIV01(Q11_02,i, lth, dx);
            double Q11_31 = DERIV11(Q11_20,i, lth, dx);
            double Q11_13 = DERIV11(Q11_02,i, lth, dx);
            
            double Q21_10 = DERIV10(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q21_01 = DERIV01(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q21_11 = DERIV11(Q[ind2(2,1,NU,K)],i, lth, dx);
            
            double Q12_10 = DERIV10(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q12_01 = DERIV01(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q12_11 = DERIV11(Q[ind2(1,2,NU,K)],i, lth, dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,3,NU,K)][I] = Q[ind2(1,2,NU,K)][I]+coef[0][cI]*((1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d60[I])+coef[4][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d15+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*d33+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d51)+coef[1][cI]*((1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d06[I])+coef[2][cI]*((1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d50)+coef[3][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d05);
            Q[ind2(2,2,NU,K)][I] = coef[0][cI]*(Q12_20[I]+((-1.0/12.0)*(dx[0]*dx[0])*Q11_40[I]))+coef[4][cI]*(Q12_11+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q11_13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_31))+coef[1][cI]*(Q12_02[I]+((-1.0/12.0)*(dx[1]*dx[1])*Q11_04[I]))+coef[2][cI]*(Q12_10+((-1.0/6.0)*(dx[0]*dx[0])*Q11_30))+coef[3][cI]*(Q12_01+((-1.0/6.0)*(dx[1]*dx[1])*Q11_03));
            Q[ind2(3,1,NU,K)][I] = coef[0][cI]*(Q21_20[I])+coef[4][cI]*(Q21_11)+coef[1][cI]*(Q21_02[I])+coef[2][cI]*(Q21_10)+coef[3][cI]*(Q21_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 4; i[1]<(lth[1]-4); i[1]++){
        for(i[0] = 4; i[0]<(lth[0]-4); i[0]++){
            int I = ind(i,lth);
            
            double d80 = DERIV20(d60,i,lth,dx);
            double d08 = DERIV02(d06,i,lth,dx);

            double d53 = DERIV11(d42,i,lth,dx);
            double d35 = DERIV11(d24,i,lth,dx);
            double d70 = DERIV10(d60,i,lth,dx);
            double d07 = DERIV01(d06,i,lth,dx);
            double d71 = DERIV11(d60,i,lth,dx);
            double d17 = DERIV11(d06,i,lth,dx);

            double Q11_60 = DERIV20(Q11_40,i,lth,dx);
            double Q11_06 = DERIV02(Q11_04,i,lth,dx);
            double Q21_40 = DERIV20(Q21_20,i,lth,dx);
            double Q21_04 = DERIV02(Q21_02,i,lth,dx);
            double Q12_40 = DERIV20(Q12_20,i,lth,dx);
            double Q12_04 = DERIV02(Q12_02,i,lth,dx);
            double Q22_20 = DERIV20(Q[ind2(2,2,NU,K)],i,lth,dx);
            double Q22_02 = DERIV02(Q[ind2(2,2,NU,K)],i,lth,dx);
            double Q31_20 = DERIV20(Q[ind2(3,1,NU,K)],i,lth,dx);
            double Q31_02 = DERIV02(Q[ind2(3,1,NU,K)],i,lth,dx);
            double Q13_20 = DERIV20(Q[ind2(1,3,NU,K)],i,lth,dx);
            double Q13_02 = DERIV02(Q[ind2(1,3,NU,K)],i,lth,dx);
            double Q11_33 = DERIV11(Q11_22,i,lth,dx);
            double Q11_50 = DERIV10(Q11_40,i,lth,dx);
            double Q11_05 = DERIV01(Q11_04,i,lth,dx);
            double Q11_51 = DERIV11(Q11_40,i,lth,dx);
            double Q11_15 = DERIV11(Q11_04,i,lth,dx);

            double Q21_30 = DERIV10(Q21_20,i,lth,dx);
            double Q21_03 = DERIV01(Q21_02,i,lth,dx);
            double Q21_31 = DERIV11(Q21_20,i,lth,dx);
            double Q21_13 = DERIV11(Q21_02,i,lth,dx);

            double Q12_30 = DERIV10(Q12_20,i,lth,dx);
            double Q12_03 = DERIV01(Q12_02,i,lth,dx);
            double Q12_31 = DERIV11(Q12_20,i,lth,dx);
            double Q12_13 = DERIV11(Q12_02,i,lth,dx);

            double Q22_10 = DERIV10(Q[ind2(2,2,NU,K)],i,lth,dx);
            double Q22_01 = DERIV01(Q[ind2(2,2,NU,K)],i,lth,dx);
            double Q22_11 = DERIV11(Q[ind2(2,2,NU,K)],i,lth,dx);

            double Q31_10 = DERIV10(Q[ind2(3,1,NU,K)],i,lth,dx);
            double Q31_01 = DERIV01(Q[ind2(3,1,NU,K)],i,lth,dx);
            double Q31_11 = DERIV11(Q[ind2(3,1,NU,K)],i,lth,dx);

            double Q13_10 = DERIV10(Q[ind2(1,3,NU,K)],i,lth,dx);
            double Q13_01 = DERIV01(Q[ind2(1,3,NU,K)],i,lth,dx);
            double Q13_11 = DERIV11(Q[ind2(1,3,NU,K)],i,lth,dx);

            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,4,NU,K)][I] = Q[ind2(1,3,NU,K)][I]+coef[0][cI]*((-1.0/560.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d80)+coef[4][cI]*(1.0*(-1.0/140.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d17+(-1.0/6.0)*(1.0/30.0)*(dx[0]*dx[0])*(dx[1]*dx[1]*dx[1]*dx[1])*d35+(1.0/30.0)*(-1.0/6.0)*(dx[0]*dx[0]*dx[0]*dx[0])*(dx[1]*dx[1])*d53+(-1.0/140.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d71)+coef[1][cI]*((-1.0/560.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d08)+coef[2][cI]*((-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d70)+coef[3][cI]*((-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d07);
            Q[ind2(2,3,NU,K)][I] = coef[0][cI]*(Q13_20+((-1.0/12.0)*(dx[0]*dx[0])*Q12_40+(1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_60))+coef[4][cI]*(Q13_11+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q12_13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q12_31+1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_15+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*Q11_33+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*Q11_51))+coef[1][cI]*(Q13_02+((-1.0/12.0)*(dx[1]*dx[1])*Q12_04+(1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_06))+coef[2][cI]*(Q13_10+((-1.0/6.0)*(dx[0]*dx[0])*Q12_30+(1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_50))+coef[3][cI]*(Q13_01+((-1.0/6.0)*(dx[1]*dx[1])*Q12_03+(1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_05));
            Q[ind2(3,2,NU,K)][I] = coef[0][cI]*(Q22_20+((-1.0/12.0)*(dx[0]*dx[0])*Q21_40))+coef[4][cI]*(Q22_11+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q21_13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q21_31))+coef[1][cI]*(Q22_02+((-1.0/12.0)*(dx[1]*dx[1])*Q21_04))+coef[2][cI]*(Q22_10+((-1.0/6.0)*(dx[0]*dx[0])*Q21_30))+coef[3][cI]*(Q22_01+((-1.0/6.0)*(dx[1]*dx[1])*Q21_03));
            Q[ind2(4,1,NU,K)][I] = coef[0][cI]*(Q31_20)+coef[4][cI]*(Q31_11)+coef[1][cI]*(Q31_02)+coef[2][cI]*(Q31_10)+coef[3][cI]*(Q31_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind2(0,0,NU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind2(1,0,NU,MU)] = Q[ind2(1,4,NU,K)][I];
    Z[ind2(2,0,NU,MU)] = Q[ind2(2,3,NU,K)][I];
    Z[ind2(3,0,NU,MU)] = Q[ind2(3,2,NU,K)][I];
    Z[ind2(4,0,NU,MU)] = Q[ind2(4,1,NU,K)][I];
    
    double d2[Lth];
    double d4[Lth];
    double d6[Lth];
    int j[3] = {wth[0],wth[1],wth[2]};

    int nu = 0;
    double d10 = DERIV10(Q[ind2(nu,1,NU,K)], wth, lth, dx);
    double d30 = DERIV10(d20, wth, lth, dx);
    double d50 = DERIV10(d40, wth, lth, dx);
    double d70 = DERIV10(d60, wth, lth, dx);
    double d80 = DERIV20(d60, wth, lth, dx);
    Z[ind2(nu,1,NU,MU)] = d10 - ((dx[0]*dx[0])/6.0)*d30 + ((dx[0]*dx[0]*dx[0]*dx[0])/30.0)*d50 - ((dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])/140.0)*d70 ;
    Z[ind2(nu,2,NU,MU)] = d20[I] - ((dx[0]*dx[0])/12.0)*d40[I] + ((dx[0]*dx[0]*dx[0]*dx[0])/90.0)*d60[I] - ((dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])/560.0)*d80;
    Z[ind2(nu,3,NU,MU)] = d30 - ((dx[0]*dx[0])/4.0)*d50 + ((7.0*dx[0]*dx[0]*dx[0]*dx[0])/120.0)*d70;
    Z[ind2(nu,4,NU,MU)] = d40[I] - ((dx[0]*dx[0])/6.0)*d60[I] + ((7.0*dx[0]*dx[0]*dx[0]*dx[0])/240.0)*d80;
    Z[ind2(nu,5,NU,MU)] = d50 - ((dx[0]*dx[0])/3.0)*d70; // to order 4
    Z[ind2(nu,6,NU,MU)] = d60[I] - ((dx[0]*dx[0])/4.0)*d80; // to order 4
    Z[ind2(nu,7,NU,MU)] = d70; // to order 2
    Z[ind2(nu,8,NU,MU)] = d80; // to order 2

    int kappa1[] = {1, 1, 1, 0};
    int kappa2[] = {1, 1, 0, 0};
    int kappa3[] = {1, 0, 0, 0};

    for(int nu = 1; nu<=p; nu++){
        int k = (p+1-nu);
        j[1] = wth[1]; j[2] = wth[2];
        for(j[0] = (wth[0]-3); j[0]<=(wth[0]+3); j[0]++){
            int I = ind(j,lth);
            d2[I] = DERIV20(Q[ind2(nu,k,NU,K)], j, lth, dx);
        }
        for(j[0] = (wth[0]-2); j[0]<=(wth[0]+2); j[0]++){
            int I = ind(j,lth);
            d4[I] = DERIV20(d2, j, lth, dx);
        }
        for(j[0] = (wth[0]-1); j[0]<=(wth[0]+1); j[0]++){
            int I = ind(j,lth);
            d6[I] = DERIV20(d4, j, lth, dx);
        }
        int I = ind(wth,lth);
        double d1 = DERIV10(Q[ind2(nu,k,NU,K)], wth, lth, dx);
        double d3 = DERIV10(d2, wth, lth, dx);
        double d5 = DERIV10(d4, wth, lth, dx);
        double d7 = DERIV10(d6, wth, lth, dx);
        double d8 = DERIV20(d6, wth, lth, dx);

        Z[ind2(nu,1,NU,MU)] = d1 - ((dx[0]*dx[0])/6.0)*d3*kappa1[(nu-1)] +((dx[0]*dx[0]*dx[0]*dx[0])/30.0)*d5*kappa2[(nu-1)]- ((dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])/140.0)*d7*kappa3[(nu-1)];
        Z[ind2(nu,2,NU,MU)] = d2[I] - ((dx[0]*dx[0])/12.0)*d4[I]*kappa1[(nu-1)] +((dx[0]*dx[0]*dx[0]*dx[0])/90.0)*d6[I]*kappa2[(nu-1)] - ((dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])/560.0)*d8*kappa3[(nu-1)];
        Z[ind2(nu,3,NU,MU)] = d3 - ((dx[0]*dx[0])/4.0)*d5*kappa2[(nu-1)] + ((7.0*dx[0]*dx[0]*dx[0]*dx[0])/120.0)*d7*kappa3[(nu-1)];
        Z[ind2(nu,4,NU,MU)] = d4[I] - ((dx[0]*dx[0])/6.0)*d6[I]*kappa2[(nu-1)]+ ((7.0*dx[0]*dx[0]*dx[0]*dx[0])/240.0)*d8*kappa3[(nu-1)];
        Z[ind2(nu,5,NU,MU)] = d5 - ((dx[0]*dx[0])/3.0)*d7*kappa3[(nu-1)];
        Z[ind2(nu,6,NU,MU)] = d6[I] - ((dx[0]*dx[0])/4.0)*d8*kappa3[(nu-1)];
        Z[ind2(nu,7,NU,MU)] = d7;
        Z[ind2(nu,8,NU,MU)] = d8;
    }
}


void LagrangeDeriv_2D_6_0(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    
    int p = 3;
    int wth[3] = {2*p,2*p,0}; wth[axis] = p;
    int lth[3] = {(2*wth[0]+1), (2*wth[1]+1), 1};
    
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 1; // other axis 1
    int v1 = 2; // other axis 2
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d20[Lth], d02[Lth], d22[Lth];
    double Q11_20[Lth], Q11_02[Lth];
    double d40[Lth], d04[Lth];
    
    int i[3]; i[2] = 0;
    for(i[1] = 0; i[1]<lth[1]; i[1]++){
        for(i[0] = 0; i[0]<=lth[0]; i[0]++){
            Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)];
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
        for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
            int I = ind(i,lth);
            d20[I] = DERIV20(Q[ind2(0,1,NU,K)],i,lth,dx);
            d02[I] = DERIV02(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            double d10 = DERIV10(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d01 = DERIV01(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d11 = DERIV11(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d20[I])+coef[4][cI]*(d11)+coef[1][cI]*(d02[I])+coef[2][cI]*(d10)+coef[3][cI]*(d01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
        for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
            int I = ind(i,lth);
            d40[I] = DERIV20(d20, i, lth, dx);
            d04[I] = DERIV02(d02, i, lth, dx);
            d22[I] = DERIV20(d02, i, lth, dx);
            
            double d13 = DERIV11(d02, i, lth, dx);
            double d31 = DERIV11(d20, i, lth, dx);
            double d03 = DERIV01(d02, i, lth, dx);
            double d30 = DERIV10(d20, i, lth, dx);
            
            Q11_20[I] = DERIV20(Q[ind2(1,1,NU,K)],i, lth, dx);
            Q11_02[I] = DERIV02(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_11 = DERIV11(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_10 = DERIV10(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_01 = DERIV01(Q[ind2(1,1,NU,K)],i, lth, dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d40[I])+coef[4][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d31)+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d04[I])+coef[2][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d30)+coef[3][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d03);
            Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_20[I])+coef[4][cI]*(Q11_11)+coef[1][cI]*(Q11_02[I])+coef[2][cI]*(Q11_10)+coef[3][cI]*(Q11_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 3; i[1]<(lth[1]-3); i[1]++){
        for(i[0] = 3; i[0]<(lth[0]-3); i[0]++){
            int I = ind(i,lth);
            
            double d60 = DERIV20(d40, i, lth, dx);
            double d06 = DERIV02(d04, i, lth, dx);
            
            double d33 = DERIV11(d22, i, lth, dx);
            double d50 = DERIV10(d40, i, lth, dx);
            double d05 = DERIV01(d04, i, lth, dx);
            double d51 = DERIV11(d40, i, lth, dx);
            double d15 = DERIV11(d04, i, lth, dx);
            
            double Q11_40 = DERIV20(Q11_20,i, lth, dx);
            double Q11_04 = DERIV02(Q11_02,i, lth, dx);
            double Q21_20 = DERIV20(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q21_02 = DERIV02(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q12_20 = DERIV20(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q12_02 = DERIV02(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q11_30 = DERIV10(Q11_20,i, lth, dx);
            double Q11_03 = DERIV01(Q11_02,i, lth, dx);
            double Q11_31 = DERIV11(Q11_20,i, lth, dx);
            double Q11_13 = DERIV11(Q11_02,i, lth, dx);
            
            double Q21_10 = DERIV10(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q21_01 = DERIV01(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q21_11 = DERIV11(Q[ind2(2,1,NU,K)],i, lth, dx);
            
            double Q12_10 = DERIV10(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q12_01 = DERIV01(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q12_11 = DERIV11(Q[ind2(1,2,NU,K)],i, lth, dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,3,NU,K)][I] = Q[ind2(1,2,NU,K)][I]+coef[0][cI]*((1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d60)+coef[4][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d15+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*d33+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d51)+coef[1][cI]*((1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d06)+coef[2][cI]*((1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d50)+coef[3][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d05);
            Q[ind2(2,2,NU,K)][I] = coef[0][cI]*(Q12_20+((-1.0/12.0)*(dx[0]*dx[0])*Q11_40))+coef[4][cI]*(Q12_11+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q11_13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_31))+coef[1][cI]*(Q12_02+((-1.0/12.0)*(dx[1]*dx[1])*Q11_04))+coef[2][cI]*(Q12_10+((-1.0/6.0)*(dx[0]*dx[0])*Q11_30))+coef[3][cI]*(Q12_01+((-1.0/6.0)*(dx[1]*dx[1])*Q11_03));
            Q[ind2(3,1,NU,K)][I] = coef[0][cI]*(Q21_20)+coef[4][cI]*(Q21_11)+coef[1][cI]*(Q21_02)+coef[2][cI]*(Q21_10)+coef[3][cI]*(Q21_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind2(0,0,NU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind2(1,0,NU,MU)] = Q[ind2(1,3,NU,K)][I];
    Z[ind2(2,0,NU,MU)] = Q[ind2(2,2,NU,K)][I];
    Z[ind2(3,0,NU,MU)] = Q[ind2(3,1,NU,K)][I];
    
    double d2[Lth];
    double d4[Lth];
    int j[3] = {wth[0],wth[1],wth[2]};
    
    int nu = 0;
    double d01 = DERIV01(Q[ind2(nu,1,NU,K)], wth, lth, dx);
    double d03 = DERIV01(d02, wth, lth, dx);
    double d05 = DERIV01(d04, wth, lth, dx);
    double d06 = DERIV02(d04, wth, lth, dx);
    Z[ind2(nu,1,NU,MU)] = d01 - ((dx[1]*dx[1])/6.0)*d03 + ((dx[1]*dx[1]*dx[1]*dx[1])/30.0)*d05;
    Z[ind2(nu,2,NU,MU)] = d02[I] - ((dx[1]*dx[1])/12.0)*d04[I] + ((dx[1]*dx[1]*dx[1]*dx[1])/90.0)*d06;
    Z[ind2(nu,3,NU,MU)] = d03 - ((dx[1]*dx[1])/4.0)*d05;
    Z[ind2(nu,4,NU,MU)] = d04[I] - ((dx[1]*dx[1])/6.0)*d06;
    Z[ind2(nu,5,NU,MU)] = d05;
    Z[ind2(nu,6,NU,MU)] = d06;
    
    int kappa1[] = {1, 1, 0};
    int kappa2[] = {1, 0, 0};
    
    for(int nu = 1; nu<=p; nu++){
            int k = (p+1-nu);
            j[0] = wth[0]; j[2] = wth[2];
            for(j[1] = (wth[1]-2); j[1]<=(wth[1]+2); j[1]++){
                int I = ind(j,lth);
                d2[I] = DERIV02(Q[ind2(nu,k,NU,K)], j, lth, dx);
            }
            for(j[1] = (wth[1]-1); j[1]<=(wth[1]+1); j[1]++){
                int I = ind(j,lth);
                d4[I] = DERIV02(d2, j, lth, dx);
            }
            int I = ind(wth,lth);
            double d1 = DERIV01(Q[ind2(nu,k,NU,K)], wth, lth, dx);
            double d3 = DERIV01(d2, wth, lth, dx);
            double d5 = DERIV01(d4, wth, lth, dx);
            double d6 = DERIV02(d4, wth, lth, dx);
            Z[ind2(nu,1,NU,MU)] = d1 - ((dx[1]*dx[1])/6.0)*d3*kappa1[(nu-1)] +((dx[1]*dx[1]*dx[1]*dx[1])/30.0)*d5*kappa2[(nu-1)];
            Z[ind2(nu,2,NU,MU)] = d2[I] - ((dx[1]*dx[1])/12.0)*d4[I]*kappa1[(nu-1)] +((dx[1]*dx[1]*dx[1]*dx[1])/90.0)*d6*kappa2[(nu-1)];
            Z[ind2(nu,3,NU,MU)] = d3 - ((dx[1]*dx[1])/4.0)*d5*kappa2[(nu-1)];
            Z[ind2(nu,4,NU,MU)] = d4[I] - ((dx[1]*dx[1])/6.0)*d6*kappa2[(nu-1)];
            Z[ind2(nu,5,NU,MU)] = d5;
            Z[ind2(nu,6,NU,MU)] = d6;
        }
}

void LagrangeDeriv_2D_6_1(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    
    int p = 3;
    int wth[3] = {2*p,2*p,0}; wth[axis] = p;
    int lth[3] = {(2*wth[0]+1), (2*wth[1]+1), 1};
    
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 1-axis; // other axis 1
    int v1 = 2; // other axis 2
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d20[Lth], d02[Lth], d22[Lth];
    double Q11_20[Lth], Q11_02[Lth];
    double d40[Lth], d04[Lth];
    
    int i[3]; i[2] = 0;
    for(i[1] = 0; i[1]<lth[1]; i[1]++){
        for(i[0] = 0; i[0]<=lth[0]; i[0]++){
            Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)];
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
        for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
            int I = ind(i,lth);
            d20[I] = DERIV20(Q[ind2(0,1,NU,K)],i,lth,dx);
            d02[I] = DERIV02(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            double d10 = DERIV10(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d01 = DERIV01(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d11 = DERIV11(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d20[I])+coef[4][cI]*(d11)+coef[1][cI]*(d02[I])+coef[2][cI]*(d10)+coef[3][cI]*(d01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
        for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
            int I = ind(i,lth);
            d40[I] = DERIV20(d20, i, lth, dx);
            d04[I] = DERIV02(d02, i, lth, dx);
            d22[I] = DERIV20(d02, i, lth, dx);
            
            double d13 = DERIV11(d02, i, lth, dx);
            double d31 = DERIV11(d20, i, lth, dx);
            double d03 = DERIV01(d02, i, lth, dx);
            double d30 = DERIV10(d20, i, lth, dx);
            
            Q11_20[I] = DERIV20(Q[ind2(1,1,NU,K)],i, lth, dx);
            Q11_02[I] = DERIV02(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_11 = DERIV11(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_10 = DERIV10(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_01 = DERIV01(Q[ind2(1,1,NU,K)],i, lth, dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d40[I])+coef[4][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d31)+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d04[I])+coef[2][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d30)+coef[3][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d03);
            Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_20[I])+coef[4][cI]*(Q11_11)+coef[1][cI]*(Q11_02[I])+coef[2][cI]*(Q11_10)+coef[3][cI]*(Q11_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 3; i[1]<(lth[1]-3); i[1]++){
        for(i[0] = 3; i[0]<(lth[0]-3); i[0]++){
            int I = ind(i,lth);
            
            double d60 = DERIV20(d40, i, lth, dx);
            double d06 = DERIV02(d04, i, lth, dx);
            
            double d33 = DERIV11(d22, i, lth, dx);
            double d50 = DERIV10(d40, i, lth, dx);
            double d05 = DERIV01(d04, i, lth, dx);
            double d51 = DERIV11(d40, i, lth, dx);
            double d15 = DERIV11(d04, i, lth, dx);
            
            double Q11_40 = DERIV20(Q11_20,i, lth, dx);
            double Q11_04 = DERIV02(Q11_02,i, lth, dx);
            double Q21_20 = DERIV20(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q21_02 = DERIV02(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q12_20 = DERIV20(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q12_02 = DERIV02(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q11_30 = DERIV10(Q11_20,i, lth, dx);
            double Q11_03 = DERIV01(Q11_02,i, lth, dx);
            double Q11_31 = DERIV11(Q11_20,i, lth, dx);
            double Q11_13 = DERIV11(Q11_02,i, lth, dx);
            
            double Q21_10 = DERIV10(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q21_01 = DERIV01(Q[ind2(2,1,NU,K)],i, lth, dx);
            double Q21_11 = DERIV11(Q[ind2(2,1,NU,K)],i, lth, dx);
            
            double Q12_10 = DERIV10(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q12_01 = DERIV01(Q[ind2(1,2,NU,K)],i, lth, dx);
            double Q12_11 = DERIV11(Q[ind2(1,2,NU,K)],i, lth, dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,3,NU,K)][I] = Q[ind2(1,2,NU,K)][I]+coef[0][cI]*((1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d60)+coef[4][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d15+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*d33+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d51)+coef[1][cI]*((1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d06)+coef[2][cI]*((1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d50)+coef[3][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d05);
            Q[ind2(2,2,NU,K)][I] = coef[0][cI]*(Q12_20+((-1.0/12.0)*(dx[0]*dx[0])*Q11_40))+coef[4][cI]*(Q12_11+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q11_13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_31))+coef[1][cI]*(Q12_02+((-1.0/12.0)*(dx[1]*dx[1])*Q11_04))+coef[2][cI]*(Q12_10+((-1.0/6.0)*(dx[0]*dx[0])*Q11_30))+coef[3][cI]*(Q12_01+((-1.0/6.0)*(dx[1]*dx[1])*Q11_03));
            Q[ind2(3,1,NU,K)][I] = coef[0][cI]*(Q21_20)+coef[4][cI]*(Q21_11)+coef[1][cI]*(Q21_02)+coef[2][cI]*(Q21_10)+coef[3][cI]*(Q21_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind2(0,0,NU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind2(1,0,NU,MU)] = Q[ind2(1,3,NU,K)][I];
    Z[ind2(2,0,NU,MU)] = Q[ind2(2,2,NU,K)][I];
    Z[ind2(3,0,NU,MU)] = Q[ind2(3,1,NU,K)][I];
    
    double d2[Lth];
    double d4[Lth];
    int j[3] = {wth[0],wth[1],wth[2]};
    
    int nu = 0;
    double d1 = DERIV10(Q[ind2(nu,1,NU,K)], wth, lth, dx);
    double d3 = DERIV10(d20, wth, lth, dx);
    double d5 = DERIV10(d40, wth, lth, dx);
    double d6 = DERIV20(d40, wth, lth, dx);
    Z[ind2(nu,1,NU,MU)] = d1 - ((dx[0]*dx[0])/6.0)*d3 + ((dx[0]*dx[0]*dx[0]*dx[0])/30.0)*d5;
    Z[ind2(nu,2,NU,MU)] = d20[I] - ((dx[0]*dx[0])/12.0)*d40[I] + ((dx[0]*dx[0]*dx[0]*dx[0])/90.0)*d6;
    Z[ind2(nu,3,NU,MU)] = d3 - ((dx[0]*dx[0])/4.0)*d5;
    Z[ind2(nu,4,NU,MU)] = d40[I] - ((dx[0]*dx[0])/6.0)*d6;
    Z[ind2(nu,5,NU,MU)] = d5;
    Z[ind2(nu,6,NU,MU)] = d6;
    
    int kappa1[] = {1, 1, 0};
    int kappa2[] = {1, 0, 0};
    
    for(int nu = 1; nu<=p; nu++){
            int k = (p+1-nu);
            j[axis] = wth[axis]; j[2] = wth[2];
            for(j[v0] = (wth[v0]-2); j[v0]<=(wth[v0]+2); j[v0]++){
                int I = ind(j,lth);
                d2[I] = DERIV20(Q[ind2(nu,k,NU,K)], j, lth, dx);
            }
            for(j[v0] = (wth[v0]-1); j[v0]<=(wth[v0]+1); j[v0]++){
                int I = ind(j,lth);
                d4[I] = DERIV20(d2, j, lth, dx);
            }
            int I = ind(wth,lth);
            double d1 = DERIV10(Q[ind2(nu,k,NU,K)], wth, lth, dx);
            double d3 = DERIV10(d2, wth, lth, dx);
            double d5 = DERIV10(d4, wth, lth, dx);
            double d6 = DERIV20(d4, wth, lth, dx);
            Z[ind2(nu,1,NU,MU)] = d1 - ((dx[0]*dx[0])/6.0)*d3*kappa1[(nu-1)] +((dx[0]*dx[0]*dx[0]*dx[0])/30.0)*d5*kappa2[(nu-1)];
            Z[ind2(nu,2,NU,MU)] = d2[I] - ((dx[0]*dx[0])/12.0)*d4[I]*kappa1[(nu-1)] +((dx[0]*dx[0]*dx[0]*dx[0])/90.0)*d6*kappa2[(nu-1)];
            Z[ind2(nu,3,NU,MU)] = d3 - ((dx[0]*dx[0])/4.0)*d5*kappa2[(nu-1)];
            Z[ind2(nu,4,NU,MU)] = d4[I] - ((dx[0]*dx[0])/6.0)*d6*kappa2[(nu-1)];
            Z[ind2(nu,5,NU,MU)] = d5;
            Z[ind2(nu,6,NU,MU)] = d6;
        }
}

void LagrangeDeriv_2D_4_0(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    
    int p = 2;
    int wth[3] = {p,2*p,0}, lth[3] = {(2*p+1), (4*p+1), 1};
    
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 1; // other axis 1
    int v1 = 2; // other axis 2
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d20[Lth], d02[Lth];
    double d40[Lth], d04[Lth];
    
    int i[3]; i[2] = 0;
    for(i[1] = 0; i[1]<lth[1]; i[1]++){
        for(i[0] = 0; i[0]<=lth[0]; i[0]++){
            Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)];
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
        for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
            int I = ind(i,lth);
            d20[I] = DERIV20(Q[ind2(0,1,NU,K)],i,lth,dx);
            d02[I] = DERIV02(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            double d10 = DERIV10(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d01 = DERIV01(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d11 = DERIV11(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d20[I])+coef[4][cI]*(d11)+coef[1][cI]*(d02[I])+coef[2][cI]*(d10)+coef[3][cI]*(d01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
        for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
            int I = ind(i,lth);
            d40[I] = DERIV20(d20, i, lth, dx);
            d04[I] = DERIV02(d02, i, lth, dx);
            
            double d13 = DERIV11(d02, i, lth, dx);
            double d31 = DERIV11(d20, i, lth, dx);
            double d03 = DERIV01(d02, i, lth, dx);
            double d30 = DERIV10(d20, i, lth, dx);
            
            double Q11_20 = DERIV20(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_02 = DERIV02(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_11 = DERIV11(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_10 = DERIV10(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_01 = DERIV01(Q[ind2(1,1,NU,K)],i, lth, dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d40[I])+coef[4][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d31)+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d04[I])+coef[2][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d30)+coef[3][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d03);
            Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_20)+coef[4][cI]*(Q11_11)+coef[1][cI]*(Q11_02)+coef[2][cI]*(Q11_10)+coef[3][cI]*(Q11_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind2(0,0,NU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind2(1,0,NU,MU)] = Q[ind2(1,2,NU,K)][I];
    Z[ind2(2,0,NU,MU)] = Q[ind2(2,1,NU,K)][I];
    
    double d2[lth[0]*lth[1]*lth[2]];
    int j[3] = {wth[0],wth[1],wth[2]};
    
    int nu = 0;
    double d01 = DERIV01(Q[ind2(nu,1,NU,K)], wth, lth, dx);
    double d03 = DERIV01(d02, wth, lth, dx);
    Z[ind2(nu,1,NU,MU)] = d01 - ((dx[1]*dx[1])/6.0)*d03;
    Z[ind2(nu,2,NU,MU)] = d02[I] - ((dx[1]*dx[1])/12.0)*d04[I];
    Z[ind2(nu,3,NU,MU)] = d03;
    Z[ind2(nu,4,NU,MU)] = d04[I];
    
    for(int nu = 1; nu<=2; nu++){
        int k = (p+1-nu);
        j[0] = wth[0]; j[2] = wth[2];
        for(j[1] = (wth[1]-1); j[1]<=(wth[1]+1); j[1]++){
            int I = ind(j,lth);
            d2[I] = DERIV02(Q[ind2(nu,k,NU,K)], j, lth, dx);
        }
        int I = ind(wth,lth);
        double d1 = DERIV01(Q[ind2(nu,k,NU,K)], wth, lth, dx);
        double d3 = DERIV01(d2, wth, lth, dx);
        double d4 = DERIV02(d2, wth, lth, dx);
        Z[ind2(nu,1,NU,MU)] = d1 - ((dx[1]*dx[1])/6.0)*d3*(k-1);
        Z[ind2(nu,2,NU,MU)] = d2[I] - ((dx[1]*dx[1])/12.0)*d4*(k-1);
        Z[ind2(nu,3,NU,MU)] = d3;
        Z[ind2(nu,4,NU,MU)] = d4;
    }
}

void LagrangeDeriv_2D_2_0(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    
    int p = 1, dim = 2;
    int wth[3] = {p,2*p,0}, lth[3] = {(2*p+1), (4*p+1), 1};
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 1;
    int v1 = 2;
    int cI = 0;
    int cInd[3];
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double V[Lth];
    double V11[Lth];
    double d20[Lth], d02[Lth];
    
    int i[3]; i[2] = 0;
    for(i[1] = 0; i[1]<lth[1]; i[1]++){
        for(i[0] = 0; i[0]<=lth[0]; i[0]++){
            V[ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)];
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
        for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
            int I = ind(i,lth);
            d20[I] = DERIV20(V,i,lth,dx);
            d02[I] = DERIV02(V,i,lth,dx);
            
            double d10 = DERIV10(V,i,lth,dx);
            double d01 = DERIV01(V,i,lth,dx);
            double d11 = DERIV11(V,i,lth,dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            V11[I] = coef[0][cI]*(d20[I])+coef[4][cI]*(d11)+coef[1][cI]*(d02[I])+coef[2][cI]*(d10)+coef[3][cI]*(d01);
        } // end i[0] loop
    }// end of i[1] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    int MU = 2*p+1;
    
    int I = ind(wth,lth);
    Z[ind2(0,0,NU,MU)] = V[I];
    Z[ind2(1,0,NU,MU)] = V11[I];
    Z[ind2(0,1,NU,MU)] = DERIV01(V, wth, lth, dx);
    Z[ind2(0,2,NU,MU)] = DERIV02(V, wth, lth, dx);
    Z[ind2(1,1,NU,MU)] = DERIV01(V11, wth, lth, dx);
    Z[ind2(1,2,NU,MU)] = DERIV02(V11, wth, lth, dx);
}

void LagrangeDeriv_2D_2_1(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    
    int p = 1, dim = 2;
    int wth[3] = {2*p,p,0}, lth[3] = {(4*p+1), (2*p+1), 1};
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 0;
    int v1 = 2;
    int cI = 0;
    int cInd[3];
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double V[Lth];
    double V11[Lth];
    double d20[Lth], d02[Lth];
    
    int i[3]; i[2] = 0;
    for(i[1] = 0; i[1]<lth[1]; i[1]++){
        for(i[0] = 0; i[0]<=lth[0]; i[0]++){
            V[ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)];
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
        for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
            int I = ind(i,lth);
            d20[I] = DERIV20(V,i,lth,dx);
            d02[I] = DERIV02(V,i,lth,dx);
            
            double d10 = DERIV10(V,i,lth,dx);
            double d01 = DERIV01(V,i,lth,dx);
            double d11 = DERIV11(V,i,lth,dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            V11[I] = coef[0][cI]*(d20[I])+coef[4][cI]*(d11)+coef[1][cI]*(d02[I])+coef[2][cI]*(d10)+coef[3][cI]*(d01);
        } // end i[0] loop
    }// end of i[1] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    int MU = 2*p+1;
    
    int I = ind(wth,lth);
    Z[ind2(0,0,NU,MU)] = V[I];
    Z[ind2(1,0,NU,MU)] = V11[I];
    Z[ind2(0,1,NU,MU)] = DERIV10(V, wth, lth, dx);
    Z[ind2(0,2,NU,MU)] = DERIV20(V, wth, lth, dx);
    Z[ind2(1,1,NU,MU)] = DERIV10(V11, wth, lth, dx);
    Z[ind2(1,2,NU,MU)] = DERIV20(V11, wth, lth, dx);
}

void LagrangeDeriv_2D_4_1(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    
    int p = 2;
    int wth[3] = {2*p,p,0}, lth[3] = {(4*p+1), (2*p+1), 1};
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 0;
    int v1 = 2;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d20[Lth], d02[Lth];
    double d40[Lth], d04[Lth];
    
    int i[3]; i[2] = 0;
    for(i[1] = 0; i[1]<lth[1]; i[1]++){
        for(i[0] = 0; i[0]<=lth[0]; i[0]++){
            Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)];
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
        for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
            int I = ind(i,lth);
            d20[I] = DERIV20(Q[ind2(0,1,NU,K)],i,lth,dx);
            d02[I] = DERIV02(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            double d10 = DERIV10(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d01 = DERIV01(Q[ind2(0,1,NU,K)],i,lth,dx);
            double d11 = DERIV11(Q[ind2(0,1,NU,K)],i,lth,dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d20[I])+coef[4][cI]*(d11)+coef[1][cI]*(d02[I])+coef[2][cI]*(d10)+coef[3][cI]*(d01);
        } // end i[0] loop
    }// end of i[1] loop
    
    for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
        for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
            int I = ind(i,lth);
            d40[I] = DERIV20(d20, i, lth, dx);
            d04[I] = DERIV02(d02, i, lth, dx);
            
            double d13 = DERIV11(d02, i, lth, dx);
            double d31 = DERIV11(d20, i, lth, dx);
            double d03 = DERIV01(d02, i, lth, dx);
            double d30 = DERIV10(d20, i, lth, dx);
            
            double Q11_20 = DERIV20(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_02 = DERIV02(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_11 = DERIV11(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_10 = DERIV10(Q[ind2(1,1,NU,K)],i, lth, dx);
            double Q11_01 = DERIV01(Q[ind2(1,1,NU,K)],i, lth, dx);
            
            if(!cstCoef){
                cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                cInd[v1] = 0;
                cI = ind(cInd, coefLth);
            }
            
            Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d40[I])+coef[4][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d31)+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d04[I])+coef[2][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d30)+coef[3][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d03);
            Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_20)+coef[4][cI]*(Q11_11)+coef[1][cI]*(Q11_02)+coef[2][cI]*(Q11_10)+coef[3][cI]*(Q11_01);
        } // end i[0] loop
    }// end of i[1] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind2(0,0,NU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind2(1,0,NU,MU)] = Q[ind2(1,2,NU,K)][I];
    Z[ind2(2,0,NU,MU)] = Q[ind2(2,1,NU,K)][I];
    
    double d2[lth[0]*lth[1]*lth[2]];
    int j[3] = {wth[0],wth[1],wth[2]};
    
    int nu = 0;
    double d10 = DERIV10(Q[ind2(nu,1,NU,K)], wth, lth, dx);
    double d30 = DERIV10(d20, wth, lth, dx);
    Z[ind2(nu,1,NU,MU)] = d10 - ((dx[0]*dx[0])/6.0)*d30;
    Z[ind2(nu,2,NU,MU)] = d20[I] - ((dx[0]*dx[0])/12.0)*d40[I];
    Z[ind2(nu,3,NU,MU)] = d30;
    Z[ind2(nu,4,NU,MU)] = d40[I];
    
    for(int nu = 1; nu<=2; nu++){
        int k = (p+1-nu);
        for(j[0] = (wth[0]-1); j[0]<=(wth[0]+1); j[0]++){
            int I = ind(j,lth);
            d2[I] = DERIV20(Q[ind2(nu,k,NU,K)], j, lth, dx);
        }
        int I = ind(wth,lth);
        double d1 = DERIV10(Q[ind2(nu,k,NU,K)], wth, lth, dx);
        double d3 = DERIV10(d2, wth, lth, dx);
        double d4 = DERIV20(d2, wth, lth, dx);
        Z[ind2(nu,1,NU,MU)] = d1 - ((dx[0]*dx[0])/6.0)*d3*(k-1);
        Z[ind2(nu,2,NU,MU)] = d2[I] - ((dx[0]*dx[0])/12.0)*d4*(k-1);
        Z[ind2(nu,3,NU,MU)] = d3;
        Z[ind2(nu,4,NU,MU)] = d4;
    }
}

void LagrangeDeriv_3D_2_0(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 1;
    int wth[3] = {(2*p),(2*p),(2*p)}; wth[axis] = p;
    int lth[3]= {(2*wth[0]+1),(2*wth[1]+1),(2*wth[2]+1)};
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 1;
    int v1 = 2;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d200[Lth], d020[Lth], d002[Lth];
    
    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_2_0.C"
    
}

void LagrangeDeriv_3D_2_1(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 1;
    int wth[3] = {(2*p),(2*p),(2*p)}; wth[axis] = p;
    int lth[3]= {(2*wth[0]+1),(2*wth[1]+1),(2*wth[2]+1)};
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 0;
    int v1 = 2;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d200[Lth], d020[Lth], d002[Lth];
    
    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_2_1.C"
    
}

void LagrangeDeriv_3D_2_2(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 1;
    int wth[3] = {(2*p),(2*p),(2*p)}; wth[axis] = p;
    int lth[3]= {(2*wth[0]+1),(2*wth[1]+1),(2*wth[2]+1)};
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 0;
    int v1 = 1;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d200[Lth], d020[Lth], d002[Lth];
    
    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_2_2.C"
    
}


void LagrangeDeriv_3D_4_0(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 2;
    int wth[3] = {p,(2*p),(2*p)}, lth[3]= {(2*p+1),(4*p+1),(4*p+1)};
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 1;
    int v1 = 2;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d200[Lth], d020[Lth], d002[Lth];
    
    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 2; i[2]<(lth[2]-2); i[2]++){
        for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
            for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
                int I = ind(i,lth);
                
                double d004 = DERIV002(d002,i,lth,dx);
                double d040 = DERIV020(d020,i,lth,dx);
                double d400 = DERIV200(d200,i,lth,dx);
                
                double d003 = DERIV001(d002,i,lth,dx);
                double d030 = DERIV010(d020,i,lth,dx);
                double d300 = DERIV100(d200,i,lth,dx);
                double d013 = DERIV011(d002,i,lth,dx);
                double d031 = DERIV011(d020,i,lth,dx);
                double d103 = DERIV101(d002,i,lth,dx);
                double d130 = DERIV110(d020,i,lth,dx);
                double d301 = DERIV101(d200,i,lth,dx);
                double d310 = DERIV110(d200,i,lth,dx);
                
                double Q11_002 = DERIV002(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_020 = DERIV020(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_200 = DERIV200(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_001 = DERIV001(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_010 = DERIV010(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_100 = DERIV100(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_011 = DERIV011(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_101 = DERIV101(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_110 = DERIV110(Q[ind2(1,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                
                Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d400)+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d040)+coef[2][cI]*((-1.0/12.0)*(dx[2]*dx[2])*d004)+coef[6][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d310)+coef[7][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d301)+coef[8][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*d031)+coef[3][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d300)+coef[4][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d030)+coef[5][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d003);
                Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_200)+coef[1][cI]*(Q11_020)+coef[2][cI]*(Q11_002)+coef[6][cI]*(Q11_110)+coef[7][cI]*(Q11_101)+coef[8][cI]*(Q11_011)+coef[3][cI]*(Q11_100)+coef[4][cI]*(Q11_010)+coef[5][cI]*(Q11_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,2,NU,K)][I];
    Z[ind3(2,0,0,NU,MU,MU)] = Q[ind2(2,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_4_0.C"
    
}

void LagrangeDeriv_3D_4_1(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 2;
    int wth[3] = {(2*p),p,(2*p)}, lth[3]= {(4*p+1),(2*p+1),(4*p+1)};
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 0;
    int v1 = 2;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d200[Lth], d020[Lth], d002[Lth];
    
    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 2; i[2]<(lth[2]-2); i[2]++){
        for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
            for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
                int I = ind(i,lth);
                
                double d004 = DERIV002(d002,i,lth,dx);
                double d040 = DERIV020(d020,i,lth,dx);
                double d400 = DERIV200(d200,i,lth,dx);
                
                double d003 = DERIV001(d002,i,lth,dx);
                double d030 = DERIV010(d020,i,lth,dx);
                double d300 = DERIV100(d200,i,lth,dx);
                double d013 = DERIV011(d002,i,lth,dx);
                double d031 = DERIV011(d020,i,lth,dx);
                double d103 = DERIV101(d002,i,lth,dx);
                double d130 = DERIV110(d020,i,lth,dx);
                double d301 = DERIV101(d200,i,lth,dx);
                double d310 = DERIV110(d200,i,lth,dx);
                
                double Q11_002 = DERIV002(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_020 = DERIV020(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_200 = DERIV200(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_001 = DERIV001(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_010 = DERIV010(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_100 = DERIV100(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_011 = DERIV011(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_101 = DERIV101(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_110 = DERIV110(Q[ind2(1,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                
                Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d400)+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d040)+coef[2][cI]*((-1.0/12.0)*(dx[2]*dx[2])*d004)+coef[6][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d310)+coef[7][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d301)+coef[8][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*d031)+coef[3][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d300)+coef[4][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d030)+coef[5][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d003);
                Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_200)+coef[1][cI]*(Q11_020)+coef[2][cI]*(Q11_002)+coef[6][cI]*(Q11_110)+coef[7][cI]*(Q11_101)+coef[8][cI]*(Q11_011)+coef[3][cI]*(Q11_100)+coef[4][cI]*(Q11_010)+coef[5][cI]*(Q11_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,2,NU,K)][I];
    Z[ind3(2,0,0,NU,MU,MU)] = Q[ind2(2,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_4_1.C"
}

void LagrangeDeriv_3D_4_2(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 2;
    int wth[3] = {(2*p),(2*p),p}, lth[3]= {(4*p+1),(4*p+1),(2*p+1)};
    int L_center = LagrangeData_center;
    int L_wth    = 2*L_center + 1;
    int v0 = 0;
    int v1 = 1;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d200[Lth], d020[Lth], d002[Lth];
    
    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 2; i[2]<(lth[2]-2); i[2]++){
        for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
            for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
                int I = ind(i,lth);
                
                double d004 = DERIV002(d002,i,lth,dx);
                double d040 = DERIV020(d020,i,lth,dx);
                double d400 = DERIV200(d200,i,lth,dx);
                
                double d003 = DERIV001(d002,i,lth,dx);
                double d030 = DERIV010(d020,i,lth,dx);
                double d300 = DERIV100(d200,i,lth,dx);
                double d013 = DERIV011(d002,i,lth,dx);
                double d031 = DERIV011(d020,i,lth,dx);
                double d103 = DERIV101(d002,i,lth,dx);
                double d130 = DERIV110(d020,i,lth,dx);
                double d301 = DERIV101(d200,i,lth,dx);
                double d310 = DERIV110(d200,i,lth,dx);
                
                double Q11_002 = DERIV002(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_020 = DERIV020(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_200 = DERIV200(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_001 = DERIV001(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_010 = DERIV010(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_100 = DERIV100(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_011 = DERIV011(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_101 = DERIV101(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_110 = DERIV110(Q[ind2(1,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                
                Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d400)+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d040)+coef[2][cI]*((-1.0/12.0)*(dx[2]*dx[2])*d004)+coef[6][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d310)+coef[7][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d301)+coef[8][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*d031)+coef[3][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d300)+coef[4][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d030)+coef[5][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d003);
                Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_200)+coef[1][cI]*(Q11_020)+coef[2][cI]*(Q11_002)+coef[6][cI]*(Q11_110)+coef[7][cI]*(Q11_101)+coef[8][cI]*(Q11_011)+coef[3][cI]*(Q11_100)+coef[4][cI]*(Q11_010)+coef[5][cI]*(Q11_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,2,NU,K)][I];
    Z[ind3(2,0,0,NU,MU,MU)] = Q[ind2(2,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_4_2.C"
}

void LagrangeDeriv_3D_6_0(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 3;
    int wth[3] = {(2*p),(2*p),(2*p)}; wth[axis] = p;
    int lth[3]= {(2*wth[0]+1),(2*wth[1]+1),(2*wth[2]+1)};
    int L_center = LagrangeData_center;
//    int L_wth    = 2*L_center + 1;
    int v0 = 1;
    int v1 = 2;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d200[Lth], d020[Lth], d002[Lth];
    double d004[Lth], d022[Lth], d040[Lth], d202[Lth], d220[Lth], d400[Lth];
    double Q11_002[Lth], Q11_020[Lth], Q11_200[Lth];

    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] - p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 2; i[2]<(lth[2]-2); i[2]++){
        for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
            for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
                int I = ind(i,lth);
                
                d004[I] = DERIV002(d002,i,lth,dx);
                d022[I] = DERIV002(d020,i,lth,dx);
                d040[I] = DERIV020(d020,i,lth,dx);
                d202[I] = DERIV200(d002,i,lth,dx);
                d220[I] = DERIV200(d020,i,lth,dx);
                d400[I] = DERIV200(d200,i,lth,dx);
                
                double d003 = DERIV001(d002,i,lth,dx);
                double d030 = DERIV010(d020,i,lth,dx);
                double d300 = DERIV100(d200,i,lth,dx);
                double d013 = DERIV011(d002,i,lth,dx);
                double d031 = DERIV011(d020,i,lth,dx);
                double d103 = DERIV101(d002,i,lth,dx);
                double d130 = DERIV110(d020,i,lth,dx);
                double d301 = DERIV101(d200,i,lth,dx);
                double d310 = DERIV110(d200,i,lth,dx);
                
                Q11_002[I] = DERIV002(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_020[I] = DERIV020(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_200[I] = DERIV200(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_001 = DERIV001(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_010 = DERIV010(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_100 = DERIV100(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_011 = DERIV011(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_101 = DERIV101(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_110 = DERIV110(Q[ind2(1,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] - p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                
                Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d400[I])+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d040[I])+coef[2][cI]*((-1.0/12.0)*(dx[2]*dx[2])*d004[I])+coef[6][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d310)+coef[7][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d301)+coef[8][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*d031)+coef[3][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d300)+coef[4][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d030)+coef[5][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d003);
                Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_200[I])+coef[1][cI]*(Q11_020[I])+coef[2][cI]*(Q11_002[I])+coef[6][cI]*(Q11_110)+coef[7][cI]*(Q11_101)+coef[8][cI]*(Q11_011)+coef[3][cI]*(Q11_100)+coef[4][cI]*(Q11_010)+coef[5][cI]*(Q11_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 3; i[2]<(lth[2]-3); i[2]++){
        for(i[1] = 3; i[1]<(lth[1]-3); i[1]++){
            for(i[0] = 3; i[0]<(lth[0]-3); i[0]++){
                int I = ind(i,lth);
                
                double d006 = DERIV002(d004,i,lth,dx);
                double d060 = DERIV020(d040,i,lth,dx);
                double d600 = DERIV200(d400,i,lth,dx);
                
                double d005 = DERIV001(d004,i,lth,dx);
                double d050 = DERIV010(d040,i,lth,dx);
                double d500 = DERIV100(d400,i,lth,dx);
                double d015 = DERIV011(d004,i,lth,dx);
                double d033 = DERIV011(d022,i,lth,dx);
                double d051 = DERIV011(d040,i,lth,dx);
                double d105 = DERIV101(d004,i,lth,dx);
                double d150 = DERIV110(d040,i,lth,dx);
                double d303 = DERIV101(d202,i,lth,dx);
                double d330 = DERIV110(d220,i,lth,dx);
                double d501 = DERIV101(d400,i,lth,dx);
                double d510 = DERIV110(d400,i,lth,dx);
                
                double Q11_004 = DERIV002(Q11_002,i,lth,dx);
                double Q11_040 = DERIV020(Q11_020,i,lth,dx);
                double Q11_400 = DERIV200(Q11_200,i,lth,dx);
                double Q21_002 = DERIV002(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_020 = DERIV020(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_200 = DERIV200(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q12_002 = DERIV002(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_020 = DERIV020(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_200 = DERIV200(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q11_003 = DERIV001(Q11_002,i,lth,dx);
                double Q11_030 = DERIV010(Q11_020,i,lth,dx);
                double Q11_300 = DERIV100(Q11_200,i,lth,dx);
                double Q11_013 = DERIV011(Q11_002,i,lth,dx);
                double Q11_031 = DERIV011(Q11_020,i,lth,dx);
                double Q11_103 = DERIV101(Q11_002,i,lth,dx);
                double Q11_130 = DERIV110(Q11_020,i,lth,dx);
                double Q11_301 = DERIV101(Q11_200,i,lth,dx);
                double Q11_310 = DERIV110(Q11_200,i,lth,dx);
                
                double Q21_001 = DERIV001(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_010 = DERIV010(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_100 = DERIV100(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_011 = DERIV011(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_101 = DERIV101(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_110 = DERIV110(Q[ind2(2,1,NU,K)],i,lth,dx);
                
                double Q12_001 = DERIV001(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_010 = DERIV010(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_100 = DERIV100(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_011 = DERIV011(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_101 = DERIV101(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_110 = DERIV110(Q[ind2(1,2,NU,K)],i,lth,dx);

                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                
                Q[ind2(1,3,NU,K)][I] = Q[ind2(1,2,NU,K)][I]+coef[0][cI]*((1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d600)+coef[1][cI]*((1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d060)+coef[2][cI]*((1.0/90.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d006)+coef[6][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d150+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*d330+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d510)+coef[7][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d105+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[2]*dx[2])*d303+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d501)+coef[8][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d015+(-1.0/6.0)*(-1.0/6.0)*(dx[1]*dx[1])*(dx[2]*dx[2])*d033+(1.0/30.0)*1.0*(dx[1]*dx[1]*dx[1]*dx[1])*(1.0)*d051)+coef[3][cI]*((1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d500)+coef[4][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d050)+coef[5][cI]*((1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d005);
                Q[ind2(2,2,NU,K)][I] = coef[0][cI]*(Q12_200+((-1.0/12.0)*(dx[0]*dx[0])*Q11_400))+coef[1][cI]*(Q12_020+((-1.0/12.0)*(dx[1]*dx[1])*Q11_040))+coef[2][cI]*(Q12_002+((-1.0/12.0)*(dx[2]*dx[2])*Q11_004))+coef[6][cI]*(Q12_110+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q11_130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_310))+coef[7][cI]*(Q12_101+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_301))+coef[8][cI]*(Q12_011+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*Q11_031))+coef[3][cI]*(Q12_100+((-1.0/6.0)*(dx[0]*dx[0])*Q11_300))+coef[4][cI]*(Q12_010+((-1.0/6.0)*(dx[1]*dx[1])*Q11_030))+coef[5][cI]*(Q12_001+((-1.0/6.0)*(dx[2]*dx[2])*Q11_003));
                Q[ind2(3,1,NU,K)][I] = coef[0][cI]*(Q21_200)+coef[1][cI]*(Q21_020)+coef[2][cI]*(Q21_002)+coef[6][cI]*(Q21_110)+coef[7][cI]*(Q21_101)+coef[8][cI]*(Q21_011)+coef[3][cI]*(Q21_100)+coef[4][cI]*(Q21_010)+coef[5][cI]*(Q21_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,3,NU,K)][I];
    Z[ind3(2,0,0,NU,MU,MU)] = Q[ind2(2,2,NU,K)][I];
    Z[ind3(3,0,0,NU,MU,MU)] = Q[ind2(3,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_6_0.C"
}

void LagrangeDeriv_3D_6_1(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 3;
    int wth[3] = {(2*p),(2*p),(2*p)}; wth[axis] = p;
    int lth[3]= {(2*wth[0]+1),(2*wth[1]+1),(2*wth[2]+1)};
    int L_center = LagrangeData_center;
//    int L_wth    = 2*L_center + 1;
    int v0 = 0;
    int v1 = 2;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d200[Lth], d020[Lth], d002[Lth];
    double d004[Lth], d022[Lth], d040[Lth], d202[Lth], d220[Lth], d400[Lth];
    double Q11_002[Lth], Q11_020[Lth], Q11_200[Lth];

    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 2; i[2]<(lth[2]-2); i[2]++){
        for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
            for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
                int I = ind(i,lth);
                
                d004[I] = DERIV002(d002,i,lth,dx);
                d022[I] = DERIV002(d020,i,lth,dx);
                d040[I] = DERIV020(d020,i,lth,dx);
                d202[I] = DERIV200(d002,i,lth,dx);
                d220[I] = DERIV200(d020,i,lth,dx);
                d400[I] = DERIV200(d200,i,lth,dx);
                
                double d003 = DERIV001(d002,i,lth,dx);
                double d030 = DERIV010(d020,i,lth,dx);
                double d300 = DERIV100(d200,i,lth,dx);
                double d013 = DERIV011(d002,i,lth,dx);
                double d031 = DERIV011(d020,i,lth,dx);
                double d103 = DERIV101(d002,i,lth,dx);
                double d130 = DERIV110(d020,i,lth,dx);
                double d301 = DERIV101(d200,i,lth,dx);
                double d310 = DERIV110(d200,i,lth,dx);
                
                Q11_002[I] = DERIV002(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_020[I] = DERIV020(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_200[I] = DERIV200(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_001 = DERIV001(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_010 = DERIV010(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_100 = DERIV100(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_011 = DERIV011(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_101 = DERIV101(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_110 = DERIV110(Q[ind2(1,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                
                Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d400[I])+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d040[I])+coef[2][cI]*((-1.0/12.0)*(dx[2]*dx[2])*d004[I])+coef[6][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d310)+coef[7][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d301)+coef[8][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*d031)+coef[3][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d300)+coef[4][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d030)+coef[5][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d003);
                Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_200[I])+coef[1][cI]*(Q11_020[I])+coef[2][cI]*(Q11_002[I])+coef[6][cI]*(Q11_110)+coef[7][cI]*(Q11_101)+coef[8][cI]*(Q11_011)+coef[3][cI]*(Q11_100)+coef[4][cI]*(Q11_010)+coef[5][cI]*(Q11_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 3; i[2]<(lth[2]-3); i[2]++){
        for(i[1] = 3; i[1]<(lth[1]-3); i[1]++){
            for(i[0] = 3; i[0]<(lth[0]-3); i[0]++){
                int I = ind(i,lth);
                
                double d006 = DERIV002(d004,i,lth,dx);
                double d060 = DERIV020(d040,i,lth,dx);
                double d600 = DERIV200(d400,i,lth,dx);
                
                double d005 = DERIV001(d004,i,lth,dx);
                double d050 = DERIV010(d040,i,lth,dx);
                double d500 = DERIV100(d400,i,lth,dx);
                double d015 = DERIV011(d004,i,lth,dx);
                double d033 = DERIV011(d022,i,lth,dx);
                double d051 = DERIV011(d040,i,lth,dx);
                double d105 = DERIV101(d004,i,lth,dx);
                double d150 = DERIV110(d040,i,lth,dx);
                double d303 = DERIV101(d202,i,lth,dx);
                double d330 = DERIV110(d220,i,lth,dx);
                double d501 = DERIV101(d400,i,lth,dx);
                double d510 = DERIV110(d400,i,lth,dx);
                
                double Q11_004 = DERIV002(Q11_002,i,lth,dx);
                double Q11_040 = DERIV020(Q11_020,i,lth,dx);
                double Q11_400 = DERIV200(Q11_200,i,lth,dx);
                double Q21_002 = DERIV002(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_020 = DERIV020(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_200 = DERIV200(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q12_002 = DERIV002(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_020 = DERIV020(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_200 = DERIV200(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q11_003 = DERIV001(Q11_002,i,lth,dx);
                double Q11_030 = DERIV010(Q11_020,i,lth,dx);
                double Q11_300 = DERIV100(Q11_200,i,lth,dx);
                double Q11_013 = DERIV011(Q11_002,i,lth,dx);
                double Q11_031 = DERIV011(Q11_020,i,lth,dx);
                double Q11_103 = DERIV101(Q11_002,i,lth,dx);
                double Q11_130 = DERIV110(Q11_020,i,lth,dx);
                double Q11_301 = DERIV101(Q11_200,i,lth,dx);
                double Q11_310 = DERIV110(Q11_200,i,lth,dx);
                
                double Q21_001 = DERIV001(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_010 = DERIV010(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_100 = DERIV100(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_011 = DERIV011(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_101 = DERIV101(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_110 = DERIV110(Q[ind2(2,1,NU,K)],i,lth,dx);
                
                double Q12_001 = DERIV001(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_010 = DERIV010(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_100 = DERIV100(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_011 = DERIV011(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_101 = DERIV101(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_110 = DERIV110(Q[ind2(1,2,NU,K)],i,lth,dx);

                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                
                Q[ind2(1,3,NU,K)][I] = Q[ind2(1,2,NU,K)][I]+coef[0][cI]*((1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d600)+coef[1][cI]*((1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d060)+coef[2][cI]*((1.0/90.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d006)+coef[6][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d150+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*d330+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d510)+coef[7][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d105+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[2]*dx[2])*d303+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d501)+coef[8][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d015+(-1.0/6.0)*(-1.0/6.0)*(dx[1]*dx[1])*(dx[2]*dx[2])*d033+(1.0/30.0)*1.0*(dx[1]*dx[1]*dx[1]*dx[1])*(1.0)*d051)+coef[3][cI]*((1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d500)+coef[4][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d050)+coef[5][cI]*((1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d005);
                Q[ind2(2,2,NU,K)][I] = coef[0][cI]*(Q12_200+((-1.0/12.0)*(dx[0]*dx[0])*Q11_400))+coef[1][cI]*(Q12_020+((-1.0/12.0)*(dx[1]*dx[1])*Q11_040))+coef[2][cI]*(Q12_002+((-1.0/12.0)*(dx[2]*dx[2])*Q11_004))+coef[6][cI]*(Q12_110+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q11_130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_310))+coef[7][cI]*(Q12_101+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_301))+coef[8][cI]*(Q12_011+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*Q11_031))+coef[3][cI]*(Q12_100+((-1.0/6.0)*(dx[0]*dx[0])*Q11_300))+coef[4][cI]*(Q12_010+((-1.0/6.0)*(dx[1]*dx[1])*Q11_030))+coef[5][cI]*(Q12_001+((-1.0/6.0)*(dx[2]*dx[2])*Q11_003));
                Q[ind2(3,1,NU,K)][I] = coef[0][cI]*(Q21_200)+coef[1][cI]*(Q21_020)+coef[2][cI]*(Q21_002)+coef[6][cI]*(Q21_110)+coef[7][cI]*(Q21_101)+coef[8][cI]*(Q21_011)+coef[3][cI]*(Q21_100)+coef[4][cI]*(Q21_010)+coef[5][cI]*(Q21_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,3,NU,K)][I];
    Z[ind3(2,0,0,NU,MU,MU)] = Q[ind2(2,2,NU,K)][I];
    Z[ind3(3,0,0,NU,MU,MU)] = Q[ind2(3,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_6_1.C"
}

void LagrangeDeriv_3D_6_2(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 3;
    int wth[3] = {(2*p),(2*p),(2*p)}; wth[axis] = p;
    int lth[3]= {(2*wth[0]+1),(2*wth[1]+1),(2*wth[2]+1)};
    int L_center = LagrangeData_center;
//    int L_wth    = 2*L_center + 1;
    int v0 = 0;
    int v1 = 1;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
    double Q[NU*K][Lth];
    double d200[Lth], d020[Lth], d002[Lth];
    double d004[Lth], d022[Lth], d040[Lth], d202[Lth], d220[Lth], d400[Lth];
    double Q11_002[Lth], Q11_020[Lth], Q11_200[Lth];

    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 2; i[2]<(lth[2]-2); i[2]++){
        for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
            for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
                int I = ind(i,lth);
                
                d004[I] = DERIV002(d002,i,lth,dx);
                d022[I] = DERIV002(d020,i,lth,dx);
                d040[I] = DERIV020(d020,i,lth,dx);
                d202[I] = DERIV200(d002,i,lth,dx);
                d220[I] = DERIV200(d020,i,lth,dx);
                d400[I] = DERIV200(d200,i,lth,dx);
                
                double d003 = DERIV001(d002,i,lth,dx);
                double d030 = DERIV010(d020,i,lth,dx);
                double d300 = DERIV100(d200,i,lth,dx);
                double d013 = DERIV011(d002,i,lth,dx);
                double d031 = DERIV011(d020,i,lth,dx);
                double d103 = DERIV101(d002,i,lth,dx);
                double d130 = DERIV110(d020,i,lth,dx);
                double d301 = DERIV101(d200,i,lth,dx);
                double d310 = DERIV110(d200,i,lth,dx);
                
                Q11_002[I] = DERIV002(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_020[I] = DERIV020(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_200[I] = DERIV200(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_001 = DERIV001(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_010 = DERIV010(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_100 = DERIV100(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_011 = DERIV011(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_101 = DERIV101(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_110 = DERIV110(Q[ind2(1,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                
                Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d400[I])+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d040[I])+coef[2][cI]*((-1.0/12.0)*(dx[2]*dx[2])*d004[I])+coef[6][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d310)+coef[7][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d301)+coef[8][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*d031)+coef[3][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d300)+coef[4][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d030)+coef[5][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d003);
                Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_200[I])+coef[1][cI]*(Q11_020[I])+coef[2][cI]*(Q11_002[I])+coef[6][cI]*(Q11_110)+coef[7][cI]*(Q11_101)+coef[8][cI]*(Q11_011)+coef[3][cI]*(Q11_100)+coef[4][cI]*(Q11_010)+coef[5][cI]*(Q11_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 3; i[2]<(lth[2]-3); i[2]++){
        for(i[1] = 3; i[1]<(lth[1]-3); i[1]++){
            for(i[0] = 3; i[0]<(lth[0]-3); i[0]++){
                int I = ind(i,lth);
                
                double d006 = DERIV002(d004,i,lth,dx);
                double d060 = DERIV020(d040,i,lth,dx);
                double d600 = DERIV200(d400,i,lth,dx);
                
                double d005 = DERIV001(d004,i,lth,dx);
                double d050 = DERIV010(d040,i,lth,dx);
                double d500 = DERIV100(d400,i,lth,dx);
                double d015 = DERIV011(d004,i,lth,dx);
                double d033 = DERIV011(d022,i,lth,dx);
                double d051 = DERIV011(d040,i,lth,dx);
                double d105 = DERIV101(d004,i,lth,dx);
                double d150 = DERIV110(d040,i,lth,dx);
                double d303 = DERIV101(d202,i,lth,dx);
                double d330 = DERIV110(d220,i,lth,dx);
                double d501 = DERIV101(d400,i,lth,dx);
                double d510 = DERIV110(d400,i,lth,dx);
                
                double Q11_004 = DERIV002(Q11_002,i,lth,dx);
                double Q11_040 = DERIV020(Q11_020,i,lth,dx);
                double Q11_400 = DERIV200(Q11_200,i,lth,dx);
                double Q21_002 = DERIV002(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_020 = DERIV020(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_200 = DERIV200(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q12_002 = DERIV002(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_020 = DERIV020(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_200 = DERIV200(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q11_003 = DERIV001(Q11_002,i,lth,dx);
                double Q11_030 = DERIV010(Q11_020,i,lth,dx);
                double Q11_300 = DERIV100(Q11_200,i,lth,dx);
                double Q11_013 = DERIV011(Q11_002,i,lth,dx);
                double Q11_031 = DERIV011(Q11_020,i,lth,dx);
                double Q11_103 = DERIV101(Q11_002,i,lth,dx);
                double Q11_130 = DERIV110(Q11_020,i,lth,dx);
                double Q11_301 = DERIV101(Q11_200,i,lth,dx);
                double Q11_310 = DERIV110(Q11_200,i,lth,dx);
                
                double Q21_001 = DERIV001(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_010 = DERIV010(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_100 = DERIV100(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_011 = DERIV011(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_101 = DERIV101(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_110 = DERIV110(Q[ind2(2,1,NU,K)],i,lth,dx);
                
                double Q12_001 = DERIV001(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_010 = DERIV010(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_100 = DERIV100(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_011 = DERIV011(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_101 = DERIV101(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_110 = DERIV110(Q[ind2(1,2,NU,K)],i,lth,dx);

                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }
                
                Q[ind2(1,3,NU,K)][I] = Q[ind2(1,2,NU,K)][I]+coef[0][cI]*((1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d600)+coef[1][cI]*((1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d060)+coef[2][cI]*((1.0/90.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d006)+coef[6][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d150+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*d330+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d510)+coef[7][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d105+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[2]*dx[2])*d303+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d501)+coef[8][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d015+(-1.0/6.0)*(-1.0/6.0)*(dx[1]*dx[1])*(dx[2]*dx[2])*d033+(1.0/30.0)*1.0*(dx[1]*dx[1]*dx[1]*dx[1])*(1.0)*d051)+coef[3][cI]*((1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d500)+coef[4][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d050)+coef[5][cI]*((1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d005);
                Q[ind2(2,2,NU,K)][I] = coef[0][cI]*(Q12_200+((-1.0/12.0)*(dx[0]*dx[0])*Q11_400))+coef[1][cI]*(Q12_020+((-1.0/12.0)*(dx[1]*dx[1])*Q11_040))+coef[2][cI]*(Q12_002+((-1.0/12.0)*(dx[2]*dx[2])*Q11_004))+coef[6][cI]*(Q12_110+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q11_130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_310))+coef[7][cI]*(Q12_101+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_301))+coef[8][cI]*(Q12_011+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*Q11_031))+coef[3][cI]*(Q12_100+((-1.0/6.0)*(dx[0]*dx[0])*Q11_300))+coef[4][cI]*(Q12_010+((-1.0/6.0)*(dx[1]*dx[1])*Q11_030))+coef[5][cI]*(Q12_001+((-1.0/6.0)*(dx[2]*dx[2])*Q11_003));
                Q[ind2(3,1,NU,K)][I] = coef[0][cI]*(Q21_200)+coef[1][cI]*(Q21_020)+coef[2][cI]*(Q21_002)+coef[6][cI]*(Q21_110)+coef[7][cI]*(Q21_101)+coef[8][cI]*(Q21_011)+coef[3][cI]*(Q21_100)+coef[4][cI]*(Q21_010)+coef[5][cI]*(Q21_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,3,NU,K)][I];
    Z[ind3(2,0,0,NU,MU,MU)] = Q[ind2(2,2,NU,K)][I];
    Z[ind3(3,0,0,NU,MU,MU)] = Q[ind2(3,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_6_2.C"
}

void LagrangeDeriv_3D_8_0(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 4;
    int wth[3] = {(2*p),(2*p),(2*p)}; wth[axis] = p;
    int lth[3]= {(2*wth[0]+1),(2*wth[1]+1),(2*wth[2]+1)};
    int L_center = LagrangeData_center;
//    int L_wth    = 2*L_center + 1;
    int v0 = 1;
    int v1 = 2;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
//    double Q[NU*K][Lth];
    // Note that: must change to dynamic allocation because you run out of static memory
//    double d002[Lth];
//    double d020[Lth];
//    double d200[Lth];
//    double d004[Lth];
//    double d022[Lth];
//    double d040[Lth];
//    double d202[Lth];
//    double d220[Lth];
//    double d400[Lth];
//    double Q11_002[Lth];
//    double Q11_020[Lth];
//    double Q11_200[Lth];
//    double d006[Lth];
//    double d024[Lth];
//    double d042[Lth];
//    double d060[Lth];
//    double d204[Lth];
//    double d240[Lth];
//    double d402[Lth];
//    double d420[Lth];
//    double d600[Lth];
//    double Q11_004[Lth];
//    double Q11_022[Lth];
//    double Q11_040[Lth];
//    double Q11_202[Lth];
//    double Q11_220[Lth];
//    double Q11_400[Lth];
//    double Q21_002[Lth];
//    double Q21_020[Lth];
//    double Q21_200[Lth];
//    double Q12_002[Lth];
//    double Q12_020[Lth];
//    double Q12_200[Lth];
    
    double *Q[NU*K];
    int nuk;
    for(nuk = 0; nuk<(NU*K); nuk++)
        Q[nuk] = new(memory[nuk]) double[Lth];
    
    double * d002 = new(memory[nuk+0]) double[Lth];
    double * d020 = new(memory[nuk+1]) double[Lth];
    double * d200 = new(memory[nuk+2]) double[Lth];
    double * d004 = new(memory[nuk+3]) double[Lth];
    double * d022 = new(memory[nuk+4]) double[Lth];
    double * d040 = new(memory[nuk+5]) double[Lth];
    double * d202 = new(memory[nuk+6]) double[Lth];
    double * d220 = new(memory[nuk+7]) double[Lth];
    double * d400 = new(memory[nuk+8]) double[Lth];
    double * Q11_002 = new(memory[nuk+9]) double[Lth];
    double * Q11_020 = new(memory[nuk+10]) double[Lth];
    double * Q11_200 = new(memory[nuk+11]) double[Lth];
    double * d006 = new(memory[nuk+12]) double[Lth];
    double * d024 = new(memory[nuk+13]) double[Lth];
    double * d042 = new(memory[nuk+14]) double[Lth];
    double * d060 = new(memory[nuk+15]) double[Lth];
    double * d204 = new(memory[nuk+16]) double[Lth];
    double * d240 = new(memory[nuk+17]) double[Lth];
    double * d402 = new(memory[nuk+18]) double[Lth];
    double * d420 = new(memory[nuk+19]) double[Lth];
    double * d600 = new(memory[nuk+20]) double[Lth];
    double * Q11_004 = new(memory[nuk+21]) double[Lth];
    double * Q11_022 = new(memory[nuk+22]) double[Lth];
    double * Q11_040 = new(memory[nuk+23]) double[Lth];
    double * Q11_202 = new(memory[nuk+24]) double[Lth];
    double * Q11_220 = new(memory[nuk+25]) double[Lth];
    double * Q11_400 = new(memory[nuk+26]) double[Lth];
    double * Q21_002 = new(memory[nuk+27]) double[Lth];
    double * Q21_020 = new(memory[nuk+28]) double[Lth];
    double * Q21_200 = new(memory[nuk+29]) double[Lth];
    double * Q12_002 = new(memory[nuk+30]) double[Lth];
    double * Q12_020 = new(memory[nuk+31]) double[Lth];
    double * Q12_200 = new(memory[nuk+32]) double[Lth];

    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);

                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 2; i[2]<(lth[2]-2); i[2]++){
        for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
            for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
                int I = ind(i,lth);
                
                d004[I] = DERIV002(d002,i,lth,dx);
                d022[I] = DERIV002(d020,i,lth,dx);
                d040[I] = DERIV020(d020,i,lth,dx);
                d202[I] = DERIV200(d002,i,lth,dx);
                d220[I] = DERIV200(d020,i,lth,dx);
                d400[I] = DERIV200(d200,i,lth,dx);

                double d003 = DERIV001(d002,i,lth,dx);
                double d030 = DERIV010(d020,i,lth,dx);
                double d300 = DERIV100(d200,i,lth,dx);
                double d013 = DERIV011(d002,i,lth,dx);
                double d031 = DERIV011(d020,i,lth,dx);
                double d103 = DERIV101(d002,i,lth,dx);
                double d130 = DERIV110(d020,i,lth,dx);
                double d301 = DERIV101(d200,i,lth,dx);
                double d310 = DERIV110(d200,i,lth,dx);

                Q11_002[I] = DERIV002(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_020[I] = DERIV020(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_200[I] = DERIV200(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_001 = DERIV001(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_010 = DERIV010(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_100 = DERIV100(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_011 = DERIV011(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_101 = DERIV101(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_110 = DERIV110(Q[ind2(1,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d400[I])+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d040[I])+coef[2][cI]*((-1.0/12.0)*(dx[2]*dx[2])*d004[I])+coef[6][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d130+(-1.0/6.0)*(dx[0]*dx[0])*d310)+coef[7][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d103+(-1.0/6.0)*(dx[0]*dx[0])*d301)+coef[8][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d013+(-1.0/6.0)*(dx[1]*dx[1])*d031)+coef[3][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d300)+coef[4][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d030)+coef[5][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d003);
                Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_200[I])+coef[1][cI]*(Q11_020[I])+coef[2][cI]*(Q11_002[I])+coef[6][cI]*(Q11_110)+coef[7][cI]*(Q11_101)+coef[8][cI]*(Q11_011)+coef[3][cI]*(Q11_100)+coef[4][cI]*(Q11_010)+coef[5][cI]*(Q11_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 3; i[2]<(lth[2]-3); i[2]++){
        for(i[1] = 3; i[1]<(lth[1]-3); i[1]++){
            for(i[0] = 3; i[0]<(lth[0]-3); i[0]++){
                int I = ind(i,lth);
                
                d006[I] = DERIV002(d004,i,lth,dx);
                d024[I] = DERIV002(d022,i,lth,dx);
                d042[I] = DERIV002(d040,i,lth,dx);
                d060[I] = DERIV020(d040,i,lth,dx);
                d204[I] = DERIV002(d202,i,lth,dx);
                d240[I] = DERIV020(d220,i,lth,dx);
                d402[I] = DERIV200(d202,i,lth,dx);
                d420[I] = DERIV200(d220,i,lth,dx);
                d600[I] = DERIV200(d400,i,lth,dx);

                double d005 = DERIV001(d004,i,lth,dx);
                double d050 = DERIV010(d040,i,lth,dx);
                double d500 = DERIV100(d400,i,lth,dx);
                double d015 = DERIV011(d004,i,lth,dx);
                double d033 = DERIV011(d022,i,lth,dx);
                double d051 = DERIV011(d040,i,lth,dx);
                double d105 = DERIV101(d004,i,lth,dx);
                double d150 = DERIV110(d040,i,lth,dx);
                double d303 = DERIV101(d202,i,lth,dx);
                double d330 = DERIV110(d220,i,lth,dx);
                double d501 = DERIV101(d400,i,lth,dx);
                double d510 = DERIV110(d400,i,lth,dx);

                Q11_004[I] = DERIV002(Q11_002,i,lth,dx);
                Q11_022[I] = DERIV002(Q11_020,i,lth,dx);
                Q11_040[I] = DERIV020(Q11_020,i,lth,dx);
                Q11_202[I] = DERIV200(Q11_002,i,lth,dx);
                Q11_220[I] = DERIV200(Q11_020,i,lth,dx);
                Q11_400[I] = DERIV200(Q11_200,i,lth,dx);
                Q21_002[I] = DERIV002(Q[ind2(2,1,NU,K)],i,lth,dx);
                Q21_020[I] = DERIV020(Q[ind2(2,1,NU,K)],i,lth,dx);
                Q21_200[I] = DERIV200(Q[ind2(2,1,NU,K)],i,lth,dx);
                Q12_002[I] = DERIV002(Q[ind2(1,2,NU,K)],i,lth,dx);
                Q12_020[I] = DERIV020(Q[ind2(1,2,NU,K)],i,lth,dx);
                Q12_200[I] = DERIV200(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q11_003 = DERIV001(Q11_002,i,lth,dx);
                double Q11_030 = DERIV010(Q11_020,i,lth,dx);
                double Q11_300 = DERIV100(Q11_200,i,lth,dx);
                double Q11_013 = DERIV011(Q11_002,i,lth,dx);
                double Q11_031 = DERIV011(Q11_020,i,lth,dx);
                double Q11_103 = DERIV101(Q11_002,i,lth,dx);
                double Q11_130 = DERIV110(Q11_020,i,lth,dx);
                double Q11_301 = DERIV101(Q11_200,i,lth,dx);
                double Q11_310 = DERIV110(Q11_200,i,lth,dx);

                double Q21_001 = DERIV001(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_010 = DERIV010(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_100 = DERIV100(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_011 = DERIV011(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_101 = DERIV101(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_110 = DERIV110(Q[ind2(2,1,NU,K)],i,lth,dx);

                double Q12_001 = DERIV001(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_010 = DERIV010(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_100 = DERIV100(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_011 = DERIV011(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_101 = DERIV101(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_110 = DERIV110(Q[ind2(1,2,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,3,NU,K)][I] = Q[ind2(1,2,NU,K)][I]+coef[0][cI]*((1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d600[I])+coef[1][cI]*((1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d060[I])+coef[2][cI]*((1.0/90.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d006[I])+coef[6][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d150+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*d330+(1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d510)+coef[7][cI]*((1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d105+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[2]*dx[2])*d303+(1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d501)+coef[8][cI]*((1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d015+(-1.0/6.0)*(-1.0/6.0)*(dx[1]*dx[1])*(dx[2]*dx[2])*d033+(1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*(1.0)*d051)+coef[3][cI]*((1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d500)+coef[4][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d050)+coef[5][cI]*((1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d005);
                Q[ind2(2,2,NU,K)][I] = coef[0][cI]*(Q12_200[I]+((-1.0/12.0)*(dx[0]*dx[0])*Q11_400[I]))+coef[1][cI]*(Q12_020[I]+((-1.0/12.0)*(dx[1]*dx[1])*Q11_040[I]))+coef[2][cI]*(Q12_002[I]+((-1.0/12.0)*(dx[2]*dx[2])*Q11_004[I]))+coef[6][cI]*(Q12_110+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q11_130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_310))+coef[7][cI]*(Q12_101+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_301))+coef[8][cI]*(Q12_011+((-1.0/6.0)*(dx[2]*dx[2])*Q11_013+(-1.0/6.0)*(dx[1]*dx[1])*(1.0)*Q11_031))+coef[3][cI]*(Q12_100+((-1.0/6.0)*(dx[0]*dx[0])*Q11_300))+coef[4][cI]*(Q12_010+((-1.0/6.0)*(dx[1]*dx[1])*Q11_030))+coef[5][cI]*(Q12_001+((-1.0/6.0)*(dx[2]*dx[2])*Q11_003));
                Q[ind2(3,1,NU,K)][I] = coef[0][cI]*(Q21_200[I])+coef[1][cI]*(Q21_020[I])+coef[2][cI]*(Q21_002[I])+coef[6][cI]*(Q21_110)+coef[7][cI]*(Q21_101)+coef[8][cI]*(Q21_011)+coef[3][cI]*(Q21_100)+coef[4][cI]*(Q21_010)+coef[5][cI]*(Q21_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 4; i[2]<(lth[2]-4); i[2]++){
        for(i[1] = 4; i[1]<(lth[1]-4); i[1]++){
            for(i[0] = 4; i[0]<(lth[0]-4); i[0]++){
                int I = ind(i,lth);
                
                double d008 = DERIV002(d006,i,lth,dx);
                double d080 = DERIV020(d060,i,lth,dx);
                double d800 = DERIV200(d600,i,lth,dx);

                double d007 = DERIV001(d006,i,lth,dx);
                double d070 = DERIV010(d060,i,lth,dx);
                double d700 = DERIV100(d600,i,lth,dx);
                double d017 = DERIV011(d006,i,lth,dx);
                double d035 = DERIV011(d024,i,lth,dx);
                double d053 = DERIV011(d042,i,lth,dx);
                double d071 = DERIV011(d060,i,lth,dx);
                double d107 = DERIV101(d006,i,lth,dx);
                double d170 = DERIV110(d060,i,lth,dx);
                double d305 = DERIV101(d204,i,lth,dx);
                double d350 = DERIV110(d240,i,lth,dx);
                double d503 = DERIV101(d402,i,lth,dx);
                double d530 = DERIV110(d420,i,lth,dx);
                double d701 = DERIV101(d600,i,lth,dx);
                double d710 = DERIV110(d600,i,lth,dx);

                double Q11_006 = DERIV002(Q11_004,i,lth,dx);
                double Q11_060 = DERIV020(Q11_040,i,lth,dx);
                double Q11_600 = DERIV200(Q11_400,i,lth,dx);
                double Q21_004 = DERIV002(Q21_002,i,lth,dx);
                double Q21_040 = DERIV020(Q21_020,i,lth,dx);
                double Q21_400 = DERIV200(Q21_200,i,lth,dx);
                double Q12_004 = DERIV002(Q12_002,i,lth,dx);
                double Q12_040 = DERIV020(Q12_020,i,lth,dx);
                double Q12_400 = DERIV200(Q12_200,i,lth,dx);
                double Q22_002 = DERIV002(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_020 = DERIV020(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_200 = DERIV200(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q31_002 = DERIV002(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_020 = DERIV020(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_200 = DERIV200(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q13_002 = DERIV002(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_020 = DERIV020(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_200 = DERIV200(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q11_005 = DERIV001(Q11_004,i,lth,dx);
                double Q11_050 = DERIV010(Q11_040,i,lth,dx);
                double Q11_500 = DERIV100(Q11_400,i,lth,dx);
                double Q11_015 = DERIV011(Q11_004,i,lth,dx);
                double Q11_033 = DERIV011(Q11_022,i,lth,dx);
                double Q11_051 = DERIV011(Q11_040,i,lth,dx);
                double Q11_105 = DERIV101(Q11_004,i,lth,dx);
                double Q11_150 = DERIV110(Q11_040,i,lth,dx);
                double Q11_303 = DERIV101(Q11_202,i,lth,dx);
                double Q11_330 = DERIV110(Q11_220,i,lth,dx);
                double Q11_501 = DERIV101(Q11_400,i,lth,dx);
                double Q11_510 = DERIV110(Q11_400,i,lth,dx);

                double Q21_003 = DERIV001(Q21_002,i,lth,dx);
                double Q21_030 = DERIV010(Q21_020,i,lth,dx);
                double Q21_300 = DERIV100(Q21_200,i,lth,dx);
                double Q21_013 = DERIV011(Q21_002,i,lth,dx);
                double Q21_031 = DERIV011(Q21_020,i,lth,dx);
                double Q21_103 = DERIV101(Q21_002,i,lth,dx);
                double Q21_130 = DERIV110(Q21_020,i,lth,dx);
                double Q21_301 = DERIV101(Q21_200,i,lth,dx);
                double Q21_310 = DERIV110(Q21_200,i,lth,dx);

                double Q12_003 = DERIV001(Q12_002,i,lth,dx);
                double Q12_030 = DERIV010(Q12_020,i,lth,dx);
                double Q12_300 = DERIV100(Q12_200,i,lth,dx);
                double Q12_013 = DERIV011(Q12_002,i,lth,dx);
                double Q12_031 = DERIV011(Q12_020,i,lth,dx);
                double Q12_103 = DERIV101(Q12_002,i,lth,dx);
                double Q12_130 = DERIV110(Q12_020,i,lth,dx);
                double Q12_301 = DERIV101(Q12_200,i,lth,dx);
                double Q12_310 = DERIV110(Q12_200,i,lth,dx);

                double Q22_001 = DERIV001(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_010 = DERIV010(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_100 = DERIV100(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_011 = DERIV011(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_101 = DERIV101(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_110 = DERIV110(Q[ind2(2,2,NU,K)],i,lth,dx);

                double Q31_001 = DERIV001(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_010 = DERIV010(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_100 = DERIV100(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_011 = DERIV011(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_101 = DERIV101(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_110 = DERIV110(Q[ind2(3,1,NU,K)],i,lth,dx);

                double Q13_001 = DERIV001(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_010 = DERIV010(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_100 = DERIV100(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_011 = DERIV011(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_101 = DERIV101(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_110 = DERIV110(Q[ind2(1,3,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,4,NU,K)][I] = Q[ind2(1,3,NU,K)][I]+coef[0][cI]*((-1.0/560.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d800)+coef[1][cI]*((-1.0/560.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d080)+coef[2][cI]*((-1.0/560.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d008)+coef[6][cI]*((-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d170+(-1.0/6.0)*(1.0/30.0)*(dx[0]*dx[0])*(dx[1]*dx[1]*dx[1]*dx[1])*d350+(1.0/30.0)*(-1.0/6.0)*(dx[0]*dx[0]*dx[0]*dx[0])*(dx[1]*dx[1])*d530+(-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d710)+coef[7][cI]*((-1.0/140.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d107+(-1.0/6.0)*(1.0/30.0)*(dx[0]*dx[0])*(dx[2]*dx[2]*dx[2]*dx[2])*d305+(1.0/30.0)*(-1.0/6.0)*(dx[0]*dx[0]*dx[0]*dx[0])*(dx[2]*dx[2])*d503+(-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d701)+coef[8][cI]*(1.0*(-1.0/140.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d017+(-1.0/6.0)*(1.0/30.0)*(dx[1]*dx[1])*(dx[2]*dx[2]*dx[2]*dx[2])*d035+(1.0/30.0)*(-1.0/6.0)*(dx[1]*dx[1]*dx[1]*dx[1])*(dx[2]*dx[2])*d053+(-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d071)+coef[3][cI]*((-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d700)+coef[4][cI]*((-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d070)+coef[5][cI]*((-1.0/140.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d007);
                Q[ind2(2,3,NU,K)][I]  = coef[0][cI]*(Q13_200+((-1.0/12.0)*(dx[0]*dx[0])*Q12_400+(1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_600))+coef[1][cI]*(Q13_020+((-1.0/12.0)*(dx[1]*dx[1])*Q12_040+(1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_060))+coef[2][cI]*(Q13_002+((-1.0/12.0)*(dx[2]*dx[2])*Q12_004+(1.0/90.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_006))+coef[6][cI]*(Q13_110+((-1.0/6.0)*(dx[1]*dx[1])*Q12_130+(-1.0/6.0)*(dx[0]*dx[0])*Q12_310+(1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_150+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*Q11_330+(1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_510))+coef[7][cI]*(Q13_101+((-1.0/6.0)*(dx[2]*dx[2])*Q12_103+(-1.0/6.0)*(dx[0]*dx[0])*Q12_301+(1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_105+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[2]*dx[2])*Q11_303+(1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_501))+coef[8][cI]*(Q13_011+((-1.0/6.0)*(dx[2]*dx[2])*Q12_013+(-1.0/6.0)*(dx[1]*dx[1])*Q12_031+(1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_015+(-1.0/6.0)*(-1.0/6.0)*(dx[1]*dx[1])*(dx[2]*dx[2])*Q11_033+(1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_051))+coef[3][cI]*(Q13_100+((-1.0/6.0)*(dx[0]*dx[0])*Q12_300+(1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_500))+coef[4][cI]*(Q13_010+((-1.0/6.0)*(dx[1]*dx[1])*Q12_030+(1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_050))+coef[5][cI]*(Q13_001+((-1.0/6.0)*(dx[2]*dx[2])*Q12_003+(1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_005));
                Q[ind2(3,2,NU,K)][I]  = coef[0][cI]*(Q22_200+((-1.0/12.0)*(dx[0]*dx[0])*Q21_400))+coef[1][cI]*(Q22_020+((-1.0/12.0)*(dx[1]*dx[1])*Q21_040))+coef[2][cI]*(Q22_002+((-1.0/12.0)*(dx[2]*dx[2])*Q21_004))+coef[6][cI]*(Q22_110+((-1.0/6.0)*(dx[1]*dx[1])*Q21_130+(-1.0/6.0)*(dx[0]*dx[0])*(1.0)*Q21_310))+coef[7][cI]*(Q22_101+((-1.0/6.0)*(dx[2]*dx[2])*Q21_103+(-1.0/6.0)*(dx[0]*dx[0])*Q21_301))+coef[8][cI]*(Q22_011+((-1.0/6.0)*(dx[2]*dx[2])*Q21_013+(-1.0/6.0)*(dx[1]*dx[1])*(1.0)*Q21_031))+coef[3][cI]*(Q22_100+((-1.0/6.0)*(dx[0]*dx[0])*Q21_300))+coef[4][cI]*(Q22_010+((-1.0/6.0)*(dx[1]*dx[1])*Q21_030))+coef[5][cI]*(Q22_001+((-1.0/6.0)*(dx[2]*dx[2])*Q21_003));
                Q[ind2(4,1,NU,K)][I]  = coef[0][cI]*(Q31_200)+coef[1][cI]*(Q31_020)+coef[2][cI]*(Q31_002)+coef[6][cI]*(Q31_110)+coef[7][cI]*(Q31_101)+coef[8][cI]*(Q31_011)+coef[3][cI]*(Q31_100)+coef[4][cI]*(Q31_010)+coef[5][cI]*(Q31_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,4,NU,K)][I];
    Z[ind3(2,0,0,NU,MU,MU)] = Q[ind2(2,3,NU,K)][I];
    Z[ind3(3,0,0,NU,MU,MU)] = Q[ind2(3,2,NU,K)][I];
    Z[ind3(4,0,0,NU,MU,MU)] = Q[ind2(4,1,NU,K)][I];
    
//    printf("%f,%f,%f,%f,%f\n",Z[ind3(0,0,0,NU,MU,MU)], Z[ind3(1,0,0,NU,MU,MU)], Z[ind3(2,0,0,NU,MU,MU)], Z[ind3(3,0,0,NU,MU,MU)], Z[ind3(4,0,0,NU,MU,MU)]);
    
    /* Generate the tangential derivatives here */
    
#include "tangentialDerivs_3D_8_0.C"
}

void LagrangeDeriv_3D_8_1(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 4;
    int wth[3] = {(2*p),(2*p),(2*p)}; wth[axis] = p;
    int lth[3]= {(2*wth[0]+1),(2*wth[1]+1),(2*wth[2]+1)};
    int L_center = LagrangeData_center;
//    int L_wth    = 2*L_center + 1;
    int v0 = 0;
    int v1 = 2;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
//    double Q[NU*K][Lth];
//    double d002[Lth];
//    double d020[Lth];
//    double d200[Lth];
//    double d004[Lth];
//    double d022[Lth];
//    double d040[Lth];
//    double d202[Lth];
//    double d220[Lth];
//    double d400[Lth];
//    double Q11_002[Lth];
//    double Q11_020[Lth];
//    double Q11_200[Lth];
//    double d006[Lth];
//    double d024[Lth];
//    double d042[Lth];
//    double d060[Lth];
//    double d204[Lth];
//    double d240[Lth];
//    double d402[Lth];
//    double d420[Lth];
//    double d600[Lth];
//    double Q11_004[Lth];
//    double Q11_022[Lth];
//    double Q11_040[Lth];
//    double Q11_202[Lth];
//    double Q11_220[Lth];
//    double Q11_400[Lth];
//    double Q21_002[Lth];
//    double Q21_020[Lth];
//    double Q21_200[Lth];
//    double Q12_002[Lth];
//    double Q12_020[Lth];
//    double Q12_200[Lth];
    
    double *Q[NU*K];
    int nuk;
    for(nuk = 0; nuk<(NU*K); nuk++)
        Q[nuk] = new(memory[nuk]) double[Lth];
    
    double * d002 = new(memory[nuk+0]) double[Lth];
    double * d020 = new(memory[nuk+1]) double[Lth];
    double * d200 = new(memory[nuk+2]) double[Lth];
    double * d004 = new(memory[nuk+3]) double[Lth];
    double * d022 = new(memory[nuk+4]) double[Lth];
    double * d040 = new(memory[nuk+5]) double[Lth];
    double * d202 = new(memory[nuk+6]) double[Lth];
    double * d220 = new(memory[nuk+7]) double[Lth];
    double * d400 = new(memory[nuk+8]) double[Lth];
    double * Q11_002 = new(memory[nuk+9]) double[Lth];
    double * Q11_020 = new(memory[nuk+10]) double[Lth];
    double * Q11_200 = new(memory[nuk+11]) double[Lth];
    double * d006 = new(memory[nuk+12]) double[Lth];
    double * d024 = new(memory[nuk+13]) double[Lth];
    double * d042 = new(memory[nuk+14]) double[Lth];
    double * d060 = new(memory[nuk+15]) double[Lth];
    double * d204 = new(memory[nuk+16]) double[Lth];
    double * d240 = new(memory[nuk+17]) double[Lth];
    double * d402 = new(memory[nuk+18]) double[Lth];
    double * d420 = new(memory[nuk+19]) double[Lth];
    double * d600 = new(memory[nuk+20]) double[Lth];
    double * Q11_004 = new(memory[nuk+21]) double[Lth];
    double * Q11_022 = new(memory[nuk+22]) double[Lth];
    double * Q11_040 = new(memory[nuk+23]) double[Lth];
    double * Q11_202 = new(memory[nuk+24]) double[Lth];
    double * Q11_220 = new(memory[nuk+25]) double[Lth];
    double * Q11_400 = new(memory[nuk+26]) double[Lth];
    double * Q21_002 = new(memory[nuk+27]) double[Lth];
    double * Q21_020 = new(memory[nuk+28]) double[Lth];
    double * Q21_200 = new(memory[nuk+29]) double[Lth];
    double * Q12_002 = new(memory[nuk+30]) double[Lth];
    double * Q12_020 = new(memory[nuk+31]) double[Lth];
    double * Q12_200 = new(memory[nuk+32]) double[Lth];


    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);

                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 2; i[2]<(lth[2]-2); i[2]++){
        for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
            for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
                int I = ind(i,lth);
                
                d004[I] = DERIV002(d002,i,lth,dx);
                d022[I] = DERIV002(d020,i,lth,dx);
                d040[I] = DERIV020(d020,i,lth,dx);
                d202[I] = DERIV200(d002,i,lth,dx);
                d220[I] = DERIV200(d020,i,lth,dx);
                d400[I] = DERIV200(d200,i,lth,dx);

                double d003 = DERIV001(d002,i,lth,dx);
                double d030 = DERIV010(d020,i,lth,dx);
                double d300 = DERIV100(d200,i,lth,dx);
                double d013 = DERIV011(d002,i,lth,dx);
                double d031 = DERIV011(d020,i,lth,dx);
                double d103 = DERIV101(d002,i,lth,dx);
                double d130 = DERIV110(d020,i,lth,dx);
                double d301 = DERIV101(d200,i,lth,dx);
                double d310 = DERIV110(d200,i,lth,dx);

                Q11_002[I] = DERIV002(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_020[I] = DERIV020(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_200[I] = DERIV200(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_001 = DERIV001(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_010 = DERIV010(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_100 = DERIV100(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_011 = DERIV011(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_101 = DERIV101(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_110 = DERIV110(Q[ind2(1,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d400[I])+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d040[I])+coef[2][cI]*((-1.0/12.0)*(dx[2]*dx[2])*d004[I])+coef[6][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d310)+coef[7][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d301)+coef[8][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*d031)+coef[3][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d300)+coef[4][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d030)+coef[5][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d003);
                Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_200[I])+coef[1][cI]*(Q11_020[I])+coef[2][cI]*(Q11_002[I])+coef[6][cI]*(Q11_110)+coef[7][cI]*(Q11_101)+coef[8][cI]*(Q11_011)+coef[3][cI]*(Q11_100)+coef[4][cI]*(Q11_010)+coef[5][cI]*(Q11_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 3; i[2]<(lth[2]-3); i[2]++){
        for(i[1] = 3; i[1]<(lth[1]-3); i[1]++){
            for(i[0] = 3; i[0]<(lth[0]-3); i[0]++){
                int I = ind(i,lth);
                
                d006[I] = DERIV002(d004,i,lth,dx);
                d024[I] = DERIV002(d022,i,lth,dx);
                d042[I] = DERIV002(d040,i,lth,dx);
                d060[I] = DERIV020(d040,i,lth,dx);
                d204[I] = DERIV002(d202,i,lth,dx);
                d240[I] = DERIV020(d220,i,lth,dx);
                d402[I] = DERIV200(d202,i,lth,dx);
                d420[I] = DERIV200(d220,i,lth,dx);
                d600[I] = DERIV200(d400,i,lth,dx);

                double d005 = DERIV001(d004,i,lth,dx);
                double d050 = DERIV010(d040,i,lth,dx);
                double d500 = DERIV100(d400,i,lth,dx);
                double d015 = DERIV011(d004,i,lth,dx);
                double d033 = DERIV011(d022,i,lth,dx);
                double d051 = DERIV011(d040,i,lth,dx);
                double d105 = DERIV101(d004,i,lth,dx);
                double d150 = DERIV110(d040,i,lth,dx);
                double d303 = DERIV101(d202,i,lth,dx);
                double d330 = DERIV110(d220,i,lth,dx);
                double d501 = DERIV101(d400,i,lth,dx);
                double d510 = DERIV110(d400,i,lth,dx);

                Q11_004[I] = DERIV002(Q11_002,i,lth,dx);
                Q11_022[I] = DERIV002(Q11_020,i,lth,dx);
                Q11_040[I] = DERIV020(Q11_020,i,lth,dx);
                Q11_202[I] = DERIV200(Q11_002,i,lth,dx);
                Q11_220[I] = DERIV200(Q11_020,i,lth,dx);
                Q11_400[I] = DERIV200(Q11_200,i,lth,dx);
                Q21_002[I] = DERIV002(Q[ind2(2,1,NU,K)],i,lth,dx);
                Q21_020[I] = DERIV020(Q[ind2(2,1,NU,K)],i,lth,dx);
                Q21_200[I] = DERIV200(Q[ind2(2,1,NU,K)],i,lth,dx);
                Q12_002[I] = DERIV002(Q[ind2(1,2,NU,K)],i,lth,dx);
                Q12_020[I] = DERIV020(Q[ind2(1,2,NU,K)],i,lth,dx);
                Q12_200[I] = DERIV200(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q11_003 = DERIV001(Q11_002,i,lth,dx);
                double Q11_030 = DERIV010(Q11_020,i,lth,dx);
                double Q11_300 = DERIV100(Q11_200,i,lth,dx);
                double Q11_013 = DERIV011(Q11_002,i,lth,dx);
                double Q11_031 = DERIV011(Q11_020,i,lth,dx);
                double Q11_103 = DERIV101(Q11_002,i,lth,dx);
                double Q11_130 = DERIV110(Q11_020,i,lth,dx);
                double Q11_301 = DERIV101(Q11_200,i,lth,dx);
                double Q11_310 = DERIV110(Q11_200,i,lth,dx);

                double Q21_001 = DERIV001(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_010 = DERIV010(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_100 = DERIV100(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_011 = DERIV011(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_101 = DERIV101(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_110 = DERIV110(Q[ind2(2,1,NU,K)],i,lth,dx);

                double Q12_001 = DERIV001(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_010 = DERIV010(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_100 = DERIV100(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_011 = DERIV011(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_101 = DERIV101(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_110 = DERIV110(Q[ind2(1,2,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,3,NU,K)][I] = Q[ind2(1,2,NU,K)][I]+coef[0][cI]*((1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d600[I])+coef[1][cI]*((1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d060[I])+coef[2][cI]*((1.0/90.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d006[I])+coef[6][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d150+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*d330+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d510)+coef[7][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d105+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[2]*dx[2])*d303+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d501)+coef[8][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d015+(-1.0/6.0)*(-1.0/6.0)*(dx[1]*dx[1])*(dx[2]*dx[2])*d033+(1.0/30.0)*1.0*(dx[1]*dx[1]*dx[1]*dx[1])*(1.0)*d051)+coef[3][cI]*((1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d500)+coef[4][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d050)+coef[5][cI]*((1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d005);
                Q[ind2(2,2,NU,K)][I] = coef[0][cI]*(Q12_200[I]+((-1.0/12.0)*(dx[0]*dx[0])*Q11_400[I]))+coef[1][cI]*(Q12_020[I]+((-1.0/12.0)*(dx[1]*dx[1])*Q11_040[I]))+coef[2][cI]*(Q12_002[I]+((-1.0/12.0)*(dx[2]*dx[2])*Q11_004[I]))+coef[6][cI]*(Q12_110+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q11_130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_310))+coef[7][cI]*(Q12_101+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_301))+coef[8][cI]*(Q12_011+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*Q11_031))+coef[3][cI]*(Q12_100+((-1.0/6.0)*(dx[0]*dx[0])*Q11_300))+coef[4][cI]*(Q12_010+((-1.0/6.0)*(dx[1]*dx[1])*Q11_030))+coef[5][cI]*(Q12_001+((-1.0/6.0)*(dx[2]*dx[2])*Q11_003));
                Q[ind2(3,1,NU,K)][I] = coef[0][cI]*(Q21_200[I])+coef[1][cI]*(Q21_020[I])+coef[2][cI]*(Q21_002[I])+coef[6][cI]*(Q21_110)+coef[7][cI]*(Q21_101)+coef[8][cI]*(Q21_011)+coef[3][cI]*(Q21_100)+coef[4][cI]*(Q21_010)+coef[5][cI]*(Q21_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 4; i[2]<(lth[2]-4); i[2]++){
        for(i[1] = 4; i[1]<(lth[1]-4); i[1]++){
            for(i[0] = 4; i[0]<(lth[0]-4); i[0]++){
                int I = ind(i,lth);
                
                double d008 = DERIV002(d006,i,lth,dx);
                double d080 = DERIV020(d060,i,lth,dx);
                double d800 = DERIV200(d600,i,lth,dx);

                double d007 = DERIV001(d006,i,lth,dx);
                double d070 = DERIV010(d060,i,lth,dx);
                double d700 = DERIV100(d600,i,lth,dx);
                double d017 = DERIV011(d006,i,lth,dx);
                double d035 = DERIV011(d024,i,lth,dx);
                double d053 = DERIV011(d042,i,lth,dx);
                double d071 = DERIV011(d060,i,lth,dx);
                double d107 = DERIV101(d006,i,lth,dx);
                double d170 = DERIV110(d060,i,lth,dx);
                double d305 = DERIV101(d204,i,lth,dx);
                double d350 = DERIV110(d240,i,lth,dx);
                double d503 = DERIV101(d402,i,lth,dx);
                double d530 = DERIV110(d420,i,lth,dx);
                double d701 = DERIV101(d600,i,lth,dx);
                double d710 = DERIV110(d600,i,lth,dx);

                double Q11_006 = DERIV002(Q11_004,i,lth,dx);
                double Q11_060 = DERIV020(Q11_040,i,lth,dx);
                double Q11_600 = DERIV200(Q11_400,i,lth,dx);
                double Q21_004 = DERIV002(Q21_002,i,lth,dx);
                double Q21_040 = DERIV020(Q21_020,i,lth,dx);
                double Q21_400 = DERIV200(Q21_200,i,lth,dx);
                double Q12_004 = DERIV002(Q12_002,i,lth,dx);
                double Q12_040 = DERIV020(Q12_020,i,lth,dx);
                double Q12_400 = DERIV200(Q12_200,i,lth,dx);
                double Q22_002 = DERIV002(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_020 = DERIV020(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_200 = DERIV200(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q31_002 = DERIV002(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_020 = DERIV020(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_200 = DERIV200(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q13_002 = DERIV002(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_020 = DERIV020(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_200 = DERIV200(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q11_005 = DERIV001(Q11_004,i,lth,dx);
                double Q11_050 = DERIV010(Q11_040,i,lth,dx);
                double Q11_500 = DERIV100(Q11_400,i,lth,dx);
                double Q11_015 = DERIV011(Q11_004,i,lth,dx);
                double Q11_033 = DERIV011(Q11_022,i,lth,dx);
                double Q11_051 = DERIV011(Q11_040,i,lth,dx);
                double Q11_105 = DERIV101(Q11_004,i,lth,dx);
                double Q11_150 = DERIV110(Q11_040,i,lth,dx);
                double Q11_303 = DERIV101(Q11_202,i,lth,dx);
                double Q11_330 = DERIV110(Q11_220,i,lth,dx);
                double Q11_501 = DERIV101(Q11_400,i,lth,dx);
                double Q11_510 = DERIV110(Q11_400,i,lth,dx);

                double Q21_003 = DERIV001(Q21_002,i,lth,dx);
                double Q21_030 = DERIV010(Q21_020,i,lth,dx);
                double Q21_300 = DERIV100(Q21_200,i,lth,dx);
                double Q21_013 = DERIV011(Q21_002,i,lth,dx);
                double Q21_031 = DERIV011(Q21_020,i,lth,dx);
                double Q21_103 = DERIV101(Q21_002,i,lth,dx);
                double Q21_130 = DERIV110(Q21_020,i,lth,dx);
                double Q21_301 = DERIV101(Q21_200,i,lth,dx);
                double Q21_310 = DERIV110(Q21_200,i,lth,dx);

                double Q12_003 = DERIV001(Q12_002,i,lth,dx);
                double Q12_030 = DERIV010(Q12_020,i,lth,dx);
                double Q12_300 = DERIV100(Q12_200,i,lth,dx);
                double Q12_013 = DERIV011(Q12_002,i,lth,dx);
                double Q12_031 = DERIV011(Q12_020,i,lth,dx);
                double Q12_103 = DERIV101(Q12_002,i,lth,dx);
                double Q12_130 = DERIV110(Q12_020,i,lth,dx);
                double Q12_301 = DERIV101(Q12_200,i,lth,dx);
                double Q12_310 = DERIV110(Q12_200,i,lth,dx);

                double Q22_001 = DERIV001(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_010 = DERIV010(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_100 = DERIV100(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_011 = DERIV011(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_101 = DERIV101(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_110 = DERIV110(Q[ind2(2,2,NU,K)],i,lth,dx);

                double Q31_001 = DERIV001(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_010 = DERIV010(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_100 = DERIV100(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_011 = DERIV011(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_101 = DERIV101(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_110 = DERIV110(Q[ind2(3,1,NU,K)],i,lth,dx);

                double Q13_001 = DERIV001(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_010 = DERIV010(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_100 = DERIV100(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_011 = DERIV011(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_101 = DERIV101(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_110 = DERIV110(Q[ind2(1,3,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,4,NU,K)][I] = Q[ind2(1,3,NU,K)][I]+coef[0][cI]*((-1.0/560.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d800)+coef[1][cI]*((-1.0/560.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d080)+coef[2][cI]*((-1.0/560.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d008)+coef[6][cI]*((-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d170+(-1.0/6.0)*(1.0/30.0)*(dx[0]*dx[0])*(dx[1]*dx[1]*dx[1]*dx[1])*d350+(1.0/30.0)*(-1.0/6.0)*(dx[0]*dx[0]*dx[0]*dx[0])*(dx[1]*dx[1])*d530+(-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d710)+coef[7][cI]*((-1.0/140.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d107+(-1.0/6.0)*(1.0/30.0)*(dx[0]*dx[0])*(dx[2]*dx[2]*dx[2]*dx[2])*d305+(1.0/30.0)*(-1.0/6.0)*(dx[0]*dx[0]*dx[0]*dx[0])*(dx[2]*dx[2])*d503+(-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d701)+coef[8][cI]*(1.0*(-1.0/140.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d017+(-1.0/6.0)*(1.0/30.0)*(dx[1]*dx[1])*(dx[2]*dx[2]*dx[2]*dx[2])*d035+(1.0/30.0)*(-1.0/6.0)*(dx[1]*dx[1]*dx[1]*dx[1])*(dx[2]*dx[2])*d053+(-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d071)+coef[3][cI]*((-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d700)+coef[4][cI]*((-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d070)+coef[5][cI]*((-1.0/140.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d007);
                Q[ind2(2,3,NU,K)][I]  = coef[0][cI]*(Q13_200+((-1.0/12.0)*(dx[0]*dx[0])*Q12_400+(1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_600))+coef[1][cI]*(Q13_020+((-1.0/12.0)*(dx[1]*dx[1])*Q12_040+(1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_060))+coef[2][cI]*(Q13_002+((-1.0/12.0)*(dx[2]*dx[2])*Q12_004+(1.0/90.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_006))+coef[6][cI]*(Q13_110+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q12_130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q12_310+1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_150+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*Q11_330+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*Q11_510))+coef[7][cI]*(Q13_101+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q12_103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q12_301+1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_105+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[2]*dx[2])*Q11_303+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*Q11_501))+coef[8][cI]*(Q13_011+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q12_013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*Q12_031+1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_015+(-1.0/6.0)*(-1.0/6.0)*(dx[1]*dx[1])*(dx[2]*dx[2])*Q11_033+(1.0/30.0)*1.0*(dx[1]*dx[1]*dx[1]*dx[1])*(1.0)*Q11_051))+coef[3][cI]*(Q13_100+((-1.0/6.0)*(dx[0]*dx[0])*Q12_300+(1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_500))+coef[4][cI]*(Q13_010+((-1.0/6.0)*(dx[1]*dx[1])*Q12_030+(1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_050))+coef[5][cI]*(Q13_001+((-1.0/6.0)*(dx[2]*dx[2])*Q12_003+(1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_005));
                Q[ind2(3,2,NU,K)][I]  = coef[0][cI]*(Q22_200+((-1.0/12.0)*(dx[0]*dx[0])*Q21_400))+coef[1][cI]*(Q22_020+((-1.0/12.0)*(dx[1]*dx[1])*Q21_040))+coef[2][cI]*(Q22_002+((-1.0/12.0)*(dx[2]*dx[2])*Q21_004))+coef[6][cI]*(Q22_110+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q21_130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q21_310))+coef[7][cI]*(Q22_101+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q21_103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q21_301))+coef[8][cI]*(Q22_011+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q21_013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*Q21_031))+coef[3][cI]*(Q22_100+((-1.0/6.0)*(dx[0]*dx[0])*Q21_300))+coef[4][cI]*(Q22_010+((-1.0/6.0)*(dx[1]*dx[1])*Q21_030))+coef[5][cI]*(Q22_001+((-1.0/6.0)*(dx[2]*dx[2])*Q21_003));
                Q[ind2(4,1,NU,K)][I]  = coef[0][cI]*(Q31_200)+coef[1][cI]*(Q31_020)+coef[2][cI]*(Q31_002)+coef[6][cI]*(Q31_110)+coef[7][cI]*(Q31_101)+coef[8][cI]*(Q31_011)+coef[3][cI]*(Q31_100)+coef[4][cI]*(Q31_010)+coef[5][cI]*(Q31_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,4,NU,K)][I];
    Z[ind3(2,0,0,NU,MU,MU)] = Q[ind2(2,3,NU,K)][I];
    Z[ind3(3,0,0,NU,MU,MU)] = Q[ind2(3,2,NU,K)][I];
    Z[ind3(4,0,0,NU,MU,MU)] = Q[ind2(4,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_8_1.C"
}

void LagrangeDeriv_3D_8_2(double Z[], int *Ind, int *LInd, double *LagrangeData, int LagrangeData_center, double **coef, int coefWth[3], int coefLth[3], bool cstCoef, double dx[3], int axis, int side, int NU, char **&memory){
    int p = 4;
    int wth[3] = {(2*p),(2*p),(2*p)}; wth[axis] = p;
    int lth[3]= {(2*wth[0]+1),(2*wth[1]+1),(2*wth[2]+1)};
    int L_center = LagrangeData_center;
//    int L_wth    = 2*L_center + 1;
    int v0 = 0;
    int v1 = 1;
    int cI = 0;
    int cInd[3] = {0,0,0};
    int K = p+1;
    
    int Lth =(lth[0]*lth[1]*lth[2]);
//    double Q[NU*K][Lth];
//    double d002[Lth];
//    double d020[Lth];
//    double d200[Lth];
//    double d004[Lth];
//    double d022[Lth];
//    double d040[Lth];
//    double d202[Lth];
//    double d220[Lth];
//    double d400[Lth];
//    double Q11_002[Lth];
//    double Q11_020[Lth];
//    double Q11_200[Lth];
//    double d006[Lth];
//    double d024[Lth];
//    double d042[Lth];
//    double d060[Lth];
//    double d204[Lth];
//    double d240[Lth];
//    double d402[Lth];
//    double d420[Lth];
//    double d600[Lth];
//    double Q11_004[Lth];
//    double Q11_022[Lth];
//    double Q11_040[Lth];
//    double Q11_202[Lth];
//    double Q11_220[Lth];
//    double Q11_400[Lth];
//    double Q21_002[Lth];
//    double Q21_020[Lth];
//    double Q21_200[Lth];
//    double Q12_002[Lth];
//    double Q12_020[Lth];
//    double Q12_200[Lth];
    
    double *Q[NU*K];
    int nuk;
    for(nuk = 0; nuk<(NU*K); nuk++)
        Q[nuk] = new(memory[nuk]) double[Lth];
    
    double * d002 = new(memory[nuk+0]) double[Lth];
    double * d020 = new(memory[nuk+1]) double[Lth];
    double * d200 = new(memory[nuk+2]) double[Lth];
    double * d004 = new(memory[nuk+3]) double[Lth];
    double * d022 = new(memory[nuk+4]) double[Lth];
    double * d040 = new(memory[nuk+5]) double[Lth];
    double * d202 = new(memory[nuk+6]) double[Lth];
    double * d220 = new(memory[nuk+7]) double[Lth];
    double * d400 = new(memory[nuk+8]) double[Lth];
    double * Q11_002 = new(memory[nuk+9]) double[Lth];
    double * Q11_020 = new(memory[nuk+10]) double[Lth];
    double * Q11_200 = new(memory[nuk+11]) double[Lth];
    double * d006 = new(memory[nuk+12]) double[Lth];
    double * d024 = new(memory[nuk+13]) double[Lth];
    double * d042 = new(memory[nuk+14]) double[Lth];
    double * d060 = new(memory[nuk+15]) double[Lth];
    double * d204 = new(memory[nuk+16]) double[Lth];
    double * d240 = new(memory[nuk+17]) double[Lth];
    double * d402 = new(memory[nuk+18]) double[Lth];
    double * d420 = new(memory[nuk+19]) double[Lth];
    double * d600 = new(memory[nuk+20]) double[Lth];
    double * Q11_004 = new(memory[nuk+21]) double[Lth];
    double * Q11_022 = new(memory[nuk+22]) double[Lth];
    double * Q11_040 = new(memory[nuk+23]) double[Lth];
    double * Q11_202 = new(memory[nuk+24]) double[Lth];
    double * Q11_220 = new(memory[nuk+25]) double[Lth];
    double * Q11_400 = new(memory[nuk+26]) double[Lth];
    double * Q21_002 = new(memory[nuk+27]) double[Lth];
    double * Q21_020 = new(memory[nuk+28]) double[Lth];
    double * Q21_200 = new(memory[nuk+29]) double[Lth];
    double * Q12_002 = new(memory[nuk+30]) double[Lth];
    double * Q12_020 = new(memory[nuk+31]) double[Lth];
    double * Q12_200 = new(memory[nuk+32]) double[Lth];

    int i[3];
    for(i[2] = 0; i[2]<lth[2]; i[2]++){
        for(i[1] = 0; i[1]<lth[1]; i[1]++){
            for(i[0] = 0; i[0]<=lth[0]; i[0]++){
                Q[ind2(0,1,NU,K)][ind(i,lth)] = LagrangeData[ind2(LInd[0],(L_center+i[0]-wth[0]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[1],(L_center+i[1]-wth[1]),(2*p+1),L_wth)]*LagrangeData[ind2(LInd[2],(L_center+i[2]-wth[2]),(2*p+1),L_wth)];
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 1; i[2]<(lth[2]-1); i[2]++){
        for(i[1] = 1; i[1]<(lth[1]-1); i[1]++){
            for(i[0] = 1; i[0]<(lth[0]-1); i[0]++){
                int I = ind(i,lth);
                d002[I] = DERIV002(Q[ind2(0,1,NU,K)],i,lth,dx);
                d020[I] = DERIV020(Q[ind2(0,1,NU,K)],i,lth,dx);
                d200[I] = DERIV200(Q[ind2(0,1,NU,K)],i,lth,dx);

                double d001 = DERIV001(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d010 = DERIV010(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d100 = DERIV100(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d011 = DERIV011(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d101 = DERIV101(Q[ind2(0,1,NU,K)],i,lth,dx);
                double d110 = DERIV110(Q[ind2(0,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis] ;
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,1,NU,K)][I] = coef[0][cI]*(d200[I])+coef[1][cI]*(d020[I])+coef[2][cI]*(d002[I])+coef[6][cI]*(d110)+coef[7][cI]*(d101)+coef[8][cI]*(d011)+coef[3][cI]*(d100)+coef[4][cI]*(d010)+coef[5][cI]*(d001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 2; i[2]<(lth[2]-2); i[2]++){
        for(i[1] = 2; i[1]<(lth[1]-2); i[1]++){
            for(i[0] = 2; i[0]<(lth[0]-2); i[0]++){
                int I = ind(i,lth);
                
                d004[I] = DERIV002(d002,i,lth,dx);
                d022[I] = DERIV002(d020,i,lth,dx);
                d040[I] = DERIV020(d020,i,lth,dx);
                d202[I] = DERIV200(d002,i,lth,dx);
                d220[I] = DERIV200(d020,i,lth,dx);
                d400[I] = DERIV200(d200,i,lth,dx);

                double d003 = DERIV001(d002,i,lth,dx);
                double d030 = DERIV010(d020,i,lth,dx);
                double d300 = DERIV100(d200,i,lth,dx);
                double d013 = DERIV011(d002,i,lth,dx);
                double d031 = DERIV011(d020,i,lth,dx);
                double d103 = DERIV101(d002,i,lth,dx);
                double d130 = DERIV110(d020,i,lth,dx);
                double d301 = DERIV101(d200,i,lth,dx);
                double d310 = DERIV110(d200,i,lth,dx);

                Q11_002[I] = DERIV002(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_020[I] = DERIV020(Q[ind2(1,1,NU,K)],i,lth,dx);
                Q11_200[I] = DERIV200(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_001 = DERIV001(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_010 = DERIV010(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_100 = DERIV100(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_011 = DERIV011(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_101 = DERIV101(Q[ind2(1,1,NU,K)],i,lth,dx);
                double Q11_110 = DERIV110(Q[ind2(1,1,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,2,NU,K)][I] = Q[ind2(1,1,NU,K)][I]+coef[0][cI]*((-1.0/12.0)*(dx[0]*dx[0])*d400[I])+coef[1][cI]*((-1.0/12.0)*(dx[1]*dx[1])*d040[I])+coef[2][cI]*((-1.0/12.0)*(dx[2]*dx[2])*d004[I])+coef[6][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*d130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d310)+coef[7][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d301)+coef[8][cI]*(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*d031)+coef[3][cI]*((-1.0/6.0)*(dx[0]*dx[0])*d300)+coef[4][cI]*((-1.0/6.0)*(dx[1]*dx[1])*d030)+coef[5][cI]*((-1.0/6.0)*(dx[2]*dx[2])*d003);
                Q[ind2(2,1,NU,K)][I] = coef[0][cI]*(Q11_200[I])+coef[1][cI]*(Q11_020[I])+coef[2][cI]*(Q11_002[I])+coef[6][cI]*(Q11_110)+coef[7][cI]*(Q11_101)+coef[8][cI]*(Q11_011)+coef[3][cI]*(Q11_100)+coef[4][cI]*(Q11_010)+coef[5][cI]*(Q11_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 3; i[2]<(lth[2]-3); i[2]++){
        for(i[1] = 3; i[1]<(lth[1]-3); i[1]++){
            for(i[0] = 3; i[0]<(lth[0]-3); i[0]++){
                int I = ind(i,lth);
                
                d006[I] = DERIV002(d004,i,lth,dx);
                d024[I] = DERIV002(d022,i,lth,dx);
                d042[I] = DERIV002(d040,i,lth,dx);
                d060[I] = DERIV020(d040,i,lth,dx);
                d204[I] = DERIV002(d202,i,lth,dx);
                d240[I] = DERIV020(d220,i,lth,dx);
                d402[I] = DERIV200(d202,i,lth,dx);
                d420[I] = DERIV200(d220,i,lth,dx);
                d600[I] = DERIV200(d400,i,lth,dx);

                double d005 = DERIV001(d004,i,lth,dx);
                double d050 = DERIV010(d040,i,lth,dx);
                double d500 = DERIV100(d400,i,lth,dx);
                double d015 = DERIV011(d004,i,lth,dx);
                double d033 = DERIV011(d022,i,lth,dx);
                double d051 = DERIV011(d040,i,lth,dx);
                double d105 = DERIV101(d004,i,lth,dx);
                double d150 = DERIV110(d040,i,lth,dx);
                double d303 = DERIV101(d202,i,lth,dx);
                double d330 = DERIV110(d220,i,lth,dx);
                double d501 = DERIV101(d400,i,lth,dx);
                double d510 = DERIV110(d400,i,lth,dx);

                Q11_004[I] = DERIV002(Q11_002,i,lth,dx);
                Q11_022[I] = DERIV002(Q11_020,i,lth,dx);
                Q11_040[I] = DERIV020(Q11_020,i,lth,dx);
                Q11_202[I] = DERIV200(Q11_002,i,lth,dx);
                Q11_220[I] = DERIV200(Q11_020,i,lth,dx);
                Q11_400[I] = DERIV200(Q11_200,i,lth,dx);
                Q21_002[I] = DERIV002(Q[ind2(2,1,NU,K)],i,lth,dx);
                Q21_020[I] = DERIV020(Q[ind2(2,1,NU,K)],i,lth,dx);
                Q21_200[I] = DERIV200(Q[ind2(2,1,NU,K)],i,lth,dx);
                Q12_002[I] = DERIV002(Q[ind2(1,2,NU,K)],i,lth,dx);
                Q12_020[I] = DERIV020(Q[ind2(1,2,NU,K)],i,lth,dx);
                Q12_200[I] = DERIV200(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q11_003 = DERIV001(Q11_002,i,lth,dx);
                double Q11_030 = DERIV010(Q11_020,i,lth,dx);
                double Q11_300 = DERIV100(Q11_200,i,lth,dx);
                double Q11_013 = DERIV011(Q11_002,i,lth,dx);
                double Q11_031 = DERIV011(Q11_020,i,lth,dx);
                double Q11_103 = DERIV101(Q11_002,i,lth,dx);
                double Q11_130 = DERIV110(Q11_020,i,lth,dx);
                double Q11_301 = DERIV101(Q11_200,i,lth,dx);
                double Q11_310 = DERIV110(Q11_200,i,lth,dx);

                double Q21_001 = DERIV001(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_010 = DERIV010(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_100 = DERIV100(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_011 = DERIV011(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_101 = DERIV101(Q[ind2(2,1,NU,K)],i,lth,dx);
                double Q21_110 = DERIV110(Q[ind2(2,1,NU,K)],i,lth,dx);

                double Q12_001 = DERIV001(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_010 = DERIV010(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_100 = DERIV100(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_011 = DERIV011(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_101 = DERIV101(Q[ind2(1,2,NU,K)],i,lth,dx);
                double Q12_110 = DERIV110(Q[ind2(1,2,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,3,NU,K)][I] = Q[ind2(1,2,NU,K)][I]+coef[0][cI]*((1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d600[I])+coef[1][cI]*((1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d060[I])+coef[2][cI]*((1.0/90.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d006[I])+coef[6][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d150+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*d330+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d510)+coef[7][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d105+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[2]*dx[2])*d303+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*d501)+coef[8][cI]*(1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d015+(-1.0/6.0)*(-1.0/6.0)*(dx[1]*dx[1])*(dx[2]*dx[2])*d033+(1.0/30.0)*1.0*(dx[1]*dx[1]*dx[1]*dx[1])*(1.0)*d051)+coef[3][cI]*((1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*d500)+coef[4][cI]*((1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*d050)+coef[5][cI]*((1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*d005);
                Q[ind2(2,2,NU,K)][I] = coef[0][cI]*(Q12_200[I]+((-1.0/12.0)*(dx[0]*dx[0])*Q11_400[I]))+coef[1][cI]*(Q12_020[I]+((-1.0/12.0)*(dx[1]*dx[1])*Q11_040[I]))+coef[2][cI]*(Q12_002[I]+((-1.0/12.0)*(dx[2]*dx[2])*Q11_004[I]))+coef[6][cI]*(Q12_110+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q11_130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_310))+coef[7][cI]*(Q12_101+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q11_301))+coef[8][cI]*(Q12_011+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q11_013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*Q11_031))+coef[3][cI]*(Q12_100+((-1.0/6.0)*(dx[0]*dx[0])*Q11_300))+coef[4][cI]*(Q12_010+((-1.0/6.0)*(dx[1]*dx[1])*Q11_030))+coef[5][cI]*(Q12_001+((-1.0/6.0)*(dx[2]*dx[2])*Q11_003));
                Q[ind2(3,1,NU,K)][I] = coef[0][cI]*(Q21_200[I])+coef[1][cI]*(Q21_020[I])+coef[2][cI]*(Q21_002[I])+coef[6][cI]*(Q21_110)+coef[7][cI]*(Q21_101)+coef[8][cI]*(Q21_011)+coef[3][cI]*(Q21_100)+coef[4][cI]*(Q21_010)+coef[5][cI]*(Q21_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    for(i[2] = 4; i[2]<(lth[2]-4); i[2]++){
        for(i[1] = 4; i[1]<(lth[1]-4); i[1]++){
            for(i[0] = 4; i[0]<(lth[0]-4); i[0]++){
                int I = ind(i,lth);
                
                double d008 = DERIV002(d006,i,lth,dx);
                double d080 = DERIV020(d060,i,lth,dx);
                double d800 = DERIV200(d600,i,lth,dx);

                double d007 = DERIV001(d006,i,lth,dx);
                double d070 = DERIV010(d060,i,lth,dx);
                double d700 = DERIV100(d600,i,lth,dx);
                double d017 = DERIV011(d006,i,lth,dx);
                double d035 = DERIV011(d024,i,lth,dx);
                double d053 = DERIV011(d042,i,lth,dx);
                double d071 = DERIV011(d060,i,lth,dx);
                double d107 = DERIV101(d006,i,lth,dx);
                double d170 = DERIV110(d060,i,lth,dx);
                double d305 = DERIV101(d204,i,lth,dx);
                double d350 = DERIV110(d240,i,lth,dx);
                double d503 = DERIV101(d402,i,lth,dx);
                double d530 = DERIV110(d420,i,lth,dx);
                double d701 = DERIV101(d600,i,lth,dx);
                double d710 = DERIV110(d600,i,lth,dx);

                double Q11_006 = DERIV002(Q11_004,i,lth,dx);
                double Q11_060 = DERIV020(Q11_040,i,lth,dx);
                double Q11_600 = DERIV200(Q11_400,i,lth,dx);
                double Q21_004 = DERIV002(Q21_002,i,lth,dx);
                double Q21_040 = DERIV020(Q21_020,i,lth,dx);
                double Q21_400 = DERIV200(Q21_200,i,lth,dx);
                double Q12_004 = DERIV002(Q12_002,i,lth,dx);
                double Q12_040 = DERIV020(Q12_020,i,lth,dx);
                double Q12_400 = DERIV200(Q12_200,i,lth,dx);
                double Q22_002 = DERIV002(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_020 = DERIV020(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_200 = DERIV200(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q31_002 = DERIV002(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_020 = DERIV020(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_200 = DERIV200(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q13_002 = DERIV002(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_020 = DERIV020(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_200 = DERIV200(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q11_005 = DERIV001(Q11_004,i,lth,dx);
                double Q11_050 = DERIV010(Q11_040,i,lth,dx);
                double Q11_500 = DERIV100(Q11_400,i,lth,dx);
                double Q11_015 = DERIV011(Q11_004,i,lth,dx);
                double Q11_033 = DERIV011(Q11_022,i,lth,dx);
                double Q11_051 = DERIV011(Q11_040,i,lth,dx);
                double Q11_105 = DERIV101(Q11_004,i,lth,dx);
                double Q11_150 = DERIV110(Q11_040,i,lth,dx);
                double Q11_303 = DERIV101(Q11_202,i,lth,dx);
                double Q11_330 = DERIV110(Q11_220,i,lth,dx);
                double Q11_501 = DERIV101(Q11_400,i,lth,dx);
                double Q11_510 = DERIV110(Q11_400,i,lth,dx);

                double Q21_003 = DERIV001(Q21_002,i,lth,dx);
                double Q21_030 = DERIV010(Q21_020,i,lth,dx);
                double Q21_300 = DERIV100(Q21_200,i,lth,dx);
                double Q21_013 = DERIV011(Q21_002,i,lth,dx);
                double Q21_031 = DERIV011(Q21_020,i,lth,dx);
                double Q21_103 = DERIV101(Q21_002,i,lth,dx);
                double Q21_130 = DERIV110(Q21_020,i,lth,dx);
                double Q21_301 = DERIV101(Q21_200,i,lth,dx);
                double Q21_310 = DERIV110(Q21_200,i,lth,dx);

                double Q12_003 = DERIV001(Q12_002,i,lth,dx);
                double Q12_030 = DERIV010(Q12_020,i,lth,dx);
                double Q12_300 = DERIV100(Q12_200,i,lth,dx);
                double Q12_013 = DERIV011(Q12_002,i,lth,dx);
                double Q12_031 = DERIV011(Q12_020,i,lth,dx);
                double Q12_103 = DERIV101(Q12_002,i,lth,dx);
                double Q12_130 = DERIV110(Q12_020,i,lth,dx);
                double Q12_301 = DERIV101(Q12_200,i,lth,dx);
                double Q12_310 = DERIV110(Q12_200,i,lth,dx);

                double Q22_001 = DERIV001(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_010 = DERIV010(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_100 = DERIV100(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_011 = DERIV011(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_101 = DERIV101(Q[ind2(2,2,NU,K)],i,lth,dx);
                double Q22_110 = DERIV110(Q[ind2(2,2,NU,K)],i,lth,dx);

                double Q31_001 = DERIV001(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_010 = DERIV010(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_100 = DERIV100(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_011 = DERIV011(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_101 = DERIV101(Q[ind2(3,1,NU,K)],i,lth,dx);
                double Q31_110 = DERIV110(Q[ind2(3,1,NU,K)],i,lth,dx);

                double Q13_001 = DERIV001(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_010 = DERIV010(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_100 = DERIV100(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_011 = DERIV011(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_101 = DERIV101(Q[ind2(1,3,NU,K)],i,lth,dx);
                double Q13_110 = DERIV110(Q[ind2(1,3,NU,K)],i,lth,dx);
                
                if(!cstCoef){
                    cInd[axis] = coefWth[axis] + i[axis] - wth[axis];
                    cInd[v0] = Ind[v0] -p + coefWth[v0] + i[v0] - wth[v0];
                    cInd[v1] = Ind[v1] - p + coefWth[v1] + i[v1] - wth[v1];
                    cI = ind(cInd, coefLth);
                }

                Q[ind2(1,4,NU,K)][I] = Q[ind2(1,3,NU,K)][I]+coef[0][cI]*((-1.0/560.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d800)+coef[1][cI]*((-1.0/560.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d080)+coef[2][cI]*((-1.0/560.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d008)+coef[6][cI]*((-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d170+(-1.0/6.0)*(1.0/30.0)*(dx[0]*dx[0])*(dx[1]*dx[1]*dx[1]*dx[1])*d350+(1.0/30.0)*(-1.0/6.0)*(dx[0]*dx[0]*dx[0]*dx[0])*(dx[1]*dx[1])*d530+(-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d710)+coef[7][cI]*((-1.0/140.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d107+(-1.0/6.0)*(1.0/30.0)*(dx[0]*dx[0])*(dx[2]*dx[2]*dx[2]*dx[2])*d305+(1.0/30.0)*(-1.0/6.0)*(dx[0]*dx[0]*dx[0]*dx[0])*(dx[2]*dx[2])*d503+(-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d701)+coef[8][cI]*(1.0*(-1.0/140.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d017+(-1.0/6.0)*(1.0/30.0)*(dx[1]*dx[1])*(dx[2]*dx[2]*dx[2]*dx[2])*d035+(1.0/30.0)*(-1.0/6.0)*(dx[1]*dx[1]*dx[1]*dx[1])*(dx[2]*dx[2])*d053+(-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d071)+coef[3][cI]*((-1.0/140.0)*(dx[0]*dx[0]*dx[0]*dx[0]*dx[0]*dx[0])*d700)+coef[4][cI]*((-1.0/140.0)*(dx[1]*dx[1]*dx[1]*dx[1]*dx[1]*dx[1])*d070)+coef[5][cI]*((-1.0/140.0)*(dx[2]*dx[2]*dx[2]*dx[2]*dx[2]*dx[2])*d007);
                Q[ind2(2,3,NU,K)][I]  = coef[0][cI]*(Q13_200+((-1.0/12.0)*(dx[0]*dx[0])*Q12_400+(1.0/90.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_600))+coef[1][cI]*(Q13_020+((-1.0/12.0)*(dx[1]*dx[1])*Q12_040+(1.0/90.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_060))+coef[2][cI]*(Q13_002+((-1.0/12.0)*(dx[2]*dx[2])*Q12_004+(1.0/90.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_006))+coef[6][cI]*(Q13_110+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q12_130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q12_310+1.0*(1.0/30.0)*(1.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_150+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[1]*dx[1])*Q11_330+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*Q11_510))+coef[7][cI]*(Q13_101+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q12_103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q12_301+1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_105+(-1.0/6.0)*(-1.0/6.0)*(dx[0]*dx[0])*(dx[2]*dx[2])*Q11_303+(1.0/30.0)*1.0*(dx[0]*dx[0]*dx[0]*dx[0])*(1.0)*Q11_501))+coef[8][cI]*(Q13_011+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q12_013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*Q12_031+1.0*(1.0/30.0)*(1.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_015+(-1.0/6.0)*(-1.0/6.0)*(dx[1]*dx[1])*(dx[2]*dx[2])*Q11_033+(1.0/30.0)*1.0*(dx[1]*dx[1]*dx[1]*dx[1])*(1.0)*Q11_051))+coef[3][cI]*(Q13_100+((-1.0/6.0)*(dx[0]*dx[0])*Q12_300+(1.0/30.0)*(dx[0]*dx[0]*dx[0]*dx[0])*Q11_500))+coef[4][cI]*(Q13_010+((-1.0/6.0)*(dx[1]*dx[1])*Q12_030+(1.0/30.0)*(dx[1]*dx[1]*dx[1]*dx[1])*Q11_050))+coef[5][cI]*(Q13_001+((-1.0/6.0)*(dx[2]*dx[2])*Q12_003+(1.0/30.0)*(dx[2]*dx[2]*dx[2]*dx[2])*Q11_005));
                Q[ind2(3,2,NU,K)][I]  = coef[0][cI]*(Q22_200+((-1.0/12.0)*(dx[0]*dx[0])*Q21_400))+coef[1][cI]*(Q22_020+((-1.0/12.0)*(dx[1]*dx[1])*Q21_040))+coef[2][cI]*(Q22_002+((-1.0/12.0)*(dx[2]*dx[2])*Q21_004))+coef[6][cI]*(Q22_110+(1.0*(-1.0/6.0)*(1.0)*(dx[1]*dx[1])*Q21_130+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q21_310))+coef[7][cI]*(Q22_101+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q21_103+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*Q21_301))+coef[8][cI]*(Q22_011+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*Q21_013+(-1.0/6.0)*1.0*(dx[1]*dx[1])*(1.0)*Q21_031))+coef[3][cI]*(Q22_100+((-1.0/6.0)*(dx[0]*dx[0])*Q21_300))+coef[4][cI]*(Q22_010+((-1.0/6.0)*(dx[1]*dx[1])*Q21_030))+coef[5][cI]*(Q22_001+((-1.0/6.0)*(dx[2]*dx[2])*Q21_003));
                Q[ind2(4,1,NU,K)][I]  = coef[0][cI]*(Q31_200)+coef[1][cI]*(Q31_020)+coef[2][cI]*(Q31_002)+coef[6][cI]*(Q31_110)+coef[7][cI]*(Q31_101)+coef[8][cI]*(Q31_011)+coef[3][cI]*(Q31_100)+coef[4][cI]*(Q31_010)+coef[5][cI]*(Q31_001);
            } // end i[0] loop
        }// end of i[1] loop
    }// end of i[2] loop
    
    /* Calculate \partial_z^{mu1}\partial_y^{mu0}Q^{nu}L_iL_jL_k evaluted at boundary point */
    
    int MU = 2*p+1;
    int I = ind(wth,lth);
    Z[ind3(0,0,0,NU,MU,MU)] = Q[ind2(0,1,NU,K)][I];
    Z[ind3(1,0,0,NU,MU,MU)] = Q[ind2(1,4,NU,K)][I];
    Z[ind3(2,0,0,NU,MU,MU)] = Q[ind2(2,3,NU,K)][I];
    Z[ind3(3,0,0,NU,MU,MU)] = Q[ind2(3,2,NU,K)][I];
    Z[ind3(4,0,0,NU,MU,MU)] = Q[ind2(4,1,NU,K)][I];
    
    /* Generate the tangential derivatives here */
#include "tangentialDerivs_3D_8_2.C"
}


LagrangeDerivFun pickLagrangeDerivFun(int dim, int p, int axis){
    if(p == 1 && dim == 2 && axis == 0){
        return &LagrangeDeriv_2D_2_0;
    }
    else if(p == 1 && dim == 2 && axis == 1){
        return &LagrangeDeriv_2D_2_1;
    }
    else if(p == 2 && dim == 2 && axis == 0){
        return &LagrangeDeriv_2D_4_0;
    }
    else if(p == 2 && dim == 2 && axis == 1){
        return &LagrangeDeriv_2D_4_1;
    }
    else if(p == 3 && dim == 2 && axis == 0){
        return &LagrangeDeriv_2D_6_0;
    }
    else if(p == 3 && dim == 2 && axis == 1){
        return &LagrangeDeriv_2D_6_1;
    }
    else if(p == 4 && dim == 2 && axis == 0){
        return &LagrangeDeriv_2D_8_0;
    }
    else if(p == 4 && dim == 2 && axis == 1){
        return &LagrangeDeriv_2D_8_1;
    }
    else if(p == 1 && dim == 3 && axis == 0){
        return &LagrangeDeriv_3D_2_0;
    }
    else if(p == 1 && dim == 3 && axis == 1){
        return &LagrangeDeriv_3D_2_1;
    }
    else if(p == 1 && dim == 3 && axis == 2){
        return &LagrangeDeriv_3D_2_2;
    }
    else if(p == 2 && dim == 3 && axis == 0){
        return &LagrangeDeriv_3D_4_0;
    }
    else if(p == 2 && dim == 3 && axis == 1){
        return &LagrangeDeriv_3D_4_1;
    }
    else if(p == 2 && dim == 3 && axis == 2){
        return &LagrangeDeriv_3D_4_2;
    }
    else if(p == 3 && dim == 3 && axis == 0){
        return &LagrangeDeriv_3D_6_0;
    }
    else if(p == 3 && dim == 3 && axis == 1){
        return &LagrangeDeriv_3D_6_1;
    }
    else if(p == 3 && dim == 3 && axis == 2){
        return &LagrangeDeriv_3D_6_2;
    }
    else if(p == 4 && dim == 3 && axis == 0){
        return &LagrangeDeriv_3D_8_0;
    }
    else if(p == 4 && dim == 3 && axis == 1){
        return &LagrangeDeriv_3D_8_1;
    }
    else if(p == 4 && dim == 3 && axis == 2){
        return &LagrangeDeriv_3D_8_2;
    }
    else{
        printf("ERROR in pickLagrangeDerivFun: LagrangeDerivFun does not exist\n");
        exit(-1);
    }
}
