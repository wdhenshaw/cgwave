double d20[Lth], d02[Lth], d22[Lth], d00[Lth];
double d10,d30,d40,d01,d11,d21,d31,d41,d12,d32,d42,d03,d13,d23,d33,d43,d04,d14,d24,d34,d44;
int j[3] = {wth[0],wth[1],wth[2]};


int nu = 0, k = 1;
for(j[2] = (wth[2]-1); j[2]<=(wth[2]+1); j[2]++){
for(j[0] = (wth[0]-1); j[0]<=(wth[0]+1); j[0]++){
int J = ind(j,lth);
d00[J] = Q[ind2(nu,k,NU,K)][J];
d20[J] = DERIV200(Q[ind2(nu,k,NU,K)], j, lth, dx);
d02[J] = DERIV002(Q[ind2(nu,k,NU,K)], j, lth, dx);
d22[J] = DERIV202(Q[ind2(nu,k,NU,K)], j, lth, dx);
}
}

 int J = ind(wth,lth);

d10=DERIV100(d00,wth,lth,dx);
d30=DERIV100(d20,wth,lth,dx);
d40=DERIV200(d20,wth,lth,dx);
d01=DERIV001(d00,wth,lth,dx);
d11=DERIV101(d00,wth,lth,dx);
d21=DERIV001(d20,wth,lth,dx);
d31=DERIV101(d20,wth,lth,dx);
d41=DERIV201(d20,wth,lth,dx);
d12=DERIV100(d02,wth,lth,dx);
d32=DERIV100(d22,wth,lth,dx);
d42=DERIV200(d22,wth,lth,dx);
d03=DERIV001(d02,wth,lth,dx);
d13=DERIV101(d02,wth,lth,dx);
d23=DERIV001(d22,wth,lth,dx);
d33=DERIV101(d22,wth,lth,dx);
d43=DERIV201(d22,wth,lth,dx);
d04=DERIV002(d02,wth,lth,dx);
d14=DERIV102(d02,wth,lth,dx);
d24=DERIV002(d22,wth,lth,dx);
d34=DERIV102(d22,wth,lth,dx);
d44=DERIV202(d22,wth,lth,dx);
Z[ind3(0,1,0,NU,MU,MU)] = d10+((-1.0/6.0)*(dx[0]*dx[0])*d30);
Z[ind3(0,2,0,NU,MU,MU)] = d20[J]+((-1.0/12.0)*(dx[0]*dx[0])*d40);
Z[ind3(0,3,0,NU,MU,MU)] = d30;
Z[ind3(0,4,0,NU,MU,MU)] = d40;
Z[ind3(0,0,1,NU,MU,MU)] = d01+((-1.0/6.0)*(dx[2]*dx[2])*d03);
Z[ind3(0,1,1,NU,MU,MU)] = d11+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d31);
Z[ind3(0,2,1,NU,MU,MU)] = d21+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d23+(-1.0/12.0)*1.0*(dx[0]*dx[0])*(1.0)*d41);
Z[ind3(0,3,1,NU,MU,MU)] = d31+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d33);
Z[ind3(0,4,1,NU,MU,MU)] = d41+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d43);
Z[ind3(0,0,2,NU,MU,MU)] = d02[J]+((-1.0/12.0)*(dx[2]*dx[2])*d04);
Z[ind3(0,1,2,NU,MU,MU)] = d12+(1.0*(-1.0/12.0)*(1.0)*(dx[2]*dx[2])*d14+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d32);
Z[ind3(0,2,2,NU,MU,MU)] = d22[J]+(1.0*(-1.0/12.0)*(1.0)*(dx[2]*dx[2])*d24+(-1.0/12.0)*1.0*(dx[0]*dx[0])*(1.0)*d42);
Z[ind3(0,3,2,NU,MU,MU)] = d32+(1.0*(-1.0/12.0)*(1.0)*(dx[2]*dx[2])*d34);
Z[ind3(0,4,2,NU,MU,MU)] = d42+(1.0*(-1.0/12.0)*(1.0)*(dx[2]*dx[2])*d44);
Z[ind3(0,0,3,NU,MU,MU)] = d03;
Z[ind3(0,1,3,NU,MU,MU)] = d13+(1.0*(-1.0/6.0)*(1.0)*(dx[0]*dx[0])*d33);
Z[ind3(0,2,3,NU,MU,MU)] = d23+(1.0*(-1.0/12.0)*(1.0)*(dx[0]*dx[0])*d43);
Z[ind3(0,3,3,NU,MU,MU)] = d33;
Z[ind3(0,4,3,NU,MU,MU)] = d43;
Z[ind3(0,0,4,NU,MU,MU)] = d04;
Z[ind3(0,1,4,NU,MU,MU)] = d14+(1.0*(-1.0/6.0)*(1.0)*(dx[0]*dx[0])*d34);
Z[ind3(0,2,4,NU,MU,MU)] = d24+(1.0*(-1.0/12.0)*(1.0)*(dx[0]*dx[0])*d44);
Z[ind3(0,3,4,NU,MU,MU)] = d34;
Z[ind3(0,4,4,NU,MU,MU)] = d44;

nu = 1; k = 2;
for(j[2] = (wth[2]-1); j[2]<=(wth[2]+1); j[2]++){
for(j[0] = (wth[0]-1); j[0]<=(wth[0]+1); j[0]++){
int J = ind(j,lth);
d00[J] = Q[ind2(nu,k,NU,K)][J];
d20[J] = DERIV200(Q[ind2(nu,k,NU,K)], j, lth, dx);
d02[J] = DERIV002(Q[ind2(nu,k,NU,K)], j, lth, dx);
d22[J] = DERIV202(Q[ind2(nu,k,NU,K)], j, lth, dx);
}
}

J = ind(wth,lth);

d10=DERIV100(d00,wth,lth,dx);
d30=DERIV100(d20,wth,lth,dx);
d40=DERIV200(d20,wth,lth,dx);
d01=DERIV001(d00,wth,lth,dx);
d11=DERIV101(d00,wth,lth,dx);
d21=DERIV001(d20,wth,lth,dx);
d31=DERIV101(d20,wth,lth,dx);
d41=DERIV201(d20,wth,lth,dx);
d12=DERIV100(d02,wth,lth,dx);
d32=DERIV100(d22,wth,lth,dx);
d42=DERIV200(d22,wth,lth,dx);
d03=DERIV001(d02,wth,lth,dx);
d13=DERIV101(d02,wth,lth,dx);
d23=DERIV001(d22,wth,lth,dx);
d33=DERIV101(d22,wth,lth,dx);
d43=DERIV201(d22,wth,lth,dx);
d04=DERIV002(d02,wth,lth,dx);
d14=DERIV102(d02,wth,lth,dx);
d24=DERIV002(d22,wth,lth,dx);
d34=DERIV102(d22,wth,lth,dx);
d44=DERIV202(d22,wth,lth,dx);
Z[ind3(1,1,0,NU,MU,MU)] = d10+((-1.0/6.0)*(dx[0]*dx[0])*d30);
Z[ind3(1,2,0,NU,MU,MU)] = d20[J]+((-1.0/12.0)*(dx[0]*dx[0])*d40);
Z[ind3(1,3,0,NU,MU,MU)] = d30;
Z[ind3(1,4,0,NU,MU,MU)] = d40;
Z[ind3(1,0,1,NU,MU,MU)] = d01+((-1.0/6.0)*(dx[2]*dx[2])*d03);
Z[ind3(1,1,1,NU,MU,MU)] = d11+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d13+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d31);
Z[ind3(1,2,1,NU,MU,MU)] = d21+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d23+(-1.0/12.0)*1.0*(dx[0]*dx[0])*(1.0)*d41);
Z[ind3(1,3,1,NU,MU,MU)] = d31+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d33);
Z[ind3(1,4,1,NU,MU,MU)] = d41+(1.0*(-1.0/6.0)*(1.0)*(dx[2]*dx[2])*d43);
Z[ind3(1,0,2,NU,MU,MU)] = d02[J]+((-1.0/12.0)*(dx[2]*dx[2])*d04);
Z[ind3(1,1,2,NU,MU,MU)] = d12+(1.0*(-1.0/12.0)*(1.0)*(dx[2]*dx[2])*d14+(-1.0/6.0)*1.0*(dx[0]*dx[0])*(1.0)*d32);
Z[ind3(1,2,2,NU,MU,MU)] = d22[J]+(1.0*(-1.0/12.0)*(1.0)*(dx[2]*dx[2])*d24+(-1.0/12.0)*1.0*(dx[0]*dx[0])*(1.0)*d42);
Z[ind3(1,3,2,NU,MU,MU)] = d32+(1.0*(-1.0/12.0)*(1.0)*(dx[2]*dx[2])*d34);
Z[ind3(1,4,2,NU,MU,MU)] = d42+(1.0*(-1.0/12.0)*(1.0)*(dx[2]*dx[2])*d44);
Z[ind3(1,0,3,NU,MU,MU)] = d03;
Z[ind3(1,1,3,NU,MU,MU)] = d13+(1.0*(-1.0/6.0)*(1.0)*(dx[0]*dx[0])*d33);
Z[ind3(1,2,3,NU,MU,MU)] = d23+(1.0*(-1.0/12.0)*(1.0)*(dx[0]*dx[0])*d43);
Z[ind3(1,3,3,NU,MU,MU)] = d33;
Z[ind3(1,4,3,NU,MU,MU)] = d43;
Z[ind3(1,0,4,NU,MU,MU)] = d04;
Z[ind3(1,1,4,NU,MU,MU)] = d14+(1.0*(-1.0/6.0)*(1.0)*(dx[0]*dx[0])*d34);
Z[ind3(1,2,4,NU,MU,MU)] = d24+(1.0*(-1.0/12.0)*(1.0)*(dx[0]*dx[0])*d44);
Z[ind3(1,3,4,NU,MU,MU)] = d34;
Z[ind3(1,4,4,NU,MU,MU)] = d44;

nu = 2; k = 1;
for(j[2] = (wth[2]-1); j[2]<=(wth[2]+1); j[2]++){
for(j[0] = (wth[0]-1); j[0]<=(wth[0]+1); j[0]++){
int J = ind(j,lth);
d00[J] = Q[ind2(nu,k,NU,K)][J];
d20[J] = DERIV200(Q[ind2(nu,k,NU,K)], j, lth, dx);
d02[J] = DERIV002(Q[ind2(nu,k,NU,K)], j, lth, dx);
d22[J] = DERIV202(Q[ind2(nu,k,NU,K)], j, lth, dx);
}
}

J = ind(wth,lth);

d10=DERIV100(d00,wth,lth,dx);
d30=DERIV100(d20,wth,lth,dx);
d40=DERIV200(d20,wth,lth,dx);
d01=DERIV001(d00,wth,lth,dx);
d11=DERIV101(d00,wth,lth,dx);
d21=DERIV001(d20,wth,lth,dx);
d31=DERIV101(d20,wth,lth,dx);
d41=DERIV201(d20,wth,lth,dx);
d12=DERIV100(d02,wth,lth,dx);
d32=DERIV100(d22,wth,lth,dx);
d42=DERIV200(d22,wth,lth,dx);
d03=DERIV001(d02,wth,lth,dx);
d13=DERIV101(d02,wth,lth,dx);
d23=DERIV001(d22,wth,lth,dx);
d33=DERIV101(d22,wth,lth,dx);
d43=DERIV201(d22,wth,lth,dx);
d04=DERIV002(d02,wth,lth,dx);
d14=DERIV102(d02,wth,lth,dx);
d24=DERIV002(d22,wth,lth,dx);
d34=DERIV102(d22,wth,lth,dx);
d44=DERIV202(d22,wth,lth,dx);
Z[ind3(2,1,0,NU,MU,MU)] = d10;
Z[ind3(2,2,0,NU,MU,MU)] = d20[J];
Z[ind3(2,3,0,NU,MU,MU)] = d30;
Z[ind3(2,4,0,NU,MU,MU)] = d40;
Z[ind3(2,0,1,NU,MU,MU)] = d01;
Z[ind3(2,1,1,NU,MU,MU)] = d11;
Z[ind3(2,2,1,NU,MU,MU)] = d21;
Z[ind3(2,3,1,NU,MU,MU)] = d31;
Z[ind3(2,4,1,NU,MU,MU)] = d41;
Z[ind3(2,0,2,NU,MU,MU)] = d02[J];
Z[ind3(2,1,2,NU,MU,MU)] = d12;
Z[ind3(2,2,2,NU,MU,MU)] = d22[J];
Z[ind3(2,3,2,NU,MU,MU)] = d32;
Z[ind3(2,4,2,NU,MU,MU)] = d42;
Z[ind3(2,0,3,NU,MU,MU)] = d03;
Z[ind3(2,1,3,NU,MU,MU)] = d13;
Z[ind3(2,2,3,NU,MU,MU)] = d23;
Z[ind3(2,3,3,NU,MU,MU)] = d33;
Z[ind3(2,4,3,NU,MU,MU)] = d43;
Z[ind3(2,0,4,NU,MU,MU)] = d04;
Z[ind3(2,1,4,NU,MU,MU)] = d14;
Z[ind3(2,2,4,NU,MU,MU)] = d24;
Z[ind3(2,3,4,NU,MU,MU)] = d34;
Z[ind3(2,4,4,NU,MU,MU)] = d44;
