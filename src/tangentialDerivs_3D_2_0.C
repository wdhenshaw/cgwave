double d20[Lth], d02[Lth], d22[Lth], d00[Lth];
double d10, d01, d11, d21, d12;
int j[3] = {wth[0],wth[1],wth[2]};

int nu = 0, k = 1;
for(j[2] = (wth[2]-1); j[2]<=(wth[2]+1); j[2]++){
for(j[1] = (wth[1]-1); j[1]<=(wth[1]+1); j[1]++){
int J = ind(j,lth);
d00[J] = Q[ind2(nu,k,NU,K)][J];
d20[J] = DERIV020(Q[ind2(nu,k,NU,K)], j, lth, dx);
d02[J] = DERIV002(Q[ind2(nu,k,NU,K)], j, lth, dx);
d22[J] = DERIV022(Q[ind2(nu,k,NU,K)], j, lth, dx);
}
}

int J = ind(wth,lth);

d10=DERIV010(d00,wth,lth,dx);
d01=DERIV001(d00,wth,lth,dx);
d11=DERIV011(d00,wth,lth,dx);
d21=DERIV001(d20,wth,lth,dx);
d12=DERIV010(d02,wth,lth,dx);
Z[ind3(0,1,0,NU,MU,MU)] = d10;
Z[ind3(0,2,0,NU,MU,MU)] = d20[J];
Z[ind3(0,0,1,NU,MU,MU)] = d01;
Z[ind3(0,1,1,NU,MU,MU)] = d11;
Z[ind3(0,2,1,NU,MU,MU)] = d21;
Z[ind3(0,0,2,NU,MU,MU)] = d02[J];
Z[ind3(0,1,2,NU,MU,MU)] = d12;
Z[ind3(0,2,2,NU,MU,MU)] = d22[J];

nu = 1; k = 1;

for(j[2] = (wth[2]-1); j[2]<=(wth[2]+1); j[2]++){
for(j[1] = (wth[1]-1); j[1]<=(wth[1]+1); j[1]++){
int J = ind(j,lth);
d00[J] = Q[ind2(nu,k,NU,K)][J];
d20[J] = DERIV020(Q[ind2(nu,k,NU,K)], j, lth, dx);
d02[J] = DERIV002(Q[ind2(nu,k,NU,K)], j, lth, dx);
d22[J] = DERIV022(Q[ind2(nu,k,NU,K)], j, lth, dx);
}
}

J = ind(wth,lth);

d10=DERIV010(d00,wth,lth,dx);
d01=DERIV001(d00,wth,lth,dx);
d11=DERIV011(d00,wth,lth,dx);
d21=DERIV001(d20,wth,lth,dx);
d12=DERIV010(d02,wth,lth,dx);
Z[ind3(1,1,0,NU,MU,MU)] = d10;
Z[ind3(1,2,0,NU,MU,MU)] = d20[J];
Z[ind3(1,0,1,NU,MU,MU)] = d01;
Z[ind3(1,1,1,NU,MU,MU)] = d11;
Z[ind3(1,2,1,NU,MU,MU)] = d21;
Z[ind3(1,0,2,NU,MU,MU)] = d02[J];
Z[ind3(1,1,2,NU,MU,MU)] = d12;
Z[ind3(1,2,2,NU,MU,MU)] = d22[J];
