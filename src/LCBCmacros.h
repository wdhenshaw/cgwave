#ifndef LCBCmacros_h
#define LCBCmacros_h

#define ind2(i,j,n1,n2) (((j)*(n1))+(i))
#define ind3(i,j,k,n1,n2,n3) (i+(j*n1)+(k*n1*n2))
#define ind4(i,j,k,l,n1,n2,n3,n4) (i+(j*n1)+(k*n1*n2) + (l*n1*n2*n3))
#define ind5(i,j,k,l,m,n1,n2,n3,n4,n5) (i+(j*n1)+(k*n1*n2) + (l*n1*n2*n3) + (m*n1*n2*n3*n4))
#define ind(I,N) (I[0]+(I[1]*N[0])+(I[2]*N[0]*N[1]))

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y));
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y));

#define sumVectors(v1,v2) {(v1[0]+v2[0]),(v1[1]+v2[1]),(v1[2]+v2[2])}
#define getOrder(m,p) ((m%2==0) ? (p-((m-2)/2)) : (p-((m-1)/2)))
#define COND(k,side) ((side == 0) ? (k<p) : (k>p))
#define sideBasedValue(side,n1,n2) ((side == 0) ? (n1) : (n2))

#define dimBasedValue(dim,n1,n2) ((dim == 2) ? (n1) : (n2))


#endif /* LCBCmacros_h */
