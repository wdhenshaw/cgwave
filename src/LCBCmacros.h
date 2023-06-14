#ifndef LCBCmacros_h
#define LCBCmacros_h

/* Define the macros for indices that need to be used */
#define ind2(i,j,n1,n2) (((j)*(n1))+(i))
#define ind3(i,j,k,n1,n2,n3) (i+(j*n1)+(k*n1*n2))
#define ind4(i,j,k,l,n1,n2,n3,n4) (i+(j*n1)+(k*n1*n2) + (l*n1*n2*n3))
#define ind5(i,j,k,l,m,n1,n2,n3,n4,n5) (i+(j*n1)+(k*n1*n2) + (l*n1*n2*n3) + (m*n1*n2*n3*n4))
#define ind(I,N) (I[0]+(I[1]*N[0])+(I[2]*N[0]*N[1]))

/* Define the min and max macros */
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

/* A macro to sum two 3D vectors */
#define sumVectors(v1,v2) {(v1[0]+v2[0]),(v1[1]+v2[1]),(v1[2]+v2[2])}

/* A macro to determine the order of accuracy to discretize auxiliary equations based on derivative order d and the order in space (p = orderInSpace/2) */
#define getOrder(m,p) ((m%2==0) ? (p-((m-2)/2)) : (p-((m-1)/2)))

/* Macros to determine if a grid point in the stencil is inside or outside a given boundary. This depends on the number of known axes (fixed axes where boundary data is known).
   k runs from 0 to (2p + 1) */
#define COND(k,side) ((side == 0) ? (k<p) : (k>p))
#define CONDITION3(k,side) (COND(k[0],side[0])||COND(k[1],side[1])||COND(k[2],side[2]))
#define CONDITION2(k,axis,side) (COND(k[axis[0]],side[0])||COND(k[axis[1]],side[1]))
#define CONDITION1(k,axis,side) (COND(k[axis[0]],side[0]))
#define CONDITION(k,axis,side,axisNum) ((axisNum == 1)?(CONDITION1(k,axis,side)):((axisNum == 2)?(CONDITION2(k,axis,side)):(CONDITION3(k,side))))

/* A macro to that sets a value to some n1 if side = 0 or n2 otherwise */
#define sideBasedValue(side,n1,n2) ((side == 0) ? (n1) : (n2))

/* A macro to that sets a value to some n1 if dim = 2 or n2 otherwise */
#define dimBasedValue(dim,n1,n2) ((dim == 2) ? (n1) : (n2))

/* Compute extra ghost points used by the user if any */
#define EXTRA_GHOST (userNumGhost - numGhost)

/* Define an indexing for the solution provided by the user while ensuring it matches the grid used in the LCBC class */
#define solInd(I,N) ((I[0] + EXTRA_GHOST) + (I[1] + EXTRA_GHOST)*(N[0] + 2*EXTRA_GHOST)+ dimBasedValue(dim,0,((I[2] + EXTRA_GHOST)*(N[0] + 2*EXTRA_GHOST)*(N[1] + 2*EXTRA_GHOST))))

/* A macro to that sets a value to some n1 if axis = fixedAxis or n2 otherwise */
#define axisBasedVal(axis,fixedAxis,n1,n2) ((axis==fixedAxis)?(n1):(n2))

/* Define an indexing for the boundary grid functions if provided by the user.
   Note: I adjusted the LCBC class such that boundary information on the grid is no longer constrained by the grid provided by the user. So this is no longer necessary but still used in the code... */
#define bdryInd(I,N,axis) (ind(I,N))

/* Define a type of functions to be used when choosing a specific function to compute Lagrange derivatives.
   The new implementation uses separate functions for separate cases rather than one general function. */
typedef void (*LagrangeDerivFun)(double[], int*, int*, double*, int, double**, int[3], int[3], bool, double[3], int, int, int, char **&);

#endif /* LCBCmacros_h */
