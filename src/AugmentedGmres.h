#ifndef AUGMENTED_GMRES_H
#define AUGMENTED_GMRES_H

// #include "mpi.h"
#include "Overture.h"
// #include "Oges.h"


#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

// forward declarations
// class CgWave;
// class GraphicsParameters;

// Function prototype for the matrix-vector function used for the matrix free AuGmes
typedef void (*MatVectFunctionPtr)( const RealArray & x, RealArray & y );

class AugmentedGmres
{

public:

AugmentedGmres();
~AugmentedGmres();

int getNumberOfIterations() const;

// return the residual from the last solve
Real getResidual() const;

RealArray & getResidualVector() const;

// Use the matrix A 
Real solve( const RealArray & A, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, RealArray & x );

// use a matrix-vector multiply function
Real solve( MatVectFunctionPtr matVectFunction, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, RealArray & x );


// ----- Utility functions ----
static Real innerProduct( const RealArray & x, const RealArray & y );

static void matVect(const RealArray & A, const RealArray & x, RealArray & y, int transpose=0, int m0=-1, int n0=-1 );

static Real norm( const RealArray & x );

// supply eigenvalues if augmented vectors are eigenvectors
int setAugmentedEigenvalues( const RealArray & augEigs );

protected:

// Generic Routine 
Real solve( const RealArray & A, MatVectFunctionPtr matVectFunction, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, RealArray & x );


// The database is used to hold parameters
mutable DataBase dbase;


};



#endif  
