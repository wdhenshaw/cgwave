#ifndef AUGMENTED_KRYLOV_H
#define AUGMENTED_KRYLOV_H

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

class AugmentedKrylov
{

public:

enum KrylovTypesEnum
{
   gmres=0,
   conjugateGradient,
   biConjugateGradientStabilized
};

AugmentedKrylov();
~AugmentedKrylov();

int getNumberOfIterations() const;

int getNumberOfMatrixVectorProducts() const;

// return the residual from the last solve
Real getResidual() const;

RealArray & getResidualVector() const;

int setKrylovType( const KrylovTypesEnum & krylovType );

// Use the matrix A 
Real solve( const RealArray & A, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, RealArray & x );

// use a matrix-vector multiply function
Real solve( MatVectFunctionPtr matVectFunction, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, RealArray & x );


// ----- Utility functions ----
static Real innerProduct( const RealArray & x, const RealArray & y );

static void matVect(const RealArray & A, const RealArray & x, RealArray & y, int transpose=0, int m0=-1, int n0=-1 );

static void matVect(Real **pA, const RealArray & x, RealArray & y, int transpose =0, int m0=-1, int n0=-1 );

static Real norm( const RealArray & x );

// supply eigenvalues if augmented vectors are eigenvectors
int setAugmentedEigenvalues( const RealArray & augEigs );

protected:

// Generic Routines taking all arguments 
Real solveCG( const RealArray & A, MatVectFunctionPtr matVectFunction, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, RealArray & x );
Real solveBiCGStab( const RealArray & A, MatVectFunctionPtr matVectFunction, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, RealArray & x );
Real solveGmres( const RealArray & A, MatVectFunctionPtr matVectFunction, const RealArray & b, const RealArray & x0, const RealArray & W, const int maxit, const Real tol, RealArray & x );


public:
   
// The database is used to hold parameters
mutable DataBase dbase;


};



#endif  
