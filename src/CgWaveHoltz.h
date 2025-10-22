#ifndef CG_WAVE_HOLTZ_H
#define CG_WAVE_HOLTZ_H

// added Oct 21, 2205
#ifdef CGWAVE_USE_PETSC
#include "mpi.h"
#endif

#include "Overture.h"
#include "Oges.h"
#include "PlotStuff.h"


// // krb do not use extern "C" if PETSc is linked using BOPT=g_c++
// extern "C"
// {
// // *wdh* 2015/09/31  To avoid having PETSc include complex.h do this: 
// #include "petscconf.h"
// #undef PETSC_HAVE_CXX_COMPLEX
// #include "petscksp.h"
// }

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

// forward declarations
class CgWave;
class GraphicsParameters;

class CgWaveHoltz
{

public:

CgWaveHoltz(CompositeGrid & cg, GenericGraphicsInterface & gi);
~CgWaveHoltz();

int assignEigenSolverInitialCondition( bool smoothInitialCondition );

// Initialize time-step and forcing 
int initialize();

// Assign parameters 
int interactiveUpdate();

int outputHeader();

// output tables for eigenWave
int outputEigenTable();

// save results to a check file
int outputCheckFile( const RealArray & maxResArray, const Real errorBetweenWaveHoltzAndHelmholtz );

// save results to a Matlab file
int outputMatlabFile( const aString & matlabFileName );

int plotEigenVectors( realCompositeGridFunction & eigenVectors, const RealArray & eig, const aString & name, GL_GraphicsInterface & ps, GraphicsParameters & psp );

// plot WaveHoltz filter beta and any eigenvalues 
// int plotFilter( RealArray & eigenValues, GL_GraphicsInterface & ps, PlotStuffParameters & psp );

// Compute the residual in the current solution
real residual( RealArray & maxRes, int useAdjustedOmega = 2 );

// Compute the residual a grid function
real residual( RealCompositeGridFunction & uh, RealArray & maxRes, int useAdjustedOmega = 2 );
real residual( RealCompositeGridFunction & uh, RealCompositeGridFunction & f, RealArray & maxRes, int useAdjustedOmega = 2 );

// save check file -- used for regression and convergence tests
int saveCheckFile( int checkFileCounter, Real maxErr, Real solutionNorm );

int setNameOfGridFile( aString & nameOfOGFile );

// Setup grids and grid functions
int setup();

int smoothEigenSolverInitialCondition();

// WaveHoltz iteration:
int solve();

// Solve Helmholtz with WaveHoltz and a Krylov method augmented with eigenvectors 
int solveAugmentedKrylov(int argc,char **args);

// Solve for eigenvalues and eigenwvectors using WaveHoltz
int solveEigen(int argc,char **argv);

// Solve Helmholtz using WaveHoltz
int solveHelmholtz(int argc,char **argv);

// Solve Helmholtz using a direct or iterative solver.
int solveHelmholtzDirect( realCompositeGridFunction & u, realCompositeGridFunction & f );

// Solve the FPI equations using PETSc 
int solvePETSc(int argc,char **args);

// Solve for eigenmodes using SLEPc
int solveSLEPc(int argc,char **args);



// The database is used to hold parameters
mutable DataBase dbase;

CompositeGrid & cg;
GenericGraphicsInterface & gi;

};



#endif  
