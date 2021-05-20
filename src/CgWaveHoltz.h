#ifndef CG_WAVE_HOLTZ_H
#define CG_WAVE_HOLTZ_H

#include "mpi.h"
#include "Overture.h"
#include "Oges.h"


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

class CgWaveHoltz
{

public:

CgWaveHoltz(CompositeGrid & cg, GenericGraphicsInterface & gi);
~CgWaveHoltz();

// advance the solution over a period (or multiple periods)
// int advance( int it );

// int getTimeStep();
  
// Initialize time-step and forcing 
int initialize();

// Assign parameters 
int interactiveUpdate();

int outputHeader();

// save results to a Matlab file
int outputMatlabFile();

// Compute the residual in the current solution
real residual();

int setNameOfGridFile( aString & nameOfOGFile );

// Setup grids and grid functions
int setup();

// WaveHoltz iteration:
int solve();

// Solve Helmholtz using a direct or iterative solver.
int solveHelmholtz( realCompositeGridFunction & u, realCompositeGridFunction & f );

// PETSc solver 
int solvePETSc(int argc,char **args);

// The database is used to hold parameters
mutable DataBase dbase;

CompositeGrid & cg;
GenericGraphicsInterface & gi;

};



#endif  
