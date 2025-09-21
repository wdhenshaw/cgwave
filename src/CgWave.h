#ifndef CG_WAVE_H
#define CG_WAVE_H

// ====================== Composite Grid Wave Equation Solver Class =====================

#include "Overture.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

class LcbcData;

// Function declaration for arguments to timeIntegralByQuadrature
typedef Real (*TimeFunc)( Real omega, Real lambda, Real t );

class CgWave
{

public:

enum BoundaryConditionEnum
{
  periodic       =-1,
  interpolation  = 0,
  dirichlet      = 1,
  neumann        = 2,
  evenSymmetry   = 3,
  radiation      = 4,
  exactBC        = 5,  // Set exact values on boundary and ghost
  abcEM2         = 6,  // absorbing BC, Engquist-Majda order 2  
  characteristic = 7,  // characteristic BC
  absorbing      = 8,   // for SuperGrid
};

enum InitialConditionOptionEnum
{
  zeroInitialCondition=0,
  twilightZoneInitialCondition,
  knownSolutionInitialCondition,
  pulseInitialCondition,
  randomInitialCondition
};

enum ForcingOptionEnum
{
  noForcing=0,
  twilightZoneForcing,
  userForcing,
  helmholtzForcing
};
    
enum TwilightZoneEnum
{
  polynomial,
  trigonometric
};


enum StepOptionEnum
{
  firstStep=0,
  middleStep,
  lastStep
};

enum PlottingOptionsEnum
{
  noPlotting =0,
  plotAndWait=1,
  plotNoWait =2
};

enum TimeSteppingMethodEnum
{
  explicitTimeStepping=0,
  implicitTimeStepping=1
};

enum BoundaryConditionApproachEnum
{
  defaultBoundaryConditionApproach=0,
  useOneSidedBoundaryConditions,
  useCompatibilityBoundaryConditions,
  useLocalCompatibilityBoundaryConditions
};

enum AssignInterpolationNeighboursEnum
{
  defaultAssignInterpNeighbours=0,
  extrapolateInterpNeighbours,
  interpolateInterpNeighbours
};

enum ModifiedEquationApproachEnum
{
  standardModifiedEquation=0,
  hierarchicalModifiedEquation,
  stencilModifiedEquation
};

enum EigenSolverEnum
{
  defaultEigenSolver=0,
  KrylovSchurEigenSolver,
  ArnoldiEigenSolver,
  ArpackEigenSolver,
  fixedPointEigenSolver,
  powerEigenSolver,
  inverseIterationEigenSolver,
  JacobiDavidsonEigenSolver,
  subspaceIterationEigenSolver
};

enum EigenSolverInitialConditionEnum
{
  defaultEigenSolverInitialCondition=0, // this means let the eigensolver choose the initial condition
  randomEigenSolverInitialCondition,
  sineEigenSolverInitialCondition,
  sumOfEigenvectorsInitialCondition
};

CgWave(CompositeGrid & cg, GenericGraphicsInterface & gi);
~CgWave();

int adjustEigenWaveFrequency();

int adjustSolutionForSuperGrid( realCompositeGridFunction & q, Real superGridLayerValue = 0. );

// advance the solution 
int advance( int it );

// advance the solution -- old version --
int advanceOld( int it );

int applyBoundaryConditions( realCompositeGridFunction & u,  realCompositeGridFunction & un,  real t, 
                             bool applyExplicitBoundaryConditions = false,
                             bool fillImplicitBoundaryConditions = false );

// Apply BC's to an eigenfunction
int applyEigenFunctionBoundaryConditions( realCompositeGridFunction & u );

int assignLCBC( realMappedGridFunction & u, Real t, Real dt, int grid );

//  Evaluate the WaveHoltz beta function (single frequency)
static Real betaWaveHoltz( Real lambda, Real omega, Real T );

// Evaluate the discrete WaveHoltz beta function (single frequency)
static Real betaDiscreteWaveHoltz( Real lambda, Real omega, Real T, Real dt );

int buildRunTimeDialog();

int buildSuperGrid( );

bool adjustBoundsForAbsorbingLayer( MappedGrid & mg, Index Iv[3], int extra =0 );

int checkDeflation();

// check parameters for consistency
int checkParameters(); 

int correctEigenfunction();

static Real coscd( Real x, Real T, Real dt );

static Real coscdPrime( Real x, Real T, Real dt );

Real cosFilter( int freq, Real lam );
Real cosFilterPrime( int freq, Real lam );


// deflate the WaveHoltz solution (or forcing)
int deflateSolution( int deflateOption= 0 );

void displayBoundaryConditions( FILE *file = stdout );

Real evalBetaFunction( const Real lam, const int freq, Real dt ) const;

int formImplicitTimeSteppingMatrix();

int getAdjustedBoundaryIndex( MappedGrid & mg, int side, int axis, Index & Ib1, Index & Ib2, Index & Ib3 );

realCompositeGridFunction& getAugmentedSolution(int current, realCompositeGridFunction & v, const real t);

Real getAverageNumberOfIterationsPerImplicitSolve() const;

realCompositeGridFunction& getCurrentSolution();
  
// Compute the residual in the eigenvalue equation:  || L v + lambda^2 v ||/lambda^2 
Real getEigenPairResidual( Real lambda, realCompositeGridFunction & v, realCompositeGridFunction & res, int component /* =0 */ );

Real getEnergyNorm( int cur, Real t );

real getErrors( realCompositeGridFunction & u, real t );

// compute errors in the eigenvalues/eigenvectors (when the true values are known)
int getErrorsInEigenmodes( Real & relErrEigenvalue, Real & relErrEigenvector, bool checkRayleighRitz=false );

// compute errors in a single eigenvalue, eigenvector pair (when the true values are known)
int getErrorInEigenPair( const Real lambda, realCompositeGridFunction & eigVector, const int component, 
                         Real & lambdaTrue, Real & relErrEigenvalue, Real & relErrEigenvector, 
                         int & eigIndex, int & multipleEigIndex, bool saveErrors = false );

// Get weights for WaveHoltz filters
int getFilterWeights( int Nt, Real dt, int numFreq, const RealArray & Tv, int orderOfAccuracy, RealArray & sigma, RealArray & filterWeights );

// Fill RHS for direct Helmholtz solver
int getHelmholtzForcing( realCompositeGridFunction & f  );

// void getLcbcData(MappedGrid & mg, Real **&fn, Real **&gn, RealArray tmpGn[], RealArray tmpFn[], Real t );
void getLcbcData(MappedGrid & mg, LcbcData *&fn, LcbcData *&gn, RealArray tmpGn[], RealArray tmpFn[], Real t );

void getInitialConditions( int current, real t, bool getTimeDerivative = false );

int getIntegrationWeights( int Nt, int numFreq, const RealArray & Tv, int orderOfAccuracy, RealArray & sigma );

aString getMethodName() const;

int getMultiFrequencyWaveHoltzMatrix( RealArray & A, bool useAdjusted = true  );

Real getRayleighQuotient( realCompositeGridFunction & v, int component =0 );


int getTimeStep();
  
void getTimeSteppingLabel( real dt, aString & label ) const;

Real getUpwindDissipationCoefficient( int grid, Real dtUpwind = -1., bool adjustForTimeStep = true );

int getUserDefinedKnownSolution(real t,  int grid, realArray & ua, const Index & I1a, const Index &I2a, const Index &I3a, 
				int numberOfTimeDerivatives = 0 );

int getWaveHoltzIterationEigenvalue( RealArray & lambda, RealArray & mu, bool useAdjusted = true  );

// Check if two composite grids match
static bool compositeGridsMatch( CompositeGrid & cg, CompositeGrid & cgsf );

Real idFilter( int freq, Real lam );
Real idFilterPrime( int freq, Real lam );

// inflate WaveHoltz solution
int inflateSolution();

// Initialize time-step and forcing 
int initialize();

// Initialize deflation for WaveHoltz
int initializeDeflation();


// Initialize local compatbility boundart conditions
int initializeLCBC();

// Initialize PETSc as needed
int initializePETSc( int argc = 0 , char **args = NULL );

// Initialize the time integral used by WaveHoltz
int initializeTimeIntegral( Real dt );

// Assign parameters 
int interactiveUpdate();

int normalizeEigenvector( int eigNumber );

int optFilterParameters( const RealArray & frequencyArray, const RealArray & periodArray, const Real dt, int & Nlam, Real & muMin, Real & muMax);
void outputHeader();

// Output results (e.g. errors to the checkFile)
int outputResults( int current, real t );

int printStatistics(FILE *file = stdout );

int plot( int current, real t, real dt );

// plot the WaveHoltz filter beta and any eigenvalues 
int plotFilter( RealArray & eigenValues );

int adjustFrequencyArrays();

int resetFrequencyArrays();

// Reset CPU timings to zero:
int resetTimings();

// Re-initialize deflation (if the number of vectors to deflate has changed, for example)
int reinitializeDeflation();

int saveSequenceInfo( real t0, RealArray & sequenceData );

int saveSequencesToShowFile();

int saveShow( int current, Real t, Real dt );

int setNameOfGridFile( aString & nameOfOGFile );

// Setup grids and grid functions
int setup();

int setupUserDefinedForcing();

static Real sincd( Real x, Real T, Real dt );

static Real sincdPrime( Real x, Real T, Real dt );

Real sinFilter( int freq, Real lam );
Real sinFilterPrime( int freq, Real lam );

// utility array
static int sortArray( RealArray & eigenValues, IntegerArray & iperm );

int takeFirstStep( int cur, real t );

int takeFirstStepHelmholtz( int cur, real t );

// Old way: 
int takeFirstBackwardStep( int cur, real t );

int takeImplicitStep( Real t );

Real timeIntegralByQuadrature( TimeFunc g, Real lambda, int freq );

int updateEigenmodes();

// update time-integral for Helmholtz projection
int updateTimeIntegral( int step, StepOptionEnum stepOption, Real t, Real dt, realCompositeGridFunction& u );

int updateUserDefinedKnownSolution();

int userDefinedForcing( realArray & f, int iparam[], real rparam[] );

// ---- functions for matrix free routines ----
int getActivePointIndex( MappedGrid & mg, Index *Iv );
int initializeGlobalIndexing( bool checkMask = true );
int gridFunctionToVector( const realCompositeGridFunction & u, Real *v, int iStart, int iEnd );
int vectorToGridFunction( const Real *v, realCompositeGridFunction & u, int iStart, int iEnd );

  enum TimingEnum
  { 
    totalTime=0,
    timeForInitialize,
    timeForInitializeBCs,
    timeForInitialConditions,
    timeForAdvance,
    timeForAdvanceRectangularGrids,
    timeForAdvanceCurvilinearGrids,
    timeForImplicitSolve,
    timeForDissipation,
    timeForBoundaryConditions,
    timeForInterpolate,
    timeForUpdateGhostBoundaries,
    timeForForcing,
    timeForUserDefinedKnownSolution,
    timeForTimeIntegral,
    timeForDeflation,
    timeForGetError,
    timeForPlotting,
    timeForOutputResults,
    timeForWaiting,
    maximumNumberOfTimings      // number of entries in this list
  };
  RealArray timing;                                     // for timings, cpu time for some function
  aString timingName[maximumNumberOfTimings];           // name of the function being timed

// The database is used to hold parameters
mutable DataBase dbase;

CompositeGrid & cg;
GenericGraphicsInterface & gi;


protected:

  int buildBoundaryConditionOptionsDialog(DialogData & dialog );
  int buildEigenWaveOptionsDialog(DialogData & dialog );
  int buildWaveHoltzOptionsDialog(DialogData & dialog );

  int getBoundaryConditionOption(const aString & answer, DialogData & dialog,IntegerArray & bcOrig );
  int getEigenWaveOption(const aString & answer, DialogData & dialog );
  int getWaveHoltzOption(const aString & answer, DialogData & dialog );


};



#endif  
