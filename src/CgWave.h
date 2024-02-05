#ifndef CG_WAVE_H
#define CG_WAVE_H

// ====================== Composite Grid Wave Equation Solver Class =====================

#include "Overture.h"

#define KK_DEBUG
#include "DBase.hh"
using namespace DBase;

class LcbcData;

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
  fixedPointEigenSolver
};

CgWave(CompositeGrid & cg, GenericGraphicsInterface & gi);
~CgWave();

int adjustEigenWaveFrequency();

int adjustSolutionForSuperGrid( realCompositeGridFunction & q, Real superGridLayerValue = 0. );

// advance the solution 
int advance( int it );

int applyBoundaryConditions( realCompositeGridFunction & u,  realCompositeGridFunction & un,  real t, 
                             bool applyExplicitBoundaryConditions = false,
                             bool fillImplicitBoundaryConditions = false );

// Apply BC's to an eigenfunction
int applyEigenFunctionBoundaryConditions( realCompositeGridFunction & u );

int assignLCBC( realMappedGridFunction & u, Real t, Real dt, int grid );

//  Evaluate the WaveHoltz beta function (single frequency)
static Real betaWaveHoltz( Real lambda, Real omega, Real T );

int buildRunTimeDialog();

int buildSuperGrid( );

bool adjustBoundsForAbsorbingLayer( MappedGrid & mg, Index Iv[3], int extra =0 );

int checkDeflation();

int correctEigenfunction();

// deflate WaveHoltz solution
int deflateSolution();

void displayBoundaryConditions( FILE *file = stdout );

int formImplicitTimeSteppingMatrix();

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
int getFilterWeights( int Nt, int numFreq, const RealArray & Tv, int orderOfAccuracy, RealArray & sigma, RealArray & filterWeights );

// Fill RHS for direct Helmholtz solver
int getHelmholtzForcing( realCompositeGridFunction & f  );

// void getLcbcData(MappedGrid & mg, Real **&fn, Real **&gn, RealArray tmpGn[], RealArray tmpFn[], Real t );
void getLcbcData(MappedGrid & mg, LcbcData *&fn, LcbcData *&gn, RealArray tmpGn[], RealArray tmpFn[], Real t );

void getInitialConditions( int current, real t, bool getTimeDerivative = false );

int getIntegrationWeights( int Nt, int numFreq, const RealArray & Tv, int orderOfAccuracy, RealArray & sigma );

aString getMethodName() const;

int getMultiFrequencyWaveHoltzMatrix( RealArray & A );

Real getRayleighQuotient( realCompositeGridFunction & v, int component =0 );


int getTimeStep();
  
void getTimeSteppingLabel( real dt, aString & label ) const;

Real getUpwindDissipationCoefficient( int grid, Real dtUpwind = -1., bool adjustForTimeStep = true );

int getUserDefinedKnownSolution(real t,  int grid, realArray & ua, const Index & I1a, const Index &I2a, const Index &I3a, 
				int numberOfTimeDerivatives = 0 );

int getWaveHoltzIterationEigenvalue( RealArray & lambda, RealArray & mu );

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

// Assign parameters 
int interactiveUpdate();

int normalizeEigenvector( int eigNumber );

void outputHeader();

// Output results (e.g. errors to the checkFile)
int outputResults( int current, real t );

int printStatistics(FILE *file = stdout );

int plot( int current, real t, real dt );

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

// utility array
static int sortArray( RealArray & eigenValues, IntegerArray & iperm );

int takeFirstStep( int cur, real t );

int takeFirstStepHelmholtz( int cur, real t );

// Old way: 
int takeFirstBackwardStep( int cur, real t );

int takeImplicitStep( Real t );

int updateEigenmodes();

// update time-integral for Helmholtz projection
int updateTimeIntegral( int step, StepOptionEnum stepOption, real t, realCompositeGridFunction& u );

int updateUserDefinedKnownSolution();

int userDefinedForcing( realArray & f, int iparam[], real rparam[] );


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
    timeForTimeIntegral,
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
