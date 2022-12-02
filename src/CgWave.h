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
  periodic      =-1,
  interpolation = 0,
  dirichlet     = 1,
  neumann       = 2,
  evenSymmetry  = 3,
  radiation     = 4,
  exactBC       = 5   // Set exact values on boundary and ghost
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

CgWave(CompositeGrid & cg, GenericGraphicsInterface & gi);
~CgWave();

// advance the solution 
int advance( int it );

int applyBoundaryConditions( realCompositeGridFunction & u, real t );

int assignLCBC( realMappedGridFunction & u, Real t, Real dt, int grid );

int buildRunTimeDialog();

void displayBoundaryConditions( FILE *file = stdout );

int formImplicitTimeSteppingMatrix();

realCompositeGridFunction& getAugmentedSolution(int current, realCompositeGridFunction & v, const real t);

realCompositeGridFunction& getCurrentSolution();
  
real getErrors( realCompositeGridFunction & u, real t );

// Fill RHS for direct Helmholtz solver
int getHelmholtzForcing( realCompositeGridFunction & f  );

// void getLcbcData(MappedGrid & mg, Real **&fn, Real **&gn, RealArray tmpGn[], RealArray tmpFn[], Real t );
void getLcbcData(MappedGrid & mg, LcbcData *&fn, LcbcData *&gn, RealArray tmpGn[], RealArray tmpFn[], Real t );

void getInitialConditions( int current, real t, bool getTimeDerivative = false );

int getIntegrationWeights( int Nt, int numFreq, const RealArray & Tv, int orderOfAccuracy, RealArray & sigma );

aString getMethodName() const;

int getMultiFrequencyWaveHoltzMatrix( RealArray & A );

int getTimeStep();
  
void getTimeSteppingLabel( real dt, aString & label ) const;

int getUserDefinedKnownSolution(real t,  int grid, realArray & ua, const Index & I1a, const Index &I2a, const Index &I3a, 
				int numberOfTimeDerivatives = 0 );

// Initialize time-step and forcing 
int initialize();

// Initialize local compatbility boundart conditions
int initializeLCBC();

// Assign parameters 
int interactiveUpdate();

void outputHeader();

// Output results (e.g. errors to the checkFile)
int outputResults( int current, real t );

int printStatistics(FILE *file = stdout );

int plot( int current, real t, real dt );

int saveSequenceInfo( real t0, RealArray & sequenceData );

int saveSequencesToShowFile();

int saveShow( int current, Real t, Real dt );

int setNameOfGridFile( aString & nameOfOGFile );

// Setup grids and grid functions
int setup();

int setupUserDefinedForcing();

int takeFirstStep( int cur, real t );

int takeFirstStepHelmholtz( int cur, real t );

// Old way: 
int takeFirstBackwardStep( int cur, real t );

int takeImplicitStep( Real t );

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

};



#endif  
