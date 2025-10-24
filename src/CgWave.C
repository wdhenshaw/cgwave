// ======================================================================================
// ====================== Composite Grid Wave Equation Solver Class =====================
// ======================================================================================

#include "CgWave.h"
#include "CompositeGridOperators.h";    
#include "PlotStuff.h"
#include "display.h"
#include "ParallelOverlappingGridInterpolator.h"
#include "ParallelUtility.h"
#include "LoadBalancer.h"
#include "gridFunctionNorms.h"
#include "OGPolyFunction.h"
#include "OGTrigFunction.h"
#include "DialogData.h"
#include "Ogshow.h"
#include "LCBC.h"
#include "OgesParameters.h"

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  

#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) \
                                 for( int side=0; side<=1; side++ )


// ================================================================================================
/// \brief Constructor for the CgWave class
// ================================================================================================
CgWave::
CgWave( CompositeGrid & cgIn, GenericGraphicsInterface & giIn ) : cg(cgIn), gi(giIn) 
{
  int & np = dbase.put<int>("np"); // number of processors 
  np = max(1,Communication_Manager::numberOfProcessors());

  int & myid = dbase.put<int>("myid"); // my processor number
  myid = max(0,Communication_Manager::My_Process_Number);

  FILE *& debugFile  = dbase.put<FILE*>("debugFile");
  FILE *& pDebugFile = dbase.put<FILE*>("pDebugFile");

  // *** open multiple debug files for each processor ****
  aString buff;
  #ifndef USE_PPP
    debugFile   = fopen("cgWave.debug","w" );        // Here is the log file
    pDebugFile= debugFile;
  #else
    debugFile = fopen(sPrintF(buff,"cgWaveNP%i.debug",np),"w" );  // Here is the debug file
    pDebugFile = fopen(sPrintF(buff,"cgWaveNP%ip%i.debug",np,myid),"w");
  #endif

  FILE *& logFile = dbase.put<FILE*>("logFile");
  logFile = fopen("cgWave.out","w" );        // log file 

  FILE *& checkFile = dbase.put<FILE*>("checkFile");
  checkFile = fopen("cgWave.check","w" );        // for regression and convergence tests
  fPrintF(checkFile,"# Check file for CgWave\n"); // check file has one title line

  dbase.put<Real>("c")=1.;
  dbase.put<Real>("cfl")=.7; // .25;
  dbase.put<Real>("tFinal")=1.;
  dbase.put<Real>("tPlot")=.1;
  dbase.put<Real>("dt")=-1.;
  dbase.put<Real>("dtUsed")=-1.; // dt actually used
  dbase.put<Real>("dtMax")=-1.;  // save maximum allowable dt (before corrections to reach a given time)
  // domainSize : estimate of the domain size (used, for example, to compute estimate points-per-wavelength for Helmholtz problems)
  //            : -1 : compute from grids
  dbase.put<Real>("domainSize")=-1.;  

  dbase.put<int>("upwind")=0;                // use upwind dissipation
  dbase.put<int>("numUpwindCorrections")=1;  // number of upwind corrections

  dbase.put<Real>("ad4")=0.;   // coeff of the artificial dissipation. (*old)
  dbase.put<int>("dissipationFrequency")=1; // apply dissipation every this many steps (1= every step)
  // preComputeUpwindUt : true=precompute Ut in upwind dissipation, 
  //                      false=compute Ut inline in Gauss-Seidel fashion  
  dbase.put<int>("preComputeUpwindUt")= false;  
  dbase.put<int>("implicitUpwind") = false; // if true, include upwinding in implicit matrix when implicit time-stepping

  // For upwind schemes with wider stencils we can extrap or interp unused points next to interpolation points 
  dbase.put<AssignInterpolationNeighboursEnum>("assignInterpNeighbours")=defaultAssignInterpNeighbours;

  dbase.put<int>("computeErrors") = 1;                  // by default, compute errors for TZ or a known solution
  dbase.put<int>("computeEnergy") = 0;                  // 1= compute energy
  dbase.put<int>("saveMaxErrors") = 0;                  // save max errors over time as a sequence in the show file
  dbase.put<int>("applyKnownSolutionAtBoundaries") = 0; // by default, do NOT apply known solution at boundaries
  dbase.put<int>("bcCount") = 0;                        // count the number of times applyBC is called

  dbase.put<int>("useKnownSolutionForFirstStep") = 0;   // use known solution for first time-step

  dbase.put<int>("current")=0;              // holds current solution

  dbase.put<int>("orderOfAccuracy")=2;
  dbase.put<int>("orderOfAccuracyInTime")=-1; // -1 means make order in time equal to order in space

  dbase.put<int>("orderOfExtrapolation")=-1; // for boundary conditions, -1= default = orderOfAccuracy+1

  dbase.put<int>("secondOrderGrid")=1;
  dbase.put<int>("debug")=0;

  dbase.put<int>("initTimeIntegral")       = 1;           // initialize time integral for WaveHoltz 
  dbase.put<int>("useFilterWeights")       = 1;           //  1 = new filter weights for WaveHoltz
  dbase.put<int>("filterTimeDerivative")   = 0;           //  filter time-derivative for WaveHoltz
  dbase.put<int>("numCompWaveHoltz")       = 0;           //  number of components in WaveHoltz
  dbase.put<int>("filterD0t")              = 0;           // directly filter D0t instead of using integration by parts
  dbase.put<int>("useOptFilter")           = 0;           // use optimized filter 

  dbase.put<int>("deflateWaveHoltz") = 0;               //  set to 1 to turn on deflation for WaveHoltz
  dbase.put<int>("numToDeflate")     = 1;               //  number of eigenvectors to deflate
  dbase.put<int>("deflateForcing")   = 0;               //  1 = "deflate" forcing at start but do not adjust iterates 
  dbase.put<bool>("deflationInitialized") = false;      //  set too true if deflation has been initialized
  dbase.put<aString>("eigenVectorFile")= "none";        //  name of file holding eigs and eigenvectors for deflation
  dbase.put<Real>("eigenVectorForcingFactor")=0.;       // scale eigenvector forcing by this amount
  dbase.put<int>("numEigsToCompute")=1;                 // number of eigenpairs to compute 
  dbase.put<int>("numArnoldiVectors")=-1;               // total number of Arnolidi vectors to keep (-1 = use default)
  dbase.put<int>("useAugmentedGmres")=0;                // 1 = use augmented GMRES for deflation
  dbase.put<int>("augmentedVectorsAreEigenvectors")=1;  // 1 = augmented vectors are true discrete eigenvectors
  dbase.put<int>("eigenVectorsAreOrthogonal")=false;    // 1 = eigenvectors read from a show file have are orthogonal
  dbase.put<int>("onlyLoadDeflatedEigenVectors")=false;    // 1 = only read in eigenvectors used for deflation from the show file (*new* way)

  dbase.put<int>("plotFilterAndDeflatedEigenvalues")=false;    // 1 = plot filter function and deflated eigenvalues after they are found

  dbase.put<int>("eigenvectorsAreTrueEigenvectors")=true;  // eigenvectors for deflation may be the true ones or approximate
  dbase.put<RealArray>("betaDeflate");                     // beta values for deflated eigenvectors 

  dbase.put<int>("iterationStartRR") = 5;
  dbase.put<int>("useAccurateInnerProduct")=true; // use accurate inner product for computing the Rayleigh quotient
  dbase.put<Real>("eigenValueTolForMultiplicity")=1.e-4;  // tolerance for eignavlue multiplicities


  dbase.put<TimeSteppingMethodEnum>("timeSteppingMethod")=explicitTimeStepping;
  dbase.put<int>("takeImplicitFirstStep")=0;     // 1 = For implicit time-stepping, use an implicit first step


  // We have different versions of modified equation time-stepping
  //   0 = standard
  //   1 = hierachical (faster but more memory)
  //   2 = stencil
  dbase.put<ModifiedEquationApproachEnum>("modifiedEquationApproach")=standardModifiedEquation;

  // coefficients in implicit time-stepping  *check me*
  //  D+t D-t u =              c^2 Delta( cImp(1,0) *u^{n+1} + cImp(0,0) *u^n + cImp(-1,0)* u^{n-1} )   :  second-order coeff cImp(-1:1,0)
  //                 -(c^4*dt^2) Delta^2( cImp(1,1) *u^{n+1} + cImp(0,1) *u^n + cImp(-1,1)* u^{n-1}  )  :  fourth-order ceoff cImp(-1:1,1) 
  RealArray & bImp = dbase.put<RealArray>("bImp");
  RealArray & cImp = dbase.put<RealArray>("cImp");
  const int maxOrderOfAccuracy=12; // being rather hopeful here 
  bImp.redim(2,maxOrderOfAccuracy);             bImp=0.; 
  cImp.redim(Range(-1,1),maxOrderOfAccuracy);   cImp=0.; 
  // For accuracy the weights depend on one parameter beta2 for second-order,
  // and a second parameter beta4 for fourth-order: (See notes in research/timeStepping/implicitTaylorSchemes.pdf)
  // Full-weighting for second-order part by default: 
  Real beta2=.5, beta4=0.; 
  bImp(0)=beta2;
  bImp(1)=beta4; 
  Real alpha2 = (1.-beta2)/2.;
  Real alpha4 = (alpha2-beta4-1./12.)/2.; 
  cImp(-1,0)=alpha2;
  cImp( 0,0)= beta2;
  cImp( 1,0)=alpha2;
  cImp(-1,1)=alpha4;
  cImp( 0,1)= beta4;
  cImp( 1,1)=alpha4;  

  // Old way: By default choose the implicit dt based on the CFL: (otherwise us dtMax)
  // dbase.put<int>("TImplicitTimeStepFromCFL")=1;
  // New way: 
  //   chooseTimeStepFromExplicitGrids = 1 : choose dt from CFL and grid spacing on explicit grids only, or if all grids are implicit
  //                                   = 0 : choose dt from CFL and grid spacing on all grids
  dbase.put<int>("chooseTimeStepFromExplicitGrids")=1;

  dbase.put<bool>("hasRadiationBoundaryConditions")=false; // true if there are radition BCs

  dbase.put<int>("useSuperGrid")=0;      
  dbase.put<Real>("superGridWidth")=.2;          // superGrid layer width (in parameter space)
  dbase.put<int>("initializeSuperGrid")=1;
  dbase.put<IntegerArray>("superGrid");          // superGrid(grid) = 1 if this grid uses superGrid
  dbase.put<int>("adjustPlotsForSuperGrid")=1;    // set solution to zero in any superGridLayers
  dbase.put<int>("adjustErrorsForSuperGrid")=1;  // set errors and residuals to zero in any superGridLayers

  dbase.put<IntegerArray>("useAbsorbingLayer");  //  useAbsorbingLayer(axis,grid) 

  // interactiveMode :
  //       0 = plot intermediate results
  //       1 = run without plotting and exit advance when finished 
  dbase.put<int>("interactiveMode")=0;
  dbase.put<int>("plotFrequency")= INT_MAX; // another way to turn on plotting every this many steps

  dbase.put<int>("solveForScatteredField")= 0;  // 1=we are solving for the scattered field
  dbase.put<int>("plotScatteredField")= 0;  // 1=plot scattered field for scattering problems
  // Parameters for the plane wave defining the incident field
  //   sin( kx*x + ky*y +kz*z - omega*t + phi )
  dbase.put<Real>("ampPlaneWave")   = 1.;
  dbase.put<Real>("kxPlaneWave")    = twoPi;
  dbase.put<Real>("kyPlaneWave")    = 0.;
  dbase.put<Real>("kzPlaneWave")    = 0.;
  dbase.put<Real>("phiPlaneWave")   = 0.;       // phase 
  dbase.put<Real>("omegaPlaneWave") = twoPi*1;  // c*k 

  dbase.put<Real>("numberOfGridPoints")=0.;
  dbase.put<int>("numberOfStepsTaken")=0;            // total steps taken
  dbase.put<int>("numberOfStepsPerSolve")=0;         // number of steps take per solve
  dbase.put<int>("totalImplicitIterations")=0;       // total iterations used in implicit solves
  dbase.put<int>("totalImplicitSolves")=0;           // total number of implicit solves
  
  dbase.put<RealArray>("dxMinMax");
  dbase.put<aString>("nameOfGridFile")="unknown";

  dbase.put<Real>("maxError")=0.;      // save max-error here 
  dbase.put<Real>("solutionNorm")=1.;  // save solution norm here 
  // dbase.put<int>("computeErrors")=0;   // true of we compute errors 

  // For Helmholtz solve with CgWaveHoltz
  dbase.put<int>("solveHelmholtz")=false;                           // if true,  use the WaveHoltz algorithm to solve the Helmholtz equation
  dbase.put<int>("computeEigenmodes")=false;                        // if true, use the WaveHoltz algorithm to compute eigenvalues and eigenvectors
  dbase.put<int>("computeEigenmodesWithSLEPc")=false;               // if true, we are solving for eigenvalues and eigenvectors with SLEPc
  dbase.put<EigenSolverEnum>("eigenSolver") = defaultEigenSolver;   // eigensolver used by SLEPSc
  dbase.put<int>("initialVectorsForEigenSolver")=true;              // provide initial vectors used by SLEPSc solvers
  dbase.put<int>("initialVectorSmooths")=500;                       // number of times to smooth the initial vectors for the eigen solvers

  dbase.put<EigenSolverInitialConditionEnum>("eigenSolverInitialCondition") = defaultEigenSolverInitialCondition; 

  dbase.put<int>("computeTimeIntegral")=false; // if true compute the time-integral (for the Helmholtz solve or other reason)
  dbase.put<int>("iteration")=-1;              // iteration number of WaveHoltz
  // Save WaveHoltz "residuals" by iteration: 
  // resVector(it) = norm( v^{n+1} - v^n )
  dbase.put<RealArray>("resVector");

  dbase.put<int>("adjustOmega")=0;                    // 1 : choose omega from the symbol of D+t D-t 
  dbase.put<int>("adjustHelmholtzForUpwinding")=1;    // 1 : correct Helmholtz for upwinding (when adjustOmega=1)
  real & omega = dbase.put<Real>("omega")=30.1;
  dbase.put<Real>("Tperiod")=twoPi/omega;
  dbase.put<int>("numPeriods")=10;
  dbase.put<int>("minStepsPerPeriod")= 10; // take at least this many steps for implicit time-stepping
  dbase.put<int>("numberOfRitzVectors")=10; // maximum number of Rayleigh-Ritz vectors to use 
  dbase.put<int>("assignRitzFrequency")=5;  // set solution to latest Ritz vector every this many steps

  dbase.put<Real>("waveHoltzAsymptoticConvergenceRate")=-1.; // in some cases we can compute this 

  dbase.put<Real>("omegaSave")   = -1.;  // used when we adjust omega 
  dbase.put<Real>("TperiodSave") = -1.;
  dbase.put<Real>("dtSave")      = -1.;

  // ---- For multi-frequency WaveHoltz ----
  int & numberOfFrequencies          = dbase.put<int>("numberOfFrequencies")=1;
  RealArray & frequencyArray         = dbase.put<RealArray>("frequencyArray");
  RealArray & frequencyArrayAdjusted = dbase.put<RealArray>("frequencyArrayAdjusted");

  RealArray & periodArray            = dbase.put<RealArray>("periodArray");
  RealArray & periodArrayAdjusted    = dbase.put<RealArray>("periodArrayAdjusted");

  IntegerArray & numPeriodsArray     = dbase.put<IntegerArray>("numPeriodsArray");
  frequencyArray.redim(numberOfFrequencies);
  frequencyArray(0)=30.1;
  periodArray.redim(numberOfFrequencies);
  periodArray(0) = twoPi/frequencyArray(0);
  numPeriodsArray.redim(numberOfFrequencies);
  numPeriodsArray=1; 

  RealArray & frequencyArraySave = dbase.put<RealArray>("frequencyArraySave");
  frequencyArraySave.redim(numberOfFrequencies); frequencyArraySave=0.; 

  frequencyArrayAdjusted.redim(numberOfFrequencies); frequencyArrayAdjusted=0.; 

  RealArray & periodArraySave    = dbase.put<RealArray>("periodArraySave");
  periodArraySave.redim(numberOfFrequencies);   periodArraySave=0.;

  dbase.put<Real>("damp")                = 0.;          //  coefficient of linear damping
  dbase.put<Real>("dampSave")            = 0.;          //  saved value when adjusting
  
  Real & omegaSign = dbase.put<Real>("omegaSign")= -1.;         //  time dependence is exp( omegaSign*I*omega*t ) for complex Helmholtz problems
  dbase.put<Real>("viFactor")                    = -omegaSign;  //  factor relating vk(:,:,2) to ui in u = ur*cos(omega*t) + ui*sin(omega*t)
  
  dbase.put<Real>("tol")=1.e-4;  // tolerance for Krylov solvers 

  real & omegaSOR = dbase.put<Real>("omegaSOR")=1.;

  // // Gaussian forcing: 
  // dbase.put<Real>("beta")=100.;
  // dbase.put<Real>("x0")=0.;
  // dbase.put<Real>("y0")=0.;
  // dbase.put<Real>("z0")=0.;

  int & numberOfTimeLevelsStored = dbase.put<int>("numberOfTimeLevelsStored")=3;

  realCompositeGridFunction *& ucg = dbase.put<realCompositeGridFunction*>("ucg");
  ucg = new realCompositeGridFunction[numberOfTimeLevelsStored];
  for( int n=0; n<numberOfTimeLevelsStored; n++ )
  {
    ucg[n].setName("u");
    ucg[n].setName("u",0);
  }
 
  dbase.put<realCompositeGridFunction>("f");  // source term 

  dbase.put<CompositeGridOperators>("operators");
  dbase.put<Interpolant*>("pInterpolant")=NULL;

  // For the Helmholtz solution: 
  dbase.put<realCompositeGridFunction>("v");
  dbase.put<realCompositeGridFunction>("vOld");

  dbase.put<IntegerArray>("gridIsImplicit"); 

  // dbase.put<realCompositeGridFunction>("vOld");
  // dbase.put<realCompositeGridFunction>("residual");

  dbase.put<int>("petscIsInitialized")=false;

  dbase.put<OGFunction*>("tz")=NULL;
  dbase.put<TwilightZoneEnum>("twilightZone")=polynomial;
  dbase.put<int>("degreeInSpace")=2;
  dbase.put<int>("degreeInTime")=2;

  // frequencies for trigonometric TZ
  RealArray & trigFreq = dbase.put<RealArray>("trigFreq");
  trigFreq.redim(4); 
  trigFreq=2.; 

  dbase.put<int>("addForcing")=0;
  dbase.put<ForcingOptionEnum>("forcingOption")=noForcing;

  dbase.put<InitialConditionOptionEnum>("initialConditionOption")=zeroInitialCondition;

  dbase.put<BoundaryConditionApproachEnum>("bcApproach")=defaultBoundaryConditionApproach;


  dbase.put<realCompositeGridFunction>("error");

  dbase.put<GUIState*>("runTimeDialog")=NULL;

  dbase.put<int>("movieFrame")=-1;
  dbase.put<int>("plotOptions")=1;
  dbase.put<int>("plotChoices")=0;
  // dbase.put<int>("plotHelmholtzErrors")=0;  // 1 = compute error w.r.t the Helmholtz solution

  dbase.put<int>("checkParameters")=1;  // 1= check problem parameters for consistency
  
  dbase.put<PlotStuffParameters>("psp");

  dbase.put<aString>("movieFileName")="cgWave";

  dbase.put<aString>("knownSolutionOption")="noKnownSolution";
  dbase.put<bool>("knownSolutionIsTimeDependent")=false; 

  // These should match BoundaryConditionEnum
  int & numberOfBCNames = dbase.put<int>("numberOfBCNames")=9;
  aString *& bcNames = dbase.put<aString*>("bcNames") = new aString [numberOfBCNames];
  bcNames[0]="periodic";
  bcNames[1]="dirichlet";
  bcNames[2]="neumann";
  bcNames[3]="evenSymmetry";
  bcNames[4]="radiation";
  bcNames[5]="exact";
  bcNames[6]="abcEM2";
  bcNames[7]="characteristic";
  bcNames[8]="absorbing";


  // show file variables:
  dbase.put<bool>("saveShowFile")     =false; 
  dbase.put<aString>("nameOfShowFile")="cgWave.show"; // name of the show file
  dbase.put<Ogshow*>("showFile")      =NULL;
  dbase.put<int>("flushFrequency")    =100;  // number of solutions per sub-showFile
  dbase.put<int>("numberOfSequences") =1;   // for sequences in the show file
  dbase.put<int>("sequenceCount")     =0;   // for sequences in the show file
  dbase.put<RealArray>("timeSequence");
  dbase.put<RealArray>("sequence");

  dbase.put<Real>("timeForOgesSolve") = 0;  // record Oges time for implicit solve

  // enum TimingEnum
  // { 
  //   totalTime=0,
  //   timeForInitialize,
  //   timeForInitializeBCs,
  //   timeForInitialConditions,
  //   timeForAdvance,
  //   timeForAdvanceRectangularGrids,
  //   timeForAdvanceCurvilinearGrids,
  //   timeForImplicitSolve,
  //   timeForDissipation,
  //   timeForBoundaryConditions,
  //   timeForInterpolate,
  //   timeForDeflation,
  //   timeForUpdateGhostBoundaries,
  //   timeForForcing,
  //   timeForUserDefinedKnownSolution,
  //   timeForTimeIntegral,
  //   timeForGetError,
  //   timeForPlotting,
  //   timeForOutputResults,
  //   timeForWaiting,
  //   maximumNumberOfTimings      // number of entries in this list
  // };


  timing.redim(maximumNumberOfTimings);
  timing=0.;
  for( int i=0; i<maximumNumberOfTimings; i++ )
    timingName[i]="";

  // only name the things that will be really timed in this run
  timingName[totalTime]                          ="total time";
  timingName[timeForInitialize]                  ="setup and initialize";
  timingName[timeForInitializeBCs]               ="initialize BCs";
  timingName[timeForInitialConditions]           ="initial conditions";
  timingName[timeForAdvance]                     ="advance";
  timingName[timeForAdvanceRectangularGrids]     ="  advance rectangular grids";
  timingName[timeForAdvanceCurvilinearGrids]     ="  advance curvilinear grids";
  timingName[timeForImplicitSolve]               ="    implicit solve";
  timingName[timeForForcing]                     ="  add forcing";
  timingName[timeForDissipation]                 ="  add dissipation";
  timingName[timeForBoundaryConditions]          ="  boundary conditions";
  timingName[timeForUserDefinedKnownSolution]    ="    user known solution";
  timingName[timeForUpdateGhostBoundaries]       ="  update ghost (parallel)";
  timingName[timeForInterpolate]                 ="  interpolation";
  timingName[timeForTimeIntegral]                ="  time integral";
  timingName[timeForDeflation]                   ="  deflation";
  timingName[timeForGetError]                    ="  get errors";
  timingName[timeForPlotting]                    ="  plotting";
  timingName[timeForOutputResults]               ="output results";
  timingName[timeForWaiting]                     ="waiting (not counted)";



}

// ================================================================================================
/// \brief Destructor for the CgWave class
// ================================================================================================
CgWave::
~CgWave()
{
  fclose(dbase.get<FILE*>("debugFile"));
  fclose(dbase.get<FILE*>("logFile"));
  fclose(dbase.get<FILE*>("checkFile"));

  delete dbase.get<OGFunction*>("tz"); 

  realCompositeGridFunction *& ucg = dbase.get<realCompositeGridFunction*>("ucg");
  delete [] ucg;

  Interpolant *& pInterpolant = dbase.get<Interpolant*>("pInterpolant");
  if( pInterpolant!=NULL && pInterpolant->decrementReferenceCount()==0 )
    delete pInterpolant;

  delete [] dbase.get<aString*>("bcNames");

  delete dbase.get<GUIState*>("runTimeDialog");

  Ogshow *& showFile = dbase.get<Ogshow*>("showFile");
  delete showFile;  

  if( dbase.has_key("LCBC") )
  {
    delete [] dbase.get<Lcbc*>("LCBC");
  }

  if( dbase.has_key("lapCoeff") )
  { // coefficients in the Laplacian for HA scheme
    delete [] dbase.get<RealArray*>("lapCoeff");
  }

  if( dbase.has_key("stencilCoeff") )
  { // stencil coefficients
    delete [] dbase.get<RealArray*>("stencilCoeff");
  }  

  if( dbase.has_key("rxOriginal") )
  {
    delete [] dbase.get<RealArray*>("rxOriginal");
  }

  if( dbase.has_key("indexVector") )
  {
    delete [] dbase.get<IntegerArray*>("indexVector");
  }

}

// ================================================================================================
/// \brief Adjust the frequency arrays for WaveHoltz 
/// \note this was out here for plotFilter
// ================================================================================================
int CgWave::adjustFrequencyArrays()
{
  Real & tFinal                     = dbase.get<real>("tFinal");
  Real & omega                      = dbase.get<real>("omega");
  Real & damp                       = dbase.get<Real>("damp"); 
  Real & dampSave                   = dbase.get<Real>("dampSave"); 
  Real & Tperiod                    = dbase.get<real>("Tperiod");
  int & numPeriods                  = dbase.get<int>("numPeriods");
  Real & omegaSave                  = dbase.get<real>("omegaSave");
  Real & TperiodSave                = dbase.get<real>("TperiodSave");  
  Real & dt                         = dbase.get<real>("dt");
  Real & dtSave                     = dbase.get<real>("dtSave");

  RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
  RealArray & frequencyArrayAdjusted= dbase.get<RealArray>("frequencyArrayAdjusted");
  RealArray & periodArray           = dbase.get<RealArray>("periodArray"); 
  RealArray & periodArrayAdjusted   = dbase.get<RealArray>("periodArrayAdjusted"); 
  IntegerArray & numPeriodsArray    = dbase.get<IntegerArray>("numPeriodsArray");
  RealArray & frequencyArraySave    = dbase.get<RealArray>("frequencyArraySave");
  RealArray & periodArraySave       = dbase.get<RealArray>("periodArraySave");     

  OV_ABORT("finish me");
  // omegaSave   = omega;
  // TperiodSave = Tperiod;
  // dtSave      = dt; 
  // Tperiod     = Tperiod; 
  // dampSave    = damp;

  // frequencyArraySave = frequencyArray;  
  // periodArraySave    = periodArray;     


  return 0;
}

// ================================================================================================
/// \brief Reset frequency arrays to their un-adjusted values
/// \note this was out here for plotFilter
// ================================================================================================
int CgWave::resetFrequencyArrays()
{
  Real & tFinal                     = dbase.get<real>("tFinal");
  Real & omega                      = dbase.get<real>("omega");
  Real & damp                       = dbase.get<Real>("damp"); 
  Real & dampSave                   = dbase.get<Real>("dampSave"); 
  Real & Tperiod                    = dbase.get<real>("Tperiod");
  int & numPeriods                  = dbase.get<int>("numPeriods");
  Real & omegaSave                  = dbase.get<real>("omegaSave");
  Real & TperiodSave                = dbase.get<real>("TperiodSave");  
  Real & dt                         = dbase.get<real>("dt");
  Real & dtSave                     = dbase.get<real>("dtSave");

  RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
  RealArray & frequencyArrayAdjusted= dbase.get<RealArray>("frequencyArrayAdjusted");
  RealArray & periodArray           = dbase.get<RealArray>("periodArray"); 
  RealArray & periodArrayAdjusted   = dbase.get<RealArray>("periodArrayAdjusted"); 
  IntegerArray & numPeriodsArray    = dbase.get<IntegerArray>("numPeriodsArray");
  RealArray & frequencyArraySave    = dbase.get<RealArray>("frequencyArraySave");
  RealArray & periodArraySave       = dbase.get<RealArray>("periodArraySave");     

  omega   = omegaSave;
  Tperiod = TperiodSave;
  dt      = dtSave; 
  tFinal  = Tperiod; // reset too *wdh* Jul 25, 2021
  damp    = dampSave;

  frequencyArray = frequencyArraySave;  // reset to original values 
  periodArray    = periodArraySave;     // reset to original values

  return 0;

}


// ----------------------------------------------------------------------------------------------
/// \brief Return the adjusted boundary index
///   (1) return extended boundaries for dirichlet BCs
///    (2) for neumann/characteristic BCs, skip ends that meet adjacent Dirichlet BCs
///
/// \param Ib1,Ib2,Ib3 (output)
// ----------------------------------------------------------------------------------------------
int CgWave::getAdjustedBoundaryIndex( MappedGrid & mg, int side, int axis, Index & Ib1, Index & Ib2, Index & Ib3 )
{
  const int numberOfDimensions = mg.numberOfDimensions();

  IntegerArray abi(2,3); // adjusted boundary index

  if( mg.boundaryCondition(side,axis)==dirichlet )
    abi = mg.dimension();     // return extended boundary for a Dirichlet BC
  else
    abi = mg.gridIndexRange();

  abi(0,axis) = mg.gridIndexRange(side,axis);
  abi(1,axis) = mg.gridIndexRange(side,axis);

  // if( mg.boundaryCondition(side,axis)==dirichlet )
  // {
  //   // return extended boundary for a Dirichlet BC
  //   abi(0,axis) = mg.dimension(side,axis);
  //   abi(1,axis) = mg.dimension(side,axis);
  // }

  if( mg.boundaryCondition(side,axis)==neumann || 
      mg.boundaryCondition(side,axis)==characteristic )
  {
    // On a neumann or characteristic boundary, do not include end points on adjacent dirichlet BCs
    for( int dir=0; dir<numberOfDimensions; dir++ )
    {
      if( dir!=axis )
      {
        for( int side2=0; side2<=1; side2++ )
        {
          if( mg.boundaryCondition(side2,dir)==dirichlet )
          {
            const int is2=1-2*side2;
            abi(side2,dir) += is2;
          }
        }
      }
    }
  }

  Ib1 = Range(abi(0,0),abi(1,0));
  Ib2 = Range(abi(0,1),abi(1,1));
  Ib3 = Range(abi(0,2),abi(1,2));

  return 0;
}

// ================================================================================================
/// \brief Re-initialize deflation (if the number of vectors to deflate has changed, for example)
// ================================================================================================
int CgWave::reinitializeDeflation()
{
  dbase.get<bool>("deflationInitialized") = false;
  return 0;
}

// ================================================================================================
/// \brief Set the name of the composite grid file 
// ================================================================================================
int CgWave::
setNameOfGridFile( aString & nameOfOGFile )
{
  dbase.get<aString>("nameOfGridFile")=nameOfOGFile;
  return 0;
}

// ================================================================================================
/// \brief Get the average number of iterations per implicit solve
// ================================================================================================
Real CgWave::
getAverageNumberOfIterationsPerImplicitSolve() const
{
  const int & totalImplicitIterations       = dbase.get<int>("totalImplicitIterations"); 
  const int & totalImplicitSolves           = dbase.get<int>("totalImplicitSolves"); 
  const Real aveNumberOfImplicitIterations  = totalImplicitIterations/Real(max(1,totalImplicitSolves)); 

  return aveNumberOfImplicitIterations;
}


// ================================================================================================
/// \brief Return the current solution
// ================================================================================================
realCompositeGridFunction& CgWave::getCurrentSolution()
{
  int current=0; // FIX ME 
  realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");

  return u[current];
}

// ================================================================================================
/// \brief Initialize time-step and forcing 
// ================================================================================================
int CgWave::initialize()
{
  real cpu0=getCPU();

  const int & useSuperGrid = dbase.get<int>("useSuperGrid");
  if( useSuperGrid )
  {
    // -- construct superGrid functions ---
    buildSuperGrid();
  }  
  
  const int & numberOfFrequencies  = dbase.get<int>("numberOfFrequencies");
  printF("CgWave::initialize and assign forcing... numberOfFrequencies=%d\n",numberOfFrequencies);

  const ModifiedEquationApproachEnum & modifiedEquationApproach = dbase.get<ModifiedEquationApproachEnum>("modifiedEquationApproach");
  const int & computeEigenmodes        = dbase.get<int>("computeEigenmodes");

  const int & addForcing = dbase.get<int>("addForcing");
  const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
  const int & solveForScatteredField      = dbase.get<int>("solveForScatteredField");

  const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");
  bool twilightZone = addForcing && forcingOption==twilightZoneForcing;
  int & computeErrors = dbase.get<int>("computeErrors");
  int & computeEnergy = dbase.get<int>("computeEnergy");
  int & saveMaxErrors = dbase.get<int>("saveMaxErrors");   
  // computeErrors = twilightZone || knownSolutionOption=="userDefinedKnownSolution";
  
  const real & omega     = dbase.get<real>("omega");
  const real & cfl       = dbase.get<real>("cfl");
  const real & dt        = dbase.get<real>("dt");


  const int & orderOfAccuracy      = dbase.get<int>("orderOfAccuracy");

  int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");
  // By default, order in time = order in space : 
  if( orderOfAccuracyInTime==-1 )
    orderOfAccuracyInTime = orderOfAccuracy;

  int & orderOfExtrapolation = dbase.get<int>("orderOfExtrapolation");
  if( orderOfExtrapolation==-1 )
    orderOfExtrapolation = orderOfAccuracy+1;


  const int & upwind                  = dbase.get<int>("upwind");
  // const real & ad4            = dbase.get<real>("ad4"); // coeff of the artificial dissipation.

  bool useUpwindDissipation = upwind;
  if( useUpwindDissipation )
  {
    // --Upwind dissipation ---

    // upwind dissipation requires an extra ghost line  -- check there are enough
    int orderOfAccuracyInSpace = orderOfAccuracy;
    const int minGhostNeeded = orderOfAccuracyInSpace/2 +1;
    Range Rx=cg.numberOfDimensions();
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg=cg[grid];
      const IntegerArray & numberOfGhostPoints = mg.numberOfGhostPoints();
      int numGhost = min(numberOfGhostPoints(Range(0,1),Rx));
      if( numGhost < minGhostNeeded )
      {
        printF("--CgWave-- initialize: ERROR: the grid does not have enough ghost points for upwind dissipation.\n"
        "   orderOfAccuracy=%i requires at least %i ghost points.\n"
        "   You could remake the grid with more ghost points to fix this error.\n",
        orderOfAccuracyInSpace,minGhostNeeded);
        OV_ABORT("ERROR");
        
      } 
    }

    // --Upwind dissipation ---

    // We need to increase the maximum allowable width to extrap interp neighbours
    // -- see Maxwell -- setupGridFunctions

    int orderOfExtrapolationForInterpolationNeighbours=orderOfAccuracy+1; 
    GenericMappedGridOperators::setDefaultMaximumWidthForExtrapolateInterpolationNeighbours(orderOfExtrapolationForInterpolationNeighbours+1);
  }


  #ifdef USE_PPP
    // In parallel - check that we have enough parallel ghost points 
    const int numParGhost = MappedGrid::getMinimumNumberOfDistributedGhostLines(); 

    assert( upwind==0 || upwind==1 );

    const int numParGhostMin = orderOfAccuracy/2 + upwind;
    printF("\n @@@@@@ cgWave: numParGhost (actual) = %d, numParGhost (minimum needed) = %d @@@@@@\n\n",numParGhost,numParGhostMin);
    if( numParGhost < numParGhostMin )
    {
      printF("CgWave:ERROR: The grid was not created with enough parallel ghost points\n");
      printF("            : Use the command line option -numParGhost=%d to set the number of parallel ghost points\n",numParGhostMin );
      OV_ABORT("ERROR");
    }

  #endif  

  // else if( dialog.getToggleValue(answer,"compute energy",computeEnergy) )
  // {
  //   printF("Setting computeEnergy=%d.\n",computeEnergy);
  //   // numberOfSequences = 2;
  // }

  // else if( dialog.getToggleValue(answer,"save max errors",saveMaxErrors) )
  // {
  //   printF("Setting saveMaxErrors=%d.\n",saveMaxErrors);
  //   // numberOfSequences = 2;
  // }


  // --- sequence info ---
  int & numberOfSequences  = dbase.get<int>("numberOfSequences");

  numberOfSequences=1;
  if( computeEnergy )
  {
    computeEnergy=numberOfSequences;
    numberOfSequences++;
  }
  if( saveMaxErrors )
  {
    saveMaxErrors=numberOfSequences;
    numberOfSequences++;
  }    

  // ------------ Helmholtz -------------
  const int & solveHelmholtz = dbase.get<int>("solveHelmholtz");
  int & computeTimeIntegral = dbase.get<int>("computeTimeIntegral");
  if( solveHelmholtz )
    computeTimeIntegral=true;

  Range all;
  if( computeTimeIntegral )
  {
    // allocate grid functions for time integral and Helmholtz solver ----
    CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");

    realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");
    
    const int & numCompWaveHoltz = dbase.get<int>("numCompWaveHoltz");

    v.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz);
    v.setOperators(operators); 
    v=0.;
    v.setName("v");
    for( int freq=0; freq<numCompWaveHoltz; freq++ )
      v.setName(sPrintF("v%d",freq),freq);

    // realCompositeGridFunction & vOld = dbase.get<realCompositeGridFunction>("vOld");
    // vOld.updateToMatchGrid(cg,all,all,all);

    // realCompositeGridFunction & residual = dbase.get<realCompositeGridFunction>("residual");
    // residual.updateToMatchGrid(cg,all,all,all);
  }
  
  
  // ---- check for radiation boundary conditions ---
  bool & hasRadiationBoundaryConditions = dbase.get<bool>("hasRadiationBoundaryConditions");
  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];
    ForBoundary(side,axis) 
    {
      int bc = mg.boundaryCondition(side,axis);
      if( bc==abcEM2 || bc==characteristic || bc==absorbing ) 
      {
        hasRadiationBoundaryConditions=true;
        break;
      }
    }
  }

  // ============= INITIALIZE FORCING ================
  realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f");
  if( forcingOption != noForcing )
    f.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);
  
  Index I1,I2,I3;

  if( forcingOption==helmholtzForcing && 
      !computeEigenmodes )
  {
    printF("CgWave::initialize: ***ASSIGN HELMHOLTZ FORCING***\n");
    
    int current=0; // not used 
    real t=0.;

    int ipar[10];
    real rpar[10];
    ipar[1]=current;
    rpar[0]=t;
    rpar[1]=dt;
    // -- evaluate the forcing for a Helmholtz solve ---
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      if( !solveForScatteredField )
      {
        ipar[0]=grid;
        userDefinedForcing( f[grid], ipar,rpar );
      }
      else
      {
        f[grid]=0.; 
      }
      // if( true )
      //   ::display(f[grid],"forcing","%6.2f");
    }
    
  }

  

  if( addForcing && forcingOption==twilightZoneForcing )
  {
    const TwilightZoneEnum & twilightZone = dbase.get<TwilightZoneEnum>("twilightZone");
    OGFunction *& tz = dbase.get<OGFunction*>("tz");

    if( twilightZone==polynomial )
    {
      int & degreeInSpace = dbase.get<int>("degreeInSpace");
      const int & degreeInTime =  dbase.get<int>("degreeInTime");
      const int & numberOfDimensions = cg.numberOfDimensions();

      if( degreeInSpace>8 )
      {
        printF("CgWave:ERROR: degreeInSpace=%d is not supported by OGPolyFunction, Reducing to 8.\n",degreeInSpace);
        degreeInSpace=8;
      }
      int numberOfComponentsForTZ=1;
      tz = new OGPolyFunction(degreeInSpace,numberOfDimensions,numberOfComponentsForTZ,degreeInTime);

      const int ndp=max(max(5,degreeInSpace+1),degreeInTime+1);
    
      printF("\n $$$$$$$ setup TZ: build OGPolyFunction: numCompTz=%i degreeSpace=%i, degreeTime=%i ndp=%i $$$$\n",
             numberOfComponentsForTZ,degreeInSpace,degreeInTime,ndp);

      RealArray spatialCoefficientsForTZ(ndp,ndp,ndp,numberOfComponentsForTZ);  
      spatialCoefficientsForTZ=0.;
      RealArray timeCoefficientsForTZ(ndp,numberOfComponentsForTZ);      
      timeCoefficientsForTZ=0.;

      const int degreeInSpaceZ = numberOfDimensions==2 ? 0 : degreeInSpace;
      for( int iz=0; iz<=degreeInSpaceZ; iz++ )
      {
        for( int iy=0; iy<=degreeInSpace; iy++ )
        {
          for( int ix=0; ix<=degreeInSpace; ix++ )
          {
            for( int n=0; n<numberOfComponentsForTZ; n++ )
            {
              // coeff of x^ix * y^iy * z^iz 
              if( ix+iy+iz <= degreeInSpace )
              {
                spatialCoefficientsForTZ(ix,iy,iz,n) = 1./( 1. + ix + 1.5*iy + 1.25*iz + n );
              }
              
            }
            
          }
        }
      }
      
      for( int n=0; n<numberOfComponentsForTZ; n++ )
      {
        for( int i=0; i<ndp; i++ )
          timeCoefficientsForTZ(i,n)= i<=degreeInTime ? 1./(i+1) : 0. ;
      }
      ::display(timeCoefficientsForTZ,"timeCoefficientsForTZ","%6.3f ");

      ((OGPolyFunction*)tz)->setCoefficients( spatialCoefficientsForTZ,timeCoefficientsForTZ );       

    }
    else if( twilightZone==trigonometric )
    {

      const int numberOfComponents=1; 
      RealArray fx( numberOfComponents),fy( numberOfComponents),fz( numberOfComponents),ft( numberOfComponents);
      RealArray gx( numberOfComponents),gy( numberOfComponents),gz( numberOfComponents),gt( numberOfComponents);
      gx=0.;
      gy=0.;
      gz=0.;
      gt=0.;
      RealArray amplitude( numberOfComponents), cc( numberOfComponents);
      amplitude=1.;
      cc=0.;

      RealArray & trigFreq = dbase.get<RealArray>("trigFreq");

      // fx= dbase.get<ArraySimpleFixed<real,4,1,1,1> >("omega")[0];
      // fy =  numberOfDimensions>1 ?  dbase.get<ArraySimpleFixed<real,4,1,1,1> >("omega")[1] : 0.;
      // fz =  numberOfDimensions>2 ?  dbase.get<ArraySimpleFixed<real,4,1,1,1> >("omega")[2] : 0.;
      // ft =  dbase.get<ArraySimpleFixed<real,4,1,1,1> >("omega")[3];

      fx=trigFreq(0); 
      fy=trigFreq(1); 
      fz=trigFreq(2); 
      ft=trigFreq(3); 

      tz = new OGTrigFunction(fx,fy,fz,ft);

      OGTrigFunction & trig = *((OGTrigFunction*)tz);  // cast tz to be an OGTrigFunction
      
      trig.setShifts(gx,gy,gz,gt);
      trig.setAmplitudes(amplitude);
      trig.setConstants(cc);

    }
    else
    {
      OV_ABORT("ERROR: unknown twilightZone");
    }
    
    CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");
    operators.setTwilightZoneFlow(true);
    operators.setTwilightZoneFlowFunction(*tz);
    

  }

  const bool & saveShowFile = dbase.get<bool>("saveShowFile"); 
  if( saveShowFile )
  {
    aString & nameOfShowFile = dbase.get<aString>("nameOfShowFile"); // name of the show file
    int & flushFrequency     = dbase.get<int>("flushFrequency");  // number of solutions per show file 
    Ogshow *& showFile       = dbase.get<Ogshow*>("showFile");

    bool useStreamMode=true;  // show files will be saved compressed
    showFile = new Ogshow( nameOfShowFile,".",useStreamMode );   
    showFile->setFlushFrequency( flushFrequency ); // save this many solutions per sub-showFile
  }

  // -- compute the time-step ---
  getTimeStep(); 
  
  printF("CgWave::initialize: dt=%g\n",dt);

  Real & dtMax = dbase.get<Real>("dtMax"); 
  dtMax = dt; // save the initial timeStep

  timing(timeForInitialize) += getCPU()-cpu0;
  
  // --- build LCBC objects ----
  const BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");
  if( bcApproach==useLocalCompatibilityBoundaryConditions )
  {
    cpu0 = getCPU();
    initializeLCBC();
    timing(timeForInitializeBCs) += getCPU()-cpu0;
  }

  return 0;
}

// ======================================================================================
/// \brief Reset CPU timings to zero:
// ======================================================================================
int CgWave::resetTimings()
{
   for( int i=0; i<maximumNumberOfTimings; i++ )
    timing(i) = 0.;

  return 0;
}


// ======================================================================================
/// \brief Determine if two composite grids match. This function is used to determine if
///   a show file (used for initial conditions or eigenvectors) has the same composite grid.
/// \return value: true if the grids (appear) to match, false otherwise.
// ======================================================================================
bool CgWave::compositeGridsMatch( CompositeGrid & cg, CompositeGrid & cgsf )
{

  bool gridsMatch=true;
  if( cg.numberOfComponentGrids() != cgsf.numberOfComponentGrids() )
  {
    gridsMatch=false;
    printF("INFO: The number of component grids (=%d) in the show file does not match the number in the current grid (=%d)\n",
      cgsf.numberOfComponentGrids(),cg.numberOfComponentGrids());
    
    // OV_ABORT("ERROR");
  }
  else
  {
    // --- check if grids have the same number of points: 
    cg.update(MappedGrid::THEmask ); 
    cgsf.update(MappedGrid::THEmask ); 
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      OV_GET_SERIAL_ARRAY(int,cg[grid].mask(),maskLocal);
      OV_GET_SERIAL_ARRAY(int,cgsf[grid].mask(),sfMaskLocal);  // use mask   
      if( maskLocal.elementCount() != sfMaskLocal.elementCount() )    
      {
        printF("grid=%d: elementCount=%d : showfile elementCount=%d DO NOT MATCH\n",grid,maskLocal.elementCount(),sfMaskLocal.elementCount() );
        gridsMatch=false;
        break;
      }
      else
      {
        printF("grid=%d: elementCount=%d : showfile elementCount=%d MATCH\n",grid,maskLocal.elementCount(),sfMaskLocal.elementCount() );
      }

      const IntegerArray & dwcg   = cg[grid].discretizationWidth();
      const IntegerArray & dwcgsf = cgsf[grid].discretizationWidth();
      if( max(abs(dwcg-dwcgsf)) > 0 )
      {
        printF("grid=%d: discretizationWidths DO NOT MATCH\n",grid);
        gridsMatch=false;
        break;
      }
    }

  }

  return gridsMatch;

}
// ===================================================================================================================
/// \brief Build options dialog for boundary condition options.
/// \param dialog (input) : graphics dialog to use.
///
// ==================================================================================================================
int CgWave::
buildBoundaryConditionOptionsDialog(DialogData & dialog )
{

  BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");
  int & applyKnownSolutionAtBoundaries        = dbase.get<int>("applyKnownSolutionAtBoundaries"); // by default, do NOT apply known solution at boundaries

  const int & orderOfExtrapolation            = dbase.get<int>("orderOfExtrapolation"); // for BCs, -1= default = orderOfAccuracy+1  

  const int & useSuperGrid                    = dbase.get<int>("useSuperGrid");
  const Real & superGridWidth                 = dbase.get<real>("superGridWidth");

  AssignInterpolationNeighboursEnum & assignInterpNeighbours = dbase.get<AssignInterpolationNeighboursEnum>("assignInterpNeighbours");

  // const int & numberOfBCNames = dbase.get<int>("numberOfBCNames");
  // aString *& bcNames = dbase.get<aString*>("bcNames");

  dialog.setOptionMenuColumns(1);

  dialog.addInfoLabel("Use: bcNumber[1|2|...]=[dirichlet|neumann|...]");
  dialog.addInfoLabel("E.g. bcNumber1=dirichlet");

  aString bcApproachLabel[] = {"useDefaultApproachForBCs", "useOneSidedBCs", "useCompatibilityBCs", "useLocalCompatibilityBCs", "" };
  dialog.addOptionMenu("BC approach:", bcApproachLabel, bcApproachLabel, (int)bcApproach );

  aString assignInterpNeighboursLabel[] = {"defaultAssignInterpNeighbours", 
                                           "extrapolateInterpNeighbours", 
                                           "interpolateInterpNeighbours", "" };
  dialog.addOptionMenu("Interp neighbours:", assignInterpNeighboursLabel, assignInterpNeighboursLabel, 
                       (int)assignInterpNeighbours );

  aString tbCommands[] = {
                          "set known on boundaries",
                          "use superGrid",
                          ""};
  int tbState[15];
  tbState[ 0] = applyKnownSolutionAtBoundaries;
  tbState[ 1] = useSuperGrid;

  int numColumns=1;
  dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 

  // // ----- Text strings ------
  const int numberOfTextStrings=20;
  aString textCommands[numberOfTextStrings];
  aString textLabels[numberOfTextStrings];
  aString textStrings[numberOfTextStrings];

  int nt=0;

  textCommands[nt] = "bc=";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "bcName ([dirichlet|evenSymmetry]");  nt++; 

  textCommands[nt] = "superGrid width";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "superGrid width %g",superGridWidth);  nt++; 

  textCommands[nt] = "orderOfExtrapolation";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "order of extrapolation %g",superGridWidth);  nt++;     

  // null strings terminal list
  textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
  dialog.setTextBoxes(textCommands, textLabels, textStrings);




  return 0;
}

//================================================================================
/// \brief: Look for an boundary condition option in the string "answer"
///
/// \param answer (input) : check this command 
///
/// \return return 1 if the command was found, 0 otherwise.
//====================================================================
int CgWave::
getBoundaryConditionOption(const aString & answer, DialogData & dialog, IntegerArray & bcOrig )
{

  int & orderOfExtrapolation = dbase.get<int>("orderOfExtrapolation"); // for BCs, -1= default = orderOfAccuracy+1  
  int & useSuperGrid         = dbase.get<int>("useSuperGrid");
  Real & superGridWidth      = dbase.get<real>("superGridWidth");


  BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");
  int & applyKnownSolutionAtBoundaries = dbase.get<int>("applyKnownSolutionAtBoundaries"); // by default, do NOT apply known solution at boundaries
  AssignInterpolationNeighboursEnum & assignInterpNeighbours = 
                             dbase.get<AssignInterpolationNeighboursEnum>("assignInterpNeighbours");
  const int & numberOfBCNames = dbase.get<int>("numberOfBCNames");
  aString *& bcNames = dbase.get<aString*>("bcNames");

  int found=true; 
  char buff[180];
  aString answer2,line;
  int len=0;

  if( answer=="useDefaultApproachForBCs" || 
           answer=="useOneSidedBCs" || 
           answer=="useCompatibilityBCs" || 
           answer=="useLocalCompatibilityBCs" )
  {
    bcApproach = ( answer=="useDefaultApproachForBCs" ? defaultBoundaryConditionApproach :
                   answer=="useOneSidedBCs"           ? useOneSidedBoundaryConditions : 
                   answer=="useCompatibilityBCs"      ? useCompatibilityBoundaryConditions :
                   answer=="useLocalCompatibilityBCs" ? useLocalCompatibilityBoundaryConditions : 
                                                        defaultBoundaryConditionApproach );
    printF("Setting approach for boundary conditions to %s\n",(const char *)answer);
  }

  else if( answer=="defaultAssignInterpNeighbours" || 
           answer=="extrapolateInterpNeighbours" || 
           answer=="interpolateInterpNeighbours" )
  {
    // For upwind schemes with wider stencils we can extrap or interp unused points next to interpolation points 
    assignInterpNeighbours = ( answer=="defaultAssignInterpNeighbours" ? defaultAssignInterpNeighbours :
                               answer=="extrapolateInterpNeighbours"   ? extrapolateInterpNeighbours :
                                                                         interpolateInterpNeighbours );
    aString assignType;
    if( assignInterpNeighbours==interpolateInterpNeighbours ) 
      assignType = "interpolated";
    else if( assignInterpNeighbours==extrapolateInterpNeighbours )
      assignType = "extrapolated"; 
    else
      assignType = "interpolated"; // new default July 4, 2022

    printF("Interpolation neighbours will be %s for upwind schemes with wider stencils.\n",(const char*)assignType );
  }
 
  else if( dialog.getToggleValue(answer,"set known on boundaries",applyKnownSolutionAtBoundaries) )
  {
    printF("Setting applyKnownSolutionAtBoundaries=%d (1=apply known solution on boundaries).\n",applyKnownSolutionAtBoundaries);
  }   
  else if( dialog.getToggleValue(answer,"use superGrid",useSuperGrid) )
  {
    printF("Setting useSuperGrid=%d (1=use superGrid absorbing layers).\n",useSuperGrid);
  }
  else if( dialog.getTextValue(answer,"superGrid width","%e",superGridWidth) )
  {
    printF("Setting superGridWidth=%g (layer width in parameter space)\n",superGridWidth);
  } 
  else if( dialog.getTextValue(answer,"order of extrapolation","%i",orderOfExtrapolation) )
  {
    printF("Setting orderOfExtrapolation=%i (for boundary conditions, -1=use default=orderOfAccuracy+1)\n",orderOfExtrapolation);
  }    
  else if( len=answer.matches("bc=") )
  {
    aString bcn = answer(len,answer.length()-1);
    int bc=-1;
    for( int i=0; i<numberOfBCNames; i++ )
    {
      if( bcn==bcNames[i] )
      {
        bc=i;
      }
    }
    if( bc==-1 )
    {
      printF("Error: unknown bc=[%s]\n",(const char*)bcn);
      printF("Valid bcNames:\n");
      for( int i=0; i<numberOfBCNames; i++ )
      {
        printF(" bcNames[%i]=[%s]\n",i,(const char*)bcNames[i]);
        OV_ABORT("ERROR: fix boundary conditions");
      }
      
    }
    else
    {
      printF("Setting all boundary conditions to bc=[%s]\n",(const char*)bcNames[bc]);
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        MappedGrid & mg=cg[grid];
        ForBoundary(side,axis) 
        {
          if( mg.boundaryCondition(side,axis)>0 )
          {
            mg.setBoundaryCondition(side,axis,bc);
          }
        }
      }
    }
  
  }

  else if( len=answer.matches("bcNumber") )
  {
    // --- specify a boundary condition for BC's with given numbers ---
    //     bcNumber2=dirichlet
    //     bcNumber5=neumann

    // (1) parse the string to find the number after "bcNumber"
    int i0=len, i1=-1;
    for( int i=i0; i<answer.length(); i++ )
    {
      if( answer[i]=='=' )
      {
        i1=i-1; break; 
      }
    }
    if( i1<i0 )
    {
      printF("getBoundaryConditionOption: ERROR parsing command=[%s]\n",(const char*)answer);
      found=false;
      return found;
    }
    int bcNumber=-1;
    sScanF(answer(i0,i1),"%d",&bcNumber);
    if( bcNumber<=0 )
    {
      printF("ERROR parsing command=[%s], bcNumber=%d is invalid.\n",(const char*)answer,bcNumber);
    }

    // Here is the bc name: 
    aString bcName =answer(i1+2,answer.length()-1);

    printF("answer=[%s], i0=%d, i1=%d, bcNumber=%d, bcName=[%s]\n",(const char*)answer,i0,i1,bcNumber,(const char*)bcName);

    BoundaryConditionEnum bcType = dirichlet;
    if( bcName=="dirichlet" || bcName=="d" )
    {
      bcType=dirichlet;
    }
    else if( bcName=="neumann" || bcName=="n" )
    {
      bcType=neumann;
    }
    else if( bcName=="absorbing" || bcName=="a" )
    {
      bcType=absorbing;
    }    
    else
    {
      printF("ERROR: unknown bcName=[%s]\n",(const char*)bcName);
    }

    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg=cg[grid];
      ForBoundary(side,axis) 
      {
        if( bcOrig(side,axis,grid)==bcNumber ) // note refer to original numbering
        {
          printF("Setting boundaryCondition(side,axis,grid)=(%d,%d,%d) to [%s]\n",side,axis,grid,(const char*)bcNames[bcType]);
          mg.setBoundaryCondition(side,axis,bcType);
        }
      }
    }

    //  OV_ABORT("stop here for now");

  }
 
  else
  {
    found=false;
  }
  

  return found;
}

// ===================================================================================================================
/// \brief Build options dialog for WaveHoltz Options.
/// \param dialog (input) : graphics dialog to use.
///
// ==================================================================================================================
int CgWave::
buildWaveHoltzOptionsDialog(DialogData & dialog )
{

  int & numToDeflate                   = dbase.get<int>("numToDeflate");
  int & solveHelmholtz                 = dbase.get<int>("solveHelmholtz");
  int & adjustHelmholtzForUpwinding    = dbase.get<int>("adjustHelmholtzForUpwinding");
  int & computeTimeIntegral            = dbase.get<int>("computeTimeIntegral");
  int & deflateWaveHoltz               = dbase.get<int>("deflateWaveHoltz");
  int & deflateForcing                 = dbase.get<int>("deflateForcing");
  int & useOptFilter                   = dbase.get<int>("useOptFilter");
  const real & omega                   = dbase.get<real>("omega");

  aString tbCommands[] = {
                          "solve Helmholtz",
                          "compute time integral",
                          "adjust Helmholtz for upwinding", 
                          "deflate WaveHoltz", 
                          "deflate forcing",   
                          "use optimized filter",                         
                            ""};
  int tbState[15];
  tbState[ 0] = solveHelmholtz;
  tbState[ 1] = computeTimeIntegral;
  tbState[ 2] = adjustHelmholtzForUpwinding;
  tbState[ 3] = deflateWaveHoltz;
  tbState[ 4] = deflateForcing;
  tbState[ 5] = useOptFilter;

  int numColumns=1;
  dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 

  // // ----- Text strings ------
  const int numberOfTextStrings=20;
  aString textCommands[numberOfTextStrings];
  aString textLabels[numberOfTextStrings];
  aString textStrings[numberOfTextStrings];

  int nt=0;

  textCommands[nt] = "omega";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",omega);  nt++; 

  textCommands[nt] = "number to deflate";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numToDeflate);  nt++;

  // null strings terminal list
  textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
  dialog.setTextBoxes(textCommands, textLabels, textStrings);




  return 0;
}

//================================================================================
/// \brief: Look for an WaveHoltz option in the string "answer"
///
/// \param answer (input) : check this command 
///
/// \return return 1 if the command was found, 0 otherwise.
//====================================================================
int CgWave::
getWaveHoltzOption(const aString & answer,
                   DialogData & dialog )
{

  int & numToDeflate                   = dbase.get<int>("numToDeflate");
  int & solveHelmholtz                 = dbase.get<int>("solveHelmholtz");
  int & adjustHelmholtzForUpwinding    = dbase.get<int>("adjustHelmholtzForUpwinding");
  int & computeTimeIntegral            = dbase.get<int>("computeTimeIntegral");
  int & deflateWaveHoltz               = dbase.get<int>("deflateWaveHoltz");
  int & deflateForcing                 = dbase.get<int>("deflateForcing");
  int & useOptFilter                   = dbase.get<int>("useOptFilter");
  real & omega                         = dbase.get<real>("omega");
  RealArray & frequencyArray           = dbase.get<RealArray>("frequencyArray");
  RealArray & frequencyArraySave       = dbase.get<RealArray>("frequencyArraySave");
  RealArray & frequencyArrayAdjusted   = dbase.get<RealArray>("frequencyArrayAdjusted");

  int found=true; 
  char buff[180];
  aString answer2,line;
  int len=0;

  if( dialog.getTextValue(answer,"number to deflate","%i",numToDeflate) )
  {
    printF("Setting numToDeflate=%i (deflate this many eigenvectors for WaveHoltz)\n",numToDeflate);
  }
  else if( dialog.getToggleValue(answer,"adjust Helmholtz for upwinding",adjustHelmholtzForUpwinding) )
  {
    printF("Setting adjustHelmholtzForUpwinding=%d (1=correct Helmholtz solution for upwinding)\n",adjustHelmholtzForUpwinding);
  }
  else if( dialog.getToggleValue(answer,"compute time integral",computeTimeIntegral) )
  {
    printF("Setting computeTimeIntegral=%i\n",computeTimeIntegral);
  } 
  else if( dialog.getToggleValue(answer,"deflate WaveHoltz",deflateWaveHoltz) )
  {
    printF("Setting deflateWaveHoltz=%i (1=apply deflate to WaveHoltz iteration\n",deflateWaveHoltz);
  }   
  else if( dialog.getToggleValue(answer,"deflate forcing",deflateForcing) )
  {
    printF("Setting deflateForcing=%i (1=deflate forcing (remove components along deflation eigenvectors) but do not adjust iterates)\n",deflateForcing);
  } 
  else if( dialog.getToggleValue(answer,"use optimized filter",useOptFilter) )
  {
    printF("Setting useOptFilter=%i (1=use optimized filter weights)\n",useOptFilter);
  }        
  else if( dialog.getTextValue(answer,"omega","%e",omega) )
  {
    printF("Setting omega=%g (and frequencyArray(0))\n",omega);

    frequencyArray(0)         = omega;
    frequencyArraySave(0)     = omega;
    frequencyArrayAdjusted(0) = omega;
  }   
  else if( dialog.getToggleValue(answer,"solve Helmholtz",solveHelmholtz) )
  {
    if( solveHelmholtz )
    {
      printF(" solveHelmholtz=true: CgWave is being used to solve the Helmholtz equation (time-periodic wave equation)\n"
             " using CgWaveHoltz\n");
    }
    
  }  
  else
  {
    found=false;
  }
  

  return found;
}


// ===================================================================================================================
/// \brief Build options dialog for EigenWave Options.
/// \param dialog (input) : graphics dialog to use.
///
// ==================================================================================================================
int CgWave::
buildEigenWaveOptionsDialog(DialogData & dialog )
{

  const EigenSolverEnum & eigenSolver        = dbase.get<EigenSolverEnum>("eigenSolver"); 
  const int & computeEigenmodes              = dbase.get<int>("computeEigenmodes");
  const int & numEigsToCompute               = dbase.get<int>("numEigsToCompute"); // number of eigenpairs to compute 
  const int & numArnoldiVectors              = dbase.get<int>("numArnoldiVectors");
  const int & useAccurateInnerProduct        = dbase.get<int>("useAccurateInnerProduct"); 
  const Real & eigenValueTolForMultiplicity  = dbase.get<Real>("eigenValueTolForMultiplicity");
  const int & initialVectorsForEigenSolver   = dbase.get<int>("initialVectorsForEigenSolver");
  const aString & eigenVectorFile            = dbase.get<aString>("eigenVectorFile"); //  name of file holding eigs and eigenvectors for deflation
  const int & minStepsPerPeriod              = dbase.get<int>("minStepsPerPeriod");
  const int & numberOfRitzVectors            = dbase.get<int>("numberOfRitzVectors");
  const int & assignRitzFrequency            = dbase.get<int>("assignRitzFrequency");
  const Real & eigenVectorForcingFactor      = dbase.get<Real>("eigenVectorForcingFactor");
  const int & numPeriods                     = dbase.get<int>("numPeriods");
  const real & omega                         = dbase.get<real>("omega");

  aString eigenSolverLabel[] = {"defaultEigenSolver", "KrylovSchur", "Arnoldi", "Arpack", "fixedPoint","subspaceIteration","" };
  dialog.addOptionMenu("eigenSolver:", eigenSolverLabel, eigenSolverLabel, (int)eigenSolver );

  aString tbCommands[] = {
                          "compute eigenModes",
                          "initial vectors for eigenSolver",
                          "use accurate inner product",
                            ""};
  int tbState[15];
  tbState[0] = computeEigenmodes;
  tbState[1] = initialVectorsForEigenSolver;
  tbState[2] = useAccurateInnerProduct;

  int numColumns=1;
  dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 

  // // ----- Text strings ------
  const int numberOfTextStrings=20;
  aString textCommands[numberOfTextStrings];
  aString textLabels[numberOfTextStrings];
  aString textStrings[numberOfTextStrings];

  int nt=0;

  textCommands[nt] = "omega";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",omega);  nt++; 

  textCommands[nt] = "number of eigenvalues";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numEigsToCompute);  nt++;

  textCommands[nt] = "num Arnoldi vectors";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i (-1 = use default)",numArnoldiVectors);  nt++;

  textCommands[nt] = "number of periods";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numPeriods);  nt++; 

  textCommands[nt] = "min steps per period";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",minStepsPerPeriod);  nt++;

  textCommands[nt] = "number of Ritz vectors";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numberOfRitzVectors);  nt++;  
  
  textCommands[nt] = "assign Ritz frequency";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",assignRitzFrequency);  nt++;  

  textCommands[nt] = "eigenVectorFile";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%s",(const char*)eigenVectorFile);  nt++;

  
  textCommands[nt] = "eigenVectorForcingFactor";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",eigenVectorForcingFactor);  nt++;

  textCommands[nt] = "eig multiplicity tol";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",eigenValueTolForMultiplicity);  nt++; 

 // null strings terminal list
  textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
  dialog.setTextBoxes(textCommands, textLabels, textStrings);




  return 0;
}

//================================================================================
/// \brief: Look for an EigenWave option in the string "answer"
///
/// \param answer (input) : check this command 
///
/// \return return 1 if the command was found, 0 otherwise.
//====================================================================
int CgWave::
getEigenWaveOption(const aString & answer,
                   DialogData & dialog )
{

  EigenSolverEnum & eigenSolver        = dbase.get<EigenSolverEnum>("eigenSolver"); 
  int & computeEigenmodes              = dbase.get<int>("computeEigenmodes");

  int & numEigsToCompute               = dbase.get<int>("numEigsToCompute"); // number of eigenpairs to compute 
  int & numArnoldiVectors              = dbase.get<int>("numArnoldiVectors");
  int & useAccurateInnerProduct        = dbase.get<int>("useAccurateInnerProduct"); 
  Real & eigenValueTolForMultiplicity  = dbase.get<Real>("eigenValueTolForMultiplicity");
  int & initialVectorsForEigenSolver   = dbase.get<int>("initialVectorsForEigenSolver");
  aString & eigenVectorFile            = dbase.get<aString>("eigenVectorFile"); //  name of file holding eigs and eigenvectors for deflation
  int & minStepsPerPeriod              = dbase.get<int>("minStepsPerPeriod");
  int & numberOfRitzVectors            = dbase.get<int>("numberOfRitzVectors");
  int & assignRitzFrequency            = dbase.get<int>("assignRitzFrequency");
  Real & eigenVectorForcingFactor      = dbase.get<Real>("eigenVectorForcingFactor"); 
  int & numPeriods                     = dbase.get<int>("numPeriods");
  real & omega                         = dbase.get<real>("omega");

  int found=true; 
  char buff[180];
  aString answer2,line;
  int len=0;

  if( dialog.getToggleValue(answer,"compute eigenModes",computeEigenmodes) )
  {
    if( computeEigenmodes )
    {
      printF(" computeEigenmodes: Use the WaveHoltz algorithm to compute eigenvalues and eigenvectors.\n");
    }
    
  }  
  else if( dialog.getTextValue(answer,"omega","%e",omega) )
  {
    printF("Setting omega=%g\n",omega);
  }  
  else if( dialog.getToggleValue(answer,"initial vectors for eigenSolver",initialVectorsForEigenSolver) )
  {
    printF("Setting initialVectorsForEigenSolver=%i (vectors used in Arnolid etc.)\n",initialVectorsForEigenSolver);
  } 
  else if( dialog.getToggleValue(answer,"use accurate inner product",useAccurateInnerProduct) )
  {
    printF("Setting useAccurateInnerProduct=%i (for computing the Rayleigh Quotient)\n",useAccurateInnerProduct);
  }  
  
  else if( dialog.getTextValue(answer,"number of eigenvalues","%i",numEigsToCompute) )
  {
    printF("Setting numEigsToCompute=%i (compute this many eigenpairs when using Arnoldi)\n",numEigsToCompute);
  } 
  else if( dialog.getTextValue(answer,"num Arnoldi vectors","%i",numArnoldiVectors) )
  {
    printF("Setting numArnoldiVectors=%i (number of vectors in the Arnoldi column space)\n",numArnoldiVectors);
  }     
  else if( dialog.getTextValue(answer,"number of periods","%i",numPeriods) )
  {
    printF("Setting numPeriods=%i\n",numPeriods);
  }   

  else if( dialog.getTextValue(answer,"min steps per period","%i",minStepsPerPeriod) )
  {
    if( minStepsPerPeriod<5 )
    {
      printF("WARNING: minStepsPerPeriod should be at least 5. Setting to 5\n");
      minStepsPerPeriod=5;
    }
    printF("Setting minStepsPerPeriod=%i (for WaveHoltz)\n",minStepsPerPeriod);
  }
  else if( dialog.getTextValue(answer,"number of Ritz vectors","%i",numberOfRitzVectors) )
  {
    printF("Setting numberOfRitzVectors=%i (max. number of Rayleigh Ritz vectors, used to accelerate EigenWave)\n",numberOfRitzVectors);
  }
  else if( dialog.getTextValue(answer,"assign Ritz frequency","%i",assignRitzFrequency) )
  {
    printF("Setting assignRitzFrequency=%i (set solution to latest Ritz vector every this many steps (EigenWave)\n",assignRitzFrequency);
  }        

  else if( dialog.getTextValue(answer,"eigenVectorFile","%s",eigenVectorFile) )
  {
    printF("Setting eigenVectorFile=%s (file holding eigenvectors for deflation)\n",(const char*)eigenVectorFile);
  }    

  else if( dialog.getTextValue(answer,"eig multiplicity tol","%e",eigenValueTolForMultiplicity) )
  {
    printF("Setting eigenValueTolForMultiplicity=%g (tolerence for detecting mutiple eigenvalues)\n",eigenValueTolForMultiplicity);
  }    
  else if( answer=="defaultEigenSolver" || 
           answer=="KrylovSchur" || 
           answer=="Arnoldi" || answer=="arnoldi" || 
           answer=="Arpack"  || answer=="arpack"  || answer=="ARPACK" ||
           answer=="fixedPoint" || answer=="fixPoint" || answer=="fp" ||
           answer=="subspaceIteration" || answer=="SI" )
  {
    if( answer=="defaultEigenSolver" )
      eigenSolver = defaultEigenSolver;
    else if( answer=="KrylovSchur" )
      eigenSolver = KrylovSchurEigenSolver;
    else if( answer=="Arnoldi" || answer=="arnoldi" )
      eigenSolver = ArnoldiEigenSolver;
   else if( answer=="Arpack" || answer=="arpack"  || answer=="ARPACK" )
      eigenSolver = ArpackEigenSolver;  
   else if( answer=="fixedPoint" || answer=="fp"  || answer=="fixPoint" )
      eigenSolver = fixedPointEigenSolver; 
   else if( answer=="subspaceIteration" || answer=="SI"  )
      eigenSolver = subspaceIterationEigenSolver;                    
    else
    {
      OV_ABORT("This case should not happen");
    }
    printF("Setting eigenSolver = %s\n",(const char*)answer);
  }  
  else if( dialog.getTextValue(answer,"eigenVectorForcingFactor","%e",eigenVectorForcingFactor) )
  {
    printF("Setting eigenVectorForcingFactor=%g\n",eigenVectorForcingFactor);
  }

  else
  {
    found=false;
  }
  

  return found;
}



// ================================================================================================
/// \brief Assign parameters 
// ================================================================================================
int CgWave::interactiveUpdate()
{
//  GenericGraphicsInterface & gi = *Overture::getGraphicsInterface();
  PlotStuffParameters psp;

  real & cfl                           = dbase.get<real>("cfl");
  real & tFinal                        = dbase.get<real>("tFinal");
  real & tPlot                         = dbase.get<real>("tPlot");
  Real & dtMax                         = dbase.get<Real>("dtMax"); 
  Real & damp                          = dbase.get<Real>("damp");

  real & omega                         = dbase.get<real>("omega");
  real & Tperiod                       = dbase.get<real>("Tperiod");
  int & numPeriods                     = dbase.get<int>("numPeriods");
  int & orderOfAccuracy                = dbase.get<int>("orderOfAccuracy");
  int & orderOfAccuracyInTime          = dbase.get<int>("orderOfAccuracyInTime");
  IntegerArray & gridIsImplicit        = dbase.get<IntegerArray>("gridIsImplicit");
  RealArray & bImp                     = dbase.get<RealArray>("bImp");
  RealArray & cImp                     = dbase.get<RealArray>("cImp");
  // int & chooseImplicitTimeStepFromCFL  = dbase.get<int>("chooseImplicitTimeStepFromCFL");
  int & chooseTimeStepFromExplicitGrids= dbase.get<int>("chooseTimeStepFromExplicitGrids");

  int & numberOfFrequencies            = dbase.get<int>("numberOfFrequencies");
  RealArray & frequencyArray           = dbase.get<RealArray>("frequencyArray");
  RealArray & periodArray              = dbase.get<RealArray>("periodArray");  

  int & upwind                         = dbase.get<int>("upwind");
  int & numUpwindCorrections           = dbase.get<int>("numUpwindCorrections");
  int & implicitUpwind                 = dbase.get<int>("implicitUpwind");
  real & ad4                           = dbase.get<real>("ad4"); // coeff of the artificial dissipation. (*old*)
  int & dissipationFrequency           = dbase.get<int>("dissipationFrequency");
  int & preComputeUpwindUt             = dbase.get<int>("preComputeUpwindUt");
  int & takeImplicitFirstStep          = dbase.get<int>("takeImplicitFirstStep");

  int & adjustPlotsForSuperGrid        =  dbase.get<int>("adjustPlotsForSuperGrid");    // set solution to zero in any superGridLayers
  int & adjustErrorsForSuperGrid       =  dbase.get<int>("adjustErrorsForSuperGrid");  


  ModifiedEquationApproachEnum & modifiedEquationApproach = dbase.get<ModifiedEquationApproachEnum>("modifiedEquationApproach");
       
  int & computeErrors                  = dbase.get<int>("computeErrors");                         // by default, compute errors for TZ or a known solution
  int & computeEnergy                  = dbase.get<int>("computeEnergy");                         // by default, compute errors for TZ or a known solution
  int & saveMaxErrors                  = dbase.get<int>("saveMaxErrors");   
  // int & applyKnownSolutionAtBoundaries = dbase.get<int>("applyKnownSolutionAtBoundaries"); // by default, do NOT apply known solution at boundaries
  int & useKnownSolutionForFirstStep   = dbase.get<int>("useKnownSolutionForFirstStep"); 
  int & debug                          = dbase.get<int>("debug");
  int & interactiveMode                = dbase.get<int>("interactiveMode");
  int & checkParameters                = dbase.get<int>("checkParameters");  // 1= check problem parameters for consistency

  int & solveForScatteredField         = dbase.get<int>("solveForScatteredField");
  int & plotScatteredField             = dbase.get<int>("plotScatteredField");
         
  int & solveHelmholtz                 = dbase.get<int>("solveHelmholtz");
  real & tol                           = dbase.get<real>("tol");
       
  bool & saveShowFile                  = dbase.get<bool>("saveShowFile"); 
  aString & nameOfShowFile             = dbase.get<aString>("nameOfShowFile"); // name of the show file
  int & flushFrequency                 = dbase.get<int>("flushFrequency");     // number of solutions per show file  
  int & numberOfSequences              = dbase.get<int>("numberOfSequences");

  real & omegaSOR                      = dbase.get<real>("omegaSOR");

  TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");
  // BoundaryConditionApproachEnum & bcApproach  = dbase.get<BoundaryConditionApproachEnum>("bcApproach");

  TwilightZoneEnum & twilightZone = dbase.get<TwilightZoneEnum>("twilightZone");
  int & degreeInSpace             = dbase.get<int>("degreeInSpace");
  int & degreeInTime              = dbase.get<int>("degreeInTime");
  RealArray & trigFreq            = dbase.get<RealArray>("trigFreq");

  int & addForcing = dbase.get<int>("addForcing");
  ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");

  InitialConditionOptionEnum & initialConditionOption = dbase.get<InitialConditionOptionEnum>("initialConditionOption");

  // AssignInterpolationNeighboursEnum & assignInterpNeighbours = 
  //                            dbase.get<AssignInterpolationNeighboursEnum>("assignInterpNeighbours");
  // const int & numberOfBCNames = dbase.get<int>("numberOfBCNames");
  // aString *& bcNames = dbase.get<aString*>("bcNames");


  // Build a dialog menu for changing parameters
  GUIState gui;
  DialogData & dialog=gui;

  dialog.setWindowTitle("CgWave - Wave Equation Solver");
  dialog.setExitCommand("exit", "exit");




  dialog.setOptionMenuColumns(1);

  aString timeSteppingLabel[] = {"explicit", "implicit", "" };
  dialog.addOptionMenu("time stepping:", timeSteppingLabel, timeSteppingLabel, (int)timeSteppingMethod );

  aString initialConditionLabel[] = {"zero initial condition", 
                                     "twilightZone initial condition", 
                                     "known solution initial condition", 
                                     "pulse initial condition", 
                                     "random initial condition", 
                                     "" };
  dialog.addOptionMenu("Initial condtions:", initialConditionLabel, initialConditionLabel, (int)initialConditionOption );

  aString forcingLabel[] = {"no forcing", "twilightZoneForcing", "userForcing", "helmholtzForcing", "" };
  dialog.addOptionMenu("forcing:", forcingLabel, forcingLabel, (int)forcingOption );

  aString tzLabel[] = {"polynomial", "trigonometric", "" };
  dialog.addOptionMenu("Twilight zone:",tzLabel,tzLabel,(int)twilightZone );

  // aString bcApproachLabel[] = {"useDefaultApproachForBCs", "useOneSidedBCs", "useCompatibilityBCs", "useLocalCompatibilityBCs", "" };
  // dialog.addOptionMenu("BC approach:", bcApproachLabel, bcApproachLabel, (int)bcApproach );

  // aString assignInterpNeighboursLabel[] = {"defaultAssignInterpNeighbours", 
  //                                          "extrapolateInterpNeighbours", 
  //                                          "interpolateInterpNeighbours", "" };
  // dialog.addOptionMenu("Interp neighbours:", assignInterpNeighboursLabel, assignInterpNeighboursLabel, 
  //                      (int)assignInterpNeighbours );

 aString modifiedEquationApproachLabel[] = {"standard modified equation", 
                                            "hierarchical modified equation", 
                                            "stencil modified equation",
                                             "" };
  dialog.addOptionMenu("ME variation:", modifiedEquationApproachLabel, modifiedEquationApproachLabel, 
                       (int)modifiedEquationApproach );

  aString pbLabels[] = {
                        "Boundary Condition Options...",
                        "WaveHoltz Options...",
                        "EigenWave Options...",
                        "user defined known solution...",
                        "user defined forcing...",
                        "choose grids for implicit",
                        "grid",
                        "implicit solver parameters",
                        "erase",
                        "exit",
                        ""};
  int numRows=5;
  dialog.setPushButtons( pbLabels, pbLabels, numRows ); 

  aString tbCommands[] = {
                          "save show file",
                          "upwind dissipation",
                          "turn on forcing",
                          "compute errors",
                          "pre-compute upwind Ut",
                          // "set known on boundaries",
                          // old "choose implicit dt from cfl",
                          "choose dt from explicit grids",
                          "use known for first step",
                          "implicit upwind",
                          "take implicit first step",
                          "adjust plots for superGrid",
                          "adjust errors for superGrid",
                          "solve for scattered field",
                          "plot scattered field",
                          "compute energy",
                          "save max errors",
                          "check parameters",
                            ""};
  int tbState[20];
  tbState[ 0] = saveShowFile;
  tbState[ 1] = upwind;
  tbState[ 2] = addForcing;
  tbState[ 3] = computeErrors;
  tbState[ 4] = preComputeUpwindUt;
  // tbState[ 5] = chooseImplicitTimeStepFromCFL;
  tbState[ 5] = chooseTimeStepFromExplicitGrids;
  tbState[ 6] = useKnownSolutionForFirstStep;
  tbState[ 7] = implicitUpwind;
  tbState[ 8] = takeImplicitFirstStep;
  tbState[ 9] = adjustPlotsForSuperGrid;
  tbState[10] = adjustErrorsForSuperGrid;
  tbState[11] = solveForScatteredField;
  tbState[12] = plotScatteredField;
  tbState[13] = computeEnergy;
  tbState[14] = saveMaxErrors;
  tbState[15] = checkParameters;

  int numColumns=2;
  dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 

  // ----- Text strings ------
  const int numberOfTextStrings=40;
  aString textCommands[numberOfTextStrings];
  aString textLabels[numberOfTextStrings];
  aString textStrings[numberOfTextStrings];

  int nt=0;
  textCommands[nt] = "tFinal";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",tFinal);  nt++; 

  textCommands[nt] = "cfl";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",cfl);  nt++; 

  textCommands[nt] = "show file name";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%s",(const char*)nameOfShowFile);  nt++;

  textCommands[nt] = "flush frequency";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",flushFrequency);  nt++;   

  textCommands[nt] = "orderInTime";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",orderOfAccuracyInTime);  nt++; 

  textCommands[nt] = "degreeInSpace";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",degreeInSpace);  nt++; 

  textCommands[nt] = "degreeInTime";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",degreeInTime);  nt++; 

  textCommands[nt] = "trig frequencies";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g, %g, %g, %g (fx,fy,fz,ft)",trigFreq(0),trigFreq(1),trigFreq(2),trigFreq(3));  nt++; 

  textCommands[nt] = "tPlot";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",tPlot);  nt++; 

  textCommands[nt] = "tol";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",tol);  nt++; 

  textCommands[nt] = "debug";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",debug);  nt++; 

  textCommands[nt] = "interactiveMode";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",interactiveMode);  nt++; 

  // textCommands[nt] = "bc=";  textLabels[nt]=textCommands[nt];
  // sPrintF(textStrings[nt], "bcName ([dirichlet|evenSymmetry]");  nt++; 

  textCommands[nt] = "dissipationFrequency";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",dissipationFrequency);  nt++; 

  textCommands[nt] = "dtMax";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",dtMax);  nt++; 

 textCommands[nt] = "damp";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",damp);  nt++;   

  textCommands[nt] = "implicit weights";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g, %g, %g, %g, %g (beta2,beta4,...)",bImp(0),bImp(1),bImp(2),bImp(3),bImp(4));  nt++; 

  textCommands[nt] = "numUpwindCorrections";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numUpwindCorrections);  nt++; 

  // null strings terminal list
  textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
  dialog.setTextBoxes(textCommands, textLabels, textStrings);


  // ********************* Boundary Condition Options **************************
  DialogData & boundaryConditionOptionsDialog = gui.getDialogSibling();

  boundaryConditionOptionsDialog.setWindowTitle("Boundary Condition Options");
  boundaryConditionOptionsDialog.setExitCommand("close Boundary Condition Options", "close");
  buildBoundaryConditionOptionsDialog( boundaryConditionOptionsDialog );

  // ********************* WaveHoltz Options **************************
  DialogData & waveHoltzOptionsDialog = gui.getDialogSibling();

  waveHoltzOptionsDialog.setWindowTitle("WaveHoltz Options");
  waveHoltzOptionsDialog.setExitCommand("close WaveHoltz Options", "close");
  buildWaveHoltzOptionsDialog( waveHoltzOptionsDialog );

  
  // ********************* EigenWave Options **************************
  DialogData & eigenWaveOptionsDialog = gui.getDialogSibling();

  eigenWaveOptionsDialog.setWindowTitle("EigenWave Options");
  eigenWaveOptionsDialog.setExitCommand("close EigenWave Options", "close");
  buildEigenWaveOptionsDialog( eigenWaveOptionsDialog );

  gi.pushGUI(gui);

  // --- save the original numbering of the BC's ----
  IntegerArray bcOrig(2,3,cg.numberOfComponentGrids());
  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg=cg[grid];
    ForBoundary(side,axis) 
    {
      bcOrig(side,axis,grid)=mg.boundaryCondition(side,axis);
    }
  }

  aString answer,line;
  char buff[200];
  int len;

  psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
  psp.set(GI_TOP_LABEL,sPrintF(buff,"CgWave"));

  for(;;) 
  {

    gi.getAnswer(answer,"");      
    if( answer=="exit" || answer=="continue" )
    {
      break;
    }
    else if( answer.matches("erase") )
    {
      gi.erase();
    }
    else if( answer=="user defined known solution..." )
    {
      if( updateUserDefinedKnownSolution() )
      {
        aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");

        knownSolutionOption = "userDefinedKnownSolution";

        initialConditionOption = knownSolutionInitialCondition;
      }
      
    }
    else if( answer=="user defined forcing..." )
    {
      setupUserDefinedForcing();
    }
    
    else if( answer.matches("grid") )
    {
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      PlotIt::plot(gi,cg,psp);                          // plot the grid
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
    }

    else if( dialog.getTextValue(answer,"tFinal","%e",tFinal) )
    {
      printF("Setting tFinal=%g\n",tFinal);
    }

    else if( len=answer.matches("cfl") )
    {
      sScanF(answer(len,answer.length()-1),"%e",&cfl);
      printF("setting cfl=%g\n",cfl);
    }

    else if( dialog.getTextValue(answer,"show file name","%s",nameOfShowFile) )
    {
      printF("Setting nameOfShowFile=%s\n",(const char*)nameOfShowFile);
    }

    else if( dialog.getTextValue(answer,"debug","%i",debug) )
    {
      printF("Setting debug=%i\n",debug);
    }

    else if( dialog.getTextValue(answer,"interactiveMode","%i",interactiveMode) )
    {
      printF("Setting interactiveMode=%i (0=normal, 1=advance and finish)\n",interactiveMode);
    }

    else if( dialog.getTextValue(answer,"orderInTime","%i",orderOfAccuracyInTime) )
    {
      printF("Setting orderOfAccuracyInTime=%i (-1 means use same as order in space)\n",orderOfAccuracyInTime);
    }

    else if( dialog.getTextValue(answer,"degreeInSpace","%i",degreeInSpace) )
    {
      printF("Setting degreeInSpace=%i\n",degreeInSpace);
    }
    else if( dialog.getTextValue(answer,"degreeInTime","%i",degreeInTime) )
    {
      printF("Setting degreeInTime=%i\n",degreeInTime);
    }

    else if( dialog.getTextValue(answer,"dissipation frequency","%i",dissipationFrequency) )
    {
      printF("Setting dissipationFrequency=%i (dissipation is applied every this many steps)\n",dissipationFrequency);
    }
    else if( dialog.getTextValue(answer,"numUpwindCorrections","%i",numUpwindCorrections) )
    {
      printF("Setting numUpwindCorrections=%i (number of upwind dissipation iterations)\n",numUpwindCorrections);
    }    

    else if( answer=="Boundary Condition Options..."  )
    {
      boundaryConditionOptionsDialog.showSibling();
    }
    else if( answer=="close Boundary Condition Options" )
    {
      boundaryConditionOptionsDialog.hideSibling();  
    }
    else if( getBoundaryConditionOption(answer,boundaryConditionOptionsDialog,bcOrig ) )
    {
      printF("Answer=%s found in getBoundaryConditionOption\n",(const char*)answer);
    }


    else if( answer=="WaveHoltz Options..."  )
    {
      waveHoltzOptionsDialog.showSibling();
    }
    else if( answer=="close WaveHoltz Options" )
    {
      waveHoltzOptionsDialog.hideSibling();  
    }
    else if( getWaveHoltzOption(answer,waveHoltzOptionsDialog ) )
    {
      printF("Answer=%s found in getWaveHoltzOption\n",(const char*)answer);
    }

    else if( answer=="EigenWave Options..."  )
    {
      eigenWaveOptionsDialog.showSibling();
    }
    else if( answer=="close EigenWave Options" )
    {
      eigenWaveOptionsDialog.hideSibling();  // pop timeStepping
    }
    else if( getEigenWaveOption(answer,eigenWaveOptionsDialog ) )
    {
      printF("Answer=%s found in getEigenWaveOption\n",(const char*)answer);
    }

    else if( dialog.getTextValue(answer,"tol","%e",tol) )
    {
      printF("Setting tol=%g (tolerence for Krylov solvers)\n",tol);
    }

    else if( dialog.getTextValue(answer,"dtMax","%e",dtMax) )
    {
      printF("Setting dtMax=%g\n",dtMax);
    }
    else if( dialog.getTextValue(answer,"damp","%e",damp) )
    {
      printF("Setting damp=%g (coefficient of linear damping)\n",damp);
    }       

    else if( dialog.getTextValue(answer,"flush frequency","%i",flushFrequency) )
    {
      printF("Setting flushFrequency=%i (number of solution in each sun showFile)\n",flushFrequency);
    }

    else if( dialog.getToggleValue(answer,"check parameters",checkParameters) )
    {
      printF("Setting checkParameters=%d (1=check parmeters for consistenty)\n",checkParameters);
    }

    else if( dialog.getToggleValue(answer,"save show file",saveShowFile) )
    {
      printF("Setting saveShowFile=%d\n",saveShowFile);
    }

    else if( dialog.getToggleValue(answer,"compute errors",computeErrors) )
    {
      printF("Setting computeErrors=%d (1=compute errors for TZ or known solutions.\n",computeErrors);
    }

    else if( dialog.getToggleValue(answer,"compute energy",computeEnergy) )
    {
      printF("Setting computeEnergy=%d.\n",computeEnergy);
      // numberOfSequences = 2;
    }

    else if( dialog.getToggleValue(answer,"save max errors",saveMaxErrors) )
    {
      printF("Setting saveMaxErrors=%d.\n",saveMaxErrors);
      // numberOfSequences = 2;
    }

    // else if( dialog.getToggleValue(answer,"set known on boundaries",applyKnownSolutionAtBoundaries) )
    // {
    //   printF("Setting applyKnownSolutionAtBoundaries=%d (1=apply known solution on boundaries).\n",applyKnownSolutionAtBoundaries);
    // } 

    else if( dialog.getToggleValue(answer,"use known for first step",useKnownSolutionForFirstStep) )
    {
      printF("Setting useKnownSolutionForFirstStep (use known solution for first time step if possible).\n",
           useKnownSolutionForFirstStep);
    }

    else if( dialog.getToggleValue(answer,"upwind dissipation",upwind) )
    {
      printF("Setting upwind=%d (upwind dissipation is on or off).\n",upwind);
    }
    else if( dialog.getToggleValue(answer,"turn on forcing",addForcing) ){}//

    else if( dialog.getToggleValue(answer,"pre-compute upwind Ut",preComputeUpwindUt) )
    {
      printF("Setting preComputeUpwindUt=%i\n",preComputeUpwindUt);
      printF("        : true  = precompute Ut in upwind dissipation.\n"
            "         : false = compute Ut inline in Gauss-Seidel fashion (allows a bigger upwind coefficient).\n");
    }   
    // else if( dialog.getToggleValue(answer,"choose implicit dt from cfl",chooseImplicitTimeStepFromCFL) )
    // {
    //   printF("Setting chooseImplicitTimeStepFromCFL=%i (1=choose implicit dt from cfl, 0=choose dt from dtMax)\n",chooseImplicitTimeStepFromCFL);
    // }

    else if( dialog.getToggleValue(answer,"choose time step from explicit grids",chooseTimeStepFromExplicitGrids) ||
             dialog.getToggleValue(answer,"choose dt from explicit grids",chooseTimeStepFromExplicitGrids) )
    {
      printF("Setting chooseTimeStepFromExplicitGrids=%i\n"
             "  1=choose dt from cfl and explicit grids only (or if all grids are implicit)\n"
             "  0=choose dt from cfl and all grids\n",    chooseTimeStepFromExplicitGrids);
    }    

    else if( dialog.getToggleValue(answer,"implicit upwind",implicitUpwind) )
    {
      printF("Setting implicitUpwind=%i (1=include upwinding in implicit matrix when implicit time-stepping\n",implicitUpwind);
    } 

    else if( dialog.getToggleValue(answer,"take implicit first step",takeImplicitFirstStep) )
    {
      printF("Setting takeImplicitFirstStep=%i (1=take an implicit first step when implicit time-stepping\n",takeImplicitFirstStep);
    } 

    else if( dialog.getToggleValue(answer,"adjust plots for superGrid",adjustPlotsForSuperGrid) )
    {
      printF("Setting adjustPlotsForSuperGrid=%i (1=do not plot solution in the superGrid layers\n",adjustPlotsForSuperGrid);
    }

    else if( dialog.getToggleValue(answer,"adjust errors for superGrid",adjustErrorsForSuperGrid) )
    {
      printF("Setting adjustErrorsForSuperGrid=%i (1=do not compute errors in the superGrid layers\n",adjustErrorsForSuperGrid);
    }

    else if( dialog.getToggleValue(answer,"plot scattered field",plotScatteredField) )
    {
      printF("Setting plotScatteredField=%i.\n",plotScatteredField);
    }             
    else if( dialog.getToggleValue(answer,"solve for scattered field",solveForScatteredField) )
    {
      printF("Setting solveForScatteredField=%i.\n",solveForScatteredField);
    }    
    else if( answer=="explicit" || answer=="implicit" )
    {
      timeSteppingMethod = ( answer=="explicit" ? explicitTimeStepping :
                             answer=="implicit" ? implicitTimeStepping : explicitTimeStepping );

      if( timeSteppingMethod == explicitTimeStepping )
        gridIsImplicit = 0;
      else
        gridIsImplicit = 1;

    } 

    else if( len=answer.matches("implicit weights") )
    {
      sScanF(answer(len,answer.length()-1),"%e %e %e %e %e",&bImp(0),&bImp(1),&bImp(2),&bImp(3),&bImp(4));
      printF("Setting implicit time-stepping weights to beta2=%g, beta4=%g, beta6=%g, beta8=%g\n",bImp(0),bImp(1),bImp(2),bImp(2));
      
      Real beta2=bImp(0), beta4=bImp(1); 
      bImp(0)=beta2;
      bImp(1)=beta4; 
      Real alpha2 = (1.-beta2)/2.;
      Real alpha4 = (alpha2-beta4-1./12.)/2.; 
      cImp(-1,0)=alpha2;
      cImp( 0,0)= beta2;
      cImp( 1,0)=alpha2;

      cImp(-1,1)=alpha4;
      cImp( 0,1)= beta4;
      cImp( 1,1)=alpha4;

    } 

    // aString initialConditionLabel[] = {"zero initial condition", "twilightZone initial condition", "known solution initial condition", "pulse initial condition", "" };
    else if( answer=="zero initial condition"           || 
             answer=="twilightZone initial condition"   || 
             answer=="known solution initial condition" || 
             answer=="pulse initial condition"          ||
             answer=="random initial condition" )
    {
      initialConditionOption = ( answer=="zero initial condition"           ? zeroInitialCondition :
                                 answer=="twilightZone initial condition"   ? twilightZoneInitialCondition :
                                 answer=="known solution initial condition" ? knownSolutionInitialCondition :
                                 answer=="pulse initial condition"          ? pulseInitialCondition : 
                                 answer=="random initial condition"         ? randomInitialCondition : 
                                                                              zeroInitialCondition );
      printF("Setting initialConditionOption=%s\n",(const char*)answer);
    }


    else if( answer=="noForcing" || answer=="twilightZoneForcing" || answer=="userForcing" || answer=="helmholtzForcing" )
    {
      forcingOption = ( answer=="noForcing" ? noForcing :
                        answer=="twilightZoneForcing" ? twilightZoneForcing :
                        answer=="userForcing" ? userForcing :
                        answer=="helmholtzForcing" ? helmholtzForcing : noForcing);

      if( forcingOption==twilightZoneForcing )
        initialConditionOption = twilightZoneInitialCondition;

      if( forcingOption != noForcing )
        addForcing=true; // *wdh* March 27, 2023
      else
        addForcing=false;
    }

    else if( answer=="polynomial" || answer=="trigonometric" )
    {
      twilightZone = ( answer=="polynomial" ? polynomial :
                       answer=="trigonometric" ? trigonometric : polynomial);
    }

    else if( len=answer.matches("trig frequencies") )
    {
      sScanF(answer(len,answer.length()-1),"%e %e %e %e",&trigFreq(0),&trigFreq(1),&trigFreq(2),&trigFreq(3));
      printF("Setting trig frequencies: fx=%g, fy=%g, fz=%g, ft=%g\n",trigFreq(0),trigFreq(1),trigFreq(2),trigFreq(3));
    }

    // else if( answer=="useDefaultApproachForBCs" || 
    //          answer=="useOneSidedBCs" || 
    //          answer=="useCompatibilityBCs" || 
    //          answer=="useLocalCompatibilityBCs" )
    // {
    //   bcApproach = ( answer=="useDefaultApproachForBCs" ? defaultBoundaryConditionApproach :
    //                  answer=="useOneSidedBCs"           ? useOneSidedBoundaryConditions : 
    //                  answer=="useCompatibilityBCs"      ? useCompatibilityBoundaryConditions :
    //                  answer=="useLocalCompatibilityBCs" ? useLocalCompatibilityBoundaryConditions : 
    //                                                       defaultBoundaryConditionApproach );
    //   printF("Setting approach for boundary conditions to %s\n",(const char *)answer);
    // }

    // else if( answer=="defaultAssignInterpNeighbours" || 
    //          answer=="extrapolateInterpNeighbours" || 
    //          answer=="interpolateInterpNeighbours" )
    // {
    //   // For upwind schemes with wider stencils we can extrap or interp unused points next to interpolation points 
    //   assignInterpNeighbours = ( answer=="defaultAssignInterpNeighbours" ? defaultAssignInterpNeighbours :
    //                              answer=="extrapolateInterpNeighbours"   ? extrapolateInterpNeighbours :
    //                                                                        interpolateInterpNeighbours );
    //   aString assignType;
    //   if( assignInterpNeighbours==interpolateInterpNeighbours ) 
    //     assignType = "interpolated";
    //   else if( assignInterpNeighbours==extrapolateInterpNeighbours )
    //     assignType = "extrapolated"; 
    //   else
    //     assignType = "interpolated"; // new default July 4, 2022

    //   printF("Interpolation neighbours will be %s for upwind schemes with wider stencils.\n",(const char*)assignType );
    // }

    else if( answer=="standard modified equation" || 
             answer=="hierarchical modified equation" ||
             answer=="stencil modified equation" )
    {
      // For upwind schemes with wider stencils we can extrap or interp unused points next to interpolation points 
      modifiedEquationApproach = ( answer=="standard modified equation"     ? standardModifiedEquation :
                                   answer=="hierarchical modified equation" ? hierarchicalModifiedEquation :
                                   answer=="stencil modified equation"      ? stencilModifiedEquation :
                                                                              standardModifiedEquation );
      printF("Setting modifiedEquationApproach=%s\n",(const char*)answer );
    }

    else if( answer=="implicit solver parameters" )
    {
      printf("Set the Oges parameters for the implicit solver.\n");
      if( !dbase.has_key("implicitSolverParameters") )
      {
        OgesParameters & par = dbase.put<OgesParameters>("implicitSolverParameters");
      }

      OgesParameters & par = dbase.get<OgesParameters>("implicitSolverParameters");
      par.update( gi,cg );

    }
    else if( len=answer.matches("tPlot") )
    {
      sScanF(answer(len,answer.length()-1),"%e",&tPlot);
      printF(" tPlot=%g\n",tPlot);
      dialog.setTextLabel("tPlot",sPrintF(line, "%g",tPlot));
    }
    else if( len=answer.matches("artificial dissipation") ) // **OLD WAY**
    {
      sScanF(answer(len,answer.length()-1),"%e",&ad4);
      printF(" artificial diffusion=%g\n",ad4);
      if( ad4>0. )
        upwind=1;
      else
        upwind=0;

      // dialog.setTextLabel("artificial dissipation",sPrintF(line, "%g",ad4));
    }

    // else if( len=answer.matches("bc=") )
    // {
    //   aString bcn = answer(len,answer.length()-1);
    //   int bc=-1;
    //   for( int i=0; i<numberOfBCNames; i++ )
    //   {
    //     if( bcn==bcNames[i] )
    //     {
    //       bc=i;
    //     }
    //   }
    //   if( bc==-1 )
    //   {
    //     printF("Error: unknown bc=[%s]\n",(const char*)bcn);
    //     printF("Valid bcNames:\n");
    //     for( int i=0; i<numberOfBCNames; i++ )
    //     {
    //       printF(" bcNames[%i]=[%s]\n",i,(const char*)bcNames[i]);
    //       OV_ABORT("ERROR: fix boundary conditions");
    //     }
        
    //   }
    //   else
    //   {
    //     printF("Setting all boundary conditions to bc=[%s]\n",(const char*)bcNames[bc]);
    //     for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //     {
    //       MappedGrid & mg=cg[grid];
    //       ForBoundary(side,axis) 
    //       {
    //         if( mg.boundaryCondition(side,axis)>0 )
    //         {
    //           mg.setBoundaryCondition(side,axis,bc);
    //         }
    //       }
    //     }
    //   }
    
    // }

    // else if( len=answer.matches("bcNumber") )
    // {
    //   // --- specify a boundary condition for BC's with given numbers ---
    //   //     bcNumber2=dirichlet
    //   //     bcNumber5=neumann

    //   // (1) parse the string to find the number after "bcNumber"
    //   int i0=len, i1=-1;
    //   for( int i=i0; i<answer.length(); i++ )
    //   {
    //     if( answer[i]=='=' )
    //     {
    //       i1=i-1; break; 
    //     }
    //   }
    //   if( i1<i0 )
    //   {
    //     printF("ERROR parsing command=[%s]\n",(const char*)answer);
    //     continue;
    //   }
    //   int bcNumber=-1;
    //   sScanF(answer(i0,i1),"%d",&bcNumber);
    //   if( bcNumber<=0 )
    //   {
    //     printF("ERROR parsing command=[%s], bcNumber=%d is invalid.\n",(const char*)answer,bcNumber);
    //   }

    //   // Here is the bc name: 
    //   aString bcName =answer(i1+2,answer.length()-1);

    //   printF("answer=[%s], i0=%d, i1=%d, bcNumber=%d, bcName=[%s]\n",(const char*)answer,i0,i1,bcNumber,(const char*)bcName);

    //   BoundaryConditionEnum bcType = dirichlet;
    //   if( bcName=="dirichlet" || bcName=="d" )
    //   {
    //     bcType=dirichlet;
    //   }
    //   else if( bcName=="neumann" || bcName=="n" )
    //   {
    //     bcType=neumann;
    //   }
    //   else
    //   {
    //     printF("ERROR: unknown bcName=[%s]\n",(const char*)bcName);
    //   }

    //   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //   {
    //     MappedGrid & mg=cg[grid];
    //     ForBoundary(side,axis) 
    //     {
    //       if( bcOrig(side,axis,grid)==bcNumber ) // note refer to original numbering
    //       {
    //         printF("Setting boundaryCondition(side,axis,grid)=(%d,%d,%d) to [%s]\n",side,axis,grid,(const char*)bcNames[bcType]);
    //         mg.setBoundaryCondition(side,axis,bcType);
    //       }
    //     }
    //   }

    //   //  OV_ABORT("stop here for now");

    // }

    else if( answer=="choose grids for implicit" )
    {

      if( timeSteppingMethod != implicitTimeStepping )
      {
        printF("WARNING: This option only currently only works for timeSteppingMethod==implicit! \n");
      }
      printF("Set grids to be implicit, or explicit. Type 'done' to finish. \n"
             " Examples: (type `help' for more examples)\n"
             "  square=explicit \n"
             "  annulus=implicit \n"
             "  all=implicit \n" );

      gi.appendToTheDefaultPrompt("implicit>");
      aString gridName;
      int implicit;// implicit used to be a bool.
      bool setRectangularGrids=false;
      for( ;; )
      {
        aString answer2;
        gi.inputString(answer2,"Specify grids to be implicit, semi-implicit or explicit. (or type `help' or `done')");
        if( answer2=="done" || answer2=="exit" )
          break;
        else if( answer2=="help" )
        {
          printF("Specify grids to be implicit, semi-implicit or explicit. Type a string of the form     \n"
                 "                                                                             \n"
                 "       <grid name>=[explicit|implicit]\n"
                 "                                                                             \n"
                 " By default all grids are implicit.                                          \n"
                 " Examples: \n"
                 "     square=explicit                            \n"
                 "     annulus=implicit                           \n"
                 "     all=implicit                               \n"
                 "     rectangular=explicit   (set all rectangular grids to be explicit)  \n"
            );
          printF("Here are the names of the grids: \n");
          for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            printF(" grid %i : name=%s \n",grid,(const char*)cg[grid].getName());
        }
        else
        {
          int length=answer2.length();
          int i,mark=-1;
          for( i=0; i<length; i++ )
          {
            if( answer2[i]=='=' )
            {
              mark=i-1;
              break;
            }
          }
          if( mark<0 )
          {
            printF("unknown form of answer=[%s]. Try again or type `help' for examples.\n",(const char *)answer2);
            gi.stopReadingCommandFile();
            continue;
          }
          else
          {
            gridName=answer2(0,mark);  // this is the name of the grid or `all'
            Range G(-1,-1);
            if( gridName=="all" )
              G=Range(0,cg.numberOfComponentGrids()-1);
            else if( gridName=="rectangular" )
            {
              setRectangularGrids=true;
              G=Range(0,cg.numberOfComponentGrids()-1);
            }
            else
            {
              // search for the name of the grid
              for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
              {
                if( gridName==cg[grid].getName() )
                {
                  G=Range(grid,grid);
                  break;
                }
              }
            }
            if( G.getBase()==-1  )
            {
              printF("Unknown grid name = <%s> \n",(const char *)gridName);
              gi.stopReadingCommandFile();
              continue;
            }

            //This could probably be done in a better way! 
            if( answer2(mark+2,mark+3)=="im" )
            {
              implicit=1;
            }
            else 
            {
              implicit=0;
            }
                     
            for( int grid=G.getBase(); grid<=G.getBound(); grid++ )
            {
              if( setRectangularGrids )
              { // set all rectangular grids 
                if( !cg[grid].isRectangular() )
                  continue;
              }
              if( timeSteppingMethod == implicitTimeStepping )
              {
                gridIsImplicit(grid)=implicit;
                printF("Setting time stepping to be %s for grid %s\n",( implicit==1 ? "implicit" : "explicit" ),(const char*)cg[grid].getName());
              }
              else
              {
                printF("NOT setting time stepping to be %s for grid %s since TIME-STEPPING IS EXPLCICIT.\n",( implicit==1 ? "implicit" : "explicit" ),(const char*)cg[grid].getName());
              }
            }
          }
        }
      }
      gi.unAppendTheDefaultPrompt();
    }

    else
    {
      printF("CgWave:ERROR: unknown answer=[%s]\n",(const char*)answer);
    }
    
  }
  
  gi.popGUI();  // pop dialog

  // Initialize time-step and forcing 
  initialize();

  // ----- output the header ----
  outputHeader();
  
  return 0;
}


// ================================================================================================
/// \brief Setup grids and grid functions
///
/// This function will perform various setup steps such as to update operators, create an interpolant, update grid functions.
// ================================================================================================
int CgWave::setup()
{
  real cpu0 = getCPU();
  
  realCompositeGridFunction *& ucg = dbase.get<realCompositeGridFunction*>("ucg");

  CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");
  operators.updateToMatchGrid(cg);

  int & orderOfAccuracy = dbase.get<int>("orderOfAccuracy");
  int maxDw = max(cg[0].discretizationWidth());
  orderOfAccuracy = maxDw-1;
  assert( orderOfAccuracy==2 || orderOfAccuracy==4 || orderOfAccuracy==6 || orderOfAccuracy==8 );
  printF("CgWave::setup SETTING orderOfAccuracy=%i\n",orderOfAccuracy);

  operators.setOrderOfAccuracy(orderOfAccuracy);

        
  Interpolant *& pInterpolant = dbase.get<Interpolant*>("pInterpolant");
  if( pInterpolant==NULL )
  {
    pInterpolant = new Interpolant(cg);
    pInterpolant->incrementReferenceCount();
    // pInterpolant->updateToMatchGrid(cg);
  }
  
  const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");    

  Range all;
  for( int i=0; i<numberOfTimeLevelsStored; i++ )
  {
    ucg[i].updateToMatchGrid(cg,all,all,all);
    ucg[i]=0.; 
    ucg[i].setOperators(operators);                                 
    ucg[i].setName("u");                       // name the grid function
  }
  
  // gridIsImplicit(grid) = 1 : this grid is advanced implicitly
  IntegerArray & gridIsImplicit = dbase.get<IntegerArray>("gridIsImplicit");
  gridIsImplicit.redim(cg.numberOfComponentGrids());
  gridIsImplicit=0;

  // superGrid(grid) = 1 if this grid has superGrid layers
  IntegerArray & superGrid =  dbase.get<IntegerArray>("superGrid");
  superGrid.redim(cg.numberOfComponentGrids());
  superGrid=0;




  timing(timeForInitialize) += getCPU()-cpu0;

  return 0;
}






// // ======================================================================================================
// /// \brief Compute errors
// // ======================================================================================================
// real CgWave::
// getErrors( realCompositeGridFunction & u, real t )
// {
//   real cpu0=getCPU();
//   // printF("+++++++++++ getErrors t=%9.3e +++++++++++\n",t);

//   real & maxError     = dbase.get<real>("maxError");      // save max-error here 
//   real & solutionNorm = dbase.get<real>("solutionNorm");  // save solution norm here


//   maxError=0.;

//   const int & addForcing = dbase.get<int>("addForcing");
//   const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
//   const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");

//   bool twilightZone = addForcing && forcingOption==twilightZoneForcing;

//   int & computeErrors       = dbase.get<int>("computeErrors");
//   // int & plotHelmholtzErrors = dbase.get<int>("plotHelmholtzErrors");

//   computeErrors = computeErrors && (twilightZone || knownSolutionOption=="userDefinedKnownSolution");
  
//   //  computeErrors = twilightZone || knownSolutionOption=="userDefinedKnownSolution";
  

//   if( computeErrors )
//   {
//     realCompositeGridFunction & error = dbase.get<realCompositeGridFunction>("error");
//     error.updateToMatchGrid(cg);

//     const int numberOfDimensions = cg.numberOfDimensions();
//     Index I1,I2,I3;
//     for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
//     {
//       MappedGrid & mg = cg[grid];
    
//       // get the local serial arrays
//       OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
//       OV_GET_SERIAL_ARRAY(real,error[grid],errLocal);
//       errLocal=0.;

//       getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points.

//       if( knownSolutionOption=="userDefinedKnownSolution" )
//       {
//         // printF("+++++++++++ getErrors for userDefinedKnownSolution +++++++++++\n");
        

//         getUserDefinedKnownSolution( t, grid, error[grid], I1,I2,I3 ); // store true solution in error[grid]
//         const int includeGhost=1;
//         bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
//         if( ok )
//         {
//           errLocal(I1,I2,I3) -= uLocal(I1,I2,I3);
//         }
//       }
//       else
//       {
//         // ----- Twilight zone ------
//         assert( dbase.get<OGFunction*>("tz")!=NULL );
//         OGFunction & e = *dbase.get<OGFunction*>("tz");

//         mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
//         OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);


//         const int includeGhost=1;
//         bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
//         if( ok )
//         {
//           int numberOfComponents=1;
//           Range C=numberOfComponents;
//           int isRectangular=0;
//           RealArray ue(I1,I2,I3);
//           e.gd( ue ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,C,t);

//           errLocal(I1,I2,I3) = ue(I1,I2,I3) - uLocal(I1,I2,I3);

//         }
//       }
      
//     }
//     maxError = maxNorm(error);
//     // printF("getErrors: t=%9.3e, maxError=%9.3e\n",t,maxError);


//   }
//   else if( knownSolutionOption=="userDefinedKnownSolution" )
//   {
//   }
  
//   // compute the solution norm
//   solutionNorm = maxNorm(u);
  

//   timing(timeForGetError)+= getCPU()-cpu0;

//   return maxError;
  


// }

// //=================================================================================================
// /// \brief Take the first step using Taylor series in time (e.g. for Helmholtz solve) 
// /// THIS ASSUMES A HELMHOLTZ SOLVE OR HELMHOLTZ FORCING 
// //=================================================================================================
// int CgWave::
// takeFirstStep( int cur, real t )
// {

//   assert( t==0. );

//   const int & debug           = dbase.get<int>("debug");
//   if( debug & 4 )
//     printF("*******  CgWave::takeFirstStep GET SOLUTION at dt *************\n");
  
//   const real & c              = dbase.get<real>("c");
//   const real & dt             = dbase.get<real>("dt");
//   const real & omega          = dbase.get<real>("omega");
//   const int & orderOfAccuracy = dbase.get<int>("orderOfAccuracy");

//   const int & numberOfFrequencies   = dbase.get<int>("numberOfFrequencies");
//   const RealArray & frequencyArray  = dbase.get<RealArray>("frequencyArray");  

//   const int & solveHelmholtz        = dbase.get<int>("solveHelmholtz");
//   ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
//   const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

//   // Do Helmholtz case for now: 

//   assert( forcingOption==helmholtzForcing );

//   const bool usePeriodicFirstStep=true;


//   //  ---- NOTE: change sign of forcing for Helmholtz since we want to solve ----
//   //       ( omega^2 I + c^2 Delta) w = f  
//   const Real fSign = forcingOption==helmholtzForcing ? -1.0 : 1.0;

//   // const int & solveHelmholtz = dbase.get<int>("solveHelmholtz");
//   // const Real fSign = solveHelmholtz ? -1.0 : 1.0;  
  

//   const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");    
//   const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
//   const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;


//   realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
//   realCompositeGridFunction & uc = u[cur];     // current time 
//   realCompositeGridFunction & un = u[next];    // next time

//   // forcing: 
//   realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f");

//   CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");

//   // ---- WaveHoltz: initial condition is provided in v: ----

//   realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");

//   Index I1,I2,I3;
//   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
//   {
//     MappedGrid & mg = cg[grid];
//     OV_GET_SERIAL_ARRAY(real,un[grid],unLocal);

//     if( usePeriodicFirstStep )
//     {
//       // do all points including ghost 
//       getIndex(mg.dimension(),I1,I2,I3);
//       bool ok=ParallelUtility::getLocalArrayBounds(un[grid],unLocal,I1,I2,I3);

//       if( solveHelmholtz )
//       {
//         OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);

//         if( ok )
//         {      
//            // // ---> Set  u(dt) = u(0)*cos(omega*(dt)) if time-periodic
//            printF("takeFirstStep:solveHelmholtz: setting  u(dt) = u(0)*cos(omega*(dt)) , dt=%20.12e, omega=%g, frequencyArray(0)=%g\n",dt,omega,frequencyArray(0));

//            // Real diff = max(fabs(ucLocal(I1,I2,I3)-vLocal(I1,I2,I3,0)));
//            // printF(" >> diff | uc - v |=%9.2e\n",diff);

//            unLocal(I1,I2,I3) = vLocal(I1,I2,I3,0)*cos( frequencyArray(0)*dt );
//            for( int freq=1; freq<numberOfFrequencies; freq++ )
//            {
//              unLocal(I1,I2,I3) += vLocal(I1,I2,I3,freq)*cos( frequencyArray(freq)*dt );
//            }
//            // unLocal(I1,I2,I3)  = ucLocal(I1,I2,I3) * cos(omega*dt);        
//         }
//       }
//       else
//       {
//         OV_GET_SERIAL_ARRAY(real,uc[grid],ucLocal);

//         if( ok )
//         {      
//            // // ---> Set  u(dt) = u(0)*cos(omega*(dt)) if time-periodic
//            printF("takeFirstStep: setting  u(dt) = u(0)*cos(omega*(dt)) , dt=%20.12e, omega=%g=%g\n",dt,omega);

//            unLocal(I1,I2,I3)  = ucLocal(I1,I2,I3) * cos(omega*dt);
//         }
//       }

//     }
//     else
//     {
//       OV_GET_SERIAL_ARRAY(real,uc[grid],ucLocal);
//       OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);

//       getIndex(mg.gridIndexRange(),I1,I2,I3);
//       bool ok=ParallelUtility::getLocalArrayBounds(un[grid],unLocal,I1,I2,I3);
//       if( ok )
//       {
      
//         // Taylor series first step **CHECK ME**
//         RealArray lap(I1,I2,I3);
//         operators[grid].derivative(MappedGridOperators::laplacianOperator,ucLocal,lap,I1,I2,I3);
   
//         // -- take a FORWARD STEP ---
//         // u(t-dt) = u(t) + dt*ut + (dt^2/2)*utt + (dt^3/6)*uttt + (dt^4/4!)*utttt
//         //  utt = c^2*Delta(u) + f
//         //  uttt = c^2*Delta(ut) + ft 
//         //  utttt = c^2*Delta(utt) + ftt
//         //        = (c^2*Delta)^2 u + c^2*Delta(f) + ftt 
//         unLocal(I1,I2,I3)  = ucLocal(I1,I2,I3) + dt*(0.) + (.5*dt*dt*c*c)*lap(I1,I2,I3) + (.5*dt*dt *cos(omega*t)*fSign)*fLocal(I1,I2,I3);
//         if( orderOfAccuracy==4 )
//         {
//           // this may be good enough for 4th-order -- local error is dt^4
//           real DeltaUt =0.; // we assume ut=0
//           // upLocal(I1,I2,I3) += ( (dt*dt*dt/6.)*(-omega*sin(omega*t) )*( (c*c)*DeltaUt + fLocal(I1,I2,I3) )
//           unLocal(I1,I2,I3) += ( fSign*(dt*dt*dt/6.)*(-omega*sin(omega*t) ) )*( fLocal(I1,I2,I3) );
//         }

//       }
      
//     }
//   } // end for grid

//   if( !usePeriodicFirstStep )
//   {
//     // This will not work for implicit since only explicit conditions are done here:
//     if( timeSteppingMethod == implicitTimeStepping )
//     {
//       printF("takeFirstStep: ERROR: fix applyBC for implicit time-stepping\n");
//       OV_ABORT("error");
//     }
//     applyBoundaryConditions( u[next], t+dt );
//   }
  
//   // u[next].display(sPrintF("u[next] after first step, t=%9.3e",t+dt),"%6.2f ");


//   return 0;

// }



// //=================================================================================================
// /// \brief Take the first BACKWARD step using Taylor series in time (e.g. for Helmholtz solve) 
// /// THIS ASSUMES A HELMHOLTZ SOLVE OR HELMHOLTZ FORCING 
// //=================================================================================================
// int CgWave::
// takeFirstBackwardStep( int cur, real t )
// {

//   const int & debug           = dbase.get<int>("debug");
//   if( debug & 4 )
//     printF("*******  CgWave::takeFirstBackwardStep GET SOLUTION at -dt *************\n");
  
//   const real & c              = dbase.get<real>("c");
//   const real & dt             = dbase.get<real>("dt");
//   const real & omega          = dbase.get<real>("omega");
//   const int & orderOfAccuracy = dbase.get<int>("orderOfAccuracy");

//   ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");

//   // Do Helmholtz case for now: 
//   assert( forcingOption==helmholtzForcing );

//     //  ---- NOTE: change sign of forcing for Helmholtz since we want to solve ----
//   //       ( omega^2 I + c^2 Delta) w = f  
//   const Real fSign = forcingOption==helmholtzForcing ? -1.0 : 1.0;

//   // const int & solveHelmholtz = dbase.get<int>("solveHelmholtz");
//   // const Real fSign = solveHelmholtz ? -1.0 : 1.0;  
  

//   const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");    
//   const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
//   // const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

//   realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
//   realCompositeGridFunction & un = u[cur];     // current time 
//   realCompositeGridFunction & up = u[prev];    // previous time

//   // forcing: 
//   realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f");

//   CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");

//   Index I1,I2,I3;
//   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
//   {
//     MappedGrid & mg = cg[grid];
//     OV_GET_SERIAL_ARRAY(real,un[grid],unLocal);
//     OV_GET_SERIAL_ARRAY(real,up[grid],upLocal);
//     OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);

//     getIndex(mg.gridIndexRange(),I1,I2,I3);
//     bool ok=ParallelUtility::getLocalArrayBounds(un[grid],unLocal,I1,I2,I3);
//     if( ok )
//     {
//       bool usePeriodicFirstStep=false;
//       if( usePeriodicFirstStep )
//       {
//          // **TESTING** u(-dt) = u(0)*cos(omega*(-dt)) if time-periodic
//         printF("takeFirstBackwardStep: setting  u(-dt) = u(0)*cos(omega*(-dt)) , dt=%20.12e\n",dt);
//         upLocal(I1,I2,I3)  = unLocal(I1,I2,I3) * cos(-omega*dt);
//       }
//       else
//       {
//         RealArray lap(I1,I2,I3);
//         operators[grid].derivative(MappedGridOperators::laplacianOperator,unLocal,lap,I1,I2,I3);
   
//         // -- take a BACKWARD STEP ---
//         // u(t-dt) = u(t) - dt*ut + (dt^2/2)*utt - (dt^3/6)*uttt + (dt^4/4!)*utttt
//         //  utt = c^2*Delta(u) + f
//         //  uttt = c^2*Delta(ut) + ft 
//         //  utttt = c^2*Delta(utt) + ftt
//         //        = (c^2*Delta)^2 u + c^2*Delta(f) + ftt 
//         upLocal(I1,I2,I3)  = unLocal(I1,I2,I3) -dt*(0.) + (.5*dt*dt*c*c)*lap(I1,I2,I3) + (.5*dt*dt *cos(omega*t)*fSign)*fLocal(I1,I2,I3);
//         if( orderOfAccuracy==4 )
//         {
//           // this may be good enough for 4th-order -- local erro is dt^4
//           real DeltaUt =0.; // we assume ut=0
//           // upLocal(I1,I2,I3) += ( -(dt*dt*dt/6.)*(-omega*sin(omega*t) )*( (c*c)*DeltaUt + fLocal(I1,I2,I3) )
//           upLocal(I1,I2,I3) += ( -fSign*(dt*dt*dt/6.)*(-omega*sin(omega*t) ) )*( fLocal(I1,I2,I3) );
//         }

//       }
      
//     }
//   } // end for grid

//   applyBoundaryConditions( u[prev], t-dt );
  

//   return 0;

// }


int CgWave::
outputResults( int current, real t )
// ===================================================================================
// /Description:
//     Save any results after time intervals of tPlot
// 
// ===================================================================================
{
  // *************************************************************
  // write to the check file for regression and convergence tests
  // *************************************************************

  FILE *& checkFile = dbase.get<FILE*>("checkFile");
  assert( checkFile != NULL );

  const real & maxError          = dbase.get<real>("maxError");      // save max-error here 
  const real & solutionNorm      = dbase.get<real>("solutionNorm");  // save solution norm here 
  const int & computeErrors      = dbase.get<int>("computeErrors");
  const int & saveMaxErrors      = dbase.get<int>("saveMaxErrors");   
  
  const int & numberOfComponents= 1;

  const int & computeEnergy  = dbase.get<int>("computeEnergy"); 
  Real energyNorm=0.; 
  if( computeEnergy ) 
    energyNorm = getEnergyNorm( current,t  );

  int numberToOutput =numberOfComponents;

  fPrintF(checkFile,"%9.2e %i  ",t,numberToOutput);
  for( int i=0; i<numberToOutput; i++ )
  {
    real err = maxError;  // maxError=0. if we do not compute errors
    real uc = solutionNorm;
    fPrintF(checkFile,"%i %9.2e %10.3e  ",i,err,uc);
  }
  fPrintF(checkFile,"\n");

  const int & numberOfSequences = dbase.get<int>("numberOfSequences");
  RealArray solutionNormVector(numberOfSequences);
  solutionNormVector(0)=solutionNorm;
  if( computeEnergy )
    solutionNormVector(computeEnergy)=energyNorm;
  if( saveMaxErrors )
  {

    // printf("\n &&&&&&&&&& saveMaxErrors=%d at t=%g numberOfSequences=%d maxError=%10.2e\n\n",saveMaxErrors,t,numberOfSequences,maxError);
    solutionNormVector(saveMaxErrors)=maxError;
  }

  saveSequenceInfo( t, solutionNormVector );

  return 0;
}

// =======================================================================
/// "Sort" the array of eigenvalues. Return the permuation array iperm.
// =======================================================================
int CgWave::
sortArray( RealArray & eigenValues, IntegerArray & iperm )
{
  int numberOfEigenvectors = eigenValues.getLength(0);

  for( int ie=0; ie<numberOfEigenvectors; ie++ )
  {
    iperm(ie)=ie; // for sorting 
  }
  // -- bubble sort ---
  for( int ie=0; ie<numberOfEigenvectors-1; ie++ )
  {
    bool changed=false;
    for( int je=0; je<numberOfEigenvectors-1; je++ )
    {
      if( eigenValues(je) > eigenValues(je+1) )  // increasing order
      {
        changed=true;
        Real temp=eigenValues(je);  eigenValues(je)= eigenValues(je+1);  eigenValues(je+1)=temp;
        int itemp=      iperm(je);        iperm(je)=       iperm(je+1);        iperm(je+1)=itemp;
      }
    }
    if( !changed ) break;
  }

  return 0;
}


// =============================================================================
/// \brief : return the coefficient of upwind dissipation for a given grid
/// \param dtUpwind (input) : if positive, use this value for dt 
// =============================================================================
Real CgWave::
getUpwindDissipationCoefficient( int grid, Real dtUpwind /* = -1. */, bool adjustForTimeStep /* = true */    )
{
  Real & c                          = dbase.get<real>("c");
  const int & orderOfAccuracy       = dbase.get<int>("orderOfAccuracy");
  const int & upwind                = dbase.get<int>("upwind");
  const int & implicitUpwind        = dbase.get<int>("implicitUpwind"); 

  const int & useSuperGrid          = dbase.get<int>("useSuperGrid"); 

  RealArray & gridCFL               = dbase.get<RealArray>("gridCFL");
  const int & solveHelmholtz        = dbase.get<int>("solveHelmholtz");
  const int & adjustOmega           = dbase.get<int>("adjustOmega");
  RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
  RealArray & frequencyArrayAdjusted= dbase.get<RealArray>("frequencyArrayAdjusted");
  RealArray & frequencyArraySave    = dbase.get<RealArray>("frequencyArraySave");
  Real & dt                         = dbase.get<real>("dt");
  Real & dtUsed                     = dbase.get<real>("dtUsed");  // dt actually used

  const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");
  IntegerArray & gridIsImplicit     = dbase.get<IntegerArray>("gridIsImplicit");
  const RealArray & bImp            = dbase.get<RealArray>("bImp");
  const RealArray & cImp            = dbase.get<RealArray>("cImp");


  Real upwindCoefficient = 1.;

  const int numberOfDimensions = cg.numberOfDimensions();

  if( dtUpwind <=0. )
    dtUpwind = dt; 

  // Default coeff for explicit UPWIND Predictor corrector
  if( timeSteppingMethod==explicitTimeStepping )
  {
    // --- explicit time-stepping + Upwind predictor-corrector ---
    // NOTE: For stability we must limit the upwind coefficient : see AMP/ssmx/ssmx.pdf 

    upwindCoefficient = (c*dtUpwind)/( sqrt(1.*numberOfDimensions) * pow(2.,orderOfAccuracy+1) );
    // printF("getUpwindDissipationCoefficient: c=%g, dt=%g, upwindCoefficient=%g\n",c,dt,upwindCoefficient);
  }
  else
  {
    // --- implicit time-stepping ---

    // **FIX ME: use Allison's new formula ...
    upwindCoefficient = (c*dtUpwind)/( sqrt(1.*numberOfDimensions) * pow(2.,orderOfAccuracy+1) );
    if( false && debug & 1 )
      printF("$$$$ getUpwindDissipationCoefficient: Set upwindCoefficient = (c*dtUpwind)/( sqrt(1.*numberOfDimensions) * pow(2.,orderOfAccuracy+1) ) = %10.2e\n",upwindCoefficient); 

    if( upwind && !implicitUpwind && adjustForTimeStep )
    {
      // --- implicit time-stepping + explicit upwind PC ---

      // NOTE: For stability we must limit the upwind coefficient : see notes in AMP/wimp/wimp.pdf

      // **CHECK ME**

      Real adjustmentFactor=1.;
      
      assert( dt>0. );

      Real myGridCFL = dt*gridCFL(grid);   // gridCFL is really c/dx , i.e. does not have the factor of "dt"

      adjustmentFactor = 1./myGridCFL;

      upwindCoefficient *= adjustmentFactor;

      if( false && debug & 1 )
        printF("$$$$ getUpwindDissipationCoefficient: Adjust useUpwindDissipationCoeff by adjustmentFactor=%9.2e, 1/gridCFL = %9.2e\n\n",
         adjustmentFactor, 1./myGridCFL);
    }

    if( false && upwind && implicitUpwind && adjustForTimeStep ) 
    {
      // **TESTING:  implicit Upwind :  increase upwind with CFL 
      
      assert( dt>0. );

      Real myGridCFL = dt*gridCFL(grid);   // gridCFL is really c/dx , i.e. does not have the factor of "dt"
      // Real adjustmentFactor = max( 1., myGridCFL );

      Real adjustmentFactor = 1.; // ################# TEST ##############

      upwindCoefficient *= adjustmentFactor;      
    } 
     
  }

  
  if( adjustForTimeStep )
  {
    // -- make adjustments to correct for time-discretization errors if we are trying to match to a direct Helmholtz solver 

    if( timeSteppingMethod==implicitTimeStepping )
    {
      if( solveHelmholtz )
      // if( upwind && implicitUpwind && solveHelmholtz )
      {
        // -- adjust for implicit time-stepping and implicitUpwind (upwinding included in implicit matrix) --

        Real omega = frequencyArray(0);
        if( adjustOmega )
        {
          omega = frequencyArrayAdjusted(0);

          if( frequencyArrayAdjusted(0)==0. )
          {
            printF("\n XXXX getUpwindDissipationCoefficient:ERROR: adjustOmega *but* frequencyArrayAdjusted==0 ! XXXX\n\n");
            OV_ABORT("etUpwindDissipationCoefficient:ERROR: fix me");
          }

        }
        if( dtUsed>0. )
        {

          Real adjustmentFactor = (frequencyArraySave(0)*dt)/tan(omega*dt);
          // Real adjustmentFactor = (frequencyArraySave(0)*dt)/tan(frequencyArrayAdjusted(0)*dt);

          if( false )
            printF("\n $$$$ getUpwindDissipationCoefficient: frequencyArray(0)=%14.6e, frequencyArrayAdjusted(0)=%14.6e, frequencyArraySave(0)=%14.6e, dt=%9.2e, dtUsed=%e, adjustmentFactor=%12.4e $$$$\n\n",
              frequencyArray(0),frequencyArrayAdjusted(0),frequencyArraySave(0),dt,dtUsed,adjustmentFactor);
          upwindCoefficient *= adjustmentFactor;
        }
        else
        {
          printF("\n $$$$ WARNING: getUpwindDissipationCoefficient: no adjustment made to upwind dissipation coefficient"
                 " since dt=%e has not been set yet $$$$\n\n",dt);
        }
      }

    }
  }

  // if( par.implicitUpwind && par.upwind )
  //    if( par.useWaveHoltz )
  //       betaAdjusted = (par.frequencyArraySaved(1)*par.dt)/tan(par.frequencyArrayAdjusted(1)*par.dt);
  //     else
  //       betaAdjusted=1;
  //     end


  // Implicit scheme + explicit UPC -- limit coeff by the gridCFL

  //  adSosup = adSosup/gridCFL

  // Adjust for solveHelmholtz 

  // Implicit monolithic scheme 

  return upwindCoefficient;

}

// =============================================================================
/// \brief : Check parameters for consistency.
/// \notes: Add more checks here to avoid user (and my) mistakes 
/// \notes: WDH : started Oct 3, 2024.
// =============================================================================
int CgWave::checkParameters()
{

  const real & c                        = dbase.get<real>("c");
  const real & dt                       = dbase.get<real>("dt");
  const real & omega                    = dbase.get<real>("omega");
  const int & orderOfAccuracy           = dbase.get<int>("orderOfAccuracy");
  const int & takeImplicitFirstStep     = dbase.get<int>("takeImplicitFirstStep");
  const int & checkParameters           = dbase.get<int>("checkParameters");  // 1= check problem parameters for consistency


  const int & numberOfFrequencies       = dbase.get<int>("numberOfFrequencies");
  const RealArray & frequencyArray      = dbase.get<RealArray>("frequencyArray");  
  const RealArray & frequencyArraySave  = dbase.get<RealArray>("frequencyArraySave");  

  const int & solveHelmholtz            = dbase.get<int>("solveHelmholtz");
  const int & filterTimeDerivative      = dbase.get<int>("filterTimeDerivative");

  ForcingOptionEnum & forcingOption     = dbase.get<ForcingOptionEnum>("forcingOption");
  const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");  

  int numErrors=0; 
  // if( solveHelmholtz && takeImplicitFirstStep && filterTimeDerivative && timeSteppingMethod == implicitTimeStepping )
  // {
  //   printF("CgWave::checkParameters: ERROR: takeImplicitFirstStep should not be used with solveHelmholtz and filterTimeDerivative.\n");
  //   numErrors++;
  // }

  if( solveHelmholtz && filterTimeDerivative && numberOfFrequencies>1  )
  {
    printF("CgWave::checkParameters: ERROR: numberOfFrequencies>1 cannot yet be used with solveHelmholtz and filterTimeDerivative.\n");
    numErrors++;
  }

  if( checkParameters && numErrors>0 )
  {
    printF("CgWave::checkParameters: Fix errors (or set checkParameters=0 to continue anyway)\n");
    OV_ABORT("CgWave::checkParameters: Fix errors (or set checkParameters=0 to continue anyway)")
  }

  return 0;

}

// // =============================================================================
// /// \brief : Return Helmholtz related parameter values that have been adjusted for 
// ///    time-discretization errors
// // =============================================================================
// int CgWave::
// getHelmholtzAdjustedParameters( Real 
// {
//   Real & c                          = dbase.get<real>("c");
//   const int & orderOfAccuracy       = dbase.get<int>("orderOfAccuracy");
//   const int & upwind                = dbase.get<int>("upwind");
//   const int & implicitUpwind        = dbase.get<int>("implicitUpwind");  

//   RealArray & gridCFL               = dbase.get<RealArray>("gridCFL");
//   const int & solveHelmholtz        = dbase.get<int>("solveHelmholtz");
//   RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
//   RealArray & frequencyArrayAdjusted= dbase.get<RealArray>("frequencyArrayAdjusted");
//   RealArray & frequencyArraySave    = dbase.get<RealArray>("frequencyArraySave");
//   Real & dt                         = dbase.get<real>("dt");
//   Real & dtUsed                     = dbase.get<real>("dtUsed");  // dt actually used

//   const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");
//   IntegerArray & gridIsImplicit     = dbase.get<IntegerArray>("gridIsImplicit");
//   const RealArray & bImp            = dbase.get<RealArray>("bImp");
//   const RealArray & cImp            = dbase.get<RealArray>("cImp");


//   Real upwindCoefficient = 1.;

//   const int numberOfDimensions = cg.numberOfDimensions();
