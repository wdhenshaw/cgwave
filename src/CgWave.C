// ====================== Composite Grid Wave Equation Solver Class =====================

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

  dbase.put<real>("c")=1.;
  dbase.put<real>("cfl")=.7; // .25;
  dbase.put<real>("tFinal")=1.;
  dbase.put<real>("tPlot")=.1;
  dbase.put<real>("dt")=-1.;
  dbase.put<Real>("dtMax")=-1.;  // save maximum allowable dt (before corrections to reach a given time)

  dbase.put<real>("ad4")=0.; // coeff of the artificial dissipation.
  dbase.put<int>("dissipationFrequency")=1; // apply dissipation every this many steps (1= every step)

  dbase.put<int>("current")=0;              // holds current solution

  dbase.put<int>("orderOfAccuracy")=2;
  dbase.put<int>("orderOfAccuracyInTime")=-1; // -1 means make order in time equal to order in space

  dbase.put<int>("secondOrderGrid")=1;
  dbase.put<int>("debug")=0;

  dbase.put<TimeSteppingMethodEnum>("timeSteppingMethod")=explicitTimeStepping;

  // coefficients in implicit time-stepping  
  //  D+t D-t u = c^2 Delta( cImp(1) *u^{n+1} + cImp(0) *u^n + cImp(-1)* u^{n-1} )
  RealArray & cImp = dbase.put<RealArray>("cImp");
  cImp.redim(Range(-1,1)); 
  // Full-weighting by default: 
  cImp(-1)=.25;
  cImp( 0)=.5;
  cImp( 1)=.25;

  // interactiveMode :
  //       0 = plot intermediate results
  //       1 = run without plotting and exit advance when finished 
  dbase.put<int>("interactiveMode")=0;


  dbase.put<real>("numberOfGridPoints")=0.;
  dbase.put<int>("numberOfStepsTaken")=0;
  
  dbase.put<RealArray>("dxMinMax");
  dbase.put<aString>("nameOfGridFile")="unknown";

  dbase.put<real>("maxError")=0.;      // save max-error here 
  dbase.put<real>("solutionNorm")=1.;  // save solution norm here 
  dbase.put<int>("computeErrors")=0;   // true of we compute errors 

  // For Helmholtz solve with CgWaveHoltz
  dbase.put<int>("solveHelmholtz")=false;
  dbase.put<int>("computeTimeIntegral")=false; // if true compute the time-integral (for the Helmholtz solve or other reason)

  dbase.put<int>("adjustOmega")=0;                 // 1 : choose omega from the symbol of D+t D-t 
  real & omega = dbase.put<real>("omega")=30.1;
  dbase.put<real>("Tperiod")=twoPi/omega;
  dbase.put<int>("numPeriods")=10;

  dbase.put<Real>("omegaSave")   = -1.;  // used when we adjust omega 
  dbase.put<Real>("TperiodSave") = -1.;
  
  dbase.put<real>("tol")=1.e-4;  // tolerance for Krylov solvers 

  real & omegaSOR = dbase.put<real>("omegaSOR")=1.;

  // Gaussian forcing: 
  dbase.put<real>("beta")=100.;
  dbase.put<real>("x0")=0.;
  dbase.put<real>("y0")=0.;
  dbase.put<real>("z0")=0.;

  int & numberOfTimeLevelsStored = dbase.put<int>("numberOfTimeLevelsStored")=3;

  realCompositeGridFunction *& ucg = dbase.put<realCompositeGridFunction*>("ucg");
  ucg = new realCompositeGridFunction[numberOfTimeLevelsStored];
 
  dbase.put<realCompositeGridFunction>("f");  // source term 

  dbase.put<CompositeGridOperators>("operators");
  dbase.put<Interpolant*>("pInterpolant")=NULL;

  dbase.put<realCompositeGridFunction>("v");

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

  dbase.put<realCompositeGridFunction>("error");

  dbase.put<GUIState*>("runTimeDialog")=NULL;

  dbase.put<int>("movieFrame")=-1;
  dbase.put<int>("plotOptions")=1;
  dbase.put<int>("plotChoices")=0;

  dbase.put<int>("myid")=max(0,Communication_Manager::My_Process_Number);
  int & np = dbase.put<int>("np");
  #ifdef USE_PPP
    np= max(1,Communication_Manager::numberOfProcessors());
  #else
    np=1;
  #endif
  
  dbase.put<PlotStuffParameters>("psp");

  dbase.put<aString>("movieFileName")="cgWave";

  dbase.put<aString>("knownSolutionOption")="noKnownSolution";
  dbase.put<bool>("knownSolutionIsTimeDependent")=false; 

  // These should match BoundaryConditionEnum
  int & numberOfBCNames = dbase.put<int>("numberOfBCNames")=5;
  aString *& bcNames = dbase.put<aString*>("bcNames") = new aString [numberOfBCNames];
  bcNames[0]="periodic";
  bcNames[1]="dirichlet";
  bcNames[2]="neumann";
  bcNames[3]="evenSymmetry";
  bcNames[4]="radiation";
  

  FILE *& debugFile = dbase.put<FILE*>("debugFile");
  debugFile = fopen("cgWave.debug","w" );        // log file 

  FILE *& logFile = dbase.put<FILE*>("logFile");
  logFile = fopen("cgWave.out","w" );        // log file 

  FILE *& checkFile = dbase.put<FILE*>("checkFile");
  checkFile = fopen("cgWave.check","w" );        // for regression and convergence tests

  // enum TimingEnum
  // { 
  //   totalTime=0,
  //   timeForInitialize,
  //   timeForInitialConditions,
  //   timeForAdvance,
  //   timeForAdvanceRectangularGrids,
  //   timeForAdvanceCurvilinearGrids,
  //   timeForAdvOpt,
  //   timeForDissipation,
  //   timeForBoundaryConditions,
  //   timeForInterpolate,
  //   timeForUpdateGhostBoundaries,
  //   timeForForcing,
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
  timingName[timeForInitialConditions]           ="initial conditions";
  timingName[timeForAdvance]                     ="advance";
  timingName[timeForAdvanceRectangularGrids]     ="  advance rectangular grids";
  timingName[timeForAdvanceCurvilinearGrids]     ="  advance curvilinear grids";
  timingName[timeForForcing]                     ="  add forcing";
  timingName[timeForImplicitSolve]               ="    implicit solve";
  timingName[timeForDissipation]                 ="  add dissipation";
  timingName[timeForBoundaryConditions]          ="  boundary conditions";
  timingName[timeForUpdateGhostBoundaries]       ="  update ghost (parallel)";
  timingName[timeForInterpolate]                 ="  interpolation";
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
  
  printF("CgWave::initialize : assign forcing...\n");
  const int & addForcing = dbase.get<int>("addForcing");
  const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");

  const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");
  bool twilightZone = addForcing && forcingOption==twilightZoneForcing;
  int & computeErrors = dbase.get<int>("computeErrors");
  computeErrors = twilightZone || knownSolutionOption=="userDefinedKnownSolution";
  
  const real & omega     = dbase.get<real>("omega");
  const real & cfl       = dbase.get<real>("cfl");
  const real & dt        = dbase.get<real>("dt");


  const int & orderOfAccuracy = dbase.get<int>("orderOfAccuracy");

  int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");
  // By default, order in time = order in space : 
  if( orderOfAccuracyInTime==-1 )
    orderOfAccuracyInTime = orderOfAccuracy;

  const real & ad4            = dbase.get<real>("ad4"); // coeff of the artificial dissipation.

  bool useUpwindDissipation = ad4  > 0.;
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
        printF("--CgWave-- setupGridFunctions: ERROR: the grid does not have enough ghost points for upwind dissipation.\n"
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

  // ------------ Helmholtz -------------
  const int & solveHelmholtz = dbase.get<int>("solveHelmholtz");
  int & computeTimeIntegral = dbase.get<int>("computeTimeIntegral");
  if( solveHelmholtz )
    computeTimeIntegral=true;

  if( computeTimeIntegral )
  {
    // allocate grid functions for time integral and Helmholtz solver ----
    CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");
    Range all;

    realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");
    v.updateToMatchGrid(cg,all,all,all);
    v.setOperators(operators); 
    v=0.;

    // realCompositeGridFunction & vOld = dbase.get<realCompositeGridFunction>("vOld");
    // vOld.updateToMatchGrid(cg,all,all,all);

    // realCompositeGridFunction & residual = dbase.get<realCompositeGridFunction>("residual");
    // residual.updateToMatchGrid(cg,all,all,all);
  }
  
  
  realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f");
  if( forcingOption != noForcing )
    f.updateToMatchGrid(cg);
  
  Index I1,I2,I3;

  if( forcingOption==helmholtzForcing )
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
      ipar[0]=grid;
      userDefinedForcing( f[grid], ipar,rpar );
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

      if( degreeInSpace>6 )
      {
        printF("CgWave:ERROR: degreeInSpace=%d is not supported by OGPolyFunction, Reducing to 6\n",degreeInSpace);
        degreeInSpace=6;
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


  // -- compute the time-step ---
  getTimeStep(); 
  
  printF("CgWave::initialize: dt=%g\n",dt);

  Real & dtMax = dbase.get<Real>("dtMax"); 
  dtMax = dt; // save the initial timeStep

  timing(timeForInitialize) += getCPU()-cpu0;
  
  return 0;
}



// ================================================================================================
/// \brief Assign parameters 
// ================================================================================================
int CgWave::interactiveUpdate()
{
//  GenericGraphicsInterface & gi = *Overture::getGraphicsInterface();
  PlotStuffParameters psp;

  real & cfl                    = dbase.get<real>("cfl");
  real & tFinal                 = dbase.get<real>("tFinal");
  real & tPlot                  = dbase.get<real>("tPlot");
  Real & dtMax                  = dbase.get<Real>("dtMax"); 
  real & omega                  = dbase.get<real>("omega");
  real & Tperiod                = dbase.get<real>("Tperiod");
  int & numPeriods              = dbase.get<int>("numPeriods");
  int & orderOfAccuracy         = dbase.get<int>("orderOfAccuracy");
  int & orderOfAccuracyInTime   = dbase.get<int>("orderOfAccuracyInTime");
  IntegerArray & gridIsImplicit = dbase.get<IntegerArray>("gridIsImplicit");
  RealArray & cImp              = dbase.get<RealArray>("cImp");


  real & ad4                    = dbase.get<real>("ad4"); // coeff of the artificial dissipation.
  int & dissipationFrequency    = dbase.get<int>("dissipationFrequency");

  int & debug                   = dbase.get<int>("debug");
  int & interactiveMode         = dbase.get<int>("interactiveMode");
  
  int & solveHelmholtz          = dbase.get<int>("solveHelmholtz");
  int & computeTimeIntegral     = dbase.get<int>("computeTimeIntegral");
  real & tol                    = dbase.get<real>("tol");
          
  real & omegaSOR               = dbase.get<real>("omegaSOR");

  // Gaussian forcing 
  real & beta = dbase.get<real>("beta");
  real & x0   = dbase.get<real>("x0");
  real & y0   = dbase.get<real>("y0");
  real & z0   = dbase.get<real>("z0");

  TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

  TwilightZoneEnum & twilightZone = dbase.get<TwilightZoneEnum>("twilightZone");
  int & degreeInSpace             = dbase.get<int>("degreeInSpace");
  int & degreeInTime              = dbase.get<int>("degreeInTime");
  RealArray & trigFreq            = dbase.get<RealArray>("trigFreq");

  int & addForcing = dbase.get<int>("addForcing");
  ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");

  InitialConditionOptionEnum & initialConditionOption = dbase.get<InitialConditionOptionEnum>("initialConditionOption");

  const int & numberOfBCNames = dbase.get<int>("numberOfBCNames");
  aString *& bcNames = dbase.get<aString*>("bcNames");


  // Build a dialog menu for changing parameters
  GUIState gui;
  DialogData & dialog=gui;

  dialog.setWindowTitle("CgWave - Wave Equation Solver");
  dialog.setExitCommand("exit", "exit");

  dialog.setOptionMenuColumns(1);

  aString timeSteppingLabel[] = {"explicit", "implicit", "" };
  dialog.addOptionMenu("time stepping:", timeSteppingLabel, timeSteppingLabel, (int)timeSteppingMethod );

// enum InitialConditionOptionEnum
// {
//   zeroInitialCondition=0,
//   twilightZoneInitialCondition,
//   knownSolutionInitialCondition,
//   pulseInitialCondition
// };

  aString initialConditionLabel[] = {"zero initial condition", "twilightZone initial condition", "known solution initial condition", "pulse initial condition", "" };
  dialog.addOptionMenu("Initial conndtions:", initialConditionLabel, initialConditionLabel, (int)initialConditionOption );

  aString forcingLabel[] = {"no forcing", "twilightZoneForcing", "userForcing", "helmholtzForcing", "" };
  dialog.addOptionMenu("forcing:", forcingLabel, forcingLabel, (int)forcingOption );

  aString tzLabel[] = {"polynomial", "trigonometric", "" };
  dialog.addOptionMenu("Twilight zone::",tzLabel,tzLabel,(int)twilightZone );

  aString pbLabels[] = {
                        "user defined known solution...",
                        "user defined forcing...",
                        "choose grids for implicit",
                        "grid",
                        "erase",
                        "exit",
                        ""};
  int numRows=5;
  dialog.setPushButtons( pbLabels, pbLabels, numRows ); 

  aString tbCommands[] = {"turn on forcing",
                          "solve Helmholtz",
                          "compute time integral",
                            ""};
  int tbState[10];
  tbState[0] = addForcing;
  int numColumns=1;
  dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 

  // ----- Text strings ------
  const int numberOfTextStrings=30;
  aString textCommands[numberOfTextStrings];
  aString textLabels[numberOfTextStrings];
  aString textStrings[numberOfTextStrings];

  int nt=0;
  textCommands[nt] = "tFinal";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",tFinal);  nt++; 

  textCommands[nt] = "cfl";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",cfl);  nt++; 

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

  textCommands[nt] = "number of periods";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numPeriods);  nt++; 

  textCommands[nt] = "artificial dissipation";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",ad4);  nt++; 

  textCommands[nt] = "Gaussian params";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g %g %g %g (beta,x0,y0,z0)",beta,x0,y0,z0);  nt++; 

  textCommands[nt] = "tol";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",tol);  nt++; 

  textCommands[nt] = "debug";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",debug);  nt++; 

  textCommands[nt] = "interactiveMode";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",interactiveMode);  nt++; 

  textCommands[nt] = "bc=";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "bcName ([dirichlet|evenSymmetry]");  nt++; 

  textCommands[nt] = "dissipationFrequency";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",dissipationFrequency);  nt++; 

  textCommands[nt] = "dtMax";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",dtMax);  nt++; 

  textCommands[nt] = "implicit weights";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g, %g, %g",cImp(-1),cImp(0),cImp(1));  nt++; 
  
  // null strings terminal list
  textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
  dialog.setTextBoxes(textCommands, textLabels, textStrings);

  

  gi.pushGUI(gui);

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

    // else if( dialog.getTextValue(answer,"omegaSOR","%e",omegaSOR) )
    // {
    //   printF("Setting omegaSOR=%g\n",omegaSOR);
    // }
    
    // else if( dialog.getTextValue(answer,"omega","%e",omega) )
    // {
    //   printF("Setting omega=%g\n",omega);
    // }

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


    else if( dialog.getTextValue(answer,"tol","%e",tol) )
    {
      printF("Setting tol=%g (tolerence for Krylov solvers)\n",tol);
    }

    else if( dialog.getTextValue(answer,"dtMax","%e",dtMax) )
    {
      printF("Setting dtMax=%g\n",dtMax);
    }    
    
    else if( dialog.getToggleValue(answer,"turn on forcing",addForcing) ){}//
    else if( dialog.getToggleValue(answer,"solve Helmholtz",solveHelmholtz) )
    {
      if( solveHelmholtz )
      {
        printF(" solveHelmholtz=true: CgWave is being used to solve the Helmholtz equation (time-periodic wave equation)\n"
               " using CgWaveHoltz\n");
      }
      
    }
    else if( dialog.getToggleValue(answer,"compute time integral",computeTimeIntegral) )
    {
      printF("Setting computeTimeIntegral=%i\n",computeTimeIntegral);
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
      sScanF(answer(len,answer.length()-1),"%e %e %e",&cImp(-1),&cImp(0),&cImp(1));
      printF("Setting implicit time-stepping weights to cImp(-1)=%g, cImp(0)=%g, cImp(1)=%g\n",cImp(-1),cImp(0),cImp(1));
    } 

    // aString initialConditionLabel[] = {"zero initial condition", "twilightZone initial condition", "known solution initial condition", "pulse initial condition", "" };
    else if( answer=="zero initial condition"           || 
             answer=="twilightZone initial condition"   || 
             answer=="known solution initial condition" || 
             answer=="pulse initial condition" )
    {
      initialConditionOption = ( answer=="zero initial condition"           ? zeroInitialCondition :
                                 answer=="twilightZone initial condition"   ? twilightZoneInitialCondition :
                                 answer=="known solution initial condition" ? knownSolutionInitialCondition :
                                 answer=="pulse initial condition"          ? pulseInitialCondition : 
                                                                              zeroInitialCondition );
    }


    else if( answer=="noForcing" || answer=="twilightZoneForcing" || answer=="userForcing" || answer=="helmholtzForcing" )
    {
      forcingOption = ( answer=="noForcing" ? noForcing :
                        answer=="twilightZoneForcing" ? twilightZoneForcing :
                        answer=="userForcing" ? userForcing :
                        answer=="helmholtzForcing" ? helmholtzForcing : noForcing);

      if( forcingOption==twilightZoneForcing )
        initialConditionOption = twilightZoneInitialCondition;
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

    else if( len=answer.matches("tPlot") )
    {
      sScanF(answer(len,answer.length()-1),"%e",&tPlot);
      printF(" tPlot=%g\n",tPlot);
      dialog.setTextLabel("tPlot",sPrintF(line, "%g",tPlot));
    }
    else if( len=answer.matches("artificial dissipation") )
    {
      sScanF(answer(len,answer.length()-1),"%e",&ad4);
      printF(" artificial diffusion=%g\n",ad4);
      dialog.setTextLabel("artificial dissipation",sPrintF(line, "%g",ad4));
    }
    else if( len=answer.matches("Gaussian params") )
    {
      printF("Gaussian forcing: \n");
      sScanF(answer(len,answer.length()-1),"%e %e %e %e",&beta,&x0,&y0,&z0);
      dialog.setTextLabel("Gaussian params",sPrintF(line, "%g %g %g %g (beta,x0,y0,z0)",
                                                    beta,x0,y0,z0));
      printF("Gaussian force: beta=%g x0=%g y0=%g z0=%g \n",beta,x0,y0,z0);
    }
    else if( dialog.getTextValue(answer,"number of periods","%i",numPeriods) )
    {
      printF("Setting numPeriods=%i\n",numPeriods);
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
    ucg[i].setOperators(operators);                                 
    ucg[i].setName("u");                       // name the grid function
  }
  
  // gridIsImplicit(grid) = 1 : this grid is advanced implicitly
  IntegerArray & gridIsImplicit = dbase.get<IntegerArray>("gridIsImplicit");
  gridIsImplicit.redim(cg.numberOfComponentGrids());
  gridIsImplicit=0;

  timing(timeForInitialize) += getCPU()-cpu0;

  return 0;
}






// ======================================================================================================
/// \brief Compute errors
// ======================================================================================================
real CgWave::
getErrors( realCompositeGridFunction & u, real t )
{
  real cpu0=getCPU();
  // printF("+++++++++++ getErrors t=%9.3e +++++++++++\n",t);

  real & maxError     = dbase.get<real>("maxError");      // save max-error here 
  real & solutionNorm = dbase.get<real>("solutionNorm");  // save solution norm here 
  int & computeErrors = dbase.get<int>("computeErrors");
  
  maxError=0.;

  const int & addForcing = dbase.get<int>("addForcing");
  const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
  const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");

  bool twilightZone = addForcing && forcingOption==twilightZoneForcing;

//  computeErrors = twilightZone || knownSolutionOption=="userDefinedKnownSolution";
  

  if( computeErrors )
  {
    realCompositeGridFunction & error = dbase.get<realCompositeGridFunction>("error");
    error.updateToMatchGrid(cg);

    const int numberOfDimensions = cg.numberOfDimensions();
    Index I1,I2,I3;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg = cg[grid];
    
      // get the local serial arrays
      OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
      OV_GET_SERIAL_ARRAY(real,error[grid],errLocal);
      errLocal=0.;

      getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points.

      if( knownSolutionOption=="userDefinedKnownSolution" )
      {
        // printF("+++++++++++ getErrors for userDefinedKnownSolution +++++++++++\n");
        

        getUserDefinedKnownSolution( t, grid, error[grid], I1,I2,I3 ); // store true solution in error[grid]
        const int includeGhost=1;
        bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
        if( ok )
        {
          errLocal(I1,I2,I3) -= uLocal(I1,I2,I3);
        }
      }
      else
      {
        // ----- Twilight zone ------
        assert( dbase.get<OGFunction*>("tz")!=NULL );
        OGFunction & e = *dbase.get<OGFunction*>("tz");

        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
        OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);


        const int includeGhost=1;
        bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
        if( ok )
        {
          int numberOfComponents=1;
          Range C=numberOfComponents;
          int isRectangular=0;
          RealArray ue(I1,I2,I3);
          e.gd( ue ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,C,t);

          errLocal(I1,I2,I3) = ue(I1,I2,I3) - uLocal(I1,I2,I3);

        }
      }
      
    }
    maxError = maxNorm(error);
    // printF("getErrors: t=%9.3e, maxError=%9.3e\n",t,maxError);


  }
  else if( knownSolutionOption=="userDefinedKnownSolution" )
  {
  }
  
  // compute the solution norm
  solutionNorm = maxNorm(u);
  

  timing(timeForGetError)+= getCPU()-cpu0;

  return maxError;
  


}


//=================================================================================================
/// \brief Take the first BACKWARD step using Taylor series in time (e.g. for Helmholtz solve) 
/// THIS ASSUMES A HELMHOLTZ SOLVE OR HELMHOLTZ FORCING 
//=================================================================================================
int CgWave::
takeFirstBackwardStep( int cur, real t )
{

  const int & debug           = dbase.get<int>("debug");
  if( debug & 4 )
    printF("*******  CgWave::takeFirstBackwardStep GET SOLUTION at -dt *************\n");
  
  const real & c              = dbase.get<real>("c");
  const real & dt             = dbase.get<real>("dt");
  const real & omega          = dbase.get<real>("omega");
  const int & orderOfAccuracy = dbase.get<int>("orderOfAccuracy");

  ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");

  // Do Helmholtz case for now: 
  assert( forcingOption==helmholtzForcing );

    //  ---- NOTE: change sign of forcing for Helmholtz since we want to solve ----
  //       ( omega^2 I + c^2 Delta) w = f  
  const Real fSign = forcingOption==helmholtzForcing ? -1.0 : 1.0;

  // const int & solveHelmholtz = dbase.get<int>("solveHelmholtz");
  // const Real fSign = solveHelmholtz ? -1.0 : 1.0;  
  

  const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");    
  const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
  // const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

  realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
  realCompositeGridFunction & un = u[cur];     // current tinme 
  realCompositeGridFunction & up = u[prev];    // previous time

  // forcing: 
  realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f");

  CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");

  Index I1,I2,I3;
  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  {
    MappedGrid & mg = cg[grid];
    OV_GET_SERIAL_ARRAY(real,un[grid],unLocal);
    OV_GET_SERIAL_ARRAY(real,up[grid],upLocal);
    OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);

    getIndex(mg.gridIndexRange(),I1,I2,I3);
    bool ok=ParallelUtility::getLocalArrayBounds(un[grid],unLocal,I1,I2,I3);
    if( ok )
    {
      RealArray lap(I1,I2,I3);
      operators[grid].derivative(MappedGridOperators::laplacianOperator,unLocal,lap,I1,I2,I3);
 
      // -- take a BACKWARD STEP ---
      // u(t-dt) = u(t) - dt*ut + (dt^2/2)*utt - (dt^3/6)*uttt + (dt^4/4!)*utttt
      //  utt = c^2*Delta(u) + f
      //  uttt = c^2*Delta(ut) + ft 
      //  utttt = c^2*Delta(utt) + ftt
      //        = (c^2*Delta)^2 u + c^2*Delta(f) + ftt 
      upLocal(I1,I2,I3)  = unLocal(I1,I2,I3) -dt*(0.) + (.5*dt*dt*c*c)*lap(I1,I2,I3) + (.5*dt*dt *cos(omega*t)*fSign)*fLocal(I1,I2,I3);
      if( orderOfAccuracy==4 )
      {
        // this may be good enough for 4th-order -- local erro is dt^4
        real DeltaUt =0.; // we assume ut=0
        // upLocal(I1,I2,I3) += ( -(dt*dt*dt/6.)*(-omega*sin(omega*t) )*( (c*c)*DeltaUt + fLocal(I1,I2,I3) )
        upLocal(I1,I2,I3) += ( -fSign*(dt*dt*dt/6.)*(-omega*sin(omega*t) ) )*( fLocal(I1,I2,I3) );
      }
      

    }
  } // end for grid

  applyBoundaryConditions( u[prev], t-dt );
  

  return 0;

}


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

  const real & maxError     = dbase.get<real>("maxError");      // save max-error here 
  const real & solutionNorm = dbase.get<real>("solutionNorm");  // save solution norm here 
  const int & computeErrors = dbase.get<int>("computeErrors");
  
  const int & numberOfComponents= 1;

  int numberToOutput =numberOfComponents;

  fPrintF(checkFile,"%9.2e %i  ",t,numberToOutput);
  for( int i=0; i<numberToOutput; i++ )
  {
    real err = maxError;  // maxError=0. if we do not compute errors
    real uc = solutionNorm;
    fPrintF(checkFile,"%i %9.2e %10.3e  ",i,err,uc);
  }
  fPrintF(checkFile,"\n");


  return 0;
}
