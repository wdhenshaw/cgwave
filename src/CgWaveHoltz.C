#include "CgWaveHoltz.h"
#include "CompositeGridOperators.h"	
#include "PlotStuff.h"
#include "display.h"
#include "ParallelOverlappingGridInterpolator.h"
#include "ParallelUtility.h"
#include "LoadBalancer.h"
#include "gridFunctionNorms.h"
#include <sstream>

#include "CgWave.h"

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  


#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) \
                                 for( int side=0; side<=1; side++ )

// ================================================================================================
/// \brief Constructor for the CgWaveHoltz class
// ================================================================================================
CgWaveHoltz::
CgWaveHoltz( CompositeGrid & cgIn, GenericGraphicsInterface & giIn ) : cg(cgIn), gi(giIn)
{

  dbase.put<int>("debug")=0;
  dbase.put<int>("cgWaveDebugMode")=0;

  real & omega = dbase.put<real>("omega")=30.1;
  dbase.put<real>("Tperiod")=twoPi/omega;
  dbase.put<int>("numPeriods")=10;

  // For multi-frequency WaveHoltz
  int & numberOfFrequencies  = dbase.put<int>("numberOfFrequencies")=1;
  RealArray & frequencyArray = dbase.put<RealArray>("frequencyArray");
  RealArray & periodArray    = dbase.put<RealArray>("periodArray");
  // IntegerArray & numPeriodsArray = dbase.put<IntegerArray>("numPeriodsArray"); // use version in cgWave 
  frequencyArray.redim(numberOfFrequencies);
  frequencyArray(0)=30.1;
  periodArray.redim(numberOfFrequencies);
  periodArray(0) = twoPi/frequencyArray(0);
  // numPeriodsArray.redim(numberOfFrequencies);
  // numPeriodsArray=1; 

  dbase.put<aString>("solverName")="fixedPoint";  // fixedPoint or gmres etc.
  dbase.put<real>("tol")=1.e-4;  // tolerance for Krylov solvers 

  dbase.put<aString>("nameOfGridFile")="unknown";

  dbase.put<int>("maximumNumberOfIterations")=500;
  dbase.put<int>("numberOfIterations")=0; // actual number of iterations taken
  dbase.put<int>("useVariableTolerance")=0; // Vary implicit solver tolerance based on current WaveHoltz residual
  dbase.put<int>("numberOfMatrixVectorMultiplications")=0;

  dbase.put<aString>("krylovType") = "gmres";  // gmres, bicgstab
  dbase.put<int>("gmresRestartLength")=-1; // restart length for WaveHoltz + GMRES

  dbase.put<int>("orderOfAccuracy")=0; 

  dbase.put<int>("useFixedPoint")=0; 

  dbase.put<Real>("convergenceRate")=0.;
  dbase.put<Real>("convergenceRatePerPeriod")=0.;

  dbase.put<Real>("maxResidual")=-1.;
  dbase.put<Real>("maxRes")=-1.; 

  dbase.put<Real>("ppw")=-1.;            // holds actual points per wavelength
  dbase.put<Real>("ppwRuleOfThumb")=-1.; // holds rule-of-thumb points per wavelength

  dbase.put<Real>("cpuWaveHoltz")=0;

  real & omegaSOR = dbase.put<real>("omegaSOR")=1.;

  dbase.put<int>("adjustOmega")=0;  // 1 : choose omega from the symbol of D+t D-t 

  dbase.put<int>("monitorResiduals")=1;  // monitor the residuals at every step
  dbase.put<int>("saveMatabFile")=1;     // save matlab file with residuals etc.
  dbase.put<aString>("matlabFileName")="cgWaveHoltz";  // name of matlab file holding residuals etc. is by default cgWaveHoltz.m

  // dbase.put<realCompositeGridFunction>("vOld");
  dbase.put<realCompositeGridFunction>("residual");
  dbase.put<Real>("numberOfActivePoints")=0.;         // hold number of active points on the grid for scaling L2 norms

  dbase.put<Real>("cpuSolveEigenWave")=0.;

  // Save "residuals" by iteration: ****now in cgWave***
  // resVector(it) = norm( v^{n+1} - v^n )
  // dbase.put<RealArray>("resVector");

  dbase.put<int>("petscIsInitialized")=false;

  FILE *& debugFile = dbase.put<FILE*>("debugFile");
  debugFile = fopen("cgWaveHoltz.debug","w" );        // log file 

  FILE *& logFile = dbase.put<FILE*>("logFile");
  logFile = fopen("cgWaveHoltz.out","w" );        // log file 

  FILE *& checkFile = dbase.put<FILE*>("checkFile");
  // checkFile = fopen("cgWaveHoltz.check","w" );        // for regression and convergence tests
  // fPrintF(checkFile,"# Check file for CgWaveHoltz\n"); // check file has one title line

  // here is the CgWave solver for the time dependent wave equation
  CgWave *& cgWave = dbase.put<CgWave*>("cgWave");
  cgWave = new CgWave(cg,gi);
  
}

// ================================================================================================
/// \brief Destructor for the CgWaveHoltz class
// ================================================================================================
CgWaveHoltz::
~CgWaveHoltz()
{
  fclose(dbase.get<FILE*>("debugFile"));
  fclose(dbase.get<FILE*>("logFile"));
  // fclose(dbase.get<FILE*>("checkFile"));

  delete dbase.get<CgWave*>("cgWave");
}

// ================================================================================================
/// \brief Initialize parameters
// ================================================================================================
int CgWaveHoltz::initialize()
{
  CgWave & cgWave = *dbase.get<CgWave*>("cgWave");

  // ---- adjust periods ----
  const int & numPeriods          = dbase.get<int>("numPeriods");
  const int & numberOfFrequencies = dbase.get<int>("numberOfFrequencies");
  RealArray & frequencyArray      = dbase.get<RealArray>("frequencyArray");
  RealArray & periodArray         = dbase.get<RealArray>("periodArray");

  IntegerArray & numPeriodsArray  = cgWave.dbase.get<IntegerArray>("numPeriodsArray");



  if( frequencyArray(0) > min(frequencyArray) )
  {
    printf("CgWaveHoltz::initialize: ERROR: frequencyArray(0) must be the smallest frequency\n");
    ::display(frequencyArray,"frequencyArray");
    OV_ABORT("ERROR");
  }
  // We may be able to integrate over more periods T2 in the larger time interval T1
  numPeriodsArray.redim(numberOfFrequencies);
  numPeriodsArray(0) = numPeriods;
  for( int freq=1; freq<numberOfFrequencies; freq++ )
    numPeriodsArray(freq) = floor( frequencyArray(freq)/frequencyArray(0))*numPeriodsArray(0); // integrate over this many periods for the "T2" integral

  for( int freq=0; freq<numberOfFrequencies; freq++ )
    periodArray(freq)=numPeriodsArray(freq)*twoPi/frequencyArray(freq);

  // ---------- set values in CgWave ----------

 
  cgWave.dbase.get<int>("numberOfFrequencies") = numberOfFrequencies;  

  RealArray & cgWaveFrequencyArray = cgWave.dbase.get<RealArray>("frequencyArray");
  cgWaveFrequencyArray.redim(numberOfFrequencies);
  cgWaveFrequencyArray = frequencyArray;

  RealArray & cgWavePeriodArray = cgWave.dbase.get<RealArray>("periodArray");
  cgWavePeriodArray.redim(numberOfFrequencies);
  cgWavePeriodArray = periodArray; 

  // IntegerArray & cgWaveNumPeriodsArray = cgWave.dbase.get<IntegerArray>("numPeriodsArray");
  // cgWaveNumPeriodsArray.redim(numberOfFrequencies);
  // cgWaveNumPeriodsArray = numPeriodsArray;


  return 0;
}


// ================================================================================================
/// \brief Set the name of the composite grid file 
// ================================================================================================
int CgWaveHoltz::
setNameOfGridFile( aString & nameOfOGFile )
{
  dbase.get<aString>("nameOfGridFile")=nameOfOGFile;
  return 0;
}


// ================================================================================================
/// \brief Assign parameters 
// ================================================================================================
int CgWaveHoltz::interactiveUpdate()
{

  PlotStuffParameters psp;

  real & omega                    = dbase.get<real>("omega");
  real & Tperiod                  = dbase.get<real>("Tperiod");
  int & numPeriods                = dbase.get<int>("numPeriods");
  real & tol                      = dbase.get<real>("tol");
  int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");
  int & useVariableTolerance      = dbase.get<int>("useVariableTolerance"); // Vary implicit solver tolerance based on current WaveHoltz residual

  aString & krylovType            = dbase.get<aString>("krylovType");     // gmres, bicgstab
  int & gmresRestartLength        = dbase.get<int>("gmresRestartLength"); // restart length for WaveHoltz + GMRES  

  int & numberOfFrequencies       = dbase.get<int>("numberOfFrequencies");
  RealArray & frequencyArray      = dbase.get<RealArray>("frequencyArray");
  RealArray & periodArray         = dbase.get<RealArray>("periodArray");

  int & adjustOmega               = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

  int & monitorResiduals          = dbase.get<int>("monitorResiduals");      // montior the residuals at every step
  int & saveMatlabFile            = dbase.get<int>("saveMatabFile");         // save matlab file with residuals etc.
  aString & matlabFileName        = dbase.get<aString>("matlabFileName");    // name of matlab file holding residuals etc.

  real & omegaSOR                 = dbase.get<real>("omegaSOR");

  CgWave & cgWave                 = *dbase.get<CgWave*>("cgWave");

  int & filterTimeDerivative      = cgWave.dbase.get<int>("filterTimeDerivative");

  int & filterD0t                 = cgWave.dbase.get<int>("filterD0t"); 

  // Build a dialog menu for changing parameters
  GUIState gui;
  DialogData & dialog=gui;

  dialog.setWindowTitle("CgWaveHoltz - Helmholtz Solver");
  dialog.setExitCommand("exit", "exit");

  // dialog.setOptionMenuColumns(1);

  // aString accuracyLabel[] = {"second order", "fourth order", "" };
  // dialog.addOptionMenu("accuracy:", accuracyLabel, accuracyLabel, (orderOfAccuracy==2 ? 0 : 1) );

  aString pbLabels[] = {
                        "grid",
                        "direct solver parameters",
                        "erase",
                        "exit",
			                  ""};
  int numRows=2;
  dialog.setPushButtons( pbLabels, pbLabels, numRows ); 


  aString tbCommands[] = {"save matlab file",
                          "monitor residuals",
                          "adjust omega",
                          "filter time derivative",
                          "filter D0t",
                          "use variable tolerance",
                           ""};
  int tbState[10];
  tbState[0] = saveMatlabFile;
  tbState[1] = monitorResiduals;
  tbState[2] = adjustOmega;
  tbState[3] = filterTimeDerivative;
  tbState[4] = filterD0t;
  tbState[5] = useVariableTolerance;
  int numColumns=1;
  dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 

  // ----- Text strings ------
  const int numberOfTextStrings=20;
  aString textCommands[numberOfTextStrings];
  aString textLabels[numberOfTextStrings];
  aString textStrings[numberOfTextStrings];

  int nt=0;
  textCommands[nt] = "omega";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",omega);  nt++; 

  textCommands[nt] = "number of periods";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numPeriods);  nt++; 

  textCommands[nt] = "number of frequencies";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numberOfFrequencies);  nt++;   

  // Do this for now -- not scalable to many frequencies
  textCommands[nt] = "frequencies";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",frequencyArray(0)); 
  aString buffer; 
  for( int freq=1; freq<numberOfFrequencies; freq++ )
  {
    textStrings[nt] += sPrintF(buffer, ", %g",frequencyArray(freq));   
  }
  nt++; 


  textCommands[nt] = "tol";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",tol);  nt++; 

  textCommands[nt] = "maximum number of iterations";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numPeriods);  nt++; 

  textCommands[nt] = "krylov type";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%s",(const char*)krylovType);  nt++; 

  textCommands[nt] = "gmres restart length";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",gmresRestartLength);  nt++; 

  textCommands[nt] = "matlab filename:";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%s",(const char*)matlabFileName);  nt++; 

  // null strings terminal list
  textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
  dialog.setTextBoxes(textCommands, textLabels, textStrings);

  

  gi.pushGUI(gui);

  aString answer,line;
  char buff[200];
  int len;

  psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
  psp.set(GI_TOP_LABEL,sPrintF(buff,"CgWaveHoltz"));


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
    else if( answer.matches("grid") )
    {
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      PlotIt::plot(gi,cg,psp);                          // plot the grid
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
    }

    else if( answer=="direct solver parameters" )
    {
      printf("Set the Oges parameters for the direct Helmholtz solver\n");
      OgesParameters & par = dbase.put<OgesParameters>("helmholtzOgesParameters");
      par.update( gi,cg );

    }

    else if( dialog.getTextValue(answer,"omega","%e",omega) )
    {
      printF("Setting omega=%g\n",omega);
      cgWave.dbase.get<int>("adjustOmega")= adjustOmega;  // 1 : choose omega from the symbol of D+t D-t 
    }

    else if( dialog.getTextValue(answer,"tol","%e",tol) )
    {
      printF("Setting tol=%g (tolerence for Krylov solvers)\n",tol);
    }
    
    else if( dialog.getTextValue(answer,"number of periods","%i",numPeriods) )
    {
      printF("Setting numPeriods=%i\n",numPeriods);
    }
  
    else if( dialog.getTextValue(answer,"number of frequencies","%i",numberOfFrequencies) )
    {
      printF("Setting numberOfFrequencies=%i\n",numberOfFrequencies);
      frequencyArray.resize(numberOfFrequencies);  frequencyArray=1.;
      periodArray.resize(numberOfFrequencies);     periodArray=twoPi/frequencyArray;
    }
    else if( (len=answer.matches("frequencies")) )
    {
      // Do this for now -- **fix me**
      if( numberOfFrequencies==1 )
        sScanF(answer(len,answer.length()-1),"%e",&frequencyArray(0));
      else if( numberOfFrequencies==2 )
        sScanF(answer(len,answer.length()-1),"%e %e",&frequencyArray(0),&frequencyArray(1));
      else if( numberOfFrequencies==3 )
        sScanF(answer(len,answer.length()-1),"%e %e %e",&frequencyArray(0),&frequencyArray(1),&frequencyArray(2));    
      else if( numberOfFrequencies==4 )
        sScanF(answer(len,answer.length()-1),"%e %e %e %e",&frequencyArray(0),&frequencyArray(1),&frequencyArray(2),&frequencyArray(3)); 
      else if( numberOfFrequencies==5 )
        sScanF(answer(len,answer.length()-1),"%e %e %e %e %e",&frequencyArray(0),&frequencyArray(1),&frequencyArray(2),&frequencyArray(3),&frequencyArray(4));  
      else if( numberOfFrequencies==6 )
        sScanF(answer(len,answer.length()-1),"%e %e %e %e %e %e",
            &frequencyArray(0),&frequencyArray(1),&frequencyArray(2),&frequencyArray(3),&frequencyArray(4),&frequencyArray(5)); 
      // else if( numberOfFrequencies==7 )
      //   sScanF(answer(len,answer.length()-1),"%e %e %e %e %e %e %e",
      //       &frequencyArray(0),&frequencyArray(1),&frequencyArray(2),&frequencyArray(3),&frequencyArray(4),&frequencyArray(5),&frequencyArray(6));             
      else
      {
        // --- general case ---
        stringstream ss;
        ss << answer(len,answer.length()-1);
        string temp;
        int num=0;
        while( !ss.eof() )
        {
          ss >> temp;
          sScanF(temp,"%e",&frequencyArray(num));
          printF("READING: frequencyArray(%d)=%9.3e\n",num,frequencyArray(num));
          num++;
          if( num==numberOfFrequencies ) break; 
        }
        if( num!=numberOfFrequencies )
        {
          printF("cgWaveHoltz:ERROR scanning for frequencies, not enough found. num=%d\n",num);
          OV_ABORT("ERROR");
        }
        // OV_ABORT("STOP HERE FOR NOW");
      }
      // else
      // {
      //   printF("ERROR: two many frequencies to input. FIX ME BILL! \n");
      //   OV_ABORT("error");
      // }                        

      for( int freq=0; freq<numberOfFrequencies; freq++ )
      {
        if( frequencyArray(freq)<=0 )
        {
          printF("ERROR: frequency %d = %12.6e !\n",freq,frequencyArray(freq));
          OV_ABORT("error");
        }
        printF("Setting frequency %d = %12.6e\n",freq,frequencyArray(freq));
        periodArray(freq) = twoPi/frequencyArray(freq); 
      }      
    }

    else if( dialog.getTextValue(answer,"maximum number of iterations","%i",maximumNumberOfIterations) )
    {
      printF("Setting maximumNumberOfIterations=%i\n",maximumNumberOfIterations);
    }

    else if( dialog.getTextValue(answer,"gmres restart length","%e",gmresRestartLength) )
    {
      printF("Setting gmresRestartLength=%i (-1 means use default)\n",gmresRestartLength);
    }  

    else if( dialog.getTextValue(answer,"krylov type","%e",krylovType) )
    {
      printF("Setting krylovType=%s (for WaveHoltz+Krylov)\n",(const char*)krylovType);
    }  

    else if( dialog.getTextValue(answer,"matlab filename:","%s",matlabFileName) )
    {
      printF("Setting matlabFileName=[%s]\n",(const char*)matlabFileName);
    }
    else if( dialog.getToggleValue(answer,"save matlab file",saveMatlabFile) )
    {
      printF("Setting saveMatlabFile=%d: 1=save a matlab file with residuals after waveHoltz solve.\n",saveMatlabFile);
    }
    else if( dialog.getToggleValue(answer,"monitor residuals",monitorResiduals) )
    {
      printF("Setting monitorResiduals=%d: 1=print residuals after each waveHoltz iteration.\n",monitorResiduals);
    }
    else if( dialog.getToggleValue(answer,"adjust omega",adjustOmega) )
    {
      printF("Setting adjustOmega=%d: 1=adjust omega for the discrete symbol pf D+t D-t \n"
             " This will make WaveHoltz solution be closer to the discrete Helhmhotz problem\n",adjustOmega);
    }
    else if( dialog.getToggleValue(answer,"filter time derivative",filterTimeDerivative) )
    {
      printF("Setting filterTimeDerivative=%d: 1=filter time-derivative of the initial condition.\n",filterTimeDerivative);
    } 
    else if( dialog.getToggleValue(answer,"filter D0t",filterD0t) )
    {
      printF("Setting filterD0t=%d: 1=filter D0t directly, 0=use integration by parts formula.\n",filterD0t);
    }        
    else if( dialog.getToggleValue(answer,"use variable tolerance",useVariableTolerance) )
    {
      printF("Setting useVariableTolerance=%d: 1=adjust implicit solver tolerance based on current WaveHoltz residual.\n",useVariableTolerance);
    }   
    else
    {
      printF("CgWaveHoltz:ERROR: unknown answer=[%s]\n",(const char*)answer);
    }
    
  }
  
  gi.popGUI();  // pop dialog

  // Initialize CgWaveHoltz parameters:
  initialize();

  return 0;
}



// ================================================================================================
/// \brief Save results to a matlab file
// ================================================================================================
int CgWaveHoltz::outputMatlabFile( const aString & matlabFileName )
{
  
  const int myid=max(0,Communication_Manager::My_Process_Number);
  if( myid!=0 )
   return 0;

  const aString & solverName               = dbase.get<aString>("solverName");
  const aString & nameOfGridFile           = dbase.get<aString>("nameOfGridFile");
  // const real & omega                       = dbase.get<real>("omega");
  // const real & Tperiod                     = dbase.get<real>("Tperiod");
  // const int & numPeriods                   = dbase.get<int>("numPeriods");
  const real & tol                         = dbase.get<real>("tol");
          
  const int & numberOfFrequencies          = dbase.get<int>("numberOfFrequencies");
  

  const Real & ppw                         = dbase.get<Real>("ppw");            // holds actual points per wavelength
  const Real & ppwRuleOfThumb              = dbase.get<Real>("ppwRuleOfThumb"); // holds rule-of-thumb points per wavelength


  const int & numberOfIterations           = dbase.get<int>("numberOfIterations");
  const Real & convergenceRate             = dbase.get<Real>("convergenceRate");
  const Real & convergenceRatePerPeriod    = dbase.get<Real>("convergenceRatePerPeriod");
  const Real & cpuWaveHoltz                = dbase.get<Real>("cpuWaveHoltz");

  CgWave & cgWave                                           = *dbase.get<CgWave*>("cgWave");

  const real & omega                                        = cgWave.dbase.get<real>("omega");
  const real & Tperiod                                      = cgWave.dbase.get<real>("Tperiod");
  const int & numPeriods                                    = cgWave.dbase.get<int>("numPeriods");
  const IntegerArray & numPeriodsArray                      = cgWave.dbase.get<IntegerArray>("numPeriodsArray");
  const int & adjustOmega                                   = cgWave.dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  const int & orderOfAccuracy                               = cgWave.dbase.get<int>("orderOfAccuracy");
  const int & orderOfAccuracyInTime                         = cgWave.dbase.get<int>("orderOfAccuracyInTime");
  const real & c                                            = cgWave.dbase.get<real>("c");
  const CgWave::TimeSteppingMethodEnum & timeSteppingMethod = cgWave.dbase.get<CgWave::TimeSteppingMethodEnum>("timeSteppingMethod");
  const RealArray & frequencyArray                          = cgWave.dbase.get<RealArray>("frequencyArray");
  const RealArray & periodArray                             = cgWave.dbase.get<RealArray>("periodArray");    
  const RealArray & frequencyArrayAdjusted                  = cgWave.dbase.get<RealArray>("frequencyArrayAdjusted");
  const RealArray & periodArrayAdjusted                     = cgWave.dbase.get<RealArray>("periodArrayAdjusted");
  const RealArray & sigma                                   = cgWave.dbase.get<RealArray>("sigma");
  const Real & dt                                           = cgWave.dbase.get<Real>("dtUsed");
  const int & minStepsPerPeriod                             = cgWave.dbase.get<int>("minStepsPerPeriod");
  const int & computeEigenmodes                             = cgWave.dbase.get<int>("computeEigenmodes");
  RealArray & resVector                                     = cgWave.dbase.get<RealArray>("resVector");

  const int & numCompWaveHoltz                              = cgWave.dbase.get<int>("numCompWaveHoltz");
  const int & filterTimeDerivative                          = cgWave.dbase.get<int>("filterTimeDerivative");
  const int & useSuperGrid                                  = cgWave.dbase.get<int>("useSuperGrid");
  const RealArray & timing                                  = cgWave.timing;

  const Real & kxPlaneWave                                  = cgWave.dbase.get<Real>("kxPlaneWave");
  const Real & kyPlaneWave                                  = cgWave.dbase.get<Real>("kyPlaneWave");
  const Real & kzPlaneWave                                  = cgWave.dbase.get<Real>("kzPlaneWave");
  RealArray & bImp                                          = cgWave.dbase.get<RealArray>("bImp");
  RealArray & cImp                                          = cgWave.dbase.get<RealArray>("cImp");


  // printF("outputMatlab: numberOfIterations=%d\n",numberOfIterations);
  // ::display(resVector,"resVector");

  // aString fileName="cgWaveHoltz.m"; // allow this to be specified
  // aString & matlabFileName           = dbase.get<aString>("matlabFileName");    // name of matlab file holding residuals etc.

  // aString suffix = sPrintF("FD%d%dTS%s",orderOfAccuracyInTime,orderOfAccuracy,
  //                          timeSteppingMethod==CgWave::implicitTimeStepping ? "I" : "E");
  // if( numberOfFrequencies>1 )
  //   suffix = suffix + sPrintF("Nf%d",numberOfFrequencies);

  // aString fileName = matlabFileName + suffix + ".m";  

  aString fileName = matlabFileName + ".m";



  // gi.inputString(fileName,sPrintF(answer,"Enter the name of the matlab file (default=%s)\n",(const char*)fileName));
  // if( fileName=="" )
  //   fileName="cgWaveHoltz.m";

  FILE *matlabFile = fopen((const char*)fileName,"w" ); 

  // Get the current date
  time_t *tp= new time_t;
  time(tp);
  // tm *ptm=localtime(tp);
  const char *dateString = ctime(tp);
  fPrintF(matlabFile,"%% File created by CgWaveHoltz %s",dateString);
  fPrintF(matlabFile,"%% Residuals versus iteration.\n");
  delete tp;

  if( cgWave.dbase.has_key("userDefinedForcingData") )
  {
    DataBase & db = cgWave.dbase.get<DataBase>("userDefinedForcingData");
    aString & option= db.get<aString>("option");
    fPrintF(matlabFile,"%% User defined forcing:\nuserDefinedForcingOption='%s';\n",(const char*)option);
  }

  fPrintF(matlabFile,"solverName='%s';\n",(const char*)solverName);
  fPrintF(matlabFile,"gridName=\'%s\';\n",(const char*)nameOfGridFile);
  const real & numberOfGridPoints    = cgWave.dbase.get<real>("numberOfGridPoints");
  fPrintF(matlabFile,"numberOfGridPoints=%12.6e;\n",numberOfGridPoints);
  fPrintF(matlabFile,"timeSteppingMethod='%s';\n",(timeSteppingMethod==CgWave::explicitTimeStepping ? "explicit" : "implicit") );

  Real beta2=bImp(0), beta4=bImp(1); 
  Real alpha2 = (1.-beta2)/2.;
  Real alpha4 = (alpha2-beta4-1./12.)/2.; 
  fPrintF(matlabFile,"beta2=%g; beta4=%g; alpha2=%g; alpha4=%g;\n",beta2,beta4,alpha2,alpha4);

  fPrintF(matlabFile,"totalTime=%9.3e;\n",cpuWaveHoltz); 
  fPrintF(matlabFile,"timeForAdvance=%9.3e;\n",timing(CgWave::timeForAdvance));
  if( dbase.has_key("cpuAugmentedInitialMatVects") )
    fPrintF(matlabFile,"cpuAugmentedInitialMatVects=%9.3e;\n",dbase.get<Real>("cpuAugmentedInitialMatVects"));
  if( dbase.has_key("cpuAugmentedQR") )
    fPrintF(matlabFile,"cpuAugmentedQR=%9.3e;\n",dbase.get<Real>("cpuAugmentedQR"));

  const int & numberOfMatrixVectorMultiplications = dbase.get<int>("numberOfMatrixVectorMultiplications");
  fPrintF(matlabFile,"numberOfMatrixVectorMultiplications=%d;\n",numberOfMatrixVectorMultiplications);

  const Real cpuSolveEigenWave = computeEigenmodes ? dbase.get<Real>("cpuSolveEigenWave") : 0.;
  fPrintF(matlabFile,"cpuSolveEigenWave=%9.2e;\n",cpuSolveEigenWave);
 
  const int & totalImplicitIterations       = cgWave.dbase.get<int>("totalImplicitIterations"); 
  const int & totalImplicitSolves           = cgWave.dbase.get<int>("totalImplicitSolves"); 
  const Real aveNumberOfImplicitIterations  = totalImplicitIterations/Real(max(1,totalImplicitSolves));  
  fPrintF(matlabFile,"aveNumberOfImplicitIterations=%5.1f;\n",aveNumberOfImplicitIterations);

  fPrintF(matlabFile,"omega=%20.14e;\n",omega);
  fPrintF(matlabFile,"kxPlaneWave=%10.3e; kyPlaneWave=%10.3e; kzPlaneWave=%10.3e;  %% plane wave number\n",kxPlaneWave/twoPi,kyPlaneWave/twoPi,kzPlaneWave/twoPi);

  fPrintF(matlabFile,"ppw=%10.3e;             %% actual points-per-wavelength\n",ppw);
  fPrintF(matlabFile,"ppwRuleOfThumb=%10.3e;  %% rule-of-thumb points-per-wavelength\n",ppwRuleOfThumb);

  int numPerLine=5;

  fPrintF(matlabFile,"numberOfFrequencies=%d;\n",numberOfFrequencies);
  fPrintF(matlabFile,"frequencyArray=[");
  for( int freq=0; freq<numberOfFrequencies; freq++ )
  {
    fPrintF(matlabFile,"%14.6e ",frequencyArray(freq));
    if( (freq % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
  }
  fPrintF(matlabFile,"];\n");

  fPrintF(matlabFile,"frequencyArrayAdjusted=[");
  for( int freq=0; freq<numberOfFrequencies; freq++ )
  {
    if( adjustOmega )
      fPrintF(matlabFile,"%14.6e ",frequencyArrayAdjusted(freq));
    else
      fPrintF(matlabFile,"%14.6e ",frequencyArray(freq));
    if( (freq % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
  }
  fPrintF(matlabFile,"];\n");  

  fPrintF(matlabFile,"periodArray=[");
  for( int freq=0; freq<numberOfFrequencies; freq++ )
  {
    fPrintF(matlabFile,"%14.6e ",periodArray(freq));
    if( (freq % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
  }
  fPrintF(matlabFile,"];\n");  

  fPrintF(matlabFile,"periodArrayAdjusted=[");
  for( int freq=0; freq<numberOfFrequencies; freq++ )
  {
    if( adjustOmega )
      fPrintF(matlabFile,"%14.6e ",periodArrayAdjusted(freq));
    else
      fPrintF(matlabFile,"%14.6e ",periodArray(freq));
    if( (freq % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
  }
  fPrintF(matlabFile,"];\n"); 

  fPrintF(matlabFile,"numPeriodsArray=[");
  for( int freq=0; freq<numberOfFrequencies; freq++ )
  {
    fPrintF(matlabFile,"%d ",numPeriodsArray(freq));
    if( (freq % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
  }
  fPrintF(matlabFile,"];\n");    

  fPrintF(matlabFile,"filterTimeDerivative=%d;\n",filterTimeDerivative);
  fPrintF(matlabFile,"numCompWaveHoltz=%d;\n",numCompWaveHoltz);
  fPrintF(matlabFile,"useSuperGrid=%d;\n",useSuperGrid);

  fPrintF(matlabFile,"computeEigenmodes=%d;\n",computeEigenmodes);
  if( computeEigenmodes )
  {

    CgWave::EigenSolverEnum & eigenSolver = cgWave.dbase.get<CgWave::EigenSolverEnum>("eigenSolver"); 
    const bool useFixedPoint =  eigenSolver==CgWave::fixedPointEigenSolver;

    aString eigenSolverName = eigenSolver==CgWave::defaultEigenSolver           ? "KrylovSchur" :
                              eigenSolver==CgWave::KrylovSchurEigenSolver       ? "KrylovSchur" :
                              eigenSolver==CgWave::ArnoldiEigenSolver           ? "Arnoldi" : 
                              eigenSolver==CgWave::ArpackEigenSolver            ? "Arpack" : 
                              eigenSolver==CgWave::fixedPointEigenSolver        ? "fixedPoint" : 
                              eigenSolver==CgWave::subspaceIterationEigenSolver ? "subspaceIteration" : 
                                                                            "unknown";
    fPrintF(matlabFile,"eigenSolver='%s';\n",(const char*)eigenSolverName);

    const int & numEigenVectors = dbase.get<int>("numEigenVectors");
    fPrintF(matlabFile,"numEigenVectors=%d; %% number of computed eigenpairs\n",numEigenVectors);

    const int & numEigsToCompute                    = cgWave.dbase.get<int>("numEigsToCompute");
    const int & numArnoldiVectors                   = cgWave.dbase.get<int>("numArnoldiVectors"); 
    fPrintF(matlabFile,"numEigsRequested=%d;  %% number of requested eigenpairs\n",numEigsToCompute);
    fPrintF(matlabFile,"numArnoldiVectors=%d; %% number of Arnolidi vectors\n",numArnoldiVectors);

    const int & numberOfStepsTaken                  = cgWave.dbase.get<int>("numberOfStepsTaken");      
    const int & numberOfMatrixVectorMultiplications = dbase.get<int>("numberOfMatrixVectorMultiplications");
    const Real numWaveSolvesPerEig = numberOfMatrixVectorMultiplications/max(Real(numEigenVectors),1.);
    const Real numTimeStepsPerEig = numberOfStepsTaken/max(Real(numEigenVectors),1.);
    fPrintF(matlabFile,"numberOfMatrixVectorMultiplications=%d;\n",numberOfMatrixVectorMultiplications);
    fPrintF(matlabFile,"numWaveSolvesPerEig=%6.2f;\n",numWaveSolvesPerEig);
    fPrintF(matlabFile,"numTimeStepsPerEig=%6.2f;\n",numTimeStepsPerEig);

    const RealArray & eigenValues = dbase.get<RealArray>("eigenValues");
    fPrintF(matlabFile,"%% computed eigenvalues:\n");
    fPrintF(matlabFile,"eigenValues=[%24.16e",eigenValues(0));
    for( int ie=1; ie<numEigenVectors; ie++ )
    {
      fPrintF(matlabFile,", %24.16e",eigenValues(ie));
      if( (ie % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
    }
    fPrintF(matlabFile,"];\n");

    if( !useFixedPoint )
    {
      if( cgWave.dbase.has_key("uev") )
      {
        // --- True values are known ----  
        // --- Output errors in eigenvalues and eigenvectors ---    
        const RealArray & relErrEigenvalue  = dbase.get<RealArray>("relErrEigenvalue");
        const RealArray & relErrEigenvector = dbase.get<RealArray>("relErrEigenvector");
        const RealArray & eigenPairResidual = dbase.get<RealArray>("eigenPairResidual");

        const int numPerLine = 10;

        fPrintF(matlabFile,"%% relative error in eigenvalues:\n");
        fPrintF(matlabFile,"relErrEigenvalue=[%10.2e",relErrEigenvalue(0));
        for( int ie=1; ie<numEigenVectors; ie++ )
        {
          fPrintF(matlabFile,",%9.2e",relErrEigenvalue(ie));
          if( (ie % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
        }
        fPrintF(matlabFile,"];\n");

        fPrintF(matlabFile,"%% relative error in eigenvectors:\n");
        fPrintF(matlabFile,"relErrEigenvector=[%10.2e",relErrEigenvector(0));
        for( int ie=1; ie<numEigenVectors; ie++ )
        {
          fPrintF(matlabFile,",%9.2e",relErrEigenvector(ie));
          if( (ie % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
        }
        fPrintF(matlabFile,"];\n");
       
        fPrintF(matlabFile,"%% relative residual in eigenvalue equation:\n");
        fPrintF(matlabFile,"eigenPairResidual=[%10.2e",eigenPairResidual(0));
        for( int ie=1; ie<numEigenVectors; ie++ )
        {
          fPrintF(matlabFile,",%9.2e",eigenPairResidual(ie));
          if( (ie % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
        }
        fPrintF(matlabFile,"];\n");

      }

    }

    if( useFixedPoint )
    {
      RealArray & eigsRayleighQuotient      = cgWave.dbase.get<RealArray>("eigsRayleighQuotient");
      RealArray & eigsRayleighRitz          = cgWave.dbase.get<RealArray>("eigsRayleighRitz");
      const int & iterationRR               = cgWave.dbase.get<int>("iterationRR");
      const int & iterationStartRR          = cgWave.dbase.get<int>("iterationStartRR"); 

      fPrintF(matlabFile,"eigsRayleighQuotient=[");
      const int numIts = min(numberOfIterations,eigsRayleighQuotient.getLength(0));
      for( int it=0; it<numIts; it++ )
      {
        fPrintF(matlabFile,"%24.16e ",eigsRayleighQuotient(it));
        if( (it % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
      }
      fPrintF(matlabFile,"];\n"); 

      fPrintF(matlabFile,"iterationStartRR=%d;\n",iterationStartRR);
      fPrintF(matlabFile,"eigsRayleighRitz=[");
      for( int it=0; it<iterationRR; it++ )
      {
        fPrintF(matlabFile,"%24.16e ",eigsRayleighRitz(it));
        if( (it % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
      }
      fPrintF(matlabFile,"];\n");   

      if( cgWave.dbase.has_key("errEigenvector") )
      {
        // --- output errors in eigenvalues and eigenvectors if they have been computed ---

        RealArray & errEigenvalue  = cgWave.dbase.get<RealArray>("errEigenvalue");
        RealArray & errEigenvector = cgWave.dbase.get<RealArray>("errEigenvector");

        const int numIts = min(numberOfIterations,errEigenvalue.getLength(0));
        fPrintF(matlabFile,"errEigenvalue=[");
        for( int it=0; it<numIts; it++ )
        {
          fPrintF(matlabFile,"%9.3e ",errEigenvalue(it));
          if( (it % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
        }
        fPrintF(matlabFile,"];\n");

        fPrintF(matlabFile,"errEigenvector=[");
        for( int it=0; it<numIts; it++ )
        {
          fPrintF(matlabFile,"%9.3e ",errEigenvector(it));
          if( (it % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
        }
        fPrintF(matlabFile,"];\n");      

      }

      if( cgWave.dbase.has_key("errEigenvectorRR") )
      {
        // --- output errors in eigenvalues and eigenvectors if they have been computed ---

        RealArray & errEigenvalueRR  = cgWave.dbase.get<RealArray>("errEigenvalueRR");
        RealArray & errEigenvectorRR = cgWave.dbase.get<RealArray>("errEigenvectorRR");

        const int numIts = min(numberOfIterations,errEigenvalueRR.getLength(0));
        fPrintF(matlabFile,"errEigenvalueRR=[");
        for( int it=0; it<numIts; it++ )
        {
          fPrintF(matlabFile,"%9.3e ",errEigenvalueRR(it));
          if( (it % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
        }
        fPrintF(matlabFile,"];\n");

        fPrintF(matlabFile,"errEigenvectorRR=[");
        for( int it=0; it<numIts; it++ )
        {
          fPrintF(matlabFile,"%9.3e ",errEigenvectorRR(it));
          if( (it % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
        }
        fPrintF(matlabFile,"];\n");      

      }
    } // end if fixPoint

  }

  fPrintF(matlabFile,"c=%20.14e;\n",c);
  fPrintF(matlabFile,"numPeriods=%d; %% (integrate over this many periods per Wave-Holtz iteration.\n",numPeriods);
  fPrintF(matlabFile,"minStepsPerPeriod=%d; %% (min time-steps per period for implicit time-stepping.\n",minStepsPerPeriod);
  fPrintF(matlabFile,"convergenceRate=%12.4e;\n",convergenceRate);
  fPrintF(matlabFile,"convergenceRatePerPeriod=%12.4e;\n",convergenceRatePerPeriod);
  fPrintF(matlabFile,"adjustOmega=%d; %% (1= adjust omega to account for discrete symbol of D+t D-t).\n",adjustOmega);
  fPrintF(matlabFile,"orderOfAccuracy=%d;\n",orderOfAccuracy);
  fPrintF(matlabFile,"orderOfAccuracyInTime=%d;\n",orderOfAccuracyInTime);
  fPrintF(matlabFile,"tol=%12.4e;\n",tol);


  // printF("saveMatlab: dt=%9.3e\n",dt);
  fPrintF(matlabFile,"dt=%20.14e;\n",dt);

  // --- Save the quadrature weights ---
  const int numTimeSteps = sigma.getLength(0);  //   Note: this count includes step 0
  // fPrintF(matlabFile,"numTimeSteps=%d;\n",numTimeSteps);
  numPerLine=40;
  fPrintF(matlabFile,"%% sigma(step,freq) = quadrature weights\n");
  fPrintF(matlabFile,"%% sigma=zeros(%d,%d)\n",numTimeSteps,numberOfFrequencies);
  fPrintF(matlabFile,"sigma=[");
  for( int step=0; step<numTimeSteps; step++ )
  {
    for( int freq=0; freq<numberOfFrequencies; freq++ )  
    {
      fPrintF(matlabFile,"%12.6e",sigma(step,freq));
      if( freq<numberOfFrequencies-1 ) fPrintF(matlabFile,", ");

      if( (step % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
    }
    fPrintF(matlabFile,"; ");
  }
  fPrintF(matlabFile,"];\n");


  numPerLine=40;
  fPrintF(matlabFile,"%% itv = iteration number\n");
  fPrintF(matlabFile,"itv=[");
  for( int i=0; i<numberOfIterations; i++ )
  {
    fPrintF(matlabFile,"%i ",i);
    if( (i % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
  }
  fPrintF(matlabFile,"];\n");

  // ---- save the iteration history of the residual ----
  numPerLine=10;
  fPrintF(matlabFile,"res=[ ...\n");
  for( int i=0; i<numberOfIterations; i++ )
  {
    fPrintF(matlabFile,"%17.10e ",resVector(i));
    if( (i % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
  }
  fPrintF(matlabFile,"];\n");

  // ---- Save info about deflation ----
  fPrintF(matlabFile,"%% ----- Info about deflation ----\n");
  const int & deflateWaveHoltz = cgWave.dbase.get<int>("deflateWaveHoltz");
  const int & numToDeflate     = cgWave.dbase.get<int>("numToDeflate");

  fPrintF(matlabFile,"deflateWaveHoltz=%d;\n",deflateWaveHoltz);
  fPrintF(matlabFile,"numToDeflate=%d;\n",numToDeflate);
  if( cgWave.dbase.has_key("eig") )
  {
    const RealArray & eig                    = cgWave.dbase.get<RealArray>("eig");
    const IntegerArray & eigNumbersToDeflate = cgWave.dbase.get<IntegerArray>("eigNumbersToDeflate");
    IntegerArray & eigMultiplicity           = cgWave.dbase.get<IntegerArray>("eigMultiplicity");
    const int numberOfEigenvectors           = eig.getBound(1) - eig.getBase(1) + 1;

    fPrintF(matlabFile,"%% True (discrete) eigenvalues\n");
    fPrintF(matlabFile,"lambdav=[ ...\n");
    for( int i=0; i<numberOfEigenvectors; i++ )
    {
      fPrintF(matlabFile,"%22.16e ",eig(0,i));
      if( i<numberOfEigenvectors-1 ) fPrintF(matlabFile,", ");
      if( (i % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
    }
    fPrintF(matlabFile,"];\n");

    fPrintF(matlabFile,"eigMultiplicity=[ ...\n");
    for( int i=0; i<numberOfEigenvectors; i++ )
    {
      fPrintF(matlabFile,"%d ",eigMultiplicity(i));
      if( i<numberOfEigenvectors-1 ) fPrintF(matlabFile,", ");
      if( (i % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
    }
    fPrintF(matlabFile,"];\n");    

    fPrintF(matlabFile,"eigsToDeflate=[ ...\n");
    for( int i=0; i<numToDeflate; i++ )
    {
      fPrintF(matlabFile,"%d ",eigNumbersToDeflate(i));
      if( i<numToDeflate-1 ) fPrintF(matlabFile,", ");
      if( (i % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
    }
    fPrintF(matlabFile,"];\n");
  }
  else
  {
    fPrintF(matlabFile,"eig=-1;\n");
    fPrintF(matlabFile,"eigsToDeflate=-1;\n");
  }


  fclose(matlabFile);

  printF("CgWaveHoltz::saved results in matlab file=[%s]\n",(const char*)fileName);


  return 0;
}


// ================================================================================================
/// \brief Save results to a check file.
// ================================================================================================
int CgWaveHoltz::outputCheckFile( const RealArray & maxResArray, const Real errorBetweenWaveHoltzAndHelmholtz )
{

  CgWaveHoltz & cgWaveHoltz = *this;
  CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");

  const Real & convergenceRate             = cgWaveHoltz.dbase.get<Real>("convergenceRate");
  const Real & convergenceRatePerPeriod    = cgWaveHoltz.dbase.get<Real>("convergenceRatePerPeriod");
  const Real & omega                       = cgWaveHoltz.dbase.get<real>("omega");
  const int & numPeriods                   = cgWaveHoltz.dbase.get<int>("numPeriods");
  const Real & tol                         = cgWaveHoltz.dbase.get<real>("tol");
  const int & adjustOmega                  = cgWaveHoltz.dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  // const int & maximumNumberOfIterations    = cgWaveHoltz.dbase.get<int>("maximumNumberOfIterations");
  const int & numberOfIterations           = cgWaveHoltz.dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken
  const int & numberOfFrequencies          = cgWaveHoltz.dbase.get<int>("numberOfFrequencies"); 
  const RealArray & frequencyArray         = cgWaveHoltz.dbase.get<RealArray>("frequencyArray"); 
  const int & useFixedPoint                = cgWaveHoltz.dbase.get<int>("useFixedPoint"); 

  const Real & maxRes                      = dbase.get<Real>("maxRes");

  const int & orderOfAccuracy                 = cgWave.dbase.get<int>("orderOfAccuracy"); 
  const int & useAugmentedGmres               = cgWave.dbase.get<int>("useAugmentedGmres"); 
  const int & augmentedVectorsAreEigenvectors = cgWave.dbase.get<int>("augmentedVectorsAreEigenvectors");  // 1 = augmented vectors are true discrete eigenvectors
  const int & deflateWaveHoltz                = cgWave.dbase.get<int>("deflateWaveHoltz"); 
  const int & numToDeflate                    = cgWave.dbase.get<int>("numToDeflate");   
  const RealArray & frequencyArrayAdjusted    = cgWave.dbase.get<RealArray>("frequencyArrayAdjusted");
  const IntegerArray & numPeriodsArray        = cgWave.dbase.get<IntegerArray>("numPeriodsArray");
  const RealArray & periodArray               = cgWave.dbase.get<RealArray>("periodArray"); 
  const int & numCompWaveHoltz                = cgWave.dbase.get<int>("numCompWaveHoltz");
  const int & filterTimeDerivative            = cgWave.dbase.get<int>("filterTimeDerivative");

  const int & numberOfStepsTaken              = cgWave.dbase.get<int>("numberOfStepsTaken");      
  const int & numberOfStepsPerSolve           = cgWave.dbase.get<int>("numberOfStepsPerSolve");
  const Real & cfl                            = cgWave.dbase.get<real>("cfl");
  const Real & c                              = cgWave.dbase.get<real>("c");
  const Real & damp                           = cgWave.dbase.get<real>("damp");
  const Real & dt                             = cgWave.dbase.get<real>("dtUsed");
  const int & minStepsPerPeriod               = cgWave.dbase.get<int>("minStepsPerPeriod");

  const CgWave::TimeSteppingMethodEnum & timeSteppingMethod = cgWave.dbase.get<CgWave::TimeSteppingMethodEnum>("timeSteppingMethod");



  aString timeSteppingName = (timeSteppingMethod==CgWave::explicitTimeStepping ? "explicit" : "implicit");
  
  aString & nameOfGridFile = cgWave.dbase.get<aString>("nameOfGridFile");
  int firstChar=0; 
  for( int i=nameOfGridFile.length()-1; i>=0; i-- )
  {
    if( nameOfGridFile[i]=='/' ){ firstChar=i+1; break; } // start from end, work backwards and look for a directory symbol
  }
  int lastChar=firstChar; 
  for( int i=firstChar; i<=nameOfGridFile.length()-1; i++ )
  {
    if( nameOfGridFile[i]=='.' ){ lastChar=i-1; break; } // remove suffix: .order2.hdf
  }

  aString gridNameNoPrefix = nameOfGridFile(firstChar,lastChar);
  FILE *& checkFile = dbase.get<FILE*>("checkFile");
  checkFile = fopen("cgWaveHoltz.check","w" );      
  assert( checkFile != NULL );      

  // Get the current date
  time_t *tp= new time_t;
  time(tp);
  const char *dateString = ctime(tp);
  fPrintF(checkFile,"# Check file for cgWaveHoltz, grid=%s, %s",(const char*)gridNameNoPrefix,dateString);  // Note: dateString include newline
  delete tp; 

  fPrintF(checkFile,"grid=%s;\n",(const char*)gridNameNoPrefix);
  fPrintF(checkFile,"timeStepping=%s;\n",(const char*)timeSteppingName);
  fPrintF(checkFile,"orderOfAccuracy=%d;\n",orderOfAccuracy);
  fPrintF(checkFile,"numPeriods=%d;\n",numPeriods);

  fPrintF(checkFile,"numberOfFrequencies=%d;\n",numberOfFrequencies);
  for( int freq=0; freq<numberOfFrequencies; freq++ )
  {
    fPrintF(checkFile," freq=%2d, omega=%8.3f",freq,frequencyArray(freq));
    if( adjustOmega )
      fPrintF(checkFile,", (adjusted=%8.3f)",frequencyArrayAdjusted(freq));

    const Real T = twoPi/frequencyArray(freq);  // period
    const Real Tbar  = T*numPeriodsArray(freq);       // periods that fit in time interval
    const Real Tbar0 = periodArray(0);          // final time integrated to 
    fPrintF(checkFile," T=%8.5f, Tbar=%8.5f, numPeriods=%3d, Tbar(0)/T=%6.2f",T,Tbar,numPeriodsArray(freq),Tbar0/T);
    fPrintF(checkFile,"\n");
  }

  fPrintF(checkFile,"numCompWaveHoltz=%d;\n",numCompWaveHoltz);
  fPrintF(checkFile,"filterTimeDerivative=%d;\n",filterTimeDerivative);

  fPrintF(checkFile,"adjustOmega=%d;\n",adjustOmega);
  fPrintF(checkFile,"useFixedPoint=%d;\n",useFixedPoint);
  fPrintF(checkFile,"useAugmentedGmres=%d;\n",useAugmentedGmres);
  fPrintF(checkFile,"minStepsPerPeriod=%d;\n",minStepsPerPeriod);
  fPrintF(checkFile,"numberOfStepsPerSolve=%d;\n",numberOfStepsPerSolve);

  fPrintF(checkFile,"convergenceRate=%5.3f;\n",convergenceRate);
  fPrintF(checkFile,"numberOfIterations=%d;\n",numberOfIterations);

  fPrintF(checkFile,"deflateWaveHoltz=%d;\n",deflateWaveHoltz);
  fPrintF(checkFile,"numToDeflate=%d;\n",numToDeflate);

  if( cgWave.dbase.has_key("eig") )
  {
    // --- Output deflation info ---
    const RealArray & eig                    = cgWave.dbase.get<RealArray>("eig");
    const IntegerArray & eigNumbersToDeflate = cgWave.dbase.get<IntegerArray>("eigNumbersToDeflate");
    IntegerArray & eigMultiplicity           = cgWave.dbase.get<IntegerArray>("eigMultiplicity");
    const int numberOfEigenvectors           = eig.getBound(1) - eig.getBase(1) + 1;

    int numPerLine=10;

    fPrintF(checkFile,"%% True (discrete) eigenvalues\n");
    fPrintF(checkFile,"lambdav=[ ...\n");
    for( int i=0; i<numberOfEigenvectors; i++ )
    {
      fPrintF(checkFile,"%22.16e ",eig(0,i));
      if( i<numberOfEigenvectors-1 ) fPrintF(checkFile,", ");
      if( (i % numPerLine)==numPerLine-1 ) fPrintF(checkFile,"...\n");
    }
    fPrintF(checkFile,"];\n");

    fPrintF(checkFile,"eigMultiplicity=[ ...\n");
    for( int i=0; i<numberOfEigenvectors; i++ )
    {
      fPrintF(checkFile,"%d ",eigMultiplicity(i));
      if( i<numberOfEigenvectors-1 ) fPrintF(checkFile,", ");
      if( (i % numPerLine)==numPerLine-1 ) fPrintF(checkFile,"...\n");
    }
    fPrintF(checkFile,"];\n");    

    fPrintF(checkFile,"eigsToDeflate=[ ...\n");
    for( int i=0; i<numToDeflate; i++ )
    {
      fPrintF(checkFile,"%d ",eigNumbersToDeflate(i));
      if( i<numToDeflate-1 ) fPrintF(checkFile,", ");
      if( (i % numPerLine)==numPerLine-1 ) fPrintF(checkFile,"...\n");
    }
    fPrintF(checkFile,"];\n");
  }  

  if( filterTimeDerivative==0 )
    fPrintF(checkFile,"maxRes=%9.3e;\n",maxRes);
  else
  {
    fPrintF(checkFile,"maxResReal=%8.2e;\n",maxResArray(0));
    fPrintF(checkFile,"maxResImag=%8.2e;\n",maxResArray(1));
  }
  fPrintF(checkFile,"errorBetweenWaveHoltzAndHelmholtz=%9.3e;\n",errorBetweenWaveHoltzAndHelmholtz);

  // fPrintF(checkFile,"numberOfStepsPerSolve=%d;\n",numberOfStepsPerSolve);
  // fPrintF(checkFile,"numEigsRequested=%d;\n",numEigsToCompute);
  // fPrintF(checkFile,"numEigsComputed=%d;\n",numEigenVectors);
  // fPrintF(checkFile,"numArnoldiVectors=%d;\n",numArnoldiVectors);
  // fPrintF(checkFile,"numWaveSolves=%d;\n",numberOfMatrixVectorMultiplications);
  // fPrintF(checkFile,"maxEigErr=%9.2e;\n",maxEigErr);
  // fPrintF(checkFile,"maxEvectErr=%9.2e;\n",maxEvectErr);
  // fPrintF(checkFile,"maxEigResid=%9.2e;\n",maxEigResid);
  fclose(checkFile);

  printF("Wrote results to the check file [cgWaveHoltz.check]\n");

  return 0;
}

// ================================================================================================
/// \brief Setup grids and grid functions
// ================================================================================================
int CgWaveHoltz::setup()
{
  const int & numberOfFrequencies = dbase.get<int>("numberOfFrequencies");

  CgWave & cgWave                  = *dbase.get<CgWave*>("cgWave");
  const int & filterTimeDerivative = cgWave.dbase.get<int>("filterTimeDerivative");
  int & numCompWaveHoltz           = cgWave.dbase.get<int>("numCompWaveHoltz");
  numCompWaveHoltz = filterTimeDerivative ? 2 : numberOfFrequencies;  

  printF("cgWaveHoltz::setup: numberOfFrequencies=%d, filterTimeDerivative=%d, numCompWaveHoltz=%d\n",numberOfFrequencies,filterTimeDerivative,numCompWaveHoltz);

  // realCompositeGridFunction & vOld = dbase.get<realCompositeGridFunction>("vOld");
  // Range all;
  // vOld.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);

  realCompositeGridFunction & residual = dbase.get<realCompositeGridFunction>("residual");
  Range all;
  residual.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz);
  for( int freq=0; freq<numCompWaveHoltz; freq++ )
  {
    residual.setName(sPrintF("r%d",freq),freq);
  }

  // v=0;
  return 0;
}




// ================================================================================================
/// \brief Solve the equations over one or more periods 
///
// ================================================================================================
int CgWaveHoltz::
solve()
{
  const int myid=max(0,Communication_Manager::My_Process_Number);
  const int np = Communication_Manager::numberOfProcessors();

  GenericGraphicsInterface & ps = gi;
  PlotStuffParameters psp;


  const real & omega        = dbase.get<real>("omega");
  real & Tperiod            = dbase.get<real>("Tperiod");
  const int & numPeriods    = dbase.get<int>("numPeriods");
  const int & adjustOmega   = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  const int & debug         = dbase.get<int>("debug");
  const real & omegaSOR     = dbase.get<real>("omegaSOR");
  const real & tol          = dbase.get<real>("tol");

  const int & numberOfFrequencies = dbase.get<int>("numberOfFrequencies");
  RealArray & frequencyArray      = dbase.get<RealArray>("frequencyArray");
  RealArray & periodArray         = dbase.get<RealArray>("periodArray");
  

  // here is the CgWave solver for the time dependent wave equation
  CgWave & cgWave = *dbase.get<CgWave*>("cgWave");
  if( omega!=0.0 )
    Tperiod=numPeriods*twoPi/omega;  
  else
    Tperiod=1;

  IntegerArray & numPeriodsArray    = cgWave.dbase.get<IntegerArray>("numPeriodsArray");
  const int & numCompWaveHoltz      = cgWave.dbase.get<int>("numCompWaveHoltz");
  const int & filterTimeDerivative  = cgWave.dbase.get<int>("filterTimeDerivative"); 

  printF("CgWaveHoltz::solve: setting tFinal = Tperiod*numPeriods = %9.3e (numPeriods=%d) \n",Tperiod,numPeriods);
 
  // // ---- adjust periods ----
  // IntegerArray & numPeriodsArray = dbase.put<IntegerArray>("numPeriodsArray");
  // // We may be able to integrate over more periods T2 in the larger time interface T1
  // numPeriodsArray = numPeriods;
  // for( int freq=1; freq<numberOfFrequencies; freq++ )
  //   numPeriodsArray(freq) = floor( frequencyArray(freq)/frequencyArray(0))*numPeriodsArray(0); // integrate over this many periods for the "T2" integral

  // for( int freq=0; freq<numberOfFrequencies; freq++ )
  //   periodArray(freq)=numPeriodsArray(freq)*twoPi/frequencyArray(freq);

  // --- set values in CgWave:  *** COULD DO BETTER ***

  cgWave.dbase.get<real>("omega")     = omega;        // ** FIX ME **
  cgWave.dbase.get<real>("tFinal")    = Tperiod;      // ** FIX ME **
  cgWave.dbase.get<real>("Tperiod")   = Tperiod;      // ** FIX ME **
  cgWave.dbase.get<int>("numPeriods") = numPeriods;          // ** FIX ME **
  cgWave.dbase.get<int>("adjustOmega")= adjustOmega;  // 1 : choose omega from the symbol of D+t D-t 

  // // >>> --- Setting these values in cgWave are now done in initialize()
  // const int & numberOfFrequencies = dbase.get<int>("numberOfFrequencies");
  // cgWave.dbase.get<int>("numberOfFrequencies") = numberOfFrequencies;

  // RealArray & cgWaveFrequencyArray = cgWave.dbase.get<RealArray>("frequencyArray");
  // cgWaveFrequencyArray.redim(numberOfFrequencies);
  // cgWaveFrequencyArray = dbase.get<RealArray>("frequencyArray");

  // RealArray & cgWavePeriodArray = cgWave.dbase.get<RealArray>("periodArray");
  // cgWavePeriodArray.redim(numberOfFrequencies);
  // cgWavePeriodArray = dbase.get<RealArray>("periodArray"); 
  // // <<<< ----



  real time0=getCPU(), timeb;

  // real timeForLaplace=0, timeForBoundaryConditions=0., timeForUpdateGhostBoundaries=0.,
  //   timeForInterpolate=0., timeForAdvance=0., timeForGetLocalArray=0.,
  //   timeForFinishBoundaryConditions=0.;
      
  printF("\n =========================  WAVE EQUATION HELMHOLTZ SOLVER ==========================\n");
  if( numberOfFrequencies==1 )
  {
    printF("CgWaveHoltz::solve using omega=%12.4e Tperiod=%12.4e (before adjust) numPeriods=%d, filterTimeDerivative=%d, numCompWaveHoltz=%d\n",
          omega,Tperiod,numPeriods,filterTimeDerivative, numCompWaveHoltz);
  }
  else
  {
    printF(" -------------------------- MULTI-FREQUENCY ALGORITHM -------------------------------\n");
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
      printF(" freq=%d : omega=%12.4e, T=%12.4e, numPeriods=%d\n",freq,frequencyArray(freq),periodArray(freq),numPeriodsArray(freq));
    }
  }

  // // Gaussian forcing 
  // const real & beta = dbase.get<real>("beta");
  // const real & x0   = dbase.get<real>("x0");
  // const real & y0   = dbase.get<real>("y0");
  // const real & z0   = dbase.get<real>("z0");
  // printF("Gaussian force: beta=%g x0=%g y0=%g z0=%g \n",beta,x0,y0,z0);


  // printF("CgWaveHoltz: c=%g, omega=%g, Tperiod=%g, numPeriods=%d, tFinal=%g, plotSteps=%d\n",
  // 	 c,omega,Tperiod,numPeriods,tFinal,plotSteps);

  realCompositeGridFunction & v    = cgWave.dbase.get<realCompositeGridFunction>("v");
  realCompositeGridFunction & vOld = cgWave.dbase.get<realCompositeGridFunction>("vOld");
  Range all;
  vOld.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz);


  if( false )
  {
    printF("CgWaveHoltz::solve\n");
    v.display("v at start of solve","%4.2f ");

  }

  // v=0; // Initial condition 

  int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");
  int & numberOfIterations        = dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken
  
  // Save "residuals" by iteration: 
  // resVector(it) = norm( v^{n+1} - v^n )
  RealArray & resVector = cgWave.dbase.get<RealArray>("resVector");
  resVector.redim(maximumNumberOfIterations+1);
  resVector=0.;

  int & plotOptions = cgWave.dbase.get<int>("plotOptions");
  plotOptions= CgWave::noPlotting; // turn of plotting in cgWave

  // ========== WaveHoltz ITERATIONS ===========
  for( int it=0; it<maximumNumberOfIterations; it++ ) 
  {
    vOld=v;  // save current solution 

    // -- advance for one period (or multiple periods ) ---
    // if( it==0 )
    //   printF("#################### CgWaveHoltz: CALL cgWave : it=%d  #########################\n",it);
    // else
    //  printF("#################### CgWaveHoltz: CALL cgWave : it=%d, L2h-norm(v-vOld)=%9.2e #########################\n",it,resVector(it-1));
    
    cgWave.advance( it );
    
    if( it>0 && omegaSOR != 1. )
    {
      v = (1.-omegaSOR)*vOld + omegaSOR*v;
    }
    
    vOld = v-vOld;

    // real errMax = maxNorm(vOld);
    Real errNorm;
    if( numCompWaveHoltz==1 )
    {
      errNorm = l2Norm(vOld);    // use L2h norm to match Krylov norm
    }
    else
    {
      // --- multiple frequencies or complex case with 2 components --- *fixed* May 9, 2025
      //   || v ||_2h^2 = SUM || v_c ||_2h^2
      // Note:
      //   l2Norm(v)^2 = SUM |v(i1,i2,i3)|^2 / (number of terms )
      errNorm=0.;
      for( int comp=0; comp<numCompWaveHoltz; comp++ )
      {
        Real errComp = l2Norm(vOld,comp);
        errNorm += SQR(errComp);
      }
      errNorm = sqrt(errNorm/numCompWaveHoltz);
    }



    printF("it=%d:  l2Norm( v-vOld )=%8.2e (tol=%g)\n",it,errNorm,tol);
    resVector(it)= errNorm; // norm( v-vOld )    

    if( it>0 )
      printF("#################### CgWaveHoltz Fixed-Point Iteration: it=%d, L2h-norm(v-vOld)=%9.2e, ratio=%5.2f #########################\n",
        it,resVector(it),resVector(it)/resVector(it-1));
    else
      printF("#################### CgWaveHoltz Fixed-Point Iteration: it=%d, L2h-norm(v-vOld)=%9.2e #########################\n",it,resVector(it));

    numberOfIterations=it+1;
    if( errNorm<tol )
      break;
    

  } // end for WaveHoltz Iteration
  if( numberOfIterations >= maximumNumberOfIterations )
  {
    printF("$$$$$$$$ CgWaveHoltz:ERROR -- number of iterations reached the maximum allowed = %d $$$$$$$$\n",maximumNumberOfIterations);
    
  }
  
  // --- Compute the average convergence rate ----
  Real & convergenceRate          = dbase.get<Real>("convergenceRate");
  Real & convergenceRatePerPeriod = dbase.get<Real>("convergenceRatePerPeriod");
  convergenceRate          = pow( resVector(numberOfIterations-1)/resVector(0), 1./( numberOfIterations ) ); 
  convergenceRatePerPeriod = pow( resVector(numberOfIterations-1)/resVector(0), 1./( numberOfIterations*numPeriods) ); 


  printF("###### DONE CgWaveHoltz: CALL cgWave : number of WaveHoltz iterations =%d ######\n",numberOfIterations);

  return 0;
}

// ================================================================================================
/// \brief Output the inital header for CgWaveHoltz with option and parameter values.
// ================================================================================================
int CgWaveHoltz::outputHeader()
{
  const aString & nameOfGridFile  = dbase.get<aString>("nameOfGridFile");
  const real & omega                    = dbase.get<real>("omega");
  const real & Tperiod                  = dbase.get<real>("Tperiod");
  const int & numPeriods                = dbase.get<int>("numPeriods");
  const real & tol                      = dbase.get<real>("tol");
  const int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");
  const int & useVariableTolerance      = dbase.get<int>("useVariableTolerance"); // Vary implicit solver tolerance based on current WaveHoltz residual

  const int & adjustOmega               = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

  const int & monitorResiduals          = dbase.get<int>("monitorResiduals");      // montior the residuals at every step
  const int & saveMatlabFile            = dbase.get<int>("saveMatabFile");         // save matlab file with residuals etc.
  const aString & matlabFileName        = dbase.get<aString>("matlabFileName");    // name of matlab file holding residuals etc.

  const real & omegaSOR  = dbase.get<real>("omegaSOR");

  FILE *& logFile = dbase.get<FILE*>("logFile");

  for( int fileio=0; fileio<2; fileio++ )
  {
    FILE *file = fileio==0 ? logFile : stdout; 
    fPrintF(file,"\n"
      "*********************************************************************************\n"
      "           CgWaveHoltz : Helmholtz Equation Solver                    \n"
      "           ---------------------------------------                  \n");

    fPrintF(file," Grid name=%s \n",(const char*)nameOfGridFile);
    fPrintF(file," omega=%14.6e, Tperiod=%14.6e \n",omega,Tperiod);
    fPrintF(file," adjustOmega=%d (1= adjust omega to account for discrete symbol of D+t D-t)\n",adjustOmega);
    fPrintF(file," saveMatlabFile=%d \n",saveMatlabFile);
    fPrintF(file," maximumNumberOfIterations=%d, useVariableTolerance=%d\n",maximumNumberOfIterations,useVariableTolerance);


    fPrintF(file,"*********************************************************************************\n\n");
    
  }


  return 0;
}


// ===================================================================================
/// \brief Save check file -- used for regression and convergence tests
/// \param checkFileCounter (input) : checkFileCounter=0,1,2,3 : a counter for different values saved in the check file.
/// \param maxErr (input) : some measure of the error.
/// \param solutionNorm (input) : norm of the solution.
// ===================================================================================
int CgWaveHoltz::saveCheckFile( int checkFileCounter, Real maxErr, Real solutionNorm )
{

  // ****** OLD WAY ****

  OV_ABORT("CgWaveHoltz::saveCheckFile CALLED -- but this is the old way\n");



  // *************************************************************
  // write to the check file for regression and convergence tests
  // *************************************************************

  FILE *& checkFile = dbase.get<FILE*>("checkFile");
  assert( checkFile != NULL );

  // const real & maxError     = dbase.get<real>("maxError");      // save max-error here 
  // const real & solutionNorm = dbase.get<real>("solutionNorm");  // save solution norm here 
  // const int & computeErrors = dbase.get<int>("computeErrors");
  
  const int & numberOfComponents= 1;

  int numberToOutput = numberOfComponents;

  Real t = checkFileCounter; 
  fPrintF(checkFile,"%9.2e %i  ",t,numberToOutput);
  for( int i=0; i<numberToOutput; i++ )
  {
    fPrintF(checkFile,"%i %9.2e %10.3e  ",i,maxErr,solutionNorm);
  }
  fPrintF(checkFile,"\n");


  return 0;
}

