#include "CgWaveHoltz.h"
#include "CompositeGridOperators.h";	
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
  IntegerArray & numPeriodsArray = dbase.put<IntegerArray>("numPeriodsArray");
  frequencyArray.redim(numberOfFrequencies);
  frequencyArray(0)=30.1;
  periodArray.redim(numberOfFrequencies);
  periodArray(0) = twoPi/frequencyArray(0);
  numPeriodsArray.redim(numberOfFrequencies);
  numPeriodsArray=1; 

  dbase.put<aString>("solverName")="fixedPoint";  // fixedPoint or gmres etc.
  dbase.put<real>("tol")=1.e-4;  // tolerance for Krylov solvers 

  dbase.put<aString>("nameOfGridFile")="unknown";

  dbase.put<int>("maximumNumberOfIterations")=500;
  dbase.put<int>("numberOfIterations")=0; // actual number of iterations taken
  dbase.put<int>("numberOfMatrixVectorMultiplications")=0;

  dbase.put<int>("orderOfAccuracy")=0; 

  dbase.put<Real>("convergenceRate")=0.;
  dbase.put<Real>("convergenceRatePerPeriod")=0.;

  dbase.put<Real>("maxResidual")=-1.;

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
  checkFile = fopen("cgWaveHoltz.check","w" );        // for regression and convergence tests
  fPrintF(checkFile,"# Check file for CgWaveHoltz\n"); // check file has one title line

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
  fclose(dbase.get<FILE*>("checkFile"));

  delete dbase.get<CgWave*>("cgWave");
}

// ================================================================================================
/// \brief Initialize time-step and forcing 
// ================================================================================================
int CgWaveHoltz::initialize()
{
  // ---- adjust periods ----
  const int & numPeriods          = dbase.get<int>("numPeriods");
  const int & numberOfFrequencies = dbase.get<int>("numberOfFrequencies");
  RealArray & frequencyArray      = dbase.get<RealArray>("frequencyArray");
  RealArray & periodArray         = dbase.get<RealArray>("periodArray");
  IntegerArray & numPeriodsArray  = dbase.get<IntegerArray>("numPeriodsArray");



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

  CgWave & cgWave = *dbase.get<CgWave*>("cgWave");
  cgWave.dbase.get<int>("numberOfFrequencies") = numberOfFrequencies;  

  RealArray & cgWaveFrequencyArray = cgWave.dbase.get<RealArray>("frequencyArray");
  cgWaveFrequencyArray.redim(numberOfFrequencies);
  cgWaveFrequencyArray = frequencyArray;

  RealArray & cgWavePeriodArray = cgWave.dbase.get<RealArray>("periodArray");
  cgWavePeriodArray.redim(numberOfFrequencies);
  cgWavePeriodArray = periodArray; 

  IntegerArray & cgWaveNumPeriodsArray = cgWave.dbase.get<IntegerArray>("numPeriodsArray");
  cgWaveNumPeriodsArray.redim(numberOfFrequencies);
  cgWaveNumPeriodsArray = numPeriodsArray;


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

  int & numberOfFrequencies       = dbase.get<int>("numberOfFrequencies");
  RealArray & frequencyArray      = dbase.get<RealArray>("frequencyArray");
  RealArray & periodArray         = dbase.get<RealArray>("periodArray");


  int & adjustOmega               = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

  int & monitorResiduals          = dbase.get<int>("monitorResiduals");      // montior the residuals at every step
  int & saveMatlabFile            = dbase.get<int>("saveMatabFile");         // save matlab file with residuals etc.
  aString & matlabFileName        = dbase.get<aString>("matlabFileName");    // name of matlab file holding residuals etc.

  real & omegaSOR                 = dbase.get<real>("omegaSOR");

  CgWave & cgWave                 = *dbase.get<CgWave*>("cgWave");

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
                           ""};
  int tbState[10];
  tbState[0] = saveMatlabFile;
  tbState[1] = monitorResiduals;
  tbState[2] = adjustOmega;
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

    else
    {
      printF("CgWaveHoltz:ERROR: unknown answer=[%s]\n",(const char*)answer);
    }
    
  }
  
  gi.popGUI();  // pop dialog

  // Initialize time-step and forcing 
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


  const int & numberOfIterations           = dbase.get<int>("numberOfIterations");
  const Real & convergenceRate             = dbase.get<Real>("convergenceRate");
  const Real & convergenceRatePerPeriod    = dbase.get<Real>("convergenceRatePerPeriod");


  CgWave & cgWave                                           = *dbase.get<CgWave*>("cgWave");

  const real & omega                                        = cgWave.dbase.get<real>("omega");
  const real & Tperiod                                      = cgWave.dbase.get<real>("Tperiod");
  const int & numPeriods                                    = cgWave.dbase.get<int>("numPeriods");
  const int & adjustOmega                                   = cgWave.dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  const int & orderOfAccuracy                               = cgWave.dbase.get<int>("orderOfAccuracy");
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


  // aString fileName="cgWaveHoltz.m"; // allow this to be specified
  // aString & matlabFileName           = dbase.get<aString>("matlabFileName");    // name of matlab file holding residuals etc.

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
  fPrintF(matlabFile,"omega=%20.14e;\n",omega);

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

  fPrintF(matlabFile,"computeEigenmodes=%d;\n",computeEigenmodes);
  if( computeEigenmodes )
  {

    CgWave::EigenSolverEnum & eigenSolver = cgWave.dbase.get<CgWave::EigenSolverEnum>("eigenSolver"); 
    const bool useFixedPoint =  eigenSolver==CgWave::fixedPointEigenSolver;

    aString eigenSolverName = eigenSolver==CgWave::defaultEigenSolver     ? "KrylovSchur" :
                              eigenSolver==CgWave::KrylovSchurEigenSolver ? "KrylovSchur" :
                              eigenSolver==CgWave::ArnoldiEigenSolver     ? "Arnoldi" : 
                              eigenSolver==CgWave::ArpackEigenSolver      ? "Arpack" : 
                              eigenSolver==CgWave::fixedPointEigenSolver  ? "fixedPoint" : 
                                                                            "unknown";
    fPrintF(matlabFile,"eigenSolver='%s';\n",(const char*)eigenSolverName);

    const int & numEigenVectors = dbase.get<int>("numEigenVectors");
    fPrintF(matlabFile,"numEigenVectors=%d; %% number of computed eigenpairs\n",numEigenVectors);

    const RealArray & eigenValues = dbase.get<RealArray>("eigenValues");
    fPrintF(matlabFile,"%% computed eigenvalues:\n");
    fPrintF(matlabFile,"eigenValues=[%24.16e",eigenValues(0));
    for( int ie=1; ie<numEigenVectors; ie++ )
    {
      fPrintF(matlabFile,", %24.16e",eigenValues(ie));
      if( (ie % numPerLine)==numPerLine-1 ) fPrintF(matlabFile,"...\n");
    }
    fPrintF(matlabFile,"];\n");

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
  fPrintF(matlabFile,"tol=%12.4e;\n",tol);
  fPrintF(matlabFile,"timeSteppingMethod='%s';\n",(timeSteppingMethod==CgWave::explicitTimeStepping ? "explicit" : "implicit") );

  printF("saveMatlab: dt=%9.3e\n",dt);
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
/// \brief Setup grids and grid functions
// ================================================================================================
int CgWaveHoltz::setup()
{
  const int & numberOfFrequencies = dbase.get<int>("numberOfFrequencies");
  printF("cgWaveHoltz::setup: numberOfFrequencies=%d\n",numberOfFrequencies);

  // realCompositeGridFunction & vOld = dbase.get<realCompositeGridFunction>("vOld");
  // Range all;
  // vOld.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);

  realCompositeGridFunction & residual = dbase.get<realCompositeGridFunction>("residual");
  Range all;
  residual.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);
  for( int freq=0; freq<numberOfFrequencies; freq++ )
  {
    residual.setName(sPrintF("r%d",freq),freq);
  }

  // v=0;
  return 0;
}



// ================================================================================================
/// \brief Compute the residual in the current solution
/// \param useAdjustedOmega (input) : 0 = use standard omega
///                                   1 = compute residual using the adjusted omega, 
///                                   2 = compute resdiual for both omega and adjusted omega.
// ================================================================================================
real CgWaveHoltz::residual( int useAdjustedOmega /* = 2 */ )
{

  const int & debug                       = dbase.get<int>("debug");

  CgWave & cgWave                         = *dbase.get<CgWave*>("cgWave");
  const real & omega                      = cgWave.dbase.get<real>("omega");
  const int & adjustOmega                 = cgWave.dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
    
  const real & dt                         = cgWave.dbase.get<real>("dt");
  const real & c                          = cgWave.dbase.get<real>("c");
  const int & upwind                      = cgWave.dbase.get<int>("upwind");
  const int & adjustHelmholtzForUpwinding = cgWave.dbase.get<int>("adjustHelmholtzForUpwinding");
  const int & computeEigenmodes           = cgWave.dbase.get<int>("computeEigenmodes");

  const int & numberOfFrequencies         = cgWave.dbase.get<int>("numberOfFrequencies");
  const RealArray & frequencyArray        = cgWave.dbase.get<RealArray>("frequencyArray");
    
  realCompositeGridFunction & v           = cgWave.dbase.get<realCompositeGridFunction>("v");
  realCompositeGridFunction & f           = cgWave.dbase.get<realCompositeGridFunction>("f");
  CompositeGridOperators & operators      = cgWave.dbase.get<CompositeGridOperators>("operators");


  // Symbol of D+t D-t : 
  // D+tD-t exp(i*omega*t^n) = -4*sin^2(omega*dt/2)/dt^2 * exp(i*omega*t^n )
  // 
  //  omegaTilde = (2/dt)*sin(omega*dt/2)
  const Real omegas = (2./dt)*sin(omega*dt/2.);

  // bool computeResidualUsingDiscreteSymbol = true;
  // bool adjustOmega = true;
  // const Real omegar = adjustOmega ? omegas : omega; // omega to use for residual

  // printF("\n +++++ ENTERING CgWaveHoltz::residual +++++ \n");

  // printF("CgWaveHoltz::residual: c=%g, omega=%14.7e, omegas=%14.7e (from symbol of D+D-), dt=%12.6e, adjustOmega=%d, upwind=%d\n",c,omega,omegas,dt,adjustOmega,upwind);
  if( upwind && !adjustHelmholtzForUpwinding )
  {
    printF("\n **** WARNING: upwinding is ON! Correcting omega for the discrete-time-symbol will not fully work to give a small residual ! *****\n\n");
  }
  realCompositeGridFunction & res = dbase.get<realCompositeGridFunction>("residual");

  RealArray maxRes(numberOfFrequencies), maxResFromDiscreteSymbol(numberOfFrequencies);
  maxRes=0.;
  maxResFromDiscreteSymbol=0.;

  Index I1,I2,I3;
  Index Ib1,Ib2,Ib3;

  // If we have not adjusted omega, then final residual is with standard omega
  // If we have adjusted omega, then final residual is with the adjusted
  int numRes=1;  // number of residuals we compute 
  bool computeResidualWithAdjustedOmega[2] ={ false,false };

  if( useAdjustedOmega==0 || numberOfFrequencies>1 || computeEigenmodes )
  {
    numRes=1;  computeResidualWithAdjustedOmega[0]=false;
  }
  else if( useAdjustedOmega==1 )
  {
    numRes=1;  computeResidualWithAdjustedOmega[0]=true;
  }
  else
  {
    numRes=2; // compute redisuals in two ways
    if( adjustOmega )
    {
      computeResidualWithAdjustedOmega[0]=false;
      computeResidualWithAdjustedOmega[1]=true;
    }
    else
    {
      computeResidualWithAdjustedOmega[0]=true;
      computeResidualWithAdjustedOmega[1]=false;
    }
  }

  for( int ires=0; ires<numRes; ires++ )
  {

    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg = cg[grid];

      OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

      OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
      OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);
      OV_GET_SERIAL_ARRAY(real,res[grid],resLocal);

      int extra=0; // -1;
      getIndex(cg[grid].gridIndexRange(),I1,I2,I3,extra); 
      
      bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
      if( ok )
      {
        resLocal=0.;
        
        RealArray lap(I1,I2,I3,numberOfFrequencies);
        operators[grid].derivative(MappedGridOperators::laplacianOperator,vLocal,lap,I1,I2,I3);

        for( int freq=0; freq<numberOfFrequencies; freq++ )
        {
          // Real om = omega;
          Real om = frequencyArray(freq);
          if( computeResidualWithAdjustedOmega[ires] )
          {
            // --- Compute the residual using omegas (from discrete symbol) ---
            if( grid==0 && debug & 1 )
              printF("CgWaveHoltz::residual: compute residual with adjusted omegas=%20.12e\n",omegas);
            om = omegas;
          }

          // --- Compute the residual ---
          if( computeEigenmodes )
          {
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
              if( maskLocal(i1,i2,i3)>0 )
              {
                // relative residual : scale by 1/omega^2
                resLocal(i1,i2,i3,freq) = ( (c*c)/(om*om) )*lap(i1,i2,i3,freq) + vLocal(i1,i2,i3,freq);
              }
            }
          }
          else
          {
            where( maskLocal(I1,I2,I3)>0 )
            {
              // resLocal(I1,I2,I3) = (c*c)*lap(I1,I2,I3) + (omega*omega)*vLocal(I1,I2,I3) + fLocal(I1,I2,I3); // change sign on f for cgWave **FIX ME**
              resLocal(I1,I2,I3,freq) = (c*c)*lap(I1,I2,I3,freq) + (om*om)*vLocal(I1,I2,I3,freq) - fLocal(I1,I2,I3,freq); 
            }
          }
        }
      }

      Range all;
      ForBoundary(side,axis)
      {
         // set residual to zero on dirichlet boundaries 
         if( mg.boundaryCondition(side,axis) == CgWave::dirichlet )
         {
           getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
           resLocal(Ib1,Ib2,Ib3,all)=0.;
         }

      }

      // ::display(res[grid],"residual","%8.2e ");
    }

    const int maskOption=1;  // maskOption=1 : check points with mask>0
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
      if( computeResidualWithAdjustedOmega[ires]  )
        maxResFromDiscreteSymbol(freq) = maxNorm(res,freq,maskOption);
      else
        maxRes(freq)                   = maxNorm(res,freq,maskOption);
    }
  }


  for( int freq=0; freq<numberOfFrequencies; freq++ )
  {
    if( computeEigenmodes )
    {
       printF("CgWaveHoltz::residual: freq=%2d, omega=%8.3f, max-rel-res=%9.3e.\n",freq,frequencyArray(freq),maxRes(freq));
    }  
    else if( useAdjustedOmega==0 || numberOfFrequencies>1 || computeEigenmodes )
    {
       printF("CgWaveHoltz::residual: freq=%2d, omega=%8.3f, max-res=%9.3e.\n",freq,frequencyArray(freq),maxRes(freq));
    }        
    else if( useAdjustedOmega==2 )
    {
      printF("CgWaveHoltz::residual: freq=%2d, omega=%8.3f, max-res=%9.3e (using omega), max-res=%9.3e (using omega from discrete symbol)\n",
               freq,frequencyArray(freq),maxRes(freq),maxResFromDiscreteSymbol(freq));
    } 
    else
    {
      printF("CgWaveHoltz::residual: freq=%d, omega=%g, max-res=%9.3e (using omega from discrete symbol)\n",freq,frequencyArray(freq),maxResFromDiscreteSymbol(freq));
    }
  }

  return max(maxRes);

}


// ================================================================================================
/// \brief Solve the equations over one or more periods 
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
  IntegerArray & numPeriodsArray  = dbase.get<IntegerArray>("numPeriodsArray");


  // here is the CgWave solver for the time dependent wave equation
  CgWave & cgWave = *dbase.get<CgWave*>("cgWave");
  if( omega!=0.0 )
    Tperiod=numPeriods*twoPi/omega;  
  else
    Tperiod=1;

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
    printF("CgWaveHoltz::solve using omega=%12.4e Tperiod=%12.4e (before adjust) numPeriods=%d\n",
          omega,Tperiod,numPeriods);
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
  vOld.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);



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
    printF("#################### CgWaveHoltz: CALL cgWave : it=%d #########################\n",it);
    
    cgWave.advance( it );
    
    if( it>0 && omegaSOR != 1. )
    {
      v = (1.-omegaSOR)*vOld + omegaSOR*v;
    }
    
    vOld = v-vOld;

    // real errMax = maxNorm(vOld);
    real errNorm = l2Norm(vOld);    // use L2h norm to match Krylov norm
    printF("it=%d:  l2Norm( v-vOld )=%8.2e (tol=%g)\n",it,errNorm,tol);
    resVector(it)= errNorm; // norm( v-vOld )     

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


  printF("###### DONE CgWaveHoltz: CALL cgWave : number of WaveHoltz iteration =%d ######\n",numberOfIterations);

  return 0;
}

// ================================================================================================
/// \brief Output the inital header for CgWaveHoltz with option and parameter values.
// ================================================================================================
int CgWaveHoltz::outputHeader()
{
  const aString & nameOfGridFile  = dbase.get<aString>("nameOfGridFile");
  real & omega                    = dbase.get<real>("omega");
  real & Tperiod                  = dbase.get<real>("Tperiod");
  int & numPeriods                = dbase.get<int>("numPeriods");
  real & tol                      = dbase.get<real>("tol");
  int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");

  int & adjustOmega               = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

  int & monitorResiduals          = dbase.get<int>("monitorResiduals");      // montior the residuals at every step
  int & saveMatlabFile            = dbase.get<int>("saveMatabFile");         // save matlab file with residuals etc.
  aString & matlabFileName        = dbase.get<aString>("matlabFileName");    // name of matlab file holding residuals etc.

  real & omegaSOR  = dbase.get<real>("omegaSOR");

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

