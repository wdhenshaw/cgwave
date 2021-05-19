#include "CgWaveHoltz.h"
#include "CompositeGridOperators.h";	
#include "PlotStuff.h"
#include "display.h"
#include "ParallelOverlappingGridInterpolator.h"
#include "ParallelUtility.h"
#include "LoadBalancer.h"
#include "gridFunctionNorms.h"

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
  real & omega = dbase.put<real>("omega")=30.1;
  dbase.put<real>("Tperiod")=twoPi/omega;
  dbase.put<int>("numPeriods")=10;

  dbase.put<aString>("solverName")="fixedPoint";  // fixedPoint or gmres etc.
  dbase.put<real>("tol")=1.e-4;  // tolerance for Krylov solvers 

  dbase.put<aString>("nameOfGridFile")="unknown";

  dbase.put<int>("maximumNumberOfIterations")=500;
  dbase.put<int>("numberOfIterations")=0; // actual number of iterations taken

  dbase.put<int>("orderOfAccuracy")=0; 

  dbase.put<Real>("convergenceRate")=0.;
  dbase.put<Real>("convergenceRatePerPeriod")=0.;

  real & omegaSOR = dbase.put<real>("omegaSOR")=1.;

  dbase.put<int>("adjustOmega")=0;  // 1 : choose omega from the symbol of D+t D-t 

  dbase.put<int>("monitorResiduals")=1;  // montior the residuals at every step
  dbase.put<int>("saveMatabFile")=1;     // save matlab file with residuals etc.
  dbase.put<aString>("matlabFileName")="cgWaveHoltz.m";  // name of matlab file holding residuals etc.

  dbase.put<realCompositeGridFunction>("vOld");
  dbase.put<realCompositeGridFunction>("residual");

  // Save "residuals" by iteration: 
  // resVector(it) = norm( v^{n+1} - v^n )
  dbase.put<RealArray>("resVector");

  dbase.put<int>("petscIsInitialized")=false;

  FILE *& debugFile = dbase.put<FILE*>("debugFile");
  debugFile = fopen("cgWaveHoltz.debug","w" );        // log file 

  FILE *& logFile = dbase.put<FILE*>("logFile");
  logFile = fopen("cgWaveHoltz.out","w" );        // log file 

    FILE *& checkFile = dbase.put<FILE*>("checkFile");
  checkFile = fopen("cgWaveHoltz.check","w" );        // for regression and convergence tests

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

  real & omega     = dbase.get<real>("omega");
  real & Tperiod   = dbase.get<real>("Tperiod");
  int & numPeriods = dbase.get<int>("numPeriods");
  real & tol       = dbase.get<real>("tol");
  int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");

  int & adjustOmega = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

  int & monitorResiduals   = dbase.get<int>("monitorResiduals");      // montior the residuals at every step
  int & saveMatlabFile     = dbase.get<int>("saveMatabFile");         // save matlab file with residuals etc.
  aString & matlabFileName = dbase.get<aString>("matlabFileName");    // name of matlab file holding residuals etc.

  real & omegaSOR  = dbase.get<real>("omegaSOR");

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

    else if( dialog.getTextValue(answer,"omega","%e",omega) )
    {
      printF("Setting omega=%g\n",omega);
    }

    else if( dialog.getTextValue(answer,"tol","%e",tol) )
    {
      printF("Setting tol=%g (tolerence for Krylov solvers)\n",tol);
    }
    
    else if( dialog.getTextValue(answer,"number of periods","%i",numPeriods) )
    {
      printF("Setting numPeriods=%i\n",numPeriods);
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
int CgWaveHoltz::outputMatlabFile()
{
  
  const int myid=max(0,Communication_Manager::My_Process_Number);
  if( myid!=0 )
   return 0;

  const aString & solverName         = dbase.get<aString>("solverName");
  const aString & nameOfGridFile     = dbase.get<aString>("nameOfGridFile");
  const real & omega                 = dbase.get<real>("omega");
  const real & Tperiod               = dbase.get<real>("Tperiod");
  const int & numPeriods             = dbase.get<int>("numPeriods");
  const real & tol                   = dbase.get<real>("tol");
  const int & adjustOmega            = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

  RealArray & resVector              = dbase.get<RealArray>("resVector");
  const int & numberOfIterations     = dbase.get<int>("numberOfIterations");
  const Real & convergenceRate          = dbase.get<Real>("convergenceRate");
  const Real & convergenceRatePerPeriod = dbase.get<Real>("convergenceRatePerPeriod");

  CgWave & cgWave = *dbase.get<CgWave*>("cgWave");
  const int & orderOfAccuracy = cgWave.dbase.get<int>("orderOfAccuracy");

  // aString fileName="cgWaveHoltz.m"; // allow this to be specified
  aString & matlabFileName           = dbase.get<aString>("matlabFileName");    // name of matlab file holding residuals etc.

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

  fPrintF(matlabFile,"solverName='%s';\n",(const char*)solverName);
  fPrintF(matlabFile,"gridName=\'%s\';\n",(const char*)nameOfGridFile);
  fPrintF(matlabFile,"omega=%20.14e;\n",omega);
  fPrintF(matlabFile,"convergenceRate=%12.4e;\n",convergenceRate);
  fPrintF(matlabFile,"convergenceRatePerPeriod=%12.4e;\n",convergenceRatePerPeriod);
  fPrintF(matlabFile,"adjustOmega=%d; %% (1= adjust omega to account for discrete symbol of D+t D-t).\n",adjustOmega);
  fPrintF(matlabFile,"numPeriods=%d; %% (integrate over this many periods per Wave-Holtz iteration.\n",numPeriods);
  fPrintF(matlabFile,"orderOfAccuracy=%d;\n",orderOfAccuracy);
  fPrintF(matlabFile,"tol=%12.4e;\n",tol);

  int numPerLine=40;
  fPrintF(matlabFile,"%% itv = iteration number\n");
  fPrintF(matlabFile,"itv=[");
  for( int i=0; i<numberOfIterations; i++ )
  {
    fPrintF(matlabFile,"%i ",i);
    if( (i % numPerLine)==numPerLine-1 ) fprintf(matlabFile,"...\n");
  }
  fPrintF(matlabFile,"];\n");

  // ---- save the iteration history of the residual ----
  numPerLine=10;
  fPrintF(matlabFile,"res=[ ...\n");
  for( int i=0; i<numberOfIterations; i++ )
  {
    fPrintF(matlabFile,"%17.10e ",resVector(i));
    if( (i % numPerLine)==numPerLine-1 ) fprintf(matlabFile,"...\n");
  }
  fPrintF(matlabFile,"];\n");

  fclose(matlabFile);

  printf("CgWaveHoltz::Results saved to file %s\n",(const char*)fileName);


  return 0;
}

// ================================================================================================
/// \brief Setup grids and grid functions
// ================================================================================================
int CgWaveHoltz::setup()
{

  realCompositeGridFunction & vOld = dbase.get<realCompositeGridFunction>("vOld");
  Range all;
  vOld.updateToMatchGrid(cg,all,all,all);

  realCompositeGridFunction & residual = dbase.get<realCompositeGridFunction>("residual");
  residual.updateToMatchGrid(cg,all,all,all);

  // v=0;
  return 0;
}



// ================================================================================================
/// \brief Compute the residual in the current solution
// ================================================================================================
real CgWaveHoltz::residual()
{

  CgWave & cgWave                     = *dbase.get<CgWave*>("cgWave");
  const real & omega                  = dbase.get<real>("omega");
  const int & adjustOmega             = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

  const real & dt                     = cgWave.dbase.get<real>("dt");
  const real & c                      = cgWave.dbase.get<real>("c");
  realCompositeGridFunction & v       = cgWave.dbase.get<realCompositeGridFunction>("v");
  realCompositeGridFunction & f       = cgWave.dbase.get<realCompositeGridFunction>("f");
  CompositeGridOperators & operators  = cgWave.dbase.get<CompositeGridOperators>("operators");


  // Symbol of D+t D-t : 
  // D+tD-t exp(i*omega*t^n) = -4*sin^2(omega*dt/2)/dt^2 * exp(i*omega*t^n )
  // 
  //  omegaTilde = (2/dt)*sin(omega*dt/2)
  const Real omegas = (2./dt)*sin(omega*dt/2.);

  bool computeResidualUsingDiscreteSymbol = true;
  // bool adjustOmega = true;
  // const Real omegar = adjustOmega ? omegas : omega; // omega to use for residual

  printF("CgWaveHoltz::residual: c=%g, omega=%12.5e, omegas=%12.5e (from symbol of D+D-), dt=%12.6e, adjustOmega=%d\n",c,omega,omegas,dt,adjustOmega);

  realCompositeGridFunction & res     = dbase.get<realCompositeGridFunction>("residual");

  Real maxRes =0., maxResFromDiscreteSymbol=0.;

  Index I1,I2,I3;
  Index Ib1,Ib2,Ib3;

  int numRes=2; // compute the residual in two ways
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
        
        RealArray lap(I1,I2,I3);
        operators[grid].derivative(MappedGridOperators::laplacianOperator,vLocal,lap,I1,I2,I3);

        Real om = omega;
        if( ires==1 ) // computeResidualUsingDiscreteSymbol )
        {
          // --- Compute the residual using omegas (from discrete symbol) ---
          om = omegas;
        }

        // --- Compute the residual ---
        where( maskLocal(I1,I2,I3)>0 )
        {
          // resLocal(I1,I2,I3) = (c*c)*lap(I1,I2,I3) + (omega*omega)*vLocal(I1,I2,I3) + fLocal(I1,I2,I3); // change sign on f for cgWave **FIX ME**
          resLocal(I1,I2,I3) = (c*c)*lap(I1,I2,I3) + (om*om)*vLocal(I1,I2,I3) - fLocal(I1,I2,I3); 
        }

      }

      ForBoundary(side,axis)
      {
         // set residual to zero on dirichlet boundaries 
         if( mg.boundaryCondition(side,axis) == CgWave::dirichlet )
         {
           getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
           resLocal(Ib1,Ib2,Ib3)=0.;
         }

      }

      // ::display(res[grid],"residual","%8.2e ");
    }

    const int cc=0, maskOption=1;  // maskOption=1 : check points with mask>0
    if( ires==0 )
      maxRes                   = maxNorm(res,cc,maskOption);
    else
      maxResFromDiscreteSymbol = maxNorm(res,cc,maskOption);

  }


  
  if( computeResidualUsingDiscreteSymbol )
  {
    printF("CgWaveHoltz::residual: max-res=%9.3e (using omega), max-res=%9.3e (using omega from discrete symbol)\n",maxRes,maxResFromDiscreteSymbol);

  } 

  return maxRes;

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

 // here is the CgWave solver for the time dependent wave equation
 CgWave & cgWave = *dbase.get<CgWave*>("cgWave");
 Tperiod=numPeriods*twoPi/omega;  
 printF("CgWaveHoltz::solve: setting tFinal = Tperiod*numPeriods = %9.3e (numPeriods=%d) \n",Tperiod,numPeriods);
 
 // // --- set values in CgWave:
 // const real & dt    = cgWave.dbase.get<real>("dt");  // is this set ? .. NO **FIX ME***

 // Real omegaForWave   = omega;
 // Real TperiodForWave = Tperiod;
 // const Real omegaSymbol = (2./dt)*sin(omega*dt/2.);  // Symbol of D+tD-t 
 // if( adjustOmega )
 // {  // -- adjust omega to account for the Fourier Symbol of D+tD-t    ** THIS IS NOT CORRECT --> NEED TO DO INSIDE cgWave
 //    omegaForWave = omegaSymbol; 
 //    TperiodForWave = twoPi/omegaForWave; 
 // }

 cgWave.dbase.get<real>("omega")     = omega;        // ** FIX ME **
 cgWave.dbase.get<real>("tFinal")    = Tperiod;      // ** FIX ME **
 cgWave.dbase.get<real>("Tperiod")   = Tperiod;      // ** FIX ME **
 cgWave.dbase.get<int>("numPeriods") = numPeriods;          // ** FIX ME **
 cgWave.dbase.get<int>("adjustOmega")= adjustOmega;  // 1 : choose omega from the symbol of D+t D-t 

 


  real time0=getCPU(), timeb;

  // real timeForLaplace=0, timeForBoundaryConditions=0., timeForUpdateGhostBoundaries=0.,
  //   timeForInterpolate=0., timeForAdvance=0., timeForGetLocalArray=0.,
  //   timeForFinishBoundaryConditions=0.;
      
  printF("\n =========================  WAVE EQUATION HELMHOLTZ SOLVER ==========================\n"
         "CgWaveHoltz::solve using omega=%12.4e Tperiod=%12.4e (before adjust) numPeriods=%d\n",
          omega,Tperiod,numPeriods);


  // // Gaussian forcing 
  // const real & beta = dbase.get<real>("beta");
  // const real & x0   = dbase.get<real>("x0");
  // const real & y0   = dbase.get<real>("y0");
  // const real & z0   = dbase.get<real>("z0");
  // printF("Gaussian force: beta=%g x0=%g y0=%g z0=%g \n",beta,x0,y0,z0);


  // printF("CgWaveHoltz: c=%g, omega=%g, Tperiod=%g, numPeriods=%d, tFinal=%g, plotSteps=%d\n",
  // 	 c,omega,Tperiod,numPeriods,tFinal,plotSteps);

  realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
  realCompositeGridFunction & vOld = dbase.get<realCompositeGridFunction>("vOld");


  // v=0; // Initial condition 

  int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");
  int & numberOfIterations = dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken
  
  // Save "residuals" by iteration: 
  // resVector(it) = norm( v^{n+1} - v^n )
  RealArray & resVector = dbase.get<RealArray>("resVector");
  resVector.redim(maximumNumberOfIterations);
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

    real errMax = maxNorm(vOld);
    printF("it=%d:  max(|v-vOld|)=%8.2e (tol=%g)\n",it,errMax,tol);
    resVector(it)= errMax; // norm( v-vOld )     

    numberOfIterations=it+1;
    if( errMax<tol )
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


  printF("#################### DONE CgWaveHoltz: CALL cgWave : number of WaveHoltz iteration =%d #########################\n",numberOfIterations);

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




