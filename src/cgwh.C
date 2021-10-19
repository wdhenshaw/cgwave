// ---------------------------------------------------------------------------
//
//     cgWaveHoltz -grid=square128.order2
//     cgWaveHoltz -grid=cice8.order2    
//
//      **** Parallel version ****
//
// Solve the second order wave equation
//       u_tt - c^2 ( u_xx + u_yy ) =0 
// Notes:
//    (1) Solve to 2nd or fourth order. If the overlapping grid is made with the 
//        option "order of accuracy" set to "fourth order" then this will be a true
//        fourth order scheme. This program will also work with the fourth order
//        approximation even if the grid is made for second order accuracy.
//    (2) Add a fourth order artificial dissipation of the form 
//                     ad4 h^4 dt (u.xxxx).t
//        This artificial diffusion may not even be needed in most cases.
//    (3) We try to be a little careful with memory usage here. Rectangular grids will
//        be efficiently treated -- only the mask array will be built for a rectangular
//        grid, even the array of verticies will not need to be built.
//   
//
// Examples for running in parallel
//
//  mpirun -np 4 pwave -grid=cice16.order4.hdf
//  mpirun -np 4 pwave -noplot -g=cice64.order4.hdf -cmd=pwave1.cmd
//
//  mpirun -np 2 -all-local pwave
//
// ---------------------------------------------------------------------------

#include "CgWaveHoltz.h"
#include "CgWave.h"

// #include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"
#include "PlotStuff.h"
#include "display.h"
#include "ParallelOverlappingGridInterpolator.h"
#include "ParallelUtility.h"
#include "LoadBalancer.h"
#include "gridFunctionNorms.h"



int 
getLineFromFile( FILE *file, char s[], int lim);

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  


  // Pulse parameters:
static real alpha=50.; // 200.;
static real pulsePow=10.; // 20
static real a0=3.;
static real xPulse=-1.2;
static real yPulse=0.;

enum InitialConditionOptionEnum
{
  smoothPulse,
  pulse
};



//=============================================================================================================
//     CgWaveHoltz
//         Compute solutions to the Helmholtz equation by solving the time dependent wave equation
// 
//  Ref. Daniel Appelo -- WaveHoltz
//=============================================================================================================
int 
main(int argc, char *argv[])
{
  int debug=0;

  Overture::start(argc,argv);  // initialize Overture
  const int myid=Communication_Manager::My_Process_Number;
  const int np = Communication_Manager::numberOfProcessors();
  
  //     cgWaveHoltz -grid=square64.order2  

  printF("Usage: `mpirun -np N cgWaveHoltz [-noplot] [file.cmd] [-g=<gridName>]'\n");
  
  // Use this to avoid un-necessary communication: 
  Optimization_Manager::setForceVSG_Update(Off);


  aString nameOfOGFile= "square32.order2.hdf"; 

  aString commandFileName="";
  aString line;
  int len=0;
  bool plotOption=true;  // by default we plot interactively
  if( argc > 1 )
  {
    int i;
    for( i=1; i<argc; i++ )
    {
      line=argv[i];
      if( line=="-noplot" || line=="noplot" )
        plotOption=false; 
      else if( len=line.matches("-g=") )
      {
        nameOfOGFile=line(len,line.length()-1);
        // printf("\n$$$$ node %i : use grid=[%s]\n",myid,(const char*)nameOfOGFile);
      }
      else if( commandFileName == "" )
      {
        commandFileName=line(len,line.length()-1);
        // printf("\n$$$$ node %i : read command file %s\n",myid,(const char*)commandFileName);
      }
    }
  }

  
  GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("cgwh",plotOption,argc,argv));

  PlotStuffParameters psp;

  // By default start saving the command file called "cgWaveHoltz.cmd"
  aString logFile="cgWaveHoltz.cmd";
  ps.saveCommandFile(logFile);
  printF("User commands are being saved in the file `%s'\n",(const char *)logFile);

  if( commandFileName!="" )
    ps.readCommandFile(commandFileName);

  // int orderOfAccuracy=4;  // **** need two ghost lines in parallel ****

  aString nameOfShowFile = "cgWaveHoltz.show";

  // create and read in a CompositeGrid
  #ifdef USE_PPP
    // On Parallel machines always add at least this many ghost lines on local arrays
    const int numGhost=2;
    MappedGrid::setMinimumNumberOfDistributedGhostLines(numGhost);
  #endif
  CompositeGrid cg;
  bool loadBalance=true; // turn on or off the load balancer
  getFromADataBase(cg,nameOfOGFile,loadBalance);

  cg.update(MappedGrid::THEmask);
  const int numberOfDimensions = cg.numberOfDimensions();


  // real omega =  30.1; // 15.3; // 7.1; // 3.25; // 7.1; // 3.25; // 25.1; // 3.25; // Helmholtz frequency 
  // int numPeriods=10;  // 1 // integrate this many periods
  // real Tperiod = numPeriods*twoPi/omega;

  // real tFinal=Tperiod;
  // real tPlot=Tperiod/10.;  // plot this often
  // real tShow=Tperiod/10.;  // save to show file this often
  // bool saveShowFile=false;
  bool showFileIsOpen=false;
  
  // real t=0;
  // real c=1., cSquared=c*c;
    
  // real ad4=1.;  // coeff of the artificial dissipation.

  // real cfl=.25, nu=0.;
  // InitialConditionOptionEnum option=smoothPulse;


  // ---- Build the WaveHoltz solver and set parameters ----
  CgWaveHoltz cgWaveHoltz(cg,ps);
  cgWaveHoltz.setNameOfGridFile(nameOfOGFile);
  cgWaveHoltz.setup();
  cgWaveHoltz.interactiveUpdate();

  Real & convergenceRate          = cgWaveHoltz.dbase.get<Real>("convergenceRate");
  Real & convergenceRatePerPeriod = cgWaveHoltz.dbase.get<Real>("convergenceRatePerPeriod");

  // put there here for now 
  // here is the CgWave solver for the time dependent wave equation
  CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
  cgWave.setNameOfGridFile(nameOfOGFile);
  cgWave.setup();
  cgWave.interactiveUpdate();

  const int & orderOfAccuracy       = cgWave.dbase.get<int>("orderOfAccuracy");
  const int & orderOfAccuracyInTime = cgWave.dbase.get<int>("orderOfAccuracyInTime");
  const int & upwind                = cgWave.dbase.get<int>("upwind"); // new way
  // const Real & ad4                  = cgWave.dbase.get<Real>("ad4");   // old way
  const int & computeErrors         = cgWave.dbase.get<int>("computeErrors");

  cgWaveHoltz.dbase.get<int>("orderOfAccuracy")=orderOfAccuracy; // set value in CgWaveHoltz

  real & omega                    = cgWaveHoltz.dbase.get<real>("omega");
  int & numPeriods                = cgWaveHoltz.dbase.get<int>("numPeriods");
  real & tol                      = cgWaveHoltz.dbase.get<real>("tol");
  int & adjustOmega               = cgWaveHoltz.dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  int & maximumNumberOfIterations = cgWaveHoltz.dbase.get<int>("maximumNumberOfIterations");
  int & cgWaveDebugMode           = cgWaveHoltz.dbase.get<int>("cgWaveDebugMode");
  // Build a dialog menu for changing parameters
  GUIState gui;
  DialogData & dialog=gui;

  dialog.setWindowTitle("CgWaveHoltz - Helmholtz Solver");
  dialog.setExitCommand("exit", "exit");

  // dialog.setOptionMenuColumns(1);

  // aString accuracyLabel[] = {"second order", "fourth order", "" };
  // dialog.addOptionMenu("accuracy:", accuracyLabel, accuracyLabel, (orderOfAccuracy==2 ? 0 : 1) );

  // aString initialConditionLabel[] = {"smooth pulse", "pulse", "" };
  // dialog.addOptionMenu("Initial Condition:",initialConditionLabel,initialConditionLabel,(int)option );

  // enum TypeOfApproximation
  // {
  //   nonConservative=0,
  //   conservative=1
  // } typeOfApproximation=nonConservative;
    
  // aString approximationLabel[] = {"nonconservative", "conservative", "" };
  // dialog.addOptionMenu("Approximation:",approximationLabel,approximationLabel,(int)typeOfApproximation );

  aString pbLabels[] = {
                        "compute with fixed-point",
                        "compute with petsc",
                        "solve Helmholtz directly",
                        "zero initial condition",
                        "random initial condition",
                        "change parameters",
                        "contour",
                        "grid",
                        "plot residual",
                        // "compute errors",
                        "run cgWave and plot",
                        "save to show",
                        "plot errors",
                        "plot forcing",
                        "erase",
                        "exit",
                        ""};
  int numRows=8;
  dialog.setPushButtons( pbLabels, pbLabels, numRows ); 

  // aString tbCommands[] = {"save show file",
  //                          ""};
  // int tbState[10];
  // tbState[0] = saveShowFile==true;
  // int numColumns=1;
  // dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 

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

  textCommands[nt] = "max iterations";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",maximumNumberOfIterations);  nt++; 

  textCommands[nt] = "show file";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%s",(const char*)nameOfShowFile);  nt++; 

  // null strings terminal list
  textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
  dialog.setTextBoxes(textCommands, textLabels, textStrings);

  
  aString tbCommands[] = {
                          "run cgWave with debugging",
                            ""};
  int tbState[10];
  tbState[0] = cgWaveDebugMode;

  int numColumns=1;
  dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 


  ps.pushGUI(gui);



  psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

  aString answer;
  char buff[200];

  // keep track of what has been plotted
  int plotChoices=0; // 1=grid, 2=contour
  bool replot=false;
  bool reComputeErrors=false; 
  int checkFileCounter=0; // keeps track of how many times the checkFile is saved with results

  // uHelmholtz holds Helmholtz solution from direct solver
  //   This is optionally used in CgWave as a known solution
  if( !cgWave.dbase.has_key("uHelmholtz") )
    cgWave.dbase.put<realCompositeGridFunction>("uHelmholtz");

  realCompositeGridFunction & uHelmholtz = cgWave.dbase.get<realCompositeGridFunction>("uHelmholtz"); 
  uHelmholtz.updateToMatchGrid(cg); // save Helmholtz solution here
  bool helmholtzFromDirectSolverWasComputed=false;

  Ogshow show;  // show file for saving solutions

  int current=0;
  for(;;) 
  {
    ps.getAnswer(answer,"");      
    if( answer=="exit" || answer=="continue" )
    {
      break;
    }
    else if( answer=="change parameters" )
    {
      cgWaveHoltz.interactiveUpdate();
    }
    
    else if( dialog.getTextValue(answer,"omega","%e",omega) )
    {
      printF("Setting omega=%g\n",omega);
    }
    else if( dialog.getTextValue(answer,"number of periods","%i",numPeriods) )
    {
      printF("Setting numPeriods=%i\n",numPeriods);
    }
    else if( dialog.getTextValue(answer,"tol","%e",tol) )
    {
      printF("Setting tol=%g (tolerence for Krylov solvers)\n",tol);
    }
    else if( dialog.getTextValue(answer,"max iterations","%i",maximumNumberOfIterations) )
    {
      printF("Setting maximumNumberOfIterations=%i\n",maximumNumberOfIterations);
    }   
    else if( dialog.getTextValue(answer,"show file","%s",nameOfShowFile) )
    {
      printF("Setting nameOfShowFile=%s\n",(const char*)nameOfShowFile);
      if( showFileIsOpen )
      {
        // close any open show file
        show.close();
        showFileIsOpen=false;
      }
    }   
    else if( dialog.getToggleValue(answer,"run cgWave with debugging",cgWaveDebugMode) )
    {
      printF("Setting cgWaveDebugMode=%d (1: run cgWave plotting every step)\n",cgWaveDebugMode);
      if( cgWaveDebugMode==1 )
      {

      }
    }          

    else if( answer=="compute with fixed-point" )
    {
      const Real cpu0=getCPU();

      cgWaveHoltz.solve();

      const Real cpu = getCPU()-cpu0;
      real maxRes = cgWaveHoltz.residual();
      printF("CgWaveHoltz: omega=%9.3e, max-res=%9.3e, CR=%5.3f CR-perPeriod=%5.3f, cpu=%9.2e(s)\n",
              omega,maxRes,convergenceRate,convergenceRatePerPeriod,cpu);

      // save results to a matlab file
      cgWaveHoltz.dbase.get<aString>("solverName")="fixedPoint";
      cgWaveHoltz.outputMatlabFile();

      replot=true;
      reComputeErrors=true;

    }
    else if( answer=="compute with petsc" )
    {
      // Krylov solver 
      const Real cpu0=getCPU();

      cgWaveHoltz.solvePETSc(argc,argv);

      const Real cpu = getCPU()-cpu0;

      real maxRes = cgWaveHoltz.residual();
      printF("CgWaveHoltz: omega=%9.3e, max-res=%9.3e, CR=%5.3f CR-perPeriod=%5.3f, cpu=%9.2e(s)\n",
             omega,maxRes,convergenceRate,convergenceRatePerPeriod,cpu);

      // save results to a matlab file 
      cgWaveHoltz.dbase.get<aString>("solverName")="gmres";  // ** fix me **
      cgWaveHoltz.outputMatlabFile();

      replot=true;
      reComputeErrors=true;
    }
    else if( answer == "solve Helmholtz directly" )
    {
       // this may require cgWave to have been called to create f ?
      // uHelmholtz.updateToMatchGrid(cg); // save Helmholtz solution here

      realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
      realCompositeGridFunction & f = cgWave.dbase.get<realCompositeGridFunction>("f");

      // Make a copy of f since the direct solver will change it 
      realCompositeGridFunction fHelmholtz(cg);
      fHelmholtz = f; 

      const Real cpu0=getCPU();

      cgWaveHoltz.solveHelmholtz( v, fHelmholtz  );

      const Real cpu = getCPU()-cpu0;

      uHelmholtz = v;
      helmholtzFromDirectSolverWasComputed=true;

      if( true )
      {
        bool useAdjustedOmega=false;
        Real maxRes = cgWaveHoltz.residual( useAdjustedOmega );
        printF("CgWaveHoltz: omega=%9.3e, max-res=%9.3e, cpu=%9.2e (Direct solution of Helmholtz).\n",
               omega,maxRes,cpu);
  
        maxRes = cgWaveHoltz.residual();
        printF("CgWaveHoltz: omega=%9.3e, max-res=%9.3e (direct Helmholtz solution)\n",omega,maxRes);
      }

      reComputeErrors=true;
      replot=true;

    }
    else if( answer=="compute errors" )
    {
      if( computeErrors )
        reComputeErrors=true; // errors are computed below
      else
        printF("compute errors: WARNING: this solution has no errors to compute. Probably not a known solution.\n"); 
    }
    else if( answer == "zero initial condition" )
    {
      realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
      v=0.; 
    }
    else if( answer == "random initial condition" )
    {
      realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
      CompositeGrid & cg = *v.getCompositeGrid();
      Index I1,I2,I3;

      std::srand(12789.);

      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        MappedGrid & mg = cg[grid];
        OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
        getIndex(mg.gridIndexRange(),I1,I2,I3);
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
          // do this for now
          // vLocal(i1,i2,i3) = sin(i1)*cos(i2);

          vLocal(i1,i2,i3) = -1. + std::rand()*(2./RAND_MAX); // [-1,1]

        }

      }
      replot=true;

    }    
    else if( answer=="save to show" )
    {
      printF("Save the current solution (and errors) to the show file [%s].\n",(const char*)nameOfShowFile);
      if( !showFileIsOpen )
      {
        showFileIsOpen=true;
        show.open(nameOfShowFile);
      }
      
      show.saveGeneralComment("Solutions from CgWaveHoltz"); // save a general comment in the show file
      // show.saveGeneralComment(" file written on April 1");      // save another general comment

      show.startFrame();                       // start a new frame
      // const bool useUpwind = ad4>0.;
      show.saveComment(0,sPrintF("CgWaveHoltz: FD%i%i%s omega=%.5g",orderOfAccuracy,orderOfAccuracyInTime,(upwind ? "u" : ""),omega));   // comment 0 (shown on plot)
      show.saveComment(1,sPrintF("v=WaveHoltz, u=Helmholtz"));               // comment 1 (shown on plot)

      // We save the current solution and optionally the direct solution and or errors
      int numShowComponents=1;
      if( helmholtzFromDirectSolverWasComputed ) numShowComponents++;
      if( computeErrors || helmholtzFromDirectSolverWasComputed ) numShowComponents++;

      Range all;
      realCompositeGridFunction q(cg,all,all,all,numShowComponents);
      q.setName("q");
      int ishow=0;  // counts components in q as we fill them in

      q.setName("v",ishow); 
      realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");

      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
         q[grid](all,all,all,ishow) = v[grid];   // waveHoltz solution

      ishow++;

      if( helmholtzFromDirectSolverWasComputed )
      { 
        q.setName("u",ishow); 
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
          q[grid](all,all,all,ishow) = uHelmholtz[grid];  // direct solver solution
        ishow++; 
      }
      if( computeErrors || helmholtzFromDirectSolverWasComputed )
      { 
        realCompositeGridFunction & error = cgWave.dbase.get<realCompositeGridFunction>("error");
        if( !computeErrors && helmholtzFromDirectSolverWasComputed )
        {
          // compute difference between WaveHoltz and Direct Helmholtz solve
          realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
          error = uHelmholtz - v;
        }
        q.setName("error",ishow);
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
          q[grid](all,all,all,ishow) = error[grid];  // error or diff with direct
        ishow++;
      }
      assert( ishow==numShowComponents );


      show.saveSolution( q );              // save to show file
      show.endFrame();

      // show.close();
    }

    else if( answer=="contour" )
    {
      // -- this is done below so that we can also treat the case of replot=true ---
    }
    else if( answer=="grid" )
    {
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      PlotIt::plot(ps,cg,psp);                          // plot the grid
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

    }
    else if( answer=="plot residual" )
    {
      realCompositeGridFunction & res = cgWaveHoltz.dbase.get<realCompositeGridFunction>("residual"); 
      res.interpolate();
           
      ps.erase();
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      psp.set(GI_TOP_LABEL,sPrintF("residual O%d omega=%.5g",orderOfAccuracy,omega));
      PlotIt::contour(ps,res,psp);
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

    }
    else if( answer=="plot errors" )
    {
      if( computeErrors || helmholtzFromDirectSolverWasComputed )
      {
        realCompositeGridFunction & error = cgWave.dbase.get<realCompositeGridFunction>("error");
        if( !computeErrors && helmholtzFromDirectSolverWasComputed )
        {
          // plot difference between WaveHoltz and Direct Helmholtz solve
          realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
          error = uHelmholtz - v;

          psp.set(GI_TOP_LABEL,sPrintF("Error (disc. sol) O%d omega=%.5g",orderOfAccuracy,omega));
        }
        else
        {
          psp.set(GI_TOP_LABEL,sPrintF("Error (true soln) O%d omega=%.5g",orderOfAccuracy,omega));

        }
        ps.erase();
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
        PlotIt::contour(ps,error,psp);
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
      }
      else
      {
        printF("plot errors: WARNING: this solution has no errors to compute. Probably not a known solution.\n"); 
        printF("           : Compute the direct Helmholtz solution in order to plot the difference between it and the WaveHoltz solution.\n"); 
      }      
    }
    else if( answer=="plot forcing" )
    {
      CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
      realCompositeGridFunction & f = cgWave.dbase.get<realCompositeGridFunction>("f");      
      ps.erase();
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      psp.set(GI_TOP_LABEL,sPrintF("forcing"));
      PlotIt::contour(ps,f,psp);
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

    }
    else if( answer.matches("erase") )
    {
      replot=false;
      plotChoices=0; 
      ps.erase();
    }
    else if( answer=="run cgWave and plot" )    
    {
      printF("Run cgWave and plot the solution over time...\n");
      int & plotChoices  = cgWave.dbase.get<int>("plotChoices");
      int & plotOptions  = cgWave.dbase.get<int>("plotOptions");

      plotChoices=2; // contours
      plotOptions = CgWave::plotAndWait; // plotNoWait

      cgWave.dbase.get<real>("omega")     = omega;

      int it=0;
      cgWave.advance( it );

    }
    else
    {
      printF("cgwh: Unknown command = [%s]\n",(const char*)answer);
      ps.stopReadingCommandFile();
       
    }

    if( reComputeErrors )
    {
      reComputeErrors=false;

      if( computeErrors )
      {
        realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
        Real t=0; // compute errors with t=0 in cos(omega*t)
        Real maxErr = cgWave.getErrors( v, t );
        printF("CgWaveHoltz: max-err =%8.2e (between WaveHoltz and known solution)\n",maxErr);

        Real solutionNorm = maxNorm( v ); 
        cgWaveHoltz.saveCheckFile( checkFileCounter,maxErr,solutionNorm ); 
        checkFileCounter++; 
      }

      if( helmholtzFromDirectSolverWasComputed )
      {
        realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
        realCompositeGridFunction err(cg);
        err = uHelmholtz - v;
        Real maxDiff = maxNorm( err );
        printF("CgWaveHoltz: max-diff=%8.2e (between WaveHoltz and Direct Helmholtz solution)\n",maxDiff);

        Real solutionNorm = maxNorm ( v );
        cgWaveHoltz.saveCheckFile( checkFileCounter,maxDiff,solutionNorm ); 
        checkFileCounter++; 

      }
    }

    if( answer=="contour" || (replot && plotChoices & 2) )
    {
      CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
      printF("plot contours...\n");

      realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
      ps.erase();
      if( answer=="contour" )
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false); // wait inside contour

      // const bool useUpwind = ad4>0.;
      psp.set(GI_TOP_LABEL,sPrintF("CgWaveHoltz: FD%i%i%s omega=%.5g",orderOfAccuracy,orderOfAccuracyInTime,(upwind ? "u" : ""),omega));

      PlotIt::contour(ps,v,psp);

      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

      plotChoices |= 2; 
    }    

  }
  
  ps.popGUI();  // pop dialog

  if( showFileIsOpen )
    show.close();

  // CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
  cgWave.printStatistics();

  
  Overture::finish();          
  return 0;
}
