// ---------------------------------------------------------------------------
//
// ====================== Composite Grid Wave Equation Solver  =====================
//     cgWaveMain -grid=square128.order2
// ---------------------------------------------------------------------------

#include "CgWave.h"

#include "Ogshow.h"  
#include "CompositeGridOperators.h"
#include "PlotStuff.h"
#include "display.h"
// #include "ParallelOverlappingGridInterpolator.h"
#include "ParallelUtility.h"
#include "LoadBalancer.h"
#include "gridFunctionNorms.h"
#include "Oges.h"



int 
getLineFromFile( FILE *file, char s[], int lim);


// bool measureCPU=TRUE;
// real
// CPU()
// // In this version of getCPU we can turn off the timing
// {
//   if( measureCPU )
//     return getCPU(); // return MPI_Wtime();
//   else
//     return 0;
// }

// real
// getMaxTime(real time)
// {
//   real maxTime=0;
// //  MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
//   return maxTime;
// }


// real
// getDt(const real & cfl, 
//       const real & a, 
//       const real & b, 
//       const real & nu, 
//       MappedGrid & mg, 
//       MappedGridOperators & op,
//       const real alpha0 = -2.,
//       const real beta0  = 1. );

// // fourth order dissipation 2D:
// #define FD4_2D(u,i1,i2,i3) \
//       (    -( u(i1-2,i2,i3)+u(i1+2,i2,i3)+u(i1,i2-2,i3)+u(i1,i2+2,i3) )   \
//         +4.*( u(i1-1,i2,i3)+u(i1+1,i2,i3)+u(i1,i2-1,i3)+u(i1,i2+1,i3) ) \
//        -12.*u(i1,i2,i3) )

// // fourth order dissipation 3D:
// #define FD4_3D(u,i1,i2,i3) \
//       (    -( u(i1-2,i2,i3)+u(i1+2,i2,i3)+u(i1,i2-2,i3)+u(i1,i2+2,i3)+u(i1,i2,i3-2)+u(i1,i2,i3+2) )   \
//         +4.*( u(i1-1,i2,i3)+u(i1+1,i2,i3)+u(i1,i2-1,i3)+u(i1,i2+1,i3)+u(i1,i2,i3-1)+u(i1,i2,i3+1) ) \
//        -18.*u(i1,i2,i3) )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  


//   // Pulse parameters:
// static real alpha=50.; // 200.;
// static real pulsePow=10.; // 20
// static real a0=3.;
// static real xPulse=-1.2;
// static real yPulse=0.;



//=============================================================================================================
//     CgWaveMain
//        Solve the wave equation on composite grids
// 
//=============================================================================================================
int 
main(int argc, char *argv[])
{
  int debug=0;

  Overture::start(argc,argv);  // initialize Overture
  const int myid=Communication_Manager::My_Process_Number;
  const int np = Communication_Manager::numberOfProcessors();
  
  printF("Usage: `mpirun -np N cgWave [-noplot] [file.cmd] [-g=<gridName>]'\n");
  
  // Use this to avoid un-necessary communication: 
  Optimization_Manager::setForceVSG_Update(Off);
  Overture::turnOnMemoryChecking(true);

  // This macro will initialize the PETSc solver if OVERTURE_USE_PETSC is defined.
  INIT_PETSC_SOLVER();

  aString nameOfOGFile= "square32.order2.hdf";

  aString commandFileName="";
  aString line;
  int len=0;
  bool plotOption=true;  // by default we plot interactively

  int numParGhost=2; // number of parallel ghost : 2 for order2, 3 for order 4
  

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
      else if( len=line.matches("-numParGhost=") )
      {
        sScanF(line(len,line.length()-1),"%i",&numParGhost);
        printF("Setting numParGhost=%i\n",numParGhost);
      }

      // old way
      else if( len=line.matches("-grid=") )
      {
        nameOfOGFile=line(len,line.length()-1);
        // printf("\n$$$$ node %i : use grid=[%s]\n",myid,(const char*)nameOfOGFile);
      }
      else if( len=line.matches("-cmd=") )
      {
        commandFileName=line(len,line.length()-1);
        // printf("\n$$$$ node %i : read command file %s\n",myid,(const char*)commandFileName);
      }
    }
  }

  
//old  PlotStuff ps(plotOption,"wave equation");

  GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("cgWave",plotOption,argc,argv));

  PlotStuffParameters psp;

  // By default start saving the command file called "cgWave.cmd"
  aString logFile="cgWave.cmd";
  ps.saveCommandFile(logFile);
  printF("User commands are being saved in the file `%s'\n",(const char *)logFile);

  if( commandFileName!="" )
    ps.readCommandFile(commandFileName);

  int orderOfAccuracy=4;  // **** need two ghost lines in parallel ****

  aString nameOfShowFile="cgWave.show";

  // create and read in a CompositeGrid
  #ifdef USE_PPP
    // On Parallel machines always add at least this many ghost lines on local arrays
    // const int numGhost=2;
    MappedGrid::setMinimumNumberOfDistributedGhostLines(numParGhost); 
  #endif
  CompositeGrid cg;
  bool loadBalance=true; // turn on or off the load balancer
  getFromADataBase(cg,nameOfOGFile,loadBalance);

  cg.update(MappedGrid::THEmask);
  const int numberOfDimensions = cg.numberOfDimensions();


  // ** TEST ***
  // cg.update( MappedGrid::THEmask | MappedGrid::THEvertex | MappedGrid::THEcenter | MappedGrid::THEinverseVertexDerivative ); // May 2, 2023

  bool saveShowFile=false;
  

  // ---- Build the CgWave solver and set parameters ----
  CgWave & cgWave = *new CgWave(cg,ps);
  cgWave.setNameOfGridFile(nameOfOGFile);
  cgWave.setup();

  // --- Prompt for changes to options and parameters ---
  cgWave.interactiveUpdate();


  // real & omega     = cgWave.dbase.get<real>("omega");
  int & numPeriods = cgWave.dbase.get<int>("numPeriods");
  real & tol       = cgWave.dbase.get<real>("tol");


  // Build a dialog menu for changing parameters
  GUIState gui;
  DialogData & dialog=gui;

  dialog.setWindowTitle("CgWave");
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
                        "solve",
                        "check known solution",
                        "change parameters",
                        "contour",
                        "grid",
                        "erase",
                        "exit",
                        ""};
  int numRows=5;
  dialog.setPushButtons( pbLabels, pbLabels, numRows ); 

  aString tbCommands[] = {"save show file",
                           ""};
  int tbState[10];
  tbState[0] = saveShowFile==true;
  int numColumns=1;
  dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 

  // ----- Text strings ------
  const int numberOfTextStrings=20;
  aString textCommands[numberOfTextStrings];
  aString textLabels[numberOfTextStrings];
  aString textStrings[numberOfTextStrings];

  int nt=0;
  textCommands[nt] = "number of periods";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numPeriods);  nt++; 
  textCommands[nt] = "tol";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",tol);  nt++; 

  // null strings terminal list
  textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
  dialog.setTextBoxes(textCommands, textLabels, textStrings);

  

  ps.pushGUI(gui);



  psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

  aString answer;
  char buff[200];
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
      cgWave.interactiveUpdate();
    }
    
    else if( answer=="solve" )
    {
      int it=0;
      int & plotOptions = cgWave.dbase.get<int>("plotOptions");
      plotOptions= CgWave::plotAndWait;
       
      cgWave.advance(it);
    }

    else if( answer=="check known solution" )
    {
      printF("cgWaveMain: check known solution\n");
      // OV_ABORT("finish me");
      Range all;
      realCompositeGridFunction u(cg), res(cg), w(cg,all,all,all,2);
      w.setName("u",0);
      w.setName("res",1);

      res=0.;

      CompositeGridOperators cgop(cg);
      const int orderOfAccuracy = cgWave.dbase.get<int>("orderOfAccuracy");
      cgop.setOrderOfAccuracy( orderOfAccuracy ); 

      Index I1,I2,I3;
      int numberOfComponents=1;
      Range C(0,0); // components
      Real t=0.0;   // check this time
      Real omega= 4.4934094579090642; // FIX ME 
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        MappedGrid & mg = cg[grid];
        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter | MappedGrid::THEmask | MappedGrid::THEinverseVertexDerivative );
        MappedGridOperators & op = cgop[grid];

 
        getIndex(mg.dimension(),I1,I2,I3);
        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
        OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
        OV_GET_SERIAL_ARRAY(real,w[grid],wLocal);
        OV_GET_SERIAL_ARRAY(real,res[grid],resLocal);

        int numberOfTimeDerivatives=0;
        cgWave.getUserDefinedKnownSolution(t, grid, u[grid],I1,I2,I3,numberOfTimeDerivatives);

        int extra=-1; // skip boundaries -- *fix* me for Neumann BCs
        getIndex( mg.gridIndexRange(),I1,I2,I3,extra ); 
         
        RealArray uLap(I1,I2,I3,C);
        op.derivative( MappedGridOperators::laplacianOperator,uLocal, uLap,I1,I2,I3,C );

        where( maskLocal(I1,I2,I3)>0 )
          resLocal(I1,I2,I3,C) = SQR(omega)*uLocal(I1,I2,I3,C) + uLap(I1,I2,I3,C); 


        getIndex(mg.dimension(),I1,I2,I3);
        wLocal(I1,I2,I3,0) = uLocal(I1,I2,I3);
        wLocal(I1,I2,I3,1) = resLocal(I1,I2,I3);
      }

      printF("============= residual norms: omega=%12.6e =============\n",omega);
      for( int m=0; m<numberOfComponents; m++ )
      {
        Real maxRes = maxNorm(res,m);
        Real l2Res = l2Norm(res,m);

        printF("testExact: equation %d:  max-res=%9.3e, l2-res=%9.3e\n",m,maxRes,l2Res);

      }      


      ps.erase();
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      psp.set(GI_TOP_LABEL,sPrintF("Known solution at t=%9.3e",t));
      PlotIt::contour(ps,w,psp);
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);              


    }

    else if( answer=="contour" )
    {
      if( false )
      {
         // int & current = dbase.get<int>("current"); // hold the current solution index
         // cgWave.plot( current, t, dt );
      }
      else
      {
        // realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
        realCompositeGridFunction & v = cgWave.getCurrentSolution();
        ps.erase();
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
        psp.set(GI_TOP_LABEL,sPrintF("Wave Equation Solution"));
        PlotIt::contour(ps,v,psp);
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
      } 
    }
    else if( answer=="grid" )
    {
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      PlotIt::plot(ps,cg,psp);                          // plot the grid
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

    }
     else if( answer.matches("erase") )
    {
      ps.erase();
    }
    else
    {
      printF("Unknown command = [%s]\n",(const char*)answer);
      ps.stopReadingCommandFile();
       
    }
    
  }
  
  ps.popGUI();  // pop dialog

  delete & cgWave; // delete here so we shut-down PETSc properly
  
  Overture::finish();          
  return 0;
}
