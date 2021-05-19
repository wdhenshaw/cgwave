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



// // Assign the initial conditions
// void
// getInitialConditions( InitialConditionOptionEnum option, realCompositeGridFunction *u, real t, real dt,
//                       real c, const bool plotOption )
// {
//   const int myid=Communication_Manager::My_Process_Number;
//   printF("get initial conditions\n");
  
//   CompositeGrid & cg = *( u[0].getCompositeGrid() );

//   // Pulse parameters:
  
//   Index I1,I2,I3;

// // define U0(x,y,t) exp( - alpha*( SQR((x)-(xPulse-c*dt)) + SQR((y)-yPulse) ) )
// // #define U0(x,y,t) exp( - alpha*( SQR((x)-(xPulse+c*(t))) ) )
// // define U0(x,y,t) exp( - alpha*( pow( a0*( (x)-(xPulse+c*(t)) ),20.) ) )
// #define U0(x,y,t) exp( - alpha*( pow( a0*( (x)-(xPulse+c*(t)) ),pulsePow) ) )

//   int grid;
//   for( grid=0; grid<cg.numberOfComponentGrids(); grid++ )
//   {
//     // initial condition is a pulse, we make an approximate guess for u(-dt) 
//     // u[1] = u(x,-dt) 
//     // u[0] = u(x,t)

//     // get the local serial arrays
//     OV_GET_SERIAL_ARRAY(real,u[0][grid],u0Local);
//     OV_GET_SERIAL_ARRAY(real,u[1][grid],u1Local);

//     getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points.
//     const bool isRectangular=cg[grid].isRectangular();

//     if( option==smoothPulse )
//     {

//       // restrict the bounds (I1,I2,i3) to the local array bounds (including parallel ghost pts):
//       const int includeGhost=1;
//       bool ok=ParallelUtility::getLocalArrayBounds(u[0][grid],u0Local,I1,I2,I3,includeGhost);

//       if( isRectangular )
//       {
//         // for a rectangular grid we avoid building the array of verticies.
//         // we assign the initial conditions with C-style loops

//         if( !ok ) continue;  // nothing to do on this processor

//         real dx[3]={0.,0.,0.}, xab[2][3]={0.,0.,0.,0.,0.,0.};
//         if( cg[grid].isRectangular() )
//           cg[grid].getRectangularGridParameters( dx, xab );

//         const real xa=xab[0][0], dx0=dx[0];
//         const real ya=xab[0][1], dy0=dx[1];
//         const real za=xab[0][2], dz0=dx[2];

//         const int i0a=cg[grid].gridIndexRange(0,0);
//         const int i1a=cg[grid].gridIndexRange(0,1);
//         const int i2a=cg[grid].gridIndexRange(0,2);

//         #define VERTEX0(i0,i1,i2) xa+dx0*(i0-i0a)
//         #define VERTEX1(i0,i1,i2) ya+dy0*(i1-i1a)
//         #define VERTEX2(i0,i1,i2) za+dz0*(i2-i2a)

//         // Here we grab a pointer to the data of the array so we can index it as a C-array
//         real *upm= u1Local.Array_Descriptor.Array_View_Pointer3;
//         real *up = u0Local.Array_Descriptor.Array_View_Pointer3;
//         const int uDim0=u0Local.getRawDataSize(0);
//         const int uDim1=u0Local.getRawDataSize(1);
//         const int d1=uDim0, d2=d1*uDim1; 
//         #define U(i0,i1,i2) up[(i0)+(i1)*d1+(i2)*d2]
//         #define UM(i0,i1,i2) upm[(i0)+(i1)*d1+(i2)*d2]

//         int i1,i2,i3;
//         FOR_3(i1,i2,i3,I1,I2,I3) // loop over all points
//         {
//           UM(i1,i2,i3)=U0(VERTEX0(i1,i2,i3),VERTEX1(i1,i2,i3),-dt);
//           U(i1,i2,i3) =U0(VERTEX0(i1,i2,i3),VERTEX1(i1,i2,i3),0.);
//         }
        
//         #undef VERTEX0
//         #undef VERTEX1
//         #undef VERTEX2
//         #undef U
//         #undef UM
//       }
//       else
//       {
//         cg[grid].update(MappedGrid::THEvertex);  // build the array of vertices
//         realArray & vertex = cg[grid].vertex();

//         if( !ok ) continue;  // nothing to do on this processor

//         // const realSerialArray & xLocal = vertex.getLocalArrayWithGhostBoundaries();
//         // display(vertex,"vertex",NULL,"%4.1f ");
//         // display(xLocal,"xLocal",NULL,"%4.1f ");

// //      u[1][grid]=U0(vertex(I1,I2,I3,0),vertex(I1,I2,I3,1),-dt);
// //      u[0][grid]=U0(vertex(I1,I2,I3,0),vertex(I1,I2,I3,1),0.);

//         real *upm= u1Local.Array_Descriptor.Array_View_Pointer3;
//         real *up = u0Local.Array_Descriptor.Array_View_Pointer3;
//         const int uDim0=u0Local.getRawDataSize(0);
//         const int uDim1=u0Local.getRawDataSize(1);
//         const int d1=uDim0, d2=d1*uDim1; 
//         #define U(i0,i1,i2) up[(i0)+(i1)*d1+(i2)*d2]
//         #define UM(i0,i1,i2) upm[(i0)+(i1)*d1+(i2)*d2]

//         OV_GET_SERIAL_ARRAY(real,vertex,vertexLocal);
//         const real *vertexp = vertexLocal.Array_Descriptor.Array_View_Pointer3;
//         const int vertexDim0=vertexLocal.getRawDataSize(0);
//         const int vertexDim1=vertexLocal.getRawDataSize(1);
//         const int vertexDim2=vertexLocal.getRawDataSize(2);
//         #define VERTEX(i0,i1,i2,i3) vertexp[i0+vertexDim0*(i1+vertexDim1*(i2+vertexDim2*(i3)))]

//         int i1,i2,i3;
//         FOR_3(i1,i2,i3,I1,I2,I3) // loop over all points
//         {
//           UM(i1,i2,i3)=U0(VERTEX(i1,i2,i3,0),VERTEX(i1,i2,i3,1),-dt);
//           U(i1,i2,i3) =U0(VERTEX(i1,i2,i3,0),VERTEX(i1,i2,i3,1),0.);
//         }

//       }
//       #undef U
//       #undef UM
//       #undef VERTEX

//     }
//     else
//     {
//       // discontinuous pulse -- this doesn't work as planned since c*dt < dx and then
//       // u[1] is just the same as u[0] 
//       cg[grid].update(MappedGrid::THEvertex);  // build the array of vertices
//       const realArray & vertex = cg[grid].vertex();
//       u[1][grid]=0.;
//       where( fabs(vertex(I1,I2,I3,0)-(xPulse-c*dt))<.2 )
//       {
//         u[1][grid]=1.;
//       }
//       u[0][grid]=0.;
//       where( fabs(vertex(I1,I2,I3,0)-xPulse)<.2 )
//       {
//         u[0][grid]=1.;
//       }
//     }
//     if( !plotOption ) 
//       cg[grid].destroy(MappedGrid::THEvertex);  // vertices are no nolonger needed.

//   }

//   printF("done initial conditions\n");

// }


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

  
//old  PlotStuff ps(plotOption,"wave equation");

  GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("cgwh",plotOption,argc,argv));

  PlotStuffParameters psp;

  // By default start saving the command file called "cgWaveHoltz.cmd"
  aString logFile="cgWaveHoltz.cmd";
  ps.saveCommandFile(logFile);
  printF("User commands are being saved in the file `%s'\n",(const char *)logFile);

  if( commandFileName!="" )
    ps.readCommandFile(commandFileName);

  // int orderOfAccuracy=4;  // **** need two ghost lines in parallel ****

  aString nameOfShowFile="cgWaveHoltz.show";

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
  bool saveShowFile=false;
  
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

  const int & orderOfAccuracy = cgWave.dbase.get<int>("orderOfAccuracy");
  cgWaveHoltz.dbase.get<int>("orderOfAccuracy")=orderOfAccuracy; // set value in CgWaveHoltz

  real & omega      = cgWaveHoltz.dbase.get<real>("omega");
  int & numPeriods  = cgWaveHoltz.dbase.get<int>("numPeriods");
  real & tol        = cgWaveHoltz.dbase.get<real>("tol");
  int & adjustOmega = cgWaveHoltz.dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

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
                        "change parameters",
                        "contour",
                        "grid",
                        "plot residual",
                        "compute errors",
                        "plot errors",
                        "plot forcing",
                        "erase",
                        "exit",
                        ""};
  int numRows=6;
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
  // textCommands[nt] = "omega";  textLabels[nt]=textCommands[nt];
  // sPrintF(textStrings[nt], "%g",omega);  nt++; 
  textCommands[nt] = "number of periods";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%i",numPeriods);  nt++; 
  textCommands[nt] = "tol";  textLabels[nt]=textCommands[nt];
  sPrintF(textStrings[nt], "%g",tol);  nt++; 

  // textCommands[nt] = "tFinal";  textLabels[nt]=textCommands[nt];
  // sPrintF(textStrings[nt], "%g",tFinal);  nt++; 
  // textCommands[nt] = "tPlot";  textLabels[nt]=textCommands[nt];
  // sPrintF(textStrings[nt], "%g",tPlot);  nt++; 
  // textCommands[nt] = "tShow";  textLabels[nt]=textCommands[nt];
  // sPrintF(textStrings[nt], "%g",tShow);  nt++; 
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
      cgWaveHoltz.interactiveUpdate();
    }
    
    // else if( dialog.getTextValue(answer,"omega","%e",omega) )
    // {
    //   printF("Setting omega=%g\n",omega);
    // }
    else if( dialog.getTextValue(answer,"number of periods","%i",numPeriods) )
    {
      printF("Setting numPeriods=%i\n",numPeriods);
    }
    else if( dialog.getTextValue(answer,"tol","%e",tol) )
    {
      printF("Setting tol=%g (tolerence for Krylov solvers)\n",tol);
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

    }
    else if( answer == "solve Helmholtz directly" )
    {
       // this may require cgWave to have been called to create f ?

      realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
      realCompositeGridFunction & f = cgWave.dbase.get<realCompositeGridFunction>("f");

      const Real cpu0=getCPU();

      cgWaveHoltz.solveHelmholtz( v, f  );

      const Real cpu = getCPU()-cpu0;

      Real maxRes = cgWaveHoltz.residual();
      printF("CgWaveHoltz: omega=%9.3e, max-res=%9.3e, cpu=%9.2e (Direct solution of Helmholtz).\n",
             omega,maxRes,cpu);

    }
    else if( answer=="compute errors" )
    {
      const int & computeErrors = cgWave.dbase.get<int>("computeErrors");
      if( computeErrors )
      {
        realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
        Real t=0; // compute errors with t=0 in cos(omega*t)
        Real maxErr = cgWave.getErrors( v, t );
        printF("CgWaveHoltz: max-err=%9.3e\n",maxErr);
      }
      else
      {
        printF("compute errors: WARNING: this solution has no errors to compute. Probably not a known solution.\n"); 
      }

    }
    else if( answer=="contour" )
    {
      CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");

      realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
      ps.erase();
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      psp.set(GI_TOP_LABEL,sPrintF("Helmholtz Solution O%d omega=%.5g",orderOfAccuracy,omega));
      PlotIt::contour(ps,v,psp);
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

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
      const int & computeErrors = cgWave.dbase.get<int>("computeErrors");
      if( computeErrors )
      {
        realCompositeGridFunction & error = cgWave.dbase.get<realCompositeGridFunction>("error");
        ps.erase();
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
        psp.set(GI_TOP_LABEL,sPrintF("error O%d omega=%.5g",orderOfAccuracy,omega));
        PlotIt::contour(ps,error,psp);
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
      }
      else
      {
        printF("plot errors: WARNING: this solution has no errors to compute. Probably not a known solution.\n"); 
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
      ps.erase();
    }
    else
    {
      printF("cgwh: Unknown command = [%s]\n",(const char*)answer);
      ps.stopReadingCommandFile();
       
    }
    
  }
  
  ps.popGUI();  // pop dialog

  // CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
  cgWave.printStatistics();

  
  Overture::finish();          
  return 0;
}
