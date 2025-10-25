// This file automatically generated from cgwh.bC with bpp.
// ---------------------------------------------------------------------------
//
//     cgWaveHoltz -grid=square128.order2
//     cgWaveHoltz -grid=cice8.order2    
//
//  SOLVE THE HELMHOLTZ EQUATION USING THE WAVE HOLTZ ALGORITHM
//
//
// ---------------------------------------------------------------------------

#include "CgWaveHoltz.h"
#include "CgWave.h"

// #include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"
#include "PlotStuff.h"
#include "display.h"
// #include "ParallelOverlappingGridInterpolator.h"
#include "ParallelUtility.h"
#include "ParallelGridUtility.h"
#include "LoadBalancer.h"
// #include "gridFunctionNorms.h"
#include "Integrate.h"



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


// ===============================================================
// Macro: Compute the inner product of two grid functions          
// ===============================================================

// ===============================================================
// Macro: Add to grid functions
//     w = fact1*v + fact2*uev   
// ===============================================================

//=============================================================================================================
///   cgwh:
///     (1) Compute solutions to the Helmholtz equation by solving the time dependent wave equation
///     (2) Compute eigenvalues and eigenvectors
/// 
///  Ref. Daniel Appelo -- WaveHoltz
//=============================================================================================================
int 
main(int argc, char *argv[])
{
    int debug=0;

    Overture::start(argc,argv);  // initialize Overture
    const int myid=Communication_Manager::My_Process_Number;
    const int np = Communication_Manager::numberOfProcessors();
    
    printF("Usage: `mpirun -np N cgwh [-noplot] [file.cmd] [-g=<gridName>]'\n");
    
  // Use this to avoid un-necessary communication: 
    Optimization_Manager::setForceVSG_Update(Off);

  // This macro will initialize the PETSc solver if OVERTURE_USE_PETSC is defined.
    #ifdef CGWAVE_USE_PETSC
    INIT_PETSC_SOLVER();
    #endif

    aString nameOfOGFile= "square32.order2.hdf"; 

    aString commandFileName="";
    aString line;
    int len=0;
    bool plotOption=true;  // by default we plot interactively
    int upwind =0;      // set to 1 for upwinding 
    int numParGhost=-1; // -1 = guess this value from the grid name and upwind 

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
            else if( len=line.matches("-numParGhost=") )
            {
                sScanF(line(len,line.length()-1),"%i",&numParGhost);
                printF("Setting numParGhost=%i\n",numParGhost);
            }      
            else if( len=line.matches("-upwind=") )
            {
                sScanF(line(len,line.length()-1),"%i",&upwind);
                printF("Setting upwind=%i\n",upwind);
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

  // aString nameOfShowFile = "cgWaveHoltz.show";


  // On Parallel machines we need to set the number of parallel ghost BEFORE reading the grid 
    if( numParGhost<0 )
        numParGhost=ParallelGridUtility::setNumberOfParallelGhost( nameOfOGFile, upwind );
    else
        ParallelGridUtility::setNumberOfParallelGhost( numParGhost );

  // #ifdef USE_PPP
  //   // On Parallel machines always add at least this many ghost lines on local arrays
  //   const int numGhost=2;
  //   MappedGrid::setMinimumNumberOfDistributedGhostLines(numGhost);
  // #endif

  // --- create and read in a CompositeGrid ---
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
  // bool showFileIsOpen=false;
    
  // real t=0;
  // real c=1., cSquared=c*c;
        
  // real ad4=1.;  // coeff of the artificial dissipation.

  // real cfl=.25, nu=0.;
  // InitialConditionOptionEnum option=smoothPulse;


  // ---- Build the WaveHoltz solver and set parameters ----
    CgWaveHoltz & cgWaveHoltz = *new CgWaveHoltz(cg,ps);
    cgWaveHoltz.setNameOfGridFile(nameOfOGFile);


    cgWaveHoltz.interactiveUpdate();
    cgWaveHoltz.setup();

    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
    cgWave.setNameOfGridFile(nameOfOGFile);
    cgWave.setup();
    cgWave.interactiveUpdate();

    const int & computeEigenmodes = cgWave.dbase.get<int>("computeEigenmodes");

    if( computeEigenmodes )
    {
  // --- Solve for eigenmodes WaveHoltz ---
        cgWaveHoltz.solveEigen(argc,argv);
    }
    else
    {
    // --- Solve the Helmholtz problem using WaveHoltz ---
        cgWaveHoltz.solveHelmholtz(argc,argv);
    }


  // CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
    const int & numberOfStepsTaken = cgWave.dbase.get<int>("numberOfStepsTaken");  
    if( numberOfStepsTaken>0 )
        cgWave.printStatistics();

    delete & cgWaveHoltz; // delete here so we shut-down PETSc properly
    
    Overture::finish();          
    return 0;
}
