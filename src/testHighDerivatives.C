// ---------------------------------------------------------------------------
//
// Test routine to check evaluation of high derivatives
//
//      *** WORK IN PROGRESS***
//
// ---------------------------------------------------------------------------

#include "CompositeGridOperators.h"
#include "PlotStuff.h"
#include "CgWave.h"

// #include "SquareMapping.h"
// #include "AnnulusMapping.h"
// #include "MatrixTransform.h"
// #include "DataPointMapping.h"

#include "OGTrigFunction.h"
#include "OGPolyFunction.h"
#include "display.h"

// #include "incompressible/SmCylinderExactSolution.h"
#include "ParallelUtility.h"
// #include "Oges.h"
// #include "CgSolverUtil.h"

#include "gridFunctionNorms.h"

#define ForBoundary(side,axis)   for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) \
                                 for( int side=0; side<=1; side++ )
int 
getLineFromFile( FILE *file, char s[], int lim);

// void display(realArray & u )
// {
//   printF("u.getlength(0)=%i\n",u.getLength(0));
  
//   ::display(u,"u");
// }

#define ogf EXTERN_C_NAME(ogf)
#define ogderiv EXTERN_C_NAME(ogderiv)
extern "C"
{

   /* Here are functions for TZ flow that can be called from fortran */

  real
  ogf(OGFunction *&ep, const real &x, const real &y,const real &z, const int & c, const real & t )
  {
    return (*ep)(x,y,z,c,t);
  }
  
  
  /* return a general derivative */
  void
  ogderiv(OGFunction *&ep, const int & ntd, const int & nxd, const int & nyd, const int & nzd, 
           const real &x, const real &y, const real &z, const real & t, const int & n, real & ud )
  {
    ud=(*ep).gd(ntd,nxd,nyd,nzd,x,y,z,n,t);
  }


}

#define getWaveDerivatives EXTERN_C_NAME(getwavederivatives)
#define hierDeriv EXTERN_C_NAME(hierderiv)

extern "C"
{
  void getWaveDerivatives( const int&nd, const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                           const int&gridIndexRange, const int&dimRange, const int&isPeriodic, 
                           const real & u, const int&mask, const real&rsxy, const real& xy, const int&boundaryCondition, 
                           const int&ipar, const real &rpar, real & maxErr, const int&ierr );

  void hierDeriv( const int&nd, const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
                  const int&gridIndexRange, const int&dimRange, const int&isPeriodic, 
                  const real & u, const int&mask, const real&rsxy, const real& xy, 
                  real & ud, real & ude,
                  const int&boundaryCondition, 
                  const int&ipar, const real &rpar, real & maxErr, real & l2Err, const int&ierr );
}



// =================================================================================================
//   Create a twilightzone exact solution object
// =================================================================================================
int createTwilightZoneSolution( CgWave::TwilightZoneEnum & twilightZone, OGFunction *& tz,
             int numberOfDimensions, int degreeInSpace, int degreeInTime, RealArray & trigFreq )
{

  if( twilightZone==CgWave::polynomial )
  {
    if( degreeInSpace>6 )
    {
      printF("createTZ:ERROR: degreeInSpace=%d is not supported by OGPolyFunction, Reducing to 6\n",degreeInSpace);
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
  else if( twilightZone==CgWave::trigonometric )
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

    // RealArray & trigFreq = dbase.get<RealArray>("trigFreq");

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
    OV_ABORT("createTZ:ERROR: unknown twilightZone");
  }
  
  

  return 0;

}



int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);
  // Use this to avoid un-necessary communication: 
  Optimization_Manager::setForceVSG_Update(Off);
  const int myid=Communication_Manager::My_Process_Number;

  // This macro will initialize the PETSc solver if OVERTURE_USE_PETSC is defined.
  // INIT_PETSC_SOLVER();

  int maxDeriv=2; // max derivative to compute
  int numResolutions=4;

  int orderOfAccuracyInSpace = 4;
  int numGhost = orderOfAccuracyInSpace/2; 

  CgWave::TwilightZoneEnum twilightZone = CgWave::polynomial;
  int degreeInSpace=2, degreeInTime=2;
  Real fx=1.; 
  RealArray trigFreq(4);
  trigFreq=fx; 

  int plotOption=true;
  bool smartRelease=false;
  bool reportMemory=false;
  bool loadBalance=false;
  int numberOfParallelGhost=2;

  aString caseName = "annulus"; // or square
  // aString caseName = "hollowCylinderTT";

  aString nameOfOGFile= "square32.order6.hdf"; 

  aString commandFileName="";
  if( argc > 1 )
  { // look at arguments for "noplot" or some other name
    int len=0;
    aString line;
    for( int i=1; i<argc; i++ )
    {
      line=argv[i];
      if( line=="-noplot" || line=="noplot" )
        plotOption=false;
      else if( len=line.matches("-g=") )
      {
        nameOfOGFile=line(len,line.length()-1);
        // printf("\n$$$$ node %i : use grid=[%s]\n",myid,(const char*)nameOfOGFile);
      }      
      else if( line=="-nopause" || line=="-abortOnEnd" || line=="-nodirect" ||
               line=="-readCollective" || line=="-writeCollective" ||
               line=="nopause" || line=="abortOnEnd" || line=="nodirect" )
        continue; // these commands are processed by getGraphicsInterface below 

      else if( len=line.matches("-tz=") )
      {
        aString tzType = line(len,line.length()-1);
        if( tzType == "poly" )
        {
          twilightZone = CgWave::polynomial;
        }
        else if( tzType == "trig" )
        {
          twilightZone = CgWave::trigonometric;
        }
        else
        {
          printF("ERROR: unknown tzType=[%s]\n",(const char*)tzType);
        }
      } 
      else if( len=line.matches("-numRes=") )
      {
        sScanF(line(len,line.length()-1),"%i",&numResolutions); 
        printF("Setting numResolutions=%d\n",numResolutions);
      }

      else if( len=line.matches("-degreeInSpace=") )
      {
        sScanF(line(len,line.length()-1),"%i",&degreeInSpace); 
        printF("Setting degreeInSpace=%d\n",degreeInSpace);
      }
      else if( len=line.matches("-order=") )
      {
        sScanF(line(len,line.length()-1),"%i",&orderOfAccuracyInSpace); 
        printF("Setting orderOfAccuracyInSpace=%d\n",orderOfAccuracyInSpace);
      }  
      else if( len=line.matches("-deriv=") )
      {
        sScanF(line(len,line.length()-1),"%i",&maxDeriv); 
        printF("Setting maxDeriv=%d\n",maxDeriv);
      }              
      else if( len=line.matches("-fx=") )
      {
        sScanF(line(len,line.length()-1),"%e",&fx); 
        printF("Setting trig frequency: fx=%g\n",fx);
        trigFreq=fx; 
      }           

      else if( line=="memory" )
      {
        reportMemory=true;
        Diagnostic_Manager::setTrackArrayData(TRUE);
      }
      else if( line=="loadBalance" || line=="-loadBalance" )
      {
        loadBalance=true;
      }
      else if( len=line.matches("-numberOfParallelGhost=") )
      {
        sScanF(line(len,line.length()-1),"%i",&numberOfParallelGhost);
        if( numberOfParallelGhost<0 || numberOfParallelGhost>10 )
        {
          printF("ERROR: numberOfParallelGhost=%i is no valid!\n",numberOfParallelGhost);
          OV_ABORT("error");
        }
        printF("Setting numberOfParallelGhost=%i\n",numberOfParallelGhost);
      }
      else if( len=line.matches("-caseName=") )
      {
        caseName = line(len,line.length()-1);
        printF("Setting caseName=[%s]\n",(const char*)caseName);
      }      
      else if( commandFileName=="" )
      {
        commandFileName=line;    
        printF("testHighDerivatives: reading commands from file [%s]\n",(const char*)commandFileName);
      }
      
    }
  }
  else
    printF("Usage: `testHighDerivatives [options][file.cmd]' \n"
            "     options:                            \n" 
            "          -noplot:   run without graphics \n" 
            "          -caseName: defines the exact solution \n" 
            "          -abortOnEnd: abort if command file ends \n" 
            "          -numberOfParallelGhost=<num> : number of parallel ghost lines \n" 
            "          memory:   run with A++ memory tracking\n" 
            "          release:  run with A++ smart release of memory\n"
            "     file.cmd: read this command file \n");

  GenericGraphicsInterface & ps = *Overture::getGraphicsInterface("testHighDerivatives",false,argc,argv);
  PlotStuffParameters psp;


  // By default start saving the command file called "testHighDerivatives.cmd"
  aString logFile="testHighDerivatives.cmd";
  ps.saveCommandFile(logFile);
  printF("User commands are being saved in the file `%s'\n",(const char *)logFile);

  ps.appendToTheDefaultPrompt("testHighDerivatives>");

  // read from a command file if given
  if( commandFileName!="" )
  {
    printF("read command file =%s\n",(const char*)commandFileName);
    ps.readCommandFile(commandFileName);
  }

  const int maxResolutions=6; 
  const int maxDerivDim=20; 
  RealArray maxErr(Range(1,maxDerivDim),maxResolutions);
  maxErr=0.;
  RealArray maxErr2(Range(1,maxDerivDim),maxResolutions);
  maxErr2=0.;

  RealArray l2Err(Range(1,maxDerivDim),maxResolutions);
  l2Err=0.;
  RealArray l2Err2(Range(1,maxDerivDim),maxResolutions);
  l2Err2=0.;  

  aString gridNames[maxResolutions];
  if( caseName=="annulus" )
  {
     gridNames[0]= "annulus1.order10.ng6";
     gridNames[1]= "annulus2.order10.ng6";
     gridNames[2]= "annulus4.order10.ng6";
     gridNames[3]= "annulus8.order10.ng6";
     gridNames[4]="annulus16.order10.ng6";
  }
  else if( caseName=="square" )
  {
     gridNames[0]= "square16.order8";
     gridNames[1]= "square32.order8";
     gridNames[2]= "square64.order8";
     gridNames[3]="square128.order8";    
  }
  else if( caseName=="nonSquare" )
  {
     gridNames[0]= "nonSquare16.order8";
     gridNames[1]= "nonSquare32.order8";
     gridNames[2]= "nonSquare64.order8";
     gridNames[3]="nonSquare128.order8";    
  }  
  else if( caseName=="rotatedSquare" )
  {
     gridNames[0]= "rotatedSquare16.order8";
     gridNames[1]= "rotatedSquare32.order8";
     gridNames[2]= "rotatedSquare64.order8";
     gridNames[3]= "rotatedSquare128.order8";    
  }  
  else if( caseName=="washer" )
  {  // smoothed polygon "washer"
     gridNames[0]="washere1.order8";
     gridNames[1]="washere2.order8";
     gridNames[2]="washere4.order8";
     gridNames[3]="washere8.order8";    
  }     
  else if( caseName=="wiggley" )
  {  // TFI for a wiggley domain 
     gridNames[0]="wiggley2.order8";
     gridNames[1]="wiggley4.order8";
     gridNames[2]="wiggley8.order8";
     gridNames[3]="wiggley16.order8";    
     gridNames[4]="wiggley32.order8";    
  } 
  else if( caseName=="rhombus" )
  {  // TFI for a wiggley domain 
     gridNames[0]="rhombus2.order8";
     gridNames[1]="rhombus4.order8";
     gridNames[2]="rhombus8.order8";
     gridNames[3]="rhombus16.order8";    
     gridNames[4]="rhombus32.order8";    
  }           
  else
  {
    OV_ABORT("unknown caseName");
  }
  for( int res=0; res<numResolutions; res++ )
  {
    nameOfOGFile = gridNames[res];

    // CompositeGrid cg;
    // const int maxWidthExtrapInterpNeighbours=4;  // This means we support 3rd-order extrap, (1,-3,3,-1)
    // aString nameOfGridFile="";
    // nameOfGridFile = readOrBuildTheGrid(ps, cg, loadBalance, numberOfParallelGhost,maxWidthExtrapInterpNeighbours );

   // create and read in a CompositeGrid
    #ifdef USE_PPP
      // On Parallel machines always add at least this many ghost lines on local arrays
      const int numGhost=2;
      MappedGrid::setMinimumNumberOfDistributedGhostLines(numGhost);
    #endif
    CompositeGrid cg;
    getFromADataBase(cg,nameOfOGFile,loadBalance);  

    Range all;
    realCompositeGridFunction u(cg,all,all,all);
    u.setName("u",0);
    u=0.;


    // Derivatives are returned here 
    //   number of distinct m'th derivatives is = m+1 
    // Eg. m=4
    //   uxxxx
    //   uxxxy
    //   uxxyy
    //   uxyyy
    //   uyyyy
    int maxTotalDerivatives=0;
    for( int m=1; m<=maxDeriv; m++ )
      maxTotalDerivatives += m+1; 

    realCompositeGridFunction ud(cg,all,all,all,maxTotalDerivatives);  // holds derivatives 
    realCompositeGridFunction ude(cg,all,all,all,maxTotalDerivatives); // holds errors in derivatives
    ud=0.;
    ude=0.;
    // aString xName="xxxxxxxxxx";
    // aString yName="yyyyyyyyyy";
    aString buff;
    int d=0; 
    for( int m=1; m<=maxDeriv; m++ )
    {
      // form names:
      //   ux, uy,
      //   uxx, uxy, uyy
      for( int n=0; n<=m; n++ )
      {
        aString dName="u"; 
        for( int j=0; j<m-n; j++ ){ dName += "x"; }
        for( int j=0; j<n;   j++ ){ dName += "y"; }
        // printF("d=%d: m=%d, n=%d, dName=[%s]\n",d,m,n,(const char*)dName); 
        ud.setName(dName,d);
        ude.setName(sPrintF(buff,"%sErr",(const char*)dName),d);
        d++; 
      }
    }
    printF(" d=%d, maxTotalDerivatives=%d\n",d,maxTotalDerivatives);
    // assert( d==maxTotalDerivatives );

    // for( int d=0; d<maxTotalDerivatives; d++ )
    // {
    //   ud.setName(sPrintF("d%d",d),d);
    //   ude.setName(sPrintF("errd%d",d),d);
    // }

    // int d=0;  // counts derivatives 
    // for( int m=1; m<=maxDeriv; m++ )
    // {
    //   dName = xName(0,m-1);
    //   for( int n=1; n<=m+1; n++ )
    //   {
    //     ud.setName(dName,d); d++;
    //     dName = xName(0:m-1-n) + yName(0:n-1);
    //   }

    //   maxTotalDerivatives += m+1; 
    // }

    // --- Assign the solution from the twilightzone ----
    OGFunction *tz=NULL;
    int numberOfDimensions = cg.numberOfDimensions(); 

    
    createTwilightZoneSolution( twilightZone, tz, numberOfDimensions, degreeInSpace, degreeInTime, trigFreq );



    assert( tz!=NULL );
    OGFunction & e = *tz;
    Index I1,I2,I3;
    Real t=0.;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg = cg[grid];
      mg.update(MappedGrid::THEmask |MappedGrid::THEcenter | MappedGrid::THEvertex );

      OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
      OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);

      getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points    

      const int includeGhost=1;
      bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
      if( ok )
      {
        int numberOfComponents=1;
        Range C=numberOfComponents;
        int isRectangular=0;
        e.gd( uLocal ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,C,t);
      }
    }

    // --- plot the solution ----
    if( plotOption )
    {
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      psp.set(GI_TOP_LABEL,"Solution");
      ps.erase();
      PlotIt::contour( ps,u,psp );
    }

    // ps.erase();
    // psp.set(GI_TOP_LABEL,"Solution");
    // psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
    // PlotIt::contour( ps,u,psp );  

    // --- compute derivatives -----


    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {

      MappedGrid & mg = cg[grid];
      mg.update(MappedGrid::THEmask |MappedGrid::THEcenter | MappedGrid::THEvertex );
      OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
      OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
      OV_GET_SERIAL_ARRAY(real,ud[grid],udLocal);
      OV_GET_SERIAL_ARRAY(real,ude[grid],udeLocal);


      const bool isRectangular = mg.isRectangular();
      int gridType = isRectangular? 0 : 1;

      // printF("advance: preComputeUpwindUt=%d\n",preComputeUpwindUt); 
      int uc=0;
      int numberOfComponents=1;
      int debug=0; 
      int twilightZone=1;

      int ipar[]={uc,                            // ipar[ 0]
                  numberOfComponents,
                  grid,                          // ipar[ 2]
                  gridType,                      // ipar[ 3]
                  orderOfAccuracyInSpace,        // ipar[ 4]
                  twilightZone,                  // ipar[ 5]
                  debug,
                  numGhost,
                  maxDeriv
                 }; 

      real dx[3]={1.,1.,1.};
      if( isRectangular )
        mg.getDeltaX(dx);
              
      real rpar[20];
      rpar[ 0]=dx[0];
      rpar[ 1]=dx[1];
      rpar[ 2]=dx[2];
      rpar[ 3]=mg.gridSpacing(0);
      rpar[ 4]=mg.gridSpacing(1);
      rpar[ 5]=mg.gridSpacing(2);
      rpar[ 6]= (real &)tz;  // twilight zone pointer
      rpar[ 7]=REAL_MIN;


      OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
      int *maskptr = maskLocal.getDataPointer(); 

      // getIndex(mg.gridIndexRange(),I1,I2,I3);

      real *uptr  = uLocal.getDataPointer();

      real *rxptr;
      if( isRectangular )
        rxptr=uptr;
      else
      {
        mg.update(MappedGrid::THEinverseVertexDerivative);
        OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal);
        rxptr = rxLocal.getDataPointer();
      }
        
      real *xyptr = uptr;
      xyptr = xLocal.getDataPointer();


      int ierr=0;
      if( false )
      {
        getWaveDerivatives( mg.numberOfDimensions(), 
                uLocal.getBase(0),uLocal.getBound(0),
                uLocal.getBase(1),uLocal.getBound(1),
                uLocal.getBase(2),uLocal.getBound(2),
                mg.indexRange(0,0), mg.dimension(0,0), mg.isPeriodic(0),            
                *uptr,*maskptr,*rxptr,*xyptr,
                mg.boundaryCondition(0,0),ipar[0],rpar[0],maxErr(1,res),ierr );

        for( int d=1; d<=maxDeriv; d++ )
        {
          printF(" %s : grid=%d: order=%d, deriv=%d maxErr=%8.2e",(const char*)gridNames[res],
            grid,orderOfAccuracyInSpace, d,maxErr(d,res));
          if( res>0 )
            printF(" ratio=%4.2f rate=%4.2f",maxErr(d,res-1)/maxErr(d,res),log2(maxErr(d,res-1)/maxErr(d,res)));

          printF("\n");
        }
      }

      ierr=0;
      hierDeriv( mg.numberOfDimensions(), 
              uLocal.getBase(0),uLocal.getBound(0),
              uLocal.getBase(1),uLocal.getBound(1),
              uLocal.getBase(2),uLocal.getBound(2),
              mg.indexRange(0,0), mg.dimension(0,0), mg.isPeriodic(0),            
              *uptr,*maskptr,*rxptr,*xyptr,
              *udLocal.getDataPointer(),*udeLocal.getDataPointer(),
              mg.boundaryCondition(0,0),ipar[0],rpar[0],maxErr2(1,res),l2Err2(1,res),ierr );  

      for( int d=1; d<=maxDeriv; d++ )
      {
        printF("HD: %s : grid=%d: order=%d, deriv=%d relMaxErr=%8.2e",(const char*)gridNames[res],
          grid,orderOfAccuracyInSpace, d,maxErr2(d,res));
        if( res>0 )
          printF(" ratio=%6.2f rate=%4.2f",maxErr2(d,res-1)/maxErr2(d,res),log2(maxErr2(d,res-1)/maxErr2(d,res)));

        printF(",  relL2Err=%8.2e",l2Err2(d,res));
        if( res>0 )
          printF(" ratio=%6.2f rate=%4.2f",l2Err2(d,res-1)/l2Err2(d,res),log2(l2Err2(d,res-1)/l2Err2(d,res)));

        printF("\n");
      }  

    // --- plot the solution ----
    if( plotOption )
    {
      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      psp.set(GI_TOP_LABEL,"derivatives");
      ps.erase();
      PlotIt::contour( ps,ud,psp );

      psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      psp.set(GI_TOP_LABEL,"derivative errors");
      ps.erase();
      PlotIt::contour( ps,ude,psp );      
    }

    }


  } // end resolutions


  // Interpolant & interpolant = *new Interpolant(cg); interpolant.incrementReferenceCount();
  // interpolant.setImplicitInterpolationMethod(Interpolant::iterateToInterpolate);


  // SmCylinderExactSolution cylExact;
  // cylExact.initialize( cg, caseName );

  // Real omega;
  // cylExact.getParameter( "omega",omega );
  // printF(">>>>> omega=%12.4e for the exact solution\n",omega);


  // int numberOfComponents=4;
  // Range all;
  // realCompositeGridFunction u(cg,all,all,all,numberOfComponents);
  // u.setName("u1",0);
  // u.setName("u2",1);
  // u.setName("u3",2);
  // u.setName("p",3);
  // u=0.;

  // Index I1,I2,I3;
  // Index Ib1,Ib2,Ib3;
  // real t=0.;
  // int numberOfTimeDerivatives = 0;
  // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  // {
  //    MappedGrid & mg = cg[grid];
  //    getIndex( mg.dimension(),I1,I2,I3);

  //    cylExact.evalSolution( t, cg, grid, u[grid], I1,I2,I3, numberOfTimeDerivatives );

  // }

  // PlotStuffParameters psp;
  // psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
  // psp.set(GI_TOP_LABEL,"Solution");
  // PlotIt::contour( ps,u,psp );

  // CompositeGridOperators cgop(cg);

  // const int numberOfDimensions = cg.numberOfDimensions();
  // const int u1c=0; 
  // const int u2c=1; 
  // const int u3c=2; 
  // const int pc = u1c + numberOfDimensions;

  // Real mu = 1.; // ** FIX ME 

  // Range Uc(u1c,u1c+numberOfDimensions-1);  // displacement components
  // Range C(u1c,pc);                         // all components 

  // // Put the residual here:
  // realCompositeGridFunction res(cg,all,all,all,numberOfComponents+1);
  // res.setName("u1",0);
  // res.setName("u2",1);
  // res.setName("u3",2);
  // res.setName("p",3);
  // res.setName("div",4);
  // const int nDiv=4; // position of div in res

  // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  // {
  //   MappedGrid & mg = cg[grid];
  //   MappedGridOperators & op = cgop[grid];

  //   OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
  //   OV_GET_SERIAL_ARRAY(real,res[grid],resLocal);

  //   resLocal = 0.;

  //   int extra=-1; // skip boundaries -- *fix* me for traction
  //   getIndex( mg.gridIndexRange(),I1,I2,I3,extra ); 
     
  //   RealArray uLap(I1,I2,I3,C);
  //   op.derivative( MappedGridOperators::laplacianOperator,uLocal, uLap,I1,I2,I3,C );

  //   RealArray px(I1,I2,I3), py(I1,I2,I3), pz(I1,I2,I3);
  //   op.derivative(MappedGridOperators::xDerivative,uLocal,px  ,I1,I2,I3,pc);
  //   op.derivative(MappedGridOperators::yDerivative,uLocal,py  ,I1,I2,I3,pc);

  //   resLocal(I1,I2,I3,u1c) = SQR(omega)*uLocal(I1,I2,I3,u1c) + mu*uLap(I1,I2,I3,u1c) - px(I1,I2,I3); 
  //   resLocal(I1,I2,I3,u2c) = SQR(omega)*uLocal(I1,I2,I3,u2c) + mu*uLap(I1,I2,I3,u2c) - py(I1,I2,I3); 
  //   if( numberOfDimensions==3 )
  //   {
  //    op.derivative(MappedGridOperators::zDerivative,uLocal,pz  ,I1,I2,I3,pc);
  //    resLocal(I1,I2,I3,u3c) = SQR(omega)*uLocal(I1,I2,I3,u3c) + mu*uLap(I1,I2,I3,u3c) - pz(I1,I2,I3); 
  //   }

  //   resLocal(I1,I2,I3,pc) = uLap(I1,I2,I3,pc); 

  //   RealArray u1x(I1,I2,I3), u2y(I1,I2,I3), u3z(I1,I2,I3);
  //   op.derivative( MappedGridOperators::xDerivative,uLocal, u1x,I1,I2,I3,u1c );
  //   op.derivative( MappedGridOperators::yDerivative,uLocal, u2y,I1,I2,I3,u2c );
  //   op.derivative( MappedGridOperators::zDerivative,uLocal, u3z,I1,I2,I3,u3c );

  //   resLocal(I1,I2,I3,nDiv) = u1x + u2y + u3z;

  //   ForBoundary(side,axis) 
  //   {
  //     if( mg.boundaryCondition(side,axis)>0 )
  //     {
  //       getBoundaryIndex(mg.gridIndexRange(),side,axis,Ib1,Ib2,Ib3);
  //       OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);

  //       if( caseName == "hollowCylinderDD" || 
  //           caseName == "solidCylinderD" )
  //       {
  //         // check displacement BCs

  //         Real u1Err = max(fabs(uLocal(Ib1,Ib2,Ib3,u1c)));
  //         Real u2Err = max(fabs(uLocal(Ib1,Ib2,Ib3,u2c)));
  //         Real u3Err = max(fabs(uLocal(Ib1,Ib2,Ib3,u3c)));
  //         Real dispErr = max(fabs(uLocal(Ib1,Ib2,Ib3,Uc)));
  //         printF(" grid=%d, (side,axis)=(%d,%d), error in u=0 BC: u1=%8.2e, u2=%8.2e, u3=%8.2e, max=%8.2e\n",
  //              grid,side,axis,u1Err,u2Err,u3Err,dispErr);
  //       }
  //     }
  //   }
  // }

  // // compute the residual norms
  // printF("============= caseName = %s, grid=%s =============\n",(const char*)caseName,(const char*)nameOfGridFile);
  // printF("============= residual norms: omega=%12.6e =============\n",omega);
  // for( int m=0; m<=numberOfComponents; m++ )
  // {
  //   Real maxRes = maxNorm(res,m);
  //   Real l2Res = l2Norm(res,m);

  //   printF("testHighDerivatives: equation %d:  max-res=%9.3e, l2-res=%9.3e  (%s)\n",m,maxRes,l2Res,(const char*)res.getName(m));

  // }


  // // --- plot the residuals ----
  // ps.erase();
  // psp.set(GI_TOP_LABEL,"Residuals");
  // psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
  // PlotIt::contour( ps,res,psp );



  Overture::finish();          
  return 0;
}
