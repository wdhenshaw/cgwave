// ---------------------------------------------------------------------------
//
// Test routine for LOCAL COMPATIBILITY BOUNDARY CONDITIONS class
//
//      *** WORK IN PROGRESS***
//
// ---------------------------------------------------------------------------

#include "CompositeGridOperators.h"
#include "PlotStuff.h"
#include "CgWave.h"
#include "LCBC.h"

#include "OGTrigFunction.h"
#include "OGPolyFunction.h"
#include "display.h"

#include "ParallelUtility.h"

#include "gridFunctionNorms.h"
#include "LCBC_annulusMap.h"

#define ForBoundary(side,axis)   for( int axis=0; axis<mg.numberOfDimensions(); axis++ ) \
for( int side=0; side<=1; side++ )
int
getLineFromFile( FILE *file, char s[], int lim);

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

/* Create Global Variable for LCBC function pointers */

#define dimBasedValue(dim,n1,n2) ((dim == 2) ? (n1) : (n2))
#define ind2(i,j,N1) (i + j*(N1))

OGFunction *ogFun;
typedef double (*LCBCFn)(double *);

void getLcbcCoef(MappedGrid mg, double **&coef, RealArray coefArray[], RealArray rxLocal, aString caseName, int p, int numberOfDimensions, int faceEval[]);
void getLcbcDataGrids(MappedGrid mg, double **&fn, double **&gn, RealArray tmpGn[], RealArray tmpFn[], OGFunction &e, RealArray xLocal, double t, int p, int numberOfDimensions, int faceEval[]);

// =========================== MAIN =====================================
int
main(int argc, char *argv[])
{
    Overture::start(argc,argv);
    // Use this to avoid un-necessary communication:
    Optimization_Manager::setForceVSG_Update(Off);
    const int myid=Communication_Manager::My_Process_Number;

    int debug=0; // debug = 1 + 2 + 4 + 8 + ...
    int numResolutions=4;

    int orderOfAccuracyInSpace = 2;
    int orderOfAccuracyInTime  = orderOfAccuracyInSpace;

    CgWave::TwilightZoneEnum twilightZone = CgWave::polynomial;
    int degreeInSpace=orderOfAccuracyInSpace, degreeInTime=orderOfAccuracyInTime;
    Real fx=1.;
    RealArray trigFreq(4);
    trigFreq=fx;

    int plotOption=true;
    bool smartRelease=false;
    bool reportMemory=false;
    bool loadBalance=false;
    int numberOfParallelGhost=2;

    aString caseName = "square"; // or "annulus"

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
            else if( len=line.matches("-caseName=") )
            {
                caseName = line(len,line.length()-1);
                printF("Setting caseName=[%s]\n",(const char*)caseName);
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
            else if( len=line.matches("-degreeInTime=") )
            {
                sScanF(line(len,line.length()-1),"%i",&degreeInTime);
                printF("Setting degreeInTime=%d\n",degreeInTime);
            }
            else if( len=line.matches("-orderInSpace=") )
            {
                sScanF(line(len,line.length()-1),"%i",&orderOfAccuracyInSpace);
                printF("Setting orderOfAccuracyInSpace=%d\n",orderOfAccuracyInSpace);
            }
            else if( len=line.matches("-orderInTime=") )
            {
                sScanF(line(len,line.length()-1),"%i",&orderOfAccuracyInTime);
                printF("Setting orderOfAccuracyInSpace=%d\n",orderOfAccuracyInTime);
            }
            else if( len=line.matches("-fx=") ) // -fx=.12345
            {
                sScanF(line(len,line.length()-1),"%e",&fx);
                printF("Setting trig frequency: fx=%g\n",fx);
                trigFreq=fx;
            }
            else if( len=line.matches("-debug=") )
            {
                sScanF(line(len,line.length()-1),"%i",&debug);
                printF("Setting debug=%i\n",debug);
            }


            else if( commandFileName=="" )
            {
                commandFileName=line;
                printF("testLCBC: reading commands from file [%s]\n",(const char*)commandFileName);
            }

        }
    }
    else
        printF("Usage: `testLCBC [options][file.cmd]' \n"
               "     options:                            \n"
               "          -noplot:        : run without graphics. \n"
               "          -caseName       : defines the grids to use.\n"
               "          -orderInSpace   : [2|4|6...] order of accuracy in space. \n"
               "          -orderInTime    : [2|4|6...] order of accuracy in time. \n"
               "          -tz=[poly|trig] : exact solution type. \n"
               "          -degreeInSpace  : degree for manufactured solutions. \n"
               "          -degreeInTime   : degree for manufactured solutions. \n"
               "          -fx=<f>         : trig frequency. \n"
               "          -debug=<i>      : debug (bit flag). \n"
               "     file.cmd: read this command file. \n");

    GenericGraphicsInterface & ps = *Overture::getGraphicsInterface("testLCBC",false,argc,argv);
    PlotStuffParameters psp;


    // By default start saving the command file called "testLCBC.cmd"
    aString logFile="testLCBC.cmd";
    ps.saveCommandFile(logFile);
    printF("User commands are being saved in the file `%s'\n",(const char *)logFile);

    ps.appendToTheDefaultPrompt("testLCBC>");

    // read from a command file if given
    if( commandFileName!="" )
    {
        printF("read command file =%s\n",(const char*)commandFileName);
        ps.readCommandFile(commandFileName);
    }

    const int maxResolutions=6;
    const int maxDerivDim=20;
    RealArray maxErr(maxResolutions);
    maxErr=0.;

    RealArray l2Err(maxResolutions);
    l2Err=0.;

    RealArray cpu(maxResolutions);
    cpu=0.;

    int GridNum[maxResolutions][3];

    int numberOfDimensions=2; // default
    aString gridNames[maxResolutions];
    if( caseName=="annulus" )
    {
        gridNames[0]= "annulus1.order6.ng4";
        gridNames[1]= "annulus2.order6.ng4";
        gridNames[2]= "annulus4.order6.ng4";
        gridNames[3]= "annulus8.order6.ng4";
        gridNames[4]="annulus16.order6.ng4";
    }
    else if( caseName=="square" )
    {
        gridNames[0]= "square8.order6";
        gridNames[1]= "square16.order6";
        gridNames[2]= "square32.order6";
        gridNames[3]= "square64.order6";
        gridNames[4]="square128.order6";
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
    else if( caseName=="box" )
    {
        numberOfDimensions=3;
        gridNames[0]= "box1.order4.ng3";
        gridNames[1]= "box2.order4.ng3";
        gridNames[2]= "box4.order4.ng3";
        gridNames[3]= "box8.order4.ng3";
        gridNames[4]= "box16.order4.ng3";
    }
    else if( caseName=="nonBox" )
    {
        numberOfDimensions=3;
        gridNames[0]= "nonBox1.order4";
        gridNames[1]= "nonBox2.order4";
        gridNames[2]= "nonBox4.order4";
        gridNames[3]= "nonBox8.order4";
    }
    else if( caseName=="twoSquares" )
    {
        numberOfDimensions=2;
        gridNames[0]= "twoSquarese1.order6";
        gridNames[1]= "twoSquarese2.order6";
        gridNames[2]= "twoSquarese4.order6";
        gridNames[3]= "twoSquarese8.order6";
        gridNames[4]= "twoSquarese16.order6";

    }else if( caseName=="twoBoxes" ){
        numberOfDimensions=3;
        gridNames[0]= "twoBoxese1.order6.hdf";
        gridNames[1]= "twoBoxese2.order6.hdf";
        gridNames[2]= "twoBoxese4.order6.hdf";
        gridNames[3]= "twoBoxese8.order6.hdf";
        gridNames[4]= "twoBoxese16.order6.hdf";
    }
    else
    {
        OV_ABORT("unknown caseName");
    }

    // --- Create the twilightzone exact solution ----
    OGFunction *tz=NULL;
    createTwilightZoneSolution( twilightZone, tz, numberOfDimensions, degreeInSpace, degreeInTime, trigFreq );
    assert( tz!=NULL );
    OGFunction & e = *tz;

    /* NEW: Adding this line to create function pointers for LCBC */
    ogFun = tz;
    LCBCFn coefFnPtr, forcingFnPtr, bdryFnPtr;
    if(caseName == "annulus"){
        coefFnPtr    = annulus_Cn;
        forcingFnPtr = annulus_Fn;
        bdryFnPtr    = annulus_Gn;
    }
    else{
        coefFnPtr    = NULL;
        forcingFnPtr = NULL;
        bdryFnPtr    = NULL;
    }
    /* ----------------------------------------------*/

    // --------- LOOP OVER GRID RESOLUTIONS -----
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

        assert( numberOfDimensions == cg.numberOfDimensions() );

        /* Adding new material here */
        int p = orderOfAccuracyInSpace/2;
        int q = 1; // this is the number of time derivatives in the temporal operator of the pde

        Range all;
        realCompositeGridFunction u(cg,all,all,all);   // holds solution
        u.setName("u",0);
        u=0.;


        realCompositeGridFunction ue(cg,all,all,all);   // holds true solution
        ue.setName("ue",0);
        ue=0.;

        realCompositeGridFunction err(cg,all,all,all); // holds the error
        err.setName("err",0);
        err=0.;

        Index I1,I2,I3;    // index values for the interior
        Index Ib1,Ib2,Ib3; // index values for boundary
        Index Ig1,Ig2,Ig3; // index values for ghost
        Real t=0.;

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            mg.update(MappedGrid::THEmask | MappedGrid::THEcenter | MappedGrid::THEvertex ); // build geometry arrays

            OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
            OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
            OV_GET_SERIAL_ARRAY(real,ue[grid],ueLocal);

            getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points

            const int includeGhost = 1;
            bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
            if( ok )
            {
                int numberOfComponents=1;
                Range C=numberOfComponents;
                int isRectangular=0;
                // -- evaluate the exact solution:
                //    ntd = number of t-derivatives of the exact solution
                //    nxd = number of x-derivatives of the exact solution
                //    nyd = number of y-derivatives of the exact solution
                //    nzd = number of z-derivatives of the exact solution
                int ntd=0, nxd=0, nyd=0, nzd=0;
                e.gd( ueLocal ,xLocal,numberOfDimensions,isRectangular,ntd,nxd,nyd,nzd,I1,I2,I3,C,t); // assign exact solution at time t

                uLocal(I1,I2,I3) = ueLocal(I1,I2,I3);
            }
            
        } // end for grid

        // --- LOCAL COMPATIBILITY BOUNDARY CONDITIONS -----
        Real cpu0 = getCPU();
        
        Lcbc lcbc[cg.numberOfComponentGrids()];

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            // build geometry arrays
            mg.update(MappedGrid::THEmask | MappedGrid::THEcenter | MappedGrid::THEvertex | MappedGrid::THEinverseVertexDerivative );

            OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
            OV_GET_SERIAL_ARRAY(real,ue[grid],ueLocal);
            OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);                   // array of grid points xLocal(i1,i2,i3,0:nd-1)
            OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal); // metrics: rxLocal(i1,i2,i3,0:*) rx, ry, sx, sy
            
            getIndex(cg[grid].dimension(),I1,I2,I3);

            printF("LCBC: grid=%d orderOfAccuracyInSpace=%d, orderOfAccuracyInTime=%d\n",grid,orderOfAccuracyInSpace,orderOfAccuracyInTime);

            ::display(mg.dimension(),"dimension","%5i ");             // index bounds for all grid points including ghost
            ::display(mg.gridIndexRange(),"gridIndexRange","%5i ");   // index bounds for grid points including boundaries


            /* ===== Grid information needed for LCBC ===== */

            printf("\n======= p = %d ======\n",p);

            int userNumGhost = -mg.dimension(0,0);
            int Nx[3]={mg.gridIndexRange(1,0),mg.gridIndexRange(1,1),mg.gridIndexRange(1,2)};
            GridNum[res][0] = Nx[0]; GridNum[res][1] = Nx[1]; GridNum[res][2] = Nx[2];
            printf("userNumGhost = %d and Nx = %d,%d,%d\n", userNumGhost,  Nx[0], Nx[1], Nx[2]);
            
            
            /* ===== BUILD the coefficient grid ===== */

            int maxCoefNum = ((numberOfDimensions+1)*(numberOfDimensions + 2))/2;
            int faceNum = 2*numberOfDimensions;
            RealArray coefArray[(maxCoefNum*faceNum)];
            double **coef = new double*[(maxCoefNum*faceNum)];
            int faceEval[faceNum];

            ForBoundary(side,axis)
            {
                int face = side + 2*axis;
                if(mg.boundaryCondition(side,axis)>0){
                    faceEval[face] = 1;
                }else if(mg.boundaryCondition(side,axis)<0){
                    faceEval[face] = -1;
                }
                else{
                    faceEval[face] = 0;
                }
            }// end of ForBoundary
            
            getLcbcCoef( mg, coef, coefArray, rxLocal, caseName, p, numberOfDimensions, faceEval);
            
            /* ===================================== */
            
            /* Call the LCBC initialize */
            printF("Call the LCBC constructor...\n");
            
            lcbc[grid].initialize(numberOfDimensions,orderOfAccuracyInSpace,orderOfAccuracyInTime,Nx,userNumGhost,faceEval, coef, coefFnPtr, bdryFnPtr, forcingFnPtr);
            printF("... done LCBC constructor\n");
            /* =========== EVALUATE BOUNDARY ================ */

            int NU = (p+1);
            double **gn, **fn;
            gn = new double*[(faceNum*NU)];
            fn = new double*[(faceNum*NU)]; // needs orderInSpace as ghost points
            
            RealArray tmpFn[(faceNum*NU)];
            RealArray tmpGn[(faceNum*NU)];
            
            getLcbcDataGrids( mg, fn, gn, tmpGn, tmpFn, e, xLocal, t, p, numberOfDimensions, faceEval);

            /*========== END OF EVALUATE BOUNDARY ================ */
            /* Update the solution at the ghost points */
            double *un = uLocal.getDataPointer();
            printf("updating %d ghost points normal to each boundary point\n",(orderOfAccuracyInSpace/2));
            lcbc[grid].updateGhost(un, t, 0.1, gn, fn);

            delete [] gn;
            gn = NULL;
            delete [] fn;
            fn = NULL;
            delete [] coef;
            coef = NULL;
            
            printf("Finished grid = %d\n", grid);

        } // end for grid

        cpu(res) = getCPU() - cpu0;
        printF("Time for LCBC = %9.3e (s)\n",cpu(res));

        // ---- compute errors ----
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];

            OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);
            OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
            OV_GET_SERIAL_ARRAY(real,ue[grid],ueLocal);
            OV_GET_SERIAL_ARRAY(real,err[grid],errLocal);

            getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points

            const int includeGhost=1;
            bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
            if( ok )
            {
                errLocal(I1,I2,I3) = ueLocal(I1,I2,I3) - uLocal(I1,I2,I3);

                Real gridErrMax = max(fabs(errLocal));
                maxErr(res) = max( maxErr(res),gridErrMax );

            }
        }
        printF("LCBC: res=%d: max-err = %9.3e\n",res,maxErr(res));


        if( plotOption )
        {

            // --- plot the solution ----
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
            psp.set(GI_TOP_LABEL,"Solution");
            psp.set(GI_NUMBER_OF_GHOST_LINES_TO_PLOT,p); // plot ghost by default
            ps.erase();
            PlotIt::contour( ps,u,psp );

            // --- plot the error ----
            psp.set(GI_TOP_LABEL,"Error");
            ps.erase();
            PlotIt::contour( ps,err,psp );

        }


    } // end resolutions

    printf("\nFor ");
    std::cout << caseName << "\n";
    printf("\n======= LCBC Grid Resolution on Study Order %d ======\n",orderOfAccuracyInSpace);
    for(int res = 0; res<numResolutions; res++){
        if(res!=0){
            printf("Nx=%2d, Ny=%2d, Nz=%2d, MaxErr = %1.5e, rate = %2.2f, cpu = %1.3e\n",GridNum[res][0],GridNum[res][1],GridNum[res][2],maxErr(res),(log(maxErr(res-1)/maxErr(res))/log(2)),cpu(res));
        }
        else{            printf("Nx=%2d, Ny=%2d, Nz=%2d, MaxErr = %1.5e,              cpu = %1.3e\n",GridNum[res][0],GridNum[res][1],GridNum[res][2],maxErr(res),cpu(res));
        }
        
//        if(res!=0){
//            printf("Nx=%2d, Ny=%2d, Nz=%2d, MaxErr = %1.16e, rate = %2.2f\n",GridNum[res][0],GridNum[res][1],GridNum[res][2],maxErr(res),(log(maxErr(res-1)/maxErr(res))/log(2)));
//        }
//        else{            printf("Nx=%2d, Ny=%2d, Nz=%2d, MaxErr = %1.16e\n",GridNum[res][0],GridNum[res][1],GridNum[res][2],maxErr(res));
//        }
    }
    printf("\n");


    Overture::finish();
    return 0;
}

void getLcbcCoef(MappedGrid mg, double **&coef, RealArray coefArray[], RealArray rxLocal, aString caseName, int p, int numberOfDimensions, int faceEval[]){
    
    int maxCoefNum = ((numberOfDimensions+1)*(numberOfDimensions + 2))/2;
    int faceNum = 2*numberOfDimensions;
    Index ib[3];
    Range ibf[3];
    Range ibf_ext[3];
    int includeGhost = p;
    
    ForBoundary(side,axis)
    {
        int face = side + 2*axis;
        if(faceEval[face]>0){
        getBoundaryIndex(mg.gridIndexRange(),side,axis,ib[0],ib[1],ib[2],includeGhost);
        
        for(int ax = 0; ax<3; ax++){
            if(ax<numberOfDimensions){
                ibf_ext[ax] = Range((ib[ax].getBase() - p),(ib[ax].getBound() + p));
                if(ax == axis){
                    ibf[ax] = Range((ib[ax].getBase() - p),(ib[ax].getBound() + p));
                }else{
                    ibf[ax] = Range(ib[ax].getBase(),ib[ax].getBound());
                }
            }else{
                ibf_ext[ax] = Range(0,0);
                ibf[ax] = Range(0,0);
            }
        }
        
        if(caseName == "square" || caseName == "twoSquares"){
            
            for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
                coefArray[ind2(face,coefNum,faceNum)] = RealArray(ibf_ext[0],ibf_ext[1],ibf_ext[2]);
            }
            coefArray[ind2(face,0,faceNum)](ibf[0],ibf[1],ibf[2]) = rxLocal(ibf[0],ibf[1],ibf[2],0)*rxLocal(ibf[0],ibf[1],ibf[2],0) +rxLocal(ibf[0],ibf[1],ibf[2],1)*rxLocal(ibf[0],ibf[1],ibf[2],1);
            coefArray[ind2(face,1,faceNum)](ibf[0],ibf[1],ibf[2]) = rxLocal(ibf[0],ibf[1],ibf[2],2)*rxLocal(ibf[0],ibf[1],ibf[2],2) +rxLocal(ibf[0],ibf[1],ibf[2],3)*rxLocal(ibf[0],ibf[1],ibf[2],3);
            coefArray[ind2(face,2,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0; // c1 ux rxx + ryy
            coefArray[ind2(face,3,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
            coefArray[ind2(face,4,faceNum)](ibf[0],ibf[1],ibf[2]) = rxLocal(ibf[0],ibf[1],ibf[2],0)*rxLocal(ibf[0],ibf[1],ibf[2],2) +rxLocal(ibf[0],ibf[1],ibf[2],1)*rxLocal(ibf[0],ibf[1],ibf[2],3);
            coefArray[ind2(face,5,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
            
            for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
                coef[ind2(face,coefNum,faceNum)] = coefArray[ind2(face,coefNum,faceNum)].getDataPointer();
            }
        }// end of if square
        else if(caseName == "box" || caseName == "twoBoxes"){
            for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
                coefArray[ind2(face,coefNum,faceNum)] = RealArray(ibf_ext[0],ibf_ext[1],ibf_ext[2]);
            }
            coefArray[ind2(face,0,faceNum)](ibf[0],ibf[1],ibf[2]) = rxLocal(ibf[0],ibf[1],ibf[2],0)*rxLocal(ibf[0],ibf[1],ibf[2],0) + rxLocal(ibf[0],ibf[1],ibf[2],1)*rxLocal(ibf[0],ibf[1],ibf[2],1) + rxLocal(ibf[0],ibf[1],ibf[2],2)*rxLocal(ibf[0],ibf[1],ibf[2],2);
            coefArray[ind2(face,1,faceNum)](ibf[0],ibf[1],ibf[2]) = rxLocal(ibf[0],ibf[1],ibf[2],3)*rxLocal(ibf[0],ibf[1],ibf[2],3) + rxLocal(ibf[0],ibf[1],ibf[2],4)*rxLocal(ibf[0],ibf[1],ibf[2],4) + rxLocal(ibf[0],ibf[1],ibf[2],5)*rxLocal(ibf[0],ibf[1],ibf[2],5);
            coefArray[ind2(face,2,faceNum)](ibf[0],ibf[1],ibf[2]) = rxLocal(ibf[0],ibf[1],ibf[2],6)*rxLocal(ibf[0],ibf[1],ibf[2],6) + rxLocal(ibf[0],ibf[1],ibf[2],7)*rxLocal(ibf[0],ibf[1],ibf[2],7) + rxLocal(ibf[0],ibf[1],ibf[2],8)*rxLocal(ibf[0],ibf[1],ibf[2],8);
            coefArray[ind2(face,3,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
            coefArray[ind2(face,4,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
            coefArray[ind2(face,5,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
            coefArray[ind2(face,6,faceNum)](ibf[0],ibf[1],ibf[2]) = rxLocal(ibf[0],ibf[1],ibf[2],0)*rxLocal(ibf[0],ibf[1],ibf[2],3) + rxLocal(ibf[0],ibf[1],ibf[2],1)*rxLocal(ibf[0],ibf[1],ibf[2],4)  + rxLocal(ibf[0],ibf[1],ibf[2],2)*rxLocal(ibf[0],ibf[1],ibf[2],5);
            coefArray[ind2(face,7,faceNum)](ibf[0],ibf[1],ibf[2]) = rxLocal(ibf[0],ibf[1],ibf[2],0)*rxLocal(ibf[0],ibf[1],ibf[2],6) + rxLocal(ibf[0],ibf[1],ibf[2],1)*rxLocal(ibf[0],ibf[1],ibf[2],7)  + rxLocal(ibf[0],ibf[1],ibf[2],2)*rxLocal(ibf[0],ibf[1],ibf[2],8);
            coefArray[ind2(face,8,faceNum)](ibf[0],ibf[1],ibf[2]) = rxLocal(ibf[0],ibf[1],ibf[2],3)*rxLocal(ibf[0],ibf[1],ibf[2],6) + rxLocal(ibf[0],ibf[1],ibf[2],4)*rxLocal(ibf[0],ibf[1],ibf[2],7)  + rxLocal(ibf[0],ibf[1],ibf[2],5)*rxLocal(ibf[0],ibf[1],ibf[2],8);
            coefArray[ind2(face,9,faceNum)](ibf[0],ibf[1],ibf[2]) = 0.0;
            
            for(int coefNum = 0; coefNum<maxCoefNum; coefNum++){
                coef[ind2(face,coefNum,faceNum)] = coefArray[ind2(face,coefNum,faceNum)].getDataPointer();
            }
        }
        else{
            coef = NULL;
        }
        }// end of if faceEval
    }// end of for Boundary
}


void getLcbcDataGrids(MappedGrid mg, double **&fn, double **&gn, RealArray tmpGn[], RealArray tmpFn[], OGFunction &e, RealArray xLocal, double t, int p, int numberOfDimensions, int faceEval[]){
    
    int faceNum = 2*numberOfDimensions;
    Index ib[3];
    Range ibf[3];
    Range ibf_ext[3];
    int includeGhost = p;
    Range C = 1;

    ForBoundary(side,axis)
    {
        int face = ind2(side,axis,2);
        if( faceEval[face]>0 )
        {
            //---- prepare the boundary grid functions for LCBC ----

            getBoundaryIndex(mg.gridIndexRange(),side,axis,ib[0],ib[1],ib[2],includeGhost);

            for(int ax = 0; ax<3; ax++){
                if(ax<numberOfDimensions){
                    ibf_ext[ax] = Range((ib[ax].getBase() - p),(ib[ax].getBound() + p));
                    if(ax == axis){
                        ibf[ax] = Range((ib[ax].getBase() - p),(ib[ax].getBound() + p));
                    }else{
                        ibf[ax] = Range(ib[ax].getBase(),ib[ax].getBound());
                    }
                }else{
                    ibf_ext[ax] = Range(0,0);
                    ibf[ax] = Range(0,0);
                }
            }


            for(int nu = 0; nu<(p+1); nu++){
                RealArray sol(ib[0],ib[1],ib[2]);    sol = 0.;
                RealArray ut(ibf[0],ibf[1],ibf[2]);  ut  = 0.;
                RealArray uxx(ibf[0],ibf[1],ibf[2]); uxx = 0.;
                RealArray uyy(ibf[0],ibf[1],ibf[2]); uyy = 0.;
                RealArray uzz(ibf[0],ibf[1],ibf[2]); uzz = 0.;
                
                e.gd( sol,xLocal,numberOfDimensions,0,nu,0,0,0,ib[0],ib[1],ib[2],C,t);

                tmpGn[ind2(face,nu,faceNum)] = RealArray(ib[0],ib[1],ib[2]);
                tmpGn[ind2(face,nu,faceNum)] = sol(ib[0],ib[1],ib[2]);
                gn[ind2(face,nu,faceNum)] = tmpGn[ind2(face,nu,faceNum)].getDataPointer();

                e.gd(  ut,xLocal,numberOfDimensions,0,(nu+1),0,0,0,ibf[0],ibf[1],ibf[2],C,t);
                e.gd( uxx,xLocal,numberOfDimensions,0,nu,2,0,0,ibf[0],ibf[1],ibf[2],C,t);
                e.gd( uyy,xLocal,numberOfDimensions,0,nu,0,2,0,ibf[0],ibf[1],ibf[2],C,t);
                if(numberOfDimensions == 3){
                    e.gd( uzz,xLocal,numberOfDimensions,0,nu,0,0,2,ibf[0],ibf[1],ibf[2],C,t);
                }
                
                tmpFn[ind2(face,nu,faceNum)] = RealArray(ibf_ext[0],ibf_ext[1],ibf_ext[2]); // created over extended range
                tmpFn[ind2(face,nu,faceNum)](ibf[0],ibf[1],ibf[2]) = ut(ibf[0],ibf[1],ibf[2]) - uxx(ibf[0],ibf[1],ibf[2]) - uyy(ibf[0],ibf[1],ibf[2]) - uzz(ibf[0],ibf[1],ibf[2]);
                fn[ind2(face,nu,faceNum)] = tmpFn[ind2(face,nu,faceNum)].getDataPointer();
            }
        }// end of if statement
    }// end of ForBoundary
}
