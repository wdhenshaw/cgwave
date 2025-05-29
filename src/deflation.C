// This file automatically generated from deflation.bC with bpp.
// ==============================================
//    Deflation Routines for WaveHoltz 
// ==============================================

#include "CgWave.h"
//#include "CompositeGridOperators.h";    
#include "PlotStuff.h"
#include "display.h"
#include "ParallelUtility.h"
//#include "LoadBalancer.h"
#include "gridFunctionNorms.h"
//#include "OGPolyFunction.h"
//#include "ProcessorInfo.h"
#include "ShowFileReader.h"
// #include "InterpolatePointsOnAGrid.h"

#include "Integrate.h"
#ifndef GRID_STATISTICS_H
#define GRID_STATISTICS_H

#include "CompositeGrid.h"

// forward declaration:
class InterfaceInfo;

/// =========================================================================================================
/// This class can be use to compute and output grid statstics such as cell volumes, areas, arclengths as well as 
/// numbers of grid points, boundary conditions etc.
/// =========================================================================================================
class GridStatistics
{
public:

enum MaskOptionEnum
{
    ignoreMask,
    positiveMask,
    nonZeroMask
};


static int
checkForNegativeVolumes(GridCollection & gc, int numberOfGhost=0, FILE *file=stdout, 
                                                bool checkActivePoints=true, bool printNegativeVolumes=true );

static int
checkForNegativeVolumes(MappedGrid & mg, int numberOfGhost=0, FILE *file=stdout, int grid=0, 
                                                bool checkActivePoints=true, bool printNegativeVolumes=true );

// check for tall cells at interfaces
static real
checkForTallCells( CompositeGrid & cg, std::vector<InterfaceInfo> & interfaceInfo, real tallCellRatioBound=-1 );

// Return x,y,z min-max bounds of a component grid
static void 
getGridCoordinateBounds(MappedGrid & mg, Real xMin[3], Real xMax[3], MaskOptionEnum maskOption=positiveMask );

// Return x,y,z min-maxbounds of a grid collection
static void 
getGridCoordinateBounds(GridCollection & gc, Real xMin[3], Real xMax[3], MaskOptionEnum maskOption=positiveMask );

static void 
getGridSpacing(MappedGrid & mg, real dsMin[3], real dsAve[3], real dsMax[3], MaskOptionEnum maskOption=ignoreMask );

static void 
getGridSpacingAndNumberOfPoints(MappedGrid & mg, real dsMin[3], real dsAve[3], real dsMax[3],
                                                                int & numberOfPoints, MaskOptionEnum maskOption=ignoreMask );

static void
getNegativeVolumes( MappedGrid & mg, int & numberOfNegativeVolumes, IntegerArray & negativeVolumeList,
                                        int grid=0, int numberOfGhost=0, bool checkActivePoints=true, bool printNegativeVolumes=true );

static void 
getNumberOfPoints(GridCollection & gc, int & totalNumberOfGridPoints, MaskOptionEnum maskOption=ignoreMask );

static void 
getNumberOfPoints(MappedGrid & mg, int & numberOfPoints, MaskOptionEnum maskOption=ignoreMask );

static void 
printGridStatistics(CompositeGrid & cg, FILE *file=stdout );

static void
printGridStatistics(GridCollection & gc, FILE *file=stdout );

static void 
printGridStatistics(MappedGrid & mg, FILE *file=stdout, int grid=0, int *ipar=NULL, real *rpar=NULL,
                                        int domainNumber = -1, const aString & domainName  = nullString  );

private:

// private routine that actually computes the negative volumes
static int
computeNegativeVolumes(MappedGrid & mg, int numberOfGhost, FILE *file, int grid, 
                   		       bool checkActivePoints, bool printNegativeVolumes,
                   		       int & negativeVolumeCount, IntegerArray *negativeVolumeList = NULL );

// The constructor should never be called
GridStatistics();
~GridStatistics();

};


#endif
#include "InterpolatePointsOnAGrid.h"


#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  
#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )


// ======================================================================
// Macro: Normalize an eigenvector 
// ======================================================================

// ======================================================================
// Macro: Compute the inner product of ui and uj
// ======================================================================



// ============================================================================================
/// \brief Plot the WaveHoltz filter function and any known eigenvalues 
// ============================================================================================
int CgWave::plotFilter( RealArray & eigenValues )
{

  // CgWaveHoltz & cgWaveHoltz = *this;
    CgWave & cgWave = *this;


    bool plotOption=true;  // by default we plot interactively
  // GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("EigenWave",plotOption,argc,argv));
    GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("EigenWave",plotOption));
    PlotStuffParameters psp;  

    const int & computeEigenmodes     = cgWave.dbase.get<int>("computeEigenmodes");
    RealArray & frequencyArray        = cgWave.dbase.get<RealArray>("frequencyArray");
    RealArray & frequencyArrayAdjusted= cgWave.dbase.get<RealArray>("frequencyArrayAdjusted");
    RealArray & periodArray           = cgWave.dbase.get<RealArray>("periodArray"); 
    RealArray & periodArrayAdjusted   = cgWave.dbase.get<RealArray>("periodArrayAdjusted"); 
    IntegerArray & numPeriodsArray    = cgWave.dbase.get<IntegerArray>("numPeriodsArray");
  //RealArray & frequencyArraySave    = cgWave.dbase.get<RealArray>("frequencyArraySave");
  //RealArray & periodArraySave       = cgWave.dbase.get<RealArray>("periodArraySave");   


  // const RealArray & frequencyArray         = cgWave.dbase.get<RealArray>("frequencyArray");
  // const RealArray & frequencyArrayAdjusted = cgWave.dbase.get<RealArray>("frequencyArrayAdjusted");
    const Real & dt                          = cgWave.dbase.get<real>("dtUsed");
  // RealArray & eigenValues            = dbase.get<RealArray>("eigenValues");  // best estimate of eigenvalues
  // const int & numEigenVectors        = dbase.get<int>("numEigenVectors");



    const int numEigenVectors = eigenValues.getLength(0); //  THIS IS NUM TO DEFLATE  ********* FIX ME *******

  // --- We need to be using the adjusted values of frequency etc. 
  // --- Sometimes this routine is called when frequencyArray holds the adjusted and sometimes not
  // --- Make a local copy so we can restore:
    RealArray frequencyArraySave         = frequencyArray;
    RealArray frequencyArrayAdjustedSave = frequencyArrayAdjusted; 
    RealArray periodArraySave            = periodArray;
    const Real dtSave = cgWave.dbase.get<real>("dt"); 
    printF("\n ^^^^^^ plotFilter: numEigenVectors=%d, dt=%9.3e, dtUsed=%9.3e, BEFORE: frequencyArray(0)=%10.3e frequencyArrayAdjusted(0)=%10.3e\n",
        numEigenVectors, dtSave, dt, frequencyArray(0),frequencyArrayAdjusted(0) );

    if( !computeEigenmodes )
    {
    // -- We need to use the adjusted frequencies in getWaveHoltzIterationEigenvalue ---
        frequencyArray = frequencyArrayAdjusted;
        periodArray    = periodArrayAdjusted; 

        cgWave.dbase.get<real>("dt") = cgWave.dbase.get<real>("dtUsed");
    }
  // printF("\n ^^^^^^ plotFilter: numEigenVectors=%d, dt=%9.3e, frequencyArray(0)=%10.3e frequencyArrayAdjusted(0)=%10.3e\n",
  //   numEigenVectors, dt, frequencyArray(0),frequencyArrayAdjusted(0) );

    
    Real maxEig = max(eigenValues);

  // --- plot WaveHoltz filter mu ----
    int numLambda=201;
    Real lamMin=0., lamMax= 2*maxEig+20; // 60;
    RealArray lambda(numLambda);
    RealArray mu(numLambda);
    for( int i=0; i<numLambda; i++ )
    {
        lambda(i) = lamMin + (lamMax-lamMin)*(i-1.)/(numLambda-1.);
    }

  // -- Evaluate the beta function at all points in lambda(:) --
    bool useAdjusted= true; // !computeEigenmodes; 
    
    cgWave.getWaveHoltzIterationEigenvalue( lambda, mu, useAdjusted );


    ps.erase();

    RealArray points, value; 

  // --- Mark true eigenvalues if known ----
  // eigTrue(0:1,0:)
    const bool trueEigenPairsKnown = cgWave.dbase.has_key("uev");
    if( trueEigenPairsKnown )
    {
        RealArray & eigTrue                = cgWave.dbase.get<RealArray>("eig");
        const int numberOfEigenvectorsTrue = eigTrue.getLength(1);  
        

        if( numberOfEigenvectorsTrue>0 )
        {
      // for( int ie=0; ie<numberOfEigenvectorsTrue; ie++ )
      //   printF("plotfilter: ie=%3d eigTrue=%12.5e\n",ie,eigTrue(0,ie));

            int numToPlot=0;
            for( int ie=0; ie<numberOfEigenvectorsTrue; ie++ )
            {
                if( eigTrue(0,ie)>lamMax )
                    break;
                numToPlot=ie+1;
            }
            assert( numToPlot > 0 );

            RealArray lam(numToPlot);
            for( int ie=0; ie<numToPlot; ie++ )
            {
                lam(ie) = eigTrue(0,ie);
            }
            RealArray betaTrue(numToPlot);
            cgWave.getWaveHoltzIterationEigenvalue( lam, betaTrue, useAdjusted );  

            points.redim(numToPlot,2);
            value.redim(numToPlot); 
            for( int ie=0; ie<numToPlot; ie++ )
            {
                points(ie,0) = eigTrue(0,ie);
                points(ie,1) = betaTrue(ie);
        // printF("ie=%3d: true: eig=%12.4e adjusted=%12.4e beta=%12.4e\n",ie,eigTrue(0,ie),lam(ie),betaTrue(ie));
                value(ie)    = 0.2;                // colour point i by value(i)
            } 
      // --- plot locations of true eigenavleus 
            printF("plot numPlot=%d true values\n",numToPlot);
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);        
            psp.set(GI_POINT_SIZE,(real) 4.);  // size in pixels
            psp.set(GI_POINT_COLOUR,"BLACK");

            ps.plotPoints(points,psp); // colour point i by value(i)  

      // ps.plotPoints(points,value,psp); // colour point i by value(i)   
        } 
    }       

  // ---- mark eigenvalues found as points ----
    if( eigenValues.getLength(0)>0 )
    {
        RealArray beta(numEigenVectors);
        points.redim(numEigenVectors,2);
        value.redim(numEigenVectors);

        cgWave.getWaveHoltzIterationEigenvalue( eigenValues, beta, useAdjusted );

        value=.6; 
        for( int ie=0; ie<numEigenVectors; ie++ )
        {
            points(ie,0) = eigenValues(ie); // plot versus original
            points(ie,1) = beta(ie);

      // printF("ie=%3d eig=%12.4e beta=%12.4e omega=%12.4e\n",ie,eigenValues(ie),beta(ie),frequencyArrayAdjusted(0));
      // value(ie)    = 0.4 + .6* (ie+1.)/(numEigenVectors+1.); // point colour in [0,1]

            value(ie)    =  .2 + .8 * fabs(beta(ie)+.5)/(1.5);  // defines the colout of the point in [0,1]
        }

        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);        
        psp.set(GI_POINT_SIZE,(real) 8.);  // size in pixels

    // ps.plotPoints(points,value,psp); // colour point i by value(i)

    // psp.set(GI_POINT_COLOUR,"LIGHTBLUE");
    // psp.set(GI_POINT_COLOUR,"BLUE");
        psp.set(GI_POINT_COLOUR,"MEDIUMAQUAMARINE");
        ps.plotPoints(points,psp); 

        aString tName="omega";
        aString title=sPrintF("WaveHoltz Filter, omega=%6.2f",frequencyArray(0));
        aString xName[1];
        xName[0]="beta";

        psp.set(GI_PLOT_GRID_POINTS_ON_CURVES,false);
        psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);  
        
    // PlotIt::plot(ps,t,x,title,tName,xName,psp);
        #ifndef USE_PPP
            PlotIt::plot(ps,lambda,mu,title,tName,xName,psp);
        #else
            printF("FINISH ME FOR PARALLEL: PlotIt::plot(ps,lambda,mu,title,tName,xName,psp)\n");
        #endif

    }

    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);        

    if( !computeEigenmodes )
    {
    // -- reset ---
        frequencyArray = frequencyArraySave;
        periodArray    = periodArraySave;
        cgWave.dbase.get<real>("dt") = dtSave; 

    }
    printF("\n ^^^^^^ plotFilter: numEigenVectors=%d, dt=%9.3e, dtUsed=%9.3e, AFTER: frequencyArray(0)=%10.3e frequencyArrayAdjusted(0)=%10.3e\n",
        numEigenVectors, dtSave, dt, frequencyArray(0),frequencyArrayAdjusted(0) );


    return 0;
}

// ===============================================================================
// MACRO: read eigenvectors from the show file
// ===============================================================================


// ==============================================================================================
/// \brief Initialize Deflation for WaveHoltz 
// ==============================================================================================
int CgWave::initializeDeflation()
{
    bool & deflationInitialized = dbase.get<bool>("deflationInitialized"); //  set too true if deflation has been initialized
    if( deflationInitialized )
        return 0;

    Real cpuStart = getCPU();

    deflationInitialized=true;

    
    const int & orderOfAccuracy              = dbase.get<int>("orderOfAccuracy");
    const int & solveHelmholtz               = dbase.get<int>("solveHelmholtz");
    const int & computeTimeIntegral          = dbase.get<int>("computeTimeIntegral");
    const int & adjustOmega                  = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
    const int & adjustHelmholtzForUpwinding  = dbase.get<int>("adjustHelmholtzForUpwinding");
    const int & deflateWaveHoltz             = dbase.get<int>("deflateWaveHoltz");
    int & numToDeflate                       = dbase.get<int>("numToDeflate");
    const int & deflateForcing               = dbase.get<int>("deflateForcing");
    const int & onlyLoadDeflatedEigenVectors = dbase.get<int>("onlyLoadDeflatedEigenVectors");
    const aString & eigenVectorFile          = dbase.get<aString>("eigenVectorFile"); //  name of file holding eigs and eigenvectors for deflation
    const int & computeEigenmodes            = dbase.get<int>("computeEigenmodes");
    const int & numberOfFrequencies          = dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray         = dbase.get<RealArray>("frequencyArray");

    printF("\n ==== INITIALIZE DEFLATION (or read known eigenpairs for checking EigenWave) =====\n");
    printF("  deflateWaveHoltz=%d, deflateForcing=%d, numToDeflate=%d, computeEigenmodes=%d, eigenVectorFile=[%s]\n",
        deflateWaveHoltz,deflateForcing,numToDeflate,computeEigenmodes,(const char*)eigenVectorFile);


    if( !dbase.has_key("uev") )
    {
        dbase.put<realCompositeGridFunction>("uev");
        dbase.put<RealArray>("eig");

        dbase.put<IntegerArray>("eigenVectorIsNormalized");
        dbase.put<IntegerArray>("eigMultiplicity");
        dbase.put<IntegerArray>("eigNumbersToDeflate");

    // dbase.put<RealArray>("eigenVectorMaxNorm");

    }

    realCompositeGridFunction & uev        = dbase.get<realCompositeGridFunction>("uev");
    IntegerArray & eigenVectorIsNormalized = dbase.get<IntegerArray>("eigenVectorIsNormalized");
    IntegerArray & eigMultiplicity         = dbase.get<IntegerArray>("eigMultiplicity");
    IntegerArray & eigNumbersToDeflate     = dbase.get<IntegerArray>("eigNumbersToDeflate");
    const Real & eigTol                    = dbase.get<Real>("eigenValueTolForMultiplicity"); 
    int & eigenVectorsAreOrthogonal        = dbase.get<int>("eigenVectorsAreOrthogonal"); 

  // RealArray    & eigenVectorMaxNorm      = dbase.get<RealArray>("eigenVectorMaxNorm");


    int numberOfEigenvectors=-1; 
    RealArray & eig = dbase.get<RealArray>("eig");

  // Integrate integrate;

  // ---- Read show file with eigenvectors and process -----
  // readEigenvectors=true;

  // CompositeGrid & cgev = cg;  // ASSUME THIS FOR NOW -- NEED TO CHECK

  // -- this didn't seem to work: seg fault in cgwh when taking a maxNorm of uev ??
    if( !dbase.has_key("cgev") )
        dbase.put<CompositeGrid>("cgev"); // do this for now, save cgev in case we need the eigenvectors later for computing errors

    CompositeGrid & cgev = dbase.get<CompositeGrid>("cgev");

    printF("---- Read show file=[%s] with eigenvectors and process -----\n",(const char*)eigenVectorFile);
    ShowFileReader showFileReader;
    showFileReader.open(eigenVectorFile);

    int numFrames = showFileReader.getNumberOfFrames();
    const int numberOfSolutions = showFileReader.getNumberOfSolutions();
    printF("There are %d solutions and %d frames\n",numberOfSolutions,numFrames);
            

  // numberOfEigenvectors = uev.getComponentBound(0) - uev.getComponentBase(0) + 1;
    numberOfEigenvectors = numberOfSolutions;
    printF(">> There are %d eigenvectors\n",numberOfEigenvectors);
  // int solutionNumber=1;
  // showFileReader.getASolution(solutionNumber,cgev,uev);        // read in a grid and solution

    if( numToDeflate > numberOfEigenvectors )
    {
        printF("CgWave::initializeDeflation::ERROR: Requesting %d eigenvectors to deflate but only %d eigenvectors are in the show file.\n",numToDeflate,numberOfEigenvectors);
        OV_ABORT("ERROR");
    }


    int solutionNumber=1; // Frame 1 has some extra data 
    HDF_DataBase & db = *(showFileReader.getFrame(solutionNumber));
    db.get(eig,"eig"); 


  // -- genEigs now can orthogonalize the eigenvectors -- 
  //  *new* Sept 5, 2024 :

    eigenVectorsAreOrthogonal=false;

    IntegerArray eigStartIndex;
    int rt = db.get(eigMultiplicity,"eigMultiplicity");
    if( rt==0 )
    {
        eigenVectorsAreOrthogonal=true;
        printF("initDeflation: Eigenvectors in the show file have been orthogonalized and have L2-norm 1.\n");
        db.get(eigStartIndex,"eigStartIndex"); // currently not used : index of first of a multiple eigenvalue

    }
    else
    {
        printF("initDeflation: Eigenvectors in the show file have NOT been orthogonalized : this must be an old file. Regenerate with genEigs and orthogonalize to avoid doing this here.\n");
    
        eigMultiplicity.redim(numberOfSolutions);
        eigMultiplicity=1; // do this since we don't know 
    }

    if( false && eigenVectorsAreOrthogonal ) // done below now
    {
        if( !dbase.has_key("integrate") )
        {
      // FINISH ME -- first check if "integrate" is there
      //  db.find ...


            Integrate & integrate = dbase.put<Integrate>("integrate");
      // Integrate integrate(cg); // ************** TEST *****************

            printF("\n ###### CgWave::initDeflation: reading Integrate object from the show file with precomputed integration weights...\n");

            integrate.updateToMatchGrid(cg);  // This is correct -- see Overture/otherStuff/testIntegrate.C
            rt = integrate.get(db,"integrate");

            if( rt==1 )
            {
                printF("\n ###### CgWave::initDeflation: Integrate object WAS found in the show file.\n\n");
                Real volume;
                volume = integrate.volume();
                printf("initializeDeflation: computed volume of the domain = %9.2e\n\n",volume);
            }
            else
            {
                printF("\n ###### CgWave::initDeflation: Integrate object was NOT found in the show file.\n\n");
            }

      // if( rt==0 )
      //   printF("\n ###### CgWave::initDeflation: Integrate object WAS found in the show file.\n\n");
      // else
      //   printF("\n ###### CgWave::initDeflation: Integrate object was NOT found in the show file.\n\n");
        }
    }


  // --------- Eigenvalues for WaveHoltz are the square root of the input eigenvalues -----
    for( int i=0; i<numberOfEigenvectors; i++ )
    {
        if( fabs(eig(1,i)) < REAL_EPSILON*100*fabs(eig(0,i)) )
        {
      // Eigenvalue is real
            Real lambda = sqrt(eig(0,i));  // Take square root
            eig(0,i)=lambda; 
            if( debug>2 )
                printF("Eigenvalue %3d : k=%16.10f (after square root)\n",i,eig(0,i)); 
        }
        else
        {
            Real lambda = sqrt(eig(0,i));  // Take square root
            eig(0,i)=lambda; 
            printF("Eigenvalue %3d : k=%16.10f (after square root)\n",i,eig(0,i)); 

            printF("WARNING: Eigenvalue %3d : k=%16.10f + (%g) I, has IMAGINARY PART (IGNORING)\n",i,eig(0,i),eig(1,i));    
      // OV_ABORT("STOP HERE FOR NOW");    
        }
    }

  // ----- Find the range of eigenvalues in the show file and check that omega is contained in that range ----
    Range all;
    const Real lambdaMin = min( eig(0,all) );
    const Real lambdaMax = max( eig(0,all) );
    const Real lambdaAve = sum( eig(0,all) )/numberOfEigenvectors; 
    printF("CgWave::initializeDeflation:: Eigenvalues from showfile: lambdaMin=%9.2e, lambdaAve=%9.2e, lambdaMax=%9.2e, omega=%9.2e\n",
                    lambdaMin,lambdaAve,lambdaMax,frequencyArray(0));
    if( !computeEigenmodes && (frequencyArray(0) < lambdaMin || frequencyArray(0) > lambdaMax) )
    {
        printF("CgWave::initializeDeflation::ERROR omega=%10.2e is outside the range of eigenvalues [%9.2e,%9.2e].\n"
                      "                                   Deflation will not work very well. Stopping here for now.\n"
                      ,frequencyArray(0),lambdaMin,lambdaMax);
        OV_ABORT("ERROR");
    }


    if( FALSE && !eigenVectorsAreOrthogonal )
    {
    // ***** OLD STUFF -- this should be done in genEigs now ****
    // ---- compute multiplicities ---

        eigenVectorIsNormalized.redim(numberOfEigenvectors); eigenVectorIsNormalized=0;
        eigMultiplicity.redim(numberOfEigenvectors);         eigMultiplicity=1;
    // eigenVectorMaxNorm.redim(numberOfEigenvectors);      eigenVectorMaxNorm=-1.;  // -1 means max-norm 

    // ---- WE COULD HAVE genEigs save the multiplicity ****** FIX ME ****
        for( int i=0; i<numberOfEigenvectors; i++ )
        {  

      // bool doubleEig=false;
      // const Real eigTol = 1.e-4; // ** FIX ME 
            int ie=i; // first eig of multiple set 
            int je=i; // last eig of multiple set
            const int maxMultiplicity=10; 
            int multiplicity=1; 
            for( int m=0; m<maxMultiplicity; m++ )
            {
                bool matchFound=false;

        // if( ie>0 )
        //   printF("ie=%d, eig(0,ie-1)=%.4g\n",ie, eig(0,ie-1));

                if( ie>0 && fabs(eig(0,ie-1)-eig(0,i))< eigTol*(1.+fabs(eig(0,i))) ) 
                {
                    multiplicity++;
          // printF("MATCH FOUND: multiplicity=%d\n",multiplicity);
                    ie--; 
                    matchFound=true;
                }
        // if( je<numberOfEigenvectors-1)
        //   printF("je=%d, eig(0,je+1)=%.4g\n",je, eig(0,je+1));

                if( (je<numberOfEigenvectors-1) && fabs(eig(0,je+1)-eig(0,i)) < eigTol*(1.+fabs(eig(0,i))) )
                {
                    multiplicity++;
          // printF("MATCH FOUND: multiplicity=%d\n",multiplicity);
                    je++; 
                    matchFound=true;
          // doubleEig=true;
          // j=i+1; 
                }
                if( !matchFound ) break;
            }
            eigMultiplicity(i) = multiplicity;
            if( debug>2 )
                printF(" i=%3d : eig=%14.6e, multiplicity=%d (eigTol=%9.2e)\n",i,eig(0,i),eigMultiplicity(i),eigTol);


            if( false )
            {
                Real evNorm = maxNorm( uev,i);
                printF("maxNorm( uev[%d] ) = %9.2e\n",i,evNorm);
            }

        }
    }

  // ::display(eigMultiplicity,"eigMultiplicity","%2d ");
  // OV_ABORT("stop here for now");

 //  showFileReader.close();



  // **NEW WAY** eigenvectors are orthonormal  Sept 5, 2024

  // ----------------------------------------------------------------------
  // ---------- MAKE A LIST OF EIGENVECTORS TO USE FOR DEFLATION ----------
  // ----------------------------------------------------------------------

  // const int & numberOfFrequencies         = dbase.get<int>("numberOfFrequencies");
  // const RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
    const RealArray & frequencyArraySave    = dbase.get<RealArray>("frequencyArraySave");
    const RealArray & periodArray           = dbase.get<RealArray>("periodArray");  

    eigenVectorIsNormalized.redim(numberOfEigenvectors);   
    eigenVectorIsNormalized=1; 

    eigNumbersToDeflate.redim(numToDeflate*2+10);

    RealArray beta(numberOfEigenvectors);
    IntegerArray iperm(numberOfEigenvectors); // for sorting 

    if( !computeEigenmodes )
    {
    // **** Choose eigenvalues to deflate based on the theoretical convergence rate ****

        RealArray lamv(numberOfEigenvectors);
        for( int i=0; i<numberOfEigenvectors; i++ )
        {
            lamv(i) = eig(0,i);
        }

    // -----------------------------------------------------------------------------------------------
    // Evaluate the "beta" or "mu" function that determines the convergence of the WaveHoltz iteration
    // -----------------------------------------------------------------------------------------------
        getWaveHoltzIterationEigenvalue( lamv, beta );

        RealArray betaSave; // save for below
        betaSave= beta;

        for( int i=0; i<numberOfEigenvectors; i++ )
        {
            Real lambda = lamv(i);
            iperm(i)=i; // for sorting 

            beta(i) = fabs( beta(i) );  // NOTE ABSOLUTE VALUE 

            if( debug>2 )
                printF(" WaveHoltz filter values: i=%4d: lambda=%12.4e, |beta|=%12.4e\n",i,lambda,beta(i));
        }
    // -- bubble sort ---
        for( int i=0; i<numberOfEigenvectors-1; i++ )
        {
            bool changed=false;
            for( int j=0; j<numberOfEigenvectors-1; j++ )
            {
                if( beta(j) < beta(j+1) )
                {
                    changed=true;
                    Real temp=beta(j);   beta(j)= beta(j+1);  beta(j+1)=temp;
                    int itemp=iperm(j); iperm(j)=iperm(j+1); iperm(j+1)=itemp;
                }
            }
            if( !changed ) break;
        }

        if( debug>2 )
            printF("\n ---- SORTED by beta: freq(0)=%12.4e\n",frequencyArray(0));
        for( int j=0; j<numberOfEigenvectors; j++ )
        {
            int i = iperm(j);  // original ordering 
            Real lambda = eig(0,i);
            if( debug>2 )
                printF(" i=%4d: lambda=%12.4e, |beta|=%12.4e\n",i,lambda,beta(j));
        }

        if( debug>2 )
            printF("\n ---- CHOOSE EIGENVALUES TO DEFLATE: (freq(0)=%12.4e) \n",frequencyArray(0));

        int ied = 0; // counts eigenvectors to deflate
        int j=0; 
        for( j=0; j<numToDeflate && j<numberOfEigenvectors-1; j++ )
        {
            int eigNumber = iperm(j);  // original ordering 
            if( ied>=eigNumbersToDeflate.getLength(0) )
            {
        // increase size
                int newSize = eigNumbersToDeflate.getLength(0)*2+10;
                printF("Increasing size of eigNumbersToDeflate from %d to %d (must be large multiplicities)\n",eigNumbersToDeflate.getLength(0),newSize);

        // ::display(eigNumbersToDeflate,"eigNumbersToDeflate before resize");
                eigNumbersToDeflate.resize(newSize);
        // ::display(eigNumbersToDeflate,"eigNumbersToDeflate after resize");
            }
            eigNumbersToDeflate(ied)=eigNumber;  ied++;


            Real lambda = eig(0,eigNumber);
            
            if( true || debug>2 )
                printF("deflate ied=%3d eigeNumber=%d: lambda=%12.4e, (multiplicity=%d) beta=%12.4e\n",ied,eigNumber,lambda,eigMultiplicity(eigNumber),beta(j));       

      // normalizeEigenvector( eigNumber );

            if( eigMultiplicity(eigNumber)>1 && j<numberOfEigenvectors-1 )
            {
        // Include the multiple eigenvalue 
        // assert( eigMultiplicity(eigNumber) );
                int eigMult = eigMultiplicity(eigNumber);
                for( int k=1; k<eigMult && j<numberOfEigenvectors-1 ; k++ )
                {
                    j++; numToDeflate++;
                    eigNumber = iperm(j);  // original ordering 
                    if( ied>=eigNumbersToDeflate.getLength(0) )
                    {
            // increase size
                        int newSize = eigNumbersToDeflate.getLength(0)*2+10;
                        printF("Increasing size of eigNumbersToDeflate from %d to %d (must be large multiplicities)\n",eigNumbersToDeflate.getLength(0),newSize);
                        eigNumbersToDeflate.resize(newSize);
                    }          
                    eigNumbersToDeflate(ied)=eigNumber;  ied++;

                    Real lambda2 =  eig(0,eigNumber);
                    if( 1==1 || debug>2 )
                        printF("deflate ied=%3d eigNumber=%d: lambda=%12.4e, (multiplicity=%d) *** duplicate ***\n",ied,eigNumber,lambda2,eigMultiplicity(eigNumber)); 

                }

            }

        } // end for j
        int numToDeflateNew =min(j,numberOfEigenvectors);

        printF("ACTUAL numToDeflate=%d, numberOfEigenvectors=%d .. using numToDeflate=%d\n",numToDeflate,numberOfEigenvectors,numToDeflateNew);
        numToDeflate=numToDeflateNew;

    // --- Save beta eigenvalues for augmented GMRES ----
        RealArray & betaDeflate = dbase.get<RealArray>("betaDeflate"); 
        betaDeflate.redim(numToDeflate);
        for( int ied=0; ied<numToDeflate; ied++ )
        {
            if( eigNumbersToDeflate(ied)>=numberOfEigenvectors || eigNumbersToDeflate(ied)>=betaSave.getLength(0) )
            {
                printF("ERROR: Setting betaDeflate but ied=%d, eigNumbersToDeflate(ied)=%d, numberOfEigenvectors=%d, betSave.getLength(0)=%d\n",
                        ied,eigNumbersToDeflate(ied),numberOfEigenvectors,betaSave.getLength(0));
                OV_ABORT("ERROR");
            }

            betaDeflate(ied) = betaSave(eigNumbersToDeflate(ied)); // check me 
        }

        Real & waveHoltzAsymptoticConvergenceRate = dbase.get<real>("waveHoltzAsymptoticConvergenceRate");
        waveHoltzAsymptoticConvergenceRate = beta(numToDeflate);

    // numToDeflate=ied;
        printF(">>> number of eigenvectors to deflate (including multiplicities) = %d\n",numToDeflate);
        printF(">>> Estimated asymptotic convergence rate=%6.2f\n",waveHoltzAsymptoticConvergenceRate);



        if( dbase.get<int>("plotFilterAndDeflatedEigenvalues") ) 
        {
      // ---- plot the WaveHoltz filter beta and deflated eigs for checking ----
            RealArray eigenValues(numToDeflate);
            for( int ie=0; ie<numToDeflate; ie++ )
                eigenValues(ie) = eig(0,eigNumbersToDeflate(ie));  // array of deflated eigs

      // ::display(eigenValues,"eigs to deflate for plotFilter");
            plotFilter( eigenValues );

        }
    }

  // --------- NOW READ EIGENVECTORS FROM THE SHOW FILE -------
    // Eigenvectors are stored as time-steps (in different frames)
        printF("======= READ EIGENVECTORS FROM THE SHOW FILE onlyLoadDeflatedEigenvectors=%d ========\n",onlyLoadDeflatedEigenVectors);
        bool gridsMatch=true; // set to false if the grid in the eigenvector file does not match cg 
        realCompositeGridFunction q;
        InterpolatePointsOnAGrid interp; // for interpolating from a grid with a different resolution
        const int numToRead = onlyLoadDeflatedEigenVectors ? numToDeflate : numberOfSolutions; 
        for( int ie=0; ie<numToRead; ie++ )
        {
      // --- read in each EV and save in uev ---
            int je = onlyLoadDeflatedEigenVectors ? eigNumbersToDeflate(ie) : ie; 
            solutionNumber= je+1; 
            printF("Read solutionNumber=%d from the show file (numberOfEigenvectors=%d) \n",solutionNumber,numberOfEigenvectors);
            showFileReader.getASolution(solutionNumber,cgev,q);
            if( ie==0 )
            {
                Range all;
        // uev.updateToMatchGrid(cgev,all,all,all,numberOfSolutions);
                uev.updateToMatchGrid(cg,all,all,all,numberOfSolutions);
                gridsMatch = CgWave::compositeGridsMatch( cg, cgev );
                if( gridsMatch )
                    printF("CgWaveHoltz::initializeDeflation: read eigenvectors from a file: Grid in eigenvector file grid seems to match with current grid\n");
                else
                    printF("CgWaveHoltz::initializeDeflation: read eigenvectors from a file: Grid in eigenvector file grid DOES NOT match with current grid... will interpolate...\n");        
            } 
            Index I1,I2,I3;
            if( gridsMatch )
            {
        // --- just copy grid function from EV file ---
                for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cgev[grid];
                    getIndex(mg.dimension(),I1,I2,I3);
                    OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                    OV_GET_SERIAL_ARRAY(real,q[grid],qLocal);        
                    bool ok=ParallelUtility::getLocalArrayBounds(q[grid],qLocal,I1,I2,I3);
                    if( ok )
                        uevLocal(I1,I2,I3,ie) = qLocal(I1,I2,I3);
                } 
            }
            else
            {
        // --- Interpolate solution from the EV file : this could be a coarser grid version used for deflation with augmented GMRES
                printF("CgWaveHoltz::initializeDeflation: interpolate the eigenvector ie=%d to the current grid ...\n",ie);
                CompositeGrid & cgev = *q.getCompositeGrid();
                cgev.update(MappedGrid::THEmask );
                interp.setAssignAllPoints(true);  // assign all points -- extrap if necessary
                int width = orderOfAccuracy+1;
        // width = orderOfAccuracy+4; // *** TEST ***;
                interp.setInterpolationWidth( width ); // set interp width 
                Range C(0,0);
                interp.interpolateAllPoints( q,uev, C, Range(ie,ie) );  // interpolate v from q
            }     
        }


  // eigenvectorsAreTrueEigenvectors : used by AugmentedGmres 
    int & eigenvectorsAreTrueEigenvectors = dbase.get<int>("eigenvectorsAreTrueEigenvectors");
    if( gridsMatch )
        eigenvectorsAreTrueEigenvectors=true;
    else
        eigenvectorsAreTrueEigenvectors=false;
    
    cgev.update(MappedGrid::THEmask ); // *wdh* March 28, 2023

  if( eigenVectorsAreOrthogonal )
    {
        if( !dbase.has_key("integrate") )
        {
      // FINISH ME -- first check if "integrate" is there
      //  db.find ...
            const int solutionNumber=1; // Frame 1 has some extra data 
            HDF_DataBase & db = *(showFileReader.getFrame(solutionNumber));

            Integrate & integrate = dbase.put<Integrate>("integrate");
      // Integrate integrate(cg); // ************** TEST *****************

            printF("\n ###### CgWave::initDeflation: reading Integrate object from the show file with precomputed integration weights...\n");

            integrate.updateToMatchGrid(cg);  // This is correct -- see Overture/otherStuff/testIntegrate.C
            if( eigenvectorsAreTrueEigenvectors )
            {
        // --- we can use the Integrate object in the show file if we have fine grid eigenvectors ----
                rt = integrate.get(db,"integrate");
                if( rt==1 )
                {
                    printF("\n ###### CgWave::initDeflation: Integrate object WAS found in the show file.\n\n");
                    Real volume;
                    volume = integrate.volume();
                    printf("initializeDeflation: computed volume of the domain = %9.2e\n\n",volume);
                }
                else
                {
                    printF("\n ###### CgWave::initDeflation: Integrate object was NOT found in the show file.\n\n");
                }
            }
            else
            {
                Real volume;
                volume = integrate.volume();
                printf("\ninitializeDeflation: computed volume of the domain = %9.2e (NOT using Integrate from show file, since eigenvectors are approximate)\n\n",volume);
            }
      // if( rt==0 )
      //   printF("\n ###### CgWave::initDeflation: Integrate object WAS found in the show file.\n\n");
      // else
      //   printF("\n ###### CgWave::initDeflation: Integrate object was NOT found in the show file.\n\n");
        }

    }

    showFileReader.close();


  // ---- DEFLATE FORCING HERE -----
  // NOTE: we currently do NOT deflate the forcing when using Augmented Krylov methods
    const int & useAugmentedGmres = dbase.get<int>("useAugmentedGmres"); 
    if( !computeEigenmodes && (deflateForcing || true) && useAugmentedGmres==0 )
    {
    // deflate forcing 
        int deflateOption=1; // deflate forcing 
        deflateSolution( deflateOption );

        if( false )
        {  // -- plot deflated forcing ---
              realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f");      
              GenericGraphicsInterface & ps = gi;
              PlotStuffParameters & psp =  dbase.get<PlotStuffParameters>("psp");

              psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
              psp.set(GI_TOP_LABEL,"deflated forcing");
              ps.erase();
              PlotIt::contour(ps,f,psp);     
              psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);


        }
    }



    Real cpuDeflationSetup = getCPU() - cpuStart;
    timing(timeForDeflation) += cpuDeflationSetup;

    printF("\n @@@@@@@@@@@  DEFLATION SETUP: cpu = %9.2e(s) @@@@@@@@@@@@@\n\n",cpuDeflationSetup);
    return 0;
}

// ==============================================================================================
/// \brief Make various checks for deflation
// ==============================================================================================
int CgWave::
checkDeflation()
{
    OV_ABORT("checkDeflation IS BROKEN -- probably no longer needed");

    const int & solveHelmholtz              = dbase.get<int>("solveHelmholtz");
    const int & computeTimeIntegral         = dbase.get<int>("computeTimeIntegral");
    const int & adjustOmega                 = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
    const int & adjustHelmholtzForUpwinding = dbase.get<int>("adjustHelmholtzForUpwinding");
    const int & deflateWaveHoltz            = dbase.get<int>("deflateWaveHoltz");
    const int & numToDeflate                = dbase.get<int>("numToDeflate");
    const aString & eigenVectorFile         = dbase.get<aString>("eigenVectorFile"); //  name of file holding eigs and eigenvectors for deflation
    const Real eigTol                       = dbase.get<Real>("eigenValueTolForMultiplicity");
    printF("\n ==== CHECK DEFLATION =====\n");
    printF("  deflateWaveHoltz=%d, numToDeflate=%d, eigenVectorFile=[%s]\n",deflateWaveHoltz,numToDeflate,(const char*)eigenVectorFile);

    if( !dbase.has_key("uev") )
    {
        initializeDeflation();
    }

    realCompositeGridFunction & uev = dbase.get<realCompositeGridFunction>("uev");
    RealArray & eig = dbase.get<RealArray>("eig");
    int numberOfEigenvectors = uev.getComponentBound(0) - uev.getComponentBase(0) + 1; 

    CompositeGrid & cgev = cg;


    Integrate & integrate = dbase.get<Integrate>("integrate");

  // -- holds intermediate results: 
    realCompositeGridFunction & u = dbase.get<realCompositeGridFunction>("uDeflate");
    u=0.; 

  // u=1.;
  // Real volume;
  // volume = integrate.volumeIntegral(u);
  // printf("Volume of domain = %9.2e\n",volume);


    Real dotProduct;
    Index I1,I2,I3;


  // ***DO THIS FOR NOW***

  // -- normalize and orthogonalize eigenvectors ----

    for( int i=0; i<numberOfEigenvectors; i++ )
    {
            for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cgev[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
                bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3);
                if( ok )
                    uLocal(I1,I2,I3) = uevLocal(I1,I2,I3,i)*uevLocal(I1,I2,I3,i);
            }
            Real eNorm = sqrt( integrate.volumeIntegral(u) );
            for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cgev[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
                if( ok )
                    uevLocal(I1,I2,I3,i) *= (1./eNorm);
            }

        if( i>0 )
        {
      // --- Look for multiple eigenvalues ---
            int j=i-1; 
            const Real delta = fabs(eig(0,i)-eig(0,j))/fabs(eig(0,i));

            if( i<100 )
                printF(" i=%d: eig(i)=%9.3e eig(i-1)=%9.3e, delta=%9.3e\n",i,eig(0,i),eig(0,j),delta);

      // const Real eigTol = 1.e-4; // ** FIX ME 
            if( delta < eigTol )
            {
        // --- We have a mutiple eigenvalue -----
        // --- Orthogonalize the eigenvectors ----

                if( i<100 )
                    printF(" >>> Multiple eigenvalue found: i=%d, eig=%10.3e, will orthogonalize eigenvectors...\n",i,eig(0,i));

                    for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
                    {
                        MappedGrid & mg = cgev[grid];
                        getIndex(mg.gridIndexRange(),I1,I2,I3);
                        OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                        OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
                        bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3);
                        if( ok )
                            uLocal(I1,I2,I3) = uevLocal(I1,I2,I3,i)*uevLocal(I1,I2,I3,j);
                    }
                    dotProduct = integrate.volumeIntegral(u);

                if( i<100 )
                    printF(" Inner product: (u%d,u%d) = %9.3e\n",i,j,dotProduct);

        // -- Gram-Schmidt --
        //   ui = ui - (ui,uj)*uj 
                for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cgev[grid];
                    getIndex(mg.dimension(),I1,I2,I3);
                    OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                    bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
                    if( ok )
                        uevLocal(I1,I2,I3,i) -= dotProduct*uevLocal(I1,I2,I3,j);
                }

        // re-normalize

                    for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
                    {
                        MappedGrid & mg = cgev[grid];
                        getIndex(mg.dimension(),I1,I2,I3);
                        OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                        OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
                        bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3);
                        if( ok )
                            uLocal(I1,I2,I3) = uevLocal(I1,I2,I3,i)*uevLocal(I1,I2,I3,i);
                    }
                    Real eNorm = sqrt( integrate.volumeIntegral(u) );
                    for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
                    {
                        MappedGrid & mg = cgev[grid];
                        getIndex(mg.dimension(),I1,I2,I3);
                        OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                        bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
                        if( ok )
                            uevLocal(I1,I2,I3,i) *= (1./eNorm);
                    }
            }
        }


    }

  // OV_ABORT("checkDeflation: stop here for now");
    return 0;
}


// ==============================================================================================
/// \brief Normalize eigenvector number "eigNumber". Orthogonalize any eigenvectors corresponding to
///      multiple eigenvalues.
/// \return value: multipilicity of the eigenvector: 1=simple, 2=double, ...
// ==============================================================================================
int CgWave::
normalizeEigenvector( int eigNumber )
{
    const int & debug                       = dbase.get<int>("debug");
    const int & solveHelmholtz              = dbase.get<int>("solveHelmholtz");
    const int & computeTimeIntegral         = dbase.get<int>("computeTimeIntegral");
    const int & adjustOmega                 = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
    const int & adjustHelmholtzForUpwinding = dbase.get<int>("adjustHelmholtzForUpwinding");
    const int & deflateWaveHoltz            = dbase.get<int>("deflateWaveHoltz");
    const int & numToDeflate                = dbase.get<int>("numToDeflate");
    const aString & eigenVectorFile         = dbase.get<aString>("eigenVectorFile"); //  name of file holding eigs and eigenvectors for deflation

    IntegerArray & eigenVectorIsNormalized = dbase.get<IntegerArray>("eigenVectorIsNormalized");
    IntegerArray & eigMultiplicity         = dbase.get<IntegerArray>("eigMultiplicity");
    const Real eigTol                      = dbase.get<Real>("eigenValueTolForMultiplicity");


    if( eigenVectorIsNormalized(eigNumber)  )
    {
        printF(">> normalizeEigenvector: INFO: eigenvector=%d already normalized\n",eigNumber);
        return eigMultiplicity(eigNumber);
    }
    int multiplicity=1; 

    if( debug>2 )
    {
        printF("\n ==== NORMALIZE EIGENVECTOR eigNumber=%d =====\n",eigNumber);
        printF("  deflateWaveHoltz=%d, numToDeflate=%d, eigenVectorFile=[%s]\n",deflateWaveHoltz,numToDeflate,(const char*)eigenVectorFile);
    }

    if( !dbase.has_key("uev") )
    {
        initializeDeflation();
    }

    realCompositeGridFunction & uev = dbase.get<realCompositeGridFunction>("uev");
    RealArray & eig = dbase.get<RealArray>("eig");
    int numberOfEigenvectors = uev.getComponentBound(0) - uev.getComponentBase(0) + 1; 

    CompositeGrid & cgev = cg;


    Integrate & integrate = dbase.get<Integrate>("integrate");

  // -- holds intermediate results: 
    realCompositeGridFunction & u = dbase.get<realCompositeGridFunction>("uDeflate");

    Real dotProduct;
    Index I1,I2,I3;


  // ***DO THIS FOR NOW***

  // -- normalize and orthogonalize eigenvectors ----

    int i = eigNumber;

        for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cgev[grid];
            getIndex(mg.dimension(),I1,I2,I3);
            OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
            OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
            bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3);
            if( ok )
                uLocal(I1,I2,I3) = uevLocal(I1,I2,I3,i)*uevLocal(I1,I2,I3,i);
        }
        Real eNorm = sqrt( integrate.volumeIntegral(u) );
        for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cgev[grid];
            getIndex(mg.dimension(),I1,I2,I3);
            OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
            bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
            if( ok )
                uevLocal(I1,I2,I3,i) *= (1./eNorm);
        }

    eigenVectorIsNormalized(i)=1;

    if( eigMultiplicity(i) > 1 )
    {
    // ---  multiple eigenvalue ----
        multiplicity = eigMultiplicity(i);

        printF(" >>> Multiple eigenvalue found: i=%d, eig=%12.5e, eigMultiplicity(i)=%d will orthogonalize eigenvectors...\n",i,eig(0,i),eigMultiplicity(i));

        int j=i;
    
        i=j+1; // check for duplicate eig here
        Real delta = fabs(eig(0,i)-eig(0,j))/fabs(eig(0,i));
    // Real eigTol = 1.e-4; // ** FIX ME 
        if( delta> eigTol && i>0 )
        {
            i=j-1;
            delta = fabs(eig(0,i)-eig(0,j))/fabs(eig(0,i));
        }

        assert( eigMultiplicity(i) > 1 );
        if( eigenVectorIsNormalized(i) != 0 )
        {
            printF(">>> ERROR: eigenVectorIsNormalized(i)=%d, expecting zero. i=%d, j=%d, eigMultiplicity(i)=%d\n",eigenVectorIsNormalized(i),i,j,eigMultiplicity(i));
            OV_ABORT("error");
        }

    // Real delta = fabs(eig(0,i)-eig(0,j))/fabs(eig(0,i));
    // Real eigTol = 1.e-3; // ** FIX ME 
    // assert( delta < eigTol );

    // --- Orthogonalize the eigenvectors ----

            for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cgev[grid];
                getIndex(mg.gridIndexRange(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
                bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3);
                if( ok )
                    uLocal(I1,I2,I3) = uevLocal(I1,I2,I3,i)*uevLocal(I1,I2,I3,j);
            }
            dotProduct = integrate.volumeIntegral(u);

        if( i<100 )
            printF(" Inner product: (u%d,u%d) = %9.3e\n",i,j,dotProduct);

    // -- Gram-Schmidt --
    //   ui = ui - (ui,uj)*uj 
        for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cgev[grid];
            getIndex(mg.dimension(),I1,I2,I3);
            OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
            bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
            if( ok )
                uevLocal(I1,I2,I3,i) -= dotProduct*uevLocal(I1,I2,I3,j);
        }

    // re-normalize
        

            for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cgev[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
                bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3);
                if( ok )
                    uLocal(I1,I2,I3) = uevLocal(I1,I2,I3,i)*uevLocal(I1,I2,I3,i);
            }
            Real eNorm = sqrt( integrate.volumeIntegral(u) );
            for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cgev[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
                if( ok )
                    uevLocal(I1,I2,I3,i) *= (1./eNorm);
            }
        eigenVectorIsNormalized(i)=1;    

        printF("Eigenvector %d is now orthonormal to eigenvector %d.\n",i,j);

    }

  // if( i>0 )
  // {
  //   // --- Look for multiple eigenvalues ---
  //   int j=i-1; 
  //   const Real delta = fabs(eig(0,i)-eig(0,j))/fabs(eig(0,i));

  //   if( i<100 )
  //     printF(" i=%d: eig(i)=%9.3e eig(i-1)=%9.3e, delta=%9.3e\n",i,eig(0,i),eig(0,j),delta);

  //   const Real eigTol = 1.e-3; // ** FIX ME 
  //   if( delta < eigTol )
  //   {
  //     // --- We have a mutiple eigenvalue -----
  //     // --- Orthogonalize the eigenvectors ----

  //     multiplicity=2;

  //     if( i<100 )
  //       printF(" >>> Multiple eigenvalue found: i=%d, eig=%10.3e, will orthogonalize eigenvectors...\n",i,eig(0,i));

  //     innerProductor( i,j,dotProduct );  

  //     if( i<100 )
  //       printF(" Inner product: (u%d,u%d) = %9.3e\n",i,j,dotProduct);

  //     // -- Gram-Schmidt --
  //     //   ui = ui - (ui,uj)*uj 
  //     for( int grid=0; grid<cgev.numberOfComponentGrids(); grid++ )
  //     {
  //       MappedGrid & mg = cgev[grid];
  //       getIndex(mg.dimension(),I1,I2,I3);
  //       OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);

  //       uevLocal(I1,I2,I3,i) -= dotProduct*uevLocal(I1,I2,I3,j);
  //     }

  //     // re-normalize

  //     normalizeEigenvectorMacro(i);             
  //   }
  // }

    return multiplicity;
}


// ==============================================================================================
/// \brief Deflate the current WaveHoltz solution (or the forcing, if deflateHelmholtzForcing==1) 
///   by subtracting out the components along the
///    selected eigenvectors.
/// \param deflateOption : 
///        0 = default
///        1 = deflate forcing 
///        2 = deflate solution even if deflateForcing==1 (for initial condition)
// ==============================================================================================
int CgWave::
deflateSolution( int deflateOption /* = 0 */  )
{
    Real cpuStart = getCPU();
    const bool useOpt=true;  // use optimized loops

    bool & deflationInitialized = dbase.get<bool>("deflationInitialized");
    const int & debug           = dbase.get<int>("debug"); 
    if( !deflationInitialized )
    {
        initializeDeflation();
    // checkDeflation();   // This will normalized eigenvectors ****** FIX ME *******
    }

    const int & deflateForcing               = dbase.get<int>("deflateForcing");

  // Do not deflate the solution if we have deflated the forcing
    if( deflateForcing )
    {
        if( deflateOption==0 )
        {
            printF("CgWave::deflateSolution: do NOT deflate solution since forcing was deflated.\n");
            return 0;
        }
    // else if( deflateOption==1 )
    // {
    //   printF("\n ######### CgWave::deflateSolution: DEFLATE FORCING ##############\n\n");      
    // }
        else if( deflateOption==2 )
        {
            printF("CgWave::deflateSolution: deflateForcing==1 BUT deflate the solution (must be the initial condition)\n");  
        }
        else
        {
            OV_ABORT("CgWave::deflateSolution: ERROR - unknown deflateOption");
        }
    }
    if( deflateOption==1 )
    {
        printF("\n ######### CgWave::deflateSolution: DEFLATE FORCING ##############\n\n");      
    }  

    const int & eigenVectorsAreOrthogonal  = dbase.get<int>("eigenVectorsAreOrthogonal"); 
    if( !eigenVectorsAreOrthogonal )
    {
          printF("\nCgWave::deflateSolution:ERROR: eigenvectors for deflation have not been orthogonalized.\n"
                        "  Use 'genEigs' to orthogonalize the eigenvectors in the show file.\n"
                        "  **OR** use the augmented GMRES algorithm which does not require orthogonal eigenvectors\n\n"
                        );
          OV_ABORT("ERROR");

    }

    const int & numberOfFrequencies          = dbase.get<int>("numberOfFrequencies");
    const int & deflateWaveHoltz             = dbase.get<int>("deflateWaveHoltz");
    const int & numToDeflate                 = dbase.get<int>("numToDeflate");
    const IntegerArray & eigNumbersToDeflate = dbase.get<IntegerArray>("eigNumbersToDeflate");
    const int & onlyLoadDeflatedEigenVectors = dbase.get<int>("onlyLoadDeflatedEigenVectors");

    if( deflateWaveHoltz==0 || numToDeflate<=0 )
        return 0;

  // -- holds intermediate results: 
  // uDeflate : holds intermediate results for deflation:
    if( !dbase.has_key("uDeflate") )
    {
        realCompositeGridFunction & u = dbase.put<realCompositeGridFunction>("uDeflate");
        u.updateToMatchGrid(cg);
    }  
    realCompositeGridFunction & u = dbase.get<realCompositeGridFunction>("uDeflate");

  // -- eigenvectors:
    realCompositeGridFunction & uev = dbase.get<realCompositeGridFunction>("uev");

  // WaveHoltz solution or forcing:
    realCompositeGridFunction & v = deflateOption==1 ? dbase.get<realCompositeGridFunction>("f") : 
                                                                                                          dbase.get<realCompositeGridFunction>("v");

    if( !dbase.has_key("integrate") )
    {
        printF("==== CgWave::deflateSolution: compute integration weights...\n");
        Integrate & integrate = dbase.put<Integrate>("integrate");
        integrate.updateToMatchGrid(cg);    
    }
    Integrate & integrate = dbase.get<Integrate>("integrate");


  // int freq=0;      // adjust this frequency **FIX ME**

  

    Index I1,I2,I3;
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
        for( int ied=0; ied<numToDeflate; ied++ )
        {
            const int eigNumber = onlyLoadDeflatedEigenVectors ? ied : eigNumbersToDeflate(ied); // deflate this eigenvector 
            if( debug>1  )
                printF("deflateSolution: freq=%d : deflate eigenvector %d (numToDeflate=%d)\n",freq,eigNumber,numToDeflate);

      // --- Compute the inner product (phi,v)_h ----
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                OV_GET_SERIAL_ARRAY(real,u[grid],uLocal); // temp space
                OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
                bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
                if( ok )
                {
                    if( useOpt )
                    {
                        Real *pu         =   uLocal.getDataPointer();
                        const Real *pv   =   vLocal.getDataPointer();
                        const Real *puev = uevLocal.getDataPointer();
                        const int i1Base=uLocal.getBase(0), i2Base=uLocal.getBase(1), i3Base=uLocal.getBase(2);
                        const int nd1= uLocal.getLength(0), nd2= uLocal.getLength(1), nd3= uLocal.getLength(2);
                        #define ua(i1,i2,i3)       pu[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base          ))]
                        #define va(i1,i2,i3,n)     pv[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base + nd3*(n)))]
                        #define ueva(i1,i2,i3,n) puev[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base + nd3*(n)))]
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            ua(i1,i2,i3) = ueva(i1,i2,i3,eigNumber)*va(i1,i2,i3,freq);
                        }
                    }
                    else
                    {
                        uLocal(I1,I2,I3) = uevLocal(I1,I2,I3,eigNumber)*vLocal(I1,I2,I3,freq);
                    }
                }
            }

            Real phiDotv = integrate.volumeIntegral(u);

            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
                bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
                if( ok )
                {
                    if( useOpt )
                    {
                        Real *pv         =   vLocal.getDataPointer();
                        const Real *puev = uevLocal.getDataPointer();
                        const int i1Base=vLocal.getBase(0), i2Base=vLocal.getBase(1), i3Base=vLocal.getBase(2);
                        const int nd1= vLocal.getLength(0), nd2= vLocal.getLength(1), nd3= vLocal.getLength(2);
            // #define va(i1,i2,i3,n)     pv[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base + nd3*(n)))]
            // #define ueva(i1,i2,i3,n) puev[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base + nd3*(n)))]
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            va(i1,i2,i3,freq) -= phiDotv*ueva(i1,i2,i3,eigNumber);
                        }
                    }
                    else
                    {
                        vLocal(I1,I2,I3,freq) -= phiDotv*uevLocal(I1,I2,I3,eigNumber);
                    }

                }
            }

        }
    }
    #undef ua
    #undef va
    #undef ueva
    timing(timeForDeflation) += getCPU()-cpuStart;

    return 0;
}


// ==============================================================================================
/// \brief Inflate the WaveHoltz solution by adding in the solution along the
///    selected eigenvectors.
// ==============================================================================================
int CgWave::
inflateSolution()
{
    Real cpuStart = getCPU();
    const bool useOpt=true;  // use optimized loops

    const int & deflateWaveHoltz             = dbase.get<int>("deflateWaveHoltz");
    const int & numToDeflate                 = dbase.get<int>("numToDeflate");
    const IntegerArray & eigNumbersToDeflate = dbase.get<IntegerArray>("eigNumbersToDeflate");
    const int & onlyLoadDeflatedEigenVectors = dbase.get<int>("onlyLoadDeflatedEigenVectors");
    const int & debug                        = dbase.get<int>("debug");

    if( deflateWaveHoltz==0 || numToDeflate<=0 )
        return 0;

  // if( true )
  // {
  //   printF("### TEST: DO NOT INFLATE THE SOLUTION ####\n");
  //   return 0; 
  // }

    if( true )
        printF("### INFLATE THE SOLUTION ####\n");

    const int & numberOfFrequencies         = dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
    const RealArray & periodArray           = dbase.get<RealArray>("periodArray");  


  // -- holds intermediate results: 
    realCompositeGridFunction & u = dbase.get<realCompositeGridFunction>("uDeflate");

  // -- eigenvectors:
    realCompositeGridFunction & uev = dbase.get<realCompositeGridFunction>("uev");

  // WaveHoltz solution:
    realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");

    Integrate & integrate = dbase.get<Integrate>("integrate");

    RealArray & eig = dbase.get<RealArray>("eig");
    int numberOfEigenvectors = uev.getComponentBound(0) - uev.getComponentBase(0) + 1;


  // int freq=0;      // adjust this frequency **FIX ME**
    
  // Forcing: 
  // realCompositeGridFunction f(cg); 
    realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f"); // changed Nov 29, 2024 -- re-eval in case forcing was deflated *FIX ME**
    getHelmholtzForcing( f );

    
    Index I1,I2,I3;
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
        for( int ied=0; ied<numToDeflate; ied++ )
        {

            const int eigNumber = onlyLoadDeflatedEigenVectors ? ied : eigNumbersToDeflate(ied); // deflate this eigenvector 

            const int eigenValueNumber = eigNumbersToDeflate(ied);

      // --- Compute the inner product (phi,f)_h ----
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
                OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);
                bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
                if( ok )
                {
                    if( useOpt )
                    {
                        Real *pu         =   uLocal.getDataPointer();
                        const Real *pf   =   fLocal.getDataPointer();
                        const Real *puev = uevLocal.getDataPointer();
                        const int i1Base=uLocal.getBase(0), i2Base=uLocal.getBase(1), i3Base=uLocal.getBase(2);
                        const int nd1= uLocal.getLength(0), nd2= uLocal.getLength(1), nd3= uLocal.getLength(2);
                        #define ua(i1,i2,i3)       pu[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base          ))]
                        #define fa(i1,i2,i3,n)     pf[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base + nd3*(n)))]
                        #define ueva(i1,i2,i3,n) puev[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base + nd3*(n)))]
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            ua(i1,i2,i3) = ueva(i1,i2,i3,eigNumber)*fa(i1,i2,i3,freq);
                        }            
                          
                    }
                    else
                    {
                        uLocal(I1,I2,I3) = uevLocal(I1,I2,I3,eigNumber)*fLocal(I1,I2,I3,freq);
                    }
                }
            }

            Real phiDotF = integrate.volumeIntegral(u);  // compute discrete inner product 

            Real lambda = eig(0,eigenValueNumber);
            Real omega = frequencyArray(freq);      // this should be unadjusted value --- *** CHECK ME ***

            const Real factor = phiDotF/( omega*omega - lambda*lambda );

            if( debug>1 )
                printF("inflateSolution: freq=%d: inflate eigenvector %d, lambda=%12.4e, omega=%12.4e\n",freq,eigNumber,lambda,omega);

            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                getIndex(mg.dimension(),I1,I2,I3);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);
                OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
                bool ok=ParallelUtility::getLocalArrayBounds(uev[grid],uevLocal,I1,I2,I3);
                if( ok )
                {
                    if( useOpt )
                    {
                        Real *pv         =   vLocal.getDataPointer();
                        const Real *puev = uevLocal.getDataPointer();
                        const int i1Base=vLocal.getBase(0), i2Base=vLocal.getBase(1), i3Base=vLocal.getBase(2);
                        const int nd1= vLocal.getLength(0), nd2= vLocal.getLength(1), nd3= vLocal.getLength(2);
                        #define va(i1,i2,i3,n)     pv[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base + nd3*(n)))]
            // #define ueva(i1,i2,i3,n) puev[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base + nd3*(n)))]
                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {
                            va(i1,i2,i3,freq) += factor*ueva(i1,i2,i3,eigNumber);
                        }            
                          
                    }
                    else
                    {
                        vLocal(I1,I2,I3,freq) += factor*uevLocal(I1,I2,I3,eigNumber);
                    }
                }
            }
        }
    }
    #undef ua
    #undef va
    #undef fa
    #undef ueva
    timing(timeForDeflation) += getCPU()-cpuStart;
    return 0;
}
