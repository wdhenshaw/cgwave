// This file automatically generated from solveAugmentedKrylov.bC with bpp.
#include "CgWaveHoltz.h" 
#include "CgWave.h" 
#include "Overture.h"
//#include "gridFunctionNorms.h"
#include "display.h"
#include "ParallelUtility.h"
//#include "SparseRep.h" 
//#include "CompositeGridOperators.h"    


//#include "PlotStuff.h"
//#include "GL_GraphicsInterface.h"

#include "AugmentedKrylov.h"


// -- global variables -- do this for now 
static RealArray *pb=NULL;

static int iteration=0;
static int computeRightHandSide =-2; // = 0;
static int computeResidual      =-3; // =-1;

static CgWaveHoltz *pCgWaveHoltz; // pointer to the CgWaveHoltz solver 
  

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)       for( i3=n3a; i3<=n3b; i3++ )                        for( i2=n2a; i2<=n2b; i2++ )                      for( i1=n1a; i1<=n1b; i1++ )                    for( n=0; n<numberOfComponents; n++ )


// --------------------------------------------------------------------------------------
//   Macro: return the index's for possible active points
//            
//  NOTE: This macro appears in solveSLEPc.bC and eigenModes.bC 
// --------------------------------------------------------------------------------------

// =========================================================================================================
//     MATRIX-VECTOR MULTIPLY FOR MATRIX FREE AUGMENTED GMRES
// 
// Compute y = M*x 
// =========================================================================================================
static void matVectFunction( const RealArray & x, RealArray & y )  
{
  // PetscErrorCode ierr = 0;

    if( true || debug>2 )
    {
        if( iteration==computeRightHandSide )
            printF(" **************** MatVec for WaveHoltz (Augmented Krylov) iteration=%i : Computing the RHS *************\n",iteration);
        else
            printF(" **************** MatVec for WaveHoltz (Augmented Krylov) iteration=%i *************\n",iteration);
    }

    assert( pCgWaveHoltz!=NULL );
    CgWaveHoltz & cgWaveHoltz = *pCgWaveHoltz;
    
    const int & monitorResiduals      = cgWaveHoltz.dbase.get<int>("monitorResiduals");      // montior the residuals at every step
    const Real & numberOfActivePoints = cgWaveHoltz.dbase.get<Real>("numberOfActivePoints");
  

  // here is the CgWave solver for the time dependent wave equation
    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
    const int & numCompWaveHoltz          = cgWave.dbase.get<int>("numCompWaveHoltz");
    const int & numberOfFrequencies       = cgWave.dbase.get<int>("numberOfFrequencies");
    const int & filterTimeDerivative      = cgWave.dbase.get<int>("filterTimeDerivative");
    const int & orderOfAccuracy           = cgWave.dbase.get<int>("orderOfAccuracy");
    const int numGhost = orderOfAccuracy/2;  

  // -- cgWave solution is stored here: 
    realCompositeGridFunction & v    = cgWave.dbase.get<realCompositeGridFunction>("v");

    realCompositeGridFunction & vOld = cgWave.dbase.get<realCompositeGridFunction>("vOld");
    CompositeGrid & cg = cgWaveHoltz.cg;



  //   ---- vectorToGridFunction ----
  // --- Set v = x ---
    v=0.; // is this needed ? 

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
  // int iab[2];

    int i=0;
    for( int freq=0; freq<numCompWaveHoltz; freq++ )
    {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
      // const IntegerArray & gid = mg.gridIndexRange();

            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);

            {
                const IntegerArray & gid = mg.gridIndexRange();
                int iab[2]; // for getActivePoints 
                Iv[2]=Range(0,0);
                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        int is = 1-2*side;
                        iab[side]=gid(side,axis);
                        const int bc = mg.boundaryCondition(side,axis);
                        if( filterTimeDerivative )
                        {
              // complex valued solution: include all points : Jan 26, 2025
                            iab[side] -= is*numGhost;
                        }      
                        else if( bc==CgWave::dirichlet )
                        {
                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                        }
                        else if( bc==CgWave::neumann )
                        {
              // include boundary
                        }
                        else if( bc>0 )
                        {
                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                            OV_ABORT("error");
                        }
                        else if( bc<0 )
                        {
              // periodic -- include left end
                            if( side==1 )
                                iab[side] += is; 
                        }
                        else
                        {
              // interpolation boundary : include end 
                        }
                    }
                    Iv[axis] = Range(iab[0],iab[1]);
                }
            }

            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( maskLocal(i1,i2,i3) > 0 )
                {
                    vLocal(i1,i2,i3,freq)=x(i);

                    i++;
                }
            }
        }
    }
    if( i !=numberOfActivePoints )
    {
        printF("matVectFunction:ERROR:(I) i=%d is not equal to numberOfActivePoints=%g, numCompWaveHoltz=%d filterTimeDerivative=%d\n",
            i,numberOfActivePoints,numCompWaveHoltz,filterTimeDerivative );
        OV_ABORT("ERROR");
    }


  // *** APPLY BOUNDARY CONDITIONS to v  ****
    Real t=0.;    // is this correct ?
    if( false )
    {
        cgWave.applyBoundaryConditions( v,v, t );
    }
    else
    {
        const bool applyExplicitBoundaryConditions=true;
        const bool fillImplicitBoundaryConditions=false;
        cgWave.applyBoundaryConditions( v,v, t, applyExplicitBoundaryConditions,fillImplicitBoundaryConditions );  
    }


    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
        OV_GET_SERIAL_ARRAY(Real,vOld[grid],vOldLocal);
        vOldLocal = vLocal;  // save current iteration
    }

  // ----------------------------------------------------------------------
  // -- advance the wave equation for one period (or multiple periods ) ---
  // ----------------------------------------------------------------------
    int it = iteration+1;
    if( iteration==computeRightHandSide || iteration==-1 )
        it=0; // it =0 means this is the start of a new WaveHoltz solve so deflate the RHS

    cgWave.advance( it);


  // To apply the oprator we do this:
  //  A v^n = v^n - v^{n+1} + b

  // Set A v^n = v^{n+1}
  // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  //   vOld[grid] = v[grid]; 

  // --- gridFunctionToVector ---
  // set y = v 

    RealArray & b = *pb; 

    i=0;
    for( int freq=0; freq<numCompWaveHoltz; freq++ )
    {       
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
      // const IntegerArray & gid = mg.gridIndexRange();      

            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
            OV_GET_SERIAL_ARRAY(Real,vOld[grid],vOldLocal);
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);

            {
                const IntegerArray & gid = mg.gridIndexRange();
                int iab[2]; // for getActivePoints 
                Iv[2]=Range(0,0);
                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        int is = 1-2*side;
                        iab[side]=gid(side,axis);
                        const int bc = mg.boundaryCondition(side,axis);
                        if( filterTimeDerivative )
                        {
              // complex valued solution: include all points : Jan 26, 2025
                            iab[side] -= is*numGhost;
                        }      
                        else if( bc==CgWave::dirichlet )
                        {
                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                        }
                        else if( bc==CgWave::neumann )
                        {
              // include boundary
                        }
                        else if( bc>0 )
                        {
                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                            OV_ABORT("error");
                        }
                        else if( bc<0 )
                        {
              // periodic -- include left end
                            if( side==1 )
                                iab[side] += is; 
                        }
                        else
                        {
              // interpolation boundary : include end 
                        }
                    }
                    Iv[axis] = Range(iab[0],iab[1]);
                }
            }
            if( iteration==computeRightHandSide )
            {
        // The RHS "b" is just the solution with zero initial conditions
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
                    if( maskLocal(i1,i2,i3) > 0 )
                    {
                        y(i) = vLocal(i1,i2,i3,freq);    

                        i++;
                    }
                }
            }
            else
            {
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
                    if( maskLocal(i1,i2,i3) > 0 )
                    {
                        y(i) = vOldLocal(i1,i2,i3,freq) - vLocal(i1,i2,i3,freq) + b(i);    // y = v^n - v^{n+1} + b

                        i++;
                    }
                }        
            }
        }
    }
    if( i !=numberOfActivePoints )
    {
        printF("matVectFunction:ERROR:(end) i=%d is not equal to numberOfActivePoints=%g\n",i,numberOfActivePoints);
        OV_ABORT("ERROR");
    }  

    iteration++;
    
  // return ierr;
}


// // ============================================================
// // Copy a component of a grid function to another
// //      u[iu] <- v[iv]
// // ============================================================
// #beginMacro copyGridFunction( u,iu, v,ju )
//   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
//   {
//     OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
//     OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);  // temp space 
//     getIndex(cg[grid].dimension(),I1,I2,I3);  

//     uLocal(I1,I2,I3,iu)= vLocal(I1,I2,I3,ju);
//   }
// #endMacro



// ============================================================================
/// \brief Solve for the Helmholtz solution using WaveHoltz and 
///           Augmented KRYLOV
// ============================================================================
int CgWaveHoltz::
solveAugmentedKrylov(int argc,char **args)
{

    printF("\n @@@@@@@@@@@@@@@@ ENTERING solveAugmentedKrylov @@@@@@@@@@@@@@@@@@@@\n\n");
    Real cpuStart=getCPU();


    pCgWaveHoltz=this;

    const int myid=max(0,Communication_Manager::My_Process_Number);
    
    const Real & omega                    = dbase.get<Real>("omega");
    Real & Tperiod                        = dbase.get<Real>("Tperiod");
    const int & numPeriods                = dbase.get<int>("numPeriods");
    const int & adjustOmega               = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t   
    const int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");
  // const int & useVariableTolerance      = dbase.get<int>("useVariableTolerance"); // Vary implicit solver tolerance based on current WaveHoltz residual

    const aString & krylovType            = dbase.get<aString>("krylovType");  // gmres, bicgstab
    const int & gmresRestartLength        = dbase.get<int>("gmresRestartLength"); // restart length for WaveHoltz + GMRES  

    Real & numberOfActivePoints           = dbase.get<Real>("numberOfActivePoints");


    
  // here is the CgWave solver for the time dependent wave equation
    CgWave & cgWave = *dbase.get<CgWave*>("cgWave");
    const int & numCompWaveHoltz                = cgWave.dbase.get<int>("numCompWaveHoltz");
    const int & filterTimeDerivative            = cgWave.dbase.get<int>("filterTimeDerivative");
    const int & augmentedVectorsAreEigenvectors = cgWave.dbase.get<int>("augmentedVectorsAreEigenvectors");  // 1 = augmented vectors are true discrete eigenvectors

    const int & orderOfAccuracy                 = cgWave.dbase.get<int>("orderOfAccuracy");
    const int numGhost = orderOfAccuracy/2;  


    if( omega!=0. )
        Tperiod=numPeriods*twoPi/omega; 
    else 
        Tperiod=1.;
  
    printF("CgWaveHoltz::solveAugmentedKrylov: setting tFinal = Tperiod*numPeriods = %9.3e (numPeriods=%d) \n",Tperiod,numPeriods);
  
  // --- set values in CgWave:  *** COULD DO BETTER ***
  
    const int & computeEigenmodes       = cgWave.dbase.get<int>("computeEigenmodes");

    cgWave.dbase.get<Real>("omega")     = omega;        // ** FIX ME **
    cgWave.dbase.get<Real>("tFinal")    = Tperiod;      // ** FIX ME **
    cgWave.dbase.get<Real>("Tperiod")   = Tperiod;      // ** FIX ME **
    cgWave.dbase.get<int>("numPeriods") = numPeriods;   // ** FIX ME **
    cgWave.dbase.get<int>("adjustOmega")= adjustOmega;  // 1 : choose omega from the symbol of D+t D-t 

  // >> These next copies should now be done in initialize()
    const int & numberOfFrequencies = dbase.get<int>("numberOfFrequencies");
    cgWave.dbase.get<int>("numberOfFrequencies") = numberOfFrequencies;

    RealArray & cgWaveFrequencyArray = cgWave.dbase.get<RealArray>("frequencyArray");
    cgWaveFrequencyArray.redim(numberOfFrequencies);
    cgWaveFrequencyArray = dbase.get<RealArray>("frequencyArray");

    RealArray & cgWavePeriodArray = cgWave.dbase.get<RealArray>("periodArray");
    cgWavePeriodArray.redim(numberOfFrequencies);
    cgWavePeriodArray = dbase.get<RealArray>("periodArray"); 

    printf("\n **solveAugGmres: omega=%g, cgWaveFrequencyArray(0)=%g\n",omega,cgWaveFrequencyArray(0) );

  // Helmholtz solution is stored here:
    realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");

    realCompositeGridFunction & vOld = cgWave.dbase.get<realCompositeGridFunction>("vOld");
    Range all;
  // const int & numCompWaveHoltz = cgWave.dbase.get<int>("numCompWaveHoltz");  
    vOld.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz);


  // Save "residuals" by iteration: 
  // resVector(it) = norm( v^{n+1} - v^n )
  // RealArray & resVector = cgWave.dbase.get<RealArray>("resVector");
  // resVector.redim(maximumNumberOfIterations+10);
  // resVector=0.;
    
    int & plotOptions = cgWave.dbase.get<int>("plotOptions");
    plotOptions= CgWave::noPlotting; // turn of plotting in cgWave  


    assert( pCgWaveHoltz!=NULL );
    CgWaveHoltz & cgWaveHoltz = *pCgWaveHoltz;
    CompositeGrid & cg = cgWaveHoltz.cg;

  // int numEquations=0;
  // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  // {
  //   getIndex(cg[grid].dimension(),I1,I2,I3);
  //   numEquations += I1.getLength()*I2.getLength()*I3.getLength();
  // }
  // numEquations *= numCompWaveHoltz;

  // printF("solveAugmentedKrylov: numEquations=%d\n",numEquations);  


  // ----- compute the number of active points ---
    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];


    numberOfActivePoints = 0.;
    for( int freq=0; freq<numCompWaveHoltz; freq++ )
    {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            {
                const IntegerArray & gid = mg.gridIndexRange();
                int iab[2]; // for getActivePoints 
                Iv[2]=Range(0,0);
                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        int is = 1-2*side;
                        iab[side]=gid(side,axis);
                        const int bc = mg.boundaryCondition(side,axis);
                        if( filterTimeDerivative )
                        {
              // complex valued solution: include all points : Jan 26, 2025
                            iab[side] -= is*numGhost;
                        }      
                        else if( bc==CgWave::dirichlet )
                        {
                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                        }
                        else if( bc==CgWave::neumann )
                        {
              // include boundary
                        }
                        else if( bc>0 )
                        {
                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                            OV_ABORT("error");
                        }
                        else if( bc<0 )
                        {
              // periodic -- include left end
                            if( side==1 )
                                iab[side] += is; 
                        }
                        else
                        {
              // interpolation boundary : include end 
                        }
                    }
                    Iv[axis] = Range(iab[0],iab[1]);
                }
            }

            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( maskLocal(i1,i2,i3) > 0 )
                    numberOfActivePoints++;
            }
        }
    }
    numberOfActivePoints = ParallelUtility::getSum(numberOfActivePoints);

    printF("solveAugmentedKrylov: numberOfActivePoints=%g\n",numberOfActivePoints);

    RealArray x(numberOfActivePoints), b(numberOfActivePoints);
    pb = &b; // set global variable pb

  // --- Set initial guess to be current v ---
    RealArray x0(numberOfActivePoints);
    x0=0; 

    Real normV=0;

    int i=0;
    for( int freq=0; freq<numCompWaveHoltz; freq++ )
    {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            {
                const IntegerArray & gid = mg.gridIndexRange();
                int iab[2]; // for getActivePoints 
                Iv[2]=Range(0,0);
                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        int is = 1-2*side;
                        iab[side]=gid(side,axis);
                        const int bc = mg.boundaryCondition(side,axis);
                        if( filterTimeDerivative )
                        {
              // complex valued solution: include all points : Jan 26, 2025
                            iab[side] -= is*numGhost;
                        }      
                        else if( bc==CgWave::dirichlet )
                        {
                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                        }
                        else if( bc==CgWave::neumann )
                        {
              // include boundary
                        }
                        else if( bc>0 )
                        {
                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                            OV_ABORT("error");
                        }
                        else if( bc<0 )
                        {
              // periodic -- include left end
                            if( side==1 )
                                iab[side] += is; 
                        }
                        else
                        {
              // interpolation boundary : include end 
                        }
                    }
                    Iv[axis] = Range(iab[0],iab[1]);
                }
            }

            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( maskLocal(i1,i2,i3) > 0 )
                {
                    x0(i) = vLocal(i1,i2,i3,freq);
                    normV = max( normV,x0(i));
                    i++;
                }
            }
        }
    }
    assert( i==numberOfActivePoints );

    printF("solveAugmentedKrylov: Set initial guess to v, max-norm(v)=%8.2e, numberOfActivePoints=%g\n",normV,numberOfActivePoints);



  // ------ Compute the RHS b : solve with zero initial conditions --------
    Real bNorm;

    x=0; // initial condition for WaveHoltz when computing b 
  // b=0; 
    iteration=computeRightHandSide;

    matVectFunction( x, b ); 

    bNorm = AugmentedKrylov::norm( b ); 
    Real bNorm2h = bNorm/sqrt(numberOfActivePoints); 
    printF("solveAugmentedKrylov: RHS is b: l2-norm(b)=%9.3e, L2h-norm(b)=%9.2e\n",bNorm,bNorm2h);

    

  // ----------- Set the augmented vectors ------------
    RealArray W; // holds augmented vectors 

    const int & numToDeflate                 = cgWave.dbase.get<int>("numToDeflate");
    const int & deflateWaveHoltz             = cgWave.dbase.get<int>("deflateWaveHoltz");
    if( deflateWaveHoltz )
    {
    // **NOTE** WE must have alrady called cgWave.advance before getting here so that  cgWave.initializeDeflation() has been called
    // to read in the known eigenvectors 

        const int & onlyLoadDeflatedEigenVectors = cgWave.dbase.get<int>("onlyLoadDeflatedEigenVectors");
        const IntegerArray & eigNumbersToDeflate = cgWave.dbase.get<IntegerArray>("eigNumbersToDeflate");

        realCompositeGridFunction & uev          = cgWave.dbase.get<realCompositeGridFunction>("uev");
        W.redim(numberOfActivePoints,numToDeflate);

        for( int ied=0; ied<numToDeflate; ied++ )
        {
            const int eigNumber = onlyLoadDeflatedEigenVectors ? ied : eigNumbersToDeflate(ied); // deflate this eigenvector
            int i=0;
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);

                {
                    const IntegerArray & gid = mg.gridIndexRange();
                    int iab[2]; // for getActivePoints 
                    Iv[2]=Range(0,0);
                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                    {
                        for( int side=0; side<=1; side++ )
                        {
                            int is = 1-2*side;
                            iab[side]=gid(side,axis);
                            const int bc = mg.boundaryCondition(side,axis);
                            if( filterTimeDerivative )
                            {
                // complex valued solution: include all points : Jan 26, 2025
                                iab[side] -= is*numGhost;
                            }      
                            else if( bc==CgWave::dirichlet )
                            {
                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                            }
                            else if( bc==CgWave::neumann )
                            {
                // include boundary
                            }
                            else if( bc>0 )
                            {
                                printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                                OV_ABORT("error");
                            }
                            else if( bc<0 )
                            {
                // periodic -- include left end
                                if( side==1 )
                                    iab[side] += is; 
                            }
                            else
                            {
                // interpolation boundary : include end 
                            }
                        }
                        Iv[axis] = Range(iab[0],iab[1]);
                    }
                }
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
                    if( maskLocal(i1,i2,i3) > 0 )
                    {
                        W(i,ied) = uevLocal(i1,i2,i3,eigNumber);
                        i++;
                    }
                }
            }

        }
        assert( i==numberOfActivePoints );

    }


  // -----------------------------------------------------------------------
  //              Augmented GMRES solve
  // -----------------------------------------------------------------------
    const Real & tol = dbase.get<Real>("tol");
    iteration=0;    // Start true iterations

    AugmentedKrylov auKrylov;

  // Choose the type of Augmented Krylov Solver 
    AugmentedKrylov::KrylovTypesEnum augmentedKrylovType = AugmentedKrylov::gmres;
    if( krylovType=="gmres" ) 
    {
          augmentedKrylovType = AugmentedKrylov::gmres;
          printF("CgWaveHoltz::solveAugmentedKrylov: setting AUGMENTED Krylov solver GMRES.\n");
    }
    else if( krylovType=="bicgstab" || krylovType=="bcgs" || krylovType=="bicgs" )
    {
        augmentedKrylovType = AugmentedKrylov::biConjugateGradientStabilized;
        printF("CgWaveHoltz::solveAugmentedKrylov: setting AUGMENTED Krylov solver to bi-conjugate-gradient stabilized.\n");
    }
    else if( krylovType=="cg" )
    {
        augmentedKrylovType = AugmentedKrylov::conjugateGradient;
        printF("CgWaveHoltz::solveAugmentedKrylov: setting AUGMENTED Krylov solver to conjugate-gradient.\n");
    }  
    else
    {
        printf("CgWaveHoltz::solveAugmentedKrylov:WARNING: unknown krylovType=[%s] for augmented methods.\n"
                      "  Valid options: gmres, bicgstab, cg\n",(const char*)krylovType);
        OV_ABORT("error");
    }    

    auKrylov.setKrylovType( augmentedKrylovType );

    const int eigenvectorsAreTrueEigenvectors = cgWave.dbase.get<int>("eigenvectorsAreTrueEigenvectors");
    printF("\n $$$$$$$$$$$$$ AuKrylov:eigenvectorsAreTrueEigenvectors=%d, augmentedVectorsAreEigenvectors=%d $$$$$$$$ (both must be true to assume this) \n\n",
        eigenvectorsAreTrueEigenvectors,augmentedVectorsAreEigenvectors);
    if( eigenvectorsAreTrueEigenvectors && augmentedVectorsAreEigenvectors)
    {
    // Supply beta eigevalues for the augmented eigenvectors -- these are used to avoid multiplying the augmented vectors by "A"
        const RealArray & betaDeflate = cgWave.dbase.get<RealArray>("betaDeflate");
        auKrylov.setAugmentedEigenvalues( betaDeflate );

    }

    auKrylov.solve( matVectFunction, b, x0, W, maximumNumberOfIterations, tol, x );


  // save some timings so we can estimate the cost by pre-computing somethings
    if( auKrylov.dbase.has_key("cpuAugmentedInitialMatVects") )
    {
        if( !dbase.has_key("cpuAugmentedInitialMatVects") )
            dbase.put<Real>("cpuAugmentedInitialMatVects")=0.;
        dbase.get<Real>("cpuAugmentedInitialMatVects")=auKrylov.dbase.get<Real>("cpuAugmentedInitialMatVects");
    }
    if( auKrylov.dbase.has_key("cpuAugmentedQR") )
    {
        if( !dbase.has_key("cpuAugmentedQR") )
            dbase.put<Real>("cpuAugmentedQR");
        dbase.get<Real>("cpuAugmentedQR")=auKrylov.dbase.get<Real>("cpuAugmentedQR"); 
    }   

  // ---- Form the solution x as a grid function v -----
    i=0;
    for( int freq=0; freq<numCompWaveHoltz; freq++ )
    {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);

            {
                const IntegerArray & gid = mg.gridIndexRange();
                int iab[2]; // for getActivePoints 
                Iv[2]=Range(0,0);
                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        int is = 1-2*side;
                        iab[side]=gid(side,axis);
                        const int bc = mg.boundaryCondition(side,axis);
                        if( filterTimeDerivative )
                        {
              // complex valued solution: include all points : Jan 26, 2025
                            iab[side] -= is*numGhost;
                        }      
                        else if( bc==CgWave::dirichlet )
                        {
                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                        }
                        else if( bc==CgWave::neumann )
                        {
              // include boundary
                        }
                        else if( bc>0 )
                        {
                            printF("getActivePointIndex:ERROR: unknown bc=%d for grid=%d\n",bc,grid);
                            OV_ABORT("error");
                        }
                        else if( bc<0 )
                        {
              // periodic -- include left end
                            if( side==1 )
                                iab[side] += is; 
                        }
                        else
                        {
              // interpolation boundary : include end 
                        }
                    }
                    Iv[axis] = Range(iab[0],iab[1]);
                }
            }
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( maskLocal(i1,i2,i3) > 0 )
                {
                    vLocal(i1,i2,i3,freq) = x(i);
                    i++;
                }
            }
        }
    }
    assert( i==numberOfActivePoints );

  // *** APPLY BOUNDARY CONDITIONS to v  ****
    Real t=0.;  // is this correct ?
    const bool applyExplicitBoundaryConditions=true;
    const bool fillImplicitBoundaryConditions=false;
    cgWave.applyBoundaryConditions( v,v, t, applyExplicitBoundaryConditions,fillImplicitBoundaryConditions );
    if( true )
        cgWave.applyBoundaryConditions( v,v, t, applyExplicitBoundaryConditions,fillImplicitBoundaryConditions );


    int & numberOfIterations = dbase.get<int>("numberOfIterations");
    numberOfIterations = auKrylov.getNumberOfIterations();
    printF("number of AugmentedKrylov iterations=%d\n",numberOfIterations);

  // save num mat-vects: 
    int & numberOfMatrixVectorMultiplications = dbase.get<int>("numberOfMatrixVectorMultiplications");
    numberOfMatrixVectorMultiplications = auKrylov.getNumberOfMatrixVectorProducts();

    RealArray & resVector = cgWave.dbase.get<RealArray>("resVector");
    resVector.redim(numberOfIterations);
    Range Ni = numberOfIterations;

  // auKrylov return the L2h norm of the residual = || r ||_2 / sqrt(numberOfPoints)
    resVector(Ni) = auKrylov.getResidualVector()(Ni);

  // // auKrylov returns || residual ||/|| b \||
  // // We save the L2h norm of the residual
  // resVector(Ni) *= bNorm/sqrt(numberOfActivePoints); 

    ::display(resVector, "resVector","%9.2e ");

  // --- Compute the average convergence rate ----
    Real & convergenceRate          = dbase.get<Real>("convergenceRate");
    Real & convergenceRatePerPeriod = dbase.get<Real>("convergenceRatePerPeriod");
    
    if( numberOfIterations>0 )
    {
        convergenceRate          = pow( resVector(numberOfIterations-1)/resVector(0), 1./( numberOfIterations ) ); 
        convergenceRatePerPeriod = pow( resVector(numberOfIterations-1)/resVector(0), 1./( numberOfIterations*numPeriods) ); 
    }
    else
    {
    // No iterations we used -- the current solution must have met the tolerance
        convergenceRate          = 1.; 
        convergenceRatePerPeriod = 1.;
    }

    dbase.get<real>("maxResidual") = auKrylov.getResidual(); // really scaled L2 residual

    Real cpuTotal = getCPU() - cpuStart;
    printF("solveAugmentedKrylov: cpuTotal=%9.2e(s)\n",cpuTotal);
    return 0;
}




