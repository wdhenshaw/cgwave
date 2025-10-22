// This file automatically generated from solvePETSc.bC with bpp.
#include "CgWaveHoltz.h" 
#include "CgWave.h" 
#include "Overture.h"
#include "gridFunctionNorms.h"
#include "display.h"
#include "ParallelUtility.h"
#include "SparseRep.h" 

#include "PlotStuff.h"
#include "GL_GraphicsInterface.h"

// krb do not use extern "C" if PETSc is linked using BOPT=g_c++
// extern "C"
// {
// // *wdh* 2015/09/31  To avoid having PETSc include complex.h do this: 
// #include "petscconf.h"
// #undef PETSC_HAVE_CXX_COMPLEX
// #include "petscksp.h"
// }

// for petsc 3.18.2
#include "petscksp.h"

// SLEPc 
#ifdef CGWAVE_USE_SLEPC
#include <slepceps.h>
#endif

static char help[] = "CgWaveHoltz test of PETSc\n";

static bool useMatrixUtilities=true; // true = use new matrix utilities (needed for parallel)

// testCase : 
// 0 = test solving laplace equation
// 1 = Real run
static int testCase =1; 
static int iteration=0;

static CgWaveHoltz *pCgWaveHoltz; // pointer to the CgWaveHoltz solver 
  

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)       for( i3=n3a; i3<=n3b; i3++ )                        for( i2=n2a; i2<=n2b; i2++ )                      for( i1=n1a; i1<=n1b; i1++ )                    for( n=0; n<numberOfComponents; n++ )


// -- global variables -- do this for now 
static KSP  ksp=NULL;     /* linear solver context */
static Mat *pA=NULL;
static Vec *pb=NULL;
static int mA,nA;

static int computeRightHandSide =-2; // = 0;
static int computeResidual      =-3; // =-1;
static int assignSolution       =-4;

// *wdh* Jan 5, 2022 --> changed cgWave to zero out unused points => thus we do not need to check the mask
static bool checkMask = false;  // if true, do not use values of the solution where mask==0 

// -- see eig/src/fillInterpolationCoefficients.bC
// #define nab(side,axis,p,grid) pnab[(side)+2*( (axis) + 3*( (p) + numberOfProcessors*( (grid) ) ) )]
// #define ndab(axis,p,grid) (nab(1,axis,p,grid)-nab(0,axis,p,grid)+1)
// #define noffset(p,grid) pnoffset[(p)+numberOfProcessors*(grid)]


// int Ogev::
// getGlobalIndex( int n, int *iv, int grid, int p ) 
// // ===============================================================================
// /// \brief: Return the global index (equation number) given the point, grid and processor
// /// \note: These equation numbers are base=0.
// // ===============================================================================
// {
//   // printF("getGlobalIndex: n=%d, (i1,i2,i3)=(%d,%d,%d) grid=%d, p=%d\n",n,iv[0],iv[1],iv[2],grid,p);

//   return  n + numberOfComponents*(
//          (iv[0]-nab(0,axis1,p,grid))+ndab(0,p,grid)*(
//           iv[1]-nab(0,axis2,p,grid) +ndab(1,p,grid)*(
//           iv[2]-nab(0,axis3,p,grid))) + noffset(p,grid) );
// }

// ======================================================================================
/// \brief Initialize PETSC if it hasn't already been initialized
// ======================================================================================
int CgWave::initializePETSc( int argc /* = 0 */, char **args  /* =NULL */ )
{
    // Do this in case we use PETSc  *wdh* Jan 7, 2023
    int & petscIsInitialized = dbase.get<int>("petscIsInitialized");
    if( !petscIsInitialized )
    {
        printF("\n######## INIT PETSc ########\n");
        static char help[] = "CgWaveHoltz test of PETSc\n";
        petscIsInitialized=true;
        PetscInitialize(&argc,&args,(char *)0,help);
    } 

    return 0; 
}

// --------------------------------------------------------------------------------------
//   Macro: return the index's for possible active points
//            
//  NOTE: This macro appears in solveSLEPc.bC and eigenModes.bC 
// --------------------------------------------------------------------------------------

// =========================================================================================================
//     MATRIX-VECTOR MULTIPLY FOR MATRIX FREE KRYLOV SOLVERS
// 
// Compute y = M*x 
//
//  NEW VERSION USING ACTIVE POINTS ONLY
//
// =========================================================================================================
extern PetscErrorCode waveHoltzMatrixVectorMultiply(Mat m ,Vec x, Vec y)
{
    PetscErrorCode ierr = 0;

    if( true )
    {
        if( iteration==computeRightHandSide )
            printF("\n ++++++++ WaveHoltz Matrix vector multiply routine: EVALUATE THE RHS\n");
        else if( iteration==computeResidual )
            printF("\n ++++++++ WaveHoltz Matrix vector multiply routine: CHECK RESIDUAL\n");
        else if( iteration==assignSolution )
            printF("\n ++++++++ WaveHoltz Matrix vector multiply routine: ASSIGN THE SOLUTION\n");    
        else
            printF("\n ++++++++ WaveHoltz Matrix vector multiply routine: waveHoltzMatrixVectorMultiply (activePoints) called iteration=%d\n",iteration);
    }


  // ---- Helmholtz solver ----

  // printF(" **************** MatVec for PETSc iteration=%i (-2=rhs, -3=residual) *************\n\n",iteration);
    if( iteration==computeRightHandSide )
    {
        printF(" >>>> Call cgWave to COMPUTE the Right-hand-side to Ax=b : b = Pi * v(0) \n");    
    }

    assert( pCgWaveHoltz!=NULL );
    CgWaveHoltz & cgWaveHoltz = *pCgWaveHoltz;
    
    const int & monitorResiduals      = cgWaveHoltz.dbase.get<int>("monitorResiduals");      // montior the residuals at every step
    const Real & numberOfActivePoints = cgWaveHoltz.dbase.get<Real>("numberOfActivePoints");
    const int & useVariableTolerance  = cgWaveHoltz.dbase.get<int>("useVariableTolerance"); // Vary implicit solver tolerance based on current WaveHoltz residual

  // here is the CgWave solver for the time dependent wave equation
    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");

    const int & numCompWaveHoltz          = cgWave.dbase.get<int>("numCompWaveHoltz");
    const int & upwind                    = cgWave.dbase.get<int>("upwind");
    const int & filterTimeDerivative      = cgWave.dbase.get<int>("filterTimeDerivative");
    const int & orderOfAccuracy           = cgWave.dbase.get<int>("orderOfAccuracy");
    const int numGhost = orderOfAccuracy/2;

  // const int & numberOfFrequencies = cgWave.dbase.get<int>("numberOfFrequencies");

  // -- cgWave solution is stored here: 
    realCompositeGridFunction & v    = cgWave.dbase.get<realCompositeGridFunction>("v");

    realCompositeGridFunction & vOld = cgWave.dbase.get<realCompositeGridFunction>("vOld");
    CompositeGrid & cg = cgWaveHoltz.cg;

  // move to solvePETSc for parallel
  // if( !cgWaveHoltz.dbase.has_key("bcg") )
  // {
  //   realCompositeGridFunction & bcg = cgWaveHoltz.dbase.put<realCompositeGridFunction>("bcg");
  //   Range all;
  //   bcg.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz);
  // }
    realCompositeGridFunction & bcg = cgWaveHoltz.dbase.get<realCompositeGridFunction>("bcg");    


    PetscScalar *xl, *yl;
    VecGetArray(x,&xl);  // get the local array from Petsc
    VecGetArray(y,&yl);  // get the local array from Petsc
    int iStart,iEnd;
    ierr = VecGetOwnershipRange(x,&iStart,&iEnd);CHKERRQ(ierr);
  // printf("  [iStart,iEnd]=[%i,%i]\n",iStart,iEnd);


  // --- Set v = x ---
    v=0.; // is this needed ? 

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    int iab[2];

    if( useMatrixUtilities )
    {
        cgWave.vectorToGridFunction( xl, v, iStart,iEnd );
    }
    else
    {
        int i=0;
        for( int freq=0; freq<numCompWaveHoltz; freq++ )
        {
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const IntegerArray & gid = mg.gridIndexRange();
                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

        // getIndex(cg[grid].dimension(),I1,I2,I3);
                {
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
                            else if( upwind || bc==CgWave::neumann )
                            {
                // include ghost ??
                // iab[side] -= is;
                            }
                            else if( bc==CgWave::abcEM2  || bc==CgWave::absorbing || bc==CgWave::radiation )
                            {
                // include ghost 
                                iab[side] -= is*numGhost;
                            }      
                            else if( bc==CgWave::dirichlet )
                            {
                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                 // iab[side] -= is*numGhost; // ************************************** TEMP ***********
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

                OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
                    assert( i<iEnd );

                    vLocal(i1,i2,i3,freq)=xl[i]; // new way  Jan 5, 2022

          // // What about upwinding -- it may use unused points
          // // if( !checkMask || maskLocal(i1,i2,i3)!=0 )  // this may not be needed if xl[i] is always zero
          // if( true )  // Jan 5, 2022
          //   vLocal(i1,i2,i3,freq)=xl[i];
          // else
          //   vLocal(i1,i2,i3,freq)=0.; 
                    i++;
                }
            }
        }
        assert( i==iEnd );
    }

  // We are just assigning "x" into "v"
    if( iteration==assignSolution )
    {    
        const bool applyExplicitBoundaryConditions=true;
        const bool fillImplicitBoundaryConditions=false;
        Real t = 0.;
        cgWave.applyBoundaryConditions( v,v, t, applyExplicitBoundaryConditions,fillImplicitBoundaryConditions ); 

        const int & numberOfFrequencies = cgWave.dbase.get<int>("numberOfFrequencies");
        if( numberOfFrequencies>1 ) 
        {
            Index I1,I2,I3;
            for( int freq=1; freq<numberOfFrequencies; freq++ )
            {
        /// In order to use applyBoundaryConditions we need to copy to a temporary
        // save component "freq" into component vOld(I1,I2,I2,0)
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                    OV_GET_SERIAL_ARRAY(Real,vOld[grid],vOldLocal);

                    MappedGrid & mg = cg[grid];
                    getIndex(mg.dimension(),I1,I2,I3);
                    const int includeParallelGhost=1;
                    bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeParallelGhost); 
                    if( ok )
                        vOldLocal(I1,I2,I3,0) = vLocal(I1,I2,I3,freq);  // save component "freq" into component vOld(I1,I2,I2,0)
                } 

        // --- apply BC's to component freq ---
        // printF("XXXXXXXXX solvePETSc : apply BCs to freq=%d after solve is complete XXXXXXXXX\n",freq);

                cgWave.applyBoundaryConditions( vOld,vOld, t, applyExplicitBoundaryConditions,fillImplicitBoundaryConditions );

        // -- copy back ---
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                    OV_GET_SERIAL_ARRAY(Real,vOld[grid],vOldLocal);

                    MappedGrid & mg = cg[grid];
                    getIndex(mg.dimension(),I1,I2,I3);
                    const int includeParallelGhost=1;
                    bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeParallelGhost); 
                    if( ok )          
                        vLocal(I1,I2,I3,freq) = vOldLocal(I1,I2,I3,0);  
                } 


            }
        }

        return 0;
    }


    if( false && iteration==1 )
    {
        ::display(v[0],"v (iteration 1)","%5.2f ");
        
    }
    
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
        OV_GET_SERIAL_ARRAY(Real,vOld[grid],vOldLocal);
        vOldLocal = vLocal;  // save current guess

    // vOld[grid] = v[grid];  // save current guess
    }

  // *** APPLY BOUNDARY CONDITIONS to v 
  // cgWave.applyBoundaryConditions( realCompositeGridFunction & u, realCompositeGridFunction & un, real t,
  //                        bool applyExplicitBoundaryConditions /* = false */,
  //                        bool fillImplicitBoundaryConditions  /* = false */ )

  // ----------------------------------------------------------------------
  // -- advance the wave equation for one period (or multiple periods ) ---
  // ----------------------------------------------------------------------
  // printF("\n\n###### solvePETSc - call advance, iteration=%d, initializeTimeStepping=%d\n\n",iteration, initializeTimeStepping);
    int it = iteration+1;
    if( iteration==computeRightHandSide || iteration==-1 )
        it=0; // it =0 means this is the start of a new WaveHoltz solve so deflate the RHS
    else if( iteration==computeResidual )
        it=1; 
    cgWave.advance( it );

  // cgWave.advance( iteration );
    

    if( iteration==computeRightHandSide )
    {
    // ---------- INITIALIZE ITERATIONS -----

        printF("COMPUTE b = Pi * v(0) \n");

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,bcg[grid],bcgLocal);

            bcgLocal = vLocal; 
        }

        if( false )
        {
            ::display(bcg[0],"b: iteration 0 ","%5.2f ");
        }

    //compute y = Pi * x 
    // set y = v 
        if( useMatrixUtilities )
        {
            cgWave.gridFunctionToVector( v,yl, iStart,iEnd );
        }
        else
        {
            int i=0;
            for( int freq=0; freq<numCompWaveHoltz; freq++ )
            {      
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
          // getIndex(cg[grid].dimension(),I1,I2,I3);
                    MappedGrid & mg = cg[grid];
                    const IntegerArray & gid = mg.gridIndexRange();
                    {
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
                                else if( upwind || bc==CgWave::neumann )
                                {
                  // include ghost ??
                  // iab[side] -= is;
                                }
                                else if( bc==CgWave::abcEM2  || bc==CgWave::absorbing || bc==CgWave::radiation )
                                {
                  // include ghost 
                                    iab[side] -= is*numGhost;
                                }      
                                else if( bc==CgWave::dirichlet )
                                {
                                      iab[side] += is;  // Dirichlet BC -- ignore the boundary
                   // iab[side] -= is*numGhost; // ************************************** TEMP ***********
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

                    OV_GET_SERIAL_ARRAY(int,cg[grid].mask(),maskLocal);

                    OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        assert( i<iEnd );
                        yl[i]= vLocal(i1,i2,i3,freq); // new way Jan 5, 2022

                        i++;
                    }
                }
            }
            assert( i==iEnd );
        }
    }
    else if( iteration==computeResidual )
    {
    // ---- iteration=computeResidual: (after final step usually) check the residual in the computed solution  : v^n - v^{n+1}
    // Compute :
    //       A v^n = v^n - v^{n+1} + b 

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            vOld[grid] -= v[grid]; // 
    
        const Real & tol = cgWaveHoltz.dbase.get<Real>("tol");
        Real & maxResidual = cgWaveHoltz.dbase.get<Real>("maxResidual");

        maxResidual = maxNorm(vOld);
        printF("it=%d:  max(residual) = max(|v-vOld|)=%8.2e, tol=%g\n",iteration,maxResidual,tol);

    }
    else
    {
    // ----------- MATRIX VECTOR MULTIPLY FOR Krylov solver -------

        const int & computeEigenmodes       = cgWave.dbase.get<int>("computeEigenmodes");

        if( !computeEigenmodes )
        {
      // Compute :
      //       A v^n = v^n - v^{n+1} + b 
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                vOld[grid] -= v[grid]; // 


            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                vOld[grid] += bcg[grid];

        }
        else
        {
      // -- for eigenmodes -- assume b=0
      // Set A v^n = v^{n+1}
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                vOld[grid] = v[grid]; // 
        }

    // should we interpolate vOld?
        if( false ) //  Jan 5, 2022 : change to false : FIXES SIC
        {
            vOld.interpolate();
        }

        assert( pb !=NULL );
        Vec & b = *pb;
        PetscScalar *bl;
        VecGetArray(b,&bl);  // get the local array from Petsc

    // ---- set y = vOld ---- 
        if( useMatrixUtilities )
        {
            cgWave.gridFunctionToVector( vOld,yl, iStart,iEnd );
        }
        else 
        {   
            int i=0;
            for( int freq=0; freq<numCompWaveHoltz; freq++ )
            {       
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
          //  getIndex(cg[grid].dimension(),I1,I2,I3);
                    MappedGrid & mg = cg[grid];
                    const IntegerArray & gid = mg.gridIndexRange();        
                    {
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
                                else if( upwind || bc==CgWave::neumann )
                                {
                  // include ghost ??
                  // iab[side] -= is;
                                }
                                else if( bc==CgWave::abcEM2  || bc==CgWave::absorbing || bc==CgWave::radiation )
                                {
                  // include ghost 
                                    iab[side] -= is*numGhost;
                                }      
                                else if( bc==CgWave::dirichlet )
                                {
                                      iab[side] += is;  // Dirichlet BC -- ignore the boundary
                   // iab[side] -= is*numGhost; // ************************************** TEMP ***********
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

                    OV_GET_SERIAL_ARRAY(int,cg[grid].mask(),maskLocal);
                    OV_GET_SERIAL_ARRAY(Real,vOld[grid],vOldLocal);
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        assert( i<iEnd );
                        yl[i]= vOldLocal(i1,i2,i3,freq);    // y = A*x 
                          
                        i++;
                    }
                }
            }
            assert( i==iEnd );
        }

    }

    const int & computeEigenmodes = cgWave.dbase.get<int>("computeEigenmodes");
    if( monitorResiduals && iteration>=0 && !computeEigenmodes )
    {
    // ----- Optionally save the current residual -----
    // Save "residuals" by iteration: 
    // resVector(it) = norm( v^{n+1} - v^n )
        RealArray & resVector = cgWave.dbase.get<RealArray>("resVector");  

        if( iteration<0 )
        {
            printF("solvePETSC: ERROR: INVALID iteration=%d.\n",iteration);
            OV_ABORT("ERROR");
        }
        const int & maximumNumberOfIterations = cgWaveHoltz.dbase.get<int>("maximumNumberOfIterations");
        if( iteration>=resVector.getLength(0) ) 
        {
            printF("solvePETSC: Increasing the size of the vector that holds residuals, iteration=%d, maximumNumberOfIterations=%d\n",iteration,maximumNumberOfIterations);
            resVector.resize(resVector.getLength(0)*2); 
        }

    // resVector(iteration)= cgWaveHoltz.residual();   // this is not correct since solution is not the right one

        Real kspResidual;
        assert( ksp !=NULL );
        KSPGetResidualNorm(ksp,&kspResidual); 
        kspResidual /= sqrt(numberOfActivePoints);  // make an approximate L2h norm

        resVector(iteration)= kspResidual;
        printF("\n ##### SAVE KRYLOV RESIDUAL: iteration=%d: L2h-residual=%9.2e, ratio=%5.2f \n\n",kspResidual,resVector(iteration)/resVector(max(iteration-1,0)));
        
  
        if( useVariableTolerance ) // ** THIS DOES NOT SEEM TO WORK -- KRYLOV NEEDS ACCURATE Matrix-vector multiplies ***
        {
            const CgWave::TimeSteppingMethodEnum & timeSteppingMethod = cgWave.dbase.get<CgWave::TimeSteppingMethodEnum>("timeSteppingMethod");
            if( timeSteppingMethod == CgWave::implicitTimeStepping )
            {
                Oges & impSolver = cgWave.dbase.get<Oges>("impSolver");
                if( impSolver.isSolverIterative() )
                {
                    Real rtol= max(1e-12,min(1.e-3,kspResidual/10));
          // rtol=1e-10;

                    Real atol=rtol;

                    printF("CgWaveHoltz:waveHoltzMatrixVectorMultiply: set implicit solver tolerance: rtol=%9.2e, atol=%9.2e (based on residual=%9.3e)\n",atol,rtol,kspResidual);
            
                    impSolver.set(OgesParameters::THErelativeTolerance,rtol);        
                    impSolver.set(OgesParameters::THEabsoluteTolerance,atol);   
                } 
            }    
        }

    // // There is no residual for iteration=0 since the Krylov solver is just computing the first A*x
    // if( iteration==1 )
    //   resVector(0) = resVector(iteration); // do this for now -- ksp is not built yet for iteration=0
    } 

    iteration++;
    
    return ierr;
}



// =========================================================================================================
//     MATRIX-VECTOR MULTIPLY FOR MATRIX FREE KRYLOV SOLVERS
// 
// Compute y = M*x 
// =========================================================================================================
extern PetscErrorCode waveHoltzMatrixVectorMultiplyOld(Mat m ,Vec x, Vec y)
{
    PetscErrorCode ierr = 0;

    if( true )
    {
        if( iteration==computeRightHandSide )
            printF("\n ++++++++ WaveHoltz Matrix vector multiply routine: EVALUATE THE RHS\n");
        else if( iteration==computeResidual )
            printF("\n ++++++++ WaveHoltz Matrix vector multiply routine: CHECK RESIDUAL\n");
      else if( iteration==assignSolution )
            printF("\n ++++++++ WaveHoltz Matrix vector multiply routine: ASSIGN THE SOLUTION\n");     
        else
            printF("\n ++++++++ WaveHoltz Matrix vector multiply routine: waveHoltzMatrixVectorMultiply called iteration=%d\n",iteration);
    }


  // ---- Helmholtz solver ----

  // printF(" **************** MatVec for PETSc iteration=%i (-2=rhs, -3=residual) *************\n\n",iteration);
    if( iteration==computeRightHandSide )
    {
        printF(" >>>> Call cgWave to COMPUTE the Right-hand-side to Ax=b : b = Pi * v(0) \n");    
    }

    assert( pCgWaveHoltz!=NULL );
    CgWaveHoltz & cgWaveHoltz = *pCgWaveHoltz;
    
    const int & monitorResiduals      = cgWaveHoltz.dbase.get<int>("monitorResiduals");      // montior the residuals at every step
    const Real & numberOfActivePoints = cgWaveHoltz.dbase.get<Real>("numberOfActivePoints");
    const int & useVariableTolerance  = cgWaveHoltz.dbase.get<int>("useVariableTolerance"); // Vary implicit solver tolerance based on current WaveHoltz residual

  // here is the CgWave solver for the time dependent wave equation
    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");

    const int & numCompWaveHoltz          = cgWave.dbase.get<int>("numCompWaveHoltz");
  // const int & numberOfFrequencies = cgWave.dbase.get<int>("numberOfFrequencies");

  // -- cgWave solution is stored here: 
    realCompositeGridFunction & v    = cgWave.dbase.get<realCompositeGridFunction>("v");

    realCompositeGridFunction & vOld = cgWave.dbase.get<realCompositeGridFunction>("vOld");
    CompositeGrid & cg = cgWaveHoltz.cg;

  // if( !cgWaveHoltz.dbase.has_key("bcg") )
  // {
  //   realCompositeGridFunction & bcg = cgWaveHoltz.dbase.put<realCompositeGridFunction>("bcg");
  //   Range all;
  //   bcg.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz);
  // }
    realCompositeGridFunction & bcg = cgWaveHoltz.dbase.get<realCompositeGridFunction>("bcg");    


    PetscScalar *xl, *yl;
    VecGetArray(x,&xl);  // get the local array from Petsc
    VecGetArray(y,&yl);  // get the local array from Petsc
    int iStart,iEnd;
    ierr = VecGetOwnershipRange(x,&iStart,&iEnd);CHKERRQ(ierr);
  // printf("  [iStart,iEnd]=[%i,%i]\n",iStart,iEnd);


  // --- Set v = x ---
    v=0.; // is this needed ? 

    Index I1,I2,I3;

    if( useMatrixUtilities )
    {
        cgWave.vectorToGridFunction( xl,v, iStart,iEnd );
    }
    else
    {
        int i=0;
        for( int freq=0; freq<numCompWaveHoltz; freq++ )
        {
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

                getIndex(cg[grid].dimension(),I1,I2,I3);

                OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
                    assert( i<iEnd );

                    vLocal(i1,i2,i3,freq)=xl[i]; // new way  Jan 5, 2022

          // // What about upwinding -- it may use unused points
          // // if( !checkMask || maskLocal(i1,i2,i3)!=0 )  // this may not be needed if xl[i] is always zero
          // if( true )  // Jan 5, 2022
          //   vLocal(i1,i2,i3,freq)=xl[i];
          // else
          //   vLocal(i1,i2,i3,freq)=0.; 
                    i++;
                }
            }
        }
        assert( i==iEnd );
    }

  // // We are just assigning "x" into "v"
  // if( iteration==assignSolution )
  // We are just assigning "x" into "v"
    if( iteration==assignSolution )
    {    
        const bool applyExplicitBoundaryConditions=true;
        const bool fillImplicitBoundaryConditions=false;
        Real t = 0.;
        cgWave.applyBoundaryConditions( v,v, t, applyExplicitBoundaryConditions,fillImplicitBoundaryConditions ); 

        return 0;
    }    
        
    if( false && iteration==1 )
    {
        ::display(v[0],"v (iteration 1)","%5.2f ");
        
    }
    
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
        OV_GET_SERIAL_ARRAY(Real,vOld[grid],vOldLocal);
        vOldLocal = vLocal;  // save current guess

    // vOld[grid] = v[grid];  // save current guess
    }

  // ----------------------------------------------------------------------
  // -- advance the wave equation for one period (or multiple periods ) ---
  // ----------------------------------------------------------------------
    int it = iteration+1;
    if( iteration==computeRightHandSide || iteration==-1 )
        it=0; // it =0 means this is the start of a new WaveHoltz solve so deflate the RHS
    else if( iteration==computeResidual )
        it=1;  

    cgWave.advance( it );  
  // cgWave.advance( iteration );
    

    if( iteration==computeRightHandSide )
    {
    // ---------- INITIALIZE ITERATIONS -----

        printF("COMPUTE b = Pi * v(0) \n");

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(Real,bcg[grid],bcgLocal);

            bcgLocal = vLocal; 

      // else
      // {
      //   // ** Jan 5, 2020 -> testing: 
      //   OV_GET_SERIAL_ARRAY(int,cg[grid].mask(),maskLocal);
      //   Range Nf = numberOfFrequencies;
      //   getIndex(cg[grid].dimension(),I1,I2,I3);
      //   FOR_3D(i1,i2,i3,I1,I2,I3)
      //   {
      //     if( maskLocal(i1,i2,i3) > 0  )
      //     {
      //       bcgLocal(i1,i2,i3,Nf)=vLocal(i1,i2,i3,Nf);
      //     }
      //     else
      //     {
      //       bcgLocal(i1,i2,i3,Nf)=0.; 
      //     }
      //   }

      // }
      // bcg[grid]= v[grid];
        }

        if( false )
        {
            ::display(bcg[0],"b: iteration 0 ","%5.2f ");
        }

    //compute y = Pi * x 
    // set y = v 
        if( useMatrixUtilities )
        {
            cgWave.gridFunctionToVector( v,yl, iStart,iEnd );
        }
        else
        {    
            int i=0;
            for( int freq=0; freq<numCompWaveHoltz; freq++ )
            {      
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    getIndex(cg[grid].dimension(),I1,I2,I3);
                    OV_GET_SERIAL_ARRAY(int,cg[grid].mask(),maskLocal);

                    OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        assert( i<iEnd );
                        yl[i]= vLocal(i1,i2,i3,freq); // new way Jan 5, 2022

            // // if( true || maskLocal(i1,i2,i3)>0 )
            // // if( !checkMask || maskLocal(i1,i2,i3) !=0 )
            // if( true ) // Jan 5, 2022 
            // {
            //   yl[i]= vLocal(i1,i2,i3,freq);
            // }
            // else
            // {
            //   // un-used points: set [Ax]_i = x_i   -- is this right ??
            //   // yl[i]=xl[i];
            //   yl[i]=.0; // Make Ax=0 at unused points
            // }
                        
                        i++;
                    }
                }
            }
            assert( i==iEnd );
        }
    }
    else if( iteration==computeResidual )
    {
    // ---- iteration=computeResidual: (after final step usually) check the residual in the computed solution  : v^n - v^{n+1}
    // Compute :
    //       A v^n = v^n - v^{n+1} + b 

        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            vOld[grid] -= v[grid]; // 
    
        const Real & tol = cgWaveHoltz.dbase.get<Real>("tol");
        Real & maxResidual = cgWaveHoltz.dbase.get<Real>("maxResidual");

        maxResidual = maxNorm(vOld);
        printF("it=%d:  max(residual) = max(|v-vOld|)=%8.2e, tol=%g\n",iteration,maxResidual,tol);

    }
    else
    {
    // ----------- MATRIX VECTOR MULTIPLY FOR Krylov solver -------

        const int & computeEigenmodes       = cgWave.dbase.get<int>("computeEigenmodes");

        if( !computeEigenmodes )
        {
      // Compute :
      //       A v^n = v^n - v^{n+1} + b 
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                vOld[grid] -= v[grid]; // 


            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                vOld[grid] += bcg[grid];

        }
        else
        {
      // -- for eigenmodes -- assume b=0
      // Set A v^n = v^{n+1}
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                vOld[grid] = v[grid]; // 
        }

    // should we interpolate vOld?
        if( false ) //  Jan 5, 2022 : change to false : FIXES SIC
        {
            vOld.interpolate();
        }

        assert( pb !=NULL );
        Vec & b = *pb;
        PetscScalar *bl;
        VecGetArray(b,&bl);  // get the local array from Petsc

    // --- set y = vOld ---
        if( useMatrixUtilities )
        {
            cgWave.gridFunctionToVector( vOld,yl, iStart,iEnd );
        }  
        else
        {  
            int i=0;
            for( int freq=0; freq<numCompWaveHoltz; freq++ )
            {       
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    getIndex(cg[grid].dimension(),I1,I2,I3);
                    OV_GET_SERIAL_ARRAY(int,cg[grid].mask(),maskLocal);
                    OV_GET_SERIAL_ARRAY(Real,vOld[grid],vOldLocal);
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        assert( i<iEnd );
                        yl[i]= vOldLocal(i1,i2,i3,freq);    // y = A*x 

            // // yl[i]= vOldLocal(i1,i2,i3,0) + bl[i]; 
            // // if( checkMask && maskLocal(i1,i2,i3)!=0 )          // ************************* check me ***************
            // if( true ) // Jan 5, 2022
            // {
            //   yl[i]= vOldLocal(i1,i2,i3,freq);    // y = A*x 
            // }
            // else
            // {
            //   // yl[i]=xl[i];
            //   yl[i]=0.;  // unused point 
            // }
                        
                        i++;
                    }
                }
            }
            assert( i==iEnd );
        }


    
    }

    const int & computeEigenmodes = cgWave.dbase.get<int>("computeEigenmodes");
    if( monitorResiduals && iteration>=0 && !computeEigenmodes )
    {
    // ----- Optionally save the current residual -----
    // Save "residuals" by iteration: 
    // resVector(it) = norm( v^{n+1} - v^n )
        RealArray & resVector = cgWave.dbase.get<RealArray>("resVector");  

        if( iteration<0 )
        {
            printF("solvePETSC: ERROR: INVALID iteration=%d.\n",iteration);
            OV_ABORT("ERROR");
        }
        const int & maximumNumberOfIterations = cgWaveHoltz.dbase.get<int>("maximumNumberOfIterations");
        if( iteration>=resVector.getLength(0) ) 
        {
            printF("solvePETSC: Increasing the size of the vector that holds residuals, iteration=%d, maximumNumberOfIterations=%d\n",iteration,maximumNumberOfIterations);
            resVector.resize(resVector.getLength(0)*2); 
        }

    // resVector(iteration)= cgWaveHoltz.residual();   // this is not correct since solution is not the right one

        Real kspResidual;
        assert( ksp !=NULL );
        KSPGetResidualNorm(ksp,&kspResidual); 
        kspResidual /= sqrt(numberOfActivePoints);  // make an approximate L2h norm

        resVector(iteration)= kspResidual;
        printF("\n ##### SAVE KRYLOV RESIDUAL: iteration=%d: L2h-residual=%9.2e \n\n",kspResidual);
        
    // if( iteration <= maximumNumberOfIterations )
    // {
    //   resVector(iteration)= kspResidual;
    //   printF("\n ##### SAVE KRYLOV RESIDUAL: iteration=%d: L2h-residual=%9.2e \n\n",kspResidual);
    // }
    // else
    // {
    //    printF("\n ##### KRYLOV RESIDUAL: iteration=%d: L2h-residual=%9.2e (NOT SAVED since iteration=%d > maximumNumberOfIterations=%d)\n\n",
    //     kspResidual,iteration,maximumNumberOfIterations);
    // }

        if( useVariableTolerance ) // ** THIS DOESNT SEEM TO WORK -- KRYLOV NEEDS ACCURATE Matrix-vector multiplies ***
        {
            const CgWave::TimeSteppingMethodEnum & timeSteppingMethod = cgWave.dbase.get<CgWave::TimeSteppingMethodEnum>("timeSteppingMethod");
            if( timeSteppingMethod == CgWave::implicitTimeStepping )
            {
                Oges & impSolver = cgWave.dbase.get<Oges>("impSolver");
                if( impSolver.isSolverIterative() )
                {
                    Real rtol= max(1e-12,min(1.e-3,kspResidual/10));
          // rtol=1e-10;

                    Real atol=rtol;

                    printF("CgWaveHoltz:waveHoltzMatrixVectorMultiply: set implicit solver tolerance: rtol=%9.2e, atol=%9.2e (based on residual=%9.3e)\n",atol,rtol,kspResidual);
            
                    impSolver.set(OgesParameters::THErelativeTolerance,rtol);        
                    impSolver.set(OgesParameters::THEabsoluteTolerance,atol);   
                } 
            }    
        }

    // // There is no residual for iteration=0 since the Krylov solver is just computing the first A*x
    // if( iteration==1 )
    //   resVector(0) = resVector(iteration); // do this for now -- ksp is not built yet for iteration=0
    } 

    iteration++;
    
    return ierr;
}


// ============================================================================
/// \brief Solve for the Helmholtz solution using PETSc
// ============================================================================
int CgWaveHoltz::
solvePETSc(int argc,char **args)
{
    pCgWaveHoltz=this;

    const int myid=max(0,Communication_Manager::My_Process_Number);
    
    const Real & omega                    = dbase.get<Real>("omega");
    Real & Tperiod                        = dbase.get<Real>("Tperiod");
    const int & numPeriods                = dbase.get<int>("numPeriods");
    const int & adjustOmega               = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t   
    const int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");
    const int & useVariableTolerance      = dbase.get<int>("useVariableTolerance"); // Vary implicit solver tolerance based on current WaveHoltz residual

    const aString & krylovType            = dbase.get<aString>("krylovType");  // gmres, bicgstab, cg, ...
    const int & gmresRestartLength        = dbase.get<int>("gmresRestartLength"); // restart length for WaveHoltz + GMRES  

    Real & numberOfActivePoints           = dbase.get<Real>("numberOfActivePoints");


    
  // here is the CgWave solver for the time dependent wave equation
    CgWave & cgWave = *dbase.get<CgWave*>("cgWave");
    const int & numCompWaveHoltz          = cgWave.dbase.get<int>("numCompWaveHoltz");
    const int & filterTimeDerivative      = cgWave.dbase.get<int>("filterTimeDerivative");
    const int & upwind                    = cgWave.dbase.get<int>("upwind");
    const int & orderOfAccuracy           = cgWave.dbase.get<int>("orderOfAccuracy");
    const int numGhost = orderOfAccuracy/2;

    if( omega!=0. )
        Tperiod=numPeriods*twoPi/omega; 
    else 
        Tperiod=1.;
  
    printF("CgWaveHoltz::solvePETSc: setting tFinal = Tperiod*numPeriods = %9.3e (numPeriods=%d) \n",Tperiod,numPeriods);
  
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

    realCompositeGridFunction & vOld = cgWave.dbase.get<realCompositeGridFunction>("vOld");
    Range all;
  // const int & numCompWaveHoltz = cgWave.dbase.get<int>("numCompWaveHoltz");  
    vOld.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz);

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    int iab[2];  

    if( !dbase.has_key("bcg") )
    {
        realCompositeGridFunction & bcg = dbase.put<realCompositeGridFunction>("bcg");
        Range all;
        bcg.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz);
    }  

  // <<<<
    if( computeEigenmodes )
    {
        int & computeEigenmodesWithSLEPc = cgWave.dbase.get<int>("computeEigenmodesWithSLEPc");
        computeEigenmodesWithSLEPc = 1; // tell cgWave we are solving with SLEPc
    }

  // Save "residuals" by iteration: 
  // resVector(it) = norm( v^{n+1} - v^n )
    RealArray & resVector = cgWave.dbase.get<RealArray>("resVector");
    resVector.redim(maximumNumberOfIterations+10);
    resVector=0.;
    
    int & plotOptions = cgWave.dbase.get<int>("plotOptions");
    plotOptions= CgWave::noPlotting; // turn of plotting in cgWave  

    Vec            x,b,u;  /* approx solution, RHS, exact solution */
    Mat            A;        /* linear system matrix */
  //  KSP            ksp;     /* linear solver context */
    PetscRandom    rctx;     /* random number generator context */
    PetscReal      norm;     /* norm of solution error */
    PetscInt       i,j,Ii,J,Istart,Iend,m = 8,n = 7,its;
    PetscErrorCode ierr;
    PetscBool      flg = PETSC_FALSE;
    PetscScalar    v;
    #if defined(PETSC_USE_LOG)
        PetscLogStage  stage;
    #endif

    int & petscIsInitialized = dbase.get<int>("petscIsInitialized");
    
    if( !petscIsInitialized )
    {
        petscIsInitialized=true;
        if( !computeEigenmodes )
        {
            PetscInitialize(&argc,&args,(char *)0,help);
        }
        else
        {
            #ifdef CGWAVE_USE_SLEPC
                SlepcInitialize(&argc,&args,(char *)0,help);
            #else
                OV_ABORT("You need SLEPc to solve for eigenvalues");
            #endif
        }

    // 3.4.5 ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
    // 3.4.5 ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
    // 3.18.2 : PetscErrorCode PetscOptionsGetInt(PetscOptions options, const char pre[], const char name[], PetscInt *ivalue, PetscBool *set)
    // ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
    // ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
    }
    

  // ============== Create a shell object for matrix free method ================
  // double mycontext; // passed to waveHoltzMatrixVectorMultiply 
    Mat *mycontext=&A;
    
  // ******  global variables
    pA = &A; mA=m; nA=n;


  // *new* way useActivePoints only for CG and biCGStab : Nov 15, 2024
    bool useActivePoints=true; 
    if( upwind )
        useActivePoints=false; 
    
    int numEquations     =0;
    int numEquationsLocal=0;  // number of equations local to this processor 

    assert( pCgWaveHoltz!=NULL );
    CgWaveHoltz & cgWaveHoltz = *pCgWaveHoltz;
    CompositeGrid & cg = cgWaveHoltz.cg;

    if( useMatrixUtilities )    
    {
        bool checkMask=false; // do this for now to match old way 
        cgWave.initializeGlobalIndexing( checkMask );
        const int & totalActive    = cgWave.dbase.get<int>("totalActive");
        const int & numActiveLocal = cgWave.dbase.get<int>("numActiveLocal");

        numEquations = totalActive;
        numEquations *= numCompWaveHoltz;

        numEquationsLocal = numActiveLocal*numCompWaveHoltz;

        numberOfActivePoints = totalActive;

    }
    else
    {
    // Index I1,I2,I3;
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            const IntegerArray & gid = mg.gridIndexRange();      
            if( useActivePoints )
                {
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
                            else if( upwind || bc==CgWave::neumann )
                            {
                // include ghost ??
                // iab[side] -= is;
                            }
                            else if( bc==CgWave::abcEM2  || bc==CgWave::absorbing || bc==CgWave::radiation )
                            {
                // include ghost 
                                iab[side] -= is*numGhost;
                            }      
                            else if( bc==CgWave::dirichlet )
                            {
                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                 // iab[side] -= is*numGhost; // ************************************** TEMP ***********
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
            else
                getIndex(mg.dimension(),I1,I2,I3);

            numEquations += I1.getLength()*I2.getLength()*I3.getLength();
        }
        numEquations *= numCompWaveHoltz;
        numEquationsLocal = numEquations;
    }

    printf("Make a Matrix Free Shell: myid=%d: numEquationsLocal=%d, numEquations=%d\n",myid, numEquationsLocal, numEquations);

    
  // PetscErrorCode MatCreateShell(MPI_Comm comm, PetscInt m, PetscInt n, PetscInt M, PetscInt N, void *ctx, Mat *A)
  // Input Parameters
  // comm - MPI communicator

  // m - number of local rows (or PETSC_DECIDE to have calculated if M is given)
  // n - number of local columns (or PETSC_DECIDE to have calculated if N is given)
  // M - number of global rows (may be PETSC_DETERMINE to have calculated if m is given)
  // N - number of global columns (may be PETSC_DETERMINE to have calculated if n is given)
  // ctx - pointer to data needed by the shell matrix routines

  // Output Parameter
  // A - the matrix

    Mat Amf;

  // ierr = MatCreateShell(PETSC_COMM_WORLD, numEquations, numEquations, PETSC_DECIDE,  PETSC_DECIDE, mycontext, &Amf);
    ierr = MatCreateShell(PETSC_COMM_WORLD, numEquationsLocal, numEquationsLocal, numEquations,  numEquations, mycontext, &Amf);

    if( useActivePoints )
        ierr = MatShellSetOperation(Amf, MATOP_MULT, (void(*)(void))waveHoltzMatrixVectorMultiply);
    else
        ierr = MatShellSetOperation(Amf, MATOP_MULT, (void(*)(void))waveHoltzMatrixVectorMultiplyOld);
    CHKERRQ(ierr);

    /* 
          Create parallel vectors.
          - We form 1 vector from scratch and then duplicate as needed.
          - When using VecCreate(), VecSetSizes and VecSetFromOptions()
          in this example, we specify only the
          vector's global dimension; the parallel partitioning is determined
          at runtime. 
          - When solving a linear system, the vectors and matrices MUST
          be partitioned accordingly.  PETSc automatically generates
          appropriately partitioned matrices and vectors when MatCreate()
          and VecCreate() are used with the same communicator.  
          - The user can alternatively specify the local vector and matrix
          dimensions when more sophisticated partitioning is needed
          (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
          below).
    */

    ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
    ierr = VecSetSizes(u,numEquationsLocal,numEquations);CHKERRQ(ierr);
    ierr = VecSetFromOptions(u);CHKERRQ(ierr);
    ierr = VecDuplicate(u,&b);CHKERRQ(ierr); 
    ierr = VecDuplicate(b,&x);CHKERRQ(ierr);

    pb = &b;  // global variable for now
    
    PetscReal bNorm; // save l2-norm of RHS b

  
    if( false )
    {
    // -- zero initial guess 
        ierr = VecSet(x,0.0);CHKERRQ(ierr);

    }
    else 
    {
    // --- Set initial guess to be current v ---

        PetscScalar *xl;
        VecGetArray(x,&xl);  // get the local array from Petsc
        int iStart,iEnd;
        ierr = VecGetOwnershipRange(x,&iStart,&iEnd); CHKERRQ(ierr);
    // printf("  [iStart,iEnd]=[%i,%i]\n",iStart,iEnd);
    
    // Helmholtz solution is stored here:
        realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
        CompositeGrid & cg = *v.getCompositeGrid();

    // Index I1,I2,I3;

        if( useMatrixUtilities )
        {
            cgWave.gridFunctionToVector( v,xl, iStart,iEnd );
        }
        else
        {
            numberOfActivePoints = 0.; // count active points for scaling norm.
            Real normV=0; 
            int i=0;
            for( int freq=0; freq<numCompWaveHoltz; freq++ )
            {
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cg[grid];
                    const IntegerArray & gid = mg.gridIndexRange();          
                    if( useActivePoints )
                        {
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
                                    else if( upwind || bc==CgWave::neumann )
                                    {
                    // include ghost ??
                    // iab[side] -= is;
                                    }
                                    else if( bc==CgWave::abcEM2  || bc==CgWave::absorbing || bc==CgWave::radiation )
                                    {
                    // include ghost 
                                        iab[side] -= is*numGhost;
                                    }      
                                    else if( bc==CgWave::dirichlet )
                                    {
                                          iab[side] += is;  // Dirichlet BC -- ignore the boundary
                     // iab[side] -= is*numGhost; // ************************************** TEMP ***********
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
                    else
                        getIndex(mg.dimension(),I1,I2,I3);

                    OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                    OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
                        assert( i<iEnd );
                        if( maskLocal(i1,i2,i3) !=0 )
                            numberOfActivePoints++;

                        xl[i] = vLocal(i1,i2,i3,freq);

            // if( true || (checkMask && maskLocal(i1,i2,i3)!=0) )
            //   xl[i] = vLocal(i1,i2,i3,freq);
            // else
            //   xl[i]=0.;

                        normV = max( normV,xl[i]);
                        i++;
                    }
                }
            }
            assert( i==iEnd );
            numberOfActivePoints = ParallelUtility::getSum(numberOfActivePoints);

            printF("solvePETSc: Set initial guess to v, max-norm(v)=%8.2e, numberOfActivePoints=%g\n",normV,numberOfActivePoints);
        }

    } // END set initial guess 

  // ---- set RHS b = A*u -----
    ierr = VecSet(u,0.0);CHKERRQ(ierr);

  // This next call to MatMult will eventually call CgWave with initial condition u=0
  // Warning: the next call will over-write v, so we need to save v in PETSc vector x before we get to here

    iteration=computeRightHandSide;
    ierr = MatMult(Amf,u,b);CHKERRQ(ierr); 

    VecNorm(b,NORM_2,&bNorm);
    Real bNorm2h = bNorm/sqrt(numberOfActivePoints); 
    printF("solvePETSc: RHS is b: l2-norm(b)=%9.3e, L2h-norm(b)=%9.2e\n",bNorm,bNorm2h);

  // *wdh* Nov 29, 2024 -- make iteration==0 be the first call
  // iteration=0; 
    iteration=-1; // Set to -1 to start the true iterations (first call to A*x does not generate a residual)


    /*
          View the exact solution vector if desired
    */
    flg  = PETSC_FALSE;
  // 3.4.5 ierr = PetscOptionsGetBool(PETSC_NULL,"-view_exact_sol",&flg,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetBool(PETSC_NULLPTR,PETSC_NULLPTR,"-view_exact_sol",&flg,PETSC_NULLPTR);CHKERRQ(ierr);
    if (flg) {ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                Create the linear solver and set various options
          - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /* 
          Create linear solver context
    */
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

    KSPType krylovSpaceMethod = KSPGMRES; // default Krylov;
    if( krylovType=="gmres" ) 
    {
          krylovSpaceMethod = KSPGMRES;
    }
    else if( krylovType=="fgmres") 
    {
    // flexible GMRES 
        krylovSpaceMethod = KSPFGMRES;
    }  
    else if( krylovType=="bicgstab" || krylovType=="bcgs" || krylovType=="bicgs" )
    {
        krylovSpaceMethod = KSPBCGS;
        printF("CghWaveHoltz::solvePETSc: setting Krylov solver to bi-conjugate-gradient stabilized.\n");
    }
    else if( krylovType=="cg" )
    {
        krylovSpaceMethod = KSPCG;
        printF("CghWaveHoltz::solvePETSc: setting Krylov solver to conjugate-gradient.\n");
    }  
    else
    {
        printf("CgWaveHoltz::solvePETSc:WARNING: unknown krylovType=[%s].\n"
                      "  Valid options: gmres, bicgstab\n"
                      "  Continuing with gmres...\n",(const char*)krylovType);
        OV_ABORT("error");
    }

  // See Oges/PETScEquationSolver.C line 390
      // krylovSpaceMethod=KSPRICHARDSON;
      // krylovSpaceMethod=KSPCHEBYSHEV;
      // krylovSpaceMethod=KSPCG;
      // krylovSpaceMethod=KSPGMRES;
      // krylovSpaceMethod=KSPBCGS;
      // krylovSpaceMethod=KSPTCQMR;
      // krylovSpaceMethod=KSPTFQMR;
      // krylovSpaceMethod=KSPCR;
      // krylovSpaceMethod=KSPLSQR;
      // krylovSpaceMethod=KSPPREONLY;
      // krylovSpaceMethod=KSPQCG;
      // krylovSpaceMethod=KSPBICG;
      // krylovSpaceMethod=KSPCGS;

  // Set the Kylov solver such as gmres, etc.
    ierr = KSPSetType(ksp, krylovSpaceMethod); CHKERRQ(ierr);

    if( gmresRestartLength>0 )
    {
        printF("CghWaveHoltz::PETScSolver: setting gmres restart length to %i\n",gmresRestartLength);
        KSPGMRESSetRestart(ksp, gmresRestartLength); CHKERRQ(ierr);
    }

    /* 
          Set operators. Here the matrix that defines the linear system
          also serves as the preconditioning matrix.
    */
  //   ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  // 3.4.5 ierr = KSPSetOperators(ksp,Amf,Amf,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  // 3.18.2: 
    ierr = KSPSetOperators(ksp,Amf,Amf);CHKERRQ(ierr);


    /* 
          Set linear solver defaults for this problem (optional).
          - By extracting the KSP and PC contexts from the KSP context,
              we can then directly call any KSP and PC routines to set
              various options.
          - The following two statements are optional; all of these
              parameters could alternatively be specified at runtime via
              KSPSetFromOptions().  All of these defaults can be
              overridden at runtime, as indicated below.
    */
    const Real & tol = dbase.get<Real>("tol");
    
  // PETSc relative tolerance is based on the l2-norms:
  //        l2Norm(res)/l2Norm(b) < relativeTol
  // We want
  //        L2hNorm(res) = l2Norm(res)/sqrt(N) < tol
  // so 
  //     relativeTol = tol*sqrt(N)/l2Norm(b)

    assert( numberOfActivePoints>0. );
    const Real relativeTol = tol*sqrt(numberOfActivePoints)/max(REAL_MIN*1000.,bNorm);

  // KSPSetTolerances(KSP ksp,PetscReal rtol,PetscReal abstol,PetscReal dtol,PetscInt maxits)
    ierr = KSPSetTolerances(ksp,relativeTol,1.e-50,PETSC_DEFAULT,maximumNumberOfIterations); CHKERRQ(ierr);

    /* 
        Set runtime options, e.g.,
                -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
        These options will override those specified above as long as
        KSPSetFromOptions() is called _after_ any other customization
        routines.
    */
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    PetscBool trueFlag=PETSC_TRUE;  // the initial guess is non-zero
    ierr = KSPSetInitialGuessNonzero(ksp,trueFlag); CHKERRQ(ierr);       

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                            Solve the linear system
          - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = KSPSolve(ksp,b,x); CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                                            Check solution and clean up
          - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    printF("***** DONE KSP SOLVE ******\n");




    int & numberOfIterations = dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken
    numberOfIterations = iteration;

  // save num mat-vects: 
    int & numberOfMatrixVectorMultiplications = dbase.get<int>("numberOfMatrixVectorMultiplications");
    numberOfMatrixVectorMultiplications = numberOfIterations;

    if( false )
    { // OLD WAY -- requires one additional call to cgWave.advance -- eliminate this Nov 30, 2024.
        printF("***** CALL AGAIN TO ASSIGN V and CHECK RESIDUAL ******\n");


        iteration=computeResidual; // This tells the matrix-vector multiply routine we are computing the residual
    // -- Here we pass the final solution ---
        ierr = MatMult(Amf,x,b);CHKERRQ(ierr);     


        Real kspResidual;
        KSPGetResidualNorm(ksp,&kspResidual);
        kspResidual /= sqrt(numberOfActivePoints);  // make an approximate L2h norm

    // add final residual
        resVector(numberOfIterations) = kspResidual;
        numberOfIterations++; 


        const Real & maxResidual = dbase.get<Real>("maxResidual");
        printF("\n ######## DONE KYRLOV ITERATIONS -- KSP residual=%8.2e (max-res=%8.2e) (tol=%8.2e) numIts=%i #######\n",
                      kspResidual,maxResidual,tol,numberOfIterations);
    }
    else
    {
        printF("***** CALL AGAIN TO ASSIGN V  ******\n");


        iteration=assignSolution; // This tells the matrix-vector multiply routine we are just assigning the solution
    // -- Here we pass the final solution ---
        ierr = MatMult(Amf,x,b);CHKERRQ(ierr);     


        Real kspResidual;
        KSPGetResidualNorm(ksp,&kspResidual);
        kspResidual /= sqrt(numberOfActivePoints);  // make an approximate L2h norm

        const Real & maxResidual = dbase.get<Real>("maxResidual");
        printF("\n ######## DONE KYRLOV ITERATIONS -- KSP residual=%8.2e (max-res=%8.2e) (tol=%8.2e) numIts=%i #######\n",
                      kspResidual,maxResidual,tol,numberOfIterations);


    }

    KSPType kspTypeName;
    KSPGetType( ksp, &kspTypeName );
    PC pc;
    ierr = KSPGetPC( ksp, &pc);    CHKERRQ(ierr);
    PCType  pcTypeName;
    PCGetType( pc, &pcTypeName );
    printf(" @@@@ INFO from PETSc: kspType=[%s], pcType=[%s]\n",kspTypeName,pcTypeName);

    /*
          Free work space.  All PETSc objects should be destroyed when they
          are no longer needed.
    */
    ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    ierr = VecDestroy(&u);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr);
    ierr = VecDestroy(&b);CHKERRQ(ierr);
    if( false )
    {
        ierr = MatDestroy(&A);CHKERRQ(ierr);
    }
    
  // --- Compute the average convergence rate ----
    Real & convergenceRate          = dbase.get<Real>("convergenceRate");
    Real & convergenceRatePerPeriod = dbase.get<Real>("convergenceRatePerPeriod");
    
  // const int & numberOfIterations = dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken  
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

    /*
          Always call PetscFinalize() before exiting a program.  This routine
              - finalizes the PETSc libraries as well as MPI
              - provides summary and diagnostic information if certain runtime
                  options are chosen (e.g., -log_summary). 
    */

  // *** DO NOT CALL PetscFinalize HERE -- must be done at very end 
//  ierr = PetscFinalize();
    return 0;
}
