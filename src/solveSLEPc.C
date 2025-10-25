// This file automatically generated from solveSLEPc.bC with bpp.
#include "CgWaveHoltz.h" 
#include "CgWave.h" 
#include "Overture.h"
#include "gridFunctionNorms.h"
#include "display.h"
#include "ParallelUtility.h"
#include "SparseRep.h" 
#include "CompositeGridOperators.h"    


#include "PlotStuff.h"
#include "GL_GraphicsInterface.h"


// for petsc 3.18.2
#include "petscksp.h"

// SLEPc 
#include <slepceps.h>

static char help[] = "Solve for eigemodes usig SLEPc\n";

static bool useMatrixUtilities=true; // true = use new matrix utilities (needed for parallel)


static int iteration=0;

static CgWaveHoltz *pCgWaveHoltz; // pointer to the CgWaveHoltz solver 
  

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR3N(i1,i2,i3,n,n1a,n1b,n2a,n2b,n3a,n3b)       for( i3=n3a; i3<=n3b; i3++ )                        for( i2=n2a; i2<=n2b; i2++ )                      for( i1=n1a; i1<=n1b; i1++ )                    for( n=0; n<numberOfComponents; n++ )


// -- global variables -- do this for now 
// static KSP  ksp=NULL;     /* linear solver context */
// static Mat *pA=NULL;
// static Vec *pb=NULL;
// static int mA,nA;

// static int computeRightHandSide =-2; // = 0;
// static int computeResidual      =-3; // =-1;

// *wdh* Jan 5, 2022 --> changed cgWave to zero out unused points => thus we do not need to check the mask
// static bool checkMask = false;  // if true, do not use values of the solution where mask==0 

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



// --------------------------------------------------------------------------------------
//   Macro: return the index's for possible active points
//            
//  NOTE: This macro appears in solveSLEPc.bC and eigenModes.bC 
// --------------------------------------------------------------------------------------


// EPSMonitorSet(EPS eps,PetscErrorCode (*monitor)(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *mctx),void *mctx,PetscErrorCode (*monitordestroy)(void**))

// User defined monitor for SLEPC
PetscErrorCode eigenWaveMonitor(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *mctx)
{
    printF("\n"
              " *************************************************************************************************\n"
              " *****************  eigenWaveMonitor: Outer iteration=%d, Number converged=%d  *******************\n"
              " *************************************************************************************************\n"
              ,its,nconv);

    return 0;
}




// =========================================================================================================
//     MATRIX-VECTOR MULTIPLY FOR MATRIX FREE KRYLOV SOLVERS
// 
// Compute y = M*x 
// =========================================================================================================
extern PetscErrorCode eigenWaveMatrixVectorMultiply(Mat m ,Vec x, Vec y)
{
    PetscErrorCode ierr = 0;

    printF("\n ++++++++ EigenWave Matrix vector multiply routine: eigenWaveMatrixVectorMultiply called iteration=%d +++\n",iteration);
    if( debug>2 )
        printF(" **************** MatVec for SLEPc iteration=%i *************\n",iteration);

    assert( pCgWaveHoltz!=NULL );
    CgWaveHoltz & cgWaveHoltz = *pCgWaveHoltz;
    
    const int & monitorResiduals      = cgWaveHoltz.dbase.get<int>("monitorResiduals");      // montior the residuals at every step
    const Real & numberOfActivePoints = cgWaveHoltz.dbase.get<Real>("numberOfActivePoints");
  

  // here is the CgWave solver for the time dependent wave equation
    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
    const int & numberOfFrequencies = cgWave.dbase.get<int>("numberOfFrequencies");
    const int & upwind               = cgWave.dbase.get<int>("upwind");

    const int & computeEigenmodes = cgWave.dbase.get<int>("computeEigenmodes");
    assert( computeEigenmodes );

  // -- cgWave solution is stored here: 
    realCompositeGridFunction & v    = cgWave.dbase.get<realCompositeGridFunction>("v");

    realCompositeGridFunction & vOld = cgWave.dbase.get<realCompositeGridFunction>("vOld");
    CompositeGrid & cg = cgWaveHoltz.cg;

  // if( !cgWaveHoltz.dbase.has_key("bcg") )
  // {
  //   realCompositeGridFunction & bcg = cgWaveHoltz.dbase.put<realCompositeGridFunction>("bcg");
  //   Range all;
  //   bcg.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);
  // }
  // realCompositeGridFunction & bcg = cgWaveHoltz.dbase.get<realCompositeGridFunction>("bcg");    


    PetscScalar *xl, *yl;
    VecGetArray(x,&xl);  // get the local array from Petsc
    VecGetArray(y,&yl);  // get the local array from Petsc
    int iStart,iEnd;
    ierr = VecGetOwnershipRange(x,&iStart,&iEnd);CHKERRQ(ierr);
  // printf("  [iStart,iEnd]=[%i,%i]\n",iStart,iEnd);


  //   ---- vectorToGridFunction ----
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
        for( int freq=0; freq<numberOfFrequencies; freq++ )
        {
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const IntegerArray & gid = mg.gridIndexRange();

                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);


        // // // getIndex(cg[grid].dimension(),I1,I2,I3);
        // int extra=-1; // *********************************************** FIX ME
        // getIndex(cg[grid].gridIndexRange(),I1,I2,I3,extra);

                    Iv[2]=Range(0,0);
                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                    {
                        for( int side=0; side<=1; side++ )
                        {
                            int is = 1-2*side;
                            iab[side]=gid(side,axis);
                            const int bc = mg.boundaryCondition(side,axis);
                            if( bc==CgWave::dirichlet )
                            {
                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                            }
                            else if( bc==CgWave::neumann || bc==CgWave::abcEM2  || bc==CgWave::absorbing || bc==CgWave::radiation )
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

                OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
                    if( maskLocal(i1,i2,i3) > 0 )
                    {
                        assert( i<iEnd );

                        vLocal(i1,i2,i3,freq)=xl[i]; // new way  Jan 5, 2022

                        i++;
                    }
                }
            }
        }
        assert( i==iEnd );
    }

  // *** APPLY BOUNDARY CONDITIONS to v  **************** IS THIS NEEDED ?? Doesn't advance do this ???
    Real t=0.; // should not matter for eigenvalue problem
    cgWave.applyEigenFunctionBoundaryConditions( v );

  // Is this needed? 
  // for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  //   v[grid].updateGhostBoundaries();  

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
    cgWave.advance( iteration );


  // Set A v^n = v^{n+1}
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        vOld[grid] = v[grid]; 

  // // should we interpolate vOld?
  // if( false ) //  Jan 5, 2022 : change to false : FIXES SIC
  // {
  //   vOld.interpolate();
  // }

  // --- gridFunctionToVector ---
  // set y = v 
    if( useMatrixUtilities )
    {
        cgWave.gridFunctionToVector( vOld,yl, iStart,iEnd );
    }
    else
    {  
        int i=0;
        for( int freq=0; freq<numberOfFrequencies; freq++ )
        {       
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const IntegerArray & gid = mg.gridIndexRange();      

                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                OV_GET_SERIAL_ARRAY(Real,vOld[grid],vOldLocal);

        // int extra=-1; // *********************************************** FIX ME
        // getIndex(cg[grid].gridIndexRange(),I1,I2,I3,extra);

                    Iv[2]=Range(0,0);
                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                    {
                        for( int side=0; side<=1; side++ )
                        {
                            int is = 1-2*side;
                            iab[side]=gid(side,axis);
                            const int bc = mg.boundaryCondition(side,axis);
                            if( bc==CgWave::dirichlet )
                            {
                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                            }
                            else if( bc==CgWave::neumann || bc==CgWave::abcEM2  || bc==CgWave::absorbing || bc==CgWave::radiation )
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

                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
                    if( maskLocal(i1,i2,i3) > 0 )
                    {
                        assert( i<iEnd );
                        yl[i]= vOldLocal(i1,i2,i3,freq);    // y = A*x 

                        i++;
                    }
                }
            }
        }
        assert( i==iEnd );
    }


    iteration++;
    
    return ierr;
}


// ============================================================
// Copy a component of a grid function to another
//      u[iu] <- v[iv]
// ============================================================


// ============================================================================
/// \brief Solve for eigenmodes using SLEPc
// ============================================================================
int CgWaveHoltz::
solveSLEPc(int argc,char **args)
{
    pCgWaveHoltz=this;

    const int myid=max(0,Communication_Manager::My_Process_Number);
    
    const int & debug                     = dbase.get<int>("debug");
    const int & cgWaveDebugMode           = dbase.get<int>("cgWaveDebugMode");

    const Real & tol                      = dbase.get<Real>("tol");
    const Real & omega                    = dbase.get<Real>("omega");
    Real & Tperiod                        = dbase.get<Real>("Tperiod");
    const int & numPeriods                = dbase.get<int>("numPeriods");
    const int & adjustOmega               = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t   
    const int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");
    Real & numberOfActivePoints           = dbase.get<Real>("numberOfActivePoints");

    const int numberOfDimensions = cg.numberOfDimensions();
    
  // here is the CgWave solver for the time dependent wave equation
    CgWave & cgWave = *dbase.get<CgWave*>("cgWave");
    if( omega!=0. )
        Tperiod=numPeriods*twoPi/omega; 
    else 
        Tperiod=1.;
  

    FILE *& debugFile  = cgWave.dbase.get<FILE*>("debugFile");
    FILE *& pDebugFile = cgWave.dbase.get<FILE*>("pDebugFile");

    printF("CgWaveHoltz::solveSLEPC: setting tFinal = Tperiod*numPeriods = %9.3e (numPeriods=%d) \n",Tperiod,numPeriods);
  
  
    const int & computeEigenmodes       = cgWave.dbase.get<int>("computeEigenmodes");
    const CgWave::EigenSolverInitialConditionEnum & eigenSolverInitialCondition 
                                                                            = cgWave.dbase.get<CgWave::EigenSolverInitialConditionEnum>("eigenSolverInitialCondition");

  // --- set values in CgWave:  *** COULD DO BETTER ***

  // SOME OF THIS IS CURRENTLY NEEDED:

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
    const int & numCompWaveHoltz = cgWave.dbase.get<int>("numCompWaveHoltz");  
    vOld.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz);

  // vOld.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);
    CompositeGridOperators & operators = cgWave.dbase.get<CompositeGridOperators>("operators");
    vOld.setOperators(operators);


    int & computeEigenmodesWithSLEPc = cgWave.dbase.get<int>("computeEigenmodesWithSLEPc");
    computeEigenmodesWithSLEPc = 1; // tell cgWave we are solving with SLEPc

  // Save "residuals" by iteration: 
  // resVector(it) = norm( v^{n+1} - v^n )
  // RealArray & resVector = cgWave.dbase.get<RealArray>("resVector");
  // resVector.redim(maximumNumberOfIterations+10);
  // resVector=0.;
    
    int & plotOptions = cgWave.dbase.get<int>("plotOptions");
    plotOptions= CgWave::noPlotting; // turn of plotting in cgWave  

    if( cgWaveDebugMode )
        plotOptions= CgWave::plotAndWait;  

  // Vec            x,b,u;  /* approx solution, RHS, exact solution */
  // Mat            A;        /* linear system matrix */
  //  KSP            ksp;     /* linear solver context */

  // PetscRandom    rctx;     /* random number generator context */
  // PetscReal      norm;     /* norm of solution error */

  // PetscInt       i,j,Ii,J,Istart,Iend,m = 8,n = 7,its;
    PetscInt       i,j,Istart,Iend,its;

    PetscErrorCode ierr;
  //PetscScalar    v;
  //#if defined(PETSC_USE_LOG)
  //  PetscLogStage  stage;
  //#endif

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
            SlepcInitialize(&argc,&args,(char *)0,help);
        }

    // 3.4.5 ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
    // 3.4.5 ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
    // 3.18.2 : PetscErrorCode PetscOptionsGetInt(PetscOptions options, const char pre[], const char name[], PetscInt *ivalue, PetscBool *set)
    // ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
    // ierr = PetscOptionsGetInt(PETSC_NULL,PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
    }
    
    
    assert( pCgWaveHoltz!=NULL );
    CgWaveHoltz & cgWaveHoltz = *pCgWaveHoltz;
    CompositeGrid & cg = cgWaveHoltz.cg;

  // --- count the number of ACTIVE equations ---
  //  These are equations that NOT constraint equations, i.e.
  //     Lap_h U = k U ,  
  //   k = eigenvalue = - lam^2 
  // 
    int iv[3], &i1=iv[0], &i2=iv[1], &i3=iv[2]; 

    Index Iv[3], &I1=Iv[0], &I2=Iv[1], &I3=Iv[2];
    int iab[2];

    int numEquations=0;
    int numEquationsLocal=0;  // number of equations local to this processor 

    if( useMatrixUtilities )    
    {
        bool checkMask=true; 
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
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            const IntegerArray & gid = mg.gridIndexRange();

            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

      // // getIndex(cg[grid].dimension(),I1,I2,I3);
      // int extra=-1;  // ******** DO THIS FOR NOW --- *********** FIX ME : not true for Neumann BC
      // getIndex(mg.gridIndexRange(),I1,I2,I3,extra);

      // printF("OLD: I1=[%d,%d] I2=[%d,%d] I3=[%d,%d]\n",
      //       Iv[0].getBase(),Iv[0].getBound(),
      //       Iv[1].getBase(),Iv[1].getBound(),
      //       Iv[2].getBase(),Iv[2].getBound());

                Iv[2]=Range(0,0);
                for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                {
                    for( int side=0; side<=1; side++ )
                    {
                        int is = 1-2*side;
                        iab[side]=gid(side,axis);
                        const int bc = mg.boundaryCondition(side,axis);
                        if( bc==CgWave::dirichlet )
                        {
                              iab[side] += is;  // Dirichlet BC -- ignore the boundary
                        }
                        else if( bc==CgWave::neumann || bc==CgWave::abcEM2  || bc==CgWave::absorbing || bc==CgWave::radiation )
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
            printF("getActivePointIndex: I1=[%d,%d] I2=[%d,%d] I3=[%d,%d]\n",
                        Iv[0].getBase(),Iv[0].getBound(),
                        Iv[1].getBase(),Iv[1].getBound(),
                        Iv[2].getBase(),Iv[2].getBound());    

            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
                if( maskLocal(i1,i2,i3) > 0 )
                {
                    numEquations++;
                }
            }
        }
        numEquations *= numberOfFrequencies;
    }

    printF("solveSLEPC:Make a Matrix Free Shell: numEquations=%d\n",numEquations);
  // OV_ABORT("stop here for now");

  // ============== Create a shell object for matrix free method ================
    Mat Amf;
    Mat *mycontext=&Amf;   // passed to waveHoltzMatrixVectorMultiply 

    PetscCall(MatCreateShell(PETSC_COMM_WORLD, numEquationsLocal, numEquationsLocal, numEquations, numEquations, mycontext, &Amf)); // destroy me 

    PetscCall(MatShellSetOperation(Amf, MATOP_MULT, (void(*)(void))eigenWaveMatrixVectorMultiply));

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
  // PetscCall(VecCreate(PETSC_COMM_WORLD,&u));                 // destroy me 
  // PetscCall(VecSetSizes(u,PETSC_DECIDE,numEquations));
  // PetscCall(VecSetFromOptions(u));
  // PetscCall(VecDuplicate(u,&b));                              // destroy me 
  // PetscCall(VecDuplicate(u,&x));                              // destroy me -- is this needed?
  
  // pb = &b;  // global variable for now                                  // destroy me , is this needed
    
    PetscReal bNorm; // save l2-norm of RHS b

    /* 
          Set exact solution; then compute right-hand-side vector.
          By default we use an exact solution of a vector with all
          elements of 1.0;  Alternatively, using the runtime option
          -random_sol forms a solution vector with random components.
    */
  


  // iteration=-1; // Set to -1 to start the true iterations (first call to A*x does not generate a residual)

    iteration=0; // is this correct ?


  // ====================================================
  // =============== COMPUTE EIGENMODES =================
  // ====================================================

    int & numEigsToCompute           = cgWave.dbase.get<int>("numEigsToCompute"); // number of eigenpairs to compute

  // -- turn OFF adjsuements for upwinding :
    int & adjustHelmholtzForUpwinding= cgWave.dbase.get<int>("adjustHelmholtzForUpwinding");
    int adjustHelmholtzForUpwindingSave = adjustHelmholtzForUpwinding;
    adjustHelmholtzForUpwinding=0; 
  // adjustHelmholtzForUpwinding=1;  // ** TESTING**

  // int & computeEigenmodesWithSLEPc = cgWave.dbase.get<int>("computeEigenmodesWithSLEPc");
  // computeEigenmodesWithSLEPc = 1; // tell cgWave we are solving with SLEPc

    EPS            eps;             /* eigenproblem solver context */

  //    Create eigensolver context

    PetscCall(EPSCreate(PETSC_COMM_WORLD,&eps));                              // destroy me 


  // Set operators. In this case, it is a standard eigenvalue problem

    PetscCall(EPSSetOperators(eps,Amf,NULL));

  // HEP = Hermitian
  // NHEP = non-Hermitian
  // PetscCall(EPSSetProblemType(eps,EPS_HEP)); 
    PetscCall(EPSSetProblemType(eps,EPS_NHEP)); 

  // Set any command line options:
    PetscCall(EPSSetFromOptions(eps));

    PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE));

    ST st; // spectral transform
    EPSGetST(eps,&st); 

    const CgWave::EigenSolverEnum & eigenSolver  = cgWave.dbase.get<CgWave::EigenSolverEnum>("eigenSolver"); 
  // -- default method is Krylov-Schur
    if( eigenSolver==CgWave::defaultEigenSolver || eigenSolver==CgWave::KrylovSchurEigenSolver )
    {
        PetscCall(EPSSetType(eps,EPSKRYLOVSCHUR));
    }
    else if( eigenSolver==CgWave::powerEigenSolver )
    {
    // --- plain power method ---
        PetscCall(EPSSetType(eps,EPSPOWER));
    } 
    else if( eigenSolver==CgWave::subspaceIterationEigenSolver )
    {
    // --- subspace iteration (block power) ---
        PetscCall(EPSSetType(eps,EPSSUBSPACE));
    }   
    else if( eigenSolver==CgWave::inverseIterationEigenSolver )
    {
    // -- inverse iteration -- but need to us an iterative solver to invert the shifted matrix (matric free)
        PetscCall(EPSSetType(eps,EPSPOWER));
        PetscCall(STSetType(st,STSINVERT));   // shift and invert

        if( false )
        {
            PetscScalar target=1.; // *************************************************** FIX ME 
            PetscCall( EPSSetTarget(eps,target ) );
        }

        PetscCall(EPSSetWhichEigenpairs(eps,EPS_TARGET_MAGNITUDE) );  // DO THIS ---

    // PetscScalar shift=1.; 
    // PetscCall(STSetShift(st,shift) ) ;

    // printF("\n >>>>>>>> Set SLEPSc Spectral transform = shit-and-invert, target=%9.2e <<<<<<<<<<<< \n",target);

    }     
    else if( eigenSolver==CgWave::JacobiDavidsonEigenSolver )
    {
        PetscCall(EPSSetType(eps,EPSJD));
    }    
    else if( eigenSolver==CgWave::ArnoldiEigenSolver )
    {
    // Arnoldi: with explicit restart and deflation.
        PetscCall(EPSSetType(eps,EPSARNOLDI));
    }
    else if( eigenSolver==CgWave::ArpackEigenSolver )
    {
    // AARPACK
        PetscCall(EPSSetType(eps,EPSARPACK));
    }
    else
    {
        OV_ABORT("error -- unknown eigenSolver");
    }

  // Arpack: --> unknown type? May need to install arpack too
  // PetscCall(EPSSetType(eps,EPSARPACK));

  // Subspace Iteration with Rayleigh-Ritz projection and locking. *BAD*
  // PetscCall(EPSSetType(eps,EPSSUBSPACE));
    

  // PetscCall(EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE));

    PetscInt mpd = PETSC_DEFAULT; // max size of the sub-space
    int numEigenValues = numEigsToCompute;
    int & numArnoldiVectors = cgWave.dbase.get<int>("numArnoldiVectors");
  // PetscInt ncv = PETSC_DEFAULT; // numEigenValues*2+1; // size of column space 

    PetscInt ncv;
    if( numArnoldiVectors<=0 )
    {
        ncv= numEigenValues*2+1; // number of column vectors to used in the space, at least 2*numEigenValues
        numArnoldiVectors = ncv;
    }
    else
    {
        ncv=numArnoldiVectors;
    }

  // mpd = 100; // ** TEST ***
    PetscCall(EPSSetDimensions(eps,numEigenValues,ncv,mpd)); 

  // const int & maximumNumberOfIterations = cgWaveHoltz.dbase.get<int>("maximumNumberOfIterations");

    PetscInt maxIt = maximumNumberOfIterations; // 500;  
    printF("\n !!!!!!!!!!!!!!!   solveSLEPC: set maxIt = %d numEigenValues=%d, nc=%d!!!!!!!!!!!!!!!!\n",maxIt,numEigenValues,ncv);

  // Arnoldi iterations may be scaled by the number of request eigenvalues  
  //  nev=3 : maxit=50 --> 349 matVecs
  // if( eigenSolver==CgWave::ArnoldiEigenSolver )
  //   maxIt = 50;         // This is for some sort of outer iteration for Arnoldi, actually mat-vecs are 10 times or so

  // PetscScalar tol=1e-12;
    PetscCall(EPSSetTolerances(eps,tol,maxIt));  


  // Set the user defined monitor 
  // EPSMonitorSet(EPS eps,PetscErrorCode (*monitor)(EPS eps,PetscInt its,PetscInt nconv,PetscScalar *eigr,PetscScalar *eigi,PetscReal *errest,PetscInt nest,void *mctx),void *mctx,PetscErrorCode (*monitordestroy)(void**))

    EPSMonitorSet(eps,eigenWaveMonitor,NULL,NULL);

  // ----------------------------------------------------
  // -------------- INITIAL CONDITIONS ------------------
  // ----------------------------------------------------
    const int & initialVectorsForEigenSolver = cgWave.dbase.get<int>("initialVectorsForEigenSolver");
    const int maxInitialVectors=5;
    Vec vInit[maxInitialVectors]; 

    PetscInt numInitialConditions=0; 
    bool setInitialConditions= false;

    if( eigenSolverInitialCondition==CgWave::defaultEigenSolverInitialCondition )
    {
          printF("\n >>>> solveSLEPC : use default initial condition provided by SLEPc <<<<<<\n\n");
    }
    else if( eigenSolverInitialCondition==CgWave::randomEigenSolverInitialCondition ||
                      eigenSolverInitialCondition==CgWave::sineEigenSolverInitialCondition   ||
                      eigenSolverInitialCondition==CgWave::sumOfEigenvectorsInitialCondition )
    {
        if( eigenSolverInitialCondition==CgWave::randomEigenSolverInitialCondition )
            printF("\n >>>> solveSLEPC : use random initial condition <<<<< \n\n");
        else if( eigenSolverInitialCondition==CgWave::sineEigenSolverInitialCondition )
            printF("\n >>>> solveSLEPC : use sine initial condition <<<<< \n\n");
        else
            printF("\n >>>> solveSLEPC : use sum of eigenvector solutions initial condition <<<<< \n\n");


        setInitialConditions=true;
        numInitialConditions=1; // we supply one initial condition -- this is all Krylov Schur seems to take
    // initial conditions are stored here: 
        realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");

  
        const int ii=0; 
    // PetscCall(MatGetVecs(Amf,NULL,&vInit[ii]));  // create a Petsc vec    // *** DELETE ME *******************
    // For PETSc 3.6 and after use:
        PetscCall(MatCreateVecs(Amf,NULL,&vInit[ii]));  // create a Petsc vec    // *** DELETE ME *******************

        PetscScalar *vil;
        PetscCall(VecGetArray( vInit[ii],&vil ));  // get the local array from Petsc
        int iStart,iEnd;
        PetscCall(VecGetOwnershipRange(vInit[ii],&iStart,&iEnd));

        if( useMatrixUtilities )
        {
            cgWave.gridFunctionToVector( v,vil, iStart,iEnd );
        }
        else    
        {
            int count=0; // counts points 
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                const IntegerArray & gid = mg.gridIndexRange();      

                OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);

                    Iv[2]=Range(0,0);
                    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                    {
                        for( int side=0; side<=1; side++ )
                        {
                            int is = 1-2*side;
                            iab[side]=gid(side,axis);
                            const int bc = mg.boundaryCondition(side,axis);
                            if( bc==CgWave::dirichlet )
                            {
                                  iab[side] += is;  // Dirichlet BC -- ignore the boundary
                            }
                            else if( bc==CgWave::neumann || bc==CgWave::abcEM2  || bc==CgWave::absorbing || bc==CgWave::radiation )
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

                FOR_3D(i1,i2,i3,I1,I2,I3)
                {
                    if( maskLocal(i1,i2,i3) > 0 )
                    {
                        assert( count<iEnd );
                        vil[count]= vLocal(i1,i2,i3,0);   

                        count++;
                    }
                }
            }
            assert( count==iEnd );
        }

    // ------- Give SLEPc initial conditions --------------
        int numInitialConditions=1; 
        printF("solveSLEPc: supply %d initial conditions to SLEPc\n",numInitialConditions);
        EPSSetInitialSpace(eps, numInitialConditions, vInit ); 

    }
    else
    {
        
        printF("\n >>>> solveSLEPC : UNKNOWN initial condition option = %d ??? Using default.<<<<< \n\n",(int)eigenSolverInitialCondition);

    }


  // *************************************************************************************
  // ***************************** SOLVE FOR EIGENPAIRS **********************************
  // *************************************************************************************
    
    PetscCall(EPSSolve(eps));


    printF("\n !!!!!!!!!!!!!!!   solveSLEPC: set maxIt = %d !!!!!!!!!!!!!!!!\n",maxIt);

  // Optional: Get some information from the solver and display it
    EPSType        type;
    PetscInt       nev;
    PetscCall(EPSGetType(eps,&type));
    printF(">>Solution method: %s\n",type);
    PetscCall(EPSGetDimensions(eps,&nev,NULL,NULL));

    printF(">>Number of requested eigenvalues: %d\n",nev);
    PetscInt nconv;
    PetscCall(EPSGetConverged(eps,&nconv));
    printF(">>Number of converged eigenvalues= %d\n",nconv);

    PetscCall(EPSGetIterationNumber(eps,&its));
    printF(">>Number of iterations of the method: %d\n",its);

    int & numberOfMatrixVectorMultiplications = dbase.get<int>("numberOfMatrixVectorMultiplications");
    numberOfMatrixVectorMultiplications = iteration;

  // PetscInt lits;
  // PetscCall(KSPGetTotalIterations(ksp,&lits));
  // printF(">>Number of linear iterations of the method: %d\n",lits);


  // Create PETSc vectors that have the same distribution as Amf
    Vec  xr,xi;
  // PetscCall(MatGetVecs(Amf,NULL,&xr));                  // destroy me 
  // PetscCall(MatGetVecs(Amf,NULL,&xi));                  // destroy me 
  // For PETSc 3.6 and after use:  
    PetscCall(MatCreateVecs(Amf,NULL,&xr));                  // destroy me 
    PetscCall(MatCreateVecs(Amf,NULL,&xi));                  // destroy me 
  // PetscCall(EPSPrintSolution(eps,NULL));
    PetscScalar kr, ki;
  // for( int i=0; i<nconv; i++ )
  // {
  //   PetscCall(EPSGetEigenpair(eps,i,&kr,&ki,xr,xi)); 
  //   // PetscCall(EPSGetEigenvalue(eps,i,&kr,&ki)); 
  //   printF("i=%d: kr=%g, ki=%g\n",i,kr,ki);
  // }

    if( !dbase.has_key("numEigenVectors") )
    {
        dbase.put<int>("numEigenVectors");
        dbase.put<RealArray>("eigenValues");
    }
    int & numEigenVectors    = dbase.get<int>("numEigenVectors");
    RealArray & eigenValues  = dbase.get<RealArray>("eigenValues");

    numEigenVectors=nconv; // min(nev,nconv); // could use nconv


  // realCompositeGridFunction ucg(cg,all,all,all,numEigenVectors);
  // for( int i=0; i<numEigenVectors; i++ )
  //   ucg.setName(sPrintF("phi%d",i),i);

    if( !dbase.has_key("eigenVectorRayleighRitz") )
    {
        dbase.put<realCompositeGridFunction>("eigenVectorRayleighRitz");
    }

    realCompositeGridFunction & ucg = dbase.get<realCompositeGridFunction>("eigenVectorRayleighRitz");
    ucg.updateToMatchGrid(cg,all,all,all,numEigenVectors);
      for( int i=0; i<numEigenVectors; i++ )
        ucg.setName(sPrintF("phi%d",i),i);

    ucg=0.; // do this for now


    int ige;

    if( nconv>0 && nev>0 ) 
    {
    // eig.redim(2,numEigenVectors);

        for( int i=0; i<numEigenVectors; i++ ) 
        {

            PetscScalar kr, ki;
            PetscCall(EPSGetEigenpair(eps,i,&kr,&ki,xr,xi)); 
      // eig(0,i) = kr;
      // eig(1,i) = ki;

            printF("WaveHoltz eigenvalue beta %d : k=%18.14e + %18.14e I \n",i,kr,ki);

      // ierr = EPSGetEigenvector(eps,i,xr,xi);CHKERRQ(ierr);

      // ierr = VecView(xr,viewer);CHKERRQ(ierr);
      // if (!ishermitian) { ierr = VecView(xi,viewer);CHKERRQ(ierr); }

            if( i<numEigenVectors )
            {
        // ---- Save the eigenvector ----
        // printF("Save eigenvector %d to the grid function.\n",i);

                const int nComp = i;            // fill in this component


                PetscScalar *xrv;
                VecGetArray(xr,&xrv);  // get the local array from Petsc
                PetscCall(VecGetOwnershipRange(xr,&Istart,&Iend));

        // FOR NOW -- FIRST SAVE in vOld so we cn apply BC's .. then copy to ucg
                vOld=0.;

                const int numberOfComponents=1;
        // *new* way (from PETScSolver.bC)

                if( useMatrixUtilities )
                {
                    cgWave.vectorToGridFunction( xrv, vOld, Istart,Iend );
                }
                else        
                {
                    int ig=0; // *** DO THIS FOR NOW **********************
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {  
                        MappedGrid & mg = cg[grid];
                        const IntegerArray & gid = mg.gridIndexRange();

            // realArray & ug= ucg[grid];
            // realSerialArray uLocal; getLocalArrayWithGhostBoundaries(ug,uLocal);
            // realSerialArray vOldLocal; getLocalArrayWithGhostBoundaries(vOld,vOldLocal);
                        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
                        OV_GET_SERIAL_ARRAY(real,vOld[grid],vOldLocal);


            // int extra=-1;  // ******** DO THIS FOR NOW --- *********** FIX ME **********
            // getIndex(cg[grid].gridIndexRange(),I1,I2,I3,extra);   

                            Iv[2]=Range(0,0);
                            for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
                            {
                                for( int side=0; side<=1; side++ )
                                {
                                    int is = 1-2*side;
                                    iab[side]=gid(side,axis);
                                    const int bc = mg.boundaryCondition(side,axis);
                                    if( bc==CgWave::dirichlet )
                                    {
                                          iab[side] += is;  // Dirichlet BC -- ignore the boundary
                                    }
                                    else if( bc==CgWave::neumann || bc==CgWave::abcEM2  || bc==CgWave::absorbing || bc==CgWave::radiation )
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

            // int ig=getGlobalIndex( n, iv, grid, myid );  // get the global index for the first point

                        FOR_3D(i1,i2,i3,I1,I2,I3)
                        {

              // ******** NOTE we can probably just increment ig by 1 if we start correctly
              // int ig=getGlobalIndex( iv, grid, myid );  // get the global index
                            if( maskLocal(i1,i2,i3) > 0 )
                            {
                                if( ig>=Istart && ig<=Iend )
                                {
                                    if( false ) printf("SP:: myid=%i: i1,i2=%i,%i, ig=%i xrv[ig]=%6.4f\n",myid,i1,i2,ig,xrv[ig-Istart]);
                                    vOldLocal(i1,i2,i3)=xrv[ig-Istart];
                                }
                                else
                                {
                                    int p=myid;
                                    printf("SP::ERROR: myid=%i, i1,i2=%i,%i, ig=%i Istart,Iend=[%i,%i]\n", myid,i1,i2,ig,Istart,Iend);
                                }
                                ig++;
                            }
                        }
                        
                    }
                }

                if( false ) // Is this needed for pipe? maybe not 
                {
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        ucg[grid].periodicUpdate();
                        ucg[grid].updateGhostBoundaries();
                    }
                }

                if( 1==1 )
                { 
                    Real t=0.;
                    cgWave.applyEigenFunctionBoundaryConditions( vOld );

                    Real vNorm = maxNorm( vOld );

          // Now copy the result into the vector grid function ucg, and normalize
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        OV_GET_SERIAL_ARRAY(real,ucg[grid],ucgLocal);
                        OV_GET_SERIAL_ARRAY(real,vOld[grid],vOldLocal);

                        getIndex(cg[grid].dimension(),I1,I2,I3);  
                        int includeParallelGhost=1;
                        bool ok=ParallelUtility::getLocalArrayBounds(vOld[grid],vOldLocal,I1,I2,I3,includeParallelGhost);
                        if( ok )             
                            ucgLocal(I1,I2,I3,nComp) = vOldLocal(I1,I2,I3)*(1./vNorm);     
                    }


                }
                else
                {
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {  
                        ucg[grid].updateGhostBoundaries();
                        if( debug & 8 )
                            display(ucg[grid],sPrintF("Eigenvectors: ucg[%i]",grid),"%6.3f ");
                    }          
                }
            }
        }

        if( false )
        {
            fflush(0);
            PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
            PetscCall(EPSConvergedReasonView(eps,PETSC_VIEWER_STDOUT_WORLD));
            PetscCall(EPSErrorView(eps,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
            PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
        }

        const CgWave::TimeSteppingMethodEnum & timeSteppingMethod = cgWave.dbase.get<CgWave::TimeSteppingMethodEnum>("timeSteppingMethod");
        int & minStepsPerPeriod              = cgWave.dbase.get<int>("minStepsPerPeriod");
        aString & nameOfGridFile             = cgWave.dbase.get<aString>("nameOfGridFile");
        printF("\n");
        printF(" -----------------------------------------------------------------------------\n");
        printF(" ---------------------- EigenWave SUMMARY (SLEPc solver) ---------------------\n");
        printF(" omega=%12.4e, grid=%s\n",omega,(const char*)nameOfGridFile);
        printF(" eigenSolver = %s\n",(eigenSolver==CgWave::defaultEigenSolver ? "KrylovSchur" :
                                                                    eigenSolver==CgWave::KrylovSchurEigenSolver ? "KrylovSchur" :
                                                                    eigenSolver==CgWave::ArnoldiEigenSolver ? "Arnoldi" : 
                                                                    eigenSolver==CgWave::ArpackEigenSolver ? "Arpack" : 
                                                                    "unknown" ) );
        printF(" provide initial vector(s) for eigenSolver = %d\n",initialVectorsForEigenSolver);
        printF(" timeSteppingMethod = %s, minStepsPerPeriod=%d\n",
                                                                                (timeSteppingMethod==CgWave::explicitTimeStepping ? "explicit (modified equation)" :
                                                                                  timeSteppingMethod==CgWave::implicitTimeStepping ? "implicit" : "unknown"),
                                                                                  minStepsPerPeriod );    
        printF(" tol=%9.2e, numIterations=%d, provide initial vectors=%d\n",tol,iteration,(int)setInitialConditions);
        printF(" num eigs requested=%d, number eigs converged=%d, numArnoldiVectors=%d\n",nev,nconv,numArnoldiVectors);
        printF(" -----------------------------------------------------------------------------\n");
    // -- compute eigenvalues of Lh from the Rayleigh Quotient ----
        eigenValues.redim(numEigenVectors);
        for( int i=0; i<numEigenVectors; i++ ) 
        {
            Real lamRQ = cgWave.getRayleighQuotient( ucg, i );
            eigenValues(i) = lamRQ;

      // if( true )
      // {
      //   realCompositeGridFunction & res = dbase.get<realCompositeGridFunction>("residual");
      //   res.updateToMatchGrid(cg,all,all,all,numEigenVectors); 

      //   Real maxRes = cgWave.getEigenPairResidual( lamRQ, ucg, res, i );
      //   printF(" i=%d: lamRQ = %16.10e,  rel-resid = || L v + lamRQ^2 v ||/lamRQ^2 = %9.2e\n",i,lamRQ,maxRes);
      // }
        }

    // ------------------- SORT -------------------------

        bool sortEigenPairs=true;
        if( sortEigenPairs )
        {
      // -- sort the eigenvalues:
            IntegerArray iperm(numEigenVectors);
            cgWave.sortArray( eigenValues, iperm );

            if( 1==1 || debug > 0 )
            {
                fprintf(pDebugFile,"debug=%d\n",debug);
                ::display(eigenValues,"sorted,eigenvalues",pDebugFile,"%20.13e ");
                ::display(iperm,"iperm : sort eigenpairs",pDebugFile,"%3d");
                fflush(pDebugFile);
            }

      // -- reorder the eigenvectors -----
            int numAssigns=0; 
            for( int i=0; i<numEigenVectors; i++ )
            {
        // Delay this copy until it is actually needed
                bool entrySaved=false; 
        // Real x = a(i);
        // copyGridFunction( vOld,0, ucg,i ); numAssigns++;

                int j=i;
                for( int ii=0; ii<numEigenVectors; ii++ )
                {
                    int k = iperm(j);
                    iperm(j)=j;
                    if( k==i )
                        break;
                    
                    if( j!=k )
                    {
                        if( !entrySaved )
                        { // now we need to save ucg[i]
                                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                                {
                                    OV_GET_SERIAL_ARRAY(real,vOld[grid],uLocal);
                                    OV_GET_SERIAL_ARRAY(real,ucg[grid],vLocal);  // temp space 
                                    getIndex(cg[grid].dimension(),I1,I2,I3);
                                    int includeParallelGhost=1;  
                                    bool ok=ParallelUtility::getLocalArrayBounds(ucg[grid],vLocal,I1,I2,I3,includeParallelGhost);
                                    if( ok )      
                                        uLocal(I1,I2,I3,0)= vLocal(I1,I2,I3,i);
                                }
                            numAssigns++; entrySaved=true;
                        }
            // a(j)=a(k); 
                            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                            {
                                OV_GET_SERIAL_ARRAY(real,ucg[grid],uLocal);
                                OV_GET_SERIAL_ARRAY(real,ucg[grid],vLocal);  // temp space 
                                getIndex(cg[grid].dimension(),I1,I2,I3);
                                int includeParallelGhost=1;  
                                bool ok=ParallelUtility::getLocalArrayBounds(ucg[grid],vLocal,I1,I2,I3,includeParallelGhost);
                                if( ok )      
                                    uLocal(I1,I2,I3,j)= vLocal(I1,I2,I3,k);
                            }
                        numAssigns++;
                    }
                    j=k;
                }
                if( i!=j )
                {
          // a(j)=x; 
                    assert( entrySaved );
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            OV_GET_SERIAL_ARRAY(real,ucg[grid],uLocal);
                            OV_GET_SERIAL_ARRAY(real,vOld[grid],vLocal);  // temp space 
                            getIndex(cg[grid].dimension(),I1,I2,I3);
                            int includeParallelGhost=1;  
                            bool ok=ParallelUtility::getLocalArrayBounds(vOld[grid],vLocal,I1,I2,I3,includeParallelGhost);
                            if( ok )      
                                uLocal(I1,I2,I3,j)= vLocal(I1,I2,I3,0);
                        }
                    numAssigns++;
                }
            }

              printF(">>> Reorder eigenVectors: number of assignments to reorder=%d\n",numAssigns);

      // // To reorder EV's we need another permuation matrix
      // IntegerArray jperm(numEigenVectors);
      // for( int ie=0; ie<numEigenVectors; ie++ )
      //   jperm(ie)=ie;


      // reorder an array according to given indexes

      // // re-order the eigenvectors 
      // for( int i=0; i<numEigenVectors; i++ )
      // {
      //   x = a[i];
      //   int j = i;
      //   for( ii=0; ii<numEigenVectors; ii++ )
      //   {
      //     int k = iperm(j);
      //     iperm(j)=j;
      //     if( k==i )
      //       break;
      //     a[j] = a[k];
      //     j=k;
      //   }
      //   a[j] = x;
      // }

      // for( int ie=0; ie<numEigenVectors; ie++ )
      // {
      //   int ip = iperm(ie);  // this element goes next;

      //   int je=  jperm( ip ); // this is where it is located 
      //   printF(" ie=%d: next element ip=iperm(ie)=%d is located here = jperm(ie) ) = %d\n",ie,iperm(ie), je);
      //   if( ie!=je )
      //   {
      //     int jtemp = jperm(ie); jperm(ie)=jperm(je); jperm(je)=jtemp;  // swap jperm
      //     for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      //     {
      //       OV_GET_SERIAL_ARRAY(real,ucg[grid],ucgLocal);
      //       OV_GET_SERIAL_ARRAY(real,vOld[grid],vOldLocal);  // temp space 
      //       getIndex(cg[grid].dimension(),I1,I2,I3);  

      //       vOldLocal(I1,I2,I3)  = ucgLocal(I1,I2,I3,ie); 
      //       ucgLocal(I1,I2,I3,ie)= ucgLocal(I1,I2,I3,je);
      //       ucgLocal(I1,I2,I3,je)= vOldLocal(I1,I2,I3);
      //     }
      //   }

      // }
      // OV_ABORT("stop here for now");
        }


        bool computeResiduals=true;
        if( computeResiduals )
        {
            for( int i=0; i<numEigenVectors; i++ ) 
            {
                realCompositeGridFunction & res = dbase.get<realCompositeGridFunction>("residual");
                res.updateToMatchGrid(cg,all,all,all,numEigenVectors); 

                Real lamRQ = eigenValues(i);
                Real maxRes = cgWave.getEigenPairResidual( lamRQ, ucg, res, i );
                printF(" i=%d: lamRQ = %16.10e,  rel-resid = || L v + lamRQ^2 v ||/lamRQ^2 = %9.2e\n",i,lamRQ,maxRes);
            }
        }    
    // ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);




    } // end nconv >0 ...

  // GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("test"));
  // ps.erase();
  // PlotStuffParameters psp;
  // psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
  // PlotIt::contour(ps,ucg,psp);
  // psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true); 



    computeEigenmodesWithSLEPc=0; // reset 
    adjustHelmholtzForUpwinding = adjustHelmholtzForUpwindingSave; // reset

  // *** DESTROY STUFF ****
    PetscCall(EPSDestroy(&eps));
    PetscCall(MatDestroy(&Amf));

  // PetscCall(VecDestroy(&u));
  // PetscCall(VecDestroy(&b));
  // PetscCall(VecDestroy(&x)); 

    PetscCall(VecDestroy(&xr)); 
    PetscCall(VecDestroy(&xi)); 
    for( int ii=0; ii<numInitialConditions; ii++ )
    {
          PetscCall(VecDestroy(&vInit[ii]));
    } 
    

  // OV_ABORT("petscSolver: computeEigenmodes - STOP HERE FOR NOW");


    return 0;
}
