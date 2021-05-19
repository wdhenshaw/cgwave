#include "CgWaveHoltz.h" 
#include "CgWave.h" 
#include "Overture.h"
#include "gridFunctionNorms.h"
#include "display.h"

// krb do not use extern "C" if PETSc is linked using BOPT=g_c++
extern "C"
{
// *wdh* 2015/09/31  To avoid having PETSc include complex.h do this: 
#include "petscconf.h"
#undef PETSC_HAVE_CXX_COMPLEX
#include "petscksp.h"
}

static char help[] = "CgWaveHoltz test of PETSc\n";

// testCase : 
// 0 = test solving laplace equation
// 1 = real run
static int testCase =1; 
static int iteration=0;

static CgWaveHoltz *pCgWaveHoltz; // pointer to the CgWaveHoltz solver 
 

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

// -- global variables -- do this for now 
Mat *pA=NULL;
Vec *pb=NULL;
int mA,nA;

// =========================================================================================================
//     MATRIX-VECTOR MULTIPLY FOR MATRIX FREE KRYLOV SOLVERS
// 
// Compute y = M*x 
// =========================================================================================================
extern PetscErrorCode usermult(Mat m ,Vec x, Vec y)
{
  PetscErrorCode ierr = 0;

  printf("+++ usermult called iteration=%d\n",iteration);
  if( testCase==0 )
  {
    // solve a Poisson equation **TEST CASE***

    Mat *matrix;
    MatShellGetContext(m, (void**)&matrix);
    // printf(" matrix=%p\n",matrix);
    // printf("     pA=%p\n",pA);
  
    PetscScalar *xl, *yl;
    VecGetArray(x,&xl);  // get the local array from Petsc
    VecGetArray(y,&yl);  // get the local array from Petsc
    int iStart,iEnd;
    ierr = VecGetOwnershipRange(x,&iStart,&iEnd);CHKERRQ(ierr);
    printf("  [iStart,iEnd]=[%i,%i] mA=%d nA=%d\n",iStart,iEnd,mA,nA);

    // column major ordering to match PETSc code 
    #define v(i,j) xl[j + nA*(i)]
    #define w(i,j) yl[j + nA*(i)]

    if( iteration<=1 )
    {
      printf("iteration=%d: MF : x=[",iteration);
      for( int i=iStart; i<iEnd; i++ ){ printf("%5.2f,",xl[i]); } // 
      printf("];\n");
    }


    for( int j=0; j<nA; j++ )
      for( int i=0; i<mA; i++ )
      {
        /// *** FIX THIS FOR PARALLEL ***
        int k = j + nA*(i);
        assert( k>=iStart && k<iEnd );   
    
        double vij = v(i,j);
        double vimj =  i>0    ? v(i-1,j) : 0.;
        double vipj =  i<mA-1 ? v(i+1,j) : 0.;
        double vijm =  j>0    ? v(i,j-1) : 0.;
        double vijp =  j<nA-1 ? v(i,j+1) : 0.;

        // w(i,j) =  4.*v(i,j) -( v(i+1,j) + v(i-1,j) + v(i,j-1) + v(i,j+1) );
        w(i,j) =  4.*vij -( vipj + vimj + vijm + vijp );
    
      }
    if( iteration<=1 )
    {
      printf("iteration=%d: MF : y=[",iteration);
      for( int i=iStart; i<iEnd; i++ ){ printf("%5.2f,",yl[i]); } // 
      printf("];\n");
    }
  
    if( 0==1 )
    {
      printf("+++ usermult call MatMult...\n");
      ierr = MatMult(*matrix, x, y);

      printf("A*x: y=[");
      for( int i=iStart; i<iEnd; i++ ){ printf("%5.2f,",yl[i]); } // 
      printf("];\n");

    }
  
    // ierr = MatMult(*pA, x, y);
    printf("... done\n");

    // assert( 1==0 );
  }
  else
  {
    // ---- Helmholtz solver ----

    printF("\n **************** MatVec for PETSc iteration=%i *************\n\n",iteration);

    assert( pCgWaveHoltz!=NULL );
    CgWaveHoltz & cgWaveHoltz = *pCgWaveHoltz;
    
    const int & monitorResiduals = cgWaveHoltz.dbase.get<int>("monitorResiduals");      // montior the residuals at every step

    // here is the CgWave solver for the time dependent wave equation
    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");

    // -- cgWave solution is stored here: 
    realCompositeGridFunction & v    = cgWave.dbase.get<realCompositeGridFunction>("v");

    realCompositeGridFunction & vOld = cgWaveHoltz.dbase.get<realCompositeGridFunction>("vOld");
    CompositeGrid & cg = cgWaveHoltz.cg;

    if( !cgWaveHoltz.dbase.has_key("bcg") )
    {
      realCompositeGridFunction & bcg = cgWaveHoltz.dbase.put<realCompositeGridFunction>("bcg");
      bcg.updateToMatchGrid(cg);
    }
    realCompositeGridFunction & bcg = cgWaveHoltz.dbase.get<realCompositeGridFunction>("bcg");    


    PetscScalar *xl, *yl;
    VecGetArray(x,&xl);  // get the local array from Petsc
    VecGetArray(y,&yl);  // get the local array from Petsc
    int iStart,iEnd;
    ierr = VecGetOwnershipRange(x,&iStart,&iEnd);CHKERRQ(ierr);
    // printf("  [iStart,iEnd]=[%i,%i]\n",iStart,iEnd);


    // --- Set v = x ---
    v=0.;

    Index I1,I2,I3;
    int i=0;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      getIndex(cg[grid].dimension(),I1,I2,I3);

      OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
      FOR_3D(i1,i2,i3,I1,I2,I3)
      {
        assert( i<iEnd );
        vLocal(i1,i2,i3)=xl[i];
        i++;
      }
    }
    assert( i==iEnd );
      
    if( false && iteration==1 )
    {
      ::display(v[0],"v (iteration 1)","%5.2f ");
      
    }
    
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      vOld[grid] = v[grid];  // save current guess
  
    // -- advance for one period (or multiple periods ) ---
    cgWave.advance( iteration );
    

    if( iteration==0 )
    {
      // ---------- INITIALIZE ITERATIONS -----

      printf("COMPUTE b = P*v(0) \n");
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        bcg[grid]= v[grid];

      if( false )
      {
        ::display(bcg[0],"b: iteration 0 ","%5.2f ");
      }

      //compute y = P*x 
      // set y = v 
      int i=0;
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        getIndex(cg[grid].dimension(),I1,I2,I3);
        OV_GET_SERIAL_ARRAY(int,cg[grid].mask(),maskLocal);

        OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
          assert( i<iEnd );
          if( true || maskLocal(i1,i2,i3)>0 )
          {
            yl[i]= vLocal(i1,i2,i3);
          }
          else
          {
            // un-used points: set [Ax]_i = x_i   -- is this right ??
            yl[i]=xl[i];
          }
          
          i++;
        }
      }
      assert( i==iEnd );
    }
    else if( iteration==-1 )
    {
      // ---- iteration=-1: (after final step usually) check the residual in the computed solution  : v^n - v^{n+1}
      // Compute :
      //       A v^n = v^n - v^{n+1} + b 

      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        vOld[grid] -= v[grid]; // 
    
      const real & tol = cgWaveHoltz.dbase.get<real>("tol");

      real errMax = maxNorm(vOld);
      printF("it=%d:  max(|v-vOld|)=%8.2e, tol=%g\n",iteration,errMax,tol);

    }
    else
    {
      // ----------- MATRIX VECTOR MULTIPLY FOR Krylov solver -------
      // Compute :
      //       A v^n = v^n - v^{n+1} + b 
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        vOld[grid] -= v[grid]; // 


      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        vOld[grid] += bcg[grid];

      // should we interpolate vOld?
      vOld.interpolate();

      assert( pb !=NULL );
      Vec & b = *pb;
      PetscScalar *bl;
      VecGetArray(b,&bl);  // get the local array from Petsc

      // set y = v 
      int i=0;
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        getIndex(cg[grid].dimension(),I1,I2,I3);
        OV_GET_SERIAL_ARRAY(int,cg[grid].mask(),maskLocal);
        OV_GET_SERIAL_ARRAY(real,vOld[grid],vOldLocal);
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
          assert( i<iEnd );
          // yl[i]= vOldLocal(i1,i2,i3,0) + bl[i]; 
          if( true || maskLocal(i1,i2,i3)>0 )          // ************************* check me ***************
          {
            yl[i]= vOldLocal(i1,i2,i3,0);
          }
          else
          {
            yl[i]=xl[i];
          }
          
          i++;
        }
      }
      assert( i==iEnd );

      // ----- Optionally save the current residual -----
      if( monitorResiduals )
      {
        if( iteration<0 || iteration>=resVector.getLength(0) )
        {
          printF("solvePETSC: ERROR: INVALID iteration=%d. resVector.getLength(0)=%d\n",iteration,resVector.getLength(0));
          OV_ABORT("ERROR");
        }
        resVector(iteration)= cgWaveHoltz.residual(); 
      } 
    
    }


  }

  iteration++;
  
  return ierr;
}


// ============================================================================
/// \brief Solve for the Helholtz solution using PETSc
// ============================================================================
int CgWaveHoltz::
solvePETSc(int argc,char **args)
{
  pCgWaveHoltz=this;

  
  const real & omega                    = dbase.get<real>("omega");
  real & Tperiod                        = dbase.get<real>("Tperiod");
  const int & numPeriods                = dbase.get<int>("numPeriods");
  const int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");

 // here is the CgWave solver for the time dependent wave equation
 CgWave & cgWave = *dbase.get<CgWave*>("cgWave");
 Tperiod=numPeriods*twoPi/omega;  
 printF("CgWaveHoltz::solvePETSc: setting tFinal = Tperiod*numPeriods = %9.3e (numPeriods=%d) \n",Tperiod,numPeriods);
 
 cgWave.dbase.get<real>("omega")=omega; // ** FIX ME **
 cgWave.dbase.get<real>("tFinal")=Tperiod; // ** FIX ME **
 cgWave.dbase.get<real>("Tperiod")=Tperiod; // ** FIX ME **
 cgWave.dbase.get<int>("numPeriods")=numPeriods; // ** FIX ME **

   // Save "residuals" by iteration: 
  // resVector(it) = norm( v^{n+1} - v^n )
  RealArray & resVector = dbase.get<RealArray>("resVector");
  resVector.redim(maximumNumberOfIterations);
  resVector=0.;
  

  Vec            x,b,u;  /* approx solution, RHS, exact solution */
  Mat            A;        /* linear system matrix */
  KSP            ksp;     /* linear solver context */
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
    PetscInitialize(&argc,&args,(char *)0,help);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
  }
  
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Compute the matrix and right-hand-side vector that define
     the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* 
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good 
     performance. See the matrix chapter of the users manual for details.
  */
  if( false )
  {

    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,m*n,m*n);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A,5,PETSC_NULL,5,PETSC_NULL);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(A,5,PETSC_NULL);CHKERRQ(ierr);
    ierr = MatSetUp(A);CHKERRQ(ierr);
  
    /* 
       Currently, all PETSc parallel matrix formats are partitioned by
       contiguous chunks of rows across the processors.  Determine which
       rows of the matrix are locally owned. 
    */
    ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

    /* 
       Set matrix elements for the 2-D, five-point stencil in parallel.
       - Each processor needs to insert only elements that it owns
       locally (but any non-local elements will be sent to the
       appropriate processor during matrix assembly). 
       - Always specify global rows and columns of matrix entries.

       Note: this uses the less common natural ordering that orders first
       all the unknowns for x = h then for x = 2h etc; Hence you see J = Ii +- n
       instead of J = I +- m as you might expect. The more standard ordering
       would first do all variables for y = h, then y = 2h etc.

    */
    ierr = PetscLogStageRegister("Assembly", &stage);CHKERRQ(ierr);
    ierr = PetscLogStagePush(stage);CHKERRQ(ierr);
    for (Ii=Istart; Ii<Iend; Ii++) { 
      v = -1.0; i = Ii/n; j = Ii - i*n;  
      if (i>0)   {J = Ii - n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
      if (i<m-1) {J = Ii + n; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
      if (j>0)   {J = Ii - 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
      if (j<n-1) {J = Ii + 1; ierr = MatSetValues(A,1,&Ii,1,&J,&v,INSERT_VALUES);CHKERRQ(ierr);}
      v = 4.0; ierr = MatSetValues(A,1,&Ii,1,&Ii,&v,INSERT_VALUES);CHKERRQ(ierr);
    }

    /* 
       Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd()
       Computations can be done while messages are in transition
       by placing code between these two statements.
    */
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = PetscLogStagePop();CHKERRQ(ierr);

    /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
    ierr = MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
  }
  

  // ============== Create a shell object for matrix free method ================
  // double mycontext; // passed to usermult 
  Mat *mycontext=&A;
  
  // ******  global variables
  pA = &A; mA=m; nA=n;
  
  int numEquations=0;
  if( testCase==0 )
  {
    numEquations=m*n;
  }
  else
  {
    assert( pCgWaveHoltz!=NULL );
    CgWaveHoltz & cgWaveHoltz = *pCgWaveHoltz;
    CompositeGrid & cg = cgWaveHoltz.cg;
    Index I1,I2,I3;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      getIndex(cg[grid].dimension(),I1,I2,I3);
      numEquations += I1.getLength()*I2.getLength()*I3.getLength();
    }
  }
  printF("Make a Matrix Free Shell: numEquations=%d\n",numEquations);

  
  Mat Amf;
  ierr = MatCreateShell(PETSC_COMM_WORLD, numEquations, numEquations, PETSC_DECIDE,  PETSC_DECIDE, mycontext, &Amf);
  ierr = MatShellSetOperation(Amf, MATOP_MULT, (void(*)(void))usermult);
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
  ierr = VecSetSizes(u,PETSC_DECIDE,numEquations);CHKERRQ(ierr);
  ierr = VecSetFromOptions(u);CHKERRQ(ierr);
  ierr = VecDuplicate(u,&b);CHKERRQ(ierr); 
  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);

  pb = &b;  // global variable for now
  

  /* 
     Set exact solution; then compute right-hand-side vector.
     By default we use an exact solution of a vector with all
     elements of 1.0;  Alternatively, using the runtime option
     -random_sol forms a solution vector with random components.
  */
  if( testCase==0 )
  {
    ierr = PetscOptionsGetBool(PETSC_NULL,"-random_exact_sol",&flg,PETSC_NULL);CHKERRQ(ierr);
    if (flg) {
      ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);
      ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);
      ierr = VecSetRandom(u,rctx);CHKERRQ(ierr);
      ierr = PetscRandomDestroy(&rctx);CHKERRQ(ierr);
    } else {
      ierr = VecSet(u,1.0);CHKERRQ(ierr);
    }
    ierr = MatMult(Amf,u,b);CHKERRQ(ierr);

  }
  else
  {
    // Initial Guess is zero:
    ierr = VecSet(x,0.0);CHKERRQ(ierr);
    ierr = VecSet(u,0.0);CHKERRQ(ierr);
    // ierr = VecSet(b,0.0);CHKERRQ(ierr);
    // set RHS: 
    ierr = MatMult(Amf,u,b);CHKERRQ(ierr);
  }
  
  /*
     View the exact solution vector if desired
  */
  flg  = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-view_exact_sol",&flg,PETSC_NULL);CHKERRQ(ierr);
  if (flg) {ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);}

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /* 
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /* 
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
//   ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,Amf,Amf,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);


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
  const real & tol = dbase.get<real>("tol");
  
  ierr = KSPSetTolerances(ksp,tol,1.e-50,PETSC_DEFAULT,
                          PETSC_DEFAULT);CHKERRQ(ierr);

  /* 
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  if( testCase == 1)
  {
    printF("***** DONE KSP SOLVE ******\n");

    if( false )
    {
      // --- Set v = x ---
      PetscScalar *xl;
      VecGetArray(x,&xl);  // get the local array from Petsc
      int iStart,iEnd;
      ierr = VecGetOwnershipRange(x,&iStart,&iEnd);CHKERRQ(ierr);
      // printf("  [iStart,iEnd]=[%i,%i]\n",iStart,iEnd);
    
      realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");
      Index I1,I2,I3;
      int i=0;
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        getIndex(cg[grid].dimension(),I1,I2,I3);

        OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
        FOR_3D(i1,i2,i3,I1,I2,I3)
        {
          assert( i<iEnd );
          vLocal(i1,i2,i3)=xl[i];
          i++;
        }
      }
      assert( i==iEnd );      

      v.interpolate();
      
    }
    else
    {

      printF("***** CALL AGAIN TO CHECK SOLUTION ******\n");
      int & numberOfIterations = dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken
      numberOfIterations = iteration;


      iteration=-1;  // **fix me**
      ierr = MatMult(Amf,x,b);CHKERRQ(ierr);     


      real maxKspResidal;
      KSPGetResidualNorm(ksp,&maxKspResidal);

      printF("\n ################ DONE KYRLOV ITERATIONS -- KSP residual=%8.2e (tol=%8.2e) numberOfIterations=%i #############\n",
             maxKspResidal,tol,numberOfIterations);

      // real maxRes = residual();
      // printF("Maximum residual = %9.3e\n",maxRes);
    }
    
  }
  else
  {

    /* 
       Check the error
    */
    ierr = VecAXPY(x,-1.0,u);CHKERRQ(ierr);
    ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
    /* Scale the norm */
    /*  norm *= sqrt(1.0/((m+1)*(n+1))); */

    /*
      Print convergence information.  PetscPrintf() produces a single 
      print statement from all processes that share a communicator.
      An alternative is PetscFPrintf(), which prints to a file.
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %G iterations %D\n",
                       norm,its);CHKERRQ(ierr);
  }
  
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
  
  convergenceRate          = pow( resVector(numberOfIterations-1)/resVector(0), 1./( numberOfIterations ) ); 
  convergenceRatePerPeriod = pow( resVector(numberOfIterations-1)/resVector(0), 1./( numberOfIterations*numPeriods) ); 

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
