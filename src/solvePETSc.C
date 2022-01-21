#include "CgWaveHoltz.h" 
#include "CgWave.h" 
#include "Overture.h"
#include "gridFunctionNorms.h"
#include "display.h"
#include "ParallelUtility.h"

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
static KSP  ksp=NULL;     /* linear solver context */
static Mat *pA=NULL;
static Vec *pb=NULL;
static int mA,nA;

static int computeRightHandSide =-2; // = 0;
static int computeResidual      =-3; // =-1;

// *wdh* Jan 5, 2022 --> changed cgWave to zero out unused points => thus we do not need to check the mask
static bool checkMask = false;  // if true, do not use values of the solution where mask==0 

// =========================================================================================================
//     MATRIX-VECTOR MULTIPLY FOR MATRIX FREE KRYLOV SOLVERS
// 
// Compute y = M*x 
// =========================================================================================================
extern PetscErrorCode waveHoltzMatrixVectorMultiply(Mat m ,Vec x, Vec y)
{
  PetscErrorCode ierr = 0;

  printF("\n ++++++++ WaveHoltz Matrix vector multiply routine: waveHoltzMatrixVectorMultiply called iteration=%d\n",iteration);

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
      printF("iteration=%d: MF : x=[",iteration);
      for( int i=iStart; i<iEnd; i++ ){ printF("%5.2f,",xl[i]); } // 
      printF("];\n");
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
      printf("+++ waveHoltzMatrixVectorMultiply call MatMult...\n");
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

    printF(" **************** MatVec for PETSc iteration=%i *************\n\n",iteration);
    if( iteration==computeRightHandSide )
    {
      printF(" >>>> Call cgWave to COMPUTE the Right-hand-side to Ax=b : b = Pi * v(0) \n");    
    }

    assert( pCgWaveHoltz!=NULL );
    CgWaveHoltz & cgWaveHoltz = *pCgWaveHoltz;
    
    const int & monitorResiduals      = cgWaveHoltz.dbase.get<int>("monitorResiduals");      // montior the residuals at every step
    const Real & numberOfActivePoints = cgWaveHoltz.dbase.get<Real>("numberOfActivePoints");
   

    // here is the CgWave solver for the time dependent wave equation
    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
    const int & numberOfFrequencies = cgWave.dbase.get<int>("numberOfFrequencies");

    // -- cgWave solution is stored here: 
    realCompositeGridFunction & v    = cgWave.dbase.get<realCompositeGridFunction>("v");

    realCompositeGridFunction & vOld = cgWave.dbase.get<realCompositeGridFunction>("vOld");
    CompositeGrid & cg = cgWaveHoltz.cg;

    if( !cgWaveHoltz.dbase.has_key("bcg") )
    {
      realCompositeGridFunction & bcg = cgWaveHoltz.dbase.put<realCompositeGridFunction>("bcg");
      Range all;
      bcg.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);
    }
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
    int i=0;
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        MappedGrid & mg = cg[grid];
        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

        getIndex(cg[grid].dimension(),I1,I2,I3);

        OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
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
      
    if( false && iteration==1 )
    {
      ::display(v[0],"v (iteration 1)","%5.2f ");
      
    }
    
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
      OV_GET_SERIAL_ARRAY(real,vOld[grid],vOldLocal);
      vOldLocal = vLocal;  // save current guess

      // vOld[grid] = v[grid];  // save current guess
    }
  
    // -- advance for one period (or multiple periods ) ---
    cgWave.advance( iteration );
    

    if( iteration==computeRightHandSide )
    {
      // ---------- INITIALIZE ITERATIONS -----

      printF("COMPUTE b = Pi * v(0) \n");

      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
        OV_GET_SERIAL_ARRAY(real,bcg[grid],bcgLocal);

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
      int i=0;
      for( int freq=0; freq<numberOfFrequencies; freq++ )
      {      
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
          getIndex(cg[grid].dimension(),I1,I2,I3);
          OV_GET_SERIAL_ARRAY(int,cg[grid].mask(),maskLocal);

          OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
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
    else if( iteration==computeResidual )
    {
      // ---- iteration=computeResidual: (after final step usually) check the residual in the computed solution  : v^n - v^{n+1}
      // Compute :
      //       A v^n = v^n - v^{n+1} + b 

      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        vOld[grid] -= v[grid]; // 
    
      const real & tol = cgWaveHoltz.dbase.get<real>("tol");
      real & maxResidual = cgWaveHoltz.dbase.get<real>("maxResidual");

      maxResidual = maxNorm(vOld);
      printF("it=%d:  max(residual) = max(|v-vOld|)=%8.2e, tol=%g\n",iteration,maxResidual,tol);

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
      if( false ) //  Jan 5, 2022 : change to false : FIXES SIC
      {
        vOld.interpolate();
      }

      assert( pb !=NULL );
      Vec & b = *pb;
      PetscScalar *bl;
      VecGetArray(b,&bl);  // get the local array from Petsc

      // set y = v 
      int i=0;
      for( int freq=0; freq<numberOfFrequencies; freq++ )
      {       
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
          getIndex(cg[grid].dimension(),I1,I2,I3);
          OV_GET_SERIAL_ARRAY(int,cg[grid].mask(),maskLocal);
          OV_GET_SERIAL_ARRAY(real,vOld[grid],vOldLocal);
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

    if( monitorResiduals && iteration>=0 )
    {
      // ----- Optionally save the current residual -----
      // Save "residuals" by iteration: 
      // resVector(it) = norm( v^{n+1} - v^n )
      RealArray & resVector = cgWaveHoltz.dbase.get<RealArray>("resVector");  

      if( iteration<0 || iteration>=resVector.getLength(0) )
      {
        printF("solvePETSC: ERROR: INVALID iteration=%d. resVector.getLength(0)=%d\n",iteration,resVector.getLength(0));
        OV_ABORT("ERROR");
      }

      // resVector(iteration)= cgWaveHoltz.residual();   // this is not correct since solution is not the right one

      Real kspResidual;
      assert( ksp !=NULL );
      KSPGetResidualNorm(ksp,&kspResidual); 
      kspResidual /= sqrt(numberOfActivePoints);  // make an approximate L2h norm

      const int & maximumNumberOfIterations = cgWaveHoltz.dbase.get<int>("maximumNumberOfIterations");
      if( iteration <= maximumNumberOfIterations )
      {
        resVector(iteration)= kspResidual;
        printF("\n ##### SAVE KRYLOV RESIDUAL: iteration=%d: L2h-residual=%9.2e \n\n",kspResidual);
      }
      else
      {
         printF("\n ##### KRYLOV RESIDUAL: iteration=%d: L2h-residual=%9.2e (NOT SAVED since iteration=%d > maximumNumberOfIterations=%d)\n\n",
          kspResidual,iteration,maximumNumberOfIterations);
      }

      // // There is no residual for iteration=0 since the Krylov solver is just computing the first A*x
      // if( iteration==1 )
      //   resVector(0) = resVector(iteration); // do this for now -- ksp is not built yet for iteration=0
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
  const int & adjustOmega               = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t   
  const int & maximumNumberOfIterations = dbase.get<int>("maximumNumberOfIterations");
  Real & numberOfActivePoints           = dbase.get<Real>("numberOfActivePoints");

  // here is the CgWave solver for the time dependent wave equation
  CgWave & cgWave = *dbase.get<CgWave*>("cgWave");
  if( omega!=0. )
    Tperiod=numPeriods*twoPi/omega; 
  else 
    Tperiod=1.;
 
  printF("CgWaveHoltz::solvePETSc: setting tFinal = Tperiod*numPeriods = %9.3e (numPeriods=%d) \n",Tperiod,numPeriods);
 
  // --- set values in CgWave:  *** COULD DO BETTER ***
 
  cgWave.dbase.get<real>("omega")     = omega;        // ** FIX ME **
  cgWave.dbase.get<real>("tFinal")    = Tperiod;      // ** FIX ME **
  cgWave.dbase.get<real>("Tperiod")   = Tperiod;      // ** FIX ME **
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
  vOld.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);

  // <<<<


  // Save "residuals" by iteration: 
  // resVector(it) = norm( v^{n+1} - v^n )
  RealArray & resVector = dbase.get<RealArray>("resVector");
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
    PetscInitialize(&argc,&args,(char *)0,help);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);
  }
  

  // ============== Create a shell object for matrix free method ================
  // double mycontext; // passed to waveHoltzMatrixVectorMultiply 
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
    numEquations *= numberOfFrequencies;
  }
  printF("Make a Matrix Free Shell: numEquations=%d\n",numEquations);

  
  Mat Amf;
  ierr = MatCreateShell(PETSC_COMM_WORLD, numEquations, numEquations, PETSC_DECIDE,  PETSC_DECIDE, mycontext, &Amf);
  ierr = MatShellSetOperation(Amf, MATOP_MULT, (void(*)(void))waveHoltzMatrixVectorMultiply);
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
  
  PetscReal bNorm; // save l2-norm of RHS b

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



    if( false )
    {
      // -- zero initial guess 
      ierr = VecSet(x,0.0);CHKERRQ(ierr);

    }
    else 
    {
       // **** FIX ME ****

      // --- Set initial guess to be current v ---

      PetscScalar *xl;
      VecGetArray(x,&xl);  // get the local array from Petsc
      int iStart,iEnd;
      ierr = VecGetOwnershipRange(x,&iStart,&iEnd); CHKERRQ(ierr);
      // printf("  [iStart,iEnd]=[%i,%i]\n",iStart,iEnd);
    
      // Helmholtz solution is stored here:
      realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
      CompositeGrid & cg = *v.getCompositeGrid();

      Index I1,I2,I3;
      Real normV=0; 
      numberOfActivePoints = 0.; // count active points for scaling norm.
      int i=0;
      for( int freq=0; freq<numberOfFrequencies; freq++ )
      {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
          MappedGrid & mg = cg[grid];
          getIndex(mg.dimension(),I1,I2,I3);

          OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
          OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
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


    } // send set initial guess 

    // ---- set RHS b = A*u -----
    ierr = VecSet(u,0.0);CHKERRQ(ierr);
    // This next call will to MatMult will eventually call CgWave with initial condition u=0
    // Warning: the next call will over-write v, so we need to save v in PETSc vector x before we get to here

    iteration=computeRightHandSide;
    ierr = MatMult(Amf,u,b);CHKERRQ(ierr); 

    VecNorm(b,NORM_2,&bNorm);
    Real bNorm2h = bNorm/sqrt(numberOfActivePoints); 
    printF("solvePETSc: RHS is b: l2-norm(b)=%9.3e, L2h-norm(b)=%9.2e\n",bNorm,bNorm2h);

    iteration=-1; // Set to -1 to start the true iterations (first call to A*x does not generate a residual)


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

  if( testCase == 1)
  {
    printF("***** DONE KSP SOLVE ******\n");

    // if( false )
    // {
    //   // --- Set v = x ---
    //   PetscScalar *xl;
    //   VecGetArray(x,&xl);  // get the local array from Petsc
    //   int iStart,iEnd;
    //   ierr = VecGetOwnershipRange(x,&iStart,&iEnd); CHKERRQ(ierr);
    //   // printf("  [iStart,iEnd]=[%i,%i]\n",iStart,iEnd);
    
    //   realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");
    //   Index I1,I2,I3;
    //   int i=0;
    //   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //   {
    //     getIndex(cg[grid].dimension(),I1,I2,I3);

    //     OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
    //     FOR_3D(i1,i2,i3,I1,I2,I3)
    //     {
    //       assert( i<iEnd );
    //       vLocal(i1,i2,i3)=xl[i];
    //       i++;
    //     }
    //   }
    //   assert( i==iEnd );      

    //   v.interpolate();
      
    // }

    // IS THIS NEEDED?  Maybe if we don't monitor the residual

    printF("***** CALL AGAIN TO CHECK SOLUTION ******\n");
    int & numberOfIterations = dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken
    numberOfIterations = iteration;


    iteration=computeResidual; // This tells the matrix-vector multiply routine we are computing the residual
    ierr = MatMult(Amf,x,b);CHKERRQ(ierr);     


    Real kspResidual;
    KSPGetResidualNorm(ksp,&kspResidual);
    kspResidual /= sqrt(numberOfActivePoints);  // make an approximate L2h norm

    // add final residual
    resVector(numberOfIterations) = kspResidual;
    numberOfIterations++;    

    const real & maxResidual = dbase.get<real>("maxResidual");
    printF("\n ######## DONE KYRLOV ITERATIONS -- KSP residual=%8.2e (max-res=%8.2e) (tol=%8.2e) numIts=%i #######\n",
           kspResidual,maxResidual,tol,numberOfIterations);

  
    // if( FALSE )
    // {
    //   // *** CHECK ME SOMETHING IS WRONG HERE ****

    //   // compute the maximum residual
    //   printF(" ***solvePETSc::compute max residual...\n");

    //   Vec res;
    //   //     numberOfVects++;
    //   //     printF("VecDup res, vect object %i.\n",numberOfVects);
    //   ierr = VecDuplicate(b,&res);CHKERRQ(ierr);

    //   ierr = MatMult(Amf,x,res);CHKERRQ(ierr);   // res = A*x

    //   //     if( myid==0 ) printf("PETScSolver::res=A*x\n");
    //   //     ierr = VecView(res,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    //   //     if( myid==0 ) printf("PETScSolver::b\n");
    //   //     ierr = VecView(res,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    //   // Computes y = x + alpha y.
    //   PetscScalar alpha; alpha=-1.;
    //   // 2.2.1 VecAYPX(&alpha,b,res);  // res = b - res
    //   VecAYPX(res,alpha,b);  // res = b - res

    //   //     if( myid==0 ) printf("PETScSolver::res=b-A*x\n");
    //   //     ierr = VecView(res,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
          
    //   //     if( myid==0 ) printf("PETScSolver::b\n");
    //   //     ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      
   
    //   PetscReal maxNorm;
    //   VecNorm(res, NORM_INFINITY, &maxNorm);

    //   ierr = VecDestroy(&res);CHKERRQ(ierr);    

    //   printF(" ***solvePETSc: KSP residual=%8.2e, max residual = %8.2e \n",kspResidual,maxNorm);

    // }  

    // real maxRes = residual();
    // printF("Maximum residual = %9.3e\n",maxRes);
    
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
  
  const int & numberOfIterations = dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken  
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
