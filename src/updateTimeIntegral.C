#include "CgWave.h"
#include "ParallelUtility.h"
#include "gridFunctionNorms.h"
#include "Integrate.h"
#include "CompositeGridOperators.h"    

// lapack routines
#ifdef OV_USE_DOUBLE
  #define GETRF EXTERN_C_NAME(dgetrf)
  #define GETRI EXTERN_C_NAME(dgetri)
  #define GETRS EXTERN_C_NAME(dgetrs)
  // #define GECON EXTERN_C_NAME(dgecon)
  // #define LANGE EXTERN_C_NAME(dlange)
  // #define GEEV  EXTERN_C_NAME(dgeev)
#else
  #define GETRF EXTERN_C_NAME(sgetrf)
  #define GETRI EXTERN_C_NAME(sgetri)
  #define GETRS EXTERN_C_NAME(sgetrs)
  // #define GECON EXTERN_C_NAME(sgecon)
  // #define LANGE EXTERN_C_NAME(slange)
  // #define GEEV  EXTERN_C_NAME(sgeev)
#endif

extern "C"
{
  // PA = LU factor
  void GETRF( int & m, int & n, real & a, const int & lda, int & ipvt, int & info );
  // Solve given LU
  void GETRS( const char *trans, int & n, int & nhrs, real & a, const int & lda, int & ipvt, real & b, const int & ldb, int & info );

  // compute inverse:
  void GETRI( int & n, real & a, const int & lda, const int & ipvt, real & work, const int & iwork, int & info );

  void GECON( const char *norm, int & n, real & a, const int & lda, real & anorm, real & rcond, real & work, int & iwork, int & info );
  real LANGE( const char *norm, int & m, int & n, real & a, const int & lda, real & work );

   void GEEV( const char *jobvl, const char* jobvr, int & n, real & a, const int & lda,
              real & wr, real & wi, real &vl, int & ldvl, real & vr, int & ldvr, real & work, int & lwork, int & info );

}

#define FOR_3D(i1,i2,i3,I1,I2,I3)                                       \
int I1Base =I1.getBase(),   I2Base =I2.getBase(),  I3Base =I3.getBase(); \
int I1Bound=I1.getBound(),  I2Bound=I2.getBound(), I3Bound=I3.getBound(); \
for(int i3=I3Base; i3<=I3Bound; i3++)                                       \
  for(int i2=I2Base; i2<=I2Bound; i2++)                                     \
    for(int i1=I1Base; i1<=I1Bound; i1++)


// -------------------------------------------------------
// Evaluate sinc(x) = sin(x)/x 
// -------------------------------------------------------
Real mySinc( Real x )
{
  const Real epsilon=sqrt(REAL_EPSILON); //  what should this be ?

  Real y;
  if( fabs(x)>epsilon )
  {
    y = sin(x)/x;
  }
  else
  {
    y = 1. - x*x/6; //  first two terms in Taylor series
  }
  return y;
}


// -------------------------------------------------------
/// \brief Evaluate the discrete sinc function
/// \notes See the waveHoltz paper for a derivation 
// -------------------------------------------------------
Real sincd( Real x, Real T, Real dt )
{
  const Real epsilon=sqrt(REAL_EPSILON); //  what should this be ?

  Real z = x*T;
  Real xdtb2= x*dt*.5;

  Real y;
  if( fabs(z)>epsilon )
  {
     y=sin(z)/(T*tan(xdtb2)/(.5*dt));
  }
  else
  {
    y = 1. - z*z/6 - SQR(xdtb2)/3.;  //  first three terms in Taylor series
  }
  return y;
}

// -------------------------------------------------------------------
/// \brief Evaluate the WaveHoltz beta function (single frequency)
// -------------------------------------------------------------------
Real 
CgWave::betaWaveHoltz( Real lambda, Real omega, Real T )
{
  //  General form of 
  //        beta(lambda; omega,T) = (2/T) int0^T (cos(omega*t)-.25)*cos(lambda*t) dt 
  //  No assumptions on omega or T (i.e. we do NOT assume that T=2*pi/omega))

  Real beta = mySinc( (omega-lambda)*T )  + mySinc( (omega+lambda)*T )- .5*mySinc( lambda*T );

  return beta;
}

// -------------------------------------------------------------------------
/// \brief Evaluate the discret WaveHoltz beta function (single frequency)
/// \notes See the waveHoltz paper for a derivation 
// -------------------------------------------------------------------------
Real 
CgWave::betaDiscreteWaveHoltz( Real lambda, Real omega, Real T, Real dt )
{
  const Real alphad = tan(omega*dt*.5)/tan(omega*dt);  // adjusted alpha 
  // printF("betaDiscreteWaveHoltz : lambda=%14.8e omega=%14.8e dt=%14.8e alphad=%14.8e\n",lambda,omega,dt,alphad);
  Real beta   = sincd(omega-lambda,T,dt)  + sincd(omega+lambda,T,dt)- alphad*sincd(lambda,T,dt);  

  return beta;
}


// -----------------------------------------------------------------------------------
/// \brief Return the value of the non-zero eigenvalue of the
///     the multi-frequency WaveHoltz iteration matrix (sometimes called mu)
/// 
/// \details For a single frequency this reduces to the WaveHoltz beta function.
///  For multiple frequencies, the single non-zero eignvalue of the iteration matrix is returned.
/// See the waho.pf document for derivation of the formula used here.
/// 
/// \param lambda (input) : array of lambda values to evaluate the mu function at.
/// \param mu (output) : array of function values
/// \param useAdjusted (input) : adjust lambda for finite dt to match theory
// ---------------------------------------------------------------------------------
int CgWave::
getWaveHoltzIterationEigenvalue( RealArray & lambda, RealArray & mu,  bool useAdjusted /* = true  */ )
{
  // printF(">>Entering getWaveHoltzIterationEigenvalue.\n");

  int & numberOfFrequencies         = dbase.get<int>("numberOfFrequencies");
  const RealArray & frequencyArray  = dbase.get<RealArray>("frequencyArray");
  const RealArray & periodArray     = dbase.get<RealArray>("periodArray"); 
  const real & dt                   = dbase.get<real>("dt");  
  const int & computeEigenmodes     = dbase.get<int>("computeEigenmodes");
  const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

  RealArray A;
  getMultiFrequencyWaveHoltzMatrix( A );

  IntegerArray ipiv(numberOfFrequencies);
  int info; 
  // Factor Matrix, PA = LU 
  GETRF( numberOfFrequencies,numberOfFrequencies, A(0,0), numberOfFrequencies, ipiv(0), info );
  if( info!=0 )
  {
    printF("getWaveHoltzIterationEigenvalue:ERROR return from LU factor getrf: info=%d\n",info );
    OV_ABORT("error");
  }

  // Form A^{-1}
  int lwork=numberOfFrequencies;
  RealArray w(lwork);
  GETRI( numberOfFrequencies, A(0,0), numberOfFrequencies, ipiv(0), w(0), lwork, info );
  if( info!=0 )
  {
    printF("getWaveHoltzIterationEigenvalue:ERROR return from computing the inverse of A using GETRI: info=%d\n",info);
    OV_ABORT("error");
  }
  // -- compute w(i) = column sum of A^{-1} --
  for( int j=0; j<numberOfFrequencies; j++ )
  {
    Real temp=0.;
    for( int i=0; i<numberOfFrequencies; i++ )
      temp += A(i,j);
    w(j)=temp;
  }

  // RealArray B(numberOfFrequencies,numberOfFrequencies);
  // RealArray M(numberOfFrequencies,numberOfFrequencies);

  assert( dt>0 );

  int numLambda = lambda.getLength(0);
  for( int ilam=0; ilam<numLambda; ilam++ )
  {

    Real lam = lambda(ilam);
    if( useAdjusted ) // *wdh* added for plot filter with EigenWave
    {
      if( timeSteppingMethod==implicitTimeStepping )
      {
       // Adjust lambda for finite time-step -- important for implicit time-stepping at large CFL
       // I think we should adjust for eigenmodes too -- the beta_j will lie on the shifted curve I think
       // 
       // See EigenWave paper for derivation: 
       // D+tD-t W^n = -lam^2 (1/2) ( W^{n+1} + W^{n-1} )
       //   W^n = e^(I lamTilde dt*n ) W^0
       //  lamTilde =  acos( 1./( 1.+.5*SQR(lam*dt) ) )/dt
       lam = acos( 1./( 1.+.5*SQR(lam*dt) ) )/dt;
      }
      else
      {
        // adjust for explicit time-stepping:
        // D+tD-t W^n = - lam^2 W^n

        lam = asin(lam*dt*.5)*2./dt; 
        // printF("getWaveHoltzIterationEigenvalue: EXPLICIT time-stepping dt=%18.12e\n",dt);
      }
    }

     // --- Form mu (see formulae in waho.pdf)
     //   mu(lambda) = SUM_I beta(i)*w(i) 
     Real muTemp=0.;
     Real betac; 
     for( int i=0; i<numberOfFrequencies; i++ )
     {
       betac = betaWaveHoltz( lam, frequencyArray(i), periodArray(i) ); // use continuous beta
       Real beta = betaDiscreteWaveHoltz( lam, frequencyArray(i), periodArray(i), dt );
       muTemp += w(i)*beta;
     }
     mu(ilam)=muTemp;

     // printF(">>>> ilam=%d: lam=%12.4e mu=%12.4e\n",ilam,lam,mu(ilam));
     if( false )
       printF("getWaveHoltzIterationEigenvalue: ilam=%3d omega=%14.7e, 2*pi/T=%14.7e lambda=%14.7e lam(adjusted)=%14.7e mu=%14.7e betac=%14.7e (dt=%12.5e)\n",
            ilam,frequencyArray(0),2.*Pi/periodArray(0),lambda(ilam),lam,mu(ilam),betac,dt);

  }


  return 0;
}


// ----------------------------------------------------------------------------
/// \brief Return the multi-frequency Wave-Holtz matrix "A"
//----------------------------------------------------------------------------
int CgWave::getMultiFrequencyWaveHoltzMatrix( RealArray & A, bool useAdjusted /* = true  */  )
{
  const int & numberOfFrequencies   = dbase.get<int>("numberOfFrequencies");
  const RealArray & frequencyArray  = dbase.get<RealArray>("frequencyArray");
  const RealArray & periodArray     = dbase.get<RealArray>("periodArray");   

  const real & tFinal               = dbase.get<real>("tFinal");
  const real & dt                   = dbase.get<real>("dt");
  const int & adjustOmega           = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  RealArray & sigma                 = dbase.get<RealArray>("sigma");  

  // ---- Compute the entries using quadrature rather than the exact formulae, in some cases----
  // const bool useQuadrature = adjustOmega && numberOfFrequencies>1;
  const bool useQuadrature = numberOfFrequencies>1;

  // form matrix A
  A.redim(numberOfFrequencies,numberOfFrequencies);
  A=0.;

  // const Real alpha=.5;
  for( int j=0; j<numberOfFrequencies; j++ )
  {
    for( int i=0; i<numberOfFrequencies; i++ )
    {
      if( !useQuadrature )
      { 
        if( i==j )
          A(i,j)=1.;  
        else
          A(i,j) = betaWaveHoltz( frequencyArray(j), frequencyArray(i), periodArray(i) );
      }
      else
      {
        // --- use quadrature to evaluate the beta function ---
        const int Nt = round( tFinal/dt ); // number of time-steps
        RealArray tv(Nt+1);
        for( int i=0; i<=Nt; i++ )
          tv(i)=i*dt;

        // beta(lambda; omega,T) = (2/T) int0^T (cos(omega*t)-.25)*cos(lambda*t) dt 
        const Real lambda=frequencyArray(j), omega=frequencyArray(i), T=periodArray(i);
        const Real alphad = tan(omega*dt*.5)/tan(omega*dt);  // adjusted alpha 
        A(i,j) = sum( sigma(Range(Nt+1),i)*( (cos(omega*tv)-.5*alphad)*cos(lambda*tv) ) )* (2./T); // qudarture approximation

        Real beta = betaWaveHoltz( frequencyArray(j), frequencyArray(i), periodArray(i) );
        if( debug & 2 )
          printF("getMultiFrequencyWaveHoltzMatrix: A(%d,%d)=%16.8e (exact)   =%16.8e (quadrature)  diff=%8.2e\n",i,j,beta,A(i,j),beta-A(i,j));

      }

    }
  }
  if( debug & 2 )
  {
    printF("getMultiFrequencyWaveHoltzMatrix: A:\n");
    for( int i=0; i<numberOfFrequencies; i++ )
    {
      printF("[");
      for( int j=0; j<numberOfFrequencies; j++ )      
        printF(" (%d,%d) %12.4e, ",i,j,A(i,j));
      printF("]\n"); 
    }
    // ::display(A," getMultiFrequencyWaveHoltzMatrix ");
  }

  // OV_ABORT("stop here for now");
  // A = inv(A);
  return 0;
}



// ================================================================================================
/// \brief Initialize the time integral used by the Helmholtz solver
// ================================================================================================
int CgWave::
initializeTimeIntegral( Real dt )
{

  int & initTimeIntegral = dbase.get<int>("initTimeIntegral");
  if( !initTimeIntegral )
  {
    printF("CgWave::initializeTimeIntegral: time integral has already been initialized\n");
    return 0;
  }
  
  Real cpu0 = getCPU();
  // realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");

  const int & debug                 = dbase.get<int>("debug");
  const real & tFinal               = dbase.get<real>("tFinal");
  // const real & dt                   = dbase.get<real>("dt");
  const int & orderOfAccuracy       = dbase.get<int>("orderOfAccuracy");
  const int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");  
  const int & solveHelmholtz        = dbase.get<int>("solveHelmholtz");
  // const int & adjustOmega           = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  const int & computeEigenmodes     = dbase.get<int>("computeEigenmodes");

  const real & omega                = dbase.get<real>("omega");
  const real & Tperiod              = dbase.get<real>("Tperiod");
  const int & numPeriods            = dbase.get<int>("numPeriods");

  int & numberOfFrequencies         = dbase.get<int>("numberOfFrequencies");
  RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
  RealArray & periodArray           = dbase.get<RealArray>("periodArray");  
  // const int & deflateWaveHoltz      = dbase.get<int>("deflateWaveHoltz");
  // const int & useAugmentedGmres     = dbase.get<int>("useAugmentedGmres"); 

  const int & filterTimeDerivative  = dbase.get<int>("filterTimeDerivative");
  const int & useFilterWeights      = dbase.get<int>("useFilterWeights"); 
  const int & filterD0t             = dbase.get<int>("filterD0t"); 


  // const int includeParallelGhost=1;

  Index I1,I2,I3;

  const int & numCompWaveHoltz = dbase.get<int>("numCompWaveHoltz");
  if( filterTimeDerivative )
    assert( numCompWaveHoltz==2 );
  else
    assert( numCompWaveHoltz==numberOfFrequencies );

  // numCompWaveHoltz = filterTimeDerivative ? 2 : numberOfFrequencies;   // ******* MOVE THIS *****

  if( !solveHelmholtz )
  {
    // -- we must be just testing the integral ---
    numberOfFrequencies=1;
    frequencyArray.redim(numberOfFrequencies);
    frequencyArray(0)=omega;
    periodArray.redim(numberOfFrequencies);
    periodArray(0)=twoPi/frequencyArray(0);
  }

  if( fabs(tFinal -periodArray(0)) > 100.*REAL_EPSILON*tFinal )
  {
    printF("CgWave::initializeTimeIntegral: tFinal != periodArray(0) : tFinal=%12.4e, periodArray(0)=%12.4e, diff=%12.4e\n",
         tFinal,periodArray(0),tFinal-periodArray(0));

    OV_ABORT("error");
  }

  // if( false )
  // {
  //   // For testing over-write the computed solution with the exact solution
  //   const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");
  //   if( knownSolutionOption == "userDefinedKnownSolution" )
  //   {
  //     for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
  //     {
  //       MappedGrid & mg = cg[grid];
  //       getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points.              
  //       // -- User defined known solution ---
  //       getUserDefinedKnownSolution( t, grid, u[grid], I1,I2,I3 );
  //     }
  //   }
  // }

  // const bool firstStep = stepOption==0;
  // const bool lastStep  = stepOption==2;

  const int Nt = round( tFinal/dt ); // number of time-steps

  // sigma(0:Nt,0:numFreq-1) : integration weights
  if( !dbase.has_key("sigma") )
  {
    dbase.put<RealArray>("sigma");         // holds integration weights
    dbase.put<RealArray>("filterWeights"); // holds *new* integration weights
  }
  RealArray & sigma = dbase.get<RealArray>("sigma");  
  RealArray & filterWeights = dbase.get<RealArray>("filterWeights");  

  if( true )
    printF("\n #### initializeTimeIntegral: omega=%16.12e, dt=%16.12e, useFilterWeights=%d filterD0t=%d ####\n\n",omega,dt,useFilterWeights,filterD0t);

  int orderOfQuadrature = orderOfAccuracy; 
  if( true && ( computeEigenmodes || numberOfFrequencies==1 ) )
  {
     orderOfQuadrature = orderOfAccuracyInTime; // *wdh* Sept 24, 2024 : Trapezoidal rule should be fine if orderInTime=2, also matches theory 
  }

  // --- Evaluate the integration weights for each of the frequencies ---
  if( useFilterWeights )
    getFilterWeights( Nt, dt, numberOfFrequencies, periodArray, orderOfQuadrature, sigma, filterWeights );
  else 
    getIntegrationWeights( Nt, numberOfFrequencies, periodArray, orderOfQuadrature, sigma );


  // // WHY IS THIS DONE HERE ????
  // if( deflateWaveHoltz && useAugmentedGmres==1 )
  // { // initialize deflation for Augmented Gmres
  //   initializeDeflation();
  // }

  initTimeIntegral=false;

  timing(timeForTimeIntegral) += getCPU()- cpu0;


  return 0;
}


// ================================================================================================
/// \brief Update the time integral used by the Helmholtz solver
///
/// The WaveHoltz algorithm computes the time integral of the solution. 
/// This routine increments the time integral using the current solution and stores 
/// the result in "v" (from the dbase).
///
/// \param step (input) : current step 
/// \param stepOption : firstStep, middleStep, lastStep. When stepOption==firstStep, the integration weights will be computed and the
///     time integral will be initialized.
/// \param t (input) : current time
/// \param u (input) : current solution
// ================================================================================================
int CgWave::
updateTimeIntegral( int step, StepOptionEnum stepOption, Real t, Real dt, realCompositeGridFunction& u )
{
  Real cpu0 = getCPU();
  
  realCompositeGridFunction & v     = dbase.get<realCompositeGridFunction>("v");

  const int & debug                 = dbase.get<int>("debug");
  // const real & tFinal               = dbase.get<real>("tFinal");
  // const real & dt                   = dbase.get<real>("dt");
  const int & orderOfAccuracy       = dbase.get<int>("orderOfAccuracy");
  const int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");  
  const int & solveHelmholtz        = dbase.get<int>("solveHelmholtz");
  const int & adjustOmega           = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  const int & computeEigenmodes     = dbase.get<int>("computeEigenmodes");

  const real & omega                = dbase.get<real>("omega");
  const real & Tperiod              = dbase.get<real>("Tperiod");
  const int & numPeriods            = dbase.get<int>("numPeriods");

  int & numberOfFrequencies         = dbase.get<int>("numberOfFrequencies");
  RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
  RealArray & periodArray           = dbase.get<RealArray>("periodArray");  
  const int & deflateWaveHoltz      = dbase.get<int>("deflateWaveHoltz");
  const int & useAugmentedGmres     = dbase.get<int>("useAugmentedGmres"); 

  const int & filterTimeDerivative = dbase.get<int>("filterTimeDerivative");
  const int & useFilterWeights      = dbase.get<int>("useFilterWeights"); 
  const int & filterD0t             = dbase.get<int>("filterD0t"); 

  const int includeParallelGhost=1;

  Index I1,I2,I3;
  Index J1,J2,J3;

  const int & numCompWaveHoltz = dbase.get<int>("numCompWaveHoltz");
  if( filterTimeDerivative )
    assert( numCompWaveHoltz==2 );
  else
    assert( numCompWaveHoltz==numberOfFrequencies );

  // numCompWaveHoltz = filterTimeDerivative ? 2 : numberOfFrequencies;   // ******* MOVE THIS *****

  if( !solveHelmholtz )
  {
    // -- we must be just testing the integral ---
    numberOfFrequencies=1;
    frequencyArray.redim(numberOfFrequencies);
    frequencyArray(0)=omega;
    periodArray.redim(numberOfFrequencies);
    periodArray(0)=twoPi/frequencyArray(0);
  }

  // if( stepOption==firstStep && adjustOmega && numberOfFrequencies==1 )
  // {
  //   // ***  DO THIS FOR NOW *****
  //   printF("\n updateTimeIntegral: WARNING: changing frequencyArray and periodArray for adjustOmega\n\n");
  //   frequencyArray(0) = omega;
  //   periodArray(0)    = twoPi/omega; 
  // }

  if( false)
  {
   printF("updateTimeIntegral: stepOption=%d omega=%16.12e, t=%12.4e, dt=%16.12e, t/dt=%.3g\n",stepOption,omega,t,dt,t/dt);
   printF("updateTimeIntegral: numberOfFrequencies=%d\n",numberOfFrequencies);
   // u.display(sPrintF("u for updateTimeIntegral, t=%9.3e",t),"%6.2f ");
  }

  // if( fabs(tFinal -periodArray(0)) > 100.*REAL_EPSILON*tFinal )
  // {
  //   printF("CgWave::updateTimeIntegral: tFinal != periodArray(0) : tFinal=%12.4e, periodArray(0)=%12.4e, diff=%12.4e\n",
  //        tFinal,periodArray(0),tFinal-periodArray(0));

  //   OV_ABORT("error");
  // }

  if( false )
  {
    // For testing over-write the computed solution with the exact solution
    const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");
    if( knownSolutionOption == "userDefinedKnownSolution" )
    {
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        MappedGrid & mg = cg[grid];
        getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points.              
        // -- User defined known solution ---
        getUserDefinedKnownSolution( t, grid, u[grid], I1,I2,I3 );
      }
    }
  }

  // const bool firstStep = stepOption==0;
  // const bool lastStep  = stepOption==2;

  // const int Nt = round( tFinal/dt ); // number of time-steps
  const int Nt = round( periodArray(0)/dt ); // number of time-steps

  // sigma(0:Nt,0:numFreq-1) : integration weights
  if( !dbase.has_key("sigma") )
  {
    dbase.put<RealArray>("sigma");         // holds integration weights
    dbase.put<RealArray>("filterWeights"); // holds *new* integration weights
  }
  RealArray & sigma = dbase.get<RealArray>("sigma");  
  RealArray & filterWeights = dbase.get<RealArray>("filterWeights");  

  if( sigma.getLength(0)!= Nt+1 )
  {
    printF("updateTimeIntegral:ERROR: sigma.getLength(0)!= Nt+1, sigma.getLength(0)=%d Nt=%d\n");
    OV_ABORT("error");
  }

  const bool useOpt=true;

  // if( FALSE && stepOption==firstStep )
  // {
  //   assert( t==0. && step==0 );
  //   if( false )
  //     printF("updateTimeIntegral: stepOption=%d omega=%16.12e, t=%12.4e, dt=%16.12e, t/dt=%.3g useFilterWeights=%d\n",stepOption,omega,t,dt,t/dt,useFilterWeights);

  //   int orderOfQuadrature = orderOfAccuracy; 
  //   if( true && ( computeEigenmodes || numberOfFrequencies==1 ) )
  //   {
  //      orderOfQuadrature = orderOfAccuracyInTime; // *wdh* Sept 24, 2024 : Trapezoidal rule should be fine if orderInTime=2, also matches theory 
  //   }

  //   // --- Evaluate the integration weights for each of the frequencies ---
  //   if( useFilterWeights )
  //     getFilterWeights( Nt, dt, numberOfFrequencies, periodArray, orderOfQuadrature, sigma, filterWeights );
  //   else 
  //     getIntegrationWeights( Nt, numberOfFrequencies, periodArray, orderOfQuadrature, sigma );


  //   // WHY IS THIS DONE HERE ????
  //   if( deflateWaveHoltz && useAugmentedGmres==1 )
  //   { // initialize deflation for Augmented Gmres
  //     initializeDeflation();
  //   }

  // }


  if( stepOption==lastStep )
  {
    int numTimeSteps=Nt;
    if( useFilterWeights && filterTimeDerivative && filterD0t )
      numTimeSteps++; //  we take one extra step

    if( step!=numTimeSteps )
    {
      const real & tFinal = dbase.get<real>("tFinal");
      printF("updateTimeIntegral:WARNING: lastStep but step=%d != numTimeSteps=%d\n",step,numTimeSteps);
      printF("tFinal=%14.6e, dt=%12.4e (Nt = round(tFinal/dt)\n",tFinal,dt);
      OV_ABORT("error");
    }

    // const Real tol = REAL_EPSILON*10000.*tFinal; 
    // if( !( fabs(t-tFinal)< tol ) )
    // {
    //   printF("CgWave:ERROR: t != tFinal, t=%14.6e tFinal=%14.6e diff=%9.3e, tol=%9.3e\n",t,tFinal,fabs(t-tFinal),tol);
      
    //   OV_ABORT("ERROR");
    // }
  }

  
  if( stepOption==firstStep )
  {
    // When solving the Helmholtz problem with CgWaveHoltz we need to evaluate an integral 
    //      v  = (1/(2*T)* Int_0^T [  ( cos(omega*t)-.25)*u(x,t) dt ] 


    realCompositeGridFunction & vOld = dbase.get<realCompositeGridFunction>("vOld");
    if( false )
      vOld.display("First Step: vOld","%5.2f ");


    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
      OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
      OV_GET_SERIAL_ARRAY(real,vOld[grid],vOldLocal);
      getIndex(cg[grid].dimension(),I1,I2,I3);

      getIndex(cg[grid].gridIndexRange(),J1,J2,J3);

      bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeParallelGhost);
      if( ok )
      {
        if( useFilterWeights && filterTimeDerivative && filterD0t )
        {
          const Real & viFactor = dbase.get<Real>("viFactor");
          // printF("updateTimeIntegral: firstStep: step=%d : filterD0t omega=%9.2e viFactor=%9.2e fw(0)=%9.2e fw(1)=%9.2ef\n",
          //     step,omega,viFactor,filterWeights(0,1),filterWeights(1,1));

          // printF("updateTimeIntegral: max(fabs(v(:,:,:,0))=%9.2e max(fabs(v(:,:,:,1))=%9.2e max(fabs(vOld(:,:,:,1))=%9.2e\n",
          //     max(fabs(vLocal(J1,J2,J3,0))),max(fabs(vLocal(J1,J2,J3,1))),max(fabs(vOldLocal(J1,J2,J3,1))));

          vLocal(I1,I2,I3,0) = filterWeights(step,0)*uLocal(I1,I2,I3);

          // CHECK ME: At the start vLocal(I1,I2,I3,1) should hold D0t W^0 
          const int nct=1;
           
          // // NOTE: currently viFactor=1 which means vLocal(I1,I2,I3,1) holds ui and NOT the time derivative
          // // THUS we need to scale by the symbol of D0t 
          // const Real scale = (sin(omega*dt)/dt)/viFactor; // 

          vLocal(I1,I2,I3,1) = filterWeights(0,nct)*vLocal(I1,I2,I3,1)      // D0t W^0 = vDot
                              +filterWeights(1,nct)*uLocal(I1,I2,I3);       // W^0

          // printF("  after first step:  |v(:,:,:,1)|=%9.2e uLocal=[%9.2e,%9.2e]\n",max(fabs(vLocal(J1,J2,J3,1))),
          //     min(uLocal(J1,J2,J3)), max(uLocal(J1,J2,J3)));                    

        }
        else
        {
          for( int freq=0; freq<numCompWaveHoltz; freq++ )
          {
            // vLocal = ( .5*( cos(omega*(t))-.25 ) )*uLocal;  // Trapezoidal first term (.5)
            if( useFilterWeights )
            {
              vLocal(I1,I2,I3,freq) = filterWeights(step,freq)*uLocal(I1,I2,I3);
            }
            else
            {
              const Real omegaFreq = frequencyArray(freq); 
              const Real alphad = tan(omegaFreq*dt*.5)/tan(omegaFreq*dt);  // adjusted alpha 
              vLocal(I1,I2,I3,freq) = ( sigma(step,freq)*( cos(omegaFreq*(t))-.5*alphad ) )*uLocal(I1,I2,I3);  // Trapezoidal first term (.5)
              // vLocal(I1,I2,I3,freq) = ( sigma(step,freq)*( cos(omegaFreq*(t))-.25 ) )*uLocal(I1,I2,I3);  // Trapezoidal first term (.5)
            }
          }
        }
      }
  
    }
  }
  else
  {
    // Integral in time -- trapezoidal rule
    // assert( i<=(numberOfTimeSteps-1) );
      
    // On the last step we scale by .5 for trap, and normalize by (2./Tperiod)*dt
    // bool lastStep = i== (numberOfTimeSteps-1);

    const real trapFactor= stepOption==lastStep ? .5 : 1.;
      
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);   // solution at new time 
      OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
      getIndex(cg[grid].dimension(),I1,I2,I3);
      getIndex(cg[grid].gridIndexRange(),J1,J2,J3);
      bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeParallelGhost);
      if( ok )
      {

        const Real *pu  = uLocal.getDataPointer();
        Real *pv        = vLocal.getDataPointer();
        const int i1Base=uLocal.getBase(0), i2Base=uLocal.getBase(1), i3Base=uLocal.getBase(2);
        const int nd1= uLocal.getLength(0), nd2= uLocal.getLength(1), nd3= uLocal.getLength(2);
        #define ua(i1,i2,i3)       pu[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base          ))]
        #define va(i1,i2,i3,m)     pv[((i1)-i1Base)+nd1*((i2)-i2Base + nd2*((i3)-i3Base +nd3*(m) ))]

        if( useFilterWeights && filterTimeDerivative && filterD0t )
        {

          vLocal(I1,I2,I3,0) += filterWeights(step  ,0)*uLocal(I1,I2,I3);  // filter W 
          vLocal(I1,I2,I3,1) += filterWeights(step+1,1)*uLocal(I1,I2,I3);  // add contribution for D0t W

          // printF("updateTimeIntegral: step=%d : filterD0t t=%9.3e fw(%d,1)=%10.3e |v(:,:,:,1)|=%9.2e\n",step,t,step+1,filterWeights(step+1,1), max(fabs(vLocal(J1,J2,J3,1))));
        }
        else
        {
          for( int freq=0; freq<numCompWaveHoltz; freq++ )
          {
            if( useFilterWeights )
            {
              if( useOpt )
              {
                const Real fw=filterWeights(step,freq);
                FOR_3D(i1,i2,i3,I1,I2,I3)
                  va(i1,i2,i3,freq) += fw*ua(i1,i2,i3);
              }
              else
              {
                vLocal(I1,I2,I3,freq) += filterWeights(step,freq)*uLocal(I1,I2,I3);
              }
            }
            else
            {
              const Real omegaFreq = frequencyArray(freq);
              const Real alphad = tan(omegaFreq*dt*.5)/tan(omegaFreq*dt);  // adjusted alpha  
              if( useOpt )
              {
                const Real factor=( sigma(step,freq)*( cos(omegaFreq*(t))-.5*alphad ) );
                FOR_3D(i1,i2,i3,I1,I2,I3)
                  va(i1,i2,i3,freq) += factor*ua(i1,i2,i3);
              }
              else
              {
                vLocal(I1,I2,I3,freq) += ( sigma(step,freq)*( cos(omegaFreq*(t))-.5*alphad ) )*uLocal(I1,I2,I3);
              }
              // vLocal(I1,I2,I3,freq) += ( sigma(step,freq)*( cos(omegaFreq*(t))-.25 ) )*uLocal(I1,I2,I3);
            }
          }
        }

        if( stepOption==lastStep && useFilterWeights==0 )
        {
          for( int freq=0; freq<numberOfFrequencies; freq++ )
          {
            vLocal(I1,I2,I3,freq) *= (2./periodArray(freq));  //  normalize the time integral  
          }        
        }
  
      }
       
    } // end for grid 
    if( stepOption==lastStep )
    {
      
      // if( false || (false && numberOfFrequencies==1) )  // ***************** Oct 31, 2021 NOT valid for multi-freq
      // {
      //   bool applyExplicitBoundaryConditions=true;
      //   applyBoundaryConditions( v,v, t, applyExplicitBoundaryConditions );
      // }

      if( false && filterTimeDerivative )
      {
        // ** testing ***
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
          getIndex(cg[grid].dimension(),I1,I2,I3);

          OV_GET_SERIAL_ARRAY(real,v[grid],vLocal); 

          Real vMax    = max(abs(vLocal(I1,I2,I3,0)));
          Real vDotMax = max(abs(vLocal(I1,I2,I3,1)));
          printF(">>>>>> INFO: max(v)=%9.2e, max(vDot) = %9.2e (useFilterWeights=%d)\n",vMax,vDotMax,useFilterWeights);

        }       
      }

      // ---- MULTI-FREQUENCY PROJECTION -----
      if( numberOfFrequencies > 1 )
      {
        RealArray A;
        getMultiFrequencyWaveHoltzMatrix( A );

        IntegerArray ipiv(numberOfFrequencies);
        int info; 
        // Factor Matrix, PA = LU 
        GETRF( numberOfFrequencies,numberOfFrequencies, A(0,0), numberOfFrequencies, ipiv(0), info);
        if( info!=0 )
        {
          printF("updateTimeIntegral:ERROR return from LU factor getrf: info=%d\n",info);
          OV_ABORT("error");
        }

        int ldb = numberOfFrequencies;
        int nrhs=1; // number of right-hand sides
        RealArray b(numberOfFrequencies);
        // b=0.;

        // --- loop over grid points and project each point ---
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
          getIndex(cg[grid].dimension(),I1,I2,I3);

          OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
          bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeParallelGhost);
          if( ok )
          {
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
              // --- fill in the right-hand-side --
              for( int freq=0; freq<numberOfFrequencies; freq++ )
              {
                b(freq) = vLocal(i1,i2,i3,freq);
              }
              // Solve A x = b 
              // TRANS : "N", "T" "C"
              // Could solve many at once
              // Could form inverse and do multiplication inline (?)
              GETRS( "N", numberOfFrequencies, nrhs, A(0,0), numberOfFrequencies, ipiv(0), b(0), ldb, info );
              if( info!=0 )
              {
                printF("updateTimeIntegral:ERROR return from Triangular solve: getrs: info=%d\n",info);
                OV_ABORT("error");
              }
              // --- over-write v with projected value ---
              for( int freq=0; freq<numberOfFrequencies; freq++ )
              {
                vLocal(i1,i2,i3,freq) = b(freq);
              }

              // if( filterTimeDerivative )
              // {
              //   const int nwt=1;   // component number for wt 
              //   b(0) = vLocal(i1,i2,i3,nwt);
              //   GETRS( "N", numberOfFrequencies, nrhs, A(0,0), numberOfFrequencies, ipiv(0), b(0), ldb, info );
              //   assert( info==0 );

              //   vLocal(i1,i2,i3,nwt) = b(0);

              //   vLocal(i1,i2,i3,nwt) = 0.; // test *********************************

              // }


            }
          }
        }
      }

      const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");

      // u.display(sPrintF("u after integral, t=%9.3e",t),"%6.2f ");
      if( false )
        v.display(sPrintF("v after WaveHoltz projection integral, t=%9.3e",t),"%6.2f ");

      if( ( debug & 2 ) && knownSolutionOption=="userDefinedKnownSolution" ) 
      {
        real maxErr = getErrors( v, t );
        printF("\n ** cgWave:updateTimeIntegral: t=%9.3e, error in WaveHoltz v:  maxErr=%9.2e ** \n\n",t,maxErr);
      }
        
      
      const int & computeEigenmodesWithSLEPc = dbase.get<int>("computeEigenmodesWithSLEPc");
      if( computeEigenmodes && !computeEigenmodesWithSLEPc ) 
      {
        updateEigenmodes();

        Real relErrEigenvalue, relErrEigenvector;
        getErrorsInEigenmodes( relErrEigenvalue, relErrEigenvector );

      }

      if( deflateWaveHoltz && useAugmentedGmres==0 )
      { // Deflate the solution by projecting out selected eigenvectors
        deflateSolution();
      }
    }
      
  }

  timing(timeForTimeIntegral) += getCPU()- cpu0;

  return 0;
}

// =================================================================================
/// \brief Compute the filter weights in time for the WaveHoltz integrals
///
/// \param Nt (input) : number of time steps 
/// \param numFreq (input) : number of frequencies
/// \param Tv(0:numFreq-1) (input) : periods
/// \param orderOfAccuracy (input) : order of accuracy of the quadrature
///
/// \param sigma(0:Nt,0:numFreq-1) (output) : integration weights
/// \param filterWeights(0:Nt,0:numFreq-1) (output) : filter weights
///
// =================================================================================
int CgWave::getFilterWeights( int Nt, Real dt, int numFreq, const RealArray & Tv, int orderOfAccuracy, RealArray & sigma, RealArray & filterWeights  )
{

  const int & useFilterWeights = dbase.get<int>("useFilterWeights"); 
  const int & filterD0t        = dbase.get<int>("filterD0t"); 

  getIntegrationWeights( Nt, numFreq, Tv, orderOfAccuracy,sigma );

  if( useFilterWeights )
  {
    // alpha=.5;
    // for( freq=1:numFreq )
    //   omega = par.frequencyArray(freq);
    //   T     = par.periodArray(freq);     
    //   for n=1:Nt+1
    //     t = (n-1)*dt; 
    //     par.fw(n,freq) = par.sigma(n,freq)*( (cos(omega*t) - .5*alpha ) )*(2/T);
    //   end
    // end

    const int & filterTimeDerivative = dbase.get<int>("filterTimeDerivative"); 
    const int & filterD0t            = dbase.get<int>("filterD0t");
    RealArray & frequencyArray       = dbase.get<RealArray>("frequencyArray");
    RealArray & frequencyArraySave   = dbase.get<RealArray>("frequencyArraySave");

    if( true ) printF("\n ##### getFilterWeights filterTimeDerivative=%d filterD0t=%d #####\n\n",filterTimeDerivative,filterD0t);


    const Real dt = Tv(0)/Nt; 

    int & numCompWaveHoltz = dbase.get<int>("numCompWaveHoltz");
    // numCompWaveHoltz = filterTimeDerivative ? 2 : numFreq;   // ******* MOVE THIS *****

    int Ntfw =Nt+1; // number of filter weights
    if( filterTimeDerivative && filterD0t ) 
    {
      //if we filter D0t W we need to take and extra time-step and save the weights for D0t W^0
      Ntfw = Nt+3;  //  we need weights for Vdot, W^n n=0,1,2,...,Nt+1
    }

    filterWeights.redim(Ntfw,numCompWaveHoltz); 
    filterWeights=0.;

    // Real alpha = 0.5;

    // if( numFreq==1 )
    // {
    //   // Adjust alpha to make discrete beta function have a max value of 1: (see overHoltz paper)
    //   // *wdh* JUly 20, 2024
    //   alpha = tan(frequencyArray(0)*dt*.5)/tan(frequencyArray(0)*dt);  
    //   // printF("\n *** Adjust filter constant for dt: alpha=%6.3f (usual=.5) ***\n\n",alpha);
    //   // OV_ABORT("STOP here for now");
    // }
    for( int freq=0; freq<numFreq; freq++ )
    {
      const Real omega = frequencyArray(freq);
      const Real alphad = tan(omega*dt*.5)/tan(omega*dt);  // adjusted alpha *wdh* Sept 10, 2023
      if( false )
        printF("getFilterWeights: omega=%14.8e, dt=%14.8e, alphad=%14.8e\n",omega,dt,alphad);

      for( int n=0; n<=Nt; n++ )
      {
        Real t = n*dt;
        filterWeights(n,freq) = sigma(n,freq)*( cos(omega*t) - .5*alphad )*(2./Tv(freq));
      }
    }

    if( filterTimeDerivative )
    {
      // -- integration weights for the time derivative ---
      if( filterD0t )
      {
        // --------------------------------------------------------------------------------------
        //  vDot^{k+1} = (2/T) SUM_{n=0}^Nt sigma_n ( cos(omega*t^n) - alpha/2) D0t W^n   (base 0 sigma_n)
        //             = SUM_{n=1}^Nt S_n D0t W^n *(2*dt) 
        //      = S_0*2*dt*vDot + S_1*( W^2-W^0) + S_2*( W^3-W^1) + S_3*(W^4-W^2) + ...
        // 
        //      = S_0*2*dt*vDot - S_1*W^0 -S_2*W^1 + (S_1-S_3)*W^2 + (S_2-S_4)*W^3 + ...
        //      = S_0*2*dt*vDot - S_1*W^0 -S_2*W^1 
        //        +  SUM_{n=2}^{Nt-1} (S_{n-1} - S_{n+1} ) W^n  
        //        + S_{Nt-1}*W^{Nt} + S_{Nt}*W^{Nt+1}
        // 
        // where 
        //     S_n = sigma_n  * (2/Tbar) * ( cos(omega*t^n) - alpha/2)  * (1/(2*dt))
        // with
        //     Dot W^0 = vDot 
        // 
        // -------------------------------------------------------------------------------------- 

        // printF("**** Form filter weights for D0t W...\n");
        int freq=0;
        const Real omega  = frequencyArray(freq);
        const Real alphad = tan(omega*dt*.5)/tan(omega*dt);  // adjusted alpha *wdh* Sept 10, 2023        
        const Real T      = Tv(freq);

        // NOTE: currently viFactor=1 which means vLocal(I1,I2,I3,1) holds ui and NOT the time derivative
        // THUS we need to scale by the symbol of D0t 
        const Real & viFactor = dbase.get<Real>("viFactor");
        const Real scaleD0t = viFactor/(sin(omega*dt)/dt);
        const Real scale = scaleD0t;

        const int nct = 1; // component number for time derivative in filterWeights


        filterWeights(0,nct) =  sigma(0,freq)*( cos(omega*(0.*dt))-.5*alphad )*(2./T);               // coeff of vDot NOTE: no scale here
        filterWeights(1,nct) = -sigma(1,freq)*( cos(omega*(1.*dt))-.5*alphad )*scale*(2./T)/(2.*dt); // coeff of W^n = W^0
        filterWeights(2,nct) = -sigma(2,freq)*( cos(omega*(2.*dt))-.5*alphad )*scale*(2./T)/(2.*dt); // coeff of W^n = W^1
        for( int n=2; n<=Nt-1; n++ )
        {
          Real t = n*dt;
          filterWeights(n+1,nct) = (  sigma(n-1,freq)*( cos(omega*(t-dt)) - .5*alphad )*(2./T)
                                     -sigma(n+1,freq)*( cos(omega*(t+dt)) - .5*alphad )*(2./T)
                                   )*scale/(2.*dt);
        }
        filterWeights(Nt+1,nct) = sigma(Nt-1,freq)*( cos(omega*((Nt-1)*dt))-.5*alphad )*scale*(2./T)/(2.*dt); // coeff of W^{n} = W^{Nt}
        filterWeights(Nt+2,nct) = sigma(Nt  ,freq)*( cos(omega*((Nt  )*dt))-.5*alphad )*scale*(2./T)/(2.*dt); // coeff of W^{n} = W^{Nt+1}

        // n=0; t=n*dt; par.fw(n+nBase,ncD0t) =  par.sigma(n   +nBase,1)*( cos(omega*(t   ))-.5*alpha )*(2/T)/(2*dt); % coeff of vDot

        // n=0; t=n*dt; par.fw(n+nBase1,ncD0t) = -par.sigma(n+1+nBase,1)*( cos(omega*(t+dt))-.5*alpha )*(2/T)/(2*dt); % coeff of W^n = W^0
        // n=1; t=n*dt; par.fw(n+nBase1,ncD0t) = -par.sigma(n+1+nBase,1)*( cos(omega*(t+dt))-.5*alpha )*(2/T)/(2*dt); % coeff of W^n = W^1
        // for n=2:Nt-1
        //   % coeff of W^{n} 
        //   t = n*dt;  
        //   coeff = (  par.sigma(n-1+nBase,1)*( cos(omega*(t-dt))-.5*alpha )*(2/T) ...
        //            - par.sigma(n+1+nBase,1)*( cos(omega*(t+dt))-.5*alpha )*(2/T) )/(2*dt); 
        //   par.fw(n+nBase1,ncD0t) = coeff;
        // end 
        // n=Nt;   t=n*dt; par.fw(n+nBase1,ncD0t) = par.sigma(n-1+nBase,1)*( cos(omega*(t-dt))-.5*alpha )*(2/T)/(2*dt); % coeff of W^{n} = W^{Nt}
        // n=Nt+1; t=n*dt; par.fw(n+nBase1,ncD0t) = par.sigma(n-1+nBase,1)*( cos(omega*(t-dt))-.5*alpha )*(2/T)/(2*dt); % coeff of W^{n} = W^{Nt+1}




      }
      else
      {
        //  ---- Use integration by parts ----
        //  
        //   INT (cos(omega*t)-1/4) w_t dt = [(cos(omega*t)-1/4) w] + omega INT sin(omega*t) w dt
        //   
        const Real omega = frequencyArray(0);
        const Real T = Tv(0);

        const Real & damp     = dbase.get<real>("damp"); // coeff of linear damping 
        const Real & viFactor = dbase.get<Real>("viFactor"); 

        // **NOTE**
        // viFactor must match in
        //     getFilterWeights
        //     takeFirstStep
        //     solveHelmholtz : definition of Im(u)
        // omegaTilde must match 1/viFactor in the definition of ui from vDot,  in u = ur*cos(omega*t) + ui*sin(omega*t)
        Real omegaTilde = omega; 
        if( true )
        {
          omegaTilde = viFactor; // *new way* June 19, 2023
        }
        else if( true )
        {
          omegaTilde = omega; // adjusted 
        }
        else if( true || damp != 0. ) 
        {
          //  --- adjust omega to match D0t used in damping: ----
          const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

          if( true || timeSteppingMethod==CgWave::implicitTimeStepping )
            omegaTilde = sin(omega*dt)/dt; 
          else
            omegaTilde = frequencyArraySave(0);  // unadjusted
        }

        const int nct = 1; // component number for time derivative in filterWeights

        Real t;
        for( int n=0; n<=Nt; n++ )
        {
          t = n*dt;
          filterWeights(n,nct) = sigma(n,0)*( omegaTilde* sin(omega*t) )*(2./T);
          // filterWeights(n,nct) = sigma(n,0)*( omegaTilde*omega sin(omega*t) )*(2./T);
        }
        // add boundary terms from IBP's
        const Real alphad = tan(omega*dt*.5)/tan(omega*dt);  // adjusted alpha *wdh* Sept 10, 2023 **CHECK ME***
        t=0; 
        filterWeights( 0,nct) += -(omegaTilde/omega)* (cos(omega*t) - .5*alphad )*(2./T);
        t = Nt*dt; 
        filterWeights(Nt,nct) +=  (omegaTilde/omega)* (cos(omega*t) - .5*alphad )*(2./T);   

      }
    }


    if( 1==1 && filterD0t )
    {
      // CHECK FILTER WEIGHTS 
      Real omega = frequencyArray(0);
      const Real & viFactor = dbase.get<Real>("viFactor");
      const Real scaleD0t = viFactor/(sin(omega*dt)/dt);
      const Real scale = scaleD0t;      

      Real tf = Tv(0); 
      RealArray tv(Nt+2), f(Nt+2), ft(Nt+2);
     
      const int nlam=6; 

      for( int ilam=0; ilam<nlam; ilam++ )
      {
        Real lambda;
        lambda = omega*(1. + .123*ilam);
        for( int i=0; i<=Nt+1; i++ )
        {
          tv(i) = i*dt; // tf/Real(Nt); 

          if( ilam<=2 )
          {
            f(i)  = cos(lambda*tv(i));
            ft(i) = ( cos(lambda*(tv(i)+dt)) - cos(lambda*(tv(i)-dt)) )/(2.*dt);  // D0t f 
          }
          else
          {
            f(i)  = sin(lambda*tv(i));
            ft(i) = ( sin(lambda*(tv(i)+dt)) - sin(lambda*(tv(i)-dt)) )/(2.*dt);  // D0t f             
          }
        }

        // ---- scale the time derivative ---
        ft *= scale; 

        Real sum1=0., sum2=0., sum1e=0., sum2e=0.;

        sum2 = filterWeights(0,1)*ft(0);   
        for( int i=0; i<=Nt+1; i++ )
        {
          sum1  += filterWeights(i  ,0)*f(i);  // filter f

          sum2  += filterWeights(i+1,1)*f(i);  // filter D0t(f) 

          sum1e += filterWeights(i  ,0)*f(i);        // do this for now
          sum2e += filterWeights(i  ,0)*ft(i); // filter D0t f for checking -- this should be the exact answer
        }

        if( ilam<=2 )
          sum1e = betaDiscreteWaveHoltz( lambda, frequencyArray(0), Tv(0), dt );

        Real err1 = sum1 -sum1e;
        Real err2 = sum2-sum2e;
        printF("Check filter weights: omega=%10.3e omegaTilde=%10.3e, lambda=%10.3e sum1=%10.3e (err=%9.2e) sum2=%10.3e (err=%9.2e)\n",
            frequencyArraySave(0),frequencyArray(0), lambda,sum1,err1,sum2,err2);
      }
      // OV_ABORT("stop here for now");

    }

  } // end useFilterWeights
  return 0;

}


// =================================================================================
/// \brief Compute the integration weights in time for the WaveHoltz integrals
///
/// \param Nt (input) : number of time steps 
/// \param numFreq (input) : number of frequencies
/// \param Tv(0:numFreq-1) (input) : periods
/// \param orderOfAccuracy (input) : order of accuracy of the quadrature
///
/// \param sigma(0:Nt,0:numFreq-1) (output) : integration weights
///
/// \authors WDH and Cassandra Carrick -- Nov, 2021
// =================================================================================
int CgWave::getIntegrationWeights( int Nt, int numFreq, const RealArray & Tv, int orderOfAccuracy, RealArray & sigma )
{
  if( false )
     printF("\n ##### getIntegrationWeights orderOfAccuracy=%d #####\n\n",orderOfAccuracy);

  if( orderOfAccuracy!=2 && orderOfAccuracy!=4 )
  {
    printF("getIntegrationWeights:ERROR: orderOfAccuracy=%d is not supported\n",orderOfAccuracy);
    OV_ABORT("error");
  }

  if( Tv(0) < max(abs(Tv)) )
  {
    printF("getIntegrationWeights:ERROR: Tv(0)=%8.2e is not the largest!\n",Tv(0));
    ::display(Tv,"Tv","%9.2e ");
    OV_ABORT("error");
  }

  const int & debug = dbase.get<int>("debug");

  Real dt = Tv(0)/Nt; 

  sigma.redim(Nt+1,numFreq); // note: base 0 

  sigma = dt;  // set all weights to dt by default, this is not really needed

  // trapezoidal rule for Tv(0) -- this is the largest period 
  sigma( 0,0)=.5*dt; 
  sigma(Nt,0)=.5*dt; 

  if( debug >3 )
    printF("getIntegrationWeights: Nt=%i, dt=%9.3e\n",Nt,dt);

  if( orderOfAccuracy == 2 )
  {
    // ------- assign weights for each frequency ------
      for( int i=1; i<numFreq; i++ )
      {
        //                    ...          tf
        //  +-----+ ...   +------+------+---X---+------+  ...  -----+
        //  0                          im      ip                   Nt
        //                  dt          |   |
        //                               dt1

        Real tf = Tv(i);                        //  final time for frequency i
        int im = min( 1.*Nt,floor( tf/dt ) );   // point to left of tf
        int ip = im+1;  
        Real dt1 = tf- im*dt;                   // size of last interval
        Real alpha = dt1/dt; 
        if( alpha<0. || alpha>1. )
        {
          printF("getIntegrationWeights:ERROR: alpha<0 or alpha>1 -- something is wrong\n");
          printF(" tf=%12.5e, dt=%10.2e, im=%d, alpha=dt1/dt = %8.3f\n",tf,dt,im,alpha);
          OV_ABORT("error");
        }

        //  u(tf) = (1-alpha)*u(im) + alpha*u(ip)

        Real wim = .5*dt + .5*dt1 + .5*(1-alpha)*dt1;  // weight for point im 
        Real wip =                  .5*(  alpha)*dt1;  // weight for point ip
        if( debug >3 )
          printF("getIntegrationWeights: freq i=%i T=%9.3e, im=%i, dt1/dt=%9.3e, wim/dt=%9.3e wip/dt=%9.3e\n",i,tf,im,dt1/dt,wim/dt,wip/dt);

        assert( im!=0 );
        assert( ip<=Nt );

        sigma( 0,i) = .5*dt;  // ** assumes im!=0 
        sigma(im,i) = wim;
        sigma(ip,i) = wip;

        if( ip<Nt)
        {
          sigma(Range(ip+1,Nt),i)=0.; // zero out remaining weights
        }
      }
    }

  if( orderOfAccuracy == 4 )
  {
    // ------- assign 4th order quadrature weights for each frequency ------
    for( int i=1; i<numFreq; i++ )
    {
      //                    ...          tf
      //  +-----+ ...   +------+------+---X---+------+  ...  -----+
      //  0                  im-1    im      im+1    im+2        Nt
      //                  dt          |   |
      //                               dt1

      assert( Nt+1>9 );   // 4th order quad rule requires at least 9 time steps        

      Real tf = Tv(i);                        // final time for frequency i

      // int im = min( 1.*Nt,floor( tf/dt ) );   // point to left of tf

      int im0 = floor( tf/dt );
      int im = min( Nt, im0 );   // point to left of tf

      // // *wdh* Should be this since we fill in sigma(im+2,i) below   **CHECK ME**
      // int im0 = floor( tf/dt );    // point to left of tf
      // int im = min( Nt-2 , im0 );  // Use one-sided at right side if we are too close.

      Real dt1 = tf- im*dt;                   // size of last interval
      Real alpha = dt1/dt; 
      if( alpha<0. || alpha>1. )
      {
        printF("getIntegrationWeights:ERROR: alpha<0 or alpha>1 -- something is wrong\n");
        printF(" tf=%12.5e, dt=%10.2e, im=%d, alpha=dt1/dt = %8.3f\n",tf,dt,im,alpha);
        OV_ABORT("error");
      }

      Real t1 = im*dt;   // time at index im
      Real sf = ( tf - t1 ) / dt;   // right end of integration interval for subdomain dt1

      Real sf4 = pow(sf,4);
      Real sf3 = pow(sf,3);
      Real sf2 = pow(sf,2);

      sigma(0,i) = (1.0  / 3.0) *dt;
      sigma(1,i) = (31.0 / 24.0)*dt;
      sigma(2,i) = (5.0  / 6.0) *dt;
      sigma(3,i) = (25.0 / 24.0)*dt;

      // add periodic wrap *wdh* Jan 15, 2022
      const int imm1 = (im-1+Nt) % Nt;
      const int imp1 = (im+1   ) % Nt; 
      const int imp2 = (im+2   ) % Nt; 

      sigma(imm1,i) = ( dt/6.0) * ( 25.0/4.0  - sf4/4.0 +     sf3     - sf2              );
      sigma(im,  i) = ( dt/2.0) * (  1.0      + sf4/4.0 - 2.0*sf3/3.0 - sf2/2.0 + 2.0*sf );
      sigma(imp1,i) = (-dt/2.0) * (  1.0/12.0 + sf4/4.0 -     sf3/3.0 - sf2              );
      sigma(imp2,i) = ( dt/6.0) * (             sf4/4.0               - sf2/2.0          );

      // sigma(im-1,i) = ( dt/6.0) * ( 25.0/4.0  - sf4/4.0 +     sf3     - sf2              );
      // sigma(im,  i) = ( dt/2.0) * (  1.0      + sf4/4.0 - 2.0*sf3/3.0 - sf2/2.0 + 2.0*sf );
      // sigma(im+1,i) = (-dt/2.0) * (  1.0/12.0 + sf4/4.0 -     sf3/3.0 - sf2              );
      // sigma(im+2,i) = ( dt/6.0) * (             sf4/4.0               - sf2/2.0          );      
      
      if( (im+2)<Nt )
      {
        sigma(Range(im+3,Nt),i)=0.; // zero out remaining weights
      }
    }
  }

  // make a higher degree poly for accuracy check of quad4, eg uv=1.+ct1*tv +ct2*tv*tv for quadratic
  
  // --- test accuracy ----
  if( debug>3 )
  {
    Real tf = Tv(0); 
    RealArray tv(Nt+1);
    for( int i=0; i<=Nt; i++ )
    {
      tv(i) = i*dt; // tf/Real(Nt); 
    }

    if( orderOfAccuracy == 2 )
    {
      RealArray uv(Nt+1);
      Real ct1=1; 
      uv = 1. + ct1*tv; // integrate 1 + t 

      for( int i=0; i<numFreq; i++ )
      {
        tf = Tv(i); //  final time for frequency i
        Real uIntTrue = tf + ct1*tf*tf/2.; 

        Range N=Nt+1;
        Real uInt = sum( sigma(N,i)*uv(N) );  //  quadrature

        Real err = fabs( uInt - uIntTrue );
        printF("getIntegrationWeights: freq i=%d, T=%9.3e, error in integral(1+%g*t) = %8.2e\n",i,tf,ct1,err);

      }
    }

    if( orderOfAccuracy == 4 )
    {
      RealArray uv(Nt+1);
      Real ct1=1;
      Real ct2=0;
      Real ct3=0;
      
      uv = 1. + ct1*tv + ct2*tv*tv + ct3*tv*tv*tv; // integrate 1 + c1*t + c2*t^2 + c3*t^3 

      for( int i=0; i<numFreq; i++ )
      {
        tf = Tv(i); //  final time for frequency i
        Real uIntTrue = tf + ct1*tf*tf/2. + ct2*tf*tf*tf/3. + ct3*tf*tf*tf*tf/4.; 

        Range N=Nt+1;
        Real uInt = sum( sigma(N,i)*uv(N) );  //  quadrature

        Real err = fabs( uInt - uIntTrue );
        printF("getIntegrationWeights: freq i=%d, T=%9.3e, error in integral(1+%g*t+%g*t^2+%g*t^3) = %8.2e\n",i,tf,ct1,ct2,ct3,err);

        ct2=1.1;
        ct3=1.3;
        uv = 1. + ct1*tv + ct2*tv*tv + ct3*tv*tv*tv;
        
      }
    }
  }

  return 0;

}  


