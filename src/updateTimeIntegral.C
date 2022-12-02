#include "CgWave.h"
#include "ParallelUtility.h"

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
  void GETRS( char *trans, int & n, int & nhrs, real & a, const int & lda, int & ipvt, real & b, const int & ldb, int & info );

  // compute inverse:
  void GETRI( int & n, real & a, const int & lda, const int & ipvt, real & work, const int & iwork, int & info );

  void GECON( char *norm, int & n, real & a, const int & lda, real & anorm, real & rcond, real & work, int & iwork, int & info );
  real LANGE( char *norm, int & m, int & n, real & a, const int & lda, real & work );

   void GEEV( char *jobvl, char* jobvr, int & n, real & a, const int & lda,
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
// Evaluate the WaveHoltz beta function
// -------------------------------------------------------
Real betaWaveHoltz( Real lambda, Real omega, Real T )
{
  //  General form of 
  //        beta(lambda; omega,T) = (2/T) int0^T (cos(omega*t)-.25)*cos(lambda*t) dt 
  //  No assumptions on omega or T (i.e. we do NOT assume that T=2*pi/omega))

  Real beta = mySinc( (omega-lambda)*T )  + mySinc( (omega+lambda)*T )- .5*mySinc( lambda*T );

  return beta;
}

// ----------------------------------------------------------------------------
/// \brief Return the multi-frequency Wave-Holtz matrix "A"
//----------------------------------------------------------------------------
int CgWave::getMultiFrequencyWaveHoltzMatrix( RealArray & A )
{
  const int & numberOfFrequencies   = dbase.get<int>("numberOfFrequencies");
  const RealArray & frequencyArray  = dbase.get<RealArray>("frequencyArray");
  const RealArray & periodArray     = dbase.get<RealArray>("periodArray");   

  const real & tFinal               = dbase.get<real>("tFinal");
  const real & dt                   = dbase.get<real>("dt");
  const int & adjustOmega           = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  RealArray & sigma = dbase.get<realArray>("sigma");  

  // ---- Compute the entries using quadrature rather than the exact formulae, in some cases----
  // const bool useQuadrature = adjustOmega && numberOfFrequencies>1;
  const bool useQuadrature = numberOfFrequencies>1;

  // form matrix A
  A.redim(numberOfFrequencies,numberOfFrequencies);
  A=0.;

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
        A(i,j) = sum( sigma(Range(Nt+1),i)*( (cos(omega*tv)-.25)*cos(lambda*tv) ) )* (2./T); // qudarture approximation

        Real beta = betaWaveHoltz( frequencyArray(j), frequencyArray(i), periodArray(i) );
        printF("getMultiFrequencyWaveHoltzMatrix: A(%d,%d)=%12.4e (exact)   =%12.4e (quadrature)  diff=%8.2e\n",i,j,beta,A(i,j),beta-A(i,j));

      }

    }
  }
  if( true )
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
/// \brief Update the time integral used by the Helmholtz solver
///
/// The WaveHoltz algorithm computes the time integral of the solution. 
/// This routine increments the time integral using the current solution and stores 
/// the result in "v" (from the dbase).
///
/// \param step (input) : current step 
/// \param stepOption : firstStep, middleStep, lastStep. When stepOption==firstStep, the integration weights will be comnputed and the
///     time integral will be initialized.
/// \param t (input) : current time
/// \param u (input) : current solution
// ================================================================================================
int CgWave::
updateTimeIntegral( int step, StepOptionEnum stepOption, real t, realCompositeGridFunction& u )
{

  realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");

  const int & debug                 = dbase.get<int>("debug");
  const real & tFinal               = dbase.get<real>("tFinal");
  const real & dt                   = dbase.get<real>("dt");
  const int & orderOfAccuracy       = dbase.get<int>("orderOfAccuracy");
  const int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");  
  const int & solveHelmholtz        = dbase.get<int>("solveHelmholtz");
  const int & adjustOmega           = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

  const real & omega                = dbase.get<real>("omega");
  const real & Tperiod              = dbase.get<real>("Tperiod");
  const int & numPeriods            = dbase.get<int>("numPeriods");

  int & numberOfFrequencies         = dbase.get<int>("numberOfFrequencies");
  RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
  RealArray & periodArray           = dbase.get<RealArray>("periodArray");  

  Index I1,I2,I3;


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

  if( false )
  {
   printF("updateTimeIntegral: stepOption=%d omega=%16.12e, t=%12.4e, dt=%16.12e, t/dt=%.3g\n",stepOption,omega,t,dt,t/dt);
   printF("updateTimeIntegral: numberOfFrequencies=%d\n",numberOfFrequencies);
   // u.display(sPrintF("u for updateTimeIntegral, t=%9.3e",t),"%6.2f ");
  }

  if( fabs(tFinal -periodArray(0)) > 100.*REAL_EPSILON*tFinal )
  {
    printF("CgWave::updateTimeIntegral: tFinal != periodArray(0) : tFinal=%12.4e, periodArray(0)=%12.4e, diff=%12.4e\n",
         tFinal,periodArray(0),tFinal-periodArray(0));

    OV_ABORT("error");
  }

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

  const int Nt = round( tFinal/dt ); // number of time-steps

  // sigma(0:Nt,0:numFreq-1) : integration weights
  if( !dbase.has_key("sigma") )
  {
    dbase.put<RealArray>("sigma"); // holds integration weights
  }
  RealArray & sigma = dbase.get<realArray>("sigma");  

  if( stepOption==firstStep )
  {
    assert( t==0. && step==0 );

    // --- Evaluate the integration weights for each of the frequencies ---
    getIntegrationWeights( Nt, numberOfFrequencies, periodArray, orderOfAccuracy, sigma );

  }


  if( stepOption==lastStep )
  {
    assert( step==Nt );

    if( !( fabs(t-tFinal)< REAL_EPSILON*1000.*tFinal ) )
    {
      printF("CgWave:ERROR: t != tFinal, t=%14.6e tFinal=%14.6e diff=%9.3e\n",t,tFinal,fabs(t-tFinal));
      
      OV_ABORT("ERROR");
    }
    
  }
  
  if( stepOption==firstStep )
  {
    // When solving the Helmholtz problem with CgWaveHoltz we need to evaluate an integral 
    //      v  = (1/(2*T)* Int_0^T [  ( cos(omega*t)-.25)*u(x,t) dt ] 
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
      OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
      getIndex(cg[grid].dimension(),I1,I2,I3);
      bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
      if( ok )
      {
        for( int freq=0; freq<numberOfFrequencies; freq++ )
        {
          // vLocal = ( .5*( cos(omega*(t))-.25 ) )*uLocal;  // Trapezoidal first term (.5)
          const Real omegaFreq = frequencyArray(freq); 
          vLocal(I1,I2,I3,freq) = ( sigma(step,freq)*( cos(omegaFreq*(t))-.25 ) )*uLocal;  // Trapezoidal first term (.5)
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
      bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
      if( ok )
      {
        if( true )
        {
          for( int freq=0; freq<numberOfFrequencies; freq++ )
          {
            const Real omegaFreq = frequencyArray(freq); 
            vLocal(I1,I2,I3,freq) += ( sigma(step,freq)*( cos(omegaFreq*(t))-.25 ) )*uLocal(I1,I2,I3);
          }
          if( stepOption==lastStep )
          {
            for( int freq=0; freq<numberOfFrequencies; freq++ )
            {
              vLocal(I1,I2,I3,freq) *= (2./periodArray(freq));  //  normalize the time integral  
            }        
          }
        }
        else
        {
          // -- old way 
          vLocal += ( trapFactor*( cos(omega*(t))-.25 ) )*uLocal;

          if( stepOption==lastStep )
            vLocal *= (2./tFinal)*dt;  // multiply by dt and normalize the time integral
        } 
      }
       
    } // end for grid 
    if( stepOption==lastStep )
    {
      
      if( false || (false && numberOfFrequencies==1) )  // ***************** Oct 31, 2021 NOT valid for multi-freq
        applyBoundaryConditions( v, t );


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
          bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
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
        
    }
      
  }


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
  // if( numFreq>1 && orderOfAccuracy!=2 )
  // {
  //   printF("CgWave::getIntegrationWeights:ERROR: orderOfAccuracy=%d is not supported yet.\n",orderOfAccuracy);
  //   OV_ABORT("error");
  // }

  // if( Tv(0) < max(abs(Tv)) )
  // {
  //   printF("getIntegrationWeights:ERROR: Tv(0)=%8.2e is not the largest!\n",Tv(0));
  //   ::display(Tv,"Tv","%9.2e ");
  //   OV_ABORT("error");
  // }

  // Real dt = Tv(0)/Nt; 

  // sigma.redim(Nt+1,numFreq); // note: base 0 

  // sigma = dt;  // set all weights to dt by default, this is not really needed

  // // trapezoidal rule for Tv(0) -- this is the largest period 
  // sigma( 0,0)=.5*dt; 
  // sigma(Nt,0)=.5*dt; 

  // printF("getIntegrationWeights: Nt=%i, dt=%9.3e\n",Nt,dt);

  // // ------- assign weights for each frequency ------
  // for( int i=1; i<numFreq; i++ )
  // {
  //   //                    ...          tf
  //   //  +-----+ ...   +------+------+---X---+------+  ...  -----+
  //   //  0                          im      ip                   Nt
  //   //                  dt          |   |
  //   //                               dt1

  //   Real tf = Tv(i);                        //  final time for frequency i
  //   int im = min( 1.*Nt,floor( tf/dt ) );   // point to left of tf
  //   int ip = im+1;  
  //   Real dt1 = tf- im*dt;                   // size of last interval
  //   Real alpha = dt1/dt; 
  //   if( alpha<0. || alpha>1. )
  //   {
  //     printF("getIntegrationWeights:ERROR: alpha<0 or alpha>1 -- something is wrong\n");
  //     printF(" tf=%12.5e, dt=%10.2e, im=%d, alpha=dt1/dt = %8.3f\n",tf,dt,im,alpha);
  //     OV_ABORT("error");
  //   }

  //   //  u(tf) = (1-alpha)*u(im) + alpha*u(ip)

  //   Real wim = .5*dt + .5*dt1 + .5*(1-alpha)*dt1;  // weight for point im 
  //   Real wip =                  .5*(  alpha)*dt1;  // weight for point ip
  //   printF("getIntegrationWeights: freq i=%i T=%9.3e, im=%i, dt1/dt=%9.3e, wim/dt=%9.3e wip/dt=%9.3e\n",i,tf,im,dt1/dt,wim/dt,wip/dt);

  //   assert( im!=0 );
  //   assert( ip<=Nt );

  //   sigma( 0,i) = .5*dt;  // ** assumes im!=0 
  //   sigma(im,i) = wim;
  //   sigma(ip,i) = wip;

  //   if( ip<Nt)
  //   {
  //     sigma(Range(ip+1,Nt),i)=0.; // zero out remaining weights
  //   }
  // }

  // // --- test accuracy ----
  // if( 1==1 )
  // {
  //   Real tf = Tv(0); 
  //   RealArray tv(Nt+1);
  //   for( int i=0; i<=Nt; i++ )
  //   {
  //     tv(i) = i*dt; // tf/Real(Nt); 
  //   }


  //   RealArray uv(Nt+1);
  //   Real ct1=1; 
  //   uv = 1. + ct1*tv; // integrate 1 + t 

  //   for( int i=0; i<numFreq; i++ )
  //   {
  //     tf = Tv(i); //  final time for frequency i
  //     Real uIntTrue = tf + ct1*tf*tf/2.; 

  //     Range N=Nt+1;
  //     Real uInt = sum( sigma(N,i)*uv(N) );  //  quadrature

  //     Real err = fabs( uInt - uIntTrue );
  //     printF("getIntegrationWeights: freq i=%d, T=%9.3e, error in integral(1+%g*t) = %8.2e\n",i,tf,ct1,err);

  //   }
  // }

  // return 0;

// }


