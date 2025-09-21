//
// Compute the optimized filter parameters for the multi-frequency WaveHoltz Algorithm
//
#include "Overture.h"
#include "display.h"
#include "ParallelUtility.h"

#include "CgWave.h"

// lapack routines
#ifdef OV_USE_DOUBLE
  #define GETRF EXTERN_C_NAME(dgetrf)
  #define GETRI EXTERN_C_NAME(dgetri)
  #define GETRS EXTERN_C_NAME(dgetrs)
  #define GELSY EXTERN_C_NAME(dgelsy)
  // #define GECON EXTERN_C_NAME(dgecon)
  // #define LANGE EXTERN_C_NAME(dlange)
  // #define GEEV  EXTERN_C_NAME(dgeev)
#else
  #define GETRF EXTERN_C_NAME(sgetrf)
  #define GETRI EXTERN_C_NAME(sgetri)
  #define GETRS EXTERN_C_NAME(sgetrs)
  #define GELSY EXTERN_C_NAME(sgelsy)
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

   // find minimum norm solution to a least squares problem
   void GELSY( const int & M, const int & N, const int & NRHS, Real & A, const int & LDA, Real & B, const int & LDB, int & JPVT, Real & RCOND, int & RANK, Real & WORK, int & LWORK, int & INFO );
}


#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR_3(i1,i2,i3,I1,I2,I3) \
  for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  \
  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  \
  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  \

// Here are the versions that allow for general quadrature rules:
#define Fone(j,lambda) CgWave::idFilter( j,lambda) 
#define Fsin(j,lambda) CgWave::sinFilter( j,lambda )
#define Fcos(j,lambda) CgWave::cosFilter( j,lambda )
#define FonePrime(j,lambda) CgWave::idFilterPrime( j,lambda) 
#define FsinPrime(j,lambda) CgWave::sinFilterPrime( j,lambda )
#define FcosPrime(j,lambda) CgWave::cosFilterPrime( j,lambda )

// here are versions that use a Trapezoidal quadrature rule
#define Fone(j,lambda) 2*CgWave::sincd(lambda,T(j),dt)
#define Fsin(j,lambda) ( CgWave::coscd((omega(j)-lambda),T(j),dt) + CgWave::coscd((omega(j)+lambda),T(j),dt) )
#define Fcos(j,lambda) ( CgWave::sincd((omega(j)-lambda),T(j),dt) + CgWave::sincd((omega(j)+lambda),T(j),dt) )
#define FonePrime(j,lambda) 2*T(j)*CgWave::sincdPrime(lambda,T(j),dt)
#define FsinPrime(j,lambda) ( T(j)*( -CgWave::coscdPrime((omega(j)-lambda),T(j),dt) + CgWave::coscdPrime((omega(j)+lambda),T(j),dt) ) )
#define FcosPrime(j,lambda) ( T(j)*( -CgWave::sincdPrime((omega(j)-lambda),T(j),dt) + CgWave::sincdPrime((omega(j)+lambda),T(j),dt) ) )

// Choose functions so that free variables are a0 
#define F0(j,lambda) Fsin(j,lambda)
#define F2(j,lambda) Fcos(j,lambda)
#define F1(j,lambda) Fone(j,lambda)
#define F0Prime(j,lambda) FsinPrime(j,lambda)
#define F2Prime(j,lambda) FcosPrime(j,lambda)
#define F1Prime(j,lambda) FonePrime(j,lambda)



/* -----
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

// -------------------------------------------------------
///
/// Evaluate the z=x*T DERIVATIVE of the discrete sinc function that comes from the trapezoidal rule quadrature
///
///  sincd(x,T,dt) = sin( x*T )/( T*tan(x*dt*.5)/(.5*dt))
///                = sin( z   )/( T*tan((z/T)*dt*.5)/(.5*dt))
/// 
/// Return
///  d( sincd(x,T,dt) )/dz = (1/T)* d( sincd(x,T,dt) )/dx 
///                        = 
/// \notes See the waveHoltz paper for a derivation 
// -------------------------------------------------------
Real sincdPrime( Real x, Real T, Real dt )
{
  const Real epsilon=sqrt(REAL_EPSILON); //  what should this be ?

  Real z = x*T;
  Real xdtb2= x*dt*.5;

  Real y;
  if( fabs(z)>epsilon )
  {
    Real tanxdt = tan(xdtb2); 
    Real denom = T*tanxdt/(.5*dt); 
    Real denomPrime= T*( 1+ SQR(tanxdt) );
    y = ( ( cos(z)*T*denom - sin(z)*denomPrime )/( SQR(denom) ) )/T;

  }
  else
  {
    y = ( -(z*T)/3. - 2.*xdtb2*(dt*.5)/3. )/T; // Taylor series
  }
  return y;
}



// -------------------------------------------------------
/// \brief Evaluate the discrete cosc function
/// \notes See the waveHoltz paper for a derivation 
// -------------------------------------------------------
Real coscd( Real x, Real T, Real dt )
{
  const Real epsilon=sqrt(REAL_EPSILON); //  what should this be ?

  Real z = x*T;
  Real xdtb2= x*dt*.5;

  Real y;
  if( fabs(z)>epsilon )
  {
     y=(1.-cos(z))/(T*tan(xdtb2)/(.5*dt));
  }
  else
  {
    y = (z/2)*( 1. - z*z/12 - SQR(xdtb2)/3.);   // Taylor series
  }
  return y;
}


// -------------------------------------------------------
///
/// Evaluate the z=x*T DERIVATIVE of the discrete cosc function that comes from the trapezoidal rule quadrature
///
///  coscd(x,T,dt) = (1-cos(x*T))/( T*tan(x*dt*.5)/(.5*dt))
/// 
/// Return
///  d( coscd(x,T,dt) )/dz = (1/T)* d( coscd(x,T,dt) )/dx 
///
/// \notes See the waveHoltz paper for a derivation 
// -------------------------------------------------------
Real coscdPrime( Real x, Real T, Real dt )
{
  const Real epsilon=sqrt(REAL_EPSILON); //  what should this be ?

  Real z = x*T;
  Real xdtb2= x*dt*.5;

  Real y;
  if( fabs(z)>epsilon )
  {
    Real tanxdt = tan(xdtb2); 
    Real denom = T*tanxdt/(.5*dt);
    Real denomPrime= T*( 1. + SQR(tanxdt) );       
    y = ( ( sin(z)*T*denom - (1.-cos(z))*denomPrime )/( SQR(denom) )  )/T;

  }
  else
  {
    Real zSq = z*z; 
    y =  .5 - zSq/8. - zSq*SQR( (dt/2.) )/( SQR(T) * 2.);   // Taylor series
  }

  // printF("cosdPrime: x=%g, T=%g, z=%g, dt=%g, y=%g\n",x,T,z,dt,y);

  return y;
}
---- */

// -----------------------------------------------------------------------------------------------------
///
/// \brief Evaluate the component filter function beta for a given frequency 
///
// -----------------------------------------------------------------------------------------------------    
Real CgWave::evalBetaFunction( const Real lam, const int freq, Real dt ) const
{

  const RealArray & omega     = dbase.get<RealArray>("frequencyArrayAdjusted");
  const RealArray & T         = dbase.get<RealArray>("periodArrayAdjusted");
  // const Real & dt             = dbase.get<real>("dt"); 
  const RealArray & filterPar = dbase.get<RealArray>("filterPar"); 

  Real beta = filterPar(0,freq)*Fone(freq,lam)
            + filterPar(1,freq)*Fsin(freq,lam)
            + filterPar(2,freq)*Fcos(freq,lam);

  return beta;
}


// -----------------------------------------------------------------------------------------------------
/// \brief Compute the optimized filter parameters for the multi-frequency WaveHoltz Algorithm
///
/// \param Nlam (input/output) : input the number of lambda values to use in the least squares problem (-1=use default)
///    output the number actually used.
/// 
/// \param filterPar(0:2,numFreq) (output in dbase) :
/// \param muMin,muMax (output) : min and max value of mu(lambda) 
/// 
/// \note See research/WaveHoltz/optFilter.m for Matlab version            
// -----------------------------------------------------------------------------------------------------
int CgWave::optFilterParameters( const RealArray & frequencyArray, const RealArray & periodArray, const Real dt, int & Nlam, Real & muMin, Real & muMax )
{

  const RealArray & omega = frequencyArray;
  const RealArray & T     = periodArray;

  const int nf = frequencyArray.getLength(0);

  printF("Entering optFilterParameters, nf=%d\n",nf);
  printF("omega=[");
  for( int i=0; i<nf; i++ )
    printF("%6.2f,",omega(i));
  printF("]\n");

  int nf3 = nf*3; // total number of unknowns a0(i), a1(i), a2(i) i=1,2,...,nf
  int nf2 = nf*2;
  RealArray C1( nf2,nf2 );
  RealArray C2( nf2,nf );
  RealArray rhs(nf2);

  // --- Fill in the constraint equations ---

  // --- set mu(omega_i)=1 ---
  for( int i=0; i<nf; i++ )
  {
    rhs(i)    = 1;
    rhs(i+nf) = 0;
    for( int j=0; j<nf; j++ )
    {
      C1(i,j   ) = F0(j,omega(i));
      C1(i,j+nf) = F2(j,omega(i));
    }
    for( int j=0; j<nf; j++ )
      C2(i,j   ) = F1(j,omega(i));
       
  }
  // ---derivative equations, mu'(omega_i)=0 ---
  for( int i=0; i<nf; i++ )
  {
    for( int j=0; j<nf; j++ )
    {
      C1(i+nf,j   ) =  F0Prime(j,omega(i));
      C1(i+nf,j+nf) =  F2Prime(j,omega(i));
      // printF("F0Prime(j,omega(i))=%9.2e\n",F0Prime(j,omega(i)));
      // printF("F2Prime(j,omega(i))=%9.2e\n",F2Prime(j,omega(i)));
    }
    for( int j=0; j<nf; j++ )
    {
      C2(i+nf,j   ) = F1Prime(j,omega(i));
      // printF("F1Prime(j,omega(i))=%9.2e\n",F1Prime(j,omega(i)));
    }
  }

  // ::display(C1,"C1","%9.2e ");
  RealArray C1Factor(nf2,nf2); // holds LU factors of C1 after call below
  C1Factor=C1;  // save C1 for use below

  IntegerArray ipiv(nf2);
  int info; 
  // Factor Matrix, PA = LU 
  GETRF( nf2,nf2, C1Factor(0,0), nf2, ipiv(0), info );
  if( info!=0 )
  {
    printF("optFilterParameters:ERROR return from LU factor getrf: info=%d\n",info );
    OV_ABORT("error");
  }

  // a0 = B1*[a1,a2]^T + d
  RealArray B1( nf2,nf );
  RealArray d( nf2,1 );
  // B1 = -C1\C2;
  // C1i = inv(C1);
  // B1 = -C1i*C2;

  // Solve A x = b 
  // TRANS : "N", "T" "C"
  // Could solve many at once
  // Could form inverse and do multiplication inline (?)
  B1 = -C2;
  int nrhs=nf; // number of right-hand sides
  GETRS( "N", nf2, nrhs, C1Factor(0,0), nf2, ipiv(0), B1(0,0), nf2, info );
  if( info!=0 )
  {
    printF("optFilterParameters:ERROR return from Triangular solve (1): getrs: info=%d\n",info);
    OV_ABORT("error");
  }

  // d  =  C1\rhs;
  //  d  =  C1i*rhs;
  nrhs=1;
  d=rhs;  
  GETRS( "N", nf2, nrhs, C1Factor(0,0), nf2, ipiv(0), d(0,0), nf2, info );
  if( info!=0 )
  {
    printF("optFilterParameters:ERROR return from Triangular solve (2): getrs: info=%d\n",info);
    OV_ABORT("error");
  }


  // ::display(B1,"B1=-inv(C1)*C2","%9.2e ");
  // ::display(d,"d=inv(C1)*rhs","%9.2e");

  RealArray a0(nf);
  RealArray a1(nf);
  RealArray a2(nf);

  if( 1==0 )
  {
    // ------ CHECKS -----

    // Set values for a0(i) and then solve for a1(i) and a2(i)

    a0=-.25; 

    //  a1(:) = [-.3,-.2]; % for omega=[2,5] Npv=[1,1]

    for( int i=0; i<nf; i++ )
    {
      // a0(i) = B1(i,:)*a1 + d(i);        // a0(i) 
      // a0(i) = B1(i,:)*a1 + d(i);        // a0(i) 
      Real tmp=0;
      for( int j=0; j<nf; j++ ) tmp+= B1(i,j)*a0(j); 
      a1(i) = tmp + d(i);        
      tmp=0; 
      for( int j=0; j<nf; j++ ) tmp+= B1(i+nf,j)*a0(j); 
      a2(i) = tmp + d(i+nf);
    }

    // check residuals:
    RealArray res(nf2);
    for( int i=0; i<nf2; i++ )
    {
      res(i) = -rhs(i);
      for( int j=0; j<nf; j++ )
        res(i) = res(i) + C1(i,j)*a1(j) + C1(i,j+nf)*a2(j) + C2(i,j)*a0(j);
      printF("res(%d)=%9.2e\n",i,res(i));
    }
  } // end checks 



  // ======= LEAST SQUARE FIT FOR a0 =======
  // 
  //  mu(lam) = mv(lam)^T a0 + m0(lam)
  //
  // Form the matrix mv(lam,0:nf-1) and vector m0(lam)
  // ---------------------------------------

  Nlam = Nlam>0 ? Nlam : 1000; // default number of lambda points 
  const int N= Nlam;


  Real lambdaMin=0; 
  Real lambdaMax= omega(nf-1)*3.;   // **** DO THIS FOR NOW ****
  RealArray lamv(N);
  for( int ilam=0; ilam<N; ilam++ )
    lamv(ilam) = lambdaMin + (lambdaMax-lambdaMin)*(ilam)/(N-1.);

  RealArray mv(N,nf);
  RealArray m0(N);
  mv=0.;
  m0=0.;

  RealArray ev(nf); // holds unit vectors
  ev=0.;
  // this can be optimized 
  for( int ilam=0; ilam<N; ilam++ )
  {
    Real lam = lamv(ilam); 
    for( int i=0; i<nf; i++ )
    {
      ev(i)=1.;  // unit vector in direction i
      for( int j=0; j<nf; j++ )
      {
        mv(ilam,j)= mv(ilam,j) 
                  + B1(i   ,j)*F0(i,lam) 
                  + B1(i+nf,j)*F2(i,lam) 
                  +      ev(j)*F1(i,lam);
      }
      m0(ilam)  = m0(ilam)
                + d(i   )*F0(i,lam) 
                + d(i+nf)*F2(i,lam);
      ev(i)=0.;
    }
  }

  // ::display(m0,"m0","%8.1e ");

  // --- Solve the least squares problem:
  //   mv a0 = - m0 

  // SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )
  // 
  //  JPVT    (input/output) INTEGER array, dimension (N)
  //          On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
  //          to the front of AP, otherwise column i is a free column.
  //          On exit, if JPVT(i) = k, then the i-th column of AP
  //          was the k-th column of A.
  //  RCOND   (input) DOUBLE PRECISION
  //          RCOND is used to determine the effective rank of A, which
  //          is defined as the order of the largest leading triangular
  //          submatrix R11 in the QR factorization with pivoting of A,
  //          whose estimated condition number < 1/RCOND.


  IntegerArray jpvt(nf);
  jpvt=0;  // used as input to possible permute columns if nonzero
  Real rcond = sqrt(REAL_EPSILON)*N; // input -- used to determine RANK  -- what should this be 
  // Real rcond = REAL_EPSILON*N; // input -- used to determine RANK  -- what should this be 
  int rank;
  int MN = min(N,nf);
  int NB=2; // block size ?
  nrhs=1; 
  int lwork=max( MN+2*nf+NB*(nf+1), 2*MN+NB*nrhs );
  info=0;
  RealArray work(lwork);

  RealArray rhsm(N);
  rhsm=-m0;
  RealArray mvFactor(N,nf);
  mvFactor=mv;  // save mv for below 
  GELSY( N, nf, nrhs, mvFactor(0,0), N, rhsm(0), N, jpvt(0), rcond, rank, work(0), lwork, info );
  for( int i=0; i<nf; i++ )
    a0(i)=rhsm(i);

  printF("Return from least squares routine GELSY: rank=%d, info=%d (rcond=%9.2e)\n",rank,info,rcond);

  // ::display(a0,"a0 (opt)","%12.5e ");

  // for( int i=0; i<nf; i++ )
  //   printF("%6.2f,",omega(i));
  // printF("]\n");

  //  --- eval mu given a1(1:nf), mv(ilam,1:nf), m0(ilam) ---
  RealArray mu(N);
  mu = m0;
  for( int i=0; i<nf; i++ )
  {
    for( int ilam=0; ilam<N; ilam++ )
      mu(ilam) += mv(ilam,i)*a0(i);
  }

  muMax= max(mu);
  muMin= min(mu);
  // printF("mu[min,max] = [%8.2e,%14.7e]\n",muMin,muMax);

  // ---- Return parameters ----

  if( !dbase.has_key("filterPar") )
  {
    dbase.put<RealArray>("filterPar");
  }
  RealArray & filterPar = dbase.get<RealArray>("filterPar");

  int numPar=3;
  filterPar.redim(numPar,nf); 
  for( int i=0; i<nf; i++ )
  {
    // filterPar(i,1) = B1(i,:)*a0 + d(i);       
    // filterPar(i,1) = B1(i+nf,:)*a0 + d(i+nf);     
    Real a1=d(i), a2=d(i+nf);
    for( int j=0; j<nf; j++ )
    { 
      a1+= B1(i   ,j)*a0(j); 
      a2+= B1(i+nf,j)*a0(j); 
    }
    filterPar(0,i) = a0(i);
    filterPar(1,i) = a1;       
    filterPar(2,i) = a2; 
  }

  if( true )
  {
    int numFreq = nf;
    printF("=================================================================================\n"
           " --- optFilterPar : compute optimized filter parameters ---  \n"
           "   numFreq=%d, dt=%9.3e, Nlam=%d, muMin=%10.2e, muMax=%16.8e\n",
          numFreq,dt,Nlam,muMin,muMax );

    printF(" Filter kernel: a0 + a1 sin(omega t)+ a2 cos(omega t)\n");
    printF("  omega=[");
    for( int i=0; i<numFreq; i++ )
      printF("%6.2f,",frequencyArray(i));
    printF("]\n");  

    // const int numPar = filterPar.getLength(1);
    for( int i=0; i<numPar; i++ )
    {
      printF("  a%d(freq)=[%12.4e",i,filterPar(i,0)); 
      for( int j=1; j<numFreq; j++ ) printF(",%14.6e",filterPar(i,j)); 
      printF("];\n");
    }
    printF("=================================================================================\n");


  }

  return 0;
}