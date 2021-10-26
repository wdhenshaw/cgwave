//===============================================================================
//
//  Test quadrature formula 
//
//==============================================================================
#include "Overture.h"  
#include "display.h"



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
// =================================================================================
int getIntegrationWeights( int Nt, int numFreq, RealArray & Tv, int orderOfAccuracy, RealArray & sigma )
{
  if( orderOfAccuracy!=2 )
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

  Real dt = Tv(0)/Nt; 

  sigma.redim(Nt+1,numFreq); // note: base 0 

  sigma = dt;  // set all weights to dt by default, this is not really needed

  // trapezoidal rule for Tv(0) -- this is the largest period 
  sigma( 0,0)=.5*dt; 
  sigma(Nt,0)=.5*dt; 

  printF("getIntegrationWeights: Nt=%i, dt=%9.3e\n",Nt,dt);

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

  // --- test accuracy ----
  if( 1==1 )
  {
    Real tf = Tv(0); 
    RealArray tv(Nt+1);
    for( int i=0; i<=Nt; i++ )
    {
      tv(i) = i*dt; // tf/Real(Nt); 
    }


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

  return 0;

}


// ================================== MAIN =================================================
int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

   // This macro will initialize the PETSc solver if OVERTURE_USE_PETSC is defined.
  // *** INIT_PETSC_SOLVER(); // ** TURN OFF FOR NOW -- May 17, 2021

  printF("Usage: testQuad -numFreq=<i> -orderOfAccuracy=<i>\n");



  int debug=0; 
  int numFreq=2; 
  int orderOfAccuracy=2;

  RealArray freq(100);
  freq=1.;

  freq(1)=1.5;
  freq(2)=2.;
  freq(3)=3.5;
  freq(4)=3.7;
  freq(5)=4.1;
  freq(6)=8.234;
  freq(7)=15.78;
  
  // real fx=2., fy=2., fz=2.; // frequencies for trig TZ
  
  int len=0;
  if( argc >= 1 )
  { 
    for( int i=1; i<argc; i++ )
    {
      aString arg = argv[i];
      if( (len=arg.matches("-debug=")) )
      {
        sScanF(arg(len,arg.length()-1),"%i",&debug);
        printF("Setting debug=%i\n",debug);
      }
      else if( (len=arg.matches("-numFreq=")) )
      {
        sScanF(arg(len,arg.length()-1),"%i",&numFreq);
        printF("Setting numFreq=%i\n",numFreq);
      }
      else if( (len=arg.matches("-orderOfAccuracy=")) )
      {
        sScanF(arg(len,arg.length()-1),"%i",&orderOfAccuracy);
        printF("Setting orderOfAccuracy=%i\n",orderOfAccuracy);
      }      
      else
      {
        printF("Unknown command line argument=[%s]\n",(const char*)arg);
      }
    }
  }

  printF("=================================================================================\n"
         " --- testQuad : test the quadrature formula ---  \n"
         "   numFreq=%d, orderOfAccuracy=%d \n",
        numFreq,orderOfAccuracy);


  RealArray Tv(numFreq);  // periods
  Tv = twoPi/freq(Range(numFreq));

  int Nt=100;
  // Real dt=Tv(0)/Nt;

  RealArray sigma;

  getIntegrationWeights( Nt, numFreq, Tv, orderOfAccuracy, sigma );

  Overture::finish();          

  return(0);
}




