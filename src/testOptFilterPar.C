//===============================================================================
//
//  Test routine for computing the optimized filter parameters
//
//==============================================================================
#include "Overture.h"
#include "GL_GraphicsInterface.h"
#include "display.h"

#include "CgWave.h"

// int optFilterParameters( const RealArray & frequencyArray, const RealArray & periodArray, const Real dt, int & Nlam, 
//                          RealArray & filterPar, Real & muMin, Real & muMax );

// ================================== MAIN =================================================
int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

   // This macro will initialize the PETSc solver if OVERTURE_USE_PETSC is defined.
  // *** INIT_PETSC_SOLVER(); // ** TURN OFF FOR NOW -- May 17, 2021

  printF("Usage: testOptFilterPar -numFreq=<i> -orderOfAccuracy=<i>\n");


  
  int debug=0; 
  int numFreq=2; 
  // int orderOfAccuracy=2;

  RealArray freq(100);
  freq=1.;

  freq(0)=5.;
  freq(1)=9.;
  freq(2)=11.;
  freq(3)=17.;
  freq(4)=23.;
  
  // real fx=2., fy=2., fz=2.; // frequencies for trig TZ
  
  aString nameOfOGFile="cice2.order2.hdf";
  aString commandFileName="";
  int plotOption=false;

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
      else if( (len=arg.matches("-omega1=")) )
      {
        sScanF(arg(len,arg.length()-1),"%e",&freq(0));
        printF("Setting omega1=%12.4e\n",freq(0));
      }  
      else if( (len=arg.matches("-omega2=")) )
      {
        sScanF(arg(len,arg.length()-1),"%e",&freq(1));
        printF("Setting omega2=%12.4e\n",freq(1));
      }   
      else if( (len=arg.matches("-omega3=")) )
      {
        sScanF(arg(len,arg.length()-1),"%e",&freq(2));
        printF("Setting omega3=%12.4e\n",freq(2));
      }              
      else
      {
        printF("Unknown command line argument=[%s]\n",(const char*)arg);
      }
    }
  }

  GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("cgWave",plotOption,argc,argv));

  GraphicsParameters psp; 
  CompositeGrid cg;

  bool loadBalance=true; // turn on or off the load balancer
  getFromADataBase(cg,nameOfOGFile,loadBalance);

  CgWave cgWave(cg,ps);

  // CgWave cgWave;

  RealArray Tv(numFreq);  // periods
  Tv = twoPi/freq(Range(numFreq));

  int Nt=10;
  Real dt=Tv(0)/Nt;

  dt = 0.05; // ************* TESTING TO MATCH MATLAB 


  RealArray frequencyArray(numFreq), periodArray(numFreq);
  for( int i=0; i<numFreq; i++ )
  {
    frequencyArray(i) = freq(i);
    periodArray(i)    = Tv(i);
  }

 
  RealArray filterPar;
  Real muMin, muMax;
  int Nlam=-1; // -1 = use default 
  cgWave.optFilterParameters( frequencyArray, periodArray, dt, Nlam, muMin, muMax );



  Overture::finish();          

  return(0);
}




