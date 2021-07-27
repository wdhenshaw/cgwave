#include "CgWave.h"
#include "ParallelUtility.h"


// ================================================================================================
/// \brief Update the time integral used by the Helmholtz solver
///
/// \param stepOption : firstStep, middleStep, lastStep
// ================================================================================================
int CgWave::
updateTimeIntegral( StepOptionEnum stepOption, real t, realCompositeGridFunction& u )
{

  realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");

  const int & debug      = dbase.get<int>("debug");
  const real & tFinal    = dbase.get<real>("tFinal");
  const real & dt        = dbase.get<real>("dt");

  const real & omega     = dbase.get<real>("omega");
  const real & Tperiod   = dbase.get<real>("Tperiod");
  const int & numPeriods = dbase.get<int>("numPeriods");

  Index I1,I2,I3;


  if( false )
  {
   printF("updateTimeIntegral: stepOption=%d omega=%16.12e, t=%12.4e, dt=%16.12e, t/dt=%.3g\n",stepOption,omega,t,dt,t/dt);
   u.display(sPrintF("u for updateTimeIntegral, t=%9.3e",t),"%6.2f ");
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

  if( stepOption==firstStep )
    assert( t==0. );
  if( stepOption==lastStep )
  {
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
        vLocal = ( .5*( cos(omega*(t))-.25 ) )*uLocal;  // Trapezoidal first term (.5)
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
        vLocal += ( trapFactor*( cos(omega*(t))-.25 ) )*uLocal;

        if( stepOption==lastStep )
          vLocal *= (2./tFinal)*dt;  // multiply by dt and normalize the time integral 
      }
       
    } // end for grid 
    if( stepOption==lastStep )
    {
      if( true )
        applyBoundaryConditions( v, t );

      const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");

      // u.display(sPrintF("u after integral, t=%9.3e",t),"%6.2f ");
      // v.display(sPrintF("v after integral, t=%9.3e",t),"%6.2f ");

      if( ( debug & 2 ) && knownSolutionOption=="userDefinedKnownSolution" ) 
      {
        real maxErr = getErrors( v, t );
        printF("\n ** cgWave:updateTimeIntegral: t=%9.3e, error in WaveHoltz v:  maxErr=%9.2e ** \n\n",t,maxErr);
      }
        
    }
      
  }


  return 0;
}
