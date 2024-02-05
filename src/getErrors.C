#include "CgWave.h"
#include "OGFunction.h"
#include "ParallelUtility.h"
#include "gridFunctionNorms.h"

// ======================================================================================================
/// \brief Compute errors
// ======================================================================================================
real CgWave::
getErrors( realCompositeGridFunction & u, real t )
{
  real cpu0=getCPU();
  // printF("+++++++++++ getErrors t=%9.3e +++++++++++\n",t);

  const int & debug  = dbase.get<int>("debug");
  FILE *& debugFile  = dbase.get<FILE*>("debugFile");
  FILE *& pDebugFile = dbase.get<FILE*>("pDebugFile");

  const int & numberOfFrequencies   = dbase.get<int>("numberOfFrequencies");
  
  const RealArray & frequencyArray  = dbase.get<RealArray>("frequencyArray");
  // printF("\n >>> CgWave::getErrors: numberOfFrequencies=%d\n",numberOfFrequencies);

  real & maxError     = dbase.get<real>("maxError");      // save max-error here 
  real & solutionNorm = dbase.get<real>("solutionNorm");  // save solution norm here


  maxError=0.;

  const int & addForcing = dbase.get<int>("addForcing");
  const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
  const aString & knownSolutionOption = dbase.get<aString>("knownSolutionOption");

  bool twilightZone = addForcing && forcingOption==twilightZoneForcing;

  int & computeErrors       = dbase.get<int>("computeErrors");
  // int & plotHelmholtzErrors = dbase.get<int>("plotHelmholtzErrors");

  computeErrors = computeErrors && (twilightZone || knownSolutionOption=="userDefinedKnownSolution");
  
  //  computeErrors = twilightZone || knownSolutionOption=="userDefinedKnownSolution";
  

  if( computeErrors )
  {
    Range all;
    realCompositeGridFunction & error = dbase.get<realCompositeGridFunction>("error");
    error.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);
    for( int freq=0; freq<numberOfFrequencies; freq++ )
      error.setName(sPrintF("err%d",freq),freq);

    const int numberOfDimensions = cg.numberOfDimensions();
    Index I1,I2,I3;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg = cg[grid];
    
      // get the local serial arrays
      OV_GET_SERIAL_ARRAY(real,u[grid],uLocal);
      OV_GET_SERIAL_ARRAY(real,error[grid],errLocal);
      errLocal=0.;

      getIndex(cg[grid].dimension(),I1,I2,I3); // assign all points including ghost points.

      if( knownSolutionOption=="userDefinedKnownSolution" )
      {
        // printF("+++++++++++ getErrors for userDefinedKnownSolution +++++++++++\n");

        getUserDefinedKnownSolution( t, grid, error[grid], I1,I2,I3 ); // store true solution in error[grid]

        // ::display(error[grid],"known solution","%9.2e ");
        // ::display(uLocal,"computed solution","%9.2e ");

        const int includeGhost=1;
        bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
        if( ok )
        {
          errLocal(I1,I2,I3,all) -= uLocal(I1,I2,I3,all);
        }
      }
      else
      {
        // ----- Twilight zone ------
        assert( numberOfFrequencies==1 ); // **FIX ME**

        assert( dbase.get<OGFunction*>("tz")!=NULL );
        OGFunction & e = *dbase.get<OGFunction*>("tz");

        mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
        OV_GET_SERIAL_ARRAY(real,mg.vertex(),xLocal);


        const int includeGhost=1;
        bool ok=ParallelUtility::getLocalArrayBounds(u[grid],uLocal,I1,I2,I3,includeGhost);
        if( ok )
        {
          int numberOfComponents=1;
          Range C=numberOfComponents;
          int isRectangular=0;
          RealArray ue(I1,I2,I3);
          e.gd( ue ,xLocal,numberOfDimensions,isRectangular,0,0,0,0,I1,I2,I3,C,t);

          errLocal(I1,I2,I3) = ue(I1,I2,I3) - uLocal(I1,I2,I3);

        }
      }

      if( debug & 16 )
      {
        OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
        where( maskLocal(I1,I2,I3)==0 )
        {
          errLocal(I1,I2,I3)=0.; 
        }
        ::display(errLocal,sPrintF("\nError on grid=%d t=%9.3e",grid,t),debugFile,"%10.2e ");
      }
      
    }
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
      Real maxErrFreq = maxNorm(error,freq);
      maxError = max(maxError,maxErrFreq);
      if( numberOfFrequencies>1  )
        printF("getErrors: t=%9.3e, freq=%3d, omega=%8.3f, maxErr=%9.3e\n",t,freq,frequencyArray(freq),maxErrFreq);
    }
    // maxError = maxNorm(error);
    // printF("getErrors: t=%9.3e, maxError=%9.3e\n",t,maxError);


  }
  else if( knownSolutionOption=="userDefinedKnownSolution" )
  {
  }
  
  // compute the solution norm
  solutionNorm = maxNorm(u);
  

  timing(timeForGetError)+= getCPU()-cpu0;

  return maxError;
  


}