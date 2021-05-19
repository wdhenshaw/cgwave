#include "CgWave.h"
#include "ParallelUtility.h"

int CgWave::
printStatistics(FILE *file /* = stdout */)
//===================================================================================
// /Description:
// Output timing statistics
//
//===================================================================================
{
  fflush(0);
  Communication_Manager::Sync();
  const int & np = dbase.get<int>("np");

  const int & debug               = dbase.get<int>("debug");
  const real & tFinal             = dbase.get<real>("tFinal");
  const real & numberOfGridPoints = dbase.get<real>("numberOfGridPoints");
  const int & numberOfStepsTaken  = dbase.get<int>("numberOfStepsTaken");

  int numberOfInterpolationPoints=0;
  for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    numberOfInterpolationPoints+=cg.numberOfInterpolationPoints(grid);

  const aString & nameOfGridFile= dbase.get<aString>("nameOfGridFile");

  FILE *& logFile = dbase.get<FILE*>("logFile");
  assert( logFile!=NULL );
  
  // // output timings from Interpolant.
  // if( cg.numberOfComponentGrids()>1  && pinterpolant!=NULL )
  //   pinterpolant->printMyStatistics(logFile);
  // printMemoryUsage(logFile);
  // // printMemoryUsage(file);

  // GenericMappedGridOperators::printBoundaryConditionStatistics(logFile);
  
  for( int i=0; i<maximumNumberOfTimings; i++ )
    timing(i) =ParallelUtility::getMaxValue(timing(i),0);  // get max over processors -- results only go to processor=0

  // adjust times for waiting
  real timeWaiting=timing(timeForWaiting);
  timing(timeForAdvance) -=timeWaiting;
  timing(timeForPlotting)-=timeWaiting;
  timing(totalTime)-=timeWaiting;

  timing(totalTime)=max(timing(totalTime),REAL_MIN*10.);

  // Get max/ave times 
  RealArray maxTiming(timing.dimension(0)),minTiming(timing.dimension(0)),aveTiming(timing.dimension(0));

//   for( int i=0; i<Parameters::maximumNumberOfTimings; i++ )
//     maxTiming(i)=ParallelUtility::getMaxValue(timing(i));   // max over all processors  -- is this the right thing to do?
  
  ParallelUtility::getMaxValues(&timing(0),&maxTiming(0),maximumNumberOfTimings);
  ParallelUtility::getMinValues(&timing(0),&minTiming(0),maximumNumberOfTimings);
  ParallelUtility::getSums(&timing(0),&aveTiming(0),maximumNumberOfTimings);
  aveTiming/=np;

  real mem=Overture::getCurrentMemoryUsage();
  real maxMem=ParallelUtility::getMaxValue(mem);  // max over all processors
  real minMem=ParallelUtility::getMinValue(mem);  // min over all processors
  real totalMem=ParallelUtility::getSum(mem);  // min over all processors
  real aveMem=totalMem/np;
  real maxMemRecorded=ParallelUtility::getMaxValue(Overture::getMaximumMemoryUsage());

  const real realsPerGridPoint = (totalMem*1024.*1024.)/numberOfGridPoints/sizeof(real);

  // Get the current date
  time_t *tp= new time_t;
  time(tp);
  // tm *ptm=localtime(tp);
  const char *dateString = ctime(tp);

  for( int fileio=0; fileio<2; fileio++ )
  {
    FILE *output = fileio==0 ? logFile : file;

    fPrintF(output,"\n"
            "              -------------------CgWave Summary----------------- \n"
            "                       %s" 
            "               Grid:   %s \n" 
            "  ==== final time=%9.2e, numberOfStepsTaken =%9i, grids=%i, gridpts =%g, interp pts=%i, processors=%i ==== \n"
            "  ==== memory per-proc: [min=%g,ave=%g,max=%g](Mb), max-recorded=%g (Mb), total=%g (Mb)\n"
            "   Timings:         (ave-sec/proc:)   seconds    sec/step   sec/step/pt     %%     [max-s/proc] [min-s/proc]\n",
            dateString,(const char*)nameOfGridFile,
            tFinal,numberOfStepsTaken,cg.numberOfComponentGrids(),numberOfGridPoints,numberOfInterpolationPoints,
            np,minMem,aveMem,maxMem,maxMemRecorded,totalMem);
    
  
    int nSpace=35;
    aString dots="........................................................................";
    if( timing(0)==0. )
      timing(0)=REAL_MIN;
    for( int i=0; i<maximumNumberOfTimings; i++ )
    {
      // if( timingName[i]!="" && timing(i)>0. )    
      if( timingName[i]!="" )
        fPrintF(output,"%s%s%10.2e  %10.2e  %10.2e   %7.3f  %10.3e  %10.3e\n",(const char*)timingName[i],
                (const char*)dots(0,max(0,nSpace-timingName[i].length())),
                aveTiming(i),aveTiming(i)/numberOfStepsTaken,aveTiming(i)/numberOfStepsTaken/numberOfGridPoints,
                100.*aveTiming(i)/aveTiming(0),maxTiming(i),minTiming(i));
      
    }
  

    fPrintF(output,"--------------------------------------------------------------------------------------------------------\n");

    fPrintF(output," Memory usage: reals/grid-point = %6.2f.\n",realsPerGridPoint);
    fPrintF(output,"--------------------------------------------------------------------------------------------------------\n");

    cg.displayDistribution("CgWave",output);

    if( fileio==0 )
    {
      // Output results as a LaTeX table
      fPrintF(output,"\n\n%% ------------ Table for LaTeX -------------------------------\n");
      fPrintF(output,
              "\\begin{table}[hbt]\n"
              "\\begin{center}\\footnotesize\n"
              "\\begin{tabular}{|l|r|r|r|r|} \\hline\n"
              "  Timings:   &  seconds &    sec/step  &  sec/step/pt &  \\%%    \\\\ \\hline\n");
      for( int i=0; i<maximumNumberOfTimings; i++ )
      {
        aString name=timingName[i];
        int len=name.length();
        for( int j=0; j<len; j++)
        {
          if( name[j]==' ' ) name[j]='~';   // replace indentation space with ~
        }
        if( timingName[i]!="" && timing(i)>0. )    
          fPrintF(output,"%s\\dotfill & %10.2e & %10.2e & %10.2e & %7.3f \\\\ \n",(const char*)name,
                  timing(i),timing(i)/numberOfStepsTaken,timing(i)/numberOfStepsTaken/numberOfGridPoints,
                  100.*timing(i)/timing(0));
      
      }     
      fPrintF(output,
              " \\hline \n"
              "\\end{tabular}\n"
              "\\end{center}\n"
              "\\caption{grid=%s, %g grid points, %i interp points, %i steps taken, %i processors.}\n"
              "\\label{tab:%s}\n"
              "\\end{table}\n", (const char*)nameOfGridFile,numberOfGridPoints,numberOfInterpolationPoints,
              numberOfStepsTaken,Communication_Manager::numberOfProcessors(),(const char*)nameOfGridFile );
    }
    
    
  }

  delete tp;

  if( np>1 )
  {
    // In parallel we print some timings for each processor
    for( int fileio=0; fileio<2; fileio++ )
    {
      FILE *output = fileio==0 ? logFile : file;
      fflush(output);
      fPrintF(output,"\n"
              " ------- Summary: Timings per processor -----------\n"
              "   p   ");
      for( int i=0; i<maximumNumberOfTimings; i++ )
      { // output a short-form name (7 chars)
        if( timingName[i]!="" && maxTiming(i)!=0. )  
        {
          aString shortName="       ";
          int m=0;
          for( int s=0; m<7 && s<timingName[i].length(); s++ ) 
          { // strip off blanks
            if( timingName[i][s]!=' ' ) {shortName[m]=timingName[i][s]; m++;} //
          }
          fPrintF(output,"%7.7s ",(const char*)shortName);
        }
      }
      fPrintF(output,"\n");
      fflush(output);
      RealArray timingLocal(timing.dimension(0));
      for( int p=0; p<np; p++ )
      {
        // Note -- it did not work very well to have processor p try to write results, so instead
        // we copy results to processor 0 to print 
        timingLocal=timing;
        broadCast(timingLocal,p);  // send timing info from processor p   -- don't need a broad cast here **fix**
        fPrintF(output,"%4i : ",p);
        for( int i=0; i<maximumNumberOfTimings; i++ )
        {
          if( timingName[i]!="" && maxTiming(i)!=0. )    
            fPrintF(output,"%7.1e ",timingLocal(i));
        }
        fflush(output);
        fPrintF(output,"\n");
      }
      fPrintF(output,"\n");
      fflush(output);
    }
  
  }
  
  if( debug & 4 )
  {
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid ++ )
      cg[grid].displayComputedGeometry();
  }

  printF("\n >>>> See the file mx.log for further timings, memory usage and other statistics <<<< \n\n");
  
  
  // reset times
  timing(timeForPlotting)+=timeWaiting;
  timing(totalTime)+=timeWaiting;

  return 0;
}
