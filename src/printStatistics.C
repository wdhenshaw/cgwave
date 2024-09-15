#include "CgWave.h"
#include "ParallelUtility.h"
#include "ProcessorInfo.h"
#include "Oges.h"

int CgWave::
printStatistics(FILE *file /* = stdout */)
//===================================================================================
/// \brief Output timing and memory usage statistics
///
//===================================================================================
{
  fflush(0);
  Communication_Manager::Sync();
  const int & np = dbase.get<int>("np");

  const int & debug                  = dbase.get<int>("debug");
  const real & tFinal                = dbase.get<real>("tFinal");
  const real & numberOfGridPoints    = dbase.get<real>("numberOfGridPoints");
  const int & numberOfStepsTaken     = dbase.get<int>("numberOfStepsTaken");
  const int & orderOfAccuracy        = dbase.get<int>("orderOfAccuracy");
  const int & orderOfAccuracyInTime  = dbase.get<int>("orderOfAccuracyInTime");
  const int & bcCount                = dbase.get<int>("bcCount");

  const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");
  const ModifiedEquationApproachEnum & modifiedEquationApproach  
                           = dbase.get<ModifiedEquationApproachEnum>("modifiedEquationApproach");

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
  
  const Real clockSpeed = ProcessorInfo::getClockSpeed()/1000.; // clock speed in GHz

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

  // implicit time-stepping: 
  const int & totalImplicitIterations       = dbase.get<int>("totalImplicitIterations"); 
  const int & totalImplicitSolves           = dbase.get<int>("totalImplicitSolves"); 
  const Real aveNumberOfImplicitIterations  = totalImplicitIterations/Real(max(1,totalImplicitSolves));

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
            "  ==== final time=%9.2e, numberOfStepsTaken =%9i, grids=%i, gridpts =%g, interp pts=%i, processors=%i, clock=%.2f GHz ==== \n"
            "  ==== memory per-proc: [min=%g,ave=%g,max=%g](Mb), max-recorded=%g (Mb), total=%g (Mb)\n",
            dateString,(const char*)nameOfGridFile,
            tFinal,numberOfStepsTaken,cg.numberOfComponentGrids(),numberOfGridPoints,numberOfInterpolationPoints,
            np,clockSpeed,minMem,aveMem,maxMem,maxMemRecorded,totalMem);

    if( 1==1 )
    {
      fPrintF(output," applyBoundaryConditions called %d times\n",bcCount);  // for debugging 
    }

    if( aveNumberOfImplicitIterations>0 )
    {
      fPrintF(output,"  ==== implicit: ave-iterations per solve = %5.1f (number of implicit solves=%d)\n",
        aveNumberOfImplicitIterations,totalImplicitSolves);
    }

    fPrintF(output,    
            "   Timings:         (ave-sec/proc:)   seconds    sec/step   sec/step/pt     %%     [max-s/proc] [min-s/proc]\n");
    
  
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

    // ---- output stats from the implicit solver ---
    if( timeSteppingMethod == implicitTimeStepping && dbase.has_key("impSolver") )
    {
      Oges & impSolver = dbase.get<Oges>("impSolver");
      if( impSolver.isSolverIterative() ) 
      {
        if( fileio==0 )
        { // logFile: 
          fPrintF(output,"\n");
          fPrintF(output,"-------------------------------------------------------------------------------------------------------\n");        
          fPrintF(output,"---------------------------- IMPLICIT SOLVER STATISTICS -----------------------------------------------\n");
          fPrintF(output," ==== implicitSolves: ave-iterations per solve = %5.1f (number of implicit solves=%d)\n",
                  aveNumberOfImplicitIterations,totalImplicitSolves);          
          fPrintF(output,"-------------------------------------------------------------------------------------------------------\n");        
          impSolver.printStatistics(output);
        }
      }     
    }

    // ----- Output CYCLES -------
    if( fileio==0 )
    { 
      fPrintF(output,"\n CYCLES : clock=%.2f GHz, gridpts =%.3g, numberOfStepsTaken =%g, \n",clockSpeed,numberOfGridPoints,(Real)numberOfStepsTaken);
      fPrintF(output,"   Timings:                        CYCLES/step/pt    %% \n",clockSpeed);
      for( int i=0; i<maximumNumberOfTimings; i++ )
      {
        // if( timingName[i]!="" && timing(i)>0. )    
        if( timingName[i]!="" )
        {
          Real myCycles = ( aveTiming(i)/numberOfStepsTaken/numberOfGridPoints ) *( clockSpeed * 1e9 ); // Ghz = 10^9 
          fPrintF(output,"%s%s",(const char*)timingName[i],
                  (const char*)dots(0,max(0,nSpace-timingName[i].length())));
          if( myCycles<100000 )
            fPrintF(output,"%9.1f",myCycles);
          else
            fPrintF(output,"%10.2e",myCycles);
          fPrintF(output,"      %6.2f\n",100.*aveTiming(i)/aveTiming(0));
        }
        
      }
    }
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

      int firstChar=0; 
      for( int i=nameOfGridFile.length()-1; i>=0; i-- )
      {
        if( nameOfGridFile[i]=='/' ){ firstChar=i+1; break; } // start from end, work backwards and look for a directory symbol
      }
      aString gridNameNoPrefix = nameOfGridFile(firstChar,nameOfGridFile.length()-1);

      const Real timeForAdvGrids = timing(timeForAdvanceRectangularGrids) + timing(timeForAdvanceCurvilinearGrids);
      // real timeForAdvGrids = timing(timeForAdvanceRectangularGrids) + timing(timeForAdvanceCurvilinearGrids)+ timing(timeForUp); 

      const Real nanoSecondsPerSecond=1e9; // 1 ns = 1e-9 s 

      aString blanks="                                                                 ";
      int numBlanks = max(0,min(30,gridNameNoPrefix.length()-1 -7));

      aString method;                              // modified equation time stepping 
      if( modifiedEquationApproach==standardModifiedEquation ) 
        method="ME"; 
      else if( modifiedEquationApproach==hierarchicalModifiedEquation ) 
        method = "FAME"; // factored (hierachical)
      else if( modifiedEquationApproach==stencilModifiedEquation)
        method = "stencil";
      else
        method = "unknown";

      if( orderOfAccuracyInTime==2 )    method + "T2";  // 2nd-order in time

      fPrintF(output,"\n\nSummary info for a performance table see cgwave/doc/performance.tex:\n");
      for( int ii=0; ii<2; ii++ ) // output type base on CPU or CYCLES 
      {
        aString myLabel;
        real scaleFactor = 1./numberOfStepsTaken/numberOfGridPoints*nanoSecondsPerSecond;
        if( ii==0 )
        {
          myLabel="PerfInfo";
          fPrintF(output,"            &       &     &       & \\multicolumn{3}{|c|}{CPU ns/step/pt}  \\\\\n",(const char*)blanks(0,numBlanks)); 
        }
        else
        {
          myLabel="CyclesPerfInfo";
          scaleFactor *= clockSpeed; 
          fPrintF(output,"            &       &     &       & \\multicolumn{3}{|c|}{CYCLES/step/pt}  \\\\\n",(const char*)blanks(0,numBlanks)); 
        }
        fPrintF(output,"   Grid  %s  & ord  & scheme & pts   & total   &  adv    & adv-r-c  &   up  &  bc  & interp \\\\\n",(const char*)blanks(0,numBlanks)); 
        fPrintF(output,"  %s & %d & %s  & %.2gM  & %7.1f & %7.1f & %7.1f & %7.1f & %7.1f & %7.1f   \\\\ %% %s\n",
             (const char*)gridNameNoPrefix,
             orderOfAccuracy,
             (const char*)method,
             numberOfGridPoints/1e6,
             timing(0)*scaleFactor,               // total time 
             timing(timeForAdvance)*scaleFactor,
             timeForAdvGrids*scaleFactor,         // time to advance rectangular and curvilinear grids 
             timing(timeForDissipation)*scaleFactor,
             timing(timeForBoundaryConditions)*scaleFactor,
             timing(timeForInterpolate)*scaleFactor,
             (const char*)myLabel
              );
      }

      // fPrintF(output,"   Grid    & ord &   Solver          &    pts  & its  & $\\|{\\rm res}\\|_\\infty$ & total        & setup     & solve & reals/pt \\\\ \n");

      // fPrintF(output," %s &  %d  &  %s   &   %.2gM  & %4d  &   %8.1e       & %7.2f  & %7.2f  & %7.2f  & %7.1f \\\\ % PerfInfo\n",
      //    (const char*)gridNameNoPrefix,
      //    orderOfAccuracy,
      //    (const char*)solverLabel,
      //    numberOfGridPoints/1e6,
      //    numIts,
      //    maximumResidual,
      //    tm[timeForSolve]+tm[timeForInitialize],
      //    tm[timeForInitialize],
      //    tm[timeForSolve],
      //    realsPerGridPoint
      //    );

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

  printF("\n >>>> See the file cgWave.out for further timings, memory usage and other statistics <<<< \n\n");
  
  
  // reset times
  timing(timeForPlotting)+=timeWaiting;
  timing(totalTime)+=timeWaiting;

  return 0;
}
