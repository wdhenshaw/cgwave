// This file automatically generated from solveEigen.bC with bpp.
// ----------------------------------------------
// -------------- EIGENWAVE ---------------------
// ----  Solve for eigenpairs using WaveHoltz ---
// ----------------------------------------------


#include "CgWaveHoltz.h"
#include "CgWave.h"
// #include "Oges.h"
#include "ParallelUtility.h"
// #include "CompositeGridOperators.h"
#include "PlotStuff.h"
#include "display.h"
#include "gridFunctionNorms.h"
#include "Ogshow.h"  
#include "HDF_DataBase.h"

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  


#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )


// ============================================================================================
/// \brief Open a dialog to conveniently plot eigenvectors, residuals or errors
/// \param name (input) : name of vectors such as "EigenVector", "Error", "Residual" 
// ============================================================================================
int
CgWaveHoltz::plotEigenVectors( realCompositeGridFunction & eigenVector, const RealArray & eig, const aString & name, GL_GraphicsInterface & ps, GraphicsParameters & psp )
{
    int vector=0; // plot this eigenvector 

    const int numberOfEigenvectors = eigenVector.getComponentBound(0) - eigenVector.getComponentBase(0) + 1;

    GUIState dialog;

    dialog.setWindowTitle(name);
    dialog.setExitCommand("done", "done");

    aString cmds[] = {"next",
                                        "previous",
                                        "first",
                                        "contour",
                                        ""};

    int numberOfPushButtons=4;  // number of entries in cmds
    int numRows=(numberOfPushButtons+1)/2;
    dialog.setPushButtons( cmds, cmds, numRows ); 

    const int numberOfTextStrings=8;
    aString textLabels[numberOfTextStrings];
    aString textStrings[numberOfTextStrings];


    int nt=0;
    textLabels[nt] = "vector";  sPrintF(textStrings[nt], "%d",vector);  nt++; 

  // null strings terminal list
    textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
    dialog.setTextBoxes(textLabels, textLabels, textStrings);
    int numberOfTextBoxes=nt;

    gi.pushGUI(dialog);      

    aString answer2, buff; 
    Real t=0.; 
    for( int it=0; ; it++ )
    {
        if( it==0 )
            answer2="plot";
        else
            gi.getAnswer(answer2,"");

        if( answer2=="done" || answer2=="exit" || answer2=="finish" )
        {
            break;
        }
        else if( dialog.getTextValue(answer2,"vector","%d",vector) )
        {
            vector = max(0,min(numberOfEigenvectors-1,vector));
            printF("Setting vector=%d\n",vector);
        }       
        else if( answer2=="first" )
        {
            vector=0;
        }
        else if( answer2=="next" )
        {
            vector = ( vector +1 ) % numberOfEigenvectors;
        }
        else if( answer2=="previous" )
        {
            vector = ( vector - 1 + numberOfEigenvectors) % numberOfEigenvectors;

        }        
        gi.erase();
        if( eig.getLength(0)==numberOfEigenvectors )
        {
            psp.set(GI_TOP_LABEL,sPrintF(buff,"%s %d, eig=%.5g",(const char*)name,vector,eig(vector)));
        }
        else
        {
      // --- possible complex eigenvalue ---
            if( abs(eig(1,vector)) < 1.e-10*abs(eig(0,vector)) )
            { // Real eigenvalue
                psp.set(GI_TOP_LABEL,sPrintF(buff,"%s %d, eig=%.5g",(const char*)name,vector,eig(0,vector)));
            }
            else
            { // complex eigenvalue 
                psp.set(GI_TOP_LABEL,sPrintF(buff,"%s %d, eig=%.5g + (%.5g) I",(const char*)name,vector,eig(0,vector),eig(1,vector)));
            }
        }

        psp.set(GI_COMPONENT_FOR_CONTOURS,vector);
        if( answer2 == "contour" )
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
        else
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);


        PlotIt::contour(gi,eigenVector,psp);

    }

    gi.popGUI(); // restore the previous GUI

    return 0;
}



// ============================================================================================
/// \brief Output a table with results for EigenWave
// ============================================================================================
int CgWaveHoltz::outputEigenTable()
{
    CgWaveHoltz & cgWaveHoltz = *this;
    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");

  // -- base latex file name on matlab file name ---
    aString eigenWaveLatexTableName = "eigenWaveTable.tex";
    const aString & matlabFileName = cgWaveHoltz.dbase.get<aString>("matlabFileName"); 
    const CgWave::EigenSolverEnum & eigenSolver  = cgWave.dbase.get<CgWave::EigenSolverEnum>("eigenSolver");
    const bool useFixedPoint = eigenSolver==CgWave::fixedPointEigenSolver; 
    if( useFixedPoint )
    {
        cgWaveHoltz.dbase.get<aString>("solverName")="FixedPoint";
        eigenWaveLatexTableName = matlabFileName + "FP"; 
    }
    else
    {
        cgWaveHoltz.dbase.get<aString>("solverName")="Krylov";
        eigenWaveLatexTableName = matlabFileName + "Krylov"; 
    }
    eigenWaveLatexTableName = eigenWaveLatexTableName + ".tex";


    FILE *output = fopen((const char*)eigenWaveLatexTableName,"w" );     // Save some tex output here 

    const int & orderOfAccuracy              = cgWave.dbase.get<int>("orderOfAccuracy");
    const int & numberOfRitzVectors          = cgWave.dbase.get<int>("numberOfRitzVectors");
    const int & initialVectorsForEigenSolver = cgWave.dbase.get<int>("initialVectorsForEigenSolver");
    const int & initialVectorSmooths         = cgWave.dbase.get<int>("initialVectorSmooths");
    const int & useAccurateInnerProduct      = cgWave.dbase.get<int>("useAccurateInnerProduct");  

  // total number of steps taken:
    const int & numberOfStepsTaken           = cgWave.dbase.get<int>("numberOfStepsTaken");      
    const int & numberOfStepsPerSolve        = cgWave.dbase.get<int>("numberOfStepsPerSolve");
    const Real & cfl                         = cgWave.dbase.get<real>("cfl");
    const Real & dt                          = cgWave.dbase.get<real>("dtUsed");
    const int & minStepsPerPeriod            = cgWave.dbase.get<int>("minStepsPerPeriod");
    const RealArray & timing                 = cgWave.timing;
    const int & assignRitzFrequency          = cgWave.dbase.get<int>("assignRitzFrequency");
    const int & iterationStartRR             = cgWave.dbase.get<int>("iterationStartRR");      // start at this WH iteration

    const int & numPeriods                   = cgWaveHoltz.dbase.get<int>("numPeriods");
    const Real & maxResidual                 = cgWaveHoltz.dbase.get<real>("maxResidual");
    const real & omega                       = cgWaveHoltz.dbase.get<real>("omega");
    const int & numberOfIterations           = cgWaveHoltz.dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken

    const CgWave::TimeSteppingMethodEnum & timeSteppingMethod = cgWave.dbase.get<CgWave::TimeSteppingMethodEnum>("timeSteppingMethod");
    aString implicitSolverName =  timeSteppingMethod==CgWave::explicitTimeStepping ? 
                                                                "none" : 
                                                                cgWave.dbase.get<aString>("implicitSolverName");

    const Real aveItsPerImplicitSolve               = cgWave.getAverageNumberOfIterationsPerImplicitSolve();
    const Real & waveHoltzAsymptoticConvergenceRate = cgWave.dbase.get<real>("waveHoltzAsymptoticConvergenceRate");

    const Real cpuWave   = timing(CgWave::totalTime);
    const Real cpuAdvance = timing(CgWave::timeForAdvance);

    const Real cpuSolveEigenWave = cgWaveHoltz.dbase.get<Real>("cpuSolveEigenWave");


    aString & nameOfGridFile = cgWave.dbase.get<aString>("nameOfGridFile");
    int firstChar=0; 
    for( int i=nameOfGridFile.length()-1; i>=0; i-- )
    {
        if( nameOfGridFile[i]=='/' ){ firstChar=i+1; break; } // start from end, work backwards and look for a directory symbol
    }
    int lastChar=firstChar; 
    for( int i=firstChar; i<=nameOfGridFile.length()-1; i++ )
    {
        if( nameOfGridFile[i]=='.' ){ lastChar=i-1; break; } // remove suffix: .order2.hdf
    }

    aString gridNameNoPrefix = nameOfGridFile(firstChar,lastChar);

  // const CgWave::TimeSteppingMethodEnum & timeSteppingMethod = cgWave.dbase.get<CgWave::TimeSteppingMethodEnum>("timeSteppingMethod");
  // int & minStepsPerPeriod              = cgWave.dbase.get<int>("minStepsPerPeriod");
  // const CgWave::EigenSolverEnum & eigenSolver  = cgWave.dbase.get<CgWave::EigenSolverEnum>("eigenSolver"); 

    const int & numEigenVectors = dbase.get<int>("numEigenVectors");

  // int setInitialConditions = -1;
    const Real & tol                                = dbase.get<Real>("tol");
    const int & numberOfMatrixVectorMultiplications = dbase.get<int>("numberOfMatrixVectorMultiplications");
    const int & numEigsToCompute                    = cgWave.dbase.get<int>("numEigsToCompute");
    const int & numArnoldiVectors                   = cgWave.dbase.get<int>("numArnoldiVectors");
    const Real & eigenValueTolForMultiplicity       = cgWave.dbase.get<Real>("eigenValueTolForMultiplicity");
  // int iteration=-1;
  // int nev=-1, nconv=-1;

    printF("\n");
    printF(" -----------------------------------------------------------------------------\n");
    printF(" ------------------------- EigenWave SUMMARY ---------------------------------\n");
    printF(" omega=%12.4e, orderOfAccuracy=%d, grid=%s\n",omega,orderOfAccuracy,(const char*)gridNameNoPrefix);
    aString eigenSolverName = (eigenSolver==CgWave::defaultEigenSolver           ? "KrylovSchur" :
                                                          eigenSolver==CgWave::KrylovSchurEigenSolver       ? "KrylovSchur" :
                                                          eigenSolver==CgWave::ArnoldiEigenSolver           ? "Arnoldi" : 
                                                          eigenSolver==CgWave::ArpackEigenSolver            ? "Arpack" : 
                                                          eigenSolver==CgWave::fixedPointEigenSolver        ? "fixedPoint" : 
                                                          eigenSolver==CgWave::powerEigenSolver             ? "power" :
                                                          eigenSolver==CgWave::inverseIterationEigenSolver  ? "inverseIteration" :
                                                          eigenSolver==CgWave::JacobiDavidsonEigenSolver    ? "JacobiDavidson" :
                                                          eigenSolver==CgWave::subspaceIterationEigenSolver ? "subspaceIteration" :
                                                          "unknown" );
    aString eigenSolverNameLong = eigenSolverName; 
    if( eigenSolver !=CgWave::fixedPointEigenSolver )
    {
        if( eigenSolver==CgWave::ArpackEigenSolver )
            eigenSolverNameLong += " (ARPACK solver)";
        else
            eigenSolverNameLong +=" (SLEPc solver)";
    }

    printF(" eigenSolver = %s\n",(const char*)eigenSolverNameLong);
  
    printF(" Provide initial vector(s) for eigenSolver = %d, eigenValueTolForMultiplicity=%g\n",initialVectorsForEigenSolver,eigenValueTolForMultiplicity);
    printF(" Use accurate inner product = %d.\n",useAccurateInnerProduct);
  // printF(" timeSteppingMethod = %s, minStepsPerPeriod =%d, numPeriods=%d\n",
  //                                     (timeSteppingMethod==CgWave::explicitTimeStepping ? "explicit (modified equation)" :
  //                                      timeSteppingMethod==CgWave::implicitTimeStepping ? "implicit" : "unknown"),
  //                                      minStepsPerPeriod,numPeriods );   
 // const Real aveItsPerImplicitSolve = cgWave.getAverageNumberOfIterationsPerImplicitSolve();
 // const int & numberOfStepsTaken    = cgWave.dbase.get<int>("numberOfStepsTaken");      
 // const int & numberOfStepsPerSolve = cgWave.dbase.get<int>("numberOfStepsPerSolve");
 // const Real & cfl                  = cgWave.dbase.get<real>("cfl");
 // const Real & dt                   = cgWave.dbase.get<real>("dtUsed");
    aString timeSteppingName = (timeSteppingMethod==CgWave::explicitTimeStepping ? "explicit" : "implicit");
    printF(" Time stepping = %s, cfl=%6.2f, dt=%9.3e, minStepsPerPeriod=%d, Np=%d, numberOfStepsPerSolve=%d, \n  totalTimeStepsStaken=%d, aveItsPerImplicitSolve=%5.1f\n",
                (const char*)timeSteppingName,
                cfl,dt,minStepsPerPeriod,numPeriods,numberOfStepsPerSolve,numberOfStepsTaken,aveItsPerImplicitSolve);   
    if( timeSteppingMethod==CgWave::implicitTimeStepping )
    {
        if( cgWave.dbase.has_key("implicitSolverParameters") )
        {
            OgesParameters & par = cgWave.dbase.get<OgesParameters>("implicitSolverParameters");
            printF(" implict solver=%s\n",(const char*)par.getSolverName());
        
            Real atol, rtol;
            par.get(OgesParameters::THEabsoluteTolerance,atol);
            par.get(OgesParameters::THErelativeTolerance,rtol);
            printF(" implicit solver: rtol=%9.2e, atol=%9.2e\n",rtol,atol);
        } 
    }
    const Real numTimeStepsPerEig = numberOfStepsTaken/max(Real(numEigenVectors),1.);

    printF(" tol=%9.2e, number of wave-solves (mat-vecs) = %d, cpu=%9.2e(s)\n",tol,numberOfMatrixVectorMultiplications,cpuSolveEigenWave);
    const Real numWaveSolvesPerEig = numberOfMatrixVectorMultiplications/max(Real(numEigenVectors),1.);
    printF(" num eigs requested=%d, number eigs converged =%d, numArnoldiVectors=%d, wave-solves/eig=%5.1f, timeSteps/eig=%5.1f\n",
        numEigsToCompute,numEigenVectors,numArnoldiVectors,numWaveSolvesPerEig,numTimeStepsPerEig);
    printF(" numberOfIterations=%d, iterationStartRR=%d, assignRitzFrequency=%d, numberOfRitzVectors=%d (for fixed-point)\n",
                numberOfIterations,iterationStartRR,assignRitzFrequency,numberOfRitzVectors);
    printF(" -----------------------------------------------------------------------------\n");

    if( numEigenVectors>0 )
    {

        realCompositeGridFunction & eigenVector = dbase.get<realCompositeGridFunction>("eigenVectorRayleighRitz");

        RealArray & eigenValues = dbase.get<RealArray>("eigenValues");  // best estimate of eigenvalues
        
        realCompositeGridFunction & res = cgWaveHoltz.dbase.get<realCompositeGridFunction>("residual");
        Range all;
        res.updateToMatchGrid(cg,all,all,all,numEigenVectors);
        res.setName("evectRes");

    // // -- sort the eigenvalues:
    // IntegerArray iperm(numEigenVectors);
    // RealArray eigenValuesSorted;
    // eigenValuesSorted = eigenValues;
    // sortArray( eigenValuesSorted, iperm );


        bool saveErrors=false;
        RealArray lambda(numEigenVectors), lambdaTrue(numEigenVectors);

        if( !dbase.has_key("relErrEigenvalue") )
        {
            dbase.put<RealArray>("relErrEigenvalue");
            dbase.put<RealArray>("relErrEigenvector");
            dbase.put<RealArray>("eigenPairResidual");

        }
        RealArray & relErrEigenvalue  = dbase.get<RealArray>("relErrEigenvalue");
        RealArray & relErrEigenvector = dbase.get<RealArray>("relErrEigenvector");
        RealArray & eigenPairResidual = dbase.get<RealArray>("eigenPairResidual");

        relErrEigenvalue.redim(numEigenVectors);
        relErrEigenvector.redim(numEigenVectors);
        eigenPairResidual.redim(numEigenVectors);

        for( int ie=0; ie<numEigenVectors; ie++ ) 
        {
            res.setName(sPrintF("res%d",ie),ie);

            Real lamRQ = eigenValues(ie);
            eigenPairResidual(ie) = cgWave.getEigenPairResidual( lamRQ, eigenVector, res, ie );
            lambda(ie) = lamRQ;  // may be over-written below 

      // printF(" i=%d: lamRQ = %16.10e,  rel-resid = || L v + lamRQ^2 v ||/lamRQ^2 = %9.2e\n",ie,lamRQ,eigenPairResidual(ie));
        }

        if( cgWave.dbase.has_key("uev") )
        {
      // --- True values are known ----
            IntegerArray multipleEigIndex(numEigenVectors);

            if( !dbase.has_key("eigIndex") )
            {
                dbase.put<IntegerArray>("eigIndex");
            }
            IntegerArray & eigIndex = dbase.get<IntegerArray>("eigIndex");
            eigIndex.redim(numEigenVectors);

      // ------ check errrors in with "true" discrete eigenvalues -----
            saveErrors=true; // save errors in cgWave.dbase.get<realCompositeGridFunction>("error")

            RealArray eigRelErr(numEigenVectors);

      // Here are the "true" eigenvectors and eigenvalues:
            realCompositeGridFunction & uev = cgWave.dbase.get<realCompositeGridFunction>("uev");
            RealArray & eigTrue             = cgWave.dbase.get<RealArray>("eig");
            IntegerArray & eigMultiplicity  = cgWave.dbase.get<IntegerArray>("eigMultiplicity");

            const int numberOfEigenvectorsTrue = uev.getComponentBound(0) - uev.getComponentBase(0) + 1;
      // printF(">> There are %d eigenvectors (true, from file).\n",numberOfEigenvectorsTrue);

            if( saveErrors )
            {
                realCompositeGridFunction & error = cgWave.dbase.get<realCompositeGridFunction>("error");
                int numErrComp = error.getComponentBound(0)- error.getComponentBase(0)+1;
                if( numErrComp != numEigenVectors || 
                        error.numberOfComponentGrids()<1 )
                {
          // printF("\n ++++++ REDIMENSION ERROR GF ++++++++\n\n");
                    error.updateToMatchGrid( cg,all,all,all,numEigenVectors );
                }
                for( int ie=0; ie<numEigenVectors; ie++ )
                    error.setName(sPrintF("err%d",ie),ie);    
            }

            printF("\n");
            printF("  ie     lambda      [eig]   lamTrue     mult  eig-err  evect-err  eig-res\n");
            printF(" .............................................................................\n");
            for( int ie=0; ie<numEigenVectors; ie++ )
            {
        // int ie = iperm(je);    // sorted order
        // int ie = je; // not sorted 
                lambda(ie) = eigenValues(ie);

        // Real diffMin = REAL_MAX*.1;
                eigIndex(ie)=0;  

        // Real lambdaTrue, relErrEigenvalue, relErrEigenvector;

                cgWave.getErrorInEigenPair( lambda(ie), eigenVector, ie, lambdaTrue(ie), relErrEigenvalue(ie), relErrEigenvector(ie), 
                                                                        eigIndex(ie), multipleEigIndex(ie), saveErrors ); 

                printF(" %2d  %14.8e  [%2d]=%14.8e   %d %9.2e %9.2e %9.2e\n",
                                    ie,lambda(ie),eigIndex(ie),lambdaTrue(ie),eigMultiplicity(eigIndex(ie)), relErrEigenvalue(ie), relErrEigenvector(ie),eigenPairResidual(ie));
        // printF(" ie=%2d: lamRQ=%16.10e, true[%2d]=%16.10e, mult=%d, eig-relerr = %9.2e, eig-vector-relerr = %9.2e\n",
        //           ie,lambda,eigIndex,lambdaTrue,eigMultiplicity(eigIndex), relErrEigenvalue, relErrEigenvector);
            }
            printf(" .............................................................................\n");
            printF(" Note:  relative-errors, eig-res = || Lv+lamRQ^2 v ||/lamRQ^2\n");

        }
        else
        {
            for( int ie=0; ie<numEigenVectors; ie++ ) 
            {
                Real lamRQ = eigenValues(ie);
        // eigenPairResidual(ie) = cgWave.getEigenPairResidual( lamRQ, eigenVector, ie );

                printF(" i=%d: lamRQ = %16.10e,  rel-resid = || L v + lamRQ^2 v ||/lamRQ^2 = %9.2e\n",ie,lamRQ,eigenPairResidual(ie));
            }
        }



    // Output results as a LaTeX table
        fPrintF(output,"\n\n%% ------------ EigenWave Table for LaTeX -------------------------------\n");

        const Real maxEigErr    = max(relErrEigenvalue);
        const Real maxEvectErr  = max(relErrEigenvector);
        const Real maxEigResid  = max(eigenPairResidual);

        if( saveErrors )
        {  
            IntegerArray & eigMultiplicity  = cgWave.dbase.get<IntegerArray>("eigMultiplicity");
            IntegerArray & eigIndex         = dbase.get<IntegerArray>("eigIndex");
            fPrintF(output,"\\newcommand{\\eigenWaveLongTable}{%% start of LongTable for %s\n",(const char*)gridNameNoPrefix);
            fPrintF(output,
                        "\\begin{table}[hbt]\n"
                        "\\begin{center}\\tableFontSize\n"
                        "\\begin{tabular}{|c|r|r|c|c|c|c|c|} \\hline\n"
                        " \\multicolumn{8}{|c|}{EigenWave: grid=%s, ts=%s, order=%d, $\\omega=%6.1f$, $N_p=%d$, %s } \\\\ \\hline \n"
                        "   $j$  & \\multicolumn{1}{c|}{$\\lambda_j$} &  \\multicolumn{1}{c|}{$\\lambda_{h,k}^{\\rm true}$}  & $k$ &  mult &  eig-err & evect-err & eig-res  \\\\ \\hline\n",
                        (const char*)gridNameNoPrefix,(const char*)timeSteppingName, orderOfAccuracy, omega, numPeriods, (const char*)eigenSolverName );
            for( int ie=0; ie<numEigenVectors; ie++ )
            {
                fPrintF(output," %3d   & %10.6f & %10.6f &  %4d & %4d  & %9.2e & %9.2e & %9.2e \\\\ \n",
                                    ie,lambda(ie),lambdaTrue(ie),eigIndex(ie),eigMultiplicity(eigIndex(ie)), relErrEigenvalue(ie), relErrEigenvector(ie),eigenPairResidual(ie));
            }

            fPrintF(output,
                        " \\hline \n"
                        "\\end{tabular}\n"
                        "\\end{center}\n"
                        "\\caption{grid=%s, method=%s, ts=%s, order=%d, $\\omega=%6.1f$, $N_p=%d$.}\n"
                        "\\label{tab:%sOrder%d}\n"
                        "\\end{table}\n", (const char*)gridNameNoPrefix,(const char*)eigenSolverName, (const char*)timeSteppingName, 
                        orderOfAccuracy, omega, numPeriods, (const char*)gridNameNoPrefix,orderOfAccuracy );

            fPrintF(output,"}%% end of LongTable for %s\n",(const char*)gridNameNoPrefix);

      // ---- summary table ----
      // Real maxEigErr    = max(relErrEigenvalue);
      // Real maxEvectErr  = max(relErrEigenvector);
      // Real maxEigResid  = max(eigenPairResidual);

            fPrintF(output,"\n");
      // fPrintF(output,
      //      "\\begin{table}[hbt]\n"
      //       "\\begin{center}\\tableFontSize\n"
      //       "\\begin{tabular}{|c|c|c|c|c|c|c|} \\hline\n"
      //       "   eig-err  & evect-err & eig-res   &  num      &  total        & wave-solve &  time-steps \\\\\n"
      //       "            &           &           &  eigs     &  wave-solves  &  per eig    &  per-eig \\\\ \\hline\n");
      // fPrintF(output," %9.2e & % 9.2e & %9.2e &   %d  &  %d    &  %d    &  %d \\\\\n",
      //      maxEigErr,maxEvectErr,maxEigResid,numEigenVectors,numberOfMatrixVectorMultiplications,
      //      (int)round(numberOfMatrixVectorMultiplications/numEigenVectors), 
      //      (int)round(numberOfStepsTaken/numEigenVectors) );

            const Real waveSolvesPerEig = 1.*numberOfMatrixVectorMultiplications/numEigenVectors;
            fPrintF(output,"\\newcommand{\\eigenWaveSummaryTable}{%% start of summary table for %s\n",(const char*)gridNameNoPrefix);
            fPrintF(output,
                      "\\begin{table}[hbt]\n"
                        "\\begin{center}\\tableFontSize\n"
                        "\\begin{tabular}{|c|c|c|c|c|c|c|c|} \\hline\n"
                        " \\multicolumn{8}{|c|}{EigenWave: grid=%s, ts=%s, order=%d, $\\omega=%6.1f$, $N_p=%d$, %s } \\\\ \\hline \n"
                        "   num   &  wave       & time-steps  & wave-solves &  time-steps &   max      &  max      &  max       \\\\ \n"
                        "   eigs  &  solves     & per period  &  per eig    &  per-eig    &   eig-err  & evect-err & eig-res    \\\\ \n \\hline\n",
                          (const char*)gridNameNoPrefix,(const char*)timeSteppingName, orderOfAccuracy, omega, numPeriods, (const char*)eigenSolverName );
            fPrintF(output,"    %d    &       %d     &      %d     &    %3.1f    &     %d     & %9.2e & % 9.2e & %9.2e \\\\\n",
                      numEigenVectors,
                      numberOfMatrixVectorMultiplications,
                      minStepsPerPeriod,
           // (int)round(numberOfMatrixVectorMultiplications/numEigenVectors), 
                      waveSolvesPerEig,
                      (int)round(numberOfStepsTaken/numEigenVectors),
                      maxEigErr,maxEvectErr,maxEigResid
                      );

            fPrintF(output,
                        " \\hline \n"
                        "\\end{tabular}\n"
                        "\\end{center}\n"
                        "\\vspace*{-1\\baselineskip}\n"
                        "\\caption{EigenWave: grid=%s, method=%s, ts=%s, order=%d, $\\omega=%6.1f$, $N_p=%d$.}\n"
                        "\\label{tab:%sOrder%dSummary}\n"
                        "\\end{table}\n",(const char*)gridNameNoPrefix,(const char*)eigenSolverName, (const char*)timeSteppingName, orderOfAccuracy, 
                        omega, numPeriods,(const char*)gridNameNoPrefix,orderOfAccuracy );    
            fPrintF(output,"}%% end of summary table for %s\n",(const char*)gridNameNoPrefix);

        }
        else
        {
      // ------------ just save residuals in the table ---------

            fPrintF(output,"\\newcommand{\\eigenWaveLongTable}{%% start of LongTable for %s\n",(const char*)gridNameNoPrefix);
            fPrintF(output,
                        "\\begin{table}[hbt]\n"
                        "\\begin{center}\\tableFontSize\n"
                        "\\begin{tabular}{|c|r|c|} \\hline\n"
                        " \\multicolumn{3}{|c|}{EigenWave: grid=%s, ts=%s, order=%d, $\\omega=%6.1f$, $N_p=%d$, %s } \\\\ \\hline \n"
                        "   $j$  & \\multicolumn{1}{c|}{$\\lambda_j$} &  eig-res  \\\\ \\hline\n",
                        (const char*)gridNameNoPrefix,(const char*)timeSteppingName, orderOfAccuracy, omega, numPeriods, (const char*)eigenSolverName );
            for( int ie=0; ie<numEigenVectors; ie++ )
            {
        // fPrintF(output," %3d   & %10.6f & %10.6f &  %4d & %4d  & %9.2e & %9.2e & %9.2e \\\\ \n",
        //           ie,lambda(ie),lambdaTrue(ie),eigIndex(ie),eigMultiplicity(eigIndex(ie)), relErrEigenvalue(ie), relErrEigenvector(ie),eigenPairResidual(ie));
                fPrintF(output," %4d & %10.6f  & %9.2e \\\\ \n",
                                    ie,lambda(ie),eigenPairResidual(ie));      
            }

            fPrintF(output,
                        " \\hline \n"
                        "\\end{tabular}\n"
                        "\\end{center}\n"
                        "\\caption{grid=%s, method=%s, ts=%s, order=%d, $\\omega=%6.1f$, $N_p=%d$.}\n"
                        "\\label{tab:%sOrder%d}\n"
                        "\\end{table}\n", (const char*)gridNameNoPrefix,(const char*)eigenSolverName, (const char*)timeSteppingName, 
                        orderOfAccuracy, omega, numPeriods, (const char*)gridNameNoPrefix,orderOfAccuracy );

            fPrintF(output,"}%% end of LongTable for %s\n",(const char*)gridNameNoPrefix);

      // -----------------------------------------------
      // ---------------- summary table ----------------
      // -----------------------------------------------

            fPrintF(output,"\n");

            const Real waveSolvesPerEig = 1.*numberOfMatrixVectorMultiplications/numEigenVectors;
            fPrintF(output,"\\newcommand{\\eigenWaveSummaryTable}{%% start of summary table for %s\n",(const char*)gridNameNoPrefix);
            fPrintF(output,
                      "\\begin{table}[hbt]\n"
                        "\\begin{center}\\tableFontSize\n"
                        "\\begin{tabular}{|c|c|c|c|c|c|} \\hline\n"
                        " \\multicolumn{6}{|c|}{EigenWave: grid=%s, ts=%s, order=%d, $\\omega=%6.1f$, $N_p=%d$, %s } \\\\ \\hline \n"
                        "   num   &  wave       & time-steps  & wave-solves &  time-steps  &  max       \\\\ \n"
                        "   eigs  &  solves     & per period  &  per eig    &  per-eig     & eig-res    \\\\ \n \\hline\n",
                          (const char*)gridNameNoPrefix,(const char*)timeSteppingName, orderOfAccuracy, omega, numPeriods, (const char*)eigenSolverName );
            fPrintF(output,"    %d   &   %d    &      %d     &    %3.1f    &     %d     &  %9.2e \\\\\n",
                      numEigenVectors,
                      numberOfMatrixVectorMultiplications,
                      minStepsPerPeriod,
                      waveSolvesPerEig,
                      (int)round(numberOfStepsTaken/numEigenVectors),
                      maxEigResid
                      );

            fPrintF(output,
                        " \\hline \n"
                        "\\end{tabular}\n"
                        "\\end{center}\n"
                        "\\vspace*{-1\\baselineskip}\n"
                        "\\caption{EigenWave: grid=%s, method=%s, ts=%s, order=%d, $\\omega=%6.1f$, $N_p=%d$.}\n"
                        "\\label{tab:%sOrder%dSummary}\n"
                        "\\end{table}\n",(const char*)gridNameNoPrefix,(const char*)eigenSolverName, (const char*)timeSteppingName, orderOfAccuracy, 
                        omega, numPeriods,(const char*)gridNameNoPrefix,orderOfAccuracy );    
            fPrintF(output,"}%% end of summary table for %s\n",(const char*)gridNameNoPrefix);



        }
        fclose(output);
        printF("Wrote results to LaTeX file [%s]\n",(const char*)eigenWaveLatexTableName);


  // --- save a check file for EigenWave ---
        FILE *checkFile=NULL;
        checkFile = fopen("eigenWave.check","w" );        // for regression and convergence tests

        assert( checkFile != NULL );

    // Get the current date
        time_t *tp= new time_t;
        time(tp);
        const char *dateString = ctime(tp);
        fPrintF(checkFile,"# Check file for EigenWave, grid=%s, %s",(const char*)gridNameNoPrefix,dateString);  // Note: dateString include newline
        delete tp; 

        fPrintF(checkFile,"grid=%s;\n",(const char*)gridNameNoPrefix);
        fPrintF(checkFile,"eigenSolver=%s;\n",(const char*)eigenSolverName);
        fPrintF(checkFile,"timeStepping=%s;\n",(const char*)timeSteppingName);
        fPrintF(checkFile,"orderOfAccuracy=%d;\n",orderOfAccuracy);
        fPrintF(checkFile,"numPeriods=%d;\n",numPeriods);
        fPrintF(checkFile,"numberOfStepsPerSolve=%d;\n",numberOfStepsPerSolve);
        fPrintF(checkFile,"numEigsRequested=%d;\n",numEigsToCompute);
        fPrintF(checkFile,"numEigsComputed=%d;\n",numEigenVectors);
        fPrintF(checkFile,"numArnoldiVectors=%d;\n",numArnoldiVectors);
        fPrintF(checkFile,"numWaveSolves=%d;\n",numberOfMatrixVectorMultiplications);
        fPrintF(checkFile,"maxEigErr=%9.2e;\n",maxEigErr);
        fPrintF(checkFile,"maxEvectErr=%9.2e;\n",maxEvectErr);
        fPrintF(checkFile,"maxEigResid=%9.2e;\n",maxEigResid);
        fclose(checkFile);

        printF("Wrote results to the check file [eigenWave.check]\n");

    } // end numEigenVectors>0 

    return 0;

}


// ============================================================================================
/// \brief Assign the initial condition used by SLEPc eigen solvers
// ============================================================================================
int CgWaveHoltz::assignEigenSolverInitialCondition( bool smoothInitialCondition )
{
    CgWaveHoltz & cgWaveHoltz = *this;
    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
    CgWave::EigenSolverInitialConditionEnum & eigenSolverInitialCondition 
                                                                                      = cgWave.dbase.get<CgWave::EigenSolverInitialConditionEnum>("eigenSolverInitialCondition");
    const int numberOfDimensions = cg.numberOfDimensions();

    if( eigenSolverInitialCondition == CgWave::defaultEigenSolverInitialCondition )
    {
        printF("assignEigenSolverInitialCondition: Setting eigen solver initial condition to use the default in SLEPSc.\n");
    }    
    else if( eigenSolverInitialCondition == CgWave::randomEigenSolverInitialCondition )
    {
        printF("assignEigenSolverInitialCondition: Set RANDOM initial conditions for SLEPc\n");

        realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
        CompositeGrid & cg = *v.getCompositeGrid();
        Index I1,I2,I3;

        std::srand(12789.);

        const Real scale=1.; // 1.e-5; // TEST
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

            vLocal=0.;
            getIndex(mg.gridIndexRange(),I1,I2,I3);
            FOR_3D(i1,i2,i3,I1,I2,I3)
            {
        // do this for now
        // vLocal(i1,i2,i3) = sin(i1)*cos(i2);
                if( maskLocal(i1,i2,i3)!=0 )
                {
                    vLocal(i1,i2,i3) = (-1. + std::rand()*(2./RAND_MAX))*scale; // [-1,1]
                }

            }
        }
    // replot=true;
    // cgWave.resetTimings(); // reset CPU timings to zero 
    }  

    else if( eigenSolverInitialCondition == CgWave::sineEigenSolverInitialCondition )
    {
    // eigenSolverInitialCondition = CgWave::sineEigenSolverInitialCondition;
        printF("assignEigenSolverInitialCondition: Set SINE initial conditions for SLEPc\n");

        realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
        CompositeGrid & cg = *v.getCompositeGrid();
        Index I1,I2,I3;

    // ---- Base frequency of the sine solution on the target frequency ----
        const Real & c     = cgWave.dbase.get<Real>("c");
        const Real & omega = cgWaveHoltz.dbase.get<real>("omega");
        const Real k = omega/c; 
        const Real theta = .5*(twoPi/4.); // Pi/4 

        const Real kxHat = cos(theta);
        const Real kyHat = sin(theta);

        const Real kx = kxHat*k;
        const Real ky = kyHat*k;
        const Real kz = kx;      
        
        const Real scale=1.; // should not matter
        printF("solveEigen: Set initial condition sin(kx*x)*sin(ky*y) [*sin(kz*z)] kx=%9.2e ky=%9.2e [kz=%8.2e]\n",kx,ky,kz);
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
            OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
            mg .update(MappedGrid::THEvertex);
            OV_GET_SERIAL_ARRAY( Real, mg.vertex(),xLocal);        

            vLocal=0.;
            getIndex(mg.gridIndexRange(),I1,I2,I3);
            int includeParallelGhost=1;
            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeParallelGhost);      
            if( ok )
            {
                if( numberOfDimensions==2 )
                    vLocal(I1,I2,I3) = scale*sin( kx*xLocal(I1,I2,I3,0) ) * sin( ky*xLocal(I1,I2,I3,1) );
                else
                    vLocal(I1,I2,I3) = scale*sin( kx*xLocal(I1,I2,I3,0) ) * sin( ky*xLocal(I1,I2,I3,1) )* sin( kz*xLocal(I1,I2,I3,2) );
            }


        }


    }

    else if( eigenSolverInitialCondition == CgWave::sumOfEigenvectorsInitialCondition )
    {

        realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
        CompositeGrid & cg = *v.getCompositeGrid();

    // --- get eigenvectors
        cgWave.initializeDeflation();

    // Here are the eigenvectors and eigenvalues:
        realCompositeGridFunction & uev = cgWave.dbase.get<realCompositeGridFunction>("uev");
        RealArray & eig                 = cgWave.dbase.get<RealArray>("eig");
        const Real & omega = cgWaveHoltz.dbase.get<real>("omega");

        const int numKnownEigs = eig.getLength(1); 
        const int & numEigsToCompute = cgWave.dbase.get<int>("numEigsToCompute"); // number of eigenpairs to compute
        const int numInSum = min(2*numEigsToCompute+1,numKnownEigs);
        printF("Set initial condition to be the sum of %d known eigenvectors. There are %d known eigenvectors\n",numInSum,numKnownEigs);

        printF("*** TO DO : find eigenvectors closest to target*****\n");
        for( int ie=0; ie<numInSum; ie++ )
        {
            const Real lambda=eig(0,ie); 
            printF("Adding eigenvector=%d: eig=%12.4e to the sum...\n",ie,lambda);

            Index I1,I2,I3;
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);

                getIndex(mg.dimension(),I1,I2,I3);
                int includeParallelGhost=1;
                bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeParallelGhost);      
                if( ok )  
                {      
                    Real scale = 1./(lambda*lambda);  // does this matter ? 
                    if( ie==0 )
                        vLocal(I1,I2,I3) = uevLocal(I1,I2,I3,ie)*scale;
                    else
                        vLocal(I1,I2,I3) += uevLocal(I1,I2,I3,ie)*scale;
                }
            }
        }
    // replot=true;      

    // printF("solveEigen:TESTING: call cgWave.updateEigenmodes ...\n");
    // cgWave.updateEigenmodes();
    // OV_ABORT("stop here for now");
    }   

    if(  smoothInitialCondition &&  eigenSolverInitialCondition != CgWave::defaultEigenSolverInitialCondition )
    {
        smoothEigenSolverInitialCondition();
    }
    return 0;
}

// ============================================================================================
/// \brief Smooth the initial condition for eigen solves
// ============================================================================================
int CgWaveHoltz::smoothEigenSolverInitialCondition()
{
    CgWaveHoltz & cgWaveHoltz = *this;
    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
    CgWave::EigenSolverInitialConditionEnum & eigenSolverInitialCondition 
                                                                                      = cgWave.dbase.get<CgWave::EigenSolverInitialConditionEnum>("eigenSolverInitialCondition");
    const int numberOfDimensions = cg.numberOfDimensions();

  // ---- SMOOTH THE INITIAL CONDITIONS ---
    realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
    CompositeGrid & cg = *v.getCompositeGrid();
    Index I1,I2,I3;      

    const int & initialVectorSmooths = cgWave.dbase.get<int>("initialVectorSmooths");
    const int numSmooths= initialVectorSmooths; 

    printF("smoothEigenSolverInitialCondition: SMOOTH (EXISTING) INITIAL CONDITIONS with %d smooths\n",initialVectorSmooths);

    Real omega = 4./5.; // optimal value for smoothing in 2D is omega=4/5 => smoothing rate mu = 3/5
    if( numberOfDimensions==3 )
        omega= 6./7.;   // Owl book p 73 -- optimal smoothing parameter in 3D,  mu=5/7 smoothing fractor

  //   (.6)^10 = .00604...
  //   (.6)^20 = 3.6 e -5 
    for( int itsm=0; itsm<numSmooths; itsm++ )
    {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
            MappedGrid & mg = cg[grid];
            getIndex(mg.gridIndexRange(),I1,I2,I3);
            OV_GET_SERIAL_ARRAY( Real,v[grid],vLocal);
            int includeParallelGhost=0;
            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeParallelGhost);      
            if( ok )  
            {        
                if( numberOfDimensions==2 )
                {
                    vLocal(I1,I2,I3) =  (1.-omega)*(vLocal(I1,I2,I3) ) +
                                                          (.25*omega)*(vLocal(I1+1,I2  ,I3) + 
                                                                                    vLocal(I1-1,I2  ,I3) + 
                                                                                    vLocal(I1  ,I2+1,I3) +
                                                                                    vLocal(I1  ,I2-1,I3) );
                }
                else
                {
                    vLocal(I1,I2,I3) =  (1.-omega)*(vLocal(I1,I2,I3) ) +
                                                            (omega/6.)*(vLocal(I1+1,I2  ,I3) + 
                                                                                    vLocal(I1-1,I2  ,I3) + 
                                                                                    vLocal(I1  ,I2+1,I3) +
                                                                                    vLocal(I1  ,I2-1,I3) +
                                                                                    vLocal(I1  ,I2,I3+1) +
                                                                                    vLocal(I1  ,I2,I3+1)                                             
                                                                                    );            
                }
            }
        }
        cgWave.applyEigenFunctionBoundaryConditions( v );
    }


    return 0;
}



// // ============================================================================================
// /// \brief Plot the WaveHoltz filter function and any known eigenvalues 
// // ============================================================================================
// int CgWaveHoltz::plotFilter( RealArray & eigenValues, GL_GraphicsInterface & ps, PlotStuffParameters & psp )
// {

//   CgWaveHoltz & cgWaveHoltz = *this;
//   CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");

//   // Plot WaveHoltz beta function 


//   const int & computeEigenmodes     = cgWave.dbase.get<int>("computeEigenmodes");
//   RealArray & frequencyArray        = cgWave.dbase.get<RealArray>("frequencyArray");
//   RealArray & frequencyArrayAdjusted= cgWave.dbase.get<RealArray>("frequencyArrayAdjusted");
//   RealArray & periodArray           = cgWave.dbase.get<RealArray>("periodArray"); 
//   RealArray & periodArrayAdjusted   = cgWave.dbase.get<RealArray>("periodArrayAdjusted"); 
//   IntegerArray & numPeriodsArray    = cgWave.dbase.get<IntegerArray>("numPeriodsArray");
//   RealArray & frequencyArraySave    = cgWave.dbase.get<RealArray>("frequencyArraySave");
//   RealArray & periodArraySave       = cgWave.dbase.get<RealArray>("periodArraySave");   


//   // const RealArray & frequencyArray         = cgWave.dbase.get<RealArray>("frequencyArray");
//   // const RealArray & frequencyArrayAdjusted = cgWave.dbase.get<RealArray>("frequencyArrayAdjusted");
//   const Real & dt                          = cgWave.dbase.get<real>("dtUsed");
//   // RealArray & eigenValues            = dbase.get<RealArray>("eigenValues");  // best estimate of eigenvalues
//   // const int & numEigenVectors        = dbase.get<int>("numEigenVectors");



//   const int numEigenVectors = eigenValues.getLength(0); //  THIS IS NUM TO DEFLATE  ********* FIX ME *******


//   const Real dtSaved = cgWave.dbase.get<real>("dt"); 
//   if( !computeEigenmodes )
//   {
//     // -- We need to use the adjusted frequencies in getWaveHoltzIterationEigenvalue ---
//     frequencyArray = frequencyArrayAdjusted;
//     periodArray    = periodArrayAdjusted; 

//     cgWave.dbase.get<real>("dt") = cgWave.dbase.get<real>("dtUsed");
//   }
//   printF("\n ^^^^^^ plotFilter: numEigenVectors=%d, dt=%9.3e, frequencyArray(0)=%10.3e frequencyArrayAdjusted(0)=%10.3e\n",
//     numEigenVectors, dt, frequencyArray(0),frequencyArrayAdjusted(0) );

    
//   Real maxEig = max(eigenValues);

//   // --- plot WaveHoltz filter mu ----
//   int numLambda=201;
//   Real lamMin=0., lamMax= 2*maxEig+20; // 60;
//   RealArray lambda(numLambda);
//   RealArray mu(numLambda);
//   for( int i=0; i<numLambda; i++ )
//   {
//     lambda(i) = lamMin + (lamMax-lamMin)*(i-1.)/(numLambda-1.);
//   }

//   // -- Evaluate the beta function at all points in lambda(:) --
//   bool useAdjusted= true; // !computeEigenmodes; 
    
//   cgWave.getWaveHoltzIterationEigenvalue( lambda, mu, useAdjusted );


//   ps.erase();

//   RealArray points, value; 

//   // --- Mark true eigenvalues if known ----
//   // eigTrue(0:1,0:)
//   const bool trueEigenPairsKnown = cgWave.dbase.has_key("uev");
//   if( trueEigenPairsKnown )
//   {
//     RealArray & eigTrue                = cgWave.dbase.get<RealArray>("eig");
//     const int numberOfEigenvectorsTrue = eigTrue.getLength(1);  
        

//     if( numberOfEigenvectorsTrue>0 )
//     {
//       // for( int ie=0; ie<numberOfEigenvectorsTrue; ie++ )
//       //   printF("plotfilter: ie=%3d eigTrue=%12.5e\n",ie,eigTrue(0,ie));

//       int numToPlot=0;
//       for( int ie=0; ie<numberOfEigenvectorsTrue; ie++ )
//       {
//         if( eigTrue(0,ie)>lamMax )
//           break;
//         numToPlot=ie+1;
//       }
//       assert( numToPlot > 0 );

//       RealArray lam(numToPlot);
//       for( int ie=0; ie<numToPlot; ie++ )
//       {
//         lam(ie) = eigTrue(0,ie);
//       }
//       RealArray betaTrue(numToPlot);
//       cgWave.getWaveHoltzIterationEigenvalue( lam, betaTrue, useAdjusted );  

//       points.redim(numToPlot,2);
//       value.redim(numToPlot); 
//       for( int ie=0; ie<numToPlot; ie++ )
//       {
//         points(ie,0) = eigTrue(0,ie);
//         points(ie,1) = betaTrue(ie);
//         // printF("ie=%3d: true: eig=%12.4e adjusted=%12.4e beta=%12.4e\n",ie,eigTrue(0,ie),lam(ie),betaTrue(ie));
//         value(ie)    = 0.2;                // colour point i by value(i)
//       } 
//       // --- plot locations of true eigenavleus 
//       printF("plot numPlot=%d true values\n",numToPlot);
//       psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);        
//       psp.set(GI_POINT_SIZE,(real) 4.);  // size in pixels
//       psp.set(GI_POINT_COLOUR,"BLACK");

//       ps.plotPoints(points,psp); // colour point i by value(i)  

//       // ps.plotPoints(points,value,psp); // colour point i by value(i)   
//     } 
//   }       

//   // ---- mark eigenvalues found as points ----
//   if( eigenValues.getLength(0)>0 )
//   {
//     RealArray beta(numEigenVectors);
//     points.redim(numEigenVectors,2);
//     value.redim(numEigenVectors);

//     cgWave.getWaveHoltzIterationEigenvalue( eigenValues, beta, useAdjusted );

//     value=.6; 
//     for( int ie=0; ie<numEigenVectors; ie++ )
//     {
//       points(ie,0) = eigenValues(ie); // plot versus original
//       points(ie,1) = beta(ie);

//       // printF("ie=%3d eig=%12.4e beta=%12.4e omega=%12.4e\n",ie,eigenValues(ie),beta(ie),frequencyArrayAdjusted(0));
//       // value(ie)    = 0.4 + .6* (ie+1.)/(numEigenVectors+1.); // point colour in [0,1]

//       value(ie)    =  .2 + .8 * fabs(beta(ie)+.5)/(1.5);  // defines the colout of the point in [0,1]
//     }

//     psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);        
//     psp.set(GI_POINT_SIZE,(real) 8.);  // size in pixels

//     // ps.plotPoints(points,value,psp); // colour point i by value(i)

//     // psp.set(GI_POINT_COLOUR,"LIGHTBLUE");
//     // psp.set(GI_POINT_COLOUR,"BLUE");
//     psp.set(GI_POINT_COLOUR,"MEDIUMAQUAMARINE");
//     ps.plotPoints(points,psp); 

//     aString tName="omega";
//     aString title=sPrintF("WaveHoltz Filter, omega=%6.2f",frequencyArray(0));
//     aString xName[1];
//     xName[0]="beta";

//     psp.set(GI_PLOT_GRID_POINTS_ON_CURVES,false);
//     psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);  
        
//     // PlotIt::plot(ps,t,x,title,tName,xName,psp);
//     #ifndef USE_PPP
//       PlotIt::plot(ps,lambda,mu,title,tName,xName,psp);
//     #else
//       printF("FINISH ME FOR PARALLEL: PlotIt::plot(ps,lambda,mu,title,tName,xName,psp)\n");
//     #endif

//   }

//   psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);        

//   if( !computeEigenmodes )
//   {
//     // -- reset ---
//     frequencyArray = frequencyArraySave;
//     periodArray    = periodArraySave;
//     cgWave.dbase.get<real>("dt") = dtSaved; 

//   }

//   return 0;
// }

// ============================================================================================
/// \brief Solve for eigenpairs using WaveHoltz
// ============================================================================================
int CgWaveHoltz::solveEigen(int argc,char **argv)
{

    CgWaveHoltz & cgWaveHoltz = *this;
    CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");

    Real & convergenceRate          = cgWaveHoltz.dbase.get<Real>("convergenceRate");
    Real & convergenceRatePerPeriod = cgWaveHoltz.dbase.get<Real>("convergenceRatePerPeriod");

    const int & orderOfAccuracy              = cgWave.dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInTime        = cgWave.dbase.get<int>("orderOfAccuracyInTime");
    const int & upwind                       = cgWave.dbase.get<int>("upwind"); // new way
    const int & adjustHelmholtzForUpwinding  = cgWave.dbase.get<int>("adjustHelmholtzForUpwinding"); 
    const int & computeErrors                = cgWave.dbase.get<int>("computeErrors");
    const RealArray & frequencyArrayAdjusted = cgWave.dbase.get<RealArray>("frequencyArrayAdjusted");
    const IntegerArray & numPeriodsArray     = cgWave.dbase.get<IntegerArray>("numPeriodsArray");
    const RealArray & periodArray            = cgWave.dbase.get<RealArray>("periodArray"); 
    int & deflateWaveHoltz                   = cgWave.dbase.get<int>("deflateWaveHoltz"); 
    int & numToDeflate                       = cgWave.dbase.get<int>("numToDeflate");
    aString & eigenVectorFile                = cgWave.dbase.get<aString>("eigenVectorFile");  
    CgWave::EigenSolverEnum & eigenSolver    = cgWave.dbase.get<CgWave::EigenSolverEnum>("eigenSolver"); 
    int & numberOfRitzVectors                = cgWave.dbase.get<int>("numberOfRitzVectors");
    int & initialVectorsForEigenSolver       = cgWave.dbase.get<int>("initialVectorsForEigenSolver");
    int & initialVectorSmooths               = cgWave.dbase.get<int>("initialVectorSmooths");
    const int & useAccurateInnerProduct      = cgWave.dbase.get<int>("useAccurateInnerProduct");

    CgWave::EigenSolverInitialConditionEnum & eigenSolverInitialCondition 
                                                                                      = cgWave.dbase.get<CgWave::EigenSolverInitialConditionEnum>("eigenSolverInitialCondition");

    cgWaveHoltz.dbase.get<int>("orderOfAccuracy")=orderOfAccuracy; // set value in CgWaveHoltz

    const int & numberOfFrequencies          = cgWaveHoltz.dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray         = cgWaveHoltz.dbase.get<RealArray>("frequencyArray");

    real & omega                             = cgWaveHoltz.dbase.get<real>("omega");
    int & numPeriods                         = cgWaveHoltz.dbase.get<int>("numPeriods");
    real & tol                               = cgWaveHoltz.dbase.get<real>("tol");
    int & adjustOmega                        = cgWaveHoltz.dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
    int & maximumNumberOfIterations          = cgWaveHoltz.dbase.get<int>("maximumNumberOfIterations");
    int & numberOfIterations                 = cgWaveHoltz.dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken
    int & cgWaveDebugMode                    = cgWaveHoltz.dbase.get<int>("cgWaveDebugMode");

    aString nameOfShowFile = "eigenWave.show";
    bool showFileIsOpen=false;
    bool smoothInitialCondition=false;

    bool plotOption=true;  // by default we plot interactively
    GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("EigenWave",plotOption,argc,argv));

    PlotStuffParameters psp;  

    const int numberOfDimensions = cg.numberOfDimensions();

  // Build a dialog menu for changing parameters
    GUIState gui;
    DialogData & dialog=gui;

    dialog.setWindowTitle("EigenWave - EigenMode Solver");
    dialog.setExitCommand("exit", "exit");

    aString pbLabels[] = {
                                                "compute",
                        // "zero initial condition",
                        // "random initial condition",
                        // "sine initial condition",
                        // "eigenvector initial condition",
                        // "smooth initial condition",
                                                "change parameters",
                                                "plot initial condition",
                                                "grid",
                                                "plot residuals",
                                                "plot eigenVectors",
                                                "run cgWave and plot",
                                                "save show file",
                                                "plot errors",
                                                "plot filter",
                        // "plot sequences",
                                                "erase",
                                                "exit",
                                                ""};
    int numRows=9;
    dialog.setPushButtons( pbLabels, pbLabels, numRows ); 

    dialog.setOptionMenuColumns(1);
    aString eigenSolverLabel[] = { "default", "KrylovSchur", "Arnoldi", "Arpack", "fixedPoint", "power", "inverseIteration", "JacobiDavidson", "" };
    dialog.addOptionMenu("EigenSolver:", eigenSolverLabel, eigenSolverLabel, (int)eigenSolver );

    aString eigenSolverInitialConditionLabel[] = { "defaultIC", "randomIC", "sineIC", "eigenvectorIC", "" };
    dialog.addOptionMenu("Initial Condition:", eigenSolverInitialConditionLabel, eigenSolverInitialConditionLabel, (int)eigenSolverInitialCondition );


  // ----- Text strings ------
    const int numberOfTextStrings=20;
    aString textCommands[numberOfTextStrings];
    aString textLabels[numberOfTextStrings];
    aString textStrings[numberOfTextStrings];

    int nt=0;
    textCommands[nt] = "omega";  textLabels[nt]=textCommands[nt];
    sPrintF(textStrings[nt], "%g",omega);  nt++; 

    textCommands[nt] = "number of periods";  textLabels[nt]=textCommands[nt];
    sPrintF(textStrings[nt], "%i",numPeriods);  nt++; 

    textCommands[nt] = "tol";  textLabels[nt]=textCommands[nt];
    sPrintF(textStrings[nt], "%g",tol);  nt++; 

    textCommands[nt] = "max iterations";  textLabels[nt]=textCommands[nt];
    sPrintF(textStrings[nt], "%i",maximumNumberOfIterations);  nt++; 

    textCommands[nt] = "show file";  textLabels[nt]=textCommands[nt];
    sPrintF(textStrings[nt], "%s",(const char*)nameOfShowFile);  nt++; 

    textCommands[nt] = "number to deflate";  textLabels[nt]=textCommands[nt];
    sPrintF(textStrings[nt], "%i",numToDeflate);  nt++; 

    textCommands[nt] = "max Ritz vectors";  textLabels[nt]=textCommands[nt];
    sPrintF(textStrings[nt], "%i",numberOfRitzVectors);  nt++;  

    textCommands[nt] = "init vector smooths";  textLabels[nt]=textCommands[nt];
    sPrintF(textStrings[nt], "%i",initialVectorSmooths);  nt++;  

  // null strings terminal list
    textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
    dialog.setTextBoxes(textCommands, textLabels, textStrings);

    
    aString tbCommands[] = {
                                                    "run cgWave with debugging",
                                                    "deflate",
                                                    "set initial vectors",
                                                    "smooth initial condition",
                                                        ""};
    int tbState[10];
    
    tbState[0] = cgWaveDebugMode;
    tbState[1] = deflateWaveHoltz;
    tbState[2] = initialVectorsForEigenSolver;
    tbState[3] = smoothInitialCondition;

    int numColumns=1;
    dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 


    ps.pushGUI(gui);



    psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

    aString answer;
    char buff[200];

  // keep track of what has been plotted
    int plotChoices=0; // 1=grid, 2=contour
    bool replot=false;
    bool reComputeErrors=false; 
  // bool saveCheckFile=false; // set to true when errors should be saved to the check file.

  // int checkFileCounter=0; // keeps track of how many times the checkFile is saved with results

  // uHelmholtz holds Helmholtz solution from direct solver
  //   This is optionally used in CgWave as a known solution
    if( !cgWave.dbase.has_key("uHelmholtz") )
        cgWave.dbase.put<realCompositeGridFunction>("uHelmholtz");

    Range all;
  // realCompositeGridFunction & uHelmholtz = cgWave.dbase.get<realCompositeGridFunction>("uHelmholtz"); 
  // uHelmholtz.updateToMatchGrid(cg,all,all,all,numberOfFrequencies); // save Helmholtz solution here

    bool helmholtzFromDirectSolverWasComputed=false;
    Real cpuSolveHelmholtz=-1.;
    bool initialConditionsHaveBeenAssigned=false;

    Ogshow show;  // show file for saving solutions

    int current=0;
    for(;;) 
    {
        ps.getAnswer(answer,"");      
        if( answer=="exit" || answer=="continue" )
        {
            break;
        }
        else if( answer=="default"          ||
                          answer=="KrylovSchur"      ||
                          answer=="Arnoldi"          ||
                          answer=="Arpack"           ||
                          answer=="fixedPoint"       ||
                          answer=="power"            ||
                          answer=="inverseIteration" ||
                          answer=="JacobiDavidson"  
                          )
        {
            if( answer=="default" )
                eigenSolver=CgWave::defaultEigenSolver;
            else if( answer=="KrylovSchur" )
                eigenSolver=CgWave::KrylovSchurEigenSolver;
            else if( answer=="Arnoldi" )
                eigenSolver=CgWave::ArnoldiEigenSolver;
            else if( answer=="Arpack" )
                eigenSolver=CgWave::ArpackEigenSolver;      
            else if( answer=="fixedPoint" )
                eigenSolver=CgWave::fixedPointEigenSolver;  
            else if( answer=="power" )
                eigenSolver=CgWave::powerEigenSolver;  
            else if( answer=="inverseIteration" )
                eigenSolver=CgWave::inverseIterationEigenSolver;  
            else if( answer=="JacobiDavidson" )
                eigenSolver=CgWave::JacobiDavidsonEigenSolver;  
            else
            {
                printF("Unknown eigenSolver =[%s]\n",(const char*)answer);
                OV_ABORT("this should not happen");
            }            
            printF("Setting eigenSolver=[%s]\n",(const char*)answer);
        }

        else if( answer=="change parameters" )
        {
            cgWaveHoltz.interactiveUpdate();
        }
        
        else if( dialog.getTextValue(answer,"omega","%e",omega) )
        {
            printF("Setting omega=%g\n",omega);

      // --- DO THIS FOR NOW ** FIX ME ****************************************************
            cgWave.dbase.get<Real>("omega")  = omega;   

            RealArray & cgWaveFrequencyArray = cgWave.dbase.get<RealArray>("frequencyArray");
            RealArray & cgWavePeriodArray    = cgWave.dbase.get<RealArray>("periodArray");
            Real & Tperiod                   = cgWave.dbase.get<real>("Tperiod");
            Real & tFinal                    = cgWave.dbase.get<real>("tFinal");
            const int & numPeriods           = cgWave.dbase.get<int>("numPeriods");

            RealArray & frequencyArray = dbase.get<RealArray>("frequencyArray");
            RealArray & periodArray    = dbase.get<RealArray>("periodArray");

            cgWaveFrequencyArray(0) = omega;
            dbase.get<Real>("omega")  = omega;   
            frequencyArray=omega;

            Tperiod=numPeriods*twoPi/omega;  
        
            tFinal                     = Tperiod; 
            cgWavePeriodArray(0)       = Tperiod; 
            periodArray                = Tperiod;
            dbase.get<real>("Tperiod") = Tperiod;

            printF("*** SETTING tFinal = %9.3e\n",tFinal);   

            cgWave.initialize();
        }
        else if( dialog.getTextValue(answer,"number of periods","%i",numPeriods) )
        {
            printF("Setting numPeriods=%i\n",numPeriods);
        }
        else if( dialog.getTextValue(answer,"tol","%e",tol) )
        {
            printF("Setting tol=%g (tolerence for Krylov solvers)\n",tol);
        }
        else if( dialog.getTextValue(answer,"max iterations","%i",maximumNumberOfIterations) )
        {
            printF("Setting maximumNumberOfIterations=%i. **NOTE*** These are OUTER iterations for Arnoldi algorithms.\n",maximumNumberOfIterations);

        }
        else if( dialog.getTextValue(answer,"number to deflate","%i",numToDeflate) )
        {
            printF("Setting numToDeflate=%i\n",numToDeflate);
            cgWave.reinitializeDeflation();
        }  
        else if( dialog.getTextValue(answer,"max Ritz vectors","%i",numberOfRitzVectors) )
        {
            printF("Setting numberOfRitzVectors=%i\n",numberOfRitzVectors);
        }  

        else if( dialog.getTextValue(answer,"init vector smooths","%i",initialVectorSmooths) )
        {
            printF("Setting initialVectorSmooths=%i\n",initialVectorSmooths);
        }  

        else if( dialog.getTextValue(answer,"show file","%s",nameOfShowFile) )
        {
            printF("Setting nameOfShowFile=%s\n",(const char*)nameOfShowFile);
            if( showFileIsOpen )
            {
        // close any open show file
                show.close();
                showFileIsOpen=false;
            }
        }   
        else if( dialog.getToggleValue(answer,"deflate",deflateWaveHoltz) )
        {
            printF("Setting deflateWaveHoltz=%d\n",deflateWaveHoltz);
        }

        else if( dialog.getToggleValue(answer,"run cgWave with debugging",cgWaveDebugMode) )
        {
            printF("Setting cgWaveDebugMode=%d (1: run cgWave plotting every step)\n",cgWaveDebugMode);
            if( cgWaveDebugMode==1 )
            {

            }
        } 
        else if( dialog.getToggleValue(answer,"set initial vectors",initialVectorsForEigenSolver) )
        {
            printF("Setting initialVectorsForEigenSolver=%d\n",initialVectorsForEigenSolver);
        }  
        else if( dialog.getToggleValue(answer,"smooth initial condition",smoothInitialCondition) )
        {
            printF("Setting smoothInitialCondition=%d. Smoothing will happen when the initial conditions are next set.\n",
                (int)smoothInitialCondition);
        }                 

        else if( answer=="compute" )
        {
            const bool useFixedPoint =  eigenSolver==CgWave::fixedPointEigenSolver;

            const int & computeEigenmodes     = cgWave.dbase.get<int>("computeEigenmodes");
            const int readTrueEigenPairs = eigenVectorFile != "none";
            if( computeEigenmodes && readTrueEigenPairs )
            {
        // -- This next call will read in any known eigenmodes ---
                cgWave.initializeDeflation();
            }

            printF("\n\n )))))))))))  initialConditionsHaveBeenAssigned=%d ((((((((((\n\n",(int)initialConditionsHaveBeenAssigned);
            if( !initialConditionsHaveBeenAssigned )
            {
        // --- initial conditions may have been assigned yet, e.g. if we re-compute ---
                assignEigenSolverInitialCondition( smoothInitialCondition );
            }
            initialConditionsHaveBeenAssigned=false;  // this means we need to re-assign the IC's next time,  since the grid function v is changed

            const Real cpu0=getCPU();
            Real cpuSolve; 
            

      // -----------------------------------------------------------------
      // ------------------WAVE HOLTZ SOLVE-------------------------------
      // -----------------------------------------------------------------
            if( !useFixedPoint )
            {
                cgWaveHoltz.solveSLEPc(argc,argv);

                cpuSolve = getCPU()-cpu0;
      
            }
            else
            {
        // -- My fixedPoint method:
                cgWaveHoltz.solve(); 

                cpuSolve = getCPU()-cpu0;

        // **********  TEMP *************** FIX ME
                if( !dbase.has_key("numEigenVectors" ))
                    dbase.put<int>("numEigenVectors")=0;

                int & numEigenVectors = dbase.get<int>("numEigenVectors");

                numEigenVectors = cgWave.dbase.get<int>("numEigenVectors"); // from CgWave::updateEigenmodes()

                if( !dbase.has_key("eigenValues") )
                    dbase.put<RealArray>("eigenValues");

                dbase.get<RealArray>("eigenValues") = cgWave.dbase.get<RealArray>("eigenValues");

                if( !dbase.has_key("eigenVectorRayleighRitz") )
                {
                    dbase.put<realCompositeGridFunction>("eigenVectorRayleighRitz");
                }

                if( !dbase.has_key("eigenVectorRayleighRitz") )
                    dbase.put<realCompositeGridFunction>("eigenVectorRayleighRitz");

                realCompositeGridFunction & eigenVector = dbase.get<realCompositeGridFunction>("eigenVectorRayleighRitz");
                eigenVector.updateToMatchGrid(cg,all,all,all,numEigenVectors);
                  for( int i=0; i<numEigenVectors; i++ )
                    eigenVector.setName(sPrintF("phi%d",i),i);

                realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");  // here is where the eigvect is stored
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    eigenVector[grid] = v[grid];
                }

            }


      // if( deflateWaveHoltz )
      //   cgWave.inflateSolution(); // un-deflate the solution

  

      // The period array may have changed since numPeriods(freq) may have changed
            cgWaveHoltz.dbase.get<RealArray>("periodArray") = cgWave.dbase.get<RealArray>("periodArray");

            cgWaveHoltz.dbase.get<Real>("cpuSolveEigenWave")=cpuSolve;
            

            const int & numEigenVectors = dbase.get<int>("numEigenVectors"); // from CgWave::updateEigenmodes()

            if( numEigenVectors>0 )
            {
        // ------------------------------------
        // ----- output a table of results ----
        // ------------------------------------

                outputEigenTable( );



                replot=true;
                reComputeErrors=true;

        // save results to a matlab file
                const aString & matlabFileName = cgWaveHoltz.dbase.get<aString>("matlabFileName");  
                aString localName; 
                if( useFixedPoint )
                {
                    cgWaveHoltz.dbase.get<aString>("solverName")="FixedPoint";
                    localName = matlabFileName + "FP"; 
                }
                else
                {
                    cgWaveHoltz.dbase.get<aString>("solverName")="Krylov";
                    localName = matlabFileName + "Krylov"; 
                }
                cgWaveHoltz.outputMatlabFile( localName );

            }

        }


        else if( answer=="compute errors" )
        {
            if( computeErrors )
                reComputeErrors=true; // errors are computed below
            else
                printF("compute errors: WARNING: this solution has no errors to compute. Probably not a known solution.\n"); 
        }
        else if( answer == "zero initial condition" )
        {
      // Keep for backward compat
      //   realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
      //   v=0.;

      //   cgWave.resetTimings(); // reset CPU timings to zero 
        }
        else if( answer=="defaultIC"      ||
                          answer=="randomIC"       || answer == "random initial condition"      ||
                          answer=="sineIC"         || answer == "sine initial condition"        ||
                          answer== "eigenvectorIC" || answer == "eigenvector initial condition"
                          )
        {
            if( answer=="defaultIC" )
                eigenSolverInitialCondition = CgWave::defaultEigenSolverInitialCondition;
            else if( answer=="randomIC"  ||  answer == "random initial condition" )
                eigenSolverInitialCondition = CgWave::randomEigenSolverInitialCondition;
            else if( answer=="sineIC"  || answer == "sine initial condition" )
                eigenSolverInitialCondition = CgWave::sineEigenSolverInitialCondition;      
            else if( answer== "eigenvectorIC" || answer == "eigenvector initial condition" )
                eigenSolverInitialCondition = CgWave::sumOfEigenvectorsInitialCondition;
            else
            {
                OV_ABORT("ERROR - this should not happen");
            }

            assignEigenSolverInitialCondition(smoothInitialCondition);

            initialConditionsHaveBeenAssigned=true;
            replot=true;

            cgWave.resetTimings(); // reset CPU timings to zero       

        }    
    // else if( answer=="randomIC"  ||  answer == "random initial condition" )
    // {
    //   printF("solveEigen: Set RANDOM initial conditions for SLEPc\n");
    //   eigenSolverInitialCondition = CgWave::randomEigenSolverInitialCondition;

    //   realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
    //   CompositeGrid & cg = *v.getCompositeGrid();
    //   Index I1,I2,I3;

    //   std::srand(12789.);

    //   const Real scale=1.; // 1.e-5; // TEST
    //   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //   {
    //     MappedGrid & mg = cg[grid];
    //     OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
    //     OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

    //     vLocal=0.;
    //     getIndex(mg.gridIndexRange(),I1,I2,I3);
    //     for( int freq=0; freq<numberOfFrequencies; freq++ )
    //     {
    //       FOR_3D(i1,i2,i3,I1,I2,I3)
    //       {
    //         // do this for now
    //         // vLocal(i1,i2,i3) = sin(i1)*cos(i2);
    //         if( maskLocal(i1,i2,i3)!=0 )
    //         {
    //           vLocal(i1,i2,i3,freq) = (-1. + std::rand()*(2./RAND_MAX))*scale; // [-1,1]
    //         }

    //       }
    //     }
    //   }
    //   replot=true;

    //   cgWave.resetTimings(); // reset CPU timings to zero 
    // }  

    // else if( answer=="sineIC"  || answer == "sine initial condition" )
    // {
    //   eigenSolverInitialCondition = CgWave::sineEigenSolverInitialCondition;
            
    //   realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
    //   CompositeGrid & cg = *v.getCompositeGrid();
    //   Index I1,I2,I3;

    //   // ---- Base frequency of the sine solution on the target frequency ----
    //   const Real & c  = cgWave.dbase.get<Real>("c");
    //   const Real k = omega/c; 
    //   const Real theta = .5*(twoPi/4.); // Pi/4 

    //   const Real kxHat = cos(theta);
    //   const Real kyHat = sin(theta);

    //   const Real kx = kxHat*k;
    //   const Real ky = kyHat*k;
    //   const Real kz = kx;      
            
    //   const Real scale=1.; // should not matter
    //   printF("solveEigen: Set initial condition sin(kx*x)*sin(ky*y) [*sin(kz*z)] kx=%9.2e ky=%9.2e [kz=%8.2e]\n",kx,ky,kz);
    //   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //   {
    //     MappedGrid & mg = cg[grid];
    //     OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
    //     OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);
    //     mg .update(MappedGrid::THEvertex);
    //     OV_GET_SERIAL_ARRAY( Real, mg.vertex(),xLocal);        

    //     vLocal=0.;
    //     getIndex(mg.gridIndexRange(),I1,I2,I3);
    //     FOR_3D(i1,i2,i3,I1,I2,I3)
    //     {
    //       if( maskLocal(i1,i2,i3)!=0 )
    //       {
    //         if( numberOfDimensions==2 )
    //           vLocal(I1,I2,I3) = scale*sin( kx*xLocal(I1,I2,I3,0) ) * sin( ky*xLocal(I1,I2,I3,1) );
    //         else
    //           vLocal(I1,I2,I3) = scale*sin( kx*xLocal(I1,I2,I3,0) ) * sin( ky*xLocal(I1,I2,I3,1) )* sin( kz*xLocal(I1,I2,I3,2) );

    //       }

    //     }
    //   }


    // }

    // else if( answer=="smooth initial condition" )
    // {
    //   // ---- SMOOTH THE INTIAL CONDITIONS ---
    //   realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
    //   CompositeGrid & cg = *v.getCompositeGrid();
    //   Index I1,I2,I3;      

    //   const int & initialVectorSmooths = cgWave.dbase.get<int>("initialVectorSmooths");
    //   const int numSmooths= initialVectorSmooths; 

    //   printF("solveEigen: SMOOTH (EXISTING) INITIAL CONDITIONS with %d smooths\n",initialVectorSmooths);

    //   Real omega = 4./5.; // optimal value for smoothing in 2D is omega=4/5 => smoothing rate mu = 3/5
    //   if( numberOfDimensions==3 )
    //     omega= 6./7.;   // Owl book p 73 -- optimal smoothing parameter in 3D,  mu=5/7 smoothing fractor

    //   //   (.6)^10 = .00604...
    //   //   (.6)^20 = 3.6 e -5 
    //   for( int itsm=0; itsm<numSmooths; itsm++ )
    //   {
    //     for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //     {
    //       MappedGrid & mg = cg[grid];
    //       getIndex(mg.gridIndexRange(),I1,I2,I3);
    //       OV_GET_SERIAL_ARRAY( Real,v[grid],vLocal);

    //       if( numberOfDimensions==2 )
    //       {
    //         vLocal(I1,I2,I3) = (1.-omega)*(vLocal(I1,I2,I3) ) +
    //                            (.25*omega)*(vLocal(I1+1,I2  ,I3) + 
    //                                         vLocal(I1-1,I2  ,I3) + 
    //                                         vLocal(I1  ,I2+1,I3) +
    //                                         vLocal(I1  ,I2-1,I3) );
    //       }
    //       else
    //       {
    //         vLocal(I1,I2,I3) = (1.-omega)*(vLocal(I1,I2,I3) ) +
    //                             (omega/6.)*(vLocal(I1+1,I2  ,I3) + 
    //                                         vLocal(I1-1,I2  ,I3) + 
    //                                         vLocal(I1  ,I2+1,I3) +
    //                                         vLocal(I1  ,I2-1,I3) +
    //                                         vLocal(I1  ,I2,I3+1) +
    //                                         vLocal(I1  ,I2,I3+1)                                             
    //                                         );            
    //       }
    //     }
    //     cgWave.applyEigenFunctionBoundaryConditions( v );
    //   }


    //   replot=true;

    //   cgWave.resetTimings(); // reset CPU timings to zero 
    // }  

    // else if( answer== "eigenvectorIC" || answer == "eigenvector initial condition" )
    // {
    //   eigenSolverInitialCondition = CgWave::sumOfEigenvectorsInitialCondition;

    //   realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
    //   CompositeGrid & cg = *v.getCompositeGrid();

    //   // --- get eigenvectors
    //   cgWave.initializeDeflation();

    //   // Here are the eigenvectors and eigenvalues:
    //   realCompositeGridFunction & uev = cgWave.dbase.get<realCompositeGridFunction>("uev");
    //   RealArray & eig                 = cgWave.dbase.get<RealArray>("eig");

    //   const int numKnownEigs = eig.getLength(1); 
    //   const int & numEigsToCompute = cgWave.dbase.get<int>("numEigsToCompute"); // number of eigenpairs to compute
    //   const int numInSum = min(2*numEigsToCompute+1,numKnownEigs);
    //   printF("Set initial condition to be the sum of %d known eigenvectors. There are %d known eigenvectors\n",numInSum,numKnownEigs);

    //   printF("*** TO DO : find eigenvectors closest to target*****\n");
    //   for( int ie=0; ie<numInSum; ie++ )
    //   {
    //     const Real lambda=eig(0,ie); 
    //     printF("Adding eigenvector=%d: eig=%12.4e to the sum...\n",ie,lambda);

    //     Index I1,I2,I3;
    //     for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    //     {
    //       MappedGrid & mg = cg[grid];
    //       OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
    //       OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);

    //       getIndex(mg.dimension(),I1,I2,I3);
    //       Real scale = 1./(lambda*lambda);  // does this matter ? 
    //       if( ie==0 )
    //         vLocal(I1,I2,I3) = uevLocal(I1,I2,I3,ie)*scale;
    //       else
    //         vLocal(I1,I2,I3) += uevLocal(I1,I2,I3,ie)*scale;
    //     }
    //   }
    //   replot=true;      

    //   // printF("solveEigen:TESTING: call cgWave.updateEigenmodes ...\n");
    //   // cgWave.updateEigenmodes();
    //   // OV_ABORT("stop here for now");
    // }    

        else if( answer=="save show file" || answer=="save to show" )
        {
            printF("Save the eigenvectors to the show file [%s].\n",(const char*)nameOfShowFile);
            if( !showFileIsOpen )
            {
                showFileIsOpen=true;
                show.open(nameOfShowFile);

                const int & flushFrequency  = cgWave.dbase.get<int>("flushFrequency");  // number of solutions per show file 
                show.setFlushFrequency( flushFrequency ); // save this many solutions per sub-showFile

            }
            
            show.saveGeneralComment("Solutions from EigenWave"); // save a general comment in the show file

            realCompositeGridFunction & eigenVector = dbase.get<realCompositeGridFunction>("eigenVectorRayleighRitz");
            RealArray & eigenValues = dbase.get<RealArray>("eigenValues");  
            const int & numEigenVectors = dbase.get<int>("numEigenVectors");

            const bool trueEigenPairsKnown = cgWave.dbase.has_key("uev");
            const IntegerArray & eigIndex = trueEigenPairsKnown? dbase.get<IntegerArray>("eigIndex") : numPeriodsArray;


            show.startFrame(); 

      // -- save eigenvalues --
            HDF_DataBase *dbp=NULL;

            dbp = show.getFrame();
            assert( dbp!=NULL );

      // save parameters that go in this frame
            HDF_DataBase & db = *dbp;
      // We save the eigs of -Delta 
            RealArray eig(2,numEigenVectors); // holds real and imag parts

            for( int ie=0; ie<numEigenVectors; ie++ ) 
            {
                eig(0,ie) = SQR(eigenValues(ie));   // NOTE SQUARE 
                eig(1,ie) = 0.; // imaginary part
            }

            db.put(eig,"eig");   

      // -- save eigenvectors as separate solutions --

            realCompositeGridFunction q(cg,all,all,all);
            q.setName("phi",0);
            Index I1,I2,I3;
            for( int ie=0; ie<numEigenVectors; ie++ )
            {
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    OV_GET_SERIAL_ARRAY(real,q[grid],qLocal);
                    OV_GET_SERIAL_ARRAY(Real,eigenVector[grid],eigenVectorLocal);
                    getIndex(cg[grid].dimension(),I1,I2,I3);
                    int includeParallelGhost=1;
                    bool ok=ParallelUtility::getLocalArrayBounds(q[grid],qLocal,I1,I2,I3,includeParallelGhost);      
                    if( ok )  
                    {                
                        qLocal(I1,I2,I3)=eigenVectorLocal(I1,I2,I3,ie);
                    }
                }

                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    q[grid].updateGhostBoundaries();
                
                if( ie>0 ) 
                    show.startFrame();                       // start a new frame
                if( trueEigenPairsKnown )
                    show.saveComment(0,sPrintF("EigenWave: FD%i, eig=%d, lambda=%.5g",orderOfAccuracy,eigIndex(ie),eigenValues(ie)));   
                else
                    show.saveComment(0,sPrintF("EigenWave: FD%i, lambda[%d]=%.5g",orderOfAccuracy,ie,eigenValues(ie)));   
                show.saveSolution( q );              // save to show file
                show.endFrame();

            }

            show.close();
            showFileIsOpen=false;
            printF("Wrote show file=[%s]\n",(const char*)nameOfShowFile);
        }

        else if( answer=="contour" )
        {
      // -- this is done below so that we can also treat the case of replot=true ---
        }
        else if( answer=="grid" )
        {
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
            PlotIt::plot(ps,cg,psp);                          // plot the grid
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

        }
        else if( answer=="plot eigenVectors" )
        {
            realCompositeGridFunction & eigenVector = dbase.get<realCompositeGridFunction>("eigenVectorRayleighRitz");
            RealArray & eig                         = dbase.get<RealArray>("eigenValues");
      // RealArray & eigTrue                = cgWave.dbase.get<RealArray>("eig");

            if( cg.numberOfDimensions()==3 )
                cg.update(MappedGrid::THEmask | MappedGrid::THEcenter | MappedGrid::THEvertex ); // This is needed in 3D for some reason 

            plotEigenVectors( eigenVector, eig, "EigenVector", ps,psp );

        }
        else if( answer=="plot residuals" )
        {
            realCompositeGridFunction & res = cgWaveHoltz.dbase.get<realCompositeGridFunction>("residual");

      // RealArray & eig                       = cgWave.dbase.get<RealArray>("eig");
            RealArray & eig                       = dbase.get<RealArray>("eigenValues");

            if( cg.numberOfDimensions()==3 )
                cg.update(MappedGrid::THEmask | MappedGrid::THEcenter | MappedGrid::THEvertex ); // This is needed in 3D for some reason 

            if( true )
            {
                plotEigenVectors( res, eig, "Residual", ps,psp );

            }
            else
            {
        // interpolate residual to make plots nicer: 
                if( false )
                { 
                    printF("INFO: interpolate the residual\n");
                    res.interpolate();
                }
                ps.erase();
                psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
                psp.set(GI_TOP_LABEL,sPrintF("residual O%d omega=%.5g",orderOfAccuracy,omega));
                PlotIt::contour(ps,res,psp);
                psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);
            }

        }
        else if( answer=="plot initial condition" )
        {
              realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");      
              ps.erase();
              psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
              psp.set(GI_TOP_LABEL,sPrintF("Initial conditions for omega=%.5g",omega));
              PlotIt::contour(ps,v,psp);
              psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);      
        }
        else if( answer=="plot errors" )
        {
            const bool trueEigenPairsKnown = cgWave.dbase.has_key("uev");
            if( !trueEigenPairsKnown )
            {
                printF("INFO: Cannot plot errors since true eigenPairs are not known\n");
                continue;
            }

            const int & computeEigenmodes     = cgWave.dbase.get<int>("computeEigenmodes");
            realCompositeGridFunction & error = cgWave.dbase.get<realCompositeGridFunction>("error");
            bool plotErrors=true;
            psp.set(GI_TOP_LABEL,sPrintF("Eigenvector Errors O%d",orderOfAccuracy));

            RealArray & eig = cgWave.dbase.get<RealArray>("eig");
            plotEigenVectors( error, eig, "Error", ps,psp );

        }

        else if( answer=="plot filter" || answer=="plot beta" )
        {
      // Plot WaveHoltz beta function
            RealArray & eigenValues = dbase.get<RealArray>("eigenValues");  // best estimate of eigenvalues
            cgWave.plotFilter( eigenValues );

        }

        else if( answer=="plot sequences" )
        {
      // Plot residual history
            const int & computeEigenmodes     = cgWave.dbase.get<int>("computeEigenmodes");
            
            const RealArray & resVector       = cgWave.dbase.get<RealArray>("resVector");
            if( numberOfIterations<=0 )
            {
                printF("There are no sequences to plot yet\n");
                continue;
            }

            const int nd=numberOfIterations;

            int numberOfComponents=1;
            if( computeEigenmodes && cgWave.dbase.has_key("errEigenvector") )
                numberOfComponents+=2;

            if( computeEigenmodes && cgWave.dbase.has_key("errEigenvectorRR") )
                  numberOfComponents+=2;
            aString xName[10];
            RealArray t(nd),x(nd,numberOfComponents);
            int m=0; // counts components
            xName[0]="log10(res)";
            for( int i=0; i<nd; i++ )
            {
                t(i)=i;
                x(i,m) = log10(max(resVector(i),1e-20));
            }
            if( computeEigenmodes && cgWave.dbase.has_key("errEigenvector") )
            {
        // --- output errors in eigenvalues and eigenvectors if they have been computed ---

                RealArray & errEigenvalue  = cgWave.dbase.get<RealArray>("errEigenvalue");
                RealArray & errEigenvector = cgWave.dbase.get<RealArray>("errEigenvector");
                xName[m+1]="log10(lamRQ-err)";
                xName[m+2]="log10(phi-err)";
                for( int i=0; i<nd; i++ )
                {
                    x(i,m+1) = log10(max(errEigenvalue(i),1e-20));
                    x(i,m+2) = log10(max(errEigenvector(i),1e-20));
                }
                m+=2;
            }
            if( computeEigenmodes &&  cgWave.dbase.has_key("errEigenvectorRR") )
            {
        // --- output errors in eigenvalues and eigenvectors if they have been computed ---
                const int & iterationStartRR = cgWave.dbase.get<int>("iterationStartRR");      // start at this WH iteration

                RealArray & errEigenvalueRR  = cgWave.dbase.get<RealArray>("errEigenvalueRR");
                RealArray & errEigenvectorRR = cgWave.dbase.get<RealArray>("errEigenvectorRR");
                const int & iterationRR      = cgWave.dbase.get<int>("iterationRR");
                const int numIts = min(iterationRR,errEigenvalueRR.getLength(0));
                const int iLast = numIts-1;
                printF("plot sequences: iterationStartRR=%d, iterationRR=%d, iLast=%d\n",iterationStartRR,iterationRR,iLast);
                xName[m+1]="log10(lamRR-err)";
                xName[m+2]="log10(phiRR-err)";
                for( int i=0; i<nd; i++ )
                {
                    int irr = i-iterationStartRR;
                    if( i<iterationStartRR )
                    {
            // Rayleigh-Ritz iterations have not started yet
                        x(i,m+1) = log10(1.);
                        x(i,m+2) = log10(1.);
                    }
          // else
          // {
          //   // These are the only valid values
          //   x(i,m+1) = log10(max(errEigenvalueRR(i),1e-20));  
          //   x(i,m+2) = log10(max(errEigenvectorRR(i),1e-20));
          // }

                    else if( i<iterationStartRR + numIts )
                    {
            // These are the only valid values
                        x(i,m+1) = log10(max(errEigenvalueRR(irr),1e-20));  
                        x(i,m+2) = log10(max(errEigenvectorRR(irr),1e-20));
                    }
                    else
                    {
            // Rayleigh-Ritz iterations ahev stoped.
                        x(i,m+1) = log10(max(errEigenvalueRR(iLast),1e-20));  
                        x(i,m+2) = log10(max(errEigenvectorRR(iLast),1e-20));
                    }
                } 
                m += 2;       
            }

            aString title="WaveHoltz";
            if( computeEigenmodes )
            {
                sPrintF(title,"EigenWave: lambda=%.4g",frequencyArray(0));
            }
            aString tName="iteration";

      // -- plot all components by default ---
            IntegerArray componentsToPlot(numberOfComponents);
            for( int i=0; i<numberOfComponents; i++ )
                componentsToPlot(i)=i;

            psp.set(GI_COMPONENTS_TO_PLOT,componentsToPlot);
  
            ps.erase();
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
            #ifndef USE_PPP
                PlotIt::plot(ps,t,x,title,tName,xName,psp);
            #else
                printF("FINISH ME FOR PARALLEL: PlotIt::plot(ps,t,x,...)\n");        
            #endif
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

            replot=true;

        }    
        else if( answer.matches("erase") )
        {
            replot=false;
            plotChoices=0; 
            ps.erase();
        }
        else if( answer=="run cgWave and plot" )    
        {
            printF("Run cgWave and plot the solution over time...\n");
            int & plotChoices  = cgWave.dbase.get<int>("plotChoices");
            int & plotOptions  = cgWave.dbase.get<int>("plotOptions");

            plotChoices=2; // contours
            plotOptions = CgWave::plotAndWait; // plotNoWait

            cgWave.dbase.get<real>("omega")     = omega;

            int it=0;
            cgWave.advance( it );

        }
        else
        {
            printF("solveEigen: Unknown command = [%s]\n",(const char*)answer);
            ps.stopReadingCommandFile();
              
        }

    // if( reComputeErrors )
    // {
    //   reComputeErrors=false;

    //   if( computeErrors )
    //   {
    //     realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
    //     Real t=0; // compute errors with t=0 in cos(omega*t)
    //     Real maxErr=1.;
    //     const int & computeEigenmodes = cgWave.dbase.get<int>("computeEigenmodes");
    //     // **FIX ME FOR computeEigenmodes
    //     if( !computeEigenmodes ) 
    //     {
    //       Real maxErr = cgWave.getErrors( v, t );
    //       printF("solveEigen: max-err =%8.2e (between WaveHoltz and known solution, all frequencies)\n",maxErr);

    //       if( saveCheckFile )
    //       {
    //         Real solutionNorm = maxNorm( v ); 
    //         cgWaveHoltz.saveCheckFile( checkFileCounter,maxErr,solutionNorm ); 
    //         checkFileCounter++; 
    //       }
    //     }
    //  }

    //   saveCheckFile=false;

    // }

        if( answer=="contour" || (replot && plotChoices & 2) )
        {
      // CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
      // printF("plot contours...\n");

      // realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
      // ps.erase();
      // if( answer=="contour" )
      //   psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false); // wait inside contour

      // // const bool useUpwind = ad4>0.;
      // if( numberOfFrequencies==1 )
      //   psp.set(GI_TOP_LABEL,sPrintF("CgWaveHoltz: FD%i%i%s omega=%.5g",orderOfAccuracy,orderOfAccuracyInTime,(upwind ? "u" : ""),omega));
      // else
      //   psp.set(GI_TOP_LABEL,sPrintF("CgWaveHoltz: FD%i%i%s numFreq=%d",orderOfAccuracy,orderOfAccuracyInTime,(upwind ? "u" : ""),numberOfFrequencies));

      // PlotIt::contour(ps,v,psp);

      // psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

      // plotChoices |= 2; 
        }    

    }
    
    ps.popGUI();  // pop dialog

    if( showFileIsOpen )
        show.close();

    return 0;
}



