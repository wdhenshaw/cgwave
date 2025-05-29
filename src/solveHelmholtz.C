// This file automatically generated from solveHelmholtz.bC with bpp.
// -------------------------------------------------------
// ----  Solve the Helmholtz equation using WaveHoltz ----
// -------------------------------------------------------

#include "CgWaveHoltz.h"
#include "CgWave.h"
// #include "Oges.h"
#include "ParallelUtility.h"
// #include "CompositeGridOperators.h"
#include "PlotStuff.h"
#include "display.h"
#include "gridFunctionNorms.h"
#include "Ogshow.h"  
#include "ShowFileReader.h"
#include "InterpolatePointsOnAGrid.h"
#include "GridStatistics.h"


#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  


#define ForBoundary(side,axis)   for( int axis=0; axis<cg.numberOfDimensions(); axis++ ) for( int side=0; side<=1; side++ )
// ============================================================================================
/// Macro: compute errors in the Helmholtz solution
// ============================================================================================


static Real factorial( int n )
{
    Real val=1.;
    for( int k=2; k<=n; k++ )
    {
        val*=k;
    }
    return val;
}

static Real nChooseK( int n, int k )
{
    Real val;

    val = factorial(n)/( factorial(k) * factorial(n-k) );
    return val;
}

// ============================================================================================
/// \brief Solve the Helmholtz equation using WaveHoltz
// ============================================================================================
int CgWaveHoltz::solveHelmholtz(int argc,char **argv)
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

    int & useAugmentedGmres                  = cgWave.dbase.get<int>("useAugmentedGmres"); 
    int & augmentedVectorsAreEigenvectors    = cgWave.dbase.get<int>("augmentedVectorsAreEigenvectors");  // 1 = augmented vectors are true discrete eigenvectors
    int & deflateWaveHoltz                   = cgWave.dbase.get<int>("deflateWaveHoltz"); 
    int & numToDeflate                       = cgWave.dbase.get<int>("numToDeflate");
    aString & eigenVectorFile                = cgWave.dbase.get<aString>("eigenVectorFile"); 
    const int & eigenVectorsAreOrthogonal    = cgWave.dbase.get<int>("eigenVectorsAreOrthogonal"); 
    int & onlyLoadDeflatedEigenVectors       = cgWave.dbase.get<int>("onlyLoadDeflatedEigenVectors");
    int & plotFilterAndDeflatedEigenvalues   = cgWave.dbase.get<int>("plotFilterAndDeflatedEigenvalues");

    const int & numCompWaveHoltz             = cgWave.dbase.get<int>("numCompWaveHoltz");
    const int & filterTimeDerivative         = cgWave.dbase.get<int>("filterTimeDerivative");

    const int & useSuperGrid                 = cgWave.dbase.get<int>("useSuperGrid");
    const int & solveForScatteredField       = cgWave.dbase.get<int>("solveForScatteredField");

    cgWaveHoltz.dbase.get<int>("orderOfAccuracy")=orderOfAccuracy; // set value in CgWaveHoltz

    const int & numberOfFrequencies          = cgWaveHoltz.dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray         = cgWaveHoltz.dbase.get<RealArray>("frequencyArray");

    real & omega                             = cgWaveHoltz.dbase.get<real>("omega");
    Real & cpuWaveHoltz                      = cgWaveHoltz.dbase.get<Real>("cpuWaveHoltz");
    int & numPeriods                         = cgWaveHoltz.dbase.get<int>("numPeriods");
    real & tol                               = cgWaveHoltz.dbase.get<real>("tol");
    int & adjustOmega                        = cgWaveHoltz.dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
    int & maximumNumberOfIterations          = cgWaveHoltz.dbase.get<int>("maximumNumberOfIterations");
    int & numberOfIterations                 = cgWaveHoltz.dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken
    int & cgWaveDebugMode                    = cgWaveHoltz.dbase.get<int>("cgWaveDebugMode");
    

    onlyLoadDeflatedEigenVectors = true; // *new way* Sept 26, 2024


    int plotTotalField = false; 


    Index I1,I2,I3;

    aString nameOfShowFile = "cgWaveHoltz.show";
    bool showFileIsOpen=false;

    aString initialConditionShowFile = "initialConditionShowFile";

    bool plotOption=true;  // by default we plot interactively
    GL_GraphicsInterface & ps = (GL_GraphicsInterface&)(*Overture::getGraphicsInterface("cgwh",plotOption,argc,argv));

    PlotStuffParameters psp;  

    int & adjustErrorsForSuperGrid = cgWave.dbase.get<int>("adjustErrorsForSuperGrid");
    int & adjustPlotsForSuperGrid  = cgWave.dbase.get<int>("adjustPlotsForSuperGrid");
  // adjustPlotsForSuperGrid  = 1;
  // adjustErrorsForSuperGrid = 1;

  // Build a dialog menu for changing parameters
    GUIState gui;
    DialogData & dialog=gui;

    dialog.setWindowTitle("CgWaveHoltz - Helmholtz Solver");
    dialog.setExitCommand("exit", "exit");

  // dialog.setOptionMenuColumns(1);

  // aString accuracyLabel[] = {"second order", "fourth order", "" };
  // dialog.addOptionMenu("accuracy:", accuracyLabel, accuracyLabel, (orderOfAccuracy==2 ? 0 : 1) );

  // aString initialConditionLabel[] = {"smooth pulse", "pulse", "" };
  // dialog.addOptionMenu("Initial Condition:",initialConditionLabel,initialConditionLabel,(int)option );

  // enum TypeOfApproximation
  // {
  //   nonConservative=0,
  //   conservative=1
  // } typeOfApproximation=nonConservative;
        
  // aString approximationLabel[] = {"nonconservative", "conservative", "" };
  // dialog.addOptionMenu("Approximation:",approximationLabel,approximationLabel,(int)typeOfApproximation );

    aString pbLabels[] = {
                                                "compute with fixed-point",
                                                "compute with krylov",
                                                "solve Helmholtz directly",
                                                "zero initial condition",
                                                "show file initial condition",
                                                "random initial condition",
                                                "eigenvector initial condition",
                                                "helmholtz initial condition",
                                                "change parameters",
                                                "contour",
                                                "grid",
                                                "plot residual",
                        // "compute errors",
                                                "run cgWave and plot",
                        // "save to show",
                                                "save show file",
                                                "plot errors",
                                                "plot forcing",
                                                "plot sequences",
                                                "plot filter",
                                                "animate",
                                                "erase",
                                                "exit",
                                                ""};
    int numRows=12;
    dialog.setPushButtons( pbLabels, pbLabels, numRows ); 

  // aString tbCommands[] = {"save show file",
  //                          ""};
  // int tbState[10];
  // tbState[0] = saveShowFile==true;
  // int numColumns=1;
  // dialog.setToggleButtons(tbCommands, tbCommands, tbState, numColumns); 

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

    textCommands[nt] = "initial condition show file";  textLabels[nt]=textCommands[nt];
    sPrintF(textStrings[nt], "%s",(const char*)initialConditionShowFile);  nt++; 


  // null strings terminal list
    textCommands[nt]="";   textLabels[nt]="";   textStrings[nt]="";  assert( nt<numberOfTextStrings );
    dialog.setTextBoxes(textCommands, textLabels, textStrings);

    
    aString tbCommands[] = {
                                                    "run cgWave with debugging",
                                                    "deflate",
                                                    "adjust plots for superGrid",
                                                    "adjust errors for superGrid",
                                                    "plot total field",
                                                    "use augmented gmres",
                                                    "augmented are eigenvectors",
                                                    "plot filter and deflated eigs",
                                                    "only load deflated eigs",
                                                        ""};
    int tbState[10];
    tbState[0] = cgWaveDebugMode;
    tbState[1] = deflateWaveHoltz;
    tbState[2] = adjustPlotsForSuperGrid;
    tbState[3] = adjustErrorsForSuperGrid;
    tbState[4] = plotTotalField;
    tbState[5] = useAugmentedGmres;
    tbState[6] = augmentedVectorsAreEigenvectors;
    tbState[7] = plotFilterAndDeflatedEigenvalues;
    tbState[8] = onlyLoadDeflatedEigenVectors; 

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
    bool saveCheckFile=false; // set to true when errors should be saved to the check file.
    int showFileSolution=-1;  // solutions from different iterations can go in the show file 

  // int checkFileCounter=0; // keeps track of how many times the checkFile is saved with results

  // uHelmholtz holds Helmholtz solution from direct solver
  //   This is optionally used in CgWave as a known solution
    if( !cgWave.dbase.has_key("uHelmholtz") )
        cgWave.dbase.put<realCompositeGridFunction>("uHelmholtz");

    realCompositeGridFunction & uHelmholtz = cgWave.dbase.get<realCompositeGridFunction>("uHelmholtz"); 
    Range all;
    uHelmholtz.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz); // save Helmholtz solution here

    realCompositeGridFunction uh;  // holds complex WaveHoltz solution 
    if( filterTimeDerivative )
    {
        uh.updateToMatchGrid(cg,all,all,all,2);  // holds complex solution
        uh.setName("ur",0);   // real part 
        uh.setName("ui",1);   // imag part
        uh=0.;
    }  

    bool helmholtzFromDirectSolverWasComputed=false;
    bool waveHoltzSolutionWasComputed=false;
    Real cpuSolveHelmholtz=-1.;

    Ogshow show;  // show file for saving solutions

    int current=0;
    for(;;) 
    {
        ps.getAnswer(answer,"");      
        if( answer=="exit" || answer=="continue" )
        {
            break;
        }
        else if( answer=="change parameters" )
        {
            cgWaveHoltz.interactiveUpdate();
        }
        
        else if( dialog.getTextValue(answer,"omega","%e",omega) )
        {
            printF("Setting omega=%g\n",omega);
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
            printF("Setting maximumNumberOfIterations=%i\n",maximumNumberOfIterations);
        }
        else if( dialog.getTextValue(answer,"number to deflate","%i",numToDeflate) )
        {
            printF("Setting numToDeflate=%i\n",numToDeflate);
            cgWave.reinitializeDeflation();
        }  
        else if( answer == "show file initial condition" ) // put before "show file"
        {   
            realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
            CompositeGrid & cg = *v.getCompositeGrid();

      // aString initialConditionShowFile;
      // gi.inputString(initialConditionShowFile,"Enter the name of the show file for the initial condition)");

            printF("---- Read show file=[%s] for the initial condition -----\n",(const char*)initialConditionShowFile);
            ShowFileReader showFileReader;
            showFileReader.open(initialConditionShowFile);

            int numFrames = showFileReader.getNumberOfFrames();
            const int numberOfSolutions = showFileReader.getNumberOfSolutions();
            printF("There are %d solutions and %d frames\n",numberOfSolutions,numFrames);     

            realCompositeGridFunction q;
            int solutionNumber  = 1; // use the first solution
            CompositeGrid cgsf;
            showFileReader.getASolution(solutionNumber,cgsf,q); 

      // new way: 
            bool gridsMatch = CgWave::compositeGridsMatch( cg, cgsf );

            if( gridsMatch )
                printF("CgWaveHoltz:read IC from show file: Show file grid seems to match with current grid\n");
            else
                printF("CgWaveHoltz:read IC from show file: Show file grid DOES NOT match with current grid... will interpolate...\n");

            if( numberOfFrequencies>1 )
            {
                printF("ERROR: read IC from a show file: finish me for numberOfFrequencies>1\n");
                OV_ABORT("ERROR");
            }

            if( gridsMatch )
            {
                Index I1,I2,I3;
                int icomp=0; // use this component
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cg[grid];
                    getIndex(mg.dimension(),I1,I2,I3);
                    OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
                    OV_GET_SERIAL_ARRAY(real,q[grid],qLocal);        
                    bool ok=ParallelUtility::getLocalArrayBounds(q[grid],qLocal,I1,I2,I3);
                    if( ok )
                        vLocal(I1,I2,I3) = qLocal(I1,I2,I3,icomp);
                } 
            }
            else 
            { 
        // This way which should handle all cases and parallel 

                Range C=numberOfFrequencies;

                InterpolatePointsOnAGrid interp;

                CompositeGrid & cgsf = *q.getCompositeGrid();
        // cgsf.update(MappedGrid::THEmask | MappedGrid::THEvertex | MappedGrid::THEinverseVertexDerivative );
                cgsf.update(MappedGrid::THEmask );
        // cg.update(MappedGrid::THEmask | MappedGrid::THEvertex | MappedGrid::THEinverseVertexDerivative );        

                interp.setAssignAllPoints(true);  // assign all points -- extrap if necessary
                int width = orderOfAccuracy+1;
                interp.setInterpolationWidth( width ); // set interp width 
                interp.interpolateAllPoints( q,v, C, C );  // interpolate v from q
            }


            replot=true;

            cgWave.resetTimings(); // reset CPU timings to zero 

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
        else if( dialog.getTextValue(answer,"initial condition show file","%s",initialConditionShowFile) )
        {
            printF("Setting initialConditionShowFile=[%s]\n",(const char*)initialConditionShowFile);
        }         
        else if( dialog.getToggleValue(answer,"deflate",deflateWaveHoltz) )
        {
            printF("Setting deflateWaveHoltz=%d\n",deflateWaveHoltz);
        }
        else if( dialog.getToggleValue(answer,"use augmented gmres",useAugmentedGmres) )
        {
            printF("Setting useAugmentedGmres=%d\n",useAugmentedGmres);
        }  
        else if( dialog.getToggleValue(answer,"augmented are eigenvectors",augmentedVectorsAreEigenvectors) )
        {
            printF("Setting augmentedVectorsAreEigenvectors=%d ( 1 = augmented vectors are true discrete eigenvectors.)\n",augmentedVectorsAreEigenvectors);
        }         
        
        else if( dialog.getToggleValue(answer,"plot filter and deflated eigs",plotFilterAndDeflatedEigenvalues) )
        {
            printF("Setting plotFilterAndDeflatedEigenvalues=%d ( 1 = plot WaveHoltz filter beta and eigenvalues after delation is initialized.)\n",plotFilterAndDeflatedEigenvalues);
        }

        else if( dialog.getToggleValue(answer,"only load deflated eigs",onlyLoadDeflatedEigenVectors) )
        {
            printF("Setting onlyLoadDeflatedEigenVectors=%d ( 1 = only read in eigenvectors used for deflation.)\n",onlyLoadDeflatedEigenVectors);
        }    

        else if( dialog.getToggleValue(answer,"adjust plots for superGrid",adjustPlotsForSuperGrid) )
        {
            printF("Setting adjustPlotsForSuperGrid=%d (zero out solution in SG layers when plotting)\n",adjustPlotsForSuperGrid);
            replot=true;
        }

        else if( dialog.getToggleValue(answer,"plot total field",plotTotalField) )
        {
            printF("Setting plotTotalField=%d (1= plot total field for scattering problems)\n",plotTotalField);
            replot=true;
        }

        else if( dialog.getToggleValue(answer,"run cgWave with debugging",cgWaveDebugMode) )
        {
            printF("Setting cgWaveDebugMode=%d (1: run cgWave plotting every step)\n",cgWaveDebugMode);
            if( cgWaveDebugMode==1 )
            {

            }
        }          

        else if( answer=="compute with fixed-point" ||
                          answer=="compute with krylov" || 
                          answer=="compute with petsc"       // backward compatible
                      )
        {
            waveHoltzSolutionWasComputed=true; 

            int & useFixedPoint = dbase.get<int>("useFixedPoint"); 
            useFixedPoint = answer=="compute with fixed-point";

            const int & computeEigenmodes   = cgWave.dbase.get<int>("computeEigenmodes");
            const int readTrueEigenPairs    = eigenVectorFile != "none";
            if( computeEigenmodes && readTrueEigenPairs )
            {
        // -- This next call will read in any known eigenmodes ---
                cgWave.initializeDeflation();
            }

            const Real cpu0=getCPU();


      // -----------------------------------------------------------------
      // ------------------WAVE HOLTZ SOLVE-------------------------------
      // -----------------------------------------------------------------
            if( useFixedPoint )
            {
                cgWaveHoltz.solve();
            }
            else
            {
                if( useAugmentedGmres )
                    cgWaveHoltz.solveAugmentedKrylov(argc,argv);
                else 
                    cgWaveHoltz.solvePETSc(argc,argv);
            }

      // {
      //   if( computeEigenmodes )
      //     cgWaveHoltz.solveSLEPc(argc,argv);
      //   else
      //     cgWaveHoltz.solvePETSc(argc,argv);
      // }


            if( deflateWaveHoltz && useAugmentedGmres==0 )
                cgWave.inflateSolution(); // un-deflate the solution

      // const Real cpu = getCPU()-cpu0;
            cpuWaveHoltz = getCPU()-cpu0;

            Real viFactor; // set below
            if( filterTimeDerivative )
            {
        // -- compute the complex solution ---
                realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");

                Real om = frequencyArrayAdjusted(0);

                const Real & dt = cgWave.dbase.get<real>("dtUsed");  // check me
                Real sigma = -1.;       // SIGN FOR exp( sigma*I*omega*t)    ** FIX ME ***   MUST MATCH VALUE IN solveHelmholtzDirect

        // **NOTE**
        // viFactor must match in
        //     getFilterWeights
        //     takeFirstStep
        //     solveHelmholtz : definition of Im(u)
        // viFactor must match 1/omegaTilde in the filterWeight definition AND TAKE FIRST STEP


                const CgWave::TimeSteppingMethodEnum & timeSteppingMethod = cgWave.dbase.get<CgWave::TimeSteppingMethodEnum>("timeSteppingMethod");
                if( true )
                {
                    viFactor = cgWave.dbase.get<Real>("viFactor");  // *new way* June 19, 2023
                }
                else if( true )
                {
                    viFactor = -sigma/om;  // use adjusted omega
                }
                else if( true || timeSteppingMethod==CgWave::implicitTimeStepping )
                {
                    viFactor = -sigma*dt/sin(om*dt); // adjust for D0t -- *check me*
                }
                else
                {
                    viFactor = -sigma/frequencyArray(0);        // unadjusted
                }

                printF("\n ###### solveHelmholtz: Set complex solution using viFactor=%9.2e#### \n\n",viFactor);

                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cg[grid];
                    getIndex(mg.dimension(),I1,I2,I3);
                    OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                    OV_GET_SERIAL_ARRAY(Real,uh[grid],uhLocal);

                    uhLocal(I1,I2,I3,0) =          vLocal(I1,I2,I3,0);
                    uhLocal(I1,I2,I3,1) = viFactor*vLocal(I1,I2,I3,1);
                }        
            }

  

      // The period array may have changed since numPeriods(freq) may have changed
            cgWaveHoltz.dbase.get<RealArray>("periodArray") = cgWave.dbase.get<RealArray>("periodArray");       

      // total number of steps taken:
            const int & numberOfStepsTaken     = cgWave.dbase.get<int>("numberOfStepsTaken");      
            const int & numberOfStepsPerSolve  = cgWave.dbase.get<int>("numberOfStepsPerSolve");
            const Real & cfl                   = cgWave.dbase.get<real>("cfl");
            const Real & c                     = cgWave.dbase.get<real>("c");
            const Real & damp                  = cgWave.dbase.get<real>("damp");
            const Real & dt                    = cgWave.dbase.get<real>("dtUsed");
            const int & minStepsPerPeriod      = cgWave.dbase.get<int>("minStepsPerPeriod");

            const RealArray & timing           = cgWave.timing;


            const Real & maxResidual          = cgWaveHoltz.dbase.get<real>("maxResidual");
            const aString & krylovType        = cgWaveHoltz.dbase.get<aString>("krylovType");  // gmres, bicgstab, cg, ...
            const int & gmresRestartLength    = cgWaveHoltz.dbase.get<int>("gmresRestartLength"); // restart length for WaveHoltz + GMRES  

            const CgWave::TimeSteppingMethodEnum & timeSteppingMethod = cgWave.dbase.get<CgWave::TimeSteppingMethodEnum>("timeSteppingMethod");
            aString implicitSolverName =  timeSteppingMethod==CgWave::explicitTimeStepping ? 
                                                                        "none" : 
                                                                        cgWave.dbase.get<aString>("implicitSolverName");

            const Real aveItsPerImplicitSolve               = cgWave.getAverageNumberOfIterationsPerImplicitSolve();
            const Real & waveHoltzAsymptoticConvergenceRate = cgWave.dbase.get<real>("waveHoltzAsymptoticConvergenceRate");

            const Real cpuWave   = timing(CgWave::totalTime);
            const Real cpuAdvance = timing(CgWave::timeForAdvance);

      // -- Compute points-per-wavelength --- 
            RealArray & dxMinMax = cgWave.dbase.get<RealArray>("dxMinMax");
            Real  kWaveNumber = frequencyArray(0)/c;     // k = omega/c
            Real lambdaWaveLength = twoPi/kWaveNumber;   // lambda = 2*pi/k
            const int gridBackGround=0; 
            Real dx = dxMinMax(gridBackGround,0);        // grid spacing from backGround grid (hopefully)
            Real & ppw            = dbase.get<Real>("ppw");            // holds actual points per wavelength
            Real & ppwRuleOfThumb = dbase.get<Real>("ppwRuleOfThumb"); // holds rule-of-thumb points per wavelength
            ppw = lambdaWaveLength/dx;        // PPW = lambda/dx 

            Real xMin[3],xMax[3];
            GridStatistics::getGridCoordinateBounds(cg,xMin,xMax);
            Real domainSize=0.;
            for( int axis=0; axis<cg.numberOfDimensions(); axis++ )
                domainSize = max(domainSize,xMax[axis]-xMin[axis]);
            const Real NLambda = domainSize/lambdaWaveLength; // size of domain in wavelengths
            const Real epsr=1e-2; // relative error tolerance
      //  bp = 2/( (mu+1)^2 *nchoosek(2*mu+2,mu+1) );
            const int mu = orderOfAccuracy/2; 
            Real bp = 2./( (mu+1)*(mu+1) * nChooseK(2*mu+2,mu+1) );
            ppwRuleOfThumb = 2.*Pi*pow( Pi*bp* NLambda/epsr,1./orderOfAccuracy); 
            Real ppwRuleOfThumbOld = Pi*pow( NLambda/epsr,1./orderOfAccuracy); // ** FIX ME
            printF("\n >>>omega=%10.2e, k=%10.2e, lambda=%10.2e, dx=%10.2e, domainSize=%10.2e, NLambda=%g, ppw=lambda/dx=%5.1f, ppw(rule-of-thumb)=%5.1f (old=%5.1f) (epsr=%g)<<<\n\n",
                  kWaveNumber,frequencyArray(0),lambdaWaveLength,dx,domainSize,NLambda,ppw,ppwRuleOfThumb,ppwRuleOfThumbOld,epsr);

            printF("\n------------------------ CgWaveHoltz SUMMARY ------------------------------------\n");
            if( filterTimeDerivative==0 )
                printF("                          Real Solution \n");
            else
                printF("                          Complex Solution \n");

            const aString & nameOfOGFile = cgWave.dbase.get<aString>("nameOfGridFile");
            printF("\n grid=%s\n",(const char*)nameOfOGFile);

            printF(" Convergence rate CR=%5.3f, CR-perPeriod=%5.3f, Time to solve =%9.2e(s) : CgWave=%9.2e(s) %5.2f%%, advance=%9.2e(s) %5.2f%%\n",
                          convergenceRate,convergenceRatePerPeriod,cpuWaveHoltz,cpuWave,(cpuWave/cpuWaveHoltz)*100.,cpuAdvance,(cpuAdvance/cpuWaveHoltz)*100.);
            printF(" max-residual=%8.2e (WaveHoltz its), numberOfIterations = %d (%s) (maximumNumberOfIterations=%d)\n",
                              maxResidual,numberOfIterations,(useFixedPoint? "FP" : (useAugmentedGmres ? "Augmented-GMRES" : (const char*)krylovType)  ),maximumNumberOfIterations);
            
            printF(" omega=%10.2e, k=%10.2e, lambda=%10.2e, dx=%10.2e, domainSize=%10.2e, NLambda=%g, ppw=lambda/dx=%5.1f, ppw(rule-of-thumb)=%5.1f (epsr=%g)<<<\n\n",
                kWaveNumber,frequencyArray(0),lambdaWaveLength,dx,domainSize,NLambda,ppw,ppwRuleOfThumb,epsr);      
            printF(" numberOfFrequencies= %d, adjustOmega=%d, adjustHelmholtzForUpwinding=%d\n",
                        numberOfFrequencies,adjustOmega,adjustHelmholtzForUpwinding);
            printF(" solveForScatteredField=%d\n",solveForScatteredField);
            if( solveForScatteredField )
            {
                const Real & amp   = cgWave.dbase.get<Real>("ampPlaneWave");
                const Real & kx    = cgWave.dbase.get<Real>("kxPlaneWave");
                const Real & ky    = cgWave.dbase.get<Real>("kyPlaneWave");
                const Real & kz    = cgWave.dbase.get<Real>("kzPlaneWave");
                const Real & phi   = cgWave.dbase.get<Real>("phiPlaneWave");
                const Real & omega = cgWave.dbase.get<Real>("omegaPlaneWave");
                printF(" Scattered field computed with plane-wave: amp=%g, [kx,ky,kz]=[%g,%g,%g]*2*pi, phi=%g\n",amp,kx/twoPi,ky/twoPi,kz/twoPi,phi);        
            }
            printF(" Time stepping = %s, \n",
                        (timeSteppingMethod==CgWave::explicitTimeStepping ? "explicit" : "implicit") );
            printF(" cfl=%6.2f, dt=%9.3e, minStepsPerPeriod=%d, numberOfStepsPerSolve=%d, totalTimeStepsStaken=%d, aveItsPerImplicitSolve=%5.1f\n",
                        (timeSteppingMethod==CgWave::explicitTimeStepping ? "explicit" : "implicit"),cfl,dt,minStepsPerPeriod,numberOfStepsPerSolve,numberOfStepsTaken,
                        aveItsPerImplicitSolve);

            const int & takeImplicitFirstStep         = cgWave.dbase.get<int>("takeImplicitFirstStep");
            printF("  takeImplicitFirstStep=%d\n",takeImplicitFirstStep);
            const int & implicitUpwind                = cgWave.dbase.get<int>("implicitUpwind");
            printF(" implicitUpwind=%d (1= add upwind diss. to implicit matrix)\n",implicitUpwind);
            if( timeSteppingMethod==CgWave::implicitTimeStepping )
            {
                RealArray & gridCFL = cgWave.dbase.get<RealArray>("gridCFL"); // really c/dx 
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    printF(" grid-CFL=%9.2e for grid=%d (%s)\n",dt*gridCFL(grid),grid,(const char*)cg[grid].getName());
                }
            }

            printF(" c=%g, useSuperGrid=%d, damp=%g\n",c,useSuperGrid,damp);
            printF(" gmres-restart-length=%d\n",gmresRestartLength);

            const int & filterTimeDerivative  = cgWave.dbase.get<int>("filterTimeDerivative");
            const int & useFilterWeights      = cgWave.dbase.get<int>("useFilterWeights");
            const int & filterD0t             = cgWave.dbase.get<int>("filterD0t");
            const int & deflateForcing        = cgWave.dbase.get<int>("deflateForcing");
  
            printF(" filterTimeDerivative=%d, useFilterWeights=%d, viFactor=%g, filterD0t=%d\n",filterTimeDerivative,useFilterWeights,viFactor,filterD0t);

            printF("    implicit solver : %s\n =====\n",(const char*)implicitSolverName);
            printF(" deflateWaveHoltz=%d, deflateForcing=%d, numToDeflate=%d, eigenVectorFile=[%s]\n",deflateWaveHoltz,deflateForcing,numToDeflate,
                    (const char*)eigenVectorFile);
            if( waveHoltzAsymptoticConvergenceRate>0 )
                printF(" Estimated asymptotic convergence rate=%6.3f.\n",waveHoltzAsymptoticConvergenceRate);

            for( int freq=0; freq<numberOfFrequencies; freq++ )
            {
                printF(" freq=%2d, omega=%8.3f",freq,frequencyArray(freq));
                if( adjustOmega )
                    printF(", (adjusted=%8.3f)",frequencyArrayAdjusted(freq));

                const Real T = twoPi/frequencyArray(freq);  // period
                const Real Tbar  = T*numPeriodsArray(freq);       // periods that fit in time interval
                const Real Tbar0 = periodArray(0);          // final time integrated to 
                printF(" T=%8.5f, Tbar=%8.5f, numPeriods=%3d, Tbar(0)/T=%6.2f",T,Tbar,numPeriodsArray(freq),Tbar0/T);
                printF("\n");
            }

            if( cpuSolveHelmholtz>0 )
            {
                int & helmholtzSolverIterations = cgWaveHoltz.dbase.get<int>("helmholtzSolverIterations");
                aString & helmholtzSolverName =  cgWaveHoltz.dbase.get<aString>("helmholtzSolverName");
                printF("Direct solve Helmholtz: cpu=%9.2e(s), iterations=%d (sum all freq)\n",cpuSolveHelmholtz,helmholtzSolverIterations);
                printF("   solver : %s\n",(const char*)helmholtzSolverName);
            }   
            printF("\n");

            Real & maxRes = dbase.get<Real>("maxRes");

            RealArray maxResArray;
            const int useAdjustedOmega=0; // compute residual using original (unadjusted) omega
            if( filterTimeDerivative==0 )
            {
                maxRes = cgWaveHoltz.residual( maxResArray,useAdjustedOmega );
            }
            else
            {
        // maxRes = cgWaveHoltz.residual( uHelmholtz );
                maxRes = cgWaveHoltz.residual( uh,maxResArray,useAdjustedOmega );
            }

            aString label = deflateWaveHoltz ? sPrintF("deflate=%d",numToDeflate) : "";
            aString schemeName;
            if( useFixedPoint )
                schemeName="FP";
            else
            {
                if( useAugmentedGmres )
                    schemeName = "a" + krylovType; // Augmented 
                else
                    schemeName = krylovType;
            }
            if( filterTimeDerivative==0 )
            {
                if( numberOfFrequencies==1 )
                {
                    printF("solveHelmholtz: max-res=%9.3e (%s, numIts=%d, ave-CR=%5.2f %s, ACR=%5.2f)\n",
                        maxRes,(const char*)schemeName,numberOfIterations,
                        convergenceRate,(const char*)label,waveHoltzAsymptoticConvergenceRate);
          // printF("solveHelmholtz: max-res=%9.3e (%s, numIts=%d, ave-CR=%5.2f %s, ACR=%5.2f)\n",
          //   maxRes,(useFixedPoint? "FP" : (useAugmentedGmres ? "Augmented-GMRES" : (const char*)krylovType)),numberOfIterations,
          //   convergenceRate,(const char*)label,waveHoltzAsymptoticConvergenceRate);          
                }

                else
                    printF("solveHelmholtz: max-res=%9.3e (all freq, %s, numIts=%d, ave-CR=%5.2f %s)\n",
                          maxRes,(const char*)schemeName,numberOfIterations,convergenceRate,(const char*)label );
            }
            else
            {
                printF("solveHelmholtz: max-residual [Re,Im]=[%8.2e,%8.2e] (%s, numIts=%d, ave-CR=%5.2f %s)\n",maxResArray(0),maxResArray(1),
                                (const char*)schemeName,numberOfIterations,convergenceRate,(const char*)label );
        // printF("solveHelmholtz: max-residual [Re,Im]=[%8.2e,%8.2e] (%s, numIts=%d, ave-CR=%5.2f %s)\n",maxResArray(0),maxResArray(1),
        //         (useFixedPoint? "FP" : (useAugmentedGmres ? "Augmented-GMRES" : (const char*)krylovType)),numberOfIterations,convergenceRate,(const char*)label );        
            }

      // --- COMPUTE ERROR COMPARED TO DIRECT HELMHOLTZ ---
            Real errorBetweenWaveHoltzAndHelmholtz=0.;
                if( helmholtzFromDirectSolverWasComputed )
                {
          // --- COMPUTE ERROR COMPARED TO DIRECT HELMHOLTZ ---
          // ** make this a separate function maybe ***
                    realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
                    realCompositeGridFunction & err = cgWave.dbase.get<realCompositeGridFunction>("error");
                    if( err.getComponentBound(0)!= (numCompWaveHoltz-1) )
                    {
                        err.updateToMatchGrid(cg,all,all,all,numCompWaveHoltz);
                        if( filterTimeDerivative==0 )
                        {
                            for( int freq=0; freq<numberOfFrequencies; freq++ )
                                err.setName(sPrintF("v%d-err",freq),freq);
                        }
                        else
                        {
                            err.setName("ur-err",0); // real part of error 
                            err.setName("ui-err",1); // imaginary part of the error
                        }
                    }
          // realCompositeGridFunction err(cg,all,all,all,numCompWaveHoltz);
                    Real maxDiff=0., vNormMax=0., solutionNorm=1.; 
                    if( filterTimeDerivative==0 )
                    {
                        RealArray vNorm(numberOfFrequencies);
                        if( useSuperGrid && adjustErrorsForSuperGrid )
                        {
              // stick v into err so it can be adjusted for superGrid without changing v
                            err = v;
                            cgWave.adjustSolutionForSuperGrid( err );
                            for( int freq=0; freq<numberOfFrequencies; freq++ )
                                vNorm(freq) = max( maxNorm( err,freq ), REAL_MIN*100);
                        } 
                        else
                        {
                            for( int freq=0; freq<numberOfFrequencies; freq++ )
                                vNorm(freq) = max( maxNorm( v,freq ), REAL_MIN*100);
                        }
                        err = uHelmholtz - v;
                        if( useSuperGrid && adjustErrorsForSuperGrid )
                        {
                            printF("INFO: Set errors to zero in the SuperGrid layers\n");
                            cgWave.adjustSolutionForSuperGrid( err );
                        }      
                        for( int freq=0; freq<numberOfFrequencies; freq++ )
                        {
              // const Real vNorm = max( maxNorm( v,freq ), REAL_MIN*100);
                            const Real diff  = maxNorm( err,freq ) / vNorm(freq);
                            if( numberOfFrequencies>1 )
                                printF("solveHelmholtz: rel-max-diff=%8.2e for freq%d = %12.4e, |v|_max=%9.2e (between WaveHoltz and Direct Helmholtz solution)\n",
                                          diff,freq,frequencyArray(freq),vNorm(freq));
                            maxDiff = max(maxDiff,diff);
                            vNormMax = max(vNormMax,vNorm(freq));
                        }
                        solutionNorm = vNormMax; // for check file
                        printF("solveHelmholtz: rel-max-diff=%8.2e, |v|_max=%9.2e (between WaveHoltz and Direct Helmholtz solution, all frequencies)\n",maxDiff,vNormMax);
                    }
                    else
                    {
            // --- complex case -----
                        const int numberOfComponents=2;
                        Range R2=numberOfComponents;
                        RealArray uNorm(R2);
                        if( useSuperGrid && adjustErrorsForSuperGrid )
                        {
              // stick uh into err so it can be adjusted for superGrid without changing uh
                            err = uh;
                            cgWave.adjustSolutionForSuperGrid( err );
                            for( int ic=0; ic<numberOfComponents; ic++ )
                                uNorm(ic) = max( maxNorm( err,ic ), REAL_MIN*100);
                        } 
                        else
                        {
                            for( int ic=0; ic<numberOfComponents; ic++ )
                                uNorm(ic) = max( maxNorm( uh,ic ), REAL_MIN*100);
                        }
                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            MappedGrid & mg = cg[grid];
                            getIndex(mg.dimension(),I1,I2,I3);
                            OV_GET_SERIAL_ARRAY(Real,uh[grid],uhLocal);
                            OV_GET_SERIAL_ARRAY(Real,err[grid],errLocal);
                            OV_GET_SERIAL_ARRAY(Real,uHelmholtz[grid],uHelmholtzLocal);
                            errLocal(I1,I2,I3,R2) = uhLocal(I1,I2,I3,R2) - uHelmholtzLocal(I1,I2,I3,R2);
                        }
                        if( useSuperGrid && adjustErrorsForSuperGrid )
                        {
                            printF("INFO: Set errors to zero in the SuperGrid layers\n");
                            cgWave.adjustSolutionForSuperGrid( err );
                        }
            // const Real urNorm = max( maxNorm( uh,0 ), REAL_MIN*100);
            // const Real uiNorm = max( maxNorm( uh,1 ), REAL_MIN*100);
                        const Real uNormTotal = max(uNorm(0),uNorm(1));
                        const Real urDiff = maxNorm( err,0 ) / uNormTotal;
                        const Real uiDiff = maxNorm( err,1 ) / uNormTotal;
                        printF("solveHelmholtz: rel-max-diff [Re,Im]=[%8.2e,%8.2e],  norm [Re,Im]=[%8.2e,%8.2e] (between WaveHoltz and Direct Helmholtz)\n",
                                      urDiff,uiDiff, uNorm(0),uNorm(1));
                        maxDiff = max(urDiff,uiDiff);  // for check file 
                        solutionNorm = uNormTotal;     // for check file
                    }
                    errorBetweenWaveHoltzAndHelmholtz = maxDiff;
          // if( saveCheckFile )
          // {
          //   // Real solutionNorm = maxNorm ( v );
          //   cgWaveHoltz.saveCheckFile( checkFileCounter,maxDiff,solutionNorm ); 
          //   checkFileCounter++; 
          // }
                } // end compute Helmholtz errors 
            printF("\n");


      // if( 1==1 )  
      //   cgWave.printStatistics();  // @@@@@@@@@@@@@@@ TEST @@@@@@@@@@@@@@@@@@@

      // ---------------------------------------
      // --------- SAVE CHECK FILE -------------
      // ---------------------------------------

            if( true )
            {
                outputCheckFile( maxResArray, errorBetweenWaveHoltzAndHelmholtz );
            }
            else
            {
        // ** OLD WAY **

                aString timeSteppingName = (timeSteppingMethod==CgWave::explicitTimeStepping ? "explicit" : "implicit");
                
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
                FILE *& checkFile = dbase.get<FILE*>("checkFile");
                checkFile = fopen("cgWaveHoltz.check","w" );      
                assert( checkFile != NULL );      

        // Get the current date
                time_t *tp= new time_t;
                time(tp);
                const char *dateString = ctime(tp);
                fPrintF(checkFile,"# Check file for cgWaveHoltz, grid=%s, %s",(const char*)gridNameNoPrefix,dateString);  // Note: dateString include newline
                delete tp; 

                fPrintF(checkFile,"grid=%s;\n",(const char*)gridNameNoPrefix);
                fPrintF(checkFile,"timeStepping=%s;\n",(const char*)timeSteppingName);
                fPrintF(checkFile,"orderOfAccuracy=%d;\n",orderOfAccuracy);
                fPrintF(checkFile,"numPeriods=%d;\n",numPeriods);

                fPrintF(checkFile,"numberOfFrequencies=%d;\n",numberOfFrequencies);
                for( int freq=0; freq<numberOfFrequencies; freq++ )
                {
                    fPrintF(checkFile," freq=%2d, omega=%8.3f",freq,frequencyArray(freq));
                    if( adjustOmega )
                        fPrintF(checkFile,", (adjusted=%8.3f)",frequencyArrayAdjusted(freq));

                    const Real T = twoPi/frequencyArray(freq);  // period
                    const Real Tbar  = T*numPeriodsArray(freq);       // periods that fit in time interval
                    const Real Tbar0 = periodArray(0);          // final time integrated to 
                    fPrintF(checkFile," T=%8.5f, Tbar=%8.5f, numPeriods=%3d, Tbar(0)/T=%6.2f",T,Tbar,numPeriodsArray(freq),Tbar0/T);
                    fPrintF(checkFile,"\n");
                }

                fPrintF(checkFile,"filterTimeDerivative=%d;\n",filterTimeDerivative);

                fPrintF(checkFile,"useFixedPoint=%d;\n",useFixedPoint);
                fPrintF(checkFile,"useAugmentedGmres=%d;\n",useAugmentedGmres);
                fPrintF(checkFile,"minStepsPerPeriod=%d;\n",minStepsPerPeriod);
                fPrintF(checkFile,"numberOfStepsPerSolve=%d;\n",numberOfStepsPerSolve);

                fPrintF(checkFile,"convergenceRate=%5.3f;\n",convergenceRate);
                fPrintF(checkFile,"numberOfIterations=%d;\n",numberOfIterations);

                fPrintF(checkFile,"numToDeflate=%d;\n",numToDeflate);

                if( filterTimeDerivative==0 )
                    fPrintF(checkFile,"maxRes=%9.3e;\n",maxRes);
                else
                {
                    fPrintF(checkFile,"maxResReal=%8.2e;\n",maxResArray(0));
                    fPrintF(checkFile,"maxResImag=%8.2e;\n",maxResArray(1));
                }
                fPrintF(checkFile,"errorBetweenWaveHoltzAndHelmholtz=%9.3e;\n",errorBetweenWaveHoltzAndHelmholtz);

        // fPrintF(checkFile,"numberOfStepsPerSolve=%d;\n",numberOfStepsPerSolve);
        // fPrintF(checkFile,"numEigsRequested=%d;\n",numEigsToCompute);
        // fPrintF(checkFile,"numEigsComputed=%d;\n",numEigenVectors);
        // fPrintF(checkFile,"numArnoldiVectors=%d;\n",numArnoldiVectors);
        // fPrintF(checkFile,"numWaveSolves=%d;\n",numberOfMatrixVectorMultiplications);
        // fPrintF(checkFile,"maxEigErr=%9.2e;\n",maxEigErr);
        // fPrintF(checkFile,"maxEvectErr=%9.2e;\n",maxEvectErr);
        // fPrintF(checkFile,"maxEigResid=%9.2e;\n",maxEigResid);
                fclose(checkFile);

                printF("Wrote results to the check file [cgWaveHoltz.check]\n");
            }

      // save results to a matlab file
            const aString & matlabFileName = cgWaveHoltz.dbase.get<aString>("matlabFileName");  
            aString localName; 
            aString suffix = sPrintF("FD%d%dTS%s",orderOfAccuracyInTime,orderOfAccuracy,
                                                              timeSteppingMethod==CgWave::implicitTimeStepping ? "I" : "E");
            if( numberOfFrequencies>1 )
                suffix = suffix + sPrintF("Nf%d",numberOfFrequencies);

            suffix = suffix + sPrintF("Np%d",numPeriods);

            localName = matlabFileName + suffix;

            if( useFixedPoint )
            {
                cgWaveHoltz.dbase.get<aString>("solverName")="FixedPoint";
                localName = localName + "FP"; 
            }
            else
            {
        // aString krylovTypeName = useAugmentedGmres ? "AGmres" : krylovType; 
                aString krylovTypeName = krylovType; 
                if( useAugmentedGmres ) 
                    krylovTypeName = "a" + krylovTypeName; // Augmented

                cgWaveHoltz.dbase.get<aString>("solverName")=krylovTypeName;
                localName = localName + krylovTypeName;
            }
            cgWaveHoltz.outputMatlabFile( localName );


            replot=true;
            reComputeErrors=true;
        }

        else if( answer == "solve Helmholtz directly" )
        {
       // this may require cgWave to have been called to create f ?
      // uHelmholtz.updateToMatchGrid(cg); // save Helmholtz solution here

  

            realCompositeGridFunction & q = filterTimeDerivative ? uh : cgWave.dbase.get<realCompositeGridFunction>("v");
            realCompositeGridFunction & f = cgWave.dbase.get<realCompositeGridFunction>("f");

      // q.updateToMatchGrid(cg,all,all,all,numberOfFrequencies); // *WDH* Aug 20. 2024 WHY IS THIS NOW NEEDED ********************************

      // f.updateToMatchGrid(cg,all,all,all,numberOfFrequencies);

      // **NOTE** Make a copy of f since the direct solver will change it 
            realCompositeGridFunction fHelmholtz(cg,all,all,all,numberOfFrequencies);
            fHelmholtz = f; 

            const Real cpu0=getCPU();

            cgWaveHoltz.solveHelmholtzDirect( q, fHelmholtz  );

            cpuSolveHelmholtz = getCPU()-cpu0;

            uHelmholtz = q;
            helmholtzFromDirectSolverWasComputed=true;

            if( true && numberOfFrequencies==1 )
            {
        // bool useAdjustedOmega=false;
        // Real maxRes = cgWaveHoltz.residual( useAdjustedOmega );
        // printF("CgWaveHoltz: omega=%9.3e, max-res=%9.3e, cpu=%9.2e (direct solution of Helmholtz).\n",
        //        omega,maxRes,cpu);
    
                int useAdjustedOmega=0; // do not adjust for omega
        // int useAdjustedOmega=2; // do not adjust for omega

        // Fill in the forcing and boundary conditions: 
        // cgWave.getHelmholtzForcing( fHelmholtz );
                printF("solveHelmholtz: call residual to check the solution from the direct Helmholtz solver...\n");
                RealArray maxResArray;
                Real maxRes = cgWaveHoltz.residual( uHelmholtz, fHelmholtz, maxResArray, useAdjustedOmega );
                printF("CgWaveHoltz: omega=%9.3e, max-res=%9.3e (direct Helmholtz solution)\n",omega,maxRes);
            }

            reComputeErrors=true;
            replot=true;

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
            realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
            v=0.;

      // v.display("v after v=0","%5.2f ");

            cgWave.resetTimings(); // reset CPU timings to zero 
        }
        else if( answer == "helmholtz initial condition" )
        {
            if( helmholtzFromDirectSolverWasComputed )
            {
                printF("Setting initial condition to the solution from the direct Helmholtz solver.\n");
                realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");

                v=uHelmholtz;  // this is ok for complex if viFactor=1 (now default)

        // if( filterTimeDerivative==0 )
        // {
        //   v=uHelmholtz;
        // }
        // else
        // {
        //   // Complex case:
        //   //  u = ur*cos(omega*t) + ui*sin(omega*t)
        //   // v(0) = ur
        //   // v(1) = D_0t u(0) = omega*ui 
        //   Real omega=frequencyArray(0);
        //   for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        //   {
        //     MappedGrid & mg = cg[grid];
        //     OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
        //     OV_GET_SERIAL_ARRAY(real,uHelmholtz[grid],uLocal);
        //     getIndex(mg.dimension(),I1,I2,I3);
        //     vLocal(I1,I2,I3,0) = uLocal(I1,I2,I3,0);
        //     vLocal(I1,I2,I3,1) = uLocal(I1,I2,I3,1)*omega; // could do better 
        //   }
        // }


                cgWave.resetTimings(); // reset CPU timings to zero 
            }
            else
            {
                printF("The direct Helmholtz solution has not been computed yet.\n");
            }
        }    

        else if( answer == "random initial condition" )
        {
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
                bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
                for( int freq=0; freq<numCompWaveHoltz; freq++ )
                {
                    FOR_3D(i1,i2,i3,I1,I2,I3)
                    {
            // do this for now
            // vLocal(i1,i2,i3) = sin(i1)*cos(i2);
                        if( maskLocal(i1,i2,i3)!=0 )
                        {
                            vLocal(i1,i2,i3,freq) = (-1. + std::rand()*(2./RAND_MAX))*scale; // [-1,1]
                        }

                    }
                }
            }
            replot=true;

            cgWave.resetTimings(); // reset CPU timings to zero 
        }  
        else if( answer == "eigenvector initial condition" )
        {
            realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
            CompositeGrid & cg = *v.getCompositeGrid();

      // --- get eigenvectors
            cgWave.initializeDeflation();

      // Here are the eigenvectors and eigenvalues:
            realCompositeGridFunction & uev          = cgWave.dbase.get<realCompositeGridFunction>("uev");
            RealArray & eig                          = cgWave.dbase.get<RealArray>("eig");
            IntegerArray & eigNumbersToDeflate       = cgWave.dbase.get<IntegerArray>("eigNumbersToDeflate");
            const int & onlyLoadDeflatedEigenVectors = cgWave.dbase.get<int>("onlyLoadDeflatedEigenVectors");

            const int ie=0;
            const int je = onlyLoadDeflatedEigenVectors ? eigNumbersToDeflate(ie) : ie; 
            printF("Choosing eigenvector=%d: eig=%12.4e\n",ie,eig(0,je));

  
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                MappedGrid & mg = cg[grid];
                OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
                OV_GET_SERIAL_ARRAY(real,uev[grid],uevLocal);

                getIndex(mg.dimension(),I1,I2,I3);
                bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
                if( ok )
                {
                    for( int freq=0; freq<numberOfFrequencies; freq++ )
                    {
                        vLocal(I1,I2,I3,freq) = uevLocal(I1,I2,I3,ie);
                    }
                }
            }
            replot=true;      

            printF("cgwh:TESTING: call cgWave.updateEigenmodes ...\n");
            cgWave.updateEigenmodes();
      // OV_ABORT("stop here for now");
        }      
        else if( answer=="save show file" || 
                          answer=="save to show" )
        {
            printF("Save the current solution, errors and forcings to the show file [%s].\n",(const char*)nameOfShowFile);
            if( !showFileIsOpen )
            {
                showFileIsOpen=true;
                show.open(nameOfShowFile);
            }
            
            showFileSolution++; 
            if( showFileSolution==0 )
                show.saveGeneralComment("Solutions from CgWaveHoltz"); // save a general comment in the show file
      // show.saveGeneralComment(" file written on April 1");      // save another general comment

            show.startFrame();                       // start a new frame
      // const bool useUpwind = ad4>0.;
            show.saveComment(0,sPrintF("CgWaveHoltz: FD%i%i%s omega=%.5g",orderOfAccuracy,orderOfAccuracyInTime,(upwind ? "u" : ""),omega));   // comment 0 (shown on plot)
            if( useSuperGrid )
                show.saveComment(1,sPrintF("SuperGrid"));               // comment 1 (shown on plot)

      // We save the current solution and optionally the direct solution and or errors
      // Also save the forcing function

      // const int & computeEigenmodes     = cgWave.dbase.get<int>("computeEigenmodes");

            int numberOfComponents = filterTimeDerivative ? 2 : numberOfFrequencies;

            int numShowComponents= numberOfComponents;

            if( !filterTimeDerivative )
                numShowComponents += numberOfFrequencies;       // save forcings

            if( helmholtzFromDirectSolverWasComputed ) 
                    numShowComponents += numberOfComponents;     // save solution from direct Helmholtz 

            int saveErrors=0;
            if( computeErrors || helmholtzFromDirectSolverWasComputed ) 
            {
                saveErrors=true;
                numShowComponents += numberOfComponents;       // save errors, either true errors or difference with direct Helmholtz solution 
            }

            if( solveForScatteredField )
            {
                numShowComponents += numberOfComponents;       // save total field
            }

            if( !waveHoltzSolutionWasComputed && helmholtzFromDirectSolverWasComputed )
            {
        // save save direct solve and forcing
                numShowComponents = 2*numberOfComponents; 
            }

            realCompositeGridFunction q(cg,all,all,all,numShowComponents);
            q.setName("q");
            int ishow=0;  // counts components in q as we fill them in

            if( waveHoltzSolutionWasComputed )
            {
                if( !filterTimeDerivative )
                {
                    realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
                    for( int freq=0; freq<numberOfFrequencies; freq++ )
                    {
                        q.setName(sPrintF("v%d",freq),ishow); 

                        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        {
                            OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                            OV_GET_SERIAL_ARRAY(Real,q[grid],qLocal);
                            getIndex(cg[grid].dimension(),I1,I2,I3);
                            bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3);
                            if( ok )
                              qLocal(I1,I2,I3,ishow) = vLocal(I1,I2,I3,freq);   // waveHoltz solution

              // q[grid](all,all,all,ishow) = v[grid](all,all,all,freq);   // waveHoltz solution
                        }

                        ishow++;
                    }
                }
                else
                {
          // ---- complex solution ---
                    int iur=ishow; q.setName("ur",iur); ishow++;
                    int iui=ishow; q.setName("ui",iui); ishow++;

                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        OV_GET_SERIAL_ARRAY(Real,uh[grid],uhLocal);
                        OV_GET_SERIAL_ARRAY(Real,q[grid],qLocal);
                        getIndex(cg[grid].dimension(),I1,I2,I3);
                        bool ok=ParallelUtility::getLocalArrayBounds(uh[grid],uhLocal,I1,I2,I3);
                        if( ok )
                        {
                            qLocal(I1,I2,I3,iur) = uhLocal(I1,I2,I3,0);   // waveHoltz solution
                            qLocal(I1,I2,I3,iui) = uhLocal(I1,I2,I3,1);   // waveHoltz solution
                        }
                    }

                }
                if( solveForScatteredField )
                {
          // -- plot TOTAL field = incident + scattered ---

                    const Real & amp   = cgWave.dbase.get<Real>("ampPlaneWave");
                    const Real & kx    = cgWave.dbase.get<Real>("kxPlaneWave");
                    const Real & ky    = cgWave.dbase.get<Real>("kyPlaneWave");
                    const Real & kz    = cgWave.dbase.get<Real>("kzPlaneWave");
                    const Real & phi   = cgWave.dbase.get<Real>("phiPlaneWave");
                    const Real & omega = cgWave.dbase.get<Real>("omegaPlaneWave");

                    int iur=ishow; q.setName("urTotal",iur); ishow++;
                    int iui=ishow; q.setName("uiTotal",iui); ishow++;

                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        MappedGrid & mg = cg[grid];

                        OV_GET_SERIAL_ARRAY(Real,uh[grid],uhLocal);
                        OV_GET_SERIAL_ARRAY(Real,q[grid],qLocal);
                        OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);

                        getIndex(cg[grid].dimension(),I1,I2,I3);
                        bool ok=ParallelUtility::getLocalArrayBounds(uh[grid],uhLocal,I1,I2,I3);
                        if( ok )
                        {
                            if( cg.numberOfDimensions()==2 )
                            {
                                qLocal(I1,I2,I3,iur) = uhLocal(I1,I2,I3,0) + amp*sin( kx*xLocal(I1,I2,I3,0) + ky*xLocal(I1,I2,I3,1) + phi );
                                qLocal(I1,I2,I3,iui) = uhLocal(I1,I2,I3,1) - amp*cos( kx*xLocal(I1,I2,I3,0) + ky*xLocal(I1,I2,I3,1) + phi );
                            }
                            else
                            {
                                qLocal(I1,I2,I3,iur) = uhLocal(I1,I2,I3,0) + amp*sin( kx*xLocal(I1,I2,I3,0) + ky*xLocal(I1,I2,I3,1) + kz*xLocal(I1,I2,I3,2) + phi );
                                qLocal(I1,I2,I3,iui) = uhLocal(I1,I2,I3,1) - amp*cos( kx*xLocal(I1,I2,I3,0) + ky*xLocal(I1,I2,I3,1) + kz*xLocal(I1,I2,I3,2) + phi );
                            }

                        }            
                    }

                }
            }

            if( helmholtzFromDirectSolverWasComputed )
            { 
        // -- save solution from direct Helmholtz solve --
                for( int freq=0; freq<numberOfComponents; freq++ )
                {
                    q.setName(sPrintF("uDHS%d",freq),ishow);

                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        OV_GET_SERIAL_ARRAY(Real,uHelmholtz[grid],uHelmholtzLocal);
                        OV_GET_SERIAL_ARRAY(Real,q[grid],qLocal);
                        getIndex(cg[grid].dimension(),I1,I2,I3);
                        bool ok=ParallelUtility::getLocalArrayBounds(q[grid],qLocal,I1,I2,I3);
                        if( ok )
                            qLocal(I1,I2,I3,ishow) = uHelmholtzLocal(I1,I2,I3,freq);  // direct solver solution

                    }
                    ishow++; 
                }

            }

            if( waveHoltzSolutionWasComputed && ( saveErrors || helmholtzFromDirectSolverWasComputed )  )
            { 
                realCompositeGridFunction & error = cgWave.dbase.get<realCompositeGridFunction>("error");

                if( false && !saveErrors && helmholtzFromDirectSolverWasComputed )
                {
          // ** IS THIS NEEDED OR IS THE ERROR COMPUTED ALREADY ???   +++ TURN OFF FOR NOW : July 16, 2023

          // compute difference between WaveHoltz and Direct Helmholtz solve
                    if( !filterTimeDerivative )
                    {
                        realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
                        error = uHelmholtz - v;
                    }
                    else
                    {
                        error = uHelmholtz - uh;
                    }
                }

                for( int freq=0; freq<numberOfComponents; freq++ )
                {        
                    q.setName(sPrintF("err%d",freq),ishow);

                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        OV_GET_SERIAL_ARRAY(Real,error[grid],errorLocal);
                        OV_GET_SERIAL_ARRAY(Real,q[grid],qLocal);
                        getIndex(cg[grid].dimension(),I1,I2,I3);
                        bool ok=ParallelUtility::getLocalArrayBounds(q[grid],qLocal,I1,I2,I3);
                        if( ok )
                            qLocal(I1,I2,I3,ishow) = errorLocal(I1,I2,I3,freq);  // error or diff with direct

            // q[grid](all,all,all,ishow) = error[grid](all,all,all,freq);  // error or diff with direct
                    }
                    ishow++;
                }
            }

      // --- Save the forcings: ---
            if( !filterTimeDerivative )
            {
                realCompositeGridFunction & f = cgWave.dbase.get<realCompositeGridFunction>("f");
                for( int freq=0; freq<numberOfFrequencies; freq++ )
                {        
                    q.setName(sPrintF("f%d",freq),ishow);
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                    {
                        OV_GET_SERIAL_ARRAY(Real,f[grid],fLocal);
                        OV_GET_SERIAL_ARRAY(Real,q[grid],qLocal);
                        getIndex(cg[grid].dimension(),I1,I2,I3);
                        bool ok=ParallelUtility::getLocalArrayBounds(q[grid],qLocal,I1,I2,I3);
                        if( ok )
                            qLocal(I1,I2,I3,ishow) = fLocal(I1,I2,I3,freq);  // forcing

            // q[grid](all,all,all,ishow) = f[grid](all,all,all,freq);  // forcing
                    }
                    ishow++;
                }
            }

            if( ishow!=numShowComponents )
            {
                printF("cgwf:ERROR: ishow=%d != numShowComponents=%d\n",ishow,numShowComponents);
            }

            if( adjustPlotsForSuperGrid && useSuperGrid )
            {
                cgWave.adjustSolutionForSuperGrid( q );
            }

            show.saveSolution( q );              // save to show file
            show.endFrame();

      // show.close();
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
        else if( answer=="plot residual" )
        {
            realCompositeGridFunction & res = cgWaveHoltz.dbase.get<realCompositeGridFunction>("residual");

              if( false )
              {
                  res.display("res","%9.2e ");

              }
      // interpolate residual to make plots nicer: 
            if( FALSE )
            { 
                printF("\n XXXXXXXXXXX  INFO: interpolate the residual  XXXXXXXXXX  \n");
                res.interpolate();
            }

            ps.erase();
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
            psp.set(GI_TOP_LABEL,sPrintF("residual O%d omega=%.5g",orderOfAccuracy,omega));
            PlotIt::contour(ps,res,psp);
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

        }
        else if( answer=="plot errors" )
        {
            const int & computeEigenmodes     = cgWave.dbase.get<int>("computeEigenmodes");
            realCompositeGridFunction & error = cgWave.dbase.get<realCompositeGridFunction>("error");
            bool plotErrors=true;
            if( computeEigenmodes )
            {
                RealArray & eigenValues = cgWave.dbase.get<RealArray>("eigenValues");
                Real lambda = eigenValues(0);
                psp.set(GI_TOP_LABEL,sPrintF("Eigenvector Error O%d lambda=%.5g",orderOfAccuracy,lambda));

            }
            else if( computeErrors || helmholtzFromDirectSolverWasComputed )
            {
        // realCompositeGridFunction & error = cgWave.dbase.get<realCompositeGridFunction>("error");
                if( !computeErrors && helmholtzFromDirectSolverWasComputed )
                {
          // plot difference between WaveHoltz and Direct Helmholtz solve
                    realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");

          // This error should have been computed already I think: June 18,
          // error = uHelmholtz - v;          

                    if( numberOfFrequencies==1 )
                        psp.set(GI_TOP_LABEL,sPrintF("Error (disc. sol) O%d omega=%.5g",orderOfAccuracy,omega));
                    else
                        psp.set(GI_TOP_LABEL,sPrintF("Error (disc. sol) O%d numFreq=%d",orderOfAccuracy,numberOfFrequencies));
                }
                else
                {
                    if( numberOfFrequencies==1 )
                        psp.set(GI_TOP_LABEL,sPrintF("Error (true soln) O%d omega=%.5g",orderOfAccuracy,omega));
                    else
                        psp.set(GI_TOP_LABEL,sPrintF("Error (true soln) O%d, numFreq=%d",orderOfAccuracy,numberOfFrequencies));

                }

            }
            else
            {
                plotErrors=false;
                printF("plot errors: WARNING: this solution has no errors to compute. Probably not a known solution.\n"); 
                printF("           : Compute the direct Helmholtz solution in order to plot the difference between it and the WaveHoltz solution.\n"); 
            } 

            if( plotErrors )  
            {
                ps.erase();
                psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
                PlotIt::contour(ps,error,psp);
                psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true); 
            }        
        }
        else if( answer=="plot forcing" )
        {
            CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
            realCompositeGridFunction & f = cgWave.dbase.get<realCompositeGridFunction>("f");
            for( int freq=0; freq<numberOfFrequencies; freq++ )    
                f.setName(sPrintF("f%d",freq),freq);  
            ps.erase();
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
            psp.set(GI_TOP_LABEL,sPrintF("forcing"));
            PlotIt::contour(ps,f,psp);
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

        }
        else if( answer== "plot filter" )
        {
      // plot the WaveHoltz filter function beta and any eigenvalues use for deflation
            if( numToDeflate>0 && cgWave.dbase.has_key("uev") )
            {


                IntegerArray & eigNumbersToDeflate = cgWave.dbase.get<IntegerArray>("eigNumbersToDeflate");

                RealArray & eigTrue                = cgWave.dbase.get<RealArray>("eig"); // known true eigenvalues 

                RealArray eigenValues(numToDeflate);
                for( int ie=0; ie<numToDeflate; ie++ )
                    eigenValues(ie) = eigTrue(0,eigNumbersToDeflate(ie));

        // ::display(eigenValues,"eigs to deflate for plotFilter");
        // plotFilter( eigenValues, ps,psp );
                cgWave.plotFilter( eigenValues );

            }
            else
            {
                printF("There are no eigenvalues to deflate to plot with the filter function\n");
            }

        }
        else if( answer=="plot sequences" )
        {
      // Plot residual history
            const int & computeEigenmodes     = cgWave.dbase.get<int>("computeEigenmodes");
            const int & numberOfIterations    = cgWaveHoltz.dbase.get<int>("numberOfIterations");  // holds actual number of iterations taken
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
            realArray t(nd),x(nd,numberOfComponents);
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
            PlotIt::plot(ps,t,x,title,tName,xName,psp);
            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

            replot=true;

      // static void
      // plot(GenericGraphicsInterface &gi, 
      //      const realArray & t, 
      //      const realArray & x, 
      //      const aString & title = nullString, 
      //      const aString & tName       = nullString,
      //      const aString *xName        = NULL,
      //      GraphicsParameters & parameters=Overture::defaultGraphicsParameters()  );



      // CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
      // realCompositeGridFunction & f = cgWave.dbase.get<realCompositeGridFunction>("f");
      // for( int freq=0; freq<numberOfFrequencies; freq++ )    
      //   f.setName(sPrintF("f%d",freq),freq);  
      // ps.erase();
      // psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
      // psp.set(GI_TOP_LABEL,sPrintF("forcing"));
      // PlotIt::contour(ps,f,psp);
      // psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

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

            const Real Tperiod=numPeriods*twoPi/omega;  
            cgWave.dbase.get<real>("omega")     = omega;        // ** FIX ME **
            cgWave.dbase.get<real>("tFinal")    = Tperiod;      // ** FIX ME **
            cgWave.dbase.get<real>("Tperiod")   = Tperiod;      // ** FIX ME **
            cgWave.dbase.get<int>("numPeriods") = numPeriods;          // ** FIX ME **
            cgWave.dbase.get<int>("adjustOmega")= adjustOmega;  // 1 : choose omega from the symbol of D+t D-t 

      // const int & numPeriods            = cgWave.dbase.get<int>("numPeriods");
      // cgWave.dbase.get<real>("tFinal")  = numPeriods*twoPi/omega;

            int it=0;
            cgWave.advance( it );

        }
        else if( answer=="animate" )    
        {
            printF("Animate the Helmholtz (WaveHoltz) solution over time...\n");
            assert( filterTimeDerivative==0 );

            const int numberOfComponents = filterTimeDerivative ? 2 : numberOfFrequencies;
            realCompositeGridFunction q(cg,all,all,all,numberOfComponents);
            q.setName("v",0);

            realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");

            const int numPeriodsToPlot=2;
            Real omega = frequencyArray(0);
            Real tf = periodArray(0)*numPeriodsToPlot;
            Real dt = 1./(5.*omega);
            int numTimeSteps = max(2,ceil(tf/dt));
            dt = tf/numTimeSteps;
            printF("omega=%g, tf=%g, numTimeSteps=%d\n",omega,tf,numTimeSteps);

      // --- NEED TO SET MIN AND MAX ---
            
            for( int n=0; n<numTimeSteps; n++ )
            {
                Real t = n*dt;
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cg[grid];
                    OV_GET_SERIAL_ARRAY(Real,q[grid],qLocal);
                    OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);

                    qLocal(I1,I2,I3,0)= vLocal(I1,I3,I3)*cos(omega*t);
                }

                if( numberOfFrequencies==1 )
                    psp.set(GI_TOP_LABEL,sPrintF("CgWaveHoltz: FD%i%i%s omega=%.5g, t=%9.2e",orderOfAccuracy,orderOfAccuracyInTime,(upwind ? "u" : ""),omega,t));
                else
                    psp.set(GI_TOP_LABEL,sPrintF("CgWaveHoltz: FD%i%i%s numFreq=%d, t=%9.2e",orderOfAccuracy,orderOfAccuracyInTime,(upwind ? "u" : ""),numberOfFrequencies,t));

                if( cg.numberOfDimensions()==2 )
                    psp.set(GI_PLOT_CONTOUR_LINES,false);
              
                printF("plot step=%d, time t=%g\n",n,t);
                ps.erase();
                psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false); // wait inside contour
                PlotIt::contour(ps,q,psp);
                ps.redraw();

            } 


        }

        else
        {
            printF("cgwh: Unknown command = [%s]\n",(const char*)answer);
            ps.stopReadingCommandFile();
              
        }

        if( reComputeErrors )
        {
            reComputeErrors=false;

            if( computeErrors )
            {
                realCompositeGridFunction & v = cgWave.dbase.get<realCompositeGridFunction>("v");
                Real t=0; // compute errors with t=0 in cos(omega*t)

                Real maxErr = cgWave.getErrors( v, t );

                printF("cgwh: max-err =%8.2e (between WaveHoltz and known solution, all frequencies)\n",maxErr);

        // if( saveCheckFile )
        // {
        //   Real solutionNorm = maxNorm( v ); 
        //   cgWaveHoltz.saveCheckFile( checkFileCounter,maxErr,solutionNorm ); 
        //   checkFileCounter++; 
        // }

            }

            saveCheckFile=false;

        }

        if( answer=="contour" || (replot && plotChoices & 2) )
        {
            CgWave & cgWave = *cgWaveHoltz.dbase.get<CgWave*>("cgWave");
      // printF("plot contours...\n");

      // realCompositeGridFunction & q = filterTimeDerivative ? uh : cgWave.dbase.get<realCompositeGridFunction>("v");

            const int numberOfComponents = filterTimeDerivative ? 2 : numberOfFrequencies;
            realCompositeGridFunction q(cg,all,all,all,numberOfComponents);
            if( filterTimeDerivative )
            {
                q = uh;
                q.setName("ur",0);
                q.setName("ui",1);
            }
            else
            {
                q = cgWave.dbase.get<realCompositeGridFunction>("v");
                for( int freq=0; freq<numberOfFrequencies; freq++ )
                    q.setName(sPrintF("v%d",freq),freq);
            }

            if( solveForScatteredField && plotTotalField )
            {
                const Real & amp   = cgWave.dbase.get<Real>("ampPlaneWave");
                const Real & kx    = cgWave.dbase.get<Real>("kxPlaneWave");
                const Real & ky    = cgWave.dbase.get<Real>("kyPlaneWave");
                const Real & kz    = cgWave.dbase.get<Real>("kzPlaneWave");
                const Real & phi   = cgWave.dbase.get<Real>("phiPlaneWave");
                const Real & omega = cgWave.dbase.get<Real>("omegaPlaneWave");

                printF(" Plot total field: add plane-wave: amp=%g, [kx,ky,kz]=[%g,%g,%g]*2*pi, phi=%g\n",amp,kx/twoPi,ky/twoPi,kz/twoPi,phi);

                Index I1,I2,I3;
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    MappedGrid & mg = cg[grid];
                    OV_GET_SERIAL_ARRAY(Real,q[grid],qLocal);

                    mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
                    OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);
                    getIndex(mg.dimension(),I1,I2,I3);

          // NOTE: 
          //   We have forced with sin( k.x - omega* t)
                    if( cg.numberOfDimensions()==2 )
                    {
                        qLocal(I1,I2,I3,0) += amp*sin( kx*xLocal(I1,I2,I3,0) + ky*xLocal(I1,I2,I3,1) + phi );
                        qLocal(I1,I2,I3,1) -= amp*cos( kx*xLocal(I1,I2,I3,0) + ky*xLocal(I1,I2,I3,1) + phi );
                    }
                    else
                    {
                        qLocal(I1,I2,I3,0) += amp*sin( kx*xLocal(I1,I2,I3,0) + ky*xLocal(I1,I2,I3,1) + kz*xLocal(I1,I2,I3,2) + phi );
                        qLocal(I1,I2,I3,1) -= amp*cos( kx*xLocal(I1,I2,I3,0) + ky*xLocal(I1,I2,I3,1) + kz*xLocal(I1,I2,I3,2) + phi );
                    }
                }
            }
            const int & useSuperGrid = cgWave.dbase.get<int>("useSuperGrid");
            if( adjustPlotsForSuperGrid && useSuperGrid )
            {
        // printF("adjustPlotsForSuperGrid : WARNING: This will change the computed solution !!\n");

                cgWave.adjustSolutionForSuperGrid( q );
            }

            ps.erase();
            if( answer=="contour" )
                psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false); // wait inside contour

      // const bool useUpwind = ad4>0.;
            if( numberOfFrequencies==1 )
                psp.set(GI_TOP_LABEL,sPrintF("CgWaveHoltz: FD%i%i%s omega=%.5g",orderOfAccuracy,orderOfAccuracyInTime,(upwind ? "u" : ""),omega));
            else
                psp.set(GI_TOP_LABEL,sPrintF("CgWaveHoltz: FD%i%i%s numFreq=%d",orderOfAccuracy,orderOfAccuracyInTime,(upwind ? "u" : ""),numberOfFrequencies));

            if( cg.numberOfDimensions()==2 )
                psp.set(GI_PLOT_CONTOUR_LINES,false);
          
            PlotIt::contour(ps,q,psp);

            psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,true);

            plotChoices |= 2; 
        }    

    }
    
    ps.popGUI();  // pop dialog

    if( showFileIsOpen )
        show.close();

    return 0;
}