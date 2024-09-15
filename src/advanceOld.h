// ================================================================================================
/// \brief Advance the solution  
///    ***OLD VERSION***
/// This is the main advance routine for cgWave.
/// \param it (input) : iteration number when solving a Helmholtz problem.
// ================================================================================================
int CgWave::
advanceOld( int it )
{
  printf("\n *********  ENTERING ADVANCE OLD *****************\n");

  real time0=getCPU();
  const int & myid = dbase.get<int>("myid");
  const int & np = dbase.get<int>("np");

  int & iteration = dbase.get<int>("iteration");
  iteration=it;

  const int & debug  = dbase.get<int>("debug");
  FILE *& debugFile  = dbase.get<FILE*>("debugFile");
  FILE *& pDebugFile = dbase.get<FILE*>("pDebugFile");

  GenericGraphicsInterface & ps = gi;
  PlotStuffParameters psp;

  Real & c                          = dbase.get<real>("c");
  Real & cfl                        = dbase.get<real>("cfl");
  Real & tFinal                     = dbase.get<real>("tFinal");
  Real & tPlot                      = dbase.get<real>("tPlot");
  Real & dtMax                      = dbase.get<Real>("dtMax"); 
  Real & damp                       = dbase.get<Real>("damp"); 
  Real & dampSave                   = dbase.get<Real>("dampSave"); 

  const int & upwind                = dbase.get<int>("upwind");
  const int & numUpwindCorrections  = dbase.get<int>("numUpwindCorrections");
  const int & implicitUpwind        = dbase.get<int>("implicitUpwind");

  const int & dissipationFrequency  = dbase.get<int>("dissipationFrequency");
  // preComputeUpwindUt : true=precompute Ut in upwind dissipation, 
  //                      false=compute Ut inline in Gauss-Seidel fashion  
  int & preComputeUpwindUt = dbase.get<int>("preComputeUpwindUt");


  real & dt                         = dbase.get<real>("dt");
  real & dtUsed                     = dbase.get<real>("dtUsed");  // dt actually used
  RealArray & gridCFL               = dbase.get<RealArray>("gridCFL");
  const int & orderOfAccuracy       = dbase.get<int>("orderOfAccuracy");
  const int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");
  int & interactiveMode             = dbase.get<int>("interactiveMode");
  int & computeErrors               = dbase.get<int>("computeErrors"); 

  int & plotOptions                 = dbase.get<int>("plotOptions");
  int & plotChoices                 = dbase.get<int>("plotChoices");
  int & plotFrequency               = dbase.get<int>("plotFrequency");

  const int & useSuperGrid          = dbase.get<int>("useSuperGrid");
  IntegerArray & superGrid          = dbase.get<IntegerArray>("superGrid");
  
  const int & solveHelmholtz              = dbase.get<int>("solveHelmholtz");
  const int & computeTimeIntegral         = dbase.get<int>("computeTimeIntegral");
  const int & adjustOmega                 = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 
  const int & adjustHelmholtzForUpwinding = dbase.get<int>("adjustHelmholtzForUpwinding");
  const int & deflateWaveHoltz            = dbase.get<int>("deflateWaveHoltz");
  const int & numToDeflate                = dbase.get<int>("numToDeflate");
  const aString & eigenVectorFile         = dbase.get<aString>("eigenVectorFile"); //  name of file holding eigs and eigenvectors for deflation
  const int & computeEigenmodes           = dbase.get<int>("computeEigenmodes");
  const int & computeEigenmodesWithSLEPc  = dbase.get<int>("computeEigenmodesWithSLEPc");

  real & omega                      = dbase.get<real>("omega");
  real & Tperiod                    = dbase.get<real>("Tperiod");
  int & numPeriods                  = dbase.get<int>("numPeriods");
  real & omegaSave                  = dbase.get<real>("omegaSave");
  real & TperiodSave                = dbase.get<real>("TperiodSave");  
  real & dtSave                     = dbase.get<real>("dtSave");

  const int & numberOfFrequencies   = dbase.get<int>("numberOfFrequencies");
  RealArray & frequencyArray        = dbase.get<RealArray>("frequencyArray");
  RealArray & frequencyArrayAdjusted= dbase.get<RealArray>("frequencyArrayAdjusted");
  RealArray & periodArray           = dbase.get<RealArray>("periodArray"); 
  RealArray & periodArrayAdjusted   = dbase.get<RealArray>("periodArrayAdjusted"); 
  IntegerArray & numPeriodsArray    = dbase.get<IntegerArray>("numPeriodsArray");
  RealArray & frequencyArraySave    = dbase.get<RealArray>("frequencyArraySave");
  RealArray & periodArraySave       = dbase.get<RealArray>("periodArraySave");     

  int & addForcing                  = dbase.get<int>("addForcing");
  ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");

  const int & useKnownSolutionForFirstStep  = dbase.get<int>("useKnownSolutionForFirstStep");
  const int & takeImplicitFirstStep         = dbase.get<int>("takeImplicitFirstStep");

  const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

  const ModifiedEquationApproachEnum & modifiedEquationApproach = dbase.get<ModifiedEquationApproachEnum>("modifiedEquationApproach");

  IntegerArray & gridIsImplicit     = dbase.get<IntegerArray>("gridIsImplicit");
  const RealArray & bImp            = dbase.get<RealArray>("bImp");
  const RealArray & cImp            = dbase.get<RealArray>("cImp");

  real & solutionNorm               = dbase.get<real>("solutionNorm");  // save solution norm here 

  int & current = dbase.get<int>("current"); // hold the current solution index

  realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
  CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");
  Interpolant *& pInterpolant = dbase.get<Interpolant*>("pInterpolant");
  assert( pInterpolant!=NULL );
  Interpolant & interpolant = *pInterpolant;

  realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f");
  
  BoundaryConditionParameters bcParams;
  
  BCTypes::BCNames boundaryCondition=BCTypes::evenSymmetry; 
  // BCTypes::BCNames boundaryCondition=BCTypes::dirichlet;

  Index I1,I2,I3;

  if( modifiedEquationApproach==hierarchicalModifiedEquation ||
      modifiedEquationApproach==stencilModifiedEquation )
  {
    // -- We compute and save the coefficients in the Laplacian for the hierachical approach --

    if( orderOfAccuracy!=orderOfAccuracyInTime )
    {
      printF("ERROR: modified equation hierachical scheme currently needs orderOfAccuracy=orderOfAccuracyInTime,\n");
      printF("   *but*  orderOfAccuracy=%d, orderOfAccuracyInTime=%d\n",orderOfAccuracy,orderOfAccuracyInTime);
      OV_ABORT("error");
    }
    
    if( !dbase.has_key("lapCoeff") )
    {
      RealArray *& lapCoeff = dbase.put<RealArray*>("lapCoeff");
      lapCoeff = new RealArray[cg.numberOfComponentGrids()];           
      // 2D: c200, c110, c020, c100, c010
      // 3D: c200, c020, c002, c110, c101, c011, c100, c010, c001
      const int numCoeff = cg.numberOfDimensions()==2 ? 5 : 9;    
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        lapCoeff[grid]=NULL;
        const bool isRectangular = cg[grid].isRectangular();
        if( !isRectangular )
        {
          getIndex(cg[grid].dimension(),I1,I2,I3);
          lapCoeff[grid].redim(I1,I2,I3,numCoeff);
          lapCoeff[grid]=-1.; // this means the coefficients have not yet been set
        }
      }
    }
    if( modifiedEquationApproach==stencilModifiedEquation )
    {
      if( !dbase.has_key("stencilCoeff") )
      {
        RealArray *& stencilCoeff = dbase.put<RealArray*>("stencilCoeff");
        stencilCoeff = new RealArray[cg.numberOfComponentGrids()];        

        int stencilWidth = orderOfAccuracy+1;
        int numStencilCoeff;
        if( cg.numberOfDimensions()==2 )
          numStencilCoeff= stencilWidth*stencilWidth;
        else
          numStencilCoeff= stencilWidth*stencilWidth*stencilWidth;
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
          stencilCoeff[grid]=NULL;
          const bool isRectangular = cg[grid].isRectangular();
          if( !isRectangular )
          {
            getIndex(cg[grid].dimension(),I1,I2,I3);
            stencilCoeff[grid].redim(numStencilCoeff,I1,I2,I3);
            stencilCoeff[grid]=-1.; // this means the coefficients have not yet been set
          }
        }        

      }
    }
  }
  RealArray *lapCoeff = modifiedEquationApproach==standardModifiedEquation ? NULL : dbase.get<RealArray*>("lapCoeff");
  RealArray *stencilCoeff = modifiedEquationApproach!=stencilModifiedEquation ? NULL : dbase.get<RealArray*>("stencilCoeff");

  real t=0.;
  // real dtSquared = SQR(dt);
  real cSquared=c*c;
    
  // if( false )
  // {
  //   Tperiod = numPeriods*twoPi/omega;
  //   tFinal=Tperiod;
  //   tPlot=Tperiod/10.;  // plot this often
  // }
  
  if( !(plotOptions==noPlotting || plotOptions==plotAndWait || plotOptions==plotNoWait) )
  {
    printF("ERROR: invalid plotOptions=%d \n",plotOptions);
    OV_ABORT("ERROR");
  }
  GUIState *& runTimeDialog = dbase.get<GUIState*>("runTimeDialog");
  if( plotOptions != noPlotting )
  {
     // ---- Build the run-time dialog if plotting is on ----
    if( runTimeDialog==NULL )
    {
      buildRunTimeDialog();
    }
    GUIState & dialog = *runTimeDialog;
    ps.pushGUI(dialog);
  }

  // estimated number of time steps: (for output)
  // int numberOfTimeSteps=int( tFinal/dt+.5);

  const int maxNumberOfTimeSteps=1e10; 
  
  // adjust dt to reach tFinal
  // dt = tFinal/numberOfTimeSteps;

  // --- choose dt so we reach the next time to plot ---
  Real nextTimeToPlot = min(tPlot,tFinal);

  // Currenly we need to make sure that dt is constant when solving Helmholtz (for the accuracy of the time integral)
  if( solveHelmholtz )
      nextTimeToPlot=tFinal;

  // For now we do not adjust dt when using implicit time-stepping:
  // ** TURN OFF ADJUST TIME STEP *wdh* Nov 23, 2023
  bool adjustTimeStep = FALSE; // timeSteppingMethod==explicitTimeStepping;

  // NOTE: dtMax is the initial time-step determined by getTimeStep
  int numPlotSteps    =  ceil( nextTimeToPlot/dtMax );  // max number of steps 

  if( solveHelmholtz && timeSteppingMethod == implicitTimeStepping )
  { 
    // For implicit time-stepping we must take at least 5 steps per period when we adjust omega and dt
    // numPlotSteps = max(5*numPeriods,numPlotSteps);

    const int & minStepsPerPeriod = dbase.get<int>("minStepsPerPeriod");

    // Note: numPlotsSteps is per period
    if( numberOfFrequencies==1 )
    {
      numPlotSteps = max(minStepsPerPeriod,numPlotSteps);
    }
    else
    {
      // we need at least 5 steps for the smallest period
      const Real minFreq = min(frequencyArray);  // should be frequencyArray(0)
      const Real maxFreq = max(frequencyArray);
      const Real Tmax = twoPi/minFreq;  // longest period
      const Real Tmin = twoPi/maxFreq;  // shortest period
      int minSteps = ceil( minStepsPerPeriod*Tmax/Tmin ); // is this correct
      // printF("### Tmax=%12.4e, Tmin=%12.4e, Tmax/Tmin=%12.4e, minSteps=%d\n",Tmax,Tmin,Tmax/Tmin,minSteps);
      numPlotSteps = max(minSteps,numPlotSteps);
    }

    // printF("###(1) numPlotSteps=%d\n",numPlotSteps);
  }


  // --------------------------------------------------------------------------------
  // --- Adjust the time-step and frequencies/periods for WaveHoltz and EigenWave ---
  // --------------------------------------------------------------------------------
  adjustTimeStepAndFrequenciesMacro(); 


  int i1,i2,i3;

  // const bool takeForwardInitialStep=true; // false; // true; // *new* way : July 26, 2021

  int step=-1;
  t=0.;

  const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");    
  const int cur = 0; // current time level
  const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
  const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

  realCompositeGridFunction & u1 = u[cur];     // current time 
  realCompositeGridFunction & u2 = u[prev];    // previous time

  // u1=0.;
  // u2=0.;
  // u[next]=0.; // initialize so there are valid values in all ghost (for WaveHoltz PETSc solver)

  // -------------------------------------------------------------
  // -------------- Initial Conditions : u1 = u(t=0) -------------
  // -------------------------------------------------------------
  assignInitialSolutionMacro();

  // --------------------------------------------
  // --------------- FIRST STEP -----------------  
  // --------------------------------------------
  if( forcingOption==helmholtzForcing )  
  {
    takeFirstStep( cur, t ); 
  }
  else
  {
    if( !useKnownSolutionForFirstStep 
        // *wdh* Dec 23, 2023 && ( timeSteppingMethod==explicitTimeStepping || takeImplicitFirstStep )
        && orderOfAccuracy < 8 ) // do this for now 
    { 
      takeFirstStep( cur, t ); 
    }
    else
    { // Just get solution at t=dt from IC function
      if( debug>0 )
        printF("\n $$$$$$ CgWave:advance: Take First Step using KNOWN solution $$$$$$$$\n\n");   

      getInitialConditions( next,t+dt );
    }
  }


  if( computeTimeIntegral )
  {
    // When solving the Helmholtz problem with CgWaveHoltz we need to evaluate an integral 
    //      v  = (1/(2*T)* Int_0^T [  ( cos(omega*t)-.25)*u(x,t) dt ] 
    updateTimeIntegral( 0, firstStep, t, dt, u[cur] );
  }


  
  real cpua=getCPU();

  dtUsed = dt; // save dt actually used

  // =======================================================
  // ================== TIME STEPS =========================
  // =======================================================
  for( int i=0; i<=maxNumberOfTimeSteps; i++ )                    // take some time steps
  {
    step++;

    const int cur = (step +numberOfTimeLevelsStored) % numberOfTimeLevelsStored; // current time level
    const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
    const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

    current = cur; 

    plotAndOutputResults(); 

    // vOld : current guess for Helmholtz solution 
    realCompositeGridFunction & vOld = solveHelmholtz ?  dbase.get<realCompositeGridFunction>("vOld") : u[prev];

    if( step==0 && !takeImplicitFirstStep )
    {
      // --- first step has already been taken ----
       t+=dt;
       if( debug & 2 )
         printF("Skip first step since set to exact, or used time-periodic, or used Taylor series method.\n");
    }
    else
    { 

      // ----------------------------------------------
      // ---------------- take a time-step ------------
      // ----------------------------------------------
      for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
      {
        if( debug & 4 ) printF("Advance grid %i, step=%i, ...\n",grid,step);
            

        if( step==0 && takeImplicitFirstStep && !gridIsImplicit(grid) )
        {  // skip explicit grids, since these have been assigned in takeFirstStep
           continue;
        }

        MappedGrid & mg = cg[grid];
        
        // realArray & upg = up[grid];  // previous 
        // realArray & ucg = uc[grid];  // current 
        // realArray & ung = un[grid];  // next 
        // realArray & vg = v[grid];
            
  
        if( debug & 16 )
        {
          ::display(u[prev][grid],sPrintF("BEFORE: u=prev grid=%d t-dt=%9.3e",grid,t-dt),debugFile,"%11.3e ");
          ::display(u[cur ][grid],sPrintF("BEFORE: u=Cur  grid=%d t   =%9.3e",grid,t),debugFile,"%11.3e ");
        }
  
       const int useUpwindDissipation = upwind;  
       // const int useUpwindDissipation = ad4>0.;  // ** OLD **
  
       // useUpwindDissipation : true if grid is implicit AND upwind dissipation in ON AND implicitUpwind=1
       int useImplicitUpwindDissipation=false;
       if( gridIsImplicit(grid) && upwind && implicitUpwind )
         useImplicitUpwindDissipation=true;
  
  
        //     ----- UPWIND PREDICTOR CORRECTOR (UPC) ----
        // UPC scheme: 
        // Predict: 
        //   (P)  u^p = 2 u^n - u^{n-1} + dt^2{ c^2 Delta_h( u^n) + ... }      (normal update) 
        //        applyBC( u^p, t+dt )
        // Correct
        //   (C)  u^{n+1} = u^p + dt^2 * [ cupx*(D+xD-x)^q + cupy*(D+yD-y)^q ] ( u^p - u^{n-1} )/(2*dt)  ( add upwinding )
        //        applyBC( u^{n+1}, t+dt )
        // 
        // NOTE: To avoid an extra application of the BC's after (C) we instead delay stage (C) to the BEGINNING of the next step,
        // and actualy do step (C) FIRST on the current solution u^n:
        //   (C)  u^{n}  <-  u^{n} +  dt^2 * [ cupx*(D+xD-x)^q + cupy*(D+yD-y)^q ] ( u^n - u^{n-2} )/(2*dt)  ( add upwinding to u^n )
        //         Skip BC routine since doesn't appear to be necessary.
        //   (P)  u^{n+1} = 2 u^n - u^{n-1} + dt^2{ Stuff }                                            (normal update for u^{n+1}) 
        //         applyBC( u^{n+1}, t+dt )
        // 
        // Note: We need step>0 below since we need 3 levels of the solution to apply the upwind dissipation.
        // Note: Add upwinding at least 2 steps in a row 
        // bool startUpwinding = (takeForwardInitialStep && step>1) || (!takeForwardInitialStep && step>0); 
        bool startUpwinding = step>1;
        if( startUpwinding && useUpwindDissipation &&
               !useImplicitUpwindDissipation &&           // implicit upwind dissipation is added with option=0 below
               ( (step % dissipationFrequency)==0 ) ||
               ( (step % dissipationFrequency)==1 ) )
        {
          // ---- upwind dissipation ------
  
          const int prev2= (cur-2+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
  
          OV_GET_SERIAL_ARRAY(real,u[prev2][grid],umLocal);  // u(t-2*dt)
          OV_GET_SERIAL_ARRAY(real,u[prev ][grid],uLocal);   // u(t-dt)
          OV_GET_SERIAL_ARRAY(real,u[cur  ][grid],unLocal);  // u(t)

          OV_GET_SERIAL_ARRAY_CONDITIONAL(real,vOld[grid],vOldLocal,solveHelmholtz);
  
          if( debug & 4 )
          {
            printF("advance: ADD DISSIPATION: step=%d, t=%9.3e using [prev2,prev,cur]=[%d,%d,%d]\n",step,t,prev2,prev,cur);
          }

          int option = 1;    // 0=advance, 1=add dissipation
          callOptimizedAdvance();
  
          // --- THIS NEXT SECTION IS NOT CORRECT SINCE applyBC does **ALL** grids  **FIX ME**
          if( FALSE ) // ** TURN THIS OFF June 13, 2023  **CHECK ME**
          {
            if( TRUE )  // turned on for cic order 8 -- byt may just need to avoid extrapolating last ghost point ++++++++++++++++++++++ FIX
            {
              // testing: apply BC's after dissipation is added

              applyBoundaryConditions( u[cur],u[prev],t);  // 
            }   
            else
            {
              u[cur].periodicUpdate(); // *wdh* Aug 14, 2021 -- this may be needed 
            } 
          }
          else
          {
            u[cur][grid].periodicUpdate();  // this may be needed?
          }   
        }

  
        // -- optimized advance ---
        OV_GET_SERIAL_ARRAY(real,u[prev][grid],umLocal);
        OV_GET_SERIAL_ARRAY(real,u[cur ][grid],uLocal);
        OV_GET_SERIAL_ARRAY(real,u[next][grid],unLocal);
        OV_GET_SERIAL_ARRAY_CONDITIONAL(real,vOld[grid],vOldLocal,solveHelmholtz);
  
        if( debug & 4 )
        {
          printF("advance: UPDATE SOLUTION: step=%d, t=%9.3e using [prev,cur,nex]=[%d,%d,%d]\n",step,t,prev,cur,next);
        }      
  
        int option = 0;    // 0=advance, 1=add dissipation
  
        // =================== TAKE A STEP ==========================
        callOptimizedAdvance();
        
       
        if( timeSteppingMethod == explicitTimeStepping && debug & 16 )
        {
          ::display(u[next][grid],sPrintF("AFTER advOpt: uNext grid=%d t+dt=%9.3e",grid,t+dt),debugFile,"%11.3e ");
        }
  
        
      }  // end for grid
  
      if( timeSteppingMethod == implicitTimeStepping )
      {
        // boundary conditions for implicit time-stepping are added here:
        takeImplicitStep( t+dt );

        if( debug & 4 )
        {
          real maxErr = getErrors( u[next], t+dt );
          printF("\n ++++ AFTER takeImplicitStep (before explicit BCs): t+dt=%9.3e, maxErr=%9.2e +++\n\n",t+dt,maxErr);
        }

        if( debug & 16 )
        {
          for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            ::display(u[next][grid],sPrintF("AFTER takeImplicitStep (before explicit BCs): uNext grid=%d t+dt=%9.3e",grid,t+dt),debugFile,"%11.3e ");
        }

      }
  
      t+=dt;
  
      if( debug & 4 ) printf("...done advance, now interpolate...\n");
  
      // ----- apply boundary conditions ------
      if( true || (timeSteppingMethod != implicitTimeStepping) ) // ** TEST ********
      {
        if( debug & 1 && t <= 2.*dt )
          printF("Call applyBoundaryConditions : explicit BC's for step=%d, t(new)=%10.3e, dt=%10.3e\n",step,t,dt);

        bool applyExplicitBoundaryConditions=true;
        if( true )
          applyBoundaryConditions( u[next],u[cur],t, applyExplicitBoundaryConditions ); 
        else
          printF(" advance: do NOT Apply explicitBCs after implicit solve ############################################################ TEMP \n");
      }
  
      if( debug & 16 )
      {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
        {
          ::display(u[next][grid],sPrintF("after applyBC uNext grid=%d t=%9.3e",grid,t),debugFile,"%11.3e ");
        }
      }

     if( debug & 4 )
      printF("end time step=%d: t=%10.3e, dt=%10.3e\n",step,t,dt);      
  
    } // end take a time-step

    if( computeTimeIntegral )
    {
      // Integral in time -- trapezoidal rule
      // assert( i<=(numberOfTimeSteps-1) );
      
      // StepOptionEnum stepOption = i== (numberOfTimeSteps-1) ? lastStep : middleStep;
      StepOptionEnum stepOption = (t >tFinal-.5*dt) ? lastStep : middleStep;
      // printF(" computeTimeIntegral : t=%9.3e, stepOption=%d\n",t,stepOption);
      updateTimeIntegral( step+1, stepOption, t, dt, u[next] );
    }
    
  } // end for i .. number of steps

  timing(timeForAdvance) += getCPU()- cpua;

  current = next; // curent solution


  if( fabs(t-tFinal)/tFinal > dt*1e-3 )
  {
    printF("CgWave::advance:ERROR: done time-stepping but t is NOT equal to tFinal! Something is wrong.\n");
    printF("      t=%14.6e tFinal=%14.6e relative-diff=%9.3e\n",t,tFinal,(t-tFinal)/tFinal);
    OV_ABORT("FIX ME");
  }

  if( solveHelmholtz && adjustOmega )
  {
     // Reset adjusted omega and Tperiod
     omega   = omegaSave;
     Tperiod = TperiodSave;
     dt      = dtSave; 
     tFinal  = Tperiod; // reset too *wdh* Jul 25, 2021
     damp = dampSave;

    frequencyArray = frequencyArraySave;  // reset to original values 
    periodArray    = periodArraySave;     // reset to original values
    printF("RESET periodArray\n");   
    for( int freq=0; freq<numberOfFrequencies; freq++ )
    {
      printF("freq=%d: periodArray=%12.4e\n",freq,periodArray(freq));
    }   
    // OV_ABORT("stop here for now");  
  }
  if( solveHelmholtz )
  {
    // correct eigenfunction for Eigenwave
 
    if( !computeEigenmodesWithSLEPc )
    {
      correctEigenfunction();
      
      // adjust frequency for EigenWave 
      adjustEigenWaveFrequency();  
    }  
  }   



  if( false && solveHelmholtz )
  {
    // --- Zero out unused points for Kyrlov solver -- *wdh* Jan 5, 2020
    // The Krylov solver will measure the residual at all grid points
    // We st used points to zero so the residual will always be zero there
    // Note: for upwinding we use some unused points; thus it is important
    //  that we apply the boundary conditions to the initial conditions to
    //  restore these values.
    printF("cgWave:INFO zero unused points at finish.\n");
    realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");
    Range all; 
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
      MappedGrid & mg = cg[grid];
      OV_GET_SERIAL_ARRAY(real,u[current][grid],uLocal);
      OV_GET_SERIAL_ARRAY(real,v[grid],vLocal);
      OV_GET_SERIAL_ARRAY(int,mg.mask(),maskLocal);

      getIndex(mg.dimension(),I1,I2,I3); 
      FOR_3D(i1,i2,i3,I1,I2,I3)
      {
        if( maskLocal(i1,i2,i3) ==0 )
        {
          // uLocal(i1,i2,i3,all)=0.; // is this needed?
          // vLocal(i1,i2,i3,all)=0.; 
        }    
      }  
    }

  }

  // total number of steps taken:
  int & numberOfStepsTaken = dbase.get<int>("numberOfStepsTaken");
  numberOfStepsTaken += step;

  // number of steps this solve
  int & numberOfStepsPerSolve = dbase.get<int>("numberOfStepsPerSolve");
  numberOfStepsPerSolve = step;

  // -- pop run time dialog 
  if( plotOptions != noPlotting )
    ps.popGUI();

  delete runTimeDialog;

  runTimeDialog=NULL;
  
  timing(totalTime) += getCPU()-time0;
      
  saveSequencesToShowFile(); // June 10, 2020

  if( !solveHelmholtz )
    printStatistics();
  
  return 0;
}
