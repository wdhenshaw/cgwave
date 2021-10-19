// This file automatically generated from advance.bC with bpp.
// ====================== CgWave: advance  =====================

#include "CgWave.h"
#include "CompositeGridOperators.h";    
#include "PlotStuff.h"
#include "display.h"
#include "ParallelOverlappingGridInterpolator.h"
#include "ParallelUtility.h"
#include "LoadBalancer.h"
#include "gridFunctionNorms.h"
#include "OGPolyFunction.h"


#define FOR_3(i1,i2,i3,I1,I2,I3) for( i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( i1=I1.getBase(); i1<=I1.getBound(); i1++ )  

#define advWave EXTERN_C_NAME(advwave)
#define ogf EXTERN_C_NAME(ogf)
#define ogderiv EXTERN_C_NAME(ogderiv)

extern "C"
{

      /* Here are functions for TZ flow that can be called from fortran */

    real
    ogf(OGFunction *&ep, const real &x, const real &y,const real &z, const int & c, const real & t )
    {
        return (*ep)(x,y,z,c,t);
    }
    
    
    /* return a general derivative */
    void
    ogderiv(OGFunction *&ep, const int & ntd, const int & nxd, const int & nyd, const int & nzd, 
                      const real &x, const real &y, const real &z, const real & t, const int & n, real & ud )
    {
        ud=(*ep).gd(ntd,nxd,nyd,nzd,x,y,z,n,t);
    }

}

extern "C"
{
    /*  Optimized Advance Routine  */
    void advWave(const int&nd,
            const int&n1a,const int&n1b,const int&n2a,const int&n2b,const int&n3a,const int&n3b,
            const int&nd1a,const int&nd1b,const int&nd2a,const int&nd2b,const int&nd3a,const int&nd3b,
            const int&nd4a,const int&nd4b,
            const int&mask,const real&xy, const real&rx,    
            const real&um, const real&u, real&un, const real&f, const real&fa, const real& v,
            const int&bc, const int&ipar, const real&rpar, int&ierr );

}


// ==============================================================================================
// Macro: call the optimized advance routine 
// ==============================================================================================

// ==============================================================================================
// Macro: Plot solution, compute errors, output results
// ==============================================================================================


// ================================================================================================
/// \brief Advance the solution 
// ================================================================================================
int CgWave::
advance( int it )
{
    real time0=getCPU();
    const int & myid = dbase.get<int>("myid");
    const int & np = dbase.get<int>("np");

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

    const int & upwind                = dbase.get<int>("upwind");
  // real & ad4                        = dbase.get<real>("ad4"); // coeff of the artificial dissipation. **OLD**
    const int & dissipationFrequency  = dbase.get<int>("dissipationFrequency");
  // preComputeUpwindUt : true=precompute Ut in upwind dissipation, 
  //                      false=compute Ut inline in Gauss-Seidel fashion  
    const int & preComputeUpwindUt = dbase.get<int>("preComputeUpwindUt");


    real & dt                         = dbase.get<real>("dt");
    const int & orderOfAccuracy       = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");
    int & interactiveMode             = dbase.get<int>("interactiveMode");
    int & computeErrors               = dbase.get<int>("computeErrors"); 

    int & plotOptions                 = dbase.get<int>("plotOptions");
    int & plotChoices                 = dbase.get<int>("plotChoices");
    int & plotFrequency               = dbase.get<int>("plotFrequency");
    
    const int & solveHelmholtz        = dbase.get<int>("solveHelmholtz");
    const int & computeTimeIntegral   = dbase.get<int>("computeTimeIntegral");
    const int & adjustOmega           = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

    real & omega                      = dbase.get<real>("omega");
    real & Tperiod                    = dbase.get<real>("Tperiod");
    int & numPeriods                  = dbase.get<int>("numPeriods");
    real & omegaSave                  = dbase.get<real>("omegaSave");
    real & TperiodSave                = dbase.get<real>("TperiodSave");  
    real & dtSave                     = dbase.get<real>("dtSave");  

    int & addForcing                  = dbase.get<int>("addForcing");
    ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");

    const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

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

    
    Index I1,I2,I3;

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
    bool adjustTimeStep = timeSteppingMethod==explicitTimeStepping;

  // NOTE: dtMax is the initial time-step determined by getTimeStep
    int numPlotSteps    =  ceil( nextTimeToPlot/dtMax );  // max number of steps 
    if( true || adjustTimeStep )
    {
        dt = nextTimeToPlot/numPlotSteps;  // we adjust dt first time only 

        if( solveHelmholtz && adjustOmega )
        {
       // ---- adjust the scheme or omega to remove time-discretization errors -----

              Real dts=dt, omegas=omega; 
              if( timeSteppingMethod == implicitTimeStepping )
              {
         // --- adjust omega for implicit time-stepping ---

         // There are at least two ways to adjust the scheme (see notes)----

         // IMP-FIX 2 : adjust omega
                  dts = (1./omega)*sqrt( 2.*( 1./cos(2.*Pi/numPlotSteps) -1. ) );
                  omegas= (2.*Pi/numPlotSteps)*1/dts;

              }
              else
              {
         // --- adjust omega for explicit time stepping ---
         // Adjust omega for Helmholtz problems so that
         //    4 sin^2( omegas*dt/2 )/dt^2 = omega^2
         //    dt = Ts/N
         //    Ts = 2*pi/omegas       : adjusted period
         // where  N = number of steps we want to take

         // dt = Ts/N = 2*pi/( omegas*N ) -> omegas*dt/2 = pi/N
         // sin(pi/N) = omega*dt/2  
         //  -> dt = sin(pi/N)*(2/omega)
          //    omegas = (pi/N)*(2/dt)
                  dts = sin(Pi/numPlotSteps) * (2./omega);
                  omegas = (Pi/numPlotSteps) * (2./dts);
              }
              printF("\n ##### CgWave:adjust omega and dt for Helmholtz: dts=%18.10e (dt=%18.10e) omegas=%14.6e (omega=%14.6e) ####\n\n",dts,dt,omegas,omega);

              omegaSave   = omega;    // save original omega ( omega is reset below)
              TperiodSave = Tperiod;  // save original Tperiod ( Tperiod is reset below)
              dtSave      = dt;

              dt      = dts;
              omega   = omegas;
              Tperiod = twoPi/omegas;
              tFinal  = Tperiod; 
        }
    }

    if( debug & 2 )
        printF("CgWave:advance: nextTimeToPlot=%9.3e, numPlotSteps=%d, new dt=%9.3e (dtMax=%9.3e)\n",
                      nextTimeToPlot,numPlotSteps, dt,dtMax);

  // int plotSteps = (int)max(1.,tPlot/dt+.5);
  // plotSteps=10;
    
  //  const int showSteps = (int)max(1.,tShow/dt+.5);
  // Put these in class:
  // real timeForLaplace=0, timeForBoundaryConditions=0., timeForUpdateGhostBoundaries=0.,
  //   timeForInterpolate=0., timeForAdvance=0., timeForGetLocalArray=0.,
  //   timeForFinishBoundaryConditions=0.;
            
    printF("CgWave::advance: c=%g, omega=%g, Tperiod=%g, numPeriods=%d, tFinal=%g, plotOptions=%d\n",
                  c,omega,Tperiod,numPeriods,tFinal,plotOptions);


    int i1,i2,i3;

    const bool takeForwardInitialStep=true; // false; // true; // *new* way : July 26, 2021

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

  // ----- Initial Conditions -----
    if( !solveHelmholtz )
    {
        getInitialConditions( cur,t );
    }
    else
    {
    // WaveHoltz: initial condition is provided in v: 
        realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");
        u1=v;

        applyBoundaryConditions( u1, t );  // *wdh* Mar 6, 2020 -- try this 
    }

    if( computeTimeIntegral )
    {
    // When solving the Helmholtz problem with CgWaveHoltz we need to evaluate an integral 
    //      v  = (1/(2*T)* Int_0^T [  ( cos(omega*t)-.25)*u(x,t) dt ] 
        updateTimeIntegral( firstStep, t, u[cur] );
    }



  // ----- FIRST STEP ------  

    if( takeForwardInitialStep )
    {
        if( forcingOption==helmholtzForcing )
            takeFirstStep( cur, t ); // -- NOTE: for now initial step is only called for Helmholtz problems *fix me*
        else
            getInitialConditions( next,t+dt );
    }
    else
    {
   // test backward step 
        if( forcingOption==helmholtzForcing )
            takeFirstBackwardStep( cur, t );
        else
            getInitialConditions( prev,t-dt );
    }

    

    
    real cpua=getCPU();

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

      // if( (i % plotSteps) == 0 || (t >tFinal-.5*dt) )  // plot solution every 'plotSteps' steps
            if( i==0 || (t >nextTimeToPlot-.5*dt) || (i % plotFrequency)==0 )  // plot solution every 'plotSteps' steps
            {
         // ---- plot first step and last step and every plotSteps ---
                if( debug & 4 )
                    printF("it=%d: completed step %i, t=%8.2e, nextTimeToPlot=%9.3e\n",it,i,t,nextTimeToPlot);
        // getErrors( u2, t-dt );
                real maxErr = getErrors( u[cur], t );
        // ** fix this: 
        // estimated number of time steps: (for output)
                int numberOfTimeSteps=int( tFinal/dt+.5);
        // -- output with a nice format ---
        // sPrintF(buff,"%%4i: %%%is   ([%%2i:%%5i],[%%2i:%%5i],[%%2i:%%5i])  %%12g   %%8.2e %%8.2e    %%s\n",maxNameLength);
                const Real cpuTime = getCPU()- cpua;
                const int numDigits = ceil( log10(numberOfTimeSteps) );
                char myFormat[180];
                if( computeErrors )
                {
                    sPrintF(myFormat,"cgWave:FD%i%i t=%%9.3e (%%%ii steps) dt=%%9.3e maxErr=%%9.2e, ||u||=%%9.2e, cpu=%%9.2e(s)\n",orderOfAccuracyInTime,orderOfAccuracy,numDigits);
                    printF(myFormat,t,step,dt,maxErr,solutionNorm,cpuTime);
                }
                else
                {
                    sPrintF(myFormat,"cgWave:FD%i%i t=%%9.3e (%%%ii steps) dt=%%9.3e ||u||=%%9.2e, cpu=%%9.2e(s)\n",orderOfAccuracyInTime,orderOfAccuracy,numDigits);
                    printF(myFormat,t,step,dt,solutionNorm,cpuTime);
                }
        // output results (e.g. print errors to the check file)
                outputResults( cur, t );
        // Optionally save results to a show file
                saveShow( current,t,dt );
        // int current = step % 2;      // ********* FIX ME
        // ------ plot the solution -----
                int finished = plot(current, t, dt );
                if( finished ) break;
                if( t >tFinal-.5*dt )
                {
          // we are done (unless tFinal is increased in the next call to plot). plot solution at final time
                    if( interactiveMode==1 && plotOptions!=noPlotting )
                        plotOptions=plotNoWait; 
                    plot(current, t, dt );
                    if( t >tFinal-.5*dt ) // tFinal may have been increased, so check again
                    { 
                        finished=true;
                        break;
                    }
                }
                if( t >nextTimeToPlot-.5*dt )
                {
           // --- Recompute nextTimeToPlot and dt ----
                    const Real prevTimeToPlot=nextTimeToPlot; 
                    nextTimeToPlot    = min(nextTimeToPlot+tPlot,tFinal);
                    Real timeInterval = nextTimeToPlot-prevTimeToPlot;
                    int numPlotSteps  = ceil( timeInterval/dtMax ); 
                    if( adjustTimeStep )
                        dt = timeInterval/numPlotSteps;
                    if( debug & 4 ) 
                    {
                        printF("CgWave:advance: nextTimeToPlot=%9.3e, numPlotSteps=%d, new dt=%9.3e (dtMax=%9.3e) ratio=%5.3f\n",
                                nextTimeToPlot,numPlotSteps, dt,dtMax,dt/dtMax);
                    }
                }
            }

        if( step==0 && takeForwardInitialStep )
        {
      // --- first step has already been taken ----
              t+=dt;
              printF("Skip first step since set to exact\n");
        }
        else
        { // ------ take a time-step -----
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
                if( debug & 4 ) printf("Advance grid %i, step=%i...\n",grid,step);
                        
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
    
       // useUpwindDissipation : true if grid is implicit AND upwind dissipation in ON 
              int useImplicitUpwindDissipation=false;
              if( gridIsImplicit(grid) && upwind )
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
                bool startUpwinding = (takeForwardInitialStep && step>1) || (!takeForwardInitialStep && step>0); 
                if( startUpwinding && useUpwindDissipation &&
               !useImplicitUpwindDissipation &&           // implicit upwind dissipation is added added with option=0 below
                              ( (step % dissipationFrequency)==0 ) ||
                              ( (step % dissipationFrequency)==1 ) )
                {
          // ---- upwind dissipation ------
    
                    const int prev2= (cur-2+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
    
                    OV_GET_SERIAL_ARRAY(real,u[prev2][grid],umLocal);  // u(t-2*dt)
                    OV_GET_SERIAL_ARRAY(real,u[prev ][grid],uLocal);   // u(t-dt)
                    OV_GET_SERIAL_ARRAY(real,u[cur  ][grid],unLocal);  // u(t)
    
                    if( debug & 4 )
                    {
                        printF("advance: ADD DISSIPATION: step=%d, t=%9.3e using [prev2,prev,cur]=[%d,%d,%d]\n",step,t,prev2,prev,cur);
                    }
                    int option = 1;    // 0=advance, 1=add dissipation
                        real cpuOpt=getCPU();
                        const int orderOfAccuracyInSpace=orderOfAccuracy;
                        const int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");
            // ** FIX ME **
                        int numberOfForcingFunctions=1;
                        int fCurrent=0;
                        real sosupParameter=1.;
                        OGFunction *tz = dbase.get<OGFunction*>("tz");
                        const bool isRectangular = mg.isRectangular();
                        int gridType = isRectangular? 0 : 1;
                        int ipar[]={option,                        // ipar[ 0]
                                                grid,                          // ipar[ 1]
                                                gridType,                      // ipar[ 2]
                                                orderOfAccuracyInSpace,        // ipar[ 3]
                                                orderOfAccuracyInTime,         // ipar[ 4]
                                                (int)addForcing,               // ipar[ 5]
                                                (int)forcingOption,            // ipar[ 6]
                                                numberOfForcingFunctions,      // ipar[ 7]
                                                fCurrent,                      // ipar[ 8] 
                                                debug,                         // ipar[ 9]
                                                gridIsImplicit(grid),          // ipar[10]
                                                useUpwindDissipation,          // ipar[11]
                                                useImplicitUpwindDissipation,  // ipar[12]
                                                preComputeUpwindUt             // ipar[13]
                                                                        };  //
                        real dx[3]={1.,1.,1.};
                        if( isRectangular )
                            mg.getDeltaX(dx);
                        real rpar[20];
                        rpar[ 0]=c;
                        rpar[ 1]=dt;
                        rpar[ 2]=dx[0];
                        rpar[ 3]=dx[1];
                        rpar[ 4]=dx[2];
                        rpar[ 5]=mg.gridSpacing(0);
                        rpar[ 6]=mg.gridSpacing(1);
                        rpar[ 7]=mg.gridSpacing(2);
                        rpar[ 8]=t;
                        rpar[ 9]= (real &)tz;  // twilight zone pointer
                        rpar[10]=sosupParameter;
                        rpar[11]=dbase.get<real>("omega");
                        rpar[12]=bImp(0);
                        rpar[13]=bImp(1);
                        rpar[14]=bImp(2);
                        rpar[15]=bImp(3);
                        intArray & mask = mg.mask();
                        OV_GET_SERIAL_ARRAY(int,mask,maskLocal);
                        int *maskptr = maskLocal.getDataPointer(); 
                        getIndex(mg.gridIndexRange(),I1,I2,I3);
                        bool ok=ParallelUtility::getLocalArrayBounds(mask,maskLocal,I1,I2,I3);
                        if( ok )
                        {
                            real *umptr = umLocal.getDataPointer();
                            real *uptr  = uLocal.getDataPointer();
                            real *unptr = unLocal.getDataPointer();
                            real *vptr = umptr;  // This may be OK for dissipation stage
              // We may need v for upwinding or 4th-order curvilinear to hold Lap_2h(u) 
                            RealArray uDot;
                            if( ( option==1 && preComputeUpwindUt ) || 
                                    ( orderOfAccuracy==4 && !isRectangular) )
                            {
                                uDot.redim(uLocal);
                                vptr = uDot.getDataPointer();  // do this for now 
                            }
                            real *fptr  = umptr;  // use this if there is nor forcing 
                            if( forcingOption != noForcing )
                            {
                                OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);
                                fptr  = fLocal.getDataPointer();
                            }
                            real *faptr  = fptr;   // do this for now 
                            real *rxptr;
                            if( isRectangular )
                            {
                                rxptr=uptr;
                            }
                            else
                            {
                                OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal);
                                rxptr = rxLocal.getDataPointer();
                            }
                            real *xyptr = uptr;
                            const bool centerNeeded = forcingOption == twilightZoneForcing; // **** IS THIS TRUE ?
                            if( centerNeeded )
                            {
                // Pass the xy array for twilight-zone 
                                mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
                                OV_GET_SERIAL_ARRAY(real,mg.center(),xLocal);
                                xyptr = xLocal.getDataPointer();
                            }
                            int ierr;
                            advWave(mg.numberOfDimensions(),
                                            I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                            uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),
                                            uLocal.getBase(2),uLocal.getBound(2),
                                            uLocal.getBase(3),uLocal.getBound(3),
                                            *maskptr,*xyptr, *rxptr,  
                                            *umptr,*uptr,*unptr, *fptr,
                                            *faptr,  // forcing at multiple time levels 
                                            *vptr,   // to hold uDot 
                                            mg.boundaryCondition(0,0), ipar[0], rpar[0], ierr );
                            if( false && option==0 )
                            {
                                ::display(uDot,"v=Lap2h(u)","%8.2e ");
                            }
                            if( false && grid==1 && option==1)
                            {
                                ::display(uDot,"uDot","%8.2e ");
                                ::display(umLocal,"umLocal","%8.2e ");
                                ::display(uLocal,"uLocal","%8.2e ");
                                ::display(unLocal,"unLocal","%8.2e ");
                            }
                        }
                        real cpu1=getCPU();
                        if( option==0 )
                        {
                            if( isRectangular )
                                timing(timeForAdvanceRectangularGrids) += cpu1-cpuOpt;
                            else
                                timing(timeForAdvanceCurvilinearGrids) += cpu1-cpuOpt;
                        }
                        else
                        {
                            timing(timeForDissipation) += cpu1-cpuOpt;
                        }
    
                    if( FALSE )
                    {
            // testing: apply BC's after dissipation is added
                        applyBoundaryConditions( u[cur],t); 
                    }   
                    else
                    {
                        u[cur].periodicUpdate(); // *wdh* Aug 14, 2021 -- this may be needed (cgsm)
                    }    
                }
    
        // -- optimized advance ---
                OV_GET_SERIAL_ARRAY(real,u[prev][grid],umLocal);
                OV_GET_SERIAL_ARRAY(real,u[cur ][grid],uLocal);
                OV_GET_SERIAL_ARRAY(real,u[next][grid],unLocal);
    
                if( debug & 4 )
                {
                    printF("advance: UPDATE SOLUTION: step=%d, t=%9.3e using [prev,cur,nex]=[%d,%d,%d]\n",step,t,prev,cur,next);
                }      
    
                int option = 0;    // 0=advance, 1=add dissipation
    
        // =================== TAKE A STEP ==========================
                    real cpuOpt=getCPU();
                    const int orderOfAccuracyInSpace=orderOfAccuracy;
                    const int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");
          // ** FIX ME **
                    int numberOfForcingFunctions=1;
                    int fCurrent=0;
                    real sosupParameter=1.;
                    OGFunction *tz = dbase.get<OGFunction*>("tz");
                    const bool isRectangular = mg.isRectangular();
                    int gridType = isRectangular? 0 : 1;
                    int ipar[]={option,                        // ipar[ 0]
                                            grid,                          // ipar[ 1]
                                            gridType,                      // ipar[ 2]
                                            orderOfAccuracyInSpace,        // ipar[ 3]
                                            orderOfAccuracyInTime,         // ipar[ 4]
                                            (int)addForcing,               // ipar[ 5]
                                            (int)forcingOption,            // ipar[ 6]
                                            numberOfForcingFunctions,      // ipar[ 7]
                                            fCurrent,                      // ipar[ 8] 
                                            debug,                         // ipar[ 9]
                                            gridIsImplicit(grid),          // ipar[10]
                                            useUpwindDissipation,          // ipar[11]
                                            useImplicitUpwindDissipation,  // ipar[12]
                                            preComputeUpwindUt             // ipar[13]
                                                                    };  //
                    real dx[3]={1.,1.,1.};
                    if( isRectangular )
                        mg.getDeltaX(dx);
                    real rpar[20];
                    rpar[ 0]=c;
                    rpar[ 1]=dt;
                    rpar[ 2]=dx[0];
                    rpar[ 3]=dx[1];
                    rpar[ 4]=dx[2];
                    rpar[ 5]=mg.gridSpacing(0);
                    rpar[ 6]=mg.gridSpacing(1);
                    rpar[ 7]=mg.gridSpacing(2);
                    rpar[ 8]=t;
                    rpar[ 9]= (real &)tz;  // twilight zone pointer
                    rpar[10]=sosupParameter;
                    rpar[11]=dbase.get<real>("omega");
                    rpar[12]=bImp(0);
                    rpar[13]=bImp(1);
                    rpar[14]=bImp(2);
                    rpar[15]=bImp(3);
                    intArray & mask = mg.mask();
                    OV_GET_SERIAL_ARRAY(int,mask,maskLocal);
                    int *maskptr = maskLocal.getDataPointer(); 
                    getIndex(mg.gridIndexRange(),I1,I2,I3);
                    bool ok=ParallelUtility::getLocalArrayBounds(mask,maskLocal,I1,I2,I3);
                    if( ok )
                    {
                        real *umptr = umLocal.getDataPointer();
                        real *uptr  = uLocal.getDataPointer();
                        real *unptr = unLocal.getDataPointer();
                        real *vptr = umptr;  // This may be OK for dissipation stage
            // We may need v for upwinding or 4th-order curvilinear to hold Lap_2h(u) 
                        RealArray uDot;
                        if( ( option==1 && preComputeUpwindUt ) || 
                                ( orderOfAccuracy==4 && !isRectangular) )
                        {
                            uDot.redim(uLocal);
                            vptr = uDot.getDataPointer();  // do this for now 
                        }
                        real *fptr  = umptr;  // use this if there is nor forcing 
                        if( forcingOption != noForcing )
                        {
                            OV_GET_SERIAL_ARRAY(real,f[grid],fLocal);
                            fptr  = fLocal.getDataPointer();
                        }
                        real *faptr  = fptr;   // do this for now 
                        real *rxptr;
                        if( isRectangular )
                        {
                            rxptr=uptr;
                        }
                        else
                        {
                            OV_GET_SERIAL_ARRAY(real,mg.inverseVertexDerivative(),rxLocal);
                            rxptr = rxLocal.getDataPointer();
                        }
                        real *xyptr = uptr;
                        const bool centerNeeded = forcingOption == twilightZoneForcing; // **** IS THIS TRUE ?
                        if( centerNeeded )
                        {
              // Pass the xy array for twilight-zone 
                            mg.update(MappedGrid::THEvertex | MappedGrid::THEcenter );
                            OV_GET_SERIAL_ARRAY(real,mg.center(),xLocal);
                            xyptr = xLocal.getDataPointer();
                        }
                        int ierr;
                        advWave(mg.numberOfDimensions(),
                                        I1.getBase(),I1.getBound(),I2.getBase(),I2.getBound(),I3.getBase(),I3.getBound(),
                                        uLocal.getBase(0),uLocal.getBound(0),uLocal.getBase(1),uLocal.getBound(1),
                                        uLocal.getBase(2),uLocal.getBound(2),
                                        uLocal.getBase(3),uLocal.getBound(3),
                                        *maskptr,*xyptr, *rxptr,  
                                        *umptr,*uptr,*unptr, *fptr,
                                        *faptr,  // forcing at multiple time levels 
                                        *vptr,   // to hold uDot 
                                        mg.boundaryCondition(0,0), ipar[0], rpar[0], ierr );
                        if( false && option==0 )
                        {
                            ::display(uDot,"v=Lap2h(u)","%8.2e ");
                        }
                        if( false && grid==1 && option==1)
                        {
                            ::display(uDot,"uDot","%8.2e ");
                            ::display(umLocal,"umLocal","%8.2e ");
                            ::display(uLocal,"uLocal","%8.2e ");
                            ::display(unLocal,"unLocal","%8.2e ");
                        }
                    }
                    real cpu1=getCPU();
                    if( option==0 )
                    {
                        if( isRectangular )
                            timing(timeForAdvanceRectangularGrids) += cpu1-cpuOpt;
                        else
                            timing(timeForAdvanceCurvilinearGrids) += cpu1-cpuOpt;
                    }
                    else
                    {
                        timing(timeForDissipation) += cpu1-cpuOpt;
                    }
                
    
                if( timeSteppingMethod == explicitTimeStepping && debug & 16 )
                {
                    ::display(u[next][grid],sPrintF("AFTER advOpt: uNext grid=%d t+dt=%9.3e",grid,t+dt),debugFile,"%11.3e ");
                }
    
                
            }  // end for grid
    
            if( timeSteppingMethod == implicitTimeStepping )
            {
                takeImplicitStep( t+dt );

                if( debug & 16 )
                {
                    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                        ::display(u[next][grid],sPrintF("AFTER takeImplicitStep: uNext grid=%d t+dt=%9.3e",grid,t+dt),debugFile,"%11.3e ");
                }

            }
    
            t+=dt;
    
            if( debug & 4 ) printf("...done advance, now interpolate...\n");
    
      // ----- apply boundary conditions ------
            applyBoundaryConditions( u[next],t);
            
    
            if( debug & 16 )
            {
                for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                {
                    ::display(u[next][grid],sPrintF("after applyBC uNext grid=%d t=%9.3e",grid,t),debugFile,"%11.3e ");
                }
            }
    
        } // end take a time-step

        if( computeTimeIntegral )
        {
      // Integral in time -- trapezoidal rule
      // assert( i<=(numberOfTimeSteps-1) );
            
      // StepOptionEnum stepOption = i== (numberOfTimeSteps-1) ? lastStep : middleStep;
            StepOptionEnum stepOption = (t >tFinal-.5*dt) ? lastStep : middleStep;
      // printF(" computeTimeIntegral : t=%9.3e, stepOption=%d\n",t,stepOption);
            updateTimeIntegral( stepOption, t, u[next] );
        }
        
    } // end for i .. number of steps

    timing(timeForAdvance) += getCPU()- cpua;

    current = next; // curent solution


    if( fabs(t-tFinal) > REAL_EPSILON*tFinal*1000. )
    {
        printF("CgWave::advance:ERROR: done time-stepping but t is NOT equal to tFinal! Something is wrong.\n");
          printF("      t=%14.6e tFinal=%14.6e diff=%9.3e\n",t,tFinal,t-tFinal);
          OV_ABORT("FIX ME");
    }

    if( solveHelmholtz && adjustOmega )
    {
     // Reset adjusted omega and Tperiod
          omega   = omegaSave;
          Tperiod = TperiodSave;
          dt      = dtSave; 
          tFinal  = Tperiod; // reset too *wdh* Jul 25, 2021
    }   

    int & numberOfStepsTaken = dbase.get<int>("numberOfStepsTaken");
    numberOfStepsTaken += step;

  // -- pop run time dialog 
    if( plotOptions != noPlotting )
        ps.popGUI();

    delete runTimeDialog;

    runTimeDialog=NULL;
    
    timing(totalTime) += getCPU()-time0;
            
    if( !solveHelmholtz )
        printStatistics();
    
    return 0;
}



