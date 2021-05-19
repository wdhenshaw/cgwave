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

#define FOR_3D(i1,i2,i3,I1,I2,I3) for( int i3=I3.getBase(); i3<=I3.getBound(); i3++ )  for( int i2=I2.getBase(); i2<=I2.getBound(); i2++ )  for( int i1=I1.getBase(); i1<=I1.getBound(); i1++ )

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

// ================================================================================================
/// \brief Advance the solution 
// ================================================================================================
int CgWave::
advance( int it )
{
    real time0=getCPU();
    const int & myid = dbase.get<int>("myid");
    const int & np = dbase.get<int>("np");

    FILE *& debugFile = dbase.get<FILE*>("debugFile");

    GenericGraphicsInterface & ps = gi;
    PlotStuffParameters psp;

    Real & c                          = dbase.get<real>("c");
    Real & cfl                        = dbase.get<real>("cfl");
    Real & tFinal                     = dbase.get<real>("tFinal");
    Real & tPlot                      = dbase.get<real>("tPlot");
    Real & dtMax                      = dbase.get<Real>("dtMax"); 

    real & ad4                        = dbase.get<real>("ad4"); // coeff of the artificial dissipation.
    const int & dissipationFrequency  = dbase.get<int>("dissipationFrequency");

    real & dt                         = dbase.get<real>("dt");
    const int & orderOfAccuracy       = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInTime = dbase.get<int>("orderOfAccuracyInTime");
    int & debug                       = dbase.get<int>("debug");
    int & interactiveMode             = dbase.get<int>("interactiveMode");
    int & computeErrors               = dbase.get<int>("computeErrors"); 

    int & plotOptions                 = dbase.get<int>("plotOptions");
    int & plotChoices                 = dbase.get<int>("plotChoices");
    
    const int & solveHelmholtz        = dbase.get<int>("solveHelmholtz");
    const int & computeTimeIntegral   = dbase.get<int>("computeTimeIntegral");
    const int & adjustOmega           = dbase.get<int>("adjustOmega");  // 1 : choose omega from the symbol of D+t D-t 

    real & omega                      = dbase.get<real>("omega");
    real & Tperiod                    = dbase.get<real>("Tperiod");
    int & numPeriods                  = dbase.get<int>("numPeriods");
    real & omegaSave                  = dbase.get<real>("omegaSave");
    real & TperiodSave                = dbase.get<real>("TperiodSave");  

    int & addForcing                  = dbase.get<int>("addForcing");
    ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");

    const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

    IntegerArray & gridIsImplicit     = dbase.get<IntegerArray>("gridIsImplicit");
    RealArray & cImp                  = dbase.get<RealArray>("cImp");

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
       // Adjust omega for Helmholtz problems so that
       //    4 sin^2( omegas*dt/2 )/dt^2 = omega^2
       //    dt = Ts/N
       //    Ts = 2*pi/omegas       : adjusted period
       // where  N = number of steps we want to take

       // dt = Ts/N = 2*pi/( omegas*N ) -> omegas*dt/2 = pi/N
       // sin(pi/N) = omega*dt/2  
       //  -> dt = sin(pi/N)*(2/omega)
        //    omegas = (pi/N)*(2/dt)
              Real dts = sin(Pi/numPlotSteps) * (2./omega);
              Real omegas = (Pi/numPlotSteps) * (2./dts);
              printF("\n ##### CgWave:adjust omega and dt for Helmholtz: dts=%12.4e (dt=%12.4e) omegas=%12.4e (omega=%12.4e) ####\n\n",dts,dt,omegas,omega);

              omegaSave   = omega;    // save original omega ( omega is reset below)
              TperiodSave = Tperiod;  // save original Tperiod ( Tperiod is reset below)

              dt      = dts;
              omega   = omegas;
              Tperiod = twoPi/omegas;
              tFinal  = Tperiod; 
        }
    }


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


    int step=-1;
    t=0.;

    const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");    
    const int cur = 0; // current time level
    const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
    const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

    realCompositeGridFunction & u1 = u[cur];     // current time 
    realCompositeGridFunction & u2 = u[prev];    // previous time

    if( !solveHelmholtz )
    {
    // ----- Initial Conditions -----
        getInitialConditions( t );
        if( forcingOption==helmholtzForcing )
        { // test backward step 
            takeFirstBackwardStep( cur, t );
        }
        
    }
    else
    {
        realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");
        u1=v;

        applyBoundaryConditions( u1, t );  // *wdh* Mar 6, 2020 -- try this 
          
    // -- we assume v has been set ---
        takeFirstBackwardStep( cur, t );
    }
    
    if( computeTimeIntegral )
    {
    // When solving the Helmholtz problem with CgWaveHoltz we need to evaluate an integral 
    //      v  = (1/(2*T)* Int_0^T [  ( cos(omega*t)-.25)*u(x,t) dt ] 
        updateTimeIntegral( firstStep, t, u[cur] );
    }
    
    real cpua=getCPU();

  // ================== TIME STEPS =========================
    for( int i=0; i<=maxNumberOfTimeSteps; i++ )                    // take some time steps
    {
        step++;

        const int cur = (step +numberOfTimeLevelsStored) % numberOfTimeLevelsStored; // current time level
        const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
        const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

        current = cur; 

    // if( (i % plotSteps) == 0 || (t >tFinal-.5*dt) )  // plot solution every 'plotSteps' steps
        if( i==0 || (t >nextTimeToPlot-.5*dt) )  // plot solution every 'plotSteps' steps
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
                sPrintF(myFormat,"cgWave:FD%i%i t=%%9.3e (%%%ii steps) maxErr=%%9.2e, ||u||=%%9.2e, cpu=%%9.2e(s)\n",orderOfAccuracyInTime,orderOfAccuracy,numDigits);
                printF(myFormat,t,step,maxErr,solutionNorm,cpuTime);
            }
            else
            {
                sPrintF(myFormat,"cgWave:FD%i%i t=%%9.3e (%%%ii steps) ||u||=%%9.2e, cpu=%%9.2e(s)\n",orderOfAccuracyInTime,orderOfAccuracy,numDigits);
                printF(myFormat,t,step,solutionNorm,cpuTime);
            }

      //   if( numberOfTimeSteps < 1e4 )
      //   {
      //     if( computeErrors )
      //       printF("cgWave: t=%9.3e (%4i steps) maxErr=%9.2e, ||u||=%9.2e \n",t,step,maxErr,solutionNorm);
      //     else
      //       printF("cgWave: t=%9.3e (%4i steps) ||u||=%9.2e \n",t,step,solutionNorm);
      //   }
      //   else if( numberOfTimeSteps < 1e7 )
      //   {
      //     if( computeErrors )
      //       printF("cgWave: t=%9.3e (%6i steps) maxErr=%9.2e, ||u||=%9.2e \n",t,step,maxErr,solutionNorm);
      //     else
      //       printF("cgWave: t=%9.3e (%6i steps) ||u||=%9.2e \n",t,step,solutionNorm);
      //   }
      //   else
      //   {
      //     if( computeErrors )        
      //       printF("cgWave: t=%9.3e (%9i steps) maxErr=%9.2e, ||u||=%9.2e \n",t,step,maxErr,solutionNorm);
      //     else
      //       printF("cgWave: t=%9.3e (%9i steps) ||u||=%9.2e \n",t,step,solutionNorm);
      //   }
      // }

      // output results (e.g. print errors to the check file)
            outputResults( cur, t );

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


    // if( saveShowFile && (i % showSteps == 0) )  // save solution every 'showSteps' steps
    // {
    //  show.startFrame();                                         // start a new frame
    //  show.saveComment(0,sPrintF(buff,"Wave equation"));
    //  show.saveComment(1,sPrintF(buff,"t=%5.2f c=%3.1f ad4=%3.1f",t,c,ad4));
    //  show.saveSolution( u1 );                                        // save the current grid function
    // }


        
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
                ::display(u[prev][grid],sPrintF("BEFORE: u=prev grid=%d t=%9.3e",grid,t+dt),debugFile,"%6.2f ");
                ::display(u[cur ][grid],sPrintF("BEFORE: u=Cur grid=%d t=%9.3e",grid,t+dt),debugFile,"%6.2f ");
            }

          const int useUpwindDissipation = ad4>0.;  // ** FIX ME**

     // useUpwindDissipation : true if grid is implicit AND upwind dissipation in ON 
          int useImplicitUpwindDissipation=false;
          if( gridIsImplicit(grid) && ad4>0. )
              useImplicitUpwindDissipation=true;


      // Add upwinding at least 2 steps in a row 
            if( step>0 && useUpwindDissipation &&
                          ( (step % dissipationFrequency)==0 ) ||
                          ( (step % dissipationFrequency)==1 ) )
            {
        // ---- upwind dissipation ------

                const int prev2= (cur-2+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

                OV_GET_SERIAL_ARRAY(real,u[prev2][grid],umLocal);
                OV_GET_SERIAL_ARRAY(real,u[prev ][grid],uLocal);
                OV_GET_SERIAL_ARRAY(real,u[cur  ][grid],unLocal);

                int option= 1;    // 0=advance, 1=add dissipation
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
                                            useImplicitUpwindDissipation   // ipar[11]
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
                    rpar[12]=cImp(-1);
                    rpar[13]=cImp( 0);
                    rpar[14]=cImp( 1);
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
            // We need v for upwinding or 4th-order curvilinear to hold Lap_2h(u) 
                        RealArray uDot;
                        if( option==1 || ( orderOfAccuracy==4 && !isRectangular) )
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
            }

      // -- optimized advance ---
            OV_GET_SERIAL_ARRAY(real,u[prev][grid],umLocal);
            OV_GET_SERIAL_ARRAY(real,u[cur ][grid],uLocal);
            OV_GET_SERIAL_ARRAY(real,u[next][grid],unLocal);

            int option= 0;    // 0=advance, 1=add dissipation

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
                                        useImplicitUpwindDissipation   // ipar[11]
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
                rpar[12]=cImp(-1);
                rpar[13]=cImp( 0);
                rpar[14]=cImp( 1);
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
          // We need v for upwinding or 4th-order curvilinear to hold Lap_2h(u) 
                    RealArray uDot;
                    if( option==1 || ( orderOfAccuracy==4 && !isRectangular) )
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
            

            if( debug & 16 )
            {
                ::display(u[next][grid],sPrintF("AFTER advOpt: uNext grid=%d t=%9.3e",grid,t+dt),debugFile,"%6.2f ");
            }

            
        }  // end for grid

        if( timeSteppingMethod == implicitTimeStepping )
        {
            takeImplictStep( t+dt );
        }

        t+=dt;

        if( debug & 4 ) printf("...done advance, now interpolate...\n");

                
    // int next = (step+1) %2;
        applyBoundaryConditions( u[next],t);
        

        if( debug & 16 )
        {
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
                ::display(u[next][grid],sPrintF("after applyBC uNext grid=%d t=%9.3e",grid,t),debugFile,"%6.2f ");
        }

        if( computeTimeIntegral )
        {
      // Integral in time -- trapezoidal rule
      // assert( i<=(numberOfTimeSteps-1) );
            
      // StepOptionEnum stepOption = i== (numberOfTimeSteps-1) ? lastStep : middleStep;
            StepOptionEnum stepOption = (t >tFinal-.5*dt) ? lastStep : middleStep;
      // printF(" computeTimeIntegral : t=%9.3e, stepOption=%d\n",t,stepOption);
            updateTimeIntegral( stepOption, t, u[next] );

              

        }
        
    } // for i .. number of steps
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



