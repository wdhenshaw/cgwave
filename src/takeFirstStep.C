// This file automatically generated from takeFirstStep.bC with bpp.
#include "CgWave.h"
#include "CompositeGridOperators.h"  
#include "display.h"
#include "ParallelUtility.h"



//=================================================================================================
// Macro: Add forcing to the first step update
//=================================================================================================


//=================================================================================================
/// \brief Take the first step using Taylor series in time, or other approach
//=================================================================================================
int CgWave::
takeFirstStep( int cur, real t )
{
    const int & solveHelmholtz              = dbase.get<int>("solveHelmholtz");
    const ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
    const int & takeImplicitFirstStep         = dbase.get<int>("takeImplicitFirstStep");
    const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

    bool takeFirstStepImplicit = timeSteppingMethod==implicitTimeStepping &&  takeImplicitFirstStep; 

    const bool useTakeFirstStepHelmholtz = solveHelmholtz || forcingOption==helmholtzForcing;
  // const bool useTakeFirstStepHelmholtz = solveHelmholtz;

    if( useTakeFirstStepHelmholtz )
    {
        return takeFirstStepHelmholtz( cur, t );
    }


    assert( t==0. );

    const int & debug                   = dbase.get<int>("debug");
    const real & c                      = dbase.get<real>("c");
    const real & dt                     = dbase.get<real>("dt");
    const real & omega                  = dbase.get<real>("omega");
    const int & orderOfAccuracy         = dbase.get<int>("orderOfAccuracy");
    const Real & damp                   = dbase.get<Real>("damp");
    const IntegerArray & gridIsImplicit = dbase.get<IntegerArray>("gridIsImplicit");

    FILE *& debugFile  = dbase.get<FILE*>("debugFile");
    FILE *& pDebugFile = dbase.get<FILE*>("pDebugFile");

    if( debug & 4 )
        printF("\n *******  CgWave::takeFirstStep t=%9.3e, cur=%d, dt=%9.3e *************\n",t,cur,dt);

  
    const int & addForcing = dbase.get<int>("addForcing");
    const bool twilightZone = addForcing && forcingOption==twilightZoneForcing ; 

    const int numberOfDimensions = cg.numberOfDimensions();

    const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");    
    const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
    const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;


    realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
    realCompositeGridFunction & up = u[prev];    // previous time 
    realCompositeGridFunction & uc = u[cur];     // current time 
    realCompositeGridFunction & un = u[next];    // next time

  // forcing: 
    realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f");

    CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");

  // --- Get the time-derivative at t=0 and save in up ---
    if( debug & 1 )
        printF("*******  CgWave::takeFirstStep; get TIME DERIVATIVE OF INITIAL CONDITIONS t=%9.3e, cur=%d, dt=%9.3e *************\n",t,cur,dt);
    bool getTimeDerivative=true; 
    getInitialConditions( prev,t,getTimeDerivative );


  // if( takeImplicitFirstStep )
  // {
  //   // --- implicit first step:
  //   //   Save the time derivative in up, this is used in advWave

  //   printF("\n **** takeFirstStep: IMPLICIT FIRST STEP : save time-derivative at t=0 *****\n\n");
  //   return 0;
  // }

      if( debug>0 )
      {
        if( takeImplicitFirstStep )
            printF("\n $$$$$$$ CgWave:takeFirstStep: take implicit first step -- eval u.t(x,0) to use with implicit time step $$$$$$$\n\n");
        else
            printF("\n $$$$$$$ CgWave:takeFirstStep: Take First Step using Taylor Series approach $$$$$$$\n\n");
      }
    const Real c2=c*c, c3=c2*c, c4=c2*c2, c6=c4*c2;
    const Real dt2=dt*dt, dt3=dt2*dt, dt4=dt3*dt, dt5=dt4*dt, dt6=dt5*dt; 
    const Real cdt = c*dt, cdt2=cdt*cdt, cdt3=cdt2*cdt, cdt4=cdt3*cdt, cdt6=cdt3*cdt3; // powers of c*dt 

    Index I1,I2,I3;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        if( gridIsImplicit(grid) && takeImplicitFirstStep )
        {
      // --- implicit first step:
      //   Save the time derivative in up, this is used in advWave
            if( debug>2 )
                printF("\n **** takeFirstStep: IMPLICIT FIRST STEP : save time-derivative at t=0 for grid=%d *****\n\n",grid);

            continue;
        }

        MappedGrid & mg = cg[grid];
        OV_GET_SERIAL_ARRAY(Real,up[grid],upLocal);
        OV_GET_SERIAL_ARRAY(Real,uc[grid],ucLocal);
        OV_GET_SERIAL_ARRAY(Real,un[grid],unLocal);

        getIndex(mg.gridIndexRange(),I1,I2,I3);
        const int includeGhost=1;
        bool ok=ParallelUtility::getLocalArrayBounds(un[grid],unLocal,I1,I2,I3,includeGhost);

        Index J1,J2,J3; // add extra so we can compute Delta^2 
        int extra=orderOfAccuracy/2; 
        getIndex(mg.gridIndexRange(),J1,J2,J3,extra);
        ok=ParallelUtility::getLocalArrayBounds(un[grid],unLocal,J1,J2,J3,includeGhost);

        if( ok )
        {
        
      // --- Taylor series in time for the first step ---- **CHECK ME**
      // -- take a FORWARD STEP ---
      // u(t-dt) = u(t) + dt*ut + (dt^2/2)*utt + (dt^3/6)*uttt + (dt^4/4!)*utttt + (dt^5/5!)*u[5] + dt^6/6!*u[6]


      //  utt = c^2*Delta(u) + f
      //  uttt = c^2*Delta(ut) + ft 
      //  utttt = c^2*Delta(utt) + ftt
      //        = (c^2*Delta)^2 u + c^2*Delta(f) + ftt 

      // TZ:
      //     f = ue.tt - c^2*Delta( ue )     

      // FINISH ME FOR WHEN THERE IS FORCING
      // if( forcingOption ~
            if( forcingOption==userForcing || forcingOption==helmholtzForcing )
                printF("**** takeFirstStep: EXPLICIT FIRST STEP -- FIX ME : ADD FORCING forcingOption=%d\n",(int)forcingOption);

            RealArray lap(J1,J2,J3);
            operators[grid].derivative(MappedGridOperators::laplacianOperator,ucLocal,lap,J1,J2,J3);
  
            unLocal(I1,I2,I3) = ucLocal(I1,I2,I3) + dt*upLocal(I1,I2,I3) + (.5*cdt2)*lap(I1,I2,I3);
            if( damp!=0. )
            {
        // u.tt = c^2 Delta(u) - damp* u.t 
        // u.t = up 
                unLocal(I1,I2,I3) += (-.5*dt*dt*damp)*upLocal(I1,I2,I3);
            }

      // if( false )
      // {
      //   printF("takeFirstStep **TEMP** t+dt=%9.3e \n",t+dt);
      //   assert( dbase.get<OGFunction*>("tz")!=NULL );
      //   OGFunction & e = *dbase.get<OGFunction*>("tz");
      //   OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);
      //   Range C=Range(0,0);
      //   int isRectangular=0;
      //   e.gd( ucLocal ,xLocal,mg.numberOfDimensions(),isRectangular,0,0,0,0,I1,I2,I3,C,t ); // u
      //   e.gd( upLocal ,xLocal,mg.numberOfDimensions(),isRectangular,1,0,0,0,I1,I2,I3,C,t ); // ut 
      //   unLocal(I1,I2,I3) = ucLocal(I1,I2,I3) + dt*upLocal(I1,I2,I3) + (.5*c*c*dt*dt)*lap(I1,I2,I3);
      //   // unLocal(I1,I2,I3) = ucLocal(I1,I2,I3) + dt*upLocal(I1,I2,I3) + (.5*cdt2)*lap(I1,I2,I3);
      // }

      // if( false )
      // {
      //   printF("takeFirstStep **TEMP** set to TZ exact at t+dt=%9.3e \n",t+dt);
      //   assert( dbase.get<OGFunction*>("tz")!=NULL );
      //   OGFunction & e = *dbase.get<OGFunction*>("tz");
      //   OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);
      //   Range C=Range(0,0);
      //   int isRectangular=0;
      //   e.gd( unLocal ,xLocal,mg.numberOfDimensions(),isRectangular,0,0,0,0,I1,I2,I3,C,t+dt);
      // }

            if( orderOfAccuracy>=4 )
            { 

                if( damp!=0. )
                {
                    OV_ABORT("takeFirstStep: ORDER>=4: FINISH ME FOR damping");
                }

                operators[grid].derivative(MappedGridOperators::laplacianOperator,upLocal,lap,I1,I2,I3); // Delta( ut )
                unLocal(I1,I2,I3) += ( c2*dt3/6.)*lap(I1,I2,I3); 

        // --  compute Delta^2 to lower order order ---
                RealArray lapSq(I1,I2,I3);
                operators.setOrderOfAccuracy(orderOfAccuracy-2);

                operators[grid].derivative(MappedGridOperators::laplacianOperator,ucLocal,lap  ,J1,J2,J3);   // Delta(uc) 
                operators[grid].derivative(MappedGridOperators::laplacianOperator,lap    ,lapSq,I1,I2,I3);   // Delta^2( uc )
  
                unLocal(I1,I2,I3) += ( cdt4/24. )*lapSq(I1,I2,I3); 


                operators.setOrderOfAccuracy(orderOfAccuracy); // reset 
            }
            
            if( orderOfAccuracy>=6 )
            { 
                if( damp!=0. )
                {
                    OV_ABORT("takeFirstStep: ORDER>=6: FINISH ME FOR damping");
                }

        // --  compute Delta^2 ( ut ) to fourth order ---
        // *check me* is this stencil too wide?

                RealArray lapSq(J1,J2,J3);
                operators.setOrderOfAccuracy(orderOfAccuracy-2);

                operators[grid].derivative(MappedGridOperators::laplacianOperator,upLocal,lap  ,J1,J2,J3);   // Delta( ut )
                operators[grid].derivative(MappedGridOperators::laplacianOperator,lap    ,lapSq,I1,I2,I3);   // Delta^2( ut )
  
                unLocal(I1,I2,I3) += ( c2*dt5/120. )*lapSq(I1,I2,I3); 

        // --  compute Delta^3( u )  to second order ---
                operators.setOrderOfAccuracy(orderOfAccuracy-4);

                operators[grid].derivative(MappedGridOperators::laplacianOperator,ucLocal,lap  ,J1,J2,J3);   // Delta(uc) 
                operators[grid].derivative(MappedGridOperators::laplacianOperator,lap    ,lapSq,J1,J2,J3);   // Delta^2( uc )
                operators[grid].derivative(MappedGridOperators::laplacianOperator,lapSq  ,lap  ,I1,I2,I3);   // Delta^3( uc )
  
                unLocal(I1,I2,I3) += ( cdt6/720. )*lap(I1,I2,I3);         


                operators.setOrderOfAccuracy(orderOfAccuracy); // reset 
            }


            if( orderOfAccuracy>6 )
            {
                printF("CgWave::takeFirstStep:WARNING: orderOfAccuracy=%d not implemented. Fix me!\n",orderOfAccuracy);
        // OV_ABORT("error");
            }

      // --- adjust update for any forcing ----
                if( twilightZone )
                {
          // --- add forcing for twilight zone ---
          // TZ:
          //     f = ue.tt - c^2*Delta( ue )   
          // 
          //  utt = c^2*Delta(u) + f
          //  uttt = c^2*Delta(ut) + ft 
          //  utttt = c^2*Delta(utt) + ftt
          //        = (c^2*Delta)^2 u + c^2*Delta(f) + ftt          
                    assert( dbase.get<OGFunction*>("tz")!=NULL );
                    OGFunction & e = *dbase.get<OGFunction*>("tz");
                    OV_GET_SERIAL_ARRAY(Real,mg.vertex(),xLocal);
                    Range C=Range(0,0);
                    int isRectangular=0;
                    RealArray f(I1,I2,I3), utt(I1,I2,I3), ut(I1,I2,I3), uxx(I1,I2,I3), uyy(I1,I2,I3), uzz(I1,I2,I3);
                    e.gd( utt,xLocal,numberOfDimensions,isRectangular,2,0,0,0,I1,I2,I3,C,t ); 
                    e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,2,0,0,I1,I2,I3,C,t ); 
                    e.gd( uyy,xLocal,numberOfDimensions,isRectangular,0,0,2,0,I1,I2,I3,C,t );
          // ----------- CORRECT u_tt with f = ue_ttt - c^2*Lap(ue )
                    if( numberOfDimensions==2 )
                    {
                        f = (.5*dt2)*( utt - (c*c)*( uxx + uyy ) );
                    }
                    else
                    {
                        e.gd( uzz,xLocal,numberOfDimensions,isRectangular,0,0,0,2,I1,I2,I3,C,t );
                        f = (.5*dt2)*( utt - (c*c)*( uxx + uyy + uzz ) );
                    }
                    if( damp!=0. )
                    {
                        e.gd( ut,xLocal,numberOfDimensions,isRectangular,1,0,0,0,I1,I2,I3,C,t ); 
                        f += (.5*dt2*damp)*ut; 
                    }    
                    if( orderOfAccuracy>=4 )
                    {
            // ----------- CORRECT u_ttt with ft = ue_tttt - c^2*Lap( uet )
            // Compute ft  = ( utt - (c*c)*( uxx + uyy + uzz ) ).t
                        e.gd( utt,xLocal,numberOfDimensions,isRectangular,3,0,0,0,I1,I2,I3,C,t ); // u.ttt
                        e.gd( uxx,xLocal,numberOfDimensions,isRectangular,1,2,0,0,I1,I2,I3,C,t ); // u.txx
                        e.gd( uyy,xLocal,numberOfDimensions,isRectangular,1,0,2,0,I1,I2,I3,C,t ); // u.tyy  
                        if( numberOfDimensions==2 )
                        {
                            f += (dt3/6.)*( utt - (c*c)*( uxx + uyy ) );
                        }
                        else
                        {
                            e.gd( uzz,xLocal,numberOfDimensions,isRectangular,1,0,0,2,I1,I2,I3,C,t ); // u.tzz
                            f += (dt3/6.)*( utt - (c*c)*( uxx + uyy + uzz ) );
                        }  
            // ----------- CORRECT u_tttt with  
            //       f_tt + c^2*Lap( f ) = ue_tttt - c^4 Delta^2( ue )
            // where 
            //     ftt  = ue_ttttt - c^2*Lap( uett )
            //     Lap( f )  =  Delta( ue.tt - c^2*Delta( ue ) )
            // *new way* Feb 5, 2022 -- account for terms cancelling to simplify
                        e.gd( utt,xLocal,numberOfDimensions,isRectangular,4,0,0,0,I1,I2,I3,C,t ); // really utttt
                        RealArray w(I1,I2,I3);
                        e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,4,0,0,I1,I2,I3,C,t );  // u.xxxx
                        e.gd( uyy,xLocal,numberOfDimensions,isRectangular,0,0,4,0,I1,I2,I3,C,t );  // u.yyyy
                        w = uxx + uyy;
                        e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,2,2,0,I1,I2,I3,C,t );  // u.xxyy
                        w += 2.*uxx; 
                        if( numberOfDimensions==3 )
                        {
                            e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,0,0,4,I1,I2,I3,C,t ); // u.zzzz
                            w += uxx;
                            e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,2,0,2,I1,I2,I3,C,t ); // u.xxzz
                            e.gd( uyy,xLocal,numberOfDimensions,isRectangular,0,0,2,2,I1,I2,I3,C,t ); // u.yyzz
                            w += 2.*( uxx + uyy );
                        } 
                        f += (dt4/24.)*( utt - c4*w );  // f += uetttt - c^4* Delta^2 ( ue )
            // else
            // {
            //   // **OLD WAY** DID not account for terms cancelling
            //   e.gd( utt,xLocal,numberOfDimensions,isRectangular,4,0,0,0,I1,I2,I3,C,t ); 
            //   e.gd( uxx,xLocal,numberOfDimensions,isRectangular,2,2,0,0,I1,I2,I3,C,t ); 
            //   e.gd( uyy,xLocal,numberOfDimensions,isRectangular,2,0,2,0,I1,I2,I3,C,t );  
            //   if( numberOfDimensions==2 )
            //   {
            //     f += (dt4/24.)*( utt - (c*c)*( uxx + uyy ) );
            //   }
            //   else
            //   {
            //     e.gd( uzz,xLocal,numberOfDimensions,isRectangular,2,0,0,2,I1,I2,I3,C,t );
            //     f += (dt4/24.)*( utt - (c*c)*( uxx + uyy + uzz ) );
            //   } 
            //   // compute Delta( f ) = Delta( ut.tt - c^2*Delta( ue ) )
            //   RealArray w(I1,I2,I3);
            //   // utt <- Delta( utt )
            //   e.gd( utt,xLocal,numberOfDimensions,isRectangular,2,2,0,0,I1,I2,I3,C,t ); // u.ttxx
            //   e.gd(   w,xLocal,numberOfDimensions,isRectangular,2,0,2,0,I1,I2,I3,C,t ); // u.ttyy
            //   utt += w;
            //   if( numberOfDimensions==3 )
            //   {
            //     e.gd(   w,xLocal,numberOfDimensions,isRectangular,2,0,0,2,I1,I2,I3,C,t ); // uttzz
            //     utt += w;
            //   }
            //   // uxx <- Delta( uxx )
            //   e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,4,0,0,I1,I2,I3,C,t );  // u.xxxx
            //   e.gd(   w,xLocal,numberOfDimensions,isRectangular,0,2,2,0,I1,I2,I3,C,t );  // u.xxyy
            //   uxx += w;
            //   if( numberOfDimensions==3 )
            //   {
            //     e.gd( w,xLocal,numberOfDimensions,isRectangular,0,2,0,2,I1,I2,I3,C,t );
            //     uxx += w;
            //   } 
            //   // uyy <- Delta( uyy )
            //   e.gd( uyy,xLocal,numberOfDimensions,isRectangular,0,2,2,0,I1,I2,I3,C,t );  // u.xxyy
            //   e.gd(   w,xLocal,numberOfDimensions,isRectangular,0,0,4,0,I1,I2,I3,C,t );  // u.yyyy
            //   uyy += w;
            //   if( numberOfDimensions==3 )
            //   {
            //     e.gd( w,xLocal,numberOfDimensions,isRectangular,0,0,2,2,I1,I2,I3,C,t );
            //     uyy += w;
            //   }  
            //   if( numberOfDimensions==3 )
            //   {
            //     // uzz <- Delta( uzz )
            //     e.gd( uzz,xLocal,numberOfDimensions,isRectangular,0,2,0,2,I1,I2,I3,C,t );  // u.xxzz
            //     e.gd(   w,xLocal,numberOfDimensions,isRectangular,0,0,2,2,I1,I2,I3,C,t );  // u.yyzz
            //     uzz += w;
            //     e.gd(   w,xLocal,numberOfDimensions,isRectangular,0,0,0,4,I1,I2,I3,C,t );  // u.zzzz
            //     uzz += w;
            //   }                              
            //   if( numberOfDimensions==2 )
            //   {
            //     f += (c*c*dt4/24.)*( utt - (c*c)*( uxx + uyy ) );  //  += (dt^4/24)*( c^2*Delta( f ) )
            //   }
            //   else
            //   {
            //     f += (c*c*dt4/24.)*( utt - (c*c)*( uxx + uyy + uzz ) );
            //   } 
            // }
                    }
                    if( orderOfAccuracy>=6 )
                    {
            // ----- correct D_t^5 u  ------
            // with  D_t^5 ue - c^4 Delta^2( ue_t )
                        e.gd( utt,xLocal,numberOfDimensions,isRectangular,5,0,0,0,I1,I2,I3,C,t ); // really uttttt
                        RealArray w(I1,I2,I3);
                        e.gd( uxx,xLocal,numberOfDimensions,isRectangular,1,4,0,0,I1,I2,I3,C,t );  // u.xxxxt
                        e.gd( uyy,xLocal,numberOfDimensions,isRectangular,1,0,4,0,I1,I2,I3,C,t );  // u.yyyyt
                        w = uxx + uyy;
                        e.gd( uxx,xLocal,numberOfDimensions,isRectangular,1,2,2,0,I1,I2,I3,C,t );  // u.xxyyt
                        w += 2.*uxx; 
                        if( numberOfDimensions==3 )
                        {
                            e.gd( uxx,xLocal,numberOfDimensions,isRectangular,1,0,0,4,I1,I2,I3,C,t ); // u.zzzzt
                            w += uxx;
                            e.gd( uxx,xLocal,numberOfDimensions,isRectangular,1,2,0,2,I1,I2,I3,C,t ); // u.xxzzt
                            e.gd( uyy,xLocal,numberOfDimensions,isRectangular,1,0,2,2,I1,I2,I3,C,t ); // u.yyzzt
                            w += 2.*( uxx + uyy );
                        } 
                        f += (dt5/120.)*( utt - c4*w );  // f += uettttt - c^4* Delta^2 ( uet )
            // ------- correct D_t^6 u with ------
            //      D_t^6 ue - c^6 Delta^3( ue )
                        e.gd( utt,xLocal,numberOfDimensions,isRectangular,6,0,0,0,I1,I2,I3,C,t ); // really D_t^6 u
                        e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,6,0,0,I1,I2,I3,C,t );  // u.xxxxxx
                        e.gd( uyy,xLocal,numberOfDimensions,isRectangular,0,0,6,0,I1,I2,I3,C,t );  // u.yyyyyy
                        w = uxx + uyy;
                        e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,4,2,0,I1,I2,I3,C,t );  // u.xxxxyy
                        e.gd( uyy,xLocal,numberOfDimensions,isRectangular,0,2,4,0,I1,I2,I3,C,t );  // u.xxyyyy
                        w += 3.*( uxx + uyy ); 
                        if( numberOfDimensions==3 )
                        {
                            e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,0,0,6,I1,I2,I3,C,t ); // u.zzzzzz
                            w += uxx;
                            e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,4,0,2,I1,I2,I3,C,t ); // u.xxxxzz
                            e.gd( uyy,xLocal,numberOfDimensions,isRectangular,0,2,0,4,I1,I2,I3,C,t ); // u.xxzzzz
                            w += 3.*( uxx + uyy );
                            e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,0,4,2,I1,I2,I3,C,t ); // u.yyyyzz
                            e.gd( uyy,xLocal,numberOfDimensions,isRectangular,0,0,2,4,I1,I2,I3,C,t ); // u.yyzzzz
                            w += 3.*( uxx + uyy );  
                            e.gd( uxx,xLocal,numberOfDimensions,isRectangular,0,2,2,2,I1,I2,I3,C,t ); // u.xxyyzz
                            w += 6.*( uxx );
                        } 
                        f += (dt6/720.)*( utt - c6*w );  // f +=   D_t^6 ue - c^6 Delta^3( ue )
            // // f = ue.tt - c^2*Delta( ue )         
            // // Compute fttt 
            // e.gd( utt,xLocal,numberOfDimensions,isRectangular,5,0,0,0,I1,I2,I3,C,t );  // ue.ttttt
            // e.gd( uxx,xLocal,numberOfDimensions,isRectangular,3,2,0,0,I1,I2,I3,C,t );  // ue.tttxx
            // e.gd( uyy,xLocal,numberOfDimensions,isRectangular,3,0,2,0,I1,I2,I3,C,t );  
            // if( numberOfDimensions==2 )
            // {
            //   f += (dt5/120.)*( uttt - (c*c)*( uxx + uyy ) );
            // }
            // else
            // {
            //   e.gd( uzz,xLocal,numberOfDimensions,isRectangular,3,0,0,2,I1,I2,I3,C,t );
            //   f += (dt5/120.)*( uttt - (c*c)*( uxx + uyy + uzz ) );
            // } 
            // // f = ue.tt - c^2*Delta( ue )         
            // // Compute ftttt 
            // e.gd( utt,xLocal,numberOfDimensions,isRectangular,6,0,0,0,I1,I2,I3,C,t );  // ue.tttttt
            // e.gd( uxx,xLocal,numberOfDimensions,isRectangular,4,2,0,0,I1,I2,I3,C,t );  // ue.ttttxx
            // e.gd( uyy,xLocal,numberOfDimensions,isRectangular,4,0,2,0,I1,I2,I3,C,t );  
            // if( numberOfDimensions==2 )
            // {
            //   f += (dt6/720.)*( uttt - (c*c)*( uxx + uyy ) );
            // }
            // else
            // {
            //   e.gd( uzz,xLocal,numberOfDimensions,isRectangular,4,0,0,2,I1,I2,I3,C,t );
            //   f += (dt6/720.)*( uttt - (c*c)*( uxx + uyy + uzz ) );
            // } 
                    }
                    if( orderOfAccuracy>=8 )
                    {
                        printF("\n ***** TakeFirstStep : add forcing: finish me for orderOfAccuracy=%d ***** \n\n",orderOfAccuracy);
                    }
                    unLocal(I1,I2,I3) += f;
                }  
  

        }

    } // end for grid

  if( debug & 8 )
    {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            ::display(u[next][grid],sPrintF("\nAFTER FIRST-STEP, BEFORE BCs: u  on grid=%d t=%9.3e",grid,t+dt),debugFile,"%10.2e ");

        fprintf(debugFile,"\n *** ERRORS AFTER TAKE FIRST STEP, BEFORE APPLY BCs\n\n");
        getErrors( u[next], t+dt );
    }  

    if( debug & 4 )
        printF("*******  CgWave::takeFirstStep; APPLY BOUNDARY CONDTIONS t=%9.3e, cur=%d, dt=%9.3e *************\n",t,cur,dt);


    bool applyExplicitBoundaryConditions=true;
    applyBoundaryConditions( u[next],u[cur], t+dt,applyExplicitBoundaryConditions  );
    

    if( debug & 8 )
    {
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            ::display(u[next][grid],sPrintF("\nAFTER FIRST-STEP: u  on grid=%d t=%9.3e",grid,t+dt),debugFile,"%10.2e ");

        fprintf(debugFile,"\n *** ERRORS AFTER TAKE FIRST STEP AND APPLY BCs\n\n");
        getErrors( u[next], t+dt );
    }

  // u[next].display(sPrintF("u[next] after first step, t=%9.3e",t+dt),"%6.2f ");


    return 0;

}





//=================================================================================================
/// \brief Take the first step using Taylor series in time (e.g. for Helmholtz solve) 
/// THIS ASSUMES A HELMHOLTZ SOLVE OR HELMHOLTZ FORCING 
//=================================================================================================
int CgWave::
takeFirstStepHelmholtz( int cur, real t )
{

    assert( t==0. );

    const int & debug           = dbase.get<int>("debug");
    if( debug & 4 )
        printF("*******  CgWave::takeFirstStepHelmholtz GET SOLUTION at t=dt *************\n");
    
    const real & c                        = dbase.get<real>("c");
    const real & dt                       = dbase.get<real>("dt");
  // const real & omega                    = dbase.get<real>("omega");
    const int & orderOfAccuracy           = dbase.get<int>("orderOfAccuracy");
    const int & orderOfAccuracyInTime     = dbase.get<int>("orderOfAccuracyInTime");

    const int & numberOfFrequencies       = dbase.get<int>("numberOfFrequencies");
    const RealArray & frequencyArray      = dbase.get<RealArray>("frequencyArray");  
    const RealArray & frequencyArraySave  = dbase.get<RealArray>("frequencyArraySave");  

    const int & solveHelmholtz            = dbase.get<int>("solveHelmholtz");
    const int & filterTimeDerivative      = dbase.get<int>("filterTimeDerivative");
    const int & takeImplicitFirstStep     = dbase.get<int>("takeImplicitFirstStep");
    const int & useSuperGrid              = dbase.get<int>("useSuperGrid");

    ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");
    const TimeSteppingMethodEnum & timeSteppingMethod = dbase.get<TimeSteppingMethodEnum>("timeSteppingMethod");

  // Do Helmholtz case for now: 

    assert( forcingOption==helmholtzForcing );

  //  ---- NOTE: change sign of forcing for Helmholtz since we want to solve ----
  //       ( omega^2 I + c^2 Delta) w = f  
    Real fSign = forcingOption==helmholtzForcing ? -1.0 : 1.0;

  // bool usePeriodicFirstStep= true;

    bool usePeriodicFirstStep = false; // *wdh* CHANGED : Sept 10, 2024 -- now this matches the theory in the overHoltz paper

  // *wdh* March 25, 2023 -- we cannot assume omega-periodic solution if we are computing eigen-modes
  // But there seems to be trouble using this !?
    const int & computeEigenmodes = dbase.get<int>("computeEigenmodes");  
    if( computeEigenmodes &&
            timeSteppingMethod == explicitTimeStepping)
    {
        usePeriodicFirstStep=false;
        fSign=0; // there is no forcing when computing eigenmodes.
        printF("+++++++++++++ takeFirstStepHelmholtz: do NOT use periodic first step ++++++++++++\n");
    }


    if( ( computeEigenmodes || numberOfFrequencies>1 ) && timeSteppingMethod == implicitTimeStepping )
    {
        usePeriodicFirstStep=true; 
    }

    if( timeSteppingMethod == implicitTimeStepping && takeImplicitFirstStep )
    {
        usePeriodicFirstStep=false; // added April 29, 2025 -- we can do explicit/implicit first steps now
        printF("+++++++++++++ takeFirstStepHelmholtz: TAKE IMPLICIT FIRST STEP ++++++++++++\n");
    }
    else if( filterTimeDerivative ) 
    {
        usePeriodicFirstStep=true; // DO this for COMPLEX solutions **Sept 28, 2024 **
    }

    

    if( filterTimeDerivative && useSuperGrid && takeImplicitFirstStep ) 
    {
        printF(">>> takeFirstStepHelmholtz:ERROR: takeImplicitFirstStep with SUPERGRID not working yet **FIX ME** ++++++++++++\n");
        OV_ABORT("error");

    }

    if( usePeriodicFirstStep )
      printF("+++++++++++++ takeFirstStepHelmholtz: USE periodic first step ++++++++++++\n");

    const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");    
    const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
    const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;


    realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
    realCompositeGridFunction & up = u[prev];    // previous time 
    realCompositeGridFunction & uc = u[cur];     // current time 
    realCompositeGridFunction & un = u[next];    // next time

  // forcing: 
    realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f");

    CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");

  // ---- WaveHoltz: initial condition is provided in v: ----

    realCompositeGridFunction & v = dbase.get<realCompositeGridFunction>("v");

    const Real & viFactor = dbase.get<Real>("viFactor");
    const Real omega  = frequencyArray(0);
    const Real scaleD0t = (sin(omega*dt)/dt)/viFactor;      

    Index I1,I2,I3;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        OV_GET_SERIAL_ARRAY(Real,un[grid],unLocal);

        if( usePeriodicFirstStep )
        {
      // do all points including ghost 
            getIndex(mg.dimension(),I1,I2,I3);
            int includeParallelGhost=1;
            bool ok=ParallelUtility::getLocalArrayBounds(un[grid],unLocal,I1,I2,I3,includeParallelGhost);

            if( solveHelmholtz )
            {
                OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);

                if( ok )
                {      
           // // ---> Set  u(dt) = u(0)*cos(omega*(dt)) if time-periodic
                    if( debug & 2 )
                        printF("takeFirstStepHelmholtz: grid=%d, setting  u(dt) = u(0)*cos(omega*(dt)) , dt=%20.12e, omega=%g, frequencyArray(0)=%g\n",
                                grid,dt,omega,frequencyArray(0));

           // Real diff = max(fabs(ucLocal(I1,I2,I3)-vLocal(I1,I2,I3,0)));
           // printF(" >> diff | uc - v |=%9.2e\n",diff);

                    if( filterTimeDerivative==0 )
                    {
                        unLocal(I1,I2,I3) = vLocal(I1,I2,I3,0)*cos( frequencyArray(0)*dt );
                        for( int freq=1; freq<numberOfFrequencies; freq++ )
                        {
                          unLocal(I1,I2,I3) += vLocal(I1,I2,I3,freq)*cos( frequencyArray(freq)*dt );
                        }
                    }
                    else
                    {
            // Periodic solution is 
            //   w(x,t) = vk*cos(omega*t) + (1/omega)*vDotk*sin(omega*t)
            // 
            // This includes the fix for the time-discretization error
            //   un(J) = vk(J,1)*cos(par.frequencyArray(1)*dt) + isign*vk(J,2)*dt; 
            // const Real isign=-1.; // FIX ME : for exp( isign*I*omega*t )

                        const Real isign=+1.; // FIX ME : for exp( isign*I*omega*t )

            // **NOTE**
            // viFactor must match in
            //     getFilterWeights
            //     takeFirstStep
            //     solveHelmholtz : definition of Im(u)
                        if( true )
                        { // *new way* June 19, 2023
                            const Real & viFactor = dbase.get<Real>("viFactor");
                            if( true || debug & 4 )
                            {
                                printF("\n @@@@@@ takeFirstStepHelmholtz: use periodic first step for filterTimeDerivative=%d frequencyArray(0)=%12.4e viFactor=%9.2e @@@@@@@@ \n\n",
                                    filterTimeDerivative,frequencyArray(0),viFactor);

                // ::display(vLocal,"vLocal","%5.2f ");

                // Real vMax = max(fabs(vLocal(I1,I2,I3,1)));
                // printF(" max(fabs(v(:,:,:,1))=%9.2e\n",vMax);
                            }

                            unLocal(I1,I2,I3) = vLocal(I1,I2,I3,0)*cos(frequencyArray(0)*dt ) + viFactor*vLocal(I1,I2,I3,1)*sin(frequencyArray(0)*dt);
                        }
                        else if( true || timeSteppingMethod==explicitTimeStepping )
                        {
                            Real viFactor = isign/frequencyArray(0);     
                            unLocal(I1,I2,I3) = vLocal(I1,I2,I3,0)*cos(frequencyArray(0)*dt ) + viFactor*vLocal(I1,I2,I3,1)*sin(frequencyArray(0)*dt);
                        }
                        else
                        {
              // viFactor = dt/sin(om*dt) from filter weights
                            unLocal(I1,I2,I3) = vLocal(I1,I2,I3,0)*cos(frequencyArray(0)*dt ) + isign*vLocal(I1,I2,I3,1)*dt;
                        }
                    }


                }
            }
            else
            {
                OV_GET_SERIAL_ARRAY(Real,uc[grid],ucLocal);

                if( ok )
                {      
           // // ---> Set  u(dt) = u(0)*cos(omega*(dt)) if time-periodic
                      printF("takeFirstStepHelmholtz: setting  u(dt) = u(0)*cos(omega*(dt)) , dt=%20.12e, omega=%g\n",dt,omega);

                      unLocal(I1,I2,I3)  = ucLocal(I1,I2,I3) * cos(omega*dt);
                }
            }

        }
        else if( timeSteppingMethod == implicitTimeStepping )
        {
      // -----IMPLICIT FIRST STEP : COMPLEX CASE ----



            if( debug>2 )
                printF("\n **** takeFirstStep: IMPLICIT FIRST STEP : save v(:,:,:,1) in up(:,:,:) for use in implicit first step *****\n\n");   

            Index I1,I2,I3;
            for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            {
        //   Save v(:,:,:,1) in up(:,:,:), this is used in advWave
  
                MappedGrid & mg = cg[grid];
                OV_GET_SERIAL_ARRAY(Real,up[grid],upLocal);
                OV_GET_SERIAL_ARRAY(Real,v[grid],vLocal);
                getIndex(mg.dimension(),I1,I2,I3);
                int includeParallelGhost=1;
                bool ok=ParallelUtility::getLocalArrayBounds(v[grid],vLocal,I1,I2,I3,includeParallelGhost);
                if( ok )
                {
                    if( filterTimeDerivative )
                    {
                        upLocal(I1,I2,I3) = vLocal(I1,I2,I3,1)*scaleD0t; // this is really D0t W^0 
                    }
                    else
                    {
                        upLocal(I1,I2,I3) = 0.; // assume initial time derivative is zero
                    }
                }

            }


        }
        else
        {
      // ---- EXPLICIT FIRST STEP -----

            OV_GET_SERIAL_ARRAY(Real,uc[grid],ucLocal);
            OV_GET_SERIAL_ARRAY(Real,f[grid],fLocal);

            RealDistributedArray & vg = filterTimeDerivative ? v[grid] : uc[grid];
            OV_GET_SERIAL_ARRAY(Real,vg,vLocal);

            getIndex(mg.gridIndexRange(),I1,I2,I3);
            bool ok=ParallelUtility::getLocalArrayBounds(un[grid],unLocal,I1,I2,I3);
            if( ok )
            {
            
        // ----- Taylor series first step -------

        // *NOTE* This same first step occurs if we impose 
        //    D+tD-t U^0 = Lph U^0 + force
        //    Dzt U^0 = 0 

                RealArray lap(I1,I2,I3);
                operators[grid].derivative(MappedGridOperators::laplacianOperator,ucLocal,lap,I1,I2,I3);
      
        // -- take a FORWARD STEP ---
        // u(t-dt) = u(t) + dt*ut + (dt^2/2)*utt + (dt^3/6)*uttt + (dt^4/4!)*utttt
        //  utt = c^2*Delta(u) + f
        //  uttt = c^2*Delta(ut) + ft 
        //  utttt = c^2*Delta(utt) + ftt
        //        = (c^2*Delta)^2 u + c^2*Delta(f) + ftt 

                if( filterTimeDerivative )
                {
                    if( debug>2 )
                        printF("~~~~~~~~~~~~~~ Take EXPLICIT first step Complex case : omega=%14.6e\n",omega);
                    unLocal(I1,I2,I3)  = ucLocal(I1,I2,I3) + dt*vLocal(I1,I2,I3,1)*scaleD0t + (.5*dt*dt*c*c)*lap(I1,I2,I3) + (.5*dt*dt *cos(omega*t)*fSign)*fLocal(I1,I2,I3);
                }
                else
                {

                    if( numberOfFrequencies==1 )
                    {
                        if( debug>2 )
                            printF("~~~~~~~~~~~~~~ Take EXPLICIT first step assuming D0t U  = 0 : omega=%14.6e\n",omega);
                        unLocal(I1,I2,I3)  = ucLocal(I1,I2,I3) + (.5*dt*dt*c*c)*lap(I1,I2,I3) + (.5*dt*dt *cos(omega*t)*fSign)*fLocal(I1,I2,I3);
                    }
                    else
                    {  // *WDH* May 15, 2025
                        if( debug>2 )
                            printF("~~~~~~~~~~~~~~ Take EXPLICIT first step assuming D0t U  = 0 ~~~~~~~~~~~~~\n");
                        unLocal(I1,I2,I3)  = ucLocal(I1,I2,I3) + (.5*dt*dt*c*c)*lap(I1,I2,I3);
                        for( int freq=0; freq<numberOfFrequencies; freq++ )
                        {
                            unLocal(I1,I2,I3) += (.5*dt*dt *cos(frequencyArray(freq)*t)*fSign)*fLocal(I1,I2,I3,freq);
                        }
                    }
                }

        // if( orderOfAccuracy==4 ) // *WDH* May 15, 2025
                if( orderOfAccuracyInTime==4 )
                {
          // this may be good enough for 4th-order -- local error is dt^4
                    real DeltaUt =0.; // we assume ut=0
          // upLocal(I1,I2,I3) += ( (dt*dt*dt/6.)*(-omega*sin(omega*t) )*( (c*c)*DeltaUt + fLocal(I1,I2,I3) )
                    unLocal(I1,I2,I3) += ( fSign*(dt*dt*dt/6.)*(-omega*sin(omega*t) ) )*( fLocal(I1,I2,I3) );
                }

            }
            
        }
    } // end for grid

    if( !usePeriodicFirstStep && 
            timeSteppingMethod == explicitTimeStepping ) // May 2, 2025 added after implicit first step added for complex case
    {

        bool applyExplicitBoundaryConditions=true;
        applyBoundaryConditions( u[next],u[cur], t+dt, applyExplicitBoundaryConditions );
    }
    else
    {
    // This next is needed for implicit time-stepping and eigenWave 
        for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
            u[next][grid].updateGhostBoundaries();
    }
    
  // u[next].display(sPrintF("u[next] after first step, t=%9.3e",t+dt),"%6.2f ");


    return 0;

}



//=================================================================================================
/// \brief Take the first BACKWARD step using Taylor series in time (e.g. for Helmholtz solve) 
/// THIS ASSUMES A HELMHOLTZ SOLVE OR HELMHOLTZ FORCING 
//=================================================================================================
int CgWave::
takeFirstBackwardStep( int cur, real t )
{

    const int & debug           = dbase.get<int>("debug");
    if( debug & 4 )
        printF("*******  CgWave::takeFirstBackwardStep GET SOLUTION at -dt *************\n");
    
    const real & c              = dbase.get<real>("c");
    const real & dt             = dbase.get<real>("dt");
    const real & omega          = dbase.get<real>("omega");
    const int & orderOfAccuracy = dbase.get<int>("orderOfAccuracy");

    ForcingOptionEnum & forcingOption = dbase.get<ForcingOptionEnum>("forcingOption");

  // Do Helmholtz case for now: 
    assert( forcingOption==helmholtzForcing );

    //  ---- NOTE: change sign of forcing for Helmholtz since we want to solve ----
  //       ( omega^2 I + c^2 Delta) w = f  
    const Real fSign = forcingOption==helmholtzForcing ? -1.0 : 1.0;

  // const int & solveHelmholtz = dbase.get<int>("solveHelmholtz");
  // const Real fSign = solveHelmholtz ? -1.0 : 1.0;  
    

    const int & numberOfTimeLevelsStored = dbase.get<int>("numberOfTimeLevelsStored");    
    const int prev= (cur-1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;
  // const int next= (cur+1+numberOfTimeLevelsStored) % numberOfTimeLevelsStored;

    realCompositeGridFunction *& u = dbase.get<realCompositeGridFunction*>("ucg");
    realCompositeGridFunction & un = u[cur];     // current time 
    realCompositeGridFunction & up = u[prev];    // previous time

  // forcing: 
    realCompositeGridFunction & f = dbase.get<realCompositeGridFunction>("f");

    CompositeGridOperators & operators = dbase.get<CompositeGridOperators>("operators");

    Index I1,I2,I3;
    for( int grid=0; grid<cg.numberOfComponentGrids(); grid++ )
    {
        MappedGrid & mg = cg[grid];
        OV_GET_SERIAL_ARRAY(Real,un[grid],unLocal);
        OV_GET_SERIAL_ARRAY(Real,up[grid],upLocal);
        OV_GET_SERIAL_ARRAY(Real,f[grid],fLocal);

        getIndex(mg.gridIndexRange(),I1,I2,I3);
        bool ok=ParallelUtility::getLocalArrayBounds(un[grid],unLocal,I1,I2,I3);
        if( ok )
        {
            bool usePeriodicFirstStep=false;
            if( usePeriodicFirstStep )
            {
         // **TESTING** u(-dt) = u(0)*cos(omega*(-dt)) if time-periodic
                printF("takeFirstBackwardStep: setting  u(-dt) = u(0)*cos(omega*(-dt)) , dt=%20.12e\n",dt);
                upLocal(I1,I2,I3)  = unLocal(I1,I2,I3) * cos(-omega*dt);
            }
            else
            {
                RealArray lap(I1,I2,I3);
                operators[grid].derivative(MappedGridOperators::laplacianOperator,unLocal,lap,I1,I2,I3);
      
        // -- take a BACKWARD STEP ---
        // u(t-dt) = u(t) - dt*ut + (dt^2/2)*utt - (dt^3/6)*uttt + (dt^4/4!)*utttt
        //  utt = c^2*Delta(u) + f
        //  uttt = c^2*Delta(ut) + ft 
        //  utttt = c^2*Delta(utt) + ftt
        //        = (c^2*Delta)^2 u + c^2*Delta(f) + ftt 
                upLocal(I1,I2,I3)  = unLocal(I1,I2,I3) -dt*(0.) + (.5*dt*dt*c*c)*lap(I1,I2,I3) + (.5*dt*dt *cos(omega*t)*fSign)*fLocal(I1,I2,I3);
                if( orderOfAccuracy==4 )
                {
          // this may be good enough for 4th-order -- local erro is dt^4
                    real DeltaUt =0.; // we assume ut=0
          // upLocal(I1,I2,I3) += ( -(dt*dt*dt/6.)*(-omega*sin(omega*t) )*( (c*c)*DeltaUt + fLocal(I1,I2,I3) )
                    upLocal(I1,I2,I3) += ( -fSign*(dt*dt*dt/6.)*(-omega*sin(omega*t) ) )*( fLocal(I1,I2,I3) );
                }

            }
            
        }
    } // end for grid

    bool applyExplicitBoundaryConditions=true;
    applyBoundaryConditions( u[prev],u[cur], t-dt, applyExplicitBoundaryConditions );
    

    return 0;

}