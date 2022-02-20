! This file automatically generated from advWaveME.bf90 with bpp.
        subroutine advWaveME2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
    ! subroutine advWaveME2dOrder8r(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,!                 mask,xy,rsxy,  um,u,un, f,fa, v, vh,  bc, frequencyArray, ipar, rpar, ierr )
   !======================================================================
   !   Advance a time step for Waves equations
   !
   ! nd : number of space dimensions
   ! um,u,un : u(t-dt), u(t), u(t+dt)
   !
   ! ipar(0)  = option : option=0 - advance wave equation
   !                           =1 - add upwind dissipation (predictor corrector mode)
   !
   !======================================================================
        implicit none
        integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b
        real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
        real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
        real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
        real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
        real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
        real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1) 
        real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
        real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
        real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
        integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
        integer bc(0:1,0:2),ierr
        integer gridIndexRange(0:1,0:2)
        real frequencyArray(0:*)
        integer ipar(0:*)
        real rpar(0:*)
   !     ---- local variables -----
        integer m1a,m1b,m2a,m2b,m3a,m3b,numGhost,nStart,nEnd,mt
        integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime,axis,dir,grid,freq
        integer addForcing,orderOfDissipation,option,gridIsImplicit,preComputeUpwindUt
        integer useNewForcingMethod,numberOfForcingFunctions,fcur,fnext,fprev,numberOfFrequencies
        real t,tm,cc,dt,dy,dz,cdt,cdtdx,cdtdy,cdtdz
    ! ,adc,adcdt,add,adddt
        real dt4by12
    ! logical addDissipation
        integer debug
        integer adjustHelmholtzForUpwinding
        real dx(0:2),dr(0:2)
    ! real dx2i,dy2i,dz2i,dxsqi,dysqi,dzsqi,dxi,dyi,dzi
    ! real dx12i,dy12i,dz12i,dxsq12i,dysq12i,dzsq12i,dxy4i,dxz4i,dyz4,time0,time1
    ! real dxi4,dyi4,dzi4,dxdyi2,dxdzi2,dydzi2
        real c0,c1,csq,dtsq,cdtsq,cdtsq12,cdtSqBy12
        integer maxOrderOfAccuracy
        parameter( maxOrderOfAccuracy=12 )
    ! Coefficients in the implicit scheme
        real bImp(0:maxOrderOfAccuracy-1)
        real cImp(-1:1,0:maxOrderOfAccuracy-1)
        real alpha2,alpha4,alpha6,alpha8, beta2,beta4,beta6,beta8
        integer rectangular,curvilinear
        parameter( rectangular=0, curvilinear=1 )
        integer timeSteppingMethod
        integer defaultTimeStepping,adamsSymmetricOrder3,rungeKuttaFourthOrder,stoermerTimeStepping,modifiedEquationTimeStepping
        parameter(defaultTimeStepping=0,adamsSymmetricOrder3=1,rungeKuttaFourthOrder=2,stoermerTimeStepping=3,modifiedEquationTimeStepping=4)
   !- ! Dispersion models
   !- integer noDispersion,drude
   !- parameter( noDispersion=0, drude=1 )
   !...........start statement function
        integer kd,m
    ! real rx,ry,rz,sx,sy,sz,tx,ty,tz
     ! --- declare variables used in the difference approximations (defined below) ---
     ! declareDifferenceOrder2(u,RX)
     ! declareDifferenceOrder2(um,none)
     ! declareDifferenceOrder2(v,none)
     ! declareDifferenceOrder2(f,none)
     ! declareDifferenceOrder4(u,RX)
     ! declareDifferenceOrder4(um,none)
     ! declareDifferenceOrder4(v,none)
     ! declareDifferenceOrder6(u,RX)
     ! declareDifferenceOrder6(um,none)
     ! declareDifferenceOrder6(v,none)
     ! declareDifferenceOrder8(u,RX)
     ! declareDifferenceOrder8(um,none)
     ! declareDifferenceOrder8(v,none)
    ! define variables for getDerivatives macros
    ! #Include "../maple/declareGetSixthDerivativesMacrosVariables.h"
        real cdtPow2,cdtPow4By12,cdtPow6By360,cdtPow8By20160
    ! real lap2d2,lap3d2,lap2d4,lap3d4,lap2d6,lap3d6,lap2d8,lap3d8,lap2d2Pow2,lap3d2Pow2,lap2d2Pow3,lap3d2Pow3,!      lap2d2Pow4,lap3d2Pow4,lap2d4Pow2,lap3d4Pow2,lap2d4Pow3,lap3d4Pow3,lap2d6Pow2,lap3d6Pow2
    ! real lap2d2m,lap3d2m
    ! real du,fd22d,fd23d,fd42d,fd43d,fd62d,fd63d,fd82d,fd83d
    ! real DztU
    ! forcing correction functions: 
    ! real lap2d2f,f2drme44, lap3d2f, f3drme44, f2dcme44, f3dcme44, 
        real ff
    ! real cdSosupx,cdSosupy,cdSosupz
        real adSosup,sosupParameter, uDotFactor, adxSosup(0:2)
        integer useSosupDissipation,sosupDissipationOption
        integer updateSolution,updateDissipation,computeUt
        integer ec 
        real ep 
        real fv(0:1) , ev(0:1), evtt(0:1), evxx(0:1), evyy(0:1), evzz(0:1)
        real evxxxx(0:1), evxxyy(0:1), evyyyy(0:1), evxxzz(0:1), evyyzz(0:1), evzzzz(0:1), evtttt(0:1)
        real evtttttt(0:1)
        real evxxxxxx(0:1)
        real evyyyyyy(0:1)
        real evzzzzzz(0:1)       
        real evxxyyyy(0:1)
        real evxxxxyy(0:1)
        real evxxxxzz(0:1)
        real evxxzzzz(0:1)
        real evyyyyzz(0:1)
        real evyyzzzz(0:1)
        real evxxyyzz(0:1)
        real omega, coswt
        integer maxFreq
        parameter( maxFreq=500 )
        real cosFreqt(0:maxFreq), coswtAve(0:maxFreq), cosineFactor(0:maxFreq)
        integer idv(0:2),j1,j2,j3
        integer iStencil,upwCase,upwindHalfStencilWidth,i1l,i2l,i3l, i1r,i2r,i3r
        integer useUpwindDissipation,useImplicitUpwindDissipation,adjustOmega,solveHelmholtz
        real upw,maxDiff,umj
    ! real upwindCoeff(-3:3,0:3) 
        integer forcingOption
    ! forcingOptions -- these should match ForcingEnum in CgWave.h 
    ! enum ForcingOptionEnum
    ! {
    !   noForcing=0,
    !   twilightZoneForcing,
    !   userForcing,
    !   helmholtzForcing
    ! };
        integer noForcing,twilightZoneForcing,userForcing,helmholtzForcing
        parameter(noForcing           =0,twilightZoneForcing =1,userForcing         =2,helmholtzForcing    =3 )
     ! ******* Variables for hierarchical method *******
     ! real ud(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:30)
     ! real ude(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:30)
          real maxErr(1:30), l2Err(1:30)
          real maxSol(30)
          real ue
          real uet8 
          real uex8 
          real uey8 
          real uez8 
          real uex6y2
          real uex4y4
          real uex2y6
          real uex6z2
          real uex4z4
          real uex2z6
          real uey6z2
          real uey4z4
          real uey2z6
          real uex4y2z2
          real uex2y4z2
          real uex2y2z4
     ! --- LOCAL arrays ---
          real d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d600(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d060(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d420(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d240(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d400i,d040i
          real d600i,d060i, d240i, d420i , d510i,d150i,d330i
          real d800i,d080i, d620i, d260i, d440i
          real d100i, d010i, d300i,d030i, d500i, d050i, d700i, d070i
          real d110i, d310i,d130i, d710i,d170i,d350i,d530i
          real d120i,d210i, d140i,d410i, d320i,d230i, d160i, d610i, d520i, d250i, d340i, d430i
     ! end if gridtype == curvilinear
          real uex,uey
          real uett,uexx,uexy,ueyy,uezz,ueLap
          real uexxx,ueyyy,uezzz,uexxy,uexyy
          real uexxxx,ueyyyy,uezzzz,uexxyy, uexxxy, uexyyy
          real uexxxxx,   uexxxxy,   uexxxyy,   uexxyyy,   uexyyyy,   ueyyyyy
          real uexxxxxx,  uexxxxxy,  uexxxxyy,  uexxxyyy,  uexxyyyy,  uexyyyyy,  ueyyyyyy
          real uexxxxxxx, uexxxxxxy, uexxxxxyy, uexxxxyyy, uexxxyyyy, uexxyyyyy, uexyyyyyy, ueyyyyyyy
          real uexxxxxxxx,uexxxxxxxy,uexxxxxxyy,uexxxxxyyy,uexxxxyyyy,uexxxyyyyy,uexxyyyyyy,uexyyyyyyy,ueyyyyyyyy
          real ux,uy,uz,uxx,uxy,uyy,uxz,uyz,uzz,uLap
          real uxxx,     uyyy,     uxxy,     uxyy
          real uxxxx,    uyyyy,    uxxyy,    uxxxy,    uxyyy
          real uxxxxx,   uxxxxy,   uxxxyy,   uxxyyy,   uxyyyy,   uyyyyy
          real uxxxxxx,  uxxxxxy,  uxxxxyy,  uxxxyyy,  uxxyyyy,  uxyyyyy,  uyyyyyy
          real uxxxxxxx ,uxxxxxxy, uxxxxxyy, uxxxxyyy, uxxxyyyy, uxxyyyyy, uxyyyyyy, uyyyyyyy
          real uxxxxxxxx,uxxxxxxxy,uxxxxxxyy,uxxxxxyyy,uxxxxyyyy,uxxxyyyyy,uxxyyyyyy,uxyyyyyyy,uyyyyyyyy
          integer maxDeriv,d,uc,count,numGhost1,m1,m2,m3
! Coefficients from ??
    real cux100
    real cux010
    real cuy100
    real cuy010

    real cuxx100 
    real cuxx200 
    real cuxx010 
    real cuxx110 
    real cuxx020 
    real cuxy100 
    real cuxy200 
    real cuxy010 
    real cuxy110 
    real cuxy020 
    real cuyy100 
    real cuyy200 
    real cuyy010 
    real cuyy110 
    real cuyy020 

    real cuxxx100
    real cuxxx200
    real cuxxx300
    real cuxxx010
    real cuxxx110
    real cuxxx210
    real cuxxx020
    real cuxxx120
    real cuxxx030
    real cuxxy100
    real cuxxy200
    real cuxxy300
    real cuxxy010
    real cuxxy110
    real cuxxy210
    real cuxxy020
    real cuxxy120
    real cuxxy030
    real cuxyy100
    real cuxyy200
    real cuxyy300
    real cuxyy010
    real cuxyy110
    real cuxyy210
    real cuxyy020
    real cuxyy120
    real cuxyy030
    real cuyyy100
    real cuyyy200
    real cuyyy300
    real cuyyy010
    real cuyyy110
    real cuyyy210
    real cuyyy020
    real cuyyy120
    real cuyyy030

    real cuxxxx100 
    real cuxxxx200 
    real cuxxxx300 
    real cuxxxx400 
    real cuxxxx010 
    real cuxxxx110 
    real cuxxxx210 
    real cuxxxx310 
    real cuxxxx020 
    real cuxxxx120 
    real cuxxxx220 
    real cuxxxx030 
    real cuxxxx130 
    real cuxxxx040 
    real cuxxxy100 
    real cuxxxy200 
    real cuxxxy300 
    real cuxxxy400 
    real cuxxxy010 
    real cuxxxy110 
    real cuxxxy210 
    real cuxxxy310 
    real cuxxxy020 
    real cuxxxy120 
    real cuxxxy220 
    real cuxxxy030 
    real cuxxxy130 
    real cuxxxy040 
    real cuxxyy100 
    real cuxxyy200 
    real cuxxyy300 
    real cuxxyy400 
    real cuxxyy010 
    real cuxxyy110 
    real cuxxyy210 
    real cuxxyy310 
    real cuxxyy020 
    real cuxxyy120 
    real cuxxyy220 
    real cuxxyy030 
    real cuxxyy130 
    real cuxxyy040 
    real cuxyyy100 
    real cuxyyy200 
    real cuxyyy300 
    real cuxyyy400 
    real cuxyyy010 
    real cuxyyy110 
    real cuxyyy210 
    real cuxyyy310 
    real cuxyyy020 
    real cuxyyy120 
    real cuxyyy220 
    real cuxyyy030 
    real cuxyyy130 
    real cuxyyy040 
    real cuyyyy100 
    real cuyyyy200 
    real cuyyyy300 
    real cuyyyy400 
    real cuyyyy010 
    real cuyyyy110 
    real cuyyyy210 
    real cuyyyy310 
    real cuyyyy020 
    real cuyyyy120 
    real cuyyyy220 
    real cuyyyy030 
    real cuyyyy130
    real cuyyyy040

    real cuxxxxx100
    real cuxxxxx200
    real cuxxxxx300
    real cuxxxxx400
    real cuxxxxx500
    real cuxxxxx010
    real cuxxxxx110
    real cuxxxxx210
    real cuxxxxx310
    real cuxxxxx410
    real cuxxxxx020
    real cuxxxxx120
    real cuxxxxx220
    real cuxxxxx320
    real cuxxxxx030
    real cuxxxxx130
    real cuxxxxx230
    real cuxxxxx040
    real cuxxxxx140
    real cuxxxxx050
    real cuxxxxy100
    real cuxxxxy200
    real cuxxxxy300
    real cuxxxxy400
    real cuxxxxy500
    real cuxxxxy010
    real cuxxxxy110
    real cuxxxxy210
    real cuxxxxy310
    real cuxxxxy410
    real cuxxxxy020
    real cuxxxxy120
    real cuxxxxy220
    real cuxxxxy320
    real cuxxxxy030
    real cuxxxxy130
    real cuxxxxy230
    real cuxxxxy040
    real cuxxxxy140
    real cuxxxxy050
    real cuxxxyy100
    real cuxxxyy200
    real cuxxxyy300
    real cuxxxyy400
    real cuxxxyy500
    real cuxxxyy010
    real cuxxxyy110
    real cuxxxyy210
    real cuxxxyy310
    real cuxxxyy410
    real cuxxxyy020
    real cuxxxyy120
    real cuxxxyy220
    real cuxxxyy320
    real cuxxxyy030
    real cuxxxyy130
    real cuxxxyy230
    real cuxxxyy040
    real cuxxxyy140
    real cuxxxyy050
    real cuxxyyy100
    real cuxxyyy200
    real cuxxyyy300
    real cuxxyyy400
    real cuxxyyy500
    real cuxxyyy010
    real cuxxyyy110
    real cuxxyyy210
    real cuxxyyy310
    real cuxxyyy410
    real cuxxyyy020
    real cuxxyyy120
    real cuxxyyy220
    real cuxxyyy320
    real cuxxyyy030
    real cuxxyyy130
    real cuxxyyy230
    real cuxxyyy040
    real cuxxyyy140
    real cuxxyyy050
    real cuxyyyy100
    real cuxyyyy200
    real cuxyyyy300
    real cuxyyyy400
    real cuxyyyy500
    real cuxyyyy010
    real cuxyyyy110
    real cuxyyyy210
    real cuxyyyy310
    real cuxyyyy410
    real cuxyyyy020
    real cuxyyyy120
    real cuxyyyy220
    real cuxyyyy320
    real cuxyyyy030
    real cuxyyyy130
    real cuxyyyy230
    real cuxyyyy040
    real cuxyyyy140
    real cuxyyyy050
    real cuyyyyy100
    real cuyyyyy200
    real cuyyyyy300
    real cuyyyyy400
    real cuyyyyy500
    real cuyyyyy010
    real cuyyyyy110
    real cuyyyyy210
    real cuyyyyy310
    real cuyyyyy410
    real cuyyyyy020
    real cuyyyyy120
    real cuyyyyy220
    real cuyyyyy320
    real cuyyyyy030
    real cuyyyyy130
    real cuyyyyy230
    real cuyyyyy040
    real cuyyyyy140
    real cuyyyyy050

    real cuxxxxxx100
    real cuxxxxxx200
    real cuxxxxxx300
    real cuxxxxxx400
    real cuxxxxxx500
    real cuxxxxxx600
    real cuxxxxxx010
    real cuxxxxxx110
    real cuxxxxxx210
    real cuxxxxxx310
    real cuxxxxxx410
    real cuxxxxxx510
    real cuxxxxxx020
    real cuxxxxxx120
    real cuxxxxxx220
    real cuxxxxxx320
    real cuxxxxxx420
    real cuxxxxxx030
    real cuxxxxxx130
    real cuxxxxxx230
    real cuxxxxxx330
    real cuxxxxxx040
    real cuxxxxxx140
    real cuxxxxxx240
    real cuxxxxxx050
    real cuxxxxxx150
    real cuxxxxxx060
    real cuxxxxxy100
    real cuxxxxxy200
    real cuxxxxxy300
    real cuxxxxxy400
    real cuxxxxxy500
    real cuxxxxxy600
    real cuxxxxxy010
    real cuxxxxxy110
    real cuxxxxxy210
    real cuxxxxxy310
    real cuxxxxxy410
    real cuxxxxxy510
    real cuxxxxxy020
    real cuxxxxxy120
    real cuxxxxxy220
    real cuxxxxxy320
    real cuxxxxxy420
    real cuxxxxxy030
    real cuxxxxxy130
    real cuxxxxxy230
    real cuxxxxxy330
    real cuxxxxxy040
    real cuxxxxxy140
    real cuxxxxxy240
    real cuxxxxxy050
    real cuxxxxxy150
    real cuxxxxxy060
    real cuxxxxyy100
    real cuxxxxyy200
    real cuxxxxyy300
    real cuxxxxyy400
    real cuxxxxyy500
    real cuxxxxyy600
    real cuxxxxyy010
    real cuxxxxyy110
    real cuxxxxyy210
    real cuxxxxyy310
    real cuxxxxyy410
    real cuxxxxyy510
    real cuxxxxyy020
    real cuxxxxyy120
    real cuxxxxyy220
    real cuxxxxyy320
    real cuxxxxyy420
    real cuxxxxyy030
    real cuxxxxyy130
    real cuxxxxyy230
    real cuxxxxyy330
    real cuxxxxyy040
    real cuxxxxyy140
    real cuxxxxyy240
    real cuxxxxyy050
    real cuxxxxyy150
    real cuxxxxyy060
    real cuxxxyyy100
    real cuxxxyyy200
    real cuxxxyyy300
    real cuxxxyyy400
    real cuxxxyyy500
    real cuxxxyyy600
    real cuxxxyyy010
    real cuxxxyyy110
    real cuxxxyyy210
    real cuxxxyyy310
    real cuxxxyyy410
    real cuxxxyyy510
    real cuxxxyyy020
    real cuxxxyyy120
    real cuxxxyyy220
    real cuxxxyyy320
    real cuxxxyyy420
    real cuxxxyyy030
    real cuxxxyyy130
    real cuxxxyyy230
    real cuxxxyyy330
    real cuxxxyyy040
    real cuxxxyyy140
    real cuxxxyyy240
    real cuxxxyyy050
    real cuxxxyyy150
    real cuxxxyyy060
    real cuxxyyyy100
    real cuxxyyyy200
    real cuxxyyyy300
    real cuxxyyyy400
    real cuxxyyyy500
    real cuxxyyyy600
    real cuxxyyyy010
    real cuxxyyyy110
    real cuxxyyyy210
    real cuxxyyyy310
    real cuxxyyyy410
    real cuxxyyyy510
    real cuxxyyyy020
    real cuxxyyyy120
    real cuxxyyyy220
    real cuxxyyyy320
    real cuxxyyyy420
    real cuxxyyyy030
    real cuxxyyyy130
    real cuxxyyyy230
    real cuxxyyyy330
    real cuxxyyyy040
    real cuxxyyyy140
    real cuxxyyyy240
    real cuxxyyyy050
    real cuxxyyyy150
    real cuxxyyyy060
    real cuxyyyyy100
    real cuxyyyyy200
    real cuxyyyyy300
    real cuxyyyyy400
    real cuxyyyyy500
    real cuxyyyyy600
    real cuxyyyyy010
    real cuxyyyyy110
    real cuxyyyyy210
    real cuxyyyyy310
    real cuxyyyyy410
    real cuxyyyyy510
    real cuxyyyyy020
    real cuxyyyyy120
    real cuxyyyyy220
    real cuxyyyyy320
    real cuxyyyyy420
    real cuxyyyyy030
    real cuxyyyyy130
    real cuxyyyyy230
    real cuxyyyyy330
    real cuxyyyyy040
    real cuxyyyyy140
    real cuxyyyyy240
    real cuxyyyyy050
    real cuxyyyyy150
    real cuxyyyyy060
    real cuyyyyyy100
    real cuyyyyyy200
    real cuyyyyyy300
    real cuyyyyyy400
    real cuyyyyyy500
    real cuyyyyyy600
    real cuyyyyyy010
    real cuyyyyyy110
    real cuyyyyyy210
    real cuyyyyyy310
    real cuyyyyyy410
    real cuyyyyyy510
    real cuyyyyyy020
    real cuyyyyyy120
    real cuyyyyyy220
    real cuyyyyyy320
    real cuyyyyyy420
    real cuyyyyyy030
    real cuyyyyyy130
    real cuyyyyyy230
    real cuyyyyyy330
    real cuyyyyyy040
    real cuyyyyyy140
    real cuyyyyyy240
    real cuyyyyyy050
    real cuyyyyyy150
    real cuyyyyyy060

    real cuxxxxxxx100
    real cuxxxxxxx200
    real cuxxxxxxx300
    real cuxxxxxxx400
    real cuxxxxxxx500
    real cuxxxxxxx600
    real cuxxxxxxx700
    real cuxxxxxxx010
    real cuxxxxxxx110
    real cuxxxxxxx210
    real cuxxxxxxx310
    real cuxxxxxxx410
    real cuxxxxxxx510
    real cuxxxxxxx610
    real cuxxxxxxx020
    real cuxxxxxxx120
    real cuxxxxxxx220
    real cuxxxxxxx320
    real cuxxxxxxx420
    real cuxxxxxxx520
    real cuxxxxxxx030
    real cuxxxxxxx130
    real cuxxxxxxx230
    real cuxxxxxxx330
    real cuxxxxxxx430
    real cuxxxxxxx040
    real cuxxxxxxx140
    real cuxxxxxxx240
    real cuxxxxxxx340
    real cuxxxxxxx050
    real cuxxxxxxx150
    real cuxxxxxxx250
    real cuxxxxxxx060
    real cuxxxxxxx160
    real cuxxxxxxx070
    real cuxxxxxxy100
    real cuxxxxxxy200
    real cuxxxxxxy300
    real cuxxxxxxy400
    real cuxxxxxxy500
    real cuxxxxxxy600
    real cuxxxxxxy700
    real cuxxxxxxy010
    real cuxxxxxxy110
    real cuxxxxxxy210
    real cuxxxxxxy310
    real cuxxxxxxy410
    real cuxxxxxxy510
    real cuxxxxxxy610
    real cuxxxxxxy020
    real cuxxxxxxy120
    real cuxxxxxxy220
    real cuxxxxxxy320
    real cuxxxxxxy420
    real cuxxxxxxy520
    real cuxxxxxxy030
    real cuxxxxxxy130
    real cuxxxxxxy230
    real cuxxxxxxy330
    real cuxxxxxxy430
    real cuxxxxxxy040
    real cuxxxxxxy140
    real cuxxxxxxy240
    real cuxxxxxxy340
    real cuxxxxxxy050
    real cuxxxxxxy150
    real cuxxxxxxy250
    real cuxxxxxxy060
    real cuxxxxxxy160
    real cuxxxxxxy070
    real cuxxxxxyy100
    real cuxxxxxyy200
    real cuxxxxxyy300
    real cuxxxxxyy400
    real cuxxxxxyy500
    real cuxxxxxyy600
    real cuxxxxxyy700
    real cuxxxxxyy010
    real cuxxxxxyy110
    real cuxxxxxyy210
    real cuxxxxxyy310
    real cuxxxxxyy410
    real cuxxxxxyy510
    real cuxxxxxyy610
    real cuxxxxxyy020
    real cuxxxxxyy120
    real cuxxxxxyy220
    real cuxxxxxyy320
    real cuxxxxxyy420
    real cuxxxxxyy520
    real cuxxxxxyy030
    real cuxxxxxyy130
    real cuxxxxxyy230
    real cuxxxxxyy330
    real cuxxxxxyy430
    real cuxxxxxyy040
    real cuxxxxxyy140
    real cuxxxxxyy240
    real cuxxxxxyy340
    real cuxxxxxyy050
    real cuxxxxxyy150
    real cuxxxxxyy250
    real cuxxxxxyy060
    real cuxxxxxyy160
    real cuxxxxxyy070
    real cuxxxxyyy100
    real cuxxxxyyy200
    real cuxxxxyyy300
    real cuxxxxyyy400
    real cuxxxxyyy500
    real cuxxxxyyy600
    real cuxxxxyyy700
    real cuxxxxyyy010
    real cuxxxxyyy110
    real cuxxxxyyy210
    real cuxxxxyyy310
    real cuxxxxyyy410
    real cuxxxxyyy510
    real cuxxxxyyy610
    real cuxxxxyyy020
    real cuxxxxyyy120
    real cuxxxxyyy220
    real cuxxxxyyy320
    real cuxxxxyyy420
    real cuxxxxyyy520
    real cuxxxxyyy030
    real cuxxxxyyy130
    real cuxxxxyyy230
    real cuxxxxyyy330
    real cuxxxxyyy430
    real cuxxxxyyy040
    real cuxxxxyyy140
    real cuxxxxyyy240
    real cuxxxxyyy340
    real cuxxxxyyy050
    real cuxxxxyyy150
    real cuxxxxyyy250
    real cuxxxxyyy060
    real cuxxxxyyy160
    real cuxxxxyyy070
    real cuxxxyyyy100
    real cuxxxyyyy200
    real cuxxxyyyy300
    real cuxxxyyyy400
    real cuxxxyyyy500
    real cuxxxyyyy600
    real cuxxxyyyy700
    real cuxxxyyyy010
    real cuxxxyyyy110
    real cuxxxyyyy210
    real cuxxxyyyy310
    real cuxxxyyyy410
    real cuxxxyyyy510
    real cuxxxyyyy610
    real cuxxxyyyy020
    real cuxxxyyyy120
    real cuxxxyyyy220
    real cuxxxyyyy320
    real cuxxxyyyy420
    real cuxxxyyyy520
    real cuxxxyyyy030
    real cuxxxyyyy130
    real cuxxxyyyy230
    real cuxxxyyyy330
    real cuxxxyyyy430
    real cuxxxyyyy040
    real cuxxxyyyy140
    real cuxxxyyyy240
    real cuxxxyyyy340
    real cuxxxyyyy050
    real cuxxxyyyy150
    real cuxxxyyyy250
    real cuxxxyyyy060
    real cuxxxyyyy160
    real cuxxxyyyy070
    real cuxxyyyyy100
    real cuxxyyyyy200
    real cuxxyyyyy300
    real cuxxyyyyy400
    real cuxxyyyyy500
    real cuxxyyyyy600
    real cuxxyyyyy700
    real cuxxyyyyy010
    real cuxxyyyyy110
    real cuxxyyyyy210
    real cuxxyyyyy310
    real cuxxyyyyy410
    real cuxxyyyyy510
    real cuxxyyyyy610
    real cuxxyyyyy020
    real cuxxyyyyy120
    real cuxxyyyyy220
    real cuxxyyyyy320
    real cuxxyyyyy420
    real cuxxyyyyy520
    real cuxxyyyyy030
    real cuxxyyyyy130
    real cuxxyyyyy230
    real cuxxyyyyy330
    real cuxxyyyyy430
    real cuxxyyyyy040
    real cuxxyyyyy140
    real cuxxyyyyy240
    real cuxxyyyyy340
    real cuxxyyyyy050
    real cuxxyyyyy150
    real cuxxyyyyy250
    real cuxxyyyyy060
    real cuxxyyyyy160
    real cuxxyyyyy070
    real cuxyyyyyy100
    real cuxyyyyyy200
    real cuxyyyyyy300
    real cuxyyyyyy400
    real cuxyyyyyy500
    real cuxyyyyyy600
    real cuxyyyyyy700
    real cuxyyyyyy010
    real cuxyyyyyy110
    real cuxyyyyyy210
    real cuxyyyyyy310
    real cuxyyyyyy410
    real cuxyyyyyy510
    real cuxyyyyyy610
    real cuxyyyyyy020
    real cuxyyyyyy120
    real cuxyyyyyy220
    real cuxyyyyyy320
    real cuxyyyyyy420
    real cuxyyyyyy520
    real cuxyyyyyy030
    real cuxyyyyyy130
    real cuxyyyyyy230
    real cuxyyyyyy330
    real cuxyyyyyy430
    real cuxyyyyyy040
    real cuxyyyyyy140
    real cuxyyyyyy240
    real cuxyyyyyy340
    real cuxyyyyyy050
    real cuxyyyyyy150
    real cuxyyyyyy250
    real cuxyyyyyy060
    real cuxyyyyyy160
    real cuxyyyyyy070
    real cuyyyyyyy100
    real cuyyyyyyy200
    real cuyyyyyyy300
    real cuyyyyyyy400
    real cuyyyyyyy500
    real cuyyyyyyy600
    real cuyyyyyyy700
    real cuyyyyyyy010
    real cuyyyyyyy110
    real cuyyyyyyy210
    real cuyyyyyyy310
    real cuyyyyyyy410
    real cuyyyyyyy510
    real cuyyyyyyy610
    real cuyyyyyyy020
    real cuyyyyyyy120
    real cuyyyyyyy220
    real cuyyyyyyy320
    real cuyyyyyyy420
    real cuyyyyyyy520
    real cuyyyyyyy030
    real cuyyyyyyy130
    real cuyyyyyyy230
    real cuyyyyyyy330
    real cuyyyyyyy430
    real cuyyyyyyy040
    real cuyyyyyyy140
    real cuyyyyyyy240
    real cuyyyyyyy340
    real cuyyyyyyy050
    real cuyyyyyyy150
    real cuyyyyyyy250
    real cuyyyyyyy060
    real cuyyyyyyy160
    real cuyyyyyyy070


    real cuxxxxxxxx100
    real cuxxxxxxxx200
    real cuxxxxxxxx300
    real cuxxxxxxxx400
    real cuxxxxxxxx500
    real cuxxxxxxxx600
    real cuxxxxxxxx700
    real cuxxxxxxxx800
    real cuxxxxxxxx010
    real cuxxxxxxxx110
    real cuxxxxxxxx210
    real cuxxxxxxxx310
    real cuxxxxxxxx410
    real cuxxxxxxxx510
    real cuxxxxxxxx610
    real cuxxxxxxxx710
    real cuxxxxxxxx020
    real cuxxxxxxxx120
    real cuxxxxxxxx220
    real cuxxxxxxxx320
    real cuxxxxxxxx420
    real cuxxxxxxxx520
    real cuxxxxxxxx620
    real cuxxxxxxxx030
    real cuxxxxxxxx130
    real cuxxxxxxxx230
    real cuxxxxxxxx330
    real cuxxxxxxxx430
    real cuxxxxxxxx530
    real cuxxxxxxxx040
    real cuxxxxxxxx140
    real cuxxxxxxxx240
    real cuxxxxxxxx340
    real cuxxxxxxxx440
    real cuxxxxxxxx050
    real cuxxxxxxxx150
    real cuxxxxxxxx250
    real cuxxxxxxxx350
    real cuxxxxxxxx060
    real cuxxxxxxxx160
    real cuxxxxxxxx260
    real cuxxxxxxxx070
    real cuxxxxxxxx170
    real cuxxxxxxxx080
    real cuxxxxxxxy100
    real cuxxxxxxxy200
    real cuxxxxxxxy300
    real cuxxxxxxxy400
    real cuxxxxxxxy500
    real cuxxxxxxxy600
    real cuxxxxxxxy700
    real cuxxxxxxxy800
    real cuxxxxxxxy010
    real cuxxxxxxxy110
    real cuxxxxxxxy210
    real cuxxxxxxxy310
    real cuxxxxxxxy410
    real cuxxxxxxxy510
    real cuxxxxxxxy610
    real cuxxxxxxxy710
    real cuxxxxxxxy020
    real cuxxxxxxxy120
    real cuxxxxxxxy220
    real cuxxxxxxxy320
    real cuxxxxxxxy420
    real cuxxxxxxxy520
    real cuxxxxxxxy620
    real cuxxxxxxxy030
    real cuxxxxxxxy130
    real cuxxxxxxxy230
    real cuxxxxxxxy330
    real cuxxxxxxxy430
    real cuxxxxxxxy530
    real cuxxxxxxxy040
    real cuxxxxxxxy140
    real cuxxxxxxxy240
    real cuxxxxxxxy340
    real cuxxxxxxxy440
    real cuxxxxxxxy050
    real cuxxxxxxxy150
    real cuxxxxxxxy250
    real cuxxxxxxxy350
    real cuxxxxxxxy060
    real cuxxxxxxxy160
    real cuxxxxxxxy260
    real cuxxxxxxxy070
    real cuxxxxxxxy170
    real cuxxxxxxxy080
    real cuxxxxxxyy100
    real cuxxxxxxyy200
    real cuxxxxxxyy300
    real cuxxxxxxyy400
    real cuxxxxxxyy500
    real cuxxxxxxyy600
    real cuxxxxxxyy700
    real cuxxxxxxyy800
    real cuxxxxxxyy010
    real cuxxxxxxyy110
    real cuxxxxxxyy210
    real cuxxxxxxyy310
    real cuxxxxxxyy410
    real cuxxxxxxyy510
    real cuxxxxxxyy610
    real cuxxxxxxyy710
    real cuxxxxxxyy020
    real cuxxxxxxyy120
    real cuxxxxxxyy220
    real cuxxxxxxyy320
    real cuxxxxxxyy420
    real cuxxxxxxyy520
    real cuxxxxxxyy620
    real cuxxxxxxyy030
    real cuxxxxxxyy130
    real cuxxxxxxyy230
    real cuxxxxxxyy330
    real cuxxxxxxyy430
    real cuxxxxxxyy530
    real cuxxxxxxyy040
    real cuxxxxxxyy140
    real cuxxxxxxyy240
    real cuxxxxxxyy340
    real cuxxxxxxyy440
    real cuxxxxxxyy050
    real cuxxxxxxyy150
    real cuxxxxxxyy250
    real cuxxxxxxyy350
    real cuxxxxxxyy060
    real cuxxxxxxyy160
    real cuxxxxxxyy260
    real cuxxxxxxyy070
    real cuxxxxxxyy170
    real cuxxxxxxyy080
    real cuxxxxxyyy100
    real cuxxxxxyyy200
    real cuxxxxxyyy300
    real cuxxxxxyyy400
    real cuxxxxxyyy500
    real cuxxxxxyyy600
    real cuxxxxxyyy700
    real cuxxxxxyyy800
    real cuxxxxxyyy010
    real cuxxxxxyyy110
    real cuxxxxxyyy210
    real cuxxxxxyyy310
    real cuxxxxxyyy410
    real cuxxxxxyyy510
    real cuxxxxxyyy610
    real cuxxxxxyyy710
    real cuxxxxxyyy020
    real cuxxxxxyyy120
    real cuxxxxxyyy220
    real cuxxxxxyyy320
    real cuxxxxxyyy420
    real cuxxxxxyyy520
    real cuxxxxxyyy620
    real cuxxxxxyyy030
    real cuxxxxxyyy130
    real cuxxxxxyyy230
    real cuxxxxxyyy330
    real cuxxxxxyyy430
    real cuxxxxxyyy530
    real cuxxxxxyyy040
    real cuxxxxxyyy140
    real cuxxxxxyyy240
    real cuxxxxxyyy340
    real cuxxxxxyyy440
    real cuxxxxxyyy050
    real cuxxxxxyyy150
    real cuxxxxxyyy250
    real cuxxxxxyyy350
    real cuxxxxxyyy060
    real cuxxxxxyyy160
    real cuxxxxxyyy260
    real cuxxxxxyyy070
    real cuxxxxxyyy170
    real cuxxxxxyyy080
    real cuxxxxyyyy100
    real cuxxxxyyyy200
    real cuxxxxyyyy300
    real cuxxxxyyyy400
    real cuxxxxyyyy500
    real cuxxxxyyyy600
    real cuxxxxyyyy700
    real cuxxxxyyyy800
    real cuxxxxyyyy010
    real cuxxxxyyyy110
    real cuxxxxyyyy210
    real cuxxxxyyyy310
    real cuxxxxyyyy410
    real cuxxxxyyyy510
    real cuxxxxyyyy610
    real cuxxxxyyyy710
    real cuxxxxyyyy020
    real cuxxxxyyyy120
    real cuxxxxyyyy220
    real cuxxxxyyyy320
    real cuxxxxyyyy420
    real cuxxxxyyyy520
    real cuxxxxyyyy620
    real cuxxxxyyyy030
    real cuxxxxyyyy130
    real cuxxxxyyyy230
    real cuxxxxyyyy330
    real cuxxxxyyyy430
    real cuxxxxyyyy530
    real cuxxxxyyyy040
    real cuxxxxyyyy140
    real cuxxxxyyyy240
    real cuxxxxyyyy340
    real cuxxxxyyyy440
    real cuxxxxyyyy050
    real cuxxxxyyyy150
    real cuxxxxyyyy250
    real cuxxxxyyyy350
    real cuxxxxyyyy060
    real cuxxxxyyyy160
    real cuxxxxyyyy260
    real cuxxxxyyyy070
    real cuxxxxyyyy170
    real cuxxxxyyyy080
    real cuxxxyyyyy100
    real cuxxxyyyyy200
    real cuxxxyyyyy300
    real cuxxxyyyyy400
    real cuxxxyyyyy500
    real cuxxxyyyyy600
    real cuxxxyyyyy700
    real cuxxxyyyyy800
    real cuxxxyyyyy010
    real cuxxxyyyyy110
    real cuxxxyyyyy210
    real cuxxxyyyyy310
    real cuxxxyyyyy410
    real cuxxxyyyyy510
    real cuxxxyyyyy610
    real cuxxxyyyyy710
    real cuxxxyyyyy020
    real cuxxxyyyyy120
    real cuxxxyyyyy220
    real cuxxxyyyyy320
    real cuxxxyyyyy420
    real cuxxxyyyyy520
    real cuxxxyyyyy620
    real cuxxxyyyyy030
    real cuxxxyyyyy130
    real cuxxxyyyyy230
    real cuxxxyyyyy330
    real cuxxxyyyyy430
    real cuxxxyyyyy530
    real cuxxxyyyyy040
    real cuxxxyyyyy140
    real cuxxxyyyyy240
    real cuxxxyyyyy340
    real cuxxxyyyyy440
    real cuxxxyyyyy050
    real cuxxxyyyyy150
    real cuxxxyyyyy250
    real cuxxxyyyyy350
    real cuxxxyyyyy060
    real cuxxxyyyyy160
    real cuxxxyyyyy260
    real cuxxxyyyyy070
    real cuxxxyyyyy170
    real cuxxxyyyyy080
    real cuxxyyyyyy100
    real cuxxyyyyyy200
    real cuxxyyyyyy300
    real cuxxyyyyyy400
    real cuxxyyyyyy500
    real cuxxyyyyyy600
    real cuxxyyyyyy700
    real cuxxyyyyyy800
    real cuxxyyyyyy010
    real cuxxyyyyyy110
    real cuxxyyyyyy210
    real cuxxyyyyyy310
    real cuxxyyyyyy410
    real cuxxyyyyyy510
    real cuxxyyyyyy610
    real cuxxyyyyyy710
    real cuxxyyyyyy020
    real cuxxyyyyyy120
    real cuxxyyyyyy220
    real cuxxyyyyyy320
    real cuxxyyyyyy420
    real cuxxyyyyyy520
    real cuxxyyyyyy620
    real cuxxyyyyyy030
    real cuxxyyyyyy130
    real cuxxyyyyyy230
    real cuxxyyyyyy330
    real cuxxyyyyyy430
    real cuxxyyyyyy530
    real cuxxyyyyyy040
    real cuxxyyyyyy140
    real cuxxyyyyyy240
    real cuxxyyyyyy340
    real cuxxyyyyyy440
    real cuxxyyyyyy050
    real cuxxyyyyyy150
    real cuxxyyyyyy250
    real cuxxyyyyyy350
    real cuxxyyyyyy060
    real cuxxyyyyyy160
    real cuxxyyyyyy260
    real cuxxyyyyyy070
    real cuxxyyyyyy170
    real cuxxyyyyyy080
    real cuxyyyyyyy100
    real cuxyyyyyyy200
    real cuxyyyyyyy300
    real cuxyyyyyyy400
    real cuxyyyyyyy500
    real cuxyyyyyyy600
    real cuxyyyyyyy700
    real cuxyyyyyyy800
    real cuxyyyyyyy010
    real cuxyyyyyyy110
    real cuxyyyyyyy210
    real cuxyyyyyyy310
    real cuxyyyyyyy410
    real cuxyyyyyyy510
    real cuxyyyyyyy610
    real cuxyyyyyyy710
    real cuxyyyyyyy020
    real cuxyyyyyyy120
    real cuxyyyyyyy220
    real cuxyyyyyyy320
    real cuxyyyyyyy420
    real cuxyyyyyyy520
    real cuxyyyyyyy620
    real cuxyyyyyyy030
    real cuxyyyyyyy130
    real cuxyyyyyyy230
    real cuxyyyyyyy330
    real cuxyyyyyyy430
    real cuxyyyyyyy530
    real cuxyyyyyyy040
    real cuxyyyyyyy140
    real cuxyyyyyyy240
    real cuxyyyyyyy340
    real cuxyyyyyyy440
    real cuxyyyyyyy050
    real cuxyyyyyyy150
    real cuxyyyyyyy250
    real cuxyyyyyyy350
    real cuxyyyyyyy060
    real cuxyyyyyyy160
    real cuxyyyyyyy260
    real cuxyyyyyyy070
    real cuxyyyyyyy170
    real cuxyyyyyyy080
    real cuyyyyyyyy100
    real cuyyyyyyyy200
    real cuyyyyyyyy300
    real cuyyyyyyyy400
    real cuyyyyyyyy500
    real cuyyyyyyyy600
    real cuyyyyyyyy700
    real cuyyyyyyyy800
    real cuyyyyyyyy010
    real cuyyyyyyyy110
    real cuyyyyyyyy210
    real cuyyyyyyyy310
    real cuyyyyyyyy410
    real cuyyyyyyyy510
    real cuyyyyyyyy610
    real cuyyyyyyyy710
    real cuyyyyyyyy020
    real cuyyyyyyyy120
    real cuyyyyyyyy220
    real cuyyyyyyyy320
    real cuyyyyyyyy420
    real cuyyyyyyyy520
    real cuyyyyyyyy620
    real cuyyyyyyyy030
    real cuyyyyyyyy130
    real cuyyyyyyyy230
    real cuyyyyyyyy330
    real cuyyyyyyyy430
    real cuyyyyyyyy530
    real cuyyyyyyyy040
    real cuyyyyyyyy140
    real cuyyyyyyyy240
    real cuyyyyyyyy340
    real cuyyyyyyyy440
    real cuyyyyyyyy050
    real cuyyyyyyyy150
    real cuyyyyyyyy250
    real cuyyyyyyyy350
    real cuyyyyyyyy060
    real cuyyyyyyyy160
    real cuyyyyyyyy260
    real cuyyyyyyyy070
    real cuyyyyyyyy170
    real cuyyyyyyyy080
   ! real unxx22r,unyy22r,unxy22r,unx22r
   !.......statement functions for jacobian
    ! rx(i1,i2,i3)=rsxy(i1,i2,i3,0,0)
    ! ry(i1,i2,i3)=rsxy(i1,i2,i3,0,1)
    ! rz(i1,i2,i3)=rsxy(i1,i2,i3,0,2)
    ! sx(i1,i2,i3)=rsxy(i1,i2,i3,1,0)
    ! sy(i1,i2,i3)=rsxy(i1,i2,i3,1,1)
    ! sz(i1,i2,i3)=rsxy(i1,i2,i3,1,2)
    ! tx(i1,i2,i3)=rsxy(i1,i2,i3,2,0)
    ! ty(i1,i2,i3)=rsxy(i1,i2,i3,2,1)
    ! tz(i1,i2,i3)=rsxy(i1,i2,i3,2,2)
    ! !  The next macro will define the difference approximation statement functions for u
    ! defineDifferenceOrder2Components1(u,RX)
    ! defineDifferenceOrder4Components1(u,RX)
    ! defineDifferenceOrder6Components1(u,RX)
    ! defineDifferenceOrder8Components1(u,RX)
    ! ! Difference approximations um (old time)
    ! defineDifferenceOrder2Components1(um,none)
    ! defineDifferenceOrder4Components1(um,none)
    ! ! Define difference approximations for v
    ! defineDifferenceOrder2Components1(v,none)
    ! defineDifferenceOrder4Components1(v,none)
    ! defineDifferenceOrder2Components1(f,none)
     !...........end   statement functions
     ! write(*,*) 'Inside advWave...'
          cc             = rpar( 0)  ! this is c
          dt             = rpar( 1)
          dx(0)          = rpar( 2)
          dx(1)          = rpar( 3)
          dx(2)          = rpar( 4)
          dr(0)          = rpar( 5)
          dr(1)          = rpar( 6)
          dr(2)          = rpar( 7)
          t              = rpar( 8)
          ep             = rpar( 9)
          sosupParameter = rpar(10)
          omega          = rpar(11) ! for helmholtz 
          bImp( 0)       = rpar(12) ! beta2 : coefficient for implicit time-stepping
          bImp( 1)       = rpar(13) ! beta4 : coefficient for implicit time-stepping
          bImp( 2)       = rpar(14) ! beta6 (for future)
          bImp( 3)       = rpar(15) ! beta8 (for future)
          dy=dx(1)  ! Are these needed?
          dz=dx(2)
          option                       = ipar( 0)
          grid                         = ipar( 1)
          gridType                     = ipar( 2)
          orderOfAccuracy              = ipar( 3)
          orderInTime                  = ipar( 4)
          addForcing                   = ipar( 5)
          forcingOption                = ipar( 6)
          numberOfForcingFunctions     = ipar( 7)
          fcur                         = ipar( 8) 
          debug                        = ipar( 9)
          gridIsImplicit               = ipar(10)
          useUpwindDissipation         = ipar(11)  ! explicit upwind dissipation
          useImplicitUpwindDissipation = ipar(12)  ! true if upwind-dissipation is on for impliciit time-stepping
          preComputeUpwindUt           = ipar(13)
          numberOfFrequencies          = ipar(14)
          adjustOmega                  = ipar(15)
          solveHelmholtz               = ipar(16)
          adjustHelmholtzForUpwinding  = ipar(17)
          fprev = mod(fcur-1+numberOfForcingFunctions,max(1,numberOfForcingFunctions))
          fnext = mod(fcur+1                         ,max(1,numberOfForcingFunctions))
     ! ** fix me ***
          timeSteppingMethod=modifiedEquationTimeStepping
     ! Set dr(:) = dx(:) for 6th-order derivatives
          if( gridType.eq.rectangular )then
              do axis=0,2
                  dr(axis)=dx(axis)
              end do
          else
              do axis=0,2
                  dx(axis)=dr(axis)
              end do
          end if  
     ! Do this for now: 
          maxDeriv=6
          uc=0
          gridIndexRange(0,0)=n1a
          gridIndexRange(1,0)=n1b
          gridIndexRange(0,1)=n2a
          gridIndexRange(1,1)=n2b
          gridIndexRange(0,2)=n3a
          gridIndexRange(1,2)=n3b    
     ! ---- Compute the coefficients in the implicit time-stepping scheme ----
          beta2=bImp(0)
          beta4=bImp(1)
          alpha2 = (1.-beta2)/2.
          alpha4 = (alpha2-beta4-1./12.)/2. 
          cImp(-1,0)=alpha2
          cImp( 0,0)= beta2
          cImp( 1,0)=alpha2
          cImp(-1,1)=alpha4
          cImp( 0,1)= beta4
          cImp( 1,1)=alpha4  
     ! ! addDissipation=.true. if we add the dissipation in the dis(i1,i2,i3,c) array
     ! !  if combineDissipationWithAdvance.ne.0 we compute the dissipation on the fly in the time step
     ! !  rather than pre-computing it in diss(i1,i2,i3,c)
     ! addDissipation = adc.gt.0.
     ! adcdt=adc*dt
          csq=cc**2
          dtsq=dt**2
          cdtsq=(cc**2)*(dt**2)
          cdt=cc*dt
     ! new: 
          cdtPow2        = cdt**2
          cdtPow4By12    = cdt**4/12.
          cdtPow6By360   = cdt**6/360. 
          cdtPow8By20160 = cdt**8/20160.  
          cdtsq12=cdtsq*cdtsq/12.  ! c^4 dt^4 /12 
     ! cdt4by360=(cdt)**4/360.  ! (c*dt)^4/360 
     ! cdt6by20160=cdt**6/(8.*7.*6.*5.*4.*3.)
          cdtSqBy12= cdtsq/12.   ! c^2*dt*2/12
          dt4by12=dtsq*dtsq/12.
          cdtdx = (cc*dt/dx(0))**2
          cdtdy = (cc*dt/dy)**2
          cdtdz = (cc*dt/dz)**2
     ! dxsqi=1./(dx(0)**2)
     ! dysqi=1./(dy**2)
     ! dzsqi=1./(dz**2)
     ! dxsq12i=1./(12.*dx(0)**2)
     ! dysq12i=1./(12.*dy**2)
     ! dzsq12i=1./(12.*dz**2)
     ! dxi4=1./(dx(0)**4)
     ! dyi4=1./(dy**4)
     ! dxdyi2=1./(dx(0)*dx(0)*dy*dy)
     ! dzi4=1./(dz**4)
     ! dxdzi2=1./(dx(0)*dx(0)*dz*dz)
     ! dydzi2=1./(dy*dy*dz*dz)
          if( option.eq.1 )then 
            useSosupDissipation = 1
          else
            useSosupDissipation = 0
          end if
          if( (.false. .or. debug.gt.1) .and. t.le.dt )then
              write(*,'("advWaveME: option=",i4," grid=",i4)') option,grid
              write(*,'("advWaveME: orderOfAccuracy=",i2," orderInTime=",i2  )') orderOfAccuracy,orderInTime
              write(*,'("advWaveME: addForcing=",i2," forcingOption=",i2)') addForcing,forcingOption
              write(*,'("advWaveME: useUpwindDissipation=",i2,"(explicit), useImplicitUpwindDissipation=",i2," (implicit)")') useUpwindDissipation,useImplicitUpwindDissipation
              write(*,'("advWaveME: useSosupDissipation=",i2,"(1= add upwind dissipation in this stage)")') useSosupDissipation
              write(*,'("advWaveME: t,dt,c,omega=",4e10.2)') t,dt,cc,omega 
              write(*,'("advWaveME: gridIsImplicit=",i2," adjustOmega=",i2," solveHelmholtz=",i2)') gridIsImplicit,adjustOmega,solveHelmholtz
              if( forcingOption.eq.helmholtzForcing )then
                  write(*,'("advWaveME: numberOfFrequencies=",i2)') numberOfFrequencies
                  write(*,'("advWaveME: frequencyArray=",(1pe12.4,1x))') (frequencyArray(freq),freq=0,numberOfFrequencies-1)
              end if
              if( gridIsImplicit.eq.1 )then
                  write(*,'("  Implicit coeff: cImp(-1:1,0) = ",3(1pe10.2,1x), "(for 2nd-order)")') cImp(-1,0),cImp(0,0),cImp(1,0)
                  write(*,'("  Implicit coeff: cImp(-1:1,1) = ",3(1pe10.2,1x), "(for 4th-order)")') cImp(-1,1),cImp(0,1),cImp(1,1)
              end if
          end if
          if( forcingOption.eq.helmholtzForcing )then
       ! --- solving the Helmholtz problem ---
              if( t.le.dt .and. debug.gt.1 )then
                  write(*,'("advWaveME: numberOfFrequencies=",i6," omega=",1pe12.4," frequencyArray(0)=",1pe12.4)') numberOfFrequencies,omega,frequencyArray(0)
              end if
              if( numberOfFrequencies.le.0 )then
                  write(*,'("advWaveME: ERROR: numberOfFrequencies=",i6," is <= 0")') numberOfFrequencies
                  stop 0123
              end if
              if( numberOfFrequencies.eq.1  .and. frequencyArray(0) .ne. omega )then
                  write(*,'("advWaveME: ERROR: frequencyArray(0)=",1pe12.4," is not equal to omega=",1pe12.4)') frequencyArray(0),omega
                  stop 1234
              end if
              if( numberOfFrequencies.gt.maxFreq )then
                  write(*,'("advWaveME: ERROR: numberOfFrequencies > maxFreq=",i6," .. FIX ME")') maxFreq
                  stop 2345
              end if
       ! if( numberOfFrequencies.gt.1 .and. gridIsImplicit.eq.1 )then
       !   write(*,'("advWave: ERROR: numberOfFrequencies > 1 and implicit time-stepping : FINISH ME")') 
       !   stop 3456  
       ! end if
              do freq=0,numberOfFrequencies-1
                  cosFreqt(freq) = cos(frequencyArray(freq)*t)
              end do
          end if
     ! write(*,'(" advWave: timeSteppingMethod=",i2)') timeSteppingMethod
          if( timeSteppingMethod.eq.defaultTimeStepping )then
            write(*,'(" advWaveME:ERROR: timeSteppingMethod=defaultTimeStepping -- this should be set")')
        ! '
            stop 83322
          end if
          if( gridIsImplicit.eq.0 )then 
       ! ------- EXPLIICT update the solution ---------
                  if( orderInTime.eq.2 )then
          ! FD24 : second-order in time and fourth-order in space
          ! FD26 : second-order in time and sixth-order in space
                        if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
                            write(*,'("advWaveME: ADVANCE dim=2 order=8 orderInTime=2, grid=rectangular... t=",e10.2)') t
                        end if
                        m=0 ! component number 
                        ec = 0 ! component number  
                                    numGhost1=3; ! should depend on the orderOfAccuracy
                                    n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                    n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                    n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                      do i3=n3a,n3b
                                      do i2=n2a,n2b
                                      do i1=n1a,n1b
                                        if( mask(i1,i2,i3).ne.0 )then
                                            d200(i1,i2,i3,0) = (u(i1+1,i2,i3,0) - 2.*u(i1,i2,i3,0) + u(i1-1,i2,i3,0))/(dx(0)**2)
                                            d020(i1,i2,i3,0) = (u(i1,i2+1,i3,0) - 2.*u(i1,i2,i3,0) + u(i1,i2-1,i3,0))/(dx(1)**2)
                                        end if ! mask .ne. 0
                                      end do
                                      end do
                                      end do
                                    numGhost1=2; 
                                    n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                    n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                    n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                      do i3=n3a,n3b
                                      do i2=n2a,n2b
                                      do i1=n1a,n1b
                                        if( mask(i1,i2,i3).ne.0 )then
                                            d400(i1,i2,i3,0) = (d200(i1+1,i2,i3,0) - 2.*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0))/(dx(0)**2)
                                            d220(i1,i2,i3,0) = (d200(i1,i2+1,i3,0) - 2.*d200(i1,i2,i3,0) + d200(i1,i2-1,i3,0))/(dx(1)**2)
                                            d040(i1,i2,i3,0) = (d020(i1,i2+1,i3,0) - 2.*d020(i1,i2,i3,0) + d020(i1,i2-1,i3,0))/(dx(1)**2)
                                        end if ! mask .ne. 0
                                      end do
                                      end do
                                      end do
                                    numGhost1=1; 
                                    n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                    n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                    n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                      do i3=n3a,n3b
                                      do i2=n2a,n2b
                                      do i1=n1a,n1b
                                        if( mask(i1,i2,i3).ne.0 )then
                                            d600(i1,i2,i3,0) = (d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0))/(dx(0)**2)
                                            d060(i1,i2,i3,0) = (d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0))/(dx(1)**2)
                                            d420(i1,i2,i3,0) = (d220(i1+1,i2,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0))/(dx(0)**2)
                                            d240(i1,i2,i3,0) = (d220(i1,i2+1,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1,i2-1,i3,0))/(dx(1)**2)
                                        end if ! mask .ne. 0
                                      end do
                                      end do
                                      end do
                                    numGhost1=0; 
                                    n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                    n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                    if( nd.eq.2 )then
                                        n3a=gridIndexRange(0,2); n3b=gridIndexRange(1,2);
                                    else
                                        n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                    end if
                  ! ------------ MAIN LOOP, 8 EIGHTH ----------------
                                    fv(m)=0.
                                      do i3=n3a,n3b
                                      do i2=n2a,n2b
                                      do i1=n1a,n1b
                                        if( mask(i1,i2,i3).ne.0 )then
                                            d600i = d600(i1,i2,i3,0)  !  (d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0))/(dx(0)**2)
                                            d060i = d060(i1,i2,i3,0)  !  (d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0))/(dx(1)**2)
                                            d240i = d240(i1,i2,i3,0)    
                                            d420i = d420(i1,i2,i3,0)    
                                            d800i = (d600(i1+1,i2,i3,0) - 2.*d600(i1,i2,i3,0) + d600(i1-1,i2,i3,0))/(dx(0)**2)
                                            d080i = (d060(i1,i2+1,i3,0) - 2.*d060(i1,i2,i3,0) + d060(i1,i2-1,i3,0))/(dx(1)**2)
                                            d620i = (d420(i1+1,i2,i3,0) - 2.*d420(i1,i2,i3,0) + d420(i1-1,i2,i3,0))/(dx(0)**2)  
                                            d260i = (d240(i1,i2+1,i3,0) - 2.*d240(i1,i2,i3,0) + d240(i1,i2-1,i3,0))/(dx(1)**2)
                                            d440i = (d240(i1+1,i2,i3,0) - 2.*d240(i1,i2,i3,0) + d240(i1-1,i2,i3,0))/(dx(0)**2)  
                      ! ! ----- for odd derivatives ----
                      ! d100i = (u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0))/(2.*dx(0))   ! D0x 
                      ! d010i = (u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0))/(2.*dx(1))   ! D0y 
                      ! d300i = (d200(i1+1,i2,i3,0) - d200(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x D_x D-x 
                      ! d030i = (d020(i1,i2+1,i3,0) - d020(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y D+y D-y 
                      ! d120i = (d020(i1+1,i2,i3,0) - d020(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x D_y D-y 
                      ! d210i = (d200(i1,i2+1,i3,0) - d200(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y D+x D-x       
                      ! d500i = (d400(i1+1,i2,i3,0) - d400(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_x D-x)^2 
                      ! d050i = (d040(i1,i2+1,i3,0) - d040(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+y D-y)^2 
                      ! d140i = (d040(i1+1,i2,i3,0) - d040(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_y D-y)^2 
                      ! d410i = (d400(i1,i2+1,i3,0) - d400(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)^2
                      ! d320i = (d220(i1+1,i2,i3,0) - d220(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D+x D-x)(D_y D-y)
                      ! d230i = (d220(i1,i2+1,i3,0) - d220(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)
                      ! d700i = (d600(i1+1,i2,i3,0) - d600(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_x D-x)^3 
                      ! d070i = (d060(i1,i2+1,i3,0) - d060(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+y D-y)^3 
                      ! d160i = (d060(i1+1,i2,i3,0) - d060(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_y D-y)^3 
                      ! d610i = (d600(i1,i2+1,i3,0) - d600(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)^3
                      ! d520i = (d420(i1+1,i2,i3,0) - d420(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D+x D-x)(D_y D-y)
                      ! d250i = (d240(i1,i2+1,i3,0) - d240(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)      
                      ! d340i = (d240(i1+1,i2,i3,0) - d240(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D+x D-x)^2(D_y D-y)
                      ! d430i = (d420(i1,i2+1,i3,0) - d420(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)^2
                      ! d900i = (d800(i1+1,i2,i3,0) - d800(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_x D-x)^4 
                      ! d090i = (d080(i1,i2+1,i3,0) - d080(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+y D-y)^4 
                      ! ! --- for mixed derivatives
                      ! d110i = (u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))                ! D0x D0y
                      ! d310i = (d200(i1+1,i2+1,i3,0) - d200(i1-1,i2+1,i3,0) - d200(i1+1,i2-1,i3,0) + d200(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) 
                      ! d130i = (d020(i1+1,i2+1,i3,0) - d020(i1-1,i2+1,i3,0) - d020(i1+1,i2-1,i3,0) + d020(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)          
                      ! d510i = (d400(i1+1,i2+1,i3,0) - d400(i1-1,i2+1,i3,0) - d400(i1+1,i2-1,i3,0) + d400(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2
                      ! d150i = (d040(i1+1,i2+1,i3,0) - d040(i1-1,i2+1,i3,0) - d040(i1+1,i2-1,i3,0) + d040(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^2     
                      ! d330i = (d220(i1+1,i2+1,i3,0) - d220(i1-1,i2+1,i3,0) - d220(i1+1,i2-1,i3,0) + d220(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)       
                      ! d710i = (d600(i1+1,i2+1,i3,0) - d600(i1-1,i2+1,i3,0) - d600(i1+1,i2-1,i3,0) + d600(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^3
                      ! d170i = (d060(i1+1,i2+1,i3,0) - d060(i1-1,i2+1,i3,0) - d060(i1+1,i2-1,i3,0) + d060(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^3  
                      ! d530i = (d420(i1+1,i2+1,i3,0) - d420(i1-1,i2+1,i3,0) - d420(i1+1,i2-1,i3,0) + d420(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2 (D+yD-y)         
                      ! d350i = (d240(i1+1,i2+1,i3,0) - d240(i1-1,i2+1,i3,0) - d240(i1+1,i2-1,i3,0) + d240(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)^2         
                      ! ! ---- first derivatives ----
                      ! !   1/6, 1/30, 1/140, 1/630 
                      ! !   ux = D0 ( 1 - (1/6)*dx^2(D+D-) + (1/30)*dx^4 (D+D-)^2 -(1/140)*(D+D-)^3 + ... )
                      ! ux = d100i - (dx(0)**2/6.)*d300i + (dx(0)**4/30.)*d500i - (dx(0)**6/140.)*d700i 
                      ! uy = d010i - (dx(1)**2/6.)*d030i + (dx(1)**4/30.)*d050i - (dx(1)**6/140.)*d070i  
                      ! ---- 2nd derivatives ----
                      ! 1/12 , 1/90, 1/560 
                      !   uxx = D+D- ( 1 - (1/12)*dx^2(D+D-) + (1/90)*dx^4 (D+D-)^2 -(1/560)*(D+D-)^3 + ... )
                                            uxx = d200(i1,i2,i3,0) - (dx(0)**2/12.)*d400(i1,i2,i3,0) + (dx(0)**4/90.)*d600i - (dx(0)**6/560.)*d800i
                                            uyy = d020(i1,i2,i3,0) - (dx(1)**2/12.)*d040(i1,i2,i3,0) + (dx(1)**4/90.)*d060i - (dx(1)**6/560.)*d080i 
                      ! uxy = D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                      !       D0y*( 1 - dy^2/6 D+yD-y + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+yD-y)^3 ) u 
                      ! uxy = d110i - (dx(0)**2/6.)*d310i - (dx(1)**2/6.)*d130i + (dx(0)**4/30.)*d510i + (dx(1)**4/30.)*d150i + (dx(0)**2 * dx(1)**2/36. )*d330i      !       - ( (1./140.)*dx(0)**6 )*d710i - ( (1./140.)*dx(1)**6 )*d170i - (dx(0)**4*dx(1)**2/(6.*30.))*d530i - (dx(0)**2*dx(1)**4/(6.*30.))*d350i 
                      ! ! ---- third derivatives ----
                      ! !   1/4, 7/120, 41/3024
                      ! !    uxxxx = D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ... 
                      ! uxxx = d300i - (dx(0)**2/4.)*d500i + (dx(0)**4 *7./120.)*d700i
                      ! uyyy = d030i - (dx(1)**2/4.)*d050i + (dx(1)**4 *7./120.)*d070i
                      ! ! uxxy = D+xD-x*( 1 - dx^2/12 D+xD-x + dx^4/90*(D+xD-x)^2 - (1/560)*dx^6*(D+xD-x)^3 ) X
                      ! !           D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                      ! uxxy = d210i - (dx(0)**2/12.)*d410i - (dx(1)**2/6.)*d230i + (dx(0)**4/90.)*d610i + (dx(1)**4/30.)*d250i + (dx(0)**2 * dx(1)**2/72. )*d430i
                      ! ! uxyy =    D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                      ! !        D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                      ! uxyy = d120i - (dx(1)**2/12.)*d140i - (dx(0)**2/6.)*d320i + (dx(1)**4/90.)*d160i + (dx(0)**4/30.)*d520i + (dx(0)**2 * dx(1)**2/72. )*d340i      
                      ! ----- fourth derivs ------
                      !  1/6, 7/240, 41/7560
                      !   uxxxx = (D+xD-x)^2*( 1 - (1/6)*D+xD-x + (7/240)*(D+xD-x)^2 + ... )
                                            uxxxx = d400(i1,i2,i3,0) - (dx(0)**2/6.)*d600(i1,i2,i3,0) + (dx(0)**4 * (7./240))*d800i
                                            uyyyy = d040(i1,i2,i3,0) - (dx(1)**2/6.)*d060(i1,i2,i3,0) + (dx(1)**4 * (7./240))*d080i
                      ! ! uxyyy = D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                      ! !         D0y D+yD-y*( 1 - (1/4)*dy^2 * D+yD-y + (7/120)*dy^4 (D+yD-y)^2 + ... 
                      ! uxyyy = d130i - (dx(0)**2/6.)*d330i - (dx(1)**2/4.)*d150i + (dx(0)**4/30.)*d530i + (dx(1)**4*7./120.)*d170i + (dx(0)**2 * dx(1)**2/24. )*d350i 
                      ! ! uxxxy
                      ! uxxxy = d310i - (dx(1)**2/6.)*d330i - (dx(0)**2/4.)*d510i + (dx(1)**4/30.)*d350i + (dx(0)**4*7./120.)*d710i + (dx(0)**2 * dx(1)**2/24. )*d530i       
                      ! uxxyy = D+xD-x*( 1 - dx^2/12 D+xD-x + dx^4/90*(D+xD-x)^2 - (1/560)*dx^6*(D+xD-x)^3 ) X 
                      !         D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                            uxxyy = d220(i1,i2,i3,0) - (dx(0)**2/12.)*d420(i1,i2,i3,0) - (dx(1)**2/12.)*d240(i1,i2,i3,0) + (dx(0)**4 * (1./90))*d620i + (dx(1)**4 * (1./90))*d260i + (dx(0)**2 * dx(1)**2 * 1./144.)*d440i 
                      ! ! ---- fixth derivatives ----
                      ! !   1/3, 13/144, 139/6048
                      ! !  uxxxxx = D0x(D+xD-x)^2 *( 1 - (1/3)*dx^2 D+xD-x + ... )
                      ! uxxxxx = d500i - (dx(0)**2/3.)*d700i
                      ! uyyyyy = d050i - (dx(1)**2/3.)*d070i 
                      ! ! uxxxxy = (D+xD-x)^2 *( 1 - dx^2/6 D+xD-x + ...) X
                      ! !                  D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                      ! uxxxxy = d410i - (dx(0)**2/6.)*d610i - (dx(1)**2/6.)*d430i
                      ! uxyyyy = d140i - (dx(1)**2/6.)*d160i - (dx(0)**2/6.)*d340i    
                      ! ! uxxxyy =  D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ...
                      ! !               D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                      ! uxxxyy = d320i - (dx(0)**2/4.)*d520i - (dx(1)**2/12.)*d340i
                      ! uxxyyy = d230i - (dx(1)**2/4.)*d250i - (dx(0)**2/12.)*d430i
                      ! ---- sixth derivs----
                      !  1/4 , 13/240 
                                            uxxxxxx = d600i - (dx(0)**2/4.)*d800i
                                            uyyyyyy = d060i - (dx(1)**2/4.)*d080i  
                      ! ! uxxxxxy = D0x(D+xD-x)^2 *( 1 - (1/3)*dx^2 D+xD-x + ... ) X 
                      ! !                      D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                      ! uxxxxxy = d510i - (dx(0)**2/3.)*d710i - (dx(1)**2/6.)*d530i
                      ! uxyyyyy = d150i - (dx(1)**2/3.)*d170i - (dx(0)**2/6.)*d350i
                      ! uxxxxyy = (D+xD-x)^2*( 1 - (1/6)*D+xD-x + (7/240)*(D+xD-x)^2 + ... )
                      !               D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                            uxxxxyy = d420i - (dx(0)**2/6.)*d620i - (dx(1)**2/12.)*d440i
                                            uxxyyyy = d240i - (dx(1)**2/6.)*d260i - (dx(0)**2/12.)*d440i
                      ! ! uxxxyyy = D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ... ) X 
                      ! !           D0y D+yD-y*( 1 - (1/4)*dy^2 * D+yD-y + (7/120)*dy^4 (D+yD-y)^2 + ...  ) u 
                      ! uxxxyyy = d330i - (dx(0)**2/4.)*d530i - (dx(1)**2/4.)*d350i
                      ! ! ---- seventh derivatives ----
                      ! uxxxxxxx = d700i
                      ! uxxxxxxy = d610i
                      ! uxxxxxyy = d520i
                      ! uxxxxyyy = d430i
                      ! uxxxyyyy = d340i
                      ! uxxyyyyy = d250i
                      ! uxyyyyyy = d160i
                      ! uyyyyyyy = d070i 
                      ! ---- eighth derivatives ----
                                            uxxxxxxxx = d800i
                      ! uxxxxxxxy = d710i
                                            uxxxxxxyy = d620i
                      ! uxxxxxyyy = d530i
                                            uxxxxyyyy = d440i
                      ! uxxxyyyyy = d350i
                                            uxxyyyyyy = d260i
                      ! uxyyyyyyy = d170i
                                            uyyyyyyyy = d080i                  
                                                if( forcingOption.eq.twilightZoneForcing )then
                                                            call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                                            call ogDeriv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evtt(m) )
                                                            call ogDeriv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxx(m) )
                                                            call ogDeriv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evyy(m) )
                                                        fv(m) = evtt(m) - csq*( evxx(m) + evyy(m) )
                                              else if( forcingOption.eq.helmholtzForcing )then
                          ! forcing for solving the Helmholtz equation   
                          ! NOTE: change sign of forcing since for Helholtz we want to solve
                          !      ( omega^2 I + c^2 Delta) w = f 
                          ! fv(m) = -f(i1,i2,i3,0)*coswt  
                                                    fv(m)=0.
                                                    do freq=0,numberOfFrequencies-1 
                                                        omega = frequencyArray(freq)
                                                        coswt = cosFreqt(freq)    
                             ! if( i1.eq.2 .and. i2.eq.2 )then 
                             !   write(*,'(" adv: forcing f(i1,i2,i3)=",1pe12.4," coswt=",1pe12.4," t=",1pe12.4," omega=",1pe12.4)') f(i1,i2,i3,0),coswt,t,omega
                             ! end if
                             ! fv(m) = -f(i1,i2,i3,0)*coswt  
                                                          fv(m) = fv(m) - f(i1,i2,i3,freq)*coswt
                                                    end do ! do freq  
                                              else if( addForcing.ne.0 )then  
                                                    fv(m) = f(i1,i2,i3,0)
                                              end if
                      ! Here is the ME update 
                                            un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx + uyy ) + cdtPow4By12*( uxxxx + uyyyy + 2.*uxxyy )  + cdtPow6By360*( uxxxxxx + uyyyyyy +3.*( uxxxxyy + uxxyyyy) ) + cdtPow8By20160*( uxxxxxxxx + uyyyyyyyy + 4.*( uxxxxxxyy +uxxyyyyyy ) +6.*uxxxxyyyy ) + dtSq*fv(m)      
                    ! if( i1.eq.5 .and. i2.eq.6 )then
                    !   call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t+dt,uc,ue )
                    !   write(*,'("ME8: (i1,i2)=(",2i3,") u=",1pe12.4," ue=",1pe12.4," err=",1pe8.2)') i1,i2,un(i1,i2,i3,0),ue,abs(un(i1,i2,i3,0)-ue)
                    ! end if
                                        end if ! mask .ne. 0
                                      end do
                                      end do
                                      end do
            ! #If 8 == 6 || 8 == 8 
            !   evalDerivativesRectangular()
            !   write(*,*) ' Stop here for now'
            !   stop 666
            ! #End
          !   ! --- TAYLOR TIME-STEPPING --- 
          !   m=0 ! component number 
          !   ec = 0 ! component number
          !   ! #If "rectangular" eq "curvilinear"
          !   !   #If "8" eq "4" && "2" eq "4"
          !   !     computeLaplacianOrder2(2)
          !   !   #End
          !   ! #End
          !   if( forcingOption.eq.helmholtzForcing )then
          !     coswt = cos(omega*t)
          !   end if 
          !   fv(m)=0.
          !   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
          !     getForcing(2,8,2,rectangular) 
          !    #If "8" eq "2"
          !      ! --- SECOND 8 ---
          !      #If "2" eq "2"
          !        ! --- TWO DIMENSIONS ---
          !        #If "rectangular" eq "rectangular"
          !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx22r(i1,i2,i3,0) + uyy22r(i1,i2,i3,0) ) + dtSq*fv(m)
          !         ! write(*,'(" adv: i1,i2=",2i4," un,u,um=",3e12.2," cdtSq,fv=",2e12.2)') i1,i2,un(i1,i2,i3,m),u(i1,i2,i3,m),um(i1,i2,i3,m),cdtSq,fv(m)
          !        #Else
          !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx22(i1,i2,i3,0)  + uyy22(i1,i2,i3,0) ) + dtSq*fv(m)
          !        #End
          !      #Else
          !        ! --- THREE DIMENSIONS ---
          !        #If "rectangular" eq "rectangular"
          !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx23r(i1,i2,i3,0) + uyy23r(i1,i2,i3,0) + uzz23r(i1,i2,i3,0) ) + dtSq*fv(m)
          !        #Else
          !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx23(i1,i2,i3,0)  + uyy23(i1,i2,i3,0)  + uzz23(i1,i2,i3,0)  ) + dtSq*fv(m)
          !        #End
          !      #End
          !    #Elif "8" eq "4"
          !      ! --- -FOURTH 8 ---
          !      #If "2" eq "2"
          !        ! --- FOUTH-8 TWO DIMENSIONS ---
          !        #If "2" eq "4"
          !          ! orderInSpace=4 and orderInTime=4 
          !          #If "rectangular" eq "rectangular"
          !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + cdtsq12*lap2d2pow2(i1,i2,i3,m) + dtSq*fv(m)
          !          #Else
          !            ! v is assumed to hold Lap(u) to 2nd-order
          !            ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
          !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !                                                            + cdtsq12*( vxx22(i1,i2,i3,m) + vyy22(i1,i2,i3,m) ) + dtSq*fv(m)
          !          #End
          !        #Else
          !          ! orderInSpace==4 and orderInTime==2                                                   
          !          #If "rectangular" eq "rectangular"
          !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + dtSq*fv(m)
          !          #Else
          !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End
          !        #End                                                          
          !      #Else
          !        ! --- FOURTH-8 THREE DIMENSIONS ---
          !        #If "2" eq "4"
          !          ! orderInSpace=4 and orderInTime=4 
          !          #If "rectangular" eq "rectangular"
          !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + cdtsq12*lap3d2pow2(i1,i2,i3,m) + dtSq*fv(m)
          !          #Else
          !            ! v is assumed to hold Lap(u) to 2nd-order
          !            ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
          !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !                                                            + cdtsq12*( vxx23(i1,i2,i3,m) + vyy23(i1,i2,i3,m) + vzz23(i1,i2,i3,m) ) + dtSq*fv(m)
          !          #End
          !        #Else
          !          ! orderInSpace==4 and orderInTime==2                                                   
          !          #If "rectangular" eq "rectangular"
          !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + dtSq*fv(m)
          !          #Else
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End
          !        #End     
          !      #End
          !    #Elif "8" eq "6"
          !      ! ---- SIXTH 8 ---
          !      #If "2" eq "2"
          !        ! --- SIXTH-8 TWO DIMENSIONS ---
          !        #If "2" eq "2"
          !          ! orderInSpace==6 and orderInTime==2                                                   
          !          #If "rectangular" eq "rectangular"
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx62r(i1,i2,i3,m) + uyy62r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #Else
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx62(i1,i2,i3,m)  +  uyy62(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End       
          !        #Else
          !          stop 6666                                                           
          ! !          ! ---- MODIFIED EQUATION 8=6 2D -----
          ! !          getSixthDerivatives2d(8,rectangular,evalMetrics,i1,i2,i3)
          ! !          if( i1.eq.21 .and. i2.eq.11 )then
          ! !            OGDERIV2D( 0,2,0,0,i1,i2,i3,t, ec, evxx(0) )
          ! !            OGDERIV2D( 0,0,2,0,i1,i2,i3,t, ec, evyy(0) )
          ! !            OGDERIV2D( 0,4,0,0,i1,i2,i3,t, ec, evxxxx(0) )
          ! !            OGDERIV2D( 0,0,4,0,i1,i2,i3,t, ec, evyyyy(0) )
          ! !            OGDERIV2D( 0,2,2,0,i1,i2,i3,t, ec, evxxyy(0) )
          ! !            write(*,'("i1,i2=",2i4)') i1,i2
          ! !            write(*,'(" uxx  =",1pe12.4," true=",1pe12.4," err=",1pe8.2)') uxx  ,evxx(0)  ,abs(uxx  -evxx(0)  )
          ! !            write(*,'(" uyy  =",1pe12.4," true=",1pe12.4," err=",1pe8.2)') uyy  ,evyy(0)  ,abs(uyy  -evyy(0)  )
          ! !            write(*,'(" uxxxx=",1pe12.4," true=",1pe12.4," err=",1pe8.2)') uxxxx,evxxxx(0),abs(uxxxx-evxxxx(0))
          ! !            write(*,'(" uxxyy=",1pe12.4," true=",1pe12.4," err=",1pe8.2)') uxxyy,evxxyy(0),abs(uxxyy-evxxyy(0))
          ! !            write(*,'(" uyyyy=",1pe12.4," true=",1pe12.4," err=",1pe8.2)') uyyyy,evyyyy(0),abs(uyyyy-evyyyy(0))
          ! ! ! uxxxx    = rxi**4*urrrr+4.*rxi**3*sxi*urrrs+6.*rxi**2*sxi**2*urrss+4.*rxi*sxi**3*ursss+sxi**4*ussss+6.*rxi**2*rxx*urrr+(7.*sxi*rxi*rxx+sxx*rxi**2+rxi*(3.*rxi*sxx+3.*rxx*sxi)+rxi*(2.*rxi*sxx+2.*rxx*sxi))*urrs+(sxi*(3.*rxi*sxx+3.*rxx*sxi)+7.*rxi*sxx*sxi+rxx*sxi**2+sxi*(2.*rxi*sxx+2.*rxx*sxi))*urss+6.*sxi**2*sxx*usss+(4.*rxi*rxxx+3.*rxx**2)*urr+(4.*rxi*sxxx+6.*rxx*sxx+4.*rxxx*sxi)*urs+(4.*sxi*sxxx+3.*sxx**2)*uss+rxxxx*ur+sxxxx*us
          ! !            stop 1111
          ! !          end if
          ! !          ! cdtPow4By12  = (c*dt)^4/12 
          ! !          ! cdtPow6By360 = (c*dt)^6/360
          ! !          un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) ! !                          + cdtsq*( uxx + uyy ) ! !                          + cdtPow4By12*( uxxxx + uyyyy + 2.*uxxyy )  ! !                          + cdtPow6By360*( uxxxxxx + uyyyyyy + 3.*(uxxxxyy + uxxyyyy) ) ! !                          + dtSq*fv(m)
          !       #End                                                          
          !     #Else
          !        ! --- SIXTH-8 THREE DIMENSIONS ---
          !        #If "2" eq "2
          !          ! orderInSpace==6 and orderInTime==2                                                   
          !          #If "rectangular" eq "rectangular"
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx63r(i1,i2,i3,m) + uyy63r(i1,i2,i3,m) + uzz63r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #Else
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx63(i1,i2,i3,m)  + uyy63(i1,i2,i3,m)  + uzz63(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End       
          !        #Else
          !          ! MODIFIED EQUATION 8=6 3D
          !          ! Turn off for now: 
          !          ! getSixthDerivatives3d(8,rectangular,evalMetrics,i1,i2,i3)
          !          write(*,'("advWave: order=8, orderInTime=2 FINISH ME")')
          !          stop 7777
          !          ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) !          !                 + cdtsq*( uxx + uyy +uzz ) !          !                 + cdtPow4By12*( uxxxx + uyyyy + uzzzz + 2.*( uxxyy +uxxzz + uyyzz ) )  !          !                 + cdtPow6By360*( uxxxxxx +  uyyyyyy + uzzzzzz + 3.*(uxxxxyy + uxxyyyy + uxxxxzz + uyyyyzz + uxxzzzz + uyyzzzz ) + 6.*uxxyyzz ) !          !                 + dtSq*fv(m)
          !        #End     
          !      #End
          !    #Elif "8" eq "8"
          !      ! ---- EIGTH 8 ---
          !      #If "2" eq "2"
          !        ! --- EIGTH-8 TWO DIMENSIONS ---
          !        #If "2" eq "2"
          !          ! orderInSpace==8 and orderInTime==2                                                   
          !          #If "rectangular" eq "rectangular"
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx82r(i1,i2,i3,m) + uyy82r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #Else
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx82(i1,i2,i3,m)  +  uyy82(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End       
          !        #Else
          !          write(*,'("advWave: order=8, orderInTime=2 FINISH ME")')
          !          stop 7777
          !          ! ! orderInSpace=4 and orderInTime=4 
          !          ! #If "rectangular" eq "rectangular"
          !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + cdtsq12*lap2d2pow2(i1,i2,i3,m) + dtSq*fv(m)
          !          ! #Else
          !          !   ! v is assumed to hold Lap(u) to 2nd-order
          !          !   ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
          !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !          !                                                   + cdtsq12*( vxx22(i1,i2,i3,m) + vyy22(i1,i2,i3,m) ) + dtSq*fv(m)
          !          ! #End
          !        #End                                                          
          !     #Else
          !        ! --- EIGTH-8 THREE DIMENSIONS ---
          !        #If "2" eq "2
          !          ! orderInSpace==8 and orderInTime==2                                                   
          !          #If "rectangular" eq "rectangular"
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx83r(i1,i2,i3,m) + uyy83r(i1,i2,i3,m) + uzz83r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #Else
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx83(i1,i2,i3,m)  + uyy83(i1,i2,i3,m)  + uzz83(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End       
          !        #Else
          !          write(*,'("advWave: order=8, orderInTime=2 FINISH ME")')
          !          stop 7777
          !          ! #If "rectangular" eq "rectangular"
          !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + cdtsq12*lap3d2pow2(i1,i2,i3,m) + dtSq*fv(m)
          !          ! #Else
          !          !   ! v is assumed to hold Lap(u) to 2nd-order
          !          !   ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
          !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !          !                                                   + cdtsq12*( vxx23(i1,i2,i3,m) + vyy23(i1,i2,i3,m) + vzz23(i1,i2,i3,m) ) + dtSq*fv(m)
          !          ! #End
          !        #End     
          !      #End
          !    #Else
          !      write(*,'("advWave: UNKNOWN order=8")')
          !      stop 7777
          !    #End         
          !    ! write(*,'("i1,i2=",2i3," u-ue=",e10.2)') i1,i2,u(i1,i2,i3,m)-ev(m)
          !    ! write(*,'(" uxx-uxxe =",e10.2)') uxx22r(i1,i2,i3,0)-evxx(m)
          !    ! OGDERIV2D( 0,0,0,0,i1,i2,i3,t+dt, ec, ev(m)  )
          !    ! write(*,'(" un-ue=",e10.2)') un(i1,i2,i3,m)-ev(m)
          !   endLoopsMask()
                  else
                          if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
                              write(*,'("advWaveME: ADVANCE dim=2 order=8 orderInTime=8, grid=rectangular... t=",e10.2)') t
                          end if
                          m=0 ! component number 
                          ec = 0 ! component number  
                                      numGhost1=3; ! should depend on the orderOfAccuracy
                                      n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                      n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                      n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                        do i3=n3a,n3b
                                        do i2=n2a,n2b
                                        do i1=n1a,n1b
                                          if( mask(i1,i2,i3).ne.0 )then
                                              d200(i1,i2,i3,0) = (u(i1+1,i2,i3,0) - 2.*u(i1,i2,i3,0) + u(i1-1,i2,i3,0))/(dx(0)**2)
                                              d020(i1,i2,i3,0) = (u(i1,i2+1,i3,0) - 2.*u(i1,i2,i3,0) + u(i1,i2-1,i3,0))/(dx(1)**2)
                                          end if ! mask .ne. 0
                                        end do
                                        end do
                                        end do
                                      numGhost1=2; 
                                      n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                      n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                      n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                        do i3=n3a,n3b
                                        do i2=n2a,n2b
                                        do i1=n1a,n1b
                                          if( mask(i1,i2,i3).ne.0 )then
                                              d400(i1,i2,i3,0) = (d200(i1+1,i2,i3,0) - 2.*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0))/(dx(0)**2)
                                              d220(i1,i2,i3,0) = (d200(i1,i2+1,i3,0) - 2.*d200(i1,i2,i3,0) + d200(i1,i2-1,i3,0))/(dx(1)**2)
                                              d040(i1,i2,i3,0) = (d020(i1,i2+1,i3,0) - 2.*d020(i1,i2,i3,0) + d020(i1,i2-1,i3,0))/(dx(1)**2)
                                          end if ! mask .ne. 0
                                        end do
                                        end do
                                        end do
                                      numGhost1=1; 
                                      n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                      n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                      n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                        do i3=n3a,n3b
                                        do i2=n2a,n2b
                                        do i1=n1a,n1b
                                          if( mask(i1,i2,i3).ne.0 )then
                                              d600(i1,i2,i3,0) = (d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0))/(dx(0)**2)
                                              d060(i1,i2,i3,0) = (d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0))/(dx(1)**2)
                                              d420(i1,i2,i3,0) = (d220(i1+1,i2,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0))/(dx(0)**2)
                                              d240(i1,i2,i3,0) = (d220(i1,i2+1,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1,i2-1,i3,0))/(dx(1)**2)
                                          end if ! mask .ne. 0
                                        end do
                                        end do
                                        end do
                                      numGhost1=0; 
                                      n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                      n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                      if( nd.eq.2 )then
                                          n3a=gridIndexRange(0,2); n3b=gridIndexRange(1,2);
                                      else
                                          n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                      end if
                   ! ------------ MAIN LOOP, 8 EIGHTH ----------------
                                      fv(m)=0.
                                        do i3=n3a,n3b
                                        do i2=n2a,n2b
                                        do i1=n1a,n1b
                                          if( mask(i1,i2,i3).ne.0 )then
                                              d600i = d600(i1,i2,i3,0)  !  (d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0))/(dx(0)**2)
                                              d060i = d060(i1,i2,i3,0)  !  (d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0))/(dx(1)**2)
                                              d240i = d240(i1,i2,i3,0)    
                                              d420i = d420(i1,i2,i3,0)    
                                              d800i = (d600(i1+1,i2,i3,0) - 2.*d600(i1,i2,i3,0) + d600(i1-1,i2,i3,0))/(dx(0)**2)
                                              d080i = (d060(i1,i2+1,i3,0) - 2.*d060(i1,i2,i3,0) + d060(i1,i2-1,i3,0))/(dx(1)**2)
                                              d620i = (d420(i1+1,i2,i3,0) - 2.*d420(i1,i2,i3,0) + d420(i1-1,i2,i3,0))/(dx(0)**2)  
                                              d260i = (d240(i1,i2+1,i3,0) - 2.*d240(i1,i2,i3,0) + d240(i1,i2-1,i3,0))/(dx(1)**2)
                                              d440i = (d240(i1+1,i2,i3,0) - 2.*d240(i1,i2,i3,0) + d240(i1-1,i2,i3,0))/(dx(0)**2)  
                       ! ! ----- for odd derivatives ----
                       ! d100i = (u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0))/(2.*dx(0))   ! D0x 
                       ! d010i = (u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0))/(2.*dx(1))   ! D0y 
                       ! d300i = (d200(i1+1,i2,i3,0) - d200(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x D_x D-x 
                       ! d030i = (d020(i1,i2+1,i3,0) - d020(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y D+y D-y 
                       ! d120i = (d020(i1+1,i2,i3,0) - d020(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x D_y D-y 
                       ! d210i = (d200(i1,i2+1,i3,0) - d200(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y D+x D-x       
                       ! d500i = (d400(i1+1,i2,i3,0) - d400(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_x D-x)^2 
                       ! d050i = (d040(i1,i2+1,i3,0) - d040(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+y D-y)^2 
                       ! d140i = (d040(i1+1,i2,i3,0) - d040(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_y D-y)^2 
                       ! d410i = (d400(i1,i2+1,i3,0) - d400(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)^2
                       ! d320i = (d220(i1+1,i2,i3,0) - d220(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D+x D-x)(D_y D-y)
                       ! d230i = (d220(i1,i2+1,i3,0) - d220(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)
                       ! d700i = (d600(i1+1,i2,i3,0) - d600(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_x D-x)^3 
                       ! d070i = (d060(i1,i2+1,i3,0) - d060(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+y D-y)^3 
                       ! d160i = (d060(i1+1,i2,i3,0) - d060(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_y D-y)^3 
                       ! d610i = (d600(i1,i2+1,i3,0) - d600(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)^3
                       ! d520i = (d420(i1+1,i2,i3,0) - d420(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D+x D-x)(D_y D-y)
                       ! d250i = (d240(i1,i2+1,i3,0) - d240(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)      
                       ! d340i = (d240(i1+1,i2,i3,0) - d240(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D+x D-x)^2(D_y D-y)
                       ! d430i = (d420(i1,i2+1,i3,0) - d420(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)^2
                       ! d900i = (d800(i1+1,i2,i3,0) - d800(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_x D-x)^4 
                       ! d090i = (d080(i1,i2+1,i3,0) - d080(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+y D-y)^4 
                       ! ! --- for mixed derivatives
                       ! d110i = (u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))                ! D0x D0y
                       ! d310i = (d200(i1+1,i2+1,i3,0) - d200(i1-1,i2+1,i3,0) - d200(i1+1,i2-1,i3,0) + d200(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) 
                       ! d130i = (d020(i1+1,i2+1,i3,0) - d020(i1-1,i2+1,i3,0) - d020(i1+1,i2-1,i3,0) + d020(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)          
                       ! d510i = (d400(i1+1,i2+1,i3,0) - d400(i1-1,i2+1,i3,0) - d400(i1+1,i2-1,i3,0) + d400(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2
                       ! d150i = (d040(i1+1,i2+1,i3,0) - d040(i1-1,i2+1,i3,0) - d040(i1+1,i2-1,i3,0) + d040(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^2     
                       ! d330i = (d220(i1+1,i2+1,i3,0) - d220(i1-1,i2+1,i3,0) - d220(i1+1,i2-1,i3,0) + d220(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)       
                       ! d710i = (d600(i1+1,i2+1,i3,0) - d600(i1-1,i2+1,i3,0) - d600(i1+1,i2-1,i3,0) + d600(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^3
                       ! d170i = (d060(i1+1,i2+1,i3,0) - d060(i1-1,i2+1,i3,0) - d060(i1+1,i2-1,i3,0) + d060(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^3  
                       ! d530i = (d420(i1+1,i2+1,i3,0) - d420(i1-1,i2+1,i3,0) - d420(i1+1,i2-1,i3,0) + d420(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2 (D+yD-y)         
                       ! d350i = (d240(i1+1,i2+1,i3,0) - d240(i1-1,i2+1,i3,0) - d240(i1+1,i2-1,i3,0) + d240(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)^2         
                       ! ! ---- first derivatives ----
                       ! !   1/6, 1/30, 1/140, 1/630 
                       ! !   ux = D0 ( 1 - (1/6)*dx^2(D+D-) + (1/30)*dx^4 (D+D-)^2 -(1/140)*(D+D-)^3 + ... )
                       ! ux = d100i - (dx(0)**2/6.)*d300i + (dx(0)**4/30.)*d500i - (dx(0)**6/140.)*d700i 
                       ! uy = d010i - (dx(1)**2/6.)*d030i + (dx(1)**4/30.)*d050i - (dx(1)**6/140.)*d070i  
                       ! ---- 2nd derivatives ----
                       ! 1/12 , 1/90, 1/560 
                       !   uxx = D+D- ( 1 - (1/12)*dx^2(D+D-) + (1/90)*dx^4 (D+D-)^2 -(1/560)*(D+D-)^3 + ... )
                                              uxx = d200(i1,i2,i3,0) - (dx(0)**2/12.)*d400(i1,i2,i3,0) + (dx(0)**4/90.)*d600i - (dx(0)**6/560.)*d800i
                                              uyy = d020(i1,i2,i3,0) - (dx(1)**2/12.)*d040(i1,i2,i3,0) + (dx(1)**4/90.)*d060i - (dx(1)**6/560.)*d080i 
                       ! uxy = D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                       !       D0y*( 1 - dy^2/6 D+yD-y + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+yD-y)^3 ) u 
                       ! uxy = d110i - (dx(0)**2/6.)*d310i - (dx(1)**2/6.)*d130i + (dx(0)**4/30.)*d510i + (dx(1)**4/30.)*d150i + (dx(0)**2 * dx(1)**2/36. )*d330i      !       - ( (1./140.)*dx(0)**6 )*d710i - ( (1./140.)*dx(1)**6 )*d170i - (dx(0)**4*dx(1)**2/(6.*30.))*d530i - (dx(0)**2*dx(1)**4/(6.*30.))*d350i 
                       ! ! ---- third derivatives ----
                       ! !   1/4, 7/120, 41/3024
                       ! !    uxxxx = D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ... 
                       ! uxxx = d300i - (dx(0)**2/4.)*d500i + (dx(0)**4 *7./120.)*d700i
                       ! uyyy = d030i - (dx(1)**2/4.)*d050i + (dx(1)**4 *7./120.)*d070i
                       ! ! uxxy = D+xD-x*( 1 - dx^2/12 D+xD-x + dx^4/90*(D+xD-x)^2 - (1/560)*dx^6*(D+xD-x)^3 ) X
                       ! !           D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                       ! uxxy = d210i - (dx(0)**2/12.)*d410i - (dx(1)**2/6.)*d230i + (dx(0)**4/90.)*d610i + (dx(1)**4/30.)*d250i + (dx(0)**2 * dx(1)**2/72. )*d430i
                       ! ! uxyy =    D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                       ! !        D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                       ! uxyy = d120i - (dx(1)**2/12.)*d140i - (dx(0)**2/6.)*d320i + (dx(1)**4/90.)*d160i + (dx(0)**4/30.)*d520i + (dx(0)**2 * dx(1)**2/72. )*d340i      
                       ! ----- fourth derivs ------
                       !  1/6, 7/240, 41/7560
                       !   uxxxx = (D+xD-x)^2*( 1 - (1/6)*D+xD-x + (7/240)*(D+xD-x)^2 + ... )
                                              uxxxx = d400(i1,i2,i3,0) - (dx(0)**2/6.)*d600(i1,i2,i3,0) + (dx(0)**4 * (7./240))*d800i
                                              uyyyy = d040(i1,i2,i3,0) - (dx(1)**2/6.)*d060(i1,i2,i3,0) + (dx(1)**4 * (7./240))*d080i
                       ! ! uxyyy = D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                       ! !         D0y D+yD-y*( 1 - (1/4)*dy^2 * D+yD-y + (7/120)*dy^4 (D+yD-y)^2 + ... 
                       ! uxyyy = d130i - (dx(0)**2/6.)*d330i - (dx(1)**2/4.)*d150i + (dx(0)**4/30.)*d530i + (dx(1)**4*7./120.)*d170i + (dx(0)**2 * dx(1)**2/24. )*d350i 
                       ! ! uxxxy
                       ! uxxxy = d310i - (dx(1)**2/6.)*d330i - (dx(0)**2/4.)*d510i + (dx(1)**4/30.)*d350i + (dx(0)**4*7./120.)*d710i + (dx(0)**2 * dx(1)**2/24. )*d530i       
                       ! uxxyy = D+xD-x*( 1 - dx^2/12 D+xD-x + dx^4/90*(D+xD-x)^2 - (1/560)*dx^6*(D+xD-x)^3 ) X 
                       !         D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                              uxxyy = d220(i1,i2,i3,0) - (dx(0)**2/12.)*d420(i1,i2,i3,0) - (dx(1)**2/12.)*d240(i1,i2,i3,0) + (dx(0)**4 * (1./90))*d620i + (dx(1)**4 * (1./90))*d260i + (dx(0)**2 * dx(1)**2 * 1./144.)*d440i 
                       ! ! ---- fixth derivatives ----
                       ! !   1/3, 13/144, 139/6048
                       ! !  uxxxxx = D0x(D+xD-x)^2 *( 1 - (1/3)*dx^2 D+xD-x + ... )
                       ! uxxxxx = d500i - (dx(0)**2/3.)*d700i
                       ! uyyyyy = d050i - (dx(1)**2/3.)*d070i 
                       ! ! uxxxxy = (D+xD-x)^2 *( 1 - dx^2/6 D+xD-x + ...) X
                       ! !                  D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                       ! uxxxxy = d410i - (dx(0)**2/6.)*d610i - (dx(1)**2/6.)*d430i
                       ! uxyyyy = d140i - (dx(1)**2/6.)*d160i - (dx(0)**2/6.)*d340i    
                       ! ! uxxxyy =  D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ...
                       ! !               D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                       ! uxxxyy = d320i - (dx(0)**2/4.)*d520i - (dx(1)**2/12.)*d340i
                       ! uxxyyy = d230i - (dx(1)**2/4.)*d250i - (dx(0)**2/12.)*d430i
                       ! ---- sixth derivs----
                       !  1/4 , 13/240 
                                              uxxxxxx = d600i - (dx(0)**2/4.)*d800i
                                              uyyyyyy = d060i - (dx(1)**2/4.)*d080i  
                       ! ! uxxxxxy = D0x(D+xD-x)^2 *( 1 - (1/3)*dx^2 D+xD-x + ... ) X 
                       ! !                      D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                       ! uxxxxxy = d510i - (dx(0)**2/3.)*d710i - (dx(1)**2/6.)*d530i
                       ! uxyyyyy = d150i - (dx(1)**2/3.)*d170i - (dx(0)**2/6.)*d350i
                       ! uxxxxyy = (D+xD-x)^2*( 1 - (1/6)*D+xD-x + (7/240)*(D+xD-x)^2 + ... )
                       !               D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                              uxxxxyy = d420i - (dx(0)**2/6.)*d620i - (dx(1)**2/12.)*d440i
                                              uxxyyyy = d240i - (dx(1)**2/6.)*d260i - (dx(0)**2/12.)*d440i
                       ! ! uxxxyyy = D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ... ) X 
                       ! !           D0y D+yD-y*( 1 - (1/4)*dy^2 * D+yD-y + (7/120)*dy^4 (D+yD-y)^2 + ...  ) u 
                       ! uxxxyyy = d330i - (dx(0)**2/4.)*d530i - (dx(1)**2/4.)*d350i
                       ! ! ---- seventh derivatives ----
                       ! uxxxxxxx = d700i
                       ! uxxxxxxy = d610i
                       ! uxxxxxyy = d520i
                       ! uxxxxyyy = d430i
                       ! uxxxyyyy = d340i
                       ! uxxyyyyy = d250i
                       ! uxyyyyyy = d160i
                       ! uyyyyyyy = d070i 
                       ! ---- eighth derivatives ----
                                              uxxxxxxxx = d800i
                       ! uxxxxxxxy = d710i
                                              uxxxxxxyy = d620i
                       ! uxxxxxyyy = d530i
                                              uxxxxyyyy = d440i
                       ! uxxxyyyyy = d350i
                                              uxxyyyyyy = d260i
                       ! uxyyyyyyy = d170i
                                              uyyyyyyyy = d080i                  
                                                  if( forcingOption.eq.twilightZoneForcing )then
                                                              call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                                              call ogDeriv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evtt(m) )
                                                              call ogDeriv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxx(m) )
                                                              call ogDeriv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evyy(m) )
                                                          fv(m) = evtt(m) - csq*( evxx(m) + evyy(m) )
                              ! Correct forcing for fourth-order ME in2D
                                                                call ogDeriv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evtttt(m) )
                                                                call ogDeriv(ep, 0,4,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxxxx(m) )
                                                                call ogDeriv(ep, 0,2,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxxyy(m) )
                                                                call ogDeriv(ep, 0,0,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evyyyy(m) )
                                                            fv(m) = fv(m) + (dtSq/12.)*evtttt(m) - (cdtsq12/dtSq)*( evxxxx(m) + 2.*evxxyy(m) + evyyyy(m) )
                              ! Correct forcing for sixth-order ME in 2D
                                                                call ogDeriv(ep, 6,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evtttttt(m) )
                                                                call ogDeriv(ep, 0,6,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxxxxxx(m) )
                                                                call ogDeriv(ep, 0,4,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxxxxyy(m) )
                                                                call ogDeriv(ep, 0,2,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxxyyyy(m) )
                                                                call ogDeriv(ep, 0,0,6,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evyyyyyy(m) )
                                                            fv(m) = fv(m) + (dtSq**2/360.)*evtttttt(m) - (cdtPow6By360/dtSq)*( evxxxxxx(m) + evyyyyyy(m) + 3.*(evxxxxyy(m) + evxxyyyy(m) )  )
                              ! Correct forcing for eighth-order ME in 2D
                                                                call ogDeriv(ep, 8,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,uet8 )
                                                                call ogDeriv(ep, 0,8,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,uex8 )
                                                                call ogDeriv(ep, 0,6,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,uex6y2 )
                                                                call ogDeriv(ep, 0,4,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,uex4y4 )
                                                                call ogDeriv(ep, 0,2,6,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,uex2y6 )
                                                                call ogDeriv(ep, 0,0,8,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,uey8 )
                              ! (x*x + y*y)^4 = x^8 + 4*x^6*y^2 + 6*x^4*y^4 + 4*x^2*y^6 + y^8
                                                            fv(m) = fv(m) + (dtSq**3/20160.)*uet8 - (cdtPow8By20160/dtSq)*( uex8 + uey8   + 4.*(uex6y2 + uex2y6 ) +6.*uex4y4 )
                                                else if( forcingOption.eq.helmholtzForcing )then
                           ! forcing for solving the Helmholtz equation   
                           ! NOTE: change sign of forcing since for Helholtz we want to solve
                           !      ( omega^2 I + c^2 Delta) w = f 
                           ! fv(m) = -f(i1,i2,i3,0)*coswt  
                                                      fv(m)=0.
                                                      do freq=0,numberOfFrequencies-1 
                                                          omega = frequencyArray(freq)
                                                          coswt = cosFreqt(freq)    
                               ! Add corrections for 4th order modified equation 
                               !  fv = f + (dt^2/12)*( c^2 Delta(u) + ftt )
                                                              write(*,*) 'fix me'
                                                              stop 4444
                                   !fv(m) = fv(m) -( f(i1,i2,i3,freq) + cdtSqBy12*( cSq*(fxx22r(i1,i2,i3,freq) + fyy22r(i1,i2,i3,freq)) - omega*omega*f(i1,i2,i3,freq)) )*coswt 
                                                      end do ! do freq  
                                                else if( addForcing.ne.0 )then  
                                                      fv(m) = f(i1,i2,i3,0)
                                                end if
                       ! Here is the ME update 
                                              un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx + uyy ) + cdtPow4By12*( uxxxx + uyyyy + 2.*uxxyy )  + cdtPow6By360*( uxxxxxx + uyyyyyy +3.*( uxxxxyy + uxxyyyy) ) + cdtPow8By20160*( uxxxxxxxx + uyyyyyyyy + 4.*( uxxxxxxyy +uxxyyyyyy ) +6.*uxxxxyyyy ) + dtSq*fv(m)      
                     ! if( i1.eq.5 .and. i2.eq.6 )then
                     !   call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t+dt,uc,ue )
                     !   write(*,'("ME8: (i1,i2)=(",2i3,") u=",1pe12.4," ue=",1pe12.4," err=",1pe8.2)') i1,i2,un(i1,i2,i3,0),ue,abs(un(i1,i2,i3,0)-ue)
                     ! end if
                                          end if ! mask .ne. 0
                                        end do
                                        end do
                                        end do
             ! #If 8 == 6 || 8 == 8 
             !   evalDerivativesRectangular()
             !   write(*,*) ' Stop here for now'
             !   stop 666
             ! #End
           !   ! --- TAYLOR TIME-STEPPING --- 
           !   m=0 ! component number 
           !   ec = 0 ! component number
           !   ! #If "rectangular" eq "curvilinear"
           !   !   #If "8" eq "4" && "8" eq "4"
           !   !     computeLaplacianOrder2(2)
           !   !   #End
           !   ! #End
           !   if( forcingOption.eq.helmholtzForcing )then
           !     coswt = cos(omega*t)
           !   end if 
           !   fv(m)=0.
           !   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
           !     getForcing(2,8,8,rectangular) 
           !    #If "8" eq "2"
           !      ! --- SECOND 8 ---
           !      #If "2" eq "2"
           !        ! --- TWO DIMENSIONS ---
           !        #If "rectangular" eq "rectangular"
           !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx22r(i1,i2,i3,0) + uyy22r(i1,i2,i3,0) ) + dtSq*fv(m)
           !         ! write(*,'(" adv: i1,i2=",2i4," un,u,um=",3e12.2," cdtSq,fv=",2e12.2)') i1,i2,un(i1,i2,i3,m),u(i1,i2,i3,m),um(i1,i2,i3,m),cdtSq,fv(m)
           !        #Else
           !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx22(i1,i2,i3,0)  + uyy22(i1,i2,i3,0) ) + dtSq*fv(m)
           !        #End
           !      #Else
           !        ! --- THREE DIMENSIONS ---
           !        #If "rectangular" eq "rectangular"
           !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx23r(i1,i2,i3,0) + uyy23r(i1,i2,i3,0) + uzz23r(i1,i2,i3,0) ) + dtSq*fv(m)
           !        #Else
           !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx23(i1,i2,i3,0)  + uyy23(i1,i2,i3,0)  + uzz23(i1,i2,i3,0)  ) + dtSq*fv(m)
           !        #End
           !      #End
           !    #Elif "8" eq "4"
           !      ! --- -FOURTH 8 ---
           !      #If "2" eq "2"
           !        ! --- FOUTH-8 TWO DIMENSIONS ---
           !        #If "8" eq "4"
           !          ! orderInSpace=4 and orderInTime=4 
           !          #If "rectangular" eq "rectangular"
           !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + cdtsq12*lap2d2pow2(i1,i2,i3,m) + dtSq*fv(m)
           !          #Else
           !            ! v is assumed to hold Lap(u) to 2nd-order
           !            ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
           !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !                                                            + cdtsq12*( vxx22(i1,i2,i3,m) + vyy22(i1,i2,i3,m) ) + dtSq*fv(m)
           !          #End
           !        #Else
           !          ! orderInSpace==4 and orderInTime==2                                                   
           !          #If "rectangular" eq "rectangular"
           !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + dtSq*fv(m)
           !          #Else
           !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End
           !        #End                                                          
           !      #Else
           !        ! --- FOURTH-8 THREE DIMENSIONS ---
           !        #If "8" eq "4"
           !          ! orderInSpace=4 and orderInTime=4 
           !          #If "rectangular" eq "rectangular"
           !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + cdtsq12*lap3d2pow2(i1,i2,i3,m) + dtSq*fv(m)
           !          #Else
           !            ! v is assumed to hold Lap(u) to 2nd-order
           !            ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
           !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !                                                            + cdtsq12*( vxx23(i1,i2,i3,m) + vyy23(i1,i2,i3,m) + vzz23(i1,i2,i3,m) ) + dtSq*fv(m)
           !          #End
           !        #Else
           !          ! orderInSpace==4 and orderInTime==2                                                   
           !          #If "rectangular" eq "rectangular"
           !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + dtSq*fv(m)
           !          #Else
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End
           !        #End     
           !      #End
           !    #Elif "8" eq "6"
           !      ! ---- SIXTH 8 ---
           !      #If "2" eq "2"
           !        ! --- SIXTH-8 TWO DIMENSIONS ---
           !        #If "8" eq "2"
           !          ! orderInSpace==6 and orderInTime==2                                                   
           !          #If "rectangular" eq "rectangular"
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx62r(i1,i2,i3,m) + uyy62r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #Else
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx62(i1,i2,i3,m)  +  uyy62(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End       
           !        #Else
           !          stop 6666                                                           
           ! !          ! ---- MODIFIED EQUATION 8=6 2D -----
           ! !          getSixthDerivatives2d(8,rectangular,evalMetrics,i1,i2,i3)
           ! !          if( i1.eq.21 .and. i2.eq.11 )then
           ! !            OGDERIV2D( 0,2,0,0,i1,i2,i3,t, ec, evxx(0) )
           ! !            OGDERIV2D( 0,0,2,0,i1,i2,i3,t, ec, evyy(0) )
           ! !            OGDERIV2D( 0,4,0,0,i1,i2,i3,t, ec, evxxxx(0) )
           ! !            OGDERIV2D( 0,0,4,0,i1,i2,i3,t, ec, evyyyy(0) )
           ! !            OGDERIV2D( 0,2,2,0,i1,i2,i3,t, ec, evxxyy(0) )
           ! !            write(*,'("i1,i2=",2i4)') i1,i2
           ! !            write(*,'(" uxx  =",1pe12.4," true=",1pe12.4," err=",1pe8.2)') uxx  ,evxx(0)  ,abs(uxx  -evxx(0)  )
           ! !            write(*,'(" uyy  =",1pe12.4," true=",1pe12.4," err=",1pe8.2)') uyy  ,evyy(0)  ,abs(uyy  -evyy(0)  )
           ! !            write(*,'(" uxxxx=",1pe12.4," true=",1pe12.4," err=",1pe8.2)') uxxxx,evxxxx(0),abs(uxxxx-evxxxx(0))
           ! !            write(*,'(" uxxyy=",1pe12.4," true=",1pe12.4," err=",1pe8.2)') uxxyy,evxxyy(0),abs(uxxyy-evxxyy(0))
           ! !            write(*,'(" uyyyy=",1pe12.4," true=",1pe12.4," err=",1pe8.2)') uyyyy,evyyyy(0),abs(uyyyy-evyyyy(0))
           ! ! ! uxxxx    = rxi**4*urrrr+4.*rxi**3*sxi*urrrs+6.*rxi**2*sxi**2*urrss+4.*rxi*sxi**3*ursss+sxi**4*ussss+6.*rxi**2*rxx*urrr+(7.*sxi*rxi*rxx+sxx*rxi**2+rxi*(3.*rxi*sxx+3.*rxx*sxi)+rxi*(2.*rxi*sxx+2.*rxx*sxi))*urrs+(sxi*(3.*rxi*sxx+3.*rxx*sxi)+7.*rxi*sxx*sxi+rxx*sxi**2+sxi*(2.*rxi*sxx+2.*rxx*sxi))*urss+6.*sxi**2*sxx*usss+(4.*rxi*rxxx+3.*rxx**2)*urr+(4.*rxi*sxxx+6.*rxx*sxx+4.*rxxx*sxi)*urs+(4.*sxi*sxxx+3.*sxx**2)*uss+rxxxx*ur+sxxxx*us
           ! !            stop 1111
           ! !          end if
           ! !          ! cdtPow4By12  = (c*dt)^4/12 
           ! !          ! cdtPow6By360 = (c*dt)^6/360
           ! !          un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) ! !                          + cdtsq*( uxx + uyy ) ! !                          + cdtPow4By12*( uxxxx + uyyyy + 2.*uxxyy )  ! !                          + cdtPow6By360*( uxxxxxx + uyyyyyy + 3.*(uxxxxyy + uxxyyyy) ) ! !                          + dtSq*fv(m)
           !       #End                                                          
           !     #Else
           !        ! --- SIXTH-8 THREE DIMENSIONS ---
           !        #If "8" eq "2
           !          ! orderInSpace==6 and orderInTime==2                                                   
           !          #If "rectangular" eq "rectangular"
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx63r(i1,i2,i3,m) + uyy63r(i1,i2,i3,m) + uzz63r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #Else
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx63(i1,i2,i3,m)  + uyy63(i1,i2,i3,m)  + uzz63(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End       
           !        #Else
           !          ! MODIFIED EQUATION 8=6 3D
           !          ! Turn off for now: 
           !          ! getSixthDerivatives3d(8,rectangular,evalMetrics,i1,i2,i3)
           !          write(*,'("advWave: order=8, orderInTime=8 FINISH ME")')
           !          stop 7777
           !          ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) !          !                 + cdtsq*( uxx + uyy +uzz ) !          !                 + cdtPow4By12*( uxxxx + uyyyy + uzzzz + 2.*( uxxyy +uxxzz + uyyzz ) )  !          !                 + cdtPow6By360*( uxxxxxx +  uyyyyyy + uzzzzzz + 3.*(uxxxxyy + uxxyyyy + uxxxxzz + uyyyyzz + uxxzzzz + uyyzzzz ) + 6.*uxxyyzz ) !          !                 + dtSq*fv(m)
           !        #End     
           !      #End
           !    #Elif "8" eq "8"
           !      ! ---- EIGTH 8 ---
           !      #If "2" eq "2"
           !        ! --- EIGTH-8 TWO DIMENSIONS ---
           !        #If "8" eq "2"
           !          ! orderInSpace==8 and orderInTime==2                                                   
           !          #If "rectangular" eq "rectangular"
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx82r(i1,i2,i3,m) + uyy82r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #Else
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx82(i1,i2,i3,m)  +  uyy82(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End       
           !        #Else
           !          write(*,'("advWave: order=8, orderInTime=8 FINISH ME")')
           !          stop 7777
           !          ! ! orderInSpace=4 and orderInTime=4 
           !          ! #If "rectangular" eq "rectangular"
           !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + cdtsq12*lap2d2pow2(i1,i2,i3,m) + dtSq*fv(m)
           !          ! #Else
           !          !   ! v is assumed to hold Lap(u) to 2nd-order
           !          !   ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
           !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !          !                                                   + cdtsq12*( vxx22(i1,i2,i3,m) + vyy22(i1,i2,i3,m) ) + dtSq*fv(m)
           !          ! #End
           !        #End                                                          
           !     #Else
           !        ! --- EIGTH-8 THREE DIMENSIONS ---
           !        #If "8" eq "2
           !          ! orderInSpace==8 and orderInTime==2                                                   
           !          #If "rectangular" eq "rectangular"
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx83r(i1,i2,i3,m) + uyy83r(i1,i2,i3,m) + uzz83r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #Else
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx83(i1,i2,i3,m)  + uyy83(i1,i2,i3,m)  + uzz83(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End       
           !        #Else
           !          write(*,'("advWave: order=8, orderInTime=8 FINISH ME")')
           !          stop 7777
           !          ! #If "rectangular" eq "rectangular"
           !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + cdtsq12*lap3d2pow2(i1,i2,i3,m) + dtSq*fv(m)
           !          ! #Else
           !          !   ! v is assumed to hold Lap(u) to 2nd-order
           !          !   ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
           !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !          !                                                   + cdtsq12*( vxx23(i1,i2,i3,m) + vyy23(i1,i2,i3,m) + vzz23(i1,i2,i3,m) ) + dtSq*fv(m)
           !          ! #End
           !        #End     
           !      #End
           !    #Else
           !      write(*,'("advWave: UNKNOWN order=8")')
           !      stop 7777
           !    #End         
           !    ! write(*,'("i1,i2=",2i3," u-ue=",e10.2)') i1,i2,u(i1,i2,i3,m)-ev(m)
           !    ! write(*,'(" uxx-uxxe =",e10.2)') uxx22r(i1,i2,i3,0)-evxx(m)
           !    ! OGDERIV2D( 0,0,0,0,i1,i2,i3,t+dt, ec, ev(m)  )
           !    ! write(*,'(" un-ue=",e10.2)') un(i1,i2,i3,m)-ev(m)
           !   endLoopsMask()
                  end if 
          else
       ! --- IMPLICIT: Fill in RHS to implicit time-stepping -----
              stop 1111
          end if
          return
          end
