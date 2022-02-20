! This file automatically generated from advWaveME.bf90 with bpp.
        subroutine advWaveME2dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
    ! subroutine advWaveME2dOrder6c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,!                 mask,xy,rsxy,  um,u,un, f,fa, v, vh,  bc, frequencyArray, ipar, rpar, ierr )
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
          real d400i,d040i
          real d600i,d060i, d240i, d420i , d510i,d150i,d330i
          real d800i,d080i, d620i, d260i, d440i
          real d100i, d010i, d300i,d030i, d500i, d050i, d700i, d070i
          real d110i, d310i,d130i, d710i,d170i,d350i,d530i
          real d120i,d210i, d140i,d410i, d320i,d230i, d160i, d610i, d520i, d250i, d340i, d430i
          real rx200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
          real rx020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
          real rx400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
          real rx220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
          real rx040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
          real rxxa(0:2,0:2,0:2)
          real rx000i(0:2,0:2)
          real rx100i(0:2,0:2)
          real rx010i(0:2,0:2)
          real rx110i(0:2,0:2)
          real rx300i(0:2,0:2)
          real rx030i(0:2,0:2)
          real rx210i(0:2,0:2)
          real rx120i(0:2,0:2)  
          real rx310i(0:2,0:2)
          real rx130i(0:2,0:2) 
          real rx500i(0:2,0:2)
          real rx050i(0:2,0:2)
          real rx410i(0:2,0:2)
          real rx140i(0:2,0:2)  
          real rx320i(0:2,0:2)
          real rx230i(0:2,0:2)  
          real rx510i(0:2,0:2)
          real rx150i(0:2,0:2)
          real rx330i(0:2,0:2)
          real rx600i(0:2,0:2)
          real rx060i(0:2,0:2)
          real rx420i(0:2,0:2)
          real rx240i(0:2,0:2)  
          real rxr(0:2,0:2,0:2)
          real rxrra(0:2,0:2),rxssa(0:2,0:2),rxrsa(0:2,0:2),rxxxa(0:2,0:2),rxxya(0:2,0:2),rxyya(0:2,0:2)
          real rxrrra(0:2,0:2),rxrrsa(0:2,0:2),rxrssa(0:2,0:2),rxsssa(0:2,0:2)
          real rxxxxa(0:2,0:2),rxxxya(0:2,0:2),rxxyya(0:2,0:2),rxyyya(0:2,0:2)
          real rxrrrra(0:2,0:2),rxssssa(0:2,0:2),rxrsssa(0:2,0:2),rxrrrsa(0:2,0:2),rxrrssa(0:2,0:2)
          real rxxxxxa(0:2,0:2),rxxxxya(0:2,0:2),rxxxyya(0:2,0:2),rxxyyya(0:2,0:2),rxyyyya(0:2,0:2)
          real rxrrrrra(0:2,0:2),rxsssssa(0:2,0:2),rxrrrrsa(0:2,0:2),rxrssssa(0:2,0:2),rxrrrssa(0:2,0:2),rxrrsssa(0:2,0:2)
          real rxxxxxxa(0:2,0:2),rxxxxxya(0:2,0:2),rxxxxyya(0:2,0:2),rxxxyyya(0:2,0:2),rxxyyyya(0:2,0:2),rxyyyyya(0:2,0:2)
          real urr,urs,uss,ur,us
          real urrr,urss,urrs,usss
          real urrrr,urrrs,urrss,ursss,ussss
          real urrrrr,urrrrs,urrrss,urrsss,urssss,usssss
          real urrrrrr,urrrrrs,urrrrss,urrrsss,urrssss,ursssss,ussssss
          real urrrrrrr,urrrrrrs,urrrrrss,urrrrsss,urrrssss,urrsssss,urssssss,usssssss
          real urrrrrrrr,urrrrrrrs,urrrrrrss,urrrrrsss,urrrrssss,urrrsssss,urrssssss,ursssssss,ussssssss
          real rx,ry,sx,sy
          real rxx,ryy,sxx,syy,rxy,sxy
          real rxxx,ryxx,sxxx,syxx
          real rxxy,ryxy,sxxy,syxy
          real rxyy,ryyy,sxyy,syyy
          real rxxxx, ryxxx, sxxxx, syxxx
          real rxxxy, ryxxy, sxxxy, syxxy
          real rxxyy, ryxyy, sxxyy, syxyy
          real rxyyy, ryyyy, sxyyy, syyyy
          real rxxxxx, ryxxxx, sxxxxx, syxxxx
          real rxxxxy, ryxxxy, sxxxxy, syxxxy
          real rxxxyy, ryxxyy, sxxxyy, syxxyy
          real rxxyyy, ryxyyy, sxxyyy, syxyyy
          real rxyyyy , ryyyyy , sxyyyy , syyyyy 
          real rxxxxxx, ryxxxxx, sxxxxxx, syxxxxx
          real rxxxxxy, ryxxxxy, sxxxxxy, syxxxxy
          real rxxxxyy, ryxxxyy, sxxxxyy, syxxxyy
          real rxxxyyy, ryxxyyy, sxxxyyy, syxxyyy
          real rxxyyyy, ryxyyyy, sxxyyyy, syxyyyy
          real rxyyyyy, ryyyyyy, sxyyyyy, syyyyyy
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
                            write(*,'("advWaveME: ADVANCE dim=2 order=6 orderInTime=2, grid=curvilinear... t=",e10.2)') t
                        end if
                        m=0 ! component number 
                        ec = 0 ! component number  
                                    numGhost1=2; ! should depend on the orderOfAccuracy
                                    n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                    n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                    n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                      do i3=n3a,n3b
                                      do i2=n2a,n2b
                                      do i1=n1a,n1b
                                        if( mask(i1,i2,i3).ne.0 )then
                                            d200(i1,i2,i3,0) = (u(i1+1,i2,i3,0) - 2.*u(i1,i2,i3,0) + u(i1-1,i2,i3,0))/(dx(0)**2)
                                            d020(i1,i2,i3,0) = (u(i1,i2+1,i3,0) - 2.*u(i1,i2,i3,0) + u(i1,i2-1,i3,0))/(dx(1)**2)
                      ! We really only need this for m2>=m1  ** FIX ME **
                                            do m1=0,nd-1
                                                do m2=0,nd-1
                                                    rx200(i1,i2,i3,m1,m2) = (rsxy(i1+1,i2,i3,m1,m2)- 2.*rsxy(i1,i2,i3,m1,m2) + rsxy(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                                                    rx020(i1,i2,i3,m1,m2) = (rsxy(i1,i2+1,i3,m1,m2)- 2.*rsxy(i1,i2,i3,m1,m2) + rsxy(i1,i2-1,i3,m1,m2))/(dx(1)**2)          
                                                end do
                                            end do
                                        end if ! mask .ne. 0
                                      end do
                                      end do
                                      end do
                                    numGhost1=1; ! should depend on the orderOfAccuracy
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
                      ! We really only need this for m2>=m1  ** FIX ME **
                                            do m1=0,nd-1
                                                do m2=0,nd-1
                                                    rx400(i1,i2,i3,m1,m2) = (rx200(i1+1,i2,i3,m1,m2) - 2.*rx200(i1,i2,i3,m1,m2) + rx200(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                                                    rx220(i1,i2,i3,m1,m2) = (rx200(i1,i2+1,i3,m1,m2) - 2.*rx200(i1,i2,i3,m1,m2) + rx200(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                                                    rx040(i1,i2,i3,m1,m2) = (rx020(i1,i2+1,i3,m1,m2) - 2.*rx020(i1,i2,i3,m1,m2) + rx020(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                                                end do
                                            end do
                                        end if ! mask .ne. 0
                                      end do
                                      end do
                                      end do
                  ! numGhost1=1; ! should depend on the orderOfAccuracy
                  ! n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                  ! n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                  ! n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                  ! beginLoops3d()
                  !   if( mask(i1,i2,i3).ne.0 )then
                  !     d600(i1,i2,i3,0) = (d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0))/(dx(0)**2)
                  !     d060(i1,i2,i3,0) = (d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0))/(dx(1)**2)
                  !     d420(i1,i2,i3,0) = (d220(i1+1,i2,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0))/(dx(0)**2)
                  !     d240(i1,i2,i3,0) = (d220(i1,i2+1,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1,i2-1,i3,0))/(dx(1)**2)
                  !     ! We really only need this for m2>=m1  ** FIX ME **
                  !     do m1=0,nd-1
                  !       do m2=0,nd-1
                  !         rx600(i1,i2,i3,m1,m2) = (rx400(i1+1,i2,i3,m1,m2) - 2.*rx400(i1,i2,i3,m1,m2) + rx400(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                  !         rx060(i1,i2,i3,m1,m2) = (rx040(i1,i2+1,i3,m1,m2) - 2.*rx040(i1,i2,i3,m1,m2) + rx040(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                  !         rx420(i1,i2,i3,m1,m2) = (rx220(i1+1,i2,i3,m1,m2) - 2.*rx220(i1,i2,i3,m1,m2) + rx220(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                  !         rx240(i1,i2,i3,m1,m2) = (rx220(i1,i2+1,i3,m1,m2) - 2.*rx220(i1,i2,i3,m1,m2) + rx220(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                  !       end do
                  !     end do      
                  !   end if ! mask .ne. 0
                  ! endLoops3d()     
                                    numGhost1=0; 
                                    n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                    n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                    if( nd.eq.2 )then
                                        n3a=gridIndexRange(0,2); n3b=gridIndexRange(1,2);
                                    else
                                        n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                    end if
                  ! See highOrderDiff.maple
                  ! u.xx = D+D-( I - b1*dx^2*D+D- + b2*(dx^2*D+D-)^2 -b3*(dx^2*D+D-)^2 + ...)
                  ! b1 = 1/12
                  ! b2 = 1/90
                  ! b3= 1/560
                  ! b4= 1/3150
                  ! ----------- MAIN LOOP : 6 6 CURVILINEAR --------
                                    fv(m)=0.
                                      do i3=n3a,n3b
                                      do i2=n2a,n2b
                                      do i1=n1a,n1b
                                        if( mask(i1,i2,i3).ne.0 )then
                                            d600i = (d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0))/(dx(0)**2)
                                            d060i = (d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0))/(dx(1)**2)
                                            d420i = (d220(i1+1,i2,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0))/(dx(0)**2)
                                            d240i = (d220(i1,i2+1,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1,i2-1,i3,0))/(dx(1)**2)
                      ! We really only need this for m2>=m1  ** FIX ME **
                                            do m1=0,nd-1
                                                do m2=0,nd-1
                                                    rx600i(m1,m2) = (rx400(i1+1,i2,i3,m1,m2) - 2.*rx400(i1,i2,i3,m1,m2) + rx400(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                                                    rx060i(m1,m2) = (rx040(i1,i2+1,i3,m1,m2) - 2.*rx040(i1,i2,i3,m1,m2) + rx040(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                                                    rx420i(m1,m2) = (rx220(i1+1,i2,i3,m1,m2) - 2.*rx220(i1,i2,i3,m1,m2) + rx220(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                                                    rx240i(m1,m2) = (rx220(i1,i2+1,i3,m1,m2) - 2.*rx220(i1,i2,i3,m1,m2) + rx220(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                                                end do
                                            end do       
                      ! d600i = (d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0))/(dx(0)**2)  
                      ! d060i = (d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0))/(dx(1)**2)
                      ! d800i = (d600(i1+1,i2,i3,0) - 2.*d600(i1,i2,i3,0) + d600(i1-1,i2,i3,0))/(dx(0)**2)
                      ! d080i = (d060(i1,i2+1,i3,0) - 2.*d060(i1,i2,i3,0) + d060(i1,i2-1,i3,0))/(dx(1)**2)
                      ! d620i = (d420(i1+1,i2,i3,0) - 2.*d420(i1,i2,i3,0) + d420(i1-1,i2,i3,0))/(dx(0)**2)  
                      ! d260i = (d240(i1,i2+1,i3,0) - 2.*d240(i1,i2,i3,0) + d240(i1,i2-1,i3,0))/(dx(1)**2)
                      ! d440i = (d240(i1+1,i2,i3,0) - 2.*d240(i1,i2,i3,0) + d240(i1-1,i2,i3,0))/(dx(0)**2)  
                      ! d240i = d240(i1,i2,i3,0)    
                      ! d420i = d420(i1,i2,i3,0)    
                      ! ----- for odd derivatives ----
                                            d100i = (u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0))/(2.*dx(0))   ! D0x 
                                            d010i = (u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0))/(2.*dx(1))   ! D0y 
                                            d300i = (d200(i1+1,i2,i3,0) - d200(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x D_x D-x 
                                            d030i = (d020(i1,i2+1,i3,0) - d020(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y D+y D-y 
                                            d120i = (d020(i1+1,i2,i3,0) - d020(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x D_y D-y 
                                            d210i = (d200(i1,i2+1,i3,0) - d200(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y D+x D-x       
                                            d500i = (d400(i1+1,i2,i3,0) - d400(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_x D-x)^2 
                                            d050i = (d040(i1,i2+1,i3,0) - d040(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+y D-y)^2 
                                            d140i = (d040(i1+1,i2,i3,0) - d040(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_y D-y)^2 
                                            d410i = (d400(i1,i2+1,i3,0) - d400(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)^2
                                            d320i = (d220(i1+1,i2,i3,0) - d220(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D+x D-x)(D_y D-y)
                                            d230i = (d220(i1,i2+1,i3,0) - d220(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)
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
                      ! --- for mixed derivatives
                                            d110i = (u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))                ! D0x D0y
                                            d310i = (d200(i1+1,i2+1,i3,0) - d200(i1-1,i2+1,i3,0) - d200(i1+1,i2-1,i3,0) + d200(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) 
                                            d130i = (d020(i1+1,i2+1,i3,0) - d020(i1-1,i2+1,i3,0) - d020(i1+1,i2-1,i3,0) + d020(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)          
                                            d510i = (d400(i1+1,i2+1,i3,0) - d400(i1-1,i2+1,i3,0) - d400(i1+1,i2-1,i3,0) + d400(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2
                                            d150i = (d040(i1+1,i2+1,i3,0) - d040(i1-1,i2+1,i3,0) - d040(i1+1,i2-1,i3,0) + d040(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^2     
                                            d330i = (d220(i1+1,i2+1,i3,0) - d220(i1-1,i2+1,i3,0) - d220(i1+1,i2-1,i3,0) + d220(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)       
                      ! d710i = (d600(i1+1,i2+1,i3,0) - d600(i1-1,i2+1,i3,0) - d600(i1+1,i2-1,i3,0) + d600(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^3
                      ! d170i = (d060(i1+1,i2+1,i3,0) - d060(i1-1,i2+1,i3,0) - d060(i1+1,i2-1,i3,0) + d060(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^3  
                      ! d530i = (d420(i1+1,i2+1,i3,0) - d420(i1-1,i2+1,i3,0) - d420(i1+1,i2-1,i3,0) + d420(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2 (D+yD-y)         
                      ! d350i = (d240(i1+1,i2+1,i3,0) - d240(i1-1,i2+1,i3,0) - d240(i1+1,i2-1,i3,0) + d240(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)^2         
                                            do m1=0,nd-1
                                                do m2=0,nd-1
                                                    rx000i(m1,m2) = rsxy(i1,i2,i3,m1,m2)
                                                    rx100i(m1,m2) = (rsxy(i1+1,i2,i3,m1,m2)   - rsxy(i1-1,i2,i3,m1,m2))/(2.*dx(0))   ! D0x 
                                                    rx010i(m1,m2) = (rsxy(i1,i2+1,i3,m1,m2)   - rsxy(i1,i2-1,i3,m1,m2))/(2.*dx(1))   ! D0y 
                                                    rx110i(m1,m2) = (rsxy(i1+1,i2+1,i3,m1,m2) - rsxy(i1-1,i2+1,i3,m1,m2) - rsxy(i1+1,i2-1,i3,m1,m2) + rsxy(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1)) 
                                                    rx300i(m1,m2) = (rx200(i1+1,i2,i3,m1,m2) - rx200(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x D_x D-x 
                                                    rx030i(m1,m2) = (rx020(i1,i2+1,i3,m1,m2) - rx020(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y D+y D-y  
                                                    rx120i(m1,m2) = (rx020(i1+1,i2,i3,m1,m2) - rx020(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x D_y D-y 
                                                    rx210i(m1,m2) = (rx200(i1,i2+1,i3,m1,m2) - rx200(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y D+x D-x             
                                                    rx310i(m1,m2) = (rx200(i1+1,i2+1,i3,m1,m2) - rx200(i1-1,i2+1,i3,m1,m2) - rx200(i1+1,i2-1,i3,m1,m2) + rx200(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) 
                                                    rx130i(m1,m2) = (rx020(i1+1,i2+1,i3,m1,m2) - rx020(i1-1,i2+1,i3,m1,m2) - rx020(i1+1,i2-1,i3,m1,m2) + rx020(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)       
                                                    rx500i(m1,m2) = (rx400(i1+1,i2,i3,m1,m2) - rx400(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D_x D-x)^2 
                                                    rx050i(m1,m2) = (rx040(i1,i2+1,i3,m1,m2) - rx040(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+y D-y)^2 
                                                    rx140i(m1,m2) = (rx040(i1+1,i2,i3,m1,m2) - rx040(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D_y D-y)^2 
                                                    rx410i(m1,m2) = (rx400(i1,i2+1,i3,m1,m2) - rx400(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+x D-x)^2
                                                    rx320i(m1,m2) = (rx220(i1+1,i2,i3,m1,m2) - rx220(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D+x D-x)(D_y D-y)
                                                    rx230i(m1,m2) = (rx220(i1,i2+1,i3,m1,m2) - rx220(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)
                                                    rx510i(m1,m2) = (rx400(i1+1,i2+1,i3,m1,m2) - rx400(i1-1,i2+1,i3,m1,m2) - rx400(i1+1,i2-1,i3,m1,m2) + rx400(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2
                                                    rx150i(m1,m2) = (rx040(i1+1,i2+1,i3,m1,m2) - rx040(i1-1,i2+1,i3,m1,m2) - rx040(i1+1,i2-1,i3,m1,m2) + rx040(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^2   
                                                    rx330i(m1,m2) = (rx220(i1+1,i2+1,i3,m1,m2) - rx220(i1-1,i2+1,i3,m1,m2) - rx220(i1+1,i2-1,i3,m1,m2) + rx220(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)   
                          ! rx240i(m1,m2) = rx240(i1,i2,i3,m1,m2)    
                          ! rx420i(m1,m2) = rx420(i1,i2,i3,m1,m2) 
                          ! rx700i(m1,m2) = (rx600(i1+1,i2,i3,m1,m2) - rx600(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D_x D-x)^3 
                          ! rx070i(m1,m2) = (rx060(i1,i2+1,i3,m1,m2) - rx060(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+y D-y)^3 
                          ! rx160i(m1,m2) = (rx060(i1+1,i2,i3,m1,m2) - rx060(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D_y D-y)^3 
                          ! rx610i(m1,m2) = (rx600(i1,i2+1,i3,m1,m2) - rx600(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+x D-x)^3
                          ! rx520i(m1,m2) = (rx420(i1+1,i2,i3,m1,m2) - rx420(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D+x D-x)(D_y D-y)
                          ! rx250i(m1,m2) = (rx240(i1,i2+1,i3,m1,m2) - rx240(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)      
                          ! rx340i(m1,m2) = (rx240(i1+1,i2,i3,m1,m2) - rx240(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D+x D-x)^2(D_y D-y)
                          ! rx430i(m1,m2) = (rx420(i1,i2+1,i3,m1,m2) - rx420(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)^2
                          ! rx800i(m1,m2) = (rx600(i1+1,i2,i3,m1,m2) - 2.*rx600(i1,i2,i3,m1,m2) + rx600(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                          ! rx080i(m1,m2) = (rx060(i1,i2+1,i3,m1,m2) - 2.*rx060(i1,i2,i3,m1,m2) + rx060(i1,i2-1,i3,m1,m2))/(dx(1)**2) 
                          ! rx620i(m1,m2) = (rx420(i1+1,i2,i3,m1,m2) - 2.*rx420(i1,i2,i3,m1,m2) + rx420(i1-1,i2,i3,m1,m2))/(dx(0)**2)  
                          ! rx260i(m1,m2) = (rx240(i1,i2+1,i3,m1,m2) - 2.*rx240(i1,i2,i3,m1,m2) + rx240(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                          ! rx440i(m1,m2) = (rx240(i1+1,i2,i3,m1,m2) - 2.*rx240(i1,i2,i3,m1,m2) + rx240(i1-1,i2,i3,m1,m2))/(dx(0)**2)          
                          ! rx530i(m1,m2) = (rx420(i1+1,i2+1,i3,m1,m2) - rx420(i1-1,i2+1,i3,m1,m2) - rx420(i1+1,i2-1,i3,m1,m2) + rx420(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2 (D+yD-y)         
                          ! rx350i(m1,m2) = (rx240(i1+1,i2+1,i3,m1,m2) - rx240(i1-1,i2+1,i3,m1,m2) - rx240(i1+1,i2-1,i3,m1,m2) + rx240(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)^2             
                          ! rx710i(m1,m2) = (rx600(i1+1,i2+1,i3,m1,m2) - rx600(i1-1,i2+1,i3,m1,m2) - rx600(i1+1,i2-1,i3,m1,m2) + rx600(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^3
                          ! rx170i(m1,m2) = (rx060(i1+1,i2+1,i3,m1,m2) - rx060(i1-1,i2+1,i3,m1,m2) - rx060(i1+1,i2-1,i3,m1,m2) + rx060(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^3  
                                                end do
                                            end do  
                      ! ---- first derivatives ----
                      !   1/6, 1/30, 1/140, 1/630 
                      !   ux = D0 ( 1 - (1/6)*dx^2(D+D-) + (1/30)*dx^4 (D+D-)^2 -(1/140)*(D+D-)^3 + ... )
                                            ur = d100i - (dx(0)**2/6.)*d300i + (dx(0)**4/30.)*d500i ! - (dx(0)**6/140.)*d700i 
                                            us = d010i - (dx(1)**2/6.)*d030i + (dx(1)**4/30.)*d050i ! - (dx(1)**6/140.)*d070i  
                      ! ---- 2nd derivatives ----
                      ! 1/12 , 1/90, 1/560 
                      !   uxx = D+D- ( 1 - (1/12)*dx^2(D+D-) + (1/90)*dx^4 (D+D-)^2 -(1/560)*(D+D-)^3 + ... )
                                            urr = d200(i1,i2,i3,0) - (dx(0)**2/12.)*d400(i1,i2,i3,0) + (dx(0)**4/90.)*d600i ! - (dx(0)**6/560.)*d800i
                                            uss = d020(i1,i2,i3,0) - (dx(1)**2/12.)*d040(i1,i2,i3,0) + (dx(1)**4/90.)*d060i ! - (dx(1)**6/560.)*d080i 
                      ! uxy = D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                      !       D0y*( 1 - dy^2/6 D+yD-y + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+yD-y)^3 ) u 
                                            urs = d110i - (dx(0)**2/6.)*d310i - (dx(1)**2/6.)*d130i + (dx(0)**4/30.)*d510i + (dx(1)**4/30.)*d150i + (dx(0)**2 * dx(1)**2/36. )*d330i 
                          !  - ( (1./140.)*dx(0)**6 )*d710i - ( (1./140.)*dx(1)**6 )*d170i - (dx(0)**4*dx(1)**2/(6.*30.))*d530i - (dx(0)**2*dx(1)**4/(6.*30.))*d350i 
                      ! !   uxx = D+D- ( 1 - (1/12)*dx^2(D+D-) + (1/90)*dx^4 (D+D-)^2 -(1/560)*(D+D-)^3 + ... )
                      ! urr = d200(i1,i2,i3,0) - (dx(0)**2/12.)*d400(i1,i2,i3,0)
                      ! uss = d020(i1,i2,i3,0) - (dx(1)**2/12.)*d040(i1,i2,i3,0)
                      ! !   ux = D0 ( 1 - (1/6)*dx^2(D+D-) + (1/30)*dx^4 (D+D-)^2 -(1/140)*(D+D-)^3 + ... )
                      ! urs = d110i - (dx(0)**2/6.)*d310i - (dx(1)**2/6.)*d130i
                      ! ur = d100i  - (dx(0)**2/6.)*d300i
                      ! us = d010i  - (dx(1)**2/6.)*d030i
                      ! spatial derivatives of metrics:
                      !   rxxi(m1,m2,m3) = Dx_m3 rsxy(m1,m2)
                                            do m1=0,nd-1
                                                do m2=0,nd-1
                                                    rxr(m1,m2,0) = rx100i(m1,m2) - (dx(0)**2/6.)*rx300i(m1,m2) + (dx(0)**4/30.)*rx500i(m1,m2) !- (dx(0)**6/140.)*rx700i(m1,m2)   
                                                    rxr(m1,m2,1) = rx010i(m1,m2) - (dx(1)**2/6.)*rx030i(m1,m2) + (dx(1)**4/30.)*rx050i(m1,m2) !- (dx(1)**6/140.)*rx070i(m1,m2)   
                                                    do m3 = 0,nd-1
                                                        rxxa(m1,m2,m3) = rx000i(0,m3)*rxr(m1,m2,0) + rx000i(1,m3)*rxr(m1,m2,1)  ! (rx).x = rx *(rx)_r + sx * (rx)_s 
                                                    end do
                                                end do
                                            end do
                                            rx = rsxy(i1,i2,i3,0,0)
                                            ry = rsxy(i1,i2,i3,0,1)
                                            sx = rsxy(i1,i2,i3,1,0)
                                            sy = rsxy(i1,i2,i3,1,1)
                                            rxx = rxxa(0,0,0)
                                            sxx = rxxa(1,0,0)
                                            ryy = rxxa(0,1,1)
                                            syy = rxxa(1,1,1)
                                            rxy = rxxa(0,0,1) ! or rxx(0,1,0)
                                            sxy = rxxa(1,0,1) ! or rxx(1,1,0)
                          ! ------ Coefficients in expansion for ux ------
                                                    cux100 = rx
                                                    cux010 = sx
                          ! ux = cux100*ur+cux010*us
                          ! ------ Coefficients in expansion for uy ------
                                                    cuy100 = ry
                                                    cuy010 = sy
                          ! uy = cuy100*ur+cuy010*us
                      ! ux = cux100*ur+cux010*us
                      ! uy = cuy100*ur+cuy010*us      
                          ! ------ Coefficients in expansion for uxx ------
                                                    cuxx100 = rxx
                                                    cuxx200 = rx**2
                                                    cuxx010 = sxx
                                                    cuxx110 = 2.*rx*sx
                                                    cuxx020 = sx**2
                          ! uxx = cuxx100*ur+cuxx200*urr+cuxx010*us+cuxx110*urs+cuxx020*uss
                          ! ------ Coefficients in expansion for uxy ------
                                                    cuxy100 = rxy
                                                    cuxy200 = rx*ry
                                                    cuxy010 = sxy
                                                    cuxy110 = rx*sy+ry*sx
                                                    cuxy020 = sx*sy
                          ! uxy = cuxy100*ur+cuxy200*urr+cuxy010*us+cuxy110*urs+cuxy020*uss
                          ! ------ Coefficients in expansion for uyy ------
                                                    cuyy100 = ryy
                                                    cuyy200 = ry**2
                                                    cuyy010 = syy
                                                    cuyy110 = 2.*ry*sy
                                                    cuyy020 = sy**2
                          ! uyy = cuyy100*ur+cuyy200*urr+cuyy010*us+cuyy110*urs+cuyy020*uss
                                            uxx = cuxx100*ur+cuxx200*urr+cuxx010*us+cuxx110*urs+cuxx020*uss
                      ! uxy = cuxy100*ur+cuxy200*urr+cuxy010*us+cuxy110*urs+cuxy020*uss
                                            uyy = cuyy100*ur+cuyy200*urr+cuyy010*us+cuyy110*urs+cuyy020*uss      
                      ! ----- START THIRD DERIVATIVES -----
                                            do m1=0,nd-1
                                            do m2=0,nd-1
                        ! ---- 2nd parameteric derivatives of the metrics ----
                                                rxrra(m1,m2) = rx200(i1,i2,i3,m1,m2) - (dx(0)**2/12.)*rx400(i1,i2,i3,m1,m2) + (dx(0)**4/90.)*rx600i(m1,m2) ! - (dx(0)**6/560.)*rx800i(m1,m2)
                                                rxssa(m1,m2) = rx020(i1,i2,i3,m1,m2) - (dx(1)**2/12.)*rx040(i1,i2,i3,m1,m2) + (dx(1)**4/90.)*rx060i(m1,m2) ! - (dx(1)**6/560.)*rx080i(m1,m2) 
                                                rxrsa(m1,m2) = rx110i(m1,m2) - (dx(0)**2/6.)*rx310i(m1,m2) - (dx(1)**2/6.)*rx130i(m1,m2) + (dx(0)**4/30.)*rx510i(m1,m2) + (dx(1)**4/30.)*rx150i(m1,m2) + (dx(0)**2 * dx(1)**2/36. )*rx330i(m1,m2) 
                              ! - ( (1./140.)*dx(0)**6 )*rx710i(m1,m2) - ( (1./140.)*dx(1)**6 )*rx170i(m1,m2) - (dx(0)**4*dx(1)**2/(6.*30.))*rx530i(m1,m2) - (dx(0)**2*dx(1)**4/(6.*30.))*rx350i(m1,m2)   
                        ! ---- 2nd spatial derivatives of the metrics ----
                                                rxxxa(m1,m2) = cuxx100*rxr(m1,m2,0)+cuxx200*rxrra(m1,m2)+cuxx010*rxr(m1,m2,1)+cuxx110*rxrsa(m1,m2)+cuxx020*rxssa(m1,m2)
                                                rxxya(m1,m2) = cuxy100*rxr(m1,m2,0)+cuxy200*rxrra(m1,m2)+cuxy010*rxr(m1,m2,1)+cuxy110*rxrsa(m1,m2)+cuxy020*rxssa(m1,m2)
                                                rxyya(m1,m2) = cuyy100*rxr(m1,m2,0)+cuyy200*rxrra(m1,m2)+cuyy010*rxr(m1,m2,1)+cuyy110*rxrsa(m1,m2)+cuyy020*rxssa(m1,m2)
                                            end do
                                            end do
                                            rxxx=rxxxa(0,0); ryxx=rxxxa(0,1); sxxx=rxxxa(1,0); syxx=rxxxa(1,1); 
                                            rxxy=rxxya(0,0); ryxy=rxxya(0,1); sxxy=rxxya(1,0); syxy=rxxya(1,1); 
                                            rxyy=rxyya(0,0); ryyy=rxyya(0,1); sxyy=rxyya(1,0); syyy=rxyya(1,1); 
                      ! ---- third parametric derivatives ----
                      !   1/4, 7/120, 41/3024
                      !    uxxxx = D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ... 
                                            urrr = d300i - (dx(0)**2/4.)*d500i ! + (dx(0)**4 *7./120.)*d700i
                                            usss = d030i - (dx(1)**2/4.)*d050i ! + (dx(1)**4 *7./120.)*d070i
                      ! uxxy = D+xD-x*( 1 - dx^2/12 D+xD-x + dx^4/90*(D+xD-x)^2 - (1/560)*dx^6*(D+xD-x)^3 ) X
                      !           D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                                            urrs = d210i - (dx(0)**2/12.)*d410i - (dx(1)**2/6.)*d230i ! + (dx(0)**4/90.)*d610i  + (dx(1)**4/30.)*d250i + (dx(0)**2 * dx(1)**2/72. )*d430i
                      ! uxyy =    D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                      !        D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                            urss = d120i - (dx(1)**2/12.)*d140i - (dx(0)**2/6.)*d320i ! + (dx(1)**4/90.)*d160i  + (dx(0)**4/30.)*d520i + (dx(0)**2 * dx(1)**2/72. )*d340i 
                      ! ----- THIRD SPATIAL DERIVATIVES -----
                      ! ------ Coefficients in expansion for uxxx ------
cuxxx100 = rxxx
cuxxx200 = 3.*rx*rxx
cuxxx300 = rx**3
cuxxx010 = sxxx
cuxxx110 = 3.*rx*sxx+3.*rxx*sx
cuxxx210 = 3.*rx**2.*sx
cuxxx020 = 3.*sx*sxx
cuxxx120 = 3.*rx*sx**2
cuxxx030 = sx**3
                      ! uxxx = cuxxx100*ur+cuxxx200*urr+cuxxx300*urrr+cuxxx010*us+cuxxx110*urs+cuxxx210*urrs+cuxxx020*uss+cuxxx120*urss+cuxxx030*usss
                      ! ------ Coefficients in expansion for uxxy ------
cuxxy100 = rxxy
cuxxy200 = 2.*rx*rxy+rxx*ry
cuxxy300 = rx**2.*ry
cuxxy010 = sxxy
cuxxy110 = 2.*rx*sxy+rxx*sy+2.*rxy*sx+ry*sxx
cuxxy210 = rx**2.*sy+2.*rx*ry*sx
cuxxy020 = 2.*sx*sxy+sxx*sy
cuxxy120 = 2.*rx*sx*sy+ry*sx**2
cuxxy030 = sx**2.*sy
                      ! uxxy = cuxxy100*ur+cuxxy200*urr+cuxxy300*urrr+cuxxy010*us+cuxxy110*urs+cuxxy210*urrs+cuxxy020*uss+cuxxy120*urss+cuxxy030*usss
                      ! ------ Coefficients in expansion for uxyy ------
cuxyy100 = rxyy
cuxyy200 = rx*ryy+2.*rxy*ry
cuxyy300 = rx*ry**2
cuxyy010 = sxyy
cuxyy110 = rx*syy+2.*rxy*sy+2.*ry*sxy+ryy*sx
cuxyy210 = 2.*rx*ry*sy+ry**2.*sx
cuxyy020 = sx*syy+2.*sxy*sy
cuxyy120 = rx*sy**2+2.*ry*sx*sy
cuxyy030 = sx*sy**2
                      ! uxyy = cuxyy100*ur+cuxyy200*urr+cuxyy300*urrr+cuxyy010*us+cuxyy110*urs+cuxyy210*urrs+cuxyy020*uss+cuxyy120*urss+cuxyy030*usss
                      ! ------ Coefficients in expansion for uyyy ------
cuyyy100 = ryyy
cuyyy200 = 3.*ry*ryy
cuyyy300 = ry**3
cuyyy010 = syyy
cuyyy110 = 3.*ry*syy+3.*ryy*sy
cuyyy210 = 3.*ry**2.*sy
cuyyy020 = 3.*sy*syy
cuyyy120 = 3.*ry*sy**2
cuyyy030 = sy**3
                      ! uyyy = cuyyy100*ur+cuyyy200*urr+cuyyy300*urrr+cuyyy010*us+cuyyy110*urs+cuyyy210*urrs+cuyyy020*uss+cuyyy120*urss+cuyyy030*usss
                      ! uxxx = cuxxx100*ur+cuxxx200*urr+cuxxx300*urrr+cuxxx010*us+cuxxx110*urs+cuxxx210*urrs+cuxxx020*uss+cuxxx120*urss+cuxxx030*usss
                      ! uxxy = cuxxy100*ur+cuxxy200*urr+cuxxy300*urrr+cuxxy010*us+cuxxy110*urs+cuxxy210*urrs+cuxxy020*uss+cuxxy120*urss+cuxxy030*usss
                      ! uxyy = cuxyy100*ur+cuxyy200*urr+cuxyy300*urrr+cuxyy010*us+cuxyy110*urs+cuxyy210*urrs+cuxyy020*uss+cuxyy120*urss+cuxyy030*usss
                      ! uyyy = cuyyy100*ur+cuyyy200*urr+cuyyy300*urrr+cuyyy010*us+cuyyy110*urs+cuyyy210*urrs+cuyyy020*uss+cuyyy120*urss+cuyyy030*usss
                      ! ----- START FOURTH DERIVATIVES -----
                                            do m1=0,nd-1
                                            do m2=0,nd-1
                        ! ---- 3rd parameteric derivatives of the metrics ----
                                                rxrrra(m1,m2) = rx300i(m1,m2) - (dx(0)**2/4.)*rx500i(m1,m2)  ! + (dx(0)**4 *7./120.)*rx700i(m1,m2)
                                                rxsssa(m1,m2) = rx030i(m1,m2) - (dx(1)**2/4.)*rx050i(m1,m2)  ! + (dx(1)**4 *7./120.)*rx070i(m1,m2)
                                                rxrrsa(m1,m2) = rx210i(m1,m2) - (dx(0)**2/12.)*rx410i(m1,m2) - (dx(1)**2/6.)*rx230i(m1,m2) ! + (dx(0)**4/90.)*rx610i(m1,m2) + (dx(1)**4/30.)*rx250i(m1,m2) + (dx(0)**2 * dx(1)**2/72. )*rx430i(m1,m2)
                                                rxrssa(m1,m2) = rx120i(m1,m2) - (dx(1)**2/12.)*rx140i(m1,m2) - (dx(0)**2/6.)*rx320i(m1,m2) ! + (dx(1)**4/90.)*rx160i(m1,m2) + (dx(0)**4/30.)*rx520i(m1,m2) + (dx(0)**2 * dx(1)**2/72. )*rx340i(m1,m2) 
                        ! ---- third spatial derivatives of the metrics ----
                                                rxxxxa(m1,m2) = cuxxx100*rxr(m1,m2,0)+cuxxx200*rxrra(m1,m2)+cuxxx300*rxrrra(m1,m2)+cuxxx010*rxr(m1,m2,1)+cuxxx110*rxrsa(m1,m2)+cuxxx210*rxrrsa(m1,m2)+cuxxx020*rxssa(m1,m2)+cuxxx120*rxrssa(m1,m2)+cuxxx030*rxsssa(m1,m2)
                                                rxxxya(m1,m2) = cuxxy100*rxr(m1,m2,0)+cuxxy200*rxrra(m1,m2)+cuxxy300*rxrrra(m1,m2)+cuxxy010*rxr(m1,m2,1)+cuxxy110*rxrsa(m1,m2)+cuxxy210*rxrrsa(m1,m2)+cuxxy020*rxssa(m1,m2)+cuxxy120*rxrssa(m1,m2)+cuxxy030*rxsssa(m1,m2)
                                                rxxyya(m1,m2) = cuxyy100*rxr(m1,m2,0)+cuxyy200*rxrra(m1,m2)+cuxyy300*rxrrra(m1,m2)+cuxyy010*rxr(m1,m2,1)+cuxyy110*rxrsa(m1,m2)+cuxyy210*rxrrsa(m1,m2)+cuxyy020*rxssa(m1,m2)+cuxyy120*rxrssa(m1,m2)+cuxyy030*rxsssa(m1,m2)
                                                rxyyya(m1,m2) = cuyyy100*rxr(m1,m2,0)+cuyyy200*rxrra(m1,m2)+cuyyy300*rxrrra(m1,m2)+cuyyy010*rxr(m1,m2,1)+cuyyy110*rxrsa(m1,m2)+cuyyy210*rxrrsa(m1,m2)+cuyyy020*rxssa(m1,m2)+cuyyy120*rxrssa(m1,m2)+cuyyy030*rxsssa(m1,m2)
                                            end do
                                            end do
                                            rxxxx=rxxxxa(0,0); ryxxx=rxxxxa(0,1); sxxxx=rxxxxa(1,0); syxxx=rxxxxa(1,1); 
                                            rxxxy=rxxxya(0,0); ryxxy=rxxxya(0,1); sxxxy=rxxxya(1,0); syxxy=rxxxya(1,1); 
                                            rxxyy=rxxyya(0,0); ryxyy=rxxyya(0,1); sxxyy=rxxyya(1,0); syxyy=rxxyya(1,1); 
                                            rxyyy=rxyyya(0,0); ryyyy=rxyyya(0,1); sxyyy=rxyyya(1,0); syyyy=rxyyya(1,1); 
                      ! ----- fourth parametric derivs ------
                      !  1/6, 7/240, 41/7560
                      !   uxxxx = (D+xD-x)^2*( 1 - (1/6)*D+xD-x + (7/240)*(D+xD-x)^2 + ... )
                                            urrrr = d400(i1,i2,i3,0) - (dx(0)**2/6.)*d600i  !  + (dx(0)**4 * (7./240))*d800i
                                            ussss = d040(i1,i2,i3,0) - (dx(1)**2/6.)*d060i  !  + (dx(1)**4 * (7./240))*d080i
                      ! uxyyy = D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                      !         D0y D+yD-y*( 1 - (1/4)*dy^2 * D+yD-y + (7/120)*dy^4 (D+yD-y)^2 + ... 
                                            ursss = d130i - (dx(0)**2/6.)*d330i - (dx(1)**2/4.)*d150i ! + (dx(0)**4/30.)*d530i + (dx(1)**4*7./120.)*d170i + (dx(0)**2 * dx(1)**2/24. )*d350i 
                      ! uxxxy
                                            urrrs = d310i - (dx(1)**2/6.)*d330i - (dx(0)**2/4.)*d510i ! + (dx(1)**4/30.)*d350i + (dx(0)**4*7./120.)*d710i + (dx(0)**2 * dx(1)**2/24. )*d530i       
                      ! uxxyy = D+xD-x*( 1 - dx^2/12 D+xD-x + dx^4/90*(D+xD-x)^2 - (1/560)*dx^6*(D+xD-x)^3 ) X 
                      !         D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                            urrss = d220(i1,i2,i3,0) - (dx(0)**2/12.)*d420i - (dx(1)**2/12.)*d240i 
                                            !  (dx(0)**4 * (1./90))*d620i + (dx(1)**4 * (1./90))*d260i + (dx(0)**2 * dx(1)**2 * 1./144.)*d440i 
                      ! ----- FOURTH SPATIAL DERIVATIVES -----
                      ! ------ Coefficients in expansion for uxxxx ------
cuxxxx100 = rxxxx
cuxxxx200 = 4.*rx*rxxx+3.*rxx**2
cuxxxx300 = 6.*rx**2.*rxx
cuxxxx400 = rx**4
cuxxxx010 = sxxxx
cuxxxx110 = 4.*rx*sxxx+6.*rxx*sxx+4.*rxxx*sx
cuxxxx210 = 6.*rx**2.*sxx+12.*rx*rxx*sx
cuxxxx310 = 4.*rx**3.*sx
cuxxxx020 = 4.*sx*sxxx+3.*sxx**2
cuxxxx120 = 12.*rx*sx*sxx+6.*rxx*sx**2
cuxxxx220 = 6.*rx**2.*sx**2
cuxxxx030 = 6.*sx**2.*sxx
cuxxxx130 = 4.*rx*sx**3
cuxxxx040 = sx**4
                      !uxxxx = cuxxxx100*ur+cuxxxx200*urr+cuxxxx300*urrr+cuxxxx400*urrrr+cuxxxx010*us+cuxxxx110*urs+cuxxxx210*urrs+cuxxxx310*urrrs+cuxxxx020*uss+cuxxxx120*urss+cuxxxx220*urrss+cuxxxx030*usss+cuxxxx130*ursss+cuxxxx040*ussss
                      ! ------ Coefficients in expansion for uxxxy ------
cuxxxy100 = rxxxy
cuxxxy200 = 3.*rx*rxxy+3.*rxx*rxy+rxxx*ry
cuxxxy300 = 3.*rx**2.*rxy+3.*rx*rxx*ry
cuxxxy400 = rx**3.*ry
cuxxxy010 = sxxxy
cuxxxy110 = 3.*rx*sxxy+3.*rxx*sxy+rxxx*sy+3.*rxxy*sx+3.*rxy*sxx+ry*sxxx
cuxxxy210 = 3.*rx**2.*sxy+(3.*rxx*sy+6.*rxy*sx+3.*ry*sxx)*rx+3.*rxx*ry*sx
cuxxxy310 = rx**3.*sy+3.*rx**2.*ry*sx
cuxxxy020 = 3.*sx*sxxy+3.*sxx*sxy+sxxx*sy
cuxxxy120 = 3.*rxy*sx**2+(6.*rx*sxy+3.*rxx*sy+3.*ry*sxx)*sx+3.*rx*sxx*sy
cuxxxy220 = 3.*rx**2.*sx*sy+3.*rx*ry*sx**2
cuxxxy030 = 3.*sx**2.*sxy+3.*sx*sxx*sy
cuxxxy130 = 3.*rx*sx**2.*sy+ry*sx**3
cuxxxy040 = sx**3.*sy
                      !uxxxy = cuxxxy100*ur+cuxxxy200*urr+cuxxxy300*urrr+cuxxxy400*urrrr+cuxxxy010*us+cuxxxy110*urs+cuxxxy210*urrs+cuxxxy310*urrrs+cuxxxy020*uss+cuxxxy120*urss+cuxxxy220*urrss+cuxxxy030*usss+cuxxxy130*ursss+cuxxxy040*ussss
                      ! ------ Coefficients in expansion for uxxyy ------
cuxxyy100 = rxxyy
cuxxyy200 = 2.*rx*rxyy+rxx*ryy+2.*rxxy*ry+2.*rxy**2
cuxxyy300 = rx**2.*ryy+4.*rx*rxy*ry+rxx*ry**2
cuxxyy400 = rx**2.*ry**2
cuxxyy010 = sxxyy
cuxxyy110 = 2.*rx*sxyy+rxx*syy+2.*rxxy*sy+4.*rxy*sxy+2.*rxyy*sx+2.*ry*sxxy+ryy*sxx
cuxxyy210 = rx**2.*syy+(4.*rxy*sy+4.*ry*sxy+2.*ryy*sx)*rx+ry*(2.*rxx*sy+4.*rxy*sx+ry*sxx)
cuxxyy310 = 2.*rx**2.*ry*sy+2.*rx*ry**2.*sx
cuxxyy020 = 2.*sx*sxyy+sxx*syy+2.*sxxy*sy+2.*sxy**2
cuxxyy120 = ryy*sx**2+(2.*rx*syy+4.*rxy*sy+4.*ry*sxy)*sx+rxx*sy**2+(4.*rx*sxy+2.*ry*sxx)*sy
cuxxyy220 = rx**2.*sy**2+4.*rx*ry*sx*sy+ry**2.*sx**2
cuxxyy030 = sx**2.*syy+4.*sx*sxy*sy+sxx*sy**2
cuxxyy130 = 2.*rx*sx*sy**2+2.*ry*sx**2.*sy
cuxxyy040 = sx**2.*sy**2
                      !uxxyy = cuxxyy100*ur+cuxxyy200*urr+cuxxyy300*urrr+cuxxyy400*urrrr+cuxxyy010*us+cuxxyy110*urs+cuxxyy210*urrs+cuxxyy310*urrrs+cuxxyy020*uss+cuxxyy120*urss+cuxxyy220*urrss+cuxxyy030*usss+cuxxyy130*ursss+cuxxyy040*ussss
                      ! ------ Coefficients in expansion for uxyyy ------
cuxyyy100 = rxyyy
cuxyyy200 = rx*ryyy+3.*rxy*ryy+3.*rxyy*ry
cuxyyy300 = 3.*rx*ry*ryy+3.*rxy*ry**2
cuxyyy400 = rx*ry**3
cuxyyy010 = sxyyy
cuxyyy110 = rx*syyy+3.*rxy*syy+3.*rxyy*sy+3.*ry*sxyy+3.*ryy*sxy+ryyy*sx
cuxyyy210 = 3.*ry**2.*sxy+(3.*rx*syy+6.*rxy*sy+3.*ryy*sx)*ry+3.*ryy*rx*sy
cuxyyy310 = 3.*rx*ry**2.*sy+ry**3.*sx
cuxyyy020 = sx*syyy+3.*sxy*syy+3.*sxyy*sy
cuxyyy120 = 3.*rxy*sy**2+(3.*rx*syy+6.*ry*sxy+3.*ryy*sx)*sy+3.*ry*sx*syy
cuxyyy220 = 3.*rx*ry*sy**2+3.*ry**2.*sx*sy
cuxyyy030 = 3.*sx*sy*syy+3.*sxy*sy**2
cuxyyy130 = rx*sy**3+3.*ry*sx*sy**2
cuxyyy040 = sx*sy**3
                      !uxyyy = cuxyyy100*ur+cuxyyy200*urr+cuxyyy300*urrr+cuxyyy400*urrrr+cuxyyy010*us+cuxyyy110*urs+cuxyyy210*urrs+cuxyyy310*urrrs+cuxyyy020*uss+cuxyyy120*urss+cuxyyy220*urrss+cuxyyy030*usss+cuxyyy130*ursss+cuxyyy040*ussss
                      ! ------ Coefficients in expansion for uyyyy ------
cuyyyy100 = ryyyy
cuyyyy200 = 4.*ry*ryyy+3.*ryy**2
cuyyyy300 = 6.*ry**2.*ryy
cuyyyy400 = ry**4
cuyyyy010 = syyyy
cuyyyy110 = 4.*ry*syyy+6.*ryy*syy+4.*ryyy*sy
cuyyyy210 = 6.*ry**2.*syy+12.*ry*ryy*sy
cuyyyy310 = 4.*ry**3.*sy
cuyyyy020 = 4.*sy*syyy+3.*syy**2
cuyyyy120 = 12.*ry*sy*syy+6.*ryy*sy**2
cuyyyy220 = 6.*ry**2.*sy**2
cuyyyy030 = 6.*sy**2.*syy
cuyyyy130 = 4.*ry*sy**3
cuyyyy040 = sy**4
                      !uyyyy = cuyyyy100*ur+cuyyyy200*urr+cuyyyy300*urrr+cuyyyy400*urrrr+cuyyyy010*us+cuyyyy110*urs+cuyyyy210*urrs+cuyyyy310*urrrs+cuyyyy020*uss+cuyyyy120*urss+cuyyyy220*urrss+cuyyyy030*usss+cuyyyy130*ursss+cuyyyy040*ussss
                                            uxxxx = cuxxxx100*ur+cuxxxx200*urr+cuxxxx300*urrr+cuxxxx400*urrrr+cuxxxx010*us+cuxxxx110*urs+cuxxxx210*urrs+cuxxxx310*urrrs+cuxxxx020*uss+cuxxxx120*urss+cuxxxx220*urrss+cuxxxx030*usss+cuxxxx130*ursss+cuxxxx040*ussss
                      ! uxxxy = cuxxxy100*ur+cuxxxy200*urr+cuxxxy300*urrr+cuxxxy400*urrrr+cuxxxy010*us+cuxxxy110*urs+cuxxxy210*urrs+cuxxxy310*urrrs+cuxxxy020*uss+cuxxxy120*urss+cuxxxy220*urrss+cuxxxy030*usss+cuxxxy130*ursss+cuxxxy040*ussss
                                            uxxyy = cuxxyy100*ur+cuxxyy200*urr+cuxxyy300*urrr+cuxxyy400*urrrr+cuxxyy010*us+cuxxyy110*urs+cuxxyy210*urrs+cuxxyy310*urrrs+cuxxyy020*uss+cuxxyy120*urss+cuxxyy220*urrss+cuxxyy030*usss+cuxxyy130*ursss+cuxxyy040*ussss
                      ! uxyyy = cuxyyy100*ur+cuxyyy200*urr+cuxyyy300*urrr+cuxyyy400*urrrr+cuxyyy010*us+cuxyyy110*urs+cuxyyy210*urrs+cuxyyy310*urrrs+cuxyyy020*uss+cuxyyy120*urss+cuxyyy220*urrss+cuxyyy030*usss+cuxyyy130*ursss+cuxyyy040*ussss
                                            uyyyy = cuyyyy100*ur+cuyyyy200*urr+cuyyyy300*urrr+cuyyyy400*urrrr+cuyyyy010*us+cuyyyy110*urs+cuyyyy210*urrs+cuyyyy310*urrrs+cuyyyy020*uss+cuyyyy120*urss+cuyyyy220*urrss+cuyyyy030*usss+cuyyyy130*ursss+cuyyyy040*ussss
                      ! ----- START FIFTH DERIVATIVES -----
                                            do m1=0,nd-1
                                            do m2=0,nd-1
                        ! ---- 4th parameteric derivatives of the metrics ----
                                                rxrrrra(m1,m2) = rx400(i1,i2,i3,m1,m2) - (dx(0)**2/6.)*rx600i(m1,m2) ! + (dx(0)**4 * (7./240))*rx800i(m1,m2)
                                                rxssssa(m1,m2) = rx040(i1,i2,i3,m1,m2) - (dx(1)**2/6.)*rx060i(m1,m2) ! + (dx(1)**4 * (7./240))*rx080i(m1,m2)
                                                rxrsssa(m1,m2) = rx130i(m1,m2) - (dx(0)**2/6.)*rx330i(m1,m2) - (dx(1)**2/4.)*rx150i(m1,m2)! + (dx(0)**4/30.)*rx530i(m1,m2) + (dx(1)**4*7./120.)*rx170i(m1,m2) + (dx(0)**2 * dx(1)**2/24. )*rx350i(m1,m2) 
                                                rxrrrsa(m1,m2) = rx310i(m1,m2) - (dx(1)**2/6.)*rx330i(m1,m2) - (dx(0)**2/4.)*rx510i(m1,m2)! + (dx(1)**4/30.)*rx350i(m1,m2) + (dx(0)**4*7./120.)*rx710i(m1,m2) + (dx(0)**2 * dx(1)**2/24. )*rx530i(m1,m2)       
                                                rxrrssa(m1,m2) = rx220(i1,i2,i3,m1,m2) - (dx(0)**2/12.)*rx420i(m1,m2) - (dx(1)**2/12.)*rx240i(m1,m2) 
                                             ! (dx(0)**4 * (1./90))*rx620i(m1,m2) + (dx(1)**4 * (1./90))*rx260i(m1,m2) + (dx(0)**2 * dx(1)**2 * 1./144.)*rx440i(m1,m2)         
                        ! ---- fourth spatial derivatives of the metrics ----
                                                rxxxxxa(m1,m2) = cuxxxx100*rxr(m1,m2,0)+cuxxxx200*rxrra(m1,m2)+cuxxxx300*rxrrra(m1,m2)+cuxxxx400*rxrrrra(m1,m2)+cuxxxx010*rxr(m1,m2,1)+cuxxxx110*rxrsa(m1,m2)+cuxxxx210*rxrrsa(m1,m2)+cuxxxx310*rxrrrsa(m1,m2)+cuxxxx020*rxssa(m1,m2)+cuxxxx120*rxrssa(m1,m2)+cuxxxx220*rxrrssa(m1,m2)+cuxxxx030*rxsssa(m1,m2)+cuxxxx130*rxrsssa(m1,m2)+cuxxxx040*rxssssa(m1,m2)
                                                rxxxxya(m1,m2) = cuxxxy100*rxr(m1,m2,0)+cuxxxy200*rxrra(m1,m2)+cuxxxy300*rxrrra(m1,m2)+cuxxxy400*rxrrrra(m1,m2)+cuxxxy010*rxr(m1,m2,1)+cuxxxy110*rxrsa(m1,m2)+cuxxxy210*rxrrsa(m1,m2)+cuxxxy310*rxrrrsa(m1,m2)+cuxxxy020*rxssa(m1,m2)+cuxxxy120*rxrssa(m1,m2)+cuxxxy220*rxrrssa(m1,m2)+cuxxxy030*rxsssa(m1,m2)+cuxxxy130*rxrsssa(m1,m2)+cuxxxy040*rxssssa(m1,m2)
                                                rxxxyya(m1,m2) = cuxxyy100*rxr(m1,m2,0)+cuxxyy200*rxrra(m1,m2)+cuxxyy300*rxrrra(m1,m2)+cuxxyy400*rxrrrra(m1,m2)+cuxxyy010*rxr(m1,m2,1)+cuxxyy110*rxrsa(m1,m2)+cuxxyy210*rxrrsa(m1,m2)+cuxxyy310*rxrrrsa(m1,m2)+cuxxyy020*rxssa(m1,m2)+cuxxyy120*rxrssa(m1,m2)+cuxxyy220*rxrrssa(m1,m2)+cuxxyy030*rxsssa(m1,m2)+cuxxyy130*rxrsssa(m1,m2)+cuxxyy040*rxssssa(m1,m2)
                                                rxxyyya(m1,m2) = cuxyyy100*rxr(m1,m2,0)+cuxyyy200*rxrra(m1,m2)+cuxyyy300*rxrrra(m1,m2)+cuxyyy400*rxrrrra(m1,m2)+cuxyyy010*rxr(m1,m2,1)+cuxyyy110*rxrsa(m1,m2)+cuxyyy210*rxrrsa(m1,m2)+cuxyyy310*rxrrrsa(m1,m2)+cuxyyy020*rxssa(m1,m2)+cuxyyy120*rxrssa(m1,m2)+cuxyyy220*rxrrssa(m1,m2)+cuxyyy030*rxsssa(m1,m2)+cuxyyy130*rxrsssa(m1,m2)+cuxyyy040*rxssssa(m1,m2)
                                                rxyyyya(m1,m2) = cuyyyy100*rxr(m1,m2,0)+cuyyyy200*rxrra(m1,m2)+cuyyyy300*rxrrra(m1,m2)+cuyyyy400*rxrrrra(m1,m2)+cuyyyy010*rxr(m1,m2,1)+cuyyyy110*rxrsa(m1,m2)+cuyyyy210*rxrrsa(m1,m2)+cuyyyy310*rxrrrsa(m1,m2)+cuyyyy020*rxssa(m1,m2)+cuyyyy120*rxrssa(m1,m2)+cuyyyy220*rxrrssa(m1,m2)+cuyyyy030*rxsssa(m1,m2)+cuyyyy130*rxrsssa(m1,m2)+cuyyyy040*rxssssa(m1,m2)
                                            end do
                                            end do
                                            rxxxxx=rxxxxxa(0,0); ryxxxx=rxxxxxa(0,1); sxxxxx=rxxxxxa(1,0); syxxxx=rxxxxxa(1,1); 
                                            rxxxxy=rxxxxya(0,0); ryxxxy=rxxxxya(0,1); sxxxxy=rxxxxya(1,0); syxxxy=rxxxxya(1,1); 
                                            rxxxyy=rxxxyya(0,0); ryxxyy=rxxxyya(0,1); sxxxyy=rxxxyya(1,0); syxxyy=rxxxyya(1,1); 
                                            rxxyyy=rxxyyya(0,0); ryxyyy=rxxyyya(0,1); sxxyyy=rxxyyya(1,0); syxyyy=rxxyyya(1,1); 
                                            rxyyyy=rxyyyya(0,0); ryyyyy=rxyyyya(0,1); sxyyyy=rxyyyya(1,0); syyyyy=rxyyyya(1,1); 
                      ! ---- FIXTH parametric derivatives ----
                      !   1/3, 13/144, 139/6048
                      !  uxxxxx = D0x(D+xD-x)^2 *( 1 - (1/3)*dx^2 D+xD-x + ... )
                                            urrrrr = d500i ! - (dx(0)**2/3.)*d700i
                                            usssss = d050i ! - (dx(1)**2/3.)*d070i 
                      ! uxxxxy = (D+xD-x)^2 *( 1 - dx^2/6 D+xD-x + ...) X
                      !                  D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                                            urrrrs = d410i ! - (dx(0)**2/6.)*d610i - (dx(1)**2/6.)*d430i
                                            urssss = d140i ! - (dx(1)**2/6.)*d160i - (dx(0)**2/6.)*d340i    
                      ! uxxxyy =  D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ...
                      !               D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                            urrrss = d320i ! - (dx(0)**2/4.)*d520i - (dx(1)**2/12.)*d340i
                                            urrsss = d230i ! - (dx(1)**2/4.)*d250i - (dx(0)**2/12.)*d430i
                       ! ----- FIXTH SPATIAL DERIVATIVES -----
                      ! ------ Coefficients in expansion for uxxxxx ------
cuxxxxx100 = rxxxxx
cuxxxxx200 = 5.*rx*rxxxx+10.*rxx*rxxx
cuxxxxx300 = 10.*rx**2.*rxxx+15.*rx*rxx**2
cuxxxxx400 = 10.*rx**3.*rxx
cuxxxxx500 = rx**5
cuxxxxx010 = sxxxxx
cuxxxxx110 = 5.*rx*sxxxx+10.*rxx*sxxx+10.*rxxx*sxx+5.*rxxxx*sx
cuxxxxx210 = 10.*rx**2.*sxxx+30.*rx*rxx*sxx+20.*rx*rxxx*sx+15.*rxx**2.*sx
cuxxxxx310 = 10.*rx**3.*sxx+30.*rx**2.*rxx*sx
cuxxxxx410 = 5.*sx*rx**4
cuxxxxx020 = 5.*sx*sxxxx+10.*sxx*sxxx
cuxxxxx120 = 20.*rx*sx*sxxx+15.*rx*sxx**2+30.*rxx*sx*sxx+10.*rxxx*sx**2
cuxxxxx220 = 30.*rx**2.*sx*sxx+30.*rx*rxx*sx**2
cuxxxxx320 = 10.*rx**3.*sx**2
cuxxxxx030 = 10.*sx**2.*sxxx+15.*sx*sxx**2
cuxxxxx130 = 30.*rx*sx**2.*sxx+10.*rxx*sx**3
cuxxxxx230 = 10.*rx**2.*sx**3
cuxxxxx040 = 10.*sx**3.*sxx
cuxxxxx140 = 5.*rx*sx**4
cuxxxxx050 = sx**5
                      ! uxxxxx = cuxxxxx100*ur+cuxxxxx200*urr+cuxxxxx300*urrr+cuxxxxx400*urrrr+cuxxxxx500*urrrrr+cuxxxxx010*us+cuxxxxx110*urs+cuxxxxx210*urrs+cuxxxxx310*urrrs+cuxxxxx410*urrrrs+cuxxxxx020*uss+cuxxxxx120*urss+cuxxxxx220*urrss+cuxxxxx320*urrrss+cuxxxxx030*usss+cuxxxxx130*ursss+cuxxxxx230*urrsss+cuxxxxx040*ussss+cuxxxxx140*urssss+cuxxxxx050*usssss
                      ! ------ Coefficients in expansion for uxxxxy ------
cuxxxxy100 = rxxxxy
cuxxxxy200 = 4.*rx*rxxxy+6.*rxx*rxxy+4.*rxxx*rxy+rxxxx*ry
cuxxxxy300 = 6.*rx**2.*rxxy+12.*rx*rxx*rxy+4.*rx*rxxx*ry+3.*rxx**2.*ry
cuxxxxy400 = 4.*(rxy*rx+3/2.*rxx*ry)*rx**2
cuxxxxy500 = ry*rx**4
cuxxxxy010 = sxxxxy
cuxxxxy110 = 4.*rx*sxxxy+6.*rxx*sxxy+4.*rxxx*sxy+rxxxx*sy+4.*rxxxy*sx+6.*rxxy*sxx+4.*rxy*sxxx+ry*sxxxx
cuxxxxy210 = 6.*rx**2.*sxxy+(12.*rxx*sxy+4.*rxxx*sy+12.*rxxy*sx+12.*rxy*sxx+4.*ry*sxxx)*rx+3.*rxx**2.*sy+(12.*rxy*sx+6.*ry*sxx)*rxx+4.*rxxx*ry*sx
cuxxxxy310 = 4.*rx**3.*sxy+6.*rx**2.*rxx*sy+12.*rx**2.*rxy*sx+6.*rx**2.*ry*sxx+12.*rx*rxx*ry*sx
cuxxxxy410 = rx**4.*sy+4.*rx**3.*ry*sx
cuxxxxy020 = 4.*sx*sxxxy+6.*sxx*sxxy+4.*sxxx*sxy+sxxxx*sy
cuxxxxy120 = 6.*rxxy*sx**2+(12.*rx*sxxy+12.*rxx*sxy+4.*rxxx*sy+12.*rxy*sxx+4.*ry*sxxx)*sx+3.*sxx**2.*ry+(12.*rx*sxy+6.*rxx*sy)*sxx+4.*rx*sxxx*sy
cuxxxxy220 = (12.*sx*sxy+6.*sxx*sy)*rx**2+12.*sx*(rxx*sy+rxy*sx+ry*sxx)*rx+6.*rxx*ry*sx**2
cuxxxxy320 = 4.*(sy*rx+3/2.*ry*sx)*rx**2.*sx
cuxxxxy030 = 6.*sx**2.*sxxy+12.*sx*sxx*sxy+4.*sx*sxxx*sy+3.*sxx**2.*sy
cuxxxxy130 = 12.*rx*sx**2.*sxy+12.*rx*sx*sxx*sy+6.*rxx*sx**2.*sy+4.*rxy*sx**3+6.*ry*sx**2.*sxx
cuxxxxy230 = 6.*rx*sx**2.*(sy*rx+2/3.*ry*sx)
cuxxxxy040 = 4.*(sx*sxy+3/2.*sxx*sy)*sx**2
cuxxxxy140 = 4.*rx*sx**3.*sy+ry*sx**4
cuxxxxy050 = sy*sx**4
                      ! uxxxxy = cuxxxxy100*ur+cuxxxxy200*urr+cuxxxxy300*urrr+cuxxxxy400*urrrr+cuxxxxy500*urrrrr+cuxxxxy010*us+cuxxxxy110*urs+cuxxxxy210*urrs+cuxxxxy310*urrrs+cuxxxxy410*urrrrs+cuxxxxy020*uss+cuxxxxy120*urss+cuxxxxy220*urrss+cuxxxxy320*urrrss+cuxxxxy030*usss+cuxxxxy130*ursss+cuxxxxy230*urrsss+cuxxxxy040*ussss+cuxxxxy140*urssss+cuxxxxy050*usssss
                      ! ------ Coefficients in expansion for uxxxyy ------
cuxxxyy100 = rxxxyy
cuxxxyy200 = 3.*rx*rxxyy+3.*rxx*rxyy+rxxx*ryy+2.*rxxxy*ry+6.*rxxy*rxy
cuxxxyy300 = 3.*rxyy*rx**2+(3.*rxx*ryy+6.*rxxy*ry+6.*rxy**2)*rx+rxxx*ry**2+6.*rxy*rxx*ry
cuxxxyy400 = rx**3.*ryy+6.*rx**2.*rxy*ry+3.*rx*rxx*ry**2
cuxxxyy500 = rx**3.*ry**2
cuxxxyy010 = sxxxyy
cuxxxyy110 = 3.*rx*sxxyy+3.*rxx*sxyy+rxxx*syy+2.*rxxxy*sy+6.*rxxy*sxy+3.*rxxyy*sx+6.*rxy*sxxy+3.*rxyy*sxx+2.*ry*sxxxy+ryy*sxxx
cuxxxyy210 = 3.*rx**2.*sxyy+(3.*rxx*syy+6.*rxxy*sy+12.*rxy*sxy+6.*rxyy*sx+6.*ry*sxxy+3.*ryy*sxx)*rx+ry**2.*sxxx+(6.*rxx*sxy+2.*rxxx*sy+6.*rxxy*sx+6.*rxy*sxx)*ry+(6.*rxy*sy+3.*ryy*sx)*rxx+6.*rxy**2.*sx
cuxxxyy310 = rx**3.*syy+(6.*rxy*sy+6.*ry*sxy+3.*ryy*sx)*rx**2+3.*ry*(2.*rxx*sy+4.*rxy*sx+ry*sxx)*rx+3.*rxx*ry**2.*sx
cuxxxyy410 = 2.*rx**3.*ry*sy+3.*rx**2.*ry**2.*sx
cuxxxyy020 = 3.*sx*sxxyy+3.*sxx*sxyy+sxxx*syy+2.*sxxxy*sy+6.*sxxy*sxy
cuxxxyy120 = 3.*rxyy*sx**2+(6.*rx*sxyy+3.*rxx*syy+6.*rxxy*sy+12.*rxy*sxy+6.*ry*sxxy+3.*ryy*sxx)*sx+rxxx*sy**2+(6.*rx*sxxy+6.*rxx*sxy+6.*rxy*sxx+2.*ry*sxxx)*sy+(3.*rx*syy+6.*ry*sxy)*sxx+6.*rx*sxy**2
cuxxxyy220 = (3.*sx*syy+6.*sxy*sy)*rx**2+(3.*ryy*sx**2+(12.*rxy*sy+12.*ry*sxy)*sx+6.*sxx*sy*ry+3.*rxx*sy**2)*rx+3.*ry*sx*(2.*rxx*sy+2.*rxy*sx+ry*sxx)
cuxxxyy320 = rx**3.*sy**2+6.*rx**2.*ry*sx*sy+3.*rx*ry**2.*sx**2
cuxxxyy030 = 3.*sx**2.*sxyy+(3.*sxx*syy+6.*sxxy*sy+6.*sxy**2)*sx+sxxx*sy**2+6.*sxx*sxy*sy
cuxxxyy130 = ryy*sx**3+(3.*rx*syy+6.*rxy*sy+6.*ry*sxy)*sx**2+(3.*rxx*sy**2+(12.*rx*sxy+6.*ry*sxx)*sy)*sx+3.*rx*sxx*sy**2
cuxxxyy230 = 3.*rx**2.*sx*sy**2+6.*rx*ry*sx**2.*sy+ry**2.*sx**3
cuxxxyy040 = sx**3.*syy+6.*sx**2.*sxy*sy+3.*sx*sxx*sy**2
cuxxxyy140 = 3.*rx*sx**2.*sy**2+2.*ry*sx**3.*sy
cuxxxyy050 = sx**3.*sy**2
                      ! uxxxyy = cuxxxyy100*ur+cuxxxyy200*urr+cuxxxyy300*urrr+cuxxxyy400*urrrr+cuxxxyy500*urrrrr+cuxxxyy010*us+cuxxxyy110*urs+cuxxxyy210*urrs+cuxxxyy310*urrrs+cuxxxyy410*urrrrs+cuxxxyy020*uss+cuxxxyy120*urss+cuxxxyy220*urrss+cuxxxyy320*urrrss+cuxxxyy030*usss+cuxxxyy130*ursss+cuxxxyy230*urrsss+cuxxxyy040*ussss+cuxxxyy140*urssss+cuxxxyy050*usssss
                      ! ------ Coefficients in expansion for uxxyyy ------
cuxxyyy100 = rxxyyy
cuxxyyy200 = 2.*rx*rxyyy+rxx*ryyy+3.*rxxy*ryy+3.*rxxyy*ry+6.*rxy*rxyy
cuxxyyy300 = 3.*rxxy*ry**2+(6.*rx*rxyy+3.*rxx*ryy+6.*rxy**2)*ry+ryyy*rx**2+6.*rxy*ryy*rx
cuxxyyy400 = 3.*rx**2.*ry*ryy+6.*rx*rxy*ry**2+rxx*ry**3
cuxxyyy500 = rx**2.*ry**3
cuxxyyy010 = sxxyyy
cuxxyyy110 = 2.*rx*sxyyy+rxx*syyy+3.*rxxy*syy+3.*rxxyy*sy+6.*rxy*sxyy+6.*rxyy*sxy+2.*rxyyy*sx+3.*ry*sxxyy+3.*ryy*sxxy+ryyy*sxx
cuxxyyy210 = 3.*ry**2.*sxxy+(6.*rx*sxyy+3.*rxx*syy+6.*rxxy*sy+12.*rxy*sxy+6.*rxyy*sx+3.*ryy*sxx)*ry+rx**2.*syyy+(6.*rxy*syy+6.*rxyy*sy+6.*ryy*sxy+2.*ryyy*sx)*rx+6.*rxy*ryy*sx+3.*rxx*ryy*sy+6.*rxy**2.*sy
cuxxyyy310 = ry**3.*sxx+(6.*rx*sxy+3.*rxx*sy+6.*rxy*sx)*ry**2+3.*rx*(rx*syy+4.*rxy*sy+2.*ryy*sx)*ry+3.*ryy*rx**2.*sy
cuxxyyy410 = 3.*rx**2.*ry**2.*sy+2.*rx*ry**3.*sx
cuxxyyy020 = 2.*sx*sxyyy+sxx*syyy+3.*sxxy*syy+3.*sxxyy*sy+6.*sxy*sxyy
cuxxyyy120 = 3.*rxxy*sy**2+(6.*rx*sxyy+3.*rxx*syy+12.*rxy*sxy+6.*rxyy*sx+6.*ry*sxxy+3.*ryy*sxx)*sy+ryyy*sx**2+(2.*rx*syyy+6.*rxy*syy+6.*ry*sxyy+6.*ryy*sxy)*sx+6.*rx*sxy*syy+3.*ry*sxx*syy+6.*ry*sxy**2
cuxxyyy220 = (6.*sx*sxy+3.*sxx*sy)*ry**2+(3.*rxx*sy**2+(12.*rx*sxy+12.*rxy*sx)*sy+6.*sx*syy*rx+3.*ryy*sx**2)*ry+3.*rx*sy*(rx*syy+2.*rxy*sy+2.*ryy*sx)
cuxxyyy320 = 3.*rx**2.*ry*sy**2+6.*rx*ry**2.*sx*sy+ry**3.*sx**2
cuxxyyy030 = 3.*sxxy*sy**2+(6.*sx*sxyy+3.*sxx*syy+6.*sxy**2)*sy+sx**2.*syyy+6.*sxy*syy*sx
cuxxyyy130 = rxx*sy**3+(6.*rx*sxy+6.*rxy*sx+3.*ry*sxx)*sy**2+6.*(rx*syy+2.*ry*sxy+1/2.*ryy*sx)*sx*sy+3.*ry*sx**2.*syy
cuxxyyy230 = rx**2.*sy**3+6.*rx*ry*sx*sy**2+3.*ry**2.*sx**2.*sy
cuxxyyy040 = 3.*sx**2.*sy*syy+6.*sx*sxy*sy**2+sxx*sy**3
cuxxyyy140 = 2.*rx*sx*sy**3+3.*ry*sx**2.*sy**2
cuxxyyy050 = sx**2.*sy**3
                      ! uxxyyy = cuxxyyy100*ur+cuxxyyy200*urr+cuxxyyy300*urrr+cuxxyyy400*urrrr+cuxxyyy500*urrrrr+cuxxyyy010*us+cuxxyyy110*urs+cuxxyyy210*urrs+cuxxyyy310*urrrs+cuxxyyy410*urrrrs+cuxxyyy020*uss+cuxxyyy120*urss+cuxxyyy220*urrss+cuxxyyy320*urrrss+cuxxyyy030*usss+cuxxyyy130*ursss+cuxxyyy230*urrsss+cuxxyyy040*ussss+cuxxyyy140*urssss+cuxxyyy050*usssss
                      ! ------ Coefficients in expansion for uxyyyy ------
cuxyyyy100 = rxyyyy
cuxyyyy200 = rx*ryyyy+4.*rxy*ryyy+6.*rxyy*ryy+4.*rxyyy*ry
cuxyyyy300 = 4.*rx*ry*ryyy+3.*rx*ryy**2+12.*rxy*ry*ryy+6.*rxyy*ry**2
cuxyyyy400 = 6.*ry**2.*(rx*ryy+2/3.*ry*rxy)
cuxyyyy500 = rx*ry**4
cuxyyyy010 = sxyyyy
cuxyyyy110 = rx*syyyy+4.*rxy*syyy+6.*rxyy*syy+4.*rxyyy*sy+4.*ry*sxyyy+6.*ryy*sxyy+4.*ryyy*sxy+ryyyy*sx
cuxyyyy210 = 6.*ry**2.*sxyy+(4.*rx*syyy+12.*rxy*syy+12.*rxyy*sy+12.*ryy*sxy+4.*ryyy*sx)*ry+3.*ryy**2.*sx+(6.*rx*syy+12.*rxy*sy)*ryy+4.*rx*ryyy*sy
cuxyyyy310 = 6.*rx*ry**2.*syy+12.*rx*ry*ryy*sy+12.*rxy*ry**2.*sy+4.*ry**3.*sxy+6.*ry**2.*ryy*sx
cuxyyyy410 = 4.*rx*ry**3.*sy+ry**4.*sx
cuxyyyy020 = sx*syyyy+4.*sxy*syyy+6.*sxyy*syy+4.*sxyyy*sy
cuxyyyy120 = 6.*rxyy*sy**2+(4.*rx*syyy+12.*rxy*syy+12.*ry*sxyy+12.*ryy*sxy+4.*ryyy*sx)*sy+3.*syy**2.*rx+(12.*ry*sxy+6.*ryy*sx)*syy+4.*ry*sx*syyy
cuxyyyy220 = (6.*sx*syy+12.*sxy*sy)*ry**2+12.*sy*(rx*syy+rxy*sy+ryy*sx)*ry+6.*ryy*rx*sy**2
cuxyyyy320 = 6.*ry**2.*sy*(rx*sy+2/3.*ry*sx)
cuxyyyy030 = 4.*sx*sy*syyy+3.*sx*syy**2+12.*sxy*sy*syy+6.*sxyy*sy**2
cuxyyyy130 = 6.*rx*sy**2.*syy+4.*rxy*sy**3+12.*ry*sx*sy*syy+12.*ry*sxy*sy**2+6.*ryy*sx*sy**2
cuxyyyy230 = 4.*ry*sy**2.*(rx*sy+3/2.*ry*sx)
cuxyyyy040 = 6.*sy**2.*(sx*syy+2/3.*sxy*sy)
cuxyyyy140 = rx*sy**4+4.*ry*sx*sy**3
cuxyyyy050 = sx*sy**4
                      ! uxyyyy = cuxyyyy100*ur+cuxyyyy200*urr+cuxyyyy300*urrr+cuxyyyy400*urrrr+cuxyyyy500*urrrrr+cuxyyyy010*us+cuxyyyy110*urs+cuxyyyy210*urrs+cuxyyyy310*urrrs+cuxyyyy410*urrrrs+cuxyyyy020*uss+cuxyyyy120*urss+cuxyyyy220*urrss+cuxyyyy320*urrrss+cuxyyyy030*usss+cuxyyyy130*ursss+cuxyyyy230*urrsss+cuxyyyy040*ussss+cuxyyyy140*urssss+cuxyyyy050*usssss
                      ! ------ Coefficients in expansion for uyyyyy ------
cuyyyyy100 = ryyyyy
cuyyyyy200 = 5.*ry*ryyyy+10.*ryy*ryyy
cuyyyyy300 = 10.*ry**2.*ryyy+15.*ry*ryy**2
cuyyyyy400 = 10.*ry**3.*ryy
cuyyyyy500 = ry**5
cuyyyyy010 = syyyyy
cuyyyyy110 = 5.*ry*syyyy+10.*ryy*syyy+10.*ryyy*syy+5.*ryyyy*sy
cuyyyyy210 = 10.*ry**2.*syyy+30.*ry*ryy*syy+20.*ry*ryyy*sy+15.*ryy**2.*sy
cuyyyyy310 = 10.*ry**3.*syy+30.*ry**2.*ryy*sy
cuyyyyy410 = 5.*sy*ry**4
cuyyyyy020 = 5.*sy*syyyy+10.*syy*syyy
cuyyyyy120 = 20.*ry*sy*syyy+15.*ry*syy**2+30.*ryy*sy*syy+10.*ryyy*sy**2
cuyyyyy220 = 30.*ry**2.*sy*syy+30.*ry*ryy*sy**2
cuyyyyy320 = 10.*ry**3.*sy**2
cuyyyyy030 = 10.*sy**2.*syyy+15.*sy*syy**2
cuyyyyy130 = 30.*ry*sy**2.*syy+10.*ryy*sy**3
cuyyyyy230 = 10.*ry**2.*sy**3
cuyyyyy040 = 10.*sy**3.*syy
cuyyyyy140 = 5.*ry*sy**4
cuyyyyy050 = sy**5
                      ! uyyyyy = cuyyyyy100*ur+cuyyyyy200*urr+cuyyyyy300*urrr+cuyyyyy400*urrrr+cuyyyyy500*urrrrr+cuyyyyy010*us+cuyyyyy110*urs+cuyyyyy210*urrs+cuyyyyy310*urrrs+cuyyyyy410*urrrrs+cuyyyyy020*uss+cuyyyyy120*urss+cuyyyyy220*urrss+cuyyyyy320*urrrss+cuyyyyy030*usss+cuyyyyy130*ursss+cuyyyyy230*urrsss+cuyyyyy040*ussss+cuyyyyy140*urssss+cuyyyyy050*usssss
                      ! uxxxxx = cuxxxxx100*ur+cuxxxxx200*urr+cuxxxxx300*urrr+cuxxxxx400*urrrr+cuxxxxx500*urrrrr+cuxxxxx010*us+cuxxxxx110*urs+cuxxxxx210*urrs+cuxxxxx310*urrrs+cuxxxxx410*urrrrs+cuxxxxx020*uss+cuxxxxx120*urss+cuxxxxx220*urrss+cuxxxxx320*urrrss+cuxxxxx030*usss+cuxxxxx130*ursss+cuxxxxx230*urrsss+cuxxxxx040*ussss+cuxxxxx140*urssss+cuxxxxx050*usssss
                      ! uxxxxy = cuxxxxy100*ur+cuxxxxy200*urr+cuxxxxy300*urrr+cuxxxxy400*urrrr+cuxxxxy500*urrrrr+cuxxxxy010*us+cuxxxxy110*urs+cuxxxxy210*urrs+cuxxxxy310*urrrs+cuxxxxy410*urrrrs+cuxxxxy020*uss+cuxxxxy120*urss+cuxxxxy220*urrss+cuxxxxy320*urrrss+cuxxxxy030*usss+cuxxxxy130*ursss+cuxxxxy230*urrsss+cuxxxxy040*ussss+cuxxxxy140*urssss+cuxxxxy050*usssss
                      ! uxxxyy = cuxxxyy100*ur+cuxxxyy200*urr+cuxxxyy300*urrr+cuxxxyy400*urrrr+cuxxxyy500*urrrrr+cuxxxyy010*us+cuxxxyy110*urs+cuxxxyy210*urrs+cuxxxyy310*urrrs+cuxxxyy410*urrrrs+cuxxxyy020*uss+cuxxxyy120*urss+cuxxxyy220*urrss+cuxxxyy320*urrrss+cuxxxyy030*usss+cuxxxyy130*ursss+cuxxxyy230*urrsss+cuxxxyy040*ussss+cuxxxyy140*urssss+cuxxxyy050*usssss
                      ! uxxyyy = cuxxyyy100*ur+cuxxyyy200*urr+cuxxyyy300*urrr+cuxxyyy400*urrrr+cuxxyyy500*urrrrr+cuxxyyy010*us+cuxxyyy110*urs+cuxxyyy210*urrs+cuxxyyy310*urrrs+cuxxyyy410*urrrrs+cuxxyyy020*uss+cuxxyyy120*urss+cuxxyyy220*urrss+cuxxyyy320*urrrss+cuxxyyy030*usss+cuxxyyy130*ursss+cuxxyyy230*urrsss+cuxxyyy040*ussss+cuxxyyy140*urssss+cuxxyyy050*usssss
                      ! uxyyyy = cuxyyyy100*ur+cuxyyyy200*urr+cuxyyyy300*urrr+cuxyyyy400*urrrr+cuxyyyy500*urrrrr+cuxyyyy010*us+cuxyyyy110*urs+cuxyyyy210*urrs+cuxyyyy310*urrrs+cuxyyyy410*urrrrs+cuxyyyy020*uss+cuxyyyy120*urss+cuxyyyy220*urrss+cuxyyyy320*urrrss+cuxyyyy030*usss+cuxyyyy130*ursss+cuxyyyy230*urrsss+cuxyyyy040*ussss+cuxyyyy140*urssss+cuxyyyy050*usssss
                      ! uyyyyy = cuyyyyy100*ur+cuyyyyy200*urr+cuyyyyy300*urrr+cuyyyyy400*urrrr+cuyyyyy500*urrrrr+cuyyyyy010*us+cuyyyyy110*urs+cuyyyyy210*urrs+cuyyyyy310*urrrs+cuyyyyy410*urrrrs+cuyyyyy020*uss+cuyyyyy120*urss+cuyyyyy220*urrss+cuyyyyy320*urrrss+cuyyyyy030*usss+cuyyyyy130*ursss+cuyyyyy230*urrsss+cuyyyyy040*ussss+cuyyyyy140*urssss+cuyyyyy050*usssss
                      ! ----- START SIXTH DERIVATIVES -----
                                            do m1=0,nd-1
                                            do m2=0,nd-1
                        ! ---- 5th parameteric derivatives of the metrics ----
                                                rxrrrrra(m1,m2) = rx500i(m1,m2) !- (dx(0)**2/3.)*rx700i(m1,m2)
                                                rxsssssa(m1,m2) = rx050i(m1,m2) !- (dx(1)**2/3.)*rx070i(m1,m2) 
                                                rxrrrrsa(m1,m2) = rx410i(m1,m2) !- (dx(0)**2/6.)*rx610i(m1,m2) -  (dx(1)**2/6.)*rx430i(m1,m2)
                                                rxrssssa(m1,m2) = rx140i(m1,m2) !- (dx(1)**2/6.)*rx160i(m1,m2) -  (dx(0)**2/6.)*rx340i(m1,m2)    
                                                rxrrrssa(m1,m2) = rx320i(m1,m2) !- (dx(0)**2/4.)*rx520i(m1,m2) - (dx(1)**2/12.)*rx340i(m1,m2)
                                                rxrrsssa(m1,m2) = rx230i(m1,m2) !- (dx(1)**2/4.)*rx250i(m1,m2) - (dx(0)**2/12.)*rx430i(m1,m2)
                        ! ---- fixth spatial derivatives of the metrics ----
                                                rxxxxxxa(m1,m2) = cuxxxxx100*rxr(m1,m2,0)+cuxxxxx200*rxrra(m1,m2)+cuxxxxx300*rxrrra(m1,m2)+cuxxxxx400*rxrrrra(m1,m2)+cuxxxxx500*rxrrrrra(m1,m2)+cuxxxxx010*rxr(m1,m2,1)+cuxxxxx110*rxrsa(m1,m2)+cuxxxxx210*rxrrsa(m1,m2)+cuxxxxx310*rxrrrsa(m1,m2)+cuxxxxx410*rxrrrrsa(m1,m2)+cuxxxxx020*rxssa(m1,m2)+cuxxxxx120*rxrssa(m1,m2)+cuxxxxx220*rxrrssa(m1,m2)+cuxxxxx320*rxrrrssa(m1,m2)+cuxxxxx030*rxsssa(m1,m2)+cuxxxxx130*rxrsssa(m1,m2)+cuxxxxx230*rxrrsssa(m1,m2)+cuxxxxx040*rxssssa(m1,m2)+cuxxxxx140*rxrssssa(m1,m2)+cuxxxxx050*rxsssssa(m1,m2)
                                                rxxxxxya(m1,m2) = cuxxxxy100*rxr(m1,m2,0)+cuxxxxy200*rxrra(m1,m2)+cuxxxxy300*rxrrra(m1,m2)+cuxxxxy400*rxrrrra(m1,m2)+cuxxxxy500*rxrrrrra(m1,m2)+cuxxxxy010*rxr(m1,m2,1)+cuxxxxy110*rxrsa(m1,m2)+cuxxxxy210*rxrrsa(m1,m2)+cuxxxxy310*rxrrrsa(m1,m2)+cuxxxxy410*rxrrrrsa(m1,m2)+cuxxxxy020*rxssa(m1,m2)+cuxxxxy120*rxrssa(m1,m2)+cuxxxxy220*rxrrssa(m1,m2)+cuxxxxy320*rxrrrssa(m1,m2)+cuxxxxy030*rxsssa(m1,m2)+cuxxxxy130*rxrsssa(m1,m2)+cuxxxxy230*rxrrsssa(m1,m2)+cuxxxxy040*rxssssa(m1,m2)+cuxxxxy140*rxrssssa(m1,m2)+cuxxxxy050*rxsssssa(m1,m2)
                                                rxxxxyya(m1,m2) = cuxxxyy100*rxr(m1,m2,0)+cuxxxyy200*rxrra(m1,m2)+cuxxxyy300*rxrrra(m1,m2)+cuxxxyy400*rxrrrra(m1,m2)+cuxxxyy500*rxrrrrra(m1,m2)+cuxxxyy010*rxr(m1,m2,1)+cuxxxyy110*rxrsa(m1,m2)+cuxxxyy210*rxrrsa(m1,m2)+cuxxxyy310*rxrrrsa(m1,m2)+cuxxxyy410*rxrrrrsa(m1,m2)+cuxxxyy020*rxssa(m1,m2)+cuxxxyy120*rxrssa(m1,m2)+cuxxxyy220*rxrrssa(m1,m2)+cuxxxyy320*rxrrrssa(m1,m2)+cuxxxyy030*rxsssa(m1,m2)+cuxxxyy130*rxrsssa(m1,m2)+cuxxxyy230*rxrrsssa(m1,m2)+cuxxxyy040*rxssssa(m1,m2)+cuxxxyy140*rxrssssa(m1,m2)+cuxxxyy050*rxsssssa(m1,m2)
                                                rxxxyyya(m1,m2) = cuxxyyy100*rxr(m1,m2,0)+cuxxyyy200*rxrra(m1,m2)+cuxxyyy300*rxrrra(m1,m2)+cuxxyyy400*rxrrrra(m1,m2)+cuxxyyy500*rxrrrrra(m1,m2)+cuxxyyy010*rxr(m1,m2,1)+cuxxyyy110*rxrsa(m1,m2)+cuxxyyy210*rxrrsa(m1,m2)+cuxxyyy310*rxrrrsa(m1,m2)+cuxxyyy410*rxrrrrsa(m1,m2)+cuxxyyy020*rxssa(m1,m2)+cuxxyyy120*rxrssa(m1,m2)+cuxxyyy220*rxrrssa(m1,m2)+cuxxyyy320*rxrrrssa(m1,m2)+cuxxyyy030*rxsssa(m1,m2)+cuxxyyy130*rxrsssa(m1,m2)+cuxxyyy230*rxrrsssa(m1,m2)+cuxxyyy040*rxssssa(m1,m2)+cuxxyyy140*rxrssssa(m1,m2)+cuxxyyy050*rxsssssa(m1,m2)
                                                rxxyyyya(m1,m2) = cuxyyyy100*rxr(m1,m2,0)+cuxyyyy200*rxrra(m1,m2)+cuxyyyy300*rxrrra(m1,m2)+cuxyyyy400*rxrrrra(m1,m2)+cuxyyyy500*rxrrrrra(m1,m2)+cuxyyyy010*rxr(m1,m2,1)+cuxyyyy110*rxrsa(m1,m2)+cuxyyyy210*rxrrsa(m1,m2)+cuxyyyy310*rxrrrsa(m1,m2)+cuxyyyy410*rxrrrrsa(m1,m2)+cuxyyyy020*rxssa(m1,m2)+cuxyyyy120*rxrssa(m1,m2)+cuxyyyy220*rxrrssa(m1,m2)+cuxyyyy320*rxrrrssa(m1,m2)+cuxyyyy030*rxsssa(m1,m2)+cuxyyyy130*rxrsssa(m1,m2)+cuxyyyy230*rxrrsssa(m1,m2)+cuxyyyy040*rxssssa(m1,m2)+cuxyyyy140*rxrssssa(m1,m2)+cuxyyyy050*rxsssssa(m1,m2)
                                                rxyyyyya(m1,m2) = cuyyyyy100*rxr(m1,m2,0)+cuyyyyy200*rxrra(m1,m2)+cuyyyyy300*rxrrra(m1,m2)+cuyyyyy400*rxrrrra(m1,m2)+cuyyyyy500*rxrrrrra(m1,m2)+cuyyyyy010*rxr(m1,m2,1)+cuyyyyy110*rxrsa(m1,m2)+cuyyyyy210*rxrrsa(m1,m2)+cuyyyyy310*rxrrrsa(m1,m2)+cuyyyyy410*rxrrrrsa(m1,m2)+cuyyyyy020*rxssa(m1,m2)+cuyyyyy120*rxrssa(m1,m2)+cuyyyyy220*rxrrssa(m1,m2)+cuyyyyy320*rxrrrssa(m1,m2)+cuyyyyy030*rxsssa(m1,m2)+cuyyyyy130*rxrsssa(m1,m2)+cuyyyyy230*rxrrsssa(m1,m2)+cuyyyyy040*rxssssa(m1,m2)+cuyyyyy140*rxrssssa(m1,m2)+cuyyyyy050*rxsssssa(m1,m2)
                                            end do
                                            end do
                                            rxxxxxx=rxxxxxxa(0,0); ryxxxxx=rxxxxxxa(0,1); sxxxxxx=rxxxxxxa(1,0); syxxxxx=rxxxxxxa(1,1); 
                                            rxxxxxy=rxxxxxya(0,0); ryxxxxy=rxxxxxya(0,1); sxxxxxy=rxxxxxya(1,0); syxxxxy=rxxxxxya(1,1); 
                                            rxxxxyy=rxxxxyya(0,0); ryxxxyy=rxxxxyya(0,1); sxxxxyy=rxxxxyya(1,0); syxxxyy=rxxxxyya(1,1); 
                                            rxxxyyy=rxxxyyya(0,0); ryxxyyy=rxxxyyya(0,1); sxxxyyy=rxxxyyya(1,0); syxxyyy=rxxxyyya(1,1); 
                                            rxxyyyy=rxxyyyya(0,0); ryxyyyy=rxxyyyya(0,1); sxxyyyy=rxxyyyya(1,0); syxyyyy=rxxyyyya(1,1); 
                                            rxyyyyy=rxyyyyya(0,0); ryyyyyy=rxyyyyya(0,1); sxyyyyy=rxyyyyya(1,0); syyyyyy=rxyyyyya(1,1); 
                       ! if( i1.eq.5 .and. i2.eq.6 )then  
                       !    write(*,'(" (i1,i2)=(",2i3,") rxxxxxx,ryxxxxx,sxxxxxx,syxxxxx=",4(1pe12.4,1x))') i1,i2,rxxxxxx,ryxxxxx,sxxxxxx,syxxxxx
                       !    write(*,'(" (i1,i2)=(",2i3,") rxxxxxy,ryxxxxy,sxxxxxy,syxxxxy=",4(1pe12.4,1x))') i1,i2,rxxxxxy,ryxxxxy,sxxxxxy,syxxxxy
                       !    write(*,'(" (i1,i2)=(",2i3,") rxxxxyy,ryxxxyy,sxxxxyy,syxxxyy=",4(1pe12.4,1x))') i1,i2,rxxxxyy,ryxxxyy,sxxxxyy,syxxxyy
                       !    write(*,'(" (i1,i2)=(",2i3,") rxxxyyy,ryxxyyy,sxxxyyy,syxxyyy=",4(1pe12.4,1x))') i1,i2,rxxxyyy,ryxxyyy,sxxxyyy,syxxyyy
                       !    write(*,'(" (i1,i2)=(",2i3,") rxxyyyy,ryxyyyy,sxxyyyy,syxyyyy=",4(1pe12.4,1x))') i1,i2,rxxyyyy,ryxyyyy,sxxyyyy,syxyyyy
                       !    write(*,'(" (i1,i2)=(",2i3,") rxyyyyy,ryyyyyy,sxyyyyy,syyyyyy=",4(1pe12.4,1x))') i1,i2,rxyyyyy,ryyyyyy,sxyyyyy,syyyyyy
                       !  end if      
                      ! ---- SIXTH parametric derivs----
                      !  1/4 , 13/240 
                                            urrrrrr = d600i ! - (dx(0)**2/4.)*d800i
                                            ussssss = d060i ! - (dx(1)**2/4.)*d080i  
                      ! uxxxxxy = D0x(D+xD-x)^2 *( 1 - (1/3)*dx^2 D+xD-x + ... ) X 
                      !                      D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                                            urrrrrs = d510i ! - (dx(0)**2/3.)*d710i - (dx(1)**2/6.)*d530i
                                            ursssss = d150i ! - (dx(1)**2/3.)*d170i - (dx(0)**2/6.)*d350i
                      ! uxxxxyy = (D+xD-x)^2*( 1 - dx^2*(1/6)*D+xD-x + dx^4*(7/240)*(D+xD-x)^2 + ... )
                      !               D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                            urrrrss = d420i ! - (dx(0)**2/6.)*d620i - (dx(1)**2/12.)*d440i
                                            urrssss = d240i ! - (dx(1)**2/6.)*d260i - (dx(0)**2/12.)*d440i
                      ! uxxxyyy = D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ... ) X 
                      !           D0y D+yD-y*( 1 - (1/4)*dy^2 * D+yD-y + (7/120)*dy^4 (D+yD-y)^2 + ...  ) u 
                                            urrrsss = d330i ! - (dx(0)**2/4.)*d530i - (dx(1)**2/4.)*d350i
                      ! ----- SIXTH SPATIAL DERIVATIVES -----
                      ! ------ Coefficients in expansion for uxxxxxx ------
cuxxxxxx100 = rxxxxxx
cuxxxxxx200 = 6.*rx*rxxxxx+15.*rxx*rxxxx+10.*rxxx**2
cuxxxxxx300 = 15.*rx**2.*rxxxx+60.*rx*rxx*rxxx+15.*rxx**3
cuxxxxxx400 = 20.*rx**3.*rxxx+45.*rx**2.*rxx**2
cuxxxxxx500 = 15.*rx**4.*rxx
cuxxxxxx600 = rx**6
cuxxxxxx010 = sxxxxxx
cuxxxxxx110 = 6.*rx*sxxxxx+15.*rxx*sxxxx+20.*rxxx*sxxx+15.*rxxxx*sxx+6.*rxxxxx*sx
cuxxxxxx210 = 15.*rx**2.*sxxxx+(60.*rxx*sxxx+60.*rxxx*sxx+30.*rxxxx*sx)*rx+60.*rxxx*rxx*sx+45.*sxx*rxx**2
cuxxxxxx310 = 20.*rx**3.*sxxx+(90.*rxx*sxx+60.*rxxx*sx)*rx**2+90.*rxx**2.*sx*rx
cuxxxxxx410 = 15.*rx**4.*sxx+60.*rx**3.*rxx*sx
cuxxxxxx510 = 6.*rx**5.*sx
cuxxxxxx020 = 6.*sx*sxxxxx+15.*sxx*sxxxx+10.*sxxx**2
cuxxxxxx120 = 15.*rxxxx*sx**2+(30.*rx*sxxxx+60.*rxx*sxxx+60.*rxxx*sxx)*sx+60.*rx*sxx*sxxx+45.*sxx**2.*rxx
cuxxxxxx220 = (60.*rx*rxxx+45.*rxx**2)*sx**2+(60.*rx**2.*sxxx+180.*rx*rxx*sxx)*sx+45.*rx**2.*sxx**2
cuxxxxxx320 = 60.*sx*rx**2.*(rx*sxx+3./2.*rxx*sx)
cuxxxxxx420 = 15.*rx**4.*sx**2
cuxxxxxx030 = 15.*sx**2.*sxxxx+60.*sx*sxx*sxxx+15.*sxx**3
cuxxxxxx130 = 20.*rxxx*sx**3+(60.*rx*sxxx+90.*rxx*sxx)*sx**2+90.*rx*sx*sxx**2
cuxxxxxx230 = 90.*sx**2.*rx*(rx*sxx+2./3.*rxx*sx)
cuxxxxxx330 = 20.*rx**3.*sx**3
cuxxxxxx040 = 20.*sx**3.*sxxx+45.*sx**2.*sxx**2
cuxxxxxx140 = 60.*rx*sx**3.*sxx+15.*rxx*sx**4
cuxxxxxx240 = 15.*rx**2.*sx**4
cuxxxxxx050 = 15.*sxx*sx**4
cuxxxxxx150 = 6.*rx*sx**5
cuxxxxxx060 = sx**6
                      ! uxxxxxx = cuxxxxxx100*ur+cuxxxxxx200*urr+cuxxxxxx300*urrr+cuxxxxxx400*urrrr+cuxxxxxx500*urrrrr+cuxxxxxx600*urrrrrr+cuxxxxxx010*us+cuxxxxxx110*urs+cuxxxxxx210*urrs+cuxxxxxx310*urrrs+cuxxxxxx410*urrrrs+cuxxxxxx510*urrrrrs+cuxxxxxx020*uss+cuxxxxxx120*urss+cuxxxxxx220*urrss+cuxxxxxx320*urrrss+cuxxxxxx420*urrrrss+cuxxxxxx030*usss+cuxxxxxx130*ursss+cuxxxxxx230*urrsss+cuxxxxxx330*urrrsss+cuxxxxxx040*ussss+cuxxxxxx140*urssss+cuxxxxxx240*urrssss+cuxxxxxx050*usssss+cuxxxxxx150*ursssss+cuxxxxxx060*ussssss
                      ! ------ Coefficients in expansion for uxxxxxy ------
cuxxxxxy100 = rxxxxxy
cuxxxxxy200 = 5.*rx*rxxxxy+10.*rxx*rxxxy+10.*rxxx*rxxy+5.*rxxxx*rxy+rxxxxx*ry
cuxxxxxy300 = 10.*rxxxy*rx**2+(30.*rxx*rxxy+20.*rxxx*rxy+5.*rxxxx*ry)*rx+15.*rxy*rxx**2+10.*rxxx*ry*rxx
cuxxxxxy400 = 10.*rx**3.*rxxy+30.*rx**2.*rxx*rxy+10.*rx**2.*rxxx*ry+15.*rx*rxx**2.*ry
cuxxxxxy500 = 5.*rx**4.*rxy+10.*rx**3.*rxx*ry
cuxxxxxy600 = ry*rx**5
cuxxxxxy010 = sxxxxxy
cuxxxxxy110 = 5.*rx*sxxxxy+10.*rxx*sxxxy+10.*rxxx*sxxy+5.*rxxxx*sxy+rxxxxx*sy+5.*rxxxxy*sx+10.*rxxxy*sxx+10.*rxxy*sxxx+5.*rxy*sxxxx+ry*sxxxxx
cuxxxxxy210 = 10.*rx**2.*sxxxy+(30.*rxx*sxxy+20.*rxxx*sxy+5.*rxxxx*sy+20.*rxxxy*sx+30.*rxxy*sxx+20.*rxy*sxxx+5.*ry*sxxxx)*rx+15.*sxy*rxx**2+(10.*rxxx*sy+30.*rxxy*sx+30.*rxy*sxx+10.*ry*sxxx)*rxx+(20.*rxy*sx+10.*ry*sxx)*rxxx+5.*rxxxx*ry*sx
cuxxxxxy310 = 10.*rx**3.*sxxy+(30.*rxx*sxy+10.*rxxx*sy+30.*rxxy*sx+30.*rxy*sxx+10.*ry*sxxx)*rx**2+(15.*rxx**2.*sy+(60.*rxy*sx+30.*ry*sxx)*rxx+20.*ry*sx*rxxx)*rx+15.*rxx**2.*ry*sx
cuxxxxxy410 = 5.*rx**2.*(rx**2.*sxy+(2.*rxx*sy+4.*rxy*sx+2.*ry*sxx)*rx+6.*rxx*ry*sx)
cuxxxxxy510 = rx**5.*sy+5.*rx**4.*ry*sx
cuxxxxxy020 = 5.*sx*sxxxxy+10.*sxx*sxxxy+10.*sxxx*sxxy+5.*sxxxx*sxy+sxxxxx*sy
cuxxxxxy120 = 10.*rxxxy*sx**2+(20.*rx*sxxxy+30.*rxx*sxxy+20.*rxxx*sxy+5.*rxxxx*sy+30.*rxxy*sxx+20.*rxy*sxxx+5.*ry*sxxxx)*sx+15.*rxy*sxx**2+(30.*rx*sxxy+30.*rxx*sxy+10.*rxxx*sy+10.*ry*sxxx)*sxx+(20.*rx*sxy+10.*rxx*sy)*sxxx+5.*rx*sxxxx*sy
cuxxxxxy220 = (30.*sx*sxxy+30.*sxx*sxy+10.*sxxx*sy)*rx**2+(30.*rxxy*sx**2+(60.*rxx*sxy+20.*rxxx*sy+60.*rxy*sxx+20.*ry*sxxx)*sx+30.*rxx*sxx*sy+15.*sxx**2.*ry)*rx+10.*((3.*rxx*rxy+rxxx*ry)*sx+3.*rxx*sxx*ry+3/2.*rxx**2.*sy)*sx
cuxxxxxy320 = 20.*rx**3.*sx*sxy+10.*rx**3.*sxx*sy+30.*rx**2.*rxx*sx*sy+30.*rx**2.*rxy*sx**2+30.*rx**2.*ry*sx*sxx+30.*rx*rxx*ry*sx**2
cuxxxxxy420 = 5.*rx**4.*sx*sy+10.*rx**3.*ry*sx**2
cuxxxxxy030 = 10.*sx**2.*sxxxy+(30.*sxx*sxxy+20.*sxxx*sxy+5.*sxxxx*sy)*sx+15.*sxy*sxx**2+10.*sxxx*sy*sxx
cuxxxxxy130 = 10.*rxxy*sx**3+(30.*rx*sxxy+30.*rxx*sxy+10.*rxxx*sy+30.*rxy*sxx+10.*ry*sxxx)*sx**2+(15.*sxx**2.*ry+(60.*rx*sxy+30.*rxx*sy)*sxx+20.*rx*sxxx*sy)*sx+15.*rx*sxx**2.*sy
cuxxxxxy230 = 30.*rx**2.*sx**2.*sxy+30.*rx**2.*sx*sxx*sy+30.*rx*rxx*sx**2.*sy+20.*rx*rxy*sx**3+30.*rx*ry*sx**2.*sxx+10.*rxx*ry*sx**3
cuxxxxxy330 = 10.*rx**3.*sx**2.*sy+10.*rx**2.*ry*sx**3
cuxxxxxy040 = 10.*sx**3.*sxxy+30.*sx**2.*sxx*sxy+10.*sx**2.*sxxx*sy+15.*sx*sxx**2.*sy
cuxxxxxy140 = 5.*rxy*sx**4+(20.*rx*sxy+10.*rxx*sy+10.*ry*sxx)*sx**3+30.*rx*sx**2.*sxx*sy
cuxxxxxy240 = 10.*rx**2.*sx**3.*sy+5.*rx*ry*sx**4
cuxxxxxy050 = 5.*sx**4.*sxy+10.*sx**3.*sxx*sy
cuxxxxxy150 = 5.*rx*sx**4.*sy+ry*sx**5
cuxxxxxy060 = sy*sx**5
                      ! uxxxxxy = cuxxxxxy100*ur+cuxxxxxy200*urr+cuxxxxxy300*urrr+cuxxxxxy400*urrrr+cuxxxxxy500*urrrrr+cuxxxxxy600*urrrrrr+cuxxxxxy010*us+cuxxxxxy110*urs+cuxxxxxy210*urrs+cuxxxxxy310*urrrs+cuxxxxxy410*urrrrs+cuxxxxxy510*urrrrrs+cuxxxxxy020*uss+cuxxxxxy120*urss+cuxxxxxy220*urrss+cuxxxxxy320*urrrss+cuxxxxxy420*urrrrss+cuxxxxxy030*usss+cuxxxxxy130*ursss+cuxxxxxy230*urrsss+cuxxxxxy330*urrrsss+cuxxxxxy040*ussss+cuxxxxxy140*urssss+cuxxxxxy240*urrssss+cuxxxxxy050*usssss+cuxxxxxy150*ursssss+cuxxxxxy060*ussssss
                      ! ------ Coefficients in expansion for uxxxxyy ------
cuxxxxyy100 = rxxxxyy
cuxxxxyy200 = 4.*rx*rxxxyy+6.*rxx*rxxyy+4.*rxxx*rxyy+rxxxx*ryy+2.*rxxxxy*ry+8.*rxxxy*rxy+6.*rxxy**2
cuxxxxyy300 = 6.*rxxyy*rx**2+(12.*rxx*rxyy+4.*rxxx*ryy+8.*rxxxy*ry+24.*rxxy*rxy)*rx+3.*ryy*rxx**2+(12.*rxxy*ry+12.*rxy**2)*rxx+rxxxx*ry**2+8.*rxxx*rxy*ry
cuxxxxyy400 = (24.*rxx*rxy*ry+4.*rxxx*ry**2)*rx+(6.*rxx*ryy+12.*rxxy*ry+12.*rxy**2)*rx**2+4.*rxyy*rx**3+3.*rxx**2.*ry**2
cuxxxxyy500 = rx**4.*ryy+8.*rx**3.*rxy*ry+6.*rx**2.*rxx*ry**2
cuxxxxyy600 = ry**2.*rx**4
cuxxxxyy010 = sxxxxyy
cuxxxxyy110 = 4.*rx*sxxxyy+6.*rxx*sxxyy+4.*rxxx*sxyy+rxxxx*syy+2.*rxxxxy*sy+8.*rxxxy*sxy+4.*rxxxyy*sx+12.*rxxy*sxxy+6.*rxxyy*sxx+8.*rxy*sxxxy+4.*rxyy*sxxx+2.*ry*sxxxxy+ryy*sxxxx
cuxxxxyy210 = 6.*rx**2.*sxxyy+(12.*rxx*sxyy+4.*rxxx*syy+8.*rxxxy*sy+24.*rxxy*sxy+12.*rxxyy*sx+24.*rxy*sxxy+12.*rxyy*sxx+8.*ry*sxxxy+4.*ryy*sxxx)*rx+ry**2.*sxxxx+(12.*rxx*sxxy+8.*rxxx*sxy+2.*rxxxx*sy+8.*rxxxy*sx+12.*rxxy*sxx+8.*rxy*sxxx)*ry+3.*syy*rxx**2+(12.*rxxy*sy+24.*rxy*sxy+12.*rxyy*sx+6.*ryy*sxx)*rxx+12.*rxy**2.*sxx+(8.*rxxx*sy+24.*rxxy*sx)*rxy+4.*rxxx*ryy*sx
cuxxxxyy310 = 4.*rx**3.*sxyy+(6.*rxx*syy+12.*rxxy*sy+24.*rxy*sxy+12.*rxyy*sx+12.*ry*sxxy+6.*ryy*sxx)*rx**2+(4.*ry**2.*sxxx+(24.*rxx*sxy+8.*rxxx*sy+24.*rxxy*sx+24.*rxy*sxx)*ry+(24.*rxy*sy+12.*ryy*sx)*rxx+24.*sx*rxy**2)*rx+(6.*rxx*sxx+4.*rxxx*sx)*ry**2+(6.*rxx**2.*sy+24.*rxx*rxy*sx)*ry
cuxxxxyy410 = rx*(rx**3.*syy+(8.*rxy*sy+8.*ry*sxy+4.*ryy*sx)*rx**2+6.*ry*(2.*rxx*sy+4.*rxy*sx+ry*sxx)*rx+12.*rxx*ry**2.*sx)
cuxxxxyy510 = 2.*rx**4.*ry*sy+4.*rx**3.*ry**2.*sx
cuxxxxyy020 = 4.*sx*sxxxyy+6.*sxx*sxxyy+4.*sxxx*sxyy+sxxxx*syy+2.*sxxxxy*sy+8.*sxxxy*sxy+6.*sxxy**2
cuxxxxyy120 = 6.*rxxyy*sx**2+(12.*rx*sxxyy+12.*rxx*sxyy+4.*rxxx*syy+8.*rxxxy*sy+24.*rxxy*sxy+24.*rxy*sxxy+12.*rxyy*sxx+8.*ry*sxxxy+4.*ryy*sxxx)*sx+rxxxx*sy**2+(8.*rx*sxxxy+12.*rxx*sxxy+8.*rxxx*sxy+12.*rxxy*sxx+8.*rxy*sxxx+2.*ry*sxxxx)*sy+3.*ryy*sxx**2+(12.*rx*sxyy+6.*rxx*syy+24.*rxy*sxy+12.*ry*sxxy)*sxx+12.*rxx*sxy**2+(24.*rx*sxxy+8.*ry*sxxx)*sxy+4.*rx*sxxx*syy
cuxxxxyy220 = (12.*sx*sxyy+6.*sxx*syy+12.*sxxy*sy+12.*sxy**2)*rx**2+(12.*rxyy*sx**2+(12.*rxx*syy+24.*rxxy*sy+48.*rxy*sxy+24.*ry*sxxy+12.*ryy*sxx)*sx+(24.*sxx*sxy+8.*sxxx*sy)*ry+4.*sy*(6.*rxx*sxy+rxxx*sy+6.*rxy*sxx))*rx+(6.*rxx*ryy+12.*rxxy*ry+12.*rxy**2)*sx**2+(4.*ry**2.*sxxx+(24.*rxx*sxy+8.*rxxx*sy+24.*rxy*sxx)*ry+24.*rxy*rxx*sy)*sx+3.*sxx**2.*ry**2+12.*rxx*sxx*sy*ry+3.*rxx**2.*sy**2
cuxxxxyy320 = (4.*sx*syy+8.*sxy*sy)*rx**3+(6.*ryy*sx**2+(24.*rxy*sy+24.*ry*sxy)*sx+12.*sxx*sy*ry+6.*rxx*sy**2)*rx**2+12.*ry*sx*(2.*rxx*sy+2.*rxy*sx+ry*sxx)*rx+6.*rxx*ry**2.*sx**2
cuxxxxyy420 = rx**4.*sy**2+8.*rx**3.*ry*sx*sy+6.*rx**2.*ry**2.*sx**2
cuxxxxyy030 = 6.*sx**2.*sxxyy+(12.*sxx*sxyy+4.*sxxx*syy+8.*sxxxy*sy+24.*sxxy*sxy)*sx+3.*syy*sxx**2+(12.*sxxy*sy+12.*sxy**2)*sxx+sxxxx*sy**2+8.*sxxx*sxy*sy
cuxxxxyy130 = 4.*rxyy*sx**3+(12.*rx*sxyy+6.*rxx*syy+12.*rxxy*sy+24.*rxy*sxy+12.*ry*sxxy+6.*ryy*sxx)*sx**2+(4.*rxxx*sy**2+(24.*rx*sxxy+24.*rxx*sxy+24.*rxy*sxx+8.*ry*sxxx)*sy+(12.*rx*syy+24.*ry*sxy)*sxx+24.*rx*sxy**2)*sx+(4.*rx*sxxx+6.*rxx*sxx)*sy**2+(24.*rx*sxx*sxy+6.*ry*sxx**2)*sy
cuxxxxyy230 = (4.*rx*ryy+8.*rxy*ry)*sx**3+(6.*rx**2.*syy+(24.*rxy*sy+24.*ry*sxy)*rx+6.*ry**2.*sxx+12.*ry*sy*rxx)*sx**2+24.*sy*(rx*sxy+ry*sxx+1/2.*rxx*sy)*rx*sx+6.*rx**2.*sxx*sy**2
cuxxxxyy330 = 4.*rx**3.*sx*sy**2+12.*rx**2.*ry*sx**2.*sy+4.*rx*ry**2.*sx**3
cuxxxxyy040 = 4.*sx**3.*sxyy+(6.*sxx*syy+12.*sxxy*sy+12.*sxy**2)*sx**2+3.*sxx**2.*sy**2+(24.*sxx*sxy*sy+4.*sxxx*sy**2)*sx
cuxxxxyy140 = 4.*(1/4.*ryy*sx**3+(rx*syy+2.*rxy*sy+2.*ry*sxy)*sx**2+6.*(rx*sxy+1/2.*ry*sxx+1/4.*rxx*sy)*sy*sx+3.*rx*sy**2.*sxx)*sx
cuxxxxyy240 = 6.*rx**2.*sx**2.*sy**2+8.*rx*ry*sx**3.*sy+ry**2.*sx**4
cuxxxxyy050 = sx**4.*syy+8.*sx**3.*sxy*sy+6.*sx**2.*sxx*sy**2
cuxxxxyy150 = 4.*rx*sx**3.*sy**2+2.*ry*sx**4.*sy
cuxxxxyy060 = sy**2.*sx**4
                      ! uxxxxyy = cuxxxxyy100*ur+cuxxxxyy200*urr+cuxxxxyy300*urrr+cuxxxxyy400*urrrr+cuxxxxyy500*urrrrr+cuxxxxyy600*urrrrrr+cuxxxxyy010*us+cuxxxxyy110*urs+cuxxxxyy210*urrs+cuxxxxyy310*urrrs+cuxxxxyy410*urrrrs+cuxxxxyy510*urrrrrs+cuxxxxyy020*uss+cuxxxxyy120*urss+cuxxxxyy220*urrss+cuxxxxyy320*urrrss+cuxxxxyy420*urrrrss+cuxxxxyy030*usss+cuxxxxyy130*ursss+cuxxxxyy230*urrsss+cuxxxxyy330*urrrsss+cuxxxxyy040*ussss+cuxxxxyy140*urssss+cuxxxxyy240*urrssss+cuxxxxyy050*usssss+cuxxxxyy150*ursssss+cuxxxxyy060*ussssss
                      ! ------ Coefficients in expansion for uxxxyyy ------
cuxxxyyy100 = rxxxyyy
cuxxxyyy200 = 3.*rx*rxxyyy+3.*rxx*rxyyy+rxxx*ryyy+3.*rxxxy*ryy+3.*rxxxyy*ry+9.*rxxy*rxyy+9.*rxxyy*rxy
cuxxxyyy300 = (9.*rxx*rxyy+3.*rxxx*ryy+18.*rxxy*rxy)*ry+(3.*rxx*ryyy+9.*rxxy*ryy+9.*rxxyy*ry+18.*rxy*rxyy)*rx+3.*rxxxy*ry**2+3.*rxyyy*rx**2+9.*rxx*ryy*rxy+6.*rxy**3
cuxxxyyy400 = ryyy*rx**3+(9.*rxy*ryy+9.*rxyy*ry)*rx**2+9.*ry*(rxx*ryy+rxxy*ry+2.*rxy**2)*rx+rxxx*ry**3+9.*rxy*rxx*ry**2
cuxxxyyy500 = 3.*rx**3.*ry*ryy+9.*rx**2.*rxy*ry**2+3.*rx*rxx*ry**3
cuxxxyyy600 = rx**3.*ry**3
cuxxxyyy010 = sxxxyyy
cuxxxyyy110 = 3.*rx*sxxyyy+3.*rxx*sxyyy+rxxx*syyy+3.*rxxxy*syy+3.*rxxxyy*sy+9.*rxxy*sxyy+9.*rxxyy*sxy+3.*rxxyyy*sx+9.*rxy*sxxyy+9.*rxyy*sxxy+3.*rxyyy*sxx+3.*ry*sxxxyy+3.*ryy*sxxxy+ryyy*sxxx
cuxxxyyy210 = 3.*rx**2.*sxyyy+(3.*rxx*syyy+9.*rxxy*syy+9.*rxxyy*sy+18.*rxy*sxyy+18.*rxyy*sxy+6.*rxyyy*sx+9.*ry*sxxyy+9.*ryy*sxxy+3.*ryyy*sxx)*rx+3.*ry**2.*sxxxy+(9.*rxx*sxyy+3.*rxxx*syy+6.*rxxxy*sy+18.*rxxy*sxy+9.*rxxyy*sx+18.*rxy*sxxy+9.*rxyy*sxx+3.*ryy*sxxx)*ry+18.*rxy**2.*sxy+(9.*rxx*syy+18.*rxxy*sy+18.*rxyy*sx+9.*ryy*sxx)*rxy+(9.*rxyy*sy+9.*ryy*sxy+3.*ryyy*sx)*rxx+3.*ryy*(rxxx*sy+3.*rxxy*sx)
cuxxxyyy310 = rx**3.*syyy+(9.*rxy*syy+9.*rxyy*sy+9.*ry*sxyy+9.*ryy*sxy+3.*ryyy*sx)*rx**2+(9.*ry**2.*sxxy+(9.*rxx*syy+18.*rxxy*sy+36.*rxy*sxy+18.*rxyy*sx+9.*ryy*sxx)*ry+18.*rxy*ryy*sx+9.*sy*rxx*ryy+18.*sy*rxy**2)*rx+3.*ry*(1/3.*ry**2.*sxxx+(3.*rxx*sxy+rxxx*sy+3.*rxxy*sx+3.*rxy*sxx)*ry+3.*sx*rxx*ryy+6.*sx*rxy**2+6.*rxy*rxx*sy)
cuxxxyyy410 = (3.*ry*syy+3.*ryy*sy)*rx**3+9.*ry*(2.*rxy*sy+ry*sxy+ryy*sx)*rx**2+3.*ry**2.*(3.*rxx*sy+6.*rxy*sx+ry*sxx)*rx+3.*rxx*sx*ry**3
cuxxxyyy510 = 3.*rx**3.*ry**2.*sy+3.*rx**2.*ry**3.*sx
cuxxxyyy020 = 3.*sx*sxxyyy+3.*sxx*sxyyy+sxxx*syyy+3.*sxxxy*syy+3.*sxxxyy*sy+9.*sxxy*sxyy+9.*sxxyy*sxy
cuxxxyyy120 = 3.*rxyyy*sx**2+(6.*rx*sxyyy+3.*rxx*syyy+9.*rxxy*syy+9.*rxxyy*sy+18.*rxy*sxyy+18.*rxyy*sxy+9.*ry*sxxyy+9.*ryy*sxxy+3.*ryyy*sxx)*sx+3.*rxxxy*sy**2+(9.*rx*sxxyy+9.*rxx*sxyy+3.*rxxx*syy+18.*rxxy*sxy+18.*rxy*sxxy+9.*rxyy*sxx+6.*ry*sxxxy+3.*ryy*sxxx)*sy+18.*sxy**2.*rxy+(18.*rx*sxyy+9.*rxx*syy+18.*ry*sxxy+9.*ryy*sxx)*sxy+(3.*rx*syyy+9.*rxy*syy+9.*ry*sxyy)*sxx+(9.*rx*sxxy+3.*ry*sxxx)*syy
cuxxxyyy220 = (3.*sx*syyy+9.*sxy*syy+9.*sxyy*sy)*rx**2+((18.*sx*sxyy+9.*sxx*syy+18.*sxxy*sy+18.*sxy**2)*ry+3.*ryyy*sx**2+(18.*rxy*syy+18.*rxyy*sy+18.*ryy*sxy)*sx+9.*sy*(rxx*syy+rxxy*sy+4.*rxy*sxy+ryy*sxx))*rx+(9.*sx*sxxy+9.*sxx*sxy+3.*sxxx*sy)*ry**2+(9.*rxyy*sx**2+(9.*rxx*syy+18.*rxxy*sy+36.*rxy*sxy+9.*ryy*sxx)*sx+3.*sy*(6.*rxx*sxy+rxxx*sy+6.*rxy*sxx))*ry+9.*rxy*ryy*sx**2+9.*sy*(rxx*ryy+2.*rxy**2)*sx+9.*rxy*rxx*sy**2
cuxxxyyy320 = 3.*sy*syy*rx**3+((9.*sx*syy+18.*sxy*sy)*ry+9.*ryy*sy*sx+9.*rxy*sy**2)*rx**2+18.*ry*((sxy*sx+1/2.*sy*sxx)*ry+1/2.*ryy*sx**2+2.*rxy*sy*sx+1/2.*rxx*sy**2)*rx+3.*ry**2.*sx*(3.*rxx*sy+3.*rxy*sx+ry*sxx)
cuxxxyyy420 = 3.*rx**3.*ry*sy**2+9.*rx**2.*ry**2.*sx*sy+3.*rx*ry**3.*sx**2
cuxxxyyy030 = 3.*sxxxy*sy**2+(3.*sxx*syyy+9.*sxxy*syy+9.*sxxyy*sy+18.*sxy*sxyy)*sx+3.*sx**2.*sxyyy+9.*sxx*syy*sxy+(9.*sxx*sxyy+3.*sxxx*syy+18.*sxxy*sxy)*sy+6.*sxy**3
cuxxxyyy130 = ryyy*sx**3+(3.*rx*syyy+9.*rxy*syy+9.*rxyy*sy+9.*ry*sxyy+9.*ryy*sxy)*sx**2+(9.*rxxy*sy**2+(18.*rx*sxyy+9.*rxx*syy+36.*rxy*sxy+18.*ry*sxxy+9.*ryy*sxx)*sy+18.*rx*syy*sxy+9.*syy*ry*sxx+18.*ry*sxy**2)*sx+9.*(1/9.*rxxx*sy**2+(rx*sxxy+1/3.*ry*sxxx+rxx*sxy+rxy*sxx)*sy+rx*syy*sxx+2.*rx*sxy**2+2.*ry*sxx*sxy)*sy
cuxxxyyy230 = 3.*ryy*ry*sx**3+((9.*rx*ryy+18.*rxy*ry)*sy+9.*rx*syy*ry+9.*ry**2.*sxy)*sx**2+9.*sy*((2.*rx*rxy+rxx*ry)*sy+rx**2.*syy+4.*ry*sxy*rx+ry**2.*sxx)*sx+9.*(rx*sxy+ry*sxx+1/3.*rxx*sy)*sy**2.*rx
cuxxxyyy330 = rx**3.*sy**3+9.*rx**2.*ry*sx*sy**2+9.*rx*ry**2.*sx**2.*sy+ry**3.*sx**3
cuxxxyyy040 = sx**3.*syyy+(9.*sxy*syy+9.*sxyy*sy)*sx**2+9.*sy*(sxx*syy+sxxy*sy+2.*sxy**2)*sx+sxxx*sy**3+9.*sxx*sxy*sy**2
cuxxxyyy140 = (3.*ry*syy+3.*ryy*sy)*sx**3+9.*sy*(rx*syy+rxy*sy+2.*ry*sxy)*sx**2+(3.*rxx*sy**3+(18.*rx*sxy+9.*ry*sxx)*sy**2)*sx+3.*rx*sxx*sy**3
cuxxxyyy240 = 3.*rx**2.*sx*sy**3+9.*rx*ry*sx**2.*sy**2+3.*ry**2.*sx**3.*sy
cuxxxyyy050 = 3.*sx**3.*sy*syy+9.*sx**2.*sxy*sy**2+3.*sx*sxx*sy**3
cuxxxyyy150 = 3.*rx*sx**2.*sy**3+3.*ry*sx**3.*sy**2
cuxxxyyy060 = sx**3.*sy**3
                      ! uxxxyyy = cuxxxyyy100*ur+cuxxxyyy200*urr+cuxxxyyy300*urrr+cuxxxyyy400*urrrr+cuxxxyyy500*urrrrr+cuxxxyyy600*urrrrrr+cuxxxyyy010*us+cuxxxyyy110*urs+cuxxxyyy210*urrs+cuxxxyyy310*urrrs+cuxxxyyy410*urrrrs+cuxxxyyy510*urrrrrs+cuxxxyyy020*uss+cuxxxyyy120*urss+cuxxxyyy220*urrss+cuxxxyyy320*urrrss+cuxxxyyy420*urrrrss+cuxxxyyy030*usss+cuxxxyyy130*ursss+cuxxxyyy230*urrsss+cuxxxyyy330*urrrsss+cuxxxyyy040*ussss+cuxxxyyy140*urssss+cuxxxyyy240*urrssss+cuxxxyyy050*usssss+cuxxxyyy150*ursssss+cuxxxyyy060*ussssss
                      ! ------ Coefficients in expansion for uxxyyyy ------
cuxxyyyy100 = rxxyyyy
cuxxyyyy200 = 2.*rx*rxyyyy+rxx*ryyyy+4.*rxxy*ryyy+6.*rxxyy*ryy+4.*rxxyyy*ry+8.*rxy*rxyyy+6.*rxyy**2
cuxxyyyy300 = 6.*rxxyy*ry**2+(8.*rx*rxyyy+4.*rxx*ryyy+12.*rxxy*ryy+24.*rxy*rxyy)*ry+3.*rxx*ryy**2+(12.*rx*rxyy+12.*rxy**2)*ryy+ryyyy*rx**2+8.*rx*rxy*ryyy
cuxxyyyy400 = (12.*rx*rxyy+6.*rxx*ryy+12.*rxy**2)*ry**2+(4.*rx**2.*ryyy+24.*rx*rxy*ryy)*ry+4.*rxxy*ry**3+3.*ryy**2.*rx**2
cuxxyyyy500 = 6.*rx**2.*ry**2.*ryy+8.*rx*rxy*ry**3+rxx*ry**4
cuxxyyyy600 = rx**2.*ry**4
cuxxyyyy010 = sxxyyyy
cuxxyyyy110 = 2.*rx*sxyyyy+rxx*syyyy+4.*rxxy*syyy+6.*rxxyy*syy+4.*rxxyyy*sy+8.*rxy*sxyyy+12.*rxyy*sxyy+8.*rxyyy*sxy+2.*rxyyyy*sx+4.*ry*sxxyyy+6.*ryy*sxxyy+4.*ryyy*sxxy+ryyyy*sxx
cuxxyyyy210 = 6.*ry**2.*sxxyy+(8.*rx*sxyyy+4.*rxx*syyy+12.*rxxy*syy+12.*rxxyy*sy+24.*rxy*sxyy+24.*rxyy*sxy+8.*rxyyy*sx+12.*ryy*sxxy+4.*ryyy*sxx)*ry+rx**2.*syyyy+(8.*rxy*syyy+12.*rxyy*syy+8.*rxyyy*sy+12.*ryy*sxyy+8.*ryyy*sxy+2.*ryyyy*sx)*rx+3.*ryy**2.*sxx+(6.*rxx*syy+12.*rxxy*sy+24.*rxy*sxy+12.*rxyy*sx)*ryy+12.*rxy**2.*syy+(24.*rxyy*sy+8.*ryyy*sx)*rxy+4.*rxx*ryyy*sy
cuxxyyyy310 = 4.*ry**3.*sxxy+(12.*rx*sxyy+6.*rxx*syy+12.*rxxy*sy+24.*rxy*sxy+12.*rxyy*sx+6.*ryy*sxx)*ry**2+(4.*rx**2.*syyy+(24.*rxy*syy+24.*rxyy*sy+24.*ryy*sxy+8.*ryyy*sx)*rx+(12.*rxx*sy+24.*rxy*sx)*ryy+24.*sy*rxy**2)*ry+6.*((syy*ryy+2/3.*ryyy*sy)*rx+ryy**2.*sx+4.*rxy*ryy*sy)*rx
cuxxyyyy410 = 6.*rx**2.*ry**2.*syy+12.*rx**2.*ry*ryy*sy+24.*rx*rxy*ry**2.*sy+8.*rx*ry**3.*sxy+12.*rx*ry**2.*ryy*sx+4.*rxx*ry**3.*sy+8.*rxy*ry**3.*sx+ry**4.*sxx
cuxxyyyy510 = 4.*rx**2.*ry**3.*sy+2.*rx*ry**4.*sx
cuxxyyyy020 = 2.*sx*sxyyyy+sxx*syyyy+4.*sxxy*syyy+6.*sxxyy*syy+4.*sxxyyy*sy+8.*sxy*sxyyy+6.*sxyy**2
cuxxyyyy120 = 6.*rxxyy*sy**2+(8.*rx*sxyyy+4.*rxx*syyy+12.*rxxy*syy+24.*rxy*sxyy+24.*rxyy*sxy+8.*rxyyy*sx+12.*ry*sxxyy+12.*ryy*sxxy+4.*ryyy*sxx)*sy+ryyyy*sx**2+(2.*rx*syyyy+8.*rxy*syyy+12.*rxyy*syy+8.*ry*sxyyy+12.*ryy*sxyy+8.*ryyy*sxy)*sx+3.*rxx*syy**2+(12.*rx*sxyy+24.*rxy*sxy+12.*ry*sxxy+6.*ryy*sxx)*syy+12.*ryy*sxy**2+(8.*rx*syyy+24.*ry*sxyy)*sxy+4.*ry*sxx*syyy
cuxxyyyy220 = (12.*sx*sxyy+6.*sxx*syy+12.*sxxy*sy+12.*sxy**2)*ry**2+(12.*rxxy*sy**2+(24.*rx*sxyy+12.*rxx*syy+48.*rxy*sxy+24.*rxyy*sx+12.*ryy*sxx)*sy+(8.*sx*syyy+24.*sxy*syy)*rx+24.*(rxy*syy+1/6.*ryyy*sx+ryy*sxy)*sx)*ry+(12.*rx*rxyy+6.*rxx*ryy+12.*rxy**2)*sy**2+(4.*rx**2.*syyy+(24.*rxy*syy+24.*ryy*sxy+8.*ryyy*sx)*rx+24.*rxy*ryy*sx)*sy+3.*syy**2.*rx**2+12.*ryy*syy*sx*rx+3.*ryy**2.*sx**2
cuxxyyyy320 = (8.*sx*sxy+4.*sxx*sy)*ry**3+(6.*rxx*sy**2+(24.*rx*sxy+24.*rxy*sx)*sy+12.*sx*syy*rx+6.*ryy*sx**2)*ry**2+12.*rx*sy*(rx*syy+2.*rxy*sy+2.*ryy*sx)*ry+6.*ryy*rx**2.*sy**2
cuxxyyyy420 = 6.*rx**2.*ry**2.*sy**2+8.*rx*ry**3.*sx*sy+ry**4.*sx**2
cuxxyyyy030 = 6.*sxxyy*sy**2+(8.*sx*sxyyy+4.*sxx*syyy+12.*sxxy*syy+24.*sxy*sxyy)*sy+3.*sxx*syy**2+(12.*sx*sxyy+12.*sxy**2)*syy+sx**2.*syyyy+8.*sx*sxy*syyy
cuxxyyyy130 = 4.*rxxy*sy**3+(12.*rx*sxyy+6.*rxx*syy+24.*rxy*sxy+12.*rxyy*sx+12.*ry*sxxy+6.*ryy*sxx)*sy**2+(4.*ryyy*sx**2+(8.*rx*syyy+24.*rxy*syy+24.*ry*sxyy+24.*ryy*sxy)*sx+(24.*rx*sxy+12.*ry*sxx)*syy+24.*ry*sxy**2)*sy+(4.*ry*syyy+6.*ryy*syy)*sx**2+(6.*rx*syy**2+24.*ry*sxy*syy)*sx
cuxxyyyy230 = (8.*rx*rxy+4.*rxx*ry)*sy**3+(6.*ry**2.*sxx+(24.*rx*sxy+24.*rxy*sx)*ry+6.*rx**2.*syy+12.*rx*sx*ryy)*sy**2+24.*ry*sx*(rx*syy+ry*sxy+1/2.*ryy*sx)*sy+6.*ry**2.*sx**2.*syy
cuxxyyyy330 = 4.*rx**2.*ry*sy**3+12.*rx*ry**2.*sx*sy**2+4.*ry**3.*sx**2.*sy
cuxxyyyy040 = (12.*sx*sxyy+6.*sxx*syy+12.*sxy**2)*sy**2+4.*sxxy*sy**3+3.*syy**2.*sx**2+(4.*sx**2.*syyy+24.*sx*sxy*syy)*sy
cuxxyyyy140 = 12.*rx*sx*sy**2.*syy+8.*rx*sxy*sy**3+rxx*sy**4+8.*rxy*sx*sy**3+12.*ry*sx**2.*sy*syy+24.*ry*sx*sxy*sy**2+4.*ry*sxx*sy**3+6.*ryy*sx**2.*sy**2
cuxxyyyy240 = rx**2.*sy**4+8.*rx*ry*sx*sy**3+6.*ry**2.*sx**2.*sy**2
cuxxyyyy050 = 6.*sx**2.*sy**2.*syy+8.*sx*sxy*sy**3+sxx*sy**4
cuxxyyyy150 = 2.*rx*sx*sy**4+4.*ry*sx**2.*sy**3
cuxxyyyy060 = sx**2.*sy**4
                      ! uxxyyyy = cuxxyyyy100*ur+cuxxyyyy200*urr+cuxxyyyy300*urrr+cuxxyyyy400*urrrr+cuxxyyyy500*urrrrr+cuxxyyyy600*urrrrrr+cuxxyyyy010*us+cuxxyyyy110*urs+cuxxyyyy210*urrs+cuxxyyyy310*urrrs+cuxxyyyy410*urrrrs+cuxxyyyy510*urrrrrs+cuxxyyyy020*uss+cuxxyyyy120*urss+cuxxyyyy220*urrss+cuxxyyyy320*urrrss+cuxxyyyy420*urrrrss+cuxxyyyy030*usss+cuxxyyyy130*ursss+cuxxyyyy230*urrsss+cuxxyyyy330*urrrsss+cuxxyyyy040*ussss+cuxxyyyy140*urssss+cuxxyyyy240*urrssss+cuxxyyyy050*usssss+cuxxyyyy150*ursssss+cuxxyyyy060*ussssss
                      ! ------ Coefficients in expansion for uxyyyyy ------
cuxyyyyy100 = rxyyyyy
cuxyyyyy200 = rx*ryyyyy+5.*rxy*ryyyy+10.*rxyy*ryyy+10.*rxyyy*ryy+5.*rxyyyy*ry
cuxyyyyy300 = 10.*rxyyy*ry**2+(5.*rx*ryyyy+20.*rxy*ryyy+30.*rxyy*ryy)*ry+15.*rxy*ryy**2+10.*ryyy*rx*ryy
cuxyyyyy400 = 10.*rx*ry**2.*ryyy+15.*rx*ry*ryy**2+30.*rxy*ry**2.*ryy+10.*rxyy*ry**3
cuxyyyyy500 = 10.*rx*ry**3.*ryy+5.*rxy*ry**4
cuxyyyyy600 = rx*ry**5
cuxyyyyy010 = sxyyyyy
cuxyyyyy110 = rx*syyyyy+5.*rxy*syyyy+10.*rxyy*syyy+10.*rxyyy*syy+5.*rxyyyy*sy+5.*ry*sxyyyy+10.*ryy*sxyyy+10.*ryyy*sxyy+5.*ryyyy*sxy+ryyyyy*sx
cuxyyyyy210 = 10.*ry**2.*sxyyy+(5.*rx*syyyy+20.*rxy*syyy+30.*rxyy*syy+20.*rxyyy*sy+30.*ryy*sxyy+20.*ryyy*sxy+5.*ryyyy*sx)*ry+15.*ryy**2.*sxy+(10.*rx*syyy+30.*rxy*syy+30.*rxyy*sy+10.*ryyy*sx)*ryy+(10.*rx*syy+20.*rxy*sy)*ryyy+5.*rx*ryyyy*sy
cuxyyyyy310 = 10.*ry**3.*sxyy+(10.*rx*syyy+30.*rxy*syy+30.*rxyy*sy+30.*ryy*sxy+10.*ryyy*sx)*ry**2+(15.*ryy**2.*sx+(30.*rx*syy+60.*rxy*sy)*ryy+20.*rx*ryyy*sy)*ry+15.*ryy**2.*rx*sy
cuxyyyyy410 = 10.*rx*ry**3.*syy+30.*rx*ry**2.*ryy*sy+20.*rxy*ry**3.*sy+5.*ry**4.*sxy+10.*ry**3.*ryy*sx
cuxyyyyy510 = 5.*rx*ry**4.*sy+ry**5.*sx
cuxyyyyy020 = sx*syyyyy+5.*sxy*syyyy+10.*sxyy*syyy+10.*sxyyy*syy+5.*sxyyyy*sy
cuxyyyyy120 = 10.*rxyyy*sy**2+(5.*rx*syyyy+20.*rxy*syyy+30.*rxyy*syy+20.*ry*sxyyy+30.*ryy*sxyy+20.*ryyy*sxy+5.*ryyyy*sx)*sy+15.*rxy*syy**2+(10.*rx*syyy+30.*ry*sxyy+30.*ryy*sxy+10.*ryyy*sx)*syy+(20.*ry*sxy+10.*ryy*sx)*syyy+5.*ry*sx*syyyy
cuxyyyyy220 = (10.*sx*syyy+30.*sxy*syy+30.*sxyy*sy)*ry**2+(30.*rxyy*sy**2+(20.*rx*syyy+60.*rxy*syy+60.*ryy*sxy+20.*ryyy*sx)*sy+30.*ryy*syy*sx+15.*syy**2.*rx)*ry+10.*sy*((rx*ryyy+3.*rxy*ryy)*sy+3.*ryy*syy*rx+3/2.*ryy**2.*sx)
cuxyyyyy320 = 30.*rx*ry**2.*sy*syy+30.*rx*ry*ryy*sy**2+30.*rxy*ry**2.*sy**2+10.*ry**3.*sx*syy+20.*ry**3.*sxy*sy+30.*ry**2.*ryy*sx*sy
cuxyyyyy420 = 10.*rx*ry**3.*sy**2+5.*ry**4.*sx*sy
cuxyyyyy030 = 10.*sxyyy*sy**2+(5.*sx*syyyy+20.*sxy*syyy+30.*sxyy*syy)*sy+15.*sxy*syy**2+10.*sx*syyy*syy
cuxyyyyy130 = 10.*rxyy*sy**3+(10.*rx*syyy+30.*rxy*syy+30.*ry*sxyy+30.*ryy*sxy+10.*ryyy*sx)*sy**2+(15.*syy**2.*rx+(60.*ry*sxy+30.*ryy*sx)*syy+20.*syyy*ry*sx)*sy+15.*ry*sx*syy**2
cuxyyyyy230 = 30.*rx*ry*sy**2.*syy+10.*rx*ryy*sy**3+20.*rxy*ry*sy**3+30.*ry**2.*sx*sy*syy+30.*ry**2.*sxy*sy**2+30.*ry*ryy*sx*sy**2
cuxyyyyy330 = 10.*rx*ry**2.*sy**3+10.*ry**3.*sx*sy**2
cuxyyyyy040 = 10.*sx*sy**2.*syyy+15.*sx*sy*syy**2+30.*sxy*sy**2.*syy+10.*sxyy*sy**3
cuxyyyyy140 = 10.*rx*sy**3.*syy+5.*rxy*sy**4+30.*ry*sx*sy**2.*syy+20.*ry*sxy*sy**3+10.*ryy*sx*sy**3
cuxyyyyy240 = 5.*rx*ry*sy**4+10.*ry**2.*sx*sy**3
cuxyyyyy050 = 10.*sx*sy**3.*syy+5.*sxy*sy**4
cuxyyyyy150 = rx*sy**5+5.*ry*sx*sy**4
cuxyyyyy060 = sx*sy**5
                      ! uxyyyyy = cuxyyyyy100*ur+cuxyyyyy200*urr+cuxyyyyy300*urrr+cuxyyyyy400*urrrr+cuxyyyyy500*urrrrr+cuxyyyyy600*urrrrrr+cuxyyyyy010*us+cuxyyyyy110*urs+cuxyyyyy210*urrs+cuxyyyyy310*urrrs+cuxyyyyy410*urrrrs+cuxyyyyy510*urrrrrs+cuxyyyyy020*uss+cuxyyyyy120*urss+cuxyyyyy220*urrss+cuxyyyyy320*urrrss+cuxyyyyy420*urrrrss+cuxyyyyy030*usss+cuxyyyyy130*ursss+cuxyyyyy230*urrsss+cuxyyyyy330*urrrsss+cuxyyyyy040*ussss+cuxyyyyy140*urssss+cuxyyyyy240*urrssss+cuxyyyyy050*usssss+cuxyyyyy150*ursssss+cuxyyyyy060*ussssss
                      ! ------ Coefficients in expansion for uyyyyyy ------
cuyyyyyy100 = ryyyyyy
cuyyyyyy200 = 6.*ry*ryyyyy+15.*ryy*ryyyy+10.*ryyy**2
cuyyyyyy300 = 15.*ry**2.*ryyyy+60.*ry*ryy*ryyy+15.*ryy**3
cuyyyyyy400 = 20.*ry**3.*ryyy+45.*ry**2.*ryy**2
cuyyyyyy500 = 15.*ry**4.*ryy
cuyyyyyy600 = ry**6
cuyyyyyy010 = syyyyyy
cuyyyyyy110 = 6.*ry*syyyyy+15.*ryy*syyyy+20.*ryyy*syyy+15.*ryyyy*syy+6.*ryyyyy*sy
cuyyyyyy210 = 15.*ry**2.*syyyy+(60.*ryy*syyy+60.*ryyy*syy+30.*ryyyy*sy)*ry+60.*ryyy*ryy*sy+45.*syy*ryy**2
cuyyyyyy310 = 20.*ry**3.*syyy+(90.*ryy*syy+60.*ryyy*sy)*ry**2+90.*ryy**2.*sy*ry
cuyyyyyy410 = 15.*ry**4.*syy+60.*ry**3.*ryy*sy
cuyyyyyy510 = 6.*sy*ry**5
cuyyyyyy020 = 6.*sy*syyyyy+15.*syy*syyyy+10.*syyy**2
cuyyyyyy120 = 15.*ryyyy*sy**2+(30.*ry*syyyy+60.*ryy*syyy+60.*ryyy*syy)*sy+60.*ry*syy*syyy+45.*syy**2.*ryy
cuyyyyyy220 = (60.*ry*ryyy+45.*ryy**2)*sy**2+(60.*ry**2.*syyy+180.*ry*ryy*syy)*sy+45.*ry**2.*syy**2
cuyyyyyy320 = 60.*sy*(ry*syy+3/2.*ryy*sy)*ry**2
cuyyyyyy420 = 15.*ry**4.*sy**2
cuyyyyyy030 = 15.*sy**2.*syyyy+60.*sy*syy*syyy+15.*syy**3
cuyyyyyy130 = 20.*ryyy*sy**3+(60.*ry*syyy+90.*ryy*syy)*sy**2+90.*ry*sy*syy**2
cuyyyyyy230 = 90.*sy**2.*ry*(ry*syy+2/3.*ryy*sy)
cuyyyyyy330 = 20.*ry**3.*sy**3
cuyyyyyy040 = 20.*sy**3.*syyy+45.*sy**2.*syy**2
cuyyyyyy140 = 60.*ry*sy**3.*syy+15.*ryy*sy**4
cuyyyyyy240 = 15.*ry**2.*sy**4
cuyyyyyy050 = 15.*sy**4.*syy
cuyyyyyy150 = 6.*ry*sy**5
cuyyyyyy060 = sy**6
                      ! uyyyyyy = cuyyyyyy100*ur+cuyyyyyy200*urr+cuyyyyyy300*urrr+cuyyyyyy400*urrrr+cuyyyyyy500*urrrrr+cuyyyyyy600*urrrrrr+cuyyyyyy010*us+cuyyyyyy110*urs+cuyyyyyy210*urrs+cuyyyyyy310*urrrs+cuyyyyyy410*urrrrs+cuyyyyyy510*urrrrrs+cuyyyyyy020*uss+cuyyyyyy120*urss+cuyyyyyy220*urrss+cuyyyyyy320*urrrss+cuyyyyyy420*urrrrss+cuyyyyyy030*usss+cuyyyyyy130*ursss+cuyyyyyy230*urrsss+cuyyyyyy330*urrrsss+cuyyyyyy040*ussss+cuyyyyyy140*urssss+cuyyyyyy240*urrssss+cuyyyyyy050*usssss+cuyyyyyy150*ursssss+cuyyyyyy060*ussssss
                                            uxxxxxx = cuxxxxxx100*ur+cuxxxxxx200*urr+cuxxxxxx300*urrr+cuxxxxxx400*urrrr+cuxxxxxx500*urrrrr+cuxxxxxx600*urrrrrr+cuxxxxxx010*us+cuxxxxxx110*urs+cuxxxxxx210*urrs+cuxxxxxx310*urrrs+cuxxxxxx410*urrrrs+cuxxxxxx510*urrrrrs+cuxxxxxx020*uss+cuxxxxxx120*urss+cuxxxxxx220*urrss+cuxxxxxx320*urrrss+cuxxxxxx420*urrrrss+cuxxxxxx030*usss+cuxxxxxx130*ursss+cuxxxxxx230*urrsss+cuxxxxxx330*urrrsss+cuxxxxxx040*ussss+cuxxxxxx140*urssss+cuxxxxxx240*urrssss+cuxxxxxx050*usssss+cuxxxxxx150*ursssss+cuxxxxxx060*ussssss
                      ! uxxxxxy = cuxxxxxy100*ur+cuxxxxxy200*urr+cuxxxxxy300*urrr+cuxxxxxy400*urrrr+cuxxxxxy500*urrrrr+cuxxxxxy600*urrrrrr+cuxxxxxy010*us+cuxxxxxy110*urs+cuxxxxxy210*urrs+cuxxxxxy310*urrrs+cuxxxxxy410*urrrrs+cuxxxxxy510*urrrrrs+cuxxxxxy020*uss+cuxxxxxy120*urss+cuxxxxxy220*urrss+cuxxxxxy320*urrrss+cuxxxxxy420*urrrrss+cuxxxxxy030*usss+cuxxxxxy130*ursss+cuxxxxxy230*urrsss+cuxxxxxy330*urrrsss+cuxxxxxy040*ussss+cuxxxxxy140*urssss+cuxxxxxy240*urrssss+cuxxxxxy050*usssss+cuxxxxxy150*ursssss+cuxxxxxy060*ussssss
                                            uxxxxyy = cuxxxxyy100*ur+cuxxxxyy200*urr+cuxxxxyy300*urrr+cuxxxxyy400*urrrr+cuxxxxyy500*urrrrr+cuxxxxyy600*urrrrrr+cuxxxxyy010*us+cuxxxxyy110*urs+cuxxxxyy210*urrs+cuxxxxyy310*urrrs+cuxxxxyy410*urrrrs+cuxxxxyy510*urrrrrs+cuxxxxyy020*uss+cuxxxxyy120*urss+cuxxxxyy220*urrss+cuxxxxyy320*urrrss+cuxxxxyy420*urrrrss+cuxxxxyy030*usss+cuxxxxyy130*ursss+cuxxxxyy230*urrsss+cuxxxxyy330*urrrsss+cuxxxxyy040*ussss+cuxxxxyy140*urssss+cuxxxxyy240*urrssss+cuxxxxyy050*usssss+cuxxxxyy150*ursssss+cuxxxxyy060*ussssss
                      ! uxxxyyy = cuxxxyyy100*ur+cuxxxyyy200*urr+cuxxxyyy300*urrr+cuxxxyyy400*urrrr+cuxxxyyy500*urrrrr+cuxxxyyy600*urrrrrr+cuxxxyyy010*us+cuxxxyyy110*urs+cuxxxyyy210*urrs+cuxxxyyy310*urrrs+cuxxxyyy410*urrrrs+cuxxxyyy510*urrrrrs+cuxxxyyy020*uss+cuxxxyyy120*urss+cuxxxyyy220*urrss+cuxxxyyy320*urrrss+cuxxxyyy420*urrrrss+cuxxxyyy030*usss+cuxxxyyy130*ursss+cuxxxyyy230*urrsss+cuxxxyyy330*urrrsss+cuxxxyyy040*ussss+cuxxxyyy140*urssss+cuxxxyyy240*urrssss+cuxxxyyy050*usssss+cuxxxyyy150*ursssss+cuxxxyyy060*ussssss
                                            uxxyyyy = cuxxyyyy100*ur+cuxxyyyy200*urr+cuxxyyyy300*urrr+cuxxyyyy400*urrrr+cuxxyyyy500*urrrrr+cuxxyyyy600*urrrrrr+cuxxyyyy010*us+cuxxyyyy110*urs+cuxxyyyy210*urrs+cuxxyyyy310*urrrs+cuxxyyyy410*urrrrs+cuxxyyyy510*urrrrrs+cuxxyyyy020*uss+cuxxyyyy120*urss+cuxxyyyy220*urrss+cuxxyyyy320*urrrss+cuxxyyyy420*urrrrss+cuxxyyyy030*usss+cuxxyyyy130*ursss+cuxxyyyy230*urrsss+cuxxyyyy330*urrrsss+cuxxyyyy040*ussss+cuxxyyyy140*urssss+cuxxyyyy240*urrssss+cuxxyyyy050*usssss+cuxxyyyy150*ursssss+cuxxyyyy060*ussssss
                      ! uxyyyyy = cuxyyyyy100*ur+cuxyyyyy200*urr+cuxyyyyy300*urrr+cuxyyyyy400*urrrr+cuxyyyyy500*urrrrr+cuxyyyyy600*urrrrrr+cuxyyyyy010*us+cuxyyyyy110*urs+cuxyyyyy210*urrs+cuxyyyyy310*urrrs+cuxyyyyy410*urrrrs+cuxyyyyy510*urrrrrs+cuxyyyyy020*uss+cuxyyyyy120*urss+cuxyyyyy220*urrss+cuxyyyyy320*urrrss+cuxyyyyy420*urrrrss+cuxyyyyy030*usss+cuxyyyyy130*ursss+cuxyyyyy230*urrsss+cuxyyyyy330*urrrsss+cuxyyyyy040*ussss+cuxyyyyy140*urssss+cuxyyyyy240*urrssss+cuxyyyyy050*usssss+cuxyyyyy150*ursssss+cuxyyyyy060*ussssss
                                            uyyyyyy = cuyyyyyy100*ur+cuyyyyyy200*urr+cuyyyyyy300*urrr+cuyyyyyy400*urrrr+cuyyyyyy500*urrrrr+cuyyyyyy600*urrrrrr+cuyyyyyy010*us+cuyyyyyy110*urs+cuyyyyyy210*urrs+cuyyyyyy310*urrrs+cuyyyyyy410*urrrrs+cuyyyyyy510*urrrrrs+cuyyyyyy020*uss+cuyyyyyy120*urss+cuyyyyyy220*urrss+cuyyyyyy320*urrrss+cuyyyyyy420*urrrrss+cuyyyyyy030*usss+cuyyyyyy130*ursss+cuyyyyyy230*urrsss+cuyyyyyy330*urrrsss+cuyyyyyy040*ussss+cuyyyyyy140*urssss+cuyyyyyy240*urrssss+cuyyyyyy050*usssss+cuyyyyyy150*ursssss+cuyyyyyy060*ussssss
                      ! if( i1.eq.5 .and. i2.eq.6 )then  
                      !   write(*,'(" (i1,i2)=(",2i3,") d420i,d620i,d440i=",4(1pe12.4,1x))') i1,i2,d420i,d620i,d440i
                      !   write(*,'(" (i1,i2)=(",2i3,") urrrrss,cuxxxxy420=",4(1pe12.4,1x))') i1,i2,urrrrss,cuxxxxyy420
                      ! end if 
                      ! ! ----- START SEVENTH DERIVATIVES -----
                      ! do m1=0,nd-1
                      ! do m2=0,nd-1
                      !   ! ---- 6th parameteric derivatives of the metrics ----
                      !   rxrrrrrra(m1,m2) = rx600(i1,i2,i3,m1,m2) - (dx(0)**2/4.)*rx800i(m1,m2)
                      !   rxssssssa(m1,m2) = rx060(i1,i2,i3,m1,m2) - (dx(1)**2/4.)*rx080i(m1,m2)  
                      !   rxrrrrrsa(m1,m2) = rx510i(m1,m2) - (dx(0)**2/3.)*rx710i(m1,m2) -  (dx(1)**2/6.)*rx530i(m1,m2)
                      !   rxrsssssa(m1,m2) = rx150i(m1,m2) - (dx(1)**2/3.)*rx170i(m1,m2) -  (dx(0)**2/6.)*rx350i(m1,m2)
                      !   rxrrrrssa(m1,m2) = rx420i(m1,m2) - (dx(0)**2/6.)*rx620i(m1,m2) - (dx(1)**2/12.)*rx440i(m1,m2)
                      !   rxrrssssa(m1,m2) = rx240i(m1,m2) - (dx(1)**2/6.)*rx260i(m1,m2) - (dx(0)**2/12.)*rx440i(m1,m2)
                      !   rxrrrsssa(m1,m2) = rx330i(m1,m2) - (dx(0)**2/4.)*rx530i(m1,m2) -  (dx(1)**2/4.)*rx350i(m1,m2)
                      !   ! ---- sixth spatial derivatives of the metrics ----
                      !   rxxxxxxxa(m1,m2) = cuxxxxxx100*rxr(m1,m2,0)+cuxxxxxx200*rxrra(m1,m2)+cuxxxxxx300*rxrrra(m1,m2)+cuxxxxxx400*rxrrrra(m1,m2)+cuxxxxxx500*rxrrrrra(m1,m2)+cuxxxxxx600*rxrrrrrra(m1,m2)+cuxxxxxx010*rxr(m1,m2,1)+cuxxxxxx110*rxrsa(m1,m2)+cuxxxxxx210*rxrrsa(m1,m2)+cuxxxxxx310*rxrrrsa(m1,m2)+cuxxxxxx410*rxrrrrsa(m1,m2)+cuxxxxxx510*rxrrrrrsa(m1,m2)+cuxxxxxx020*rxssa(m1,m2)+cuxxxxxx120*rxrssa(m1,m2)+cuxxxxxx220*rxrrssa(m1,m2)+cuxxxxxx320*rxrrrssa(m1,m2)+cuxxxxxx420*rxrrrrssa(m1,m2)+cuxxxxxx030*rxsssa(m1,m2)+cuxxxxxx130*rxrsssa(m1,m2)+cuxxxxxx230*rxrrsssa(m1,m2)+cuxxxxxx330*rxrrrsssa(m1,m2)+cuxxxxxx040*rxssssa(m1,m2)+cuxxxxxx140*rxrssssa(m1,m2)+cuxxxxxx240*rxrrssssa(m1,m2)+cuxxxxxx050*rxsssssa(m1,m2)+cuxxxxxx150*rxrsssssa(m1,m2)+cuxxxxxx060*rxssssssa(m1,m2)
                      !   rxxxxxxya(m1,m2) = cuxxxxxy100*rxr(m1,m2,0)+cuxxxxxy200*rxrra(m1,m2)+cuxxxxxy300*rxrrra(m1,m2)+cuxxxxxy400*rxrrrra(m1,m2)+cuxxxxxy500*rxrrrrra(m1,m2)+cuxxxxxy600*rxrrrrrra(m1,m2)+cuxxxxxy010*rxr(m1,m2,1)+cuxxxxxy110*rxrsa(m1,m2)+cuxxxxxy210*rxrrsa(m1,m2)+cuxxxxxy310*rxrrrsa(m1,m2)+cuxxxxxy410*rxrrrrsa(m1,m2)+cuxxxxxy510*rxrrrrrsa(m1,m2)+cuxxxxxy020*rxssa(m1,m2)+cuxxxxxy120*rxrssa(m1,m2)+cuxxxxxy220*rxrrssa(m1,m2)+cuxxxxxy320*rxrrrssa(m1,m2)+cuxxxxxy420*rxrrrrssa(m1,m2)+cuxxxxxy030*rxsssa(m1,m2)+cuxxxxxy130*rxrsssa(m1,m2)+cuxxxxxy230*rxrrsssa(m1,m2)+cuxxxxxy330*rxrrrsssa(m1,m2)+cuxxxxxy040*rxssssa(m1,m2)+cuxxxxxy140*rxrssssa(m1,m2)+cuxxxxxy240*rxrrssssa(m1,m2)+cuxxxxxy050*rxsssssa(m1,m2)+cuxxxxxy150*rxrsssssa(m1,m2)+cuxxxxxy060*rxssssssa(m1,m2)
                      !   rxxxxxyya(m1,m2) = cuxxxxyy100*rxr(m1,m2,0)+cuxxxxyy200*rxrra(m1,m2)+cuxxxxyy300*rxrrra(m1,m2)+cuxxxxyy400*rxrrrra(m1,m2)+cuxxxxyy500*rxrrrrra(m1,m2)+cuxxxxyy600*rxrrrrrra(m1,m2)+cuxxxxyy010*rxr(m1,m2,1)+cuxxxxyy110*rxrsa(m1,m2)+cuxxxxyy210*rxrrsa(m1,m2)+cuxxxxyy310*rxrrrsa(m1,m2)+cuxxxxyy410*rxrrrrsa(m1,m2)+cuxxxxyy510*rxrrrrrsa(m1,m2)+cuxxxxyy020*rxssa(m1,m2)+cuxxxxyy120*rxrssa(m1,m2)+cuxxxxyy220*rxrrssa(m1,m2)+cuxxxxyy320*rxrrrssa(m1,m2)+cuxxxxyy420*rxrrrrssa(m1,m2)+cuxxxxyy030*rxsssa(m1,m2)+cuxxxxyy130*rxrsssa(m1,m2)+cuxxxxyy230*rxrrsssa(m1,m2)+cuxxxxyy330*rxrrrsssa(m1,m2)+cuxxxxyy040*rxssssa(m1,m2)+cuxxxxyy140*rxrssssa(m1,m2)+cuxxxxyy240*rxrrssssa(m1,m2)+cuxxxxyy050*rxsssssa(m1,m2)+cuxxxxyy150*rxrsssssa(m1,m2)+cuxxxxyy060*rxssssssa(m1,m2)
                      !   rxxxxyyya(m1,m2) = cuxxxyyy100*rxr(m1,m2,0)+cuxxxyyy200*rxrra(m1,m2)+cuxxxyyy300*rxrrra(m1,m2)+cuxxxyyy400*rxrrrra(m1,m2)+cuxxxyyy500*rxrrrrra(m1,m2)+cuxxxyyy600*rxrrrrrra(m1,m2)+cuxxxyyy010*rxr(m1,m2,1)+cuxxxyyy110*rxrsa(m1,m2)+cuxxxyyy210*rxrrsa(m1,m2)+cuxxxyyy310*rxrrrsa(m1,m2)+cuxxxyyy410*rxrrrrsa(m1,m2)+cuxxxyyy510*rxrrrrrsa(m1,m2)+cuxxxyyy020*rxssa(m1,m2)+cuxxxyyy120*rxrssa(m1,m2)+cuxxxyyy220*rxrrssa(m1,m2)+cuxxxyyy320*rxrrrssa(m1,m2)+cuxxxyyy420*rxrrrrssa(m1,m2)+cuxxxyyy030*rxsssa(m1,m2)+cuxxxyyy130*rxrsssa(m1,m2)+cuxxxyyy230*rxrrsssa(m1,m2)+cuxxxyyy330*rxrrrsssa(m1,m2)+cuxxxyyy040*rxssssa(m1,m2)+cuxxxyyy140*rxrssssa(m1,m2)+cuxxxyyy240*rxrrssssa(m1,m2)+cuxxxyyy050*rxsssssa(m1,m2)+cuxxxyyy150*rxrsssssa(m1,m2)+cuxxxyyy060*rxssssssa(m1,m2)
                      !   rxxxyyyya(m1,m2) = cuxxyyyy100*rxr(m1,m2,0)+cuxxyyyy200*rxrra(m1,m2)+cuxxyyyy300*rxrrra(m1,m2)+cuxxyyyy400*rxrrrra(m1,m2)+cuxxyyyy500*rxrrrrra(m1,m2)+cuxxyyyy600*rxrrrrrra(m1,m2)+cuxxyyyy010*rxr(m1,m2,1)+cuxxyyyy110*rxrsa(m1,m2)+cuxxyyyy210*rxrrsa(m1,m2)+cuxxyyyy310*rxrrrsa(m1,m2)+cuxxyyyy410*rxrrrrsa(m1,m2)+cuxxyyyy510*rxrrrrrsa(m1,m2)+cuxxyyyy020*rxssa(m1,m2)+cuxxyyyy120*rxrssa(m1,m2)+cuxxyyyy220*rxrrssa(m1,m2)+cuxxyyyy320*rxrrrssa(m1,m2)+cuxxyyyy420*rxrrrrssa(m1,m2)+cuxxyyyy030*rxsssa(m1,m2)+cuxxyyyy130*rxrsssa(m1,m2)+cuxxyyyy230*rxrrsssa(m1,m2)+cuxxyyyy330*rxrrrsssa(m1,m2)+cuxxyyyy040*rxssssa(m1,m2)+cuxxyyyy140*rxrssssa(m1,m2)+cuxxyyyy240*rxrrssssa(m1,m2)+cuxxyyyy050*rxsssssa(m1,m2)+cuxxyyyy150*rxrsssssa(m1,m2)+cuxxyyyy060*rxssssssa(m1,m2)
                      !   rxxyyyyya(m1,m2) = cuxyyyyy100*rxr(m1,m2,0)+cuxyyyyy200*rxrra(m1,m2)+cuxyyyyy300*rxrrra(m1,m2)+cuxyyyyy400*rxrrrra(m1,m2)+cuxyyyyy500*rxrrrrra(m1,m2)+cuxyyyyy600*rxrrrrrra(m1,m2)+cuxyyyyy010*rxr(m1,m2,1)+cuxyyyyy110*rxrsa(m1,m2)+cuxyyyyy210*rxrrsa(m1,m2)+cuxyyyyy310*rxrrrsa(m1,m2)+cuxyyyyy410*rxrrrrsa(m1,m2)+cuxyyyyy510*rxrrrrrsa(m1,m2)+cuxyyyyy020*rxssa(m1,m2)+cuxyyyyy120*rxrssa(m1,m2)+cuxyyyyy220*rxrrssa(m1,m2)+cuxyyyyy320*rxrrrssa(m1,m2)+cuxyyyyy420*rxrrrrssa(m1,m2)+cuxyyyyy030*rxsssa(m1,m2)+cuxyyyyy130*rxrsssa(m1,m2)+cuxyyyyy230*rxrrsssa(m1,m2)+cuxyyyyy330*rxrrrsssa(m1,m2)+cuxyyyyy040*rxssssa(m1,m2)+cuxyyyyy140*rxrssssa(m1,m2)+cuxyyyyy240*rxrrssssa(m1,m2)+cuxyyyyy050*rxsssssa(m1,m2)+cuxyyyyy150*rxrsssssa(m1,m2)+cuxyyyyy060*rxssssssa(m1,m2)
                      !   rxyyyyyya(m1,m2) = cuyyyyyy100*rxr(m1,m2,0)+cuyyyyyy200*rxrra(m1,m2)+cuyyyyyy300*rxrrra(m1,m2)+cuyyyyyy400*rxrrrra(m1,m2)+cuyyyyyy500*rxrrrrra(m1,m2)+cuyyyyyy600*rxrrrrrra(m1,m2)+cuyyyyyy010*rxr(m1,m2,1)+cuyyyyyy110*rxrsa(m1,m2)+cuyyyyyy210*rxrrsa(m1,m2)+cuyyyyyy310*rxrrrsa(m1,m2)+cuyyyyyy410*rxrrrrsa(m1,m2)+cuyyyyyy510*rxrrrrrsa(m1,m2)+cuyyyyyy020*rxssa(m1,m2)+cuyyyyyy120*rxrssa(m1,m2)+cuyyyyyy220*rxrrssa(m1,m2)+cuyyyyyy320*rxrrrssa(m1,m2)+cuyyyyyy420*rxrrrrssa(m1,m2)+cuyyyyyy030*rxsssa(m1,m2)+cuyyyyyy130*rxrsssa(m1,m2)+cuyyyyyy230*rxrrsssa(m1,m2)+cuyyyyyy330*rxrrrsssa(m1,m2)+cuyyyyyy040*rxssssa(m1,m2)+cuyyyyyy140*rxrssssa(m1,m2)+cuyyyyyy240*rxrrssssa(m1,m2)+cuyyyyyy050*rxsssssa(m1,m2)+cuyyyyyy150*rxrsssssa(m1,m2)+cuyyyyyy060*rxssssssa(m1,m2)
                      ! end do
                      ! end do
                      ! rxxxxxxx=rxxxxxxxa(0,0); ryxxxxxx=rxxxxxxxa(0,1); sxxxxxxx=rxxxxxxxa(1,0); syxxxxxx=rxxxxxxxa(1,1); 
                      ! rxxxxxxy=rxxxxxxya(0,0); ryxxxxxy=rxxxxxxya(0,1); sxxxxxxy=rxxxxxxya(1,0); syxxxxxy=rxxxxxxya(1,1); 
                      ! rxxxxxyy=rxxxxxyya(0,0); ryxxxxyy=rxxxxxyya(0,1); sxxxxxyy=rxxxxxyya(1,0); syxxxxyy=rxxxxxyya(1,1); 
                      ! rxxxxyyy=rxxxxyyya(0,0); ryxxxyyy=rxxxxyyya(0,1); sxxxxyyy=rxxxxyyya(1,0); syxxxyyy=rxxxxyyya(1,1); 
                      ! rxxxyyyy=rxxxyyyya(0,0); ryxxyyyy=rxxxyyyya(0,1); sxxxyyyy=rxxxyyyya(1,0); syxxyyyy=rxxxyyyya(1,1); 
                      ! rxxyyyyy=rxxyyyyya(0,0); ryxyyyyy=rxxyyyyya(0,1); sxxyyyyy=rxxyyyyya(1,0); syxyyyyy=rxxyyyyya(1,1); 
                      ! rxyyyyyy=rxyyyyyya(0,0); ryyyyyyy=rxyyyyyya(0,1); sxyyyyyy=rxyyyyyya(1,0); syyyyyyy=rxyyyyyya(1,1); 
                      ! ! ---- seventh parametric derivatives ----
                      ! urrrrrrr = d700i
                      ! urrrrrrs = d610i
                      ! urrrrrss = d520i
                      ! urrrrsss = d430i
                      ! urrrssss = d340i
                      ! urrsssss = d250i
                      ! urssssss = d160i
                      ! usssssss = d070i 
                      ! ! ----- SEVENTH SPATIAL DERIVATIVES -----
                      ! getDerivCoeff2d(7)
                      ! ! uxxxxxxx = cuxxxxxxx100*ur+cuxxxxxxx200*urr+cuxxxxxxx300*urrr+cuxxxxxxx400*urrrr+cuxxxxxxx500*urrrrr+cuxxxxxxx600*urrrrrr+cuxxxxxxx700*urrrrrrr+cuxxxxxxx010*us+cuxxxxxxx110*urs+cuxxxxxxx210*urrs+cuxxxxxxx310*urrrs+cuxxxxxxx410*urrrrs+cuxxxxxxx510*urrrrrs+cuxxxxxxx610*urrrrrrs+cuxxxxxxx020*uss+cuxxxxxxx120*urss+cuxxxxxxx220*urrss+cuxxxxxxx320*urrrss+cuxxxxxxx420*urrrrss+cuxxxxxxx520*urrrrrss+cuxxxxxxx030*usss+cuxxxxxxx130*ursss+cuxxxxxxx230*urrsss+cuxxxxxxx330*urrrsss+cuxxxxxxx430*urrrrsss+cuxxxxxxx040*ussss+cuxxxxxxx140*urssss+cuxxxxxxx240*urrssss+cuxxxxxxx340*urrrssss+cuxxxxxxx050*usssss+cuxxxxxxx150*ursssss+cuxxxxxxx250*urrsssss+cuxxxxxxx060*ussssss+cuxxxxxxx160*urssssss+cuxxxxxxx070*usssssss
                      ! ! uxxxxxxy = cuxxxxxxy100*ur+cuxxxxxxy200*urr+cuxxxxxxy300*urrr+cuxxxxxxy400*urrrr+cuxxxxxxy500*urrrrr+cuxxxxxxy600*urrrrrr+cuxxxxxxy700*urrrrrrr+cuxxxxxxy010*us+cuxxxxxxy110*urs+cuxxxxxxy210*urrs+cuxxxxxxy310*urrrs+cuxxxxxxy410*urrrrs+cuxxxxxxy510*urrrrrs+cuxxxxxxy610*urrrrrrs+cuxxxxxxy020*uss+cuxxxxxxy120*urss+cuxxxxxxy220*urrss+cuxxxxxxy320*urrrss+cuxxxxxxy420*urrrrss+cuxxxxxxy520*urrrrrss+cuxxxxxxy030*usss+cuxxxxxxy130*ursss+cuxxxxxxy230*urrsss+cuxxxxxxy330*urrrsss+cuxxxxxxy430*urrrrsss+cuxxxxxxy040*ussss+cuxxxxxxy140*urssss+cuxxxxxxy240*urrssss+cuxxxxxxy340*urrrssss+cuxxxxxxy050*usssss+cuxxxxxxy150*ursssss+cuxxxxxxy250*urrsssss+cuxxxxxxy060*ussssss+cuxxxxxxy160*urssssss+cuxxxxxxy070*usssssss
                      ! ! uxxxxxyy = cuxxxxxyy100*ur+cuxxxxxyy200*urr+cuxxxxxyy300*urrr+cuxxxxxyy400*urrrr+cuxxxxxyy500*urrrrr+cuxxxxxyy600*urrrrrr+cuxxxxxyy700*urrrrrrr+cuxxxxxyy010*us+cuxxxxxyy110*urs+cuxxxxxyy210*urrs+cuxxxxxyy310*urrrs+cuxxxxxyy410*urrrrs+cuxxxxxyy510*urrrrrs+cuxxxxxyy610*urrrrrrs+cuxxxxxyy020*uss+cuxxxxxyy120*urss+cuxxxxxyy220*urrss+cuxxxxxyy320*urrrss+cuxxxxxyy420*urrrrss+cuxxxxxyy520*urrrrrss+cuxxxxxyy030*usss+cuxxxxxyy130*ursss+cuxxxxxyy230*urrsss+cuxxxxxyy330*urrrsss+cuxxxxxyy430*urrrrsss+cuxxxxxyy040*ussss+cuxxxxxyy140*urssss+cuxxxxxyy240*urrssss+cuxxxxxyy340*urrrssss+cuxxxxxyy050*usssss+cuxxxxxyy150*ursssss+cuxxxxxyy250*urrsssss+cuxxxxxyy060*ussssss+cuxxxxxyy160*urssssss+cuxxxxxyy070*usssssss
                      ! ! uxxxxyyy = cuxxxxyyy100*ur+cuxxxxyyy200*urr+cuxxxxyyy300*urrr+cuxxxxyyy400*urrrr+cuxxxxyyy500*urrrrr+cuxxxxyyy600*urrrrrr+cuxxxxyyy700*urrrrrrr+cuxxxxyyy010*us+cuxxxxyyy110*urs+cuxxxxyyy210*urrs+cuxxxxyyy310*urrrs+cuxxxxyyy410*urrrrs+cuxxxxyyy510*urrrrrs+cuxxxxyyy610*urrrrrrs+cuxxxxyyy020*uss+cuxxxxyyy120*urss+cuxxxxyyy220*urrss+cuxxxxyyy320*urrrss+cuxxxxyyy420*urrrrss+cuxxxxyyy520*urrrrrss+cuxxxxyyy030*usss+cuxxxxyyy130*ursss+cuxxxxyyy230*urrsss+cuxxxxyyy330*urrrsss+cuxxxxyyy430*urrrrsss+cuxxxxyyy040*ussss+cuxxxxyyy140*urssss+cuxxxxyyy240*urrssss+cuxxxxyyy340*urrrssss+cuxxxxyyy050*usssss+cuxxxxyyy150*ursssss+cuxxxxyyy250*urrsssss+cuxxxxyyy060*ussssss+cuxxxxyyy160*urssssss+cuxxxxyyy070*usssssss
                      ! ! uxxxyyyy = cuxxxyyyy100*ur+cuxxxyyyy200*urr+cuxxxyyyy300*urrr+cuxxxyyyy400*urrrr+cuxxxyyyy500*urrrrr+cuxxxyyyy600*urrrrrr+cuxxxyyyy700*urrrrrrr+cuxxxyyyy010*us+cuxxxyyyy110*urs+cuxxxyyyy210*urrs+cuxxxyyyy310*urrrs+cuxxxyyyy410*urrrrs+cuxxxyyyy510*urrrrrs+cuxxxyyyy610*urrrrrrs+cuxxxyyyy020*uss+cuxxxyyyy120*urss+cuxxxyyyy220*urrss+cuxxxyyyy320*urrrss+cuxxxyyyy420*urrrrss+cuxxxyyyy520*urrrrrss+cuxxxyyyy030*usss+cuxxxyyyy130*ursss+cuxxxyyyy230*urrsss+cuxxxyyyy330*urrrsss+cuxxxyyyy430*urrrrsss+cuxxxyyyy040*ussss+cuxxxyyyy140*urssss+cuxxxyyyy240*urrssss+cuxxxyyyy340*urrrssss+cuxxxyyyy050*usssss+cuxxxyyyy150*ursssss+cuxxxyyyy250*urrsssss+cuxxxyyyy060*ussssss+cuxxxyyyy160*urssssss+cuxxxyyyy070*usssssss
                      ! ! uxxyyyyy = cuxxyyyyy100*ur+cuxxyyyyy200*urr+cuxxyyyyy300*urrr+cuxxyyyyy400*urrrr+cuxxyyyyy500*urrrrr+cuxxyyyyy600*urrrrrr+cuxxyyyyy700*urrrrrrr+cuxxyyyyy010*us+cuxxyyyyy110*urs+cuxxyyyyy210*urrs+cuxxyyyyy310*urrrs+cuxxyyyyy410*urrrrs+cuxxyyyyy510*urrrrrs+cuxxyyyyy610*urrrrrrs+cuxxyyyyy020*uss+cuxxyyyyy120*urss+cuxxyyyyy220*urrss+cuxxyyyyy320*urrrss+cuxxyyyyy420*urrrrss+cuxxyyyyy520*urrrrrss+cuxxyyyyy030*usss+cuxxyyyyy130*ursss+cuxxyyyyy230*urrsss+cuxxyyyyy330*urrrsss+cuxxyyyyy430*urrrrsss+cuxxyyyyy040*ussss+cuxxyyyyy140*urssss+cuxxyyyyy240*urrssss+cuxxyyyyy340*urrrssss+cuxxyyyyy050*usssss+cuxxyyyyy150*ursssss+cuxxyyyyy250*urrsssss+cuxxyyyyy060*ussssss+cuxxyyyyy160*urssssss+cuxxyyyyy070*usssssss
                      ! ! uxyyyyyy = cuxyyyyyy100*ur+cuxyyyyyy200*urr+cuxyyyyyy300*urrr+cuxyyyyyy400*urrrr+cuxyyyyyy500*urrrrr+cuxyyyyyy600*urrrrrr+cuxyyyyyy700*urrrrrrr+cuxyyyyyy010*us+cuxyyyyyy110*urs+cuxyyyyyy210*urrs+cuxyyyyyy310*urrrs+cuxyyyyyy410*urrrrs+cuxyyyyyy510*urrrrrs+cuxyyyyyy610*urrrrrrs+cuxyyyyyy020*uss+cuxyyyyyy120*urss+cuxyyyyyy220*urrss+cuxyyyyyy320*urrrss+cuxyyyyyy420*urrrrss+cuxyyyyyy520*urrrrrss+cuxyyyyyy030*usss+cuxyyyyyy130*ursss+cuxyyyyyy230*urrsss+cuxyyyyyy330*urrrsss+cuxyyyyyy430*urrrrsss+cuxyyyyyy040*ussss+cuxyyyyyy140*urssss+cuxyyyyyy240*urrssss+cuxyyyyyy340*urrrssss+cuxyyyyyy050*usssss+cuxyyyyyy150*ursssss+cuxyyyyyy250*urrsssss+cuxyyyyyy060*ussssss+cuxyyyyyy160*urssssss+cuxyyyyyy070*usssssss
                      ! ! uyyyyyyy = cuyyyyyyy100*ur+cuyyyyyyy200*urr+cuyyyyyyy300*urrr+cuyyyyyyy400*urrrr+cuyyyyyyy500*urrrrr+cuyyyyyyy600*urrrrrr+cuyyyyyyy700*urrrrrrr+cuyyyyyyy010*us+cuyyyyyyy110*urs+cuyyyyyyy210*urrs+cuyyyyyyy310*urrrs+cuyyyyyyy410*urrrrs+cuyyyyyyy510*urrrrrs+cuyyyyyyy610*urrrrrrs+cuyyyyyyy020*uss+cuyyyyyyy120*urss+cuyyyyyyy220*urrss+cuyyyyyyy320*urrrss+cuyyyyyyy420*urrrrss+cuyyyyyyy520*urrrrrss+cuyyyyyyy030*usss+cuyyyyyyy130*ursss+cuyyyyyyy230*urrsss+cuyyyyyyy330*urrrsss+cuyyyyyyy430*urrrrsss+cuyyyyyyy040*ussss+cuyyyyyyy140*urssss+cuyyyyyyy240*urrssss+cuyyyyyyy340*urrrssss+cuyyyyyyy050*usssss+cuyyyyyyy150*ursssss+cuyyyyyyy250*urrsssss+cuyyyyyyy060*ussssss+cuyyyyyyy160*urssssss+cuyyyyyyy070*usssssss
                      ! ! ----- START EIGHTH DERIVATIVES -----
                      ! do m1=0,nd-1
                      ! do m2=0,nd-1
                      !   ! ---- 7th parameteric derivatives of the metrics ----
                      !   rxrrrrrrra(m1,m2) = rx700i(m1,m2)
                      !   rxrrrrrrsa(m1,m2) = rx610i(m1,m2)
                      !   rxrrrrrssa(m1,m2) = rx520i(m1,m2)
                      !   rxrrrrsssa(m1,m2) = rx430i(m1,m2)
                      !   rxrrrssssa(m1,m2) = rx340i(m1,m2)
                      !   rxrrsssssa(m1,m2) = rx250i(m1,m2)
                      !   rxrssssssa(m1,m2) = rx160i(m1,m2)
                      !   rxsssssssa(m1,m2) = rx070i(m1,m2) 
                      !   ! ---- seventh spatial derivatives of the metrics ----
                      !   rxxxxxxxxa(m1,m2) = cuxxxxxxx100*rxr(m1,m2,0)+cuxxxxxxx200*rxrra(m1,m2)+cuxxxxxxx300*rxrrra(m1,m2)+cuxxxxxxx400*rxrrrra(m1,m2)+cuxxxxxxx500*rxrrrrra(m1,m2)+cuxxxxxxx600*rxrrrrrra(m1,m2)+cuxxxxxxx700*rxrrrrrrra(m1,m2)+cuxxxxxxx010*rxr(m1,m2,1)+cuxxxxxxx110*rxrsa(m1,m2)+cuxxxxxxx210*rxrrsa(m1,m2)+cuxxxxxxx310*rxrrrsa(m1,m2)+cuxxxxxxx410*rxrrrrsa(m1,m2)+cuxxxxxxx510*rxrrrrrsa(m1,m2)+cuxxxxxxx610*rxrrrrrrsa(m1,m2)+cuxxxxxxx020*rxssa(m1,m2)+cuxxxxxxx120*rxrssa(m1,m2)+cuxxxxxxx220*rxrrssa(m1,m2)+cuxxxxxxx320*rxrrrssa(m1,m2)+cuxxxxxxx420*rxrrrrssa(m1,m2)+cuxxxxxxx520*rxrrrrrssa(m1,m2)+cuxxxxxxx030*rxsssa(m1,m2)+cuxxxxxxx130*rxrsssa(m1,m2)+cuxxxxxxx230*rxrrsssa(m1,m2)+cuxxxxxxx330*rxrrrsssa(m1,m2)+cuxxxxxxx430*rxrrrrsssa(m1,m2)+cuxxxxxxx040*rxssssa(m1,m2)+cuxxxxxxx140*rxrssssa(m1,m2)+cuxxxxxxx240*rxrrssssa(m1,m2)+cuxxxxxxx340*rxrrrssssa(m1,m2)+cuxxxxxxx050*rxsssssa(m1,m2)+cuxxxxxxx150*rxrsssssa(m1,m2)+cuxxxxxxx250*rxrrsssssa(m1,m2)+cuxxxxxxx060*rxssssssa(m1,m2)+cuxxxxxxx160*rxrssssssa(m1,m2)+cuxxxxxxx070*rxsssssssa(m1,m2)
                      !   rxxxxxxxya(m1,m2) = cuxxxxxxy100*rxr(m1,m2,0)+cuxxxxxxy200*rxrra(m1,m2)+cuxxxxxxy300*rxrrra(m1,m2)+cuxxxxxxy400*rxrrrra(m1,m2)+cuxxxxxxy500*rxrrrrra(m1,m2)+cuxxxxxxy600*rxrrrrrra(m1,m2)+cuxxxxxxy700*rxrrrrrrra(m1,m2)+cuxxxxxxy010*rxr(m1,m2,1)+cuxxxxxxy110*rxrsa(m1,m2)+cuxxxxxxy210*rxrrsa(m1,m2)+cuxxxxxxy310*rxrrrsa(m1,m2)+cuxxxxxxy410*rxrrrrsa(m1,m2)+cuxxxxxxy510*rxrrrrrsa(m1,m2)+cuxxxxxxy610*rxrrrrrrsa(m1,m2)+cuxxxxxxy020*rxssa(m1,m2)+cuxxxxxxy120*rxrssa(m1,m2)+cuxxxxxxy220*rxrrssa(m1,m2)+cuxxxxxxy320*rxrrrssa(m1,m2)+cuxxxxxxy420*rxrrrrssa(m1,m2)+cuxxxxxxy520*rxrrrrrssa(m1,m2)+cuxxxxxxy030*rxsssa(m1,m2)+cuxxxxxxy130*rxrsssa(m1,m2)+cuxxxxxxy230*rxrrsssa(m1,m2)+cuxxxxxxy330*rxrrrsssa(m1,m2)+cuxxxxxxy430*rxrrrrsssa(m1,m2)+cuxxxxxxy040*rxssssa(m1,m2)+cuxxxxxxy140*rxrssssa(m1,m2)+cuxxxxxxy240*rxrrssssa(m1,m2)+cuxxxxxxy340*rxrrrssssa(m1,m2)+cuxxxxxxy050*rxsssssa(m1,m2)+cuxxxxxxy150*rxrsssssa(m1,m2)+cuxxxxxxy250*rxrrsssssa(m1,m2)+cuxxxxxxy060*rxssssssa(m1,m2)+cuxxxxxxy160*rxrssssssa(m1,m2)+cuxxxxxxy070*rxsssssssa(m1,m2)
                      !   rxxxxxxyya(m1,m2) = cuxxxxxyy100*rxr(m1,m2,0)+cuxxxxxyy200*rxrra(m1,m2)+cuxxxxxyy300*rxrrra(m1,m2)+cuxxxxxyy400*rxrrrra(m1,m2)+cuxxxxxyy500*rxrrrrra(m1,m2)+cuxxxxxyy600*rxrrrrrra(m1,m2)+cuxxxxxyy700*rxrrrrrrra(m1,m2)+cuxxxxxyy010*rxr(m1,m2,1)+cuxxxxxyy110*rxrsa(m1,m2)+cuxxxxxyy210*rxrrsa(m1,m2)+cuxxxxxyy310*rxrrrsa(m1,m2)+cuxxxxxyy410*rxrrrrsa(m1,m2)+cuxxxxxyy510*rxrrrrrsa(m1,m2)+cuxxxxxyy610*rxrrrrrrsa(m1,m2)+cuxxxxxyy020*rxssa(m1,m2)+cuxxxxxyy120*rxrssa(m1,m2)+cuxxxxxyy220*rxrrssa(m1,m2)+cuxxxxxyy320*rxrrrssa(m1,m2)+cuxxxxxyy420*rxrrrrssa(m1,m2)+cuxxxxxyy520*rxrrrrrssa(m1,m2)+cuxxxxxyy030*rxsssa(m1,m2)+cuxxxxxyy130*rxrsssa(m1,m2)+cuxxxxxyy230*rxrrsssa(m1,m2)+cuxxxxxyy330*rxrrrsssa(m1,m2)+cuxxxxxyy430*rxrrrrsssa(m1,m2)+cuxxxxxyy040*rxssssa(m1,m2)+cuxxxxxyy140*rxrssssa(m1,m2)+cuxxxxxyy240*rxrrssssa(m1,m2)+cuxxxxxyy340*rxrrrssssa(m1,m2)+cuxxxxxyy050*rxsssssa(m1,m2)+cuxxxxxyy150*rxrsssssa(m1,m2)+cuxxxxxyy250*rxrrsssssa(m1,m2)+cuxxxxxyy060*rxssssssa(m1,m2)+cuxxxxxyy160*rxrssssssa(m1,m2)+cuxxxxxyy070*rxsssssssa(m1,m2)
                      !   rxxxxxyyya(m1,m2) = cuxxxxyyy100*rxr(m1,m2,0)+cuxxxxyyy200*rxrra(m1,m2)+cuxxxxyyy300*rxrrra(m1,m2)+cuxxxxyyy400*rxrrrra(m1,m2)+cuxxxxyyy500*rxrrrrra(m1,m2)+cuxxxxyyy600*rxrrrrrra(m1,m2)+cuxxxxyyy700*rxrrrrrrra(m1,m2)+cuxxxxyyy010*rxr(m1,m2,1)+cuxxxxyyy110*rxrsa(m1,m2)+cuxxxxyyy210*rxrrsa(m1,m2)+cuxxxxyyy310*rxrrrsa(m1,m2)+cuxxxxyyy410*rxrrrrsa(m1,m2)+cuxxxxyyy510*rxrrrrrsa(m1,m2)+cuxxxxyyy610*rxrrrrrrsa(m1,m2)+cuxxxxyyy020*rxssa(m1,m2)+cuxxxxyyy120*rxrssa(m1,m2)+cuxxxxyyy220*rxrrssa(m1,m2)+cuxxxxyyy320*rxrrrssa(m1,m2)+cuxxxxyyy420*rxrrrrssa(m1,m2)+cuxxxxyyy520*rxrrrrrssa(m1,m2)+cuxxxxyyy030*rxsssa(m1,m2)+cuxxxxyyy130*rxrsssa(m1,m2)+cuxxxxyyy230*rxrrsssa(m1,m2)+cuxxxxyyy330*rxrrrsssa(m1,m2)+cuxxxxyyy430*rxrrrrsssa(m1,m2)+cuxxxxyyy040*rxssssa(m1,m2)+cuxxxxyyy140*rxrssssa(m1,m2)+cuxxxxyyy240*rxrrssssa(m1,m2)+cuxxxxyyy340*rxrrrssssa(m1,m2)+cuxxxxyyy050*rxsssssa(m1,m2)+cuxxxxyyy150*rxrsssssa(m1,m2)+cuxxxxyyy250*rxrrsssssa(m1,m2)+cuxxxxyyy060*rxssssssa(m1,m2)+cuxxxxyyy160*rxrssssssa(m1,m2)+cuxxxxyyy070*rxsssssssa(m1,m2)
                      !   rxxxxyyyya(m1,m2) = cuxxxyyyy100*rxr(m1,m2,0)+cuxxxyyyy200*rxrra(m1,m2)+cuxxxyyyy300*rxrrra(m1,m2)+cuxxxyyyy400*rxrrrra(m1,m2)+cuxxxyyyy500*rxrrrrra(m1,m2)+cuxxxyyyy600*rxrrrrrra(m1,m2)+cuxxxyyyy700*rxrrrrrrra(m1,m2)+cuxxxyyyy010*rxr(m1,m2,1)+cuxxxyyyy110*rxrsa(m1,m2)+cuxxxyyyy210*rxrrsa(m1,m2)+cuxxxyyyy310*rxrrrsa(m1,m2)+cuxxxyyyy410*rxrrrrsa(m1,m2)+cuxxxyyyy510*rxrrrrrsa(m1,m2)+cuxxxyyyy610*rxrrrrrrsa(m1,m2)+cuxxxyyyy020*rxssa(m1,m2)+cuxxxyyyy120*rxrssa(m1,m2)+cuxxxyyyy220*rxrrssa(m1,m2)+cuxxxyyyy320*rxrrrssa(m1,m2)+cuxxxyyyy420*rxrrrrssa(m1,m2)+cuxxxyyyy520*rxrrrrrssa(m1,m2)+cuxxxyyyy030*rxsssa(m1,m2)+cuxxxyyyy130*rxrsssa(m1,m2)+cuxxxyyyy230*rxrrsssa(m1,m2)+cuxxxyyyy330*rxrrrsssa(m1,m2)+cuxxxyyyy430*rxrrrrsssa(m1,m2)+cuxxxyyyy040*rxssssa(m1,m2)+cuxxxyyyy140*rxrssssa(m1,m2)+cuxxxyyyy240*rxrrssssa(m1,m2)+cuxxxyyyy340*rxrrrssssa(m1,m2)+cuxxxyyyy050*rxsssssa(m1,m2)+cuxxxyyyy150*rxrsssssa(m1,m2)+cuxxxyyyy250*rxrrsssssa(m1,m2)+cuxxxyyyy060*rxssssssa(m1,m2)+cuxxxyyyy160*rxrssssssa(m1,m2)+cuxxxyyyy070*rxsssssssa(m1,m2)
                      !   rxxxyyyyya(m1,m2) = cuxxyyyyy100*rxr(m1,m2,0)+cuxxyyyyy200*rxrra(m1,m2)+cuxxyyyyy300*rxrrra(m1,m2)+cuxxyyyyy400*rxrrrra(m1,m2)+cuxxyyyyy500*rxrrrrra(m1,m2)+cuxxyyyyy600*rxrrrrrra(m1,m2)+cuxxyyyyy700*rxrrrrrrra(m1,m2)+cuxxyyyyy010*rxr(m1,m2,1)+cuxxyyyyy110*rxrsa(m1,m2)+cuxxyyyyy210*rxrrsa(m1,m2)+cuxxyyyyy310*rxrrrsa(m1,m2)+cuxxyyyyy410*rxrrrrsa(m1,m2)+cuxxyyyyy510*rxrrrrrsa(m1,m2)+cuxxyyyyy610*rxrrrrrrsa(m1,m2)+cuxxyyyyy020*rxssa(m1,m2)+cuxxyyyyy120*rxrssa(m1,m2)+cuxxyyyyy220*rxrrssa(m1,m2)+cuxxyyyyy320*rxrrrssa(m1,m2)+cuxxyyyyy420*rxrrrrssa(m1,m2)+cuxxyyyyy520*rxrrrrrssa(m1,m2)+cuxxyyyyy030*rxsssa(m1,m2)+cuxxyyyyy130*rxrsssa(m1,m2)+cuxxyyyyy230*rxrrsssa(m1,m2)+cuxxyyyyy330*rxrrrsssa(m1,m2)+cuxxyyyyy430*rxrrrrsssa(m1,m2)+cuxxyyyyy040*rxssssa(m1,m2)+cuxxyyyyy140*rxrssssa(m1,m2)+cuxxyyyyy240*rxrrssssa(m1,m2)+cuxxyyyyy340*rxrrrssssa(m1,m2)+cuxxyyyyy050*rxsssssa(m1,m2)+cuxxyyyyy150*rxrsssssa(m1,m2)+cuxxyyyyy250*rxrrsssssa(m1,m2)+cuxxyyyyy060*rxssssssa(m1,m2)+cuxxyyyyy160*rxrssssssa(m1,m2)+cuxxyyyyy070*rxsssssssa(m1,m2)
                      !   rxxyyyyyya(m1,m2) = cuxyyyyyy100*rxr(m1,m2,0)+cuxyyyyyy200*rxrra(m1,m2)+cuxyyyyyy300*rxrrra(m1,m2)+cuxyyyyyy400*rxrrrra(m1,m2)+cuxyyyyyy500*rxrrrrra(m1,m2)+cuxyyyyyy600*rxrrrrrra(m1,m2)+cuxyyyyyy700*rxrrrrrrra(m1,m2)+cuxyyyyyy010*rxr(m1,m2,1)+cuxyyyyyy110*rxrsa(m1,m2)+cuxyyyyyy210*rxrrsa(m1,m2)+cuxyyyyyy310*rxrrrsa(m1,m2)+cuxyyyyyy410*rxrrrrsa(m1,m2)+cuxyyyyyy510*rxrrrrrsa(m1,m2)+cuxyyyyyy610*rxrrrrrrsa(m1,m2)+cuxyyyyyy020*rxssa(m1,m2)+cuxyyyyyy120*rxrssa(m1,m2)+cuxyyyyyy220*rxrrssa(m1,m2)+cuxyyyyyy320*rxrrrssa(m1,m2)+cuxyyyyyy420*rxrrrrssa(m1,m2)+cuxyyyyyy520*rxrrrrrssa(m1,m2)+cuxyyyyyy030*rxsssa(m1,m2)+cuxyyyyyy130*rxrsssa(m1,m2)+cuxyyyyyy230*rxrrsssa(m1,m2)+cuxyyyyyy330*rxrrrsssa(m1,m2)+cuxyyyyyy430*rxrrrrsssa(m1,m2)+cuxyyyyyy040*rxssssa(m1,m2)+cuxyyyyyy140*rxrssssa(m1,m2)+cuxyyyyyy240*rxrrssssa(m1,m2)+cuxyyyyyy340*rxrrrssssa(m1,m2)+cuxyyyyyy050*rxsssssa(m1,m2)+cuxyyyyyy150*rxrsssssa(m1,m2)+cuxyyyyyy250*rxrrsssssa(m1,m2)+cuxyyyyyy060*rxssssssa(m1,m2)+cuxyyyyyy160*rxrssssssa(m1,m2)+cuxyyyyyy070*rxsssssssa(m1,m2)
                      !   rxyyyyyyya(m1,m2) = cuyyyyyyy100*rxr(m1,m2,0)+cuyyyyyyy200*rxrra(m1,m2)+cuyyyyyyy300*rxrrra(m1,m2)+cuyyyyyyy400*rxrrrra(m1,m2)+cuyyyyyyy500*rxrrrrra(m1,m2)+cuyyyyyyy600*rxrrrrrra(m1,m2)+cuyyyyyyy700*rxrrrrrrra(m1,m2)+cuyyyyyyy010*rxr(m1,m2,1)+cuyyyyyyy110*rxrsa(m1,m2)+cuyyyyyyy210*rxrrsa(m1,m2)+cuyyyyyyy310*rxrrrsa(m1,m2)+cuyyyyyyy410*rxrrrrsa(m1,m2)+cuyyyyyyy510*rxrrrrrsa(m1,m2)+cuyyyyyyy610*rxrrrrrrsa(m1,m2)+cuyyyyyyy020*rxssa(m1,m2)+cuyyyyyyy120*rxrssa(m1,m2)+cuyyyyyyy220*rxrrssa(m1,m2)+cuyyyyyyy320*rxrrrssa(m1,m2)+cuyyyyyyy420*rxrrrrssa(m1,m2)+cuyyyyyyy520*rxrrrrrssa(m1,m2)+cuyyyyyyy030*rxsssa(m1,m2)+cuyyyyyyy130*rxrsssa(m1,m2)+cuyyyyyyy230*rxrrsssa(m1,m2)+cuyyyyyyy330*rxrrrsssa(m1,m2)+cuyyyyyyy430*rxrrrrsssa(m1,m2)+cuyyyyyyy040*rxssssa(m1,m2)+cuyyyyyyy140*rxrssssa(m1,m2)+cuyyyyyyy240*rxrrssssa(m1,m2)+cuyyyyyyy340*rxrrrssssa(m1,m2)+cuyyyyyyy050*rxsssssa(m1,m2)+cuyyyyyyy150*rxrsssssa(m1,m2)+cuyyyyyyy250*rxrrsssssa(m1,m2)+cuyyyyyyy060*rxssssssa(m1,m2)+cuyyyyyyy160*rxrssssssa(m1,m2)+cuyyyyyyy070*rxsssssssa(m1,m2)
                      ! end do
                      ! end do
                      ! rxxxxxxxx=rxxxxxxxxa(0,0); ryxxxxxxx=rxxxxxxxxa(0,1); sxxxxxxxx=rxxxxxxxxa(1,0); syxxxxxxx=rxxxxxxxxa(1,1); 
                      ! rxxxxxxxy=rxxxxxxxya(0,0); ryxxxxxxy=rxxxxxxxya(0,1); sxxxxxxxy=rxxxxxxxya(1,0); syxxxxxxy=rxxxxxxxya(1,1); 
                      ! rxxxxxxyy=rxxxxxxyya(0,0); ryxxxxxyy=rxxxxxxyya(0,1); sxxxxxxyy=rxxxxxxyya(1,0); syxxxxxyy=rxxxxxxyya(1,1); 
                      ! rxxxxxyyy=rxxxxxyyya(0,0); ryxxxxyyy=rxxxxxyyya(0,1); sxxxxxyyy=rxxxxxyyya(1,0); syxxxxyyy=rxxxxxyyya(1,1); 
                      ! rxxxxyyyy=rxxxxyyyya(0,0); ryxxxyyyy=rxxxxyyyya(0,1); sxxxxyyyy=rxxxxyyyya(1,0); syxxxyyyy=rxxxxyyyya(1,1); 
                      ! rxxxyyyyy=rxxxyyyyya(0,0); ryxxyyyyy=rxxxyyyyya(0,1); sxxxyyyyy=rxxxyyyyya(1,0); syxxyyyyy=rxxxyyyyya(1,1); 
                      ! rxxyyyyyy=rxxyyyyyya(0,0); ryxyyyyyy=rxxyyyyyya(0,1); sxxyyyyyy=rxxyyyyyya(1,0); syxyyyyyy=rxxyyyyyya(1,1); 
                      ! rxyyyyyyy=rxyyyyyyya(0,0); ryyyyyyyy=rxyyyyyyya(0,1); sxyyyyyyy=rxyyyyyyya(1,0); syyyyyyyy=rxyyyyyyya(1,1); 
                      ! ! ---- eighth parametric derivatives ----
                      ! urrrrrrrr = d800i
                      ! urrrrrrrs = d710i
                      ! urrrrrrss = d620i
                      ! urrrrrsss = d530i
                      ! urrrrssss = d440i
                      ! urrrsssss = d350i
                      ! urrssssss = d260i
                      ! ursssssss = d170i
                      ! ussssssss = d080i
                      ! ! ----- EIGHTH SPATIAL DERIVATIVES -----
                      ! getDerivCoeff2d(8)  
                      ! uxxxxxxxx = cuxxxxxxxx100*ur+cuxxxxxxxx200*urr+cuxxxxxxxx300*urrr+cuxxxxxxxx400*urrrr+cuxxxxxxxx500*urrrrr+cuxxxxxxxx600*urrrrrr+cuxxxxxxxx700*urrrrrrr+cuxxxxxxxx800*urrrrrrrr+cuxxxxxxxx010*us+cuxxxxxxxx110*urs+cuxxxxxxxx210*urrs+cuxxxxxxxx310*urrrs+cuxxxxxxxx410*urrrrs+cuxxxxxxxx510*urrrrrs+cuxxxxxxxx610*urrrrrrs+cuxxxxxxxx710*urrrrrrrs+cuxxxxxxxx020*uss+cuxxxxxxxx120*urss+cuxxxxxxxx220*urrss+cuxxxxxxxx320*urrrss+cuxxxxxxxx420*urrrrss+cuxxxxxxxx520*urrrrrss+cuxxxxxxxx620*urrrrrrss+cuxxxxxxxx030*usss+cuxxxxxxxx130*ursss+cuxxxxxxxx230*urrsss+cuxxxxxxxx330*urrrsss+cuxxxxxxxx430*urrrrsss+cuxxxxxxxx530*urrrrrsss+cuxxxxxxxx040*ussss+cuxxxxxxxx140*urssss+cuxxxxxxxx240*urrssss+cuxxxxxxxx340*urrrssss+cuxxxxxxxx440*urrrrssss+cuxxxxxxxx050*usssss+cuxxxxxxxx150*ursssss+cuxxxxxxxx250*urrsssss+cuxxxxxxxx350*urrrsssss+cuxxxxxxxx060*ussssss+cuxxxxxxxx160*urssssss+cuxxxxxxxx260*urrssssss+cuxxxxxxxx070*usssssss+cuxxxxxxxx170*ursssssss+cuxxxxxxxx080*ussssssss
                      ! ! uxxxxxxxy = cuxxxxxxxy100*ur+cuxxxxxxxy200*urr+cuxxxxxxxy300*urrr+cuxxxxxxxy400*urrrr+cuxxxxxxxy500*urrrrr+cuxxxxxxxy600*urrrrrr+cuxxxxxxxy700*urrrrrrr+cuxxxxxxxy800*urrrrrrrr+cuxxxxxxxy010*us+cuxxxxxxxy110*urs+cuxxxxxxxy210*urrs+cuxxxxxxxy310*urrrs+cuxxxxxxxy410*urrrrs+cuxxxxxxxy510*urrrrrs+cuxxxxxxxy610*urrrrrrs+cuxxxxxxxy710*urrrrrrrs+cuxxxxxxxy020*uss+cuxxxxxxxy120*urss+cuxxxxxxxy220*urrss+cuxxxxxxxy320*urrrss+cuxxxxxxxy420*urrrrss+cuxxxxxxxy520*urrrrrss+cuxxxxxxxy620*urrrrrrss+cuxxxxxxxy030*usss+cuxxxxxxxy130*ursss+cuxxxxxxxy230*urrsss+cuxxxxxxxy330*urrrsss+cuxxxxxxxy430*urrrrsss+cuxxxxxxxy530*urrrrrsss+cuxxxxxxxy040*ussss+cuxxxxxxxy140*urssss+cuxxxxxxxy240*urrssss+cuxxxxxxxy340*urrrssss+cuxxxxxxxy440*urrrrssss+cuxxxxxxxy050*usssss+cuxxxxxxxy150*ursssss+cuxxxxxxxy250*urrsssss+cuxxxxxxxy350*urrrsssss+cuxxxxxxxy060*ussssss+cuxxxxxxxy160*urssssss+cuxxxxxxxy260*urrssssss+cuxxxxxxxy070*usssssss+cuxxxxxxxy170*ursssssss+cuxxxxxxxy080*ussssssss
                      ! uxxxxxxyy = cuxxxxxxyy100*ur+cuxxxxxxyy200*urr+cuxxxxxxyy300*urrr+cuxxxxxxyy400*urrrr+cuxxxxxxyy500*urrrrr+cuxxxxxxyy600*urrrrrr+cuxxxxxxyy700*urrrrrrr+cuxxxxxxyy800*urrrrrrrr+cuxxxxxxyy010*us+cuxxxxxxyy110*urs+cuxxxxxxyy210*urrs+cuxxxxxxyy310*urrrs+cuxxxxxxyy410*urrrrs+cuxxxxxxyy510*urrrrrs+cuxxxxxxyy610*urrrrrrs+cuxxxxxxyy710*urrrrrrrs+cuxxxxxxyy020*uss+cuxxxxxxyy120*urss+cuxxxxxxyy220*urrss+cuxxxxxxyy320*urrrss+cuxxxxxxyy420*urrrrss+cuxxxxxxyy520*urrrrrss+cuxxxxxxyy620*urrrrrrss+cuxxxxxxyy030*usss+cuxxxxxxyy130*ursss+cuxxxxxxyy230*urrsss+cuxxxxxxyy330*urrrsss+cuxxxxxxyy430*urrrrsss+cuxxxxxxyy530*urrrrrsss+cuxxxxxxyy040*ussss+cuxxxxxxyy140*urssss+cuxxxxxxyy240*urrssss+cuxxxxxxyy340*urrrssss+cuxxxxxxyy440*urrrrssss+cuxxxxxxyy050*usssss+cuxxxxxxyy150*ursssss+cuxxxxxxyy250*urrsssss+cuxxxxxxyy350*urrrsssss+cuxxxxxxyy060*ussssss+cuxxxxxxyy160*urssssss+cuxxxxxxyy260*urrssssss+cuxxxxxxyy070*usssssss+cuxxxxxxyy170*ursssssss+cuxxxxxxyy080*ussssssss
                      ! ! uxxxxxyyy = cuxxxxxyyy100*ur+cuxxxxxyyy200*urr+cuxxxxxyyy300*urrr+cuxxxxxyyy400*urrrr+cuxxxxxyyy500*urrrrr+cuxxxxxyyy600*urrrrrr+cuxxxxxyyy700*urrrrrrr+cuxxxxxyyy800*urrrrrrrr+cuxxxxxyyy010*us+cuxxxxxyyy110*urs+cuxxxxxyyy210*urrs+cuxxxxxyyy310*urrrs+cuxxxxxyyy410*urrrrs+cuxxxxxyyy510*urrrrrs+cuxxxxxyyy610*urrrrrrs+cuxxxxxyyy710*urrrrrrrs+cuxxxxxyyy020*uss+cuxxxxxyyy120*urss+cuxxxxxyyy220*urrss+cuxxxxxyyy320*urrrss+cuxxxxxyyy420*urrrrss+cuxxxxxyyy520*urrrrrss+cuxxxxxyyy620*urrrrrrss+cuxxxxxyyy030*usss+cuxxxxxyyy130*ursss+cuxxxxxyyy230*urrsss+cuxxxxxyyy330*urrrsss+cuxxxxxyyy430*urrrrsss+cuxxxxxyyy530*urrrrrsss+cuxxxxxyyy040*ussss+cuxxxxxyyy140*urssss+cuxxxxxyyy240*urrssss+cuxxxxxyyy340*urrrssss+cuxxxxxyyy440*urrrrssss+cuxxxxxyyy050*usssss+cuxxxxxyyy150*ursssss+cuxxxxxyyy250*urrsssss+cuxxxxxyyy350*urrrsssss+cuxxxxxyyy060*ussssss+cuxxxxxyyy160*urssssss+cuxxxxxyyy260*urrssssss+cuxxxxxyyy070*usssssss+cuxxxxxyyy170*ursssssss+cuxxxxxyyy080*ussssssss
                      ! uxxxxyyyy = cuxxxxyyyy100*ur+cuxxxxyyyy200*urr+cuxxxxyyyy300*urrr+cuxxxxyyyy400*urrrr+cuxxxxyyyy500*urrrrr+cuxxxxyyyy600*urrrrrr+cuxxxxyyyy700*urrrrrrr+cuxxxxyyyy800*urrrrrrrr+cuxxxxyyyy010*us+cuxxxxyyyy110*urs+cuxxxxyyyy210*urrs+cuxxxxyyyy310*urrrs+cuxxxxyyyy410*urrrrs+cuxxxxyyyy510*urrrrrs+cuxxxxyyyy610*urrrrrrs+cuxxxxyyyy710*urrrrrrrs+cuxxxxyyyy020*uss+cuxxxxyyyy120*urss+cuxxxxyyyy220*urrss+cuxxxxyyyy320*urrrss+cuxxxxyyyy420*urrrrss+cuxxxxyyyy520*urrrrrss+cuxxxxyyyy620*urrrrrrss+cuxxxxyyyy030*usss+cuxxxxyyyy130*ursss+cuxxxxyyyy230*urrsss+cuxxxxyyyy330*urrrsss+cuxxxxyyyy430*urrrrsss+cuxxxxyyyy530*urrrrrsss+cuxxxxyyyy040*ussss+cuxxxxyyyy140*urssss+cuxxxxyyyy240*urrssss+cuxxxxyyyy340*urrrssss+cuxxxxyyyy440*urrrrssss+cuxxxxyyyy050*usssss+cuxxxxyyyy150*ursssss+cuxxxxyyyy250*urrsssss+cuxxxxyyyy350*urrrsssss+cuxxxxyyyy060*ussssss+cuxxxxyyyy160*urssssss+cuxxxxyyyy260*urrssssss+cuxxxxyyyy070*usssssss+cuxxxxyyyy170*ursssssss+cuxxxxyyyy080*ussssssss
                      ! ! uxxxyyyyy = cuxxxyyyyy100*ur+cuxxxyyyyy200*urr+cuxxxyyyyy300*urrr+cuxxxyyyyy400*urrrr+cuxxxyyyyy500*urrrrr+cuxxxyyyyy600*urrrrrr+cuxxxyyyyy700*urrrrrrr+cuxxxyyyyy800*urrrrrrrr+cuxxxyyyyy010*us+cuxxxyyyyy110*urs+cuxxxyyyyy210*urrs+cuxxxyyyyy310*urrrs+cuxxxyyyyy410*urrrrs+cuxxxyyyyy510*urrrrrs+cuxxxyyyyy610*urrrrrrs+cuxxxyyyyy710*urrrrrrrs+cuxxxyyyyy020*uss+cuxxxyyyyy120*urss+cuxxxyyyyy220*urrss+cuxxxyyyyy320*urrrss+cuxxxyyyyy420*urrrrss+cuxxxyyyyy520*urrrrrss+cuxxxyyyyy620*urrrrrrss+cuxxxyyyyy030*usss+cuxxxyyyyy130*ursss+cuxxxyyyyy230*urrsss+cuxxxyyyyy330*urrrsss+cuxxxyyyyy430*urrrrsss+cuxxxyyyyy530*urrrrrsss+cuxxxyyyyy040*ussss+cuxxxyyyyy140*urssss+cuxxxyyyyy240*urrssss+cuxxxyyyyy340*urrrssss+cuxxxyyyyy440*urrrrssss+cuxxxyyyyy050*usssss+cuxxxyyyyy150*ursssss+cuxxxyyyyy250*urrsssss+cuxxxyyyyy350*urrrsssss+cuxxxyyyyy060*ussssss+cuxxxyyyyy160*urssssss+cuxxxyyyyy260*urrssssss+cuxxxyyyyy070*usssssss+cuxxxyyyyy170*ursssssss+cuxxxyyyyy080*ussssssss
                      ! uxxyyyyyy = cuxxyyyyyy100*ur+cuxxyyyyyy200*urr+cuxxyyyyyy300*urrr+cuxxyyyyyy400*urrrr+cuxxyyyyyy500*urrrrr+cuxxyyyyyy600*urrrrrr+cuxxyyyyyy700*urrrrrrr+cuxxyyyyyy800*urrrrrrrr+cuxxyyyyyy010*us+cuxxyyyyyy110*urs+cuxxyyyyyy210*urrs+cuxxyyyyyy310*urrrs+cuxxyyyyyy410*urrrrs+cuxxyyyyyy510*urrrrrs+cuxxyyyyyy610*urrrrrrs+cuxxyyyyyy710*urrrrrrrs+cuxxyyyyyy020*uss+cuxxyyyyyy120*urss+cuxxyyyyyy220*urrss+cuxxyyyyyy320*urrrss+cuxxyyyyyy420*urrrrss+cuxxyyyyyy520*urrrrrss+cuxxyyyyyy620*urrrrrrss+cuxxyyyyyy030*usss+cuxxyyyyyy130*ursss+cuxxyyyyyy230*urrsss+cuxxyyyyyy330*urrrsss+cuxxyyyyyy430*urrrrsss+cuxxyyyyyy530*urrrrrsss+cuxxyyyyyy040*ussss+cuxxyyyyyy140*urssss+cuxxyyyyyy240*urrssss+cuxxyyyyyy340*urrrssss+cuxxyyyyyy440*urrrrssss+cuxxyyyyyy050*usssss+cuxxyyyyyy150*ursssss+cuxxyyyyyy250*urrsssss+cuxxyyyyyy350*urrrsssss+cuxxyyyyyy060*ussssss+cuxxyyyyyy160*urssssss+cuxxyyyyyy260*urrssssss+cuxxyyyyyy070*usssssss+cuxxyyyyyy170*ursssssss+cuxxyyyyyy080*ussssssss
                      ! ! uxyyyyyyy = cuxyyyyyyy100*ur+cuxyyyyyyy200*urr+cuxyyyyyyy300*urrr+cuxyyyyyyy400*urrrr+cuxyyyyyyy500*urrrrr+cuxyyyyyyy600*urrrrrr+cuxyyyyyyy700*urrrrrrr+cuxyyyyyyy800*urrrrrrrr+cuxyyyyyyy010*us+cuxyyyyyyy110*urs+cuxyyyyyyy210*urrs+cuxyyyyyyy310*urrrs+cuxyyyyyyy410*urrrrs+cuxyyyyyyy510*urrrrrs+cuxyyyyyyy610*urrrrrrs+cuxyyyyyyy710*urrrrrrrs+cuxyyyyyyy020*uss+cuxyyyyyyy120*urss+cuxyyyyyyy220*urrss+cuxyyyyyyy320*urrrss+cuxyyyyyyy420*urrrrss+cuxyyyyyyy520*urrrrrss+cuxyyyyyyy620*urrrrrrss+cuxyyyyyyy030*usss+cuxyyyyyyy130*ursss+cuxyyyyyyy230*urrsss+cuxyyyyyyy330*urrrsss+cuxyyyyyyy430*urrrrsss+cuxyyyyyyy530*urrrrrsss+cuxyyyyyyy040*ussss+cuxyyyyyyy140*urssss+cuxyyyyyyy240*urrssss+cuxyyyyyyy340*urrrssss+cuxyyyyyyy440*urrrrssss+cuxyyyyyyy050*usssss+cuxyyyyyyy150*ursssss+cuxyyyyyyy250*urrsssss+cuxyyyyyyy350*urrrsssss+cuxyyyyyyy060*ussssss+cuxyyyyyyy160*urssssss+cuxyyyyyyy260*urrssssss+cuxyyyyyyy070*usssssss+cuxyyyyyyy170*ursssssss+cuxyyyyyyy080*ussssssss
                      ! uyyyyyyyy = cuyyyyyyyy100*ur+cuyyyyyyyy200*urr+cuyyyyyyyy300*urrr+cuyyyyyyyy400*urrrr+cuyyyyyyyy500*urrrrr+cuyyyyyyyy600*urrrrrr+cuyyyyyyyy700*urrrrrrr+cuyyyyyyyy800*urrrrrrrr+cuyyyyyyyy010*us+cuyyyyyyyy110*urs+cuyyyyyyyy210*urrs+cuyyyyyyyy310*urrrs+cuyyyyyyyy410*urrrrs+cuyyyyyyyy510*urrrrrs+cuyyyyyyyy610*urrrrrrs+cuyyyyyyyy710*urrrrrrrs+cuyyyyyyyy020*uss+cuyyyyyyyy120*urss+cuyyyyyyyy220*urrss+cuyyyyyyyy320*urrrss+cuyyyyyyyy420*urrrrss+cuyyyyyyyy520*urrrrrss+cuyyyyyyyy620*urrrrrrss+cuyyyyyyyy030*usss+cuyyyyyyyy130*ursss+cuyyyyyyyy230*urrsss+cuyyyyyyyy330*urrrsss+cuyyyyyyyy430*urrrrsss+cuyyyyyyyy530*urrrrrsss+cuyyyyyyyy040*ussss+cuyyyyyyyy140*urssss+cuyyyyyyyy240*urrssss+cuyyyyyyyy340*urrrssss+cuyyyyyyyy440*urrrrssss+cuyyyyyyyy050*usssss+cuyyyyyyyy150*ursssss+cuyyyyyyyy250*urrsssss+cuyyyyyyyy350*urrrsssss+cuyyyyyyyy060*ussssss+cuyyyyyyyy160*urssssss+cuyyyyyyyy260*urrssssss+cuyyyyyyyy070*usssssss+cuyyyyyyyy170*ursssssss+cuyyyyyyyy080*ussssssss
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
                                            un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx + uyy ) + cdtPow4By12*( uxxxx + uyyyy + 2.*uxxyy )  + cdtPow6By360*( uxxxxxx +  uyyyyyy +  3.*( uxxxxyy +  uxxyyyy) ) + dtSq*fv(m)      
                    ! if( i1.eq.5 .and. i2.eq.6 )then
                    !   call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t+dt,uc,ue )
                    !   write(*,'("ME8C: (i1,i2)=(",2i3,") u=",1pe12.4," ue=",1pe12.4," err=",1pe8.2)') i1,i2,un(i1,i2,i3,0),ue,abs(un(i1,i2,i3,0)-ue)
                    ! end if
                                        end if ! mask 
                                      end do
                                      end do
                                      end do
            ! #If 6 == 6 || 6 == 8 
            !   evalDerivativesRectangular()
            !   write(*,*) ' Stop here for now'
            !   stop 666
            ! #End
          !   ! --- TAYLOR TIME-STEPPING --- 
          !   m=0 ! component number 
          !   ec = 0 ! component number
          !   ! #If "curvilinear" eq "curvilinear"
          !   !   #If "6" eq "4" && "2" eq "4"
          !   !     computeLaplacianOrder2(2)
          !   !   #End
          !   ! #End
          !   if( forcingOption.eq.helmholtzForcing )then
          !     coswt = cos(omega*t)
          !   end if 
          !   fv(m)=0.
          !   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
          !     getForcing(2,6,2,curvilinear) 
          !    #If "6" eq "2"
          !      ! --- SECOND 6 ---
          !      #If "2" eq "2"
          !        ! --- TWO DIMENSIONS ---
          !        #If "curvilinear" eq "rectangular"
          !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx22r(i1,i2,i3,0) + uyy22r(i1,i2,i3,0) ) + dtSq*fv(m)
          !         ! write(*,'(" adv: i1,i2=",2i4," un,u,um=",3e12.2," cdtSq,fv=",2e12.2)') i1,i2,un(i1,i2,i3,m),u(i1,i2,i3,m),um(i1,i2,i3,m),cdtSq,fv(m)
          !        #Else
          !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx22(i1,i2,i3,0)  + uyy22(i1,i2,i3,0) ) + dtSq*fv(m)
          !        #End
          !      #Else
          !        ! --- THREE DIMENSIONS ---
          !        #If "curvilinear" eq "rectangular"
          !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx23r(i1,i2,i3,0) + uyy23r(i1,i2,i3,0) + uzz23r(i1,i2,i3,0) ) + dtSq*fv(m)
          !        #Else
          !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx23(i1,i2,i3,0)  + uyy23(i1,i2,i3,0)  + uzz23(i1,i2,i3,0)  ) + dtSq*fv(m)
          !        #End
          !      #End
          !    #Elif "6" eq "4"
          !      ! --- -FOURTH 6 ---
          !      #If "2" eq "2"
          !        ! --- FOUTH-6 TWO DIMENSIONS ---
          !        #If "2" eq "4"
          !          ! orderInSpace=4 and orderInTime=4 
          !          #If "curvilinear" eq "rectangular"
          !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + cdtsq12*lap2d2pow2(i1,i2,i3,m) + dtSq*fv(m)
          !          #Else
          !            ! v is assumed to hold Lap(u) to 2nd-order
          !            ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
          !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !                                                            + cdtsq12*( vxx22(i1,i2,i3,m) + vyy22(i1,i2,i3,m) ) + dtSq*fv(m)
          !          #End
          !        #Else
          !          ! orderInSpace==4 and orderInTime==2                                                   
          !          #If "curvilinear" eq "rectangular"
          !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + dtSq*fv(m)
          !          #Else
          !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End
          !        #End                                                          
          !      #Else
          !        ! --- FOURTH-6 THREE DIMENSIONS ---
          !        #If "2" eq "4"
          !          ! orderInSpace=4 and orderInTime=4 
          !          #If "curvilinear" eq "rectangular"
          !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + cdtsq12*lap3d2pow2(i1,i2,i3,m) + dtSq*fv(m)
          !          #Else
          !            ! v is assumed to hold Lap(u) to 2nd-order
          !            ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
          !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !                                                            + cdtsq12*( vxx23(i1,i2,i3,m) + vyy23(i1,i2,i3,m) + vzz23(i1,i2,i3,m) ) + dtSq*fv(m)
          !          #End
          !        #Else
          !          ! orderInSpace==4 and orderInTime==2                                                   
          !          #If "curvilinear" eq "rectangular"
          !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + dtSq*fv(m)
          !          #Else
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End
          !        #End     
          !      #End
          !    #Elif "6" eq "6"
          !      ! ---- SIXTH 6 ---
          !      #If "2" eq "2"
          !        ! --- SIXTH-6 TWO DIMENSIONS ---
          !        #If "2" eq "2"
          !          ! orderInSpace==6 and orderInTime==2                                                   
          !          #If "curvilinear" eq "rectangular"
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx62r(i1,i2,i3,m) + uyy62r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #Else
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx62(i1,i2,i3,m)  +  uyy62(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End       
          !        #Else
          !          stop 6666                                                           
          ! !          ! ---- MODIFIED EQUATION 6=6 2D -----
          ! !          getSixthDerivatives2d(6,curvilinear,evalMetrics,i1,i2,i3)
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
          !        ! --- SIXTH-6 THREE DIMENSIONS ---
          !        #If "2" eq "2
          !          ! orderInSpace==6 and orderInTime==2                                                   
          !          #If "curvilinear" eq "rectangular"
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx63r(i1,i2,i3,m) + uyy63r(i1,i2,i3,m) + uzz63r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #Else
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx63(i1,i2,i3,m)  + uyy63(i1,i2,i3,m)  + uzz63(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End       
          !        #Else
          !          ! MODIFIED EQUATION 6=6 3D
          !          ! Turn off for now: 
          !          ! getSixthDerivatives3d(6,curvilinear,evalMetrics,i1,i2,i3)
          !          write(*,'("advWave: order=6, orderInTime=2 FINISH ME")')
          !          stop 7777
          !          ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) !          !                 + cdtsq*( uxx + uyy +uzz ) !          !                 + cdtPow4By12*( uxxxx + uyyyy + uzzzz + 2.*( uxxyy +uxxzz + uyyzz ) )  !          !                 + cdtPow6By360*( uxxxxxx +  uyyyyyy + uzzzzzz + 3.*(uxxxxyy + uxxyyyy + uxxxxzz + uyyyyzz + uxxzzzz + uyyzzzz ) + 6.*uxxyyzz ) !          !                 + dtSq*fv(m)
          !        #End     
          !      #End
          !    #Elif "6" eq "8"
          !      ! ---- EIGTH 6 ---
          !      #If "2" eq "2"
          !        ! --- EIGTH-6 TWO DIMENSIONS ---
          !        #If "2" eq "2"
          !          ! orderInSpace==8 and orderInTime==2                                                   
          !          #If "curvilinear" eq "rectangular"
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx82r(i1,i2,i3,m) + uyy82r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #Else
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx82(i1,i2,i3,m)  +  uyy82(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End       
          !        #Else
          !          write(*,'("advWave: order=6, orderInTime=2 FINISH ME")')
          !          stop 7777
          !          ! ! orderInSpace=4 and orderInTime=4 
          !          ! #If "curvilinear" eq "rectangular"
          !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + cdtsq12*lap2d2pow2(i1,i2,i3,m) + dtSq*fv(m)
          !          ! #Else
          !          !   ! v is assumed to hold Lap(u) to 2nd-order
          !          !   ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
          !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !          !                                                   + cdtsq12*( vxx22(i1,i2,i3,m) + vyy22(i1,i2,i3,m) ) + dtSq*fv(m)
          !          ! #End
          !        #End                                                          
          !     #Else
          !        ! --- EIGTH-6 THREE DIMENSIONS ---
          !        #If "2" eq "2
          !          ! orderInSpace==8 and orderInTime==2                                                   
          !          #If "curvilinear" eq "rectangular"
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx83r(i1,i2,i3,m) + uyy83r(i1,i2,i3,m) + uzz83r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #Else
          !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx83(i1,i2,i3,m)  + uyy83(i1,i2,i3,m)  + uzz83(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
          !          #End       
          !        #Else
          !          write(*,'("advWave: order=6, orderInTime=2 FINISH ME")')
          !          stop 7777
          !          ! #If "curvilinear" eq "rectangular"
          !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + cdtsq12*lap3d2pow2(i1,i2,i3,m) + dtSq*fv(m)
          !          ! #Else
          !          !   ! v is assumed to hold Lap(u) to 2nd-order
          !          !   ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
          !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !          !                                                   + cdtsq12*( vxx23(i1,i2,i3,m) + vyy23(i1,i2,i3,m) + vzz23(i1,i2,i3,m) ) + dtSq*fv(m)
          !          ! #End
          !        #End     
          !      #End
          !    #Else
          !      write(*,'("advWave: UNKNOWN order=6")')
          !      stop 7777
          !    #End         
          !    ! write(*,'("i1,i2=",2i3," u-ue=",e10.2)') i1,i2,u(i1,i2,i3,m)-ev(m)
          !    ! write(*,'(" uxx-uxxe =",e10.2)') uxx22r(i1,i2,i3,0)-evxx(m)
          !    ! OGDERIV2D( 0,0,0,0,i1,i2,i3,t+dt, ec, ev(m)  )
          !    ! write(*,'(" un-ue=",e10.2)') un(i1,i2,i3,m)-ev(m)
          !   endLoopsMask()
                  else
                          if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
                              write(*,'("advWaveME: ADVANCE dim=2 order=6 orderInTime=6, grid=curvilinear... t=",e10.2)') t
                          end if
                          m=0 ! component number 
                          ec = 0 ! component number  
                                      numGhost1=2; ! should depend on the orderOfAccuracy
                                      n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                      n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                      n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                        do i3=n3a,n3b
                                        do i2=n2a,n2b
                                        do i1=n1a,n1b
                                          if( mask(i1,i2,i3).ne.0 )then
                                              d200(i1,i2,i3,0) = (u(i1+1,i2,i3,0) - 2.*u(i1,i2,i3,0) + u(i1-1,i2,i3,0))/(dx(0)**2)
                                              d020(i1,i2,i3,0) = (u(i1,i2+1,i3,0) - 2.*u(i1,i2,i3,0) + u(i1,i2-1,i3,0))/(dx(1)**2)
                       ! We really only need this for m2>=m1  ** FIX ME **
                                              do m1=0,nd-1
                                                  do m2=0,nd-1
                                                      rx200(i1,i2,i3,m1,m2) = (rsxy(i1+1,i2,i3,m1,m2)- 2.*rsxy(i1,i2,i3,m1,m2) + rsxy(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                                                      rx020(i1,i2,i3,m1,m2) = (rsxy(i1,i2+1,i3,m1,m2)- 2.*rsxy(i1,i2,i3,m1,m2) + rsxy(i1,i2-1,i3,m1,m2))/(dx(1)**2)          
                                                  end do
                                              end do
                                          end if ! mask .ne. 0
                                        end do
                                        end do
                                        end do
                                      numGhost1=1; ! should depend on the orderOfAccuracy
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
                       ! We really only need this for m2>=m1  ** FIX ME **
                                              do m1=0,nd-1
                                                  do m2=0,nd-1
                                                      rx400(i1,i2,i3,m1,m2) = (rx200(i1+1,i2,i3,m1,m2) - 2.*rx200(i1,i2,i3,m1,m2) + rx200(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                                                      rx220(i1,i2,i3,m1,m2) = (rx200(i1,i2+1,i3,m1,m2) - 2.*rx200(i1,i2,i3,m1,m2) + rx200(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                                                      rx040(i1,i2,i3,m1,m2) = (rx020(i1,i2+1,i3,m1,m2) - 2.*rx020(i1,i2,i3,m1,m2) + rx020(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                                                  end do
                                              end do
                                          end if ! mask .ne. 0
                                        end do
                                        end do
                                        end do
                   ! numGhost1=1; ! should depend on the orderOfAccuracy
                   ! n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                   ! n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                   ! n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                   ! beginLoops3d()
                   !   if( mask(i1,i2,i3).ne.0 )then
                   !     d600(i1,i2,i3,0) = (d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0))/(dx(0)**2)
                   !     d060(i1,i2,i3,0) = (d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0))/(dx(1)**2)
                   !     d420(i1,i2,i3,0) = (d220(i1+1,i2,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0))/(dx(0)**2)
                   !     d240(i1,i2,i3,0) = (d220(i1,i2+1,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1,i2-1,i3,0))/(dx(1)**2)
                   !     ! We really only need this for m2>=m1  ** FIX ME **
                   !     do m1=0,nd-1
                   !       do m2=0,nd-1
                   !         rx600(i1,i2,i3,m1,m2) = (rx400(i1+1,i2,i3,m1,m2) - 2.*rx400(i1,i2,i3,m1,m2) + rx400(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                   !         rx060(i1,i2,i3,m1,m2) = (rx040(i1,i2+1,i3,m1,m2) - 2.*rx040(i1,i2,i3,m1,m2) + rx040(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                   !         rx420(i1,i2,i3,m1,m2) = (rx220(i1+1,i2,i3,m1,m2) - 2.*rx220(i1,i2,i3,m1,m2) + rx220(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                   !         rx240(i1,i2,i3,m1,m2) = (rx220(i1,i2+1,i3,m1,m2) - 2.*rx220(i1,i2,i3,m1,m2) + rx220(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                   !       end do
                   !     end do      
                   !   end if ! mask .ne. 0
                   ! endLoops3d()     
                                      numGhost1=0; 
                                      n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                      n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                      if( nd.eq.2 )then
                                          n3a=gridIndexRange(0,2); n3b=gridIndexRange(1,2);
                                      else
                                          n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                      end if
                   ! See highOrderDiff.maple
                   ! u.xx = D+D-( I - b1*dx^2*D+D- + b2*(dx^2*D+D-)^2 -b3*(dx^2*D+D-)^2 + ...)
                   ! b1 = 1/12
                   ! b2 = 1/90
                   ! b3= 1/560
                   ! b4= 1/3150
                   ! ----------- MAIN LOOP : 6 6 CURVILINEAR --------
                                      fv(m)=0.
                                        do i3=n3a,n3b
                                        do i2=n2a,n2b
                                        do i1=n1a,n1b
                                          if( mask(i1,i2,i3).ne.0 )then
                                              d600i = (d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0))/(dx(0)**2)
                                              d060i = (d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0))/(dx(1)**2)
                                              d420i = (d220(i1+1,i2,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0))/(dx(0)**2)
                                              d240i = (d220(i1,i2+1,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1,i2-1,i3,0))/(dx(1)**2)
                       ! We really only need this for m2>=m1  ** FIX ME **
                                              do m1=0,nd-1
                                                  do m2=0,nd-1
                                                      rx600i(m1,m2) = (rx400(i1+1,i2,i3,m1,m2) - 2.*rx400(i1,i2,i3,m1,m2) + rx400(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                                                      rx060i(m1,m2) = (rx040(i1,i2+1,i3,m1,m2) - 2.*rx040(i1,i2,i3,m1,m2) + rx040(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                                                      rx420i(m1,m2) = (rx220(i1+1,i2,i3,m1,m2) - 2.*rx220(i1,i2,i3,m1,m2) + rx220(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                                                      rx240i(m1,m2) = (rx220(i1,i2+1,i3,m1,m2) - 2.*rx220(i1,i2,i3,m1,m2) + rx220(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                                                  end do
                                              end do       
                       ! d600i = (d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0))/(dx(0)**2)  
                       ! d060i = (d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0))/(dx(1)**2)
                       ! d800i = (d600(i1+1,i2,i3,0) - 2.*d600(i1,i2,i3,0) + d600(i1-1,i2,i3,0))/(dx(0)**2)
                       ! d080i = (d060(i1,i2+1,i3,0) - 2.*d060(i1,i2,i3,0) + d060(i1,i2-1,i3,0))/(dx(1)**2)
                       ! d620i = (d420(i1+1,i2,i3,0) - 2.*d420(i1,i2,i3,0) + d420(i1-1,i2,i3,0))/(dx(0)**2)  
                       ! d260i = (d240(i1,i2+1,i3,0) - 2.*d240(i1,i2,i3,0) + d240(i1,i2-1,i3,0))/(dx(1)**2)
                       ! d440i = (d240(i1+1,i2,i3,0) - 2.*d240(i1,i2,i3,0) + d240(i1-1,i2,i3,0))/(dx(0)**2)  
                       ! d240i = d240(i1,i2,i3,0)    
                       ! d420i = d420(i1,i2,i3,0)    
                       ! ----- for odd derivatives ----
                                              d100i = (u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0))/(2.*dx(0))   ! D0x 
                                              d010i = (u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0))/(2.*dx(1))   ! D0y 
                                              d300i = (d200(i1+1,i2,i3,0) - d200(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x D_x D-x 
                                              d030i = (d020(i1,i2+1,i3,0) - d020(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y D+y D-y 
                                              d120i = (d020(i1+1,i2,i3,0) - d020(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x D_y D-y 
                                              d210i = (d200(i1,i2+1,i3,0) - d200(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y D+x D-x       
                                              d500i = (d400(i1+1,i2,i3,0) - d400(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_x D-x)^2 
                                              d050i = (d040(i1,i2+1,i3,0) - d040(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+y D-y)^2 
                                              d140i = (d040(i1+1,i2,i3,0) - d040(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D_y D-y)^2 
                                              d410i = (d400(i1,i2+1,i3,0) - d400(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)^2
                                              d320i = (d220(i1+1,i2,i3,0) - d220(i1-1,i2,i3,0))/(2.*dx(0))  ! D0x (D+x D-x)(D_y D-y)
                                              d230i = (d220(i1,i2+1,i3,0) - d220(i1,i2-1,i3,0))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)
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
                       ! --- for mixed derivatives
                                              d110i = (u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))                ! D0x D0y
                                              d310i = (d200(i1+1,i2+1,i3,0) - d200(i1-1,i2+1,i3,0) - d200(i1+1,i2-1,i3,0) + d200(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) 
                                              d130i = (d020(i1+1,i2+1,i3,0) - d020(i1-1,i2+1,i3,0) - d020(i1+1,i2-1,i3,0) + d020(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)          
                                              d510i = (d400(i1+1,i2+1,i3,0) - d400(i1-1,i2+1,i3,0) - d400(i1+1,i2-1,i3,0) + d400(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2
                                              d150i = (d040(i1+1,i2+1,i3,0) - d040(i1-1,i2+1,i3,0) - d040(i1+1,i2-1,i3,0) + d040(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^2     
                                              d330i = (d220(i1+1,i2+1,i3,0) - d220(i1-1,i2+1,i3,0) - d220(i1+1,i2-1,i3,0) + d220(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)       
                       ! d710i = (d600(i1+1,i2+1,i3,0) - d600(i1-1,i2+1,i3,0) - d600(i1+1,i2-1,i3,0) + d600(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^3
                       ! d170i = (d060(i1+1,i2+1,i3,0) - d060(i1-1,i2+1,i3,0) - d060(i1+1,i2-1,i3,0) + d060(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^3  
                       ! d530i = (d420(i1+1,i2+1,i3,0) - d420(i1-1,i2+1,i3,0) - d420(i1+1,i2-1,i3,0) + d420(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2 (D+yD-y)         
                       ! d350i = (d240(i1+1,i2+1,i3,0) - d240(i1-1,i2+1,i3,0) - d240(i1+1,i2-1,i3,0) + d240(i1-1,i2-1,i3,0))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)^2         
                                              do m1=0,nd-1
                                                  do m2=0,nd-1
                                                      rx000i(m1,m2) = rsxy(i1,i2,i3,m1,m2)
                                                      rx100i(m1,m2) = (rsxy(i1+1,i2,i3,m1,m2)   - rsxy(i1-1,i2,i3,m1,m2))/(2.*dx(0))   ! D0x 
                                                      rx010i(m1,m2) = (rsxy(i1,i2+1,i3,m1,m2)   - rsxy(i1,i2-1,i3,m1,m2))/(2.*dx(1))   ! D0y 
                                                      rx110i(m1,m2) = (rsxy(i1+1,i2+1,i3,m1,m2) - rsxy(i1-1,i2+1,i3,m1,m2) - rsxy(i1+1,i2-1,i3,m1,m2) + rsxy(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1)) 
                                                      rx300i(m1,m2) = (rx200(i1+1,i2,i3,m1,m2) - rx200(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x D_x D-x 
                                                      rx030i(m1,m2) = (rx020(i1,i2+1,i3,m1,m2) - rx020(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y D+y D-y  
                                                      rx120i(m1,m2) = (rx020(i1+1,i2,i3,m1,m2) - rx020(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x D_y D-y 
                                                      rx210i(m1,m2) = (rx200(i1,i2+1,i3,m1,m2) - rx200(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y D+x D-x             
                                                      rx310i(m1,m2) = (rx200(i1+1,i2+1,i3,m1,m2) - rx200(i1-1,i2+1,i3,m1,m2) - rx200(i1+1,i2-1,i3,m1,m2) + rx200(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) 
                                                      rx130i(m1,m2) = (rx020(i1+1,i2+1,i3,m1,m2) - rx020(i1-1,i2+1,i3,m1,m2) - rx020(i1+1,i2-1,i3,m1,m2) + rx020(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)       
                                                      rx500i(m1,m2) = (rx400(i1+1,i2,i3,m1,m2) - rx400(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D_x D-x)^2 
                                                      rx050i(m1,m2) = (rx040(i1,i2+1,i3,m1,m2) - rx040(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+y D-y)^2 
                                                      rx140i(m1,m2) = (rx040(i1+1,i2,i3,m1,m2) - rx040(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D_y D-y)^2 
                                                      rx410i(m1,m2) = (rx400(i1,i2+1,i3,m1,m2) - rx400(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+x D-x)^2
                                                      rx320i(m1,m2) = (rx220(i1+1,i2,i3,m1,m2) - rx220(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D+x D-x)(D_y D-y)
                                                      rx230i(m1,m2) = (rx220(i1,i2+1,i3,m1,m2) - rx220(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)
                                                      rx510i(m1,m2) = (rx400(i1+1,i2+1,i3,m1,m2) - rx400(i1-1,i2+1,i3,m1,m2) - rx400(i1+1,i2-1,i3,m1,m2) + rx400(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2
                                                      rx150i(m1,m2) = (rx040(i1+1,i2+1,i3,m1,m2) - rx040(i1-1,i2+1,i3,m1,m2) - rx040(i1+1,i2-1,i3,m1,m2) + rx040(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^2   
                                                      rx330i(m1,m2) = (rx220(i1+1,i2+1,i3,m1,m2) - rx220(i1-1,i2+1,i3,m1,m2) - rx220(i1+1,i2-1,i3,m1,m2) + rx220(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)   
                           ! rx240i(m1,m2) = rx240(i1,i2,i3,m1,m2)    
                           ! rx420i(m1,m2) = rx420(i1,i2,i3,m1,m2) 
                           ! rx700i(m1,m2) = (rx600(i1+1,i2,i3,m1,m2) - rx600(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D_x D-x)^3 
                           ! rx070i(m1,m2) = (rx060(i1,i2+1,i3,m1,m2) - rx060(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+y D-y)^3 
                           ! rx160i(m1,m2) = (rx060(i1+1,i2,i3,m1,m2) - rx060(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D_y D-y)^3 
                           ! rx610i(m1,m2) = (rx600(i1,i2+1,i3,m1,m2) - rx600(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+x D-x)^3
                           ! rx520i(m1,m2) = (rx420(i1+1,i2,i3,m1,m2) - rx420(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D+x D-x)(D_y D-y)
                           ! rx250i(m1,m2) = (rx240(i1,i2+1,i3,m1,m2) - rx240(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)      
                           ! rx340i(m1,m2) = (rx240(i1+1,i2,i3,m1,m2) - rx240(i1-1,i2,i3,m1,m2))/(2.*dx(0))  ! D0x (D+x D-x)^2(D_y D-y)
                           ! rx430i(m1,m2) = (rx420(i1,i2+1,i3,m1,m2) - rx420(i1,i2-1,i3,m1,m2))/(2.*dx(1))  ! D0y (D+x D-x)(D_y D-y)^2
                           ! rx800i(m1,m2) = (rx600(i1+1,i2,i3,m1,m2) - 2.*rx600(i1,i2,i3,m1,m2) + rx600(i1-1,i2,i3,m1,m2))/(dx(0)**2)
                           ! rx080i(m1,m2) = (rx060(i1,i2+1,i3,m1,m2) - 2.*rx060(i1,i2,i3,m1,m2) + rx060(i1,i2-1,i3,m1,m2))/(dx(1)**2) 
                           ! rx620i(m1,m2) = (rx420(i1+1,i2,i3,m1,m2) - 2.*rx420(i1,i2,i3,m1,m2) + rx420(i1-1,i2,i3,m1,m2))/(dx(0)**2)  
                           ! rx260i(m1,m2) = (rx240(i1,i2+1,i3,m1,m2) - 2.*rx240(i1,i2,i3,m1,m2) + rx240(i1,i2-1,i3,m1,m2))/(dx(1)**2)
                           ! rx440i(m1,m2) = (rx240(i1+1,i2,i3,m1,m2) - 2.*rx240(i1,i2,i3,m1,m2) + rx240(i1-1,i2,i3,m1,m2))/(dx(0)**2)          
                           ! rx530i(m1,m2) = (rx420(i1+1,i2+1,i3,m1,m2) - rx420(i1-1,i2+1,i3,m1,m2) - rx420(i1+1,i2-1,i3,m1,m2) + rx420(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^2 (D+yD-y)         
                           ! rx350i(m1,m2) = (rx240(i1+1,i2+1,i3,m1,m2) - rx240(i1-1,i2+1,i3,m1,m2) - rx240(i1+1,i2-1,i3,m1,m2) + rx240(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x) (D+yD-y)^2             
                           ! rx710i(m1,m2) = (rx600(i1+1,i2+1,i3,m1,m2) - rx600(i1-1,i2+1,i3,m1,m2) - rx600(i1+1,i2-1,i3,m1,m2) + rx600(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y (D+xD-x)^3
                           ! rx170i(m1,m2) = (rx060(i1+1,i2+1,i3,m1,m2) - rx060(i1-1,i2+1,i3,m1,m2) - rx060(i1+1,i2-1,i3,m1,m2) + rx060(i1-1,i2-1,i3,m1,m2))/(4.*dx(0)*dx(1))    ! D0x D0y          (D+yD-y)^3  
                                                  end do
                                              end do  
                       ! ---- first derivatives ----
                       !   1/6, 1/30, 1/140, 1/630 
                       !   ux = D0 ( 1 - (1/6)*dx^2(D+D-) + (1/30)*dx^4 (D+D-)^2 -(1/140)*(D+D-)^3 + ... )
                                              ur = d100i - (dx(0)**2/6.)*d300i + (dx(0)**4/30.)*d500i ! - (dx(0)**6/140.)*d700i 
                                              us = d010i - (dx(1)**2/6.)*d030i + (dx(1)**4/30.)*d050i ! - (dx(1)**6/140.)*d070i  
                       ! ---- 2nd derivatives ----
                       ! 1/12 , 1/90, 1/560 
                       !   uxx = D+D- ( 1 - (1/12)*dx^2(D+D-) + (1/90)*dx^4 (D+D-)^2 -(1/560)*(D+D-)^3 + ... )
                                              urr = d200(i1,i2,i3,0) - (dx(0)**2/12.)*d400(i1,i2,i3,0) + (dx(0)**4/90.)*d600i ! - (dx(0)**6/560.)*d800i
                                              uss = d020(i1,i2,i3,0) - (dx(1)**2/12.)*d040(i1,i2,i3,0) + (dx(1)**4/90.)*d060i ! - (dx(1)**6/560.)*d080i 
                       ! uxy = D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                       !       D0y*( 1 - dy^2/6 D+yD-y + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+yD-y)^3 ) u 
                                              urs = d110i - (dx(0)**2/6.)*d310i - (dx(1)**2/6.)*d130i + (dx(0)**4/30.)*d510i + (dx(1)**4/30.)*d150i + (dx(0)**2 * dx(1)**2/36. )*d330i 
                           !  - ( (1./140.)*dx(0)**6 )*d710i - ( (1./140.)*dx(1)**6 )*d170i - (dx(0)**4*dx(1)**2/(6.*30.))*d530i - (dx(0)**2*dx(1)**4/(6.*30.))*d350i 
                       ! !   uxx = D+D- ( 1 - (1/12)*dx^2(D+D-) + (1/90)*dx^4 (D+D-)^2 -(1/560)*(D+D-)^3 + ... )
                       ! urr = d200(i1,i2,i3,0) - (dx(0)**2/12.)*d400(i1,i2,i3,0)
                       ! uss = d020(i1,i2,i3,0) - (dx(1)**2/12.)*d040(i1,i2,i3,0)
                       ! !   ux = D0 ( 1 - (1/6)*dx^2(D+D-) + (1/30)*dx^4 (D+D-)^2 -(1/140)*(D+D-)^3 + ... )
                       ! urs = d110i - (dx(0)**2/6.)*d310i - (dx(1)**2/6.)*d130i
                       ! ur = d100i  - (dx(0)**2/6.)*d300i
                       ! us = d010i  - (dx(1)**2/6.)*d030i
                       ! spatial derivatives of metrics:
                       !   rxxi(m1,m2,m3) = Dx_m3 rsxy(m1,m2)
                                              do m1=0,nd-1
                                                  do m2=0,nd-1
                                                      rxr(m1,m2,0) = rx100i(m1,m2) - (dx(0)**2/6.)*rx300i(m1,m2) + (dx(0)**4/30.)*rx500i(m1,m2) !- (dx(0)**6/140.)*rx700i(m1,m2)   
                                                      rxr(m1,m2,1) = rx010i(m1,m2) - (dx(1)**2/6.)*rx030i(m1,m2) + (dx(1)**4/30.)*rx050i(m1,m2) !- (dx(1)**6/140.)*rx070i(m1,m2)   
                                                      do m3 = 0,nd-1
                                                          rxxa(m1,m2,m3) = rx000i(0,m3)*rxr(m1,m2,0) + rx000i(1,m3)*rxr(m1,m2,1)  ! (rx).x = rx *(rx)_r + sx * (rx)_s 
                                                      end do
                                                  end do
                                              end do
                                              rx = rsxy(i1,i2,i3,0,0)
                                              ry = rsxy(i1,i2,i3,0,1)
                                              sx = rsxy(i1,i2,i3,1,0)
                                              sy = rsxy(i1,i2,i3,1,1)
                                              rxx = rxxa(0,0,0)
                                              sxx = rxxa(1,0,0)
                                              ryy = rxxa(0,1,1)
                                              syy = rxxa(1,1,1)
                                              rxy = rxxa(0,0,1) ! or rxx(0,1,0)
                                              sxy = rxxa(1,0,1) ! or rxx(1,1,0)
                           ! ------ Coefficients in expansion for ux ------
                                                      cux100 = rx
                                                      cux010 = sx
                           ! ux = cux100*ur+cux010*us
                           ! ------ Coefficients in expansion for uy ------
                                                      cuy100 = ry
                                                      cuy010 = sy
                           ! uy = cuy100*ur+cuy010*us
                       ! ux = cux100*ur+cux010*us
                       ! uy = cuy100*ur+cuy010*us      
                           ! ------ Coefficients in expansion for uxx ------
                                                      cuxx100 = rxx
                                                      cuxx200 = rx**2
                                                      cuxx010 = sxx
                                                      cuxx110 = 2.*rx*sx
                                                      cuxx020 = sx**2
                           ! uxx = cuxx100*ur+cuxx200*urr+cuxx010*us+cuxx110*urs+cuxx020*uss
                           ! ------ Coefficients in expansion for uxy ------
                                                      cuxy100 = rxy
                                                      cuxy200 = rx*ry
                                                      cuxy010 = sxy
                                                      cuxy110 = rx*sy+ry*sx
                                                      cuxy020 = sx*sy
                           ! uxy = cuxy100*ur+cuxy200*urr+cuxy010*us+cuxy110*urs+cuxy020*uss
                           ! ------ Coefficients in expansion for uyy ------
                                                      cuyy100 = ryy
                                                      cuyy200 = ry**2
                                                      cuyy010 = syy
                                                      cuyy110 = 2.*ry*sy
                                                      cuyy020 = sy**2
                           ! uyy = cuyy100*ur+cuyy200*urr+cuyy010*us+cuyy110*urs+cuyy020*uss
                                              uxx = cuxx100*ur+cuxx200*urr+cuxx010*us+cuxx110*urs+cuxx020*uss
                       ! uxy = cuxy100*ur+cuxy200*urr+cuxy010*us+cuxy110*urs+cuxy020*uss
                                              uyy = cuyy100*ur+cuyy200*urr+cuyy010*us+cuyy110*urs+cuyy020*uss      
                       ! ----- START THIRD DERIVATIVES -----
                                              do m1=0,nd-1
                                              do m2=0,nd-1
                         ! ---- 2nd parameteric derivatives of the metrics ----
                                                  rxrra(m1,m2) = rx200(i1,i2,i3,m1,m2) - (dx(0)**2/12.)*rx400(i1,i2,i3,m1,m2) + (dx(0)**4/90.)*rx600i(m1,m2) ! - (dx(0)**6/560.)*rx800i(m1,m2)
                                                  rxssa(m1,m2) = rx020(i1,i2,i3,m1,m2) - (dx(1)**2/12.)*rx040(i1,i2,i3,m1,m2) + (dx(1)**4/90.)*rx060i(m1,m2) ! - (dx(1)**6/560.)*rx080i(m1,m2) 
                                                  rxrsa(m1,m2) = rx110i(m1,m2) - (dx(0)**2/6.)*rx310i(m1,m2) - (dx(1)**2/6.)*rx130i(m1,m2) + (dx(0)**4/30.)*rx510i(m1,m2) + (dx(1)**4/30.)*rx150i(m1,m2) + (dx(0)**2 * dx(1)**2/36. )*rx330i(m1,m2) 
                               ! - ( (1./140.)*dx(0)**6 )*rx710i(m1,m2) - ( (1./140.)*dx(1)**6 )*rx170i(m1,m2) - (dx(0)**4*dx(1)**2/(6.*30.))*rx530i(m1,m2) - (dx(0)**2*dx(1)**4/(6.*30.))*rx350i(m1,m2)   
                         ! ---- 2nd spatial derivatives of the metrics ----
                                                  rxxxa(m1,m2) = cuxx100*rxr(m1,m2,0)+cuxx200*rxrra(m1,m2)+cuxx010*rxr(m1,m2,1)+cuxx110*rxrsa(m1,m2)+cuxx020*rxssa(m1,m2)
                                                  rxxya(m1,m2) = cuxy100*rxr(m1,m2,0)+cuxy200*rxrra(m1,m2)+cuxy010*rxr(m1,m2,1)+cuxy110*rxrsa(m1,m2)+cuxy020*rxssa(m1,m2)
                                                  rxyya(m1,m2) = cuyy100*rxr(m1,m2,0)+cuyy200*rxrra(m1,m2)+cuyy010*rxr(m1,m2,1)+cuyy110*rxrsa(m1,m2)+cuyy020*rxssa(m1,m2)
                                              end do
                                              end do
                                              rxxx=rxxxa(0,0); ryxx=rxxxa(0,1); sxxx=rxxxa(1,0); syxx=rxxxa(1,1); 
                                              rxxy=rxxya(0,0); ryxy=rxxya(0,1); sxxy=rxxya(1,0); syxy=rxxya(1,1); 
                                              rxyy=rxyya(0,0); ryyy=rxyya(0,1); sxyy=rxyya(1,0); syyy=rxyya(1,1); 
                       ! ---- third parametric derivatives ----
                       !   1/4, 7/120, 41/3024
                       !    uxxxx = D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ... 
                                              urrr = d300i - (dx(0)**2/4.)*d500i ! + (dx(0)**4 *7./120.)*d700i
                                              usss = d030i - (dx(1)**2/4.)*d050i ! + (dx(1)**4 *7./120.)*d070i
                       ! uxxy = D+xD-x*( 1 - dx^2/12 D+xD-x + dx^4/90*(D+xD-x)^2 - (1/560)*dx^6*(D+xD-x)^3 ) X
                       !           D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                                              urrs = d210i - (dx(0)**2/12.)*d410i - (dx(1)**2/6.)*d230i ! + (dx(0)**4/90.)*d610i  + (dx(1)**4/30.)*d250i + (dx(0)**2 * dx(1)**2/72. )*d430i
                       ! uxyy =    D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                       !        D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                              urss = d120i - (dx(1)**2/12.)*d140i - (dx(0)**2/6.)*d320i ! + (dx(1)**4/90.)*d160i  + (dx(0)**4/30.)*d520i + (dx(0)**2 * dx(1)**2/72. )*d340i 
                       ! ----- THIRD SPATIAL DERIVATIVES -----
                       ! ------ Coefficients in expansion for uxxx ------
cuxxx100 = rxxx
cuxxx200 = 3.*rx*rxx
cuxxx300 = rx**3
cuxxx010 = sxxx
cuxxx110 = 3.*rx*sxx+3.*rxx*sx
cuxxx210 = 3.*rx**2.*sx
cuxxx020 = 3.*sx*sxx
cuxxx120 = 3.*rx*sx**2
cuxxx030 = sx**3
                       ! uxxx = cuxxx100*ur+cuxxx200*urr+cuxxx300*urrr+cuxxx010*us+cuxxx110*urs+cuxxx210*urrs+cuxxx020*uss+cuxxx120*urss+cuxxx030*usss
                       ! ------ Coefficients in expansion for uxxy ------
cuxxy100 = rxxy
cuxxy200 = 2.*rx*rxy+rxx*ry
cuxxy300 = rx**2.*ry
cuxxy010 = sxxy
cuxxy110 = 2.*rx*sxy+rxx*sy+2.*rxy*sx+ry*sxx
cuxxy210 = rx**2.*sy+2.*rx*ry*sx
cuxxy020 = 2.*sx*sxy+sxx*sy
cuxxy120 = 2.*rx*sx*sy+ry*sx**2
cuxxy030 = sx**2.*sy
                       ! uxxy = cuxxy100*ur+cuxxy200*urr+cuxxy300*urrr+cuxxy010*us+cuxxy110*urs+cuxxy210*urrs+cuxxy020*uss+cuxxy120*urss+cuxxy030*usss
                       ! ------ Coefficients in expansion for uxyy ------
cuxyy100 = rxyy
cuxyy200 = rx*ryy+2.*rxy*ry
cuxyy300 = rx*ry**2
cuxyy010 = sxyy
cuxyy110 = rx*syy+2.*rxy*sy+2.*ry*sxy+ryy*sx
cuxyy210 = 2.*rx*ry*sy+ry**2.*sx
cuxyy020 = sx*syy+2.*sxy*sy
cuxyy120 = rx*sy**2+2.*ry*sx*sy
cuxyy030 = sx*sy**2
                       ! uxyy = cuxyy100*ur+cuxyy200*urr+cuxyy300*urrr+cuxyy010*us+cuxyy110*urs+cuxyy210*urrs+cuxyy020*uss+cuxyy120*urss+cuxyy030*usss
                       ! ------ Coefficients in expansion for uyyy ------
cuyyy100 = ryyy
cuyyy200 = 3.*ry*ryy
cuyyy300 = ry**3
cuyyy010 = syyy
cuyyy110 = 3.*ry*syy+3.*ryy*sy
cuyyy210 = 3.*ry**2.*sy
cuyyy020 = 3.*sy*syy
cuyyy120 = 3.*ry*sy**2
cuyyy030 = sy**3
                       ! uyyy = cuyyy100*ur+cuyyy200*urr+cuyyy300*urrr+cuyyy010*us+cuyyy110*urs+cuyyy210*urrs+cuyyy020*uss+cuyyy120*urss+cuyyy030*usss
                       ! uxxx = cuxxx100*ur+cuxxx200*urr+cuxxx300*urrr+cuxxx010*us+cuxxx110*urs+cuxxx210*urrs+cuxxx020*uss+cuxxx120*urss+cuxxx030*usss
                       ! uxxy = cuxxy100*ur+cuxxy200*urr+cuxxy300*urrr+cuxxy010*us+cuxxy110*urs+cuxxy210*urrs+cuxxy020*uss+cuxxy120*urss+cuxxy030*usss
                       ! uxyy = cuxyy100*ur+cuxyy200*urr+cuxyy300*urrr+cuxyy010*us+cuxyy110*urs+cuxyy210*urrs+cuxyy020*uss+cuxyy120*urss+cuxyy030*usss
                       ! uyyy = cuyyy100*ur+cuyyy200*urr+cuyyy300*urrr+cuyyy010*us+cuyyy110*urs+cuyyy210*urrs+cuyyy020*uss+cuyyy120*urss+cuyyy030*usss
                       ! ----- START FOURTH DERIVATIVES -----
                                              do m1=0,nd-1
                                              do m2=0,nd-1
                         ! ---- 3rd parameteric derivatives of the metrics ----
                                                  rxrrra(m1,m2) = rx300i(m1,m2) - (dx(0)**2/4.)*rx500i(m1,m2)  ! + (dx(0)**4 *7./120.)*rx700i(m1,m2)
                                                  rxsssa(m1,m2) = rx030i(m1,m2) - (dx(1)**2/4.)*rx050i(m1,m2)  ! + (dx(1)**4 *7./120.)*rx070i(m1,m2)
                                                  rxrrsa(m1,m2) = rx210i(m1,m2) - (dx(0)**2/12.)*rx410i(m1,m2) - (dx(1)**2/6.)*rx230i(m1,m2) ! + (dx(0)**4/90.)*rx610i(m1,m2) + (dx(1)**4/30.)*rx250i(m1,m2) + (dx(0)**2 * dx(1)**2/72. )*rx430i(m1,m2)
                                                  rxrssa(m1,m2) = rx120i(m1,m2) - (dx(1)**2/12.)*rx140i(m1,m2) - (dx(0)**2/6.)*rx320i(m1,m2) ! + (dx(1)**4/90.)*rx160i(m1,m2) + (dx(0)**4/30.)*rx520i(m1,m2) + (dx(0)**2 * dx(1)**2/72. )*rx340i(m1,m2) 
                         ! ---- third spatial derivatives of the metrics ----
                                                  rxxxxa(m1,m2) = cuxxx100*rxr(m1,m2,0)+cuxxx200*rxrra(m1,m2)+cuxxx300*rxrrra(m1,m2)+cuxxx010*rxr(m1,m2,1)+cuxxx110*rxrsa(m1,m2)+cuxxx210*rxrrsa(m1,m2)+cuxxx020*rxssa(m1,m2)+cuxxx120*rxrssa(m1,m2)+cuxxx030*rxsssa(m1,m2)
                                                  rxxxya(m1,m2) = cuxxy100*rxr(m1,m2,0)+cuxxy200*rxrra(m1,m2)+cuxxy300*rxrrra(m1,m2)+cuxxy010*rxr(m1,m2,1)+cuxxy110*rxrsa(m1,m2)+cuxxy210*rxrrsa(m1,m2)+cuxxy020*rxssa(m1,m2)+cuxxy120*rxrssa(m1,m2)+cuxxy030*rxsssa(m1,m2)
                                                  rxxyya(m1,m2) = cuxyy100*rxr(m1,m2,0)+cuxyy200*rxrra(m1,m2)+cuxyy300*rxrrra(m1,m2)+cuxyy010*rxr(m1,m2,1)+cuxyy110*rxrsa(m1,m2)+cuxyy210*rxrrsa(m1,m2)+cuxyy020*rxssa(m1,m2)+cuxyy120*rxrssa(m1,m2)+cuxyy030*rxsssa(m1,m2)
                                                  rxyyya(m1,m2) = cuyyy100*rxr(m1,m2,0)+cuyyy200*rxrra(m1,m2)+cuyyy300*rxrrra(m1,m2)+cuyyy010*rxr(m1,m2,1)+cuyyy110*rxrsa(m1,m2)+cuyyy210*rxrrsa(m1,m2)+cuyyy020*rxssa(m1,m2)+cuyyy120*rxrssa(m1,m2)+cuyyy030*rxsssa(m1,m2)
                                              end do
                                              end do
                                              rxxxx=rxxxxa(0,0); ryxxx=rxxxxa(0,1); sxxxx=rxxxxa(1,0); syxxx=rxxxxa(1,1); 
                                              rxxxy=rxxxya(0,0); ryxxy=rxxxya(0,1); sxxxy=rxxxya(1,0); syxxy=rxxxya(1,1); 
                                              rxxyy=rxxyya(0,0); ryxyy=rxxyya(0,1); sxxyy=rxxyya(1,0); syxyy=rxxyya(1,1); 
                                              rxyyy=rxyyya(0,0); ryyyy=rxyyya(0,1); sxyyy=rxyyya(1,0); syyyy=rxyyya(1,1); 
                       ! ----- fourth parametric derivs ------
                       !  1/6, 7/240, 41/7560
                       !   uxxxx = (D+xD-x)^2*( 1 - (1/6)*D+xD-x + (7/240)*(D+xD-x)^2 + ... )
                                              urrrr = d400(i1,i2,i3,0) - (dx(0)**2/6.)*d600i  !  + (dx(0)**4 * (7./240))*d800i
                                              ussss = d040(i1,i2,i3,0) - (dx(1)**2/6.)*d060i  !  + (dx(1)**4 * (7./240))*d080i
                       ! uxyyy = D0x*( 1 - dx^2/6 D+xD-x + dx^4/30*(D+xD-x)^2 - (1/140)*dx^6*(D+xD-x)^3 ) X
                       !         D0y D+yD-y*( 1 - (1/4)*dy^2 * D+yD-y + (7/120)*dy^4 (D+yD-y)^2 + ... 
                                              ursss = d130i - (dx(0)**2/6.)*d330i - (dx(1)**2/4.)*d150i ! + (dx(0)**4/30.)*d530i + (dx(1)**4*7./120.)*d170i + (dx(0)**2 * dx(1)**2/24. )*d350i 
                       ! uxxxy
                                              urrrs = d310i - (dx(1)**2/6.)*d330i - (dx(0)**2/4.)*d510i ! + (dx(1)**4/30.)*d350i + (dx(0)**4*7./120.)*d710i + (dx(0)**2 * dx(1)**2/24. )*d530i       
                       ! uxxyy = D+xD-x*( 1 - dx^2/12 D+xD-x + dx^4/90*(D+xD-x)^2 - (1/560)*dx^6*(D+xD-x)^3 ) X 
                       !         D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                              urrss = d220(i1,i2,i3,0) - (dx(0)**2/12.)*d420i - (dx(1)**2/12.)*d240i 
                                             !  (dx(0)**4 * (1./90))*d620i + (dx(1)**4 * (1./90))*d260i + (dx(0)**2 * dx(1)**2 * 1./144.)*d440i 
                       ! ----- FOURTH SPATIAL DERIVATIVES -----
                       ! ------ Coefficients in expansion for uxxxx ------
cuxxxx100 = rxxxx
cuxxxx200 = 4.*rx*rxxx+3.*rxx**2
cuxxxx300 = 6.*rx**2.*rxx
cuxxxx400 = rx**4
cuxxxx010 = sxxxx
cuxxxx110 = 4.*rx*sxxx+6.*rxx*sxx+4.*rxxx*sx
cuxxxx210 = 6.*rx**2.*sxx+12.*rx*rxx*sx
cuxxxx310 = 4.*rx**3.*sx
cuxxxx020 = 4.*sx*sxxx+3.*sxx**2
cuxxxx120 = 12.*rx*sx*sxx+6.*rxx*sx**2
cuxxxx220 = 6.*rx**2.*sx**2
cuxxxx030 = 6.*sx**2.*sxx
cuxxxx130 = 4.*rx*sx**3
cuxxxx040 = sx**4
                       !uxxxx = cuxxxx100*ur+cuxxxx200*urr+cuxxxx300*urrr+cuxxxx400*urrrr+cuxxxx010*us+cuxxxx110*urs+cuxxxx210*urrs+cuxxxx310*urrrs+cuxxxx020*uss+cuxxxx120*urss+cuxxxx220*urrss+cuxxxx030*usss+cuxxxx130*ursss+cuxxxx040*ussss
                       ! ------ Coefficients in expansion for uxxxy ------
cuxxxy100 = rxxxy
cuxxxy200 = 3.*rx*rxxy+3.*rxx*rxy+rxxx*ry
cuxxxy300 = 3.*rx**2.*rxy+3.*rx*rxx*ry
cuxxxy400 = rx**3.*ry
cuxxxy010 = sxxxy
cuxxxy110 = 3.*rx*sxxy+3.*rxx*sxy+rxxx*sy+3.*rxxy*sx+3.*rxy*sxx+ry*sxxx
cuxxxy210 = 3.*rx**2.*sxy+(3.*rxx*sy+6.*rxy*sx+3.*ry*sxx)*rx+3.*rxx*ry*sx
cuxxxy310 = rx**3.*sy+3.*rx**2.*ry*sx
cuxxxy020 = 3.*sx*sxxy+3.*sxx*sxy+sxxx*sy
cuxxxy120 = 3.*rxy*sx**2+(6.*rx*sxy+3.*rxx*sy+3.*ry*sxx)*sx+3.*rx*sxx*sy
cuxxxy220 = 3.*rx**2.*sx*sy+3.*rx*ry*sx**2
cuxxxy030 = 3.*sx**2.*sxy+3.*sx*sxx*sy
cuxxxy130 = 3.*rx*sx**2.*sy+ry*sx**3
cuxxxy040 = sx**3.*sy
                       !uxxxy = cuxxxy100*ur+cuxxxy200*urr+cuxxxy300*urrr+cuxxxy400*urrrr+cuxxxy010*us+cuxxxy110*urs+cuxxxy210*urrs+cuxxxy310*urrrs+cuxxxy020*uss+cuxxxy120*urss+cuxxxy220*urrss+cuxxxy030*usss+cuxxxy130*ursss+cuxxxy040*ussss
                       ! ------ Coefficients in expansion for uxxyy ------
cuxxyy100 = rxxyy
cuxxyy200 = 2.*rx*rxyy+rxx*ryy+2.*rxxy*ry+2.*rxy**2
cuxxyy300 = rx**2.*ryy+4.*rx*rxy*ry+rxx*ry**2
cuxxyy400 = rx**2.*ry**2
cuxxyy010 = sxxyy
cuxxyy110 = 2.*rx*sxyy+rxx*syy+2.*rxxy*sy+4.*rxy*sxy+2.*rxyy*sx+2.*ry*sxxy+ryy*sxx
cuxxyy210 = rx**2.*syy+(4.*rxy*sy+4.*ry*sxy+2.*ryy*sx)*rx+ry*(2.*rxx*sy+4.*rxy*sx+ry*sxx)
cuxxyy310 = 2.*rx**2.*ry*sy+2.*rx*ry**2.*sx
cuxxyy020 = 2.*sx*sxyy+sxx*syy+2.*sxxy*sy+2.*sxy**2
cuxxyy120 = ryy*sx**2+(2.*rx*syy+4.*rxy*sy+4.*ry*sxy)*sx+rxx*sy**2+(4.*rx*sxy+2.*ry*sxx)*sy
cuxxyy220 = rx**2.*sy**2+4.*rx*ry*sx*sy+ry**2.*sx**2
cuxxyy030 = sx**2.*syy+4.*sx*sxy*sy+sxx*sy**2
cuxxyy130 = 2.*rx*sx*sy**2+2.*ry*sx**2.*sy
cuxxyy040 = sx**2.*sy**2
                       !uxxyy = cuxxyy100*ur+cuxxyy200*urr+cuxxyy300*urrr+cuxxyy400*urrrr+cuxxyy010*us+cuxxyy110*urs+cuxxyy210*urrs+cuxxyy310*urrrs+cuxxyy020*uss+cuxxyy120*urss+cuxxyy220*urrss+cuxxyy030*usss+cuxxyy130*ursss+cuxxyy040*ussss
                       ! ------ Coefficients in expansion for uxyyy ------
cuxyyy100 = rxyyy
cuxyyy200 = rx*ryyy+3.*rxy*ryy+3.*rxyy*ry
cuxyyy300 = 3.*rx*ry*ryy+3.*rxy*ry**2
cuxyyy400 = rx*ry**3
cuxyyy010 = sxyyy
cuxyyy110 = rx*syyy+3.*rxy*syy+3.*rxyy*sy+3.*ry*sxyy+3.*ryy*sxy+ryyy*sx
cuxyyy210 = 3.*ry**2.*sxy+(3.*rx*syy+6.*rxy*sy+3.*ryy*sx)*ry+3.*ryy*rx*sy
cuxyyy310 = 3.*rx*ry**2.*sy+ry**3.*sx
cuxyyy020 = sx*syyy+3.*sxy*syy+3.*sxyy*sy
cuxyyy120 = 3.*rxy*sy**2+(3.*rx*syy+6.*ry*sxy+3.*ryy*sx)*sy+3.*ry*sx*syy
cuxyyy220 = 3.*rx*ry*sy**2+3.*ry**2.*sx*sy
cuxyyy030 = 3.*sx*sy*syy+3.*sxy*sy**2
cuxyyy130 = rx*sy**3+3.*ry*sx*sy**2
cuxyyy040 = sx*sy**3
                       !uxyyy = cuxyyy100*ur+cuxyyy200*urr+cuxyyy300*urrr+cuxyyy400*urrrr+cuxyyy010*us+cuxyyy110*urs+cuxyyy210*urrs+cuxyyy310*urrrs+cuxyyy020*uss+cuxyyy120*urss+cuxyyy220*urrss+cuxyyy030*usss+cuxyyy130*ursss+cuxyyy040*ussss
                       ! ------ Coefficients in expansion for uyyyy ------
cuyyyy100 = ryyyy
cuyyyy200 = 4.*ry*ryyy+3.*ryy**2
cuyyyy300 = 6.*ry**2.*ryy
cuyyyy400 = ry**4
cuyyyy010 = syyyy
cuyyyy110 = 4.*ry*syyy+6.*ryy*syy+4.*ryyy*sy
cuyyyy210 = 6.*ry**2.*syy+12.*ry*ryy*sy
cuyyyy310 = 4.*ry**3.*sy
cuyyyy020 = 4.*sy*syyy+3.*syy**2
cuyyyy120 = 12.*ry*sy*syy+6.*ryy*sy**2
cuyyyy220 = 6.*ry**2.*sy**2
cuyyyy030 = 6.*sy**2.*syy
cuyyyy130 = 4.*ry*sy**3
cuyyyy040 = sy**4
                       !uyyyy = cuyyyy100*ur+cuyyyy200*urr+cuyyyy300*urrr+cuyyyy400*urrrr+cuyyyy010*us+cuyyyy110*urs+cuyyyy210*urrs+cuyyyy310*urrrs+cuyyyy020*uss+cuyyyy120*urss+cuyyyy220*urrss+cuyyyy030*usss+cuyyyy130*ursss+cuyyyy040*ussss
                                              uxxxx = cuxxxx100*ur+cuxxxx200*urr+cuxxxx300*urrr+cuxxxx400*urrrr+cuxxxx010*us+cuxxxx110*urs+cuxxxx210*urrs+cuxxxx310*urrrs+cuxxxx020*uss+cuxxxx120*urss+cuxxxx220*urrss+cuxxxx030*usss+cuxxxx130*ursss+cuxxxx040*ussss
                       ! uxxxy = cuxxxy100*ur+cuxxxy200*urr+cuxxxy300*urrr+cuxxxy400*urrrr+cuxxxy010*us+cuxxxy110*urs+cuxxxy210*urrs+cuxxxy310*urrrs+cuxxxy020*uss+cuxxxy120*urss+cuxxxy220*urrss+cuxxxy030*usss+cuxxxy130*ursss+cuxxxy040*ussss
                                              uxxyy = cuxxyy100*ur+cuxxyy200*urr+cuxxyy300*urrr+cuxxyy400*urrrr+cuxxyy010*us+cuxxyy110*urs+cuxxyy210*urrs+cuxxyy310*urrrs+cuxxyy020*uss+cuxxyy120*urss+cuxxyy220*urrss+cuxxyy030*usss+cuxxyy130*ursss+cuxxyy040*ussss
                       ! uxyyy = cuxyyy100*ur+cuxyyy200*urr+cuxyyy300*urrr+cuxyyy400*urrrr+cuxyyy010*us+cuxyyy110*urs+cuxyyy210*urrs+cuxyyy310*urrrs+cuxyyy020*uss+cuxyyy120*urss+cuxyyy220*urrss+cuxyyy030*usss+cuxyyy130*ursss+cuxyyy040*ussss
                                              uyyyy = cuyyyy100*ur+cuyyyy200*urr+cuyyyy300*urrr+cuyyyy400*urrrr+cuyyyy010*us+cuyyyy110*urs+cuyyyy210*urrs+cuyyyy310*urrrs+cuyyyy020*uss+cuyyyy120*urss+cuyyyy220*urrss+cuyyyy030*usss+cuyyyy130*ursss+cuyyyy040*ussss
                       ! ----- START FIFTH DERIVATIVES -----
                                              do m1=0,nd-1
                                              do m2=0,nd-1
                         ! ---- 4th parameteric derivatives of the metrics ----
                                                  rxrrrra(m1,m2) = rx400(i1,i2,i3,m1,m2) - (dx(0)**2/6.)*rx600i(m1,m2) ! + (dx(0)**4 * (7./240))*rx800i(m1,m2)
                                                  rxssssa(m1,m2) = rx040(i1,i2,i3,m1,m2) - (dx(1)**2/6.)*rx060i(m1,m2) ! + (dx(1)**4 * (7./240))*rx080i(m1,m2)
                                                  rxrsssa(m1,m2) = rx130i(m1,m2) - (dx(0)**2/6.)*rx330i(m1,m2) - (dx(1)**2/4.)*rx150i(m1,m2)! + (dx(0)**4/30.)*rx530i(m1,m2) + (dx(1)**4*7./120.)*rx170i(m1,m2) + (dx(0)**2 * dx(1)**2/24. )*rx350i(m1,m2) 
                                                  rxrrrsa(m1,m2) = rx310i(m1,m2) - (dx(1)**2/6.)*rx330i(m1,m2) - (dx(0)**2/4.)*rx510i(m1,m2)! + (dx(1)**4/30.)*rx350i(m1,m2) + (dx(0)**4*7./120.)*rx710i(m1,m2) + (dx(0)**2 * dx(1)**2/24. )*rx530i(m1,m2)       
                                                  rxrrssa(m1,m2) = rx220(i1,i2,i3,m1,m2) - (dx(0)**2/12.)*rx420i(m1,m2) - (dx(1)**2/12.)*rx240i(m1,m2) 
                                              ! (dx(0)**4 * (1./90))*rx620i(m1,m2) + (dx(1)**4 * (1./90))*rx260i(m1,m2) + (dx(0)**2 * dx(1)**2 * 1./144.)*rx440i(m1,m2)         
                         ! ---- fourth spatial derivatives of the metrics ----
                                                  rxxxxxa(m1,m2) = cuxxxx100*rxr(m1,m2,0)+cuxxxx200*rxrra(m1,m2)+cuxxxx300*rxrrra(m1,m2)+cuxxxx400*rxrrrra(m1,m2)+cuxxxx010*rxr(m1,m2,1)+cuxxxx110*rxrsa(m1,m2)+cuxxxx210*rxrrsa(m1,m2)+cuxxxx310*rxrrrsa(m1,m2)+cuxxxx020*rxssa(m1,m2)+cuxxxx120*rxrssa(m1,m2)+cuxxxx220*rxrrssa(m1,m2)+cuxxxx030*rxsssa(m1,m2)+cuxxxx130*rxrsssa(m1,m2)+cuxxxx040*rxssssa(m1,m2)
                                                  rxxxxya(m1,m2) = cuxxxy100*rxr(m1,m2,0)+cuxxxy200*rxrra(m1,m2)+cuxxxy300*rxrrra(m1,m2)+cuxxxy400*rxrrrra(m1,m2)+cuxxxy010*rxr(m1,m2,1)+cuxxxy110*rxrsa(m1,m2)+cuxxxy210*rxrrsa(m1,m2)+cuxxxy310*rxrrrsa(m1,m2)+cuxxxy020*rxssa(m1,m2)+cuxxxy120*rxrssa(m1,m2)+cuxxxy220*rxrrssa(m1,m2)+cuxxxy030*rxsssa(m1,m2)+cuxxxy130*rxrsssa(m1,m2)+cuxxxy040*rxssssa(m1,m2)
                                                  rxxxyya(m1,m2) = cuxxyy100*rxr(m1,m2,0)+cuxxyy200*rxrra(m1,m2)+cuxxyy300*rxrrra(m1,m2)+cuxxyy400*rxrrrra(m1,m2)+cuxxyy010*rxr(m1,m2,1)+cuxxyy110*rxrsa(m1,m2)+cuxxyy210*rxrrsa(m1,m2)+cuxxyy310*rxrrrsa(m1,m2)+cuxxyy020*rxssa(m1,m2)+cuxxyy120*rxrssa(m1,m2)+cuxxyy220*rxrrssa(m1,m2)+cuxxyy030*rxsssa(m1,m2)+cuxxyy130*rxrsssa(m1,m2)+cuxxyy040*rxssssa(m1,m2)
                                                  rxxyyya(m1,m2) = cuxyyy100*rxr(m1,m2,0)+cuxyyy200*rxrra(m1,m2)+cuxyyy300*rxrrra(m1,m2)+cuxyyy400*rxrrrra(m1,m2)+cuxyyy010*rxr(m1,m2,1)+cuxyyy110*rxrsa(m1,m2)+cuxyyy210*rxrrsa(m1,m2)+cuxyyy310*rxrrrsa(m1,m2)+cuxyyy020*rxssa(m1,m2)+cuxyyy120*rxrssa(m1,m2)+cuxyyy220*rxrrssa(m1,m2)+cuxyyy030*rxsssa(m1,m2)+cuxyyy130*rxrsssa(m1,m2)+cuxyyy040*rxssssa(m1,m2)
                                                  rxyyyya(m1,m2) = cuyyyy100*rxr(m1,m2,0)+cuyyyy200*rxrra(m1,m2)+cuyyyy300*rxrrra(m1,m2)+cuyyyy400*rxrrrra(m1,m2)+cuyyyy010*rxr(m1,m2,1)+cuyyyy110*rxrsa(m1,m2)+cuyyyy210*rxrrsa(m1,m2)+cuyyyy310*rxrrrsa(m1,m2)+cuyyyy020*rxssa(m1,m2)+cuyyyy120*rxrssa(m1,m2)+cuyyyy220*rxrrssa(m1,m2)+cuyyyy030*rxsssa(m1,m2)+cuyyyy130*rxrsssa(m1,m2)+cuyyyy040*rxssssa(m1,m2)
                                              end do
                                              end do
                                              rxxxxx=rxxxxxa(0,0); ryxxxx=rxxxxxa(0,1); sxxxxx=rxxxxxa(1,0); syxxxx=rxxxxxa(1,1); 
                                              rxxxxy=rxxxxya(0,0); ryxxxy=rxxxxya(0,1); sxxxxy=rxxxxya(1,0); syxxxy=rxxxxya(1,1); 
                                              rxxxyy=rxxxyya(0,0); ryxxyy=rxxxyya(0,1); sxxxyy=rxxxyya(1,0); syxxyy=rxxxyya(1,1); 
                                              rxxyyy=rxxyyya(0,0); ryxyyy=rxxyyya(0,1); sxxyyy=rxxyyya(1,0); syxyyy=rxxyyya(1,1); 
                                              rxyyyy=rxyyyya(0,0); ryyyyy=rxyyyya(0,1); sxyyyy=rxyyyya(1,0); syyyyy=rxyyyya(1,1); 
                       ! ---- FIXTH parametric derivatives ----
                       !   1/3, 13/144, 139/6048
                       !  uxxxxx = D0x(D+xD-x)^2 *( 1 - (1/3)*dx^2 D+xD-x + ... )
                                              urrrrr = d500i ! - (dx(0)**2/3.)*d700i
                                              usssss = d050i ! - (dx(1)**2/3.)*d070i 
                       ! uxxxxy = (D+xD-x)^2 *( 1 - dx^2/6 D+xD-x + ...) X
                       !                  D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                                              urrrrs = d410i ! - (dx(0)**2/6.)*d610i - (dx(1)**2/6.)*d430i
                                              urssss = d140i ! - (dx(1)**2/6.)*d160i - (dx(0)**2/6.)*d340i    
                       ! uxxxyy =  D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ...
                       !               D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                              urrrss = d320i ! - (dx(0)**2/4.)*d520i - (dx(1)**2/12.)*d340i
                                              urrsss = d230i ! - (dx(1)**2/4.)*d250i - (dx(0)**2/12.)*d430i
                        ! ----- FIXTH SPATIAL DERIVATIVES -----
                       ! ------ Coefficients in expansion for uxxxxx ------
cuxxxxx100 = rxxxxx
cuxxxxx200 = 5.*rx*rxxxx+10.*rxx*rxxx
cuxxxxx300 = 10.*rx**2.*rxxx+15.*rx*rxx**2
cuxxxxx400 = 10.*rx**3.*rxx
cuxxxxx500 = rx**5
cuxxxxx010 = sxxxxx
cuxxxxx110 = 5.*rx*sxxxx+10.*rxx*sxxx+10.*rxxx*sxx+5.*rxxxx*sx
cuxxxxx210 = 10.*rx**2.*sxxx+30.*rx*rxx*sxx+20.*rx*rxxx*sx+15.*rxx**2.*sx
cuxxxxx310 = 10.*rx**3.*sxx+30.*rx**2.*rxx*sx
cuxxxxx410 = 5.*sx*rx**4
cuxxxxx020 = 5.*sx*sxxxx+10.*sxx*sxxx
cuxxxxx120 = 20.*rx*sx*sxxx+15.*rx*sxx**2+30.*rxx*sx*sxx+10.*rxxx*sx**2
cuxxxxx220 = 30.*rx**2.*sx*sxx+30.*rx*rxx*sx**2
cuxxxxx320 = 10.*rx**3.*sx**2
cuxxxxx030 = 10.*sx**2.*sxxx+15.*sx*sxx**2
cuxxxxx130 = 30.*rx*sx**2.*sxx+10.*rxx*sx**3
cuxxxxx230 = 10.*rx**2.*sx**3
cuxxxxx040 = 10.*sx**3.*sxx
cuxxxxx140 = 5.*rx*sx**4
cuxxxxx050 = sx**5
                       ! uxxxxx = cuxxxxx100*ur+cuxxxxx200*urr+cuxxxxx300*urrr+cuxxxxx400*urrrr+cuxxxxx500*urrrrr+cuxxxxx010*us+cuxxxxx110*urs+cuxxxxx210*urrs+cuxxxxx310*urrrs+cuxxxxx410*urrrrs+cuxxxxx020*uss+cuxxxxx120*urss+cuxxxxx220*urrss+cuxxxxx320*urrrss+cuxxxxx030*usss+cuxxxxx130*ursss+cuxxxxx230*urrsss+cuxxxxx040*ussss+cuxxxxx140*urssss+cuxxxxx050*usssss
                       ! ------ Coefficients in expansion for uxxxxy ------
cuxxxxy100 = rxxxxy
cuxxxxy200 = 4.*rx*rxxxy+6.*rxx*rxxy+4.*rxxx*rxy+rxxxx*ry
cuxxxxy300 = 6.*rx**2.*rxxy+12.*rx*rxx*rxy+4.*rx*rxxx*ry+3.*rxx**2.*ry
cuxxxxy400 = 4.*(rxy*rx+3/2.*rxx*ry)*rx**2
cuxxxxy500 = ry*rx**4
cuxxxxy010 = sxxxxy
cuxxxxy110 = 4.*rx*sxxxy+6.*rxx*sxxy+4.*rxxx*sxy+rxxxx*sy+4.*rxxxy*sx+6.*rxxy*sxx+4.*rxy*sxxx+ry*sxxxx
cuxxxxy210 = 6.*rx**2.*sxxy+(12.*rxx*sxy+4.*rxxx*sy+12.*rxxy*sx+12.*rxy*sxx+4.*ry*sxxx)*rx+3.*rxx**2.*sy+(12.*rxy*sx+6.*ry*sxx)*rxx+4.*rxxx*ry*sx
cuxxxxy310 = 4.*rx**3.*sxy+6.*rx**2.*rxx*sy+12.*rx**2.*rxy*sx+6.*rx**2.*ry*sxx+12.*rx*rxx*ry*sx
cuxxxxy410 = rx**4.*sy+4.*rx**3.*ry*sx
cuxxxxy020 = 4.*sx*sxxxy+6.*sxx*sxxy+4.*sxxx*sxy+sxxxx*sy
cuxxxxy120 = 6.*rxxy*sx**2+(12.*rx*sxxy+12.*rxx*sxy+4.*rxxx*sy+12.*rxy*sxx+4.*ry*sxxx)*sx+3.*sxx**2.*ry+(12.*rx*sxy+6.*rxx*sy)*sxx+4.*rx*sxxx*sy
cuxxxxy220 = (12.*sx*sxy+6.*sxx*sy)*rx**2+12.*sx*(rxx*sy+rxy*sx+ry*sxx)*rx+6.*rxx*ry*sx**2
cuxxxxy320 = 4.*(sy*rx+3/2.*ry*sx)*rx**2.*sx
cuxxxxy030 = 6.*sx**2.*sxxy+12.*sx*sxx*sxy+4.*sx*sxxx*sy+3.*sxx**2.*sy
cuxxxxy130 = 12.*rx*sx**2.*sxy+12.*rx*sx*sxx*sy+6.*rxx*sx**2.*sy+4.*rxy*sx**3+6.*ry*sx**2.*sxx
cuxxxxy230 = 6.*rx*sx**2.*(sy*rx+2/3.*ry*sx)
cuxxxxy040 = 4.*(sx*sxy+3/2.*sxx*sy)*sx**2
cuxxxxy140 = 4.*rx*sx**3.*sy+ry*sx**4
cuxxxxy050 = sy*sx**4
                       ! uxxxxy = cuxxxxy100*ur+cuxxxxy200*urr+cuxxxxy300*urrr+cuxxxxy400*urrrr+cuxxxxy500*urrrrr+cuxxxxy010*us+cuxxxxy110*urs+cuxxxxy210*urrs+cuxxxxy310*urrrs+cuxxxxy410*urrrrs+cuxxxxy020*uss+cuxxxxy120*urss+cuxxxxy220*urrss+cuxxxxy320*urrrss+cuxxxxy030*usss+cuxxxxy130*ursss+cuxxxxy230*urrsss+cuxxxxy040*ussss+cuxxxxy140*urssss+cuxxxxy050*usssss
                       ! ------ Coefficients in expansion for uxxxyy ------
cuxxxyy100 = rxxxyy
cuxxxyy200 = 3.*rx*rxxyy+3.*rxx*rxyy+rxxx*ryy+2.*rxxxy*ry+6.*rxxy*rxy
cuxxxyy300 = 3.*rxyy*rx**2+(3.*rxx*ryy+6.*rxxy*ry+6.*rxy**2)*rx+rxxx*ry**2+6.*rxy*rxx*ry
cuxxxyy400 = rx**3.*ryy+6.*rx**2.*rxy*ry+3.*rx*rxx*ry**2
cuxxxyy500 = rx**3.*ry**2
cuxxxyy010 = sxxxyy
cuxxxyy110 = 3.*rx*sxxyy+3.*rxx*sxyy+rxxx*syy+2.*rxxxy*sy+6.*rxxy*sxy+3.*rxxyy*sx+6.*rxy*sxxy+3.*rxyy*sxx+2.*ry*sxxxy+ryy*sxxx
cuxxxyy210 = 3.*rx**2.*sxyy+(3.*rxx*syy+6.*rxxy*sy+12.*rxy*sxy+6.*rxyy*sx+6.*ry*sxxy+3.*ryy*sxx)*rx+ry**2.*sxxx+(6.*rxx*sxy+2.*rxxx*sy+6.*rxxy*sx+6.*rxy*sxx)*ry+(6.*rxy*sy+3.*ryy*sx)*rxx+6.*rxy**2.*sx
cuxxxyy310 = rx**3.*syy+(6.*rxy*sy+6.*ry*sxy+3.*ryy*sx)*rx**2+3.*ry*(2.*rxx*sy+4.*rxy*sx+ry*sxx)*rx+3.*rxx*ry**2.*sx
cuxxxyy410 = 2.*rx**3.*ry*sy+3.*rx**2.*ry**2.*sx
cuxxxyy020 = 3.*sx*sxxyy+3.*sxx*sxyy+sxxx*syy+2.*sxxxy*sy+6.*sxxy*sxy
cuxxxyy120 = 3.*rxyy*sx**2+(6.*rx*sxyy+3.*rxx*syy+6.*rxxy*sy+12.*rxy*sxy+6.*ry*sxxy+3.*ryy*sxx)*sx+rxxx*sy**2+(6.*rx*sxxy+6.*rxx*sxy+6.*rxy*sxx+2.*ry*sxxx)*sy+(3.*rx*syy+6.*ry*sxy)*sxx+6.*rx*sxy**2
cuxxxyy220 = (3.*sx*syy+6.*sxy*sy)*rx**2+(3.*ryy*sx**2+(12.*rxy*sy+12.*ry*sxy)*sx+6.*sxx*sy*ry+3.*rxx*sy**2)*rx+3.*ry*sx*(2.*rxx*sy+2.*rxy*sx+ry*sxx)
cuxxxyy320 = rx**3.*sy**2+6.*rx**2.*ry*sx*sy+3.*rx*ry**2.*sx**2
cuxxxyy030 = 3.*sx**2.*sxyy+(3.*sxx*syy+6.*sxxy*sy+6.*sxy**2)*sx+sxxx*sy**2+6.*sxx*sxy*sy
cuxxxyy130 = ryy*sx**3+(3.*rx*syy+6.*rxy*sy+6.*ry*sxy)*sx**2+(3.*rxx*sy**2+(12.*rx*sxy+6.*ry*sxx)*sy)*sx+3.*rx*sxx*sy**2
cuxxxyy230 = 3.*rx**2.*sx*sy**2+6.*rx*ry*sx**2.*sy+ry**2.*sx**3
cuxxxyy040 = sx**3.*syy+6.*sx**2.*sxy*sy+3.*sx*sxx*sy**2
cuxxxyy140 = 3.*rx*sx**2.*sy**2+2.*ry*sx**3.*sy
cuxxxyy050 = sx**3.*sy**2
                       ! uxxxyy = cuxxxyy100*ur+cuxxxyy200*urr+cuxxxyy300*urrr+cuxxxyy400*urrrr+cuxxxyy500*urrrrr+cuxxxyy010*us+cuxxxyy110*urs+cuxxxyy210*urrs+cuxxxyy310*urrrs+cuxxxyy410*urrrrs+cuxxxyy020*uss+cuxxxyy120*urss+cuxxxyy220*urrss+cuxxxyy320*urrrss+cuxxxyy030*usss+cuxxxyy130*ursss+cuxxxyy230*urrsss+cuxxxyy040*ussss+cuxxxyy140*urssss+cuxxxyy050*usssss
                       ! ------ Coefficients in expansion for uxxyyy ------
cuxxyyy100 = rxxyyy
cuxxyyy200 = 2.*rx*rxyyy+rxx*ryyy+3.*rxxy*ryy+3.*rxxyy*ry+6.*rxy*rxyy
cuxxyyy300 = 3.*rxxy*ry**2+(6.*rx*rxyy+3.*rxx*ryy+6.*rxy**2)*ry+ryyy*rx**2+6.*rxy*ryy*rx
cuxxyyy400 = 3.*rx**2.*ry*ryy+6.*rx*rxy*ry**2+rxx*ry**3
cuxxyyy500 = rx**2.*ry**3
cuxxyyy010 = sxxyyy
cuxxyyy110 = 2.*rx*sxyyy+rxx*syyy+3.*rxxy*syy+3.*rxxyy*sy+6.*rxy*sxyy+6.*rxyy*sxy+2.*rxyyy*sx+3.*ry*sxxyy+3.*ryy*sxxy+ryyy*sxx
cuxxyyy210 = 3.*ry**2.*sxxy+(6.*rx*sxyy+3.*rxx*syy+6.*rxxy*sy+12.*rxy*sxy+6.*rxyy*sx+3.*ryy*sxx)*ry+rx**2.*syyy+(6.*rxy*syy+6.*rxyy*sy+6.*ryy*sxy+2.*ryyy*sx)*rx+6.*rxy*ryy*sx+3.*rxx*ryy*sy+6.*rxy**2.*sy
cuxxyyy310 = ry**3.*sxx+(6.*rx*sxy+3.*rxx*sy+6.*rxy*sx)*ry**2+3.*rx*(rx*syy+4.*rxy*sy+2.*ryy*sx)*ry+3.*ryy*rx**2.*sy
cuxxyyy410 = 3.*rx**2.*ry**2.*sy+2.*rx*ry**3.*sx
cuxxyyy020 = 2.*sx*sxyyy+sxx*syyy+3.*sxxy*syy+3.*sxxyy*sy+6.*sxy*sxyy
cuxxyyy120 = 3.*rxxy*sy**2+(6.*rx*sxyy+3.*rxx*syy+12.*rxy*sxy+6.*rxyy*sx+6.*ry*sxxy+3.*ryy*sxx)*sy+ryyy*sx**2+(2.*rx*syyy+6.*rxy*syy+6.*ry*sxyy+6.*ryy*sxy)*sx+6.*rx*sxy*syy+3.*ry*sxx*syy+6.*ry*sxy**2
cuxxyyy220 = (6.*sx*sxy+3.*sxx*sy)*ry**2+(3.*rxx*sy**2+(12.*rx*sxy+12.*rxy*sx)*sy+6.*sx*syy*rx+3.*ryy*sx**2)*ry+3.*rx*sy*(rx*syy+2.*rxy*sy+2.*ryy*sx)
cuxxyyy320 = 3.*rx**2.*ry*sy**2+6.*rx*ry**2.*sx*sy+ry**3.*sx**2
cuxxyyy030 = 3.*sxxy*sy**2+(6.*sx*sxyy+3.*sxx*syy+6.*sxy**2)*sy+sx**2.*syyy+6.*sxy*syy*sx
cuxxyyy130 = rxx*sy**3+(6.*rx*sxy+6.*rxy*sx+3.*ry*sxx)*sy**2+6.*(rx*syy+2.*ry*sxy+1/2.*ryy*sx)*sx*sy+3.*ry*sx**2.*syy
cuxxyyy230 = rx**2.*sy**3+6.*rx*ry*sx*sy**2+3.*ry**2.*sx**2.*sy
cuxxyyy040 = 3.*sx**2.*sy*syy+6.*sx*sxy*sy**2+sxx*sy**3
cuxxyyy140 = 2.*rx*sx*sy**3+3.*ry*sx**2.*sy**2
cuxxyyy050 = sx**2.*sy**3
                       ! uxxyyy = cuxxyyy100*ur+cuxxyyy200*urr+cuxxyyy300*urrr+cuxxyyy400*urrrr+cuxxyyy500*urrrrr+cuxxyyy010*us+cuxxyyy110*urs+cuxxyyy210*urrs+cuxxyyy310*urrrs+cuxxyyy410*urrrrs+cuxxyyy020*uss+cuxxyyy120*urss+cuxxyyy220*urrss+cuxxyyy320*urrrss+cuxxyyy030*usss+cuxxyyy130*ursss+cuxxyyy230*urrsss+cuxxyyy040*ussss+cuxxyyy140*urssss+cuxxyyy050*usssss
                       ! ------ Coefficients in expansion for uxyyyy ------
cuxyyyy100 = rxyyyy
cuxyyyy200 = rx*ryyyy+4.*rxy*ryyy+6.*rxyy*ryy+4.*rxyyy*ry
cuxyyyy300 = 4.*rx*ry*ryyy+3.*rx*ryy**2+12.*rxy*ry*ryy+6.*rxyy*ry**2
cuxyyyy400 = 6.*ry**2.*(rx*ryy+2/3.*ry*rxy)
cuxyyyy500 = rx*ry**4
cuxyyyy010 = sxyyyy
cuxyyyy110 = rx*syyyy+4.*rxy*syyy+6.*rxyy*syy+4.*rxyyy*sy+4.*ry*sxyyy+6.*ryy*sxyy+4.*ryyy*sxy+ryyyy*sx
cuxyyyy210 = 6.*ry**2.*sxyy+(4.*rx*syyy+12.*rxy*syy+12.*rxyy*sy+12.*ryy*sxy+4.*ryyy*sx)*ry+3.*ryy**2.*sx+(6.*rx*syy+12.*rxy*sy)*ryy+4.*rx*ryyy*sy
cuxyyyy310 = 6.*rx*ry**2.*syy+12.*rx*ry*ryy*sy+12.*rxy*ry**2.*sy+4.*ry**3.*sxy+6.*ry**2.*ryy*sx
cuxyyyy410 = 4.*rx*ry**3.*sy+ry**4.*sx
cuxyyyy020 = sx*syyyy+4.*sxy*syyy+6.*sxyy*syy+4.*sxyyy*sy
cuxyyyy120 = 6.*rxyy*sy**2+(4.*rx*syyy+12.*rxy*syy+12.*ry*sxyy+12.*ryy*sxy+4.*ryyy*sx)*sy+3.*syy**2.*rx+(12.*ry*sxy+6.*ryy*sx)*syy+4.*ry*sx*syyy
cuxyyyy220 = (6.*sx*syy+12.*sxy*sy)*ry**2+12.*sy*(rx*syy+rxy*sy+ryy*sx)*ry+6.*ryy*rx*sy**2
cuxyyyy320 = 6.*ry**2.*sy*(rx*sy+2/3.*ry*sx)
cuxyyyy030 = 4.*sx*sy*syyy+3.*sx*syy**2+12.*sxy*sy*syy+6.*sxyy*sy**2
cuxyyyy130 = 6.*rx*sy**2.*syy+4.*rxy*sy**3+12.*ry*sx*sy*syy+12.*ry*sxy*sy**2+6.*ryy*sx*sy**2
cuxyyyy230 = 4.*ry*sy**2.*(rx*sy+3/2.*ry*sx)
cuxyyyy040 = 6.*sy**2.*(sx*syy+2/3.*sxy*sy)
cuxyyyy140 = rx*sy**4+4.*ry*sx*sy**3
cuxyyyy050 = sx*sy**4
                       ! uxyyyy = cuxyyyy100*ur+cuxyyyy200*urr+cuxyyyy300*urrr+cuxyyyy400*urrrr+cuxyyyy500*urrrrr+cuxyyyy010*us+cuxyyyy110*urs+cuxyyyy210*urrs+cuxyyyy310*urrrs+cuxyyyy410*urrrrs+cuxyyyy020*uss+cuxyyyy120*urss+cuxyyyy220*urrss+cuxyyyy320*urrrss+cuxyyyy030*usss+cuxyyyy130*ursss+cuxyyyy230*urrsss+cuxyyyy040*ussss+cuxyyyy140*urssss+cuxyyyy050*usssss
                       ! ------ Coefficients in expansion for uyyyyy ------
cuyyyyy100 = ryyyyy
cuyyyyy200 = 5.*ry*ryyyy+10.*ryy*ryyy
cuyyyyy300 = 10.*ry**2.*ryyy+15.*ry*ryy**2
cuyyyyy400 = 10.*ry**3.*ryy
cuyyyyy500 = ry**5
cuyyyyy010 = syyyyy
cuyyyyy110 = 5.*ry*syyyy+10.*ryy*syyy+10.*ryyy*syy+5.*ryyyy*sy
cuyyyyy210 = 10.*ry**2.*syyy+30.*ry*ryy*syy+20.*ry*ryyy*sy+15.*ryy**2.*sy
cuyyyyy310 = 10.*ry**3.*syy+30.*ry**2.*ryy*sy
cuyyyyy410 = 5.*sy*ry**4
cuyyyyy020 = 5.*sy*syyyy+10.*syy*syyy
cuyyyyy120 = 20.*ry*sy*syyy+15.*ry*syy**2+30.*ryy*sy*syy+10.*ryyy*sy**2
cuyyyyy220 = 30.*ry**2.*sy*syy+30.*ry*ryy*sy**2
cuyyyyy320 = 10.*ry**3.*sy**2
cuyyyyy030 = 10.*sy**2.*syyy+15.*sy*syy**2
cuyyyyy130 = 30.*ry*sy**2.*syy+10.*ryy*sy**3
cuyyyyy230 = 10.*ry**2.*sy**3
cuyyyyy040 = 10.*sy**3.*syy
cuyyyyy140 = 5.*ry*sy**4
cuyyyyy050 = sy**5
                       ! uyyyyy = cuyyyyy100*ur+cuyyyyy200*urr+cuyyyyy300*urrr+cuyyyyy400*urrrr+cuyyyyy500*urrrrr+cuyyyyy010*us+cuyyyyy110*urs+cuyyyyy210*urrs+cuyyyyy310*urrrs+cuyyyyy410*urrrrs+cuyyyyy020*uss+cuyyyyy120*urss+cuyyyyy220*urrss+cuyyyyy320*urrrss+cuyyyyy030*usss+cuyyyyy130*ursss+cuyyyyy230*urrsss+cuyyyyy040*ussss+cuyyyyy140*urssss+cuyyyyy050*usssss
                       ! uxxxxx = cuxxxxx100*ur+cuxxxxx200*urr+cuxxxxx300*urrr+cuxxxxx400*urrrr+cuxxxxx500*urrrrr+cuxxxxx010*us+cuxxxxx110*urs+cuxxxxx210*urrs+cuxxxxx310*urrrs+cuxxxxx410*urrrrs+cuxxxxx020*uss+cuxxxxx120*urss+cuxxxxx220*urrss+cuxxxxx320*urrrss+cuxxxxx030*usss+cuxxxxx130*ursss+cuxxxxx230*urrsss+cuxxxxx040*ussss+cuxxxxx140*urssss+cuxxxxx050*usssss
                       ! uxxxxy = cuxxxxy100*ur+cuxxxxy200*urr+cuxxxxy300*urrr+cuxxxxy400*urrrr+cuxxxxy500*urrrrr+cuxxxxy010*us+cuxxxxy110*urs+cuxxxxy210*urrs+cuxxxxy310*urrrs+cuxxxxy410*urrrrs+cuxxxxy020*uss+cuxxxxy120*urss+cuxxxxy220*urrss+cuxxxxy320*urrrss+cuxxxxy030*usss+cuxxxxy130*ursss+cuxxxxy230*urrsss+cuxxxxy040*ussss+cuxxxxy140*urssss+cuxxxxy050*usssss
                       ! uxxxyy = cuxxxyy100*ur+cuxxxyy200*urr+cuxxxyy300*urrr+cuxxxyy400*urrrr+cuxxxyy500*urrrrr+cuxxxyy010*us+cuxxxyy110*urs+cuxxxyy210*urrs+cuxxxyy310*urrrs+cuxxxyy410*urrrrs+cuxxxyy020*uss+cuxxxyy120*urss+cuxxxyy220*urrss+cuxxxyy320*urrrss+cuxxxyy030*usss+cuxxxyy130*ursss+cuxxxyy230*urrsss+cuxxxyy040*ussss+cuxxxyy140*urssss+cuxxxyy050*usssss
                       ! uxxyyy = cuxxyyy100*ur+cuxxyyy200*urr+cuxxyyy300*urrr+cuxxyyy400*urrrr+cuxxyyy500*urrrrr+cuxxyyy010*us+cuxxyyy110*urs+cuxxyyy210*urrs+cuxxyyy310*urrrs+cuxxyyy410*urrrrs+cuxxyyy020*uss+cuxxyyy120*urss+cuxxyyy220*urrss+cuxxyyy320*urrrss+cuxxyyy030*usss+cuxxyyy130*ursss+cuxxyyy230*urrsss+cuxxyyy040*ussss+cuxxyyy140*urssss+cuxxyyy050*usssss
                       ! uxyyyy = cuxyyyy100*ur+cuxyyyy200*urr+cuxyyyy300*urrr+cuxyyyy400*urrrr+cuxyyyy500*urrrrr+cuxyyyy010*us+cuxyyyy110*urs+cuxyyyy210*urrs+cuxyyyy310*urrrs+cuxyyyy410*urrrrs+cuxyyyy020*uss+cuxyyyy120*urss+cuxyyyy220*urrss+cuxyyyy320*urrrss+cuxyyyy030*usss+cuxyyyy130*ursss+cuxyyyy230*urrsss+cuxyyyy040*ussss+cuxyyyy140*urssss+cuxyyyy050*usssss
                       ! uyyyyy = cuyyyyy100*ur+cuyyyyy200*urr+cuyyyyy300*urrr+cuyyyyy400*urrrr+cuyyyyy500*urrrrr+cuyyyyy010*us+cuyyyyy110*urs+cuyyyyy210*urrs+cuyyyyy310*urrrs+cuyyyyy410*urrrrs+cuyyyyy020*uss+cuyyyyy120*urss+cuyyyyy220*urrss+cuyyyyy320*urrrss+cuyyyyy030*usss+cuyyyyy130*ursss+cuyyyyy230*urrsss+cuyyyyy040*ussss+cuyyyyy140*urssss+cuyyyyy050*usssss
                       ! ----- START SIXTH DERIVATIVES -----
                                              do m1=0,nd-1
                                              do m2=0,nd-1
                         ! ---- 5th parameteric derivatives of the metrics ----
                                                  rxrrrrra(m1,m2) = rx500i(m1,m2) !- (dx(0)**2/3.)*rx700i(m1,m2)
                                                  rxsssssa(m1,m2) = rx050i(m1,m2) !- (dx(1)**2/3.)*rx070i(m1,m2) 
                                                  rxrrrrsa(m1,m2) = rx410i(m1,m2) !- (dx(0)**2/6.)*rx610i(m1,m2) -  (dx(1)**2/6.)*rx430i(m1,m2)
                                                  rxrssssa(m1,m2) = rx140i(m1,m2) !- (dx(1)**2/6.)*rx160i(m1,m2) -  (dx(0)**2/6.)*rx340i(m1,m2)    
                                                  rxrrrssa(m1,m2) = rx320i(m1,m2) !- (dx(0)**2/4.)*rx520i(m1,m2) - (dx(1)**2/12.)*rx340i(m1,m2)
                                                  rxrrsssa(m1,m2) = rx230i(m1,m2) !- (dx(1)**2/4.)*rx250i(m1,m2) - (dx(0)**2/12.)*rx430i(m1,m2)
                         ! ---- fixth spatial derivatives of the metrics ----
                                                  rxxxxxxa(m1,m2) = cuxxxxx100*rxr(m1,m2,0)+cuxxxxx200*rxrra(m1,m2)+cuxxxxx300*rxrrra(m1,m2)+cuxxxxx400*rxrrrra(m1,m2)+cuxxxxx500*rxrrrrra(m1,m2)+cuxxxxx010*rxr(m1,m2,1)+cuxxxxx110*rxrsa(m1,m2)+cuxxxxx210*rxrrsa(m1,m2)+cuxxxxx310*rxrrrsa(m1,m2)+cuxxxxx410*rxrrrrsa(m1,m2)+cuxxxxx020*rxssa(m1,m2)+cuxxxxx120*rxrssa(m1,m2)+cuxxxxx220*rxrrssa(m1,m2)+cuxxxxx320*rxrrrssa(m1,m2)+cuxxxxx030*rxsssa(m1,m2)+cuxxxxx130*rxrsssa(m1,m2)+cuxxxxx230*rxrrsssa(m1,m2)+cuxxxxx040*rxssssa(m1,m2)+cuxxxxx140*rxrssssa(m1,m2)+cuxxxxx050*rxsssssa(m1,m2)
                                                  rxxxxxya(m1,m2) = cuxxxxy100*rxr(m1,m2,0)+cuxxxxy200*rxrra(m1,m2)+cuxxxxy300*rxrrra(m1,m2)+cuxxxxy400*rxrrrra(m1,m2)+cuxxxxy500*rxrrrrra(m1,m2)+cuxxxxy010*rxr(m1,m2,1)+cuxxxxy110*rxrsa(m1,m2)+cuxxxxy210*rxrrsa(m1,m2)+cuxxxxy310*rxrrrsa(m1,m2)+cuxxxxy410*rxrrrrsa(m1,m2)+cuxxxxy020*rxssa(m1,m2)+cuxxxxy120*rxrssa(m1,m2)+cuxxxxy220*rxrrssa(m1,m2)+cuxxxxy320*rxrrrssa(m1,m2)+cuxxxxy030*rxsssa(m1,m2)+cuxxxxy130*rxrsssa(m1,m2)+cuxxxxy230*rxrrsssa(m1,m2)+cuxxxxy040*rxssssa(m1,m2)+cuxxxxy140*rxrssssa(m1,m2)+cuxxxxy050*rxsssssa(m1,m2)
                                                  rxxxxyya(m1,m2) = cuxxxyy100*rxr(m1,m2,0)+cuxxxyy200*rxrra(m1,m2)+cuxxxyy300*rxrrra(m1,m2)+cuxxxyy400*rxrrrra(m1,m2)+cuxxxyy500*rxrrrrra(m1,m2)+cuxxxyy010*rxr(m1,m2,1)+cuxxxyy110*rxrsa(m1,m2)+cuxxxyy210*rxrrsa(m1,m2)+cuxxxyy310*rxrrrsa(m1,m2)+cuxxxyy410*rxrrrrsa(m1,m2)+cuxxxyy020*rxssa(m1,m2)+cuxxxyy120*rxrssa(m1,m2)+cuxxxyy220*rxrrssa(m1,m2)+cuxxxyy320*rxrrrssa(m1,m2)+cuxxxyy030*rxsssa(m1,m2)+cuxxxyy130*rxrsssa(m1,m2)+cuxxxyy230*rxrrsssa(m1,m2)+cuxxxyy040*rxssssa(m1,m2)+cuxxxyy140*rxrssssa(m1,m2)+cuxxxyy050*rxsssssa(m1,m2)
                                                  rxxxyyya(m1,m2) = cuxxyyy100*rxr(m1,m2,0)+cuxxyyy200*rxrra(m1,m2)+cuxxyyy300*rxrrra(m1,m2)+cuxxyyy400*rxrrrra(m1,m2)+cuxxyyy500*rxrrrrra(m1,m2)+cuxxyyy010*rxr(m1,m2,1)+cuxxyyy110*rxrsa(m1,m2)+cuxxyyy210*rxrrsa(m1,m2)+cuxxyyy310*rxrrrsa(m1,m2)+cuxxyyy410*rxrrrrsa(m1,m2)+cuxxyyy020*rxssa(m1,m2)+cuxxyyy120*rxrssa(m1,m2)+cuxxyyy220*rxrrssa(m1,m2)+cuxxyyy320*rxrrrssa(m1,m2)+cuxxyyy030*rxsssa(m1,m2)+cuxxyyy130*rxrsssa(m1,m2)+cuxxyyy230*rxrrsssa(m1,m2)+cuxxyyy040*rxssssa(m1,m2)+cuxxyyy140*rxrssssa(m1,m2)+cuxxyyy050*rxsssssa(m1,m2)
                                                  rxxyyyya(m1,m2) = cuxyyyy100*rxr(m1,m2,0)+cuxyyyy200*rxrra(m1,m2)+cuxyyyy300*rxrrra(m1,m2)+cuxyyyy400*rxrrrra(m1,m2)+cuxyyyy500*rxrrrrra(m1,m2)+cuxyyyy010*rxr(m1,m2,1)+cuxyyyy110*rxrsa(m1,m2)+cuxyyyy210*rxrrsa(m1,m2)+cuxyyyy310*rxrrrsa(m1,m2)+cuxyyyy410*rxrrrrsa(m1,m2)+cuxyyyy020*rxssa(m1,m2)+cuxyyyy120*rxrssa(m1,m2)+cuxyyyy220*rxrrssa(m1,m2)+cuxyyyy320*rxrrrssa(m1,m2)+cuxyyyy030*rxsssa(m1,m2)+cuxyyyy130*rxrsssa(m1,m2)+cuxyyyy230*rxrrsssa(m1,m2)+cuxyyyy040*rxssssa(m1,m2)+cuxyyyy140*rxrssssa(m1,m2)+cuxyyyy050*rxsssssa(m1,m2)
                                                  rxyyyyya(m1,m2) = cuyyyyy100*rxr(m1,m2,0)+cuyyyyy200*rxrra(m1,m2)+cuyyyyy300*rxrrra(m1,m2)+cuyyyyy400*rxrrrra(m1,m2)+cuyyyyy500*rxrrrrra(m1,m2)+cuyyyyy010*rxr(m1,m2,1)+cuyyyyy110*rxrsa(m1,m2)+cuyyyyy210*rxrrsa(m1,m2)+cuyyyyy310*rxrrrsa(m1,m2)+cuyyyyy410*rxrrrrsa(m1,m2)+cuyyyyy020*rxssa(m1,m2)+cuyyyyy120*rxrssa(m1,m2)+cuyyyyy220*rxrrssa(m1,m2)+cuyyyyy320*rxrrrssa(m1,m2)+cuyyyyy030*rxsssa(m1,m2)+cuyyyyy130*rxrsssa(m1,m2)+cuyyyyy230*rxrrsssa(m1,m2)+cuyyyyy040*rxssssa(m1,m2)+cuyyyyy140*rxrssssa(m1,m2)+cuyyyyy050*rxsssssa(m1,m2)
                                              end do
                                              end do
                                              rxxxxxx=rxxxxxxa(0,0); ryxxxxx=rxxxxxxa(0,1); sxxxxxx=rxxxxxxa(1,0); syxxxxx=rxxxxxxa(1,1); 
                                              rxxxxxy=rxxxxxya(0,0); ryxxxxy=rxxxxxya(0,1); sxxxxxy=rxxxxxya(1,0); syxxxxy=rxxxxxya(1,1); 
                                              rxxxxyy=rxxxxyya(0,0); ryxxxyy=rxxxxyya(0,1); sxxxxyy=rxxxxyya(1,0); syxxxyy=rxxxxyya(1,1); 
                                              rxxxyyy=rxxxyyya(0,0); ryxxyyy=rxxxyyya(0,1); sxxxyyy=rxxxyyya(1,0); syxxyyy=rxxxyyya(1,1); 
                                              rxxyyyy=rxxyyyya(0,0); ryxyyyy=rxxyyyya(0,1); sxxyyyy=rxxyyyya(1,0); syxyyyy=rxxyyyya(1,1); 
                                              rxyyyyy=rxyyyyya(0,0); ryyyyyy=rxyyyyya(0,1); sxyyyyy=rxyyyyya(1,0); syyyyyy=rxyyyyya(1,1); 
                        ! if( i1.eq.5 .and. i2.eq.6 )then  
                        !    write(*,'(" (i1,i2)=(",2i3,") rxxxxxx,ryxxxxx,sxxxxxx,syxxxxx=",4(1pe12.4,1x))') i1,i2,rxxxxxx,ryxxxxx,sxxxxxx,syxxxxx
                        !    write(*,'(" (i1,i2)=(",2i3,") rxxxxxy,ryxxxxy,sxxxxxy,syxxxxy=",4(1pe12.4,1x))') i1,i2,rxxxxxy,ryxxxxy,sxxxxxy,syxxxxy
                        !    write(*,'(" (i1,i2)=(",2i3,") rxxxxyy,ryxxxyy,sxxxxyy,syxxxyy=",4(1pe12.4,1x))') i1,i2,rxxxxyy,ryxxxyy,sxxxxyy,syxxxyy
                        !    write(*,'(" (i1,i2)=(",2i3,") rxxxyyy,ryxxyyy,sxxxyyy,syxxyyy=",4(1pe12.4,1x))') i1,i2,rxxxyyy,ryxxyyy,sxxxyyy,syxxyyy
                        !    write(*,'(" (i1,i2)=(",2i3,") rxxyyyy,ryxyyyy,sxxyyyy,syxyyyy=",4(1pe12.4,1x))') i1,i2,rxxyyyy,ryxyyyy,sxxyyyy,syxyyyy
                        !    write(*,'(" (i1,i2)=(",2i3,") rxyyyyy,ryyyyyy,sxyyyyy,syyyyyy=",4(1pe12.4,1x))') i1,i2,rxyyyyy,ryyyyyy,sxyyyyy,syyyyyy
                        !  end if      
                       ! ---- SIXTH parametric derivs----
                       !  1/4 , 13/240 
                                              urrrrrr = d600i ! - (dx(0)**2/4.)*d800i
                                              ussssss = d060i ! - (dx(1)**2/4.)*d080i  
                       ! uxxxxxy = D0x(D+xD-x)^2 *( 1 - (1/3)*dx^2 D+xD-x + ... ) X 
                       !                      D0y*( 1 - dy^2/6 D+yD-y  + dy^4/30*(D+yD-y)^2 - (1/140)*dy^6*(D+xD-x)^3 ) u 
                                              urrrrrs = d510i ! - (dx(0)**2/3.)*d710i - (dx(1)**2/6.)*d530i
                                              ursssss = d150i ! - (dx(1)**2/3.)*d170i - (dx(0)**2/6.)*d350i
                       ! uxxxxyy = (D+xD-x)^2*( 1 - dx^2*(1/6)*D+xD-x + dx^4*(7/240)*(D+xD-x)^2 + ... )
                       !               D+yD-y*( 1 - dy^2/12 D+yD-y + dy^4/90*(D+yD-y)^2 - (1/560)*dy^6*(D+yD-y)^3 ) u 
                                              urrrrss = d420i ! - (dx(0)**2/6.)*d620i - (dx(1)**2/12.)*d440i
                                              urrssss = d240i ! - (dx(1)**2/6.)*d260i - (dx(0)**2/12.)*d440i
                       ! uxxxyyy = D0x D+xD-x*( 1 - (1/4)*dx^2 * D+xD-x + (7/120)*dx^4 (D+xD-x)^2 + ... ) X 
                       !           D0y D+yD-y*( 1 - (1/4)*dy^2 * D+yD-y + (7/120)*dy^4 (D+yD-y)^2 + ...  ) u 
                                              urrrsss = d330i ! - (dx(0)**2/4.)*d530i - (dx(1)**2/4.)*d350i
                       ! ----- SIXTH SPATIAL DERIVATIVES -----
                       ! ------ Coefficients in expansion for uxxxxxx ------
cuxxxxxx100 = rxxxxxx
cuxxxxxx200 = 6.*rx*rxxxxx+15.*rxx*rxxxx+10.*rxxx**2
cuxxxxxx300 = 15.*rx**2.*rxxxx+60.*rx*rxx*rxxx+15.*rxx**3
cuxxxxxx400 = 20.*rx**3.*rxxx+45.*rx**2.*rxx**2
cuxxxxxx500 = 15.*rx**4.*rxx
cuxxxxxx600 = rx**6
cuxxxxxx010 = sxxxxxx
cuxxxxxx110 = 6.*rx*sxxxxx+15.*rxx*sxxxx+20.*rxxx*sxxx+15.*rxxxx*sxx+6.*rxxxxx*sx
cuxxxxxx210 = 15.*rx**2.*sxxxx+(60.*rxx*sxxx+60.*rxxx*sxx+30.*rxxxx*sx)*rx+60.*rxxx*rxx*sx+45.*sxx*rxx**2
cuxxxxxx310 = 20.*rx**3.*sxxx+(90.*rxx*sxx+60.*rxxx*sx)*rx**2+90.*rxx**2.*sx*rx
cuxxxxxx410 = 15.*rx**4.*sxx+60.*rx**3.*rxx*sx
cuxxxxxx510 = 6.*rx**5.*sx
cuxxxxxx020 = 6.*sx*sxxxxx+15.*sxx*sxxxx+10.*sxxx**2
cuxxxxxx120 = 15.*rxxxx*sx**2+(30.*rx*sxxxx+60.*rxx*sxxx+60.*rxxx*sxx)*sx+60.*rx*sxx*sxxx+45.*sxx**2.*rxx
cuxxxxxx220 = (60.*rx*rxxx+45.*rxx**2)*sx**2+(60.*rx**2.*sxxx+180.*rx*rxx*sxx)*sx+45.*rx**2.*sxx**2
cuxxxxxx320 = 60.*sx*rx**2.*(rx*sxx+3./2.*rxx*sx)
cuxxxxxx420 = 15.*rx**4.*sx**2
cuxxxxxx030 = 15.*sx**2.*sxxxx+60.*sx*sxx*sxxx+15.*sxx**3
cuxxxxxx130 = 20.*rxxx*sx**3+(60.*rx*sxxx+90.*rxx*sxx)*sx**2+90.*rx*sx*sxx**2
cuxxxxxx230 = 90.*sx**2.*rx*(rx*sxx+2./3.*rxx*sx)
cuxxxxxx330 = 20.*rx**3.*sx**3
cuxxxxxx040 = 20.*sx**3.*sxxx+45.*sx**2.*sxx**2
cuxxxxxx140 = 60.*rx*sx**3.*sxx+15.*rxx*sx**4
cuxxxxxx240 = 15.*rx**2.*sx**4
cuxxxxxx050 = 15.*sxx*sx**4
cuxxxxxx150 = 6.*rx*sx**5
cuxxxxxx060 = sx**6
                       ! uxxxxxx = cuxxxxxx100*ur+cuxxxxxx200*urr+cuxxxxxx300*urrr+cuxxxxxx400*urrrr+cuxxxxxx500*urrrrr+cuxxxxxx600*urrrrrr+cuxxxxxx010*us+cuxxxxxx110*urs+cuxxxxxx210*urrs+cuxxxxxx310*urrrs+cuxxxxxx410*urrrrs+cuxxxxxx510*urrrrrs+cuxxxxxx020*uss+cuxxxxxx120*urss+cuxxxxxx220*urrss+cuxxxxxx320*urrrss+cuxxxxxx420*urrrrss+cuxxxxxx030*usss+cuxxxxxx130*ursss+cuxxxxxx230*urrsss+cuxxxxxx330*urrrsss+cuxxxxxx040*ussss+cuxxxxxx140*urssss+cuxxxxxx240*urrssss+cuxxxxxx050*usssss+cuxxxxxx150*ursssss+cuxxxxxx060*ussssss
                       ! ------ Coefficients in expansion for uxxxxxy ------
cuxxxxxy100 = rxxxxxy
cuxxxxxy200 = 5.*rx*rxxxxy+10.*rxx*rxxxy+10.*rxxx*rxxy+5.*rxxxx*rxy+rxxxxx*ry
cuxxxxxy300 = 10.*rxxxy*rx**2+(30.*rxx*rxxy+20.*rxxx*rxy+5.*rxxxx*ry)*rx+15.*rxy*rxx**2+10.*rxxx*ry*rxx
cuxxxxxy400 = 10.*rx**3.*rxxy+30.*rx**2.*rxx*rxy+10.*rx**2.*rxxx*ry+15.*rx*rxx**2.*ry
cuxxxxxy500 = 5.*rx**4.*rxy+10.*rx**3.*rxx*ry
cuxxxxxy600 = ry*rx**5
cuxxxxxy010 = sxxxxxy
cuxxxxxy110 = 5.*rx*sxxxxy+10.*rxx*sxxxy+10.*rxxx*sxxy+5.*rxxxx*sxy+rxxxxx*sy+5.*rxxxxy*sx+10.*rxxxy*sxx+10.*rxxy*sxxx+5.*rxy*sxxxx+ry*sxxxxx
cuxxxxxy210 = 10.*rx**2.*sxxxy+(30.*rxx*sxxy+20.*rxxx*sxy+5.*rxxxx*sy+20.*rxxxy*sx+30.*rxxy*sxx+20.*rxy*sxxx+5.*ry*sxxxx)*rx+15.*sxy*rxx**2+(10.*rxxx*sy+30.*rxxy*sx+30.*rxy*sxx+10.*ry*sxxx)*rxx+(20.*rxy*sx+10.*ry*sxx)*rxxx+5.*rxxxx*ry*sx
cuxxxxxy310 = 10.*rx**3.*sxxy+(30.*rxx*sxy+10.*rxxx*sy+30.*rxxy*sx+30.*rxy*sxx+10.*ry*sxxx)*rx**2+(15.*rxx**2.*sy+(60.*rxy*sx+30.*ry*sxx)*rxx+20.*ry*sx*rxxx)*rx+15.*rxx**2.*ry*sx
cuxxxxxy410 = 5.*rx**2.*(rx**2.*sxy+(2.*rxx*sy+4.*rxy*sx+2.*ry*sxx)*rx+6.*rxx*ry*sx)
cuxxxxxy510 = rx**5.*sy+5.*rx**4.*ry*sx
cuxxxxxy020 = 5.*sx*sxxxxy+10.*sxx*sxxxy+10.*sxxx*sxxy+5.*sxxxx*sxy+sxxxxx*sy
cuxxxxxy120 = 10.*rxxxy*sx**2+(20.*rx*sxxxy+30.*rxx*sxxy+20.*rxxx*sxy+5.*rxxxx*sy+30.*rxxy*sxx+20.*rxy*sxxx+5.*ry*sxxxx)*sx+15.*rxy*sxx**2+(30.*rx*sxxy+30.*rxx*sxy+10.*rxxx*sy+10.*ry*sxxx)*sxx+(20.*rx*sxy+10.*rxx*sy)*sxxx+5.*rx*sxxxx*sy
cuxxxxxy220 = (30.*sx*sxxy+30.*sxx*sxy+10.*sxxx*sy)*rx**2+(30.*rxxy*sx**2+(60.*rxx*sxy+20.*rxxx*sy+60.*rxy*sxx+20.*ry*sxxx)*sx+30.*rxx*sxx*sy+15.*sxx**2.*ry)*rx+10.*((3.*rxx*rxy+rxxx*ry)*sx+3.*rxx*sxx*ry+3/2.*rxx**2.*sy)*sx
cuxxxxxy320 = 20.*rx**3.*sx*sxy+10.*rx**3.*sxx*sy+30.*rx**2.*rxx*sx*sy+30.*rx**2.*rxy*sx**2+30.*rx**2.*ry*sx*sxx+30.*rx*rxx*ry*sx**2
cuxxxxxy420 = 5.*rx**4.*sx*sy+10.*rx**3.*ry*sx**2
cuxxxxxy030 = 10.*sx**2.*sxxxy+(30.*sxx*sxxy+20.*sxxx*sxy+5.*sxxxx*sy)*sx+15.*sxy*sxx**2+10.*sxxx*sy*sxx
cuxxxxxy130 = 10.*rxxy*sx**3+(30.*rx*sxxy+30.*rxx*sxy+10.*rxxx*sy+30.*rxy*sxx+10.*ry*sxxx)*sx**2+(15.*sxx**2.*ry+(60.*rx*sxy+30.*rxx*sy)*sxx+20.*rx*sxxx*sy)*sx+15.*rx*sxx**2.*sy
cuxxxxxy230 = 30.*rx**2.*sx**2.*sxy+30.*rx**2.*sx*sxx*sy+30.*rx*rxx*sx**2.*sy+20.*rx*rxy*sx**3+30.*rx*ry*sx**2.*sxx+10.*rxx*ry*sx**3
cuxxxxxy330 = 10.*rx**3.*sx**2.*sy+10.*rx**2.*ry*sx**3
cuxxxxxy040 = 10.*sx**3.*sxxy+30.*sx**2.*sxx*sxy+10.*sx**2.*sxxx*sy+15.*sx*sxx**2.*sy
cuxxxxxy140 = 5.*rxy*sx**4+(20.*rx*sxy+10.*rxx*sy+10.*ry*sxx)*sx**3+30.*rx*sx**2.*sxx*sy
cuxxxxxy240 = 10.*rx**2.*sx**3.*sy+5.*rx*ry*sx**4
cuxxxxxy050 = 5.*sx**4.*sxy+10.*sx**3.*sxx*sy
cuxxxxxy150 = 5.*rx*sx**4.*sy+ry*sx**5
cuxxxxxy060 = sy*sx**5
                       ! uxxxxxy = cuxxxxxy100*ur+cuxxxxxy200*urr+cuxxxxxy300*urrr+cuxxxxxy400*urrrr+cuxxxxxy500*urrrrr+cuxxxxxy600*urrrrrr+cuxxxxxy010*us+cuxxxxxy110*urs+cuxxxxxy210*urrs+cuxxxxxy310*urrrs+cuxxxxxy410*urrrrs+cuxxxxxy510*urrrrrs+cuxxxxxy020*uss+cuxxxxxy120*urss+cuxxxxxy220*urrss+cuxxxxxy320*urrrss+cuxxxxxy420*urrrrss+cuxxxxxy030*usss+cuxxxxxy130*ursss+cuxxxxxy230*urrsss+cuxxxxxy330*urrrsss+cuxxxxxy040*ussss+cuxxxxxy140*urssss+cuxxxxxy240*urrssss+cuxxxxxy050*usssss+cuxxxxxy150*ursssss+cuxxxxxy060*ussssss
                       ! ------ Coefficients in expansion for uxxxxyy ------
cuxxxxyy100 = rxxxxyy
cuxxxxyy200 = 4.*rx*rxxxyy+6.*rxx*rxxyy+4.*rxxx*rxyy+rxxxx*ryy+2.*rxxxxy*ry+8.*rxxxy*rxy+6.*rxxy**2
cuxxxxyy300 = 6.*rxxyy*rx**2+(12.*rxx*rxyy+4.*rxxx*ryy+8.*rxxxy*ry+24.*rxxy*rxy)*rx+3.*ryy*rxx**2+(12.*rxxy*ry+12.*rxy**2)*rxx+rxxxx*ry**2+8.*rxxx*rxy*ry
cuxxxxyy400 = (24.*rxx*rxy*ry+4.*rxxx*ry**2)*rx+(6.*rxx*ryy+12.*rxxy*ry+12.*rxy**2)*rx**2+4.*rxyy*rx**3+3.*rxx**2.*ry**2
cuxxxxyy500 = rx**4.*ryy+8.*rx**3.*rxy*ry+6.*rx**2.*rxx*ry**2
cuxxxxyy600 = ry**2.*rx**4
cuxxxxyy010 = sxxxxyy
cuxxxxyy110 = 4.*rx*sxxxyy+6.*rxx*sxxyy+4.*rxxx*sxyy+rxxxx*syy+2.*rxxxxy*sy+8.*rxxxy*sxy+4.*rxxxyy*sx+12.*rxxy*sxxy+6.*rxxyy*sxx+8.*rxy*sxxxy+4.*rxyy*sxxx+2.*ry*sxxxxy+ryy*sxxxx
cuxxxxyy210 = 6.*rx**2.*sxxyy+(12.*rxx*sxyy+4.*rxxx*syy+8.*rxxxy*sy+24.*rxxy*sxy+12.*rxxyy*sx+24.*rxy*sxxy+12.*rxyy*sxx+8.*ry*sxxxy+4.*ryy*sxxx)*rx+ry**2.*sxxxx+(12.*rxx*sxxy+8.*rxxx*sxy+2.*rxxxx*sy+8.*rxxxy*sx+12.*rxxy*sxx+8.*rxy*sxxx)*ry+3.*syy*rxx**2+(12.*rxxy*sy+24.*rxy*sxy+12.*rxyy*sx+6.*ryy*sxx)*rxx+12.*rxy**2.*sxx+(8.*rxxx*sy+24.*rxxy*sx)*rxy+4.*rxxx*ryy*sx
cuxxxxyy310 = 4.*rx**3.*sxyy+(6.*rxx*syy+12.*rxxy*sy+24.*rxy*sxy+12.*rxyy*sx+12.*ry*sxxy+6.*ryy*sxx)*rx**2+(4.*ry**2.*sxxx+(24.*rxx*sxy+8.*rxxx*sy+24.*rxxy*sx+24.*rxy*sxx)*ry+(24.*rxy*sy+12.*ryy*sx)*rxx+24.*sx*rxy**2)*rx+(6.*rxx*sxx+4.*rxxx*sx)*ry**2+(6.*rxx**2.*sy+24.*rxx*rxy*sx)*ry
cuxxxxyy410 = rx*(rx**3.*syy+(8.*rxy*sy+8.*ry*sxy+4.*ryy*sx)*rx**2+6.*ry*(2.*rxx*sy+4.*rxy*sx+ry*sxx)*rx+12.*rxx*ry**2.*sx)
cuxxxxyy510 = 2.*rx**4.*ry*sy+4.*rx**3.*ry**2.*sx
cuxxxxyy020 = 4.*sx*sxxxyy+6.*sxx*sxxyy+4.*sxxx*sxyy+sxxxx*syy+2.*sxxxxy*sy+8.*sxxxy*sxy+6.*sxxy**2
cuxxxxyy120 = 6.*rxxyy*sx**2+(12.*rx*sxxyy+12.*rxx*sxyy+4.*rxxx*syy+8.*rxxxy*sy+24.*rxxy*sxy+24.*rxy*sxxy+12.*rxyy*sxx+8.*ry*sxxxy+4.*ryy*sxxx)*sx+rxxxx*sy**2+(8.*rx*sxxxy+12.*rxx*sxxy+8.*rxxx*sxy+12.*rxxy*sxx+8.*rxy*sxxx+2.*ry*sxxxx)*sy+3.*ryy*sxx**2+(12.*rx*sxyy+6.*rxx*syy+24.*rxy*sxy+12.*ry*sxxy)*sxx+12.*rxx*sxy**2+(24.*rx*sxxy+8.*ry*sxxx)*sxy+4.*rx*sxxx*syy
cuxxxxyy220 = (12.*sx*sxyy+6.*sxx*syy+12.*sxxy*sy+12.*sxy**2)*rx**2+(12.*rxyy*sx**2+(12.*rxx*syy+24.*rxxy*sy+48.*rxy*sxy+24.*ry*sxxy+12.*ryy*sxx)*sx+(24.*sxx*sxy+8.*sxxx*sy)*ry+4.*sy*(6.*rxx*sxy+rxxx*sy+6.*rxy*sxx))*rx+(6.*rxx*ryy+12.*rxxy*ry+12.*rxy**2)*sx**2+(4.*ry**2.*sxxx+(24.*rxx*sxy+8.*rxxx*sy+24.*rxy*sxx)*ry+24.*rxy*rxx*sy)*sx+3.*sxx**2.*ry**2+12.*rxx*sxx*sy*ry+3.*rxx**2.*sy**2
cuxxxxyy320 = (4.*sx*syy+8.*sxy*sy)*rx**3+(6.*ryy*sx**2+(24.*rxy*sy+24.*ry*sxy)*sx+12.*sxx*sy*ry+6.*rxx*sy**2)*rx**2+12.*ry*sx*(2.*rxx*sy+2.*rxy*sx+ry*sxx)*rx+6.*rxx*ry**2.*sx**2
cuxxxxyy420 = rx**4.*sy**2+8.*rx**3.*ry*sx*sy+6.*rx**2.*ry**2.*sx**2
cuxxxxyy030 = 6.*sx**2.*sxxyy+(12.*sxx*sxyy+4.*sxxx*syy+8.*sxxxy*sy+24.*sxxy*sxy)*sx+3.*syy*sxx**2+(12.*sxxy*sy+12.*sxy**2)*sxx+sxxxx*sy**2+8.*sxxx*sxy*sy
cuxxxxyy130 = 4.*rxyy*sx**3+(12.*rx*sxyy+6.*rxx*syy+12.*rxxy*sy+24.*rxy*sxy+12.*ry*sxxy+6.*ryy*sxx)*sx**2+(4.*rxxx*sy**2+(24.*rx*sxxy+24.*rxx*sxy+24.*rxy*sxx+8.*ry*sxxx)*sy+(12.*rx*syy+24.*ry*sxy)*sxx+24.*rx*sxy**2)*sx+(4.*rx*sxxx+6.*rxx*sxx)*sy**2+(24.*rx*sxx*sxy+6.*ry*sxx**2)*sy
cuxxxxyy230 = (4.*rx*ryy+8.*rxy*ry)*sx**3+(6.*rx**2.*syy+(24.*rxy*sy+24.*ry*sxy)*rx+6.*ry**2.*sxx+12.*ry*sy*rxx)*sx**2+24.*sy*(rx*sxy+ry*sxx+1/2.*rxx*sy)*rx*sx+6.*rx**2.*sxx*sy**2
cuxxxxyy330 = 4.*rx**3.*sx*sy**2+12.*rx**2.*ry*sx**2.*sy+4.*rx*ry**2.*sx**3
cuxxxxyy040 = 4.*sx**3.*sxyy+(6.*sxx*syy+12.*sxxy*sy+12.*sxy**2)*sx**2+3.*sxx**2.*sy**2+(24.*sxx*sxy*sy+4.*sxxx*sy**2)*sx
cuxxxxyy140 = 4.*(1/4.*ryy*sx**3+(rx*syy+2.*rxy*sy+2.*ry*sxy)*sx**2+6.*(rx*sxy+1/2.*ry*sxx+1/4.*rxx*sy)*sy*sx+3.*rx*sy**2.*sxx)*sx
cuxxxxyy240 = 6.*rx**2.*sx**2.*sy**2+8.*rx*ry*sx**3.*sy+ry**2.*sx**4
cuxxxxyy050 = sx**4.*syy+8.*sx**3.*sxy*sy+6.*sx**2.*sxx*sy**2
cuxxxxyy150 = 4.*rx*sx**3.*sy**2+2.*ry*sx**4.*sy
cuxxxxyy060 = sy**2.*sx**4
                       ! uxxxxyy = cuxxxxyy100*ur+cuxxxxyy200*urr+cuxxxxyy300*urrr+cuxxxxyy400*urrrr+cuxxxxyy500*urrrrr+cuxxxxyy600*urrrrrr+cuxxxxyy010*us+cuxxxxyy110*urs+cuxxxxyy210*urrs+cuxxxxyy310*urrrs+cuxxxxyy410*urrrrs+cuxxxxyy510*urrrrrs+cuxxxxyy020*uss+cuxxxxyy120*urss+cuxxxxyy220*urrss+cuxxxxyy320*urrrss+cuxxxxyy420*urrrrss+cuxxxxyy030*usss+cuxxxxyy130*ursss+cuxxxxyy230*urrsss+cuxxxxyy330*urrrsss+cuxxxxyy040*ussss+cuxxxxyy140*urssss+cuxxxxyy240*urrssss+cuxxxxyy050*usssss+cuxxxxyy150*ursssss+cuxxxxyy060*ussssss
                       ! ------ Coefficients in expansion for uxxxyyy ------
cuxxxyyy100 = rxxxyyy
cuxxxyyy200 = 3.*rx*rxxyyy+3.*rxx*rxyyy+rxxx*ryyy+3.*rxxxy*ryy+3.*rxxxyy*ry+9.*rxxy*rxyy+9.*rxxyy*rxy
cuxxxyyy300 = (9.*rxx*rxyy+3.*rxxx*ryy+18.*rxxy*rxy)*ry+(3.*rxx*ryyy+9.*rxxy*ryy+9.*rxxyy*ry+18.*rxy*rxyy)*rx+3.*rxxxy*ry**2+3.*rxyyy*rx**2+9.*rxx*ryy*rxy+6.*rxy**3
cuxxxyyy400 = ryyy*rx**3+(9.*rxy*ryy+9.*rxyy*ry)*rx**2+9.*ry*(rxx*ryy+rxxy*ry+2.*rxy**2)*rx+rxxx*ry**3+9.*rxy*rxx*ry**2
cuxxxyyy500 = 3.*rx**3.*ry*ryy+9.*rx**2.*rxy*ry**2+3.*rx*rxx*ry**3
cuxxxyyy600 = rx**3.*ry**3
cuxxxyyy010 = sxxxyyy
cuxxxyyy110 = 3.*rx*sxxyyy+3.*rxx*sxyyy+rxxx*syyy+3.*rxxxy*syy+3.*rxxxyy*sy+9.*rxxy*sxyy+9.*rxxyy*sxy+3.*rxxyyy*sx+9.*rxy*sxxyy+9.*rxyy*sxxy+3.*rxyyy*sxx+3.*ry*sxxxyy+3.*ryy*sxxxy+ryyy*sxxx
cuxxxyyy210 = 3.*rx**2.*sxyyy+(3.*rxx*syyy+9.*rxxy*syy+9.*rxxyy*sy+18.*rxy*sxyy+18.*rxyy*sxy+6.*rxyyy*sx+9.*ry*sxxyy+9.*ryy*sxxy+3.*ryyy*sxx)*rx+3.*ry**2.*sxxxy+(9.*rxx*sxyy+3.*rxxx*syy+6.*rxxxy*sy+18.*rxxy*sxy+9.*rxxyy*sx+18.*rxy*sxxy+9.*rxyy*sxx+3.*ryy*sxxx)*ry+18.*rxy**2.*sxy+(9.*rxx*syy+18.*rxxy*sy+18.*rxyy*sx+9.*ryy*sxx)*rxy+(9.*rxyy*sy+9.*ryy*sxy+3.*ryyy*sx)*rxx+3.*ryy*(rxxx*sy+3.*rxxy*sx)
cuxxxyyy310 = rx**3.*syyy+(9.*rxy*syy+9.*rxyy*sy+9.*ry*sxyy+9.*ryy*sxy+3.*ryyy*sx)*rx**2+(9.*ry**2.*sxxy+(9.*rxx*syy+18.*rxxy*sy+36.*rxy*sxy+18.*rxyy*sx+9.*ryy*sxx)*ry+18.*rxy*ryy*sx+9.*sy*rxx*ryy+18.*sy*rxy**2)*rx+3.*ry*(1/3.*ry**2.*sxxx+(3.*rxx*sxy+rxxx*sy+3.*rxxy*sx+3.*rxy*sxx)*ry+3.*sx*rxx*ryy+6.*sx*rxy**2+6.*rxy*rxx*sy)
cuxxxyyy410 = (3.*ry*syy+3.*ryy*sy)*rx**3+9.*ry*(2.*rxy*sy+ry*sxy+ryy*sx)*rx**2+3.*ry**2.*(3.*rxx*sy+6.*rxy*sx+ry*sxx)*rx+3.*rxx*sx*ry**3
cuxxxyyy510 = 3.*rx**3.*ry**2.*sy+3.*rx**2.*ry**3.*sx
cuxxxyyy020 = 3.*sx*sxxyyy+3.*sxx*sxyyy+sxxx*syyy+3.*sxxxy*syy+3.*sxxxyy*sy+9.*sxxy*sxyy+9.*sxxyy*sxy
cuxxxyyy120 = 3.*rxyyy*sx**2+(6.*rx*sxyyy+3.*rxx*syyy+9.*rxxy*syy+9.*rxxyy*sy+18.*rxy*sxyy+18.*rxyy*sxy+9.*ry*sxxyy+9.*ryy*sxxy+3.*ryyy*sxx)*sx+3.*rxxxy*sy**2+(9.*rx*sxxyy+9.*rxx*sxyy+3.*rxxx*syy+18.*rxxy*sxy+18.*rxy*sxxy+9.*rxyy*sxx+6.*ry*sxxxy+3.*ryy*sxxx)*sy+18.*sxy**2.*rxy+(18.*rx*sxyy+9.*rxx*syy+18.*ry*sxxy+9.*ryy*sxx)*sxy+(3.*rx*syyy+9.*rxy*syy+9.*ry*sxyy)*sxx+(9.*rx*sxxy+3.*ry*sxxx)*syy
cuxxxyyy220 = (3.*sx*syyy+9.*sxy*syy+9.*sxyy*sy)*rx**2+((18.*sx*sxyy+9.*sxx*syy+18.*sxxy*sy+18.*sxy**2)*ry+3.*ryyy*sx**2+(18.*rxy*syy+18.*rxyy*sy+18.*ryy*sxy)*sx+9.*sy*(rxx*syy+rxxy*sy+4.*rxy*sxy+ryy*sxx))*rx+(9.*sx*sxxy+9.*sxx*sxy+3.*sxxx*sy)*ry**2+(9.*rxyy*sx**2+(9.*rxx*syy+18.*rxxy*sy+36.*rxy*sxy+9.*ryy*sxx)*sx+3.*sy*(6.*rxx*sxy+rxxx*sy+6.*rxy*sxx))*ry+9.*rxy*ryy*sx**2+9.*sy*(rxx*ryy+2.*rxy**2)*sx+9.*rxy*rxx*sy**2
cuxxxyyy320 = 3.*sy*syy*rx**3+((9.*sx*syy+18.*sxy*sy)*ry+9.*ryy*sy*sx+9.*rxy*sy**2)*rx**2+18.*ry*((sxy*sx+1/2.*sy*sxx)*ry+1/2.*ryy*sx**2+2.*rxy*sy*sx+1/2.*rxx*sy**2)*rx+3.*ry**2.*sx*(3.*rxx*sy+3.*rxy*sx+ry*sxx)
cuxxxyyy420 = 3.*rx**3.*ry*sy**2+9.*rx**2.*ry**2.*sx*sy+3.*rx*ry**3.*sx**2
cuxxxyyy030 = 3.*sxxxy*sy**2+(3.*sxx*syyy+9.*sxxy*syy+9.*sxxyy*sy+18.*sxy*sxyy)*sx+3.*sx**2.*sxyyy+9.*sxx*syy*sxy+(9.*sxx*sxyy+3.*sxxx*syy+18.*sxxy*sxy)*sy+6.*sxy**3
cuxxxyyy130 = ryyy*sx**3+(3.*rx*syyy+9.*rxy*syy+9.*rxyy*sy+9.*ry*sxyy+9.*ryy*sxy)*sx**2+(9.*rxxy*sy**2+(18.*rx*sxyy+9.*rxx*syy+36.*rxy*sxy+18.*ry*sxxy+9.*ryy*sxx)*sy+18.*rx*syy*sxy+9.*syy*ry*sxx+18.*ry*sxy**2)*sx+9.*(1/9.*rxxx*sy**2+(rx*sxxy+1/3.*ry*sxxx+rxx*sxy+rxy*sxx)*sy+rx*syy*sxx+2.*rx*sxy**2+2.*ry*sxx*sxy)*sy
cuxxxyyy230 = 3.*ryy*ry*sx**3+((9.*rx*ryy+18.*rxy*ry)*sy+9.*rx*syy*ry+9.*ry**2.*sxy)*sx**2+9.*sy*((2.*rx*rxy+rxx*ry)*sy+rx**2.*syy+4.*ry*sxy*rx+ry**2.*sxx)*sx+9.*(rx*sxy+ry*sxx+1/3.*rxx*sy)*sy**2.*rx
cuxxxyyy330 = rx**3.*sy**3+9.*rx**2.*ry*sx*sy**2+9.*rx*ry**2.*sx**2.*sy+ry**3.*sx**3
cuxxxyyy040 = sx**3.*syyy+(9.*sxy*syy+9.*sxyy*sy)*sx**2+9.*sy*(sxx*syy+sxxy*sy+2.*sxy**2)*sx+sxxx*sy**3+9.*sxx*sxy*sy**2
cuxxxyyy140 = (3.*ry*syy+3.*ryy*sy)*sx**3+9.*sy*(rx*syy+rxy*sy+2.*ry*sxy)*sx**2+(3.*rxx*sy**3+(18.*rx*sxy+9.*ry*sxx)*sy**2)*sx+3.*rx*sxx*sy**3
cuxxxyyy240 = 3.*rx**2.*sx*sy**3+9.*rx*ry*sx**2.*sy**2+3.*ry**2.*sx**3.*sy
cuxxxyyy050 = 3.*sx**3.*sy*syy+9.*sx**2.*sxy*sy**2+3.*sx*sxx*sy**3
cuxxxyyy150 = 3.*rx*sx**2.*sy**3+3.*ry*sx**3.*sy**2
cuxxxyyy060 = sx**3.*sy**3
                       ! uxxxyyy = cuxxxyyy100*ur+cuxxxyyy200*urr+cuxxxyyy300*urrr+cuxxxyyy400*urrrr+cuxxxyyy500*urrrrr+cuxxxyyy600*urrrrrr+cuxxxyyy010*us+cuxxxyyy110*urs+cuxxxyyy210*urrs+cuxxxyyy310*urrrs+cuxxxyyy410*urrrrs+cuxxxyyy510*urrrrrs+cuxxxyyy020*uss+cuxxxyyy120*urss+cuxxxyyy220*urrss+cuxxxyyy320*urrrss+cuxxxyyy420*urrrrss+cuxxxyyy030*usss+cuxxxyyy130*ursss+cuxxxyyy230*urrsss+cuxxxyyy330*urrrsss+cuxxxyyy040*ussss+cuxxxyyy140*urssss+cuxxxyyy240*urrssss+cuxxxyyy050*usssss+cuxxxyyy150*ursssss+cuxxxyyy060*ussssss
                       ! ------ Coefficients in expansion for uxxyyyy ------
cuxxyyyy100 = rxxyyyy
cuxxyyyy200 = 2.*rx*rxyyyy+rxx*ryyyy+4.*rxxy*ryyy+6.*rxxyy*ryy+4.*rxxyyy*ry+8.*rxy*rxyyy+6.*rxyy**2
cuxxyyyy300 = 6.*rxxyy*ry**2+(8.*rx*rxyyy+4.*rxx*ryyy+12.*rxxy*ryy+24.*rxy*rxyy)*ry+3.*rxx*ryy**2+(12.*rx*rxyy+12.*rxy**2)*ryy+ryyyy*rx**2+8.*rx*rxy*ryyy
cuxxyyyy400 = (12.*rx*rxyy+6.*rxx*ryy+12.*rxy**2)*ry**2+(4.*rx**2.*ryyy+24.*rx*rxy*ryy)*ry+4.*rxxy*ry**3+3.*ryy**2.*rx**2
cuxxyyyy500 = 6.*rx**2.*ry**2.*ryy+8.*rx*rxy*ry**3+rxx*ry**4
cuxxyyyy600 = rx**2.*ry**4
cuxxyyyy010 = sxxyyyy
cuxxyyyy110 = 2.*rx*sxyyyy+rxx*syyyy+4.*rxxy*syyy+6.*rxxyy*syy+4.*rxxyyy*sy+8.*rxy*sxyyy+12.*rxyy*sxyy+8.*rxyyy*sxy+2.*rxyyyy*sx+4.*ry*sxxyyy+6.*ryy*sxxyy+4.*ryyy*sxxy+ryyyy*sxx
cuxxyyyy210 = 6.*ry**2.*sxxyy+(8.*rx*sxyyy+4.*rxx*syyy+12.*rxxy*syy+12.*rxxyy*sy+24.*rxy*sxyy+24.*rxyy*sxy+8.*rxyyy*sx+12.*ryy*sxxy+4.*ryyy*sxx)*ry+rx**2.*syyyy+(8.*rxy*syyy+12.*rxyy*syy+8.*rxyyy*sy+12.*ryy*sxyy+8.*ryyy*sxy+2.*ryyyy*sx)*rx+3.*ryy**2.*sxx+(6.*rxx*syy+12.*rxxy*sy+24.*rxy*sxy+12.*rxyy*sx)*ryy+12.*rxy**2.*syy+(24.*rxyy*sy+8.*ryyy*sx)*rxy+4.*rxx*ryyy*sy
cuxxyyyy310 = 4.*ry**3.*sxxy+(12.*rx*sxyy+6.*rxx*syy+12.*rxxy*sy+24.*rxy*sxy+12.*rxyy*sx+6.*ryy*sxx)*ry**2+(4.*rx**2.*syyy+(24.*rxy*syy+24.*rxyy*sy+24.*ryy*sxy+8.*ryyy*sx)*rx+(12.*rxx*sy+24.*rxy*sx)*ryy+24.*sy*rxy**2)*ry+6.*((syy*ryy+2/3.*ryyy*sy)*rx+ryy**2.*sx+4.*rxy*ryy*sy)*rx
cuxxyyyy410 = 6.*rx**2.*ry**2.*syy+12.*rx**2.*ry*ryy*sy+24.*rx*rxy*ry**2.*sy+8.*rx*ry**3.*sxy+12.*rx*ry**2.*ryy*sx+4.*rxx*ry**3.*sy+8.*rxy*ry**3.*sx+ry**4.*sxx
cuxxyyyy510 = 4.*rx**2.*ry**3.*sy+2.*rx*ry**4.*sx
cuxxyyyy020 = 2.*sx*sxyyyy+sxx*syyyy+4.*sxxy*syyy+6.*sxxyy*syy+4.*sxxyyy*sy+8.*sxy*sxyyy+6.*sxyy**2
cuxxyyyy120 = 6.*rxxyy*sy**2+(8.*rx*sxyyy+4.*rxx*syyy+12.*rxxy*syy+24.*rxy*sxyy+24.*rxyy*sxy+8.*rxyyy*sx+12.*ry*sxxyy+12.*ryy*sxxy+4.*ryyy*sxx)*sy+ryyyy*sx**2+(2.*rx*syyyy+8.*rxy*syyy+12.*rxyy*syy+8.*ry*sxyyy+12.*ryy*sxyy+8.*ryyy*sxy)*sx+3.*rxx*syy**2+(12.*rx*sxyy+24.*rxy*sxy+12.*ry*sxxy+6.*ryy*sxx)*syy+12.*ryy*sxy**2+(8.*rx*syyy+24.*ry*sxyy)*sxy+4.*ry*sxx*syyy
cuxxyyyy220 = (12.*sx*sxyy+6.*sxx*syy+12.*sxxy*sy+12.*sxy**2)*ry**2+(12.*rxxy*sy**2+(24.*rx*sxyy+12.*rxx*syy+48.*rxy*sxy+24.*rxyy*sx+12.*ryy*sxx)*sy+(8.*sx*syyy+24.*sxy*syy)*rx+24.*(rxy*syy+1/6.*ryyy*sx+ryy*sxy)*sx)*ry+(12.*rx*rxyy+6.*rxx*ryy+12.*rxy**2)*sy**2+(4.*rx**2.*syyy+(24.*rxy*syy+24.*ryy*sxy+8.*ryyy*sx)*rx+24.*rxy*ryy*sx)*sy+3.*syy**2.*rx**2+12.*ryy*syy*sx*rx+3.*ryy**2.*sx**2
cuxxyyyy320 = (8.*sx*sxy+4.*sxx*sy)*ry**3+(6.*rxx*sy**2+(24.*rx*sxy+24.*rxy*sx)*sy+12.*sx*syy*rx+6.*ryy*sx**2)*ry**2+12.*rx*sy*(rx*syy+2.*rxy*sy+2.*ryy*sx)*ry+6.*ryy*rx**2.*sy**2
cuxxyyyy420 = 6.*rx**2.*ry**2.*sy**2+8.*rx*ry**3.*sx*sy+ry**4.*sx**2
cuxxyyyy030 = 6.*sxxyy*sy**2+(8.*sx*sxyyy+4.*sxx*syyy+12.*sxxy*syy+24.*sxy*sxyy)*sy+3.*sxx*syy**2+(12.*sx*sxyy+12.*sxy**2)*syy+sx**2.*syyyy+8.*sx*sxy*syyy
cuxxyyyy130 = 4.*rxxy*sy**3+(12.*rx*sxyy+6.*rxx*syy+24.*rxy*sxy+12.*rxyy*sx+12.*ry*sxxy+6.*ryy*sxx)*sy**2+(4.*ryyy*sx**2+(8.*rx*syyy+24.*rxy*syy+24.*ry*sxyy+24.*ryy*sxy)*sx+(24.*rx*sxy+12.*ry*sxx)*syy+24.*ry*sxy**2)*sy+(4.*ry*syyy+6.*ryy*syy)*sx**2+(6.*rx*syy**2+24.*ry*sxy*syy)*sx
cuxxyyyy230 = (8.*rx*rxy+4.*rxx*ry)*sy**3+(6.*ry**2.*sxx+(24.*rx*sxy+24.*rxy*sx)*ry+6.*rx**2.*syy+12.*rx*sx*ryy)*sy**2+24.*ry*sx*(rx*syy+ry*sxy+1/2.*ryy*sx)*sy+6.*ry**2.*sx**2.*syy
cuxxyyyy330 = 4.*rx**2.*ry*sy**3+12.*rx*ry**2.*sx*sy**2+4.*ry**3.*sx**2.*sy
cuxxyyyy040 = (12.*sx*sxyy+6.*sxx*syy+12.*sxy**2)*sy**2+4.*sxxy*sy**3+3.*syy**2.*sx**2+(4.*sx**2.*syyy+24.*sx*sxy*syy)*sy
cuxxyyyy140 = 12.*rx*sx*sy**2.*syy+8.*rx*sxy*sy**3+rxx*sy**4+8.*rxy*sx*sy**3+12.*ry*sx**2.*sy*syy+24.*ry*sx*sxy*sy**2+4.*ry*sxx*sy**3+6.*ryy*sx**2.*sy**2
cuxxyyyy240 = rx**2.*sy**4+8.*rx*ry*sx*sy**3+6.*ry**2.*sx**2.*sy**2
cuxxyyyy050 = 6.*sx**2.*sy**2.*syy+8.*sx*sxy*sy**3+sxx*sy**4
cuxxyyyy150 = 2.*rx*sx*sy**4+4.*ry*sx**2.*sy**3
cuxxyyyy060 = sx**2.*sy**4
                       ! uxxyyyy = cuxxyyyy100*ur+cuxxyyyy200*urr+cuxxyyyy300*urrr+cuxxyyyy400*urrrr+cuxxyyyy500*urrrrr+cuxxyyyy600*urrrrrr+cuxxyyyy010*us+cuxxyyyy110*urs+cuxxyyyy210*urrs+cuxxyyyy310*urrrs+cuxxyyyy410*urrrrs+cuxxyyyy510*urrrrrs+cuxxyyyy020*uss+cuxxyyyy120*urss+cuxxyyyy220*urrss+cuxxyyyy320*urrrss+cuxxyyyy420*urrrrss+cuxxyyyy030*usss+cuxxyyyy130*ursss+cuxxyyyy230*urrsss+cuxxyyyy330*urrrsss+cuxxyyyy040*ussss+cuxxyyyy140*urssss+cuxxyyyy240*urrssss+cuxxyyyy050*usssss+cuxxyyyy150*ursssss+cuxxyyyy060*ussssss
                       ! ------ Coefficients in expansion for uxyyyyy ------
cuxyyyyy100 = rxyyyyy
cuxyyyyy200 = rx*ryyyyy+5.*rxy*ryyyy+10.*rxyy*ryyy+10.*rxyyy*ryy+5.*rxyyyy*ry
cuxyyyyy300 = 10.*rxyyy*ry**2+(5.*rx*ryyyy+20.*rxy*ryyy+30.*rxyy*ryy)*ry+15.*rxy*ryy**2+10.*ryyy*rx*ryy
cuxyyyyy400 = 10.*rx*ry**2.*ryyy+15.*rx*ry*ryy**2+30.*rxy*ry**2.*ryy+10.*rxyy*ry**3
cuxyyyyy500 = 10.*rx*ry**3.*ryy+5.*rxy*ry**4
cuxyyyyy600 = rx*ry**5
cuxyyyyy010 = sxyyyyy
cuxyyyyy110 = rx*syyyyy+5.*rxy*syyyy+10.*rxyy*syyy+10.*rxyyy*syy+5.*rxyyyy*sy+5.*ry*sxyyyy+10.*ryy*sxyyy+10.*ryyy*sxyy+5.*ryyyy*sxy+ryyyyy*sx
cuxyyyyy210 = 10.*ry**2.*sxyyy+(5.*rx*syyyy+20.*rxy*syyy+30.*rxyy*syy+20.*rxyyy*sy+30.*ryy*sxyy+20.*ryyy*sxy+5.*ryyyy*sx)*ry+15.*ryy**2.*sxy+(10.*rx*syyy+30.*rxy*syy+30.*rxyy*sy+10.*ryyy*sx)*ryy+(10.*rx*syy+20.*rxy*sy)*ryyy+5.*rx*ryyyy*sy
cuxyyyyy310 = 10.*ry**3.*sxyy+(10.*rx*syyy+30.*rxy*syy+30.*rxyy*sy+30.*ryy*sxy+10.*ryyy*sx)*ry**2+(15.*ryy**2.*sx+(30.*rx*syy+60.*rxy*sy)*ryy+20.*rx*ryyy*sy)*ry+15.*ryy**2.*rx*sy
cuxyyyyy410 = 10.*rx*ry**3.*syy+30.*rx*ry**2.*ryy*sy+20.*rxy*ry**3.*sy+5.*ry**4.*sxy+10.*ry**3.*ryy*sx
cuxyyyyy510 = 5.*rx*ry**4.*sy+ry**5.*sx
cuxyyyyy020 = sx*syyyyy+5.*sxy*syyyy+10.*sxyy*syyy+10.*sxyyy*syy+5.*sxyyyy*sy
cuxyyyyy120 = 10.*rxyyy*sy**2+(5.*rx*syyyy+20.*rxy*syyy+30.*rxyy*syy+20.*ry*sxyyy+30.*ryy*sxyy+20.*ryyy*sxy+5.*ryyyy*sx)*sy+15.*rxy*syy**2+(10.*rx*syyy+30.*ry*sxyy+30.*ryy*sxy+10.*ryyy*sx)*syy+(20.*ry*sxy+10.*ryy*sx)*syyy+5.*ry*sx*syyyy
cuxyyyyy220 = (10.*sx*syyy+30.*sxy*syy+30.*sxyy*sy)*ry**2+(30.*rxyy*sy**2+(20.*rx*syyy+60.*rxy*syy+60.*ryy*sxy+20.*ryyy*sx)*sy+30.*ryy*syy*sx+15.*syy**2.*rx)*ry+10.*sy*((rx*ryyy+3.*rxy*ryy)*sy+3.*ryy*syy*rx+3/2.*ryy**2.*sx)
cuxyyyyy320 = 30.*rx*ry**2.*sy*syy+30.*rx*ry*ryy*sy**2+30.*rxy*ry**2.*sy**2+10.*ry**3.*sx*syy+20.*ry**3.*sxy*sy+30.*ry**2.*ryy*sx*sy
cuxyyyyy420 = 10.*rx*ry**3.*sy**2+5.*ry**4.*sx*sy
cuxyyyyy030 = 10.*sxyyy*sy**2+(5.*sx*syyyy+20.*sxy*syyy+30.*sxyy*syy)*sy+15.*sxy*syy**2+10.*sx*syyy*syy
cuxyyyyy130 = 10.*rxyy*sy**3+(10.*rx*syyy+30.*rxy*syy+30.*ry*sxyy+30.*ryy*sxy+10.*ryyy*sx)*sy**2+(15.*syy**2.*rx+(60.*ry*sxy+30.*ryy*sx)*syy+20.*syyy*ry*sx)*sy+15.*ry*sx*syy**2
cuxyyyyy230 = 30.*rx*ry*sy**2.*syy+10.*rx*ryy*sy**3+20.*rxy*ry*sy**3+30.*ry**2.*sx*sy*syy+30.*ry**2.*sxy*sy**2+30.*ry*ryy*sx*sy**2
cuxyyyyy330 = 10.*rx*ry**2.*sy**3+10.*ry**3.*sx*sy**2
cuxyyyyy040 = 10.*sx*sy**2.*syyy+15.*sx*sy*syy**2+30.*sxy*sy**2.*syy+10.*sxyy*sy**3
cuxyyyyy140 = 10.*rx*sy**3.*syy+5.*rxy*sy**4+30.*ry*sx*sy**2.*syy+20.*ry*sxy*sy**3+10.*ryy*sx*sy**3
cuxyyyyy240 = 5.*rx*ry*sy**4+10.*ry**2.*sx*sy**3
cuxyyyyy050 = 10.*sx*sy**3.*syy+5.*sxy*sy**4
cuxyyyyy150 = rx*sy**5+5.*ry*sx*sy**4
cuxyyyyy060 = sx*sy**5
                       ! uxyyyyy = cuxyyyyy100*ur+cuxyyyyy200*urr+cuxyyyyy300*urrr+cuxyyyyy400*urrrr+cuxyyyyy500*urrrrr+cuxyyyyy600*urrrrrr+cuxyyyyy010*us+cuxyyyyy110*urs+cuxyyyyy210*urrs+cuxyyyyy310*urrrs+cuxyyyyy410*urrrrs+cuxyyyyy510*urrrrrs+cuxyyyyy020*uss+cuxyyyyy120*urss+cuxyyyyy220*urrss+cuxyyyyy320*urrrss+cuxyyyyy420*urrrrss+cuxyyyyy030*usss+cuxyyyyy130*ursss+cuxyyyyy230*urrsss+cuxyyyyy330*urrrsss+cuxyyyyy040*ussss+cuxyyyyy140*urssss+cuxyyyyy240*urrssss+cuxyyyyy050*usssss+cuxyyyyy150*ursssss+cuxyyyyy060*ussssss
                       ! ------ Coefficients in expansion for uyyyyyy ------
cuyyyyyy100 = ryyyyyy
cuyyyyyy200 = 6.*ry*ryyyyy+15.*ryy*ryyyy+10.*ryyy**2
cuyyyyyy300 = 15.*ry**2.*ryyyy+60.*ry*ryy*ryyy+15.*ryy**3
cuyyyyyy400 = 20.*ry**3.*ryyy+45.*ry**2.*ryy**2
cuyyyyyy500 = 15.*ry**4.*ryy
cuyyyyyy600 = ry**6
cuyyyyyy010 = syyyyyy
cuyyyyyy110 = 6.*ry*syyyyy+15.*ryy*syyyy+20.*ryyy*syyy+15.*ryyyy*syy+6.*ryyyyy*sy
cuyyyyyy210 = 15.*ry**2.*syyyy+(60.*ryy*syyy+60.*ryyy*syy+30.*ryyyy*sy)*ry+60.*ryyy*ryy*sy+45.*syy*ryy**2
cuyyyyyy310 = 20.*ry**3.*syyy+(90.*ryy*syy+60.*ryyy*sy)*ry**2+90.*ryy**2.*sy*ry
cuyyyyyy410 = 15.*ry**4.*syy+60.*ry**3.*ryy*sy
cuyyyyyy510 = 6.*sy*ry**5
cuyyyyyy020 = 6.*sy*syyyyy+15.*syy*syyyy+10.*syyy**2
cuyyyyyy120 = 15.*ryyyy*sy**2+(30.*ry*syyyy+60.*ryy*syyy+60.*ryyy*syy)*sy+60.*ry*syy*syyy+45.*syy**2.*ryy
cuyyyyyy220 = (60.*ry*ryyy+45.*ryy**2)*sy**2+(60.*ry**2.*syyy+180.*ry*ryy*syy)*sy+45.*ry**2.*syy**2
cuyyyyyy320 = 60.*sy*(ry*syy+3/2.*ryy*sy)*ry**2
cuyyyyyy420 = 15.*ry**4.*sy**2
cuyyyyyy030 = 15.*sy**2.*syyyy+60.*sy*syy*syyy+15.*syy**3
cuyyyyyy130 = 20.*ryyy*sy**3+(60.*ry*syyy+90.*ryy*syy)*sy**2+90.*ry*sy*syy**2
cuyyyyyy230 = 90.*sy**2.*ry*(ry*syy+2/3.*ryy*sy)
cuyyyyyy330 = 20.*ry**3.*sy**3
cuyyyyyy040 = 20.*sy**3.*syyy+45.*sy**2.*syy**2
cuyyyyyy140 = 60.*ry*sy**3.*syy+15.*ryy*sy**4
cuyyyyyy240 = 15.*ry**2.*sy**4
cuyyyyyy050 = 15.*sy**4.*syy
cuyyyyyy150 = 6.*ry*sy**5
cuyyyyyy060 = sy**6
                       ! uyyyyyy = cuyyyyyy100*ur+cuyyyyyy200*urr+cuyyyyyy300*urrr+cuyyyyyy400*urrrr+cuyyyyyy500*urrrrr+cuyyyyyy600*urrrrrr+cuyyyyyy010*us+cuyyyyyy110*urs+cuyyyyyy210*urrs+cuyyyyyy310*urrrs+cuyyyyyy410*urrrrs+cuyyyyyy510*urrrrrs+cuyyyyyy020*uss+cuyyyyyy120*urss+cuyyyyyy220*urrss+cuyyyyyy320*urrrss+cuyyyyyy420*urrrrss+cuyyyyyy030*usss+cuyyyyyy130*ursss+cuyyyyyy230*urrsss+cuyyyyyy330*urrrsss+cuyyyyyy040*ussss+cuyyyyyy140*urssss+cuyyyyyy240*urrssss+cuyyyyyy050*usssss+cuyyyyyy150*ursssss+cuyyyyyy060*ussssss
                                              uxxxxxx = cuxxxxxx100*ur+cuxxxxxx200*urr+cuxxxxxx300*urrr+cuxxxxxx400*urrrr+cuxxxxxx500*urrrrr+cuxxxxxx600*urrrrrr+cuxxxxxx010*us+cuxxxxxx110*urs+cuxxxxxx210*urrs+cuxxxxxx310*urrrs+cuxxxxxx410*urrrrs+cuxxxxxx510*urrrrrs+cuxxxxxx020*uss+cuxxxxxx120*urss+cuxxxxxx220*urrss+cuxxxxxx320*urrrss+cuxxxxxx420*urrrrss+cuxxxxxx030*usss+cuxxxxxx130*ursss+cuxxxxxx230*urrsss+cuxxxxxx330*urrrsss+cuxxxxxx040*ussss+cuxxxxxx140*urssss+cuxxxxxx240*urrssss+cuxxxxxx050*usssss+cuxxxxxx150*ursssss+cuxxxxxx060*ussssss
                       ! uxxxxxy = cuxxxxxy100*ur+cuxxxxxy200*urr+cuxxxxxy300*urrr+cuxxxxxy400*urrrr+cuxxxxxy500*urrrrr+cuxxxxxy600*urrrrrr+cuxxxxxy010*us+cuxxxxxy110*urs+cuxxxxxy210*urrs+cuxxxxxy310*urrrs+cuxxxxxy410*urrrrs+cuxxxxxy510*urrrrrs+cuxxxxxy020*uss+cuxxxxxy120*urss+cuxxxxxy220*urrss+cuxxxxxy320*urrrss+cuxxxxxy420*urrrrss+cuxxxxxy030*usss+cuxxxxxy130*ursss+cuxxxxxy230*urrsss+cuxxxxxy330*urrrsss+cuxxxxxy040*ussss+cuxxxxxy140*urssss+cuxxxxxy240*urrssss+cuxxxxxy050*usssss+cuxxxxxy150*ursssss+cuxxxxxy060*ussssss
                                              uxxxxyy = cuxxxxyy100*ur+cuxxxxyy200*urr+cuxxxxyy300*urrr+cuxxxxyy400*urrrr+cuxxxxyy500*urrrrr+cuxxxxyy600*urrrrrr+cuxxxxyy010*us+cuxxxxyy110*urs+cuxxxxyy210*urrs+cuxxxxyy310*urrrs+cuxxxxyy410*urrrrs+cuxxxxyy510*urrrrrs+cuxxxxyy020*uss+cuxxxxyy120*urss+cuxxxxyy220*urrss+cuxxxxyy320*urrrss+cuxxxxyy420*urrrrss+cuxxxxyy030*usss+cuxxxxyy130*ursss+cuxxxxyy230*urrsss+cuxxxxyy330*urrrsss+cuxxxxyy040*ussss+cuxxxxyy140*urssss+cuxxxxyy240*urrssss+cuxxxxyy050*usssss+cuxxxxyy150*ursssss+cuxxxxyy060*ussssss
                       ! uxxxyyy = cuxxxyyy100*ur+cuxxxyyy200*urr+cuxxxyyy300*urrr+cuxxxyyy400*urrrr+cuxxxyyy500*urrrrr+cuxxxyyy600*urrrrrr+cuxxxyyy010*us+cuxxxyyy110*urs+cuxxxyyy210*urrs+cuxxxyyy310*urrrs+cuxxxyyy410*urrrrs+cuxxxyyy510*urrrrrs+cuxxxyyy020*uss+cuxxxyyy120*urss+cuxxxyyy220*urrss+cuxxxyyy320*urrrss+cuxxxyyy420*urrrrss+cuxxxyyy030*usss+cuxxxyyy130*ursss+cuxxxyyy230*urrsss+cuxxxyyy330*urrrsss+cuxxxyyy040*ussss+cuxxxyyy140*urssss+cuxxxyyy240*urrssss+cuxxxyyy050*usssss+cuxxxyyy150*ursssss+cuxxxyyy060*ussssss
                                              uxxyyyy = cuxxyyyy100*ur+cuxxyyyy200*urr+cuxxyyyy300*urrr+cuxxyyyy400*urrrr+cuxxyyyy500*urrrrr+cuxxyyyy600*urrrrrr+cuxxyyyy010*us+cuxxyyyy110*urs+cuxxyyyy210*urrs+cuxxyyyy310*urrrs+cuxxyyyy410*urrrrs+cuxxyyyy510*urrrrrs+cuxxyyyy020*uss+cuxxyyyy120*urss+cuxxyyyy220*urrss+cuxxyyyy320*urrrss+cuxxyyyy420*urrrrss+cuxxyyyy030*usss+cuxxyyyy130*ursss+cuxxyyyy230*urrsss+cuxxyyyy330*urrrsss+cuxxyyyy040*ussss+cuxxyyyy140*urssss+cuxxyyyy240*urrssss+cuxxyyyy050*usssss+cuxxyyyy150*ursssss+cuxxyyyy060*ussssss
                       ! uxyyyyy = cuxyyyyy100*ur+cuxyyyyy200*urr+cuxyyyyy300*urrr+cuxyyyyy400*urrrr+cuxyyyyy500*urrrrr+cuxyyyyy600*urrrrrr+cuxyyyyy010*us+cuxyyyyy110*urs+cuxyyyyy210*urrs+cuxyyyyy310*urrrs+cuxyyyyy410*urrrrs+cuxyyyyy510*urrrrrs+cuxyyyyy020*uss+cuxyyyyy120*urss+cuxyyyyy220*urrss+cuxyyyyy320*urrrss+cuxyyyyy420*urrrrss+cuxyyyyy030*usss+cuxyyyyy130*ursss+cuxyyyyy230*urrsss+cuxyyyyy330*urrrsss+cuxyyyyy040*ussss+cuxyyyyy140*urssss+cuxyyyyy240*urrssss+cuxyyyyy050*usssss+cuxyyyyy150*ursssss+cuxyyyyy060*ussssss
                                              uyyyyyy = cuyyyyyy100*ur+cuyyyyyy200*urr+cuyyyyyy300*urrr+cuyyyyyy400*urrrr+cuyyyyyy500*urrrrr+cuyyyyyy600*urrrrrr+cuyyyyyy010*us+cuyyyyyy110*urs+cuyyyyyy210*urrs+cuyyyyyy310*urrrs+cuyyyyyy410*urrrrs+cuyyyyyy510*urrrrrs+cuyyyyyy020*uss+cuyyyyyy120*urss+cuyyyyyy220*urrss+cuyyyyyy320*urrrss+cuyyyyyy420*urrrrss+cuyyyyyy030*usss+cuyyyyyy130*ursss+cuyyyyyy230*urrsss+cuyyyyyy330*urrrsss+cuyyyyyy040*ussss+cuyyyyyy140*urssss+cuyyyyyy240*urrssss+cuyyyyyy050*usssss+cuyyyyyy150*ursssss+cuyyyyyy060*ussssss
                       ! if( i1.eq.5 .and. i2.eq.6 )then  
                       !   write(*,'(" (i1,i2)=(",2i3,") d420i,d620i,d440i=",4(1pe12.4,1x))') i1,i2,d420i,d620i,d440i
                       !   write(*,'(" (i1,i2)=(",2i3,") urrrrss,cuxxxxy420=",4(1pe12.4,1x))') i1,i2,urrrrss,cuxxxxyy420
                       ! end if 
                       ! ! ----- START SEVENTH DERIVATIVES -----
                       ! do m1=0,nd-1
                       ! do m2=0,nd-1
                       !   ! ---- 6th parameteric derivatives of the metrics ----
                       !   rxrrrrrra(m1,m2) = rx600(i1,i2,i3,m1,m2) - (dx(0)**2/4.)*rx800i(m1,m2)
                       !   rxssssssa(m1,m2) = rx060(i1,i2,i3,m1,m2) - (dx(1)**2/4.)*rx080i(m1,m2)  
                       !   rxrrrrrsa(m1,m2) = rx510i(m1,m2) - (dx(0)**2/3.)*rx710i(m1,m2) -  (dx(1)**2/6.)*rx530i(m1,m2)
                       !   rxrsssssa(m1,m2) = rx150i(m1,m2) - (dx(1)**2/3.)*rx170i(m1,m2) -  (dx(0)**2/6.)*rx350i(m1,m2)
                       !   rxrrrrssa(m1,m2) = rx420i(m1,m2) - (dx(0)**2/6.)*rx620i(m1,m2) - (dx(1)**2/12.)*rx440i(m1,m2)
                       !   rxrrssssa(m1,m2) = rx240i(m1,m2) - (dx(1)**2/6.)*rx260i(m1,m2) - (dx(0)**2/12.)*rx440i(m1,m2)
                       !   rxrrrsssa(m1,m2) = rx330i(m1,m2) - (dx(0)**2/4.)*rx530i(m1,m2) -  (dx(1)**2/4.)*rx350i(m1,m2)
                       !   ! ---- sixth spatial derivatives of the metrics ----
                       !   rxxxxxxxa(m1,m2) = cuxxxxxx100*rxr(m1,m2,0)+cuxxxxxx200*rxrra(m1,m2)+cuxxxxxx300*rxrrra(m1,m2)+cuxxxxxx400*rxrrrra(m1,m2)+cuxxxxxx500*rxrrrrra(m1,m2)+cuxxxxxx600*rxrrrrrra(m1,m2)+cuxxxxxx010*rxr(m1,m2,1)+cuxxxxxx110*rxrsa(m1,m2)+cuxxxxxx210*rxrrsa(m1,m2)+cuxxxxxx310*rxrrrsa(m1,m2)+cuxxxxxx410*rxrrrrsa(m1,m2)+cuxxxxxx510*rxrrrrrsa(m1,m2)+cuxxxxxx020*rxssa(m1,m2)+cuxxxxxx120*rxrssa(m1,m2)+cuxxxxxx220*rxrrssa(m1,m2)+cuxxxxxx320*rxrrrssa(m1,m2)+cuxxxxxx420*rxrrrrssa(m1,m2)+cuxxxxxx030*rxsssa(m1,m2)+cuxxxxxx130*rxrsssa(m1,m2)+cuxxxxxx230*rxrrsssa(m1,m2)+cuxxxxxx330*rxrrrsssa(m1,m2)+cuxxxxxx040*rxssssa(m1,m2)+cuxxxxxx140*rxrssssa(m1,m2)+cuxxxxxx240*rxrrssssa(m1,m2)+cuxxxxxx050*rxsssssa(m1,m2)+cuxxxxxx150*rxrsssssa(m1,m2)+cuxxxxxx060*rxssssssa(m1,m2)
                       !   rxxxxxxya(m1,m2) = cuxxxxxy100*rxr(m1,m2,0)+cuxxxxxy200*rxrra(m1,m2)+cuxxxxxy300*rxrrra(m1,m2)+cuxxxxxy400*rxrrrra(m1,m2)+cuxxxxxy500*rxrrrrra(m1,m2)+cuxxxxxy600*rxrrrrrra(m1,m2)+cuxxxxxy010*rxr(m1,m2,1)+cuxxxxxy110*rxrsa(m1,m2)+cuxxxxxy210*rxrrsa(m1,m2)+cuxxxxxy310*rxrrrsa(m1,m2)+cuxxxxxy410*rxrrrrsa(m1,m2)+cuxxxxxy510*rxrrrrrsa(m1,m2)+cuxxxxxy020*rxssa(m1,m2)+cuxxxxxy120*rxrssa(m1,m2)+cuxxxxxy220*rxrrssa(m1,m2)+cuxxxxxy320*rxrrrssa(m1,m2)+cuxxxxxy420*rxrrrrssa(m1,m2)+cuxxxxxy030*rxsssa(m1,m2)+cuxxxxxy130*rxrsssa(m1,m2)+cuxxxxxy230*rxrrsssa(m1,m2)+cuxxxxxy330*rxrrrsssa(m1,m2)+cuxxxxxy040*rxssssa(m1,m2)+cuxxxxxy140*rxrssssa(m1,m2)+cuxxxxxy240*rxrrssssa(m1,m2)+cuxxxxxy050*rxsssssa(m1,m2)+cuxxxxxy150*rxrsssssa(m1,m2)+cuxxxxxy060*rxssssssa(m1,m2)
                       !   rxxxxxyya(m1,m2) = cuxxxxyy100*rxr(m1,m2,0)+cuxxxxyy200*rxrra(m1,m2)+cuxxxxyy300*rxrrra(m1,m2)+cuxxxxyy400*rxrrrra(m1,m2)+cuxxxxyy500*rxrrrrra(m1,m2)+cuxxxxyy600*rxrrrrrra(m1,m2)+cuxxxxyy010*rxr(m1,m2,1)+cuxxxxyy110*rxrsa(m1,m2)+cuxxxxyy210*rxrrsa(m1,m2)+cuxxxxyy310*rxrrrsa(m1,m2)+cuxxxxyy410*rxrrrrsa(m1,m2)+cuxxxxyy510*rxrrrrrsa(m1,m2)+cuxxxxyy020*rxssa(m1,m2)+cuxxxxyy120*rxrssa(m1,m2)+cuxxxxyy220*rxrrssa(m1,m2)+cuxxxxyy320*rxrrrssa(m1,m2)+cuxxxxyy420*rxrrrrssa(m1,m2)+cuxxxxyy030*rxsssa(m1,m2)+cuxxxxyy130*rxrsssa(m1,m2)+cuxxxxyy230*rxrrsssa(m1,m2)+cuxxxxyy330*rxrrrsssa(m1,m2)+cuxxxxyy040*rxssssa(m1,m2)+cuxxxxyy140*rxrssssa(m1,m2)+cuxxxxyy240*rxrrssssa(m1,m2)+cuxxxxyy050*rxsssssa(m1,m2)+cuxxxxyy150*rxrsssssa(m1,m2)+cuxxxxyy060*rxssssssa(m1,m2)
                       !   rxxxxyyya(m1,m2) = cuxxxyyy100*rxr(m1,m2,0)+cuxxxyyy200*rxrra(m1,m2)+cuxxxyyy300*rxrrra(m1,m2)+cuxxxyyy400*rxrrrra(m1,m2)+cuxxxyyy500*rxrrrrra(m1,m2)+cuxxxyyy600*rxrrrrrra(m1,m2)+cuxxxyyy010*rxr(m1,m2,1)+cuxxxyyy110*rxrsa(m1,m2)+cuxxxyyy210*rxrrsa(m1,m2)+cuxxxyyy310*rxrrrsa(m1,m2)+cuxxxyyy410*rxrrrrsa(m1,m2)+cuxxxyyy510*rxrrrrrsa(m1,m2)+cuxxxyyy020*rxssa(m1,m2)+cuxxxyyy120*rxrssa(m1,m2)+cuxxxyyy220*rxrrssa(m1,m2)+cuxxxyyy320*rxrrrssa(m1,m2)+cuxxxyyy420*rxrrrrssa(m1,m2)+cuxxxyyy030*rxsssa(m1,m2)+cuxxxyyy130*rxrsssa(m1,m2)+cuxxxyyy230*rxrrsssa(m1,m2)+cuxxxyyy330*rxrrrsssa(m1,m2)+cuxxxyyy040*rxssssa(m1,m2)+cuxxxyyy140*rxrssssa(m1,m2)+cuxxxyyy240*rxrrssssa(m1,m2)+cuxxxyyy050*rxsssssa(m1,m2)+cuxxxyyy150*rxrsssssa(m1,m2)+cuxxxyyy060*rxssssssa(m1,m2)
                       !   rxxxyyyya(m1,m2) = cuxxyyyy100*rxr(m1,m2,0)+cuxxyyyy200*rxrra(m1,m2)+cuxxyyyy300*rxrrra(m1,m2)+cuxxyyyy400*rxrrrra(m1,m2)+cuxxyyyy500*rxrrrrra(m1,m2)+cuxxyyyy600*rxrrrrrra(m1,m2)+cuxxyyyy010*rxr(m1,m2,1)+cuxxyyyy110*rxrsa(m1,m2)+cuxxyyyy210*rxrrsa(m1,m2)+cuxxyyyy310*rxrrrsa(m1,m2)+cuxxyyyy410*rxrrrrsa(m1,m2)+cuxxyyyy510*rxrrrrrsa(m1,m2)+cuxxyyyy020*rxssa(m1,m2)+cuxxyyyy120*rxrssa(m1,m2)+cuxxyyyy220*rxrrssa(m1,m2)+cuxxyyyy320*rxrrrssa(m1,m2)+cuxxyyyy420*rxrrrrssa(m1,m2)+cuxxyyyy030*rxsssa(m1,m2)+cuxxyyyy130*rxrsssa(m1,m2)+cuxxyyyy230*rxrrsssa(m1,m2)+cuxxyyyy330*rxrrrsssa(m1,m2)+cuxxyyyy040*rxssssa(m1,m2)+cuxxyyyy140*rxrssssa(m1,m2)+cuxxyyyy240*rxrrssssa(m1,m2)+cuxxyyyy050*rxsssssa(m1,m2)+cuxxyyyy150*rxrsssssa(m1,m2)+cuxxyyyy060*rxssssssa(m1,m2)
                       !   rxxyyyyya(m1,m2) = cuxyyyyy100*rxr(m1,m2,0)+cuxyyyyy200*rxrra(m1,m2)+cuxyyyyy300*rxrrra(m1,m2)+cuxyyyyy400*rxrrrra(m1,m2)+cuxyyyyy500*rxrrrrra(m1,m2)+cuxyyyyy600*rxrrrrrra(m1,m2)+cuxyyyyy010*rxr(m1,m2,1)+cuxyyyyy110*rxrsa(m1,m2)+cuxyyyyy210*rxrrsa(m1,m2)+cuxyyyyy310*rxrrrsa(m1,m2)+cuxyyyyy410*rxrrrrsa(m1,m2)+cuxyyyyy510*rxrrrrrsa(m1,m2)+cuxyyyyy020*rxssa(m1,m2)+cuxyyyyy120*rxrssa(m1,m2)+cuxyyyyy220*rxrrssa(m1,m2)+cuxyyyyy320*rxrrrssa(m1,m2)+cuxyyyyy420*rxrrrrssa(m1,m2)+cuxyyyyy030*rxsssa(m1,m2)+cuxyyyyy130*rxrsssa(m1,m2)+cuxyyyyy230*rxrrsssa(m1,m2)+cuxyyyyy330*rxrrrsssa(m1,m2)+cuxyyyyy040*rxssssa(m1,m2)+cuxyyyyy140*rxrssssa(m1,m2)+cuxyyyyy240*rxrrssssa(m1,m2)+cuxyyyyy050*rxsssssa(m1,m2)+cuxyyyyy150*rxrsssssa(m1,m2)+cuxyyyyy060*rxssssssa(m1,m2)
                       !   rxyyyyyya(m1,m2) = cuyyyyyy100*rxr(m1,m2,0)+cuyyyyyy200*rxrra(m1,m2)+cuyyyyyy300*rxrrra(m1,m2)+cuyyyyyy400*rxrrrra(m1,m2)+cuyyyyyy500*rxrrrrra(m1,m2)+cuyyyyyy600*rxrrrrrra(m1,m2)+cuyyyyyy010*rxr(m1,m2,1)+cuyyyyyy110*rxrsa(m1,m2)+cuyyyyyy210*rxrrsa(m1,m2)+cuyyyyyy310*rxrrrsa(m1,m2)+cuyyyyyy410*rxrrrrsa(m1,m2)+cuyyyyyy510*rxrrrrrsa(m1,m2)+cuyyyyyy020*rxssa(m1,m2)+cuyyyyyy120*rxrssa(m1,m2)+cuyyyyyy220*rxrrssa(m1,m2)+cuyyyyyy320*rxrrrssa(m1,m2)+cuyyyyyy420*rxrrrrssa(m1,m2)+cuyyyyyy030*rxsssa(m1,m2)+cuyyyyyy130*rxrsssa(m1,m2)+cuyyyyyy230*rxrrsssa(m1,m2)+cuyyyyyy330*rxrrrsssa(m1,m2)+cuyyyyyy040*rxssssa(m1,m2)+cuyyyyyy140*rxrssssa(m1,m2)+cuyyyyyy240*rxrrssssa(m1,m2)+cuyyyyyy050*rxsssssa(m1,m2)+cuyyyyyy150*rxrsssssa(m1,m2)+cuyyyyyy060*rxssssssa(m1,m2)
                       ! end do
                       ! end do
                       ! rxxxxxxx=rxxxxxxxa(0,0); ryxxxxxx=rxxxxxxxa(0,1); sxxxxxxx=rxxxxxxxa(1,0); syxxxxxx=rxxxxxxxa(1,1); 
                       ! rxxxxxxy=rxxxxxxya(0,0); ryxxxxxy=rxxxxxxya(0,1); sxxxxxxy=rxxxxxxya(1,0); syxxxxxy=rxxxxxxya(1,1); 
                       ! rxxxxxyy=rxxxxxyya(0,0); ryxxxxyy=rxxxxxyya(0,1); sxxxxxyy=rxxxxxyya(1,0); syxxxxyy=rxxxxxyya(1,1); 
                       ! rxxxxyyy=rxxxxyyya(0,0); ryxxxyyy=rxxxxyyya(0,1); sxxxxyyy=rxxxxyyya(1,0); syxxxyyy=rxxxxyyya(1,1); 
                       ! rxxxyyyy=rxxxyyyya(0,0); ryxxyyyy=rxxxyyyya(0,1); sxxxyyyy=rxxxyyyya(1,0); syxxyyyy=rxxxyyyya(1,1); 
                       ! rxxyyyyy=rxxyyyyya(0,0); ryxyyyyy=rxxyyyyya(0,1); sxxyyyyy=rxxyyyyya(1,0); syxyyyyy=rxxyyyyya(1,1); 
                       ! rxyyyyyy=rxyyyyyya(0,0); ryyyyyyy=rxyyyyyya(0,1); sxyyyyyy=rxyyyyyya(1,0); syyyyyyy=rxyyyyyya(1,1); 
                       ! ! ---- seventh parametric derivatives ----
                       ! urrrrrrr = d700i
                       ! urrrrrrs = d610i
                       ! urrrrrss = d520i
                       ! urrrrsss = d430i
                       ! urrrssss = d340i
                       ! urrsssss = d250i
                       ! urssssss = d160i
                       ! usssssss = d070i 
                       ! ! ----- SEVENTH SPATIAL DERIVATIVES -----
                       ! getDerivCoeff2d(7)
                       ! ! uxxxxxxx = cuxxxxxxx100*ur+cuxxxxxxx200*urr+cuxxxxxxx300*urrr+cuxxxxxxx400*urrrr+cuxxxxxxx500*urrrrr+cuxxxxxxx600*urrrrrr+cuxxxxxxx700*urrrrrrr+cuxxxxxxx010*us+cuxxxxxxx110*urs+cuxxxxxxx210*urrs+cuxxxxxxx310*urrrs+cuxxxxxxx410*urrrrs+cuxxxxxxx510*urrrrrs+cuxxxxxxx610*urrrrrrs+cuxxxxxxx020*uss+cuxxxxxxx120*urss+cuxxxxxxx220*urrss+cuxxxxxxx320*urrrss+cuxxxxxxx420*urrrrss+cuxxxxxxx520*urrrrrss+cuxxxxxxx030*usss+cuxxxxxxx130*ursss+cuxxxxxxx230*urrsss+cuxxxxxxx330*urrrsss+cuxxxxxxx430*urrrrsss+cuxxxxxxx040*ussss+cuxxxxxxx140*urssss+cuxxxxxxx240*urrssss+cuxxxxxxx340*urrrssss+cuxxxxxxx050*usssss+cuxxxxxxx150*ursssss+cuxxxxxxx250*urrsssss+cuxxxxxxx060*ussssss+cuxxxxxxx160*urssssss+cuxxxxxxx070*usssssss
                       ! ! uxxxxxxy = cuxxxxxxy100*ur+cuxxxxxxy200*urr+cuxxxxxxy300*urrr+cuxxxxxxy400*urrrr+cuxxxxxxy500*urrrrr+cuxxxxxxy600*urrrrrr+cuxxxxxxy700*urrrrrrr+cuxxxxxxy010*us+cuxxxxxxy110*urs+cuxxxxxxy210*urrs+cuxxxxxxy310*urrrs+cuxxxxxxy410*urrrrs+cuxxxxxxy510*urrrrrs+cuxxxxxxy610*urrrrrrs+cuxxxxxxy020*uss+cuxxxxxxy120*urss+cuxxxxxxy220*urrss+cuxxxxxxy320*urrrss+cuxxxxxxy420*urrrrss+cuxxxxxxy520*urrrrrss+cuxxxxxxy030*usss+cuxxxxxxy130*ursss+cuxxxxxxy230*urrsss+cuxxxxxxy330*urrrsss+cuxxxxxxy430*urrrrsss+cuxxxxxxy040*ussss+cuxxxxxxy140*urssss+cuxxxxxxy240*urrssss+cuxxxxxxy340*urrrssss+cuxxxxxxy050*usssss+cuxxxxxxy150*ursssss+cuxxxxxxy250*urrsssss+cuxxxxxxy060*ussssss+cuxxxxxxy160*urssssss+cuxxxxxxy070*usssssss
                       ! ! uxxxxxyy = cuxxxxxyy100*ur+cuxxxxxyy200*urr+cuxxxxxyy300*urrr+cuxxxxxyy400*urrrr+cuxxxxxyy500*urrrrr+cuxxxxxyy600*urrrrrr+cuxxxxxyy700*urrrrrrr+cuxxxxxyy010*us+cuxxxxxyy110*urs+cuxxxxxyy210*urrs+cuxxxxxyy310*urrrs+cuxxxxxyy410*urrrrs+cuxxxxxyy510*urrrrrs+cuxxxxxyy610*urrrrrrs+cuxxxxxyy020*uss+cuxxxxxyy120*urss+cuxxxxxyy220*urrss+cuxxxxxyy320*urrrss+cuxxxxxyy420*urrrrss+cuxxxxxyy520*urrrrrss+cuxxxxxyy030*usss+cuxxxxxyy130*ursss+cuxxxxxyy230*urrsss+cuxxxxxyy330*urrrsss+cuxxxxxyy430*urrrrsss+cuxxxxxyy040*ussss+cuxxxxxyy140*urssss+cuxxxxxyy240*urrssss+cuxxxxxyy340*urrrssss+cuxxxxxyy050*usssss+cuxxxxxyy150*ursssss+cuxxxxxyy250*urrsssss+cuxxxxxyy060*ussssss+cuxxxxxyy160*urssssss+cuxxxxxyy070*usssssss
                       ! ! uxxxxyyy = cuxxxxyyy100*ur+cuxxxxyyy200*urr+cuxxxxyyy300*urrr+cuxxxxyyy400*urrrr+cuxxxxyyy500*urrrrr+cuxxxxyyy600*urrrrrr+cuxxxxyyy700*urrrrrrr+cuxxxxyyy010*us+cuxxxxyyy110*urs+cuxxxxyyy210*urrs+cuxxxxyyy310*urrrs+cuxxxxyyy410*urrrrs+cuxxxxyyy510*urrrrrs+cuxxxxyyy610*urrrrrrs+cuxxxxyyy020*uss+cuxxxxyyy120*urss+cuxxxxyyy220*urrss+cuxxxxyyy320*urrrss+cuxxxxyyy420*urrrrss+cuxxxxyyy520*urrrrrss+cuxxxxyyy030*usss+cuxxxxyyy130*ursss+cuxxxxyyy230*urrsss+cuxxxxyyy330*urrrsss+cuxxxxyyy430*urrrrsss+cuxxxxyyy040*ussss+cuxxxxyyy140*urssss+cuxxxxyyy240*urrssss+cuxxxxyyy340*urrrssss+cuxxxxyyy050*usssss+cuxxxxyyy150*ursssss+cuxxxxyyy250*urrsssss+cuxxxxyyy060*ussssss+cuxxxxyyy160*urssssss+cuxxxxyyy070*usssssss
                       ! ! uxxxyyyy = cuxxxyyyy100*ur+cuxxxyyyy200*urr+cuxxxyyyy300*urrr+cuxxxyyyy400*urrrr+cuxxxyyyy500*urrrrr+cuxxxyyyy600*urrrrrr+cuxxxyyyy700*urrrrrrr+cuxxxyyyy010*us+cuxxxyyyy110*urs+cuxxxyyyy210*urrs+cuxxxyyyy310*urrrs+cuxxxyyyy410*urrrrs+cuxxxyyyy510*urrrrrs+cuxxxyyyy610*urrrrrrs+cuxxxyyyy020*uss+cuxxxyyyy120*urss+cuxxxyyyy220*urrss+cuxxxyyyy320*urrrss+cuxxxyyyy420*urrrrss+cuxxxyyyy520*urrrrrss+cuxxxyyyy030*usss+cuxxxyyyy130*ursss+cuxxxyyyy230*urrsss+cuxxxyyyy330*urrrsss+cuxxxyyyy430*urrrrsss+cuxxxyyyy040*ussss+cuxxxyyyy140*urssss+cuxxxyyyy240*urrssss+cuxxxyyyy340*urrrssss+cuxxxyyyy050*usssss+cuxxxyyyy150*ursssss+cuxxxyyyy250*urrsssss+cuxxxyyyy060*ussssss+cuxxxyyyy160*urssssss+cuxxxyyyy070*usssssss
                       ! ! uxxyyyyy = cuxxyyyyy100*ur+cuxxyyyyy200*urr+cuxxyyyyy300*urrr+cuxxyyyyy400*urrrr+cuxxyyyyy500*urrrrr+cuxxyyyyy600*urrrrrr+cuxxyyyyy700*urrrrrrr+cuxxyyyyy010*us+cuxxyyyyy110*urs+cuxxyyyyy210*urrs+cuxxyyyyy310*urrrs+cuxxyyyyy410*urrrrs+cuxxyyyyy510*urrrrrs+cuxxyyyyy610*urrrrrrs+cuxxyyyyy020*uss+cuxxyyyyy120*urss+cuxxyyyyy220*urrss+cuxxyyyyy320*urrrss+cuxxyyyyy420*urrrrss+cuxxyyyyy520*urrrrrss+cuxxyyyyy030*usss+cuxxyyyyy130*ursss+cuxxyyyyy230*urrsss+cuxxyyyyy330*urrrsss+cuxxyyyyy430*urrrrsss+cuxxyyyyy040*ussss+cuxxyyyyy140*urssss+cuxxyyyyy240*urrssss+cuxxyyyyy340*urrrssss+cuxxyyyyy050*usssss+cuxxyyyyy150*ursssss+cuxxyyyyy250*urrsssss+cuxxyyyyy060*ussssss+cuxxyyyyy160*urssssss+cuxxyyyyy070*usssssss
                       ! ! uxyyyyyy = cuxyyyyyy100*ur+cuxyyyyyy200*urr+cuxyyyyyy300*urrr+cuxyyyyyy400*urrrr+cuxyyyyyy500*urrrrr+cuxyyyyyy600*urrrrrr+cuxyyyyyy700*urrrrrrr+cuxyyyyyy010*us+cuxyyyyyy110*urs+cuxyyyyyy210*urrs+cuxyyyyyy310*urrrs+cuxyyyyyy410*urrrrs+cuxyyyyyy510*urrrrrs+cuxyyyyyy610*urrrrrrs+cuxyyyyyy020*uss+cuxyyyyyy120*urss+cuxyyyyyy220*urrss+cuxyyyyyy320*urrrss+cuxyyyyyy420*urrrrss+cuxyyyyyy520*urrrrrss+cuxyyyyyy030*usss+cuxyyyyyy130*ursss+cuxyyyyyy230*urrsss+cuxyyyyyy330*urrrsss+cuxyyyyyy430*urrrrsss+cuxyyyyyy040*ussss+cuxyyyyyy140*urssss+cuxyyyyyy240*urrssss+cuxyyyyyy340*urrrssss+cuxyyyyyy050*usssss+cuxyyyyyy150*ursssss+cuxyyyyyy250*urrsssss+cuxyyyyyy060*ussssss+cuxyyyyyy160*urssssss+cuxyyyyyy070*usssssss
                       ! ! uyyyyyyy = cuyyyyyyy100*ur+cuyyyyyyy200*urr+cuyyyyyyy300*urrr+cuyyyyyyy400*urrrr+cuyyyyyyy500*urrrrr+cuyyyyyyy600*urrrrrr+cuyyyyyyy700*urrrrrrr+cuyyyyyyy010*us+cuyyyyyyy110*urs+cuyyyyyyy210*urrs+cuyyyyyyy310*urrrs+cuyyyyyyy410*urrrrs+cuyyyyyyy510*urrrrrs+cuyyyyyyy610*urrrrrrs+cuyyyyyyy020*uss+cuyyyyyyy120*urss+cuyyyyyyy220*urrss+cuyyyyyyy320*urrrss+cuyyyyyyy420*urrrrss+cuyyyyyyy520*urrrrrss+cuyyyyyyy030*usss+cuyyyyyyy130*ursss+cuyyyyyyy230*urrsss+cuyyyyyyy330*urrrsss+cuyyyyyyy430*urrrrsss+cuyyyyyyy040*ussss+cuyyyyyyy140*urssss+cuyyyyyyy240*urrssss+cuyyyyyyy340*urrrssss+cuyyyyyyy050*usssss+cuyyyyyyy150*ursssss+cuyyyyyyy250*urrsssss+cuyyyyyyy060*ussssss+cuyyyyyyy160*urssssss+cuyyyyyyy070*usssssss
                       ! ! ----- START EIGHTH DERIVATIVES -----
                       ! do m1=0,nd-1
                       ! do m2=0,nd-1
                       !   ! ---- 7th parameteric derivatives of the metrics ----
                       !   rxrrrrrrra(m1,m2) = rx700i(m1,m2)
                       !   rxrrrrrrsa(m1,m2) = rx610i(m1,m2)
                       !   rxrrrrrssa(m1,m2) = rx520i(m1,m2)
                       !   rxrrrrsssa(m1,m2) = rx430i(m1,m2)
                       !   rxrrrssssa(m1,m2) = rx340i(m1,m2)
                       !   rxrrsssssa(m1,m2) = rx250i(m1,m2)
                       !   rxrssssssa(m1,m2) = rx160i(m1,m2)
                       !   rxsssssssa(m1,m2) = rx070i(m1,m2) 
                       !   ! ---- seventh spatial derivatives of the metrics ----
                       !   rxxxxxxxxa(m1,m2) = cuxxxxxxx100*rxr(m1,m2,0)+cuxxxxxxx200*rxrra(m1,m2)+cuxxxxxxx300*rxrrra(m1,m2)+cuxxxxxxx400*rxrrrra(m1,m2)+cuxxxxxxx500*rxrrrrra(m1,m2)+cuxxxxxxx600*rxrrrrrra(m1,m2)+cuxxxxxxx700*rxrrrrrrra(m1,m2)+cuxxxxxxx010*rxr(m1,m2,1)+cuxxxxxxx110*rxrsa(m1,m2)+cuxxxxxxx210*rxrrsa(m1,m2)+cuxxxxxxx310*rxrrrsa(m1,m2)+cuxxxxxxx410*rxrrrrsa(m1,m2)+cuxxxxxxx510*rxrrrrrsa(m1,m2)+cuxxxxxxx610*rxrrrrrrsa(m1,m2)+cuxxxxxxx020*rxssa(m1,m2)+cuxxxxxxx120*rxrssa(m1,m2)+cuxxxxxxx220*rxrrssa(m1,m2)+cuxxxxxxx320*rxrrrssa(m1,m2)+cuxxxxxxx420*rxrrrrssa(m1,m2)+cuxxxxxxx520*rxrrrrrssa(m1,m2)+cuxxxxxxx030*rxsssa(m1,m2)+cuxxxxxxx130*rxrsssa(m1,m2)+cuxxxxxxx230*rxrrsssa(m1,m2)+cuxxxxxxx330*rxrrrsssa(m1,m2)+cuxxxxxxx430*rxrrrrsssa(m1,m2)+cuxxxxxxx040*rxssssa(m1,m2)+cuxxxxxxx140*rxrssssa(m1,m2)+cuxxxxxxx240*rxrrssssa(m1,m2)+cuxxxxxxx340*rxrrrssssa(m1,m2)+cuxxxxxxx050*rxsssssa(m1,m2)+cuxxxxxxx150*rxrsssssa(m1,m2)+cuxxxxxxx250*rxrrsssssa(m1,m2)+cuxxxxxxx060*rxssssssa(m1,m2)+cuxxxxxxx160*rxrssssssa(m1,m2)+cuxxxxxxx070*rxsssssssa(m1,m2)
                       !   rxxxxxxxya(m1,m2) = cuxxxxxxy100*rxr(m1,m2,0)+cuxxxxxxy200*rxrra(m1,m2)+cuxxxxxxy300*rxrrra(m1,m2)+cuxxxxxxy400*rxrrrra(m1,m2)+cuxxxxxxy500*rxrrrrra(m1,m2)+cuxxxxxxy600*rxrrrrrra(m1,m2)+cuxxxxxxy700*rxrrrrrrra(m1,m2)+cuxxxxxxy010*rxr(m1,m2,1)+cuxxxxxxy110*rxrsa(m1,m2)+cuxxxxxxy210*rxrrsa(m1,m2)+cuxxxxxxy310*rxrrrsa(m1,m2)+cuxxxxxxy410*rxrrrrsa(m1,m2)+cuxxxxxxy510*rxrrrrrsa(m1,m2)+cuxxxxxxy610*rxrrrrrrsa(m1,m2)+cuxxxxxxy020*rxssa(m1,m2)+cuxxxxxxy120*rxrssa(m1,m2)+cuxxxxxxy220*rxrrssa(m1,m2)+cuxxxxxxy320*rxrrrssa(m1,m2)+cuxxxxxxy420*rxrrrrssa(m1,m2)+cuxxxxxxy520*rxrrrrrssa(m1,m2)+cuxxxxxxy030*rxsssa(m1,m2)+cuxxxxxxy130*rxrsssa(m1,m2)+cuxxxxxxy230*rxrrsssa(m1,m2)+cuxxxxxxy330*rxrrrsssa(m1,m2)+cuxxxxxxy430*rxrrrrsssa(m1,m2)+cuxxxxxxy040*rxssssa(m1,m2)+cuxxxxxxy140*rxrssssa(m1,m2)+cuxxxxxxy240*rxrrssssa(m1,m2)+cuxxxxxxy340*rxrrrssssa(m1,m2)+cuxxxxxxy050*rxsssssa(m1,m2)+cuxxxxxxy150*rxrsssssa(m1,m2)+cuxxxxxxy250*rxrrsssssa(m1,m2)+cuxxxxxxy060*rxssssssa(m1,m2)+cuxxxxxxy160*rxrssssssa(m1,m2)+cuxxxxxxy070*rxsssssssa(m1,m2)
                       !   rxxxxxxyya(m1,m2) = cuxxxxxyy100*rxr(m1,m2,0)+cuxxxxxyy200*rxrra(m1,m2)+cuxxxxxyy300*rxrrra(m1,m2)+cuxxxxxyy400*rxrrrra(m1,m2)+cuxxxxxyy500*rxrrrrra(m1,m2)+cuxxxxxyy600*rxrrrrrra(m1,m2)+cuxxxxxyy700*rxrrrrrrra(m1,m2)+cuxxxxxyy010*rxr(m1,m2,1)+cuxxxxxyy110*rxrsa(m1,m2)+cuxxxxxyy210*rxrrsa(m1,m2)+cuxxxxxyy310*rxrrrsa(m1,m2)+cuxxxxxyy410*rxrrrrsa(m1,m2)+cuxxxxxyy510*rxrrrrrsa(m1,m2)+cuxxxxxyy610*rxrrrrrrsa(m1,m2)+cuxxxxxyy020*rxssa(m1,m2)+cuxxxxxyy120*rxrssa(m1,m2)+cuxxxxxyy220*rxrrssa(m1,m2)+cuxxxxxyy320*rxrrrssa(m1,m2)+cuxxxxxyy420*rxrrrrssa(m1,m2)+cuxxxxxyy520*rxrrrrrssa(m1,m2)+cuxxxxxyy030*rxsssa(m1,m2)+cuxxxxxyy130*rxrsssa(m1,m2)+cuxxxxxyy230*rxrrsssa(m1,m2)+cuxxxxxyy330*rxrrrsssa(m1,m2)+cuxxxxxyy430*rxrrrrsssa(m1,m2)+cuxxxxxyy040*rxssssa(m1,m2)+cuxxxxxyy140*rxrssssa(m1,m2)+cuxxxxxyy240*rxrrssssa(m1,m2)+cuxxxxxyy340*rxrrrssssa(m1,m2)+cuxxxxxyy050*rxsssssa(m1,m2)+cuxxxxxyy150*rxrsssssa(m1,m2)+cuxxxxxyy250*rxrrsssssa(m1,m2)+cuxxxxxyy060*rxssssssa(m1,m2)+cuxxxxxyy160*rxrssssssa(m1,m2)+cuxxxxxyy070*rxsssssssa(m1,m2)
                       !   rxxxxxyyya(m1,m2) = cuxxxxyyy100*rxr(m1,m2,0)+cuxxxxyyy200*rxrra(m1,m2)+cuxxxxyyy300*rxrrra(m1,m2)+cuxxxxyyy400*rxrrrra(m1,m2)+cuxxxxyyy500*rxrrrrra(m1,m2)+cuxxxxyyy600*rxrrrrrra(m1,m2)+cuxxxxyyy700*rxrrrrrrra(m1,m2)+cuxxxxyyy010*rxr(m1,m2,1)+cuxxxxyyy110*rxrsa(m1,m2)+cuxxxxyyy210*rxrrsa(m1,m2)+cuxxxxyyy310*rxrrrsa(m1,m2)+cuxxxxyyy410*rxrrrrsa(m1,m2)+cuxxxxyyy510*rxrrrrrsa(m1,m2)+cuxxxxyyy610*rxrrrrrrsa(m1,m2)+cuxxxxyyy020*rxssa(m1,m2)+cuxxxxyyy120*rxrssa(m1,m2)+cuxxxxyyy220*rxrrssa(m1,m2)+cuxxxxyyy320*rxrrrssa(m1,m2)+cuxxxxyyy420*rxrrrrssa(m1,m2)+cuxxxxyyy520*rxrrrrrssa(m1,m2)+cuxxxxyyy030*rxsssa(m1,m2)+cuxxxxyyy130*rxrsssa(m1,m2)+cuxxxxyyy230*rxrrsssa(m1,m2)+cuxxxxyyy330*rxrrrsssa(m1,m2)+cuxxxxyyy430*rxrrrrsssa(m1,m2)+cuxxxxyyy040*rxssssa(m1,m2)+cuxxxxyyy140*rxrssssa(m1,m2)+cuxxxxyyy240*rxrrssssa(m1,m2)+cuxxxxyyy340*rxrrrssssa(m1,m2)+cuxxxxyyy050*rxsssssa(m1,m2)+cuxxxxyyy150*rxrsssssa(m1,m2)+cuxxxxyyy250*rxrrsssssa(m1,m2)+cuxxxxyyy060*rxssssssa(m1,m2)+cuxxxxyyy160*rxrssssssa(m1,m2)+cuxxxxyyy070*rxsssssssa(m1,m2)
                       !   rxxxxyyyya(m1,m2) = cuxxxyyyy100*rxr(m1,m2,0)+cuxxxyyyy200*rxrra(m1,m2)+cuxxxyyyy300*rxrrra(m1,m2)+cuxxxyyyy400*rxrrrra(m1,m2)+cuxxxyyyy500*rxrrrrra(m1,m2)+cuxxxyyyy600*rxrrrrrra(m1,m2)+cuxxxyyyy700*rxrrrrrrra(m1,m2)+cuxxxyyyy010*rxr(m1,m2,1)+cuxxxyyyy110*rxrsa(m1,m2)+cuxxxyyyy210*rxrrsa(m1,m2)+cuxxxyyyy310*rxrrrsa(m1,m2)+cuxxxyyyy410*rxrrrrsa(m1,m2)+cuxxxyyyy510*rxrrrrrsa(m1,m2)+cuxxxyyyy610*rxrrrrrrsa(m1,m2)+cuxxxyyyy020*rxssa(m1,m2)+cuxxxyyyy120*rxrssa(m1,m2)+cuxxxyyyy220*rxrrssa(m1,m2)+cuxxxyyyy320*rxrrrssa(m1,m2)+cuxxxyyyy420*rxrrrrssa(m1,m2)+cuxxxyyyy520*rxrrrrrssa(m1,m2)+cuxxxyyyy030*rxsssa(m1,m2)+cuxxxyyyy130*rxrsssa(m1,m2)+cuxxxyyyy230*rxrrsssa(m1,m2)+cuxxxyyyy330*rxrrrsssa(m1,m2)+cuxxxyyyy430*rxrrrrsssa(m1,m2)+cuxxxyyyy040*rxssssa(m1,m2)+cuxxxyyyy140*rxrssssa(m1,m2)+cuxxxyyyy240*rxrrssssa(m1,m2)+cuxxxyyyy340*rxrrrssssa(m1,m2)+cuxxxyyyy050*rxsssssa(m1,m2)+cuxxxyyyy150*rxrsssssa(m1,m2)+cuxxxyyyy250*rxrrsssssa(m1,m2)+cuxxxyyyy060*rxssssssa(m1,m2)+cuxxxyyyy160*rxrssssssa(m1,m2)+cuxxxyyyy070*rxsssssssa(m1,m2)
                       !   rxxxyyyyya(m1,m2) = cuxxyyyyy100*rxr(m1,m2,0)+cuxxyyyyy200*rxrra(m1,m2)+cuxxyyyyy300*rxrrra(m1,m2)+cuxxyyyyy400*rxrrrra(m1,m2)+cuxxyyyyy500*rxrrrrra(m1,m2)+cuxxyyyyy600*rxrrrrrra(m1,m2)+cuxxyyyyy700*rxrrrrrrra(m1,m2)+cuxxyyyyy010*rxr(m1,m2,1)+cuxxyyyyy110*rxrsa(m1,m2)+cuxxyyyyy210*rxrrsa(m1,m2)+cuxxyyyyy310*rxrrrsa(m1,m2)+cuxxyyyyy410*rxrrrrsa(m1,m2)+cuxxyyyyy510*rxrrrrrsa(m1,m2)+cuxxyyyyy610*rxrrrrrrsa(m1,m2)+cuxxyyyyy020*rxssa(m1,m2)+cuxxyyyyy120*rxrssa(m1,m2)+cuxxyyyyy220*rxrrssa(m1,m2)+cuxxyyyyy320*rxrrrssa(m1,m2)+cuxxyyyyy420*rxrrrrssa(m1,m2)+cuxxyyyyy520*rxrrrrrssa(m1,m2)+cuxxyyyyy030*rxsssa(m1,m2)+cuxxyyyyy130*rxrsssa(m1,m2)+cuxxyyyyy230*rxrrsssa(m1,m2)+cuxxyyyyy330*rxrrrsssa(m1,m2)+cuxxyyyyy430*rxrrrrsssa(m1,m2)+cuxxyyyyy040*rxssssa(m1,m2)+cuxxyyyyy140*rxrssssa(m1,m2)+cuxxyyyyy240*rxrrssssa(m1,m2)+cuxxyyyyy340*rxrrrssssa(m1,m2)+cuxxyyyyy050*rxsssssa(m1,m2)+cuxxyyyyy150*rxrsssssa(m1,m2)+cuxxyyyyy250*rxrrsssssa(m1,m2)+cuxxyyyyy060*rxssssssa(m1,m2)+cuxxyyyyy160*rxrssssssa(m1,m2)+cuxxyyyyy070*rxsssssssa(m1,m2)
                       !   rxxyyyyyya(m1,m2) = cuxyyyyyy100*rxr(m1,m2,0)+cuxyyyyyy200*rxrra(m1,m2)+cuxyyyyyy300*rxrrra(m1,m2)+cuxyyyyyy400*rxrrrra(m1,m2)+cuxyyyyyy500*rxrrrrra(m1,m2)+cuxyyyyyy600*rxrrrrrra(m1,m2)+cuxyyyyyy700*rxrrrrrrra(m1,m2)+cuxyyyyyy010*rxr(m1,m2,1)+cuxyyyyyy110*rxrsa(m1,m2)+cuxyyyyyy210*rxrrsa(m1,m2)+cuxyyyyyy310*rxrrrsa(m1,m2)+cuxyyyyyy410*rxrrrrsa(m1,m2)+cuxyyyyyy510*rxrrrrrsa(m1,m2)+cuxyyyyyy610*rxrrrrrrsa(m1,m2)+cuxyyyyyy020*rxssa(m1,m2)+cuxyyyyyy120*rxrssa(m1,m2)+cuxyyyyyy220*rxrrssa(m1,m2)+cuxyyyyyy320*rxrrrssa(m1,m2)+cuxyyyyyy420*rxrrrrssa(m1,m2)+cuxyyyyyy520*rxrrrrrssa(m1,m2)+cuxyyyyyy030*rxsssa(m1,m2)+cuxyyyyyy130*rxrsssa(m1,m2)+cuxyyyyyy230*rxrrsssa(m1,m2)+cuxyyyyyy330*rxrrrsssa(m1,m2)+cuxyyyyyy430*rxrrrrsssa(m1,m2)+cuxyyyyyy040*rxssssa(m1,m2)+cuxyyyyyy140*rxrssssa(m1,m2)+cuxyyyyyy240*rxrrssssa(m1,m2)+cuxyyyyyy340*rxrrrssssa(m1,m2)+cuxyyyyyy050*rxsssssa(m1,m2)+cuxyyyyyy150*rxrsssssa(m1,m2)+cuxyyyyyy250*rxrrsssssa(m1,m2)+cuxyyyyyy060*rxssssssa(m1,m2)+cuxyyyyyy160*rxrssssssa(m1,m2)+cuxyyyyyy070*rxsssssssa(m1,m2)
                       !   rxyyyyyyya(m1,m2) = cuyyyyyyy100*rxr(m1,m2,0)+cuyyyyyyy200*rxrra(m1,m2)+cuyyyyyyy300*rxrrra(m1,m2)+cuyyyyyyy400*rxrrrra(m1,m2)+cuyyyyyyy500*rxrrrrra(m1,m2)+cuyyyyyyy600*rxrrrrrra(m1,m2)+cuyyyyyyy700*rxrrrrrrra(m1,m2)+cuyyyyyyy010*rxr(m1,m2,1)+cuyyyyyyy110*rxrsa(m1,m2)+cuyyyyyyy210*rxrrsa(m1,m2)+cuyyyyyyy310*rxrrrsa(m1,m2)+cuyyyyyyy410*rxrrrrsa(m1,m2)+cuyyyyyyy510*rxrrrrrsa(m1,m2)+cuyyyyyyy610*rxrrrrrrsa(m1,m2)+cuyyyyyyy020*rxssa(m1,m2)+cuyyyyyyy120*rxrssa(m1,m2)+cuyyyyyyy220*rxrrssa(m1,m2)+cuyyyyyyy320*rxrrrssa(m1,m2)+cuyyyyyyy420*rxrrrrssa(m1,m2)+cuyyyyyyy520*rxrrrrrssa(m1,m2)+cuyyyyyyy030*rxsssa(m1,m2)+cuyyyyyyy130*rxrsssa(m1,m2)+cuyyyyyyy230*rxrrsssa(m1,m2)+cuyyyyyyy330*rxrrrsssa(m1,m2)+cuyyyyyyy430*rxrrrrsssa(m1,m2)+cuyyyyyyy040*rxssssa(m1,m2)+cuyyyyyyy140*rxrssssa(m1,m2)+cuyyyyyyy240*rxrrssssa(m1,m2)+cuyyyyyyy340*rxrrrssssa(m1,m2)+cuyyyyyyy050*rxsssssa(m1,m2)+cuyyyyyyy150*rxrsssssa(m1,m2)+cuyyyyyyy250*rxrrsssssa(m1,m2)+cuyyyyyyy060*rxssssssa(m1,m2)+cuyyyyyyy160*rxrssssssa(m1,m2)+cuyyyyyyy070*rxsssssssa(m1,m2)
                       ! end do
                       ! end do
                       ! rxxxxxxxx=rxxxxxxxxa(0,0); ryxxxxxxx=rxxxxxxxxa(0,1); sxxxxxxxx=rxxxxxxxxa(1,0); syxxxxxxx=rxxxxxxxxa(1,1); 
                       ! rxxxxxxxy=rxxxxxxxya(0,0); ryxxxxxxy=rxxxxxxxya(0,1); sxxxxxxxy=rxxxxxxxya(1,0); syxxxxxxy=rxxxxxxxya(1,1); 
                       ! rxxxxxxyy=rxxxxxxyya(0,0); ryxxxxxyy=rxxxxxxyya(0,1); sxxxxxxyy=rxxxxxxyya(1,0); syxxxxxyy=rxxxxxxyya(1,1); 
                       ! rxxxxxyyy=rxxxxxyyya(0,0); ryxxxxyyy=rxxxxxyyya(0,1); sxxxxxyyy=rxxxxxyyya(1,0); syxxxxyyy=rxxxxxyyya(1,1); 
                       ! rxxxxyyyy=rxxxxyyyya(0,0); ryxxxyyyy=rxxxxyyyya(0,1); sxxxxyyyy=rxxxxyyyya(1,0); syxxxyyyy=rxxxxyyyya(1,1); 
                       ! rxxxyyyyy=rxxxyyyyya(0,0); ryxxyyyyy=rxxxyyyyya(0,1); sxxxyyyyy=rxxxyyyyya(1,0); syxxyyyyy=rxxxyyyyya(1,1); 
                       ! rxxyyyyyy=rxxyyyyyya(0,0); ryxyyyyyy=rxxyyyyyya(0,1); sxxyyyyyy=rxxyyyyyya(1,0); syxyyyyyy=rxxyyyyyya(1,1); 
                       ! rxyyyyyyy=rxyyyyyyya(0,0); ryyyyyyyy=rxyyyyyyya(0,1); sxyyyyyyy=rxyyyyyyya(1,0); syyyyyyyy=rxyyyyyyya(1,1); 
                       ! ! ---- eighth parametric derivatives ----
                       ! urrrrrrrr = d800i
                       ! urrrrrrrs = d710i
                       ! urrrrrrss = d620i
                       ! urrrrrsss = d530i
                       ! urrrrssss = d440i
                       ! urrrsssss = d350i
                       ! urrssssss = d260i
                       ! ursssssss = d170i
                       ! ussssssss = d080i
                       ! ! ----- EIGHTH SPATIAL DERIVATIVES -----
                       ! getDerivCoeff2d(8)  
                       ! uxxxxxxxx = cuxxxxxxxx100*ur+cuxxxxxxxx200*urr+cuxxxxxxxx300*urrr+cuxxxxxxxx400*urrrr+cuxxxxxxxx500*urrrrr+cuxxxxxxxx600*urrrrrr+cuxxxxxxxx700*urrrrrrr+cuxxxxxxxx800*urrrrrrrr+cuxxxxxxxx010*us+cuxxxxxxxx110*urs+cuxxxxxxxx210*urrs+cuxxxxxxxx310*urrrs+cuxxxxxxxx410*urrrrs+cuxxxxxxxx510*urrrrrs+cuxxxxxxxx610*urrrrrrs+cuxxxxxxxx710*urrrrrrrs+cuxxxxxxxx020*uss+cuxxxxxxxx120*urss+cuxxxxxxxx220*urrss+cuxxxxxxxx320*urrrss+cuxxxxxxxx420*urrrrss+cuxxxxxxxx520*urrrrrss+cuxxxxxxxx620*urrrrrrss+cuxxxxxxxx030*usss+cuxxxxxxxx130*ursss+cuxxxxxxxx230*urrsss+cuxxxxxxxx330*urrrsss+cuxxxxxxxx430*urrrrsss+cuxxxxxxxx530*urrrrrsss+cuxxxxxxxx040*ussss+cuxxxxxxxx140*urssss+cuxxxxxxxx240*urrssss+cuxxxxxxxx340*urrrssss+cuxxxxxxxx440*urrrrssss+cuxxxxxxxx050*usssss+cuxxxxxxxx150*ursssss+cuxxxxxxxx250*urrsssss+cuxxxxxxxx350*urrrsssss+cuxxxxxxxx060*ussssss+cuxxxxxxxx160*urssssss+cuxxxxxxxx260*urrssssss+cuxxxxxxxx070*usssssss+cuxxxxxxxx170*ursssssss+cuxxxxxxxx080*ussssssss
                       ! ! uxxxxxxxy = cuxxxxxxxy100*ur+cuxxxxxxxy200*urr+cuxxxxxxxy300*urrr+cuxxxxxxxy400*urrrr+cuxxxxxxxy500*urrrrr+cuxxxxxxxy600*urrrrrr+cuxxxxxxxy700*urrrrrrr+cuxxxxxxxy800*urrrrrrrr+cuxxxxxxxy010*us+cuxxxxxxxy110*urs+cuxxxxxxxy210*urrs+cuxxxxxxxy310*urrrs+cuxxxxxxxy410*urrrrs+cuxxxxxxxy510*urrrrrs+cuxxxxxxxy610*urrrrrrs+cuxxxxxxxy710*urrrrrrrs+cuxxxxxxxy020*uss+cuxxxxxxxy120*urss+cuxxxxxxxy220*urrss+cuxxxxxxxy320*urrrss+cuxxxxxxxy420*urrrrss+cuxxxxxxxy520*urrrrrss+cuxxxxxxxy620*urrrrrrss+cuxxxxxxxy030*usss+cuxxxxxxxy130*ursss+cuxxxxxxxy230*urrsss+cuxxxxxxxy330*urrrsss+cuxxxxxxxy430*urrrrsss+cuxxxxxxxy530*urrrrrsss+cuxxxxxxxy040*ussss+cuxxxxxxxy140*urssss+cuxxxxxxxy240*urrssss+cuxxxxxxxy340*urrrssss+cuxxxxxxxy440*urrrrssss+cuxxxxxxxy050*usssss+cuxxxxxxxy150*ursssss+cuxxxxxxxy250*urrsssss+cuxxxxxxxy350*urrrsssss+cuxxxxxxxy060*ussssss+cuxxxxxxxy160*urssssss+cuxxxxxxxy260*urrssssss+cuxxxxxxxy070*usssssss+cuxxxxxxxy170*ursssssss+cuxxxxxxxy080*ussssssss
                       ! uxxxxxxyy = cuxxxxxxyy100*ur+cuxxxxxxyy200*urr+cuxxxxxxyy300*urrr+cuxxxxxxyy400*urrrr+cuxxxxxxyy500*urrrrr+cuxxxxxxyy600*urrrrrr+cuxxxxxxyy700*urrrrrrr+cuxxxxxxyy800*urrrrrrrr+cuxxxxxxyy010*us+cuxxxxxxyy110*urs+cuxxxxxxyy210*urrs+cuxxxxxxyy310*urrrs+cuxxxxxxyy410*urrrrs+cuxxxxxxyy510*urrrrrs+cuxxxxxxyy610*urrrrrrs+cuxxxxxxyy710*urrrrrrrs+cuxxxxxxyy020*uss+cuxxxxxxyy120*urss+cuxxxxxxyy220*urrss+cuxxxxxxyy320*urrrss+cuxxxxxxyy420*urrrrss+cuxxxxxxyy520*urrrrrss+cuxxxxxxyy620*urrrrrrss+cuxxxxxxyy030*usss+cuxxxxxxyy130*ursss+cuxxxxxxyy230*urrsss+cuxxxxxxyy330*urrrsss+cuxxxxxxyy430*urrrrsss+cuxxxxxxyy530*urrrrrsss+cuxxxxxxyy040*ussss+cuxxxxxxyy140*urssss+cuxxxxxxyy240*urrssss+cuxxxxxxyy340*urrrssss+cuxxxxxxyy440*urrrrssss+cuxxxxxxyy050*usssss+cuxxxxxxyy150*ursssss+cuxxxxxxyy250*urrsssss+cuxxxxxxyy350*urrrsssss+cuxxxxxxyy060*ussssss+cuxxxxxxyy160*urssssss+cuxxxxxxyy260*urrssssss+cuxxxxxxyy070*usssssss+cuxxxxxxyy170*ursssssss+cuxxxxxxyy080*ussssssss
                       ! ! uxxxxxyyy = cuxxxxxyyy100*ur+cuxxxxxyyy200*urr+cuxxxxxyyy300*urrr+cuxxxxxyyy400*urrrr+cuxxxxxyyy500*urrrrr+cuxxxxxyyy600*urrrrrr+cuxxxxxyyy700*urrrrrrr+cuxxxxxyyy800*urrrrrrrr+cuxxxxxyyy010*us+cuxxxxxyyy110*urs+cuxxxxxyyy210*urrs+cuxxxxxyyy310*urrrs+cuxxxxxyyy410*urrrrs+cuxxxxxyyy510*urrrrrs+cuxxxxxyyy610*urrrrrrs+cuxxxxxyyy710*urrrrrrrs+cuxxxxxyyy020*uss+cuxxxxxyyy120*urss+cuxxxxxyyy220*urrss+cuxxxxxyyy320*urrrss+cuxxxxxyyy420*urrrrss+cuxxxxxyyy520*urrrrrss+cuxxxxxyyy620*urrrrrrss+cuxxxxxyyy030*usss+cuxxxxxyyy130*ursss+cuxxxxxyyy230*urrsss+cuxxxxxyyy330*urrrsss+cuxxxxxyyy430*urrrrsss+cuxxxxxyyy530*urrrrrsss+cuxxxxxyyy040*ussss+cuxxxxxyyy140*urssss+cuxxxxxyyy240*urrssss+cuxxxxxyyy340*urrrssss+cuxxxxxyyy440*urrrrssss+cuxxxxxyyy050*usssss+cuxxxxxyyy150*ursssss+cuxxxxxyyy250*urrsssss+cuxxxxxyyy350*urrrsssss+cuxxxxxyyy060*ussssss+cuxxxxxyyy160*urssssss+cuxxxxxyyy260*urrssssss+cuxxxxxyyy070*usssssss+cuxxxxxyyy170*ursssssss+cuxxxxxyyy080*ussssssss
                       ! uxxxxyyyy = cuxxxxyyyy100*ur+cuxxxxyyyy200*urr+cuxxxxyyyy300*urrr+cuxxxxyyyy400*urrrr+cuxxxxyyyy500*urrrrr+cuxxxxyyyy600*urrrrrr+cuxxxxyyyy700*urrrrrrr+cuxxxxyyyy800*urrrrrrrr+cuxxxxyyyy010*us+cuxxxxyyyy110*urs+cuxxxxyyyy210*urrs+cuxxxxyyyy310*urrrs+cuxxxxyyyy410*urrrrs+cuxxxxyyyy510*urrrrrs+cuxxxxyyyy610*urrrrrrs+cuxxxxyyyy710*urrrrrrrs+cuxxxxyyyy020*uss+cuxxxxyyyy120*urss+cuxxxxyyyy220*urrss+cuxxxxyyyy320*urrrss+cuxxxxyyyy420*urrrrss+cuxxxxyyyy520*urrrrrss+cuxxxxyyyy620*urrrrrrss+cuxxxxyyyy030*usss+cuxxxxyyyy130*ursss+cuxxxxyyyy230*urrsss+cuxxxxyyyy330*urrrsss+cuxxxxyyyy430*urrrrsss+cuxxxxyyyy530*urrrrrsss+cuxxxxyyyy040*ussss+cuxxxxyyyy140*urssss+cuxxxxyyyy240*urrssss+cuxxxxyyyy340*urrrssss+cuxxxxyyyy440*urrrrssss+cuxxxxyyyy050*usssss+cuxxxxyyyy150*ursssss+cuxxxxyyyy250*urrsssss+cuxxxxyyyy350*urrrsssss+cuxxxxyyyy060*ussssss+cuxxxxyyyy160*urssssss+cuxxxxyyyy260*urrssssss+cuxxxxyyyy070*usssssss+cuxxxxyyyy170*ursssssss+cuxxxxyyyy080*ussssssss
                       ! ! uxxxyyyyy = cuxxxyyyyy100*ur+cuxxxyyyyy200*urr+cuxxxyyyyy300*urrr+cuxxxyyyyy400*urrrr+cuxxxyyyyy500*urrrrr+cuxxxyyyyy600*urrrrrr+cuxxxyyyyy700*urrrrrrr+cuxxxyyyyy800*urrrrrrrr+cuxxxyyyyy010*us+cuxxxyyyyy110*urs+cuxxxyyyyy210*urrs+cuxxxyyyyy310*urrrs+cuxxxyyyyy410*urrrrs+cuxxxyyyyy510*urrrrrs+cuxxxyyyyy610*urrrrrrs+cuxxxyyyyy710*urrrrrrrs+cuxxxyyyyy020*uss+cuxxxyyyyy120*urss+cuxxxyyyyy220*urrss+cuxxxyyyyy320*urrrss+cuxxxyyyyy420*urrrrss+cuxxxyyyyy520*urrrrrss+cuxxxyyyyy620*urrrrrrss+cuxxxyyyyy030*usss+cuxxxyyyyy130*ursss+cuxxxyyyyy230*urrsss+cuxxxyyyyy330*urrrsss+cuxxxyyyyy430*urrrrsss+cuxxxyyyyy530*urrrrrsss+cuxxxyyyyy040*ussss+cuxxxyyyyy140*urssss+cuxxxyyyyy240*urrssss+cuxxxyyyyy340*urrrssss+cuxxxyyyyy440*urrrrssss+cuxxxyyyyy050*usssss+cuxxxyyyyy150*ursssss+cuxxxyyyyy250*urrsssss+cuxxxyyyyy350*urrrsssss+cuxxxyyyyy060*ussssss+cuxxxyyyyy160*urssssss+cuxxxyyyyy260*urrssssss+cuxxxyyyyy070*usssssss+cuxxxyyyyy170*ursssssss+cuxxxyyyyy080*ussssssss
                       ! uxxyyyyyy = cuxxyyyyyy100*ur+cuxxyyyyyy200*urr+cuxxyyyyyy300*urrr+cuxxyyyyyy400*urrrr+cuxxyyyyyy500*urrrrr+cuxxyyyyyy600*urrrrrr+cuxxyyyyyy700*urrrrrrr+cuxxyyyyyy800*urrrrrrrr+cuxxyyyyyy010*us+cuxxyyyyyy110*urs+cuxxyyyyyy210*urrs+cuxxyyyyyy310*urrrs+cuxxyyyyyy410*urrrrs+cuxxyyyyyy510*urrrrrs+cuxxyyyyyy610*urrrrrrs+cuxxyyyyyy710*urrrrrrrs+cuxxyyyyyy020*uss+cuxxyyyyyy120*urss+cuxxyyyyyy220*urrss+cuxxyyyyyy320*urrrss+cuxxyyyyyy420*urrrrss+cuxxyyyyyy520*urrrrrss+cuxxyyyyyy620*urrrrrrss+cuxxyyyyyy030*usss+cuxxyyyyyy130*ursss+cuxxyyyyyy230*urrsss+cuxxyyyyyy330*urrrsss+cuxxyyyyyy430*urrrrsss+cuxxyyyyyy530*urrrrrsss+cuxxyyyyyy040*ussss+cuxxyyyyyy140*urssss+cuxxyyyyyy240*urrssss+cuxxyyyyyy340*urrrssss+cuxxyyyyyy440*urrrrssss+cuxxyyyyyy050*usssss+cuxxyyyyyy150*ursssss+cuxxyyyyyy250*urrsssss+cuxxyyyyyy350*urrrsssss+cuxxyyyyyy060*ussssss+cuxxyyyyyy160*urssssss+cuxxyyyyyy260*urrssssss+cuxxyyyyyy070*usssssss+cuxxyyyyyy170*ursssssss+cuxxyyyyyy080*ussssssss
                       ! ! uxyyyyyyy = cuxyyyyyyy100*ur+cuxyyyyyyy200*urr+cuxyyyyyyy300*urrr+cuxyyyyyyy400*urrrr+cuxyyyyyyy500*urrrrr+cuxyyyyyyy600*urrrrrr+cuxyyyyyyy700*urrrrrrr+cuxyyyyyyy800*urrrrrrrr+cuxyyyyyyy010*us+cuxyyyyyyy110*urs+cuxyyyyyyy210*urrs+cuxyyyyyyy310*urrrs+cuxyyyyyyy410*urrrrs+cuxyyyyyyy510*urrrrrs+cuxyyyyyyy610*urrrrrrs+cuxyyyyyyy710*urrrrrrrs+cuxyyyyyyy020*uss+cuxyyyyyyy120*urss+cuxyyyyyyy220*urrss+cuxyyyyyyy320*urrrss+cuxyyyyyyy420*urrrrss+cuxyyyyyyy520*urrrrrss+cuxyyyyyyy620*urrrrrrss+cuxyyyyyyy030*usss+cuxyyyyyyy130*ursss+cuxyyyyyyy230*urrsss+cuxyyyyyyy330*urrrsss+cuxyyyyyyy430*urrrrsss+cuxyyyyyyy530*urrrrrsss+cuxyyyyyyy040*ussss+cuxyyyyyyy140*urssss+cuxyyyyyyy240*urrssss+cuxyyyyyyy340*urrrssss+cuxyyyyyyy440*urrrrssss+cuxyyyyyyy050*usssss+cuxyyyyyyy150*ursssss+cuxyyyyyyy250*urrsssss+cuxyyyyyyy350*urrrsssss+cuxyyyyyyy060*ussssss+cuxyyyyyyy160*urssssss+cuxyyyyyyy260*urrssssss+cuxyyyyyyy070*usssssss+cuxyyyyyyy170*ursssssss+cuxyyyyyyy080*ussssssss
                       ! uyyyyyyyy = cuyyyyyyyy100*ur+cuyyyyyyyy200*urr+cuyyyyyyyy300*urrr+cuyyyyyyyy400*urrrr+cuyyyyyyyy500*urrrrr+cuyyyyyyyy600*urrrrrr+cuyyyyyyyy700*urrrrrrr+cuyyyyyyyy800*urrrrrrrr+cuyyyyyyyy010*us+cuyyyyyyyy110*urs+cuyyyyyyyy210*urrs+cuyyyyyyyy310*urrrs+cuyyyyyyyy410*urrrrs+cuyyyyyyyy510*urrrrrs+cuyyyyyyyy610*urrrrrrs+cuyyyyyyyy710*urrrrrrrs+cuyyyyyyyy020*uss+cuyyyyyyyy120*urss+cuyyyyyyyy220*urrss+cuyyyyyyyy320*urrrss+cuyyyyyyyy420*urrrrss+cuyyyyyyyy520*urrrrrss+cuyyyyyyyy620*urrrrrrss+cuyyyyyyyy030*usss+cuyyyyyyyy130*ursss+cuyyyyyyyy230*urrsss+cuyyyyyyyy330*urrrsss+cuyyyyyyyy430*urrrrsss+cuyyyyyyyy530*urrrrrsss+cuyyyyyyyy040*ussss+cuyyyyyyyy140*urssss+cuyyyyyyyy240*urrssss+cuyyyyyyyy340*urrrssss+cuyyyyyyyy440*urrrrssss+cuyyyyyyyy050*usssss+cuyyyyyyyy150*ursssss+cuyyyyyyyy250*urrsssss+cuyyyyyyyy350*urrrsssss+cuyyyyyyyy060*ussssss+cuyyyyyyyy160*urssssss+cuyyyyyyyy260*urrssssss+cuyyyyyyyy070*usssssss+cuyyyyyyyy170*ursssssss+cuyyyyyyyy080*ussssssss
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
                                   !fv(m) = fv(m) -( f(i1,i2,i3,freq) + cdtSqBy12*( cSq*(fxx22(i1,i2,i3,freq) + fyy22(i1,i2,i3,freq)) - omega*omega*f(i1,i2,i3,freq)) )*coswt 
                                                      end do ! do freq  
                                                else if( addForcing.ne.0 )then  
                                                      fv(m) = f(i1,i2,i3,0)
                                                end if
                       ! Here is the ME update 
                                              un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx + uyy ) + cdtPow4By12*( uxxxx + uyyyy + 2.*uxxyy )  + cdtPow6By360*( uxxxxxx +  uyyyyyy +  3.*( uxxxxyy +  uxxyyyy) ) + dtSq*fv(m)      
                     ! if( i1.eq.5 .and. i2.eq.6 )then
                     !   call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t+dt,uc,ue )
                     !   write(*,'("ME8C: (i1,i2)=(",2i3,") u=",1pe12.4," ue=",1pe12.4," err=",1pe8.2)') i1,i2,un(i1,i2,i3,0),ue,abs(un(i1,i2,i3,0)-ue)
                     ! end if
                                          end if ! mask 
                                        end do
                                        end do
                                        end do
             ! #If 6 == 6 || 6 == 8 
             !   evalDerivativesRectangular()
             !   write(*,*) ' Stop here for now'
             !   stop 666
             ! #End
           !   ! --- TAYLOR TIME-STEPPING --- 
           !   m=0 ! component number 
           !   ec = 0 ! component number
           !   ! #If "curvilinear" eq "curvilinear"
           !   !   #If "6" eq "4" && "6" eq "4"
           !   !     computeLaplacianOrder2(2)
           !   !   #End
           !   ! #End
           !   if( forcingOption.eq.helmholtzForcing )then
           !     coswt = cos(omega*t)
           !   end if 
           !   fv(m)=0.
           !   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
           !     getForcing(2,6,6,curvilinear) 
           !    #If "6" eq "2"
           !      ! --- SECOND 6 ---
           !      #If "2" eq "2"
           !        ! --- TWO DIMENSIONS ---
           !        #If "curvilinear" eq "rectangular"
           !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx22r(i1,i2,i3,0) + uyy22r(i1,i2,i3,0) ) + dtSq*fv(m)
           !         ! write(*,'(" adv: i1,i2=",2i4," un,u,um=",3e12.2," cdtSq,fv=",2e12.2)') i1,i2,un(i1,i2,i3,m),u(i1,i2,i3,m),um(i1,i2,i3,m),cdtSq,fv(m)
           !        #Else
           !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx22(i1,i2,i3,0)  + uyy22(i1,i2,i3,0) ) + dtSq*fv(m)
           !        #End
           !      #Else
           !        ! --- THREE DIMENSIONS ---
           !        #If "curvilinear" eq "rectangular"
           !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx23r(i1,i2,i3,0) + uyy23r(i1,i2,i3,0) + uzz23r(i1,i2,i3,0) ) + dtSq*fv(m)
           !        #Else
           !         ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx23(i1,i2,i3,0)  + uyy23(i1,i2,i3,0)  + uzz23(i1,i2,i3,0)  ) + dtSq*fv(m)
           !        #End
           !      #End
           !    #Elif "6" eq "4"
           !      ! --- -FOURTH 6 ---
           !      #If "2" eq "2"
           !        ! --- FOUTH-6 TWO DIMENSIONS ---
           !        #If "6" eq "4"
           !          ! orderInSpace=4 and orderInTime=4 
           !          #If "curvilinear" eq "rectangular"
           !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + cdtsq12*lap2d2pow2(i1,i2,i3,m) + dtSq*fv(m)
           !          #Else
           !            ! v is assumed to hold Lap(u) to 2nd-order
           !            ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
           !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !                                                            + cdtsq12*( vxx22(i1,i2,i3,m) + vyy22(i1,i2,i3,m) ) + dtSq*fv(m)
           !          #End
           !        #Else
           !          ! orderInSpace==4 and orderInTime==2                                                   
           !          #If "curvilinear" eq "rectangular"
           !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + dtSq*fv(m)
           !          #Else
           !            ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End
           !        #End                                                          
           !      #Else
           !        ! --- FOURTH-6 THREE DIMENSIONS ---
           !        #If "6" eq "4"
           !          ! orderInSpace=4 and orderInTime=4 
           !          #If "curvilinear" eq "rectangular"
           !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + cdtsq12*lap3d2pow2(i1,i2,i3,m) + dtSq*fv(m)
           !          #Else
           !            ! v is assumed to hold Lap(u) to 2nd-order
           !            ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
           !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !                                                            + cdtsq12*( vxx23(i1,i2,i3,m) + vyy23(i1,i2,i3,m) + vzz23(i1,i2,i3,m) ) + dtSq*fv(m)
           !          #End
           !        #Else
           !          ! orderInSpace==4 and orderInTime==2                                                   
           !          #If "curvilinear" eq "rectangular"
           !            !un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + dtSq*fv(m)
           !          #Else
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End
           !        #End     
           !      #End
           !    #Elif "6" eq "6"
           !      ! ---- SIXTH 6 ---
           !      #If "2" eq "2"
           !        ! --- SIXTH-6 TWO DIMENSIONS ---
           !        #If "6" eq "2"
           !          ! orderInSpace==6 and orderInTime==2                                                   
           !          #If "curvilinear" eq "rectangular"
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx62r(i1,i2,i3,m) + uyy62r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #Else
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx62(i1,i2,i3,m)  +  uyy62(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End       
           !        #Else
           !          stop 6666                                                           
           ! !          ! ---- MODIFIED EQUATION 6=6 2D -----
           ! !          getSixthDerivatives2d(6,curvilinear,evalMetrics,i1,i2,i3)
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
           !        ! --- SIXTH-6 THREE DIMENSIONS ---
           !        #If "6" eq "2
           !          ! orderInSpace==6 and orderInTime==2                                                   
           !          #If "curvilinear" eq "rectangular"
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx63r(i1,i2,i3,m) + uyy63r(i1,i2,i3,m) + uzz63r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #Else
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx63(i1,i2,i3,m)  + uyy63(i1,i2,i3,m)  + uzz63(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End       
           !        #Else
           !          ! MODIFIED EQUATION 6=6 3D
           !          ! Turn off for now: 
           !          ! getSixthDerivatives3d(6,curvilinear,evalMetrics,i1,i2,i3)
           !          write(*,'("advWave: order=6, orderInTime=6 FINISH ME")')
           !          stop 7777
           !          ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) !          !                 + cdtsq*( uxx + uyy +uzz ) !          !                 + cdtPow4By12*( uxxxx + uyyyy + uzzzz + 2.*( uxxyy +uxxzz + uyyzz ) )  !          !                 + cdtPow6By360*( uxxxxxx +  uyyyyyy + uzzzzzz + 3.*(uxxxxyy + uxxyyyy + uxxxxzz + uyyyyzz + uxxzzzz + uyyzzzz ) + 6.*uxxyyzz ) !          !                 + dtSq*fv(m)
           !        #End     
           !      #End
           !    #Elif "6" eq "8"
           !      ! ---- EIGTH 6 ---
           !      #If "2" eq "2"
           !        ! --- EIGTH-6 TWO DIMENSIONS ---
           !        #If "6" eq "2"
           !          ! orderInSpace==8 and orderInTime==2                                                   
           !          #If "curvilinear" eq "rectangular"
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx82r(i1,i2,i3,m) + uyy82r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #Else
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx82(i1,i2,i3,m)  +  uyy82(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End       
           !        #Else
           !          write(*,'("advWave: order=6, orderInTime=6 FINISH ME")')
           !          stop 7777
           !          ! ! orderInSpace=4 and orderInTime=4 
           !          ! #If "curvilinear" eq "rectangular"
           !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap2d4(i1,i2,i3,m) + cdtsq12*lap2d2pow2(i1,i2,i3,m) + dtSq*fv(m)
           !          ! #Else
           !          !   ! v is assumed to hold Lap(u) to 2nd-order
           !          !   ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
           !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx42(i1,i2,i3,m) + uyy42(i1,i2,i3,m) ) !          !                                                   + cdtsq12*( vxx22(i1,i2,i3,m) + vyy22(i1,i2,i3,m) ) + dtSq*fv(m)
           !          ! #End
           !        #End                                                          
           !     #Else
           !        ! --- EIGTH-6 THREE DIMENSIONS ---
           !        #If "6" eq "2
           !          ! orderInSpace==8 and orderInTime==2                                                   
           !          #If "curvilinear" eq "rectangular"
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx83r(i1,i2,i3,m) + uyy83r(i1,i2,i3,m) + uzz83r(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #Else
           !           ! un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( uxx83(i1,i2,i3,m)  + uyy83(i1,i2,i3,m)  + uzz83(i1,i2,i3,m) ) !                                                            + dtSq*fv(m)
           !          #End       
           !        #Else
           !          write(*,'("advWave: order=6, orderInTime=6 FINISH ME")')
           !          stop 7777
           !          ! #If "curvilinear" eq "rectangular"
           !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*lap3d4(i1,i2,i3,m) + cdtsq12*lap3d2pow2(i1,i2,i3,m) + dtSq*fv(m)
           !          ! #Else
           !          !   ! v is assumed to hold Lap(u) to 2nd-order
           !          !   ! write(*,'(" i1,i2=",2i4," uxx4=",e10.2," true=",e10.2)') i1,i2,uxx42(i1,i2,i3,m),evxx(m)
           !          !   un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*  ( uxx43(i1,i2,i3,m) + uyy43(i1,i2,i3,m) + uzz43(i1,i2,i3,m) ) !          !                                                   + cdtsq12*( vxx23(i1,i2,i3,m) + vyy23(i1,i2,i3,m) + vzz23(i1,i2,i3,m) ) + dtSq*fv(m)
           !          ! #End
           !        #End     
           !      #End
           !    #Else
           !      write(*,'("advWave: UNKNOWN order=6")')
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
