! This file automatically generated from advWaveStencil.bf90 with bpp.
    subroutine advWaveStencil3dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,sc,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
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
  ! real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
    real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1) 
    real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
    real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
    real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds coeff of Laplacian for HA scheme
    real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
          real sc(1:125,nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
          real scr(1:125)
    integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
    integer bc(0:1,0:2),ierr
    integer gridIndexRange(0:1,0:2)
    real frequencyArray(0:*)
    integer ipar(0:*)
    real rpar(0:*)
  !     ---- local variables -----
    integer m1a,m1b,m2a,m2b,m3a,m3b,numGhost,nStart,nEnd,mt,ig,useMask
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
 !...........start statement function
    integer kd,m
  ! real rx,ry,rz,sx,sy,sz,tx,ty,tz
    real cdtPow2,cdtPow4By12,cdtPow6By360,cdtPow8By20160
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
      integer maxDeriv,d,uc,count,numGhost1,m1,m2,m3
   ! ====== variables for stencils ========
   ! real czm,czp,czz,cmz,cpz
      real dt2,dt4,dt6,dt8
      real dx2i,dy2i,dz2i
      real cdx2i,cdy2i,cdz2i
      integer stencilWidth,numStencilCoeff
      real dr1, dr2, dr3, dr1i, dr2i, dr3i, rx, ry, rz, sx, sy
      real sz, tx, ty, tz, diffOrder1, diffOrder2, diffOrder3, rxr, rxs, rxt, ryr
      real rys, ryt, rzr, rzs, rzt, sxr, sxs, sxt, syr, sys, syt
      real szr, szs, szt, txr, txs, txt, tyr, tys, tyt, tzr, tzs
      real tzt, rxx, ryy, rzz, sxx, syy, szz, txx, tyy, tzz
   ! *** The next include files were generated by cgWave/maple/writeStencilFiles.mpl ***
! Define variables to valuate stencil coefficients, dim=3, order=4, gridType=Curvilinear
! File generated by cgWave/maple/writeStencilFiles.mpl
integer i1m3,i1m2,i1m1,i1p1,i1p2,i1p3
integer i2m3,i2m2,i2m1,i2p1,i2p2,i2p3
integer i3m3,i3m2,i3m1,i3p1,i3p2,i3p3
real t0,t1,t2,t6,t7,t8,t9,t10,t14,t15,t17,t18,t20,t21,t22,t23,t24,t30,t31,t33,t34,t38,t42,t44,t46,t47,t48,t49,t50,t56,t57,t61,t62,t63,t66,t67,t70,t73,t80,t81,t82,t83,t93,t94,t98,t99,t100,t101,t114,t115,t116,t117,t118,t122,t123,t125,t126,t128,t129,t130,t131,t132,t139,t140,t144,t145,t146,t147,t151,t152,t154,t155,t157,t158,t160,t162,t164,t165,t167,t168,t169,t170,t171,t174,t175,t176,t177,t180,t181,t182,t183,t189,t191,t192,t193,t194,t195,t198,t199,t200,t205,t206,t209,t210,t213,t215,t216,t218,t223,t225,t227,t228,t229,t230,t231,t234,t235,t236,t237,t245,t246,t250,t252,t254,t255,t256,t257,t258,t265,t267,t268,t269,t270,t271,t276,t277,t278,t281,t282,t285,t287,t288,t290,t291,t292,t299,t300,t301,t302,t305,t306,t307,t308,t313,t315,t316,t318,t319,t320,t321,t324,t325,t336,t337,t341,t342,t343,t350,t351,t354,t356,t364,t366,t373,t374,t378,t380,t382,t383,t384,t385,t386,t391,t392,t393,t394,t400,t401,t402,t409,t410,t415,t417,t422,t423,t424,t425,t435,t439,t440,t445,t447,t454,t458,t462,t464,t466,t467,t468,t469,t470,t476,t477,t481,t484,t485,t486,t487,t490,t493,t500,t501,t502,t503,t513,t515,t517,t518,t519,t520,t521,t528,t530,t533,t534,t535,t536,t539,t540,t541,t544,t546,t547,t548,t551,t553,t554,t555,t562,t563,t564,t565,t568,t569,t570,t571,t576,t578,t579,t581,t582,t583,t584,t587,t588,t599,t600,t604,t605,t606,t613,t615,t616,t617,t627,t629,t638,t641,t644,t645,t646,t647,t650,t657,t658,t659,t660,t663,t664,t665,t666,t671,t673,t674,t676,t677,t678,t679,t682,t683,t694,t695,t699,t702,t703,t704,t705,t708,t709,t710,t711,t714,t717,t720,t721,t722,t723,t726,t727,t728,t729,t732,t735,t746,t747,t748,t749,t752,t753,t754,t755,t787,t788,t789,t790,t799,t800,t801,t808,t810,t811,t812,t821,t822,t823,t824,t827,t828,t829,t830,t872,t874,t904,t905,t909,t910,t911,t912,t922,t923,t927,t929,t931,t934,t935,t936,t937,t940,t941,t942,t943,t951,t952,t953,t960,t961,t964,t966,t973,t974,t975,t976,t984,t988,t989,t990,t991,t1000,t1001,t1002,t1009,t1011,t1012,t1013,t1022,t1023,t1024,t1025,t1028,t1029,t1030,t1031,t1076,t1077,t1078,t1079,t1111,t1127,t1128,t1133,t1135,t1142,t1150,t1152,t1182;
    ! #Include "../include/defineStencilVariables2dOrder2Curvilinear.h"
   !...........end   statement functions
   ! write(*,*) 'Inside advWaveStencil...'
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
      useMask=0  ! do this for now -- do not check mask in loops, these seems faster
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
      dt2 = dt*dt;
      dt4 = dt2*dt2;
      dt6 = dt4*dt2;
      dt8 = dt6*dt2;
      dx2i = 1./dx(0)**2
      dy2i = 1./dx(1)**2
      dz2i = 1./dx(2)**2
      cdx2i = csq/dx(0)**2
      cdy2i = csq/dx(1)**2
      cdz2i = csq/dx(2)**2  
   ! if( option.eq.1 )then 
   !  useSosupDissipation = 1
   ! else
   !  useSosupDissipation = 0
   ! end if
      if( (.true. .or. debug.gt.1) .and. t.le.dt )then
          write(*,'("advWaveStencil: option=",i4," grid=",i4)') option,grid
          write(*,'("advWaveStencil: orderOfAccuracy=",i2," orderInTime=",i2  )') orderOfAccuracy,orderInTime
          write(*,'("advWaveStencil: addForcing=",i2," forcingOption=",i2)') addForcing,forcingOption
     ! write(*,'("advWaveStencil: useUpwindDissipation=",i2,"(explicit), useImplicitUpwindDissipation=",i2," (implicit)")') useUpwindDissipation,useImplicitUpwindDissipation
          write(*,'("advWaveStencil: t,dt,c,omega=",4e10.2)') t,dt,cc,omega 
          write(*,'("advWaveStencil: gridIsImplicit=",i2," adjustOmega=",i2," solveHelmholtz=",i2)') gridIsImplicit,adjustOmega,solveHelmholtz
          if( forcingOption.eq.helmholtzForcing )then
              write(*,'("advWaveStencil: numberOfFrequencies=",i2)') numberOfFrequencies
              write(*,'("advWaveStencil: frequencyArray=",(1pe12.4,1x))') (frequencyArray(freq),freq=0,numberOfFrequencies-1)
          end if
          if( gridIsImplicit.eq.1 )then
              write(*,'("  Implicit coeff: cImp(-1:1,0) = ",3(1pe10.2,1x), "(for 2nd-order)")') cImp(-1,0),cImp(0,0),cImp(1,0)
              write(*,'("  Implicit coeff: cImp(-1:1,1) = ",3(1pe10.2,1x), "(for 4th-order)")') cImp(-1,1),cImp(0,1),cImp(1,1)
          end if
      end if
      if( forcingOption.eq.helmholtzForcing )then
     ! --- solving the Helmholtz problem ---
          if( t.le.dt .and. debug.gt.1 )then
              write(*,'("advWaveStencil: numberOfFrequencies=",i6," omega=",1pe12.4," frequencyArray(0)=",1pe12.4)') numberOfFrequencies,omega,frequencyArray(0)
          end if
          if( numberOfFrequencies.le.0 )then
              write(*,'("advWaveStencil: ERROR: numberOfFrequencies=",i6," is <= 0")') numberOfFrequencies
              stop 0123
          end if
          if( numberOfFrequencies.eq.1  .and. frequencyArray(0) .ne. omega )then
              write(*,'("advWaveStencil: ERROR: frequencyArray(0)=",1pe12.4," is not equal to omega=",1pe12.4)') frequencyArray(0),omega
              stop 1234
          end if
          if( numberOfFrequencies.gt.maxFreq )then
              write(*,'("advWaveStencil: ERROR: numberOfFrequencies > maxFreq=",i6," .. FIX ME")') maxFreq
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
        write(*,'(" advWaveStencil:ERROR: timeSteppingMethod=defaultTimeStepping -- this should be set")')
      ! '
        stop 83322
      end if
   ! -- first time through, eval lapCoeff and the stencil coefficients --
          if( lapCoeff(0,0,0,0).lt.0. )then
                    dr1=dr(0); dr1i=1./dr1;
                    dr2=dr(1); dr2i=1./dr2;
                    dr3=dr(2); dr3i=1./dr3;
          ! --- Evaluate and store coefficients in Laplacian ---
                    write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF 3D'
                    numGhost1=orderOfAccuracy/2 -1; ! check me 
                    n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                    n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                    n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                      do i3=n3a,n3b
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                        rx = rsxy(i1,i2,i3,0,0)
                        ry = rsxy(i1,i2,i3,0,1)
                        rz = rsxy(i1,i2,i3,0,2)
                        sx = rsxy(i1,i2,i3,1,0)
                        sy = rsxy(i1,i2,i3,1,1)
                        sz = rsxy(i1,i2,i3,1,2)
                        tx = rsxy(i1,i2,i3,2,0)
                        ty = rsxy(i1,i2,i3,2,1)
                        tz = rsxy(i1,i2,i3,2,2)
            ! --- choose order for (r,s,t) derivatives based on available ghost points, less accuracy is needed in ghost points  ---
                        if( (i1-4).ge.nd1a .and. (i1+4).le.nd1b )then
                            diffOrder1=min(8,4)
                        elseif( (i1-3).ge.nd1a .and. (i1+3).le.nd1b )then
                            diffOrder1=min(6,4)
                        elseif( (i1-2).ge.nd1a .and. (i1+2).le.nd1b )then
                            diffOrder1=min(4,4)
                        elseif( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
                            diffOrder1=min(2,4)
                        else
                            write(*,*) "i1,nd1a,nd1b=",i1,nd1a,nd1b
                            stop 999
                        end if
                        if( (i2-4).ge.nd2a .and. (i2+4).le.nd2b )then
                            diffOrder2=min(8,4)
                        elseif( (i2-3).ge.nd2a .and. (i2+3).le.nd2b )then
                            diffOrder2=min(6,4)
                        elseif( (i2-2).ge.nd2a .and. (i2+2).le.nd2b )then
                            diffOrder2=min(4,4)
                        elseif( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
                            diffOrder2=min(2,4)
                        else
                            write(*,*) "i2,nd2a,nd2b=",i2,nd2a,nd2b
                            stop 999
                        end if
                        if( (i3-4).ge.nd3a .and. (i3+4).le.nd3b )then
                            diffOrder3=min(8,4)
                        elseif( (i3-3).ge.nd3a .and. (i3+3).le.nd3b )then
                            diffOrder3=min(6,4)
                        elseif( (i3-2).ge.nd3a .and. (i3+2).le.nd3b )then
                            diffOrder3=min(4,4)
                        elseif( (i3-1).ge.nd3a .and. (i3+1).le.nd3b )then
                            diffOrder3=min(2,4)
                        else
                            write(*,*) "i3,nd3a,nd3b=",i3,nd3a,nd3b
                            stop 999
                        end if
                        if( diffOrder1.eq.2 )then
                            rxr = (rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0))*(.5*dr1i) 
                            ryr = (rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1))*(.5*dr1i) 
                            rzr = (rsxy(i1+1,i2,i3,0,2)-rsxy(i1-1,i2,i3,0,2))*(.5*dr1i) 
                            sxr = (rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0))*(.5*dr1i) 
                            syr = (rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1))*(.5*dr1i) 
                            szr = (rsxy(i1+1,i2,i3,1,2)-rsxy(i1-1,i2,i3,1,2))*(.5*dr1i) 
                            txr = (rsxy(i1+1,i2,i3,2,0)-rsxy(i1-1,i2,i3,2,0))*(.5*dr1i) 
                            tyr = (rsxy(i1+1,i2,i3,2,1)-rsxy(i1-1,i2,i3,2,1))*(.5*dr1i) 
                            tzr = (rsxy(i1+1,i2,i3,2,2)-rsxy(i1-1,i2,i3,2,2))*(.5*dr1i) 
                        elseif( diffOrder1.eq.4 )then
                            rxr = ( 8*(rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0)) -(rsxy(i1+2,i2,i3,0,0)-rsxy(i1-2,i2,i3,0,0)) )*(dr1i/12.) 
                            ryr = ( 8*(rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1)) -(rsxy(i1+2,i2,i3,0,1)-rsxy(i1-2,i2,i3,0,1)) )*(dr1i/12.) 
                            rzr = ( 8*(rsxy(i1+1,i2,i3,0,2)-rsxy(i1-1,i2,i3,0,2)) -(rsxy(i1+2,i2,i3,0,2)-rsxy(i1-2,i2,i3,0,2)) )*(dr1i/12.) 
                            sxr = ( 8*(rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0)) -(rsxy(i1+2,i2,i3,1,0)-rsxy(i1-2,i2,i3,1,0)) )*(dr1i/12.) 
                            syr = ( 8*(rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1)) -(rsxy(i1+2,i2,i3,1,1)-rsxy(i1-2,i2,i3,1,1)) )*(dr1i/12.) 
                            szr = ( 8*(rsxy(i1+1,i2,i3,1,2)-rsxy(i1-1,i2,i3,1,2)) -(rsxy(i1+2,i2,i3,1,2)-rsxy(i1-2,i2,i3,1,2)) )*(dr1i/12.) 
                            txr = ( 8*(rsxy(i1+1,i2,i3,2,0)-rsxy(i1-1,i2,i3,2,0)) -(rsxy(i1+2,i2,i3,2,0)-rsxy(i1-2,i2,i3,2,0)) )*(dr1i/12.) 
                            tyr = ( 8*(rsxy(i1+1,i2,i3,2,1)-rsxy(i1-1,i2,i3,2,1)) -(rsxy(i1+2,i2,i3,2,1)-rsxy(i1-2,i2,i3,2,1)) )*(dr1i/12.) 
                            tzr = ( 8*(rsxy(i1+1,i2,i3,2,2)-rsxy(i1-1,i2,i3,2,2)) -(rsxy(i1+2,i2,i3,2,2)-rsxy(i1-2,i2,i3,2,2)) )*(dr1i/12.) 
                        elseif( diffOrder1.eq.6 )then
                            rxr = ( 45.*(rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0)) -9.*(rsxy(i1+2,i2,i3,0,0)-rsxy(i1-2,i2,i3,0,0)) +(rsxy(i1+3,i2,i3,0,0)-rsxy(i1-3,i2,i3,0,0)) )*(dr1i/60.) 
                            ryr = ( 45.*(rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1)) -9.*(rsxy(i1+2,i2,i3,0,1)-rsxy(i1-2,i2,i3,0,1)) +(rsxy(i1+3,i2,i3,0,1)-rsxy(i1-3,i2,i3,0,1)) )*(dr1i/60.) 
                            rzr = ( 45.*(rsxy(i1+1,i2,i3,0,2)-rsxy(i1-1,i2,i3,0,2)) -9.*(rsxy(i1+2,i2,i3,0,2)-rsxy(i1-2,i2,i3,0,2)) +(rsxy(i1+3,i2,i3,0,2)-rsxy(i1-3,i2,i3,0,2)) )*(dr1i/60.) 
                            sxr = ( 45.*(rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0)) -9.*(rsxy(i1+2,i2,i3,1,0)-rsxy(i1-2,i2,i3,1,0)) +(rsxy(i1+3,i2,i3,1,0)-rsxy(i1-3,i2,i3,1,0)) )*(dr1i/60.) 
                            syr = ( 45.*(rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1)) -9.*(rsxy(i1+2,i2,i3,1,1)-rsxy(i1-2,i2,i3,1,1)) +(rsxy(i1+3,i2,i3,1,1)-rsxy(i1-3,i2,i3,1,1)) )*(dr1i/60.) 
                            szr = ( 45.*(rsxy(i1+1,i2,i3,1,2)-rsxy(i1-1,i2,i3,1,2)) -9.*(rsxy(i1+2,i2,i3,1,2)-rsxy(i1-2,i2,i3,1,2)) +(rsxy(i1+3,i2,i3,1,2)-rsxy(i1-3,i2,i3,1,2)) )*(dr1i/60.) 
                            txr = ( 45.*(rsxy(i1+1,i2,i3,2,0)-rsxy(i1-1,i2,i3,2,0)) -9.*(rsxy(i1+2,i2,i3,2,0)-rsxy(i1-2,i2,i3,2,0)) +(rsxy(i1+3,i2,i3,2,0)-rsxy(i1-3,i2,i3,2,0)) )*(dr1i/60.) 
                            tyr = ( 45.*(rsxy(i1+1,i2,i3,2,1)-rsxy(i1-1,i2,i3,2,1)) -9.*(rsxy(i1+2,i2,i3,2,1)-rsxy(i1-2,i2,i3,2,1)) +(rsxy(i1+3,i2,i3,2,1)-rsxy(i1-3,i2,i3,2,1)) )*(dr1i/60.) 
                            tzr = ( 45.*(rsxy(i1+1,i2,i3,2,2)-rsxy(i1-1,i2,i3,2,2)) -9.*(rsxy(i1+2,i2,i3,2,2)-rsxy(i1-2,i2,i3,2,2)) +(rsxy(i1+3,i2,i3,2,2)-rsxy(i1-3,i2,i3,2,2)) )*(dr1i/60.) 
                        elseif( diffOrder1.eq.8 )then
                            rxr = ( 672.*(rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0)) -168.*(rsxy(i1+2,i2,i3,0,0)-rsxy(i1-2,i2,i3,0,0)) +32*(rsxy(i1+3,i2,i3,0,0)-rsxy(i1-3,i2,i3,0,0)) -3.*(rsxy(i1+4,i2,i3,0,0)-rsxy(i1-4,i2,i3,0,0)) )*(dr1i/840.) 
                            ryr = ( 672.*(rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1)) -168.*(rsxy(i1+2,i2,i3,0,1)-rsxy(i1-2,i2,i3,0,1)) +32*(rsxy(i1+3,i2,i3,0,1)-rsxy(i1-3,i2,i3,0,1)) -3.*(rsxy(i1+4,i2,i3,0,1)-rsxy(i1-4,i2,i3,0,1)) )*(dr1i/840.) 
                            rzr = ( 672.*(rsxy(i1+1,i2,i3,0,2)-rsxy(i1-1,i2,i3,0,2)) -168.*(rsxy(i1+2,i2,i3,0,2)-rsxy(i1-2,i2,i3,0,2)) +32*(rsxy(i1+3,i2,i3,0,2)-rsxy(i1-3,i2,i3,0,2)) -3.*(rsxy(i1+4,i2,i3,0,2)-rsxy(i1-4,i2,i3,0,2)) )*(dr1i/840.) 
                            sxr = ( 672.*(rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0)) -168.*(rsxy(i1+2,i2,i3,1,0)-rsxy(i1-2,i2,i3,1,0)) +32*(rsxy(i1+3,i2,i3,1,0)-rsxy(i1-3,i2,i3,1,0)) -3.*(rsxy(i1+4,i2,i3,1,0)-rsxy(i1-4,i2,i3,1,0)) )*(dr1i/840.) 
                            syr = ( 672.*(rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1)) -168.*(rsxy(i1+2,i2,i3,1,1)-rsxy(i1-2,i2,i3,1,1)) +32*(rsxy(i1+3,i2,i3,1,1)-rsxy(i1-3,i2,i3,1,1)) -3.*(rsxy(i1+4,i2,i3,1,1)-rsxy(i1-4,i2,i3,1,1)) )*(dr1i/840.) 
                            szr = ( 672.*(rsxy(i1+1,i2,i3,1,2)-rsxy(i1-1,i2,i3,1,2)) -168.*(rsxy(i1+2,i2,i3,1,2)-rsxy(i1-2,i2,i3,1,2)) +32*(rsxy(i1+3,i2,i3,1,2)-rsxy(i1-3,i2,i3,1,2)) -3.*(rsxy(i1+4,i2,i3,1,2)-rsxy(i1-4,i2,i3,1,2)) )*(dr1i/840.) 
                            txr = ( 672.*(rsxy(i1+1,i2,i3,2,0)-rsxy(i1-1,i2,i3,2,0)) -168.*(rsxy(i1+2,i2,i3,2,0)-rsxy(i1-2,i2,i3,2,0)) +32*(rsxy(i1+3,i2,i3,2,0)-rsxy(i1-3,i2,i3,2,0)) -3.*(rsxy(i1+4,i2,i3,2,0)-rsxy(i1-4,i2,i3,2,0)) )*(dr1i/840.) 
                            tyr = ( 672.*(rsxy(i1+1,i2,i3,2,1)-rsxy(i1-1,i2,i3,2,1)) -168.*(rsxy(i1+2,i2,i3,2,1)-rsxy(i1-2,i2,i3,2,1)) +32*(rsxy(i1+3,i2,i3,2,1)-rsxy(i1-3,i2,i3,2,1)) -3.*(rsxy(i1+4,i2,i3,2,1)-rsxy(i1-4,i2,i3,2,1)) )*(dr1i/840.) 
                            tzr = ( 672.*(rsxy(i1+1,i2,i3,2,2)-rsxy(i1-1,i2,i3,2,2)) -168.*(rsxy(i1+2,i2,i3,2,2)-rsxy(i1-2,i2,i3,2,2)) +32*(rsxy(i1+3,i2,i3,2,2)-rsxy(i1-3,i2,i3,2,2)) -3.*(rsxy(i1+4,i2,i3,2,2)-rsxy(i1-4,i2,i3,2,2)) )*(dr1i/840.) 
                        end if
                        if( diffOrder2.eq.2 )then
                            rxs = (rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0))*(.5*dr2i) 
                            rys = (rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1))*(.5*dr2i) 
                            rzs = (rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2))*(.5*dr2i) 
                            sxs = (rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0))*(.5*dr2i) 
                            sys = (rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1))*(.5*dr2i) 
                            szs = (rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2))*(.5*dr2i) 
                            txs = (rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0))*(.5*dr2i) 
                            tys = (rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1))*(.5*dr2i) 
                            tzs = (rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2))*(.5*dr2i) 
                        elseif( diffOrder2.eq.4 )then
                            rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
                            rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
                            rzs = ( 8*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) )*(dr2i/12.) 
                            sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
                            sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
                            szs = ( 8*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) )*(dr2i/12.) 
                            txs = ( 8*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) )*(dr2i/12.) 
                            tys = ( 8*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) )*(dr2i/12.) 
                            tzs = ( 8*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) )*(dr2i/12.) 
                        elseif( diffOrder2.eq.6 )then
                            rxs = ( 45.*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -9.*(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) +(rsxy(i1,i2+3,i3,0,0)-rsxy(i1,i2-3,i3,0,0)) )*(dr2i/60.) 
                            rys = ( 45.*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -9.*(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) +(rsxy(i1,i2+3,i3,0,1)-rsxy(i1,i2-3,i3,0,1)) )*(dr2i/60.) 
                            rzs = ( 45.*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -9.*(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) +(rsxy(i1,i2+3,i3,0,2)-rsxy(i1,i2-3,i3,0,2)) )*(dr2i/60.) 
                            sxs = ( 45.*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -9.*(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) +(rsxy(i1,i2+3,i3,1,0)-rsxy(i1,i2-3,i3,1,0)) )*(dr2i/60.) 
                            sys = ( 45.*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -9.*(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) +(rsxy(i1,i2+3,i3,1,1)-rsxy(i1,i2-3,i3,1,1)) )*(dr2i/60.) 
                            szs = ( 45.*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -9.*(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) +(rsxy(i1,i2+3,i3,1,2)-rsxy(i1,i2-3,i3,1,2)) )*(dr2i/60.) 
                            txs = ( 45.*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -9.*(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) +(rsxy(i1,i2+3,i3,2,0)-rsxy(i1,i2-3,i3,2,0)) )*(dr2i/60.) 
                            tys = ( 45.*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -9.*(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) +(rsxy(i1,i2+3,i3,2,1)-rsxy(i1,i2-3,i3,2,1)) )*(dr2i/60.) 
                            tzs = ( 45.*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -9.*(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) +(rsxy(i1,i2+3,i3,2,2)-rsxy(i1,i2-3,i3,2,2)) )*(dr2i/60.) 
                        elseif( diffOrder2.eq.8 )then
                            rxs = ( 672.*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -168.*(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) +32*(rsxy(i1,i2+3,i3,0,0)-rsxy(i1,i2-3,i3,0,0)) -3.*(rsxy(i1,i2+4,i3,0,0)-rsxy(i1,i2-4,i3,0,0)) )*(dr2i/840.) 
                            rys = ( 672.*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -168.*(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) +32*(rsxy(i1,i2+3,i3,0,1)-rsxy(i1,i2-3,i3,0,1)) -3.*(rsxy(i1,i2+4,i3,0,1)-rsxy(i1,i2-4,i3,0,1)) )*(dr2i/840.) 
                            rzs = ( 672.*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -168.*(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) +32*(rsxy(i1,i2+3,i3,0,2)-rsxy(i1,i2-3,i3,0,2)) -3.*(rsxy(i1,i2+4,i3,0,2)-rsxy(i1,i2-4,i3,0,2)) )*(dr2i/840.) 
                            sxs = ( 672.*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -168.*(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) +32*(rsxy(i1,i2+3,i3,1,0)-rsxy(i1,i2-3,i3,1,0)) -3.*(rsxy(i1,i2+4,i3,1,0)-rsxy(i1,i2-4,i3,1,0)) )*(dr2i/840.) 
                            sys = ( 672.*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -168.*(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) +32*(rsxy(i1,i2+3,i3,1,1)-rsxy(i1,i2-3,i3,1,1)) -3.*(rsxy(i1,i2+4,i3,1,1)-rsxy(i1,i2-4,i3,1,1)) )*(dr2i/840.) 
                            szs = ( 672.*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -168.*(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) +32*(rsxy(i1,i2+3,i3,1,2)-rsxy(i1,i2-3,i3,1,2)) -3.*(rsxy(i1,i2+4,i3,1,2)-rsxy(i1,i2-4,i3,1,2)) )*(dr2i/840.) 
                            txs = ( 672.*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -168.*(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) +32*(rsxy(i1,i2+3,i3,2,0)-rsxy(i1,i2-3,i3,2,0)) -3.*(rsxy(i1,i2+4,i3,2,0)-rsxy(i1,i2-4,i3,2,0)) )*(dr2i/840.) 
                            tys = ( 672.*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -168.*(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) +32*(rsxy(i1,i2+3,i3,2,1)-rsxy(i1,i2-3,i3,2,1)) -3.*(rsxy(i1,i2+4,i3,2,1)-rsxy(i1,i2-4,i3,2,1)) )*(dr2i/840.) 
                            tzs = ( 672.*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -168.*(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) +32*(rsxy(i1,i2+3,i3,2,2)-rsxy(i1,i2-3,i3,2,2)) -3.*(rsxy(i1,i2+4,i3,2,2)-rsxy(i1,i2-4,i3,2,2)) )*(dr2i/840.) 
                        end if
                        if( diffOrder3.eq.2 )then
                            rxt = (rsxy(i1,i2,i3+1,0,0)-rsxy(i1,i2,i3-1,0,0))*(.5*dr3i) 
                            ryt = (rsxy(i1,i2,i3+1,0,1)-rsxy(i1,i2,i3-1,0,1))*(.5*dr3i) 
                            rzt = (rsxy(i1,i2,i3+1,0,2)-rsxy(i1,i2,i3-1,0,2))*(.5*dr3i) 
                            sxt = (rsxy(i1,i2,i3+1,1,0)-rsxy(i1,i2,i3-1,1,0))*(.5*dr3i) 
                            syt = (rsxy(i1,i2,i3+1,1,1)-rsxy(i1,i2,i3-1,1,1))*(.5*dr3i) 
                            szt = (rsxy(i1,i2,i3+1,1,2)-rsxy(i1,i2,i3-1,1,2))*(.5*dr3i) 
                            txt = (rsxy(i1,i2,i3+1,2,0)-rsxy(i1,i2,i3-1,2,0))*(.5*dr3i) 
                            tyt = (rsxy(i1,i2,i3+1,2,1)-rsxy(i1,i2,i3-1,2,1))*(.5*dr3i) 
                            tzt = (rsxy(i1,i2,i3+1,2,2)-rsxy(i1,i2,i3-1,2,2))*(.5*dr3i) 
                        elseif( diffOrder3.eq.4 )then
                            rxt = ( 8*(rsxy(i1,i2,i3+1,0,0)-rsxy(i1,i2,i3-1,0,0)) -(rsxy(i1,i2,i3+2,0,0)-rsxy(i1,i2,i3-2,0,0)) )*(dr3i/12.) 
                            ryt = ( 8*(rsxy(i1,i2,i3+1,0,1)-rsxy(i1,i2,i3-1,0,1)) -(rsxy(i1,i2,i3+2,0,1)-rsxy(i1,i2,i3-2,0,1)) )*(dr3i/12.) 
                            rzt = ( 8*(rsxy(i1,i2,i3+1,0,2)-rsxy(i1,i2,i3-1,0,2)) -(rsxy(i1,i2,i3+2,0,2)-rsxy(i1,i2,i3-2,0,2)) )*(dr3i/12.) 
                            sxt = ( 8*(rsxy(i1,i2,i3+1,1,0)-rsxy(i1,i2,i3-1,1,0)) -(rsxy(i1,i2,i3+2,1,0)-rsxy(i1,i2,i3-2,1,0)) )*(dr3i/12.) 
                            syt = ( 8*(rsxy(i1,i2,i3+1,1,1)-rsxy(i1,i2,i3-1,1,1)) -(rsxy(i1,i2,i3+2,1,1)-rsxy(i1,i2,i3-2,1,1)) )*(dr3i/12.) 
                            szt = ( 8*(rsxy(i1,i2,i3+1,1,2)-rsxy(i1,i2,i3-1,1,2)) -(rsxy(i1,i2,i3+2,1,2)-rsxy(i1,i2,i3-2,1,2)) )*(dr3i/12.) 
                            txt = ( 8*(rsxy(i1,i2,i3+1,2,0)-rsxy(i1,i2,i3-1,2,0)) -(rsxy(i1,i2,i3+2,2,0)-rsxy(i1,i2,i3-2,2,0)) )*(dr3i/12.) 
                            tyt = ( 8*(rsxy(i1,i2,i3+1,2,1)-rsxy(i1,i2,i3-1,2,1)) -(rsxy(i1,i2,i3+2,2,1)-rsxy(i1,i2,i3-2,2,1)) )*(dr3i/12.) 
                            tzt = ( 8*(rsxy(i1,i2,i3+1,2,2)-rsxy(i1,i2,i3-1,2,2)) -(rsxy(i1,i2,i3+2,2,2)-rsxy(i1,i2,i3-2,2,2)) )*(dr3i/12.) 
                        elseif( diffOrder3.eq.6 )then
                            rxt = ( 45.*(rsxy(i1,i2,i3+1,0,0)-rsxy(i1,i2,i3-1,0,0)) -9.*(rsxy(i1,i2,i3+2,0,0)-rsxy(i1,i2,i3-2,0,0)) +(rsxy(i1,i2,i3+3,0,0)-rsxy(i1,i2,i3-3,0,0)))*(dr3i/60.) 
                            ryt = ( 45.*(rsxy(i1,i2,i3+1,0,1)-rsxy(i1,i2,i3-1,0,1)) -9.*(rsxy(i1,i2,i3+2,0,1)-rsxy(i1,i2,i3-2,0,1)) +(rsxy(i1,i2,i3+3,0,1)-rsxy(i1,i2,i3-3,0,1)))*(dr3i/60.) 
                            rzt = ( 45.*(rsxy(i1,i2,i3+1,0,2)-rsxy(i1,i2,i3-1,0,2)) -9.*(rsxy(i1,i2,i3+2,0,2)-rsxy(i1,i2,i3-2,0,2)) +(rsxy(i1,i2,i3+3,0,2)-rsxy(i1,i2,i3-3,0,2)))*(dr3i/60.) 
                            sxt = ( 45.*(rsxy(i1,i2,i3+1,1,0)-rsxy(i1,i2,i3-1,1,0)) -9.*(rsxy(i1,i2,i3+2,1,0)-rsxy(i1,i2,i3-2,1,0)) +(rsxy(i1,i2,i3+3,1,0)-rsxy(i1,i2,i3-3,1,0)))*(dr3i/60.) 
                            syt = ( 45.*(rsxy(i1,i2,i3+1,1,1)-rsxy(i1,i2,i3-1,1,1)) -9.*(rsxy(i1,i2,i3+2,1,1)-rsxy(i1,i2,i3-2,1,1)) +(rsxy(i1,i2,i3+3,1,1)-rsxy(i1,i2,i3-3,1,1)))*(dr3i/60.) 
                            szt = ( 45.*(rsxy(i1,i2,i3+1,1,2)-rsxy(i1,i2,i3-1,1,2)) -9.*(rsxy(i1,i2,i3+2,1,2)-rsxy(i1,i2,i3-2,1,2)) +(rsxy(i1,i2,i3+3,1,2)-rsxy(i1,i2,i3-3,1,2)))*(dr3i/60.) 
                            txt = ( 45.*(rsxy(i1,i2,i3+1,2,0)-rsxy(i1,i2,i3-1,2,0)) -9.*(rsxy(i1,i2,i3+2,2,0)-rsxy(i1,i2,i3-2,2,0)) +(rsxy(i1,i2,i3+3,2,0)-rsxy(i1,i2,i3-3,2,0)))*(dr3i/60.) 
                            tyt = ( 45.*(rsxy(i1,i2,i3+1,2,1)-rsxy(i1,i2,i3-1,2,1)) -9.*(rsxy(i1,i2,i3+2,2,1)-rsxy(i1,i2,i3-2,2,1)) +(rsxy(i1,i2,i3+3,2,1)-rsxy(i1,i2,i3-3,2,1)))*(dr3i/60.) 
                            tzt = ( 45.*(rsxy(i1,i2,i3+1,2,2)-rsxy(i1,i2,i3-1,2,2)) -9.*(rsxy(i1,i2,i3+2,2,2)-rsxy(i1,i2,i3-2,2,2)) +(rsxy(i1,i2,i3+3,2,2)-rsxy(i1,i2,i3-3,2,2)))*(dr3i/60.) 
                        elseif( diffOrder3.eq.8 )then
                            rxt = ( 672.*(rsxy(i1,i2,i3+1,0,0)-rsxy(i1,i2,i3-1,0,0)) -168.*(rsxy(i1,i2,i3+2,0,0)-rsxy(i1,i2,i3-2,0,0)) +32*(rsxy(i1,i2,i3+3,0,0)-rsxy(i1,i2,i3-3,0,0)) -3.*(rsxy(i1,i2,i3+4,0,0)-rsxy(i1,i2,i3-4,0,0)) )*(dr3i/840.) 
                            ryt = ( 672.*(rsxy(i1,i2,i3+1,0,1)-rsxy(i1,i2,i3-1,0,1)) -168.*(rsxy(i1,i2,i3+2,0,1)-rsxy(i1,i2,i3-2,0,1)) +32*(rsxy(i1,i2,i3+3,0,1)-rsxy(i1,i2,i3-3,0,1)) -3.*(rsxy(i1,i2,i3+4,0,1)-rsxy(i1,i2,i3-4,0,1)) )*(dr3i/840.) 
                            rzt = ( 672.*(rsxy(i1,i2,i3+1,0,2)-rsxy(i1,i2,i3-1,0,2)) -168.*(rsxy(i1,i2,i3+2,0,2)-rsxy(i1,i2,i3-2,0,2)) +32*(rsxy(i1,i2,i3+3,0,2)-rsxy(i1,i2,i3-3,0,2)) -3.*(rsxy(i1,i2,i3+4,0,2)-rsxy(i1,i2,i3-4,0,2)) )*(dr3i/840.) 
                            sxt = ( 672.*(rsxy(i1,i2,i3+1,1,0)-rsxy(i1,i2,i3-1,1,0)) -168.*(rsxy(i1,i2,i3+2,1,0)-rsxy(i1,i2,i3-2,1,0)) +32*(rsxy(i1,i2,i3+3,1,0)-rsxy(i1,i2,i3-3,1,0)) -3.*(rsxy(i1,i2,i3+4,1,0)-rsxy(i1,i2,i3-4,1,0)) )*(dr3i/840.) 
                            syt = ( 672.*(rsxy(i1,i2,i3+1,1,1)-rsxy(i1,i2,i3-1,1,1)) -168.*(rsxy(i1,i2,i3+2,1,1)-rsxy(i1,i2,i3-2,1,1)) +32*(rsxy(i1,i2,i3+3,1,1)-rsxy(i1,i2,i3-3,1,1)) -3.*(rsxy(i1,i2,i3+4,1,1)-rsxy(i1,i2,i3-4,1,1)) )*(dr3i/840.) 
                            szt = ( 672.*(rsxy(i1,i2,i3+1,1,2)-rsxy(i1,i2,i3-1,1,2)) -168.*(rsxy(i1,i2,i3+2,1,2)-rsxy(i1,i2,i3-2,1,2)) +32*(rsxy(i1,i2,i3+3,1,2)-rsxy(i1,i2,i3-3,1,2)) -3.*(rsxy(i1,i2,i3+4,1,2)-rsxy(i1,i2,i3-4,1,2)) )*(dr3i/840.) 
                            txt = ( 672.*(rsxy(i1,i2,i3+1,2,0)-rsxy(i1,i2,i3-1,2,0)) -168.*(rsxy(i1,i2,i3+2,2,0)-rsxy(i1,i2,i3-2,2,0)) +32*(rsxy(i1,i2,i3+3,2,0)-rsxy(i1,i2,i3-3,2,0)) -3.*(rsxy(i1,i2,i3+4,2,0)-rsxy(i1,i2,i3-4,2,0)) )*(dr3i/840.) 
                            tyt = ( 672.*(rsxy(i1,i2,i3+1,2,1)-rsxy(i1,i2,i3-1,2,1)) -168.*(rsxy(i1,i2,i3+2,2,1)-rsxy(i1,i2,i3-2,2,1)) +32*(rsxy(i1,i2,i3+3,2,1)-rsxy(i1,i2,i3-3,2,1)) -3.*(rsxy(i1,i2,i3+4,2,1)-rsxy(i1,i2,i3-4,2,1)) )*(dr3i/840.) 
                            tzt = ( 672.*(rsxy(i1,i2,i3+1,2,2)-rsxy(i1,i2,i3-1,2,2)) -168.*(rsxy(i1,i2,i3+2,2,2)-rsxy(i1,i2,i3-2,2,2)) +32*(rsxy(i1,i2,i3+3,2,2)-rsxy(i1,i2,i3-3,2,2)) -3.*(rsxy(i1,i2,i3+4,2,2)-rsxy(i1,i2,i3-4,2,2)) )*(dr3i/840.) 
                        end if
                        rxx = rx*rxr + sx*rxs + tx*rxt
                        ryy = ry*ryr + sy*rys + ty*ryt
                        rzz = rz*rzr + sz*rzs + tz*rzt
                        sxx = rx*sxr + sx*sxs + tx*sxt
                        syy = ry*syr + sy*sys + ty*syt
                        szz = rz*szr + sz*szs + tz*szt
                        txx = rx*txr + sx*txs + tx*txt
                        tyy = ry*tyr + sy*tys + ty*tyt
                        tzz = rz*tzr + sz*tzs + tz*tzt
            ! -- Coefficients in the Laplacian (scaled)
                        lapCoeff(i1,i2,i3,0) = (rx**2 + ry**2 + rz**2 )*dr1i**2
                        lapCoeff(i1,i2,i3,1) = (sx**2 + sy**2 + sz**2 )*dr2i**2
                        lapCoeff(i1,i2,i3,2) = (tx**2 + ty**2 + tz**2 )*dr3i**2
                        lapCoeff(i1,i2,i3,3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i
                        lapCoeff(i1,i2,i3,4) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr3i
                        lapCoeff(i1,i2,i3,5) = 2.*(sx*tx + sy*ty + sz*tz )*dr2i*dr3i
                        lapCoeff(i1,i2,i3,6) = (rxx + ryy + rzz)*dr1i
                        lapCoeff(i1,i2,i3,7) = (sxx + syy + szz)*dr2i 
                        lapCoeff(i1,i2,i3,8) = (txx + tyy + tzz)*dr3i 
                      end do
                      end do
                      end do
                  write(*,*) 'EVAL STENCIL COEFF'
           ! --- THREE DIMENSIONS ---
                  numGhost1=0;
                  n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                  n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                  n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
         ! ----- START LOOPS ----
                    do i3=n3a,n3b
                    do i2=n2a,n2b
                    do i1=n1a,n1b
! Evaluate stencil coefficients, dim=3, order=4, gridType=Curvilinear
! File generated by cgWave/maple/writeStencilFiles.mpl
i1m1=i1-1; i1p1=i1+1;
i2m1=i2-1; i2p1=i2+1;
i3m1=i3-1; i3p1=i3+1;
t1 = lapCoeff(i1,i2m1,i3m1,5);
t2 = lapCoeff(i1,i2,i3,5);
t6 = lapCoeff(i1m1,i2,i3m1,5);
t7 = lapCoeff(i1,i2,i3,4);
t8 = t6 * t7;
t9 = lapCoeff(i1,i2m1,i3m1,4);
t10 = t2 * t9;
t14 = lapCoeff(i1,i2,i3,8);
t15 = lapCoeff(i1,i2,i3m1,5);
t17 = t14 * t15 / 0.8e1;
t18 = lapCoeff(i1,i2,i3,2);
t20 = t18 * t15 / 0.4e1;
t21 = lapCoeff(i1,i2m1,i3m1,8);
t22 = t21 / 0.8e1;
t23 = lapCoeff(i1,i2m1,i3m1,2);
t24 = t23 / 0.4e1;
t30 = t2 * dt2;
t31 = t30 / 0.24e2;
t33 = lapCoeff(i1p1,i2,i3m1,5);
t34 = t33 * t7;
t38 = lapCoeff(i1m1,i2,i3m1,4);
t42 = lapCoeff(i1,i2,i3m1,4);
t44 = t14 * t42 / 0.8e1;
t46 = t18 * t42 / 0.4e1;
t47 = lapCoeff(i1m1,i2,i3m1,8);
t48 = t47 / 0.8e1;
t49 = lapCoeff(i1m1,i2,i3m1,2);
t50 = t49 / 0.4e1;
t56 = t7 * dt2;
t57 = t56 / 0.24e2;
t61 = lapCoeff(i1,i2,i3m1,8);
t62 = t61 / 0.2e1;
t63 = lapCoeff(i1,i2,i3m1,2);
t66 = t61 / 0.4e1;
t67 = t63 / 0.2e1;
t70 = lapCoeff(i1,i2p1,i3m1,5);
t73 = lapCoeff(i1p1,i2,i3m1,4);
t80 = lapCoeff(i1p1,i2,i3m1,8);
t81 = t80 / 0.8e1;
t82 = lapCoeff(i1p1,i2,i3m1,2);
t83 = t82 / 0.4e1;
t93 = lapCoeff(i1,i2p1,i3m1,4);
t94 = t2 * t93;
t98 = lapCoeff(i1,i2p1,i3m1,8);
t99 = t98 / 0.8e1;
t100 = lapCoeff(i1,i2p1,i3m1,2);
t101 = t100 / 0.4e1;
t114 = lapCoeff(i1m1,i2m1,i3,5);
t115 = lapCoeff(i1,i2,i3,3);
t116 = t114 * t115;
t117 = lapCoeff(i1,i2m1,i3m1,3);
t118 = t2 * t117;
t122 = lapCoeff(i1,i2,i3,7);
t123 = lapCoeff(i1,i2m1,i3,5);
t125 = t122 * t123 / 0.8e1;
t126 = lapCoeff(i1,i2,i3,1);
t128 = t123 * t126 / 0.4e1;
t129 = lapCoeff(i1,i2m1,i3m1,7);
t130 = t129 / 0.8e1;
t131 = lapCoeff(i1,i2m1,i3m1,1);
t132 = t131 / 0.4e1;
t139 = lapCoeff(i1p1,i2m1,i3,5);
t140 = t139 * t115;
t144 = lapCoeff(i1m1,i2m1,i3,4);
t145 = t144 * t115;
t146 = lapCoeff(i1m1,i2,i3m1,3);
t147 = t7 * t146;
t151 = lapCoeff(i1m1,i2,i3,5);
t152 = lapCoeff(i1,i2,i3,6);
t154 = t151 * t152 / 0.8e1;
t155 = lapCoeff(i1,i2m1,i3,4);
t157 = t122 * t155 / 0.8e1;
t158 = lapCoeff(i1,i2,i3m1,3);
t160 = t14 * t158 / 0.8e1;
t162 = t18 * t158 / 0.4e1;
t164 = t126 * t155 / 0.4e1;
t165 = lapCoeff(i1,i2,i3,0);
t167 = t151 * t165 / 0.4e1;
t168 = lapCoeff(i1m1,i2m1,i3,8);
t169 = t168 / 0.8e1;
t170 = lapCoeff(i1m1,i2m1,i3,2);
t171 = t170 / 0.4e1;
t174 = lapCoeff(i1m1,i2,i3m1,7);
t175 = t174 / 0.8e1;
t176 = lapCoeff(i1m1,i2,i3m1,1);
t177 = t176 / 0.4e1;
t180 = lapCoeff(i1,i2m1,i3m1,6);
t181 = t180 / 0.8e1;
t182 = lapCoeff(i1,i2m1,i3m1,0);
t183 = t182 / 0.4e1;
t189 = 0.5e1 / 0.12e2 * t30;
t191 = t2 * t165 / 0.2e1;
t192 = lapCoeff(i1,i2m1,i3,8);
t193 = t192 / 0.2e1;
t194 = lapCoeff(i1,i2m1,i3,2);
t195 = t2 / 0.2e1;
t198 = lapCoeff(i1,i2,i3m1,7);
t199 = t198 / 0.2e1;
t200 = lapCoeff(i1,i2,i3m1,1);
t205 = t192 / 0.4e1;
t206 = t194 / 0.2e1;
t209 = t198 / 0.4e1;
t210 = t200 / 0.2e1;
t213 = lapCoeff(i1p1,i2m1,i3,4);
t215 = t115 * (t144 + t213) / 0.16e2;
t216 = lapCoeff(i1p1,i2,i3m1,3);
t218 = t7 * (t146 + t216) / 0.16e2;
t223 = lapCoeff(i1p1,i2,i3,5);
t225 = t223 * t152 / 0.8e1;
t227 = t223 * t165 / 0.4e1;
t228 = lapCoeff(i1p1,i2m1,i3,8);
t229 = t228 / 0.8e1;
t230 = lapCoeff(i1p1,i2m1,i3,2);
t231 = t230 / 0.4e1;
t234 = lapCoeff(i1p1,i2,i3m1,7);
t235 = t234 / 0.8e1;
t236 = lapCoeff(i1p1,i2,i3m1,1);
t237 = t236 / 0.4e1;
t245 = t213 * t115;
t246 = t7 * t216;
t250 = lapCoeff(i1m1,i2,i3,4);
t252 = t152 * t250 / 0.8e1;
t254 = t250 * t165 / 0.4e1;
t255 = lapCoeff(i1m1,i2,i3m1,6);
t256 = t255 / 0.8e1;
t257 = lapCoeff(i1m1,i2,i3m1,0);
t258 = t257 / 0.4e1;
t265 = 0.5e1 / 0.12e2 * t56;
t267 = t126 * t7 / 0.2e1;
t268 = lapCoeff(i1m1,i2,i3,8);
t269 = t268 / 0.2e1;
t270 = lapCoeff(i1m1,i2,i3,2);
t271 = t7 / 0.2e1;
t276 = lapCoeff(i1,i2,i3m1,6);
t277 = t276 / 0.2e1;
t278 = lapCoeff(i1,i2,i3m1,0);
t281 = t268 / 0.4e1;
t282 = t270 / 0.2e1;
t285 = lapCoeff(i1m1,i2p1,i3,5);
t287 = t115 * (t114 + t285) / 0.16e2;
t288 = lapCoeff(i1,i2p1,i3m1,3);
t290 = t2 * (t117 + t288) / 0.16e2;
t291 = t276 / 0.4e1;
t292 = t278 / 0.2e1;
t299 = lapCoeff(i1,i2p1,i3m1,7);
t300 = t299 / 0.8e1;
t301 = lapCoeff(i1,i2p1,i3m1,1);
t302 = t301 / 0.4e1;
t305 = lapCoeff(i1p1,i2,i3m1,6);
t306 = t305 / 0.8e1;
t307 = lapCoeff(i1p1,i2,i3m1,0);
t308 = t307 / 0.4e1;
t313 = lapCoeff(i1,i2p1,i3,5);
t315 = t122 * (t123 + t313) / 0.8e1;
t316 = lapCoeff(i1p1,i2,i3,4);
t318 = t152 * (t250 + t316) / 0.8e1;
t319 = 0.2e1 * t18;
t320 = t123 / 0.4e1;
t321 = t313 / 0.4e1;
t324 = t250 / 0.4e1;
t325 = t316 / 0.4e1;
t336 = 0.2e1 / 0.3e1 * t14;
t337 = 0.4e1 / 0.3e1 * t18;
t341 = lapCoeff(i1p1,i2,i3,2);
t342 = lapCoeff(i1p1,i2,i3,8);
t343 = t342 / 0.2e1;
t350 = t342 / 0.4e1;
t351 = t341 / 0.2e1;
t354 = lapCoeff(i1p1,i2p1,i3,5);
t356 = t115 * (t139 + t354) / 0.16e2;
t364 = t152 * t316 / 0.8e1;
t366 = t316 * t165 / 0.4e1;
t373 = lapCoeff(i1m1,i2p1,i3,4);
t374 = t373 * t115;
t378 = lapCoeff(i1,i2p1,i3,4);
t380 = t122 * t378 / 0.8e1;
t382 = t126 * t378 / 0.4e1;
t383 = lapCoeff(i1m1,i2p1,i3,8);
t384 = t383 / 0.8e1;
t385 = lapCoeff(i1m1,i2p1,i3,2);
t386 = t385 / 0.4e1;
t391 = lapCoeff(i1,i2p1,i3m1,6);
t392 = t391 / 0.8e1;
t393 = lapCoeff(i1,i2p1,i3m1,0);
t394 = t393 / 0.4e1;
t400 = lapCoeff(i1,i2p1,i3,2);
t401 = lapCoeff(i1,i2p1,i3,8);
t402 = t401 / 0.2e1;
t409 = t401 / 0.4e1;
t410 = t400 / 0.2e1;
t415 = lapCoeff(i1p1,i2p1,i3,4);
t417 = t115 * (t373 + t415) / 0.16e2;
t422 = lapCoeff(i1p1,i2p1,i3,8);
t423 = t422 / 0.8e1;
t424 = lapCoeff(i1p1,i2p1,i3,2);
t425 = t424 / 0.4e1;
t435 = t415 * t115;
t439 = t285 * t115;
t440 = t2 * t288;
t445 = t122 * t313 / 0.8e1;
t447 = t313 * t126 / 0.4e1;
t454 = t354 * t115;
t458 = lapCoeff(i1m1,i2m1,i3,3);
t462 = lapCoeff(i1,i2m1,i3,3);
t464 = t122 * t462 / 0.8e1;
t466 = t126 * t462 / 0.4e1;
t467 = lapCoeff(i1m1,i2m1,i3,7);
t468 = t467 / 0.8e1;
t469 = lapCoeff(i1m1,i2m1,i3,1);
t470 = t469 / 0.4e1;
t476 = t115 * dt2;
t477 = t476 / 0.24e2;
t481 = lapCoeff(i1,i2m1,i3p1,5);
t484 = lapCoeff(i1,i2m1,i3,7);
t485 = t484 / 0.4e1;
t486 = lapCoeff(i1,i2m1,i3,1);
t487 = t486 / 0.2e1;
t490 = t484 / 0.2e1;
t493 = lapCoeff(i1p1,i2m1,i3,3);
t500 = lapCoeff(i1p1,i2m1,i3,7);
t501 = t500 / 0.8e1;
t502 = lapCoeff(i1p1,i2m1,i3,1);
t503 = t502 / 0.4e1;
t513 = lapCoeff(i1m1,i2,i3,3);
t515 = t152 * t513 / 0.8e1;
t517 = t513 * t165 / 0.4e1;
t518 = lapCoeff(i1m1,i2m1,i3,6);
t519 = t518 / 0.8e1;
t520 = lapCoeff(i1m1,i2m1,i3,0);
t521 = t520 / 0.4e1;
t528 = 0.5e1 / 0.12e2 * t476;
t530 = t18 * t115 / 0.2e1;
t533 = lapCoeff(i1m1,i2,i3,7);
t534 = t533 / 0.2e1;
t535 = lapCoeff(i1m1,i2,i3,1);
t536 = t115 / 0.2e1;
t539 = lapCoeff(i1,i2m1,i3,6);
t540 = t539 / 0.2e1;
t541 = lapCoeff(i1,i2m1,i3,0);
t544 = lapCoeff(i1m1,i2,i3p1,5);
t546 = t7 * (t6 + t544) / 0.16e2;
t547 = t533 / 0.4e1;
t548 = t535 / 0.2e1;
t551 = lapCoeff(i1,i2m1,i3p1,4);
t553 = t2 * (t9 + t551) / 0.16e2;
t554 = t539 / 0.4e1;
t555 = t541 / 0.2e1;
t562 = lapCoeff(i1,i2m1,i3p1,8);
t563 = t562 / 0.8e1;
t564 = lapCoeff(i1,i2m1,i3p1,2);
t565 = t564 / 0.4e1;
t568 = lapCoeff(i1p1,i2m1,i3,6);
t569 = t568 / 0.8e1;
t570 = lapCoeff(i1p1,i2m1,i3,0);
t571 = t570 / 0.4e1;
t576 = lapCoeff(i1,i2,i3p1,5);
t578 = t14 * (t15 + t576) / 0.8e1;
t579 = lapCoeff(i1p1,i2,i3,3);
t581 = t152 * (t513 + t579) / 0.8e1;
t582 = t15 / 0.4e1;
t583 = t576 / 0.4e1;
t584 = 0.2e1 * t126;
t587 = t513 / 0.4e1;
t588 = t579 / 0.4e1;
t599 = 0.2e1 / 0.3e1 * t122;
t600 = 0.4e1 / 0.3e1 * t126;
t604 = lapCoeff(i1p1,i2,i3,1);
t605 = lapCoeff(i1p1,i2,i3,7);
t606 = t605 / 0.2e1;
t613 = lapCoeff(i1p1,i2,i3p1,5);
t615 = t7 * (t33 + t613) / 0.16e2;
t616 = t605 / 0.4e1;
t617 = t604 / 0.2e1;
t627 = t152 * t579 / 0.8e1;
t629 = t579 * t165 / 0.4e1;
t638 = lapCoeff(i1m1,i2,i3p1,4);
t641 = lapCoeff(i1m1,i2p1,i3,3);
t644 = lapCoeff(i1m1,i2,i3,6);
t645 = t644 / 0.4e1;
t646 = lapCoeff(i1m1,i2,i3,0);
t647 = t646 / 0.2e1;
t650 = t644 / 0.2e1;
t657 = lapCoeff(i1m1,i2,i3p1,8);
t658 = t657 / 0.8e1;
t659 = lapCoeff(i1m1,i2,i3p1,2);
t660 = t659 / 0.4e1;
t663 = lapCoeff(i1m1,i2p1,i3,7);
t664 = t663 / 0.8e1;
t665 = lapCoeff(i1m1,i2p1,i3,1);
t666 = t665 / 0.4e1;
t671 = lapCoeff(i1,i2,i3p1,4);
t673 = t14 * (t42 + t671) / 0.8e1;
t674 = lapCoeff(i1,i2p1,i3,3);
t676 = t122 * (t462 + t674) / 0.8e1;
t677 = t42 / 0.4e1;
t678 = t671 / 0.4e1;
t679 = 0.2e1 * t165;
t682 = t462 / 0.4e1;
t683 = t674 / 0.4e1;
t694 = 0.2e1 / 0.3e1 * t152;
t695 = 0.4e1 / 0.3e1 * t165;
t699 = lapCoeff(i1,i2p1,i3p1,5);
t702 = lapCoeff(i1,i2,i3p1,8);
t703 = t702 / 0.4e1;
t704 = lapCoeff(i1,i2,i3p1,2);
t705 = t704 / 0.2e1;
t708 = lapCoeff(i1,i2p1,i3,7);
t709 = t708 / 0.4e1;
t710 = lapCoeff(i1,i2p1,i3,1);
t711 = t710 / 0.2e1;
t714 = lapCoeff(i1p1,i2,i3p1,4);
t717 = lapCoeff(i1p1,i2p1,i3,3);
t720 = lapCoeff(i1p1,i2,i3,6);
t721 = t720 / 0.4e1;
t722 = lapCoeff(i1p1,i2,i3,0);
t723 = t722 / 0.2e1;
t726 = t702 / 0.2e1;
t727 = 0.4e1 * t18;
t728 = 0.4e1 * t126;
t729 = 0.4e1 * t165;
t732 = t708 / 0.2e1;
t735 = t720 / 0.2e1;
t746 = lapCoeff(i1p1,i2,i3p1,8);
t747 = t746 / 0.8e1;
t748 = lapCoeff(i1p1,i2,i3p1,2);
t749 = t748 / 0.4e1;
t752 = lapCoeff(i1p1,i2p1,i3,7);
t753 = t752 / 0.8e1;
t754 = lapCoeff(i1p1,i2p1,i3,1);
t755 = t754 / 0.4e1;
t787 = lapCoeff(i1m1,i2p1,i3,6);
t788 = t787 / 0.8e1;
t789 = lapCoeff(i1m1,i2p1,i3,0);
t790 = t789 / 0.4e1;
t799 = lapCoeff(i1,i2p1,i3,6);
t800 = t799 / 0.2e1;
t801 = lapCoeff(i1,i2p1,i3,0);
t808 = lapCoeff(i1,i2p1,i3p1,4);
t810 = t2 * (t93 + t808) / 0.16e2;
t811 = t799 / 0.4e1;
t812 = t801 / 0.2e1;
t821 = lapCoeff(i1,i2p1,i3p1,8);
t822 = t821 / 0.8e1;
t823 = lapCoeff(i1,i2p1,i3p1,2);
t824 = t823 / 0.4e1;
t827 = lapCoeff(i1p1,i2p1,i3,6);
t828 = t827 / 0.8e1;
t829 = lapCoeff(i1p1,i2p1,i3,0);
t830 = t829 / 0.4e1;
t872 = t122 * t674 / 0.8e1;
t874 = t126 * t674 / 0.4e1;
t904 = lapCoeff(i1,i2m1,i3p1,3);
t905 = t2 * t904;
t909 = lapCoeff(i1,i2m1,i3p1,7);
t910 = t909 / 0.8e1;
t911 = lapCoeff(i1,i2m1,i3p1,1);
t912 = t911 / 0.4e1;
t922 = lapCoeff(i1m1,i2,i3p1,3);
t923 = t7 * t922;
t927 = lapCoeff(i1,i2,i3p1,3);
t929 = t14 * t927 / 0.8e1;
t931 = t18 * t927 / 0.4e1;
t934 = lapCoeff(i1m1,i2,i3p1,7);
t935 = t934 / 0.8e1;
t936 = lapCoeff(i1m1,i2,i3p1,1);
t937 = t936 / 0.4e1;
t940 = lapCoeff(i1,i2m1,i3p1,6);
t941 = t940 / 0.8e1;
t942 = lapCoeff(i1,i2m1,i3p1,0);
t943 = t942 / 0.4e1;
t951 = lapCoeff(i1,i2,i3p1,7);
t952 = t951 / 0.2e1;
t953 = lapCoeff(i1,i2,i3p1,1);
t960 = t951 / 0.4e1;
t961 = t953 / 0.2e1;
t964 = lapCoeff(i1p1,i2,i3p1,3);
t966 = t7 * (t922 + t964) / 0.16e2;
t973 = lapCoeff(i1p1,i2,i3p1,7);
t974 = t973 / 0.8e1;
t975 = lapCoeff(i1p1,i2,i3p1,1);
t976 = t975 / 0.4e1;
t984 = t7 * t964;
t988 = lapCoeff(i1m1,i2,i3p1,6);
t989 = t988 / 0.8e1;
t990 = lapCoeff(i1m1,i2,i3p1,0);
t991 = t990 / 0.4e1;
t1000 = lapCoeff(i1,i2,i3p1,6);
t1001 = t1000 / 0.2e1;
t1002 = lapCoeff(i1,i2,i3p1,0);
t1009 = lapCoeff(i1,i2p1,i3p1,3);
t1011 = t2 * (t904 + t1009) / 0.16e2;
t1012 = t1000 / 0.4e1;
t1013 = t1002 / 0.2e1;
t1022 = lapCoeff(i1,i2p1,i3p1,7);
t1023 = t1022 / 0.8e1;
t1024 = lapCoeff(i1,i2p1,i3p1,1);
t1025 = t1024 / 0.4e1;
t1028 = lapCoeff(i1p1,i2,i3p1,6);
t1029 = t1028 / 0.8e1;
t1030 = lapCoeff(i1p1,i2,i3p1,0);
t1031 = t1030 / 0.4e1;
t1076 = lapCoeff(i1,i2p1,i3p1,6);
t1077 = t1076 / 0.8e1;
t1078 = lapCoeff(i1,i2p1,i3p1,0);
t1079 = t1078 / 0.4e1;
t1111 = t2 * t1009;
t1127 = t544 * t7;
t1128 = t2 * t551;
t1133 = t14 * t576 / 0.8e1;
t1135 = t18 * t576 / 0.4e1;
t1142 = t613 * t7;
t1150 = t14 * t671 / 0.8e1;
t1152 = t18 * t671 / 0.4e1;
t1182 = t2 * t808;
sc(1,i1,i2,i3) = 0;
sc(2,i1,i2,i3) = 0;
sc(3,i1,i2,i3) = (t1 * t2 * dt4 / 0.192e3);
sc(4,i1,i2,i3) = 0;
sc(5,i1,i2,i3) = 0;
sc(6,i1,i2,i3) = 0;
sc(7,i1,i2,i3) = (dt4 * (t8 + t10) / 0.192e3);
sc(8,i1,i2,i3) = (-dt4 * (t17 - t20 + t2 * (t22 - t24)) / 0.12e2 - t31);
sc(9,i1,i2,i3) = -(dt4 * (t10 + t34) / 0.192e3);
sc(10,i1,i2,i3) = 0;
sc(11,i1,i2,i3) = (t38 * t7 * dt4 / 0.192e3);
sc(12,i1,i2,i3) = (-dt4 * (t44 - t46 + t7 * (t48 - t50)) / 0.12e2 - t57);
sc(13,i1,i2,i3) = (dt2 * (t14 - t18) / 0.12e2 - dt4 * (t18 * (t62 - t63) - t14 * (t66 - t67) + t2 * (t1 + t70) / 0.16e2 + t7 * (t38 + t73) / 0.16e2) / 0.12e2);
sc(14,i1,i2,i3) = (dt4 * (t44 - t46 + t7 * (t81 - t83)) / 0.12e2 + t57);
sc(15,i1,i2,i3) = (t7 * t73 * dt4 / 0.192e3);
sc(16,i1,i2,i3) = 0;
sc(17,i1,i2,i3) = -(dt4 * (t8 + t94) / 0.192e3);
sc(18,i1,i2,i3) = (dt4 * (t17 - t20 + t2 * (t99 - t101)) / 0.12e2 + t31);
sc(19,i1,i2,i3) = (dt4 * (t94 + t34) / 0.192e3);
sc(20,i1,i2,i3) = 0;
sc(21,i1,i2,i3) = 0;
sc(22,i1,i2,i3) = 0;
sc(23,i1,i2,i3) = (t2 * t70 * dt4 / 0.192e3);
sc(24,i1,i2,i3) = 0;
sc(25,i1,i2,i3) = 0;
sc(26,i1,i2,i3) = 0;
sc(27,i1,i2,i3) = (dt4 * (t116 + t118) / 0.192e3);
sc(28,i1,i2,i3) = (-dt4 * (t125 - t128 + t2 * (t130 - t132)) / 0.12e2 - t31);
sc(29,i1,i2,i3) = -(dt4 * (t118 + t140) / 0.192e3);
sc(30,i1,i2,i3) = 0;
sc(31,i1,i2,i3) = (dt4 * (t145 + t147) / 0.192e3);
sc(32,i1,i2,i3) = -(dt4 * (t154 + t157 + t160 - t162 - t164 - t167 + t115 * (t169 - t171) + t7 * (t175 - t177) + t2 * (t181 - t183)) / 0.12e2);
sc(33,i1,i2,i3) = (t189 - dt4 * (t191 + t126 * (t193 - t194 + t195) + t18 * (t199 + t195 - t200) + t2 * (t23 + t131 + t182) / 0.2e1 - t122 * (t205 - t206) - t14 * (t209 - t210) + t215 + t218) / 0.12e2);
sc(34,i1,i2,i3) = (dt4 * (t157 + t160 + t225 - t162 - t164 + t227 + t115 * (t229 - t231) + t7 * (t235 - t237) + t2 * (t181 + t183)) / 0.12e2);
sc(35,i1,i2,i3) = (dt4 * (t245 + t246) / 0.192e3);
sc(36,i1,i2,i3) = (-dt4 * (t252 - t254 + t7 * (t256 - t258)) / 0.12e2 - t57);
sc(37,i1,i2,i3) = (t265 - dt4 * (t267 + t165 * (t269 - t270 + t271) + t7 * (t49 + t176 + t257) / 0.2e1 + t18 * (t277 + t271 - t278) - t152 * (t281 - t282) + t287 + t290 - t14 * (t291 - t292)) / 0.12e2);
sc(38,i1,i2,i3) = (dt4 * (t2 * (t130 + t300 + t132 - t302) + t7 * (t256 + t306 + t258 - t308) + t14 * (t63 + t200 + t278) + t315 + t318 + t126 * (t14 - t319 - t320 + t321) + t165 * (t14 - t319 - t324 + t325) - t18 * (0.2e1 * t63 - t14 + t319 + 0.2e1 * t200 + 0.2e1 * t278)) / 0.12e2 - dt2 * (t336 - t337));
sc(39,i1,i2,i3) = (dt4 * (t165 * (t341 - t343 + t271) + t18 * (t277 + t271 + t278) + t267 + t7 * (t82 + t236 + t307) / 0.2e1 - t152 * (t350 - t351) + t356 + t290 - t14 * (t291 + t292)) / 0.12e2 - t265);
sc(40,i1,i2,i3) = (t57 - dt4 * (t364 + t366 + t7 * (t306 + t308)) / 0.12e2);
sc(41,i1,i2,i3) = -(dt4 * (t147 + t374) / 0.192e3);
sc(42,i1,i2,i3) = (dt4 * (t154 + t160 + t380 - t162 + t382 - t167 + t115 * (t384 - t386) + t7 * (t175 + t177) + t2 * (t392 - t394)) / 0.12e2);
sc(43,i1,i2,i3) = (dt4 * (t126 * (t400 - t402 + t195) + t18 * (t199 + t195 + t200) + t191 + t2 * (t100 + t301 + t393) / 0.2e1 - t122 * (t409 - t410) - t14 * (t209 + t210) + t417 + t218) / 0.12e2 - t189);
sc(44,i1,i2,i3) = -(dt4 * (t160 + t380 + t225 - t162 + t382 + t227 + t115 * (t423 - t425) + t7 * (t235 + t237) + t2 * (t392 + t394)) / 0.12e2);
sc(45,i1,i2,i3) = -(dt4 * (t246 + t435) / 0.192e3);
sc(46,i1,i2,i3) = 0;
sc(47,i1,i2,i3) = (dt4 * (t439 + t440) / 0.192e3);
sc(48,i1,i2,i3) = (t31 - dt4 * (t445 + t447 + t2 * (t300 + t302)) / 0.12e2);
sc(49,i1,i2,i3) = -(dt4 * (t440 + t454) / 0.192e3);
sc(50,i1,i2,i3) = 0;
sc(51,i1,i2,i3) = (t458 * t115 * dt4 / 0.192e3);
sc(52,i1,i2,i3) = (-dt4 * (t464 - t466 + t115 * (t468 - t470)) / 0.12e2 - t477);
sc(53,i1,i2,i3) = (dt2 * (t122 - t126) / 0.12e2 - dt4 * (t2 * (t1 + t481) / 0.16e2 - t122 * (t485 - t487) + t126 * (t490 - t486) + t115 * (t458 + t493) / 0.16e2) / 0.12e2);
sc(54,i1,i2,i3) = (dt4 * (t464 - t466 + t115 * (t501 - t503)) / 0.12e2 + t477);
sc(55,i1,i2,i3) = (t115 * t493 * dt4 / 0.192e3);
sc(56,i1,i2,i3) = (-dt4 * (t515 - t517 + t115 * (t519 - t521)) / 0.12e2 - t477);
sc(57,i1,i2,i3) = (t528 - dt4 * (t530 + t115 * (t170 + t469 + t520) / 0.2e1 + t165 * (t534 - t535 + t536) + t126 * (t540 + t536 - t541) + t546 - t152 * (t547 - t548) + t553 - t122 * (t554 - t555)) / 0.12e2);
sc(58,i1,i2,i3) = (dt4 * (t2 * (t22 + t563 + t24 - t565) + t115 * (t519 + t569 + t521 - t571) + t122 * (t194 + t486 + t541) + t578 + t581 + t18 * (t122 - t582 + t583 - t584) + t165 * (t122 - t584 - t587 + t588) - t126 * (0.2e1 * t194 - t122 + 0.2e1 * t486 + t584 + 0.2e1 * t541)) / 0.12e2 - dt2 * (t599 - t600));
sc(59,i1,i2,i3) = (dt4 * (t165 * (t604 - t606 + t536) + t126 * (t540 + t536 + t541) + t530 + t115 * (t230 + t502 + t570) / 0.2e1 + t615 - t152 * (t616 - t617) + t553 - t122 * (t554 + t555)) / 0.12e2 - t528);
sc(60,i1,i2,i3) = (t477 - dt4 * (t627 + t629 + t115 * (t569 + t571)) / 0.12e2);
sc(61,i1,i2,i3) = (dt2 * (t152 - t165) / 0.12e2 - dt4 * (t7 * (t38 + t638) / 0.16e2 + t115 * (t458 + t641) / 0.16e2 - t152 * (t645 - t647) + t165 * (t650 - t646)) / 0.12e2);
sc(62,i1,i2,i3) = (dt4 * (t7 * (t48 + t658 + t50 - t660) + t115 * (t468 + t664 + t470 - t666) + t152 * (t270 + t535 + t646) + t673 + t676 + t18 * (t152 - t677 + t678 - t679) + t126 * (t152 - t682 + t683 - t679) - t165 * (0.2e1 * t270 + 0.2e1 * t535 - t152 + 0.2e1 * t646 + t679)) / 0.12e2 - dt2 * (t694 - t695));
sc(63,i1,i2,i3) = (dt4 * (t2 * (t1 + t481 + t70 + t699) / 0.16e2 - t14 * (t66 + t703 + t67 - t705) - t122 * (t485 + t709 + t487 - t711) + t7 * (t38 + t638 + t73 + t714) / 0.16e2 + t115 * (t458 + t641 + t493 + t717) / 0.16e2 - t152 * (t645 + t721 + t647 - t723) + t18 * (t62 - t726 + t63 + t727 + t704 + t728 + t729) + t126 * (t727 + t490 - t732 + t486 + t728 + t710 + t729) + t165 * (t727 + t728 + t650 - t735 + t646 + t729 + t722)) / 0.12e2 - 0.5e1 / 0.2e1 * dt2 * (t18 + t126 + t165) + 0.2e1);
sc(64,i1,i2,i3) = (dt2 * (t694 + t695) - dt4 * (t7 * (t81 + t747 + t83 - t749) + t115 * (t501 + t753 + t503 - t755) + t152 * (t341 + t604 + t722) + t165 * (0.2e1 * t341 + 0.2e1 * t604 + t152 + t679 + 0.2e1 * t722) + t673 + t676 + t18 * (t152 - t677 + t678 + t679) + t126 * (t152 - t682 + t683 + t679)) / 0.12e2);
sc(65,i1,i2,i3) = (dt4 * (t165 * (t735 + t722) - t7 * (t73 + t714) / 0.16e2 - t115 * (t493 + t717) / 0.16e2 + t152 * (t721 + t723)) / 0.12e2 - dt2 * (t152 + t165) / 0.12e2);
sc(66,i1,i2,i3) = (dt4 * (t515 - t517 + t115 * (t788 - t790)) / 0.12e2 + t477);
sc(67,i1,i2,i3) = (dt4 * (t165 * (t534 + t535 + t536) + t126 * (t536 - t800 + t801) + t530 + t115 * (t385 + t665 + t789) / 0.2e1 + t546 - t152 * (t547 + t548) + t810 - t122 * (t811 - t812)) / 0.12e2 - t528);
sc(68,i1,i2,i3) = (dt2 * (t599 + t600) - dt4 * (t2 * (t99 + t822 + t101 - t824) + t115 * (t788 + t828 + t790 - t830) + t122 * (t400 + t710 + t801) + t126 * (0.2e1 * t400 + t122 + t584 + 0.2e1 * t710 + 0.2e1 * t801) + t578 + t581 + t18 * (t122 - t582 + t583 + t584) + t165 * (t122 + t584 - t587 + t588)) / 0.12e2);
sc(69,i1,i2,i3) = (t528 + dt4 * (t165 * (t606 + t604 - t536) + t126 * (t800 - t536 + t801) - t530 - t115 * (t424 + t754 + t829) / 0.2e1 - t615 + t152 * (t616 + t617) - t810 + t122 * (t811 + t812)) / 0.12e2);
sc(70,i1,i2,i3) = (dt4 * (t627 + t629 + t115 * (t828 + t830)) / 0.12e2 - t477);
sc(71,i1,i2,i3) = (t641 * t115 * dt4 / 0.192e3);
sc(72,i1,i2,i3) = (t477 - dt4 * (t872 + t874 + t115 * (t664 + t666)) / 0.12e2);
sc(73,i1,i2,i3) = (dt4 * (t126 * (t732 + t710) - t2 * (t70 + t699) / 0.16e2 + t122 * (t709 + t711) - t115 * (t641 + t717) / 0.16e2) / 0.12e2 - dt2 * (t122 + t126) / 0.12e2);
sc(74,i1,i2,i3) = (dt4 * (t872 + t874 + t115 * (t753 + t755)) / 0.12e2 - t477);
sc(75,i1,i2,i3) = (t115 * t717 * dt4 / 0.192e3);
sc(76,i1,i2,i3) = 0;
sc(77,i1,i2,i3) = -(dt4 * (t116 + t905) / 0.192e3);
sc(78,i1,i2,i3) = (dt4 * (t125 - t128 + t2 * (t910 - t912)) / 0.12e2 + t31);
sc(79,i1,i2,i3) = (dt4 * (t905 + t140) / 0.192e3);
sc(80,i1,i2,i3) = 0;
sc(81,i1,i2,i3) = -(dt4 * (t145 + t923) / 0.192e3);
sc(82,i1,i2,i3) = (dt4 * (t154 + t157 + t929 + t931 - t164 - t167 + t115 * (t169 + t171) + t7 * (t935 - t937) + t2 * (t941 - t943)) / 0.12e2);
sc(83,i1,i2,i3) = (dt4 * (t126 * (t193 + t194 + t195) + t18 * (t195 - t952 + t953) + t191 + t2 * (t564 + t911 + t942) / 0.2e1 - t122 * (t205 + t206) - t14 * (t960 - t961) + t215 + t966) / 0.12e2 - t189);
sc(84,i1,i2,i3) = -(dt4 * (t157 + t929 + t225 + t931 - t164 + t227 + t115 * (t229 + t231) + t7 * (t974 - t976) + t2 * (t941 + t943)) / 0.12e2);
sc(85,i1,i2,i3) = -(dt4 * (t245 + t984) / 0.192e3);
sc(86,i1,i2,i3) = (dt4 * (t252 - t254 + t7 * (t989 - t991)) / 0.12e2 + t57);
sc(87,i1,i2,i3) = (dt4 * (t165 * (t269 + t270 + t271) + t18 * (t271 - t1001 + t1002) + t267 + t7 * (t659 + t936 + t990) / 0.2e1 - t152 * (t281 + t282) + t287 + t1011 - t14 * (t1012 - t1013)) / 0.12e2 - t265);
sc(88,i1,i2,i3) = (dt2 * (t336 + t337) - dt4 * (t2 * (t910 + t1023 + t912 - t1025) + t7 * (t989 + t1029 + t991 - t1031) + t14 * (t704 + t953 + t1002) + t18 * (t14 + t319 + 0.2e1 * t704 + 0.2e1 * t953 + 0.2e1 * t1002) + t315 + t318 + t126 * (t14 + t319 - t320 + t321) + t165 * (t14 + t319 - t324 + t325)) / 0.12e2);
sc(89,i1,i2,i3) = (t265 + dt4 * (t165 * (t343 + t341 - t271) + t18 * (t1001 - t271 + t1002) - t267 - t7 * (t748 + t975 + t1030) / 0.2e1 + t152 * (t350 + t351) - t356 - t1011 + t14 * (t1012 + t1013)) / 0.12e2);
sc(90,i1,i2,i3) = (dt4 * (t364 + t366 + t7 * (t1029 + t1031)) / 0.12e2 - t57);
sc(91,i1,i2,i3) = (dt4 * (t923 + t374) / 0.192e3);
sc(92,i1,i2,i3) = -(dt4 * (t154 + t929 + t380 + t931 + t382 - t167 + t115 * (t384 + t386) + t7 * (t935 + t937) + t2 * (t1077 - t1079)) / 0.12e2);
sc(93,i1,i2,i3) = (t189 + dt4 * (t126 * (t402 + t400 - t195) + t18 * (t952 - t195 + t953) - t191 - t2 * (t823 + t1024 + t1078) / 0.2e1 + t122 * (t409 + t410) + t14 * (t960 + t961) - t417 - t966) / 0.12e2);
sc(94,i1,i2,i3) = (dt4 * (t929 + t380 + t225 + t931 + t382 + t227 + t115 * (t423 + t425) + t7 * (t974 + t976) + t2 * (t1077 + t1079)) / 0.12e2);
sc(95,i1,i2,i3) = (dt4 * (t984 + t435) / 0.192e3);
sc(96,i1,i2,i3) = 0;
sc(97,i1,i2,i3) = -(dt4 * (t439 + t1111) / 0.192e3);
sc(98,i1,i2,i3) = (dt4 * (t445 + t447 + t2 * (t1023 + t1025)) / 0.12e2 - t31);
sc(99,i1,i2,i3) = (dt4 * (t1111 + t454) / 0.192e3);
sc(100,i1,i2,i3) = 0;
sc(101,i1,i2,i3) = 0;
sc(102,i1,i2,i3) = 0;
sc(103,i1,i2,i3) = (t481 * t2 * dt4 / 0.192e3);
sc(104,i1,i2,i3) = 0;
sc(105,i1,i2,i3) = 0;
sc(106,i1,i2,i3) = 0;
sc(107,i1,i2,i3) = (dt4 * (t1127 + t1128) / 0.192e3);
sc(108,i1,i2,i3) = (t31 - dt4 * (t1133 + t1135 + t2 * (t563 + t565)) / 0.12e2);
sc(109,i1,i2,i3) = -(dt4 * (t1128 + t1142) / 0.192e3);
sc(110,i1,i2,i3) = 0;
sc(111,i1,i2,i3) = (t638 * t7 * dt4 / 0.192e3);
sc(112,i1,i2,i3) = (t57 - dt4 * (t1150 + t1152 + t7 * (t658 + t660)) / 0.12e2);
sc(113,i1,i2,i3) = (dt4 * (t18 * (t726 + t704) + t14 * (t703 + t705) - t2 * (t481 + t699) / 0.16e2 - t7 * (t638 + t714) / 0.16e2) / 0.12e2 - dt2 * (t14 + t18) / 0.12e2);
sc(114,i1,i2,i3) = (dt4 * (t1150 + t1152 + t7 * (t747 + t749)) / 0.12e2 - t57);
sc(115,i1,i2,i3) = (t7 * t714 * dt4 / 0.192e3);
sc(116,i1,i2,i3) = 0;
sc(117,i1,i2,i3) = -(dt4 * (t1127 + t1182) / 0.192e3);
sc(118,i1,i2,i3) = (dt4 * (t1133 + t1135 + t2 * (t822 + t824)) / 0.12e2 - t31);
sc(119,i1,i2,i3) = (dt4 * (t1182 + t1142) / 0.192e3);
sc(120,i1,i2,i3) = 0;
sc(121,i1,i2,i3) = 0;
sc(122,i1,i2,i3) = 0;
sc(123,i1,i2,i3) = (t2 * t699 * dt4 / 0.192e3);
sc(124,i1,i2,i3) = 0;
sc(125,i1,i2,i3) = 0;
             ! #Include "../include/defineStencilVariables2dOrder2Curvilinear.h"
                    end do
                    end do
                    end do
         ! ----- END LOOPS ----
          end if
      if( gridIsImplicit.eq.0 )then 
     ! ------- EXPLICIT update the solution ---------
              if( orderInTime.eq.2 .and. orderOfAccuracy.gt.2 )then
                  stop 2222
         ! if( addForcing.eq.0 )then
         !   updateWaveOpt(3,4,2,curvilinear,NOFORCING)
         ! else 
         !   updateWaveOpt(3,4,2,curvilinear,FORCING)
         ! end if
              else
                  if( addForcing.eq.0 )then
                          if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
                              write(*,'("advWaveStencil: ADVANCE dim=3 order=4 orderInTime=4, grid=curvilinear... t=",e10.2)') t
                          end if
             ! --- TAYLOR TIME-STEPPING --- 
                          m=0 ! component number 
                          ec = 0 ! component number
                          if( forcingOption.eq.helmholtzForcing )then
                              coswt = cos(omega*t)
                          end if 
                          fv(m)=0.
             ! >>>> NOTE: NO-MASK IS FASTER for square128
             ! beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
! Stencil: nd=3, orderOfAccuracy=4, gridType=Curvilinear
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ sc(  3,i1,i2,i3)*u(i1+0,i2-2,i3-2,m)                                                                              + sc(  7,i1,i2,i3)*u(i1-1,i2-1,i3-2,m) + sc(  8,i1,i2,i3)*u(i1+0,i2-1,i3-2,m) + sc(  9,i1,i2,i3)*u(i1+1,i2-1,i3-2,m)                                       + sc( 11,i1,i2,i3)*u(i1-2,i2+0,i3-2,m) + sc( 12,i1,i2,i3)*u(i1-1,i2+0,i3-2,m) + sc( 13,i1,i2,i3)*u(i1+0,i2+0,i3-2,m) + sc( 14,i1,i2,i3)*u(i1+1,i2+0,i3-2,m) + sc( 15,i1,i2,i3)*u(i1+2,i2+0,i3-2,m)+ sc( 17,i1,i2,i3)*u(i1-1,i2+1,i3-2,m) + sc( 18,i1,i2,i3)*u(i1+0,i2+1,i3-2,m) + sc( 19,i1,i2,i3)*u(i1+1,i2+1,i3-2,m)                                       + sc( 23,i1,i2,i3)*u(i1+0,i2+2,i3-2,m)                                                                              + sc( 27,i1,i2,i3)*u(i1-1,i2-2,i3-1,m) + sc( 28,i1,i2,i3)*u(i1+0,i2-2,i3-1,m) + sc( 29,i1,i2,i3)*u(i1+1,i2-2,i3-1,m)                                       + sc( 31,i1,i2,i3)*u(i1-2,i2-1,i3-1,m) + sc( 32,i1,i2,i3)*u(i1-1,i2-1,i3-1,m) + sc( 33,i1,i2,i3)*u(i1+0,i2-1,i3-1,m) + sc( 34,i1,i2,i3)*u(i1+1,i2-1,i3-1,m) + sc( 35,i1,i2,i3)*u(i1+2,i2-1,i3-1,m)+ sc( 36,i1,i2,i3)*u(i1-2,i2+0,i3-1,m) + sc( 37,i1,i2,i3)*u(i1-1,i2+0,i3-1,m) + sc( 38,i1,i2,i3)*u(i1+0,i2+0,i3-1,m) + sc( 39,i1,i2,i3)*u(i1+1,i2+0,i3-1,m) + sc( 40,i1,i2,i3)*u(i1+2,i2+0,i3-1,m)+ sc( 41,i1,i2,i3)*u(i1-2,i2+1,i3-1,m) + sc( 42,i1,i2,i3)*u(i1-1,i2+1,i3-1,m) + sc( 43,i1,i2,i3)*u(i1+0,i2+1,i3-1,m) + sc( 44,i1,i2,i3)*u(i1+1,i2+1,i3-1,m) + sc( 45,i1,i2,i3)*u(i1+2,i2+1,i3-1,m)+ sc( 47,i1,i2,i3)*u(i1-1,i2+2,i3-1,m) + sc( 48,i1,i2,i3)*u(i1+0,i2+2,i3-1,m) + sc( 49,i1,i2,i3)*u(i1+1,i2+2,i3-1,m)                                       + sc( 51,i1,i2,i3)*u(i1-2,i2-2,i3+0,m) + sc( 52,i1,i2,i3)*u(i1-1,i2-2,i3+0,m) + sc( 53,i1,i2,i3)*u(i1+0,i2-2,i3+0,m) + sc( 54,i1,i2,i3)*u(i1+1,i2-2,i3+0,m) + sc( 55,i1,i2,i3)*u(i1+2,i2-2,i3+0,m)+ sc( 56,i1,i2,i3)*u(i1-2,i2-1,i3+0,m) + sc( 57,i1,i2,i3)*u(i1-1,i2-1,i3+0,m) + sc( 58,i1,i2,i3)*u(i1+0,i2-1,i3+0,m) + sc( 59,i1,i2,i3)*u(i1+1,i2-1,i3+0,m) + sc( 60,i1,i2,i3)*u(i1+2,i2-1,i3+0,m)+ sc( 61,i1,i2,i3)*u(i1-2,i2+0,i3+0,m) + sc( 62,i1,i2,i3)*u(i1-1,i2+0,i3+0,m) + sc( 63,i1,i2,i3)*u(i1+0,i2+0,i3+0,m) + sc( 64,i1,i2,i3)*u(i1+1,i2+0,i3+0,m) + sc( 65,i1,i2,i3)*u(i1+2,i2+0,i3+0,m)+ sc( 66,i1,i2,i3)*u(i1-2,i2+1,i3+0,m) + sc( 67,i1,i2,i3)*u(i1-1,i2+1,i3+0,m) + sc( 68,i1,i2,i3)*u(i1+0,i2+1,i3+0,m) + sc( 69,i1,i2,i3)*u(i1+1,i2+1,i3+0,m) + sc( 70,i1,i2,i3)*u(i1+2,i2+1,i3+0,m)+ sc( 71,i1,i2,i3)*u(i1-2,i2+2,i3+0,m) + sc( 72,i1,i2,i3)*u(i1-1,i2+2,i3+0,m) + sc( 73,i1,i2,i3)*u(i1+0,i2+2,i3+0,m) + sc( 74,i1,i2,i3)*u(i1+1,i2+2,i3+0,m) + sc( 75,i1,i2,i3)*u(i1+2,i2+2,i3+0,m)+ sc( 77,i1,i2,i3)*u(i1-1,i2-2,i3+1,m) + sc( 78,i1,i2,i3)*u(i1+0,i2-2,i3+1,m) + sc( 79,i1,i2,i3)*u(i1+1,i2-2,i3+1,m)                                       + sc( 81,i1,i2,i3)*u(i1-2,i2-1,i3+1,m) + sc( 82,i1,i2,i3)*u(i1-1,i2-1,i3+1,m) + sc( 83,i1,i2,i3)*u(i1+0,i2-1,i3+1,m) + sc( 84,i1,i2,i3)*u(i1+1,i2-1,i3+1,m) + sc( 85,i1,i2,i3)*u(i1+2,i2-1,i3+1,m)+ sc( 86,i1,i2,i3)*u(i1-2,i2+0,i3+1,m) + sc( 87,i1,i2,i3)*u(i1-1,i2+0,i3+1,m) + sc( 88,i1,i2,i3)*u(i1+0,i2+0,i3+1,m) + sc( 89,i1,i2,i3)*u(i1+1,i2+0,i3+1,m) + sc( 90,i1,i2,i3)*u(i1+2,i2+0,i3+1,m)+ sc( 91,i1,i2,i3)*u(i1-2,i2+1,i3+1,m) + sc( 92,i1,i2,i3)*u(i1-1,i2+1,i3+1,m) + sc( 93,i1,i2,i3)*u(i1+0,i2+1,i3+1,m) + sc( 94,i1,i2,i3)*u(i1+1,i2+1,i3+1,m) + sc( 95,i1,i2,i3)*u(i1+2,i2+1,i3+1,m)+ sc( 97,i1,i2,i3)*u(i1-1,i2+2,i3+1,m) + sc( 98,i1,i2,i3)*u(i1+0,i2+2,i3+1,m) + sc( 99,i1,i2,i3)*u(i1+1,i2+2,i3+1,m)                                       + sc(103,i1,i2,i3)*u(i1+0,i2-2,i3+2,m)                                                                              + sc(107,i1,i2,i3)*u(i1-1,i2-1,i3+2,m) + sc(108,i1,i2,i3)*u(i1+0,i2-1,i3+2,m) + sc(109,i1,i2,i3)*u(i1+1,i2-1,i3+2,m)                                       + sc(111,i1,i2,i3)*u(i1-2,i2+0,i3+2,m) + sc(112,i1,i2,i3)*u(i1-1,i2+0,i3+2,m) + sc(113,i1,i2,i3)*u(i1+0,i2+0,i3+2,m) + sc(114,i1,i2,i3)*u(i1+1,i2+0,i3+2,m) + sc(115,i1,i2,i3)*u(i1+2,i2+0,i3+2,m)+ sc(117,i1,i2,i3)*u(i1-1,i2+1,i3+2,m) + sc(118,i1,i2,i3)*u(i1+0,i2+1,i3+2,m) + sc(119,i1,i2,i3)*u(i1+1,i2+1,i3+2,m)                                       + sc(123,i1,i2,i3)*u(i1+0,i2+2,i3+2,m)                                                                              
                              end do
                              end do
                              end do
             ! endLoopsMask()
                  else
                          if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
                              write(*,'("advWaveStencil: ADVANCE dim=3 order=4 orderInTime=4, grid=curvilinear... t=",e10.2)') t
                          end if
             ! --- TAYLOR TIME-STEPPING --- 
                          m=0 ! component number 
                          ec = 0 ! component number
                          if( forcingOption.eq.helmholtzForcing )then
                              coswt = cos(omega*t)
                          end if 
                          fv(m)=0.
             ! >>>> NOTE: NO-MASK IS FASTER for square128
             ! beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                      if( forcingOption.eq.twilightZoneForcing )then
                                                  call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                                  call ogDeriv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evtt(m) )
                                                  call ogDeriv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxx(m) )
                                                  call ogDeriv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evyy(m) )
                                                  call ogDeriv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evzz(m) )
                                              fv(m) = evtt(m) - csq*( evxx(m) + evyy(m)  + evzz(m) )
                        ! Correct forcing for fourth-order ME in 3D
                                                    call ogDeriv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evtttt(m) )
                                                    call ogDeriv(ep, 0,4,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxxxx(m) )
                                                    call ogDeriv(ep, 0,2,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxxyy(m) )
                                                    call ogDeriv(ep, 0,2,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxxzz(m) )
                                                    call ogDeriv(ep, 0,0,2,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evyyzz(m) )
                                                    call ogDeriv(ep, 0,0,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evyyyy(m) )
                                                    call ogDeriv(ep, 0,0,0,4, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evzzzz(m) )
                                                fv(m) = fv(m) + (dtSq/12.)*evtttt(m) - (cdtsq12/dtSq)*( evxxxx(m) + 2.*( evxxyy(m) + evxxzz(m) + evyyzz(m) ) + evyyyy(m) + evzzzz(m) )       
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
                             !fv(m) = fv(m) -( f(i1,i2,i3,freq) + cdtSqBy12*( cSq*(fxx23(i1,i2,i3,freq) + fyy23(i1,i2,i3,freq) + fzz23(i1,i2,i3,freq)) - omega*omega*f(i1,i2,i3,freq)) )*coswt 
                                          end do ! do freq  
                                    else if( addForcing.ne.0 )then  
                                          fv(m) = f(i1,i2,i3,0)
                                    end if
! Stencil: nd=3, orderOfAccuracy=4, gridType=Curvilinear
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ sc(  3,i1,i2,i3)*u(i1+0,i2-2,i3-2,m)                                                                              + sc(  7,i1,i2,i3)*u(i1-1,i2-1,i3-2,m) + sc(  8,i1,i2,i3)*u(i1+0,i2-1,i3-2,m) + sc(  9,i1,i2,i3)*u(i1+1,i2-1,i3-2,m)                                       + sc( 11,i1,i2,i3)*u(i1-2,i2+0,i3-2,m) + sc( 12,i1,i2,i3)*u(i1-1,i2+0,i3-2,m) + sc( 13,i1,i2,i3)*u(i1+0,i2+0,i3-2,m) + sc( 14,i1,i2,i3)*u(i1+1,i2+0,i3-2,m) + sc( 15,i1,i2,i3)*u(i1+2,i2+0,i3-2,m)+ sc( 17,i1,i2,i3)*u(i1-1,i2+1,i3-2,m) + sc( 18,i1,i2,i3)*u(i1+0,i2+1,i3-2,m) + sc( 19,i1,i2,i3)*u(i1+1,i2+1,i3-2,m)                                       + sc( 23,i1,i2,i3)*u(i1+0,i2+2,i3-2,m)                                                                              + sc( 27,i1,i2,i3)*u(i1-1,i2-2,i3-1,m) + sc( 28,i1,i2,i3)*u(i1+0,i2-2,i3-1,m) + sc( 29,i1,i2,i3)*u(i1+1,i2-2,i3-1,m)                                       + sc( 31,i1,i2,i3)*u(i1-2,i2-1,i3-1,m) + sc( 32,i1,i2,i3)*u(i1-1,i2-1,i3-1,m) + sc( 33,i1,i2,i3)*u(i1+0,i2-1,i3-1,m) + sc( 34,i1,i2,i3)*u(i1+1,i2-1,i3-1,m) + sc( 35,i1,i2,i3)*u(i1+2,i2-1,i3-1,m)+ sc( 36,i1,i2,i3)*u(i1-2,i2+0,i3-1,m) + sc( 37,i1,i2,i3)*u(i1-1,i2+0,i3-1,m) + sc( 38,i1,i2,i3)*u(i1+0,i2+0,i3-1,m) + sc( 39,i1,i2,i3)*u(i1+1,i2+0,i3-1,m) + sc( 40,i1,i2,i3)*u(i1+2,i2+0,i3-1,m)+ sc( 41,i1,i2,i3)*u(i1-2,i2+1,i3-1,m) + sc( 42,i1,i2,i3)*u(i1-1,i2+1,i3-1,m) + sc( 43,i1,i2,i3)*u(i1+0,i2+1,i3-1,m) + sc( 44,i1,i2,i3)*u(i1+1,i2+1,i3-1,m) + sc( 45,i1,i2,i3)*u(i1+2,i2+1,i3-1,m)+ sc( 47,i1,i2,i3)*u(i1-1,i2+2,i3-1,m) + sc( 48,i1,i2,i3)*u(i1+0,i2+2,i3-1,m) + sc( 49,i1,i2,i3)*u(i1+1,i2+2,i3-1,m)                                       + sc( 51,i1,i2,i3)*u(i1-2,i2-2,i3+0,m) + sc( 52,i1,i2,i3)*u(i1-1,i2-2,i3+0,m) + sc( 53,i1,i2,i3)*u(i1+0,i2-2,i3+0,m) + sc( 54,i1,i2,i3)*u(i1+1,i2-2,i3+0,m) + sc( 55,i1,i2,i3)*u(i1+2,i2-2,i3+0,m)+ sc( 56,i1,i2,i3)*u(i1-2,i2-1,i3+0,m) + sc( 57,i1,i2,i3)*u(i1-1,i2-1,i3+0,m) + sc( 58,i1,i2,i3)*u(i1+0,i2-1,i3+0,m) + sc( 59,i1,i2,i3)*u(i1+1,i2-1,i3+0,m) + sc( 60,i1,i2,i3)*u(i1+2,i2-1,i3+0,m)+ sc( 61,i1,i2,i3)*u(i1-2,i2+0,i3+0,m) + sc( 62,i1,i2,i3)*u(i1-1,i2+0,i3+0,m) + sc( 63,i1,i2,i3)*u(i1+0,i2+0,i3+0,m) + sc( 64,i1,i2,i3)*u(i1+1,i2+0,i3+0,m) + sc( 65,i1,i2,i3)*u(i1+2,i2+0,i3+0,m)+ sc( 66,i1,i2,i3)*u(i1-2,i2+1,i3+0,m) + sc( 67,i1,i2,i3)*u(i1-1,i2+1,i3+0,m) + sc( 68,i1,i2,i3)*u(i1+0,i2+1,i3+0,m) + sc( 69,i1,i2,i3)*u(i1+1,i2+1,i3+0,m) + sc( 70,i1,i2,i3)*u(i1+2,i2+1,i3+0,m)+ sc( 71,i1,i2,i3)*u(i1-2,i2+2,i3+0,m) + sc( 72,i1,i2,i3)*u(i1-1,i2+2,i3+0,m) + sc( 73,i1,i2,i3)*u(i1+0,i2+2,i3+0,m) + sc( 74,i1,i2,i3)*u(i1+1,i2+2,i3+0,m) + sc( 75,i1,i2,i3)*u(i1+2,i2+2,i3+0,m)+ sc( 77,i1,i2,i3)*u(i1-1,i2-2,i3+1,m) + sc( 78,i1,i2,i3)*u(i1+0,i2-2,i3+1,m) + sc( 79,i1,i2,i3)*u(i1+1,i2-2,i3+1,m)                                       + sc( 81,i1,i2,i3)*u(i1-2,i2-1,i3+1,m) + sc( 82,i1,i2,i3)*u(i1-1,i2-1,i3+1,m) + sc( 83,i1,i2,i3)*u(i1+0,i2-1,i3+1,m) + sc( 84,i1,i2,i3)*u(i1+1,i2-1,i3+1,m) + sc( 85,i1,i2,i3)*u(i1+2,i2-1,i3+1,m)+ sc( 86,i1,i2,i3)*u(i1-2,i2+0,i3+1,m) + sc( 87,i1,i2,i3)*u(i1-1,i2+0,i3+1,m) + sc( 88,i1,i2,i3)*u(i1+0,i2+0,i3+1,m) + sc( 89,i1,i2,i3)*u(i1+1,i2+0,i3+1,m) + sc( 90,i1,i2,i3)*u(i1+2,i2+0,i3+1,m)+ sc( 91,i1,i2,i3)*u(i1-2,i2+1,i3+1,m) + sc( 92,i1,i2,i3)*u(i1-1,i2+1,i3+1,m) + sc( 93,i1,i2,i3)*u(i1+0,i2+1,i3+1,m) + sc( 94,i1,i2,i3)*u(i1+1,i2+1,i3+1,m) + sc( 95,i1,i2,i3)*u(i1+2,i2+1,i3+1,m)+ sc( 97,i1,i2,i3)*u(i1-1,i2+2,i3+1,m) + sc( 98,i1,i2,i3)*u(i1+0,i2+2,i3+1,m) + sc( 99,i1,i2,i3)*u(i1+1,i2+2,i3+1,m)                                       + sc(103,i1,i2,i3)*u(i1+0,i2-2,i3+2,m)                                                                              + sc(107,i1,i2,i3)*u(i1-1,i2-1,i3+2,m) + sc(108,i1,i2,i3)*u(i1+0,i2-1,i3+2,m) + sc(109,i1,i2,i3)*u(i1+1,i2-1,i3+2,m)                                       + sc(111,i1,i2,i3)*u(i1-2,i2+0,i3+2,m) + sc(112,i1,i2,i3)*u(i1-1,i2+0,i3+2,m) + sc(113,i1,i2,i3)*u(i1+0,i2+0,i3+2,m) + sc(114,i1,i2,i3)*u(i1+1,i2+0,i3+2,m) + sc(115,i1,i2,i3)*u(i1+2,i2+0,i3+2,m)+ sc(117,i1,i2,i3)*u(i1-1,i2+1,i3+2,m) + sc(118,i1,i2,i3)*u(i1+0,i2+1,i3+2,m) + sc(119,i1,i2,i3)*u(i1+1,i2+1,i3+2,m)                                       + sc(123,i1,i2,i3)*u(i1+0,i2+2,i3+2,m)                                                                              +dtSq*fv(m)
                              end do
                              end do
                              end do
             ! endLoopsMask()
                  end if
              end if 
      else
     ! --- IMPLICIT: Fill in RHS to implicit time-stepping -----
          stop 1111
      end if
      return
      end
