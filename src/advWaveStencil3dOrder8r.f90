! This file automatically generated from advWaveStencil.bf90 with bpp.
  subroutine advWaveStencil3dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,sc,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
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
   ! real stencilCoeff(0:*)   ! holds stencil coeffs
   real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1) 
   real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
   real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
   real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds coeff of Laplacian for HA scheme
   real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
   real etax(nd1a:nd1b)  ! superGrid functions
   real etay(nd2a:nd2b)
   real etaz(nd3a:nd3b)
   integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
   integer bc(0:1,0:2),ierr  
   real frequencyArray(0:*)
   double precision pdb  ! pointer to the data base
   integer ipar(0:*)
   real rpar(0:*)
  ! integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b
  ! real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  ! real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  ! real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  ! real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  ! ! real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
  ! real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1) 
  ! real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  ! real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
  ! real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds coeff of Laplacian for HA scheme
  ! real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
  ! integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
  ! integer bc(0:1,0:2),ierr
  ! real frequencyArray(0:*)
  ! integer ipar(0:*)
  ! real rpar(0:*)
     real sc(1:729,nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
     real scr(1:729)    
  !     ---- local variables -----
  integer gridIndexRange(0:1,0:2)
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
  real omega, coswt, damp
  integer ok,getInt,getReal
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
! Define variables to valuate stencil coefficients, dim=3, order=8, gridType=Rectangular
! File generated by cgWave/maple/writeStencilFiles.mpl
integer i1m3,i1m2,i1m1,i1p1,i1p2,i1p3
integer i2m3,i2m2,i2m1,i2p1,i2p2,i2p3
integer i3m3,i3m2,i3m1,i3p1,i3p2,i3p3
real t0,t1,t3,t4,t9,t12,t14,t15,t16,t19,t22,t23,t24,t25,t30,t31,t32,t33,t34,t36,t37,t40,t43,t44,t45,t46,t47,t48,t50,t51,t52,t53,t55,t56,t57,t59,t60,t61,t62,t67,t69,t71,t72,t73,t74,t75,t76,t77,t83,t85,t86,t87,t88,t90,t91,t92,t94,t97,t99,t102,t106,t107,t108,t109,t110,t111,t112,t113,t114,t115,t116,t117,t119,t120,t121,t123,t124,t128,t129,t131,t132,t133,t134,t135,t136,t137,t138,t142,t147,t148,t150,t151,t152,t154,t157,t159,t161,t162,t163,t164,t165,t166,t168,t169,t170,t175,t176,t178,t179,t180,t181,t182,t183,t190,t191,t192,t194,t195,t198,t201,t202,t204,t205,t207,t208,t210,t211,t212,t213,t214,t215,t216,t217,t218,t219,t220,t221,t222,t223,t226,t229,t235,t236,t237,t238,t240,t241,t243,t244,t245,t246,t247,t248,t250,t251,t253,t254,t260,t261,t262,t263,t268,t278,t279,t280,t281,t282,t287,t297,t298,t301,t302,t305,t308,t309,t310,t312,t313,t315,t320,t321,t322,t323,t325,t326,t327,t329,t333,t340,t341,t342,t344,t346,t347,t348,t349,t350,t351,t356,t358,t359,t360,t361,t363,t364,t365,t368,t369,t374,t375,t376,t377,t378,t381,t382,t383,t384,t385,t386,t387,t388,t389,t390,t393,t401,t402,t403,t404,t405,t406,t411,t414,t415,t418,t420,t421,t422,t424,t429,t430,t431,t432,t434,t435,t436,t438,t442,t447,t448,t449,t453,t454,t455,t460,t461,t462,t463,t466,t467,t468,t469,t470,t471,t472,t473,t482,t483,t485,t486,t489,t492,t493,t494,t495,t496,t497,t499,t500,t501,t504,t505,t506,t507,t508,t511,t514,t517,t518,t523,t524,t525,t526,t528,t530,t532,t533,t534,t545,t546,t550,t551,t552,t553,t555,t563,t565,t567,t568,t575,t577,t580,t581,t586,t587,t591,t592,t593,t594,t596,t604,t612,t613,t615,t616,t623,t625,t627,t629,t636,t637,t638,t643,t644,t645,t649,t652,t653,t654,t656,t657,t662,t664,t665,t666,t672,t674,t675,t677,t679,t683,t685,t686,t687,t689,t690,t695,t696,t698,t705,t706,t707,t711,t714,t715,t716,t719,t722,t728,t729,t730,t732,t733,t735,t736,t742,t743,t750,t759,t764,t765,t767,t772,t773,t774,t776,t783,t784,t786,t791,t792,t793,t794,t797,t806,t807,t811,t814,t815,t816,t817,t818,t820,t821,t823,t825,t830,t831,t832,t843,t844,t848,t850,t858,t860,t862,t869,t871,t873,t875,t882,t886,t889,t890,t892,t893,t898,t900,t901,t907,t909,t910,t911,t913,t915,t920,t924,t927,t928,t929,t932,t935,t945,t948,t949,t950,t951,t953,t956,t961,t962,t963,t975,t976,t979,t993,t994,t995,t996,t997,t998,t999,t1003,t1007,t1012,t1013,t1014,t1015,t1016,t1017,t1030,t1031,t1032,t1051,t1080,t1091,t1098,t1099,t1100,t1127,t1133,t1137,t1174,t1190,t1195,t1206,t1213,t1214,t1215,t1216,t1256,t1264,t1267,t1275,t1279;
    ! #Include "../include/defineStencilVariables2dOrder2Rectangular.h"
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
   ! new way: extract parameters from the dataBase
    ok=getReal(pdb,'damp',damp) 
    if( ok.eq.0 )then
      write(*,'("*** advWaveStencil:getReal:ERROR: unable to find damp")') 
      stop 1133
    end if
   if( damp.ne.0. )then
     write(*,*) "advWaveStencil: FINISH ME FOR damping"
     stop 1010
   end if 
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
     ! ---- CARTESIAN GRID -----
     i1=0; i2=0; i3=0;
! Evaluate stencil coefficients, dim=3, order=8, gridType=Rectangular
! File generated by cgWave/maple/writeStencilFiles.mpl
i1m1=i1-1; i1p1=i1+1;
i2m1=i2-1; i2p1=i2+1;
i3m1=i3-1; i3p1=i3+1;
i1m2=i1-2; i1p2=i1+2;
i2m2=i2-2; i2p2=i2+2;
i3m2=i3-2; i3p2=i3+2;
i1m3=i1-3; i1p3=i1+3;
i2m3=i2-3; i2p3=i2+3;
i3m3=i3-3; i3p3=i3+3;
t1 = dt2 * cdz2i;
t3 = cdz2i ** 2;
t4 = t3 ** 2;
t9 = t3 * cdz2i;
t12 = -t1 / 0.560e3 + dt8 * t4 / 0.20160e5 + 0.7e1 / 0.2880e4 * dt4 * t3 - dt6 * t9 / 0.1440e4;
t14 = dt4 * cdy2i * cdz2i;
t15 = t14 / 0.540e3;
t16 = dt6 * t3;
t19 = dt8 * t9;
t22 = t15 - t16 * cdy2i / 0.720e3 + t19 * cdy2i / 0.5040e4;
t23 = dt4 * cdx2i;
t24 = t23 * cdz2i;
t25 = t24 / 0.540e3;
t30 = t25 - t16 * cdx2i / 0.720e3 + t19 * cdx2i / 0.5040e4;
t31 = cdy2i * cdz2i;
t32 = t31 / 0.45e2;
t33 = cdx2i * cdz2i;
t34 = t33 / 0.45e2;
t36 = cdy2i / 0.45e2;
t37 = cdx2i / 0.45e2;
t40 = -t32 - t34 - cdz2i * (0.2e1 / 0.5e1 * cdz2i + t36 + t37);
t43 = 0.2e1 / 0.3e1 * t3;
t44 = t31 / 0.3e1;
t45 = t33 / 0.3e1;
t46 = 0.2e1 * cdz2i;
t47 = cdy2i / 0.6e1;
t48 = cdx2i / 0.6e1;
t50 = cdz2i * (t46 + t47 + t48);
t51 = 0.4e1 * cdz2i;
t52 = 0.2e1 * cdy2i;
t53 = 0.2e1 * cdx2i;
t55 = cdz2i * (t51 + t52 + t53);
t56 = t55 / 0.12e2;
t57 = -t43 - t44 - t45 - t50 - t56;
t59 = t3 * cdy2i;
t60 = t59 / 0.3e1;
t61 = t3 * cdx2i;
t62 = t61 / 0.3e1;
t67 = 0.2e1 * t9 * cdy2i;
t69 = 0.2e1 * t9 * cdx2i;
t71 = 0.2e1 * t31;
t72 = 0.2e1 * t33;
t73 = -0.2e1 * t3 - t71 - t72 - t55;
t74 = cdz2i * t73;
t75 = 0.2e1 * t59;
t76 = 0.2e1 * t61;
t77 = 0.2e1 * t9;
t83 = 0.8e1 / 0.315e3 * t1;
t85 = t14 / 0.864e3;
t86 = t31 / 0.6e1;
t87 = t3 / 0.12e2;
t88 = -t86 - t87;
t90 = cdy2i ** 2;
t91 = t90 / 0.12e2;
t92 = -t86 - t91;
t94 = -cdy2i * t88 - cdz2i * t92;
t97 = dt8 * t90;
t99 = t97 * t3 / 0.3360e4;
t102 = cdx2i * cdy2i;
t106 = dt6 * cdx2i * t31;
t107 = t106 / 0.720e3;
t108 = dt8 * t3 * t102 / 0.1680e4 - t107;
t109 = 0.4e1 / 0.135e3 * t14;
t110 = t31 / 0.2e1;
t111 = t3 / 0.3e1;
t112 = t33 / 0.6e1;
t113 = -t110 - t111 - t112 - t50;
t114 = cdy2i * t113;
t115 = 0.9e1 / 0.2e1 * t31;
t116 = t102 / 0.6e1;
t117 = 0.4e1 * cdy2i;
t119 = cdy2i * (t46 + t117 + t53);
t120 = t119 / 0.12e2;
t121 = -t115 - t116 - t120;
t123 = t102 * cdz2i;
t124 = t123 / 0.3e1;
t128 = 0.8e1 * t59;
t129 = t128 + t76 - t74;
t131 = 0.6e1 * t59;
t132 = 0.6e1 * t31;
t133 = -t132 - t72 - t55;
t134 = cdy2i * t133;
t135 = 0.2e1 * t102;
t136 = -t132 - t135 - t119;
t137 = cdz2i * t136;
t138 = 0.4e1 * t123;
t142 = 0.6e1 * t61 * cdy2i;
t147 = t24 / 0.864e3;
t148 = -t112 - t87;
t150 = cdx2i ** 2;
t151 = t150 / 0.12e2;
t152 = -t112 - t151;
t154 = -cdx2i * t148 - cdz2i * t152;
t157 = dt8 * t150;
t159 = t157 * t3 / 0.3360e4;
t161 = 0.4e1 / 0.135e3 * t24;
t162 = t33 / 0.2e1;
t163 = -t162 - t111 - t86 - t50;
t164 = cdx2i * t163;
t165 = 0.9e1 / 0.2e1 * t33;
t166 = 0.4e1 * cdx2i;
t168 = cdx2i * (t46 + t52 + t166);
t169 = t168 / 0.12e2;
t170 = -t116 - t165 - t169;
t175 = 0.8e1 * t61;
t176 = t175 + t75 - t74;
t178 = 0.6e1 * t61;
t179 = 0.6e1 * t33;
t180 = -t179 - t71 - t55;
t181 = cdx2i * t180;
t182 = -t179 - t135 - t168;
t183 = cdz2i * t182;
t190 = t1 / 0.5e1;
t191 = 0.41e2 / 0.120e3 * t31;
t192 = 0.41e2 / 0.120e3 * t33;
t194 = 0.41e2 / 0.120e3 * cdy2i;
t195 = 0.41e2 / 0.120e3 * cdx2i;
t198 = -t191 - t192 - cdz2i * (0.169e3 / 0.60e2 * cdz2i + t194 + t195);
t201 = 0.2e1 * t74;
t202 = 0.6e1 * cdz2i;
t204 = cdz2i * (t202 + t117 + t166);
t205 = 0.6e1 * cdy2i;
t207 = cdy2i * (t51 + t205 + t166);
t208 = 0.6e1 * cdx2i;
t210 = cdx2i * (t51 + t117 + t208);
t211 = 0.4e1 * t31;
t212 = 0.4e1 * t33;
t213 = 0.2e1 * t55;
t214 = -t3 - t204 - t207 - t210 - t211 - t212 - t213;
t215 = cdz2i * t214;
t216 = 0.4e1 * t59;
t217 = 0.4e1 * t61;
t218 = 0.8e1 * t31;
t219 = -t218 - t212 - t213;
t220 = cdy2i * t219;
t221 = 0.8e1 * t33;
t222 = -t221 - t211 - t213;
t223 = cdx2i * t222;
t226 = 0.10e2 * t59;
t229 = 0.10e2 * t61;
t235 = 0.2e1 / 0.3e1 * t31;
t236 = 0.2e1 * t50;
t237 = t3 / 0.2e1;
t238 = -t235 - t45 - t236 - t237;
t240 = 0.2e1 / 0.3e1 * t33;
t241 = -t44 - t240 - t236 - t237;
t243 = t207 / 0.12e2;
t244 = t210 / 0.12e2;
t245 = t204 / 0.12e2;
t246 = 0.25e2 / 0.6e1 * t31;
t247 = 0.25e2 / 0.6e1 * t33;
t248 = t55 / 0.3e1;
t250 = 0.19e2 / 0.6e1 * cdy2i;
t251 = 0.19e2 / 0.6e1 * cdx2i;
t253 = cdz2i * (0.13e2 / 0.2e1 * cdz2i + t250 + t251);
t254 = -t243 - t244 - t245 - t246 - t247 - t248 - t253 - t43 - t236;
t260 = -cdz2i * t182;
t261 = -cdx2i * t180;
t262 = -t138 - t178 - t260 - t261;
t263 = cdz2i * t262;
t268 = -cdz2i * t170;
t278 = t159 + t147 - dt6 * (-cdx2i * t148 - cdz2i * t152) / 0.360e3;
t279 = -cdz2i * t136;
t280 = -cdy2i * t133;
t281 = -t138 - t131 - t279 - t280;
t282 = cdz2i * t281;
t287 = -cdz2i * t121;
t297 = t99 + t85 - dt6 * (-cdy2i * t88 - cdz2i * t92) / 0.360e3;
t298 = dt6 * t90;
t301 = t90 * cdy2i;
t302 = dt8 * t301;
t305 = t15 - t298 * cdz2i / 0.720e3 + t302 * cdz2i / 0.5040e4;
t308 = t97 * t33 / 0.1680e4 - t107;
t309 = t90 / 0.3e1;
t310 = cdz2i / 0.6e1;
t312 = cdy2i * (t310 + t52 + t48);
t313 = -t110 - t309 - t116 - t312;
t315 = -t115 - t112 - t56;
t320 = t90 * cdz2i;
t321 = 0.8e1 * t320;
t322 = t90 * cdx2i;
t323 = 0.2e1 * t322;
t325 = -0.2e1 * t90 - t71 - t135 - t119;
t326 = cdy2i * t325;
t327 = t321 + t323 - t326;
t329 = 0.6e1 * t320;
t333 = 0.6e1 * t322 * cdz2i;
t340 = t157 * t31 / 0.1680e4 - t107;
t341 = 0.16e2 * t123;
t342 = t341 - t134 - t137;
t344 = t341 - t181 - t183;
t346 = 0.6e1 * t102;
t347 = -t346 - t71 - t119;
t348 = cdx2i * t347;
t349 = -t346 - t72 - t168;
t350 = cdy2i * t349;
t351 = t341 - t348 - t350;
t356 = t106 / 0.30e2;
t358 = 0.23e2 / 0.3e1 * t123;
t359 = 0.71e2 / 0.6e1 * t31;
t360 = 0.23e2 / 0.6e1 * t33;
t361 = -t359 - t360 - t248 - t253;
t363 = t119 / 0.3e1;
t364 = 0.23e2 / 0.6e1 * t102;
t365 = 0.19e2 / 0.6e1 * cdz2i;
t368 = cdy2i * (t365 + 0.13e2 / 0.2e1 * cdy2i + t251);
t369 = -t363 - t359 - t364 - t368;
t374 = 0.19e2 / 0.54e2 * t14;
t375 = 0.8e1 * t123;
t376 = 0.2e1 * t134;
t377 = 0.2e1 * t137;
t378 = 0.3e1 * t320;
t381 = 0.3e1 * t59;
t382 = 0.8e1 * t102;
t383 = 0.2e1 * t119;
t384 = -t382 - t211 - t383;
t385 = cdx2i * t384;
t386 = 0.4e1 * t102;
t387 = -t218 - t386 - t383;
t388 = cdz2i * t387;
t389 = -t90 - t204 - t207 - t210 - t211 - t386 - t383;
t390 = cdy2i * t389;
t393 = 0.20e2 * t123;
t401 = -t341 - t261 - t260;
t402 = cdy2i * t401;
t403 = -cdx2i * t347;
t404 = -cdy2i * t349;
t405 = -t341 - t403 - t404;
t406 = cdz2i * t405;
t411 = dt6 * t150;
t414 = t150 * cdx2i;
t415 = dt8 * t414;
t418 = t25 - t411 * cdz2i / 0.720e3 + t415 * cdz2i / 0.5040e4;
t420 = cdx2i * (t310 + t47 + t53);
t421 = t150 / 0.3e1;
t422 = -t420 - t162 - t116 - t421;
t424 = -t165 - t86 - t56;
t429 = t150 * cdy2i;
t430 = 0.2e1 * t429;
t431 = t150 * cdz2i;
t432 = 0.8e1 * t431;
t434 = -0.2e1 * t150 - t72 - t135 - t168;
t435 = cdx2i * t434;
t436 = t430 + t432 - t435;
t438 = 0.6e1 * t431;
t442 = 0.6e1 * t429 * cdz2i;
t447 = 0.71e2 / 0.6e1 * t33;
t448 = 0.23e2 / 0.6e1 * t31;
t449 = -t253 - t447 - t448 - t248;
t453 = cdx2i * (t365 + t250 + 0.13e2 / 0.2e1 * cdx2i);
t454 = t168 / 0.3e1;
t455 = -t453 - t447 - t364 - t454;
t460 = 0.19e2 / 0.54e2 * t24;
t461 = 0.2e1 * t183;
t462 = 0.2e1 * t181;
t463 = 0.3e1 * t431;
t466 = 0.2e1 * t168;
t467 = -t382 - t212 - t466;
t468 = cdy2i * t467;
t469 = 0.3e1 * t61;
t470 = -t221 - t386 - t466;
t471 = cdz2i * t470;
t472 = -t150 - t204 - t207 - t210 - t212 - t386 - t466;
t473 = cdx2i * t472;
t482 = 0.35e2 / 0.9e1 * t31;
t483 = 0.35e2 / 0.9e1 * t33;
t485 = 0.35e2 / 0.9e1 * cdy2i;
t486 = 0.35e2 / 0.9e1 * cdx2i;
t489 = -t482 - t483 - cdz2i * (0.122e3 / 0.15e2 * cdz2i + t485 + t486);
t492 = 0.8e1 / 0.5e1 * t1;
t493 = 0.46e2 / 0.3e1 * t31;
t494 = 0.22e2 / 0.3e1 * t33;
t495 = t55 / 0.2e1;
t496 = 0.2e1 * t253;
t497 = -t493 - t494 - t495 - t496;
t499 = 0.22e2 / 0.3e1 * t31;
t500 = 0.46e2 / 0.3e1 * t33;
t501 = -t499 - t500 - t495 - t496;
t504 = t204 / 0.3e1;
t505 = t207 / 0.3e1;
t506 = t210 / 0.3e1;
t507 = 0.23e2 / 0.3e1 * t31;
t508 = 0.23e2 / 0.3e1 * t33;
t511 = cdz2i * (0.28e2 / 0.3e1 * cdz2i + t205 + t208);
t514 = cdy2i * (t202 + 0.28e2 / 0.3e1 * cdy2i + t208);
t517 = cdx2i * (t202 + t205 + 0.28e2 / 0.3e1 * cdx2i);
t518 = -0.7e1 / 0.12e2 * t55 - t111 - t504 - t505 - t506 - t507 - t508 - t496 - t511 - t50 - t514 - t517;
t523 = 0.2e1 * t215;
t524 = 0.2e1 * t204;
t525 = 0.2e1 * t207;
t526 = 0.2e1 * t210;
t528 = cdz2i * (t211 + t212 + t213 + t524 + t525 + t526);
t530 = cdy2i * (t383 + t211 + t386 + t524 + t525 + t526);
t532 = cdx2i * (t466 + t212 + t386 + t524 + t525 + t526);
t533 = 0.2e1 * t220;
t534 = 0.2e1 * t223;
t545 = 0.2e1 * t261;
t546 = 0.2e1 * t260;
t550 = cdy2i * (t545 + t546 + t393);
t551 = -cdx2i * t472;
t552 = -cdz2i * t470;
t553 = -cdy2i * t467;
t555 = cdz2i * (t469 + t551 + t552 + t553 + t545 + t546 + t375);
t563 = dt6 * (-cdx2i * t449 - cdz2i * t455 + t358) / 0.360e3;
t565 = t138 + t438 + t260 + t261;
t567 = -cdx2i * t434;
t568 = t432 + t430 + t567;
t575 = -cdx2i * t424 - cdz2i * t422 + t124;
t577 = dt6 * t575 / 0.360e3;
t580 = -t341 - t280 - t279;
t581 = cdx2i * t580;
t586 = 0.2e1 * t280;
t587 = 0.2e1 * t279;
t591 = cdx2i * (t586 + t587 + t393);
t592 = -cdy2i * t389;
t593 = -cdz2i * t387;
t594 = -cdx2i * t384;
t596 = cdz2i * (t381 + t586 + t587 + t375 + t592 + t593 + t594);
t604 = dt6 * (-cdy2i * t361 - cdz2i * t369 + t358) / 0.360e3;
t612 = t356 - dt8 * (-cdx2i * t580 - cdy2i * t401 - cdz2i * t405) / 0.20160e5;
t613 = t138 + t329 + t279 + t280;
t615 = -cdy2i * t325;
t616 = t321 + t323 + t615;
t623 = -cdy2i * t315 - cdz2i * t313 + t124;
t625 = dt6 * t623 / 0.360e3;
t627 = dt2 * cdy2i;
t629 = t90 ** 2;
t636 = -t627 / 0.560e3 + dt8 * t629 / 0.20160e5 + 0.7e1 / 0.2880e4 * dt4 * t90 - dt6 * t301 / 0.1440e4;
t637 = t23 * cdy2i;
t638 = t637 / 0.540e3;
t643 = t638 - t298 * cdx2i / 0.720e3 + t302 * cdx2i / 0.5040e4;
t644 = t102 / 0.45e2;
t645 = cdz2i / 0.45e2;
t649 = -t32 - t644 - cdy2i * (t645 + 0.2e1 / 0.5e1 * cdy2i + t37);
t652 = 0.2e1 / 0.3e1 * t90;
t653 = t102 / 0.3e1;
t654 = -t44 - t652 - t653 - t312 - t120;
t656 = t320 / 0.3e1;
t657 = t322 / 0.3e1;
t662 = 0.2e1 * t301 * cdz2i;
t664 = 0.2e1 * t301 * cdx2i;
t665 = 0.2e1 * t320;
t666 = 0.2e1 * t301;
t672 = 0.8e1 / 0.315e3 * t627;
t674 = t637 / 0.864e3;
t675 = -t116 - t91;
t677 = -t116 - t151;
t679 = -cdx2i * t675 - cdy2i * t677;
t683 = t157 * t90 / 0.3360e4;
t685 = 0.4e1 / 0.135e3 * t637;
t686 = t102 / 0.2e1;
t687 = -t686 - t309 - t312 - t86;
t689 = 0.9e1 / 0.2e1 * t102;
t690 = -t689 - t112 - t169;
t695 = 0.8e1 * t322;
t696 = t695 + t665 - t326;
t698 = 0.6e1 * t322;
t705 = t627 / 0.5e1;
t706 = 0.41e2 / 0.120e3 * t102;
t707 = 0.41e2 / 0.120e3 * cdz2i;
t711 = -t191 - t706 - cdy2i * (t707 + 0.169e3 / 0.60e2 * cdy2i + t195);
t714 = 0.2e1 * t326;
t715 = 0.4e1 * t320;
t716 = 0.4e1 * t322;
t719 = 0.10e2 * t320;
t722 = 0.10e2 * t322;
t728 = 0.2e1 * t312;
t729 = t90 / 0.2e1;
t730 = -t235 - t653 - t728 - t729;
t732 = 0.2e1 / 0.3e1 * t102;
t733 = -t732 - t44 - t728 - t729;
t735 = 0.25e2 / 0.6e1 * t102;
t736 = -t246 - t652 - t735 - t368 - t728 - t244 - t245 - t243 - t363;
t742 = -t138 - t698 - t404 - t403;
t743 = cdy2i * t742;
t750 = cdx2i * t687 + cdy2i * t690 - t124;
t759 = t683 + t674 - dt6 * (-cdx2i * t675 - cdy2i * t677) / 0.360e3;
t764 = t638 - t411 * cdy2i / 0.720e3 + t415 * cdy2i / 0.5040e4;
t765 = -t686 - t421 - t420 - t112;
t767 = -t689 - t86 - t120;
t772 = 0.8e1 * t429;
t773 = 0.2e1 * t431;
t774 = t772 + t773 - t435;
t776 = 0.6e1 * t429;
t783 = 0.71e2 / 0.6e1 * t102;
t784 = -t783 - t448 - t363 - t368;
t786 = -t360 - t783 - t454 - t453;
t791 = 0.19e2 / 0.54e2 * t637;
t792 = 0.3e1 * t429;
t793 = 0.2e1 * t348;
t794 = 0.2e1 * t350;
t797 = 0.3e1 * t322;
t806 = 0.35e2 / 0.9e1 * t102;
t807 = 0.35e2 / 0.9e1 * cdz2i;
t811 = -t482 - t806 - cdy2i * (t807 + 0.122e3 / 0.15e2 * cdy2i + t486);
t814 = 0.8e1 / 0.5e1 * t627;
t815 = 0.22e2 / 0.3e1 * t102;
t816 = t119 / 0.2e1;
t817 = 0.2e1 * t368;
t818 = -t493 - t815 - t816 - t817;
t820 = 0.46e2 / 0.3e1 * t102;
t821 = -t820 - t499 - t816 - t817;
t823 = 0.23e2 / 0.3e1 * t102;
t825 = -t511 - t514 - t312 - t517 - t507 - t823 - 0.7e1 / 0.12e2 * t119 - t817 - t309 - t504 - t505 - t506;
t830 = 0.2e1 * t388;
t831 = 0.2e1 * t385;
t832 = 0.2e1 * t390;
t843 = 0.2e1 * t403;
t844 = 0.2e1 * t404;
t848 = cdz2i * (t843 + t844 + t393);
t850 = cdy2i * (t552 + t553 + t797 + t551 + t843 + t844 + t375);
t858 = dt6 * (-cdx2i * t784 - cdy2i * t786 + t358) / 0.360e3;
t860 = t138 + t776 + t404 + t403;
t862 = t567 + t773 + t772;
t869 = -cdx2i * t767 - cdy2i * t765 + t124;
t871 = dt6 * t869 / 0.360e3;
t873 = dt2 * cdx2i;
t875 = t150 ** 2;
t882 = -t873 / 0.560e3 + dt8 * t875 / 0.20160e5 + 0.7e1 / 0.2880e4 * dt4 * t150 - dt6 * t414 / 0.1440e4;
t886 = -t34 - t644 - cdx2i * (t645 + t36 + 0.2e1 / 0.5e1 * cdx2i);
t889 = 0.2e1 / 0.3e1 * t150;
t890 = -t45 - t653 - t889 - t420 - t169;
t892 = t431 / 0.3e1;
t893 = t429 / 0.3e1;
t898 = 0.2e1 * t414 * cdz2i;
t900 = 0.2e1 * t414 * cdy2i;
t901 = 0.2e1 * t414;
t907 = 0.8e1 / 0.315e3 * t873;
t909 = 0.2e1 * t420;
t910 = t150 / 0.2e1;
t911 = -t240 - t653 - t909 - t910;
t913 = -t732 - t45 - t909 - t910;
t915 = -t247 - t735 - t889 - t245 - t243 - t244 - t454 - t453 - t909;
t920 = t873 / 0.5e1;
t924 = -t192 - t706 - cdx2i * (t707 + t194 + 0.169e3 / 0.60e2 * cdx2i);
t927 = 0.2e1 * t435;
t928 = 0.4e1 * t431;
t929 = 0.4e1 * t429;
t932 = 0.10e2 * t431;
t935 = 0.10e2 * t429;
t945 = -t483 - t806 - cdx2i * (t807 + t485 + 0.122e3 / 0.15e2 * cdx2i);
t948 = 0.8e1 / 0.5e1 * t873;
t949 = t168 / 0.2e1;
t950 = 0.2e1 * t453;
t951 = -t500 - t815 - t949 - t950;
t953 = -t820 - t494 - t949 - t950;
t956 = -t511 - t514 - t517 - t420 - t823 - t508 - 0.7e1 / 0.12e2 * t168 - t421 - t504 - t505 - t506 - t950;
t961 = 0.2e1 * t471;
t962 = 0.2e1 * t468;
t963 = 0.2e1 * t473;
t975 = 0.257e3 / 0.36e2 * cdy2i;
t976 = 0.257e3 / 0.36e2 * cdx2i;
t979 = 0.257e3 / 0.36e2 * cdz2i;
t993 = 0.2e1 * t511;
t994 = 0.2e1 * t514;
t995 = 0.2e1 * t517;
t996 = t204 / 0.2e1;
t997 = t207 / 0.2e1;
t998 = t210 / 0.2e1;
t999 = t507 + t508 + t3 / 0.6e1 + 0.2e1 / 0.3e1 * t55 + t496 + t993 + t994 + t995 + t996 + t997 + t998;
t1003 = t507 + t823 + t90 / 0.6e1 + 0.2e1 / 0.3e1 * t119 + t817 + t993 + t994 + t995 + t996 + t997 + t998;
t1007 = t508 + t823 + t150 / 0.6e1 + 0.2e1 / 0.3e1 * t168 + t950 + t993 + t994 + t995 + t996 + t997 + t998;
t1012 = -cdx2i * t222;
t1013 = -cdy2i * t219;
t1014 = -cdz2i * t214;
t1015 = 0.2e1 * t528;
t1016 = 0.2e1 * t530;
t1017 = 0.2e1 * t532;
t1030 = 0.2e1 * t553;
t1031 = 0.2e1 * t552;
t1032 = 0.2e1 * t551;
t1051 = 0.2e1 * t567;
t1080 = -cdx2i * t860;
t1091 = cdx2i * (t792 + t593 + t594 + t592 + t844 + t843 + t375);
t1098 = 0.2e1 * t592;
t1099 = 0.2e1 * t593;
t1100 = 0.2e1 * t594;
t1127 = t665 + t695 + t615;
t1133 = -dt6 * t750 / 0.360e3;
t1137 = 0.2e1 * t615;
t1174 = -cdy2i * t613;
t1190 = cdy2i * (t1013 + t1012 + t378 + t1014 + t587 + t586 + t375);
t1195 = -cdx2i * t565;
t1206 = cdx2i * (t1012 + t463 + t1014 + t1013 + t546 + t545 + t375);
t1213 = -cdz2i * t73;
t1214 = 0.2e1 * t1013;
t1215 = 0.2e1 * t1012;
t1216 = 0.2e1 * t1014;
t1256 = t128 + t76 + t1213;
t1264 = dt6 * (-cdy2i * t113 + t124 + t287) / 0.360e3;
t1267 = t75 + t175 + t1213;
t1275 = dt6 * (-cdx2i * t163 + t124 + t268) / 0.360e3;
t1279 = 0.2e1 * t1213;
scr(1) = 0;
scr(2) = 0;
scr(3) = 0;
scr(4) = 0;
scr(5) = 0;
scr(6) = 0;
scr(7) = 0;
scr(8) = 0;
scr(9) = 0;
scr(10) = 0;
scr(11) = 0;
scr(12) = 0;
scr(13) = 0;
scr(14) = 0;
scr(15) = 0;
scr(16) = 0;
scr(17) = 0;
scr(18) = 0;
scr(19) = 0;
scr(20) = 0;
scr(21) = 0;
scr(22) = 0;
scr(23) = 0;
scr(24) = 0;
scr(25) = 0;
scr(26) = 0;
scr(27) = 0;
scr(28) = 0;
scr(29) = 0;
scr(30) = 0;
scr(31) = 0;
scr(32) = 0;
scr(33) = 0;
scr(34) = 0;
scr(35) = 0;
scr(36) = 0;
scr(37) = 0;
scr(38) = 0;
scr(39) = 0;
scr(40) = 0;
scr(41) = t12;
scr(42) = 0;
scr(43) = 0;
scr(44) = 0;
scr(45) = 0;
scr(46) = 0;
scr(47) = 0;
scr(48) = 0;
scr(49) = 0;
scr(50) = 0;
scr(51) = 0;
scr(52) = 0;
scr(53) = 0;
scr(54) = 0;
scr(55) = 0;
scr(56) = 0;
scr(57) = 0;
scr(58) = 0;
scr(59) = 0;
scr(60) = 0;
scr(61) = 0;
scr(62) = 0;
scr(63) = 0;
scr(64) = 0;
scr(65) = 0;
scr(66) = 0;
scr(67) = 0;
scr(68) = 0;
scr(69) = 0;
scr(70) = 0;
scr(71) = 0;
scr(72) = 0;
scr(73) = 0;
scr(74) = 0;
scr(75) = 0;
scr(76) = 0;
scr(77) = 0;
scr(78) = 0;
scr(79) = 0;
scr(80) = 0;
scr(81) = 0;
scr(82) = 0;
scr(83) = 0;
scr(84) = 0;
scr(85) = 0;
scr(86) = 0;
scr(87) = 0;
scr(88) = 0;
scr(89) = 0;
scr(90) = 0;
scr(91) = 0;
scr(92) = 0;
scr(93) = 0;
scr(94) = 0;
scr(95) = 0;
scr(96) = 0;
scr(97) = 0;
scr(98) = 0;
scr(99) = 0;
scr(100) = 0;
scr(101) = 0;
scr(102) = 0;
scr(103) = 0;
scr(104) = 0;
scr(105) = 0;
scr(106) = 0;
scr(107) = 0;
scr(108) = 0;
scr(109) = 0;
scr(110) = 0;
scr(111) = 0;
scr(112) = 0;
scr(113) = t22;
scr(114) = 0;
scr(115) = 0;
scr(116) = 0;
scr(117) = 0;
scr(118) = 0;
scr(119) = 0;
scr(120) = 0;
scr(121) = t30;
scr(122) = (dt4 * t40 / 0.12e2 - dt6 * (cdz2i * t57 - t60 - t62) / 0.360e3 - dt8 * (t67 + t69 - cdz2i * (t74 - t75 - t76 - t77)) / 0.20160e5 + t83);
scr(123) = t30;
scr(124) = 0;
scr(125) = 0;
scr(126) = 0;
scr(127) = 0;
scr(128) = 0;
scr(129) = 0;
scr(130) = 0;
scr(131) = t22;
scr(132) = 0;
scr(133) = 0;
scr(134) = 0;
scr(135) = 0;
scr(136) = 0;
scr(137) = 0;
scr(138) = 0;
scr(139) = 0;
scr(140) = 0;
scr(141) = 0;
scr(142) = 0;
scr(143) = 0;
scr(144) = 0;
scr(145) = 0;
scr(146) = 0;
scr(147) = 0;
scr(148) = 0;
scr(149) = 0;
scr(150) = 0;
scr(151) = 0;
scr(152) = 0;
scr(153) = 0;
scr(154) = 0;
scr(155) = 0;
scr(156) = 0;
scr(157) = 0;
scr(158) = 0;
scr(159) = 0;
scr(160) = 0;
scr(161) = 0;
scr(162) = 0;
scr(163) = 0;
scr(164) = 0;
scr(165) = 0;
scr(166) = 0;
scr(167) = 0;
scr(168) = 0;
scr(169) = 0;
scr(170) = 0;
scr(171) = 0;
scr(172) = 0;
scr(173) = 0;
scr(174) = 0;
scr(175) = 0;
scr(176) = 0;
scr(177) = 0;
scr(178) = 0;
scr(179) = 0;
scr(180) = 0;
scr(181) = 0;
scr(182) = 0;
scr(183) = 0;
scr(184) = 0;
scr(185) = (t85 - dt6 * t94 / 0.360e3 + t99);
scr(186) = 0;
scr(187) = 0;
scr(188) = 0;
scr(189) = 0;
scr(190) = 0;
scr(191) = 0;
scr(192) = 0;
scr(193) = t108;
scr(194) = (-t109 - dt6 * (cdz2i * t121 + t114 - t124) / 0.360e3 - dt8 * (cdy2i * t129 + cdz2i * (t131 - t134 - t137 + t138) + t142) / 0.20160e5);
scr(195) = t108;
scr(196) = 0;
scr(197) = 0;
scr(198) = 0;
scr(199) = 0;
scr(200) = 0;
scr(201) = (t147 - dt6 * t154 / 0.360e3 + t159);
scr(202) = (-t161 - dt6 * (cdz2i * t170 - t124 + t164) / 0.360e3 - dt8 * (cdx2i * t176 + cdz2i * (t178 - t181 - t183 + t138) + t142) / 0.20160e5);
scr(203) = (-t190 - dt4 * t198 / 0.12e2 - dt8 * (cdz2i * (t201 + t215 - t216 - t217 - t9 + t220 + t223) + cdy2i * (-t226 - t217 + t201) + cdx2i * (-t229 - t216 + t201)) / 0.20160e5 + dt6 * (cdx2i * t241 + cdy2i * t238 + cdz2i * t254) / 0.360e3);
scr(204) = (-t161 + dt8 * (-cdx2i * t176 - t142 + t263) / 0.20160e5 - dt6 * (t164 - t268 - t124) / 0.360e3);
scr(205) = t278;
scr(206) = 0;
scr(207) = 0;
scr(208) = 0;
scr(209) = 0;
scr(210) = 0;
scr(211) = t108;
scr(212) = (-t109 + dt8 * (-cdy2i * t129 - t142 + t282) / 0.20160e5 - dt6 * (t114 - t287 - t124) / 0.360e3);
scr(213) = t108;
scr(214) = 0;
scr(215) = 0;
scr(216) = 0;
scr(217) = 0;
scr(218) = 0;
scr(219) = 0;
scr(220) = 0;
scr(221) = t297;
scr(222) = 0;
scr(223) = 0;
scr(224) = 0;
scr(225) = 0;
scr(226) = 0;
scr(227) = 0;
scr(228) = 0;
scr(229) = 0;
scr(230) = 0;
scr(231) = 0;
scr(232) = 0;
scr(233) = 0;
scr(234) = 0;
scr(235) = 0;
scr(236) = 0;
scr(237) = 0;
scr(238) = 0;
scr(239) = 0;
scr(240) = 0;
scr(241) = 0;
scr(242) = 0;
scr(243) = 0;
scr(244) = 0;
scr(245) = 0;
scr(246) = 0;
scr(247) = 0;
scr(248) = 0;
scr(249) = 0;
scr(250) = 0;
scr(251) = 0;
scr(252) = 0;
scr(253) = 0;
scr(254) = 0;
scr(255) = 0;
scr(256) = 0;
scr(257) = t305;
scr(258) = 0;
scr(259) = 0;
scr(260) = 0;
scr(261) = 0;
scr(262) = 0;
scr(263) = 0;
scr(264) = 0;
scr(265) = t308;
scr(266) = (-t109 - dt6 * (cdy2i * t315 + cdz2i * t313 - t124) / 0.360e3 - dt8 * (cdz2i * t327 + cdy2i * (t329 - t134 - t137 + t138) + t333) / 0.20160e5);
scr(267) = t308;
scr(268) = 0;
scr(269) = 0;
scr(270) = 0;
scr(271) = 0;
scr(272) = 0;
scr(273) = t340;
scr(274) = (-dt8 * (cdx2i * t342 + cdy2i * t344 + cdz2i * t351) / 0.20160e5 + t356);
scr(275) = (dt6 * (cdy2i * t361 + cdz2i * t369 - t358) / 0.360e3 + t374 - dt8 * (cdy2i * (-t375 + t376 + t377 + t223 - t378 + t220 + t215) + cdz2i * (-t375 + t376 + t377 - t381 + t385 + t388 + t390) + cdx2i * (-t393 + t376 + t377)) / 0.20160e5);
scr(276) = (dt8 * (-cdx2i * t342 + t402 + t406) / 0.20160e5 + t356);
scr(277) = t340;
scr(278) = 0;
scr(279) = 0;
scr(280) = 0;
scr(281) = t418;
scr(282) = (-t161 - dt6 * (cdx2i * t424 + cdz2i * t422 - t124) / 0.360e3 - dt8 * (cdz2i * t436 + cdx2i * (t438 - t181 - t183 + t138) + t442) / 0.20160e5);
scr(283) = (dt6 * (cdx2i * t449 + cdz2i * t455 - t358) / 0.360e3 + t460 - dt8 * (cdx2i * (-t375 + t461 + t462 - t463 + t223 + t220 + t215) + cdz2i * (-t375 + t461 + t462 + t468 - t469 + t471 + t473) + cdy2i * (-t393 + t461 + t462)) / 0.20160e5);
scr(284) = (dt4 * t489 / 0.12e2 + t492 - dt6 * (cdx2i * t501 + cdy2i * t497 + cdz2i * t518) / 0.360e3 + dt8 * (cdz2i * (t523 - t528 - t530 - t532 - t75 - t76 + t74 + t533 + t534) + cdy2i * (t523 + t533 + t534 - t375 - t279 + t134 + t137 - t280) + cdx2i * (t523 + t533 + t534 - t375 - t260 + t181 + t183 - t261)) / 0.20160e5);
scr(285) = (dt8 * (cdx2i * (t545 + t546 + t375 - t220 - t223 + t463 - t215) + t550 + t555) / 0.20160e5 - t563 + t460);
scr(286) = (-t161 - dt8 * (cdx2i * t565 + cdz2i * t568 + t442) / 0.20160e5 + t577);
scr(287) = t418;
scr(288) = 0;
scr(289) = 0;
scr(290) = 0;
scr(291) = t340;
scr(292) = (dt8 * (-cdy2i * t344 + t406 + t581) / 0.20160e5 + t356);
scr(293) = (dt8 * (cdy2i * (t586 + t587 + t375 - t223 - t220 + t378 - t215) + t591 + t596) / 0.20160e5 - t604 + t374);
scr(294) = t612;
scr(295) = t340;
scr(296) = 0;
scr(297) = 0;
scr(298) = 0;
scr(299) = 0;
scr(300) = 0;
scr(301) = t308;
scr(302) = (-t109 - dt8 * (cdy2i * t613 + cdz2i * t616 + t333) / 0.20160e5 + t625);
scr(303) = t308;
scr(304) = 0;
scr(305) = 0;
scr(306) = 0;
scr(307) = 0;
scr(308) = 0;
scr(309) = 0;
scr(310) = 0;
scr(311) = t305;
scr(312) = 0;
scr(313) = 0;
scr(314) = 0;
scr(315) = 0;
scr(316) = 0;
scr(317) = 0;
scr(318) = 0;
scr(319) = 0;
scr(320) = 0;
scr(321) = 0;
scr(322) = 0;
scr(323) = 0;
scr(324) = 0;
scr(325) = 0;
scr(326) = 0;
scr(327) = 0;
scr(328) = 0;
scr(329) = t636;
scr(330) = 0;
scr(331) = 0;
scr(332) = 0;
scr(333) = 0;
scr(334) = 0;
scr(335) = 0;
scr(336) = 0;
scr(337) = t643;
scr(338) = (dt4 * t649 / 0.12e2 - dt6 * (cdy2i * t654 - t656 - t657) / 0.360e3 - dt8 * (t662 + t664 - cdy2i * (t326 - t665 - t323 - t666)) / 0.20160e5 + t672);
scr(339) = t643;
scr(340) = 0;
scr(341) = 0;
scr(342) = 0;
scr(343) = 0;
scr(344) = 0;
scr(345) = (t674 - dt6 * t679 / 0.360e3 + t683);
scr(346) = (-t685 - dt6 * (cdx2i * t687 + cdy2i * t690 - t124) / 0.360e3 - dt8 * (cdx2i * t696 + cdy2i * (t698 - t348 - t350 + t138) + t333) / 0.20160e5);
scr(347) = (-t705 - dt4 * t711 / 0.12e2 - dt8 * (cdy2i * (t714 + t390 - t715 - t716 - t301 + t388 + t385) + cdz2i * (-t719 - t716 + t714) + cdx2i * (-t722 - t715 + t714)) / 0.20160e5 + dt6 * (cdx2i * t733 + cdy2i * t736 + cdz2i * t730) / 0.360e3);
scr(348) = (-t685 + dt8 * (-cdx2i * t696 - t333 + t743) / 0.20160e5 - dt6 * t750 / 0.360e3);
scr(349) = t759;
scr(350) = 0;
scr(351) = 0;
scr(352) = 0;
scr(353) = t764;
scr(354) = (-t685 - dt6 * (cdx2i * t767 + cdy2i * t765 - t124) / 0.360e3 - dt8 * (cdy2i * t774 + cdx2i * (t776 - t348 - t350 + t138) + t442) / 0.20160e5);
scr(355) = (dt6 * (cdx2i * t784 + cdy2i * t786 - t358) / 0.360e3 + t791 - dt8 * (cdx2i * (t390 - t792 + t385 - t375 + t793 + t794 + t388) + cdy2i * (-t375 + t793 + t794 - t797 + t468 + t471 + t473) + cdz2i * (-t393 + t793 + t794)) / 0.20160e5);
scr(356) = (dt4 * t811 / 0.12e2 + t814 - dt6 * (cdx2i * t821 + cdy2i * t825 + cdz2i * t818) / 0.360e3 + dt8 * (cdy2i * (t830 + t831 + t832 - t528 - t530 - t532 + t326 - t665 - t323) + cdz2i * (-t279 - t375 - t280 + t134 + t137 + t830 + t831 + t832) + cdx2i * (-t403 - t375 - t404 + t348 + t350 + t830 + t831 + t832)) / 0.20160e5);
scr(357) = (dt8 * (cdx2i * (t843 + t844 + t375 - t388 - t385 + t792 - t390) + t848 + t850) / 0.20160e5 - t858 + t791);
scr(358) = (-t685 - dt8 * (cdx2i * t860 + cdy2i * t862 + t442) / 0.20160e5 + t871);
scr(359) = t764;
scr(360) = 0;
scr(361) = t882;
scr(362) = (dt4 * t886 / 0.12e2 - dt6 * (cdx2i * t890 - t892 - t893) / 0.360e3 - dt8 * (t898 + t900 - cdx2i * (t435 - t773 - t430 - t901)) / 0.20160e5 + t907);
scr(363) = (dt6 * (cdx2i * t915 + cdy2i * t913 + cdz2i * t911) / 0.360e3 - t920 - dt4 * t924 / 0.12e2 - dt8 * (cdx2i * (t927 + t473 - t928 - t929 - t414 + t471 + t468) + cdz2i * (-t932 - t929 + t927) + cdy2i * (-t935 - t928 + t927)) / 0.20160e5);
scr(364) = (dt4 * t945 / 0.12e2 + t948 - dt6 * (cdx2i * t956 + cdy2i * t953 + cdz2i * t951) / 0.360e3 + dt8 * (cdx2i * (t961 + t962 + t963 - t532 - t773 - t430 - t528 - t530 + t435) + cdz2i * (-t260 - t375 - t261 + t181 + t183 + t961 + t962 + t963) + cdy2i * (-t404 - t375 - t403 + t348 + t350 + t961 + t962 + t963)) / 0.20160e5);
scr(365) = (dt4 * (cdz2i * (0.91e2 / 0.8e1 * cdz2i + t975 + t976) + cdy2i * (t979 + 0.91e2 / 0.8e1 * cdy2i + t976) + cdx2i * (t979 + t975 + 0.91e2 / 0.8e1 * cdx2i)) / 0.12e2 - 0.205e3 / 0.72e2 * dt2 * (cdz2i + cdy2i + cdx2i) - dt6 * (cdx2i * t1007 + cdy2i * t1003 + cdz2i * t999) / 0.360e3 + dt8 * (cdz2i * (-t215 - t223 - t220 + t1012 + t1013 + t1014 + t1015 + t1016 + t1017) + cdy2i * (t592 - t390 - t385 - t388 + t594 + t593 + t1015 + t1016 + t1017) + cdx2i * (-t468 - t471 + t553 + t552 + t551 - t473 + t1015 + t1016 + t1017)) / 0.20160e5 + 0.2e1);
scr(366) = (t948 + dt4 * t945 / 0.12e2 - dt8 * (cdx2i * (t1030 + t1031 + t1032 + t532 + t530 + t528 + t567 + t773 + t430) + cdz2i * (t546 + t545 + t375 + t1030 + t1031 + t1032) + cdy2i * (t844 + t843 + t375 + t1030 + t1031 + t1032)) / 0.20160e5 + dt6 * (-cdx2i * t956 - cdy2i * t953 - cdz2i * t951) / 0.360e3);
scr(367) = (-dt4 * t924 / 0.12e2 - t920 + dt8 * (cdz2i * (t932 + t929 + t1051) + cdy2i * (t928 + t935 + t1051) + cdx2i * (t928 + t929 + t414 + t1051 + t551 + t552 + t553)) / 0.20160e5 - dt6 * (-cdx2i * t915 - cdy2i * t913 - cdz2i * t911) / 0.360e3);
scr(368) = (t907 + dt8 * (cdx2i * (-t773 - t430 - t901 - t567) - t898 - t900) / 0.20160e5 - dt6 * (cdx2i * t890 - t892 - t893) / 0.360e3 + dt4 * t886 / 0.12e2);
scr(369) = t882;
scr(370) = 0;
scr(371) = t764;
scr(372) = (-t685 + dt8 * (-cdy2i * t774 + t1080 - t442) / 0.20160e5 + dt6 * t869 / 0.360e3);
scr(373) = (dt8 * (cdy2i * (t844 + t843 + t375 - t468 - t471 + t797 - t473) + t848 + t1091) / 0.20160e5 - t858 + t791);
scr(374) = (t814 + dt4 * t811 / 0.12e2 - dt8 * (cdy2i * (t1098 + t532 + t530 + t528 + t615 + t665 + t323 + t1099 + t1100) + cdz2i * (t1098 + t586 + t587 + t375 + t1099 + t1100) + cdx2i * (t1098 + t844 + t843 + t375 + t1099 + t1100)) / 0.20160e5 + dt6 * (-cdx2i * t821 - cdy2i * t825 - cdz2i * t818) / 0.360e3);
scr(375) = (-t858 + t791 + dt8 * (t848 + t1091 + t850) / 0.20160e5);
scr(376) = (-t685 + dt8 * (-cdy2i * t862 + t1080 - t442) / 0.20160e5 + t871);
scr(377) = t764;
scr(378) = 0;
scr(379) = 0;
scr(380) = 0;
scr(381) = t759;
scr(382) = (-t685 - dt8 * (cdx2i * t1127 - cdy2i * t742 + t333) / 0.20160e5 + t1133);
scr(383) = (-dt4 * t711 / 0.12e2 - t705 + dt8 * (cdz2i * (t719 + t716 + t1137) + cdx2i * (t715 + t722 + t1137) + cdy2i * (t715 + t716 + t301 + t1137 + t592 + t593 + t594)) / 0.20160e5 - dt6 * (-cdx2i * t733 - cdy2i * t736 - cdz2i * t730) / 0.360e3);
scr(384) = (-t685 + dt8 * (-cdx2i * t1127 - t333 + t743) / 0.20160e5 + t1133);
scr(385) = (t674 - dt6 * t679 / 0.360e3 + t683);
scr(386) = 0;
scr(387) = 0;
scr(388) = 0;
scr(389) = 0;
scr(390) = 0;
scr(391) = t643;
scr(392) = (t672 + dt8 * (cdy2i * (-t665 - t323 - t666 - t615) - t662 - t664) / 0.20160e5 - dt6 * (cdy2i * t654 - t656 - t657) / 0.360e3 + dt4 * t649 / 0.12e2);
scr(393) = t643;
scr(394) = 0;
scr(395) = 0;
scr(396) = 0;
scr(397) = 0;
scr(398) = 0;
scr(399) = 0;
scr(400) = 0;
scr(401) = t636;
scr(402) = 0;
scr(403) = 0;
scr(404) = 0;
scr(405) = 0;
scr(406) = 0;
scr(407) = 0;
scr(408) = 0;
scr(409) = 0;
scr(410) = 0;
scr(411) = 0;
scr(412) = 0;
scr(413) = 0;
scr(414) = 0;
scr(415) = 0;
scr(416) = 0;
scr(417) = 0;
scr(418) = 0;
scr(419) = t305;
scr(420) = 0;
scr(421) = 0;
scr(422) = 0;
scr(423) = 0;
scr(424) = 0;
scr(425) = 0;
scr(426) = 0;
scr(427) = t308;
scr(428) = (-t109 + dt8 * (-cdz2i * t327 + t1174 - t333) / 0.20160e5 + dt6 * t623 / 0.360e3);
scr(429) = t308;
scr(430) = 0;
scr(431) = 0;
scr(432) = 0;
scr(433) = 0;
scr(434) = 0;
scr(435) = t340;
scr(436) = (dt8 * (-cdz2i * t351 + t402 + t581) / 0.20160e5 + t356);
scr(437) = (dt8 * (cdz2i * (t587 + t586 + t375 - t388 - t385 + t381 - t390) + t591 + t1190) / 0.20160e5 - t604 + t374);
scr(438) = t612;
scr(439) = t340;
scr(440) = 0;
scr(441) = 0;
scr(442) = 0;
scr(443) = t418;
scr(444) = (-t161 + dt8 * (-cdz2i * t436 + t1195 - t442) / 0.20160e5 + dt6 * t575 / 0.360e3);
scr(445) = (dt8 * (cdz2i * (t546 + t545 + t375 - t471 - t468 + t469 - t473) + t550 + t1206) / 0.20160e5 - t563 + t460);
scr(446) = (t492 + dt4 * t489 / 0.12e2 - dt8 * (cdz2i * (t532 + t530 + t528 + t1213 + t75 + t76 + t1214 + t1215 + t1216) + cdy2i * (t586 + t587 + t375 + t1214 + t1215 + t1216) + cdx2i * (t545 + t546 + t375 + t1214 + t1215 + t1216)) / 0.20160e5 + dt6 * (-cdx2i * t501 - cdy2i * t497 - cdz2i * t518) / 0.360e3);
scr(447) = (-t563 + t460 + dt8 * (t550 + t1206 + t555) / 0.20160e5);
scr(448) = (-t161 + dt8 * (-cdz2i * t568 + t1195 - t442) / 0.20160e5 + t577);
scr(449) = t418;
scr(450) = 0;
scr(451) = 0;
scr(452) = 0;
scr(453) = t340;
scr(454) = t612;
scr(455) = (-t604 + t374 + dt8 * (t591 + t1190 + t596) / 0.20160e5);
scr(456) = (t356 + dt8 * (t581 + t402 + t406) / 0.20160e5);
scr(457) = t340;
scr(458) = 0;
scr(459) = 0;
scr(460) = 0;
scr(461) = 0;
scr(462) = 0;
scr(463) = t308;
scr(464) = (-t109 + dt8 * (-cdz2i * t616 + t1174 - t333) / 0.20160e5 + t625);
scr(465) = t308;
scr(466) = 0;
scr(467) = 0;
scr(468) = 0;
scr(469) = 0;
scr(470) = 0;
scr(471) = 0;
scr(472) = 0;
scr(473) = t305;
scr(474) = 0;
scr(475) = 0;
scr(476) = 0;
scr(477) = 0;
scr(478) = 0;
scr(479) = 0;
scr(480) = 0;
scr(481) = 0;
scr(482) = 0;
scr(483) = 0;
scr(484) = 0;
scr(485) = 0;
scr(486) = 0;
scr(487) = 0;
scr(488) = 0;
scr(489) = 0;
scr(490) = 0;
scr(491) = 0;
scr(492) = 0;
scr(493) = 0;
scr(494) = 0;
scr(495) = 0;
scr(496) = 0;
scr(497) = 0;
scr(498) = 0;
scr(499) = 0;
scr(500) = 0;
scr(501) = 0;
scr(502) = 0;
scr(503) = 0;
scr(504) = 0;
scr(505) = 0;
scr(506) = 0;
scr(507) = 0;
scr(508) = 0;
scr(509) = t297;
scr(510) = 0;
scr(511) = 0;
scr(512) = 0;
scr(513) = 0;
scr(514) = 0;
scr(515) = 0;
scr(516) = 0;
scr(517) = t108;
scr(518) = (-t109 - dt8 * (cdy2i * t1256 - cdz2i * t281 + t142) / 0.20160e5 + t1264);
scr(519) = t108;
scr(520) = 0;
scr(521) = 0;
scr(522) = 0;
scr(523) = 0;
scr(524) = 0;
scr(525) = t278;
scr(526) = (-t161 - dt8 * (cdx2i * t1267 - cdz2i * t262 + t142) / 0.20160e5 + t1275);
scr(527) = (-dt4 * t198 / 0.12e2 - t190 + dt8 * (cdy2i * (t1279 + t226 + t217) + cdx2i * (t1279 + t216 + t229) + cdz2i * (t216 + t217 + t9 + t1279 + t1014 + t1013 + t1012)) / 0.20160e5 - dt6 * (-cdx2i * t241 - cdy2i * t238 - cdz2i * t254) / 0.360e3);
scr(528) = (-t161 + dt8 * (-cdx2i * t1267 - t142 + t263) / 0.20160e5 + t1275);
scr(529) = (t147 - dt6 * t154 / 0.360e3 + t159);
scr(530) = 0;
scr(531) = 0;
scr(532) = 0;
scr(533) = 0;
scr(534) = 0;
scr(535) = t108;
scr(536) = (-t109 + dt8 * (-cdy2i * t1256 - t142 + t282) / 0.20160e5 + t1264);
scr(537) = t108;
scr(538) = 0;
scr(539) = 0;
scr(540) = 0;
scr(541) = 0;
scr(542) = 0;
scr(543) = 0;
scr(544) = 0;
scr(545) = (t85 - dt6 * t94 / 0.360e3 + t99);
scr(546) = 0;
scr(547) = 0;
scr(548) = 0;
scr(549) = 0;
scr(550) = 0;
scr(551) = 0;
scr(552) = 0;
scr(553) = 0;
scr(554) = 0;
scr(555) = 0;
scr(556) = 0;
scr(557) = 0;
scr(558) = 0;
scr(559) = 0;
scr(560) = 0;
scr(561) = 0;
scr(562) = 0;
scr(563) = 0;
scr(564) = 0;
scr(565) = 0;
scr(566) = 0;
scr(567) = 0;
scr(568) = 0;
scr(569) = 0;
scr(570) = 0;
scr(571) = 0;
scr(572) = 0;
scr(573) = 0;
scr(574) = 0;
scr(575) = 0;
scr(576) = 0;
scr(577) = 0;
scr(578) = 0;
scr(579) = 0;
scr(580) = 0;
scr(581) = 0;
scr(582) = 0;
scr(583) = 0;
scr(584) = 0;
scr(585) = 0;
scr(586) = 0;
scr(587) = 0;
scr(588) = 0;
scr(589) = 0;
scr(590) = 0;
scr(591) = 0;
scr(592) = 0;
scr(593) = 0;
scr(594) = 0;
scr(595) = 0;
scr(596) = 0;
scr(597) = 0;
scr(598) = 0;
scr(599) = t22;
scr(600) = 0;
scr(601) = 0;
scr(602) = 0;
scr(603) = 0;
scr(604) = 0;
scr(605) = 0;
scr(606) = 0;
scr(607) = t30;
scr(608) = (t83 + dt8 * (cdz2i * (-t75 - t76 - t77 - t1213) - t67 - t69) / 0.20160e5 + dt6 * (-cdz2i * t57 + t60 + t62) / 0.360e3 + dt4 * t40 / 0.12e2);
scr(609) = t30;
scr(610) = 0;
scr(611) = 0;
scr(612) = 0;
scr(613) = 0;
scr(614) = 0;
scr(615) = 0;
scr(616) = 0;
scr(617) = t22;
scr(618) = 0;
scr(619) = 0;
scr(620) = 0;
scr(621) = 0;
scr(622) = 0;
scr(623) = 0;
scr(624) = 0;
scr(625) = 0;
scr(626) = 0;
scr(627) = 0;
scr(628) = 0;
scr(629) = 0;
scr(630) = 0;
scr(631) = 0;
scr(632) = 0;
scr(633) = 0;
scr(634) = 0;
scr(635) = 0;
scr(636) = 0;
scr(637) = 0;
scr(638) = 0;
scr(639) = 0;
scr(640) = 0;
scr(641) = 0;
scr(642) = 0;
scr(643) = 0;
scr(644) = 0;
scr(645) = 0;
scr(646) = 0;
scr(647) = 0;
scr(648) = 0;
scr(649) = 0;
scr(650) = 0;
scr(651) = 0;
scr(652) = 0;
scr(653) = 0;
scr(654) = 0;
scr(655) = 0;
scr(656) = 0;
scr(657) = 0;
scr(658) = 0;
scr(659) = 0;
scr(660) = 0;
scr(661) = 0;
scr(662) = 0;
scr(663) = 0;
scr(664) = 0;
scr(665) = 0;
scr(666) = 0;
scr(667) = 0;
scr(668) = 0;
scr(669) = 0;
scr(670) = 0;
scr(671) = 0;
scr(672) = 0;
scr(673) = 0;
scr(674) = 0;
scr(675) = 0;
scr(676) = 0;
scr(677) = 0;
scr(678) = 0;
scr(679) = 0;
scr(680) = 0;
scr(681) = 0;
scr(682) = 0;
scr(683) = 0;
scr(684) = 0;
scr(685) = 0;
scr(686) = 0;
scr(687) = 0;
scr(688) = 0;
scr(689) = t12;
scr(690) = 0;
scr(691) = 0;
scr(692) = 0;
scr(693) = 0;
scr(694) = 0;
scr(695) = 0;
scr(696) = 0;
scr(697) = 0;
scr(698) = 0;
scr(699) = 0;
scr(700) = 0;
scr(701) = 0;
scr(702) = 0;
scr(703) = 0;
scr(704) = 0;
scr(705) = 0;
scr(706) = 0;
scr(707) = 0;
scr(708) = 0;
scr(709) = 0;
scr(710) = 0;
scr(711) = 0;
scr(712) = 0;
scr(713) = 0;
scr(714) = 0;
scr(715) = 0;
scr(716) = 0;
scr(717) = 0;
scr(718) = 0;
scr(719) = 0;
scr(720) = 0;
scr(721) = 0;
scr(722) = 0;
scr(723) = 0;
scr(724) = 0;
scr(725) = 0;
scr(726) = 0;
scr(727) = 0;
scr(728) = 0;
scr(729) = 0;
        ! #Include "../include/defineStencilVariables2dOrder2Rectangular.h"
     ! ! now copy results to small arrays
     ! stencilWidth = orderOfAccuracy+1;
     ! numStencilCoeff = stencilWidth**nd;
     ! do m=1,numStencilCoeff
     !   scr(m) = sc(m,0,0)
     ! end do
     ! if( t.le.dt )then
     !   do m=1,numStencilCoeff
     !     if( scr(m).eq.0 )then
     !       write(*,'("c(",i3,")=",e10.2)') m,scr(m)
     !     end if
     !   end do
     ! end if
   if( gridIsImplicit.eq.0 )then 
     ! ------- EXPLICIT update the solution ---------
       if( orderInTime.eq.2 .and. orderOfAccuracy.gt.2 )then
         stop 2222
         ! if( addForcing.eq.0 )then
         !   updateWaveOpt(3,8,2,rectangular,NOFORCING)
         ! else 
         !   updateWaveOpt(3,8,2,rectangular,FORCING)
         ! end if
       else
         if( addForcing.eq.0 )then
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=3 order=8 orderInTime=8, grid=rectangular... t=",e10.2)') t
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
! Stencil: nd=3, orderOfAccuracy=8, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ scr( 41)*u(i1+0,i2+0,i3-4,m)                                                                                                                                + scr(113)*u(i1+0,i2-1,i3-3,m)                                                                                                                                + scr(121)*u(i1-1,i2+0,i3-3,m) + scr(122)*u(i1+0,i2+0,i3-3,m) + scr(123)*u(i1+1,i2+0,i3-3,m)                                                                                                + scr(131)*u(i1+0,i2+1,i3-3,m)                                                                                                                                + scr(185)*u(i1+0,i2-2,i3-2,m)                                                                                                                                + scr(193)*u(i1-1,i2-1,i3-2,m) + scr(194)*u(i1+0,i2-1,i3-2,m) + scr(195)*u(i1+1,i2-1,i3-2,m)                                                                                                + scr(201)*u(i1-2,i2+0,i3-2,m) + scr(202)*u(i1-1,i2+0,i3-2,m) + scr(203)*u(i1+0,i2+0,i3-2,m) + scr(204)*u(i1+1,i2+0,i3-2,m) + scr(205)*u(i1+2,i2+0,i3-2,m)                                                                + scr(211)*u(i1-1,i2+1,i3-2,m) + scr(212)*u(i1+0,i2+1,i3-2,m) + scr(213)*u(i1+1,i2+1,i3-2,m)                                                                                                + scr(221)*u(i1+0,i2+2,i3-2,m)                                                                                                                                + scr(257)*u(i1+0,i2-3,i3-1,m)                                                                                                                                + scr(265)*u(i1-1,i2-2,i3-1,m) + scr(266)*u(i1+0,i2-2,i3-1,m) + scr(267)*u(i1+1,i2-2,i3-1,m)                                                                                                + scr(273)*u(i1-2,i2-1,i3-1,m) + scr(274)*u(i1-1,i2-1,i3-1,m) + scr(275)*u(i1+0,i2-1,i3-1,m) + scr(276)*u(i1+1,i2-1,i3-1,m) + scr(277)*u(i1+2,i2-1,i3-1,m)                                                                + scr(281)*u(i1-3,i2+0,i3-1,m) + scr(282)*u(i1-2,i2+0,i3-1,m) + scr(283)*u(i1-1,i2+0,i3-1,m) + scr(284)*u(i1+0,i2+0,i3-1,m) + scr(285)*u(i1+1,i2+0,i3-1,m) + scr(286)*u(i1+2,i2+0,i3-1,m) + scr(287)*u(i1+3,i2+0,i3-1,m)                                + scr(291)*u(i1-2,i2+1,i3-1,m) + scr(292)*u(i1-1,i2+1,i3-1,m) + scr(293)*u(i1+0,i2+1,i3-1,m) + scr(294)*u(i1+1,i2+1,i3-1,m) + scr(295)*u(i1+2,i2+1,i3-1,m)                                                                + scr(301)*u(i1-1,i2+2,i3-1,m) + scr(302)*u(i1+0,i2+2,i3-1,m) + scr(303)*u(i1+1,i2+2,i3-1,m)                                                                                                + scr(311)*u(i1+0,i2+3,i3-1,m)                                                                                                                                + scr(329)*u(i1+0,i2-4,i3+0,m)                                                                                                                                + scr(337)*u(i1-1,i2-3,i3+0,m) + scr(338)*u(i1+0,i2-3,i3+0,m) + scr(339)*u(i1+1,i2-3,i3+0,m)                                                                                                + scr(345)*u(i1-2,i2-2,i3+0,m) + scr(346)*u(i1-1,i2-2,i3+0,m) + scr(347)*u(i1+0,i2-2,i3+0,m) + scr(348)*u(i1+1,i2-2,i3+0,m) + scr(349)*u(i1+2,i2-2,i3+0,m)                                                                + scr(353)*u(i1-3,i2-1,i3+0,m) + scr(354)*u(i1-2,i2-1,i3+0,m) + scr(355)*u(i1-1,i2-1,i3+0,m) + scr(356)*u(i1+0,i2-1,i3+0,m) + scr(357)*u(i1+1,i2-1,i3+0,m) + scr(358)*u(i1+2,i2-1,i3+0,m) + scr(359)*u(i1+3,i2-1,i3+0,m)                                + scr(361)*u(i1-4,i2+0,i3+0,m) + scr(362)*u(i1-3,i2+0,i3+0,m) + scr(363)*u(i1-2,i2+0,i3+0,m) + scr(364)*u(i1-1,i2+0,i3+0,m) + scr(365)*u(i1+0,i2+0,i3+0,m) + scr(366)*u(i1+1,i2+0,i3+0,m) + scr(367)*u(i1+2,i2+0,i3+0,m) + scr(368)*u(i1+3,i2+0,i3+0,m) + scr(369)*u(i1+4,i2+0,i3+0,m)+ scr(371)*u(i1-3,i2+1,i3+0,m) + scr(372)*u(i1-2,i2+1,i3+0,m) + scr(373)*u(i1-1,i2+1,i3+0,m) + scr(374)*u(i1+0,i2+1,i3+0,m) + scr(375)*u(i1+1,i2+1,i3+0,m) + scr(376)*u(i1+2,i2+1,i3+0,m) + scr(377)*u(i1+3,i2+1,i3+0,m)                                + scr(381)*u(i1-2,i2+2,i3+0,m) + scr(382)*u(i1-1,i2+2,i3+0,m) + scr(383)*u(i1+0,i2+2,i3+0,m) + scr(384)*u(i1+1,i2+2,i3+0,m) + scr(385)*u(i1+2,i2+2,i3+0,m)                                                                + scr(391)*u(i1-1,i2+3,i3+0,m) + scr(392)*u(i1+0,i2+3,i3+0,m) + scr(393)*u(i1+1,i2+3,i3+0,m)                                                                                                + scr(401)*u(i1+0,i2+4,i3+0,m)                                                                                                                                + scr(419)*u(i1+0,i2-3,i3+1,m)                                                                                                                                + scr(427)*u(i1-1,i2-2,i3+1,m) + scr(428)*u(i1+0,i2-2,i3+1,m) + scr(429)*u(i1+1,i2-2,i3+1,m)                                                                                                + scr(435)*u(i1-2,i2-1,i3+1,m) + scr(436)*u(i1-1,i2-1,i3+1,m) + scr(437)*u(i1+0,i2-1,i3+1,m) + scr(438)*u(i1+1,i2-1,i3+1,m) + scr(439)*u(i1+2,i2-1,i3+1,m)                                                                + scr(443)*u(i1-3,i2+0,i3+1,m) + scr(444)*u(i1-2,i2+0,i3+1,m) + scr(445)*u(i1-1,i2+0,i3+1,m) + scr(446)*u(i1+0,i2+0,i3+1,m) + scr(447)*u(i1+1,i2+0,i3+1,m) + scr(448)*u(i1+2,i2+0,i3+1,m) + scr(449)*u(i1+3,i2+0,i3+1,m)                                + scr(453)*u(i1-2,i2+1,i3+1,m) + scr(454)*u(i1-1,i2+1,i3+1,m) + scr(455)*u(i1+0,i2+1,i3+1,m) + scr(456)*u(i1+1,i2+1,i3+1,m) + scr(457)*u(i1+2,i2+1,i3+1,m)                                                                + scr(463)*u(i1-1,i2+2,i3+1,m) + scr(464)*u(i1+0,i2+2,i3+1,m) + scr(465)*u(i1+1,i2+2,i3+1,m)                                                                                                + scr(473)*u(i1+0,i2+3,i3+1,m)                                                                                                                                + scr(509)*u(i1+0,i2-2,i3+2,m)                                                                                                                                + scr(517)*u(i1-1,i2-1,i3+2,m) + scr(518)*u(i1+0,i2-1,i3+2,m) + scr(519)*u(i1+1,i2-1,i3+2,m)                                                                                                + scr(525)*u(i1-2,i2+0,i3+2,m) + scr(526)*u(i1-1,i2+0,i3+2,m) + scr(527)*u(i1+0,i2+0,i3+2,m) + scr(528)*u(i1+1,i2+0,i3+2,m) + scr(529)*u(i1+2,i2+0,i3+2,m)                                                                + scr(535)*u(i1-1,i2+1,i3+2,m) + scr(536)*u(i1+0,i2+1,i3+2,m) + scr(537)*u(i1+1,i2+1,i3+2,m)                                                                                                + scr(545)*u(i1+0,i2+2,i3+2,m)                                                                                                                                + scr(599)*u(i1+0,i2-1,i3+3,m)                                                                                                                                + scr(607)*u(i1-1,i2+0,i3+3,m) + scr(608)*u(i1+0,i2+0,i3+3,m) + scr(609)*u(i1+1,i2+0,i3+3,m)                                                                                                + scr(617)*u(i1+0,i2+1,i3+3,m)                                                                                                                                + scr(689)*u(i1+0,i2+0,i3+4,m)                                                                                                                                
               end do
               end do
               end do
             ! endLoopsMask()
         else
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=3 order=8 orderInTime=8, grid=rectangular... t=",e10.2)') t
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
                        ! Correct forcing for sixth-order ME in 3D
                          call ogDeriv(ep, 6,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evtttttt(m) )
                          call ogDeriv(ep, 0,6,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxxxxxx(m) )
                          call ogDeriv(ep, 0,0,6,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evyyyyyy(m) )
                          call ogDeriv(ep, 0,0,0,6, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evzzzzzz(m) )
                          call ogDeriv(ep, 0,4,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxxxxyy(m) )
                          call ogDeriv(ep, 0,2,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxxyyyy(m) )
                          call ogDeriv(ep, 0,4,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxxxxzz(m) )
                          call ogDeriv(ep, 0,2,0,4, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxxzzzz(m) )
                          call ogDeriv(ep, 0,0,4,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evyyyyzz(m) )
                          call ogDeriv(ep, 0,0,2,4, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evyyzzzz(m) )
                          call ogDeriv(ep, 0,2,2,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxxyyzz(m) )
                        fv(m) = fv(m) + (dtSq**2/360.)*evtttttt(m) - (cdtPow6By360/dtSq)*( evxxxxxx(m) + evyyyyyy(m) + evzzzzzz(m) + 3.*(evxxxxyy(m) + evxxyyyy(m) + evxxxxzz(m) + evxxzzzz(m) + evyyyyzz(m) + evyyzzzz(m) ) + 6.*evxxyyzz(m)  )
                        ! Correct forcing for eighth-order ME in 3D
                          call ogDeriv(ep, 8,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uet8 )
                          call ogDeriv(ep, 0,8,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uex8 )
                          call ogDeriv(ep, 0,0,8,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uey8 )
                          call ogDeriv(ep, 0,0,0,8, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uez8 )
                          call ogDeriv(ep, 0,6,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uex6y2 )
                          call ogDeriv(ep, 0,4,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uex4y4 )
                          call ogDeriv(ep, 0,2,6,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uex2y6 )
                          call ogDeriv(ep, 0,6,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uex6z2 )
                          call ogDeriv(ep, 0,4,0,4, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uex4z4 )
                          call ogDeriv(ep, 0,2,0,6, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uex2z6 )
                          call ogDeriv(ep, 0,0,6,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uey6z2 )
                          call ogDeriv(ep, 0,0,4,4, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uey4z4 )
                          call ogDeriv(ep, 0,0,2,6, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uey2z6 )
                          call ogDeriv(ep, 0,4,2,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uex4y2z2 )
                          call ogDeriv(ep, 0,2,4,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uex2y4z2 )
                          call ogDeriv(ep, 0,2,2,4, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,uex2y2z4 )
                        ! ( x*x _ y*y + z*z )^4 = x^8 + 4*x^6*y^2 + 4*x^6*z^2 + 6*x^4*y^4 + 12*x^4*y^2*z^2 + 6*x^4*z^4 + 4*x^2*y^6 + 12*x^2*y^4*z^2 + 12*x^2*y^2*z^4 + 4*x^2*z^6 + y^8 + 4*y^6*z^2 + 6*y^4*z^4 + 4*y^2*z^6 + z^8 
                        fv(m) = fv(m) + (dtSq**3/20160.)*uet8       - (cdtPow8By20160/dtSq)*( uex8 +       uey8 +       uez8         +  4.*(uex6y2 +    uex2y6 +    uex6z2 +    uex2z6 +    uey6z2 +    uey2z6 )    +  6.*( uex4y4 +   uex4z4 +   uey4z4 )   + 12.*( uex4y2z2 + uex2y4z2 + uex2y2z4 ) )
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
                             !fv(m) = fv(m) -( f(i1,i2,i3,freq) + cdtSqBy12*( cSq*(fxx23r(i1,i2,i3,freq) + fyy23r(i1,i2,i3,freq) + fzz23r(i1,i2,i3,freq)) - omega*omega*f(i1,i2,i3,freq)) )*coswt 
                     end do ! do freq  
                  else if( addForcing.ne.0 )then  
                     fv(m) = f(i1,i2,i3,0)
                  end if
! Stencil: nd=3, orderOfAccuracy=8, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ scr( 41)*u(i1+0,i2+0,i3-4,m)                                                                                                                                + scr(113)*u(i1+0,i2-1,i3-3,m)                                                                                                                                + scr(121)*u(i1-1,i2+0,i3-3,m) + scr(122)*u(i1+0,i2+0,i3-3,m) + scr(123)*u(i1+1,i2+0,i3-3,m)                                                                                                + scr(131)*u(i1+0,i2+1,i3-3,m)                                                                                                                                + scr(185)*u(i1+0,i2-2,i3-2,m)                                                                                                                                + scr(193)*u(i1-1,i2-1,i3-2,m) + scr(194)*u(i1+0,i2-1,i3-2,m) + scr(195)*u(i1+1,i2-1,i3-2,m)                                                                                                + scr(201)*u(i1-2,i2+0,i3-2,m) + scr(202)*u(i1-1,i2+0,i3-2,m) + scr(203)*u(i1+0,i2+0,i3-2,m) + scr(204)*u(i1+1,i2+0,i3-2,m) + scr(205)*u(i1+2,i2+0,i3-2,m)                                                                + scr(211)*u(i1-1,i2+1,i3-2,m) + scr(212)*u(i1+0,i2+1,i3-2,m) + scr(213)*u(i1+1,i2+1,i3-2,m)                                                                                                + scr(221)*u(i1+0,i2+2,i3-2,m)                                                                                                                                + scr(257)*u(i1+0,i2-3,i3-1,m)                                                                                                                                + scr(265)*u(i1-1,i2-2,i3-1,m) + scr(266)*u(i1+0,i2-2,i3-1,m) + scr(267)*u(i1+1,i2-2,i3-1,m)                                                                                                + scr(273)*u(i1-2,i2-1,i3-1,m) + scr(274)*u(i1-1,i2-1,i3-1,m) + scr(275)*u(i1+0,i2-1,i3-1,m) + scr(276)*u(i1+1,i2-1,i3-1,m) + scr(277)*u(i1+2,i2-1,i3-1,m)                                                                + scr(281)*u(i1-3,i2+0,i3-1,m) + scr(282)*u(i1-2,i2+0,i3-1,m) + scr(283)*u(i1-1,i2+0,i3-1,m) + scr(284)*u(i1+0,i2+0,i3-1,m) + scr(285)*u(i1+1,i2+0,i3-1,m) + scr(286)*u(i1+2,i2+0,i3-1,m) + scr(287)*u(i1+3,i2+0,i3-1,m)                                + scr(291)*u(i1-2,i2+1,i3-1,m) + scr(292)*u(i1-1,i2+1,i3-1,m) + scr(293)*u(i1+0,i2+1,i3-1,m) + scr(294)*u(i1+1,i2+1,i3-1,m) + scr(295)*u(i1+2,i2+1,i3-1,m)                                                                + scr(301)*u(i1-1,i2+2,i3-1,m) + scr(302)*u(i1+0,i2+2,i3-1,m) + scr(303)*u(i1+1,i2+2,i3-1,m)                                                                                                + scr(311)*u(i1+0,i2+3,i3-1,m)                                                                                                                                + scr(329)*u(i1+0,i2-4,i3+0,m)                                                                                                                                + scr(337)*u(i1-1,i2-3,i3+0,m) + scr(338)*u(i1+0,i2-3,i3+0,m) + scr(339)*u(i1+1,i2-3,i3+0,m)                                                                                                + scr(345)*u(i1-2,i2-2,i3+0,m) + scr(346)*u(i1-1,i2-2,i3+0,m) + scr(347)*u(i1+0,i2-2,i3+0,m) + scr(348)*u(i1+1,i2-2,i3+0,m) + scr(349)*u(i1+2,i2-2,i3+0,m)                                                                + scr(353)*u(i1-3,i2-1,i3+0,m) + scr(354)*u(i1-2,i2-1,i3+0,m) + scr(355)*u(i1-1,i2-1,i3+0,m) + scr(356)*u(i1+0,i2-1,i3+0,m) + scr(357)*u(i1+1,i2-1,i3+0,m) + scr(358)*u(i1+2,i2-1,i3+0,m) + scr(359)*u(i1+3,i2-1,i3+0,m)                                + scr(361)*u(i1-4,i2+0,i3+0,m) + scr(362)*u(i1-3,i2+0,i3+0,m) + scr(363)*u(i1-2,i2+0,i3+0,m) + scr(364)*u(i1-1,i2+0,i3+0,m) + scr(365)*u(i1+0,i2+0,i3+0,m) + scr(366)*u(i1+1,i2+0,i3+0,m) + scr(367)*u(i1+2,i2+0,i3+0,m) + scr(368)*u(i1+3,i2+0,i3+0,m) + scr(369)*u(i1+4,i2+0,i3+0,m)+ scr(371)*u(i1-3,i2+1,i3+0,m) + scr(372)*u(i1-2,i2+1,i3+0,m) + scr(373)*u(i1-1,i2+1,i3+0,m) + scr(374)*u(i1+0,i2+1,i3+0,m) + scr(375)*u(i1+1,i2+1,i3+0,m) + scr(376)*u(i1+2,i2+1,i3+0,m) + scr(377)*u(i1+3,i2+1,i3+0,m)                                + scr(381)*u(i1-2,i2+2,i3+0,m) + scr(382)*u(i1-1,i2+2,i3+0,m) + scr(383)*u(i1+0,i2+2,i3+0,m) + scr(384)*u(i1+1,i2+2,i3+0,m) + scr(385)*u(i1+2,i2+2,i3+0,m)                                                                + scr(391)*u(i1-1,i2+3,i3+0,m) + scr(392)*u(i1+0,i2+3,i3+0,m) + scr(393)*u(i1+1,i2+3,i3+0,m)                                                                                                + scr(401)*u(i1+0,i2+4,i3+0,m)                                                                                                                                + scr(419)*u(i1+0,i2-3,i3+1,m)                                                                                                                                + scr(427)*u(i1-1,i2-2,i3+1,m) + scr(428)*u(i1+0,i2-2,i3+1,m) + scr(429)*u(i1+1,i2-2,i3+1,m)                                                                                                + scr(435)*u(i1-2,i2-1,i3+1,m) + scr(436)*u(i1-1,i2-1,i3+1,m) + scr(437)*u(i1+0,i2-1,i3+1,m) + scr(438)*u(i1+1,i2-1,i3+1,m) + scr(439)*u(i1+2,i2-1,i3+1,m)                                                                + scr(443)*u(i1-3,i2+0,i3+1,m) + scr(444)*u(i1-2,i2+0,i3+1,m) + scr(445)*u(i1-1,i2+0,i3+1,m) + scr(446)*u(i1+0,i2+0,i3+1,m) + scr(447)*u(i1+1,i2+0,i3+1,m) + scr(448)*u(i1+2,i2+0,i3+1,m) + scr(449)*u(i1+3,i2+0,i3+1,m)                                + scr(453)*u(i1-2,i2+1,i3+1,m) + scr(454)*u(i1-1,i2+1,i3+1,m) + scr(455)*u(i1+0,i2+1,i3+1,m) + scr(456)*u(i1+1,i2+1,i3+1,m) + scr(457)*u(i1+2,i2+1,i3+1,m)                                                                + scr(463)*u(i1-1,i2+2,i3+1,m) + scr(464)*u(i1+0,i2+2,i3+1,m) + scr(465)*u(i1+1,i2+2,i3+1,m)                                                                                                + scr(473)*u(i1+0,i2+3,i3+1,m)                                                                                                                                + scr(509)*u(i1+0,i2-2,i3+2,m)                                                                                                                                + scr(517)*u(i1-1,i2-1,i3+2,m) + scr(518)*u(i1+0,i2-1,i3+2,m) + scr(519)*u(i1+1,i2-1,i3+2,m)                                                                                                + scr(525)*u(i1-2,i2+0,i3+2,m) + scr(526)*u(i1-1,i2+0,i3+2,m) + scr(527)*u(i1+0,i2+0,i3+2,m) + scr(528)*u(i1+1,i2+0,i3+2,m) + scr(529)*u(i1+2,i2+0,i3+2,m)                                                                + scr(535)*u(i1-1,i2+1,i3+2,m) + scr(536)*u(i1+0,i2+1,i3+2,m) + scr(537)*u(i1+1,i2+1,i3+2,m)                                                                                                + scr(545)*u(i1+0,i2+2,i3+2,m)                                                                                                                                + scr(599)*u(i1+0,i2-1,i3+3,m)                                                                                                                                + scr(607)*u(i1-1,i2+0,i3+3,m) + scr(608)*u(i1+0,i2+0,i3+3,m) + scr(609)*u(i1+1,i2+0,i3+3,m)                                                                                                + scr(617)*u(i1+0,i2+1,i3+3,m)                                                                                                                                + scr(689)*u(i1+0,i2+0,i3+4,m)                                                                                                                                +dtSq*fv(m)
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
