! This file automatically generated from advWaveStencil.bf90 with bpp.
  subroutine advWaveStencil2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,sc,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
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
     real sc(1:81,nd1a:nd1b,nd2a:nd2b)
     real scr(1:81)    
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
! Define variables to valuate stencil coefficients, dim=2, order=8, gridType=Rectangular
! File generated by cgWave/maple/writeStencilFiles.mpl
integer i1m3,i1m2,i1m1,i1p1,i1p2,i1p3
integer i2m3,i2m2,i2m1,i2p1,i2p2,i2p3
real t0,t1,t3,t4,t9,t12,t14,t15,t22,t23,t24,t29,t32,t34,t35,t36,t37,t39,t41,t42,t43,t44,t45,t51,t52,t53,t54,t55,t58,t59,t65,t68,t69,t71,t74,t79,t80,t81,t82,t83,t85,t86,t87,t89,t90,t91,t92,t96,t97,t99,t100,t101,t102,t103,t104,t111,t112,t114,t115,t117,t118,t122,t123,t125,t126,t127,t129,t130,t131,t136,t141,t144,t145,t148,t149,t150,t151,t152,t153,t154,t155,t162,t163,t164,t170,t178,t182,t183,t185,t188,t189,t194,t195,t197,t198,t199,t201,t208,t209,t211,t215,t216,t221,t222,t223,t226,t227,t228,t229,t230,t231,t237,t239,t240,t243,t244,t246,t248,t254,t259,t260,t263,t264,t266,t268,t271,t274,t276,t277,t278,t284,t285,t288,t289,t290,t299,t303,t305,t306,t307,t308,t309,t315,t317,t324,t329,t332,t334,t335,t336,t342,t343,t344,t350,t352,t353,t355,t364,t367,t368,t371,t378,t379,t387,t392,t395,t397,t400,t417,t418,t419,t420,t428,t429,t431,t432,t465,t501,t514,t539,t540,t541,t549;
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
! Evaluate stencil coefficients, dim=2, order=8, gridType=Rectangular
! File generated by cgWave/maple/writeStencilFiles.mpl
i1m1=i1-1; i1p1=i1+1;
i2m1=i2-1; i2p1=i2+1;
i1m2=i1-2; i1p2=i1+2;
i2m2=i2-2; i2p2=i2+2;
i1m3=i1-3; i1p3=i1+3;
i2m3=i2-3; i2p3=i2+3;
t1 = dt2 * cdy2i;
t3 = cdy2i ** 2;
t4 = t3 ** 2;
t9 = t3 * cdy2i;
t12 = -t1 / 0.560e3 + dt8 * t4 / 0.20160e5 + 0.7e1 / 0.2880e4 * dt4 * t3 - dt6 * t9 / 0.1440e4;
t14 = dt4 * cdx2i * cdy2i;
t15 = t14 / 0.540e3;
t22 = t15 - dt6 * cdx2i * t3 / 0.720e3 + dt8 * t9 * cdx2i / 0.5040e4;
t23 = cdx2i * cdy2i;
t24 = t23 / 0.45e2;
t29 = -t24 - cdy2i * (0.2e1 / 0.5e1 * cdy2i + cdx2i / 0.45e2);
t32 = 0.8e1 / 0.315e3 * t1;
t34 = 0.2e1 * t9 * cdx2i;
t35 = 0.2e1 * t9;
t36 = 0.4e1 * cdy2i;
t37 = 0.2e1 * cdx2i;
t39 = cdy2i * (t36 + t37);
t41 = 0.2e1 * t23;
t42 = -t39 - 0.2e1 * t3 - t41;
t43 = cdy2i * t42;
t44 = cdx2i * t3;
t45 = 0.2e1 * t44;
t51 = t44 / 0.3e1;
t52 = 0.2e1 / 0.3e1 * t3;
t53 = t39 / 0.12e2;
t54 = t23 / 0.3e1;
t55 = 0.2e1 * cdy2i;
t58 = cdy2i * (t55 + cdx2i / 0.6e1);
t59 = -t52 - t53 - t54 - t58;
t65 = cdx2i ** 2;
t68 = dt8 * t65 * t3 / 0.3360e4;
t69 = t23 / 0.6e1;
t71 = t69 + t3 / 0.12e2;
t74 = t69 + t65 / 0.12e2;
t79 = t14 / 0.864e3;
t80 = t68 - dt6 * (cdx2i * t71 + cdy2i * t74) / 0.360e3 + t79;
t81 = 0.4e1 / 0.135e3 * t14;
t82 = 0.9e1 / 0.2e1 * t23;
t83 = 0.4e1 * cdx2i;
t85 = cdx2i * (t55 + t83);
t86 = t85 / 0.12e2;
t87 = -t82 - t86;
t89 = t23 / 0.2e1;
t90 = t3 / 0.3e1;
t91 = -t89 - t90 - t58;
t92 = cdx2i * t91;
t96 = 0.8e1 * t44;
t97 = t96 - t43;
t99 = 0.6e1 * t44;
t100 = 0.6e1 * t23;
t101 = -t85 - t100;
t102 = cdy2i * t101;
t103 = -t39 - t100;
t104 = cdx2i * t103;
t111 = t1 / 0.5e1;
t112 = 0.2e1 / 0.3e1 * t23;
t114 = 0.2e1 * t58;
t115 = -t112 - t3 / 0.2e1 - t114;
t117 = 0.25e2 / 0.6e1 * t23;
t118 = t39 / 0.3e1;
t122 = cdy2i * (0.13e2 / 0.2e1 * cdy2i + 0.19e2 / 0.6e1 * cdx2i);
t123 = 0.6e1 * cdy2i;
t125 = cdy2i * (t123 + t83);
t126 = t125 / 0.12e2;
t127 = 0.6e1 * cdx2i;
t129 = cdx2i * (t36 + t127);
t130 = t129 / 0.12e2;
t131 = -t117 - t52 - t114 - t118 - t122 - t126 - t130;
t136 = 0.41e2 / 0.120e3 * t23;
t141 = -t136 - cdy2i * (0.169e3 / 0.60e2 * cdy2i + 0.41e2 / 0.120e3 * cdx2i);
t144 = 0.10e2 * t44;
t145 = 0.2e1 * t43;
t148 = 0.2e1 * t39;
t149 = 0.4e1 * t23;
t150 = -t148 - t3 - t149 - t125 - t129;
t151 = cdy2i * t150;
t152 = 0.8e1 * t23;
t153 = -t152 - t148;
t154 = cdx2i * t153;
t155 = 0.4e1 * t44;
t162 = -cdy2i * t101;
t163 = -cdx2i * t103;
t164 = -t99 - t162 - t163;
t170 = -cdy2i * t87;
t178 = t65 * cdx2i;
t182 = t15 - dt6 * cdy2i * t65 / 0.720e3 + dt8 * t178 * cdy2i / 0.5040e4;
t183 = -t53 - t82;
t185 = t65 / 0.3e1;
t188 = cdx2i * (cdy2i / 0.6e1 + t37);
t189 = -t89 - t185 - t188;
t194 = cdy2i * t65;
t195 = 0.8e1 * t194;
t197 = -t85 - 0.2e1 * t65 - t41;
t198 = cdx2i * t197;
t199 = t195 - t198;
t201 = 0.6e1 * t194;
t208 = 0.71e2 / 0.6e1 * t23;
t209 = -t208 - t118 - t122;
t211 = t85 / 0.3e1;
t215 = cdx2i * (0.19e2 / 0.6e1 * cdy2i + 0.13e2 / 0.2e1 * cdx2i);
t216 = -t208 - t211 - t215;
t221 = 0.3e1 * t194;
t222 = 0.2e1 * t102;
t223 = 0.2e1 * t104;
t226 = 0.3e1 * t44;
t227 = 0.2e1 * t85;
t228 = -t227 - t65 - t149 - t125 - t129;
t229 = cdx2i * t228;
t230 = -t152 - t227;
t231 = cdy2i * t230;
t237 = 0.19e2 / 0.54e2 * t14;
t239 = 0.2e1 * t151;
t240 = 0.2e1 * t154;
t243 = 0.2e1 * t125;
t244 = 0.2e1 * t129;
t246 = cdy2i * (t148 + t149 + t243 + t244);
t248 = cdx2i * (t227 + t149 + t243 + t244);
t254 = 0.8e1 / 0.5e1 * t1;
t259 = 0.35e2 / 0.9e1 * t23;
t260 = -cdy2i * (0.122e3 / 0.15e2 * cdy2i + 0.35e2 / 0.9e1 * cdx2i) - t259;
t263 = 0.46e2 / 0.3e1 * t23;
t264 = 0.2e1 * t122;
t266 = -t263 - t264 - t39 / 0.2e1;
t268 = 0.23e2 / 0.3e1 * t23;
t271 = cdy2i * (0.28e2 / 0.3e1 * cdy2i + t127);
t274 = cdx2i * (t123 + 0.28e2 / 0.3e1 * cdx2i);
t276 = t125 / 0.3e1;
t277 = t129 / 0.3e1;
t278 = -t268 - t271 - t58 - t274 - t264 - 0.7e1 / 0.12e2 * t39 - t276 - t277 - t90;
t284 = 0.2e1 * t162;
t285 = 0.2e1 * t163;
t288 = -cdx2i * t228;
t289 = -cdy2i * t230;
t290 = -t226 - t284 - t288 - t289 - t285;
t299 = dt6 * (-cdx2i * t209 - cdy2i * t216) / 0.360e3;
t303 = -cdx2i * t183 - cdy2i * t189;
t305 = dt6 * t303 / 0.360e3;
t306 = t201 + t162 + t163;
t307 = cdx2i * t306;
t308 = -cdx2i * t197;
t309 = t308 + t195;
t315 = dt2 * cdx2i;
t317 = t65 ** 2;
t324 = -t315 / 0.560e3 + dt8 * t317 / 0.20160e5 + 0.7e1 / 0.2880e4 * dt4 * t65 - dt6 * t178 / 0.1440e4;
t329 = -t24 - cdx2i * (cdy2i / 0.45e2 + 0.2e1 / 0.5e1 * cdx2i);
t332 = 0.8e1 / 0.315e3 * t315;
t334 = 0.2e1 * t178 * cdy2i;
t335 = 0.2e1 * t178;
t336 = 0.2e1 * t194;
t342 = t194 / 0.3e1;
t343 = 0.2e1 / 0.3e1 * t65;
t344 = -t54 - t86 - t343 - t188;
t350 = t315 / 0.5e1;
t352 = 0.2e1 * t188;
t353 = -t112 - t65 / 0.2e1 - t352;
t355 = -t117 - t211 - t343 - t215 - t126 - t130 - t352;
t364 = -t136 - cdx2i * (0.41e2 / 0.120e3 * cdy2i + 0.169e3 / 0.60e2 * cdx2i);
t367 = 0.10e2 * t194;
t368 = 0.2e1 * t198;
t371 = 0.4e1 * t194;
t378 = 0.2e1 * t229;
t379 = 0.2e1 * t231;
t387 = 0.8e1 / 0.5e1 * t315;
t392 = -cdx2i * (0.35e2 / 0.9e1 * cdy2i + 0.122e3 / 0.15e2 * cdx2i) - t259;
t395 = 0.2e1 * t215;
t397 = -t263 - t395 - t85 / 0.2e1;
t400 = -t268 - t395 - t271 - t274 - t188 - t276 - t277 - 0.7e1 / 0.12e2 * t85 - t185;
t417 = -cdy2i * t150;
t418 = -cdx2i * t153;
t419 = 0.2e1 * t246;
t420 = 0.2e1 * t248;
t428 = 0.2e1 * t271;
t429 = 0.2e1 * t274;
t431 = t125 / 0.2e1;
t432 = t129 / 0.2e1;
t465 = 0.2e1 * t308;
t501 = -t221 - t285 - t417 - t418 - t284;
t514 = -cdy2i * t42;
t539 = dt6 * (-cdx2i * t91 + t170) / 0.360e3;
t540 = -cdy2i * t164;
t541 = t514 + t96;
t549 = 0.2e1 * t514;
scr(1) = 0;
scr(2) = 0;
scr(3) = 0;
scr(4) = 0;
scr(5) = t12;
scr(6) = 0;
scr(7) = 0;
scr(8) = 0;
scr(9) = 0;
scr(10) = 0;
scr(11) = 0;
scr(12) = 0;
scr(13) = t22;
scr(14) = (dt4 * t29 / 0.12e2 + t32 - dt8 * (t34 + cdy2i * (t35 - t43 + t45)) / 0.20160e5 - dt6 * (cdy2i * t59 - t51) / 0.360e3);
scr(15) = t22;
scr(16) = 0;
scr(17) = 0;
scr(18) = 0;
scr(19) = 0;
scr(20) = 0;
scr(21) = t80;
scr(22) = (-t81 - dt6 * (cdy2i * t87 + t92) / 0.360e3 - dt8 * (cdx2i * t97 + cdy2i * (t99 - t102 - t104)) / 0.20160e5);
scr(23) = (-t111 + dt6 * (cdx2i * t115 + cdy2i * t131) / 0.360e3 - dt4 * t141 / 0.12e2 + dt8 * (cdx2i * (t144 - t145) - cdy2i * (-t9 + t151 + t154 + t145 - t155)) / 0.20160e5);
scr(24) = (dt8 * (-cdx2i * t97 + cdy2i * t164) / 0.20160e5 - dt6 * (-t170 + t92) / 0.360e3 - t81);
scr(25) = t80;
scr(26) = 0;
scr(27) = 0;
scr(28) = 0;
scr(29) = t182;
scr(30) = (-t81 - dt6 * (cdx2i * t183 + cdy2i * t189) / 0.360e3 - dt8 * (cdy2i * t199 + cdx2i * (t201 - t102 - t104)) / 0.20160e5);
scr(31) = (dt6 * (cdx2i * t209 + cdy2i * t216) / 0.360e3 - dt8 * (cdx2i * (-t221 + t151 + t222 + t223 + t154) + cdy2i * (-t226 + t229 + t222 + t223 + t231)) / 0.20160e5 + t237);
scr(32) = (dt8 * (cdx2i * (t239 - t162 + t102 + t104 + t240 - t163) + cdy2i * (t239 - t246 - t248 + t240 + t43 - t45)) / 0.20160e5 + t254 + dt4 * t260 / 0.12e2 - dt6 * (cdx2i * t266 + cdy2i * t278) / 0.360e3);
scr(33) = (t237 - dt8 * (-cdx2i * (t221 - t151 + t284 - t154 + t285) + cdy2i * t290) / 0.20160e5 - t299);
scr(34) = (t305 - dt8 * (cdy2i * t309 + t307) / 0.20160e5 - t81);
scr(35) = t182;
scr(36) = 0;
scr(37) = t324;
scr(38) = (dt4 * t329 / 0.12e2 + t332 - dt8 * (t334 + cdx2i * (t335 - t198 + t336)) / 0.20160e5 - dt6 * (cdx2i * t344 - t342) / 0.360e3);
scr(39) = (-t350 + dt6 * (cdx2i * t355 + cdy2i * t353) / 0.360e3 - dt4 * t364 / 0.12e2 + dt8 * (cdy2i * (t367 - t368) - cdx2i * (-t178 + t229 + t231 + t368 - t371)) / 0.20160e5);
scr(40) = (dt8 * (cdy2i * (t378 - t163 + t102 + t104 + t379 - t162) + cdx2i * (t378 - t246 - t248 + t379 + t198 - t336)) / 0.20160e5 + t387 + dt4 * t392 / 0.12e2 - dt6 * (cdx2i * t400 + cdy2i * t397) / 0.360e3);
scr(41) = (dt4 * (cdy2i * (0.91e2 / 0.8e1 * cdy2i + 0.257e3 / 0.36e2 * cdx2i) + cdx2i * (0.257e3 / 0.36e2 * cdy2i + 0.91e2 / 0.8e1 * cdx2i)) / 0.12e2 - dt8 * (cdy2i * (t151 - t417 - t418 - t419 - t420 + t154) + cdx2i * (t229 - t288 - t289 - t419 - t420 + t231)) / 0.20160e5 - dt6 * (cdy2i * (t268 + t264 + t428 + t429 + t3 / 0.6e1 + t431 + t432 + 0.2e1 / 0.3e1 * t39) + cdx2i * (t268 + t395 + t428 + t429 + t65 / 0.6e1 + 0.2e1 / 0.3e1 * t85 + t431 + t432)) / 0.360e3 - 0.205e3 / 0.72e2 * dt2 * (cdy2i + cdx2i) + 0.2e1);
scr(42) = (t387 + dt6 * (-cdx2i * t400 - cdy2i * t397) / 0.360e3 + dt4 * t392 / 0.12e2 - dt8 * (cdx2i * (t308 + 0.2e1 * t288 + 0.2e1 * t289 + t246 + t248 + t336) + 0.2e1 * cdy2i * (t162 + t163 + t288 + t289)) / 0.20160e5);
scr(43) = (-dt4 * t364 / 0.12e2 - t350 + dt8 * (cdy2i * (t465 + t367) + cdx2i * (t465 + t288 + t289 + t371 + t178)) / 0.20160e5 - dt6 * (-cdx2i * t355 - cdy2i * t353) / 0.360e3);
scr(44) = (t332 - dt8 * (t334 + cdx2i * (t308 + t336 + t335)) / 0.20160e5 + dt4 * t329 / 0.12e2 + dt6 * (-cdx2i * t344 + t342) / 0.360e3);
scr(45) = t324;
scr(46) = 0;
scr(47) = t182;
scr(48) = (dt8 * (-cdx2i * t306 - cdy2i * t199) / 0.20160e5 + dt6 * t303 / 0.360e3 - t81);
scr(49) = (t237 - dt8 * (-cdy2i * (t226 - t229 + t285 - t231 + t284) + cdx2i * t501) / 0.20160e5 - t299);
scr(50) = (t254 + dt6 * (-cdx2i * t266 - cdy2i * t278) / 0.360e3 + dt4 * t260 / 0.12e2 - dt8 * (cdy2i * (t514 + 0.2e1 * t417 + 0.2e1 * t418 + t246 + t248 + t45) + 0.2e1 * cdx2i * (t162 + t163 + t417 + t418)) / 0.20160e5);
scr(51) = (dt8 * (-cdx2i * t501 - cdy2i * t290) / 0.20160e5 - t299 + t237);
scr(52) = (t305 - dt8 * (cdy2i * t309 + t307) / 0.20160e5 - t81);
scr(53) = t182;
scr(54) = 0;
scr(55) = 0;
scr(56) = 0;
scr(57) = t80;
scr(58) = (t539 - dt8 * (cdx2i * t541 + t540) / 0.20160e5 - t81);
scr(59) = (-dt4 * t141 / 0.12e2 - t111 + dt8 * (cdx2i * (t549 + t144) + cdy2i * (t549 + t417 + t418 + t155 + t9)) / 0.20160e5 - dt6 * (-cdx2i * t115 - cdy2i * t131) / 0.360e3);
scr(60) = (t539 - dt8 * (cdx2i * t541 + t540) / 0.20160e5 - t81);
scr(61) = (dt6 * (-cdx2i * t71 - cdy2i * t74) / 0.360e3 + t79 + t68);
scr(62) = 0;
scr(63) = 0;
scr(64) = 0;
scr(65) = 0;
scr(66) = 0;
scr(67) = t22;
scr(68) = (t32 - dt8 * (t34 + cdy2i * (t514 + t45 + t35)) / 0.20160e5 + dt4 * t29 / 0.12e2 + dt6 * (-cdy2i * t59 + t51) / 0.360e3);
scr(69) = t22;
scr(70) = 0;
scr(71) = 0;
scr(72) = 0;
scr(73) = 0;
scr(74) = 0;
scr(75) = 0;
scr(76) = 0;
scr(77) = t12;
scr(78) = 0;
scr(79) = 0;
scr(80) = 0;
scr(81) = 0;
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
         !   updateWaveOpt(2,8,2,rectangular,NOFORCING)
         ! else 
         !   updateWaveOpt(2,8,2,rectangular,FORCING)
         ! end if
       else
         if( addForcing.eq.0 )then
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=2 order=8 orderInTime=8, grid=rectangular... t=",e10.2)') t
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
! Stencil: nd=2, orderOfAccuracy=8, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ scr(  5)*u(i1+0,i2-4,i3,m)                                                                                                                    + scr( 13)*u(i1-1,i2-3,i3,m) + scr( 14)*u(i1+0,i2-3,i3,m) + scr( 15)*u(i1+1,i2-3,i3,m)                                                                                       + scr( 21)*u(i1-2,i2-2,i3,m) + scr( 22)*u(i1-1,i2-2,i3,m) + scr( 23)*u(i1+0,i2-2,i3,m) + scr( 24)*u(i1+1,i2-2,i3,m) + scr( 25)*u(i1+2,i2-2,i3,m)                                                          + scr( 29)*u(i1-3,i2-1,i3,m) + scr( 30)*u(i1-2,i2-1,i3,m) + scr( 31)*u(i1-1,i2-1,i3,m) + scr( 32)*u(i1+0,i2-1,i3,m) + scr( 33)*u(i1+1,i2-1,i3,m) + scr( 34)*u(i1+2,i2-1,i3,m) + scr( 35)*u(i1+3,i2-1,i3,m)                             + scr( 37)*u(i1-4,i2+0,i3,m) + scr( 38)*u(i1-3,i2+0,i3,m) + scr( 39)*u(i1-2,i2+0,i3,m) + scr( 40)*u(i1-1,i2+0,i3,m) + scr( 41)*u(i1+0,i2+0,i3,m) + scr( 42)*u(i1+1,i2+0,i3,m) + scr( 43)*u(i1+2,i2+0,i3,m) + scr( 44)*u(i1+3,i2+0,i3,m) + scr( 45)*u(i1+4,i2+0,i3,m)+ scr( 47)*u(i1-3,i2+1,i3,m) + scr( 48)*u(i1-2,i2+1,i3,m) + scr( 49)*u(i1-1,i2+1,i3,m) + scr( 50)*u(i1+0,i2+1,i3,m) + scr( 51)*u(i1+1,i2+1,i3,m) + scr( 52)*u(i1+2,i2+1,i3,m) + scr( 53)*u(i1+3,i2+1,i3,m)                             + scr( 57)*u(i1-2,i2+2,i3,m) + scr( 58)*u(i1-1,i2+2,i3,m) + scr( 59)*u(i1+0,i2+2,i3,m) + scr( 60)*u(i1+1,i2+2,i3,m) + scr( 61)*u(i1+2,i2+2,i3,m)                                                          + scr( 67)*u(i1-1,i2+3,i3,m) + scr( 68)*u(i1+0,i2+3,i3,m) + scr( 69)*u(i1+1,i2+3,i3,m)                                                                                       + scr( 77)*u(i1+0,i2+4,i3,m)                                                                                                                    
               end do
               end do
               end do
             ! endLoopsMask()
         else
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=2 order=8 orderInTime=8, grid=rectangular... t=",e10.2)') t
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
! Stencil: nd=2, orderOfAccuracy=8, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ scr(  5)*u(i1+0,i2-4,i3,m)                                                                                                                    + scr( 13)*u(i1-1,i2-3,i3,m) + scr( 14)*u(i1+0,i2-3,i3,m) + scr( 15)*u(i1+1,i2-3,i3,m)                                                                                       + scr( 21)*u(i1-2,i2-2,i3,m) + scr( 22)*u(i1-1,i2-2,i3,m) + scr( 23)*u(i1+0,i2-2,i3,m) + scr( 24)*u(i1+1,i2-2,i3,m) + scr( 25)*u(i1+2,i2-2,i3,m)                                                          + scr( 29)*u(i1-3,i2-1,i3,m) + scr( 30)*u(i1-2,i2-1,i3,m) + scr( 31)*u(i1-1,i2-1,i3,m) + scr( 32)*u(i1+0,i2-1,i3,m) + scr( 33)*u(i1+1,i2-1,i3,m) + scr( 34)*u(i1+2,i2-1,i3,m) + scr( 35)*u(i1+3,i2-1,i3,m)                             + scr( 37)*u(i1-4,i2+0,i3,m) + scr( 38)*u(i1-3,i2+0,i3,m) + scr( 39)*u(i1-2,i2+0,i3,m) + scr( 40)*u(i1-1,i2+0,i3,m) + scr( 41)*u(i1+0,i2+0,i3,m) + scr( 42)*u(i1+1,i2+0,i3,m) + scr( 43)*u(i1+2,i2+0,i3,m) + scr( 44)*u(i1+3,i2+0,i3,m) + scr( 45)*u(i1+4,i2+0,i3,m)+ scr( 47)*u(i1-3,i2+1,i3,m) + scr( 48)*u(i1-2,i2+1,i3,m) + scr( 49)*u(i1-1,i2+1,i3,m) + scr( 50)*u(i1+0,i2+1,i3,m) + scr( 51)*u(i1+1,i2+1,i3,m) + scr( 52)*u(i1+2,i2+1,i3,m) + scr( 53)*u(i1+3,i2+1,i3,m)                             + scr( 57)*u(i1-2,i2+2,i3,m) + scr( 58)*u(i1-1,i2+2,i3,m) + scr( 59)*u(i1+0,i2+2,i3,m) + scr( 60)*u(i1+1,i2+2,i3,m) + scr( 61)*u(i1+2,i2+2,i3,m)                                                          + scr( 67)*u(i1-1,i2+3,i3,m) + scr( 68)*u(i1+0,i2+3,i3,m) + scr( 69)*u(i1+1,i2+3,i3,m)                                                                                       + scr( 77)*u(i1+0,i2+4,i3,m)                                                                                                                    +dtSq*fv(m)
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
