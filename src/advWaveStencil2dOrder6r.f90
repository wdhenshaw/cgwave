! This file automatically generated from advWaveStencil.bf90 with bpp.
  subroutine advWaveStencil2dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,sc,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
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
     real sc(1:49,nd1a:nd1b,nd2a:nd2b)
     real scr(1:49)
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
! Define variables to valuate stencil coefficients, dim=2, order=6, gridType=Rectangular
! File generated by cgWave/maple/writeStencilFiles.mpl
integer i1m3,i1m2,i1m1,i1p1,i1p2,i1p3
integer i2m3,i2m2,i2m1,i2p1,i2p2,i2p3
real t0,t1,t4,t9,t14,t15,t16,t17,t18,t19,t21,t23,t24,t25,t28,t32,t33,t37,t42,t45,t46,t47,t49,t50,t51,t53,t59,t64,t67,t68,t69,t70,t72,t73,t75,t76,t78,t79,t90,t93,t98,t99,t101,t104,t111,t119,t122,t123,t124,t126,t143,t144;
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
! Evaluate stencil coefficients, dim=2, order=6, gridType=Rectangular
! File generated by cgWave/maple/writeStencilFiles.mpl
i1m1=i1-1; i1p1=i1+1;
i2m1=i2-1; i2p1=i2+1;
i1m2=i1-2; i1p2=i1+2;
i2m2=i2-2; i2p2=i2+2;
t1 = cdy2i ** 2;
t4 = dt2 * cdy2i;
t9 = -dt4 * t1 / 0.72e2 + t4 / 0.90e2 + dt6 * t1 * cdy2i / 0.360e3;
t14 = dt4 * cdx2i * cdy2i;
t15 = t14 / 0.72e2;
t16 = dt6 * cdx2i * t1 / 0.120e3 - t15;
t17 = 0.3e1 / 0.20e2 * t4;
t18 = 0.4e1 * cdy2i;
t19 = 0.2e1 * cdx2i;
t21 = cdy2i * (t18 + t19);
t23 = cdx2i * cdy2i;
t24 = 0.2e1 * t23;
t25 = -t21 - 0.2e1 * t1 - t24;
t28 = 0.2e1 * cdx2i * t1;
t32 = t23 / 0.6e1;
t33 = 0.2e1 * cdy2i;
t37 = -t32 - cdy2i * (t33 + cdx2i / 0.6e1);
t42 = cdx2i ** 2;
t45 = dt6 * cdy2i * t42 / 0.120e3 - t15;
t46 = 0.5e1 / 0.18e2 * t14;
t47 = 0.4e1 * cdx2i;
t49 = cdx2i * (t33 + t47);
t50 = 0.6e1 * t23;
t51 = -t49 - t50;
t53 = -t21 - t50;
t59 = 0.19e2 / 0.6e1 * t23;
t64 = -t59 - cdy2i * (0.13e2 / 0.2e1 * cdy2i + 0.19e2 / 0.6e1 * cdx2i);
t67 = 0.3e1 / 0.2e1 * t4;
t68 = 0.2e1 * t21;
t69 = 0.4e1 * t23;
t70 = 0.6e1 * cdy2i;
t72 = cdy2i * (t70 + t47);
t73 = 0.6e1 * cdx2i;
t75 = cdx2i * (t18 + t73);
t76 = -t68 - t1 - t69 - t72 - t75;
t78 = 0.8e1 * t23;
t79 = -t78 - t68;
t90 = t46 - dt6 * (-cdx2i * t53 - cdy2i * t51) / 0.360e3;
t93 = dt2 * cdx2i;
t98 = -dt4 * t42 / 0.72e2 + t93 / 0.90e2 + dt6 * t42 * cdx2i / 0.360e3;
t99 = 0.3e1 / 0.20e2 * t93;
t101 = -t49 - 0.2e1 * t42 - t24;
t104 = 0.2e1 * cdy2i * t42;
t111 = -t32 - cdx2i * (cdy2i / 0.6e1 + t19);
t119 = -t59 - cdx2i * (0.19e2 / 0.6e1 * cdy2i + 0.13e2 / 0.2e1 * cdx2i);
t122 = 0.3e1 / 0.2e1 * t93;
t123 = 0.2e1 * t49;
t124 = -t123 - t42 - t69 - t72 - t75;
t126 = -t78 - t123;
t143 = 0.2e1 * t72;
t144 = 0.2e1 * t75;
scr(1) = 0;
scr(2) = 0;
scr(3) = 0;
scr(4) = t9;
scr(5) = 0;
scr(6) = 0;
scr(7) = 0;
scr(8) = 0;
scr(9) = 0;
scr(10) = t16;
scr(11) = (-t17 + dt6 * (cdy2i * t25 - t28) / 0.360e3 - dt4 * t37 / 0.12e2);
scr(12) = t16;
scr(13) = 0;
scr(14) = 0;
scr(15) = 0;
scr(16) = t45;
scr(17) = (t46 + dt6 * (cdx2i * t53 + cdy2i * t51) / 0.360e3);
scr(18) = (dt4 * t64 / 0.12e2 + t67 - dt6 * (cdx2i * t79 + cdy2i * t76) / 0.360e3);
scr(19) = t90;
scr(20) = t45;
scr(21) = 0;
scr(22) = t98;
scr(23) = (-t99 + dt6 * (cdx2i * t101 - t104) / 0.360e3 - dt4 * t111 / 0.12e2);
scr(24) = (dt4 * t119 / 0.12e2 + t122 - dt6 * (cdx2i * t124 + cdy2i * t126) / 0.360e3);
scr(25) = (dt4 * (cdy2i * (0.28e2 / 0.3e1 * cdy2i + t73) + cdx2i * (t70 + 0.28e2 / 0.3e1 * cdx2i)) / 0.12e2 - 0.49e2 / 0.18e2 * dt2 * (cdy2i + cdx2i) + dt6 * (-cdy2i * (t68 + t69 + t143 + t144) - cdx2i * (t123 + t69 + t143 + t144)) / 0.360e3 + 0.2e1);
scr(26) = (t122 + dt4 * t119 / 0.12e2 + dt6 * (-cdx2i * t124 - cdy2i * t126) / 0.360e3);
scr(27) = (-dt4 * t111 / 0.12e2 - dt6 * (-cdx2i * t101 + t104) / 0.360e3 - t99);
scr(28) = t98;
scr(29) = 0;
scr(30) = t45;
scr(31) = t90;
scr(32) = (t67 + dt4 * t64 / 0.12e2 + dt6 * (-cdx2i * t79 - cdy2i * t76) / 0.360e3);
scr(33) = t90;
scr(34) = t45;
scr(35) = 0;
scr(36) = 0;
scr(37) = 0;
scr(38) = t16;
scr(39) = (-dt4 * t37 / 0.12e2 - dt6 * (-cdy2i * t25 + t28) / 0.360e3 - t17);
scr(40) = t16;
scr(41) = 0;
scr(42) = 0;
scr(43) = 0;
scr(44) = 0;
scr(45) = 0;
scr(46) = t9;
scr(47) = 0;
scr(48) = 0;
scr(49) = 0;
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
         !   updateWaveOpt(2,6,2,rectangular,NOFORCING)
         ! else 
         !   updateWaveOpt(2,6,2,rectangular,FORCING)
         ! end if
       else
         if( addForcing.eq.0 )then
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=2 order=6 orderInTime=6, grid=rectangular... t=",e10.2)') t
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
! Stencil: nd=2, orderOfAccuracy=6, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ scr(  4)*u(i1+0,i2-3,i3,m)                                                                                       + scr( 10)*u(i1-1,i2-2,i3,m) + scr( 11)*u(i1+0,i2-2,i3,m) + scr( 12)*u(i1+1,i2-2,i3,m)                                                          + scr( 16)*u(i1-2,i2-1,i3,m) + scr( 17)*u(i1-1,i2-1,i3,m) + scr( 18)*u(i1+0,i2-1,i3,m) + scr( 19)*u(i1+1,i2-1,i3,m) + scr( 20)*u(i1+2,i2-1,i3,m)                             + scr( 22)*u(i1-3,i2+0,i3,m) + scr( 23)*u(i1-2,i2+0,i3,m) + scr( 24)*u(i1-1,i2+0,i3,m) + scr( 25)*u(i1+0,i2+0,i3,m) + scr( 26)*u(i1+1,i2+0,i3,m) + scr( 27)*u(i1+2,i2+0,i3,m) + scr( 28)*u(i1+3,i2+0,i3,m)+ scr( 30)*u(i1-2,i2+1,i3,m) + scr( 31)*u(i1-1,i2+1,i3,m) + scr( 32)*u(i1+0,i2+1,i3,m) + scr( 33)*u(i1+1,i2+1,i3,m) + scr( 34)*u(i1+2,i2+1,i3,m)                             + scr( 38)*u(i1-1,i2+2,i3,m) + scr( 39)*u(i1+0,i2+2,i3,m) + scr( 40)*u(i1+1,i2+2,i3,m)                                                          + scr( 46)*u(i1+0,i2+3,i3,m)                                                                                       
               end do
               end do
               end do
             ! endLoopsMask()
         else
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=2 order=6 orderInTime=6, grid=rectangular... t=",e10.2)') t
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
! Stencil: nd=2, orderOfAccuracy=6, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ scr(  4)*u(i1+0,i2-3,i3,m)                                                                                       + scr( 10)*u(i1-1,i2-2,i3,m) + scr( 11)*u(i1+0,i2-2,i3,m) + scr( 12)*u(i1+1,i2-2,i3,m)                                                          + scr( 16)*u(i1-2,i2-1,i3,m) + scr( 17)*u(i1-1,i2-1,i3,m) + scr( 18)*u(i1+0,i2-1,i3,m) + scr( 19)*u(i1+1,i2-1,i3,m) + scr( 20)*u(i1+2,i2-1,i3,m)                             + scr( 22)*u(i1-3,i2+0,i3,m) + scr( 23)*u(i1-2,i2+0,i3,m) + scr( 24)*u(i1-1,i2+0,i3,m) + scr( 25)*u(i1+0,i2+0,i3,m) + scr( 26)*u(i1+1,i2+0,i3,m) + scr( 27)*u(i1+2,i2+0,i3,m) + scr( 28)*u(i1+3,i2+0,i3,m)+ scr( 30)*u(i1-2,i2+1,i3,m) + scr( 31)*u(i1-1,i2+1,i3,m) + scr( 32)*u(i1+0,i2+1,i3,m) + scr( 33)*u(i1+1,i2+1,i3,m) + scr( 34)*u(i1+2,i2+1,i3,m)                             + scr( 38)*u(i1-1,i2+2,i3,m) + scr( 39)*u(i1+0,i2+2,i3,m) + scr( 40)*u(i1+1,i2+2,i3,m)                                                          + scr( 46)*u(i1+0,i2+3,i3,m)                                                                                       +dtSq*fv(m)
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
