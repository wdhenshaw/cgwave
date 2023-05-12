! This file automatically generated from advWaveStencil.bf90 with bpp.
  subroutine advWaveStencil3dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,sc,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
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
! Define variables to valuate stencil coefficients, dim=3, order=4, gridType=Rectangular
! File generated by cgWave/maple/writeStencilFiles.mpl
integer i1m3,i1m2,i1m1,i1p1,i1p2,i1p3
integer i2m3,i2m2,i2m1,i2p1,i2p2,i2p3
integer i3m3,i3m2,i3m1,i3p1,i3p2,i3p3
real t0,t1,t3,t4,t7,t8,t10,t12,t14,t15,t16,t17,t20,t23,t25,t27,t28,t30,t32,t33,t34,t37,t40,t42,t44,t45,t46,t49,t52;
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
! Evaluate stencil coefficients, dim=3, order=4, gridType=Rectangular
! File generated by cgWave/maple/writeStencilFiles.mpl
i1m1=i1-1; i1p1=i1+1;
i2m1=i2-1; i2p1=i2+1;
i3m1=i3-1; i3p1=i3+1;
t1 = cdz2i ** 2;
t3 = dt2 * cdz2i;
t4 = dt4 * t1 - t3;
t7 = dt4 * cdy2i * cdz2i / 0.6e1;
t8 = dt4 * cdx2i;
t10 = t8 * cdz2i / 0.6e1;
t12 = 0.2e1 * cdy2i * cdz2i;
t14 = 0.2e1 * cdx2i * cdz2i;
t15 = 0.4e1 * cdz2i;
t16 = 0.2e1 * cdy2i;
t17 = 0.2e1 * cdx2i;
t20 = -t12 - t14 - cdz2i * (t15 + t16 + t17);
t23 = 0.4e1 / 0.3e1 * t3;
t25 = cdy2i ** 2;
t27 = dt2 * cdy2i;
t28 = dt4 * t25 - t27;
t30 = t8 * cdy2i / 0.6e1;
t32 = 0.2e1 * cdx2i * cdy2i;
t33 = 0.2e1 * cdz2i;
t34 = 0.4e1 * cdy2i;
t37 = -t12 - t32 - cdy2i * (t33 + t34 + t17);
t40 = 0.4e1 / 0.3e1 * t27;
t42 = cdx2i ** 2;
t44 = dt2 * cdx2i;
t45 = dt4 * t42 - t44;
t46 = 0.4e1 * cdx2i;
t49 = -t14 - t32 - cdx2i * (t33 + t16 + t46);
t52 = 0.4e1 / 0.3e1 * t44;
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
scr(13) = (t4 / 0.12e2);
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
scr(33) = t7;
scr(34) = 0;
scr(35) = 0;
scr(36) = 0;
scr(37) = t10;
scr(38) = (dt4 * t20 / 0.12e2 + t23);
scr(39) = t10;
scr(40) = 0;
scr(41) = 0;
scr(42) = 0;
scr(43) = t7;
scr(44) = 0;
scr(45) = 0;
scr(46) = 0;
scr(47) = 0;
scr(48) = 0;
scr(49) = 0;
scr(50) = 0;
scr(51) = 0;
scr(52) = 0;
scr(53) = (t28 / 0.12e2);
scr(54) = 0;
scr(55) = 0;
scr(56) = 0;
scr(57) = t30;
scr(58) = (dt4 * t37 / 0.12e2 + t40);
scr(59) = t30;
scr(60) = 0;
scr(61) = (t45 / 0.12e2);
scr(62) = (dt4 * t49 / 0.12e2 + t52);
scr(63) = (dt4 * (cdz2i * (0.6e1 * cdz2i + t34 + t46) + cdy2i * (t15 + 0.6e1 * cdy2i + t46) + cdx2i * (t15 + t34 + 0.6e1 * cdx2i)) / 0.12e2 - 0.5e1 / 0.2e1 * dt2 * (cdz2i + cdy2i + cdx2i) + 0.2e1);
scr(64) = (dt4 * t49 / 0.12e2 + t52);
scr(65) = (t45 / 0.12e2);
scr(66) = 0;
scr(67) = t30;
scr(68) = (dt4 * t37 / 0.12e2 + t40);
scr(69) = t30;
scr(70) = 0;
scr(71) = 0;
scr(72) = 0;
scr(73) = (t28 / 0.12e2);
scr(74) = 0;
scr(75) = 0;
scr(76) = 0;
scr(77) = 0;
scr(78) = 0;
scr(79) = 0;
scr(80) = 0;
scr(81) = 0;
scr(82) = 0;
scr(83) = t7;
scr(84) = 0;
scr(85) = 0;
scr(86) = 0;
scr(87) = t10;
scr(88) = (dt4 * t20 / 0.12e2 + t23);
scr(89) = t10;
scr(90) = 0;
scr(91) = 0;
scr(92) = 0;
scr(93) = t7;
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
scr(113) = (t4 / 0.12e2);
scr(114) = 0;
scr(115) = 0;
scr(116) = 0;
scr(117) = 0;
scr(118) = 0;
scr(119) = 0;
scr(120) = 0;
scr(121) = 0;
scr(122) = 0;
scr(123) = 0;
scr(124) = 0;
scr(125) = 0;
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
         !   updateWaveOpt(3,4,2,rectangular,NOFORCING)
         ! else 
         !   updateWaveOpt(3,4,2,rectangular,FORCING)
         ! end if
       else
         if( addForcing.eq.0 )then
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=3 order=4 orderInTime=4, grid=rectangular... t=",e10.2)') t
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
! Stencil: nd=3, orderOfAccuracy=4, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ scr( 13)*u(i1+0,i2+0,i3-2,m)                                                                + scr( 33)*u(i1+0,i2-1,i3-1,m)                                                                + scr( 37)*u(i1-1,i2+0,i3-1,m) + scr( 38)*u(i1+0,i2+0,i3-1,m) + scr( 39)*u(i1+1,i2+0,i3-1,m)                                + scr( 43)*u(i1+0,i2+1,i3-1,m)                                                                + scr( 53)*u(i1+0,i2-2,i3+0,m)                                                                + scr( 57)*u(i1-1,i2-1,i3+0,m) + scr( 58)*u(i1+0,i2-1,i3+0,m) + scr( 59)*u(i1+1,i2-1,i3+0,m)                                + scr( 61)*u(i1-2,i2+0,i3+0,m) + scr( 62)*u(i1-1,i2+0,i3+0,m) + scr( 63)*u(i1+0,i2+0,i3+0,m) + scr( 64)*u(i1+1,i2+0,i3+0,m) + scr( 65)*u(i1+2,i2+0,i3+0,m)+ scr( 67)*u(i1-1,i2+1,i3+0,m) + scr( 68)*u(i1+0,i2+1,i3+0,m) + scr( 69)*u(i1+1,i2+1,i3+0,m)                                + scr( 73)*u(i1+0,i2+2,i3+0,m)                                                                + scr( 83)*u(i1+0,i2-1,i3+1,m)                                                                + scr( 87)*u(i1-1,i2+0,i3+1,m) + scr( 88)*u(i1+0,i2+0,i3+1,m) + scr( 89)*u(i1+1,i2+0,i3+1,m)                                + scr( 93)*u(i1+0,i2+1,i3+1,m)                                                                + scr(113)*u(i1+0,i2+0,i3+2,m)                                                                
               end do
               end do
               end do
             ! endLoopsMask()
         else
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=3 order=4 orderInTime=4, grid=rectangular... t=",e10.2)') t
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
                             !fv(m) = fv(m) -( f(i1,i2,i3,freq) + cdtSqBy12*( cSq*(fxx23r(i1,i2,i3,freq) + fyy23r(i1,i2,i3,freq) + fzz23r(i1,i2,i3,freq)) - omega*omega*f(i1,i2,i3,freq)) )*coswt 
                     end do ! do freq  
                  else if( addForcing.ne.0 )then  
                     fv(m) = f(i1,i2,i3,0)
                  end if
! Stencil: nd=3, orderOfAccuracy=4, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ scr( 13)*u(i1+0,i2+0,i3-2,m)                                                                + scr( 33)*u(i1+0,i2-1,i3-1,m)                                                                + scr( 37)*u(i1-1,i2+0,i3-1,m) + scr( 38)*u(i1+0,i2+0,i3-1,m) + scr( 39)*u(i1+1,i2+0,i3-1,m)                                + scr( 43)*u(i1+0,i2+1,i3-1,m)                                                                + scr( 53)*u(i1+0,i2-2,i3+0,m)                                                                + scr( 57)*u(i1-1,i2-1,i3+0,m) + scr( 58)*u(i1+0,i2-1,i3+0,m) + scr( 59)*u(i1+1,i2-1,i3+0,m)                                + scr( 61)*u(i1-2,i2+0,i3+0,m) + scr( 62)*u(i1-1,i2+0,i3+0,m) + scr( 63)*u(i1+0,i2+0,i3+0,m) + scr( 64)*u(i1+1,i2+0,i3+0,m) + scr( 65)*u(i1+2,i2+0,i3+0,m)+ scr( 67)*u(i1-1,i2+1,i3+0,m) + scr( 68)*u(i1+0,i2+1,i3+0,m) + scr( 69)*u(i1+1,i2+1,i3+0,m)                                + scr( 73)*u(i1+0,i2+2,i3+0,m)                                                                + scr( 83)*u(i1+0,i2-1,i3+1,m)                                                                + scr( 87)*u(i1-1,i2+0,i3+1,m) + scr( 88)*u(i1+0,i2+0,i3+1,m) + scr( 89)*u(i1+1,i2+0,i3+1,m)                                + scr( 93)*u(i1+0,i2+1,i3+1,m)                                                                + scr(113)*u(i1+0,i2+0,i3+2,m)                                                                +dtSq*fv(m)
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
