! This file automatically generated from advWaveStencil.bf90 with bpp.
    subroutine advWaveStencil3dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,sc,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
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
          real sc(1:343,nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
          real scr(1:343)
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
! Define variables to valuate stencil coefficients, dim=3, order=6, gridType=Rectangular
! File generated by cgWave/maple/writeStencilFiles.mpl
integer i1m3,i1m2,i1m1,i1p1,i1p2,i1p3
integer i2m3,i2m2,i2m1,i2p1,i2p2,i2p3
integer i3m3,i3m2,i3m1,i3p1,i3p2,i3p3
real t0,t1,t5,t9,t11,t12,t13,t16,t17,t18,t19,t22,t23,t25,t26,t27,t28,t29,t30,t31,t33,t34,t37,t39,t43,t44,t45,t46,t47,t50,t54,t55,t58,t61,t62,t63,t64,t66,t67,t68,t70,t71,t74,t79,t80,t83,t84,t85,t86,t88,t90,t91,t97,t98,t100,t101,t104,t107,t109,t110,t112,t113,t115,t116,t117,t118,t119,t121,t122,t124,t125,t130,t137,t143,t147,t151,t152,t153,t156,t157,t159,t162,t164,t168,t169,t172,t178,t179,t180,t181,t183,t189,t190,t194,t197,t198,t199,t201,t203,t204,t209,t216,t220,t224,t225,t227,t230,t232,t238,t245,t248,t249,t251,t253,t258,t274,t275,t276;
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
! Evaluate stencil coefficients, dim=3, order=6, gridType=Rectangular
! File generated by cgWave/maple/writeStencilFiles.mpl
i1m1=i1-1; i1p1=i1+1;
i2m1=i2-1; i2p1=i2+1;
i3m1=i3-1; i3p1=i3+1;
i1m2=i1-2; i1p2=i1+2;
i2m2=i2-2; i2p2=i2+2;
i3m2=i3-2; i3p2=i3+2;
t1 = cdz2i ** 2;
t5 = dt2 * cdz2i;
t9 = dt6 * t1 * cdz2i / 0.360e3 + t5 / 0.90e2 - dt4 * t1 / 0.72e2;
t11 = dt4 * cdy2i * cdz2i;
t12 = t11 / 0.72e2;
t13 = dt6 * t1;
t16 = -t12 + t13 * cdy2i / 0.120e3;
t17 = dt4 * cdx2i;
t18 = t17 * cdz2i;
t19 = t18 / 0.72e2;
t22 = -t19 + t13 * cdx2i / 0.120e3;
t23 = 0.3e1 / 0.20e2 * t5;
t25 = cdy2i * cdz2i;
t26 = 0.2e1 * t25;
t27 = cdx2i * cdz2i;
t28 = 0.2e1 * t27;
t29 = 0.4e1 * cdz2i;
t30 = 0.2e1 * cdy2i;
t31 = 0.2e1 * cdx2i;
t33 = cdz2i * (t29 + t30 + t31);
t34 = -0.2e1 * t1 - t26 - t28 - t33;
t37 = 0.2e1 * t1 * cdy2i;
t39 = 0.2e1 * t1 * cdx2i;
t43 = t25 / 0.6e1;
t44 = t27 / 0.6e1;
t45 = 0.2e1 * cdz2i;
t46 = cdy2i / 0.6e1;
t47 = cdx2i / 0.6e1;
t50 = -t43 - t44 - cdz2i * (t45 + t46 + t47);
t54 = cdy2i ** 2;
t55 = dt6 * t54;
t58 = -t12 + t55 * cdz2i / 0.120e3;
t61 = dt6 * cdx2i * t25 / 0.60e2;
t62 = 0.5e1 / 0.18e2 * t11;
t63 = 0.6e1 * t25;
t64 = -t63 - t28 - t33;
t66 = cdx2i * cdy2i;
t67 = 0.2e1 * t66;
t68 = 0.4e1 * cdy2i;
t70 = cdy2i * (t45 + t68 + t31);
t71 = -t63 - t67 - t70;
t74 = 0.4e1 * t66 * cdz2i;
t79 = cdx2i ** 2;
t80 = dt6 * t79;
t83 = -t19 + t80 * cdz2i / 0.120e3;
t84 = 0.5e1 / 0.18e2 * t18;
t85 = 0.6e1 * t27;
t86 = -t85 - t26 - t33;
t88 = 0.4e1 * cdx2i;
t90 = cdx2i * (t45 + t30 + t88);
t91 = -t85 - t67 - t90;
t97 = 0.19e2 / 0.6e1 * t25;
t98 = 0.19e2 / 0.6e1 * t27;
t100 = 0.19e2 / 0.6e1 * cdy2i;
t101 = 0.19e2 / 0.6e1 * cdx2i;
t104 = -t97 - t98 - cdz2i * (0.13e2 / 0.2e1 * cdz2i + t100 + t101);
t107 = 0.6e1 * cdz2i;
t109 = cdz2i * (t107 + t68 + t88);
t110 = 0.6e1 * cdy2i;
t112 = cdy2i * (t29 + t110 + t88);
t113 = 0.6e1 * cdx2i;
t115 = cdx2i * (t29 + t68 + t113);
t116 = 0.4e1 * t25;
t117 = 0.4e1 * t27;
t118 = 0.2e1 * t33;
t119 = -t1 - t109 - t112 - t115 - t116 - t117 - t118;
t121 = 0.8e1 * t25;
t122 = -t121 - t117 - t118;
t124 = 0.8e1 * t27;
t125 = -t124 - t116 - t118;
t130 = 0.3e1 / 0.2e1 * t5;
t137 = t84 - dt6 * (-cdx2i * t86 - cdz2i * t91 + t74) / 0.360e3;
t143 = t62 - dt6 * (-cdy2i * t64 - cdz2i * t71 + t74) / 0.360e3;
t147 = dt2 * cdy2i;
t151 = dt6 * t54 * cdy2i / 0.360e3 + t147 / 0.90e2 - dt4 * t54 / 0.72e2;
t152 = t17 * cdy2i;
t153 = t152 / 0.72e2;
t156 = -t153 + t55 * cdx2i / 0.120e3;
t157 = 0.3e1 / 0.20e2 * t147;
t159 = -0.2e1 * t54 - t26 - t67 - t70;
t162 = 0.2e1 * t54 * cdz2i;
t164 = 0.2e1 * t54 * cdx2i;
t168 = t66 / 0.6e1;
t169 = cdz2i / 0.6e1;
t172 = -t43 - t168 - cdy2i * (t169 + t30 + t47);
t178 = -t153 + t80 * cdy2i / 0.120e3;
t179 = 0.5e1 / 0.18e2 * t152;
t180 = 0.6e1 * t66;
t181 = -t180 - t26 - t70;
t183 = -t180 - t28 - t90;
t189 = 0.19e2 / 0.6e1 * t66;
t190 = 0.19e2 / 0.6e1 * cdz2i;
t194 = -t97 - t189 - cdy2i * (t190 + 0.13e2 / 0.2e1 * cdy2i + t101);
t197 = 0.4e1 * t66;
t198 = 0.2e1 * t70;
t199 = -t54 - t109 - t112 - t115 - t116 - t197 - t198;
t201 = -t121 - t197 - t198;
t203 = 0.8e1 * t66;
t204 = -t203 - t116 - t198;
t209 = 0.3e1 / 0.2e1 * t147;
t216 = t179 - dt6 * (-cdx2i * t181 - cdy2i * t183 + t74) / 0.360e3;
t220 = dt2 * cdx2i;
t224 = dt6 * t79 * cdx2i / 0.360e3 + t220 / 0.90e2 - dt4 * t79 / 0.72e2;
t225 = 0.3e1 / 0.20e2 * t220;
t227 = -0.2e1 * t79 - t28 - t67 - t90;
t230 = 0.2e1 * t79 * cdz2i;
t232 = 0.2e1 * t79 * cdy2i;
t238 = -t44 - t168 - cdx2i * (t169 + t46 + t31);
t245 = -t98 - t189 - cdx2i * (t190 + t100 + 0.13e2 / 0.2e1 * cdx2i);
t248 = 0.2e1 * t90;
t249 = -t79 - t109 - t112 - t115 - t117 - t197 - t248;
t251 = -t124 - t197 - t248;
t253 = -t203 - t117 - t248;
t258 = 0.3e1 / 0.2e1 * t220;
t274 = 0.2e1 * t109;
t275 = 0.2e1 * t112;
t276 = 0.2e1 * t115;
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
scr(25) = t9;
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
scr(41) = 0;
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
scr(67) = t16;
scr(68) = 0;
scr(69) = 0;
scr(70) = 0;
scr(71) = 0;
scr(72) = 0;
scr(73) = t22;
scr(74) = (-t23 + dt6 * (cdz2i * t34 - t37 - t39) / 0.360e3 - dt4 * t50 / 0.12e2);
scr(75) = t22;
scr(76) = 0;
scr(77) = 0;
scr(78) = 0;
scr(79) = 0;
scr(80) = 0;
scr(81) = t16;
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
scr(109) = t58;
scr(110) = 0;
scr(111) = 0;
scr(112) = 0;
scr(113) = 0;
scr(114) = 0;
scr(115) = t61;
scr(116) = (t62 + dt6 * (cdy2i * t64 + cdz2i * t71 - t74) / 0.360e3);
scr(117) = t61;
scr(118) = 0;
scr(119) = 0;
scr(120) = 0;
scr(121) = t83;
scr(122) = (t84 + dt6 * (cdx2i * t86 + cdz2i * t91 - t74) / 0.360e3);
scr(123) = (dt4 * t104 / 0.12e2 - dt6 * (cdx2i * t125 + cdy2i * t122 + cdz2i * t119) / 0.360e3 + t130);
scr(124) = t137;
scr(125) = t83;
scr(126) = 0;
scr(127) = 0;
scr(128) = 0;
scr(129) = t61;
scr(130) = t143;
scr(131) = t61;
scr(132) = 0;
scr(133) = 0;
scr(134) = 0;
scr(135) = 0;
scr(136) = 0;
scr(137) = t58;
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
scr(151) = t151;
scr(152) = 0;
scr(153) = 0;
scr(154) = 0;
scr(155) = 0;
scr(156) = 0;
scr(157) = t156;
scr(158) = (-t157 + dt6 * (cdy2i * t159 - t162 - t164) / 0.360e3 - dt4 * t172 / 0.12e2);
scr(159) = t156;
scr(160) = 0;
scr(161) = 0;
scr(162) = 0;
scr(163) = t178;
scr(164) = (t179 + dt6 * (cdx2i * t181 + cdy2i * t183 - t74) / 0.360e3);
scr(165) = (dt4 * t194 / 0.12e2 - dt6 * (cdx2i * t204 + cdy2i * t199 + cdz2i * t201) / 0.360e3 + t209);
scr(166) = t216;
scr(167) = t178;
scr(168) = 0;
scr(169) = t224;
scr(170) = (-t225 + dt6 * (cdx2i * t227 - t230 - t232) / 0.360e3 - dt4 * t238 / 0.12e2);
scr(171) = (dt4 * t245 / 0.12e2 - dt6 * (cdx2i * t249 + cdy2i * t253 + cdz2i * t251) / 0.360e3 + t258);
scr(172) = (dt4 * (cdz2i * (0.28e2 / 0.3e1 * cdz2i + t110 + t113) + cdy2i * (t107 + 0.28e2 / 0.3e1 * cdy2i + t113) + cdx2i * (t107 + t110 + 0.28e2 / 0.3e1 * cdx2i)) / 0.12e2 - 0.49e2 / 0.18e2 * dt2 * (cdz2i + cdy2i + cdx2i) - dt6 * (cdz2i * (t274 + t275 + t276 + t116 + t117 + t118) + cdy2i * (t274 + t275 + t276 + t198 + t116 + t197) + cdx2i * (t274 + t275 + t276 + t248 + t117 + t197)) / 0.360e3 + 0.2e1);
scr(173) = (t258 + dt6 * (-cdx2i * t249 - cdy2i * t253 - cdz2i * t251) / 0.360e3 + dt4 * t245 / 0.12e2);
scr(174) = (-t225 + dt6 * (cdx2i * t227 - t230 - t232) / 0.360e3 - dt4 * t238 / 0.12e2);
scr(175) = t224;
scr(176) = 0;
scr(177) = t178;
scr(178) = t216;
scr(179) = (t209 + dt6 * (-cdx2i * t204 - cdy2i * t199 - cdz2i * t201) / 0.360e3 + dt4 * t194 / 0.12e2);
scr(180) = t216;
scr(181) = t178;
scr(182) = 0;
scr(183) = 0;
scr(184) = 0;
scr(185) = t156;
scr(186) = (-t157 + dt6 * (cdy2i * t159 - t162 - t164) / 0.360e3 - dt4 * t172 / 0.12e2);
scr(187) = t156;
scr(188) = 0;
scr(189) = 0;
scr(190) = 0;
scr(191) = 0;
scr(192) = 0;
scr(193) = t151;
scr(194) = 0;
scr(195) = 0;
scr(196) = 0;
scr(197) = 0;
scr(198) = 0;
scr(199) = 0;
scr(200) = 0;
scr(201) = 0;
scr(202) = 0;
scr(203) = 0;
scr(204) = 0;
scr(205) = 0;
scr(206) = 0;
scr(207) = t58;
scr(208) = 0;
scr(209) = 0;
scr(210) = 0;
scr(211) = 0;
scr(212) = 0;
scr(213) = t61;
scr(214) = t143;
scr(215) = t61;
scr(216) = 0;
scr(217) = 0;
scr(218) = 0;
scr(219) = t83;
scr(220) = t137;
scr(221) = (t130 + dt6 * (-cdx2i * t125 - cdy2i * t122 - cdz2i * t119) / 0.360e3 + dt4 * t104 / 0.12e2);
scr(222) = t137;
scr(223) = t83;
scr(224) = 0;
scr(225) = 0;
scr(226) = 0;
scr(227) = t61;
scr(228) = t143;
scr(229) = t61;
scr(230) = 0;
scr(231) = 0;
scr(232) = 0;
scr(233) = 0;
scr(234) = 0;
scr(235) = t58;
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
scr(257) = 0;
scr(258) = 0;
scr(259) = 0;
scr(260) = 0;
scr(261) = 0;
scr(262) = 0;
scr(263) = t16;
scr(264) = 0;
scr(265) = 0;
scr(266) = 0;
scr(267) = 0;
scr(268) = 0;
scr(269) = t22;
scr(270) = (-t23 + dt6 * (cdz2i * t34 - t37 - t39) / 0.360e3 - dt4 * t50 / 0.12e2);
scr(271) = t22;
scr(272) = 0;
scr(273) = 0;
scr(274) = 0;
scr(275) = 0;
scr(276) = 0;
scr(277) = t16;
scr(278) = 0;
scr(279) = 0;
scr(280) = 0;
scr(281) = 0;
scr(282) = 0;
scr(283) = 0;
scr(284) = 0;
scr(285) = 0;
scr(286) = 0;
scr(287) = 0;
scr(288) = 0;
scr(289) = 0;
scr(290) = 0;
scr(291) = 0;
scr(292) = 0;
scr(293) = 0;
scr(294) = 0;
scr(295) = 0;
scr(296) = 0;
scr(297) = 0;
scr(298) = 0;
scr(299) = 0;
scr(300) = 0;
scr(301) = 0;
scr(302) = 0;
scr(303) = 0;
scr(304) = 0;
scr(305) = 0;
scr(306) = 0;
scr(307) = 0;
scr(308) = 0;
scr(309) = 0;
scr(310) = 0;
scr(311) = 0;
scr(312) = 0;
scr(313) = 0;
scr(314) = 0;
scr(315) = 0;
scr(316) = 0;
scr(317) = 0;
scr(318) = 0;
scr(319) = t9;
scr(320) = 0;
scr(321) = 0;
scr(322) = 0;
scr(323) = 0;
scr(324) = 0;
scr(325) = 0;
scr(326) = 0;
scr(327) = 0;
scr(328) = 0;
scr(329) = 0;
scr(330) = 0;
scr(331) = 0;
scr(332) = 0;
scr(333) = 0;
scr(334) = 0;
scr(335) = 0;
scr(336) = 0;
scr(337) = 0;
scr(338) = 0;
scr(339) = 0;
scr(340) = 0;
scr(341) = 0;
scr(342) = 0;
scr(343) = 0;
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
         !   updateWaveOpt(3,6,2,rectangular,NOFORCING)
         ! else 
         !   updateWaveOpt(3,6,2,rectangular,FORCING)
         ! end if
              else
                  if( addForcing.eq.0 )then
                          if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
                              write(*,'("advWaveStencil: ADVANCE dim=3 order=6 orderInTime=6, grid=rectangular... t=",e10.2)') t
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
! Stencil: nd=3, orderOfAccuracy=6, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ scr( 25)*u(i1+0,i2+0,i3-3,m)                                                                                                + scr( 67)*u(i1+0,i2-1,i3-2,m)                                                                                                + scr( 73)*u(i1-1,i2+0,i3-2,m) + scr( 74)*u(i1+0,i2+0,i3-2,m) + scr( 75)*u(i1+1,i2+0,i3-2,m)                                                                + scr( 81)*u(i1+0,i2+1,i3-2,m)                                                                                                + scr(109)*u(i1+0,i2-2,i3-1,m)                                                                                                + scr(115)*u(i1-1,i2-1,i3-1,m) + scr(116)*u(i1+0,i2-1,i3-1,m) + scr(117)*u(i1+1,i2-1,i3-1,m)                                                                + scr(121)*u(i1-2,i2+0,i3-1,m) + scr(122)*u(i1-1,i2+0,i3-1,m) + scr(123)*u(i1+0,i2+0,i3-1,m) + scr(124)*u(i1+1,i2+0,i3-1,m) + scr(125)*u(i1+2,i2+0,i3-1,m)                                + scr(129)*u(i1-1,i2+1,i3-1,m) + scr(130)*u(i1+0,i2+1,i3-1,m) + scr(131)*u(i1+1,i2+1,i3-1,m)                                                                + scr(137)*u(i1+0,i2+2,i3-1,m)                                                                                                + scr(151)*u(i1+0,i2-3,i3+0,m)                                                                                                + scr(157)*u(i1-1,i2-2,i3+0,m) + scr(158)*u(i1+0,i2-2,i3+0,m) + scr(159)*u(i1+1,i2-2,i3+0,m)                                                                + scr(163)*u(i1-2,i2-1,i3+0,m) + scr(164)*u(i1-1,i2-1,i3+0,m) + scr(165)*u(i1+0,i2-1,i3+0,m) + scr(166)*u(i1+1,i2-1,i3+0,m) + scr(167)*u(i1+2,i2-1,i3+0,m)                                + scr(169)*u(i1-3,i2+0,i3+0,m) + scr(170)*u(i1-2,i2+0,i3+0,m) + scr(171)*u(i1-1,i2+0,i3+0,m) + scr(172)*u(i1+0,i2+0,i3+0,m) + scr(173)*u(i1+1,i2+0,i3+0,m) + scr(174)*u(i1+2,i2+0,i3+0,m) + scr(175)*u(i1+3,i2+0,i3+0,m)+ scr(177)*u(i1-2,i2+1,i3+0,m) + scr(178)*u(i1-1,i2+1,i3+0,m) + scr(179)*u(i1+0,i2+1,i3+0,m) + scr(180)*u(i1+1,i2+1,i3+0,m) + scr(181)*u(i1+2,i2+1,i3+0,m)                                + scr(185)*u(i1-1,i2+2,i3+0,m) + scr(186)*u(i1+0,i2+2,i3+0,m) + scr(187)*u(i1+1,i2+2,i3+0,m)                                                                + scr(193)*u(i1+0,i2+3,i3+0,m)                                                                                                + scr(207)*u(i1+0,i2-2,i3+1,m)                                                                                                + scr(213)*u(i1-1,i2-1,i3+1,m) + scr(214)*u(i1+0,i2-1,i3+1,m) + scr(215)*u(i1+1,i2-1,i3+1,m)                                                                + scr(219)*u(i1-2,i2+0,i3+1,m) + scr(220)*u(i1-1,i2+0,i3+1,m) + scr(221)*u(i1+0,i2+0,i3+1,m) + scr(222)*u(i1+1,i2+0,i3+1,m) + scr(223)*u(i1+2,i2+0,i3+1,m)                                + scr(227)*u(i1-1,i2+1,i3+1,m) + scr(228)*u(i1+0,i2+1,i3+1,m) + scr(229)*u(i1+1,i2+1,i3+1,m)                                                                + scr(235)*u(i1+0,i2+2,i3+1,m)                                                                                                + scr(263)*u(i1+0,i2-1,i3+2,m)                                                                                                + scr(269)*u(i1-1,i2+0,i3+2,m) + scr(270)*u(i1+0,i2+0,i3+2,m) + scr(271)*u(i1+1,i2+0,i3+2,m)                                                                + scr(277)*u(i1+0,i2+1,i3+2,m)                                                                                                + scr(319)*u(i1+0,i2+0,i3+3,m)                                                                                                
                              end do
                              end do
                              end do
             ! endLoopsMask()
                  else
                          if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
                              write(*,'("advWaveStencil: ADVANCE dim=3 order=6 orderInTime=6, grid=rectangular... t=",e10.2)') t
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
! Stencil: nd=3, orderOfAccuracy=6, gridType=Rectangular
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ scr( 25)*u(i1+0,i2+0,i3-3,m)                                                                                                + scr( 67)*u(i1+0,i2-1,i3-2,m)                                                                                                + scr( 73)*u(i1-1,i2+0,i3-2,m) + scr( 74)*u(i1+0,i2+0,i3-2,m) + scr( 75)*u(i1+1,i2+0,i3-2,m)                                                                + scr( 81)*u(i1+0,i2+1,i3-2,m)                                                                                                + scr(109)*u(i1+0,i2-2,i3-1,m)                                                                                                + scr(115)*u(i1-1,i2-1,i3-1,m) + scr(116)*u(i1+0,i2-1,i3-1,m) + scr(117)*u(i1+1,i2-1,i3-1,m)                                                                + scr(121)*u(i1-2,i2+0,i3-1,m) + scr(122)*u(i1-1,i2+0,i3-1,m) + scr(123)*u(i1+0,i2+0,i3-1,m) + scr(124)*u(i1+1,i2+0,i3-1,m) + scr(125)*u(i1+2,i2+0,i3-1,m)                                + scr(129)*u(i1-1,i2+1,i3-1,m) + scr(130)*u(i1+0,i2+1,i3-1,m) + scr(131)*u(i1+1,i2+1,i3-1,m)                                                                + scr(137)*u(i1+0,i2+2,i3-1,m)                                                                                                + scr(151)*u(i1+0,i2-3,i3+0,m)                                                                                                + scr(157)*u(i1-1,i2-2,i3+0,m) + scr(158)*u(i1+0,i2-2,i3+0,m) + scr(159)*u(i1+1,i2-2,i3+0,m)                                                                + scr(163)*u(i1-2,i2-1,i3+0,m) + scr(164)*u(i1-1,i2-1,i3+0,m) + scr(165)*u(i1+0,i2-1,i3+0,m) + scr(166)*u(i1+1,i2-1,i3+0,m) + scr(167)*u(i1+2,i2-1,i3+0,m)                                + scr(169)*u(i1-3,i2+0,i3+0,m) + scr(170)*u(i1-2,i2+0,i3+0,m) + scr(171)*u(i1-1,i2+0,i3+0,m) + scr(172)*u(i1+0,i2+0,i3+0,m) + scr(173)*u(i1+1,i2+0,i3+0,m) + scr(174)*u(i1+2,i2+0,i3+0,m) + scr(175)*u(i1+3,i2+0,i3+0,m)+ scr(177)*u(i1-2,i2+1,i3+0,m) + scr(178)*u(i1-1,i2+1,i3+0,m) + scr(179)*u(i1+0,i2+1,i3+0,m) + scr(180)*u(i1+1,i2+1,i3+0,m) + scr(181)*u(i1+2,i2+1,i3+0,m)                                + scr(185)*u(i1-1,i2+2,i3+0,m) + scr(186)*u(i1+0,i2+2,i3+0,m) + scr(187)*u(i1+1,i2+2,i3+0,m)                                                                + scr(193)*u(i1+0,i2+3,i3+0,m)                                                                                                + scr(207)*u(i1+0,i2-2,i3+1,m)                                                                                                + scr(213)*u(i1-1,i2-1,i3+1,m) + scr(214)*u(i1+0,i2-1,i3+1,m) + scr(215)*u(i1+1,i2-1,i3+1,m)                                                                + scr(219)*u(i1-2,i2+0,i3+1,m) + scr(220)*u(i1-1,i2+0,i3+1,m) + scr(221)*u(i1+0,i2+0,i3+1,m) + scr(222)*u(i1+1,i2+0,i3+1,m) + scr(223)*u(i1+2,i2+0,i3+1,m)                                + scr(227)*u(i1-1,i2+1,i3+1,m) + scr(228)*u(i1+0,i2+1,i3+1,m) + scr(229)*u(i1+1,i2+1,i3+1,m)                                                                + scr(235)*u(i1+0,i2+2,i3+1,m)                                                                                                + scr(263)*u(i1+0,i2-1,i3+2,m)                                                                                                + scr(269)*u(i1-1,i2+0,i3+2,m) + scr(270)*u(i1+0,i2+0,i3+2,m) + scr(271)*u(i1+1,i2+0,i3+2,m)                                                                + scr(277)*u(i1+0,i2+1,i3+2,m)                                                                                                + scr(319)*u(i1+0,i2+0,i3+3,m)                                                                                                +dtSq*fv(m)
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
