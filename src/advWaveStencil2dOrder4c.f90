! This file automatically generated from advWaveStencil.bf90 with bpp.
  subroutine advWaveStencil2dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,sc,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
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
     real sc(1:25,nd1a:nd1b,nd2a:nd2b)
     real scr(1:25)
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
! Define variables to valuate stencil coefficients, dim=2, order=4, gridType=Curvilinear
! File generated by cgWave/maple/writeStencilFiles.mpl
integer i1m3,i1m2,i1m1,i1p1,i1p2,i1p3
integer i2m3,i2m2,i2m1,i2p1,i2p2,i2p3
real t0,t1,t2,t6,t7,t9,t10,t12,t13,t14,t15,t16,t22,t23,t27,t28,t29,t32,t33,t36,t43,t44,t45,t46,t56,t57,t59,t60,t62,t63,t64,t65,t66,t73,t74,t75,t76,t77,t80,t81,t82,t85,t86,t89,t90,t99,t100,t101,t102,t106,t110,t112,t113,t114,t122,t123,t127,t128,t129,t134,t135,t147,t149,t158,t161,t162,t163,t164,t167,t174,t175,t176,t177,t182,t185,t187,t188,t189,t197,t198,t202,t205,t206,t207,t208,t211,t212,t213,t214,t217,t218,t219,t222,t233,t234,t235,t236,t263,t264,t265,t266,t275,t276,t277,t282,t283,t294,t295,t296,t297,t336,t338;
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
          write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF 2D'
          numGhost1=orderOfAccuracy/2 -1; ! check me 
          n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
          n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
          n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            rx = rsxy(i1,i2,i3,0,0)
            ry = rsxy(i1,i2,i3,0,1)
            sx = rsxy(i1,i2,i3,1,0)
            sy = rsxy(i1,i2,i3,1,1)
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
            if( diffOrder1.eq.2 )then
              rxr = (rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0))*(.5*dr1i) 
              ryr = (rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1))*(.5*dr1i) 
              sxr = (rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0))*(.5*dr1i) 
              syr = (rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1))*(.5*dr1i) 
            elseif( diffOrder1.eq.4 )then
              rxr = ( 8*(rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0)) -(rsxy(i1+2,i2,i3,0,0)-rsxy(i1-2,i2,i3,0,0)) )*(dr1i/12.) 
              ryr = ( 8*(rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1)) -(rsxy(i1+2,i2,i3,0,1)-rsxy(i1-2,i2,i3,0,1)) )*(dr1i/12.) 
              sxr = ( 8*(rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0)) -(rsxy(i1+2,i2,i3,1,0)-rsxy(i1-2,i2,i3,1,0)) )*(dr1i/12.) 
              syr = ( 8*(rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1)) -(rsxy(i1+2,i2,i3,1,1)-rsxy(i1-2,i2,i3,1,1)) )*(dr1i/12.) 
            elseif( diffOrder1.eq.6 )then
              rxr = ( 45.*(rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0)) -9.*(rsxy(i1+2,i2,i3,0,0)-rsxy(i1-2,i2,i3,0,0)) +(rsxy(i1+3,i2,i3,0,0)-rsxy(i1-3,i2,i3,0,0)) )*(dr1i/60.) 
              ryr = ( 45.*(rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1)) -9.*(rsxy(i1+2,i2,i3,0,1)-rsxy(i1-2,i2,i3,0,1)) +(rsxy(i1+3,i2,i3,0,1)-rsxy(i1-3,i2,i3,0,1)) )*(dr1i/60.) 
              sxr = ( 45.*(rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0)) -9.*(rsxy(i1+2,i2,i3,1,0)-rsxy(i1-2,i2,i3,1,0)) +(rsxy(i1+3,i2,i3,1,0)-rsxy(i1-3,i2,i3,1,0)) )*(dr1i/60.) 
              syr = ( 45.*(rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1)) -9.*(rsxy(i1+2,i2,i3,1,1)-rsxy(i1-2,i2,i3,1,1)) +(rsxy(i1+3,i2,i3,1,1)-rsxy(i1-3,i2,i3,1,1)) )*(dr1i/60.) 
            elseif( diffOrder1.eq.8 )then
              rxr = ( 672.*(rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0)) -168.*(rsxy(i1+2,i2,i3,0,0)-rsxy(i1-2,i2,i3,0,0)) +32*(rsxy(i1+3,i2,i3,0,0)-rsxy(i1-3,i2,i3,0,0)) -3.*(rsxy(i1+4,i2,i3,0,0)-rsxy(i1-4,i2,i3,0,0)) )*(dr1i/840.) 
              ryr = ( 672.*(rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1)) -168.*(rsxy(i1+2,i2,i3,0,1)-rsxy(i1-2,i2,i3,0,1)) +32*(rsxy(i1+3,i2,i3,0,1)-rsxy(i1-3,i2,i3,0,1)) -3.*(rsxy(i1+4,i2,i3,0,1)-rsxy(i1-4,i2,i3,0,1)) )*(dr1i/840.) 
              sxr = ( 672.*(rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0)) -168.*(rsxy(i1+2,i2,i3,1,0)-rsxy(i1-2,i2,i3,1,0)) +32*(rsxy(i1+3,i2,i3,1,0)-rsxy(i1-3,i2,i3,1,0)) -3.*(rsxy(i1+4,i2,i3,1,0)-rsxy(i1-4,i2,i3,1,0)) )*(dr1i/840.) 
              syr = ( 672.*(rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1)) -168.*(rsxy(i1+2,i2,i3,1,1)-rsxy(i1-2,i2,i3,1,1)) +32*(rsxy(i1+3,i2,i3,1,1)-rsxy(i1-3,i2,i3,1,1)) -3.*(rsxy(i1+4,i2,i3,1,1)-rsxy(i1-4,i2,i3,1,1)) )*(dr1i/840.) 
            end if
            if( diffOrder2.eq.2 )then
              rxs = (rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0))*(.5*dr2i) 
              rys = (rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1))*(.5*dr2i) 
              sxs = (rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0))*(.5*dr2i) 
              sys = (rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1))*(.5*dr2i) 
            elseif( diffOrder2.eq.4 )then
              rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
              rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
              sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
              sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
            elseif( diffOrder2.eq.6 )then
              rxs = ( 45.*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -9.*(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) +(rsxy(i1,i2+3,i3,0,0)-rsxy(i1,i2-3,i3,0,0)) )*(dr2i/60.) 
              rys = ( 45.*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -9.*(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) +(rsxy(i1,i2+3,i3,0,1)-rsxy(i1,i2-3,i3,0,1)) )*(dr2i/60.) 
              sxs = ( 45.*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -9.*(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) +(rsxy(i1,i2+3,i3,1,0)-rsxy(i1,i2-3,i3,1,0)) )*(dr2i/60.) 
              sys = ( 45.*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -9.*(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) +(rsxy(i1,i2+3,i3,1,1)-rsxy(i1,i2-3,i3,1,1)) )*(dr2i/60.) 
            elseif( diffOrder2.eq.8 )then
              rxs = ( 672.*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -168.*(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) +32*(rsxy(i1,i2+3,i3,0,0)-rsxy(i1,i2-3,i3,0,0)) -3.*(rsxy(i1,i2+4,i3,0,0)-rsxy(i1,i2-4,i3,0,0)) )*(dr2i/840.) 
              rys = ( 672.*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -168.*(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) +32*(rsxy(i1,i2+3,i3,0,1)-rsxy(i1,i2-3,i3,0,1)) -3.*(rsxy(i1,i2+4,i3,0,1)-rsxy(i1,i2-4,i3,0,1)) )*(dr2i/840.) 
              sxs = ( 672.*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -168.*(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) +32*(rsxy(i1,i2+3,i3,1,0)-rsxy(i1,i2-3,i3,1,0)) -3.*(rsxy(i1,i2+4,i3,1,0)-rsxy(i1,i2-4,i3,1,0)) )*(dr2i/840.) 
              sys = ( 672.*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -168.*(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) +32*(rsxy(i1,i2+3,i3,1,1)-rsxy(i1,i2-3,i3,1,1)) -3.*(rsxy(i1,i2+4,i3,1,1)-rsxy(i1,i2-4,i3,1,1)) )*(dr2i/840.) 
            end if
            rxx = rx*rxr + sx*rxs 
            ryy = ry*ryr + sy*rys 
            sxx = rx*sxr + sx*sxs 
            syy = ry*syr + sy*sys 
            ! -- Coefficients in the Laplacian (scaled)
            lapCoeff(i1,i2,i3,0) = (rx**2 + ry**2   )*dr1i**2
            lapCoeff(i1,i2,i3,2) = 2.*(rx*sx + ry*sy)*dr1i*dr2i
            lapCoeff(i1,i2,i3,1) = (sx**2 + sy**2   )*dr2i**2
            lapCoeff(i1,i2,i3,3) = (rxx + ryy       )*dr1i
            lapCoeff(i1,i2,i3,4) = (sxx + syy       )*dr2i 
           end do
           end do
           end do
          ! if( c200(0,0,0).le.0. )then
          !   ! --- Evaluate and store coefficients in Laplacian ---
          !   write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'
          !   numGhost1=orderOfAccuracy/2 -1; ! check me 
          !   n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
          !   n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
          !   n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
          !   ! ----- START LOOPS ----
          !   beginLoops3d()
          !   #If "noMask" eq "USEMASK" 
          !     if( mask(i1,i2,i3).ne.0 )then
          !   #End 
          !     rx = rsxy(i1,i2,i3,0,0)
          !     ry = rsxy(i1,i2,i3,0,1)
          !     sx = rsxy(i1,i2,i3,1,0)
          !     sy = rsxy(i1,i2,i3,1,1)
          !     ! --- choose order for (r,s,t) derivatives based on available ghost points, less accuracy is needed in ghost points  ---
          !     if( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
          !       diffOrder1=2
          !     else
          !       stop 999
          !     end if
          !     if( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
          !       diffOrder2=2
          !     else
          !       stop 999
          !     end if
          !     if( diffOrder1.eq.2 )then
          !       rxr = (rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0))*(.5*dr1i) 
          !       ryr = (rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1))*(.5*dr1i) 
          !       sxr = (rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0))*(.5*dr1i) 
          !       syr = (rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1))*(.5*dr1i) 
          !     end if
          !     if( diffOrder2.eq.2 )then
          !       rxs = (rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0))*(.5*dr2i) 
          !       rys = (rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1))*(.5*dr2i) 
          !       sxs = (rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0))*(.5*dr2i) 
          !       sys = (rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1))*(.5*dr2i) 
          !     end if
          !     rxx = rx*rxr + sx*rxs 
          !     ryy = ry*ryr + sy*rys 
          !     sxx = rx*sxr + sx*sxs 
          !     syy = ry*syr + sy*sys 
          !     ! -- Coefficients in the Laplacian (scaled)
          !     c200(i1,i2,i3) = (rx**2 + ry**2   )*dr1i**2
          !     c110(i1,i2,i3) = 2.*(rx*sx + ry*sy)*dr1i*dr2i*.25
          !     c020(i1,i2,i3) = (sx**2 + sy**2   )*dr2i**2
          !     c100(i1,i2,i3) = (rxx + ryy       )*dr1i*.5
          !     c010(i1,i2,i3) = (sxx + syy       )*dr2i*.5 
          !   #If "noMask" eq "USEMASK" 
          !     end if ! mask .ne. 0
          !   #End 
          !   endLoops3d() 
          !   ! ----- END LOOPS ----
          ! end if ! end assignLapCoeff
         write(*,*) 'EVAL STENCIL COEFF'
          ! ---- 2D -----
         numGhost1=0;
         n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
         n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
         n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
         ! ----- START LOOPS ----
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
! Evaluate stencil coefficients, dim=2, order=4, gridType=Curvilinear
! File generated by cgWave/maple/writeStencilFiles.mpl
i1m1=i1-1; i1p1=i1+1;
i2m1=i2-1; i2p1=i2+1;
t1 = lapCoeff(i1m1,i2m1,i3,2);
t2 = lapCoeff(i1,i2,i3,2);
t6 = lapCoeff(i1,i2,i3,4);
t7 = lapCoeff(i1,i2m1,i3,2);
t9 = t6 * t7 / 0.8e1;
t10 = lapCoeff(i1,i2,i3,1);
t12 = t10 * t7 / 0.4e1;
t13 = lapCoeff(i1m1,i2m1,i3,4);
t14 = t13 / 0.8e1;
t15 = lapCoeff(i1m1,i2m1,i3,1);
t16 = t15 / 0.4e1;
t22 = t2 * dt2;
t23 = t22 / 0.24e2;
t27 = lapCoeff(i1,i2m1,i3,4);
t28 = t27 / 0.2e1;
t29 = lapCoeff(i1,i2m1,i3,1);
t32 = t27 / 0.4e1;
t33 = t29 / 0.2e1;
t36 = lapCoeff(i1p1,i2m1,i3,2);
t43 = lapCoeff(i1p1,i2m1,i3,4);
t44 = t43 / 0.8e1;
t45 = lapCoeff(i1p1,i2m1,i3,1);
t46 = t45 / 0.4e1;
t56 = lapCoeff(i1,i2,i3,3);
t57 = lapCoeff(i1m1,i2,i3,2);
t59 = t56 * t57 / 0.8e1;
t60 = lapCoeff(i1,i2,i3,0);
t62 = t57 * t60 / 0.4e1;
t63 = lapCoeff(i1m1,i2m1,i3,3);
t64 = t63 / 0.8e1;
t65 = lapCoeff(i1m1,i2m1,i3,0);
t66 = t65 / 0.4e1;
t73 = 0.5e1 / 0.12e2 * t22;
t74 = lapCoeff(i1m1,i2,i3,4);
t75 = t74 / 0.2e1;
t76 = lapCoeff(i1m1,i2,i3,1);
t77 = t2 / 0.2e1;
t80 = lapCoeff(i1,i2m1,i3,3);
t81 = t80 / 0.2e1;
t82 = lapCoeff(i1,i2m1,i3,0);
t85 = t74 / 0.4e1;
t86 = t76 / 0.2e1;
t89 = t80 / 0.4e1;
t90 = t82 / 0.2e1;
t99 = lapCoeff(i1p1,i2m1,i3,3);
t100 = t99 / 0.8e1;
t101 = lapCoeff(i1p1,i2m1,i3,0);
t102 = t101 / 0.4e1;
t106 = 0.2e1 * t10;
t110 = lapCoeff(i1p1,i2,i3,2);
t112 = t56 * (t57 + t110) / 0.8e1;
t113 = t57 / 0.4e1;
t114 = t110 / 0.4e1;
t122 = 0.2e1 / 0.3e1 * t6;
t123 = 0.4e1 / 0.3e1 * t10;
t127 = lapCoeff(i1p1,i2,i3,1);
t128 = lapCoeff(i1p1,i2,i3,4);
t129 = t128 / 0.2e1;
t134 = t128 / 0.4e1;
t135 = t127 / 0.2e1;
t147 = t56 * t110 / 0.8e1;
t149 = t110 * t60 / 0.4e1;
t158 = lapCoeff(i1m1,i2p1,i3,2);
t161 = lapCoeff(i1m1,i2,i3,3);
t162 = t161 / 0.4e1;
t163 = lapCoeff(i1m1,i2,i3,0);
t164 = t163 / 0.2e1;
t167 = t161 / 0.2e1;
t174 = lapCoeff(i1m1,i2p1,i3,4);
t175 = t174 / 0.8e1;
t176 = lapCoeff(i1m1,i2p1,i3,1);
t177 = t176 / 0.4e1;
t182 = 0.2e1 * t60;
t185 = lapCoeff(i1,i2p1,i3,2);
t187 = t6 * (t7 + t185) / 0.8e1;
t188 = t7 / 0.4e1;
t189 = t185 / 0.4e1;
t197 = 0.2e1 / 0.3e1 * t56;
t198 = 0.4e1 / 0.3e1 * t60;
t202 = lapCoeff(i1p1,i2p1,i3,2);
t205 = lapCoeff(i1,i2p1,i3,4);
t206 = t205 / 0.4e1;
t207 = lapCoeff(i1,i2p1,i3,1);
t208 = t207 / 0.2e1;
t211 = lapCoeff(i1p1,i2,i3,3);
t212 = t211 / 0.4e1;
t213 = lapCoeff(i1p1,i2,i3,0);
t214 = t213 / 0.2e1;
t217 = t205 / 0.2e1;
t218 = 0.4e1 * t10;
t219 = 0.4e1 * t60;
t222 = t211 / 0.2e1;
t233 = lapCoeff(i1p1,i2p1,i3,4);
t234 = t233 / 0.8e1;
t235 = lapCoeff(i1p1,i2p1,i3,1);
t236 = t235 / 0.4e1;
t263 = lapCoeff(i1m1,i2p1,i3,3);
t264 = t263 / 0.8e1;
t265 = lapCoeff(i1m1,i2p1,i3,0);
t266 = t265 / 0.4e1;
t275 = lapCoeff(i1,i2p1,i3,3);
t276 = t275 / 0.2e1;
t277 = lapCoeff(i1,i2p1,i3,0);
t282 = t275 / 0.4e1;
t283 = t277 / 0.2e1;
t294 = lapCoeff(i1p1,i2p1,i3,3);
t295 = t294 / 0.8e1;
t296 = lapCoeff(i1p1,i2p1,i3,0);
t297 = t296 / 0.4e1;
t336 = t6 * t185 / 0.8e1;
t338 = t10 * t185 / 0.4e1;
sc(1,i1,i2) = t1 * t2 * dt4 / 0.192e3;
sc(2,i1,i2) = -dt4 * (t9 - t12 + t2 * (t14 - t16)) / 0.12e2 - t23;
sc(3,i1,i2) = dt2 * (t6 - t10) / 0.12e2 - dt4 * (t10 * (t28 - t29) - t6 * (t32 - t33) + t2 * (t1 + t36) / 0.16e2) / 0.12e2;
sc(4,i1,i2) = dt4 * (t9 - t12 + t2 * (t44 - t46)) / 0.12e2 + t23;
sc(5,i1,i2) = t2 * t36 * dt4 / 0.192e3;
sc(6,i1,i2) = -dt4 * (t59 - t62 + t2 * (t64 - t66)) / 0.12e2 - t23;
sc(7,i1,i2) = t73 - dt4 * (t60 * (t75 - t76 + t77) + t10 * (t81 + t77 - t82) - t56 * (t85 - t86) - t6 * (t89 - t90) + t2 * (t15 + t65) / 0.2e1) / 0.12e2;
sc(8,i1,i2) = dt4 * (t2 * (t64 + t100 + t66 - t102) - t10 * (0.2e1 * t29 - t6 + t106 + 0.2e1 * t82) + t112 + t60 * (t6 - t106 - t113 + t114) + t6 * (t29 + t82)) / 0.12e2 - dt2 * (t122 - t123);
sc(9,i1,i2) = dt4 * (t60 * (t127 - t129 + t77) + t10 * (t81 + t77 + t82) - t56 * (t134 - t135) - t6 * (t89 + t90) + t2 * (t45 + t101) / 0.2e1) / 0.12e2 - t73;
sc(10,i1,i2) = t23 - dt4 * (t147 + t149 + t2 * (t100 + t102)) / 0.12e2;
sc(11,i1,i2) = dt2 * (t56 - t60) / 0.12e2 - dt4 * (t2 * (t1 + t158) / 0.16e2 - t56 * (t162 - t164) + t60 * (t167 - t163)) / 0.12e2;
sc(12,i1,i2) = dt4 * (t2 * (t14 + t175 + t16 - t177) - t60 * (0.2e1 * t76 - t56 + 0.2e1 * t163 + t182) + t187 + t10 * (t56 - t188 + t189 - t182) + t56 * (t76 + t163)) / 0.12e2 - dt2 * (t197 - t198);
sc(13,i1,i2) = dt4 * (t2 * (t1 + t158 + t36 + t202) / 0.16e2 - t6 * (t32 + t206 + t33 - t208) - t56 * (t162 + t212 + t164 - t214) + t10 * (t28 - t217 + t29 + t218 + t207 + t219) + t60 * (t218 + t167 - t222 + t163 + t219 + t213)) / 0.12e2 - 0.5e1 / 0.2e1 * dt2 * (t10 + t60) + 0.2e1;
sc(14,i1,i2) = dt2 * (t197 + t198) - dt4 * (t2 * (t44 + t234 + t46 - t236) + t187 + t10 * (t56 - t188 + t189 + t182) + t60 * (0.2e1 * t127 + t56 + t182 + 0.2e1 * t213) + t56 * (t127 + t213)) / 0.12e2;
sc(15,i1,i2) = dt4 * (t60 * (t222 + t213) - t2 * (t36 + t202) / 0.16e2 + t56 * (t212 + t214)) / 0.12e2 - dt2 * (t56 + t60) / 0.12e2;
sc(16,i1,i2) = dt4 * (t59 - t62 + t2 * (t264 - t266)) / 0.12e2 + t23;
sc(17,i1,i2) = dt4 * (t60 * (t75 + t76 + t77) + t10 * (t77 - t276 + t277) - t56 * (t85 + t86) - t6 * (t282 - t283) + t2 * (t176 + t265) / 0.2e1) / 0.12e2 - t73;
sc(18,i1,i2) = dt2 * (t122 + t123) - dt4 * (t2 * (t264 + t295 + t266 - t297) + t112 + t10 * (t6 + t106 + 0.2e1 * t207 + 0.2e1 * t277) + t60 * (t6 + t106 - t113 + t114) + t6 * (t207 + t277)) / 0.12e2;
sc(19,i1,i2) = t73 + dt4 * (t60 * (t129 + t127 - t77) + t10 * (t276 - t77 + t277) + t56 * (t134 + t135) + t6 * (t282 + t283) - t2 * (t235 + t296) / 0.2e1) / 0.12e2;
sc(20,i1,i2) = dt4 * (t147 + t149 + t2 * (t295 + t297)) / 0.12e2 - t23;
sc(21,i1,i2) = t158 * t2 * dt4 / 0.192e3;
sc(22,i1,i2) = t23 - dt4 * (t336 + t338 + t2 * (t175 + t177)) / 0.12e2;
sc(23,i1,i2) = dt4 * (t10 * (t217 + t207) + t6 * (t206 + t208) - t2 * (t158 + t202) / 0.16e2) / 0.12e2 - dt2 * (t6 + t10) / 0.12e2;
sc(24,i1,i2) = dt4 * (t336 + t338 + t2 * (t234 + t236)) / 0.12e2 - t23;
sc(25,i1,i2) = t2 * t202 * dt4 / 0.192e3;
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
         !   updateWaveOpt(2,4,2,curvilinear,NOFORCING)
         ! else 
         !   updateWaveOpt(2,4,2,curvilinear,FORCING)
         ! end if
       else
         if( addForcing.eq.0 )then
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=2 order=4 orderInTime=4, grid=curvilinear... t=",e10.2)') t
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
! Stencil: nd=2, orderOfAccuracy=4, gridType=Curvilinear
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ sc(  1,i1,i2)*u(i1-2,i2-2,i3,m) + sc(  2,i1,i2)*u(i1-1,i2-2,i3,m) + sc(  3,i1,i2)*u(i1+0,i2-2,i3,m) + sc(  4,i1,i2)*u(i1+1,i2-2,i3,m) + sc(  5,i1,i2)*u(i1+2,i2-2,i3,m)+ sc(  6,i1,i2)*u(i1-2,i2-1,i3,m) + sc(  7,i1,i2)*u(i1-1,i2-1,i3,m) + sc(  8,i1,i2)*u(i1+0,i2-1,i3,m) + sc(  9,i1,i2)*u(i1+1,i2-1,i3,m) + sc( 10,i1,i2)*u(i1+2,i2-1,i3,m)+ sc( 11,i1,i2)*u(i1-2,i2+0,i3,m) + sc( 12,i1,i2)*u(i1-1,i2+0,i3,m) + sc( 13,i1,i2)*u(i1+0,i2+0,i3,m) + sc( 14,i1,i2)*u(i1+1,i2+0,i3,m) + sc( 15,i1,i2)*u(i1+2,i2+0,i3,m)+ sc( 16,i1,i2)*u(i1-2,i2+1,i3,m) + sc( 17,i1,i2)*u(i1-1,i2+1,i3,m) + sc( 18,i1,i2)*u(i1+0,i2+1,i3,m) + sc( 19,i1,i2)*u(i1+1,i2+1,i3,m) + sc( 20,i1,i2)*u(i1+2,i2+1,i3,m)+ sc( 21,i1,i2)*u(i1-2,i2+2,i3,m) + sc( 22,i1,i2)*u(i1-1,i2+2,i3,m) + sc( 23,i1,i2)*u(i1+0,i2+2,i3,m) + sc( 24,i1,i2)*u(i1+1,i2+2,i3,m) + sc( 25,i1,i2)*u(i1+2,i2+2,i3,m)
               end do
               end do
               end do
             ! endLoopsMask()
         else
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=2 order=4 orderInTime=4, grid=curvilinear... t=",e10.2)') t
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
! Stencil: nd=2, orderOfAccuracy=4, gridType=Curvilinear
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ sc(  1,i1,i2)*u(i1-2,i2-2,i3,m) + sc(  2,i1,i2)*u(i1-1,i2-2,i3,m) + sc(  3,i1,i2)*u(i1+0,i2-2,i3,m) + sc(  4,i1,i2)*u(i1+1,i2-2,i3,m) + sc(  5,i1,i2)*u(i1+2,i2-2,i3,m)+ sc(  6,i1,i2)*u(i1-2,i2-1,i3,m) + sc(  7,i1,i2)*u(i1-1,i2-1,i3,m) + sc(  8,i1,i2)*u(i1+0,i2-1,i3,m) + sc(  9,i1,i2)*u(i1+1,i2-1,i3,m) + sc( 10,i1,i2)*u(i1+2,i2-1,i3,m)+ sc( 11,i1,i2)*u(i1-2,i2+0,i3,m) + sc( 12,i1,i2)*u(i1-1,i2+0,i3,m) + sc( 13,i1,i2)*u(i1+0,i2+0,i3,m) + sc( 14,i1,i2)*u(i1+1,i2+0,i3,m) + sc( 15,i1,i2)*u(i1+2,i2+0,i3,m)+ sc( 16,i1,i2)*u(i1-2,i2+1,i3,m) + sc( 17,i1,i2)*u(i1-1,i2+1,i3,m) + sc( 18,i1,i2)*u(i1+0,i2+1,i3,m) + sc( 19,i1,i2)*u(i1+1,i2+1,i3,m) + sc( 20,i1,i2)*u(i1+2,i2+1,i3,m)+ sc( 21,i1,i2)*u(i1-2,i2+2,i3,m) + sc( 22,i1,i2)*u(i1-1,i2+2,i3,m) + sc( 23,i1,i2)*u(i1+0,i2+2,i3,m) + sc( 24,i1,i2)*u(i1+1,i2+2,i3,m) + sc( 25,i1,i2)*u(i1+2,i2+2,i3,m)+dtSq*fv(m)
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
