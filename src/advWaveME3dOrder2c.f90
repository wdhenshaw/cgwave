! This file automatically generated from advWaveME.bf90 with bpp.
    subroutine advWaveME3dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
  ! subroutine advWaveME3dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,!                 mask,xy,rsxy,  um,u,un, f,fa, v, vh,  bc, frequencyArray, ipar, rpar, ierr )
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
    real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds coeff of Laplacian for HA scheme
    real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
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
          real dr1, dr2, dr3, dr1i, dr2i, dr3i, rx, ry, rz, sx, sy
          real sz, tx, ty, tz, diffOrder1, diffOrder2, diffOrder3, rxr, rxs, rxt, ryr
          real rys, ryt, rzr, rzs, rzt, sxr, sxs, sxt, syr, sys, syt
          real szr, szs, szt, txr, txs, txt, tyr, tys, tyt, tzr, tzs
          real tzt, rxx, ryy, rzz, sxx, syy, szz, txx, tyy, tzz, d200i
          real d020i, d002i, lap2h, d100i, d010i, d110i, d001i, d101i, d011i
      integer maxDeriv,d,uc,count,numGhost1,m1,m2,m3
 ! declare coefficients in the chain rule for curvilinear grids (from cgwave/maple/chainRuleCoefficients.mw)
 ! #If "curvilinear" eq "curvilinear"
 !   #If 3 == 2
 !     #Include "../maple/declareChainRuleCoefficients2d.h"
 !   #Else
 !     #Include "../maple/declareChainRuleCoefficients3d.h"
 !   #End
 ! #End
  ! statement functions for coefficients
  ! real c200,c020,c002, c110, c101, c011, c100, c010, cux001
  ! #If 3 == 2
  !   c200(i1,i2,i3) = lapCoeff(i1,i2,i3,0)
  !   c020(i1,i2,i3) = lapCoeff(i1,i2,i3,1)
  !   c110(i1,i2,i3) = lapCoeff(i1,i2,i3,2)
  !   c100(i1,i2,i3) = lapCoeff(i1,i2,i3,3)
  !   c010(i1,i2,i3) = lapCoeff(i1,i2,i3,4)
  ! #Else
  !   c200(i1,i2,i3) = lapCoeff(i1,i2,i3,0)
  !   c020(i1,i2,i3) = lapCoeff(i1,i2,i3,1)
  !   c002(i1,i2,i3) = lapCoeff(i1,i2,i3,2)
  !   c110(i1,i2,i3) = lapCoeff(i1,i2,i3,3)
  !   c101(i1,i2,i3) = lapCoeff(i1,i2,i3,4)
  !   c011(i1,i2,i3) = lapCoeff(i1,i2,i3,5)
  !   c100(i1,i2,i3) = lapCoeff(i1,i2,i3,6)
  !   c010(i1,i2,i3) = lapCoeff(i1,i2,i3,7)   
  !   c001(i1,i2,i3) = lapCoeff(i1,i2,i3,8)   
  ! #End
 ! OLD  #Include "derivativesByChainRuleCoefficients.h"
   !...........end   statement functions
   ! write(*,*) 'Inside advWaveME...'
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
                        write(*,'("advWaveME: ADVANCE dim=3 order=2 orderInTime=2, grid=curvilinear... t=",e10.2)') t
                    end if
                    m=0 ! component number 
                    ec = 0 ! component number 
          ! -- call the appropriate macro:
          !  update2dOrder2Rectangular(3,2,2,curvilinear)
          !  update3dOrder6Curvilinear(3,2,2,curvilinear)
                    if( useMask.eq.0 .and. addForcing.eq.0 )then
            ! No-mask, no-forcing
              ! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
              ! Example: 
              ! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...
                            dr1=dr(0); dr1i=1./dr1;
                            dr2=dr(1); dr2i=1./dr2;
                            dr3=dr(2); dr3i=1./dr3;
                            fv(m)=0.
                            if( lapCoeff(0,0,0,0).le.0. )then
                ! --- Evaluate and store coefficients in Laplacian ---
                                write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'
                                numGhost1=0;
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
                                    if( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
                                        diffOrder1=2
                                    else
                                        stop 999
                                    end if
                                    if( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
                                        diffOrder2=2
                                    else
                                        stop 999
                                    end if
                                    if( (i3-1).ge.nd3a .and. (i3+1).le.nd3b )then
                                        diffOrder3=2
                                    else
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
                                    lapCoeff(i1,i2,i3,3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i*.25
                                    lapCoeff(i1,i2,i3,4) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr3i*.25
                                    lapCoeff(i1,i2,i3,5) = 2.*(sx*tx + sy*ty + sz*tz )*dr2i*dr3i*.25
                                    lapCoeff(i1,i2,i3,6) = (rxx + ryy + rzz)*dr1i*.5
                                    lapCoeff(i1,i2,i3,7) = (sxx + syy + szz)*dr2i*.5 
                                    lapCoeff(i1,i2,i3,8) = (txx + tyy + tzz)*dr3i*.5 
                                  end do
                                  end do
                                  end do
                            end if ! end assignLapCoeff
              ! ===========  FINAL LOOP TO FILL IN THE SOLUTION ============
                            numGhost1=0;
                            n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                            n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                            n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                    d200i = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
                                    d020i = u(i1,i2+1,i3,0) - 2*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
                                    d002i = u(i1,i2,i3+1,0) - 2*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
                                    d100i = u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0)
                                    d010i = u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0)
                                    d110i = u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0)
                                    d001i = u(i1,i2,i3+1,0) - u(i1,i2,i3-1,0)
                                    d101i = u(i1+1,i2,i3+1,0) - u(i1-1,i2,i3+1,0) - u(i1+1,i2,i3-1,0) + u(i1-1,i2,i3-1,0)
                                    d011i = u(i1,i2+1,i3+1,0) - u(i1,i2-1,i3+1,0) - u(i1,i2+1,i3-1,0) + u(i1,i2-1,i3-1,0)
                                    lap2h = lapCoeff(i1,i2,i3,0)*d200i +lapCoeff(i1,i2,i3,1)*d020i +lapCoeff(i1,i2,i3,2)*d002i +lapCoeff(i1,i2,i3,3)*d110i + lapCoeff(i1,i2,i3,4)*d101i + lapCoeff(i1,i2,i3,5)*d011i + lapCoeff(i1,i2,i3,6)*d100i + lapCoeff(i1,i2,i3,7)*d010i + lapCoeff(i1,i2,i3,8)*d001i
                  ! --- Modified equation space-time update ----
                                    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( lap2h )                         +dtSq*fv(m)                             
                                  end do
                                  end do
                                  end do
                    else
              ! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
              ! Example: 
              ! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...
                            dr1=dr(0); dr1i=1./dr1;
                            dr2=dr(1); dr2i=1./dr2;
                            dr3=dr(2); dr3i=1./dr3;
                            fv(m)=0.
                            if( lapCoeff(0,0,0,0).le.0. )then
                ! --- Evaluate and store coefficients in Laplacian ---
                                write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'
                                numGhost1=0;
                                n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                if( mask(i1,i2,i3).ne.0 )then
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
                                    if( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
                                        diffOrder1=2
                                    else
                                        stop 999
                                    end if
                                    if( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
                                        diffOrder2=2
                                    else
                                        stop 999
                                    end if
                                    if( (i3-1).ge.nd3a .and. (i3+1).le.nd3b )then
                                        diffOrder3=2
                                    else
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
                                    lapCoeff(i1,i2,i3,3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i*.25
                                    lapCoeff(i1,i2,i3,4) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr3i*.25
                                    lapCoeff(i1,i2,i3,5) = 2.*(sx*tx + sy*ty + sz*tz )*dr2i*dr3i*.25
                                    lapCoeff(i1,i2,i3,6) = (rxx + ryy + rzz)*dr1i*.5
                                    lapCoeff(i1,i2,i3,7) = (sxx + syy + szz)*dr2i*.5 
                                    lapCoeff(i1,i2,i3,8) = (txx + tyy + tzz)*dr3i*.5 
                                    end if ! mask .ne. 0
                                  end do
                                  end do
                                  end do
                            end if ! end assignLapCoeff
              ! ===========  FINAL LOOP TO FILL IN THE SOLUTION ============
                            numGhost1=0;
                            n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                            n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                            n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                if( mask(i1,i2,i3).ne.0 )then
                                    d200i = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
                                    d020i = u(i1,i2+1,i3,0) - 2*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
                                    d002i = u(i1,i2,i3+1,0) - 2*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
                                    d100i = u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0)
                                    d010i = u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0)
                                    d110i = u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0)
                                    d001i = u(i1,i2,i3+1,0) - u(i1,i2,i3-1,0)
                                    d101i = u(i1+1,i2,i3+1,0) - u(i1-1,i2,i3+1,0) - u(i1+1,i2,i3-1,0) + u(i1-1,i2,i3-1,0)
                                    d011i = u(i1,i2+1,i3+1,0) - u(i1,i2-1,i3+1,0) - u(i1,i2+1,i3-1,0) + u(i1,i2-1,i3-1,0)
                                    lap2h = lapCoeff(i1,i2,i3,0)*d200i +lapCoeff(i1,i2,i3,1)*d020i +lapCoeff(i1,i2,i3,2)*d002i +lapCoeff(i1,i2,i3,3)*d110i + lapCoeff(i1,i2,i3,4)*d101i + lapCoeff(i1,i2,i3,5)*d011i + lapCoeff(i1,i2,i3,6)*d100i + lapCoeff(i1,i2,i3,7)*d010i + lapCoeff(i1,i2,i3,8)*d001i
                                            if( forcingOption.eq.twilightZoneForcing )then
                                                        call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                                        call ogDeriv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evtt(m) )
                                                        call ogDeriv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxx(m) )
                                                        call ogDeriv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evyy(m) )
                                                        call ogDeriv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evzz(m) )
                                                    fv(m) = evtt(m) - csq*( evxx(m) + evyy(m)  + evzz(m) )
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
                  ! --- Modified equation space-time update ----
                                    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( lap2h )                         +dtSq*fv(m)                             
                                    end if ! mask .ne. 0
                                  end do
                                  end do
                                  end do
                    end if
              else
                      if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
                          write(*,'("advWaveME: ADVANCE dim=3 order=2 orderInTime=2, grid=curvilinear... t=",e10.2)') t
                      end if
                      m=0 ! component number 
                      ec = 0 ! component number 
           ! -- call the appropriate macro:
           !  update2dOrder2Rectangular(3,2,2,curvilinear)
           !  update3dOrder6Curvilinear(3,2,2,curvilinear)
                      if( useMask.eq.0 .and. addForcing.eq.0 )then
             ! No-mask, no-forcing
               ! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
               ! Example: 
               ! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...
                              dr1=dr(0); dr1i=1./dr1;
                              dr2=dr(1); dr2i=1./dr2;
                              dr3=dr(2); dr3i=1./dr3;
                              fv(m)=0.
                              if( lapCoeff(0,0,0,0).le.0. )then
                 ! --- Evaluate and store coefficients in Laplacian ---
                                  write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'
                                  numGhost1=0;
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
                                      if( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
                                          diffOrder1=2
                                      else
                                          stop 999
                                      end if
                                      if( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
                                          diffOrder2=2
                                      else
                                          stop 999
                                      end if
                                      if( (i3-1).ge.nd3a .and. (i3+1).le.nd3b )then
                                          diffOrder3=2
                                      else
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
                                      lapCoeff(i1,i2,i3,3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i*.25
                                      lapCoeff(i1,i2,i3,4) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr3i*.25
                                      lapCoeff(i1,i2,i3,5) = 2.*(sx*tx + sy*ty + sz*tz )*dr2i*dr3i*.25
                                      lapCoeff(i1,i2,i3,6) = (rxx + ryy + rzz)*dr1i*.5
                                      lapCoeff(i1,i2,i3,7) = (sxx + syy + szz)*dr2i*.5 
                                      lapCoeff(i1,i2,i3,8) = (txx + tyy + tzz)*dr3i*.5 
                                    end do
                                    end do
                                    end do
                              end if ! end assignLapCoeff
               ! ===========  FINAL LOOP TO FILL IN THE SOLUTION ============
                              numGhost1=0;
                              n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                              n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                              n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                do i3=n3a,n3b
                                do i2=n2a,n2b
                                do i1=n1a,n1b
                                      d200i = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
                                      d020i = u(i1,i2+1,i3,0) - 2*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
                                      d002i = u(i1,i2,i3+1,0) - 2*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
                                      d100i = u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0)
                                      d010i = u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0)
                                      d110i = u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0)
                                      d001i = u(i1,i2,i3+1,0) - u(i1,i2,i3-1,0)
                                      d101i = u(i1+1,i2,i3+1,0) - u(i1-1,i2,i3+1,0) - u(i1+1,i2,i3-1,0) + u(i1-1,i2,i3-1,0)
                                      d011i = u(i1,i2+1,i3+1,0) - u(i1,i2-1,i3+1,0) - u(i1,i2+1,i3-1,0) + u(i1,i2-1,i3-1,0)
                                      lap2h = lapCoeff(i1,i2,i3,0)*d200i +lapCoeff(i1,i2,i3,1)*d020i +lapCoeff(i1,i2,i3,2)*d002i +lapCoeff(i1,i2,i3,3)*d110i + lapCoeff(i1,i2,i3,4)*d101i + lapCoeff(i1,i2,i3,5)*d011i + lapCoeff(i1,i2,i3,6)*d100i + lapCoeff(i1,i2,i3,7)*d010i + lapCoeff(i1,i2,i3,8)*d001i
                   ! --- Modified equation space-time update ----
                                      un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( lap2h )                         +dtSq*fv(m)                             
                                    end do
                                    end do
                                    end do
                      else
               ! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
               ! Example: 
               ! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...
                              dr1=dr(0); dr1i=1./dr1;
                              dr2=dr(1); dr2i=1./dr2;
                              dr3=dr(2); dr3i=1./dr3;
                              fv(m)=0.
                              if( lapCoeff(0,0,0,0).le.0. )then
                 ! --- Evaluate and store coefficients in Laplacian ---
                                  write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'
                                  numGhost1=0;
                                  n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                                  n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                                  n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                do i3=n3a,n3b
                                do i2=n2a,n2b
                                do i1=n1a,n1b
                                  if( mask(i1,i2,i3).ne.0 )then
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
                                      if( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
                                          diffOrder1=2
                                      else
                                          stop 999
                                      end if
                                      if( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
                                          diffOrder2=2
                                      else
                                          stop 999
                                      end if
                                      if( (i3-1).ge.nd3a .and. (i3+1).le.nd3b )then
                                          diffOrder3=2
                                      else
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
                                      lapCoeff(i1,i2,i3,3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i*.25
                                      lapCoeff(i1,i2,i3,4) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr3i*.25
                                      lapCoeff(i1,i2,i3,5) = 2.*(sx*tx + sy*ty + sz*tz )*dr2i*dr3i*.25
                                      lapCoeff(i1,i2,i3,6) = (rxx + ryy + rzz)*dr1i*.5
                                      lapCoeff(i1,i2,i3,7) = (sxx + syy + szz)*dr2i*.5 
                                      lapCoeff(i1,i2,i3,8) = (txx + tyy + tzz)*dr3i*.5 
                                      end if ! mask .ne. 0
                                    end do
                                    end do
                                    end do
                              end if ! end assignLapCoeff
               ! ===========  FINAL LOOP TO FILL IN THE SOLUTION ============
                              numGhost1=0;
                              n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                              n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                              n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                do i3=n3a,n3b
                                do i2=n2a,n2b
                                do i1=n1a,n1b
                                  if( mask(i1,i2,i3).ne.0 )then
                                      d200i = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
                                      d020i = u(i1,i2+1,i3,0) - 2*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
                                      d002i = u(i1,i2,i3+1,0) - 2*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
                                      d100i = u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0)
                                      d010i = u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0)
                                      d110i = u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0)
                                      d001i = u(i1,i2,i3+1,0) - u(i1,i2,i3-1,0)
                                      d101i = u(i1+1,i2,i3+1,0) - u(i1-1,i2,i3+1,0) - u(i1+1,i2,i3-1,0) + u(i1-1,i2,i3-1,0)
                                      d011i = u(i1,i2+1,i3+1,0) - u(i1,i2-1,i3+1,0) - u(i1,i2+1,i3-1,0) + u(i1,i2-1,i3-1,0)
                                      lap2h = lapCoeff(i1,i2,i3,0)*d200i +lapCoeff(i1,i2,i3,1)*d020i +lapCoeff(i1,i2,i3,2)*d002i +lapCoeff(i1,i2,i3,3)*d110i + lapCoeff(i1,i2,i3,4)*d101i + lapCoeff(i1,i2,i3,5)*d011i + lapCoeff(i1,i2,i3,6)*d100i + lapCoeff(i1,i2,i3,7)*d010i + lapCoeff(i1,i2,i3,8)*d001i
                                              if( forcingOption.eq.twilightZoneForcing )then
                                                          call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                                          call ogDeriv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evtt(m) )
                                                          call ogDeriv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxx(m) )
                                                          call ogDeriv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evyy(m) )
                                                          call ogDeriv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evzz(m) )
                                                      fv(m) = evtt(m) - csq*( evxx(m) + evyy(m)  + evzz(m) )
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
                   ! --- Modified equation space-time update ----
                                      un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) + cdtsq*( lap2h )                         +dtSq*fv(m)                             
                                      end if ! mask .ne. 0
                                    end do
                                    end do
                                    end do
                      end if
              end if 
      else
     ! --- IMPLICIT: Fill in RHS to implicit time-stepping -----
          stop 1111
      end if
      return
      end
