! This file automatically generated from advWaveME.bf90 with bpp.
    subroutine advWaveME3dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
  ! subroutine advWaveME3dOrder8c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,!                 mask,xy,rsxy,  um,u,un, f,fa, v, vh,  bc, frequencyArray, ipar, rpar, ierr )
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
          real cr0, cr1, cr2, cr3, cs0, cs1, cs2, cs3, ct0, ct1, ct2
          real ct3, crr0, crr1, crr2, crr3, css0, css1, css2, css3, ctt0, ctt1
          real ctt2, ctt3, crrr0, crrr1, crrr2, crrr3, csss0, csss1, csss2, csss3, cttt0
          real cttt1, cttt2, cttt3, crrrr0, crrrr1, crrrr2, crrrr3, cssss0, cssss1, cssss2, cssss3
          real ctttt0, ctttt1, ctttt2, ctttt3, crrrrr0, crrrrr1, crrrrr2, crrrrr3, csssss0, csssss1, csssss2
          real csssss3, cttttt0, cttttt1, cttttt2, cttttt3, crrrrrr0, crrrrrr1, crrrrrr2, crrrrrr3, cssssss0, cssssss1
          real cssssss2, cssssss3, ctttttt0, ctttttt1, ctttttt2, ctttttt3, dr1, dr2, dr3, dr1i, dr2i
          real dr3i, rx, ry, rz, sx, sy, sz, tx, ty, tz, diffOrder1
          real diffOrder2, diffOrder3, rxr, rxs, rxt, ryr, rys, ryt, rzr, rzs, rzt
          real sxr, sxs, sxt, syr, sys, syt, szr, szs, szt, txr, txs
          real txt, tyr, tys, tyt, tzr, tzs, tzt, rxx, ryy, rzz, sxx
          real syy, szz, txx, tyy, tzz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d100i, d010i
          real d110i, d001i, d101i, d011i, d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d004(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d202(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d022(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d300i, d030i, d003i, d310i, d130i, d301i, d103i, d031i, d013i, lap2h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real lap2h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h100i, lap2h010i, lap2h110i, lap2h001i, lap2h101i, lap2h011i, lap6h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hCubed(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real d600(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d060(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d006(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d420(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d240(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d402(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d204(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d042(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d024(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d500i, d050i
          real d005i, d510i, d150i, d330i, d501i, d105i, d051i, d015i, d303i, d033i, lap2hSq200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real lap2hSq020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq100i, lap2hSq010i, lap2hSq001i, lap2hSq110i, lap2hSq101i, lap2hSq011i, lap4h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
          real lap2h400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h004(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h202(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h022(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h100i, lap4h010i, lap4h001i, lap4h110i, lap4h101i
          real lap4h011i, lap2h300i, lap2h030i, lap2h003i, lap2h310i, lap2h130i, lap2h301i, lap2h103i, lap2h031i, lap2h013i, d800i
          real d080i, d008i, d700i, d070i, d007i, d710i, d170i, d701i, d107i, d071i, d017i
          real d530i, d350i, d503i, d305i, d053i, d035i, lap8h, lap2hCubed200i, lap2hCubed020i, lap2hCubed002i, lap2hCubed100i
          real lap2hCubed010i, lap2hCubed001i, lap2hCubed110i, lap2hCubed101i, lap2hCubed011i, lap2h4p, lap6h200i, lap6h020i, lap6h002i, lap6h100i, lap6h010i
          real lap6h001i, lap6h110i, lap6h101i, lap6h011i, lap4h400i, lap4h040i, lap4h004i, lap4h300i, lap4h030i, lap4h003i, lap4h310i
          real lap4h130i, lap4h301i, lap4h103i, lap4h031i, lap4h013i, lap2h600i, lap2h060i, lap2h006i, lap2h500i, lap2h050i, lap2h005i
          real lap2h510i, lap2h150i, lap2h330i, lap2h501i, lap2h105i, lap2h051i, lap2h015i, lap2h303i, lap2h033i, lap6hSq, lap4hSq200i
          real lap4hSq020i, lap4hSq002i, lap4hSq100i, lap4hSq010i, lap4hSq001i, lap4hSq110i, lap4hSq101i, lap4hSq011i, lap2hSq400i, lap2hSq040i, lap2hSq004i
          real lap2hSq300i, lap2hSq030i, lap2hSq003i, lap2hSq310i, lap2hSq130i, lap2hSq301i, lap2hSq103i, lap2hSq031i, lap2hSq013i, lap4hCubed
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
                        write(*,'("advWaveME: ADVANCE dim=3 order=8 orderInTime=2, grid=curvilinear... t=",e10.2)') t
                    end if
                    m=0 ! component number 
                    ec = 0 ! component number 
          ! -- call the appropriate macro:
          !  update2dOrder2Rectangular(3,8,2,curvilinear)
          !  update3dOrder6Curvilinear(3,8,2,curvilinear)
                    if( useMask.eq.0 .and. addForcing.eq.0 )then
            ! No-mask, no-forcing
              ! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
              ! Example: 
              ! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...
cr0 = 1.; cr1 = -1/6.; cr2 = 1/30.; cr3 = -1/140.; 
cs0 = 1.; cs1 = -1/6.; cs2 = 1/30.; cs3 = -1/140.; 
ct0 = 1.; ct1 = -1/6.; ct2 = 1/30.; ct3 = -1/140.; 
crr0 = 1.; crr1 = -1/12.; crr2 = 1/90.; crr3 = -1/560.; 
css0 = 1.; css1 = -1/12.; css2 = 1/90.; css3 = -1/560.; 
ctt0 = 1.; ctt1 = -1/12.; ctt2 = 1/90.; ctt3 = -1/560.; 
crrr0 = 1.; crrr1 = -1/4.; crrr2 = 7/120.; crrr3 = -41/3024.; 
csss0 = 1.; csss1 = -1/4.; csss2 = 7/120.; csss3 = -41/3024.; 
cttt0 = 1.; cttt1 = -1/4.; cttt2 = 7/120.; cttt3 = -41/3024.; 
crrrr0 = 1.; crrrr1 = -1/6.; crrrr2 = 7/240.; crrrr3 = -41/7560.; 
cssss0 = 1.; cssss1 = -1/6.; cssss2 = 7/240.; cssss3 = -41/7560.; 
ctttt0 = 1.; ctttt1 = -1/6.; ctttt2 = 7/240.; ctttt3 = -41/7560.; 
crrrrr0 = 1.; crrrrr1 = -1/3.; crrrrr2 = 13/144.; crrrrr3 = -139/6048.; 
csssss0 = 1.; csssss1 = -1/3.; csssss2 = 13/144.; csssss3 = -139/6048.; 
cttttt0 = 1.; cttttt1 = -1/3.; cttttt2 = 13/144.; cttttt3 = -139/6048.; 
crrrrrr0 = 1.; crrrrrr1 = -1/4.; crrrrrr2 = 13/240.; crrrrrr3 = -139/12096.; 
cssssss0 = 1.; cssssss1 = -1/4.; cssssss2 = 13/240.; cssssss3 = -139/12096.; 
ctttttt0 = 1.; ctttttt1 = -1/4.; ctttttt2 = 13/240.; ctttttt3 = -139/12096.; 
                            dr1=dr(0); dr1i=1./dr1;
                            dr2=dr(1); dr2i=1./dr2;
                            dr3=dr(2); dr3i=1./dr3;
                            fv(m)=0.
                            if( lapCoeff(0,0,0,0).le.0. )then
                ! --- Evaluate and store coefficients in Laplacian ---
                                write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'
                                numGhost1=3;
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
                                        diffOrder1=8
                                    elseif( (i1-3).ge.nd1a .and. (i1+3).le.nd1b )then
                                        diffOrder1=6
                                    elseif( (i1-2).ge.nd1a .and. (i1+2).le.nd1b )then
                                        diffOrder1=4
                                    elseif( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
                                        diffOrder1=2
                                    else
                                        stop 999
                                    end if
                                    if( (i2-4).ge.nd2a .and. (i2+4).le.nd2b )then
                                        diffOrder2=8
                                    elseif( (i2-3).ge.nd2a .and. (i2+3).le.nd2b )then
                                        diffOrder2=6
                                    elseif( (i2-2).ge.nd2a .and. (i2+2).le.nd2b )then
                                        diffOrder2=4
                                    elseif( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
                                        diffOrder2=2
                                    else
                                        stop 999
                                    end if
                                    if( (i3-4).ge.nd3a .and. (i3+4).le.nd3b )then
                                        diffOrder3=8
                                    elseif( (i3-3).ge.nd3a .and. (i3+3).le.nd3b )then
                                        diffOrder3=6
                                    elseif( (i3-2).ge.nd3a .and. (i3+2).le.nd3b )then
                                        diffOrder3=4
                                    elseif( (i3-1).ge.nd3a .and. (i3+1).le.nd3b )then
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
                                        rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
                                        rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
                                        rzs = ( 8*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) )*(dr2i/12.) 
                                        sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
                                        sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
                                        szs = ( 8*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) )*(dr2i/12.) 
                                        txs = ( 8*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) )*(dr2i/12.) 
                                        tys = ( 8*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) )*(dr2i/12.) 
                                        tzs = ( 8*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) )*(dr2i/12.) 
                                    elseif( diffOrder2.eq.8 )then
                                        rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
                                        rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
                                        rzs = ( 8*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) )*(dr2i/12.) 
                                        sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
                                        sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
                                        szs = ( 8*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) )*(dr2i/12.) 
                                        txs = ( 8*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) )*(dr2i/12.) 
                                        tys = ( 8*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) )*(dr2i/12.) 
                                        tzs = ( 8*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) )*(dr2i/12.) 
                                    end if
                                    if( diffOrder3.eq.2 )then
                                        rxt = (rsxy(i1,i2,i3+1,0,0)-rsxy(i1,i2,i3-1,0,0))*(.5*dr2i) 
                                        ryt = (rsxy(i1,i2,i3+1,0,1)-rsxy(i1,i2,i3-1,0,1))*(.5*dr2i) 
                                        rzt = (rsxy(i1,i2,i3+1,0,2)-rsxy(i1,i2,i3-1,0,2))*(.5*dr2i) 
                                        sxt = (rsxy(i1,i2,i3+1,1,0)-rsxy(i1,i2,i3-1,1,0))*(.5*dr2i) 
                                        syt = (rsxy(i1,i2,i3+1,1,1)-rsxy(i1,i2,i3-1,1,1))*(.5*dr2i) 
                                        szt = (rsxy(i1,i2,i3+1,1,2)-rsxy(i1,i2,i3-1,1,2))*(.5*dr2i) 
                                        txt = (rsxy(i1,i2,i3+1,2,0)-rsxy(i1,i2,i3-1,2,0))*(.5*dr2i) 
                                        tyt = (rsxy(i1,i2,i3+1,2,1)-rsxy(i1,i2,i3-1,2,1))*(.5*dr2i) 
                                        tzt = (rsxy(i1,i2,i3+1,2,2)-rsxy(i1,i2,i3-1,2,2))*(.5*dr2i) 
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
                                    lapCoeff(i1,i2,i3,3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i*.25
                                    lapCoeff(i1,i2,i3,4) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr2i*.25
                                    lapCoeff(i1,i2,i3,5) = 2.*(sx*tx + sy*ty + sz*tz )*dr1i*dr2i*.25
                                    lapCoeff(i1,i2,i3,6) = (rxx + ryy + rzz)*dr1i*.5
                                    lapCoeff(i1,i2,i3,7) = (sxx + syy + tyy)*dr2i*.5 
                                    lapCoeff(i1,i2,i3,8) = (txx + tyy + tzz)*dr3i*.5 
                                  end do
                                  end do
                                  end do
                            end if ! end assignLapCoeff
                            numGhost1=3;
                            n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                            n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                            n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                    d200(i1,i2,i3,0) = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
                                    d020(i1,i2,i3,0) = u(i1,i2+1,i3,0) - 2*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
                                    d002(i1,i2,i3,0) = u(i1,i2,i3+1,0) - 2*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
                                    d100i = u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0)
                                    d010i = u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0)
                                    d110i = u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0)
                                    d001i = u(i1,i2,i3+1,0) - u(i1,i2,i3-1,0)
                                    d101i = u(i1+1,i2,i3+1,0) - u(i1-1,i2,i3+1,0) - u(i1+1,i2,i3-1,0) + u(i1-1,i2,i3-1,0)
                                    d011i = u(i1,i2+1,i3+1,0) - u(i1,i2-1,i3+1,0) - u(i1,i2+1,i3-1,0) + u(i1,i2-1,i3-1,0)
                                    lap2h(i1,i2,i3,0) = lapCoeff(i1,i2,i3,0)*d200(i1,i2,i3,0) +lapCoeff(i1,i2,i3,1)*d020(i1,i2,i3,0) +lapCoeff(i1,i2,i3,2)*d002(i1,i2,i3,0) +lapCoeff(i1,i2,i3,3)*d110i + lapCoeff(i1,i2,i3,4)*d101i + lapCoeff(i1,i2,i3,5)*d011i + lapCoeff(i1,i2,i3,6)*d100i + lapCoeff(i1,i2,i3,7)*d010i + lapCoeff(i1,i2,i3,8)*d001i
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
                                    d400(i1,i2,i3,0) = d200(i1+1,i2,i3,0) - 2*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0)
                                    d040(i1,i2,i3,0) = d020(i1,i2+1,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1,i2-1,i3,0)
                                    d004(i1,i2,i3,0) = d002(i1,i2,i3+1,0) - 2*d002(i1,i2,i3,0) + d002(i1,i2,i3-1,0)
                                    d220(i1,i2,i3,0) = d020(i1+1,i2,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1-1,i2,i3,0)
                                    d202(i1,i2,i3,0) = d002(i1+1,i2,i3,0) - 2*d002(i1,i2,i3,0) + d002(i1-1,i2,i3,0)
                                    d022(i1,i2,i3,0) = d002(i1,i2+1,i3,0) - 2*d002(i1,i2,i3,0) + d002(i1,i2-1,i3,0)
                                    d300i = d200(i1+1,i2,i3,0) - d200(i1-1,i2,i3,0)
                                    d030i = d020(i1,i2+1,i3,0) - d020(i1,i2-1,i3,0)
                                    d003i = d002(i1,i2,i3+1,0) - d002(i1,i2,i3-1,0)
                                    d310i = d200(i1+1,i2+1,i3,0) - d200(i1-1,i2+1,i3,0) - d200(i1+1,i2-1,i3,0) + d200(i1-1,i2-1,i3,0)
                                    d130i = d020(i1+1,i2+1,i3,0) - d020(i1-1,i2+1,i3,0) - d020(i1+1,i2-1,i3,0) + d020(i1-1,i2-1,i3,0)
                                    d301i = d200(i1+1,i2,i3+1,0) - d200(i1-1,i2,i3+1,0) - d200(i1+1,i2,i3-1,0) + d200(i1-1,i2,i3-1,0)
                                    d103i = d002(i1+1,i2,i3+1,0) - d002(i1-1,i2,i3+1,0) - d002(i1+1,i2,i3-1,0) + d002(i1-1,i2,i3-1,0)
                                    d031i = d020(i1,i2+1,i3+1,0) - d020(i1,i2-1,i3+1,0) - d020(i1,i2+1,i3-1,0) + d020(i1,i2-1,i3-1,0)
                                    d013i = d002(i1,i2+1,i3+1,0) - d002(i1,i2-1,i3+1,0) - d002(i1,i2+1,i3-1,0) + d002(i1,i2-1,i3-1,0)
                  ! --- Laplacian to order 4 = lap2h + corrections 
                                    lap4h(i1,i2,i3,0) = lap2h(i1,i2,i3,0) + lapCoeff(i1,i2,i3,0)*crr1*d400(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*css1*d040(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*ctt1*d004(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*(cr1*d310i + cs1*d130i) + lapCoeff(i1,i2,i3,4)*(cr1*d301i + ct1*d103i) + lapCoeff(i1,i2,i3,5)*(cs1*d031i + ct1*d013i) + lapCoeff(i1,i2,i3,6)*cr1 *d300i + lapCoeff(i1,i2,i3,7)*cs1 *d030i + lapCoeff(i1,i2,i3,8)*ct1 *d003i 
                  ! --- Laplacian squared to order 2:
                                    lap2h200(i1,i2,i3,0) = lap2h(i1+1,i2,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1-1,i2,i3,0)
                                    lap2h020(i1,i2,i3,0) = lap2h(i1,i2+1,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2-1,i3,0)
                                    lap2h002(i1,i2,i3,0) = lap2h(i1,i2,i3+1,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2,i3-1,0)
                                    lap2h100i = lap2h(i1+1,i2,i3,0) - lap2h(i1-1,i2,i3,0)
                                    lap2h010i = lap2h(i1,i2+1,i3,0) - lap2h(i1,i2-1,i3,0)
                                    lap2h110i = lap2h(i1+1,i2+1,i3,0) - lap2h(i1-1,i2+1,i3,0) - lap2h(i1+1,i2-1,i3,0) + lap2h(i1-1,i2-1,i3,0)
                                    lap2h001i = lap2h(i1,i2,i3+1,0) - lap2h(i1,i2,i3-1,0)
                                    lap2h101i = lap2h(i1+1,i2,i3+1,0) - lap2h(i1-1,i2,i3+1,0) - lap2h(i1+1,i2,i3-1,0) + lap2h(i1-1,i2,i3-1,0)
                                    lap2h011i = lap2h(i1,i2+1,i3+1,0) - lap2h(i1,i2-1,i3+1,0) - lap2h(i1,i2+1,i3-1,0) + lap2h(i1,i2-1,i3-1,0)
                                    lap2hSq(i1,i2,i3,0) =  lapCoeff(i1,i2,i3,0)*lap2h200(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*lap2h020(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*lap2h002(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*lap2h110i  + lapCoeff(i1,i2,i3,4)*lap2h101i  + lapCoeff(i1,i2,i3,5)*lap2h011i  + lapCoeff(i1,i2,i3,6)*lap2h100i  + lapCoeff(i1,i2,i3,7)*lap2h010i  + lapCoeff(i1,i2,i3,8)*lap2h001i    
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
                                    d600(i1,i2,i3,0) = d400(i1+1,i2,i3,0) - 2*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0)
                                    d060(i1,i2,i3,0) = d040(i1,i2+1,i3,0) - 2*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0)
                                    d006(i1,i2,i3,0) = d004(i1,i2,i3+1,0) - 2*d004(i1,i2,i3,0) + d004(i1,i2,i3-1,0)
                                    d420(i1,i2,i3,0) = d220(i1+1,i2,i3,0) - 2*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0)
                                    d240(i1,i2,i3,0) = d040(i1+1,i2,i3,0) - 2*d040(i1,i2,i3,0) + d040(i1-1,i2,i3,0)
                                    d402(i1,i2,i3,0) = d202(i1+1,i2,i3,0) - 2*d202(i1,i2,i3,0) + d202(i1-1,i2,i3,0)
                                    d204(i1,i2,i3,0) = d004(i1+1,i2,i3,0) - 2*d004(i1,i2,i3,0) + d004(i1-1,i2,i3,0)
                                    d042(i1,i2,i3,0) = d022(i1,i2+1,i3,0) - 2*d022(i1,i2,i3,0) + d022(i1,i2-1,i3,0)
                                    d024(i1,i2,i3,0) = d004(i1,i2+1,i3,0) - 2*d004(i1,i2,i3,0) + d004(i1,i2-1,i3,0)
                                    d500i = d400(i1+1,i2,i3,0) - d400(i1-1,i2,i3,0)
                                    d050i = d040(i1,i2+1,i3,0) - d040(i1,i2-1,i3,0)
                                    d005i = d004(i1,i2,i3+1,0) - d004(i1,i2,i3-1,0)
                                    d510i = d400(i1+1,i2+1,i3,0) - d400(i1-1,i2+1,i3,0) - d400(i1+1,i2-1,i3,0) + d400(i1-1,i2-1,i3,0)
                                    d150i = d040(i1+1,i2+1,i3,0) - d040(i1-1,i2+1,i3,0) - d040(i1+1,i2-1,i3,0) + d040(i1-1,i2-1,i3,0)
                                    d330i = d220(i1+1,i2+1,i3,0) - d220(i1-1,i2+1,i3,0) - d220(i1+1,i2-1,i3,0) + d220(i1-1,i2-1,i3,0)
                                    d501i = d400(i1+1,i2,i3+1,0) - d400(i1-1,i2,i3+1,0) - d400(i1+1,i2,i3-1,0) + d400(i1-1,i2,i3-1,0)
                                    d105i = d004(i1+1,i2,i3+1,0) - d004(i1-1,i2,i3+1,0) - d004(i1+1,i2,i3-1,0) + d004(i1-1,i2,i3-1,0)
                                    d051i = d040(i1,i2+1,i3+1,0) - d040(i1,i2-1,i3+1,0) - d040(i1,i2+1,i3-1,0) + d040(i1,i2-1,i3-1,0)
                                    d015i = d004(i1,i2+1,i3+1,0) - d004(i1,i2-1,i3+1,0) - d004(i1,i2+1,i3-1,0) + d004(i1,i2-1,i3-1,0)
                                    d303i = d202(i1+1,i2,i3+1,0) - d202(i1-1,i2,i3+1,0) - d202(i1+1,i2,i3-1,0) + d202(i1-1,i2,i3-1,0)
                                    d033i = d022(i1,i2+1,i3+1,0) - d022(i1,i2-1,i3+1,0) - d022(i1,i2+1,i3-1,0) + d022(i1,i2-1,i3-1,0)
                  ! --- Laplacian to order 6 = lap4h + corrections 
                                    lap6h(i1,i2,i3,0) = lap4h(i1,i2,i3,0) + lapCoeff(i1,i2,i3,0)*crr2*d600(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*css2*d060(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*ctt2*d006(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*(cr2*d510i + cs2*d150i + cr1*cs1*d330i ) + lapCoeff(i1,i2,i3,4)*(cr2*d501i + ct2*d105i + cr1*ct1*d303i ) + lapCoeff(i1,i2,i3,5)*(cs2*d051i + ct2*d015i + cs1*ct1*d033i ) + lapCoeff(i1,i2,i3,6)*cr2 *d500i + lapCoeff(i1,i2,i3,7)*cs2 *d050i + lapCoeff(i1,i2,i3,8)*ct2 *d005i 
                                    lap2hSq200(i1,i2,i3,0) = lap2hSq(i1+1,i2,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1-1,i2,i3,0)
                                    lap2hSq020(i1,i2,i3,0) = lap2hSq(i1,i2+1,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2-1,i3,0)
                                    lap2hSq002(i1,i2,i3,0) = lap2hSq(i1,i2,i3+1,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2,i3-1,0)
                                    lap2hSq100i = lap2hSq(i1+1,i2,i3,0) - lap2hSq(i1-1,i2,i3,0)
                                    lap2hSq010i = lap2hSq(i1,i2+1,i3,0) - lap2hSq(i1,i2-1,i3,0)
                                    lap2hSq001i = lap2hSq(i1,i2,i3+1,0) - lap2hSq(i1,i2,i3-1,0)
                                    lap2hSq110i = lap2hSq(i1+1,i2+1,i3,0) - lap2hSq(i1-1,i2+1,i3,0) - lap2hSq(i1+1,i2-1,i3,0) + lap2hSq(i1-1,i2-1,i3,0)
                                    lap2hSq101i = lap2hSq(i1+1,i2,i3+1,0) - lap2hSq(i1-1,i2,i3+1,0) - lap2hSq(i1+1,i2,i3-1,0) + lap2hSq(i1-1,i2,i3-1,0)
                                    lap2hSq011i = lap2hSq(i1,i2+1,i3+1,0) - lap2hSq(i1,i2-1,i3+1,0) - lap2hSq(i1,i2+1,i3-1,0) + lap2hSq(i1,i2-1,i3-1,0)
                                    lap2hCubed(i1,i2,i3,0) =  + lapCoeff(i1,i2,i3,0)*lap2hSq200(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*lap2hSq020(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*lap2hSq002(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*lap2hSq110i  + lapCoeff(i1,i2,i3,4)*lap2hSq101i  + lapCoeff(i1,i2,i3,5)*lap2hSq011i  + lapCoeff(i1,i2,i3,6)*lap2hSq100i  + lapCoeff(i1,i2,i3,7)*lap2hSq010i  + lapCoeff(i1,i2,i3,8)*lap2hSq001i   
                                  end do
                                  end do
                                  end do
               ! --- SPLIT LOOPS
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                  ! --- Laplacian squared to order 4 = 
                  !  lap2h*( lap4h ) + corrections*( Lap2h )
                                    lap4h200(i1,i2,i3,0) = lap4h(i1+1,i2,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1-1,i2,i3,0)
                                    lap4h020(i1,i2,i3,0) = lap4h(i1,i2+1,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2-1,i3,0)
                                    lap4h002(i1,i2,i3,0) = lap4h(i1,i2,i3+1,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2,i3-1,0)
                                    lap2h400(i1,i2,i3,0) = lap2h200(i1+1,i2,i3,0) - 2*lap2h200(i1,i2,i3,0) + lap2h200(i1-1,i2,i3,0)
                                    lap2h040(i1,i2,i3,0) = lap2h020(i1,i2+1,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1,i2-1,i3,0)
                                    lap2h004(i1,i2,i3,0) = lap2h002(i1,i2,i3+1,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1,i2,i3-1,0)
                                    lap2h220(i1,i2,i3,0) = lap2h020(i1+1,i2,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1-1,i2,i3,0)
                                    lap2h202(i1,i2,i3,0) = lap2h002(i1+1,i2,i3,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1-1,i2,i3,0)
                                    lap2h022(i1,i2,i3,0) = lap2h002(i1,i2+1,i3,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1,i2-1,i3,0)
                                    lap4h100i = lap4h(i1+1,i2,i3,0) - lap4h(i1-1,i2,i3,0)
                                    lap4h010i = lap4h(i1,i2+1,i3,0) - lap4h(i1,i2-1,i3,0)
                                    lap4h001i = lap4h(i1,i2,i3+1,0) - lap4h(i1,i2,i3-1,0)
                                    lap4h110i = lap4h(i1+1,i2+1,i3,0) - lap4h(i1-1,i2+1,i3,0) - lap4h(i1+1,i2-1,i3,0) + lap4h(i1-1,i2-1,i3,0)
                                    lap4h101i = lap4h(i1+1,i2,i3+1,0) - lap4h(i1-1,i2,i3+1,0) - lap4h(i1+1,i2,i3-1,0) + lap4h(i1-1,i2,i3-1,0)
                                    lap4h011i = lap4h(i1,i2+1,i3+1,0) - lap4h(i1,i2-1,i3+1,0) - lap4h(i1,i2+1,i3-1,0) + lap4h(i1,i2-1,i3-1,0)
                                    lap2h300i = lap2h200(i1+1,i2,i3,0) - lap2h200(i1-1,i2,i3,0)
                                    lap2h030i = lap2h020(i1,i2+1,i3,0) - lap2h020(i1,i2-1,i3,0)
                                    lap2h003i = lap2h002(i1,i2,i3+1,0) - lap2h002(i1,i2,i3-1,0)
                                    lap2h310i = lap2h200(i1+1,i2+1,i3,0) - lap2h200(i1-1,i2+1,i3,0) - lap2h200(i1+1,i2-1,i3,0) + lap2h200(i1-1,i2-1,i3,0)
                                    lap2h130i = lap2h020(i1+1,i2+1,i3,0) - lap2h020(i1-1,i2+1,i3,0) - lap2h020(i1+1,i2-1,i3,0) + lap2h020(i1-1,i2-1,i3,0)
                                    lap2h301i = lap2h200(i1+1,i2,i3+1,0) - lap2h200(i1-1,i2,i3+1,0) - lap2h200(i1+1,i2,i3-1,0) + lap2h200(i1-1,i2,i3-1,0)
                                    lap2h103i = lap2h002(i1+1,i2,i3+1,0) - lap2h002(i1-1,i2,i3+1,0) - lap2h002(i1+1,i2,i3-1,0) + lap2h002(i1-1,i2,i3-1,0)
                                    lap2h031i = lap2h020(i1,i2+1,i3+1,0) - lap2h020(i1,i2-1,i3+1,0) - lap2h020(i1,i2+1,i3-1,0) + lap2h020(i1,i2-1,i3-1,0)
                                    lap2h013i = lap2h002(i1,i2+1,i3+1,0) - lap2h002(i1,i2-1,i3+1,0) - lap2h002(i1,i2+1,i3-1,0) + lap2h002(i1,i2-1,i3-1,0)
                                    lap4hSq(i1,i2,i3,0) =     lapCoeff(i1,i2,i3,0)*( lap4h200(i1,i2,i3,0) + crr1*lap2h400(i1,i2,i3,0) )    + lapCoeff(i1,i2,i3,1)*( lap4h020(i1,i2,i3,0) + css1*lap2h040(i1,i2,i3,0) )     + lapCoeff(i1,i2,i3,2)*( lap4h002(i1,i2,i3,0) + ctt1*lap2h004(i1,i2,i3,0) )     + lapCoeff(i1,i2,i3,3)*( lap4h110i + cr1*lap2h310i + cs1*lap2h130i ) + lapCoeff(i1,i2,i3,4)*( lap4h101i + cr1*lap2h301i + ct1*lap2h103i ) + lapCoeff(i1,i2,i3,5)*( lap4h011i + cs1*lap2h031i + ct1*lap2h013i ) + lapCoeff(i1,i2,i3,6)*( lap4h100i + cr1 *lap2h300i )    + lapCoeff(i1,i2,i3,7)*( lap4h010i + cs1 *lap2h030i )    + lapCoeff(i1,i2,i3,8)*( lap4h001i + ct1 *lap2h003i )      
                                  end do
                                  end do
                                  end do
              ! ===========  FINAL LOOP TO FILL IN THE SOLUTION ============
                            numGhost1=0;
                            n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                            n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                            n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                    d800i = d600(i1+1,i2,i3,0) - 2*d600(i1,i2,i3,0) + d600(i1-1,i2,i3,0)
                                    d080i = d060(i1,i2+1,i3,0) - 2*d060(i1,i2,i3,0) + d060(i1,i2-1,i3,0)
                                    d008i = d006(i1,i2,i3+1,0) - 2*d006(i1,i2,i3,0) + d006(i1,i2,i3-1,0)
                                    d700i = d600(i1+1,i2,i3,0) - d600(i1-1,i2,i3,0)
                                    d070i = d060(i1,i2+1,i3,0) - d060(i1,i2-1,i3,0)
                                    d007i = d006(i1,i2,i3+1,0) - d006(i1,i2,i3-1,0)
                                    d710i = d600(i1+1,i2+1,i3,0) - d600(i1-1,i2+1,i3,0) - d600(i1+1,i2-1,i3,0) + d600(i1-1,i2-1,i3,0)
                                    d170i = d060(i1+1,i2+1,i3,0) - d060(i1-1,i2+1,i3,0) - d060(i1+1,i2-1,i3,0) + d060(i1-1,i2-1,i3,0)
                                    d701i = d600(i1+1,i2,i3+1,0) - d600(i1-1,i2,i3+1,0) - d600(i1+1,i2,i3-1,0) + d600(i1-1,i2,i3-1,0)
                                    d107i = d006(i1+1,i2,i3+1,0) - d006(i1-1,i2,i3+1,0) - d006(i1+1,i2,i3-1,0) + d006(i1-1,i2,i3-1,0)
                                    d071i = d060(i1,i2+1,i3+1,0) - d060(i1,i2-1,i3+1,0) - d060(i1,i2+1,i3-1,0) + d060(i1,i2-1,i3-1,0)
                                    d017i = d006(i1,i2+1,i3+1,0) - d006(i1,i2-1,i3+1,0) - d006(i1,i2+1,i3-1,0) + d006(i1,i2-1,i3-1,0)
                                    d530i = d420(i1+1,i2+1,i3,0) - d420(i1-1,i2+1,i3,0) - d420(i1+1,i2-1,i3,0) + d420(i1-1,i2-1,i3,0)
                                    d350i = d240(i1+1,i2+1,i3,0) - d240(i1-1,i2+1,i3,0) - d240(i1+1,i2-1,i3,0) + d240(i1-1,i2-1,i3,0)
                                    d503i = d402(i1+1,i2,i3+1,0) - d402(i1-1,i2,i3+1,0) - d402(i1+1,i2,i3-1,0) + d402(i1-1,i2,i3-1,0)
                                    d305i = d204(i1+1,i2,i3+1,0) - d204(i1-1,i2,i3+1,0) - d204(i1+1,i2,i3-1,0) + d204(i1-1,i2,i3-1,0)
                                    d053i = d042(i1,i2+1,i3+1,0) - d042(i1,i2-1,i3+1,0) - d042(i1,i2+1,i3-1,0) + d042(i1,i2-1,i3-1,0)
                                    d035i = d024(i1,i2+1,i3+1,0) - d024(i1,i2-1,i3+1,0) - d024(i1,i2+1,i3-1,0) + d024(i1,i2-1,i3-1,0)
                  ! --- Laplacian to order 8 = lap6h + corrections 
                                    lap8h = lap6h(i1,i2,i3,0)                                                         + lapCoeff(i1,i2,i3,0)*crr3*d800i                                               + lapCoeff(i1,i2,i3,1)*css3*d080i                                               + lapCoeff(i1,i2,i3,2)*ctt3*d008i                                               + lapCoeff(i1,i2,i3,3)*(cr3*d710i + cs3*d170i + cr2*cs1*d530i + cr1*cs2*d350i ) + lapCoeff(i1,i2,i3,4)*(cr3*d701i + ct3*d107i + cr2*ct1*d503i + cr1*ct2*d305i ) + lapCoeff(i1,i2,i3,5)*(cs3*d071i + ct3*d017i + cs2*ct1*d053i + cs1*ct2*d035i ) + lapCoeff(i1,i2,i3,6)* cr3*d700i                                               + lapCoeff(i1,i2,i3,7)* cs3*d070i                                               + lapCoeff(i1,i2,i3,8)* ct3*d007i 
                  ! --- Laplacian^4 4p (4th power) order 2: 
                                    lap2hCubed200i = lap2hCubed(i1+1,i2,i3,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1-1,i2,i3,0)
                                    lap2hCubed020i = lap2hCubed(i1,i2+1,i3,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1,i2-1,i3,0)
                                    lap2hCubed002i = lap2hCubed(i1,i2,i3+1,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1,i2,i3-1,0)
                                    lap2hCubed100i = lap2hCubed(i1+1,i2,i3,0) - lap2hCubed(i1-1,i2,i3,0)
                                    lap2hCubed010i = lap2hCubed(i1,i2+1,i3,0) - lap2hCubed(i1,i2-1,i3,0)
                                    lap2hCubed001i = lap2hCubed(i1,i2,i3+1,0) - lap2hCubed(i1,i2,i3-1,0)
                                    lap2hCubed110i = lap2hCubed(i1+1,i2+1,i3,0) - lap2hCubed(i1-1,i2+1,i3,0) - lap2hCubed(i1+1,i2-1,i3,0) + lap2hCubed(i1-1,i2-1,i3,0)
                                    lap2hCubed101i = lap2hCubed(i1+1,i2,i3+1,0) - lap2hCubed(i1-1,i2,i3+1,0) - lap2hCubed(i1+1,i2,i3-1,0) + lap2hCubed(i1-1,i2,i3-1,0)
                                    lap2hCubed011i = lap2hCubed(i1,i2+1,i3+1,0) - lap2hCubed(i1,i2-1,i3+1,0) - lap2hCubed(i1,i2+1,i3-1,0) + lap2hCubed(i1,i2-1,i3-1,0)
                                    lap2h4p  =                             + lapCoeff(i1,i2,i3,0)*lap2hCubed200i  + lapCoeff(i1,i2,i3,1)*lap2hCubed020i  + lapCoeff(i1,i2,i3,2)*lap2hCubed002i  + lapCoeff(i1,i2,i3,3)*lap2hCubed110i  + lapCoeff(i1,i2,i3,4)*lap2hCubed101i  + lapCoeff(i1,i2,i3,5)*lap2hCubed011i  + lapCoeff(i1,i2,i3,6)*lap2hCubed100i  + lapCoeff(i1,i2,i3,7)*lap2hCubed010i  + lapCoeff(i1,i2,i3,8)*lap2hCubed001i    
                  ! --- Laplacian squared to order 6 :
                  !   Lap6h = Lap4h + M4  = (Lap2h) + M2 + M4 
                  !   Lap6h*Lap6h = [ (Lap2h) + M2 + M4 ] [ (Lap2h) + M2 + M4 ]
                  !               = Lap2h*Lap6h + M2*Lap4h + M4*Lap2h + O(h^6)
                                    lap6h200i = lap6h(i1+1,i2,i3,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1-1,i2,i3,0)
                                    lap6h020i = lap6h(i1,i2+1,i3,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1,i2-1,i3,0)
                                    lap6h002i = lap6h(i1,i2,i3+1,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1,i2,i3-1,0)
                                    lap6h100i = lap6h(i1+1,i2,i3,0) - lap6h(i1-1,i2,i3,0)
                                    lap6h010i = lap6h(i1,i2+1,i3,0) - lap6h(i1,i2-1,i3,0)
                                    lap6h001i = lap6h(i1,i2,i3+1,0) - lap6h(i1,i2,i3-1,0)
                                    lap6h110i = lap6h(i1+1,i2+1,i3,0) - lap6h(i1-1,i2+1,i3,0) - lap6h(i1+1,i2-1,i3,0) + lap6h(i1-1,i2-1,i3,0)
                                    lap6h101i = lap6h(i1+1,i2,i3+1,0) - lap6h(i1-1,i2,i3+1,0) - lap6h(i1+1,i2,i3-1,0) + lap6h(i1-1,i2,i3-1,0)
                                    lap6h011i = lap6h(i1,i2+1,i3+1,0) - lap6h(i1,i2-1,i3+1,0) - lap6h(i1,i2+1,i3-1,0) + lap6h(i1,i2-1,i3-1,0)
                                    lap4h400i = lap4h200(i1+1,i2,i3,0) - 2*lap4h200(i1,i2,i3,0) + lap4h200(i1-1,i2,i3,0)
                                    lap4h040i = lap4h020(i1,i2+1,i3,0) - 2*lap4h020(i1,i2,i3,0) + lap4h020(i1,i2-1,i3,0)
                                    lap4h004i = lap4h002(i1,i2,i3+1,0) - 2*lap4h002(i1,i2,i3,0) + lap4h002(i1,i2,i3-1,0)
                                    lap4h300i = lap4h200(i1+1,i2,i3,0) - lap4h200(i1-1,i2,i3,0)
                                    lap4h030i = lap4h020(i1,i2+1,i3,0) - lap4h020(i1,i2-1,i3,0)
                                    lap4h003i = lap4h002(i1,i2,i3+1,0) - lap4h002(i1,i2,i3-1,0)
                                    lap4h310i = lap4h200(i1+1,i2+1,i3,0) - lap4h200(i1-1,i2+1,i3,0) - lap4h200(i1+1,i2-1,i3,0) + lap4h200(i1-1,i2-1,i3,0)
                                    lap4h130i = lap4h020(i1+1,i2+1,i3,0) - lap4h020(i1-1,i2+1,i3,0) - lap4h020(i1+1,i2-1,i3,0) + lap4h020(i1-1,i2-1,i3,0)
                                    lap4h301i = lap4h200(i1+1,i2,i3+1,0) - lap4h200(i1-1,i2,i3+1,0) - lap4h200(i1+1,i2,i3-1,0) + lap4h200(i1-1,i2,i3-1,0)
                                    lap4h103i = lap4h002(i1+1,i2,i3+1,0) - lap4h002(i1-1,i2,i3+1,0) - lap4h002(i1+1,i2,i3-1,0) + lap4h002(i1-1,i2,i3-1,0)
                                    lap4h031i = lap4h020(i1,i2+1,i3+1,0) - lap4h020(i1,i2-1,i3+1,0) - lap4h020(i1,i2+1,i3-1,0) + lap4h020(i1,i2-1,i3-1,0)
                                    lap4h013i = lap4h002(i1,i2+1,i3+1,0) - lap4h002(i1,i2-1,i3+1,0) - lap4h002(i1,i2+1,i3-1,0) + lap4h002(i1,i2-1,i3-1,0)
                                    lap2h600i = lap2h400(i1+1,i2,i3,0) - 2*lap2h400(i1,i2,i3,0) + lap2h400(i1-1,i2,i3,0)
                                    lap2h060i = lap2h040(i1,i2+1,i3,0) - 2*lap2h040(i1,i2,i3,0) + lap2h040(i1,i2-1,i3,0)
                                    lap2h006i = lap2h004(i1,i2,i3+1,0) - 2*lap2h004(i1,i2,i3,0) + lap2h004(i1,i2,i3-1,0)
                                    lap2h500i = lap2h400(i1+1,i2,i3,0) - lap2h400(i1-1,i2,i3,0)
                                    lap2h050i = lap2h040(i1,i2+1,i3,0) - lap2h040(i1,i2-1,i3,0)
                                    lap2h005i = lap2h004(i1,i2,i3+1,0) - lap2h004(i1,i2,i3-1,0)
                                    lap2h510i = lap2h400(i1+1,i2+1,i3,0) - lap2h400(i1-1,i2+1,i3,0) - lap2h400(i1+1,i2-1,i3,0) + lap2h400(i1-1,i2-1,i3,0)
                                    lap2h150i = lap2h040(i1+1,i2+1,i3,0) - lap2h040(i1-1,i2+1,i3,0) - lap2h040(i1+1,i2-1,i3,0) + lap2h040(i1-1,i2-1,i3,0)
                                    lap2h330i = lap2h220(i1+1,i2+1,i3,0) - lap2h220(i1-1,i2+1,i3,0) - lap2h220(i1+1,i2-1,i3,0) + lap2h220(i1-1,i2-1,i3,0)
                                    lap2h501i = lap2h400(i1+1,i2,i3+1,0) - lap2h400(i1-1,i2,i3+1,0) - lap2h400(i1+1,i2,i3-1,0) + lap2h400(i1-1,i2,i3-1,0)
                                    lap2h105i = lap2h004(i1+1,i2,i3+1,0) - lap2h004(i1-1,i2,i3+1,0) - lap2h004(i1+1,i2,i3-1,0) + lap2h004(i1-1,i2,i3-1,0)
                                    lap2h051i = lap2h040(i1,i2+1,i3+1,0) - lap2h040(i1,i2-1,i3+1,0) - lap2h040(i1,i2+1,i3-1,0) + lap2h040(i1,i2-1,i3-1,0)
                                    lap2h015i = lap2h004(i1,i2+1,i3+1,0) - lap2h004(i1,i2-1,i3+1,0) - lap2h004(i1,i2+1,i3-1,0) + lap2h004(i1,i2-1,i3-1,0)
                                    lap2h303i = lap2h202(i1+1,i2,i3+1,0) - lap2h202(i1-1,i2,i3+1,0) - lap2h202(i1+1,i2,i3-1,0) + lap2h202(i1-1,i2,i3-1,0)
                                    lap2h033i = lap2h022(i1,i2+1,i3+1,0) - lap2h022(i1,i2-1,i3+1,0) - lap2h022(i1,i2+1,i3-1,0) + lap2h022(i1,i2-1,i3-1,0)
                                    lap6hSq =                                                                                     lapCoeff(i1,i2,i3,0)*(lap6h200i + crr1*lap4h400i + crr2*lap2h600i )                       + lapCoeff(i1,i2,i3,1)*(lap6h020i + css1*lap4h040i + css2*lap2h060i )                       + lapCoeff(i1,i2,i3,2)*(lap6h002i + ctt1*lap4h004i + ctt2*lap2h006i )                       + lapCoeff(i1,i2,i3,3)*(lap6h110i +  cr1*lap4h310i +  cr2*lap2h510i                         +  cs1*lap4h130i +  cs2*lap2h150i + cr1*cs1*lap2h330i )   + lapCoeff(i1,i2,i3,4)*(lap6h101i +  cr1*lap4h301i +  cr2*lap2h501i                         +  ct1*lap4h103i +  ct2*lap2h105i + cr1*ct1*lap2h303i )   + lapCoeff(i1,i2,i3,5)*(lap6h011i +  cs1*lap4h031i +  cs2*lap2h051i                         +  ct1*lap4h013i +  ct2*lap2h015i + cs1*ct1*lap2h033i )   + lapCoeff(i1,i2,i3,6)*(lap6h100i +  cr1*lap4h300i +  cr2*lap2h500i )                       + lapCoeff(i1,i2,i3,7)*(lap6h010i +  cs1*lap4h030i +  cs2*lap2h050i )                       + lapCoeff(i1,i2,i3,8)*(lap6h010i +  ct1*lap4h003i +  ct2*lap2h005i )                         
                  ! --- Laplacian CUBED to order 4 
                                    lap4hSq200i = lap4hSq(i1+1,i2,i3,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1-1,i2,i3,0)
                                    lap4hSq020i = lap4hSq(i1,i2+1,i3,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1,i2-1,i3,0)
                                    lap4hSq002i = lap4hSq(i1,i2,i3+1,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1,i2,i3-1,0)
                                    lap4hSq100i = lap4hSq(i1+1,i2,i3,0) - lap4hSq(i1-1,i2,i3,0)
                                    lap4hSq010i = lap4hSq(i1,i2+1,i3,0) - lap4hSq(i1,i2-1,i3,0)
                                    lap4hSq001i = lap4hSq(i1,i2,i3+1,0) - lap4hSq(i1,i2,i3-1,0)
                                    lap4hSq110i = lap4hSq(i1+1,i2+1,i3,0) - lap4hSq(i1-1,i2+1,i3,0) - lap4hSq(i1+1,i2-1,i3,0) + lap4hSq(i1-1,i2-1,i3,0)
                                    lap4hSq101i = lap4hSq(i1+1,i2,i3+1,0) - lap4hSq(i1-1,i2,i3+1,0) - lap4hSq(i1+1,i2,i3-1,0) + lap4hSq(i1-1,i2,i3-1,0)
                                    lap4hSq011i = lap4hSq(i1,i2+1,i3+1,0) - lap4hSq(i1,i2-1,i3+1,0) - lap4hSq(i1,i2+1,i3-1,0) + lap4hSq(i1,i2-1,i3-1,0)
                                    lap2hSq400i = lap2hSq200(i1+1,i2,i3,0) - 2*lap2hSq200(i1,i2,i3,0) + lap2hSq200(i1-1,i2,i3,0)
                                    lap2hSq040i = lap2hSq020(i1,i2+1,i3,0) - 2*lap2hSq020(i1,i2,i3,0) + lap2hSq020(i1,i2-1,i3,0)
                                    lap2hSq004i = lap2hSq002(i1,i2,i3+1,0) - 2*lap2hSq002(i1,i2,i3,0) + lap2hSq002(i1,i2,i3-1,0)
                                    lap2hSq300i = lap2hSq200(i1+1,i2,i3,0) - lap2hSq200(i1-1,i2,i3,0)
                                    lap2hSq030i = lap2hSq020(i1,i2+1,i3,0) - lap2hSq020(i1,i2-1,i3,0)
                                    lap2hSq003i = lap2hSq002(i1,i2,i3+1,0) - lap2hSq002(i1,i2,i3-1,0)
                                    lap2hSq310i = lap2hSq200(i1+1,i2+1,i3,0) - lap2hSq200(i1-1,i2+1,i3,0) - lap2hSq200(i1+1,i2-1,i3,0) + lap2hSq200(i1-1,i2-1,i3,0)
                                    lap2hSq130i = lap2hSq020(i1+1,i2+1,i3,0) - lap2hSq020(i1-1,i2+1,i3,0) - lap2hSq020(i1+1,i2-1,i3,0) + lap2hSq020(i1-1,i2-1,i3,0)
                                    lap2hSq301i = lap2hSq200(i1+1,i2,i3+1,0) - lap2hSq200(i1-1,i2,i3+1,0) - lap2hSq200(i1+1,i2,i3-1,0) + lap2hSq200(i1-1,i2,i3-1,0)
                                    lap2hSq103i = lap2hSq002(i1+1,i2,i3+1,0) - lap2hSq002(i1-1,i2,i3+1,0) - lap2hSq002(i1+1,i2,i3-1,0) + lap2hSq002(i1-1,i2,i3-1,0)
                                    lap2hSq031i = lap2hSq020(i1,i2+1,i3+1,0) - lap2hSq020(i1,i2-1,i3+1,0) - lap2hSq020(i1,i2+1,i3-1,0) + lap2hSq020(i1,i2-1,i3-1,0)
                                    lap2hSq013i = lap2hSq002(i1,i2+1,i3+1,0) - lap2hSq002(i1,i2-1,i3+1,0) - lap2hSq002(i1,i2+1,i3-1,0) + lap2hSq002(i1,i2-1,i3-1,0)
                                    lap4hCubed =                                                    lapCoeff(i1,i2,i3,0)*(lap4hSq200i + crr1*lap2hSq400i )   + lapCoeff(i1,i2,i3,1)*(lap4hSq020i + css1*lap2hSq040i )   + lapCoeff(i1,i2,i3,2)*(lap4hSq002i + ctt1*lap2hSq004i )   + lapCoeff(i1,i2,i3,3)*(lap4hSq110i +  cr1*lap2hSq310i     +  cs1*lap2hSq130i )   + lapCoeff(i1,i2,i3,4)*(lap4hSq101i +  cr1*lap2hSq301i     +  ct1*lap2hSq103i )   + lapCoeff(i1,i2,i3,5)*(lap4hSq011i +  cs1*lap2hSq031i     +  ct1*lap2hSq013i )   + lapCoeff(i1,i2,i3,6)*(lap4hSq100i + cr1 *lap2hSq300i )   + lapCoeff(i1,i2,i3,7)*(lap4hSq010i + cs1 *lap2hSq030i )   + lapCoeff(i1,i2,i3,8)*(lap4hSq001i + ct1 *lap2hSq003i )     
                  ! --- Modified equation space-time update ----
                                    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m)  + cdtsq*( lap8h )               + cdtPow4By12*( lap6hSq )       + cdtPow6By360*( lap4hCubed )   + cdtPow8By20160*( lap2h4p )    +dtSq*fv(m)                    
                                  end do
                                  end do
                                  end do
                    else
              ! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
              ! Example: 
              ! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...
cr0 = 1.; cr1 = -1/6.; cr2 = 1/30.; cr3 = -1/140.; 
cs0 = 1.; cs1 = -1/6.; cs2 = 1/30.; cs3 = -1/140.; 
ct0 = 1.; ct1 = -1/6.; ct2 = 1/30.; ct3 = -1/140.; 
crr0 = 1.; crr1 = -1/12.; crr2 = 1/90.; crr3 = -1/560.; 
css0 = 1.; css1 = -1/12.; css2 = 1/90.; css3 = -1/560.; 
ctt0 = 1.; ctt1 = -1/12.; ctt2 = 1/90.; ctt3 = -1/560.; 
crrr0 = 1.; crrr1 = -1/4.; crrr2 = 7/120.; crrr3 = -41/3024.; 
csss0 = 1.; csss1 = -1/4.; csss2 = 7/120.; csss3 = -41/3024.; 
cttt0 = 1.; cttt1 = -1/4.; cttt2 = 7/120.; cttt3 = -41/3024.; 
crrrr0 = 1.; crrrr1 = -1/6.; crrrr2 = 7/240.; crrrr3 = -41/7560.; 
cssss0 = 1.; cssss1 = -1/6.; cssss2 = 7/240.; cssss3 = -41/7560.; 
ctttt0 = 1.; ctttt1 = -1/6.; ctttt2 = 7/240.; ctttt3 = -41/7560.; 
crrrrr0 = 1.; crrrrr1 = -1/3.; crrrrr2 = 13/144.; crrrrr3 = -139/6048.; 
csssss0 = 1.; csssss1 = -1/3.; csssss2 = 13/144.; csssss3 = -139/6048.; 
cttttt0 = 1.; cttttt1 = -1/3.; cttttt2 = 13/144.; cttttt3 = -139/6048.; 
crrrrrr0 = 1.; crrrrrr1 = -1/4.; crrrrrr2 = 13/240.; crrrrrr3 = -139/12096.; 
cssssss0 = 1.; cssssss1 = -1/4.; cssssss2 = 13/240.; cssssss3 = -139/12096.; 
ctttttt0 = 1.; ctttttt1 = -1/4.; ctttttt2 = 13/240.; ctttttt3 = -139/12096.; 
                            dr1=dr(0); dr1i=1./dr1;
                            dr2=dr(1); dr2i=1./dr2;
                            dr3=dr(2); dr3i=1./dr3;
                            fv(m)=0.
                            if( lapCoeff(0,0,0,0).le.0. )then
                ! --- Evaluate and store coefficients in Laplacian ---
                                write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'
                                numGhost1=3;
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
                                    if( (i1-4).ge.nd1a .and. (i1+4).le.nd1b )then
                                        diffOrder1=8
                                    elseif( (i1-3).ge.nd1a .and. (i1+3).le.nd1b )then
                                        diffOrder1=6
                                    elseif( (i1-2).ge.nd1a .and. (i1+2).le.nd1b )then
                                        diffOrder1=4
                                    elseif( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
                                        diffOrder1=2
                                    else
                                        stop 999
                                    end if
                                    if( (i2-4).ge.nd2a .and. (i2+4).le.nd2b )then
                                        diffOrder2=8
                                    elseif( (i2-3).ge.nd2a .and. (i2+3).le.nd2b )then
                                        diffOrder2=6
                                    elseif( (i2-2).ge.nd2a .and. (i2+2).le.nd2b )then
                                        diffOrder2=4
                                    elseif( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
                                        diffOrder2=2
                                    else
                                        stop 999
                                    end if
                                    if( (i3-4).ge.nd3a .and. (i3+4).le.nd3b )then
                                        diffOrder3=8
                                    elseif( (i3-3).ge.nd3a .and. (i3+3).le.nd3b )then
                                        diffOrder3=6
                                    elseif( (i3-2).ge.nd3a .and. (i3+2).le.nd3b )then
                                        diffOrder3=4
                                    elseif( (i3-1).ge.nd3a .and. (i3+1).le.nd3b )then
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
                                        rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
                                        rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
                                        rzs = ( 8*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) )*(dr2i/12.) 
                                        sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
                                        sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
                                        szs = ( 8*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) )*(dr2i/12.) 
                                        txs = ( 8*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) )*(dr2i/12.) 
                                        tys = ( 8*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) )*(dr2i/12.) 
                                        tzs = ( 8*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) )*(dr2i/12.) 
                                    elseif( diffOrder2.eq.8 )then
                                        rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
                                        rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
                                        rzs = ( 8*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) )*(dr2i/12.) 
                                        sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
                                        sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
                                        szs = ( 8*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) )*(dr2i/12.) 
                                        txs = ( 8*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) )*(dr2i/12.) 
                                        tys = ( 8*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) )*(dr2i/12.) 
                                        tzs = ( 8*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) )*(dr2i/12.) 
                                    end if
                                    if( diffOrder3.eq.2 )then
                                        rxt = (rsxy(i1,i2,i3+1,0,0)-rsxy(i1,i2,i3-1,0,0))*(.5*dr2i) 
                                        ryt = (rsxy(i1,i2,i3+1,0,1)-rsxy(i1,i2,i3-1,0,1))*(.5*dr2i) 
                                        rzt = (rsxy(i1,i2,i3+1,0,2)-rsxy(i1,i2,i3-1,0,2))*(.5*dr2i) 
                                        sxt = (rsxy(i1,i2,i3+1,1,0)-rsxy(i1,i2,i3-1,1,0))*(.5*dr2i) 
                                        syt = (rsxy(i1,i2,i3+1,1,1)-rsxy(i1,i2,i3-1,1,1))*(.5*dr2i) 
                                        szt = (rsxy(i1,i2,i3+1,1,2)-rsxy(i1,i2,i3-1,1,2))*(.5*dr2i) 
                                        txt = (rsxy(i1,i2,i3+1,2,0)-rsxy(i1,i2,i3-1,2,0))*(.5*dr2i) 
                                        tyt = (rsxy(i1,i2,i3+1,2,1)-rsxy(i1,i2,i3-1,2,1))*(.5*dr2i) 
                                        tzt = (rsxy(i1,i2,i3+1,2,2)-rsxy(i1,i2,i3-1,2,2))*(.5*dr2i) 
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
                                    lapCoeff(i1,i2,i3,3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i*.25
                                    lapCoeff(i1,i2,i3,4) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr2i*.25
                                    lapCoeff(i1,i2,i3,5) = 2.*(sx*tx + sy*ty + sz*tz )*dr1i*dr2i*.25
                                    lapCoeff(i1,i2,i3,6) = (rxx + ryy + rzz)*dr1i*.5
                                    lapCoeff(i1,i2,i3,7) = (sxx + syy + tyy)*dr2i*.5 
                                    lapCoeff(i1,i2,i3,8) = (txx + tyy + tzz)*dr3i*.5 
                                    end if ! mask .ne. 0
                                  end do
                                  end do
                                  end do
                            end if ! end assignLapCoeff
                            numGhost1=3;
                            n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                            n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                            n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                if( mask(i1,i2,i3).ne.0 )then
                                    d200(i1,i2,i3,0) = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
                                    d020(i1,i2,i3,0) = u(i1,i2+1,i3,0) - 2*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
                                    d002(i1,i2,i3,0) = u(i1,i2,i3+1,0) - 2*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
                                    d100i = u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0)
                                    d010i = u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0)
                                    d110i = u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0)
                                    d001i = u(i1,i2,i3+1,0) - u(i1,i2,i3-1,0)
                                    d101i = u(i1+1,i2,i3+1,0) - u(i1-1,i2,i3+1,0) - u(i1+1,i2,i3-1,0) + u(i1-1,i2,i3-1,0)
                                    d011i = u(i1,i2+1,i3+1,0) - u(i1,i2-1,i3+1,0) - u(i1,i2+1,i3-1,0) + u(i1,i2-1,i3-1,0)
                                    lap2h(i1,i2,i3,0) = lapCoeff(i1,i2,i3,0)*d200(i1,i2,i3,0) +lapCoeff(i1,i2,i3,1)*d020(i1,i2,i3,0) +lapCoeff(i1,i2,i3,2)*d002(i1,i2,i3,0) +lapCoeff(i1,i2,i3,3)*d110i + lapCoeff(i1,i2,i3,4)*d101i + lapCoeff(i1,i2,i3,5)*d011i + lapCoeff(i1,i2,i3,6)*d100i + lapCoeff(i1,i2,i3,7)*d010i + lapCoeff(i1,i2,i3,8)*d001i
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
                                    d400(i1,i2,i3,0) = d200(i1+1,i2,i3,0) - 2*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0)
                                    d040(i1,i2,i3,0) = d020(i1,i2+1,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1,i2-1,i3,0)
                                    d004(i1,i2,i3,0) = d002(i1,i2,i3+1,0) - 2*d002(i1,i2,i3,0) + d002(i1,i2,i3-1,0)
                                    d220(i1,i2,i3,0) = d020(i1+1,i2,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1-1,i2,i3,0)
                                    d202(i1,i2,i3,0) = d002(i1+1,i2,i3,0) - 2*d002(i1,i2,i3,0) + d002(i1-1,i2,i3,0)
                                    d022(i1,i2,i3,0) = d002(i1,i2+1,i3,0) - 2*d002(i1,i2,i3,0) + d002(i1,i2-1,i3,0)
                                    d300i = d200(i1+1,i2,i3,0) - d200(i1-1,i2,i3,0)
                                    d030i = d020(i1,i2+1,i3,0) - d020(i1,i2-1,i3,0)
                                    d003i = d002(i1,i2,i3+1,0) - d002(i1,i2,i3-1,0)
                                    d310i = d200(i1+1,i2+1,i3,0) - d200(i1-1,i2+1,i3,0) - d200(i1+1,i2-1,i3,0) + d200(i1-1,i2-1,i3,0)
                                    d130i = d020(i1+1,i2+1,i3,0) - d020(i1-1,i2+1,i3,0) - d020(i1+1,i2-1,i3,0) + d020(i1-1,i2-1,i3,0)
                                    d301i = d200(i1+1,i2,i3+1,0) - d200(i1-1,i2,i3+1,0) - d200(i1+1,i2,i3-1,0) + d200(i1-1,i2,i3-1,0)
                                    d103i = d002(i1+1,i2,i3+1,0) - d002(i1-1,i2,i3+1,0) - d002(i1+1,i2,i3-1,0) + d002(i1-1,i2,i3-1,0)
                                    d031i = d020(i1,i2+1,i3+1,0) - d020(i1,i2-1,i3+1,0) - d020(i1,i2+1,i3-1,0) + d020(i1,i2-1,i3-1,0)
                                    d013i = d002(i1,i2+1,i3+1,0) - d002(i1,i2-1,i3+1,0) - d002(i1,i2+1,i3-1,0) + d002(i1,i2-1,i3-1,0)
                  ! --- Laplacian to order 4 = lap2h + corrections 
                                    lap4h(i1,i2,i3,0) = lap2h(i1,i2,i3,0) + lapCoeff(i1,i2,i3,0)*crr1*d400(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*css1*d040(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*ctt1*d004(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*(cr1*d310i + cs1*d130i) + lapCoeff(i1,i2,i3,4)*(cr1*d301i + ct1*d103i) + lapCoeff(i1,i2,i3,5)*(cs1*d031i + ct1*d013i) + lapCoeff(i1,i2,i3,6)*cr1 *d300i + lapCoeff(i1,i2,i3,7)*cs1 *d030i + lapCoeff(i1,i2,i3,8)*ct1 *d003i 
                  ! --- Laplacian squared to order 2:
                                    lap2h200(i1,i2,i3,0) = lap2h(i1+1,i2,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1-1,i2,i3,0)
                                    lap2h020(i1,i2,i3,0) = lap2h(i1,i2+1,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2-1,i3,0)
                                    lap2h002(i1,i2,i3,0) = lap2h(i1,i2,i3+1,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2,i3-1,0)
                                    lap2h100i = lap2h(i1+1,i2,i3,0) - lap2h(i1-1,i2,i3,0)
                                    lap2h010i = lap2h(i1,i2+1,i3,0) - lap2h(i1,i2-1,i3,0)
                                    lap2h110i = lap2h(i1+1,i2+1,i3,0) - lap2h(i1-1,i2+1,i3,0) - lap2h(i1+1,i2-1,i3,0) + lap2h(i1-1,i2-1,i3,0)
                                    lap2h001i = lap2h(i1,i2,i3+1,0) - lap2h(i1,i2,i3-1,0)
                                    lap2h101i = lap2h(i1+1,i2,i3+1,0) - lap2h(i1-1,i2,i3+1,0) - lap2h(i1+1,i2,i3-1,0) + lap2h(i1-1,i2,i3-1,0)
                                    lap2h011i = lap2h(i1,i2+1,i3+1,0) - lap2h(i1,i2-1,i3+1,0) - lap2h(i1,i2+1,i3-1,0) + lap2h(i1,i2-1,i3-1,0)
                                    lap2hSq(i1,i2,i3,0) =  lapCoeff(i1,i2,i3,0)*lap2h200(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*lap2h020(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*lap2h002(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*lap2h110i  + lapCoeff(i1,i2,i3,4)*lap2h101i  + lapCoeff(i1,i2,i3,5)*lap2h011i  + lapCoeff(i1,i2,i3,6)*lap2h100i  + lapCoeff(i1,i2,i3,7)*lap2h010i  + lapCoeff(i1,i2,i3,8)*lap2h001i    
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
                                    d600(i1,i2,i3,0) = d400(i1+1,i2,i3,0) - 2*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0)
                                    d060(i1,i2,i3,0) = d040(i1,i2+1,i3,0) - 2*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0)
                                    d006(i1,i2,i3,0) = d004(i1,i2,i3+1,0) - 2*d004(i1,i2,i3,0) + d004(i1,i2,i3-1,0)
                                    d420(i1,i2,i3,0) = d220(i1+1,i2,i3,0) - 2*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0)
                                    d240(i1,i2,i3,0) = d040(i1+1,i2,i3,0) - 2*d040(i1,i2,i3,0) + d040(i1-1,i2,i3,0)
                                    d402(i1,i2,i3,0) = d202(i1+1,i2,i3,0) - 2*d202(i1,i2,i3,0) + d202(i1-1,i2,i3,0)
                                    d204(i1,i2,i3,0) = d004(i1+1,i2,i3,0) - 2*d004(i1,i2,i3,0) + d004(i1-1,i2,i3,0)
                                    d042(i1,i2,i3,0) = d022(i1,i2+1,i3,0) - 2*d022(i1,i2,i3,0) + d022(i1,i2-1,i3,0)
                                    d024(i1,i2,i3,0) = d004(i1,i2+1,i3,0) - 2*d004(i1,i2,i3,0) + d004(i1,i2-1,i3,0)
                                    d500i = d400(i1+1,i2,i3,0) - d400(i1-1,i2,i3,0)
                                    d050i = d040(i1,i2+1,i3,0) - d040(i1,i2-1,i3,0)
                                    d005i = d004(i1,i2,i3+1,0) - d004(i1,i2,i3-1,0)
                                    d510i = d400(i1+1,i2+1,i3,0) - d400(i1-1,i2+1,i3,0) - d400(i1+1,i2-1,i3,0) + d400(i1-1,i2-1,i3,0)
                                    d150i = d040(i1+1,i2+1,i3,0) - d040(i1-1,i2+1,i3,0) - d040(i1+1,i2-1,i3,0) + d040(i1-1,i2-1,i3,0)
                                    d330i = d220(i1+1,i2+1,i3,0) - d220(i1-1,i2+1,i3,0) - d220(i1+1,i2-1,i3,0) + d220(i1-1,i2-1,i3,0)
                                    d501i = d400(i1+1,i2,i3+1,0) - d400(i1-1,i2,i3+1,0) - d400(i1+1,i2,i3-1,0) + d400(i1-1,i2,i3-1,0)
                                    d105i = d004(i1+1,i2,i3+1,0) - d004(i1-1,i2,i3+1,0) - d004(i1+1,i2,i3-1,0) + d004(i1-1,i2,i3-1,0)
                                    d051i = d040(i1,i2+1,i3+1,0) - d040(i1,i2-1,i3+1,0) - d040(i1,i2+1,i3-1,0) + d040(i1,i2-1,i3-1,0)
                                    d015i = d004(i1,i2+1,i3+1,0) - d004(i1,i2-1,i3+1,0) - d004(i1,i2+1,i3-1,0) + d004(i1,i2-1,i3-1,0)
                                    d303i = d202(i1+1,i2,i3+1,0) - d202(i1-1,i2,i3+1,0) - d202(i1+1,i2,i3-1,0) + d202(i1-1,i2,i3-1,0)
                                    d033i = d022(i1,i2+1,i3+1,0) - d022(i1,i2-1,i3+1,0) - d022(i1,i2+1,i3-1,0) + d022(i1,i2-1,i3-1,0)
                  ! --- Laplacian to order 6 = lap4h + corrections 
                                    lap6h(i1,i2,i3,0) = lap4h(i1,i2,i3,0) + lapCoeff(i1,i2,i3,0)*crr2*d600(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*css2*d060(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*ctt2*d006(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*(cr2*d510i + cs2*d150i + cr1*cs1*d330i ) + lapCoeff(i1,i2,i3,4)*(cr2*d501i + ct2*d105i + cr1*ct1*d303i ) + lapCoeff(i1,i2,i3,5)*(cs2*d051i + ct2*d015i + cs1*ct1*d033i ) + lapCoeff(i1,i2,i3,6)*cr2 *d500i + lapCoeff(i1,i2,i3,7)*cs2 *d050i + lapCoeff(i1,i2,i3,8)*ct2 *d005i 
                                    lap2hSq200(i1,i2,i3,0) = lap2hSq(i1+1,i2,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1-1,i2,i3,0)
                                    lap2hSq020(i1,i2,i3,0) = lap2hSq(i1,i2+1,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2-1,i3,0)
                                    lap2hSq002(i1,i2,i3,0) = lap2hSq(i1,i2,i3+1,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2,i3-1,0)
                                    lap2hSq100i = lap2hSq(i1+1,i2,i3,0) - lap2hSq(i1-1,i2,i3,0)
                                    lap2hSq010i = lap2hSq(i1,i2+1,i3,0) - lap2hSq(i1,i2-1,i3,0)
                                    lap2hSq001i = lap2hSq(i1,i2,i3+1,0) - lap2hSq(i1,i2,i3-1,0)
                                    lap2hSq110i = lap2hSq(i1+1,i2+1,i3,0) - lap2hSq(i1-1,i2+1,i3,0) - lap2hSq(i1+1,i2-1,i3,0) + lap2hSq(i1-1,i2-1,i3,0)
                                    lap2hSq101i = lap2hSq(i1+1,i2,i3+1,0) - lap2hSq(i1-1,i2,i3+1,0) - lap2hSq(i1+1,i2,i3-1,0) + lap2hSq(i1-1,i2,i3-1,0)
                                    lap2hSq011i = lap2hSq(i1,i2+1,i3+1,0) - lap2hSq(i1,i2-1,i3+1,0) - lap2hSq(i1,i2+1,i3-1,0) + lap2hSq(i1,i2-1,i3-1,0)
                                    lap2hCubed(i1,i2,i3,0) =  + lapCoeff(i1,i2,i3,0)*lap2hSq200(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*lap2hSq020(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*lap2hSq002(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*lap2hSq110i  + lapCoeff(i1,i2,i3,4)*lap2hSq101i  + lapCoeff(i1,i2,i3,5)*lap2hSq011i  + lapCoeff(i1,i2,i3,6)*lap2hSq100i  + lapCoeff(i1,i2,i3,7)*lap2hSq010i  + lapCoeff(i1,i2,i3,8)*lap2hSq001i   
                                    end if ! mask .ne. 0
                                  end do
                                  end do
                                  end do
               ! --- SPLIT LOOPS
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                if( mask(i1,i2,i3).ne.0 )then
                  ! --- Laplacian squared to order 4 = 
                  !  lap2h*( lap4h ) + corrections*( Lap2h )
                                    lap4h200(i1,i2,i3,0) = lap4h(i1+1,i2,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1-1,i2,i3,0)
                                    lap4h020(i1,i2,i3,0) = lap4h(i1,i2+1,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2-1,i3,0)
                                    lap4h002(i1,i2,i3,0) = lap4h(i1,i2,i3+1,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2,i3-1,0)
                                    lap2h400(i1,i2,i3,0) = lap2h200(i1+1,i2,i3,0) - 2*lap2h200(i1,i2,i3,0) + lap2h200(i1-1,i2,i3,0)
                                    lap2h040(i1,i2,i3,0) = lap2h020(i1,i2+1,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1,i2-1,i3,0)
                                    lap2h004(i1,i2,i3,0) = lap2h002(i1,i2,i3+1,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1,i2,i3-1,0)
                                    lap2h220(i1,i2,i3,0) = lap2h020(i1+1,i2,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1-1,i2,i3,0)
                                    lap2h202(i1,i2,i3,0) = lap2h002(i1+1,i2,i3,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1-1,i2,i3,0)
                                    lap2h022(i1,i2,i3,0) = lap2h002(i1,i2+1,i3,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1,i2-1,i3,0)
                                    lap4h100i = lap4h(i1+1,i2,i3,0) - lap4h(i1-1,i2,i3,0)
                                    lap4h010i = lap4h(i1,i2+1,i3,0) - lap4h(i1,i2-1,i3,0)
                                    lap4h001i = lap4h(i1,i2,i3+1,0) - lap4h(i1,i2,i3-1,0)
                                    lap4h110i = lap4h(i1+1,i2+1,i3,0) - lap4h(i1-1,i2+1,i3,0) - lap4h(i1+1,i2-1,i3,0) + lap4h(i1-1,i2-1,i3,0)
                                    lap4h101i = lap4h(i1+1,i2,i3+1,0) - lap4h(i1-1,i2,i3+1,0) - lap4h(i1+1,i2,i3-1,0) + lap4h(i1-1,i2,i3-1,0)
                                    lap4h011i = lap4h(i1,i2+1,i3+1,0) - lap4h(i1,i2-1,i3+1,0) - lap4h(i1,i2+1,i3-1,0) + lap4h(i1,i2-1,i3-1,0)
                                    lap2h300i = lap2h200(i1+1,i2,i3,0) - lap2h200(i1-1,i2,i3,0)
                                    lap2h030i = lap2h020(i1,i2+1,i3,0) - lap2h020(i1,i2-1,i3,0)
                                    lap2h003i = lap2h002(i1,i2,i3+1,0) - lap2h002(i1,i2,i3-1,0)
                                    lap2h310i = lap2h200(i1+1,i2+1,i3,0) - lap2h200(i1-1,i2+1,i3,0) - lap2h200(i1+1,i2-1,i3,0) + lap2h200(i1-1,i2-1,i3,0)
                                    lap2h130i = lap2h020(i1+1,i2+1,i3,0) - lap2h020(i1-1,i2+1,i3,0) - lap2h020(i1+1,i2-1,i3,0) + lap2h020(i1-1,i2-1,i3,0)
                                    lap2h301i = lap2h200(i1+1,i2,i3+1,0) - lap2h200(i1-1,i2,i3+1,0) - lap2h200(i1+1,i2,i3-1,0) + lap2h200(i1-1,i2,i3-1,0)
                                    lap2h103i = lap2h002(i1+1,i2,i3+1,0) - lap2h002(i1-1,i2,i3+1,0) - lap2h002(i1+1,i2,i3-1,0) + lap2h002(i1-1,i2,i3-1,0)
                                    lap2h031i = lap2h020(i1,i2+1,i3+1,0) - lap2h020(i1,i2-1,i3+1,0) - lap2h020(i1,i2+1,i3-1,0) + lap2h020(i1,i2-1,i3-1,0)
                                    lap2h013i = lap2h002(i1,i2+1,i3+1,0) - lap2h002(i1,i2-1,i3+1,0) - lap2h002(i1,i2+1,i3-1,0) + lap2h002(i1,i2-1,i3-1,0)
                                    lap4hSq(i1,i2,i3,0) =     lapCoeff(i1,i2,i3,0)*( lap4h200(i1,i2,i3,0) + crr1*lap2h400(i1,i2,i3,0) )    + lapCoeff(i1,i2,i3,1)*( lap4h020(i1,i2,i3,0) + css1*lap2h040(i1,i2,i3,0) )     + lapCoeff(i1,i2,i3,2)*( lap4h002(i1,i2,i3,0) + ctt1*lap2h004(i1,i2,i3,0) )     + lapCoeff(i1,i2,i3,3)*( lap4h110i + cr1*lap2h310i + cs1*lap2h130i ) + lapCoeff(i1,i2,i3,4)*( lap4h101i + cr1*lap2h301i + ct1*lap2h103i ) + lapCoeff(i1,i2,i3,5)*( lap4h011i + cs1*lap2h031i + ct1*lap2h013i ) + lapCoeff(i1,i2,i3,6)*( lap4h100i + cr1 *lap2h300i )    + lapCoeff(i1,i2,i3,7)*( lap4h010i + cs1 *lap2h030i )    + lapCoeff(i1,i2,i3,8)*( lap4h001i + ct1 *lap2h003i )      
                                    end if ! mask .ne. 0
                                  end do
                                  end do
                                  end do
              ! ===========  FINAL LOOP TO FILL IN THE SOLUTION ============
                            numGhost1=0;
                            n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                            n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                            n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                if( mask(i1,i2,i3).ne.0 )then
                                    d800i = d600(i1+1,i2,i3,0) - 2*d600(i1,i2,i3,0) + d600(i1-1,i2,i3,0)
                                    d080i = d060(i1,i2+1,i3,0) - 2*d060(i1,i2,i3,0) + d060(i1,i2-1,i3,0)
                                    d008i = d006(i1,i2,i3+1,0) - 2*d006(i1,i2,i3,0) + d006(i1,i2,i3-1,0)
                                    d700i = d600(i1+1,i2,i3,0) - d600(i1-1,i2,i3,0)
                                    d070i = d060(i1,i2+1,i3,0) - d060(i1,i2-1,i3,0)
                                    d007i = d006(i1,i2,i3+1,0) - d006(i1,i2,i3-1,0)
                                    d710i = d600(i1+1,i2+1,i3,0) - d600(i1-1,i2+1,i3,0) - d600(i1+1,i2-1,i3,0) + d600(i1-1,i2-1,i3,0)
                                    d170i = d060(i1+1,i2+1,i3,0) - d060(i1-1,i2+1,i3,0) - d060(i1+1,i2-1,i3,0) + d060(i1-1,i2-1,i3,0)
                                    d701i = d600(i1+1,i2,i3+1,0) - d600(i1-1,i2,i3+1,0) - d600(i1+1,i2,i3-1,0) + d600(i1-1,i2,i3-1,0)
                                    d107i = d006(i1+1,i2,i3+1,0) - d006(i1-1,i2,i3+1,0) - d006(i1+1,i2,i3-1,0) + d006(i1-1,i2,i3-1,0)
                                    d071i = d060(i1,i2+1,i3+1,0) - d060(i1,i2-1,i3+1,0) - d060(i1,i2+1,i3-1,0) + d060(i1,i2-1,i3-1,0)
                                    d017i = d006(i1,i2+1,i3+1,0) - d006(i1,i2-1,i3+1,0) - d006(i1,i2+1,i3-1,0) + d006(i1,i2-1,i3-1,0)
                                    d530i = d420(i1+1,i2+1,i3,0) - d420(i1-1,i2+1,i3,0) - d420(i1+1,i2-1,i3,0) + d420(i1-1,i2-1,i3,0)
                                    d350i = d240(i1+1,i2+1,i3,0) - d240(i1-1,i2+1,i3,0) - d240(i1+1,i2-1,i3,0) + d240(i1-1,i2-1,i3,0)
                                    d503i = d402(i1+1,i2,i3+1,0) - d402(i1-1,i2,i3+1,0) - d402(i1+1,i2,i3-1,0) + d402(i1-1,i2,i3-1,0)
                                    d305i = d204(i1+1,i2,i3+1,0) - d204(i1-1,i2,i3+1,0) - d204(i1+1,i2,i3-1,0) + d204(i1-1,i2,i3-1,0)
                                    d053i = d042(i1,i2+1,i3+1,0) - d042(i1,i2-1,i3+1,0) - d042(i1,i2+1,i3-1,0) + d042(i1,i2-1,i3-1,0)
                                    d035i = d024(i1,i2+1,i3+1,0) - d024(i1,i2-1,i3+1,0) - d024(i1,i2+1,i3-1,0) + d024(i1,i2-1,i3-1,0)
                  ! --- Laplacian to order 8 = lap6h + corrections 
                                    lap8h = lap6h(i1,i2,i3,0)                                                         + lapCoeff(i1,i2,i3,0)*crr3*d800i                                               + lapCoeff(i1,i2,i3,1)*css3*d080i                                               + lapCoeff(i1,i2,i3,2)*ctt3*d008i                                               + lapCoeff(i1,i2,i3,3)*(cr3*d710i + cs3*d170i + cr2*cs1*d530i + cr1*cs2*d350i ) + lapCoeff(i1,i2,i3,4)*(cr3*d701i + ct3*d107i + cr2*ct1*d503i + cr1*ct2*d305i ) + lapCoeff(i1,i2,i3,5)*(cs3*d071i + ct3*d017i + cs2*ct1*d053i + cs1*ct2*d035i ) + lapCoeff(i1,i2,i3,6)* cr3*d700i                                               + lapCoeff(i1,i2,i3,7)* cs3*d070i                                               + lapCoeff(i1,i2,i3,8)* ct3*d007i 
                  ! --- Laplacian^4 4p (4th power) order 2: 
                                    lap2hCubed200i = lap2hCubed(i1+1,i2,i3,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1-1,i2,i3,0)
                                    lap2hCubed020i = lap2hCubed(i1,i2+1,i3,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1,i2-1,i3,0)
                                    lap2hCubed002i = lap2hCubed(i1,i2,i3+1,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1,i2,i3-1,0)
                                    lap2hCubed100i = lap2hCubed(i1+1,i2,i3,0) - lap2hCubed(i1-1,i2,i3,0)
                                    lap2hCubed010i = lap2hCubed(i1,i2+1,i3,0) - lap2hCubed(i1,i2-1,i3,0)
                                    lap2hCubed001i = lap2hCubed(i1,i2,i3+1,0) - lap2hCubed(i1,i2,i3-1,0)
                                    lap2hCubed110i = lap2hCubed(i1+1,i2+1,i3,0) - lap2hCubed(i1-1,i2+1,i3,0) - lap2hCubed(i1+1,i2-1,i3,0) + lap2hCubed(i1-1,i2-1,i3,0)
                                    lap2hCubed101i = lap2hCubed(i1+1,i2,i3+1,0) - lap2hCubed(i1-1,i2,i3+1,0) - lap2hCubed(i1+1,i2,i3-1,0) + lap2hCubed(i1-1,i2,i3-1,0)
                                    lap2hCubed011i = lap2hCubed(i1,i2+1,i3+1,0) - lap2hCubed(i1,i2-1,i3+1,0) - lap2hCubed(i1,i2+1,i3-1,0) + lap2hCubed(i1,i2-1,i3-1,0)
                                    lap2h4p  =                             + lapCoeff(i1,i2,i3,0)*lap2hCubed200i  + lapCoeff(i1,i2,i3,1)*lap2hCubed020i  + lapCoeff(i1,i2,i3,2)*lap2hCubed002i  + lapCoeff(i1,i2,i3,3)*lap2hCubed110i  + lapCoeff(i1,i2,i3,4)*lap2hCubed101i  + lapCoeff(i1,i2,i3,5)*lap2hCubed011i  + lapCoeff(i1,i2,i3,6)*lap2hCubed100i  + lapCoeff(i1,i2,i3,7)*lap2hCubed010i  + lapCoeff(i1,i2,i3,8)*lap2hCubed001i    
                  ! --- Laplacian squared to order 6 :
                  !   Lap6h = Lap4h + M4  = (Lap2h) + M2 + M4 
                  !   Lap6h*Lap6h = [ (Lap2h) + M2 + M4 ] [ (Lap2h) + M2 + M4 ]
                  !               = Lap2h*Lap6h + M2*Lap4h + M4*Lap2h + O(h^6)
                                    lap6h200i = lap6h(i1+1,i2,i3,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1-1,i2,i3,0)
                                    lap6h020i = lap6h(i1,i2+1,i3,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1,i2-1,i3,0)
                                    lap6h002i = lap6h(i1,i2,i3+1,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1,i2,i3-1,0)
                                    lap6h100i = lap6h(i1+1,i2,i3,0) - lap6h(i1-1,i2,i3,0)
                                    lap6h010i = lap6h(i1,i2+1,i3,0) - lap6h(i1,i2-1,i3,0)
                                    lap6h001i = lap6h(i1,i2,i3+1,0) - lap6h(i1,i2,i3-1,0)
                                    lap6h110i = lap6h(i1+1,i2+1,i3,0) - lap6h(i1-1,i2+1,i3,0) - lap6h(i1+1,i2-1,i3,0) + lap6h(i1-1,i2-1,i3,0)
                                    lap6h101i = lap6h(i1+1,i2,i3+1,0) - lap6h(i1-1,i2,i3+1,0) - lap6h(i1+1,i2,i3-1,0) + lap6h(i1-1,i2,i3-1,0)
                                    lap6h011i = lap6h(i1,i2+1,i3+1,0) - lap6h(i1,i2-1,i3+1,0) - lap6h(i1,i2+1,i3-1,0) + lap6h(i1,i2-1,i3-1,0)
                                    lap4h400i = lap4h200(i1+1,i2,i3,0) - 2*lap4h200(i1,i2,i3,0) + lap4h200(i1-1,i2,i3,0)
                                    lap4h040i = lap4h020(i1,i2+1,i3,0) - 2*lap4h020(i1,i2,i3,0) + lap4h020(i1,i2-1,i3,0)
                                    lap4h004i = lap4h002(i1,i2,i3+1,0) - 2*lap4h002(i1,i2,i3,0) + lap4h002(i1,i2,i3-1,0)
                                    lap4h300i = lap4h200(i1+1,i2,i3,0) - lap4h200(i1-1,i2,i3,0)
                                    lap4h030i = lap4h020(i1,i2+1,i3,0) - lap4h020(i1,i2-1,i3,0)
                                    lap4h003i = lap4h002(i1,i2,i3+1,0) - lap4h002(i1,i2,i3-1,0)
                                    lap4h310i = lap4h200(i1+1,i2+1,i3,0) - lap4h200(i1-1,i2+1,i3,0) - lap4h200(i1+1,i2-1,i3,0) + lap4h200(i1-1,i2-1,i3,0)
                                    lap4h130i = lap4h020(i1+1,i2+1,i3,0) - lap4h020(i1-1,i2+1,i3,0) - lap4h020(i1+1,i2-1,i3,0) + lap4h020(i1-1,i2-1,i3,0)
                                    lap4h301i = lap4h200(i1+1,i2,i3+1,0) - lap4h200(i1-1,i2,i3+1,0) - lap4h200(i1+1,i2,i3-1,0) + lap4h200(i1-1,i2,i3-1,0)
                                    lap4h103i = lap4h002(i1+1,i2,i3+1,0) - lap4h002(i1-1,i2,i3+1,0) - lap4h002(i1+1,i2,i3-1,0) + lap4h002(i1-1,i2,i3-1,0)
                                    lap4h031i = lap4h020(i1,i2+1,i3+1,0) - lap4h020(i1,i2-1,i3+1,0) - lap4h020(i1,i2+1,i3-1,0) + lap4h020(i1,i2-1,i3-1,0)
                                    lap4h013i = lap4h002(i1,i2+1,i3+1,0) - lap4h002(i1,i2-1,i3+1,0) - lap4h002(i1,i2+1,i3-1,0) + lap4h002(i1,i2-1,i3-1,0)
                                    lap2h600i = lap2h400(i1+1,i2,i3,0) - 2*lap2h400(i1,i2,i3,0) + lap2h400(i1-1,i2,i3,0)
                                    lap2h060i = lap2h040(i1,i2+1,i3,0) - 2*lap2h040(i1,i2,i3,0) + lap2h040(i1,i2-1,i3,0)
                                    lap2h006i = lap2h004(i1,i2,i3+1,0) - 2*lap2h004(i1,i2,i3,0) + lap2h004(i1,i2,i3-1,0)
                                    lap2h500i = lap2h400(i1+1,i2,i3,0) - lap2h400(i1-1,i2,i3,0)
                                    lap2h050i = lap2h040(i1,i2+1,i3,0) - lap2h040(i1,i2-1,i3,0)
                                    lap2h005i = lap2h004(i1,i2,i3+1,0) - lap2h004(i1,i2,i3-1,0)
                                    lap2h510i = lap2h400(i1+1,i2+1,i3,0) - lap2h400(i1-1,i2+1,i3,0) - lap2h400(i1+1,i2-1,i3,0) + lap2h400(i1-1,i2-1,i3,0)
                                    lap2h150i = lap2h040(i1+1,i2+1,i3,0) - lap2h040(i1-1,i2+1,i3,0) - lap2h040(i1+1,i2-1,i3,0) + lap2h040(i1-1,i2-1,i3,0)
                                    lap2h330i = lap2h220(i1+1,i2+1,i3,0) - lap2h220(i1-1,i2+1,i3,0) - lap2h220(i1+1,i2-1,i3,0) + lap2h220(i1-1,i2-1,i3,0)
                                    lap2h501i = lap2h400(i1+1,i2,i3+1,0) - lap2h400(i1-1,i2,i3+1,0) - lap2h400(i1+1,i2,i3-1,0) + lap2h400(i1-1,i2,i3-1,0)
                                    lap2h105i = lap2h004(i1+1,i2,i3+1,0) - lap2h004(i1-1,i2,i3+1,0) - lap2h004(i1+1,i2,i3-1,0) + lap2h004(i1-1,i2,i3-1,0)
                                    lap2h051i = lap2h040(i1,i2+1,i3+1,0) - lap2h040(i1,i2-1,i3+1,0) - lap2h040(i1,i2+1,i3-1,0) + lap2h040(i1,i2-1,i3-1,0)
                                    lap2h015i = lap2h004(i1,i2+1,i3+1,0) - lap2h004(i1,i2-1,i3+1,0) - lap2h004(i1,i2+1,i3-1,0) + lap2h004(i1,i2-1,i3-1,0)
                                    lap2h303i = lap2h202(i1+1,i2,i3+1,0) - lap2h202(i1-1,i2,i3+1,0) - lap2h202(i1+1,i2,i3-1,0) + lap2h202(i1-1,i2,i3-1,0)
                                    lap2h033i = lap2h022(i1,i2+1,i3+1,0) - lap2h022(i1,i2-1,i3+1,0) - lap2h022(i1,i2+1,i3-1,0) + lap2h022(i1,i2-1,i3-1,0)
                                    lap6hSq =                                                                                     lapCoeff(i1,i2,i3,0)*(lap6h200i + crr1*lap4h400i + crr2*lap2h600i )                       + lapCoeff(i1,i2,i3,1)*(lap6h020i + css1*lap4h040i + css2*lap2h060i )                       + lapCoeff(i1,i2,i3,2)*(lap6h002i + ctt1*lap4h004i + ctt2*lap2h006i )                       + lapCoeff(i1,i2,i3,3)*(lap6h110i +  cr1*lap4h310i +  cr2*lap2h510i                         +  cs1*lap4h130i +  cs2*lap2h150i + cr1*cs1*lap2h330i )   + lapCoeff(i1,i2,i3,4)*(lap6h101i +  cr1*lap4h301i +  cr2*lap2h501i                         +  ct1*lap4h103i +  ct2*lap2h105i + cr1*ct1*lap2h303i )   + lapCoeff(i1,i2,i3,5)*(lap6h011i +  cs1*lap4h031i +  cs2*lap2h051i                         +  ct1*lap4h013i +  ct2*lap2h015i + cs1*ct1*lap2h033i )   + lapCoeff(i1,i2,i3,6)*(lap6h100i +  cr1*lap4h300i +  cr2*lap2h500i )                       + lapCoeff(i1,i2,i3,7)*(lap6h010i +  cs1*lap4h030i +  cs2*lap2h050i )                       + lapCoeff(i1,i2,i3,8)*(lap6h010i +  ct1*lap4h003i +  ct2*lap2h005i )                         
                  ! --- Laplacian CUBED to order 4 
                                    lap4hSq200i = lap4hSq(i1+1,i2,i3,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1-1,i2,i3,0)
                                    lap4hSq020i = lap4hSq(i1,i2+1,i3,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1,i2-1,i3,0)
                                    lap4hSq002i = lap4hSq(i1,i2,i3+1,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1,i2,i3-1,0)
                                    lap4hSq100i = lap4hSq(i1+1,i2,i3,0) - lap4hSq(i1-1,i2,i3,0)
                                    lap4hSq010i = lap4hSq(i1,i2+1,i3,0) - lap4hSq(i1,i2-1,i3,0)
                                    lap4hSq001i = lap4hSq(i1,i2,i3+1,0) - lap4hSq(i1,i2,i3-1,0)
                                    lap4hSq110i = lap4hSq(i1+1,i2+1,i3,0) - lap4hSq(i1-1,i2+1,i3,0) - lap4hSq(i1+1,i2-1,i3,0) + lap4hSq(i1-1,i2-1,i3,0)
                                    lap4hSq101i = lap4hSq(i1+1,i2,i3+1,0) - lap4hSq(i1-1,i2,i3+1,0) - lap4hSq(i1+1,i2,i3-1,0) + lap4hSq(i1-1,i2,i3-1,0)
                                    lap4hSq011i = lap4hSq(i1,i2+1,i3+1,0) - lap4hSq(i1,i2-1,i3+1,0) - lap4hSq(i1,i2+1,i3-1,0) + lap4hSq(i1,i2-1,i3-1,0)
                                    lap2hSq400i = lap2hSq200(i1+1,i2,i3,0) - 2*lap2hSq200(i1,i2,i3,0) + lap2hSq200(i1-1,i2,i3,0)
                                    lap2hSq040i = lap2hSq020(i1,i2+1,i3,0) - 2*lap2hSq020(i1,i2,i3,0) + lap2hSq020(i1,i2-1,i3,0)
                                    lap2hSq004i = lap2hSq002(i1,i2,i3+1,0) - 2*lap2hSq002(i1,i2,i3,0) + lap2hSq002(i1,i2,i3-1,0)
                                    lap2hSq300i = lap2hSq200(i1+1,i2,i3,0) - lap2hSq200(i1-1,i2,i3,0)
                                    lap2hSq030i = lap2hSq020(i1,i2+1,i3,0) - lap2hSq020(i1,i2-1,i3,0)
                                    lap2hSq003i = lap2hSq002(i1,i2,i3+1,0) - lap2hSq002(i1,i2,i3-1,0)
                                    lap2hSq310i = lap2hSq200(i1+1,i2+1,i3,0) - lap2hSq200(i1-1,i2+1,i3,0) - lap2hSq200(i1+1,i2-1,i3,0) + lap2hSq200(i1-1,i2-1,i3,0)
                                    lap2hSq130i = lap2hSq020(i1+1,i2+1,i3,0) - lap2hSq020(i1-1,i2+1,i3,0) - lap2hSq020(i1+1,i2-1,i3,0) + lap2hSq020(i1-1,i2-1,i3,0)
                                    lap2hSq301i = lap2hSq200(i1+1,i2,i3+1,0) - lap2hSq200(i1-1,i2,i3+1,0) - lap2hSq200(i1+1,i2,i3-1,0) + lap2hSq200(i1-1,i2,i3-1,0)
                                    lap2hSq103i = lap2hSq002(i1+1,i2,i3+1,0) - lap2hSq002(i1-1,i2,i3+1,0) - lap2hSq002(i1+1,i2,i3-1,0) + lap2hSq002(i1-1,i2,i3-1,0)
                                    lap2hSq031i = lap2hSq020(i1,i2+1,i3+1,0) - lap2hSq020(i1,i2-1,i3+1,0) - lap2hSq020(i1,i2+1,i3-1,0) + lap2hSq020(i1,i2-1,i3-1,0)
                                    lap2hSq013i = lap2hSq002(i1,i2+1,i3+1,0) - lap2hSq002(i1,i2-1,i3+1,0) - lap2hSq002(i1,i2+1,i3-1,0) + lap2hSq002(i1,i2-1,i3-1,0)
                                    lap4hCubed =                                                    lapCoeff(i1,i2,i3,0)*(lap4hSq200i + crr1*lap2hSq400i )   + lapCoeff(i1,i2,i3,1)*(lap4hSq020i + css1*lap2hSq040i )   + lapCoeff(i1,i2,i3,2)*(lap4hSq002i + ctt1*lap2hSq004i )   + lapCoeff(i1,i2,i3,3)*(lap4hSq110i +  cr1*lap2hSq310i     +  cs1*lap2hSq130i )   + lapCoeff(i1,i2,i3,4)*(lap4hSq101i +  cr1*lap2hSq301i     +  ct1*lap2hSq103i )   + lapCoeff(i1,i2,i3,5)*(lap4hSq011i +  cs1*lap2hSq031i     +  ct1*lap2hSq013i )   + lapCoeff(i1,i2,i3,6)*(lap4hSq100i + cr1 *lap2hSq300i )   + lapCoeff(i1,i2,i3,7)*(lap4hSq010i + cs1 *lap2hSq030i )   + lapCoeff(i1,i2,i3,8)*(lap4hSq001i + ct1 *lap2hSq003i )     
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
                                    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m)  + cdtsq*( lap8h )               + cdtPow4By12*( lap6hSq )       + cdtPow6By360*( lap4hCubed )   + cdtPow8By20160*( lap2h4p )    +dtSq*fv(m)                    
                                    end if ! mask .ne. 0
                                  end do
                                  end do
                                  end do
                    end if
              else
                      if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
                          write(*,'("advWaveME: ADVANCE dim=3 order=8 orderInTime=8, grid=curvilinear... t=",e10.2)') t
                      end if
                      m=0 ! component number 
                      ec = 0 ! component number 
           ! -- call the appropriate macro:
           !  update2dOrder2Rectangular(3,8,8,curvilinear)
           !  update3dOrder6Curvilinear(3,8,8,curvilinear)
                      if( useMask.eq.0 .and. addForcing.eq.0 )then
             ! No-mask, no-forcing
               ! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
               ! Example: 
               ! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...
cr0 = 1.; cr1 = -1/6.; cr2 = 1/30.; cr3 = -1/140.; 
cs0 = 1.; cs1 = -1/6.; cs2 = 1/30.; cs3 = -1/140.; 
ct0 = 1.; ct1 = -1/6.; ct2 = 1/30.; ct3 = -1/140.; 
crr0 = 1.; crr1 = -1/12.; crr2 = 1/90.; crr3 = -1/560.; 
css0 = 1.; css1 = -1/12.; css2 = 1/90.; css3 = -1/560.; 
ctt0 = 1.; ctt1 = -1/12.; ctt2 = 1/90.; ctt3 = -1/560.; 
crrr0 = 1.; crrr1 = -1/4.; crrr2 = 7/120.; crrr3 = -41/3024.; 
csss0 = 1.; csss1 = -1/4.; csss2 = 7/120.; csss3 = -41/3024.; 
cttt0 = 1.; cttt1 = -1/4.; cttt2 = 7/120.; cttt3 = -41/3024.; 
crrrr0 = 1.; crrrr1 = -1/6.; crrrr2 = 7/240.; crrrr3 = -41/7560.; 
cssss0 = 1.; cssss1 = -1/6.; cssss2 = 7/240.; cssss3 = -41/7560.; 
ctttt0 = 1.; ctttt1 = -1/6.; ctttt2 = 7/240.; ctttt3 = -41/7560.; 
crrrrr0 = 1.; crrrrr1 = -1/3.; crrrrr2 = 13/144.; crrrrr3 = -139/6048.; 
csssss0 = 1.; csssss1 = -1/3.; csssss2 = 13/144.; csssss3 = -139/6048.; 
cttttt0 = 1.; cttttt1 = -1/3.; cttttt2 = 13/144.; cttttt3 = -139/6048.; 
crrrrrr0 = 1.; crrrrrr1 = -1/4.; crrrrrr2 = 13/240.; crrrrrr3 = -139/12096.; 
cssssss0 = 1.; cssssss1 = -1/4.; cssssss2 = 13/240.; cssssss3 = -139/12096.; 
ctttttt0 = 1.; ctttttt1 = -1/4.; ctttttt2 = 13/240.; ctttttt3 = -139/12096.; 
                              dr1=dr(0); dr1i=1./dr1;
                              dr2=dr(1); dr2i=1./dr2;
                              dr3=dr(2); dr3i=1./dr3;
                              fv(m)=0.
                              if( lapCoeff(0,0,0,0).le.0. )then
                 ! --- Evaluate and store coefficients in Laplacian ---
                                  write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'
                                  numGhost1=3;
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
                                          diffOrder1=8
                                      elseif( (i1-3).ge.nd1a .and. (i1+3).le.nd1b )then
                                          diffOrder1=6
                                      elseif( (i1-2).ge.nd1a .and. (i1+2).le.nd1b )then
                                          diffOrder1=4
                                      elseif( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
                                          diffOrder1=2
                                      else
                                          stop 999
                                      end if
                                      if( (i2-4).ge.nd2a .and. (i2+4).le.nd2b )then
                                          diffOrder2=8
                                      elseif( (i2-3).ge.nd2a .and. (i2+3).le.nd2b )then
                                          diffOrder2=6
                                      elseif( (i2-2).ge.nd2a .and. (i2+2).le.nd2b )then
                                          diffOrder2=4
                                      elseif( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
                                          diffOrder2=2
                                      else
                                          stop 999
                                      end if
                                      if( (i3-4).ge.nd3a .and. (i3+4).le.nd3b )then
                                          diffOrder3=8
                                      elseif( (i3-3).ge.nd3a .and. (i3+3).le.nd3b )then
                                          diffOrder3=6
                                      elseif( (i3-2).ge.nd3a .and. (i3+2).le.nd3b )then
                                          diffOrder3=4
                                      elseif( (i3-1).ge.nd3a .and. (i3+1).le.nd3b )then
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
                                          rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
                                          rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
                                          rzs = ( 8*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) )*(dr2i/12.) 
                                          sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
                                          sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
                                          szs = ( 8*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) )*(dr2i/12.) 
                                          txs = ( 8*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) )*(dr2i/12.) 
                                          tys = ( 8*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) )*(dr2i/12.) 
                                          tzs = ( 8*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) )*(dr2i/12.) 
                                      elseif( diffOrder2.eq.8 )then
                                          rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
                                          rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
                                          rzs = ( 8*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) )*(dr2i/12.) 
                                          sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
                                          sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
                                          szs = ( 8*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) )*(dr2i/12.) 
                                          txs = ( 8*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) )*(dr2i/12.) 
                                          tys = ( 8*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) )*(dr2i/12.) 
                                          tzs = ( 8*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) )*(dr2i/12.) 
                                      end if
                                      if( diffOrder3.eq.2 )then
                                          rxt = (rsxy(i1,i2,i3+1,0,0)-rsxy(i1,i2,i3-1,0,0))*(.5*dr2i) 
                                          ryt = (rsxy(i1,i2,i3+1,0,1)-rsxy(i1,i2,i3-1,0,1))*(.5*dr2i) 
                                          rzt = (rsxy(i1,i2,i3+1,0,2)-rsxy(i1,i2,i3-1,0,2))*(.5*dr2i) 
                                          sxt = (rsxy(i1,i2,i3+1,1,0)-rsxy(i1,i2,i3-1,1,0))*(.5*dr2i) 
                                          syt = (rsxy(i1,i2,i3+1,1,1)-rsxy(i1,i2,i3-1,1,1))*(.5*dr2i) 
                                          szt = (rsxy(i1,i2,i3+1,1,2)-rsxy(i1,i2,i3-1,1,2))*(.5*dr2i) 
                                          txt = (rsxy(i1,i2,i3+1,2,0)-rsxy(i1,i2,i3-1,2,0))*(.5*dr2i) 
                                          tyt = (rsxy(i1,i2,i3+1,2,1)-rsxy(i1,i2,i3-1,2,1))*(.5*dr2i) 
                                          tzt = (rsxy(i1,i2,i3+1,2,2)-rsxy(i1,i2,i3-1,2,2))*(.5*dr2i) 
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
                                      lapCoeff(i1,i2,i3,3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i*.25
                                      lapCoeff(i1,i2,i3,4) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr2i*.25
                                      lapCoeff(i1,i2,i3,5) = 2.*(sx*tx + sy*ty + sz*tz )*dr1i*dr2i*.25
                                      lapCoeff(i1,i2,i3,6) = (rxx + ryy + rzz)*dr1i*.5
                                      lapCoeff(i1,i2,i3,7) = (sxx + syy + tyy)*dr2i*.5 
                                      lapCoeff(i1,i2,i3,8) = (txx + tyy + tzz)*dr3i*.5 
                                    end do
                                    end do
                                    end do
                              end if ! end assignLapCoeff
                              numGhost1=3;
                              n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                              n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                              n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                do i3=n3a,n3b
                                do i2=n2a,n2b
                                do i1=n1a,n1b
                                      d200(i1,i2,i3,0) = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
                                      d020(i1,i2,i3,0) = u(i1,i2+1,i3,0) - 2*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
                                      d002(i1,i2,i3,0) = u(i1,i2,i3+1,0) - 2*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
                                      d100i = u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0)
                                      d010i = u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0)
                                      d110i = u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0)
                                      d001i = u(i1,i2,i3+1,0) - u(i1,i2,i3-1,0)
                                      d101i = u(i1+1,i2,i3+1,0) - u(i1-1,i2,i3+1,0) - u(i1+1,i2,i3-1,0) + u(i1-1,i2,i3-1,0)
                                      d011i = u(i1,i2+1,i3+1,0) - u(i1,i2-1,i3+1,0) - u(i1,i2+1,i3-1,0) + u(i1,i2-1,i3-1,0)
                                      lap2h(i1,i2,i3,0) = lapCoeff(i1,i2,i3,0)*d200(i1,i2,i3,0) +lapCoeff(i1,i2,i3,1)*d020(i1,i2,i3,0) +lapCoeff(i1,i2,i3,2)*d002(i1,i2,i3,0) +lapCoeff(i1,i2,i3,3)*d110i + lapCoeff(i1,i2,i3,4)*d101i + lapCoeff(i1,i2,i3,5)*d011i + lapCoeff(i1,i2,i3,6)*d100i + lapCoeff(i1,i2,i3,7)*d010i + lapCoeff(i1,i2,i3,8)*d001i
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
                                      d400(i1,i2,i3,0) = d200(i1+1,i2,i3,0) - 2*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0)
                                      d040(i1,i2,i3,0) = d020(i1,i2+1,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1,i2-1,i3,0)
                                      d004(i1,i2,i3,0) = d002(i1,i2,i3+1,0) - 2*d002(i1,i2,i3,0) + d002(i1,i2,i3-1,0)
                                      d220(i1,i2,i3,0) = d020(i1+1,i2,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1-1,i2,i3,0)
                                      d202(i1,i2,i3,0) = d002(i1+1,i2,i3,0) - 2*d002(i1,i2,i3,0) + d002(i1-1,i2,i3,0)
                                      d022(i1,i2,i3,0) = d002(i1,i2+1,i3,0) - 2*d002(i1,i2,i3,0) + d002(i1,i2-1,i3,0)
                                      d300i = d200(i1+1,i2,i3,0) - d200(i1-1,i2,i3,0)
                                      d030i = d020(i1,i2+1,i3,0) - d020(i1,i2-1,i3,0)
                                      d003i = d002(i1,i2,i3+1,0) - d002(i1,i2,i3-1,0)
                                      d310i = d200(i1+1,i2+1,i3,0) - d200(i1-1,i2+1,i3,0) - d200(i1+1,i2-1,i3,0) + d200(i1-1,i2-1,i3,0)
                                      d130i = d020(i1+1,i2+1,i3,0) - d020(i1-1,i2+1,i3,0) - d020(i1+1,i2-1,i3,0) + d020(i1-1,i2-1,i3,0)
                                      d301i = d200(i1+1,i2,i3+1,0) - d200(i1-1,i2,i3+1,0) - d200(i1+1,i2,i3-1,0) + d200(i1-1,i2,i3-1,0)
                                      d103i = d002(i1+1,i2,i3+1,0) - d002(i1-1,i2,i3+1,0) - d002(i1+1,i2,i3-1,0) + d002(i1-1,i2,i3-1,0)
                                      d031i = d020(i1,i2+1,i3+1,0) - d020(i1,i2-1,i3+1,0) - d020(i1,i2+1,i3-1,0) + d020(i1,i2-1,i3-1,0)
                                      d013i = d002(i1,i2+1,i3+1,0) - d002(i1,i2-1,i3+1,0) - d002(i1,i2+1,i3-1,0) + d002(i1,i2-1,i3-1,0)
                   ! --- Laplacian to order 4 = lap2h + corrections 
                                      lap4h(i1,i2,i3,0) = lap2h(i1,i2,i3,0) + lapCoeff(i1,i2,i3,0)*crr1*d400(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*css1*d040(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*ctt1*d004(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*(cr1*d310i + cs1*d130i) + lapCoeff(i1,i2,i3,4)*(cr1*d301i + ct1*d103i) + lapCoeff(i1,i2,i3,5)*(cs1*d031i + ct1*d013i) + lapCoeff(i1,i2,i3,6)*cr1 *d300i + lapCoeff(i1,i2,i3,7)*cs1 *d030i + lapCoeff(i1,i2,i3,8)*ct1 *d003i 
                   ! --- Laplacian squared to order 2:
                                      lap2h200(i1,i2,i3,0) = lap2h(i1+1,i2,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1-1,i2,i3,0)
                                      lap2h020(i1,i2,i3,0) = lap2h(i1,i2+1,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2-1,i3,0)
                                      lap2h002(i1,i2,i3,0) = lap2h(i1,i2,i3+1,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2,i3-1,0)
                                      lap2h100i = lap2h(i1+1,i2,i3,0) - lap2h(i1-1,i2,i3,0)
                                      lap2h010i = lap2h(i1,i2+1,i3,0) - lap2h(i1,i2-1,i3,0)
                                      lap2h110i = lap2h(i1+1,i2+1,i3,0) - lap2h(i1-1,i2+1,i3,0) - lap2h(i1+1,i2-1,i3,0) + lap2h(i1-1,i2-1,i3,0)
                                      lap2h001i = lap2h(i1,i2,i3+1,0) - lap2h(i1,i2,i3-1,0)
                                      lap2h101i = lap2h(i1+1,i2,i3+1,0) - lap2h(i1-1,i2,i3+1,0) - lap2h(i1+1,i2,i3-1,0) + lap2h(i1-1,i2,i3-1,0)
                                      lap2h011i = lap2h(i1,i2+1,i3+1,0) - lap2h(i1,i2-1,i3+1,0) - lap2h(i1,i2+1,i3-1,0) + lap2h(i1,i2-1,i3-1,0)
                                      lap2hSq(i1,i2,i3,0) =  lapCoeff(i1,i2,i3,0)*lap2h200(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*lap2h020(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*lap2h002(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*lap2h110i  + lapCoeff(i1,i2,i3,4)*lap2h101i  + lapCoeff(i1,i2,i3,5)*lap2h011i  + lapCoeff(i1,i2,i3,6)*lap2h100i  + lapCoeff(i1,i2,i3,7)*lap2h010i  + lapCoeff(i1,i2,i3,8)*lap2h001i    
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
                                      d600(i1,i2,i3,0) = d400(i1+1,i2,i3,0) - 2*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0)
                                      d060(i1,i2,i3,0) = d040(i1,i2+1,i3,0) - 2*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0)
                                      d006(i1,i2,i3,0) = d004(i1,i2,i3+1,0) - 2*d004(i1,i2,i3,0) + d004(i1,i2,i3-1,0)
                                      d420(i1,i2,i3,0) = d220(i1+1,i2,i3,0) - 2*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0)
                                      d240(i1,i2,i3,0) = d040(i1+1,i2,i3,0) - 2*d040(i1,i2,i3,0) + d040(i1-1,i2,i3,0)
                                      d402(i1,i2,i3,0) = d202(i1+1,i2,i3,0) - 2*d202(i1,i2,i3,0) + d202(i1-1,i2,i3,0)
                                      d204(i1,i2,i3,0) = d004(i1+1,i2,i3,0) - 2*d004(i1,i2,i3,0) + d004(i1-1,i2,i3,0)
                                      d042(i1,i2,i3,0) = d022(i1,i2+1,i3,0) - 2*d022(i1,i2,i3,0) + d022(i1,i2-1,i3,0)
                                      d024(i1,i2,i3,0) = d004(i1,i2+1,i3,0) - 2*d004(i1,i2,i3,0) + d004(i1,i2-1,i3,0)
                                      d500i = d400(i1+1,i2,i3,0) - d400(i1-1,i2,i3,0)
                                      d050i = d040(i1,i2+1,i3,0) - d040(i1,i2-1,i3,0)
                                      d005i = d004(i1,i2,i3+1,0) - d004(i1,i2,i3-1,0)
                                      d510i = d400(i1+1,i2+1,i3,0) - d400(i1-1,i2+1,i3,0) - d400(i1+1,i2-1,i3,0) + d400(i1-1,i2-1,i3,0)
                                      d150i = d040(i1+1,i2+1,i3,0) - d040(i1-1,i2+1,i3,0) - d040(i1+1,i2-1,i3,0) + d040(i1-1,i2-1,i3,0)
                                      d330i = d220(i1+1,i2+1,i3,0) - d220(i1-1,i2+1,i3,0) - d220(i1+1,i2-1,i3,0) + d220(i1-1,i2-1,i3,0)
                                      d501i = d400(i1+1,i2,i3+1,0) - d400(i1-1,i2,i3+1,0) - d400(i1+1,i2,i3-1,0) + d400(i1-1,i2,i3-1,0)
                                      d105i = d004(i1+1,i2,i3+1,0) - d004(i1-1,i2,i3+1,0) - d004(i1+1,i2,i3-1,0) + d004(i1-1,i2,i3-1,0)
                                      d051i = d040(i1,i2+1,i3+1,0) - d040(i1,i2-1,i3+1,0) - d040(i1,i2+1,i3-1,0) + d040(i1,i2-1,i3-1,0)
                                      d015i = d004(i1,i2+1,i3+1,0) - d004(i1,i2-1,i3+1,0) - d004(i1,i2+1,i3-1,0) + d004(i1,i2-1,i3-1,0)
                                      d303i = d202(i1+1,i2,i3+1,0) - d202(i1-1,i2,i3+1,0) - d202(i1+1,i2,i3-1,0) + d202(i1-1,i2,i3-1,0)
                                      d033i = d022(i1,i2+1,i3+1,0) - d022(i1,i2-1,i3+1,0) - d022(i1,i2+1,i3-1,0) + d022(i1,i2-1,i3-1,0)
                   ! --- Laplacian to order 6 = lap4h + corrections 
                                      lap6h(i1,i2,i3,0) = lap4h(i1,i2,i3,0) + lapCoeff(i1,i2,i3,0)*crr2*d600(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*css2*d060(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*ctt2*d006(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*(cr2*d510i + cs2*d150i + cr1*cs1*d330i ) + lapCoeff(i1,i2,i3,4)*(cr2*d501i + ct2*d105i + cr1*ct1*d303i ) + lapCoeff(i1,i2,i3,5)*(cs2*d051i + ct2*d015i + cs1*ct1*d033i ) + lapCoeff(i1,i2,i3,6)*cr2 *d500i + lapCoeff(i1,i2,i3,7)*cs2 *d050i + lapCoeff(i1,i2,i3,8)*ct2 *d005i 
                                      lap2hSq200(i1,i2,i3,0) = lap2hSq(i1+1,i2,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1-1,i2,i3,0)
                                      lap2hSq020(i1,i2,i3,0) = lap2hSq(i1,i2+1,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2-1,i3,0)
                                      lap2hSq002(i1,i2,i3,0) = lap2hSq(i1,i2,i3+1,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2,i3-1,0)
                                      lap2hSq100i = lap2hSq(i1+1,i2,i3,0) - lap2hSq(i1-1,i2,i3,0)
                                      lap2hSq010i = lap2hSq(i1,i2+1,i3,0) - lap2hSq(i1,i2-1,i3,0)
                                      lap2hSq001i = lap2hSq(i1,i2,i3+1,0) - lap2hSq(i1,i2,i3-1,0)
                                      lap2hSq110i = lap2hSq(i1+1,i2+1,i3,0) - lap2hSq(i1-1,i2+1,i3,0) - lap2hSq(i1+1,i2-1,i3,0) + lap2hSq(i1-1,i2-1,i3,0)
                                      lap2hSq101i = lap2hSq(i1+1,i2,i3+1,0) - lap2hSq(i1-1,i2,i3+1,0) - lap2hSq(i1+1,i2,i3-1,0) + lap2hSq(i1-1,i2,i3-1,0)
                                      lap2hSq011i = lap2hSq(i1,i2+1,i3+1,0) - lap2hSq(i1,i2-1,i3+1,0) - lap2hSq(i1,i2+1,i3-1,0) + lap2hSq(i1,i2-1,i3-1,0)
                                      lap2hCubed(i1,i2,i3,0) =  + lapCoeff(i1,i2,i3,0)*lap2hSq200(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*lap2hSq020(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*lap2hSq002(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*lap2hSq110i  + lapCoeff(i1,i2,i3,4)*lap2hSq101i  + lapCoeff(i1,i2,i3,5)*lap2hSq011i  + lapCoeff(i1,i2,i3,6)*lap2hSq100i  + lapCoeff(i1,i2,i3,7)*lap2hSq010i  + lapCoeff(i1,i2,i3,8)*lap2hSq001i   
                                    end do
                                    end do
                                    end do
                ! --- SPLIT LOOPS
                                do i3=n3a,n3b
                                do i2=n2a,n2b
                                do i1=n1a,n1b
                   ! --- Laplacian squared to order 4 = 
                   !  lap2h*( lap4h ) + corrections*( Lap2h )
                                      lap4h200(i1,i2,i3,0) = lap4h(i1+1,i2,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1-1,i2,i3,0)
                                      lap4h020(i1,i2,i3,0) = lap4h(i1,i2+1,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2-1,i3,0)
                                      lap4h002(i1,i2,i3,0) = lap4h(i1,i2,i3+1,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2,i3-1,0)
                                      lap2h400(i1,i2,i3,0) = lap2h200(i1+1,i2,i3,0) - 2*lap2h200(i1,i2,i3,0) + lap2h200(i1-1,i2,i3,0)
                                      lap2h040(i1,i2,i3,0) = lap2h020(i1,i2+1,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1,i2-1,i3,0)
                                      lap2h004(i1,i2,i3,0) = lap2h002(i1,i2,i3+1,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1,i2,i3-1,0)
                                      lap2h220(i1,i2,i3,0) = lap2h020(i1+1,i2,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1-1,i2,i3,0)
                                      lap2h202(i1,i2,i3,0) = lap2h002(i1+1,i2,i3,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1-1,i2,i3,0)
                                      lap2h022(i1,i2,i3,0) = lap2h002(i1,i2+1,i3,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1,i2-1,i3,0)
                                      lap4h100i = lap4h(i1+1,i2,i3,0) - lap4h(i1-1,i2,i3,0)
                                      lap4h010i = lap4h(i1,i2+1,i3,0) - lap4h(i1,i2-1,i3,0)
                                      lap4h001i = lap4h(i1,i2,i3+1,0) - lap4h(i1,i2,i3-1,0)
                                      lap4h110i = lap4h(i1+1,i2+1,i3,0) - lap4h(i1-1,i2+1,i3,0) - lap4h(i1+1,i2-1,i3,0) + lap4h(i1-1,i2-1,i3,0)
                                      lap4h101i = lap4h(i1+1,i2,i3+1,0) - lap4h(i1-1,i2,i3+1,0) - lap4h(i1+1,i2,i3-1,0) + lap4h(i1-1,i2,i3-1,0)
                                      lap4h011i = lap4h(i1,i2+1,i3+1,0) - lap4h(i1,i2-1,i3+1,0) - lap4h(i1,i2+1,i3-1,0) + lap4h(i1,i2-1,i3-1,0)
                                      lap2h300i = lap2h200(i1+1,i2,i3,0) - lap2h200(i1-1,i2,i3,0)
                                      lap2h030i = lap2h020(i1,i2+1,i3,0) - lap2h020(i1,i2-1,i3,0)
                                      lap2h003i = lap2h002(i1,i2,i3+1,0) - lap2h002(i1,i2,i3-1,0)
                                      lap2h310i = lap2h200(i1+1,i2+1,i3,0) - lap2h200(i1-1,i2+1,i3,0) - lap2h200(i1+1,i2-1,i3,0) + lap2h200(i1-1,i2-1,i3,0)
                                      lap2h130i = lap2h020(i1+1,i2+1,i3,0) - lap2h020(i1-1,i2+1,i3,0) - lap2h020(i1+1,i2-1,i3,0) + lap2h020(i1-1,i2-1,i3,0)
                                      lap2h301i = lap2h200(i1+1,i2,i3+1,0) - lap2h200(i1-1,i2,i3+1,0) - lap2h200(i1+1,i2,i3-1,0) + lap2h200(i1-1,i2,i3-1,0)
                                      lap2h103i = lap2h002(i1+1,i2,i3+1,0) - lap2h002(i1-1,i2,i3+1,0) - lap2h002(i1+1,i2,i3-1,0) + lap2h002(i1-1,i2,i3-1,0)
                                      lap2h031i = lap2h020(i1,i2+1,i3+1,0) - lap2h020(i1,i2-1,i3+1,0) - lap2h020(i1,i2+1,i3-1,0) + lap2h020(i1,i2-1,i3-1,0)
                                      lap2h013i = lap2h002(i1,i2+1,i3+1,0) - lap2h002(i1,i2-1,i3+1,0) - lap2h002(i1,i2+1,i3-1,0) + lap2h002(i1,i2-1,i3-1,0)
                                      lap4hSq(i1,i2,i3,0) =     lapCoeff(i1,i2,i3,0)*( lap4h200(i1,i2,i3,0) + crr1*lap2h400(i1,i2,i3,0) )    + lapCoeff(i1,i2,i3,1)*( lap4h020(i1,i2,i3,0) + css1*lap2h040(i1,i2,i3,0) )     + lapCoeff(i1,i2,i3,2)*( lap4h002(i1,i2,i3,0) + ctt1*lap2h004(i1,i2,i3,0) )     + lapCoeff(i1,i2,i3,3)*( lap4h110i + cr1*lap2h310i + cs1*lap2h130i ) + lapCoeff(i1,i2,i3,4)*( lap4h101i + cr1*lap2h301i + ct1*lap2h103i ) + lapCoeff(i1,i2,i3,5)*( lap4h011i + cs1*lap2h031i + ct1*lap2h013i ) + lapCoeff(i1,i2,i3,6)*( lap4h100i + cr1 *lap2h300i )    + lapCoeff(i1,i2,i3,7)*( lap4h010i + cs1 *lap2h030i )    + lapCoeff(i1,i2,i3,8)*( lap4h001i + ct1 *lap2h003i )      
                                    end do
                                    end do
                                    end do
               ! ===========  FINAL LOOP TO FILL IN THE SOLUTION ============
                              numGhost1=0;
                              n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                              n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                              n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                do i3=n3a,n3b
                                do i2=n2a,n2b
                                do i1=n1a,n1b
                                      d800i = d600(i1+1,i2,i3,0) - 2*d600(i1,i2,i3,0) + d600(i1-1,i2,i3,0)
                                      d080i = d060(i1,i2+1,i3,0) - 2*d060(i1,i2,i3,0) + d060(i1,i2-1,i3,0)
                                      d008i = d006(i1,i2,i3+1,0) - 2*d006(i1,i2,i3,0) + d006(i1,i2,i3-1,0)
                                      d700i = d600(i1+1,i2,i3,0) - d600(i1-1,i2,i3,0)
                                      d070i = d060(i1,i2+1,i3,0) - d060(i1,i2-1,i3,0)
                                      d007i = d006(i1,i2,i3+1,0) - d006(i1,i2,i3-1,0)
                                      d710i = d600(i1+1,i2+1,i3,0) - d600(i1-1,i2+1,i3,0) - d600(i1+1,i2-1,i3,0) + d600(i1-1,i2-1,i3,0)
                                      d170i = d060(i1+1,i2+1,i3,0) - d060(i1-1,i2+1,i3,0) - d060(i1+1,i2-1,i3,0) + d060(i1-1,i2-1,i3,0)
                                      d701i = d600(i1+1,i2,i3+1,0) - d600(i1-1,i2,i3+1,0) - d600(i1+1,i2,i3-1,0) + d600(i1-1,i2,i3-1,0)
                                      d107i = d006(i1+1,i2,i3+1,0) - d006(i1-1,i2,i3+1,0) - d006(i1+1,i2,i3-1,0) + d006(i1-1,i2,i3-1,0)
                                      d071i = d060(i1,i2+1,i3+1,0) - d060(i1,i2-1,i3+1,0) - d060(i1,i2+1,i3-1,0) + d060(i1,i2-1,i3-1,0)
                                      d017i = d006(i1,i2+1,i3+1,0) - d006(i1,i2-1,i3+1,0) - d006(i1,i2+1,i3-1,0) + d006(i1,i2-1,i3-1,0)
                                      d530i = d420(i1+1,i2+1,i3,0) - d420(i1-1,i2+1,i3,0) - d420(i1+1,i2-1,i3,0) + d420(i1-1,i2-1,i3,0)
                                      d350i = d240(i1+1,i2+1,i3,0) - d240(i1-1,i2+1,i3,0) - d240(i1+1,i2-1,i3,0) + d240(i1-1,i2-1,i3,0)
                                      d503i = d402(i1+1,i2,i3+1,0) - d402(i1-1,i2,i3+1,0) - d402(i1+1,i2,i3-1,0) + d402(i1-1,i2,i3-1,0)
                                      d305i = d204(i1+1,i2,i3+1,0) - d204(i1-1,i2,i3+1,0) - d204(i1+1,i2,i3-1,0) + d204(i1-1,i2,i3-1,0)
                                      d053i = d042(i1,i2+1,i3+1,0) - d042(i1,i2-1,i3+1,0) - d042(i1,i2+1,i3-1,0) + d042(i1,i2-1,i3-1,0)
                                      d035i = d024(i1,i2+1,i3+1,0) - d024(i1,i2-1,i3+1,0) - d024(i1,i2+1,i3-1,0) + d024(i1,i2-1,i3-1,0)
                   ! --- Laplacian to order 8 = lap6h + corrections 
                                      lap8h = lap6h(i1,i2,i3,0)                                                         + lapCoeff(i1,i2,i3,0)*crr3*d800i                                               + lapCoeff(i1,i2,i3,1)*css3*d080i                                               + lapCoeff(i1,i2,i3,2)*ctt3*d008i                                               + lapCoeff(i1,i2,i3,3)*(cr3*d710i + cs3*d170i + cr2*cs1*d530i + cr1*cs2*d350i ) + lapCoeff(i1,i2,i3,4)*(cr3*d701i + ct3*d107i + cr2*ct1*d503i + cr1*ct2*d305i ) + lapCoeff(i1,i2,i3,5)*(cs3*d071i + ct3*d017i + cs2*ct1*d053i + cs1*ct2*d035i ) + lapCoeff(i1,i2,i3,6)* cr3*d700i                                               + lapCoeff(i1,i2,i3,7)* cs3*d070i                                               + lapCoeff(i1,i2,i3,8)* ct3*d007i 
                   ! --- Laplacian^4 4p (4th power) order 2: 
                                      lap2hCubed200i = lap2hCubed(i1+1,i2,i3,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1-1,i2,i3,0)
                                      lap2hCubed020i = lap2hCubed(i1,i2+1,i3,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1,i2-1,i3,0)
                                      lap2hCubed002i = lap2hCubed(i1,i2,i3+1,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1,i2,i3-1,0)
                                      lap2hCubed100i = lap2hCubed(i1+1,i2,i3,0) - lap2hCubed(i1-1,i2,i3,0)
                                      lap2hCubed010i = lap2hCubed(i1,i2+1,i3,0) - lap2hCubed(i1,i2-1,i3,0)
                                      lap2hCubed001i = lap2hCubed(i1,i2,i3+1,0) - lap2hCubed(i1,i2,i3-1,0)
                                      lap2hCubed110i = lap2hCubed(i1+1,i2+1,i3,0) - lap2hCubed(i1-1,i2+1,i3,0) - lap2hCubed(i1+1,i2-1,i3,0) + lap2hCubed(i1-1,i2-1,i3,0)
                                      lap2hCubed101i = lap2hCubed(i1+1,i2,i3+1,0) - lap2hCubed(i1-1,i2,i3+1,0) - lap2hCubed(i1+1,i2,i3-1,0) + lap2hCubed(i1-1,i2,i3-1,0)
                                      lap2hCubed011i = lap2hCubed(i1,i2+1,i3+1,0) - lap2hCubed(i1,i2-1,i3+1,0) - lap2hCubed(i1,i2+1,i3-1,0) + lap2hCubed(i1,i2-1,i3-1,0)
                                      lap2h4p  =                             + lapCoeff(i1,i2,i3,0)*lap2hCubed200i  + lapCoeff(i1,i2,i3,1)*lap2hCubed020i  + lapCoeff(i1,i2,i3,2)*lap2hCubed002i  + lapCoeff(i1,i2,i3,3)*lap2hCubed110i  + lapCoeff(i1,i2,i3,4)*lap2hCubed101i  + lapCoeff(i1,i2,i3,5)*lap2hCubed011i  + lapCoeff(i1,i2,i3,6)*lap2hCubed100i  + lapCoeff(i1,i2,i3,7)*lap2hCubed010i  + lapCoeff(i1,i2,i3,8)*lap2hCubed001i    
                   ! --- Laplacian squared to order 6 :
                   !   Lap6h = Lap4h + M4  = (Lap2h) + M2 + M4 
                   !   Lap6h*Lap6h = [ (Lap2h) + M2 + M4 ] [ (Lap2h) + M2 + M4 ]
                   !               = Lap2h*Lap6h + M2*Lap4h + M4*Lap2h + O(h^6)
                                      lap6h200i = lap6h(i1+1,i2,i3,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1-1,i2,i3,0)
                                      lap6h020i = lap6h(i1,i2+1,i3,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1,i2-1,i3,0)
                                      lap6h002i = lap6h(i1,i2,i3+1,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1,i2,i3-1,0)
                                      lap6h100i = lap6h(i1+1,i2,i3,0) - lap6h(i1-1,i2,i3,0)
                                      lap6h010i = lap6h(i1,i2+1,i3,0) - lap6h(i1,i2-1,i3,0)
                                      lap6h001i = lap6h(i1,i2,i3+1,0) - lap6h(i1,i2,i3-1,0)
                                      lap6h110i = lap6h(i1+1,i2+1,i3,0) - lap6h(i1-1,i2+1,i3,0) - lap6h(i1+1,i2-1,i3,0) + lap6h(i1-1,i2-1,i3,0)
                                      lap6h101i = lap6h(i1+1,i2,i3+1,0) - lap6h(i1-1,i2,i3+1,0) - lap6h(i1+1,i2,i3-1,0) + lap6h(i1-1,i2,i3-1,0)
                                      lap6h011i = lap6h(i1,i2+1,i3+1,0) - lap6h(i1,i2-1,i3+1,0) - lap6h(i1,i2+1,i3-1,0) + lap6h(i1,i2-1,i3-1,0)
                                      lap4h400i = lap4h200(i1+1,i2,i3,0) - 2*lap4h200(i1,i2,i3,0) + lap4h200(i1-1,i2,i3,0)
                                      lap4h040i = lap4h020(i1,i2+1,i3,0) - 2*lap4h020(i1,i2,i3,0) + lap4h020(i1,i2-1,i3,0)
                                      lap4h004i = lap4h002(i1,i2,i3+1,0) - 2*lap4h002(i1,i2,i3,0) + lap4h002(i1,i2,i3-1,0)
                                      lap4h300i = lap4h200(i1+1,i2,i3,0) - lap4h200(i1-1,i2,i3,0)
                                      lap4h030i = lap4h020(i1,i2+1,i3,0) - lap4h020(i1,i2-1,i3,0)
                                      lap4h003i = lap4h002(i1,i2,i3+1,0) - lap4h002(i1,i2,i3-1,0)
                                      lap4h310i = lap4h200(i1+1,i2+1,i3,0) - lap4h200(i1-1,i2+1,i3,0) - lap4h200(i1+1,i2-1,i3,0) + lap4h200(i1-1,i2-1,i3,0)
                                      lap4h130i = lap4h020(i1+1,i2+1,i3,0) - lap4h020(i1-1,i2+1,i3,0) - lap4h020(i1+1,i2-1,i3,0) + lap4h020(i1-1,i2-1,i3,0)
                                      lap4h301i = lap4h200(i1+1,i2,i3+1,0) - lap4h200(i1-1,i2,i3+1,0) - lap4h200(i1+1,i2,i3-1,0) + lap4h200(i1-1,i2,i3-1,0)
                                      lap4h103i = lap4h002(i1+1,i2,i3+1,0) - lap4h002(i1-1,i2,i3+1,0) - lap4h002(i1+1,i2,i3-1,0) + lap4h002(i1-1,i2,i3-1,0)
                                      lap4h031i = lap4h020(i1,i2+1,i3+1,0) - lap4h020(i1,i2-1,i3+1,0) - lap4h020(i1,i2+1,i3-1,0) + lap4h020(i1,i2-1,i3-1,0)
                                      lap4h013i = lap4h002(i1,i2+1,i3+1,0) - lap4h002(i1,i2-1,i3+1,0) - lap4h002(i1,i2+1,i3-1,0) + lap4h002(i1,i2-1,i3-1,0)
                                      lap2h600i = lap2h400(i1+1,i2,i3,0) - 2*lap2h400(i1,i2,i3,0) + lap2h400(i1-1,i2,i3,0)
                                      lap2h060i = lap2h040(i1,i2+1,i3,0) - 2*lap2h040(i1,i2,i3,0) + lap2h040(i1,i2-1,i3,0)
                                      lap2h006i = lap2h004(i1,i2,i3+1,0) - 2*lap2h004(i1,i2,i3,0) + lap2h004(i1,i2,i3-1,0)
                                      lap2h500i = lap2h400(i1+1,i2,i3,0) - lap2h400(i1-1,i2,i3,0)
                                      lap2h050i = lap2h040(i1,i2+1,i3,0) - lap2h040(i1,i2-1,i3,0)
                                      lap2h005i = lap2h004(i1,i2,i3+1,0) - lap2h004(i1,i2,i3-1,0)
                                      lap2h510i = lap2h400(i1+1,i2+1,i3,0) - lap2h400(i1-1,i2+1,i3,0) - lap2h400(i1+1,i2-1,i3,0) + lap2h400(i1-1,i2-1,i3,0)
                                      lap2h150i = lap2h040(i1+1,i2+1,i3,0) - lap2h040(i1-1,i2+1,i3,0) - lap2h040(i1+1,i2-1,i3,0) + lap2h040(i1-1,i2-1,i3,0)
                                      lap2h330i = lap2h220(i1+1,i2+1,i3,0) - lap2h220(i1-1,i2+1,i3,0) - lap2h220(i1+1,i2-1,i3,0) + lap2h220(i1-1,i2-1,i3,0)
                                      lap2h501i = lap2h400(i1+1,i2,i3+1,0) - lap2h400(i1-1,i2,i3+1,0) - lap2h400(i1+1,i2,i3-1,0) + lap2h400(i1-1,i2,i3-1,0)
                                      lap2h105i = lap2h004(i1+1,i2,i3+1,0) - lap2h004(i1-1,i2,i3+1,0) - lap2h004(i1+1,i2,i3-1,0) + lap2h004(i1-1,i2,i3-1,0)
                                      lap2h051i = lap2h040(i1,i2+1,i3+1,0) - lap2h040(i1,i2-1,i3+1,0) - lap2h040(i1,i2+1,i3-1,0) + lap2h040(i1,i2-1,i3-1,0)
                                      lap2h015i = lap2h004(i1,i2+1,i3+1,0) - lap2h004(i1,i2-1,i3+1,0) - lap2h004(i1,i2+1,i3-1,0) + lap2h004(i1,i2-1,i3-1,0)
                                      lap2h303i = lap2h202(i1+1,i2,i3+1,0) - lap2h202(i1-1,i2,i3+1,0) - lap2h202(i1+1,i2,i3-1,0) + lap2h202(i1-1,i2,i3-1,0)
                                      lap2h033i = lap2h022(i1,i2+1,i3+1,0) - lap2h022(i1,i2-1,i3+1,0) - lap2h022(i1,i2+1,i3-1,0) + lap2h022(i1,i2-1,i3-1,0)
                                      lap6hSq =                                                                                     lapCoeff(i1,i2,i3,0)*(lap6h200i + crr1*lap4h400i + crr2*lap2h600i )                       + lapCoeff(i1,i2,i3,1)*(lap6h020i + css1*lap4h040i + css2*lap2h060i )                       + lapCoeff(i1,i2,i3,2)*(lap6h002i + ctt1*lap4h004i + ctt2*lap2h006i )                       + lapCoeff(i1,i2,i3,3)*(lap6h110i +  cr1*lap4h310i +  cr2*lap2h510i                         +  cs1*lap4h130i +  cs2*lap2h150i + cr1*cs1*lap2h330i )   + lapCoeff(i1,i2,i3,4)*(lap6h101i +  cr1*lap4h301i +  cr2*lap2h501i                         +  ct1*lap4h103i +  ct2*lap2h105i + cr1*ct1*lap2h303i )   + lapCoeff(i1,i2,i3,5)*(lap6h011i +  cs1*lap4h031i +  cs2*lap2h051i                         +  ct1*lap4h013i +  ct2*lap2h015i + cs1*ct1*lap2h033i )   + lapCoeff(i1,i2,i3,6)*(lap6h100i +  cr1*lap4h300i +  cr2*lap2h500i )                       + lapCoeff(i1,i2,i3,7)*(lap6h010i +  cs1*lap4h030i +  cs2*lap2h050i )                       + lapCoeff(i1,i2,i3,8)*(lap6h010i +  ct1*lap4h003i +  ct2*lap2h005i )                         
                   ! --- Laplacian CUBED to order 4 
                                      lap4hSq200i = lap4hSq(i1+1,i2,i3,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1-1,i2,i3,0)
                                      lap4hSq020i = lap4hSq(i1,i2+1,i3,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1,i2-1,i3,0)
                                      lap4hSq002i = lap4hSq(i1,i2,i3+1,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1,i2,i3-1,0)
                                      lap4hSq100i = lap4hSq(i1+1,i2,i3,0) - lap4hSq(i1-1,i2,i3,0)
                                      lap4hSq010i = lap4hSq(i1,i2+1,i3,0) - lap4hSq(i1,i2-1,i3,0)
                                      lap4hSq001i = lap4hSq(i1,i2,i3+1,0) - lap4hSq(i1,i2,i3-1,0)
                                      lap4hSq110i = lap4hSq(i1+1,i2+1,i3,0) - lap4hSq(i1-1,i2+1,i3,0) - lap4hSq(i1+1,i2-1,i3,0) + lap4hSq(i1-1,i2-1,i3,0)
                                      lap4hSq101i = lap4hSq(i1+1,i2,i3+1,0) - lap4hSq(i1-1,i2,i3+1,0) - lap4hSq(i1+1,i2,i3-1,0) + lap4hSq(i1-1,i2,i3-1,0)
                                      lap4hSq011i = lap4hSq(i1,i2+1,i3+1,0) - lap4hSq(i1,i2-1,i3+1,0) - lap4hSq(i1,i2+1,i3-1,0) + lap4hSq(i1,i2-1,i3-1,0)
                                      lap2hSq400i = lap2hSq200(i1+1,i2,i3,0) - 2*lap2hSq200(i1,i2,i3,0) + lap2hSq200(i1-1,i2,i3,0)
                                      lap2hSq040i = lap2hSq020(i1,i2+1,i3,0) - 2*lap2hSq020(i1,i2,i3,0) + lap2hSq020(i1,i2-1,i3,0)
                                      lap2hSq004i = lap2hSq002(i1,i2,i3+1,0) - 2*lap2hSq002(i1,i2,i3,0) + lap2hSq002(i1,i2,i3-1,0)
                                      lap2hSq300i = lap2hSq200(i1+1,i2,i3,0) - lap2hSq200(i1-1,i2,i3,0)
                                      lap2hSq030i = lap2hSq020(i1,i2+1,i3,0) - lap2hSq020(i1,i2-1,i3,0)
                                      lap2hSq003i = lap2hSq002(i1,i2,i3+1,0) - lap2hSq002(i1,i2,i3-1,0)
                                      lap2hSq310i = lap2hSq200(i1+1,i2+1,i3,0) - lap2hSq200(i1-1,i2+1,i3,0) - lap2hSq200(i1+1,i2-1,i3,0) + lap2hSq200(i1-1,i2-1,i3,0)
                                      lap2hSq130i = lap2hSq020(i1+1,i2+1,i3,0) - lap2hSq020(i1-1,i2+1,i3,0) - lap2hSq020(i1+1,i2-1,i3,0) + lap2hSq020(i1-1,i2-1,i3,0)
                                      lap2hSq301i = lap2hSq200(i1+1,i2,i3+1,0) - lap2hSq200(i1-1,i2,i3+1,0) - lap2hSq200(i1+1,i2,i3-1,0) + lap2hSq200(i1-1,i2,i3-1,0)
                                      lap2hSq103i = lap2hSq002(i1+1,i2,i3+1,0) - lap2hSq002(i1-1,i2,i3+1,0) - lap2hSq002(i1+1,i2,i3-1,0) + lap2hSq002(i1-1,i2,i3-1,0)
                                      lap2hSq031i = lap2hSq020(i1,i2+1,i3+1,0) - lap2hSq020(i1,i2-1,i3+1,0) - lap2hSq020(i1,i2+1,i3-1,0) + lap2hSq020(i1,i2-1,i3-1,0)
                                      lap2hSq013i = lap2hSq002(i1,i2+1,i3+1,0) - lap2hSq002(i1,i2-1,i3+1,0) - lap2hSq002(i1,i2+1,i3-1,0) + lap2hSq002(i1,i2-1,i3-1,0)
                                      lap4hCubed =                                                    lapCoeff(i1,i2,i3,0)*(lap4hSq200i + crr1*lap2hSq400i )   + lapCoeff(i1,i2,i3,1)*(lap4hSq020i + css1*lap2hSq040i )   + lapCoeff(i1,i2,i3,2)*(lap4hSq002i + ctt1*lap2hSq004i )   + lapCoeff(i1,i2,i3,3)*(lap4hSq110i +  cr1*lap2hSq310i     +  cs1*lap2hSq130i )   + lapCoeff(i1,i2,i3,4)*(lap4hSq101i +  cr1*lap2hSq301i     +  ct1*lap2hSq103i )   + lapCoeff(i1,i2,i3,5)*(lap4hSq011i +  cs1*lap2hSq031i     +  ct1*lap2hSq013i )   + lapCoeff(i1,i2,i3,6)*(lap4hSq100i + cr1 *lap2hSq300i )   + lapCoeff(i1,i2,i3,7)*(lap4hSq010i + cs1 *lap2hSq030i )   + lapCoeff(i1,i2,i3,8)*(lap4hSq001i + ct1 *lap2hSq003i )     
                   ! --- Modified equation space-time update ----
                                      un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m)  + cdtsq*( lap8h )               + cdtPow4By12*( lap6hSq )       + cdtPow6By360*( lap4hCubed )   + cdtPow8By20160*( lap2h4p )    +dtSq*fv(m)                    
                                    end do
                                    end do
                                    end do
                      else
               ! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
               ! Example: 
               ! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...
cr0 = 1.; cr1 = -1/6.; cr2 = 1/30.; cr3 = -1/140.; 
cs0 = 1.; cs1 = -1/6.; cs2 = 1/30.; cs3 = -1/140.; 
ct0 = 1.; ct1 = -1/6.; ct2 = 1/30.; ct3 = -1/140.; 
crr0 = 1.; crr1 = -1/12.; crr2 = 1/90.; crr3 = -1/560.; 
css0 = 1.; css1 = -1/12.; css2 = 1/90.; css3 = -1/560.; 
ctt0 = 1.; ctt1 = -1/12.; ctt2 = 1/90.; ctt3 = -1/560.; 
crrr0 = 1.; crrr1 = -1/4.; crrr2 = 7/120.; crrr3 = -41/3024.; 
csss0 = 1.; csss1 = -1/4.; csss2 = 7/120.; csss3 = -41/3024.; 
cttt0 = 1.; cttt1 = -1/4.; cttt2 = 7/120.; cttt3 = -41/3024.; 
crrrr0 = 1.; crrrr1 = -1/6.; crrrr2 = 7/240.; crrrr3 = -41/7560.; 
cssss0 = 1.; cssss1 = -1/6.; cssss2 = 7/240.; cssss3 = -41/7560.; 
ctttt0 = 1.; ctttt1 = -1/6.; ctttt2 = 7/240.; ctttt3 = -41/7560.; 
crrrrr0 = 1.; crrrrr1 = -1/3.; crrrrr2 = 13/144.; crrrrr3 = -139/6048.; 
csssss0 = 1.; csssss1 = -1/3.; csssss2 = 13/144.; csssss3 = -139/6048.; 
cttttt0 = 1.; cttttt1 = -1/3.; cttttt2 = 13/144.; cttttt3 = -139/6048.; 
crrrrrr0 = 1.; crrrrrr1 = -1/4.; crrrrrr2 = 13/240.; crrrrrr3 = -139/12096.; 
cssssss0 = 1.; cssssss1 = -1/4.; cssssss2 = 13/240.; cssssss3 = -139/12096.; 
ctttttt0 = 1.; ctttttt1 = -1/4.; ctttttt2 = 13/240.; ctttttt3 = -139/12096.; 
                              dr1=dr(0); dr1i=1./dr1;
                              dr2=dr(1); dr2i=1./dr2;
                              dr3=dr(2); dr3i=1./dr3;
                              fv(m)=0.
                              if( lapCoeff(0,0,0,0).le.0. )then
                 ! --- Evaluate and store coefficients in Laplacian ---
                                  write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'
                                  numGhost1=3;
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
                                      if( (i1-4).ge.nd1a .and. (i1+4).le.nd1b )then
                                          diffOrder1=8
                                      elseif( (i1-3).ge.nd1a .and. (i1+3).le.nd1b )then
                                          diffOrder1=6
                                      elseif( (i1-2).ge.nd1a .and. (i1+2).le.nd1b )then
                                          diffOrder1=4
                                      elseif( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
                                          diffOrder1=2
                                      else
                                          stop 999
                                      end if
                                      if( (i2-4).ge.nd2a .and. (i2+4).le.nd2b )then
                                          diffOrder2=8
                                      elseif( (i2-3).ge.nd2a .and. (i2+3).le.nd2b )then
                                          diffOrder2=6
                                      elseif( (i2-2).ge.nd2a .and. (i2+2).le.nd2b )then
                                          diffOrder2=4
                                      elseif( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
                                          diffOrder2=2
                                      else
                                          stop 999
                                      end if
                                      if( (i3-4).ge.nd3a .and. (i3+4).le.nd3b )then
                                          diffOrder3=8
                                      elseif( (i3-3).ge.nd3a .and. (i3+3).le.nd3b )then
                                          diffOrder3=6
                                      elseif( (i3-2).ge.nd3a .and. (i3+2).le.nd3b )then
                                          diffOrder3=4
                                      elseif( (i3-1).ge.nd3a .and. (i3+1).le.nd3b )then
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
                                          rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
                                          rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
                                          rzs = ( 8*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) )*(dr2i/12.) 
                                          sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
                                          sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
                                          szs = ( 8*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) )*(dr2i/12.) 
                                          txs = ( 8*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) )*(dr2i/12.) 
                                          tys = ( 8*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) )*(dr2i/12.) 
                                          tzs = ( 8*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) )*(dr2i/12.) 
                                      elseif( diffOrder2.eq.8 )then
                                          rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
                                          rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
                                          rzs = ( 8*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) )*(dr2i/12.) 
                                          sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
                                          sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
                                          szs = ( 8*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) )*(dr2i/12.) 
                                          txs = ( 8*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) )*(dr2i/12.) 
                                          tys = ( 8*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) )*(dr2i/12.) 
                                          tzs = ( 8*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) )*(dr2i/12.) 
                                      end if
                                      if( diffOrder3.eq.2 )then
                                          rxt = (rsxy(i1,i2,i3+1,0,0)-rsxy(i1,i2,i3-1,0,0))*(.5*dr2i) 
                                          ryt = (rsxy(i1,i2,i3+1,0,1)-rsxy(i1,i2,i3-1,0,1))*(.5*dr2i) 
                                          rzt = (rsxy(i1,i2,i3+1,0,2)-rsxy(i1,i2,i3-1,0,2))*(.5*dr2i) 
                                          sxt = (rsxy(i1,i2,i3+1,1,0)-rsxy(i1,i2,i3-1,1,0))*(.5*dr2i) 
                                          syt = (rsxy(i1,i2,i3+1,1,1)-rsxy(i1,i2,i3-1,1,1))*(.5*dr2i) 
                                          szt = (rsxy(i1,i2,i3+1,1,2)-rsxy(i1,i2,i3-1,1,2))*(.5*dr2i) 
                                          txt = (rsxy(i1,i2,i3+1,2,0)-rsxy(i1,i2,i3-1,2,0))*(.5*dr2i) 
                                          tyt = (rsxy(i1,i2,i3+1,2,1)-rsxy(i1,i2,i3-1,2,1))*(.5*dr2i) 
                                          tzt = (rsxy(i1,i2,i3+1,2,2)-rsxy(i1,i2,i3-1,2,2))*(.5*dr2i) 
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
                                      lapCoeff(i1,i2,i3,3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i*.25
                                      lapCoeff(i1,i2,i3,4) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr2i*.25
                                      lapCoeff(i1,i2,i3,5) = 2.*(sx*tx + sy*ty + sz*tz )*dr1i*dr2i*.25
                                      lapCoeff(i1,i2,i3,6) = (rxx + ryy + rzz)*dr1i*.5
                                      lapCoeff(i1,i2,i3,7) = (sxx + syy + tyy)*dr2i*.5 
                                      lapCoeff(i1,i2,i3,8) = (txx + tyy + tzz)*dr3i*.5 
                                      end if ! mask .ne. 0
                                    end do
                                    end do
                                    end do
                              end if ! end assignLapCoeff
                              numGhost1=3;
                              n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                              n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                              n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                do i3=n3a,n3b
                                do i2=n2a,n2b
                                do i1=n1a,n1b
                                  if( mask(i1,i2,i3).ne.0 )then
                                      d200(i1,i2,i3,0) = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
                                      d020(i1,i2,i3,0) = u(i1,i2+1,i3,0) - 2*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
                                      d002(i1,i2,i3,0) = u(i1,i2,i3+1,0) - 2*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
                                      d100i = u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0)
                                      d010i = u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0)
                                      d110i = u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0)
                                      d001i = u(i1,i2,i3+1,0) - u(i1,i2,i3-1,0)
                                      d101i = u(i1+1,i2,i3+1,0) - u(i1-1,i2,i3+1,0) - u(i1+1,i2,i3-1,0) + u(i1-1,i2,i3-1,0)
                                      d011i = u(i1,i2+1,i3+1,0) - u(i1,i2-1,i3+1,0) - u(i1,i2+1,i3-1,0) + u(i1,i2-1,i3-1,0)
                                      lap2h(i1,i2,i3,0) = lapCoeff(i1,i2,i3,0)*d200(i1,i2,i3,0) +lapCoeff(i1,i2,i3,1)*d020(i1,i2,i3,0) +lapCoeff(i1,i2,i3,2)*d002(i1,i2,i3,0) +lapCoeff(i1,i2,i3,3)*d110i + lapCoeff(i1,i2,i3,4)*d101i + lapCoeff(i1,i2,i3,5)*d011i + lapCoeff(i1,i2,i3,6)*d100i + lapCoeff(i1,i2,i3,7)*d010i + lapCoeff(i1,i2,i3,8)*d001i
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
                                      d400(i1,i2,i3,0) = d200(i1+1,i2,i3,0) - 2*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0)
                                      d040(i1,i2,i3,0) = d020(i1,i2+1,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1,i2-1,i3,0)
                                      d004(i1,i2,i3,0) = d002(i1,i2,i3+1,0) - 2*d002(i1,i2,i3,0) + d002(i1,i2,i3-1,0)
                                      d220(i1,i2,i3,0) = d020(i1+1,i2,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1-1,i2,i3,0)
                                      d202(i1,i2,i3,0) = d002(i1+1,i2,i3,0) - 2*d002(i1,i2,i3,0) + d002(i1-1,i2,i3,0)
                                      d022(i1,i2,i3,0) = d002(i1,i2+1,i3,0) - 2*d002(i1,i2,i3,0) + d002(i1,i2-1,i3,0)
                                      d300i = d200(i1+1,i2,i3,0) - d200(i1-1,i2,i3,0)
                                      d030i = d020(i1,i2+1,i3,0) - d020(i1,i2-1,i3,0)
                                      d003i = d002(i1,i2,i3+1,0) - d002(i1,i2,i3-1,0)
                                      d310i = d200(i1+1,i2+1,i3,0) - d200(i1-1,i2+1,i3,0) - d200(i1+1,i2-1,i3,0) + d200(i1-1,i2-1,i3,0)
                                      d130i = d020(i1+1,i2+1,i3,0) - d020(i1-1,i2+1,i3,0) - d020(i1+1,i2-1,i3,0) + d020(i1-1,i2-1,i3,0)
                                      d301i = d200(i1+1,i2,i3+1,0) - d200(i1-1,i2,i3+1,0) - d200(i1+1,i2,i3-1,0) + d200(i1-1,i2,i3-1,0)
                                      d103i = d002(i1+1,i2,i3+1,0) - d002(i1-1,i2,i3+1,0) - d002(i1+1,i2,i3-1,0) + d002(i1-1,i2,i3-1,0)
                                      d031i = d020(i1,i2+1,i3+1,0) - d020(i1,i2-1,i3+1,0) - d020(i1,i2+1,i3-1,0) + d020(i1,i2-1,i3-1,0)
                                      d013i = d002(i1,i2+1,i3+1,0) - d002(i1,i2-1,i3+1,0) - d002(i1,i2+1,i3-1,0) + d002(i1,i2-1,i3-1,0)
                   ! --- Laplacian to order 4 = lap2h + corrections 
                                      lap4h(i1,i2,i3,0) = lap2h(i1,i2,i3,0) + lapCoeff(i1,i2,i3,0)*crr1*d400(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*css1*d040(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*ctt1*d004(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*(cr1*d310i + cs1*d130i) + lapCoeff(i1,i2,i3,4)*(cr1*d301i + ct1*d103i) + lapCoeff(i1,i2,i3,5)*(cs1*d031i + ct1*d013i) + lapCoeff(i1,i2,i3,6)*cr1 *d300i + lapCoeff(i1,i2,i3,7)*cs1 *d030i + lapCoeff(i1,i2,i3,8)*ct1 *d003i 
                   ! --- Laplacian squared to order 2:
                                      lap2h200(i1,i2,i3,0) = lap2h(i1+1,i2,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1-1,i2,i3,0)
                                      lap2h020(i1,i2,i3,0) = lap2h(i1,i2+1,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2-1,i3,0)
                                      lap2h002(i1,i2,i3,0) = lap2h(i1,i2,i3+1,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2,i3-1,0)
                                      lap2h100i = lap2h(i1+1,i2,i3,0) - lap2h(i1-1,i2,i3,0)
                                      lap2h010i = lap2h(i1,i2+1,i3,0) - lap2h(i1,i2-1,i3,0)
                                      lap2h110i = lap2h(i1+1,i2+1,i3,0) - lap2h(i1-1,i2+1,i3,0) - lap2h(i1+1,i2-1,i3,0) + lap2h(i1-1,i2-1,i3,0)
                                      lap2h001i = lap2h(i1,i2,i3+1,0) - lap2h(i1,i2,i3-1,0)
                                      lap2h101i = lap2h(i1+1,i2,i3+1,0) - lap2h(i1-1,i2,i3+1,0) - lap2h(i1+1,i2,i3-1,0) + lap2h(i1-1,i2,i3-1,0)
                                      lap2h011i = lap2h(i1,i2+1,i3+1,0) - lap2h(i1,i2-1,i3+1,0) - lap2h(i1,i2+1,i3-1,0) + lap2h(i1,i2-1,i3-1,0)
                                      lap2hSq(i1,i2,i3,0) =  lapCoeff(i1,i2,i3,0)*lap2h200(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*lap2h020(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*lap2h002(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*lap2h110i  + lapCoeff(i1,i2,i3,4)*lap2h101i  + lapCoeff(i1,i2,i3,5)*lap2h011i  + lapCoeff(i1,i2,i3,6)*lap2h100i  + lapCoeff(i1,i2,i3,7)*lap2h010i  + lapCoeff(i1,i2,i3,8)*lap2h001i    
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
                                      d600(i1,i2,i3,0) = d400(i1+1,i2,i3,0) - 2*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0)
                                      d060(i1,i2,i3,0) = d040(i1,i2+1,i3,0) - 2*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0)
                                      d006(i1,i2,i3,0) = d004(i1,i2,i3+1,0) - 2*d004(i1,i2,i3,0) + d004(i1,i2,i3-1,0)
                                      d420(i1,i2,i3,0) = d220(i1+1,i2,i3,0) - 2*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0)
                                      d240(i1,i2,i3,0) = d040(i1+1,i2,i3,0) - 2*d040(i1,i2,i3,0) + d040(i1-1,i2,i3,0)
                                      d402(i1,i2,i3,0) = d202(i1+1,i2,i3,0) - 2*d202(i1,i2,i3,0) + d202(i1-1,i2,i3,0)
                                      d204(i1,i2,i3,0) = d004(i1+1,i2,i3,0) - 2*d004(i1,i2,i3,0) + d004(i1-1,i2,i3,0)
                                      d042(i1,i2,i3,0) = d022(i1,i2+1,i3,0) - 2*d022(i1,i2,i3,0) + d022(i1,i2-1,i3,0)
                                      d024(i1,i2,i3,0) = d004(i1,i2+1,i3,0) - 2*d004(i1,i2,i3,0) + d004(i1,i2-1,i3,0)
                                      d500i = d400(i1+1,i2,i3,0) - d400(i1-1,i2,i3,0)
                                      d050i = d040(i1,i2+1,i3,0) - d040(i1,i2-1,i3,0)
                                      d005i = d004(i1,i2,i3+1,0) - d004(i1,i2,i3-1,0)
                                      d510i = d400(i1+1,i2+1,i3,0) - d400(i1-1,i2+1,i3,0) - d400(i1+1,i2-1,i3,0) + d400(i1-1,i2-1,i3,0)
                                      d150i = d040(i1+1,i2+1,i3,0) - d040(i1-1,i2+1,i3,0) - d040(i1+1,i2-1,i3,0) + d040(i1-1,i2-1,i3,0)
                                      d330i = d220(i1+1,i2+1,i3,0) - d220(i1-1,i2+1,i3,0) - d220(i1+1,i2-1,i3,0) + d220(i1-1,i2-1,i3,0)
                                      d501i = d400(i1+1,i2,i3+1,0) - d400(i1-1,i2,i3+1,0) - d400(i1+1,i2,i3-1,0) + d400(i1-1,i2,i3-1,0)
                                      d105i = d004(i1+1,i2,i3+1,0) - d004(i1-1,i2,i3+1,0) - d004(i1+1,i2,i3-1,0) + d004(i1-1,i2,i3-1,0)
                                      d051i = d040(i1,i2+1,i3+1,0) - d040(i1,i2-1,i3+1,0) - d040(i1,i2+1,i3-1,0) + d040(i1,i2-1,i3-1,0)
                                      d015i = d004(i1,i2+1,i3+1,0) - d004(i1,i2-1,i3+1,0) - d004(i1,i2+1,i3-1,0) + d004(i1,i2-1,i3-1,0)
                                      d303i = d202(i1+1,i2,i3+1,0) - d202(i1-1,i2,i3+1,0) - d202(i1+1,i2,i3-1,0) + d202(i1-1,i2,i3-1,0)
                                      d033i = d022(i1,i2+1,i3+1,0) - d022(i1,i2-1,i3+1,0) - d022(i1,i2+1,i3-1,0) + d022(i1,i2-1,i3-1,0)
                   ! --- Laplacian to order 6 = lap4h + corrections 
                                      lap6h(i1,i2,i3,0) = lap4h(i1,i2,i3,0) + lapCoeff(i1,i2,i3,0)*crr2*d600(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*css2*d060(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*ctt2*d006(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*(cr2*d510i + cs2*d150i + cr1*cs1*d330i ) + lapCoeff(i1,i2,i3,4)*(cr2*d501i + ct2*d105i + cr1*ct1*d303i ) + lapCoeff(i1,i2,i3,5)*(cs2*d051i + ct2*d015i + cs1*ct1*d033i ) + lapCoeff(i1,i2,i3,6)*cr2 *d500i + lapCoeff(i1,i2,i3,7)*cs2 *d050i + lapCoeff(i1,i2,i3,8)*ct2 *d005i 
                                      lap2hSq200(i1,i2,i3,0) = lap2hSq(i1+1,i2,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1-1,i2,i3,0)
                                      lap2hSq020(i1,i2,i3,0) = lap2hSq(i1,i2+1,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2-1,i3,0)
                                      lap2hSq002(i1,i2,i3,0) = lap2hSq(i1,i2,i3+1,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2,i3-1,0)
                                      lap2hSq100i = lap2hSq(i1+1,i2,i3,0) - lap2hSq(i1-1,i2,i3,0)
                                      lap2hSq010i = lap2hSq(i1,i2+1,i3,0) - lap2hSq(i1,i2-1,i3,0)
                                      lap2hSq001i = lap2hSq(i1,i2,i3+1,0) - lap2hSq(i1,i2,i3-1,0)
                                      lap2hSq110i = lap2hSq(i1+1,i2+1,i3,0) - lap2hSq(i1-1,i2+1,i3,0) - lap2hSq(i1+1,i2-1,i3,0) + lap2hSq(i1-1,i2-1,i3,0)
                                      lap2hSq101i = lap2hSq(i1+1,i2,i3+1,0) - lap2hSq(i1-1,i2,i3+1,0) - lap2hSq(i1+1,i2,i3-1,0) + lap2hSq(i1-1,i2,i3-1,0)
                                      lap2hSq011i = lap2hSq(i1,i2+1,i3+1,0) - lap2hSq(i1,i2-1,i3+1,0) - lap2hSq(i1,i2+1,i3-1,0) + lap2hSq(i1,i2-1,i3-1,0)
                                      lap2hCubed(i1,i2,i3,0) =  + lapCoeff(i1,i2,i3,0)*lap2hSq200(i1,i2,i3,0) + lapCoeff(i1,i2,i3,1)*lap2hSq020(i1,i2,i3,0) + lapCoeff(i1,i2,i3,2)*lap2hSq002(i1,i2,i3,0) + lapCoeff(i1,i2,i3,3)*lap2hSq110i  + lapCoeff(i1,i2,i3,4)*lap2hSq101i  + lapCoeff(i1,i2,i3,5)*lap2hSq011i  + lapCoeff(i1,i2,i3,6)*lap2hSq100i  + lapCoeff(i1,i2,i3,7)*lap2hSq010i  + lapCoeff(i1,i2,i3,8)*lap2hSq001i   
                                      end if ! mask .ne. 0
                                    end do
                                    end do
                                    end do
                ! --- SPLIT LOOPS
                                do i3=n3a,n3b
                                do i2=n2a,n2b
                                do i1=n1a,n1b
                                  if( mask(i1,i2,i3).ne.0 )then
                   ! --- Laplacian squared to order 4 = 
                   !  lap2h*( lap4h ) + corrections*( Lap2h )
                                      lap4h200(i1,i2,i3,0) = lap4h(i1+1,i2,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1-1,i2,i3,0)
                                      lap4h020(i1,i2,i3,0) = lap4h(i1,i2+1,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2-1,i3,0)
                                      lap4h002(i1,i2,i3,0) = lap4h(i1,i2,i3+1,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2,i3-1,0)
                                      lap2h400(i1,i2,i3,0) = lap2h200(i1+1,i2,i3,0) - 2*lap2h200(i1,i2,i3,0) + lap2h200(i1-1,i2,i3,0)
                                      lap2h040(i1,i2,i3,0) = lap2h020(i1,i2+1,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1,i2-1,i3,0)
                                      lap2h004(i1,i2,i3,0) = lap2h002(i1,i2,i3+1,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1,i2,i3-1,0)
                                      lap2h220(i1,i2,i3,0) = lap2h020(i1+1,i2,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1-1,i2,i3,0)
                                      lap2h202(i1,i2,i3,0) = lap2h002(i1+1,i2,i3,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1-1,i2,i3,0)
                                      lap2h022(i1,i2,i3,0) = lap2h002(i1,i2+1,i3,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1,i2-1,i3,0)
                                      lap4h100i = lap4h(i1+1,i2,i3,0) - lap4h(i1-1,i2,i3,0)
                                      lap4h010i = lap4h(i1,i2+1,i3,0) - lap4h(i1,i2-1,i3,0)
                                      lap4h001i = lap4h(i1,i2,i3+1,0) - lap4h(i1,i2,i3-1,0)
                                      lap4h110i = lap4h(i1+1,i2+1,i3,0) - lap4h(i1-1,i2+1,i3,0) - lap4h(i1+1,i2-1,i3,0) + lap4h(i1-1,i2-1,i3,0)
                                      lap4h101i = lap4h(i1+1,i2,i3+1,0) - lap4h(i1-1,i2,i3+1,0) - lap4h(i1+1,i2,i3-1,0) + lap4h(i1-1,i2,i3-1,0)
                                      lap4h011i = lap4h(i1,i2+1,i3+1,0) - lap4h(i1,i2-1,i3+1,0) - lap4h(i1,i2+1,i3-1,0) + lap4h(i1,i2-1,i3-1,0)
                                      lap2h300i = lap2h200(i1+1,i2,i3,0) - lap2h200(i1-1,i2,i3,0)
                                      lap2h030i = lap2h020(i1,i2+1,i3,0) - lap2h020(i1,i2-1,i3,0)
                                      lap2h003i = lap2h002(i1,i2,i3+1,0) - lap2h002(i1,i2,i3-1,0)
                                      lap2h310i = lap2h200(i1+1,i2+1,i3,0) - lap2h200(i1-1,i2+1,i3,0) - lap2h200(i1+1,i2-1,i3,0) + lap2h200(i1-1,i2-1,i3,0)
                                      lap2h130i = lap2h020(i1+1,i2+1,i3,0) - lap2h020(i1-1,i2+1,i3,0) - lap2h020(i1+1,i2-1,i3,0) + lap2h020(i1-1,i2-1,i3,0)
                                      lap2h301i = lap2h200(i1+1,i2,i3+1,0) - lap2h200(i1-1,i2,i3+1,0) - lap2h200(i1+1,i2,i3-1,0) + lap2h200(i1-1,i2,i3-1,0)
                                      lap2h103i = lap2h002(i1+1,i2,i3+1,0) - lap2h002(i1-1,i2,i3+1,0) - lap2h002(i1+1,i2,i3-1,0) + lap2h002(i1-1,i2,i3-1,0)
                                      lap2h031i = lap2h020(i1,i2+1,i3+1,0) - lap2h020(i1,i2-1,i3+1,0) - lap2h020(i1,i2+1,i3-1,0) + lap2h020(i1,i2-1,i3-1,0)
                                      lap2h013i = lap2h002(i1,i2+1,i3+1,0) - lap2h002(i1,i2-1,i3+1,0) - lap2h002(i1,i2+1,i3-1,0) + lap2h002(i1,i2-1,i3-1,0)
                                      lap4hSq(i1,i2,i3,0) =     lapCoeff(i1,i2,i3,0)*( lap4h200(i1,i2,i3,0) + crr1*lap2h400(i1,i2,i3,0) )    + lapCoeff(i1,i2,i3,1)*( lap4h020(i1,i2,i3,0) + css1*lap2h040(i1,i2,i3,0) )     + lapCoeff(i1,i2,i3,2)*( lap4h002(i1,i2,i3,0) + ctt1*lap2h004(i1,i2,i3,0) )     + lapCoeff(i1,i2,i3,3)*( lap4h110i + cr1*lap2h310i + cs1*lap2h130i ) + lapCoeff(i1,i2,i3,4)*( lap4h101i + cr1*lap2h301i + ct1*lap2h103i ) + lapCoeff(i1,i2,i3,5)*( lap4h011i + cs1*lap2h031i + ct1*lap2h013i ) + lapCoeff(i1,i2,i3,6)*( lap4h100i + cr1 *lap2h300i )    + lapCoeff(i1,i2,i3,7)*( lap4h010i + cs1 *lap2h030i )    + lapCoeff(i1,i2,i3,8)*( lap4h001i + ct1 *lap2h003i )      
                                      end if ! mask .ne. 0
                                    end do
                                    end do
                                    end do
               ! ===========  FINAL LOOP TO FILL IN THE SOLUTION ============
                              numGhost1=0;
                              n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
                              n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
                              n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
                                do i3=n3a,n3b
                                do i2=n2a,n2b
                                do i1=n1a,n1b
                                  if( mask(i1,i2,i3).ne.0 )then
                                      d800i = d600(i1+1,i2,i3,0) - 2*d600(i1,i2,i3,0) + d600(i1-1,i2,i3,0)
                                      d080i = d060(i1,i2+1,i3,0) - 2*d060(i1,i2,i3,0) + d060(i1,i2-1,i3,0)
                                      d008i = d006(i1,i2,i3+1,0) - 2*d006(i1,i2,i3,0) + d006(i1,i2,i3-1,0)
                                      d700i = d600(i1+1,i2,i3,0) - d600(i1-1,i2,i3,0)
                                      d070i = d060(i1,i2+1,i3,0) - d060(i1,i2-1,i3,0)
                                      d007i = d006(i1,i2,i3+1,0) - d006(i1,i2,i3-1,0)
                                      d710i = d600(i1+1,i2+1,i3,0) - d600(i1-1,i2+1,i3,0) - d600(i1+1,i2-1,i3,0) + d600(i1-1,i2-1,i3,0)
                                      d170i = d060(i1+1,i2+1,i3,0) - d060(i1-1,i2+1,i3,0) - d060(i1+1,i2-1,i3,0) + d060(i1-1,i2-1,i3,0)
                                      d701i = d600(i1+1,i2,i3+1,0) - d600(i1-1,i2,i3+1,0) - d600(i1+1,i2,i3-1,0) + d600(i1-1,i2,i3-1,0)
                                      d107i = d006(i1+1,i2,i3+1,0) - d006(i1-1,i2,i3+1,0) - d006(i1+1,i2,i3-1,0) + d006(i1-1,i2,i3-1,0)
                                      d071i = d060(i1,i2+1,i3+1,0) - d060(i1,i2-1,i3+1,0) - d060(i1,i2+1,i3-1,0) + d060(i1,i2-1,i3-1,0)
                                      d017i = d006(i1,i2+1,i3+1,0) - d006(i1,i2-1,i3+1,0) - d006(i1,i2+1,i3-1,0) + d006(i1,i2-1,i3-1,0)
                                      d530i = d420(i1+1,i2+1,i3,0) - d420(i1-1,i2+1,i3,0) - d420(i1+1,i2-1,i3,0) + d420(i1-1,i2-1,i3,0)
                                      d350i = d240(i1+1,i2+1,i3,0) - d240(i1-1,i2+1,i3,0) - d240(i1+1,i2-1,i3,0) + d240(i1-1,i2-1,i3,0)
                                      d503i = d402(i1+1,i2,i3+1,0) - d402(i1-1,i2,i3+1,0) - d402(i1+1,i2,i3-1,0) + d402(i1-1,i2,i3-1,0)
                                      d305i = d204(i1+1,i2,i3+1,0) - d204(i1-1,i2,i3+1,0) - d204(i1+1,i2,i3-1,0) + d204(i1-1,i2,i3-1,0)
                                      d053i = d042(i1,i2+1,i3+1,0) - d042(i1,i2-1,i3+1,0) - d042(i1,i2+1,i3-1,0) + d042(i1,i2-1,i3-1,0)
                                      d035i = d024(i1,i2+1,i3+1,0) - d024(i1,i2-1,i3+1,0) - d024(i1,i2+1,i3-1,0) + d024(i1,i2-1,i3-1,0)
                   ! --- Laplacian to order 8 = lap6h + corrections 
                                      lap8h = lap6h(i1,i2,i3,0)                                                         + lapCoeff(i1,i2,i3,0)*crr3*d800i                                               + lapCoeff(i1,i2,i3,1)*css3*d080i                                               + lapCoeff(i1,i2,i3,2)*ctt3*d008i                                               + lapCoeff(i1,i2,i3,3)*(cr3*d710i + cs3*d170i + cr2*cs1*d530i + cr1*cs2*d350i ) + lapCoeff(i1,i2,i3,4)*(cr3*d701i + ct3*d107i + cr2*ct1*d503i + cr1*ct2*d305i ) + lapCoeff(i1,i2,i3,5)*(cs3*d071i + ct3*d017i + cs2*ct1*d053i + cs1*ct2*d035i ) + lapCoeff(i1,i2,i3,6)* cr3*d700i                                               + lapCoeff(i1,i2,i3,7)* cs3*d070i                                               + lapCoeff(i1,i2,i3,8)* ct3*d007i 
                   ! --- Laplacian^4 4p (4th power) order 2: 
                                      lap2hCubed200i = lap2hCubed(i1+1,i2,i3,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1-1,i2,i3,0)
                                      lap2hCubed020i = lap2hCubed(i1,i2+1,i3,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1,i2-1,i3,0)
                                      lap2hCubed002i = lap2hCubed(i1,i2,i3+1,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1,i2,i3-1,0)
                                      lap2hCubed100i = lap2hCubed(i1+1,i2,i3,0) - lap2hCubed(i1-1,i2,i3,0)
                                      lap2hCubed010i = lap2hCubed(i1,i2+1,i3,0) - lap2hCubed(i1,i2-1,i3,0)
                                      lap2hCubed001i = lap2hCubed(i1,i2,i3+1,0) - lap2hCubed(i1,i2,i3-1,0)
                                      lap2hCubed110i = lap2hCubed(i1+1,i2+1,i3,0) - lap2hCubed(i1-1,i2+1,i3,0) - lap2hCubed(i1+1,i2-1,i3,0) + lap2hCubed(i1-1,i2-1,i3,0)
                                      lap2hCubed101i = lap2hCubed(i1+1,i2,i3+1,0) - lap2hCubed(i1-1,i2,i3+1,0) - lap2hCubed(i1+1,i2,i3-1,0) + lap2hCubed(i1-1,i2,i3-1,0)
                                      lap2hCubed011i = lap2hCubed(i1,i2+1,i3+1,0) - lap2hCubed(i1,i2-1,i3+1,0) - lap2hCubed(i1,i2+1,i3-1,0) + lap2hCubed(i1,i2-1,i3-1,0)
                                      lap2h4p  =                             + lapCoeff(i1,i2,i3,0)*lap2hCubed200i  + lapCoeff(i1,i2,i3,1)*lap2hCubed020i  + lapCoeff(i1,i2,i3,2)*lap2hCubed002i  + lapCoeff(i1,i2,i3,3)*lap2hCubed110i  + lapCoeff(i1,i2,i3,4)*lap2hCubed101i  + lapCoeff(i1,i2,i3,5)*lap2hCubed011i  + lapCoeff(i1,i2,i3,6)*lap2hCubed100i  + lapCoeff(i1,i2,i3,7)*lap2hCubed010i  + lapCoeff(i1,i2,i3,8)*lap2hCubed001i    
                   ! --- Laplacian squared to order 6 :
                   !   Lap6h = Lap4h + M4  = (Lap2h) + M2 + M4 
                   !   Lap6h*Lap6h = [ (Lap2h) + M2 + M4 ] [ (Lap2h) + M2 + M4 ]
                   !               = Lap2h*Lap6h + M2*Lap4h + M4*Lap2h + O(h^6)
                                      lap6h200i = lap6h(i1+1,i2,i3,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1-1,i2,i3,0)
                                      lap6h020i = lap6h(i1,i2+1,i3,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1,i2-1,i3,0)
                                      lap6h002i = lap6h(i1,i2,i3+1,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1,i2,i3-1,0)
                                      lap6h100i = lap6h(i1+1,i2,i3,0) - lap6h(i1-1,i2,i3,0)
                                      lap6h010i = lap6h(i1,i2+1,i3,0) - lap6h(i1,i2-1,i3,0)
                                      lap6h001i = lap6h(i1,i2,i3+1,0) - lap6h(i1,i2,i3-1,0)
                                      lap6h110i = lap6h(i1+1,i2+1,i3,0) - lap6h(i1-1,i2+1,i3,0) - lap6h(i1+1,i2-1,i3,0) + lap6h(i1-1,i2-1,i3,0)
                                      lap6h101i = lap6h(i1+1,i2,i3+1,0) - lap6h(i1-1,i2,i3+1,0) - lap6h(i1+1,i2,i3-1,0) + lap6h(i1-1,i2,i3-1,0)
                                      lap6h011i = lap6h(i1,i2+1,i3+1,0) - lap6h(i1,i2-1,i3+1,0) - lap6h(i1,i2+1,i3-1,0) + lap6h(i1,i2-1,i3-1,0)
                                      lap4h400i = lap4h200(i1+1,i2,i3,0) - 2*lap4h200(i1,i2,i3,0) + lap4h200(i1-1,i2,i3,0)
                                      lap4h040i = lap4h020(i1,i2+1,i3,0) - 2*lap4h020(i1,i2,i3,0) + lap4h020(i1,i2-1,i3,0)
                                      lap4h004i = lap4h002(i1,i2,i3+1,0) - 2*lap4h002(i1,i2,i3,0) + lap4h002(i1,i2,i3-1,0)
                                      lap4h300i = lap4h200(i1+1,i2,i3,0) - lap4h200(i1-1,i2,i3,0)
                                      lap4h030i = lap4h020(i1,i2+1,i3,0) - lap4h020(i1,i2-1,i3,0)
                                      lap4h003i = lap4h002(i1,i2,i3+1,0) - lap4h002(i1,i2,i3-1,0)
                                      lap4h310i = lap4h200(i1+1,i2+1,i3,0) - lap4h200(i1-1,i2+1,i3,0) - lap4h200(i1+1,i2-1,i3,0) + lap4h200(i1-1,i2-1,i3,0)
                                      lap4h130i = lap4h020(i1+1,i2+1,i3,0) - lap4h020(i1-1,i2+1,i3,0) - lap4h020(i1+1,i2-1,i3,0) + lap4h020(i1-1,i2-1,i3,0)
                                      lap4h301i = lap4h200(i1+1,i2,i3+1,0) - lap4h200(i1-1,i2,i3+1,0) - lap4h200(i1+1,i2,i3-1,0) + lap4h200(i1-1,i2,i3-1,0)
                                      lap4h103i = lap4h002(i1+1,i2,i3+1,0) - lap4h002(i1-1,i2,i3+1,0) - lap4h002(i1+1,i2,i3-1,0) + lap4h002(i1-1,i2,i3-1,0)
                                      lap4h031i = lap4h020(i1,i2+1,i3+1,0) - lap4h020(i1,i2-1,i3+1,0) - lap4h020(i1,i2+1,i3-1,0) + lap4h020(i1,i2-1,i3-1,0)
                                      lap4h013i = lap4h002(i1,i2+1,i3+1,0) - lap4h002(i1,i2-1,i3+1,0) - lap4h002(i1,i2+1,i3-1,0) + lap4h002(i1,i2-1,i3-1,0)
                                      lap2h600i = lap2h400(i1+1,i2,i3,0) - 2*lap2h400(i1,i2,i3,0) + lap2h400(i1-1,i2,i3,0)
                                      lap2h060i = lap2h040(i1,i2+1,i3,0) - 2*lap2h040(i1,i2,i3,0) + lap2h040(i1,i2-1,i3,0)
                                      lap2h006i = lap2h004(i1,i2,i3+1,0) - 2*lap2h004(i1,i2,i3,0) + lap2h004(i1,i2,i3-1,0)
                                      lap2h500i = lap2h400(i1+1,i2,i3,0) - lap2h400(i1-1,i2,i3,0)
                                      lap2h050i = lap2h040(i1,i2+1,i3,0) - lap2h040(i1,i2-1,i3,0)
                                      lap2h005i = lap2h004(i1,i2,i3+1,0) - lap2h004(i1,i2,i3-1,0)
                                      lap2h510i = lap2h400(i1+1,i2+1,i3,0) - lap2h400(i1-1,i2+1,i3,0) - lap2h400(i1+1,i2-1,i3,0) + lap2h400(i1-1,i2-1,i3,0)
                                      lap2h150i = lap2h040(i1+1,i2+1,i3,0) - lap2h040(i1-1,i2+1,i3,0) - lap2h040(i1+1,i2-1,i3,0) + lap2h040(i1-1,i2-1,i3,0)
                                      lap2h330i = lap2h220(i1+1,i2+1,i3,0) - lap2h220(i1-1,i2+1,i3,0) - lap2h220(i1+1,i2-1,i3,0) + lap2h220(i1-1,i2-1,i3,0)
                                      lap2h501i = lap2h400(i1+1,i2,i3+1,0) - lap2h400(i1-1,i2,i3+1,0) - lap2h400(i1+1,i2,i3-1,0) + lap2h400(i1-1,i2,i3-1,0)
                                      lap2h105i = lap2h004(i1+1,i2,i3+1,0) - lap2h004(i1-1,i2,i3+1,0) - lap2h004(i1+1,i2,i3-1,0) + lap2h004(i1-1,i2,i3-1,0)
                                      lap2h051i = lap2h040(i1,i2+1,i3+1,0) - lap2h040(i1,i2-1,i3+1,0) - lap2h040(i1,i2+1,i3-1,0) + lap2h040(i1,i2-1,i3-1,0)
                                      lap2h015i = lap2h004(i1,i2+1,i3+1,0) - lap2h004(i1,i2-1,i3+1,0) - lap2h004(i1,i2+1,i3-1,0) + lap2h004(i1,i2-1,i3-1,0)
                                      lap2h303i = lap2h202(i1+1,i2,i3+1,0) - lap2h202(i1-1,i2,i3+1,0) - lap2h202(i1+1,i2,i3-1,0) + lap2h202(i1-1,i2,i3-1,0)
                                      lap2h033i = lap2h022(i1,i2+1,i3+1,0) - lap2h022(i1,i2-1,i3+1,0) - lap2h022(i1,i2+1,i3-1,0) + lap2h022(i1,i2-1,i3-1,0)
                                      lap6hSq =                                                                                     lapCoeff(i1,i2,i3,0)*(lap6h200i + crr1*lap4h400i + crr2*lap2h600i )                       + lapCoeff(i1,i2,i3,1)*(lap6h020i + css1*lap4h040i + css2*lap2h060i )                       + lapCoeff(i1,i2,i3,2)*(lap6h002i + ctt1*lap4h004i + ctt2*lap2h006i )                       + lapCoeff(i1,i2,i3,3)*(lap6h110i +  cr1*lap4h310i +  cr2*lap2h510i                         +  cs1*lap4h130i +  cs2*lap2h150i + cr1*cs1*lap2h330i )   + lapCoeff(i1,i2,i3,4)*(lap6h101i +  cr1*lap4h301i +  cr2*lap2h501i                         +  ct1*lap4h103i +  ct2*lap2h105i + cr1*ct1*lap2h303i )   + lapCoeff(i1,i2,i3,5)*(lap6h011i +  cs1*lap4h031i +  cs2*lap2h051i                         +  ct1*lap4h013i +  ct2*lap2h015i + cs1*ct1*lap2h033i )   + lapCoeff(i1,i2,i3,6)*(lap6h100i +  cr1*lap4h300i +  cr2*lap2h500i )                       + lapCoeff(i1,i2,i3,7)*(lap6h010i +  cs1*lap4h030i +  cs2*lap2h050i )                       + lapCoeff(i1,i2,i3,8)*(lap6h010i +  ct1*lap4h003i +  ct2*lap2h005i )                         
                   ! --- Laplacian CUBED to order 4 
                                      lap4hSq200i = lap4hSq(i1+1,i2,i3,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1-1,i2,i3,0)
                                      lap4hSq020i = lap4hSq(i1,i2+1,i3,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1,i2-1,i3,0)
                                      lap4hSq002i = lap4hSq(i1,i2,i3+1,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1,i2,i3-1,0)
                                      lap4hSq100i = lap4hSq(i1+1,i2,i3,0) - lap4hSq(i1-1,i2,i3,0)
                                      lap4hSq010i = lap4hSq(i1,i2+1,i3,0) - lap4hSq(i1,i2-1,i3,0)
                                      lap4hSq001i = lap4hSq(i1,i2,i3+1,0) - lap4hSq(i1,i2,i3-1,0)
                                      lap4hSq110i = lap4hSq(i1+1,i2+1,i3,0) - lap4hSq(i1-1,i2+1,i3,0) - lap4hSq(i1+1,i2-1,i3,0) + lap4hSq(i1-1,i2-1,i3,0)
                                      lap4hSq101i = lap4hSq(i1+1,i2,i3+1,0) - lap4hSq(i1-1,i2,i3+1,0) - lap4hSq(i1+1,i2,i3-1,0) + lap4hSq(i1-1,i2,i3-1,0)
                                      lap4hSq011i = lap4hSq(i1,i2+1,i3+1,0) - lap4hSq(i1,i2-1,i3+1,0) - lap4hSq(i1,i2+1,i3-1,0) + lap4hSq(i1,i2-1,i3-1,0)
                                      lap2hSq400i = lap2hSq200(i1+1,i2,i3,0) - 2*lap2hSq200(i1,i2,i3,0) + lap2hSq200(i1-1,i2,i3,0)
                                      lap2hSq040i = lap2hSq020(i1,i2+1,i3,0) - 2*lap2hSq020(i1,i2,i3,0) + lap2hSq020(i1,i2-1,i3,0)
                                      lap2hSq004i = lap2hSq002(i1,i2,i3+1,0) - 2*lap2hSq002(i1,i2,i3,0) + lap2hSq002(i1,i2,i3-1,0)
                                      lap2hSq300i = lap2hSq200(i1+1,i2,i3,0) - lap2hSq200(i1-1,i2,i3,0)
                                      lap2hSq030i = lap2hSq020(i1,i2+1,i3,0) - lap2hSq020(i1,i2-1,i3,0)
                                      lap2hSq003i = lap2hSq002(i1,i2,i3+1,0) - lap2hSq002(i1,i2,i3-1,0)
                                      lap2hSq310i = lap2hSq200(i1+1,i2+1,i3,0) - lap2hSq200(i1-1,i2+1,i3,0) - lap2hSq200(i1+1,i2-1,i3,0) + lap2hSq200(i1-1,i2-1,i3,0)
                                      lap2hSq130i = lap2hSq020(i1+1,i2+1,i3,0) - lap2hSq020(i1-1,i2+1,i3,0) - lap2hSq020(i1+1,i2-1,i3,0) + lap2hSq020(i1-1,i2-1,i3,0)
                                      lap2hSq301i = lap2hSq200(i1+1,i2,i3+1,0) - lap2hSq200(i1-1,i2,i3+1,0) - lap2hSq200(i1+1,i2,i3-1,0) + lap2hSq200(i1-1,i2,i3-1,0)
                                      lap2hSq103i = lap2hSq002(i1+1,i2,i3+1,0) - lap2hSq002(i1-1,i2,i3+1,0) - lap2hSq002(i1+1,i2,i3-1,0) + lap2hSq002(i1-1,i2,i3-1,0)
                                      lap2hSq031i = lap2hSq020(i1,i2+1,i3+1,0) - lap2hSq020(i1,i2-1,i3+1,0) - lap2hSq020(i1,i2+1,i3-1,0) + lap2hSq020(i1,i2-1,i3-1,0)
                                      lap2hSq013i = lap2hSq002(i1,i2+1,i3+1,0) - lap2hSq002(i1,i2-1,i3+1,0) - lap2hSq002(i1,i2+1,i3-1,0) + lap2hSq002(i1,i2-1,i3-1,0)
                                      lap4hCubed =                                                    lapCoeff(i1,i2,i3,0)*(lap4hSq200i + crr1*lap2hSq400i )   + lapCoeff(i1,i2,i3,1)*(lap4hSq020i + css1*lap2hSq040i )   + lapCoeff(i1,i2,i3,2)*(lap4hSq002i + ctt1*lap2hSq004i )   + lapCoeff(i1,i2,i3,3)*(lap4hSq110i +  cr1*lap2hSq310i     +  cs1*lap2hSq130i )   + lapCoeff(i1,i2,i3,4)*(lap4hSq101i +  cr1*lap2hSq301i     +  ct1*lap2hSq103i )   + lapCoeff(i1,i2,i3,5)*(lap4hSq011i +  cs1*lap2hSq031i     +  ct1*lap2hSq013i )   + lapCoeff(i1,i2,i3,6)*(lap4hSq100i + cr1 *lap2hSq300i )   + lapCoeff(i1,i2,i3,7)*(lap4hSq010i + cs1 *lap2hSq030i )   + lapCoeff(i1,i2,i3,8)*(lap4hSq001i + ct1 *lap2hSq003i )     
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
                                 !fv(m) = fv(m) -( f(i1,i2,i3,freq) + cdtSqBy12*( cSq*(fxx23(i1,i2,i3,freq) + fyy23(i1,i2,i3,freq) + fzz23(i1,i2,i3,freq)) - omega*omega*f(i1,i2,i3,freq)) )*coswt 
                                                  end do ! do freq  
                                            else if( addForcing.ne.0 )then  
                                                  fv(m) = f(i1,i2,i3,0)
                                            end if
                   ! --- Modified equation space-time update ----
                                      un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m)  + cdtsq*( lap8h )               + cdtPow4By12*( lap6hSq )       + cdtPow6By360*( lap4hCubed )   + cdtPow8By20160*( lap2h4p )    +dtSq*fv(m)                    
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
