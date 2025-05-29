! This file automatically generated from abcWave.bf90 with bpp.
! *******************************************************************************
!   Absorbing boundary conditions
!
!  THIS FILE STARTED FROM THE VERSION in cg/mx/src
!
! *******************************************************************************

! These next include files will define the macros that will define the difference approximations
! The actual macro is called below
! Use this next macro to declare the statement functions that are defined below
! To include derivatives of rx use OPTION=RX


! Define statement functions for difference approximations of order 2 
! To include derivatives of rx use OPTION=RX
! To include derivatives of rx use OPTION=RX



! Use this next macro to declare the statement functions that are defined below
! To include derivatives of rx use OPTION=RX


! Define statement functions for difference approximations of order 4 
! To include derivatives of rx use OPTION=RX
! To include derivatives of rx use OPTION=RX

! Here are macros that define the planeWave solution
! -*- mode: f90; -*-

! **************************************************
! Here are macros that define the:
!      planeWave solution 
! **************************************************

! ======================================================================
!  Slow start function 
!    tba = length of slow start interval (<0 mean no slow start)
! ======================================================================

! cubic ramp
! tba=max(REAL_EPSILON,tb-ta);
! dta=t-ta;
	  
! This (cubic) ramp has 1-derivative zero at t=0 and t=tba

! This ramp has 3-derivatives zero at t=0 and t=1
! This is from ramp.maple
! r=-84*t**5+35*t**4-20*t**7+70*t**6
! rt=-420*t**4+140*t**3-140*t**6+420*t**5
! rtt=-1680*t**3+420*t**2-840*t**5+2100*t**4
! rttt=-5040*t**2+840*t-4200*t**4+8400*t**3


! This ramp has 4-derivatives zero at t=0 and t=1
! This is from ramp.maple
! r=126*(t)**5-315*(t)**8+70*(t)**9-420*(t)**6+540*(t)**7
! rt=630*(t)**4-2520*(t)**7+630*(t)**8-2520*(t)**5+3780*(t)**6
! rtt=2520*(t)**3-17640*(t)**6+5040*(t)**7-12600*(t)**4+22680*(t)**5
! rttt=7560*(t)**2-105840*(t)**5+35280*(t)**6-50400*(t)**3+113400*(t)**4


! ============================================================
!  Initialize parameters for the boundary forcing
!   tba: slow start time interval -- no slow start if this is negative
! ===========================================================

! **************** Here is the new generic plane wave solution *******************

! component n=ex,ey,ez, hx,hy,hz (assumes ex=0)
! one time derivative:
! two time derivatives:
! three time derivatives:

! *************** Here is the 2D planeWave solution ******************************


! one time derivative:

! two time derivatives:

! three time derivatives:

! four time derivatives:

! Here are the slow start versions

! one time derivative:

! two time derivatives:

! three time derivatives:

! four time derivatives:


! **************** Here is the 3D planeWave solution ***************************************



! one time derivative:


! two time derivatives:


! three time derivatives:


! four time derivatives:


! Here are the slow start versions


! one time derivative:


! two time derivatives:

! three time derivatives:

! four time derivatives:


! -------------------------------------------------------------------
! Helper function: Return minus the second time derivative
! -------------------------------------------------------------------


! --------------------------------------------------------------------
! Evaluate the plane wave in 2D
! 
!  x,y,t (input) : point to evaluate at 
!  numberOfTimeDerivatives : evaluate this time derivative
!  ubc(.)  (output) : ubc(ex), etc. 
! --------------------------------------------------------------------


! --------------------------------------------------------------------
! Evaluate the plane wave in 3D
! 
!  x,y,z,t (input) : point to evaluate at 
!  numberOfTimeDerivatives : evaluate this time derivative
!  ubc(.)  (output) : ubc(ex), etc. 
! --------------------------------------------------------------------

! Evaluate the twilight-zone forcing 





! ************************************************************************************************
!  This macro is used for looping over the faces of a grid to assign booundary conditions
!
! extra: extra points to assign
!          Case 1: extra=numberOfGhostPoints -- for assigning extended boundaries
!          Case 2: extra=-1 -- for assigning ghost points but not including extended boundaries
! numberOfGhostPoints : number of ghost points (1 for 2nd order, 2 for fourth-order ...)
! ***********************************************************************************************


! ========================================================================
! Begin loop over edges in 3D
! ========================================================================




! ABC - Engquist Majda order 2
! This is only a first order in time approx.
! Generalized form:
! u.xt = c1abcem2*u.xx + c2abcem2*( u.yy + u.zz )
!   Taylor: p0=1 p2=-1/2
!   Cheby:  p0=1.00023, p2=-.515555

! -------------------- CARTESIAN GRID ---------------------

! Here are first-order-in-time formula that do not require other ghost point at new time (un)
!  Solve for ghost value from:
!      D+t D0x ( u^n ) = c1abcem2 * D+xD-x u^n + c2abcem2 D+yD-y u^n  + f(t^n+dt/2)
! These are used at corners.


! Here are 2nd-order in time approximations -- centered in space-time, solve for ghost at new time: 
!   D+t D0x ( u^n ) = A+t[ c1abcem2 * D+xD-x u^n + c2abcem2 D+yD-y u^n ] + f(t^n+dt/2)
!   Average in time operator:  A+t u^n = .5*( u^(n+1) + u^n )





! --------------------------------------------------------------
! Macro: 
!     ------ 2nd-order accurate corner approximations ----
!
! Parameters: 
!   side1,side2 : 0,1 to denote which corner in 2D
! --------------------------------------------------------------














! ======================================================================================
! Setup Macro to apply the 2nd-order accurate Engquist-Majda ABC on a curvilinear grid
! 
! On a Curvilinear grid we write:
!        u_tt = L u  
!        u_tt = D_n^2 u + (L-D_n^2) u 
! where the "normal" derivative is 
!        D_n = sqrt( rx^2 + ry^2) D_r 
! 
!    sqrt( rx^2 + ry^2) u_{rt} = c1abcem2*( (rx^2 + ry^2) u_{rr} ) + c2abcem2*( L - (rx^2 + ry^2) u_{rr} )
! 
! ======================================================================================


! ======================================================================================
! Macro to apply the 2nd-order accurate Engquist-Majda ABC on a curvilinear grid
! ======================================================================================

! First-order in time explicit version: 




! ============================================================================
! Macro to extrapolate points adjacent to an edge or corner
! EXCEPT for the first ghost point on the extended boundary
! ============================================================================

! --------------------------------------------------------------------------
! Macro: Evaluate the forcing for the ABC EM2 for twilight zone   
!        Cartesian grid version
!    DIR = X or Y or Z
!    DIM = 2 or 3
!    tf : evaluate the forcing at this time
! Output: 
!    force(0:2) 
! --------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! Macro: Get the forcing for the trapezoid rule by averaging times tp amd tf
!   This will make the scheme exact for degree 2 polynomials in time 
! --------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! Macro: Extrapolate the ghost point along a given direction
! --------------------------------------------------------------------------


! ----------------------------------------------------------------------------
! Macro: Assign edges and corners : ORDER OF ACCURACY 2
! --------------------------------------------------------------------------



! ----------------------------------------------------------------------------
! Macro: Evaluate the residuals in the EM conditions at a corner
! --------------------------------------------------------------------------

! ----------------------------------------------------------------------------
! Macro: Assign edges and corners : ORDER OF ACCURACY 4
! --------------------------------------------------------------------------


      subroutine abcWave( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,ndf1a,ndf1b,ndf2a,ndf2b,ndf3a,ndf3b,gridIndexRange, u, un, f,mask,rsxy, xy,bc, boundaryCondition, ipar, rpar, ierr )
! ===================================================================================
!  Absorbing boundary conditions for the Wave Equation
!
!  gridType : 0=rectangular, 1=curvilinear
!  useForcing : 1=use f for RHS to BC
!  side,axis : 0:1 and 0:2
!
!  u : solution at time t-dt
!  un : solution at time t (apply BC to this solution)
!
! ===================================================================================

  implicit none

  integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, ndf1a,ndf1b,ndf2a,ndf2b,ndf3a,ndf3b,n1a,n1b,n2a,n2b,n3a,n3b, ndc, bc,ierr

  real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
  real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
  real f(ndf1a:ndf1b,ndf2a:ndf2b,ndf3a:ndf3b,0:*)
  integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
  real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
  real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
  integer gridIndexRange(0:1,0:2)

  integer ipar(0:*),boundaryCondition(0:1,0:2)
  real rpar(0:*)

  !     --- local variables ----
      
  integer side,axis,gridType,orderOfAccuracy,orderOfExtrapolation,useForcing,ex,ey,ez,hx,hy,hz,useWhereMask,grid,debug,side1,side2,side3,forcingOption,myid
  real dx(0:2),dr(0:2),t,ep,dt,c,cEM2     
  real dxa,dya,dza
  integer axisp1,axisp2,i1,i2,i3,is1,is2,is3,js1,js2,js3,ks1,ks2,ks3,is,j1,j2,j3,m1,m2,m3,mSum
  integer ip1,ip2,ip3,ig1,ig2,ig3,ghost1,ghost2,ghost3,mm
  integer extra,extra1a,extra1b,extra2a,extra2b,extra3a,extra3b,numberOfGhostPoints
  integer edgeDirection,sidea,sideb,sidec,bc1,bc2,bc3

  real p0,p2,q0,q2,c1abcem2,c2abcem2
  real an1,an2,an3,aNorm,epsX

  real rx0,ry0,rz0 , rxx0,ryy0, rzz0 
  real dr0,cxt,cxx,cyy,czz,cm1,g,bxx,byy,bzz
  real rxNorm, rxNormSq, Dn2, Lu, ur0,urr0, unr0, unrr0
  real ux0,uy0,uz0, uxx0,uyy0,uzz0
  real unx0,uny0,unz0, unxx0,unyy0,unzz0
  real t0,t1,t2
  real x,y,z,eyTrue
  real forcex(0:11),forcey(0:11),forcez(0:11),forcep(0:11),forcef(0:11)
  real tp,tm,tf,utx,uty,utz,uxx,uyy,uzz

  integer ksv(0:2)
  integer isign1,isign2,idir,extrapOrder
  real r1,r2,f1,f2,a11,a12,a21,a22,uA,uB,det

  real ux,vy,vxy,alpha,uGhost,aGhost

  real eps,mu,kx,ky,kz,slowStartInterval,twoPi,cc

  real ax,ay,az,aSq,div,divCoeff,res
  logical adjacentFaceIsABC, applyABC
  integer projectDivLine

  ! For solves with 4 unknowns
  real a4(4,4), b4(4), work4(4), rcond, resv(4), uv(10)
  integer job, ipvt4(4), numberOfEquations,axis1,axis2,axis3

  ! ! boundary conditions parameters
  ! #Include "../include/bcDefineFortranInclude.h"

  integer dirichlet,neumann,evenSymmetry,radiation,exactBC,abcEM2,characteristic,absorbing
  parameter( dirichlet=1, neumann=2, evenSymmetry=3, radiation=4, exactBC=5, abcEM2=6, characteristic=7, absorbing=8  )    

  ! DO THIS FOR NOW 
  integer symmetryBoundaryCondition
  parameter( symmetryBoundaryCondition=evenSymmetry )  

  ! forcing options
      ! forcingOptions -- these should match ForcingEnum in Maxwell.h 
      integer noForcing,magneticSinusoidalPointSource,gaussianSource,twilightZoneForcing,	gaussianChargeSource, userDefinedForcingOption
	integer noBoundaryForcing,planeWaveBoundaryForcing,chirpedPlaneWaveBoundaryForcing
      parameter(noForcing                =0,magneticSinusoidalPointSource =1,gaussianSource                =2,twilightZoneForcing           =3,	   gaussianChargeSource          =4,userDefinedForcingOption      =5 )
      ! boundary forcing options when solved directly for the scattered field:
      parameter( noBoundaryForcing              =0,		 planeWaveBoundaryForcing       =1,chirpedPlaneWaveBoundaryForcing=2 )

  integer method,nfdtd,bamx
  parameter( nfdtd=5,bamx=7 ) 

  integer rectangular,curvilinear
  parameter(rectangular=0,curvilinear=1)


  ! ---bcOpt variables ---
  integer addForcingBC,assignBCForImplicit,assignKnownSolutionAtBoundaries,bcApproach,gridIsImplicit,knownSolutionOption
  integer numberOfComponents,numberOfFrequencies,numberOfProcessors,numGhost,twilightZone,useUpwindDissipation
  integer uc,ghost
  real REAL_MIN
  logical checkResiduals

!     --- start statement function ----
  integer kd,m,n
  real rx,ry,rz,sx,sy,sz,tx,ty,tz

   real d12
   real d22
   real h12
   real h22
   real rxr2
   real rxs2
   real rxt2
   real rxrr2
   real rxss2
   real rxrs2
   real ryr2
   real rys2
   real ryt2
   real ryrr2
   real ryss2
   real ryrs2
   real rzr2
   real rzs2
   real rzt2
   real rzrr2
   real rzss2
   real rzrs2
   real sxr2
   real sxs2
   real sxt2
   real sxrr2
   real sxss2
   real sxrs2
   real syr2
   real sys2
   real syt2
   real syrr2
   real syss2
   real syrs2
   real szr2
   real szs2
   real szt2
   real szrr2
   real szss2
   real szrs2
   real txr2
   real txs2
   real txt2
   real txrr2
   real txss2
   real txrs2
   real tyr2
   real tys2
   real tyt2
   real tyrr2
   real tyss2
   real tyrs2
   real tzr2
   real tzs2
   real tzt2
   real tzrr2
   real tzss2
   real tzrs2
   real rxx21
   real rxx22
   real rxy22
   real rxx23
   real rxy23
   real rxz23
   real ryx22
   real ryy22
   real ryx23
   real ryy23
   real ryz23
   real rzx22
   real rzy22
   real rzx23
   real rzy23
   real rzz23
   real sxx22
   real sxy22
   real sxx23
   real sxy23
   real sxz23
   real syx22
   real syy22
   real syx23
   real syy23
   real syz23
   real szx22
   real szy22
   real szx23
   real szy23
   real szz23
   real txx22
   real txy22
   real txx23
   real txy23
   real txz23
   real tyx22
   real tyy22
   real tyx23
   real tyy23
   real tyz23
   real tzx22
   real tzy22
   real tzx23
   real tzy23
   real tzz23
   real ur2
   real us2
   real ut2
   real urr2
   real uss2
   real urs2
   real utt2
   real urt2
   real ust2
   real urrr2
   real usss2
   real uttt2
   real ux21
   real uy21
   real uz21
   real ux22
   real uy22
   real uz22
   real ux23
   real uy23
   real uz23
   real uxx21
   real uyy21
   real uxy21
   real uxz21
   real uyz21
   real uzz21
   real ulaplacian21
   real uxx22
   real uyy22
   real uxy22
   real uxz22
   real uyz22
   real uzz22
   real ulaplacian22
   real uxx23
   real uyy23
   real uzz23
   real uxy23
   real uxz23
   real uyz23
   real ulaplacian23
   real ux23r
   real uy23r
   real uz23r
   real uxx23r
   real uyy23r
   real uxy23r
   real uzz23r
   real uxz23r
   real uyz23r
   real ux21r
   real uy21r
   real uz21r
   real uxx21r
   real uyy21r
   real uzz21r
   real uxy21r
   real uxz21r
   real uyz21r
   real ulaplacian21r
   real ux22r
   real uy22r
   real uz22r
   real uxx22r
   real uyy22r
   real uzz22r
   real uxy22r
   real uxz22r
   real uyz22r
   real ulaplacian22r
   real ulaplacian23r
   real uxxx22r
   real uyyy22r
   real uxxy22r
   real uxyy22r
   real uxxxx22r
   real uyyyy22r
   real uxxyy22r
   real uxxx23r
   real uyyy23r
   real uzzz23r
   real uxxy23r
   real uxxz23r
   real uxyy23r
   real uyyz23r
   real uxzz23r
   real uyzz23r
   real uxxxx23r
   real uyyyy23r
   real uzzzz23r
   real uxxyy23r
   real uxxzz23r
   real uyyzz23r
   real uLapSq22r
   real uLapSq23r
   real unr2
   real uns2
   real unt2
   real unrr2
   real unss2
   real unrs2
   real untt2
   real unrt2
   real unst2
   real unrrr2
   real unsss2
   real unttt2
   real unx21
   real uny21
   real unz21
   real unx22
   real uny22
   real unz22
   real unx23
   real uny23
   real unz23
   real unxx21
   real unyy21
   real unxy21
   real unxz21
   real unyz21
   real unzz21
   real unlaplacian21
   real unxx22
   real unyy22
   real unxy22
   real unxz22
   real unyz22
   real unzz22
   real unlaplacian22
   real unxx23
   real unyy23
   real unzz23
   real unxy23
   real unxz23
   real unyz23
   real unlaplacian23
   real unx23r
   real uny23r
   real unz23r
   real unxx23r
   real unyy23r
   real unxy23r
   real unzz23r
   real unxz23r
   real unyz23r
   real unx21r
   real uny21r
   real unz21r
   real unxx21r
   real unyy21r
   real unzz21r
   real unxy21r
   real unxz21r
   real unyz21r
   real unlaplacian21r
   real unx22r
   real uny22r
   real unz22r
   real unxx22r
   real unyy22r
   real unzz22r
   real unxy22r
   real unxz22r
   real unyz22r
   real unlaplacian22r
   real unlaplacian23r
   real unxxx22r
   real unyyy22r
   real unxxy22r
   real unxyy22r
   real unxxxx22r
   real unyyyy22r
   real unxxyy22r
   real unxxx23r
   real unyyy23r
   real unzzz23r
   real unxxy23r
   real unxxz23r
   real unxyy23r
   real unyyz23r
   real unxzz23r
   real unyzz23r
   real unxxxx23r
   real unyyyy23r
   real unzzzz23r
   real unxxyy23r
   real unxxzz23r
   real unyyzz23r
   real unLapSq22r
   real unLapSq23r

   real d14
   real d24
   real h41
   real h42
   real rxr4
   real rxs4
   real rxt4
   real ryr4
   real rys4
   real ryt4
   real rzr4
   real rzs4
   real rzt4
   real sxr4
   real sxs4
   real sxt4
   real syr4
   real sys4
   real syt4
   real szr4
   real szs4
   real szt4
   real txr4
   real txs4
   real txt4
   real tyr4
   real tys4
   real tyt4
   real tzr4
   real tzs4
   real tzt4
   real rxx41
   real rxx42
   real rxy42
   real rxx43
   real rxy43
   real rxz43
   real ryx42
   real ryy42
   real ryx43
   real ryy43
   real ryz43
   real rzx42
   real rzy42
   real rzx43
   real rzy43
   real rzz43
   real sxx42
   real sxy42
   real sxx43
   real sxy43
   real sxz43
   real syx42
   real syy42
   real syx43
   real syy43
   real syz43
   real szx42
   real szy42
   real szx43
   real szy43
   real szz43
   real txx42
   real txy42
   real txx43
   real txy43
   real txz43
   real tyx42
   real tyy42
   real tyx43
   real tyy43
   real tyz43
   real tzx42
   real tzy42
   real tzx43
   real tzy43
   real tzz43
   real ur4
   real us4
   real ut4
   real urr4
   real uss4
   real utt4
   real urs4
   real urt4
   real ust4
   real ux41
   real uy41
   real uz41
   real ux42
   real uy42
   real uz42
   real ux43
   real uy43
   real uz43
   real uxx41
   real uyy41
   real uxy41
   real uxz41
   real uyz41
   real uzz41
   real ulaplacian41
   real uxx42
   real uyy42
   real uxy42
   real uxz42
   real uyz42
   real uzz42
   real ulaplacian42
   real uxx43
   real uyy43
   real uzz43
   real uxy43
   real uxz43
   real uyz43
   real ulaplacian43
   real ux43r
   real uy43r
   real uz43r
   real uxx43r
   real uyy43r
   real uzz43r
   real uxy43r
   real uxz43r
   real uyz43r
   real ux41r
   real uy41r
   real uz41r
   real uxx41r
   real uyy41r
   real uzz41r
   real uxy41r
   real uxz41r
   real uyz41r
   real ulaplacian41r
   real ux42r
   real uy42r
   real uz42r
   real uxx42r
   real uyy42r
   real uzz42r
   real uxy42r
   real uxz42r
   real uyz42r
   real ulaplacian42r
   real ulaplacian43r
   real unr4
   real uns4
   real unt4
   real unrr4
   real unss4
   real untt4
   real unrs4
   real unrt4
   real unst4
   real unx41
   real uny41
   real unz41
   real unx42
   real uny42
   real unz42
   real unx43
   real uny43
   real unz43
   real unxx41
   real unyy41
   real unxy41
   real unxz41
   real unyz41
   real unzz41
   real unlaplacian41
   real unxx42
   real unyy42
   real unxy42
   real unxz42
   real unyz42
   real unzz42
   real unlaplacian42
   real unxx43
   real unyy43
   real unzz43
   real unxy43
   real unxz43
   real unyz43
   real unlaplacian43
   real unx43r
   real uny43r
   real unz43r
   real unxx43r
   real unyy43r
   real unzz43r
   real unxy43r
   real unxz43r
   real unyz43r
   real unx41r
   real uny41r
   real unz41r
   real unxx41r
   real unyy41r
   real unzz41r
   real unxy41r
   real unxz41r
   real unyz41r
   real unlaplacian41r
   real unx42r
   real uny42r
   real unz42r
   real unxx42r
   real unyy42r
   real unzz42r
   real unxy42r
   real unxz42r
   real unyz42r
   real unlaplacian42r
   real unlaplacian43r
      real rsxyr2,rsxys2,rsxyt2,rsxyx22,rsxyy22,rsxyr4,rsxys4,rsxyx42,rsxyy42
      real rsxyxs42, rsxyys42, rsxyxr42, rsxyyr42
      real rsxyrr2,rsxyss2,rsxyrs2, rsxyrr4,rsxyss4,rsxyrs4
     
      real rsxyx43,rsxyy43,rsxyz43,rsxyt4,rsxytt4,rsxyrt4,rsxyst4
      real rsxyxr43,rsxyxs43,rsxyxt43
      real rsxyyr43,rsxyys43,rsxyyt43
      real rsxyzr43,rsxyzs43,rsxyzt43
     
      real rsxyxr22,rsxyxs22,rsxyyr22,rsxyys22
      real rsxyx23,rsxyy23,rsxyz23,rsxytt2,rsxyrt2,rsxyst2
      real rsxyxr23,rsxyxs23,rsxyxt23
      real rsxyyr23,rsxyys23,rsxyyt23
      real rsxyzr23,rsxyzs23,rsxyzt23

!.......statement functions for jacobian
  rx(i1,i2,i3)=rsxy(i1,i2,i3,0,0)
  ry(i1,i2,i3)=rsxy(i1,i2,i3,0,1)
  rz(i1,i2,i3)=rsxy(i1,i2,i3,0,2)
  sx(i1,i2,i3)=rsxy(i1,i2,i3,1,0)
  sy(i1,i2,i3)=rsxy(i1,i2,i3,1,1)
  sz(i1,i2,i3)=rsxy(i1,i2,i3,1,2)
  tx(i1,i2,i3)=rsxy(i1,i2,i3,2,0)
  ty(i1,i2,i3)=rsxy(i1,i2,i3,2,1)
  tz(i1,i2,i3)=rsxy(i1,i2,i3,2,2)


!     The next macro call will define the difference approximation statement functions
  d12(kd) = 1./(2.*dr(kd))
  d22(kd) = 1./(dr(kd)**2)
  ur2(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*d12(0)
  us2(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*d12(1)
  ut2(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*d12(2)
  urr2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) )*d22(0)
  uss2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) )*d22(1)
  urs2(i1,i2,i3,kd)=(ur2(i1,i2+1,i3,kd)-ur2(i1,i2-1,i3,kd))*d12(1)
  utt2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd)) )*d22(2)
  urt2(i1,i2,i3,kd)=(ur2(i1,i2,i3+1,kd)-ur2(i1,i2,i3-1,kd))*d12(2)
  ust2(i1,i2,i3,kd)=(us2(i1,i2,i3+1,kd)-us2(i1,i2,i3-1,kd))*d12(2)
  urrr2(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
  usss2(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
  uttt2(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
  rxr2(i1,i2,i3)=(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))*d12(0)
  rxs2(i1,i2,i3)=(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))*d12(1)
  rxt2(i1,i2,i3)=(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))*d12(2)
  rxrr2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1+1,i2,i3)+rx(i1-1,i2,i3)) )*d22(0)
  rxss2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1,i2+1,i3)+rx(i1,i2-1,i3)) )*d22(1)
  rxrs2(i1,i2,i3)=(rxr2(i1,i2+1,i3)-rxr2(i1,i2-1,i3))*d12(1)
  ryr2(i1,i2,i3)=(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))*d12(0)
  rys2(i1,i2,i3)=(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))*d12(1)
  ryt2(i1,i2,i3)=(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))*d12(2)
  ryrr2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1+1,i2,i3)+ry(i1-1,i2,i3)) )*d22(0)
  ryss2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1,i2+1,i3)+ry(i1,i2-1,i3)) )*d22(1)
  ryrs2(i1,i2,i3)=(ryr2(i1,i2+1,i3)-ryr2(i1,i2-1,i3))*d12(1)
  rzr2(i1,i2,i3)=(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))*d12(0)
  rzs2(i1,i2,i3)=(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))*d12(1)
  rzt2(i1,i2,i3)=(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))*d12(2)
  rzrr2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1+1,i2,i3)+rz(i1-1,i2,i3)) )*d22(0)
  rzss2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1,i2+1,i3)+rz(i1,i2-1,i3)) )*d22(1)
  rzrs2(i1,i2,i3)=(rzr2(i1,i2+1,i3)-rzr2(i1,i2-1,i3))*d12(1)
  sxr2(i1,i2,i3)=(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))*d12(0)
  sxs2(i1,i2,i3)=(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))*d12(1)
  sxt2(i1,i2,i3)=(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))*d12(2)
  sxrr2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1+1,i2,i3)+sx(i1-1,i2,i3)) )*d22(0)
  sxss2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1,i2+1,i3)+sx(i1,i2-1,i3)) )*d22(1)
  sxrs2(i1,i2,i3)=(sxr2(i1,i2+1,i3)-sxr2(i1,i2-1,i3))*d12(1)
  syr2(i1,i2,i3)=(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))*d12(0)
  sys2(i1,i2,i3)=(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))*d12(1)
  syt2(i1,i2,i3)=(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))*d12(2)
  syrr2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1+1,i2,i3)+sy(i1-1,i2,i3)) )*d22(0)
  syss2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1,i2+1,i3)+sy(i1,i2-1,i3)) )*d22(1)
  syrs2(i1,i2,i3)=(syr2(i1,i2+1,i3)-syr2(i1,i2-1,i3))*d12(1)
  szr2(i1,i2,i3)=(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))*d12(0)
  szs2(i1,i2,i3)=(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))*d12(1)
  szt2(i1,i2,i3)=(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))*d12(2)
  szrr2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1+1,i2,i3)+sz(i1-1,i2,i3)) )*d22(0)
  szss2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1,i2+1,i3)+sz(i1,i2-1,i3)) )*d22(1)
  szrs2(i1,i2,i3)=(szr2(i1,i2+1,i3)-szr2(i1,i2-1,i3))*d12(1)
  txr2(i1,i2,i3)=(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))*d12(0)
  txs2(i1,i2,i3)=(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))*d12(1)
  txt2(i1,i2,i3)=(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))*d12(2)
  txrr2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1+1,i2,i3)+tx(i1-1,i2,i3)) )*d22(0)
  txss2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1,i2+1,i3)+tx(i1,i2-1,i3)) )*d22(1)
  txrs2(i1,i2,i3)=(txr2(i1,i2+1,i3)-txr2(i1,i2-1,i3))*d12(1)
  tyr2(i1,i2,i3)=(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))*d12(0)
  tys2(i1,i2,i3)=(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))*d12(1)
  tyt2(i1,i2,i3)=(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))*d12(2)
  tyrr2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1+1,i2,i3)+ty(i1-1,i2,i3)) )*d22(0)
  tyss2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1,i2+1,i3)+ty(i1,i2-1,i3)) )*d22(1)
  tyrs2(i1,i2,i3)=(tyr2(i1,i2+1,i3)-tyr2(i1,i2-1,i3))*d12(1)
  tzr2(i1,i2,i3)=(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))*d12(0)
  tzs2(i1,i2,i3)=(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))*d12(1)
  tzt2(i1,i2,i3)=(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))*d12(2)
  tzrr2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1+1,i2,i3)+tz(i1-1,i2,i3)) )*d22(0)
  tzss2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1,i2+1,i3)+tz(i1,i2-1,i3)) )*d22(1)
  tzrs2(i1,i2,i3)=(tzr2(i1,i2+1,i3)-tzr2(i1,i2-1,i3))*d12(1)
  ux21(i1,i2,i3,kd)= rx(i1,i2,i3)*ur2(i1,i2,i3,kd)
  uy21(i1,i2,i3,kd)=0
  uz21(i1,i2,i3,kd)=0
  ux22(i1,i2,i3,kd)= rx(i1,i2,i3)*ur2(i1,i2,i3,kd)+sx(i1,i2,i3)*us2(i1,i2,i3,kd)
  uy22(i1,i2,i3,kd)= ry(i1,i2,i3)*ur2(i1,i2,i3,kd)+sy(i1,i2,i3)*us2(i1,i2,i3,kd)
  uz22(i1,i2,i3,kd)=0
  ux23(i1,i2,i3,kd)=rx(i1,i2,i3)*ur2(i1,i2,i3,kd)+sx(i1,i2,i3)*us2(i1,i2,i3,kd)+tx(i1,i2,i3)*ut2(i1,i2,i3,kd)
  uy23(i1,i2,i3,kd)=ry(i1,i2,i3)*ur2(i1,i2,i3,kd)+sy(i1,i2,i3)*us2(i1,i2,i3,kd)+ty(i1,i2,i3)*ut2(i1,i2,i3,kd)
  uz23(i1,i2,i3,kd)=rz(i1,i2,i3)*ur2(i1,i2,i3,kd)+sz(i1,i2,i3)*us2(i1,i2,i3,kd)+tz(i1,i2,i3)*ut2(i1,i2,i3,kd)
  rxx21(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)
  rxx22(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(i1,i2,i3)
  rxy22(i1,i2,i3)= ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(i1,i2,i3)
  rxx23(i1,i2,i3)=rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(i1,i2,i3)+tx(i1,i2,i3)*rxt2(i1,i2,i3)
  rxy23(i1,i2,i3)=ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(i1,i2,i3)+ty(i1,i2,i3)*rxt2(i1,i2,i3)
  rxz23(i1,i2,i3)=rz(i1,i2,i3)*rxr2(i1,i2,i3)+sz(i1,i2,i3)*rxs2(i1,i2,i3)+tz(i1,i2,i3)*rxt2(i1,i2,i3)
  ryx22(i1,i2,i3)= rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(i1,i2,i3)
  ryy22(i1,i2,i3)= ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(i1,i2,i3)
  ryx23(i1,i2,i3)=rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(i1,i2,i3)+tx(i1,i2,i3)*ryt2(i1,i2,i3)
  ryy23(i1,i2,i3)=ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(i1,i2,i3)+ty(i1,i2,i3)*ryt2(i1,i2,i3)
  ryz23(i1,i2,i3)=rz(i1,i2,i3)*ryr2(i1,i2,i3)+sz(i1,i2,i3)*rys2(i1,i2,i3)+tz(i1,i2,i3)*ryt2(i1,i2,i3)
  rzx22(i1,i2,i3)= rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(i1,i2,i3)
  rzy22(i1,i2,i3)= ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(i1,i2,i3)
  rzx23(i1,i2,i3)=rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(i1,i2,i3)+tx(i1,i2,i3)*rzt2(i1,i2,i3)
  rzy23(i1,i2,i3)=ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(i1,i2,i3)+ty(i1,i2,i3)*rzt2(i1,i2,i3)
  rzz23(i1,i2,i3)=rz(i1,i2,i3)*rzr2(i1,i2,i3)+sz(i1,i2,i3)*rzs2(i1,i2,i3)+tz(i1,i2,i3)*rzt2(i1,i2,i3)
  sxx22(i1,i2,i3)= rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(i1,i2,i3)
  sxy22(i1,i2,i3)= ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(i1,i2,i3)
  sxx23(i1,i2,i3)=rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(i1,i2,i3)+tx(i1,i2,i3)*sxt2(i1,i2,i3)
  sxy23(i1,i2,i3)=ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(i1,i2,i3)+ty(i1,i2,i3)*sxt2(i1,i2,i3)
  sxz23(i1,i2,i3)=rz(i1,i2,i3)*sxr2(i1,i2,i3)+sz(i1,i2,i3)*sxs2(i1,i2,i3)+tz(i1,i2,i3)*sxt2(i1,i2,i3)
  syx22(i1,i2,i3)= rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(i1,i2,i3)
  syy22(i1,i2,i3)= ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(i1,i2,i3)
  syx23(i1,i2,i3)=rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(i1,i2,i3)+tx(i1,i2,i3)*syt2(i1,i2,i3)
  syy23(i1,i2,i3)=ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(i1,i2,i3)+ty(i1,i2,i3)*syt2(i1,i2,i3)
  syz23(i1,i2,i3)=rz(i1,i2,i3)*syr2(i1,i2,i3)+sz(i1,i2,i3)*sys2(i1,i2,i3)+tz(i1,i2,i3)*syt2(i1,i2,i3)
  szx22(i1,i2,i3)= rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(i1,i2,i3)
  szy22(i1,i2,i3)= ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(i1,i2,i3)
  szx23(i1,i2,i3)=rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(i1,i2,i3)+tx(i1,i2,i3)*szt2(i1,i2,i3)
  szy23(i1,i2,i3)=ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(i1,i2,i3)+ty(i1,i2,i3)*szt2(i1,i2,i3)
  szz23(i1,i2,i3)=rz(i1,i2,i3)*szr2(i1,i2,i3)+sz(i1,i2,i3)*szs2(i1,i2,i3)+tz(i1,i2,i3)*szt2(i1,i2,i3)
  txx22(i1,i2,i3)= rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(i1,i2,i3)
  txy22(i1,i2,i3)= ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(i1,i2,i3)
  txx23(i1,i2,i3)=rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(i1,i2,i3)+tx(i1,i2,i3)*txt2(i1,i2,i3)
  txy23(i1,i2,i3)=ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(i1,i2,i3)+ty(i1,i2,i3)*txt2(i1,i2,i3)
  txz23(i1,i2,i3)=rz(i1,i2,i3)*txr2(i1,i2,i3)+sz(i1,i2,i3)*txs2(i1,i2,i3)+tz(i1,i2,i3)*txt2(i1,i2,i3)
  tyx22(i1,i2,i3)= rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(i1,i2,i3)
  tyy22(i1,i2,i3)= ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(i1,i2,i3)
  tyx23(i1,i2,i3)=rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(i1,i2,i3)+tx(i1,i2,i3)*tyt2(i1,i2,i3)
  tyy23(i1,i2,i3)=ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(i1,i2,i3)+ty(i1,i2,i3)*tyt2(i1,i2,i3)
  tyz23(i1,i2,i3)=rz(i1,i2,i3)*tyr2(i1,i2,i3)+sz(i1,i2,i3)*tys2(i1,i2,i3)+tz(i1,i2,i3)*tyt2(i1,i2,i3)
  tzx22(i1,i2,i3)= rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*tzs2(i1,i2,i3)
  tzy22(i1,i2,i3)= ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*tzs2(i1,i2,i3)
  tzx23(i1,i2,i3)=rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*tzs2(i1,i2,i3)+tx(i1,i2,i3)*tzt2(i1,i2,i3)
  tzy23(i1,i2,i3)=ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*tzs2(i1,i2,i3)+ty(i1,i2,i3)*tzt2(i1,i2,i3)
  tzz23(i1,i2,i3)=rz(i1,i2,i3)*tzr2(i1,i2,i3)+sz(i1,i2,i3)*tzs2(i1,i2,i3)+tz(i1,i2,i3)*tzt2(i1,i2,i3)
  uxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*ur2(i1,i2,i3,kd)
  uyy21(i1,i2,i3,kd)=0
  uxy21(i1,i2,i3,kd)=0
  uxz21(i1,i2,i3,kd)=0
  uyz21(i1,i2,i3,kd)=0
  uzz21(i1,i2,i3,kd)=0
  ulaplacian21(i1,i2,i3,kd)=uxx21(i1,i2,i3,kd)
  uxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx22(i1,i2,i3))*us2(i1,i2,i3,kd)
  uyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(syy22(i1,i2,i3))*us2(i1,i2,i3,kd)
  uxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss2(i1,i2,i3,kd)+rxy22(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*us2(i1,i2,i3,kd)
  uxz22(i1,i2,i3,kd)=0
  uyz22(i1,i2,i3,kd)=0
  uzz22(i1,i2,i3,kd)=0
  ulaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*us2(i1,i2,i3,kd)
  uxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sx(i1,i2,i3)**2*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*us2(i1,i2,i3,kd)+txx23(i1,i2,i3)*ut2(i1,i2,i3,kd)
  uyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sy(i1,i2,i3)**2*uss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*ur2(i1,i2,i3,kd)+syy23(i1,i2,i3)*us2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*ut2(i1,i2,i3,kd)
  uzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sz(i1,i2,i3)**2*uss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+szz23(i1,i2,i3)*us2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
  uxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*utt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust2(i1,i2,i3,kd)+rxy23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*us2(i1,i2,i3,kd)+txy23(i1,i2,i3)*ut2(i1,i2,i3,kd)
  uxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*utt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust2(i1,i2,i3,kd)+rxz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*us2(i1,i2,i3,kd)+txz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
  uyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr2(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*uss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*utt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust2(i1,i2,i3,kd)+ryz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+syz23(i1,i2,i3)*us2(i1,i2,i3,kd)+tyz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
  ulaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*utt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*urs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*urt2(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*ust2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*us2(i1,i2,i3,kd)+(txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*ut2(i1,i2,i3,kd)
  !============================================================================================
  ! Define derivatives for a rectangular grid
  !
  !============================================================================================
  h12(kd) = 1./(2.*dx(kd))
  h22(kd) = 1./(dx(kd)**2)
  ux23r(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*h12(0)
  uy23r(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*h12(1)
  uz23r(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*h12(2)
  uxx23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) )*h22(0)
  uyy23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) )*h22(1)
  uxy23r(i1,i2,i3,kd)=(ux23r(i1,i2+1,i3,kd)-ux23r(i1,i2-1,i3,kd))*h12(1)
  uzz23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd)) )*h22(2)
  uxz23r(i1,i2,i3,kd)=(ux23r(i1,i2,i3+1,kd)-ux23r(i1,i2,i3-1,kd))*h12(2)
  uyz23r(i1,i2,i3,kd)=(uy23r(i1,i2,i3+1,kd)-uy23r(i1,i2,i3-1,kd))*h12(2)
  ux21r(i1,i2,i3,kd)= ux23r(i1,i2,i3,kd)
  uy21r(i1,i2,i3,kd)= uy23r(i1,i2,i3,kd)
  uz21r(i1,i2,i3,kd)= uz23r(i1,i2,i3,kd)
  uxx21r(i1,i2,i3,kd)= uxx23r(i1,i2,i3,kd)
  uyy21r(i1,i2,i3,kd)= uyy23r(i1,i2,i3,kd)
  uzz21r(i1,i2,i3,kd)= uzz23r(i1,i2,i3,kd)
  uxy21r(i1,i2,i3,kd)= uxy23r(i1,i2,i3,kd)
  uxz21r(i1,i2,i3,kd)= uxz23r(i1,i2,i3,kd)
  uyz21r(i1,i2,i3,kd)= uyz23r(i1,i2,i3,kd)
  ulaplacian21r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)
  ux22r(i1,i2,i3,kd)= ux23r(i1,i2,i3,kd)
  uy22r(i1,i2,i3,kd)= uy23r(i1,i2,i3,kd)
  uz22r(i1,i2,i3,kd)= uz23r(i1,i2,i3,kd)
  uxx22r(i1,i2,i3,kd)= uxx23r(i1,i2,i3,kd)
  uyy22r(i1,i2,i3,kd)= uyy23r(i1,i2,i3,kd)
  uzz22r(i1,i2,i3,kd)= uzz23r(i1,i2,i3,kd)
  uxy22r(i1,i2,i3,kd)= uxy23r(i1,i2,i3,kd)
  uxz22r(i1,i2,i3,kd)= uxz23r(i1,i2,i3,kd)
  uyz22r(i1,i2,i3,kd)= uyz23r(i1,i2,i3,kd)
  ulaplacian22r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)+uyy23r(i1,i2,i3,kd)
  ulaplacian23r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)+uyy23r(i1,i2,i3,kd)+uzz23r(i1,i2,i3,kd)
  uxxx22r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
  uyyy22r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
  uxxy22r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
  uxyy22r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
  uxxxx22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4)
  uyyyy22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)
  uxxyy22r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
  ! 2D laplacian squared = u.xxxx + 2 u.xxyy + u.yyyy
  uLapSq22r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +2.*(u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
  uxxx23r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
  uyyy23r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
  uzzz23r(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
  uxxy23r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
  uxyy23r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
  uxxz23r(i1,i2,i3,kd)=( uxx22r(i1,i2,i3+1,kd)-uxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
  uyyz23r(i1,i2,i3,kd)=( uyy22r(i1,i2,i3+1,kd)-uyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
  uxzz23r(i1,i2,i3,kd)=( uzz22r(i1+1,i2,i3,kd)-uzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
  uyzz23r(i1,i2,i3,kd)=( uzz22r(i1,i2+1,i3,kd)-uzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
  uxxxx23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4)
  uyyyy23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)
  uzzzz23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)**4)
  uxxyy23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
  uxxzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))   +   (u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
  uyyzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1,i2+1,i3,kd)  +u(i1,i2-1,i3,kd)+  u(i1,i2  ,i3+1,kd)+u(i1,i2  ,i3-1,kd))   +   (u(i1,i2+1,i3+1,kd)+u(i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
  ! 3D laplacian squared = u.xxxx + u.yyyy + u.zzzz + 2 (u.xxyy + u.xxzz + u.yyzz )
  uLapSq23r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))    +(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)  +u(i1-1,i2,i3,kd)  +u(i1  ,i2+1,i3,kd)+u(i1  ,i2-1,i3,kd))   +2.*(u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)  +u(i1-1,i2,i3,kd)  +u(i1  ,i2,i3+1,kd)+u(i1  ,i2,i3-1,kd))   +2.*(u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*u(i1,i2,i3,kd)     -4.*(u(i1,i2+1,i3,kd)  +u(i1,i2-1,i3,kd)  +u(i1,i2  ,i3+1,kd)+u(i1,i2  ,i3-1,kd))   +2.*(u(i1,i2+1,i3+1,kd)+u(i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
  unr2(i1,i2,i3,kd)=(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))*d12(0)
  uns2(i1,i2,i3,kd)=(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))*d12(1)
  unt2(i1,i2,i3,kd)=(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))*d12(2)
  unrr2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)) )*d22(0)
  unss2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd)) )*d22(1)
  unrs2(i1,i2,i3,kd)=(unr2(i1,i2+1,i3,kd)-unr2(i1,i2-1,i3,kd))*d12(1)
  untt2(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd)) )*d22(2)
  unrt2(i1,i2,i3,kd)=(unr2(i1,i2,i3+1,kd)-unr2(i1,i2,i3-1,kd))*d12(2)
  unst2(i1,i2,i3,kd)=(uns2(i1,i2,i3+1,kd)-uns2(i1,i2,i3-1,kd))*d12(2)
  unrrr2(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
  unsss2(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
  unttt2(i1,i2,i3,kd)=(-2.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))+(un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
  unx21(i1,i2,i3,kd)= rx(i1,i2,i3)*unr2(i1,i2,i3,kd)
  uny21(i1,i2,i3,kd)=0
  unz21(i1,i2,i3,kd)=0
  unx22(i1,i2,i3,kd)= rx(i1,i2,i3)*unr2(i1,i2,i3,kd)+sx(i1,i2,i3)*uns2(i1,i2,i3,kd)
  uny22(i1,i2,i3,kd)= ry(i1,i2,i3)*unr2(i1,i2,i3,kd)+sy(i1,i2,i3)*uns2(i1,i2,i3,kd)
  unz22(i1,i2,i3,kd)=0
  unx23(i1,i2,i3,kd)=rx(i1,i2,i3)*unr2(i1,i2,i3,kd)+sx(i1,i2,i3)*uns2(i1,i2,i3,kd)+tx(i1,i2,i3)*unt2(i1,i2,i3,kd)
  uny23(i1,i2,i3,kd)=ry(i1,i2,i3)*unr2(i1,i2,i3,kd)+sy(i1,i2,i3)*uns2(i1,i2,i3,kd)+ty(i1,i2,i3)*unt2(i1,i2,i3,kd)
  unz23(i1,i2,i3,kd)=rz(i1,i2,i3)*unr2(i1,i2,i3,kd)+sz(i1,i2,i3)*uns2(i1,i2,i3,kd)+tz(i1,i2,i3)*unt2(i1,i2,i3,kd)
  unxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*unr2(i1,i2,i3,kd)
  unyy21(i1,i2,i3,kd)=0
  unxy21(i1,i2,i3,kd)=0
  unxz21(i1,i2,i3,kd)=0
  unyz21(i1,i2,i3,kd)=0
  unzz21(i1,i2,i3,kd)=0
  unlaplacian21(i1,i2,i3,kd)=unxx21(i1,i2,i3,kd)
  unxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*unss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*unr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3))*uns2(i1,i2,i3,kd)
  unyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*unss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*unr2(i1,i2,i3,kd)+(syy22(i1,i2,i3))*uns2(i1,i2,i3,kd)
  unxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss2(i1,i2,i3,kd)+rxy22(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*uns2(i1,i2,i3,kd)
  unxz22(i1,i2,i3,kd)=0
  unyz22(i1,i2,i3,kd)=0
  unzz22(i1,i2,i3,kd)=0
  unlaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*unss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*unr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*uns2(i1,i2,i3,kd)
  unxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*unrr2(i1,i2,i3,kd)+sx(i1,i2,i3)**2*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*untt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*unrs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*unrt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*unst2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*uns2(i1,i2,i3,kd)+txx23(i1,i2,i3)*unt2(i1,i2,i3,kd)
  unyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*unrr2(i1,i2,i3,kd)+sy(i1,i2,i3)**2*unss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*untt2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*unrs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*unrt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*unst2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*unr2(i1,i2,i3,kd)+syy23(i1,i2,i3)*uns2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*unt2(i1,i2,i3,kd)
  unzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*unrr2(i1,i2,i3,kd)+sz(i1,i2,i3)**2*unss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*untt2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*unrs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*unrt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*unst2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*unr2(i1,i2,i3,kd)+szz23(i1,i2,i3)*uns2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*unt2(i1,i2,i3,kd)
  unxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*untt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*unst2(i1,i2,i3,kd)+rxy23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*uns2(i1,i2,i3,kd)+txy23(i1,i2,i3)*unt2(i1,i2,i3,kd)
  unxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*unrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*unss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*untt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*unst2(i1,i2,i3,kd)+rxz23(i1,i2,i3)*unr2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*uns2(i1,i2,i3,kd)+txz23(i1,i2,i3)*unt2(i1,i2,i3,kd)
  unyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*unrr2(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*unss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*untt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*unrs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*unrt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*unst2(i1,i2,i3,kd)+ryz23(i1,i2,i3)*unr2(i1,i2,i3,kd)+syz23(i1,i2,i3)*uns2(i1,i2,i3,kd)+tyz23(i1,i2,i3)*unt2(i1,i2,i3,kd)
  unlaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*unrr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*unss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*untt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*unrs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*unrt2(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*unst2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*unr2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*uns2(i1,i2,i3,kd)+(txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*unt2(i1,i2,i3,kd)
  !============================================================================================
  ! Define derivatives for a rectangular grid
  !
  !============================================================================================
  unx23r(i1,i2,i3,kd)=(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))*h12(0)
  uny23r(i1,i2,i3,kd)=(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))*h12(1)
  unz23r(i1,i2,i3,kd)=(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))*h12(2)
  unxx23r(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)) )*h22(0)
  unyy23r(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd)) )*h22(1)
  unxy23r(i1,i2,i3,kd)=(unx23r(i1,i2+1,i3,kd)-unx23r(i1,i2-1,i3,kd))*h12(1)
  unzz23r(i1,i2,i3,kd)=(-2.*un(i1,i2,i3,kd)+(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd)) )*h22(2)
  unxz23r(i1,i2,i3,kd)=(unx23r(i1,i2,i3+1,kd)-unx23r(i1,i2,i3-1,kd))*h12(2)
  unyz23r(i1,i2,i3,kd)=(uny23r(i1,i2,i3+1,kd)-uny23r(i1,i2,i3-1,kd))*h12(2)
  unx21r(i1,i2,i3,kd)= unx23r(i1,i2,i3,kd)
  uny21r(i1,i2,i3,kd)= uny23r(i1,i2,i3,kd)
  unz21r(i1,i2,i3,kd)= unz23r(i1,i2,i3,kd)
  unxx21r(i1,i2,i3,kd)= unxx23r(i1,i2,i3,kd)
  unyy21r(i1,i2,i3,kd)= unyy23r(i1,i2,i3,kd)
  unzz21r(i1,i2,i3,kd)= unzz23r(i1,i2,i3,kd)
  unxy21r(i1,i2,i3,kd)= unxy23r(i1,i2,i3,kd)
  unxz21r(i1,i2,i3,kd)= unxz23r(i1,i2,i3,kd)
  unyz21r(i1,i2,i3,kd)= unyz23r(i1,i2,i3,kd)
  unlaplacian21r(i1,i2,i3,kd)=unxx23r(i1,i2,i3,kd)
  unx22r(i1,i2,i3,kd)= unx23r(i1,i2,i3,kd)
  uny22r(i1,i2,i3,kd)= uny23r(i1,i2,i3,kd)
  unz22r(i1,i2,i3,kd)= unz23r(i1,i2,i3,kd)
  unxx22r(i1,i2,i3,kd)= unxx23r(i1,i2,i3,kd)
  unyy22r(i1,i2,i3,kd)= unyy23r(i1,i2,i3,kd)
  unzz22r(i1,i2,i3,kd)= unzz23r(i1,i2,i3,kd)
  unxy22r(i1,i2,i3,kd)= unxy23r(i1,i2,i3,kd)
  unxz22r(i1,i2,i3,kd)= unxz23r(i1,i2,i3,kd)
  unyz22r(i1,i2,i3,kd)= unyz23r(i1,i2,i3,kd)
  unlaplacian22r(i1,i2,i3,kd)=unxx23r(i1,i2,i3,kd)+unyy23r(i1,i2,i3,kd)
  unlaplacian23r(i1,i2,i3,kd)=unxx23r(i1,i2,i3,kd)+unyy23r(i1,i2,i3,kd)+unzz23r(i1,i2,i3,kd)
  unxxx22r(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
  unyyy22r(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
  unxxy22r(i1,i2,i3,kd)=( unxx22r(i1,i2+1,i3,kd)-unxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
  unxyy22r(i1,i2,i3,kd)=( unyy22r(i1+1,i2,i3,kd)-unyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
  unxxxx22r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(dx(0)**4)
  unyyyy22r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(dx(1)**4)
  unxxyy22r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))   +   (un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
  ! 2D laplacian squared = un.xxxx + 2 un.xxyy + un.yyyy
  unLapSq22r(i1,i2,i3,kd)= ( 6.*un(i1,i2,i3,kd)   - 4.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))    +(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))    +(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 8.*un(i1,i2,i3,kd)     -4.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))   +2.*(un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
  unxxx23r(i1,i2,i3,kd)=(-2.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
  unyyy23r(i1,i2,i3,kd)=(-2.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
  unzzz23r(i1,i2,i3,kd)=(-2.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))+(un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
  unxxy23r(i1,i2,i3,kd)=( unxx22r(i1,i2+1,i3,kd)-unxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
  unxyy23r(i1,i2,i3,kd)=( unyy22r(i1+1,i2,i3,kd)-unyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
  unxxz23r(i1,i2,i3,kd)=( unxx22r(i1,i2,i3+1,kd)-unxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
  unyyz23r(i1,i2,i3,kd)=( unyy22r(i1,i2,i3+1,kd)-unyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
  unxzz23r(i1,i2,i3,kd)=( unzz22r(i1+1,i2,i3,kd)-unzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
  unyzz23r(i1,i2,i3,kd)=( unzz22r(i1,i2+1,i3,kd)-unzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
  unxxxx23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))+(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(dx(0)**4)
  unyyyy23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))+(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(dx(1)**4)
  unzzzz23r(i1,i2,i3,kd)=(6.*un(i1,i2,i3,kd)-4.*(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))+(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )/(dx(2)**4)
  unxxyy23r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))   +   (un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
  unxxzz23r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd)+un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))   +   (un(i1+1,i2,i3+1,kd)+un(i1-1,i2,i3+1,kd)+un(i1+1,i2,i3-1,kd)+un(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
  unyyzz23r(i1,i2,i3,kd)=( 4.*un(i1,i2,i3,kd)     -2.*(un(i1,i2+1,i3,kd)  +un(i1,i2-1,i3,kd)+  un(i1,i2  ,i3+1,kd)+un(i1,i2  ,i3-1,kd))   +   (un(i1,i2+1,i3+1,kd)+un(i1,i2-1,i3+1,kd)+un(i1,i2+1,i3-1,kd)+un(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
  ! 3D laplacian squared = un.xxxx + un.yyyy + un.zzzz + 2 (un.xxyy + un.xxzz + un.yyzz )
  unLapSq23r(i1,i2,i3,kd)= ( 6.*un(i1,i2,i3,kd)   - 4.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))    +(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))    +(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 6.*un(i1,i2,i3,kd)    -4.*(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))    +(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 8.*un(i1,i2,i3,kd)     -4.*(un(i1+1,i2,i3,kd)  +un(i1-1,i2,i3,kd)  +un(i1  ,i2+1,i3,kd)+un(i1  ,i2-1,i3,kd))   +2.*(un(i1+1,i2+1,i3,kd)+un(i1-1,i2+1,i3,kd)+un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*un(i1,i2,i3,kd)     -4.*(un(i1+1,i2,i3,kd)  +un(i1-1,i2,i3,kd)  +un(i1  ,i2,i3+1,kd)+un(i1  ,i2,i3-1,kd))   +2.*(un(i1+1,i2,i3+1,kd)+un(i1-1,i2,i3+1,kd)+un(i1+1,i2,i3-1,kd)+un(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*un(i1,i2,i3,kd)     -4.*(un(i1,i2+1,i3,kd)  +un(i1,i2-1,i3,kd)  +un(i1,i2  ,i3+1,kd)+un(i1,i2  ,i3-1,kd))   +2.*(un(i1,i2+1,i3+1,kd)+un(i1,i2-1,i3+1,kd)+un(i1,i2+1,i3-1,kd)+un(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
  d14(kd) = 1./(12.*dr(kd))
  d24(kd) = 1./(12.*dr(kd)**2)
  ur4(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*d14(0)
  us4(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*d14(1)
  ut4(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*d14(2)
  urr4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*d24(0)
  uss4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*d24(1)
  utt4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*d24(2)
  urs4(i1,i2,i3,kd)=(8.*(ur4(i1,i2+1,i3,kd)-ur4(i1,i2-1,i3,kd))-(ur4(i1,i2+2,i3,kd)-ur4(i1,i2-2,i3,kd)))*d14(1)
  urt4(i1,i2,i3,kd)=(8.*(ur4(i1,i2,i3+1,kd)-ur4(i1,i2,i3-1,kd))-(ur4(i1,i2,i3+2,kd)-ur4(i1,i2,i3-2,kd)))*d14(2)
  ust4(i1,i2,i3,kd)=(8.*(us4(i1,i2,i3+1,kd)-us4(i1,i2,i3-1,kd))-(us4(i1,i2,i3+2,kd)-us4(i1,i2,i3-2,kd)))*d14(2)
  rxr4(i1,i2,i3)=(8.*(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))-(rx(i1+2,i2,i3)-rx(i1-2,i2,i3)))*d14(0)
  rxs4(i1,i2,i3)=(8.*(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))-(rx(i1,i2+2,i3)-rx(i1,i2-2,i3)))*d14(1)
  rxt4(i1,i2,i3)=(8.*(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))-(rx(i1,i2,i3+2)-rx(i1,i2,i3-2)))*d14(2)
  ryr4(i1,i2,i3)=(8.*(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))-(ry(i1+2,i2,i3)-ry(i1-2,i2,i3)))*d14(0)
  rys4(i1,i2,i3)=(8.*(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))-(ry(i1,i2+2,i3)-ry(i1,i2-2,i3)))*d14(1)
  ryt4(i1,i2,i3)=(8.*(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))-(ry(i1,i2,i3+2)-ry(i1,i2,i3-2)))*d14(2)
  rzr4(i1,i2,i3)=(8.*(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))-(rz(i1+2,i2,i3)-rz(i1-2,i2,i3)))*d14(0)
  rzs4(i1,i2,i3)=(8.*(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))-(rz(i1,i2+2,i3)-rz(i1,i2-2,i3)))*d14(1)
  rzt4(i1,i2,i3)=(8.*(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))-(rz(i1,i2,i3+2)-rz(i1,i2,i3-2)))*d14(2)
  sxr4(i1,i2,i3)=(8.*(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))-(sx(i1+2,i2,i3)-sx(i1-2,i2,i3)))*d14(0)
  sxs4(i1,i2,i3)=(8.*(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))-(sx(i1,i2+2,i3)-sx(i1,i2-2,i3)))*d14(1)
  sxt4(i1,i2,i3)=(8.*(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))-(sx(i1,i2,i3+2)-sx(i1,i2,i3-2)))*d14(2)
  syr4(i1,i2,i3)=(8.*(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))-(sy(i1+2,i2,i3)-sy(i1-2,i2,i3)))*d14(0)
  sys4(i1,i2,i3)=(8.*(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))-(sy(i1,i2+2,i3)-sy(i1,i2-2,i3)))*d14(1)
  syt4(i1,i2,i3)=(8.*(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))-(sy(i1,i2,i3+2)-sy(i1,i2,i3-2)))*d14(2)
  szr4(i1,i2,i3)=(8.*(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))-(sz(i1+2,i2,i3)-sz(i1-2,i2,i3)))*d14(0)
  szs4(i1,i2,i3)=(8.*(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))-(sz(i1,i2+2,i3)-sz(i1,i2-2,i3)))*d14(1)
  szt4(i1,i2,i3)=(8.*(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))-(sz(i1,i2,i3+2)-sz(i1,i2,i3-2)))*d14(2)
  txr4(i1,i2,i3)=(8.*(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))-(tx(i1+2,i2,i3)-tx(i1-2,i2,i3)))*d14(0)
  txs4(i1,i2,i3)=(8.*(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))-(tx(i1,i2+2,i3)-tx(i1,i2-2,i3)))*d14(1)
  txt4(i1,i2,i3)=(8.*(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))-(tx(i1,i2,i3+2)-tx(i1,i2,i3-2)))*d14(2)
  tyr4(i1,i2,i3)=(8.*(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))-(ty(i1+2,i2,i3)-ty(i1-2,i2,i3)))*d14(0)
  tys4(i1,i2,i3)=(8.*(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))-(ty(i1,i2+2,i3)-ty(i1,i2-2,i3)))*d14(1)
  tyt4(i1,i2,i3)=(8.*(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))-(ty(i1,i2,i3+2)-ty(i1,i2,i3-2)))*d14(2)
  tzr4(i1,i2,i3)=(8.*(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))-(tz(i1+2,i2,i3)-tz(i1-2,i2,i3)))*d14(0)
  tzs4(i1,i2,i3)=(8.*(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))-(tz(i1,i2+2,i3)-tz(i1,i2-2,i3)))*d14(1)
  tzt4(i1,i2,i3)=(8.*(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))-(tz(i1,i2,i3+2)-tz(i1,i2,i3-2)))*d14(2)
  ux41(i1,i2,i3,kd)= rx(i1,i2,i3)*ur4(i1,i2,i3,kd)
  uy41(i1,i2,i3,kd)=0
  uz41(i1,i2,i3,kd)=0
  ux42(i1,i2,i3,kd)= rx(i1,i2,i3)*ur4(i1,i2,i3,kd)+sx(i1,i2,i3)*us4(i1,i2,i3,kd)
  uy42(i1,i2,i3,kd)= ry(i1,i2,i3)*ur4(i1,i2,i3,kd)+sy(i1,i2,i3)*us4(i1,i2,i3,kd)
  uz42(i1,i2,i3,kd)=0
  ux43(i1,i2,i3,kd)=rx(i1,i2,i3)*ur4(i1,i2,i3,kd)+sx(i1,i2,i3)*us4(i1,i2,i3,kd)+tx(i1,i2,i3)*ut4(i1,i2,i3,kd)
  uy43(i1,i2,i3,kd)=ry(i1,i2,i3)*ur4(i1,i2,i3,kd)+sy(i1,i2,i3)*us4(i1,i2,i3,kd)+ty(i1,i2,i3)*ut4(i1,i2,i3,kd)
  uz43(i1,i2,i3,kd)=rz(i1,i2,i3)*ur4(i1,i2,i3,kd)+sz(i1,i2,i3)*us4(i1,i2,i3,kd)+tz(i1,i2,i3)*ut4(i1,i2,i3,kd)
  rxx41(i1,i2,i3)= rx(i1,i2,i3)*rxr4(i1,i2,i3)
  rxx42(i1,i2,i3)= rx(i1,i2,i3)*rxr4(i1,i2,i3)+sx(i1,i2,i3)*rxs4(i1,i2,i3)
  rxy42(i1,i2,i3)= ry(i1,i2,i3)*rxr4(i1,i2,i3)+sy(i1,i2,i3)*rxs4(i1,i2,i3)
  rxx43(i1,i2,i3)=rx(i1,i2,i3)*rxr4(i1,i2,i3)+sx(i1,i2,i3)*rxs4(i1,i2,i3)+tx(i1,i2,i3)*rxt4(i1,i2,i3)
  rxy43(i1,i2,i3)=ry(i1,i2,i3)*rxr4(i1,i2,i3)+sy(i1,i2,i3)*rxs4(i1,i2,i3)+ty(i1,i2,i3)*rxt4(i1,i2,i3)
  rxz43(i1,i2,i3)=rz(i1,i2,i3)*rxr4(i1,i2,i3)+sz(i1,i2,i3)*rxs4(i1,i2,i3)+tz(i1,i2,i3)*rxt4(i1,i2,i3)
  ryx42(i1,i2,i3)= rx(i1,i2,i3)*ryr4(i1,i2,i3)+sx(i1,i2,i3)*rys4(i1,i2,i3)
  ryy42(i1,i2,i3)= ry(i1,i2,i3)*ryr4(i1,i2,i3)+sy(i1,i2,i3)*rys4(i1,i2,i3)
  ryx43(i1,i2,i3)=rx(i1,i2,i3)*ryr4(i1,i2,i3)+sx(i1,i2,i3)*rys4(i1,i2,i3)+tx(i1,i2,i3)*ryt4(i1,i2,i3)
  ryy43(i1,i2,i3)=ry(i1,i2,i3)*ryr4(i1,i2,i3)+sy(i1,i2,i3)*rys4(i1,i2,i3)+ty(i1,i2,i3)*ryt4(i1,i2,i3)
  ryz43(i1,i2,i3)=rz(i1,i2,i3)*ryr4(i1,i2,i3)+sz(i1,i2,i3)*rys4(i1,i2,i3)+tz(i1,i2,i3)*ryt4(i1,i2,i3)
  rzx42(i1,i2,i3)= rx(i1,i2,i3)*rzr4(i1,i2,i3)+sx(i1,i2,i3)*rzs4(i1,i2,i3)
  rzy42(i1,i2,i3)= ry(i1,i2,i3)*rzr4(i1,i2,i3)+sy(i1,i2,i3)*rzs4(i1,i2,i3)
  rzx43(i1,i2,i3)=rx(i1,i2,i3)*rzr4(i1,i2,i3)+sx(i1,i2,i3)*rzs4(i1,i2,i3)+tx(i1,i2,i3)*rzt4(i1,i2,i3)
  rzy43(i1,i2,i3)=ry(i1,i2,i3)*rzr4(i1,i2,i3)+sy(i1,i2,i3)*rzs4(i1,i2,i3)+ty(i1,i2,i3)*rzt4(i1,i2,i3)
  rzz43(i1,i2,i3)=rz(i1,i2,i3)*rzr4(i1,i2,i3)+sz(i1,i2,i3)*rzs4(i1,i2,i3)+tz(i1,i2,i3)*rzt4(i1,i2,i3)
  sxx42(i1,i2,i3)= rx(i1,i2,i3)*sxr4(i1,i2,i3)+sx(i1,i2,i3)*sxs4(i1,i2,i3)
  sxy42(i1,i2,i3)= ry(i1,i2,i3)*sxr4(i1,i2,i3)+sy(i1,i2,i3)*sxs4(i1,i2,i3)
  sxx43(i1,i2,i3)=rx(i1,i2,i3)*sxr4(i1,i2,i3)+sx(i1,i2,i3)*sxs4(i1,i2,i3)+tx(i1,i2,i3)*sxt4(i1,i2,i3)
  sxy43(i1,i2,i3)=ry(i1,i2,i3)*sxr4(i1,i2,i3)+sy(i1,i2,i3)*sxs4(i1,i2,i3)+ty(i1,i2,i3)*sxt4(i1,i2,i3)
  sxz43(i1,i2,i3)=rz(i1,i2,i3)*sxr4(i1,i2,i3)+sz(i1,i2,i3)*sxs4(i1,i2,i3)+tz(i1,i2,i3)*sxt4(i1,i2,i3)
  syx42(i1,i2,i3)= rx(i1,i2,i3)*syr4(i1,i2,i3)+sx(i1,i2,i3)*sys4(i1,i2,i3)
  syy42(i1,i2,i3)= ry(i1,i2,i3)*syr4(i1,i2,i3)+sy(i1,i2,i3)*sys4(i1,i2,i3)
  syx43(i1,i2,i3)=rx(i1,i2,i3)*syr4(i1,i2,i3)+sx(i1,i2,i3)*sys4(i1,i2,i3)+tx(i1,i2,i3)*syt4(i1,i2,i3)
  syy43(i1,i2,i3)=ry(i1,i2,i3)*syr4(i1,i2,i3)+sy(i1,i2,i3)*sys4(i1,i2,i3)+ty(i1,i2,i3)*syt4(i1,i2,i3)
  syz43(i1,i2,i3)=rz(i1,i2,i3)*syr4(i1,i2,i3)+sz(i1,i2,i3)*sys4(i1,i2,i3)+tz(i1,i2,i3)*syt4(i1,i2,i3)
  szx42(i1,i2,i3)= rx(i1,i2,i3)*szr4(i1,i2,i3)+sx(i1,i2,i3)*szs4(i1,i2,i3)
  szy42(i1,i2,i3)= ry(i1,i2,i3)*szr4(i1,i2,i3)+sy(i1,i2,i3)*szs4(i1,i2,i3)
  szx43(i1,i2,i3)=rx(i1,i2,i3)*szr4(i1,i2,i3)+sx(i1,i2,i3)*szs4(i1,i2,i3)+tx(i1,i2,i3)*szt4(i1,i2,i3)
  szy43(i1,i2,i3)=ry(i1,i2,i3)*szr4(i1,i2,i3)+sy(i1,i2,i3)*szs4(i1,i2,i3)+ty(i1,i2,i3)*szt4(i1,i2,i3)
  szz43(i1,i2,i3)=rz(i1,i2,i3)*szr4(i1,i2,i3)+sz(i1,i2,i3)*szs4(i1,i2,i3)+tz(i1,i2,i3)*szt4(i1,i2,i3)
  txx42(i1,i2,i3)= rx(i1,i2,i3)*txr4(i1,i2,i3)+sx(i1,i2,i3)*txs4(i1,i2,i3)
  txy42(i1,i2,i3)= ry(i1,i2,i3)*txr4(i1,i2,i3)+sy(i1,i2,i3)*txs4(i1,i2,i3)
  txx43(i1,i2,i3)=rx(i1,i2,i3)*txr4(i1,i2,i3)+sx(i1,i2,i3)*txs4(i1,i2,i3)+tx(i1,i2,i3)*txt4(i1,i2,i3)
  txy43(i1,i2,i3)=ry(i1,i2,i3)*txr4(i1,i2,i3)+sy(i1,i2,i3)*txs4(i1,i2,i3)+ty(i1,i2,i3)*txt4(i1,i2,i3)
  txz43(i1,i2,i3)=rz(i1,i2,i3)*txr4(i1,i2,i3)+sz(i1,i2,i3)*txs4(i1,i2,i3)+tz(i1,i2,i3)*txt4(i1,i2,i3)
  tyx42(i1,i2,i3)= rx(i1,i2,i3)*tyr4(i1,i2,i3)+sx(i1,i2,i3)*tys4(i1,i2,i3)
  tyy42(i1,i2,i3)= ry(i1,i2,i3)*tyr4(i1,i2,i3)+sy(i1,i2,i3)*tys4(i1,i2,i3)
  tyx43(i1,i2,i3)=rx(i1,i2,i3)*tyr4(i1,i2,i3)+sx(i1,i2,i3)*tys4(i1,i2,i3)+tx(i1,i2,i3)*tyt4(i1,i2,i3)
  tyy43(i1,i2,i3)=ry(i1,i2,i3)*tyr4(i1,i2,i3)+sy(i1,i2,i3)*tys4(i1,i2,i3)+ty(i1,i2,i3)*tyt4(i1,i2,i3)
  tyz43(i1,i2,i3)=rz(i1,i2,i3)*tyr4(i1,i2,i3)+sz(i1,i2,i3)*tys4(i1,i2,i3)+tz(i1,i2,i3)*tyt4(i1,i2,i3)
  tzx42(i1,i2,i3)= rx(i1,i2,i3)*tzr4(i1,i2,i3)+sx(i1,i2,i3)*tzs4(i1,i2,i3)
  tzy42(i1,i2,i3)= ry(i1,i2,i3)*tzr4(i1,i2,i3)+sy(i1,i2,i3)*tzs4(i1,i2,i3)
  tzx43(i1,i2,i3)=rx(i1,i2,i3)*tzr4(i1,i2,i3)+sx(i1,i2,i3)*tzs4(i1,i2,i3)+tx(i1,i2,i3)*tzt4(i1,i2,i3)
  tzy43(i1,i2,i3)=ry(i1,i2,i3)*tzr4(i1,i2,i3)+sy(i1,i2,i3)*tzs4(i1,i2,i3)+ty(i1,i2,i3)*tzt4(i1,i2,i3)
  tzz43(i1,i2,i3)=rz(i1,i2,i3)*tzr4(i1,i2,i3)+sz(i1,i2,i3)*tzs4(i1,i2,i3)+tz(i1,i2,i3)*tzt4(i1,i2,i3)
  uxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*ur4(i1,i2,i3,kd)
  uyy41(i1,i2,i3,kd)=0
  uxy41(i1,i2,i3,kd)=0
  uxz41(i1,i2,i3,kd)=0
  uyz41(i1,i2,i3,kd)=0
  uzz41(i1,i2,i3,kd)=0
  ulaplacian41(i1,i2,i3,kd)=uxx41(i1,i2,i3,kd)
  uxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx42(i1,i2,i3))*us4(i1,i2,i3,kd)
  uyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(syy42(i1,i2,i3))*us4(i1,i2,i3,kd)
  uxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss4(i1,i2,i3,kd)+rxy42(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*us4(i1,i2,i3,kd)
  uxz42(i1,i2,i3,kd)=0
  uyz42(i1,i2,i3,kd)=0
  uzz42(i1,i2,i3,kd)=0
  ulaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*us4(i1,i2,i3,kd)
  uxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sx(i1,i2,i3)**2*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust4(i1,i2,i3,kd)+rxx43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxx43(i1,i2,i3)*us4(i1,i2,i3,kd)+txx43(i1,i2,i3)*ut4(i1,i2,i3,kd)
  uyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sy(i1,i2,i3)**2*uss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust4(i1,i2,i3,kd)+ryy43(i1,i2,i3)*ur4(i1,i2,i3,kd)+syy43(i1,i2,i3)*us4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*ut4(i1,i2,i3,kd)
  uzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sz(i1,i2,i3)**2*uss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust4(i1,i2,i3,kd)+rzz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+szz43(i1,i2,i3)*us4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
  uxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*utt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust4(i1,i2,i3,kd)+rxy43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*us4(i1,i2,i3,kd)+txy43(i1,i2,i3)*ut4(i1,i2,i3,kd)
  uxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*utt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust4(i1,i2,i3,kd)+rxz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*us4(i1,i2,i3,kd)+txz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
  uyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr4(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*uss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*utt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust4(i1,i2,i3,kd)+ryz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+syz43(i1,i2,i3)*us4(i1,i2,i3,kd)+tyz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
  ulaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*utt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*urs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*urt4(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*ust4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+ryy43(i1,i2,i3)+rzz43(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx43(i1,i2,i3)+syy43(i1,i2,i3)+szz43(i1,i2,i3))*us4(i1,i2,i3,kd)+(txx43(i1,i2,i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*ut4(i1,i2,i3,kd)
  !============================================================================================
  ! Define derivatives for a rectangular grid
  !
  !============================================================================================
  h41(kd) = 1./(12.*dx(kd))
  h42(kd) = 1./(12.*dx(kd)**2)
  ux43r(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*h41(0)
  uy43r(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*h41(1)
  uz43r(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*h41(2)
  uxx43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*h42(0) 
  uyy43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*h42(1) 
  uzz43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*h42(2)
  uxy43r(i1,i2,i3,kd)=( (u(i1+2,i2+2,i3,kd)-u(i1-2,i2+2,i3,kd)- u(i1+2,i2-2,i3,kd)+u(i1-2,i2-2,i3,kd)) +8.*(u(i1-1,i2+2,i3,kd)-u(i1-1,i2-2,i3,kd)-u(i1+1,i2+2,i3,kd)+u(i1+1,i2-2,i3,kd) +u(i1+2,i2-1,i3,kd)-u(i1-2,i2-1,i3,kd)-u(i1+2,i2+1,i3,kd)+u(i1-2,i2+1,i3,kd))+64.*(u(i1+1,i2+1,i3,kd)-u(i1-1,i2+1,i3,kd)- u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
  uxz43r(i1,i2,i3,kd)=( (u(i1+2,i2,i3+2,kd)-u(i1-2,i2,i3+2,kd)-u(i1+2,i2,i3-2,kd)+u(i1-2,i2,i3-2,kd)) +8.*(u(i1-1,i2,i3+2,kd)-u(i1-1,i2,i3-2,kd)-u(i1+1,i2,i3+2,kd)+u(i1+1,i2,i3-2,kd) +u(i1+2,i2,i3-1,kd)-u(i1-2,i2,i3-1,kd)- u(i1+2,i2,i3+1,kd)+u(i1-2,i2,i3+1,kd)) +64.*(u(i1+1,i2,i3+1,kd)-u(i1-1,i2,i3+1,kd)-u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
  uyz43r(i1,i2,i3,kd)=( (u(i1,i2+2,i3+2,kd)-u(i1,i2-2,i3+2,kd)-u(i1,i2+2,i3-2,kd)+u(i1,i2-2,i3-2,kd)) +8.*(u(i1,i2-1,i3+2,kd)-u(i1,i2-1,i3-2,kd)-u(i1,i2+1,i3+2,kd)+u(i1,i2+1,i3-2,kd) +u(i1,i2+2,i3-1,kd)-u(i1,i2-2,i3-1,kd)-u(i1,i2+2,i3+1,kd)+u(i1,i2-2,i3+1,kd)) +64.*(u(i1,i2+1,i3+1,kd)-u(i1,i2-1,i3+1,kd)-u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
  ux41r(i1,i2,i3,kd)= ux43r(i1,i2,i3,kd)
  uy41r(i1,i2,i3,kd)= uy43r(i1,i2,i3,kd)
  uz41r(i1,i2,i3,kd)= uz43r(i1,i2,i3,kd)
  uxx41r(i1,i2,i3,kd)= uxx43r(i1,i2,i3,kd)
  uyy41r(i1,i2,i3,kd)= uyy43r(i1,i2,i3,kd)
  uzz41r(i1,i2,i3,kd)= uzz43r(i1,i2,i3,kd)
  uxy41r(i1,i2,i3,kd)= uxy43r(i1,i2,i3,kd)
  uxz41r(i1,i2,i3,kd)= uxz43r(i1,i2,i3,kd)
  uyz41r(i1,i2,i3,kd)= uyz43r(i1,i2,i3,kd)
  ulaplacian41r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)
  ux42r(i1,i2,i3,kd)= ux43r(i1,i2,i3,kd)
  uy42r(i1,i2,i3,kd)= uy43r(i1,i2,i3,kd)
  uz42r(i1,i2,i3,kd)= uz43r(i1,i2,i3,kd)
  uxx42r(i1,i2,i3,kd)= uxx43r(i1,i2,i3,kd)
  uyy42r(i1,i2,i3,kd)= uyy43r(i1,i2,i3,kd)
  uzz42r(i1,i2,i3,kd)= uzz43r(i1,i2,i3,kd)
  uxy42r(i1,i2,i3,kd)= uxy43r(i1,i2,i3,kd)
  uxz42r(i1,i2,i3,kd)= uxz43r(i1,i2,i3,kd)
  uyz42r(i1,i2,i3,kd)= uyz43r(i1,i2,i3,kd)
  ulaplacian42r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)+uyy43r(i1,i2,i3,kd)
  ulaplacian43r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)+uyy43r(i1,i2,i3,kd)+uzz43r(i1,i2,i3,kd)
  unr4(i1,i2,i3,kd)=(8.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))-(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)))*d14(0)
  uns4(i1,i2,i3,kd)=(8.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))-(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)))*d14(1)
  unt4(i1,i2,i3,kd)=(8.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))-(un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)))*d14(2)
  unrr4(i1,i2,i3,kd)=(-30.*un(i1,i2,i3,kd)+16.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))-(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )*d24(0)
  unss4(i1,i2,i3,kd)=(-30.*un(i1,i2,i3,kd)+16.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))-(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )*d24(1)
  untt4(i1,i2,i3,kd)=(-30.*un(i1,i2,i3,kd)+16.*(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))-(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )*d24(2)
  unrs4(i1,i2,i3,kd)=(8.*(unr4(i1,i2+1,i3,kd)-unr4(i1,i2-1,i3,kd))-(unr4(i1,i2+2,i3,kd)-unr4(i1,i2-2,i3,kd)))*d14(1)
  unrt4(i1,i2,i3,kd)=(8.*(unr4(i1,i2,i3+1,kd)-unr4(i1,i2,i3-1,kd))-(unr4(i1,i2,i3+2,kd)-unr4(i1,i2,i3-2,kd)))*d14(2)
  unst4(i1,i2,i3,kd)=(8.*(uns4(i1,i2,i3+1,kd)-uns4(i1,i2,i3-1,kd))-(uns4(i1,i2,i3+2,kd)-uns4(i1,i2,i3-2,kd)))*d14(2)
  unx41(i1,i2,i3,kd)= rx(i1,i2,i3)*unr4(i1,i2,i3,kd)
  uny41(i1,i2,i3,kd)=0
  unz41(i1,i2,i3,kd)=0
  unx42(i1,i2,i3,kd)= rx(i1,i2,i3)*unr4(i1,i2,i3,kd)+sx(i1,i2,i3)*uns4(i1,i2,i3,kd)
  uny42(i1,i2,i3,kd)= ry(i1,i2,i3)*unr4(i1,i2,i3,kd)+sy(i1,i2,i3)*uns4(i1,i2,i3,kd)
  unz42(i1,i2,i3,kd)=0
  unx43(i1,i2,i3,kd)=rx(i1,i2,i3)*unr4(i1,i2,i3,kd)+sx(i1,i2,i3)*uns4(i1,i2,i3,kd)+tx(i1,i2,i3)*unt4(i1,i2,i3,kd)
  uny43(i1,i2,i3,kd)=ry(i1,i2,i3)*unr4(i1,i2,i3,kd)+sy(i1,i2,i3)*uns4(i1,i2,i3,kd)+ty(i1,i2,i3)*unt4(i1,i2,i3,kd)
  unz43(i1,i2,i3,kd)=rz(i1,i2,i3)*unr4(i1,i2,i3,kd)+sz(i1,i2,i3)*uns4(i1,i2,i3,kd)+tz(i1,i2,i3)*unt4(i1,i2,i3,kd)
  unxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*unr4(i1,i2,i3,kd)
  unyy41(i1,i2,i3,kd)=0
  unxy41(i1,i2,i3,kd)=0
  unxz41(i1,i2,i3,kd)=0
  unyz41(i1,i2,i3,kd)=0
  unzz41(i1,i2,i3,kd)=0
  unlaplacian41(i1,i2,i3,kd)=unxx41(i1,i2,i3,kd)
  unxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*unss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*unr4(i1,i2,i3,kd)+(sxx42(i1,i2,i3))*uns4(i1,i2,i3,kd)
  unyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*unss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*unr4(i1,i2,i3,kd)+(syy42(i1,i2,i3))*uns4(i1,i2,i3,kd)
  unxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss4(i1,i2,i3,kd)+rxy42(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*uns4(i1,i2,i3,kd)
  unxz42(i1,i2,i3,kd)=0
  unyz42(i1,i2,i3,kd)=0
  unzz42(i1,i2,i3,kd)=0
  unlaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*unss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*unr4(i1,i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*uns4(i1,i2,i3,kd)
  unxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*unrr4(i1,i2,i3,kd)+sx(i1,i2,i3)**2*unss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*untt4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*unrs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*unrt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*unst4(i1,i2,i3,kd)+rxx43(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxx43(i1,i2,i3)*uns4(i1,i2,i3,kd)+txx43(i1,i2,i3)*unt4(i1,i2,i3,kd)
  unyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*unrr4(i1,i2,i3,kd)+sy(i1,i2,i3)**2*unss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*untt4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*unrs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*unrt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*unst4(i1,i2,i3,kd)+ryy43(i1,i2,i3)*unr4(i1,i2,i3,kd)+syy43(i1,i2,i3)*uns4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*unt4(i1,i2,i3,kd)
  unzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*unrr4(i1,i2,i3,kd)+sz(i1,i2,i3)**2*unss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*untt4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*unrs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*unrt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*unst4(i1,i2,i3,kd)+rzz43(i1,i2,i3)*unr4(i1,i2,i3,kd)+szz43(i1,i2,i3)*uns4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*unt4(i1,i2,i3,kd)
  unxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*unrr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*unss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*untt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*unrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*unst4(i1,i2,i3,kd)+rxy43(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*uns4(i1,i2,i3,kd)+txy43(i1,i2,i3)*unt4(i1,i2,i3,kd)
  unxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*unrr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*unss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*untt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*unrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*unst4(i1,i2,i3,kd)+rxz43(i1,i2,i3)*unr4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*uns4(i1,i2,i3,kd)+txz43(i1,i2,i3)*unt4(i1,i2,i3,kd)
  unyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*unrr4(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*unss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*untt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*unrs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*unrt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*unst4(i1,i2,i3,kd)+ryz43(i1,i2,i3)*unr4(i1,i2,i3,kd)+syz43(i1,i2,i3)*uns4(i1,i2,i3,kd)+tyz43(i1,i2,i3)*unt4(i1,i2,i3,kd)
  unlaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*unrr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*unss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*untt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*unrs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*unrt4(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*unst4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+ryy43(i1,i2,i3)+rzz43(i1,i2,i3))*unr4(i1,i2,i3,kd)+(sxx43(i1,i2,i3)+syy43(i1,i2,i3)+szz43(i1,i2,i3))*uns4(i1,i2,i3,kd)+(txx43(i1,i2,i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*unt4(i1,i2,i3,kd)
  !============================================================================================
  ! Define derivatives for a rectangular grid
  !
  !============================================================================================
  unx43r(i1,i2,i3,kd)=(8.*(un(i1+1,i2,i3,kd)-un(i1-1,i2,i3,kd))-(un(i1+2,i2,i3,kd)-un(i1-2,i2,i3,kd)))*h41(0)
  uny43r(i1,i2,i3,kd)=(8.*(un(i1,i2+1,i3,kd)-un(i1,i2-1,i3,kd))-(un(i1,i2+2,i3,kd)-un(i1,i2-2,i3,kd)))*h41(1)
  unz43r(i1,i2,i3,kd)=(8.*(un(i1,i2,i3+1,kd)-un(i1,i2,i3-1,kd))-(un(i1,i2,i3+2,kd)-un(i1,i2,i3-2,kd)))*h41(2)
  unxx43r(i1,i2,i3,kd)=( -30.*un(i1,i2,i3,kd)+16.*(un(i1+1,i2,i3,kd)+un(i1-1,i2,i3,kd))-(un(i1+2,i2,i3,kd)+un(i1-2,i2,i3,kd)) )*h42(0) 
  unyy43r(i1,i2,i3,kd)=( -30.*un(i1,i2,i3,kd)+16.*(un(i1,i2+1,i3,kd)+un(i1,i2-1,i3,kd))-(un(i1,i2+2,i3,kd)+un(i1,i2-2,i3,kd)) )*h42(1) 
  unzz43r(i1,i2,i3,kd)=( -30.*un(i1,i2,i3,kd)+16.*(un(i1,i2,i3+1,kd)+un(i1,i2,i3-1,kd))-(un(i1,i2,i3+2,kd)+un(i1,i2,i3-2,kd)) )*h42(2)
  unxy43r(i1,i2,i3,kd)=( (un(i1+2,i2+2,i3,kd)-un(i1-2,i2+2,i3,kd)- un(i1+2,i2-2,i3,kd)+un(i1-2,i2-2,i3,kd)) +8.*(un(i1-1,i2+2,i3,kd)-un(i1-1,i2-2,i3,kd)-un(i1+1,i2+2,i3,kd)+un(i1+1,i2-2,i3,kd) +un(i1+2,i2-1,i3,kd)-un(i1-2,i2-1,i3,kd)-un(i1+2,i2+1,i3,kd)+un(i1-2,i2+1,i3,kd))+64.*(un(i1+1,i2+1,i3,kd)-un(i1-1,i2+1,i3,kd)- un(i1+1,i2-1,i3,kd)+un(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
  unxz43r(i1,i2,i3,kd)=( (un(i1+2,i2,i3+2,kd)-un(i1-2,i2,i3+2,kd)-un(i1+2,i2,i3-2,kd)+un(i1-2,i2,i3-2,kd)) +8.*(un(i1-1,i2,i3+2,kd)-un(i1-1,i2,i3-2,kd)-un(i1+1,i2,i3+2,kd)+un(i1+1,i2,i3-2,kd) +un(i1+2,i2,i3-1,kd)-un(i1-2,i2,i3-1,kd)- un(i1+2,i2,i3+1,kd)+un(i1-2,i2,i3+1,kd)) +64.*(un(i1+1,i2,i3+1,kd)-un(i1-1,i2,i3+1,kd)-un(i1+1,i2,i3-1,kd)+un(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
  unyz43r(i1,i2,i3,kd)=( (un(i1,i2+2,i3+2,kd)-un(i1,i2-2,i3+2,kd)-un(i1,i2+2,i3-2,kd)+un(i1,i2-2,i3-2,kd)) +8.*(un(i1,i2-1,i3+2,kd)-un(i1,i2-1,i3-2,kd)-un(i1,i2+1,i3+2,kd)+un(i1,i2+1,i3-2,kd) +un(i1,i2+2,i3-1,kd)-un(i1,i2-2,i3-1,kd)-un(i1,i2+2,i3+1,kd)+un(i1,i2-2,i3+1,kd)) +64.*(un(i1,i2+1,i3+1,kd)-un(i1,i2-1,i3+1,kd)-un(i1,i2+1,i3-1,kd)+un(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
  unx41r(i1,i2,i3,kd)= unx43r(i1,i2,i3,kd)
  uny41r(i1,i2,i3,kd)= uny43r(i1,i2,i3,kd)
  unz41r(i1,i2,i3,kd)= unz43r(i1,i2,i3,kd)
  unxx41r(i1,i2,i3,kd)= unxx43r(i1,i2,i3,kd)
  unyy41r(i1,i2,i3,kd)= unyy43r(i1,i2,i3,kd)
  unzz41r(i1,i2,i3,kd)= unzz43r(i1,i2,i3,kd)
  unxy41r(i1,i2,i3,kd)= unxy43r(i1,i2,i3,kd)
  unxz41r(i1,i2,i3,kd)= unxz43r(i1,i2,i3,kd)
  unyz41r(i1,i2,i3,kd)= unyz43r(i1,i2,i3,kd)
  unlaplacian41r(i1,i2,i3,kd)=unxx43r(i1,i2,i3,kd)
  unx42r(i1,i2,i3,kd)= unx43r(i1,i2,i3,kd)
  uny42r(i1,i2,i3,kd)= uny43r(i1,i2,i3,kd)
  unz42r(i1,i2,i3,kd)= unz43r(i1,i2,i3,kd)
  unxx42r(i1,i2,i3,kd)= unxx43r(i1,i2,i3,kd)
  unyy42r(i1,i2,i3,kd)= unyy43r(i1,i2,i3,kd)
  unzz42r(i1,i2,i3,kd)= unzz43r(i1,i2,i3,kd)
  unxy42r(i1,i2,i3,kd)= unxy43r(i1,i2,i3,kd)
  unxz42r(i1,i2,i3,kd)= unxz43r(i1,i2,i3,kd)
  unyz42r(i1,i2,i3,kd)= unyz43r(i1,i2,i3,kd)
  unlaplacian42r(i1,i2,i3,kd)=unxx43r(i1,i2,i3,kd)+unyy43r(i1,i2,i3,kd)
  unlaplacian43r(i1,i2,i3,kd)=unxx43r(i1,i2,i3,kd)+unyy43r(i1,i2,i3,kd)+unzz43r(i1,i2,i3,kd)

      rsxyr2(i1,i2,i3,m,n)=(rsxy(i1+1,i2,i3,m,n)-rsxy(i1-1,i2,i3,m,n))*d12(0)
      rsxys2(i1,i2,i3,m,n)=(rsxy(i1,i2+1,i3,m,n)-rsxy(i1,i2-1,i3,m,n))*d12(1)
      rsxyt2(i1,i2,i3,m,n)=(rsxy(i1,i2,i3+1,m,n)-rsxy(i1,i2,i3-1,m,n))*d12(2)
     
      rsxyrr2(i1,i2,i3,m,n)=(rsxy(i1+1,i2,i3,m,n)-2.*rsxy(i1,i2,i3,m,n)+rsxy(i1-1,i2,i3,m,n))*d22(0)
      rsxyss2(i1,i2,i3,m,n)=(rsxy(i1,i2+1,i3,m,n)-2.*rsxy(i1,i2,i3,m,n)+rsxy(i1,i2-1,i3,m,n))*d22(1)
      rsxytt2(i1,i2,i3,m,n)=(rsxy(i1,i2,i3+1,m,n)-2.*rsxy(i1,i2,i3,m,n)+rsxy(i1,i2,i3-1,m,n))*d22(2)
     
      rsxyx22(i1,i2,i3,m,n)= rx(i1,i2,i3)*rsxyr2(i1,i2,i3,m,n)+sx(i1,i2,i3)*rsxys2(i1,i2,i3,m,n)
      rsxyy22(i1,i2,i3,m,n)= ry(i1,i2,i3)*rsxyr2(i1,i2,i3,m,n)+sy(i1,i2,i3)*rsxys2(i1,i2,i3,m,n)
     
     
      ! check these again:
      !  -- 2nd -order ---
     
      rsxyrs2(i1,i2,i3,m,n)=(rsxyr2(i1,i2+1,i3,m,n)-rsxyr2(i1,i2-1,i3,m,n))*d12(1)
      rsxyrt2(i1,i2,i3,m,n)=(rsxyr2(i1,i2,i3+1,m,n)-rsxyr2(i1,i2,i3-1,m,n))*d12(2)
      rsxyst2(i1,i2,i3,m,n)=(rsxys2(i1,i2,i3+1,m,n)-rsxys2(i1,i2,i3-1,m,n))*d12(2)
     
      rsxyxr22(i1,i2,i3,m,n)= rsxyr2(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n) + rx(i1,i2,i3)*rsxyrr2(i1,i2,i3,m,n)+rsxyr2(i1,i2,i3,1,0)*rsxys2(i1,i2,i3,m,n) + sx(i1,i2,i3)*rsxyrs2(i1,i2,i3,m,n)
      rsxyxs22(i1,i2,i3,m,n)= rsxys2(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n) + rx(i1,i2,i3)*rsxyrs2(i1,i2,i3,m,n)+rsxys2(i1,i2,i3,1,0)*rsxys2(i1,i2,i3,m,n) + sx(i1,i2,i3)*rsxyss2(i1,i2,i3,m,n)
     
      rsxyyr22(i1,i2,i3,m,n)= rsxyr2(i1,i2,i3,0,1)*rsxyr2(i1,i2,i3,m,n) + ry(i1,i2,i3)*rsxyrr2(i1,i2,i3,m,n)+rsxyr2(i1,i2,i3,1,1)*rsxys2(i1,i2,i3,m,n) + sy(i1,i2,i3)*rsxyrs2(i1,i2,i3,m,n)
      rsxyys22(i1,i2,i3,m,n)= rsxys2(i1,i2,i3,0,1)*rsxyr2(i1,i2,i3,m,n) + ry(i1,i2,i3)*rsxyrs2(i1,i2,i3,m,n)+rsxys2(i1,i2,i3,1,1)*rsxys2(i1,i2,i3,m,n) + sy(i1,i2,i3)*rsxyss2(i1,i2,i3,m,n)
     
      ! 3d versions -- check these again
      rsxyx23(i1,i2,i3,m,n)= rx(i1,i2,i3)*rsxyr2(i1,i2,i3,m,n)+sx(i1,i2,i3)*rsxys2(i1,i2,i3,m,n)+tx(i1,i2,i3)*rsxyt2(i1,i2,i3,m,n)
      rsxyy23(i1,i2,i3,m,n)= ry(i1,i2,i3)*rsxyr2(i1,i2,i3,m,n)+sy(i1,i2,i3)*rsxys2(i1,i2,i3,m,n)+ty(i1,i2,i3)*rsxyt2(i1,i2,i3,m,n)
      rsxyz23(i1,i2,i3,m,n)= rz(i1,i2,i3)*rsxyr2(i1,i2,i3,m,n)+sz(i1,i2,i3)*rsxys2(i1,i2,i3,m,n)+tz(i1,i2,i3)*rsxyt2(i1,i2,i3,m,n)
     
      rsxyxr23(i1,i2,i3,m,n)= rsxyr2(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n) + rx(i1,i2,i3)*rsxyrr2(i1,i2,i3,m,n)+rsxyr2(i1,i2,i3,1,0)*rsxys2(i1,i2,i3,m,n) + sx(i1,i2,i3)*rsxyrs2(i1,i2,i3,m,n)+rsxyr2(i1,i2,i3,2,0)*rsxyt2(i1,i2,i3,m,n) + tx(i1,i2,i3)*rsxyrt2(i1,i2,i3,m,n)
     
      rsxyxs23(i1,i2,i3,m,n)= rsxys2(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n) + rx(i1,i2,i3)*rsxyrs2(i1,i2,i3,m,n)+rsxys2(i1,i2,i3,1,0)*rsxys2(i1,i2,i3,m,n) + sx(i1,i2,i3)*rsxyss2(i1,i2,i3,m,n)+rsxys2(i1,i2,i3,2,0)*rsxyt2(i1,i2,i3,m,n) + tx(i1,i2,i3)*rsxyst2(i1,i2,i3,m,n)
     
      rsxyxt23(i1,i2,i3,m,n)= rsxyt2(i1,i2,i3,0,0)*rsxyr2(i1,i2,i3,m,n) + rx(i1,i2,i3)*rsxyrt2(i1,i2,i3,m,n)+rsxyt2(i1,i2,i3,1,0)*rsxys2(i1,i2,i3,m,n) + sx(i1,i2,i3)*rsxyst2(i1,i2,i3,m,n)+rsxyt2(i1,i2,i3,2,0)*rsxyt2(i1,i2,i3,m,n) + tx(i1,i2,i3)*rsxytt2(i1,i2,i3,m,n)
     
      rsxyyr23(i1,i2,i3,m,n)= rsxyr2(i1,i2,i3,0,1)*rsxyr2(i1,i2,i3,m,n) + ry(i1,i2,i3)*rsxyrr2(i1,i2,i3,m,n)+rsxyr2(i1,i2,i3,1,1)*rsxys2(i1,i2,i3,m,n) + sy(i1,i2,i3)*rsxyrs2(i1,i2,i3,m,n)+rsxyr2(i1,i2,i3,2,1)*rsxyt2(i1,i2,i3,m,n) + ty(i1,i2,i3)*rsxyrt2(i1,i2,i3,m,n)
     
      rsxyys23(i1,i2,i3,m,n)= rsxys2(i1,i2,i3,0,1)*rsxyr2(i1,i2,i3,m,n) + ry(i1,i2,i3)*rsxyrs2(i1,i2,i3,m,n)+rsxys2(i1,i2,i3,1,1)*rsxys2(i1,i2,i3,m,n) + sy(i1,i2,i3)*rsxyss2(i1,i2,i3,m,n)+rsxys2(i1,i2,i3,2,1)*rsxyt2(i1,i2,i3,m,n) + ty(i1,i2,i3)*rsxyst2(i1,i2,i3,m,n)
     
      rsxyyt23(i1,i2,i3,m,n)= rsxyt2(i1,i2,i3,0,1)*rsxyr2(i1,i2,i3,m,n) + ry(i1,i2,i3)*rsxyrt2(i1,i2,i3,m,n)+rsxyt2(i1,i2,i3,1,1)*rsxys2(i1,i2,i3,m,n) + sy(i1,i2,i3)*rsxyst2(i1,i2,i3,m,n)+rsxyt2(i1,i2,i3,2,1)*rsxyt2(i1,i2,i3,m,n) + ty(i1,i2,i3)*rsxytt2(i1,i2,i3,m,n)
     
      rsxyzr23(i1,i2,i3,m,n)= rsxyr2(i1,i2,i3,0,2)*rsxyr2(i1,i2,i3,m,n) + rz(i1,i2,i3)*rsxyrr2(i1,i2,i3,m,n)+rsxyr2(i1,i2,i3,1,2)*rsxys2(i1,i2,i3,m,n) + sz(i1,i2,i3)*rsxyrs2(i1,i2,i3,m,n)+rsxyr2(i1,i2,i3,2,2)*rsxyt2(i1,i2,i3,m,n) + tz(i1,i2,i3)*rsxyrt2(i1,i2,i3,m,n)
     
      rsxyzs23(i1,i2,i3,m,n)= rsxys2(i1,i2,i3,0,2)*rsxyr2(i1,i2,i3,m,n) + rz(i1,i2,i3)*rsxyrs2(i1,i2,i3,m,n)+rsxys2(i1,i2,i3,1,2)*rsxys2(i1,i2,i3,m,n) + sz(i1,i2,i3)*rsxyss2(i1,i2,i3,m,n)+rsxys2(i1,i2,i3,2,2)*rsxyt2(i1,i2,i3,m,n) + tz(i1,i2,i3)*rsxyst2(i1,i2,i3,m,n)
     
      rsxyzt23(i1,i2,i3,m,n)= rsxyt2(i1,i2,i3,0,2)*rsxyr2(i1,i2,i3,m,n) + rz(i1,i2,i3)*rsxyrt2(i1,i2,i3,m,n)+rsxyt2(i1,i2,i3,1,2)*rsxys2(i1,i2,i3,m,n) + sz(i1,i2,i3)*rsxyst2(i1,i2,i3,m,n)+rsxyt2(i1,i2,i3,2,2)*rsxyt2(i1,i2,i3,m,n) + tz(i1,i2,i3)*rsxytt2(i1,i2,i3,m,n)
     
      ! ---- 4th order ---
     
      rsxyr4(i1,i2,i3,m,n)=(8.*(rsxy(i1+1,i2,i3,m,n)-rsxy(i1-1,i2,i3,m,n))-(rsxy(i1+2,i2,i3,m,n)-rsxy(i1-2,i2,i3,m,n)))*d14(0)
      rsxys4(i1,i2,i3,m,n)=(8.*(rsxy(i1,i2+1,i3,m,n)-rsxy(i1,i2-1,i3,m,n))-(rsxy(i1,i2+2,i3,m,n)-rsxy(i1,i2-2,i3,m,n)))*d14(1)
      rsxyt4(i1,i2,i3,m,n)=(8.*(rsxy(i1,i2,i3+1,m,n)-rsxy(i1,i2,i3-1,m,n))-(rsxy(i1,i2,i3+2,m,n)-rsxy(i1,i2,i3-2,m,n)))*d14(2)
     
      rsxyrr4(i1,i2,i3,m,n)=(-30.*rsxy(i1,i2,i3,m,n)+16.*(rsxy(i1+1,i2,i3,m,n)+rsxy(i1-1,i2,i3,m,n))-(rsxy(i1+2,i2,i3,m,n)+rsxy(i1-2,i2,i3,m,n)) )*d24(0)
     
      rsxyss4(i1,i2,i3,m,n)=(-30.*rsxy(i1,i2,i3,m,n)+16.*(rsxy(i1,i2+1,i3,m,n)+rsxy(i1,i2-1,i3,m,n))-(rsxy(i1,i2+2,i3,m,n)+rsxy(i1,i2-2,i3,m,n)) )*d24(1)
     
      rsxytt4(i1,i2,i3,m,n)=(-30.*rsxy(i1,i2,i3,m,n)+16.*(rsxy(i1,i2,i3+1,m,n)+rsxy(i1,i2,i3-1,m,n))-(rsxy(i1,i2,i3+2,m,n)+rsxy(i1,i2,i3-2,m,n)) )*d24(2)
     
      rsxyrs4(i1,i2,i3,m,n)=(8.*(rsxyr4(i1,i2+1,i3,m,n)-rsxyr4(i1,i2-1,i3,m,n))-(rsxyr4(i1,i2+2,i3,m,n)-rsxyr4(i1,i2-2,i3,m,n)))*d14(1)
     
      rsxyrt4(i1,i2,i3,m,n)=(8.*(rsxyr4(i1,i2,i3+1,m,n)-rsxyr4(i1,i2,i3-1,m,n))-(rsxyr4(i1,i2,i3+2,m,n)-rsxyr4(i1,i2,i3-2,m,n)))*d14(2)
     
      rsxyst4(i1,i2,i3,m,n)=(8.*(rsxys4(i1,i2,i3+1,m,n)-rsxys4(i1,i2,i3-1,m,n))-(rsxys4(i1,i2,i3+2,m,n)-rsxys4(i1,i2,i3-2,m,n)))*d14(2)
     
      rsxyx42(i1,i2,i3,m,n)= rx(i1,i2,i3)*rsxyr4(i1,i2,i3,m,n)+sx(i1,i2,i3)*rsxys4(i1,i2,i3,m,n)
      rsxyy42(i1,i2,i3,m,n)= ry(i1,i2,i3)*rsxyr4(i1,i2,i3,m,n)+sy(i1,i2,i3)*rsxys4(i1,i2,i3,m,n)
     
      rsxyxr42(i1,i2,i3,m,n)= rsxyr4(i1,i2,i3,0,0)*rsxyr4(i1,i2,i3,m,n) + rx(i1,i2,i3)*rsxyrr4(i1,i2,i3,m,n)+rsxyr4(i1,i2,i3,1,0)*rsxys4(i1,i2,i3,m,n) + sx(i1,i2,i3)*rsxyrs4(i1,i2,i3,m,n)
      rsxyxs42(i1,i2,i3,m,n)= rsxys4(i1,i2,i3,0,0)*rsxyr4(i1,i2,i3,m,n) + rx(i1,i2,i3)*rsxyrs4(i1,i2,i3,m,n)+rsxys4(i1,i2,i3,1,0)*rsxys4(i1,i2,i3,m,n) + sx(i1,i2,i3)*rsxyss4(i1,i2,i3,m,n)
     
      rsxyyr42(i1,i2,i3,m,n)= rsxyr4(i1,i2,i3,0,1)*rsxyr4(i1,i2,i3,m,n) + ry(i1,i2,i3)*rsxyrr4(i1,i2,i3,m,n)+rsxyr4(i1,i2,i3,1,1)*rsxys4(i1,i2,i3,m,n) + sy(i1,i2,i3)*rsxyrs4(i1,i2,i3,m,n)
      rsxyys42(i1,i2,i3,m,n)= rsxys4(i1,i2,i3,0,1)*rsxyr4(i1,i2,i3,m,n) + ry(i1,i2,i3)*rsxyrs4(i1,i2,i3,m,n)+rsxys4(i1,i2,i3,1,1)*rsxys4(i1,i2,i3,m,n) + sy(i1,i2,i3)*rsxyss4(i1,i2,i3,m,n)

      rsxyx43(i1,i2,i3,m,n)= rx(i1,i2,i3)*rsxyr4(i1,i2,i3,m,n)+sx(i1,i2,i3)*rsxys4(i1,i2,i3,m,n)+tx(i1,i2,i3)*rsxyt4(i1,i2,i3,m,n)
      rsxyy43(i1,i2,i3,m,n)= ry(i1,i2,i3)*rsxyr4(i1,i2,i3,m,n)+sy(i1,i2,i3)*rsxys4(i1,i2,i3,m,n)+ty(i1,i2,i3)*rsxyt4(i1,i2,i3,m,n)
      rsxyz43(i1,i2,i3,m,n)= rz(i1,i2,i3)*rsxyr4(i1,i2,i3,m,n)+sz(i1,i2,i3)*rsxys4(i1,i2,i3,m,n)+tz(i1,i2,i3)*rsxyt4(i1,i2,i3,m,n)
     
      rsxyxr43(i1,i2,i3,m,n)= rsxyr4(i1,i2,i3,0,0)*rsxyr4(i1,i2,i3,m,n) + rx(i1,i2,i3)*rsxyrr4(i1,i2,i3,m,n)+rsxyr4(i1,i2,i3,1,0)*rsxys4(i1,i2,i3,m,n) + sx(i1,i2,i3)*rsxyrs4(i1,i2,i3,m,n)+rsxyr4(i1,i2,i3,2,0)*rsxyt4(i1,i2,i3,m,n) + tx(i1,i2,i3)*rsxyrt4(i1,i2,i3,m,n)
     
      rsxyxs43(i1,i2,i3,m,n)= rsxys4(i1,i2,i3,0,0)*rsxyr4(i1,i2,i3,m,n) + rx(i1,i2,i3)*rsxyrs4(i1,i2,i3,m,n)+rsxys4(i1,i2,i3,1,0)*rsxys4(i1,i2,i3,m,n) + sx(i1,i2,i3)*rsxyss4(i1,i2,i3,m,n)+rsxys4(i1,i2,i3,2,0)*rsxyt4(i1,i2,i3,m,n) + tx(i1,i2,i3)*rsxyst4(i1,i2,i3,m,n)
     
      rsxyxt43(i1,i2,i3,m,n)= rsxyt4(i1,i2,i3,0,0)*rsxyr4(i1,i2,i3,m,n) + rx(i1,i2,i3)*rsxyrt4(i1,i2,i3,m,n)+rsxyt4(i1,i2,i3,1,0)*rsxys4(i1,i2,i3,m,n) + sx(i1,i2,i3)*rsxyst4(i1,i2,i3,m,n)+rsxyt4(i1,i2,i3,2,0)*rsxyt4(i1,i2,i3,m,n) + tx(i1,i2,i3)*rsxytt4(i1,i2,i3,m,n)
     
      rsxyyr43(i1,i2,i3,m,n)= rsxyr4(i1,i2,i3,0,1)*rsxyr4(i1,i2,i3,m,n) + ry(i1,i2,i3)*rsxyrr4(i1,i2,i3,m,n)+rsxyr4(i1,i2,i3,1,1)*rsxys4(i1,i2,i3,m,n) + sy(i1,i2,i3)*rsxyrs4(i1,i2,i3,m,n)+rsxyr4(i1,i2,i3,2,1)*rsxyt4(i1,i2,i3,m,n) + ty(i1,i2,i3)*rsxyrt4(i1,i2,i3,m,n)
     
      rsxyys43(i1,i2,i3,m,n)= rsxys4(i1,i2,i3,0,1)*rsxyr4(i1,i2,i3,m,n) + ry(i1,i2,i3)*rsxyrs4(i1,i2,i3,m,n)+rsxys4(i1,i2,i3,1,1)*rsxys4(i1,i2,i3,m,n) + sy(i1,i2,i3)*rsxyss4(i1,i2,i3,m,n)+rsxys4(i1,i2,i3,2,1)*rsxyt4(i1,i2,i3,m,n) + ty(i1,i2,i3)*rsxyst4(i1,i2,i3,m,n)
     
      rsxyyt43(i1,i2,i3,m,n)= rsxyt4(i1,i2,i3,0,1)*rsxyr4(i1,i2,i3,m,n) + ry(i1,i2,i3)*rsxyrt4(i1,i2,i3,m,n)+rsxyt4(i1,i2,i3,1,1)*rsxys4(i1,i2,i3,m,n) + sy(i1,i2,i3)*rsxyst4(i1,i2,i3,m,n)+rsxyt4(i1,i2,i3,2,1)*rsxyt4(i1,i2,i3,m,n) + ty(i1,i2,i3)*rsxytt4(i1,i2,i3,m,n)
     
      rsxyzr43(i1,i2,i3,m,n)= rsxyr4(i1,i2,i3,0,2)*rsxyr4(i1,i2,i3,m,n) + rz(i1,i2,i3)*rsxyrr4(i1,i2,i3,m,n)+rsxyr4(i1,i2,i3,1,2)*rsxys4(i1,i2,i3,m,n) + sz(i1,i2,i3)*rsxyrs4(i1,i2,i3,m,n)+rsxyr4(i1,i2,i3,2,2)*rsxyt4(i1,i2,i3,m,n) + tz(i1,i2,i3)*rsxyrt4(i1,i2,i3,m,n)
     
      rsxyzs43(i1,i2,i3,m,n)= rsxys4(i1,i2,i3,0,2)*rsxyr4(i1,i2,i3,m,n) + rz(i1,i2,i3)*rsxyrs4(i1,i2,i3,m,n)+rsxys4(i1,i2,i3,1,2)*rsxys4(i1,i2,i3,m,n) + sz(i1,i2,i3)*rsxyss4(i1,i2,i3,m,n)+rsxys4(i1,i2,i3,2,2)*rsxyt4(i1,i2,i3,m,n) + tz(i1,i2,i3)*rsxyst4(i1,i2,i3,m,n)
     
      rsxyzt43(i1,i2,i3,m,n)= rsxyt4(i1,i2,i3,0,2)*rsxyr4(i1,i2,i3,m,n) + rz(i1,i2,i3)*rsxyrt4(i1,i2,i3,m,n)+rsxyt4(i1,i2,i3,1,2)*rsxys4(i1,i2,i3,m,n) + sz(i1,i2,i3)*rsxyst4(i1,i2,i3,m,n)+rsxyt4(i1,i2,i3,2,2)*rsxyt4(i1,i2,i3,m,n) + tz(i1,i2,i3)*rsxytt4(i1,i2,i3,m,n)
     

!............... end statement functions

  ! write(*,'("ENTERING abcWave ")') 

  ierr=0

  uc                              = ipar( 0)
  numberOfComponents              = ipar( 1)
  grid                            = ipar( 2)
  gridType                        = ipar( 3)
  orderOfAccuracy                 = ipar( 4)
  gridIsImplicit                  = ipar( 5)
  twilightZone                    = ipar( 6)
  numberOfProcessors              = ipar( 7)
  debug                           = ipar( 8)
  myid                            = ipar( 9)
  assignKnownSolutionAtBoundaries = ipar(10)
  knownSolutionOption             = ipar(11)
  addForcingBC                    = ipar(12)
  forcingOption                   = ipar(13)
  useUpwindDissipation            = ipar(14)
  numGhost                        = ipar(15)  
  assignBCForImplicit  = ipar(16)
  bcApproach                      = ipar(17)
  numberOfFrequencies             = ipar(18)
 
  t         = rpar( 0)
  dt        = rpar( 1)
  dx(0)     = rpar( 2)
  dx(1)     = rpar( 3)
  dx(2)     = rpar( 4)
  dr(0)     = rpar( 5)
  dr(1)     = rpar( 6)
  dr(2)     = rpar( 7)
  ep        = rpar( 8) ! pointer for exact solution -- new : 110311 
  REAL_MIN  = rpar( 9)
  c         = rpar(10)
  cEM2      = rpar(11)

  if( t==0. )then
    ! write(*,'("abcWave: t=0: DO NOTHING since past time is not known ")') 
    return
  end if



  ierr=0
  axis1=0
  axis2=1
  axis3=2

  ! ** Fill in old parameters ***

  side=0; axis=0; 
  n1a=gridIndexRange(0,0); n1b=gridIndexRange(1,0);
  n2a=gridIndexRange(0,1); n2b=gridIndexRange(1,1);
  n3a=gridIndexRange(0,2); n3b=gridIndexRange(1,2);

  orderOfExtrapolation=3;     ! not used ? 
  useForcing = 0              ! not used ?
  ex=uc; ey=uc; ez=uc; 
  hx=uc; hy=uc; hz=uc;
  useWhereMask=1;

  if( twilightZone==1 )then
    forcingOption=twilightZoneForcing
  else
    forcingOption=0  ! CHECK ME
  end if

  method=nfdtd
  myid=0;

  eps=1; mu=1; 
  kx=1; ky=0; kz=0;
  slowStartInterval=.5; 



  tp=t-dt ! previous time
  tm=t-.5*dt ! midpoint in time

  ! -----------------------------------------
  ! ------- In 3D we just set hz=ez ---------
  ! -----------------------------------------
  if( nd.eq.3 )then
    hz=ez
  end if

  ! Return if there are no ABC's *wdh* Aug 16, 2018
  applyABC=.false.
  do axis=0,nd-1
  do side=0,1
    if( boundaryCondition(side,axis).eq.abcEM2 .or. boundaryCondition(side,axis).eq.absorbing )then
      applyABC=.true.
      exit
    end if
  end do
  end do
  if( .not.applyABC )then
    return
  end if


  if( twilightZone==1 .and. gridType==curvilinear )then
    write(*,'("abcWave: TZ NOT IMPLEMENTED FOR CURVILINEAR GRIDS")') 
    stop 1111
  end if

  if( debug.gt.1 .and. t.le.3*dt )then
    write(*,'("abcWave: grid=",i4," order=",i2," gridType=",i2," t=",e9.2", dt=",e9.2)') grid,orderOfAccuracy,gridType,t,dt
    write(*,'("abcWave: c=",e14.6," cEM2=",e14.6)') c,cEM2
    write(*,'("abcWave: useForcing=",i2," forcingOption=",i2)') useForcing,forcingOption
    write(*,'("abcWave: method=",i2," (5=nfdtd,7=bamx)")') method
    write(*,'("abcWave: gridIsImplicit=",i2," numGhost=",i4)') gridIsImplicit,numGhost
  end if
  if( method .ne. nfdtd .and. method.ne.bamx )then
     write(*,'("abcWave:ERROR: unknown method")')
     stop 12234
  end if 

  ! if( gridIsImplicit==1 )then
  !   write(*,'("abcWave: gridIsImplicit: SKIPPING explicit ABCs")')
  !   return
  ! end if


  ! if( debug.gt.1 .ab)then
  !   write(*,'(" abcMaxwell: **START** grid=",i4," side,axis=",2i2," projectDivLine=",i2)') grid,side,axis,projectDivLine
  ! end if
     


  ! for plane wave forcing 
  twoPi=8.*atan2(1.,1.)
  cc= c*sqrt( kx*kx+ky*ky+kz*kz )

  epsX=1.e-30 ! fix this ***

  ! Engquist-Majda 2nd-order
  !    u.xt = (1/c)*u.tt - c/2 * (u.yy + u.zz)   at x=0
  !         = c*( u.xx + .5*( u.yy + u.zz ) 
 
  ! We need un : u(t+dt) 
  !         u  : u(t)

  ! Generalized form:
  ! u.xt = c1abcem2*u.xx + c2abcem2*( u.yy + u.zz )
  !   Taylor: p0=1 p2=-1/2
  !   Cheby:  p0=1.00023, p2=-.515555
  p0=1.  
  p2=-.5
  ! p0=1.00023   !   Cheby on a subinterval
  ! p2=-.515555  !   Cheby on a subinterval
  c1abcem2=cEM2*p0
  c2abcem2=cEM2*(p0+p2)

  if( method.eq.bamx )then
    ! bamx: turn off tangential derivatives for now 
    p2=-p0             
    c2abcem2=cEM2*(p0+p2) 
  end if
     
  ! ***************************************************
  ! write(*,'("abcWave -- stop here for now ")') 
  ! stop 1234

  extra=-1  ! no need to do corners -- these are already done in another way
  extra=0   ! re-compute corners
  numberOfGhostPoints=orderOfAccuracy/2

  if( gridType.eq.curvilinear )then
    ! do this for testing:
    dx(0)=dr(0)
    dx(1)=dr(1)
    dx(2)=dr(2)
  end if

  dxa=dx(0) 
  dya=dx(1) 
  dza=dx(2) 

  if( gridType.eq.curvilinear )then
    ! On a curvilinear grid we need to make sure that there are valid values on the
    ! first ghost line -- these may be used on non-orthogonal grids (cross terms in uxx, uyy, ...)

    ! Note: This next loop only applies to boundaries that are ABCs : 
     extra1a=extra
     extra1b=extra
     extra2a=extra
     extra2b=extra
     if( nd.eq.3 )then
       extra3a=extra
       extra3b=extra
     else
       extra3a=0
       extra3b=0
     end if
     if( boundaryCondition(0,0).lt.0 )then
       extra1a=max(0,extra1a) ! over-ride extra=-1 : assign ends in periodic directions
       extra1b=extra1a
     else
       if( boundaryCondition(0,0).eq.0 )then
         extra1a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
       end if
       if( boundaryCondition(1,0).eq.0 )then
         extra1b=numberOfGhostPoints
       end if
     end if
     if( boundaryCondition(0,1).lt.0 )then
      extra2a=max(0,extra2a) ! over-ride extra=-1 : assign ends in periodic directions
      extra2b=extra2a
     else 
       if( boundaryCondition(0,1).eq.0 )then
         extra2a=numberOfGhostPoints
       end if
       if( boundaryCondition(1,1).eq.0 )then
         extra2b=numberOfGhostPoints
       end if
     end if
     if(  nd.eq.3 .and. boundaryCondition(0,2).lt.0 )then
      extra3a=max(0,extra3a) ! over-ride extra=-1 : assign ends in periodic directions
      extra3b=extra3a
     else 
       if( boundaryCondition(0,2).eq.0 )then
         extra3a=numberOfGhostPoints
       end if
       if( boundaryCondition(1,2).eq.0 )then
         extra3b=numberOfGhostPoints
       end if
     end if
     do axis=0,nd-1
     do side=0,1
       if( boundaryCondition(side,axis).eq.abcEM2 .or. boundaryCondition(side,axis).eq.absorbing )then
         ! write(*,'(" abc: grid,side,axis,bc=",3i2)') grid,side,axis,boundaryCondition(side,axis)
         n1a=gridIndexRange(0,0)-extra1a
         n1b=gridIndexRange(1,0)+extra1b
         n2a=gridIndexRange(0,1)-extra2a
         n2b=gridIndexRange(1,1)+extra2b
         n3a=gridIndexRange(0,2)-extra3a
         n3b=gridIndexRange(1,2)+extra3b
         if( axis.eq.0 )then
           n1a=gridIndexRange(side,axis)
           n1b=gridIndexRange(side,axis)
         else if( axis.eq.1 )then
           n2a=gridIndexRange(side,axis)
           n2b=gridIndexRange(side,axis)
         else
           n3a=gridIndexRange(side,axis)
           n3b=gridIndexRange(side,axis)
         end if
         is1=0
         is2=0
         is3=0
         if( axis.eq.0 )then
           is1=1-2*side
         else if( axis.eq.1 )then
           is2=1-2*side
         else if( axis.eq.2 )then
           is3=1-2*side
         else
           stop 5
         end if
         is = 1 - 2*side
         axisp1=mod(axis+1,nd)
         axisp2=mod(axis+2,nd)
         ! (js1,js2,js3) used to compute tangential derivatives
         js1=0
         js2=0
         js3=0
         if( axisp1.eq.0 )then
           js1=1-2*side
         else if( axisp1.eq.1 )then
           js2=1-2*side
         else if( axisp1.eq.2 )then
           js3=1-2*side
         else
           stop 5
         end if
         ! (ks1,ks2,ks3) used to compute second tangential derivative
         ks1=0
         ks2=0
         ks3=0
         if( axisp2.eq.0 )then
           ks1=1-2*side
         else if( axisp2.eq.1 )then
           ks2=1-2*side
         else if( axisp2.eq.2 )then
           ks3=1-2*side
         else
           stop 5
         end if
      if( orderOfAccuracy.eq.2 )then

       ! ** we could also impose the first-order in time explicit formula **
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b
          un(i1-is1,i2-is2,i3-is3,ex)=2.*un(i1,i2,i3,ex)-un(i1+is1,i2+is2,i3+is3,ex)
       end do
       end do
       end do

      else if( orderOfAccuracy.eq.4 )then

       ! extrap to order 3 in case we adjust for incident fields 
       do i3=n3a,n3b
       do i2=n2a,n2b
       do i1=n1a,n1b
          un(i1-is1,i2-is2,i3-is3,ex)=3.*un(i1,i2,i3,ex)-3.*un(i1+is1,i2+is2,i3+is3,ex)+un(i1+2*is1,i2+2*is2,i3+2*is3,ex)
       end do
       end do
       end do

      else
        stop 8822 ! unknown orderOfAccuracy
      end if

       end if
     end do
     end do

  end if

  ! -- initialize for forcing:
  z=0.
  do mm=0,11
    forcex(mm)=0.
    forcey(mm)=0.
    forcez(mm)=0.
    forcep(mm)=0.
    forcef(mm)=0.
  end do

  ! ------------------------------------------------------------------------
  ! ------------------ EDGES AND CORNERS------------------------------------
  ! ------------------------------------------------------------------------
  if( orderOfAccuracy==2 )then
    
    ! write(*,*) 'abcWave: assign edges and corners, order 2'
      ! ------------------------------------------------------------------------
      ! ------------------Corners-----------------------------------------------
      ! ------------------------------------------------------------------------
      ! We need to assign points in the corner region:
      !
      !           |  |  |
      !           +--+--+--
      !           |  |  |
      !     D--A--X--+--+--
      !     |  |  |
      !     D--C--B
      !     |  |  |
      !     D--D--D
      ! 
      if( nd.eq.2 )then
       ! **** 2D ****
        i3=gridIndexRange(0,2)
        do side1=0,1
        do side2=0,1
         bc1=boundaryCondition(side1,0)
         bc2=boundaryCondition(side2,1)
         if( ((bc1.eq.abcEM2 .or. bc1.eq.absorbing) .and. bc2.gt.0 ) .or. ((bc2.eq.abcEM2 .or. bc2.eq.absorbing) .and. bc1.gt.0 ) )then
           ! --- One of the faces at this corner is an ABC and the other has bc>0 ---         
           ! Adjacent side is also an ABC: 
           adjacentFaceIsABC = bc1.eq.abcEM2 .or. bc1.eq.absorbing .and. bc2.eq.abcEM2 .or. bc2.eq.absorbing
           i1=gridIndexRange(side1,0) ! (i1,i2,i3)=corner point
           i2=gridIndexRange(side2,1)
           if( mask(i1,i2,i3).gt.0 )then ! *wdh* 090712
            ! write(*,'(" ABC:set corner: grid,side1,side2,i1,i2=",3i3,2i5)') grid,side1,side2,i1,i2
            ! --- start by extrapolating all points on the extended boundary and adjacent to the corner ---
            !* is1=1-2*side1
            !* is2=1-2*side2
            !* is3=0
            !* j3=0
            !* extrapolateGhost(ex,ey,hz,numberOfGhostPoints,numberOfGhostPoints,0)
            ! --- Assign points 'A' and 'B" on the extended boundary ---
            if( adjacentFaceIsABC )then
              is1=1-2*side1
              is2=0
              is3=0
              if( gridType.eq.rectangular )then
               if( .true. )then
                ! *new* way
                is1=1-2*side1
                is2=1-2*side2
                  ! getForcingEM2(X,2,tm,is1,is2,forcex)
                   if( forcingOption.eq.twilightZoneForcing )then
                     ! Test: set to exact solution at time t:
                     ! x=xy(i1-is1,i2,i3,0)
                     ! y=xy(i1-is1,i2,i3,1)
                     ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                     ! un(i1-is1,i2,i3,ey)=eyTrue
                     ! add TZ forcing *wdh* Sept 17, 2016
                     ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                     x=xy(i1,i2,i3,0)
                     y=xy(i1,i2,i3,1)
                        ! Values for forcep(ex) are currently needed at corners:
                          call ogDeriv(ep, 1,1,0,0,x,y,z,tp,ex,utx)
                          call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                          call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                        forcep(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                        ! OGDERIV(1,1,0,0,x,y,z,tp,ey,utx)
                        ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                        ! forcep(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                        ! OGDERIV(1,1,0,0,x,y,z,tp,hz,utx)
                        ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
                        ! forcep(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
                     ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                   end if
                   if( forcingOption.eq.twilightZoneForcing )then
                     ! Test: set to exact solution at time t:
                     ! x=xy(i1-is1,i2,i3,0)
                     ! y=xy(i1-is1,i2,i3,1)
                     ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                     ! un(i1-is1,i2,i3,ey)=eyTrue
                     ! add TZ forcing *wdh* Sept 17, 2016
                     ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                     x=xy(i1,i2,i3,0)
                     y=xy(i1,i2,i3,1)
                        ! Values for forcef(ex) are currently needed at corners:
                          call ogDeriv(ep, 1,1,0,0,x,y,z,t,ex,utx)
                          call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                          call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                        forcef(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                        ! OGDERIV(1,1,0,0,x,y,z,t,ey,utx)
                        ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
                        ! forcef(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                        ! OGDERIV(1,1,0,0,x,y,z,t,hz,utx)
                        ! OGDERIV(0,2,0,0,x,y,z,t,hz,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,t,hz,uyy)
                        ! forcef(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
                     ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                   end if
                  ! do idir=0,2
                  do idir=0,0 ! ex only 
                    forcex(idir)=.5*(forcep(idir)+forcef(idir))
                  end do
                  ! getForcingEM2(Y,2,tm,is1,is2,forcey)
                   if( forcingOption.eq.twilightZoneForcing )then
                     ! Test: set to exact solution at time t:
                     ! x=xy(i1-is1,i2,i3,0)
                     ! y=xy(i1-is1,i2,i3,1)
                     ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                     ! un(i1-is1,i2,i3,ey)=eyTrue
                     ! add TZ forcing *wdh* Sept 17, 2016
                     ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                     x=xy(i1,i2,i3,0)
                     y=xy(i1,i2,i3,1)
                          call ogDeriv(ep, 1,0,1,0,x,y,z,tp,ex,uty)
                          call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                          call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                        forcep(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                        ! ! Values for forcep(ey) are currently needed at corners:
                        ! OGDERIV(1,0,1,0,x,y,z,tp,ey,uty)
                        ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                        ! forcep(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                        ! OGDERIV(1,0,1,0,x,y,z,tp,hz,uty)
                        ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
                        ! forcep(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
                     ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                   end if
                   if( forcingOption.eq.twilightZoneForcing )then
                     ! Test: set to exact solution at time t:
                     ! x=xy(i1-is1,i2,i3,0)
                     ! y=xy(i1-is1,i2,i3,1)
                     ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                     ! un(i1-is1,i2,i3,ey)=eyTrue
                     ! add TZ forcing *wdh* Sept 17, 2016
                     ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                     x=xy(i1,i2,i3,0)
                     y=xy(i1,i2,i3,1)
                          call ogDeriv(ep, 1,0,1,0,x,y,z,t,ex,uty)
                          call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                          call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                        forcef(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                        ! ! Values for forcef(ey) are currently needed at corners:
                        ! OGDERIV(1,0,1,0,x,y,z,t,ey,uty)
                        ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
                        ! forcef(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                        ! OGDERIV(1,0,1,0,x,y,z,t,hz,uty)
                        ! OGDERIV(0,2,0,0,x,y,z,t,hz,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,t,hz,uyy)
                        ! forcef(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
                     ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                   end if
                  ! do idir=0,2
                  do idir=0,0 ! ex only 
                    forcey(idir)=.5*(forcep(idir)+forcef(idir))
                  end do
                  ! At a corner there are two coupled equations we need to solve for ghost points A,B below
                  !                   |
                  !                 A +---+----
                  !                   B
                  !  f(u)  = [ f(u_old) - A (u_old) ] + A u = 0 
                  !  [ a11 a12 ][ uA ] = [ a11 a12 ][ uA_old ] - [ f1(u_old) ]
                  !  [ a21 a22 ][ uB ]   [ a21 a22 ][ uB_old ]   [ f2(u_old) ]
                  isign1=1-2*side1
                  isign2=1-2*side2
                  ! first evaluate residuals in equations given current (wrong) values at A, B
                  r1 = isign1*(unx22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ex))/(dt)- .5*( c1abcem2*unxx22r(i1,i2,i3,ex) + c2abcem2*unyy22r(i1,i2,i3,ex) +c1abcem2* uxx22r(i1,i2,i3,ex) + c2abcem2* uyy22r(i1,i2,i3,ex) ) - forcex(ex) 
                  r2 = isign2*(uny22r(i1,i2,i3,ex)-uy22r(i1,i2,i3,ex))/(dt)- .5*( c1abcem2*unyy22r(i1,i2,i3,ex) + c2abcem2*unxx22r(i1,i2,i3,ex) +c1abcem2* uyy22r(i1,i2,i3,ex) + c2abcem2* uxx22r(i1,i2,i3,ex) ) - forcey(ex)
                  a11 = -1./(2.*dt*dx(0))  - .5*c1abcem2/(dx(0)**2)
                  a12 = -.5*c2abcem2/(dx(1)**2)
                  a21 = -.5*c2abcem2/(dx(0)**2)
                  a22 = -1./(2.*dt*dx(1))  - .5*c1abcem2/(dx(1)**2)
                  det = a11*a22-a21*a12
                  uA = un(i1-isign1,i2,i3,ex)
                  uB = un(i1,i2-isign2,i3,ex)
                  f1 = a11*uA + a12*uB - r1 
                  f2 = a21*uA + a22*uB - r2 
                  ! Solve for A, B
                  un(i1-isign1,i2,i3,ex) = ( f1*a22 - f2*a12)/det
                  un(i1,i2-isign2,i3,ex) = (-f1*a21 + f2*a11)/det
               else           
                ! *old* way          
                ! --- Assign point 'A' on the extended boundary ---
                ! Use first-order-in-time formula since it doesn't require other ghost point 'B' at new time (un)
                 if( forcingOption.eq.twilightZoneForcing )then
                   ! Test: set to exact solution at time t:
                   ! x=xy(i1-is1,i2,i3,0)
                   ! y=xy(i1-is1,i2,i3,1)
                   ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                   ! un(i1-is1,i2,i3,ey)=eyTrue
                   ! add TZ forcing *wdh* Sept 17, 2016
                   ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                      ! Values for forcex(ex) are currently needed at corners:
                        call ogDeriv(ep, 1,1,0,0,x,y,z,tp,ex,utx)
                        call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                        call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                      forcex(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                      ! OGDERIV(1,1,0,0,x,y,z,tp,ey,utx)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                      ! forcex(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                      ! OGDERIV(1,1,0,0,x,y,z,tp,hz,utx)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
                      ! forcex(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
                   ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                 end if
                un(i1-is1,i2-is2,i3-is3,ex)=(un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(2.*dxa*dt)*(c1abcem2*uxx22r(i1,i2,i3,ex)+c2abcem2*uyy22r(i1,i2,i3,ex)+forcex(ex)))
               end if
              else
               ! curvilinear grid 
               side=side1
               axis=0 
                 is =1-2*side
                 dr0=dr(axis)
                 rx0 = rsxy(i1,i2,i3,axis,0)
                 ry0 = rsxy(i1,i2,i3,axis,1)
                 rxNormSq = rx0**2 + ry0**2 
                 rxNorm = max( epsX, sqrt(rxNormSq) )
                 rxx0 = rsxyx22(i1,i2,i3,axis,0)
                 ryy0 = rsxyy22(i1,i2,i3,axis,1)
                 ! cm1 : coeff of u(i1-is1,i2-is2,i3-is3,cc) in g (given below): 
                 cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0)/(2.*dr0) ) )
                 ! u: derivatives at time t: 
                 ! un: derivatives at time t+dt : evaluate using the incorrect ghost values 
                 if( axis.eq.0 )then
                   ur0   =   ur2(i1,i2,i3,ex)
                   urr0  =  urr2(i1,i2,i3,ex)
                   unr0  =  unr2(i1,i2,i3,ex)
                   unrr0 = unrr2(i1,i2,i3,ex)
                 else
                   ur0   =   us2(i1,i2,i3,ex)
                   urr0  =  uss2(i1,i2,i3,ex)
                   unr0  =  uns2(i1,i2,i3,ex)
                   unrr0 = unss2(i1,i2,i3,ex)
                 end if
                 uxx0 = uxx22(i1,i2,i3,ex)
                 uyy0 = uyy22(i1,i2,i3,ex)
                 unxx0 = unxx22(i1,i2,i3,ex)
                 unyy0 = unyy22(i1,i2,i3,ex)
                 ! first evaluate the BC using the incorrect ghost values 
                 Dn2 = rxNormSq*(unrr0+urr0) 
                 Lu  = unxx0+unyy0 + uxx0+uyy0 
                 g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
                 ! note: this assumes an orthogonal grid -- we should make sure that the 
                 !       ghost values have an initial guess in them (extrapolate ?)
                 un(i1-is1,i2-is2,i3-is3,ex) = -(g - cm1*un(i1-is1,i2-is2,i3-is3,ex) )/cm1 
              end if
              ! --- Assign point 'B'  on the extended boundary --
              is1=0
              is2=1-2*side2
              if( gridType.eq.rectangular )then
               if( .false. )then
                ! *old* way
                 if( forcingOption.eq.twilightZoneForcing )then
                   ! Test: set to exact solution at time t:
                   ! x=xy(i1-is1,i2,i3,0)
                   ! y=xy(i1-is1,i2,i3,1)
                   ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                   ! un(i1-is1,i2,i3,ey)=eyTrue
                   ! add TZ forcing *wdh* Sept 17, 2016
                   ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                        call ogDeriv(ep, 1,0,1,0,x,y,z,tp,ex,uty)
                        call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                        call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                      forcey(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                      ! ! Values for forcey(ey) are currently needed at corners:
                      ! OGDERIV(1,0,1,0,x,y,z,tp,ey,uty)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                      ! forcey(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                      ! OGDERIV(1,0,1,0,x,y,z,tp,hz,uty)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
                      ! forcey(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
                   ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                 end if
                un(i1-is1,i2-is2,i3-is3,ex)=(un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(2.*dya*dt)*(c1abcem2*uyy22r(i1,i2,i3,ex)+c2abcem2*uxx22r(i1,i2,i3,ex)+forcey(ex)))
               end if
              else
               ! curvilinear grid 
               side=side2
               axis=1 
                 is =1-2*side
                 dr0=dr(axis)
                 rx0 = rsxy(i1,i2,i3,axis,0)
                 ry0 = rsxy(i1,i2,i3,axis,1)
                 rxNormSq = rx0**2 + ry0**2 
                 rxNorm = max( epsX, sqrt(rxNormSq) )
                 rxx0 = rsxyx22(i1,i2,i3,axis,0)
                 ryy0 = rsxyy22(i1,i2,i3,axis,1)
                 ! cm1 : coeff of u(i1-is1,i2-is2,i3-is3,cc) in g (given below): 
                 cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0)/(2.*dr0) ) )
                 ! u: derivatives at time t: 
                 ! un: derivatives at time t+dt : evaluate using the incorrect ghost values 
                 if( axis.eq.0 )then
                   ur0   =   ur2(i1,i2,i3,ex)
                   urr0  =  urr2(i1,i2,i3,ex)
                   unr0  =  unr2(i1,i2,i3,ex)
                   unrr0 = unrr2(i1,i2,i3,ex)
                 else
                   ur0   =   us2(i1,i2,i3,ex)
                   urr0  =  uss2(i1,i2,i3,ex)
                   unr0  =  uns2(i1,i2,i3,ex)
                   unrr0 = unss2(i1,i2,i3,ex)
                 end if
                 uxx0 = uxx22(i1,i2,i3,ex)
                 uyy0 = uyy22(i1,i2,i3,ex)
                 unxx0 = unxx22(i1,i2,i3,ex)
                 unyy0 = unyy22(i1,i2,i3,ex)
                 ! first evaluate the BC using the incorrect ghost values 
                 Dn2 = rxNormSq*(unrr0+urr0) 
                 Lu  = unxx0+unyy0 + uxx0+uyy0 
                 g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
                 ! note: this assumes an orthogonal grid -- we should make sure that the 
                 !       ghost values have an initial guess in them (extrapolate ?)
                 un(i1-is1,i2-is2,i3-is3,ex) = -(g - cm1*un(i1-is1,i2-is2,i3-is3,ex) )/cm1 
              end if
            else if( .true. )then
              ! .false. .and. bc1.ne.symmetryBoundaryCondition .and. bc2.ne.symmetryBoundaryCondition )then
              ! --- adjacent face is NOT another ABC ---
              ! Do this for now *wdh* Sept 19, 2016
              if( orderOfAccuracy.eq.2 )then
                extrapOrder=3
              else if( orderOfAccuracy.eq.4 )then
                extrapOrder=5
              else
                stop 4114
              end if
              if( bc1.ne.symmetryBoundaryCondition )then
                  ksv(0)=0
                  ksv(1)=0
                  ksv(2)=0
                  ksv(0)=1-2*side1
                  ks1=ksv(0)
                  ks2=ksv(1)
                  ks3=ksv(2)
                  if( extrapOrder.eq.3 )then
                    un(i1-ks1,i2-ks2,i3-ks3,ex)=(3.*un(i1,i2,i3,ex)-3.*un(i1+ks1,i2+ks2,i3+ks3,ex)+un(i1+2*ks1,i2+2*ks2,i3+2*ks3,ex))
                    ! un(i1-ks1,i2-ks2,i3-ks3,ey)=extrap3(un,i1,i2,i3,ey,ks1,ks2,ks3)
                    ! un(i1-ks1,i2-ks2,i3-ks3,hz)=extrap3(un,i1,i2,i3,hz,ks1,ks2,ks3)
                  else if( extrapOrder.eq.5 )then
                    un(i1-ks1,i2-ks2,i3-ks3,ex)=(5.*un(i1,i2,i3,ex)-10.*un(i1+ks1,i2+ks2,i3+ks3,ex)+10.*un(i1+2*ks1,i2+2*ks2,i3+2*ks3,ex)-5.*un(i1+3*ks1,i2+3*ks2,i3+3*ks3,ex)+un(i1+4*ks1,i2+4*ks2,i3+4*ks3,ex))
                    ! un(i1-ks1,i2-ks2,i3-ks3,ey)=extrap5(un,i1,i2,i3,ey,ks1,ks2,ks3)
                    ! un(i1-ks1,i2-ks2,i3-ks3,hz)=extrap5(un,i1,i2,i3,hz,ks1,ks2,ks3)
                  else
                    stop 1782
                  end if
              end if
              if( bc2.ne.symmetryBoundaryCondition )then
                  ksv(0)=0
                  ksv(1)=0
                  ksv(2)=0
                  ksv(1)=1-2*side2
                  ks1=ksv(0)
                  ks2=ksv(1)
                  ks3=ksv(2)
                  if( extrapOrder.eq.3 )then
                    un(i1-ks1,i2-ks2,i3-ks3,ex)=(3.*un(i1,i2,i3,ex)-3.*un(i1+ks1,i2+ks2,i3+ks3,ex)+un(i1+2*ks1,i2+2*ks2,i3+2*ks3,ex))
                    ! un(i1-ks1,i2-ks2,i3-ks3,ey)=extrap3(un,i1,i2,i3,ey,ks1,ks2,ks3)
                    ! un(i1-ks1,i2-ks2,i3-ks3,hz)=extrap3(un,i1,i2,i3,hz,ks1,ks2,ks3)
                  else if( extrapOrder.eq.5 )then
                    un(i1-ks1,i2-ks2,i3-ks3,ex)=(5.*un(i1,i2,i3,ex)-10.*un(i1+ks1,i2+ks2,i3+ks3,ex)+10.*un(i1+2*ks1,i2+2*ks2,i3+2*ks3,ex)-5.*un(i1+3*ks1,i2+3*ks2,i3+3*ks3,ex)+un(i1+4*ks1,i2+4*ks2,i3+4*ks3,ex))
                    ! un(i1-ks1,i2-ks2,i3-ks3,ey)=extrap5(un,i1,i2,i3,ey,ks1,ks2,ks3)
                    ! un(i1-ks1,i2-ks2,i3-ks3,hz)=extrap5(un,i1,i2,i3,hz,ks1,ks2,ks3)
                  else
                    stop 1782
                  end if
              end if 
            end if ! end if adjacentFaceIsABC
            ! --- Now extrapolate all other points on the extended boundary and adjacent to the corner ---
            is1=1-2*side1
            is2=1-2*side2
            is3=0
            j3=0
             do m3=0,0
             do m2=0,numberOfGhostPoints
             do m1=0,numberOfGhostPoints
              mSum = m1+m2+m3
              if( mSum.gt.0 .and. mSum.ne.1 )then ! mSum=1 : these points have already been assigned
               ! extrap ghost point (j1,j2,j3)
               j1=i1-is1*m1
               j2=i2-is2*m2
               j3=i3-is3*m3
               ! js1=0 if m1=0 and js1=is1 if m1>0
               js1 = is1*min(m1,1)
               js2 = is2*min(m2,1)
               js3 = is3*min(m3,1)
               if( orderOfAccuracy.eq.2 )then
                 ! Changed to third-order extra *wdh* Sept 18, 2016
                 un(j1,j2,j3,ex)=(3.*un(j1+js1,j2+js2,j3+js3,ex)-3.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex))
                 ! un(j1,j2,j3,ey)=extrap3(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                 ! un(j1,j2,j3,hz)=extrap3(un,j1+js1,j2+js2,j3+js3,hz,js1,js2,js3)
                 ! un(j1,j2,j3,ex)=extrap2(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                 ! un(j1,j2,j3,ey)=extrap2(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                 ! un(j1,j2,j3,hz)=extrap2(un,j1+js1,j2+js2,j3+js3,hz,js1,js2,js3)
               else                                                                           
                 ! Note: adjust for incident fields should take into account the width of extrapolation: 
                 if( m1.le.1 .and. m2.le.1 .and. m3.le.1 )then
                   ! increased extrapolation to order=5 *wdh* June 20, 2016
                   un(j1,j2,j3,ex)=(5.*un(j1+js1,j2+js2,j3+js3,ex)-10.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+10.*un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex)-5.*un(j1+js1+3*js1,j2+js2+3*js2,j3+js3+3*js3,ex)+un(j1+js1+4*js1,j2+js2+4*js2,j3+js3+4*js3,ex))
                   ! un(j1,j2,j3,ey)=extrap5(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                   ! un(j1,j2,j3,hz)=extrap5(un,j1+js1,j2+js2,j3+js3,hz,js1,js2,js3)
                   !un(j1,j2,j3,ex)=extrap3(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                   !un(j1,j2,j3,ey)=extrap3(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                   !un(j1,j2,j3,hz)=extrap3(un,j1+js1,j2+js2,j3+js3,hz,js1,js2,js3)
                 else
                   ! 2nd-ghost line 
                   un(j1,j2,j3,ex)=(5.*un(j1+js1,j2+js2,j3+js3,ex)-10.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+10.*un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex)-5.*un(j1+js1+3*js1,j2+js2+3*js2,j3+js3+3*js3,ex)+un(j1+js1+4*js1,j2+js2+4*js2,j3+js3+4*js3,ex))
                   ! un(j1,j2,j3,ey)=extrap5(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                   ! un(j1,j2,j3,hz)=extrap5(un,j1+js1,j2+js2,j3+js3,hz,js1,js2,js3)
                   ! un(j1,j2,j3,ex)=extrap4(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                   ! un(j1,j2,j3,ey)=extrap4(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                   ! un(j1,j2,j3,hz)=extrap4(un,j1+js1,j2+js2,j3+js3,hz,js1,js2,js3)
                 end if
              end if
              end if
             end do
             end do
             end do
           end if ! mask 
          end if ! if one face is ABC
        end do          
        end do          
      else 
        ! ***** 3D *****
         do edgeDirection=0,2 ! direction parallel to the edge
         do sidea=0,1
         do sideb=0,1
          if( edgeDirection.eq.0 )then
            side1=0
            side2=sidea
            side3=sideb
          else if( edgeDirection.eq.1 )then
            side1=sideb 
            side2=0
            side3=sidea
          else
            side1=sidea
            side2=sideb
            side3=0
          end if
         is1=1-2*(side1)
         is2=1-2*(side2)
         is3=1-2*(side3)
         if( edgeDirection.eq.2 )then
          is3=0
          n1a=gridIndexRange(side1,0)
          n1b=gridIndexRange(side1,0)
          n2a=gridIndexRange(side2,1)
          n2b=gridIndexRange(side2,1)
          n3a=gridIndexRange(0,2)
          n3b=gridIndexRange(1,2)
          bc1=boundaryCondition(side1,0)
          bc2=boundaryCondition(side2,1)
         else if( edgeDirection.eq.1 )then
          is2=0
          n1a=gridIndexRange(side1,0)
          n1b=gridIndexRange(side1,0)
          n2a=gridIndexRange(    0,1)
          n2b=gridIndexRange(    1,1)
          n3a=gridIndexRange(side3,2)
          n3b=gridIndexRange(side3,2)
          bc1=boundaryCondition(side1,0)
          bc2=boundaryCondition(side3,2)
         else 
          is1=0  
          n1a=gridIndexRange(    0,0)
          n1b=gridIndexRange(    1,0)
          n2a=gridIndexRange(side2,1)
          n2b=gridIndexRange(side2,1)
          n3a=gridIndexRange(side3,2)
          n3b=gridIndexRange(side3,2)
          bc1=boundaryCondition(side2,1)
          bc2=boundaryCondition(side3,2)
         end if
         if( ( (bc1.eq.abcEM2 .or. bc1.eq.absorbing) .and. bc2.gt.0 ) .or. ( (bc2.eq.abcEM2 .or. bc2.eq.absorbing) .and. bc1.gt.0 ) )then
          ! --- One face is an ABC  and the other has bc>0 ---
          ! Adjacent side is also an ABC: 
          adjacentFaceIsABC = bc1.eq.abcEM2 .or. bc1.eq.absorbing .and. bc2.eq.abcEM2 .or. bc2.eq.absorbing
          if( edgeDirection.eq.0 )then
           i2=n2a
           i3=n3a
           do i1=n1a,n1b
           if( mask(i1,i2,i3).gt.0 )then ! *wdh* 090712
             ! Use first-order-in-time formula since it doesn't require other ghost point at new time (un)
             if( adjacentFaceIsABC )then
              if( gridType.eq.rectangular )then
                if( forcingOption.eq.twilightZoneForcing )then
                  ! Test: set to exact solution at time t:
                  ! x=xy(i1-is1,i2,i3,0)
                  ! y=xy(i1-is1,i2,i3,1)
                  ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                  ! un(i1-is1,i2,i3,ey)=eyTrue
                  ! add TZ forcing *wdh* Sept 17, 2016
                  ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                  x=xy(i1,i2,i3,0)
                  y=xy(i1,i2,i3,1)
                    ! ------ Cartesian Grid 3d forcing ----------
                    z=xy(i1,i2,i3,2)
                       call ogDeriv(ep, 1,0,0,1,x,y,z,tp,ex,utz)
                       call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                       call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                       call ogDeriv(ep, 0,0,0,2,x,y,z,tp,ex,uzz)
                     forcez(ex) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) ) 
                     ! OGDERIV(1,0,0,1,x,y,z,tp,ey,utz)
                     ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                     ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                     ! OGDERIV(0,0,0,2,x,y,z,tp,ey,uzz)
                     ! forcez(ey) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) ) 
                     ! OGDERIV(1,0,0,1,x,y,z,tp,ez,utz)
                     ! OGDERIV(0,2,0,0,x,y,z,tp,ez,uxx)
                     ! OGDERIV(0,0,2,0,x,y,z,tp,ez,uyy)
                     ! OGDERIV(0,0,0,2,x,y,z,tp,ez,uzz)
                     ! forcez(ez) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) )
                  ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                end if
               un(i1,i2,i3-is3,ex)=(un(i1+0,i2+0,i3+is3,ex)-(u(i1+0,i2+0,i3+is3,ex)-u(i1-0,i2-0,i3-is3,ex))-(2.*dza*dt)*(c1abcem2*uzz23r(i1,i2,i3,ex)+c2abcem2*(uxx23r(i1,i2,i3,ex)+uyy23r(i1,i2,i3,ex))+forcez(ex)))
                if( forcingOption.eq.twilightZoneForcing )then
                  ! Test: set to exact solution at time t:
                  ! x=xy(i1-is1,i2,i3,0)
                  ! y=xy(i1-is1,i2,i3,1)
                  ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                  ! un(i1-is1,i2,i3,ey)=eyTrue
                  ! add TZ forcing *wdh* Sept 17, 2016
                  ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                  x=xy(i1,i2,i3,0)
                  y=xy(i1,i2,i3,1)
                    ! ------ Cartesian Grid 3d forcing ----------
                    z=xy(i1,i2,i3,2)
                       call ogDeriv(ep, 1,0,1,0,x,y,z,tp,ex,uty)
                       call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                       call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                       call ogDeriv(ep, 0,0,0,2,x,y,z,tp,ex,uzz)
                     forcey(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) ) 
                     ! OGDERIV(1,0,1,0,x,y,z,tp,ey,uty)
                     ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                     ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                     ! OGDERIV(0,0,0,2,x,y,z,tp,ey,uzz)
                     ! forcey(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) ) 
                     ! OGDERIV(1,0,1,0,x,y,z,tp,ez,uty)
                     ! OGDERIV(0,2,0,0,x,y,z,tp,ez,uxx)
                     ! OGDERIV(0,0,2,0,x,y,z,tp,ez,uyy)
                     ! OGDERIV(0,0,0,2,x,y,z,tp,ez,uzz)
                     ! forcey(ez) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) )
                  ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                end if
               un(i1,i2-is2,i3,ex)=(un(i1+0,i2+is2,i3+0,ex)-(u(i1+0,i2+is2,i3+0,ex)-u(i1-0,i2-is2,i3-0,ex))-(2.*dya*dt)*(c1abcem2*uyy23r(i1,i2,i3,ex)+c2abcem2*(uxx23r(i1,i2,i3,ex)+uzz23r(i1,i2,i3,ex))+forcey(ex)))
              else
                 is =1-2*side3
                 dr0=dr(2)
                 rx0 = rsxy(i1,i2,i3,2,0)
                 ry0 = rsxy(i1,i2,i3,2,1)
                 rz0 = rsxy(i1,i2,i3,2,2)
                 rxNormSq = rx0**2 + ry0**2 + rz0**2
                 rxNorm = max( epsX, sqrt(rxNormSq) )
                 rxx0 = rsxyx23(i1,i2,i3,2,0)
                 ryy0 = rsxyy23(i1,i2,i3,2,1)
                 rzz0 = rsxyz23(i1,i2,i3,2,2)
                 ! cm1 : coeff of u(i1-0,i2-0,i3-is3,cc) in g (given below): 
                 cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0+rzz0)/(2.*dr0) ) )
                 ! derivatives at time t: 
                 ! derivatives at time t+dt : evaluate using the incorrect ghost values 
                 if( 2.eq.0 )then
                   ur0   =   ur2(i1,i2,i3,ex)
                   urr0  =  urr2(i1,i2,i3,ex)
                   unr0  =  unr2(i1,i2,i3,ex)
                   unrr0 = unrr2(i1,i2,i3,ex)
                 else if( 2.eq.1 )then
                   ur0   =   us2(i1,i2,i3,ex)
                   urr0  =  uss2(i1,i2,i3,ex)
                   unr0  =  uns2(i1,i2,i3,ex)
                   unrr0 = unss2(i1,i2,i3,ex)
                 else
                   ur0   =   ut2(i1,i2,i3,ex)
                   urr0  =  utt2(i1,i2,i3,ex)
                   unr0  =  unt2(i1,i2,i3,ex)
                   unrr0 = untt2(i1,i2,i3,ex)
                 end if
                 uxx0 = uxx23(i1,i2,i3,ex)
                 uyy0 = uyy23(i1,i2,i3,ex)
                 uzz0 = uzz23(i1,i2,i3,ex)
                 unxx0 = unxx23(i1,i2,i3,ex)
                 unyy0 = unyy23(i1,i2,i3,ex)
                 unzz0 = unzz23(i1,i2,i3,ex)
                 ! first evaluate the BC using the incorrect ghost values 
                 Dn2 = rxNormSq*(unrr0+urr0) 
                 Lu  = unxx0 + unyy0 + unzz0 + uxx0 + uyy0 + uzz0 
                 g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
                 ! note: this assumes an orthogonal grid -- we should make sure that the 
                 !       ghost values have an initial guess in them (extrapolate ?)
                 un(i1-0,i2-0,i3-is3,ex) = -(g - cm1*un(i1-0,i2-0,i3-is3,ex) )/cm1 
                 is =1-2*side2
                 dr0=dr(1)
                 rx0 = rsxy(i1,i2,i3,1,0)
                 ry0 = rsxy(i1,i2,i3,1,1)
                 rz0 = rsxy(i1,i2,i3,1,2)
                 rxNormSq = rx0**2 + ry0**2 + rz0**2
                 rxNorm = max( epsX, sqrt(rxNormSq) )
                 rxx0 = rsxyx23(i1,i2,i3,1,0)
                 ryy0 = rsxyy23(i1,i2,i3,1,1)
                 rzz0 = rsxyz23(i1,i2,i3,1,2)
                 ! cm1 : coeff of u(i1-0,i2-is2,i3-0,cc) in g (given below): 
                 cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0+rzz0)/(2.*dr0) ) )
                 ! derivatives at time t: 
                 ! derivatives at time t+dt : evaluate using the incorrect ghost values 
                 if( 1.eq.0 )then
                   ur0   =   ur2(i1,i2,i3,ex)
                   urr0  =  urr2(i1,i2,i3,ex)
                   unr0  =  unr2(i1,i2,i3,ex)
                   unrr0 = unrr2(i1,i2,i3,ex)
                 else if( 1.eq.1 )then
                   ur0   =   us2(i1,i2,i3,ex)
                   urr0  =  uss2(i1,i2,i3,ex)
                   unr0  =  uns2(i1,i2,i3,ex)
                   unrr0 = unss2(i1,i2,i3,ex)
                 else
                   ur0   =   ut2(i1,i2,i3,ex)
                   urr0  =  utt2(i1,i2,i3,ex)
                   unr0  =  unt2(i1,i2,i3,ex)
                   unrr0 = untt2(i1,i2,i3,ex)
                 end if
                 uxx0 = uxx23(i1,i2,i3,ex)
                 uyy0 = uyy23(i1,i2,i3,ex)
                 uzz0 = uzz23(i1,i2,i3,ex)
                 unxx0 = unxx23(i1,i2,i3,ex)
                 unyy0 = unyy23(i1,i2,i3,ex)
                 unzz0 = unzz23(i1,i2,i3,ex)
                 ! first evaluate the BC using the incorrect ghost values 
                 Dn2 = rxNormSq*(unrr0+urr0) 
                 Lu  = unxx0 + unyy0 + unzz0 + uxx0 + uyy0 + uzz0 
                 g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
                 ! note: this assumes an orthogonal grid -- we should make sure that the 
                 !       ghost values have an initial guess in them (extrapolate ?)
                 un(i1-0,i2-is2,i3-0,ex) = -(g - cm1*un(i1-0,i2-is2,i3-0,ex) )/cm1 
              end if
             else
              ! --- adjacent face is NOT another ABC ---
              ! *CHECK ME*
              ! Do this for now *wdh* Sept 20, 2016
              js1=0
              js2=0
              js3=is3
              un(i1-js1,i2-js2,i3-js3,ex)=(3.*un(i1,i2,i3,ex)-3.*un(i1+js1,i2+js2,i3+js3,ex)+un(i1+2*js1,i2+2*js2,i3+2*js3,ex))
              js1=0
              js2=is2
              js3=0
              un(i1-js1,i2-js2,i3-js3,ex)=(3.*un(i1,i2,i3,ex)-3.*un(i1+js1,i2+js2,i3+js3,ex)+un(i1+2*js1,i2+2*js2,i3+2*js3,ex))
             end if ! end adjacentFace
              do m3=0,numberOfGhostPoints
              do m2=0,numberOfGhostPoints
              do m1=0,0
               mSum = m1+m2+m3
               if( mSum.gt.0 .and. mSum.ne.1 )then ! mSum=1 : these points have already been assigned
                ! extrap ghost point (j1,j2,j3)
                j1=i1-is1*m1
                j2=i2-is2*m2
                j3=i3-is3*m3
                ! js1=0 if m1=0 and js1=is1 if m1>0
                js1 = is1*min(m1,1)
                js2 = is2*min(m2,1)
                js3 = is3*min(m3,1)
                if( orderOfAccuracy.eq.2 )then
                  ! Changed to third-order extra *wdh* Sept 18, 2016
                  un(j1,j2,j3,ex)=(3.*un(j1+js1,j2+js2,j3+js3,ex)-3.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex))
                  ! un(j1,j2,j3,ey)=extrap3(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                  ! un(j1,j2,j3,ez)=extrap3(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  ! un(j1,j2,j3,ex)=extrap2(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                  ! un(j1,j2,j3,ey)=extrap2(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                  ! un(j1,j2,j3,ez)=extrap2(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                else                                                                           
                  ! Note: adjust for incident fields should take into account the width of extrapolation: 
                  if( m1.le.1 .and. m2.le.1 .and. m3.le.1 )then
                    ! increased extrapolation to order=5 *wdh* June 20, 2016
                    un(j1,j2,j3,ex)=(5.*un(j1+js1,j2+js2,j3+js3,ex)-10.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+10.*un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex)-5.*un(j1+js1+3*js1,j2+js2+3*js2,j3+js3+3*js3,ex)+un(j1+js1+4*js1,j2+js2+4*js2,j3+js3+4*js3,ex))
                    ! un(j1,j2,j3,ey)=extrap5(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    ! un(j1,j2,j3,ez)=extrap5(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                    !un(j1,j2,j3,ex)=extrap3(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                    !un(j1,j2,j3,ey)=extrap3(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    !un(j1,j2,j3,ez)=extrap3(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  else
                    ! 2nd-ghost line 
                    un(j1,j2,j3,ex)=(5.*un(j1+js1,j2+js2,j3+js3,ex)-10.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+10.*un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex)-5.*un(j1+js1+3*js1,j2+js2+3*js2,j3+js3+3*js3,ex)+un(j1+js1+4*js1,j2+js2+4*js2,j3+js3+4*js3,ex))
                    ! un(j1,j2,j3,ey)=extrap5(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    ! un(j1,j2,j3,ez)=extrap5(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                    ! un(j1,j2,j3,ex)=extrap4(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                    ! un(j1,j2,j3,ey)=extrap4(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    ! un(j1,j2,j3,ez)=extrap4(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  end if
               end if
               end if
              end do
              end do
              end do
           end if
           end do
          else if( edgeDirection.eq.1 )then
           i1=n1a
           i3=n3a
           do i2=n2a,n2b
           if( mask(i1,i2,i3).gt.0 )then ! *wdh* 090712
             if( adjacentFaceIsABC )then
               ! Use first-order-in-time formula since it doesn't require other ghost point at new time (un)
               if( gridType.eq.rectangular )then
                 if( forcingOption.eq.twilightZoneForcing )then
                   ! Test: set to exact solution at time t:
                   ! x=xy(i1-is1,i2,i3,0)
                   ! y=xy(i1-is1,i2,i3,1)
                   ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                   ! un(i1-is1,i2,i3,ey)=eyTrue
                   ! add TZ forcing *wdh* Sept 17, 2016
                   ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                     ! ------ Cartesian Grid 3d forcing ----------
                     z=xy(i1,i2,i3,2)
                      ! Values for forcex(ex) are currently needed at corners:
                        call ogDeriv(ep, 1,1,0,0,x,y,z,tp,ex,utx)
                        call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                        call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                        call ogDeriv(ep, 0,0,0,2,x,y,z,tp,ex,uzz)
                      forcex(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) ) 
                      ! OGDERIV(1,1,0,0,x,y,z,tp,ey,utx)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                      ! OGDERIV(0,0,0,2,x,y,z,tp,ey,uzz)
                      ! forcex(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) ) 
                      ! OGDERIV(1,1,0,0,x,y,z,tp,ez,utx)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,ez,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,ez,uyy)
                      ! OGDERIV(0,0,0,2,x,y,z,tp,ez,uzz)
                      ! forcex(ez) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) )
                   ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                 end if
                un(i1-is1,i2,i3,ex)=(un(i1+is1,i2+0,i3+0,ex)-(u(i1+is1,i2+0,i3+0,ex)-u(i1-is1,i2-0,i3-0,ex))-(2.*dxa*dt)*(c1abcem2*uxx23r(i1,i2,i3,ex)+c2abcem2*(uyy23r(i1,i2,i3,ex)+uzz23r(i1,i2,i3,ex))+forcex(ex)))
                 if( forcingOption.eq.twilightZoneForcing )then
                   ! Test: set to exact solution at time t:
                   ! x=xy(i1-is1,i2,i3,0)
                   ! y=xy(i1-is1,i2,i3,1)
                   ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                   ! un(i1-is1,i2,i3,ey)=eyTrue
                   ! add TZ forcing *wdh* Sept 17, 2016
                   ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                     ! ------ Cartesian Grid 3d forcing ----------
                     z=xy(i1,i2,i3,2)
                        call ogDeriv(ep, 1,0,0,1,x,y,z,tp,ex,utz)
                        call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                        call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                        call ogDeriv(ep, 0,0,0,2,x,y,z,tp,ex,uzz)
                      forcez(ex) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) ) 
                      ! OGDERIV(1,0,0,1,x,y,z,tp,ey,utz)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                      ! OGDERIV(0,0,0,2,x,y,z,tp,ey,uzz)
                      ! forcez(ey) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) ) 
                      ! OGDERIV(1,0,0,1,x,y,z,tp,ez,utz)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,ez,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,ez,uyy)
                      ! OGDERIV(0,0,0,2,x,y,z,tp,ez,uzz)
                      ! forcez(ez) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) )
                   ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                 end if
                un(i1,i2,i3-is3,ex)=(un(i1+0,i2+0,i3+is3,ex)-(u(i1+0,i2+0,i3+is3,ex)-u(i1-0,i2-0,i3-is3,ex))-(2.*dza*dt)*(c1abcem2*uzz23r(i1,i2,i3,ex)+c2abcem2*(uxx23r(i1,i2,i3,ex)+uyy23r(i1,i2,i3,ex))+forcez(ex)))
               else
                  is =1-2*side1
                  dr0=dr(0)
                  rx0 = rsxy(i1,i2,i3,0,0)
                  ry0 = rsxy(i1,i2,i3,0,1)
                  rz0 = rsxy(i1,i2,i3,0,2)
                  rxNormSq = rx0**2 + ry0**2 + rz0**2
                  rxNorm = max( epsX, sqrt(rxNormSq) )
                  rxx0 = rsxyx23(i1,i2,i3,0,0)
                  ryy0 = rsxyy23(i1,i2,i3,0,1)
                  rzz0 = rsxyz23(i1,i2,i3,0,2)
                  ! cm1 : coeff of u(i1-is1,i2-0,i3-0,cc) in g (given below): 
                  cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0+rzz0)/(2.*dr0) ) )
                  ! derivatives at time t: 
                  ! derivatives at time t+dt : evaluate using the incorrect ghost values 
                  if( 0.eq.0 )then
                    ur0   =   ur2(i1,i2,i3,ex)
                    urr0  =  urr2(i1,i2,i3,ex)
                    unr0  =  unr2(i1,i2,i3,ex)
                    unrr0 = unrr2(i1,i2,i3,ex)
                  else if( 0.eq.1 )then
                    ur0   =   us2(i1,i2,i3,ex)
                    urr0  =  uss2(i1,i2,i3,ex)
                    unr0  =  uns2(i1,i2,i3,ex)
                    unrr0 = unss2(i1,i2,i3,ex)
                  else
                    ur0   =   ut2(i1,i2,i3,ex)
                    urr0  =  utt2(i1,i2,i3,ex)
                    unr0  =  unt2(i1,i2,i3,ex)
                    unrr0 = untt2(i1,i2,i3,ex)
                  end if
                  uxx0 = uxx23(i1,i2,i3,ex)
                  uyy0 = uyy23(i1,i2,i3,ex)
                  uzz0 = uzz23(i1,i2,i3,ex)
                  unxx0 = unxx23(i1,i2,i3,ex)
                  unyy0 = unyy23(i1,i2,i3,ex)
                  unzz0 = unzz23(i1,i2,i3,ex)
                  ! first evaluate the BC using the incorrect ghost values 
                  Dn2 = rxNormSq*(unrr0+urr0) 
                  Lu  = unxx0 + unyy0 + unzz0 + uxx0 + uyy0 + uzz0 
                  g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
                  ! note: this assumes an orthogonal grid -- we should make sure that the 
                  !       ghost values have an initial guess in them (extrapolate ?)
                  un(i1-is1,i2-0,i3-0,ex) = -(g - cm1*un(i1-is1,i2-0,i3-0,ex) )/cm1 
                  is =1-2*side3
                  dr0=dr(2)
                  rx0 = rsxy(i1,i2,i3,2,0)
                  ry0 = rsxy(i1,i2,i3,2,1)
                  rz0 = rsxy(i1,i2,i3,2,2)
                  rxNormSq = rx0**2 + ry0**2 + rz0**2
                  rxNorm = max( epsX, sqrt(rxNormSq) )
                  rxx0 = rsxyx23(i1,i2,i3,2,0)
                  ryy0 = rsxyy23(i1,i2,i3,2,1)
                  rzz0 = rsxyz23(i1,i2,i3,2,2)
                  ! cm1 : coeff of u(i1-0,i2-0,i3-is3,cc) in g (given below): 
                  cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0+rzz0)/(2.*dr0) ) )
                  ! derivatives at time t: 
                  ! derivatives at time t+dt : evaluate using the incorrect ghost values 
                  if( 2.eq.0 )then
                    ur0   =   ur2(i1,i2,i3,ex)
                    urr0  =  urr2(i1,i2,i3,ex)
                    unr0  =  unr2(i1,i2,i3,ex)
                    unrr0 = unrr2(i1,i2,i3,ex)
                  else if( 2.eq.1 )then
                    ur0   =   us2(i1,i2,i3,ex)
                    urr0  =  uss2(i1,i2,i3,ex)
                    unr0  =  uns2(i1,i2,i3,ex)
                    unrr0 = unss2(i1,i2,i3,ex)
                  else
                    ur0   =   ut2(i1,i2,i3,ex)
                    urr0  =  utt2(i1,i2,i3,ex)
                    unr0  =  unt2(i1,i2,i3,ex)
                    unrr0 = untt2(i1,i2,i3,ex)
                  end if
                  uxx0 = uxx23(i1,i2,i3,ex)
                  uyy0 = uyy23(i1,i2,i3,ex)
                  uzz0 = uzz23(i1,i2,i3,ex)
                  unxx0 = unxx23(i1,i2,i3,ex)
                  unyy0 = unyy23(i1,i2,i3,ex)
                  unzz0 = unzz23(i1,i2,i3,ex)
                  ! first evaluate the BC using the incorrect ghost values 
                  Dn2 = rxNormSq*(unrr0+urr0) 
                  Lu  = unxx0 + unyy0 + unzz0 + uxx0 + uyy0 + uzz0 
                  g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
                  ! note: this assumes an orthogonal grid -- we should make sure that the 
                  !       ghost values have an initial guess in them (extrapolate ?)
                  un(i1-0,i2-0,i3-is3,ex) = -(g - cm1*un(i1-0,i2-0,i3-is3,ex) )/cm1 
               end if
             else
              ! --- adjacent face is NOT another ABC ---
              ! *CHECK ME*
              ! Do this for now *wdh* Sept 20, 2016
              js1=is1
              js2=0
              js3=0
              un(i1-js1,i2-js2,i3-js3,ex)=(3.*un(i1,i2,i3,ex)-3.*un(i1+js1,i2+js2,i3+js3,ex)+un(i1+2*js1,i2+2*js2,i3+2*js3,ex))
              js1=0
              js2=0
              js3=is3
              un(i1-js1,i2-js2,i3-js3,ex)=(3.*un(i1,i2,i3,ex)-3.*un(i1+js1,i2+js2,i3+js3,ex)+un(i1+2*js1,i2+2*js2,i3+2*js3,ex))
             end if ! end adjacentFace
              do m3=0,numberOfGhostPoints
              do m2=0,0
              do m1=0,numberOfGhostPoints
               mSum = m1+m2+m3
               if( mSum.gt.0 .and. mSum.ne.1 )then ! mSum=1 : these points have already been assigned
                ! extrap ghost point (j1,j2,j3)
                j1=i1-is1*m1
                j2=i2-is2*m2
                j3=i3-is3*m3
                ! js1=0 if m1=0 and js1=is1 if m1>0
                js1 = is1*min(m1,1)
                js2 = is2*min(m2,1)
                js3 = is3*min(m3,1)
                if( orderOfAccuracy.eq.2 )then
                  ! Changed to third-order extra *wdh* Sept 18, 2016
                  un(j1,j2,j3,ex)=(3.*un(j1+js1,j2+js2,j3+js3,ex)-3.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex))
                  ! un(j1,j2,j3,ey)=extrap3(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                  ! un(j1,j2,j3,ez)=extrap3(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  ! un(j1,j2,j3,ex)=extrap2(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                  ! un(j1,j2,j3,ey)=extrap2(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                  ! un(j1,j2,j3,ez)=extrap2(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                else                                                                           
                  ! Note: adjust for incident fields should take into account the width of extrapolation: 
                  if( m1.le.1 .and. m2.le.1 .and. m3.le.1 )then
                    ! increased extrapolation to order=5 *wdh* June 20, 2016
                    un(j1,j2,j3,ex)=(5.*un(j1+js1,j2+js2,j3+js3,ex)-10.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+10.*un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex)-5.*un(j1+js1+3*js1,j2+js2+3*js2,j3+js3+3*js3,ex)+un(j1+js1+4*js1,j2+js2+4*js2,j3+js3+4*js3,ex))
                    ! un(j1,j2,j3,ey)=extrap5(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    ! un(j1,j2,j3,ez)=extrap5(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                    !un(j1,j2,j3,ex)=extrap3(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                    !un(j1,j2,j3,ey)=extrap3(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    !un(j1,j2,j3,ez)=extrap3(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  else
                    ! 2nd-ghost line 
                    un(j1,j2,j3,ex)=(5.*un(j1+js1,j2+js2,j3+js3,ex)-10.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+10.*un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex)-5.*un(j1+js1+3*js1,j2+js2+3*js2,j3+js3+3*js3,ex)+un(j1+js1+4*js1,j2+js2+4*js2,j3+js3+4*js3,ex))
                    ! un(j1,j2,j3,ey)=extrap5(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    ! un(j1,j2,j3,ez)=extrap5(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                    ! un(j1,j2,j3,ex)=extrap4(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                    ! un(j1,j2,j3,ey)=extrap4(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    ! un(j1,j2,j3,ez)=extrap4(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  end if
               end if
               end if
              end do
              end do
              end do
           end if
           end do
          else if( edgeDirection.eq.2 )then
           ! write(*,'(" ABC:set corner: grid,side1,side2,i1,i2=",3i3,2i5)') grid,side1,side2,i1,i2
           i1=n1a
           i2=n2a
           do i3=n3a,n3b
           if( mask(i1,i2,i3).gt.0 )then ! *wdh* 090712
             if( adjacentFaceIsABC )then
               ! Use first-order-in-time formula since it doesn't require other ghost point at new time (un)
               if( gridType.eq.rectangular )then
                 if( forcingOption.eq.twilightZoneForcing )then
                   ! Test: set to exact solution at time t:
                   ! x=xy(i1-is1,i2,i3,0)
                   ! y=xy(i1-is1,i2,i3,1)
                   ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                   ! un(i1-is1,i2,i3,ey)=eyTrue
                   ! add TZ forcing *wdh* Sept 17, 2016
                   ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                     ! ------ Cartesian Grid 3d forcing ----------
                     z=xy(i1,i2,i3,2)
                      ! Values for forcex(ex) are currently needed at corners:
                        call ogDeriv(ep, 1,1,0,0,x,y,z,tp,ex,utx)
                        call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                        call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                        call ogDeriv(ep, 0,0,0,2,x,y,z,tp,ex,uzz)
                      forcex(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) ) 
                      ! OGDERIV(1,1,0,0,x,y,z,tp,ey,utx)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                      ! OGDERIV(0,0,0,2,x,y,z,tp,ey,uzz)
                      ! forcex(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) ) 
                      ! OGDERIV(1,1,0,0,x,y,z,tp,ez,utx)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,ez,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,ez,uyy)
                      ! OGDERIV(0,0,0,2,x,y,z,tp,ez,uzz)
                      ! forcex(ez) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) )
                   ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                 end if
                un(i1-is1,i2,i3,ex)=(un(i1+is1,i2+0,i3+0,ex)-(u(i1+is1,i2+0,i3+0,ex)-u(i1-is1,i2-0,i3-0,ex))-(2.*dxa*dt)*(c1abcem2*uxx23r(i1,i2,i3,ex)+c2abcem2*(uyy23r(i1,i2,i3,ex)+uzz23r(i1,i2,i3,ex))+forcex(ex)))
                 if( forcingOption.eq.twilightZoneForcing )then
                   ! Test: set to exact solution at time t:
                   ! x=xy(i1-is1,i2,i3,0)
                   ! y=xy(i1-is1,i2,i3,1)
                   ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                   ! un(i1-is1,i2,i3,ey)=eyTrue
                   ! add TZ forcing *wdh* Sept 17, 2016
                   ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                   x=xy(i1,i2,i3,0)
                   y=xy(i1,i2,i3,1)
                     ! ------ Cartesian Grid 3d forcing ----------
                     z=xy(i1,i2,i3,2)
                        call ogDeriv(ep, 1,0,1,0,x,y,z,tp,ex,uty)
                        call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                        call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                        call ogDeriv(ep, 0,0,0,2,x,y,z,tp,ex,uzz)
                      forcey(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) ) 
                      ! OGDERIV(1,0,1,0,x,y,z,tp,ey,uty)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                      ! OGDERIV(0,0,0,2,x,y,z,tp,ey,uzz)
                      ! forcey(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) ) 
                      ! OGDERIV(1,0,1,0,x,y,z,tp,ez,uty)
                      ! OGDERIV(0,2,0,0,x,y,z,tp,ez,uxx)
                      ! OGDERIV(0,0,2,0,x,y,z,tp,ez,uyy)
                      ! OGDERIV(0,0,0,2,x,y,z,tp,ez,uzz)
                      ! forcey(ez) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) )
                   ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                 end if
                un(i1,i2-is2,i3,ex)=(un(i1+0,i2+is2,i3+0,ex)-(u(i1+0,i2+is2,i3+0,ex)-u(i1-0,i2-is2,i3-0,ex))-(2.*dya*dt)*(c1abcem2*uyy23r(i1,i2,i3,ex)+c2abcem2*(uxx23r(i1,i2,i3,ex)+uzz23r(i1,i2,i3,ex))+forcey(ex)))
               else
                  is =1-2*side1
                  dr0=dr(0)
                  rx0 = rsxy(i1,i2,i3,0,0)
                  ry0 = rsxy(i1,i2,i3,0,1)
                  rz0 = rsxy(i1,i2,i3,0,2)
                  rxNormSq = rx0**2 + ry0**2 + rz0**2
                  rxNorm = max( epsX, sqrt(rxNormSq) )
                  rxx0 = rsxyx23(i1,i2,i3,0,0)
                  ryy0 = rsxyy23(i1,i2,i3,0,1)
                  rzz0 = rsxyz23(i1,i2,i3,0,2)
                  ! cm1 : coeff of u(i1-is1,i2-0,i3-0,cc) in g (given below): 
                  cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0+rzz0)/(2.*dr0) ) )
                  ! derivatives at time t: 
                  ! derivatives at time t+dt : evaluate using the incorrect ghost values 
                  if( 0.eq.0 )then
                    ur0   =   ur2(i1,i2,i3,ex)
                    urr0  =  urr2(i1,i2,i3,ex)
                    unr0  =  unr2(i1,i2,i3,ex)
                    unrr0 = unrr2(i1,i2,i3,ex)
                  else if( 0.eq.1 )then
                    ur0   =   us2(i1,i2,i3,ex)
                    urr0  =  uss2(i1,i2,i3,ex)
                    unr0  =  uns2(i1,i2,i3,ex)
                    unrr0 = unss2(i1,i2,i3,ex)
                  else
                    ur0   =   ut2(i1,i2,i3,ex)
                    urr0  =  utt2(i1,i2,i3,ex)
                    unr0  =  unt2(i1,i2,i3,ex)
                    unrr0 = untt2(i1,i2,i3,ex)
                  end if
                  uxx0 = uxx23(i1,i2,i3,ex)
                  uyy0 = uyy23(i1,i2,i3,ex)
                  uzz0 = uzz23(i1,i2,i3,ex)
                  unxx0 = unxx23(i1,i2,i3,ex)
                  unyy0 = unyy23(i1,i2,i3,ex)
                  unzz0 = unzz23(i1,i2,i3,ex)
                  ! first evaluate the BC using the incorrect ghost values 
                  Dn2 = rxNormSq*(unrr0+urr0) 
                  Lu  = unxx0 + unyy0 + unzz0 + uxx0 + uyy0 + uzz0 
                  g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
                  ! note: this assumes an orthogonal grid -- we should make sure that the 
                  !       ghost values have an initial guess in them (extrapolate ?)
                  un(i1-is1,i2-0,i3-0,ex) = -(g - cm1*un(i1-is1,i2-0,i3-0,ex) )/cm1 
                  is =1-2*side2
                  dr0=dr(1)
                  rx0 = rsxy(i1,i2,i3,1,0)
                  ry0 = rsxy(i1,i2,i3,1,1)
                  rz0 = rsxy(i1,i2,i3,1,2)
                  rxNormSq = rx0**2 + ry0**2 + rz0**2
                  rxNorm = max( epsX, sqrt(rxNormSq) )
                  rxx0 = rsxyx23(i1,i2,i3,1,0)
                  ryy0 = rsxyy23(i1,i2,i3,1,1)
                  rzz0 = rsxyz23(i1,i2,i3,1,2)
                  ! cm1 : coeff of u(i1-0,i2-is2,i3-0,cc) in g (given below): 
                  cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0+rzz0)/(2.*dr0) ) )
                  ! derivatives at time t: 
                  ! derivatives at time t+dt : evaluate using the incorrect ghost values 
                  if( 1.eq.0 )then
                    ur0   =   ur2(i1,i2,i3,ex)
                    urr0  =  urr2(i1,i2,i3,ex)
                    unr0  =  unr2(i1,i2,i3,ex)
                    unrr0 = unrr2(i1,i2,i3,ex)
                  else if( 1.eq.1 )then
                    ur0   =   us2(i1,i2,i3,ex)
                    urr0  =  uss2(i1,i2,i3,ex)
                    unr0  =  uns2(i1,i2,i3,ex)
                    unrr0 = unss2(i1,i2,i3,ex)
                  else
                    ur0   =   ut2(i1,i2,i3,ex)
                    urr0  =  utt2(i1,i2,i3,ex)
                    unr0  =  unt2(i1,i2,i3,ex)
                    unrr0 = untt2(i1,i2,i3,ex)
                  end if
                  uxx0 = uxx23(i1,i2,i3,ex)
                  uyy0 = uyy23(i1,i2,i3,ex)
                  uzz0 = uzz23(i1,i2,i3,ex)
                  unxx0 = unxx23(i1,i2,i3,ex)
                  unyy0 = unyy23(i1,i2,i3,ex)
                  unzz0 = unzz23(i1,i2,i3,ex)
                  ! first evaluate the BC using the incorrect ghost values 
                  Dn2 = rxNormSq*(unrr0+urr0) 
                  Lu  = unxx0 + unyy0 + unzz0 + uxx0 + uyy0 + uzz0 
                  g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
                  ! note: this assumes an orthogonal grid -- we should make sure that the 
                  !       ghost values have an initial guess in them (extrapolate ?)
                  un(i1-0,i2-is2,i3-0,ex) = -(g - cm1*un(i1-0,i2-is2,i3-0,ex) )/cm1 
               end if
             else
              ! --- adjacent face is NOT another ABC ---
              ! *CHECK ME*
              ! Do this for now *wdh* Sept 20, 2016
              js1=is1
              js2=0
              js3=0
              un(i1-js1,i2-js2,i3-js3,ex)=(3.*un(i1,i2,i3,ex)-3.*un(i1+js1,i2+js2,i3+js3,ex)+un(i1+2*js1,i2+2*js2,i3+2*js3,ex))
              js1=0
              js2=is2
              js3=0
              un(i1-js1,i2-js2,i3-js3,ex)=(3.*un(i1,i2,i3,ex)-3.*un(i1+js1,i2+js2,i3+js3,ex)+un(i1+2*js1,i2+2*js2,i3+2*js3,ex))
             end if ! end adjacentFace
              do m3=0,0
              do m2=0,numberOfGhostPoints
              do m1=0,numberOfGhostPoints
               mSum = m1+m2+m3
               if( mSum.gt.0 .and. mSum.ne.1 )then ! mSum=1 : these points have already been assigned
                ! extrap ghost point (j1,j2,j3)
                j1=i1-is1*m1
                j2=i2-is2*m2
                j3=i3-is3*m3
                ! js1=0 if m1=0 and js1=is1 if m1>0
                js1 = is1*min(m1,1)
                js2 = is2*min(m2,1)
                js3 = is3*min(m3,1)
                if( orderOfAccuracy.eq.2 )then
                  ! Changed to third-order extra *wdh* Sept 18, 2016
                  un(j1,j2,j3,ex)=(3.*un(j1+js1,j2+js2,j3+js3,ex)-3.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex))
                  ! un(j1,j2,j3,ey)=extrap3(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                  ! un(j1,j2,j3,ez)=extrap3(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  ! un(j1,j2,j3,ex)=extrap2(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                  ! un(j1,j2,j3,ey)=extrap2(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                  ! un(j1,j2,j3,ez)=extrap2(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                else                                                                           
                  ! Note: adjust for incident fields should take into account the width of extrapolation: 
                  if( m1.le.1 .and. m2.le.1 .and. m3.le.1 )then
                    ! increased extrapolation to order=5 *wdh* June 20, 2016
                    un(j1,j2,j3,ex)=(5.*un(j1+js1,j2+js2,j3+js3,ex)-10.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+10.*un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex)-5.*un(j1+js1+3*js1,j2+js2+3*js2,j3+js3+3*js3,ex)+un(j1+js1+4*js1,j2+js2+4*js2,j3+js3+4*js3,ex))
                    ! un(j1,j2,j3,ey)=extrap5(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    ! un(j1,j2,j3,ez)=extrap5(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                    !un(j1,j2,j3,ex)=extrap3(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                    !un(j1,j2,j3,ey)=extrap3(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    !un(j1,j2,j3,ez)=extrap3(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  else
                    ! 2nd-ghost line 
                    un(j1,j2,j3,ex)=(5.*un(j1+js1,j2+js2,j3+js3,ex)-10.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+10.*un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex)-5.*un(j1+js1+3*js1,j2+js2+3*js2,j3+js3+3*js3,ex)+un(j1+js1+4*js1,j2+js2+4*js2,j3+js3+4*js3,ex))
                    ! un(j1,j2,j3,ey)=extrap5(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    ! un(j1,j2,j3,ez)=extrap5(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                    ! un(j1,j2,j3,ex)=extrap4(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                    ! un(j1,j2,j3,ey)=extrap4(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                    ! un(j1,j2,j3,ez)=extrap4(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  end if
               end if
               end if
              end do
              end do
              end do
           end if
           end do
          end if ! end if edgeDirection
         end if ! bc
         end do ! end sideb
         end do ! end sidea
         end do ! end edgeDirection
        ! ***** vertices in 3D  *****
        do side3=0,1
        do side2=0,1
        do side1=0,1
         bc1=boundaryCondition(side1,0)
         bc2=boundaryCondition(side2,1)
         bc3=boundaryCondition(side3,2)
         if( ( (bc1.eq.abcEM2 .or. bc1.eq.absorbing) .or. (bc2.eq.abcEM2 .or. bc2.eq.absorbing) .or. (bc3.eq.abcEM2 .or. bc3.eq.absorbing) ) .and.  ( bc1.gt.0 .and. bc2.gt.0 .and. bc3.gt.0 ) )then
          ! Three physical faces meet at this corner and at least one face is an ABC
          i1=gridIndexRange(side1,0)
          i2=gridIndexRange(side2,1)
          i3=gridIndexRange(side3,2)
          if( mask(i1,i2,i3).gt.0 )then 
           is1=1-2*side1
           is2=1-2*side2
           is3=1-2*side3
            do m3=0,numberOfGhostPoints
            do m2=0,numberOfGhostPoints
            do m1=0,numberOfGhostPoints
             mSum = m1+m2+m3
             if( mSum.gt.0 .and. mSum.ne.1 )then ! mSum=1 : these points have already been assigned
              ! extrap ghost point (j1,j2,j3)
              j1=i1-is1*m1
              j2=i2-is2*m2
              j3=i3-is3*m3
              ! js1=0 if m1=0 and js1=is1 if m1>0
              js1 = is1*min(m1,1)
              js2 = is2*min(m2,1)
              js3 = is3*min(m3,1)
              if( orderOfAccuracy.eq.2 )then
                ! Changed to third-order extra *wdh* Sept 18, 2016
                un(j1,j2,j3,ex)=(3.*un(j1+js1,j2+js2,j3+js3,ex)-3.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex))
                ! un(j1,j2,j3,ey)=extrap3(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                ! un(j1,j2,j3,ez)=extrap3(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                ! un(j1,j2,j3,ex)=extrap2(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                ! un(j1,j2,j3,ey)=extrap2(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                ! un(j1,j2,j3,ez)=extrap2(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
              else                                                                           
                ! Note: adjust for incident fields should take into account the width of extrapolation: 
                if( m1.le.1 .and. m2.le.1 .and. m3.le.1 )then
                  ! increased extrapolation to order=5 *wdh* June 20, 2016
                  un(j1,j2,j3,ex)=(5.*un(j1+js1,j2+js2,j3+js3,ex)-10.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+10.*un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex)-5.*un(j1+js1+3*js1,j2+js2+3*js2,j3+js3+3*js3,ex)+un(j1+js1+4*js1,j2+js2+4*js2,j3+js3+4*js3,ex))
                  ! un(j1,j2,j3,ey)=extrap5(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                  ! un(j1,j2,j3,ez)=extrap5(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  !un(j1,j2,j3,ex)=extrap3(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                  !un(j1,j2,j3,ey)=extrap3(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                  !un(j1,j2,j3,ez)=extrap3(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                else
                  ! 2nd-ghost line 
                  un(j1,j2,j3,ex)=(5.*un(j1+js1,j2+js2,j3+js3,ex)-10.*un(j1+js1+js1,j2+js2+js2,j3+js3+js3,ex)+10.*un(j1+js1+2*js1,j2+js2+2*js2,j3+js3+2*js3,ex)-5.*un(j1+js1+3*js1,j2+js2+3*js2,j3+js3+3*js3,ex)+un(j1+js1+4*js1,j2+js2+4*js2,j3+js3+4*js3,ex))
                  ! un(j1,j2,j3,ey)=extrap5(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                  ! un(j1,j2,j3,ez)=extrap5(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                  ! un(j1,j2,j3,ex)=extrap4(un,j1+js1,j2+js2,j3+js3,ex,js1,js2,js3)
                  ! un(j1,j2,j3,ey)=extrap4(un,j1+js1,j2+js2,j3+js3,ey,js1,js2,js3)
                  ! un(j1,j2,j3,ez)=extrap4(un,j1+js1,j2+js2,j3+js3,ez,js1,js2,js3)
                end if
             end if
             end if
            end do
            end do
            end do
          end if ! mask
         end if
        end do
        end do
        end do
       !  end if
      end if

  else if( orderOfAccuracy==4 )then

    ! write(*,*) 'abcWave: assign edges and corners, order 4'
      ! ------------------------------------------------------------------------
      ! ------------------Corners-----------------------------------------------
      ! ------------------------------------------------------------------------
      ! We need to assign points in the corner region:
      !
      !           |  |  |
      !           +--+--+--
      !           |  |  |
      !     D--A--X--+--+--
      !     |  |  |
      !     D--C--B
      !     |  |  |
      !     D--D--D
      ! 
      if( nd.eq.2 )then
        ! **** 2D ****
        i3=gridIndexRange(0,2)
        is3=0
        do side1=0,1
        do side2=0,1
         bc1=boundaryCondition(side1,0)
         bc2=boundaryCondition(side2,1)
         if( ((bc1.eq.abcEM2 .or. bc1.eq.absorbing) .and. bc2.gt.0 ) .or. ((bc2.eq.abcEM2 .or. bc2.eq.absorbing) .and. bc1.gt.0 ) )then
           ! --- One of the faces at this corner is an ABC and the other has bc>0 ---         
           ! Adjacent side is also an ABC: 
           adjacentFaceIsABC = bc1.eq.abcEM2 .or. bc1.eq.absorbing .and. bc2.eq.abcEM2 .or. bc2.eq.absorbing
           i1=gridIndexRange(side1,0) ! (i1,i2,i3)=corner point
           i2=gridIndexRange(side2,1)
           if( mask(i1,i2,i3).gt.0 )then ! *wdh* 090712
            ! write(*,'(" ABC:set corner order 4: grid,side1,side2,i1,i2=",3i3,2i5)') grid,side1,side2,i1,i2
            ! --- start by extrapolating all points on the extended boundary and adjacent to the corner ---
            !* is1=1-2*side1
            !* is2=1-2*side2
            !* is3=0
            !* j3=0
            !* extrapolateGhost(ex,ey,hz,numberOfGhostPoints,numberOfGhostPoints,0)
            ! Solve: 
            !   f1 = is1* u_xt - c( u_xx + .5*u_yy )
            !   f2 = (D+x)^5 U_{-2} = 0 
            !   f3 = is2* u_yt - c( u_yy + .5*u_xx )
            !   f4 = (D+y)^5 U_{-2} = 0 
            !
            ! Evaluate the equations with current ghost 
            ! There four coupled equations we need to solve for ghost points A,B,C,D below
            !                   |
            !                   +   X = corner 
            !                   |
            !             B--A--X---+---+
            !                   |
            !             G  G  C
            !                   |
            !             G  G  D 
            !
            !  f(u)  = [ f(u_old) - A (u_old) ] + A u = 0 
            !
            !  [ a11 a12 a13 a14 ][ uA ]   [ a11 a12 a13 a14 ][ uA_old ]   [ f1(u_old) ]
            !  [ a21 a22 a23 a24 ][ uB ] = [ a21 a22 a23 a24 ][ uB_old ] - [ f2(u_old) ]  
            !  [ a31 a32 a33 a34 ][ uC ]   [ a31 a32 a33 a34 ][ uC_old ]   [ f3(u_old) ]   
            !  [ a41 a42 a43 a44 ][ uD ]   [ a41 a42 a43 a44 ][ uD_old ]   [ f4(u_old) ]   
            if( adjacentFaceIsABC )then
              is1=1-2*side1
              is2=1-2*side2
              is3 = 0 
              if( gridType.eq.rectangular )then
                  ! getForcingEM2(X,2,tm,is1,is2,forcex)
                   if( forcingOption.eq.twilightZoneForcing )then
                     ! Test: set to exact solution at time t:
                     ! x=xy(i1-is1,i2,i3,0)
                     ! y=xy(i1-is1,i2,i3,1)
                     ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                     ! un(i1-is1,i2,i3,ey)=eyTrue
                     ! add TZ forcing *wdh* Sept 17, 2016
                     ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                     x=xy(i1,i2,i3,0)
                     y=xy(i1,i2,i3,1)
                        ! Values for forcep(ex) are currently needed at corners:
                          call ogDeriv(ep, 1,1,0,0,x,y,z,tp,ex,utx)
                          call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                          call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                        forcep(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                        ! OGDERIV(1,1,0,0,x,y,z,tp,ey,utx)
                        ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                        ! forcep(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                        ! OGDERIV(1,1,0,0,x,y,z,tp,hz,utx)
                        ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
                        ! forcep(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
                     ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                   end if
                   if( forcingOption.eq.twilightZoneForcing )then
                     ! Test: set to exact solution at time t:
                     ! x=xy(i1-is1,i2,i3,0)
                     ! y=xy(i1-is1,i2,i3,1)
                     ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                     ! un(i1-is1,i2,i3,ey)=eyTrue
                     ! add TZ forcing *wdh* Sept 17, 2016
                     ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                     x=xy(i1,i2,i3,0)
                     y=xy(i1,i2,i3,1)
                        ! Values for forcef(ex) are currently needed at corners:
                          call ogDeriv(ep, 1,1,0,0,x,y,z,t,ex,utx)
                          call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                          call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                        forcef(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                        ! OGDERIV(1,1,0,0,x,y,z,t,ey,utx)
                        ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
                        ! forcef(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                        ! OGDERIV(1,1,0,0,x,y,z,t,hz,utx)
                        ! OGDERIV(0,2,0,0,x,y,z,t,hz,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,t,hz,uyy)
                        ! forcef(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
                     ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                   end if
                  ! do idir=0,2
                  do idir=0,0 ! ex only 
                    forcex(idir)=.5*(forcep(idir)+forcef(idir))
                  end do
                  ! getForcingEM2(Y,2,tm,is1,is2,forcey)
                   if( forcingOption.eq.twilightZoneForcing )then
                     ! Test: set to exact solution at time t:
                     ! x=xy(i1-is1,i2,i3,0)
                     ! y=xy(i1-is1,i2,i3,1)
                     ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                     ! un(i1-is1,i2,i3,ey)=eyTrue
                     ! add TZ forcing *wdh* Sept 17, 2016
                     ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                     x=xy(i1,i2,i3,0)
                     y=xy(i1,i2,i3,1)
                          call ogDeriv(ep, 1,0,1,0,x,y,z,tp,ex,uty)
                          call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                          call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                        forcep(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                        ! ! Values for forcep(ey) are currently needed at corners:
                        ! OGDERIV(1,0,1,0,x,y,z,tp,ey,uty)
                        ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                        ! forcep(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                        ! OGDERIV(1,0,1,0,x,y,z,tp,hz,uty)
                        ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
                        ! forcep(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
                     ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                   end if
                   if( forcingOption.eq.twilightZoneForcing )then
                     ! Test: set to exact solution at time t:
                     ! x=xy(i1-is1,i2,i3,0)
                     ! y=xy(i1-is1,i2,i3,1)
                     ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
                     ! un(i1-is1,i2,i3,ey)=eyTrue
                     ! add TZ forcing *wdh* Sept 17, 2016
                     ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
                     x=xy(i1,i2,i3,0)
                     y=xy(i1,i2,i3,1)
                          call ogDeriv(ep, 1,0,1,0,x,y,z,t,ex,uty)
                          call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                          call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                        forcef(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                        ! ! Values for forcef(ey) are currently needed at corners:
                        ! OGDERIV(1,0,1,0,x,y,z,t,ey,uty)
                        ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
                        ! forcef(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                        ! OGDERIV(1,0,1,0,x,y,z,t,hz,uty)
                        ! OGDERIV(0,2,0,0,x,y,z,t,hz,uxx)
                        ! OGDERIV(0,0,2,0,x,y,z,t,hz,uyy)
                        ! forcef(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
                     ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
                   end if
                  ! do idir=0,2
                  do idir=0,0 ! ex only 
                    forcey(idir)=.5*(forcep(idir)+forcef(idir))
                  end do
                ! Evaluate the residuals in the corner equations using current ghos values 
                ! evalCornerResidual(ORDER,DIM)
                  resv(1) = is1*(unx42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ex))/(dt)- .5*( c1abcem2*unxx42r(i1,i2,i3,ex) + c2abcem2*unyy42r(i1,i2,i3,ex) + c1abcem2* uxx42r(i1,i2,i3,ex) + c2abcem2* uyy42r(i1,i2,i3,ex) ) - forcex(ex) 
                  resv(2) = is2*(uny42r(i1,i2,i3,ex)-uy42r(i1,i2,i3,ex))/(dt)- .5*( c1abcem2*unyy42r(i1,i2,i3,ex) + c2abcem2*unxx42r(i1,i2,i3,ex) + c1abcem2* uyy42r(i1,i2,i3,ex) + c2abcem2* uxx42r(i1,i2,i3,ex) ) - forcey(ex)  
                  resv(3) =    un(i1-2*is1,i2,i3,ex) -5.*un(i1-1*is1,i2,i3,ex) +10.*un(i1+0*is1,i2,i3,ex) -10.*un(i1+1*is1,i2,i3,ex) +5.*un(i1+2*is1,i2,i3,ex) -un(i1+3*is1,i2,i3,ex)
                  resv(4) =    un(i1,i2-2*is2,i3,ex) -5.*un(i1,i2-1*is2,i3,ex) +10.*un(i1,i2+0*is2,i3,ex) -10.*un(i1,i2+1*is2,i3,ex) +5.*un(i1,i2+2*is2,i3,ex) -un(i1,i2+3*is2,i3,ex)
                a4(1,1) = -8./(12.*dt*dx(axis1))  - .5*c1abcem2*(16.)/(12.*dx(axis1)**2) ! coeff of uA in r1
                a4(1,2) =  1./(12.*dt*dx(axis1))  - .5*c1abcem2*(-1.)/(12.*dx(axis1)**2) ! coeff of uB in r1
                a4(1,3) =                         - .5*c2abcem2*(16.)/(12.*dx(axis2)**2) ! coeff of uC in r1
                a4(1,4) =                         - .5*c2abcem2*(-1.)/(12.*dx(axis2)**2) ! coeff of uD in r1
                a4(2,1) =                         - .5*c2abcem2*(16.)/(12.*dx(axis1)**2) ! coeff of uC in r2
                a4(2,2) =                         - .5*c2abcem2*(-1.)/(12.*dx(axis1)**2) ! coeff of uD in r2
                a4(2,3) = -8./(12.*dt*dx(axis2))  - .5*c1abcem2*(16.)/(12.*dx(axis2)**2) ! coeff of uA in r2
                a4(2,4) =  1./(12.*dt*dx(axis2))  - .5*c1abcem2*(-1.)/(12.*dx(axis2)**2) ! coeff of uB in r2
                a4(3,1) = -5.                                                    ! coeff of uA in r1
                a4(3,2) =  1.                                                    ! coeff of uB in r2
                a4(3,3)  = 0.
                a4(3,4)  = 0.
                a4(4,1) =  0.                                                    ! coeff of uA in r1
                a4(4,2) =  0.                                                    ! coeff of uB in r2
                a4(4,3) = -5.
                a4(4,4) =  1.
                uv(1) = un(i1-1*is1,i2      ,i3,ex)   ! uA
                uv(2) = un(i1-2*is1,i2      ,i3,ex)   ! uB
                uv(3) = un(i1      ,i2-1*is2,i3,ex)   ! uC
                uv(4) = un(i1      ,i2-2*is2,i3,ex)   ! uD
                do m=1,4
                  b4(m) = a4(m,1)*uv(1) + a4(m,2)*uv(2) + a4(m,3)*uv(3) + a4(m,4)*uv(4) - resv(m)
                end do
                ! Solve A4 * x = b4
                numberOfEquations=4
                call dgeco( a4(1,1), numberOfEquations, numberOfEquations, ipvt4(1), rcond, work4(1) )
                job=0
                call dgesl( a4(1,1), numberOfEquations, numberOfEquations, ipvt4(1), b4(1), job )
                un(i1-1*is1,i2      ,i3,ex) = b4(1)            
                un(i1-2*is1,i2      ,i3,ex) = b4(2)
                un(i1      ,i2-1*is2,i3,ex) = b4(3)
                un(i1      ,i2-2*is2,i3,ex) = b4(4)
                if( .false. )then
                  ! check the residuals after solve
                    resv(1) = is1*(unx42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ex))/(dt)- .5*( c1abcem2*unxx42r(i1,i2,i3,ex) + c2abcem2*unyy42r(i1,i2,i3,ex) + c1abcem2* uxx42r(i1,i2,i3,ex) + c2abcem2* uyy42r(i1,i2,i3,ex) ) - forcex(ex) 
                    resv(2) = is2*(uny42r(i1,i2,i3,ex)-uy42r(i1,i2,i3,ex))/(dt)- .5*( c1abcem2*unyy42r(i1,i2,i3,ex) + c2abcem2*unxx42r(i1,i2,i3,ex) + c1abcem2* uyy42r(i1,i2,i3,ex) + c2abcem2* uxx42r(i1,i2,i3,ex) ) - forcey(ex)  
                    resv(3) =    un(i1-2*is1,i2,i3,ex) -5.*un(i1-1*is1,i2,i3,ex) +10.*un(i1+0*is1,i2,i3,ex) -10.*un(i1+1*is1,i2,i3,ex) +5.*un(i1+2*is1,i2,i3,ex) -un(i1+3*is1,i2,i3,ex)
                    resv(4) =    un(i1,i2-2*is2,i3,ex) -5.*un(i1,i2-1*is2,i3,ex) +10.*un(i1,i2+0*is2,i3,ex) -10.*un(i1,i2+1*is2,i3,ex) +5.*un(i1,i2+2*is2,i3,ex) -un(i1,i2+3*is2,i3,ex)
                  write(*,'("abcWave: EM order 4: i1,i2=",2i4," corner residuals=",4(1pe10.2,1x))') i1,i2,(resv(m),m=1,4)
                end if
              else
               ! --- curvilinear grid ----
               stop 4449
               side=side1
               axis=0 
                 is =1-2*side
                 dr0=dr(axis)
                 rx0 = rsxy(i1,i2,i3,axis,0)
                 ry0 = rsxy(i1,i2,i3,axis,1)
                 rxNormSq = rx0**2 + ry0**2 
                 rxNorm = max( epsX, sqrt(rxNormSq) )
                 rxx0 = rsxyx22(i1,i2,i3,axis,0)
                 ryy0 = rsxyy22(i1,i2,i3,axis,1)
                 ! cm1 : coeff of u(i1-is1,i2-is2,i3-is3,cc) in g (given below): 
                 cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0)/(2.*dr0) ) )
                 ! u: derivatives at time t: 
                 ! un: derivatives at time t+dt : evaluate using the incorrect ghost values 
                 if( axis.eq.0 )then
                   ur0   =   ur2(i1,i2,i3,ex)
                   urr0  =  urr2(i1,i2,i3,ex)
                   unr0  =  unr2(i1,i2,i3,ex)
                   unrr0 = unrr2(i1,i2,i3,ex)
                 else
                   ur0   =   us2(i1,i2,i3,ex)
                   urr0  =  uss2(i1,i2,i3,ex)
                   unr0  =  uns2(i1,i2,i3,ex)
                   unrr0 = unss2(i1,i2,i3,ex)
                 end if
                 uxx0 = uxx22(i1,i2,i3,ex)
                 uyy0 = uyy22(i1,i2,i3,ex)
                 unxx0 = unxx22(i1,i2,i3,ex)
                 unyy0 = unyy22(i1,i2,i3,ex)
                 ! first evaluate the BC using the incorrect ghost values 
                 Dn2 = rxNormSq*(unrr0+urr0) 
                 Lu  = unxx0+unyy0 + uxx0+uyy0 
                 g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
                 ! note: this assumes an orthogonal grid -- we should make sure that the 
                 !       ghost values have an initial guess in them (extrapolate ?)
                 un(i1-is1,i2-is2,i3-is3,ex) = -(g - cm1*un(i1-is1,i2-is2,i3-is3,ex) )/cm1 
              end if
            else 
              ! --- adjacent face is NOT another ABC ---
              ! Extrapolate -- do this for now
              extrapOrder=5
              if( bc1.ne.symmetryBoundaryCondition )then
                  ksv(0)=0
                  ksv(1)=0
                  ksv(2)=0
                  ksv(0)=1-2*side1
                  ks1=ksv(0)
                  ks2=ksv(1)
                  ks3=ksv(2)
                  if( extrapOrder.eq.3 )then
                    un(i1-ks1,i2-ks2,i3-ks3,ex)=(3.*un(i1,i2,i3,ex)-3.*un(i1+ks1,i2+ks2,i3+ks3,ex)+un(i1+2*ks1,i2+2*ks2,i3+2*ks3,ex))
                    ! un(i1-ks1,i2-ks2,i3-ks3,ey)=extrap3(un,i1,i2,i3,ey,ks1,ks2,ks3)
                    ! un(i1-ks1,i2-ks2,i3-ks3,hz)=extrap3(un,i1,i2,i3,hz,ks1,ks2,ks3)
                  else if( extrapOrder.eq.5 )then
                    un(i1-ks1,i2-ks2,i3-ks3,ex)=(5.*un(i1,i2,i3,ex)-10.*un(i1+ks1,i2+ks2,i3+ks3,ex)+10.*un(i1+2*ks1,i2+2*ks2,i3+2*ks3,ex)-5.*un(i1+3*ks1,i2+3*ks2,i3+3*ks3,ex)+un(i1+4*ks1,i2+4*ks2,i3+4*ks3,ex))
                    ! un(i1-ks1,i2-ks2,i3-ks3,ey)=extrap5(un,i1,i2,i3,ey,ks1,ks2,ks3)
                    ! un(i1-ks1,i2-ks2,i3-ks3,hz)=extrap5(un,i1,i2,i3,hz,ks1,ks2,ks3)
                  else
                    stop 1782
                  end if
              end if
              if( bc2.ne.symmetryBoundaryCondition )then
                  ksv(0)=0
                  ksv(1)=0
                  ksv(2)=0
                  ksv(1)=1-2*side2
                  ks1=ksv(0)
                  ks2=ksv(1)
                  ks3=ksv(2)
                  if( extrapOrder.eq.3 )then
                    un(i1-ks1,i2-ks2,i3-ks3,ex)=(3.*un(i1,i2,i3,ex)-3.*un(i1+ks1,i2+ks2,i3+ks3,ex)+un(i1+2*ks1,i2+2*ks2,i3+2*ks3,ex))
                    ! un(i1-ks1,i2-ks2,i3-ks3,ey)=extrap3(un,i1,i2,i3,ey,ks1,ks2,ks3)
                    ! un(i1-ks1,i2-ks2,i3-ks3,hz)=extrap3(un,i1,i2,i3,hz,ks1,ks2,ks3)
                  else if( extrapOrder.eq.5 )then
                    un(i1-ks1,i2-ks2,i3-ks3,ex)=(5.*un(i1,i2,i3,ex)-10.*un(i1+ks1,i2+ks2,i3+ks3,ex)+10.*un(i1+2*ks1,i2+2*ks2,i3+2*ks3,ex)-5.*un(i1+3*ks1,i2+3*ks2,i3+3*ks3,ex)+un(i1+4*ks1,i2+4*ks2,i3+4*ks3,ex))
                    ! un(i1-ks1,i2-ks2,i3-ks3,ey)=extrap5(un,i1,i2,i3,ey,ks1,ks2,ks3)
                    ! un(i1-ks1,i2-ks2,i3-ks3,hz)=extrap5(un,i1,i2,i3,hz,ks1,ks2,ks3)
                  else
                    stop 1782
                  end if
              end if 
            end if ! end if adjacentFaceIsABC
            ! --- Now extrapolate the 3 ghost points "G" adjacent to the corner ---
            j1=i1-is1; j2=i2-is2; j3=0; 
            un(j1,j2,j3,ex)=(5.*un(j1+is1,j2+is2,j3+is3,ex)-10.*un(j1+is1+is1,j2+is2+is2,j3+is3+is3,ex)+10.*un(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,ex)-5.*un(j1+is1+3*is1,j2+is2+3*is2,j3+is3+3*is3,ex)+un(j1+is1+4*is1,j2+is2+4*is2,j3+is3+4*is3,ex))
            j1=i1-2*is1; j2=i2-is2; 
            un(j1,j2,j3,ex)=(5.*un(j1+is1,j2+is2,j3+is3,ex)-10.*un(j1+is1+is1,j2+is2+is2,j3+is3+is3,ex)+10.*un(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,ex)-5.*un(j1+is1+3*is1,j2+is2+3*is2,j3+is3+3*is3,ex)+un(j1+is1+4*is1,j2+is2+4*is2,j3+is3+4*is3,ex))
            j1=i1-is1; j2=i2-2*is2; 
            un(j1,j2,j3,ex)=(5.*un(j1+is1,j2+is2,j3+is3,ex)-10.*un(j1+is1+is1,j2+is2+is2,j3+is3+is3,ex)+10.*un(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,ex)-5.*un(j1+is1+3*is1,j2+is2+3*is2,j3+is3+3*is3,ex)+un(j1+is1+4*is1,j2+is2+4*is2,j3+is3+4*is3,ex))
           end if ! mask 
          end if ! if one face is ABC
        end do          
        end do          
      else 
        ! ***** 3D *****
        write(*,*) 'abcWave: corners and edges FINISH ME - 3D'
        stop 3434
      end if  

  else
    write(*,*) 'abcWave: assign edges and corners, finish me for orderOfAccuracy=',orderOfAccuracy
    stop 6666
  end if



  ! -------------------------------------------------------------------------
  ! ------------------Loop over Sides----------------------------------------
  ! -------------------------------------------------------------------------
   extra1a=extra
   extra1b=extra
   extra2a=extra
   extra2b=extra
   if( nd.eq.3 )then
     extra3a=extra
     extra3b=extra
   else
     extra3a=0
     extra3b=0
   end if
   if( boundaryCondition(0,0).lt.0 )then
     extra1a=max(0,extra1a) ! over-ride extra=-1 : assign ends in periodic directions
     extra1b=extra1a
   else
     if( boundaryCondition(0,0).eq.0 )then
       extra1a=numberOfGhostPoints  ! include interpolation points since we assign ghost points outside these
     end if
     if( boundaryCondition(1,0).eq.0 )then
       extra1b=numberOfGhostPoints
     end if
   end if
   if( boundaryCondition(0,1).lt.0 )then
    extra2a=max(0,extra2a) ! over-ride extra=-1 : assign ends in periodic directions
    extra2b=extra2a
   else 
     if( boundaryCondition(0,1).eq.0 )then
       extra2a=numberOfGhostPoints
     end if
     if( boundaryCondition(1,1).eq.0 )then
       extra2b=numberOfGhostPoints
     end if
   end if
   if(  nd.eq.3 .and. boundaryCondition(0,2).lt.0 )then
    extra3a=max(0,extra3a) ! over-ride extra=-1 : assign ends in periodic directions
    extra3b=extra3a
   else 
     if( boundaryCondition(0,2).eq.0 )then
       extra3a=numberOfGhostPoints
     end if
     if( boundaryCondition(1,2).eq.0 )then
       extra3b=numberOfGhostPoints
     end if
   end if
   do axis=0,nd-1
   do side=0,1
     if( boundaryCondition(side,axis).eq.abcEM2 .or. boundaryCondition(side,axis).eq.absorbing )then
       ! write(*,'(" abc: grid,side,axis,bc=",3i2)') grid,side,axis,boundaryCondition(side,axis)
       n1a=gridIndexRange(0,0)-extra1a
       n1b=gridIndexRange(1,0)+extra1b
       n2a=gridIndexRange(0,1)-extra2a
       n2b=gridIndexRange(1,1)+extra2b
       n3a=gridIndexRange(0,2)-extra3a
       n3b=gridIndexRange(1,2)+extra3b
       if( axis.eq.0 )then
         n1a=gridIndexRange(side,axis)
         n1b=gridIndexRange(side,axis)
       else if( axis.eq.1 )then
         n2a=gridIndexRange(side,axis)
         n2b=gridIndexRange(side,axis)
       else
         n3a=gridIndexRange(side,axis)
         n3b=gridIndexRange(side,axis)
       end if
       is1=0
       is2=0
       is3=0
       if( axis.eq.0 )then
         is1=1-2*side
       else if( axis.eq.1 )then
         is2=1-2*side
       else if( axis.eq.2 )then
         is3=1-2*side
       else
         stop 5
       end if
       is = 1 - 2*side
       axisp1=mod(axis+1,nd)
       axisp2=mod(axis+2,nd)
       ! (js1,js2,js3) used to compute tangential derivatives
       js1=0
       js2=0
       js3=0
       if( axisp1.eq.0 )then
         js1=1-2*side
       else if( axisp1.eq.1 )then
         js2=1-2*side
       else if( axisp1.eq.2 )then
         js3=1-2*side
       else
         stop 5
       end if
       ! (ks1,ks2,ks3) used to compute second tangential derivative
       ks1=0
       ks2=0
       ks3=0
       if( axisp2.eq.0 )then
         ks1=1-2*side
       else if( axisp2.eq.1 )then
         ks2=1-2*side
       else if( axisp2.eq.2 )then
         ks3=1-2*side
       else
         stop 5
       end if


   if( gridType.eq.rectangular .and. orderOfAccuracy.eq.2 )then
    ! ***********************************************
    ! ************rectangular grid*******************
    ! ***********************************************

   
    ! write(*,'(" Apply abcEM2 rectangular: grid,side,axis=",3i3," dt,c=",2e12.3," is1,is2=",2i2)') grid,side,axis,dt,c,is1,is2
    if( nd.eq.2 )then
     if( axis.eq.0 )then

      checkResiduals=.false.
      if( checkResiduals )then
        do i3=n3a,n3b
        do i2=n2a,n2b
        do i1=n1a,n1b
        if( mask(i1,i2,i3).gt.0 )then
              ! macro: 
           ! getForcingEM2(X,2,tm,is1,is2,forcex)
            if( forcingOption.eq.twilightZoneForcing )then
              ! Test: set to exact solution at time t:
              ! x=xy(i1-is1,i2,i3,0)
              ! y=xy(i1-is1,i2,i3,1)
              ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
              ! un(i1-is1,i2,i3,ey)=eyTrue
              ! add TZ forcing *wdh* Sept 17, 2016
              ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
              x=xy(i1,i2,i3,0)
              y=xy(i1,i2,i3,1)
                 ! Values for forcep(ex) are currently needed at corners:
                   call ogDeriv(ep, 1,1,0,0,x,y,z,tp,ex,utx)
                   call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                   call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                 forcep(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                 ! OGDERIV(1,1,0,0,x,y,z,tp,ey,utx)
                 ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                 ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                 ! forcep(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                 ! OGDERIV(1,1,0,0,x,y,z,tp,hz,utx)
                 ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
                 ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
                 ! forcep(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
              ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
            end if
            if( forcingOption.eq.twilightZoneForcing )then
              ! Test: set to exact solution at time t:
              ! x=xy(i1-is1,i2,i3,0)
              ! y=xy(i1-is1,i2,i3,1)
              ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
              ! un(i1-is1,i2,i3,ey)=eyTrue
              ! add TZ forcing *wdh* Sept 17, 2016
              ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
              x=xy(i1,i2,i3,0)
              y=xy(i1,i2,i3,1)
                 ! Values for forcef(ex) are currently needed at corners:
                   call ogDeriv(ep, 1,1,0,0,x,y,z,t,ex,utx)
                   call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                   call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                 forcef(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                 ! OGDERIV(1,1,0,0,x,y,z,t,ey,utx)
                 ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
                 ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
                 ! forcef(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                 ! OGDERIV(1,1,0,0,x,y,z,t,hz,utx)
                 ! OGDERIV(0,2,0,0,x,y,z,t,hz,uxx)
                 ! OGDERIV(0,0,2,0,x,y,z,t,hz,uyy)
                 ! forcef(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
              ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
            end if
           ! do idir=0,2
           do idir=0,0 ! ex only 
             forcex(idir)=.5*(forcep(idir)+forcef(idir))
           end do

         ! check residual on input (use to check implicit solve)
         res = un(i1-is1,i2-is2,i3-is3,ex) - ((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dxa*dt)*(c1abcem2*uxx22r(i1,i2,i3,ex)+c2abcem2*uyy22r(i1,i2,i3,ex)+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1+is1,i2,i3,ex))/dxa**2+c2abcem2*(un(i1,i2-1,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2+1,i3,ex))/dx(1)**2+2.*forcex(ex)))/(1.+c1abcem2*dt/dxa))

         ! res -= -is*(unx-ucx)/dt  + (.5*ca)*( unxx + ucxx ) + (.25*ca)*( unyy + ucyy );
         r1 = is1*(unx22r(i1,i2,i3,ex)-ux22r(i1,i2,i3,ex))/(dt)- .5*( (cEM2)*unxx22r(i1,i2,i3,ex) + (.5*cEM2)*unyy22r(i1,i2,i3,ex) +(cEM2)* uxx22r(i1,i2,i3,ex) + (.5*cEM2)* uyy22r(i1,i2,i3,ex) ) - forcex(ex) 

         ! res2 = -is1*(unx-ucx)/dt  + (.5*cEM2)*( unxx + ucxx ) + (.25*cEM2)*( unyy + ucyy );         
         write(*,'(">EM2: t=",1pe10.2," cEM2=",1pe12.4," i1,i2,forcex,un,res,r1=",2i3,1pe10.2,1pe10.2,1pe10.2,1pe10.2)') t,cEM2,i1,i2,forcex(ex),un(i1-is1,i2-is2,i3-is3,ex),res,r1
        end if
        end do
        end do
        end do
      end if

      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
   
       ! macro: 
         ! getForcingEM2(X,2,tm,is1,is2,forcex)
          if( forcingOption.eq.twilightZoneForcing )then
            ! Test: set to exact solution at time t:
            ! x=xy(i1-is1,i2,i3,0)
            ! y=xy(i1-is1,i2,i3,1)
            ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
            ! un(i1-is1,i2,i3,ey)=eyTrue
            ! add TZ forcing *wdh* Sept 17, 2016
            ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
            x=xy(i1,i2,i3,0)
            y=xy(i1,i2,i3,1)
               ! Values for forcep(ex) are currently needed at corners:
                 call ogDeriv(ep, 1,1,0,0,x,y,z,tp,ex,utx)
                 call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                 call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
               forcep(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
               ! OGDERIV(1,1,0,0,x,y,z,tp,ey,utx)
               ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
               ! forcep(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
               ! OGDERIV(1,1,0,0,x,y,z,tp,hz,utx)
               ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
               ! forcep(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
            ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
          end if
          if( forcingOption.eq.twilightZoneForcing )then
            ! Test: set to exact solution at time t:
            ! x=xy(i1-is1,i2,i3,0)
            ! y=xy(i1-is1,i2,i3,1)
            ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
            ! un(i1-is1,i2,i3,ey)=eyTrue
            ! add TZ forcing *wdh* Sept 17, 2016
            ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
            x=xy(i1,i2,i3,0)
            y=xy(i1,i2,i3,1)
               ! Values for forcef(ex) are currently needed at corners:
                 call ogDeriv(ep, 1,1,0,0,x,y,z,t,ex,utx)
                 call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                 call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
               forcef(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
               ! OGDERIV(1,1,0,0,x,y,z,t,ey,utx)
               ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
               ! forcef(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
               ! OGDERIV(1,1,0,0,x,y,z,t,hz,utx)
               ! OGDERIV(0,2,0,0,x,y,z,t,hz,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,t,hz,uyy)
               ! forcef(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
            ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
          end if
         ! do idir=0,2
         do idir=0,0 ! ex only 
           forcex(idir)=.5*(forcep(idir)+forcef(idir))
         end do

       un(i1-is1,i2-is2,i3-is3,ex)=((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dxa*dt)*(c1abcem2*uxx22r(i1,i2,i3,ex)+c2abcem2*uyy22r(i1,i2,i3,ex)+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1+is1,i2,i3,ex))/dxa**2+c2abcem2*(un(i1,i2-1,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2+1,i3,ex))/dx(1)**2+2.*forcex(ex)))/(1.+c1abcem2*dt/dxa))

      end if
      end do
      end do
      end do


     else if( axis.eq.1 )then

      if( checkResiduals )then
        do i3=n3a,n3b
        do i2=n2a,n2b
        do i1=n1a,n1b
        if( mask(i1,i2,i3).gt.0 )then
          ! macro: 
            ! getForcingEM2(Y,2,tm,is1,is2,forcey)
             if( forcingOption.eq.twilightZoneForcing )then
               ! Test: set to exact solution at time t:
               ! x=xy(i1-is1,i2,i3,0)
               ! y=xy(i1-is1,i2,i3,1)
               ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
               ! un(i1-is1,i2,i3,ey)=eyTrue
               ! add TZ forcing *wdh* Sept 17, 2016
               ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
               x=xy(i1,i2,i3,0)
               y=xy(i1,i2,i3,1)
                    call ogDeriv(ep, 1,0,1,0,x,y,z,tp,ex,uty)
                    call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                    call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                  forcep(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                  ! ! Values for forcep(ey) are currently needed at corners:
                  ! OGDERIV(1,0,1,0,x,y,z,tp,ey,uty)
                  ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                  ! forcep(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                  ! OGDERIV(1,0,1,0,x,y,z,tp,hz,uty)
                  ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
                  ! forcep(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
               ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
             end if
             if( forcingOption.eq.twilightZoneForcing )then
               ! Test: set to exact solution at time t:
               ! x=xy(i1-is1,i2,i3,0)
               ! y=xy(i1-is1,i2,i3,1)
               ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
               ! un(i1-is1,i2,i3,ey)=eyTrue
               ! add TZ forcing *wdh* Sept 17, 2016
               ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
               x=xy(i1,i2,i3,0)
               y=xy(i1,i2,i3,1)
                    call ogDeriv(ep, 1,0,1,0,x,y,z,t,ex,uty)
                    call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                    call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                  forcef(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                  ! ! Values for forcef(ey) are currently needed at corners:
                  ! OGDERIV(1,0,1,0,x,y,z,t,ey,uty)
                  ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
                  ! forcef(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                  ! OGDERIV(1,0,1,0,x,y,z,t,hz,uty)
                  ! OGDERIV(0,2,0,0,x,y,z,t,hz,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,t,hz,uyy)
                  ! forcef(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
               ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
             end if
            ! do idir=0,2
            do idir=0,0 ! ex only 
              forcey(idir)=.5*(forcep(idir)+forcef(idir))
            end do

         ! check residual on input (use to check implicit solve)
         res = un(i1-is1,i2-is2,i3-is3,ex) - ((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dya*dt)*(c1abcem2*uyy22r(i1,i2,i3,ex)+c2abcem2*uxx22r(i1,i2,i3,ex)+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1,i2+is2,i3,ex))/dya**2+c2abcem2*(un(i1-1,i2,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1+1,i2,i3,ex))/dx(0)**2+2.*forcey(ex)))/(1.+c1abcem2*dt/dya))
         write(*,'(">EM2: t=",1pe10.2," cEM2=",1pe12.4," i1,i2,forcey,residual=",2i3,1pe10.2,1pe10.2)') t,cEM2,i1,i2,forcey(ex),res

        end if
        end do
        end do
        end do
      end if

      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then

       ! macro: 
         ! getForcingEM2(Y,2,tm,is1,is2,forcey)
          if( forcingOption.eq.twilightZoneForcing )then
            ! Test: set to exact solution at time t:
            ! x=xy(i1-is1,i2,i3,0)
            ! y=xy(i1-is1,i2,i3,1)
            ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
            ! un(i1-is1,i2,i3,ey)=eyTrue
            ! add TZ forcing *wdh* Sept 17, 2016
            ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
            x=xy(i1,i2,i3,0)
            y=xy(i1,i2,i3,1)
                 call ogDeriv(ep, 1,0,1,0,x,y,z,tp,ex,uty)
                 call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                 call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
               forcep(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
               ! ! Values for forcep(ey) are currently needed at corners:
               ! OGDERIV(1,0,1,0,x,y,z,tp,ey,uty)
               ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
               ! forcep(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
               ! OGDERIV(1,0,1,0,x,y,z,tp,hz,uty)
               ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
               ! forcep(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
            ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
          end if
          if( forcingOption.eq.twilightZoneForcing )then
            ! Test: set to exact solution at time t:
            ! x=xy(i1-is1,i2,i3,0)
            ! y=xy(i1-is1,i2,i3,1)
            ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
            ! un(i1-is1,i2,i3,ey)=eyTrue
            ! add TZ forcing *wdh* Sept 17, 2016
            ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
            x=xy(i1,i2,i3,0)
            y=xy(i1,i2,i3,1)
                 call ogDeriv(ep, 1,0,1,0,x,y,z,t,ex,uty)
                 call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                 call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
               forcef(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
               ! ! Values for forcef(ey) are currently needed at corners:
               ! OGDERIV(1,0,1,0,x,y,z,t,ey,uty)
               ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
               ! forcef(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
               ! OGDERIV(1,0,1,0,x,y,z,t,hz,uty)
               ! OGDERIV(0,2,0,0,x,y,z,t,hz,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,t,hz,uyy)
               ! forcef(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
            ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
          end if
         ! do idir=0,2
         do idir=0,0 ! ex only 
           forcey(idir)=.5*(forcep(idir)+forcef(idir))
         end do

       un(i1-is1,i2-is2,i3-is3,ex)=((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dya*dt)*(c1abcem2*uyy22r(i1,i2,i3,ex)+c2abcem2*uxx22r(i1,i2,i3,ex)+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1,i2+is2,i3,ex))/dya**2+c2abcem2*(un(i1-1,i2,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1+1,i2,i3,ex))/dx(0)**2+2.*forcey(ex)))/(1.+c1abcem2*dt/dya))
 
      end if
      end do
      end do
      end do
     else
      stop 9467
     end if

    else ! ***** 3D *****

     ! -------- Cartesian Grid 3D --------
     if( axis.eq.0 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then

         ! getForcingEM2(X,3,tm,is1,is2,forcex)
          if( forcingOption.eq.twilightZoneForcing )then
            ! Test: set to exact solution at time t:
            ! x=xy(i1-is1,i2,i3,0)
            ! y=xy(i1-is1,i2,i3,1)
            ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
            ! un(i1-is1,i2,i3,ey)=eyTrue
            ! add TZ forcing *wdh* Sept 17, 2016
            ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
            x=xy(i1,i2,i3,0)
            y=xy(i1,i2,i3,1)
              ! ------ Cartesian Grid 3d forcing ----------
              z=xy(i1,i2,i3,2)
               ! Values for forcep(ex) are currently needed at corners:
                 call ogDeriv(ep, 1,1,0,0,x,y,z,tp,ex,utx)
                 call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                 call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                 call ogDeriv(ep, 0,0,0,2,x,y,z,tp,ex,uzz)
               forcep(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) ) 
               ! OGDERIV(1,1,0,0,x,y,z,tp,ey,utx)
               ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,tp,ey,uzz)
               ! forcep(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) ) 
               ! OGDERIV(1,1,0,0,x,y,z,tp,ez,utx)
               ! OGDERIV(0,2,0,0,x,y,z,tp,ez,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,tp,ez,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,tp,ez,uzz)
               ! forcep(ez) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) )
            ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
          end if
          if( forcingOption.eq.twilightZoneForcing )then
            ! Test: set to exact solution at time t:
            ! x=xy(i1-is1,i2,i3,0)
            ! y=xy(i1-is1,i2,i3,1)
            ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
            ! un(i1-is1,i2,i3,ey)=eyTrue
            ! add TZ forcing *wdh* Sept 17, 2016
            ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
            x=xy(i1,i2,i3,0)
            y=xy(i1,i2,i3,1)
              ! ------ Cartesian Grid 3d forcing ----------
              z=xy(i1,i2,i3,2)
               ! Values for forcef(ex) are currently needed at corners:
                 call ogDeriv(ep, 1,1,0,0,x,y,z,t,ex,utx)
                 call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                 call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                 call ogDeriv(ep, 0,0,0,2,x,y,z,t,ex,uzz)
               forcef(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) ) 
               ! OGDERIV(1,1,0,0,x,y,z,t,ey,utx)
               ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,t,ey,uzz)
               ! forcef(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) ) 
               ! OGDERIV(1,1,0,0,x,y,z,t,ez,utx)
               ! OGDERIV(0,2,0,0,x,y,z,t,ez,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,t,ez,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,t,ez,uzz)
               ! forcef(ez) = is1*utx - ( c1abcem2*uxx + c2abcem2*(uyy+uzz) )
            ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
          end if
         ! do idir=0,2
         do idir=0,0 ! ex only 
           forcex(idir)=.5*(forcep(idir)+forcef(idir))
         end do

       un(i1-is1,i2-is2,i3-is3,ex)=((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dxa*dt)*(c1abcem2*uxx23r(i1,i2,i3,ex)+c2abcem2*(uyy23r(i1,i2,i3,ex)+uzz23r(i1,i2,i3,ex))+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1+is1,i2,i3,ex))/dxa**2+c2abcem2*(un(i1,i2-1,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2+1,i3,ex))/dx(1)**2+c2abcem2*(un(i1,i2,i3-1,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2,i3+1,ex))/dx(2)**2+2.*forcex(ex)))/(1.+c1abcem2*dt/dxa))

      end if
      end do
      end do
      end do
     else if( axis.eq.1 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then

         ! getForcingEM2(Y,3,tm,is1,is2,forcey)
          if( forcingOption.eq.twilightZoneForcing )then
            ! Test: set to exact solution at time t:
            ! x=xy(i1-is1,i2,i3,0)
            ! y=xy(i1-is1,i2,i3,1)
            ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
            ! un(i1-is1,i2,i3,ey)=eyTrue
            ! add TZ forcing *wdh* Sept 17, 2016
            ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
            x=xy(i1,i2,i3,0)
            y=xy(i1,i2,i3,1)
              ! ------ Cartesian Grid 3d forcing ----------
              z=xy(i1,i2,i3,2)
                 call ogDeriv(ep, 1,0,1,0,x,y,z,tp,ex,uty)
                 call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                 call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                 call ogDeriv(ep, 0,0,0,2,x,y,z,tp,ex,uzz)
               forcep(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) ) 
               ! OGDERIV(1,0,1,0,x,y,z,tp,ey,uty)
               ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,tp,ey,uzz)
               ! forcep(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) ) 
               ! OGDERIV(1,0,1,0,x,y,z,tp,ez,uty)
               ! OGDERIV(0,2,0,0,x,y,z,tp,ez,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,tp,ez,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,tp,ez,uzz)
               ! forcep(ez) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) )
            ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
          end if
          if( forcingOption.eq.twilightZoneForcing )then
            ! Test: set to exact solution at time t:
            ! x=xy(i1-is1,i2,i3,0)
            ! y=xy(i1-is1,i2,i3,1)
            ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
            ! un(i1-is1,i2,i3,ey)=eyTrue
            ! add TZ forcing *wdh* Sept 17, 2016
            ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
            x=xy(i1,i2,i3,0)
            y=xy(i1,i2,i3,1)
              ! ------ Cartesian Grid 3d forcing ----------
              z=xy(i1,i2,i3,2)
                 call ogDeriv(ep, 1,0,1,0,x,y,z,t,ex,uty)
                 call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                 call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                 call ogDeriv(ep, 0,0,0,2,x,y,z,t,ex,uzz)
               forcef(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) ) 
               ! OGDERIV(1,0,1,0,x,y,z,t,ey,uty)
               ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,t,ey,uzz)
               ! forcef(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) ) 
               ! OGDERIV(1,0,1,0,x,y,z,t,ez,uty)
               ! OGDERIV(0,2,0,0,x,y,z,t,ez,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,t,ez,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,t,ez,uzz)
               ! forcef(ez) = is2*uty - ( c1abcem2*uyy + c2abcem2*(uxx+uzz) )
            ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
          end if
         ! do idir=0,2
         do idir=0,0 ! ex only 
           forcey(idir)=.5*(forcep(idir)+forcef(idir))
         end do

       un(i1-is1,i2-is2,i3-is3,ex)=((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dya*dt)*(c1abcem2*uyy23r(i1,i2,i3,ex)+c2abcem2*(uxx23r(i1,i2,i3,ex)+uzz23r(i1,i2,i3,ex))+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1,i2+is2,i3,ex))/dya**2+c2abcem2*(un(i1-1,i2,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1+1,i2,i3,ex))/dx(0)**2+c2abcem2*(un(i1,i2,i3-1,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2,i3+1,ex))/dx(2)**2+2.*forcey(ex)))/(1.+c1abcem2*dt/dya))

      end if
      end do
      end do
      end do
     else if( axis.eq.2 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then

         ! getForcingEM2(Z,3,tm,is1,is2,forcez)
          if( forcingOption.eq.twilightZoneForcing )then
            ! Test: set to exact solution at time t:
            ! x=xy(i1-is1,i2,i3,0)
            ! y=xy(i1-is1,i2,i3,1)
            ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
            ! un(i1-is1,i2,i3,ey)=eyTrue
            ! add TZ forcing *wdh* Sept 17, 2016
            ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
            x=xy(i1,i2,i3,0)
            y=xy(i1,i2,i3,1)
              ! ------ Cartesian Grid 3d forcing ----------
              z=xy(i1,i2,i3,2)
                 call ogDeriv(ep, 1,0,0,1,x,y,z,tp,ex,utz)
                 call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                 call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                 call ogDeriv(ep, 0,0,0,2,x,y,z,tp,ex,uzz)
               forcep(ex) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) ) 
               ! OGDERIV(1,0,0,1,x,y,z,tp,ey,utz)
               ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,tp,ey,uzz)
               ! forcep(ey) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) ) 
               ! OGDERIV(1,0,0,1,x,y,z,tp,ez,utz)
               ! OGDERIV(0,2,0,0,x,y,z,tp,ez,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,tp,ez,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,tp,ez,uzz)
               ! forcep(ez) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) )
            ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
          end if
          if( forcingOption.eq.twilightZoneForcing )then
            ! Test: set to exact solution at time t:
            ! x=xy(i1-is1,i2,i3,0)
            ! y=xy(i1-is1,i2,i3,1)
            ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
            ! un(i1-is1,i2,i3,ey)=eyTrue
            ! add TZ forcing *wdh* Sept 17, 2016
            ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
            x=xy(i1,i2,i3,0)
            y=xy(i1,i2,i3,1)
              ! ------ Cartesian Grid 3d forcing ----------
              z=xy(i1,i2,i3,2)
                 call ogDeriv(ep, 1,0,0,1,x,y,z,t,ex,utz)
                 call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                 call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                 call ogDeriv(ep, 0,0,0,2,x,y,z,t,ex,uzz)
               forcef(ex) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) ) 
               ! OGDERIV(1,0,0,1,x,y,z,t,ey,utz)
               ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,t,ey,uzz)
               ! forcef(ey) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) ) 
               ! OGDERIV(1,0,0,1,x,y,z,t,ez,utz)
               ! OGDERIV(0,2,0,0,x,y,z,t,ez,uxx)
               ! OGDERIV(0,0,2,0,x,y,z,t,ez,uyy)
               ! OGDERIV(0,0,0,2,x,y,z,t,ez,uzz)
               ! forcef(ez) = is3*utz - ( c1abcem2*uzz + c2abcem2*(uxx+uyy) )
            ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
          end if
         ! do idir=0,2
         do idir=0,0 ! ex only 
           forcez(idir)=.5*(forcep(idir)+forcef(idir))
         end do

       un(i1-is1,i2-is2,i3-is3,ex)=((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dza*dt)*(c1abcem2*uzz23r(i1,i2,i3,ex)+c2abcem2*(uxx23r(i1,i2,i3,ex)+uyy23r(i1,i2,i3,ex))+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1,i2,i3+is3,ex))/dza**2+c2abcem2*(un(i1-1,i2,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1+1,i2,i3,ex))/dx(0)**2+c2abcem2*(un(i1,i2-1,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2+1,i3,ex))/dx(1)**2+2.*forcez(ex)))/(1.+c1abcem2*dt/dza))

      end if
      end do
      end do
      end do
     else
      stop 9468
     end if

    end if

    ! added: July 3, 2023 
    if( numGhost>1 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
        do ghost=2,numGhost
          j1= i1-ghost*is1; j2=i2-ghost*is2; j3=i3-ghost*is3;
          un(j1,j2,j3,ex)=(3.*un(j1+is1,j2+is2,j3+is3,ex)-3.*un(j1+is1+is1,j2+is2+is2,j3+is3+is3,ex)+un(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,ex))
        end do
      end if
      end do
      end do
      end do
    end if

 
   else if( gridType.eq.rectangular .and. orderOfAccuracy.eq.4 )then

    ! ***********************************************
    ! ************rectangular grid*******************
    ! ************ fourth-order   *******************
    ! ***********************************************

   
    ! write(*,'(" Apply abcEM2 CARTESIAN order=4: grid,side,axis=",3i3," dt,c=",2e12.3)') grid,side,axis,dt,c

    if( nd.eq.2 )then

     if( .true. )then

      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
        ! Solve: 
        !   f1 = is* u_xt - c( u_xx + .5*u_yy )
        !   f2 = D^5 U_{-2} = 0 
        !
        !   f1 = a11*U(-1) + a12*U(-2) + ...

        ! Evaluate the equations with current ghost 
        ! There two coupled equations we need to solve for ghost points A,B below
        !                   +
        !                   |
        !             B--A--+---+---+
        !                   |
        !                   +
        !  f(u)  = [ f(u_old) - A (u_old) ] + A u = 0 
        !
        !  [ a11 a12 ][ uA ] = [ a11 a12 ][ uA_old ] - [ f1(u_old) ]
        !  [ a21 a22 ][ uB ]   [ a21 a22 ][ uB_old ]   [ f2(u_old) ]

        ! first evaluate residuals in equations given current (wrong) values at A, B

        if( axis==0 )then
            ! getForcingEM2(X,2,tm,is1,is2,forcex)
             if( forcingOption.eq.twilightZoneForcing )then
               ! Test: set to exact solution at time t:
               ! x=xy(i1-is1,i2,i3,0)
               ! y=xy(i1-is1,i2,i3,1)
               ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
               ! un(i1-is1,i2,i3,ey)=eyTrue
               ! add TZ forcing *wdh* Sept 17, 2016
               ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
               x=xy(i1,i2,i3,0)
               y=xy(i1,i2,i3,1)
                  ! Values for forcep(ex) are currently needed at corners:
                    call ogDeriv(ep, 1,1,0,0,x,y,z,tp,ex,utx)
                    call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                    call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                  forcep(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                  ! OGDERIV(1,1,0,0,x,y,z,tp,ey,utx)
                  ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                  ! forcep(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                  ! OGDERIV(1,1,0,0,x,y,z,tp,hz,utx)
                  ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
                  ! forcep(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
               ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
             end if
             if( forcingOption.eq.twilightZoneForcing )then
               ! Test: set to exact solution at time t:
               ! x=xy(i1-is1,i2,i3,0)
               ! y=xy(i1-is1,i2,i3,1)
               ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
               ! un(i1-is1,i2,i3,ey)=eyTrue
               ! add TZ forcing *wdh* Sept 17, 2016
               ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
               x=xy(i1,i2,i3,0)
               y=xy(i1,i2,i3,1)
                  ! Values for forcef(ex) are currently needed at corners:
                    call ogDeriv(ep, 1,1,0,0,x,y,z,t,ex,utx)
                    call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                    call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                  forcef(ex) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                  ! OGDERIV(1,1,0,0,x,y,z,t,ey,utx)
                  ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
                  ! forcef(ey) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy ) 
                  ! OGDERIV(1,1,0,0,x,y,z,t,hz,utx)
                  ! OGDERIV(0,2,0,0,x,y,z,t,hz,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,t,hz,uyy)
                  ! forcef(hz) = is1*utx - ( c1abcem2*uxx + c2abcem2*uyy )
               ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
             end if
            ! do idir=0,2
            do idir=0,0 ! ex only 
              forcex(idir)=.5*(forcep(idir)+forcef(idir))
            end do
          r1 = is1*(unx42r(i1,i2,i3,ex)-ux42r(i1,i2,i3,ex))/(dt)- .5*( c1abcem2*unxx42r(i1,i2,i3,ex) + c2abcem2*unyy42r(i1,i2,i3,ex) + c1abcem2* uxx42r(i1,i2,i3,ex) + c2abcem2* uyy42r(i1,i2,i3,ex) ) - forcex(ex) 
        else
            ! getForcingEM2(Y,2,tm,is1,is2,forcey)
             if( forcingOption.eq.twilightZoneForcing )then
               ! Test: set to exact solution at time t:
               ! x=xy(i1-is1,i2,i3,0)
               ! y=xy(i1-is1,i2,i3,1)
               ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
               ! un(i1-is1,i2,i3,ey)=eyTrue
               ! add TZ forcing *wdh* Sept 17, 2016
               ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
               x=xy(i1,i2,i3,0)
               y=xy(i1,i2,i3,1)
                    call ogDeriv(ep, 1,0,1,0,x,y,z,tp,ex,uty)
                    call ogDeriv(ep, 0,2,0,0,x,y,z,tp,ex,uxx)
                    call ogDeriv(ep, 0,0,2,0,x,y,z,tp,ex,uyy)
                  forcep(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                  ! ! Values for forcep(ey) are currently needed at corners:
                  ! OGDERIV(1,0,1,0,x,y,z,tp,ey,uty)
                  ! OGDERIV(0,2,0,0,x,y,z,tp,ey,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,tp,ey,uyy)
                  ! forcep(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                  ! OGDERIV(1,0,1,0,x,y,z,tp,hz,uty)
                  ! OGDERIV(0,2,0,0,x,y,z,tp,hz,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,tp,hz,uyy)
                  ! forcep(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
               ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
             end if
             if( forcingOption.eq.twilightZoneForcing )then
               ! Test: set to exact solution at time t:
               ! x=xy(i1-is1,i2,i3,0)
               ! y=xy(i1-is1,i2,i3,1)
               ! OGDERIV(0,0,0,0,x,y,z,t,ey,eyTrue)
               ! un(i1-is1,i2,i3,ey)=eyTrue
               ! add TZ forcing *wdh* Sept 17, 2016
               ! OGDERIV(ntd,nxd,nyd,nzd,x,y,z,t,n,ud)
               x=xy(i1,i2,i3,0)
               y=xy(i1,i2,i3,1)
                    call ogDeriv(ep, 1,0,1,0,x,y,z,t,ex,uty)
                    call ogDeriv(ep, 0,2,0,0,x,y,z,t,ex,uxx)
                    call ogDeriv(ep, 0,0,2,0,x,y,z,t,ex,uyy)
                  forcef(ex) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                  ! ! Values for forcef(ey) are currently needed at corners:
                  ! OGDERIV(1,0,1,0,x,y,z,t,ey,uty)
                  ! OGDERIV(0,2,0,0,x,y,z,t,ey,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,t,ey,uyy)
                  ! forcef(ey) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx ) 
                  ! OGDERIV(1,0,1,0,x,y,z,t,hz,uty)
                  ! OGDERIV(0,2,0,0,x,y,z,t,hz,uxx)
                  ! OGDERIV(0,0,2,0,x,y,z,t,hz,uyy)
                  ! forcef(hz) = is2*uty - ( c1abcem2*uyy + c2abcem2*uxx )
               ! write(*,'(" Apply abcEM2: add TZ forcing t,dt,utx,uxx,uyy=",5e10.3)') t,dt,utx,uxx,uyy
             end if
            ! do idir=0,2
            do idir=0,0 ! ex only 
              forcey(idir)=.5*(forcep(idir)+forcef(idir))
            end do
          r1 = is2*(uny42r(i1,i2,i3,ex)-uy42r(i1,i2,i3,ex))/(dt)- .5*( c1abcem2*unyy42r(i1,i2,i3,ex) + c2abcem2*unxx42r(i1,i2,i3,ex) + c1abcem2* uyy42r(i1,i2,i3,ex) + c2abcem2* uxx42r(i1,i2,i3,ex) ) - forcey(ex)           
        end if

        r2 =    un(i1-2*is1,i2-2*is2,i3,ex) -5.*un(i1-1*is1,i2-1*is2,i3,ex) +10.*un(i1+0*is1,i2+0*is2,i3,ex) -10.*un(i1+1*is1,i2+1*is2,i3,ex) +5.*un(i1+2*is1,i2+2*is2,i3,ex) -un(i1+3*is1,i2+3*is2,i3,ex)


        a11 = -8./(12.*dt*dx(axis))  - .5*c1abcem2*(16.)/(12.*dx(axis)**2) ! coeff of uA in r1
        a12 =  1./(12.*dt*dx(axis))  - .5*c1abcem2*(-1.)/(12.*dx(axis)**2) ! coeff of uB in r2
        a21 = -5.                                                    ! coeff of uA in r1
        a22 =  1.                                                    ! coeff of uB in r2

        det = a11*a22-a21*a12

        uA = un(i1-1*is1,i2-1*is2,i3,ex)
        uB = un(i1-2*is1,i2-2*is2,i3,ex)
        f1 = a11*uA + a12*uB - r1 
        f2 = a21*uA + a22*uB - r2 

        ! Solve for A, B
        un(i1-1*is1,i2-1*is2,i3,ex) = ( f1*a22 - f2*a12)/det
        un(i1-2*is1,i2-2*is2,i3,ex) = (-f1*a21 + f2*a11)/det         

      end if
      end do
      end do
      end do

     else if( axis.eq.0 )then

        ! ----- OLD WAY -----

        ! Version 1: 2nd-order
        do i3=n3a,n3b
        do i2=n2a,n2b
        do i1=n1a,n1b
        if( mask(i1,i2,i3).gt.0 )then
     
         un(i1-is1,i2-is2,i3-is3,ex)=((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dxa*dt)*(c1abcem2*uxx22r(i1,i2,i3,ex)+c2abcem2*uyy22r(i1,i2,i3,ex)+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1+is1,i2,i3,ex))/dxa**2+c2abcem2*(un(i1,i2-1,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2+1,i3,ex))/dx(1)**2+2.*forcex(ex)))/(1.+c1abcem2*dt/dxa))
     
           un(i1-2*is1,i2-2*is2,i3-2*is3,ex)=4.*un(i1-is1,i2-is2,i3-is3,ex)-6.*un(i1,i2,i3,ex)+4.*un(i1+is1,i2+is2,i3+is3,ex)-un(i1+2*is1,i2+2*is2,i3+2*is3,ex)

        end if
        end do
        end do
        end do
      

     else if( axis.eq.1 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
       un(i1-is1,i2-is2,i3-is3,ex)=((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dya*dt)*(c1abcem2*uyy22r(i1,i2,i3,ex)+c2abcem2*uxx22r(i1,i2,i3,ex)+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1,i2+is2,i3,ex))/dya**2+c2abcem2*(un(i1-1,i2,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1+1,i2,i3,ex))/dx(0)**2+2.*forcey(ex)))/(1.+c1abcem2*dt/dya))

         un(i1-2*is1,i2-2*is2,i3-2*is3,ex)=4.*un(i1-is1,i2-is2,i3-is3,ex)-6.*un(i1,i2,i3,ex)+4.*un(i1+is1,i2+is2,i3+is3,ex)-un(i1+2*is1,i2+2*is2,i3+2*is3,ex)

      end if
      end do
      end do
      end do
     else
      stop 9469
     end if

    else ! ***** 3D *****
     if( axis.eq.0 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
       un(i1-is1,i2-is2,i3-is3,ex)=((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dxa*dt)*(c1abcem2*uxx23r(i1,i2,i3,ex)+c2abcem2*(uyy23r(i1,i2,i3,ex)+uzz23r(i1,i2,i3,ex))+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1+is1,i2,i3,ex))/dxa**2+c2abcem2*(un(i1,i2-1,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2+1,i3,ex))/dx(1)**2+c2abcem2*(un(i1,i2,i3-1,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2,i3+1,ex))/dx(2)**2+2.*forcex(ex)))/(1.+c1abcem2*dt/dxa))

         un(i1-2*is1,i2-2*is2,i3-2*is3,ex)=4.*un(i1-is1,i2-is2,i3-is3,ex)-6.*un(i1,i2,i3,ex)+4.*un(i1+is1,i2+is2,i3+is3,ex)-un(i1+2*is1,i2+2*is2,i3+2*is3,ex)

       ! write(1,'("t=",e10.2,", i1,i2,i3=",3i3," div43d=",e10.2)') t,i1,i2,i3,unx43r(i1,i2,i3,ex)+uny43r(i1,i2,i3,ey)+unz43r(i1,i2,i3,ez)

      end if
      end do
      end do
      end do
     else if( axis.eq.1 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
       un(i1-is1,i2-is2,i3-is3,ex)=((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dya*dt)*(c1abcem2*uyy23r(i1,i2,i3,ex)+c2abcem2*(uxx23r(i1,i2,i3,ex)+uzz23r(i1,i2,i3,ex))+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1,i2+is2,i3,ex))/dya**2+c2abcem2*(un(i1-1,i2,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1+1,i2,i3,ex))/dx(0)**2+c2abcem2*(un(i1,i2,i3-1,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2,i3+1,ex))/dx(2)**2+2.*forcey(ex)))/(1.+c1abcem2*dt/dya))
 
      end if
      end do
      end do
      end do
     else if( axis.eq.2 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
       un(i1-is1,i2-is2,i3-is3,ex)=((un(i1+is1,i2+is2,i3+is3,ex)-(u(i1+is1,i2+is2,i3+is3,ex)-u(i1-is1,i2-is2,i3-is3,ex))-(dza*dt)*(c1abcem2*uzz23r(i1,i2,i3,ex)+c2abcem2*(uxx23r(i1,i2,i3,ex)+uyy23r(i1,i2,i3,ex))+c1abcem2*(-2.*un(i1,i2,i3,ex)+un(i1,i2,i3+is3,ex))/dza**2+c2abcem2*(un(i1-1,i2,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1+1,i2,i3,ex))/dx(0)**2+c2abcem2*(un(i1,i2-1,i3,ex)-2.*un(i1,i2,i3,ex)+un(i1,i2+1,i3,ex))/dx(1)**2+2.*forcez(ex)))/(1.+c1abcem2*dt/dza))

         un(i1-2*is1,i2-2*is2,i3-2*is3,ex)=4.*un(i1-is1,i2-is2,i3-is3,ex)-6.*un(i1,i2,i3,ex)+4.*un(i1+is1,i2+is2,i3+is3,ex)-un(i1+2*is1,i2+2*is2,i3+2*is3,ex)

      ! write(1,'("t=",e10.2,", i1,i2,i3=",3i3," div43d=",e10.2)') t,i1,i2,i3,unx43r(i1,i2,i3,ex)+uny43r(i1,i2,i3,ey)+unz43r(i1,i2,i3,ez)

      end if
      end do
      end do
      end do
     else
      stop 9470
     end if

    end if

    ! added: July 3, 2023 
    if( numGhost>2 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
        do ghost=3,numGhost
          j1= i1-ghost*is1; j2=i2-ghost*is2; j3=i3-ghost*is3
          un(j1,j2,j3,ex)=(5.*un(j1+is1,j2+is2,j3+is3,ex)-10.*un(j1+is1+is1,j2+is2+is2,j3+is3+is3,ex)+10.*un(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,ex)-5.*un(j1+is1+3*is1,j2+is2+3*is2,j3+is3+3*is3,ex)+un(j1+is1+4*is1,j2+is2+4*is2,j3+is3+4*is3,ex))
      end do
      end if
      end do
      end do
      end do
    end if

 
   else if( gridType.eq.curvilinear .and. orderOfAccuracy.eq.2 )then
    ! ***********************************************
    ! ************curvilinear grid*******************
    ! ***********************************************

    ! --- this really assumes that the boundary is not curved : do this for now ----

    ! write(*,'(" Apply abcEM2 CURV ORDER 2: grid,side,axis=",3i3," dt,c=",2e12.3)') grid,side,axis,dt,c

    if( nd.eq.2 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then

         is =1-2*side
         dr0=dr(axis)
         rx0 = rsxy(i1,i2,i3,axis,0)
         ry0 = rsxy(i1,i2,i3,axis,1)
         rxNormSq = rx0**2 + ry0**2 
         rxNorm = max( epsX, sqrt(rxNormSq) )
         rxx0 = rsxyx22(i1,i2,i3,axis,0)
         ryy0 = rsxyy22(i1,i2,i3,axis,1)
         ! cm1 : coeff of u(i1-is1,i2-is2,i3-is3,cc) in g (given below): 
         cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0)/(2.*dr0) ) )
         ! u: derivatives at time t: 
         ! un: derivatives at time t+dt : evaluate using the incorrect ghost values 
         if( axis.eq.0 )then
           ur0   =   ur2(i1,i2,i3,ex)
           urr0  =  urr2(i1,i2,i3,ex)
           unr0  =  unr2(i1,i2,i3,ex)
           unrr0 = unrr2(i1,i2,i3,ex)
         else
           ur0   =   us2(i1,i2,i3,ex)
           urr0  =  uss2(i1,i2,i3,ex)
           unr0  =  uns2(i1,i2,i3,ex)
           unrr0 = unss2(i1,i2,i3,ex)
         end if
         uxx0 = uxx22(i1,i2,i3,ex)
         uyy0 = uyy22(i1,i2,i3,ex)
         unxx0 = unxx22(i1,i2,i3,ex)
         unyy0 = unyy22(i1,i2,i3,ex)
         ! first evaluate the BC using the incorrect ghost values 
         Dn2 = rxNormSq*(unrr0+urr0) 
         Lu  = unxx0+unyy0 + uxx0+uyy0 
         g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
         ! note: this assumes an orthogonal grid -- we should make sure that the 
         !       ghost values have an initial guess in them (extrapolate ?)
         un(i1-is1,i2-is2,i3-is3,ex) = -(g - cm1*un(i1-is1,i2-is2,i3-is3,ex) )/cm1 

      end if
      end do
      end do
      end do

    else ! ***** 3D *****

      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
         is =1-2*side
         dr0=dr(axis)
         rx0 = rsxy(i1,i2,i3,axis,0)
         ry0 = rsxy(i1,i2,i3,axis,1)
         rz0 = rsxy(i1,i2,i3,axis,2)
         rxNormSq = rx0**2 + ry0**2 + rz0**2
         rxNorm = max( epsX, sqrt(rxNormSq) )
         rxx0 = rsxyx23(i1,i2,i3,axis,0)
         ryy0 = rsxyy23(i1,i2,i3,axis,1)
         rzz0 = rsxyz23(i1,i2,i3,axis,2)
         ! cm1 : coeff of u(i1-is1,i2-is2,i3-is3,cc) in g (given below): 
         cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0+rzz0)/(2.*dr0) ) )
         ! derivatives at time t: 
         ! derivatives at time t+dt : evaluate using the incorrect ghost values 
         if( axis.eq.0 )then
           ur0   =   ur2(i1,i2,i3,ex)
           urr0  =  urr2(i1,i2,i3,ex)
           unr0  =  unr2(i1,i2,i3,ex)
           unrr0 = unrr2(i1,i2,i3,ex)
         else if( axis.eq.1 )then
           ur0   =   us2(i1,i2,i3,ex)
           urr0  =  uss2(i1,i2,i3,ex)
           unr0  =  uns2(i1,i2,i3,ex)
           unrr0 = unss2(i1,i2,i3,ex)
         else
           ur0   =   ut2(i1,i2,i3,ex)
           urr0  =  utt2(i1,i2,i3,ex)
           unr0  =  unt2(i1,i2,i3,ex)
           unrr0 = untt2(i1,i2,i3,ex)
         end if
         uxx0 = uxx23(i1,i2,i3,ex)
         uyy0 = uyy23(i1,i2,i3,ex)
         uzz0 = uzz23(i1,i2,i3,ex)
         unxx0 = unxx23(i1,i2,i3,ex)
         unyy0 = unyy23(i1,i2,i3,ex)
         unzz0 = unzz23(i1,i2,i3,ex)
         ! first evaluate the BC using the incorrect ghost values 
         Dn2 = rxNormSq*(unrr0+urr0) 
         Lu  = unxx0 + unyy0 + unzz0 + uxx0 + uyy0 + uzz0 
         g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
         ! note: this assumes an orthogonal grid -- we should make sure that the 
         !       ghost values have an initial guess in them (extrapolate ?)
         un(i1-is1,i2-is2,i3-is3,ex) = -(g - cm1*un(i1-is1,i2-is2,i3-is3,ex) )/cm1 

      end if
      end do
      end do
      end do

    end if

    ! added: July 3, 2023 
    if( numGhost>1 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
        do ghost=2,numGhost
          j1= i1-ghost*is1; j2=i2-ghost*is2; j3=i3-ghost*is3
          un(j1,j2,j3,ex)=(3.*un(j1+is1,j2+is2,j3+is3,ex)-3.*un(j1+is1+is1,j2+is2+is2,j3+is3+is3,ex)+un(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,ex))
        end do
      end if
      end do
      end do
      end do
    end if

 

   else if( gridType.eq.curvilinear .and. orderOfAccuracy.eq.4 )then

    ! ***********************************************
    ! ************curvilinear grid*******************
    ! ***********************************************

   
     ! >>>>>>>>>>>>>>>> this is only second-order ---- fix this <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       

    ! write(*,'(" Apply abcEM2: grid,side,axis=",3i3," dt,c=",2e12.3)') grid,side,axis,dt,c
    if( nd.eq.2 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
   
         is =1-2*side
         dr0=dr(axis)
         rx0 = rsxy(i1,i2,i3,axis,0)
         ry0 = rsxy(i1,i2,i3,axis,1)
         rxNormSq = rx0**2 + ry0**2 
         rxNorm = max( epsX, sqrt(rxNormSq) )
         rxx0 = rsxyx22(i1,i2,i3,axis,0)
         ryy0 = rsxyy22(i1,i2,i3,axis,1)
         ! cm1 : coeff of u(i1-is1,i2-is2,i3-is3,cc) in g (given below): 
         cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0)/(2.*dr0) ) )
         ! u: derivatives at time t: 
         ! un: derivatives at time t+dt : evaluate using the incorrect ghost values 
         if( axis.eq.0 )then
           ur0   =   ur2(i1,i2,i3,ex)
           urr0  =  urr2(i1,i2,i3,ex)
           unr0  =  unr2(i1,i2,i3,ex)
           unrr0 = unrr2(i1,i2,i3,ex)
         else
           ur0   =   us2(i1,i2,i3,ex)
           urr0  =  uss2(i1,i2,i3,ex)
           unr0  =  uns2(i1,i2,i3,ex)
           unrr0 = unss2(i1,i2,i3,ex)
         end if
         uxx0 = uxx22(i1,i2,i3,ex)
         uyy0 = uyy22(i1,i2,i3,ex)
         unxx0 = unxx22(i1,i2,i3,ex)
         unyy0 = unyy22(i1,i2,i3,ex)
         ! first evaluate the BC using the incorrect ghost values 
         Dn2 = rxNormSq*(unrr0+urr0) 
         Lu  = unxx0+unyy0 + uxx0+uyy0 
         g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
         ! note: this assumes an orthogonal grid -- we should make sure that the 
         !       ghost values have an initial guess in them (extrapolate ?)
         un(i1-is1,i2-is2,i3-is3,ex) = -(g - cm1*un(i1-is1,i2-is2,i3-is3,ex) )/cm1 
   
         un(i1-2*is1,i2-2*is2,i3-2*is3,ex)=4.*un(i1-is1,i2-is2,i3-is3,ex)-6.*un(i1,i2,i3,ex)+4.*un(i1+is1,i2+is2,i3+is3,ex)-un(i1+2*is1,i2+2*is2,i3+2*is3,ex)

      end if
      end do
      end do
      end do

    else ! ***** 3D *****

      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
         is =1-2*side
         dr0=dr(axis)
         rx0 = rsxy(i1,i2,i3,axis,0)
         ry0 = rsxy(i1,i2,i3,axis,1)
         rz0 = rsxy(i1,i2,i3,axis,2)
         rxNormSq = rx0**2 + ry0**2 + rz0**2
         rxNorm = max( epsX, sqrt(rxNormSq) )
         rxx0 = rsxyx23(i1,i2,i3,axis,0)
         ryy0 = rsxyy23(i1,i2,i3,axis,1)
         rzz0 = rsxyz23(i1,i2,i3,axis,2)
         ! cm1 : coeff of u(i1-is1,i2-is2,i3-is3,cc) in g (given below): 
         cm1 = -rxNorm/(2.*dr0*dt) -.5*( c1abcem2*( rxNormSq/dr0**2 ) + c2abcem2*( -is*(rxx0+ryy0+rzz0)/(2.*dr0) ) )
         ! derivatives at time t: 
         ! derivatives at time t+dt : evaluate using the incorrect ghost values 
         if( axis.eq.0 )then
           ur0   =   ur2(i1,i2,i3,ex)
           urr0  =  urr2(i1,i2,i3,ex)
           unr0  =  unr2(i1,i2,i3,ex)
           unrr0 = unrr2(i1,i2,i3,ex)
         else if( axis.eq.1 )then
           ur0   =   us2(i1,i2,i3,ex)
           urr0  =  uss2(i1,i2,i3,ex)
           unr0  =  uns2(i1,i2,i3,ex)
           unrr0 = unss2(i1,i2,i3,ex)
         else
           ur0   =   ut2(i1,i2,i3,ex)
           urr0  =  utt2(i1,i2,i3,ex)
           unr0  =  unt2(i1,i2,i3,ex)
           unrr0 = untt2(i1,i2,i3,ex)
         end if
         uxx0 = uxx23(i1,i2,i3,ex)
         uyy0 = uyy23(i1,i2,i3,ex)
         uzz0 = uzz23(i1,i2,i3,ex)
         unxx0 = unxx23(i1,i2,i3,ex)
         unyy0 = unyy23(i1,i2,i3,ex)
         unzz0 = unzz23(i1,i2,i3,ex)
         ! first evaluate the BC using the incorrect ghost values 
         Dn2 = rxNormSq*(unrr0+urr0) 
         Lu  = unxx0 + unyy0 + unzz0 + uxx0 + uyy0 + uzz0 
         g = is*rxNorm*(unr0-ur0)/dt -.5*( c1abcem2*( Dn2 ) + c2abcem2*( Lu - Dn2 ) )
         ! note: this assumes an orthogonal grid -- we should make sure that the 
         !       ghost values have an initial guess in them (extrapolate ?)
         un(i1-is1,i2-is2,i3-is3,ex) = -(g - cm1*un(i1-is1,i2-is2,i3-is3,ex) )/cm1 

         un(i1-2*is1,i2-2*is2,i3-2*is3,ex)=4.*un(i1-is1,i2-is2,i3-is3,ex)-6.*un(i1,i2,i3,ex)+4.*un(i1+is1,i2+is2,i3+is3,ex)-un(i1+2*is1,i2+2*is2,i3+2*is3,ex)

       ! write(1,'("t=",e10.2,", i1,i2,i3=",3i3," div43dc=",e10.2)') t,i1,i2,i3,unx43(i1,i2,i3,ex)+uny43(i1,i2,i3,ey)+unz43(i1,i2,i3,ez)

      end if
      end do
      end do
      end do

    end if
 

    ! added: July 3, 2023 
    if( numGhost>2 )then
      do i3=n3a,n3b
      do i2=n2a,n2b
      do i1=n1a,n1b
      if( mask(i1,i2,i3).gt.0 )then
        do ghost=3,numGhost
          j1= i1-ghost*is1; j2=i2-ghost*is2; j3=i3-ghost*is3
          un(j1,j2,j3,ex)=(5.*un(j1+is1,j2+is2,j3+is3,ex)-10.*un(j1+is1+is1,j2+is2+is2,j3+is3+is3,ex)+10.*un(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,ex)-5.*un(j1+is1+3*is1,j2+is2+3*is2,j3+is3+3*is3,ex)+un(j1+is1+4*is1,j2+is2+4*is2,j3+is3+4*is3,ex))
        end do
      end if
      end do
      end do
      end do
    end if


   else
     write(*,'(" ABC: ERROR gridType=",i2," orderOfAccuracy=",i4," unexpected")') gridType,orderOfAccuracy

     stop 2255
   end if

     end if
   end do
   end do

  return
  end
