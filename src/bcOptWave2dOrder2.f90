! This file automatically generated from bcOptWave.bf90 with bpp.
 subroutine bcOptWave2dOrder2( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
 ! ===================================================================================
 !  Boundary conditions for CgWave
 !
 !  gridType : 0=rectangular, 1=curvilinear
 !
 ! The forcing for the boundary conditions can be accessed using the statement function:
 !         bcf(side,axis,i1,i2,i3,m)
 ! which is defined below. 
 ! ===================================================================================
   implicit none
   integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, ndb, ierr
   real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
   integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
   real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
   real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
   integer gridIndexRange(0:1,0:2),boundaryCondition(0:1,0:2), dimRange(0:1,0:2), isPeriodic(0:*)
   real frequencyArray(0:*)
   ! temp space for CBC order 4 -- fix me : just make a stencil
   real uTemp(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
   real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
   ! ! *** TEMP ARRAYS FOR WORK SPACE --> THIS IS SLOW!!!
   ! #If "2" eq "4"
   !   real uTemp(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
   !   real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
   ! #Else
   ! real v(2,2,2,0:0)
   ! #End
   double precision pdb  ! pointer to data base
   ! integer addBoundaryForcing(0:1,0:2)
   ! integer interfaceType(0:1,0:2,0:*)
   ! integer dim(0:1,0:2,0:1,0:2)
   ! real bcf0(0:*)
   ! integer*8 bcOffset(0:1,0:2)
   ! real bcData(0:ndb-1,0:1,0:nd-1,0:*)
   integer ipar(0:*)
   real rpar(0:*)
   !     --- local variables ----
   integer bc(0:1,0:2) ! local version, normally equal to boundaryCondition
   integer uc,numberOfComponents,assignTwilightZone,assignKnownSolutionAtBoundaries,freq
   integer grid,gridType,orderOfAccuracy,useWhereMask,gridIsImplicit,useUpwindDissipation
   integer twilightZone,numberOfProcessors,addForcingBC,assignBCForImplicit
   integer debug,myid,ghost
   integer checkCoeff
   real maxDiff
   integer ok,getInt,getReal
   real omega,cfl,c
   real kx,ky,kz,twoPi
   real t,dt,epsx,REAL_MIN 
   real ep
   real a0,a1,an1,an2,an3,aNormi, t1,t2,t3
   real dx(0:2),dr(0:2),gravity(0:2)
   real r1,r2,r3
   real dxn,b0,b1,ue,uex,uey,uez,ff,urv(0:2),ur0,cosPW
   real c2,c4,c6,c8
   real gtt,rFactor,uLap,vLap
   real a11,a12,a21,a22
   real r1a,r2a, r1b,r2b, a11c,a12c,a21c,a22c  
   real u1Save,u2Save  
   real uett,uexx,ueyy,uezz,ueLap
   real ue1,ue2,ue3,ue4,f1,f2,f3,det
   real uettxx,uettyy,uettzz, uexxxx,ueyyyy,uezzzz,uexxyy,uexxzz,ueyyzz
   real uetttt,uettLap,ueLap2,lap3d2Pow2
   real uettx,uetty,uettz
   real uexxx,uexxy,uexxz,uexyy,uexzz,ueyyy,ueyyz,ueyzz,uezzz
   real uexxxxxx,uexxxxyy,uexxxxzz,ueyyyyyy,uexxyyyy,uexxzzzz,ueyyyyzz,ueyyzzzz,uezzzzzz,uexxyyzz
   real fLap,ftt,gtttt
   real gg,nDotGradF,crv(0:3)
   integer side,axis,axisp1,axisp2,i1,i2,i3,is1,is2,is3,j1,j2,j3,js1,js2,js3,k1,k2,k3,ks1,ks2,ks3,is,js
   integer i1p,i2p,i3p
   integer l1,l2,l3
   integer numGhost,numGhost3,numberOfGhostPoints,extraForNeumann,extraForDirichlet,numberOfFrequencies
   integer side1,side2,side3
   integer n1a,n1b,n2a,n2b,n3a,n3b
   integer m1a,m1b,m2a,m2b,m3a,m3b
   integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
   integer j1a,j1b,j2a,j2b,ja,j3b
   integer extra1a,extra1b,extra2a,extra2b,extra3a,extra3b,extram
   integer maxExtrapWidth,extrapWidth
   integer cornerBC(0:2,0:2,0:2), iparc(0:10), orderOfExtrapolationForCorners, orderOfExtrapolation
   real rparc(0:10)
   real ca,cEM2,rhs
   ! boundary conditions parameters and interfaceType values
   ! #Include "bcDefineFortran.h"
   ! These should mauch the values in Parameters.h
   integer dirichletBoundaryCondition,neumannBoundaryCondition,dirichletInterface,neumannInterface,mixedBoundaryCondition
   parameter( dirichletBoundaryCondition=12, neumannBoundaryCondition=18, dirichletInterface=21, neumannInterface=22, mixedBoundaryCondition=30 )
   integer rectangular,curvilinear
   parameter(rectangular=0,curvilinear=1)
   ! Boundary conditions: These must mauch the values in CgWave.h
   ! periodic       =-1,
   ! interpolation  = 0,
   ! dirichlet      = 1,
   ! neumann        = 2,
   ! evenSymmetry   = 3,
   ! radiation      = 4   
   ! exactBC        = 5 
   ! abcEM2         = 6,  // absorbing BC, Engquist-Majda order 2  
   ! characteristic = 7,  // characteristic BC
   ! absorbing      = 8,   // for SuperGrid  
   integer dirichlet,neumann,evenSymmetry,radiation,exactBC,abcEM2,characteristic,absorbing
   parameter( dirichlet=1, neumann=2, evenSymmetry=3, radiation=4, exactBC=5, abcEM2=6, characteristic=7, absorbing=8  )
   ! Corner conditions (from op/fortranDeriv/assignCornersOpt.bf)
   integer doNothingCorner,extrapolateCorner,symmetryCorner,taylor2ndOrder
   integer evenSymmetryCorner,oddSymmetryCorner,taylor2ndOrderEvenCorner,taylor4thOrderEvenCorner,vectorSymmetryAxis1Corner,vectorSymmetryAxis2Corner,vectorSymmetryAxis3Corner
   parameter(doNothingCorner=-1,extrapolateCorner=0,symmetryCorner=1,taylor2ndOrder=2, evenSymmetryCorner=3,oddSymmetryCorner=4,taylor2ndOrderEvenCorner=5,taylor4thOrderEvenCorner=6, vectorSymmetryAxis1Corner=7,vectorSymmetryAxis2Corner=8,vectorSymmetryAxis3Corner=9 )      
   ! known solutions
   integer knownSolutionOption
   integer planeWave, gaussianPlaneWave, boxHelmHoltz, polyPeriodic, otherKnownSolution
   parameter( planeWave=1, gaussianPlaneWave=2, boxHelmHoltz=3, polyPeriodic=4, otherKnownSolution=1000 )
   ! parameters for plane wave known solution
   real ampPlaneWave, kxPlaneWave,kyPlaneWave,kzPlaneWave, omegaPlaneWave, omegaTol
   integer solveHelmholtz, solveForScatteredField
   ! parameters for Gaussian plane wave
   real kxGPW,kyGPW,kzGPW, x0GPW,y0GPW,z0GPW, k0GPW, betaGPW
   real xi
   ! parameters for boxHelmholtz known solution
   real kxBoxHelmholtz,kyBoxHelmholtz,kzBoxHelmholtz,omegaBoxHelmholtz,coswt
   ! parameters for polyPeriodic known solution
   real omegaPolyPeriodic,a0PolyPeriodic, a1PolyPeriodic, b1PolyPeriodic, c1PolyPeriodic
   ! --- forcing options ----
   ! These must match the values in CgWave.h: 
   ! enum ForcingOptionEnum
   ! {
   !   noForcing=0,
   !   twilightZoneForcing,
   !   userForcing,
   !   helmholtzForcing
   ! };  
   integer forcingOption
   integer noForcing,twilightZoneForcing,userForcing,helmholtzForcing
   parameter( noForcing=0, twilightZoneForcing=1, userForcing=2, helmholtzForcing=2 )
   ! BC APPROACH -- these must match the values in CgWave.h 
   ! enum BoundaryConditionApproachEnum
   ! {
   !   defaultBoundaryConditionApproach,
   !   useOneSidedBoundaryConditions,
   !   useCompatibilityBoundaryConditions,
   !   useLocalCompatibilityBoundaryConditions
   ! };  
   integer bcApproach
   integer defaultBoundaryConditionApproach
   integer useOneSidedBoundaryConditions
   integer useCompatibilityBoundaryConditions
   integer useLocalCompatibilityBoundaryConditions  
   parameter( defaultBoundaryConditionApproach       =0, useOneSidedBoundaryConditions          =1, useCompatibilityBoundaryConditions     =2, useLocalCompatibilityBoundaryConditions=3 )
   real r3v(0:2),a3(0:2,0:2),a3i(0:2,0:2),f3v(0:2)
   real scale1,scale2,scale3
   integer m1,m2,m3
   integer assignCornerGhostPoints
   logical firstTimeForCBC6
   real symSign
   integer sidea, iab(0:1,0:2)
   !     --- start statement function ----
   real bcf,mixedRHS,mixedCoeff,mixedNormalCoeff
   integer kd,m,n,component
   real uxOneSided,lap2d2Pow2
   ! real uxxx,uxxy,uxxz,uxyy,uxzz, uyyy,uyyz,uyzz, uzzz 
   real rx,ry,rz,sx,sy,sz,tx,ty,tz
   ! define variables for getDerivatives macros
   ! #Include "../maple/declareGetDerivativesMacrosVariables.h"
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
    real vr2
    real vs2
    real vt2
    real vrr2
    real vss2
    real vrs2
    real vtt2
    real vrt2
    real vst2
    real vrrr2
    real vsss2
    real vttt2
    real vx21
    real vy21
    real vz21
    real vx22
    real vy22
    real vz22
    real vx23
    real vy23
    real vz23
    real vxx21
    real vyy21
    real vxy21
    real vxz21
    real vyz21
    real vzz21
    real vlaplacian21
    real vxx22
    real vyy22
    real vxy22
    real vxz22
    real vyz22
    real vzz22
    real vlaplacian22
    real vxx23
    real vyy23
    real vzz23
    real vxy23
    real vxz23
    real vyz23
    real vlaplacian23
    real vx23r
    real vy23r
    real vz23r
    real vxx23r
    real vyy23r
    real vxy23r
    real vzz23r
    real vxz23r
    real vyz23r
    real vx21r
    real vy21r
    real vz21r
    real vxx21r
    real vyy21r
    real vzz21r
    real vxy21r
    real vxz21r
    real vyz21r
    real vlaplacian21r
    real vx22r
    real vy22r
    real vz22r
    real vxx22r
    real vyy22r
    real vzz22r
    real vxy22r
    real vxz22r
    real vyz22r
    real vlaplacian22r
    real vlaplacian23r
    real vxxx22r
    real vyyy22r
    real vxxy22r
    real vxyy22r
    real vxxxx22r
    real vyyyy22r
    real vxxyy22r
    real vxxx23r
    real vyyy23r
    real vzzz23r
    real vxxy23r
    real vxxz23r
    real vxyy23r
    real vyyz23r
    real vxzz23r
    real vyzz23r
    real vxxxx23r
    real vyyyy23r
    real vzzzz23r
    real vxxyy23r
    real vxxzz23r
    real vyyzz23r
    real vLapSq22r
    real vLapSq23r
   !  The next macro call will define the difference approximation statement functions
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
   vr2(i1,i2,i3,kd)=(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))*d12(0)
   vs2(i1,i2,i3,kd)=(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))*d12(1)
   vt2(i1,i2,i3,kd)=(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))*d12(2)
   vrr2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)) )*d22(0)
   vss2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd)) )*d22(1)
   vrs2(i1,i2,i3,kd)=(vr2(i1,i2+1,i3,kd)-vr2(i1,i2-1,i3,kd))*d12(1)
   vtt2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd)) )*d22(2)
   vrt2(i1,i2,i3,kd)=(vr2(i1,i2,i3+1,kd)-vr2(i1,i2,i3-1,kd))*d12(2)
   vst2(i1,i2,i3,kd)=(vs2(i1,i2,i3+1,kd)-vs2(i1,i2,i3-1,kd))*d12(2)
   vrrr2(i1,i2,i3,kd)=(-2.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
   vsss2(i1,i2,i3,kd)=(-2.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
   vttt2(i1,i2,i3,kd)=(-2.*(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))+(v(i1,i2,i3+2,kd)-v(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
   vx21(i1,i2,i3,kd)= rx(i1,i2,i3)*vr2(i1,i2,i3,kd)
   vy21(i1,i2,i3,kd)=0
   vz21(i1,i2,i3,kd)=0
   vx22(i1,i2,i3,kd)= rx(i1,i2,i3)*vr2(i1,i2,i3,kd)+sx(i1,i2,i3)*vs2(i1,i2,i3,kd)
   vy22(i1,i2,i3,kd)= ry(i1,i2,i3)*vr2(i1,i2,i3,kd)+sy(i1,i2,i3)*vs2(i1,i2,i3,kd)
   vz22(i1,i2,i3,kd)=0
   vx23(i1,i2,i3,kd)=rx(i1,i2,i3)*vr2(i1,i2,i3,kd)+sx(i1,i2,i3)*vs2(i1,i2,i3,kd)+tx(i1,i2,i3)*vt2(i1,i2,i3,kd)
   vy23(i1,i2,i3,kd)=ry(i1,i2,i3)*vr2(i1,i2,i3,kd)+sy(i1,i2,i3)*vs2(i1,i2,i3,kd)+ty(i1,i2,i3)*vt2(i1,i2,i3,kd)
   vz23(i1,i2,i3,kd)=rz(i1,i2,i3)*vr2(i1,i2,i3,kd)+sz(i1,i2,i3)*vs2(i1,i2,i3,kd)+tz(i1,i2,i3)*vt2(i1,i2,i3,kd)
   vxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*vr2(i1,i2,i3,kd)
   vyy21(i1,i2,i3,kd)=0
   vxy21(i1,i2,i3,kd)=0
   vxz21(i1,i2,i3,kd)=0
   vyz21(i1,i2,i3,kd)=0
   vzz21(i1,i2,i3,kd)=0
   vlaplacian21(i1,i2,i3,kd)=vxx21(i1,i2,i3,kd)
   vxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*vss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*vr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3))*vs2(i1,i2,i3,kd)
   vyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*vss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*vr2(i1,i2,i3,kd)+(syy22(i1,i2,i3))*vs2(i1,i2,i3,kd)
   vxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*vrr2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*vss2(i1,i2,i3,kd)+rxy22(i1,i2,i3)*vr2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*vs2(i1,i2,i3,kd)
   vxz22(i1,i2,i3,kd)=0
   vyz22(i1,i2,i3,kd)=0
   vzz22(i1,i2,i3,kd)=0
   vlaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*vss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*vr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*vs2(i1,i2,i3,kd)
   vxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*vrr2(i1,i2,i3,kd)+sx(i1,i2,i3)**2*vss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*vtt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*vrs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*vrt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*vst2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*vr2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*vs2(i1,i2,i3,kd)+txx23(i1,i2,i3)*vt2(i1,i2,i3,kd)
   vyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*vrr2(i1,i2,i3,kd)+sy(i1,i2,i3)**2*vss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*vtt2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*vrs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*vrt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*vst2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*vr2(i1,i2,i3,kd)+syy23(i1,i2,i3)*vs2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*vt2(i1,i2,i3,kd)
   vzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*vrr2(i1,i2,i3,kd)+sz(i1,i2,i3)**2*vss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*vtt2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*vrs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*vrt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*vst2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*vr2(i1,i2,i3,kd)+szz23(i1,i2,i3)*vs2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*vt2(i1,i2,i3,kd)
   vxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*vrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*vss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*vtt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*vrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*vst2(i1,i2,i3,kd)+rxy23(i1,i2,i3)*vr2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*vs2(i1,i2,i3,kd)+txy23(i1,i2,i3)*vt2(i1,i2,i3,kd)
   vxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*vrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*vss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*vtt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*vrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*vst2(i1,i2,i3,kd)+rxz23(i1,i2,i3)*vr2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*vs2(i1,i2,i3,kd)+txz23(i1,i2,i3)*vt2(i1,i2,i3,kd)
   vyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*vrr2(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*vss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*vtt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*vrt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*vst2(i1,i2,i3,kd)+ryz23(i1,i2,i3)*vr2(i1,i2,i3,kd)+syz23(i1,i2,i3)*vs2(i1,i2,i3,kd)+tyz23(i1,i2,i3)*vt2(i1,i2,i3,kd)
   vlaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*vss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*vtt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*vrs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*vrt2(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*vst2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*vr2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*vs2(i1,i2,i3,kd)+(txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*vt2(i1,i2,i3,kd)
   !============================================================================================
   ! Define derivatives for a rectangular grid
   !
   !============================================================================================
   vx23r(i1,i2,i3,kd)=(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))*h12(0)
   vy23r(i1,i2,i3,kd)=(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))*h12(1)
   vz23r(i1,i2,i3,kd)=(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))*h12(2)
   vxx23r(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)) )*h22(0)
   vyy23r(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd)) )*h22(1)
   vxy23r(i1,i2,i3,kd)=(vx23r(i1,i2+1,i3,kd)-vx23r(i1,i2-1,i3,kd))*h12(1)
   vzz23r(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd)) )*h22(2)
   vxz23r(i1,i2,i3,kd)=(vx23r(i1,i2,i3+1,kd)-vx23r(i1,i2,i3-1,kd))*h12(2)
   vyz23r(i1,i2,i3,kd)=(vy23r(i1,i2,i3+1,kd)-vy23r(i1,i2,i3-1,kd))*h12(2)
   vx21r(i1,i2,i3,kd)= vx23r(i1,i2,i3,kd)
   vy21r(i1,i2,i3,kd)= vy23r(i1,i2,i3,kd)
   vz21r(i1,i2,i3,kd)= vz23r(i1,i2,i3,kd)
   vxx21r(i1,i2,i3,kd)= vxx23r(i1,i2,i3,kd)
   vyy21r(i1,i2,i3,kd)= vyy23r(i1,i2,i3,kd)
   vzz21r(i1,i2,i3,kd)= vzz23r(i1,i2,i3,kd)
   vxy21r(i1,i2,i3,kd)= vxy23r(i1,i2,i3,kd)
   vxz21r(i1,i2,i3,kd)= vxz23r(i1,i2,i3,kd)
   vyz21r(i1,i2,i3,kd)= vyz23r(i1,i2,i3,kd)
   vlaplacian21r(i1,i2,i3,kd)=vxx23r(i1,i2,i3,kd)
   vx22r(i1,i2,i3,kd)= vx23r(i1,i2,i3,kd)
   vy22r(i1,i2,i3,kd)= vy23r(i1,i2,i3,kd)
   vz22r(i1,i2,i3,kd)= vz23r(i1,i2,i3,kd)
   vxx22r(i1,i2,i3,kd)= vxx23r(i1,i2,i3,kd)
   vyy22r(i1,i2,i3,kd)= vyy23r(i1,i2,i3,kd)
   vzz22r(i1,i2,i3,kd)= vzz23r(i1,i2,i3,kd)
   vxy22r(i1,i2,i3,kd)= vxy23r(i1,i2,i3,kd)
   vxz22r(i1,i2,i3,kd)= vxz23r(i1,i2,i3,kd)
   vyz22r(i1,i2,i3,kd)= vyz23r(i1,i2,i3,kd)
   vlaplacian22r(i1,i2,i3,kd)=vxx23r(i1,i2,i3,kd)+vyy23r(i1,i2,i3,kd)
   vlaplacian23r(i1,i2,i3,kd)=vxx23r(i1,i2,i3,kd)+vyy23r(i1,i2,i3,kd)+vzz23r(i1,i2,i3,kd)
   vxxx22r(i1,i2,i3,kd)=(-2.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
   vyyy22r(i1,i2,i3,kd)=(-2.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
   vxxy22r(i1,i2,i3,kd)=( vxx22r(i1,i2+1,i3,kd)-vxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
   vxyy22r(i1,i2,i3,kd)=( vyy22r(i1+1,i2,i3,kd)-vyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
   vxxxx22r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )/(dx(0)**4)
   vyyyy22r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(1)**4)
   vxxyy22r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)+v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))   +   (v(i1+1,i2+1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
   ! 2D laplacian squared = v.xxxx + 2 v.xxyy + v.yyyy
   vLapSq22r(i1,i2,i3,kd)= ( 6.*v(i1,i2,i3,kd)   - 4.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd))    +(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*v(i1,i2,i3,kd)    -4.*(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))    +(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 8.*v(i1,i2,i3,kd)     -4.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)+v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))   +2.*(v(i1+1,i2+1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
   vxxx23r(i1,i2,i3,kd)=(-2.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
   vyyy23r(i1,i2,i3,kd)=(-2.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
   vzzz23r(i1,i2,i3,kd)=(-2.*(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))+(v(i1,i2,i3+2,kd)-v(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
   vxxy23r(i1,i2,i3,kd)=( vxx22r(i1,i2+1,i3,kd)-vxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
   vxyy23r(i1,i2,i3,kd)=( vyy22r(i1+1,i2,i3,kd)-vyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
   vxxz23r(i1,i2,i3,kd)=( vxx22r(i1,i2,i3+1,kd)-vxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
   vyyz23r(i1,i2,i3,kd)=( vyy22r(i1,i2,i3+1,kd)-vyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
   vxzz23r(i1,i2,i3,kd)=( vzz22r(i1+1,i2,i3,kd)-vzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
   vyzz23r(i1,i2,i3,kd)=( vzz22r(i1,i2+1,i3,kd)-vzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
   vxxxx23r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )/(dx(0)**4)
   vyyyy23r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(1)**4)
   vzzzz23r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd))+(v(i1,i2,i3+2,kd)+v(i1,i2,i3-2,kd)) )/(dx(2)**4)
   vxxyy23r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)+v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))   +   (v(i1+1,i2+1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
   vxxzz23r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)+v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd))   +   (v(i1+1,i2,i3+1,kd)+v(i1-1,i2,i3+1,kd)+v(i1+1,i2,i3-1,kd)+v(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
   vyyzz23r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1,i2+1,i3,kd)  +v(i1,i2-1,i3,kd)+  v(i1,i2  ,i3+1,kd)+v(i1,i2  ,i3-1,kd))   +   (v(i1,i2+1,i3+1,kd)+v(i1,i2-1,i3+1,kd)+v(i1,i2+1,i3-1,kd)+v(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
   ! 3D laplacian squared = v.xxxx + v.yyyy + v.zzzz + 2 (v.xxyy + v.xxzz + v.yyzz )
   vLapSq23r(i1,i2,i3,kd)= ( 6.*v(i1,i2,i3,kd)   - 4.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd))    +(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*v(i1,i2,i3,kd)    -4.*(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))    +(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 6.*v(i1,i2,i3,kd)    -4.*(v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd))    +(v(i1,i2,i3+2,kd)+v(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 8.*v(i1,i2,i3,kd)     -4.*(v(i1+1,i2,i3,kd)  +v(i1-1,i2,i3,kd)  +v(i1  ,i2+1,i3,kd)+v(i1  ,i2-1,i3,kd))   +2.*(v(i1+1,i2+1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*v(i1,i2,i3,kd)     -4.*(v(i1+1,i2,i3,kd)  +v(i1-1,i2,i3,kd)  +v(i1  ,i2,i3+1,kd)+v(i1  ,i2,i3-1,kd))   +2.*(v(i1+1,i2,i3+1,kd)+v(i1-1,i2,i3+1,kd)+v(i1+1,i2,i3-1,kd)+v(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*v(i1,i2,i3,kd)     -4.*(v(i1,i2+1,i3,kd)  +v(i1,i2-1,i3,kd)  +v(i1,i2  ,i3+1,kd)+v(i1,i2  ,i3-1,kd))   +2.*(v(i1,i2+1,i3+1,kd)+v(i1,i2-1,i3+1,kd)+v(i1,i2+1,i3-1,kd)+v(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
 ! declare variables for getDerivatives macros
 !! turned off May 4, 2023
 !! #Include "../include/declareGetSixthDerivativesMacrosVariables.h"
 ! instead: 
! ****** File written by makeGetDerivativesMacros.maple  ******
real ur
real urr
real urrr
real urrrr
real us
real urs
real urrs
real urrrs
real uss
real urss
real urrss
real usss
real ursss
real ussss
real ut
real urt
real urrt
real urrrt
real ust
real urst
real urrst
real usst
real ursst
real ussst
real utt
real urtt
real urrtt
real ustt
real urstt
real usstt
real uttt
real urttt
real usttt
real utttt
real rxr
real rxrr
real rxrrr
real rxs
real rxrs
real rxrrs
real rxss
real rxrss
real rxsss
real rxt
real rxrt
real rxrrt
real rxst
real rxrst
real rxsst
real rxtt
real rxrtt
real rxstt
real rxttt
real ryr
real ryrr
real ryrrr
real rys
real ryrs
real ryrrs
real ryss
real ryrss
real rysss
real ryt
real ryrt
real ryrrt
real ryst
real ryrst
real rysst
real rytt
real ryrtt
real rystt
real ryttt
real sxr
real sxrr
real sxrrr
real sxs
real sxrs
real sxrrs
real sxss
real sxrss
real sxsss
real sxt
real sxrt
real sxrrt
real sxst
real sxrst
real sxsst
real sxtt
real sxrtt
real sxstt
real sxttt
real syr
real syrr
real syrrr
real sys
real syrs
real syrrs
real syss
real syrss
real sysss
real syt
real syrt
real syrrt
real syst
real syrst
real sysst
real sytt
real syrtt
real systt
real syttt
real rzr
real rzrr
real rzrrr
real rzs
real rzrs
real rzrrs
real rzss
real rzrss
real rzsss
real rzt
real rzrt
real rzrrt
real rzst
real rzrst
real rzsst
real rztt
real rzrtt
real rzstt
real rzttt
real szr
real szrr
real szrrr
real szs
real szrs
real szrrs
real szss
real szrss
real szsss
real szt
real szrt
real szrrt
real szst
real szrst
real szsst
real sztt
real szrtt
real szstt
real szttt
real txr
real txrr
real txrrr
real txs
real txrs
real txrrs
real txss
real txrss
real txsss
real txt
real txrt
real txrrt
real txst
real txrst
real txsst
real txtt
real txrtt
real txstt
real txttt
real tyr
real tyrr
real tyrrr
real tys
real tyrs
real tyrrs
real tyss
real tyrss
real tysss
real tyt
real tyrt
real tyrrt
real tyst
real tyrst
real tysst
real tytt
real tyrtt
real tystt
real tyttt
real tzr
real tzrr
real tzrrr
real tzs
real tzrs
real tzrrs
real tzss
real tzrss
real tzsss
real tzt
real tzrt
real tzrrt
real tzst
real tzrst
real tzsst
real tztt
real tzrtt
real tzstt
real tzttt
real rxi
real ryi
real sxi
real syi
real rzi
real szi
real txi
real tyi
real tzi
real rxx
real rxy
real rxz
real ryy
real ryz
real rzz
real sxx
real sxy
real sxz
real syy
real syz
real szz
real txx
real txy
real txz
real tyy
real tyz
real tzz
real rxxx
real rxxy
real rxyy
real rxxz
real rxyz
real rxzz
real ryyy
real ryyz
real ryzz
real rzzz
real sxxx
real sxxy
real sxyy
real sxxz
real sxyz
real sxzz
real syyy
real syyz
real syzz
real szzz
real txxx
real txxy
real txyy
real txxz
real txyz
real txzz
real tyyy
real tyyz
real tyzz
real tzzz
real rxxxx
real rxxxy
real rxxyy
real rxyyy
real rxxxz
real rxxyz
real rxyyz
real rxxzz
real rxyzz
real rxzzz
real ryyyy
real ryyyz
real ryyzz
real ryzzz
real rzzzz
real sxxxx
real sxxxy
real sxxyy
real sxyyy
real sxxxz
real sxxyz
real sxyyz
real sxxzz
real sxyzz
real sxzzz
real syyyy
real syyyz
real syyzz
real syzzz
real szzzz
real txxxx
real txxxy
real txxyy
real txyyy
real txxxz
real txxyz
real txyyz
real txxzz
real txyzz
real txzzz
real tyyyy
real tyyyz
real tyyzz
real tyzzz
real tzzzz
real uxx
real uxxxx
real uyy
real uxxyy
real uyyyy
real uzz
real uxxzz
real uyyzz
real uzzzz
  real ux,uy,uz
  real uxxx,uxxy,uxyy,uyyy,uxxz,uxzz,uzzz,uyyz,uyzz,uxyz
   ! 4th-order 1 sided derivative  extrap=(1 5 10 10 5 1)
   uxOneSided(i1,i2,i3,m)=-(10./3.)*u(i1,i2,i3,m)+6.*u(i1+is1,i2+is2,i3+is3,m)-2.*u(i1+2*is1,i2+2*is2,i3+2*is3,m)+(1./3.)*u(i1+3*is1,i2+3*is2,i3+3*is3,m)
   ! 2D laplacian squared = u.xxxx + 2 u.xxyy + u.yyyy
   lap2d2Pow2(i1,i2,i3,m)= ( 6.*u(i1,i2,i3,m)   - 4.*(u(i1+1,i2,i3,m)+u(i1-1,i2,i3,m))    +(u(i1+2,i2,i3,m)+u(i1-2,i2,i3,m)) )/(dx(0)**4) +( 6.*u(i1,i2,i3,m)    -4.*(u(i1,i2+1,i3,m)+u(i1,i2-1,i3,m))    +(u(i1,i2+2,i3,m)+u(i1,i2-2,i3,m)) )/(dx(1)**4)  +( 8.*u(i1,i2,i3,m)     -4.*(u(i1+1,i2,i3,m)+u(i1-1,i2,i3,m)+u(i1,i2+1,i3,m)+u(i1,i2-1,i3,m))   +2.*(u(i1+1,i2+1,i3,m)+u(i1-1,i2+1,i3,m)+u(i1+1,i2-1,i3,m)+u(i1-1,i2-1,i3,m)) )/( (dx(0)*dx(1))**2 )
   ! 3D laplacian squared = u.xxxx + u.yyyy + u.zzzz + 2 (u.xxyy + u.xxzz + u.yyzz )
   lap3d2Pow2(i1,i2,i3,m)= ( 6.*u(i1,i2,i3,m)   - 4.*(u(i1+1,i2,i3,m)+u(i1-1,i2,i3,m))    +(u(i1+2,i2,i3,m)+u(i1-2,i2,i3,m)) )/(dx(0)**4) +(  +6.*u(i1,i2,i3,m)    -4.*(u(i1,i2+1,i3,m)+u(i1,i2-1,i3,m))    +(u(i1,i2+2,i3,m)+u(i1,i2-2,i3,m)) )/(dx(1)**4)+(  +6.*u(i1,i2,i3,m)    -4.*(u(i1,i2,i3+1,m)+u(i1,i2,i3-1,m))    +(u(i1,i2,i3+2,m)+u(i1,i2,i3-2,m)) )/(dx(2)**4)+(8.*u(i1,i2,i3,m)     -4.*(u(i1+1,i2,i3,m)+u(i1-1,i2,i3,m)+u(i1,i2+1,i3,m)+u(i1,i2-1,i3,m))   +2.*(u(i1+1,i2+1,i3,m)+u(i1-1,i2+1,i3,m)+u(i1+1,i2-1,i3,m)+u(i1-1,i2-1,i3,m)) )/( (dx(0)*dx(1))**2 ) +(8.*u(i1,i2,i3,m)     -4.*(u(i1+1,i2,i3,m)+u(i1-1,i2,i3,m)+u(i1,i2,i3+1,m)+u(i1,i2,i3-1,m))   +2.*(u(i1+1,i2,i3+1,m)+u(i1-1,i2,i3+1,m)+u(i1+1,i2,i3-1,m)+u(i1-1,i2,i3-1,m)) )/( (dx(0)*dx(2))**2 ) +(8.*u(i1,i2,i3,m)     -4.*(u(i1,i2+1,i3,m)+u(i1,i2-1,i3,m)+u(i1,i2,i3+1,m)+u(i1,i2,i3-1,m))   +2.*(u(i1,i2+1,i3+1,m)+u(i1,i2-1,i3+1,m)+u(i1,i2+1,i3-1,m)+u(i1,i2-1,i3-1,m)) )/( (dx(1)*dx(2))**2 ) 
   ! ! Here is the the generic boundary condition forcing array. It uses the bcOffset(side,axis) values as an
   ! ! an offset from the bcf0 array to access the bcf10, bcf01, bcf11, ... arrays
   ! bcf(side,axis,i1,i2,i3,m) = bcf0(bcOffset(side,axis) + !     (i1-dim(0,0,side,axis)+(dim(1,0,side,axis)-dim(0,0,side,axis)+1)* !     (i2-dim(0,1,side,axis)+(dim(1,1,side,axis)-dim(0,1,side,axis)+1)* !     (i3-dim(0,2,side,axis)+(dim(1,2,side,axis)-dim(0,2,side,axis)+1)*(m)))))
   ! mixedRHS(component,side,axis,grid) = bcData(component+numberOfComponents*(0),side,axis,grid)
   ! mixedCoeff(component,side,axis,grid) = bcData(component+numberOfComponents*(1),side,axis,grid)
   ! mixedNormalCoeff(component,side,axis,grid) =  bcData(component+numberOfComponents*(2),side,axis,grid)
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
   !............... end statement functions
   ! if( .true. )then ! ********************* TESTING FOR TIMING
   !   return
   ! end if
   checkCoeff=0 ! 1 ! set to 1 to check coefficients in CBCs using discrete delta approach
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
   assignBCForImplicit             = ipar(16)
   bcApproach                      = ipar(17)
   numberOfFrequencies             = ipar(18)
   assignCornerGhostPoints         = ipar(19)
   orderOfExtrapolation            = ipar(20)
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
   c2 = c**2
   c4 = c**4
   c6 = c**6
   c8 = c**8
   twoPi = atan2(1.,1.)*8.; ! atan2(1,1)=pi/4
   if( gridType.eq.rectangular )then
     ! some macros want dr=dx for rectangular grids
     do axis=0,2
       dr(axis)=dx(axis)
     end do
   end if
   ! numberOfGhostPoints=orderOfAccuracy/2
   numberOfGhostPoints=numGhost ! now passed in 
   numGhost3          =numGhost ! num ghost in 3rd direction (i3)
   if( nd.eq.2 )then
     numGhost3=0
   end if
   if( assignCornerGhostPoints.ne.0 .and. assignCornerGhostPoints.ne.1 )then
     write(*,'("bcOptWave: ERROR: assignCornerGhostPoints=",i6," is unexpected")') assignCornerGhostPoints
     stop 4321
   end if
   ! write(*,'(" bcOptWave2dOrder2: dim=2, order=2")')
   ! *wdh* Nov 22, 2023 try turning of explicit BC's for implicit time-stepping
   if( .false. .and. gridIsImplicit.ne.0 .and. bcApproach==useCompatibilityBoundaryConditions .and. assignBCForImplicit==0 )then
     write(*,'(" bcOptWave: Skip explicit CBCs for implicit grid=",i4)') grid
     return
   end if
   if( t.le.3*dt .and. debug.gt.1 )then
   ! if( .true. )then
     write(*,'(" bcOptWave: nd=",i2," grid=",i4," gridType=",i2," orderOfAccuracy=",i2," uc=",i3," twilightZone=",i2)') nd,grid,gridType,orderOfAccuracy,uc,twilightZone
     write(*,'("  addForcingBC=",i4," forcingOption=",i4," assignKnownSolutionAtBoundaries=",i4)') addForcingBC, forcingOption, assignKnownSolutionAtBoundaries
     write(*,'("  t=",e10.2," dt=",e10.2," knownSolutionOption=",i4," REAL_MIN=",e10.2)') t,dt,knownSolutionOption,REAL_MIN
     write(*,'("  abcWave: c=",e14.6," cEM2=",e14.6)') c,cEM2
     write(*,'("  useUpwindDissipation=",i2," numGhost=",i2," assignCornerGhostPoints=",i2)') useUpwindDissipation,numGhost,assignCornerGhostPoints
     write(*,'("  assignBCForImplicit=",i4," bcApproach=",i4," gridIsImplicit=",i2)') assignBCForImplicit,bcApproach,gridIsImplicit
     write(*,'("  boundaryCondition=",6i4)') ((boundaryCondition(side,axis),side=0,1),axis=0,2)
   end if
   ! if( bcApproach.eq.useCompatibilityBoundaryConditions )then
   !   write(*,'("bcOptWave: ERROR: useCompatibilityBoundaryConditions not implemented yet.")') 
   !   stop 1010
   ! end if
   if( bcApproach.eq.useLocalCompatibilityBoundaryConditions )then
     write(*,'("bcOptWave: ERROR: useLocalCompatibilityBoundaryConditions not implemented yet.")') 
     stop 2020    
   end if
   ! ---- Make a local version of the boundaryCondition array ----
   do side=0,1
     do axis=0,2
       bc(side,axis)=boundaryCondition(side,axis)
       ! if( bc(side,axis).eq.absorbing )then
       !   ! DO THIS FOR NOW : use dirichlet for absorbing boundaries 
       !   bc(side,axis)=dirichlet
       ! end if
     end do
   end do
   if( assignKnownSolutionAtBoundaries.eq.1 ) then
     if( knownSolutionOption.eq.planeWave )then
       ! get parameter values from the C++ data-base
        ok=getReal(pdb,'ampPlaneWave',ampPlaneWave) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find ampPlaneWave")') 
          stop 1133
        end if
        ok=getReal(pdb,'kxPlaneWave',kxPlaneWave) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find kxPlaneWave")') 
          stop 1133
        end if
        ok=getReal(pdb,'kyPlaneWave',kyPlaneWave) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find kyPlaneWave")') 
          stop 1133
        end if
        ok=getReal(pdb,'kzPlaneWave',kzPlaneWave) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find kzPlaneWave")') 
          stop 1133
        end if
        ok=getReal(pdb,'omegaPlaneWave',omegaPlaneWave) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find omegaPlaneWave")') 
          stop 1133
        end if
        ok=getInt(pdb,'solveForScatteredField',solveForScatteredField) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getInt:ERROR: unable to find solveForScatteredField")') 
          stop 1122
        end if
       if( solveForScatteredField==1 )then
          ! If we solve for the scattered field then we flip the sign of the plane wave since this has been subtracted out
          ampPlaneWave = -ampPlaneWave
           ok=getInt(pdb,'solveHelmholtz',solveHelmholtz) 
           if( ok.eq.0 )then
             write(*,'("*** bcOptWave:getInt:ERROR: unable to find solveHelmholtz")') 
             stop 1122
           end if
          if( solveHelmholtz==1 )then 
             ! Get adjusted omega:
              ok=getReal(pdb,'omega',omega) 
              if( ok.eq.0 )then
                write(*,'("*** bcOptWave:getReal:ERROR: unable to find omega")') 
                stop 1133
              end if
             if( t.le.2*dt )then
               write(*,'(" bcOptWave:solveHelmholtz:scattering: Use adjusted omega=",1pe15.8," in place of omegaPlaneWave=",1pe15.8)') omega,omegaPlaneWave
             end if
             omegaPlaneWave = omega
          end if
          ! getIntParameter(solveHelmholtz)  
          ! if( solveHelmholtz==1 )then
          !    ! We are solving a Helmholtz problem : check that omega in the plane wave solution matches the omega for Helmholtz
          !    getRealParameter(omega)
          !    omegaTol = 1e-10  ! **FIX ME** use a multiple for REAL_EPSILON 
          !    if( abs(omega-omegaPlaneWave) .gt. omegaTol * omega )then
          !      write(*,'(" bcOptWave: solveForScatteredField=1, boundary forcing is a plane wave")') 
          !      write(*,'(" bcOptWave: ERROR: omegaPlaneWave=",e18.8," is not equal to omega(Helmholtz)=",e18.8)') omegaPlaneWave,omega
          !      stop 1234
          !    end if
          ! end if
       end if
       if(  t.le.2*dt .and. debug.gt.1  )then
         write(*,'(" bcOptWave:  knownSolutionOption=planeWave: solveForScatteredField=",i2," ampPlaneWave=",e10.2," kxPlaneWave=",e10.2," kyPlaneWave=",e10.2," omegaPlaneWave=",e14.4)') solveForScatteredField, ampPlaneWave,kxPlaneWave,kyPlaneWave,omegaPlaneWave
       end if 
       ! write(*,'(" bcOptWave:  knownSolutionOption=planeWave: solveForScatteredField=",i2," omegaPlaneWave=",e14.4)') solveForScatteredField, omegaPlaneWave
     else if( knownSolutionOption.eq.gaussianPlaneWave )then
       ! Get the parameters in the Gaussian plane wave (Set in userDefinedKnownSolution)
        ok=getReal(pdb,'kxGPW',kxGPW) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find kxGPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'kyGPW',kyGPW) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find kyGPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'kzGPW',kzGPW) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find kzGPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'x0GPW',x0GPW) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find x0GPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'y0GPW',y0GPW) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find y0GPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'z0GPW',z0GPW) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find z0GPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'k0GPW',k0GPW) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find k0GPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'betaGPW',betaGPW) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find betaGPW")') 
          stop 1133
        end if
       if(  t.le.dt .and. debug.ge.0  )then
         write(*,'(" bcOptWave:  knownSolutionOption=gaussianPlaneWave: kx,ky,kz=",3(1pe10.2)," x0,y0,z0=",3(1pe10.2)," k0,beta=",2(1pe10.2))') kxGPW,kyGPW,kzGPW,x0GPW,y0GPW,z0GPW,k0GPW,betaGPW
       end if           
     else if( knownSolutionOption.eq.boxHelmholtz )then
       ! get parameter values from the C++ data-base
        ok=getReal(pdb,'kxBoxHelmholtz',kxBoxHelmholtz) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find kxBoxHelmholtz")') 
          stop 1133
        end if
        ok=getReal(pdb,'kyBoxHelmholtz',kyBoxHelmholtz) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find kyBoxHelmholtz")') 
          stop 1133
        end if
        ok=getReal(pdb,'kzBoxHelmholtz',kzBoxHelmholtz) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find kzBoxHelmholtz")') 
          stop 1133
        end if
        ok=getReal(pdb,'omegaBoxHelmholtz',omegaBoxHelmholtz) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find omegaBoxHelmholtz")') 
          stop 1133
        end if
       coswt = cos(omegaBoxHelmholtz*t)
       assignKnownSolutionAtBoundaries=1  ! for inhomogeneous BCs
       if(  t.le.dt .and. debug.ge.1   )then
         write(*,'(" bcOptWave:  assignKnownSolutionAtBoundaries=",i4)') assignKnownSolutionAtBoundaries
         write(*,'(" bcOptWave:  numberOfFrequencies=",i4)') numberOfFrequencies
         write(*,'(" bcOptWave:  frequencyArray=",10(1pe12.4,1x))') (frequencyArray(freq),freq=0,numberOfFrequencies-1)
         write(*,'(" bcOptWave:  knownSolutionOption=boxHelmholtz: kx,ky,kz,omega=",4e10.2)') kxBoxHelmholtz,kyBoxHelmholtz,kzBoxHelmholtz,omegaBoxHelmholtz
       end if
     else if( knownSolutionOption.eq.polyPeriodic )then
       ! get parameter values from the C++ data-base
        ok=getReal(pdb,'omegaPolyPeriodic',omegaPolyPeriodic) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find omegaPolyPeriodic")') 
          stop 1133
        end if
        ok=getReal(pdb,'a0PolyPeriodic',a0PolyPeriodic) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find a0PolyPeriodic")') 
          stop 1133
        end if
        ok=getReal(pdb,'a1PolyPeriodic',a1PolyPeriodic) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find a1PolyPeriodic")') 
          stop 1133
        end if
        ok=getReal(pdb,'b1PolyPeriodic',b1PolyPeriodic) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find b1PolyPeriodic")') 
          stop 1133
        end if
        ok=getReal(pdb,'c1PolyPeriodic',c1PolyPeriodic) 
        if( ok.eq.0 )then
          write(*,'("*** bcOptWave:getReal:ERROR: unable to find c1PolyPeriodic")') 
          stop 1133
        end if
       coswt = cos(omegaPolyPeriodic*t)
       if(  t.le.dt .and. debug.gt.1   )then
         write(*,'(" bcOptWave:  knownSolutionOption=polyPeriodic: a0,a1,b1,c1,omega=",5e10.2)') a0PolyPeriodic,a1PolyPeriodic,b1PolyPeriodic,c1PolyPeriodic,omegaPolyPeriodic
       end if
     else if( knownSolutionOption.ne.0 .and. knownSolutionOption.ne.otherKnownSolution )then
       write(*,'("bcOptWave:ERROR: unknown knownSolutionOption=",i6)') knownSolutionOption
       stop 1111
     end if 
   end if
   ! TEST: 
   ! getRealParameter(omega)
   ! getRealParameter(cfl)
   ! write(*,'(" bcOptWave:  cfl=",e10.2)') cfl
   if( uc.lt.0 .or. uc.ge.numberOfComponents )then
     write(*,'("bcOptWave:ERROR: invalid uc=",i6," but numberOfComponents=",i3)')  uc,numberOfComponents
     stop 1111
   end if
   epsx=REAL_MIN*100.  ! for normal
   if( orderOfAccuracy.ne.2 .and. orderOfAccuracy.ne.4 .and. orderOfAccuracy.ne.6 .and. orderOfAccuracy.ne.8 )then
     write(*,'("bcOptWave:ERROR: orderOfAccuracy is not 2, 4 or 6, orderOfAccuracy=",i4)') orderOfAccuracy
     stop 1111
   end if
   ! Now passed in: 
   ! numGhost=orderOfAccuracy/2
   if( assignBCForImplicit.eq.1 .or. assignBCForImplicit.eq.2 )then
     ! -------- IMPLICIT BoundaryConditions --------
       ! -------- IMPLICIT BoundaryConditions --------
       ! 
       !   assignBCForImplicit = 1 : BCs for implicit time stepping
       !   assignBCForImplicit = 2 : BCs for direct Helmholtz solve
       ! if( .true. )then
       !   write(*,'("bcOptWave: fill BCs into RHS for implicit solver")')
       ! end if
       ! write(*,'("FINISH ME")')
       ! stop 6789
       an3=0.
        extra1a=numGhost
        extra1b=numGhost
        extra2a=numGhost
        extra2b=numGhost
        if( nd.eq.3 )then
          extra3a=numGhost
          extra3b=numGhost
        else
          extra3a=0
          extra3b=0
        end if
        if( bc(0,0).lt.0 )then
          extra1a=max(0,extra1a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
        else if( bc(0,0).eq.0 )then
          extra1a=numGhost  ! include interpolation points since we assign ghost points outside these
        end if
        ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
        if( bc(1,0).lt.0 )then
          extra1b=max(0,extra1b) ! over-ride numGhost=-1 : assign ends in periodic directions
        else if( bc(1,0).eq.0 )then
          extra1b=numGhost
        end if
        if( bc(0,1).lt.0 )then
          extra2a=max(0,extra2a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
        else if( bc(0,1).eq.0 )then
          extra2a=numGhost  ! include interpolation points since we assign ghost points outside these
        end if
        ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
        if( bc(1,1).lt.0 )then
          extra2b=max(0,extra2b) ! over-ride numGhost=-1 : assign ends in periodic directions
        else if( bc(1,1).eq.0 )then
          extra2b=numGhost
        end if
        if(  nd.eq.3 )then
         if( bc(0,2).lt.0 )then
           extra3a=max(0,extra3a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
         else if( bc(0,2).eq.0 )then
           extra3a=numGhost  ! include interpolation points since we assign ghost points outside these
         end if
         ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
         if( bc(1,2).lt.0 )then
           extra3b=max(0,extra3b) ! over-ride numGhost=-1 : assign ends in periodic directions
         else if( bc(1,2).eq.0 )then
           extra3b=numGhost
         end if
        end if
        do axis=0,nd-1
        do side=0,1
          if( .true. .or. bc(side,axis).gt.0 )then ! we may set ghost outside interp for implicit
            ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,bc(side,axis)
            n1a=gridIndexRange(0,0)
            n1b=gridIndexRange(1,0)
            n2a=gridIndexRange(0,1)
            n2b=gridIndexRange(1,1)
            n3a=gridIndexRange(0,2)
            n3b=gridIndexRange(1,2)
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
            nn1a=gridIndexRange(0,0)-extra1a
            nn1b=gridIndexRange(1,0)+extra1b
            nn2a=gridIndexRange(0,1)-extra2a
            nn2b=gridIndexRange(1,1)+extra2b
            nn3a=gridIndexRange(0,2)-extra3a
            nn3b=gridIndexRange(1,2)+extra3b
            if( axis.eq.0 )then
              nn1a=gridIndexRange(side,axis)
              nn1b=gridIndexRange(side,axis)
            else if( axis.eq.1 )then
              nn2a=gridIndexRange(side,axis)
              nn2b=gridIndexRange(side,axis)
            else
              nn3a=gridIndexRange(side,axis)
              nn3b=gridIndexRange(side,axis)
            end if
            is=1-2*side
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
            axisp1=mod(axis+1,nd)
            axisp2=mod(axis+2,nd)
            i3=n3a
            if( debug.gt.31 )then
              write(*,'(" bcOptWave: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i5)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
              write(*,'(" bcOptWave: numGhost,extra1a,extra2a,extra3a=",4i4,", loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i5)') numGhost,extra1a,extra2a,extra3a,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
            end if
          end if ! if bc>0 
          assignTwilightZone=twilightZone
         if( bc(side,axis) == dirichlet )then
            ! ============== IMPLICIT DIRICHLET =============
           ff=0.
           if( bcApproach.eq.useOneSidedBoundaryConditions .or. orderOfAccuracy.eq.2 )then
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
               if( mask(i1,i2,i3).ne.0 )then
                   if( assignTwilightZone.eq.1 )then
                     ! compute RHS from TZ
                     if( nd.eq.2 )then
                       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                     else
                       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                     end if
                     ff = ue
                   else if( assignKnownSolutionAtBoundaries.eq.1 )then
                     ! -- we set inhomogeneous Dirichlet values for some known solutions 
                     if( knownSolutionOption.eq.planeWave )then
                       ! --- evaluate the plane wave solution ---
                       if( nd.eq.2 )then
                         ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                       else
                         ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                       end if 
                     else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                       ! Eval the Gaussian plane wave solution
                       !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                       !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                       !  
                       if( nd.eq.2 )then
                         xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) - c*t
                       else
                         xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) + kzGPW*(xy(i1,i2,i3,2)-z0GPW) - c*t
                       end if 
                       ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                     else if( knownSolutionOption.eq.boxHelmholtz ) then
                       ! --- evaluate the boxHelmholtz solution ---
                       ! For multi-freq we add up all the component frequencies
                       ff = 0. 
                       do freq=0,numberOfFrequencies-1
                         ! kx = kxBoxHelmholtz + twoPi*freq
                         ! ky = kyBoxHelmholtz + twoPi*freq
                         ! kz = kzBoxHelmholtz + twoPi*freq
                         ! coswt = cos( frequencyArray(freq)*t )
                             ! This macro is used in bcOptWave.bf90
                             omega = frequencyArray(freq);
                             kx = kxBoxHelmholtz*(freq*.5+1.)
                             ky = kyBoxHelmholtz*(freq*.5+1.)
                             kz = kzBoxHelmholtz*(freq*.5+1.)
                             ! kx = kxBoxHelmholtz + twoPi*freq.
                             ! ky = kyBoxHelmholtz + twoPi*freq.
                             ! kz = kzBoxHelmholtz + twoPi*freq.    
                         coswt = cos( omega*t )
                         if( nd.eq.2 )then
                           ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) *coswt
                         else
                           ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) * sin( kz*xy(i1,i2,i3,2) ) *coswt
                         end if 
                       end do
                     else if( knownSolutionOption.eq.polyPeriodic ) then
                       ! --- evaluate the polyPeriodic solution ---
                       if( nd.eq.2 )then
                         ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1)                                 ) *coswt
                       else
                         ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                       end if 
                     else
                       stop 9876
                     end if 
                   end if
                 ! fill in boundary value: 
                 u(i1,i2,i3,uc)=ff
                 ! -- Set ghost to zero (RHS to extrapolation conditions) ---
                 ! Is this necessary ?
                 do ghost=1,numGhost
                   j1=i1-is1*ghost
                   j2=i2-is2*ghost
                   j3=i3-is3*ghost
                   u(j1,j2,j3,uc) = 0.
                 end do  
               end if 
              end do
              end do
              end do
           else if( bcApproach.eq.useCompatibilityBoundaryConditions )then
             ! if( t.le.5*dt )then
             !   write(*,'("bcOpt: implicit: finish me for CBC ***")') 
             ! end if
             ! assert orderOfAccuracy==4 ************************************
             if( orderOfAccuracy.ne.4 )then
               write(*,'("bcOpt:ERROR:Dirichlet: implicit, CBC, but orderOfAccuracy=",i4)') orderOfAccuracy
               stop 1234
             end if
             ! At a Dirichlet-Dirichlet Corner we cannot use a CBCs on both sides since this is a duplicate
             !         |
             !   E--C--+
             !         |
             !   D--D--X--+--+--+--
             !         |  |  |
             !         D  C  C
             !         |  |  |
             !         D  E  E
             ! D = use dirichlet BC on extended ghost line(s)
             ! C = use CBC
             ! E = extrap
             ! Assign extended boundary 
             m1a=n1a-numberOfGhostPoints; m1b=n1b+numberOfGhostPoints; 
             m2a=n2a-numberOfGhostPoints; m2b=n2b+numberOfGhostPoints; 
             m3a=n3a-numberOfGhostPoints; m3b=n3b+numberOfGhostPoints; 
             if( nd.eq.2 )then
               m3a=n3a; m3b=n3b;
             end if
             if( axis.eq.0 )then
               m1a=n1a; m1b=n1b;
             else if( axis.eq.1 )then
               m2a=n2a; m2b=n2b;
             else
               m3a=n3a; m3b=n3b;
             end if
             ! write(*,'("bcOpt:Imp:Dirichlet:RHS side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i5)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
             do i3=m3a,m3b
             do i2=m2a,m2b
             do i1=m1a,m1b
               if( mask(i1,i2,i3).ne.0 )then
                   if( assignTwilightZone.eq.1 )then
                     ! compute RHS from TZ
                     if( nd.eq.2 )then
                       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                     else
                       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                     end if
                     ff = ue
                   else if( assignKnownSolutionAtBoundaries.eq.1 )then
                     ! -- we set inhomogeneous Dirichlet values for some known solutions 
                     if( knownSolutionOption.eq.planeWave )then
                       ! --- evaluate the plane wave solution ---
                       if( nd.eq.2 )then
                         ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                       else
                         ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                       end if 
                     else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                       ! Eval the Gaussian plane wave solution
                       !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                       !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                       !  
                       if( nd.eq.2 )then
                         xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) - c*t
                       else
                         xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) + kzGPW*(xy(i1,i2,i3,2)-z0GPW) - c*t
                       end if 
                       ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                     else if( knownSolutionOption.eq.boxHelmholtz ) then
                       ! --- evaluate the boxHelmholtz solution ---
                       ! For multi-freq we add up all the component frequencies
                       ff = 0. 
                       do freq=0,numberOfFrequencies-1
                         ! kx = kxBoxHelmholtz + twoPi*freq
                         ! ky = kyBoxHelmholtz + twoPi*freq
                         ! kz = kzBoxHelmholtz + twoPi*freq
                         ! coswt = cos( frequencyArray(freq)*t )
                             ! This macro is used in bcOptWave.bf90
                             omega = frequencyArray(freq);
                             kx = kxBoxHelmholtz*(freq*.5+1.)
                             ky = kyBoxHelmholtz*(freq*.5+1.)
                             kz = kzBoxHelmholtz*(freq*.5+1.)
                             ! kx = kxBoxHelmholtz + twoPi*freq.
                             ! ky = kyBoxHelmholtz + twoPi*freq.
                             ! kz = kzBoxHelmholtz + twoPi*freq.    
                         coswt = cos( omega*t )
                         if( nd.eq.2 )then
                           ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) *coswt
                         else
                           ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) * sin( kz*xy(i1,i2,i3,2) ) *coswt
                         end if 
                       end do
                     else if( knownSolutionOption.eq.polyPeriodic ) then
                       ! --- evaluate the polyPeriodic solution ---
                       if( nd.eq.2 )then
                         ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1)                                 ) *coswt
                       else
                         ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                       end if 
                     else
                       stop 9876
                     end if 
                   end if
                 ! fill in boundary value: 
                 u(i1,i2,i3,uc)=ff
                 ! write(*,'("IMP-RHS CBC: i1,i2=",2i4," u=",e10.2)') i1,i2,ff 
               end if
             end do
             end do
             end do
             ! --- Adjust loops for CBC RHS at Dirichlet-Dirichlet corners ---
             iab(0,0)=n1a; iab(1,0)=n1b; 
             iab(0,1)=n2a; iab(1,1)=n2b; 
             iab(0,2)=n3a; iab(1,2)=n3b;
             do sidea=0,1  ! adjacent side
               if( bc(sidea,axisp1) == dirichlet )then
                 iab(sidea,axisp1)=iab(sidea,axisp1) + 1-2*sidea  ! avoid end-point
               end if 
             end do
             if( nd.eq.3 )then
               do sidea=0,1  ! adjacent side
                 if( bc(sidea,axisp2) == dirichlet )then
                   iab(sidea,axisp2)=iab(sidea,axisp2) + 1-2*sidea  ! avoid end-point
                 end if 
               end do          
             endif
             m1a=iab(0,0); m1b=iab(1,0); 
             m2a=iab(0,1); m2b=iab(1,1); 
             m3a=iab(0,2); m3b=iab(1,2);
             ! write(*,'("bcOpt:Imp:Dirichlet:CBC2:RHS side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i5)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
             ! -- assign RHS to CBC and extrapolation --
             do i3=m3a,m3b
             do i2=m2a,m2b
             do i1=m1a,m1b
               if( mask(i1,i2,i3).ne.0 )then
                 ! getDirichletForcing(ff,i1,i2,i3)
                 ! ! fill in boundary value: 
                 ! u(i1,i2,i3,uc)=ff
                 j1=i1-is1
                 j2=i2-is2
                 j3=i3-is3 
                 if( assignTwilightZone.eq.1 )then
                   ! compute RHS from TZ
                   if( nd.eq.2 )then
                     call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexx )
                     call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyy )
                     ueLap = uexx + ueyy                  
                   else
                    call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexx )
                    call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyy )
                    call ogDeriv(ep,0,0,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezz )
                    ueLap = uexx + ueyy + uezz                  
                   end if 
                   ! write(*,'("IMP-RHS CBC: j1,j2=",2i4," ueLap=",e10.2)') j1,j2,ueLap 
                   u(j1,j2,j3,uc)=c2*ueLap
                   do ghost=2,numGhost
                     j1=i1-is1*ghost
                     j2=i2-is2*ghost
                     j3=i3-is3*ghost
                     u(j1,j2,j3,uc) = 0.
                   end do                 
                 else
                   ! -- FIX ME FOR PLANE WAVE ETC ---
                   ! Is this necessary ?
                   do ghost=1,numGhost
                     j1=i1-is1*ghost
                     j2=i2-is2*ghost
                     j3=i3-is3*ghost
                     u(j1,j2,j3,uc) = 0.
                   end do 
                 end if 
               end if 
              end do
              end do
              end do
           else
             write(*,'("bcOpt: implicit: unexpected bcApproach=",i6)') bcApproach
             stop 1111
           end if
         else if( bc(side,axis) == exactBC )then
           ! ------------------------------------------------------------------------
           ! NOTE: KnownSolution values are assigned in **applyBoundaryConditions**
           ! ------------------------------------------------------------------------
           if( assignTwilightZone==1 )then
             if( t.le.3*dt )then
               write(*,'("bcOptWave: Assign IMP EXACT BC : Twilightzone grid,side,axis=",3i3)') grid,side,axis 
             end if        
             ! --- Assign extended boundary AND ghost points ----
             m1a=n1a-numberOfGhostPoints; m1b=n1b+numberOfGhostPoints; 
             m2a=n2a-numberOfGhostPoints; m2b=n2b+numberOfGhostPoints; 
             m3a=n3a-numberOfGhostPoints; m3b=n3b+numberOfGhostPoints; 
             if( nd.eq.2 )then
               m3a=n3a; m3b=n3b;
             end if
             if( axis.eq.0 )then
               if( side==0 )then
                 m1a=n1a-numberOfGhostPoints 
                 m1b=n1a
               else
                 m1a=n1b
                 m1b=n1b+numberOfGhostPoints
               end if
             else if( axis.eq.1 )then
              if( side==0 )then
                 m2a=n2a-numberOfGhostPoints 
                 m2b=n2a
               else
                 m2a=n2b
                 m2b=n2b+numberOfGhostPoints
               end if          
             else
              if( side==0 )then
                 m3a=n3a-numberOfGhostPoints 
                 m3b=n3a
               else
                 m3a=n3b
                 m3b=n3b+numberOfGhostPoints
               end if          
             end if
             do i3=m3a,m3b
             do i2=m2a,m2b
             do i1=m1a,m1b
               if( mask(i1,i2,i3).ne.0 )then
                   if( assignTwilightZone.eq.1 )then
                     ! compute RHS from TZ
                     if( nd.eq.2 )then
                       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                     else
                       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                     end if
                     ff = ue
                   else if( assignKnownSolutionAtBoundaries.eq.1 )then
                     ! -- we set inhomogeneous Dirichlet values for some known solutions 
                     if( knownSolutionOption.eq.planeWave )then
                       ! --- evaluate the plane wave solution ---
                       if( nd.eq.2 )then
                         ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                       else
                         ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                       end if 
                     else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                       ! Eval the Gaussian plane wave solution
                       !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                       !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                       !  
                       if( nd.eq.2 )then
                         xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) - c*t
                       else
                         xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) + kzGPW*(xy(i1,i2,i3,2)-z0GPW) - c*t
                       end if 
                       ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                     else if( knownSolutionOption.eq.boxHelmholtz ) then
                       ! --- evaluate the boxHelmholtz solution ---
                       ! For multi-freq we add up all the component frequencies
                       ff = 0. 
                       do freq=0,numberOfFrequencies-1
                         ! kx = kxBoxHelmholtz + twoPi*freq
                         ! ky = kyBoxHelmholtz + twoPi*freq
                         ! kz = kzBoxHelmholtz + twoPi*freq
                         ! coswt = cos( frequencyArray(freq)*t )
                             ! This macro is used in bcOptWave.bf90
                             omega = frequencyArray(freq);
                             kx = kxBoxHelmholtz*(freq*.5+1.)
                             ky = kyBoxHelmholtz*(freq*.5+1.)
                             kz = kzBoxHelmholtz*(freq*.5+1.)
                             ! kx = kxBoxHelmholtz + twoPi*freq.
                             ! ky = kyBoxHelmholtz + twoPi*freq.
                             ! kz = kzBoxHelmholtz + twoPi*freq.    
                         coswt = cos( omega*t )
                         if( nd.eq.2 )then
                           ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) *coswt
                         else
                           ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) * sin( kz*xy(i1,i2,i3,2) ) *coswt
                         end if 
                       end do
                     else if( knownSolutionOption.eq.polyPeriodic ) then
                       ! --- evaluate the polyPeriodic solution ---
                       if( nd.eq.2 )then
                         ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1)                                 ) *coswt
                       else
                         ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                       end if 
                     else
                       stop 9876
                     end if 
                   end if
                 ! fill in boundary value: 
                 u(i1,i2,i3,uc)=ff
                 ! write(*,'("IMP-RHS CBC: i1,i2=",2i4," u=",e10.2)') i1,i2,ff 
               end if
             end do
             end do
             end do
           end if ! end if assignTwilightZone
         else if( bc(side,axis) == neumann )then
           if( gridType.eq.rectangular )then
             ! compute the outward normal (an1,an2,an3)
             an1 = 0.
             an2 = 0.
             an3 = 0.
             if( axis.eq.0 )then
              an1=-is
             else if( axis.eq.1 )then
              an2=-is
             else
              an3=-is
             end if
           end if        
           ! BC is a0*u + a1*u.n = 
           a0=0.
           a1=1.
           ff=0.
           if( bcApproach==useCompatibilityBoundaryConditions .and. orderOfAccuracy==4 )then
             ! Exclude the ends when adjacent boundaries are dirichlet or exact: *wdh* Dec 4, 2023
             extram=0
               ! Assume sides with axis=0 are adjacent sides (fixed later)
               if( bc(0,0)==dirichlet .or. bc(0,0)==exactBC )then
                 m1a= gridIndexRange(0,0)+1
               else
                 m1a=gridIndexRange(0,0)-extram
               end if
               if( bc(1,0)==dirichlet .or. bc(1,0)==exactBC )then
                 m1b = gridIndexRange(1,0)-1
               else
                 m1b = gridIndexRange(1,0)+extram
               end if
               ! Assume sides with axis=1 are adjacent sides (fixed later)
               if( bc(0,1)==dirichlet .or. bc(0,1)==exactBC )then
                 m2a= gridIndexRange(0,1)+1
               else
                 m2a=gridIndexRange(0,1)-extram
               end if
               if( bc(1,1)==dirichlet .or. bc(1,1)==exactBC )then
                 m2b = gridIndexRange(1,1)-1
               else
                 m2b = gridIndexRange(1,1)+extram
               end if  
               if( nd.eq.2 )then
                 m3a=gridIndexRange(0,2)
                 m3b=gridIndexRange(1,2)
               else
                 if( bc(0,2)==dirichlet .or. bc(0,2)==exactBC )then
                   m3a= gridIndexRange(0,2)+1
                 else
                   m3a=gridIndexRange(0,2)-extram
                 end if
                 if( bc(1,2)==dirichlet .or. bc(1,2)==exactBC )then
                   m3b = gridIndexRange(1,2)-1
                 else
                   m3b = gridIndexRange(1,2)+extram
                 end if
               end if
               if( axis.eq.0 )then
                 m1a=gridIndexRange(side,axis)
                 m1b=gridIndexRange(side,axis)
               else if( axis.eq.1 )then
                 m2a=gridIndexRange(side,axis)
                 m2b=gridIndexRange(side,axis)
               else
                 m3a=gridIndexRange(side,axis)
                 m3b=gridIndexRange(side,axis)
               end if
           else
             m1a=n1a; m1b=n1b; m2a=n2a; m2b=n2b; m3a=n3a; m3b=n3b; 
           end if
           ! write(*,'("bcOpt:Imp:Neumann:CBC:RHS side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i5)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
           do i3=m3a,m3b
           do i2=m2a,m2b
           do i1=m1a,m1b
           ! beginLoops3d()
             if( mask(i1,i2,i3).ne.0 )then
               if( gridType.eq.curvilinear )then
                 ! compute the outward normal (an1,an2,an3)
                     an1 = rsxy(i1,i2,i3,axis,0)
                     an2 = rsxy(i1,i2,i3,axis,1)
                     if( nd.eq.2 )then
                      aNormi = (-is)/sqrt(an1**2+an2**2)
                      an1=an1*aNormi
                      an2=an2*aNormi
                     else
                      an3 = rsxy(i1,i2,i3,axis,2)
                      aNormi = (-is)/sqrt(an1**2+an2**2+an3**2)
                      an1=an1*aNormi
                      an2=an2*aNormi
                      an3=an3*aNormi
                     end if
               end if            
                 if( assignTwilightZone.eq.1 )then
                   ! compute RHS from TZ
                   if( nd.eq.2 )then
                     call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                     call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uex)
                     call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uey)
                     ff = a0*ue + a1*( an1*uex + an2*uey )
                   else
                     call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                     call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uex)
                     call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uey)
                     call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uez)
                     ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
                   end if
                 else if( assignKnownSolutionAtBoundaries.eq.1 )then
                   ! -- we set inhomogeneous Neumann values for some known solutions 
                   if( knownSolutionOption.eq.planeWave )then
                     ! --- evaluate RHS for the plane wave solution ---
                     if( nd.eq.2 )then
                       ue    = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                       cosPW = ampPlaneWave*cos( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                       uex   = kxPlaneWave*cosPW
                       uey   = kyPlaneWave*cosPw
                       ff = a0*ue + a1*( an1*uex + an2*uey )
                     else
                       ue    = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                       cosPW = ampPlaneWave*cos( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                       uex   = kxPlaneWave*cosPW
                       uey   = kyPlaneWave*cosPw
                       uez   = kzPlaneWave*cosPw
                       ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
                     end if 
                   else if( knownSolutionOption.eq.gaussianPlaneWave )then
                     ! Do nothing for Gaussian plane wave solution for now
                     ff = 0.
                   else if( knownSolutionOption.eq.boxHelmholtz ) then
                     ! --- evaluate RHS the boxHelmholtz solution ---
                     if( nd.eq.2 )then
                       ue  = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) * sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) *coswt
                       uex = cos( kxBoxHelmholtz*xy(i1,i2,i3,0) ) * sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) *coswt * kxBoxHelmholtz
                       uey = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) * cos( kyBoxHelmholtz*xy(i1,i2,i3,1) ) *coswt * kyBoxHelmholtz
                       ff = a0*ue + a1*( an1*uex + an2*uey )
                     else
                       ue  = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * sin( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt
                       uex = cos( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * sin( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt * kxBoxHelmholtz
                       uey = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *cos( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * sin( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt * kyBoxHelmholtz
                       uez = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * cos( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt * kzBoxHelmholtz
                       ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
                     end if
                   else if( knownSolutionOption.eq.polyPeriodic ) then
                     ! --- evaluate RHS the polyPeriodic solution ---
                     if( nd.eq.2 )then
                       ue  = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) ) *coswt
                       uex = (      a1PolyPeriodic                                                            ) *coswt
                       uey = (                          b1PolyPeriodic                                        ) *coswt
                       ff = a0*ue + a1*( an1*uex + an2*uey )
                     else
                       ue  = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                       uex = (      a1PolyPeriodic                                                                                            ) *coswt
                       uey = (                          b1PolyPeriodic                                                                        ) *coswt
                       uez = (                                              c1PolyPeriodic                                                    ) *coswt
                       ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
                     end if
                   else
                     stop 9876
                   end if 
                 end if
               ! fill in first ghost:
               j1=i1-is1
               j2=i2-is2
               j3=i3-is3
               u(j1,j2,j3,uc)=ff
               ! write(*,'("IMP-RHS NEUMANN: CBC1: j1,j2,j3=",3i4," an1,an2,an3=",3(1pe10.2)," rhs=",e10.2)') j1,j2,j3,an1,an2,an3,ff        
             end if 
            end do
            end do
            end do
           if( bcApproach.eq.useCompatibilityBoundaryConditions .and. orderOfAccuracy.ge.4 )then   
             if( orderOfAccuracy.ne.4 )then
               write(*,'("bcOpt:ERROR:Neumann implicit, CBC, but orderOfAccuracy=",i4)') orderOfAccuracy
               stop 3456
             end if      
             ! if( t.le.5*dt )then
             !   write(*,'("bcOpt: implicit:Neumann: Assign RHS for CBC ***")') 
             !   ! stop 8888
             ! end if
             ! -- assign RHS to CBC  --
             !   Neumann:   u.n = g 
             !   CBC        (u_tt).n = g_tt ->   c^2 Delta( u.n ) = g.tt 
             !
             an1 = -is1; an2 = -is2; an3 = -is3; ! outward normal
             ! beginLoops3d()
             do i3=m3a,m3b
             do i2=m2a,m2b
             do i1=m1a,m1b
               if( mask(i1,i2,i3).ne.0 )then
                 j1=i1-2*is1 ! 2nd ghost 
                 j2=i2-2*is2
                 j3=i3-2*is3 
                 if( gridType.eq.curvilinear )then
                   ! compute the outward normal (an1,an2,an3)
                       an1 = rsxy(i1,i2,i3,axis,0)
                       an2 = rsxy(i1,i2,i3,axis,1)
                       if( nd.eq.2 )then
                        aNormi = (-is)/sqrt(an1**2+an2**2)
                        an1=an1*aNormi
                        an2=an2*aNormi
                       else
                        an3 = rsxy(i1,i2,i3,axis,2)
                        aNormi = (-is)/sqrt(an1**2+an2**2+an3**2)
                        an1=an1*aNormi
                        an2=an2*aNormi
                        an3=an3*aNormi
                       end if
                 end if                    
                 if( assignTwilightZone.eq.1 )then
                   ! compute RHS from TZ
                   if( nd.eq.2 )then
                     call ogDeriv(ep,0,3,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxx )
                     call ogDeriv(ep,0,1,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexyy )
                     call ogDeriv(ep,0,2,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxy )
                     call ogDeriv(ep,0,0,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyy )
                     rhs = an1*(uexxx + uexyy) + an2*(uexxy + ueyyy)
                     ! test: rhs = an1*(uexxx) + an2*(ueyyy)
                   else
                     call ogDeriv(ep,0,3,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxx )
                     call ogDeriv(ep,0,1,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexyy )
                     call ogDeriv(ep,0,1,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexzz )
                     call ogDeriv(ep,0,2,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxy )
                     call ogDeriv(ep,0,0,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyy )
                     call ogDeriv(ep,0,0,1,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyzz )
                     call ogDeriv(ep,0,2,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxz )
                     call ogDeriv(ep,0,0,2,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyz )
                     call ogDeriv(ep,0,0,0,3,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezzz )
                     rhs = an1*(uexxx+uexyy+uexzz) + an2*(uexxy+ueyyy+ueyzz) + an3*(uexxz+ueyyz+uezzz)
                     ! test: rhs = an1*(uexxx) + an2*(ueyyy) + an3*(uezzz)
                   end if 
                   ! write(*,'("IMP-RHS NEUMANN: CBC2: j1,j2,j3=",3i4," an1,an2,an3=",3(1pe10.2)," rhs=",e10.2)') j1,j2,j3,an1,an2,an3,rhs
                   u(j1,j2,j3,uc)=c2*rhs
                 else
                   ! -- FIX ME FOR PLANE WAVE ETC ---
                   u(j1,j2,j3,uc) = 0.
                 end if 
               end if 
              end do
              end do
              end do
           end if       
         else if( bc(side,axis) == absorbing .or. bc(side,axis) == abcEM2 )then
           ! if( t<=4.*dt )then
           !   write(*,*) "bcOpt: implicit BC for absorbing/EM2 "
           ! end if
           ! stop 4444
           !  Use adjusted c for the EM2 absorbing BC to account for time-discretization errors
           !    D+t (Dx ) w + A+( ... )
           ! cEM2 = c*tan(frequencyArray(0)*dt/2.)/(frequencyArraySave(0)*dt/2.);
           ca = cEM2;  ! Adjusted c 
           if( gridType.eq.rectangular )then
              do i3=n3a,n3b
              do i2=n2a,n2b
              do i1=n1a,n1b
               if( mask(i1,i2,i3).ne.0 )then
                 j1  = i1-is1; j2  = i2-is2; j3  = i3-is3;     ! ghost 
                 i1p = i1+is1; i2p = i2+is2; i3p = i3+is3;     ! first line inside
                 ! We need current solution un here
                 ! res = -is*(unx-ucx)/dt + (.5*c)*( unxx + ucxx) + (.25*c)*( unyy + ucyy );
                 if( assignBCForImplicit.eq.1 )then
                   if( axis==0 )then
                     u(j1,j2,j3,uc) = (un(j1,j2,j3,uc)-un(i1p,i2p,i3p,uc))/(2.*dx(axis)*dt)                    - .5*(   ca*(un(i1+1,i2,i3,uc)-2.*un(i1,i2,i3,uc)+un(i1-1,i2,i3,uc))/(dx(0)**2) )   - .5*(.5*ca*(un(i1,i2+1,i3,uc)-2.*un(i1,i2,i3,uc)+un(i1,i2-1,i3,uc))/(dx(1)**2) ) 
                   else
                     u(j1,j2,j3,uc) = (un(j1,j2,j3,uc)-un(i1p,i2p,i3p,uc))/(2.*dx(axis)*dt)                     - .5*( .5*ca*(un(i1+1,i2,i3,uc)-2.*un(i1,i2,i3,uc)+un(i1-1,i2,i3,uc))/(dx(0)**2) )   - .5*(    ca*(un(i1,i2+1,i3,uc)-2.*un(i1,i2,i3,uc)+un(i1,i2-1,i3,uc))/(dx(1)**2) )                   
                   end if 
                 else
                    ! Bc for direct Helmholtz solve
                    u(j1,j2,j3,uc) = 0.
                 end if
                 do ghost=2,numGhost
                   j1=i1-is1*ghost
                   j2=i2-is2*ghost
                   j3=i3-is3*ghost
                   u(j1,j2,j3,uc) = 0.
                 end do                
               end if
              end do
              end do
              end do
           else
            write(*,*) "bcOpt: implicit BC for absorbing/EM2 -- CURVILINEAR : finish me"
             stop 4444          
           end if
         else if( bc(side,axis) > 0 )then
           write(*,'("bcOptWave:fill RHS for direct Helmholtz solver, unexpected boundaryCondition=",i4)') bc(side,axis)
           stop 6666
         else if( bc(side,axis) == 0 )then
           !write(*,'("QQQQ Set ghost outside interp to zero: side,axis,grid,numGhost=",4i4)') side,axis,grid,numGhost
           !write(*,'(" bcOptWave: grid,side,axis=",3i3,", ! loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i6)') grid,side,axis,! n1a,n1b,n2a,n2b,n3a,n3b        
            do i3=n3a,n3b
            do i2=n2a,n2b
            do i1=n1a,n1b
            ! -- Set ghost outside interpolation boundaries to zero (for active unused points)
             do ghost=1,numGhost
               j1=i1-is1*ghost
               j2=i2-is2*ghost
               j3=i3-is3*ghost
               u(j1,j2,j3,uc) = 0.
             end do  
            end do
            end do
            end do
         end if
        end do ! end side
        end do ! end axis
     if( assignCornerGhostPoints.eq.1 )then
       ! ---- Assign RHS for corner and edge points ----
       call cornersWave( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange, dimRange, isPeriodic, u, un, mask,rsxy, xy, uTemp, v, boundaryCondition, frequencyArray, pdb, ipar, rpar, ierr )
     end if
     ! ---------------- RETURN ---------------
     return
   end if
   ! if( .true. .and. assignCornerGhostPoints.eq.1 )then  ! ************************** TESTING **************
   !   ! call new routine: nov 28, 2024 -- 
   !   call cornersWave( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,!                     gridIndexRange, dimRange, isPeriodic, u, un, mask,rsxy, xy, uTemp, v, boundaryCondition, !                     frequencyArray, pdb, ipar, rpar, ierr )
   ! end if
   ! ---------------------------------------------------------------
   ! ----------- STAGE I : Assign Dirichlet Conditions -------------
   ! ---------------------------------------------------------------
     ! NOTE: the numGhost args are used in ghost loops
     extraForDirichlet=numGhost
     ff =0. ! default value 
      extra1a=extraForDirichlet
      extra1b=extraForDirichlet
      extra2a=extraForDirichlet
      extra2b=extraForDirichlet
      if( nd.eq.3 )then
        extra3a=extraForDirichlet
        extra3b=extraForDirichlet
      else
        extra3a=0
        extra3b=0
      end if
      if( bc(0,0).lt.0 )then
        extra1a=max(0,extra1a) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions (or internal parallel boundaries)
      else if( bc(0,0).eq.0 )then
        extra1a=numGhost  ! include interpolation points since we assign ghost points outside these
      end if
      ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
      if( bc(1,0).lt.0 )then
        extra1b=max(0,extra1b) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions
      else if( bc(1,0).eq.0 )then
        extra1b=numGhost
      end if
      if( bc(0,1).lt.0 )then
        extra2a=max(0,extra2a) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions (or internal parallel boundaries)
      else if( bc(0,1).eq.0 )then
        extra2a=numGhost  ! include interpolation points since we assign ghost points outside these
      end if
      ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
      if( bc(1,1).lt.0 )then
        extra2b=max(0,extra2b) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions
      else if( bc(1,1).eq.0 )then
        extra2b=numGhost
      end if
      if(  nd.eq.3 )then
       if( bc(0,2).lt.0 )then
         extra3a=max(0,extra3a) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions (or internal parallel boundaries)
       else if( bc(0,2).eq.0 )then
         extra3a=numGhost  ! include interpolation points since we assign ghost points outside these
       end if
       ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
       if( bc(1,2).lt.0 )then
         extra3b=max(0,extra3b) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions
       else if( bc(1,2).eq.0 )then
         extra3b=numGhost
       end if
      end if
      do axis=0,nd-1
      do side=0,1
        if( .true. .or. bc(side,axis).gt.0 )then ! we may set ghost outside interp for implicit
          ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,bc(side,axis)
          n1a=gridIndexRange(0,0)
          n1b=gridIndexRange(1,0)
          n2a=gridIndexRange(0,1)
          n2b=gridIndexRange(1,1)
          n3a=gridIndexRange(0,2)
          n3b=gridIndexRange(1,2)
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
          nn1a=gridIndexRange(0,0)-extra1a
          nn1b=gridIndexRange(1,0)+extra1b
          nn2a=gridIndexRange(0,1)-extra2a
          nn2b=gridIndexRange(1,1)+extra2b
          nn3a=gridIndexRange(0,2)-extra3a
          nn3b=gridIndexRange(1,2)+extra3b
          if( axis.eq.0 )then
            nn1a=gridIndexRange(side,axis)
            nn1b=gridIndexRange(side,axis)
          else if( axis.eq.1 )then
            nn2a=gridIndexRange(side,axis)
            nn2b=gridIndexRange(side,axis)
          else
            nn3a=gridIndexRange(side,axis)
            nn3b=gridIndexRange(side,axis)
          end if
          is=1-2*side
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
          axisp1=mod(axis+1,nd)
          axisp2=mod(axis+2,nd)
          i3=n3a
          if( debug.gt.31 )then
            write(*,'(" bcOptWave: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i5)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
            write(*,'(" bcOptWave: extraForDirichlet,extra1a,extra2a,extra3a=",4i4,", loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i5)') extraForDirichlet,extra1a,extra2a,extra3a,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
          end if
        end if ! if bc>0 
        assignTwilightZone=twilightZone
       if( bc(side,axis).eq.exactBC )then
         ! ==== Set the boundary and ghost with the exact solution ====
         ! NOTE: known solutions are now done in applyBoundaryCondtions.bC 
         if( assignTwilightZone.eq.1 )then
             ! assign extram points in the tangential directions
             extram = numGhost 
               m1a=gridIndexRange(0,0)-extram
               m1b=gridIndexRange(1,0)+extram
               m2a=gridIndexRange(0,1)-extram
               m2b=gridIndexRange(1,1)+extram
               if( nd.eq.2 )then
                 m3a=gridIndexRange(0,2)
                 m3b=gridIndexRange(1,2)
               else
                 m3a=gridIndexRange(0,2)-extram
                 m3b=gridIndexRange(1,2)+extram
               end if
               if( axis.eq.0 )then
                m1a=gridIndexRange(side,axis)
                m1b=gridIndexRange(side,axis)
               else if( axis.eq.1 )then
                m2a=gridIndexRange(side,axis)
                m2b=gridIndexRange(side,axis)
               else
                m3a=gridIndexRange(side,axis)
                m3b=gridIndexRange(side,axis)
               end if
             ff=0.
             do i3=m3a,m3b
             do i2=m2a,m2b
             do i1=m1a,m1b
               if( mask(i1,i2,i3).ne.0 )then
                 do ghost=0,numGhost
                   j1=i1-is1*ghost
                   j2=i2-is2*ghost
                   j3=i3-is3*ghost
                     if( assignTwilightZone.eq.1 )then
                       ! compute RHS from TZ
                       if( nd.eq.2 )then
                         call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue )
                       else
                         call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),xy(j1,j2,j3,2),t,uc,ue )
                       end if
                       ff = ue
                     else if( assignKnownSolutionAtBoundaries.eq.1 )then
                       ! -- we set inhomogeneous Dirichlet values for some known solutions 
                       if( knownSolutionOption.eq.planeWave )then
                         ! --- evaluate the plane wave solution ---
                         if( nd.eq.2 )then
                           ff = ampPlaneWave*sin( kxPlaneWave*xy(j1,j2,j3,0) + kyPlaneWave*xy(j1,j2,j3,1) - omegaPlaneWave*t )
                         else
                           ff = ampPlaneWave*sin( kxPlaneWave*xy(j1,j2,j3,0) + kyPlaneWave*xy(j1,j2,j3,1) + kzPlaneWave*xy(j1,j2,j3,2) - omegaPlaneWave*t )
                         end if 
                       else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                         ! Eval the Gaussian plane wave solution
                         !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                         !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                         !  
                         if( nd.eq.2 )then
                           xi = kxGPW*(xy(j1,j2,j3,0)-x0GPW) + kyGPW*(xy(j1,j2,j3,1)-y0GPW) - c*t
                         else
                           xi = kxGPW*(xy(j1,j2,j3,0)-x0GPW) + kyGPW*(xy(j1,j2,j3,1)-y0GPW) + kzGPW*(xy(j1,j2,j3,2)-z0GPW) - c*t
                         end if 
                         ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                       else if( knownSolutionOption.eq.boxHelmholtz ) then
                         ! --- evaluate the boxHelmholtz solution ---
                         ! For multi-freq we add up all the component frequencies
                         ff = 0. 
                         do freq=0,numberOfFrequencies-1
                           ! kx = kxBoxHelmholtz + twoPi*freq
                           ! ky = kyBoxHelmholtz + twoPi*freq
                           ! kz = kzBoxHelmholtz + twoPi*freq
                           ! coswt = cos( frequencyArray(freq)*t )
                               ! This macro is used in bcOptWave.bf90
                               omega = frequencyArray(freq);
                               kx = kxBoxHelmholtz*(freq*.5+1.)
                               ky = kyBoxHelmholtz*(freq*.5+1.)
                               kz = kzBoxHelmholtz*(freq*.5+1.)
                               ! kx = kxBoxHelmholtz + twoPi*freq.
                               ! ky = kyBoxHelmholtz + twoPi*freq.
                               ! kz = kzBoxHelmholtz + twoPi*freq.    
                           coswt = cos( omega*t )
                           if( nd.eq.2 )then
                             ff = ff + sin( kx*xy(j1,j2,j3,0) ) * sin( ky*xy(j1,j2,j3,1) ) *coswt
                           else
                             ff = ff + sin( kx*xy(j1,j2,j3,0) ) * sin( ky*xy(j1,j2,j3,1) ) * sin( kz*xy(j1,j2,j3,2) ) *coswt
                           end if 
                         end do
                       else if( knownSolutionOption.eq.polyPeriodic ) then
                         ! --- evaluate the polyPeriodic solution ---
                         if( nd.eq.2 )then
                           ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(j1,j2,j3,0) + b1PolyPeriodic*xy(j1,j2,j3,1)                                 ) *coswt
                         else
                           ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(j1,j2,j3,0) + b1PolyPeriodic*xy(j1,j2,j3,1) + c1PolyPeriodic*xy(j1,j2,j3,2) ) *coswt
                         end if 
                       else
                         stop 9876
                       end if 
                     end if
                   u(j1,j2,j3,uc) = ff
                 end do
               end if ! mask .ne. 0
              end do
              end do
              end do
         end if
       else if( bc(side,axis).eq.dirichlet )then
         if( bcApproach.eq.useCompatibilityBoundaryConditions )then
           ! --- Assign values on the boundary for CBCs ---
           if( orderOfAccuracy.eq.2 )then
             if( gridType.eq.rectangular )then
                 ! if( t.le.2*dt )then
                 !   write(*,'("assign extended Dirichlet BC: side,axis=",2i2," nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                 ! end if
                 ff=0.
                 if( addForcingBC.eq.1 )then
                   ! beginLoops3d()
                   ! *wdh* Dec 1, 2023: assign extended bndry
                   do i3=nn3a,nn3b
                   do i2=nn2a,nn2b
                   do i1=nn1a,nn1b
                     if( mask(i1,i2,i3).ne.0 )then
                       ! --- get the RHS to the Dirichlet BC ---
                       ! #If #FORCING eq "USEFORCING" *wdh* This was a bug, turned off Nov 22, 2023
                         if( assignTwilightZone.eq.1 )then
                           ! compute RHS from TZ
                           if( nd.eq.2 )then
                             call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                           else
                             call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                           end if
                           ff = ue
                         else if( assignKnownSolutionAtBoundaries.eq.1 )then
                           ! -- we set inhomogeneous Dirichlet values for some known solutions 
                           if( knownSolutionOption.eq.planeWave )then
                             ! --- evaluate the plane wave solution ---
                             if( nd.eq.2 )then
                               ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                             else
                               ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                             end if 
                           else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                             ! Eval the Gaussian plane wave solution
                             !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                             !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                             !  
                             if( nd.eq.2 )then
                               xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) - c*t
                             else
                               xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) + kzGPW*(xy(i1,i2,i3,2)-z0GPW) - c*t
                             end if 
                             ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                           else if( knownSolutionOption.eq.boxHelmholtz ) then
                             ! --- evaluate the boxHelmholtz solution ---
                             ! For multi-freq we add up all the component frequencies
                             ff = 0. 
                             do freq=0,numberOfFrequencies-1
                               ! kx = kxBoxHelmholtz + twoPi*freq
                               ! ky = kyBoxHelmholtz + twoPi*freq
                               ! kz = kzBoxHelmholtz + twoPi*freq
                               ! coswt = cos( frequencyArray(freq)*t )
                                   ! This macro is used in bcOptWave.bf90
                                   omega = frequencyArray(freq);
                                   kx = kxBoxHelmholtz*(freq*.5+1.)
                                   ky = kyBoxHelmholtz*(freq*.5+1.)
                                   kz = kzBoxHelmholtz*(freq*.5+1.)
                                   ! kx = kxBoxHelmholtz + twoPi*freq.
                                   ! ky = kyBoxHelmholtz + twoPi*freq.
                                   ! kz = kzBoxHelmholtz + twoPi*freq.    
                               coswt = cos( omega*t )
                               if( nd.eq.2 )then
                                 ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) *coswt
                               else
                                 ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) * sin( kz*xy(i1,i2,i3,2) ) *coswt
                               end if 
                             end do
                           else if( knownSolutionOption.eq.polyPeriodic ) then
                             ! --- evaluate the polyPeriodic solution ---
                             if( nd.eq.2 )then
                               ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1)                                 ) *coswt
                             else
                               ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                             end if 
                           else
                             stop 9876
                           end if 
                         end if
                       ! #End
                       ! --- Dirichlet BC ---
                       u(i1,i2,i3,uc)=ff
                       ! write(*,'("assignDirichletBndry: side,axis=",2i2," i1,i2=",2i3," ff=",e10.2)') side,axis,i1,i2,ff
                     end if ! mask .ne. 0
                    end do
                    end do
                    end do
                 else
                   ! Homogenous BCs 
                   ! beginLoops3d()
                   do i3=nn3a,nn3b
                   do i2=nn2a,nn2b
                   do i1=nn1a,nn1b
                     u(i1,i2,i3,uc)=0.
                    end do
                    end do
                    end do
                 end if
             else
                 ! if( t.le.2*dt )then
                 !   write(*,'("assign extended Dirichlet BC: side,axis=",2i2," nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                 ! end if
                 ff=0.
                 if( addForcingBC.eq.1 )then
                   ! beginLoops3d()
                   ! *wdh* Dec 1, 2023: assign extended bndry
                   do i3=nn3a,nn3b
                   do i2=nn2a,nn2b
                   do i1=nn1a,nn1b
                     if( mask(i1,i2,i3).ne.0 )then
                       ! --- get the RHS to the Dirichlet BC ---
                       ! #If #FORCING eq "USEFORCING" *wdh* This was a bug, turned off Nov 22, 2023
                         if( assignTwilightZone.eq.1 )then
                           ! compute RHS from TZ
                           if( nd.eq.2 )then
                             call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                           else
                             call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                           end if
                           ff = ue
                         else if( assignKnownSolutionAtBoundaries.eq.1 )then
                           ! -- we set inhomogeneous Dirichlet values for some known solutions 
                           if( knownSolutionOption.eq.planeWave )then
                             ! --- evaluate the plane wave solution ---
                             if( nd.eq.2 )then
                               ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                             else
                               ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                             end if 
                           else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                             ! Eval the Gaussian plane wave solution
                             !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                             !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                             !  
                             if( nd.eq.2 )then
                               xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) - c*t
                             else
                               xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) + kzGPW*(xy(i1,i2,i3,2)-z0GPW) - c*t
                             end if 
                             ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                           else if( knownSolutionOption.eq.boxHelmholtz ) then
                             ! --- evaluate the boxHelmholtz solution ---
                             ! For multi-freq we add up all the component frequencies
                             ff = 0. 
                             do freq=0,numberOfFrequencies-1
                               ! kx = kxBoxHelmholtz + twoPi*freq
                               ! ky = kyBoxHelmholtz + twoPi*freq
                               ! kz = kzBoxHelmholtz + twoPi*freq
                               ! coswt = cos( frequencyArray(freq)*t )
                                   ! This macro is used in bcOptWave.bf90
                                   omega = frequencyArray(freq);
                                   kx = kxBoxHelmholtz*(freq*.5+1.)
                                   ky = kyBoxHelmholtz*(freq*.5+1.)
                                   kz = kzBoxHelmholtz*(freq*.5+1.)
                                   ! kx = kxBoxHelmholtz + twoPi*freq.
                                   ! ky = kyBoxHelmholtz + twoPi*freq.
                                   ! kz = kzBoxHelmholtz + twoPi*freq.    
                               coswt = cos( omega*t )
                               if( nd.eq.2 )then
                                 ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) *coswt
                               else
                                 ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) * sin( kz*xy(i1,i2,i3,2) ) *coswt
                               end if 
                             end do
                           else if( knownSolutionOption.eq.polyPeriodic ) then
                             ! --- evaluate the polyPeriodic solution ---
                             if( nd.eq.2 )then
                               ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1)                                 ) *coswt
                             else
                               ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                             end if 
                           else
                             stop 9876
                           end if 
                         end if
                       ! #End
                       ! --- Dirichlet BC ---
                       u(i1,i2,i3,uc)=ff
                       ! write(*,'("assignDirichletBndry: side,axis=",2i2," i1,i2=",2i3," ff=",e10.2)') side,axis,i1,i2,ff
                     end if ! mask .ne. 0
                    end do
                    end do
                    end do
                 else
                   ! Homogenous BCs 
                   ! beginLoops3d()
                   do i3=nn3a,nn3b
                   do i2=nn2a,nn2b
                   do i1=nn1a,nn1b
                     u(i1,i2,i3,uc)=0.
                    end do
                    end do
                    end do
                 end if
             end if 
           else if( orderOfAccuracy.eq.4 )then
             if( t.le.3*dt .and. debug.gt.1 )then
               write(*,'("APPLY BC order =4 useCompatibilityBoundaryConditions, side,axis=",2i4)') side,axis
             end if
               stop 444 
           else if( orderOfAccuracy.eq.6 )then 
             write(*,'("CgWave::bcOpt:ERROR: CBC order=6 turned off, use LCBC")') 
             stop 666
           else if( orderOfAccuracy.eq.8 )then   
             stop 888
           else
             write(*,'("CgWave::bcOpt:ERROR: unexpected orderOfAccuracy=",i6)') orderOfAccuracy
             stop 8888
           end if
         else 
           ! ----- assign boundary and then ghost by extrapolation ----
           ! This was a temporary method until other approaches were implemented
           if( orderOfAccuracy.eq.2 )then
                 ff=0.
                 if( addForcingBC.eq.1 )then  
                    do i3=n3a,n3b
                    do i2=n2a,n2b
                    do i1=n1a,n1b
                     if( mask(i1,i2,i3).ne.0 )then
                       ! --- get the RHS to the Dirichlet BC ---
                         if( assignTwilightZone.eq.1 )then
                           ! compute RHS from TZ
                           if( nd.eq.2 )then
                             call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                           else
                             call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                           end if
                           ff = ue
                         else if( assignKnownSolutionAtBoundaries.eq.1 )then
                           ! -- we set inhomogeneous Dirichlet values for some known solutions 
                           if( knownSolutionOption.eq.planeWave )then
                             ! --- evaluate the plane wave solution ---
                             if( nd.eq.2 )then
                               ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                             else
                               ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                             end if 
                           else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                             ! Eval the Gaussian plane wave solution
                             !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                             !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                             !  
                             if( nd.eq.2 )then
                               xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) - c*t
                             else
                               xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) + kzGPW*(xy(i1,i2,i3,2)-z0GPW) - c*t
                             end if 
                             ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                           else if( knownSolutionOption.eq.boxHelmholtz ) then
                             ! --- evaluate the boxHelmholtz solution ---
                             ! For multi-freq we add up all the component frequencies
                             ff = 0. 
                             do freq=0,numberOfFrequencies-1
                               ! kx = kxBoxHelmholtz + twoPi*freq
                               ! ky = kyBoxHelmholtz + twoPi*freq
                               ! kz = kzBoxHelmholtz + twoPi*freq
                               ! coswt = cos( frequencyArray(freq)*t )
                                   ! This macro is used in bcOptWave.bf90
                                   omega = frequencyArray(freq);
                                   kx = kxBoxHelmholtz*(freq*.5+1.)
                                   ky = kyBoxHelmholtz*(freq*.5+1.)
                                   kz = kzBoxHelmholtz*(freq*.5+1.)
                                   ! kx = kxBoxHelmholtz + twoPi*freq.
                                   ! ky = kyBoxHelmholtz + twoPi*freq.
                                   ! kz = kzBoxHelmholtz + twoPi*freq.    
                               coswt = cos( omega*t )
                               if( nd.eq.2 )then
                                 ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) *coswt
                               else
                                 ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) * sin( kz*xy(i1,i2,i3,2) ) *coswt
                               end if 
                             end do
                           else if( knownSolutionOption.eq.polyPeriodic ) then
                             ! --- evaluate the polyPeriodic solution ---
                             if( nd.eq.2 )then
                               ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1)                                 ) *coswt
                             else
                               ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                             end if 
                           else
                             stop 9876
                           end if 
                         end if
                       ! --- Dirichlet BC ---
                       u(i1,i2,i3,uc)=ff
                       ! -- extrapolate ghost ---
                       do ghost=1,numGhost
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,uc) = (3.*u(j1+is1,j2+is2,j3+is3,uc)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)) 
                       end do
                     end if ! mask .ne. 0
                    end do
                    end do
                    end do
                 else
                   ! --- no forcing ----
                   if( numGhost.eq.1 )then
                     ghost=1
                      do i3=n3a,n3b
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                       ! --- Dirichlet BC ---
                       u(i1,i2,i3,uc)=0.
                       ! -- extrapolate ghost ---
                       j1=i1-is1*ghost
                       j2=i2-is2*ghost
                       j3=i3-is3*ghost         
                       u(j1,j2,j3,uc) = (3.*u(j1+is1,j2+is2,j3+is3,uc)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)) 
                      end do
                      end do
                      end do
                   else 
                      do i3=n3a,n3b
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                       ! --- Dirichlet BC ---
                       u(i1,i2,i3,uc)=0.
                       ! -- extrapolate ghost ---
                       do ghost=1,numGhost
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,uc) = (3.*u(j1+is1,j2+is2,j3+is3,uc)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)) 
                       end do
                      end do
                      end do
                      end do
                   end if
                   if( .false. )then
                     ! -- two loops ---
                     ! --- no forcing ----
                      do i3=n3a,n3b
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                       ! --- Dirichlet BC ---
                       u(i1,i2,i3,uc)=0.
                      end do
                      end do
                      end do
                     ! -- extrapolate ghost ---
                     do ghost=1,numGhost
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost               
                         u(j1,j2,j3,uc) = (3.*u(j1+is1,j2+is2,j3+is3,uc)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)) 
                        end do
                        end do
                        end do
                     end do
                   end if
                 endif
           else if( orderOfAccuracy.eq.4 )then
           else if( orderOfAccuracy.eq.6 )then   
           else if( orderOfAccuracy.eq.8 )then   
             stop 888
           else
             write(*,'("CgWave::bcOpt:ERROR: unexpected orderOfAccuracy=",i6)') orderOfAccuracy
             stop 8888
           end if
         end if
       end if ! end if dirichlet 
      end do ! end side
      end do ! end axis
   if( .false. ) then
    n1a=gridIndexRange(0,0)
    n1b=gridIndexRange(1,0)
    n2a=gridIndexRange(0,1)
    n2b=gridIndexRange(1,1)
    n3a=gridIndexRange(0,2)
    n3b=gridIndexRange(1,2)    
     write(*,'(/,"bcOpt: After dirichlet BC, n1a,n1b,n2a,n2b,n3a,n3b=",6i5, " numGhost3=",i3)') n1a,n1b,n2a,n2b,n3a,n3b,numGhost3
     do i3=n3a-numGhost3,n3b+numGhost3
       if( nd.eq.3 )then
         write(*,'("i3=",i4)') i3
       end if
       do i2=n2a-numGhost,n2b+numGhost
         write(*,'("i2=",i4,1x,100(1pe10.2))') i2,(u(i1,i2,i3,0),i1=n1a-numGhost,n1b+numGhost)
       end do
     end do
   end if
    ! if( .true. )then
   !   return   ! ************************ TESTING TEMP **********************
   ! end if
   ! --  Extrap values on remaining sides to give initial values 
   !     --> maybe we only need to do this along extended boundaries on
   !         curvilinear grids so we have values for the Neumann BC
    extra1a=numGhost
    extra1b=numGhost
    extra2a=numGhost
    extra2b=numGhost
    if( nd.eq.3 )then
      extra3a=numGhost
      extra3b=numGhost
    else
      extra3a=0
      extra3b=0
    end if
    if( bc(0,0).lt.0 )then
      extra1a=max(0,extra1a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
    else if( bc(0,0).eq.0 )then
      extra1a=numGhost  ! include interpolation points since we assign ghost points outside these
    end if
    ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
    if( bc(1,0).lt.0 )then
      extra1b=max(0,extra1b) ! over-ride numGhost=-1 : assign ends in periodic directions
    else if( bc(1,0).eq.0 )then
      extra1b=numGhost
    end if
    if( bc(0,1).lt.0 )then
      extra2a=max(0,extra2a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
    else if( bc(0,1).eq.0 )then
      extra2a=numGhost  ! include interpolation points since we assign ghost points outside these
    end if
    ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
    if( bc(1,1).lt.0 )then
      extra2b=max(0,extra2b) ! over-ride numGhost=-1 : assign ends in periodic directions
    else if( bc(1,1).eq.0 )then
      extra2b=numGhost
    end if
    if(  nd.eq.3 )then
     if( bc(0,2).lt.0 )then
       extra3a=max(0,extra3a) ! over-ride numGhost=-1 : assign ends in periodic directions (or internal parallel boundaries)
     else if( bc(0,2).eq.0 )then
       extra3a=numGhost  ! include interpolation points since we assign ghost points outside these
     end if
     ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
     if( bc(1,2).lt.0 )then
       extra3b=max(0,extra3b) ! over-ride numGhost=-1 : assign ends in periodic directions
     else if( bc(1,2).eq.0 )then
       extra3b=numGhost
     end if
    end if
    do axis=0,nd-1
    do side=0,1
      if( .true. .or. bc(side,axis).gt.0 )then ! we may set ghost outside interp for implicit
        ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,bc(side,axis)
        n1a=gridIndexRange(0,0)
        n1b=gridIndexRange(1,0)
        n2a=gridIndexRange(0,1)
        n2b=gridIndexRange(1,1)
        n3a=gridIndexRange(0,2)
        n3b=gridIndexRange(1,2)
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
        nn1a=gridIndexRange(0,0)-extra1a
        nn1b=gridIndexRange(1,0)+extra1b
        nn2a=gridIndexRange(0,1)-extra2a
        nn2b=gridIndexRange(1,1)+extra2b
        nn3a=gridIndexRange(0,2)-extra3a
        nn3b=gridIndexRange(1,2)+extra3b
        if( axis.eq.0 )then
          nn1a=gridIndexRange(side,axis)
          nn1b=gridIndexRange(side,axis)
        else if( axis.eq.1 )then
          nn2a=gridIndexRange(side,axis)
          nn2b=gridIndexRange(side,axis)
        else
          nn3a=gridIndexRange(side,axis)
          nn3b=gridIndexRange(side,axis)
        end if
        is=1-2*side
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
        axisp1=mod(axis+1,nd)
        axisp2=mod(axis+2,nd)
        i3=n3a
        if( debug.gt.31 )then
          write(*,'(" bcOptWave: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i5)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
          write(*,'(" bcOptWave: numGhost,extra1a,extra2a,extra3a=",4i4,", loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i5)') numGhost,extra1a,extra2a,extra3a,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
        end if
      end if ! if bc>0 
      assignTwilightZone=twilightZone
     ! if( ( (bc(side,axis).ne.dirichlet .and. bc(side,axis).ne.exactBC .and. bc(side,axis).ne.absorbing .and. bc(side,axis).ne.abcEM2 ) !      .or. bcApproach.eq.useCompatibilityBoundaryConditions ) !      .and. bc(side,axis).gt.0 )then
     if(  bc(side,axis).gt.0 .and. bc(side,axis).ne.exactBC )then ! *wdh* Dec 2, 2023 -- always do this at least for CBCs
       ! If the grid is too coarse then we can only extrapolate using point from n1a to n1b 
       !         E--+--+--+--+
       !           n1a       n1b
       maxExtrapWidth = gridIndexRange(1,axis)-gridIndexRange(0,axis)+1
       ! extrapWidth = min(orderOfAccuracy+1,maxExtrapWidth) ! *wdh* July 21, 2024
       extrapWidth = min(orderOfExtrapolation,maxExtrapWidth)
       if( extrapWidth .lt. orderOfExtrapolation )then
         write(*,'("bcOptWave:WARNING: reducing extrapolation width to ",i2," since there are not enough grid points")') extrapWidth
       end if
       if( extrapWidth==2 )then
           ! OLD WAY: extrap all ghost, this may overwite ghost points on extended dirichlet boundaries
           ! ghost-loops will assign extra points in tangential directions
           ! if( t.le. 2*dt )then
           !    write(*,'("extrapolateGhost: side,axis=",2i3," nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i5)') side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
           ! end if 
           ! beginGhostLoops3d()  
           ! New way: Exclude the ends when adjacent boundaries are dirichlet or exact: *wdh* Dec 2, 2023
           extram=numGhost
             ! Assume sides with axis=0 are adjacent sides (fixed later)
             if( bc(0,0)==dirichlet .or. bc(0,0)==exactBC )then
               m1a= gridIndexRange(0,0)+1
             else
               m1a=gridIndexRange(0,0)-extram
             end if
             if( bc(1,0)==dirichlet .or. bc(1,0)==exactBC )then
               m1b = gridIndexRange(1,0)-1
             else
               m1b = gridIndexRange(1,0)+extram
             end if
             ! Assume sides with axis=1 are adjacent sides (fixed later)
             if( bc(0,1)==dirichlet .or. bc(0,1)==exactBC )then
               m2a= gridIndexRange(0,1)+1
             else
               m2a=gridIndexRange(0,1)-extram
             end if
             if( bc(1,1)==dirichlet .or. bc(1,1)==exactBC )then
               m2b = gridIndexRange(1,1)-1
             else
               m2b = gridIndexRange(1,1)+extram
             end if  
             if( nd.eq.2 )then
               m3a=gridIndexRange(0,2)
               m3b=gridIndexRange(1,2)
             else
               if( bc(0,2)==dirichlet .or. bc(0,2)==exactBC )then
                 m3a= gridIndexRange(0,2)+1
               else
                 m3a=gridIndexRange(0,2)-extram
               end if
               if( bc(1,2)==dirichlet .or. bc(1,2)==exactBC )then
                 m3b = gridIndexRange(1,2)-1
               else
                 m3b = gridIndexRange(1,2)+extram
               end if
             end if
             if( axis.eq.0 )then
               m1a=gridIndexRange(side,axis)
               m1b=gridIndexRange(side,axis)
             else if( axis.eq.1 )then
               m2a=gridIndexRange(side,axis)
               m2b=gridIndexRange(side,axis)
             else
               m3a=gridIndexRange(side,axis)
               m3b=gridIndexRange(side,axis)
             end if
           if( .false. .and. t.le. 2*dt )then
             write(*,'("extrapolateGhost: grid,side,axis=",3i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i5)') grid,side,axis,m1a,m1b,m2a,m2b,m3a,m3b
           end if  
           do i3=m3a,m3b
           do i2=m2a,m2b
           do i1=m1a,m1b
             if( mask(i1,i2,i3).ne.0 )then
               ! -- extrapolate ghost ---
               do ghost=1,numGhost
                 j1=i1-is1*ghost
                 j2=i2-is2*ghost
                 j3=i3-is3*ghost
                 u(j1,j2,j3,uc) = (2.*u(j1+is1,j2+is2,j3+is3,uc)-u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)) 
                 ! write(*,'("extrap: u(",3i4,")=",1pe10.2)') j1,j2,j3,u(j1,j2,j3,uc)
               end do
             end if ! mask .ne. 0
            end do
            end do
            end do
       else if( extrapWidth==3 )then
           ! OLD WAY: extrap all ghost, this may overwite ghost points on extended dirichlet boundaries
           ! ghost-loops will assign extra points in tangential directions
           ! if( t.le. 2*dt )then
           !    write(*,'("extrapolateGhost: side,axis=",2i3," nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i5)') side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
           ! end if 
           ! beginGhostLoops3d()  
           ! New way: Exclude the ends when adjacent boundaries are dirichlet or exact: *wdh* Dec 2, 2023
           extram=numGhost
             ! Assume sides with axis=0 are adjacent sides (fixed later)
             if( bc(0,0)==dirichlet .or. bc(0,0)==exactBC )then
               m1a= gridIndexRange(0,0)+1
             else
               m1a=gridIndexRange(0,0)-extram
             end if
             if( bc(1,0)==dirichlet .or. bc(1,0)==exactBC )then
               m1b = gridIndexRange(1,0)-1
             else
               m1b = gridIndexRange(1,0)+extram
             end if
             ! Assume sides with axis=1 are adjacent sides (fixed later)
             if( bc(0,1)==dirichlet .or. bc(0,1)==exactBC )then
               m2a= gridIndexRange(0,1)+1
             else
               m2a=gridIndexRange(0,1)-extram
             end if
             if( bc(1,1)==dirichlet .or. bc(1,1)==exactBC )then
               m2b = gridIndexRange(1,1)-1
             else
               m2b = gridIndexRange(1,1)+extram
             end if  
             if( nd.eq.2 )then
               m3a=gridIndexRange(0,2)
               m3b=gridIndexRange(1,2)
             else
               if( bc(0,2)==dirichlet .or. bc(0,2)==exactBC )then
                 m3a= gridIndexRange(0,2)+1
               else
                 m3a=gridIndexRange(0,2)-extram
               end if
               if( bc(1,2)==dirichlet .or. bc(1,2)==exactBC )then
                 m3b = gridIndexRange(1,2)-1
               else
                 m3b = gridIndexRange(1,2)+extram
               end if
             end if
             if( axis.eq.0 )then
               m1a=gridIndexRange(side,axis)
               m1b=gridIndexRange(side,axis)
             else if( axis.eq.1 )then
               m2a=gridIndexRange(side,axis)
               m2b=gridIndexRange(side,axis)
             else
               m3a=gridIndexRange(side,axis)
               m3b=gridIndexRange(side,axis)
             end if
           if( .false. .and. t.le. 2*dt )then
             write(*,'("extrapolateGhost: grid,side,axis=",3i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i5)') grid,side,axis,m1a,m1b,m2a,m2b,m3a,m3b
           end if  
           do i3=m3a,m3b
           do i2=m2a,m2b
           do i1=m1a,m1b
             if( mask(i1,i2,i3).ne.0 )then
               ! -- extrapolate ghost ---
               do ghost=1,numGhost
                 j1=i1-is1*ghost
                 j2=i2-is2*ghost
                 j3=i3-is3*ghost
                 u(j1,j2,j3,uc) = (3.*u(j1+is1,j2+is2,j3+is3,uc)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)) 
                 ! write(*,'("extrap: u(",3i4,")=",1pe10.2)') j1,j2,j3,u(j1,j2,j3,uc)
               end do
             end if ! mask .ne. 0
            end do
            end do
            end do
       else
         write(*,'("CgWave::bcOpt:ERROR: unexpected extrapWidth=",i3," for orderOfAccuracy=",i6)') extrapWidth,orderOfAccuracy
         stop 8888
       end if
     end if
    end do ! end side
    end do ! end axis
   if( .false.  ) then
     n1a=gridIndexRange(0,0)
     n1b=gridIndexRange(1,0)
     n2a=gridIndexRange(0,1)
     n2b=gridIndexRange(1,1)
     n3a=gridIndexRange(0,2)
     n3b=gridIndexRange(1,2)
     write(*,'(/,"bcOpt: After dirichlet BC and extrapolate ghost, grid=",i3," n1a,n1b,n2a,n2b,n3a,n3b=",6i5)') grid,n1a,n1b,n2a,n2b,n3a,n3b
     do i3=n3a-numGhost3,n3b+numGhost3
       if( nd.eq.3 )then
         write(*,'("i3=",i4)') i3
       else
         write(*,'("i1=[",i4," : ",i4,"]")') n1a-numGhost,n1b+numGhost
       end if
       do i2=n2a-numGhost,n2b+numGhost
         write(*,'("i2=",i4,1x,100(1pe10.2))') i2,(u(i1,i2,i3,0),i1=n1a-numGhost,n1b+numGhost)
       end do
     end do
   end if
   ! if( .true. )then
   !   write(*,'("bcOpt:TEST: return after Stage Ia extrapolate ghost")')
   !   return
   ! end if
   ! TESTING: Is this needed now ?? --> Needed for Dirichlet-Neumann Corners order=4
   if( assignCornerGhostPoints.eq.1 .and. bcApproach.eq.useCompatibilityBoundaryConditions )then
     ! Fill in corners and edges -- this assumes the corners and edges do not need other ghost points!!
     call cornersWave( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange, dimRange, isPeriodic, u, un, mask,rsxy, xy, uTemp, v, boundaryCondition, frequencyArray, pdb, ipar, rpar, ierr )
   end if  
   ! ---------------------------------------------------------------------
   ! ----------- STAGE II : Neumann-like Boundary Conditions -------------
   ! -----------            Ghost values for CBCs            -------------
   ! ---------------------------------------------------------------------
   ! CHECK ME --> numGhost here ??
   extraForNeumann=0 ! only assign Neumann conditions to the adjacent boundaries
    extra1a=extraForNeumann
    extra1b=extraForNeumann
    extra2a=extraForNeumann
    extra2b=extraForNeumann
    if( nd.eq.3 )then
      extra3a=extraForNeumann
      extra3b=extraForNeumann
    else
      extra3a=0
      extra3b=0
    end if
    if( bc(0,0).lt.0 )then
      extra1a=max(0,extra1a) ! over-ride extraForNeumann=-1 : assign ends in periodic directions (or internal parallel boundaries)
    else if( bc(0,0).eq.0 )then
      extra1a=numGhost  ! include interpolation points since we assign ghost points outside these
    end if
    ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
    if( bc(1,0).lt.0 )then
      extra1b=max(0,extra1b) ! over-ride extraForNeumann=-1 : assign ends in periodic directions
    else if( bc(1,0).eq.0 )then
      extra1b=numGhost
    end if
    if( bc(0,1).lt.0 )then
      extra2a=max(0,extra2a) ! over-ride extraForNeumann=-1 : assign ends in periodic directions (or internal parallel boundaries)
    else if( bc(0,1).eq.0 )then
      extra2a=numGhost  ! include interpolation points since we assign ghost points outside these
    end if
    ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
    if( bc(1,1).lt.0 )then
      extra2b=max(0,extra2b) ! over-ride extraForNeumann=-1 : assign ends in periodic directions
    else if( bc(1,1).eq.0 )then
      extra2b=numGhost
    end if
    if(  nd.eq.3 )then
     if( bc(0,2).lt.0 )then
       extra3a=max(0,extra3a) ! over-ride extraForNeumann=-1 : assign ends in periodic directions (or internal parallel boundaries)
     else if( bc(0,2).eq.0 )then
       extra3a=numGhost  ! include interpolation points since we assign ghost points outside these
     end if
     ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
     if( bc(1,2).lt.0 )then
       extra3b=max(0,extra3b) ! over-ride extraForNeumann=-1 : assign ends in periodic directions
     else if( bc(1,2).eq.0 )then
       extra3b=numGhost
     end if
    end if
    do axis=0,nd-1
    do side=0,1
      if( .true. .or. bc(side,axis).gt.0 )then ! we may set ghost outside interp for implicit
        ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,bc(side,axis)
        n1a=gridIndexRange(0,0)
        n1b=gridIndexRange(1,0)
        n2a=gridIndexRange(0,1)
        n2b=gridIndexRange(1,1)
        n3a=gridIndexRange(0,2)
        n3b=gridIndexRange(1,2)
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
        nn1a=gridIndexRange(0,0)-extra1a
        nn1b=gridIndexRange(1,0)+extra1b
        nn2a=gridIndexRange(0,1)-extra2a
        nn2b=gridIndexRange(1,1)+extra2b
        nn3a=gridIndexRange(0,2)-extra3a
        nn3b=gridIndexRange(1,2)+extra3b
        if( axis.eq.0 )then
          nn1a=gridIndexRange(side,axis)
          nn1b=gridIndexRange(side,axis)
        else if( axis.eq.1 )then
          nn2a=gridIndexRange(side,axis)
          nn2b=gridIndexRange(side,axis)
        else
          nn3a=gridIndexRange(side,axis)
          nn3b=gridIndexRange(side,axis)
        end if
        is=1-2*side
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
        axisp1=mod(axis+1,nd)
        axisp2=mod(axis+2,nd)
        i3=n3a
        if( debug.gt.31 )then
          write(*,'(" bcOptWave: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i5)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
          write(*,'(" bcOptWave: extraForNeumann,extra1a,extra2a,extra3a=",4i4,", loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i5)') extraForNeumann,extra1a,extra2a,extra3a,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
        end if
      end if ! if bc>0 
      assignTwilightZone=twilightZone
     if( bc(side,axis).eq.dirichlet .and. bcApproach.eq.useCompatibilityBoundaryConditions )then
       ! -- fill ghost using CBCs ----
       if( orderOfAccuracy.eq.2 )then
             if( forcingOption.eq.noForcing )then
               if( gridType.eq.rectangular )then
                 if( nd.eq.2 )then
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     ! Try this: *wdh* Dec 2, 2023 : Exclude the ends when adjacent boundaries are dirichlet or exact:
                     ! getRestrictedLoopBounds(side,axis,extram, m1a,m1b,m2a,m2b,m3a,m3b)
                     ! if( t.le. 2*dt )then
                     !    write(*,'("dirichlet CBC: Stage I: side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i3)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
                     ! end if
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ghost =1 
                         j1=i1-is1*ghost; j2=i2-is2*ghost; j3=i3-is3*ghost
                         ! --- get the compatibility forcings for order=2 ---
                             ! No forcing, do nothing 
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22r or ulaplacian23r 
                           uLap = ulaplacian22r(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2            ! residual in equation using current ghost value
                           a11 = 1./( dx(axis)**2 )                                  ! coeff of u(-1) in r1 
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                           ! if( j1==7 .and. j2==-1 )then
                           !   write(*,'(" BC2: j1,j2=",2i3," set u=",e10.2)') j1,j2,u(j1,j2,j3,0)
                           ! end if 
                         ! if( .false. )then
                         !   ! check the error
                         !   call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                         !   write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         ! end if
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ! This should have been done already
                         ! ghost = 1 
                         ! j1=i1-is1*ghost
                         ! j2=i2-is2*ghost
                         ! j3=i3-is3*ghost
                         ! u(j1,j2,j3,uc)=extrap3(u,i1,i2,i3,uc,is1,is2,is3)
                         ! #If "2" eq "2" 
                         !   if( numGhost.gt.1 )then
                         !     ! extrap second ghost (UPW)
                         !     ghost = 2
                         !     k1=i1-is1*ghost
                         !     k2=i2-is2*ghost
                         !     k3=i3-is3*ghost              
                         !     u(k1,k2,k3,uc)=extrap3(u,j1,j2,j3,uc,is1,is2,is3)
                         !   end if
                         ! #End
                       end if ! mask 
                     end do
                     end do
                     end do
                 else
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     ! Try this: *wdh* Dec 2, 2023 : Exclude the ends when adjacent boundaries are dirichlet or exact:
                     ! getRestrictedLoopBounds(side,axis,extram, m1a,m1b,m2a,m2b,m3a,m3b)
                     ! if( t.le. 2*dt )then
                     !    write(*,'("dirichlet CBC: Stage I: side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i3)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
                     ! end if
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ghost =1 
                         j1=i1-is1*ghost; j2=i2-is2*ghost; j3=i3-is3*ghost
                         ! --- get the compatibility forcings for order=2 ---
                             ! No forcing, do nothing 
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22r or ulaplacian23r 
                           uLap = ulaplacian23r(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2            ! residual in equation using current ghost value
                           a11 = 1./( dx(axis)**2 )                                  ! coeff of u(-1) in r1 
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                           ! if( j1==7 .and. j2==-1 )then
                           !   write(*,'(" BC2: j1,j2=",2i3," set u=",e10.2)') j1,j2,u(j1,j2,j3,0)
                           ! end if 
                         ! if( .false. )then
                         !   ! check the error
                         !   call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                         !   write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         ! end if
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ! This should have been done already
                         ! ghost = 1 
                         ! j1=i1-is1*ghost
                         ! j2=i2-is2*ghost
                         ! j3=i3-is3*ghost
                         ! u(j1,j2,j3,uc)=extrap3(u,i1,i2,i3,uc,is1,is2,is3)
                         ! #If "2" eq "2" 
                         !   if( numGhost.gt.1 )then
                         !     ! extrap second ghost (UPW)
                         !     ghost = 2
                         !     k1=i1-is1*ghost
                         !     k2=i2-is2*ghost
                         !     k3=i3-is3*ghost              
                         !     u(k1,k2,k3,uc)=extrap3(u,j1,j2,j3,uc,is1,is2,is3)
                         !   end if
                         ! #End
                       end if ! mask 
                     end do
                     end do
                     end do
                 end if
               else
                 if( nd.eq.2 )then
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     ! Try this: *wdh* Dec 2, 2023 : Exclude the ends when adjacent boundaries are dirichlet or exact:
                     ! getRestrictedLoopBounds(side,axis,extram, m1a,m1b,m2a,m2b,m3a,m3b)
                     ! if( t.le. 2*dt )then
                     !    write(*,'("dirichlet CBC: Stage I: side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i3)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
                     ! end if
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ghost =1 
                         j1=i1-is1*ghost; j2=i2-is2*ghost; j3=i3-is3*ghost
                         ! --- get the compatibility forcings for order=2 ---
                             ! No forcing, do nothing 
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22 or ulaplacian23 
                           uLap = ulaplacian22(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2
                             a11 = ( rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2 )/( dr(axis)**2 )  ! coeff of u(-1) in r1 
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                           ! if( j1==7 .and. j2==-1 )then
                           !   write(*,'(" BC2: j1,j2=",2i3," set u=",e10.2)') j1,j2,u(j1,j2,j3,0)
                           ! end if 
                         ! if( .false. )then
                         !   ! check the error
                         !   call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                         !   write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         ! end if
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ! This should have been done already
                         ! ghost = 1 
                         ! j1=i1-is1*ghost
                         ! j2=i2-is2*ghost
                         ! j3=i3-is3*ghost
                         ! u(j1,j2,j3,uc)=extrap3(u,i1,i2,i3,uc,is1,is2,is3)
                         ! #If "2" eq "2" 
                         !   if( numGhost.gt.1 )then
                         !     ! extrap second ghost (UPW)
                         !     ghost = 2
                         !     k1=i1-is1*ghost
                         !     k2=i2-is2*ghost
                         !     k3=i3-is3*ghost              
                         !     u(k1,k2,k3,uc)=extrap3(u,j1,j2,j3,uc,is1,is2,is3)
                         !   end if
                         ! #End
                       end if ! mask 
                     end do
                     end do
                     end do
                 else
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     ! Try this: *wdh* Dec 2, 2023 : Exclude the ends when adjacent boundaries are dirichlet or exact:
                     ! getRestrictedLoopBounds(side,axis,extram, m1a,m1b,m2a,m2b,m3a,m3b)
                     ! if( t.le. 2*dt )then
                     !    write(*,'("dirichlet CBC: Stage I: side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i3)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
                     ! end if
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ghost =1 
                         j1=i1-is1*ghost; j2=i2-is2*ghost; j3=i3-is3*ghost
                         ! --- get the compatibility forcings for order=2 ---
                             ! No forcing, do nothing 
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22 or ulaplacian23 
                           uLap = ulaplacian23(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2
                             a11 = ( rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2 + rsxy(i1,i2,i3,axis,2)**2 )/( dr(axis)**2 )  ! coeff of u(-1) in r1 
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                           ! if( j1==7 .and. j2==-1 )then
                           !   write(*,'(" BC2: j1,j2=",2i3," set u=",e10.2)') j1,j2,u(j1,j2,j3,0)
                           ! end if 
                         ! if( .false. )then
                         !   ! check the error
                         !   call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                         !   write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         ! end if
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ! This should have been done already
                         ! ghost = 1 
                         ! j1=i1-is1*ghost
                         ! j2=i2-is2*ghost
                         ! j3=i3-is3*ghost
                         ! u(j1,j2,j3,uc)=extrap3(u,i1,i2,i3,uc,is1,is2,is3)
                         ! #If "2" eq "2" 
                         !   if( numGhost.gt.1 )then
                         !     ! extrap second ghost (UPW)
                         !     ghost = 2
                         !     k1=i1-is1*ghost
                         !     k2=i2-is2*ghost
                         !     k3=i3-is3*ghost              
                         !     u(k1,k2,k3,uc)=extrap3(u,j1,j2,j3,uc,is1,is2,is3)
                         !   end if
                         ! #End
                       end if ! mask 
                     end do
                     end do
                     end do
                 end if
               end if 
             else
               if( gridType.eq.rectangular )then
                 if( nd.eq.2 )then
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     ! Try this: *wdh* Dec 2, 2023 : Exclude the ends when adjacent boundaries are dirichlet or exact:
                     ! getRestrictedLoopBounds(side,axis,extram, m1a,m1b,m2a,m2b,m3a,m3b)
                     ! if( t.le. 2*dt )then
                     !    write(*,'("dirichlet CBC: Stage I: side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i3)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
                     ! end if
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ghost =1 
                         j1=i1-is1*ghost; j2=i2-is2*ghost; j3=i3-is3*ghost
                         ! --- get the compatibility forcings for order=2 ---
                             if( assignTwilightZone.eq.1 )then
                               ! compute RHS from TZ
                                 call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexx )
                                 call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyy )
                                 ueLap = uexx + ueyy
                                 call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uett )
                               ff = uett - c2*ueLap
                               gtt = uett
                             else
                               ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.;
                             end if
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22r or ulaplacian23r 
                           uLap = ulaplacian22r(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2            ! residual in equation using current ghost value
                           a11 = 1./( dx(axis)**2 )                                  ! coeff of u(-1) in r1 
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                           ! if( j1==7 .and. j2==-1 )then
                           !   write(*,'(" BC2: j1,j2=",2i3," set u=",e10.2)') j1,j2,u(j1,j2,j3,0)
                           ! end if 
                         ! if( .false. )then
                         !   ! check the error
                         !   call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                         !   write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         ! end if
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ! This should have been done already
                         ! ghost = 1 
                         ! j1=i1-is1*ghost
                         ! j2=i2-is2*ghost
                         ! j3=i3-is3*ghost
                         ! u(j1,j2,j3,uc)=extrap3(u,i1,i2,i3,uc,is1,is2,is3)
                         ! #If "2" eq "2" 
                         !   if( numGhost.gt.1 )then
                         !     ! extrap second ghost (UPW)
                         !     ghost = 2
                         !     k1=i1-is1*ghost
                         !     k2=i2-is2*ghost
                         !     k3=i3-is3*ghost              
                         !     u(k1,k2,k3,uc)=extrap3(u,j1,j2,j3,uc,is1,is2,is3)
                         !   end if
                         ! #End
                       end if ! mask 
                     end do
                     end do
                     end do
                 else
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     ! Try this: *wdh* Dec 2, 2023 : Exclude the ends when adjacent boundaries are dirichlet or exact:
                     ! getRestrictedLoopBounds(side,axis,extram, m1a,m1b,m2a,m2b,m3a,m3b)
                     ! if( t.le. 2*dt )then
                     !    write(*,'("dirichlet CBC: Stage I: side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i3)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
                     ! end if
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ghost =1 
                         j1=i1-is1*ghost; j2=i2-is2*ghost; j3=i3-is3*ghost
                         ! --- get the compatibility forcings for order=2 ---
                             if( assignTwilightZone.eq.1 )then
                               ! compute RHS from TZ
                                 ! 3D 
                                 call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexx )
                                 call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyy )
                                 call ogDeriv(ep,0,0,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezz )
                                 ueLap = uexx + ueyy + uezz
                                 call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uett )
                               ff = uett - c2*ueLap
                               gtt = uett
                             else
                               ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.;
                             end if
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22r or ulaplacian23r 
                           uLap = ulaplacian23r(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2            ! residual in equation using current ghost value
                           a11 = 1./( dx(axis)**2 )                                  ! coeff of u(-1) in r1 
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                           ! if( j1==7 .and. j2==-1 )then
                           !   write(*,'(" BC2: j1,j2=",2i3," set u=",e10.2)') j1,j2,u(j1,j2,j3,0)
                           ! end if 
                         ! if( .false. )then
                         !   ! check the error
                         !   call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                         !   write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         ! end if
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ! This should have been done already
                         ! ghost = 1 
                         ! j1=i1-is1*ghost
                         ! j2=i2-is2*ghost
                         ! j3=i3-is3*ghost
                         ! u(j1,j2,j3,uc)=extrap3(u,i1,i2,i3,uc,is1,is2,is3)
                         ! #If "2" eq "2" 
                         !   if( numGhost.gt.1 )then
                         !     ! extrap second ghost (UPW)
                         !     ghost = 2
                         !     k1=i1-is1*ghost
                         !     k2=i2-is2*ghost
                         !     k3=i3-is3*ghost              
                         !     u(k1,k2,k3,uc)=extrap3(u,j1,j2,j3,uc,is1,is2,is3)
                         !   end if
                         ! #End
                       end if ! mask 
                     end do
                     end do
                     end do
                 end if
               else
                 if( nd.eq.2 )then
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     ! Try this: *wdh* Dec 2, 2023 : Exclude the ends when adjacent boundaries are dirichlet or exact:
                     ! getRestrictedLoopBounds(side,axis,extram, m1a,m1b,m2a,m2b,m3a,m3b)
                     ! if( t.le. 2*dt )then
                     !    write(*,'("dirichlet CBC: Stage I: side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i3)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
                     ! end if
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ghost =1 
                         j1=i1-is1*ghost; j2=i2-is2*ghost; j3=i3-is3*ghost
                         ! --- get the compatibility forcings for order=2 ---
                             if( assignTwilightZone.eq.1 )then
                               ! compute RHS from TZ
                                 call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexx )
                                 call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyy )
                                 ueLap = uexx + ueyy
                                 call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uett )
                               ff = uett - c2*ueLap
                               gtt = uett
                             else
                               ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.;
                             end if
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22 or ulaplacian23 
                           uLap = ulaplacian22(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2
                             a11 = ( rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2 )/( dr(axis)**2 )  ! coeff of u(-1) in r1 
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                           ! if( j1==7 .and. j2==-1 )then
                           !   write(*,'(" BC2: j1,j2=",2i3," set u=",e10.2)') j1,j2,u(j1,j2,j3,0)
                           ! end if 
                         ! if( .false. )then
                         !   ! check the error
                         !   call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                         !   write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         ! end if
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ! This should have been done already
                         ! ghost = 1 
                         ! j1=i1-is1*ghost
                         ! j2=i2-is2*ghost
                         ! j3=i3-is3*ghost
                         ! u(j1,j2,j3,uc)=extrap3(u,i1,i2,i3,uc,is1,is2,is3)
                         ! #If "2" eq "2" 
                         !   if( numGhost.gt.1 )then
                         !     ! extrap second ghost (UPW)
                         !     ghost = 2
                         !     k1=i1-is1*ghost
                         !     k2=i2-is2*ghost
                         !     k3=i3-is3*ghost              
                         !     u(k1,k2,k3,uc)=extrap3(u,j1,j2,j3,uc,is1,is2,is3)
                         !   end if
                         ! #End
                       end if ! mask 
                     end do
                     end do
                     end do
                 else
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     ! Try this: *wdh* Dec 2, 2023 : Exclude the ends when adjacent boundaries are dirichlet or exact:
                     ! getRestrictedLoopBounds(side,axis,extram, m1a,m1b,m2a,m2b,m3a,m3b)
                     ! if( t.le. 2*dt )then
                     !    write(*,'("dirichlet CBC: Stage I: side,axis=",2i3," m1a,m1b,m2a,m2b,m3a,m3b=",6i3)') side,axis,m1a,m1b,m2a,m2b,m3a,m3b
                     ! end if
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ghost =1 
                         j1=i1-is1*ghost; j2=i2-is2*ghost; j3=i3-is3*ghost
                         ! --- get the compatibility forcings for order=2 ---
                             if( assignTwilightZone.eq.1 )then
                               ! compute RHS from TZ
                                 ! 3D 
                                 call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexx )
                                 call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyy )
                                 call ogDeriv(ep,0,0,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezz )
                                 ueLap = uexx + ueyy + uezz
                                 call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uett )
                               ff = uett - c2*ueLap
                               gtt = uett
                             else
                               ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.;
                             end if
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22 or ulaplacian23 
                           uLap = ulaplacian23(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2
                             a11 = ( rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2 + rsxy(i1,i2,i3,axis,2)**2 )/( dr(axis)**2 )  ! coeff of u(-1) in r1 
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                           ! if( j1==7 .and. j2==-1 )then
                           !   write(*,'(" BC2: j1,j2=",2i3," set u=",e10.2)') j1,j2,u(j1,j2,j3,0)
                           ! end if 
                         ! if( .false. )then
                         !   ! check the error
                         !   call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                         !   write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         ! end if
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ! This should have been done already
                         ! ghost = 1 
                         ! j1=i1-is1*ghost
                         ! j2=i2-is2*ghost
                         ! j3=i3-is3*ghost
                         ! u(j1,j2,j3,uc)=extrap3(u,i1,i2,i3,uc,is1,is2,is3)
                         ! #If "2" eq "2" 
                         !   if( numGhost.gt.1 )then
                         !     ! extrap second ghost (UPW)
                         !     ghost = 2
                         !     k1=i1-is1*ghost
                         !     k2=i2-is2*ghost
                         !     k3=i3-is3*ghost              
                         !     u(k1,k2,k3,uc)=extrap3(u,j1,j2,j3,uc,is1,is2,is3)
                         !   end if
                         ! #End
                       end if ! mask 
                     end do
                     end do
                     end do
                 end if
               end if     
             end if
       else if( orderOfAccuracy.eq.4 )then
       else if( orderOfAccuracy.eq.6 )then  
         ! turned off May 4, 2023
         write(*,'("CgWave::bcOpt:ERROR: CBC order=6 turned off, use LCBC")')
         stop 666
         ! callDirichletGhostCompatibility(6)
       else if( orderOfAccuracy.eq.8 )then   
         stop 888
       else
         write(*,'("CgWave::bcOpt:ERROR:Dirichlet CBC unexpected orderOfAccuracy=",i6)') orderOfAccuracy
         stop 8888
       end if
     else if( bc(side,axis).eq.neumann )then
       ! ------ NEUMANN ----------
       if( bcApproach.eq.useCompatibilityBoundaryConditions )then
         if( orderOfAccuracy.eq.2 )then
             if( forcingOption.eq.noForcing )then
               if( gridType.eq.rectangular )then
                 if( nd.eq.2 )then
                     ! write(*,'("START ASSIGN NEUMANN GHOST COMPATIBILITY dim=[2] gridType=[rectangular] forcing=[noForcing] ======")') 
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                       ! compute the outward normal (an1,an2,an3)
                       an1 = 0.
                       an2 = 0.
                       an3 = 0.
                       if( axis.eq.0 )then
                        an1=-is
                       else if( axis.eq.1 )then
                        an2=-is
                       else
                        an3=-is
                       end if
                       dxn=dx(axis)
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                   ! *wdh* Dec 7 : Apply 2nd order conditions for fourth order scheme too 
                   ! *wdh* Dec 7, 2023 #If "2" eq "2"
                     if( orderOfAccuracy==4 )then
                       if( t.le.3*dt )then
                         write(*,'("Neumann CBC order=4: Stage I: apply 2nd order conditions first")') 
                       end if
                     end if
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     gg=0.; nDotGradF=0.; gtt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       ! first ghost pt:
                       j1=i1-is1; j2=i2-is2; j3=i3-is3
                       if( mask(i1,i2,i3).gt.0 )then
                           ! write(*,'("getNeumannCompatForcing forcing=noForcing dim=2 order=2 assignTwilightZone=",i2)') assignTwilightZone
                             ! No forcing, do nothing 
                             gg=0.;  gtt=0.; nDotGradF=0.; 
                           ! --- NEUMANN 2=2 rectangular ---
                           !    (+/-)*a1*[-u(i-1)  + u(i+1)] + 2*dxn*a0*u(i) = 2*dxn*gg 
                           !  *check me* 
                           b0 = -2.*dxn*a0/a1 
                           b1 =  2.*dxn/a1 
                           u(j1,j2,j3,uc)=  b0*u(j1+  is1,j2+  is2,j3+  is3,uc)+    u(j1+2*is1,j2+2*is2,j3+2*is3,uc)+ b1*gg
                           ! ----- Assign extra ghost ----
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                           if( numGhost.gt.1 )then
                             ! extrap second ghost (UPW)
                             ghost = 2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost              
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))
                           end if
                       end if
                      end do
                      end do
                      end do
                   ! Dec 7, 2023 #End
                 else
                     ! write(*,'("START ASSIGN NEUMANN GHOST COMPATIBILITY dim=[3] gridType=[rectangular] forcing=[noForcing] ======")') 
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                       ! compute the outward normal (an1,an2,an3)
                       an1 = 0.
                       an2 = 0.
                       an3 = 0.
                       if( axis.eq.0 )then
                        an1=-is
                       else if( axis.eq.1 )then
                        an2=-is
                       else
                        an3=-is
                       end if
                       dxn=dx(axis)
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                   ! *wdh* Dec 7 : Apply 2nd order conditions for fourth order scheme too 
                   ! *wdh* Dec 7, 2023 #If "2" eq "2"
                     if( orderOfAccuracy==4 )then
                       if( t.le.3*dt )then
                         write(*,'("Neumann CBC order=4: Stage I: apply 2nd order conditions first")') 
                       end if
                     end if
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     gg=0.; nDotGradF=0.; gtt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       ! first ghost pt:
                       j1=i1-is1; j2=i2-is2; j3=i3-is3
                       if( mask(i1,i2,i3).gt.0 )then
                           ! write(*,'("getNeumannCompatForcing forcing=noForcing dim=3 order=2 assignTwilightZone=",i2)') assignTwilightZone
                             ! No forcing, do nothing 
                             gg=0.;  gtt=0.; nDotGradF=0.; 
                           ! --- NEUMANN 2=2 rectangular ---
                           !    (+/-)*a1*[-u(i-1)  + u(i+1)] + 2*dxn*a0*u(i) = 2*dxn*gg 
                           !  *check me* 
                           b0 = -2.*dxn*a0/a1 
                           b1 =  2.*dxn/a1 
                           u(j1,j2,j3,uc)=  b0*u(j1+  is1,j2+  is2,j3+  is3,uc)+    u(j1+2*is1,j2+2*is2,j3+2*is3,uc)+ b1*gg
                           ! ----- Assign extra ghost ----
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                           if( numGhost.gt.1 )then
                             ! extrap second ghost (UPW)
                             ghost = 2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost              
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))
                           end if
                       end if
                      end do
                      end do
                      end do
                   ! Dec 7, 2023 #End
                 end if
               else
                 if( nd.eq.2 )then
                     ! write(*,'("START ASSIGN NEUMANN GHOST COMPATIBILITY dim=[2] gridType=[curvilinear] forcing=[noForcing] ======")') 
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                   ! *wdh* Dec 7 : Apply 2nd order conditions for fourth order scheme too 
                   ! *wdh* Dec 7, 2023 #If "2" eq "2"
                     if( orderOfAccuracy==4 )then
                       if( t.le.3*dt )then
                         write(*,'("Neumann CBC order=4: Stage I: apply 2nd order conditions first")') 
                       end if
                     end if
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     gg=0.; nDotGradF=0.; gtt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       ! first ghost pt:
                       j1=i1-is1; j2=i2-is2; j3=i3-is3
                       if( mask(i1,i2,i3).gt.0 )then
                           ! compute the outward normal (an1,an2,an3)
                               an1 = rsxy(i1,i2,i3,axis,0)
                               an2 = rsxy(i1,i2,i3,axis,1)
                               if( nd.eq.2 )then
                                aNormi = (-is)/sqrt(an1**2+an2**2)
                                an1=an1*aNormi
                                an2=an2*aNormi
                               else
                                an3 = rsxy(i1,i2,i3,axis,2)
                                aNormi = (-is)/sqrt(an1**2+an2**2+an3**2)
                                an1=an1*aNormi
                                an2=an2*aNormi
                                an3=an3*aNormi
                               end if
                           ! write(*,'("getNeumannCompatForcing forcing=noForcing dim=2 order=2 assignTwilightZone=",i2)') assignTwilightZone
                             ! No forcing, do nothing 
                             gg=0.;  gtt=0.; nDotGradF=0.; 
                           ! ------ curvilinear grid: -------
                           ! a1*( n1*ux + n2*ux + n3*uz ) + a0*u = f 
                           ! a1*( (n1*rx+n2*ry+n3*rz)*ur + (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ) + a0*u = f 
                           ! =>
                           !  ur = [ f - (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ]/( a1*( (n1*rx+n2*ry+n3*rz) ) 
                           ! ----- NEUMANN 2=2 curvilinear ----
                           ! ur = ( u(i+1) - u(i-1) )/2*dr
                           ! ur = ur0 -> 
                           ! u(i-1) = u(i+1) - 2*dr*( ur0 )   (left)
                           ! u(i+1) = u(i-1) + 2*dr*( ur0 )   (right) 
                           urv(0) = ur2(i1,i2,i3,uc)
                           urv(1) = us2(i1,i2,i3,uc)
                             t1=a1*( an1*rsxy(i1,i2,i3,axis  ,0)+an2*rsxy(i1,i2,i3,axis  ,1) )
                             t2=a1*( an1*rsxy(i1,i2,i3,axisp1,0)+an2*rsxy(i1,i2,i3,axisp1,1) ) 
                             ur0 = (gg - ( t2*urv(axisp1) + a0*u(i1,i2,i3,uc) ) )/t1
                           u(j1,j2,j3,uc) =  u(j1+2*is1,j2+2*is2,j3+2*is3,uc) -2.*is*dr(axis)*ur0
                           ! write(*,'("neumann CBC j1,j2=",2i4," u=",1pe12.4," gg=",1pe10.2)') j1,j2,u(j1,j2,j3,uc),gg
                           ! ----- Assign extra ghost ----
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                           if( numGhost.gt.1 )then
                             ! extrap second ghost (UPW)
                             ghost = 2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost              
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))
                           end if
                       end if
                      end do
                      end do
                      end do
                   ! Dec 7, 2023 #End
                 else
                     ! write(*,'("START ASSIGN NEUMANN GHOST COMPATIBILITY dim=[3] gridType=[curvilinear] forcing=[noForcing] ======")') 
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                   ! *wdh* Dec 7 : Apply 2nd order conditions for fourth order scheme too 
                   ! *wdh* Dec 7, 2023 #If "2" eq "2"
                     if( orderOfAccuracy==4 )then
                       if( t.le.3*dt )then
                         write(*,'("Neumann CBC order=4: Stage I: apply 2nd order conditions first")') 
                       end if
                     end if
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     gg=0.; nDotGradF=0.; gtt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       ! first ghost pt:
                       j1=i1-is1; j2=i2-is2; j3=i3-is3
                       if( mask(i1,i2,i3).gt.0 )then
                           ! compute the outward normal (an1,an2,an3)
                               an1 = rsxy(i1,i2,i3,axis,0)
                               an2 = rsxy(i1,i2,i3,axis,1)
                               if( nd.eq.2 )then
                                aNormi = (-is)/sqrt(an1**2+an2**2)
                                an1=an1*aNormi
                                an2=an2*aNormi
                               else
                                an3 = rsxy(i1,i2,i3,axis,2)
                                aNormi = (-is)/sqrt(an1**2+an2**2+an3**2)
                                an1=an1*aNormi
                                an2=an2*aNormi
                                an3=an3*aNormi
                               end if
                           ! write(*,'("getNeumannCompatForcing forcing=noForcing dim=3 order=2 assignTwilightZone=",i2)') assignTwilightZone
                             ! No forcing, do nothing 
                             gg=0.;  gtt=0.; nDotGradF=0.; 
                           ! ------ curvilinear grid: -------
                           ! a1*( n1*ux + n2*ux + n3*uz ) + a0*u = f 
                           ! a1*( (n1*rx+n2*ry+n3*rz)*ur + (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ) + a0*u = f 
                           ! =>
                           !  ur = [ f - (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ]/( a1*( (n1*rx+n2*ry+n3*rz) ) 
                           ! ----- NEUMANN 2=2 curvilinear ----
                           ! ur = ( u(i+1) - u(i-1) )/2*dr
                           ! ur = ur0 -> 
                           ! u(i-1) = u(i+1) - 2*dr*( ur0 )   (left)
                           ! u(i+1) = u(i-1) + 2*dr*( ur0 )   (right) 
                           urv(0) = ur2(i1,i2,i3,uc)
                           urv(1) = us2(i1,i2,i3,uc)
                             urv(2) = ut2(i1,i2,i3,uc)
                             t1=a1*( an1*rsxy(i1,i2,i3,axis  ,0)+an2*rsxy(i1,i2,i3,axis  ,1)+an3*rsxy(i1,i2,i3,axis  ,2) )
                             t2=a1*( an1*rsxy(i1,i2,i3,axisp1,0)+an2*rsxy(i1,i2,i3,axisp1,1)+an3*rsxy(i1,i2,i3,axisp1,2) )
                             t3=a1*( an1*rsxy(i1,i2,i3,axisp2,0)+an2*rsxy(i1,i2,i3,axisp2,1)+an3*rsxy(i1,i2,i3,axisp2,2) )
                             ur0 = ( gg - ( t2*urv(axisp1) + t3*urv(axisp2) + a0*u(i1,i2,i3,uc) ) )/t1
                           u(j1,j2,j3,uc) =  u(j1+2*is1,j2+2*is2,j3+2*is3,uc) -2.*is*dr(axis)*ur0
                           ! write(*,'("neumann CBC j1,j2=",2i4," u=",1pe12.4," gg=",1pe10.2)') j1,j2,u(j1,j2,j3,uc),gg
                           ! ----- Assign extra ghost ----
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                           if( numGhost.gt.1 )then
                             ! extrap second ghost (UPW)
                             ghost = 2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost              
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))
                           end if
                       end if
                      end do
                      end do
                      end do
                   ! Dec 7, 2023 #End
                 end if
               end if 
             else
               if( gridType.eq.rectangular )then
                 if( nd.eq.2 )then
                     ! write(*,'("START ASSIGN NEUMANN GHOST COMPATIBILITY dim=[2] gridType=[rectangular] forcing=[forcing] ======")') 
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                       ! compute the outward normal (an1,an2,an3)
                       an1 = 0.
                       an2 = 0.
                       an3 = 0.
                       if( axis.eq.0 )then
                        an1=-is
                       else if( axis.eq.1 )then
                        an2=-is
                       else
                        an3=-is
                       end if
                       dxn=dx(axis)
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                   ! *wdh* Dec 7 : Apply 2nd order conditions for fourth order scheme too 
                   ! *wdh* Dec 7, 2023 #If "2" eq "2"
                     if( orderOfAccuracy==4 )then
                       if( t.le.3*dt )then
                         write(*,'("Neumann CBC order=4: Stage I: apply 2nd order conditions first")') 
                       end if
                     end if
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     gg=0.; nDotGradF=0.; gtt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       ! first ghost pt:
                       j1=i1-is1; j2=i2-is2; j3=i3-is3
                       if( mask(i1,i2,i3).gt.0 )then
                           ! write(*,'("getNeumannCompatForcing forcing=forcing dim=2 order=2 assignTwilightZone=",i2)') assignTwilightZone
                             if( assignTwilightZone.eq.1 )then
                               ! compute RHS from TZ
                                 call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uex )
                                 call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uey )
                                 gg = an1*uex + an2*uey
                             else
                               gg=0.;  gtt=0.; nDotGradF=0.; 
                             end if
                           ! --- NEUMANN 2=2 rectangular ---
                           !    (+/-)*a1*[-u(i-1)  + u(i+1)] + 2*dxn*a0*u(i) = 2*dxn*gg 
                           !  *check me* 
                           b0 = -2.*dxn*a0/a1 
                           b1 =  2.*dxn/a1 
                           u(j1,j2,j3,uc)=  b0*u(j1+  is1,j2+  is2,j3+  is3,uc)+    u(j1+2*is1,j2+2*is2,j3+2*is3,uc)+ b1*gg
                           ! ----- Assign extra ghost ----
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                           if( numGhost.gt.1 )then
                             ! extrap second ghost (UPW)
                             ghost = 2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost              
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))
                           end if
                       end if
                      end do
                      end do
                      end do
                   ! Dec 7, 2023 #End
                 else
                     ! write(*,'("START ASSIGN NEUMANN GHOST COMPATIBILITY dim=[3] gridType=[rectangular] forcing=[forcing] ======")') 
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                       ! compute the outward normal (an1,an2,an3)
                       an1 = 0.
                       an2 = 0.
                       an3 = 0.
                       if( axis.eq.0 )then
                        an1=-is
                       else if( axis.eq.1 )then
                        an2=-is
                       else
                        an3=-is
                       end if
                       dxn=dx(axis)
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                   ! *wdh* Dec 7 : Apply 2nd order conditions for fourth order scheme too 
                   ! *wdh* Dec 7, 2023 #If "2" eq "2"
                     if( orderOfAccuracy==4 )then
                       if( t.le.3*dt )then
                         write(*,'("Neumann CBC order=4: Stage I: apply 2nd order conditions first")') 
                       end if
                     end if
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     gg=0.; nDotGradF=0.; gtt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       ! first ghost pt:
                       j1=i1-is1; j2=i2-is2; j3=i3-is3
                       if( mask(i1,i2,i3).gt.0 )then
                           ! write(*,'("getNeumannCompatForcing forcing=forcing dim=3 order=2 assignTwilightZone=",i2)') assignTwilightZone
                             if( assignTwilightZone.eq.1 )then
                               ! compute RHS from TZ
                                 ! ----- 3D  -----
                                 call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uex )
                                 call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uey )
                                 call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uez )
                                 gg = an1*uex + an2*uey + an3*uez
                             else
                               gg=0.;  gtt=0.; nDotGradF=0.; 
                             end if
                           ! --- NEUMANN 2=2 rectangular ---
                           !    (+/-)*a1*[-u(i-1)  + u(i+1)] + 2*dxn*a0*u(i) = 2*dxn*gg 
                           !  *check me* 
                           b0 = -2.*dxn*a0/a1 
                           b1 =  2.*dxn/a1 
                           u(j1,j2,j3,uc)=  b0*u(j1+  is1,j2+  is2,j3+  is3,uc)+    u(j1+2*is1,j2+2*is2,j3+2*is3,uc)+ b1*gg
                           ! ----- Assign extra ghost ----
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                           if( numGhost.gt.1 )then
                             ! extrap second ghost (UPW)
                             ghost = 2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost              
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))
                           end if
                       end if
                      end do
                      end do
                      end do
                   ! Dec 7, 2023 #End
                 end if
               else
                 if( nd.eq.2 )then
                     ! write(*,'("START ASSIGN NEUMANN GHOST COMPATIBILITY dim=[2] gridType=[curvilinear] forcing=[forcing] ======")') 
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                   ! *wdh* Dec 7 : Apply 2nd order conditions for fourth order scheme too 
                   ! *wdh* Dec 7, 2023 #If "2" eq "2"
                     if( orderOfAccuracy==4 )then
                       if( t.le.3*dt )then
                         write(*,'("Neumann CBC order=4: Stage I: apply 2nd order conditions first")') 
                       end if
                     end if
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     gg=0.; nDotGradF=0.; gtt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       ! first ghost pt:
                       j1=i1-is1; j2=i2-is2; j3=i3-is3
                       if( mask(i1,i2,i3).gt.0 )then
                           ! compute the outward normal (an1,an2,an3)
                               an1 = rsxy(i1,i2,i3,axis,0)
                               an2 = rsxy(i1,i2,i3,axis,1)
                               if( nd.eq.2 )then
                                aNormi = (-is)/sqrt(an1**2+an2**2)
                                an1=an1*aNormi
                                an2=an2*aNormi
                               else
                                an3 = rsxy(i1,i2,i3,axis,2)
                                aNormi = (-is)/sqrt(an1**2+an2**2+an3**2)
                                an1=an1*aNormi
                                an2=an2*aNormi
                                an3=an3*aNormi
                               end if
                           ! write(*,'("getNeumannCompatForcing forcing=forcing dim=2 order=2 assignTwilightZone=",i2)') assignTwilightZone
                             if( assignTwilightZone.eq.1 )then
                               ! compute RHS from TZ
                                 call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uex )
                                 call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uey )
                                 gg = an1*uex + an2*uey
                             else
                               gg=0.;  gtt=0.; nDotGradF=0.; 
                             end if
                           ! ------ curvilinear grid: -------
                           ! a1*( n1*ux + n2*ux + n3*uz ) + a0*u = f 
                           ! a1*( (n1*rx+n2*ry+n3*rz)*ur + (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ) + a0*u = f 
                           ! =>
                           !  ur = [ f - (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ]/( a1*( (n1*rx+n2*ry+n3*rz) ) 
                           ! ----- NEUMANN 2=2 curvilinear ----
                           ! ur = ( u(i+1) - u(i-1) )/2*dr
                           ! ur = ur0 -> 
                           ! u(i-1) = u(i+1) - 2*dr*( ur0 )   (left)
                           ! u(i+1) = u(i-1) + 2*dr*( ur0 )   (right) 
                           urv(0) = ur2(i1,i2,i3,uc)
                           urv(1) = us2(i1,i2,i3,uc)
                             t1=a1*( an1*rsxy(i1,i2,i3,axis  ,0)+an2*rsxy(i1,i2,i3,axis  ,1) )
                             t2=a1*( an1*rsxy(i1,i2,i3,axisp1,0)+an2*rsxy(i1,i2,i3,axisp1,1) ) 
                             ur0 = (gg - ( t2*urv(axisp1) + a0*u(i1,i2,i3,uc) ) )/t1
                           u(j1,j2,j3,uc) =  u(j1+2*is1,j2+2*is2,j3+2*is3,uc) -2.*is*dr(axis)*ur0
                           ! write(*,'("neumann CBC j1,j2=",2i4," u=",1pe12.4," gg=",1pe10.2)') j1,j2,u(j1,j2,j3,uc),gg
                           ! ----- Assign extra ghost ----
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                           if( numGhost.gt.1 )then
                             ! extrap second ghost (UPW)
                             ghost = 2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost              
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))
                           end if
                       end if
                      end do
                      end do
                      end do
                   ! Dec 7, 2023 #End
                 else
                     ! write(*,'("START ASSIGN NEUMANN GHOST COMPATIBILITY dim=[3] gridType=[curvilinear] forcing=[forcing] ======")') 
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                   ! *wdh* Dec 7 : Apply 2nd order conditions for fourth order scheme too 
                   ! *wdh* Dec 7, 2023 #If "2" eq "2"
                     if( orderOfAccuracy==4 )then
                       if( t.le.3*dt )then
                         write(*,'("Neumann CBC order=4: Stage I: apply 2nd order conditions first")') 
                       end if
                     end if
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       m1a=gridIndexRange(0,0)-extram
                       m1b=gridIndexRange(1,0)+extram
                       m2a=gridIndexRange(0,1)-extram
                       m2b=gridIndexRange(1,1)+extram
                       if( nd.eq.2 )then
                         m3a=gridIndexRange(0,2)
                         m3b=gridIndexRange(1,2)
                       else
                         m3a=gridIndexRange(0,2)-extram
                         m3b=gridIndexRange(1,2)+extram
                       end if
                       if( axis.eq.0 )then
                        m1a=gridIndexRange(side,axis)
                        m1b=gridIndexRange(side,axis)
                       else if( axis.eq.1 )then
                        m2a=gridIndexRange(side,axis)
                        m2b=gridIndexRange(side,axis)
                       else
                        m3a=gridIndexRange(side,axis)
                        m3b=gridIndexRange(side,axis)
                       end if
                     gg=0.; nDotGradF=0.; gtt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       ! first ghost pt:
                       j1=i1-is1; j2=i2-is2; j3=i3-is3
                       if( mask(i1,i2,i3).gt.0 )then
                           ! compute the outward normal (an1,an2,an3)
                               an1 = rsxy(i1,i2,i3,axis,0)
                               an2 = rsxy(i1,i2,i3,axis,1)
                               if( nd.eq.2 )then
                                aNormi = (-is)/sqrt(an1**2+an2**2)
                                an1=an1*aNormi
                                an2=an2*aNormi
                               else
                                an3 = rsxy(i1,i2,i3,axis,2)
                                aNormi = (-is)/sqrt(an1**2+an2**2+an3**2)
                                an1=an1*aNormi
                                an2=an2*aNormi
                                an3=an3*aNormi
                               end if
                           ! write(*,'("getNeumannCompatForcing forcing=forcing dim=3 order=2 assignTwilightZone=",i2)') assignTwilightZone
                             if( assignTwilightZone.eq.1 )then
                               ! compute RHS from TZ
                                 ! ----- 3D  -----
                                 call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uex )
                                 call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uey )
                                 call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uez )
                                 gg = an1*uex + an2*uey + an3*uez
                             else
                               gg=0.;  gtt=0.; nDotGradF=0.; 
                             end if
                           ! ------ curvilinear grid: -------
                           ! a1*( n1*ux + n2*ux + n3*uz ) + a0*u = f 
                           ! a1*( (n1*rx+n2*ry+n3*rz)*ur + (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ) + a0*u = f 
                           ! =>
                           !  ur = [ f - (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ]/( a1*( (n1*rx+n2*ry+n3*rz) ) 
                           ! ----- NEUMANN 2=2 curvilinear ----
                           ! ur = ( u(i+1) - u(i-1) )/2*dr
                           ! ur = ur0 -> 
                           ! u(i-1) = u(i+1) - 2*dr*( ur0 )   (left)
                           ! u(i+1) = u(i-1) + 2*dr*( ur0 )   (right) 
                           urv(0) = ur2(i1,i2,i3,uc)
                           urv(1) = us2(i1,i2,i3,uc)
                             urv(2) = ut2(i1,i2,i3,uc)
                             t1=a1*( an1*rsxy(i1,i2,i3,axis  ,0)+an2*rsxy(i1,i2,i3,axis  ,1)+an3*rsxy(i1,i2,i3,axis  ,2) )
                             t2=a1*( an1*rsxy(i1,i2,i3,axisp1,0)+an2*rsxy(i1,i2,i3,axisp1,1)+an3*rsxy(i1,i2,i3,axisp1,2) )
                             t3=a1*( an1*rsxy(i1,i2,i3,axisp2,0)+an2*rsxy(i1,i2,i3,axisp2,1)+an3*rsxy(i1,i2,i3,axisp2,2) )
                             ur0 = ( gg - ( t2*urv(axisp1) + t3*urv(axisp2) + a0*u(i1,i2,i3,uc) ) )/t1
                           u(j1,j2,j3,uc) =  u(j1+2*is1,j2+2*is2,j3+2*is3,uc) -2.*is*dr(axis)*ur0
                           ! write(*,'("neumann CBC j1,j2=",2i4," u=",1pe12.4," gg=",1pe10.2)') j1,j2,u(j1,j2,j3,uc),gg
                           ! ----- Assign extra ghost ----
                           if( numGhost.gt.1 )then
                             ghost =2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost   
                             ! extrap second ghost (UPW)
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                           end if
                       else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                           if( numGhost.gt.1 )then
                             ! extrap second ghost (UPW)
                             ghost = 2
                             k1=i1-is1*ghost
                             k2=i2-is2*ghost
                             k3=i3-is3*ghost              
                             u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))
                           end if
                       end if
                      end do
                      end do
                      end do
                   ! Dec 7, 2023 #End
                 end if
               end if     
             end if
         else if( orderOfAccuracy.eq.4 )then
         else if( orderOfAccuracy.eq.6 )then   
           stop 666
         else if( orderOfAccuracy.eq.8 )then   
           stop 888
         else
           write(*,'("CgWave::bcOpt:ERROR:neumann CBC unexpected orderOfAccuracy=",i6)') orderOfAccuracy
          stop 8888
         end if
       else ! one-sided 
         if( orderOfAccuracy.eq.2 )then
              ! BC: a0*T + a1*T.n = 
              ! a0=mixedCoeff(uc,side,axis,grid)
              ! a1=mixedNormalCoeff(uc,side,axis,grid)
              a0=0.
              a1=1.
              ! rectangular case:
              if( gridType.eq.rectangular )then
                ! compute the outward normal (an1,an2,an3)
                an1 = 0.
                an2 = 0.
                an3 = 0.
                if( axis.eq.0 )then
                 an1=-is
                else if( axis.eq.1 )then
                 an2=-is
                else
                 an3=-is
                end if
                dxn=dx(axis)
                b0=-4.*dxn*a0/a1-10./3.
                b1=4.*(dxn/a1)
              end if
              ff=0.
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
               ! first ghost pt:
               j1=i1-is1
               j2=i2-is2
               j3=i3-is3
               ! 2nd ghost:
               k1=j1-is1
               k2=j2-is2
               k3=j3-is3
               ! 3rd ghost:
               l1=k1-is1
               l2=k2-is2
               l3=k3-is3    
               if( mask(i1,i2,i3).gt.0 )then
                 if( gridType.eq.curvilinear )then
                   ! compute the outward normal (an1,an2,an3)
                       an1 = rsxy(i1,i2,i3,axis,0)
                       an2 = rsxy(i1,i2,i3,axis,1)
                       if( nd.eq.2 )then
                        aNormi = (-is)/sqrt(an1**2+an2**2)
                        an1=an1*aNormi
                        an2=an2*aNormi
                       else
                        an3 = rsxy(i1,i2,i3,axis,2)
                        aNormi = (-is)/sqrt(an1**2+an2**2+an3**2)
                        an1=an1*aNormi
                        an2=an2*aNormi
                        an3=an3*aNormi
                       end if
                 end if
                   if( assignTwilightZone.eq.1 )then
                     ! compute RHS from TZ
                     if( nd.eq.2 )then
                       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                       call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uex)
                       call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uey)
                       ff = a0*ue + a1*( an1*uex + an2*uey )
                     else
                       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                       call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uex)
                       call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uey)
                       call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uez)
                       ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
                     end if
                   else if( assignKnownSolutionAtBoundaries.eq.1 )then
                     ! -- we set inhomogeneous Neumann values for some known solutions 
                     if( knownSolutionOption.eq.planeWave )then
                       ! --- evaluate RHS for the plane wave solution ---
                       if( nd.eq.2 )then
                         ue    = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                         cosPW = ampPlaneWave*cos( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                         uex   = kxPlaneWave*cosPW
                         uey   = kyPlaneWave*cosPw
                         ff = a0*ue + a1*( an1*uex + an2*uey )
                       else
                         ue    = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                         cosPW = ampPlaneWave*cos( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                         uex   = kxPlaneWave*cosPW
                         uey   = kyPlaneWave*cosPw
                         uez   = kzPlaneWave*cosPw
                         ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
                       end if 
                     else if( knownSolutionOption.eq.gaussianPlaneWave )then
                       ! Do nothing for Gaussian plane wave solution for now
                       ff = 0.
                     else if( knownSolutionOption.eq.boxHelmholtz ) then
                       ! --- evaluate RHS the boxHelmholtz solution ---
                       if( nd.eq.2 )then
                         ue  = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) * sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) *coswt
                         uex = cos( kxBoxHelmholtz*xy(i1,i2,i3,0) ) * sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) *coswt * kxBoxHelmholtz
                         uey = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) * cos( kyBoxHelmholtz*xy(i1,i2,i3,1) ) *coswt * kyBoxHelmholtz
                         ff = a0*ue + a1*( an1*uex + an2*uey )
                       else
                         ue  = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * sin( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt
                         uex = cos( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * sin( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt * kxBoxHelmholtz
                         uey = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *cos( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * sin( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt * kyBoxHelmholtz
                         uez = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * cos( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt * kzBoxHelmholtz
                         ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
                       end if
                     else if( knownSolutionOption.eq.polyPeriodic ) then
                       ! --- evaluate RHS the polyPeriodic solution ---
                       if( nd.eq.2 )then
                         ue  = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) ) *coswt
                         uex = (      a1PolyPeriodic                                                            ) *coswt
                         uey = (                          b1PolyPeriodic                                        ) *coswt
                         ff = a0*ue + a1*( an1*uex + an2*uey )
                       else
                         ue  = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                         uex = (      a1PolyPeriodic                                                                                            ) *coswt
                         uey = (                          b1PolyPeriodic                                                                        ) *coswt
                         uez = (                                              c1PolyPeriodic                                                    ) *coswt
                         ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
                       end if
                     else
                       stop 9876
                     end if 
                   end if
                 ! 2=4: 
                 ! --- assign 2 ghost points using:
                 !  (1) Apply Neumann BC to 4th order
                 !  (2) Extrap. 2nd ghost to 5th order
                 if( gridType.eq.rectangular )then
                   ! write(*,'(" TBC: j1,j2=",2i3," u,ff=",2e12.2)') j1,j2,ff,u(j1,j2,j3,uc)
                   !if( orderOfAccuracy.eq.2 )then
                     ! --- NEUMANN 2=2 rectangular ---
                     !    (+/-)*a1*[-u(i-1)  + u(i+1)] + 2*dxn*a0*u(i) = 2*dxn*ff 
                     !  *check me* 
                     b0 = -2.*dxn*a0/a1 
                     b1 =  2.*dxn/a1 
                     u(j1,j2,j3,uc)=  b0*u(j1+  is1,j2+  is2,j3+  is3,uc)+   u(j1+2*is1,j2+2*is2,j3+2*is3,uc)+ b1*ff
                 else 
                   ! ------ curvilinear grid: -------
                   ! a1*( n1*ux + n2*ux + n3*uz ) + a0*u = f 
                   ! a1*( (n1*rx+n2*ry+n3*rz)*ur + (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ) + a0*u = f 
                   ! =>
                   !  ur = [ f - (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ]/( a1*( (n1*rx+n2*ry+n3*rz) ) 
                     ! ----- NEUMANN 2=2 curvilinear ----
                     ! ur = ( u(i+1) - u(i-1) )/2*dr
                     ! ur = ur0 -> 
                     ! u(i-1) = u(i+1) - 2*dr*( ur0 )   (left)
                     ! u(i+1) = u(i-1) + 2*dr*( ur0 )   (right) 
                     urv(0) = ur2(i1,i2,i3,uc)
                     urv(1) = us2(i1,i2,i3,uc)
                     if( nd.eq.2 )then
                       t1=a1*( an1*rsxy(i1,i2,i3,axis  ,0)+an2*rsxy(i1,i2,i3,axis  ,1) )
                       t2=a1*( an1*rsxy(i1,i2,i3,axisp1,0)+an2*rsxy(i1,i2,i3,axisp1,1) ) 
                       ur0 = (ff - ( t2*urv(axisp1) + a0*u(i1,i2,i3,uc) ) )/t1
                     else
                       urv(2) = ut2(i1,i2,i3,uc)
                       t1=a1*( an1*rsxy(i1,i2,i3,axis  ,0)+an2*rsxy(i1,i2,i3,axis  ,1)+an3*rsxy(i1,i2,i3,axis  ,2) )
                       t2=a1*( an1*rsxy(i1,i2,i3,axisp1,0)+an2*rsxy(i1,i2,i3,axisp1,1)+an3*rsxy(i1,i2,i3,axisp1,2) )
                       t3=a1*( an1*rsxy(i1,i2,i3,axisp2,0)+an2*rsxy(i1,i2,i3,axisp2,1)+an3*rsxy(i1,i2,i3,axisp2,2) )
                      ur0 = ( ff - ( t2*urv(axisp1) + t3*urv(axisp2) + a0*u(i1,i2,i3,uc) ) )/t1
                    end if
                    u(j1,j2,j3,uc) =  u(j1+2*is1,j2+2*is2,j3+2*is3,uc) -2.*is*dr(axis)*ur0
                 end if ! curvilinear grid 
                 ! ----- Assign extra ghost ----
                   if( numGhost.gt.1 )then
                     ! extrap second ghost (UPW)
                     u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))          
                   end if
               else if( mask(i1,i2,i3).lt.0 )then
                 ! ----- extrap ghost outside interp. pts on physical boundaries ------
                   u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                   if( numGhost.gt.1 )then
                     u(k1,k2,k3,uc)=(3.*u(j1,j2,j3,uc)-3.*u(j1+is1,j2+is2,j3+is3,uc)+u(j1+2*is1,j2+2*is2,j3+2*is3,uc))
                   end if
               end if
               end do
               end do
               end do
         else if( orderOfAccuracy.eq.4 )then
         else if( orderOfAccuracy.eq.6 )then
         else if( orderOfAccuracy.eq.8 )then   
         else
           write(*,'("CgWave::bcOpt:ERROR: unexpected orderOfAccuracy=",i6)') orderOfAccuracy
           stop 8888
         end if
       end if
     else if( bc(side,axis) == absorbing .or. bc(side,axis) == abcEM2 )then
       ! --- ABC's are done elsewhere ---
       ! write(*,'("bcOptWave: bc=absorbing/EM2 called??")')
       ! stop 3434
     else if(  bc(side,axis).eq.dirichlet .or. bc(side,axis).eq.exactBC .or. bc(side,axis).le.0 )then
       ! do nothing
     else
       write(*,'("bcOptWave: unexpected boundaryCondition=",i4)') bc(side,axis)
       stop 5151
     end if 
    end do ! end side
    end do ! end axis
  if( .false.  ) then
     n1a=gridIndexRange(0,0)
     n1b=gridIndexRange(1,0)
     n2a=gridIndexRange(0,1)
     n2b=gridIndexRange(1,1)
     n3a=gridIndexRange(0,2)
     n3b=gridIndexRange(1,2)
     write(*,'(/,"bcOpt: After dirichlet COMPAT, grid=",i3," n1a,n1b,n2a,n2b,n3a,n3b=",6i5)') grid,n1a,n1b,n2a,n2b,n3a,n3b
     do i3=n3a-numGhost3,n3b+numGhost3
       if( nd.eq.3 )then
         write(*,'("i3=",i4)') i3
       else
         write(*,'("i1=[",i4," : ",i4,"]")') n1a-numGhost,n1b+numGhost
       end if
       do i2=n2a-numGhost,n2b+numGhost
         write(*,'("i2=",i4,1x,100(1pe10.2))') i2,(u(i1,i2,i3,0),i1=n1a-numGhost,n1b+numGhost)
       end do
     end do
   end if
   ! **TEST** DEC 2, 2023
   ! Dec 8, 2023 : this should no longer be needed as we now restore the Dirichlet values
   if( .false. .and. nd==3 .and. bcApproach==useCompatibilityBoundaryConditions )then
     if( t.le.2*dt )then
       write(*,'(/,"*** bcOptWave: RE-ASSIGN DIRICHLET BCS ***",/)') 
     end if
       ! NOTE: the numGhost args are used in ghost loops
       extraForDirichlet=numGhost
       ff =0. ! default value 
        extra1a=extraForDirichlet
        extra1b=extraForDirichlet
        extra2a=extraForDirichlet
        extra2b=extraForDirichlet
        if( nd.eq.3 )then
          extra3a=extraForDirichlet
          extra3b=extraForDirichlet
        else
          extra3a=0
          extra3b=0
        end if
        if( bc(0,0).lt.0 )then
          extra1a=max(0,extra1a) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions (or internal parallel boundaries)
        else if( bc(0,0).eq.0 )then
          extra1a=numGhost  ! include interpolation points since we assign ghost points outside these
        end if
        ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
        if( bc(1,0).lt.0 )then
          extra1b=max(0,extra1b) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions
        else if( bc(1,0).eq.0 )then
          extra1b=numGhost
        end if
        if( bc(0,1).lt.0 )then
          extra2a=max(0,extra2a) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions (or internal parallel boundaries)
        else if( bc(0,1).eq.0 )then
          extra2a=numGhost  ! include interpolation points since we assign ghost points outside these
        end if
        ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
        if( bc(1,1).lt.0 )then
          extra2b=max(0,extra2b) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions
        else if( bc(1,1).eq.0 )then
          extra2b=numGhost
        end if
        if(  nd.eq.3 )then
         if( bc(0,2).lt.0 )then
           extra3a=max(0,extra3a) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions (or internal parallel boundaries)
         else if( bc(0,2).eq.0 )then
           extra3a=numGhost  ! include interpolation points since we assign ghost points outside these
         end if
         ! **NOTE** the bc on the right may be negative even it is not on the left (for parallel)
         if( bc(1,2).lt.0 )then
           extra3b=max(0,extra3b) ! over-ride extraForDirichlet=-1 : assign ends in periodic directions
         else if( bc(1,2).eq.0 )then
           extra3b=numGhost
         end if
        end if
        do axis=0,nd-1
        do side=0,1
          if( .true. .or. bc(side,axis).gt.0 )then ! we may set ghost outside interp for implicit
            ! write(*,'(" bcOpt: side,axis,bc=",3i2)') side,axis,bc(side,axis)
            n1a=gridIndexRange(0,0)
            n1b=gridIndexRange(1,0)
            n2a=gridIndexRange(0,1)
            n2b=gridIndexRange(1,1)
            n3a=gridIndexRange(0,2)
            n3b=gridIndexRange(1,2)
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
            nn1a=gridIndexRange(0,0)-extra1a
            nn1b=gridIndexRange(1,0)+extra1b
            nn2a=gridIndexRange(0,1)-extra2a
            nn2b=gridIndexRange(1,1)+extra2b
            nn3a=gridIndexRange(0,2)-extra3a
            nn3b=gridIndexRange(1,2)+extra3b
            if( axis.eq.0 )then
              nn1a=gridIndexRange(side,axis)
              nn1b=gridIndexRange(side,axis)
            else if( axis.eq.1 )then
              nn2a=gridIndexRange(side,axis)
              nn2b=gridIndexRange(side,axis)
            else
              nn3a=gridIndexRange(side,axis)
              nn3b=gridIndexRange(side,axis)
            end if
            is=1-2*side
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
            axisp1=mod(axis+1,nd)
            axisp2=mod(axis+2,nd)
            i3=n3a
            if( debug.gt.31 )then
              write(*,'(" bcOptWave: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i5)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
              write(*,'(" bcOptWave: extraForDirichlet,extra1a,extra2a,extra3a=",4i4,", loop bounds: nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i5)') extraForDirichlet,extra1a,extra2a,extra3a,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
            end if
          end if ! if bc>0 
          assignTwilightZone=twilightZone
         if( bc(side,axis).eq.exactBC )then
           ! ==== Set the boundary and ghost with the exact solution ====
           ! NOTE: known solutions are now done in applyBoundaryCondtions.bC 
           if( assignTwilightZone.eq.1 )then
               ! assign extram points in the tangential directions
               extram = numGhost 
                 m1a=gridIndexRange(0,0)-extram
                 m1b=gridIndexRange(1,0)+extram
                 m2a=gridIndexRange(0,1)-extram
                 m2b=gridIndexRange(1,1)+extram
                 if( nd.eq.2 )then
                   m3a=gridIndexRange(0,2)
                   m3b=gridIndexRange(1,2)
                 else
                   m3a=gridIndexRange(0,2)-extram
                   m3b=gridIndexRange(1,2)+extram
                 end if
                 if( axis.eq.0 )then
                  m1a=gridIndexRange(side,axis)
                  m1b=gridIndexRange(side,axis)
                 else if( axis.eq.1 )then
                  m2a=gridIndexRange(side,axis)
                  m2b=gridIndexRange(side,axis)
                 else
                  m3a=gridIndexRange(side,axis)
                  m3b=gridIndexRange(side,axis)
                 end if
               ff=0.
               do i3=m3a,m3b
               do i2=m2a,m2b
               do i1=m1a,m1b
                 if( mask(i1,i2,i3).ne.0 )then
                   do ghost=0,numGhost
                     j1=i1-is1*ghost
                     j2=i2-is2*ghost
                     j3=i3-is3*ghost
                       if( assignTwilightZone.eq.1 )then
                         ! compute RHS from TZ
                         if( nd.eq.2 )then
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue )
                         else
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),xy(j1,j2,j3,2),t,uc,ue )
                         end if
                         ff = ue
                       else if( assignKnownSolutionAtBoundaries.eq.1 )then
                         ! -- we set inhomogeneous Dirichlet values for some known solutions 
                         if( knownSolutionOption.eq.planeWave )then
                           ! --- evaluate the plane wave solution ---
                           if( nd.eq.2 )then
                             ff = ampPlaneWave*sin( kxPlaneWave*xy(j1,j2,j3,0) + kyPlaneWave*xy(j1,j2,j3,1) - omegaPlaneWave*t )
                           else
                             ff = ampPlaneWave*sin( kxPlaneWave*xy(j1,j2,j3,0) + kyPlaneWave*xy(j1,j2,j3,1) + kzPlaneWave*xy(j1,j2,j3,2) - omegaPlaneWave*t )
                           end if 
                         else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                           ! Eval the Gaussian plane wave solution
                           !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                           !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                           !  
                           if( nd.eq.2 )then
                             xi = kxGPW*(xy(j1,j2,j3,0)-x0GPW) + kyGPW*(xy(j1,j2,j3,1)-y0GPW) - c*t
                           else
                             xi = kxGPW*(xy(j1,j2,j3,0)-x0GPW) + kyGPW*(xy(j1,j2,j3,1)-y0GPW) + kzGPW*(xy(j1,j2,j3,2)-z0GPW) - c*t
                           end if 
                           ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                         else if( knownSolutionOption.eq.boxHelmholtz ) then
                           ! --- evaluate the boxHelmholtz solution ---
                           ! For multi-freq we add up all the component frequencies
                           ff = 0. 
                           do freq=0,numberOfFrequencies-1
                             ! kx = kxBoxHelmholtz + twoPi*freq
                             ! ky = kyBoxHelmholtz + twoPi*freq
                             ! kz = kzBoxHelmholtz + twoPi*freq
                             ! coswt = cos( frequencyArray(freq)*t )
                                 ! This macro is used in bcOptWave.bf90
                                 omega = frequencyArray(freq);
                                 kx = kxBoxHelmholtz*(freq*.5+1.)
                                 ky = kyBoxHelmholtz*(freq*.5+1.)
                                 kz = kzBoxHelmholtz*(freq*.5+1.)
                                 ! kx = kxBoxHelmholtz + twoPi*freq.
                                 ! ky = kyBoxHelmholtz + twoPi*freq.
                                 ! kz = kzBoxHelmholtz + twoPi*freq.    
                             coswt = cos( omega*t )
                             if( nd.eq.2 )then
                               ff = ff + sin( kx*xy(j1,j2,j3,0) ) * sin( ky*xy(j1,j2,j3,1) ) *coswt
                             else
                               ff = ff + sin( kx*xy(j1,j2,j3,0) ) * sin( ky*xy(j1,j2,j3,1) ) * sin( kz*xy(j1,j2,j3,2) ) *coswt
                             end if 
                           end do
                         else if( knownSolutionOption.eq.polyPeriodic ) then
                           ! --- evaluate the polyPeriodic solution ---
                           if( nd.eq.2 )then
                             ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(j1,j2,j3,0) + b1PolyPeriodic*xy(j1,j2,j3,1)                                 ) *coswt
                           else
                             ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(j1,j2,j3,0) + b1PolyPeriodic*xy(j1,j2,j3,1) + c1PolyPeriodic*xy(j1,j2,j3,2) ) *coswt
                           end if 
                         else
                           stop 9876
                         end if 
                       end if
                     u(j1,j2,j3,uc) = ff
                   end do
                 end if ! mask .ne. 0
                end do
                end do
                end do
           end if
         else if( bc(side,axis).eq.dirichlet )then
           if( bcApproach.eq.useCompatibilityBoundaryConditions )then
             ! --- Assign values on the boundary for CBCs ---
             if( orderOfAccuracy.eq.2 )then
               if( gridType.eq.rectangular )then
                   ! if( t.le.2*dt )then
                   !   write(*,'("assign extended Dirichlet BC: side,axis=",2i2," nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                   ! end if
                   ff=0.
                   if( addForcingBC.eq.1 )then
                     ! beginLoops3d()
                     ! *wdh* Dec 1, 2023: assign extended bndry
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                       if( mask(i1,i2,i3).ne.0 )then
                         ! --- get the RHS to the Dirichlet BC ---
                         ! #If #FORCING eq "USEFORCING" *wdh* This was a bug, turned off Nov 22, 2023
                           if( assignTwilightZone.eq.1 )then
                             ! compute RHS from TZ
                             if( nd.eq.2 )then
                               call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                             else
                               call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                             end if
                             ff = ue
                           else if( assignKnownSolutionAtBoundaries.eq.1 )then
                             ! -- we set inhomogeneous Dirichlet values for some known solutions 
                             if( knownSolutionOption.eq.planeWave )then
                               ! --- evaluate the plane wave solution ---
                               if( nd.eq.2 )then
                                 ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                               else
                                 ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                               end if 
                             else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                               ! Eval the Gaussian plane wave solution
                               !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                               !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                               !  
                               if( nd.eq.2 )then
                                 xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) - c*t
                               else
                                 xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) + kzGPW*(xy(i1,i2,i3,2)-z0GPW) - c*t
                               end if 
                               ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                             else if( knownSolutionOption.eq.boxHelmholtz ) then
                               ! --- evaluate the boxHelmholtz solution ---
                               ! For multi-freq we add up all the component frequencies
                               ff = 0. 
                               do freq=0,numberOfFrequencies-1
                                 ! kx = kxBoxHelmholtz + twoPi*freq
                                 ! ky = kyBoxHelmholtz + twoPi*freq
                                 ! kz = kzBoxHelmholtz + twoPi*freq
                                 ! coswt = cos( frequencyArray(freq)*t )
                                     ! This macro is used in bcOptWave.bf90
                                     omega = frequencyArray(freq);
                                     kx = kxBoxHelmholtz*(freq*.5+1.)
                                     ky = kyBoxHelmholtz*(freq*.5+1.)
                                     kz = kzBoxHelmholtz*(freq*.5+1.)
                                     ! kx = kxBoxHelmholtz + twoPi*freq.
                                     ! ky = kyBoxHelmholtz + twoPi*freq.
                                     ! kz = kzBoxHelmholtz + twoPi*freq.    
                                 coswt = cos( omega*t )
                                 if( nd.eq.2 )then
                                   ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) *coswt
                                 else
                                   ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) * sin( kz*xy(i1,i2,i3,2) ) *coswt
                                 end if 
                               end do
                             else if( knownSolutionOption.eq.polyPeriodic ) then
                               ! --- evaluate the polyPeriodic solution ---
                               if( nd.eq.2 )then
                                 ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1)                                 ) *coswt
                               else
                                 ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                               end if 
                             else
                               stop 9876
                             end if 
                           end if
                         ! #End
                         ! --- Dirichlet BC ---
                         u(i1,i2,i3,uc)=ff
                         ! write(*,'("assignDirichletBndry: side,axis=",2i2," i1,i2=",2i3," ff=",e10.2)') side,axis,i1,i2,ff
                       end if ! mask .ne. 0
                      end do
                      end do
                      end do
                   else
                     ! Homogenous BCs 
                     ! beginLoops3d()
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                       u(i1,i2,i3,uc)=0.
                      end do
                      end do
                      end do
                   end if
               else
                   ! if( t.le.2*dt )then
                   !   write(*,'("assign extended Dirichlet BC: side,axis=",2i2," nn1a,nn1b,nn2a,nn2b,nn3a,nn3b=",6i3)') side,axis,nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
                   ! end if
                   ff=0.
                   if( addForcingBC.eq.1 )then
                     ! beginLoops3d()
                     ! *wdh* Dec 1, 2023: assign extended bndry
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                       if( mask(i1,i2,i3).ne.0 )then
                         ! --- get the RHS to the Dirichlet BC ---
                         ! #If #FORCING eq "USEFORCING" *wdh* This was a bug, turned off Nov 22, 2023
                           if( assignTwilightZone.eq.1 )then
                             ! compute RHS from TZ
                             if( nd.eq.2 )then
                               call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                             else
                               call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                             end if
                             ff = ue
                           else if( assignKnownSolutionAtBoundaries.eq.1 )then
                             ! -- we set inhomogeneous Dirichlet values for some known solutions 
                             if( knownSolutionOption.eq.planeWave )then
                               ! --- evaluate the plane wave solution ---
                               if( nd.eq.2 )then
                                 ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                               else
                                 ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                               end if 
                             else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                               ! Eval the Gaussian plane wave solution
                               !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                               !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                               !  
                               if( nd.eq.2 )then
                                 xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) - c*t
                               else
                                 xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) + kzGPW*(xy(i1,i2,i3,2)-z0GPW) - c*t
                               end if 
                               ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                             else if( knownSolutionOption.eq.boxHelmholtz ) then
                               ! --- evaluate the boxHelmholtz solution ---
                               ! For multi-freq we add up all the component frequencies
                               ff = 0. 
                               do freq=0,numberOfFrequencies-1
                                 ! kx = kxBoxHelmholtz + twoPi*freq
                                 ! ky = kyBoxHelmholtz + twoPi*freq
                                 ! kz = kzBoxHelmholtz + twoPi*freq
                                 ! coswt = cos( frequencyArray(freq)*t )
                                     ! This macro is used in bcOptWave.bf90
                                     omega = frequencyArray(freq);
                                     kx = kxBoxHelmholtz*(freq*.5+1.)
                                     ky = kyBoxHelmholtz*(freq*.5+1.)
                                     kz = kzBoxHelmholtz*(freq*.5+1.)
                                     ! kx = kxBoxHelmholtz + twoPi*freq.
                                     ! ky = kyBoxHelmholtz + twoPi*freq.
                                     ! kz = kzBoxHelmholtz + twoPi*freq.    
                                 coswt = cos( omega*t )
                                 if( nd.eq.2 )then
                                   ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) *coswt
                                 else
                                   ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) * sin( kz*xy(i1,i2,i3,2) ) *coswt
                                 end if 
                               end do
                             else if( knownSolutionOption.eq.polyPeriodic ) then
                               ! --- evaluate the polyPeriodic solution ---
                               if( nd.eq.2 )then
                                 ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1)                                 ) *coswt
                               else
                                 ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                               end if 
                             else
                               stop 9876
                             end if 
                           end if
                         ! #End
                         ! --- Dirichlet BC ---
                         u(i1,i2,i3,uc)=ff
                         ! write(*,'("assignDirichletBndry: side,axis=",2i2," i1,i2=",2i3," ff=",e10.2)') side,axis,i1,i2,ff
                       end if ! mask .ne. 0
                      end do
                      end do
                      end do
                   else
                     ! Homogenous BCs 
                     ! beginLoops3d()
                     do i3=nn3a,nn3b
                     do i2=nn2a,nn2b
                     do i1=nn1a,nn1b
                       u(i1,i2,i3,uc)=0.
                      end do
                      end do
                      end do
                   end if
               end if 
             else if( orderOfAccuracy.eq.4 )then
               if( t.le.3*dt .and. debug.gt.1 )then
                 write(*,'("APPLY BC order =4 useCompatibilityBoundaryConditions, side,axis=",2i4)') side,axis
               end if
                 stop 444 
             else if( orderOfAccuracy.eq.6 )then 
               write(*,'("CgWave::bcOpt:ERROR: CBC order=6 turned off, use LCBC")') 
               stop 666
             else if( orderOfAccuracy.eq.8 )then   
               stop 888
             else
               write(*,'("CgWave::bcOpt:ERROR: unexpected orderOfAccuracy=",i6)') orderOfAccuracy
               stop 8888
             end if
           else 
             ! ----- assign boundary and then ghost by extrapolation ----
             ! This was a temporary method until other approaches were implemented
             if( orderOfAccuracy.eq.2 )then
                   ff=0.
                   if( addForcingBC.eq.1 )then  
                      do i3=n3a,n3b
                      do i2=n2a,n2b
                      do i1=n1a,n1b
                       if( mask(i1,i2,i3).ne.0 )then
                         ! --- get the RHS to the Dirichlet BC ---
                           if( assignTwilightZone.eq.1 )then
                             ! compute RHS from TZ
                             if( nd.eq.2 )then
                               call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
                             else
                               call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
                             end if
                             ff = ue
                           else if( assignKnownSolutionAtBoundaries.eq.1 )then
                             ! -- we set inhomogeneous Dirichlet values for some known solutions 
                             if( knownSolutionOption.eq.planeWave )then
                               ! --- evaluate the plane wave solution ---
                               if( nd.eq.2 )then
                                 ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
                               else
                                 ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
                               end if 
                             else if( knownSolutionOption.eq.gaussianPlaneWave ) then
                               ! Eval the Gaussian plane wave solution
                               !    u = exp( -beta*(xi^2) )*cos( k0*xi )
                               !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
                               !  
                               if( nd.eq.2 )then
                                 xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) - c*t
                               else
                                 xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) + kzGPW*(xy(i1,i2,i3,2)-z0GPW) - c*t
                               end if 
                               ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      
                             else if( knownSolutionOption.eq.boxHelmholtz ) then
                               ! --- evaluate the boxHelmholtz solution ---
                               ! For multi-freq we add up all the component frequencies
                               ff = 0. 
                               do freq=0,numberOfFrequencies-1
                                 ! kx = kxBoxHelmholtz + twoPi*freq
                                 ! ky = kyBoxHelmholtz + twoPi*freq
                                 ! kz = kzBoxHelmholtz + twoPi*freq
                                 ! coswt = cos( frequencyArray(freq)*t )
                                     ! This macro is used in bcOptWave.bf90
                                     omega = frequencyArray(freq);
                                     kx = kxBoxHelmholtz*(freq*.5+1.)
                                     ky = kyBoxHelmholtz*(freq*.5+1.)
                                     kz = kzBoxHelmholtz*(freq*.5+1.)
                                     ! kx = kxBoxHelmholtz + twoPi*freq.
                                     ! ky = kyBoxHelmholtz + twoPi*freq.
                                     ! kz = kzBoxHelmholtz + twoPi*freq.    
                                 coswt = cos( omega*t )
                                 if( nd.eq.2 )then
                                   ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) *coswt
                                 else
                                   ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) * sin( kz*xy(i1,i2,i3,2) ) *coswt
                                 end if 
                               end do
                             else if( knownSolutionOption.eq.polyPeriodic ) then
                               ! --- evaluate the polyPeriodic solution ---
                               if( nd.eq.2 )then
                                 ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1)                                 ) *coswt
                               else
                                 ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
                               end if 
                             else
                               stop 9876
                             end if 
                           end if
                         ! --- Dirichlet BC ---
                         u(i1,i2,i3,uc)=ff
                         ! -- extrapolate ghost ---
                         do ghost=1,numGhost
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc) = (3.*u(j1+is1,j2+is2,j3+is3,uc)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)) 
                         end do
                       end if ! mask .ne. 0
                      end do
                      end do
                      end do
                   else
                     ! --- no forcing ----
                     if( numGhost.eq.1 )then
                       ghost=1
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ! --- Dirichlet BC ---
                         u(i1,i2,i3,uc)=0.
                         ! -- extrapolate ghost ---
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost         
                         u(j1,j2,j3,uc) = (3.*u(j1+is1,j2+is2,j3+is3,uc)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)) 
                        end do
                        end do
                        end do
                     else 
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ! --- Dirichlet BC ---
                         u(i1,i2,i3,uc)=0.
                         ! -- extrapolate ghost ---
                         do ghost=1,numGhost
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc) = (3.*u(j1+is1,j2+is2,j3+is3,uc)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)) 
                         end do
                        end do
                        end do
                        end do
                     end if
                     if( .false. )then
                       ! -- two loops ---
                       ! --- no forcing ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ! --- Dirichlet BC ---
                         u(i1,i2,i3,uc)=0.
                        end do
                        end do
                        end do
                       ! -- extrapolate ghost ---
                       do ghost=1,numGhost
                          do i3=n3a,n3b
                          do i2=n2a,n2b
                          do i1=n1a,n1b
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost               
                           u(j1,j2,j3,uc) = (3.*u(j1+is1,j2+is2,j3+is3,uc)-3.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)) 
                          end do
                          end do
                          end do
                       end do
                     end if
                   endif
             else if( orderOfAccuracy.eq.4 )then
             else if( orderOfAccuracy.eq.6 )then   
             else if( orderOfAccuracy.eq.8 )then   
               stop 888
             else
               write(*,'("CgWave::bcOpt:ERROR: unexpected orderOfAccuracy=",i6)') orderOfAccuracy
               stop 8888
             end if
           end if
         end if ! end if dirichlet 
        end do ! end side
        end do ! end axis
   end if
   if( assignCornerGhostPoints.eq.1 )then
     ! call new routine: nov 28, 2024 -- 
     ! assignCornerGhostPoints : this variable is set in the calling routine
     call cornersWave( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange, dimRange, isPeriodic, u, un, mask,rsxy, xy, uTemp, v, boundaryCondition, frequencyArray, pdb, ipar, rpar, ierr )
   else
     !  --- Assign ghost points outside corners ---
     if( .false. .and. orderOfAccuracy==4 .and. gridType==rectangular .and. bcApproach==useCompatibilityBoundaryConditions )then 
       ! *wdh* Nov 28, 2023
       if( t.le.3*dt .and. debug.gt.1 )then
         write(*,'("bcOpt: Assign symmetry corner ghost for CBC")')
         stop 1234
       end if
       ! if( forcingOption.eq.noForcing )then
       !   assignSymmetryCornerGhost(noForcing)
       ! else
       !   assignSymmetryCornerGhost(forcing)
       ! end if
     else
       ! This is broken, at least in 3D: May 2, 2023 
       ! 
       if( nd.eq.2 )then
         ! Turn this back on for 2D -- Nov 23, 2023
         if( t.le.3*dt .and. debug.gt.1 )then
           write(*,'("bcOpt: Assign corner ghost")')
         end if
           ! ---------------------------------
           ! --- assign corners and edges: ---
           ! ---------------------------------
          ! From op/src/BoundaryConditionParameters.h
          ! enum CornerBoundaryConditionEnum
          !  {
          !    doNothingCorner=-1,  
          !    extrapolateCorner=0,
          !    symmetryCorner,  // should be replaced by the one of the odd,even below -- keep for compatibility
          !    taylor2ndOrder,  // should be replaced by the taylor2ndOrderOddCorner below -- keep for compatibility
          !    evenSymmetryCorner,
          !    oddSymmetryCorner,
          !    taylor2ndOrderEvenCorner,
          !    taylor4thOrderEvenCorner,
          !    vectorSymmetryAxis1Corner,       // even symmetry on all variables except normal component of the "velocity"
          !    vectorSymmetryAxis2Corner, 
          !    vectorSymmetryAxis3Corner
          !  };
           if( nd.eq.3 )then
             write(*,'("assignCornerGhost -- *check me for 3D : stop here for now")')
             stop 666
           end if
           do side3=0,2
           do side2=0,2
           do side1=0,2
             cornerBC(side1,side2,side3)=0
           end do
           end do
           end do
           ! ** FIX ME :
           !   cornerBC(0:2,0:2,0:2) : 2=edge in 3D 
           do side3=0,1
           do side2=0,1
           do side1=0,1
             if( orderOfAccuracy.ge.4. .and. bc(side1,0).eq.neumann .and. bc(side2,1).eq.neumann .and. ( nd.eq.2 .or. bc(side3,2).eq.neumann ) )then
               ! ---- This is a Neumann-Neumann corner ----
               cornerBC(side1,side2,side3)=extrapolateCorner
               ! if( t.le.2*dt  )then
               !   write(*,'("Assign special Neumann corners conditions ")')
               ! end if
               ! cornerBC(side1,side2,side3)=taylor4thOrderEvenCorner
               ! cornerBC(side1,side2,side3)=evenSymmetryCorner
               ! cornerBC(side1,side2,side3)=0 
             else if(            bc(side1,0).eq.exactBC .and. bc(side2,1).eq.exactBC .and. ( nd.eq.2 .or. bc(side3,2).eq.exactBC ) )then
               ! ---- Do nothing at this exact corner 
               cornerBC(side1,side2,side3)=-1
             else 
               cornerBC(side1,side2,side3)=0         ! extrapolateCorner=0, (BoundaryConditionParameters)
             end if 
           end do
           end do
           end do
           ! orderOfExtrapolationForCorners=5
           ! orderOfExtrapolationForCorners= orderOfAccuracy+1 ! *wdh* July 21, 2024
           orderOfExtrapolationForCorners= orderOfExtrapolation;
           iparc(0)=uc
           iparc(1)=uc
           iparc(2)=0                              ! useWhereMask;
           iparc(3)=orderOfExtrapolationForCorners
           iparc(4)=numGhost                       ! numberOfCornerGhostLinesToAssign
           iparc(5)=0                              ! cornerExtrapolationOption : 0=extrap along diagonals
           iparc(6)=0                              ! vectorSymmetryCornerComponent
           iparc(7)=gridType
           rparc(0)=epsx ! normEps
           ! Note: is it ok to use gridIndexRange instead of indexRange here: ??
           call fixBoundaryCornersOpt( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,0,uc,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, u,mask,rsxy, gridIndexRange, dimRange, isPeriodic, boundaryCondition, cornerBC, iparc, rparc )
       end if
     end if
   end if
   return
   end
