! This file automatically generated from bcOptWave.bf90 with bpp.
 subroutine bcOptWave2dOrder4( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
 ! subroutine bcOptWave( nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,!                       gridIndexRange, dimRange, isPeriodic, u, mask,rsxy, xy, boundaryCondition, !                       frequencyArray, pdb, ipar, rpar, ierr )
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
   integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
   real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
   real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
   integer gridIndexRange(0:1,0:2),boundaryCondition(0:1,0:2), dimRange(0:1,0:2), isPeriodic(0:*)
   real frequencyArray(0:*)
   ! temp space for CBC order 4 -- fix me : just make a stencil
   real uTemp(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
   real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
   ! ! *** TEMP ARRAYS FOR WORK SPACE --> THIS IS SLOW!!!
   ! #If "4" eq "4"
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
   integer twilightZone,numberOfProcessors,addForcingBC,assignBCForImplicitForImplicit
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
   integer l1,l2,l3
   integer numGhost,numberOfGhostPoints,extraForNeumann,extraForDirichlet,numberOfFrequencies
   integer side1,side2,side3
   integer n1a,n1b,n2a,n2b,n3a,n3b
   integer m1a,m1b,m2a,m2b,m3a,m3b
   integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
   integer j1a,j1b,j2a,j2b,ja,j3b
   integer extra1a,extra1b,extra2a,extra2b,extra3a,extra3b,extram
   integer cornerBC(0:2,0:2,0:2), iparc(0:10), orderOfExtrapolationForCorners
   real rparc(0:10)
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
   real ampPlaneWave, kxPlaneWave,kyPlaneWave,kzPlaneWave, omegaPlaneWave
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
   logical firstTimeForCBC6
   real symSign
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
   checkCoeff=0 ! set to 1 to check coefficients in CBCs using discrete delta approach
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
   assignBCForImplicitForImplicit  = ipar(16)
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
   c2 = c**2
   c4 = c**4
   c6 = c**6
   c8 = c**8
   twoPi = atan2(1.,1.)*8.; ! atan2(1,1)=pi/4
   ! numberOfGhostPoints=orderOfAccuracy/2
   numberOfGhostPoints=numGhost ! now passed in 
   ! write(*,'(" bcOptWave2dOrder4: dim=2, order=4")')
   if( t.le.dt .and. debug.gt.1 )then
   ! if( .true. )then
     write(*,'(" bcOptWave: grid=",i4," gridType=",i2," orderOfAccuracy=",i2," uc=",i3," twilightZone=",i2)') grid,gridType,orderOfAccuracy,uc,twilightZone
     write(*,'("  addForcingBC=",i4," forcingOption=",i4," assignKnownSolutionAtBoundaries=",i4)') addForcingBC, forcingOption, assignKnownSolutionAtBoundaries
     write(*,'("  t=",e10.2," dt=",e10.2," knownSolutionOption=",i4," REAL_MIN=",e10.2)') t,dt,knownSolutionOption,REAL_MIN
     write(*,'("  useUpwindDissipation=",i2," numGhost=",i2)') useUpwindDissipation,numGhost
     write(*,'("  assignBCForImplicitForImplicit=",i4," bcApproach=",i4)') assignBCForImplicitForImplicit,bcApproach
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
       if( bc(side,axis).eq.absorbing )then
         ! DO THIS FOR NOW : use dirichlet for absorbing boundaries 
         bc(side,axis)=dirichlet
       end if
     end do
   end do
   if( .true. ) then
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
       if(  t.le.dt .and. debug.gt.1  )then
         write(*,'(" bcOptWave:  knownSolutionOption=planeWave: ampPlaneWave=",e10.2," kxPlaneWave=",e10.2," kyPlaneWave=",e10.2)') ampPlaneWave,kxPlaneWave,kyPlaneWave
       end if 
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
   if( assignBCForImplicitForImplicit.eq.1 )then
     ! -------- IMPLICIT BoundaryConditions --------
     ! if( .true. )then
     !   write(*,'("bcOptWave: fill BCs into RHS for implicit solver")')
     ! end if
     ! write(*,'("FINISH ME")')
     ! stop 6789
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
        if( bc(side,axis).gt.0 )then
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
          if( debug.gt.7 )then
            write(*,'(" bcOptWave: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
          end if
        end if ! if bc>0 
        assignTwilightZone=twilightZone
       if( bc(side,axis) == dirichlet )then
         ff=0.
         if( bcApproach.eq.useOneSidedBoundaryConditions )then
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
           write(*,'("bcOpt: implicit: finish me for CBC")') 
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
         else
           write(*,'("bcOpt: implicit: unexpected bcApproach=",i6)') bcApproach
           stop 1111
         end if
       else if( bc(side,axis) == exactBC )then
         write(*,*) "bcOpt: implicit BC for exactBc -- finish me"
         stop 4444
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
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
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
           end if 
          end do
          end do
          end do
       else if( bc(side,axis) > 0 )then
         write(*,'("bcOptWave:fill RHS for direct Helmholtz solver, unexpected boundaryCondition=",i4)') bc(side,axis)
         stop 6666
       end if
      end do ! end side
      end do ! end axis
     ! ---------------- RETURN ---------------
     return
   end if
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
      if( bc(side,axis).gt.0 )then
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
        if( debug.gt.7 )then
          write(*,'(" bcOptWave: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
        end if
      end if ! if bc>0 
      assignTwilightZone=twilightZone
     if( bc(side,axis).eq.exactBC )then
       ! ==== Set the boundary and ghost with the exact solution ====
       ! *wdh* June 13, 2023 
       ! ***** this case is now done in applyBoundaryCondtions.bC *******
       !! assignExactBoundary()
     else if( bc(side,axis).eq.dirichlet )then
       if( bcApproach.eq.useCompatibilityBoundaryConditions )then
         ! --- Assign values on the boundary for CBCs ---
         if( orderOfAccuracy.eq.2 )then
             stop 222 
         else if( orderOfAccuracy.eq.4 )then
           if( gridType.eq.rectangular )then
               ff=0.
               if( addForcingBC.eq.1 )then
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                   if( mask(i1,i2,i3).ne.0 )then
                     ! --- get the RHS to the Dirichlet BC ---
                     ! --- Dirichlet BC ---
                     u(i1,i2,i3,uc)=ff
                   end if ! mask .ne. 0
                  end do
                  end do
                  end do
               else
                 ! Homogenous BCs 
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                   u(i1,i2,i3,uc)=0.
                  end do
                  end do
                  end do
               end if
           else
               ff=0.
               if( addForcingBC.eq.1 )then
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                   if( mask(i1,i2,i3).ne.0 )then
                     ! --- get the RHS to the Dirichlet BC ---
                     ! --- Dirichlet BC ---
                     u(i1,i2,i3,uc)=ff
                   end if ! mask .ne. 0
                  end do
                  end do
                  end do
               else
                 ! Homogenous BCs 
                  do i3=n3a,n3b
                  do i2=n2a,n2b
                  do i1=n1a,n1b
                   u(i1,i2,i3,uc)=0.
                  end do
                  end do
                  end do
               end if
           end if  
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
         else if( orderOfAccuracy.eq.4 )then
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
                     u(j1,j2,j3,uc) = (5.*u(j1+is1,j2+is2,j3+is3,uc)-10.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+10.*u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)-5.*u(j1+is1+3*is1,j2+is2+3*is2,j3+is3+3*is3,uc)+u(j1+is1+4*is1,j2+is2+4*is2,j3+is3+4*is3,uc)) 
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
                   u(j1,j2,j3,uc) = (5.*u(j1+is1,j2+is2,j3+is3,uc)-10.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+10.*u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)-5.*u(j1+is1+3*is1,j2+is2+3*is2,j3+is3+3*is3,uc)+u(j1+is1+4*is1,j2+is2+4*is2,j3+is3+4*is3,uc)) 
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
                     u(j1,j2,j3,uc) = (5.*u(j1+is1,j2+is2,j3+is3,uc)-10.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+10.*u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)-5.*u(j1+is1+3*is1,j2+is2+3*is2,j3+is3+3*is3,uc)+u(j1+is1+4*is1,j2+is2+4*is2,j3+is3+4*is3,uc)) 
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
                     u(j1,j2,j3,uc) = (5.*u(j1+is1,j2+is2,j3+is3,uc)-10.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+10.*u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)-5.*u(j1+is1+3*is1,j2+is2+3*is2,j3+is3+3*is3,uc)+u(j1+is1+4*is1,j2+is2+4*is2,j3+is3+4*is3,uc)) 
                    end do
                    end do
                    end do
                 end do
               end if
             endif
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
      if( bc(side,axis).gt.0 )then
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
        if( debug.gt.7 )then
          write(*,'(" bcOptWave: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
        end if
      end if ! if bc>0 
      assignTwilightZone=twilightZone
     if( ( (bc(side,axis).ne.dirichlet .and. bc(side,axis).ne.exactBC ) .or. bcApproach.eq.useCompatibilityBoundaryConditions ) .and. bc(side,axis).gt.0 )then
       if( orderOfAccuracy.eq.2 )then
       else if( orderOfAccuracy.eq.4 )then
           ! ghost-loops will assign extra points in tangential directions
            do i3=nn3a,nn3b
            do i2=nn2a,nn2b
            do i1=nn1a,nn1b
             if( mask(i1,i2,i3).ne.0 )then
               ! -- extrapolate ghost ---
               do ghost=1,numGhost
                 j1=i1-is1*ghost
                 j2=i2-is2*ghost
                 j3=i3-is3*ghost
                 u(j1,j2,j3,uc) = (5.*u(j1+is1,j2+is2,j3+is3,uc)-10.*u(j1+is1+is1,j2+is2+is2,j3+is3+is3,uc)+10.*u(j1+is1+2*is1,j2+is2+2*is2,j3+is3+2*is3,uc)-5.*u(j1+is1+3*is1,j2+is2+3*is2,j3+is3+3*is3,uc)+u(j1+is1+4*is1,j2+is2+4*is2,j3+is3+4*is3,uc)) 
               end do
             end if ! mask .ne. 0
            end do
            end do
            end do
       else if( orderOfAccuracy.eq.6 )then   
       else if( orderOfAccuracy.eq.8 )then   
       else
         write(*,'("CgWave::bcOpt:ERROR: unexpected orderOfAccuracy=",i6)') orderOfAccuracy
         stop 8888
       end if
     end if
    end do ! end side
    end do ! end axis
   !  --- Assign ghost points outside corners ---
   ! assignCornerGhostsMacro()
   ! Try this for order==6
   ! if( .false. .and. orderOfAccuracy.eq.6 .and. nd.eq.2 .and. gridType.eq.rectangular )then
   !   if( forcingOption.eq.noForcing )then
   !     assignSymmetryCornerGhost(noForcing)
   !   else
   !     assignSymmetryCornerGhost(forcing)
   !   end if
   ! end if
   ! ---------------------------------------------------------------------
   ! ----------- STAGE II : Neumann-like Boundary Conditions -------------
   ! -----------            Ghost values for CBCs            -------------
   ! ---------------------------------------------------------------------
   ! CHECK ME --> numGhost here ??
   extraForNeumann=0 ! only assign Neumann conditions to the boundary
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
      if( bc(side,axis).gt.0 )then
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
        if( debug.gt.7 )then
          write(*,'(" bcOptWave: grid,side,axis=",3i3,", loop bounds: n1a,n1b,n2a,n2b,n3a,n3b=",6i3)') grid,side,axis,n1a,n1b,n2a,n2b,n3a,n3b
        end if
      end if ! if bc>0 
      assignTwilightZone=twilightZone
     if( bc(side,axis).eq.dirichlet .and. bcApproach.eq.useCompatibilityBoundaryConditions )then
       ! -- fill ghost using CBCs ----
       if( orderOfAccuracy.eq.2 )then
       else if( orderOfAccuracy.eq.4 )then
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
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ! --- get the compatibility forcings for order=2 ---
                             ! No forcing, do nothing 
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22r or ulaplacian23r 
                           uLap = ulaplacian22r(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2            ! residual in equation using current ghost value
                           a11 = 1./( dx(axis)**2 )                                  ! coeff of u(-1) in r1 
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                         if( .false. )then
                           ! check the error
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ghost = 1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                       end if ! mask 
                     end do
                     end do
                     end do
                       !----------------------------------------------------------------
                       ! --- STAGE II fill in first ghost by 4th-order compatibility ---
                       !----------------------------------------------------------------
                       ! --- fill in two ghost using 4th- order compatibility ---
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ! --- get the compatibility forcings at order=4 ---
                               ! No forcing, do nothing 
                           ! u_tt = c^2*Lap(u) + f 
                           ! u_tttt = c^2*Lap(u_tt) + f_tt 
                           !        = c^2*Lap( c^2*Lap(u) + f ) + f_tt 
                           !        = c^4*Lap^2(u) + c^2*Lap(f) ) + f_tt 
                             ! uxx43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*dx42(0) 
                             uLap = ulaplacian42r(i1,i2,i3,0)
                             r1 =  c2*uLap + ff - gtt            ! residual in equation using current ghost value
                             a11 = c2*16./(12.*dx(axis)**2)                           ! coeff of u(-1) in r1 
                             a12 =    -c2/(12.*dx(axis)**2)                           ! coeff of u(-2) in r1
                             ! uxxxx = 1 -4 6 -4 1
                             ! vLap = Lap^2( u )
                             vLap = lap2d2Pow2(i1,i2,i3,0)
                             r2 =  c4*vLap + c2*fLap + ftt - gtttt
                             a21 =  -c4*4./(dx(axis)**4)                              ! coeff of u(-1) in r2
                             a22 =   c4*1./(dx(axis)**4)
                         if( .false. )then
                           write(*,'("(i1,i2)=(",2i3,"),  r1,r2,a11,a12,a21,a22,vLap=",7(1pe9.2,1x))') i1,i2,r1,r2,a11,a12,a21,a22,vLap
                         end if
                           ghost =1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost =2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost        
                           ! Solve
                           !  [ a11 a12 ][ u(-1)] = [ a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1 ]
                           !  [ a21 a22 ][ u(-2)]   [ a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2 ]
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2)/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1)/det
                           if( .false. )then
                             ! check the error
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                             write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2," u(-1)=",e10.3," err=",e8.2)') i1,i2,uTemp(j1,j2,j3,0),abs(uTemp(j1,j2,j3,0)-ue1),uTemp(k1,k2,k3,0),abs(uTemp(k1,k2,k3,0)-ue2)
                           end if
                         end if ! mask .gt. 0
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
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
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ! --- get the compatibility forcings for order=2 ---
                             ! No forcing, do nothing 
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22r or ulaplacian23r 
                           uLap = ulaplacian23r(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2            ! residual in equation using current ghost value
                           a11 = 1./( dx(axis)**2 )                                  ! coeff of u(-1) in r1 
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                         if( .false. )then
                           ! check the error
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ghost = 1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                       end if ! mask 
                     end do
                     end do
                     end do
                       !----------------------------------------------------------------
                       ! --- STAGE II fill in first ghost by 4th-order compatibility ---
                       !----------------------------------------------------------------
                       ! --- fill in two ghost using 4th- order compatibility ---
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ! --- get the compatibility forcings at order=4 ---
                               ! No forcing, do nothing 
                           ! u_tt = c^2*Lap(u) + f 
                           ! u_tttt = c^2*Lap(u_tt) + f_tt 
                           !        = c^2*Lap( c^2*Lap(u) + f ) + f_tt 
                           !        = c^4*Lap^2(u) + c^2*Lap(f) ) + f_tt 
                             ! uxx43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*dx42(0) 
                             uLap = ulaplacian43r(i1,i2,i3,0)
                             r1 =  c2*uLap + ff - gtt            ! residual in equation using current ghost value
                             a11 = c2*16./(12.*dx(axis)**2)                           ! coeff of u(-1) in r1 
                             a12 =    -c2/(12.*dx(axis)**2)                           ! coeff of u(-2) in r1
                             ! uxxxx = 1 -4 6 -4 1
                             ! vLap = Lap^2( u )
                             vLap = lap3d2Pow2(i1,i2,i3,0)
                             r2 =  c4*vLap + c2*fLap + ftt - gtttt
                             a21 =  -c4*4./(dx(axis)**4)                              ! coeff of u(-1) in r2
                             a22 =   c4*1./(dx(axis)**4)
                         if( .false. )then
                           write(*,'("(i1,i2)=(",2i3,"),  r1,r2,a11,a12,a21,a22,vLap=",7(1pe9.2,1x))') i1,i2,r1,r2,a11,a12,a21,a22,vLap
                         end if
                           ghost =1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost =2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost        
                           ! Solve
                           !  [ a11 a12 ][ u(-1)] = [ a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1 ]
                           !  [ a21 a22 ][ u(-2)]   [ a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2 ]
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2)/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1)/det
                           if( .false. )then
                             ! check the error
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                             write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2," u(-1)=",e10.3," err=",e8.2)') i1,i2,uTemp(j1,j2,j3,0),abs(uTemp(j1,j2,j3,0)-ue1),uTemp(k1,k2,k3,0),abs(uTemp(k1,k2,k3,0)-ue2)
                           end if
                         end if ! mask .gt. 0
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
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
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ! --- get the compatibility forcings for order=2 ---
                             ! No forcing, do nothing 
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22 or ulaplacian23 
                           uLap = ulaplacian22(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2
                             a11 = ( rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2 )/( dr(axis)**2 )  ! coeff of u(-1) in r1 
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                         if( .false. )then
                           ! check the error
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ghost = 1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                       end if ! mask 
                     end do
                     end do
                     end do
                       !----------------------------------------------------------------
                       ! --- STAGE II fill in first ghost by 4th-order compatibility ---
                       !----------------------------------------------------------------
                         ! Compute and save v = Lap(u) at some points near the boundary
                         !  for use below to compute Lap^2 (u)
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !         r=0 
                         do i3=m3a,m3b
                         do i2=m2a,m2b
                         do i1=m1a,m1b
                           ! eval Lap(u) on ghost, boundary and first line in:
                           if( mask(i1,i2,i3).ne.0 )then 
                             v(i1-is1,i2-is2,i3-is3,0) = ulaplacian22(i1-is1,i2-is2,i3-is3,0)
                             v(i1    ,i2    ,i3    ,0) = ulaplacian22(i1    ,i2    ,i3    ,0)
                             v(i1+is1,i2+is2,i3+is3,0) = ulaplacian22(i1+is1,i2+is2,i3+is3,0)
                           end if
                         end do
                         end do
                         end do
                       ! --- fill in two ghost using 4th- order compatibility ---
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ! --- get the compatibility forcings at order=4 ---
                               ! No forcing, do nothing 
                           ! u_tt = c^2*Lap(u) + f 
                           ! u_tttt = c^2*Lap(u_tt) + f_tt 
                           !        = c^2*Lap( c^2*Lap(u) + f ) + f_tt 
                           !        = c^4*Lap^2(u) + c^2*Lap(f) ) + f_tt 
                             ! curvilinear 
                             ! uxx42(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2)*urr4(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0))*urs4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)**2)*uss4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,0,0))*ur4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,1,0))*us4(i1,i2,i3,kd)
                             uLap = ulaplacian42(i1,i2,i3,0)
                             r1 =  c2*uLap + ff - gtt            ! residual in equation using current ghost value
                               rFactor = rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2
                             a11 =  c2*rFactor*16./(12.*dr(axis)**2)                           ! coeff of u(-1) in r1 
                             a12 = -c2*rFactor    /(12.*dr(axis)**2)                           ! coeff of u(-2) in r1
                             ! vLap = Lap^2( u) to second-order
                             ! vLap = vlaplacian22 or vlaplacian23 
                             vLap = vlaplacian22(i1,i2,i3,0)
                             r2 =  c4*vLap + c2*fLap + ftt - gtttt
                             ! Here is the leading order term in a21, a22   ** this may not be good enough for stability but should remain accurate ***
                             a21 =  -c4*( rFactor**2*4./(dr(axis)**4) )                             ! coeff of u(-1) in r2
                             a22 =   c4*( rFactor**2   /(dr(axis)**4) )         
                         if( .false. )then
                           write(*,'("(i1,i2)=(",2i3,"),  r1,r2,a11,a12,a21,a22,vLap=",7(1pe9.2,1x))') i1,i2,r1,r2,a11,a12,a21,a22,vLap
                         end if
                           ghost =1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost =2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost        
                           ! Solve
                           !  [ a11 a12 ][ u(-1)] = [ a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1 ]
                           !  [ a21 a22 ][ u(-2)]   [ a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2 ]
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2)/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1)/det
                           if( .false. )then
                             ! check the error
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                             write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2," u(-1)=",e10.3," err=",e8.2)') i1,i2,uTemp(j1,j2,j3,0),abs(uTemp(j1,j2,j3,0)-ue1),uTemp(k1,k2,k3,0),abs(uTemp(k1,k2,k3,0)-ue2)
                           end if
                         end if ! mask .gt. 0
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
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
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
                         ! --- get the compatibility forcings for order=2 ---
                             ! No forcing, do nothing 
                         ! u_tt = c^2*Lap(u) + f 
                           ! uLap = ulaplacian22 or ulaplacian23 
                           uLap = ulaplacian23(i1,i2,i3,0)
                           r1 = uLap + (ff - gtt)/c2
                             a11 = ( rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2 + rsxy(i1,i2,i3,axis,2)**2 )/( dr(axis)**2 )  ! coeff of u(-1) in r1 
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                         if( .false. )then
                           ! check the error
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ghost = 1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                       end if ! mask 
                     end do
                     end do
                     end do
                       !----------------------------------------------------------------
                       ! --- STAGE II fill in first ghost by 4th-order compatibility ---
                       !----------------------------------------------------------------
                         ! Compute and save v = Lap(u) at some points near the boundary
                         !  for use below to compute Lap^2 (u)
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !         r=0 
                         do i3=m3a,m3b
                         do i2=m2a,m2b
                         do i1=m1a,m1b
                           ! eval Lap(u) on ghost, boundary and first line in:
                           if( mask(i1,i2,i3).ne.0 )then 
                             v(i1-is1,i2-is2,i3-is3,0) = ulaplacian23(i1-is1,i2-is2,i3-is3,0)
                             v(i1    ,i2    ,i3    ,0) = ulaplacian23(i1    ,i2    ,i3    ,0)
                             v(i1+is1,i2+is2,i3+is3,0) = ulaplacian23(i1+is1,i2+is2,i3+is3,0)
                           end if
                         end do
                         end do
                         end do
                       ! --- fill in two ghost using 4th- order compatibility ---
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ! --- get the compatibility forcings at order=4 ---
                               ! No forcing, do nothing 
                           ! u_tt = c^2*Lap(u) + f 
                           ! u_tttt = c^2*Lap(u_tt) + f_tt 
                           !        = c^2*Lap( c^2*Lap(u) + f ) + f_tt 
                           !        = c^4*Lap^2(u) + c^2*Lap(f) ) + f_tt 
                             ! curvilinear 
                             ! uxx42(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2)*urr4(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0))*urs4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)**2)*uss4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,0,0))*ur4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,1,0))*us4(i1,i2,i3,kd)
                             uLap = ulaplacian43(i1,i2,i3,0)
                             r1 =  c2*uLap + ff - gtt            ! residual in equation using current ghost value
                               rFactor = rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2 + rsxy(i1,i2,i3,axis,2)**2
                             a11 =  c2*rFactor*16./(12.*dr(axis)**2)                           ! coeff of u(-1) in r1 
                             a12 = -c2*rFactor    /(12.*dr(axis)**2)                           ! coeff of u(-2) in r1
                             ! vLap = Lap^2( u) to second-order
                             ! vLap = vlaplacian22 or vlaplacian23 
                             vLap = vlaplacian23(i1,i2,i3,0)
                             r2 =  c4*vLap + c2*fLap + ftt - gtttt
                             ! Here is the leading order term in a21, a22   ** this may not be good enough for stability but should remain accurate ***
                             a21 =  -c4*( rFactor**2*4./(dr(axis)**4) )                             ! coeff of u(-1) in r2
                             a22 =   c4*( rFactor**2   /(dr(axis)**4) )         
                         if( .false. )then
                           write(*,'("(i1,i2)=(",2i3,"),  r1,r2,a11,a12,a21,a22,vLap=",7(1pe9.2,1x))') i1,i2,r1,r2,a11,a12,a21,a22,vLap
                         end if
                           ghost =1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost =2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost        
                           ! Solve
                           !  [ a11 a12 ][ u(-1)] = [ a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1 ]
                           !  [ a21 a22 ][ u(-2)]   [ a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2 ]
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2)/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1)/det
                           if( .false. )then
                             ! check the error
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                             write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2," u(-1)=",e10.3," err=",e8.2)') i1,i2,uTemp(j1,j2,j3,0),abs(uTemp(j1,j2,j3,0)-ue1),uTemp(k1,k2,k3,0),abs(uTemp(k1,k2,k3,0)-ue2)
                           end if
                         end if ! mask .gt. 0
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
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
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
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
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                         if( .false. )then
                           ! check the error
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ghost = 1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                       end if ! mask 
                     end do
                     end do
                     end do
                       !----------------------------------------------------------------
                       ! --- STAGE II fill in first ghost by 4th-order compatibility ---
                       !----------------------------------------------------------------
                       ! --- fill in two ghost using 4th- order compatibility ---
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ! --- get the compatibility forcings at order=4 ---
                               if( assignTwilightZone.eq.1 )then
                                 ! compute RHS from TZ
                                   call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexx )
                                   call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyy )
                                   ueLap = uexx + ueyy
                                   call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uett )
                                     call ogDeriv(ep,2,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uettxx )
                                     call ogDeriv(ep,2,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uettyy )
                                     uettLap = uettxx + uettyy
                                     call ogDeriv(ep,0,4,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxx )
                                     call ogDeriv(ep,0,2,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxyy )
                                     call ogDeriv(ep,0,0,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyyy )
                                     ueLap2 = uexxxx + 2.*uexxyy + ueyyyy ! Lap^2( u )
                                     call ogDeriv(ep,4,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uetttt )
                                 ff = uett - c2*ueLap
                                 gtt = uett
                                   fLap = uettLap - c2*ueLap2 
                                   ftt  = uetttt - c2*uettLap
                                   gtttt = uetttt
                               else
                                 ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.;
                               end if
                           ! u_tt = c^2*Lap(u) + f 
                           ! u_tttt = c^2*Lap(u_tt) + f_tt 
                           !        = c^2*Lap( c^2*Lap(u) + f ) + f_tt 
                           !        = c^4*Lap^2(u) + c^2*Lap(f) ) + f_tt 
                             ! uxx43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*dx42(0) 
                             uLap = ulaplacian42r(i1,i2,i3,0)
                             r1 =  c2*uLap + ff - gtt            ! residual in equation using current ghost value
                             a11 = c2*16./(12.*dx(axis)**2)                           ! coeff of u(-1) in r1 
                             a12 =    -c2/(12.*dx(axis)**2)                           ! coeff of u(-2) in r1
                             ! uxxxx = 1 -4 6 -4 1
                             ! vLap = Lap^2( u )
                             vLap = lap2d2Pow2(i1,i2,i3,0)
                             r2 =  c4*vLap + c2*fLap + ftt - gtttt
                             a21 =  -c4*4./(dx(axis)**4)                              ! coeff of u(-1) in r2
                             a22 =   c4*1./(dx(axis)**4)
                         if( .false. )then
                           write(*,'("(i1,i2)=(",2i3,"),  r1,r2,a11,a12,a21,a22,vLap=",7(1pe9.2,1x))') i1,i2,r1,r2,a11,a12,a21,a22,vLap
                         end if
                           ghost =1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost =2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost        
                           ! Solve
                           !  [ a11 a12 ][ u(-1)] = [ a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1 ]
                           !  [ a21 a22 ][ u(-2)]   [ a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2 ]
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2)/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1)/det
                           if( .false. )then
                             ! check the error
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                             write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2," u(-1)=",e10.3," err=",e8.2)') i1,i2,uTemp(j1,j2,j3,0),abs(uTemp(j1,j2,j3,0)-ue1),uTemp(k1,k2,k3,0),abs(uTemp(k1,k2,k3,0)-ue2)
                           end if
                         end if ! mask .gt. 0
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
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
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
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
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                         if( .false. )then
                           ! check the error
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ghost = 1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                       end if ! mask 
                     end do
                     end do
                     end do
                       !----------------------------------------------------------------
                       ! --- STAGE II fill in first ghost by 4th-order compatibility ---
                       !----------------------------------------------------------------
                       ! --- fill in two ghost using 4th- order compatibility ---
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ! --- get the compatibility forcings at order=4 ---
                               if( assignTwilightZone.eq.1 )then
                                 ! compute RHS from TZ
                                   ! 3D 
                                   call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexx )
                                   call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyy )
                                   call ogDeriv(ep,0,0,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezz )
                                   ueLap = uexx + ueyy + uezz
                                   call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uett )
                                     call ogDeriv(ep,2,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettxx )
                                     call ogDeriv(ep,2,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettyy )
                                     call ogDeriv(ep,2,0,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettzz )
                                     uettLap = uettxx + uettyy + uettzz
                                     call ogDeriv(ep,0,4,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxxx )
                                     call ogDeriv(ep,0,0,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyyy )
                                     call ogDeriv(ep,0,0,0,4,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezzzz )
                                     call ogDeriv(ep,0,2,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxyy )
                                     call ogDeriv(ep,0,2,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxzz )
                                     call ogDeriv(ep,0,0,2,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyzz )
                                     ueLap2 = uexxxx + 2.*( uexxyy + uexxzz + ueyyzz)  + ueyyyy + uezzzz ! Lap^2( u )
                                     call ogDeriv(ep,4,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uetttt )
                                 ff = uett - c2*ueLap
                                 gtt = uett
                                   fLap = uettLap - c2*ueLap2 
                                   ftt  = uetttt - c2*uettLap
                                   gtttt = uetttt
                               else
                                 ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.;
                               end if
                           ! u_tt = c^2*Lap(u) + f 
                           ! u_tttt = c^2*Lap(u_tt) + f_tt 
                           !        = c^2*Lap( c^2*Lap(u) + f ) + f_tt 
                           !        = c^4*Lap^2(u) + c^2*Lap(f) ) + f_tt 
                             ! uxx43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*dx42(0) 
                             uLap = ulaplacian43r(i1,i2,i3,0)
                             r1 =  c2*uLap + ff - gtt            ! residual in equation using current ghost value
                             a11 = c2*16./(12.*dx(axis)**2)                           ! coeff of u(-1) in r1 
                             a12 =    -c2/(12.*dx(axis)**2)                           ! coeff of u(-2) in r1
                             ! uxxxx = 1 -4 6 -4 1
                             ! vLap = Lap^2( u )
                             vLap = lap3d2Pow2(i1,i2,i3,0)
                             r2 =  c4*vLap + c2*fLap + ftt - gtttt
                             a21 =  -c4*4./(dx(axis)**4)                              ! coeff of u(-1) in r2
                             a22 =   c4*1./(dx(axis)**4)
                         if( .false. )then
                           write(*,'("(i1,i2)=(",2i3,"),  r1,r2,a11,a12,a21,a22,vLap=",7(1pe9.2,1x))') i1,i2,r1,r2,a11,a12,a21,a22,vLap
                         end if
                           ghost =1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost =2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost        
                           ! Solve
                           !  [ a11 a12 ][ u(-1)] = [ a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1 ]
                           !  [ a21 a22 ][ u(-2)]   [ a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2 ]
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2)/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1)/det
                           if( .false. )then
                             ! check the error
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                             write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2," u(-1)=",e10.3," err=",e8.2)') i1,i2,uTemp(j1,j2,j3,0),abs(uTemp(j1,j2,j3,0)-ue1),uTemp(k1,k2,k3,0),abs(uTemp(k1,k2,k3,0)-ue2)
                           end if
                         end if ! mask .gt. 0
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
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
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
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
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                         if( .false. )then
                           ! check the error
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ghost = 1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                       end if ! mask 
                     end do
                     end do
                     end do
                       !----------------------------------------------------------------
                       ! --- STAGE II fill in first ghost by 4th-order compatibility ---
                       !----------------------------------------------------------------
                         ! Compute and save v = Lap(u) at some points near the boundary
                         !  for use below to compute Lap^2 (u)
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !         r=0 
                         do i3=m3a,m3b
                         do i2=m2a,m2b
                         do i1=m1a,m1b
                           ! eval Lap(u) on ghost, boundary and first line in:
                           if( mask(i1,i2,i3).ne.0 )then 
                             v(i1-is1,i2-is2,i3-is3,0) = ulaplacian22(i1-is1,i2-is2,i3-is3,0)
                             v(i1    ,i2    ,i3    ,0) = ulaplacian22(i1    ,i2    ,i3    ,0)
                             v(i1+is1,i2+is2,i3+is3,0) = ulaplacian22(i1+is1,i2+is2,i3+is3,0)
                           end if
                         end do
                         end do
                         end do
                       ! --- fill in two ghost using 4th- order compatibility ---
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ! --- get the compatibility forcings at order=4 ---
                               if( assignTwilightZone.eq.1 )then
                                 ! compute RHS from TZ
                                   call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexx )
                                   call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyy )
                                   ueLap = uexx + ueyy
                                   call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uett )
                                     call ogDeriv(ep,2,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uettxx )
                                     call ogDeriv(ep,2,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uettyy )
                                     uettLap = uettxx + uettyy
                                     call ogDeriv(ep,0,4,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxx )
                                     call ogDeriv(ep,0,2,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxyy )
                                     call ogDeriv(ep,0,0,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyyy )
                                     ueLap2 = uexxxx + 2.*uexxyy + ueyyyy ! Lap^2( u )
                                     call ogDeriv(ep,4,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uetttt )
                                 ff = uett - c2*ueLap
                                 gtt = uett
                                   fLap = uettLap - c2*ueLap2 
                                   ftt  = uetttt - c2*uettLap
                                   gtttt = uetttt
                               else
                                 ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.;
                               end if
                           ! u_tt = c^2*Lap(u) + f 
                           ! u_tttt = c^2*Lap(u_tt) + f_tt 
                           !        = c^2*Lap( c^2*Lap(u) + f ) + f_tt 
                           !        = c^4*Lap^2(u) + c^2*Lap(f) ) + f_tt 
                             ! curvilinear 
                             ! uxx42(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2)*urr4(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0))*urs4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)**2)*uss4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,0,0))*ur4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,1,0))*us4(i1,i2,i3,kd)
                             uLap = ulaplacian42(i1,i2,i3,0)
                             r1 =  c2*uLap + ff - gtt            ! residual in equation using current ghost value
                               rFactor = rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2
                             a11 =  c2*rFactor*16./(12.*dr(axis)**2)                           ! coeff of u(-1) in r1 
                             a12 = -c2*rFactor    /(12.*dr(axis)**2)                           ! coeff of u(-2) in r1
                             ! vLap = Lap^2( u) to second-order
                             ! vLap = vlaplacian22 or vlaplacian23 
                             vLap = vlaplacian22(i1,i2,i3,0)
                             r2 =  c4*vLap + c2*fLap + ftt - gtttt
                             ! Here is the leading order term in a21, a22   ** this may not be good enough for stability but should remain accurate ***
                             a21 =  -c4*( rFactor**2*4./(dr(axis)**4) )                             ! coeff of u(-1) in r2
                             a22 =   c4*( rFactor**2   /(dr(axis)**4) )         
                         if( .false. )then
                           write(*,'("(i1,i2)=(",2i3,"),  r1,r2,a11,a12,a21,a22,vLap=",7(1pe9.2,1x))') i1,i2,r1,r2,a11,a12,a21,a22,vLap
                         end if
                           ghost =1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost =2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost        
                           ! Solve
                           !  [ a11 a12 ][ u(-1)] = [ a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1 ]
                           !  [ a21 a22 ][ u(-2)]   [ a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2 ]
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2)/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1)/det
                           if( .false. )then
                             ! check the error
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                             write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2," u(-1)=",e10.3," err=",e8.2)') i1,i2,uTemp(j1,j2,j3,0),abs(uTemp(j1,j2,j3,0)-ue1),uTemp(k1,k2,k3,0),abs(uTemp(k1,k2,k3,0)-ue2)
                           end if
                         end if ! mask .gt. 0
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
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
                     ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.; 
                     do i3=m3a,m3b
                     do i2=m2a,m2b
                     do i1=m1a,m1b
                       if( mask(i1,i2,i3).gt.0 )then ! valid discretization point on the boundary 
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
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,0) = u(j1,j2,j3,0) - r1/a11 
                         if( .false. )then
                           ! check the error
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2, "  (order 2 update)")') i1,i2,u(j1,j2,j3,0),abs(u(j1,j2,j3,0)-ue1)
                         end if
                       else if( mask(i1,i2,i3).lt.0 )then ! interpolation point on the boundary 
                         ! ----- extrap ghost outside interp. pts on physical boundaries ------
                         ghost = 1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         u(j1,j2,j3,uc)=(3.*u(i1,i2,i3,uc)-3.*u(i1+is1,i2+is2,i3+is3,uc)+u(i1+2*is1,i2+2*is2,i3+2*is3,uc))
                       end if ! mask 
                     end do
                     end do
                     end do
                       !----------------------------------------------------------------
                       ! --- STAGE II fill in first ghost by 4th-order compatibility ---
                       !----------------------------------------------------------------
                         ! Compute and save v = Lap(u) at some points near the boundary
                         !  for use below to compute Lap^2 (u)
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !        +-X-+
                         !          |
                         !         r=0 
                         do i3=m3a,m3b
                         do i2=m2a,m2b
                         do i1=m1a,m1b
                           ! eval Lap(u) on ghost, boundary and first line in:
                           if( mask(i1,i2,i3).ne.0 )then 
                             v(i1-is1,i2-is2,i3-is3,0) = ulaplacian23(i1-is1,i2-is2,i3-is3,0)
                             v(i1    ,i2    ,i3    ,0) = ulaplacian23(i1    ,i2    ,i3    ,0)
                             v(i1+is1,i2+is2,i3+is3,0) = ulaplacian23(i1+is1,i2+is2,i3+is3,0)
                           end if
                         end do
                         end do
                         end do
                       ! --- fill in two ghost using 4th- order compatibility ---
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ! --- get the compatibility forcings at order=4 ---
                               if( assignTwilightZone.eq.1 )then
                                 ! compute RHS from TZ
                                   ! 3D 
                                   call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexx )
                                   call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyy )
                                   call ogDeriv(ep,0,0,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezz )
                                   ueLap = uexx + ueyy + uezz
                                   call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uett )
                                     call ogDeriv(ep,2,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettxx )
                                     call ogDeriv(ep,2,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettyy )
                                     call ogDeriv(ep,2,0,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettzz )
                                     uettLap = uettxx + uettyy + uettzz
                                     call ogDeriv(ep,0,4,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxxx )
                                     call ogDeriv(ep,0,0,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyyy )
                                     call ogDeriv(ep,0,0,0,4,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezzzz )
                                     call ogDeriv(ep,0,2,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxyy )
                                     call ogDeriv(ep,0,2,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxzz )
                                     call ogDeriv(ep,0,0,2,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyzz )
                                     ueLap2 = uexxxx + 2.*( uexxyy + uexxzz + ueyyzz)  + ueyyyy + uezzzz ! Lap^2( u )
                                     call ogDeriv(ep,4,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uetttt )
                                 ff = uett - c2*ueLap
                                 gtt = uett
                                   fLap = uettLap - c2*ueLap2 
                                   ftt  = uetttt - c2*uettLap
                                   gtttt = uetttt
                               else
                                 ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.;
                               end if
                           ! u_tt = c^2*Lap(u) + f 
                           ! u_tttt = c^2*Lap(u_tt) + f_tt 
                           !        = c^2*Lap( c^2*Lap(u) + f ) + f_tt 
                           !        = c^4*Lap^2(u) + c^2*Lap(f) ) + f_tt 
                             ! curvilinear 
                             ! uxx42(i1,i2,i3,kd)=(rsxy(i1,i2,i3,0,0)**2)*urr4(i1,i2,i3,kd)+2.*(rsxy(i1,i2,i3,0,0)*rsxy(i1,i2,i3,1,0))*urs4(i1,i2,i3,kd)+(rsxy(i1,i2,i3,1,0)**2)*uss4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,0,0))*ur4(i1,i2,i3,kd)+(rsxyx42(i1,i2,i3,1,0))*us4(i1,i2,i3,kd)
                             uLap = ulaplacian43(i1,i2,i3,0)
                             r1 =  c2*uLap + ff - gtt            ! residual in equation using current ghost value
                               rFactor = rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2 + rsxy(i1,i2,i3,axis,2)**2
                             a11 =  c2*rFactor*16./(12.*dr(axis)**2)                           ! coeff of u(-1) in r1 
                             a12 = -c2*rFactor    /(12.*dr(axis)**2)                           ! coeff of u(-2) in r1
                             ! vLap = Lap^2( u) to second-order
                             ! vLap = vlaplacian22 or vlaplacian23 
                             vLap = vlaplacian23(i1,i2,i3,0)
                             r2 =  c4*vLap + c2*fLap + ftt - gtttt
                             ! Here is the leading order term in a21, a22   ** this may not be good enough for stability but should remain accurate ***
                             a21 =  -c4*( rFactor**2*4./(dr(axis)**4) )                             ! coeff of u(-1) in r2
                             a22 =   c4*( rFactor**2   /(dr(axis)**4) )         
                         if( .false. )then
                           write(*,'("(i1,i2)=(",2i3,"),  r1,r2,a11,a12,a21,a22,vLap=",7(1pe9.2,1x))') i1,i2,r1,r2,a11,a12,a21,a22,vLap
                         end if
                           ghost =1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost =2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost        
                           ! Solve
                           !  [ a11 a12 ][ u(-1)] = [ a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1 ]
                           !  [ a21 a22 ][ u(-2)]   [ a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2 ]
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2)/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1)/det
                           if( .false. )then
                             ! check the error
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                             write(*,'("(i1,i2)=(",2i3,") u(-1)=",e10.3," err=",e8.2," u(-1)=",e10.3," err=",e8.2)') i1,i2,uTemp(j1,j2,j3,0),abs(uTemp(j1,j2,j3,0)-ue1),uTemp(k1,k2,k3,0),abs(uTemp(k1,k2,k3,0)-ue2)
                           end if
                         end if ! mask .gt. 0
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
                        end do
                        end do
                        end do
                 end if
               end if     
             end if
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
         else if( orderOfAccuracy.eq.4 )then
             if( forcingOption.eq.noForcing )then
               if( gridType.eq.rectangular )then
                 if( nd.eq.2 )then
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
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       ! ------------------- 4 4 NEUMANN CBC --------------------
                       maxDiff=0. ! for checkCoeff
                       gg=0.
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         ghost =2
                         k1=i1-is1*ghost
                         k2=i2-is2*ghost
                         k3=i3-is3*ghost   
                         if( mask(i1,i2,i3).gt.0 )then
                               ! No forcing, do nothing 
                             ! ---------- RECTANGULAR  ---------
                             ! This assumes dr(0:2) = dx(0:2)
                             ! --- NEUMANN 4=4 rectangular ---
                             ! u_tt = c^2*Lap(u) + f 
                             !   u.n = g
                             ! g_tt = c^2 n.grad( Lap(u) ) + n.grad(f)
                             ! ux = [ u(-2) - 8*u(-1) + 8*u(1) - u(2) ]/(12*h) 
                             ! uxxx = (-u(-1) +2*u(-1) - 2*u(1) + u(2) ]/(2*h^3)
                               ! eval equation with wrong values at ghost: 
                               r1 =  an1*ux42r(i1,i2,i3,0) + an2*uy42r(i1,i2,i3,0) - gg
                               ! note: an1=-1 on left side ans an1=+1 on right side
                               ! **CHECK ME**
                               a11 =  8./(12.*dx(axis))  ! coeff of u(-1)
                               a12 = -1./(12.*dx(axis))  ! coeff of u(-2)
                               r2 = c2*( an1*( uxxx + uxyy ) + an2*( uxxy + uyyy ) ) + nDotGradF - gtt
                               ! **CHECK ME**
                               a21 = -c2*( 2./(2.*dx(axis)**3) )
                               a22 =  c2*( 1./(2.*dx(axis)**3) )
                               ! define the residual functions for the discrete delta method
                                 if( checkCoeff.eq.1 )then
                                   ! --- Check coefficients a11,a12,... by discrete delta ---
                                   u1Save = u(j1,j2,j3,0)
                                   u2Save = u(k1,k2,k3,0)
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=0.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   r1a = (an1*ux42r(i1,i2,i3,0)+an2*uy42r(i1,i2,i3,0))
                                   r2a = (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   u(j1,j2,j3,0)=1.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   r1b =  (an1*ux42r(i1,i2,i3,0)+an2*uy42r(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   a11c = r1b - r1a
                                   a21c = r2b - r2a
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=1.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   r1b =  (an1*ux42r(i1,i2,i3,0)+an2*uy42r(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   a12c = r1b - r1a
                                   a22c = r2b - r2a
                                   u(j1,j2,j3,0) = u1Save  
                                   u(k1,k2,k3,0) = u2Save      
                                   maxDiff=max(maxDiff,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c),abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c))
                                   ! #If "rectangular" eq "curvilinear"
                                   ! write(*,'("++ i1,i2=",2i4," a11,a11c=",2(1pe12.4)," a12,a12c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a11,a11c,a12,a12c,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c)
                                   ! write(*,'("++ i1,i2=",2i4," a21,a21c=",2(1pe12.4)," a22,a22c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a21,a21c,a22,a22c,abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c)
                                   ! #End
                                 end if 
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2 )/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1 )/det
                           if( .false. )then
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,uex )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,uey )
                             write(*,'(" (i1,i2)=(",2i3,") u(-1),ue(-1)=",2(1pe12.4)," u(-2),ue(-2)=",2(1pe12.4))') i1,i2,uTemp(j1,j2,j3,0),uex,uTemp(k1,k2,k3,0),uey
                           end if
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! this case done below 
                         end if
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
                        end do
                        end do
                        end do
                       if( checkCoeff.eq.1 )then
                         write(*,'("bcOptWave: t=",e12.3," checkCoeff: (side,axis)=(",2i2,") rel-maxDiff=",e9.2," rectangular 4")') t,side,axis,maxDiff
                       end if
                 else
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
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       ! ------------------- 4 4 NEUMANN CBC --------------------
                       maxDiff=0. ! for checkCoeff
                       gg=0.
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         ghost =2
                         k1=i1-is1*ghost
                         k2=i2-is2*ghost
                         k3=i3-is3*ghost   
                         if( mask(i1,i2,i3).gt.0 )then
                               ! No forcing, do nothing 
                             ! evaluate 3rd derivatives : uxxx,uxxy,uxxz,uxyy,uxzz, uyyy,uyyz,uyzz, uzzz 
                             ! ---------- RECTANGULAR  ---------
                             ! This assumes dr(0:2) = dx(0:2)
                             uxyz     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                             ! --- NEUMANN 4=4 rectangular ---
                             ! u_tt = c^2*Lap(u) + f 
                             !   u.n = g
                             ! g_tt = c^2 n.grad( Lap(u) ) + n.grad(f)
                             ! ux = [ u(-2) - 8*u(-1) + 8*u(1) - u(2) ]/(12*h) 
                             ! uxxx = (-u(-1) +2*u(-1) - 2*u(1) + u(2) ]/(2*h^3)
                               ! 3D 
                               r1 =  an1*ux43r(i1,i2,i3,0) + an2*uy43r(i1,i2,i3,0) + an3*uz43r(i1,i2,i3,0) - gg
                               ! **CHECK ME**
                               a11 =  8./(12.*dx(axis))  ! coeff of u(-1)
                               a12 = -1./(12.*dx(axis))  ! coeff of u(-2)  
                               r2 = c2*( an1*( uxxx + uxyy + uxzz ) + an2*( uxxy + uyyy + uyzz ) + an3*( uxxz + uyyz + uzzz ) ) + nDotGradF - gtt
                               ! **CHECK ME**
                               a21 = -c2*( 2./(2.*dx(axis)**3) )
                               a22 =  c2*( 1./(2.*dx(axis)**3) )
                               ! define the residual functions for the discrete delta method
                                 if( checkCoeff.eq.1 )then
                                   ! --- Check coefficients a11,a12,... by discrete delta ---
                                   u1Save = u(j1,j2,j3,0)
                                   u2Save = u(k1,k2,k3,0)
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=0.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   uxyz     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   r1a = (an1*ux43r(i1,i2,i3,0)+an2*uy43r(i1,i2,i3,0)+an3*uz43r(i1,i2,i3,0))
                                   r2a = (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   u(j1,j2,j3,0)=1.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   uxyz     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   r1b =  (an1*ux43r(i1,i2,i3,0)+an2*uy43r(i1,i2,i3,0)+an3*uz43r(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   a11c = r1b - r1a
                                   a21c = r2b - r2a
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=1.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   uxyz     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   r1b =  (an1*ux43r(i1,i2,i3,0)+an2*uy43r(i1,i2,i3,0)+an3*uz43r(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   a12c = r1b - r1a
                                   a22c = r2b - r2a
                                   u(j1,j2,j3,0) = u1Save  
                                   u(k1,k2,k3,0) = u2Save      
                                   maxDiff=max(maxDiff,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c),abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c))
                                   ! #If "rectangular" eq "curvilinear"
                                   ! write(*,'("++ i1,i2=",2i4," a11,a11c=",2(1pe12.4)," a12,a12c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a11,a11c,a12,a12c,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c)
                                   ! write(*,'("++ i1,i2=",2i4," a21,a21c=",2(1pe12.4)," a22,a22c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a21,a21c,a22,a22c,abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c)
                                   ! #End
                                 end if 
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2 )/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1 )/det
                           if( .false. )then
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,uex )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,uey )
                             write(*,'(" (i1,i2)=(",2i3,") u(-1),ue(-1)=",2(1pe12.4)," u(-2),ue(-2)=",2(1pe12.4))') i1,i2,uTemp(j1,j2,j3,0),uex,uTemp(k1,k2,k3,0),uey
                           end if
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! this case done below 
                         end if
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
                        end do
                        end do
                        end do
                       if( checkCoeff.eq.1 )then
                         write(*,'("bcOptWave: t=",e12.3," checkCoeff: (side,axis)=(",2i2,") rel-maxDiff=",e9.2," rectangular 4")') t,side,axis,maxDiff
                       end if
                 end if
               else
                 if( nd.eq.2 )then
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       ! ------------------- 4 4 NEUMANN CBC --------------------
                       maxDiff=0. ! for checkCoeff
                       gg=0.
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         ghost =2
                         k1=i1-is1*ghost
                         k2=i2-is2*ghost
                         k3=i3-is3*ghost   
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
                               ! No forcing, do nothing 
                             ! ---------- START CURVILINEAR  ---------
                             ! ---------- Parametric derivatives ---------
                             ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                             urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                             urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                             us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                             urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                             urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                             uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                             urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                             usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                             rxr      = (rx(i1-2,i2,i3)-8.*rx(i1-1,i2,i3)+8.*rx(i1+1,i2,i3)-rx(i1+2,i2,i3))/(12.*dr(0))
                             rxrr     = (-rx(i1-2,i2,i3)+16.*rx(i1-1,i2,i3)-30.*rx(i1,i2,i3)+16.*rx(i1+1,i2,i3)-rx(i1+2,i2,i3))/(12.*dr(0)**2)
                             rxs      = (rx(i1,i2-2,i3)-8.*rx(i1,i2-1,i3)+8.*rx(i1,i2+1,i3)-rx(i1,i2+2,i3))/(12.*dr(1))
                             rxrs     = ((rx(i1-2,i2-2,i3)-8.*rx(i1-2,i2-1,i3)+8.*rx(i1-2,i2+1,i3)-rx(i1-2,i2+2,i3))/(12.*dr(1))-8.*(rx(i1-1,i2-2,i3)-8.*rx(i1-1,i2-1,i3)+8.*rx(i1-1,i2+1,i3)-rx(i1-1,i2+2,i3))/(12.*dr(1))+8.*(rx(i1+1,i2-2,i3)-8.*rx(i1+1,i2-1,i3)+8.*rx(i1+1,i2+1,i3)-rx(i1+1,i2+2,i3))/(12.*dr(1))-(rx(i1+2,i2-2,i3)-8.*rx(i1+2,i2-1,i3)+8.*rx(i1+2,i2+1,i3)-rx(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             rxss     = (-rx(i1,i2-2,i3)+16.*rx(i1,i2-1,i3)-30.*rx(i1,i2,i3)+16.*rx(i1,i2+1,i3)-rx(i1,i2+2,i3))/(12.*dr(1)**2)
                             ryr      = (ry(i1-2,i2,i3)-8.*ry(i1-1,i2,i3)+8.*ry(i1+1,i2,i3)-ry(i1+2,i2,i3))/(12.*dr(0))
                             ryrr     = (-ry(i1-2,i2,i3)+16.*ry(i1-1,i2,i3)-30.*ry(i1,i2,i3)+16.*ry(i1+1,i2,i3)-ry(i1+2,i2,i3))/(12.*dr(0)**2)
                             rys      = (ry(i1,i2-2,i3)-8.*ry(i1,i2-1,i3)+8.*ry(i1,i2+1,i3)-ry(i1,i2+2,i3))/(12.*dr(1))
                             ryrs     = ((ry(i1-2,i2-2,i3)-8.*ry(i1-2,i2-1,i3)+8.*ry(i1-2,i2+1,i3)-ry(i1-2,i2+2,i3))/(12.*dr(1))-8.*(ry(i1-1,i2-2,i3)-8.*ry(i1-1,i2-1,i3)+8.*ry(i1-1,i2+1,i3)-ry(i1-1,i2+2,i3))/(12.*dr(1))+8.*(ry(i1+1,i2-2,i3)-8.*ry(i1+1,i2-1,i3)+8.*ry(i1+1,i2+1,i3)-ry(i1+1,i2+2,i3))/(12.*dr(1))-(ry(i1+2,i2-2,i3)-8.*ry(i1+2,i2-1,i3)+8.*ry(i1+2,i2+1,i3)-ry(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             ryss     = (-ry(i1,i2-2,i3)+16.*ry(i1,i2-1,i3)-30.*ry(i1,i2,i3)+16.*ry(i1,i2+1,i3)-ry(i1,i2+2,i3))/(12.*dr(1)**2)
                             sxr      = (sx(i1-2,i2,i3)-8.*sx(i1-1,i2,i3)+8.*sx(i1+1,i2,i3)-sx(i1+2,i2,i3))/(12.*dr(0))
                             sxrr     = (-sx(i1-2,i2,i3)+16.*sx(i1-1,i2,i3)-30.*sx(i1,i2,i3)+16.*sx(i1+1,i2,i3)-sx(i1+2,i2,i3))/(12.*dr(0)**2)
                             sxs      = (sx(i1,i2-2,i3)-8.*sx(i1,i2-1,i3)+8.*sx(i1,i2+1,i3)-sx(i1,i2+2,i3))/(12.*dr(1))
                             sxrs     = ((sx(i1-2,i2-2,i3)-8.*sx(i1-2,i2-1,i3)+8.*sx(i1-2,i2+1,i3)-sx(i1-2,i2+2,i3))/(12.*dr(1))-8.*(sx(i1-1,i2-2,i3)-8.*sx(i1-1,i2-1,i3)+8.*sx(i1-1,i2+1,i3)-sx(i1-1,i2+2,i3))/(12.*dr(1))+8.*(sx(i1+1,i2-2,i3)-8.*sx(i1+1,i2-1,i3)+8.*sx(i1+1,i2+1,i3)-sx(i1+1,i2+2,i3))/(12.*dr(1))-(sx(i1+2,i2-2,i3)-8.*sx(i1+2,i2-1,i3)+8.*sx(i1+2,i2+1,i3)-sx(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             sxss     = (-sx(i1,i2-2,i3)+16.*sx(i1,i2-1,i3)-30.*sx(i1,i2,i3)+16.*sx(i1,i2+1,i3)-sx(i1,i2+2,i3))/(12.*dr(1)**2)
                             syr      = (sy(i1-2,i2,i3)-8.*sy(i1-1,i2,i3)+8.*sy(i1+1,i2,i3)-sy(i1+2,i2,i3))/(12.*dr(0))
                             syrr     = (-sy(i1-2,i2,i3)+16.*sy(i1-1,i2,i3)-30.*sy(i1,i2,i3)+16.*sy(i1+1,i2,i3)-sy(i1+2,i2,i3))/(12.*dr(0)**2)
                             sys      = (sy(i1,i2-2,i3)-8.*sy(i1,i2-1,i3)+8.*sy(i1,i2+1,i3)-sy(i1,i2+2,i3))/(12.*dr(1))
                             syrs     = ((sy(i1-2,i2-2,i3)-8.*sy(i1-2,i2-1,i3)+8.*sy(i1-2,i2+1,i3)-sy(i1-2,i2+2,i3))/(12.*dr(1))-8.*(sy(i1-1,i2-2,i3)-8.*sy(i1-1,i2-1,i3)+8.*sy(i1-1,i2+1,i3)-sy(i1-1,i2+2,i3))/(12.*dr(1))+8.*(sy(i1+1,i2-2,i3)-8.*sy(i1+1,i2-1,i3)+8.*sy(i1+1,i2+1,i3)-sy(i1+1,i2+2,i3))/(12.*dr(1))-(sy(i1+2,i2-2,i3)-8.*sy(i1+2,i2-1,i3)+8.*sy(i1+2,i2+1,i3)-sy(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             syss     = (-sy(i1,i2-2,i3)+16.*sy(i1,i2-1,i3)-30.*sy(i1,i2,i3)+16.*sy(i1,i2+1,i3)-sy(i1,i2+2,i3))/(12.*dr(1)**2)
                             ! ---------- Spatial derivatives of metrics rx, sx, ry, ... ---------
                             rxi = rx(i1,i2,i3)
                             ryi = ry(i1,i2,i3)
                             sxi = sx(i1,i2,i3)
                             syi = sy(i1,i2,i3)
                             rxx      = rxi*rxr+sxi*rxs
                             rxy      = ryi*rxr+syi*rxs
                             ryy      = ryi*ryr+syi*rys
                             sxx      = rxi*sxr+sxi*sxs
                             sxy      = ryi*sxr+syi*sxs
                             syy      = ryi*syr+syi*sys
                             rxxx     = rxi**2*rxrr+2.*rxi*sxi*rxrs+sxi**2*rxss+rxx*rxr+sxx*rxs
                             rxxy     = rxi*ryi*rxrr+(rxi*syi+ryi*sxi)*rxrs+sxi*syi*rxss+rxy*rxr+sxy*rxs
                             rxyy     = ryi**2*rxrr+2.*ryi*syi*rxrs+syi**2*rxss+ryy*rxr+syy*rxs
                             ryyy     = ryi**2*ryrr+2.*ryi*syi*ryrs+syi**2*ryss+ryy*ryr+syy*rys
                             sxxx     = rxi**2*sxrr+2.*rxi*sxi*sxrs+sxi**2*sxss+rxx*sxr+sxx*sxs
                             sxxy     = rxi*ryi*sxrr+(rxi*syi+ryi*sxi)*sxrs+sxi*syi*sxss+rxy*sxr+sxy*sxs
                             sxyy     = ryi**2*sxrr+2.*ryi*syi*sxrs+syi**2*sxss+ryy*sxr+syy*sxs
                             syyy     = ryi**2*syrr+2.*ryi*syi*syrs+syi**2*syss+ryy*syr+syy*sys
                             ! ---- end evalMetrics eq evalMetrics ---
                             ! ---------- Third spatial derivatives of u ---------
                             ! ---------- END CURVILINEAR  ---------
                             ! ------ curvilinear grid: -------
                               ! eval equation with wrong values at ghost:
                               r1 =  an1*ux42(i1,i2,i3,0) + an2*uy42(i1,i2,i3,0) - gg
                               crv(axis) = an1*rsxy(i1,i2,i3,axis,0) + an2*rsxy(i1,i2,i3,axis,1)
                               ! crv(1) = an1*rsxy(i1,i2,i3,1,0) + an2*rsxy(i1,i2,i3,1,1)
                               ! **CHECK ME**
                               a11 = -is*( crv(axis)*8./(12.*dr(axis)) )  ! coeff of u(-1)
                               a12 =  is*( crv(axis)*1./(12.*dr(axis)) )  ! coeff of u(-2)
                               r2 = c2*( an1*( uxxx + uxyy ) + an2*( uxxy + uyyy ) ) + nDotGradF - gtt
                   ! uxxx = rxi**3*urrr+3.*rxi**2*sxi*urrs+3.*rxi*sxi**2*urss+sxi**3*usss+3.*rxi*rxx*urr+(3.*rxi*sxx+3.*rxx*sxi)*urs+3.*sxx*sxi*uss+rxxx*ur+sxxx*us
                   ! uxyy = ryi**2*rxi*urrr+(syi*ryi*rxi+ryi*(rxi*syi+ryi*sxi))*urrs+(ryi*syi*sxi+syi*(rxi*syi+ryi*sxi))*urss+syi**2*sxi*usss+(rxi*ryy+2.*rxy*ryi)*urr+(rxi*syy+2.*rxy*syi+2.*ryi*sxy+ryy*sxi)*urs+(sxi*syy+2.*sxy*syi)*uss+rxyy*ur+sxyy*us
                   ! uxxy = ryi*rxi**2*urrr+(rxi**2*syi+2.*rxi*ryi*sxi)*urrs+(2.*rxi*sxi*syi+ryi*sxi**2)*urss+syi*sxi**2*usss+(2.*rxi*rxy+rxx*ryi)*urr+(2.*rxi*sxy+rxx*syi+2.*rxy*sxi+ryi*sxx)*urs+(2.*sxi*sxy+sxx*syi)*uss+rxxy*ur+sxxy*us
                   ! uyyy = ryi**3*urrr+3.*ryi**2*syi*urrs+3.*ryi*syi**2*urss+syi**3*usss+3.*ryi*ryy*urr+(3.*ryi*syy+3.*ryy*syi)*urs+3.*syy*syi*uss+ryyy*ur+syyy*us
                               ! uxxx = rxi**3*urrr +3.*rxi*rxx*urr + rxxx*ur
                               ! uxyy = ryi**2*rxi*urrr +(rxi*ryy+2.*rxy*ryi)*urr+rxyy*ur
                               ! uxxy = ryi*rxi**2*urrr +(2.*rxi*rxy+rxx*ryi)*urr+rxxy*ur
                               ! Coeff of terms involving urrr (do not include terms involving urr and ur as these used 2nd-order values)
                               crv(axis) = (an1*rsxy(i1,i2,i3,axis,0) + an2*rsxy(i1,i2,i3,axis,1) )*( rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2 )
                               ! if( axis.eq.0 )then
                               !   crv(axis) = (an1*rxi + an2*ryi )*( rxi**2 + ryi**2 )
                               !   ! crv(axis) = an1*( rxi**3 + ryi**2*rxi ) + !   !             an2*( ryi**3 + ryi*rxi**2 )
                               ! else
                               !   crv(axis) = (an1*sxi + an2*syi )*( sxi**2 + syi**2 )
                               !   ! crv(axis) = an1*( sxi**3 + syi**2*sxi ) + !   !             an2*( syi**3 + syi*sxi**2 )              
                               ! end if
                               ! crv(1) = an1*rsxy(i1,i2,i3,1,0)**3 + an2*rsxy(i1,i2,i3,1,1)**3
                               ! **CHECK ME**
                               a21 =  is*c2*( 2.*crv(axis)/(2.*dr(axis)**3) )
                               a22 = -is*c2*( 1.*crv(axis)/(2.*dr(axis)**3) )
                                 if( checkCoeff.eq.1 )then
                                   ! --- Check coefficients a11,a12,... by discrete delta ---
                                   u1Save = u(j1,j2,j3,0)
                                   u2Save = u(k1,k2,k3,0)
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=0.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   ! ---------- END CURVILINEAR  ---------
                                   r1a = (an1*ux42(i1,i2,i3,0)+an2*uy42(i1,i2,i3,0))
                                   r2a = (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   u(j1,j2,j3,0)=1.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   ! ---------- END CURVILINEAR  ---------
                                   r1b =  (an1*ux42(i1,i2,i3,0)+an2*uy42(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   a11c = r1b - r1a
                                   a21c = r2b - r2a
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=1.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   ! ---------- END CURVILINEAR  ---------
                                   r1b =  (an1*ux42(i1,i2,i3,0)+an2*uy42(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   a12c = r1b - r1a
                                   a22c = r2b - r2a
                                   u(j1,j2,j3,0) = u1Save  
                                   u(k1,k2,k3,0) = u2Save      
                                   maxDiff=max(maxDiff,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c),abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c))
                                   ! #If "curvilinear" eq "curvilinear"
                                   ! write(*,'("++ i1,i2=",2i4," a11,a11c=",2(1pe12.4)," a12,a12c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a11,a11c,a12,a12c,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c)
                                   ! write(*,'("++ i1,i2=",2i4," a21,a21c=",2(1pe12.4)," a22,a22c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a21,a21c,a22,a22c,abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c)
                                   ! #End
                                 end if 
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2 )/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1 )/det
                           if( .false. )then
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,uex )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,uey )
                             write(*,'(" (i1,i2)=(",2i3,") u(-1),ue(-1)=",2(1pe12.4)," u(-2),ue(-2)=",2(1pe12.4))') i1,i2,uTemp(j1,j2,j3,0),uex,uTemp(k1,k2,k3,0),uey
                           end if
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! this case done below 
                         end if
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
                        end do
                        end do
                        end do
                       if( checkCoeff.eq.1 )then
                         write(*,'("bcOptWave: t=",e12.3," checkCoeff: (side,axis)=(",2i2,") rel-maxDiff=",e9.2," curvilinear 4")') t,side,axis,maxDiff
                       end if
                 else
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       ! ------------------- 4 4 NEUMANN CBC --------------------
                       maxDiff=0. ! for checkCoeff
                       gg=0.
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         ghost =2
                         k1=i1-is1*ghost
                         k2=i2-is2*ghost
                         k3=i3-is3*ghost   
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
                               ! No forcing, do nothing 
                             ! evaluate 3rd derivatives : uxxx,uxxy,uxxz,uxyy,uxzz, uyyy,uyyz,uyzz, uzzz 
                             ! ---------- START CURVILINEAR  ---------
                             ! ---------- Parametric derivatives ---------
                             ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                             urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                             urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                             us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                             urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                             urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                             uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                             urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                             usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                             ut       = (u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))
                             urt      = ((u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0))
                             urrt     = (-(u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0)**2)
                             ust      = ((u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))
                             urst     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                             usst     = (-(u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)**2)
                             utt      = (-u(i1,i2,i3-2,0)+16.*u(i1,i2,i3-1,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2)**2)
                             urtt     = ((-u(i1-2,i2,i3-2,0)+16.*u(i1-2,i2,i3-1,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1-1,i2,i3-2,0)+16.*u(i1-1,i2,i3-1,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1+1,i2,i3-2,0)+16.*u(i1+1,i2,i3-1,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2)**2)-(-u(i1+2,i2,i3-2,0)+16.*u(i1+2,i2,i3-1,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(0))
                             ustt     = ((-u(i1,i2-2,i3-2,0)+16.*u(i1,i2-2,i3-1,0)-30.*u(i1,i2-2,i3,0)+16.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1,i2-1,i3-2,0)+16.*u(i1,i2-1,i3-1,0)-30.*u(i1,i2-1,i3,0)+16.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1,i2+1,i3-2,0)+16.*u(i1,i2+1,i3-1,0)-30.*u(i1,i2+1,i3,0)+16.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2)**2)-(-u(i1,i2+2,i3-2,0)+16.*u(i1,i2+2,i3-1,0)-30.*u(i1,i2+2,i3,0)+16.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(1))
                             uttt     = (-u(i1,i2,i3-2,0)+2.*u(i1,i2,i3-1,0)-2.*u(i1,i2,i3+1,0)+u(i1,i2,i3+2,0))/(2.*dr(2)**3)
                             rxr      = (rx(i1-2,i2,i3)-8.*rx(i1-1,i2,i3)+8.*rx(i1+1,i2,i3)-rx(i1+2,i2,i3))/(12.*dr(0))
                             rxrr     = (-rx(i1-2,i2,i3)+16.*rx(i1-1,i2,i3)-30.*rx(i1,i2,i3)+16.*rx(i1+1,i2,i3)-rx(i1+2,i2,i3))/(12.*dr(0)**2)
                             rxs      = (rx(i1,i2-2,i3)-8.*rx(i1,i2-1,i3)+8.*rx(i1,i2+1,i3)-rx(i1,i2+2,i3))/(12.*dr(1))
                             rxrs     = ((rx(i1-2,i2-2,i3)-8.*rx(i1-2,i2-1,i3)+8.*rx(i1-2,i2+1,i3)-rx(i1-2,i2+2,i3))/(12.*dr(1))-8.*(rx(i1-1,i2-2,i3)-8.*rx(i1-1,i2-1,i3)+8.*rx(i1-1,i2+1,i3)-rx(i1-1,i2+2,i3))/(12.*dr(1))+8.*(rx(i1+1,i2-2,i3)-8.*rx(i1+1,i2-1,i3)+8.*rx(i1+1,i2+1,i3)-rx(i1+1,i2+2,i3))/(12.*dr(1))-(rx(i1+2,i2-2,i3)-8.*rx(i1+2,i2-1,i3)+8.*rx(i1+2,i2+1,i3)-rx(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             rxss     = (-rx(i1,i2-2,i3)+16.*rx(i1,i2-1,i3)-30.*rx(i1,i2,i3)+16.*rx(i1,i2+1,i3)-rx(i1,i2+2,i3))/(12.*dr(1)**2)
                             rxt      = (rx(i1,i2,i3-2)-8.*rx(i1,i2,i3-1)+8.*rx(i1,i2,i3+1)-rx(i1,i2,i3+2))/(12.*dr(2))
                             rxrt     = ((rx(i1-2,i2,i3-2)-8.*rx(i1-2,i2,i3-1)+8.*rx(i1-2,i2,i3+1)-rx(i1-2,i2,i3+2))/(12.*dr(2))-8.*(rx(i1-1,i2,i3-2)-8.*rx(i1-1,i2,i3-1)+8.*rx(i1-1,i2,i3+1)-rx(i1-1,i2,i3+2))/(12.*dr(2))+8.*(rx(i1+1,i2,i3-2)-8.*rx(i1+1,i2,i3-1)+8.*rx(i1+1,i2,i3+1)-rx(i1+1,i2,i3+2))/(12.*dr(2))-(rx(i1+2,i2,i3-2)-8.*rx(i1+2,i2,i3-1)+8.*rx(i1+2,i2,i3+1)-rx(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             rxst     = ((rx(i1,i2-2,i3-2)-8.*rx(i1,i2-2,i3-1)+8.*rx(i1,i2-2,i3+1)-rx(i1,i2-2,i3+2))/(12.*dr(2))-8.*(rx(i1,i2-1,i3-2)-8.*rx(i1,i2-1,i3-1)+8.*rx(i1,i2-1,i3+1)-rx(i1,i2-1,i3+2))/(12.*dr(2))+8.*(rx(i1,i2+1,i3-2)-8.*rx(i1,i2+1,i3-1)+8.*rx(i1,i2+1,i3+1)-rx(i1,i2+1,i3+2))/(12.*dr(2))-(rx(i1,i2+2,i3-2)-8.*rx(i1,i2+2,i3-1)+8.*rx(i1,i2+2,i3+1)-rx(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             rxtt     = (-rx(i1,i2,i3-2)+16.*rx(i1,i2,i3-1)-30.*rx(i1,i2,i3)+16.*rx(i1,i2,i3+1)-rx(i1,i2,i3+2))/(12.*dr(2)**2)
                             ryr      = (ry(i1-2,i2,i3)-8.*ry(i1-1,i2,i3)+8.*ry(i1+1,i2,i3)-ry(i1+2,i2,i3))/(12.*dr(0))
                             ryrr     = (-ry(i1-2,i2,i3)+16.*ry(i1-1,i2,i3)-30.*ry(i1,i2,i3)+16.*ry(i1+1,i2,i3)-ry(i1+2,i2,i3))/(12.*dr(0)**2)
                             rys      = (ry(i1,i2-2,i3)-8.*ry(i1,i2-1,i3)+8.*ry(i1,i2+1,i3)-ry(i1,i2+2,i3))/(12.*dr(1))
                             ryrs     = ((ry(i1-2,i2-2,i3)-8.*ry(i1-2,i2-1,i3)+8.*ry(i1-2,i2+1,i3)-ry(i1-2,i2+2,i3))/(12.*dr(1))-8.*(ry(i1-1,i2-2,i3)-8.*ry(i1-1,i2-1,i3)+8.*ry(i1-1,i2+1,i3)-ry(i1-1,i2+2,i3))/(12.*dr(1))+8.*(ry(i1+1,i2-2,i3)-8.*ry(i1+1,i2-1,i3)+8.*ry(i1+1,i2+1,i3)-ry(i1+1,i2+2,i3))/(12.*dr(1))-(ry(i1+2,i2-2,i3)-8.*ry(i1+2,i2-1,i3)+8.*ry(i1+2,i2+1,i3)-ry(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             ryss     = (-ry(i1,i2-2,i3)+16.*ry(i1,i2-1,i3)-30.*ry(i1,i2,i3)+16.*ry(i1,i2+1,i3)-ry(i1,i2+2,i3))/(12.*dr(1)**2)
                             ryt      = (ry(i1,i2,i3-2)-8.*ry(i1,i2,i3-1)+8.*ry(i1,i2,i3+1)-ry(i1,i2,i3+2))/(12.*dr(2))
                             ryrt     = ((ry(i1-2,i2,i3-2)-8.*ry(i1-2,i2,i3-1)+8.*ry(i1-2,i2,i3+1)-ry(i1-2,i2,i3+2))/(12.*dr(2))-8.*(ry(i1-1,i2,i3-2)-8.*ry(i1-1,i2,i3-1)+8.*ry(i1-1,i2,i3+1)-ry(i1-1,i2,i3+2))/(12.*dr(2))+8.*(ry(i1+1,i2,i3-2)-8.*ry(i1+1,i2,i3-1)+8.*ry(i1+1,i2,i3+1)-ry(i1+1,i2,i3+2))/(12.*dr(2))-(ry(i1+2,i2,i3-2)-8.*ry(i1+2,i2,i3-1)+8.*ry(i1+2,i2,i3+1)-ry(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             ryst     = ((ry(i1,i2-2,i3-2)-8.*ry(i1,i2-2,i3-1)+8.*ry(i1,i2-2,i3+1)-ry(i1,i2-2,i3+2))/(12.*dr(2))-8.*(ry(i1,i2-1,i3-2)-8.*ry(i1,i2-1,i3-1)+8.*ry(i1,i2-1,i3+1)-ry(i1,i2-1,i3+2))/(12.*dr(2))+8.*(ry(i1,i2+1,i3-2)-8.*ry(i1,i2+1,i3-1)+8.*ry(i1,i2+1,i3+1)-ry(i1,i2+1,i3+2))/(12.*dr(2))-(ry(i1,i2+2,i3-2)-8.*ry(i1,i2+2,i3-1)+8.*ry(i1,i2+2,i3+1)-ry(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             rytt     = (-ry(i1,i2,i3-2)+16.*ry(i1,i2,i3-1)-30.*ry(i1,i2,i3)+16.*ry(i1,i2,i3+1)-ry(i1,i2,i3+2))/(12.*dr(2)**2)
                             sxr      = (sx(i1-2,i2,i3)-8.*sx(i1-1,i2,i3)+8.*sx(i1+1,i2,i3)-sx(i1+2,i2,i3))/(12.*dr(0))
                             sxrr     = (-sx(i1-2,i2,i3)+16.*sx(i1-1,i2,i3)-30.*sx(i1,i2,i3)+16.*sx(i1+1,i2,i3)-sx(i1+2,i2,i3))/(12.*dr(0)**2)
                             sxs      = (sx(i1,i2-2,i3)-8.*sx(i1,i2-1,i3)+8.*sx(i1,i2+1,i3)-sx(i1,i2+2,i3))/(12.*dr(1))
                             sxrs     = ((sx(i1-2,i2-2,i3)-8.*sx(i1-2,i2-1,i3)+8.*sx(i1-2,i2+1,i3)-sx(i1-2,i2+2,i3))/(12.*dr(1))-8.*(sx(i1-1,i2-2,i3)-8.*sx(i1-1,i2-1,i3)+8.*sx(i1-1,i2+1,i3)-sx(i1-1,i2+2,i3))/(12.*dr(1))+8.*(sx(i1+1,i2-2,i3)-8.*sx(i1+1,i2-1,i3)+8.*sx(i1+1,i2+1,i3)-sx(i1+1,i2+2,i3))/(12.*dr(1))-(sx(i1+2,i2-2,i3)-8.*sx(i1+2,i2-1,i3)+8.*sx(i1+2,i2+1,i3)-sx(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             sxss     = (-sx(i1,i2-2,i3)+16.*sx(i1,i2-1,i3)-30.*sx(i1,i2,i3)+16.*sx(i1,i2+1,i3)-sx(i1,i2+2,i3))/(12.*dr(1)**2)
                             sxt      = (sx(i1,i2,i3-2)-8.*sx(i1,i2,i3-1)+8.*sx(i1,i2,i3+1)-sx(i1,i2,i3+2))/(12.*dr(2))
                             sxrt     = ((sx(i1-2,i2,i3-2)-8.*sx(i1-2,i2,i3-1)+8.*sx(i1-2,i2,i3+1)-sx(i1-2,i2,i3+2))/(12.*dr(2))-8.*(sx(i1-1,i2,i3-2)-8.*sx(i1-1,i2,i3-1)+8.*sx(i1-1,i2,i3+1)-sx(i1-1,i2,i3+2))/(12.*dr(2))+8.*(sx(i1+1,i2,i3-2)-8.*sx(i1+1,i2,i3-1)+8.*sx(i1+1,i2,i3+1)-sx(i1+1,i2,i3+2))/(12.*dr(2))-(sx(i1+2,i2,i3-2)-8.*sx(i1+2,i2,i3-1)+8.*sx(i1+2,i2,i3+1)-sx(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             sxst     = ((sx(i1,i2-2,i3-2)-8.*sx(i1,i2-2,i3-1)+8.*sx(i1,i2-2,i3+1)-sx(i1,i2-2,i3+2))/(12.*dr(2))-8.*(sx(i1,i2-1,i3-2)-8.*sx(i1,i2-1,i3-1)+8.*sx(i1,i2-1,i3+1)-sx(i1,i2-1,i3+2))/(12.*dr(2))+8.*(sx(i1,i2+1,i3-2)-8.*sx(i1,i2+1,i3-1)+8.*sx(i1,i2+1,i3+1)-sx(i1,i2+1,i3+2))/(12.*dr(2))-(sx(i1,i2+2,i3-2)-8.*sx(i1,i2+2,i3-1)+8.*sx(i1,i2+2,i3+1)-sx(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             sxtt     = (-sx(i1,i2,i3-2)+16.*sx(i1,i2,i3-1)-30.*sx(i1,i2,i3)+16.*sx(i1,i2,i3+1)-sx(i1,i2,i3+2))/(12.*dr(2)**2)
                             syr      = (sy(i1-2,i2,i3)-8.*sy(i1-1,i2,i3)+8.*sy(i1+1,i2,i3)-sy(i1+2,i2,i3))/(12.*dr(0))
                             syrr     = (-sy(i1-2,i2,i3)+16.*sy(i1-1,i2,i3)-30.*sy(i1,i2,i3)+16.*sy(i1+1,i2,i3)-sy(i1+2,i2,i3))/(12.*dr(0)**2)
                             sys      = (sy(i1,i2-2,i3)-8.*sy(i1,i2-1,i3)+8.*sy(i1,i2+1,i3)-sy(i1,i2+2,i3))/(12.*dr(1))
                             syrs     = ((sy(i1-2,i2-2,i3)-8.*sy(i1-2,i2-1,i3)+8.*sy(i1-2,i2+1,i3)-sy(i1-2,i2+2,i3))/(12.*dr(1))-8.*(sy(i1-1,i2-2,i3)-8.*sy(i1-1,i2-1,i3)+8.*sy(i1-1,i2+1,i3)-sy(i1-1,i2+2,i3))/(12.*dr(1))+8.*(sy(i1+1,i2-2,i3)-8.*sy(i1+1,i2-1,i3)+8.*sy(i1+1,i2+1,i3)-sy(i1+1,i2+2,i3))/(12.*dr(1))-(sy(i1+2,i2-2,i3)-8.*sy(i1+2,i2-1,i3)+8.*sy(i1+2,i2+1,i3)-sy(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             syss     = (-sy(i1,i2-2,i3)+16.*sy(i1,i2-1,i3)-30.*sy(i1,i2,i3)+16.*sy(i1,i2+1,i3)-sy(i1,i2+2,i3))/(12.*dr(1)**2)
                             syt      = (sy(i1,i2,i3-2)-8.*sy(i1,i2,i3-1)+8.*sy(i1,i2,i3+1)-sy(i1,i2,i3+2))/(12.*dr(2))
                             syrt     = ((sy(i1-2,i2,i3-2)-8.*sy(i1-2,i2,i3-1)+8.*sy(i1-2,i2,i3+1)-sy(i1-2,i2,i3+2))/(12.*dr(2))-8.*(sy(i1-1,i2,i3-2)-8.*sy(i1-1,i2,i3-1)+8.*sy(i1-1,i2,i3+1)-sy(i1-1,i2,i3+2))/(12.*dr(2))+8.*(sy(i1+1,i2,i3-2)-8.*sy(i1+1,i2,i3-1)+8.*sy(i1+1,i2,i3+1)-sy(i1+1,i2,i3+2))/(12.*dr(2))-(sy(i1+2,i2,i3-2)-8.*sy(i1+2,i2,i3-1)+8.*sy(i1+2,i2,i3+1)-sy(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             syst     = ((sy(i1,i2-2,i3-2)-8.*sy(i1,i2-2,i3-1)+8.*sy(i1,i2-2,i3+1)-sy(i1,i2-2,i3+2))/(12.*dr(2))-8.*(sy(i1,i2-1,i3-2)-8.*sy(i1,i2-1,i3-1)+8.*sy(i1,i2-1,i3+1)-sy(i1,i2-1,i3+2))/(12.*dr(2))+8.*(sy(i1,i2+1,i3-2)-8.*sy(i1,i2+1,i3-1)+8.*sy(i1,i2+1,i3+1)-sy(i1,i2+1,i3+2))/(12.*dr(2))-(sy(i1,i2+2,i3-2)-8.*sy(i1,i2+2,i3-1)+8.*sy(i1,i2+2,i3+1)-sy(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             sytt     = (-sy(i1,i2,i3-2)+16.*sy(i1,i2,i3-1)-30.*sy(i1,i2,i3)+16.*sy(i1,i2,i3+1)-sy(i1,i2,i3+2))/(12.*dr(2)**2)
                             rzr      = (rz(i1-2,i2,i3)-8.*rz(i1-1,i2,i3)+8.*rz(i1+1,i2,i3)-rz(i1+2,i2,i3))/(12.*dr(0))
                             rzrr     = (-rz(i1-2,i2,i3)+16.*rz(i1-1,i2,i3)-30.*rz(i1,i2,i3)+16.*rz(i1+1,i2,i3)-rz(i1+2,i2,i3))/(12.*dr(0)**2)
                             rzs      = (rz(i1,i2-2,i3)-8.*rz(i1,i2-1,i3)+8.*rz(i1,i2+1,i3)-rz(i1,i2+2,i3))/(12.*dr(1))
                             rzrs     = ((rz(i1-2,i2-2,i3)-8.*rz(i1-2,i2-1,i3)+8.*rz(i1-2,i2+1,i3)-rz(i1-2,i2+2,i3))/(12.*dr(1))-8.*(rz(i1-1,i2-2,i3)-8.*rz(i1-1,i2-1,i3)+8.*rz(i1-1,i2+1,i3)-rz(i1-1,i2+2,i3))/(12.*dr(1))+8.*(rz(i1+1,i2-2,i3)-8.*rz(i1+1,i2-1,i3)+8.*rz(i1+1,i2+1,i3)-rz(i1+1,i2+2,i3))/(12.*dr(1))-(rz(i1+2,i2-2,i3)-8.*rz(i1+2,i2-1,i3)+8.*rz(i1+2,i2+1,i3)-rz(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             rzss     = (-rz(i1,i2-2,i3)+16.*rz(i1,i2-1,i3)-30.*rz(i1,i2,i3)+16.*rz(i1,i2+1,i3)-rz(i1,i2+2,i3))/(12.*dr(1)**2)
                             rzt      = (rz(i1,i2,i3-2)-8.*rz(i1,i2,i3-1)+8.*rz(i1,i2,i3+1)-rz(i1,i2,i3+2))/(12.*dr(2))
                             rzrt     = ((rz(i1-2,i2,i3-2)-8.*rz(i1-2,i2,i3-1)+8.*rz(i1-2,i2,i3+1)-rz(i1-2,i2,i3+2))/(12.*dr(2))-8.*(rz(i1-1,i2,i3-2)-8.*rz(i1-1,i2,i3-1)+8.*rz(i1-1,i2,i3+1)-rz(i1-1,i2,i3+2))/(12.*dr(2))+8.*(rz(i1+1,i2,i3-2)-8.*rz(i1+1,i2,i3-1)+8.*rz(i1+1,i2,i3+1)-rz(i1+1,i2,i3+2))/(12.*dr(2))-(rz(i1+2,i2,i3-2)-8.*rz(i1+2,i2,i3-1)+8.*rz(i1+2,i2,i3+1)-rz(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             rzst     = ((rz(i1,i2-2,i3-2)-8.*rz(i1,i2-2,i3-1)+8.*rz(i1,i2-2,i3+1)-rz(i1,i2-2,i3+2))/(12.*dr(2))-8.*(rz(i1,i2-1,i3-2)-8.*rz(i1,i2-1,i3-1)+8.*rz(i1,i2-1,i3+1)-rz(i1,i2-1,i3+2))/(12.*dr(2))+8.*(rz(i1,i2+1,i3-2)-8.*rz(i1,i2+1,i3-1)+8.*rz(i1,i2+1,i3+1)-rz(i1,i2+1,i3+2))/(12.*dr(2))-(rz(i1,i2+2,i3-2)-8.*rz(i1,i2+2,i3-1)+8.*rz(i1,i2+2,i3+1)-rz(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             rztt     = (-rz(i1,i2,i3-2)+16.*rz(i1,i2,i3-1)-30.*rz(i1,i2,i3)+16.*rz(i1,i2,i3+1)-rz(i1,i2,i3+2))/(12.*dr(2)**2)
                             szr      = (sz(i1-2,i2,i3)-8.*sz(i1-1,i2,i3)+8.*sz(i1+1,i2,i3)-sz(i1+2,i2,i3))/(12.*dr(0))
                             szrr     = (-sz(i1-2,i2,i3)+16.*sz(i1-1,i2,i3)-30.*sz(i1,i2,i3)+16.*sz(i1+1,i2,i3)-sz(i1+2,i2,i3))/(12.*dr(0)**2)
                             szs      = (sz(i1,i2-2,i3)-8.*sz(i1,i2-1,i3)+8.*sz(i1,i2+1,i3)-sz(i1,i2+2,i3))/(12.*dr(1))
                             szrs     = ((sz(i1-2,i2-2,i3)-8.*sz(i1-2,i2-1,i3)+8.*sz(i1-2,i2+1,i3)-sz(i1-2,i2+2,i3))/(12.*dr(1))-8.*(sz(i1-1,i2-2,i3)-8.*sz(i1-1,i2-1,i3)+8.*sz(i1-1,i2+1,i3)-sz(i1-1,i2+2,i3))/(12.*dr(1))+8.*(sz(i1+1,i2-2,i3)-8.*sz(i1+1,i2-1,i3)+8.*sz(i1+1,i2+1,i3)-sz(i1+1,i2+2,i3))/(12.*dr(1))-(sz(i1+2,i2-2,i3)-8.*sz(i1+2,i2-1,i3)+8.*sz(i1+2,i2+1,i3)-sz(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             szss     = (-sz(i1,i2-2,i3)+16.*sz(i1,i2-1,i3)-30.*sz(i1,i2,i3)+16.*sz(i1,i2+1,i3)-sz(i1,i2+2,i3))/(12.*dr(1)**2)
                             szt      = (sz(i1,i2,i3-2)-8.*sz(i1,i2,i3-1)+8.*sz(i1,i2,i3+1)-sz(i1,i2,i3+2))/(12.*dr(2))
                             szrt     = ((sz(i1-2,i2,i3-2)-8.*sz(i1-2,i2,i3-1)+8.*sz(i1-2,i2,i3+1)-sz(i1-2,i2,i3+2))/(12.*dr(2))-8.*(sz(i1-1,i2,i3-2)-8.*sz(i1-1,i2,i3-1)+8.*sz(i1-1,i2,i3+1)-sz(i1-1,i2,i3+2))/(12.*dr(2))+8.*(sz(i1+1,i2,i3-2)-8.*sz(i1+1,i2,i3-1)+8.*sz(i1+1,i2,i3+1)-sz(i1+1,i2,i3+2))/(12.*dr(2))-(sz(i1+2,i2,i3-2)-8.*sz(i1+2,i2,i3-1)+8.*sz(i1+2,i2,i3+1)-sz(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             szst     = ((sz(i1,i2-2,i3-2)-8.*sz(i1,i2-2,i3-1)+8.*sz(i1,i2-2,i3+1)-sz(i1,i2-2,i3+2))/(12.*dr(2))-8.*(sz(i1,i2-1,i3-2)-8.*sz(i1,i2-1,i3-1)+8.*sz(i1,i2-1,i3+1)-sz(i1,i2-1,i3+2))/(12.*dr(2))+8.*(sz(i1,i2+1,i3-2)-8.*sz(i1,i2+1,i3-1)+8.*sz(i1,i2+1,i3+1)-sz(i1,i2+1,i3+2))/(12.*dr(2))-(sz(i1,i2+2,i3-2)-8.*sz(i1,i2+2,i3-1)+8.*sz(i1,i2+2,i3+1)-sz(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             sztt     = (-sz(i1,i2,i3-2)+16.*sz(i1,i2,i3-1)-30.*sz(i1,i2,i3)+16.*sz(i1,i2,i3+1)-sz(i1,i2,i3+2))/(12.*dr(2)**2)
                             txr      = (tx(i1-2,i2,i3)-8.*tx(i1-1,i2,i3)+8.*tx(i1+1,i2,i3)-tx(i1+2,i2,i3))/(12.*dr(0))
                             txrr     = (-tx(i1-2,i2,i3)+16.*tx(i1-1,i2,i3)-30.*tx(i1,i2,i3)+16.*tx(i1+1,i2,i3)-tx(i1+2,i2,i3))/(12.*dr(0)**2)
                             txs      = (tx(i1,i2-2,i3)-8.*tx(i1,i2-1,i3)+8.*tx(i1,i2+1,i3)-tx(i1,i2+2,i3))/(12.*dr(1))
                             txrs     = ((tx(i1-2,i2-2,i3)-8.*tx(i1-2,i2-1,i3)+8.*tx(i1-2,i2+1,i3)-tx(i1-2,i2+2,i3))/(12.*dr(1))-8.*(tx(i1-1,i2-2,i3)-8.*tx(i1-1,i2-1,i3)+8.*tx(i1-1,i2+1,i3)-tx(i1-1,i2+2,i3))/(12.*dr(1))+8.*(tx(i1+1,i2-2,i3)-8.*tx(i1+1,i2-1,i3)+8.*tx(i1+1,i2+1,i3)-tx(i1+1,i2+2,i3))/(12.*dr(1))-(tx(i1+2,i2-2,i3)-8.*tx(i1+2,i2-1,i3)+8.*tx(i1+2,i2+1,i3)-tx(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             txss     = (-tx(i1,i2-2,i3)+16.*tx(i1,i2-1,i3)-30.*tx(i1,i2,i3)+16.*tx(i1,i2+1,i3)-tx(i1,i2+2,i3))/(12.*dr(1)**2)
                             txt      = (tx(i1,i2,i3-2)-8.*tx(i1,i2,i3-1)+8.*tx(i1,i2,i3+1)-tx(i1,i2,i3+2))/(12.*dr(2))
                             txrt     = ((tx(i1-2,i2,i3-2)-8.*tx(i1-2,i2,i3-1)+8.*tx(i1-2,i2,i3+1)-tx(i1-2,i2,i3+2))/(12.*dr(2))-8.*(tx(i1-1,i2,i3-2)-8.*tx(i1-1,i2,i3-1)+8.*tx(i1-1,i2,i3+1)-tx(i1-1,i2,i3+2))/(12.*dr(2))+8.*(tx(i1+1,i2,i3-2)-8.*tx(i1+1,i2,i3-1)+8.*tx(i1+1,i2,i3+1)-tx(i1+1,i2,i3+2))/(12.*dr(2))-(tx(i1+2,i2,i3-2)-8.*tx(i1+2,i2,i3-1)+8.*tx(i1+2,i2,i3+1)-tx(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             txst     = ((tx(i1,i2-2,i3-2)-8.*tx(i1,i2-2,i3-1)+8.*tx(i1,i2-2,i3+1)-tx(i1,i2-2,i3+2))/(12.*dr(2))-8.*(tx(i1,i2-1,i3-2)-8.*tx(i1,i2-1,i3-1)+8.*tx(i1,i2-1,i3+1)-tx(i1,i2-1,i3+2))/(12.*dr(2))+8.*(tx(i1,i2+1,i3-2)-8.*tx(i1,i2+1,i3-1)+8.*tx(i1,i2+1,i3+1)-tx(i1,i2+1,i3+2))/(12.*dr(2))-(tx(i1,i2+2,i3-2)-8.*tx(i1,i2+2,i3-1)+8.*tx(i1,i2+2,i3+1)-tx(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             txtt     = (-tx(i1,i2,i3-2)+16.*tx(i1,i2,i3-1)-30.*tx(i1,i2,i3)+16.*tx(i1,i2,i3+1)-tx(i1,i2,i3+2))/(12.*dr(2)**2)
                             tyr      = (ty(i1-2,i2,i3)-8.*ty(i1-1,i2,i3)+8.*ty(i1+1,i2,i3)-ty(i1+2,i2,i3))/(12.*dr(0))
                             tyrr     = (-ty(i1-2,i2,i3)+16.*ty(i1-1,i2,i3)-30.*ty(i1,i2,i3)+16.*ty(i1+1,i2,i3)-ty(i1+2,i2,i3))/(12.*dr(0)**2)
                             tys      = (ty(i1,i2-2,i3)-8.*ty(i1,i2-1,i3)+8.*ty(i1,i2+1,i3)-ty(i1,i2+2,i3))/(12.*dr(1))
                             tyrs     = ((ty(i1-2,i2-2,i3)-8.*ty(i1-2,i2-1,i3)+8.*ty(i1-2,i2+1,i3)-ty(i1-2,i2+2,i3))/(12.*dr(1))-8.*(ty(i1-1,i2-2,i3)-8.*ty(i1-1,i2-1,i3)+8.*ty(i1-1,i2+1,i3)-ty(i1-1,i2+2,i3))/(12.*dr(1))+8.*(ty(i1+1,i2-2,i3)-8.*ty(i1+1,i2-1,i3)+8.*ty(i1+1,i2+1,i3)-ty(i1+1,i2+2,i3))/(12.*dr(1))-(ty(i1+2,i2-2,i3)-8.*ty(i1+2,i2-1,i3)+8.*ty(i1+2,i2+1,i3)-ty(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             tyss     = (-ty(i1,i2-2,i3)+16.*ty(i1,i2-1,i3)-30.*ty(i1,i2,i3)+16.*ty(i1,i2+1,i3)-ty(i1,i2+2,i3))/(12.*dr(1)**2)
                             tyt      = (ty(i1,i2,i3-2)-8.*ty(i1,i2,i3-1)+8.*ty(i1,i2,i3+1)-ty(i1,i2,i3+2))/(12.*dr(2))
                             tyrt     = ((ty(i1-2,i2,i3-2)-8.*ty(i1-2,i2,i3-1)+8.*ty(i1-2,i2,i3+1)-ty(i1-2,i2,i3+2))/(12.*dr(2))-8.*(ty(i1-1,i2,i3-2)-8.*ty(i1-1,i2,i3-1)+8.*ty(i1-1,i2,i3+1)-ty(i1-1,i2,i3+2))/(12.*dr(2))+8.*(ty(i1+1,i2,i3-2)-8.*ty(i1+1,i2,i3-1)+8.*ty(i1+1,i2,i3+1)-ty(i1+1,i2,i3+2))/(12.*dr(2))-(ty(i1+2,i2,i3-2)-8.*ty(i1+2,i2,i3-1)+8.*ty(i1+2,i2,i3+1)-ty(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             tyst     = ((ty(i1,i2-2,i3-2)-8.*ty(i1,i2-2,i3-1)+8.*ty(i1,i2-2,i3+1)-ty(i1,i2-2,i3+2))/(12.*dr(2))-8.*(ty(i1,i2-1,i3-2)-8.*ty(i1,i2-1,i3-1)+8.*ty(i1,i2-1,i3+1)-ty(i1,i2-1,i3+2))/(12.*dr(2))+8.*(ty(i1,i2+1,i3-2)-8.*ty(i1,i2+1,i3-1)+8.*ty(i1,i2+1,i3+1)-ty(i1,i2+1,i3+2))/(12.*dr(2))-(ty(i1,i2+2,i3-2)-8.*ty(i1,i2+2,i3-1)+8.*ty(i1,i2+2,i3+1)-ty(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             tytt     = (-ty(i1,i2,i3-2)+16.*ty(i1,i2,i3-1)-30.*ty(i1,i2,i3)+16.*ty(i1,i2,i3+1)-ty(i1,i2,i3+2))/(12.*dr(2)**2)
                             tzr      = (tz(i1-2,i2,i3)-8.*tz(i1-1,i2,i3)+8.*tz(i1+1,i2,i3)-tz(i1+2,i2,i3))/(12.*dr(0))
                             tzrr     = (-tz(i1-2,i2,i3)+16.*tz(i1-1,i2,i3)-30.*tz(i1,i2,i3)+16.*tz(i1+1,i2,i3)-tz(i1+2,i2,i3))/(12.*dr(0)**2)
                             tzs      = (tz(i1,i2-2,i3)-8.*tz(i1,i2-1,i3)+8.*tz(i1,i2+1,i3)-tz(i1,i2+2,i3))/(12.*dr(1))
                             tzrs     = ((tz(i1-2,i2-2,i3)-8.*tz(i1-2,i2-1,i3)+8.*tz(i1-2,i2+1,i3)-tz(i1-2,i2+2,i3))/(12.*dr(1))-8.*(tz(i1-1,i2-2,i3)-8.*tz(i1-1,i2-1,i3)+8.*tz(i1-1,i2+1,i3)-tz(i1-1,i2+2,i3))/(12.*dr(1))+8.*(tz(i1+1,i2-2,i3)-8.*tz(i1+1,i2-1,i3)+8.*tz(i1+1,i2+1,i3)-tz(i1+1,i2+2,i3))/(12.*dr(1))-(tz(i1+2,i2-2,i3)-8.*tz(i1+2,i2-1,i3)+8.*tz(i1+2,i2+1,i3)-tz(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             tzss     = (-tz(i1,i2-2,i3)+16.*tz(i1,i2-1,i3)-30.*tz(i1,i2,i3)+16.*tz(i1,i2+1,i3)-tz(i1,i2+2,i3))/(12.*dr(1)**2)
                             tzt      = (tz(i1,i2,i3-2)-8.*tz(i1,i2,i3-1)+8.*tz(i1,i2,i3+1)-tz(i1,i2,i3+2))/(12.*dr(2))
                             tzrt     = ((tz(i1-2,i2,i3-2)-8.*tz(i1-2,i2,i3-1)+8.*tz(i1-2,i2,i3+1)-tz(i1-2,i2,i3+2))/(12.*dr(2))-8.*(tz(i1-1,i2,i3-2)-8.*tz(i1-1,i2,i3-1)+8.*tz(i1-1,i2,i3+1)-tz(i1-1,i2,i3+2))/(12.*dr(2))+8.*(tz(i1+1,i2,i3-2)-8.*tz(i1+1,i2,i3-1)+8.*tz(i1+1,i2,i3+1)-tz(i1+1,i2,i3+2))/(12.*dr(2))-(tz(i1+2,i2,i3-2)-8.*tz(i1+2,i2,i3-1)+8.*tz(i1+2,i2,i3+1)-tz(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             tzst     = ((tz(i1,i2-2,i3-2)-8.*tz(i1,i2-2,i3-1)+8.*tz(i1,i2-2,i3+1)-tz(i1,i2-2,i3+2))/(12.*dr(2))-8.*(tz(i1,i2-1,i3-2)-8.*tz(i1,i2-1,i3-1)+8.*tz(i1,i2-1,i3+1)-tz(i1,i2-1,i3+2))/(12.*dr(2))+8.*(tz(i1,i2+1,i3-2)-8.*tz(i1,i2+1,i3-1)+8.*tz(i1,i2+1,i3+1)-tz(i1,i2+1,i3+2))/(12.*dr(2))-(tz(i1,i2+2,i3-2)-8.*tz(i1,i2+2,i3-1)+8.*tz(i1,i2+2,i3+1)-tz(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             tztt     = (-tz(i1,i2,i3-2)+16.*tz(i1,i2,i3-1)-30.*tz(i1,i2,i3)+16.*tz(i1,i2,i3+1)-tz(i1,i2,i3+2))/(12.*dr(2)**2)
                             ! ---------- Spatial derivatives of metrics rx, sx, ry, ... ---------
                             rxi = rx(i1,i2,i3)
                             ryi = ry(i1,i2,i3)
                             sxi = sx(i1,i2,i3)
                             syi = sy(i1,i2,i3)
                             rzi = rz(i1,i2,i3)
                             szi = sz(i1,i2,i3)
                             txi = tx(i1,i2,i3)
                             tyi = ty(i1,i2,i3)
                             tzi = tz(i1,i2,i3)
                             rxx      = rxi*rxr+sxi*rxs+txi*rxt
                             rxy      = ryi*rxr+syi*rxs+tyi*rxt
                             rxz      = rzi*rxr+szi*rxs+tzi*rxt
                             ryy      = ryi*ryr+syi*rys+tyi*ryt
                             ryz      = rzi*ryr+szi*rys+tzi*ryt
                             rzz      = rzi*rzr+szi*rzs+tzi*rzt
                             sxx      = rxi*sxr+sxi*sxs+txi*sxt
                             sxy      = ryi*sxr+syi*sxs+tyi*sxt
                             sxz      = rzi*sxr+szi*sxs+tzi*sxt
                             syy      = ryi*syr+syi*sys+tyi*syt
                             syz      = rzi*syr+szi*sys+tzi*syt
                             szz      = rzi*szr+szi*szs+tzi*szt
                             txx      = rxi*txr+sxi*txs+txi*txt
                             txy      = ryi*txr+syi*txs+tyi*txt
                             txz      = rzi*txr+szi*txs+tzi*txt
                             tyy      = ryi*tyr+syi*tys+tyi*tyt
                             tyz      = rzi*tyr+szi*tys+tzi*tyt
                             tzz      = rzi*tzr+szi*tzs+tzi*tzt
                             rxxx     = rxi**2*rxrr+2.*rxi*sxi*rxrs+2.*rxi*txi*rxrt+sxi**2*rxss+2.*sxi*txi*rxst+txi**2*rxtt+rxx*rxr+sxx*rxs+txx*rxt
                             rxxy     = rxi*ryi*rxrr+(rxi*syi+ryi*sxi)*rxrs+sxi*syi*rxss+(rxi*tyi+ryi*txi)*rxrt+(sxi*tyi+syi*txi)*rxst+txi*tyi*rxtt+rxy*rxr+sxy*rxs+txy*rxt
                             rxyy     = ryi**2*rxrr+2.*ryi*syi*rxrs+2.*ryi*tyi*rxrt+syi**2*rxss+2.*syi*tyi*rxst+tyi**2*rxtt+ryy*rxr+syy*rxs+tyy*rxt
                             rxxz     = rxi*rzi*rxrr+(rxi*szi+rzi*sxi)*rxrs+sxi*szi*rxss+(rxi*tzi+rzi*txi)*rxrt+(sxi*tzi+szi*txi)*rxst+txi*tzi*rxtt+rxz*rxr+sxz*rxs+txz*rxt
                             rxyz     = ryi*rzi*rxrr+(ryi*szi+rzi*syi)*rxrs+syi*szi*rxss+(ryi*tzi+rzi*tyi)*rxrt+(syi*tzi+szi*tyi)*rxst+tyi*tzi*rxtt+ryz*rxr+syz*rxs+tyz*rxt
                             rxzz     = rzi**2*rxrr+2.*rzi*szi*rxrs+2.*rzi*tzi*rxrt+szi**2*rxss+2.*szi*tzi*rxst+tzi**2*rxtt+rzz*rxr+szz*rxs+tzz*rxt
                             ryyy     = ryi**2*ryrr+2.*ryi*syi*ryrs+2.*ryi*tyi*ryrt+syi**2*ryss+2.*syi*tyi*ryst+tyi**2*rytt+ryy*ryr+syy*rys+tyy*ryt
                             ryyz     = ryi*rzi*ryrr+(ryi*szi+rzi*syi)*ryrs+syi*szi*ryss+(ryi*tzi+rzi*tyi)*ryrt+(syi*tzi+szi*tyi)*ryst+tyi*tzi*rytt+ryz*ryr+syz*rys+tyz*ryt
                             ryzz     = rzi**2*ryrr+2.*rzi*szi*ryrs+2.*rzi*tzi*ryrt+szi**2*ryss+2.*szi*tzi*ryst+tzi**2*rytt+rzz*ryr+szz*rys+tzz*ryt
                             rzzz     = rzi**2*rzrr+2.*rzi*szi*rzrs+2.*rzi*tzi*rzrt+szi**2*rzss+2.*szi*tzi*rzst+tzi**2*rztt+rzz*rzr+szz*rzs+tzz*rzt
                             sxxx     = rxi**2*sxrr+2.*rxi*sxi*sxrs+2.*rxi*txi*sxrt+sxi**2*sxss+2.*sxi*txi*sxst+txi**2*sxtt+rxx*sxr+sxx*sxs+txx*sxt
                             sxxy     = rxi*ryi*sxrr+(rxi*syi+ryi*sxi)*sxrs+sxi*syi*sxss+(rxi*tyi+ryi*txi)*sxrt+(sxi*tyi+syi*txi)*sxst+txi*tyi*sxtt+rxy*sxr+sxy*sxs+txy*sxt
                             sxyy     = ryi**2*sxrr+2.*ryi*syi*sxrs+2.*ryi*tyi*sxrt+syi**2*sxss+2.*syi*tyi*sxst+tyi**2*sxtt+ryy*sxr+syy*sxs+tyy*sxt
                             sxxz     = rxi*rzi*sxrr+(rxi*szi+rzi*sxi)*sxrs+sxi*szi*sxss+(rxi*tzi+rzi*txi)*sxrt+(sxi*tzi+szi*txi)*sxst+txi*tzi*sxtt+rxz*sxr+sxz*sxs+txz*sxt
                             sxyz     = ryi*rzi*sxrr+(ryi*szi+rzi*syi)*sxrs+syi*szi*sxss+(ryi*tzi+rzi*tyi)*sxrt+(syi*tzi+szi*tyi)*sxst+tyi*tzi*sxtt+ryz*sxr+syz*sxs+tyz*sxt
                             sxzz     = rzi**2*sxrr+2.*rzi*szi*sxrs+2.*rzi*tzi*sxrt+szi**2*sxss+2.*szi*tzi*sxst+tzi**2*sxtt+rzz*sxr+szz*sxs+tzz*sxt
                             syyy     = ryi**2*syrr+2.*ryi*syi*syrs+2.*ryi*tyi*syrt+syi**2*syss+2.*syi*tyi*syst+tyi**2*sytt+ryy*syr+syy*sys+tyy*syt
                             syyz     = ryi*rzi*syrr+(ryi*szi+rzi*syi)*syrs+syi*szi*syss+(ryi*tzi+rzi*tyi)*syrt+(syi*tzi+szi*tyi)*syst+tyi*tzi*sytt+ryz*syr+syz*sys+tyz*syt
                             syzz     = rzi**2*syrr+2.*rzi*szi*syrs+2.*rzi*tzi*syrt+szi**2*syss+2.*szi*tzi*syst+tzi**2*sytt+rzz*syr+szz*sys+tzz*syt
                             szzz     = rzi**2*szrr+2.*rzi*szi*szrs+2.*rzi*tzi*szrt+szi**2*szss+2.*szi*tzi*szst+tzi**2*sztt+rzz*szr+szz*szs+tzz*szt
                             txxx     = rxi**2*txrr+2.*rxi*sxi*txrs+2.*rxi*txi*txrt+sxi**2*txss+2.*sxi*txi*txst+txi**2*txtt+rxx*txr+sxx*txs+txx*txt
                             txxy     = rxi*ryi*txrr+(rxi*syi+ryi*sxi)*txrs+sxi*syi*txss+(rxi*tyi+ryi*txi)*txrt+(sxi*tyi+syi*txi)*txst+txi*tyi*txtt+rxy*txr+sxy*txs+txy*txt
                             txyy     = ryi**2*txrr+2.*ryi*syi*txrs+2.*ryi*tyi*txrt+syi**2*txss+2.*syi*tyi*txst+tyi**2*txtt+ryy*txr+syy*txs+tyy*txt
                             txxz     = rxi*rzi*txrr+(rxi*szi+rzi*sxi)*txrs+sxi*szi*txss+(rxi*tzi+rzi*txi)*txrt+(sxi*tzi+szi*txi)*txst+txi*tzi*txtt+rxz*txr+sxz*txs+txz*txt
                             txyz     = ryi*rzi*txrr+(ryi*szi+rzi*syi)*txrs+syi*szi*txss+(ryi*tzi+rzi*tyi)*txrt+(syi*tzi+szi*tyi)*txst+tyi*tzi*txtt+ryz*txr+syz*txs+tyz*txt
                             txzz     = rzi**2*txrr+2.*rzi*szi*txrs+2.*rzi*tzi*txrt+szi**2*txss+2.*szi*tzi*txst+tzi**2*txtt+rzz*txr+szz*txs+tzz*txt
                             tyyy     = ryi**2*tyrr+2.*ryi*syi*tyrs+2.*ryi*tyi*tyrt+syi**2*tyss+2.*syi*tyi*tyst+tyi**2*tytt+ryy*tyr+syy*tys+tyy*tyt
                             tyyz     = ryi*rzi*tyrr+(ryi*szi+rzi*syi)*tyrs+syi*szi*tyss+(ryi*tzi+rzi*tyi)*tyrt+(syi*tzi+szi*tyi)*tyst+tyi*tzi*tytt+ryz*tyr+syz*tys+tyz*tyt
                             tyzz     = rzi**2*tyrr+2.*rzi*szi*tyrs+2.*rzi*tzi*tyrt+szi**2*tyss+2.*szi*tzi*tyst+tzi**2*tytt+rzz*tyr+szz*tys+tzz*tyt
                             tzzz     = rzi**2*tzrr+2.*rzi*szi*tzrs+2.*rzi*tzi*tzrt+szi**2*tzss+2.*szi*tzi*tzst+tzi**2*tztt+rzz*tzr+szz*tzs+tzz*tzt
                             ! ---- end evalMetrics eq evalMetrics ---
                             ! ---------- Third spatial derivatives of u ---------
                             uxyz     = rxi*ryi*rzi*urrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*urrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*urss+sxi*syi*szi*usss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*urrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*urst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*usst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*urtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*ustt+txi*tyi*tzi*uttt+(rxi*ryz+rxy*rzi+rxz*ryi)*urr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*urs+(sxi*syz+sxy*szi+sxz*syi)*uss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*urt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*ust+(txi*tyz+txy*tzi+txz*tyi)*utt+rxyz*ur+sxyz*us+txyz*ut
                             ! ---------- END CURVILINEAR  ---------
                             ! ------ curvilinear grid: -------
                               ! 3D 
                               r1 =  an1*ux43(i1,i2,i3,0) + an2*uy43(i1,i2,i3,0) + an3*uz43(i1,i2,i3,0) - gg
                               crv(axis) = an1*rsxy(i1,i2,i3,axis,0) + an2*rsxy(i1,i2,i3,axis,1) + an3*rsxy(i1,i2,i3,axis,2)
                               ! **CHECK ME**
                               a11 = -is*( crv(axis)*8./(12.*dr(axis)) )  ! coeff of u(-1)
                               a12 =  is*( crv(axis)*1./(12.*dr(axis)) )  ! coeff of u(-2)
                               r2 = c2*( an1*( uxxx + uxyy + uxzz ) + an2*( uxxy + uyyy + uyzz ) + an3*( uxxz + uyyz + uzzz ) ) + nDotGradF - gtt
                               ! crv(axis) = an1*rsxy(i1,i2,i3,axis,0)**3 + an2*rsxy(i1,i2,i3,axis,1)**3 + an3*rsxy(i1,i2,i3,axis,2)**3
                               ! Coeff of "urrr" term *check me*
                               crv(axis) = (an1*rsxy(i1,i2,i3,axis,0) + an2*rsxy(i1,i2,i3,axis,1) + an3*rsxy(i1,i2,i3,axis,2) )* (    rsxy(i1,i2,i3,axis,0)**2  + rsxy(i1,i2,i3,axis,1)**2  + rsxy(i1,i2,i3,axis,2)**2 )
                               ! if( axis.eq.0 )then
                               !   crv(axis) = ( an1*rxi + an2*ryi + an3*rzi )*( rxi**2 + ryi**2 + rzi**2 ) 
                               !   ! crv(axis) = an1*( rxi*( rxi**2 + ryi**2 + rzi**2 ) ) + !   !             an2*( ryi*( rxi**2 + ryi**2 + rzi**2 ) ) + !   !             an3*( rzi*( rxi**2 + ryi**2 + rzi**2 ) ) 
                               ! else if( axis.eq.1 )then
                               !   crv(axis) = ( an1*sxi + an2*syi + an3*szi )*( sxi**2 + syi**2 + szi**2 ) 
                               ! else           
                               !   crv(axis) = ( an1*txi + an2*tyi + an3*tzi )*( txi**2 + tyi**2 + tzi**2 ) 
                               ! end if
                               ! **CHECK ME**
                               a21 =  is*c2*( 2.*crv(axis)/(2.*dr(axis)**3) )
                               a22 = -is*c2*( 1.*crv(axis)/(2.*dr(axis)**3) )            
                               ! define the residual functions for the discrete delta method
                                 if( checkCoeff.eq.1 )then
                                   ! --- Check coefficients a11,a12,... by discrete delta ---
                                   u1Save = u(j1,j2,j3,0)
                                   u2Save = u(k1,k2,k3,0)
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=0.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ut       = (u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))
                                   urt      = ((u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0))
                                   urrt     = (-(u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0)**2)
                                   ust      = ((u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))
                                   urst     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   usst     = (-(u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)**2)
                                   utt      = (-u(i1,i2,i3-2,0)+16.*u(i1,i2,i3-1,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2)**2)
                                   urtt     = ((-u(i1-2,i2,i3-2,0)+16.*u(i1-2,i2,i3-1,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1-1,i2,i3-2,0)+16.*u(i1-1,i2,i3-1,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1+1,i2,i3-2,0)+16.*u(i1+1,i2,i3-1,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2)**2)-(-u(i1+2,i2,i3-2,0)+16.*u(i1+2,i2,i3-1,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(0))
                                   ustt     = ((-u(i1,i2-2,i3-2,0)+16.*u(i1,i2-2,i3-1,0)-30.*u(i1,i2-2,i3,0)+16.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1,i2-1,i3-2,0)+16.*u(i1,i2-1,i3-1,0)-30.*u(i1,i2-1,i3,0)+16.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1,i2+1,i3-2,0)+16.*u(i1,i2+1,i3-1,0)-30.*u(i1,i2+1,i3,0)+16.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2)**2)-(-u(i1,i2+2,i3-2,0)+16.*u(i1,i2+2,i3-1,0)-30.*u(i1,i2+2,i3,0)+16.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(1))
                                   uttt     = (-u(i1,i2,i3-2,0)+2.*u(i1,i2,i3-1,0)-2.*u(i1,i2,i3+1,0)+u(i1,i2,i3+2,0))/(2.*dr(2)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   uxyz     = rxi*ryi*rzi*urrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*urrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*urss+sxi*syi*szi*usss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*urrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*urst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*usst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*urtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*ustt+txi*tyi*tzi*uttt+(rxi*ryz+rxy*rzi+rxz*ryi)*urr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*urs+(sxi*syz+sxy*szi+sxz*syi)*uss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*urt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*ust+(txi*tyz+txy*tzi+txz*tyi)*utt+rxyz*ur+sxyz*us+txyz*ut
                                   ! ---------- END CURVILINEAR  ---------
                                   r1a = (an1*ux43(i1,i2,i3,0)+an2*uy43(i1,i2,i3,0)+an3*uz43(i1,i2,i3,0))
                                   r2a = (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   u(j1,j2,j3,0)=1.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ut       = (u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))
                                   urt      = ((u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0))
                                   urrt     = (-(u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0)**2)
                                   ust      = ((u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))
                                   urst     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   usst     = (-(u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)**2)
                                   utt      = (-u(i1,i2,i3-2,0)+16.*u(i1,i2,i3-1,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2)**2)
                                   urtt     = ((-u(i1-2,i2,i3-2,0)+16.*u(i1-2,i2,i3-1,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1-1,i2,i3-2,0)+16.*u(i1-1,i2,i3-1,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1+1,i2,i3-2,0)+16.*u(i1+1,i2,i3-1,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2)**2)-(-u(i1+2,i2,i3-2,0)+16.*u(i1+2,i2,i3-1,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(0))
                                   ustt     = ((-u(i1,i2-2,i3-2,0)+16.*u(i1,i2-2,i3-1,0)-30.*u(i1,i2-2,i3,0)+16.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1,i2-1,i3-2,0)+16.*u(i1,i2-1,i3-1,0)-30.*u(i1,i2-1,i3,0)+16.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1,i2+1,i3-2,0)+16.*u(i1,i2+1,i3-1,0)-30.*u(i1,i2+1,i3,0)+16.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2)**2)-(-u(i1,i2+2,i3-2,0)+16.*u(i1,i2+2,i3-1,0)-30.*u(i1,i2+2,i3,0)+16.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(1))
                                   uttt     = (-u(i1,i2,i3-2,0)+2.*u(i1,i2,i3-1,0)-2.*u(i1,i2,i3+1,0)+u(i1,i2,i3+2,0))/(2.*dr(2)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   uxyz     = rxi*ryi*rzi*urrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*urrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*urss+sxi*syi*szi*usss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*urrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*urst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*usst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*urtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*ustt+txi*tyi*tzi*uttt+(rxi*ryz+rxy*rzi+rxz*ryi)*urr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*urs+(sxi*syz+sxy*szi+sxz*syi)*uss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*urt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*ust+(txi*tyz+txy*tzi+txz*tyi)*utt+rxyz*ur+sxyz*us+txyz*ut
                                   ! ---------- END CURVILINEAR  ---------
                                   r1b =  (an1*ux43(i1,i2,i3,0)+an2*uy43(i1,i2,i3,0)+an3*uz43(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   a11c = r1b - r1a
                                   a21c = r2b - r2a
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=1.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ut       = (u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))
                                   urt      = ((u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0))
                                   urrt     = (-(u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0)**2)
                                   ust      = ((u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))
                                   urst     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   usst     = (-(u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)**2)
                                   utt      = (-u(i1,i2,i3-2,0)+16.*u(i1,i2,i3-1,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2)**2)
                                   urtt     = ((-u(i1-2,i2,i3-2,0)+16.*u(i1-2,i2,i3-1,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1-1,i2,i3-2,0)+16.*u(i1-1,i2,i3-1,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1+1,i2,i3-2,0)+16.*u(i1+1,i2,i3-1,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2)**2)-(-u(i1+2,i2,i3-2,0)+16.*u(i1+2,i2,i3-1,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(0))
                                   ustt     = ((-u(i1,i2-2,i3-2,0)+16.*u(i1,i2-2,i3-1,0)-30.*u(i1,i2-2,i3,0)+16.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1,i2-1,i3-2,0)+16.*u(i1,i2-1,i3-1,0)-30.*u(i1,i2-1,i3,0)+16.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1,i2+1,i3-2,0)+16.*u(i1,i2+1,i3-1,0)-30.*u(i1,i2+1,i3,0)+16.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2)**2)-(-u(i1,i2+2,i3-2,0)+16.*u(i1,i2+2,i3-1,0)-30.*u(i1,i2+2,i3,0)+16.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(1))
                                   uttt     = (-u(i1,i2,i3-2,0)+2.*u(i1,i2,i3-1,0)-2.*u(i1,i2,i3+1,0)+u(i1,i2,i3+2,0))/(2.*dr(2)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   uxyz     = rxi*ryi*rzi*urrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*urrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*urss+sxi*syi*szi*usss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*urrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*urst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*usst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*urtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*ustt+txi*tyi*tzi*uttt+(rxi*ryz+rxy*rzi+rxz*ryi)*urr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*urs+(sxi*syz+sxy*szi+sxz*syi)*uss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*urt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*ust+(txi*tyz+txy*tzi+txz*tyi)*utt+rxyz*ur+sxyz*us+txyz*ut
                                   ! ---------- END CURVILINEAR  ---------
                                   r1b =  (an1*ux43(i1,i2,i3,0)+an2*uy43(i1,i2,i3,0)+an3*uz43(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   a12c = r1b - r1a
                                   a22c = r2b - r2a
                                   u(j1,j2,j3,0) = u1Save  
                                   u(k1,k2,k3,0) = u2Save      
                                   maxDiff=max(maxDiff,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c),abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c))
                                   ! #If "curvilinear" eq "curvilinear"
                                   ! write(*,'("++ i1,i2=",2i4," a11,a11c=",2(1pe12.4)," a12,a12c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a11,a11c,a12,a12c,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c)
                                   ! write(*,'("++ i1,i2=",2i4," a21,a21c=",2(1pe12.4)," a22,a22c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a21,a21c,a22,a22c,abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c)
                                   ! #End
                                 end if 
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2 )/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1 )/det
                           if( .false. )then
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,uex )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,uey )
                             write(*,'(" (i1,i2)=(",2i3,") u(-1),ue(-1)=",2(1pe12.4)," u(-2),ue(-2)=",2(1pe12.4))') i1,i2,uTemp(j1,j2,j3,0),uex,uTemp(k1,k2,k3,0),uey
                           end if
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! this case done below 
                         end if
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
                        end do
                        end do
                        end do
                       if( checkCoeff.eq.1 )then
                         write(*,'("bcOptWave: t=",e12.3," checkCoeff: (side,axis)=(",2i2,") rel-maxDiff=",e9.2," curvilinear 4")') t,side,axis,maxDiff
                       end if
                 end if
               end if 
             else
               if( gridType.eq.rectangular )then
                 if( nd.eq.2 )then
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
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       ! ------------------- 4 4 NEUMANN CBC --------------------
                       maxDiff=0. ! for checkCoeff
                       gg=0.
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         ghost =2
                         k1=i1-is1*ghost
                         k2=i2-is2*ghost
                         k3=i3-is3*ghost   
                         if( mask(i1,i2,i3).gt.0 )then
                               if( assignTwilightZone.eq.1 )then
                                 ! compute RHS from TZ
                                   call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uex )
                                   call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uey )
                                   gg = an1*uex + an2*uey
                                     call ogDeriv(ep,2,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uettx )
                                     call ogDeriv(ep,2,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uetty )
                                     gtt = an1*( uettx ) + an2*( uetty )
                                     call ogDeriv(ep,0,3,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxx )
                                     call ogDeriv(ep,0,1,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexyy )
                                     call ogDeriv(ep,0,2,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxy )
                                     call ogDeriv(ep,0,0,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyy )
                                     nDotGradF = an1*( uettx - c2*( uexxx + uexyy ) ) + an2*( uetty - c2*( uexxy + ueyyy ) )
                               else
                                 gg=0.;  gtt=0.; nDotGradF=0.; 
                               end if
                             ! ---------- RECTANGULAR  ---------
                             ! This assumes dr(0:2) = dx(0:2)
                             ! --- NEUMANN 4=4 rectangular ---
                             ! u_tt = c^2*Lap(u) + f 
                             !   u.n = g
                             ! g_tt = c^2 n.grad( Lap(u) ) + n.grad(f)
                             ! ux = [ u(-2) - 8*u(-1) + 8*u(1) - u(2) ]/(12*h) 
                             ! uxxx = (-u(-1) +2*u(-1) - 2*u(1) + u(2) ]/(2*h^3)
                               ! eval equation with wrong values at ghost: 
                               r1 =  an1*ux42r(i1,i2,i3,0) + an2*uy42r(i1,i2,i3,0) - gg
                               ! note: an1=-1 on left side ans an1=+1 on right side
                               ! **CHECK ME**
                               a11 =  8./(12.*dx(axis))  ! coeff of u(-1)
                               a12 = -1./(12.*dx(axis))  ! coeff of u(-2)
                               r2 = c2*( an1*( uxxx + uxyy ) + an2*( uxxy + uyyy ) ) + nDotGradF - gtt
                               ! **CHECK ME**
                               a21 = -c2*( 2./(2.*dx(axis)**3) )
                               a22 =  c2*( 1./(2.*dx(axis)**3) )
                               ! define the residual functions for the discrete delta method
                                 if( checkCoeff.eq.1 )then
                                   ! --- Check coefficients a11,a12,... by discrete delta ---
                                   u1Save = u(j1,j2,j3,0)
                                   u2Save = u(k1,k2,k3,0)
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=0.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   r1a = (an1*ux42r(i1,i2,i3,0)+an2*uy42r(i1,i2,i3,0))
                                   r2a = (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   u(j1,j2,j3,0)=1.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   r1b =  (an1*ux42r(i1,i2,i3,0)+an2*uy42r(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   a11c = r1b - r1a
                                   a21c = r2b - r2a
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=1.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   r1b =  (an1*ux42r(i1,i2,i3,0)+an2*uy42r(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   a12c = r1b - r1a
                                   a22c = r2b - r2a
                                   u(j1,j2,j3,0) = u1Save  
                                   u(k1,k2,k3,0) = u2Save      
                                   maxDiff=max(maxDiff,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c),abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c))
                                   ! #If "rectangular" eq "curvilinear"
                                   ! write(*,'("++ i1,i2=",2i4," a11,a11c=",2(1pe12.4)," a12,a12c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a11,a11c,a12,a12c,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c)
                                   ! write(*,'("++ i1,i2=",2i4," a21,a21c=",2(1pe12.4)," a22,a22c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a21,a21c,a22,a22c,abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c)
                                   ! #End
                                 end if 
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2 )/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1 )/det
                           if( .false. )then
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,uex )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,uey )
                             write(*,'(" (i1,i2)=(",2i3,") u(-1),ue(-1)=",2(1pe12.4)," u(-2),ue(-2)=",2(1pe12.4))') i1,i2,uTemp(j1,j2,j3,0),uex,uTemp(k1,k2,k3,0),uey
                           end if
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! this case done below 
                         end if
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
                        end do
                        end do
                        end do
                       if( checkCoeff.eq.1 )then
                         write(*,'("bcOptWave: t=",e12.3," checkCoeff: (side,axis)=(",2i2,") rel-maxDiff=",e9.2," rectangular 4")') t,side,axis,maxDiff
                       end if
                 else
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
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       ! ------------------- 4 4 NEUMANN CBC --------------------
                       maxDiff=0. ! for checkCoeff
                       gg=0.
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         ghost =2
                         k1=i1-is1*ghost
                         k2=i2-is2*ghost
                         k3=i3-is3*ghost   
                         if( mask(i1,i2,i3).gt.0 )then
                               if( assignTwilightZone.eq.1 )then
                                 ! compute RHS from TZ
                                   ! ----- 3D  -----
                                   call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uex )
                                   call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uey )
                                   call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uez )
                                   gg = an1*uex + an2*uey + an3*uez
                                     call ogDeriv(ep,2,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettx )
                                     call ogDeriv(ep,2,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uetty )
                                     call ogDeriv(ep,2,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettz )
                                     gtt = an1*( uettx ) + an2*( uetty ) + an3*( uettz )
                                     call ogDeriv(ep,0,3,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxx )
                                     call ogDeriv(ep,0,2,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxy )
                                     call ogDeriv(ep,0,1,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexyy )
                                     call ogDeriv(ep,0,0,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyy )
                                     call ogDeriv(ep,0,2,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxz )
                                     call ogDeriv(ep,0,1,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexzz )
                                     call ogDeriv(ep,0,0,2,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyz )
                                     call ogDeriv(ep,0,0,1,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyzz )
                                     call ogDeriv(ep,0,0,0,3,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezzz )
                                     nDotGradF = an1*( uettx - c2*( uexxx + uexyy + uexzz ) ) + an2*( uetty - c2*( uexxy + ueyyy + ueyzz ) ) + an3*( uettz - c2*( uexxz + ueyyz + uezzz ) )
                               else
                                 gg=0.;  gtt=0.; nDotGradF=0.; 
                               end if
                             ! evaluate 3rd derivatives : uxxx,uxxy,uxxz,uxyy,uxzz, uyyy,uyyz,uyzz, uzzz 
                             ! ---------- RECTANGULAR  ---------
                             ! This assumes dr(0:2) = dx(0:2)
                             uxyz     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                             ! --- NEUMANN 4=4 rectangular ---
                             ! u_tt = c^2*Lap(u) + f 
                             !   u.n = g
                             ! g_tt = c^2 n.grad( Lap(u) ) + n.grad(f)
                             ! ux = [ u(-2) - 8*u(-1) + 8*u(1) - u(2) ]/(12*h) 
                             ! uxxx = (-u(-1) +2*u(-1) - 2*u(1) + u(2) ]/(2*h^3)
                               ! 3D 
                               r1 =  an1*ux43r(i1,i2,i3,0) + an2*uy43r(i1,i2,i3,0) + an3*uz43r(i1,i2,i3,0) - gg
                               ! **CHECK ME**
                               a11 =  8./(12.*dx(axis))  ! coeff of u(-1)
                               a12 = -1./(12.*dx(axis))  ! coeff of u(-2)  
                               r2 = c2*( an1*( uxxx + uxyy + uxzz ) + an2*( uxxy + uyyy + uyzz ) + an3*( uxxz + uyyz + uzzz ) ) + nDotGradF - gtt
                               ! **CHECK ME**
                               a21 = -c2*( 2./(2.*dx(axis)**3) )
                               a22 =  c2*( 1./(2.*dx(axis)**3) )
                               ! define the residual functions for the discrete delta method
                                 if( checkCoeff.eq.1 )then
                                   ! --- Check coefficients a11,a12,... by discrete delta ---
                                   u1Save = u(j1,j2,j3,0)
                                   u2Save = u(k1,k2,k3,0)
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=0.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   uxyz     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   r1a = (an1*ux43r(i1,i2,i3,0)+an2*uy43r(i1,i2,i3,0)+an3*uz43r(i1,i2,i3,0))
                                   r2a = (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   u(j1,j2,j3,0)=1.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   uxyz     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   r1b =  (an1*ux43r(i1,i2,i3,0)+an2*uy43r(i1,i2,i3,0)+an3*uz43r(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   a11c = r1b - r1a
                                   a21c = r2b - r2a
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=1.
                                   ! ---------- RECTANGULAR  ---------
                                   ! This assumes dr(0:2) = dx(0:2)
                                   uxyz     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   r1b =  (an1*ux43r(i1,i2,i3,0)+an2*uy43r(i1,i2,i3,0)+an3*uz43r(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   a12c = r1b - r1a
                                   a22c = r2b - r2a
                                   u(j1,j2,j3,0) = u1Save  
                                   u(k1,k2,k3,0) = u2Save      
                                   maxDiff=max(maxDiff,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c),abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c))
                                   ! #If "rectangular" eq "curvilinear"
                                   ! write(*,'("++ i1,i2=",2i4," a11,a11c=",2(1pe12.4)," a12,a12c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a11,a11c,a12,a12c,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c)
                                   ! write(*,'("++ i1,i2=",2i4," a21,a21c=",2(1pe12.4)," a22,a22c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a21,a21c,a22,a22c,abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c)
                                   ! #End
                                 end if 
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2 )/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1 )/det
                           if( .false. )then
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,uex )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,uey )
                             write(*,'(" (i1,i2)=(",2i3,") u(-1),ue(-1)=",2(1pe12.4)," u(-2),ue(-2)=",2(1pe12.4))') i1,i2,uTemp(j1,j2,j3,0),uex,uTemp(k1,k2,k3,0),uey
                           end if
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! this case done below 
                         end if
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
                        end do
                        end do
                        end do
                       if( checkCoeff.eq.1 )then
                         write(*,'("bcOptWave: t=",e12.3," checkCoeff: (side,axis)=(",2i2,") rel-maxDiff=",e9.2," rectangular 4")') t,side,axis,maxDiff
                       end if
                 end if
               else
                 if( nd.eq.2 )then
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       ! ------------------- 4 4 NEUMANN CBC --------------------
                       maxDiff=0. ! for checkCoeff
                       gg=0.
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         ghost =2
                         k1=i1-is1*ghost
                         k2=i2-is2*ghost
                         k3=i3-is3*ghost   
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
                               if( assignTwilightZone.eq.1 )then
                                 ! compute RHS from TZ
                                   call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uex )
                                   call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uey )
                                   gg = an1*uex + an2*uey
                                     call ogDeriv(ep,2,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uettx )
                                     call ogDeriv(ep,2,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uetty )
                                     gtt = an1*( uettx ) + an2*( uetty )
                                     call ogDeriv(ep,0,3,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxx )
                                     call ogDeriv(ep,0,1,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexyy )
                                     call ogDeriv(ep,0,2,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxy )
                                     call ogDeriv(ep,0,0,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyy )
                                     nDotGradF = an1*( uettx - c2*( uexxx + uexyy ) ) + an2*( uetty - c2*( uexxy + ueyyy ) )
                               else
                                 gg=0.;  gtt=0.; nDotGradF=0.; 
                               end if
                             ! ---------- START CURVILINEAR  ---------
                             ! ---------- Parametric derivatives ---------
                             ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                             urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                             urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                             us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                             urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                             urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                             uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                             urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                             usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                             rxr      = (rx(i1-2,i2,i3)-8.*rx(i1-1,i2,i3)+8.*rx(i1+1,i2,i3)-rx(i1+2,i2,i3))/(12.*dr(0))
                             rxrr     = (-rx(i1-2,i2,i3)+16.*rx(i1-1,i2,i3)-30.*rx(i1,i2,i3)+16.*rx(i1+1,i2,i3)-rx(i1+2,i2,i3))/(12.*dr(0)**2)
                             rxs      = (rx(i1,i2-2,i3)-8.*rx(i1,i2-1,i3)+8.*rx(i1,i2+1,i3)-rx(i1,i2+2,i3))/(12.*dr(1))
                             rxrs     = ((rx(i1-2,i2-2,i3)-8.*rx(i1-2,i2-1,i3)+8.*rx(i1-2,i2+1,i3)-rx(i1-2,i2+2,i3))/(12.*dr(1))-8.*(rx(i1-1,i2-2,i3)-8.*rx(i1-1,i2-1,i3)+8.*rx(i1-1,i2+1,i3)-rx(i1-1,i2+2,i3))/(12.*dr(1))+8.*(rx(i1+1,i2-2,i3)-8.*rx(i1+1,i2-1,i3)+8.*rx(i1+1,i2+1,i3)-rx(i1+1,i2+2,i3))/(12.*dr(1))-(rx(i1+2,i2-2,i3)-8.*rx(i1+2,i2-1,i3)+8.*rx(i1+2,i2+1,i3)-rx(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             rxss     = (-rx(i1,i2-2,i3)+16.*rx(i1,i2-1,i3)-30.*rx(i1,i2,i3)+16.*rx(i1,i2+1,i3)-rx(i1,i2+2,i3))/(12.*dr(1)**2)
                             ryr      = (ry(i1-2,i2,i3)-8.*ry(i1-1,i2,i3)+8.*ry(i1+1,i2,i3)-ry(i1+2,i2,i3))/(12.*dr(0))
                             ryrr     = (-ry(i1-2,i2,i3)+16.*ry(i1-1,i2,i3)-30.*ry(i1,i2,i3)+16.*ry(i1+1,i2,i3)-ry(i1+2,i2,i3))/(12.*dr(0)**2)
                             rys      = (ry(i1,i2-2,i3)-8.*ry(i1,i2-1,i3)+8.*ry(i1,i2+1,i3)-ry(i1,i2+2,i3))/(12.*dr(1))
                             ryrs     = ((ry(i1-2,i2-2,i3)-8.*ry(i1-2,i2-1,i3)+8.*ry(i1-2,i2+1,i3)-ry(i1-2,i2+2,i3))/(12.*dr(1))-8.*(ry(i1-1,i2-2,i3)-8.*ry(i1-1,i2-1,i3)+8.*ry(i1-1,i2+1,i3)-ry(i1-1,i2+2,i3))/(12.*dr(1))+8.*(ry(i1+1,i2-2,i3)-8.*ry(i1+1,i2-1,i3)+8.*ry(i1+1,i2+1,i3)-ry(i1+1,i2+2,i3))/(12.*dr(1))-(ry(i1+2,i2-2,i3)-8.*ry(i1+2,i2-1,i3)+8.*ry(i1+2,i2+1,i3)-ry(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             ryss     = (-ry(i1,i2-2,i3)+16.*ry(i1,i2-1,i3)-30.*ry(i1,i2,i3)+16.*ry(i1,i2+1,i3)-ry(i1,i2+2,i3))/(12.*dr(1)**2)
                             sxr      = (sx(i1-2,i2,i3)-8.*sx(i1-1,i2,i3)+8.*sx(i1+1,i2,i3)-sx(i1+2,i2,i3))/(12.*dr(0))
                             sxrr     = (-sx(i1-2,i2,i3)+16.*sx(i1-1,i2,i3)-30.*sx(i1,i2,i3)+16.*sx(i1+1,i2,i3)-sx(i1+2,i2,i3))/(12.*dr(0)**2)
                             sxs      = (sx(i1,i2-2,i3)-8.*sx(i1,i2-1,i3)+8.*sx(i1,i2+1,i3)-sx(i1,i2+2,i3))/(12.*dr(1))
                             sxrs     = ((sx(i1-2,i2-2,i3)-8.*sx(i1-2,i2-1,i3)+8.*sx(i1-2,i2+1,i3)-sx(i1-2,i2+2,i3))/(12.*dr(1))-8.*(sx(i1-1,i2-2,i3)-8.*sx(i1-1,i2-1,i3)+8.*sx(i1-1,i2+1,i3)-sx(i1-1,i2+2,i3))/(12.*dr(1))+8.*(sx(i1+1,i2-2,i3)-8.*sx(i1+1,i2-1,i3)+8.*sx(i1+1,i2+1,i3)-sx(i1+1,i2+2,i3))/(12.*dr(1))-(sx(i1+2,i2-2,i3)-8.*sx(i1+2,i2-1,i3)+8.*sx(i1+2,i2+1,i3)-sx(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             sxss     = (-sx(i1,i2-2,i3)+16.*sx(i1,i2-1,i3)-30.*sx(i1,i2,i3)+16.*sx(i1,i2+1,i3)-sx(i1,i2+2,i3))/(12.*dr(1)**2)
                             syr      = (sy(i1-2,i2,i3)-8.*sy(i1-1,i2,i3)+8.*sy(i1+1,i2,i3)-sy(i1+2,i2,i3))/(12.*dr(0))
                             syrr     = (-sy(i1-2,i2,i3)+16.*sy(i1-1,i2,i3)-30.*sy(i1,i2,i3)+16.*sy(i1+1,i2,i3)-sy(i1+2,i2,i3))/(12.*dr(0)**2)
                             sys      = (sy(i1,i2-2,i3)-8.*sy(i1,i2-1,i3)+8.*sy(i1,i2+1,i3)-sy(i1,i2+2,i3))/(12.*dr(1))
                             syrs     = ((sy(i1-2,i2-2,i3)-8.*sy(i1-2,i2-1,i3)+8.*sy(i1-2,i2+1,i3)-sy(i1-2,i2+2,i3))/(12.*dr(1))-8.*(sy(i1-1,i2-2,i3)-8.*sy(i1-1,i2-1,i3)+8.*sy(i1-1,i2+1,i3)-sy(i1-1,i2+2,i3))/(12.*dr(1))+8.*(sy(i1+1,i2-2,i3)-8.*sy(i1+1,i2-1,i3)+8.*sy(i1+1,i2+1,i3)-sy(i1+1,i2+2,i3))/(12.*dr(1))-(sy(i1+2,i2-2,i3)-8.*sy(i1+2,i2-1,i3)+8.*sy(i1+2,i2+1,i3)-sy(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             syss     = (-sy(i1,i2-2,i3)+16.*sy(i1,i2-1,i3)-30.*sy(i1,i2,i3)+16.*sy(i1,i2+1,i3)-sy(i1,i2+2,i3))/(12.*dr(1)**2)
                             ! ---------- Spatial derivatives of metrics rx, sx, ry, ... ---------
                             rxi = rx(i1,i2,i3)
                             ryi = ry(i1,i2,i3)
                             sxi = sx(i1,i2,i3)
                             syi = sy(i1,i2,i3)
                             rxx      = rxi*rxr+sxi*rxs
                             rxy      = ryi*rxr+syi*rxs
                             ryy      = ryi*ryr+syi*rys
                             sxx      = rxi*sxr+sxi*sxs
                             sxy      = ryi*sxr+syi*sxs
                             syy      = ryi*syr+syi*sys
                             rxxx     = rxi**2*rxrr+2.*rxi*sxi*rxrs+sxi**2*rxss+rxx*rxr+sxx*rxs
                             rxxy     = rxi*ryi*rxrr+(rxi*syi+ryi*sxi)*rxrs+sxi*syi*rxss+rxy*rxr+sxy*rxs
                             rxyy     = ryi**2*rxrr+2.*ryi*syi*rxrs+syi**2*rxss+ryy*rxr+syy*rxs
                             ryyy     = ryi**2*ryrr+2.*ryi*syi*ryrs+syi**2*ryss+ryy*ryr+syy*rys
                             sxxx     = rxi**2*sxrr+2.*rxi*sxi*sxrs+sxi**2*sxss+rxx*sxr+sxx*sxs
                             sxxy     = rxi*ryi*sxrr+(rxi*syi+ryi*sxi)*sxrs+sxi*syi*sxss+rxy*sxr+sxy*sxs
                             sxyy     = ryi**2*sxrr+2.*ryi*syi*sxrs+syi**2*sxss+ryy*sxr+syy*sxs
                             syyy     = ryi**2*syrr+2.*ryi*syi*syrs+syi**2*syss+ryy*syr+syy*sys
                             ! ---- end evalMetrics eq evalMetrics ---
                             ! ---------- Third spatial derivatives of u ---------
                             ! ---------- END CURVILINEAR  ---------
                             ! ------ curvilinear grid: -------
                               ! eval equation with wrong values at ghost:
                               r1 =  an1*ux42(i1,i2,i3,0) + an2*uy42(i1,i2,i3,0) - gg
                               crv(axis) = an1*rsxy(i1,i2,i3,axis,0) + an2*rsxy(i1,i2,i3,axis,1)
                               ! crv(1) = an1*rsxy(i1,i2,i3,1,0) + an2*rsxy(i1,i2,i3,1,1)
                               ! **CHECK ME**
                               a11 = -is*( crv(axis)*8./(12.*dr(axis)) )  ! coeff of u(-1)
                               a12 =  is*( crv(axis)*1./(12.*dr(axis)) )  ! coeff of u(-2)
                               r2 = c2*( an1*( uxxx + uxyy ) + an2*( uxxy + uyyy ) ) + nDotGradF - gtt
                   ! uxxx = rxi**3*urrr+3.*rxi**2*sxi*urrs+3.*rxi*sxi**2*urss+sxi**3*usss+3.*rxi*rxx*urr+(3.*rxi*sxx+3.*rxx*sxi)*urs+3.*sxx*sxi*uss+rxxx*ur+sxxx*us
                   ! uxyy = ryi**2*rxi*urrr+(syi*ryi*rxi+ryi*(rxi*syi+ryi*sxi))*urrs+(ryi*syi*sxi+syi*(rxi*syi+ryi*sxi))*urss+syi**2*sxi*usss+(rxi*ryy+2.*rxy*ryi)*urr+(rxi*syy+2.*rxy*syi+2.*ryi*sxy+ryy*sxi)*urs+(sxi*syy+2.*sxy*syi)*uss+rxyy*ur+sxyy*us
                   ! uxxy = ryi*rxi**2*urrr+(rxi**2*syi+2.*rxi*ryi*sxi)*urrs+(2.*rxi*sxi*syi+ryi*sxi**2)*urss+syi*sxi**2*usss+(2.*rxi*rxy+rxx*ryi)*urr+(2.*rxi*sxy+rxx*syi+2.*rxy*sxi+ryi*sxx)*urs+(2.*sxi*sxy+sxx*syi)*uss+rxxy*ur+sxxy*us
                   ! uyyy = ryi**3*urrr+3.*ryi**2*syi*urrs+3.*ryi*syi**2*urss+syi**3*usss+3.*ryi*ryy*urr+(3.*ryi*syy+3.*ryy*syi)*urs+3.*syy*syi*uss+ryyy*ur+syyy*us
                               ! uxxx = rxi**3*urrr +3.*rxi*rxx*urr + rxxx*ur
                               ! uxyy = ryi**2*rxi*urrr +(rxi*ryy+2.*rxy*ryi)*urr+rxyy*ur
                               ! uxxy = ryi*rxi**2*urrr +(2.*rxi*rxy+rxx*ryi)*urr+rxxy*ur
                               ! Coeff of terms involving urrr (do not include terms involving urr and ur as these used 2nd-order values)
                               crv(axis) = (an1*rsxy(i1,i2,i3,axis,0) + an2*rsxy(i1,i2,i3,axis,1) )*( rsxy(i1,i2,i3,axis,0)**2 + rsxy(i1,i2,i3,axis,1)**2 )
                               ! if( axis.eq.0 )then
                               !   crv(axis) = (an1*rxi + an2*ryi )*( rxi**2 + ryi**2 )
                               !   ! crv(axis) = an1*( rxi**3 + ryi**2*rxi ) + !   !             an2*( ryi**3 + ryi*rxi**2 )
                               ! else
                               !   crv(axis) = (an1*sxi + an2*syi )*( sxi**2 + syi**2 )
                               !   ! crv(axis) = an1*( sxi**3 + syi**2*sxi ) + !   !             an2*( syi**3 + syi*sxi**2 )              
                               ! end if
                               ! crv(1) = an1*rsxy(i1,i2,i3,1,0)**3 + an2*rsxy(i1,i2,i3,1,1)**3
                               ! **CHECK ME**
                               a21 =  is*c2*( 2.*crv(axis)/(2.*dr(axis)**3) )
                               a22 = -is*c2*( 1.*crv(axis)/(2.*dr(axis)**3) )
                                 if( checkCoeff.eq.1 )then
                                   ! --- Check coefficients a11,a12,... by discrete delta ---
                                   u1Save = u(j1,j2,j3,0)
                                   u2Save = u(k1,k2,k3,0)
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=0.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   ! ---------- END CURVILINEAR  ---------
                                   r1a = (an1*ux42(i1,i2,i3,0)+an2*uy42(i1,i2,i3,0))
                                   r2a = (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   u(j1,j2,j3,0)=1.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   ! ---------- END CURVILINEAR  ---------
                                   r1b =  (an1*ux42(i1,i2,i3,0)+an2*uy42(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   a11c = r1b - r1a
                                   a21c = r2b - r2a
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=1.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   ! ---------- END CURVILINEAR  ---------
                                   r1b =  (an1*ux42(i1,i2,i3,0)+an2*uy42(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy)+an2*(uxxy+uyyy)))
                                   a12c = r1b - r1a
                                   a22c = r2b - r2a
                                   u(j1,j2,j3,0) = u1Save  
                                   u(k1,k2,k3,0) = u2Save      
                                   maxDiff=max(maxDiff,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c),abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c))
                                   ! #If "curvilinear" eq "curvilinear"
                                   ! write(*,'("++ i1,i2=",2i4," a11,a11c=",2(1pe12.4)," a12,a12c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a11,a11c,a12,a12c,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c)
                                   ! write(*,'("++ i1,i2=",2i4," a21,a21c=",2(1pe12.4)," a22,a22c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a21,a21c,a22,a22c,abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c)
                                   ! #End
                                 end if 
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2 )/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1 )/det
                           if( .false. )then
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,uex )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,uey )
                             write(*,'(" (i1,i2)=(",2i3,") u(-1),ue(-1)=",2(1pe12.4)," u(-2),ue(-2)=",2(1pe12.4))') i1,i2,uTemp(j1,j2,j3,0),uex,uTemp(k1,k2,k3,0),uey
                           end if
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! this case done below 
                         end if
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
                        end do
                        end do
                        end do
                       if( checkCoeff.eq.1 )then
                         write(*,'("bcOptWave: t=",e12.3," checkCoeff: (side,axis)=(",2i2,") rel-maxDiff=",e9.2," curvilinear 4")') t,side,axis,maxDiff
                       end if
                 else
                     ! Mixed BC: a0*u + a1*u.n = 
                     a0=0.
                     a1=1.
                     ! STAGE I always fill in first ghost from 2nd-order scheme 
                     ! rectangular case:
                     !---------------------------------------------------------------
                     ! --- STAGE I fill in first ghost by 2nd-order compatibility ---
                     !---------------------------------------------------------------
                     ! assign extram points in the tangential directions
                     extram = numGhost-1 
                       ! ------------------- 4 4 NEUMANN CBC --------------------
                       maxDiff=0. ! for checkCoeff
                       gg=0.
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         ghost =1 
                         j1=i1-is1*ghost
                         j2=i2-is2*ghost
                         j3=i3-is3*ghost
                         ghost =2
                         k1=i1-is1*ghost
                         k2=i2-is2*ghost
                         k3=i3-is3*ghost   
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
                               if( assignTwilightZone.eq.1 )then
                                 ! compute RHS from TZ
                                   ! ----- 3D  -----
                                   call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uex )
                                   call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uey )
                                   call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uez )
                                   gg = an1*uex + an2*uey + an3*uez
                                     call ogDeriv(ep,2,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettx )
                                     call ogDeriv(ep,2,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uetty )
                                     call ogDeriv(ep,2,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettz )
                                     gtt = an1*( uettx ) + an2*( uetty ) + an3*( uettz )
                                     call ogDeriv(ep,0,3,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxx )
                                     call ogDeriv(ep,0,2,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxy )
                                     call ogDeriv(ep,0,1,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexyy )
                                     call ogDeriv(ep,0,0,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyy )
                                     call ogDeriv(ep,0,2,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxz )
                                     call ogDeriv(ep,0,1,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexzz )
                                     call ogDeriv(ep,0,0,2,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyz )
                                     call ogDeriv(ep,0,0,1,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyzz )
                                     call ogDeriv(ep,0,0,0,3,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezzz )
                                     nDotGradF = an1*( uettx - c2*( uexxx + uexyy + uexzz ) ) + an2*( uetty - c2*( uexxy + ueyyy + ueyzz ) ) + an3*( uettz - c2*( uexxz + ueyyz + uezzz ) )
                               else
                                 gg=0.;  gtt=0.; nDotGradF=0.; 
                               end if
                             ! evaluate 3rd derivatives : uxxx,uxxy,uxxz,uxyy,uxzz, uyyy,uyyz,uyzz, uzzz 
                             ! ---------- START CURVILINEAR  ---------
                             ! ---------- Parametric derivatives ---------
                             ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                             urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                             urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                             us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                             urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                             urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                             uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                             urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                             usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                             ut       = (u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))
                             urt      = ((u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0))
                             urrt     = (-(u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0)**2)
                             ust      = ((u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))
                             urst     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                             usst     = (-(u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)**2)
                             utt      = (-u(i1,i2,i3-2,0)+16.*u(i1,i2,i3-1,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2)**2)
                             urtt     = ((-u(i1-2,i2,i3-2,0)+16.*u(i1-2,i2,i3-1,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1-1,i2,i3-2,0)+16.*u(i1-1,i2,i3-1,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1+1,i2,i3-2,0)+16.*u(i1+1,i2,i3-1,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2)**2)-(-u(i1+2,i2,i3-2,0)+16.*u(i1+2,i2,i3-1,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(0))
                             ustt     = ((-u(i1,i2-2,i3-2,0)+16.*u(i1,i2-2,i3-1,0)-30.*u(i1,i2-2,i3,0)+16.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1,i2-1,i3-2,0)+16.*u(i1,i2-1,i3-1,0)-30.*u(i1,i2-1,i3,0)+16.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1,i2+1,i3-2,0)+16.*u(i1,i2+1,i3-1,0)-30.*u(i1,i2+1,i3,0)+16.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2)**2)-(-u(i1,i2+2,i3-2,0)+16.*u(i1,i2+2,i3-1,0)-30.*u(i1,i2+2,i3,0)+16.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(1))
                             uttt     = (-u(i1,i2,i3-2,0)+2.*u(i1,i2,i3-1,0)-2.*u(i1,i2,i3+1,0)+u(i1,i2,i3+2,0))/(2.*dr(2)**3)
                             rxr      = (rx(i1-2,i2,i3)-8.*rx(i1-1,i2,i3)+8.*rx(i1+1,i2,i3)-rx(i1+2,i2,i3))/(12.*dr(0))
                             rxrr     = (-rx(i1-2,i2,i3)+16.*rx(i1-1,i2,i3)-30.*rx(i1,i2,i3)+16.*rx(i1+1,i2,i3)-rx(i1+2,i2,i3))/(12.*dr(0)**2)
                             rxs      = (rx(i1,i2-2,i3)-8.*rx(i1,i2-1,i3)+8.*rx(i1,i2+1,i3)-rx(i1,i2+2,i3))/(12.*dr(1))
                             rxrs     = ((rx(i1-2,i2-2,i3)-8.*rx(i1-2,i2-1,i3)+8.*rx(i1-2,i2+1,i3)-rx(i1-2,i2+2,i3))/(12.*dr(1))-8.*(rx(i1-1,i2-2,i3)-8.*rx(i1-1,i2-1,i3)+8.*rx(i1-1,i2+1,i3)-rx(i1-1,i2+2,i3))/(12.*dr(1))+8.*(rx(i1+1,i2-2,i3)-8.*rx(i1+1,i2-1,i3)+8.*rx(i1+1,i2+1,i3)-rx(i1+1,i2+2,i3))/(12.*dr(1))-(rx(i1+2,i2-2,i3)-8.*rx(i1+2,i2-1,i3)+8.*rx(i1+2,i2+1,i3)-rx(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             rxss     = (-rx(i1,i2-2,i3)+16.*rx(i1,i2-1,i3)-30.*rx(i1,i2,i3)+16.*rx(i1,i2+1,i3)-rx(i1,i2+2,i3))/(12.*dr(1)**2)
                             rxt      = (rx(i1,i2,i3-2)-8.*rx(i1,i2,i3-1)+8.*rx(i1,i2,i3+1)-rx(i1,i2,i3+2))/(12.*dr(2))
                             rxrt     = ((rx(i1-2,i2,i3-2)-8.*rx(i1-2,i2,i3-1)+8.*rx(i1-2,i2,i3+1)-rx(i1-2,i2,i3+2))/(12.*dr(2))-8.*(rx(i1-1,i2,i3-2)-8.*rx(i1-1,i2,i3-1)+8.*rx(i1-1,i2,i3+1)-rx(i1-1,i2,i3+2))/(12.*dr(2))+8.*(rx(i1+1,i2,i3-2)-8.*rx(i1+1,i2,i3-1)+8.*rx(i1+1,i2,i3+1)-rx(i1+1,i2,i3+2))/(12.*dr(2))-(rx(i1+2,i2,i3-2)-8.*rx(i1+2,i2,i3-1)+8.*rx(i1+2,i2,i3+1)-rx(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             rxst     = ((rx(i1,i2-2,i3-2)-8.*rx(i1,i2-2,i3-1)+8.*rx(i1,i2-2,i3+1)-rx(i1,i2-2,i3+2))/(12.*dr(2))-8.*(rx(i1,i2-1,i3-2)-8.*rx(i1,i2-1,i3-1)+8.*rx(i1,i2-1,i3+1)-rx(i1,i2-1,i3+2))/(12.*dr(2))+8.*(rx(i1,i2+1,i3-2)-8.*rx(i1,i2+1,i3-1)+8.*rx(i1,i2+1,i3+1)-rx(i1,i2+1,i3+2))/(12.*dr(2))-(rx(i1,i2+2,i3-2)-8.*rx(i1,i2+2,i3-1)+8.*rx(i1,i2+2,i3+1)-rx(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             rxtt     = (-rx(i1,i2,i3-2)+16.*rx(i1,i2,i3-1)-30.*rx(i1,i2,i3)+16.*rx(i1,i2,i3+1)-rx(i1,i2,i3+2))/(12.*dr(2)**2)
                             ryr      = (ry(i1-2,i2,i3)-8.*ry(i1-1,i2,i3)+8.*ry(i1+1,i2,i3)-ry(i1+2,i2,i3))/(12.*dr(0))
                             ryrr     = (-ry(i1-2,i2,i3)+16.*ry(i1-1,i2,i3)-30.*ry(i1,i2,i3)+16.*ry(i1+1,i2,i3)-ry(i1+2,i2,i3))/(12.*dr(0)**2)
                             rys      = (ry(i1,i2-2,i3)-8.*ry(i1,i2-1,i3)+8.*ry(i1,i2+1,i3)-ry(i1,i2+2,i3))/(12.*dr(1))
                             ryrs     = ((ry(i1-2,i2-2,i3)-8.*ry(i1-2,i2-1,i3)+8.*ry(i1-2,i2+1,i3)-ry(i1-2,i2+2,i3))/(12.*dr(1))-8.*(ry(i1-1,i2-2,i3)-8.*ry(i1-1,i2-1,i3)+8.*ry(i1-1,i2+1,i3)-ry(i1-1,i2+2,i3))/(12.*dr(1))+8.*(ry(i1+1,i2-2,i3)-8.*ry(i1+1,i2-1,i3)+8.*ry(i1+1,i2+1,i3)-ry(i1+1,i2+2,i3))/(12.*dr(1))-(ry(i1+2,i2-2,i3)-8.*ry(i1+2,i2-1,i3)+8.*ry(i1+2,i2+1,i3)-ry(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             ryss     = (-ry(i1,i2-2,i3)+16.*ry(i1,i2-1,i3)-30.*ry(i1,i2,i3)+16.*ry(i1,i2+1,i3)-ry(i1,i2+2,i3))/(12.*dr(1)**2)
                             ryt      = (ry(i1,i2,i3-2)-8.*ry(i1,i2,i3-1)+8.*ry(i1,i2,i3+1)-ry(i1,i2,i3+2))/(12.*dr(2))
                             ryrt     = ((ry(i1-2,i2,i3-2)-8.*ry(i1-2,i2,i3-1)+8.*ry(i1-2,i2,i3+1)-ry(i1-2,i2,i3+2))/(12.*dr(2))-8.*(ry(i1-1,i2,i3-2)-8.*ry(i1-1,i2,i3-1)+8.*ry(i1-1,i2,i3+1)-ry(i1-1,i2,i3+2))/(12.*dr(2))+8.*(ry(i1+1,i2,i3-2)-8.*ry(i1+1,i2,i3-1)+8.*ry(i1+1,i2,i3+1)-ry(i1+1,i2,i3+2))/(12.*dr(2))-(ry(i1+2,i2,i3-2)-8.*ry(i1+2,i2,i3-1)+8.*ry(i1+2,i2,i3+1)-ry(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             ryst     = ((ry(i1,i2-2,i3-2)-8.*ry(i1,i2-2,i3-1)+8.*ry(i1,i2-2,i3+1)-ry(i1,i2-2,i3+2))/(12.*dr(2))-8.*(ry(i1,i2-1,i3-2)-8.*ry(i1,i2-1,i3-1)+8.*ry(i1,i2-1,i3+1)-ry(i1,i2-1,i3+2))/(12.*dr(2))+8.*(ry(i1,i2+1,i3-2)-8.*ry(i1,i2+1,i3-1)+8.*ry(i1,i2+1,i3+1)-ry(i1,i2+1,i3+2))/(12.*dr(2))-(ry(i1,i2+2,i3-2)-8.*ry(i1,i2+2,i3-1)+8.*ry(i1,i2+2,i3+1)-ry(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             rytt     = (-ry(i1,i2,i3-2)+16.*ry(i1,i2,i3-1)-30.*ry(i1,i2,i3)+16.*ry(i1,i2,i3+1)-ry(i1,i2,i3+2))/(12.*dr(2)**2)
                             sxr      = (sx(i1-2,i2,i3)-8.*sx(i1-1,i2,i3)+8.*sx(i1+1,i2,i3)-sx(i1+2,i2,i3))/(12.*dr(0))
                             sxrr     = (-sx(i1-2,i2,i3)+16.*sx(i1-1,i2,i3)-30.*sx(i1,i2,i3)+16.*sx(i1+1,i2,i3)-sx(i1+2,i2,i3))/(12.*dr(0)**2)
                             sxs      = (sx(i1,i2-2,i3)-8.*sx(i1,i2-1,i3)+8.*sx(i1,i2+1,i3)-sx(i1,i2+2,i3))/(12.*dr(1))
                             sxrs     = ((sx(i1-2,i2-2,i3)-8.*sx(i1-2,i2-1,i3)+8.*sx(i1-2,i2+1,i3)-sx(i1-2,i2+2,i3))/(12.*dr(1))-8.*(sx(i1-1,i2-2,i3)-8.*sx(i1-1,i2-1,i3)+8.*sx(i1-1,i2+1,i3)-sx(i1-1,i2+2,i3))/(12.*dr(1))+8.*(sx(i1+1,i2-2,i3)-8.*sx(i1+1,i2-1,i3)+8.*sx(i1+1,i2+1,i3)-sx(i1+1,i2+2,i3))/(12.*dr(1))-(sx(i1+2,i2-2,i3)-8.*sx(i1+2,i2-1,i3)+8.*sx(i1+2,i2+1,i3)-sx(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             sxss     = (-sx(i1,i2-2,i3)+16.*sx(i1,i2-1,i3)-30.*sx(i1,i2,i3)+16.*sx(i1,i2+1,i3)-sx(i1,i2+2,i3))/(12.*dr(1)**2)
                             sxt      = (sx(i1,i2,i3-2)-8.*sx(i1,i2,i3-1)+8.*sx(i1,i2,i3+1)-sx(i1,i2,i3+2))/(12.*dr(2))
                             sxrt     = ((sx(i1-2,i2,i3-2)-8.*sx(i1-2,i2,i3-1)+8.*sx(i1-2,i2,i3+1)-sx(i1-2,i2,i3+2))/(12.*dr(2))-8.*(sx(i1-1,i2,i3-2)-8.*sx(i1-1,i2,i3-1)+8.*sx(i1-1,i2,i3+1)-sx(i1-1,i2,i3+2))/(12.*dr(2))+8.*(sx(i1+1,i2,i3-2)-8.*sx(i1+1,i2,i3-1)+8.*sx(i1+1,i2,i3+1)-sx(i1+1,i2,i3+2))/(12.*dr(2))-(sx(i1+2,i2,i3-2)-8.*sx(i1+2,i2,i3-1)+8.*sx(i1+2,i2,i3+1)-sx(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             sxst     = ((sx(i1,i2-2,i3-2)-8.*sx(i1,i2-2,i3-1)+8.*sx(i1,i2-2,i3+1)-sx(i1,i2-2,i3+2))/(12.*dr(2))-8.*(sx(i1,i2-1,i3-2)-8.*sx(i1,i2-1,i3-1)+8.*sx(i1,i2-1,i3+1)-sx(i1,i2-1,i3+2))/(12.*dr(2))+8.*(sx(i1,i2+1,i3-2)-8.*sx(i1,i2+1,i3-1)+8.*sx(i1,i2+1,i3+1)-sx(i1,i2+1,i3+2))/(12.*dr(2))-(sx(i1,i2+2,i3-2)-8.*sx(i1,i2+2,i3-1)+8.*sx(i1,i2+2,i3+1)-sx(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             sxtt     = (-sx(i1,i2,i3-2)+16.*sx(i1,i2,i3-1)-30.*sx(i1,i2,i3)+16.*sx(i1,i2,i3+1)-sx(i1,i2,i3+2))/(12.*dr(2)**2)
                             syr      = (sy(i1-2,i2,i3)-8.*sy(i1-1,i2,i3)+8.*sy(i1+1,i2,i3)-sy(i1+2,i2,i3))/(12.*dr(0))
                             syrr     = (-sy(i1-2,i2,i3)+16.*sy(i1-1,i2,i3)-30.*sy(i1,i2,i3)+16.*sy(i1+1,i2,i3)-sy(i1+2,i2,i3))/(12.*dr(0)**2)
                             sys      = (sy(i1,i2-2,i3)-8.*sy(i1,i2-1,i3)+8.*sy(i1,i2+1,i3)-sy(i1,i2+2,i3))/(12.*dr(1))
                             syrs     = ((sy(i1-2,i2-2,i3)-8.*sy(i1-2,i2-1,i3)+8.*sy(i1-2,i2+1,i3)-sy(i1-2,i2+2,i3))/(12.*dr(1))-8.*(sy(i1-1,i2-2,i3)-8.*sy(i1-1,i2-1,i3)+8.*sy(i1-1,i2+1,i3)-sy(i1-1,i2+2,i3))/(12.*dr(1))+8.*(sy(i1+1,i2-2,i3)-8.*sy(i1+1,i2-1,i3)+8.*sy(i1+1,i2+1,i3)-sy(i1+1,i2+2,i3))/(12.*dr(1))-(sy(i1+2,i2-2,i3)-8.*sy(i1+2,i2-1,i3)+8.*sy(i1+2,i2+1,i3)-sy(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             syss     = (-sy(i1,i2-2,i3)+16.*sy(i1,i2-1,i3)-30.*sy(i1,i2,i3)+16.*sy(i1,i2+1,i3)-sy(i1,i2+2,i3))/(12.*dr(1)**2)
                             syt      = (sy(i1,i2,i3-2)-8.*sy(i1,i2,i3-1)+8.*sy(i1,i2,i3+1)-sy(i1,i2,i3+2))/(12.*dr(2))
                             syrt     = ((sy(i1-2,i2,i3-2)-8.*sy(i1-2,i2,i3-1)+8.*sy(i1-2,i2,i3+1)-sy(i1-2,i2,i3+2))/(12.*dr(2))-8.*(sy(i1-1,i2,i3-2)-8.*sy(i1-1,i2,i3-1)+8.*sy(i1-1,i2,i3+1)-sy(i1-1,i2,i3+2))/(12.*dr(2))+8.*(sy(i1+1,i2,i3-2)-8.*sy(i1+1,i2,i3-1)+8.*sy(i1+1,i2,i3+1)-sy(i1+1,i2,i3+2))/(12.*dr(2))-(sy(i1+2,i2,i3-2)-8.*sy(i1+2,i2,i3-1)+8.*sy(i1+2,i2,i3+1)-sy(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             syst     = ((sy(i1,i2-2,i3-2)-8.*sy(i1,i2-2,i3-1)+8.*sy(i1,i2-2,i3+1)-sy(i1,i2-2,i3+2))/(12.*dr(2))-8.*(sy(i1,i2-1,i3-2)-8.*sy(i1,i2-1,i3-1)+8.*sy(i1,i2-1,i3+1)-sy(i1,i2-1,i3+2))/(12.*dr(2))+8.*(sy(i1,i2+1,i3-2)-8.*sy(i1,i2+1,i3-1)+8.*sy(i1,i2+1,i3+1)-sy(i1,i2+1,i3+2))/(12.*dr(2))-(sy(i1,i2+2,i3-2)-8.*sy(i1,i2+2,i3-1)+8.*sy(i1,i2+2,i3+1)-sy(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             sytt     = (-sy(i1,i2,i3-2)+16.*sy(i1,i2,i3-1)-30.*sy(i1,i2,i3)+16.*sy(i1,i2,i3+1)-sy(i1,i2,i3+2))/(12.*dr(2)**2)
                             rzr      = (rz(i1-2,i2,i3)-8.*rz(i1-1,i2,i3)+8.*rz(i1+1,i2,i3)-rz(i1+2,i2,i3))/(12.*dr(0))
                             rzrr     = (-rz(i1-2,i2,i3)+16.*rz(i1-1,i2,i3)-30.*rz(i1,i2,i3)+16.*rz(i1+1,i2,i3)-rz(i1+2,i2,i3))/(12.*dr(0)**2)
                             rzs      = (rz(i1,i2-2,i3)-8.*rz(i1,i2-1,i3)+8.*rz(i1,i2+1,i3)-rz(i1,i2+2,i3))/(12.*dr(1))
                             rzrs     = ((rz(i1-2,i2-2,i3)-8.*rz(i1-2,i2-1,i3)+8.*rz(i1-2,i2+1,i3)-rz(i1-2,i2+2,i3))/(12.*dr(1))-8.*(rz(i1-1,i2-2,i3)-8.*rz(i1-1,i2-1,i3)+8.*rz(i1-1,i2+1,i3)-rz(i1-1,i2+2,i3))/(12.*dr(1))+8.*(rz(i1+1,i2-2,i3)-8.*rz(i1+1,i2-1,i3)+8.*rz(i1+1,i2+1,i3)-rz(i1+1,i2+2,i3))/(12.*dr(1))-(rz(i1+2,i2-2,i3)-8.*rz(i1+2,i2-1,i3)+8.*rz(i1+2,i2+1,i3)-rz(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             rzss     = (-rz(i1,i2-2,i3)+16.*rz(i1,i2-1,i3)-30.*rz(i1,i2,i3)+16.*rz(i1,i2+1,i3)-rz(i1,i2+2,i3))/(12.*dr(1)**2)
                             rzt      = (rz(i1,i2,i3-2)-8.*rz(i1,i2,i3-1)+8.*rz(i1,i2,i3+1)-rz(i1,i2,i3+2))/(12.*dr(2))
                             rzrt     = ((rz(i1-2,i2,i3-2)-8.*rz(i1-2,i2,i3-1)+8.*rz(i1-2,i2,i3+1)-rz(i1-2,i2,i3+2))/(12.*dr(2))-8.*(rz(i1-1,i2,i3-2)-8.*rz(i1-1,i2,i3-1)+8.*rz(i1-1,i2,i3+1)-rz(i1-1,i2,i3+2))/(12.*dr(2))+8.*(rz(i1+1,i2,i3-2)-8.*rz(i1+1,i2,i3-1)+8.*rz(i1+1,i2,i3+1)-rz(i1+1,i2,i3+2))/(12.*dr(2))-(rz(i1+2,i2,i3-2)-8.*rz(i1+2,i2,i3-1)+8.*rz(i1+2,i2,i3+1)-rz(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             rzst     = ((rz(i1,i2-2,i3-2)-8.*rz(i1,i2-2,i3-1)+8.*rz(i1,i2-2,i3+1)-rz(i1,i2-2,i3+2))/(12.*dr(2))-8.*(rz(i1,i2-1,i3-2)-8.*rz(i1,i2-1,i3-1)+8.*rz(i1,i2-1,i3+1)-rz(i1,i2-1,i3+2))/(12.*dr(2))+8.*(rz(i1,i2+1,i3-2)-8.*rz(i1,i2+1,i3-1)+8.*rz(i1,i2+1,i3+1)-rz(i1,i2+1,i3+2))/(12.*dr(2))-(rz(i1,i2+2,i3-2)-8.*rz(i1,i2+2,i3-1)+8.*rz(i1,i2+2,i3+1)-rz(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             rztt     = (-rz(i1,i2,i3-2)+16.*rz(i1,i2,i3-1)-30.*rz(i1,i2,i3)+16.*rz(i1,i2,i3+1)-rz(i1,i2,i3+2))/(12.*dr(2)**2)
                             szr      = (sz(i1-2,i2,i3)-8.*sz(i1-1,i2,i3)+8.*sz(i1+1,i2,i3)-sz(i1+2,i2,i3))/(12.*dr(0))
                             szrr     = (-sz(i1-2,i2,i3)+16.*sz(i1-1,i2,i3)-30.*sz(i1,i2,i3)+16.*sz(i1+1,i2,i3)-sz(i1+2,i2,i3))/(12.*dr(0)**2)
                             szs      = (sz(i1,i2-2,i3)-8.*sz(i1,i2-1,i3)+8.*sz(i1,i2+1,i3)-sz(i1,i2+2,i3))/(12.*dr(1))
                             szrs     = ((sz(i1-2,i2-2,i3)-8.*sz(i1-2,i2-1,i3)+8.*sz(i1-2,i2+1,i3)-sz(i1-2,i2+2,i3))/(12.*dr(1))-8.*(sz(i1-1,i2-2,i3)-8.*sz(i1-1,i2-1,i3)+8.*sz(i1-1,i2+1,i3)-sz(i1-1,i2+2,i3))/(12.*dr(1))+8.*(sz(i1+1,i2-2,i3)-8.*sz(i1+1,i2-1,i3)+8.*sz(i1+1,i2+1,i3)-sz(i1+1,i2+2,i3))/(12.*dr(1))-(sz(i1+2,i2-2,i3)-8.*sz(i1+2,i2-1,i3)+8.*sz(i1+2,i2+1,i3)-sz(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             szss     = (-sz(i1,i2-2,i3)+16.*sz(i1,i2-1,i3)-30.*sz(i1,i2,i3)+16.*sz(i1,i2+1,i3)-sz(i1,i2+2,i3))/(12.*dr(1)**2)
                             szt      = (sz(i1,i2,i3-2)-8.*sz(i1,i2,i3-1)+8.*sz(i1,i2,i3+1)-sz(i1,i2,i3+2))/(12.*dr(2))
                             szrt     = ((sz(i1-2,i2,i3-2)-8.*sz(i1-2,i2,i3-1)+8.*sz(i1-2,i2,i3+1)-sz(i1-2,i2,i3+2))/(12.*dr(2))-8.*(sz(i1-1,i2,i3-2)-8.*sz(i1-1,i2,i3-1)+8.*sz(i1-1,i2,i3+1)-sz(i1-1,i2,i3+2))/(12.*dr(2))+8.*(sz(i1+1,i2,i3-2)-8.*sz(i1+1,i2,i3-1)+8.*sz(i1+1,i2,i3+1)-sz(i1+1,i2,i3+2))/(12.*dr(2))-(sz(i1+2,i2,i3-2)-8.*sz(i1+2,i2,i3-1)+8.*sz(i1+2,i2,i3+1)-sz(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             szst     = ((sz(i1,i2-2,i3-2)-8.*sz(i1,i2-2,i3-1)+8.*sz(i1,i2-2,i3+1)-sz(i1,i2-2,i3+2))/(12.*dr(2))-8.*(sz(i1,i2-1,i3-2)-8.*sz(i1,i2-1,i3-1)+8.*sz(i1,i2-1,i3+1)-sz(i1,i2-1,i3+2))/(12.*dr(2))+8.*(sz(i1,i2+1,i3-2)-8.*sz(i1,i2+1,i3-1)+8.*sz(i1,i2+1,i3+1)-sz(i1,i2+1,i3+2))/(12.*dr(2))-(sz(i1,i2+2,i3-2)-8.*sz(i1,i2+2,i3-1)+8.*sz(i1,i2+2,i3+1)-sz(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             sztt     = (-sz(i1,i2,i3-2)+16.*sz(i1,i2,i3-1)-30.*sz(i1,i2,i3)+16.*sz(i1,i2,i3+1)-sz(i1,i2,i3+2))/(12.*dr(2)**2)
                             txr      = (tx(i1-2,i2,i3)-8.*tx(i1-1,i2,i3)+8.*tx(i1+1,i2,i3)-tx(i1+2,i2,i3))/(12.*dr(0))
                             txrr     = (-tx(i1-2,i2,i3)+16.*tx(i1-1,i2,i3)-30.*tx(i1,i2,i3)+16.*tx(i1+1,i2,i3)-tx(i1+2,i2,i3))/(12.*dr(0)**2)
                             txs      = (tx(i1,i2-2,i3)-8.*tx(i1,i2-1,i3)+8.*tx(i1,i2+1,i3)-tx(i1,i2+2,i3))/(12.*dr(1))
                             txrs     = ((tx(i1-2,i2-2,i3)-8.*tx(i1-2,i2-1,i3)+8.*tx(i1-2,i2+1,i3)-tx(i1-2,i2+2,i3))/(12.*dr(1))-8.*(tx(i1-1,i2-2,i3)-8.*tx(i1-1,i2-1,i3)+8.*tx(i1-1,i2+1,i3)-tx(i1-1,i2+2,i3))/(12.*dr(1))+8.*(tx(i1+1,i2-2,i3)-8.*tx(i1+1,i2-1,i3)+8.*tx(i1+1,i2+1,i3)-tx(i1+1,i2+2,i3))/(12.*dr(1))-(tx(i1+2,i2-2,i3)-8.*tx(i1+2,i2-1,i3)+8.*tx(i1+2,i2+1,i3)-tx(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             txss     = (-tx(i1,i2-2,i3)+16.*tx(i1,i2-1,i3)-30.*tx(i1,i2,i3)+16.*tx(i1,i2+1,i3)-tx(i1,i2+2,i3))/(12.*dr(1)**2)
                             txt      = (tx(i1,i2,i3-2)-8.*tx(i1,i2,i3-1)+8.*tx(i1,i2,i3+1)-tx(i1,i2,i3+2))/(12.*dr(2))
                             txrt     = ((tx(i1-2,i2,i3-2)-8.*tx(i1-2,i2,i3-1)+8.*tx(i1-2,i2,i3+1)-tx(i1-2,i2,i3+2))/(12.*dr(2))-8.*(tx(i1-1,i2,i3-2)-8.*tx(i1-1,i2,i3-1)+8.*tx(i1-1,i2,i3+1)-tx(i1-1,i2,i3+2))/(12.*dr(2))+8.*(tx(i1+1,i2,i3-2)-8.*tx(i1+1,i2,i3-1)+8.*tx(i1+1,i2,i3+1)-tx(i1+1,i2,i3+2))/(12.*dr(2))-(tx(i1+2,i2,i3-2)-8.*tx(i1+2,i2,i3-1)+8.*tx(i1+2,i2,i3+1)-tx(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             txst     = ((tx(i1,i2-2,i3-2)-8.*tx(i1,i2-2,i3-1)+8.*tx(i1,i2-2,i3+1)-tx(i1,i2-2,i3+2))/(12.*dr(2))-8.*(tx(i1,i2-1,i3-2)-8.*tx(i1,i2-1,i3-1)+8.*tx(i1,i2-1,i3+1)-tx(i1,i2-1,i3+2))/(12.*dr(2))+8.*(tx(i1,i2+1,i3-2)-8.*tx(i1,i2+1,i3-1)+8.*tx(i1,i2+1,i3+1)-tx(i1,i2+1,i3+2))/(12.*dr(2))-(tx(i1,i2+2,i3-2)-8.*tx(i1,i2+2,i3-1)+8.*tx(i1,i2+2,i3+1)-tx(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             txtt     = (-tx(i1,i2,i3-2)+16.*tx(i1,i2,i3-1)-30.*tx(i1,i2,i3)+16.*tx(i1,i2,i3+1)-tx(i1,i2,i3+2))/(12.*dr(2)**2)
                             tyr      = (ty(i1-2,i2,i3)-8.*ty(i1-1,i2,i3)+8.*ty(i1+1,i2,i3)-ty(i1+2,i2,i3))/(12.*dr(0))
                             tyrr     = (-ty(i1-2,i2,i3)+16.*ty(i1-1,i2,i3)-30.*ty(i1,i2,i3)+16.*ty(i1+1,i2,i3)-ty(i1+2,i2,i3))/(12.*dr(0)**2)
                             tys      = (ty(i1,i2-2,i3)-8.*ty(i1,i2-1,i3)+8.*ty(i1,i2+1,i3)-ty(i1,i2+2,i3))/(12.*dr(1))
                             tyrs     = ((ty(i1-2,i2-2,i3)-8.*ty(i1-2,i2-1,i3)+8.*ty(i1-2,i2+1,i3)-ty(i1-2,i2+2,i3))/(12.*dr(1))-8.*(ty(i1-1,i2-2,i3)-8.*ty(i1-1,i2-1,i3)+8.*ty(i1-1,i2+1,i3)-ty(i1-1,i2+2,i3))/(12.*dr(1))+8.*(ty(i1+1,i2-2,i3)-8.*ty(i1+1,i2-1,i3)+8.*ty(i1+1,i2+1,i3)-ty(i1+1,i2+2,i3))/(12.*dr(1))-(ty(i1+2,i2-2,i3)-8.*ty(i1+2,i2-1,i3)+8.*ty(i1+2,i2+1,i3)-ty(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             tyss     = (-ty(i1,i2-2,i3)+16.*ty(i1,i2-1,i3)-30.*ty(i1,i2,i3)+16.*ty(i1,i2+1,i3)-ty(i1,i2+2,i3))/(12.*dr(1)**2)
                             tyt      = (ty(i1,i2,i3-2)-8.*ty(i1,i2,i3-1)+8.*ty(i1,i2,i3+1)-ty(i1,i2,i3+2))/(12.*dr(2))
                             tyrt     = ((ty(i1-2,i2,i3-2)-8.*ty(i1-2,i2,i3-1)+8.*ty(i1-2,i2,i3+1)-ty(i1-2,i2,i3+2))/(12.*dr(2))-8.*(ty(i1-1,i2,i3-2)-8.*ty(i1-1,i2,i3-1)+8.*ty(i1-1,i2,i3+1)-ty(i1-1,i2,i3+2))/(12.*dr(2))+8.*(ty(i1+1,i2,i3-2)-8.*ty(i1+1,i2,i3-1)+8.*ty(i1+1,i2,i3+1)-ty(i1+1,i2,i3+2))/(12.*dr(2))-(ty(i1+2,i2,i3-2)-8.*ty(i1+2,i2,i3-1)+8.*ty(i1+2,i2,i3+1)-ty(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             tyst     = ((ty(i1,i2-2,i3-2)-8.*ty(i1,i2-2,i3-1)+8.*ty(i1,i2-2,i3+1)-ty(i1,i2-2,i3+2))/(12.*dr(2))-8.*(ty(i1,i2-1,i3-2)-8.*ty(i1,i2-1,i3-1)+8.*ty(i1,i2-1,i3+1)-ty(i1,i2-1,i3+2))/(12.*dr(2))+8.*(ty(i1,i2+1,i3-2)-8.*ty(i1,i2+1,i3-1)+8.*ty(i1,i2+1,i3+1)-ty(i1,i2+1,i3+2))/(12.*dr(2))-(ty(i1,i2+2,i3-2)-8.*ty(i1,i2+2,i3-1)+8.*ty(i1,i2+2,i3+1)-ty(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             tytt     = (-ty(i1,i2,i3-2)+16.*ty(i1,i2,i3-1)-30.*ty(i1,i2,i3)+16.*ty(i1,i2,i3+1)-ty(i1,i2,i3+2))/(12.*dr(2)**2)
                             tzr      = (tz(i1-2,i2,i3)-8.*tz(i1-1,i2,i3)+8.*tz(i1+1,i2,i3)-tz(i1+2,i2,i3))/(12.*dr(0))
                             tzrr     = (-tz(i1-2,i2,i3)+16.*tz(i1-1,i2,i3)-30.*tz(i1,i2,i3)+16.*tz(i1+1,i2,i3)-tz(i1+2,i2,i3))/(12.*dr(0)**2)
                             tzs      = (tz(i1,i2-2,i3)-8.*tz(i1,i2-1,i3)+8.*tz(i1,i2+1,i3)-tz(i1,i2+2,i3))/(12.*dr(1))
                             tzrs     = ((tz(i1-2,i2-2,i3)-8.*tz(i1-2,i2-1,i3)+8.*tz(i1-2,i2+1,i3)-tz(i1-2,i2+2,i3))/(12.*dr(1))-8.*(tz(i1-1,i2-2,i3)-8.*tz(i1-1,i2-1,i3)+8.*tz(i1-1,i2+1,i3)-tz(i1-1,i2+2,i3))/(12.*dr(1))+8.*(tz(i1+1,i2-2,i3)-8.*tz(i1+1,i2-1,i3)+8.*tz(i1+1,i2+1,i3)-tz(i1+1,i2+2,i3))/(12.*dr(1))-(tz(i1+2,i2-2,i3)-8.*tz(i1+2,i2-1,i3)+8.*tz(i1+2,i2+1,i3)-tz(i1+2,i2+2,i3))/(12.*dr(1)))/(12.*dr(0))
                             tzss     = (-tz(i1,i2-2,i3)+16.*tz(i1,i2-1,i3)-30.*tz(i1,i2,i3)+16.*tz(i1,i2+1,i3)-tz(i1,i2+2,i3))/(12.*dr(1)**2)
                             tzt      = (tz(i1,i2,i3-2)-8.*tz(i1,i2,i3-1)+8.*tz(i1,i2,i3+1)-tz(i1,i2,i3+2))/(12.*dr(2))
                             tzrt     = ((tz(i1-2,i2,i3-2)-8.*tz(i1-2,i2,i3-1)+8.*tz(i1-2,i2,i3+1)-tz(i1-2,i2,i3+2))/(12.*dr(2))-8.*(tz(i1-1,i2,i3-2)-8.*tz(i1-1,i2,i3-1)+8.*tz(i1-1,i2,i3+1)-tz(i1-1,i2,i3+2))/(12.*dr(2))+8.*(tz(i1+1,i2,i3-2)-8.*tz(i1+1,i2,i3-1)+8.*tz(i1+1,i2,i3+1)-tz(i1+1,i2,i3+2))/(12.*dr(2))-(tz(i1+2,i2,i3-2)-8.*tz(i1+2,i2,i3-1)+8.*tz(i1+2,i2,i3+1)-tz(i1+2,i2,i3+2))/(12.*dr(2)))/(12.*dr(0))
                             tzst     = ((tz(i1,i2-2,i3-2)-8.*tz(i1,i2-2,i3-1)+8.*tz(i1,i2-2,i3+1)-tz(i1,i2-2,i3+2))/(12.*dr(2))-8.*(tz(i1,i2-1,i3-2)-8.*tz(i1,i2-1,i3-1)+8.*tz(i1,i2-1,i3+1)-tz(i1,i2-1,i3+2))/(12.*dr(2))+8.*(tz(i1,i2+1,i3-2)-8.*tz(i1,i2+1,i3-1)+8.*tz(i1,i2+1,i3+1)-tz(i1,i2+1,i3+2))/(12.*dr(2))-(tz(i1,i2+2,i3-2)-8.*tz(i1,i2+2,i3-1)+8.*tz(i1,i2+2,i3+1)-tz(i1,i2+2,i3+2))/(12.*dr(2)))/(12.*dr(1))
                             tztt     = (-tz(i1,i2,i3-2)+16.*tz(i1,i2,i3-1)-30.*tz(i1,i2,i3)+16.*tz(i1,i2,i3+1)-tz(i1,i2,i3+2))/(12.*dr(2)**2)
                             ! ---------- Spatial derivatives of metrics rx, sx, ry, ... ---------
                             rxi = rx(i1,i2,i3)
                             ryi = ry(i1,i2,i3)
                             sxi = sx(i1,i2,i3)
                             syi = sy(i1,i2,i3)
                             rzi = rz(i1,i2,i3)
                             szi = sz(i1,i2,i3)
                             txi = tx(i1,i2,i3)
                             tyi = ty(i1,i2,i3)
                             tzi = tz(i1,i2,i3)
                             rxx      = rxi*rxr+sxi*rxs+txi*rxt
                             rxy      = ryi*rxr+syi*rxs+tyi*rxt
                             rxz      = rzi*rxr+szi*rxs+tzi*rxt
                             ryy      = ryi*ryr+syi*rys+tyi*ryt
                             ryz      = rzi*ryr+szi*rys+tzi*ryt
                             rzz      = rzi*rzr+szi*rzs+tzi*rzt
                             sxx      = rxi*sxr+sxi*sxs+txi*sxt
                             sxy      = ryi*sxr+syi*sxs+tyi*sxt
                             sxz      = rzi*sxr+szi*sxs+tzi*sxt
                             syy      = ryi*syr+syi*sys+tyi*syt
                             syz      = rzi*syr+szi*sys+tzi*syt
                             szz      = rzi*szr+szi*szs+tzi*szt
                             txx      = rxi*txr+sxi*txs+txi*txt
                             txy      = ryi*txr+syi*txs+tyi*txt
                             txz      = rzi*txr+szi*txs+tzi*txt
                             tyy      = ryi*tyr+syi*tys+tyi*tyt
                             tyz      = rzi*tyr+szi*tys+tzi*tyt
                             tzz      = rzi*tzr+szi*tzs+tzi*tzt
                             rxxx     = rxi**2*rxrr+2.*rxi*sxi*rxrs+2.*rxi*txi*rxrt+sxi**2*rxss+2.*sxi*txi*rxst+txi**2*rxtt+rxx*rxr+sxx*rxs+txx*rxt
                             rxxy     = rxi*ryi*rxrr+(rxi*syi+ryi*sxi)*rxrs+sxi*syi*rxss+(rxi*tyi+ryi*txi)*rxrt+(sxi*tyi+syi*txi)*rxst+txi*tyi*rxtt+rxy*rxr+sxy*rxs+txy*rxt
                             rxyy     = ryi**2*rxrr+2.*ryi*syi*rxrs+2.*ryi*tyi*rxrt+syi**2*rxss+2.*syi*tyi*rxst+tyi**2*rxtt+ryy*rxr+syy*rxs+tyy*rxt
                             rxxz     = rxi*rzi*rxrr+(rxi*szi+rzi*sxi)*rxrs+sxi*szi*rxss+(rxi*tzi+rzi*txi)*rxrt+(sxi*tzi+szi*txi)*rxst+txi*tzi*rxtt+rxz*rxr+sxz*rxs+txz*rxt
                             rxyz     = ryi*rzi*rxrr+(ryi*szi+rzi*syi)*rxrs+syi*szi*rxss+(ryi*tzi+rzi*tyi)*rxrt+(syi*tzi+szi*tyi)*rxst+tyi*tzi*rxtt+ryz*rxr+syz*rxs+tyz*rxt
                             rxzz     = rzi**2*rxrr+2.*rzi*szi*rxrs+2.*rzi*tzi*rxrt+szi**2*rxss+2.*szi*tzi*rxst+tzi**2*rxtt+rzz*rxr+szz*rxs+tzz*rxt
                             ryyy     = ryi**2*ryrr+2.*ryi*syi*ryrs+2.*ryi*tyi*ryrt+syi**2*ryss+2.*syi*tyi*ryst+tyi**2*rytt+ryy*ryr+syy*rys+tyy*ryt
                             ryyz     = ryi*rzi*ryrr+(ryi*szi+rzi*syi)*ryrs+syi*szi*ryss+(ryi*tzi+rzi*tyi)*ryrt+(syi*tzi+szi*tyi)*ryst+tyi*tzi*rytt+ryz*ryr+syz*rys+tyz*ryt
                             ryzz     = rzi**2*ryrr+2.*rzi*szi*ryrs+2.*rzi*tzi*ryrt+szi**2*ryss+2.*szi*tzi*ryst+tzi**2*rytt+rzz*ryr+szz*rys+tzz*ryt
                             rzzz     = rzi**2*rzrr+2.*rzi*szi*rzrs+2.*rzi*tzi*rzrt+szi**2*rzss+2.*szi*tzi*rzst+tzi**2*rztt+rzz*rzr+szz*rzs+tzz*rzt
                             sxxx     = rxi**2*sxrr+2.*rxi*sxi*sxrs+2.*rxi*txi*sxrt+sxi**2*sxss+2.*sxi*txi*sxst+txi**2*sxtt+rxx*sxr+sxx*sxs+txx*sxt
                             sxxy     = rxi*ryi*sxrr+(rxi*syi+ryi*sxi)*sxrs+sxi*syi*sxss+(rxi*tyi+ryi*txi)*sxrt+(sxi*tyi+syi*txi)*sxst+txi*tyi*sxtt+rxy*sxr+sxy*sxs+txy*sxt
                             sxyy     = ryi**2*sxrr+2.*ryi*syi*sxrs+2.*ryi*tyi*sxrt+syi**2*sxss+2.*syi*tyi*sxst+tyi**2*sxtt+ryy*sxr+syy*sxs+tyy*sxt
                             sxxz     = rxi*rzi*sxrr+(rxi*szi+rzi*sxi)*sxrs+sxi*szi*sxss+(rxi*tzi+rzi*txi)*sxrt+(sxi*tzi+szi*txi)*sxst+txi*tzi*sxtt+rxz*sxr+sxz*sxs+txz*sxt
                             sxyz     = ryi*rzi*sxrr+(ryi*szi+rzi*syi)*sxrs+syi*szi*sxss+(ryi*tzi+rzi*tyi)*sxrt+(syi*tzi+szi*tyi)*sxst+tyi*tzi*sxtt+ryz*sxr+syz*sxs+tyz*sxt
                             sxzz     = rzi**2*sxrr+2.*rzi*szi*sxrs+2.*rzi*tzi*sxrt+szi**2*sxss+2.*szi*tzi*sxst+tzi**2*sxtt+rzz*sxr+szz*sxs+tzz*sxt
                             syyy     = ryi**2*syrr+2.*ryi*syi*syrs+2.*ryi*tyi*syrt+syi**2*syss+2.*syi*tyi*syst+tyi**2*sytt+ryy*syr+syy*sys+tyy*syt
                             syyz     = ryi*rzi*syrr+(ryi*szi+rzi*syi)*syrs+syi*szi*syss+(ryi*tzi+rzi*tyi)*syrt+(syi*tzi+szi*tyi)*syst+tyi*tzi*sytt+ryz*syr+syz*sys+tyz*syt
                             syzz     = rzi**2*syrr+2.*rzi*szi*syrs+2.*rzi*tzi*syrt+szi**2*syss+2.*szi*tzi*syst+tzi**2*sytt+rzz*syr+szz*sys+tzz*syt
                             szzz     = rzi**2*szrr+2.*rzi*szi*szrs+2.*rzi*tzi*szrt+szi**2*szss+2.*szi*tzi*szst+tzi**2*sztt+rzz*szr+szz*szs+tzz*szt
                             txxx     = rxi**2*txrr+2.*rxi*sxi*txrs+2.*rxi*txi*txrt+sxi**2*txss+2.*sxi*txi*txst+txi**2*txtt+rxx*txr+sxx*txs+txx*txt
                             txxy     = rxi*ryi*txrr+(rxi*syi+ryi*sxi)*txrs+sxi*syi*txss+(rxi*tyi+ryi*txi)*txrt+(sxi*tyi+syi*txi)*txst+txi*tyi*txtt+rxy*txr+sxy*txs+txy*txt
                             txyy     = ryi**2*txrr+2.*ryi*syi*txrs+2.*ryi*tyi*txrt+syi**2*txss+2.*syi*tyi*txst+tyi**2*txtt+ryy*txr+syy*txs+tyy*txt
                             txxz     = rxi*rzi*txrr+(rxi*szi+rzi*sxi)*txrs+sxi*szi*txss+(rxi*tzi+rzi*txi)*txrt+(sxi*tzi+szi*txi)*txst+txi*tzi*txtt+rxz*txr+sxz*txs+txz*txt
                             txyz     = ryi*rzi*txrr+(ryi*szi+rzi*syi)*txrs+syi*szi*txss+(ryi*tzi+rzi*tyi)*txrt+(syi*tzi+szi*tyi)*txst+tyi*tzi*txtt+ryz*txr+syz*txs+tyz*txt
                             txzz     = rzi**2*txrr+2.*rzi*szi*txrs+2.*rzi*tzi*txrt+szi**2*txss+2.*szi*tzi*txst+tzi**2*txtt+rzz*txr+szz*txs+tzz*txt
                             tyyy     = ryi**2*tyrr+2.*ryi*syi*tyrs+2.*ryi*tyi*tyrt+syi**2*tyss+2.*syi*tyi*tyst+tyi**2*tytt+ryy*tyr+syy*tys+tyy*tyt
                             tyyz     = ryi*rzi*tyrr+(ryi*szi+rzi*syi)*tyrs+syi*szi*tyss+(ryi*tzi+rzi*tyi)*tyrt+(syi*tzi+szi*tyi)*tyst+tyi*tzi*tytt+ryz*tyr+syz*tys+tyz*tyt
                             tyzz     = rzi**2*tyrr+2.*rzi*szi*tyrs+2.*rzi*tzi*tyrt+szi**2*tyss+2.*szi*tzi*tyst+tzi**2*tytt+rzz*tyr+szz*tys+tzz*tyt
                             tzzz     = rzi**2*tzrr+2.*rzi*szi*tzrs+2.*rzi*tzi*tzrt+szi**2*tzss+2.*szi*tzi*tzst+tzi**2*tztt+rzz*tzr+szz*tzs+tzz*tzt
                             ! ---- end evalMetrics eq evalMetrics ---
                             ! ---------- Third spatial derivatives of u ---------
                             uxyz     = rxi*ryi*rzi*urrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*urrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*urss+sxi*syi*szi*usss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*urrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*urst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*usst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*urtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*ustt+txi*tyi*tzi*uttt+(rxi*ryz+rxy*rzi+rxz*ryi)*urr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*urs+(sxi*syz+sxy*szi+sxz*syi)*uss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*urt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*ust+(txi*tyz+txy*tzi+txz*tyi)*utt+rxyz*ur+sxyz*us+txyz*ut
                             ! ---------- END CURVILINEAR  ---------
                             ! ------ curvilinear grid: -------
                               ! 3D 
                               r1 =  an1*ux43(i1,i2,i3,0) + an2*uy43(i1,i2,i3,0) + an3*uz43(i1,i2,i3,0) - gg
                               crv(axis) = an1*rsxy(i1,i2,i3,axis,0) + an2*rsxy(i1,i2,i3,axis,1) + an3*rsxy(i1,i2,i3,axis,2)
                               ! **CHECK ME**
                               a11 = -is*( crv(axis)*8./(12.*dr(axis)) )  ! coeff of u(-1)
                               a12 =  is*( crv(axis)*1./(12.*dr(axis)) )  ! coeff of u(-2)
                               r2 = c2*( an1*( uxxx + uxyy + uxzz ) + an2*( uxxy + uyyy + uyzz ) + an3*( uxxz + uyyz + uzzz ) ) + nDotGradF - gtt
                               ! crv(axis) = an1*rsxy(i1,i2,i3,axis,0)**3 + an2*rsxy(i1,i2,i3,axis,1)**3 + an3*rsxy(i1,i2,i3,axis,2)**3
                               ! Coeff of "urrr" term *check me*
                               crv(axis) = (an1*rsxy(i1,i2,i3,axis,0) + an2*rsxy(i1,i2,i3,axis,1) + an3*rsxy(i1,i2,i3,axis,2) )* (    rsxy(i1,i2,i3,axis,0)**2  + rsxy(i1,i2,i3,axis,1)**2  + rsxy(i1,i2,i3,axis,2)**2 )
                               ! if( axis.eq.0 )then
                               !   crv(axis) = ( an1*rxi + an2*ryi + an3*rzi )*( rxi**2 + ryi**2 + rzi**2 ) 
                               !   ! crv(axis) = an1*( rxi*( rxi**2 + ryi**2 + rzi**2 ) ) + !   !             an2*( ryi*( rxi**2 + ryi**2 + rzi**2 ) ) + !   !             an3*( rzi*( rxi**2 + ryi**2 + rzi**2 ) ) 
                               ! else if( axis.eq.1 )then
                               !   crv(axis) = ( an1*sxi + an2*syi + an3*szi )*( sxi**2 + syi**2 + szi**2 ) 
                               ! else           
                               !   crv(axis) = ( an1*txi + an2*tyi + an3*tzi )*( txi**2 + tyi**2 + tzi**2 ) 
                               ! end if
                               ! **CHECK ME**
                               a21 =  is*c2*( 2.*crv(axis)/(2.*dr(axis)**3) )
                               a22 = -is*c2*( 1.*crv(axis)/(2.*dr(axis)**3) )            
                               ! define the residual functions for the discrete delta method
                                 if( checkCoeff.eq.1 )then
                                   ! --- Check coefficients a11,a12,... by discrete delta ---
                                   u1Save = u(j1,j2,j3,0)
                                   u2Save = u(k1,k2,k3,0)
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=0.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ut       = (u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))
                                   urt      = ((u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0))
                                   urrt     = (-(u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0)**2)
                                   ust      = ((u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))
                                   urst     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   usst     = (-(u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)**2)
                                   utt      = (-u(i1,i2,i3-2,0)+16.*u(i1,i2,i3-1,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2)**2)
                                   urtt     = ((-u(i1-2,i2,i3-2,0)+16.*u(i1-2,i2,i3-1,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1-1,i2,i3-2,0)+16.*u(i1-1,i2,i3-1,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1+1,i2,i3-2,0)+16.*u(i1+1,i2,i3-1,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2)**2)-(-u(i1+2,i2,i3-2,0)+16.*u(i1+2,i2,i3-1,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(0))
                                   ustt     = ((-u(i1,i2-2,i3-2,0)+16.*u(i1,i2-2,i3-1,0)-30.*u(i1,i2-2,i3,0)+16.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1,i2-1,i3-2,0)+16.*u(i1,i2-1,i3-1,0)-30.*u(i1,i2-1,i3,0)+16.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1,i2+1,i3-2,0)+16.*u(i1,i2+1,i3-1,0)-30.*u(i1,i2+1,i3,0)+16.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2)**2)-(-u(i1,i2+2,i3-2,0)+16.*u(i1,i2+2,i3-1,0)-30.*u(i1,i2+2,i3,0)+16.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(1))
                                   uttt     = (-u(i1,i2,i3-2,0)+2.*u(i1,i2,i3-1,0)-2.*u(i1,i2,i3+1,0)+u(i1,i2,i3+2,0))/(2.*dr(2)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   uxyz     = rxi*ryi*rzi*urrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*urrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*urss+sxi*syi*szi*usss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*urrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*urst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*usst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*urtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*ustt+txi*tyi*tzi*uttt+(rxi*ryz+rxy*rzi+rxz*ryi)*urr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*urs+(sxi*syz+sxy*szi+sxz*syi)*uss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*urt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*ust+(txi*tyz+txy*tzi+txz*tyi)*utt+rxyz*ur+sxyz*us+txyz*ut
                                   ! ---------- END CURVILINEAR  ---------
                                   r1a = (an1*ux43(i1,i2,i3,0)+an2*uy43(i1,i2,i3,0)+an3*uz43(i1,i2,i3,0))
                                   r2a = (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   u(j1,j2,j3,0)=1.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ut       = (u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))
                                   urt      = ((u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0))
                                   urrt     = (-(u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0)**2)
                                   ust      = ((u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))
                                   urst     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   usst     = (-(u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)**2)
                                   utt      = (-u(i1,i2,i3-2,0)+16.*u(i1,i2,i3-1,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2)**2)
                                   urtt     = ((-u(i1-2,i2,i3-2,0)+16.*u(i1-2,i2,i3-1,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1-1,i2,i3-2,0)+16.*u(i1-1,i2,i3-1,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1+1,i2,i3-2,0)+16.*u(i1+1,i2,i3-1,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2)**2)-(-u(i1+2,i2,i3-2,0)+16.*u(i1+2,i2,i3-1,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(0))
                                   ustt     = ((-u(i1,i2-2,i3-2,0)+16.*u(i1,i2-2,i3-1,0)-30.*u(i1,i2-2,i3,0)+16.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1,i2-1,i3-2,0)+16.*u(i1,i2-1,i3-1,0)-30.*u(i1,i2-1,i3,0)+16.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1,i2+1,i3-2,0)+16.*u(i1,i2+1,i3-1,0)-30.*u(i1,i2+1,i3,0)+16.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2)**2)-(-u(i1,i2+2,i3-2,0)+16.*u(i1,i2+2,i3-1,0)-30.*u(i1,i2+2,i3,0)+16.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(1))
                                   uttt     = (-u(i1,i2,i3-2,0)+2.*u(i1,i2,i3-1,0)-2.*u(i1,i2,i3+1,0)+u(i1,i2,i3+2,0))/(2.*dr(2)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   uxyz     = rxi*ryi*rzi*urrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*urrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*urss+sxi*syi*szi*usss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*urrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*urst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*usst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*urtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*ustt+txi*tyi*tzi*uttt+(rxi*ryz+rxy*rzi+rxz*ryi)*urr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*urs+(sxi*syz+sxy*szi+sxz*syi)*uss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*urt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*ust+(txi*tyz+txy*tzi+txz*tyi)*utt+rxyz*ur+sxyz*us+txyz*ut
                                   ! ---------- END CURVILINEAR  ---------
                                   r1b =  (an1*ux43(i1,i2,i3,0)+an2*uy43(i1,i2,i3,0)+an3*uz43(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   a11c = r1b - r1a
                                   a21c = r2b - r2a
                                   u(j1,j2,j3,0)=0.
                                   u(k1,k2,k3,0)=1.
                                   ! ---------- START CURVILINEAR  ---------
                                   ! ---------- Parametric derivatives ---------
                                   ur       = (u(i1-2,i2,i3,0)-8.*u(i1-1,i2,i3,0)+8.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0))
                                   urr      = (-u(i1-2,i2,i3,0)+16.*u(i1-1,i2,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1+1,i2,i3,0)-u(i1+2,i2,i3,0))/(12.*dr(0)**2)
                                   urrr     = (-u(i1-2,i2,i3,0)+2.*u(i1-1,i2,i3,0)-2.*u(i1+1,i2,i3,0)+u(i1+2,i2,i3,0))/(2.*dr(0)**3)
                                   us       = (u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))
                                   urs      = ((u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))-8.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))+8.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0))
                                   urrs     = (-(u(i1-2,i2-2,i3,0)-8.*u(i1-2,i2-1,i3,0)+8.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1-1,i2-2,i3,0)-8.*u(i1-1,i2-1,i3,0)+8.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1))-30.*(u(i1,i2-2,i3,0)-8.*u(i1,i2-1,i3,0)+8.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1))+16.*(u(i1+1,i2-2,i3,0)-8.*u(i1+1,i2-1,i3,0)+8.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1))-(u(i1+2,i2-2,i3,0)-8.*u(i1+2,i2-1,i3,0)+8.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)))/(12.*dr(0)**2)
                                   uss      = (-u(i1,i2-2,i3,0)+16.*u(i1,i2-1,i3,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2+1,i3,0)-u(i1,i2+2,i3,0))/(12.*dr(1)**2)
                                   urss     = ((-u(i1-2,i2-2,i3,0)+16.*u(i1-2,i2-1,i3,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2+1,i3,0)-u(i1-2,i2+2,i3,0))/(12.*dr(1)**2)-8.*(-u(i1-1,i2-2,i3,0)+16.*u(i1-1,i2-1,i3,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2+1,i3,0)-u(i1-1,i2+2,i3,0))/(12.*dr(1)**2)+8.*(-u(i1+1,i2-2,i3,0)+16.*u(i1+1,i2-1,i3,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2+1,i3,0)-u(i1+1,i2+2,i3,0))/(12.*dr(1)**2)-(-u(i1+2,i2-2,i3,0)+16.*u(i1+2,i2-1,i3,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2+1,i3,0)-u(i1+2,i2+2,i3,0))/(12.*dr(1)**2))/(12.*dr(0))
                                   usss     = (-u(i1,i2-2,i3,0)+2.*u(i1,i2-1,i3,0)-2.*u(i1,i2+1,i3,0)+u(i1,i2+2,i3,0))/(2.*dr(1)**3)
                                   ut       = (u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))
                                   urt      = ((u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0))
                                   urrt     = (-(u(i1-2,i2,i3-2,0)-8.*u(i1-2,i2,i3-1,0)+8.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1-1,i2,i3-2,0)-8.*u(i1-1,i2,i3-1,0)+8.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1+1,i2,i3-2,0)-8.*u(i1+1,i2,i3-1,0)+8.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2))-(u(i1+2,i2,i3-2,0)-8.*u(i1+2,i2,i3-1,0)+8.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)))/(12.*dr(0)**2)
                                   ust      = ((u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))
                                   urst     = (((u(i1-2,i2-2,i3-2,0)-8.*u(i1-2,i2-2,i3-1,0)+8.*u(i1-2,i2-2,i3+1,0)-u(i1-2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-2,i2-1,i3-2,0)-8.*u(i1-2,i2-1,i3-1,0)+8.*u(i1-2,i2-1,i3+1,0)-u(i1-2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-2,i2+1,i3-2,0)-8.*u(i1-2,i2+1,i3-1,0)+8.*u(i1-2,i2+1,i3+1,0)-u(i1-2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-2,i2+2,i3-2,0)-8.*u(i1-2,i2+2,i3-1,0)+8.*u(i1-2,i2+2,i3+1,0)-u(i1-2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-8.*((u(i1-1,i2-2,i3-2,0)-8.*u(i1-1,i2-2,i3-1,0)+8.*u(i1-1,i2-2,i3+1,0)-u(i1-1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1-1,i2-1,i3-2,0)-8.*u(i1-1,i2-1,i3-1,0)+8.*u(i1-1,i2-1,i3+1,0)-u(i1-1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1-1,i2+1,i3-2,0)-8.*u(i1-1,i2+1,i3-1,0)+8.*u(i1-1,i2+1,i3+1,0)-u(i1-1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1-1,i2+2,i3-2,0)-8.*u(i1-1,i2+2,i3-1,0)+8.*u(i1-1,i2+2,i3+1,0)-u(i1-1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))+8.*((u(i1+1,i2-2,i3-2,0)-8.*u(i1+1,i2-2,i3-1,0)+8.*u(i1+1,i2-2,i3+1,0)-u(i1+1,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+1,i2-1,i3-2,0)-8.*u(i1+1,i2-1,i3-1,0)+8.*u(i1+1,i2-1,i3+1,0)-u(i1+1,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+1,i2+1,i3-2,0)-8.*u(i1+1,i2+1,i3-1,0)+8.*u(i1+1,i2+1,i3+1,0)-u(i1+1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+1,i2+2,i3-2,0)-8.*u(i1+1,i2+2,i3-1,0)+8.*u(i1+1,i2+2,i3+1,0)-u(i1+1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1))-((u(i1+2,i2-2,i3-2,0)-8.*u(i1+2,i2-2,i3-1,0)+8.*u(i1+2,i2-2,i3+1,0)-u(i1+2,i2-2,i3+2,0))/(12.*dr(2))-8.*(u(i1+2,i2-1,i3-2,0)-8.*u(i1+2,i2-1,i3-1,0)+8.*u(i1+2,i2-1,i3+1,0)-u(i1+2,i2-1,i3+2,0))/(12.*dr(2))+8.*(u(i1+2,i2+1,i3-2,0)-8.*u(i1+2,i2+1,i3-1,0)+8.*u(i1+2,i2+1,i3+1,0)-u(i1+2,i2+1,i3+2,0))/(12.*dr(2))-(u(i1+2,i2+2,i3-2,0)-8.*u(i1+2,i2+2,i3-1,0)+8.*u(i1+2,i2+2,i3+1,0)-u(i1+2,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)))/(12.*dr(0))
                                   usst     = (-(u(i1,i2-2,i3-2,0)-8.*u(i1,i2-2,i3-1,0)+8.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2-1,i3-2,0)-8.*u(i1,i2-1,i3-1,0)+8.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2))-30.*(u(i1,i2,i3-2,0)-8.*u(i1,i2,i3-1,0)+8.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2))+16.*(u(i1,i2+1,i3-2,0)-8.*u(i1,i2+1,i3-1,0)+8.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2))-(u(i1,i2+2,i3-2,0)-8.*u(i1,i2+2,i3-1,0)+8.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)))/(12.*dr(1)**2)
                                   utt      = (-u(i1,i2,i3-2,0)+16.*u(i1,i2,i3-1,0)-30.*u(i1,i2,i3,0)+16.*u(i1,i2,i3+1,0)-u(i1,i2,i3+2,0))/(12.*dr(2)**2)
                                   urtt     = ((-u(i1-2,i2,i3-2,0)+16.*u(i1-2,i2,i3-1,0)-30.*u(i1-2,i2,i3,0)+16.*u(i1-2,i2,i3+1,0)-u(i1-2,i2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1-1,i2,i3-2,0)+16.*u(i1-1,i2,i3-1,0)-30.*u(i1-1,i2,i3,0)+16.*u(i1-1,i2,i3+1,0)-u(i1-1,i2,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1+1,i2,i3-2,0)+16.*u(i1+1,i2,i3-1,0)-30.*u(i1+1,i2,i3,0)+16.*u(i1+1,i2,i3+1,0)-u(i1+1,i2,i3+2,0))/(12.*dr(2)**2)-(-u(i1+2,i2,i3-2,0)+16.*u(i1+2,i2,i3-1,0)-30.*u(i1+2,i2,i3,0)+16.*u(i1+2,i2,i3+1,0)-u(i1+2,i2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(0))
                                   ustt     = ((-u(i1,i2-2,i3-2,0)+16.*u(i1,i2-2,i3-1,0)-30.*u(i1,i2-2,i3,0)+16.*u(i1,i2-2,i3+1,0)-u(i1,i2-2,i3+2,0))/(12.*dr(2)**2)-8.*(-u(i1,i2-1,i3-2,0)+16.*u(i1,i2-1,i3-1,0)-30.*u(i1,i2-1,i3,0)+16.*u(i1,i2-1,i3+1,0)-u(i1,i2-1,i3+2,0))/(12.*dr(2)**2)+8.*(-u(i1,i2+1,i3-2,0)+16.*u(i1,i2+1,i3-1,0)-30.*u(i1,i2+1,i3,0)+16.*u(i1,i2+1,i3+1,0)-u(i1,i2+1,i3+2,0))/(12.*dr(2)**2)-(-u(i1,i2+2,i3-2,0)+16.*u(i1,i2+2,i3-1,0)-30.*u(i1,i2+2,i3,0)+16.*u(i1,i2+2,i3+1,0)-u(i1,i2+2,i3+2,0))/(12.*dr(2)**2))/(12.*dr(1))
                                   uttt     = (-u(i1,i2,i3-2,0)+2.*u(i1,i2,i3-1,0)-2.*u(i1,i2,i3+1,0)+u(i1,i2,i3+2,0))/(2.*dr(2)**3)
                                   ! ---- end nothing eq evalMetrics ---
                                   ! ---------- Third spatial derivatives of u ---------
                                   uxyz     = rxi*ryi*rzi*urrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*urrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*urss+sxi*syi*szi*usss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*urrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*urst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*usst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*urtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*ustt+txi*tyi*tzi*uttt+(rxi*ryz+rxy*rzi+rxz*ryi)*urr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*urs+(sxi*syz+sxy*szi+sxz*syi)*uss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*urt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*ust+(txi*tyz+txy*tzi+txz*tyi)*utt+rxyz*ur+sxyz*us+txyz*ut
                                   ! ---------- END CURVILINEAR  ---------
                                   r1b =  (an1*ux43(i1,i2,i3,0)+an2*uy43(i1,i2,i3,0)+an3*uz43(i1,i2,i3,0))
                                   r2b =  (c2*(an1*(uxxx+uxyy+uxzz)+an2*(uxxy+uyyy+uyzz)+an3*(uxxz+uyyz+uzzz)))
                                   a12c = r1b - r1a
                                   a22c = r2b - r2a
                                   u(j1,j2,j3,0) = u1Save  
                                   u(k1,k2,k3,0) = u2Save      
                                   maxDiff=max(maxDiff,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c),abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c))
                                   ! #If "curvilinear" eq "curvilinear"
                                   ! write(*,'("++ i1,i2=",2i4," a11,a11c=",2(1pe12.4)," a12,a12c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a11,a11c,a12,a12c,abs(a11-a11c)/abs(a11c),abs(a12-a12c)/abs(a12c)
                                   ! write(*,'("++ i1,i2=",2i4," a21,a21c=",2(1pe12.4)," a22,a22c=",2(1pe12.4)," rel-diff=",2(e9.2))') i1,i2,a21,a21c,a22,a22c,abs(a21-a21c)/abs(a21c),abs(a22-a22c)/abs(a22c)
                                   ! #End
                                 end if 
                           f1 = a11*u(j1,j2,j3,0) + a12*u(k1,k2,k3,0) - r1
                           f2 = a21*u(j1,j2,j3,0) + a22*u(k1,k2,k3,0) - r2
                           det = a11*a22 - a21*a12
                           uTemp(j1,j2,j3,0) = ( a22*f1 - a12*f2 )/det
                           uTemp(k1,k2,k3,0) = ( a11*f2 - a21*f1 )/det
                           if( .false. )then
                             call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,uex )
                             call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,uey )
                             write(*,'(" (i1,i2)=(",2i3,") u(-1),ue(-1)=",2(1pe12.4)," u(-2),ue(-2)=",2(1pe12.4))') i1,i2,uTemp(j1,j2,j3,0),uex,uTemp(k1,k2,k3,0),uey
                           end if
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! this case done below 
                         end if
                        end do
                        end do
                        end do
                       ! ------ fill in ghost values from uTemp ----
                        do i3=n3a,n3b
                        do i2=n2a,n2b
                        do i1=n1a,n1b
                         if( mask(i1,i2,i3).gt.0 )then
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost 
                           u(j1,j2,j3,0) = uTemp(j1,j2,j3,0)
                           u(k1,k2,k3,0) = uTemp(k1,k2,k3,0) 
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if        
                         else if( mask(i1,i2,i3).lt.0 )then
                           ! ----- extrap ghost outside interp. pts on physical boundaries ------
                           ghost = 1 
                           j1=i1-is1*ghost
                           j2=i2-is2*ghost
                           j3=i3-is3*ghost
                           u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                           ghost = 2
                           k1=i1-is1*ghost
                           k2=i2-is2*ghost
                           k3=i3-is3*ghost           
                           u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                           if( numGhost.gt.2 )then
                             !  extrap third ghost (UPW)
                             ghost = 3
                             l1=i1-is1*ghost
                             l2=i2-is2*ghost
                             l3=i3-is3*ghost           
                             u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                           end if
                         end if    
                        end do
                        end do
                        end do
                       if( checkCoeff.eq.1 )then
                         write(*,'("bcOptWave: t=",e12.3," checkCoeff: (side,axis)=(",2i2,") rel-maxDiff=",e9.2," curvilinear 4")') t,side,axis,maxDiff
                       end if
                 end if
               end if     
             end if
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
         else if( orderOfAccuracy.eq.4 )then
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
                 ! 4=4: 
                 ! --- assign 2 ghost points using:
                 !  (1) Apply Neumann BC to 4th order
                 !  (2) Extrap. 2nd ghost to 5th order
                 if( gridType.eq.rectangular )then
                   ! write(*,'(" TBC: j1,j2=",2i3," u,ff=",2e12.2)') j1,j2,ff,u(j1,j2,j3,uc)
                   !if( orderOfAccuracy.eq.2 )then
                    !  --- NEUMANN 4=4 rectangular ---
                     u(j1,j2,j3,uc)=  b0*u(j1+  is1,j2+  is2,j3+  is3,uc)+6.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-2.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc)/3.+b1*ff
                 else 
                   ! ------ curvilinear grid: -------
                   ! a1*( n1*ux + n2*ux + n3*uz ) + a0*u = f 
                   ! a1*( (n1*rx+n2*ry+n3*rz)*ur + (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ) + a0*u = f 
                   ! =>
                   !  ur = [ f - (n1*sx+n2*sy+n3*sz)*us + (n1*tx+n2*ty+n3*st)*ut ]/( a1*( (n1*rx+n2*ry+n3*rz) ) 
                     ! ---- NEUMANN 4 4 curvilinear ----
                     !       d14(kd) = 1./(12.*dr(kd))
                     !       ur4(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(i1+2,
                     !        & i2,i3,kd)-u(i1-2,i2,i3,kd)))*d14(0)
                     ! ur = f -> 
                     ! u(-2) -8*u(-1) =            -8*u(1)   + u(2)        + 12*dr( f )    --- (A)
                     ! u(-2) -5*u(-1) = -10*u(0) + 10*u(1) - 5*u(2) + u(3)                 --- (B)
                     ! A - B = 
                     !       -3*u(-1) =  10*u(0) - 18*u(1) + 6*u(2) - u(3) + 12*dr*( f ) 
                     urv(0) = ur4(i1,i2,i3,uc)
                     urv(1) = us4(i1,i2,i3,uc)
                     if( nd.eq.2 )then
                       t1=a1*( an1*rsxy(i1,i2,i3,axis  ,0)+an2*rsxy(i1,i2,i3,axis  ,1) )
                       t2=a1*( an1*rsxy(i1,i2,i3,axisp1,0)+an2*rsxy(i1,i2,i3,axisp1,1) )
                       ur0 = (ff - ( t2*urv(axisp1) + a0*u(i1,i2,i3,uc) ) )/t1
                     else
                       urv(2) = ut4(i1,i2,i3,uc)
                       t1=a1*( an1*rsxy(i1,i2,i3,axis  ,0)+an2*rsxy(i1,i2,i3,axis  ,1)+an3*rsxy(i1,i2,i3,axis  ,2) )
                       t2=a1*( an1*rsxy(i1,i2,i3,axisp1,0)+an2*rsxy(i1,i2,i3,axisp1,1)+an3*rsxy(i1,i2,i3,axisp1,2) )
                       t3=a1*( an1*rsxy(i1,i2,i3,axisp2,0)+an2*rsxy(i1,i2,i3,axisp2,1)+an3*rsxy(i1,i2,i3,axisp2,2) )
                       ur0 = ( ff - ( t2*urv(axisp1) + t3*urv(axisp2) + a0*u(i1,i2,i3,uc) ) )/t1
                     end if
                     u(j1,j2,j3,uc) = (-10./3.)*u(j1+  is1,j2+  is2,j3+  is3,uc)+6.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-2.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+(1./3)*u(j1+4*is1,j2+4*is2,j3+4*is3,uc)-4.*is*dr(axis)*ur0
                 end if ! curvilinear grid 
                 ! ----- Assign extra ghost ----
                   ! For Neumann BC's it IS necessary to extrap to order 5 for fourth order. 
                   ! extrap second ghost
                   u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                   if( numGhost.gt.2 )then
                     !  extrap third ghost (UPW)
                     u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                   end if
               else if( mask(i1,i2,i3).lt.0 )then
                 ! ----- extrap ghost outside interp. pts on physical boundaries ------
                   u(j1,j2,j3,uc)=(5.*u(i1,i2,i3,uc)-10.*u(i1+is1,i2+is2,i3+is3,uc)+10.*u(i1+2*is1,i2+2*is2,i3+2*is3,uc)-5.*u(i1+3*is1,i2+3*is2,i3+3*is3,uc)+u(i1+4*is1,i2+4*is2,i3+4*is3,uc))
                   u(k1,k2,k3,uc)=(5.*u(j1,j2,j3,uc)-10.*u(j1+is1,j2+is2,j3+is3,uc)+10.*u(j1+2*is1,j2+2*is2,j3+2*is3,uc)-5.*u(j1+3*is1,j2+3*is2,j3+3*is3,uc)+u(j1+4*is1,j2+4*is2,j3+4*is3,uc))
                   if( numGhost.gt.2 )then
                     !  extrap third ghost 
                     u(l1,l2,l3,uc)=(5.*u(k1,k2,k3,uc)-10.*u(k1+is1,k2+is2,k3+is3,uc)+10.*u(k1+2*is1,k2+2*is2,k3+2*is3,uc)-5.*u(k1+3*is1,k2+3*is2,k3+3*is3,uc)+u(k1+4*is1,k2+4*is2,k3+4*is3,uc))            
                   end if
               end if
               end do
               end do
               end do
         else if( orderOfAccuracy.eq.6 )then
         else if( orderOfAccuracy.eq.8 )then   
         else
           write(*,'("CgWave::bcOpt:ERROR: unexpected orderOfAccuracy=",i6)') orderOfAccuracy
           stop 8888
         end if
       end if
     else if(  bc(side,axis).eq.dirichlet .or. bc(side,axis).eq.exactBC .or. bc(side,axis).le.0 )then
       ! do nothing
     else
       write(*,'("bcOptWave: unexpected boundaryCondition=",i4)') bc(side,axis)
       stop 5151
     end if 
    end do ! end side
    end do ! end axis
   !  --- Assign ghost points outside corners ---
   if( .false. .and. orderOfAccuracy.eq.6 .and. nd.eq.2 .and. gridType.eq.rectangular )then
     if( forcingOption.eq.noForcing )then
         ! ---------------------------------
         ! --- assign corners and edges: ---
         ! ---------------------------------
         if( nd.eq.2 )then
           ! --- TWO DIMENSIONS ----
           do side2=0,1
           do side1=0,1
             if( bc(side1,0).gt.0 .and. bc(side2,1).gt.0 )then
               ! **fix me** for exact BCs
               if( ( bc(side1,0).ne.dirichlet .and. bc(side1,0).ne.neumann ) .or. ( bc(side2,1).ne.dirichlet .and. bc(side2,1).ne.neumann ) )then
                 write(*,*) "Un-supported corner bcs =",bc(side1,0),bc(side2,1)
                 stop 2222
               end if
               symSign = +1. ! even symmetry for D-D or N-N corners
               if( (bc(side1,0).eq.dirichlet .and. bc(side2,1).eq.neumann   ) .or. (bc(side1,0).eq.neumann   .and. bc(side2,1).eq.dirichlet ) ) then
                 symSign=-1.;
               end if
               is1 = 1-2*side1
               is2 = 1-2*side2
               i1 = gridIndexRange(side1,0)
               i2 = gridIndexRange(side2,1)
               i3 = gridIndexRange(    0,2)
               do m2=0,numGhost
               do m1=0,numGhost
                 j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                 k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point 
                 u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0)
               end do
               end do
             end if 
           end do
           end do
         else
           ! --- THREE DIMENSIONS ----
           write(*,*) "symmetry corners -- finish me in 3D"
           stop 9999
           do side3=0,1
           do side2=0,1
           do side1=0,1
           end do
           end do
           end do
         end if
     else
         ! ---------------------------------
         ! --- assign corners and edges: ---
         ! ---------------------------------
         if( nd.eq.2 )then
           ! --- TWO DIMENSIONS ----
           do side2=0,1
           do side1=0,1
             if( bc(side1,0).gt.0 .and. bc(side2,1).gt.0 )then
               ! **fix me** for exact BCs
               if( ( bc(side1,0).ne.dirichlet .and. bc(side1,0).ne.neumann ) .or. ( bc(side2,1).ne.dirichlet .and. bc(side2,1).ne.neumann ) )then
                 write(*,*) "Un-supported corner bcs =",bc(side1,0),bc(side2,1)
                 stop 2222
               end if
               symSign = +1. ! even symmetry for D-D or N-N corners
               if( (bc(side1,0).eq.dirichlet .and. bc(side2,1).eq.neumann   ) .or. (bc(side1,0).eq.neumann   .and. bc(side2,1).eq.dirichlet ) ) then
                 symSign=-1.;
               end if
               is1 = 1-2*side1
               is2 = 1-2*side2
               i1 = gridIndexRange(side1,0)
               i2 = gridIndexRange(side2,1)
               i3 = gridIndexRange(    0,2)
               do m2=0,numGhost
               do m1=0,numGhost
                 j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                 k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point 
                 u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0)
                   if( assignTwilightZone.eq.1 )then
                         call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                         call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                     u(j1,j2,j3,0) = u(j1,j2,j3,0) + ue1 - symSign*ue2 
                   else
                     ! finish me 
                   end if
               end do
               end do
             end if 
           end do
           end do
         else
           ! --- THREE DIMENSIONS ----
           write(*,*) "symmetry corners -- finish me in 3D"
           stop 9999
           do side3=0,1
           do side2=0,1
           do side1=0,1
           end do
           end do
           end do
         end if
     end if
   else
     ! This is broken, at least in 3D: May 2, 2023 
     ! 
     ! assignCornerGhostsMacro()
   end if
   return
   end
