! This file automatically generated from cornersWave.bf90 with bpp.
 subroutine cornersWave3dOrder2( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
 ! ===================================================================================
 !  Corner and edge boundary conditions for CgWave
 !
 !  gridType : 0=rectangular, 1=curvilinear
 !
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
   integer numGhost,numberOfGhostPoints,extraForNeumann,extraForDirichlet,numberOfFrequencies
   integer side1,side2,side3
   integer n1a,n1b,n2a,n2b,n3a,n3b
   integer m1a,m1b,m2a,m2b,m3a,m3b
   integer nn1a,nn1b,nn2a,nn2b,nn3a,nn3b
   integer j1a,j1b,j2a,j2b,ja,j3b
   integer extra1a,extra1b,extra2a,extra2b,extra3a,extra3b,extram
   integer cornerBC(0:2,0:2,0:2), iparc(0:10), orderOfExtrapolationForCorners
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
   logical firstTimeForCBC6
   real symSign
   integer sidea, iab(0:1,0:2)
   !     --- start statement function ----
   real bcf,mixedRHS,mixedCoeff,mixedNormalCoeff
   integer kd,m,n,component
   real uxOneSided,lap2d2Pow2
   integer bc1,bc2,bc3,edgeDirection,sideb,extra
   ! real uxxx,uxxy,uxxz,uxyy,uxzz, uyyy,uyyz,uyzz, uzzz 
   real rx,ry,rz,sx,sy,sz,tx,ty,tz
   ! define variables for getDerivatives macros
   ! #Include "../maple/declareGetDerivativesMacrosVariables.h"
   ! declareDifferenceOrder2(u,RX)
   ! #If "2" eq "4"
   !  declareDifferenceOrder4(u,RX)
   ! #End
   ! declareDifferenceOrder2(v,none)
   ! !  The next macro call will define the difference approximation statement functions
   ! defineDifferenceOrder2Components1(u,RX)
   ! #If "2" eq "4"
   !  defineDifferenceOrder4Components1(u,RX)
   ! #End
   ! defineDifferenceOrder2Components1(v,none)
 ! declare variables for getDerivatives macros
 !! turned off May 4, 2023
 !! #Include "../include/declareGetSixthDerivativesMacrosVariables.h"
 ! instead: 
 !! #Include "../include/declareGetFourthDerivativesMacrosVariables.h"
  real ux,uy,uz
  real uxxx,uxxy,uxyy,uyyy,uxxz,uxzz,uzzz,uyyz,uyzz,uxyz
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
   assignBCForImplicit             = ipar(16)
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
   c2 = c**2
   c4 = c**4
   c6 = c**6
   c8 = c**8
   twoPi = atan2(1.,1.)*8.; ! atan2(1,1)=pi/4
   assignTwilightZone=twilightZone
   if( gridType.eq.rectangular )then
     ! some macros want dr=dx for rectangular grids
     do axis=0,2
       dr(axis)=dx(axis)
     end do
   end if
   ! numberOfGhostPoints=orderOfAccuracy/2
   numberOfGhostPoints=numGhost ! now passed in 
   ! write(*,'(" cornersWave3dOrder2: dim=3, order=2")')
   ! *wdh* Nov 22, 2023 try turning of explicit BC's for implicit time-stepping
   if( .false. .and. gridIsImplicit.ne.0 .and. bcApproach==useCompatibilityBoundaryConditions .and. assignBCForImplicit==0 )then
     write(*,'(" corners: Skip explicit CBCs for implicit grid=",i4)') grid
     return
   end if
   if( t.le.3*dt .and. debug.gt.1 )then
   ! if( .true. )then
     write(*,'(" cornerWave: grid=",i4," gridType=",i2," orderOfAccuracy=",i2," uc=",i3," twilightZone=",i2)') grid,gridType,orderOfAccuracy,uc,twilightZone
     write(*,'("  addForcingBC=",i4," forcingOption=",i4," assignKnownSolutionAtBoundaries=",i4)') addForcingBC, forcingOption, assignKnownSolutionAtBoundaries
     write(*,'("  t=",e10.2," dt=",e10.2," knownSolutionOption=",i4," REAL_MIN=",e10.2)') t,dt,knownSolutionOption,REAL_MIN
     write(*,'("  abcWave: c=",e14.6," cEM2=",e14.6)') c,cEM2
     write(*,'("  useUpwindDissipation=",i2," numGhost=",i2)') useUpwindDissipation,numGhost
     write(*,'("  assignBCForImplicit=",i4," bcApproach=",i4," gridIsImplicit=",i2)') assignBCForImplicit,bcApproach,gridIsImplicit
     write(*,'("  boundaryCondition=",6i4)') ((boundaryCondition(side,axis),side=0,1),axis=0,2)
   end if
   ! if( bcApproach.eq.useCompatibilityBoundaryConditions )then
   !   write(*,'("corners: ERROR: useCompatibilityBoundaryConditions not implemented yet.")') 
   !   stop 1010
   ! end if
   if( bcApproach.eq.useLocalCompatibilityBoundaryConditions )then
     write(*,'("corners: ERROR: useLocalCompatibilityBoundaryConditions not implemented yet.")') 
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
   if( assignKnownSolutionAtBoundaries.eq.1 )then
     if( knownSolutionOption.eq.planeWave )then
       ! get parameter values from the C++ data-base
        ok=getReal(pdb,'ampPlaneWave',ampPlaneWave) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find ampPlaneWave")') 
          stop 1133
        end if
        ok=getReal(pdb,'kxPlaneWave',kxPlaneWave) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find kxPlaneWave")') 
          stop 1133
        end if
        ok=getReal(pdb,'kyPlaneWave',kyPlaneWave) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find kyPlaneWave")') 
          stop 1133
        end if
        ok=getReal(pdb,'kzPlaneWave',kzPlaneWave) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find kzPlaneWave")') 
          stop 1133
        end if
        ok=getReal(pdb,'omegaPlaneWave',omegaPlaneWave) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find omegaPlaneWave")') 
          stop 1133
        end if
        ok=getInt(pdb,'solveForScatteredField',solveForScatteredField) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getInt:ERROR: unable to find solveForScatteredField")') 
          stop 1122
        end if
       if( solveForScatteredField==1 )then
          ! If we solve for the scattered field then we flip the sign of the plane wave since this has been subtracted out
          ampPlaneWave = -ampPlaneWave
           ok=getInt(pdb,'solveHelmholtz',solveHelmholtz) 
           if( ok.eq.0 )then
             write(*,'("*** corners:getInt:ERROR: unable to find solveHelmholtz")') 
             stop 1122
           end if
          if( solveHelmholtz==1 )then 
             ! Get adjusted omega:
              ok=getReal(pdb,'omega',omega) 
              if( ok.eq.0 )then
                write(*,'("*** corners:getReal:ERROR: unable to find omega")') 
                stop 1133
              end if
             if( t.le.2*dt )then
               write(*,'(" corners:solveHelmholtz:scattering: Use adjusted omega=",1pe15.8," in place of omegaPlaneWave=",1pe15.8)') omega,omegaPlaneWave
             end if
             omegaPlaneWave = omega
          end if
          ! getIntParameter(solveHelmholtz)  
          ! if( solveHelmholtz==1 )then
          !    ! We are solving a Helmholtz problem : check that omega in the plane wave solution matches the omega for Helmholtz
          !    getRealParameter(omega)
          !    omegaTol = 1e-10  ! **FIX ME** use a multiple for REAL_EPSILON 
          !    if( abs(omega-omegaPlaneWave) .gt. omegaTol * omega )then
          !      write(*,'(" corners: solveForScatteredField=1, boundary forcing is a plane wave")') 
          !      write(*,'(" corners: ERROR: omegaPlaneWave=",e18.8," is not equal to omega(Helmholtz)=",e18.8)') omegaPlaneWave,omega
          !      stop 1234
          !    end if
          ! end if
       end if
       if(  t.le.2*dt .and. debug.gt.1  )then
         write(*,'(" corners:  knownSolutionOption=planeWave: solveForScatteredField=",i2," ampPlaneWave=",e10.2," kxPlaneWave=",e10.2," kyPlaneWave=",e10.2," omegaPlaneWave=",e14.4)') solveForScatteredField, ampPlaneWave,kxPlaneWave,kyPlaneWave,omegaPlaneWave
       end if 
       ! write(*,'(" corners:  knownSolutionOption=planeWave: solveForScatteredField=",i2," omegaPlaneWave=",e14.4)') solveForScatteredField, omegaPlaneWave
     else if( knownSolutionOption.eq.gaussianPlaneWave )then
       ! Get the parameters in the Gaussian plane wave (Set in userDefinedKnownSolution)
        ok=getReal(pdb,'kxGPW',kxGPW) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find kxGPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'kyGPW',kyGPW) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find kyGPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'kzGPW',kzGPW) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find kzGPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'x0GPW',x0GPW) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find x0GPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'y0GPW',y0GPW) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find y0GPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'z0GPW',z0GPW) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find z0GPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'k0GPW',k0GPW) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find k0GPW")') 
          stop 1133
        end if
        ok=getReal(pdb,'betaGPW',betaGPW) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find betaGPW")') 
          stop 1133
        end if
       if(  t.le.dt .and. debug.ge.0  )then
         write(*,'(" corners:  knownSolutionOption=gaussianPlaneWave: kx,ky,kz=",3(1pe10.2)," x0,y0,z0=",3(1pe10.2)," k0,beta=",2(1pe10.2))') kxGPW,kyGPW,kzGPW,x0GPW,y0GPW,z0GPW,k0GPW,betaGPW
       end if           
     else if( knownSolutionOption.eq.boxHelmholtz )then
       ! get parameter values from the C++ data-base
        ok=getReal(pdb,'kxBoxHelmholtz',kxBoxHelmholtz) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find kxBoxHelmholtz")') 
          stop 1133
        end if
        ok=getReal(pdb,'kyBoxHelmholtz',kyBoxHelmholtz) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find kyBoxHelmholtz")') 
          stop 1133
        end if
        ok=getReal(pdb,'kzBoxHelmholtz',kzBoxHelmholtz) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find kzBoxHelmholtz")') 
          stop 1133
        end if
        ok=getReal(pdb,'omegaBoxHelmholtz',omegaBoxHelmholtz) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find omegaBoxHelmholtz")') 
          stop 1133
        end if
       coswt = cos(omegaBoxHelmholtz*t)
       assignKnownSolutionAtBoundaries=1  ! for inhomogeneous BCs
       if(  t.le.dt .and. debug.ge.1   )then
         write(*,'(" corners:  assignKnownSolutionAtBoundaries=",i4)') assignKnownSolutionAtBoundaries
         write(*,'(" corners:  numberOfFrequencies=",i4)') numberOfFrequencies
         write(*,'(" corners:  frequencyArray=",10(1pe12.4,1x))') (frequencyArray(freq),freq=0,numberOfFrequencies-1)
         write(*,'(" corners:  knownSolutionOption=boxHelmholtz: kx,ky,kz,omega=",4e10.2)') kxBoxHelmholtz,kyBoxHelmholtz,kzBoxHelmholtz,omegaBoxHelmholtz
       end if
     else if( knownSolutionOption.eq.polyPeriodic )then
       ! get parameter values from the C++ data-base
        ok=getReal(pdb,'omegaPolyPeriodic',omegaPolyPeriodic) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find omegaPolyPeriodic")') 
          stop 1133
        end if
        ok=getReal(pdb,'a0PolyPeriodic',a0PolyPeriodic) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find a0PolyPeriodic")') 
          stop 1133
        end if
        ok=getReal(pdb,'a1PolyPeriodic',a1PolyPeriodic) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find a1PolyPeriodic")') 
          stop 1133
        end if
        ok=getReal(pdb,'b1PolyPeriodic',b1PolyPeriodic) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find b1PolyPeriodic")') 
          stop 1133
        end if
        ok=getReal(pdb,'c1PolyPeriodic',c1PolyPeriodic) 
        if( ok.eq.0 )then
          write(*,'("*** corners:getReal:ERROR: unable to find c1PolyPeriodic")') 
          stop 1133
        end if
       coswt = cos(omegaPolyPeriodic*t)
       if(  t.le.dt .and. debug.gt.1   )then
         write(*,'(" corners:  knownSolutionOption=polyPeriodic: a0,a1,b1,c1,omega=",5e10.2)') a0PolyPeriodic,a1PolyPeriodic,b1PolyPeriodic,c1PolyPeriodic,omegaPolyPeriodic
       end if
     else if( knownSolutionOption.ne.0 .and. knownSolutionOption.ne.otherKnownSolution )then
       write(*,'("corners:ERROR: unknown knownSolutionOption=",i6)') knownSolutionOption
       stop 1111
     end if 
   end if
   ! TEST: 
   ! getRealParameter(omega)
   ! getRealParameter(cfl)
   ! write(*,'(" corners:  cfl=",e10.2)') cfl
   if( uc.lt.0 .or. uc.ge.numberOfComponents )then
     write(*,'("corners:ERROR: invalid uc=",i6," but numberOfComponents=",i3)')  uc,numberOfComponents
     stop 1111
   end if
   epsx=REAL_MIN*100.  ! for normal
   if( orderOfAccuracy.ne.2 .and. orderOfAccuracy.ne.4 .and. orderOfAccuracy.ne.6 .and. orderOfAccuracy.ne.8 )then
     write(*,'("corners:ERROR: orderOfAccuracy is not 2, 4 or 6, orderOfAccuracy=",i4)') orderOfAccuracy
     stop 1111
   end if
   if( assignBCForImplicit.eq.1 .or. assignBCForImplicit.eq.2 )then
     ! --- FILL IN THE RHS FOR THE IMPLICIT MATRIX ---
     extra=0
     if( forcingOption.eq.noForcing )then
         ! ---------------------------------
         ! --- assign corners and edges: ---
         ! ---------------------------------
         if( .false. .and. t.le.2.*dt )then
           write(*,'("assign GENERAL corner ghost points forcing=noForcing, method=implicit, nd=",i3, " ***FINISH ME***")') nd
         end if
         if( nd.eq.2 )then
           ! ------ TWO DIMENSIONS ----
           if( debug.gt.1 .and. t.le.2*dt )then
             write(*,'("assign corner ghost points in 2D - general and curvilinear case")')
           end if
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
               if( ( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann ) )then
                 write(*,*) "Un-supported corner bc1=",bc1
                 stop 2222
               end if
               if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                 write(*,*) "Un-supported corner bc2 =",bc2
                 stop 2222
               end if
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign = +1. 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if   
               ! if( (bc(side1,0).eq.dirichlet .and. bc(side2,1).eq.neumann   ) .or. !     (bc(side1,0).eq.neumann   .and. bc(side2,1).eq.dirichlet ) ) then
               !   symSign=-1.;
               ! end if
               is1 = 1-2*side1
               is2 = 1-2*side2
               i1 = gridIndexRange(side1,0)
               i2 = gridIndexRange(side2,1)
               i3 = gridIndexRange(    0,2)
               if( mask(i1,i2,i3).gt.0 )then
                 do m2=1,numberOfGhostPoints
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point 
                       u(j1,j2,j3,0) = 0.
                     ! write(*,'("Corner ghost: set u(",3i3,")=u(",3i3,")*",f5.2)') j1,j2,j3,k1,k2,k3,symSign
                 end do
                 end do
               end if
             end if 
           end do
           end do
         else
           ! --- THREE DIMENSIONS ----
           ! ------ TWO DIMENSIONS ----
           if( debug.gt.1 .and. t.le.2*dt )then
             write(*,'("Assign corner and edge ghost points in 3D - general and curvilinear case")')
           end if
           ! write(*,'("assign general edge ghost points in 3D forcing=noForcing")')
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
                 n3a=gridIndexRange(0,2)-extra
                 n3b=gridIndexRange(1,2)+extra
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side2,1)
               else if( edgeDirection.eq.1 )then
                 is2=0
                 n1a=gridIndexRange(side1,0)
                 n1b=gridIndexRange(side1,0)
                 n2a=gridIndexRange(    0,1)-extra
                 n2b=gridIndexRange(    1,1)+extra
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side3,2)
               else 
                 is1=0  
                 n1a=gridIndexRange(    0,0)-extra
                 n1b=gridIndexRange(    1,0)+extra
                 n2a=gridIndexRange(side2,1)
                 n2b=gridIndexRange(side2,1)
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side2,1)
                 bc2=boundaryCondition(side3,2)
               end if
               if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
                 ! -- this is an edge between two physical boundaries --
                 if( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann )then
                   write(*,*) "Un-supported edge bc1=",bc1
                   stop 2222
                 end if
                 if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                   write(*,*) "Un-supported edge bc2 =",bc2
                   stop 2222
                 end if
                 ! symSign = 1 : even symmetry
                 !          -1 : odd symmetry
                 symSign=1; 
                 if( bc1==dirichlet .or. bc1==exactBC )then
                   symSign = -symSign
                 end if
                 if( bc2==dirichlet .or. bc2==exactBC )then
                   symSign = -symSign
                 end if 
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                   if( mask(i1,i2,i3).gt.0 )then
                     do m1=1,numberOfGhostPoints
                     do m2=1,numberOfGhostPoints 
                       if( edgeDirection==0 )then
                         ! edge lies along i1=const
                         j1 = i1; j2=i2-is2*m1; j3=i3-is3*m2; ! ghost 
                         k1 = i1; k2=i2+is2*m1; k3=i3+is3*m2; ! interior point              
                       else if( edgeDirection==1 )then
                         ! edge lies along i2=const
                         j1 = i1-is1*m1; j2=i2; j3=i3-is3*m2; ! ghost 
                         k1 = i1+is1*m1; k2=i2; k3=i3+is3*m2; ! interior point  
                       else
                        ! edge lies along i3=constant
                         j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                         k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point  
                       end if
                           u(j1,j2,j3,0) = 0.
                         ! write(*,'("Edge=",i2," sidea,sideb=",2i2," set u(",3i3,")=u(",3i3,")*",f5.2)') edgeDirection,sidea,sideb,j1,j2,j3,k1,k2,k3,symSign
                     end do
                     end do     
                   end if ! if mask 
                 end do
                 end do
                 end do
               end if ! bc1>0 and bc2>0
             end do ! end do sideb
             end do ! end do sidea
           end do ! end edgeDirection
           do side3=0,1
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             bc3 = bc(side3,2)
             if( bc1.gt.0 .and.  bc2.gt.0 .and. bc3.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC .and. bc3.ne.exactBC )then
               ! --- Vertex where three physical boundaries meet ---
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign=1; 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if        
               if( bc3==dirichlet .or. bc3==exactBC )then
                 symSign = -symSign
               end if  
               ! vertex is (i1,i2,i3)
               i1=gridIndexRange(side1,0)
               i2=gridIndexRange(side2,1)
               i3=gridIndexRange(side3,2)
               if( mask(i1,i2,i3).gt.0 )then
                 is1 = 1-2*side1
                 is2 = 1-2*side2
                 is3 = 1-2*side3          
                 do m3=1,numberOfGhostPoints  
                 do m2=1,numberOfGhostPoints  
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3-is3*m3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3+is3*m3; ! interior point             
                       u(j1,j2,j3,0) = 0.
                     ! write(*,'("Vertex: set u(",3i3,")=u(",3i3,")*",f5.2)') j1,j2,j3,k1,k2,k3,symSign
                 end do
                 end do          
                 end do
               end if
             end if          
           end do
           end do
           end do
         end if
     else
         ! ---------------------------------
         ! --- assign corners and edges: ---
         ! ---------------------------------
         if( .false. .and. t.le.2.*dt )then
           write(*,'("assign GENERAL corner ghost points forcing=forcing, method=implicit, nd=",i3, " ***FINISH ME***")') nd
         end if
         if( nd.eq.2 )then
           ! ------ TWO DIMENSIONS ----
           if( debug.gt.1 .and. t.le.2*dt )then
             write(*,'("assign corner ghost points in 2D - general and curvilinear case")')
           end if
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
               if( ( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann ) )then
                 write(*,*) "Un-supported corner bc1=",bc1
                 stop 2222
               end if
               if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                 write(*,*) "Un-supported corner bc2 =",bc2
                 stop 2222
               end if
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign = +1. 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if   
               ! if( (bc(side1,0).eq.dirichlet .and. bc(side2,1).eq.neumann   ) .or. !     (bc(side1,0).eq.neumann   .and. bc(side2,1).eq.dirichlet ) ) then
               !   symSign=-1.;
               ! end if
               is1 = 1-2*side1
               is2 = 1-2*side2
               i1 = gridIndexRange(side1,0)
               i2 = gridIndexRange(side2,1)
               i3 = gridIndexRange(    0,2)
               if( mask(i1,i2,i3).gt.0 )then
                 do m2=1,numberOfGhostPoints
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point 
                     if( assignTwilightZone.eq.1 )then
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                         u(j1,j2,j3,0) =                       -symSign*ue2 + ue1 
                     else
                       ! finish me 
                     end if
                 end do
                 end do
               end if
             end if 
           end do
           end do
         else
           ! --- THREE DIMENSIONS ----
           ! ------ TWO DIMENSIONS ----
           if( debug.gt.1 .and. t.le.2*dt )then
             write(*,'("Assign corner and edge ghost points in 3D - general and curvilinear case")')
           end if
           ! write(*,'("assign general edge ghost points in 3D forcing=forcing")')
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
                 n3a=gridIndexRange(0,2)-extra
                 n3b=gridIndexRange(1,2)+extra
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side2,1)
               else if( edgeDirection.eq.1 )then
                 is2=0
                 n1a=gridIndexRange(side1,0)
                 n1b=gridIndexRange(side1,0)
                 n2a=gridIndexRange(    0,1)-extra
                 n2b=gridIndexRange(    1,1)+extra
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side3,2)
               else 
                 is1=0  
                 n1a=gridIndexRange(    0,0)-extra
                 n1b=gridIndexRange(    1,0)+extra
                 n2a=gridIndexRange(side2,1)
                 n2b=gridIndexRange(side2,1)
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side2,1)
                 bc2=boundaryCondition(side3,2)
               end if
               if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
                 ! -- this is an edge between two physical boundaries --
                 if( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann )then
                   write(*,*) "Un-supported edge bc1=",bc1
                   stop 2222
                 end if
                 if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                   write(*,*) "Un-supported edge bc2 =",bc2
                   stop 2222
                 end if
                 ! symSign = 1 : even symmetry
                 !          -1 : odd symmetry
                 symSign=1; 
                 if( bc1==dirichlet .or. bc1==exactBC )then
                   symSign = -symSign
                 end if
                 if( bc2==dirichlet .or. bc2==exactBC )then
                   symSign = -symSign
                 end if 
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                   if( mask(i1,i2,i3).gt.0 )then
                     do m1=1,numberOfGhostPoints
                     do m2=1,numberOfGhostPoints 
                       if( edgeDirection==0 )then
                         ! edge lies along i1=const
                         j1 = i1; j2=i2-is2*m1; j3=i3-is3*m2; ! ghost 
                         k1 = i1; k2=i2+is2*m1; k3=i3+is3*m2; ! interior point              
                       else if( edgeDirection==1 )then
                         ! edge lies along i2=const
                         j1 = i1-is1*m1; j2=i2; j3=i3-is3*m2; ! ghost 
                         k1 = i1+is1*m1; k2=i2; k3=i3+is3*m2; ! interior point  
                       else
                        ! edge lies along i3=constant
                         j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                         k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point  
                       end if
                         if( assignTwilightZone.eq.1 )then
                               call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),xy(j1,j2,j3,2),t,uc,ue1 )
                               call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),xy(k1,k2,k3,2),t,uc,ue2 )
                             u(j1,j2,j3,0) =                       - symSign*ue2 + ue1 
               ! write(*,'("Edge: ghost j1,j2,j3=",3i4," k1,k2,k3=",3i4," u=",e10.2," err=",e9.2)') j1,j2,j3,k1,k2,k3,u(j1,j2,j3,0),u(j1,j2,j3,0)-ue1
                         else
                           ! finish me 
                           write(*,*) "EdgeBC: finish me for forcing option"
                           stop 3333
                         end if
                     end do
                     end do     
                   end if ! if mask 
                 end do
                 end do
                 end do
               end if ! bc1>0 and bc2>0
             end do ! end do sideb
             end do ! end do sidea
           end do ! end edgeDirection
           do side3=0,1
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             bc3 = bc(side3,2)
             if( bc1.gt.0 .and.  bc2.gt.0 .and. bc3.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC .and. bc3.ne.exactBC )then
               ! --- Vertex where three physical boundaries meet ---
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign=1; 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if        
               if( bc3==dirichlet .or. bc3==exactBC )then
                 symSign = -symSign
               end if  
               ! vertex is (i1,i2,i3)
               i1=gridIndexRange(side1,0)
               i2=gridIndexRange(side2,1)
               i3=gridIndexRange(side3,2)
               if( mask(i1,i2,i3).gt.0 )then
                 is1 = 1-2*side1
                 is2 = 1-2*side2
                 is3 = 1-2*side3          
                 do m3=1,numberOfGhostPoints  
                 do m2=1,numberOfGhostPoints  
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3-is3*m3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3+is3*m3; ! interior point             
                     if( assignTwilightZone.eq.1 )then
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),xy(j1,j2,j3,2),t,uc,ue1 )
                           call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),xy(k1,k2,k3,2),t,uc,ue2 )
                         u(j1,j2,j3,0) =                       - symSign*ue2 + ue1 
                     else
                       ! finish me 
                       write(*,*) "SymBC: finish me for forcing option"
                       stop 3333
                     end if
                 end do
                 end do          
                 end do
               end if
             end if          
           end do
           end do
           end do
         end if
     end if    
     ! ---------------- RETURN ---------------
     return
   end if
   !  --- Assign ghost points outside corners ---
   ! if( .true. .and. orderOfAccuracy==4 .and. gridType==rectangular .and. bcApproach==useCompatibilityBoundaryConditions )then 
   if( gridType==rectangular  )then 
     if( t.le.3*dt .and. debug.gt.1 )then
       write(*,'("cornerWave: Assign corner ghost (symmetry conditions)")')
     end if
     extra=0
     if( forcingOption.eq.noForcing )then
         ! ---------------------------------
         ! --- assign corners and edges: ---
         ! ---------------------------------
         ! write(*,'("assign symmetry corner ghost points forcing=noForcing, nd=",i3)') nd
         if( nd.eq.2 )then
           ! ------ TWO DIMENSIONS ----
           if( debug.gt.0 .and. t.le.2*dt )then
             write(*,'("assign symmetry Ghost points in 2D")')
           end if
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
               if( ( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann ) )then
                 write(*,*) "Un-supported corner bc1=",bc1
                 stop 2222
               end if
               if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                 write(*,*) "Un-supported corner bc2 =",bc2
                 stop 2222
               end if
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign = +1. 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if   
               ! if( (bc(side1,0).eq.dirichlet .and. bc(side2,1).eq.neumann   ) .or. !     (bc(side1,0).eq.neumann   .and. bc(side2,1).eq.dirichlet ) ) then
               !   symSign=-1.;
               ! end if
               is1 = 1-2*side1
               is2 = 1-2*side2
               i1 = gridIndexRange(side1,0)
               i2 = gridIndexRange(side2,1)
               i3 = gridIndexRange(    0,2)
               if( mask(i1,i2,i3).gt.0 )then
                 do m2=1,numberOfGhostPoints
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point 
                     u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0)
                     ! write(*,'("Corner ghost: set u(",3i3,")=u(",3i3,")*",f5.2)') j1,j2,j3,k1,k2,k3,symSign
                 end do
                 end do
               end if
             end if 
           end do
           end do
         else
           ! --- THREE DIMENSIONS ----
           ! write(*,'("assign symmetry edge ghost points in 3D forcing=noForcing")')
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
                 n3a=gridIndexRange(0,2)-extra
                 n3b=gridIndexRange(1,2)+extra
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side2,1)
               else if( edgeDirection.eq.1 )then
                 is2=0
                 n1a=gridIndexRange(side1,0)
                 n1b=gridIndexRange(side1,0)
                 n2a=gridIndexRange(    0,1)-extra
                 n2b=gridIndexRange(    1,1)+extra
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side3,2)
               else 
                 is1=0  
                 n1a=gridIndexRange(    0,0)-extra
                 n1b=gridIndexRange(    1,0)+extra
                 n2a=gridIndexRange(side2,1)
                 n2b=gridIndexRange(side2,1)
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side2,1)
                 bc2=boundaryCondition(side3,2)
               end if
               if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
                 ! -- this is an edge between two physical boundaries --
                 if( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann )then
                   write(*,*) "Un-supported edge bc1=",bc1
                   stop 2222
                 end if
                 if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                   write(*,*) "Un-supported edge bc2 =",bc2
                   stop 2222
                 end if
                 ! symSign = 1 : even symmetry
                 !          -1 : odd symmetry
                 symSign=1; 
                 if( bc1==dirichlet .or. bc1==exactBC )then
                   symSign = -symSign
                 end if
                 if( bc2==dirichlet .or. bc2==exactBC )then
                   symSign = -symSign
                 end if 
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                   if( mask(i1,i2,i3).gt.0 )then
                     do m1=1,numberOfGhostPoints
                     do m2=1,numberOfGhostPoints 
                       if( edgeDirection==0 )then
                         ! edge lies along i1=const
                         j1 = i1; j2=i2-is2*m1; j3=i3-is3*m2; ! ghost 
                         k1 = i1; k2=i2+is2*m1; k3=i3+is3*m2; ! interior point              
                       else if( edgeDirection==1 )then
                         ! edge lies along i2=const
                         j1 = i1-is1*m1; j2=i2; j3=i3-is3*m2; ! ghost 
                         k1 = i1+is1*m1; k2=i2; k3=i3+is3*m2; ! interior point  
                       else
                        ! edge lies along i3=constant
                         j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                         k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point  
                       end if
                         u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0)
                         ! write(*,'("Edge=",i2," sidea,sideb=",2i2," set u(",3i3,")=u(",3i3,")*",f5.2)') edgeDirection,sidea,sideb,j1,j2,j3,k1,k2,k3,symSign
                     end do
                     end do     
                   end if ! if mask 
                 end do
                 end do
                 end do
               end if ! bc1>0 and bc2>0
             end do ! end do sideb
             end do ! end do sidea
           end do ! end edgeDirection
           ! write(*,*) "symmetry corners -- finish me in 3D"
           ! stop 9999
           do side3=0,1
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             bc3 = bc(side3,2)
             if( bc1.gt.0 .and.  bc2.gt.0 .and. bc3.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC .and. bc3.ne.exactBC )then
               ! --- Vertex where three physical boundaries meet ---
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign=1; 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if        
               if( bc3==dirichlet .or. bc3==exactBC )then
                 symSign = -symSign
               end if  
               ! vertex is (i1,i2,i3)
               i1=gridIndexRange(side1,0)
               i2=gridIndexRange(side2,1)
               i3=gridIndexRange(side3,2)
               if( mask(i1,i2,i3).gt.0 )then
                 is1 = 1-2*side1
                 is2 = 1-2*side2
                 is3 = 1-2*side3          
                 do m3=1,numberOfGhostPoints  
                 do m2=1,numberOfGhostPoints  
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3-is3*m3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3+is3*m3; ! interior point             
                     u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0)
                     ! write(*,'("Vertex: set u(",3i3,")=u(",3i3,")*",f5.2)') j1,j2,j3,k1,k2,k3,symSign
                 end do
                 end do          
                 end do
               end if
             end if          
           end do
           end do
           end do
         end if
     else
         ! ---------------------------------
         ! --- assign corners and edges: ---
         ! ---------------------------------
         ! write(*,'("assign symmetry corner ghost points forcing=forcing, nd=",i3)') nd
         if( nd.eq.2 )then
           ! ------ TWO DIMENSIONS ----
           if( debug.gt.0 .and. t.le.2*dt )then
             write(*,'("assign symmetry Ghost points in 2D")')
           end if
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
               if( ( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann ) )then
                 write(*,*) "Un-supported corner bc1=",bc1
                 stop 2222
               end if
               if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                 write(*,*) "Un-supported corner bc2 =",bc2
                 stop 2222
               end if
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign = +1. 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if   
               ! if( (bc(side1,0).eq.dirichlet .and. bc(side2,1).eq.neumann   ) .or. !     (bc(side1,0).eq.neumann   .and. bc(side2,1).eq.dirichlet ) ) then
               !   symSign=-1.;
               ! end if
               is1 = 1-2*side1
               is2 = 1-2*side2
               i1 = gridIndexRange(side1,0)
               i2 = gridIndexRange(side2,1)
               i3 = gridIndexRange(    0,2)
               if( mask(i1,i2,i3).gt.0 )then
                 do m2=1,numberOfGhostPoints
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point 
                     if( assignTwilightZone.eq.1 )then
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                       u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0) -symSign*ue2 + ue1 
                     else
                       ! finish me 
                     end if
                 end do
                 end do
               end if
             end if 
           end do
           end do
         else
           ! --- THREE DIMENSIONS ----
           ! write(*,'("assign symmetry edge ghost points in 3D forcing=forcing")')
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
                 n3a=gridIndexRange(0,2)-extra
                 n3b=gridIndexRange(1,2)+extra
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side2,1)
               else if( edgeDirection.eq.1 )then
                 is2=0
                 n1a=gridIndexRange(side1,0)
                 n1b=gridIndexRange(side1,0)
                 n2a=gridIndexRange(    0,1)-extra
                 n2b=gridIndexRange(    1,1)+extra
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side3,2)
               else 
                 is1=0  
                 n1a=gridIndexRange(    0,0)-extra
                 n1b=gridIndexRange(    1,0)+extra
                 n2a=gridIndexRange(side2,1)
                 n2b=gridIndexRange(side2,1)
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side2,1)
                 bc2=boundaryCondition(side3,2)
               end if
               if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
                 ! -- this is an edge between two physical boundaries --
                 if( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann )then
                   write(*,*) "Un-supported edge bc1=",bc1
                   stop 2222
                 end if
                 if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                   write(*,*) "Un-supported edge bc2 =",bc2
                   stop 2222
                 end if
                 ! symSign = 1 : even symmetry
                 !          -1 : odd symmetry
                 symSign=1; 
                 if( bc1==dirichlet .or. bc1==exactBC )then
                   symSign = -symSign
                 end if
                 if( bc2==dirichlet .or. bc2==exactBC )then
                   symSign = -symSign
                 end if 
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                   if( mask(i1,i2,i3).gt.0 )then
                     do m1=1,numberOfGhostPoints
                     do m2=1,numberOfGhostPoints 
                       if( edgeDirection==0 )then
                         ! edge lies along i1=const
                         j1 = i1; j2=i2-is2*m1; j3=i3-is3*m2; ! ghost 
                         k1 = i1; k2=i2+is2*m1; k3=i3+is3*m2; ! interior point              
                       else if( edgeDirection==1 )then
                         ! edge lies along i2=const
                         j1 = i1-is1*m1; j2=i2; j3=i3-is3*m2; ! ghost 
                         k1 = i1+is1*m1; k2=i2; k3=i3+is3*m2; ! interior point  
                       else
                        ! edge lies along i3=constant
                         j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                         k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point  
                       end if
                         if( assignTwilightZone.eq.1 )then
                               call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),xy(j1,j2,j3,2),t,uc,ue1 )
                               call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),xy(k1,k2,k3,2),t,uc,ue2 )
                           u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0) - symSign*ue2 + ue1 
               ! write(*,'("Edge: ghost j1,j2,j3=",3i4," k1,k2,k3=",3i4," u=",e10.2," err=",e9.2)') j1,j2,j3,k1,k2,k3,u(j1,j2,j3,0),u(j1,j2,j3,0)-ue1
                         else
                           ! finish me 
                           write(*,*) "EdgeBC: finish me for forcing option"
                           stop 3333
                         end if
                     end do
                     end do     
                   end if ! if mask 
                 end do
                 end do
                 end do
               end if ! bc1>0 and bc2>0
             end do ! end do sideb
             end do ! end do sidea
           end do ! end edgeDirection
           ! write(*,*) "symmetry corners -- finish me in 3D"
           ! stop 9999
           do side3=0,1
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             bc3 = bc(side3,2)
             if( bc1.gt.0 .and.  bc2.gt.0 .and. bc3.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC .and. bc3.ne.exactBC )then
               ! --- Vertex where three physical boundaries meet ---
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign=1; 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if        
               if( bc3==dirichlet .or. bc3==exactBC )then
                 symSign = -symSign
               end if  
               ! vertex is (i1,i2,i3)
               i1=gridIndexRange(side1,0)
               i2=gridIndexRange(side2,1)
               i3=gridIndexRange(side3,2)
               if( mask(i1,i2,i3).gt.0 )then
                 is1 = 1-2*side1
                 is2 = 1-2*side2
                 is3 = 1-2*side3          
                 do m3=1,numberOfGhostPoints  
                 do m2=1,numberOfGhostPoints  
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3-is3*m3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3+is3*m3; ! interior point             
                     if( assignTwilightZone.eq.1 )then
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),xy(j1,j2,j3,2),t,uc,ue1 )
                           call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),xy(k1,k2,k3,2),t,uc,ue2 )
                       u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0) - symSign*ue2 + ue1 
                     else
                       ! finish me 
                       write(*,*) "SymBC: finish me for forcing option"
                       stop 3333
                     end if
                 end do
                 end do          
                 end do
               end if
             end if          
           end do
           end do
           end do
         end if
     end if
   else if( .true. )then
     if( t.le.3*dt .and. debug.gt.1 )then
       write(*,'("cornerWave: Assign corner ghost (general case, curvilinear grids)")')
     end if
     extra=0
     if( forcingOption.eq.noForcing )then
         ! ---------------------------------
         ! --- assign corners and edges: ---
         ! ---------------------------------
         if( .false. .and. t.le.2.*dt )then
           write(*,'("assign GENERAL corner ghost points forcing=noForcing, method=explicit, nd=",i3, " ***FINISH ME***")') nd
         end if
         if( nd.eq.2 )then
           ! ------ TWO DIMENSIONS ----
           if( debug.gt.1 .and. t.le.2*dt )then
             write(*,'("assign corner ghost points in 2D - general and curvilinear case")')
           end if
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
               if( ( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann ) )then
                 write(*,*) "Un-supported corner bc1=",bc1
                 stop 2222
               end if
               if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                 write(*,*) "Un-supported corner bc2 =",bc2
                 stop 2222
               end if
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign = +1. 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if   
               ! if( (bc(side1,0).eq.dirichlet .and. bc(side2,1).eq.neumann   ) .or. !     (bc(side1,0).eq.neumann   .and. bc(side2,1).eq.dirichlet ) ) then
               !   symSign=-1.;
               ! end if
               is1 = 1-2*side1
               is2 = 1-2*side2
               i1 = gridIndexRange(side1,0)
               i2 = gridIndexRange(side2,1)
               i3 = gridIndexRange(    0,2)
               if( mask(i1,i2,i3).gt.0 )then
                 do m2=1,numberOfGhostPoints
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point 
                       u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0)
                     ! write(*,'("Corner ghost: set u(",3i3,")=u(",3i3,")*",f5.2)') j1,j2,j3,k1,k2,k3,symSign
                 end do
                 end do
               end if
             end if 
           end do
           end do
         else
           ! --- THREE DIMENSIONS ----
           ! ------ TWO DIMENSIONS ----
           if( debug.gt.1 .and. t.le.2*dt )then
             write(*,'("Assign corner and edge ghost points in 3D - general and curvilinear case")')
           end if
           ! write(*,'("assign general edge ghost points in 3D forcing=noForcing")')
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
                 n3a=gridIndexRange(0,2)-extra
                 n3b=gridIndexRange(1,2)+extra
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side2,1)
               else if( edgeDirection.eq.1 )then
                 is2=0
                 n1a=gridIndexRange(side1,0)
                 n1b=gridIndexRange(side1,0)
                 n2a=gridIndexRange(    0,1)-extra
                 n2b=gridIndexRange(    1,1)+extra
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side3,2)
               else 
                 is1=0  
                 n1a=gridIndexRange(    0,0)-extra
                 n1b=gridIndexRange(    1,0)+extra
                 n2a=gridIndexRange(side2,1)
                 n2b=gridIndexRange(side2,1)
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side2,1)
                 bc2=boundaryCondition(side3,2)
               end if
               if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
                 ! -- this is an edge between two physical boundaries --
                 if( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann )then
                   write(*,*) "Un-supported edge bc1=",bc1
                   stop 2222
                 end if
                 if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                   write(*,*) "Un-supported edge bc2 =",bc2
                   stop 2222
                 end if
                 ! symSign = 1 : even symmetry
                 !          -1 : odd symmetry
                 symSign=1; 
                 if( bc1==dirichlet .or. bc1==exactBC )then
                   symSign = -symSign
                 end if
                 if( bc2==dirichlet .or. bc2==exactBC )then
                   symSign = -symSign
                 end if 
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                   if( mask(i1,i2,i3).gt.0 )then
                     do m1=1,numberOfGhostPoints
                     do m2=1,numberOfGhostPoints 
                       if( edgeDirection==0 )then
                         ! edge lies along i1=const
                         j1 = i1; j2=i2-is2*m1; j3=i3-is3*m2; ! ghost 
                         k1 = i1; k2=i2+is2*m1; k3=i3+is3*m2; ! interior point              
                       else if( edgeDirection==1 )then
                         ! edge lies along i2=const
                         j1 = i1-is1*m1; j2=i2; j3=i3-is3*m2; ! ghost 
                         k1 = i1+is1*m1; k2=i2; k3=i3+is3*m2; ! interior point  
                       else
                        ! edge lies along i3=constant
                         j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                         k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point  
                       end if
                           u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0)
                         ! write(*,'("Edge=",i2," sidea,sideb=",2i2," set u(",3i3,")=u(",3i3,")*",f5.2)') edgeDirection,sidea,sideb,j1,j2,j3,k1,k2,k3,symSign
                     end do
                     end do     
                   end if ! if mask 
                 end do
                 end do
                 end do
               end if ! bc1>0 and bc2>0
             end do ! end do sideb
             end do ! end do sidea
           end do ! end edgeDirection
           do side3=0,1
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             bc3 = bc(side3,2)
             if( bc1.gt.0 .and.  bc2.gt.0 .and. bc3.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC .and. bc3.ne.exactBC )then
               ! --- Vertex where three physical boundaries meet ---
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign=1; 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if        
               if( bc3==dirichlet .or. bc3==exactBC )then
                 symSign = -symSign
               end if  
               ! vertex is (i1,i2,i3)
               i1=gridIndexRange(side1,0)
               i2=gridIndexRange(side2,1)
               i3=gridIndexRange(side3,2)
               if( mask(i1,i2,i3).gt.0 )then
                 is1 = 1-2*side1
                 is2 = 1-2*side2
                 is3 = 1-2*side3          
                 do m3=1,numberOfGhostPoints  
                 do m2=1,numberOfGhostPoints  
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3-is3*m3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3+is3*m3; ! interior point             
                       u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0)
                     ! write(*,'("Vertex: set u(",3i3,")=u(",3i3,")*",f5.2)') j1,j2,j3,k1,k2,k3,symSign
                 end do
                 end do          
                 end do
               end if
             end if          
           end do
           end do
           end do
         end if
     else
         ! ---------------------------------
         ! --- assign corners and edges: ---
         ! ---------------------------------
         if( .false. .and. t.le.2.*dt )then
           write(*,'("assign GENERAL corner ghost points forcing=forcing, method=explicit, nd=",i3, " ***FINISH ME***")') nd
         end if
         if( nd.eq.2 )then
           ! ------ TWO DIMENSIONS ----
           if( debug.gt.1 .and. t.le.2*dt )then
             write(*,'("assign corner ghost points in 2D - general and curvilinear case")')
           end if
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
               if( ( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann ) )then
                 write(*,*) "Un-supported corner bc1=",bc1
                 stop 2222
               end if
               if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                 write(*,*) "Un-supported corner bc2 =",bc2
                 stop 2222
               end if
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign = +1. 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if   
               ! if( (bc(side1,0).eq.dirichlet .and. bc(side2,1).eq.neumann   ) .or. !     (bc(side1,0).eq.neumann   .and. bc(side2,1).eq.dirichlet ) ) then
               !   symSign=-1.;
               ! end if
               is1 = 1-2*side1
               is2 = 1-2*side2
               i1 = gridIndexRange(side1,0)
               i2 = gridIndexRange(side2,1)
               i3 = gridIndexRange(    0,2)
               if( mask(i1,i2,i3).gt.0 )then
                 do m2=1,numberOfGhostPoints
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point 
                     if( assignTwilightZone.eq.1 )then
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),0.,t,uc,ue1 )
                           call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),0.,t,uc,ue2 )
                         u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0) -symSign*ue2 + ue1 
                     else
                       ! finish me 
                     end if
                 end do
                 end do
               end if
             end if 
           end do
           end do
         else
           ! --- THREE DIMENSIONS ----
           ! ------ TWO DIMENSIONS ----
           if( debug.gt.1 .and. t.le.2*dt )then
             write(*,'("Assign corner and edge ghost points in 3D - general and curvilinear case")')
           end if
           ! write(*,'("assign general edge ghost points in 3D forcing=forcing")')
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
                 n3a=gridIndexRange(0,2)-extra
                 n3b=gridIndexRange(1,2)+extra
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side2,1)
               else if( edgeDirection.eq.1 )then
                 is2=0
                 n1a=gridIndexRange(side1,0)
                 n1b=gridIndexRange(side1,0)
                 n2a=gridIndexRange(    0,1)-extra
                 n2b=gridIndexRange(    1,1)+extra
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side1,0)
                 bc2=boundaryCondition(side3,2)
               else 
                 is1=0  
                 n1a=gridIndexRange(    0,0)-extra
                 n1b=gridIndexRange(    1,0)+extra
                 n2a=gridIndexRange(side2,1)
                 n2b=gridIndexRange(side2,1)
                 n3a=gridIndexRange(side3,2)
                 n3b=gridIndexRange(side3,2)
                 bc1=boundaryCondition(side2,1)
                 bc2=boundaryCondition(side3,2)
               end if
               if( bc1.gt.0 .and. bc2.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC )then
                 ! -- this is an edge between two physical boundaries --
                 if( bc1.ne.dirichlet .and. bc1.ne.exactBC .and. bc1.ne.neumann )then
                   write(*,*) "Un-supported edge bc1=",bc1
                   stop 2222
                 end if
                 if( bc2.ne.dirichlet .and. bc2.ne.exactBC .and. bc2.ne.neumann )then
                   write(*,*) "Un-supported edge bc2 =",bc2
                   stop 2222
                 end if
                 ! symSign = 1 : even symmetry
                 !          -1 : odd symmetry
                 symSign=1; 
                 if( bc1==dirichlet .or. bc1==exactBC )then
                   symSign = -symSign
                 end if
                 if( bc2==dirichlet .or. bc2==exactBC )then
                   symSign = -symSign
                 end if 
                 do i3=n3a,n3b
                 do i2=n2a,n2b
                 do i1=n1a,n1b
                   if( mask(i1,i2,i3).gt.0 )then
                     do m1=1,numberOfGhostPoints
                     do m2=1,numberOfGhostPoints 
                       if( edgeDirection==0 )then
                         ! edge lies along i1=const
                         j1 = i1; j2=i2-is2*m1; j3=i3-is3*m2; ! ghost 
                         k1 = i1; k2=i2+is2*m1; k3=i3+is3*m2; ! interior point              
                       else if( edgeDirection==1 )then
                         ! edge lies along i2=const
                         j1 = i1-is1*m1; j2=i2; j3=i3-is3*m2; ! ghost 
                         k1 = i1+is1*m1; k2=i2; k3=i3+is3*m2; ! interior point  
                       else
                        ! edge lies along i3=constant
                         j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3; ! ghost 
                         k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3; ! interior point  
                       end if
                         if( assignTwilightZone.eq.1 )then
                               call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),xy(j1,j2,j3,2),t,uc,ue1 )
                               call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),xy(k1,k2,k3,2),t,uc,ue2 )
                             u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0) - symSign*ue2 + ue1 
               ! write(*,'("Edge: ghost j1,j2,j3=",3i4," k1,k2,k3=",3i4," u=",e10.2," err=",e9.2)') j1,j2,j3,k1,k2,k3,u(j1,j2,j3,0),u(j1,j2,j3,0)-ue1
                         else
                           ! finish me 
                           write(*,*) "EdgeBC: finish me for forcing option"
                           stop 3333
                         end if
                     end do
                     end do     
                   end if ! if mask 
                 end do
                 end do
                 end do
               end if ! bc1>0 and bc2>0
             end do ! end do sideb
             end do ! end do sidea
           end do ! end edgeDirection
           do side3=0,1
           do side2=0,1
           do side1=0,1
             bc1 = bc(side1,0)
             bc2 = bc(side2,1)
             bc3 = bc(side3,2)
             if( bc1.gt.0 .and.  bc2.gt.0 .and. bc3.gt.0 .and. bc1.ne.exactBC .and. bc2.ne.exactBC .and. bc3.ne.exactBC )then
               ! --- Vertex where three physical boundaries meet ---
               ! symSign = 1 : even symmetry
               !          -1 : odd symmetry
               symSign=1; 
               if( bc1==dirichlet .or. bc1==exactBC )then
                 symSign = -symSign
               end if
               if( bc2==dirichlet .or. bc2==exactBC )then
                 symSign = -symSign
               end if        
               if( bc3==dirichlet .or. bc3==exactBC )then
                 symSign = -symSign
               end if  
               ! vertex is (i1,i2,i3)
               i1=gridIndexRange(side1,0)
               i2=gridIndexRange(side2,1)
               i3=gridIndexRange(side3,2)
               if( mask(i1,i2,i3).gt.0 )then
                 is1 = 1-2*side1
                 is2 = 1-2*side2
                 is3 = 1-2*side3          
                 do m3=1,numberOfGhostPoints  
                 do m2=1,numberOfGhostPoints  
                 do m1=1,numberOfGhostPoints
                   j1 = i1-is1*m1; j2=i2-is2*m2; j3=i3-is3*m3; ! ghost 
                   k1 = i1+is1*m1; k2=i2+is2*m2; k3=i3+is3*m3; ! interior point             
                     if( assignTwilightZone.eq.1 )then
                           call ogDeriv(ep,0,0,0,0,xy(j1,j2,j3,0),xy(j1,j2,j3,1),xy(j1,j2,j3,2),t,uc,ue1 )
                           call ogDeriv(ep,0,0,0,0,xy(k1,k2,k3,0),xy(k1,k2,k3,1),xy(k1,k2,k3,2),t,uc,ue2 )
                         u(j1,j2,j3,0) = symSign*u(k1,k2,k3,0) - symSign*ue2 + ue1 
                     else
                       ! finish me 
                       write(*,*) "SymBC: finish me for forcing option"
                       stop 3333
                     end if
                 end do
                 end do          
                 end do
               end if
             end if          
           end do
           end do
           end do
         end if
     end if
   else 
     ! **** OLD WAY *****
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
         orderOfExtrapolationForCorners= orderOfAccuracy+1
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
   return
   end
