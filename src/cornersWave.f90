! This file automatically generated from cornersWave.bf90 with bpp.
! ==================================================================================
!
!        ASSIGN CORNER AND EDGE GHOST POINTS
!        -----------------------------------
!
! ==================================================================================


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



! define macros to evaluate higher derivatives (from maple/makeGetDerivativesMacros.maple)
!! ** June 13, 2023 : TURN OFF ??
! #Include "../maple/defineGetDerivativesMacros.h"

! NEW VERSION WITH DISTINCTIVE NAMES: *bug fixed: Nov 24, 2023
! #Include "../maple/defineGetThirdDerivativesMacros.h"




! From bcOptSmFOS.bf
! DataBase *pdb = &parameters.dbase;
! double precision pdb  ! pointer to data base
! ====================================================================
! Look up an integer parameter from the data base
! ====================================================================

! ====================================================================
! Look up a real parameter from the data base
! ====================================================================



! General begin loops macro











! ----- define extrapolation formulae ------











! ************************************************************************************************
!  This macro is used for looping over the faces of a grid to assign booundary conditions
!
! extra: extra points to assign
!          Case 1: extra=numberOfGhostPoints -- for assigning extended boundaries
!          Case 2: extra=-1 -- for assigning ghost points but not including extended boundaries
! numberOfGhostPoints : number of ghost points (1 for 2nd order, 2 for fourth-order ...)
!
!
! Output:
!  n1a,n1b,n2a,n2b,n3a,n3b : from gridIndexRange
!  nn1a,nn1b,nn2a,nn2b,nn3a,nn3b : includes "extra" points
! 
! ***********************************************************************************************









! =========================================================================
! Compute the normal on a curvilinear grid.
!
! Assumes is=1-2*side is defined. 
! =========================================================================

! ========================================================================================
!  Assign ghost points outside corners
! ========================================================================================




! ! ===================================================================================================
! ! Macro: Extrapolate Ghost Points 
! ! ORDER : 2,4,6,8
! ! ===================================================================================================
! #beginMacro extrapolateGhost(ORDER)

!   #If #ORDER eq "2"
!     #defineMacro EXTRAP(u,i1,i2,i3,n,is1,is2,is3) extrap3(u,i1,i2,i3,n,is1,is2,is3)
!   #Elif #ORDER eq "4"
!     #defineMacro EXTRAP(u,i1,i2,i3,n,is1,is2,is3) extrap5(u,i1,i2,i3,n,is1,is2,is3)
!   #Elif #ORDER eq "6"
!     ! #defineMacro EXTRAP(u,i1,i2,i3,n,is1,is2,is3) extrap6(u,i1,i2,i3,n,is1,is2,is3)
!     #defineMacro EXTRAP(u,i1,i2,i3,n,is1,is2,is3) extrap7(u,i1,i2,i3,n,is1,is2,is3)
!   #Elif #ORDER eq "8"
!     #defineMacro EXTRAP(u,i1,i2,i3,n,is1,is2,is3) extrap9(u,i1,i2,i3,n,is1,is2,is3)
!   #Else
!      write(*,'("corners: unexpected order=ORDER")') 
!     stop 1010
!   #End

!   ! ghost-loops will assign extra points in tangential directions
!   beginGhostLoops3d()
!     if( mask(i1,i2,i3).ne.0 )then
     
!       ! -- extrapolate ghost ---
!       do ghost=1,numGhost
!         j1=i1-is1*ghost
!         j2=i2-is2*ghost
!         j3=i3-is3*ghost

!         u(j1,j2,j3,uc) = EXTRAP(u,j1+is1,j2+is2,j3+is3,uc,is1,is2,is3) 

!       end do
       
!     end if ! mask .ne. 0
!   endLoops3d()
! #endMacro


! ! ============================================================================================
! ! Macro: evaluate the solution on the boundary for Dirichlet BCs
! ! ============================================================================================
! #beginMacro getDirichletForcing(ff,i1,i2,i3)
!   if( assignTwilightZone.eq.1 )then
!     ! compute RHS from TZ
!     if( nd.eq.2 )then
!       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
!     else
!       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
!     end if
!     ff = ue
!   else if( assignKnownSolutionAtBoundaries.eq.1 )then
!     ! -- we set inhomogeneous Dirichlet values for some known solutions 
!     if( knownSolutionOption.eq.planeWave )then
!       ! --- evaluate the plane wave solution ---
!       if( nd.eq.2 )then
!         ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
!       else
!         ff = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
!       end if 

!     else if( knownSolutionOption.eq.gaussianPlaneWave ) then
!       ! Eval the Gaussian plane wave solution
!       !    u = exp( -beta*(xi^2) )*cos( k0*xi )
!       !    xi = kx*( x-x0) + ky*(y-y0) + kz*(z-z0) - c*t       
!       !  
!       if( nd.eq.2 )then
!         xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) - c*t
!       else
!         xi = kxGPW*(xy(i1,i2,i3,0)-x0GPW) + kyGPW*(xy(i1,i2,i3,1)-y0GPW) + kzGPW*(xy(i1,i2,i3,2)-z0GPW) - c*t
!       end if 
!       ff = exp( -betaGPW*xi**2 ) * cos( k0GPW*xi )      

!     else if( knownSolutionOption.eq.boxHelmholtz ) then
!       ! --- evaluate the boxHelmholtz solution ---
!       ! For multi-freq we add up all the component frequencies
!       ff = 0. 
!       do freq=0,numberOfFrequencies-1
!         ! kx = kxBoxHelmholtz + twoPi*freq
!         ! ky = kyBoxHelmholtz + twoPi*freq
!         ! kz = kzBoxHelmholtz + twoPi*freq
!         ! coswt = cos( frequencyArray(freq)*t )

!         getBoxHelmholtzParameters( freq, omega, kx,ky,kz, F90 )
!         coswt = cos( omega*t )
!         if( nd.eq.2 )then
!           ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) *coswt
!         else
!           ff = ff + sin( kx*xy(i1,i2,i3,0) ) * sin( ky*xy(i1,i2,i3,1) ) * sin( kz*xy(i1,i2,i3,2) ) *coswt
!         end if 
!       end do

!     else if( knownSolutionOption.eq.polyPeriodic ) then
!       ! --- evaluate the polyPeriodic solution ---
!       if( nd.eq.2 )then
!         ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1)                                 ) *coswt
!       else
!         ff = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
!       end if 

!     else
!       stop 9876
!     end if 
!   end if
! #endMacro      




! ! ============================================================================================
! ! Macro: evaluate the solution on the boundary for Dirichlet BCs
! !    Solving
! !        u_tt = c^2 * Delta( u ) + f(x,y,z,t)
! !        u = g
! ! For TZ at order=2:
! !    ff = ue_tt - c^2*Delta(ue) 
! !    gtt = g_tt = uett
! !  order=4:
! !    fLap 
! !    ftt
! !    gtttt
! ! ============================================================================================
! #beginMacro getDirichletCompatibilityForcing(FORCING,DIM,ORDER,ff,gtt,fLap,ftt,gtttt)
!   #If #FORCING eq "noForcing"
!     ! No forcing, do nothing 
!   #Else
!     if( assignTwilightZone.eq.1 )then
!       ! compute RHS from TZ
!       #If #DIM eq 2
!         call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexx )
!         call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyy )
!         ueLap = uexx + ueyy
!         call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uett )
!         #If #ORDER eq "4"
!           call ogDeriv(ep,2,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uettxx )
!           call ogDeriv(ep,2,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uettyy )
!           uettLap = uettxx + uettyy

!           call ogDeriv(ep,0,4,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxxx )
!           call ogDeriv(ep,0,2,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxyy )
!           call ogDeriv(ep,0,0,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyyy )

!           ueLap2 = uexxxx + 2.*uexxyy + ueyyyy ! Lap^2( u )

!           call ogDeriv(ep,4,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uetttt )

!         #End
!       #Else
!         ! 3D 
!         call ogDeriv(ep,0,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexx )
!         call ogDeriv(ep,0,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyy )
!         call ogDeriv(ep,0,0,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezz )
!         ueLap = uexx + ueyy + uezz
!         call ogDeriv(ep,2,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uett )
!         #If #ORDER eq "4"
!           call ogDeriv(ep,2,2,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettxx )
!           call ogDeriv(ep,2,0,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettyy )
!           call ogDeriv(ep,2,0,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettzz )
!           uettLap = uettxx + uettyy + uettzz

!           call ogDeriv(ep,0,4,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxxx )
!           call ogDeriv(ep,0,0,4,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyyy )
!           call ogDeriv(ep,0,0,0,4,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezzzz )
!           call ogDeriv(ep,0,2,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxyy )
!           call ogDeriv(ep,0,2,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxzz )
!           call ogDeriv(ep,0,0,2,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyzz )

!           ueLap2 = uexxxx + 2.*( uexxyy + uexxzz + ueyyzz)  + ueyyyy + uezzzz ! Lap^2( u )

!           call ogDeriv(ep,4,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uetttt )

!         #End      

!       #End
!       ff = uett - c2*ueLap
!       gtt = uett

!       #If #ORDER eq "4"
!         fLap = uettLap - c2*ueLap2 
!         ftt  = uetttt - c2*uettLap
!         gtttt = uetttt
!       #End    


!     else
!       ff=0.; gtt=0.; fLap=0.; ftt=0.; gtttt=0.;
!     end if
!   #End
! #endMacro      





! ===================================================================================================
! Macro: get loop bounds over boundaries with extram points in tangential directions 
! ===================================================================================================




! ===================================================================================================
! Macro: Get the TZ solution in 2d or 3d 
! ===================================================================================================



! ! ------------------------------------------------------------------------------------
! !  Macro: evaluate the RHS to the Neumann BC
! ! ------------------------------------------------------------------------------------
! #beginMacro getNeumannForcing(ff)
!   if( assignTwilightZone.eq.1 )then
!     ! compute RHS from TZ
!     if( nd.eq.2 )then
!       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ue )
!       call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uex)
!       call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uey)
!       ff = a0*ue + a1*( an1*uex + an2*uey )
!     else
!       call ogDeriv(ep,0,0,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ue )
!       call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uex)
!       call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uey)
!       call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uez)
!       ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
!     end if

!   else if( assignKnownSolutionAtBoundaries.eq.1 )then

!     ! -- we set inhomogeneous Neumann values for some known solutions 
!     if( knownSolutionOption.eq.planeWave )then
!       ! --- evaluate RHS for the plane wave solution ---
!       if( nd.eq.2 )then
!         ue    = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
!         cosPW = ampPlaneWave*cos( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) - omegaPlaneWave*t )
!         uex   = kxPlaneWave*cosPW
!         uey   = kyPlaneWave*cosPw

!         ff = a0*ue + a1*( an1*uex + an2*uey )
!       else
!         ue    = ampPlaneWave*sin( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
!         cosPW = ampPlaneWave*cos( kxPlaneWave*xy(i1,i2,i3,0) + kyPlaneWave*xy(i1,i2,i3,1) + kzPlaneWave*xy(i1,i2,i3,2) - omegaPlaneWave*t )
!         uex   = kxPlaneWave*cosPW
!         uey   = kyPlaneWave*cosPw
!         uez   = kzPlaneWave*cosPw

!         ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
!       end if 

!     else if( knownSolutionOption.eq.gaussianPlaneWave )then
!       ! Do nothing for Gaussian plane wave solution for now
!       ff = 0.
      
!     else if( knownSolutionOption.eq.boxHelmholtz ) then
!       ! --- evaluate RHS the boxHelmholtz solution ---
!       if( nd.eq.2 )then
!         ue  = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) * sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) *coswt
!         uex = cos( kxBoxHelmholtz*xy(i1,i2,i3,0) ) * sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) *coswt * kxBoxHelmholtz
!         uey = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) * cos( kyBoxHelmholtz*xy(i1,i2,i3,1) ) *coswt * kyBoxHelmholtz

!         ff = a0*ue + a1*( an1*uex + an2*uey )
!       else
!         ue  = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * sin( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt
!         uex = cos( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * sin( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt * kxBoxHelmholtz
!         uey = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *cos( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * sin( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt * kyBoxHelmholtz
!         uez = sin( kxBoxHelmholtz*xy(i1,i2,i3,0) ) *sin( kyBoxHelmholtz*xy(i1,i2,i3,1) ) * cos( kzBoxHelmholtz*xy(i1,i2,i3,2) ) *coswt * kzBoxHelmholtz

!         ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
!       end if

!     else if( knownSolutionOption.eq.polyPeriodic ) then
!       ! --- evaluate RHS the polyPeriodic solution ---
!       if( nd.eq.2 )then
!         ue  = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) ) *coswt
!         uex = (      a1PolyPeriodic                                                            ) *coswt
!         uey = (                          b1PolyPeriodic                                        ) *coswt

!         ff = a0*ue + a1*( an1*uex + an2*uey )
!       else
!         ue  = ( a0PolyPeriodic + a1PolyPeriodic*xy(i1,i2,i3,0) + b1PolyPeriodic*xy(i1,i2,i3,1) + c1PolyPeriodic*xy(i1,i2,i3,2) ) *coswt
!         uex = (      a1PolyPeriodic                                                                                            ) *coswt
!         uey = (                          b1PolyPeriodic                                                                        ) *coswt
!         uez = (                                              c1PolyPeriodic                                                    ) *coswt

!         ff = a0*ue + a1*( an1*uex + an2*uey + an3*uez )
!       end if


!     else

!       stop 9876

!     end if 

!   end if
! #endMacro 


! ! ============================================================================================
! ! Macro: evaluate the forcings for Neumann CBCs
! !    Solving
! !        u_tt = c^2 * Delta( u ) + f(x,y,z,t)
! !        u.n = g
! ! For TZ at order=2:
! !    gg = g
! !  order=4:
! !    gg,
! !    nDotGradF = n.grad( f ), f = ue_tt - c^2*lap(ue)
! !    gtt
! ! ============================================================================================
! #beginMacro getNeumannCompatibilityForcing(FORCING,DIM,ORDER,gg,nDotGradF,gtt )
!   #If #FORCING eq "noForcing"
!     ! No forcing, do nothing 
!   #Else
!     if( assignTwilightZone.eq.1 )then
!       ! compute RHS from TZ
!       #If #DIM eq 2
!         call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uex )
!         call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uey )
!         gg = an1*uex + an2*uey

!         #If #ORDER eq "4"
!           call ogDeriv(ep,2,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uettx )
!           call ogDeriv(ep,2,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uetty )
!           gtt = an1*( uettx ) + an2*( uetty )

!           call ogDeriv(ep,0,3,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxx )
!           call ogDeriv(ep,0,1,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexyy )
!           call ogDeriv(ep,0,2,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,uexxy )
!           call ogDeriv(ep,0,0,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t,uc,ueyyy )

!           nDotGradF = an1*( uettx - c2*( uexxx + uexyy ) ) + !                       an2*( uetty - c2*( uexxy + ueyyy ) )

!         #End

!       #Else
!         ! ----- 3D  -----
!         call ogDeriv(ep,0,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uex )
!         call ogDeriv(ep,0,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uey )
!         call ogDeriv(ep,0,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uez )
!         gg = an1*uex + an2*uey + an3*uez

!         #If #ORDER eq "4"
!           call ogDeriv(ep,2,1,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettx )
!           call ogDeriv(ep,2,0,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uetty )
!           call ogDeriv(ep,2,0,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uettz )
!           gtt = an1*( uettx ) + an2*( uetty ) + an3*( uettz )

!           call ogDeriv(ep,0,3,0,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxx )
!           call ogDeriv(ep,0,2,1,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxy )
!           call ogDeriv(ep,0,1,2,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexyy )
!           call ogDeriv(ep,0,0,3,0,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyy )
!           call ogDeriv(ep,0,2,0,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexxz )
!           call ogDeriv(ep,0,1,0,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uexzz )
!           call ogDeriv(ep,0,0,2,1,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyyz )
!           call ogDeriv(ep,0,0,1,2,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,ueyzz )
!           call ogDeriv(ep,0,0,0,3,xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t,uc,uezzz )

!           nDotGradF = an1*( uettx - c2*( uexxx + uexyy + uexzz ) ) + !                       an2*( uetty - c2*( uexxy + ueyyy + ueyzz ) ) + !                       an3*( uettz - c2*( uexxz + ueyyz + uezzz ) )
!         #End      

!       #End
!     else
!       gg=0.;  gtt=0.; nDotGradF=0.; 
!     end if
!   #End
! #endMacro      


! ========================================================================================
!  Apply symmetry conditions for ghost along edges in 3D 
!
!    ASSIGN EDGES WITH DIRICHLET AND/OR NEUMANN BUT NOT EXACT BCs
! ========================================================================================


! ========================================================================================
!  Apply symmetry conditions in corner ghost for Cartesian grids
!
!    ASSIGN CORNERS WITH DIRICHLET AND/OR NEUMANN BUT NOT EXACT BCs 
! ========================================================================================



! ========================================================================================
!  Apply general conditions for ghost along edges in 3D 
!  
!    ASSIGN EDGES WITH DIRICHLET AND/OR NEUMANN BUT NOT EXACT BCs
!
! METHOD: explicit or implicit (time-stepping method)
! ========================================================================================


! ========================================================================================
!  Apply general corner and edge conditions 
!
! METHOD: explicit or implicit (time-stepping method)
!
! ========================================================================================






! Argument list

! **********************************************************************************
! Macro cornerWave...
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
! **********************************************************************************

! --- Macro to build the file for each dimension and order ---

! --- construct the different files ----



! buildFile(corners2dOrder4,2,6)
! buildFile(corners3dOrder4,3,6)




subroutine cornersWave( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr ) 
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

  real uTemp(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
  real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)

  integer gridIndexRange(0:1,0:2),boundaryCondition(0:1,0:2), dimRange(0:1,0:2), isPeriodic(0:*)
  real frequencyArray(0:*)

  double precision pdb  ! pointer to data base


  integer ipar(0:*)
  real rpar(0:*)

  integer orderOfAccuracy

  ! extract parameters we need: 
  orderOfAccuracy  = ipar( 4)

  if( nd.eq.2 )then
    if( orderOfAccuracy.eq.2 )then
      call cornersWave2dOrder2( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    elseif( orderOfAccuracy.eq.4 )then
      call cornersWave2dOrder4( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    else
      stop 6666
    end if
  else
    if( orderOfAccuracy.eq.2 )then
      call cornersWave2dOrder2( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    elseif( orderOfAccuracy.eq.4 )then
      call cornersWave2dOrder4( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    else
      stop 7777
    end if    

  end if

  return
  end

