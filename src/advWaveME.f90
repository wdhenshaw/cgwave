! This file automatically generated from advWaveME.bf90 with bpp.
!
! =======================================================================
! ============ Optimized advance routines for CgWave ====================
!              NEW MODIFIED EQUATION VERSIONS
!
! Feb 2022 : first version created from advWave.bf90 and hierDeriv.bf90
! =======================================================================
!
! These next include files will define the macros that will define the difference approximations
! The actual macro is called below
! #Include "defineDiffOrder2f.h"
! #Include "defineDiffOrder4f.h"

! ! 6th-order
! #Include "defineDiffOrder6f.h"

! ! 8th order
! #Include "defineDiffOrder8f.h"

! define macros to evaluate derivatives for the 6th order method (from maple/makeGetDerivativesMacros.maple)
! #Include "../maple/defineGetSixthDerivativesMacros.h"


! ======================================================================================
!   Evaluate the TZ exact solution in 2D
! ======================================================================================

! ======================================================================================
!   Evaluate the TZ exact solution in 3D
! ======================================================================================

    


! ---------------------------------------------------------------------------
! Macro : beginLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : beginLoops
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
! Macro : endLoopsMask
! ---------------------------------------------------------------------------




! ===========================================================================================
! Macro: Output some debug info for the first few time-steps 
! ===========================================================================================

! =============================================================================
!  Get coefficients in the derivatives on a CURVILINEAR GRID in TWO DIMENSIONS
!
! DERIV : 2,3,4,5,6,7,8
!
! These coefficients are from chainRuleCoefficients.maple
! =============================================================================








    
! =========================================================================================
! Macro: Compute the forcing for the update of u
! =========================================================================================


! =========================================================================
! Compute errors in the derivatives and save the derivatives and errors
! =========================================================================



! =============================================================================
!  Evaluate derivatives using the hierarchical approach on a rectangular grid
!     **** THIS IS A TEST ROUTINE **** 
! =============================================================================



! =============================================================================
!  SIXTH ORDER ME SCHEME
!       RECTANGULAR GRID 
! =============================================================================


! ===================================================================================================
!  SIXYTH ORDER ME 
!         CURVILINEAR grid 
! ===================================================================================================

! =============================================================================
!  EIGHTH ORDER ME SCHEME
!       RECTANGULAR GRID 
! =============================================================================


! ===================================================================================================
!  EIGHTH ORDER ME 
!         CURVILINEAR grid 
! ===================================================================================================

! =========================================================================================
!
! Macro: Advance the wave equation, EXPLICIT TIME-STEPPING
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   ORDERINTIME : 2 or 4 
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================




! Argument list

! **********************************************************************************
! Macro ADV_WAVE:
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
!  GRIDTYPE : rectangular, curvilinear
! **********************************************************************************


  


  ! buildFile(advWaveME2dOrder2r,2,2,rectangular)
  ! buildFile(advWaveME3dOrder2r,3,2,rectangular)
  ! buildFile(advWaveME2dOrder2c,2,2,curvilinear)
  ! buildFile(advWaveME3dOrder2c,3,2,curvilinear)
  ! buildFile(advWaveME2dOrder4r,2,4,rectangular)
  ! buildFile(advWaveME3dOrder4r,3,4,rectangular)
  ! buildFile(advWaveME2dOrder4c,2,4,curvilinear)
  ! buildFile(advWaveME3dOrder4c,3,4,curvilinear)
  ! buildFile(advWaveME3dOrder6r,3,6,rectangular)
  ! buildFile(advWaveME3dOrder6c,3,6,curvilinear)
  ! buildFile(advWaveME3dOrder8r,3,8,rectangular)
  ! buildFile(advWaveME3dOrder8c,3,8,curvilinear)



            subroutine advWaveME( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
!======================================================================
!   Advance a time step for Maxwells eqution
!     OPTIMIZED version for rectangular grids.
! nd : number of space dimensions
!
! ipar(0)  = option : option=0 - Maxwell+Artificial diffusion
!                           =1 - AD only
!======================================================================
            implicit none
            integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

            real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
            real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions

            real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
            real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

            integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
            integer bc(0:1,0:2),ierr

            integer ipar(0:*)
            real rpar(0:*)

            real frequencyArray(0:*)
            
!     ---- local variables -----
            integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime
            integer addForcing,orderOfDissipation,option
            integer useWhereMask,solveForE,solveForH,grid
            integer ex,ey,ez, hx,hy,hz

            integer rectangular,curvilinear
            parameter( rectangular=0, curvilinear=1 )
!...........end   statement functions


      ! write(*,*) 'Inside advWave...'

            gridType           =ipar(2)
            orderOfAccuracy    =ipar(3)

            if( orderOfAccuracy.eq.2 )then

        ! if( nd.eq.2 .and. gridType.eq.rectangular ) then
        !   call advWaveME2dOrder2r( ARGLIST() )
        ! else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
        !   call advWaveME2dOrder2c( ARGLIST() )
        ! else if( nd.eq.3 .and. gridType.eq.rectangular ) then
        !   call advWaveME3dOrder2r( ARGLIST() )
        ! else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
        !   call advWaveME3dOrder2c( ARGLIST() )
        ! else
        !   stop 2271
        ! end if
                stop 2222

            else if( orderOfAccuracy.eq.4 ) then
       !  if( nd.eq.2 .and. gridType.eq.rectangular )then
       !    call advWaveME2dOrder4r( ARGLIST() )
       !  else if(nd.eq.2 .and. gridType.eq.curvilinear )then
       !    call advWaveME2dOrder4c( ARGLIST() )
       !  else if(  nd.eq.3 .and. gridType.eq.rectangular )then
       !    call advWaveME3dOrder4r( ARGLIST() )
       !  else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       !    call advWaveME3dOrder4c( ARGLIST() )
       ! else
       !   stop 8843
       ! end if
              stop 4444

            else if( orderOfAccuracy.eq.6 ) then

                if( nd.eq.2 .and. gridType.eq.rectangular )then
                    call advWaveME2dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
        ! else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        !   !call advWaveME2dOrder6c( ARGLIST() )
        ! else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        !   !call advWaveME3dOrder6r( ARGLIST() )
        ! else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        !   !call advWaveME3dOrder6c( ARGLIST() )
              else
                  stop 8843
              end if

            else if( orderOfAccuracy.eq.8 ) then

              if( nd.eq.2 .and. gridType.eq.rectangular )then
                  call advWaveME2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
       !  else if(nd.eq.2 .and. gridType.eq.curvilinear )then
       !    call advWaveME2dOrder8c( ARGLIST() )
       !  else if(  nd.eq.3 .and. gridType.eq.rectangular )then
       !    call advWaveME3dOrder8r( ARGLIST() )
       !  else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       !    call advWaveME3dOrder8c( ARGLIST() )
              else
                  stop 8843
              end if

              stop 8888 

            else
                write(*,'(" advWaveME:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy

                stop 11122
            end if

            return
            end








