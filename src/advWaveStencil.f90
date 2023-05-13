! This file automatically generated from advWaveStencil.bf90 with bpp.
!
! =======================================================================
! ============ Optimized advance routines for CgWave ====================
!              STENCIL VERSION
!
! April 2023 -- version 1
! =======================================================================
!


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



! =========================================================================================
! Macro: Compute the forcing for the update of u
! =========================================================================================




! ===========================================================================
!   EVALUATE COEFFICIENTS OF THE LAPLACIAN 
!         TWO DIMENSIONS
! ===========================================================================



! ===========================================================================
!   EVALUATE COEFFICIENTS OF THE LAPLACIAN 
!         THREE DIMENSIONS
! ===========================================================================





! order=6 stencil: 
! #Include 'meStencil2dOrder6.h'

! order=8 stencil: 
! #Include 'meStencil2dOrder8.h'

! ===========================================================================
!   Macro: eval stencil coeff at a given order
! ===========================================================================



! ===========================================================================
!   Macro: Main routine
!    EVALUATE STENCIL COEFFICIENTS OF THE MODIFIED EQUATION SCHEME
! ===========================================================================


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



! ==============================================
! ============== CONSTRUCT A FILE ==============
! ==============================================





! buildFile(advWaveStencil3dOrder8c,3,8,curvilinear)



            subroutine advWaveStencil( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,sc,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
!======================================================================
!   Advance a time step for Maxwells eqution
!     OPTIMIZED MODIFIED EQUATION VERSIONS
! 
!       ***** STENCIL VERSION *****
!
! nd : number of space dimensions
!
!======================================================================
            implicit none
            integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

            real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      ! real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
            real sc(0:*)
            real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
            real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
            real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds coeff of Laplacian for HA scheme

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



      ! write(*,*) 'Inside advWaveStencil...'

            gridType           =ipar(2)
            orderOfAccuracy    =ipar(3)

            if( orderOfAccuracy.eq.2 )then

                if( nd.eq.2 .and. gridType.eq.rectangular ) then
                    call advWaveStencil2dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,sc,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )

        ! else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
        !   call advWaveStencil2dOrder2c( ARGLIST() )
        ! else if( nd.eq.3 .and. gridType.eq.rectangular ) then
        !   call advWaveStencil3dOrder2r( ARGLIST() )
        ! else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
        !   call advWaveStencil3dOrder2c( ARGLIST() )
        ! else
        !   stop 2271
                end if
                stop 2222

            else if( orderOfAccuracy.eq.4 ) then
       ! if( nd.eq.2 .and. gridType.eq.rectangular )then
       !   call advWaveStencil2dOrder4r( ARGLIST() )
       !  else if(nd.eq.2 .and. gridType.eq.curvilinear )then
       !    call advWaveStencil2dOrder4c( ARGLIST() )
       !  else if(  nd.eq.3 .and. gridType.eq.rectangular )then
       !    call advWaveStencil3dOrder4r( ARGLIST() )
       !  else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       !    call advWaveStencil3dOrder4c( ARGLIST() )
       ! else
       !   stop 8843
       !end if
              stop 4444

            else if( orderOfAccuracy.eq.6 ) then

                if( nd.eq.2 .and. gridType.eq.rectangular )then
        !  call advWaveStencil2dOrder6r( ARGLIST() )
        ! else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        !   !call advWaveStencil2dOrder6c( ARGLIST() )
                else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        ! call advWaveStencil3dOrder6r( ARGLIST() )
        ! else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        !   !call advWaveStencil3dOrder6c( ARGLIST() )
              else
                  stop 8843
              end if

            else if( orderOfAccuracy.eq.8 ) then

                if( nd.eq.2 .and. gridType.eq.rectangular )then
        ! call advWaveStencil2dOrder8r( ARGLIST() )
        !  else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        !    call advWaveStencil2dOrder8c( ARGLIST() )
                else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        !  call advWaveStencil3dOrder8r( ARGLIST() )
        !  else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        !    call advWaveStencil3dOrder8c( ARGLIST() )
                else
                    stop 8843
                end if

            else
                write(*,'(" advWaveStencil:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy

                stop 11122
            end if

            return
            end








