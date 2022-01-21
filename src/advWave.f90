! This file automatically generated from advWave.bf90 with bpp.
!
! =======================================================================
! ============ Optimized advance routines for CgWave ====================
! =======================================================================
!
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

! 6th-order
! Use this next macro to declare the statement functions that are defined below
! To include derivatives of rx use OPTION=RX


! Define statement functions for difference approximations of order 6 
! To include derivatives of rx use OPTION=RX
! To include derivatives of rx use OPTION=RX

! 8th order
! Use this next macro to declare the statement functions that are defined below
! To include derivatives of rx use OPTION=RX


! Define statement functions for difference approximations of order 8 
! To include derivatives of rx use OPTION=RX
! To include derivatives of rx use OPTION=RX

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
! Macro: compute the coefficients in the sosup dissipation for curvilinear grids
! ===========================================================================================


! ===========================================================================================
! Macro: Output some debug info for the first few time-steps 
! ===========================================================================================

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (4th-order difference used with 2nd-order scheme) 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (6th-order difference used with 4th-order scheme)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! =========================================================================================
! Macro: Compute v=Delta(u) to second order for fourth-order scheme on curvilinear grids
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================
    
! =========================================================================================
! Macro: Compute the forcing for the update of u
! =========================================================================================

! =========================================================================================
!
! Macro: Advance the wave equation, EXPLICIT TIME-STEPPING
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   ORDERINTIME : 2 or 4 
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================


! =========================================================================================
! Macro: Compute the forcing for the IMPLICIT update of u
! =========================================================================================

! =========================================================================================
!
! Macro: FILL IN THE RHS FOR IMPLICIT TIME-STEPPING
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   ORDERINTIME : 2 or 4 
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================



! =========================================================================================
! Macro: ADD UPWIND DISSIPATION
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================


! =========================================================================================
! Macro: ADD UPWIND DISSIPATION FOR IMPLICIT TIME-STEPPING
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================

    
! =========================================================================================
! Macro: COMPUTE uDot = (un -um )
!
! ========================================================================================


! Argument list

! **********************************************************************************
! Macro ADV_WAVE:
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
!  GRIDTYPE : rectangular, curvilinear
! **********************************************************************************


  

    ! NOTE: For now 3D versions are just null versions below 



  !   ORDER=6 : BC's not implemented yet

  ! ORDER 8 -- finish BCs



            subroutine advWave( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
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

                if( nd.eq.2 .and. gridType.eq.rectangular ) then
                    call advWave2dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
                    call advWave2dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if( nd.eq.3 .and. gridType.eq.rectangular ) then
                    call advWave3dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
                    call advWave3dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else
                    stop 2271
                end if

            else if( orderOfAccuracy.eq.4 ) then
                if( nd.eq.2 .and. gridType.eq.rectangular )then
                    call advWave2dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if(nd.eq.2 .and. gridType.eq.curvilinear )then
                    call advWave2dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if(  nd.eq.3 .and. gridType.eq.rectangular )then
                    call advWave3dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
                    call advWave3dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
              else
                  stop 8843
              end if

            else if( orderOfAccuracy.eq.6 ) then

                if( nd.eq.2 .and. gridType.eq.rectangular )then
                    call advWave2dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if(nd.eq.2 .and. gridType.eq.curvilinear )then
                    call advWave2dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if(  nd.eq.3 .and. gridType.eq.rectangular )then
                    call advWave3dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
                    call advWave3dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
              else
                  stop 8843
              end if

            else if( orderOfAccuracy.eq.8 ) then

                if( nd.eq.2 .and. gridType.eq.rectangular )then
                    call advWave2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if(nd.eq.2 .and. gridType.eq.curvilinear )then
                    call advWave2dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if(  nd.eq.3 .and. gridType.eq.rectangular )then
                    call advWave3dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
                else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
                    call advWave3dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,bc,frequencyArray,ipar,rpar,ierr )
              else
                  stop 8843
              end if



            else
                write(*,'(" advWave:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy
          ! '
                stop 11122
            end if

            return
            end








