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

! define macros to evaluate derivatives for the 6th order method (from maple/makeGetDerivativesMacros.maple)
! ****** File written by makeGetDerivativesMacros.maple  ******



! =======================================================
!  Macro to compute Sixth derivatives in 2 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================

! =======================================================
!  Macro to compute Sixth derivatives in 3 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================



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


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (8th-order difference used with 6th-order scheme)
! (x+y)^8 = x^8 + 8*x^7*y + 28*x^6*y^2 + 56*x^5*y^3 + 70*x^4*y^4 + 56*x^3*y^5 + 28*x^2*y^6 + 8*x*y^7 + y^8
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Upwind (sosup) dissipation (10th-order difference used with 8th-order scheme)
! (x+y)^10 = x^10 + 10*x^9*y + 45*x^8*y^2 + 120*x^7*y^3 + 210*x^6*y^4 + 252*x^5*y^5 + 210*x^4*y^6 + 120*x^3*y^7 + 45*x^2*y^8 + 10*x*y^9 + y^10
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
!
! Macro: Advance the wave equation, EXPLICIT TIME-STEPPING and 
!
!     +++++++++ SUPERGRID ++++++++++++++
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

! ==================================================================
! ==================================================================

! **********************************************************************************
! Macro ADV_WAVE:
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
!  GRIDTYPE : rectangular, curvilinear
! **********************************************************************************


! Macro to build separate files 


  ! ==== Build separate f90 files for different cases ====





  ! !   ORDER=6 : BC's not implemented yet -- needed for upwinding, SuperGrid

  ! ! ORDER 8 -- needed for upwinding,  SuperGrid


      subroutine advWave( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
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
  ! real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
  real stencilCoeff(0:*)   ! holds stencil coeffs
  real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1) 
  real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
  real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds coeff of Laplacian for HA scheme
  real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
  real etax(nd1a:nd1b,0:*)  ! superGrid functions
  real etay(nd2a:nd2b,0:*)
  real etaz(nd3a:nd3b,0:*)
  integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
  integer bc(0:1,0:2),ierr  
  real frequencyArray(0:*)
  integer ipar(0:*)
  real rpar(0:*)

      ! integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b

      ! real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      ! real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      ! real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      ! real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      ! ! real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
      ! real stencilCoeff(0:*)
      ! real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
      ! real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
      ! real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds current Helmholtz solutions

      ! real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)
      ! real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)

      ! integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
      ! integer bc(0:1,0:2),ierr

      ! integer ipar(0:*)
      ! real rpar(0:*)

      ! real frequencyArray(0:*)
      
!     ---- local variables -----
  integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime
  integer addForcing,orderOfDissipation,option,addUpwinding,modifiedEquationApproach
  ! integer useImplicitUpwindDissipation,useUpwindDissipation
  integer useWhereMask,solveForE,solveForH,grid
  integer ex,ey,ez, hx,hy,hz

  integer rectangular,curvilinear
  parameter( rectangular=0, curvilinear=1 )

  integer standardME, hierarchicalME, stencilME
  parameter( standardME=0, hierarchicalME=1, stencilME=2 )
!...........end   statement functions



       
  option                    = ipar( 0)
  gridType                  = ipar( 2)
  orderOfAccuracy           = ipar( 3)
  orderInTime               = ipar( 4)
  modifiedEquationApproach  = ipar(18)

  ! write(*,*) 'Inside advWave...'
  ! write(*,*) 'option, orderOfAccuracy, modifiedEquationApproach=',option, orderOfAccuracy, modifiedEquationApproach

        ! useUpwindDissipation         = ipar(11)  ! explicit upwind dissipation
  ! useImplicitUpwindDissipation = ipar(12)  ! true if upwind-dissipation is on for impliciit time-stepping 

  if( option.eq.1 )then 
   addUpwinding = 1
  else
   addUpwinding = 0
  end if


  if( orderOfAccuracy.eq.2 )then

    if( modifiedEquationApproach.eq.standardME .or. addUpwinding.ne.0 )then
      ! standard ME scheme, or upwind stage
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWave2dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWave2dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWave3dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWave3dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if
    else if( modifiedEquationApproach.eq.hierarchicalME ) then
      ! new Hierarchical scheme
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveME2dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveME2dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWaveME3dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWaveME3dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if 
    else if( modifiedEquationApproach.eq.stencilME ) then

     if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveStencil2dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveStencil2dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
       call advWaveStencil3dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       call advWaveStencil3dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if 

    else  
      write(*,'("Unknown modifiedEquationApproach=",i6)') modifiedEquationApproach
      stop 1111
    end if

    ! if( nd.eq.2 .and. gridType.eq.rectangular ) then
    !   call advWave2dOrder2r( ARGLIST() )
    ! else if( nd.eq.2 .and. gridType.eq.curvilinear ) then
    !   call advWave2dOrder2c( ARGLIST() )
    ! else if( nd.eq.3 .and. gridType.eq.rectangular ) then
    !   call advWave3dOrder2r( ARGLIST() )
    ! else if( nd.eq.3 .and. gridType.eq.curvilinear ) then
    !   call advWave3dOrder2c( ARGLIST() )
    ! else
    !   stop 2271
    ! end if

  else if( orderOfAccuracy.eq.4 ) then

    if( modifiedEquationApproach.eq.standardME .or. addUpwinding.ne.0 )then
      ! standard ME scheme, or upwind stage
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWave2dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWave2dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWave3dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWave3dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if
    else if( modifiedEquationApproach.eq.hierarchicalME ) then
      ! new Hierarchical scheme
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveME2dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveME2dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWaveME3dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWaveME3dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if 

    else if( modifiedEquationApproach.eq.stencilME ) then

     if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveStencil2dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveStencil2dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWaveStencil3dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       call advWaveStencil3dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if 

    else  
      write(*,'("Unknown modifiedEquationApproach=",i6)') modifiedEquationApproach
      stop 1111
    end if

  else if( orderOfAccuracy.eq.6 ) then

    if( modifiedEquationApproach.eq.standardME .or. addUpwinding.ne.0 )then
      ! standard ME scheme, or upwind stage
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWave2dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWave2dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWave3dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWave3dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if
    else if( modifiedEquationApproach.eq.hierarchicalME ) then
      ! new Hierarchical scheme
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveME2dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveME2dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWaveME3dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWaveME3dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if 
    else if( modifiedEquationApproach.eq.stencilME ) then

     if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveStencil2dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveStencil2dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
       call advWaveStencil3dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       call advWaveStencil3dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if 

    else  
      write(*,'("Unknown modifiedEquationApproach=",i6)') modifiedEquationApproach
      stop 1111
    end if


  else if( orderOfAccuracy.eq.8 ) then

  if( modifiedEquationApproach.eq.standardME .or. addUpwinding.ne.0 )then
      ! standard ME scheme, or upwind stage
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWave2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWave2dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWave3dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWave3dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if
    else if( modifiedEquationApproach.eq.hierarchicalME ) then
      ! new Hierarchical scheme
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveME2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveME2dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWaveME3dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWaveME3dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else
       stop 8843
      end if 

    else if( modifiedEquationApproach.eq.stencilME ) then


     if( nd.eq.2 .and. gridType.eq.rectangular )then
       call advWaveStencil2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
       call advWaveStencil2dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
       call advWaveStencil3dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       ! call advWaveStencil3dOrder8c( ARGLIST() )
      else
       stop 8843
      end if                  
    else  
      write(*,'("Unknown modifiedEquationApproach=",i6)') modifiedEquationApproach
      stop 1111
    end if


   !  if( nd.eq.2 .and. gridType.eq.rectangular )then
   !    if( addUpwinding==0 )then
   !      ! new ME scheme
   !      call advWaveME2dOrder8r( ARGLIST() )
   !    else
   !      call advWave2dOrder8r( ARGLIST() )
   !    end if

   !  else if(nd.eq.2 .and. gridType.eq.curvilinear )then
   !   if( addUpwinding==0 )then
   !      ! new ME scheme
   !      call advWaveME2dOrder8c( ARGLIST() )
   !    else
   !      call advWave2dOrder8c( ARGLIST() )
   !    end if  

   !  else if(  nd.eq.3 .and. gridType.eq.rectangular )then
   !    call advWave3dOrder8r( ARGLIST() )
   !  else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
   !    call advWave3dOrder8c( ARGLIST() )
   ! else
   !   stop 8843
   ! end if



  else
    write(*,'(" advWave:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy
      ! '
    stop 11122
  end if

  return
  end








