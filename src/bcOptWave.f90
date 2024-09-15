! This file automatically generated from bcOptWave.bf90 with bpp.
! ==================================================================================
!
!        Optimized Assign Boundary Conditions for CgWave
!        -----------------------------------------------
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
! ****** File written by makeGetDerivativesMacros.maple  ******



! =======================================================
!  Macro to compute Third derivatives in 2 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================

! =======================================================
!  Macro to compute Third derivatives in 3 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================



! define macros to evaluate derivatives for the 6th order method (from maple/makeGetDerivativesMacros.maple)
!! turned off May 4, 2023
!! #Include "../maple/defineGetSixthDerivativesMacros.h"


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




! ===================================================================================================
! Macro: Extrapolate Ghost Points 
!  WIDTH    forumala 
!   2         2 -1 
!   3         3 -3 1
!   4         4 -6 4 1
!   5         5 -10 10 -5 1
!   6
! ===================================================================================================


! ============================================================================================
! Macro: evaluate the solution on the boundary for Dirichlet BCs
! ============================================================================================

! ===================================================================================================
! Macro: Assign boundary and ghost points on Dirichlet boundaries 
! ORDER : 2,4,6,8
! ===================================================================================================


! ============================================================================================
! Macro: evaluate the solution on the boundary for Dirichlet BCs
!    Solving
!        u_tt = c^2 * Delta( u ) + f(x,y,z,t)
!        u = g
! For TZ at order=2:
!    ff = ue_tt - c^2*Delta(ue) 
!    gtt = g_tt = uett
!  order=4:
!    fLap 
!    ftt
!    gtttt
! ============================================================================================


! ===================================================================================================
! Macro: Assign the boundary values on Dirichlet boundaries
!
!  NOTE: DIM AND GRIDTYPE NOT CURRENTLY USED
!  FORCE : USEFORCING or NOFORCING 
! ===================================================================================================

! ===================================================================================================
! Macro: Assign the boundary values on "exact" boundaries
!
! ===================================================================================================



! ===================================================================================================
! Macro: get loop bounds over boundaries with extram points in tangential directions 
! ===================================================================================================

! ===================================================================================================
! Macro: get loop bounds over boundaries with extram points in tangential directions 
!   Do NOT include the end grid points when adjacent sides are Dirichlet BCs
!   since we assume the extended dirichlet boundary has been assigned
!
!                |                         |
!   BC=Neumann   +                         + BC=Dirichlet
!                |                         |
!          G--G--B--+--+--+-- ...  --+--+--B--G--G
!            m1a |                     m1b |
!                X                         X
!
! ===================================================================================================


! ===================================================================================================
! Macro: Assign ghost points on Dirichlet boundaries using: 
!          *** COMPATIBILITY BOUNDARY CONDITIONS ****
! ORDER : 2,4,6,8
! FORCING : noForcing, forcing
! ===================================================================================================

! ===================================================================================================
! Macro: Get the TZ solution in 2d or 3d 
! ===================================================================================================

! ===================================================================================================
! Macro: Evaluate the forcing for the CBCs
! ===================================================================================================

! ===================================================================================================
! Macro: Assign ghost points on Dirichlet boundaries using: 
!          *** COMPATIBILITY BOUNDARY CONDITIONS ****
! ORDER : 6 
! FORCING : noForcing, forcing
! ===================================================================================================




! ------------------------------------------------------------------------------------
!  Macro: evaluate the RHS to the Neumann BC
! ------------------------------------------------------------------------------------

! ===================================================================================================
! Macro: Assign ghost points on Neumann boundaries 
! ORDER : 2,4,6,8
! ===================================================================================================





! ============================================================================================
! Macro: evaluate the forcings for Neumann CBCs
!    Solving
!        u_tt = c^2 * Delta( u ) + f(x,y,z,t)
!        u.n = g
! For TZ at order=2:
!    gg = g
!  order=4:
!    gg,
!    nDotGradF = n.grad( f ), f = ue_tt - c^2*lap(ue)
!    gtt
! ============================================================================================

!==========================================================================
!  Check the coefficients in the ghost points of the residual equations
! using the discrete delta approach
!==========================================================================

! ===================================================================================================
! Macro: Assign ghost points on Neumann boundaries using: 
!          *** COMPATIBILITY BOUNDARY CONDITIONS ****
! ORDER : 2,4,6,8
! FORCING : noForcing, forcing
! ===================================================================================================

! ========================================================================================
!  Apply symmetry conditions for ghost along edges in 3D 
! ========================================================================================


! ========================================================================================
!  Apply symmetry conditions in corner ghost for Cartesian grids 
! ========================================================================================





! ===================================================================================
! Utility macro to call different versions of assignDirichletGhostCompatibility
! for a given ORDER
! ===================================================================================

! ===================================================================================
! Utility macro to call different versions of assignNeumannGhostCompatibility
! for a given ORDER
! ===================================================================================


! =========================================================================
! Macro: 
!   Assign the RHS for implicit boundary conditions
!
! =========================================================================

! ===============================================================================
! ----------- STAGE I : Assign Dirichlet Conditions -------------
! ===============================================================================


! Argument list

! **********************************************************************************
! Macro BC_WAVE:
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
! **********************************************************************************



! --- Macro to build the file for each dimension and order ---

! ****************************************************************
! --- construct the different files ----
! ****************************************************************



! buildFile(bcOptWave2dOrder4,2,6)
! buildFile(bcOptWave3dOrder4,3,6)


subroutine bcOptWave( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr ) 
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
      call bcOptWave2dOrder2( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    elseif( orderOfAccuracy.eq.4 )then
      call bcOptWave2dOrder4( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    else
      stop 6666
    end if
  else
    if( orderOfAccuracy.eq.2 )then
      call bcOptWave3dOrder2( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    elseif( orderOfAccuracy.eq.4 )then
      call bcOptWave3dOrder4( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,un,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    else
      stop 7777
    end if    

  end if


  return
  end

