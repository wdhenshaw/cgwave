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

! From bcOptSmFOS.bf
! DataBase *pdb = &parameters.dbase;
! double precision pdb  ! pointer to data base
! ====================================================================
! Look up an integer parameter from the data base
! ====================================================================

! ====================================================================
! Look up a real parameter from the data base
! ====================================================================

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
!   OPTION = EXPLIICT, IMPLICIT
!      EXPLICIT: compute v(i1,i2,i3,0) = uLap
!      IMPLICIT: compute v(i1,i2,i3,0) = uLap, and v(i1,i2,i3,1) = umLap
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
!
! Macro: FILL IN THE RHS FOR IMPLICIT TIME-STEPPING
!          RECTANGULAR GRID + SUPERGRID 
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
!    *** OLD ****
!
!   DIM : dimension : 2 or 3
!   ORDER : 2 or 4   
!   GRIDTYPE : rectangular or curvilinear
! ========================================================================================
! #beginMacro addUpwindDissImplicitOld(DIM,ORDER,GRIDTYPE)

!   if( (.false. .or. debug.gt.3) .and. t.lt.4*dt )then
!     write(*,'(/,"advWave:addUpwindDissImplicitOld: UPWIND DISS dim=DIM order=ORDER grid=GRIDTYPE... t=",e10.2)') t
!     if( gridType==rectangular )then
!       write(*,'(" upwind-diss-coeff: adxSosup=",3(1pe12.4,1x))') adxSosup(0), adxSosup(1),adxSosup(2)
!     else
!       write(*,'(" upwind-diss-coeff: adSosup=",1pe12.4)') adSosup
!     end if
!   end if

!   ! from formImplicitTimeSteppingMatrix
!     ! // fourth-order dissipation for 2nd-order scheme:
!     ! Real upwindCoeff4[4][5] = { 1.,-4.,6.,-4.,1.,
!     !                             1.,-3.,3.,-1.,0.,   // extrap right-most point D-^3 u(2)
!     !                             0.,-1.,3.,-3.,1.,   // extrap left -most point D+^3 u(-2)
!     !                             0.,-1.,2.,-1.,0.
!     !                           };

!     ! // sixth-order dissipation for 4th-order scheme
!     ! Real upwindCoeff6[4][7] = {1.,-6.,15.,-20.,15.,-6.,1.,
!     !                            1.,-5.,10.,-10., 5.,-1.,0.,  // extrap right-most point D-^5 u(3)
!     !                            0.,-1., 5.,-10.,10., 5.,1.,  // extrap left -most point D+^5 u(-3)
!     !                            0.,-1., 4., -6., 4.,-1.,0.
!     !                           };
  

!   ! -- Note: Could adjust loop bounds to avoid Dirichlet boundaries

!   m=0 ! component number 
!   ec = 0 ! component number
  
!   ! Compute some things needed in the loops below
!   if( adjustHelmholtzForUpwinding.eq.1 )then
!     do freq=0,numberOfFrequencies-1
!       cosineFactor(freq) = cos(frequencyArray(freq)*(t+dt)) - cos(frequencyArray(freq)*(t-dt))
!       sineFactor(freq)   = sin(frequencyArray(freq)*(t+dt)) - sin(frequencyArray(freq)*(t-dt))
!     end do
!   end if 
    
!   ! stencil width = order + 1
!   ! upwind stencil = stencilWidth + 2 = order+1 + 2
!   upwindHalfStencilWidth = (orderOfAccuracy+2)/2 

!   beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)

!     #If #GRIDTYPE eq "curvilinear"
!       ! --- Compute UPW coefficient for curvilinear grids ---
!       #If #DIM eq "2"
!         getSosupDissipationCoeff2d(adxSosup)
!       #Else
!         getSosupDissipationCoeff3d(adxSosup)
!       #End
!     #End

!     ! --- loop over directions ---
!     do dir=0,nd-1

!       idv(0)=0; idv(1)=0; idv(2)=0
!       idv(dir)=1 ! active direction

!       ! check if left-most and right-most entries in the upwind stencil are valid 
!       i1l = i1-upwindHalfStencilWidth*idv(0); i1r = i1+upwindHalfStencilWidth*idv(0);
!       i2l = i2-upwindHalfStencilWidth*idv(1); i2r = i2+upwindHalfStencilWidth*idv(1);
!       i3l = i3-upwindHalfStencilWidth*idv(2); i3r = i3+upwindHalfStencilWidth*idv(2);

!       !  Note: there are at most four cases at any order, since we have order/2 layers of interpolation points 
!       !   Example, order=2, upwind-order=4
!       !     X---X---C---X---X           C = center point = valid discretization point
!       !     X---X---C---X               missing right-most 
!       !         X---C---X---X           missing left-most
!       !         X---C---X               missing left and right-most

!       upwCase=0 
!       if( mask(i1l,i2l,i3l) .ne. 0 .and. mask(i1r,i2r,i3r) .ne. 0 ) then
!         upwCase=0  ! centred, full-width stencil
!       else if( mask(i1l,i2l,i3l) .ne. 0 ) then
!         upwCase=1  ! left biased stencil
!       else if( mask(i1r,i2r,i3r) .ne. 0 ) then
!         upwCase=2  ! right biased stencil   
!       else  
!         upwCase=3  ! centred smaller stencil     
!       end if



!       upw = 0. 
!       do iStencil=-upwindHalfStencilWidth,upwindHalfStencilWidth

!         j1 = i1 + iStencil*idv(0);  j2 = i2 + iStencil*idv(1);  j3 = i3 + iStencil*idv(2)

!         umj = um(j1,j2,j3,0)

!         if( adjustHelmholtzForUpwinding.eq.1 )then
!           do freq=0,numberOfFrequencies-1
!             umj = umj + cosineFactor(freq)*vh(j1,j2,j3,freq)
!           end do
!         end if

!         upw = upw + upwindCoeff(iStencil,upwCase)*umj

!         ! upw = upw + upwindCoeff(iStencil,upwCase)*um(j1,j2,j3,ec)  ! *** CHECK ME 


!         ! write(*,'("upw-rhs: i1,i2=",2i4," j1,j2=",2i4," upwindCoeff=",1pe9.2, " um=",1pe9.2," upw=",1pe9.2)') i1,i2,j1,j2,upwindCoeff(iStencil,upwCase),um(j1,j2,j3,ec),upw
!       end do 
!       ! if( abs(upw).gt.1e-10 )then
!       !   write(*,'(">>upw-rhs: i1,i2=",2i4," upw=",1pe9.2)') i1,i2,upw
!       ! end if 

!       ! This is the coeff of um in
!       !         + adxSosup(dir)*(UpwindStencil)( un - um )
!       ! un(i1,i2,i3,ec) = un(i1,i2,i3,ec) - adxSosup(dir)*upw 
!       ! *wdh* July 3, 2023 CHANGE SIGN
!       un(i1,i2,i3,ec) = un(i1,i2,i3,ec) + adxSosup(dir)*upw 

!       ! TEST un(i1,i2,i3,ec) = un(i1,i2,i3,ec) - adxSosup(dir)*upw 

!     end do ! end do fir 

!   endLoopsMask()



! #endMacro

  
! =========================================================================================
! Macro: ADD UPWIND DISSIPATION FOR IMPLICIT TIME-STEPPING
!
!       ** NEW VERSION : July 10, 2023 **
!
!    Assumes that interpolation neighbours ahave been filled in the the implicitSolve 
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
  ! buildFile(advWave2dOrder6r,2,6,rectangular)
  ! buildFile(advWave3dOrder6r,3,6,rectangular)
  ! buildFile(advWave2dOrder6c,2,6,curvilinear)
  ! buildFile(advWave3dOrder6c,3,6,curvilinear)

  ! ! ORDER 8 -- needed for upwinding,  SuperGrid
  ! buildFile(advWave2dOrder8r,2,8,rectangular)
  ! buildFile(advWave3dOrder8r,3,8,rectangular)
  ! buildFile(advWave2dOrder8c,2,8,curvilinear)
  ! buildFile(advWave3dOrder8c,3,8,curvilinear)


      subroutine advWave( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
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
  double precision pdb  ! pointer to the data base
  integer ipar(0:*)
  real rpar(0:*)


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

       
  option                    = ipar( 0)
  gridType                  = ipar( 2)
  orderOfAccuracy           = ipar( 3)
  orderInTime               = ipar( 4)
  modifiedEquationApproach  = ipar(18)

  ! write(*,*) 'Inside advWave...'
  ! write(*,*) 'f: '
  ! write(*,*) (((f(i1,i2,i3,0),i1=nd1a,nd1b),i2=nd2a,nd2b),i3=nd3a,nd3b)

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
        call advWave2dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWave2dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWave3dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWave3dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else
       stop 8843
      end if
    else if( modifiedEquationApproach.eq.hierarchicalME ) then
      ! new Hierarchical scheme
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveME2dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveME2dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWaveME3dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWaveME3dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else
       stop 8843
      end if 
    else if( modifiedEquationApproach.eq.stencilME ) then

     if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveStencil2dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveStencil2dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
       call advWaveStencil3dOrder2r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       call advWaveStencil3dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
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
        call advWave2dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWave2dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWave3dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWave3dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else
       stop 8843
      end if
    else if( modifiedEquationApproach.eq.hierarchicalME ) then
      ! new Hierarchical scheme
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveME2dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveME2dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWaveME3dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWaveME3dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else
       stop 8843
      end if 

    else if( modifiedEquationApproach.eq.stencilME ) then

     if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveStencil2dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveStencil2dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWaveStencil3dOrder4r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       call advWaveStencil3dOrder4c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
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
        call advWave2dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWave2dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWave3dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWave3dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else
       stop 8843
      end if
    else if( modifiedEquationApproach.eq.hierarchicalME ) then
      ! new Hierarchical scheme
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveME2dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveME2dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWaveME3dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWaveME3dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else
       stop 8843
      end if 
    else if( modifiedEquationApproach.eq.stencilME ) then

     if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveStencil2dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveStencil2dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
       call advWaveStencil3dOrder6r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       call advWaveStencil3dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
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
        call advWave2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWave2dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWave3dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWave3dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else
       stop 8843
      end if
    else if( modifiedEquationApproach.eq.hierarchicalME ) then
      ! new Hierarchical scheme
      if( nd.eq.2 .and. gridType.eq.rectangular )then
        call advWaveME2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
        call advWaveME2dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
        call advWaveME3dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
        call advWaveME3dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else
       stop 8843
      end if 

    else if( modifiedEquationApproach.eq.stencilME ) then


     if( nd.eq.2 .and. gridType.eq.rectangular )then
       call advWaveStencil2dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(nd.eq.2 .and. gridType.eq.curvilinear )then
       call advWaveStencil2dOrder8c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.rectangular )then
       call advWaveStencil3dOrder8r( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,stencilCoeff,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,pdb,ipar,rpar,ierr )
      else if(  nd.eq.3 .and. gridType.eq.curvilinear )then
       ! call advWaveStencil3dOrder8c( ARGLIST() )
      else
       stop 8843
      end if                  
    else  
      write(*,'("Unknown modifiedEquationApproach=",i6)') modifiedEquationApproach
      stop 1111
    end if

  else
    write(*,'(" advWave:ERROR: un-implemented order of accuracy =",i6)') orderOfAccuracy
      ! '
    stop 11122
  end if

  return
  end








