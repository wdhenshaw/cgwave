! ===========================================================================
!   Modified Equation : order=2, DIMENSIONS=3, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
!  MASK : USEMASK or NOMASK 
!  FORCING : NOFORCING, TZ, USEFORCING 
! ===========================================================================
#beginMacro update3dOrder2Rectangular(DIM,ORDER,ORDERINTIME,GRIDTYPE,MASK,FORCING)
#If #FORCING eq noForcing
#defineMacro FV(m) 
#Else
#defineMacro FV(m) +dtSq*fv(m)
#End
cxx=1./dx(0)**2;
cyy=1./dx(1)**2;
czz=1./dx(2)**2;

! ------------ MAIN LOOP, 3D, ORDER 2 ----------------
numGhost1=0; 
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);

! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
! Example: 
! u.xx = D+D-( I + cxx1*D+D- + cxx2*(D+D-x)^2 + ...
cx0 = 1./(dx(0)**1); 
cy0 = 1./(dx(1)**1); 
cz0 = 1./(dx(2)**1); 
cxx0 = 1./(dx(0)**2); 
cyy0 = 1./(dx(1)**2); 
czz0 = 1./(dx(2)**2); 

fv(m)=0.
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 

    ! -- Note: Highest differences do not need to be stored in an array 

    ! ----- DIFFERENCES OF U -----
    d200i = d000(i1+1,i2,i3,0) - 2.*d000(i1,i2,i3,0) + d000(i1-1,i2,i3,0)
    d020i = d000(i1,i2+1,i3,0) - 2.*d000(i1,i2,i3,0) + d000(i1,i2-1,i3,0)
    d002i = d000(i1,i2,i3+1,0) - 2.*d000(i1,i2,i3,0) + d000(i1,i2,i3-1,0)



    ! ----- spatial derivatives of order 2 ----
    uxx =cxx0*(d200i)
    uyy =cyy0*(d020i)
    uzz =czz0*(d002i)

  #If #FORCING ne "NOFORCING" 
    getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) 
  #End 

    ! ---- Modified Equation UPDATE ----- 
    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) \
                     + cdtsq*( uxx + uyy +uzz ) \
                       FV(m)  

  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 
#endMacro