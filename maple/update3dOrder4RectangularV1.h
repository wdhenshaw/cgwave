! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=3, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
!  MASK : USEMASK or NOMASK 
!  FORCING : NOFORCING, TZ, USEFORCING 
! ===========================================================================
#beginMacro update3dOrder4Rectangular(DIM,ORDER,ORDERINTIME,GRIDTYPE,MASK,FORCING)
#If #FORCING eq noForcing
#defineMacro FV(m) 
#Else
#defineMacro FV(m) +dtSq*fv(m)
#End
cxx=1./dx(0)**2;
cyy=1./dx(1)**2;
czz=1./dx(2)**2;
numGhost1=1;
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 
    d200(i1,i2,i3,0) = u(i1+1,i2,i3,0) - 2.*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
    d020(i1,i2,i3,0) = u(i1,i2+1,i3,0) - 2.*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
    d002(i1,i2,i3,0) = u(i1,i2,i3+1,0) - 2.*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 

! ------------ MAIN LOOP, 3D, ORDER 4 ----------------
numGhost1=0; 
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);

! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
! Example: 
! u.xx = D+D-( I + cxx1*D+D- + cxx2*(D+D-x)^2 + ...
cx0 = 1./(dx(0)**1); cx1 = (-1/6.); 
cy0 = 1./(dx(1)**1); cy1 = (-1/6.); 
cz0 = 1./(dx(2)**1); cz1 = (-1/6.); 
cxx0 = 1./(dx(0)**2); cxx1 = (-1/12.); 
cyy0 = 1./(dx(1)**2); cyy1 = (-1/12.); 
czz0 = 1./(dx(2)**2); czz1 = (-1/12.); 
cxxx0 = 1./(dx(0)**3); cxxx1 = (-1/4.); 
cyyy0 = 1./(dx(1)**3); cyyy1 = (-1/4.); 
czzz0 = 1./(dx(2)**3); czzz1 = (-1/4.); 
cxxxx0 = 1./(dx(0)**4); cxxxx1 = (-1/6.); 
cyyyy0 = 1./(dx(1)**4); cyyyy1 = (-1/6.); 
czzzz0 = 1./(dx(2)**4); czzzz1 = (-1/6.); 

fv(m)=0.
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 

    ! -- Note: Highest differences do not need to be stored in an array 

    ! ----- DIFFERENCES OF U -----
    d200i = d200(i1,i2,i3,0)
    d400i = d200(i1+1,i2,i3,0) - 2.*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0)
    d020i = d020(i1,i2,i3,0)
    d220i = d020(i1+1,i2,i3,0) - 2.*d020(i1,i2,i3,0) + d020(i1-1,i2,i3,0)
    d040i = d020(i1,i2+1,i3,0) - 2.*d020(i1,i2,i3,0) + d020(i1,i2-1,i3,0)
    d002i = d002(i1,i2,i3,0)
    d202i = d002(i1+1,i2,i3,0) - 2.*d002(i1,i2,i3,0) + d002(i1-1,i2,i3,0)
    d022i = d002(i1,i2+1,i3,0) - 2.*d002(i1,i2,i3,0) + d002(i1,i2-1,i3,0)
    d004i = d002(i1,i2,i3+1,0) - 2.*d002(i1,i2,i3,0) + d002(i1,i2,i3-1,0)



    ! ----- spatial derivatives of order 2 ----
    uxx =cxx0*(d200i + cxx1*d400i)
    uyy =cyy0*(d020i + cyy1*d040i)
    uzz =czz0*(d002i + czz1*d004i)

    ! ----- spatial derivatives of order 4 ----
    uxxxx =cxxxx0*(d400i)
    uxxyy =cxx0*cyy0*(d220i)
    uyyyy =cyyyy0*(d040i)
    uxxzz =cxx0*czz0*(d202i)
    uyyzz =cyy0*czz0*(d022i)
    uzzzz =czzzz0*(d004i)

  #If #FORCING ne "NOFORCING" 
    getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) 
  #End 

    ! ---- Modified Equation UPDATE ----- 
    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) \
                     + cdtsq*( uxx + uyy +uzz ) \
                     + cdtPow4By12*( uxxxx + uyyyy + uzzzz + 2.*(uxxyy+uxxzz+uyyzz) )  \
                       FV(m)  

  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 
#endMacro