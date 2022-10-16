! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=3, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
!  MASK : USEMASK or NOMASK 
!  FORCING : NOFORCING, TZ, USEFORCING 
! ===========================================================================
#beginMacro update3dOrder6Rectangular(DIM,ORDER,ORDERINTIME,GRIDTYPE,MASK,FORCING)
#If #FORCING eq noForcing
#defineMacro FV(m) 
#Else
#defineMacro FV(m) +dtSq*fv(m)
#End
cxx=1./dx(0)**2;
cyy=1./dx(1)**2;
czz=1./dx(2)**2;
numGhost1=2;
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
numGhost1=1;
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 
    d400(i1,i2,i3,0) = d200(i1+1,i2,i3,0) - 2.*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0)
    d220(i1,i2,i3,0) = d020(i1+1,i2,i3,0) - 2.*d020(i1,i2,i3,0) + d020(i1-1,i2,i3,0)
    d040(i1,i2,i3,0) = d020(i1,i2+1,i3,0) - 2.*d020(i1,i2,i3,0) + d020(i1,i2-1,i3,0)
    d202(i1,i2,i3,0) = d002(i1+1,i2,i3,0) - 2.*d002(i1,i2,i3,0) + d002(i1-1,i2,i3,0)
    d022(i1,i2,i3,0) = d002(i1,i2+1,i3,0) - 2.*d002(i1,i2,i3,0) + d002(i1,i2-1,i3,0)
    d004(i1,i2,i3,0) = d002(i1,i2,i3+1,0) - 2.*d002(i1,i2,i3,0) + d002(i1,i2,i3-1,0)
  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 

! ------------ MAIN LOOP, 3D, ORDER 6 ----------------
numGhost1=0; 
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);

! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
! Example: 
! u.xx = D+D-( I + cxx1*D+D- + cxx2*(D+D-x)^2 + ...
cx0 = 1./(dx(0)**1); cx1 = (-1/6.); cx2 = (1/30.); 
cy0 = 1./(dx(1)**1); cy1 = (-1/6.); cy2 = (1/30.); 
cz0 = 1./(dx(2)**1); cz1 = (-1/6.); cz2 = (1/30.); 
cxx0 = 1./(dx(0)**2); cxx1 = (-1/12.); cxx2 = (1/90.); 
cyy0 = 1./(dx(1)**2); cyy1 = (-1/12.); cyy2 = (1/90.); 
czz0 = 1./(dx(2)**2); czz1 = (-1/12.); czz2 = (1/90.); 
cxxx0 = 1./(dx(0)**3); cxxx1 = (-1/4.); cxxx2 = (7/120.); 
cyyy0 = 1./(dx(1)**3); cyyy1 = (-1/4.); cyyy2 = (7/120.); 
czzz0 = 1./(dx(2)**3); czzz1 = (-1/4.); czzz2 = (7/120.); 
cxxxx0 = 1./(dx(0)**4); cxxxx1 = (-1/6.); cxxxx2 = (7/240.); 
cyyyy0 = 1./(dx(1)**4); cyyyy1 = (-1/6.); cyyyy2 = (7/240.); 
czzzz0 = 1./(dx(2)**4); czzzz1 = (-1/6.); czzzz2 = (7/240.); 
cxxxxx0 = 1./(dx(0)**5); cxxxxx1 = (-1/3.); cxxxxx2 = (13/144.); 
cyyyyy0 = 1./(dx(1)**5); cyyyyy1 = (-1/3.); cyyyyy2 = (13/144.); 
czzzzz0 = 1./(dx(2)**5); czzzzz1 = (-1/3.); czzzzz2 = (13/144.); 
cxxxxxx0 = 1./(dx(0)**6); cxxxxxx1 = (-1/4.); cxxxxxx2 = (13/240.); 
cyyyyyy0 = 1./(dx(1)**6); cyyyyyy1 = (-1/4.); cyyyyyy2 = (13/240.); 
czzzzzz0 = 1./(dx(2)**6); czzzzzz1 = (-1/4.); czzzzzz2 = (13/240.); 

fv(m)=0.
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 

    ! -- Note: Highest differences do not need to be stored in an array 

    ! ----- DIFFERENCES OF U -----
    d200i = d200(i1,i2,i3,0)
    d400i = d400(i1,i2,i3,0)
    d600i = d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0)
    d020i = d020(i1,i2,i3,0)
    d220i = d220(i1,i2,i3,0)
    d420i = d220(i1+1,i2,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0)
    d040i = d040(i1,i2,i3,0)
    d240i = d040(i1+1,i2,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1-1,i2,i3,0)
    d060i = d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0)
    d002i = d002(i1,i2,i3,0)
    d202i = d202(i1,i2,i3,0)
    d402i = d202(i1+1,i2,i3,0) - 2.*d202(i1,i2,i3,0) + d202(i1-1,i2,i3,0)
    d022i = d022(i1,i2,i3,0)
    d222i = d022(i1+1,i2,i3,0) - 2.*d022(i1,i2,i3,0) + d022(i1-1,i2,i3,0)
    d042i = d022(i1,i2+1,i3,0) - 2.*d022(i1,i2,i3,0) + d022(i1,i2-1,i3,0)
    d004i = d004(i1,i2,i3,0)
    d204i = d004(i1+1,i2,i3,0) - 2.*d004(i1,i2,i3,0) + d004(i1-1,i2,i3,0)
    d024i = d004(i1,i2+1,i3,0) - 2.*d004(i1,i2,i3,0) + d004(i1,i2-1,i3,0)
    d006i = d004(i1,i2,i3+1,0) - 2.*d004(i1,i2,i3,0) + d004(i1,i2,i3-1,0)



    ! ----- spatial derivatives of order 2 ----
    uxx =cxx0*(d200i + cxx1*d400i + cxx2*d600i)
    uyy =cyy0*(d020i + cyy1*d040i + cyy2*d060i)
    uzz =czz0*(d002i + czz1*d004i + czz2*d006i)

    ! ----- spatial derivatives of order 4 ----
    uxxxx =cxxxx0*(d400i + cxxxx1*d600i)
    uxxyy =cxx0*cyy0*(d220i + cxx1*d420i + cyy1*d240i)
    uyyyy =cyyyy0*(d040i + cyyyy1*d060i)
    uxxzz =cxx0*czz0*(d202i + cxx1*d402i + czz1*d204i)
    uyyzz =cyy0*czz0*(d022i + cyy1*d042i + czz1*d024i)
    uzzzz =czzzz0*(d004i + czzzz1*d006i)

    ! ----- spatial derivatives of order 6 ----
    uxxxxxx =cxxxxxx0*(d600i)
    uxxxxyy =cxxxx0*cyy0*(d420i)
    uxxyyyy =cxx0*cyyyy0*(d240i)
    uyyyyyy =cyyyyyy0*(d060i)
    uxxxxzz =cxxxx0*czz0*(d402i)
    uxxyyzz =cxx0*cyy0*czz0*(d222i)
    uyyyyzz =cyyyy0*czz0*(d042i)
    uxxzzzz =cxx0*czzzz0*(d204i)
    uyyzzzz =cyy0*czzzz0*(d024i)
    uzzzzzz =czzzzzz0*(d006i)

  #If #FORCING ne "NOFORCING" 
    getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) 
  #End 

    ! ---- Modified Equation UPDATE ----- 
    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) \
                     + cdtsq*( uxx + uyy +uzz ) \
                     + cdtPow4By12*( uxxxx + uyyyy + uzzzz + 2.*(uxxyy+uxxzz+uyyzz) )  \
                     + cdtPow6By360*( uxxxxxx + uyyyyyy + uzzzzzz \
                                + 3.*( uxxxxyy + uxxyyyy + uxxxxzz + uxxzzzz + uyyyyzz + uyyzzzz ) \
                                + 6.*(  uxxyyzz  ) ) \
                       FV(m)  

  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 
#endMacro