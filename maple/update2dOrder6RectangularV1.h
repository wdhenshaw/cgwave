! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=2, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
!  MASK : USEMASK or NOMASK 
!  FORCING : NOFORCING, TZ, USEFORCING 
! ===========================================================================
#beginMacro update2dOrder6Rectangular(DIM,ORDER,ORDERINTIME,GRIDTYPE,MASK,FORCING)
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
  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 

! ------------ MAIN LOOP, 2D, ORDER 6 ----------------
numGhost1=0; 
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);

! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
! Example: 
! u.xx = D+D-( I + cxx1*D+D- + cxx2*(D+D-x)^2 + ...
cx0 = 1./(dx(0)**1); cx1 = (-1/6.); cx2 = (1/30.); 
cy0 = 1./(dx(1)**1); cy1 = (-1/6.); cy2 = (1/30.); 
cxx0 = 1./(dx(0)**2); cxx1 = (-1/12.); cxx2 = (1/90.); 
cyy0 = 1./(dx(1)**2); cyy1 = (-1/12.); cyy2 = (1/90.); 
cxxx0 = 1./(dx(0)**3); cxxx1 = (-1/4.); cxxx2 = (7/120.); 
cyyy0 = 1./(dx(1)**3); cyyy1 = (-1/4.); cyyy2 = (7/120.); 
cxxxx0 = 1./(dx(0)**4); cxxxx1 = (-1/6.); cxxxx2 = (7/240.); 
cyyyy0 = 1./(dx(1)**4); cyyyy1 = (-1/6.); cyyyy2 = (7/240.); 
cxxxxx0 = 1./(dx(0)**5); cxxxxx1 = (-1/3.); cxxxxx2 = (13/144.); 
cyyyyy0 = 1./(dx(1)**5); cyyyyy1 = (-1/3.); cyyyyy2 = (13/144.); 
cxxxxxx0 = 1./(dx(0)**6); cxxxxxx1 = (-1/4.); cxxxxxx2 = (13/240.); 
cyyyyyy0 = 1./(dx(1)**6); cyyyyyy1 = (-1/4.); cyyyyyy2 = (13/240.); 

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



    ! ----- spatial derivatives of order 2 ----
    uxx =cxx0*(d200i + cxx1*d400i + cxx2*d600i)
    uyy =cyy0*(d020i + cyy1*d040i + cyy2*d060i)

    ! ----- spatial derivatives of order 4 ----
    uxxxx =cxxxx0*(d400i + cxxxx1*d600i)
    uxxyy =cxx0*cyy0*(d220i + cxx1*d420i + cyy1*d240i)
    uyyyy =cyyyy0*(d040i + cyyyy1*d060i)

    ! ----- spatial derivatives of order 6 ----
    uxxxxxx =cxxxxxx0*(d600i)
    uxxxxyy =cxxxx0*cyy0*(d420i)
    uxxyyyy =cxx0*cyyyy0*(d240i)
    uyyyyyy =cyyyyyy0*(d060i)

  #If #FORCING ne "NOFORCING" 
    getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) 
  #End 

    ! ---- Modified Equation UPDATE ----- 
    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) \
                     + cdtsq*( uxx + uyy ) \
                     + cdtPow4By12*( uxxxx + uyyyy + 2.*uxxyy )  \
                     + cdtPow6By360*( uxxxxxx + uyyyyyy\
                               + 3.*( uxxxxyy + uxxyyyy) ) \
     FV(m)  

  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 
#endMacro