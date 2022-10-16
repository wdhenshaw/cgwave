! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=3, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
!  MASK : USEMASK or NOMASK 
!  FORCING : NOFORCING, TZ, USEFORCING 
! ===========================================================================
#beginMacro update3dOrder8Rectangular(DIM,ORDER,ORDERINTIME,GRIDTYPE,MASK,FORCING)
#If #FORCING eq noForcing
#defineMacro FV(m) 
#Else
#defineMacro FV(m) +dtSq*fv(m)
#End
cxx=1./dx(0)**2;
cyy=1./dx(1)**2;
czz=1./dx(2)**2;
numGhost1=3;
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
numGhost1=2;
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
numGhost1=1;
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 
    d600(i1,i2,i3,0) = d400(i1+1,i2,i3,0) - 2.*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0)
    d420(i1,i2,i3,0) = d220(i1+1,i2,i3,0) - 2.*d220(i1,i2,i3,0) + d220(i1-1,i2,i3,0)
    d240(i1,i2,i3,0) = d040(i1+1,i2,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1-1,i2,i3,0)
    d060(i1,i2,i3,0) = d040(i1,i2+1,i3,0) - 2.*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0)
    d402(i1,i2,i3,0) = d202(i1+1,i2,i3,0) - 2.*d202(i1,i2,i3,0) + d202(i1-1,i2,i3,0)
    d222(i1,i2,i3,0) = d022(i1+1,i2,i3,0) - 2.*d022(i1,i2,i3,0) + d022(i1-1,i2,i3,0)
    d042(i1,i2,i3,0) = d022(i1,i2+1,i3,0) - 2.*d022(i1,i2,i3,0) + d022(i1,i2-1,i3,0)
    d204(i1,i2,i3,0) = d004(i1+1,i2,i3,0) - 2.*d004(i1,i2,i3,0) + d004(i1-1,i2,i3,0)
    d024(i1,i2,i3,0) = d004(i1,i2+1,i3,0) - 2.*d004(i1,i2,i3,0) + d004(i1,i2-1,i3,0)
    d006(i1,i2,i3,0) = d004(i1,i2,i3+1,0) - 2.*d004(i1,i2,i3,0) + d004(i1,i2,i3-1,0)
  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 

! ------------ MAIN LOOP, 3D, ORDER 8 ----------------
numGhost1=0; 
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);

! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
! Example: 
! u.xx = D+D-( I + cxx1*D+D- + cxx2*(D+D-x)^2 + ...
cx0 = 1./(dx(0)**1); cx1 = (-1/6.); cx2 = (1/30.); cx3 = (-1/140.); 
cy0 = 1./(dx(1)**1); cy1 = (-1/6.); cy2 = (1/30.); cy3 = (-1/140.); 
cz0 = 1./(dx(2)**1); cz1 = (-1/6.); cz2 = (1/30.); cz3 = (-1/140.); 
cxx0 = 1./(dx(0)**2); cxx1 = (-1/12.); cxx2 = (1/90.); cxx3 = (-1/560.); 
cyy0 = 1./(dx(1)**2); cyy1 = (-1/12.); cyy2 = (1/90.); cyy3 = (-1/560.); 
czz0 = 1./(dx(2)**2); czz1 = (-1/12.); czz2 = (1/90.); czz3 = (-1/560.); 
cxxx0 = 1./(dx(0)**3); cxxx1 = (-1/4.); cxxx2 = (7/120.); cxxx3 = (-41/3024.); 
cyyy0 = 1./(dx(1)**3); cyyy1 = (-1/4.); cyyy2 = (7/120.); cyyy3 = (-41/3024.); 
czzz0 = 1./(dx(2)**3); czzz1 = (-1/4.); czzz2 = (7/120.); czzz3 = (-41/3024.); 
cxxxx0 = 1./(dx(0)**4); cxxxx1 = (-1/6.); cxxxx2 = (7/240.); cxxxx3 = (-41/7560.); 
cyyyy0 = 1./(dx(1)**4); cyyyy1 = (-1/6.); cyyyy2 = (7/240.); cyyyy3 = (-41/7560.); 
czzzz0 = 1./(dx(2)**4); czzzz1 = (-1/6.); czzzz2 = (7/240.); czzzz3 = (-41/7560.); 
cxxxxx0 = 1./(dx(0)**5); cxxxxx1 = (-1/3.); cxxxxx2 = (13/144.); cxxxxx3 = (-139/6048.); 
cyyyyy0 = 1./(dx(1)**5); cyyyyy1 = (-1/3.); cyyyyy2 = (13/144.); cyyyyy3 = (-139/6048.); 
czzzzz0 = 1./(dx(2)**5); czzzzz1 = (-1/3.); czzzzz2 = (13/144.); czzzzz3 = (-139/6048.); 
cxxxxxx0 = 1./(dx(0)**6); cxxxxxx1 = (-1/4.); cxxxxxx2 = (13/240.); cxxxxxx3 = (-139/12096.); 
cyyyyyy0 = 1./(dx(1)**6); cyyyyyy1 = (-1/4.); cyyyyyy2 = (13/240.); cyyyyyy3 = (-139/12096.); 
czzzzzz0 = 1./(dx(2)**6); czzzzzz1 = (-1/4.); czzzzzz2 = (13/240.); czzzzzz3 = (-139/12096.); 
cxxxxxxx0 = 1./(dx(0)**7); cxxxxxxx1 = (0.); cxxxxxxx2 = (0.); cxxxxxxx3 = (0.); 
cyyyyyyy0 = 1./(dx(1)**7); cyyyyyyy1 = (0.); cyyyyyyy2 = (0.); cyyyyyyy3 = (0.); 
czzzzzzz0 = 1./(dx(2)**7); czzzzzzz1 = (0.); czzzzzzz2 = (0.); czzzzzzz3 = (0.); 
cxxxxxxxx0 = 1./(dx(0)**8); cxxxxxxxx1 = (0.); cxxxxxxxx2 = (0.); cxxxxxxxx3 = (0.); 
cyyyyyyyy0 = 1./(dx(1)**8); cyyyyyyyy1 = (0.); cyyyyyyyy2 = (0.); cyyyyyyyy3 = (0.); 
czzzzzzzz0 = 1./(dx(2)**8); czzzzzzzz1 = (0.); czzzzzzzz2 = (0.); czzzzzzzz3 = (0.); 

fv(m)=0.
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 

    ! -- Note: Highest differences do not need to be stored in an array 

    ! ----- DIFFERENCES OF U -----
    d200i = d200(i1,i2,i3,0)
    d400i = d400(i1,i2,i3,0)
    d600i = d600(i1,i2,i3,0)
    d800i = d600(i1+1,i2,i3,0) - 2.*d600(i1,i2,i3,0) + d600(i1-1,i2,i3,0)
    d020i = d020(i1,i2,i3,0)
    d220i = d220(i1,i2,i3,0)
    d420i = d420(i1,i2,i3,0)
    d620i = d420(i1+1,i2,i3,0) - 2.*d420(i1,i2,i3,0) + d420(i1-1,i2,i3,0)
    d040i = d040(i1,i2,i3,0)
    d240i = d240(i1,i2,i3,0)
    d440i = d240(i1+1,i2,i3,0) - 2.*d240(i1,i2,i3,0) + d240(i1-1,i2,i3,0)
    d060i = d060(i1,i2,i3,0)
    d260i = d060(i1+1,i2,i3,0) - 2.*d060(i1,i2,i3,0) + d060(i1-1,i2,i3,0)
    d080i = d060(i1,i2+1,i3,0) - 2.*d060(i1,i2,i3,0) + d060(i1,i2-1,i3,0)
    d002i = d002(i1,i2,i3,0)
    d202i = d202(i1,i2,i3,0)
    d402i = d402(i1,i2,i3,0)
    d602i = d402(i1+1,i2,i3,0) - 2.*d402(i1,i2,i3,0) + d402(i1-1,i2,i3,0)
    d022i = d022(i1,i2,i3,0)
    d222i = d222(i1,i2,i3,0)
    d422i = d222(i1+1,i2,i3,0) - 2.*d222(i1,i2,i3,0) + d222(i1-1,i2,i3,0)
    d042i = d042(i1,i2,i3,0)
    d242i = d042(i1+1,i2,i3,0) - 2.*d042(i1,i2,i3,0) + d042(i1-1,i2,i3,0)
    d062i = d042(i1,i2+1,i3,0) - 2.*d042(i1,i2,i3,0) + d042(i1,i2-1,i3,0)
    d004i = d004(i1,i2,i3,0)
    d204i = d204(i1,i2,i3,0)
    d404i = d204(i1+1,i2,i3,0) - 2.*d204(i1,i2,i3,0) + d204(i1-1,i2,i3,0)
    d024i = d024(i1,i2,i3,0)
    d224i = d024(i1+1,i2,i3,0) - 2.*d024(i1,i2,i3,0) + d024(i1-1,i2,i3,0)
    d044i = d024(i1,i2+1,i3,0) - 2.*d024(i1,i2,i3,0) + d024(i1,i2-1,i3,0)
    d006i = d006(i1,i2,i3,0)
    d206i = d006(i1+1,i2,i3,0) - 2.*d006(i1,i2,i3,0) + d006(i1-1,i2,i3,0)
    d026i = d006(i1,i2+1,i3,0) - 2.*d006(i1,i2,i3,0) + d006(i1,i2-1,i3,0)
    d008i = d006(i1,i2,i3+1,0) - 2.*d006(i1,i2,i3,0) + d006(i1,i2,i3-1,0)



    ! ----- spatial derivatives of order 2 ----
    uxx =cxx0*(d200i + cxx1*d400i + cxx2*d600i + cxx3*d800i)
    uyy =cyy0*(d020i + cyy1*d040i + cyy2*d060i + cyy3*d080i)
    uzz =czz0*(d002i + czz1*d004i + czz2*d006i + czz3*d008i)

    ! ----- spatial derivatives of order 4 ----
    uxxxx =cxxxx0*(d400i + cxxxx1*d600i + cxxxx2*d800i)
    uxxyy =cxx0*cyy0*(d220i + cxx1*d420i + cxx2*d620i + cyy1*d240i + cxx1*cyy1*d440i + cyy2*d260i)
    uyyyy =cyyyy0*(d040i + cyyyy1*d060i + cyyyy2*d080i)
    uxxzz =cxx0*czz0*(d202i + cxx1*d402i + cxx2*d602i + czz1*d204i + cxx1*czz1*d404i + czz2*d206i)
    uyyzz =cyy0*czz0*(d022i + cyy1*d042i + cyy2*d062i + czz1*d024i + cyy1*czz1*d044i + czz2*d026i)
    uzzzz =czzzz0*(d004i + czzzz1*d006i + czzzz2*d008i)

    ! ----- spatial derivatives of order 6 ----
    uxxxxxx =cxxxxxx0*(d600i + cxxxxxx1*d800i)
    uxxxxyy =cxxxx0*cyy0*(d420i + cxxxx1*d620i + cyy1*d440i)
    uxxyyyy =cxx0*cyyyy0*(d240i + cxx1*d440i + cyyyy1*d260i)
    uyyyyyy =cyyyyyy0*(d060i + cyyyyyy1*d080i)
    uxxxxzz =cxxxx0*czz0*(d402i + cxxxx1*d602i + czz1*d404i)
    uxxyyzz =cxx0*cyy0*czz0*(d222i + cxx1*d422i + cyy1*d242i + czz1*d224i)
    uyyyyzz =cyyyy0*czz0*(d042i + cyyyy1*d062i + czz1*d044i)
    uxxzzzz =cxx0*czzzz0*(d204i + cxx1*d404i + czzzz1*d206i)
    uyyzzzz =cyy0*czzzz0*(d024i + cyy1*d044i + czzzz1*d026i)
    uzzzzzz =czzzzzz0*(d006i + czzzzzz1*d008i)

    ! ----- spatial derivatives of order 8 ----
    uxxxxxxxx =cxxxxxxxx0*(d800i)
    uxxxxxxyy =cxxxxxx0*cyy0*(d620i)
    uxxxxyyyy =cxxxx0*cyyyy0*(d440i)
    uxxyyyyyy =cxx0*cyyyyyy0*(d260i)
    uyyyyyyyy =cyyyyyyyy0*(d080i)
    uxxxxxxzz =cxxxxxx0*czz0*(d602i)
    uxxxxyyzz =cxxxx0*cyy0*czz0*(d422i)
    uxxyyyyzz =cxx0*cyyyy0*czz0*(d242i)
    uyyyyyyzz =cyyyyyy0*czz0*(d062i)
    uxxxxzzzz =cxxxx0*czzzz0*(d404i)
    uxxyyzzzz =cxx0*cyy0*czzzz0*(d224i)
    uyyyyzzzz =cyyyy0*czzzz0*(d044i)
    uxxzzzzzz =cxx0*czzzzzz0*(d206i)
    uyyzzzzzz =cyy0*czzzzzz0*(d026i)
    uzzzzzzzz =czzzzzzzz0*(d008i)

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
                     + cdtPow8By20160*( uxxxxxxxx + uyyyyyyyy + uzzzzzzzz \
                                 + 4.*( uxxxxxxyy + uxxyyyyyy + uxxxxxxzz + uxxzzzzzz + uyyyyyyzz + uyyzzzzzz ) \
                                 + 6.*( uxxxxyyyy + uxxxxzzzz + uyyyyzzzz ) \
                                 +12.*( uxxxxyyzz + uxxyyyyzz + uxxyyzzzz ) )\
                       FV(m)  

  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 
#endMacro