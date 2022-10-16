! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=3, gridType=Curvilinear
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
!  MASK : USEMASK or NOMASK 
!  FORCING : NOFORCING, TZ, USEFORCING 
! ===========================================================================
#beginMacro update3dOrder4Curvilinear(DIM,ORDER,ORDERINTIME,GRIDTYPE,MASK,FORCING)
#If #FORCING eq noForcing
#defineMacro FV(m) 
#Else
#defineMacro FV(m) +dtSq*fv(m)
#End

! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
! Example: 
! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...
cr0 = 1.; cr1 = -1/6.; 
cs0 = 1.; cs1 = -1/6.; 
ct0 = 1.; ct1 = -1/6.; 
crr0 = 1.; crr1 = -1/12.; 
css0 = 1.; css1 = -1/12.; 
ctt0 = 1.; ctt1 = -1/12.; 
dr1=dr(0); dr1i=1./dr1;
dr2=dr(1); dr2i=1./dr2;
dr3=dr(2); dr3i=1./dr3;
fv(m)=0.

#defineMacro c200(i1,i2,i3) lapCoeff(i1,i2,i3,0)
#defineMacro c020(i1,i2,i3) lapCoeff(i1,i2,i3,1)
#defineMacro c002(i1,i2,i3) lapCoeff(i1,i2,i3,2)
#defineMacro c110(i1,i2,i3) lapCoeff(i1,i2,i3,3)
#defineMacro c101(i1,i2,i3) lapCoeff(i1,i2,i3,4)
#defineMacro c011(i1,i2,i3) lapCoeff(i1,i2,i3,5)
#defineMacro c100(i1,i2,i3) lapCoeff(i1,i2,i3,6)
#defineMacro c010(i1,i2,i3) lapCoeff(i1,i2,i3,7)
#defineMacro c001(i1,i2,i3) lapCoeff(i1,i2,i3,8)

if( c200(0,0,0).le.0. )then

  ! --- Evaluate and store coefficients in Laplacian ---
  write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'

  numGhost1=1;
  n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
  n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
  n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
  beginLoops3d()
  #If #MASK eq "USEMASK" 
    if( mask(i1,i2,i3).ne.0 )then
  #End 
    rx = rsxy(i1,i2,i3,0,0)
    ry = rsxy(i1,i2,i3,0,1)
    rz = rsxy(i1,i2,i3,0,2)
    sx = rsxy(i1,i2,i3,1,0)
    sy = rsxy(i1,i2,i3,1,1)
    sz = rsxy(i1,i2,i3,1,2)
    tx = rsxy(i1,i2,i3,2,0)
    ty = rsxy(i1,i2,i3,2,1)
    tz = rsxy(i1,i2,i3,2,2)

    ! --- choose order for (r,s,t) derivatives based on available ghost points, less accuracy is needed in ghost points  ---
    if( (i1-2).ge.nd1a .and. (i1+2).le.nd1b )then
      diffOrder1=4
    elseif( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
      diffOrder1=2
    else
      stop 999
    end if
    if( (i2-2).ge.nd2a .and. (i2+2).le.nd2b )then
      diffOrder2=4
    elseif( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
      diffOrder2=2
    else
      stop 999
    end if
    if( (i3-2).ge.nd3a .and. (i3+2).le.nd3b )then
      diffOrder3=4
    elseif( (i3-1).ge.nd3a .and. (i3+1).le.nd3b )then
      diffOrder3=2
    else
      stop 999
    end if
    if( diffOrder1.eq.2 )then
      rxr = (rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0))*(.5*dr1i) 
      ryr = (rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1))*(.5*dr1i) 
      rzr = (rsxy(i1+1,i2,i3,0,2)-rsxy(i1-1,i2,i3,0,2))*(.5*dr1i) 
      sxr = (rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0))*(.5*dr1i) 
      syr = (rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1))*(.5*dr1i) 
      szr = (rsxy(i1+1,i2,i3,1,2)-rsxy(i1-1,i2,i3,1,2))*(.5*dr1i) 
      txr = (rsxy(i1+1,i2,i3,2,0)-rsxy(i1-1,i2,i3,2,0))*(.5*dr1i) 
      tyr = (rsxy(i1+1,i2,i3,2,1)-rsxy(i1-1,i2,i3,2,1))*(.5*dr1i) 
      tzr = (rsxy(i1+1,i2,i3,2,2)-rsxy(i1-1,i2,i3,2,2))*(.5*dr1i) 
    elseif( diffOrder1.eq.4 )then
      rxr = ( 8*(rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0)) -(rsxy(i1+2,i2,i3,0,0)-rsxy(i1-2,i2,i3,0,0)) )*(dr1i/12.) 
      ryr = ( 8*(rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1)) -(rsxy(i1+2,i2,i3,0,1)-rsxy(i1-2,i2,i3,0,1)) )*(dr1i/12.) 
      rzr = ( 8*(rsxy(i1+1,i2,i3,0,2)-rsxy(i1-1,i2,i3,0,2)) -(rsxy(i1+2,i2,i3,0,2)-rsxy(i1-2,i2,i3,0,2)) )*(dr1i/12.) 
      sxr = ( 8*(rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0)) -(rsxy(i1+2,i2,i3,1,0)-rsxy(i1-2,i2,i3,1,0)) )*(dr1i/12.) 
      syr = ( 8*(rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1)) -(rsxy(i1+2,i2,i3,1,1)-rsxy(i1-2,i2,i3,1,1)) )*(dr1i/12.) 
      szr = ( 8*(rsxy(i1+1,i2,i3,1,2)-rsxy(i1-1,i2,i3,1,2)) -(rsxy(i1+2,i2,i3,1,2)-rsxy(i1-2,i2,i3,1,2)) )*(dr1i/12.) 
      txr = ( 8*(rsxy(i1+1,i2,i3,2,0)-rsxy(i1-1,i2,i3,2,0)) -(rsxy(i1+2,i2,i3,2,0)-rsxy(i1-2,i2,i3,2,0)) )*(dr1i/12.) 
      tyr = ( 8*(rsxy(i1+1,i2,i3,2,1)-rsxy(i1-1,i2,i3,2,1)) -(rsxy(i1+2,i2,i3,2,1)-rsxy(i1-2,i2,i3,2,1)) )*(dr1i/12.) 
      tzr = ( 8*(rsxy(i1+1,i2,i3,2,2)-rsxy(i1-1,i2,i3,2,2)) -(rsxy(i1+2,i2,i3,2,2)-rsxy(i1-2,i2,i3,2,2)) )*(dr1i/12.) 
    end if
    if( diffOrder2.eq.2 )then
      rxs = (rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0))*(.5*dr2i) 
      rys = (rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1))*(.5*dr2i) 
      rzs = (rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2))*(.5*dr2i) 
      sxs = (rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0))*(.5*dr2i) 
      sys = (rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1))*(.5*dr2i) 
      szs = (rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2))*(.5*dr2i) 
      txs = (rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0))*(.5*dr2i) 
      tys = (rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1))*(.5*dr2i) 
      tzs = (rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2))*(.5*dr2i) 
    elseif( diffOrder2.eq.4 )then
      rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
      rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
      rzs = ( 8*(rsxy(i1,i2+1,i3,0,2)-rsxy(i1,i2-1,i3,0,2)) -(rsxy(i1,i2+2,i3,0,2)-rsxy(i1,i2-2,i3,0,2)) )*(dr2i/12.) 
      sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
      sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
      szs = ( 8*(rsxy(i1,i2+1,i3,1,2)-rsxy(i1,i2-1,i3,1,2)) -(rsxy(i1,i2+2,i3,1,2)-rsxy(i1,i2-2,i3,1,2)) )*(dr2i/12.) 
      txs = ( 8*(rsxy(i1,i2+1,i3,2,0)-rsxy(i1,i2-1,i3,2,0)) -(rsxy(i1,i2+2,i3,2,0)-rsxy(i1,i2-2,i3,2,0)) )*(dr2i/12.) 
      tys = ( 8*(rsxy(i1,i2+1,i3,2,1)-rsxy(i1,i2-1,i3,2,1)) -(rsxy(i1,i2+2,i3,2,1)-rsxy(i1,i2-2,i3,2,1)) )*(dr2i/12.) 
      tzs = ( 8*(rsxy(i1,i2+1,i3,2,2)-rsxy(i1,i2-1,i3,2,2)) -(rsxy(i1,i2+2,i3,2,2)-rsxy(i1,i2-2,i3,2,2)) )*(dr2i/12.) 
    end if
    if( diffOrder3.eq.2 )then
      rxt = (rsxy(i1,i2,i3+1,0,0)-rsxy(i1,i2,i3-1,0,0))*(.5*dr2i) 
      ryt = (rsxy(i1,i2,i3+1,0,1)-rsxy(i1,i2,i3-1,0,1))*(.5*dr2i) 
      rzt = (rsxy(i1,i2,i3+1,0,2)-rsxy(i1,i2,i3-1,0,2))*(.5*dr2i) 
      sxt = (rsxy(i1,i2,i3+1,1,0)-rsxy(i1,i2,i3-1,1,0))*(.5*dr2i) 
      syt = (rsxy(i1,i2,i3+1,1,1)-rsxy(i1,i2,i3-1,1,1))*(.5*dr2i) 
      szt = (rsxy(i1,i2,i3+1,1,2)-rsxy(i1,i2,i3-1,1,2))*(.5*dr2i) 
      txt = (rsxy(i1,i2,i3+1,2,0)-rsxy(i1,i2,i3-1,2,0))*(.5*dr2i) 
      tyt = (rsxy(i1,i2,i3+1,2,1)-rsxy(i1,i2,i3-1,2,1))*(.5*dr2i) 
      tzt = (rsxy(i1,i2,i3+1,2,2)-rsxy(i1,i2,i3-1,2,2))*(.5*dr2i) 
    elseif( diffOrder3.eq.4 )then
      rxt = ( 8*(rsxy(i1,i2,i3+1,0,0)-rsxy(i1,i2,i3-1,0,0)) -(rsxy(i1,i2,i3+2,0,0)-rsxy(i1,i2,i3-2,0,0)) )*(dr3i/12.) 
      ryt = ( 8*(rsxy(i1,i2,i3+1,0,1)-rsxy(i1,i2,i3-1,0,1)) -(rsxy(i1,i2,i3+2,0,1)-rsxy(i1,i2,i3-2,0,1)) )*(dr3i/12.) 
      rzt = ( 8*(rsxy(i1,i2,i3+1,0,2)-rsxy(i1,i2,i3-1,0,2)) -(rsxy(i1,i2,i3+2,0,2)-rsxy(i1,i2,i3-2,0,2)) )*(dr3i/12.) 
      sxt = ( 8*(rsxy(i1,i2,i3+1,1,0)-rsxy(i1,i2,i3-1,1,0)) -(rsxy(i1,i2,i3+2,1,0)-rsxy(i1,i2,i3-2,1,0)) )*(dr3i/12.) 
      syt = ( 8*(rsxy(i1,i2,i3+1,1,1)-rsxy(i1,i2,i3-1,1,1)) -(rsxy(i1,i2,i3+2,1,1)-rsxy(i1,i2,i3-2,1,1)) )*(dr3i/12.) 
      szt = ( 8*(rsxy(i1,i2,i3+1,1,2)-rsxy(i1,i2,i3-1,1,2)) -(rsxy(i1,i2,i3+2,1,2)-rsxy(i1,i2,i3-2,1,2)) )*(dr3i/12.) 
      txt = ( 8*(rsxy(i1,i2,i3+1,2,0)-rsxy(i1,i2,i3-1,2,0)) -(rsxy(i1,i2,i3+2,2,0)-rsxy(i1,i2,i3-2,2,0)) )*(dr3i/12.) 
      tyt = ( 8*(rsxy(i1,i2,i3+1,2,1)-rsxy(i1,i2,i3-1,2,1)) -(rsxy(i1,i2,i3+2,2,1)-rsxy(i1,i2,i3-2,2,1)) )*(dr3i/12.) 
      tzt = ( 8*(rsxy(i1,i2,i3+1,2,2)-rsxy(i1,i2,i3-1,2,2)) -(rsxy(i1,i2,i3+2,2,2)-rsxy(i1,i2,i3-2,2,2)) )*(dr3i/12.) 
    end if
    rxx = rx*rxr + sx*rxs + tx*rxt
    ryy = ry*ryr + sy*rys + ty*ryt
    rzz = rz*rzr + sz*rzs + tz*rzt
    sxx = rx*sxr + sx*sxs + tx*sxt
    syy = ry*syr + sy*sys + ty*syt
    szz = rz*szr + sz*szs + tz*szt
    txx = rx*txr + sx*txs + tx*txt
    tyy = ry*tyr + sy*tys + ty*tyt
    tzz = rz*tzr + sz*tzs + tz*tzt

    ! -- Coefficients in the Laplacian (scaled)
    c200(i1,i2,i3) = (rx**2 + ry**2 + rz**2 )*dr1i**2
    c020(i1,i2,i3) = (sx**2 + sy**2 + sz**2 )*dr2i**2
    c002(i1,i2,i3) = (tx**2 + ty**2 + tz**2 )*dr3i**2
    c110(i1,i2,i3) = 2.*(rx*sx + ry*sy + rz*sz )*dr1i*dr2i*.25
    c101(i1,i2,i3) = 2.*(rx*tx + ry*ty + rz*tz )*dr1i*dr2i*.25
    c011(i1,i2,i3) = 2.*(sx*tx + sy*ty + sz*tz )*dr1i*dr2i*.25
    c100(i1,i2,i3) = (rxx + ryy + rzz)*dr1i*.5
    c010(i1,i2,i3) = (sxx + syy + tyy)*dr2i*.5 
    c001(i1,i2,i3) = (txx + tyy + tzz)*dr3i*.5 

  #If #MASK eq "USEMASK" 
    end if ! mask .ne. 0
  #End 
  endLoops3d() 
end if ! end assignLapCoeff

numGhost1=1;
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 
    d200(i1,i2,i3,0) = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
    d020(i1,i2,i3,0) = u(i1,i2+1,i3,0) - 2*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
    d002(i1,i2,i3,0) = u(i1,i2,i3+1,0) - 2*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
    d100i = u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0)
    d010i = u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0)
    d110i = u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0)
    d001i = u(i1,i2,i3+1,0) - u(i1,i2,i3-1,0)
    d101i = u(i1+1,i2,i3+1,0) - u(i1-1,i2,i3+1,0) - u(i1+1,i2,i3-1,0) + u(i1-1,i2,i3-1,0)
    d011i = u(i1,i2+1,i3+1,0) - u(i1,i2-1,i3+1,0) - u(i1,i2+1,i3-1,0) + u(i1,i2-1,i3-1,0)
    lap2h(i1,i2,i3,0) = \
         c200(i1,i2,i3)*d200(i1,i2,i3,0) +\
         c020(i1,i2,i3)*d020(i1,i2,i3,0) +\
         c002(i1,i2,i3)*d002(i1,i2,i3,0) +\
         c110(i1,i2,i3)*d110i + \
         c101(i1,i2,i3)*d101i + \
         c011(i1,i2,i3)*d011i + \
         c100(i1,i2,i3)*d100i + \
         c010(i1,i2,i3)*d010i + \
         c001(i1,i2,i3)*d001i
  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 

! ===========  FINAL LOOP TO FILL IN THE SOLUTION ============

numGhost1=0;
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 
    d400i = d200(i1+1,i2,i3,0) - 2*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0)
    d040i = d020(i1,i2+1,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1,i2-1,i3,0)
    d004i = d002(i1,i2,i3+1,0) - 2*d002(i1,i2,i3,0) + d002(i1,i2,i3-1,0)
    d300i = d200(i1+1,i2,i3,0) - d200(i1-1,i2,i3,0)
    d030i = d020(i1,i2+1,i3,0) - d020(i1,i2-1,i3,0)
    d003i = d002(i1,i2,i3+1,0) - d002(i1,i2,i3-1,0)
    d310i = d200(i1+1,i2+1,i3,0) - d200(i1-1,i2+1,i3,0) - d200(i1+1,i2-1,i3,0) + d200(i1-1,i2-1,i3,0)
    d130i = d020(i1+1,i2+1,i3,0) - d020(i1-1,i2+1,i3,0) - d020(i1+1,i2-1,i3,0) + d020(i1-1,i2-1,i3,0)
    d301i = d200(i1+1,i2,i3+1,0) - d200(i1-1,i2,i3+1,0) - d200(i1+1,i2,i3-1,0) + d200(i1-1,i2,i3-1,0)
    d103i = d002(i1+1,i2,i3+1,0) - d002(i1-1,i2,i3+1,0) - d002(i1+1,i2,i3-1,0) + d002(i1-1,i2,i3-1,0)
    d031i = d020(i1,i2+1,i3+1,0) - d020(i1,i2-1,i3+1,0) - d020(i1,i2+1,i3-1,0) + d020(i1,i2-1,i3-1,0)
    d013i = d002(i1,i2+1,i3+1,0) - d002(i1,i2-1,i3+1,0) - d002(i1,i2+1,i3-1,0) + d002(i1,i2-1,i3-1,0)

    ! --- Laplacian to order 4 = lap2h + corrections 
    lap4h = lap2h(i1,i2,i3,0) \
         + c200(i1,i2,i3)*crr1*d400i \
         + c020(i1,i2,i3)*css1*d040i \
         + c002(i1,i2,i3)*ctt1*d004i \
         + c110(i1,i2,i3)*(cr1*d310i + cs1*d130i) \
         + c101(i1,i2,i3)*(cr1*d301i + ct1*d103i) \
         + c011(i1,i2,i3)*(cs1*d031i + ct1*d013i) \
         + c100(i1,i2,i3)*cr1 *d300i \
         + c010(i1,i2,i3)*cs1 *d030i \
         + c001(i1,i2,i3)*ct1 *d003i 

    ! --- Laplacian squared to order 2:
    lap2h200i = lap2h(i1+1,i2,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1-1,i2,i3,0)
    lap2h020i = lap2h(i1,i2+1,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2-1,i3,0)
    lap2h002i = lap2h(i1,i2,i3+1,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2,i3-1,0)
    lap2h100i = lap2h(i1+1,i2,i3,0) - lap2h(i1-1,i2,i3,0)
    lap2h010i = lap2h(i1,i2+1,i3,0) - lap2h(i1,i2-1,i3,0)
    lap2h110i = lap2h(i1+1,i2+1,i3,0) - lap2h(i1-1,i2+1,i3,0) - lap2h(i1+1,i2-1,i3,0) + lap2h(i1-1,i2-1,i3,0)
    lap2h001i = lap2h(i1,i2,i3+1,0) - lap2h(i1,i2,i3-1,0)
    lap2h101i = lap2h(i1+1,i2,i3+1,0) - lap2h(i1-1,i2,i3+1,0) - lap2h(i1+1,i2,i3-1,0) + lap2h(i1-1,i2,i3-1,0)
    lap2h011i = lap2h(i1,i2+1,i3+1,0) - lap2h(i1,i2-1,i3+1,0) - lap2h(i1,i2+1,i3-1,0) + lap2h(i1,i2-1,i3-1,0)
    lap2hSq =  \
            c200(i1,i2,i3)*lap2h200i \
          + c020(i1,i2,i3)*lap2h020i \
          + c002(i1,i2,i3)*lap2h002i \
          + c110(i1,i2,i3)*lap2h110i  \
          + c101(i1,i2,i3)*lap2h101i  \
          + c011(i1,i2,i3)*lap2h011i  \
          + c100(i1,i2,i3)*lap2h100i  \
          + c010(i1,i2,i3)*lap2h010i  \
          + c001(i1,i2,i3)*lap2h001i    
    #If #FORCING ne "NOFORCING" 
      getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) 
    #End 


    ! --- Modified equation space-time update ----
    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m) \
             + cdtsq*( lap4h )                         \
             + cdtPow4By12*( lap2hSq )                   \
             FV(m)                                    
  #If #MASK eq "USEMASK" 
  end if ! mask .ne. 0
  #End 
endLoops3d() 
#endMacro