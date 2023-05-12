! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=2, gridType=Curvilinear
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
!  MASK : USEMASK or NOMASK 
!  FORCING : NOFORCING, TZ, USEFORCING 
! ===========================================================================
#beginMacro update2dOrder6Curvilinear(DIM,ORDER,ORDERINTIME,GRIDTYPE,MASK,FORCING)
#If #FORCING eq "USEFORCING"
#defineMacro FV(m) +dtSq*fv(m)
#Else
#defineMacro FV(m) 
#End

! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
! Example: 
! u.rr = D+D-( I + crr1*D+D- + crr2*(D+D-x)^2 + ...
cr0 = 1.; cr1 = -1/6.; cr2 = 1/30.; 
cs0 = 1.; cs1 = -1/6.; cs2 = 1/30.; 
crr0 = 1.; crr1 = -1/12.; crr2 = 1/90.; 
css0 = 1.; css1 = -1/12.; css2 = 1/90.; 
crrr0 = 1.; crrr1 = -1/4.; crrr2 = 7/120.; 
csss0 = 1.; csss1 = -1/4.; csss2 = 7/120.; 
crrrr0 = 1.; crrrr1 = -1/6.; crrrr2 = 7/240.; 
cssss0 = 1.; cssss1 = -1/6.; cssss2 = 7/240.; 
dr1=dr(0); dr1i=1./dr1;
dr2=dr(1); dr2i=1./dr2;
dr3=dr(2); dr3i=1./dr3;
fv(m)=0.

#defineMacro c200(i1,i2,i3) lapCoeff(i1,i2,i3,0)
#defineMacro c020(i1,i2,i3) lapCoeff(i1,i2,i3,1)
#defineMacro c110(i1,i2,i3) lapCoeff(i1,i2,i3,2)
#defineMacro c100(i1,i2,i3) lapCoeff(i1,i2,i3,3)
#defineMacro c010(i1,i2,i3) lapCoeff(i1,i2,i3,4)

if( c200(0,0,0).le.0. )then

  ! --- Evaluate and store coefficients in Laplacian ---
  write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'

  numGhost1=2;
  n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
  n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
  n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 
    rx = rsxy(i1,i2,i3,0,0)
    ry = rsxy(i1,i2,i3,0,1)
    sx = rsxy(i1,i2,i3,1,0)
    sy = rsxy(i1,i2,i3,1,1)

    ! --- choose order for (r,s,t) derivatives based on available ghost points, less accuracy is needed in ghost points  ---
    if( (i1-3).ge.nd1a .and. (i1+3).le.nd1b )then
      diffOrder1=6
    elseif( (i1-2).ge.nd1a .and. (i1+2).le.nd1b )then
      diffOrder1=4
    elseif( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
      diffOrder1=2
    else
      stop 999
    end if
    if( (i2-3).ge.nd2a .and. (i2+3).le.nd2b )then
      diffOrder2=6
    elseif( (i2-2).ge.nd2a .and. (i2+2).le.nd2b )then
      diffOrder2=4
    elseif( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
      diffOrder2=2
    else
      stop 999
    end if
    if( diffOrder1.eq.2 )then
      rxr = (rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0))*(.5*dr1i) 
      ryr = (rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1))*(.5*dr1i) 
      sxr = (rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0))*(.5*dr1i) 
      syr = (rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1))*(.5*dr1i) 
    elseif( diffOrder1.eq.4 )then
      rxr = ( 8*(rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0)) -(rsxy(i1+2,i2,i3,0,0)-rsxy(i1-2,i2,i3,0,0)) )*(dr1i/12.) 
      ryr = ( 8*(rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1)) -(rsxy(i1+2,i2,i3,0,1)-rsxy(i1-2,i2,i3,0,1)) )*(dr1i/12.) 
      sxr = ( 8*(rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0)) -(rsxy(i1+2,i2,i3,1,0)-rsxy(i1-2,i2,i3,1,0)) )*(dr1i/12.) 
      syr = ( 8*(rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1)) -(rsxy(i1+2,i2,i3,1,1)-rsxy(i1-2,i2,i3,1,1)) )*(dr1i/12.) 
    elseif( diffOrder1.eq.6 )then
      rxr = ( 45.*(rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0)) -9.*(rsxy(i1+2,i2,i3,0,0)-rsxy(i1-2,i2,i3,0,0)) +(rsxy(i1+3,i2,i3,0,0)-rsxy(i1-3,i2,i3,0,0)) )*(dr1i/60.) 
      ryr = ( 45.*(rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1)) -9.*(rsxy(i1+2,i2,i3,0,1)-rsxy(i1-2,i2,i3,0,1)) +(rsxy(i1+3,i2,i3,0,1)-rsxy(i1-3,i2,i3,0,1)) )*(dr1i/60.) 
      sxr = ( 45.*(rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0)) -9.*(rsxy(i1+2,i2,i3,1,0)-rsxy(i1-2,i2,i3,1,0)) +(rsxy(i1+3,i2,i3,1,0)-rsxy(i1-3,i2,i3,1,0)) )*(dr1i/60.) 
      syr = ( 45.*(rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1)) -9.*(rsxy(i1+2,i2,i3,1,1)-rsxy(i1-2,i2,i3,1,1)) +(rsxy(i1+3,i2,i3,1,1)-rsxy(i1-3,i2,i3,1,1)) )*(dr1i/60.) 
    end if
    if( diffOrder2.eq.2 )then
      rxs = (rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0))*(.5*dr2i) 
      rys = (rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1))*(.5*dr2i) 
      sxs = (rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0))*(.5*dr2i) 
      sys = (rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1))*(.5*dr2i) 
    elseif( diffOrder2.eq.4 )then
      rxs = ( 8*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) )*(dr2i/12.) 
      rys = ( 8*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) )*(dr2i/12.) 
      sxs = ( 8*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) )*(dr2i/12.) 
      sys = ( 8*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) )*(dr2i/12.) 
    elseif( diffOrder2.eq.6 )then
      rxs = ( 45.*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -9.*(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) +(rsxy(i1,i2+3,i3,0,0)-rsxy(i1,i2-3,i3,0,0)) )*(dr2i/60.) 
      rys = ( 45.*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -9.*(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) +(rsxy(i1,i2+3,i3,0,1)-rsxy(i1,i2-3,i3,0,1)) )*(dr2i/60.) 
      sxs = ( 45.*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -9.*(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) +(rsxy(i1,i2+3,i3,1,0)-rsxy(i1,i2-3,i3,1,0)) )*(dr2i/60.) 
      sys = ( 45.*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -9.*(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) +(rsxy(i1,i2+3,i3,1,1)-rsxy(i1,i2-3,i3,1,1)) )*(dr2i/60.) 
    end if
    rxx = rx*rxr + sx*rxs 
    ryy = ry*ryr + sy*rys 
    sxx = rx*sxr + sx*sxs 
    syy = ry*syr + sy*sys 

    ! -- Coefficients in the Laplacian (scaled)
    c200(i1,i2,i3) = (rx**2 + ry**2   )*dr1i**2
    c110(i1,i2,i3) = 2.*(rx*sx + ry*sy)*dr1i*dr2i*.25
    c020(i1,i2,i3) = (sx**2 + sy**2   )*dr2i**2
    c100(i1,i2,i3) = (rxx + ryy       )*dr1i*.5
    c010(i1,i2,i3) = (sxx + syy       )*dr2i*.5 

  #If #MASK eq "USEMASK" 
    end if ! mask .ne. 0
  #End 
  endLoops3d() 
end if ! end assignLapCoeff

numGhost1=2;
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 
    d200(i1,i2,i3,0) = u(i1+1,i2,i3,0) - 2*u(i1,i2,i3,0) + u(i1-1,i2,i3,0)
    d020(i1,i2,i3,0) = u(i1,i2+1,i3,0) - 2*u(i1,i2,i3,0) + u(i1,i2-1,i3,0)
    d100i = u(i1+1,i2,i3,0) - u(i1-1,i2,i3,0)
    d010i = u(i1,i2+1,i3,0) - u(i1,i2-1,i3,0)
    d110i = u(i1+1,i2+1,i3,0) - u(i1-1,i2+1,i3,0) - u(i1+1,i2-1,i3,0) + u(i1-1,i2-1,i3,0)
    lap2h(i1,i2,i3,0) = \
         c200(i1,i2,i3)*d200(i1,i2,i3,0) + \
         c110(i1,i2,i3)*d110i + \
         c020(i1,i2,i3)*d020(i1,i2,i3,0) +\
         c100(i1,i2,i3)*d100i + \
         c010(i1,i2,i3)*d010i
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
    d400(i1,i2,i3,0) = d200(i1+1,i2,i3,0) - 2*d200(i1,i2,i3,0) + d200(i1-1,i2,i3,0)
    d040(i1,i2,i3,0) = d020(i1,i2+1,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1,i2-1,i3,0)
    d220(i1,i2,i3,0) = d020(i1+1,i2,i3,0) - 2*d020(i1,i2,i3,0) + d020(i1-1,i2,i3,0)
    d300i = d200(i1+1,i2,i3,0) - d200(i1-1,i2,i3,0)
    d030i = d020(i1,i2+1,i3,0) - d020(i1,i2-1,i3,0)
    d310i = d200(i1+1,i2+1,i3,0) - d200(i1-1,i2+1,i3,0) - d200(i1+1,i2-1,i3,0) + d200(i1-1,i2-1,i3,0)
    d130i = d020(i1+1,i2+1,i3,0) - d020(i1-1,i2+1,i3,0) - d020(i1+1,i2-1,i3,0) + d020(i1-1,i2-1,i3,0)

    ! --- Laplacian to order 4 = lap2h + corrections 
    lap4h(i1,i2,i3,0) = lap2h(i1,i2,i3,0) \
         + c200(i1,i2,i3)*crr1*d400(i1,i2,i3,0) \
         + c110(i1,i2,i3)*(cr1*d310i + cs1*d130i) \
         + c020(i1,i2,i3)*css1*d040(i1,i2,i3,0) \
         + c100(i1,i2,i3)*cr1 *d300i \
         + c010(i1,i2,i3)*cs1 *d030i 

    ! --- Laplacian squared to order 2:
    lap2h200(i1,i2,i3,0) = lap2h(i1+1,i2,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1-1,i2,i3,0)
    lap2h020(i1,i2,i3,0) = lap2h(i1,i2+1,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2-1,i3,0)
    lap2h100i = lap2h(i1+1,i2,i3,0) - lap2h(i1-1,i2,i3,0)
    lap2h010i = lap2h(i1,i2+1,i3,0) - lap2h(i1,i2-1,i3,0)
    lap2h110i = lap2h(i1+1,i2+1,i3,0) - lap2h(i1-1,i2+1,i3,0) - lap2h(i1+1,i2-1,i3,0) + lap2h(i1-1,i2-1,i3,0)
    lap2hSq(i1,i2,i3,0) =  \
             c200(i1,i2,i3)*lap2h200(i1,i2,i3,0) \
           + c020(i1,i2,i3)*lap2h020(i1,i2,i3,0) \
           + c110(i1,i2,i3)*lap2h110i  \
           + c100(i1,i2,i3)*lap2h100i  \
           + c010(i1,i2,i3)*lap2h010i    
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
    d600i = d400(i1+1,i2,i3,0) - 2*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0)
    d060i = d040(i1,i2+1,i3,0) - 2*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0)
    d500i = d400(i1+1,i2,i3,0) - d400(i1-1,i2,i3,0)
    d050i = d040(i1,i2+1,i3,0) - d040(i1,i2-1,i3,0)
    d510i = d400(i1+1,i2+1,i3,0) - d400(i1-1,i2+1,i3,0) - d400(i1+1,i2-1,i3,0) + d400(i1-1,i2-1,i3,0)
    d150i = d040(i1+1,i2+1,i3,0) - d040(i1-1,i2+1,i3,0) - d040(i1+1,i2-1,i3,0) + d040(i1-1,i2-1,i3,0)
    d330i = d220(i1+1,i2+1,i3,0) - d220(i1-1,i2+1,i3,0) - d220(i1+1,i2-1,i3,0) + d220(i1-1,i2-1,i3,0)
    ! --- Laplacian to order 6 = lap4h + corrections 
    lap6h = lap4h(i1,i2,i3,0) \
       + c200(i1,i2,i3)*crr2*d600i \
       + c110(i1,i2,i3)*(cr2*d510i + cs2*d150i + cr1*cs1*d330i ) \
       + c020(i1,i2,i3)*css2*d060i \
       + c100(i1,i2,i3)*cr2 *d500i \
       + c010(i1,i2,i3)*cs2 *d050i 
    lap2hSq200i = lap2hSq(i1+1,i2,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1-1,i2,i3,0)
    lap2hSq020i = lap2hSq(i1,i2+1,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2-1,i3,0)
    lap2hSq100i = lap2hSq(i1+1,i2,i3,0) - lap2hSq(i1-1,i2,i3,0)
    lap2hSq010i = lap2hSq(i1,i2+1,i3,0) - lap2hSq(i1,i2-1,i3,0)
    lap2hSq110i = lap2hSq(i1+1,i2+1,i3,0) - lap2hSq(i1-1,i2+1,i3,0) - lap2hSq(i1+1,i2-1,i3,0) + lap2hSq(i1-1,i2-1,i3,0)
    lap2hCubed =  \
       + c200(i1,i2,i3)*lap2hSq200i \
       + c110(i1,i2,i3)*lap2hSq110i  \
       + c020(i1,i2,i3)*lap2hSq020i \
       + c100(i1,i2,i3)*lap2hSq100i  \
       + c010(i1,i2,i3)*lap2hSq010i   

    ! --- Laplacian squared to order 4 = 
    !  lap2h*( lap4h ) + corrections*( Lap2h )
    lap4h200i = lap4h(i1+1,i2,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1-1,i2,i3,0)
    lap4h020i = lap4h(i1,i2+1,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2-1,i3,0)
    lap2h400i = lap2h200(i1+1,i2,i3,0) - 2*lap2h200(i1,i2,i3,0) + lap2h200(i1-1,i2,i3,0)
    lap2h040i = lap2h020(i1,i2+1,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1,i2-1,i3,0)
    lap4h100i = lap4h(i1+1,i2,i3,0) - lap4h(i1-1,i2,i3,0)
    lap4h010i = lap4h(i1,i2+1,i3,0) - lap4h(i1,i2-1,i3,0)
    lap4h110i = lap4h(i1+1,i2+1,i3,0) - lap4h(i1-1,i2+1,i3,0) - lap4h(i1+1,i2-1,i3,0) + lap4h(i1-1,i2-1,i3,0)
    lap2h300i = lap2h200(i1+1,i2,i3,0) - lap2h200(i1-1,i2,i3,0)
    lap2h030i = lap2h020(i1,i2+1,i3,0) - lap2h020(i1,i2-1,i3,0)
    lap2h310i = lap2h200(i1+1,i2+1,i3,0) - lap2h200(i1-1,i2+1,i3,0) - lap2h200(i1+1,i2-1,i3,0) + lap2h200(i1-1,i2-1,i3,0)
    lap2h130i = lap2h020(i1+1,i2+1,i3,0) - lap2h020(i1-1,i2+1,i3,0) - lap2h020(i1+1,i2-1,i3,0) + lap2h020(i1-1,i2-1,i3,0)
    lap4hSq =     \
         c200(i1,i2,i3)*( lap4h200i + crr1*lap2h400i )    \
       + c110(i1,i2,i3)*( lap4h110i + cr1*lap2h310i + cs1*lap2h130i ) \
       + c020(i1,i2,i3)*( lap4h020i + css1*lap2h040i )     \
       + c100(i1,i2,i3)*( lap4h100i + cr1 *lap2h300i )    \
       + c010(i1,i2,i3)*( lap4h010i + cs1 *lap2h030i )      
    #If #FORCING ne "NOFORCING" 
      getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) 
    #End 


    ! --- Modified equation space-time update ----
    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m)  \
             + cdtsq*( lap6h )               \
             + cdtPow4By12*( lap4hSq )       \
             + cdtPow6By360*( lap2hCubed )   \
             FV(m)                    
  #If #MASK eq "USEMASK" 
    end if ! mask .ne. 0
  #End 
  endLoops3d() 
#endMacro