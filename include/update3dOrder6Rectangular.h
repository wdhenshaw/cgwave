! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=3, gridType=Rectangular
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
!  MASK : USEMASK or NOMASK 
!  FORCING : NOFORCING, TZ, USEFORCING 
! ===========================================================================
#beginMacro update3dOrder6Rectangular(DIM,ORDER,ORDERINTIME,GRIDTYPE,MASK,FORCING)
#If #FORCING eq "USEFORCING"
#defineMacro FV(m) +dtSq*fv(m)
#Else
#defineMacro FV(m) 
#End

! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
! Example: 
! u.xx = D+D-( I + cxx1*D+D- + cxx2*(D+D-x)^2 + ...
cx0 = 1.; cx1 = -1/6.; cx2 = 1/30.; 
cy0 = 1.; cy1 = -1/6.; cy2 = 1/30.; 
cz0 = 1.; cz1 = -1/6.; cz2 = 1/30.; 
cxx0 = 1.; cxx1 = -1/12.; cxx2 = 1/90.; 
cyy0 = 1.; cyy1 = -1/12.; cyy2 = 1/90.; 
czz0 = 1.; czz1 = -1/12.; czz2 = 1/90.; 
cxxx0 = 1.; cxxx1 = -1/4.; cxxx2 = 7/120.; 
cyyy0 = 1.; cyyy1 = -1/4.; cyyy2 = 7/120.; 
czzz0 = 1.; czzz1 = -1/4.; czzz2 = 7/120.; 
cxxxx0 = 1.; cxxxx1 = -1/6.; cxxxx2 = 7/240.; 
cyyyy0 = 1.; cyyyy1 = -1/6.; cyyyy2 = 7/240.; 
czzzz0 = 1.; czzzz1 = -1/6.; czzzz2 = 7/240.; 
cxx=1./dx(0)**2;
cyy=1./dx(1)**2;
czz=1./dx(2)**2;
fv(m)=0.

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
    d002(i1,i2,i3,0) = u(i1,i2,i3+1,0) - 2*u(i1,i2,i3,0) + u(i1,i2,i3-1,0)
    lap2h(i1,i2,i3,0) = cxx*d200(i1,i2,i3,0) + cyy*d020(i1,i2,i3,0) + czz*d002(i1,i2,i3,0) 
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
    d004(i1,i2,i3,0) = d002(i1,i2,i3+1,0) - 2*d002(i1,i2,i3,0) + d002(i1,i2,i3-1,0)

    ! --- Laplacian to order 4 = lap2h + corrections 
    lap4h(i1,i2,i3,0) = lap2h(i1,i2,i3,0) \
         + cxx*cxx1*d400(i1,i2,i3,0) \
         + cyy*cyy1*d040(i1,i2,i3,0) \
         + czz*czz1*d004(i1,i2,i3,0)     

    ! --- Laplacian squared to order 2:
    lap2h200(i1,i2,i3,0) = lap2h(i1+1,i2,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1-1,i2,i3,0)
    lap2h020(i1,i2,i3,0) = lap2h(i1,i2+1,i3,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2-1,i3,0)
    lap2h002(i1,i2,i3,0) = lap2h(i1,i2,i3+1,0) - 2*lap2h(i1,i2,i3,0) + lap2h(i1,i2,i3-1,0)
    lap2hSq(i1,i2,i3,0) =                      \
             cxx*lap2h200(i1,i2,i3,0)   \
           + cyy*lap2h020(i1,i2,i3,0)   \
           + czz*lap2h002(i1,i2,i3,0)     
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
    d006i = d004(i1,i2,i3+1,0) - 2*d004(i1,i2,i3,0) + d004(i1,i2,i3-1,0)
    ! --- Laplacian to order 6 = lap4h + corrections 
    lap6h = lap4h(i1,i2,i3,0) \
       + cxx*cxx2*d600i \
       + cyy*cyy2*d060i \
       + czz*czz2*d006i   
    lap2hSq200i = lap2hSq(i1+1,i2,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1-1,i2,i3,0)
    lap2hSq020i = lap2hSq(i1,i2+1,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2-1,i3,0)
    lap2hSq002i = lap2hSq(i1,i2,i3+1,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2,i3-1,0)
    lap2hCubed =                    \
         cxx*lap2hSq200i  \
       + cyy*lap2hSq020i  \
       + czz*lap2hSq020i    

    ! --- Laplacian squared to order 4 = 
    !  lap2h*( lap4h ) + corrections*( Lap2h )
    lap4h200i = lap4h(i1+1,i2,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1-1,i2,i3,0)
    lap4h020i = lap4h(i1,i2+1,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2-1,i3,0)
    lap4h002i = lap4h(i1,i2,i3+1,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2,i3-1,0)
    lap2h400i = lap2h200(i1+1,i2,i3,0) - 2*lap2h200(i1,i2,i3,0) + lap2h200(i1-1,i2,i3,0)
    lap2h040i = lap2h020(i1,i2+1,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1,i2-1,i3,0)
    lap2h004i = lap2h002(i1,i2,i3+1,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1,i2,i3-1,0)
    lap4hSq =                                       \
         cxx*( lap4h200i + cxx1*lap2h400i )  \
       + cyy*( lap4h020i + cyy1*lap2h040i )  \
       + czz*( lap4h002i + czz1*lap2h004i )    
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