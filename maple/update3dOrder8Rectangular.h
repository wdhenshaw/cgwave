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

! ---- DEFINE CONSTANTS IN EXPANSIONS OF DERIVATIVES ----
! Example: 
! u.xx = D+D-( I + cxx1*D+D- + cxx2*(D+D-x)^2 + ...
cx0 = 1.; cx1 = -1/6.; cx2 = 1/30.; cx3 = -1/140.; 
cy0 = 1.; cy1 = -1/6.; cy2 = 1/30.; cy3 = -1/140.; 
cz0 = 1.; cz1 = -1/6.; cz2 = 1/30.; cz3 = -1/140.; 
cxx0 = 1.; cxx1 = -1/12.; cxx2 = 1/90.; cxx3 = -1/560.; 
cyy0 = 1.; cyy1 = -1/12.; cyy2 = 1/90.; cyy3 = -1/560.; 
czz0 = 1.; czz1 = -1/12.; czz2 = 1/90.; czz3 = -1/560.; 
cxxx0 = 1.; cxxx1 = -1/4.; cxxx2 = 7/120.; cxxx3 = -41/3024.; 
cyyy0 = 1.; cyyy1 = -1/4.; cyyy2 = 7/120.; cyyy3 = -41/3024.; 
czzz0 = 1.; czzz1 = -1/4.; czzz2 = 7/120.; czzz3 = -41/3024.; 
cxxxx0 = 1.; cxxxx1 = -1/6.; cxxxx2 = 7/240.; cxxxx3 = -41/7560.; 
cyyyy0 = 1.; cyyyy1 = -1/6.; cyyyy2 = 7/240.; cyyyy3 = -41/7560.; 
czzzz0 = 1.; czzzz1 = -1/6.; czzzz2 = 7/240.; czzzz3 = -41/7560.; 
cxxxxx0 = 1.; cxxxxx1 = -1/3.; cxxxxx2 = 13/144.; cxxxxx3 = -139/6048.; 
cyyyyy0 = 1.; cyyyyy1 = -1/3.; cyyyyy2 = 13/144.; cyyyyy3 = -139/6048.; 
czzzzz0 = 1.; czzzzz1 = -1/3.; czzzzz2 = 13/144.; czzzzz3 = -139/6048.; 
cxxxxxx0 = 1.; cxxxxxx1 = -1/4.; cxxxxxx2 = 13/240.; cxxxxxx3 = -139/12096.; 
cyyyyyy0 = 1.; cyyyyyy1 = -1/4.; cyyyyyy2 = 13/240.; cyyyyyy3 = -139/12096.; 
czzzzzz0 = 1.; czzzzzz1 = -1/4.; czzzzzz2 = 13/240.; czzzzzz3 = -139/12096.; 
cxx=1./dx(0)**2;
cyy=1./dx(1)**2;
czz=1./dx(2)**2;
fv(m)=0.

numGhost1=3;
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

numGhost1=2;
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

numGhost1=1;
n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
beginLoops3d()
  #If #MASK eq "USEMASK" 
  if( mask(i1,i2,i3).ne.0 )then
  #End 
    d600(i1,i2,i3,0) = d400(i1+1,i2,i3,0) - 2*d400(i1,i2,i3,0) + d400(i1-1,i2,i3,0)
    d060(i1,i2,i3,0) = d040(i1,i2+1,i3,0) - 2*d040(i1,i2,i3,0) + d040(i1,i2-1,i3,0)
    d006(i1,i2,i3,0) = d004(i1,i2,i3+1,0) - 2*d004(i1,i2,i3,0) + d004(i1,i2,i3-1,0)
    ! --- Laplacian to order 6 = lap4h + corrections 
    lap6h(i1,i2,i3,0) = lap4h(i1,i2,i3,0) \
       + cxx*cxx2*d600(i1,i2,i3,0) \
       + cyy*cyy2*d060(i1,i2,i3,0) \
       + czz*czz2*d006(i1,i2,i3,0)   
    lap2hSq200(i1,i2,i3,0) = lap2hSq(i1+1,i2,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1-1,i2,i3,0)
    lap2hSq020(i1,i2,i3,0) = lap2hSq(i1,i2+1,i3,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2-1,i3,0)
    lap2hSq002(i1,i2,i3,0) = lap2hSq(i1,i2,i3+1,0) - 2*lap2hSq(i1,i2,i3,0) + lap2hSq(i1,i2,i3-1,0)
    lap2hCubed(i1,i2,i3,0) =                    \
         cxx*lap2hSq200(i1,i2,i3,0)  \
       + cyy*lap2hSq020(i1,i2,i3,0)  \
       + czz*lap2hSq020(i1,i2,i3,0)    

    ! --- Laplacian squared to order 4 = 
    !  lap2h*( lap4h ) + corrections*( Lap2h )
    lap4h200(i1,i2,i3,0) = lap4h(i1+1,i2,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1-1,i2,i3,0)
    lap4h020(i1,i2,i3,0) = lap4h(i1,i2+1,i3,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2-1,i3,0)
    lap4h002(i1,i2,i3,0) = lap4h(i1,i2,i3+1,0) - 2*lap4h(i1,i2,i3,0) + lap4h(i1,i2,i3-1,0)
    lap2h400(i1,i2,i3,0) = lap2h200(i1+1,i2,i3,0) - 2*lap2h200(i1,i2,i3,0) + lap2h200(i1-1,i2,i3,0)
    lap2h040(i1,i2,i3,0) = lap2h020(i1,i2+1,i3,0) - 2*lap2h020(i1,i2,i3,0) + lap2h020(i1,i2-1,i3,0)
    lap2h004(i1,i2,i3,0) = lap2h002(i1,i2,i3+1,0) - 2*lap2h002(i1,i2,i3,0) + lap2h002(i1,i2,i3-1,0)
    lap4hSq(i1,i2,i3,0) =                                       \
         cxx*( lap4h200(i1,i2,i3,0) + cxx1*lap2h400(i1,i2,i3,0) )  \
       + cyy*( lap4h020(i1,i2,i3,0) + cyy1*lap2h040(i1,i2,i3,0) )  \
       + czz*( lap4h002(i1,i2,i3,0) + czz1*lap2h004(i1,i2,i3,0) )    
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
    d800i = d600(i1+1,i2,i3,0) - 2*d600(i1,i2,i3,0) + d600(i1-1,i2,i3,0)
    d080i = d060(i1,i2+1,i3,0) - 2*d060(i1,i2,i3,0) + d060(i1,i2-1,i3,0)
    d008i = d006(i1,i2,i3+1,0) - 2*d006(i1,i2,i3,0) + d006(i1,i2,i3-1,0)
    ! --- Laplacian to order 8 = lap6h + corrections 
    lap8h = lap6h(i1,i2,i3,0)  \
            + cxx*cxx3*d800i  \
            + cyy*cyy3*d080i  \
            + czz*czz3*d008i    

    ! --- Laplacian^4 4p (4th power) order 2: 
    lap2hCubed200i = lap2hCubed(i1+1,i2,i3,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1-1,i2,i3,0)
    lap2hCubed020i = lap2hCubed(i1,i2+1,i3,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1,i2-1,i3,0)
    lap2hCubed002i = lap2hCubed(i1,i2,i3+1,0) - 2*lap2hCubed(i1,i2,i3,0) + lap2hCubed(i1,i2,i3-1,0)
    lap2h4p =                   \
            + cxx*lap2hCubed200i   \
            + cyy*lap2hCubed020i   \
            + czz*lap2hCubed002i     
    ! --- Laplacian squared to order 6 :
    !   Lap6h = Lap4h + M4  = (Lap2h) + M2 + M4 
    !   Lap6h*Lap6h = [ (Lap2h) + M2 + M4 ] [ (Lap2h) + M2 + M4 ]
    !               = Lap2h*Lap6h + M2*Lap4h + M4*Lap2h + O(h^6)
    lap6h200i = lap6h(i1+1,i2,i3,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1-1,i2,i3,0)
    lap6h020i = lap6h(i1,i2+1,i3,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1,i2-1,i3,0)
    lap6h002i = lap6h(i1,i2,i3+1,0) - 2*lap6h(i1,i2,i3,0) + lap6h(i1,i2,i3-1,0)
    lap4h400i = lap4h200(i1+1,i2,i3,0) - 2*lap4h200(i1,i2,i3,0) + lap4h200(i1-1,i2,i3,0)
    lap4h040i = lap4h020(i1,i2+1,i3,0) - 2*lap4h020(i1,i2,i3,0) + lap4h020(i1,i2-1,i3,0)
    lap4h004i = lap4h002(i1,i2,i3+1,0) - 2*lap4h002(i1,i2,i3,0) + lap4h002(i1,i2,i3-1,0)
    lap2h600i = lap2h400(i1+1,i2,i3,0) - 2*lap2h400(i1,i2,i3,0) + lap2h400(i1-1,i2,i3,0)
    lap2h060i = lap2h040(i1,i2+1,i3,0) - 2*lap2h040(i1,i2,i3,0) + lap2h040(i1,i2-1,i3,0)
    lap2h006i = lap2h004(i1,i2,i3+1,0) - 2*lap2h004(i1,i2,i3,0) + lap2h004(i1,i2,i3-1,0)
    lap6hSq =                                                    \
              cxx*(lap6h200i + cxx1*lap4h400i + cxx2*lap2h600i ) \
            + cyy*(lap6h020i + cyy1*lap4h040i + cyy2*lap2h060i ) \
            + czz*(lap6h002i + czz1*lap4h004i + czz2*lap2h006i )   

    ! --- Laplacian CUBED to order 4 
    lap4hSq200i = lap4hSq(i1+1,i2,i3,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1-1,i2,i3,0)
    lap4hSq020i = lap4hSq(i1,i2+1,i3,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1,i2-1,i3,0)
    lap4hSq002i = lap4hSq(i1,i2,i3+1,0) - 2*lap4hSq(i1,i2,i3,0) + lap4hSq(i1,i2,i3-1,0)
    lap2hSq400i = lap2hSq200(i1+1,i2,i3,0) - 2*lap2hSq200(i1,i2,i3,0) + lap2hSq200(i1-1,i2,i3,0)
    lap2hSq040i = lap2hSq020(i1,i2+1,i3,0) - 2*lap2hSq020(i1,i2,i3,0) + lap2hSq020(i1,i2-1,i3,0)
    lap2hSq004i = lap2hSq002(i1,i2,i3+1,0) - 2*lap2hSq002(i1,i2,i3,0) + lap2hSq002(i1,i2,i3-1,0)
    lap4hCubed =                                         \
                 cxx*(lap4hSq200i + cxx1*lap2hSq400i )   \
               + cyy*(lap4hSq020i + cyy1*lap2hSq040i )   \
               + czz*(lap4hSq002i + czz1*lap2hSq004i )     
    #If #FORCING ne "NOFORCING" 
      getForcing(DIM,ORDER,ORDERINTIME,GRIDTYPE) 
    #End 


    ! --- Modified equation space-time update ----
    un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m)-um(i1,i2,i3,m)  \
             + cdtsq*( lap8h )               \
             + cdtPow4By12*( lap6hSq )       \
             + cdtPow6By360*( lap4hCubed )   \
             + cdtPow8By20160*( lap2h4p )    \
             FV(m)                    
  #If #MASK eq "USEMASK" 
    end if ! mask .ne. 0
  #End 
  endLoops3d() 
#endMacro