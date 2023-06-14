! This file automatically generated from bcOptWave.bf90 with bpp.
! ==================================================================================
!
!        Optimized Assign Boundary Conditions for CgWave
!        -----------------------------------------------
!
! ==================================================================================


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



! define macros to evaluate higher derivatives (from maple/makeGetDerivativesMacros.maple)
!! ** June 13, 2023 : TURN OFF ??
!! #Include "../maple/defineGetDerivativesMacros.h"

! NEW VERSION WITH DISTINCTIVE NAMES:
! ****** File written by makeGetDerivativesMacros.maple  ******



! =======================================================
!  Macro to compute Third derivatives in 2 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================

! =======================================================
!  Macro to compute Third derivatives in 3 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================


! ! ************ TEMP ********************
! #beginMacro defineParameticDerivativesComponents0(u)
!  #defineMacro u ## Ar2(i1,i2,i3) (-u(i1-1,i2,i3)+u(i1+1,i2,i3))/(2.*dr(0))
!  #defineMacro u ## As2(i1,i2,i3) (-u(i1,i2-1,i3)+u(i1,i2+1,i3))/(2.*dr(1))
!  #defineMacro u ## At2(i1,i2,i3) (-u(i1,i2,i3-1)+u(i1,i2,i3+1))/(2.*dr(2))
!  #defineMacro u ## Arr2(i1,i2,i3) (u(i1-1,i2,i3)-2.*u(i1,i2,i3)+u(i1+1,i2,i3))/(dr(0)**2)
!  #defineMacro u ## Ars2(i1,i2,i3) (-u ## s2(i1-1,i2,i3)+u ## s2(i1+1,i2,i3))/(2.*dr(0))
!  #defineMacro u ## Ass2(i1,i2,i3) (u(i1,i2-1,i3)-2.*u(i1,i2,i3)+u(i1,i2+1,i3))/(dr(1)**2)
!  #defineMacro u ## Art2(i1,i2,i3) (-u ## t2(i1-1,i2,i3)+u ## t2(i1+1,i2,i3))/(2.*dr(0))
!  #defineMacro u ## Ast2(i1,i2,i3) (-u ## t2(i1,i2-1,i3)+u ## t2(i1,i2+1,i3))/(2.*dr(1))
!  #defineMacro u ## Att2(i1,i2,i3) (u(i1,i2,i3-1)-2.*u(i1,i2,i3)+u(i1,i2,i3+1))/(dr(2)**2)
!  #defineMacro u ## Arrr2(i1,i2,i3) (-u(i1-2,i2,i3)+2.*u(i1-1,i2,i3)-2.*u(i1+1,i2,i3)+u(i1+2,i2,i3))/(2.*dr(0)**3)
!  #defineMacro u ## Arrs2(i1,i2,i3) (u ## s2(i1-1,i2,i3)-2.*u ## s2(i1,i2,i3)+u ## s2(i1+1,i2,i3))/(dr(0)**2)
!  #defineMacro u ## Arss2(i1,i2,i3) (-u ## ss2(i1-1,i2,i3)+u ## ss2(i1+1,i2,i3))/(2.*dr(0))
!  #defineMacro u ## Asss2(i1,i2,i3) (-u(i1,i2-2,i3)+2.*u(i1,i2-1,i3)-2.*u(i1,i2+1,i3)+u(i1,i2+2,i3))/(2.*dr(1)**3)
!  #defineMacro u ## Arrt2(i1,i2,i3) (u ## t2(i1-1,i2,i3)-2.*u ## t2(i1,i2,i3)+u ## t2(i1+1,i2,i3))/(dr(0)**2)
!  #defineMacro u ## Arst2(i1,i2,i3) (-u ## st2(i1-1,i2,i3)+u ## st2(i1+1,i2,i3))/(2.*dr(0))
!  #defineMacro u ## Asst2(i1,i2,i3) (u ## t2(i1,i2-1,i3)-2.*u ## t2(i1,i2,i3)+u ## t2(i1,i2+1,i3))/(dr(1)**2)
!  #defineMacro u ## Artt2(i1,i2,i3) (-u ## tt2(i1-1,i2,i3)+u ## tt2(i1+1,i2,i3))/(2.*dr(0))
!  #defineMacro u ## Astt2(i1,i2,i3) (-u ## tt2(i1,i2-1,i3)+u ## tt2(i1,i2+1,i3))/(2.*dr(1))
!  #defineMacro u ## Attt2(i1,i2,i3) (-u(i1,i2,i3-2)+2.*u(i1,i2,i3-1)-2.*u(i1,i2,i3+1)+u(i1,i2,i3+2))/(2.*dr(2)**3)
!  #defineMacro u ## Ar4(i1,i2,i3) (u(i1-2,i2,i3)-8.*u(i1-1,i2,i3)+8.*u(i1+1,i2,i3)-u(i1+2,i2,i3))/(12.*dr(0))
!  #defineMacro u ## As4(i1,i2,i3) (u(i1,i2-2,i3)-8.*u(i1,i2-1,i3)+8.*u(i1,i2+1,i3)-u(i1,i2+2,i3))/(12.*dr(1))
!  #defineMacro u ## At4(i1,i2,i3) (u(i1,i2,i3-2)-8.*u(i1,i2,i3-1)+8.*u(i1,i2,i3+1)-u(i1,i2,i3+2))/(12.*dr(2))
!  #defineMacro u ## Arr4(i1,i2,i3) (-u(i1-2,i2,i3)+16.*u(i1-1,i2,i3)-30.*u(i1,i2,i3)+16.*u(i1+1,i2,i3)-u(i1+2,i2,i3))/(12.*dr(0)**2)
!  #defineMacro u ## Ars4(i1,i2,i3) (u ## s4(i1-2,i2,i3)-8.*u ## s4(i1-1,i2,i3)+8.*u ## s4(i1+1,i2,i3)-u ## s4(i1+2,i2,i3))/(12.*dr(0))
!  #defineMacro u ## Ass4(i1,i2,i3) (-u(i1,i2-2,i3)+16.*u(i1,i2-1,i3)-30.*u(i1,i2,i3)+16.*u(i1,i2+1,i3)-u(i1,i2+2,i3))/(12.*dr(1)**2)
!  #defineMacro u ## Art4(i1,i2,i3) (u ## t4(i1-2,i2,i3)-8.*u ## t4(i1-1,i2,i3)+8.*u ## t4(i1+1,i2,i3)-u ## t4(i1+2,i2,i3))/(12.*dr(0))
!  #defineMacro u ## Ast4(i1,i2,i3) (u ## t4(i1,i2-2,i3)-8.*u ## t4(i1,i2-1,i3)+8.*u ## t4(i1,i2+1,i3)-u ## t4(i1,i2+2,i3))/(12.*dr(1))
!  #defineMacro u ## Att4(i1,i2,i3) (-u(i1,i2,i3-2)+16.*u(i1,i2,i3-1)-30.*u(i1,i2,i3)+16.*u(i1,i2,i3+1)-u(i1,i2,i3+2))/(12.*dr(2)**2)
!  #defineMacro u ## Arrr4(i1,i2,i3) (u(i1-3,i2,i3)-8.*u(i1-2,i2,i3)+13.*u(i1-1,i2,i3)-13.*u(i1+1,i2,i3)+8.*u(i1+2,i2,i3)-u(i1+3,i2,i3))/(8.*dr(0)**3)
!  #defineMacro u ## Arrs4(i1,i2,i3) (-u ## s4(i1-2,i2,i3)+16.*u ## s4(i1-1,i2,i3)-30.*u ## s4(i1,i2,i3)+16.*u ## s4(i1+1,i2,i3)-u ## s4(i1+2,i2,i3))/(12.*dr(0)**2)
!  #defineMacro u ## Arss4(i1,i2,i3) (u ## ss4(i1-2,i2,i3)-8.*u ## ss4(i1-1,i2,i3)+8.*u ## ss4(i1+1,i2,i3)-u ## ss4(i1+2,i2,i3))/(12.*dr(0))
!  #defineMacro u ## Asss4(i1,i2,i3) (u(i1,i2-3,i3)-8.*u(i1,i2-2,i3)+13.*u(i1,i2-1,i3)-13.*u(i1,i2+1,i3)+8.*u(i1,i2+2,i3)-u(i1,i2+3,i3))/(8.*dr(1)**3)
!  #defineMacro u ## Arrt4(i1,i2,i3) (-u ## t4(i1-2,i2,i3)+16.*u ## t4(i1-1,i2,i3)-30.*u ## t4(i1,i2,i3)+16.*u ## t4(i1+1,i2,i3)-u ## t4(i1+2,i2,i3))/(12.*dr(0)**2)
!  #defineMacro u ## Arst4(i1,i2,i3) (u ## st4(i1-2,i2,i3)-8.*u ## st4(i1-1,i2,i3)+8.*u ## st4(i1+1,i2,i3)-u ## st4(i1+2,i2,i3))/(12.*dr(0))
!  #defineMacro u ## Asst4(i1,i2,i3) (-u ## t4(i1,i2-2,i3)+16.*u ## t4(i1,i2-1,i3)-30.*u ## t4(i1,i2,i3)+16.*u ## t4(i1,i2+1,i3)-u ## t4(i1,i2+2,i3))/(12.*dr(1)**2)
!  #defineMacro u ## Artt4(i1,i2,i3) (u ## tt4(i1-2,i2,i3)-8.*u ## tt4(i1-1,i2,i3)+8.*u ## tt4(i1+1,i2,i3)-u ## tt4(i1+2,i2,i3))/(12.*dr(0))
!  #defineMacro u ## Astt4(i1,i2,i3) (u ## tt4(i1,i2-2,i3)-8.*u ## tt4(i1,i2-1,i3)+8.*u ## tt4(i1,i2+1,i3)-u ## tt4(i1,i2+2,i3))/(12.*dr(1))
!  #defineMacro u ## Attt4(i1,i2,i3) (u(i1,i2,i3-3)-8.*u(i1,i2,i3-2)+13.*u(i1,i2,i3-1)-13.*u(i1,i2,i3+1)+8.*u(i1,i2,i3+2)-u(i1,i2,i3+3))/(8.*dr(2)**3)
! #endMacro

! #beginMacro defineParameticDerivativesComponents1(u)
!  #defineMacro u ## Ar2(i1,i2,i3,m) (-u(i1-1,i2,i3,m)+u(i1+1,i2,i3,m))/(2.*dr(0))
!  #defineMacro u ## As2(i1,i2,i3,m) (-u(i1,i2-1,i3,m)+u(i1,i2+1,i3,m))/(2.*dr(1))
!  #defineMacro u ## At2(i1,i2,i3,m) (-u(i1,i2,i3-1,m)+u(i1,i2,i3+1,m))/(2.*dr(2))
!  #defineMacro u ## Arr2(i1,i2,i3,m) (u(i1-1,i2,i3,m)-2.*u(i1,i2,i3,m)+u(i1+1,i2,i3,m))/(dr(0)**2)
!  #defineMacro u ## Ars2(i1,i2,i3,m) (-u ## s2(i1-1,i2,i3,m)+u ## s2(i1+1,i2,i3,m))/(2.*dr(0))
!  #defineMacro u ## Ass2(i1,i2,i3,m) (u(i1,i2-1,i3,m)-2.*u(i1,i2,i3,m)+u(i1,i2+1,i3,m))/(dr(1)**2)
!  #defineMacro u ## Art2(i1,i2,i3,m) (-u ## t2(i1-1,i2,i3,m)+u ## t2(i1+1,i2,i3,m))/(2.*dr(0))
!  #defineMacro u ## Ast2(i1,i2,i3,m) (-u ## t2(i1,i2-1,i3,m)+u ## t2(i1,i2+1,i3,m))/(2.*dr(1))
!  #defineMacro u ## Att2(i1,i2,i3,m) (u(i1,i2,i3-1,m)-2.*u(i1,i2,i3,m)+u(i1,i2,i3+1,m))/(dr(2)**2)
!  #defineMacro u ## Arrr2(i1,i2,i3,m) (-u(i1-2,i2,i3,m)+2.*u(i1-1,i2,i3,m)-2.*u(i1+1,i2,i3,m)+u(i1+2,i2,i3,m))/(2.*dr(0)**3)
!  #defineMacro u ## Arrs2(i1,i2,i3,m) (u ## s2(i1-1,i2,i3,m)-2.*u ## s2(i1,i2,i3,m)+u ## s2(i1+1,i2,i3,m))/(dr(0)**2)
!  #defineMacro u ## Arss2(i1,i2,i3,m) (-u ## ss2(i1-1,i2,i3,m)+u ## ss2(i1+1,i2,i3,m))/(2.*dr(0))
!  #defineMacro u ## Asss2(i1,i2,i3,m) (-u(i1,i2-2,i3,m)+2.*u(i1,i2-1,i3,m)-2.*u(i1,i2+1,i3,m)+u(i1,i2+2,i3,m))/(2.*dr(1)**3)
!  #defineMacro u ## Arrt2(i1,i2,i3,m) (u ## t2(i1-1,i2,i3,m)-2.*u ## t2(i1,i2,i3,m)+u ## t2(i1+1,i2,i3,m))/(dr(0)**2)
!  #defineMacro u ## Arst2(i1,i2,i3,m) (-u ## st2(i1-1,i2,i3,m)+u ## st2(i1+1,i2,i3,m))/(2.*dr(0))
!  #defineMacro u ## Asst2(i1,i2,i3,m) (u ## t2(i1,i2-1,i3,m)-2.*u ## t2(i1,i2,i3,m)+u ## t2(i1,i2+1,i3,m))/(dr(1)**2)
!  #defineMacro u ## Artt2(i1,i2,i3,m) (-u ## tt2(i1-1,i2,i3,m)+u ## tt2(i1+1,i2,i3,m))/(2.*dr(0))
!  #defineMacro u ## Astt2(i1,i2,i3,m) (-u ## tt2(i1,i2-1,i3,m)+u ## tt2(i1,i2+1,i3,m))/(2.*dr(1))
!  #defineMacro u ## Attt2(i1,i2,i3,m) (-u(i1,i2,i3-2,m)+2.*u(i1,i2,i3-1,m)-2.*u(i1,i2,i3+1,m)+u(i1,i2,i3+2,m))/(2.*dr(2)**3)
!  #defineMacro u ## Ar4(i1,i2,i3,m) (u(i1-2,i2,i3,m)-8.*u(i1-1,i2,i3,m)+8.*u(i1+1,i2,i3,m)-u(i1+2,i2,i3,m))/(12.*dr(0))
!  #defineMacro u ## As4(i1,i2,i3,m) (u(i1,i2-2,i3,m)-8.*u(i1,i2-1,i3,m)+8.*u(i1,i2+1,i3,m)-u(i1,i2+2,i3,m))/(12.*dr(1))
!  #defineMacro u ## At4(i1,i2,i3,m) (u(i1,i2,i3-2,m)-8.*u(i1,i2,i3-1,m)+8.*u(i1,i2,i3+1,m)-u(i1,i2,i3+2,m))/(12.*dr(2))
!  #defineMacro u ## Arr4(i1,i2,i3,m) (-u(i1-2,i2,i3,m)+16.*u(i1-1,i2,i3,m)-30.*u(i1,i2,i3,m)+16.*u(i1+1,i2,i3,m)-u(i1+2,i2,i3,m))/(12.*dr(0)**2)
!  #defineMacro u ## Ars4(i1,i2,i3,m) (u ## s4(i1-2,i2,i3,m)-8.*u ## s4(i1-1,i2,i3,m)+8.*u ## s4(i1+1,i2,i3,m)-u ## s4(i1+2,i2,i3,m))/(12.*dr(0))
!  #defineMacro u ## Ass4(i1,i2,i3,m) (-u(i1,i2-2,i3,m)+16.*u(i1,i2-1,i3,m)-30.*u(i1,i2,i3,m)+16.*u(i1,i2+1,i3,m)-u(i1,i2+2,i3,m))/(12.*dr(1)**2)
!  #defineMacro u ## Art4(i1,i2,i3,m) (u ## t4(i1-2,i2,i3,m)-8.*u ## t4(i1-1,i2,i3,m)+8.*u ## t4(i1+1,i2,i3,m)-u ## t4(i1+2,i2,i3,m))/(12.*dr(0))
!  #defineMacro u ## Ast4(i1,i2,i3,m) (u ## t4(i1,i2-2,i3,m)-8.*u ## t4(i1,i2-1,i3,m)+8.*u ## t4(i1,i2+1,i3,m)-u ## t4(i1,i2+2,i3,m))/(12.*dr(1))
!  #defineMacro u ## Att4(i1,i2,i3,m) (-u(i1,i2,i3-2,m)+16.*u(i1,i2,i3-1,m)-30.*u(i1,i2,i3,m)+16.*u(i1,i2,i3+1,m)-u(i1,i2,i3+2,m))/(12.*dr(2)**2)
!  #defineMacro u ## Arrr4(i1,i2,i3,m) (u(i1-3,i2,i3,m)-8.*u(i1-2,i2,i3,m)+13.*u(i1-1,i2,i3,m)-13.*u(i1+1,i2,i3,m)+8.*u(i1+2,i2,i3,m)-u(i1+3,i2,i3,m))/(8.*dr(0)**3)
!  #defineMacro u ## Arrs4(i1,i2,i3,m) (-u ## s4(i1-2,i2,i3,m)+16.*u ## s4(i1-1,i2,i3,m)-30.*u ## s4(i1,i2,i3,m)+16.*u ## s4(i1+1,i2,i3,m)-u ## s4(i1+2,i2,i3,m))/(12.*dr(0)**2)
!  #defineMacro u ## Arss4(i1,i2,i3,m) (u ## ss4(i1-2,i2,i3,m)-8.*u ## ss4(i1-1,i2,i3,m)+8.*u ## ss4(i1+1,i2,i3,m)-u ## ss4(i1+2,i2,i3,m))/(12.*dr(0))
!  #defineMacro u ## Asss4(i1,i2,i3,m) (u(i1,i2-3,i3,m)-8.*u(i1,i2-2,i3,m)+13.*u(i1,i2-1,i3,m)-13.*u(i1,i2+1,i3,m)+8.*u(i1,i2+2,i3,m)-u(i1,i2+3,i3,m))/(8.*dr(1)**3)
!  #defineMacro u ## Arrt4(i1,i2,i3,m) (-u ## t4(i1-2,i2,i3,m)+16.*u ## t4(i1-1,i2,i3,m)-30.*u ## t4(i1,i2,i3,m)+16.*u ## t4(i1+1,i2,i3,m)-u ## t4(i1+2,i2,i3,m))/(12.*dr(0)**2)
!  #defineMacro u ## Arst4(i1,i2,i3,m) (u ## st4(i1-2,i2,i3,m)-8.*u ## st4(i1-1,i2,i3,m)+8.*u ## st4(i1+1,i2,i3,m)-u ## st4(i1+2,i2,i3,m))/(12.*dr(0))
!  #defineMacro u ## Asst4(i1,i2,i3,m) (-u ## t4(i1,i2-2,i3,m)+16.*u ## t4(i1,i2-1,i3,m)-30.*u ## t4(i1,i2,i3,m)+16.*u ## t4(i1,i2+1,i3,m)-u ## t4(i1,i2+2,i3,m))/(12.*dr(1)**2)
!  #defineMacro u ## Artt4(i1,i2,i3,m) (u ## tt4(i1-2,i2,i3,m)-8.*u ## tt4(i1-1,i2,i3,m)+8.*u ## tt4(i1+1,i2,i3,m)-u ## tt4(i1+2,i2,i3,m))/(12.*dr(0))
!  #defineMacro u ## Astt4(i1,i2,i3,m) (u ## tt4(i1,i2-2,i3,m)-8.*u ## tt4(i1,i2-1,i3,m)+8.*u ## tt4(i1,i2+1,i3,m)-u ## tt4(i1,i2+2,i3,m))/(12.*dr(1))
!  #defineMacro u ## Attt4(i1,i2,i3,m) (u(i1,i2,i3-3,m)-8.*u(i1,i2,i3-2,m)+13.*u(i1,i2,i3-1,m)-13.*u(i1,i2,i3+1,m)+8.*u(i1,i2,i3+2,m)-u(i1,i2,i3+3,m))/(8.*dr(2)**3)
! #endMacro


! ! =======================================================
! !  Macro to compute Third derivatives in 2 dimensions 
! !  OPTION : evalMetrics : evaluate the derivatives of the metrics
! !          (metrics need only be evaluated once when using discrete delta to get coeffs)
! ! =======================================================
! #beginMacro getThirdDerivatives2d(ORDER,GRIDTYPE,OPTION,i1,i2,i3)

! #If #GRIDTYPE eq "rectangular" 
! ! ---------- RECTANGULAR  ---------
! uxxx = uxxx22r(i1,i2,i3,0)
! uxxy = uxxy22r(i1,i2,i3,0)
! uxyy = uxyy22r(i1,i2,i3,0)
! uyyy = uyyy22r(i1,i2,i3,0)

! #Else
! ! ---------- START CURVILINEAR  ---------
! defineParameticDerivativesComponents1(u)
! #If #OPTION eq "evalMetrics"
! defineParameticDerivativesComponents0(rx)
! defineParameticDerivativesComponents0(ry)
! defineParameticDerivativesComponents0(sx)
! defineParameticDerivativesComponents0(sy)
! #End

! ! ---------- Parametric derivatives ---------
! ur     = uAr4(i1,i2,i3,0)
! urr    = uArr4(i1,i2,i3,0)
! urrr   = uArrr2(i1,i2,i3,0)
! us     = uAs4(i1,i2,i3,0)
! urs    = uArs4(i1,i2,i3,0)
! urrs   = uArrs2(i1,i2,i3,0)
! uss    = uAss4(i1,i2,i3,0)
! urss   = uArss2(i1,i2,i3,0)
! usss   = uAsss2(i1,i2,i3,0)
! #If #OPTION eq "evalMetrics"
! rxr    = rxAr4(i1,i2,i3)
! rxrr   = rxArr4(i1,i2,i3)
! rxs    = rxAs4(i1,i2,i3)
! rxrs   = rxArs4(i1,i2,i3)
! rxss   = rxAss4(i1,i2,i3)
! ryr    = ryAr4(i1,i2,i3)
! ryrr   = ryArr4(i1,i2,i3)
! rys    = ryAs4(i1,i2,i3)
! ryrs   = ryArs4(i1,i2,i3)
! ryss   = ryAss4(i1,i2,i3)
! sxr    = sxAr4(i1,i2,i3)
! sxrr   = sxArr4(i1,i2,i3)
! sxs    = sxAs4(i1,i2,i3)
! sxrs   = sxArs4(i1,i2,i3)
! sxss   = sxAss4(i1,i2,i3)
! syr    = syAr4(i1,i2,i3)
! syrr   = syArr4(i1,i2,i3)
! sys    = syAs4(i1,i2,i3)
! syrs   = syArs4(i1,i2,i3)
! syss   = syAss4(i1,i2,i3)

! ! ---------- Spatial derivatives of metrics rx, sx, ry, ... ---------
! rxi = rx(i1,i2,i3)
! ryi = ry(i1,i2,i3)
! sxi = sx(i1,i2,i3)
! syi = sy(i1,i2,i3)
! rxx = rxi*rxr+sxi*rxs
! rxxx = rxi**2*rxrr+2.*rxi*sxi*rxrs+sxi**2*rxss+rxx*rxr+sxx*rxs
! rxy = ryi*rxr+syi*rxs
! rxxy = ryi*rxi*rxrr+(rxi*syi+ryi*sxi)*rxrs+syi*sxi*rxss+rxy*rxr+sxy*rxs
! rxyy = ryi**2*rxrr+2.*ryi*syi*rxrs+syi**2*rxss+ryy*rxr+syy*rxs
! ryy = ryi*ryr+syi*rys
! ryyy = ryi**2*ryrr+2.*ryi*syi*ryrs+syi**2*ryss+ryy*ryr+syy*rys
! sxx = rxi*sxr+sxi*sxs
! sxxx = rxi**2*sxrr+2.*rxi*sxi*sxrs+sxi**2*sxss+rxx*sxr+sxx*sxs
! sxy = ryi*sxr+syi*sxs
! sxxy = ryi*rxi*sxrr+(rxi*syi+ryi*sxi)*sxrs+syi*sxi*sxss+rxy*sxr+sxy*sxs
! sxyy = ryi**2*sxrr+2.*ryi*syi*sxrs+syi**2*sxss+ryy*sxr+syy*sxs
! syy = ryi*syr+syi*sys
! syyy = ryi**2*syrr+2.*ryi*syi*syrs+syi**2*syss+ryy*syr+syy*sys
! #End
! ! ---- end OPTION eq evalMetrics ---

! ! ---------- Third spatial derivatives of u ---------
! uxxx = rxi**3*urrr+3.*rxi**2*sxi*urrs+3.*rxi*sxi**2*urss+sxi**3*usss+3.*rxi*rxx*urr+(3.*rxi*sxx+3.*rxx*sxi)*urs+3.*sxx*sxi*uss+rxxx*ur+sxxx*us
! uxxy = ryi*rxi**2*urrr+(rxi**2*syi+2.*rxi*ryi*sxi)*urrs+(2.*rxi*sxi*syi+ryi*sxi**2)*urss+syi*sxi**2*usss+(2.*rxi*rxy+rxx*ryi)*urr+(2.*rxi*sxy+rxx*syi+2.*rxy*sxi+ryi*sxx)*urs+(2.*sxi*sxy+sxx*syi)*uss+rxxy*ur+sxxy*us
! uxyy = ryi**2*rxi*urrr+(syi*ryi*rxi+ryi*(rxi*syi+ryi*sxi))*urrs+(ryi*syi*sxi+syi*(rxi*syi+ryi*sxi))*urss+syi**2*sxi*usss+(rxi*ryy+2.*rxy*ryi)*urr+(rxi*syy+2.*rxy*syi+2.*ryi*sxy+ryy*sxi)*urs+(sxi*syy+2.*sxy*syi)*uss+rxyy*ur+sxyy*us
! uyyy = ryi**3*urrr+3.*ryi**2*syi*urrs+3.*ryi*syi**2*urss+syi**3*usss+3.*ryi*ryy*urr+(3.*ryi*syy+3.*ryy*syi)*urs+3.*syy*syi*uss+ryyy*ur+syyy*us
! ! ---------- END CURVILINEAR  ---------
! #End
! #endMacro


! ! #beginMacro getThirdDerivatives2d(ORDER,GRIDTYPE,OPTION,i1,i2,i3)
! ! #endMacro

! #beginMacro getThirdDerivatives3d(ORDER,GRIDTYPE,OPTION,i1,i2,i3)
! #endMacro

! define macros to evaluate derivatives for the 6th order method (from maple/makeGetDerivativesMacros.maple)
!! turned off May 4, 2023
!! #Include "../maple/defineGetSixthDerivativesMacros.h"


! From bcOptSmFOS.bf
! DataBase *pdb = &parameters.dbase;
! double precision pdb  ! pointer to data base
! ====================================================================
! Look up an integer parameter from the data base
! ====================================================================

! ====================================================================
! Look up a real parameter from the data base
! ====================================================================



! General begin loops macro











! ----- define extrapolation formulae ------











! ************************************************************************************************
!  This macro is used for looping over the faces of a grid to assign booundary conditions
!
! extra: extra points to assign
!          Case 1: extra=numberOfGhostPoints -- for assigning extended boundaries
!          Case 2: extra=-1 -- for assigning ghost points but not including extended boundaries
! numberOfGhostPoints : number of ghost points (1 for 2nd order, 2 for fourth-order ...)
!
!
! Output:
!  n1a,n1b,n2a,n2b,n3a,n3b : from gridIndexRange
!  nn1a,nn1b,nn2a,nn2b,nn3a,nn3b : includes "extra" points
! 
! ***********************************************************************************************







! =================================================================================
!   Assign values in the corners in 2D (see bcMaxwellCorners.bf)
!
!  Set the normal component of the solution on the extended boundaries (points N in figure)
!  Set the corner points "C" 
!              |
!              X
!              |
!        N--N--X--X----
!              |
!        C  C  N
!              |
!        C  C  N
!
! ORDER: 2 or 4
! GRIDTYPE: rectangular, curvilinear
! FORCING: none, twilightZone
! =================================================================================







! =========================================================================
! Compute the normal on a curvilinear grid.
!
! Assumes is=1-2*side is defined. 
! =========================================================================

! ========================================================================================
!  Assign ghost points outside corners
! ========================================================================================




! ===================================================================================================
! Macro: Extrapolate Ghost Points 
! ORDER : 2,4,6,8
! ===================================================================================================


! ============================================================================================
! Macro: evaluate the solution on the boundary for Dirichlet BCs
! ============================================================================================

! ===================================================================================================
! Macro: Assign boundary and ghost points on Dirichlet boundaries 
! ORDER : 2,4,6,8
! ===================================================================================================


! ============================================================================================
! Macro: evaluate the solution on the boundary for Dirichlet BCs
!    Solving
!        u_tt = c^2 * Delta( u ) + f(x,y,z,t)
!        u = g
! For TZ at order=2:
!    ff = ue_tt - c^2*Delta(ue) 
!    gtt = g_tt = uett
!  order=4:
!    fLap 
!    ftt
!    gtttt
! ============================================================================================


! ===================================================================================================
! Macro: Assign the boundary values on Dirichlet boundaries
!
!  NOTE: DIM AND GRIDTYPE NOT CURRENTLY USED
!  FORCE : USEFORCING or NOFORCING 
! ===================================================================================================

! ===================================================================================================
! Macro: Assign the boundary values on "exact" boundaries
!
! **NOT USED ANY MORE**
! Exact boundaries are set in applyBoundaryConditions.bC
!
! ===================================================================================================
! #beginMacro assignExactBoundary()

!   ! assign extram points in the tangential directions
!   extram = numGhost 
!   getLoopBounds(side,axis,extram, m1a,m1b,m2a,m2b,m3a,m3b)

!   ff=0.
!   beginLoops(m1a,m1b,m2a,m2b,m3a,m3b)
!     if( mask(i1,i2,i3).ne.0 )then
     
!       do ghost=0,numGhost
!         j1=i1-is1*ghost
!         j2=i2-is2*ghost
!         j3=i3-is3*ghost

!         getDirichletForcing(ff,j1,j2,j3)

!         u(j1,j2,j3,uc) = ff

!       end do
  
!     end if ! mask .ne. 0
!   endLoops3d()

! #endMacro



! ===================================================================================================
! Macro: get loop bounds over boundaries with extram points in tangential directions 
! ===================================================================================================


! ===================================================================================================
! Macro: Assign ghost points on Dirichlet boundaries using: 
!          *** COMPATIBILITY BOUNDARY CONDITIONS ****
! ORDER : 2,4,6,8
! FORCING : noForcing, forcing
! ===================================================================================================

! ===================================================================================================
! Macro: Get the TZ solution in 2d or 3d 
! ===================================================================================================

! ===================================================================================================
! Macro: Evaluate the forcing for the CBCs
! ===================================================================================================

! ===================================================================================================
! Macro: Assign ghost points on Dirichlet boundaries using: 
!          *** COMPATIBILITY BOUNDARY CONDITIONS ****
! ORDER : 6 
! FORCING : noForcing, forcing
! ===================================================================================================




! ------------------------------------------------------------------------------------
!  Macro: evaluate the RHS to the Neumann BC
! ------------------------------------------------------------------------------------

! ===================================================================================================
! Macro: Assign ghost points on Neumann boundaries 
! ORDER : 2,4,6,8
! ===================================================================================================


! ===================================================================================================
! Macro: Assign ghost points on Neumann boundaries
!          *** COMPATIBILITY BOUNDARY CONDITIONS **** 
! ORDER : 2,4,6,8
! ===================================================================================================
! *** OLD VERSION ****


! ============================================================================================
! Macro: evaluate the forcings for Neumann CBCs
!    Solving
!        u_tt = c^2 * Delta( u ) + f(x,y,z,t)
!        u.n = g
! For TZ at order=2:
!    gg = g
!  order=4:
!    gg,
!    nDotGradF = n.grad( f ), f = ue_tt - c^2*lap(ue)
!    gtt
! ============================================================================================

!==========================================================================
!  Check the coefficients in the ghost points of the residual equations
! using the discrete delta approach
!==========================================================================

! ===================================================================================================
! Macro: Assign ghost points on Dirichlet boundaries using: 
!          *** COMPATIBILITY BOUNDARY CONDITIONS ****
! ORDER : 2,4,6,8
! FORCING : noForcing, forcing
! ===================================================================================================


! ========================================================================================
!  Apply symmetry conditions in corner ghost for Cartesian grids 
! ========================================================================================





! ===================================================================================
! Utility macro to call different versions of assignDirichletGhostCompatibility
! for a given ORDER
! ===================================================================================

! ===================================================================================
! Utility macro to call different versions of assignNeumannGhostCompatibility
! for a given ORDER
! ===================================================================================



! Argument list

! **********************************************************************************
! Macro BC_WAVE:
!  NAME: name of the subroutine
!  DIM : 2 or 3
!  ORDER : 2 ,4, 6 or 8
! **********************************************************************************

! --- Macro to build the file for each dimension and order ---

! --- construct the different files ----



! buildFile(bcOptWave2dOrder4,2,6)
! buildFile(bcOptWave3dOrder4,3,6)




subroutine bcOptWave( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr ) 
! ===================================================================================
!  Boundary conditions for CgWave
!
!  gridType : 0=rectangular, 1=curvilinear
!
! The forcing for the boundary conditions can be accessed using the statement function:
!         bcf(side,axis,i1,i2,i3,m)
! which is defined below. 
! ===================================================================================

  implicit none

  integer nd, nd1a,nd1b,nd2a,nd2b,nd3a,nd3b, ndb, ierr

  real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)
  integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
  real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
  real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1)

  real uTemp(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
  real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)

  integer gridIndexRange(0:1,0:2),boundaryCondition(0:1,0:2), dimRange(0:1,0:2), isPeriodic(0:*)
  real frequencyArray(0:*)

  double precision pdb  ! pointer to data base


  integer ipar(0:*)
  real rpar(0:*)

  integer orderOfAccuracy

  ! extract parameters we need: 
  orderOfAccuracy  = ipar( 4)

  if( nd.eq.2 )then
    if( orderOfAccuracy.eq.2 )then
      call bcOptWave2dOrder2( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    elseif( orderOfAccuracy.eq.4 )then
      call bcOptWave2dOrder4( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    else
      stop 6666
    end if
  else
    if( orderOfAccuracy.eq.2 )then
      call bcOptWave2dOrder2( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    elseif( orderOfAccuracy.eq.4 )then
      call bcOptWave2dOrder4( nd,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,gridIndexRange,dimRange,isPeriodic,u,mask,rsxy,xy,uTemp,v,boundaryCondition,frequencyArray,pdb,ipar,rpar,ierr )
    else
      stop 7777
    end if    

  end if

  return
  end

