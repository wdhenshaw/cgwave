! ****** File written by makeGetDerivativesMacros.maple  ******

#beginMacro defineParameticDerivativesComponents0(u)
 #defineMacro u ## r2(i1,i2,i3) (-u(i1-1,i2,i3)+u(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## s2(i1,i2,i3) (-u(i1,i2-1,i3)+u(i1,i2+1,i3))/(2.*dr(1))
 #defineMacro u ## t2(i1,i2,i3) (-u(i1,i2,i3-1)+u(i1,i2,i3+1))/(2.*dr(2))
 #defineMacro u ## rr2(i1,i2,i3) (u(i1-1,i2,i3)-2.*u(i1,i2,i3)+u(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rs2(i1,i2,i3) (-u ## s2(i1-1,i2,i3)+u ## s2(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## ss2(i1,i2,i3) (u(i1,i2-1,i3)-2.*u(i1,i2,i3)+u(i1,i2+1,i3))/(dr(1)**2)
 #defineMacro u ## rt2(i1,i2,i3) (-u ## t2(i1-1,i2,i3)+u ## t2(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## st2(i1,i2,i3) (-u ## t2(i1,i2-1,i3)+u ## t2(i1,i2+1,i3))/(2.*dr(1))
 #defineMacro u ## tt2(i1,i2,i3) (u(i1,i2,i3-1)-2.*u(i1,i2,i3)+u(i1,i2,i3+1))/(dr(2)**2)
 #defineMacro u ## rrr2(i1,i2,i3) (-u(i1-2,i2,i3)+2.*u(i1-1,i2,i3)-2.*u(i1+1,i2,i3)+u(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrs2(i1,i2,i3) (u ## s2(i1-1,i2,i3)-2.*u ## s2(i1,i2,i3)+u ## s2(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rss2(i1,i2,i3) (-u ## ss2(i1-1,i2,i3)+u ## ss2(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## sss2(i1,i2,i3) (-u(i1,i2-2,i3)+2.*u(i1,i2-1,i3)-2.*u(i1,i2+1,i3)+u(i1,i2+2,i3))/(2.*dr(1)**3)
 #defineMacro u ## rrt2(i1,i2,i3) (u ## t2(i1-1,i2,i3)-2.*u ## t2(i1,i2,i3)+u ## t2(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rst2(i1,i2,i3) (-u ## st2(i1-1,i2,i3)+u ## st2(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## sst2(i1,i2,i3) (u ## t2(i1,i2-1,i3)-2.*u ## t2(i1,i2,i3)+u ## t2(i1,i2+1,i3))/(dr(1)**2)
 #defineMacro u ## rtt2(i1,i2,i3) (-u ## tt2(i1-1,i2,i3)+u ## tt2(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## stt2(i1,i2,i3) (-u ## tt2(i1,i2-1,i3)+u ## tt2(i1,i2+1,i3))/(2.*dr(1))
 #defineMacro u ## ttt2(i1,i2,i3) (-u(i1,i2,i3-2)+2.*u(i1,i2,i3-1)-2.*u(i1,i2,i3+1)+u(i1,i2,i3+2))/(2.*dr(2)**3)
 #defineMacro u ## r4(i1,i2,i3) (u(i1-2,i2,i3)-8.*u(i1-1,i2,i3)+8.*u(i1+1,i2,i3)-u(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## s4(i1,i2,i3) (u(i1,i2-2,i3)-8.*u(i1,i2-1,i3)+8.*u(i1,i2+1,i3)-u(i1,i2+2,i3))/(12.*dr(1))
 #defineMacro u ## t4(i1,i2,i3) (u(i1,i2,i3-2)-8.*u(i1,i2,i3-1)+8.*u(i1,i2,i3+1)-u(i1,i2,i3+2))/(12.*dr(2))
 #defineMacro u ## rr4(i1,i2,i3) (-u(i1-2,i2,i3)+16.*u(i1-1,i2,i3)-30.*u(i1,i2,i3)+16.*u(i1+1,i2,i3)-u(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## rs4(i1,i2,i3) (u ## s4(i1-2,i2,i3)-8.*u ## s4(i1-1,i2,i3)+8.*u ## s4(i1+1,i2,i3)-u ## s4(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## ss4(i1,i2,i3) (-u(i1,i2-2,i3)+16.*u(i1,i2-1,i3)-30.*u(i1,i2,i3)+16.*u(i1,i2+1,i3)-u(i1,i2+2,i3))/(12.*dr(1)**2)
 #defineMacro u ## rt4(i1,i2,i3) (u ## t4(i1-2,i2,i3)-8.*u ## t4(i1-1,i2,i3)+8.*u ## t4(i1+1,i2,i3)-u ## t4(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## st4(i1,i2,i3) (u ## t4(i1,i2-2,i3)-8.*u ## t4(i1,i2-1,i3)+8.*u ## t4(i1,i2+1,i3)-u ## t4(i1,i2+2,i3))/(12.*dr(1))
 #defineMacro u ## tt4(i1,i2,i3) (-u(i1,i2,i3-2)+16.*u(i1,i2,i3-1)-30.*u(i1,i2,i3)+16.*u(i1,i2,i3+1)-u(i1,i2,i3+2))/(12.*dr(2)**2)
 #defineMacro u ## rrr4(i1,i2,i3) (u(i1-3,i2,i3)-8.*u(i1-2,i2,i3)+13.*u(i1-1,i2,i3)-13.*u(i1+1,i2,i3)+8.*u(i1+2,i2,i3)-u(i1+3,i2,i3))/(8.*dr(0)**3)
 #defineMacro u ## rrs4(i1,i2,i3) (-u ## s4(i1-2,i2,i3)+16.*u ## s4(i1-1,i2,i3)-30.*u ## s4(i1,i2,i3)+16.*u ## s4(i1+1,i2,i3)-u ## s4(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## rss4(i1,i2,i3) (u ## ss4(i1-2,i2,i3)-8.*u ## ss4(i1-1,i2,i3)+8.*u ## ss4(i1+1,i2,i3)-u ## ss4(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## sss4(i1,i2,i3) (u(i1,i2-3,i3)-8.*u(i1,i2-2,i3)+13.*u(i1,i2-1,i3)-13.*u(i1,i2+1,i3)+8.*u(i1,i2+2,i3)-u(i1,i2+3,i3))/(8.*dr(1)**3)
 #defineMacro u ## rrt4(i1,i2,i3) (-u ## t4(i1-2,i2,i3)+16.*u ## t4(i1-1,i2,i3)-30.*u ## t4(i1,i2,i3)+16.*u ## t4(i1+1,i2,i3)-u ## t4(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## rst4(i1,i2,i3) (u ## st4(i1-2,i2,i3)-8.*u ## st4(i1-1,i2,i3)+8.*u ## st4(i1+1,i2,i3)-u ## st4(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## sst4(i1,i2,i3) (-u ## t4(i1,i2-2,i3)+16.*u ## t4(i1,i2-1,i3)-30.*u ## t4(i1,i2,i3)+16.*u ## t4(i1,i2+1,i3)-u ## t4(i1,i2+2,i3))/(12.*dr(1)**2)
 #defineMacro u ## rtt4(i1,i2,i3) (u ## tt4(i1-2,i2,i3)-8.*u ## tt4(i1-1,i2,i3)+8.*u ## tt4(i1+1,i2,i3)-u ## tt4(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## stt4(i1,i2,i3) (u ## tt4(i1,i2-2,i3)-8.*u ## tt4(i1,i2-1,i3)+8.*u ## tt4(i1,i2+1,i3)-u ## tt4(i1,i2+2,i3))/(12.*dr(1))
 #defineMacro u ## ttt4(i1,i2,i3) (u(i1,i2,i3-3)-8.*u(i1,i2,i3-2)+13.*u(i1,i2,i3-1)-13.*u(i1,i2,i3+1)+8.*u(i1,i2,i3+2)-u(i1,i2,i3+3))/(8.*dr(2)**3)
#endMacro

#beginMacro defineParameticDerivativesComponents1(u)
 #defineMacro u ## r2(i1,i2,i3,m) (-u(i1-1,i2,i3,m)+u(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## s2(i1,i2,i3,m) (-u(i1,i2-1,i3,m)+u(i1,i2+1,i3,m))/(2.*dr(1))
 #defineMacro u ## t2(i1,i2,i3,m) (-u(i1,i2,i3-1,m)+u(i1,i2,i3+1,m))/(2.*dr(2))
 #defineMacro u ## rr2(i1,i2,i3,m) (u(i1-1,i2,i3,m)-2.*u(i1,i2,i3,m)+u(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rs2(i1,i2,i3,m) (-u ## s2(i1-1,i2,i3,m)+u ## s2(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## ss2(i1,i2,i3,m) (u(i1,i2-1,i3,m)-2.*u(i1,i2,i3,m)+u(i1,i2+1,i3,m))/(dr(1)**2)
 #defineMacro u ## rt2(i1,i2,i3,m) (-u ## t2(i1-1,i2,i3,m)+u ## t2(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## st2(i1,i2,i3,m) (-u ## t2(i1,i2-1,i3,m)+u ## t2(i1,i2+1,i3,m))/(2.*dr(1))
 #defineMacro u ## tt2(i1,i2,i3,m) (u(i1,i2,i3-1,m)-2.*u(i1,i2,i3,m)+u(i1,i2,i3+1,m))/(dr(2)**2)
 #defineMacro u ## rrr2(i1,i2,i3,m) (-u(i1-2,i2,i3,m)+2.*u(i1-1,i2,i3,m)-2.*u(i1+1,i2,i3,m)+u(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrs2(i1,i2,i3,m) (u ## s2(i1-1,i2,i3,m)-2.*u ## s2(i1,i2,i3,m)+u ## s2(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rss2(i1,i2,i3,m) (-u ## ss2(i1-1,i2,i3,m)+u ## ss2(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## sss2(i1,i2,i3,m) (-u(i1,i2-2,i3,m)+2.*u(i1,i2-1,i3,m)-2.*u(i1,i2+1,i3,m)+u(i1,i2+2,i3,m))/(2.*dr(1)**3)
 #defineMacro u ## rrt2(i1,i2,i3,m) (u ## t2(i1-1,i2,i3,m)-2.*u ## t2(i1,i2,i3,m)+u ## t2(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rst2(i1,i2,i3,m) (-u ## st2(i1-1,i2,i3,m)+u ## st2(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## sst2(i1,i2,i3,m) (u ## t2(i1,i2-1,i3,m)-2.*u ## t2(i1,i2,i3,m)+u ## t2(i1,i2+1,i3,m))/(dr(1)**2)
 #defineMacro u ## rtt2(i1,i2,i3,m) (-u ## tt2(i1-1,i2,i3,m)+u ## tt2(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## stt2(i1,i2,i3,m) (-u ## tt2(i1,i2-1,i3,m)+u ## tt2(i1,i2+1,i3,m))/(2.*dr(1))
 #defineMacro u ## ttt2(i1,i2,i3,m) (-u(i1,i2,i3-2,m)+2.*u(i1,i2,i3-1,m)-2.*u(i1,i2,i3+1,m)+u(i1,i2,i3+2,m))/(2.*dr(2)**3)
 #defineMacro u ## r4(i1,i2,i3,m) (u(i1-2,i2,i3,m)-8.*u(i1-1,i2,i3,m)+8.*u(i1+1,i2,i3,m)-u(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## s4(i1,i2,i3,m) (u(i1,i2-2,i3,m)-8.*u(i1,i2-1,i3,m)+8.*u(i1,i2+1,i3,m)-u(i1,i2+2,i3,m))/(12.*dr(1))
 #defineMacro u ## t4(i1,i2,i3,m) (u(i1,i2,i3-2,m)-8.*u(i1,i2,i3-1,m)+8.*u(i1,i2,i3+1,m)-u(i1,i2,i3+2,m))/(12.*dr(2))
 #defineMacro u ## rr4(i1,i2,i3,m) (-u(i1-2,i2,i3,m)+16.*u(i1-1,i2,i3,m)-30.*u(i1,i2,i3,m)+16.*u(i1+1,i2,i3,m)-u(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## rs4(i1,i2,i3,m) (u ## s4(i1-2,i2,i3,m)-8.*u ## s4(i1-1,i2,i3,m)+8.*u ## s4(i1+1,i2,i3,m)-u ## s4(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## ss4(i1,i2,i3,m) (-u(i1,i2-2,i3,m)+16.*u(i1,i2-1,i3,m)-30.*u(i1,i2,i3,m)+16.*u(i1,i2+1,i3,m)-u(i1,i2+2,i3,m))/(12.*dr(1)**2)
 #defineMacro u ## rt4(i1,i2,i3,m) (u ## t4(i1-2,i2,i3,m)-8.*u ## t4(i1-1,i2,i3,m)+8.*u ## t4(i1+1,i2,i3,m)-u ## t4(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## st4(i1,i2,i3,m) (u ## t4(i1,i2-2,i3,m)-8.*u ## t4(i1,i2-1,i3,m)+8.*u ## t4(i1,i2+1,i3,m)-u ## t4(i1,i2+2,i3,m))/(12.*dr(1))
 #defineMacro u ## tt4(i1,i2,i3,m) (-u(i1,i2,i3-2,m)+16.*u(i1,i2,i3-1,m)-30.*u(i1,i2,i3,m)+16.*u(i1,i2,i3+1,m)-u(i1,i2,i3+2,m))/(12.*dr(2)**2)
 #defineMacro u ## rrr4(i1,i2,i3,m) (u(i1-3,i2,i3,m)-8.*u(i1-2,i2,i3,m)+13.*u(i1-1,i2,i3,m)-13.*u(i1+1,i2,i3,m)+8.*u(i1+2,i2,i3,m)-u(i1+3,i2,i3,m))/(8.*dr(0)**3)
 #defineMacro u ## rrs4(i1,i2,i3,m) (-u ## s4(i1-2,i2,i3,m)+16.*u ## s4(i1-1,i2,i3,m)-30.*u ## s4(i1,i2,i3,m)+16.*u ## s4(i1+1,i2,i3,m)-u ## s4(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## rss4(i1,i2,i3,m) (u ## ss4(i1-2,i2,i3,m)-8.*u ## ss4(i1-1,i2,i3,m)+8.*u ## ss4(i1+1,i2,i3,m)-u ## ss4(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## sss4(i1,i2,i3,m) (u(i1,i2-3,i3,m)-8.*u(i1,i2-2,i3,m)+13.*u(i1,i2-1,i3,m)-13.*u(i1,i2+1,i3,m)+8.*u(i1,i2+2,i3,m)-u(i1,i2+3,i3,m))/(8.*dr(1)**3)
 #defineMacro u ## rrt4(i1,i2,i3,m) (-u ## t4(i1-2,i2,i3,m)+16.*u ## t4(i1-1,i2,i3,m)-30.*u ## t4(i1,i2,i3,m)+16.*u ## t4(i1+1,i2,i3,m)-u ## t4(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## rst4(i1,i2,i3,m) (u ## st4(i1-2,i2,i3,m)-8.*u ## st4(i1-1,i2,i3,m)+8.*u ## st4(i1+1,i2,i3,m)-u ## st4(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## sst4(i1,i2,i3,m) (-u ## t4(i1,i2-2,i3,m)+16.*u ## t4(i1,i2-1,i3,m)-30.*u ## t4(i1,i2,i3,m)+16.*u ## t4(i1,i2+1,i3,m)-u ## t4(i1,i2+2,i3,m))/(12.*dr(1)**2)
 #defineMacro u ## rtt4(i1,i2,i3,m) (u ## tt4(i1-2,i2,i3,m)-8.*u ## tt4(i1-1,i2,i3,m)+8.*u ## tt4(i1+1,i2,i3,m)-u ## tt4(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## stt4(i1,i2,i3,m) (u ## tt4(i1,i2-2,i3,m)-8.*u ## tt4(i1,i2-1,i3,m)+8.*u ## tt4(i1,i2+1,i3,m)-u ## tt4(i1,i2+2,i3,m))/(12.*dr(1))
 #defineMacro u ## ttt4(i1,i2,i3,m) (u(i1,i2,i3-3,m)-8.*u(i1,i2,i3-2,m)+13.*u(i1,i2,i3-1,m)-13.*u(i1,i2,i3+1,m)+8.*u(i1,i2,i3+2,m)-u(i1,i2,i3+3,m))/(8.*dr(2)**3)
#endMacro

! =======================================================
!  Macro to compute Third derivatives in 2 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================
#beginMacro getThirdDerivatives2d(ORDER,GRIDTYPE,OPTION,i1,i2,i3)

#If #GRIDTYPE eq "rectangular" 
! ---------- RECTANGULAR  ---------
uxxx = uxxx22r(i1,i2,i3,0)
uxxy = uxxy22r(i1,i2,i3,0)
uxyy = uxyy22r(i1,i2,i3,0)
uyyy = uyyy22r(i1,i2,i3,0)

#Else
! ---------- START CURVILINEAR  ---------
defineParameticDerivativesComponents1(u)
#If #OPTION eq "evalMetrics"
defineParameticDerivativesComponents0(rx)
defineParameticDerivativesComponents0(ry)
defineParameticDerivativesComponents0(sx)
defineParameticDerivativesComponents0(sy)
#End

! ---------- Parametric derivatives ---------
ur     = ur4(i1,i2,i3,0)
urr    = urr4(i1,i2,i3,0)
urrr   = urrr2(i1,i2,i3,0)
us     = us4(i1,i2,i3,0)
urs    = urs4(i1,i2,i3,0)
urrs   = urrs2(i1,i2,i3,0)
uss    = uss4(i1,i2,i3,0)
urss   = urss2(i1,i2,i3,0)
usss   = usss2(i1,i2,i3,0)
#If #OPTION eq "evalMetrics"
rxr    = rxr4(i1,i2,i3)
rxrr   = rxrr4(i1,i2,i3)
rxs    = rxs4(i1,i2,i3)
rxrs   = rxrs4(i1,i2,i3)
rxss   = rxss4(i1,i2,i3)
ryr    = ryr4(i1,i2,i3)
ryrr   = ryrr4(i1,i2,i3)
rys    = rys4(i1,i2,i3)
ryrs   = ryrs4(i1,i2,i3)
ryss   = ryss4(i1,i2,i3)
sxr    = sxr4(i1,i2,i3)
sxrr   = sxrr4(i1,i2,i3)
sxs    = sxs4(i1,i2,i3)
sxrs   = sxrs4(i1,i2,i3)
sxss   = sxss4(i1,i2,i3)
syr    = syr4(i1,i2,i3)
syrr   = syrr4(i1,i2,i3)
sys    = sys4(i1,i2,i3)
syrs   = syrs4(i1,i2,i3)
syss   = syss4(i1,i2,i3)

! ---------- Spatial derivatives of metrics rx, sx, ry, ... ---------
rxi = rx(i1,i2,i3)
ryi = ry(i1,i2,i3)
sxi = sx(i1,i2,i3)
syi = sy(i1,i2,i3)
rxx = rxi*rxr+sxi*rxs
rxxx = rxi**2*rxrr+2.*rxi*sxi*rxrs+sxi**2*rxss+rxx*rxr+sxx*rxs
rxy = ryi*rxr+syi*rxs
rxxy = ryi*rxi*rxrr+(rxi*syi+ryi*sxi)*rxrs+syi*sxi*rxss+rxy*rxr+sxy*rxs
rxyy = ryi**2*rxrr+2.*ryi*syi*rxrs+syi**2*rxss+ryy*rxr+syy*rxs
ryy = ryi*ryr+syi*rys
ryyy = ryi**2*ryrr+2.*ryi*syi*ryrs+syi**2*ryss+ryy*ryr+syy*rys
sxx = rxi*sxr+sxi*sxs
sxxx = rxi**2*sxrr+2.*rxi*sxi*sxrs+sxi**2*sxss+rxx*sxr+sxx*sxs
sxy = ryi*sxr+syi*sxs
sxxy = ryi*rxi*sxrr+(rxi*syi+ryi*sxi)*sxrs+syi*sxi*sxss+rxy*sxr+sxy*sxs
sxyy = ryi**2*sxrr+2.*ryi*syi*sxrs+syi**2*sxss+ryy*sxr+syy*sxs
syy = ryi*syr+syi*sys
syyy = ryi**2*syrr+2.*ryi*syi*syrs+syi**2*syss+ryy*syr+syy*sys
#End
! ---- end OPTION eq evalMetrics ---

! ---------- Third spatial derivatives of u ---------
uxxx = rxi**3*urrr+3.*rxi**2*sxi*urrs+3.*rxi*sxi**2*urss+sxi**3*usss+3.*rxi*rxx*urr+(3.*rxi*sxx+3.*rxx*sxi)*urs+3.*sxx*sxi*uss+rxxx*ur+sxxx*us
uxxy = ryi*rxi**2*urrr+(rxi**2*syi+2.*rxi*ryi*sxi)*urrs+(2.*rxi*sxi*syi+ryi*sxi**2)*urss+syi*sxi**2*usss+(2.*rxi*rxy+rxx*ryi)*urr+(2.*rxi*sxy+rxx*syi+2.*rxy*sxi+ryi*sxx)*urs+(2.*sxi*sxy+sxx*syi)*uss+rxxy*ur+sxxy*us
uxyy = ryi**2*rxi*urrr+(syi*ryi*rxi+ryi*(rxi*syi+ryi*sxi))*urrs+(ryi*syi*sxi+syi*(rxi*syi+ryi*sxi))*urss+syi**2*sxi*usss+(rxi*ryy+2.*rxy*ryi)*urr+(rxi*syy+2.*rxy*syi+2.*ryi*sxy+ryy*sxi)*urs+(sxi*syy+2.*sxy*syi)*uss+rxyy*ur+sxyy*us
uyyy = ryi**3*urrr+3.*ryi**2*syi*urrs+3.*ryi*syi**2*urss+syi**3*usss+3.*ryi*ryy*urr+(3.*ryi*syy+3.*ryy*syi)*urs+3.*syy*syi*uss+ryyy*ur+syyy*us
! ---------- END CURVILINEAR  ---------
#End
#endMacro

! =======================================================
!  Macro to compute Third derivatives in 3 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================
#beginMacro getThirdDerivatives3d(ORDER,GRIDTYPE,OPTION,i1,i2,i3)

#If #GRIDTYPE eq "rectangular" 
! ---------- RECTANGULAR  ---------
uxxx = uxxx23r(i1,i2,i3,0)
uxxy = uxxy23r(i1,i2,i3,0)
uxxz = uxxz23r(i1,i2,i3,0)
uxyy = uxyy23r(i1,i2,i3,0)
uxzz = uxzz23r(i1,i2,i3,0)
uyyy = uyyy23r(i1,i2,i3,0)
uyyz = uyyz23r(i1,i2,i3,0)
uyzz = uyzz23r(i1,i2,i3,0)
uzzz = uzzz23r(i1,i2,i3,0)

#Else
! ---------- START CURVILINEAR  ---------
defineParameticDerivativesComponents1(u)
#If #OPTION eq "evalMetrics"
defineParameticDerivativesComponents0(rx)
defineParameticDerivativesComponents0(ry)
defineParameticDerivativesComponents0(sx)
defineParameticDerivativesComponents0(sy)
defineParameticDerivativesComponents0(rz)
defineParameticDerivativesComponents0(sz)
defineParameticDerivativesComponents0(tx)
defineParameticDerivativesComponents0(ty)
defineParameticDerivativesComponents0(tz)
#End

! ---------- Parametric derivatives ---------
ur     = ur4(i1,i2,i3,0)
urr    = urr4(i1,i2,i3,0)
urrr   = urrr2(i1,i2,i3,0)
us     = us4(i1,i2,i3,0)
urs    = urs4(i1,i2,i3,0)
urrs   = urrs2(i1,i2,i3,0)
uss    = uss4(i1,i2,i3,0)
urss   = urss2(i1,i2,i3,0)
usss   = usss2(i1,i2,i3,0)
ut     = ut4(i1,i2,i3,0)
urt    = urt4(i1,i2,i3,0)
urrt   = urrt2(i1,i2,i3,0)
ust    = ust4(i1,i2,i3,0)
urst   = urst2(i1,i2,i3,0)
usst   = usst2(i1,i2,i3,0)
utt    = utt4(i1,i2,i3,0)
urtt   = urtt2(i1,i2,i3,0)
ustt   = ustt2(i1,i2,i3,0)
uttt   = uttt2(i1,i2,i3,0)
#If #OPTION eq "evalMetrics"
rxr    = rxr4(i1,i2,i3)
rxrr   = rxrr4(i1,i2,i3)
rxs    = rxs4(i1,i2,i3)
rxrs   = rxrs4(i1,i2,i3)
rxss   = rxss4(i1,i2,i3)
rxt    = rxt4(i1,i2,i3)
rxrt   = rxrt4(i1,i2,i3)
rxst   = rxst4(i1,i2,i3)
rxtt   = rxtt4(i1,i2,i3)
ryr    = ryr4(i1,i2,i3)
ryrr   = ryrr4(i1,i2,i3)
rys    = rys4(i1,i2,i3)
ryrs   = ryrs4(i1,i2,i3)
ryss   = ryss4(i1,i2,i3)
ryt    = ryt4(i1,i2,i3)
ryrt   = ryrt4(i1,i2,i3)
ryst   = ryst4(i1,i2,i3)
rytt   = rytt4(i1,i2,i3)
sxr    = sxr4(i1,i2,i3)
sxrr   = sxrr4(i1,i2,i3)
sxs    = sxs4(i1,i2,i3)
sxrs   = sxrs4(i1,i2,i3)
sxss   = sxss4(i1,i2,i3)
sxt    = sxt4(i1,i2,i3)
sxrt   = sxrt4(i1,i2,i3)
sxst   = sxst4(i1,i2,i3)
sxtt   = sxtt4(i1,i2,i3)
syr    = syr4(i1,i2,i3)
syrr   = syrr4(i1,i2,i3)
sys    = sys4(i1,i2,i3)
syrs   = syrs4(i1,i2,i3)
syss   = syss4(i1,i2,i3)
syt    = syt4(i1,i2,i3)
syrt   = syrt4(i1,i2,i3)
syst   = syst4(i1,i2,i3)
sytt   = sytt4(i1,i2,i3)
rzr    = rzr4(i1,i2,i3)
rzrr   = rzrr4(i1,i2,i3)
rzs    = rzs4(i1,i2,i3)
rzrs   = rzrs4(i1,i2,i3)
rzss   = rzss4(i1,i2,i3)
rzt    = rzt4(i1,i2,i3)
rzrt   = rzrt4(i1,i2,i3)
rzst   = rzst4(i1,i2,i3)
rztt   = rztt4(i1,i2,i3)
szr    = szr4(i1,i2,i3)
szrr   = szrr4(i1,i2,i3)
szs    = szs4(i1,i2,i3)
szrs   = szrs4(i1,i2,i3)
szss   = szss4(i1,i2,i3)
szt    = szt4(i1,i2,i3)
szrt   = szrt4(i1,i2,i3)
szst   = szst4(i1,i2,i3)
sztt   = sztt4(i1,i2,i3)
txr    = txr4(i1,i2,i3)
txrr   = txrr4(i1,i2,i3)
txs    = txs4(i1,i2,i3)
txrs   = txrs4(i1,i2,i3)
txss   = txss4(i1,i2,i3)
txt    = txt4(i1,i2,i3)
txrt   = txrt4(i1,i2,i3)
txst   = txst4(i1,i2,i3)
txtt   = txtt4(i1,i2,i3)
tyr    = tyr4(i1,i2,i3)
tyrr   = tyrr4(i1,i2,i3)
tys    = tys4(i1,i2,i3)
tyrs   = tyrs4(i1,i2,i3)
tyss   = tyss4(i1,i2,i3)
tyt    = tyt4(i1,i2,i3)
tyrt   = tyrt4(i1,i2,i3)
tyst   = tyst4(i1,i2,i3)
tytt   = tytt4(i1,i2,i3)
tzr    = tzr4(i1,i2,i3)
tzrr   = tzrr4(i1,i2,i3)
tzs    = tzs4(i1,i2,i3)
tzrs   = tzrs4(i1,i2,i3)
tzss   = tzss4(i1,i2,i3)
tzt    = tzt4(i1,i2,i3)
tzrt   = tzrt4(i1,i2,i3)
tzst   = tzst4(i1,i2,i3)
tztt   = tztt4(i1,i2,i3)

! ---------- Spatial derivatives of metrics rx, sx, ry, ... ---------
rxi = rx(i1,i2,i3)
ryi = ry(i1,i2,i3)
sxi = sx(i1,i2,i3)
syi = sy(i1,i2,i3)
rzi = rz(i1,i2,i3)
szi = sz(i1,i2,i3)
txi = tx(i1,i2,i3)
tyi = ty(i1,i2,i3)
tzi = tz(i1,i2,i3)
rxx = rxi*rxr+sxi*rxs+txi*rxt
rxxx = rxi**2*rxrr+2.*rxi*sxi*rxrs+2.*rxi*txi*rxrt+sxi**2*rxss+2.*sxi*txi*rxst+txi**2*rxtt+rxx*rxr+sxx*rxs+txx*rxt
rxy = ryi*rxr+syi*rxs+tyi*rxt
rxxy = ryi*rxi*rxrr+(rxi*syi+ryi*sxi)*rxrs+syi*sxi*rxss+(rxi*tyi+ryi*txi)*rxrt+(sxi*tyi+syi*txi)*rxst+tyi*txi*rxtt+rxy*rxr+sxy*rxs+txy*rxt
rxyy = ryi**2*rxrr+2.*ryi*syi*rxrs+2.*ryi*tyi*rxrt+syi**2*rxss+2.*syi*tyi*rxst+tyi**2*rxtt+ryy*rxr+syy*rxs+tyy*rxt
rxz = rzi*rxr+szi*rxs+tzi*rxt
rxxz = rzi*rxi*rxrr+(rxi*szi+rzi*sxi)*rxrs+szi*sxi*rxss+(rxi*tzi+rzi*txi)*rxrt+(sxi*tzi+szi*txi)*rxst+tzi*txi*rxtt+rxz*rxr+sxz*rxs+txz*rxt
rxyz = rzi*ryi*rxrr+(ryi*szi+rzi*syi)*rxrs+szi*syi*rxss+(ryi*tzi+rzi*tyi)*rxrt+(syi*tzi+szi*tyi)*rxst+tzi*tyi*rxtt+ryz*rxr+syz*rxs+tyz*rxt
rxzz = rzi**2*rxrr+2.*rzi*szi*rxrs+2.*rzi*tzi*rxrt+szi**2*rxss+2.*szi*tzi*rxst+tzi**2*rxtt+rzz*rxr+szz*rxs+tzz*rxt
ryy = ryi*ryr+syi*rys+tyi*ryt
ryyy = ryi**2*ryrr+2.*ryi*syi*ryrs+2.*ryi*tyi*ryrt+syi**2*ryss+2.*syi*tyi*ryst+tyi**2*rytt+ryy*ryr+syy*rys+tyy*ryt
ryz = rzi*ryr+szi*rys+tzi*ryt
ryyz = rzi*ryi*ryrr+(ryi*szi+rzi*syi)*ryrs+szi*syi*ryss+(ryi*tzi+rzi*tyi)*ryrt+(syi*tzi+szi*tyi)*ryst+tzi*tyi*rytt+ryz*ryr+syz*rys+tyz*ryt
ryzz = rzi**2*ryrr+2.*rzi*szi*ryrs+2.*rzi*tzi*ryrt+szi**2*ryss+2.*szi*tzi*ryst+tzi**2*rytt+rzz*ryr+szz*rys+tzz*ryt
rzz = rzi*rzr+szi*rzs+tzi*rzt
rzzz = rzi**2*rzrr+2.*rzi*szi*rzrs+2.*rzi*tzi*rzrt+szi**2*rzss+2.*szi*tzi*rzst+tzi**2*rztt+rzz*rzr+szz*rzs+tzz*rzt
sxx = rxi*sxr+sxi*sxs+txi*sxt
sxxx = rxi**2*sxrr+2.*rxi*sxi*sxrs+2.*rxi*txi*sxrt+sxi**2*sxss+2.*sxi*txi*sxst+txi**2*sxtt+rxx*sxr+sxx*sxs+txx*sxt
sxy = ryi*sxr+syi*sxs+tyi*sxt
sxxy = ryi*rxi*sxrr+(rxi*syi+ryi*sxi)*sxrs+syi*sxi*sxss+(rxi*tyi+ryi*txi)*sxrt+(sxi*tyi+syi*txi)*sxst+tyi*txi*sxtt+rxy*sxr+sxy*sxs+txy*sxt
sxyy = ryi**2*sxrr+2.*ryi*syi*sxrs+2.*ryi*tyi*sxrt+syi**2*sxss+2.*syi*tyi*sxst+tyi**2*sxtt+ryy*sxr+syy*sxs+tyy*sxt
sxz = rzi*sxr+szi*sxs+tzi*sxt
sxxz = rzi*rxi*sxrr+(rxi*szi+rzi*sxi)*sxrs+szi*sxi*sxss+(rxi*tzi+rzi*txi)*sxrt+(sxi*tzi+szi*txi)*sxst+tzi*txi*sxtt+rxz*sxr+sxz*sxs+txz*sxt
sxyz = rzi*ryi*sxrr+(ryi*szi+rzi*syi)*sxrs+szi*syi*sxss+(ryi*tzi+rzi*tyi)*sxrt+(syi*tzi+szi*tyi)*sxst+tzi*tyi*sxtt+ryz*sxr+syz*sxs+tyz*sxt
sxzz = rzi**2*sxrr+2.*rzi*szi*sxrs+2.*rzi*tzi*sxrt+szi**2*sxss+2.*szi*tzi*sxst+tzi**2*sxtt+rzz*sxr+szz*sxs+tzz*sxt
syy = ryi*syr+syi*sys+tyi*syt
syyy = ryi**2*syrr+2.*ryi*syi*syrs+2.*ryi*tyi*syrt+syi**2*syss+2.*syi*tyi*syst+tyi**2*sytt+ryy*syr+syy*sys+tyy*syt
syz = rzi*syr+szi*sys+tzi*syt
syyz = rzi*ryi*syrr+(ryi*szi+rzi*syi)*syrs+szi*syi*syss+(ryi*tzi+rzi*tyi)*syrt+(syi*tzi+szi*tyi)*syst+tzi*tyi*sytt+ryz*syr+syz*sys+tyz*syt
syzz = rzi**2*syrr+2.*rzi*szi*syrs+2.*rzi*tzi*syrt+szi**2*syss+2.*szi*tzi*syst+tzi**2*sytt+rzz*syr+szz*sys+tzz*syt
szz = rzi*szr+szi*szs+tzi*szt
szzz = rzi**2*szrr+2.*rzi*szi*szrs+2.*rzi*tzi*szrt+szi**2*szss+2.*szi*tzi*szst+tzi**2*sztt+rzz*szr+szz*szs+tzz*szt
txx = rxi*txr+sxi*txs+txi*txt
txxx = rxi**2*txrr+2.*rxi*sxi*txrs+2.*rxi*txi*txrt+sxi**2*txss+2.*sxi*txi*txst+txi**2*txtt+rxx*txr+sxx*txs+txx*txt
txy = ryi*txr+syi*txs+tyi*txt
txxy = ryi*rxi*txrr+(rxi*syi+ryi*sxi)*txrs+syi*sxi*txss+(rxi*tyi+ryi*txi)*txrt+(sxi*tyi+syi*txi)*txst+tyi*txi*txtt+rxy*txr+sxy*txs+txy*txt
txyy = ryi**2*txrr+2.*ryi*syi*txrs+2.*ryi*tyi*txrt+syi**2*txss+2.*syi*tyi*txst+tyi**2*txtt+ryy*txr+syy*txs+tyy*txt
txz = rzi*txr+szi*txs+tzi*txt
txxz = rzi*rxi*txrr+(rxi*szi+rzi*sxi)*txrs+szi*sxi*txss+(rxi*tzi+rzi*txi)*txrt+(sxi*tzi+szi*txi)*txst+tzi*txi*txtt+rxz*txr+sxz*txs+txz*txt
txyz = rzi*ryi*txrr+(ryi*szi+rzi*syi)*txrs+szi*syi*txss+(ryi*tzi+rzi*tyi)*txrt+(syi*tzi+szi*tyi)*txst+tzi*tyi*txtt+ryz*txr+syz*txs+tyz*txt
txzz = rzi**2*txrr+2.*rzi*szi*txrs+2.*rzi*tzi*txrt+szi**2*txss+2.*szi*tzi*txst+tzi**2*txtt+rzz*txr+szz*txs+tzz*txt
tyy = ryi*tyr+syi*tys+tyi*tyt
tyyy = ryi**2*tyrr+2.*ryi*syi*tyrs+2.*ryi*tyi*tyrt+syi**2*tyss+2.*syi*tyi*tyst+tyi**2*tytt+ryy*tyr+syy*tys+tyy*tyt
tyz = rzi*tyr+szi*tys+tzi*tyt
tyyz = rzi*ryi*tyrr+(ryi*szi+rzi*syi)*tyrs+szi*syi*tyss+(ryi*tzi+rzi*tyi)*tyrt+(syi*tzi+szi*tyi)*tyst+tzi*tyi*tytt+ryz*tyr+syz*tys+tyz*tyt
tyzz = rzi**2*tyrr+2.*rzi*szi*tyrs+2.*rzi*tzi*tyrt+szi**2*tyss+2.*szi*tzi*tyst+tzi**2*tytt+rzz*tyr+szz*tys+tzz*tyt
tzz = rzi*tzr+szi*tzs+tzi*tzt
tzzz = rzi**2*tzrr+2.*rzi*szi*tzrs+2.*rzi*tzi*tzrt+szi**2*tzss+2.*szi*tzi*tzst+tzi**2*tztt+rzz*tzr+szz*tzs+tzz*tzt
#End
! ---- end OPTION eq evalMetrics ---

! ---------- Third spatial derivatives of u ---------
uxxx = rxi**3*urrr+3.*rxi**2*sxi*urrs+3.*rxi*sxi**2*urss+sxi**3*usss+3.*rxi**2*txi*urrt+6.*rxi*sxi*txi*urst+3.*txi*sxi**2*usst+3.*rxi*txi**2*urtt+3.*txi**2*sxi*ustt+txi**3*uttt+3.*rxi*rxx*urr+(3.*rxi*sxx+3.*rxx*sxi)*urs+3.*sxx*sxi*uss+(3.*rxi*txx+3.*rxx*txi)*urt+(3.*sxi*txx+3.*sxx*txi)*ust+3.*txx*txi*utt+rxxx*ur+sxxx*us+txxx*ut
uxxy = ryi*rxi**2*urrr+(rxi**2*syi+2.*rxi*ryi*sxi)*urrs+(2.*rxi*sxi*syi+ryi*sxi**2)*urss+syi*sxi**2*usss+(rxi**2*tyi+2.*rxi*ryi*txi)*urrt+(2.*rxi*sxi*tyi+2.*rxi*syi*txi+2.*ryi*sxi*txi)*urst+(sxi**2*tyi+2.*sxi*syi*txi)*usst+(2.*rxi*txi*tyi+ryi*txi**2)*urtt+(2.*sxi*txi*tyi+syi*txi**2)*ustt+tyi*txi**2*uttt+(2.*rxi*rxy+rxx*ryi)*urr+(2.*rxi*sxy+rxx*syi+2.*rxy*sxi+ryi*sxx)*urs+(2.*sxi*sxy+sxx*syi)*uss+(2.*rxi*txy+rxx*tyi+2.*rxy*txi+ryi*txx)*urt+(2.*sxi*txy+sxx*tyi+2.*sxy*txi+syi*txx)*ust+(2.*txi*txy+txx*tyi)*utt+rxxy*ur+sxxy*us+txxy*ut
uxyy = ryi**2*rxi*urrr+(syi*ryi*rxi+ryi*(rxi*syi+ryi*sxi))*urrs+(ryi*syi*sxi+syi*(rxi*syi+ryi*sxi))*urss+syi**2*sxi*usss+(tyi*ryi*rxi+ryi*(rxi*tyi+ryi*txi))*urrt+(tyi*(rxi*syi+ryi*sxi)+syi*(rxi*tyi+ryi*txi)+ryi*(sxi*tyi+syi*txi))*urst+(syi*(sxi*tyi+syi*txi)+tyi*syi*sxi)*usst+(tyi*(rxi*tyi+ryi*txi)+ryi*tyi*txi)*urtt+(tyi*(sxi*tyi+syi*txi)+syi*tyi*txi)*ustt+tyi**2*txi*uttt+(rxi*ryy+2.*rxy*ryi)*urr+(rxi*syy+2.*rxy*syi+2.*ryi*sxy+ryy*sxi)*urs+(sxi*syy+2.*sxy*syi)*uss+(rxi*tyy+2.*rxy*tyi+2.*ryi*txy+ryy*txi)*urt+(sxi*tyy+2.*sxy*tyi+2.*syi*txy+syy*txi)*ust+(txi*tyy+2.*txy*tyi)*utt+rxyy*ur+sxyy*us+txyy*ut
uyyy = ryi**3*urrr+3.*ryi**2*syi*urrs+3.*ryi*syi**2*urss+syi**3*usss+3.*ryi**2*tyi*urrt+6.*ryi*tyi*syi*urst+3.*tyi*syi**2*usst+3.*ryi*tyi**2*urtt+3.*tyi**2*syi*ustt+tyi**3*uttt+3.*ryi*ryy*urr+(3.*ryi*syy+3.*ryy*syi)*urs+3.*syy*syi*uss+(3.*ryi*tyy+3.*ryy*tyi)*urt+(3.*syi*tyy+3.*syy*tyi)*ust+3.*tyy*tyi*utt+ryyy*ur+syyy*us+tyyy*ut
uxxz = rzi*rxi**2*urrr+(rxi**2*szi+2.*rxi*rzi*sxi)*urrs+(2.*rxi*sxi*szi+rzi*sxi**2)*urss+szi*sxi**2*usss+(rxi**2*tzi+2.*rxi*rzi*txi)*urrt+(2.*rxi*sxi*tzi+2.*rxi*szi*txi+2.*rzi*sxi*txi)*urst+(sxi**2*tzi+2.*sxi*szi*txi)*usst+(2.*rxi*txi*tzi+rzi*txi**2)*urtt+(2.*sxi*txi*tzi+szi*txi**2)*ustt+tzi*txi**2*uttt+(2.*rxi*rxz+rxx*rzi)*urr+(2.*rxi*sxz+rxx*szi+2.*rxz*sxi+rzi*sxx)*urs+(2.*sxi*sxz+sxx*szi)*uss+(2.*rxi*txz+rxx*tzi+2.*rxz*txi+rzi*txx)*urt+(2.*sxi*txz+sxx*tzi+2.*sxz*txi+szi*txx)*ust+(2.*txi*txz+txx*tzi)*utt+rxxz*ur+sxxz*us+txxz*ut
uxyz = rzi*ryi*rxi*urrr+(szi*ryi*rxi+rzi*(rxi*syi+ryi*sxi))*urrs+(rzi*syi*sxi+szi*(rxi*syi+ryi*sxi))*urss+szi*syi*sxi*usss+(tzi*ryi*rxi+rzi*(rxi*tyi+ryi*txi))*urrt+(tzi*(rxi*syi+ryi*sxi)+szi*(rxi*tyi+ryi*txi)+rzi*(sxi*tyi+syi*txi))*urst+(szi*(sxi*tyi+syi*txi)+tzi*syi*sxi)*usst+(tzi*(rxi*tyi+ryi*txi)+rzi*tyi*txi)*urtt+(tzi*(sxi*tyi+syi*txi)+szi*tyi*txi)*ustt+tzi*tyi*txi*uttt+(rxi*ryz+rxy*rzi+rxz*ryi)*urr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*urs+(sxi*syz+sxy*szi+sxz*syi)*uss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*urt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*ust+(txi*tyz+txy*tzi+txz*tyi)*utt+rxyz*ur+sxyz*us+txyz*ut
uyyz = rzi*ryi**2*urrr+(ryi**2*szi+2.*ryi*rzi*syi)*urrs+(2.*ryi*syi*szi+rzi*syi**2)*urss+szi*syi**2*usss+(ryi**2*tzi+2.*ryi*rzi*tyi)*urrt+(2.*ryi*syi*tzi+2.*ryi*szi*tyi+2.*rzi*syi*tyi)*urst+(syi**2*tzi+2.*syi*szi*tyi)*usst+(2.*ryi*tyi*tzi+rzi*tyi**2)*urtt+(2.*syi*tyi*tzi+szi*tyi**2)*ustt+tzi*tyi**2*uttt+(2.*ryi*ryz+ryy*rzi)*urr+(2.*ryi*syz+ryy*szi+2.*ryz*syi+rzi*syy)*urs+(2.*syi*syz+syy*szi)*uss+(2.*ryi*tyz+ryy*tzi+2.*ryz*tyi+rzi*tyy)*urt+(2.*syi*tyz+syy*tzi+2.*syz*tyi+szi*tyy)*ust+(2.*tyi*tyz+tyy*tzi)*utt+ryyz*ur+syyz*us+tyyz*ut
uxzz = rzi**2*rxi*urrr+(szi*rzi*rxi+rzi*(rxi*szi+rzi*sxi))*urrs+(rzi*szi*sxi+szi*(rxi*szi+rzi*sxi))*urss+szi**2*sxi*usss+(tzi*rzi*rxi+rzi*(rxi*tzi+rzi*txi))*urrt+(tzi*(rxi*szi+rzi*sxi)+szi*(rxi*tzi+rzi*txi)+rzi*(sxi*tzi+szi*txi))*urst+(szi*(sxi*tzi+szi*txi)+tzi*szi*sxi)*usst+(tzi*(rxi*tzi+rzi*txi)+rzi*tzi*txi)*urtt+(tzi*(sxi*tzi+szi*txi)+szi*tzi*txi)*ustt+tzi**2*txi*uttt+(rxi*rzz+2.*rxz*rzi)*urr+(rxi*szz+2.*rxz*szi+2.*rzi*sxz+rzz*sxi)*urs+(sxi*szz+2.*sxz*szi)*uss+(rxi*tzz+2.*rxz*tzi+2.*rzi*txz+rzz*txi)*urt+(sxi*tzz+2.*sxz*tzi+2.*szi*txz+szz*txi)*ust+(txi*tzz+2.*txz*tzi)*utt+rxzz*ur+sxzz*us+txzz*ut
uyzz = rzi**2*ryi*urrr+(szi*rzi*ryi+rzi*(ryi*szi+rzi*syi))*urrs+(rzi*szi*syi+szi*(ryi*szi+rzi*syi))*urss+szi**2*syi*usss+(tzi*rzi*ryi+rzi*(ryi*tzi+rzi*tyi))*urrt+(tzi*(ryi*szi+rzi*syi)+szi*(ryi*tzi+rzi*tyi)+rzi*(syi*tzi+szi*tyi))*urst+(szi*(syi*tzi+szi*tyi)+tzi*szi*syi)*usst+(tzi*(ryi*tzi+rzi*tyi)+rzi*tzi*tyi)*urtt+(tzi*(syi*tzi+szi*tyi)+szi*tzi*tyi)*ustt+tzi**2*tyi*uttt+(ryi*rzz+2.*ryz*rzi)*urr+(ryi*szz+2.*ryz*szi+2.*rzi*syz+rzz*syi)*urs+(syi*szz+2.*syz*szi)*uss+(ryi*tzz+2.*ryz*tzi+2.*rzi*tyz+rzz*tyi)*urt+(syi*tzz+2.*syz*tzi+2.*szi*tyz+szz*tyi)*ust+(tyi*tzz+2.*tyz*tzi)*utt+ryzz*ur+syzz*us+tyzz*ut
uzzz = rzi**3*urrr+3.*rzi**2*szi*urrs+3.*rzi*szi**2*urss+szi**3*usss+3.*rzi**2*tzi*urrt+6.*rzi*tzi*szi*urst+3.*tzi*szi**2*usst+3.*rzi*tzi**2*urtt+3.*tzi**2*szi*ustt+tzi**3*uttt+3.*rzi*rzz*urr+(3.*rzi*szz+3.*rzz*szi)*urs+3.*szz*szi*uss+(3.*rzi*tzz+3.*rzz*tzi)*urt+(3.*szi*tzz+3.*szz*tzi)*ust+3.*tzz*tzi*utt+rzzz*ur+szzz*us+tzzz*ut
! ---------- END CURVILINEAR  ---------
#End
#endMacro
