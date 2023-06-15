! ****** File written by makeGetDerivativesMacros.maple  ******

#beginMacro defineThirdParameticDerivativesComponents0(u)
 #defineMacro u ## r2Macro(i1,i2,i3) (-u(i1-1,i2,i3)+u(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rr2Macro(i1,i2,i3) (u(i1-1,i2,i3)-2.*u(i1,i2,i3)+u(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrr2Macro(i1,i2,i3) (-u(i1-2,i2,i3)+2.*u(i1-1,i2,i3)-2.*u(i1+1,i2,i3)+u(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrr2Macro(i1,i2,i3) (u(i1-2,i2,i3)-4.*u(i1-1,i2,i3)+6.*u(i1,i2,i3)-4.*u(i1+1,i2,i3)+u(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## s2Macro(i1,i2,i3) (-u(i1,i2-1,i3)+u(i1,i2+1,i3))/(2.*dr(1))
 #defineMacro u ## rs2Macro(i1,i2,i3) (-u ## s2Macro(i1-1,i2,i3)+u ## s2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrs2Macro(i1,i2,i3) (u ## s2Macro(i1-1,i2,i3)-2.*u ## s2Macro(i1,i2,i3)+u ## s2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrs2Macro(i1,i2,i3) (-u ## s2Macro(i1-2,i2,i3)+2.*u ## s2Macro(i1-1,i2,i3)-2.*u ## s2Macro(i1+1,i2,i3)+u ## s2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrs2Macro(i1,i2,i3) (u ## s2Macro(i1-2,i2,i3)-4.*u ## s2Macro(i1-1,i2,i3)+6.*u ## s2Macro(i1,i2,i3)-4.*u ## s2Macro(i1+1,i2,i3)+u ## s2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## ss2Macro(i1,i2,i3) (u(i1,i2-1,i3)-2.*u(i1,i2,i3)+u(i1,i2+1,i3))/(dr(1)**2)
 #defineMacro u ## rss2Macro(i1,i2,i3) (-u ## ss2Macro(i1-1,i2,i3)+u ## ss2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrss2Macro(i1,i2,i3) (u ## ss2Macro(i1-1,i2,i3)-2.*u ## ss2Macro(i1,i2,i3)+u ## ss2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrss2Macro(i1,i2,i3) (-u ## ss2Macro(i1-2,i2,i3)+2.*u ## ss2Macro(i1-1,i2,i3)-2.*u ## ss2Macro(i1+1,i2,i3)+u ## ss2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrss2Macro(i1,i2,i3) (u ## ss2Macro(i1-2,i2,i3)-4.*u ## ss2Macro(i1-1,i2,i3)+6.*u ## ss2Macro(i1,i2,i3)-4.*u ## ss2Macro(i1+1,i2,i3)+u ## ss2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## sss2Macro(i1,i2,i3) (-u(i1,i2-2,i3)+2.*u(i1,i2-1,i3)-2.*u(i1,i2+1,i3)+u(i1,i2+2,i3))/(2.*dr(1)**3)
 #defineMacro u ## rsss2Macro(i1,i2,i3) (-u ## sss2Macro(i1-1,i2,i3)+u ## sss2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrsss2Macro(i1,i2,i3) (u ## sss2Macro(i1-1,i2,i3)-2.*u ## sss2Macro(i1,i2,i3)+u ## sss2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrsss2Macro(i1,i2,i3) (-u ## sss2Macro(i1-2,i2,i3)+2.*u ## sss2Macro(i1-1,i2,i3)-2.*u ## sss2Macro(i1+1,i2,i3)+u ## sss2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsss2Macro(i1,i2,i3) (u ## sss2Macro(i1-2,i2,i3)-4.*u ## sss2Macro(i1-1,i2,i3)+6.*u ## sss2Macro(i1,i2,i3)-4.*u ## sss2Macro(i1+1,i2,i3)+u ## sss2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## ssss2Macro(i1,i2,i3) (u(i1,i2-2,i3)-4.*u(i1,i2-1,i3)+6.*u(i1,i2,i3)-4.*u(i1,i2+1,i3)+u(i1,i2+2,i3))/(dr(1)**4)
 #defineMacro u ## rssss2Macro(i1,i2,i3) (-u ## ssss2Macro(i1-1,i2,i3)+u ## ssss2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrssss2Macro(i1,i2,i3) (u ## ssss2Macro(i1-1,i2,i3)-2.*u ## ssss2Macro(i1,i2,i3)+u ## ssss2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrssss2Macro(i1,i2,i3) (-u ## ssss2Macro(i1-2,i2,i3)+2.*u ## ssss2Macro(i1-1,i2,i3)-2.*u ## ssss2Macro(i1+1,i2,i3)+u ## ssss2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssss2Macro(i1,i2,i3) (u ## ssss2Macro(i1-2,i2,i3)-4.*u ## ssss2Macro(i1-1,i2,i3)+6.*u ## ssss2Macro(i1,i2,i3)-4.*u ## ssss2Macro(i1+1,i2,i3)+u ## ssss2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## t2Macro(i1,i2,i3) (-u(i1,i2,i3-1)+u(i1,i2,i3+1))/(2.*dr(2))
 #defineMacro u ## rt2Macro(i1,i2,i3) (-u ## t2Macro(i1-1,i2,i3)+u ## t2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrt2Macro(i1,i2,i3) (u ## t2Macro(i1-1,i2,i3)-2.*u ## t2Macro(i1,i2,i3)+u ## t2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrt2Macro(i1,i2,i3) (-u ## t2Macro(i1-2,i2,i3)+2.*u ## t2Macro(i1-1,i2,i3)-2.*u ## t2Macro(i1+1,i2,i3)+u ## t2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrt2Macro(i1,i2,i3) (u ## t2Macro(i1-2,i2,i3)-4.*u ## t2Macro(i1-1,i2,i3)+6.*u ## t2Macro(i1,i2,i3)-4.*u ## t2Macro(i1+1,i2,i3)+u ## t2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## st2Macro(i1,i2,i3) (-u ## t2Macro(i1,i2-1,i3)+u ## t2Macro(i1,i2+1,i3))/(2.*dr(1))
 #defineMacro u ## rst2Macro(i1,i2,i3) (-u ## st2Macro(i1-1,i2,i3)+u ## st2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrst2Macro(i1,i2,i3) (u ## st2Macro(i1-1,i2,i3)-2.*u ## st2Macro(i1,i2,i3)+u ## st2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrst2Macro(i1,i2,i3) (-u ## st2Macro(i1-2,i2,i3)+2.*u ## st2Macro(i1-1,i2,i3)-2.*u ## st2Macro(i1+1,i2,i3)+u ## st2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrst2Macro(i1,i2,i3) (u ## st2Macro(i1-2,i2,i3)-4.*u ## st2Macro(i1-1,i2,i3)+6.*u ## st2Macro(i1,i2,i3)-4.*u ## st2Macro(i1+1,i2,i3)+u ## st2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## sst2Macro(i1,i2,i3) (u ## t2Macro(i1,i2-1,i3)-2.*u ## t2Macro(i1,i2,i3)+u ## t2Macro(i1,i2+1,i3))/(dr(1)**2)
 #defineMacro u ## rsst2Macro(i1,i2,i3) (-u ## sst2Macro(i1-1,i2,i3)+u ## sst2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrsst2Macro(i1,i2,i3) (u ## sst2Macro(i1-1,i2,i3)-2.*u ## sst2Macro(i1,i2,i3)+u ## sst2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrsst2Macro(i1,i2,i3) (-u ## sst2Macro(i1-2,i2,i3)+2.*u ## sst2Macro(i1-1,i2,i3)-2.*u ## sst2Macro(i1+1,i2,i3)+u ## sst2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsst2Macro(i1,i2,i3) (u ## sst2Macro(i1-2,i2,i3)-4.*u ## sst2Macro(i1-1,i2,i3)+6.*u ## sst2Macro(i1,i2,i3)-4.*u ## sst2Macro(i1+1,i2,i3)+u ## sst2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## ssst2Macro(i1,i2,i3) (-u ## t2Macro(i1,i2-2,i3)+2.*u ## t2Macro(i1,i2-1,i3)-2.*u ## t2Macro(i1,i2+1,i3)+u ## t2Macro(i1,i2+2,i3))/(2.*dr(1)**3)
 #defineMacro u ## rssst2Macro(i1,i2,i3) (-u ## ssst2Macro(i1-1,i2,i3)+u ## ssst2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrssst2Macro(i1,i2,i3) (u ## ssst2Macro(i1-1,i2,i3)-2.*u ## ssst2Macro(i1,i2,i3)+u ## ssst2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrssst2Macro(i1,i2,i3) (-u ## ssst2Macro(i1-2,i2,i3)+2.*u ## ssst2Macro(i1-1,i2,i3)-2.*u ## ssst2Macro(i1+1,i2,i3)+u ## ssst2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssst2Macro(i1,i2,i3) (u ## ssst2Macro(i1-2,i2,i3)-4.*u ## ssst2Macro(i1-1,i2,i3)+6.*u ## ssst2Macro(i1,i2,i3)-4.*u ## ssst2Macro(i1+1,i2,i3)+u ## ssst2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## sssst2Macro(i1,i2,i3) (u ## t2Macro(i1,i2-2,i3)-4.*u ## t2Macro(i1,i2-1,i3)+6.*u ## t2Macro(i1,i2,i3)-4.*u ## t2Macro(i1,i2+1,i3)+u ## t2Macro(i1,i2+2,i3))/(dr(1)**4)
 #defineMacro u ## rsssst2Macro(i1,i2,i3) (-u ## sssst2Macro(i1-1,i2,i3)+u ## sssst2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrsssst2Macro(i1,i2,i3) (u ## sssst2Macro(i1-1,i2,i3)-2.*u ## sssst2Macro(i1,i2,i3)+u ## sssst2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrsssst2Macro(i1,i2,i3) (-u ## sssst2Macro(i1-2,i2,i3)+2.*u ## sssst2Macro(i1-1,i2,i3)-2.*u ## sssst2Macro(i1+1,i2,i3)+u ## sssst2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsssst2Macro(i1,i2,i3) (u ## sssst2Macro(i1-2,i2,i3)-4.*u ## sssst2Macro(i1-1,i2,i3)+6.*u ## sssst2Macro(i1,i2,i3)-4.*u ## sssst2Macro(i1+1,i2,i3)+u ## sssst2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## tt2Macro(i1,i2,i3) (u(i1,i2,i3-1)-2.*u(i1,i2,i3)+u(i1,i2,i3+1))/(dr(2)**2)
 #defineMacro u ## rtt2Macro(i1,i2,i3) (-u ## tt2Macro(i1-1,i2,i3)+u ## tt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrtt2Macro(i1,i2,i3) (u ## tt2Macro(i1-1,i2,i3)-2.*u ## tt2Macro(i1,i2,i3)+u ## tt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrtt2Macro(i1,i2,i3) (-u ## tt2Macro(i1-2,i2,i3)+2.*u ## tt2Macro(i1-1,i2,i3)-2.*u ## tt2Macro(i1+1,i2,i3)+u ## tt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrtt2Macro(i1,i2,i3) (u ## tt2Macro(i1-2,i2,i3)-4.*u ## tt2Macro(i1-1,i2,i3)+6.*u ## tt2Macro(i1,i2,i3)-4.*u ## tt2Macro(i1+1,i2,i3)+u ## tt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## stt2Macro(i1,i2,i3) (-u ## tt2Macro(i1,i2-1,i3)+u ## tt2Macro(i1,i2+1,i3))/(2.*dr(1))
 #defineMacro u ## rstt2Macro(i1,i2,i3) (-u ## stt2Macro(i1-1,i2,i3)+u ## stt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrstt2Macro(i1,i2,i3) (u ## stt2Macro(i1-1,i2,i3)-2.*u ## stt2Macro(i1,i2,i3)+u ## stt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrstt2Macro(i1,i2,i3) (-u ## stt2Macro(i1-2,i2,i3)+2.*u ## stt2Macro(i1-1,i2,i3)-2.*u ## stt2Macro(i1+1,i2,i3)+u ## stt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrstt2Macro(i1,i2,i3) (u ## stt2Macro(i1-2,i2,i3)-4.*u ## stt2Macro(i1-1,i2,i3)+6.*u ## stt2Macro(i1,i2,i3)-4.*u ## stt2Macro(i1+1,i2,i3)+u ## stt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## sstt2Macro(i1,i2,i3) (u ## tt2Macro(i1,i2-1,i3)-2.*u ## tt2Macro(i1,i2,i3)+u ## tt2Macro(i1,i2+1,i3))/(dr(1)**2)
 #defineMacro u ## rsstt2Macro(i1,i2,i3) (-u ## sstt2Macro(i1-1,i2,i3)+u ## sstt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrsstt2Macro(i1,i2,i3) (u ## sstt2Macro(i1-1,i2,i3)-2.*u ## sstt2Macro(i1,i2,i3)+u ## sstt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrsstt2Macro(i1,i2,i3) (-u ## sstt2Macro(i1-2,i2,i3)+2.*u ## sstt2Macro(i1-1,i2,i3)-2.*u ## sstt2Macro(i1+1,i2,i3)+u ## sstt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsstt2Macro(i1,i2,i3) (u ## sstt2Macro(i1-2,i2,i3)-4.*u ## sstt2Macro(i1-1,i2,i3)+6.*u ## sstt2Macro(i1,i2,i3)-4.*u ## sstt2Macro(i1+1,i2,i3)+u ## sstt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## ssstt2Macro(i1,i2,i3) (-u ## tt2Macro(i1,i2-2,i3)+2.*u ## tt2Macro(i1,i2-1,i3)-2.*u ## tt2Macro(i1,i2+1,i3)+u ## tt2Macro(i1,i2+2,i3))/(2.*dr(1)**3)
 #defineMacro u ## rssstt2Macro(i1,i2,i3) (-u ## ssstt2Macro(i1-1,i2,i3)+u ## ssstt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrssstt2Macro(i1,i2,i3) (u ## ssstt2Macro(i1-1,i2,i3)-2.*u ## ssstt2Macro(i1,i2,i3)+u ## ssstt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrssstt2Macro(i1,i2,i3) (-u ## ssstt2Macro(i1-2,i2,i3)+2.*u ## ssstt2Macro(i1-1,i2,i3)-2.*u ## ssstt2Macro(i1+1,i2,i3)+u ## ssstt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssstt2Macro(i1,i2,i3) (u ## ssstt2Macro(i1-2,i2,i3)-4.*u ## ssstt2Macro(i1-1,i2,i3)+6.*u ## ssstt2Macro(i1,i2,i3)-4.*u ## ssstt2Macro(i1+1,i2,i3)+u ## ssstt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## sssstt2Macro(i1,i2,i3) (u ## tt2Macro(i1,i2-2,i3)-4.*u ## tt2Macro(i1,i2-1,i3)+6.*u ## tt2Macro(i1,i2,i3)-4.*u ## tt2Macro(i1,i2+1,i3)+u ## tt2Macro(i1,i2+2,i3))/(dr(1)**4)
 #defineMacro u ## rsssstt2Macro(i1,i2,i3) (-u ## sssstt2Macro(i1-1,i2,i3)+u ## sssstt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrsssstt2Macro(i1,i2,i3) (u ## sssstt2Macro(i1-1,i2,i3)-2.*u ## sssstt2Macro(i1,i2,i3)+u ## sssstt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrsssstt2Macro(i1,i2,i3) (-u ## sssstt2Macro(i1-2,i2,i3)+2.*u ## sssstt2Macro(i1-1,i2,i3)-2.*u ## sssstt2Macro(i1+1,i2,i3)+u ## sssstt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsssstt2Macro(i1,i2,i3) (u ## sssstt2Macro(i1-2,i2,i3)-4.*u ## sssstt2Macro(i1-1,i2,i3)+6.*u ## sssstt2Macro(i1,i2,i3)-4.*u ## sssstt2Macro(i1+1,i2,i3)+u ## sssstt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## ttt2Macro(i1,i2,i3) (-u(i1,i2,i3-2)+2.*u(i1,i2,i3-1)-2.*u(i1,i2,i3+1)+u(i1,i2,i3+2))/(2.*dr(2)**3)
 #defineMacro u ## rttt2Macro(i1,i2,i3) (-u ## ttt2Macro(i1-1,i2,i3)+u ## ttt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrttt2Macro(i1,i2,i3) (u ## ttt2Macro(i1-1,i2,i3)-2.*u ## ttt2Macro(i1,i2,i3)+u ## ttt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrttt2Macro(i1,i2,i3) (-u ## ttt2Macro(i1-2,i2,i3)+2.*u ## ttt2Macro(i1-1,i2,i3)-2.*u ## ttt2Macro(i1+1,i2,i3)+u ## ttt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrttt2Macro(i1,i2,i3) (u ## ttt2Macro(i1-2,i2,i3)-4.*u ## ttt2Macro(i1-1,i2,i3)+6.*u ## ttt2Macro(i1,i2,i3)-4.*u ## ttt2Macro(i1+1,i2,i3)+u ## ttt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## sttt2Macro(i1,i2,i3) (-u ## ttt2Macro(i1,i2-1,i3)+u ## ttt2Macro(i1,i2+1,i3))/(2.*dr(1))
 #defineMacro u ## rsttt2Macro(i1,i2,i3) (-u ## sttt2Macro(i1-1,i2,i3)+u ## sttt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrsttt2Macro(i1,i2,i3) (u ## sttt2Macro(i1-1,i2,i3)-2.*u ## sttt2Macro(i1,i2,i3)+u ## sttt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrsttt2Macro(i1,i2,i3) (-u ## sttt2Macro(i1-2,i2,i3)+2.*u ## sttt2Macro(i1-1,i2,i3)-2.*u ## sttt2Macro(i1+1,i2,i3)+u ## sttt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsttt2Macro(i1,i2,i3) (u ## sttt2Macro(i1-2,i2,i3)-4.*u ## sttt2Macro(i1-1,i2,i3)+6.*u ## sttt2Macro(i1,i2,i3)-4.*u ## sttt2Macro(i1+1,i2,i3)+u ## sttt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## ssttt2Macro(i1,i2,i3) (u ## ttt2Macro(i1,i2-1,i3)-2.*u ## ttt2Macro(i1,i2,i3)+u ## ttt2Macro(i1,i2+1,i3))/(dr(1)**2)
 #defineMacro u ## rssttt2Macro(i1,i2,i3) (-u ## ssttt2Macro(i1-1,i2,i3)+u ## ssttt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrssttt2Macro(i1,i2,i3) (u ## ssttt2Macro(i1-1,i2,i3)-2.*u ## ssttt2Macro(i1,i2,i3)+u ## ssttt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrssttt2Macro(i1,i2,i3) (-u ## ssttt2Macro(i1-2,i2,i3)+2.*u ## ssttt2Macro(i1-1,i2,i3)-2.*u ## ssttt2Macro(i1+1,i2,i3)+u ## ssttt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssttt2Macro(i1,i2,i3) (u ## ssttt2Macro(i1-2,i2,i3)-4.*u ## ssttt2Macro(i1-1,i2,i3)+6.*u ## ssttt2Macro(i1,i2,i3)-4.*u ## ssttt2Macro(i1+1,i2,i3)+u ## ssttt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## sssttt2Macro(i1,i2,i3) (-u ## ttt2Macro(i1,i2-2,i3)+2.*u ## ttt2Macro(i1,i2-1,i3)-2.*u ## ttt2Macro(i1,i2+1,i3)+u ## ttt2Macro(i1,i2+2,i3))/(2.*dr(1)**3)
 #defineMacro u ## rsssttt2Macro(i1,i2,i3) (-u ## sssttt2Macro(i1-1,i2,i3)+u ## sssttt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrsssttt2Macro(i1,i2,i3) (u ## sssttt2Macro(i1-1,i2,i3)-2.*u ## sssttt2Macro(i1,i2,i3)+u ## sssttt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrsssttt2Macro(i1,i2,i3) (-u ## sssttt2Macro(i1-2,i2,i3)+2.*u ## sssttt2Macro(i1-1,i2,i3)-2.*u ## sssttt2Macro(i1+1,i2,i3)+u ## sssttt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsssttt2Macro(i1,i2,i3) (u ## sssttt2Macro(i1-2,i2,i3)-4.*u ## sssttt2Macro(i1-1,i2,i3)+6.*u ## sssttt2Macro(i1,i2,i3)-4.*u ## sssttt2Macro(i1+1,i2,i3)+u ## sssttt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## ssssttt2Macro(i1,i2,i3) (u ## ttt2Macro(i1,i2-2,i3)-4.*u ## ttt2Macro(i1,i2-1,i3)+6.*u ## ttt2Macro(i1,i2,i3)-4.*u ## ttt2Macro(i1,i2+1,i3)+u ## ttt2Macro(i1,i2+2,i3))/(dr(1)**4)
 #defineMacro u ## rssssttt2Macro(i1,i2,i3) (-u ## ssssttt2Macro(i1-1,i2,i3)+u ## ssssttt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrssssttt2Macro(i1,i2,i3) (u ## ssssttt2Macro(i1-1,i2,i3)-2.*u ## ssssttt2Macro(i1,i2,i3)+u ## ssssttt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrssssttt2Macro(i1,i2,i3) (-u ## ssssttt2Macro(i1-2,i2,i3)+2.*u ## ssssttt2Macro(i1-1,i2,i3)-2.*u ## ssssttt2Macro(i1+1,i2,i3)+u ## ssssttt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssssttt2Macro(i1,i2,i3) (u ## ssssttt2Macro(i1-2,i2,i3)-4.*u ## ssssttt2Macro(i1-1,i2,i3)+6.*u ## ssssttt2Macro(i1,i2,i3)-4.*u ## ssssttt2Macro(i1+1,i2,i3)+u ## ssssttt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## tttt2Macro(i1,i2,i3) (u(i1,i2,i3-2)-4.*u(i1,i2,i3-1)+6.*u(i1,i2,i3)-4.*u(i1,i2,i3+1)+u(i1,i2,i3+2))/(dr(2)**4)
 #defineMacro u ## rtttt2Macro(i1,i2,i3) (-u ## tttt2Macro(i1-1,i2,i3)+u ## tttt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrtttt2Macro(i1,i2,i3) (u ## tttt2Macro(i1-1,i2,i3)-2.*u ## tttt2Macro(i1,i2,i3)+u ## tttt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrtttt2Macro(i1,i2,i3) (-u ## tttt2Macro(i1-2,i2,i3)+2.*u ## tttt2Macro(i1-1,i2,i3)-2.*u ## tttt2Macro(i1+1,i2,i3)+u ## tttt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrtttt2Macro(i1,i2,i3) (u ## tttt2Macro(i1-2,i2,i3)-4.*u ## tttt2Macro(i1-1,i2,i3)+6.*u ## tttt2Macro(i1,i2,i3)-4.*u ## tttt2Macro(i1+1,i2,i3)+u ## tttt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## stttt2Macro(i1,i2,i3) (-u ## tttt2Macro(i1,i2-1,i3)+u ## tttt2Macro(i1,i2+1,i3))/(2.*dr(1))
 #defineMacro u ## rstttt2Macro(i1,i2,i3) (-u ## stttt2Macro(i1-1,i2,i3)+u ## stttt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrstttt2Macro(i1,i2,i3) (u ## stttt2Macro(i1-1,i2,i3)-2.*u ## stttt2Macro(i1,i2,i3)+u ## stttt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrstttt2Macro(i1,i2,i3) (-u ## stttt2Macro(i1-2,i2,i3)+2.*u ## stttt2Macro(i1-1,i2,i3)-2.*u ## stttt2Macro(i1+1,i2,i3)+u ## stttt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrstttt2Macro(i1,i2,i3) (u ## stttt2Macro(i1-2,i2,i3)-4.*u ## stttt2Macro(i1-1,i2,i3)+6.*u ## stttt2Macro(i1,i2,i3)-4.*u ## stttt2Macro(i1+1,i2,i3)+u ## stttt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## sstttt2Macro(i1,i2,i3) (u ## tttt2Macro(i1,i2-1,i3)-2.*u ## tttt2Macro(i1,i2,i3)+u ## tttt2Macro(i1,i2+1,i3))/(dr(1)**2)
 #defineMacro u ## rsstttt2Macro(i1,i2,i3) (-u ## sstttt2Macro(i1-1,i2,i3)+u ## sstttt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrsstttt2Macro(i1,i2,i3) (u ## sstttt2Macro(i1-1,i2,i3)-2.*u ## sstttt2Macro(i1,i2,i3)+u ## sstttt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrsstttt2Macro(i1,i2,i3) (-u ## sstttt2Macro(i1-2,i2,i3)+2.*u ## sstttt2Macro(i1-1,i2,i3)-2.*u ## sstttt2Macro(i1+1,i2,i3)+u ## sstttt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsstttt2Macro(i1,i2,i3) (u ## sstttt2Macro(i1-2,i2,i3)-4.*u ## sstttt2Macro(i1-1,i2,i3)+6.*u ## sstttt2Macro(i1,i2,i3)-4.*u ## sstttt2Macro(i1+1,i2,i3)+u ## sstttt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## ssstttt2Macro(i1,i2,i3) (-u ## tttt2Macro(i1,i2-2,i3)+2.*u ## tttt2Macro(i1,i2-1,i3)-2.*u ## tttt2Macro(i1,i2+1,i3)+u ## tttt2Macro(i1,i2+2,i3))/(2.*dr(1)**3)
 #defineMacro u ## rssstttt2Macro(i1,i2,i3) (-u ## ssstttt2Macro(i1-1,i2,i3)+u ## ssstttt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrssstttt2Macro(i1,i2,i3) (u ## ssstttt2Macro(i1-1,i2,i3)-2.*u ## ssstttt2Macro(i1,i2,i3)+u ## ssstttt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrssstttt2Macro(i1,i2,i3) (-u ## ssstttt2Macro(i1-2,i2,i3)+2.*u ## ssstttt2Macro(i1-1,i2,i3)-2.*u ## ssstttt2Macro(i1+1,i2,i3)+u ## ssstttt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssstttt2Macro(i1,i2,i3) (u ## ssstttt2Macro(i1-2,i2,i3)-4.*u ## ssstttt2Macro(i1-1,i2,i3)+6.*u ## ssstttt2Macro(i1,i2,i3)-4.*u ## ssstttt2Macro(i1+1,i2,i3)+u ## ssstttt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## sssstttt2Macro(i1,i2,i3) (u ## tttt2Macro(i1,i2-2,i3)-4.*u ## tttt2Macro(i1,i2-1,i3)+6.*u ## tttt2Macro(i1,i2,i3)-4.*u ## tttt2Macro(i1,i2+1,i3)+u ## tttt2Macro(i1,i2+2,i3))/(dr(1)**4)
 #defineMacro u ## rsssstttt2Macro(i1,i2,i3) (-u ## sssstttt2Macro(i1-1,i2,i3)+u ## sssstttt2Macro(i1+1,i2,i3))/(2.*dr(0))
 #defineMacro u ## rrsssstttt2Macro(i1,i2,i3) (u ## sssstttt2Macro(i1-1,i2,i3)-2.*u ## sssstttt2Macro(i1,i2,i3)+u ## sssstttt2Macro(i1+1,i2,i3))/(dr(0)**2)
 #defineMacro u ## rrrsssstttt2Macro(i1,i2,i3) (-u ## sssstttt2Macro(i1-2,i2,i3)+2.*u ## sssstttt2Macro(i1-1,i2,i3)-2.*u ## sssstttt2Macro(i1+1,i2,i3)+u ## sssstttt2Macro(i1+2,i2,i3))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsssstttt2Macro(i1,i2,i3) (u ## sssstttt2Macro(i1-2,i2,i3)-4.*u ## sssstttt2Macro(i1-1,i2,i3)+6.*u ## sssstttt2Macro(i1,i2,i3)-4.*u ## sssstttt2Macro(i1+1,i2,i3)+u ## sssstttt2Macro(i1+2,i2,i3))/(dr(0)**4)
 #defineMacro u ## r4Macro(i1,i2,i3) (u(i1-2,i2,i3)-8.*u(i1-1,i2,i3)+8.*u(i1+1,i2,i3)-u(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## rr4Macro(i1,i2,i3) (-u(i1-2,i2,i3)+16.*u(i1-1,i2,i3)-30.*u(i1,i2,i3)+16.*u(i1+1,i2,i3)-u(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## s4Macro(i1,i2,i3) (u(i1,i2-2,i3)-8.*u(i1,i2-1,i3)+8.*u(i1,i2+1,i3)-u(i1,i2+2,i3))/(12.*dr(1))
 #defineMacro u ## rs4Macro(i1,i2,i3) (u ## s4Macro(i1-2,i2,i3)-8.*u ## s4Macro(i1-1,i2,i3)+8.*u ## s4Macro(i1+1,i2,i3)-u ## s4Macro(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## rrs4Macro(i1,i2,i3) (-u ## s4Macro(i1-2,i2,i3)+16.*u ## s4Macro(i1-1,i2,i3)-30.*u ## s4Macro(i1,i2,i3)+16.*u ## s4Macro(i1+1,i2,i3)-u ## s4Macro(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## ss4Macro(i1,i2,i3) (-u(i1,i2-2,i3)+16.*u(i1,i2-1,i3)-30.*u(i1,i2,i3)+16.*u(i1,i2+1,i3)-u(i1,i2+2,i3))/(12.*dr(1)**2)
 #defineMacro u ## rss4Macro(i1,i2,i3) (u ## ss4Macro(i1-2,i2,i3)-8.*u ## ss4Macro(i1-1,i2,i3)+8.*u ## ss4Macro(i1+1,i2,i3)-u ## ss4Macro(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## rrss4Macro(i1,i2,i3) (-u ## ss4Macro(i1-2,i2,i3)+16.*u ## ss4Macro(i1-1,i2,i3)-30.*u ## ss4Macro(i1,i2,i3)+16.*u ## ss4Macro(i1+1,i2,i3)-u ## ss4Macro(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## t4Macro(i1,i2,i3) (u(i1,i2,i3-2)-8.*u(i1,i2,i3-1)+8.*u(i1,i2,i3+1)-u(i1,i2,i3+2))/(12.*dr(2))
 #defineMacro u ## rt4Macro(i1,i2,i3) (u ## t4Macro(i1-2,i2,i3)-8.*u ## t4Macro(i1-1,i2,i3)+8.*u ## t4Macro(i1+1,i2,i3)-u ## t4Macro(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## rrt4Macro(i1,i2,i3) (-u ## t4Macro(i1-2,i2,i3)+16.*u ## t4Macro(i1-1,i2,i3)-30.*u ## t4Macro(i1,i2,i3)+16.*u ## t4Macro(i1+1,i2,i3)-u ## t4Macro(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## st4Macro(i1,i2,i3) (u ## t4Macro(i1,i2-2,i3)-8.*u ## t4Macro(i1,i2-1,i3)+8.*u ## t4Macro(i1,i2+1,i3)-u ## t4Macro(i1,i2+2,i3))/(12.*dr(1))
 #defineMacro u ## rst4Macro(i1,i2,i3) (u ## st4Macro(i1-2,i2,i3)-8.*u ## st4Macro(i1-1,i2,i3)+8.*u ## st4Macro(i1+1,i2,i3)-u ## st4Macro(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## rrst4Macro(i1,i2,i3) (-u ## st4Macro(i1-2,i2,i3)+16.*u ## st4Macro(i1-1,i2,i3)-30.*u ## st4Macro(i1,i2,i3)+16.*u ## st4Macro(i1+1,i2,i3)-u ## st4Macro(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## sst4Macro(i1,i2,i3) (-u ## t4Macro(i1,i2-2,i3)+16.*u ## t4Macro(i1,i2-1,i3)-30.*u ## t4Macro(i1,i2,i3)+16.*u ## t4Macro(i1,i2+1,i3)-u ## t4Macro(i1,i2+2,i3))/(12.*dr(1)**2)
 #defineMacro u ## rsst4Macro(i1,i2,i3) (u ## sst4Macro(i1-2,i2,i3)-8.*u ## sst4Macro(i1-1,i2,i3)+8.*u ## sst4Macro(i1+1,i2,i3)-u ## sst4Macro(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## rrsst4Macro(i1,i2,i3) (-u ## sst4Macro(i1-2,i2,i3)+16.*u ## sst4Macro(i1-1,i2,i3)-30.*u ## sst4Macro(i1,i2,i3)+16.*u ## sst4Macro(i1+1,i2,i3)-u ## sst4Macro(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## tt4Macro(i1,i2,i3) (-u(i1,i2,i3-2)+16.*u(i1,i2,i3-1)-30.*u(i1,i2,i3)+16.*u(i1,i2,i3+1)-u(i1,i2,i3+2))/(12.*dr(2)**2)
 #defineMacro u ## rtt4Macro(i1,i2,i3) (u ## tt4Macro(i1-2,i2,i3)-8.*u ## tt4Macro(i1-1,i2,i3)+8.*u ## tt4Macro(i1+1,i2,i3)-u ## tt4Macro(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## rrtt4Macro(i1,i2,i3) (-u ## tt4Macro(i1-2,i2,i3)+16.*u ## tt4Macro(i1-1,i2,i3)-30.*u ## tt4Macro(i1,i2,i3)+16.*u ## tt4Macro(i1+1,i2,i3)-u ## tt4Macro(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## stt4Macro(i1,i2,i3) (u ## tt4Macro(i1,i2-2,i3)-8.*u ## tt4Macro(i1,i2-1,i3)+8.*u ## tt4Macro(i1,i2+1,i3)-u ## tt4Macro(i1,i2+2,i3))/(12.*dr(1))
 #defineMacro u ## rstt4Macro(i1,i2,i3) (u ## stt4Macro(i1-2,i2,i3)-8.*u ## stt4Macro(i1-1,i2,i3)+8.*u ## stt4Macro(i1+1,i2,i3)-u ## stt4Macro(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## rrstt4Macro(i1,i2,i3) (-u ## stt4Macro(i1-2,i2,i3)+16.*u ## stt4Macro(i1-1,i2,i3)-30.*u ## stt4Macro(i1,i2,i3)+16.*u ## stt4Macro(i1+1,i2,i3)-u ## stt4Macro(i1+2,i2,i3))/(12.*dr(0)**2)
 #defineMacro u ## sstt4Macro(i1,i2,i3) (-u ## tt4Macro(i1,i2-2,i3)+16.*u ## tt4Macro(i1,i2-1,i3)-30.*u ## tt4Macro(i1,i2,i3)+16.*u ## tt4Macro(i1,i2+1,i3)-u ## tt4Macro(i1,i2+2,i3))/(12.*dr(1)**2)
 #defineMacro u ## rsstt4Macro(i1,i2,i3) (u ## sstt4Macro(i1-2,i2,i3)-8.*u ## sstt4Macro(i1-1,i2,i3)+8.*u ## sstt4Macro(i1+1,i2,i3)-u ## sstt4Macro(i1+2,i2,i3))/(12.*dr(0))
 #defineMacro u ## rrsstt4Macro(i1,i2,i3) (-u ## sstt4Macro(i1-2,i2,i3)+16.*u ## sstt4Macro(i1-1,i2,i3)-30.*u ## sstt4Macro(i1,i2,i3)+16.*u ## sstt4Macro(i1+1,i2,i3)-u ## sstt4Macro(i1+2,i2,i3))/(12.*dr(0)**2)
#endMacro

#beginMacro defineThirdParameticDerivativesComponents1(u)
 #defineMacro u ## r2Macro(i1,i2,i3,m) (-u(i1-1,i2,i3,m)+u(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rr2Macro(i1,i2,i3,m) (u(i1-1,i2,i3,m)-2.*u(i1,i2,i3,m)+u(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrr2Macro(i1,i2,i3,m) (-u(i1-2,i2,i3,m)+2.*u(i1-1,i2,i3,m)-2.*u(i1+1,i2,i3,m)+u(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrr2Macro(i1,i2,i3,m) (u(i1-2,i2,i3,m)-4.*u(i1-1,i2,i3,m)+6.*u(i1,i2,i3,m)-4.*u(i1+1,i2,i3,m)+u(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## s2Macro(i1,i2,i3,m) (-u(i1,i2-1,i3,m)+u(i1,i2+1,i3,m))/(2.*dr(1))
 #defineMacro u ## rs2Macro(i1,i2,i3,m) (-u ## s2Macro(i1-1,i2,i3,m)+u ## s2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrs2Macro(i1,i2,i3,m) (u ## s2Macro(i1-1,i2,i3,m)-2.*u ## s2Macro(i1,i2,i3,m)+u ## s2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrs2Macro(i1,i2,i3,m) (-u ## s2Macro(i1-2,i2,i3,m)+2.*u ## s2Macro(i1-1,i2,i3,m)-2.*u ## s2Macro(i1+1,i2,i3,m)+u ## s2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrs2Macro(i1,i2,i3,m) (u ## s2Macro(i1-2,i2,i3,m)-4.*u ## s2Macro(i1-1,i2,i3,m)+6.*u ## s2Macro(i1,i2,i3,m)-4.*u ## s2Macro(i1+1,i2,i3,m)+u ## s2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## ss2Macro(i1,i2,i3,m) (u(i1,i2-1,i3,m)-2.*u(i1,i2,i3,m)+u(i1,i2+1,i3,m))/(dr(1)**2)
 #defineMacro u ## rss2Macro(i1,i2,i3,m) (-u ## ss2Macro(i1-1,i2,i3,m)+u ## ss2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrss2Macro(i1,i2,i3,m) (u ## ss2Macro(i1-1,i2,i3,m)-2.*u ## ss2Macro(i1,i2,i3,m)+u ## ss2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrss2Macro(i1,i2,i3,m) (-u ## ss2Macro(i1-2,i2,i3,m)+2.*u ## ss2Macro(i1-1,i2,i3,m)-2.*u ## ss2Macro(i1+1,i2,i3,m)+u ## ss2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrss2Macro(i1,i2,i3,m) (u ## ss2Macro(i1-2,i2,i3,m)-4.*u ## ss2Macro(i1-1,i2,i3,m)+6.*u ## ss2Macro(i1,i2,i3,m)-4.*u ## ss2Macro(i1+1,i2,i3,m)+u ## ss2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## sss2Macro(i1,i2,i3,m) (-u(i1,i2-2,i3,m)+2.*u(i1,i2-1,i3,m)-2.*u(i1,i2+1,i3,m)+u(i1,i2+2,i3,m))/(2.*dr(1)**3)
 #defineMacro u ## rsss2Macro(i1,i2,i3,m) (-u ## sss2Macro(i1-1,i2,i3,m)+u ## sss2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrsss2Macro(i1,i2,i3,m) (u ## sss2Macro(i1-1,i2,i3,m)-2.*u ## sss2Macro(i1,i2,i3,m)+u ## sss2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrsss2Macro(i1,i2,i3,m) (-u ## sss2Macro(i1-2,i2,i3,m)+2.*u ## sss2Macro(i1-1,i2,i3,m)-2.*u ## sss2Macro(i1+1,i2,i3,m)+u ## sss2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsss2Macro(i1,i2,i3,m) (u ## sss2Macro(i1-2,i2,i3,m)-4.*u ## sss2Macro(i1-1,i2,i3,m)+6.*u ## sss2Macro(i1,i2,i3,m)-4.*u ## sss2Macro(i1+1,i2,i3,m)+u ## sss2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## ssss2Macro(i1,i2,i3,m) (u(i1,i2-2,i3,m)-4.*u(i1,i2-1,i3,m)+6.*u(i1,i2,i3,m)-4.*u(i1,i2+1,i3,m)+u(i1,i2+2,i3,m))/(dr(1)**4)
 #defineMacro u ## rssss2Macro(i1,i2,i3,m) (-u ## ssss2Macro(i1-1,i2,i3,m)+u ## ssss2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrssss2Macro(i1,i2,i3,m) (u ## ssss2Macro(i1-1,i2,i3,m)-2.*u ## ssss2Macro(i1,i2,i3,m)+u ## ssss2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrssss2Macro(i1,i2,i3,m) (-u ## ssss2Macro(i1-2,i2,i3,m)+2.*u ## ssss2Macro(i1-1,i2,i3,m)-2.*u ## ssss2Macro(i1+1,i2,i3,m)+u ## ssss2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssss2Macro(i1,i2,i3,m) (u ## ssss2Macro(i1-2,i2,i3,m)-4.*u ## ssss2Macro(i1-1,i2,i3,m)+6.*u ## ssss2Macro(i1,i2,i3,m)-4.*u ## ssss2Macro(i1+1,i2,i3,m)+u ## ssss2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## t2Macro(i1,i2,i3,m) (-u(i1,i2,i3-1,m)+u(i1,i2,i3+1,m))/(2.*dr(2))
 #defineMacro u ## rt2Macro(i1,i2,i3,m) (-u ## t2Macro(i1-1,i2,i3,m)+u ## t2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrt2Macro(i1,i2,i3,m) (u ## t2Macro(i1-1,i2,i3,m)-2.*u ## t2Macro(i1,i2,i3,m)+u ## t2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrt2Macro(i1,i2,i3,m) (-u ## t2Macro(i1-2,i2,i3,m)+2.*u ## t2Macro(i1-1,i2,i3,m)-2.*u ## t2Macro(i1+1,i2,i3,m)+u ## t2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrt2Macro(i1,i2,i3,m) (u ## t2Macro(i1-2,i2,i3,m)-4.*u ## t2Macro(i1-1,i2,i3,m)+6.*u ## t2Macro(i1,i2,i3,m)-4.*u ## t2Macro(i1+1,i2,i3,m)+u ## t2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## st2Macro(i1,i2,i3,m) (-u ## t2Macro(i1,i2-1,i3,m)+u ## t2Macro(i1,i2+1,i3,m))/(2.*dr(1))
 #defineMacro u ## rst2Macro(i1,i2,i3,m) (-u ## st2Macro(i1-1,i2,i3,m)+u ## st2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrst2Macro(i1,i2,i3,m) (u ## st2Macro(i1-1,i2,i3,m)-2.*u ## st2Macro(i1,i2,i3,m)+u ## st2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrst2Macro(i1,i2,i3,m) (-u ## st2Macro(i1-2,i2,i3,m)+2.*u ## st2Macro(i1-1,i2,i3,m)-2.*u ## st2Macro(i1+1,i2,i3,m)+u ## st2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrst2Macro(i1,i2,i3,m) (u ## st2Macro(i1-2,i2,i3,m)-4.*u ## st2Macro(i1-1,i2,i3,m)+6.*u ## st2Macro(i1,i2,i3,m)-4.*u ## st2Macro(i1+1,i2,i3,m)+u ## st2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## sst2Macro(i1,i2,i3,m) (u ## t2Macro(i1,i2-1,i3,m)-2.*u ## t2Macro(i1,i2,i3,m)+u ## t2Macro(i1,i2+1,i3,m))/(dr(1)**2)
 #defineMacro u ## rsst2Macro(i1,i2,i3,m) (-u ## sst2Macro(i1-1,i2,i3,m)+u ## sst2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrsst2Macro(i1,i2,i3,m) (u ## sst2Macro(i1-1,i2,i3,m)-2.*u ## sst2Macro(i1,i2,i3,m)+u ## sst2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrsst2Macro(i1,i2,i3,m) (-u ## sst2Macro(i1-2,i2,i3,m)+2.*u ## sst2Macro(i1-1,i2,i3,m)-2.*u ## sst2Macro(i1+1,i2,i3,m)+u ## sst2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsst2Macro(i1,i2,i3,m) (u ## sst2Macro(i1-2,i2,i3,m)-4.*u ## sst2Macro(i1-1,i2,i3,m)+6.*u ## sst2Macro(i1,i2,i3,m)-4.*u ## sst2Macro(i1+1,i2,i3,m)+u ## sst2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## ssst2Macro(i1,i2,i3,m) (-u ## t2Macro(i1,i2-2,i3,m)+2.*u ## t2Macro(i1,i2-1,i3,m)-2.*u ## t2Macro(i1,i2+1,i3,m)+u ## t2Macro(i1,i2+2,i3,m))/(2.*dr(1)**3)
 #defineMacro u ## rssst2Macro(i1,i2,i3,m) (-u ## ssst2Macro(i1-1,i2,i3,m)+u ## ssst2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrssst2Macro(i1,i2,i3,m) (u ## ssst2Macro(i1-1,i2,i3,m)-2.*u ## ssst2Macro(i1,i2,i3,m)+u ## ssst2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrssst2Macro(i1,i2,i3,m) (-u ## ssst2Macro(i1-2,i2,i3,m)+2.*u ## ssst2Macro(i1-1,i2,i3,m)-2.*u ## ssst2Macro(i1+1,i2,i3,m)+u ## ssst2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssst2Macro(i1,i2,i3,m) (u ## ssst2Macro(i1-2,i2,i3,m)-4.*u ## ssst2Macro(i1-1,i2,i3,m)+6.*u ## ssst2Macro(i1,i2,i3,m)-4.*u ## ssst2Macro(i1+1,i2,i3,m)+u ## ssst2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## sssst2Macro(i1,i2,i3,m) (u ## t2Macro(i1,i2-2,i3,m)-4.*u ## t2Macro(i1,i2-1,i3,m)+6.*u ## t2Macro(i1,i2,i3,m)-4.*u ## t2Macro(i1,i2+1,i3,m)+u ## t2Macro(i1,i2+2,i3,m))/(dr(1)**4)
 #defineMacro u ## rsssst2Macro(i1,i2,i3,m) (-u ## sssst2Macro(i1-1,i2,i3,m)+u ## sssst2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrsssst2Macro(i1,i2,i3,m) (u ## sssst2Macro(i1-1,i2,i3,m)-2.*u ## sssst2Macro(i1,i2,i3,m)+u ## sssst2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrsssst2Macro(i1,i2,i3,m) (-u ## sssst2Macro(i1-2,i2,i3,m)+2.*u ## sssst2Macro(i1-1,i2,i3,m)-2.*u ## sssst2Macro(i1+1,i2,i3,m)+u ## sssst2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsssst2Macro(i1,i2,i3,m) (u ## sssst2Macro(i1-2,i2,i3,m)-4.*u ## sssst2Macro(i1-1,i2,i3,m)+6.*u ## sssst2Macro(i1,i2,i3,m)-4.*u ## sssst2Macro(i1+1,i2,i3,m)+u ## sssst2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## tt2Macro(i1,i2,i3,m) (u(i1,i2,i3-1,m)-2.*u(i1,i2,i3,m)+u(i1,i2,i3+1,m))/(dr(2)**2)
 #defineMacro u ## rtt2Macro(i1,i2,i3,m) (-u ## tt2Macro(i1-1,i2,i3,m)+u ## tt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrtt2Macro(i1,i2,i3,m) (u ## tt2Macro(i1-1,i2,i3,m)-2.*u ## tt2Macro(i1,i2,i3,m)+u ## tt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrtt2Macro(i1,i2,i3,m) (-u ## tt2Macro(i1-2,i2,i3,m)+2.*u ## tt2Macro(i1-1,i2,i3,m)-2.*u ## tt2Macro(i1+1,i2,i3,m)+u ## tt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrtt2Macro(i1,i2,i3,m) (u ## tt2Macro(i1-2,i2,i3,m)-4.*u ## tt2Macro(i1-1,i2,i3,m)+6.*u ## tt2Macro(i1,i2,i3,m)-4.*u ## tt2Macro(i1+1,i2,i3,m)+u ## tt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## stt2Macro(i1,i2,i3,m) (-u ## tt2Macro(i1,i2-1,i3,m)+u ## tt2Macro(i1,i2+1,i3,m))/(2.*dr(1))
 #defineMacro u ## rstt2Macro(i1,i2,i3,m) (-u ## stt2Macro(i1-1,i2,i3,m)+u ## stt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrstt2Macro(i1,i2,i3,m) (u ## stt2Macro(i1-1,i2,i3,m)-2.*u ## stt2Macro(i1,i2,i3,m)+u ## stt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrstt2Macro(i1,i2,i3,m) (-u ## stt2Macro(i1-2,i2,i3,m)+2.*u ## stt2Macro(i1-1,i2,i3,m)-2.*u ## stt2Macro(i1+1,i2,i3,m)+u ## stt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrstt2Macro(i1,i2,i3,m) (u ## stt2Macro(i1-2,i2,i3,m)-4.*u ## stt2Macro(i1-1,i2,i3,m)+6.*u ## stt2Macro(i1,i2,i3,m)-4.*u ## stt2Macro(i1+1,i2,i3,m)+u ## stt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## sstt2Macro(i1,i2,i3,m) (u ## tt2Macro(i1,i2-1,i3,m)-2.*u ## tt2Macro(i1,i2,i3,m)+u ## tt2Macro(i1,i2+1,i3,m))/(dr(1)**2)
 #defineMacro u ## rsstt2Macro(i1,i2,i3,m) (-u ## sstt2Macro(i1-1,i2,i3,m)+u ## sstt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrsstt2Macro(i1,i2,i3,m) (u ## sstt2Macro(i1-1,i2,i3,m)-2.*u ## sstt2Macro(i1,i2,i3,m)+u ## sstt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrsstt2Macro(i1,i2,i3,m) (-u ## sstt2Macro(i1-2,i2,i3,m)+2.*u ## sstt2Macro(i1-1,i2,i3,m)-2.*u ## sstt2Macro(i1+1,i2,i3,m)+u ## sstt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsstt2Macro(i1,i2,i3,m) (u ## sstt2Macro(i1-2,i2,i3,m)-4.*u ## sstt2Macro(i1-1,i2,i3,m)+6.*u ## sstt2Macro(i1,i2,i3,m)-4.*u ## sstt2Macro(i1+1,i2,i3,m)+u ## sstt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## ssstt2Macro(i1,i2,i3,m) (-u ## tt2Macro(i1,i2-2,i3,m)+2.*u ## tt2Macro(i1,i2-1,i3,m)-2.*u ## tt2Macro(i1,i2+1,i3,m)+u ## tt2Macro(i1,i2+2,i3,m))/(2.*dr(1)**3)
 #defineMacro u ## rssstt2Macro(i1,i2,i3,m) (-u ## ssstt2Macro(i1-1,i2,i3,m)+u ## ssstt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrssstt2Macro(i1,i2,i3,m) (u ## ssstt2Macro(i1-1,i2,i3,m)-2.*u ## ssstt2Macro(i1,i2,i3,m)+u ## ssstt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrssstt2Macro(i1,i2,i3,m) (-u ## ssstt2Macro(i1-2,i2,i3,m)+2.*u ## ssstt2Macro(i1-1,i2,i3,m)-2.*u ## ssstt2Macro(i1+1,i2,i3,m)+u ## ssstt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssstt2Macro(i1,i2,i3,m) (u ## ssstt2Macro(i1-2,i2,i3,m)-4.*u ## ssstt2Macro(i1-1,i2,i3,m)+6.*u ## ssstt2Macro(i1,i2,i3,m)-4.*u ## ssstt2Macro(i1+1,i2,i3,m)+u ## ssstt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## sssstt2Macro(i1,i2,i3,m) (u ## tt2Macro(i1,i2-2,i3,m)-4.*u ## tt2Macro(i1,i2-1,i3,m)+6.*u ## tt2Macro(i1,i2,i3,m)-4.*u ## tt2Macro(i1,i2+1,i3,m)+u ## tt2Macro(i1,i2+2,i3,m))/(dr(1)**4)
 #defineMacro u ## rsssstt2Macro(i1,i2,i3,m) (-u ## sssstt2Macro(i1-1,i2,i3,m)+u ## sssstt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrsssstt2Macro(i1,i2,i3,m) (u ## sssstt2Macro(i1-1,i2,i3,m)-2.*u ## sssstt2Macro(i1,i2,i3,m)+u ## sssstt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrsssstt2Macro(i1,i2,i3,m) (-u ## sssstt2Macro(i1-2,i2,i3,m)+2.*u ## sssstt2Macro(i1-1,i2,i3,m)-2.*u ## sssstt2Macro(i1+1,i2,i3,m)+u ## sssstt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsssstt2Macro(i1,i2,i3,m) (u ## sssstt2Macro(i1-2,i2,i3,m)-4.*u ## sssstt2Macro(i1-1,i2,i3,m)+6.*u ## sssstt2Macro(i1,i2,i3,m)-4.*u ## sssstt2Macro(i1+1,i2,i3,m)+u ## sssstt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## ttt2Macro(i1,i2,i3,m) (-u(i1,i2,i3-2,m)+2.*u(i1,i2,i3-1,m)-2.*u(i1,i2,i3+1,m)+u(i1,i2,i3+2,m))/(2.*dr(2)**3)
 #defineMacro u ## rttt2Macro(i1,i2,i3,m) (-u ## ttt2Macro(i1-1,i2,i3,m)+u ## ttt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrttt2Macro(i1,i2,i3,m) (u ## ttt2Macro(i1-1,i2,i3,m)-2.*u ## ttt2Macro(i1,i2,i3,m)+u ## ttt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrttt2Macro(i1,i2,i3,m) (-u ## ttt2Macro(i1-2,i2,i3,m)+2.*u ## ttt2Macro(i1-1,i2,i3,m)-2.*u ## ttt2Macro(i1+1,i2,i3,m)+u ## ttt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrttt2Macro(i1,i2,i3,m) (u ## ttt2Macro(i1-2,i2,i3,m)-4.*u ## ttt2Macro(i1-1,i2,i3,m)+6.*u ## ttt2Macro(i1,i2,i3,m)-4.*u ## ttt2Macro(i1+1,i2,i3,m)+u ## ttt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## sttt2Macro(i1,i2,i3,m) (-u ## ttt2Macro(i1,i2-1,i3,m)+u ## ttt2Macro(i1,i2+1,i3,m))/(2.*dr(1))
 #defineMacro u ## rsttt2Macro(i1,i2,i3,m) (-u ## sttt2Macro(i1-1,i2,i3,m)+u ## sttt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrsttt2Macro(i1,i2,i3,m) (u ## sttt2Macro(i1-1,i2,i3,m)-2.*u ## sttt2Macro(i1,i2,i3,m)+u ## sttt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrsttt2Macro(i1,i2,i3,m) (-u ## sttt2Macro(i1-2,i2,i3,m)+2.*u ## sttt2Macro(i1-1,i2,i3,m)-2.*u ## sttt2Macro(i1+1,i2,i3,m)+u ## sttt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsttt2Macro(i1,i2,i3,m) (u ## sttt2Macro(i1-2,i2,i3,m)-4.*u ## sttt2Macro(i1-1,i2,i3,m)+6.*u ## sttt2Macro(i1,i2,i3,m)-4.*u ## sttt2Macro(i1+1,i2,i3,m)+u ## sttt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## ssttt2Macro(i1,i2,i3,m) (u ## ttt2Macro(i1,i2-1,i3,m)-2.*u ## ttt2Macro(i1,i2,i3,m)+u ## ttt2Macro(i1,i2+1,i3,m))/(dr(1)**2)
 #defineMacro u ## rssttt2Macro(i1,i2,i3,m) (-u ## ssttt2Macro(i1-1,i2,i3,m)+u ## ssttt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrssttt2Macro(i1,i2,i3,m) (u ## ssttt2Macro(i1-1,i2,i3,m)-2.*u ## ssttt2Macro(i1,i2,i3,m)+u ## ssttt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrssttt2Macro(i1,i2,i3,m) (-u ## ssttt2Macro(i1-2,i2,i3,m)+2.*u ## ssttt2Macro(i1-1,i2,i3,m)-2.*u ## ssttt2Macro(i1+1,i2,i3,m)+u ## ssttt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssttt2Macro(i1,i2,i3,m) (u ## ssttt2Macro(i1-2,i2,i3,m)-4.*u ## ssttt2Macro(i1-1,i2,i3,m)+6.*u ## ssttt2Macro(i1,i2,i3,m)-4.*u ## ssttt2Macro(i1+1,i2,i3,m)+u ## ssttt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## sssttt2Macro(i1,i2,i3,m) (-u ## ttt2Macro(i1,i2-2,i3,m)+2.*u ## ttt2Macro(i1,i2-1,i3,m)-2.*u ## ttt2Macro(i1,i2+1,i3,m)+u ## ttt2Macro(i1,i2+2,i3,m))/(2.*dr(1)**3)
 #defineMacro u ## rsssttt2Macro(i1,i2,i3,m) (-u ## sssttt2Macro(i1-1,i2,i3,m)+u ## sssttt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrsssttt2Macro(i1,i2,i3,m) (u ## sssttt2Macro(i1-1,i2,i3,m)-2.*u ## sssttt2Macro(i1,i2,i3,m)+u ## sssttt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrsssttt2Macro(i1,i2,i3,m) (-u ## sssttt2Macro(i1-2,i2,i3,m)+2.*u ## sssttt2Macro(i1-1,i2,i3,m)-2.*u ## sssttt2Macro(i1+1,i2,i3,m)+u ## sssttt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsssttt2Macro(i1,i2,i3,m) (u ## sssttt2Macro(i1-2,i2,i3,m)-4.*u ## sssttt2Macro(i1-1,i2,i3,m)+6.*u ## sssttt2Macro(i1,i2,i3,m)-4.*u ## sssttt2Macro(i1+1,i2,i3,m)+u ## sssttt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## ssssttt2Macro(i1,i2,i3,m) (u ## ttt2Macro(i1,i2-2,i3,m)-4.*u ## ttt2Macro(i1,i2-1,i3,m)+6.*u ## ttt2Macro(i1,i2,i3,m)-4.*u ## ttt2Macro(i1,i2+1,i3,m)+u ## ttt2Macro(i1,i2+2,i3,m))/(dr(1)**4)
 #defineMacro u ## rssssttt2Macro(i1,i2,i3,m) (-u ## ssssttt2Macro(i1-1,i2,i3,m)+u ## ssssttt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrssssttt2Macro(i1,i2,i3,m) (u ## ssssttt2Macro(i1-1,i2,i3,m)-2.*u ## ssssttt2Macro(i1,i2,i3,m)+u ## ssssttt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrssssttt2Macro(i1,i2,i3,m) (-u ## ssssttt2Macro(i1-2,i2,i3,m)+2.*u ## ssssttt2Macro(i1-1,i2,i3,m)-2.*u ## ssssttt2Macro(i1+1,i2,i3,m)+u ## ssssttt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssssttt2Macro(i1,i2,i3,m) (u ## ssssttt2Macro(i1-2,i2,i3,m)-4.*u ## ssssttt2Macro(i1-1,i2,i3,m)+6.*u ## ssssttt2Macro(i1,i2,i3,m)-4.*u ## ssssttt2Macro(i1+1,i2,i3,m)+u ## ssssttt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## tttt2Macro(i1,i2,i3,m) (u(i1,i2,i3-2,m)-4.*u(i1,i2,i3-1,m)+6.*u(i1,i2,i3,m)-4.*u(i1,i2,i3+1,m)+u(i1,i2,i3+2,m))/(dr(2)**4)
 #defineMacro u ## rtttt2Macro(i1,i2,i3,m) (-u ## tttt2Macro(i1-1,i2,i3,m)+u ## tttt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrtttt2Macro(i1,i2,i3,m) (u ## tttt2Macro(i1-1,i2,i3,m)-2.*u ## tttt2Macro(i1,i2,i3,m)+u ## tttt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrtttt2Macro(i1,i2,i3,m) (-u ## tttt2Macro(i1-2,i2,i3,m)+2.*u ## tttt2Macro(i1-1,i2,i3,m)-2.*u ## tttt2Macro(i1+1,i2,i3,m)+u ## tttt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrtttt2Macro(i1,i2,i3,m) (u ## tttt2Macro(i1-2,i2,i3,m)-4.*u ## tttt2Macro(i1-1,i2,i3,m)+6.*u ## tttt2Macro(i1,i2,i3,m)-4.*u ## tttt2Macro(i1+1,i2,i3,m)+u ## tttt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## stttt2Macro(i1,i2,i3,m) (-u ## tttt2Macro(i1,i2-1,i3,m)+u ## tttt2Macro(i1,i2+1,i3,m))/(2.*dr(1))
 #defineMacro u ## rstttt2Macro(i1,i2,i3,m) (-u ## stttt2Macro(i1-1,i2,i3,m)+u ## stttt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrstttt2Macro(i1,i2,i3,m) (u ## stttt2Macro(i1-1,i2,i3,m)-2.*u ## stttt2Macro(i1,i2,i3,m)+u ## stttt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrstttt2Macro(i1,i2,i3,m) (-u ## stttt2Macro(i1-2,i2,i3,m)+2.*u ## stttt2Macro(i1-1,i2,i3,m)-2.*u ## stttt2Macro(i1+1,i2,i3,m)+u ## stttt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrstttt2Macro(i1,i2,i3,m) (u ## stttt2Macro(i1-2,i2,i3,m)-4.*u ## stttt2Macro(i1-1,i2,i3,m)+6.*u ## stttt2Macro(i1,i2,i3,m)-4.*u ## stttt2Macro(i1+1,i2,i3,m)+u ## stttt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## sstttt2Macro(i1,i2,i3,m) (u ## tttt2Macro(i1,i2-1,i3,m)-2.*u ## tttt2Macro(i1,i2,i3,m)+u ## tttt2Macro(i1,i2+1,i3,m))/(dr(1)**2)
 #defineMacro u ## rsstttt2Macro(i1,i2,i3,m) (-u ## sstttt2Macro(i1-1,i2,i3,m)+u ## sstttt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrsstttt2Macro(i1,i2,i3,m) (u ## sstttt2Macro(i1-1,i2,i3,m)-2.*u ## sstttt2Macro(i1,i2,i3,m)+u ## sstttt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrsstttt2Macro(i1,i2,i3,m) (-u ## sstttt2Macro(i1-2,i2,i3,m)+2.*u ## sstttt2Macro(i1-1,i2,i3,m)-2.*u ## sstttt2Macro(i1+1,i2,i3,m)+u ## sstttt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsstttt2Macro(i1,i2,i3,m) (u ## sstttt2Macro(i1-2,i2,i3,m)-4.*u ## sstttt2Macro(i1-1,i2,i3,m)+6.*u ## sstttt2Macro(i1,i2,i3,m)-4.*u ## sstttt2Macro(i1+1,i2,i3,m)+u ## sstttt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## ssstttt2Macro(i1,i2,i3,m) (-u ## tttt2Macro(i1,i2-2,i3,m)+2.*u ## tttt2Macro(i1,i2-1,i3,m)-2.*u ## tttt2Macro(i1,i2+1,i3,m)+u ## tttt2Macro(i1,i2+2,i3,m))/(2.*dr(1)**3)
 #defineMacro u ## rssstttt2Macro(i1,i2,i3,m) (-u ## ssstttt2Macro(i1-1,i2,i3,m)+u ## ssstttt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrssstttt2Macro(i1,i2,i3,m) (u ## ssstttt2Macro(i1-1,i2,i3,m)-2.*u ## ssstttt2Macro(i1,i2,i3,m)+u ## ssstttt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrssstttt2Macro(i1,i2,i3,m) (-u ## ssstttt2Macro(i1-2,i2,i3,m)+2.*u ## ssstttt2Macro(i1-1,i2,i3,m)-2.*u ## ssstttt2Macro(i1+1,i2,i3,m)+u ## ssstttt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrssstttt2Macro(i1,i2,i3,m) (u ## ssstttt2Macro(i1-2,i2,i3,m)-4.*u ## ssstttt2Macro(i1-1,i2,i3,m)+6.*u ## ssstttt2Macro(i1,i2,i3,m)-4.*u ## ssstttt2Macro(i1+1,i2,i3,m)+u ## ssstttt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## sssstttt2Macro(i1,i2,i3,m) (u ## tttt2Macro(i1,i2-2,i3,m)-4.*u ## tttt2Macro(i1,i2-1,i3,m)+6.*u ## tttt2Macro(i1,i2,i3,m)-4.*u ## tttt2Macro(i1,i2+1,i3,m)+u ## tttt2Macro(i1,i2+2,i3,m))/(dr(1)**4)
 #defineMacro u ## rsssstttt2Macro(i1,i2,i3,m) (-u ## sssstttt2Macro(i1-1,i2,i3,m)+u ## sssstttt2Macro(i1+1,i2,i3,m))/(2.*dr(0))
 #defineMacro u ## rrsssstttt2Macro(i1,i2,i3,m) (u ## sssstttt2Macro(i1-1,i2,i3,m)-2.*u ## sssstttt2Macro(i1,i2,i3,m)+u ## sssstttt2Macro(i1+1,i2,i3,m))/(dr(0)**2)
 #defineMacro u ## rrrsssstttt2Macro(i1,i2,i3,m) (-u ## sssstttt2Macro(i1-2,i2,i3,m)+2.*u ## sssstttt2Macro(i1-1,i2,i3,m)-2.*u ## sssstttt2Macro(i1+1,i2,i3,m)+u ## sssstttt2Macro(i1+2,i2,i3,m))/(2.*dr(0)**3)
 #defineMacro u ## rrrrsssstttt2Macro(i1,i2,i3,m) (u ## sssstttt2Macro(i1-2,i2,i3,m)-4.*u ## sssstttt2Macro(i1-1,i2,i3,m)+6.*u ## sssstttt2Macro(i1,i2,i3,m)-4.*u ## sssstttt2Macro(i1+1,i2,i3,m)+u ## sssstttt2Macro(i1+2,i2,i3,m))/(dr(0)**4)
 #defineMacro u ## r4Macro(i1,i2,i3,m) (u(i1-2,i2,i3,m)-8.*u(i1-1,i2,i3,m)+8.*u(i1+1,i2,i3,m)-u(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## rr4Macro(i1,i2,i3,m) (-u(i1-2,i2,i3,m)+16.*u(i1-1,i2,i3,m)-30.*u(i1,i2,i3,m)+16.*u(i1+1,i2,i3,m)-u(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## s4Macro(i1,i2,i3,m) (u(i1,i2-2,i3,m)-8.*u(i1,i2-1,i3,m)+8.*u(i1,i2+1,i3,m)-u(i1,i2+2,i3,m))/(12.*dr(1))
 #defineMacro u ## rs4Macro(i1,i2,i3,m) (u ## s4Macro(i1-2,i2,i3,m)-8.*u ## s4Macro(i1-1,i2,i3,m)+8.*u ## s4Macro(i1+1,i2,i3,m)-u ## s4Macro(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## rrs4Macro(i1,i2,i3,m) (-u ## s4Macro(i1-2,i2,i3,m)+16.*u ## s4Macro(i1-1,i2,i3,m)-30.*u ## s4Macro(i1,i2,i3,m)+16.*u ## s4Macro(i1+1,i2,i3,m)-u ## s4Macro(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## ss4Macro(i1,i2,i3,m) (-u(i1,i2-2,i3,m)+16.*u(i1,i2-1,i3,m)-30.*u(i1,i2,i3,m)+16.*u(i1,i2+1,i3,m)-u(i1,i2+2,i3,m))/(12.*dr(1)**2)
 #defineMacro u ## rss4Macro(i1,i2,i3,m) (u ## ss4Macro(i1-2,i2,i3,m)-8.*u ## ss4Macro(i1-1,i2,i3,m)+8.*u ## ss4Macro(i1+1,i2,i3,m)-u ## ss4Macro(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## rrss4Macro(i1,i2,i3,m) (-u ## ss4Macro(i1-2,i2,i3,m)+16.*u ## ss4Macro(i1-1,i2,i3,m)-30.*u ## ss4Macro(i1,i2,i3,m)+16.*u ## ss4Macro(i1+1,i2,i3,m)-u ## ss4Macro(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## t4Macro(i1,i2,i3,m) (u(i1,i2,i3-2,m)-8.*u(i1,i2,i3-1,m)+8.*u(i1,i2,i3+1,m)-u(i1,i2,i3+2,m))/(12.*dr(2))
 #defineMacro u ## rt4Macro(i1,i2,i3,m) (u ## t4Macro(i1-2,i2,i3,m)-8.*u ## t4Macro(i1-1,i2,i3,m)+8.*u ## t4Macro(i1+1,i2,i3,m)-u ## t4Macro(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## rrt4Macro(i1,i2,i3,m) (-u ## t4Macro(i1-2,i2,i3,m)+16.*u ## t4Macro(i1-1,i2,i3,m)-30.*u ## t4Macro(i1,i2,i3,m)+16.*u ## t4Macro(i1+1,i2,i3,m)-u ## t4Macro(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## st4Macro(i1,i2,i3,m) (u ## t4Macro(i1,i2-2,i3,m)-8.*u ## t4Macro(i1,i2-1,i3,m)+8.*u ## t4Macro(i1,i2+1,i3,m)-u ## t4Macro(i1,i2+2,i3,m))/(12.*dr(1))
 #defineMacro u ## rst4Macro(i1,i2,i3,m) (u ## st4Macro(i1-2,i2,i3,m)-8.*u ## st4Macro(i1-1,i2,i3,m)+8.*u ## st4Macro(i1+1,i2,i3,m)-u ## st4Macro(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## rrst4Macro(i1,i2,i3,m) (-u ## st4Macro(i1-2,i2,i3,m)+16.*u ## st4Macro(i1-1,i2,i3,m)-30.*u ## st4Macro(i1,i2,i3,m)+16.*u ## st4Macro(i1+1,i2,i3,m)-u ## st4Macro(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## sst4Macro(i1,i2,i3,m) (-u ## t4Macro(i1,i2-2,i3,m)+16.*u ## t4Macro(i1,i2-1,i3,m)-30.*u ## t4Macro(i1,i2,i3,m)+16.*u ## t4Macro(i1,i2+1,i3,m)-u ## t4Macro(i1,i2+2,i3,m))/(12.*dr(1)**2)
 #defineMacro u ## rsst4Macro(i1,i2,i3,m) (u ## sst4Macro(i1-2,i2,i3,m)-8.*u ## sst4Macro(i1-1,i2,i3,m)+8.*u ## sst4Macro(i1+1,i2,i3,m)-u ## sst4Macro(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## rrsst4Macro(i1,i2,i3,m) (-u ## sst4Macro(i1-2,i2,i3,m)+16.*u ## sst4Macro(i1-1,i2,i3,m)-30.*u ## sst4Macro(i1,i2,i3,m)+16.*u ## sst4Macro(i1+1,i2,i3,m)-u ## sst4Macro(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## tt4Macro(i1,i2,i3,m) (-u(i1,i2,i3-2,m)+16.*u(i1,i2,i3-1,m)-30.*u(i1,i2,i3,m)+16.*u(i1,i2,i3+1,m)-u(i1,i2,i3+2,m))/(12.*dr(2)**2)
 #defineMacro u ## rtt4Macro(i1,i2,i3,m) (u ## tt4Macro(i1-2,i2,i3,m)-8.*u ## tt4Macro(i1-1,i2,i3,m)+8.*u ## tt4Macro(i1+1,i2,i3,m)-u ## tt4Macro(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## rrtt4Macro(i1,i2,i3,m) (-u ## tt4Macro(i1-2,i2,i3,m)+16.*u ## tt4Macro(i1-1,i2,i3,m)-30.*u ## tt4Macro(i1,i2,i3,m)+16.*u ## tt4Macro(i1+1,i2,i3,m)-u ## tt4Macro(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## stt4Macro(i1,i2,i3,m) (u ## tt4Macro(i1,i2-2,i3,m)-8.*u ## tt4Macro(i1,i2-1,i3,m)+8.*u ## tt4Macro(i1,i2+1,i3,m)-u ## tt4Macro(i1,i2+2,i3,m))/(12.*dr(1))
 #defineMacro u ## rstt4Macro(i1,i2,i3,m) (u ## stt4Macro(i1-2,i2,i3,m)-8.*u ## stt4Macro(i1-1,i2,i3,m)+8.*u ## stt4Macro(i1+1,i2,i3,m)-u ## stt4Macro(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## rrstt4Macro(i1,i2,i3,m) (-u ## stt4Macro(i1-2,i2,i3,m)+16.*u ## stt4Macro(i1-1,i2,i3,m)-30.*u ## stt4Macro(i1,i2,i3,m)+16.*u ## stt4Macro(i1+1,i2,i3,m)-u ## stt4Macro(i1+2,i2,i3,m))/(12.*dr(0)**2)
 #defineMacro u ## sstt4Macro(i1,i2,i3,m) (-u ## tt4Macro(i1,i2-2,i3,m)+16.*u ## tt4Macro(i1,i2-1,i3,m)-30.*u ## tt4Macro(i1,i2,i3,m)+16.*u ## tt4Macro(i1,i2+1,i3,m)-u ## tt4Macro(i1,i2+2,i3,m))/(12.*dr(1)**2)
 #defineMacro u ## rsstt4Macro(i1,i2,i3,m) (u ## sstt4Macro(i1-2,i2,i3,m)-8.*u ## sstt4Macro(i1-1,i2,i3,m)+8.*u ## sstt4Macro(i1+1,i2,i3,m)-u ## sstt4Macro(i1+2,i2,i3,m))/(12.*dr(0))
 #defineMacro u ## rrsstt4Macro(i1,i2,i3,m) (-u ## sstt4Macro(i1-2,i2,i3,m)+16.*u ## sstt4Macro(i1-1,i2,i3,m)-30.*u ## sstt4Macro(i1,i2,i3,m)+16.*u ## sstt4Macro(i1+1,i2,i3,m)-u ## sstt4Macro(i1+2,i2,i3,m))/(12.*dr(0)**2)
#endMacro

! =======================================================
!  Macro to compute Third derivatives in 2 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================
#beginMacro getThirdDerivatives2d(ORDER,GRIDTYPE,OPTION,i1,i2,i3)

#If #GRIDTYPE eq "rectangular" 
! ---------- RECTANGULAR  ---------
! This assumes dr(0:2) = dx(0:2)
defineThirdParameticDerivativesComponents1(u)

#Else
! ---------- START CURVILINEAR  ---------
defineThirdParameticDerivativesComponents1(u)
#If #OPTION eq "evalMetrics"
defineThirdParameticDerivativesComponents0(rx)
defineThirdParameticDerivativesComponents0(ry)
defineThirdParameticDerivativesComponents0(sx)
defineThirdParameticDerivativesComponents0(sy)
#End

! ---------- Parametric derivatives ---------
ur       = ur4Macro(i1,i2,i3,0)
urr      = urr4Macro(i1,i2,i3,0)
urrr     = urrr2Macro(i1,i2,i3,0)
us       = us4Macro(i1,i2,i3,0)
urs      = urs4Macro(i1,i2,i3,0)
urrs     = urrs4Macro(i1,i2,i3,0)
uss      = uss4Macro(i1,i2,i3,0)
urss     = urss4Macro(i1,i2,i3,0)
usss     = usss2Macro(i1,i2,i3,0)
#If #OPTION eq "evalMetrics"
rxr      = rxr4Macro(i1,i2,i3)
rxrr     = rxrr4Macro(i1,i2,i3)
rxs      = rxs4Macro(i1,i2,i3)
rxrs     = rxrs4Macro(i1,i2,i3)
rxss     = rxss4Macro(i1,i2,i3)
ryr      = ryr4Macro(i1,i2,i3)
ryrr     = ryrr4Macro(i1,i2,i3)
rys      = rys4Macro(i1,i2,i3)
ryrs     = ryrs4Macro(i1,i2,i3)
ryss     = ryss4Macro(i1,i2,i3)
sxr      = sxr4Macro(i1,i2,i3)
sxrr     = sxrr4Macro(i1,i2,i3)
sxs      = sxs4Macro(i1,i2,i3)
sxrs     = sxrs4Macro(i1,i2,i3)
sxss     = sxss4Macro(i1,i2,i3)
syr      = syr4Macro(i1,i2,i3)
syrr     = syrr4Macro(i1,i2,i3)
sys      = sys4Macro(i1,i2,i3)
syrs     = syrs4Macro(i1,i2,i3)
syss     = syss4Macro(i1,i2,i3)

! ---------- Spatial derivatives of metrics rx, sx, ry, ... ---------
rxi = rx(i1,i2,i3)
ryi = ry(i1,i2,i3)
sxi = sx(i1,i2,i3)
syi = sy(i1,i2,i3)
rxx      = rxi*rxr+sxi*rxs
rxy      = ryi*rxr+syi*rxs
ryy      = ryi*ryr+syi*rys
sxx      = rxi*sxr+sxi*sxs
sxy      = ryi*sxr+syi*sxs
syy      = ryi*syr+syi*sys
rxxx     = rxi**2*rxrr+2.*rxi*sxi*rxrs+sxi**2*rxss+rxx*rxr+sxx*rxs
rxxy     = rxi*ryi*rxrr+(rxi*syi+ryi*sxi)*rxrs+sxi*syi*rxss+rxy*rxr+sxy*rxs
rxyy     = ryi**2*rxrr+2.*ryi*syi*rxrs+syi**2*rxss+ryy*rxr+syy*rxs
ryyy     = ryi**2*ryrr+2.*ryi*syi*ryrs+syi**2*ryss+ryy*ryr+syy*rys
sxxx     = rxi**2*sxrr+2.*rxi*sxi*sxrs+sxi**2*sxss+rxx*sxr+sxx*sxs
sxxy     = rxi*ryi*sxrr+(rxi*syi+ryi*sxi)*sxrs+sxi*syi*sxss+rxy*sxr+sxy*sxs
sxyy     = ryi**2*sxrr+2.*ryi*syi*sxrs+syi**2*sxss+ryy*sxr+syy*sxs
syyy     = ryi**2*syrr+2.*ryi*syi*syrs+syi**2*syss+ryy*syr+syy*sys
#End
! ---- end OPTION eq evalMetrics ---

! ---------- Third spatial derivatives of u ---------
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
! This assumes dr(0:2) = dx(0:2)
defineThirdParameticDerivativesComponents1(u)
uxyz     = urst4Macro(i1,i2,i3,0)

#Else
! ---------- START CURVILINEAR  ---------
defineThirdParameticDerivativesComponents1(u)
#If #OPTION eq "evalMetrics"
defineThirdParameticDerivativesComponents0(rx)
defineThirdParameticDerivativesComponents0(ry)
defineThirdParameticDerivativesComponents0(sx)
defineThirdParameticDerivativesComponents0(sy)
defineThirdParameticDerivativesComponents0(rz)
defineThirdParameticDerivativesComponents0(sz)
defineThirdParameticDerivativesComponents0(tx)
defineThirdParameticDerivativesComponents0(ty)
defineThirdParameticDerivativesComponents0(tz)
#End

! ---------- Parametric derivatives ---------
ur       = ur4Macro(i1,i2,i3,0)
urr      = urr4Macro(i1,i2,i3,0)
urrr     = urrr2Macro(i1,i2,i3,0)
us       = us4Macro(i1,i2,i3,0)
urs      = urs4Macro(i1,i2,i3,0)
urrs     = urrs4Macro(i1,i2,i3,0)
uss      = uss4Macro(i1,i2,i3,0)
urss     = urss4Macro(i1,i2,i3,0)
usss     = usss2Macro(i1,i2,i3,0)
ut       = ut4Macro(i1,i2,i3,0)
urt      = urt4Macro(i1,i2,i3,0)
urrt     = urrt4Macro(i1,i2,i3,0)
ust      = ust4Macro(i1,i2,i3,0)
urst     = urst4Macro(i1,i2,i3,0)
usst     = usst4Macro(i1,i2,i3,0)
utt      = utt4Macro(i1,i2,i3,0)
urtt     = urtt4Macro(i1,i2,i3,0)
ustt     = ustt4Macro(i1,i2,i3,0)
uttt     = uttt2Macro(i1,i2,i3,0)
#If #OPTION eq "evalMetrics"
rxr      = rxr4Macro(i1,i2,i3)
rxrr     = rxrr4Macro(i1,i2,i3)
rxs      = rxs4Macro(i1,i2,i3)
rxrs     = rxrs4Macro(i1,i2,i3)
rxss     = rxss4Macro(i1,i2,i3)
rxt      = rxt4Macro(i1,i2,i3)
rxrt     = rxrt4Macro(i1,i2,i3)
rxst     = rxst4Macro(i1,i2,i3)
rxtt     = rxtt4Macro(i1,i2,i3)
ryr      = ryr4Macro(i1,i2,i3)
ryrr     = ryrr4Macro(i1,i2,i3)
rys      = rys4Macro(i1,i2,i3)
ryrs     = ryrs4Macro(i1,i2,i3)
ryss     = ryss4Macro(i1,i2,i3)
ryt      = ryt4Macro(i1,i2,i3)
ryrt     = ryrt4Macro(i1,i2,i3)
ryst     = ryst4Macro(i1,i2,i3)
rytt     = rytt4Macro(i1,i2,i3)
sxr      = sxr4Macro(i1,i2,i3)
sxrr     = sxrr4Macro(i1,i2,i3)
sxs      = sxs4Macro(i1,i2,i3)
sxrs     = sxrs4Macro(i1,i2,i3)
sxss     = sxss4Macro(i1,i2,i3)
sxt      = sxt4Macro(i1,i2,i3)
sxrt     = sxrt4Macro(i1,i2,i3)
sxst     = sxst4Macro(i1,i2,i3)
sxtt     = sxtt4Macro(i1,i2,i3)
syr      = syr4Macro(i1,i2,i3)
syrr     = syrr4Macro(i1,i2,i3)
sys      = sys4Macro(i1,i2,i3)
syrs     = syrs4Macro(i1,i2,i3)
syss     = syss4Macro(i1,i2,i3)
syt      = syt4Macro(i1,i2,i3)
syrt     = syrt4Macro(i1,i2,i3)
syst     = syst4Macro(i1,i2,i3)
sytt     = sytt4Macro(i1,i2,i3)
rzr      = rzr4Macro(i1,i2,i3)
rzrr     = rzrr4Macro(i1,i2,i3)
rzs      = rzs4Macro(i1,i2,i3)
rzrs     = rzrs4Macro(i1,i2,i3)
rzss     = rzss4Macro(i1,i2,i3)
rzt      = rzt4Macro(i1,i2,i3)
rzrt     = rzrt4Macro(i1,i2,i3)
rzst     = rzst4Macro(i1,i2,i3)
rztt     = rztt4Macro(i1,i2,i3)
szr      = szr4Macro(i1,i2,i3)
szrr     = szrr4Macro(i1,i2,i3)
szs      = szs4Macro(i1,i2,i3)
szrs     = szrs4Macro(i1,i2,i3)
szss     = szss4Macro(i1,i2,i3)
szt      = szt4Macro(i1,i2,i3)
szrt     = szrt4Macro(i1,i2,i3)
szst     = szst4Macro(i1,i2,i3)
sztt     = sztt4Macro(i1,i2,i3)
txr      = txr4Macro(i1,i2,i3)
txrr     = txrr4Macro(i1,i2,i3)
txs      = txs4Macro(i1,i2,i3)
txrs     = txrs4Macro(i1,i2,i3)
txss     = txss4Macro(i1,i2,i3)
txt      = txt4Macro(i1,i2,i3)
txrt     = txrt4Macro(i1,i2,i3)
txst     = txst4Macro(i1,i2,i3)
txtt     = txtt4Macro(i1,i2,i3)
tyr      = tyr4Macro(i1,i2,i3)
tyrr     = tyrr4Macro(i1,i2,i3)
tys      = tys4Macro(i1,i2,i3)
tyrs     = tyrs4Macro(i1,i2,i3)
tyss     = tyss4Macro(i1,i2,i3)
tyt      = tyt4Macro(i1,i2,i3)
tyrt     = tyrt4Macro(i1,i2,i3)
tyst     = tyst4Macro(i1,i2,i3)
tytt     = tytt4Macro(i1,i2,i3)
tzr      = tzr4Macro(i1,i2,i3)
tzrr     = tzrr4Macro(i1,i2,i3)
tzs      = tzs4Macro(i1,i2,i3)
tzrs     = tzrs4Macro(i1,i2,i3)
tzss     = tzss4Macro(i1,i2,i3)
tzt      = tzt4Macro(i1,i2,i3)
tzrt     = tzrt4Macro(i1,i2,i3)
tzst     = tzst4Macro(i1,i2,i3)
tztt     = tztt4Macro(i1,i2,i3)

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
rxx      = rxi*rxr+sxi*rxs+txi*rxt
rxy      = ryi*rxr+syi*rxs+tyi*rxt
rxz      = rzi*rxr+szi*rxs+tzi*rxt
ryy      = ryi*ryr+syi*rys+tyi*ryt
ryz      = rzi*ryr+szi*rys+tzi*ryt
rzz      = rzi*rzr+szi*rzs+tzi*rzt
sxx      = rxi*sxr+sxi*sxs+txi*sxt
sxy      = ryi*sxr+syi*sxs+tyi*sxt
sxz      = rzi*sxr+szi*sxs+tzi*sxt
syy      = ryi*syr+syi*sys+tyi*syt
syz      = rzi*syr+szi*sys+tzi*syt
szz      = rzi*szr+szi*szs+tzi*szt
txx      = rxi*txr+sxi*txs+txi*txt
txy      = ryi*txr+syi*txs+tyi*txt
txz      = rzi*txr+szi*txs+tzi*txt
tyy      = ryi*tyr+syi*tys+tyi*tyt
tyz      = rzi*tyr+szi*tys+tzi*tyt
tzz      = rzi*tzr+szi*tzs+tzi*tzt
rxxx     = rxi**2*rxrr+2.*rxi*sxi*rxrs+2.*rxi*txi*rxrt+sxi**2*rxss+2.*sxi*txi*rxst+txi**2*rxtt+rxx*rxr+sxx*rxs+txx*rxt
rxxy     = rxi*ryi*rxrr+(rxi*syi+ryi*sxi)*rxrs+sxi*syi*rxss+(rxi*tyi+ryi*txi)*rxrt+(sxi*tyi+syi*txi)*rxst+txi*tyi*rxtt+rxy*rxr+sxy*rxs+txy*rxt
rxyy     = ryi**2*rxrr+2.*ryi*syi*rxrs+2.*ryi*tyi*rxrt+syi**2*rxss+2.*syi*tyi*rxst+tyi**2*rxtt+ryy*rxr+syy*rxs+tyy*rxt
rxxz     = rxi*rzi*rxrr+(rxi*szi+rzi*sxi)*rxrs+sxi*szi*rxss+(rxi*tzi+rzi*txi)*rxrt+(sxi*tzi+szi*txi)*rxst+txi*tzi*rxtt+rxz*rxr+sxz*rxs+txz*rxt
rxyz     = ryi*rzi*rxrr+(ryi*szi+rzi*syi)*rxrs+syi*szi*rxss+(ryi*tzi+rzi*tyi)*rxrt+(syi*tzi+szi*tyi)*rxst+tyi*tzi*rxtt+ryz*rxr+syz*rxs+tyz*rxt
rxzz     = rzi**2*rxrr+2.*rzi*szi*rxrs+2.*rzi*tzi*rxrt+szi**2*rxss+2.*szi*tzi*rxst+tzi**2*rxtt+rzz*rxr+szz*rxs+tzz*rxt
ryyy     = ryi**2*ryrr+2.*ryi*syi*ryrs+2.*ryi*tyi*ryrt+syi**2*ryss+2.*syi*tyi*ryst+tyi**2*rytt+ryy*ryr+syy*rys+tyy*ryt
ryyz     = ryi*rzi*ryrr+(ryi*szi+rzi*syi)*ryrs+syi*szi*ryss+(ryi*tzi+rzi*tyi)*ryrt+(syi*tzi+szi*tyi)*ryst+tyi*tzi*rytt+ryz*ryr+syz*rys+tyz*ryt
ryzz     = rzi**2*ryrr+2.*rzi*szi*ryrs+2.*rzi*tzi*ryrt+szi**2*ryss+2.*szi*tzi*ryst+tzi**2*rytt+rzz*ryr+szz*rys+tzz*ryt
rzzz     = rzi**2*rzrr+2.*rzi*szi*rzrs+2.*rzi*tzi*rzrt+szi**2*rzss+2.*szi*tzi*rzst+tzi**2*rztt+rzz*rzr+szz*rzs+tzz*rzt
sxxx     = rxi**2*sxrr+2.*rxi*sxi*sxrs+2.*rxi*txi*sxrt+sxi**2*sxss+2.*sxi*txi*sxst+txi**2*sxtt+rxx*sxr+sxx*sxs+txx*sxt
sxxy     = rxi*ryi*sxrr+(rxi*syi+ryi*sxi)*sxrs+sxi*syi*sxss+(rxi*tyi+ryi*txi)*sxrt+(sxi*tyi+syi*txi)*sxst+txi*tyi*sxtt+rxy*sxr+sxy*sxs+txy*sxt
sxyy     = ryi**2*sxrr+2.*ryi*syi*sxrs+2.*ryi*tyi*sxrt+syi**2*sxss+2.*syi*tyi*sxst+tyi**2*sxtt+ryy*sxr+syy*sxs+tyy*sxt
sxxz     = rxi*rzi*sxrr+(rxi*szi+rzi*sxi)*sxrs+sxi*szi*sxss+(rxi*tzi+rzi*txi)*sxrt+(sxi*tzi+szi*txi)*sxst+txi*tzi*sxtt+rxz*sxr+sxz*sxs+txz*sxt
sxyz     = ryi*rzi*sxrr+(ryi*szi+rzi*syi)*sxrs+syi*szi*sxss+(ryi*tzi+rzi*tyi)*sxrt+(syi*tzi+szi*tyi)*sxst+tyi*tzi*sxtt+ryz*sxr+syz*sxs+tyz*sxt
sxzz     = rzi**2*sxrr+2.*rzi*szi*sxrs+2.*rzi*tzi*sxrt+szi**2*sxss+2.*szi*tzi*sxst+tzi**2*sxtt+rzz*sxr+szz*sxs+tzz*sxt
syyy     = ryi**2*syrr+2.*ryi*syi*syrs+2.*ryi*tyi*syrt+syi**2*syss+2.*syi*tyi*syst+tyi**2*sytt+ryy*syr+syy*sys+tyy*syt
syyz     = ryi*rzi*syrr+(ryi*szi+rzi*syi)*syrs+syi*szi*syss+(ryi*tzi+rzi*tyi)*syrt+(syi*tzi+szi*tyi)*syst+tyi*tzi*sytt+ryz*syr+syz*sys+tyz*syt
syzz     = rzi**2*syrr+2.*rzi*szi*syrs+2.*rzi*tzi*syrt+szi**2*syss+2.*szi*tzi*syst+tzi**2*sytt+rzz*syr+szz*sys+tzz*syt
szzz     = rzi**2*szrr+2.*rzi*szi*szrs+2.*rzi*tzi*szrt+szi**2*szss+2.*szi*tzi*szst+tzi**2*sztt+rzz*szr+szz*szs+tzz*szt
txxx     = rxi**2*txrr+2.*rxi*sxi*txrs+2.*rxi*txi*txrt+sxi**2*txss+2.*sxi*txi*txst+txi**2*txtt+rxx*txr+sxx*txs+txx*txt
txxy     = rxi*ryi*txrr+(rxi*syi+ryi*sxi)*txrs+sxi*syi*txss+(rxi*tyi+ryi*txi)*txrt+(sxi*tyi+syi*txi)*txst+txi*tyi*txtt+rxy*txr+sxy*txs+txy*txt
txyy     = ryi**2*txrr+2.*ryi*syi*txrs+2.*ryi*tyi*txrt+syi**2*txss+2.*syi*tyi*txst+tyi**2*txtt+ryy*txr+syy*txs+tyy*txt
txxz     = rxi*rzi*txrr+(rxi*szi+rzi*sxi)*txrs+sxi*szi*txss+(rxi*tzi+rzi*txi)*txrt+(sxi*tzi+szi*txi)*txst+txi*tzi*txtt+rxz*txr+sxz*txs+txz*txt
txyz     = ryi*rzi*txrr+(ryi*szi+rzi*syi)*txrs+syi*szi*txss+(ryi*tzi+rzi*tyi)*txrt+(syi*tzi+szi*tyi)*txst+tyi*tzi*txtt+ryz*txr+syz*txs+tyz*txt
txzz     = rzi**2*txrr+2.*rzi*szi*txrs+2.*rzi*tzi*txrt+szi**2*txss+2.*szi*tzi*txst+tzi**2*txtt+rzz*txr+szz*txs+tzz*txt
tyyy     = ryi**2*tyrr+2.*ryi*syi*tyrs+2.*ryi*tyi*tyrt+syi**2*tyss+2.*syi*tyi*tyst+tyi**2*tytt+ryy*tyr+syy*tys+tyy*tyt
tyyz     = ryi*rzi*tyrr+(ryi*szi+rzi*syi)*tyrs+syi*szi*tyss+(ryi*tzi+rzi*tyi)*tyrt+(syi*tzi+szi*tyi)*tyst+tyi*tzi*tytt+ryz*tyr+syz*tys+tyz*tyt
tyzz     = rzi**2*tyrr+2.*rzi*szi*tyrs+2.*rzi*tzi*tyrt+szi**2*tyss+2.*szi*tzi*tyst+tzi**2*tytt+rzz*tyr+szz*tys+tzz*tyt
tzzz     = rzi**2*tzrr+2.*rzi*szi*tzrs+2.*rzi*tzi*tzrt+szi**2*tzss+2.*szi*tzi*tzst+tzi**2*tztt+rzz*tzr+szz*tzs+tzz*tzt
#End
! ---- end OPTION eq evalMetrics ---

! ---------- Third spatial derivatives of u ---------
uxyz     = rxi*ryi*rzi*urrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*urrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*urss+sxi*syi*szi*usss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*urrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*urst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*usst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*urtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*ustt+txi*tyi*tzi*uttt+(rxi*ryz+rxy*rzi+rxz*ryi)*urr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*urs+(sxi*syz+sxy*szi+sxz*syi)*uss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*urt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*ust+(txi*tyz+txy*tzi+txz*tyi)*utt+rxyz*ur+sxyz*us+txyz*ut
! ---------- END CURVILINEAR  ---------
#End
#endMacro
