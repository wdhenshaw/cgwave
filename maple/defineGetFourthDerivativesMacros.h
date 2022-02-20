! ****** File written by makeGetDerivativesMacros.maple  ******

#beginMacro defineFourthParameticDerivativesComponents0(u)
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

#beginMacro defineFourthParameticDerivativesComponents1(u)
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
!  Macro to compute Fourth derivatives in 2 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================
#beginMacro getFourthDerivatives2d(ORDER,GRIDTYPE,OPTION,i1,i2,i3)

#If #GRIDTYPE eq "rectangular" 
! ---------- RECTANGULAR  ---------
! This assumes dr(0:2) = dx(0:2)
defineFourthParameticDerivativesComponents1(u)
uxx      = urr4Macro(i1,i2,i3,0)
uxxxx    = urrrr2Macro(i1,i2,i3,0)
uyy      = uss4Macro(i1,i2,i3,0)
uxxyy    = urrss4Macro(i1,i2,i3,0)
uyyyy    = ussss2Macro(i1,i2,i3,0)

#Else
! ---------- START CURVILINEAR  ---------
defineFourthParameticDerivativesComponents1(u)
#If #OPTION eq "evalMetrics"
defineFourthParameticDerivativesComponents0(rx)
defineFourthParameticDerivativesComponents0(ry)
defineFourthParameticDerivativesComponents0(sx)
defineFourthParameticDerivativesComponents0(sy)
#End

! ---------- Parametric derivatives ---------
ur       = ur4Macro(i1,i2,i3,0)
urr      = urr4Macro(i1,i2,i3,0)
urrr     = urrr2Macro(i1,i2,i3,0)
urrrr    = urrrr2Macro(i1,i2,i3,0)
us       = us4Macro(i1,i2,i3,0)
urs      = urs4Macro(i1,i2,i3,0)
urrs     = urrs4Macro(i1,i2,i3,0)
urrrs    = urrrs2Macro(i1,i2,i3,0)
uss      = uss4Macro(i1,i2,i3,0)
urss     = urss4Macro(i1,i2,i3,0)
urrss    = urrss4Macro(i1,i2,i3,0)
usss     = usss2Macro(i1,i2,i3,0)
ursss    = ursss2Macro(i1,i2,i3,0)
ussss    = ussss2Macro(i1,i2,i3,0)
#If #OPTION eq "evalMetrics"
rxr      = rxr4Macro(i1,i2,i3)
rxrr     = rxrr4Macro(i1,i2,i3)
rxrrr    = rxrrr2Macro(i1,i2,i3)
rxs      = rxs4Macro(i1,i2,i3)
rxrs     = rxrs4Macro(i1,i2,i3)
rxrrs    = rxrrs4Macro(i1,i2,i3)
rxss     = rxss4Macro(i1,i2,i3)
rxrss    = rxrss4Macro(i1,i2,i3)
rxsss    = rxsss2Macro(i1,i2,i3)
ryr      = ryr4Macro(i1,i2,i3)
ryrr     = ryrr4Macro(i1,i2,i3)
ryrrr    = ryrrr2Macro(i1,i2,i3)
rys      = rys4Macro(i1,i2,i3)
ryrs     = ryrs4Macro(i1,i2,i3)
ryrrs    = ryrrs4Macro(i1,i2,i3)
ryss     = ryss4Macro(i1,i2,i3)
ryrss    = ryrss4Macro(i1,i2,i3)
rysss    = rysss2Macro(i1,i2,i3)
sxr      = sxr4Macro(i1,i2,i3)
sxrr     = sxrr4Macro(i1,i2,i3)
sxrrr    = sxrrr2Macro(i1,i2,i3)
sxs      = sxs4Macro(i1,i2,i3)
sxrs     = sxrs4Macro(i1,i2,i3)
sxrrs    = sxrrs4Macro(i1,i2,i3)
sxss     = sxss4Macro(i1,i2,i3)
sxrss    = sxrss4Macro(i1,i2,i3)
sxsss    = sxsss2Macro(i1,i2,i3)
syr      = syr4Macro(i1,i2,i3)
syrr     = syrr4Macro(i1,i2,i3)
syrrr    = syrrr2Macro(i1,i2,i3)
sys      = sys4Macro(i1,i2,i3)
syrs     = syrs4Macro(i1,i2,i3)
syrrs    = syrrs4Macro(i1,i2,i3)
syss     = syss4Macro(i1,i2,i3)
syrss    = syrss4Macro(i1,i2,i3)
sysss    = sysss2Macro(i1,i2,i3)

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
rxxxx    = rxi**3*rxrrr+3.*rxi**2*sxi*rxrrs+3.*rxi*sxi**2*rxrss+sxi**3*rxsss+3.*rxi*rxx*rxrr+(3.*rxi*sxx+3.*rxx*sxi)*rxrs+3.*sxi*sxx*rxss+rxxx*rxr+sxxx*rxs
rxxxy    = rxi**2*ryi*rxrrr+(rxi**2*syi+2.*rxi*ryi*sxi)*rxrrs+(2.*rxi*sxi*syi+ryi*sxi**2)*rxrss+sxi**2*syi*rxsss+(2.*rxi*rxy+rxx*ryi)*rxrr+(2.*rxi*sxy+rxx*syi+2.*rxy*sxi+ryi*sxx)*rxrs+(2.*sxi*sxy+sxx*syi)*rxss+rxxy*rxr+sxxy*rxs
rxxyy    = rxi*ryi**2*rxrrr+(2.*rxi*ryi*syi+ryi**2*sxi)*rxrrs+(rxi*syi**2+2.*ryi*sxi*syi)*rxrss+sxi*syi**2*rxsss+(rxi*ryy+2.*rxy*ryi)*rxrr+(rxi*syy+2.*rxy*syi+2.*ryi*sxy+ryy*sxi)*rxrs+(sxi*syy+2.*sxy*syi)*rxss+rxyy*rxr+sxyy*rxs
rxyyy    = ryi**3*rxrrr+3.*ryi**2*syi*rxrrs+3.*ryi*syi**2*rxrss+syi**3*rxsss+3.*ryi*ryy*rxrr+(3.*ryi*syy+3.*ryy*syi)*rxrs+3.*syi*syy*rxss+ryyy*rxr+syyy*rxs
ryyyy    = ryi**3*ryrrr+3.*ryi**2*syi*ryrrs+3.*ryi*syi**2*ryrss+syi**3*rysss+3.*ryi*ryy*ryrr+(3.*ryi*syy+3.*ryy*syi)*ryrs+3.*syi*syy*ryss+ryyy*ryr+syyy*rys
sxxxx    = rxi**3*sxrrr+3.*rxi**2*sxi*sxrrs+3.*rxi*sxi**2*sxrss+sxi**3*sxsss+3.*rxi*rxx*sxrr+(3.*rxi*sxx+3.*rxx*sxi)*sxrs+3.*sxi*sxx*sxss+rxxx*sxr+sxxx*sxs
sxxxy    = rxi**2*ryi*sxrrr+(rxi**2*syi+2.*rxi*ryi*sxi)*sxrrs+(2.*rxi*sxi*syi+ryi*sxi**2)*sxrss+sxi**2*syi*sxsss+(2.*rxi*rxy+rxx*ryi)*sxrr+(2.*rxi*sxy+rxx*syi+2.*rxy*sxi+ryi*sxx)*sxrs+(2.*sxi*sxy+sxx*syi)*sxss+rxxy*sxr+sxxy*sxs
sxxyy    = rxi*ryi**2*sxrrr+(2.*rxi*ryi*syi+ryi**2*sxi)*sxrrs+(rxi*syi**2+2.*ryi*sxi*syi)*sxrss+sxi*syi**2*sxsss+(rxi*ryy+2.*rxy*ryi)*sxrr+(rxi*syy+2.*rxy*syi+2.*ryi*sxy+ryy*sxi)*sxrs+(sxi*syy+2.*sxy*syi)*sxss+rxyy*sxr+sxyy*sxs
sxyyy    = ryi**3*sxrrr+3.*ryi**2*syi*sxrrs+3.*ryi*syi**2*sxrss+syi**3*sxsss+3.*ryi*ryy*sxrr+(3.*ryi*syy+3.*ryy*syi)*sxrs+3.*syi*syy*sxss+ryyy*sxr+syyy*sxs
syyyy    = ryi**3*syrrr+3.*ryi**2*syi*syrrs+3.*ryi*syi**2*syrss+syi**3*sysss+3.*ryi*ryy*syrr+(3.*ryi*syy+3.*ryy*syi)*syrs+3.*syi*syy*syss+ryyy*syr+syyy*sys
#End
! ---- end OPTION eq evalMetrics ---

! ---------- Fourth spatial derivatives of u ---------
uxx      = rxi**2*urr+2.*rxi*sxi*urs+sxi**2*uss+rxx*ur+sxx*us
uxxxx    = rxi**4*urrrr+4.*rxi**3*sxi*urrrs+6.*rxi**2*sxi**2*urrss+4.*rxi*sxi**3*ursss+sxi**4*ussss+6.*rxi**2*rxx*urrr+(6.*rxi**2*sxx+12.*rxi*rxx*sxi)*urrs+(12.*rxi*sxi*sxx+6.*rxx*sxi**2)*urss+6.*sxi**2*sxx*usss+(4.*rxi*rxxx+3.*rxx**2)*urr+(4.*rxi*sxxx+6.*rxx*sxx+4.*rxxx*sxi)*urs+(4.*sxi*sxxx+3.*sxx**2)*uss+rxxxx*ur+sxxxx*us
uyy      = ryi**2*urr+2.*ryi*syi*urs+syi**2*uss+ryy*ur+syy*us
uxxyy    = rxi**2*ryi**2*urrrr+(2.*rxi**2*ryi*syi+2.*rxi*ryi**2*sxi)*urrrs+(rxi**2*syi**2+4.*rxi*ryi*sxi*syi+ryi**2*sxi**2)*urrss+(2.*rxi*sxi*syi**2+2.*ryi*sxi**2*syi)*ursss+sxi**2*syi**2*ussss+(rxi**2*ryy+4.*rxi*rxy*ryi+rxx*ryi**2)*urrr+(rxi**2*syy+(4.*rxy*syi+4.*ryi*sxy+2.*ryy*sxi)*rxi+ryi**2*sxx+(2.*rxx*syi+4.*rxy*sxi)*ryi)*urrs+((2.*sxi*syy+4.*sxy*syi)*rxi+(4.*sxi*sxy+2.*sxx*syi)*ryi+ryy*sxi**2+4.*rxy*syi*sxi+rxx*syi**2)*urss+(sxi**2*syy+4.*sxi*sxy*syi+sxx*syi**2)*usss+(2.*rxi*rxyy+rxx*ryy+2.*rxxy*ryi+2.*rxy**2)*urr+(2.*rxi*sxyy+rxx*syy+2.*rxxy*syi+4.*rxy*sxy+2.*rxyy*sxi+2.*ryi*sxxy+ryy*sxx)*urs+(2.*sxi*sxyy+sxx*syy+2.*sxxy*syi+2.*sxy**2)*uss+rxxyy*ur+sxxyy*us
uyyyy    = ryi**4*urrrr+4.*ryi**3*syi*urrrs+6.*ryi**2*syi**2*urrss+4.*ryi*syi**3*ursss+syi**4*ussss+6.*ryi**2*ryy*urrr+(6.*ryi**2*syy+12.*ryi*ryy*syi)*urrs+(12.*ryi*syi*syy+6.*ryy*syi**2)*urss+6.*syi**2*syy*usss+(4.*ryi*ryyy+3.*ryy**2)*urr+(4.*ryi*syyy+6.*ryy*syy+4.*ryyy*syi)*urs+(4.*syi*syyy+3.*syy**2)*uss+ryyyy*ur+syyyy*us
! ---------- END CURVILINEAR  ---------
#End
#endMacro

! =======================================================
!  Macro to compute Fourth derivatives in 3 dimensions 
!  OPTION : evalMetrics : evaluate the derivatives of the metrics
!          (metrics need only be evaluated once when using discrete delta to get coeffs)
! =======================================================
#beginMacro getFourthDerivatives3d(ORDER,GRIDTYPE,OPTION,i1,i2,i3)

#If #GRIDTYPE eq "rectangular" 
! ---------- RECTANGULAR  ---------
! This assumes dr(0:2) = dx(0:2)
defineFourthParameticDerivativesComponents1(u)
uxx      = urr4Macro(i1,i2,i3,0)
uxxxx    = urrrr2Macro(i1,i2,i3,0)
uyy      = uss4Macro(i1,i2,i3,0)
uxxyy    = urrss4Macro(i1,i2,i3,0)
uyyyy    = ussss2Macro(i1,i2,i3,0)
uzz      = utt4Macro(i1,i2,i3,0)
uxxzz    = urrtt4Macro(i1,i2,i3,0)
uyyzz    = usstt4Macro(i1,i2,i3,0)
uzzzz    = utttt2Macro(i1,i2,i3,0)

#Else
! ---------- START CURVILINEAR  ---------
defineFourthParameticDerivativesComponents1(u)
#If #OPTION eq "evalMetrics"
defineFourthParameticDerivativesComponents0(rx)
defineFourthParameticDerivativesComponents0(ry)
defineFourthParameticDerivativesComponents0(sx)
defineFourthParameticDerivativesComponents0(sy)
defineFourthParameticDerivativesComponents0(rz)
defineFourthParameticDerivativesComponents0(sz)
defineFourthParameticDerivativesComponents0(tx)
defineFourthParameticDerivativesComponents0(ty)
defineFourthParameticDerivativesComponents0(tz)
#End

! ---------- Parametric derivatives ---------
ur       = ur4Macro(i1,i2,i3,0)
urr      = urr4Macro(i1,i2,i3,0)
urrr     = urrr2Macro(i1,i2,i3,0)
urrrr    = urrrr2Macro(i1,i2,i3,0)
us       = us4Macro(i1,i2,i3,0)
urs      = urs4Macro(i1,i2,i3,0)
urrs     = urrs4Macro(i1,i2,i3,0)
urrrs    = urrrs2Macro(i1,i2,i3,0)
uss      = uss4Macro(i1,i2,i3,0)
urss     = urss4Macro(i1,i2,i3,0)
urrss    = urrss4Macro(i1,i2,i3,0)
usss     = usss2Macro(i1,i2,i3,0)
ursss    = ursss2Macro(i1,i2,i3,0)
ussss    = ussss2Macro(i1,i2,i3,0)
ut       = ut4Macro(i1,i2,i3,0)
urt      = urt4Macro(i1,i2,i3,0)
urrt     = urrt4Macro(i1,i2,i3,0)
urrrt    = urrrt2Macro(i1,i2,i3,0)
ust      = ust4Macro(i1,i2,i3,0)
urst     = urst4Macro(i1,i2,i3,0)
urrst    = urrst4Macro(i1,i2,i3,0)
usst     = usst4Macro(i1,i2,i3,0)
ursst    = ursst4Macro(i1,i2,i3,0)
ussst    = ussst2Macro(i1,i2,i3,0)
utt      = utt4Macro(i1,i2,i3,0)
urtt     = urtt4Macro(i1,i2,i3,0)
urrtt    = urrtt4Macro(i1,i2,i3,0)
ustt     = ustt4Macro(i1,i2,i3,0)
urstt    = urstt4Macro(i1,i2,i3,0)
usstt    = usstt4Macro(i1,i2,i3,0)
uttt     = uttt2Macro(i1,i2,i3,0)
urttt    = urttt2Macro(i1,i2,i3,0)
usttt    = usttt2Macro(i1,i2,i3,0)
utttt    = utttt2Macro(i1,i2,i3,0)
#If #OPTION eq "evalMetrics"
rxr      = rxr4Macro(i1,i2,i3)
rxrr     = rxrr4Macro(i1,i2,i3)
rxrrr    = rxrrr2Macro(i1,i2,i3)
rxs      = rxs4Macro(i1,i2,i3)
rxrs     = rxrs4Macro(i1,i2,i3)
rxrrs    = rxrrs4Macro(i1,i2,i3)
rxss     = rxss4Macro(i1,i2,i3)
rxrss    = rxrss4Macro(i1,i2,i3)
rxsss    = rxsss2Macro(i1,i2,i3)
rxt      = rxt4Macro(i1,i2,i3)
rxrt     = rxrt4Macro(i1,i2,i3)
rxrrt    = rxrrt4Macro(i1,i2,i3)
rxst     = rxst4Macro(i1,i2,i3)
rxrst    = rxrst4Macro(i1,i2,i3)
rxsst    = rxsst4Macro(i1,i2,i3)
rxtt     = rxtt4Macro(i1,i2,i3)
rxrtt    = rxrtt4Macro(i1,i2,i3)
rxstt    = rxstt4Macro(i1,i2,i3)
rxttt    = rxttt2Macro(i1,i2,i3)
ryr      = ryr4Macro(i1,i2,i3)
ryrr     = ryrr4Macro(i1,i2,i3)
ryrrr    = ryrrr2Macro(i1,i2,i3)
rys      = rys4Macro(i1,i2,i3)
ryrs     = ryrs4Macro(i1,i2,i3)
ryrrs    = ryrrs4Macro(i1,i2,i3)
ryss     = ryss4Macro(i1,i2,i3)
ryrss    = ryrss4Macro(i1,i2,i3)
rysss    = rysss2Macro(i1,i2,i3)
ryt      = ryt4Macro(i1,i2,i3)
ryrt     = ryrt4Macro(i1,i2,i3)
ryrrt    = ryrrt4Macro(i1,i2,i3)
ryst     = ryst4Macro(i1,i2,i3)
ryrst    = ryrst4Macro(i1,i2,i3)
rysst    = rysst4Macro(i1,i2,i3)
rytt     = rytt4Macro(i1,i2,i3)
ryrtt    = ryrtt4Macro(i1,i2,i3)
rystt    = rystt4Macro(i1,i2,i3)
ryttt    = ryttt2Macro(i1,i2,i3)
sxr      = sxr4Macro(i1,i2,i3)
sxrr     = sxrr4Macro(i1,i2,i3)
sxrrr    = sxrrr2Macro(i1,i2,i3)
sxs      = sxs4Macro(i1,i2,i3)
sxrs     = sxrs4Macro(i1,i2,i3)
sxrrs    = sxrrs4Macro(i1,i2,i3)
sxss     = sxss4Macro(i1,i2,i3)
sxrss    = sxrss4Macro(i1,i2,i3)
sxsss    = sxsss2Macro(i1,i2,i3)
sxt      = sxt4Macro(i1,i2,i3)
sxrt     = sxrt4Macro(i1,i2,i3)
sxrrt    = sxrrt4Macro(i1,i2,i3)
sxst     = sxst4Macro(i1,i2,i3)
sxrst    = sxrst4Macro(i1,i2,i3)
sxsst    = sxsst4Macro(i1,i2,i3)
sxtt     = sxtt4Macro(i1,i2,i3)
sxrtt    = sxrtt4Macro(i1,i2,i3)
sxstt    = sxstt4Macro(i1,i2,i3)
sxttt    = sxttt2Macro(i1,i2,i3)
syr      = syr4Macro(i1,i2,i3)
syrr     = syrr4Macro(i1,i2,i3)
syrrr    = syrrr2Macro(i1,i2,i3)
sys      = sys4Macro(i1,i2,i3)
syrs     = syrs4Macro(i1,i2,i3)
syrrs    = syrrs4Macro(i1,i2,i3)
syss     = syss4Macro(i1,i2,i3)
syrss    = syrss4Macro(i1,i2,i3)
sysss    = sysss2Macro(i1,i2,i3)
syt      = syt4Macro(i1,i2,i3)
syrt     = syrt4Macro(i1,i2,i3)
syrrt    = syrrt4Macro(i1,i2,i3)
syst     = syst4Macro(i1,i2,i3)
syrst    = syrst4Macro(i1,i2,i3)
sysst    = sysst4Macro(i1,i2,i3)
sytt     = sytt4Macro(i1,i2,i3)
syrtt    = syrtt4Macro(i1,i2,i3)
systt    = systt4Macro(i1,i2,i3)
syttt    = syttt2Macro(i1,i2,i3)
rzr      = rzr4Macro(i1,i2,i3)
rzrr     = rzrr4Macro(i1,i2,i3)
rzrrr    = rzrrr2Macro(i1,i2,i3)
rzs      = rzs4Macro(i1,i2,i3)
rzrs     = rzrs4Macro(i1,i2,i3)
rzrrs    = rzrrs4Macro(i1,i2,i3)
rzss     = rzss4Macro(i1,i2,i3)
rzrss    = rzrss4Macro(i1,i2,i3)
rzsss    = rzsss2Macro(i1,i2,i3)
rzt      = rzt4Macro(i1,i2,i3)
rzrt     = rzrt4Macro(i1,i2,i3)
rzrrt    = rzrrt4Macro(i1,i2,i3)
rzst     = rzst4Macro(i1,i2,i3)
rzrst    = rzrst4Macro(i1,i2,i3)
rzsst    = rzsst4Macro(i1,i2,i3)
rztt     = rztt4Macro(i1,i2,i3)
rzrtt    = rzrtt4Macro(i1,i2,i3)
rzstt    = rzstt4Macro(i1,i2,i3)
rzttt    = rzttt2Macro(i1,i2,i3)
szr      = szr4Macro(i1,i2,i3)
szrr     = szrr4Macro(i1,i2,i3)
szrrr    = szrrr2Macro(i1,i2,i3)
szs      = szs4Macro(i1,i2,i3)
szrs     = szrs4Macro(i1,i2,i3)
szrrs    = szrrs4Macro(i1,i2,i3)
szss     = szss4Macro(i1,i2,i3)
szrss    = szrss4Macro(i1,i2,i3)
szsss    = szsss2Macro(i1,i2,i3)
szt      = szt4Macro(i1,i2,i3)
szrt     = szrt4Macro(i1,i2,i3)
szrrt    = szrrt4Macro(i1,i2,i3)
szst     = szst4Macro(i1,i2,i3)
szrst    = szrst4Macro(i1,i2,i3)
szsst    = szsst4Macro(i1,i2,i3)
sztt     = sztt4Macro(i1,i2,i3)
szrtt    = szrtt4Macro(i1,i2,i3)
szstt    = szstt4Macro(i1,i2,i3)
szttt    = szttt2Macro(i1,i2,i3)
txr      = txr4Macro(i1,i2,i3)
txrr     = txrr4Macro(i1,i2,i3)
txrrr    = txrrr2Macro(i1,i2,i3)
txs      = txs4Macro(i1,i2,i3)
txrs     = txrs4Macro(i1,i2,i3)
txrrs    = txrrs4Macro(i1,i2,i3)
txss     = txss4Macro(i1,i2,i3)
txrss    = txrss4Macro(i1,i2,i3)
txsss    = txsss2Macro(i1,i2,i3)
txt      = txt4Macro(i1,i2,i3)
txrt     = txrt4Macro(i1,i2,i3)
txrrt    = txrrt4Macro(i1,i2,i3)
txst     = txst4Macro(i1,i2,i3)
txrst    = txrst4Macro(i1,i2,i3)
txsst    = txsst4Macro(i1,i2,i3)
txtt     = txtt4Macro(i1,i2,i3)
txrtt    = txrtt4Macro(i1,i2,i3)
txstt    = txstt4Macro(i1,i2,i3)
txttt    = txttt2Macro(i1,i2,i3)
tyr      = tyr4Macro(i1,i2,i3)
tyrr     = tyrr4Macro(i1,i2,i3)
tyrrr    = tyrrr2Macro(i1,i2,i3)
tys      = tys4Macro(i1,i2,i3)
tyrs     = tyrs4Macro(i1,i2,i3)
tyrrs    = tyrrs4Macro(i1,i2,i3)
tyss     = tyss4Macro(i1,i2,i3)
tyrss    = tyrss4Macro(i1,i2,i3)
tysss    = tysss2Macro(i1,i2,i3)
tyt      = tyt4Macro(i1,i2,i3)
tyrt     = tyrt4Macro(i1,i2,i3)
tyrrt    = tyrrt4Macro(i1,i2,i3)
tyst     = tyst4Macro(i1,i2,i3)
tyrst    = tyrst4Macro(i1,i2,i3)
tysst    = tysst4Macro(i1,i2,i3)
tytt     = tytt4Macro(i1,i2,i3)
tyrtt    = tyrtt4Macro(i1,i2,i3)
tystt    = tystt4Macro(i1,i2,i3)
tyttt    = tyttt2Macro(i1,i2,i3)
tzr      = tzr4Macro(i1,i2,i3)
tzrr     = tzrr4Macro(i1,i2,i3)
tzrrr    = tzrrr2Macro(i1,i2,i3)
tzs      = tzs4Macro(i1,i2,i3)
tzrs     = tzrs4Macro(i1,i2,i3)
tzrrs    = tzrrs4Macro(i1,i2,i3)
tzss     = tzss4Macro(i1,i2,i3)
tzrss    = tzrss4Macro(i1,i2,i3)
tzsss    = tzsss2Macro(i1,i2,i3)
tzt      = tzt4Macro(i1,i2,i3)
tzrt     = tzrt4Macro(i1,i2,i3)
tzrrt    = tzrrt4Macro(i1,i2,i3)
tzst     = tzst4Macro(i1,i2,i3)
tzrst    = tzrst4Macro(i1,i2,i3)
tzsst    = tzsst4Macro(i1,i2,i3)
tztt     = tztt4Macro(i1,i2,i3)
tzrtt    = tzrtt4Macro(i1,i2,i3)
tzstt    = tzstt4Macro(i1,i2,i3)
tzttt    = tzttt2Macro(i1,i2,i3)

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
rxxxx    = rxi**3*rxrrr+3.*rxi**2*sxi*rxrrs+3.*rxi*sxi**2*rxrss+sxi**3*rxsss+3.*rxi**2*txi*rxrrt+6.*rxi*sxi*txi*rxrst+3.*sxi**2*txi*rxsst+3.*rxi*txi**2*rxrtt+3.*sxi*txi**2*rxstt+txi**3*rxttt+3.*rxi*rxx*rxrr+(3.*rxi*sxx+3.*rxx*sxi)*rxrs+3.*sxi*sxx*rxss+(3.*rxi*txx+3.*rxx*txi)*rxrt+(3.*sxi*txx+3.*sxx*txi)*rxst+3.*txi*txx*rxtt+rxxx*rxr+sxxx*rxs+txxx*rxt
rxxxy    = rxi**2*ryi*rxrrr+(rxi**2*syi+2.*rxi*ryi*sxi)*rxrrs+(2.*rxi*sxi*syi+ryi*sxi**2)*rxrss+sxi**2*syi*rxsss+(rxi**2*tyi+2.*rxi*ryi*txi)*rxrrt+((2.*sxi*tyi+2.*syi*txi)*rxi+2.*ryi*txi*sxi)*rxrst+(sxi**2*tyi+2.*sxi*syi*txi)*rxsst+(2.*rxi*txi*tyi+ryi*txi**2)*rxrtt+(2.*sxi*txi*tyi+syi*txi**2)*rxstt+txi**2*tyi*rxttt+(2.*rxi*rxy+rxx*ryi)*rxrr+(2.*rxi*sxy+rxx*syi+2.*rxy*sxi+ryi*sxx)*rxrs+(2.*sxi*sxy+sxx*syi)*rxss+(2.*rxi*txy+rxx*tyi+2.*rxy*txi+ryi*txx)*rxrt+(2.*sxi*txy+sxx*tyi+2.*sxy*txi+syi*txx)*rxst+(2.*txi*txy+txx*tyi)*rxtt+rxxy*rxr+sxxy*rxs+txxy*rxt
rxxyy    = rxi*ryi**2*rxrrr+(2.*rxi*ryi*syi+ryi**2*sxi)*rxrrs+(rxi*syi**2+2.*ryi*sxi*syi)*rxrss+sxi*syi**2*rxsss+(2.*rxi*ryi*tyi+ryi**2*txi)*rxrrt+((2.*sxi*tyi+2.*syi*txi)*ryi+2.*tyi*rxi*syi)*rxrst+(2.*sxi*syi*tyi+syi**2*txi)*rxsst+(rxi*tyi**2+2.*ryi*txi*tyi)*rxrtt+(sxi*tyi**2+2.*syi*txi*tyi)*rxstt+txi*tyi**2*rxttt+(rxi*ryy+2.*rxy*ryi)*rxrr+(rxi*syy+2.*rxy*syi+2.*ryi*sxy+ryy*sxi)*rxrs+(sxi*syy+2.*sxy*syi)*rxss+(rxi*tyy+2.*rxy*tyi+2.*ryi*txy+ryy*txi)*rxrt+(sxi*tyy+2.*sxy*tyi+2.*syi*txy+syy*txi)*rxst+(txi*tyy+2.*txy*tyi)*rxtt+rxyy*rxr+sxyy*rxs+txyy*rxt
rxyyy    = ryi**3*rxrrr+3.*ryi**2*syi*rxrrs+3.*ryi*syi**2*rxrss+syi**3*rxsss+3.*ryi**2*tyi*rxrrt+6.*ryi*syi*tyi*rxrst+3.*syi**2*tyi*rxsst+3.*ryi*tyi**2*rxrtt+3.*syi*tyi**2*rxstt+tyi**3*rxttt+3.*ryi*ryy*rxrr+(3.*ryi*syy+3.*ryy*syi)*rxrs+3.*syi*syy*rxss+(3.*ryi*tyy+3.*ryy*tyi)*rxrt+(3.*syi*tyy+3.*syy*tyi)*rxst+3.*tyi*tyy*rxtt+ryyy*rxr+syyy*rxs+tyyy*rxt
rxxxz    = rxi**2*rzi*rxrrr+(rxi**2*szi+2.*rxi*rzi*sxi)*rxrrs+(2.*rxi*sxi*szi+rzi*sxi**2)*rxrss+sxi**2*szi*rxsss+(rxi**2*tzi+2.*rxi*rzi*txi)*rxrrt+((2.*sxi*tzi+2.*szi*txi)*rxi+2.*rzi*txi*sxi)*rxrst+(sxi**2*tzi+2.*sxi*szi*txi)*rxsst+(2.*rxi*txi*tzi+rzi*txi**2)*rxrtt+(2.*sxi*txi*tzi+szi*txi**2)*rxstt+txi**2*tzi*rxttt+(2.*rxi*rxz+rxx*rzi)*rxrr+(2.*rxi*sxz+rxx*szi+2.*rxz*sxi+rzi*sxx)*rxrs+(2.*sxi*sxz+sxx*szi)*rxss+(2.*rxi*txz+rxx*tzi+2.*rxz*txi+rzi*txx)*rxrt+(2.*sxi*txz+sxx*tzi+2.*sxz*txi+szi*txx)*rxst+(2.*txi*txz+txx*tzi)*rxtt+rxxz*rxr+sxxz*rxs+txxz*rxt
rxxyz    = rxi*ryi*rzi*rxrrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*rxrrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*rxrss+sxi*syi*szi*rxsss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*rxrrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*rxrst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*rxsst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*rxrtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*rxstt+txi*tyi*tzi*rxttt+(rxi*ryz+rxy*rzi+rxz*ryi)*rxrr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*rxrs+(sxi*syz+sxy*szi+sxz*syi)*rxss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*rxrt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*rxst+(txi*tyz+txy*tzi+txz*tyi)*rxtt+rxyz*rxr+sxyz*rxs+txyz*rxt
rxyyz    = ryi**2*rzi*rxrrr+(ryi**2*szi+2.*ryi*rzi*syi)*rxrrs+(2.*ryi*syi*szi+rzi*syi**2)*rxrss+syi**2*szi*rxsss+(ryi**2*tzi+2.*ryi*rzi*tyi)*rxrrt+((2.*syi*tzi+2.*szi*tyi)*ryi+2.*rzi*tyi*syi)*rxrst+(syi**2*tzi+2.*syi*szi*tyi)*rxsst+(2.*ryi*tyi*tzi+rzi*tyi**2)*rxrtt+(2.*syi*tyi*tzi+szi*tyi**2)*rxstt+tyi**2*tzi*rxttt+(2.*ryi*ryz+ryy*rzi)*rxrr+(2.*ryi*syz+ryy*szi+2.*ryz*syi+rzi*syy)*rxrs+(2.*syi*syz+syy*szi)*rxss+(2.*ryi*tyz+ryy*tzi+2.*ryz*tyi+rzi*tyy)*rxrt+(2.*syi*tyz+syy*tzi+2.*syz*tyi+szi*tyy)*rxst+(2.*tyi*tyz+tyy*tzi)*rxtt+ryyz*rxr+syyz*rxs+tyyz*rxt
rxxzz    = rxi*rzi**2*rxrrr+(2.*rxi*rzi*szi+rzi**2*sxi)*rxrrs+(rxi*szi**2+2.*rzi*sxi*szi)*rxrss+sxi*szi**2*rxsss+(2.*rxi*rzi*tzi+rzi**2*txi)*rxrrt+((2.*sxi*tzi+2.*szi*txi)*rzi+2.*tzi*rxi*szi)*rxrst+(2.*sxi*szi*tzi+szi**2*txi)*rxsst+(rxi*tzi**2+2.*rzi*txi*tzi)*rxrtt+(sxi*tzi**2+2.*szi*txi*tzi)*rxstt+txi*tzi**2*rxttt+(rxi*rzz+2.*rxz*rzi)*rxrr+(rxi*szz+2.*rxz*szi+2.*rzi*sxz+rzz*sxi)*rxrs+(sxi*szz+2.*sxz*szi)*rxss+(rxi*tzz+2.*rxz*tzi+2.*rzi*txz+rzz*txi)*rxrt+(sxi*tzz+2.*sxz*tzi+2.*szi*txz+szz*txi)*rxst+(txi*tzz+2.*txz*tzi)*rxtt+rxzz*rxr+sxzz*rxs+txzz*rxt
rxyzz    = ryi*rzi**2*rxrrr+(2.*ryi*rzi*szi+rzi**2*syi)*rxrrs+(ryi*szi**2+2.*rzi*syi*szi)*rxrss+syi*szi**2*rxsss+(2.*ryi*rzi*tzi+rzi**2*tyi)*rxrrt+((2.*syi*tzi+2.*szi*tyi)*rzi+2.*tzi*ryi*szi)*rxrst+(2.*syi*szi*tzi+szi**2*tyi)*rxsst+(ryi*tzi**2+2.*rzi*tyi*tzi)*rxrtt+(syi*tzi**2+2.*szi*tyi*tzi)*rxstt+tyi*tzi**2*rxttt+(ryi*rzz+2.*ryz*rzi)*rxrr+(ryi*szz+2.*ryz*szi+2.*rzi*syz+rzz*syi)*rxrs+(syi*szz+2.*syz*szi)*rxss+(ryi*tzz+2.*ryz*tzi+2.*rzi*tyz+rzz*tyi)*rxrt+(syi*tzz+2.*syz*tzi+2.*szi*tyz+szz*tyi)*rxst+(tyi*tzz+2.*tyz*tzi)*rxtt+ryzz*rxr+syzz*rxs+tyzz*rxt
rxzzz    = rzi**3*rxrrr+3.*rzi**2*szi*rxrrs+3.*rzi*szi**2*rxrss+szi**3*rxsss+3.*rzi**2*tzi*rxrrt+6.*rzi*szi*tzi*rxrst+3.*szi**2*tzi*rxsst+3.*rzi*tzi**2*rxrtt+3.*szi*tzi**2*rxstt+tzi**3*rxttt+3.*rzi*rzz*rxrr+(3.*rzi*szz+3.*rzz*szi)*rxrs+3.*szi*szz*rxss+(3.*rzi*tzz+3.*rzz*tzi)*rxrt+(3.*szi*tzz+3.*szz*tzi)*rxst+3.*tzi*tzz*rxtt+rzzz*rxr+szzz*rxs+tzzz*rxt
ryyyy    = ryi**3*ryrrr+3.*ryi**2*syi*ryrrs+3.*ryi*syi**2*ryrss+syi**3*rysss+3.*ryi**2*tyi*ryrrt+6.*ryi*syi*tyi*ryrst+3.*syi**2*tyi*rysst+3.*ryi*tyi**2*ryrtt+3.*syi*tyi**2*rystt+tyi**3*ryttt+3.*ryi*ryy*ryrr+(3.*ryi*syy+3.*ryy*syi)*ryrs+3.*syi*syy*ryss+(3.*ryi*tyy+3.*ryy*tyi)*ryrt+(3.*syi*tyy+3.*syy*tyi)*ryst+3.*tyi*tyy*rytt+ryyy*ryr+syyy*rys+tyyy*ryt
ryyyz    = ryi**2*rzi*ryrrr+(ryi**2*szi+2.*ryi*rzi*syi)*ryrrs+(2.*ryi*syi*szi+rzi*syi**2)*ryrss+syi**2*szi*rysss+(ryi**2*tzi+2.*ryi*rzi*tyi)*ryrrt+((2.*syi*tzi+2.*szi*tyi)*ryi+2.*rzi*tyi*syi)*ryrst+(syi**2*tzi+2.*syi*szi*tyi)*rysst+(2.*ryi*tyi*tzi+rzi*tyi**2)*ryrtt+(2.*syi*tyi*tzi+szi*tyi**2)*rystt+tyi**2*tzi*ryttt+(2.*ryi*ryz+ryy*rzi)*ryrr+(2.*ryi*syz+ryy*szi+2.*ryz*syi+rzi*syy)*ryrs+(2.*syi*syz+syy*szi)*ryss+(2.*ryi*tyz+ryy*tzi+2.*ryz*tyi+rzi*tyy)*ryrt+(2.*syi*tyz+syy*tzi+2.*syz*tyi+szi*tyy)*ryst+(2.*tyi*tyz+tyy*tzi)*rytt+ryyz*ryr+syyz*rys+tyyz*ryt
ryyzz    = ryi*rzi**2*ryrrr+(2.*ryi*rzi*szi+rzi**2*syi)*ryrrs+(ryi*szi**2+2.*rzi*syi*szi)*ryrss+syi*szi**2*rysss+(2.*ryi*rzi*tzi+rzi**2*tyi)*ryrrt+((2.*syi*tzi+2.*szi*tyi)*rzi+2.*tzi*ryi*szi)*ryrst+(2.*syi*szi*tzi+szi**2*tyi)*rysst+(ryi*tzi**2+2.*rzi*tyi*tzi)*ryrtt+(syi*tzi**2+2.*szi*tyi*tzi)*rystt+tyi*tzi**2*ryttt+(ryi*rzz+2.*ryz*rzi)*ryrr+(ryi*szz+2.*ryz*szi+2.*rzi*syz+rzz*syi)*ryrs+(syi*szz+2.*syz*szi)*ryss+(ryi*tzz+2.*ryz*tzi+2.*rzi*tyz+rzz*tyi)*ryrt+(syi*tzz+2.*syz*tzi+2.*szi*tyz+szz*tyi)*ryst+(tyi*tzz+2.*tyz*tzi)*rytt+ryzz*ryr+syzz*rys+tyzz*ryt
ryzzz    = rzi**3*ryrrr+3.*rzi**2*szi*ryrrs+3.*rzi*szi**2*ryrss+szi**3*rysss+3.*rzi**2*tzi*ryrrt+6.*rzi*szi*tzi*ryrst+3.*szi**2*tzi*rysst+3.*rzi*tzi**2*ryrtt+3.*szi*tzi**2*rystt+tzi**3*ryttt+3.*rzi*rzz*ryrr+(3.*rzi*szz+3.*rzz*szi)*ryrs+3.*szi*szz*ryss+(3.*rzi*tzz+3.*rzz*tzi)*ryrt+(3.*szi*tzz+3.*szz*tzi)*ryst+3.*tzi*tzz*rytt+rzzz*ryr+szzz*rys+tzzz*ryt
rzzzz    = rzi**3*rzrrr+3.*rzi**2*szi*rzrrs+3.*rzi*szi**2*rzrss+szi**3*rzsss+3.*rzi**2*tzi*rzrrt+6.*rzi*szi*tzi*rzrst+3.*szi**2*tzi*rzsst+3.*rzi*tzi**2*rzrtt+3.*szi*tzi**2*rzstt+tzi**3*rzttt+3.*rzi*rzz*rzrr+(3.*rzi*szz+3.*rzz*szi)*rzrs+3.*szi*szz*rzss+(3.*rzi*tzz+3.*rzz*tzi)*rzrt+(3.*szi*tzz+3.*szz*tzi)*rzst+3.*tzi*tzz*rztt+rzzz*rzr+szzz*rzs+tzzz*rzt
sxxxx    = rxi**3*sxrrr+3.*rxi**2*sxi*sxrrs+3.*rxi*sxi**2*sxrss+sxi**3*sxsss+3.*rxi**2*txi*sxrrt+6.*rxi*sxi*txi*sxrst+3.*sxi**2*txi*sxsst+3.*rxi*txi**2*sxrtt+3.*sxi*txi**2*sxstt+txi**3*sxttt+3.*rxi*rxx*sxrr+(3.*rxi*sxx+3.*rxx*sxi)*sxrs+3.*sxi*sxx*sxss+(3.*rxi*txx+3.*rxx*txi)*sxrt+(3.*sxi*txx+3.*sxx*txi)*sxst+3.*txi*txx*sxtt+rxxx*sxr+sxxx*sxs+txxx*sxt
sxxxy    = rxi**2*ryi*sxrrr+(rxi**2*syi+2.*rxi*ryi*sxi)*sxrrs+(2.*rxi*sxi*syi+ryi*sxi**2)*sxrss+sxi**2*syi*sxsss+(rxi**2*tyi+2.*rxi*ryi*txi)*sxrrt+((2.*sxi*tyi+2.*syi*txi)*rxi+2.*ryi*txi*sxi)*sxrst+(sxi**2*tyi+2.*sxi*syi*txi)*sxsst+(2.*rxi*txi*tyi+ryi*txi**2)*sxrtt+(2.*sxi*txi*tyi+syi*txi**2)*sxstt+txi**2*tyi*sxttt+(2.*rxi*rxy+rxx*ryi)*sxrr+(2.*rxi*sxy+rxx*syi+2.*rxy*sxi+ryi*sxx)*sxrs+(2.*sxi*sxy+sxx*syi)*sxss+(2.*rxi*txy+rxx*tyi+2.*rxy*txi+ryi*txx)*sxrt+(2.*sxi*txy+sxx*tyi+2.*sxy*txi+syi*txx)*sxst+(2.*txi*txy+txx*tyi)*sxtt+rxxy*sxr+sxxy*sxs+txxy*sxt
sxxyy    = rxi*ryi**2*sxrrr+(2.*rxi*ryi*syi+ryi**2*sxi)*sxrrs+(rxi*syi**2+2.*ryi*sxi*syi)*sxrss+sxi*syi**2*sxsss+(2.*rxi*ryi*tyi+ryi**2*txi)*sxrrt+((2.*sxi*tyi+2.*syi*txi)*ryi+2.*tyi*rxi*syi)*sxrst+(2.*sxi*syi*tyi+syi**2*txi)*sxsst+(rxi*tyi**2+2.*ryi*txi*tyi)*sxrtt+(sxi*tyi**2+2.*syi*txi*tyi)*sxstt+txi*tyi**2*sxttt+(rxi*ryy+2.*rxy*ryi)*sxrr+(rxi*syy+2.*rxy*syi+2.*ryi*sxy+ryy*sxi)*sxrs+(sxi*syy+2.*sxy*syi)*sxss+(rxi*tyy+2.*rxy*tyi+2.*ryi*txy+ryy*txi)*sxrt+(sxi*tyy+2.*sxy*tyi+2.*syi*txy+syy*txi)*sxst+(txi*tyy+2.*txy*tyi)*sxtt+rxyy*sxr+sxyy*sxs+txyy*sxt
sxyyy    = ryi**3*sxrrr+3.*ryi**2*syi*sxrrs+3.*ryi*syi**2*sxrss+syi**3*sxsss+3.*ryi**2*tyi*sxrrt+6.*ryi*syi*tyi*sxrst+3.*syi**2*tyi*sxsst+3.*ryi*tyi**2*sxrtt+3.*syi*tyi**2*sxstt+tyi**3*sxttt+3.*ryi*ryy*sxrr+(3.*ryi*syy+3.*ryy*syi)*sxrs+3.*syi*syy*sxss+(3.*ryi*tyy+3.*ryy*tyi)*sxrt+(3.*syi*tyy+3.*syy*tyi)*sxst+3.*tyi*tyy*sxtt+ryyy*sxr+syyy*sxs+tyyy*sxt
sxxxz    = rxi**2*rzi*sxrrr+(rxi**2*szi+2.*rxi*rzi*sxi)*sxrrs+(2.*rxi*sxi*szi+rzi*sxi**2)*sxrss+sxi**2*szi*sxsss+(rxi**2*tzi+2.*rxi*rzi*txi)*sxrrt+((2.*sxi*tzi+2.*szi*txi)*rxi+2.*rzi*txi*sxi)*sxrst+(sxi**2*tzi+2.*sxi*szi*txi)*sxsst+(2.*rxi*txi*tzi+rzi*txi**2)*sxrtt+(2.*sxi*txi*tzi+szi*txi**2)*sxstt+txi**2*tzi*sxttt+(2.*rxi*rxz+rxx*rzi)*sxrr+(2.*rxi*sxz+rxx*szi+2.*rxz*sxi+rzi*sxx)*sxrs+(2.*sxi*sxz+sxx*szi)*sxss+(2.*rxi*txz+rxx*tzi+2.*rxz*txi+rzi*txx)*sxrt+(2.*sxi*txz+sxx*tzi+2.*sxz*txi+szi*txx)*sxst+(2.*txi*txz+txx*tzi)*sxtt+rxxz*sxr+sxxz*sxs+txxz*sxt
sxxyz    = rxi*ryi*rzi*sxrrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*sxrrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*sxrss+sxi*syi*szi*sxsss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*sxrrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*sxrst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*sxsst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*sxrtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*sxstt+txi*tyi*tzi*sxttt+(rxi*ryz+rxy*rzi+rxz*ryi)*sxrr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*sxrs+(sxi*syz+sxy*szi+sxz*syi)*sxss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*sxrt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*sxst+(txi*tyz+txy*tzi+txz*tyi)*sxtt+rxyz*sxr+sxyz*sxs+txyz*sxt
sxyyz    = ryi**2*rzi*sxrrr+(ryi**2*szi+2.*ryi*rzi*syi)*sxrrs+(2.*ryi*syi*szi+rzi*syi**2)*sxrss+syi**2*szi*sxsss+(ryi**2*tzi+2.*ryi*rzi*tyi)*sxrrt+((2.*syi*tzi+2.*szi*tyi)*ryi+2.*rzi*tyi*syi)*sxrst+(syi**2*tzi+2.*syi*szi*tyi)*sxsst+(2.*ryi*tyi*tzi+rzi*tyi**2)*sxrtt+(2.*syi*tyi*tzi+szi*tyi**2)*sxstt+tyi**2*tzi*sxttt+(2.*ryi*ryz+ryy*rzi)*sxrr+(2.*ryi*syz+ryy*szi+2.*ryz*syi+rzi*syy)*sxrs+(2.*syi*syz+syy*szi)*sxss+(2.*ryi*tyz+ryy*tzi+2.*ryz*tyi+rzi*tyy)*sxrt+(2.*syi*tyz+syy*tzi+2.*syz*tyi+szi*tyy)*sxst+(2.*tyi*tyz+tyy*tzi)*sxtt+ryyz*sxr+syyz*sxs+tyyz*sxt
sxxzz    = rxi*rzi**2*sxrrr+(2.*rxi*rzi*szi+rzi**2*sxi)*sxrrs+(rxi*szi**2+2.*rzi*sxi*szi)*sxrss+sxi*szi**2*sxsss+(2.*rxi*rzi*tzi+rzi**2*txi)*sxrrt+((2.*sxi*tzi+2.*szi*txi)*rzi+2.*tzi*rxi*szi)*sxrst+(2.*sxi*szi*tzi+szi**2*txi)*sxsst+(rxi*tzi**2+2.*rzi*txi*tzi)*sxrtt+(sxi*tzi**2+2.*szi*txi*tzi)*sxstt+txi*tzi**2*sxttt+(rxi*rzz+2.*rxz*rzi)*sxrr+(rxi*szz+2.*rxz*szi+2.*rzi*sxz+rzz*sxi)*sxrs+(sxi*szz+2.*sxz*szi)*sxss+(rxi*tzz+2.*rxz*tzi+2.*rzi*txz+rzz*txi)*sxrt+(sxi*tzz+2.*sxz*tzi+2.*szi*txz+szz*txi)*sxst+(txi*tzz+2.*txz*tzi)*sxtt+rxzz*sxr+sxzz*sxs+txzz*sxt
sxyzz    = ryi*rzi**2*sxrrr+(2.*ryi*rzi*szi+rzi**2*syi)*sxrrs+(ryi*szi**2+2.*rzi*syi*szi)*sxrss+syi*szi**2*sxsss+(2.*ryi*rzi*tzi+rzi**2*tyi)*sxrrt+((2.*syi*tzi+2.*szi*tyi)*rzi+2.*tzi*ryi*szi)*sxrst+(2.*syi*szi*tzi+szi**2*tyi)*sxsst+(ryi*tzi**2+2.*rzi*tyi*tzi)*sxrtt+(syi*tzi**2+2.*szi*tyi*tzi)*sxstt+tyi*tzi**2*sxttt+(ryi*rzz+2.*ryz*rzi)*sxrr+(ryi*szz+2.*ryz*szi+2.*rzi*syz+rzz*syi)*sxrs+(syi*szz+2.*syz*szi)*sxss+(ryi*tzz+2.*ryz*tzi+2.*rzi*tyz+rzz*tyi)*sxrt+(syi*tzz+2.*syz*tzi+2.*szi*tyz+szz*tyi)*sxst+(tyi*tzz+2.*tyz*tzi)*sxtt+ryzz*sxr+syzz*sxs+tyzz*sxt
sxzzz    = rzi**3*sxrrr+3.*rzi**2*szi*sxrrs+3.*rzi*szi**2*sxrss+szi**3*sxsss+3.*rzi**2*tzi*sxrrt+6.*rzi*szi*tzi*sxrst+3.*szi**2*tzi*sxsst+3.*rzi*tzi**2*sxrtt+3.*szi*tzi**2*sxstt+tzi**3*sxttt+3.*rzi*rzz*sxrr+(3.*rzi*szz+3.*rzz*szi)*sxrs+3.*szi*szz*sxss+(3.*rzi*tzz+3.*rzz*tzi)*sxrt+(3.*szi*tzz+3.*szz*tzi)*sxst+3.*tzi*tzz*sxtt+rzzz*sxr+szzz*sxs+tzzz*sxt
syyyy    = ryi**3*syrrr+3.*ryi**2*syi*syrrs+3.*ryi*syi**2*syrss+syi**3*sysss+3.*ryi**2*tyi*syrrt+6.*ryi*syi*tyi*syrst+3.*syi**2*tyi*sysst+3.*ryi*tyi**2*syrtt+3.*syi*tyi**2*systt+tyi**3*syttt+3.*ryi*ryy*syrr+(3.*ryi*syy+3.*ryy*syi)*syrs+3.*syi*syy*syss+(3.*ryi*tyy+3.*ryy*tyi)*syrt+(3.*syi*tyy+3.*syy*tyi)*syst+3.*tyi*tyy*sytt+ryyy*syr+syyy*sys+tyyy*syt
syyyz    = ryi**2*rzi*syrrr+(ryi**2*szi+2.*ryi*rzi*syi)*syrrs+(2.*ryi*syi*szi+rzi*syi**2)*syrss+syi**2*szi*sysss+(ryi**2*tzi+2.*ryi*rzi*tyi)*syrrt+((2.*syi*tzi+2.*szi*tyi)*ryi+2.*rzi*tyi*syi)*syrst+(syi**2*tzi+2.*syi*szi*tyi)*sysst+(2.*ryi*tyi*tzi+rzi*tyi**2)*syrtt+(2.*syi*tyi*tzi+szi*tyi**2)*systt+tyi**2*tzi*syttt+(2.*ryi*ryz+ryy*rzi)*syrr+(2.*ryi*syz+ryy*szi+2.*ryz*syi+rzi*syy)*syrs+(2.*syi*syz+syy*szi)*syss+(2.*ryi*tyz+ryy*tzi+2.*ryz*tyi+rzi*tyy)*syrt+(2.*syi*tyz+syy*tzi+2.*syz*tyi+szi*tyy)*syst+(2.*tyi*tyz+tyy*tzi)*sytt+ryyz*syr+syyz*sys+tyyz*syt
syyzz    = ryi*rzi**2*syrrr+(2.*ryi*rzi*szi+rzi**2*syi)*syrrs+(ryi*szi**2+2.*rzi*syi*szi)*syrss+syi*szi**2*sysss+(2.*ryi*rzi*tzi+rzi**2*tyi)*syrrt+((2.*syi*tzi+2.*szi*tyi)*rzi+2.*tzi*ryi*szi)*syrst+(2.*syi*szi*tzi+szi**2*tyi)*sysst+(ryi*tzi**2+2.*rzi*tyi*tzi)*syrtt+(syi*tzi**2+2.*szi*tyi*tzi)*systt+tyi*tzi**2*syttt+(ryi*rzz+2.*ryz*rzi)*syrr+(ryi*szz+2.*ryz*szi+2.*rzi*syz+rzz*syi)*syrs+(syi*szz+2.*syz*szi)*syss+(ryi*tzz+2.*ryz*tzi+2.*rzi*tyz+rzz*tyi)*syrt+(syi*tzz+2.*syz*tzi+2.*szi*tyz+szz*tyi)*syst+(tyi*tzz+2.*tyz*tzi)*sytt+ryzz*syr+syzz*sys+tyzz*syt
syzzz    = rzi**3*syrrr+3.*rzi**2*szi*syrrs+3.*rzi*szi**2*syrss+szi**3*sysss+3.*rzi**2*tzi*syrrt+6.*rzi*szi*tzi*syrst+3.*szi**2*tzi*sysst+3.*rzi*tzi**2*syrtt+3.*szi*tzi**2*systt+tzi**3*syttt+3.*rzi*rzz*syrr+(3.*rzi*szz+3.*rzz*szi)*syrs+3.*szi*szz*syss+(3.*rzi*tzz+3.*rzz*tzi)*syrt+(3.*szi*tzz+3.*szz*tzi)*syst+3.*tzi*tzz*sytt+rzzz*syr+szzz*sys+tzzz*syt
szzzz    = rzi**3*szrrr+3.*rzi**2*szi*szrrs+3.*rzi*szi**2*szrss+szi**3*szsss+3.*rzi**2*tzi*szrrt+6.*rzi*szi*tzi*szrst+3.*szi**2*tzi*szsst+3.*rzi*tzi**2*szrtt+3.*szi*tzi**2*szstt+tzi**3*szttt+3.*rzi*rzz*szrr+(3.*rzi*szz+3.*rzz*szi)*szrs+3.*szi*szz*szss+(3.*rzi*tzz+3.*rzz*tzi)*szrt+(3.*szi*tzz+3.*szz*tzi)*szst+3.*tzi*tzz*sztt+rzzz*szr+szzz*szs+tzzz*szt
txxxx    = rxi**3*txrrr+3.*rxi**2*sxi*txrrs+3.*rxi*sxi**2*txrss+sxi**3*txsss+3.*rxi**2*txi*txrrt+6.*rxi*sxi*txi*txrst+3.*sxi**2*txi*txsst+3.*rxi*txi**2*txrtt+3.*sxi*txi**2*txstt+txi**3*txttt+3.*rxi*rxx*txrr+(3.*rxi*sxx+3.*rxx*sxi)*txrs+3.*sxi*sxx*txss+(3.*rxi*txx+3.*rxx*txi)*txrt+(3.*sxi*txx+3.*sxx*txi)*txst+3.*txi*txx*txtt+rxxx*txr+sxxx*txs+txxx*txt
txxxy    = rxi**2*ryi*txrrr+(rxi**2*syi+2.*rxi*ryi*sxi)*txrrs+(2.*rxi*sxi*syi+ryi*sxi**2)*txrss+sxi**2*syi*txsss+(rxi**2*tyi+2.*rxi*ryi*txi)*txrrt+((2.*sxi*tyi+2.*syi*txi)*rxi+2.*ryi*txi*sxi)*txrst+(sxi**2*tyi+2.*sxi*syi*txi)*txsst+(2.*rxi*txi*tyi+ryi*txi**2)*txrtt+(2.*sxi*txi*tyi+syi*txi**2)*txstt+txi**2*tyi*txttt+(2.*rxi*rxy+rxx*ryi)*txrr+(2.*rxi*sxy+rxx*syi+2.*rxy*sxi+ryi*sxx)*txrs+(2.*sxi*sxy+sxx*syi)*txss+(2.*rxi*txy+rxx*tyi+2.*rxy*txi+ryi*txx)*txrt+(2.*sxi*txy+sxx*tyi+2.*sxy*txi+syi*txx)*txst+(2.*txi*txy+txx*tyi)*txtt+rxxy*txr+sxxy*txs+txxy*txt
txxyy    = rxi*ryi**2*txrrr+(2.*rxi*ryi*syi+ryi**2*sxi)*txrrs+(rxi*syi**2+2.*ryi*sxi*syi)*txrss+sxi*syi**2*txsss+(2.*rxi*ryi*tyi+ryi**2*txi)*txrrt+((2.*sxi*tyi+2.*syi*txi)*ryi+2.*tyi*rxi*syi)*txrst+(2.*sxi*syi*tyi+syi**2*txi)*txsst+(rxi*tyi**2+2.*ryi*txi*tyi)*txrtt+(sxi*tyi**2+2.*syi*txi*tyi)*txstt+txi*tyi**2*txttt+(rxi*ryy+2.*rxy*ryi)*txrr+(rxi*syy+2.*rxy*syi+2.*ryi*sxy+ryy*sxi)*txrs+(sxi*syy+2.*sxy*syi)*txss+(rxi*tyy+2.*rxy*tyi+2.*ryi*txy+ryy*txi)*txrt+(sxi*tyy+2.*sxy*tyi+2.*syi*txy+syy*txi)*txst+(txi*tyy+2.*txy*tyi)*txtt+rxyy*txr+sxyy*txs+txyy*txt
txyyy    = ryi**3*txrrr+3.*ryi**2*syi*txrrs+3.*ryi*syi**2*txrss+syi**3*txsss+3.*ryi**2*tyi*txrrt+6.*ryi*syi*tyi*txrst+3.*syi**2*tyi*txsst+3.*ryi*tyi**2*txrtt+3.*syi*tyi**2*txstt+tyi**3*txttt+3.*ryi*ryy*txrr+(3.*ryi*syy+3.*ryy*syi)*txrs+3.*syi*syy*txss+(3.*ryi*tyy+3.*ryy*tyi)*txrt+(3.*syi*tyy+3.*syy*tyi)*txst+3.*tyi*tyy*txtt+ryyy*txr+syyy*txs+tyyy*txt
txxxz    = rxi**2*rzi*txrrr+(rxi**2*szi+2.*rxi*rzi*sxi)*txrrs+(2.*rxi*sxi*szi+rzi*sxi**2)*txrss+sxi**2*szi*txsss+(rxi**2*tzi+2.*rxi*rzi*txi)*txrrt+((2.*sxi*tzi+2.*szi*txi)*rxi+2.*rzi*txi*sxi)*txrst+(sxi**2*tzi+2.*sxi*szi*txi)*txsst+(2.*rxi*txi*tzi+rzi*txi**2)*txrtt+(2.*sxi*txi*tzi+szi*txi**2)*txstt+txi**2*tzi*txttt+(2.*rxi*rxz+rxx*rzi)*txrr+(2.*rxi*sxz+rxx*szi+2.*rxz*sxi+rzi*sxx)*txrs+(2.*sxi*sxz+sxx*szi)*txss+(2.*rxi*txz+rxx*tzi+2.*rxz*txi+rzi*txx)*txrt+(2.*sxi*txz+sxx*tzi+2.*sxz*txi+szi*txx)*txst+(2.*txi*txz+txx*tzi)*txtt+rxxz*txr+sxxz*txs+txxz*txt
txxyz    = rxi*ryi*rzi*txrrr+((ryi*szi+rzi*syi)*rxi+rzi*sxi*ryi)*txrrs+(rxi*syi*szi+ryi*sxi*szi+rzi*sxi*syi)*txrss+sxi*syi*szi*txsss+((ryi*tzi+rzi*tyi)*rxi+rzi*txi*ryi)*txrrt+((syi*tzi+szi*tyi)*rxi+(sxi*tzi+szi*txi)*ryi+(sxi*tyi+syi*txi)*rzi)*txrst+((syi*tzi+szi*tyi)*sxi+szi*txi*syi)*txsst+(rxi*tyi*tzi+ryi*txi*tzi+rzi*txi*tyi)*txrtt+(sxi*tyi*tzi+syi*txi*tzi+szi*txi*tyi)*txstt+txi*tyi*tzi*txttt+(rxi*ryz+rxy*rzi+rxz*ryi)*txrr+(rxi*syz+rxy*szi+rxz*syi+ryi*sxz+ryz*sxi+rzi*sxy)*txrs+(sxi*syz+sxy*szi+sxz*syi)*txss+(rxi*tyz+rxy*tzi+rxz*tyi+ryi*txz+ryz*txi+rzi*txy)*txrt+(sxi*tyz+sxy*tzi+sxz*tyi+syi*txz+syz*txi+szi*txy)*txst+(txi*tyz+txy*tzi+txz*tyi)*txtt+rxyz*txr+sxyz*txs+txyz*txt
txyyz    = ryi**2*rzi*txrrr+(ryi**2*szi+2.*ryi*rzi*syi)*txrrs+(2.*ryi*syi*szi+rzi*syi**2)*txrss+syi**2*szi*txsss+(ryi**2*tzi+2.*ryi*rzi*tyi)*txrrt+((2.*syi*tzi+2.*szi*tyi)*ryi+2.*rzi*tyi*syi)*txrst+(syi**2*tzi+2.*syi*szi*tyi)*txsst+(2.*ryi*tyi*tzi+rzi*tyi**2)*txrtt+(2.*syi*tyi*tzi+szi*tyi**2)*txstt+tyi**2*tzi*txttt+(2.*ryi*ryz+ryy*rzi)*txrr+(2.*ryi*syz+ryy*szi+2.*ryz*syi+rzi*syy)*txrs+(2.*syi*syz+syy*szi)*txss+(2.*ryi*tyz+ryy*tzi+2.*ryz*tyi+rzi*tyy)*txrt+(2.*syi*tyz+syy*tzi+2.*syz*tyi+szi*tyy)*txst+(2.*tyi*tyz+tyy*tzi)*txtt+ryyz*txr+syyz*txs+tyyz*txt
txxzz    = rxi*rzi**2*txrrr+(2.*rxi*rzi*szi+rzi**2*sxi)*txrrs+(rxi*szi**2+2.*rzi*sxi*szi)*txrss+sxi*szi**2*txsss+(2.*rxi*rzi*tzi+rzi**2*txi)*txrrt+((2.*sxi*tzi+2.*szi*txi)*rzi+2.*tzi*rxi*szi)*txrst+(2.*sxi*szi*tzi+szi**2*txi)*txsst+(rxi*tzi**2+2.*rzi*txi*tzi)*txrtt+(sxi*tzi**2+2.*szi*txi*tzi)*txstt+txi*tzi**2*txttt+(rxi*rzz+2.*rxz*rzi)*txrr+(rxi*szz+2.*rxz*szi+2.*rzi*sxz+rzz*sxi)*txrs+(sxi*szz+2.*sxz*szi)*txss+(rxi*tzz+2.*rxz*tzi+2.*rzi*txz+rzz*txi)*txrt+(sxi*tzz+2.*sxz*tzi+2.*szi*txz+szz*txi)*txst+(txi*tzz+2.*txz*tzi)*txtt+rxzz*txr+sxzz*txs+txzz*txt
txyzz    = ryi*rzi**2*txrrr+(2.*ryi*rzi*szi+rzi**2*syi)*txrrs+(ryi*szi**2+2.*rzi*syi*szi)*txrss+syi*szi**2*txsss+(2.*ryi*rzi*tzi+rzi**2*tyi)*txrrt+((2.*syi*tzi+2.*szi*tyi)*rzi+2.*tzi*ryi*szi)*txrst+(2.*syi*szi*tzi+szi**2*tyi)*txsst+(ryi*tzi**2+2.*rzi*tyi*tzi)*txrtt+(syi*tzi**2+2.*szi*tyi*tzi)*txstt+tyi*tzi**2*txttt+(ryi*rzz+2.*ryz*rzi)*txrr+(ryi*szz+2.*ryz*szi+2.*rzi*syz+rzz*syi)*txrs+(syi*szz+2.*syz*szi)*txss+(ryi*tzz+2.*ryz*tzi+2.*rzi*tyz+rzz*tyi)*txrt+(syi*tzz+2.*syz*tzi+2.*szi*tyz+szz*tyi)*txst+(tyi*tzz+2.*tyz*tzi)*txtt+ryzz*txr+syzz*txs+tyzz*txt
txzzz    = rzi**3*txrrr+3.*rzi**2*szi*txrrs+3.*rzi*szi**2*txrss+szi**3*txsss+3.*rzi**2*tzi*txrrt+6.*rzi*szi*tzi*txrst+3.*szi**2*tzi*txsst+3.*rzi*tzi**2*txrtt+3.*szi*tzi**2*txstt+tzi**3*txttt+3.*rzi*rzz*txrr+(3.*rzi*szz+3.*rzz*szi)*txrs+3.*szi*szz*txss+(3.*rzi*tzz+3.*rzz*tzi)*txrt+(3.*szi*tzz+3.*szz*tzi)*txst+3.*tzi*tzz*txtt+rzzz*txr+szzz*txs+tzzz*txt
tyyyy    = ryi**3*tyrrr+3.*ryi**2*syi*tyrrs+3.*ryi*syi**2*tyrss+syi**3*tysss+3.*ryi**2*tyi*tyrrt+6.*ryi*syi*tyi*tyrst+3.*syi**2*tyi*tysst+3.*ryi*tyi**2*tyrtt+3.*syi*tyi**2*tystt+tyi**3*tyttt+3.*ryi*ryy*tyrr+(3.*ryi*syy+3.*ryy*syi)*tyrs+3.*syi*syy*tyss+(3.*ryi*tyy+3.*ryy*tyi)*tyrt+(3.*syi*tyy+3.*syy*tyi)*tyst+3.*tyi*tyy*tytt+ryyy*tyr+syyy*tys+tyyy*tyt
tyyyz    = ryi**2*rzi*tyrrr+(ryi**2*szi+2.*ryi*rzi*syi)*tyrrs+(2.*ryi*syi*szi+rzi*syi**2)*tyrss+syi**2*szi*tysss+(ryi**2*tzi+2.*ryi*rzi*tyi)*tyrrt+((2.*syi*tzi+2.*szi*tyi)*ryi+2.*rzi*tyi*syi)*tyrst+(syi**2*tzi+2.*syi*szi*tyi)*tysst+(2.*ryi*tyi*tzi+rzi*tyi**2)*tyrtt+(2.*syi*tyi*tzi+szi*tyi**2)*tystt+tyi**2*tzi*tyttt+(2.*ryi*ryz+ryy*rzi)*tyrr+(2.*ryi*syz+ryy*szi+2.*ryz*syi+rzi*syy)*tyrs+(2.*syi*syz+syy*szi)*tyss+(2.*ryi*tyz+ryy*tzi+2.*ryz*tyi+rzi*tyy)*tyrt+(2.*syi*tyz+syy*tzi+2.*syz*tyi+szi*tyy)*tyst+(2.*tyi*tyz+tyy*tzi)*tytt+ryyz*tyr+syyz*tys+tyyz*tyt
tyyzz    = ryi*rzi**2*tyrrr+(2.*ryi*rzi*szi+rzi**2*syi)*tyrrs+(ryi*szi**2+2.*rzi*syi*szi)*tyrss+syi*szi**2*tysss+(2.*ryi*rzi*tzi+rzi**2*tyi)*tyrrt+((2.*syi*tzi+2.*szi*tyi)*rzi+2.*tzi*ryi*szi)*tyrst+(2.*syi*szi*tzi+szi**2*tyi)*tysst+(ryi*tzi**2+2.*rzi*tyi*tzi)*tyrtt+(syi*tzi**2+2.*szi*tyi*tzi)*tystt+tyi*tzi**2*tyttt+(ryi*rzz+2.*ryz*rzi)*tyrr+(ryi*szz+2.*ryz*szi+2.*rzi*syz+rzz*syi)*tyrs+(syi*szz+2.*syz*szi)*tyss+(ryi*tzz+2.*ryz*tzi+2.*rzi*tyz+rzz*tyi)*tyrt+(syi*tzz+2.*syz*tzi+2.*szi*tyz+szz*tyi)*tyst+(tyi*tzz+2.*tyz*tzi)*tytt+ryzz*tyr+syzz*tys+tyzz*tyt
tyzzz    = rzi**3*tyrrr+3.*rzi**2*szi*tyrrs+3.*rzi*szi**2*tyrss+szi**3*tysss+3.*rzi**2*tzi*tyrrt+6.*rzi*szi*tzi*tyrst+3.*szi**2*tzi*tysst+3.*rzi*tzi**2*tyrtt+3.*szi*tzi**2*tystt+tzi**3*tyttt+3.*rzi*rzz*tyrr+(3.*rzi*szz+3.*rzz*szi)*tyrs+3.*szi*szz*tyss+(3.*rzi*tzz+3.*rzz*tzi)*tyrt+(3.*szi*tzz+3.*szz*tzi)*tyst+3.*tzi*tzz*tytt+rzzz*tyr+szzz*tys+tzzz*tyt
tzzzz    = rzi**3*tzrrr+3.*rzi**2*szi*tzrrs+3.*rzi*szi**2*tzrss+szi**3*tzsss+3.*rzi**2*tzi*tzrrt+6.*rzi*szi*tzi*tzrst+3.*szi**2*tzi*tzsst+3.*rzi*tzi**2*tzrtt+3.*szi*tzi**2*tzstt+tzi**3*tzttt+3.*rzi*rzz*tzrr+(3.*rzi*szz+3.*rzz*szi)*tzrs+3.*szi*szz*tzss+(3.*rzi*tzz+3.*rzz*tzi)*tzrt+(3.*szi*tzz+3.*szz*tzi)*tzst+3.*tzi*tzz*tztt+rzzz*tzr+szzz*tzs+tzzz*tzt
#End
! ---- end OPTION eq evalMetrics ---

! ---------- Fourth spatial derivatives of u ---------
uxx      = rxi**2*urr+2.*rxi*sxi*urs+2.*rxi*txi*urt+sxi**2*uss+2.*sxi*txi*ust+txi**2*utt+rxx*ur+sxx*us+txx*ut
uxxxx    = 12.*rxi**2*sxi*txi*urrst+12.*rxi*sxi**2*txi*ursst+12.*rxi*sxi*txi**2*urstt+txi**4*utttt+(6.*rxi**2*txx+12.*rxi*rxx*txi)*urrt+((12.*sxi*txx+12.*sxx*txi)*rxi+12.*rxx*txi*sxi)*urst+(6.*sxi**2*txx+12.*sxi*sxx*txi)*usst+(12.*rxi*txi*txx+6.*rxx*txi**2)*urtt+(12.*sxi*txi*txx+6.*sxx*txi**2)*ustt+(4.*rxi*txxx+6.*rxx*txx+4.*rxxx*txi)*urt+(4.*sxi*txxx+6.*sxx*txx+4.*sxxx*txi)*ust+(4.*txi*txxx+3.*txx**2)*utt+txxxx*ut+4.*rxi**3*txi*urrrt+4.*sxi**3*txi*ussst+6.*rxi**2*txi**2*urrtt+6.*sxi**2*txi**2*usstt+4.*rxi*txi**3*urttt+4.*sxi*txi**3*usttt+6.*txi**2*txx*uttt+rxi**4*urrrr+sxi**4*ussss+(6.*rxi**2*sxx+12.*rxi*rxx*sxi)*urrs+(12.*rxi*sxi*sxx+6.*rxx*sxi**2)*urss+(4.*rxi*rxxx+3.*rxx**2)*urr+(4.*rxi*sxxx+6.*rxx*sxx+4.*rxxx*sxi)*urs+(4.*sxi*sxxx+3.*sxx**2)*uss+rxxxx*ur+sxxxx*us+4.*rxi**3*sxi*urrrs+6.*rxi**2*sxi**2*urrss+4.*rxi*sxi**3*ursss+6.*rxi**2*rxx*urrr+6.*sxi**2*sxx*usss
uyy      = ryi**2*urr+2.*ryi*syi*urs+2.*ryi*tyi*urt+syi**2*uss+2.*syi*tyi*ust+tyi**2*utt+ryy*ur+syy*us+tyy*ut
uxxyy    = txi**2*tyi**2*utttt+(rxi**2*tyy+ryi**2*txx+(4.*rxy*tyi+4.*ryi*txy+2.*ryy*txi)*rxi+(2.*rxx*tyi+4.*rxy*txi)*ryi)*urrt+((4.*rxy*tyi+2.*ryy*txi)*sxi+(2.*sxi*tyy+4.*sxy*tyi+4.*syi*txy+2.*syy*txi)*rxi+(4.*sxi*txy+2.*sxx*tyi+4.*sxy*txi+2.*syi*txx)*ryi+(2.*rxx*tyi+4.*rxy*txi)*syi)*urst+(sxi**2*tyy+txx*syi**2+(4.*sxy*tyi+4.*syi*txy+2.*syy*txi)*sxi+(2.*sxx*tyi+4.*sxy*txi)*syi)*usst+(rxx*tyi**2+4.*rxy*tyi*txi+(2.*txi*tyy+4.*txy*tyi)*rxi+(4.*txi*txy+2.*txx*tyi)*ryi+ryy*txi**2)*urtt+(sxx*tyi**2+4.*sxy*tyi*txi+(2.*txi*tyy+4.*txy*tyi)*sxi+(4.*txi*txy+2.*txx*tyi)*syi+syy*txi**2)*ustt+(txi**2*tyy+4.*txi*txy*tyi+txx*tyi**2)*uttt+(2.*rxi*txyy+rxx*tyy+2.*rxxy*tyi+4.*rxy*txy+2.*rxyy*txi+2.*ryi*txxy+ryy*txx)*urt+(2.*sxi*txyy+sxx*tyy+2.*sxxy*tyi+4.*sxy*txy+2.*sxyy*txi+2.*syi*txxy+syy*txx)*ust+(2.*txi*txyy+txx*tyy+2.*txxy*tyi+2.*txy**2)*utt+txxyy*ut+(2.*rxi**2*ryi*tyi+2.*rxi*ryi**2*txi)*urrrt+(2.*syi*tyi*rxi**2+(4.*sxi*tyi+4.*syi*txi)*ryi*rxi+2.*sxi*txi*ryi**2)*urrst+((4.*sxi*syi*tyi+2.*syi**2*txi)*rxi+(2.*sxi**2*tyi+4.*sxi*syi*txi)*ryi)*ursst+(2.*sxi**2*syi*tyi+2.*sxi*syi**2*txi)*ussst+(rxi**2*tyi**2+4.*rxi*ryi*txi*tyi+ryi**2*txi**2)*urrtt+((2.*sxi*tyi**2+4.*syi*txi*tyi)*rxi+(4.*sxi*txi*tyi+2.*syi*txi**2)*ryi)*urstt+(sxi**2*tyi**2+4.*sxi*syi*txi*tyi+syi**2*txi**2)*usstt+(2.*rxi*txi*tyi**2+2.*ryi*txi**2*tyi)*urttt+(2.*sxi*txi*tyi**2+2.*syi*txi**2*tyi)*usttt+(2.*rxi*sxyy+rxx*syy+2.*rxxy*syi+4.*rxy*sxy+2.*rxyy*sxi+2.*ryi*sxxy+ryy*sxx)*urs+(2.*sxi*sxyy+sxx*syy+2.*sxxy*syi+2.*sxy**2)*uss+rxxyy*ur+sxxyy*us+(2.*rxi**2*ryi*syi+2.*rxi*ryi**2*sxi)*urrrs+(rxi**2*syi**2+4.*rxi*ryi*sxi*syi+ryi**2*sxi**2)*urrss+(2.*rxi*sxi*syi**2+2.*ryi*sxi**2*syi)*ursss+(rxi**2*ryy+4.*rxi*rxy*ryi+rxx*ryi**2)*urrr+(rxi**2*syy+(4.*rxy*syi+4.*ryi*sxy+2.*ryy*sxi)*rxi+ryi**2*sxx+(2.*rxx*syi+4.*rxy*sxi)*ryi)*urrs+((2.*sxi*syy+4.*sxy*syi)*rxi+(4.*sxi*sxy+2.*sxx*syi)*ryi+ryy*sxi**2+4.*rxy*syi*sxi+rxx*syi**2)*urss+(sxi**2*syy+4.*sxi*sxy*syi+sxx*syi**2)*usss+(2.*rxi*rxyy+rxx*ryy+2.*rxxy*ryi+2.*rxy**2)*urr+rxi**2*ryi**2*urrrr+sxi**2*syi**2*ussss
uyyyy    = 12.*ryi**2*syi*tyi*urrst+12.*ryi*syi**2*tyi*ursst+12.*ryi*syi*tyi**2*urstt+tyi**4*utttt+(6.*ryi**2*tyy+12.*ryi*ryy*tyi)*urrt+((12.*syi*tyy+12.*syy*tyi)*ryi+12.*ryy*tyi*syi)*urst+(6.*syi**2*tyy+12.*syi*syy*tyi)*usst+(12.*ryi*tyi*tyy+6.*ryy*tyi**2)*urtt+(12.*syi*tyi*tyy+6.*syy*tyi**2)*ustt+(4.*ryi*tyyy+6.*ryy*tyy+4.*ryyy*tyi)*urt+(4.*syi*tyyy+6.*syy*tyy+4.*syyy*tyi)*ust+(4.*tyi*tyyy+3.*tyy**2)*utt+tyyyy*ut+4.*ryi**3*tyi*urrrt+4.*syi**3*tyi*ussst+6.*ryi**2*tyi**2*urrtt+6.*syi**2*tyi**2*usstt+4.*ryi*tyi**3*urttt+4.*syi*tyi**3*usttt+6.*tyi**2*tyy*uttt+ryi**4*urrrr+syi**4*ussss+(6.*ryi**2*syy+12.*ryi*ryy*syi)*urrs+(12.*ryi*syi*syy+6.*ryy*syi**2)*urss+(4.*ryi*ryyy+3.*ryy**2)*urr+(4.*ryi*syyy+6.*ryy*syy+4.*ryyy*syi)*urs+(4.*syi*syyy+3.*syy**2)*uss+ryyyy*ur+syyyy*us+4.*ryi**3*syi*urrrs+6.*ryi**2*syi**2*urrss+4.*ryi*syi**3*ursss+6.*ryi**2*ryy*urrr+6.*syi**2*syy*usss
uzz      = rzi**2*urr+2.*rzi*szi*urs+2.*rzi*tzi*urt+szi**2*uss+2.*szi*tzi*ust+tzi**2*utt+rzz*ur+szz*us+tzz*ut
uxxzz    = rxi**2*rzi**2*urrrr+sxi**2*szi**2*ussss+txi**2*tzi**2*utttt+((2.*rxx*tzi+4.*rxz*txi)*szi+(2.*sxi*tzz+4.*sxz*tzi+4.*szi*txz+2.*szz*txi)*rxi+(4.*sxi*txz+2.*sxx*tzi+4.*sxz*txi+2.*szi*txx)*rzi+(4.*rxz*tzi+2.*rzz*txi)*sxi)*urst+((2.*sxx*tzi+4.*sxz*txi)*szi+txx*szi**2+(4.*sxz*tzi+4.*szi*txz+2.*szz*txi)*sxi+sxi**2*tzz)*usst+(rzz*txi**2+(2.*txi*tzz+4.*txz*tzi)*rxi+(4.*txi*txz+2.*txx*tzi)*rzi+rxx*tzi**2+4.*rxz*tzi*txi)*urtt+((4.*txi*txz+2.*txx*tzi)*szi+szz*txi**2+(2.*txi*tzz+4.*txz*tzi)*sxi+sxx*tzi**2+4.*sxz*tzi*txi)*ustt+(txi**2*tzz+4.*txi*txz*tzi+txx*tzi**2)*uttt+(2.*rxi*rxzz+rxx*rzz+2.*rxxz*rzi+2.*rxz**2)*urr+(2.*rxi*sxzz+rxx*szz+2.*rxxz*szi+4.*rxz*sxz+2.*rxzz*sxi+2.*rzi*sxxz+rzz*sxx)*urs+(2.*sxi*sxzz+sxx*szz+2.*sxxz*szi+2.*sxz**2)*uss+(2.*rxi*txzz+rxx*tzz+2.*rxxz*tzi+4.*rxz*txz+2.*rxzz*txi+2.*rzi*txxz+rzz*txx)*urt+(2.*sxi*txzz+sxx*tzz+2.*sxxz*tzi+4.*sxz*txz+2.*sxzz*txi+2.*szi*txxz+szz*txx)*ust+(2.*txi*txzz+txx*tzz+2.*txxz*tzi+2.*txz**2)*utt+rxxzz*ur+sxxzz*us+txxzz*ut+(2.*rxi**2*rzi*szi+2.*rxi*rzi**2*sxi)*urrrs+(rxi**2*szi**2+4.*rxi*rzi*sxi*szi+rzi**2*sxi**2)*urrss+(2.*rxi*sxi*szi**2+2.*rzi*sxi**2*szi)*ursss+(2.*rxi**2*rzi*tzi+2.*rxi*rzi**2*txi)*urrrt+(2.*sxi*txi*rzi**2+2.*szi*tzi*rxi**2+(4.*sxi*tzi+4.*szi*txi)*rzi*rxi)*urrst+((4.*sxi*szi*tzi+2.*szi**2*txi)*rxi+(2.*sxi**2*tzi+4.*sxi*szi*txi)*rzi)*ursst+(2.*sxi**2*szi*tzi+2.*sxi*szi**2*txi)*ussst+(rxi**2*tzi**2+4.*rxi*rzi*txi*tzi+rzi**2*txi**2)*urrtt+((4.*sxi*txi*tzi+2.*szi*txi**2)*rzi+(2.*sxi*tzi**2+4.*szi*txi*tzi)*rxi)*urstt+(sxi**2*tzi**2+4.*sxi*szi*txi*tzi+szi**2*txi**2)*usstt+(2.*rxi*txi*tzi**2+2.*rzi*txi**2*tzi)*urttt+(2.*sxi*txi*tzi**2+2.*szi*txi**2*tzi)*usttt+(rxi**2*rzz+4.*rxi*rxz*rzi+rxx*rzi**2)*urrr+(sxx*rzi**2+szz*rxi**2+(4.*rxz*szi+4.*rzi*sxz+2.*rzz*sxi)*rxi+(2.*rxx*szi+4.*rxz*sxi)*rzi)*urrs+(4.*rxz*szi*sxi+(2.*sxi*szz+4.*sxz*szi)*rxi+(4.*sxi*sxz+2.*sxx*szi)*rzi+rxx*szi**2+rzz*sxi**2)*urss+(sxi**2*szz+4.*sxi*sxz*szi+sxx*szi**2)*usss+(rzi**2*txx+(4.*rxz*tzi+4.*rzi*txz+2.*rzz*txi)*rxi+(2.*rxx*tzi+4.*rxz*txi)*rzi+rxi**2*tzz)*urrt
uyyzz    = (2.*syi**2*szi*tzi+2.*syi*szi**2*tyi)*ussst+(ryi**2*tzi**2+4.*ryi*rzi*tyi*tzi+rzi**2*tyi**2)*urrtt+((4.*syi*tyi*tzi+2.*szi*tyi**2)*rzi+(2.*syi*tzi**2+4.*szi*tyi*tzi)*ryi)*urstt+(syi**2*tzi**2+4.*syi*szi*tyi*tzi+szi**2*tyi**2)*usstt+(2.*ryi*tyi*tzi**2+2.*rzi*tyi**2*tzi)*urttt+(2.*syi*tyi*tzi**2+2.*szi*tyi**2*tzi)*usttt+(ryi**2*rzz+4.*ryi*ryz*rzi+ryy*rzi**2)*urrr+((2.*ryy*szi+4.*ryz*syi)*rzi+syy*rzi**2+szz*ryi**2+(4.*ryz*szi+4.*rzi*syz+2.*rzz*syi)*ryi)*urrs+(rzz*syi**2+ryy*szi**2+(4.*syi*syz+2.*syy*szi)*rzi+(2.*syi*szz+4.*syz*szi)*ryi+4.*ryz*szi*syi)*urss+(syi**2*szz+4.*syi*syz*szi+syy*szi**2)*usss+((2.*ryy*tzi+4.*ryz*tyi)*rzi+(4.*ryz*tzi+4.*rzi*tyz+2.*rzz*tyi)*ryi+ryi**2*tzz+rzi**2*tyy)*urrt+((4.*syi*tyz+2.*syy*tzi+4.*syz*tyi+2.*szi*tyy)*rzi+(2.*syi*tzz+4.*syz*tzi+4.*szi*tyz+2.*szz*tyi)*ryi+(4.*ryz*tzi+2.*rzz*tyi)*syi+(2.*ryy*tzi+4.*ryz*tyi)*szi)*urst+((4.*syz*tzi+4.*szi*tyz+2.*szz*tyi)*syi+syi**2*tzz+(2.*syy*tzi+4.*syz*tyi)*szi+tyy*szi**2)*usst+((4.*tyi*tyz+2.*tyy*tzi)*rzi+(2.*tyi*tzz+4.*tyz*tzi)*ryi+ryy*tzi**2+4.*ryz*tzi*tyi+rzz*tyi**2)*urtt+((2.*tyi*tzz+4.*tyz*tzi)*syi+syy*tzi**2+4.*syz*tzi*tyi+(4.*tyi*tyz+2.*tyy*tzi)*szi+szz*tyi**2)*ustt+(tyi**2*tzz+4.*tyi*tyz*tzi+tyy*tzi**2)*uttt+(2.*ryi*ryzz+ryy*rzz+2.*ryyz*rzi+2.*ryz**2)*urr+(2.*ryi*syzz+ryy*szz+2.*ryyz*szi+4.*ryz*syz+2.*ryzz*syi+2.*rzi*syyz+rzz*syy)*urs+(2.*syi*syzz+syy*szz+2.*syyz*szi+2.*syz**2)*uss+(2.*ryi*tyzz+ryy*tzz+2.*ryyz*tzi+4.*ryz*tyz+2.*ryzz*tyi+2.*rzi*tyyz+rzz*tyy)*urt+(2.*syi*tyzz+syy*tzz+2.*syyz*tzi+4.*syz*tyz+2.*syzz*tyi+2.*szi*tyyz+szz*tyy)*ust+(2.*tyi*tyzz+tyy*tzz+2.*tyyz*tzi+2.*tyz**2)*utt+ryyzz*ur+syyzz*us+tyyzz*ut+(2.*ryi**2*rzi*szi+2.*ryi*rzi**2*syi)*urrrs+(ryi**2*szi**2+4.*ryi*rzi*syi*szi+rzi**2*syi**2)*urrss+(2.*ryi*syi*szi**2+2.*rzi*syi**2*szi)*ursss+(2.*ryi**2*rzi*tzi+2.*ryi*rzi**2*tyi)*urrrt+((4.*syi*tzi+4.*szi*tyi)*rzi*ryi+2.*szi*tzi*ryi**2+2.*syi*tyi*rzi**2)*urrst+((2.*syi**2*tzi+4.*syi*szi*tyi)*rzi+(4.*syi*szi*tzi+2.*szi**2*tyi)*ryi)*ursst+syi**2*szi**2*ussss+tyi**2*tzi**2*utttt+rzi**2*ryi**2*urrrr
uzzzz    = rzi**4*urrrr+szi**4*ussss+tzi**4*utttt+(6.*rzi**2*szz+12.*rzi*rzz*szi)*urrs+(12.*rzi*szi*szz+6.*rzz*szi**2)*urss+(6.*rzi**2*tzz+12.*rzi*rzz*tzi)*urrt+((12.*szi*tzz+12.*szz*tzi)*rzi+12.*rzz*tzi*szi)*urst+(6.*szi**2*tzz+12.*szi*szz*tzi)*usst+(12.*rzi*tzi*tzz+6.*rzz*tzi**2)*urtt+(12.*szi*tzi*tzz+6.*szz*tzi**2)*ustt+(4.*rzi*rzzz+3.*rzz**2)*urr+(4.*rzi*szzz+6.*rzz*szz+4.*rzzz*szi)*urs+(4.*szi*szzz+3.*szz**2)*uss+(4.*rzi*tzzz+6.*rzz*tzz+4.*rzzz*tzi)*urt+(4.*szi*tzzz+6.*szz*tzz+4.*szzz*tzi)*ust+(4.*tzi*tzzz+3.*tzz**2)*utt+rzzzz*ur+szzzz*us+tzzzz*ut+12.*rzi**2*szi*tzi*urrst+12.*rzi*szi**2*tzi*ursst+12.*rzi*szi*tzi**2*urstt+4.*rzi**3*szi*urrrs+6.*rzi**2*szi**2*urrss+4.*rzi*szi**3*ursss+4.*rzi**3*tzi*urrrt+4.*szi**3*tzi*ussst+6.*rzi**2*tzi**2*urrtt+6.*szi**2*tzi**2*usstt+4.*rzi*tzi**3*urttt+4.*szi*tzi**3*usttt+6.*rzi**2*rzz*urrr+6.*szi**2*szz*usss+6.*tzi**2*tzz*uttt
! ---------- END CURVILINEAR  ---------
#End
#endMacro
