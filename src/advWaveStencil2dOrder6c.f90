! This file automatically generated from advWaveStencil.bf90 with bpp.
  subroutine advWaveStencil2dOrder6c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,sc,v,vh,lapCoeff,etax,etay,etaz,bc,frequencyArray,ipar,rpar,ierr )
 !======================================================================
 !   Advance a time step for Waves equations
 !
 ! nd : number of space dimensions
 ! um,u,un : u(t-dt), u(t), u(t+dt)
 !
 ! ipar(0)  = option : option=0 - advance wave equation
 !                           =1 - add upwind dissipation (predictor corrector mode)
 !
 !======================================================================
  implicit none
  integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b
   real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
   real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
   real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
   real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
   ! real stencilCoeff(0:*)   ! holds stencil coeffs
   real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1) 
   real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
   real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
   real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds coeff of Laplacian for HA scheme
   real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
   real etax(nd1a:nd1b)  ! superGrid functions
   real etay(nd2a:nd2b)
   real etaz(nd3a:nd3b)
   integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
   integer bc(0:1,0:2),ierr  
   real frequencyArray(0:*)
   integer ipar(0:*)
   real rpar(0:*)
  ! integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b
  ! real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  ! real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  ! real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  ! real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  ! ! real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
  ! real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1) 
  ! real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
  ! real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
  ! real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds coeff of Laplacian for HA scheme
  ! real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
  ! integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
  ! integer bc(0:1,0:2),ierr
  ! real frequencyArray(0:*)
  ! integer ipar(0:*)
  ! real rpar(0:*)
     real sc(1:49,nd1a:nd1b,nd2a:nd2b)
     real scr(1:49)
  !     ---- local variables -----
  integer gridIndexRange(0:1,0:2)
  integer m1a,m1b,m2a,m2b,m3a,m3b,numGhost,nStart,nEnd,mt,ig,useMask
  integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime,axis,dir,grid,freq
  integer addForcing,orderOfDissipation,option,gridIsImplicit,preComputeUpwindUt
  integer useNewForcingMethod,numberOfForcingFunctions,fcur,fnext,fprev,numberOfFrequencies
  real t,tm,cc,dt,dy,dz,cdt,cdtdx,cdtdy,cdtdz
  ! ,adc,adcdt,add,adddt
  real dt4by12
  ! logical addDissipation
  integer debug
  integer adjustHelmholtzForUpwinding
  real dx(0:2),dr(0:2)
  real c0,c1,csq,dtsq,cdtsq,cdtsq12,cdtSqBy12
  integer maxOrderOfAccuracy
  parameter( maxOrderOfAccuracy=12 )
  ! Coefficients in the implicit scheme
  real bImp(0:maxOrderOfAccuracy-1)
  real cImp(-1:1,0:maxOrderOfAccuracy-1)
  real alpha2,alpha4,alpha6,alpha8, beta2,beta4,beta6,beta8
  integer rectangular,curvilinear
  parameter( rectangular=0, curvilinear=1 )
  integer timeSteppingMethod
  integer defaultTimeStepping,adamsSymmetricOrder3,rungeKuttaFourthOrder,stoermerTimeStepping,modifiedEquationTimeStepping
  parameter(defaultTimeStepping=0,adamsSymmetricOrder3=1,rungeKuttaFourthOrder=2,stoermerTimeStepping=3,modifiedEquationTimeStepping=4)
 !...........start statement function
  integer kd,m
  ! real rx,ry,rz,sx,sy,sz,tx,ty,tz
  real cdtPow2,cdtPow4By12,cdtPow6By360,cdtPow8By20160
  real ff
  ! real cdSosupx,cdSosupy,cdSosupz
  real adSosup,sosupParameter, uDotFactor, adxSosup(0:2)
  integer useSosupDissipation,sosupDissipationOption
  integer updateSolution,updateDissipation,computeUt
  integer ec 
  real ep 
  real fv(0:1) , ev(0:1), evtt(0:1), evxx(0:1), evyy(0:1), evzz(0:1)
  real evxxxx(0:1), evxxyy(0:1), evyyyy(0:1), evxxzz(0:1), evyyzz(0:1), evzzzz(0:1), evtttt(0:1)
  real evtttttt(0:1)
  real evxxxxxx(0:1)
  real evyyyyyy(0:1)
  real evzzzzzz(0:1)       
  real evxxyyyy(0:1)
  real evxxxxyy(0:1)
  real evxxxxzz(0:1)
  real evxxzzzz(0:1)
  real evyyyyzz(0:1)
  real evyyzzzz(0:1)
  real evxxyyzz(0:1)
  real omega, coswt
  integer maxFreq
  parameter( maxFreq=500 )
  real cosFreqt(0:maxFreq), coswtAve(0:maxFreq), cosineFactor(0:maxFreq)
  integer idv(0:2),j1,j2,j3
  integer iStencil,upwCase,upwindHalfStencilWidth,i1l,i2l,i3l, i1r,i2r,i3r
  integer useUpwindDissipation,useImplicitUpwindDissipation,adjustOmega,solveHelmholtz
  real upw,maxDiff,umj
  ! real upwindCoeff(-3:3,0:3) 
  integer forcingOption
  ! forcingOptions -- these should match ForcingEnum in CgWave.h 
  ! enum ForcingOptionEnum
  ! {
  !   noForcing=0,
  !   twilightZoneForcing,
  !   userForcing,
  !   helmholtzForcing
  ! };
  integer noForcing,twilightZoneForcing,userForcing,helmholtzForcing
  parameter(noForcing           =0,twilightZoneForcing =1,userForcing         =2,helmholtzForcing    =3 )
   real maxErr(1:30), l2Err(1:30)
   real maxSol(30)
   real ue
   real uet8 
   real uex8 
   real uey8 
   real uez8 
   real uex6y2
   real uex4y4
   real uex2y6
   real uex6z2
   real uex4z4
   real uex2z6
   real uey6z2
   real uey4z4
   real uey2z6
   real uex4y2z2
   real uex2y4z2
   real uex2y2z4
   integer maxDeriv,d,uc,count,numGhost1,m1,m2,m3
   ! ====== variables for stencils ========
   ! real czm,czp,czz,cmz,cpz
   real dt2,dt4,dt6,dt8
   real dx2i,dy2i,dz2i
   real cdx2i,cdy2i,cdz2i
   integer stencilWidth,numStencilCoeff
   real dr1, dr2, dr3, dr1i, dr2i, dr3i, rx, ry, rz, sx, sy
   real sz, tx, ty, tz, diffOrder1, diffOrder2, diffOrder3, rxr, rxs, rxt, ryr
   real rys, ryt, rzr, rzs, rzt, sxr, sxs, sxt, syr, sys, syt
   real szr, szs, szt, txr, txs, txt, tyr, tys, tyt, tzr, tzs
   real tzt, rxx, ryy, rzz, sxx, syy, szz, txx, tyy, tzz
   ! *** The next include files were generated by cgWave/maple/writeStencilFiles.mpl ***
! Define variables to valuate stencil coefficients, dim=2, order=6, gridType=Curvilinear
! File generated by cgWave/maple/writeStencilFiles.mpl
integer i1m3,i1m2,i1m1,i1p1,i1p2,i1p3
integer i2m3,i2m2,i2m1,i2p1,i2p2,i2p3
real t0,t1,t2,t4,t5,t8,t9,t11,t12,t14,t15,t16,t17,t18,t24,t26,t29,t36,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50,t51,t52,t53,t55,t58,t59,t63,t66,t67,t71,t83,t84,t85,t86,t90,t91,t92,t94,t99,t100,t101,t102,t103,t104,t107,t108,t111,t117,t118,t121,t123,t124,t126,t127,t128,t129,t136,t139,t140,t145,t152,t153,t154,t155,t157,t165,t166,t170,t174,t183,t184,t185,t186,t193,t194,t195,t196,t202,t205,t215,t219,t220,t222,t223,t225,t226,t227,t228,t229,t235,t237,t241,t251,t252,t253,t254,t255,t256,t257,t258,t259,t260,t261,t262,t263,t264,t265,t271,t272,t273,t274,t275,t276,t277,t278,t279,t280,t282,t285,t286,t287,t288,t289,t290,t291,t292,t294,t298,t299,t303,t304,t307,t308,t311,t314,t315,t318,t319,t329,t331,t332,t333,t334,t335,t336,t337,t338,t339,t340,t342,t343,t344,t347,t348,t349,t352,t358,t360,t361,t363,t365,t370,t371,t372,t373,t378,t382,t383,t387,t392,t394,t395,t397,t398,t399,t400,t402,t403,t411,t412,t413,t414,t424,t425,t435,t436,t437,t438,t440,t442,t445,t446,t447,t448,t450,t451,t454,t455,t456,t458,t459,t460,t462,t463,t464,t467,t468,t469,t470,t471,t472,t473,t474,t475,t476,t478,t488,t491,t498,t499,t506,t507,t513,t514,t515,t516,t519,t520,t521,t522,t525,t526,t527,t528,t531,t532,t533,t548,t550,t551,t552,t553,t560,t562,t564,t565,t566,t568,t569,t570,t573,t574,t577,t580,t581,t582,t584,t585,t586,t589,t596,t598,t600,t605,t606,t607,t608,t613,t617,t629,t630,t631,t632,t633,t634,t635,t636,t642,t643,t645,t649,t650,t651,t652,t654,t658,t663,t664,t673,t676,t677,t685,t687,t705,t706,t707,t708,t709,t710,t711,t712,t713,t714,t715,t717,t720,t721,t725,t731,t732,t736,t745,t746,t747,t748,t752,t753,t755,t760,t761,t762,t763,t764,t765,t766,t767,t768,t769,t770,t771,t772,t775,t776,t778,t781,t787,t789,t790,t792,t793,t794,t795,t796,t802,t810,t811,t812,t816,t821,t822,t824,t826,t827,t828,t830,t831,t832,t840,t841,t851,t852,t860,t861,t862,t863,t866,t867,t868,t871,t873,t874,t876,t879,t880,t881,t884,t885,t886,t892,t893,t894,t895,t896,t897,t898,t899,t901,t902,t904,t905,t906,t909,t910,t913,t914,t915,t916,t918,t920,t921,t922,t923,t924,t925,t927,t928,t930,t931,t932,t934,t935,t936,t938,t939,t940,t941,t942,t943,t944,t946,t947,t948,t950,t951,t952,t954,t955,t957,t958,t959,t961,t962,t987,t988,t995,t997,t998,t1002,t1007,t1008,t1011,t1021,t1022,t1023,t1026,t1028,t1032,t1033,t1034,t1035,t1038,t1039,t1041,t1048,t1049,t1052,t1054,t1055,t1057,t1058,t1059,t1060,t1065,t1069,t1074,t1076,t1080,t1081,t1088,t1093,t1095,t1096,t1098,t1099,t1102,t1104,t1106,t1107,t1108,t1111,t1112,t1113,t1114,t1116,t1118,t1119,t1122,t1124,t1125,t1126,t1128,t1129,t1131,t1132,t1137,t1139,t1141,t1143,t1144,t1145,t1147,t1149,t1150,t1151,t1152,t1160,t1162,t1163,t1164,t1166,t1167,t1172,t1173,t1174,t1175,t1176,t1179,t1180,t1188,t1190,t1191,t1198,t1199,t1206,t1209,t1215,t1218,t1221,t1222,t1226,t1230,t1232,t1233,t1235,t1238,t1241,t1244,t1246,t1247,t1253,t1256,t1258,t1260,t1262,t1263,t1265,t1266,t1267,t1268,t1269,t1271,t1273,t1274,t1276,t1278,t1279,t1281,t1283,t1285,t1287,t1291,t1292,t1293,t1299,t1300,t1301,t1303,t1305,t1306,t1307,t1308,t1309,t1311,t1312,t1313,t1314,t1315,t1318,t1323,t1324,t1343,t1344,t1346,t1347,t1348,t1351,t1360,t1361,t1381,t1383,t1387,t1388,t1389,t1390,t1392,t1396,t1397,t1417,t1420,t1421,t1422,t1423,t1424,t1425,t1428,t1429,t1435,t1436,t1439,t1441,t1442,t1444,t1445,t1446,t1447,t1456,t1457,t1460,t1463,t1472,t1473,t1474,t1475,t1477,t1480,t1482,t1483,t1484,t1486,t1487,t1488,t1489,t1491,t1492,t1495,t1496,t1498,t1499,t1500,t1503,t1504,t1505,t1506,t1507,t1508,t1509,t1510,t1511,t1512,t1514,t1524,t1525,t1532,t1533,t1542,t1543,t1549,t1550,t1553,t1554,t1555,t1556,t1559,t1560,t1561,t1562,t1563,t1566,t1567,t1574,t1575,t1576,t1579,t1580,t1581,t1585,t1587,t1590,t1592,t1594,t1600,t1601,t1604,t1606,t1607,t1609,t1610,t1614,t1615,t1616,t1621,t1623,t1624,t1635,t1636,t1640,t1645,t1648,t1649,t1650,t1653,t1655,t1656,t1657,t1659,t1660,t1663,t1665,t1666,t1667,t1668,t1669,t1671,t1674,t1676,t1677,t1678,t1679,t1681,t1683,t1684,t1689,t1691,t1693,t1695,t1697,t1699,t1700,t1701,t1702,t1703,t1704,t1712,t1714,t1715,t1718,t1719,t1722,t1723,t1724,t1725,t1728,t1729,t1735,t1737,t1738,t1739,t1740,t1741,t1743,t1744,t1745,t1746,t1748,t1750,t1751,t1753,t1754,t1755,t1758,t1759,t1760,t1761,t1762,t1763,t1765,t1766,t1767,t1769,t1771,t1773,t1774,t1775,t1778,t1779,t1782,t1784,t1785,t1786,t1789,t1790,t1791,t1792,t1795,t1796,t1799,t1800,t1802,t1816,t1821,t1823,t1827,t1830,t1837,t1847,t1848,t1849,t1850,t1854,t1855,t1857,t1858,t1869,t1871,t1872,t1876,t1877,t1878,t1883,t1885,t1886,t1897,t1901,t1906,t1909,t1911,t1912,t1913,t1915,t1917,t1918,t1919,t1920,t1921,t1922,t1923,t1925,t1926,t1928,t1932,t1934,t1936,t1938,t1940,t1942,t1943,t1944,t1945,t1946,t1949,t1950,t1953,t1956,t1959,t1961,t1962,t1964,t1972,t1973,t1976,t1977,t1978,t1979,t1980,t1988,t1989,t1990,t1991,t1992,t1993,t1995,t1999,t2000,t2002,t2005,t2007,t2008,t2010,t2012,t2016,t2017,t2018,t2023,t2032,t2033,t2044,t2045,t2071,t2074,t2076,t2077,t2089,t2090,t2091,t2092,t2094,t2102,t2106,t2107,t2111,t2120,t2121,t2122,t2123,t2139,t2140,t2141,t2143,t2150,t2151,t2153,t2155,t2156,t2157,t2158,t2159,t2162,t2163,t2166,t2169,t2170,t2171,t2173,t2182,t2184,t2185,t2186,t2187,t2188,t2194,t2202,t2203,t2216,t2217,t2218,t2219,t2220,t2222,t2225,t2226,t2227,t2228,t2229,t2236,t2237,t2244,t2245,t2252,t2253,t2260,t2265,t2268,t2271,t2275,t2276,t2277,t2284,t2286,t2287,t2289,t2292,t2295,t2298,t2300,t2301,t2307,t2310,t2312,t2314,t2316,t2317,t2319,t2320,t2321,t2322,t2323,t2324,t2326,t2328,t2330,t2333,t2335,t2336,t2337,t2339,t2341,t2345,t2346,t2347,t2357,t2358,t2360,t2361,t2365,t2366,t2367,t2368,t2379,t2381,t2382,t2383,t2384,t2389,t2393,t2398,t2400,t2404,t2411,t2416,t2418,t2419,t2421,t2423,t2425,t2427,t2428,t2429,t2430,t2431,t2432,t2433,t2435,t2436,t2438,t2442,t2444,t2446,t2448,t2449,t2450,t2452,t2454,t2455,t2456,t2459,t2460,t2463,t2466,t2469,t2471,t2472,t2474,t2484,t2485,t2487,t2488,t2489,t2490,t2492,t2493,t2503,t2504,t2506,t2509,t2512,t2515,t2517,t2518,t2522,t2523,t2525,t2528,t2531,t2534,t2536,t2542,t2545,t2546,t2548,t2549,t2551,t2553,t2554,t2556,t2559,t2560,t2564,t2568,t2569,t2587,t2588,t2589,t2593,t2597,t2603,t2604,t2605,t2606,t2621,t2627,t2629,t2631,t2632,t2633,t2634,t2639,t2643,t2644,t2645,t2673,t2701,t2702,t2703,t2704,t2723,t2724,t2725,t2726,t2727,t2728,t2729,t2730,t2736,t2737,t2738,t2739,t2741,t2745,t2746,t2748,t2754,t2759,t2760,t2767,t2768,t2769,t2779,t2780,t2781,t2782,t2784,t2785,t2786,t2787,t2788,t2791,t2792,t2793,t2795,t2796,t2797,t2802,t2803,t2823,t2824,t2825,t2836,t2837,t2854,t2855,t2862,t2863,t2864,t2865,t2871,t2872,t2873,t2874,t2875,t2876,t2878,t2882,t2883,t2885,t2888,t2889,t2891,t2893,t2897,t2899,t2900,t2901,t2910,t2917,t2918,t2927,t2928,t2937,t2938,t2939,t2940,t2955,t2961,t2963,t2964,t2965,t2967,t2968,t2973,t2977,t2978,t2979,t3007,t3012,t3028,t3063,t3065,t3087,t3089,t3093,t3094,t3095,t3096,t3098,t3102,t3103,t3140,t3143,t3145,t3146,t3159;
    ! #Include "../include/defineStencilVariables2dOrder2Curvilinear.h"
   !...........end   statement functions
   ! write(*,*) 'Inside advWaveStencil...'
   cc             = rpar( 0)  ! this is c
   dt             = rpar( 1)
   dx(0)          = rpar( 2)
   dx(1)          = rpar( 3)
   dx(2)          = rpar( 4)
   dr(0)          = rpar( 5)
   dr(1)          = rpar( 6)
   dr(2)          = rpar( 7)
   t              = rpar( 8)
   ep             = rpar( 9)
   sosupParameter = rpar(10)
   omega          = rpar(11) ! for helmholtz 
   bImp( 0)       = rpar(12) ! beta2 : coefficient for implicit time-stepping
   bImp( 1)       = rpar(13) ! beta4 : coefficient for implicit time-stepping
   bImp( 2)       = rpar(14) ! beta6 (for future)
   bImp( 3)       = rpar(15) ! beta8 (for future)
   dy=dx(1)  ! Are these needed?
   dz=dx(2)
   option                       = ipar( 0)
   grid                         = ipar( 1)
   gridType                     = ipar( 2)
   orderOfAccuracy              = ipar( 3)
   orderInTime                  = ipar( 4)
   addForcing                   = ipar( 5)
   forcingOption                = ipar( 6)
   numberOfForcingFunctions     = ipar( 7)
   fcur                         = ipar( 8) 
   debug                        = ipar( 9)
   gridIsImplicit               = ipar(10)
   useUpwindDissipation         = ipar(11)  ! explicit upwind dissipation
   useImplicitUpwindDissipation = ipar(12)  ! true if upwind-dissipation is on for impliciit time-stepping
   preComputeUpwindUt           = ipar(13)
   numberOfFrequencies          = ipar(14)
   adjustOmega                  = ipar(15)
   solveHelmholtz               = ipar(16)
   adjustHelmholtzForUpwinding  = ipar(17)
   fprev = mod(fcur-1+numberOfForcingFunctions,max(1,numberOfForcingFunctions))
   fnext = mod(fcur+1                         ,max(1,numberOfForcingFunctions))
   ! ** fix me ***
   timeSteppingMethod=modifiedEquationTimeStepping
   useMask=0  ! do this for now -- do not check mask in loops, these seems faster
   ! Set dr(:) = dx(:) for 6th-order derivatives
   if( gridType.eq.rectangular )then
     do axis=0,2
       dr(axis)=dx(axis)
     end do
   else
     do axis=0,2
       dx(axis)=dr(axis)
     end do
   end if  
   ! Do this for now: 
   maxDeriv=6
   uc=0
   gridIndexRange(0,0)=n1a
   gridIndexRange(1,0)=n1b
   gridIndexRange(0,1)=n2a
   gridIndexRange(1,1)=n2b
   gridIndexRange(0,2)=n3a
   gridIndexRange(1,2)=n3b    
   ! ---- Compute the coefficients in the implicit time-stepping scheme ----
   beta2=bImp(0)
   beta4=bImp(1)
   alpha2 = (1.-beta2)/2.
   alpha4 = (alpha2-beta4-1./12.)/2. 
   cImp(-1,0)=alpha2
   cImp( 0,0)= beta2
   cImp( 1,0)=alpha2
   cImp(-1,1)=alpha4
   cImp( 0,1)= beta4
   cImp( 1,1)=alpha4  
   csq=cc**2
   dtsq=dt**2
   cdtsq=(cc**2)*(dt**2)
   cdt=cc*dt
   ! new: 
   cdtPow2        = cdt**2
   cdtPow4By12    = cdt**4/12.
   cdtPow6By360   = cdt**6/360. 
   cdtPow8By20160 = cdt**8/20160.  
   cdtsq12=cdtsq*cdtsq/12.  ! c^4 dt^4 /12 
   ! cdt4by360=(cdt)**4/360.  ! (c*dt)^4/360 
   ! cdt6by20160=cdt**6/(8.*7.*6.*5.*4.*3.)
   cdtSqBy12= cdtsq/12.   ! c^2*dt*2/12
   dt4by12=dtsq*dtsq/12.
   cdtdx = (cc*dt/dx(0))**2
   cdtdy = (cc*dt/dy)**2
   cdtdz = (cc*dt/dz)**2
   dt2 = dt*dt;
   dt4 = dt2*dt2;
   dt6 = dt4*dt2;
   dt8 = dt6*dt2;
   dx2i = 1./dx(0)**2
   dy2i = 1./dx(1)**2
   dz2i = 1./dx(2)**2
   cdx2i = csq/dx(0)**2
   cdy2i = csq/dx(1)**2
   cdz2i = csq/dx(2)**2  
   ! if( option.eq.1 )then 
   !  useSosupDissipation = 1
   ! else
   !  useSosupDissipation = 0
   ! end if
   if( (.true. .or. debug.gt.1) .and. t.le.dt )then
     write(*,'("advWaveStencil: option=",i4," grid=",i4)') option,grid
     write(*,'("advWaveStencil: orderOfAccuracy=",i2," orderInTime=",i2  )') orderOfAccuracy,orderInTime
     write(*,'("advWaveStencil: addForcing=",i2," forcingOption=",i2)') addForcing,forcingOption
     ! write(*,'("advWaveStencil: useUpwindDissipation=",i2,"(explicit), useImplicitUpwindDissipation=",i2," (implicit)")') useUpwindDissipation,useImplicitUpwindDissipation
     write(*,'("advWaveStencil: t,dt,c,omega=",4e10.2)') t,dt,cc,omega 
     write(*,'("advWaveStencil: gridIsImplicit=",i2," adjustOmega=",i2," solveHelmholtz=",i2)') gridIsImplicit,adjustOmega,solveHelmholtz
     if( forcingOption.eq.helmholtzForcing )then
       write(*,'("advWaveStencil: numberOfFrequencies=",i2)') numberOfFrequencies
       write(*,'("advWaveStencil: frequencyArray=",(1pe12.4,1x))') (frequencyArray(freq),freq=0,numberOfFrequencies-1)
     end if
     if( gridIsImplicit.eq.1 )then
       write(*,'("  Implicit coeff: cImp(-1:1,0) = ",3(1pe10.2,1x), "(for 2nd-order)")') cImp(-1,0),cImp(0,0),cImp(1,0)
       write(*,'("  Implicit coeff: cImp(-1:1,1) = ",3(1pe10.2,1x), "(for 4th-order)")') cImp(-1,1),cImp(0,1),cImp(1,1)
     end if
   end if
   if( forcingOption.eq.helmholtzForcing )then
     ! --- solving the Helmholtz problem ---
     if( t.le.dt .and. debug.gt.1 )then
       write(*,'("advWaveStencil: numberOfFrequencies=",i6," omega=",1pe12.4," frequencyArray(0)=",1pe12.4)') numberOfFrequencies,omega,frequencyArray(0)
     end if
     if( numberOfFrequencies.le.0 )then
       write(*,'("advWaveStencil: ERROR: numberOfFrequencies=",i6," is <= 0")') numberOfFrequencies
       stop 0123
     end if
     if( numberOfFrequencies.eq.1  .and. frequencyArray(0) .ne. omega )then
       write(*,'("advWaveStencil: ERROR: frequencyArray(0)=",1pe12.4," is not equal to omega=",1pe12.4)') frequencyArray(0),omega
       stop 1234
     end if
     if( numberOfFrequencies.gt.maxFreq )then
       write(*,'("advWaveStencil: ERROR: numberOfFrequencies > maxFreq=",i6," .. FIX ME")') maxFreq
       stop 2345
     end if
     ! if( numberOfFrequencies.gt.1 .and. gridIsImplicit.eq.1 )then
     !   write(*,'("advWave: ERROR: numberOfFrequencies > 1 and implicit time-stepping : FINISH ME")') 
     !   stop 3456  
     ! end if
     do freq=0,numberOfFrequencies-1
       cosFreqt(freq) = cos(frequencyArray(freq)*t)
     end do
   end if
   ! write(*,'(" advWave: timeSteppingMethod=",i2)') timeSteppingMethod
   if( timeSteppingMethod.eq.defaultTimeStepping )then
    write(*,'(" advWaveStencil:ERROR: timeSteppingMethod=defaultTimeStepping -- this should be set")')
      ! '
    stop 83322
   end if
   ! -- first time through, eval lapCoeff and the stencil coefficients --
     if( lapCoeff(0,0,0,0).lt.0. )then
          dr1=dr(0); dr1i=1./dr1;
          dr2=dr(1); dr2i=1./dr2;
          dr3=dr(2); dr3i=1./dr3;
          ! --- Evaluate and store coefficients in Laplacian ---
          write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF 2D'
          numGhost1=orderOfAccuracy/2 -1; ! check me 
          n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
          n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
          n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
           do i3=n3a,n3b
           do i2=n2a,n2b
           do i1=n1a,n1b
            rx = rsxy(i1,i2,i3,0,0)
            ry = rsxy(i1,i2,i3,0,1)
            sx = rsxy(i1,i2,i3,1,0)
            sy = rsxy(i1,i2,i3,1,1)
            ! --- choose order for (r,s,t) derivatives based on available ghost points, less accuracy is needed in ghost points  ---
            if( (i1-4).ge.nd1a .and. (i1+4).le.nd1b )then
              diffOrder1=min(8,6)
            elseif( (i1-3).ge.nd1a .and. (i1+3).le.nd1b )then
              diffOrder1=min(6,6)
            elseif( (i1-2).ge.nd1a .and. (i1+2).le.nd1b )then
              diffOrder1=min(4,6)
            elseif( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
              diffOrder1=min(2,6)
            else
              write(*,*) "i1,nd1a,nd1b=",i1,nd1a,nd1b
              stop 999
            end if
            if( (i2-4).ge.nd2a .and. (i2+4).le.nd2b )then
              diffOrder2=min(8,6)
            elseif( (i2-3).ge.nd2a .and. (i2+3).le.nd2b )then
              diffOrder2=min(6,6)
            elseif( (i2-2).ge.nd2a .and. (i2+2).le.nd2b )then
              diffOrder2=min(4,6)
            elseif( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
              diffOrder2=min(2,6)
            else
              write(*,*) "i2,nd2a,nd2b=",i2,nd2a,nd2b
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
            elseif( diffOrder1.eq.8 )then
              rxr = ( 672.*(rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0)) -168.*(rsxy(i1+2,i2,i3,0,0)-rsxy(i1-2,i2,i3,0,0)) +32*(rsxy(i1+3,i2,i3,0,0)-rsxy(i1-3,i2,i3,0,0)) -3.*(rsxy(i1+4,i2,i3,0,0)-rsxy(i1-4,i2,i3,0,0)) )*(dr1i/840.) 
              ryr = ( 672.*(rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1)) -168.*(rsxy(i1+2,i2,i3,0,1)-rsxy(i1-2,i2,i3,0,1)) +32*(rsxy(i1+3,i2,i3,0,1)-rsxy(i1-3,i2,i3,0,1)) -3.*(rsxy(i1+4,i2,i3,0,1)-rsxy(i1-4,i2,i3,0,1)) )*(dr1i/840.) 
              sxr = ( 672.*(rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0)) -168.*(rsxy(i1+2,i2,i3,1,0)-rsxy(i1-2,i2,i3,1,0)) +32*(rsxy(i1+3,i2,i3,1,0)-rsxy(i1-3,i2,i3,1,0)) -3.*(rsxy(i1+4,i2,i3,1,0)-rsxy(i1-4,i2,i3,1,0)) )*(dr1i/840.) 
              syr = ( 672.*(rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1)) -168.*(rsxy(i1+2,i2,i3,1,1)-rsxy(i1-2,i2,i3,1,1)) +32*(rsxy(i1+3,i2,i3,1,1)-rsxy(i1-3,i2,i3,1,1)) -3.*(rsxy(i1+4,i2,i3,1,1)-rsxy(i1-4,i2,i3,1,1)) )*(dr1i/840.) 
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
            elseif( diffOrder2.eq.8 )then
              rxs = ( 672.*(rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0)) -168.*(rsxy(i1,i2+2,i3,0,0)-rsxy(i1,i2-2,i3,0,0)) +32*(rsxy(i1,i2+3,i3,0,0)-rsxy(i1,i2-3,i3,0,0)) -3.*(rsxy(i1,i2+4,i3,0,0)-rsxy(i1,i2-4,i3,0,0)) )*(dr2i/840.) 
              rys = ( 672.*(rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1)) -168.*(rsxy(i1,i2+2,i3,0,1)-rsxy(i1,i2-2,i3,0,1)) +32*(rsxy(i1,i2+3,i3,0,1)-rsxy(i1,i2-3,i3,0,1)) -3.*(rsxy(i1,i2+4,i3,0,1)-rsxy(i1,i2-4,i3,0,1)) )*(dr2i/840.) 
              sxs = ( 672.*(rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0)) -168.*(rsxy(i1,i2+2,i3,1,0)-rsxy(i1,i2-2,i3,1,0)) +32*(rsxy(i1,i2+3,i3,1,0)-rsxy(i1,i2-3,i3,1,0)) -3.*(rsxy(i1,i2+4,i3,1,0)-rsxy(i1,i2-4,i3,1,0)) )*(dr2i/840.) 
              sys = ( 672.*(rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1)) -168.*(rsxy(i1,i2+2,i3,1,1)-rsxy(i1,i2-2,i3,1,1)) +32*(rsxy(i1,i2+3,i3,1,1)-rsxy(i1,i2-3,i3,1,1)) -3.*(rsxy(i1,i2+4,i3,1,1)-rsxy(i1,i2-4,i3,1,1)) )*(dr2i/840.) 
            end if
            rxx = rx*rxr + sx*rxs 
            ryy = ry*ryr + sy*rys 
            sxx = rx*sxr + sx*sxs 
            syy = ry*syr + sy*sys 
            ! -- Coefficients in the Laplacian (scaled)
            lapCoeff(i1,i2,i3,0) = (rx**2 + ry**2   )*dr1i**2
            lapCoeff(i1,i2,i3,2) = 2.*(rx*sx + ry*sy)*dr1i*dr2i
            lapCoeff(i1,i2,i3,1) = (sx**2 + sy**2   )*dr2i**2
            lapCoeff(i1,i2,i3,3) = (rxx + ryy       )*dr1i
            lapCoeff(i1,i2,i3,4) = (sxx + syy       )*dr2i 
           end do
           end do
           end do
          ! if( c200(0,0,0).le.0. )then
          !   ! --- Evaluate and store coefficients in Laplacian ---
          !   write(*,*) 'ASSIGN SCALED LAPLACIAN COEFF'
          !   numGhost1=orderOfAccuracy/2 -1; ! check me 
          !   n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
          !   n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
          !   n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
          !   ! ----- START LOOPS ----
          !   beginLoops3d()
          !   #If "noMask" eq "USEMASK" 
          !     if( mask(i1,i2,i3).ne.0 )then
          !   #End 
          !     rx = rsxy(i1,i2,i3,0,0)
          !     ry = rsxy(i1,i2,i3,0,1)
          !     sx = rsxy(i1,i2,i3,1,0)
          !     sy = rsxy(i1,i2,i3,1,1)
          !     ! --- choose order for (r,s,t) derivatives based on available ghost points, less accuracy is needed in ghost points  ---
          !     if( (i1-1).ge.nd1a .and. (i1+1).le.nd1b )then
          !       diffOrder1=2
          !     else
          !       stop 999
          !     end if
          !     if( (i2-1).ge.nd2a .and. (i2+1).le.nd2b )then
          !       diffOrder2=2
          !     else
          !       stop 999
          !     end if
          !     if( diffOrder1.eq.2 )then
          !       rxr = (rsxy(i1+1,i2,i3,0,0)-rsxy(i1-1,i2,i3,0,0))*(.5*dr1i) 
          !       ryr = (rsxy(i1+1,i2,i3,0,1)-rsxy(i1-1,i2,i3,0,1))*(.5*dr1i) 
          !       sxr = (rsxy(i1+1,i2,i3,1,0)-rsxy(i1-1,i2,i3,1,0))*(.5*dr1i) 
          !       syr = (rsxy(i1+1,i2,i3,1,1)-rsxy(i1-1,i2,i3,1,1))*(.5*dr1i) 
          !     end if
          !     if( diffOrder2.eq.2 )then
          !       rxs = (rsxy(i1,i2+1,i3,0,0)-rsxy(i1,i2-1,i3,0,0))*(.5*dr2i) 
          !       rys = (rsxy(i1,i2+1,i3,0,1)-rsxy(i1,i2-1,i3,0,1))*(.5*dr2i) 
          !       sxs = (rsxy(i1,i2+1,i3,1,0)-rsxy(i1,i2-1,i3,1,0))*(.5*dr2i) 
          !       sys = (rsxy(i1,i2+1,i3,1,1)-rsxy(i1,i2-1,i3,1,1))*(.5*dr2i) 
          !     end if
          !     rxx = rx*rxr + sx*rxs 
          !     ryy = ry*ryr + sy*rys 
          !     sxx = rx*sxr + sx*sxs 
          !     syy = ry*syr + sy*sys 
          !     ! -- Coefficients in the Laplacian (scaled)
          !     c200(i1,i2,i3) = (rx**2 + ry**2   )*dr1i**2
          !     c110(i1,i2,i3) = 2.*(rx*sx + ry*sy)*dr1i*dr2i*.25
          !     c020(i1,i2,i3) = (sx**2 + sy**2   )*dr2i**2
          !     c100(i1,i2,i3) = (rxx + ryy       )*dr1i*.5
          !     c010(i1,i2,i3) = (sxx + syy       )*dr2i*.5 
          !   #If "noMask" eq "USEMASK" 
          !     end if ! mask .ne. 0
          !   #End 
          !   endLoops3d() 
          !   ! ----- END LOOPS ----
          ! end if ! end assignLapCoeff
         write(*,*) 'EVAL STENCIL COEFF'
          ! ---- 2D -----
         numGhost1=0;
         n1a=max(nd1a,gridIndexRange(0,0)-numGhost1);  n1b=min(nd1b,gridIndexRange(1,0)+numGhost1);
         n2a=max(nd2a,gridIndexRange(0,1)-numGhost1);  n2b=min(nd2b,gridIndexRange(1,1)+numGhost1);
         n3a=max(nd3a,gridIndexRange(0,2)-numGhost1);  n3b=min(nd3b,gridIndexRange(1,2)+numGhost1);
         ! ----- START LOOPS ----
          do i3=n3a,n3b
          do i2=n2a,n2b
          do i1=n1a,n1b
! Evaluate stencil coefficients, dim=2, order=6, gridType=Curvilinear
! File generated by cgWave/maple/writeStencilFiles.mpl
i1m1=i1-1; i1p1=i1+1;
i2m1=i2-1; i2p1=i2+1;
i1m2=i1-2; i1p2=i1+2;
i2m2=i2-2; i2p2=i2+2;
t1 = lapCoeff(i1m2,i2m2,i3,2);
t2 = lapCoeff(i1m1,i2m1,i3,2);
t4 = lapCoeff(i1,i2,i3,2);
t5 = t4 * dt6;
t8 = lapCoeff(i1m1,i2m1,i3,4);
t9 = lapCoeff(i1m1,i2m2,i3,2);
t11 = t8 * t9 / 0.32e2;
t12 = lapCoeff(i1m1,i2m1,i3,1);
t14 = t12 * t9 / 0.16e2;
t15 = lapCoeff(i1m2,i2m2,i3,4);
t16 = t15 / 0.8e1;
t17 = lapCoeff(i1m2,i2m2,i3,1);
t18 = t17 / 0.4e1;
t24 = lapCoeff(i1,i2,i3,4);
t26 = lapCoeff(i1,i2m1,i3,2);
t29 = lapCoeff(i1,i2,i3,1);
t36 = t4 * dt4;
t41 = t4 * dt2;
t42 = t41 / 0.120e3;
t43 = lapCoeff(i1,i2m1,i3,4);
t44 = lapCoeff(i1,i2m2,i3,2);
t45 = t43 * t44;
t46 = t45 / 0.8e1;
t47 = lapCoeff(i1,i2m1,i3,1);
t48 = t47 * t44;
t49 = t48 / 0.4e1;
t50 = lapCoeff(i1m1,i2m2,i3,4);
t51 = t50 / 0.8e1;
t52 = lapCoeff(i1m1,i2m2,i3,1);
t53 = t52 / 0.4e1;
t55 = t26 * (t51 - t53);
t58 = t45 / 0.16e2;
t59 = t48 / 0.8e1;
t63 = lapCoeff(i1p1,i2m1,i3,2);
t66 = t50 / 0.4e1;
t67 = t52 / 0.2e1;
t71 = t50 / 0.2e1;
t83 = t50 / 0.48e2;
t84 = t8 / 0.48e2;
t85 = t52 / 0.24e2;
t86 = t12 / 0.48e2;
t90 = t24 * (t44 + t26) / 0.48e2;
t91 = t44 / 0.48e2;
t92 = t26 / 0.24e2;
t94 = t29 * (t91 + t92);
t99 = lapCoeff(i1,i2m2,i3,4);
t100 = t99 / 0.24e2;
t101 = t43 / 0.12e2;
t102 = lapCoeff(i1,i2m2,i3,1);
t103 = t102 / 0.12e2;
t104 = t47 / 0.12e2;
t107 = t43 / 0.24e2;
t108 = t47 / 0.24e2;
t111 = lapCoeff(i1p1,i2m2,i3,2);
t117 = t24 / 0.60e2;
t118 = t29 / 0.90e2;
t121 = lapCoeff(i1p1,i2m1,i3,4);
t123 = t121 * t111 / 0.32e2;
t124 = lapCoeff(i1p1,i2m1,i3,1);
t126 = t124 * t111 / 0.16e2;
t127 = t99 / 0.8e1;
t128 = t102 / 0.4e1;
t129 = t127 - t128;
t136 = t99 / 0.2e1;
t139 = t99 / 0.4e1;
t140 = t102 / 0.2e1;
t145 = t47 * (t136 - t102) - t43 * (t139 - t140) + t26 * (t9 + t111) / 0.16e2;
t152 = lapCoeff(i1p1,i2m2,i3,4);
t153 = t152 / 0.8e1;
t154 = lapCoeff(i1p1,i2m2,i3,1);
t155 = t154 / 0.4e1;
t157 = t26 * (t153 - t155);
t165 = t152 / 0.4e1;
t166 = t154 / 0.2e1;
t170 = t152 / 0.2e1;
t174 = lapCoeff(i1p2,i2m2,i3,2);
t183 = t152 / 0.48e2;
t184 = t121 / 0.48e2;
t185 = t154 / 0.24e2;
t186 = t124 / 0.48e2;
t193 = lapCoeff(i1p2,i2m2,i3,4);
t194 = t193 / 0.8e1;
t195 = lapCoeff(i1p2,i2m2,i3,1);
t196 = t195 / 0.4e1;
t202 = t24 * t26;
t205 = t29 * t26;
t215 = t4 * t63;
t219 = lapCoeff(i1m1,i2m1,i3,3);
t220 = lapCoeff(i1m2,i2m1,i3,2);
t222 = t219 * t220 / 0.32e2;
t223 = lapCoeff(i1m1,i2m1,i3,0);
t225 = t220 * t223 / 0.16e2;
t226 = lapCoeff(i1m2,i2m2,i3,3);
t227 = t226 / 0.8e1;
t228 = lapCoeff(i1m2,i2m2,i3,0);
t229 = t228 / 0.4e1;
t235 = lapCoeff(i1,i2,i3,3);
t237 = lapCoeff(i1m1,i2,i3,2);
t241 = lapCoeff(i1,i2,i3,0);
t251 = t202 / 0.48e2;
t252 = t205 / 0.24e2;
t253 = t235 * t237;
t254 = t253 / 0.48e2;
t255 = t237 * t241;
t256 = t255 / 0.24e2;
t257 = lapCoeff(i1m2,i2m1,i3,4);
t258 = t257 / 0.48e2;
t259 = lapCoeff(i1m2,i2m1,i3,1);
t260 = t259 / 0.24e2;
t261 = lapCoeff(i1m1,i2m2,i3,3);
t262 = t261 / 0.48e2;
t263 = 0.7e1 / 0.48e2 * t2;
t264 = lapCoeff(i1m1,i2m2,i3,0);
t265 = t264 / 0.24e2;
t271 = lapCoeff(i1m1,i2,i3,4);
t272 = t271 * t2;
t273 = t272 / 0.8e1;
t274 = lapCoeff(i1m1,i2,i3,1);
t275 = t274 * t2;
t276 = t275 / 0.4e1;
t277 = t2 * t4;
t278 = t277 / 0.8e1;
t279 = t257 / 0.8e1;
t280 = t259 / 0.4e1;
t282 = t237 * (t279 - t280);
t285 = lapCoeff(i1,i2m1,i3,3);
t286 = t285 * t2;
t287 = t286 / 0.16e2;
t288 = lapCoeff(i1,i2m1,i3,0);
t289 = t2 * t288;
t290 = t289 / 0.8e1;
t291 = t261 / 0.8e1;
t292 = t264 / 0.4e1;
t294 = t26 * (t291 - t292);
t298 = t272 / 0.16e2;
t299 = t275 / 0.8e1;
t303 = t286 / 0.8e1;
t304 = t289 / 0.4e1;
t307 = t257 / 0.2e1;
t308 = t2 / 0.2e1;
t311 = t261 / 0.2e1;
t314 = t257 / 0.4e1;
t315 = t259 / 0.2e1;
t318 = t261 / 0.4e1;
t319 = t264 / 0.2e1;
t329 = t41 / 0.144e3;
t331 = t202 / 0.4e1;
t332 = t205 / 0.2e1;
t333 = lapCoeff(i1p1,i2,i3,2);
t334 = t26 * t333;
t335 = t334 / 0.16e2;
t336 = t8 / 0.4e1;
t337 = t12 / 0.2e1;
t338 = t336 - t337;
t339 = t271 * t338;
t340 = t8 / 0.2e1;
t342 = t274 * (t340 - t12);
t343 = t8 / 0.8e1;
t344 = t12 / 0.4e1;
t347 = 0.2e1 * t4 * (t343 - t344);
t348 = t220 + t26;
t349 = t237 * t348 / 0.16e2;
t352 = t334 / 0.32e2;
t358 = lapCoeff(i1p1,i2m1,i3,3);
t360 = t358 * t26 / 0.32e2;
t361 = lapCoeff(i1p1,i2m1,i3,0);
t363 = t26 * t361 / 0.16e2;
t365 = 0.2e1 * t12;
t370 = lapCoeff(i1,i2m2,i3,3);
t371 = t370 / 0.8e1;
t372 = lapCoeff(i1,i2m2,i3,0);
t373 = t372 / 0.4e1;
t378 = t219 * t348 / 0.32e2;
t382 = t220 / 0.4e1;
t383 = t26 / 0.4e1;
t387 = t52 + t264;
t392 = t26 / 0.2e1;
t394 = t288 * (t340 - t12 + t392);
t395 = t370 / 0.2e1;
t397 = t47 * (t395 + t392 - t372);
t398 = t285 * t338;
t399 = t370 / 0.4e1;
t400 = t372 / 0.2e1;
t402 = t43 * (t399 - t400);
t403 = t26 * t387 / 0.2e1;
t411 = 0.4e1 / 0.45e2 * t41;
t412 = t370 / 0.24e2;
t413 = t4 / 0.12e2;
t414 = t372 / 0.12e2;
t424 = t220 / 0.96e2;
t425 = t63 / 0.96e2;
t435 = lapCoeff(i1p1,i2m2,i3,3);
t436 = t435 / 0.8e1;
t437 = lapCoeff(i1p1,i2m2,i3,0);
t438 = t437 / 0.4e1;
t440 = t26 * (t291 + t436 + t292 - t438);
t442 = 0.2e1 * t47;
t445 = t47 * (0.2e1 * t102 - t43 + t442 + 0.2e1 * t372);
t446 = t43 / 0.4e1;
t447 = t47 / 0.2e1;
t448 = t446 - t447;
t450 = 0.2e1 * t24 * t448;
t451 = t43 / 0.2e1;
t454 = 0.2e1 * t29 * (t451 - t47);
t455 = t2 + t63;
t456 = t285 * t455 / 0.8e1;
t458 = t4 * t455 / 0.8e1;
t459 = t2 / 0.4e1;
t460 = t63 / 0.4e1;
t462 = t288 * (t43 - t442 - t459 + t460);
t463 = t102 + t372;
t464 = t43 * t463;
t467 = lapCoeff(i1p1,i2,i3,4);
t468 = t467 * t63;
t469 = t468 / 0.16e2;
t470 = lapCoeff(i1p1,i2,i3,1);
t471 = t470 * t63;
t472 = t471 / 0.8e1;
t473 = t43 / 0.8e1;
t474 = t47 / 0.4e1;
t475 = t473 - t474;
t476 = t237 * t475;
t478 = t333 * t475;
t488 = t63 / 0.2e1;
t491 = t435 / 0.2e1;
t498 = t435 / 0.4e1;
t499 = t437 / 0.2e1;
t506 = t468 / 0.8e1;
t507 = t471 / 0.4e1;
t513 = t102 / 0.6e1;
t514 = 0.5e1 / 0.12e2 * t43;
t515 = 0.5e1 / 0.6e1 * t47;
t516 = t372 / 0.6e1;
t519 = t24 / 0.6e1;
t520 = t29 / 0.6e1;
t521 = t237 / 0.24e2;
t522 = t333 / 0.24e2;
t525 = t237 + t333;
t526 = t235 * t525 / 0.48e2;
t527 = 0.5e1 / 0.6e1 * t43;
t528 = 0.5e1 / 0.3e1 * t47;
t531 = t435 / 0.48e2;
t532 = 0.7e1 / 0.48e2 * t63;
t533 = t437 / 0.24e2;
t548 = t121 / 0.4e1;
t550 = t124 / 0.2e1;
t551 = t2 / 0.96e2;
t552 = lapCoeff(i1p2,i2m1,i3,2);
t553 = t552 / 0.96e2;
t560 = t121 / 0.2e1;
t562 = t288 * (t124 - t560 + t392);
t564 = t47 * (t395 + t392 + t372);
t565 = t548 - t550;
t566 = t285 * t565;
t568 = t43 * (t399 + t400);
t569 = t154 + t437;
t570 = t26 * t569 / 0.2e1;
t573 = t121 / 0.8e1;
t574 = t124 / 0.4e1;
t577 = 0.2e1 * t4 * (t573 - t574);
t580 = t237 * t26;
t581 = t580 / 0.16e2;
t582 = t467 * t565;
t584 = t470 * (t560 - t124);
t585 = t26 + t552;
t586 = t333 * t585 / 0.16e2;
t589 = t580 / 0.32e2;
t596 = t219 * t26 / 0.32e2;
t598 = t26 * t223 / 0.16e2;
t600 = 0.2e1 * t124;
t605 = lapCoeff(i1p2,i2m2,i3,3);
t606 = t605 / 0.8e1;
t607 = lapCoeff(i1p2,i2m2,i3,0);
t608 = t607 / 0.4e1;
t613 = t358 * t585 / 0.32e2;
t617 = t552 / 0.4e1;
t629 = t235 * t333;
t630 = t629 / 0.48e2;
t631 = t333 * t241;
t632 = t631 / 0.24e2;
t633 = lapCoeff(i1p2,i2m1,i3,1);
t634 = t633 / 0.24e2;
t635 = lapCoeff(i1p2,i2m1,i3,4);
t636 = t635 / 0.48e2;
t642 = t635 / 0.8e1;
t643 = t633 / 0.4e1;
t645 = t333 * (t642 - t643);
t649 = t285 * t63;
t650 = t649 / 0.16e2;
t651 = t63 * t288;
t652 = t651 / 0.8e1;
t654 = t26 * (t436 + t438);
t658 = t635 / 0.2e1;
t663 = t635 / 0.4e1;
t664 = t633 / 0.2e1;
t673 = t215 / 0.8e1;
t676 = t649 / 0.8e1;
t677 = t651 / 0.4e1;
t685 = t358 * t552 / 0.32e2;
t687 = t552 * t361 / 0.16e2;
t705 = lapCoeff(i1m1,i2,i3,3);
t706 = lapCoeff(i1m2,i2,i3,2);
t707 = t705 * t706;
t708 = t707 / 0.8e1;
t709 = lapCoeff(i1m1,i2,i3,0);
t710 = t706 * t709;
t711 = t710 / 0.4e1;
t712 = lapCoeff(i1m2,i2m1,i3,3);
t713 = t712 / 0.8e1;
t714 = lapCoeff(i1m2,i2m1,i3,0);
t715 = t714 / 0.4e1;
t717 = t237 * (t713 - t715);
t720 = t707 / 0.16e2;
t721 = t710 / 0.8e1;
t725 = lapCoeff(i1m1,i2p1,i3,2);
t731 = t712 / 0.4e1;
t732 = t714 / 0.2e1;
t736 = t712 / 0.2e1;
t745 = t712 / 0.48e2;
t746 = t219 / 0.48e2;
t747 = t714 / 0.24e2;
t748 = t223 / 0.48e2;
t752 = t235 * (t706 + t237) / 0.48e2;
t753 = t706 / 0.48e2;
t755 = t241 * (t753 + t521);
t760 = t253 / 0.4e1;
t761 = lapCoeff(i1,i2p1,i3,2);
t762 = t237 * t761;
t763 = t762 / 0.16e2;
t764 = t255 / 0.2e1;
t765 = t9 + t237;
t766 = t26 * t765 / 0.16e2;
t767 = t219 / 0.4e1;
t768 = t223 / 0.2e1;
t769 = t767 - t768;
t770 = t285 * t769;
t771 = t219 / 0.8e1;
t772 = t223 / 0.4e1;
t775 = 0.2e1 * t4 * (t771 - t772);
t776 = t219 / 0.2e1;
t778 = t288 * (t776 - t223);
t781 = t762 / 0.32e2;
t787 = lapCoeff(i1m1,i2p1,i3,4);
t789 = t787 * t237 / 0.32e2;
t790 = lapCoeff(i1m1,i2p1,i3,1);
t792 = t790 * t237 / 0.16e2;
t793 = lapCoeff(i1m2,i2,i3,4);
t794 = t793 / 0.8e1;
t795 = lapCoeff(i1m2,i2,i3,1);
t796 = t795 / 0.4e1;
t802 = 0.2e1 * t223;
t810 = t8 * t765 / 0.32e2;
t811 = t9 / 0.4e1;
t812 = t237 / 0.4e1;
t816 = t259 + t714;
t821 = t793 / 0.2e1;
t822 = t237 / 0.2e1;
t824 = t709 * (t821 - t795 + t822);
t826 = t274 * (t776 + t822 - t223);
t827 = t793 / 0.4e1;
t828 = t795 / 0.2e1;
t830 = t705 * (t827 - t828);
t831 = t271 * t769;
t832 = t237 * t816 / 0.2e1;
t840 = t793 / 0.24e2;
t841 = t795 / 0.12e2;
t851 = t9 / 0.96e2;
t852 = t725 / 0.96e2;
t860 = 0.19e2 / 0.36e2 * t41;
t861 = 0.5e1 / 0.6e1 * t274;
t862 = 0.5e1 / 0.12e2 * t271;
t863 = t333 / 0.48e2;
t866 = 0.5e1 / 0.12e2 * t285;
t867 = t761 / 0.48e2;
t868 = 0.5e1 / 0.6e1 * t288;
t871 = t787 / 0.48e2;
t873 = t790 / 0.48e2;
t874 = t358 / 0.48e2;
t876 = t361 / 0.48e2;
t879 = 0.5e1 / 0.6e1 * t271;
t880 = 0.5e1 / 0.3e1 * t274;
t881 = 0.23e2 / 0.24e2 * t4;
t884 = 0.5e1 / 0.6e1 * t285;
t885 = t761 / 0.24e2;
t886 = 0.5e1 / 0.3e1 * t288;
t892 = lapCoeff(i1,i2p1,i3,1);
t893 = t892 * t4;
t894 = t893 / 0.4e1;
t895 = lapCoeff(i1,i2p1,i3,4);
t896 = t895 * t4;
t897 = t896 / 0.8e1;
t898 = t271 / 0.8e1;
t899 = t274 / 0.4e1;
t901 = t26 * (t51 + t898 + t53 - t899);
t902 = 0.2e1 * t288;
t904 = t288 * (t365 - t285 + t802 + t902);
t905 = t271 / 0.2e1;
t906 = t4 / 0.2e1;
t909 = 0.2e1 * t241 * (t905 - t274 + t906);
t910 = t285 / 0.2e1;
t913 = 0.2e1 * t29 * (t910 + t906 - t288);
t914 = t271 / 0.4e1;
t915 = t274 / 0.2e1;
t916 = t914 - t915;
t918 = 0.2e1 * t235 * t916;
t920 = t761 * (t898 - t899);
t921 = t44 + t4;
t922 = t43 * t921 / 0.8e1;
t923 = t285 / 0.4e1;
t924 = t288 / 0.2e1;
t925 = t923 - t924;
t927 = 0.2e1 * t24 * t925;
t928 = t12 + t223;
t930 = t4 * t928;
t931 = t44 / 0.4e1;
t932 = t4 / 0.4e1;
t934 = t47 * (t285 - t931 + t932 - t902);
t935 = t285 * t928;
t936 = t894 - t897 + t901 - t904 + t909 + t913 - t918 - t920 + t922 - t927 + t930 + t934 + t935;
t938 = lapCoeff(i1p1,i2,i3,0);
t939 = t4 * t938;
t940 = t939 / 0.4e1;
t941 = lapCoeff(i1p1,i2,i3,3);
t942 = t941 * t4;
t943 = t942 / 0.8e1;
t944 = 0.2e1 * t274;
t946 = t274 * (t365 - t271 + t944 + t802);
t947 = t285 / 0.8e1;
t948 = t288 / 0.4e1;
t950 = t237 * (t713 + t947 + t715 - t948);
t951 = t706 + t4;
t952 = t705 * t951 / 0.8e1;
t954 = t333 * (t947 - t948);
t955 = t706 / 0.4e1;
t957 = t709 * (t271 - t944 - t955 + t932);
t958 = t271 * t928;
t959 = t940 - t943 - t946 + t950 + t909 + t913 - t918 - t927 + t952 + t930 - t954 + t957 + t958;
t961 = lapCoeff(i1p1,i2p1,i3,2);
t962 = t4 * t961;
t987 = 0.4e1 * t12;
t988 = 0.4e1 * t223;
t995 = t962 / 0.64e2 - t8 * (t66 + t914 + t67 - t915) / 0.4e1 + t2 * (t1 + t706 + t44 + t4) / 0.64e2 - t219 * (t731 + t923 + t732 - t924) / 0.4e1 - t787 * t916 / 0.4e1 + t790 * (t905 - t274) / 0.4e1 + t725 * t951 / 0.64e2 + t63 * t921 / 0.64e2 - t358 * t925 / 0.4e1 + t361 * (t910 - t288) / 0.4e1 + t12 * (t71 - t905 + t52 + t987 + t274 + t988) / 0.4e1 + t223 * (t987 + t736 - t910 + t714 + t988 + t288) / 0.4e1;
t997 = t896 / 0.16e2;
t998 = t893 / 0.8e1;
t1002 = t922 / 0.2e1;
t1007 = t942 / 0.16e2;
t1008 = t939 / 0.8e1;
t1011 = t952 / 0.2e1;
t1021 = 0.19e2 / 0.12e2 * t24;
t1022 = 0.19e2 / 0.6e1 * t29;
t1023 = t333 / 0.2e1;
t1026 = t895 / 0.12e2;
t1028 = t892 / 0.12e2;
t1032 = t358 / 0.4e1;
t1033 = t111 / 0.96e2;
t1034 = t961 / 0.96e2;
t1035 = t361 / 0.2e1;
t1038 = t235 * t525 / 0.4e1;
t1039 = t895 / 0.24e2;
t1041 = t892 / 0.24e2;
t1048 = 0.3e1 / 0.4e1 * t24;
t1049 = 0.3e1 / 0.2e1 * t29;
t1052 = lapCoeff(i1p1,i2p1,i3,4);
t1054 = t1052 * t333 / 0.32e2;
t1055 = lapCoeff(i1p1,i2p1,i3,1);
t1057 = t1055 * t333 / 0.16e2;
t1058 = t24 / 0.8e1;
t1059 = t29 / 0.4e1;
t1060 = t127 + t1058 + t128 - t1059;
t1065 = 0.2e1 * t361;
t1069 = t1058 - t1059;
t1074 = t111 + t333;
t1076 = t121 * t1074 / 0.32e2;
t1080 = t111 / 0.4e1;
t1081 = t333 / 0.4e1;
t1088 = t47 + t288;
t1093 = t789 + t1054 - t792 - t1057 + t2 * t1060 / 0.4e1 + t63 * t1060 / 0.4e1 - t361 * (t442 - t358 + t902 + t1065) / 0.4e1 + t725 * t1069 / 0.4e1 + t961 * t1069 / 0.4e1 + t810 + t1076 + t12 * (t219 - t811 + t812 + t802) / 0.4e1 + t124 * (t358 - t1080 + t1081 - t1065) / 0.4e1 + t223 * (t442 + t219 + t802 + t902) / 0.4e1 + t219 * t1088 / 0.4e1 + t358 * t1088 / 0.4e1;
t1095 = t24 / 0.4e1;
t1096 = t29 / 0.2e1;
t1098 = t43 * (t139 + t1095 + t140 - t1096);
t1099 = 0.2e1 * t29;
t1102 = 0.2e1 * t29 * (t442 - t24 + t1099 + t902);
t1104 = t26 * (t9 + t237 + t111 + t333) / 0.16e2;
t1106 = t285 * (t767 + t1032 + t768 - t1035);
t1107 = t358 / 0.8e1;
t1108 = t361 / 0.4e1;
t1111 = 0.2e1 * t4 * (t771 + t1107 + t772 - t1108);
t1112 = t1095 - t1096;
t1113 = t895 * t1112;
t1114 = t24 / 0.2e1;
t1116 = t892 * (t1114 - t29);
t1118 = t235 * t525 / 0.4e1;
t1119 = t761 * t525 / 0.16e2;
t1122 = 0.2e1 * t241 * (t24 - t1099 - t812 + t1081);
t1124 = 0.2e1 * t24 * t1088;
t1125 = 0.4e1 * t47;
t1126 = 0.4e1 * t288;
t1128 = t47 * (t136 - t1114 + t102 + t1125 + t29 + t1126);
t1129 = t358 / 0.2e1;
t1131 = t288 * (t1125 + t776 - t1129 + t223 + t1126 + t361);
t1132 = t1098 - t1102 - t1104 + t1106 + t1111 - t1113 + t1116 + t1118 + t1119 + t1122 + t1124 - t1128 - t1131;
t1137 = t274 * (t776 + t822 + t223);
t1139 = t709 * (t29 - t1114 + t822);
t1141 = t938 * (t1114 - t29 + t1023);
t1143 = t470 * (t1129 + t1023 - t361);
t1144 = t705 * t1112;
t1145 = t941 * t1112;
t1147 = t271 * (t767 + t768);
t1149 = t467 * (t1032 - t1035);
t1150 = t237 * t1088 / 0.2e1;
t1151 = t333 * t1088 / 0.2e1;
t1152 = t1111 - t1137 - t1102 - t1139 + t1141 + t1143 + t1144 - t1145 + t1147 - t1149 + t1118 - t1150 + t1151 + t1122 + t1124;
t1160 = t1052 / 0.48e2;
t1162 = t1055 / 0.48e2;
t1163 = lapCoeff(i1p2,i2m1,i3,3);
t1164 = t1163 / 0.48e2;
t1166 = lapCoeff(i1p2,i2m1,i3,0);
t1167 = t1166 / 0.24e2;
t1172 = 0.5e1 / 0.12e2 * t467;
t1173 = 0.5e1 / 0.6e1 * t470;
t1174 = t237 / 0.48e2;
t1175 = lapCoeff(i1p2,i2,i3,2);
t1176 = t1175 / 0.48e2;
t1179 = 0.5e1 / 0.3e1 * t470;
t1180 = 0.5e1 / 0.6e1 * t467;
t1188 = t725 * t4;
t1190 = t467 / 0.4e1;
t1191 = t470 / 0.2e1;
t1198 = t1163 / 0.4e1;
t1199 = t1166 / 0.2e1;
t1206 = t1190 - t1191;
t1209 = t467 / 0.2e1;
t1215 = t4 + t1175;
t1218 = t923 + t924;
t1221 = 0.4e1 * t124;
t1222 = 0.4e1 * t361;
t1226 = t1163 / 0.2e1;
t1230 = t1188 / 0.64e2 - t121 * (t165 + t1190 + t166 - t1191) / 0.4e1 + t63 * (t44 + t4 + t174 + t1175) / 0.64e2 - t358 * (t923 + t1198 + t924 - t1199) / 0.4e1 - t223 * (t910 + t288) / 0.4e1 - t1052 * t1206 / 0.4e1 + t1055 * (t1209 - t470) / 0.4e1 + t2 * t921 / 0.64e2 + t961 * t1215 / 0.64e2 - t219 * t1218 / 0.4e1 + t124 * (t170 - t1209 + t154 + t1221 + t470 + t1222) / 0.4e1 + t361 * (t1221 + t910 - t1226 + t288 + t1222 + t1166) / 0.4e1;
t1232 = t467 / 0.8e1;
t1233 = t470 / 0.4e1;
t1235 = t26 * (t153 + t1232 + t155 - t1233);
t1238 = t761 * (t1232 - t1233);
t1241 = t47 * (t285 - t931 + t932 + t902);
t1244 = t288 * (t600 + t285 + t902 + t1065);
t1246 = t124 + t361;
t1247 = t285 * t1246;
t1253 = 0.2e1 * t241 * (t470 - t1209 + t906);
t1256 = 0.2e1 * t29 * (t910 + t906 + t288);
t1258 = 0.2e1 * t235 * t1206;
t1260 = 0.2e1 * t24 * t1218;
t1262 = t4 * t1246;
t1263 = t1253 + t1256 - t897 + t894 + t1235 - t1258 - t1238 + t922 - t1260 + t1262 + t1241 + t1244 + t1247;
t1265 = t705 * t4;
t1266 = t1265 / 0.16e2;
t1267 = t4 * t709;
t1268 = t1267 / 0.8e1;
t1269 = 0.2e1 * t470;
t1271 = t470 * (t600 - t467 + t1269 + t1065);
t1273 = t1163 / 0.8e1;
t1274 = t1166 / 0.4e1;
t1276 = t333 * (t947 + t1273 + t948 - t1274);
t1278 = t941 * t1215 / 0.8e1;
t1279 = t1278 / 0.2e1;
t1281 = t237 * (t947 + t948);
t1283 = t1175 / 0.4e1;
t1285 = t938 * (t467 - t1269 - t932 + t1283);
t1287 = t467 * t1246;
t1291 = t1265 / 0.8e1;
t1292 = t1267 / 0.4e1;
t1293 = t1253 + t1256 + t1291 + t1292 + t1271 - t1276 - t1258 - t1260 - t1278 + t1262 + t1281 - t1285 - t1287;
t1299 = lapCoeff(i1p2,i2,i3,1);
t1300 = lapCoeff(i1p2,i2,i3,4);
t1301 = t1300 / 0.2e1;
t1303 = t938 * (t1299 - t1301 + t1023);
t1305 = t470 * (t1129 + t1023 + t361);
t1306 = t629 / 0.4e1;
t1307 = t631 / 0.2e1;
t1308 = t1300 / 0.4e1;
t1309 = t1299 / 0.2e1;
t1311 = t941 * (t1308 - t1309);
t1312 = t1032 + t1035;
t1313 = t467 * t1312;
t1314 = t633 + t1166;
t1315 = t333 * t1314 / 0.2e1;
t1318 = 0.2e1 * t4 * (t1107 + t1108);
t1323 = t1300 / 0.8e1;
t1324 = t1299 / 0.4e1;
t1343 = t761 * t333;
t1344 = t1343 / 0.16e2;
t1346 = t288 * (t1129 + t361);
t1347 = t26 * t1074 / 0.16e2;
t1348 = t285 * t1312;
t1351 = t1343 / 0.32e2;
t1360 = t1299 / 0.12e2;
t1361 = t1300 / 0.24e2;
t1381 = t235 * (t333 + t1175) / 0.48e2;
t1383 = t241 * (t522 + t1176);
t1387 = t941 * t1175;
t1388 = t1387 / 0.16e2;
t1389 = t1175 * t938;
t1390 = t1389 / 0.8e1;
t1392 = t333 * (t1273 + t1274);
t1396 = t1387 / 0.8e1;
t1397 = t1389 / 0.4e1;
t1417 = lapCoeff(i1m2,i2p1,i3,2);
t1420 = lapCoeff(i1m2,i2,i3,3);
t1421 = t1420 / 0.24e2;
t1422 = t705 / 0.24e2;
t1423 = lapCoeff(i1m2,i2,i3,0);
t1424 = t1423 / 0.12e2;
t1425 = t709 / 0.24e2;
t1428 = t705 / 0.12e2;
t1429 = t709 / 0.12e2;
t1435 = t235 / 0.60e2;
t1436 = t241 / 0.90e2;
t1439 = lapCoeff(i1m1,i2p1,i3,3);
t1441 = t1439 * t1417 / 0.32e2;
t1442 = lapCoeff(i1m1,i2p1,i3,0);
t1444 = t1417 * t1442 / 0.16e2;
t1445 = t1420 / 0.8e1;
t1446 = t1423 / 0.4e1;
t1447 = t1445 - t1446;
t1456 = t1420 / 0.4e1;
t1457 = t1423 / 0.2e1;
t1460 = t1420 / 0.2e1;
t1463 = t237 * (t220 + t1417) / 0.16e2 - t705 * (t1456 - t1457) + t709 * (t1460 - t1423);
t1472 = lapCoeff(i1m2,i2p1,i3,4);
t1473 = t1472 / 0.8e1;
t1474 = lapCoeff(i1m2,i2p1,i3,1);
t1475 = t1474 / 0.4e1;
t1477 = t237 * (t279 + t1473 + t280 - t1475);
t1480 = 0.2e1 * t709;
t1482 = t709 * (0.2e1 * t795 - t705 + 0.2e1 * t1423 + t1480);
t1483 = t2 + t725;
t1484 = t271 * t1483 / 0.8e1;
t1486 = t4 * t1483 / 0.8e1;
t1487 = t705 / 0.4e1;
t1488 = t709 / 0.2e1;
t1489 = t1487 - t1488;
t1491 = 0.2e1 * t235 * t1489;
t1492 = t705 / 0.2e1;
t1495 = 0.2e1 * t241 * (t1492 - t709);
t1496 = t725 / 0.4e1;
t1498 = t274 * (t705 - t459 + t1496 - t1480);
t1499 = t795 + t1423;
t1500 = t705 * t1499;
t1503 = lapCoeff(i1,i2p1,i3,3);
t1504 = t1503 * t725;
t1505 = t1504 / 0.16e2;
t1506 = lapCoeff(i1,i2p1,i3,0);
t1507 = t725 * t1506;
t1508 = t1507 / 0.8e1;
t1509 = t705 / 0.8e1;
t1510 = t709 / 0.4e1;
t1511 = t1509 - t1510;
t1512 = t26 * t1511;
t1514 = t761 * t1511;
t1524 = t1472 / 0.2e1;
t1525 = t725 / 0.2e1;
t1532 = t1472 / 0.4e1;
t1533 = t1474 / 0.2e1;
t1542 = t1504 / 0.8e1;
t1543 = t1507 / 0.4e1;
t1549 = t235 / 0.6e1;
t1550 = t241 / 0.6e1;
t1553 = t795 / 0.6e1;
t1554 = 0.5e1 / 0.12e2 * t705;
t1555 = t1423 / 0.6e1;
t1556 = 0.5e1 / 0.6e1 * t709;
t1559 = t26 + t761;
t1560 = t24 * t1559 / 0.48e2;
t1561 = t1472 / 0.48e2;
t1562 = t1474 / 0.24e2;
t1563 = 0.7e1 / 0.48e2 * t725;
t1566 = 0.5e1 / 0.6e1 * t705;
t1567 = 0.5e1 / 0.3e1 * t709;
t1574 = 0.19e2 / 0.12e2 * t235;
t1575 = t761 / 0.2e1;
t1576 = 0.19e2 / 0.6e1 * t241;
t1579 = t790 / 0.2e1;
t1580 = t787 / 0.4e1;
t1581 = t1417 / 0.96e2;
t1585 = t941 / 0.12e2;
t1587 = t938 / 0.12e2;
t1590 = t24 * t1559 / 0.4e1;
t1592 = t941 / 0.24e2;
t1594 = t938 / 0.24e2;
t1600 = 0.3e1 / 0.4e1 * t235;
t1601 = 0.3e1 / 0.2e1 * t241;
t1604 = lapCoeff(i1p1,i2p1,i3,3);
t1606 = t1604 * t761 / 0.32e2;
t1607 = lapCoeff(i1p1,i2p1,i3,0);
t1609 = t761 * t1607 / 0.16e2;
t1610 = 0.2e1 * t790;
t1614 = t235 / 0.8e1;
t1615 = t241 / 0.4e1;
t1616 = t1445 + t1614 + t1446 - t1615;
t1621 = t1417 + t761;
t1623 = t1439 * t1621 / 0.32e2;
t1624 = t1614 - t1615;
t1635 = t1417 / 0.4e1;
t1636 = t761 / 0.4e1;
t1640 = t274 + t709;
t1645 = t360 + t1606 - t363 - t1609 - t790 * (t944 - t787 + t1610 + t1480) / 0.4e1 + t2 * t1616 / 0.4e1 + t725 * t1616 / 0.4e1 + t378 + t1623 + t63 * t1624 / 0.4e1 + t961 * t1624 / 0.4e1 + t12 * (t8 + t365 + t944 + t1480) / 0.4e1 + t223 * (t8 + t365 - t382 + t383) / 0.4e1 + t1442 * (t787 - t1610 - t1635 + t1636) / 0.4e1 + t8 * t1640 / 0.4e1 + t787 * t1640 / 0.4e1;
t1648 = t271 * (t336 + t1580 + t337 - t1579);
t1649 = t787 / 0.8e1;
t1650 = t790 / 0.4e1;
t1653 = 0.2e1 * t4 * (t343 + t1649 + t344 - t1650);
t1655 = t237 * (t220 + t1417 + t26 + t761) / 0.16e2;
t1656 = t235 / 0.4e1;
t1657 = t241 / 0.2e1;
t1659 = t705 * (t1456 + t1656 + t1457 - t1657);
t1660 = 0.2e1 * t241;
t1663 = 0.2e1 * t241 * (t944 - t235 + t1480 + t1660);
t1665 = t24 * t1559 / 0.4e1;
t1666 = t333 * t1559 / 0.16e2;
t1667 = t1656 - t1657;
t1668 = t941 * t1667;
t1669 = t235 / 0.2e1;
t1671 = t938 * (t1669 - t241);
t1674 = 0.2e1 * t29 * (t235 - t383 + t1636 - t1660);
t1676 = 0.2e1 * t235 * t1640;
t1677 = t787 / 0.2e1;
t1678 = 0.4e1 * t274;
t1679 = 0.4e1 * t709;
t1681 = t274 * (t340 - t1677 + t12 + t1678 + t790 + t1679);
t1683 = t709 * (t1678 + t1460 - t1669 + t1423 + t1679 + t241);
t1684 = t1648 + t1653 - t1655 + t1659 - t1663 + t1665 + t1666 - t1668 + t1671 + t1674 + t1676 - t1681 - t1683;
t1689 = t47 * (t392 - t1669 + t241);
t1691 = t288 * (t340 + t12 + t392);
t1693 = t1506 * (t1677 - t790 + t1575);
t1695 = t892 * (t1669 + t1575 - t241);
t1697 = t285 * (t336 + t337);
t1699 = t1503 * (t1580 - t1579);
t1700 = t43 * t1667;
t1701 = t895 * t1667;
t1702 = t26 * t1640 / 0.2e1;
t1703 = t761 * t1640 / 0.2e1;
t1704 = t1653 - t1689 - t1691 - t1663 + t1693 + t1695 + t1697 - t1699 + t1665 + t1700 - t1701 - t1702 + t1703 + t1674 + t1676;
t1712 = t2 + t725 + t63 + t961;
t1714 = 0.5e1 / 0.12e2 * t895;
t1715 = 0.5e1 / 0.6e1 * t892;
t1718 = 0.5e1 / 0.12e2 * t941;
t1719 = 0.5e1 / 0.6e1 * t938;
t1722 = 0.5e1 / 0.6e1 * t895;
t1723 = 0.6e1 * t29;
t1724 = 0.5e1 / 0.3e1 * t892;
t1725 = 0.6e1 * t241;
t1728 = 0.5e1 / 0.6e1 * t941;
t1729 = 0.5e1 / 0.3e1 * t938;
t1735 = t29 + t241;
t1737 = t941 / 0.8e1;
t1738 = t938 / 0.4e1;
t1739 = t1509 + t1737 + t1510 - t1738;
t1740 = t26 * t1739;
t1741 = 0.2e1 * t892;
t1743 = t892 * (t1099 - t895 + t1741 + t1660);
t1744 = t761 * t1739;
t1745 = t725 + t961;
t1746 = t1503 * t1745 / 0.8e1;
t1748 = t47 * (t43 + t442 + t1099 + t1660);
t1750 = t288 * (t43 + t442 - t459 + t460);
t1751 = t961 / 0.4e1;
t1753 = t1506 * (t895 - t1741 - t1496 + t1751);
t1754 = t43 * t1735;
t1755 = t895 * t1735;
t1758 = t895 / 0.8e1;
t1759 = t892 / 0.4e1;
t1760 = t473 + t1758 + t474 - t1759;
t1761 = t237 * t1760;
t1762 = t333 * t1760;
t1763 = 0.2e1 * t938;
t1765 = t938 * (t1099 - t941 + t1660 + t1763);
t1766 = t63 + t961;
t1767 = t467 * t1766 / 0.8e1;
t1769 = t274 * (t705 - t459 + t1496 + t1480);
t1771 = t470 * (t941 - t460 + t1751 - t1763);
t1773 = t709 * (t1099 + t705 + t1480 + t1660);
t1774 = t705 * t1735;
t1775 = t941 * t1735;
t1778 = t895 / 0.4e1;
t1779 = t892 / 0.2e1;
t1782 = 0.2e1 * t24 * (t446 + t1778 + t447 - t1779);
t1784 = t4 * t1712 / 0.8e1;
t1785 = t941 / 0.4e1;
t1786 = t938 / 0.2e1;
t1789 = 0.2e1 * t235 * (t1487 + t1785 + t1488 - t1786);
t1790 = t895 / 0.2e1;
t1791 = 0.4e1 * t29;
t1792 = 0.4e1 * t241;
t1795 = 0.2e1 * t29 * (t451 - t1790 + t47 + t1791 + t892 + t1792);
t1796 = t941 / 0.2e1;
t1799 = 0.2e1 * t241 * (t1791 + t1492 - t1796 + t709 + t1792 + t938);
t1800 = t1743 - t1782 + t1784 - t1789 + t1740 - t1744 + t456 - t1746 + t1748 + t1750 - t1753 + t1754 - t1755 + t1795 + t1799;
t1802 = t1761 - t1782 - t1762 + t1784 - t1789 + t1765 + t1484 - t1767 + t1769 - t1771 + t1773 + t1774 - t1775 + t1795 + t1799;
t1816 = t961 / 0.2e1;
t1821 = t446 + t447;
t1823 = t1778 - t1779;
t1827 = t1487 + t1488;
t1830 = t1785 - t1786;
t1837 = t223 * (t451 + t47 - t308) - t1442 * (t892 - t1790 + t1525) - t361 * (t451 + t47 + t488) + t12 * (t1492 - t308 + t709) - t790 * (t1492 + t1525 + t709) - t124 * (t488 - t1796 + t938) - t1607 * (t1790 - t892 + t1816) - t1055 * (t1796 + t1816 - t938) + t219 * t1821 + t1439 * t1823 + t358 * t1821 + t1604 * t1823 + t8 * t1827 + t787 * t1827 + t121 * t1830 + t1052 * t1830 - t2 * t1735 / 0.2e1 - t725 * t1735 / 0.2e1 - t63 * t1735 / 0.2e1 - t961 * t1735 / 0.2e1;
t1847 = t1052 / 0.4e1;
t1848 = t1055 / 0.2e1;
t1849 = lapCoeff(i1p2,i2p1,i3,2);
t1850 = t1849 / 0.96e2;
t1854 = lapCoeff(i1p2,i2,i3,3);
t1855 = t1854 / 0.24e2;
t1857 = lapCoeff(i1p2,i2,i3,0);
t1858 = t1857 / 0.12e2;
t1869 = t1439 * t761 / 0.32e2;
t1871 = t761 * t1442 / 0.16e2;
t1872 = 0.2e1 * t1055;
t1876 = t1854 / 0.8e1;
t1877 = t1857 / 0.4e1;
t1878 = t1614 + t1876 + t1615 - t1877;
t1883 = t761 + t1849;
t1885 = t1604 * t1883 / 0.32e2;
t1886 = t1614 + t1615;
t1897 = t1849 / 0.4e1;
t1901 = t470 + t938;
t1906 = t596 + t1869 + t598 + t1871 - t1055 * (t1269 - t1052 + t1872 + t1763) / 0.4e1 + t63 * t1878 / 0.4e1 + t961 * t1878 / 0.4e1 + t613 + t1885 + t2 * t1886 / 0.4e1 + t725 * t1886 / 0.4e1 + t124 * (t121 + t600 + t1269 + t1763) / 0.4e1 + t361 * (t121 + t600 - t383 + t617) / 0.4e1 + t1607 * (t1052 - t1872 - t1636 + t1897) / 0.4e1 + t121 * t1901 / 0.4e1 + t1052 * t1901 / 0.4e1;
t1909 = t467 * (t548 + t1847 + t550 - t1848);
t1911 = t333 * (t26 + t761 + t552 + t1849) / 0.16e2;
t1912 = t1854 / 0.4e1;
t1913 = t1857 / 0.2e1;
t1915 = t941 * (t1656 + t1912 + t1657 - t1913);
t1917 = t709 * (t1669 + t241);
t1918 = t237 * t1559 / 0.16e2;
t1919 = t1656 + t1657;
t1920 = t705 * t1919;
t1921 = t1052 / 0.2e1;
t1922 = 0.4e1 * t470;
t1923 = 0.4e1 * t938;
t1925 = t470 * (t560 - t1921 + t124 + t1922 + t1055 + t1923);
t1926 = t1854 / 0.2e1;
t1928 = t938 * (t1922 + t1669 - t1926 + t241 + t1923 + t1857);
t1932 = t288 * (t560 + t124 - t392);
t1934 = t1506 * (t1055 - t1921 + t1575);
t1936 = t47 * (t1669 - t392 + t241);
t1938 = t892 * (t1669 + t1575 + t241);
t1940 = t285 * (t548 + t550);
t1942 = t1503 * (t1847 - t1848);
t1943 = t43 * t1919;
t1944 = t895 * t1919;
t1945 = t26 * t1901 / 0.2e1;
t1946 = t761 * t1901 / 0.2e1;
t1949 = t1052 / 0.8e1;
t1950 = t1055 / 0.4e1;
t1953 = 0.2e1 * t4 * (t573 + t1949 + t574 - t1950);
t1956 = 0.2e1 * t29 * (t235 - t383 + t1636 + t1660);
t1959 = 0.2e1 * t241 * (t1269 + t235 + t1660 + t1763);
t1961 = 0.2e1 * t235 * t1901;
t1962 = t1953 - t1909 + t1911 - t1915 + t1917 + t1665 - t1918 + t1920 + t1956 + t1959 + t1961 + t1925 + t1928;
t1964 = t1932 + t1934 + t1936 + t1938 + t1953 + t1940 - t1942 + t1665 + t1943 - t1944 - t1945 + t1946 + t1956 + t1959 + t1961;
t1972 = t1299 / 0.6e1;
t1973 = t1857 / 0.6e1;
t1976 = lapCoeff(i1p2,i2p1,i3,4);
t1977 = t1976 / 0.48e2;
t1978 = lapCoeff(i1p2,i2p1,i3,1);
t1979 = t1978 / 0.24e2;
t1980 = 0.7e1 / 0.48e2 * t961;
t1988 = t1503 * t961;
t1989 = t1988 / 0.16e2;
t1990 = t961 * t1506;
t1991 = t1990 / 0.8e1;
t1992 = t1737 + t1738;
t1993 = t26 * t1992;
t1995 = t761 * t1992;
t1999 = t1976 / 0.8e1;
t2000 = t1978 / 0.4e1;
t2002 = t333 * (t642 + t1999 + t643 - t2000);
t2005 = 0.2e1 * t241 * (t1796 + t938);
t2007 = t4 * t1766 / 0.8e1;
t2008 = t1785 + t1786;
t2010 = 0.2e1 * t235 * t2008;
t2012 = t470 * (t941 - t460 + t1751 + t1763);
t2016 = t938 * (0.2e1 * t1299 + t941 + t1763 + 0.2e1 * t1857);
t2017 = t1299 + t1857;
t2018 = t941 * t2017;
t2023 = t1976 / 0.2e1;
t2032 = t1976 / 0.4e1;
t2033 = t1978 / 0.2e1;
t2044 = t1988 / 0.8e1;
t2045 = t1990 / 0.4e1;
t2071 = t938 * (t1926 + t1857) - t333 * (t552 + t1849) / 0.16e2 + t941 * (t1912 + t1913);
t2074 = t1604 * t1849 / 0.32e2;
t2076 = t1849 * t1607 / 0.16e2;
t2077 = t1876 + t1877;
t2089 = lapCoeff(i1m2,i2p1,i3,3);
t2090 = t2089 / 0.8e1;
t2091 = lapCoeff(i1m2,i2p1,i3,0);
t2092 = t2091 / 0.4e1;
t2094 = t237 * (t2090 - t2092);
t2102 = lapCoeff(i1m2,i2p2,i3,2);
t2106 = t2089 / 0.4e1;
t2107 = t2091 / 0.2e1;
t2111 = t2089 / 0.2e1;
t2120 = t2089 / 0.48e2;
t2121 = t1439 / 0.48e2;
t2122 = t2091 / 0.24e2;
t2123 = t1442 / 0.48e2;
t2139 = t1439 / 0.4e1;
t2140 = lapCoeff(i1m1,i2p2,i3,2);
t2141 = t2140 / 0.96e2;
t2143 = t1442 / 0.2e1;
t2150 = t709 * (t821 + t795 + t822);
t2151 = t1439 / 0.2e1;
t2153 = t274 * (t822 - t2151 + t1442);
t2155 = t705 * (t827 + t828);
t2156 = t2139 - t2143;
t2157 = t271 * t2156;
t2158 = t1474 + t2091;
t2159 = t237 * t2158 / 0.2e1;
t2162 = t1439 / 0.8e1;
t2163 = t1442 / 0.4e1;
t2166 = 0.2e1 * t4 * (t2162 - t2163);
t2169 = t237 + t2140;
t2170 = t761 * t2169 / 0.16e2;
t2171 = t1503 * t2156;
t2173 = t1506 * (t2151 - t1442);
t2182 = t8 * t237 / 0.32e2;
t2184 = t12 * t237 / 0.16e2;
t2185 = lapCoeff(i1m2,i2p2,i3,4);
t2186 = t2185 / 0.8e1;
t2187 = lapCoeff(i1m2,i2p2,i3,1);
t2188 = t2187 / 0.4e1;
t2194 = 0.2e1 * t1442;
t2202 = t787 * t2169 / 0.32e2;
t2203 = t2140 / 0.4e1;
t2216 = lapCoeff(i1m1,i2p2,i3,4);
t2217 = t2216 / 0.48e2;
t2218 = lapCoeff(i1m1,i2p2,i3,1);
t2219 = t2218 / 0.24e2;
t2220 = t1604 / 0.48e2;
t2222 = t1607 / 0.48e2;
t2225 = 0.5e1 / 0.12e2 * t1503;
t2226 = t26 / 0.48e2;
t2227 = lapCoeff(i1,i2p2,i3,2);
t2228 = t2227 / 0.48e2;
t2229 = 0.5e1 / 0.6e1 * t1506;
t2236 = 0.5e1 / 0.6e1 * t1503;
t2237 = 0.5e1 / 0.3e1 * t1506;
t2244 = t2216 / 0.4e1;
t2245 = t2218 / 0.2e1;
t2252 = t1503 / 0.4e1;
t2253 = t1506 / 0.2e1;
t2260 = t914 + t915;
t2265 = t4 + t2227;
t2268 = t2252 - t2253;
t2271 = t1503 / 0.2e1;
t2275 = t2216 / 0.2e1;
t2276 = 0.4e1 * t790;
t2277 = 0.4e1 * t1442;
t2284 = t215 / 0.64e2 - t787 * (t914 + t2244 + t915 - t2245) / 0.4e1 + t725 * (t706 + t2102 + t4 + t2227) / 0.64e2 - t1439 * (t2106 + t2252 + t2107 - t2253) / 0.4e1 - t12 * (t905 + t274) / 0.4e1 - t8 * t2260 / 0.4e1 + t2 * t951 / 0.64e2 + t961 * t2265 / 0.64e2 - t1604 * t2268 / 0.4e1 + t1607 * (t2271 - t1506) / 0.4e1 + t790 * (t905 - t2275 + t274 + t2276 + t2218 + t2277) / 0.4e1 + t1442 * (t2276 + t2111 - t2271 + t2091 + t2277 + t1506) / 0.4e1;
t2286 = t1503 / 0.8e1;
t2287 = t1506 / 0.4e1;
t2289 = t237 * (t2090 + t2286 + t2092 - t2287);
t2292 = t333 * (t2286 - t2287);
t2295 = t274 * (t271 + t944 + t1610 + t2194);
t2298 = t709 * (t271 + t944 - t955 + t932);
t2300 = t790 + t1442;
t2301 = t271 * t2300;
t2307 = 0.2e1 * t241 * (t905 + t274 + t906);
t2310 = 0.2e1 * t29 * (t906 - t2271 + t1506);
t2312 = 0.2e1 * t235 * t2260;
t2314 = 0.2e1 * t24 * t2268;
t2316 = t4 * t2300;
t2317 = t2307 + t2310 - t943 + t940 + t2289 - t2312 - t2314 + t952 + t2316 - t2292 + t2295 + t2298 + t2301;
t2319 = t43 * t4;
t2320 = t2319 / 0.16e2;
t2321 = t47 * t4;
t2322 = t2321 / 0.8e1;
t2323 = t2216 / 0.8e1;
t2324 = t2218 / 0.4e1;
t2326 = t761 * (t898 + t2323 + t899 - t2324);
t2328 = 0.2e1 * t1506;
t2330 = t1506 * (t1610 - t1503 + t2194 + t2328);
t2333 = t26 * (t898 + t899);
t2335 = t895 * t2265 / 0.8e1;
t2336 = t2335 / 0.2e1;
t2337 = t2227 / 0.4e1;
t2339 = t892 * (t1503 - t932 + t2337 - t2328);
t2341 = t1503 * t2300;
t2345 = t2319 / 0.8e1;
t2346 = t2321 / 0.4e1;
t2347 = t2307 + t2310 + t2345 + t2346 - t2326 + t2330 - t2312 + t2333 - t2335 - t2314 + t2316 - t2339 - t2341;
t2357 = lapCoeff(i1,i2p2,i3,4);
t2358 = t2357 / 0.24e2;
t2360 = lapCoeff(i1,i2p2,i3,1);
t2361 = t2360 / 0.12e2;
t2365 = t1604 / 0.4e1;
t2366 = lapCoeff(i1p1,i2p2,i3,2);
t2367 = t2366 / 0.96e2;
t2368 = t1607 / 0.2e1;
t2379 = t121 * t333 / 0.32e2;
t2381 = t124 * t333 / 0.16e2;
t2382 = t2357 / 0.8e1;
t2383 = t2360 / 0.4e1;
t2384 = t1058 + t2382 + t1059 - t2383;
t2389 = 0.2e1 * t1607;
t2393 = t1058 + t1059;
t2398 = t333 + t2366;
t2400 = t1052 * t2398 / 0.32e2;
t2404 = t2366 / 0.4e1;
t2411 = t892 + t1506;
t2416 = t2182 + t2379 + t2184 + t2381 + t725 * t2384 / 0.4e1 + t961 * t2384 / 0.4e1 - t1607 * (t1741 - t1604 + t2328 + t2389) / 0.4e1 + t2 * t2393 / 0.4e1 + t63 * t2393 / 0.4e1 + t2202 + t2400 + t790 * (t1439 - t812 + t2203 + t2194) / 0.4e1 + t1055 * (t1604 - t1081 + t2404 - t2389) / 0.4e1 + t1442 * (t1741 + t1439 + t2194 + t2328) / 0.4e1 + t1439 * t2411 / 0.4e1 + t1604 * t2411 / 0.4e1;
t2418 = t2357 / 0.4e1;
t2419 = t2360 / 0.2e1;
t2421 = t895 * (t1095 + t2418 + t1096 - t2419);
t2423 = t761 * (t237 + t2140 + t333 + t2366) / 0.16e2;
t2425 = t1503 * (t2139 + t2365 + t2143 - t2368);
t2427 = t47 * (t1114 + t29);
t2428 = t1095 + t1096;
t2429 = t43 * t2428;
t2430 = t26 * t525 / 0.16e2;
t2431 = t2357 / 0.2e1;
t2432 = 0.4e1 * t892;
t2433 = 0.4e1 * t1506;
t2435 = t892 * (t1114 - t2431 + t29 + t2432 + t2360 + t2433);
t2436 = t1604 / 0.2e1;
t2438 = t1506 * (t2432 + t2151 - t2436 + t1442 + t2433 + t1607);
t2442 = t709 * (t1114 + t29 - t822);
t2444 = t938 * (t1114 + t29 + t1023);
t2446 = t274 * (t2151 - t822 + t1442);
t2448 = t470 * (t1023 - t2436 + t1607);
t2449 = t705 * t2428;
t2450 = t941 * t2428;
t2452 = t271 * (t2139 + t2143);
t2454 = t467 * (t2365 - t2368);
t2455 = t237 * t2411 / 0.2e1;
t2456 = t333 * t2411 / 0.2e1;
t2459 = t1604 / 0.8e1;
t2460 = t1607 / 0.4e1;
t2463 = 0.2e1 * t4 * (t2162 + t2459 + t2163 - t2460);
t2466 = 0.2e1 * t29 * (t24 + t1099 + t1741 + t2328);
t2469 = 0.2e1 * t241 * (t24 + t1099 - t812 + t1081);
t2471 = 0.2e1 * t24 * t2411;
t2472 = t2423 - t2421 - t2425 + t2463 + t2427 + t2429 + t1118 - t2430 + t2466 + t2469 + t2471 + t2435 + t2438;
t2474 = t2442 + t2444 + t2446 + t2448 + t2463 + t2449 - t2450 + t2452 - t2454 + t1118 - t2455 + t2456 + t2466 + t2469 + t2471;
t2484 = lapCoeff(i1p1,i2p2,i3,4);
t2485 = t2484 / 0.48e2;
t2487 = lapCoeff(i1p1,i2p2,i3,1);
t2488 = t2487 / 0.24e2;
t2489 = lapCoeff(i1p2,i2p1,i3,3);
t2490 = t2489 / 0.48e2;
t2492 = lapCoeff(i1p2,i2p1,i3,0);
t2493 = t2492 / 0.24e2;
t2503 = t2484 / 0.8e1;
t2504 = t2487 / 0.4e1;
t2506 = t761 * (t1232 + t2503 + t1233 - t2504);
t2509 = t26 * (t1232 + t1233);
t2512 = t892 * (t1503 - t932 + t2337 + t2328);
t2515 = t1506 * (t1872 + t1503 + t2328 + t2389);
t2517 = t1055 + t1607;
t2518 = t1503 * t2517;
t2522 = t2489 / 0.8e1;
t2523 = t2492 / 0.4e1;
t2525 = t333 * (t2286 + t2522 + t2287 - t2523);
t2528 = t237 * (t2286 + t2287);
t2531 = t470 * (t467 + t1269 + t1872 + t2389);
t2534 = t938 * (t467 + t1269 - t932 + t1283);
t2536 = t467 * t2517;
t2542 = 0.2e1 * t241 * (t1209 + t470 - t906);
t2545 = 0.2e1 * t29 * (t2271 - t906 + t1506);
t2546 = t1190 + t1191;
t2548 = 0.2e1 * t235 * t2546;
t2549 = t2252 + t2253;
t2551 = 0.2e1 * t24 * t2549;
t2553 = t4 * t2517;
t2554 = t2542 + t2545 - t2345 - t2346 + t2506 + t2548 - t2509 + t2335 + t2551 - t2553 + t2512 + t2515 + t2518;
t2556 = t2542 + t2545 - t1291 - t1292 + t2525 + t2548 + t2551 + t1278 - t2553 - t2528 + t2531 + t2534 + t2536;
t2559 = t2484 / 0.4e1;
t2560 = t2487 / 0.2e1;
t2564 = lapCoeff(i1p2,i2p2,i3,2);
t2568 = t2489 / 0.4e1;
t2569 = t2492 / 0.2e1;
t2587 = t2484 / 0.2e1;
t2588 = 0.4e1 * t1055;
t2589 = 0.4e1 * t1607;
t2593 = t2489 / 0.2e1;
t2597 = t277 / 0.64e2 - t1052 * (t1190 + t2559 + t1191 - t2560) / 0.4e1 + t961 * (t4 + t2227 + t1175 + t2564) / 0.64e2 - t1604 * (t2252 + t2568 + t2253 - t2569) / 0.4e1 - t124 * (t1209 + t470) / 0.4e1 - t1442 * (t2271 + t1506) / 0.4e1 - t121 * t2546 / 0.4e1 + t725 * t2265 / 0.64e2 + t63 * t1215 / 0.64e2 - t1439 * t2549 / 0.4e1 + t1055 * (t1209 - t2587 + t470 + t2588 + t2487 + t2589) / 0.4e1 + t1607 * (t2588 + t2271 - t2593 + t1506 + t2589 + t2492) / 0.4e1;
t2603 = lapCoeff(i1p2,i2p2,i3,4);
t2604 = t2603 / 0.8e1;
t2605 = lapCoeff(i1p2,i2p2,i3,1);
t2606 = t2605 / 0.4e1;
t2621 = t1978 + t2492;
t2627 = t938 * (t1301 + t1299 - t1023);
t2629 = t470 * (t2436 - t1023 + t1607);
t2631 = t941 * (t1308 + t1309);
t2632 = t2365 + t2368;
t2633 = t467 * t2632;
t2634 = t333 * t2621 / 0.2e1;
t2639 = 0.2e1 * t4 * (t2459 + t2460);
t2643 = t1506 * (t2436 + t1607);
t2644 = t761 * t2398 / 0.16e2;
t2645 = t1503 * t2632;
t2673 = t333 * (t2522 + t2523);
t2701 = lapCoeff(i1m2,i2p2,i3,3);
t2702 = t2701 / 0.8e1;
t2703 = lapCoeff(i1m2,i2p2,i3,0);
t2704 = t2703 / 0.4e1;
t2723 = t24 * t761;
t2724 = t2723 / 0.48e2;
t2725 = t29 * t761;
t2726 = t2725 / 0.24e2;
t2727 = lapCoeff(i1m1,i2p2,i3,3);
t2728 = t2727 / 0.48e2;
t2729 = lapCoeff(i1m1,i2p2,i3,0);
t2730 = t2729 / 0.24e2;
t2736 = t271 * t725;
t2737 = t2736 / 0.16e2;
t2738 = t274 * t725;
t2739 = t2738 / 0.8e1;
t2741 = t237 * (t1473 + t1475);
t2745 = t2727 / 0.8e1;
t2746 = t2729 / 0.4e1;
t2748 = t761 * (t2745 - t2746);
t2754 = t2727 / 0.2e1;
t2759 = t2727 / 0.4e1;
t2760 = t2729 / 0.2e1;
t2767 = t2736 / 0.8e1;
t2768 = t2738 / 0.4e1;
t2769 = t1188 / 0.8e1;
t2779 = t1506 * (t1677 + t790 + t1575);
t2780 = lapCoeff(i1,i2p2,i3,3);
t2781 = t2780 / 0.2e1;
t2782 = lapCoeff(i1,i2p2,i3,0);
t2784 = t892 * (t1575 - t2781 + t2782);
t2785 = t2723 / 0.4e1;
t2786 = t2725 / 0.2e1;
t2787 = t1580 + t1579;
t2788 = t1503 * t2787;
t2791 = 0.2e1 * t4 * (t1649 + t1650);
t2792 = t2780 / 0.4e1;
t2793 = t2782 / 0.2e1;
t2795 = t895 * (t2792 - t2793);
t2796 = t2218 + t2729;
t2797 = t761 * t2796 / 0.2e1;
t2802 = t2780 / 0.8e1;
t2803 = t2782 / 0.4e1;
t2823 = t274 * (t1677 + t790);
t2824 = t271 * t2787;
t2825 = t237 * t1621 / 0.16e2;
t2836 = t2780 / 0.24e2;
t2837 = t2782 / 0.12e2;
t2854 = t2360 / 0.6e1;
t2855 = t2782 / 0.6e1;
t2862 = lapCoeff(i1p1,i2p2,i3,3);
t2863 = t2862 / 0.48e2;
t2864 = lapCoeff(i1p1,i2p2,i3,0);
t2865 = t2864 / 0.24e2;
t2871 = t467 * t961;
t2872 = t2871 / 0.16e2;
t2873 = t470 * t961;
t2874 = t2873 / 0.8e1;
t2875 = t1758 + t1759;
t2876 = t237 * t2875;
t2878 = t333 * t2875;
t2882 = t2862 / 0.8e1;
t2883 = t2864 / 0.4e1;
t2885 = t761 * (t2745 + t2882 + t2746 - t2883);
t2888 = 0.2e1 * t29 * (t1790 + t892);
t2889 = t1778 + t1779;
t2891 = 0.2e1 * t24 * t2889;
t2893 = t4 * t1745 / 0.8e1;
t2897 = t892 * (t895 + t1741 + 0.2e1 * t2360 + 0.2e1 * t2782);
t2899 = t1506 * (t895 + t1741 - t1496 + t1751);
t2900 = t2360 + t2782;
t2901 = t895 * t2900;
t2910 = t2862 / 0.2e1;
t2917 = t2862 / 0.4e1;
t2918 = t2864 / 0.2e1;
t2927 = t2871 / 0.8e1;
t2928 = t2873 / 0.4e1;
t2937 = lapCoeff(i1p2,i2p2,i3,3);
t2938 = t2937 / 0.8e1;
t2939 = lapCoeff(i1p2,i2p2,i3,0);
t2940 = t2939 / 0.4e1;
t2955 = t2487 + t2864;
t2961 = t1506 * (t1921 + t1055 - t1575);
t2963 = t892 * (t2781 - t1575 + t2782);
t2964 = t1847 + t1848;
t2965 = t1503 * t2964;
t2967 = t895 * (t2792 + t2793);
t2968 = t761 * t2955 / 0.2e1;
t2973 = 0.2e1 * t4 * (t1949 + t1950);
t2977 = t470 * (t1921 + t1055);
t2978 = t467 * t2964;
t2979 = t333 * t1883 / 0.16e2;
t3007 = t333 * (t1999 + t2000);
t3012 = t761 * (t2882 + t2883);
t3028 = t962 / 0.8e1;
t3063 = t787 * t2140 / 0.32e2;
t3065 = t790 * t2140 / 0.16e2;
t3087 = t24 * (t761 + t2227) / 0.48e2;
t3089 = t29 * (t885 + t2228);
t3093 = t895 * t2227;
t3094 = t3093 / 0.16e2;
t3095 = t892 * t2227;
t3096 = t3095 / 0.8e1;
t3098 = t761 * (t2323 + t2324);
t3102 = t3093 / 0.8e1;
t3103 = t3095 / 0.4e1;
t3140 = t892 * (t2431 + t2360) + t895 * (t2418 + t2419) - t761 * (t2140 + t2366) / 0.16e2;
t3143 = t1052 * t2366 / 0.32e2;
t3145 = t1055 * t2366 / 0.16e2;
t3146 = t2382 + t2383;
t3159 = t761 * (t2503 + t2504);
sc(1,i1,i2) = t1 * t2 * t5 / 0.23040e5;
sc(2,i1,i2) = -dt6 * (t4 * (t11 - t14 + t2 * (t16 - t18) / 0.4e1) + t24 * t9 * t26 / 0.32e2 - t29 * t9 * t26 / 0.16e2) / 0.360e3 - t36 * (t9 + t2) / 0.1152e4;
sc(3,i1,i2) = t42 - dt6 * (t29 * (t46 - t49 + t55) - t24 * (t58 - t59 + t55 / 0.2e1) + t4 * (t44 * t63 / 0.64e2 - t8 * (t66 - t67) / 0.4e1 + t12 * (t71 - t52) / 0.4e1 + t2 * (t1 + t44) / 0.64e2)) / 0.360e3 + dt4 * (t4 * (t83 + t84 - t85 - t86) + t90 - t94) / 0.12e2;
sc(4,i1,i2) = dt4 * (t29 * (t100 + t101 - t103 - t104) - t24 * (t100 + t107 - t103 - t108) + t4 * (t9 + t2 + t111 + t63) / 0.96e2) / 0.12e2 - dt2 * (t117 - t118) + dt6 * (t4 * (t11 + t123 - t14 - t126 + t2 * t129 / 0.4e1 + t63 * t129 / 0.4e1) + t24 * t145 / 0.2e1 - t29 * t145) / 0.360e3;
sc(5,i1,i2) = dt6 * (t29 * (t46 - t49 + t157) - t24 * (t58 - t59 + t157 / 0.2e1) + t4 * (t2 * t44 / 0.64e2 - t121 * (t165 - t166) / 0.4e1 + t124 * (t170 - t154) / 0.4e1 + t63 * (t44 + t174) / 0.64e2)) / 0.360e3 - t42 - dt4 * (t4 * (t183 + t184 - t185 - t186) + t90 - t94) / 0.12e2;
sc(6,i1,i2) = -dt6 * (t4 * (t123 - t126 + t63 * (t194 - t196) / 0.4e1) + t202 * t111 / 0.32e2 - t205 * t111 / 0.16e2) / 0.360e3 - t36 * (t111 + t63) / 0.1152e4;
sc(7,i1,i2) = -t215 * t174 * dt6 / 0.23040e5;
sc(8,i1,i2) = -dt6 * (t4 * (t222 - t225 + t2 * (t227 - t229) / 0.4e1) + t235 * t220 * t237 / 0.32e2 - t220 * t237 * t241 / 0.16e2) / 0.360e3 - t36 * (t220 + t2) / 0.1152e4;
sc(9,i1,i2) = dt4 * (t251 - t252 + t254 - t256 + t4 * (t258 - t260 + t262 + t263 - t265)) / 0.12e2 - dt6 * (t241 * (t273 - t276 + t278 + t282) - t24 * (t287 - t290 + t294 / 0.2e1) - t235 * (t298 - t299 + t282 / 0.2e1) + t29 * (t303 + t278 - t304 + t294) + t4 * (t223 * (t307 - t259 + t308) + t12 * (t311 + t308 - t264) - t219 * (t314 - t315) - t8 * (t318 - t319) + t2 * (t17 + t228) / 0.2e1) / 0.4e1) / 0.360e3 + t329;
sc(10,i1,i2) = dt6 * (t241 * (t331 - t332 + t335 + t339 - t342 + t347 - t349) + t235 * (t352 - t339 / 0.2e1 + t342 / 0.2e1 + t349 / 0.2e1) + t4 * (t360 - t363 - t12 * (0.2e1 * t52 - t8 + t365 + 0.2e1 * t264) / 0.4e1 + t2 * (t227 + t371 + t229 - t373) / 0.4e1 + t378 + t63 * (t371 - t373) / 0.4e1 + t223 * (t8 - t365 - t382 + t383) / 0.4e1 + t8 * t387 / 0.4e1) + t29 * (t331 - t332 - t394 - t397 + t398 + t347 + t402 - t403) + t24 * (t394 + t397 - t398 - t402 + t403) / 0.2e1) / 0.360e3 - t411 + dt4 * (t29 * (t412 + t392 + t413 - t414) + t241 * (t271 - t274 + t4) / 0.12e2 - t24 * (t412 + t383 - t414) - t235 * (t271 - t274) / 0.24e2 + t4 * (t52 / 0.12e2 - t336 + t337 + t424 + t425 + t264 / 0.12e2)) / 0.12e2;
sc(11,i1,i2) = 0.3e1 / 0.20e2 * dt2 * (t24 - t29) + dt6 * (t29 * (t440 - t445 - t450 + t454 + t456 + t458 + t462 + t464) - t235 * (t298 + t469 - t299 - t472 + t476 / 0.2e1 + t478 / 0.2e1) - t24 * (t440 - t445 + t456 + t462 + t464) / 0.2e1 + t4 * (t223 * (t47 - t451 + t308) + t12 * (t311 + t308 + t264) + t361 * (t451 - t47 + t488) + t124 * (t491 + t488 - t437) - t219 * t448 - t358 * t448 - t8 * (t318 + t319) - t121 * (t498 - t499) + t2 * t463 / 0.2e1 + t63 * t463 / 0.2e1) / 0.4e1 + t241 * (t273 - t506 - t276 + t507 - t450 + t454 + t476 - t478 + t458)) / 0.360e3 - dt4 * (t24 * (t513 - t514 + t515 + t516) + t241 * (t519 - t520 - t521 + t522) + t526 - t29 * (t513 - t519 - t527 + t528 + t520 + t516) + t4 * (t262 + t531 + t263 + t532 + t265 - t533)) / 0.12e2;
sc(12,i1,i2) = t411 + dt4 * (t24 * (t412 + t383 + t414) - t241 * (t470 - t467 + t4) / 0.12e2 - t29 * (t412 + t392 + t413 + t414) + t235 * (t467 - t470) / 0.24e2 + t4 * (t548 - t154 / 0.12e2 - t550 + t551 + t553 - t437 / 0.12e2)) / 0.12e2 - dt6 * (t24 * (t562 + t564 - t566 - t568 + t570) / 0.2e1 - t29 * (t562 + t564 - t331 + t332 - t566 - t577 - t568 + t570) + t241 * (t331 - t332 - t581 - t582 + t584 + t577 + t586) + t235 * (t589 - t582 / 0.2e1 + t584 / 0.2e1 + t586 / 0.2e1) + t4 * (t596 + t598 - t124 * (0.2e1 * t154 - t121 + t600 + 0.2e1 * t437) / 0.4e1 + t63 * (t371 + t606 + t373 - t608) / 0.4e1 + t613 + t2 * (t371 + t373) / 0.4e1 + t361 * (t121 - t600 - t383 + t617) / 0.4e1 + t121 * t569 / 0.4e1)) / 0.360e3;
sc(13,i1,i2) = dt4 * (t252 - t251 + t630 + t632 + t4 * (t634 - t636 + t531 + t532 + t533)) / 0.12e2 + dt6 * (t235 * (t469 - t472 + t645 / 0.2e1) + t24 * (t650 + t652 + t654 / 0.2e1) - t4 * (t361 * (t633 - t658 + t488) + t124 * (t491 + t488 + t437) - t358 * (t663 - t664) - t121 * (t498 + t499) + t63 * (t195 + t607) / 0.2e1) / 0.4e1 + t241 * (t506 - t507 - t673 + t645) - t29 * (t676 + t673 + t677 + t654)) / 0.360e3 - t329;
sc(14,i1,i2) = dt6 * (t4 * (t685 + t687 + t63 * (t606 + t608) / 0.4e1) + t629 * t552 / 0.32e2 + t333 * t552 * t241 / 0.16e2) / 0.360e3 - t36 * (t63 + t552) / 0.1152e4;
sc(15,i1,i2) = t42 - dt6 * (t241 * (t708 - t711 + t717) - t235 * (t720 - t721 + t717 / 0.2e1) + t4 * (t706 * t725 / 0.64e2 + t2 * (t1 + t706) / 0.64e2 - t219 * (t731 - t732) / 0.4e1 + t223 * (t736 - t714) / 0.4e1)) / 0.360e3 + dt4 * (t4 * (t745 + t746 - t747 - t748) + t752 - t755) / 0.12e2;
sc(16,i1,i2) = dt6 * (t29 * (t760 + t763 - t764 - t766 + t770 + t775 - t778) + t24 * (t781 + t766 / 0.2e1 - t770 / 0.2e1 + t778 / 0.2e1) + t4 * (t789 - t792 + t2 * (t16 + t794 + t18 - t796) / 0.4e1 - t223 * (0.2e1 * t259 - t219 + 0.2e1 * t714 + t802) / 0.4e1 + t725 * (t794 - t796) / 0.4e1 + t810 + t12 * (t219 - t811 + t812 - t802) / 0.4e1 + t219 * t816 / 0.4e1) + t241 * (t760 - t764 - t824 - t826 + t830 + t831 - t832 + t775) + t235 * (t824 + t826 - t830 - t831 + t832) / 0.2e1) / 0.360e3 - t411 + dt4 * (t241 * (t840 - t841 + t822 + t413) - t235 * (t840 - t841 + t812) + t29 * (t285 + t4 - t288) / 0.12e2 - t24 * (t285 - t288) / 0.24e2 + t4 * (t259 / 0.12e2 - t767 + t851 + t852 + t714 / 0.12e2 + t768)) / 0.12e2;
sc(17,i1,i2) = t860 - dt4 * (t235 * (t861 - t862 + t753 + t863) + t24 * (t91 - t866 + t867 + t868) + t4 * (t83 + t871 + t85 + 0.23e2 / 0.24e2 * t12 - t873 + t745 + t874 + t747 + 0.23e2 / 0.24e2 * t223 - t876) + t241 * (t879 - t880 - t753 + t881 + t522) + t29 * (t884 - t91 + t881 + t885 - t886)) / 0.12e2 + dt6 * (t29 * t936 + t241 * t959 + t4 * t995 - t24 * (t997 - t998 + t901 / 0.2e1 - t904 / 0.2e1 + t920 / 0.2e1 + t1002 + t934 / 0.2e1 + t935 / 0.2e1) - t235 * (t1007 - t1008 - t946 / 0.2e1 + t950 / 0.2e1 + t1011 + t954 / 0.2e1 + t957 / 0.2e1 + t958 / 0.2e1)) / 0.360e3;
sc(18,i1,i2) = dt4 * (t241 * (t1021 - t1022 - t822 + t1023) - t29 * (t100 - t1021 - t1026 + t103 + 0.19e2 / 0.6e1 * t47 + t1022 + t1028 + 0.19e2 / 0.6e1 * t288) - t4 * (t851 - t1032 - t767 + t852 + t1033 + t1034 - t768 + t1035) + t1038 + t24 * (t100 + t1039 + t103 + 0.19e2 / 0.12e2 * t47 - t1041 + 0.19e2 / 0.12e2 * t288)) / 0.12e2 - dt2 * (t1048 - t1049) - dt6 * (t4 * t1093 + t29 * t1132 + t24 * (t1104 - t1098 - t1106 - t1113 + t1116 + t1119 + t1128 + t1131) / 0.2e1 + t241 * t1152 + t235 * (t1139 + t1137 + t1141 + t1143 - t1144 - t1145 - t1147 - t1149 + t1150 + t1151) / 0.2e1) / 0.360e3;
sc(19,i1,i2) = dt4 * (t4 * (t183 + t1160 + t185 + 0.23e2 / 0.24e2 * t124 - t1162 - t746 - t1164 - t748 + 0.23e2 / 0.24e2 * t361 + t1167) - t24 * (t866 - t91 - t867 + t868) - t235 * (t1172 - t1173 + t1174 + t1176) + t241 * (t1179 - t1180 + t521 + t881 - t1176) + t29 * (t884 - t91 + t881 + t885 + t886)) / 0.12e2 - dt6 * (t4 * t1230 - t24 * (t997 - t998 + t1235 / 0.2e1 + t1238 / 0.2e1 + t1002 + t1241 / 0.2e1 + t1244 / 0.2e1 + t1247 / 0.2e1) + t29 * t1263 - t235 * (t1266 + t1268 - t1271 / 0.2e1 + t1276 / 0.2e1 + t1279 + t1281 / 0.2e1 + t1285 / 0.2e1 + t1287 / 0.2e1) + t241 * t1293) / 0.360e3 - t860;
sc(20,i1,i2) = t411 + dt6 * (t241 * (t1303 + t1305 + t1306 + t1307 - t1311 - t1313 + t1315 + t1318) + t235 * (t1303 + t1305 - t1311 - t1313 + t1315) / 0.2e1 + t4 * (t1054 - t1057 + t63 * (t194 + t1323 + t196 - t1324) / 0.4e1 + t961 * (t1323 - t1324) / 0.4e1 + t1076 + t124 * (t358 - t1080 + t1081 + t1065) / 0.4e1 + t361 * (0.2e1 * t633 + t358 + t1065 + 0.2e1 * t1166) / 0.4e1 + t358 * t1314 / 0.4e1) + t29 * (t1306 + t1344 + t1307 + t1346 - t1347 + t1348 + t1318) + t24 * (t1351 - t1346 / 0.2e1 + t1347 / 0.2e1 - t1348 / 0.2e1)) / 0.360e3 - dt4 * (t241 * (t1360 - t1361 + t413 + t1023) + t235 * (t1360 - t1361 + t1081) + t29 * (t285 + t4 + t288) / 0.12e2 - t24 * (t285 + t288) / 0.24e2 + t4 * (t633 / 0.12e2 + t1032 - t1033 - t1034 + t1035 + t1166 / 0.12e2)) / 0.12e2;
sc(21,i1,i2) = dt4 * (t4 * (t874 + t1164 + t876 + t1167) + t1381 + t1383) / 0.12e2 - t42 - dt6 * (t235 * (t1388 + t1390 + t1392 / 0.2e1) + t241 * (t1396 + t1397 + t1392) - t4 * (t961 * t1175 / 0.64e2 - t361 * (t1226 + t1166) / 0.4e1 + t63 * (t174 + t1175) / 0.64e2 - t358 * (t1198 + t1199) / 0.4e1)) / 0.360e3;
sc(22,i1,i2) = dt4 * (t4 * (t220 + t1417 + t2 + t725) / 0.96e2 - t235 * (t1421 + t1422 - t1424 - t1425) + t241 * (t1421 + t1428 - t1424 - t1429)) / 0.12e2 - dt2 * (t1435 - t1436) + dt6 * (t4 * (t222 + t1441 - t225 - t1444 + t2 * t1447 / 0.4e1 + t725 * t1447 / 0.4e1) + t235 * t1463 / 0.2e1 - t241 * t1463) / 0.360e3;
sc(23,i1,i2) = 0.3e1 / 0.20e2 * dt2 * (t235 - t241) + dt6 * (t241 * (t1477 - t1482 + t1484 + t1486 - t1491 + t1495 + t1498 + t1500) - t24 * (t287 + t1505 - t290 - t1508 + t1512 / 0.2e1 + t1514 / 0.2e1) - t235 * (t1477 - t1482 + t1484 + t1498 + t1500) / 0.2e1 + t4 * (t223 * (t307 + t259 + t308) + t12 * (t308 - t1492 + t709) + t1442 * (t1524 - t1474 + t1525) + t790 * (t1492 + t1525 - t709) - t219 * (t314 + t315) - t1439 * (t1532 - t1533) - t8 * t1489 - t787 * t1489 + t2 * t1499 / 0.2e1 + t725 * t1499 / 0.2e1) / 0.4e1 + t29 * (t303 - t1542 - t304 + t1543 + t1486 - t1491 + t1512 - t1514 + t1495)) / 0.360e3 - dt4 * (t29 * (t1549 - t92 + t885 - t1550) + t235 * (t1553 - t1554 + t1555 + t1556) + t1560 + t4 * (t258 + t1561 + t260 - t1562 + t263 + t1563) - t241 * (t1553 - t1566 - t1549 + t1555 + t1567 + t1550)) / 0.12e2;
sc(24,i1,i2) = dt4 * (t29 * (t1574 - t392 + t1575 - t1576) - t4 * (t1579 - t1580 - t337 - t336 + t424 + t1581 + t425 + t1034) - t241 * (0.19e2 / 0.6e1 * t274 + t1421 - t1574 - t1585 + t1424 + 0.19e2 / 0.6e1 * t709 + t1576 + t1587) + t1590 + t235 * (0.19e2 / 0.12e2 * t274 + t1421 + t1592 + t1424 + 0.19e2 / 0.12e2 * t709 - t1594)) / 0.12e2 - dt2 * (t1600 - t1601) - dt6 * (t4 * t1645 + t241 * t1684 + t235 * (t1655 - t1648 - t1659 + t1666 - t1668 + t1671 + t1681 + t1683) / 0.2e1 + t29 * t1704 + t24 * (t1691 + t1689 + t1693 + t1695 - t1697 - t1699 - t1700 - t1701 + t1702 + t1703) / 0.2e1) / 0.360e3;
sc(25,i1,i2) = dt4 * (0.7e1 / 0.48e2 * t4 * t1712 - t24 * (t514 + t1714 + t515 - t1715) - t235 * (t1554 + t1718 + t1556 - t1719) + t29 * (t527 - t1722 + t528 + t1723 + t1724 + t1725) + t241 * (t1723 + t1566 - t1728 + t1567 + t1725 + t1729)) / 0.12e2 - 0.49e2 / 0.18e2 * dt2 * t1735 + dt6 * (t24 * (t1740 - t1743 + t1744 + t456 + t1746 + t1748 + t1750 + t1753 + t1754 + t1755) / 0.2e1 + t235 * (t1761 + t1762 - t1765 + t1484 + t1767 + t1769 + t1771 + t1773 + t1774 + t1775) / 0.2e1 - t29 * t1800 - t241 * t1802 + t4 * t1837 / 0.4e1) / 0.360e3 + 0.2e1;
sc(26,i1,i2) = dt2 * (t1600 + t1601) - dt4 * (t29 * (t1574 - t392 + t1575 + t1576) + t4 * (t548 + t1847 + t550 - t1848 + t551 + t852 + t553 + t1850) + t241 * (0.19e2 / 0.6e1 * t470 + t1428 + t1574 - t1855 + t1429 + t1576 + 0.19e2 / 0.6e1 * t938 + t1858) + t1590 + t235 * (0.19e2 / 0.12e2 * t470 - t1422 - t1855 - t1425 + 0.19e2 / 0.12e2 * t938 + t1858)) / 0.12e2 + dt6 * (t4 * t1906 - t235 * (t1909 - t1911 + t1915 + t1917 - t1918 + t1920 - t1925 - t1928) / 0.2e1 - t24 * (t1932 - t1934 + t1936 - t1938 + t1940 + t1942 + t1943 + t1944 - t1945 - t1946) / 0.2e1 + t241 * t1962 + t29 * t1964) / 0.360e3;
sc(27,i1,i2) = dt4 * (t29 * (t1549 - t92 + t885 + t1550) + t235 * (t1972 + t1718 + t1719 + t1973) + t1560 + t4 * (t636 + t1977 + t634 - t1979 - t532 - t1980) + t241 * (t1972 + t1549 + t1728 + t1550 + t1729 + t1973)) / 0.12e2 - dt6 * (t24 * (t650 + t1989 + t652 + t1991 + t1993 / 0.2e1 + t1995 / 0.2e1) + t241 * (t2002 + t2005 + t1767 - t2007 + t2010 + t2012 + t2016 + t2018) + t4 * (t361 * (t658 + t633 - t488) - t1607 * (t1978 - t2023 + t1816) + t124 * (t1796 - t488 + t938) - t1055 * (t1796 + t1816 + t938) + t358 * (t663 + t664) + t1604 * (t2032 - t2033) + t121 * t2008 + t1052 * t2008 - t63 * t2017 / 0.2e1 - t961 * t2017 / 0.2e1) / 0.4e1 + t235 * (t2002 + t1767 + t2012 + t2016 + t2018) / 0.2e1 + t29 * (t2044 - t676 - t677 + t2045 + t2005 - t2007 + t2010 - t1993 + t1995)) / 0.360e3 - 0.3e1 / 0.20e2 * dt2 * (t235 + t241);
sc(28,i1,i2) = dt2 * (t1435 + t1436) - dt4 * (t235 * (t1592 + t1855 + t1594 + t1858) - t4 * (t63 + t961 + t552 + t1849) / 0.96e2 + t241 * (t1585 + t1855 + t1587 + t1858)) / 0.12e2 + dt6 * (t235 * t2071 / 0.2e1 - t4 * (t685 + t2074 + t687 + t2076 + t63 * t2077 / 0.4e1 + t961 * t2077 / 0.4e1) + t241 * t2071) / 0.360e3;
sc(29,i1,i2) = dt6 * (t241 * (t708 - t711 + t2094) - t235 * (t720 - t721 + t2094 / 0.2e1) + t4 * (t706 * t2 / 0.64e2 + t725 * (t706 + t2102) / 0.64e2 - t1439 * (t2106 - t2107) / 0.4e1 + t1442 * (t2111 - t2091) / 0.4e1)) / 0.360e3 - t42 - dt4 * (t4 * (t2120 + t2121 - t2122 - t2123) + t752 - t755) / 0.12e2;
sc(30,i1,i2) = t411 - dt4 * (t241 * (t840 + t841 + t822 + t413) - t235 * (t840 + t841 + t812) + t29 * (t4 - t1503 + t1506) / 0.12e2 - t24 * (t1503 - t1506) / 0.24e2 + t4 * (t1474 / 0.12e2 - t2139 - t551 - t2141 + t2091 / 0.12e2 + t2143)) / 0.12e2 - dt6 * (t235 * (t2150 + t2153 - t2155 - t2157 + t2159) / 0.2e1 - t241 * (t2150 + t2153 - t760 + t764 - t2155 - t2157 + t2159 - t2166) + t29 * (t760 - t581 - t764 + t2170 - t2171 + t2166 + t2173) + t24 * (t589 + t2170 / 0.2e1 - t2171 / 0.2e1 + t2173 / 0.2e1) + t4 * (t2182 + t2184 + t725 * (t794 + t2186 + t796 - t2188) / 0.4e1 - t1442 * (0.2e1 * t1474 - t1439 + 0.2e1 * t2091 + t2194) / 0.4e1 + t2 * (t794 + t796) / 0.4e1 + t2202 + t790 * (t1439 - t812 + t2203 - t2194) / 0.4e1 + t1439 * t2158 / 0.4e1)) / 0.360e3;
sc(31,i1,i2) = dt4 * (t4 * (0.23e2 / 0.24e2 * t790 - t2217 - t86 - t84 + t2219 + t2120 + t2220 + t2122 + 0.23e2 / 0.24e2 * t1442 - t2222) - t24 * (t2225 + t2226 + t2228 - t2229) - t235 * (t862 + t861 - t753 - t863) + t241 * (t879 + t880 - t753 + t881 + t522) + t29 * (t92 - t2236 + t881 - t2228 + t2237)) / 0.12e2 - dt6 * (t4 * t2284 - t235 * (t1007 - t1008 + t2289 / 0.2e1 + t1011 + t2292 / 0.2e1 + t2295 / 0.2e1 + t2298 / 0.2e1 + t2301 / 0.2e1) + t241 * t2317 - t24 * (t2320 + t2322 + t2326 / 0.2e1 - t2330 / 0.2e1 + t2333 / 0.2e1 + t2336 + t2339 / 0.2e1 + t2341 / 0.2e1) + t29 * t2347) / 0.360e3 - t860;
sc(32,i1,i2) = dt2 * (t1048 + t1049) - dt4 * (t241 * (t1021 + t1022 - t822 + t1023) + t29 * (t101 + t1021 - t2358 + t104 + t1022 + 0.19e2 / 0.6e1 * t892 + t2361 + 0.19e2 / 0.6e1 * t1506) + t4 * (t2139 + t2365 + t551 + t2141 + t425 + t2367 + t2143 - t2368) + t1038 - t24 * (t107 + t2358 + t108 - 0.19e2 / 0.12e2 * t892 - t2361 - 0.19e2 / 0.12e2 * t1506)) / 0.12e2 + dt6 * (t4 * t2416 - t24 * (t2421 - t2423 + t2425 + t2427 + t2429 - t2430 - t2435 - t2438) / 0.2e1 - t235 * (t2442 - t2444 + t2446 - t2448 + t2449 + t2450 + t2452 + t2454 - t2455 - t2456) / 0.2e1 + t29 * t2472 + t241 * t2474) / 0.360e3;
sc(33,i1,i2) = dt4 * (t235 * (t1172 + t1173 + t1174 + t1176) + t24 * (t2225 + t2226 + t2228 + t2229) + t4 * (t184 + t2485 + t186 - 0.23e2 / 0.24e2 * t1055 - t2488 + t2121 + t2490 + t2123 - 0.23e2 / 0.24e2 * t1607 - t2493) + t241 * (t1180 + t1179 - t521 - t881 + t1176) + t29 * (t2236 - t92 - t881 + t2228 + t2237)) / 0.12e2 - dt6 * (t24 * (t2320 + t2322 + t2506 / 0.2e1 + t2509 / 0.2e1 + t2336 + t2512 / 0.2e1 + t2515 / 0.2e1 + t2518 / 0.2e1) + t235 * (t1266 + t1268 + t2525 / 0.2e1 + t1279 + t2528 / 0.2e1 + t2531 / 0.2e1 + t2534 / 0.2e1 + t2536 / 0.2e1) + t29 * t2554 + t241 * t2556 - t4 * t2597) / 0.360e3 + t860;
sc(34,i1,i2) = -t411 - dt6 * (t4 * (t2379 + t2381 + t961 * (t1323 + t2604 + t1324 - t2606) / 0.4e1 + t63 * (t1323 + t1324) / 0.4e1 + t2400 + t1055 * (t1604 - t1081 + t2404 + t2389) / 0.4e1 + t1607 * (0.2e1 * t1978 + t1604 + t2389 + 0.2e1 * t2492) / 0.4e1 + t1604 * t2621 / 0.4e1) - t235 * (t2627 + t2629 + t2631 + t2633 - t2634) / 0.2e1 - t241 * (t2627 + t2629 - t1306 - t1307 + t2631 + t2633 - t2634 - t2639) + t29 * (t1306 - t335 + t1307 - t2643 + t2644 - t2645 + t2639) + t24 * (t352 - t2643 / 0.2e1 + t2644 / 0.2e1 - t2645 / 0.2e1)) / 0.360e3 - dt4 * (t241 * (t1361 + t1360 - t413 - t1023) + t235 * (t1361 + t1360 - t1081) + t29 * (t1503 - t4 + t1506) / 0.12e2 + t24 * (t1503 + t1506) / 0.24e2 - t4 * (t1978 / 0.12e2 + t2365 + t425 + t2367 + t2368 + t2492 / 0.12e2)) / 0.12e2;
sc(35,i1,i2) = dt6 * (t235 * (t1388 + t1390 + t2673 / 0.2e1) + t241 * (t1396 + t1397 + t2673) - t4 * (t63 * t1175 / 0.64e2 - t1607 * (t2593 + t2492) / 0.4e1 + t961 * (t1175 + t2564) / 0.64e2 - t1604 * (t2568 + t2569) / 0.4e1)) / 0.360e3 + t42 - dt4 * (t4 * (t2220 + t2490 + t2222 + t2493) + t1381 + t1383) / 0.12e2;
sc(36,i1,i2) = -dt6 * (t4 * (t1441 - t1444 + t725 * (t2702 - t2704) / 0.4e1) + t235 * t1417 * t237 / 0.32e2 - t1417 * t237 * t241 / 0.16e2) / 0.360e3 - t36 * (t1417 + t725) / 0.1152e4;
sc(37,i1,i2) = dt4 * (t2724 + t2726 - t254 + t256 + t4 * (t1561 + t1562 - t2728 + t1563 + t2730)) / 0.12e2 + dt6 * (t235 * (t2737 + t2739 + t2741 / 0.2e1) + t24 * (t1505 - t1508 + t2748 / 0.2e1) - t4 * (t1442 * (t1524 + t1474 + t1525) + t790 * (t1525 - t2754 + t2729) - t1439 * (t1532 + t1533) - t787 * (t2759 - t2760) + t725 * (t2187 + t2703) / 0.2e1) / 0.4e1 - t241 * (t2767 + t2768 + t2769 + t2741) + t29 * (t1542 - t2769 - t1543 + t2748)) / 0.360e3 - t329;
sc(38,i1,i2) = t411 + dt6 * (t29 * (t2779 + t2784 + t2785 + t2786 - t2788 + t2791 - t2795 + t2797) + t24 * (t2779 + t2784 - t2788 - t2795 + t2797) / 0.2e1 + t4 * (t1606 - t1609 + t725 * (t2702 + t2802 + t2704 - t2803) / 0.4e1 + t1623 + t961 * (t2802 - t2803) / 0.4e1 + t790 * (t787 + t1610 + 0.2e1 * t2218 + 0.2e1 * t2729) / 0.4e1 + t1442 * (t787 + t1610 - t1635 + t1636) / 0.4e1 + t787 * t2796 / 0.4e1) + t241 * (t2785 + t2786 + t1344 + t2823 + t2824 + t2791 - t2825) + t235 * (t1351 - t2823 / 0.2e1 - t2824 / 0.2e1 + t2825 / 0.2e1)) / 0.360e3 - dt4 * (t29 * (t413 - t2836 + t1575 + t2837) + t241 * (t271 + t274 + t4) / 0.12e2 + t24 * (t1636 - t2836 + t2837) - t235 * (t271 + t274) / 0.24e2 + t4 * (t1580 + t1579 + t2218 / 0.12e2 - t1581 - t1034 + t2729 / 0.12e2)) / 0.12e2;
sc(39,i1,i2) = dt4 * (t24 * (t1714 + t1715 + t2854 + t2855) + t241 * (t519 + t520 - t521 + t522) + t526 + t29 * (t519 + t1722 + t520 + t1724 + t2854 + t2855) + t4 * (t2728 + t2863 - t1563 - t1980 + t2730 - t2865)) / 0.12e2 - dt6 * (t235 * (t2737 + t2872 + t2739 + t2874 + t2876 / 0.2e1 + t2878 / 0.2e1) + t29 * (t2885 + t2888 + t2891 + t1746 - t2893 + t2897 + t2899 + t2901) + t4 * (t1442 * (t1790 + t892 - t1525) - t1607 * (t1790 + t892 + t1816) + t790 * (t2754 - t1525 + t2729) - t1055 * (t1816 - t2910 + t2864) + t1439 * t2889 + t1604 * t2889 + t787 * (t2759 + t2760) + t1052 * (t2917 - t2918) - t725 * t2900 / 0.2e1 - t961 * t2900 / 0.2e1) / 0.4e1 + t24 * (t2885 + t1746 + t2897 + t2899 + t2901) / 0.2e1 + t241 * (t2927 - t2767 - t2768 + t2928 + t2888 + t2891 - t2876 + t2878 - t2893)) / 0.360e3 - 0.3e1 / 0.20e2 * dt2 * (t24 + t29);
sc(40,i1,i2) = -t411 - dt6 * (t4 * (t1869 + t1871 + t961 * (t2802 + t2938 + t2803 - t2940) / 0.4e1 + t1885 + t725 * (t2802 + t2803) / 0.4e1 + t1055 * (t1052 + t1872 + 0.2e1 * t2487 + 0.2e1 * t2864) / 0.4e1 + t1607 * (t1052 + t1872 - t1636 + t1897) / 0.4e1 + t1052 * t2955 / 0.4e1) - t24 * (t2961 + t2963 + t2965 + t2967 - t2968) / 0.2e1 - t29 * (t2961 + t2963 - t2785 - t2786 + t2965 - t2973 + t2967 - t2968) + t241 * (t2785 + t2786 - t763 - t2977 - t2978 + t2973 + t2979) + t235 * (t781 - t2977 / 0.2e1 - t2978 / 0.2e1 + t2979 / 0.2e1)) / 0.360e3 - dt4 * (t29 * (t2836 - t413 - t1575 + t2837) + t241 * (t467 + t470 - t4) / 0.12e2 + t24 * (t2836 - t1636 + t2837) + t235 * (t467 + t470) / 0.24e2 - t4 * (t1847 + t1848 + t2487 / 0.12e2 + t852 + t1850 + t2864 / 0.12e2)) / 0.12e2;
sc(41,i1,i2) = dt6 * (t235 * (t2872 + t2874 + t3007 / 0.2e1) + t24 * (t1989 + t1991 + t3012 / 0.2e1) + t4 * (t1607 * (t2023 + t1978 - t1816) + t1055 * (t2910 - t1816 + t2864) + t1604 * (t2032 + t2033) + t1052 * (t2917 + t2918) - t961 * (t2605 + t2939) / 0.2e1) / 0.4e1 + t241 * (t2927 + t2928 - t3028 + t3007) + t29 * (t2044 - t3028 + t2045 + t3012)) / 0.360e3 - dt4 * (t2724 + t2726 + t630 + t632 + t4 * (t1977 + t1979 + t2863 - t1980 + t2865)) / 0.12e2 + t329;
sc(42,i1,i2) = dt6 * (t4 * (t2074 + t2076 + t961 * (t2938 + t2940) / 0.4e1) + t629 * t1849 / 0.32e2 + t333 * t1849 * t241 / 0.16e2) / 0.360e3 - t36 * (t961 + t1849) / 0.1152e4;
sc(43,i1,i2) = -t2102 * t725 * t5 / 0.23040e5;
sc(44,i1,i2) = dt6 * (t4 * (t3063 + t3065 + t725 * (t2186 + t2188) / 0.4e1) + t24 * t2140 * t761 / 0.32e2 + t29 * t2140 * t761 / 0.16e2) / 0.360e3 - t36 * (t725 + t2140) / 0.1152e4;
sc(45,i1,i2) = dt4 * (t4 * (t871 + t2217 + t873 + t2219) + t3087 + t3089) / 0.12e2 - t42 - dt6 * (t24 * (t3094 + t3096 + t3098 / 0.2e1) + t29 * (t3102 + t3103 + t3098) - t4 * (t2227 * t961 / 0.64e2 - t790 * (t2275 + t2218) / 0.4e1 - t787 * (t2244 + t2245) / 0.4e1 + t725 * (t2102 + t2227) / 0.64e2)) / 0.360e3;
sc(46,i1,i2) = dt2 * (t117 + t118) - dt4 * (t24 * (t1039 + t2358 + t1041 + t2361) + t29 * (t1026 + t2358 + t1028 + t2361) - t4 * (t725 + t2140 + t961 + t2366) / 0.96e2) / 0.12e2 + dt6 * (t24 * t3140 / 0.2e1 - t4 * (t3063 + t3143 + t3065 + t3145 + t725 * t3146 / 0.4e1 + t961 * t3146 / 0.4e1) + t29 * t3140) / 0.360e3;
sc(47,i1,i2) = dt6 * (t24 * (t3094 + t3096 + t3159 / 0.2e1) + t29 * (t3102 + t3103 + t3159) - t4 * (t725 * t2227 / 0.64e2 - t1055 * (t2587 + t2487) / 0.4e1 - t1052 * (t2559 + t2560) / 0.4e1 + t961 * (t2227 + t2564) / 0.64e2)) / 0.360e3 + t42 - dt4 * (t4 * (t1160 + t2485 + t1162 + t2488) + t3087 + t3089) / 0.12e2;
sc(48,i1,i2) = dt6 * (t4 * (t3143 + t3145 + t961 * (t2604 + t2606) / 0.4e1) + t2723 * t2366 / 0.32e2 + t2725 * t2366 / 0.16e2) / 0.360e3 - t36 * (t961 + t2366) / 0.1152e4;
sc(49,i1,i2) = t962 * t2564 * dt6 / 0.23040e5;
             ! #Include "../include/defineStencilVariables2dOrder2Curvilinear.h"
          end do
          end do
          end do
         ! ----- END LOOPS ----
     end if
   if( gridIsImplicit.eq.0 )then 
     ! ------- EXPLICIT update the solution ---------
       if( orderInTime.eq.2 .and. orderOfAccuracy.gt.2 )then
         stop 2222
         ! if( addForcing.eq.0 )then
         !   updateWaveOpt(2,6,2,curvilinear,NOFORCING)
         ! else 
         !   updateWaveOpt(2,6,2,curvilinear,FORCING)
         ! end if
       else
         if( addForcing.eq.0 )then
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=2 order=6 orderInTime=6, grid=curvilinear... t=",e10.2)') t
             end if
             ! --- TAYLOR TIME-STEPPING --- 
             m=0 ! component number 
             ec = 0 ! component number
             if( forcingOption.eq.helmholtzForcing )then
               coswt = cos(omega*t)
             end if 
             fv(m)=0.
             ! >>>> NOTE: NO-MASK IS FASTER for square128
             ! beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
! Stencil: nd=2, orderOfAccuracy=6, gridType=Curvilinear
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ sc(  1,i1,i2)*u(i1-3,i2-3,i3,m) + sc(  2,i1,i2)*u(i1-2,i2-3,i3,m) + sc(  3,i1,i2)*u(i1-1,i2-3,i3,m) + sc(  4,i1,i2)*u(i1+0,i2-3,i3,m) + sc(  5,i1,i2)*u(i1+1,i2-3,i3,m) + sc(  6,i1,i2)*u(i1+2,i2-3,i3,m) + sc(  7,i1,i2)*u(i1+3,i2-3,i3,m)+ sc(  8,i1,i2)*u(i1-3,i2-2,i3,m) + sc(  9,i1,i2)*u(i1-2,i2-2,i3,m) + sc( 10,i1,i2)*u(i1-1,i2-2,i3,m) + sc( 11,i1,i2)*u(i1+0,i2-2,i3,m) + sc( 12,i1,i2)*u(i1+1,i2-2,i3,m) + sc( 13,i1,i2)*u(i1+2,i2-2,i3,m) + sc( 14,i1,i2)*u(i1+3,i2-2,i3,m)+ sc( 15,i1,i2)*u(i1-3,i2-1,i3,m) + sc( 16,i1,i2)*u(i1-2,i2-1,i3,m) + sc( 17,i1,i2)*u(i1-1,i2-1,i3,m) + sc( 18,i1,i2)*u(i1+0,i2-1,i3,m) + sc( 19,i1,i2)*u(i1+1,i2-1,i3,m) + sc( 20,i1,i2)*u(i1+2,i2-1,i3,m) + sc( 21,i1,i2)*u(i1+3,i2-1,i3,m)+ sc( 22,i1,i2)*u(i1-3,i2+0,i3,m) + sc( 23,i1,i2)*u(i1-2,i2+0,i3,m) + sc( 24,i1,i2)*u(i1-1,i2+0,i3,m) + sc( 25,i1,i2)*u(i1+0,i2+0,i3,m) + sc( 26,i1,i2)*u(i1+1,i2+0,i3,m) + sc( 27,i1,i2)*u(i1+2,i2+0,i3,m) + sc( 28,i1,i2)*u(i1+3,i2+0,i3,m)+ sc( 29,i1,i2)*u(i1-3,i2+1,i3,m) + sc( 30,i1,i2)*u(i1-2,i2+1,i3,m) + sc( 31,i1,i2)*u(i1-1,i2+1,i3,m) + sc( 32,i1,i2)*u(i1+0,i2+1,i3,m) + sc( 33,i1,i2)*u(i1+1,i2+1,i3,m) + sc( 34,i1,i2)*u(i1+2,i2+1,i3,m) + sc( 35,i1,i2)*u(i1+3,i2+1,i3,m)+ sc( 36,i1,i2)*u(i1-3,i2+2,i3,m) + sc( 37,i1,i2)*u(i1-2,i2+2,i3,m) + sc( 38,i1,i2)*u(i1-1,i2+2,i3,m) + sc( 39,i1,i2)*u(i1+0,i2+2,i3,m) + sc( 40,i1,i2)*u(i1+1,i2+2,i3,m) + sc( 41,i1,i2)*u(i1+2,i2+2,i3,m) + sc( 42,i1,i2)*u(i1+3,i2+2,i3,m)+ sc( 43,i1,i2)*u(i1-3,i2+3,i3,m) + sc( 44,i1,i2)*u(i1-2,i2+3,i3,m) + sc( 45,i1,i2)*u(i1-1,i2+3,i3,m) + sc( 46,i1,i2)*u(i1+0,i2+3,i3,m) + sc( 47,i1,i2)*u(i1+1,i2+3,i3,m) + sc( 48,i1,i2)*u(i1+2,i2+3,i3,m) + sc( 49,i1,i2)*u(i1+3,i2+3,i3,m)
               end do
               end do
               end do
             ! endLoopsMask()
         else
             if( ( .true. .or. debug.gt.3) .and. t.lt.2*dt )then
               write(*,'("advWaveStencil: ADVANCE dim=2 order=6 orderInTime=6, grid=curvilinear... t=",e10.2)') t
             end if
             ! --- TAYLOR TIME-STEPPING --- 
             m=0 ! component number 
             ec = 0 ! component number
             if( forcingOption.eq.helmholtzForcing )then
               coswt = cos(omega*t)
             end if 
             fv(m)=0.
             ! >>>> NOTE: NO-MASK IS FASTER for square128
             ! beginLoopsMask(i1,i2,i3,n1a,n1b,n2a,n2b,n3a,n3b)
               do i3=n3a,n3b
               do i2=n2a,n2b
               do i1=n1a,n1b
                   if( forcingOption.eq.twilightZoneForcing )then
                         call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                         call ogDeriv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evtt(m) )
                         call ogDeriv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxx(m) )
                         call ogDeriv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evyy(m) )
                       fv(m) = evtt(m) - csq*( evxx(m) + evyy(m) )
                        ! Correct forcing for fourth-order ME in2D
                          call ogDeriv(ep, 4,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evtttt(m) )
                          call ogDeriv(ep, 0,4,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxxxx(m) )
                          call ogDeriv(ep, 0,2,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxxyy(m) )
                          call ogDeriv(ep, 0,0,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evyyyy(m) )
                        fv(m) = fv(m) + (dtSq/12.)*evtttt(m) - (cdtsq12/dtSq)*( evxxxx(m) + 2.*evxxyy(m) + evyyyy(m) )
                        ! Correct forcing for sixth-order ME in 2D
                          call ogDeriv(ep, 6,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evtttttt(m) )
                          call ogDeriv(ep, 0,6,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxxxxxx(m) )
                          call ogDeriv(ep, 0,4,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxxxxyy(m) )
                          call ogDeriv(ep, 0,2,4,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxxyyyy(m) )
                          call ogDeriv(ep, 0,0,6,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evyyyyyy(m) )
                        fv(m) = fv(m) + (dtSq**2/360.)*evtttttt(m) - (cdtPow6By360/dtSq)*( evxxxxxx(m) + evyyyyyy(m) + 3.*(evxxxxyy(m) + evxxyyyy(m) )  )
                  else if( forcingOption.eq.helmholtzForcing )then
                     ! forcing for solving the Helmholtz equation   
                     ! NOTE: change sign of forcing since for Helholtz we want to solve
                     !      ( omega^2 I + c^2 Delta) w = f 
                     ! fv(m) = -f(i1,i2,i3,0)*coswt  
                     fv(m)=0.
                     do freq=0,numberOfFrequencies-1 
                       omega = frequencyArray(freq)
                       coswt = cosFreqt(freq)    
                         ! Add corrections for 4th order modified equation 
                         !  fv = f + (dt^2/12)*( c^2 Delta(u) + ftt )
                         write(*,*) 'fix me'
                         stop 4444
                             !fv(m) = fv(m) -( f(i1,i2,i3,freq) + cdtSqBy12*( cSq*(fxx22(i1,i2,i3,freq) + fyy22(i1,i2,i3,freq)) - omega*omega*f(i1,i2,i3,freq)) )*coswt 
                     end do ! do freq  
                  else if( addForcing.ne.0 )then  
                     fv(m) = f(i1,i2,i3,0)
                  end if
! Stencil: nd=2, orderOfAccuracy=6, gridType=Curvilinear
un(i1,i2,i3,m)=  - um(i1,i2,i3,m)+ sc(  1,i1,i2)*u(i1-3,i2-3,i3,m) + sc(  2,i1,i2)*u(i1-2,i2-3,i3,m) + sc(  3,i1,i2)*u(i1-1,i2-3,i3,m) + sc(  4,i1,i2)*u(i1+0,i2-3,i3,m) + sc(  5,i1,i2)*u(i1+1,i2-3,i3,m) + sc(  6,i1,i2)*u(i1+2,i2-3,i3,m) + sc(  7,i1,i2)*u(i1+3,i2-3,i3,m)+ sc(  8,i1,i2)*u(i1-3,i2-2,i3,m) + sc(  9,i1,i2)*u(i1-2,i2-2,i3,m) + sc( 10,i1,i2)*u(i1-1,i2-2,i3,m) + sc( 11,i1,i2)*u(i1+0,i2-2,i3,m) + sc( 12,i1,i2)*u(i1+1,i2-2,i3,m) + sc( 13,i1,i2)*u(i1+2,i2-2,i3,m) + sc( 14,i1,i2)*u(i1+3,i2-2,i3,m)+ sc( 15,i1,i2)*u(i1-3,i2-1,i3,m) + sc( 16,i1,i2)*u(i1-2,i2-1,i3,m) + sc( 17,i1,i2)*u(i1-1,i2-1,i3,m) + sc( 18,i1,i2)*u(i1+0,i2-1,i3,m) + sc( 19,i1,i2)*u(i1+1,i2-1,i3,m) + sc( 20,i1,i2)*u(i1+2,i2-1,i3,m) + sc( 21,i1,i2)*u(i1+3,i2-1,i3,m)+ sc( 22,i1,i2)*u(i1-3,i2+0,i3,m) + sc( 23,i1,i2)*u(i1-2,i2+0,i3,m) + sc( 24,i1,i2)*u(i1-1,i2+0,i3,m) + sc( 25,i1,i2)*u(i1+0,i2+0,i3,m) + sc( 26,i1,i2)*u(i1+1,i2+0,i3,m) + sc( 27,i1,i2)*u(i1+2,i2+0,i3,m) + sc( 28,i1,i2)*u(i1+3,i2+0,i3,m)+ sc( 29,i1,i2)*u(i1-3,i2+1,i3,m) + sc( 30,i1,i2)*u(i1-2,i2+1,i3,m) + sc( 31,i1,i2)*u(i1-1,i2+1,i3,m) + sc( 32,i1,i2)*u(i1+0,i2+1,i3,m) + sc( 33,i1,i2)*u(i1+1,i2+1,i3,m) + sc( 34,i1,i2)*u(i1+2,i2+1,i3,m) + sc( 35,i1,i2)*u(i1+3,i2+1,i3,m)+ sc( 36,i1,i2)*u(i1-3,i2+2,i3,m) + sc( 37,i1,i2)*u(i1-2,i2+2,i3,m) + sc( 38,i1,i2)*u(i1-1,i2+2,i3,m) + sc( 39,i1,i2)*u(i1+0,i2+2,i3,m) + sc( 40,i1,i2)*u(i1+1,i2+2,i3,m) + sc( 41,i1,i2)*u(i1+2,i2+2,i3,m) + sc( 42,i1,i2)*u(i1+3,i2+2,i3,m)+ sc( 43,i1,i2)*u(i1-3,i2+3,i3,m) + sc( 44,i1,i2)*u(i1-2,i2+3,i3,m) + sc( 45,i1,i2)*u(i1-1,i2+3,i3,m) + sc( 46,i1,i2)*u(i1+0,i2+3,i3,m) + sc( 47,i1,i2)*u(i1+1,i2+3,i3,m) + sc( 48,i1,i2)*u(i1+2,i2+3,i3,m) + sc( 49,i1,i2)*u(i1+3,i2+3,i3,m)+dtSq*fv(m)
               end do
               end do
               end do
             ! endLoopsMask()
         end if
       end if 
   else
     ! --- IMPLICIT: Fill in RHS to implicit time-stepping -----
     stop 1111
   end if
   return
   end
