! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=3, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare3dOrder8Rectangular()
real cxx, cyy, czz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d202(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d022(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real d004(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d600(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d420(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d240(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d060(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d402(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d222(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d042(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d204(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d024(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d006(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real cx0, cx1, cx2, cx3, cy0, cy1, cy2, cy3, cz0, cz1, cz2
real cz3, cxx0, cxx1, cxx2, cxx3, cyy0, cyy1, cyy2, cyy3, czz0, czz1
real czz2, czz3, cxxx0, cxxx1, cxxx2, cxxx3, cyyy0, cyyy1, cyyy2, cyyy3, czzz0
real czzz1, czzz2, czzz3, cxxxx0, cxxxx1, cxxxx2, cxxxx3, cyyyy0, cyyyy1, cyyyy2, cyyyy3
real czzzz0, czzzz1, czzzz2, czzzz3, cxxxxx0, cxxxxx1, cxxxxx2, cxxxxx3, cyyyyy0, cyyyyy1, cyyyyy2
real cyyyyy3, czzzzz0, czzzzz1, czzzzz2, czzzzz3, cxxxxxx0, cxxxxxx1, cxxxxxx2, cxxxxxx3, cyyyyyy0, cyyyyyy1
real cyyyyyy2, cyyyyyy3, czzzzzz0, czzzzzz1, czzzzzz2, czzzzzz3, cxxxxxxx0, cxxxxxxx1, cxxxxxxx2, cxxxxxxx3, cyyyyyyy0
real cyyyyyyy1, cyyyyyyy2, cyyyyyyy3, czzzzzzz0, czzzzzzz1, czzzzzzz2, czzzzzzz3, cxxxxxxxx0, cxxxxxxxx1, cxxxxxxxx2, cxxxxxxxx3
real cyyyyyyyy0, cyyyyyyyy1, cyyyyyyyy2, cyyyyyyyy3, czzzzzzzz0, czzzzzzzz1, czzzzzzzz2, czzzzzzzz3, d200i, d400i, d600i
real d800i, d020i, d220i, d420i, d620i, d040i, d240i, d440i, d060i, d260i, d080i
real d002i, d202i, d402i, d602i, d022i, d222i, d422i, d042i, d242i, d062i, d004i
real d204i, d404i, d024i, d224i, d044i, d006i, d206i, d026i, d008i, uxx, uyy
real uzz, uxxxx, uxxyy, uyyyy, uxxzz, uyyzz, uzzzz, uxxxxxx, uxxxxyy, uxxyyyy, uyyyyyy
real uxxxxzz, uxxyyzz, uyyyyzz, uxxzzzz, uyyzzzz, uzzzzzz, uxxxxxxxx, uxxxxxxyy, uxxxxyyyy, uxxyyyyyy, uyyyyyyyy
real uxxxxxxzz, uxxxxyyzz, uxxyyyyzz, uyyyyyyzz, uxxxxzzzz, uxxyyzzzz, uyyyyzzzz, uxxzzzzzz, uyyzzzzzz, uzzzzzzzz
#endMacro