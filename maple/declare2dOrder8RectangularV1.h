! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=2, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare2dOrder8Rectangular()
real cxx, cyy, czz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d600(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d420(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d240(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real d060(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), cx0, cx1, cx2, cx3, cy0, cy1, cy2, cy3, cxx0, cxx1
real cxx2, cxx3, cyy0, cyy1, cyy2, cyy3, cxxx0, cxxx1, cxxx2, cxxx3, cyyy0
real cyyy1, cyyy2, cyyy3, cxxxx0, cxxxx1, cxxxx2, cxxxx3, cyyyy0, cyyyy1, cyyyy2, cyyyy3
real cxxxxx0, cxxxxx1, cxxxxx2, cxxxxx3, cyyyyy0, cyyyyy1, cyyyyy2, cyyyyy3, cxxxxxx0, cxxxxxx1, cxxxxxx2
real cxxxxxx3, cyyyyyy0, cyyyyyy1, cyyyyyy2, cyyyyyy3, cxxxxxxx0, cxxxxxxx1, cxxxxxxx2, cxxxxxxx3, cyyyyyyy0, cyyyyyyy1
real cyyyyyyy2, cyyyyyyy3, cxxxxxxxx0, cxxxxxxxx1, cxxxxxxxx2, cxxxxxxxx3, cyyyyyyyy0, cyyyyyyyy1, cyyyyyyyy2, cyyyyyyyy3, d200i
real d400i, d600i, d800i, d020i, d220i, d420i, d620i, d040i, d240i, d440i, d060i
real d260i, d080i, uxx, uyy, uxxxx, uxxyy, uyyyy, uxxxxxx, uxxxxyy, uxxyyyy, uyyyyyy
real uxxxxxxxx, uxxxxxxyy, uxxxxyyyy, uxxyyyyyy, uyyyyyyyy
#endMacro