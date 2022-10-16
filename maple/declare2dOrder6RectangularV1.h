! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=2, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare2dOrder6Rectangular()
real cxx, cyy, czz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), cx0, cx1, cx2
real cy0, cy1, cy2, cxx0, cxx1, cxx2, cyy0, cyy1, cyy2, cxxx0, cxxx1
real cxxx2, cyyy0, cyyy1, cyyy2, cxxxx0, cxxxx1, cxxxx2, cyyyy0, cyyyy1, cyyyy2, cxxxxx0
real cxxxxx1, cxxxxx2, cyyyyy0, cyyyyy1, cyyyyy2, cxxxxxx0, cxxxxxx1, cxxxxxx2, cyyyyyy0, cyyyyyy1, cyyyyyy2
real d200i, d400i, d600i, d020i, d220i, d420i, d040i, d240i, d060i, uxx, uyy
real uxxxx, uxxyy, uyyyy, uxxxxxx, uxxxxyy, uxxyyyy, uyyyyyy
#endMacro