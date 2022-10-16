! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=3, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare3dOrder6Rectangular()
real cxx, cyy, czz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d202(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d022(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real d004(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), cx0, cx1, cx2, cy0, cy1, cy2, cz0, cz1, cz2, cxx0
real cxx1, cxx2, cyy0, cyy1, cyy2, czz0, czz1, czz2, cxxx0, cxxx1, cxxx2
real cyyy0, cyyy1, cyyy2, czzz0, czzz1, czzz2, cxxxx0, cxxxx1, cxxxx2, cyyyy0, cyyyy1
real cyyyy2, czzzz0, czzzz1, czzzz2, cxxxxx0, cxxxxx1, cxxxxx2, cyyyyy0, cyyyyy1, cyyyyy2, czzzzz0
real czzzzz1, czzzzz2, cxxxxxx0, cxxxxxx1, cxxxxxx2, cyyyyyy0, cyyyyyy1, cyyyyyy2, czzzzzz0, czzzzzz1, czzzzzz2
real d200i, d400i, d600i, d020i, d220i, d420i, d040i, d240i, d060i, d002i, d202i
real d402i, d022i, d222i, d042i, d004i, d204i, d024i, d006i, uxx, uyy, uzz
real uxxxx, uxxyy, uyyyy, uxxzz, uyyzz, uzzzz, uxxxxxx, uxxxxyy, uxxyyyy, uyyyyyy, uxxxxzz
real uxxyyzz, uyyyyzz, uxxzzzz, uyyzzzz, uzzzzzz
#endMacro