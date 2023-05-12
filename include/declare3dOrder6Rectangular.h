! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=3, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare3dOrder6Rectangular()
real cx0, cx1, cx2, cy0, cy1, cy2, cz0, cz1, cz2, cxx0, cxx1
real cxx2, cyy0, cyy1, cyy2, czz0, czz1, czz2, cxxx0, cxxx1, cxxx2, cyyy0
real cyyy1, cyyy2, czzz0, czzz1, czzz2, cxxxx0, cxxxx1, cxxxx2, cyyyy0, cyyyy1, cyyyy2
real czzzz0, czzzz1, czzzz2, cxx, cyy, czz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d004(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap6h, lap4hSq, lap2hCubed, d600i
real d060i, d006i, lap2hSq200i, lap2hSq020i, lap2hSq002i, lap4h200i, lap4h020i, lap4h002i, lap2h400i, lap2h040i, lap2h004i

#endMacro