! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=2, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare2dOrder6Rectangular()
real cx0, cx1, cx2, cy0, cy1, cy2, cxx0, cxx1, cxx2, cyy0, cyy1
real cyy2, cxxx0, cxxx1, cxxx2, cyyy0, cyyy1, cyyy2, cxxxx0, cxxxx1, cxxxx2, cyyyy0
real cyyyy1, cyyyy2, cxx, cyy, czz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real lap2hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap6h, lap4hSq, lap2hCubed, d600i, d060i, lap2hSq200i, lap2hSq020i, lap4h200i
real lap4h020i, lap2h400i, lap2h040i
#endMacro