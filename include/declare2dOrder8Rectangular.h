! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=2, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare2dOrder8Rectangular()
real cx0, cx1, cx2, cx3, cy0, cy1, cy2, cy3, cxx0, cxx1, cxx2
real cxx3, cyy0, cyy1, cyy2, cyy3, cxxx0, cxxx1, cxxx2, cxxx3, cyyy0, cyyy1
real cyyy2, cyyy3, cxxxx0, cxxxx1, cxxxx2, cxxxx3, cyyyy0, cyyyy1, cyyyy2, cyyyy3, cxxxxx0
real cxxxxx1, cxxxxx2, cxxxxx3, cyyyyy0, cyyyyy1, cyyyyy2, cyyyyy3, cxxxxxx0, cxxxxxx1, cxxxxxx2, cxxxxxx3
real cyyyyyy0, cyyyyyy1, cyyyyyy2, cyyyyyy3, cxx, cyy, czz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap6h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hCubed(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d600(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d060(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real lap2hSq020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d800i, d080i, lap8h, lap2hCubed200i, lap2hCubed020i, lap2h4p
real lap6h200i, lap6h020i, lap4h400i, lap4h040i, lap2h600i, lap2h060i, lap6hSq, lap4hSq200i, lap4hSq020i, lap2hSq400i, lap2hSq040i
real lap4hCubed
#endMacro