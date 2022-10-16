! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=3, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare3dOrder8Rectangular()
real cx0, cx1, cx2, cx3, cy0, cy1, cy2, cy3, cz0, cz1, cz2
real cz3, cxx0, cxx1, cxx2, cxx3, cyy0, cyy1, cyy2, cyy3, czz0, czz1
real czz2, czz3, cxxx0, cxxx1, cxxx2, cxxx3, cyyy0, cyyy1, cyyy2, cyyy3, czzz0
real czzz1, czzz2, czzz3, cxxxx0, cxxxx1, cxxxx2, cxxxx3, cyyyy0, cyyyy1, cyyyy2, cyyyy3
real czzzz0, czzzz1, czzzz2, czzzz3, cxxxxx0, cxxxxx1, cxxxxx2, cxxxxx3, cyyyyy0, cyyyyy1, cyyyyy2
real cyyyyy3, czzzzz0, czzzzz1, czzzzz2, czzzzz3, cxxxxxx0, cxxxxxx1, cxxxxxx2, cxxxxxx3, cyyyyyy0, cyyyyyy1
real cyyyyyy2, cyyyyyy3, czzzzzz0, czzzzzz1, czzzzzz2, czzzzzz3, cxx, cyy, czz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real d002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d004(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap6h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real lap4hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hCubed(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d600(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d060(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d006(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real lap2h400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h004(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d800i, d080i, d008i, lap8h, lap2hCubed200i, lap2hCubed020i, lap2hCubed002i, lap2h4p
real lap6h200i, lap6h020i, lap6h002i, lap4h400i, lap4h040i, lap4h004i, lap2h600i, lap2h060i, lap2h006i, lap6hSq, lap4hSq200i
real lap4hSq020i, lap4hSq002i, lap2hSq400i, lap2hSq040i, lap2hSq004i, lap4hCubed
#endMacro