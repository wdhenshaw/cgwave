! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=2, gridType=Rectangular
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare2dOrder4Rectangular()
real cx0, cx1, cy0, cy1, cxx0, cxx1, cyy0, cyy1, cxx, cyy, czz
real d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d400i, d040i, lap4h, lap2hSq, lap2h200i, lap2h020i
#endMacro