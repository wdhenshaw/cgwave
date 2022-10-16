! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=2, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare2dOrder4Curvilinear()
real cr0, cr1, cs0, cs1, crr0, crr1, css0, css1, dr1, dr2, dr3
real dr1i, dr2i, dr3i, rx, ry, sx, sy, diffOrder1, diffOrder2, diffOrder3, rxr
real rxs, ryr, rys, sxr, sxs, syr, sys, rxx, ryy, sxx, syy
real d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d100i, d010i, d110i, d400i, d040i, lap4h, lap2hSq, d300i
real d030i, d310i, d130i, lap2h200i, lap2h020i, lap2h100i, lap2h010i, lap2h110i
#endMacro