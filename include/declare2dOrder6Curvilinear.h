! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=2, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare2dOrder6Curvilinear()
real cr0, cr1, cr2, cs0, cs1, cs2, crr0, crr1, crr2, css0, css1
real css2, crrr0, crrr1, crrr2, csss0, csss1, csss2, crrrr0, crrrr1, crrrr2, cssss0
real cssss1, cssss2, dr1, dr2, dr3, dr1i, dr2i, dr3i, rx, ry, sx
real sy, diffOrder1, diffOrder2, diffOrder3, rxr, rxs, ryr, rys, sxr, sxs, syr
real sys, rxx, ryy, sxx, syy, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d100i, d010i, d110i
real d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d300i, d030i, d310i, d130i, lap2h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real lap2h100i, lap2h010i, lap2h110i, lap6h, lap4hSq, lap2hCubed, d600i, d060i, d500i, d050i, d510i
real d150i, d330i, lap2hSq200i, lap2hSq020i, lap2hSq100i, lap2hSq010i, lap2hSq110i, lap4h200i, lap4h020i, lap2h400i, lap2h040i
real lap4h100i, lap4h010i, lap4h110i, lap2h300i, lap2h030i, lap2h310i, lap2h130i
#endMacro