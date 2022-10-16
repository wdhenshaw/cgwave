! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=2, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare2dOrder8Curvilinear()
real cr0, cr1, cr2, cr3, cs0, cs1, cs2, cs3, crr0, crr1, crr2
real crr3, css0, css1, css2, css3, crrr0, crrr1, crrr2, crrr3, csss0, csss1
real csss2, csss3, crrrr0, crrrr1, crrrr2, crrrr3, cssss0, cssss1, cssss2, cssss3, crrrrr0
real crrrrr1, crrrrr2, crrrrr3, csssss0, csssss1, csssss2, csssss3, crrrrrr0, crrrrrr1, crrrrrr2, crrrrrr3
real cssssss0, cssssss1, cssssss2, cssssss3, dr1, dr2, dr3, dr1i, dr2i, dr3i, rx
real ry, sx, sy, diffOrder1, diffOrder2, diffOrder3, rxr, rxs, ryr, rys, sxr
real sxs, syr, sys, rxx, ryy, sxx, syy, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d100i
real d010i, d110i, d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d300i, d030i, d310i, d130i
real lap2h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h100i, lap2h010i, lap2h110i, lap6h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hCubed(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d600(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d060(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d420(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real d240(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d500i, d050i, d510i, d150i, d330i, lap2hSq200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq100i, lap2hSq010i, lap2hSq110i
real lap4h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h100i, lap4h010i, lap4h110i, lap2h300i, lap2h030i, lap2h310i
real lap2h130i, d800i, d080i, d700i, d070i, d710i, d170i, d530i, d350i, lap8h, lap2hCubed200i
real lap2hCubed020i, lap2hCubed100i, lap2hCubed010i, lap2hCubed110i, lap2h4p, lap6h200i, lap6h020i, lap6h100i, lap6h010i, lap6h110i, lap4h400i
real lap4h040i, lap4h300i, lap4h030i, lap4h310i, lap4h130i, lap2h600i, lap2h060i, lap2h500i, lap2h050i, lap2h510i, lap2h150i
real lap2h330i, lap6hSq, lap4hSq200i, lap4hSq020i, lap4hSq100i, lap4hSq010i, lap4hSq110i, lap2hSq400i, lap2hSq040i, lap2hSq300i, lap2hSq030i
real lap2hSq310i, lap2hSq130i, lap4hCubed
#endMacro