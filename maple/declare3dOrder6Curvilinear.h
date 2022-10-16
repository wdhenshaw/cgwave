! ===========================================================================
!   Modified Equation : order=6, DIMENSIONS=3, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare3dOrder6Curvilinear()
real cr0, cr1, cr2, cs0, cs1, cs2, ct0, ct1, ct2, crr0, crr1
real crr2, css0, css1, css2, ctt0, ctt1, ctt2, crrr0, crrr1, crrr2, csss0
real csss1, csss2, cttt0, cttt1, cttt2, crrrr0, crrrr1, crrrr2, cssss0, cssss1, cssss2
real ctttt0, ctttt1, ctttt2, dr1, dr2, dr3, dr1i, dr2i, dr3i, rx, ry
real rz, sx, sy, sz, tx, ty, tz, diffOrder1, diffOrder2, diffOrder3, rxr
real rxs, rxt, ryr, rys, ryt, rzr, rzs, rzt, sxr, sxs, sxt
real syr, sys, syt, szr, szs, szt, txr, txs, txt, tyr, tys
real tyt, tzr, tzs, tzt, rxx, ryy, rzz, sxx, syy, szz, txx
real tyy, tzz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d100i, d010i, d110i, d001i, d101i
real d011i, d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d004(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d202(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d022(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d300i, d030i
real d003i, d310i, d130i, d301i, d103i, d031i, d013i, lap2h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h100i
real lap2h010i, lap2h110i, lap2h001i, lap2h101i, lap2h011i, lap6h, lap4hSq, lap2hCubed, d600i, d060i, d006i
real d500i, d050i, d005i, d510i, d150i, d330i, d501i, d105i, d051i, d015i, d303i
real d033i, lap2hSq200i, lap2hSq020i, lap2hSq002i, lap2hSq100i, lap2hSq010i, lap2hSq001i, lap2hSq110i, lap2hSq101i, lap2hSq011i, lap4h200i
real lap4h020i, lap4h002i, lap2h400i, lap2h040i, lap2h004i, lap4h100i, lap4h010i, lap4h001i, lap4h110i, lap4h101i, lap4h011i
real lap2h300i, lap2h030i, lap2h003i, lap2h310i, lap2h130i, lap2h301i, lap2h103i, lap2h031i, lap2h013i
#endMacro