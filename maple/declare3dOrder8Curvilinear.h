! ===========================================================================
!   Modified Equation : order=8, DIMENSIONS=3, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare3dOrder8Curvilinear()
real cr0, cr1, cr2, cr3, cs0, cs1, cs2, cs3, ct0, ct1, ct2
real ct3, crr0, crr1, crr2, crr3, css0, css1, css2, css3, ctt0, ctt1
real ctt2, ctt3, crrr0, crrr1, crrr2, crrr3, csss0, csss1, csss2, csss3, cttt0
real cttt1, cttt2, cttt3, crrrr0, crrrr1, crrrr2, crrrr3, cssss0, cssss1, cssss2, cssss3
real ctttt0, ctttt1, ctttt2, ctttt3, crrrrr0, crrrrr1, crrrrr2, crrrrr3, csssss0, csssss1, csssss2
real csssss3, cttttt0, cttttt1, cttttt2, cttttt3, crrrrrr0, crrrrrr1, crrrrrr2, crrrrrr3, cssssss0, cssssss1
real cssssss2, cssssss3, ctttttt0, ctttttt1, ctttttt2, ctttttt3, dr1, dr2, dr3, dr1i, dr2i
real dr3i, rx, ry, rz, sx, sy, sz, tx, ty, tz, diffOrder1
real diffOrder2, diffOrder3, rxr, rxs, rxt, ryr, rys, ryt, rzr, rzs, rzt
real sxr, sxs, sxt, syr, sys, syt, szr, szs, szt, txr, txs
real txt, tyr, tys, tyt, tzr, tzs, tzt, rxx, ryy, rzz, sxx
real syy, szz, txx, tyy, tzz, d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d100i, d010i
real d110i, d001i, d101i, d011i, d400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d004(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d202(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real d022(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d300i, d030i, d003i, d310i, d130i, d301i, d103i, d031i, d013i, lap2h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real lap2h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h100i, lap2h010i, lap2h110i, lap2h001i, lap2h101i, lap2h011i, lap6h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4hSq(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hCubed(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real d600(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d060(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d006(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d420(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d240(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d402(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d204(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d042(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d024(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d500i, d050i
real d005i, d510i, d150i, d330i, d501i, d105i, d051i, d015i, d303i, d033i, lap2hSq200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real lap2hSq020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2hSq100i, lap2hSq010i, lap2hSq001i, lap2hSq110i, lap2hSq101i, lap2hSq011i, lap4h200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0)
real lap2h400(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h040(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h004(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h220(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h202(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h022(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap4h100i, lap4h010i, lap4h001i, lap4h110i, lap4h101i
real lap4h011i, lap2h300i, lap2h030i, lap2h003i, lap2h310i, lap2h130i, lap2h301i, lap2h103i, lap2h031i, lap2h013i, d800i
real d080i, d008i, d700i, d070i, d007i, d710i, d170i, d701i, d107i, d071i, d017i
real d530i, d350i, d503i, d305i, d053i, d035i, lap8h, lap2hCubed200i, lap2hCubed020i, lap2hCubed002i, lap2hCubed100i
real lap2hCubed010i, lap2hCubed001i, lap2hCubed110i, lap2hCubed101i, lap2hCubed011i, lap2h4p, lap6h200i, lap6h020i, lap6h002i, lap6h100i, lap6h010i
real lap6h001i, lap6h110i, lap6h101i, lap6h011i, lap4h400i, lap4h040i, lap4h004i, lap4h300i, lap4h030i, lap4h003i, lap4h310i
real lap4h130i, lap4h301i, lap4h103i, lap4h031i, lap4h013i, lap2h600i, lap2h060i, lap2h006i, lap2h500i, lap2h050i, lap2h005i
real lap2h510i, lap2h150i, lap2h330i, lap2h501i, lap2h105i, lap2h051i, lap2h015i, lap2h303i, lap2h033i, lap6hSq, lap4hSq200i
real lap4hSq020i, lap4hSq002i, lap4hSq100i, lap4hSq010i, lap4hSq001i, lap4hSq110i, lap4hSq101i, lap4hSq011i, lap2hSq400i, lap2hSq040i, lap2hSq004i
real lap2hSq300i, lap2hSq030i, lap2hSq003i, lap2hSq310i, lap2hSq130i, lap2hSq301i, lap2hSq103i, lap2hSq031i, lap2hSq013i, lap4hCubed
#endMacro