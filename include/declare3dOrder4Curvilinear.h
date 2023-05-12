! ===========================================================================
!   Modified Equation : order=4, DIMENSIONS=3, gridType=Curvilinear
!  Macro to DECLARE VARIABLES
!  Macro created by maple code: cgwave/maple/writeModifiedEquationCode.mpl
! ===========================================================================
#beginMacro declare3dOrder4Curvilinear()
real cr0, cr1, cs0, cs1, ct0, ct1, crr0, crr1, css0, css1, ctt0
real ctt1, dr1, dr2, dr3, dr1i, dr2i, dr3i, rx, ry, rz, sx
real sy, sz, tx, ty, tz, diffOrder1, diffOrder2, diffOrder3, rxr, rxs, rxt
real ryr, rys, ryt, rzr, rzs, rzt, sxr, sxs, sxt, syr, sys
real syt, szr, szs, szt, txr, txs, txt, tyr, tys, tyt, tzr
real tzs, tzt, rxx, ryy, rzz, sxx, syy, szz, txx, tyy, tzz
real d200(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d020(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d002(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), lap2h(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:0), d100i, d010i, d110i, d001i, d101i, d011i, d400i
real d040i, d004i, lap4h, lap2hSq, d300i, d030i, d003i, d310i, d130i, d301i, d103i
real d031i, d013i, lap2h200i, lap2h020i, lap2h002i, lap2h100i, lap2h010i, lap2h110i, lap2h001i, lap2h101i, lap2h011i

#endMacro