! This file automatically generated from advWave.bf90 with bpp.
        subroutine advWave2dOrder2c( nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,mask,xy,rsxy,um,u,un,f,fa,v,vh,lapCoeff,bc,frequencyArray,ipar,rpar,ierr )
    ! subroutine advWave2dOrder2c(nd,n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b,!                 mask,xy,rsxy,  um,u,un, f,fa, v, vh,  bc, frequencyArray, ipar, rpar, ierr )
   !======================================================================
   !   Advance a time step for Waves equations
   !
   ! nd : number of space dimensions
   ! um,u,un : u(t-dt), u(t), u(t+dt)
   !
   ! ipar(0)  = option : option=0 - advance wave equation
   !                           =1 - add upwind dissipation (predictor corrector mode)
   !
   !======================================================================
        implicit none
        integer nd, n1a,n1b,n2a,n2b,n3a,n3b,nd1a,nd1b,nd2a,nd2b,nd3a,nd3b,nd4a,nd4b
        real um(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
        real u(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
        real un(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
        real f(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
        real fa(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b,0:*)  ! forcings at different times
        real xy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1) 
        real v(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)
        real vh(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,nd4a:nd4b)  ! holds current Helmholtz solutions
        real lapCoeff(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:*)  ! holds coeff of Laplacian for HA scheme
        real rsxy(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b,0:nd-1,0:nd-1)
        integer mask(nd1a:nd1b,nd2a:nd2b,nd3a:nd3b)
        integer bc(0:1,0:2),ierr
        real frequencyArray(0:*)
        integer ipar(0:*)
        real rpar(0:*)
   !     ---- local variables -----
        integer m1a,m1b,m2a,m2b,m3a,m3b,numGhost,nStart,nEnd,mt
        integer c,i1,i2,i3,n,gridType,orderOfAccuracy,orderInTime,axis,dir,grid,freq
        integer addForcing,orderOfDissipation,option,gridIsImplicit,preComputeUpwindUt,modifiedEquationApproach
        integer useNewForcingMethod,numberOfForcingFunctions,fcur,fnext,fprev,numberOfFrequencies
        real t,tm,cc,dt,dy,dz,cdt,cdtdx,cdtdy,cdtdz
    ! ,adc,adcdt,add,adddt
        real dt4by12
    ! logical addDissipation
        integer debug
        integer adjustHelmholtzForUpwinding
        real dx(0:2),dr(0:2)
        real dx2i,dy2i,dz2i,dxsqi,dysqi,dzsqi,dxi,dyi,dzi
        real dx12i,dy12i,dz12i,dxsq12i,dysq12i,dzsq12i,dxy4i,dxz4i,dyz4,time0,time1
        real dxi4,dyi4,dzi4,dxdyi2,dxdzi2,dydzi2
        real c0,c1,csq,dtsq,cdtsq,cdtsq12,cdtSqBy12
        real gridCFL
        integer maxOrderOfAccuracy
        parameter( maxOrderOfAccuracy=12 )
    ! Coefficients in the implicit scheme
        real bImp(0:maxOrderOfAccuracy-1)
        real cImp(-1:1,0:maxOrderOfAccuracy-1)
        real alpha2,alpha4,alpha6,alpha8, beta2,beta4,beta6,beta8
        integer rectangular,curvilinear
        parameter( rectangular=0, curvilinear=1 )
        integer timeSteppingMethod
        integer defaultTimeStepping,adamsSymmetricOrder3,rungeKuttaFourthOrder,stoermerTimeStepping,modifiedEquationTimeStepping
        parameter(defaultTimeStepping=0,adamsSymmetricOrder3=1,rungeKuttaFourthOrder=2,stoermerTimeStepping=3,modifiedEquationTimeStepping=4)
   !- ! Dispersion models
   !- integer noDispersion,drude
   !- parameter( noDispersion=0, drude=1 )
   !...........start statement function
        integer kd,m
        real rx,ry,rz,sx,sy,sz,tx,ty,tz
     ! --- declare variables used in the difference approximations (defined below) ---
            real d12
            real d22
            real h12
            real h22
            real rxr2
            real rxs2
            real rxt2
            real rxrr2
            real rxss2
            real rxrs2
            real ryr2
            real rys2
            real ryt2
            real ryrr2
            real ryss2
            real ryrs2
            real rzr2
            real rzs2
            real rzt2
            real rzrr2
            real rzss2
            real rzrs2
            real sxr2
            real sxs2
            real sxt2
            real sxrr2
            real sxss2
            real sxrs2
            real syr2
            real sys2
            real syt2
            real syrr2
            real syss2
            real syrs2
            real szr2
            real szs2
            real szt2
            real szrr2
            real szss2
            real szrs2
            real txr2
            real txs2
            real txt2
            real txrr2
            real txss2
            real txrs2
            real tyr2
            real tys2
            real tyt2
            real tyrr2
            real tyss2
            real tyrs2
            real tzr2
            real tzs2
            real tzt2
            real tzrr2
            real tzss2
            real tzrs2
            real rxx21
            real rxx22
            real rxy22
            real rxx23
            real rxy23
            real rxz23
            real ryx22
            real ryy22
            real ryx23
            real ryy23
            real ryz23
            real rzx22
            real rzy22
            real rzx23
            real rzy23
            real rzz23
            real sxx22
            real sxy22
            real sxx23
            real sxy23
            real sxz23
            real syx22
            real syy22
            real syx23
            real syy23
            real syz23
            real szx22
            real szy22
            real szx23
            real szy23
            real szz23
            real txx22
            real txy22
            real txx23
            real txy23
            real txz23
            real tyx22
            real tyy22
            real tyx23
            real tyy23
            real tyz23
            real tzx22
            real tzy22
            real tzx23
            real tzy23
            real tzz23
            real ur2
            real us2
            real ut2
            real urr2
            real uss2
            real urs2
            real utt2
            real urt2
            real ust2
            real urrr2
            real usss2
            real uttt2
            real ux21
            real uy21
            real uz21
            real ux22
            real uy22
            real uz22
            real ux23
            real uy23
            real uz23
            real uxx21
            real uyy21
            real uxy21
            real uxz21
            real uyz21
            real uzz21
            real ulaplacian21
            real uxx22
            real uyy22
            real uxy22
            real uxz22
            real uyz22
            real uzz22
            real ulaplacian22
            real uxx23
            real uyy23
            real uzz23
            real uxy23
            real uxz23
            real uyz23
            real ulaplacian23
            real ux23r
            real uy23r
            real uz23r
            real uxx23r
            real uyy23r
            real uxy23r
            real uzz23r
            real uxz23r
            real uyz23r
            real ux21r
            real uy21r
            real uz21r
            real uxx21r
            real uyy21r
            real uzz21r
            real uxy21r
            real uxz21r
            real uyz21r
            real ulaplacian21r
            real ux22r
            real uy22r
            real uz22r
            real uxx22r
            real uyy22r
            real uzz22r
            real uxy22r
            real uxz22r
            real uyz22r
            real ulaplacian22r
            real ulaplacian23r
            real uxxx22r
            real uyyy22r
            real uxxy22r
            real uxyy22r
            real uxxxx22r
            real uyyyy22r
            real uxxyy22r
            real uxxx23r
            real uyyy23r
            real uzzz23r
            real uxxy23r
            real uxxz23r
            real uxyy23r
            real uyyz23r
            real uxzz23r
            real uyzz23r
            real uxxxx23r
            real uyyyy23r
            real uzzzz23r
            real uxxyy23r
            real uxxzz23r
            real uyyzz23r
            real uLapSq22r
            real uLapSq23r
            real umr2
            real ums2
            real umt2
            real umrr2
            real umss2
            real umrs2
            real umtt2
            real umrt2
            real umst2
            real umrrr2
            real umsss2
            real umttt2
            real umx21
            real umy21
            real umz21
            real umx22
            real umy22
            real umz22
            real umx23
            real umy23
            real umz23
            real umxx21
            real umyy21
            real umxy21
            real umxz21
            real umyz21
            real umzz21
            real umlaplacian21
            real umxx22
            real umyy22
            real umxy22
            real umxz22
            real umyz22
            real umzz22
            real umlaplacian22
            real umxx23
            real umyy23
            real umzz23
            real umxy23
            real umxz23
            real umyz23
            real umlaplacian23
            real umx23r
            real umy23r
            real umz23r
            real umxx23r
            real umyy23r
            real umxy23r
            real umzz23r
            real umxz23r
            real umyz23r
            real umx21r
            real umy21r
            real umz21r
            real umxx21r
            real umyy21r
            real umzz21r
            real umxy21r
            real umxz21r
            real umyz21r
            real umlaplacian21r
            real umx22r
            real umy22r
            real umz22r
            real umxx22r
            real umyy22r
            real umzz22r
            real umxy22r
            real umxz22r
            real umyz22r
            real umlaplacian22r
            real umlaplacian23r
            real umxxx22r
            real umyyy22r
            real umxxy22r
            real umxyy22r
            real umxxxx22r
            real umyyyy22r
            real umxxyy22r
            real umxxx23r
            real umyyy23r
            real umzzz23r
            real umxxy23r
            real umxxz23r
            real umxyy23r
            real umyyz23r
            real umxzz23r
            real umyzz23r
            real umxxxx23r
            real umyyyy23r
            real umzzzz23r
            real umxxyy23r
            real umxxzz23r
            real umyyzz23r
            real umLapSq22r
            real umLapSq23r
            real vr2
            real vs2
            real vt2
            real vrr2
            real vss2
            real vrs2
            real vtt2
            real vrt2
            real vst2
            real vrrr2
            real vsss2
            real vttt2
            real vx21
            real vy21
            real vz21
            real vx22
            real vy22
            real vz22
            real vx23
            real vy23
            real vz23
            real vxx21
            real vyy21
            real vxy21
            real vxz21
            real vyz21
            real vzz21
            real vlaplacian21
            real vxx22
            real vyy22
            real vxy22
            real vxz22
            real vyz22
            real vzz22
            real vlaplacian22
            real vxx23
            real vyy23
            real vzz23
            real vxy23
            real vxz23
            real vyz23
            real vlaplacian23
            real vx23r
            real vy23r
            real vz23r
            real vxx23r
            real vyy23r
            real vxy23r
            real vzz23r
            real vxz23r
            real vyz23r
            real vx21r
            real vy21r
            real vz21r
            real vxx21r
            real vyy21r
            real vzz21r
            real vxy21r
            real vxz21r
            real vyz21r
            real vlaplacian21r
            real vx22r
            real vy22r
            real vz22r
            real vxx22r
            real vyy22r
            real vzz22r
            real vxy22r
            real vxz22r
            real vyz22r
            real vlaplacian22r
            real vlaplacian23r
            real vxxx22r
            real vyyy22r
            real vxxy22r
            real vxyy22r
            real vxxxx22r
            real vyyyy22r
            real vxxyy22r
            real vxxx23r
            real vyyy23r
            real vzzz23r
            real vxxy23r
            real vxxz23r
            real vxyy23r
            real vyyz23r
            real vxzz23r
            real vyzz23r
            real vxxxx23r
            real vyyyy23r
            real vzzzz23r
            real vxxyy23r
            real vxxzz23r
            real vyyzz23r
            real vLapSq22r
            real vLapSq23r
   !  declareDifferenceOrder2(um,none)
            real fr2
            real fs2
            real ft2
            real frr2
            real fss2
            real frs2
            real ftt2
            real frt2
            real fst2
            real frrr2
            real fsss2
            real fttt2
            real fx21
            real fy21
            real fz21
            real fx22
            real fy22
            real fz22
            real fx23
            real fy23
            real fz23
            real fxx21
            real fyy21
            real fxy21
            real fxz21
            real fyz21
            real fzz21
            real flaplacian21
            real fxx22
            real fyy22
            real fxy22
            real fxz22
            real fyz22
            real fzz22
            real flaplacian22
            real fxx23
            real fyy23
            real fzz23
            real fxy23
            real fxz23
            real fyz23
            real flaplacian23
            real fx23r
            real fy23r
            real fz23r
            real fxx23r
            real fyy23r
            real fxy23r
            real fzz23r
            real fxz23r
            real fyz23r
            real fx21r
            real fy21r
            real fz21r
            real fxx21r
            real fyy21r
            real fzz21r
            real fxy21r
            real fxz21r
            real fyz21r
            real flaplacian21r
            real fx22r
            real fy22r
            real fz22r
            real fxx22r
            real fyy22r
            real fzz22r
            real fxy22r
            real fxz22r
            real fyz22r
            real flaplacian22r
            real flaplacian23r
            real fxxx22r
            real fyyy22r
            real fxxy22r
            real fxyy22r
            real fxxxx22r
            real fyyyy22r
            real fxxyy22r
            real fxxx23r
            real fyyy23r
            real fzzz23r
            real fxxy23r
            real fxxz23r
            real fxyy23r
            real fyyz23r
            real fxzz23r
            real fyzz23r
            real fxxxx23r
            real fyyyy23r
            real fzzzz23r
            real fxxyy23r
            real fxxzz23r
            real fyyzz23r
            real fLapSq22r
            real fLapSq23r
            real d14
            real d24
            real h41
            real h42
            real rxr4
            real rxs4
            real rxt4
            real ryr4
            real rys4
            real ryt4
            real rzr4
            real rzs4
            real rzt4
            real sxr4
            real sxs4
            real sxt4
            real syr4
            real sys4
            real syt4
            real szr4
            real szs4
            real szt4
            real txr4
            real txs4
            real txt4
            real tyr4
            real tys4
            real tyt4
            real tzr4
            real tzs4
            real tzt4
            real rxx41
            real rxx42
            real rxy42
            real rxx43
            real rxy43
            real rxz43
            real ryx42
            real ryy42
            real ryx43
            real ryy43
            real ryz43
            real rzx42
            real rzy42
            real rzx43
            real rzy43
            real rzz43
            real sxx42
            real sxy42
            real sxx43
            real sxy43
            real sxz43
            real syx42
            real syy42
            real syx43
            real syy43
            real syz43
            real szx42
            real szy42
            real szx43
            real szy43
            real szz43
            real txx42
            real txy42
            real txx43
            real txy43
            real txz43
            real tyx42
            real tyy42
            real tyx43
            real tyy43
            real tyz43
            real tzx42
            real tzy42
            real tzx43
            real tzy43
            real tzz43
            real ur4
            real us4
            real ut4
            real urr4
            real uss4
            real utt4
            real urs4
            real urt4
            real ust4
            real ux41
            real uy41
            real uz41
            real ux42
            real uy42
            real uz42
            real ux43
            real uy43
            real uz43
            real uxx41
            real uyy41
            real uxy41
            real uxz41
            real uyz41
            real uzz41
            real ulaplacian41
            real uxx42
            real uyy42
            real uxy42
            real uxz42
            real uyz42
            real uzz42
            real ulaplacian42
            real uxx43
            real uyy43
            real uzz43
            real uxy43
            real uxz43
            real uyz43
            real ulaplacian43
            real ux43r
            real uy43r
            real uz43r
            real uxx43r
            real uyy43r
            real uzz43r
            real uxy43r
            real uxz43r
            real uyz43r
            real ux41r
            real uy41r
            real uz41r
            real uxx41r
            real uyy41r
            real uzz41r
            real uxy41r
            real uxz41r
            real uyz41r
            real ulaplacian41r
            real ux42r
            real uy42r
            real uz42r
            real uxx42r
            real uyy42r
            real uzz42r
            real uxy42r
            real uxz42r
            real uyz42r
            real ulaplacian42r
            real ulaplacian43r
            real umr4
            real ums4
            real umt4
            real umrr4
            real umss4
            real umtt4
            real umrs4
            real umrt4
            real umst4
            real umx41
            real umy41
            real umz41
            real umx42
            real umy42
            real umz42
            real umx43
            real umy43
            real umz43
            real umxx41
            real umyy41
            real umxy41
            real umxz41
            real umyz41
            real umzz41
            real umlaplacian41
            real umxx42
            real umyy42
            real umxy42
            real umxz42
            real umyz42
            real umzz42
            real umlaplacian42
            real umxx43
            real umyy43
            real umzz43
            real umxy43
            real umxz43
            real umyz43
            real umlaplacian43
            real umx43r
            real umy43r
            real umz43r
            real umxx43r
            real umyy43r
            real umzz43r
            real umxy43r
            real umxz43r
            real umyz43r
            real umx41r
            real umy41r
            real umz41r
            real umxx41r
            real umyy41r
            real umzz41r
            real umxy41r
            real umxz41r
            real umyz41r
            real umlaplacian41r
            real umx42r
            real umy42r
            real umz42r
            real umxx42r
            real umyy42r
            real umzz42r
            real umxy42r
            real umxz42r
            real umyz42r
            real umlaplacian42r
            real umlaplacian43r
            real vr4
            real vs4
            real vt4
            real vrr4
            real vss4
            real vtt4
            real vrs4
            real vrt4
            real vst4
            real vx41
            real vy41
            real vz41
            real vx42
            real vy42
            real vz42
            real vx43
            real vy43
            real vz43
            real vxx41
            real vyy41
            real vxy41
            real vxz41
            real vyz41
            real vzz41
            real vlaplacian41
            real vxx42
            real vyy42
            real vxy42
            real vxz42
            real vyz42
            real vzz42
            real vlaplacian42
            real vxx43
            real vyy43
            real vzz43
            real vxy43
            real vxz43
            real vyz43
            real vlaplacian43
            real vx43r
            real vy43r
            real vz43r
            real vxx43r
            real vyy43r
            real vzz43r
            real vxy43r
            real vxz43r
            real vyz43r
            real vx41r
            real vy41r
            real vz41r
            real vxx41r
            real vyy41r
            real vzz41r
            real vxy41r
            real vxz41r
            real vyz41r
            real vlaplacian41r
            real vx42r
            real vy42r
            real vz42r
            real vxx42r
            real vyy42r
            real vzz42r
            real vxy42r
            real vxz42r
            real vyz42r
            real vlaplacian42r
            real vlaplacian43r
            real d16
            real d26
            real h16
            real h26
            real rxr6
            real rxs6
            real rxt6
            real ryr6
            real rys6
            real ryt6
            real rzr6
            real rzs6
            real rzt6
            real sxr6
            real sxs6
            real sxt6
            real syr6
            real sys6
            real syt6
            real szr6
            real szs6
            real szt6
            real txr6
            real txs6
            real txt6
            real tyr6
            real tys6
            real tyt6
            real tzr6
            real tzs6
            real tzt6
            real rxx61
            real rxx62
            real rxy62
            real rxx63
            real rxy63
            real rxz63
            real ryx62
            real ryy62
            real ryx63
            real ryy63
            real ryz63
            real rzx62
            real rzy62
            real rzx63
            real rzy63
            real rzz63
            real sxx62
            real sxy62
            real sxx63
            real sxy63
            real sxz63
            real syx62
            real syy62
            real syx63
            real syy63
            real syz63
            real szx62
            real szy62
            real szx63
            real szy63
            real szz63
            real txx62
            real txy62
            real txx63
            real txy63
            real txz63
            real tyx62
            real tyy62
            real tyx63
            real tyy63
            real tyz63
            real tzx62
            real tzy62
            real tzx63
            real tzy63
            real tzz63
            real ur6
            real us6
            real ut6
            real urr6
            real uss6
            real utt6
            real urs6
            real urt6
            real ust6
            real ux61
            real uy61
            real uz61
            real ux62
            real uy62
            real uz62
            real ux63
            real uy63
            real uz63
            real uxx61
            real uyy61
            real uxy61
            real uxz61
            real uyz61
            real uzz61
            real ulaplacian61
            real uxx62
            real uyy62
            real uxy62
            real uxz62
            real uyz62
            real uzz62
            real ulaplacian62
            real uxx63
            real uyy63
            real uzz63
            real uxy63
            real uxz63
            real uyz63
            real ulaplacian63
            real ux63r
            real uy63r
            real uz63r
            real uxx63r
            real uyy63r
            real uzz63r
            real uxy63r
            real uxz63r
            real uyz63r
            real ux61r
            real uy61r
            real uz61r
            real uxx61r
            real uyy61r
            real uzz61r
            real uxy61r
            real uxz61r
            real uyz61r
            real ulaplacian61r
            real ux62r
            real uy62r
            real uz62r
            real uxx62r
            real uyy62r
            real uzz62r
            real uxy62r
            real uxz62r
            real uyz62r
            real ulaplacian62r
            real ulaplacian63r
            real umr6
            real ums6
            real umt6
            real umrr6
            real umss6
            real umtt6
            real umrs6
            real umrt6
            real umst6
            real umx61
            real umy61
            real umz61
            real umx62
            real umy62
            real umz62
            real umx63
            real umy63
            real umz63
            real umxx61
            real umyy61
            real umxy61
            real umxz61
            real umyz61
            real umzz61
            real umlaplacian61
            real umxx62
            real umyy62
            real umxy62
            real umxz62
            real umyz62
            real umzz62
            real umlaplacian62
            real umxx63
            real umyy63
            real umzz63
            real umxy63
            real umxz63
            real umyz63
            real umlaplacian63
            real umx63r
            real umy63r
            real umz63r
            real umxx63r
            real umyy63r
            real umzz63r
            real umxy63r
            real umxz63r
            real umyz63r
            real umx61r
            real umy61r
            real umz61r
            real umxx61r
            real umyy61r
            real umzz61r
            real umxy61r
            real umxz61r
            real umyz61r
            real umlaplacian61r
            real umx62r
            real umy62r
            real umz62r
            real umxx62r
            real umyy62r
            real umzz62r
            real umxy62r
            real umxz62r
            real umyz62r
            real umlaplacian62r
            real umlaplacian63r
            real vr6
            real vs6
            real vt6
            real vrr6
            real vss6
            real vtt6
            real vrs6
            real vrt6
            real vst6
            real vx61
            real vy61
            real vz61
            real vx62
            real vy62
            real vz62
            real vx63
            real vy63
            real vz63
            real vxx61
            real vyy61
            real vxy61
            real vxz61
            real vyz61
            real vzz61
            real vlaplacian61
            real vxx62
            real vyy62
            real vxy62
            real vxz62
            real vyz62
            real vzz62
            real vlaplacian62
            real vxx63
            real vyy63
            real vzz63
            real vxy63
            real vxz63
            real vyz63
            real vlaplacian63
            real vx63r
            real vy63r
            real vz63r
            real vxx63r
            real vyy63r
            real vzz63r
            real vxy63r
            real vxz63r
            real vyz63r
            real vx61r
            real vy61r
            real vz61r
            real vxx61r
            real vyy61r
            real vzz61r
            real vxy61r
            real vxz61r
            real vyz61r
            real vlaplacian61r
            real vx62r
            real vy62r
            real vz62r
            real vxx62r
            real vyy62r
            real vzz62r
            real vxy62r
            real vxz62r
            real vyz62r
            real vlaplacian62r
            real vlaplacian63r
            real d18
            real d28
            real h18
            real h28
            real rxr8
            real rxs8
            real rxt8
            real ryr8
            real rys8
            real ryt8
            real rzr8
            real rzs8
            real rzt8
            real sxr8
            real sxs8
            real sxt8
            real syr8
            real sys8
            real syt8
            real szr8
            real szs8
            real szt8
            real txr8
            real txs8
            real txt8
            real tyr8
            real tys8
            real tyt8
            real tzr8
            real tzs8
            real tzt8
            real rxx81
            real rxx82
            real rxy82
            real rxx83
            real rxy83
            real rxz83
            real ryx82
            real ryy82
            real ryx83
            real ryy83
            real ryz83
            real rzx82
            real rzy82
            real rzx83
            real rzy83
            real rzz83
            real sxx82
            real sxy82
            real sxx83
            real sxy83
            real sxz83
            real syx82
            real syy82
            real syx83
            real syy83
            real syz83
            real szx82
            real szy82
            real szx83
            real szy83
            real szz83
            real txx82
            real txy82
            real txx83
            real txy83
            real txz83
            real tyx82
            real tyy82
            real tyx83
            real tyy83
            real tyz83
            real tzx82
            real tzy82
            real tzx83
            real tzy83
            real tzz83
            real ur8
            real us8
            real ut8
            real urr8
            real uss8
            real utt8
            real urs8
            real urt8
            real ust8
            real ux81
            real uy81
            real uz81
            real ux82
            real uy82
            real uz82
            real ux83
            real uy83
            real uz83
            real uxx81
            real uyy81
            real uxy81
            real uxz81
            real uyz81
            real uzz81
            real ulaplacian81
            real uxx82
            real uyy82
            real uxy82
            real uxz82
            real uyz82
            real uzz82
            real ulaplacian82
            real uxx83
            real uyy83
            real uzz83
            real uxy83
            real uxz83
            real uyz83
            real ulaplacian83
            real ux83r
            real uy83r
            real uz83r
            real uxx83r
            real uyy83r
            real uzz83r
            real uxy83r
            real uxz83r
            real uyz83r
            real ux81r
            real uy81r
            real uz81r
            real uxx81r
            real uyy81r
            real uzz81r
            real uxy81r
            real uxz81r
            real uyz81r
            real ulaplacian81r
            real ux82r
            real uy82r
            real uz82r
            real uxx82r
            real uyy82r
            real uzz82r
            real uxy82r
            real uxz82r
            real uyz82r
            real ulaplacian82r
            real ulaplacian83r
            real umr8
            real ums8
            real umt8
            real umrr8
            real umss8
            real umtt8
            real umrs8
            real umrt8
            real umst8
            real umx81
            real umy81
            real umz81
            real umx82
            real umy82
            real umz82
            real umx83
            real umy83
            real umz83
            real umxx81
            real umyy81
            real umxy81
            real umxz81
            real umyz81
            real umzz81
            real umlaplacian81
            real umxx82
            real umyy82
            real umxy82
            real umxz82
            real umyz82
            real umzz82
            real umlaplacian82
            real umxx83
            real umyy83
            real umzz83
            real umxy83
            real umxz83
            real umyz83
            real umlaplacian83
            real umx83r
            real umy83r
            real umz83r
            real umxx83r
            real umyy83r
            real umzz83r
            real umxy83r
            real umxz83r
            real umyz83r
            real umx81r
            real umy81r
            real umz81r
            real umxx81r
            real umyy81r
            real umzz81r
            real umxy81r
            real umxz81r
            real umyz81r
            real umlaplacian81r
            real umx82r
            real umy82r
            real umz82r
            real umxx82r
            real umyy82r
            real umzz82r
            real umxy82r
            real umxz82r
            real umyz82r
            real umlaplacian82r
            real umlaplacian83r
            real vr8
            real vs8
            real vt8
            real vrr8
            real vss8
            real vtt8
            real vrs8
            real vrt8
            real vst8
            real vx81
            real vy81
            real vz81
            real vx82
            real vy82
            real vz82
            real vx83
            real vy83
            real vz83
            real vxx81
            real vyy81
            real vxy81
            real vxz81
            real vyz81
            real vzz81
            real vlaplacian81
            real vxx82
            real vyy82
            real vxy82
            real vxz82
            real vyz82
            real vzz82
            real vlaplacian82
            real vxx83
            real vyy83
            real vzz83
            real vxy83
            real vxz83
            real vyz83
            real vlaplacian83
            real vx83r
            real vy83r
            real vz83r
            real vxx83r
            real vyy83r
            real vzz83r
            real vxy83r
            real vxz83r
            real vyz83r
            real vx81r
            real vy81r
            real vz81r
            real vxx81r
            real vyy81r
            real vzz81r
            real vxy81r
            real vxz81r
            real vyz81r
            real vlaplacian81r
            real vx82r
            real vy82r
            real vz82r
            real vxx82r
            real vyy82r
            real vzz82r
            real vxy82r
            real vxz82r
            real vyz82r
            real vlaplacian82r
            real vlaplacian83r
    ! define variables for getDerivatives macros
! ****** File written by makeGetDerivativesMacros.maple  ******
real ur
real urr
real urrr
real urrrr
real urrrrr
real urrrrrr
real us
real urs
real urrs
real urrrs
real urrrrs
real urrrrrs
real uss
real urss
real urrss
real urrrss
real urrrrss
real usss
real ursss
real urrsss
real urrrsss
real ussss
real urssss
real urrssss
real usssss
real ursssss
real ussssss
real ut
real urt
real urrt
real urrrt
real urrrrt
real urrrrrt
real ust
real urst
real urrst
real urrrst
real urrrrst
real usst
real ursst
real urrsst
real urrrsst
real ussst
real urssst
real urrssst
real usssst
real ursssst
real ussssst
real utt
real urtt
real urrtt
real urrrtt
real urrrrtt
real ustt
real urstt
real urrstt
real urrrstt
real usstt
real ursstt
real urrsstt
real ussstt
real urssstt
real usssstt
real uttt
real urttt
real urrttt
real urrrttt
real usttt
real ursttt
real urrsttt
real ussttt
real urssttt
real usssttt
real utttt
real urtttt
real urrtttt
real ustttt
real urstttt
real usstttt
real uttttt
real urttttt
real usttttt
real utttttt
real rxr
real rxrr
real rxrrr
real rxrrrr
real rxrrrrr
real rxs
real rxrs
real rxrrs
real rxrrrs
real rxrrrrs
real rxss
real rxrss
real rxrrss
real rxrrrss
real rxsss
real rxrsss
real rxrrsss
real rxssss
real rxrssss
real rxsssss
real rxt
real rxrt
real rxrrt
real rxrrrt
real rxrrrrt
real rxst
real rxrst
real rxrrst
real rxrrrst
real rxsst
real rxrsst
real rxrrsst
real rxssst
real rxrssst
real rxsssst
real rxtt
real rxrtt
real rxrrtt
real rxrrrtt
real rxstt
real rxrstt
real rxrrstt
real rxsstt
real rxrsstt
real rxssstt
real rxttt
real rxrttt
real rxrrttt
real rxsttt
real rxrsttt
real rxssttt
real rxtttt
real rxrtttt
real rxstttt
real rxttttt
real ryr
real ryrr
real ryrrr
real ryrrrr
real ryrrrrr
real rys
real ryrs
real ryrrs
real ryrrrs
real ryrrrrs
real ryss
real ryrss
real ryrrss
real ryrrrss
real rysss
real ryrsss
real ryrrsss
real ryssss
real ryrssss
real rysssss
real ryt
real ryrt
real ryrrt
real ryrrrt
real ryrrrrt
real ryst
real ryrst
real ryrrst
real ryrrrst
real rysst
real ryrsst
real ryrrsst
real ryssst
real ryrssst
real rysssst
real rytt
real ryrtt
real ryrrtt
real ryrrrtt
real rystt
real ryrstt
real ryrrstt
real rysstt
real ryrsstt
real ryssstt
real ryttt
real ryrttt
real ryrrttt
real rysttt
real ryrsttt
real ryssttt
real rytttt
real ryrtttt
real rystttt
real ryttttt
real sxr
real sxrr
real sxrrr
real sxrrrr
real sxrrrrr
real sxs
real sxrs
real sxrrs
real sxrrrs
real sxrrrrs
real sxss
real sxrss
real sxrrss
real sxrrrss
real sxsss
real sxrsss
real sxrrsss
real sxssss
real sxrssss
real sxsssss
real sxt
real sxrt
real sxrrt
real sxrrrt
real sxrrrrt
real sxst
real sxrst
real sxrrst
real sxrrrst
real sxsst
real sxrsst
real sxrrsst
real sxssst
real sxrssst
real sxsssst
real sxtt
real sxrtt
real sxrrtt
real sxrrrtt
real sxstt
real sxrstt
real sxrrstt
real sxsstt
real sxrsstt
real sxssstt
real sxttt
real sxrttt
real sxrrttt
real sxsttt
real sxrsttt
real sxssttt
real sxtttt
real sxrtttt
real sxstttt
real sxttttt
real syr
real syrr
real syrrr
real syrrrr
real syrrrrr
real sys
real syrs
real syrrs
real syrrrs
real syrrrrs
real syss
real syrss
real syrrss
real syrrrss
real sysss
real syrsss
real syrrsss
real syssss
real syrssss
real sysssss
real syt
real syrt
real syrrt
real syrrrt
real syrrrrt
real syst
real syrst
real syrrst
real syrrrst
real sysst
real syrsst
real syrrsst
real syssst
real syrssst
real sysssst
real sytt
real syrtt
real syrrtt
real syrrrtt
real systt
real syrstt
real syrrstt
real sysstt
real syrsstt
real syssstt
real syttt
real syrttt
real syrrttt
real systtt
real syrsttt
real syssttt
real sytttt
real syrtttt
real systttt
real syttttt
real rzr
real rzrr
real rzrrr
real rzrrrr
real rzrrrrr
real rzs
real rzrs
real rzrrs
real rzrrrs
real rzrrrrs
real rzss
real rzrss
real rzrrss
real rzrrrss
real rzsss
real rzrsss
real rzrrsss
real rzssss
real rzrssss
real rzsssss
real rzt
real rzrt
real rzrrt
real rzrrrt
real rzrrrrt
real rzst
real rzrst
real rzrrst
real rzrrrst
real rzsst
real rzrsst
real rzrrsst
real rzssst
real rzrssst
real rzsssst
real rztt
real rzrtt
real rzrrtt
real rzrrrtt
real rzstt
real rzrstt
real rzrrstt
real rzsstt
real rzrsstt
real rzssstt
real rzttt
real rzrttt
real rzrrttt
real rzsttt
real rzrsttt
real rzssttt
real rztttt
real rzrtttt
real rzstttt
real rzttttt
real szr
real szrr
real szrrr
real szrrrr
real szrrrrr
real szs
real szrs
real szrrs
real szrrrs
real szrrrrs
real szss
real szrss
real szrrss
real szrrrss
real szsss
real szrsss
real szrrsss
real szssss
real szrssss
real szsssss
real szt
real szrt
real szrrt
real szrrrt
real szrrrrt
real szst
real szrst
real szrrst
real szrrrst
real szsst
real szrsst
real szrrsst
real szssst
real szrssst
real szsssst
real sztt
real szrtt
real szrrtt
real szrrrtt
real szstt
real szrstt
real szrrstt
real szsstt
real szrsstt
real szssstt
real szttt
real szrttt
real szrrttt
real szsttt
real szrsttt
real szssttt
real sztttt
real szrtttt
real szstttt
real szttttt
real txr
real txrr
real txrrr
real txrrrr
real txrrrrr
real txs
real txrs
real txrrs
real txrrrs
real txrrrrs
real txss
real txrss
real txrrss
real txrrrss
real txsss
real txrsss
real txrrsss
real txssss
real txrssss
real txsssss
real txt
real txrt
real txrrt
real txrrrt
real txrrrrt
real txst
real txrst
real txrrst
real txrrrst
real txsst
real txrsst
real txrrsst
real txssst
real txrssst
real txsssst
real txtt
real txrtt
real txrrtt
real txrrrtt
real txstt
real txrstt
real txrrstt
real txsstt
real txrsstt
real txssstt
real txttt
real txrttt
real txrrttt
real txsttt
real txrsttt
real txssttt
real txtttt
real txrtttt
real txstttt
real txttttt
real tyr
real tyrr
real tyrrr
real tyrrrr
real tyrrrrr
real tys
real tyrs
real tyrrs
real tyrrrs
real tyrrrrs
real tyss
real tyrss
real tyrrss
real tyrrrss
real tysss
real tyrsss
real tyrrsss
real tyssss
real tyrssss
real tysssss
real tyt
real tyrt
real tyrrt
real tyrrrt
real tyrrrrt
real tyst
real tyrst
real tyrrst
real tyrrrst
real tysst
real tyrsst
real tyrrsst
real tyssst
real tyrssst
real tysssst
real tytt
real tyrtt
real tyrrtt
real tyrrrtt
real tystt
real tyrstt
real tyrrstt
real tysstt
real tyrsstt
real tyssstt
real tyttt
real tyrttt
real tyrrttt
real tysttt
real tyrsttt
real tyssttt
real tytttt
real tyrtttt
real tystttt
real tyttttt
real tzr
real tzrr
real tzrrr
real tzrrrr
real tzrrrrr
real tzs
real tzrs
real tzrrs
real tzrrrs
real tzrrrrs
real tzss
real tzrss
real tzrrss
real tzrrrss
real tzsss
real tzrsss
real tzrrsss
real tzssss
real tzrssss
real tzsssss
real tzt
real tzrt
real tzrrt
real tzrrrt
real tzrrrrt
real tzst
real tzrst
real tzrrst
real tzrrrst
real tzsst
real tzrsst
real tzrrsst
real tzssst
real tzrssst
real tzsssst
real tztt
real tzrtt
real tzrrtt
real tzrrrtt
real tzstt
real tzrstt
real tzrrstt
real tzsstt
real tzrsstt
real tzssstt
real tzttt
real tzrttt
real tzrrttt
real tzsttt
real tzrsttt
real tzssttt
real tztttt
real tzrtttt
real tzstttt
real tzttttt
real rxi
real ryi
real sxi
real syi
real rzi
real szi
real txi
real tyi
real tzi
real rxx
real rxy
real rxz
real ryy
real ryz
real rzz
real sxx
real sxy
real sxz
real syy
real syz
real szz
real txx
real txy
real txz
real tyy
real tyz
real tzz
real rxxx
real rxxy
real rxyy
real rxxz
real rxyz
real rxzz
real ryyy
real ryyz
real ryzz
real rzzz
real sxxx
real sxxy
real sxyy
real sxxz
real sxyz
real sxzz
real syyy
real syyz
real syzz
real szzz
real txxx
real txxy
real txyy
real txxz
real txyz
real txzz
real tyyy
real tyyz
real tyzz
real tzzz
real rxxxx
real rxxxy
real rxxyy
real rxyyy
real rxxxz
real rxxyz
real rxyyz
real rxxzz
real rxyzz
real rxzzz
real ryyyy
real ryyyz
real ryyzz
real ryzzz
real rzzzz
real sxxxx
real sxxxy
real sxxyy
real sxyyy
real sxxxz
real sxxyz
real sxyyz
real sxxzz
real sxyzz
real sxzzz
real syyyy
real syyyz
real syyzz
real syzzz
real szzzz
real txxxx
real txxxy
real txxyy
real txyyy
real txxxz
real txxyz
real txyyz
real txxzz
real txyzz
real txzzz
real tyyyy
real tyyyz
real tyyzz
real tyzzz
real tzzzz
real rxxxxx
real rxxxxy
real rxxxyy
real rxxyyy
real rxyyyy
real rxxxxz
real rxxxyz
real rxxyyz
real rxyyyz
real rxxxzz
real rxxyzz
real rxyyzz
real rxxzzz
real rxyzzz
real rxzzzz
real ryyyyy
real ryyyyz
real ryyyzz
real ryyzzz
real ryzzzz
real rzzzzz
real sxxxxx
real sxxxxy
real sxxxyy
real sxxyyy
real sxyyyy
real sxxxxz
real sxxxyz
real sxxyyz
real sxyyyz
real sxxxzz
real sxxyzz
real sxyyzz
real sxxzzz
real sxyzzz
real sxzzzz
real syyyyy
real syyyyz
real syyyzz
real syyzzz
real syzzzz
real szzzzz
real txxxxx
real txxxxy
real txxxyy
real txxyyy
real txyyyy
real txxxxz
real txxxyz
real txxyyz
real txyyyz
real txxxzz
real txxyzz
real txyyzz
real txxzzz
real txyzzz
real txzzzz
real tyyyyy
real tyyyyz
real tyyyzz
real tyyzzz
real tyzzzz
real tzzzzz
real rxxxxxx
real rxxxxxy
real rxxxxyy
real rxxxyyy
real rxxyyyy
real rxyyyyy
real rxxxxxz
real rxxxxyz
real rxxxyyz
real rxxyyyz
real rxyyyyz
real rxxxxzz
real rxxxyzz
real rxxyyzz
real rxyyyzz
real rxxxzzz
real rxxyzzz
real rxyyzzz
real rxxzzzz
real rxyzzzz
real rxzzzzz
real ryyyyyy
real ryyyyyz
real ryyyyzz
real ryyyzzz
real ryyzzzz
real ryzzzzz
real rzzzzzz
real sxxxxxx
real sxxxxxy
real sxxxxyy
real sxxxyyy
real sxxyyyy
real sxyyyyy
real sxxxxxz
real sxxxxyz
real sxxxyyz
real sxxyyyz
real sxyyyyz
real sxxxxzz
real sxxxyzz
real sxxyyzz
real sxyyyzz
real sxxxzzz
real sxxyzzz
real sxyyzzz
real sxxzzzz
real sxyzzzz
real sxzzzzz
real syyyyyy
real syyyyyz
real syyyyzz
real syyyzzz
real syyzzzz
real syzzzzz
real szzzzzz
real txxxxxx
real txxxxxy
real txxxxyy
real txxxyyy
real txxyyyy
real txyyyyy
real txxxxxz
real txxxxyz
real txxxyyz
real txxyyyz
real txyyyyz
real txxxxzz
real txxxyzz
real txxyyzz
real txyyyzz
real txxxzzz
real txxyzzz
real txyyzzz
real txxzzzz
real txyzzzz
real txzzzzz
real tyyyyyy
real tyyyyyz
real tyyyyzz
real tyyyzzz
real tyyzzzz
real tyzzzzz
real tzzzzzz
real uxx
real uxxxx
real uxxxxxx
real uyy
real uxxyy
real uxxxxyy
real uyyyy
real uxxyyyy
real uyyyyyy
real uzz
real uxxzz
real uxxxxzz
real uyyzz
real uxxyyzz
real uyyyyzz
real uzzzz
real uxxzzzz
real uyyzzzz
real uzzzzzz
    ! real cdt4by360,cdt6by20160
        real cdtPow4By12,cdtPow6By360
        real lap2d2,lap3d2,lap2d4,lap3d4,lap2d6,lap3d6,lap2d8,lap3d8,lap2d2Pow2,lap3d2Pow2,lap2d2Pow3,lap3d2Pow3,lap2d2Pow4,lap3d2Pow4,lap2d4Pow2,lap3d4Pow2,lap2d4Pow3,lap3d4Pow3,lap2d6Pow2,lap3d6Pow2
        real lap2d2m,lap3d2m
        real du,fd22d,fd23d,fd42d,fd43d,fd62d,fd63d,fd82d,fd83d
        real DztU
    ! forcing correction functions: 
        real lap2d2f,f2drme44, lap3d2f, f3drme44, f2dcme44, f3dcme44, ff
    ! real cdSosupx,cdSosupy,cdSosupz
        real adSosup,sosupParameter, uDotFactor, adxSosup(0:2)
        integer useSosupDissipation,sosupDissipationOption
        integer updateSolution,updateDissipation,computeUt
        integer ec 
        real ep 
        real fv(0:1) , ev(0:1), evtt(0:1), evxx(0:1), evyy(0:1), evzz(0:1)
        real evxxxx(0:1), evxxyy(0:1), evyyyy(0:1), evxxzz(0:1), evyyzz(0:1), evzzzz(0:1), evtttt(0:1)
        real evtttttt(0:1)
        real evxxxxxx(0:1)
        real evyyyyyy(0:1)
        real evzzzzzz(0:1)       
        real evxxyyyy(0:1)
        real evxxxxyy(0:1)
        real evxxxxzz(0:1)
        real evxxzzzz(0:1)
        real evyyyyzz(0:1)
        real evyyzzzz(0:1)
        real evxxyyzz(0:1)
        real omega, coswt
        integer maxFreq
        parameter( maxFreq=500 )
        real cosFreqt(0:maxFreq), coswtAve(0:maxFreq), cosineFactor(0:maxFreq)
        integer idv(0:2),j1,j2,j3
        integer iStencil,upwCase,upwindHalfStencilWidth,i1l,i2l,i3l, i1r,i2r,i3r
        integer useUpwindDissipation,useImplicitUpwindDissipation,adjustOmega,solveHelmholtz
        real upw,maxDiff,umj
        real upwindCoeff(-3:3,0:3) 
        integer forcingOption
    ! forcingOptions -- these should match ForcingEnum in CgWave.h 
    ! enum ForcingOptionEnum
    ! {
    !   noForcing=0,
    !   twilightZoneForcing,
    !   userForcing,
    !   helmholtzForcing
    ! };
        integer noForcing,twilightZoneForcing,userForcing,helmholtzForcing
        parameter(noForcing           =0,twilightZoneForcing =1,userForcing         =2,helmholtzForcing    =3 )
   ! real unxx22r,unyy22r,unxy22r,unx22r
   !.......statement functions for jacobian
        rx(i1,i2,i3)=rsxy(i1,i2,i3,0,0)
        ry(i1,i2,i3)=rsxy(i1,i2,i3,0,1)
        rz(i1,i2,i3)=rsxy(i1,i2,i3,0,2)
        sx(i1,i2,i3)=rsxy(i1,i2,i3,1,0)
        sy(i1,i2,i3)=rsxy(i1,i2,i3,1,1)
        sz(i1,i2,i3)=rsxy(i1,i2,i3,1,2)
        tx(i1,i2,i3)=rsxy(i1,i2,i3,2,0)
        ty(i1,i2,i3)=rsxy(i1,i2,i3,2,1)
        tz(i1,i2,i3)=rsxy(i1,i2,i3,2,2)
    !  The next macro will define the difference approximation statement functions for u
        d12(kd) = 1./(2.*dr(kd))
        d22(kd) = 1./(dr(kd)**2)
        ur2(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*d12(0)
        us2(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*d12(1)
        ut2(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*d12(2)
        urr2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) )*d22(0)
        uss2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) )*d22(1)
        urs2(i1,i2,i3,kd)=(ur2(i1,i2+1,i3,kd)-ur2(i1,i2-1,i3,kd))*d12(1)
        utt2(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd)) )*d22(2)
        urt2(i1,i2,i3,kd)=(ur2(i1,i2,i3+1,kd)-ur2(i1,i2,i3-1,kd))*d12(2)
        ust2(i1,i2,i3,kd)=(us2(i1,i2,i3+1,kd)-us2(i1,i2,i3-1,kd))*d12(2)
        urrr2(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
        usss2(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
        uttt2(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
        rxr2(i1,i2,i3)=(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))*d12(0)
        rxs2(i1,i2,i3)=(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))*d12(1)
        rxt2(i1,i2,i3)=(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))*d12(2)
        rxrr2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1+1,i2,i3)+rx(i1-1,i2,i3)) )*d22(0)
        rxss2(i1,i2,i3)=(-2.*rx(i1,i2,i3)+(rx(i1,i2+1,i3)+rx(i1,i2-1,i3)) )*d22(1)
        rxrs2(i1,i2,i3)=(rxr2(i1,i2+1,i3)-rxr2(i1,i2-1,i3))*d12(1)
        ryr2(i1,i2,i3)=(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))*d12(0)
        rys2(i1,i2,i3)=(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))*d12(1)
        ryt2(i1,i2,i3)=(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))*d12(2)
        ryrr2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1+1,i2,i3)+ry(i1-1,i2,i3)) )*d22(0)
        ryss2(i1,i2,i3)=(-2.*ry(i1,i2,i3)+(ry(i1,i2+1,i3)+ry(i1,i2-1,i3)) )*d22(1)
        ryrs2(i1,i2,i3)=(ryr2(i1,i2+1,i3)-ryr2(i1,i2-1,i3))*d12(1)
        rzr2(i1,i2,i3)=(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))*d12(0)
        rzs2(i1,i2,i3)=(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))*d12(1)
        rzt2(i1,i2,i3)=(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))*d12(2)
        rzrr2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1+1,i2,i3)+rz(i1-1,i2,i3)) )*d22(0)
        rzss2(i1,i2,i3)=(-2.*rz(i1,i2,i3)+(rz(i1,i2+1,i3)+rz(i1,i2-1,i3)) )*d22(1)
        rzrs2(i1,i2,i3)=(rzr2(i1,i2+1,i3)-rzr2(i1,i2-1,i3))*d12(1)
        sxr2(i1,i2,i3)=(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))*d12(0)
        sxs2(i1,i2,i3)=(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))*d12(1)
        sxt2(i1,i2,i3)=(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))*d12(2)
        sxrr2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1+1,i2,i3)+sx(i1-1,i2,i3)) )*d22(0)
        sxss2(i1,i2,i3)=(-2.*sx(i1,i2,i3)+(sx(i1,i2+1,i3)+sx(i1,i2-1,i3)) )*d22(1)
        sxrs2(i1,i2,i3)=(sxr2(i1,i2+1,i3)-sxr2(i1,i2-1,i3))*d12(1)
        syr2(i1,i2,i3)=(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))*d12(0)
        sys2(i1,i2,i3)=(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))*d12(1)
        syt2(i1,i2,i3)=(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))*d12(2)
        syrr2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1+1,i2,i3)+sy(i1-1,i2,i3)) )*d22(0)
        syss2(i1,i2,i3)=(-2.*sy(i1,i2,i3)+(sy(i1,i2+1,i3)+sy(i1,i2-1,i3)) )*d22(1)
        syrs2(i1,i2,i3)=(syr2(i1,i2+1,i3)-syr2(i1,i2-1,i3))*d12(1)
        szr2(i1,i2,i3)=(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))*d12(0)
        szs2(i1,i2,i3)=(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))*d12(1)
        szt2(i1,i2,i3)=(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))*d12(2)
        szrr2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1+1,i2,i3)+sz(i1-1,i2,i3)) )*d22(0)
        szss2(i1,i2,i3)=(-2.*sz(i1,i2,i3)+(sz(i1,i2+1,i3)+sz(i1,i2-1,i3)) )*d22(1)
        szrs2(i1,i2,i3)=(szr2(i1,i2+1,i3)-szr2(i1,i2-1,i3))*d12(1)
        txr2(i1,i2,i3)=(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))*d12(0)
        txs2(i1,i2,i3)=(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))*d12(1)
        txt2(i1,i2,i3)=(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))*d12(2)
        txrr2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1+1,i2,i3)+tx(i1-1,i2,i3)) )*d22(0)
        txss2(i1,i2,i3)=(-2.*tx(i1,i2,i3)+(tx(i1,i2+1,i3)+tx(i1,i2-1,i3)) )*d22(1)
        txrs2(i1,i2,i3)=(txr2(i1,i2+1,i3)-txr2(i1,i2-1,i3))*d12(1)
        tyr2(i1,i2,i3)=(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))*d12(0)
        tys2(i1,i2,i3)=(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))*d12(1)
        tyt2(i1,i2,i3)=(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))*d12(2)
        tyrr2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1+1,i2,i3)+ty(i1-1,i2,i3)) )*d22(0)
        tyss2(i1,i2,i3)=(-2.*ty(i1,i2,i3)+(ty(i1,i2+1,i3)+ty(i1,i2-1,i3)) )*d22(1)
        tyrs2(i1,i2,i3)=(tyr2(i1,i2+1,i3)-tyr2(i1,i2-1,i3))*d12(1)
        tzr2(i1,i2,i3)=(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))*d12(0)
        tzs2(i1,i2,i3)=(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))*d12(1)
        tzt2(i1,i2,i3)=(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))*d12(2)
        tzrr2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1+1,i2,i3)+tz(i1-1,i2,i3)) )*d22(0)
        tzss2(i1,i2,i3)=(-2.*tz(i1,i2,i3)+(tz(i1,i2+1,i3)+tz(i1,i2-1,i3)) )*d22(1)
        tzrs2(i1,i2,i3)=(tzr2(i1,i2+1,i3)-tzr2(i1,i2-1,i3))*d12(1)
        ux21(i1,i2,i3,kd)= rx(i1,i2,i3)*ur2(i1,i2,i3,kd)
        uy21(i1,i2,i3,kd)=0
        uz21(i1,i2,i3,kd)=0
        ux22(i1,i2,i3,kd)= rx(i1,i2,i3)*ur2(i1,i2,i3,kd)+sx(i1,i2,i3)*us2(i1,i2,i3,kd)
        uy22(i1,i2,i3,kd)= ry(i1,i2,i3)*ur2(i1,i2,i3,kd)+sy(i1,i2,i3)*us2(i1,i2,i3,kd)
        uz22(i1,i2,i3,kd)=0
        ux23(i1,i2,i3,kd)=rx(i1,i2,i3)*ur2(i1,i2,i3,kd)+sx(i1,i2,i3)*us2(i1,i2,i3,kd)+tx(i1,i2,i3)*ut2(i1,i2,i3,kd)
        uy23(i1,i2,i3,kd)=ry(i1,i2,i3)*ur2(i1,i2,i3,kd)+sy(i1,i2,i3)*us2(i1,i2,i3,kd)+ty(i1,i2,i3)*ut2(i1,i2,i3,kd)
        uz23(i1,i2,i3,kd)=rz(i1,i2,i3)*ur2(i1,i2,i3,kd)+sz(i1,i2,i3)*us2(i1,i2,i3,kd)+tz(i1,i2,i3)*ut2(i1,i2,i3,kd)
        rxx21(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)
        rxx22(i1,i2,i3)= rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(i1,i2,i3)
        rxy22(i1,i2,i3)= ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(i1,i2,i3)
        rxx23(i1,i2,i3)=rx(i1,i2,i3)*rxr2(i1,i2,i3)+sx(i1,i2,i3)*rxs2(i1,i2,i3)+tx(i1,i2,i3)*rxt2(i1,i2,i3)
        rxy23(i1,i2,i3)=ry(i1,i2,i3)*rxr2(i1,i2,i3)+sy(i1,i2,i3)*rxs2(i1,i2,i3)+ty(i1,i2,i3)*rxt2(i1,i2,i3)
        rxz23(i1,i2,i3)=rz(i1,i2,i3)*rxr2(i1,i2,i3)+sz(i1,i2,i3)*rxs2(i1,i2,i3)+tz(i1,i2,i3)*rxt2(i1,i2,i3)
        ryx22(i1,i2,i3)= rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(i1,i2,i3)
        ryy22(i1,i2,i3)= ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(i1,i2,i3)
        ryx23(i1,i2,i3)=rx(i1,i2,i3)*ryr2(i1,i2,i3)+sx(i1,i2,i3)*rys2(i1,i2,i3)+tx(i1,i2,i3)*ryt2(i1,i2,i3)
        ryy23(i1,i2,i3)=ry(i1,i2,i3)*ryr2(i1,i2,i3)+sy(i1,i2,i3)*rys2(i1,i2,i3)+ty(i1,i2,i3)*ryt2(i1,i2,i3)
        ryz23(i1,i2,i3)=rz(i1,i2,i3)*ryr2(i1,i2,i3)+sz(i1,i2,i3)*rys2(i1,i2,i3)+tz(i1,i2,i3)*ryt2(i1,i2,i3)
        rzx22(i1,i2,i3)= rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(i1,i2,i3)
        rzy22(i1,i2,i3)= ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(i1,i2,i3)
        rzx23(i1,i2,i3)=rx(i1,i2,i3)*rzr2(i1,i2,i3)+sx(i1,i2,i3)*rzs2(i1,i2,i3)+tx(i1,i2,i3)*rzt2(i1,i2,i3)
        rzy23(i1,i2,i3)=ry(i1,i2,i3)*rzr2(i1,i2,i3)+sy(i1,i2,i3)*rzs2(i1,i2,i3)+ty(i1,i2,i3)*rzt2(i1,i2,i3)
        rzz23(i1,i2,i3)=rz(i1,i2,i3)*rzr2(i1,i2,i3)+sz(i1,i2,i3)*rzs2(i1,i2,i3)+tz(i1,i2,i3)*rzt2(i1,i2,i3)
        sxx22(i1,i2,i3)= rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(i1,i2,i3)
        sxy22(i1,i2,i3)= ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(i1,i2,i3)
        sxx23(i1,i2,i3)=rx(i1,i2,i3)*sxr2(i1,i2,i3)+sx(i1,i2,i3)*sxs2(i1,i2,i3)+tx(i1,i2,i3)*sxt2(i1,i2,i3)
        sxy23(i1,i2,i3)=ry(i1,i2,i3)*sxr2(i1,i2,i3)+sy(i1,i2,i3)*sxs2(i1,i2,i3)+ty(i1,i2,i3)*sxt2(i1,i2,i3)
        sxz23(i1,i2,i3)=rz(i1,i2,i3)*sxr2(i1,i2,i3)+sz(i1,i2,i3)*sxs2(i1,i2,i3)+tz(i1,i2,i3)*sxt2(i1,i2,i3)
        syx22(i1,i2,i3)= rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(i1,i2,i3)
        syy22(i1,i2,i3)= ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(i1,i2,i3)
        syx23(i1,i2,i3)=rx(i1,i2,i3)*syr2(i1,i2,i3)+sx(i1,i2,i3)*sys2(i1,i2,i3)+tx(i1,i2,i3)*syt2(i1,i2,i3)
        syy23(i1,i2,i3)=ry(i1,i2,i3)*syr2(i1,i2,i3)+sy(i1,i2,i3)*sys2(i1,i2,i3)+ty(i1,i2,i3)*syt2(i1,i2,i3)
        syz23(i1,i2,i3)=rz(i1,i2,i3)*syr2(i1,i2,i3)+sz(i1,i2,i3)*sys2(i1,i2,i3)+tz(i1,i2,i3)*syt2(i1,i2,i3)
        szx22(i1,i2,i3)= rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(i1,i2,i3)
        szy22(i1,i2,i3)= ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(i1,i2,i3)
        szx23(i1,i2,i3)=rx(i1,i2,i3)*szr2(i1,i2,i3)+sx(i1,i2,i3)*szs2(i1,i2,i3)+tx(i1,i2,i3)*szt2(i1,i2,i3)
        szy23(i1,i2,i3)=ry(i1,i2,i3)*szr2(i1,i2,i3)+sy(i1,i2,i3)*szs2(i1,i2,i3)+ty(i1,i2,i3)*szt2(i1,i2,i3)
        szz23(i1,i2,i3)=rz(i1,i2,i3)*szr2(i1,i2,i3)+sz(i1,i2,i3)*szs2(i1,i2,i3)+tz(i1,i2,i3)*szt2(i1,i2,i3)
        txx22(i1,i2,i3)= rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(i1,i2,i3)
        txy22(i1,i2,i3)= ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(i1,i2,i3)
        txx23(i1,i2,i3)=rx(i1,i2,i3)*txr2(i1,i2,i3)+sx(i1,i2,i3)*txs2(i1,i2,i3)+tx(i1,i2,i3)*txt2(i1,i2,i3)
        txy23(i1,i2,i3)=ry(i1,i2,i3)*txr2(i1,i2,i3)+sy(i1,i2,i3)*txs2(i1,i2,i3)+ty(i1,i2,i3)*txt2(i1,i2,i3)
        txz23(i1,i2,i3)=rz(i1,i2,i3)*txr2(i1,i2,i3)+sz(i1,i2,i3)*txs2(i1,i2,i3)+tz(i1,i2,i3)*txt2(i1,i2,i3)
        tyx22(i1,i2,i3)= rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(i1,i2,i3)
        tyy22(i1,i2,i3)= ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(i1,i2,i3)
        tyx23(i1,i2,i3)=rx(i1,i2,i3)*tyr2(i1,i2,i3)+sx(i1,i2,i3)*tys2(i1,i2,i3)+tx(i1,i2,i3)*tyt2(i1,i2,i3)
        tyy23(i1,i2,i3)=ry(i1,i2,i3)*tyr2(i1,i2,i3)+sy(i1,i2,i3)*tys2(i1,i2,i3)+ty(i1,i2,i3)*tyt2(i1,i2,i3)
        tyz23(i1,i2,i3)=rz(i1,i2,i3)*tyr2(i1,i2,i3)+sz(i1,i2,i3)*tys2(i1,i2,i3)+tz(i1,i2,i3)*tyt2(i1,i2,i3)
        tzx22(i1,i2,i3)= rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*tzs2(i1,i2,i3)
        tzy22(i1,i2,i3)= ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*tzs2(i1,i2,i3)
        tzx23(i1,i2,i3)=rx(i1,i2,i3)*tzr2(i1,i2,i3)+sx(i1,i2,i3)*tzs2(i1,i2,i3)+tx(i1,i2,i3)*tzt2(i1,i2,i3)
        tzy23(i1,i2,i3)=ry(i1,i2,i3)*tzr2(i1,i2,i3)+sy(i1,i2,i3)*tzs2(i1,i2,i3)+ty(i1,i2,i3)*tzt2(i1,i2,i3)
        tzz23(i1,i2,i3)=rz(i1,i2,i3)*tzr2(i1,i2,i3)+sz(i1,i2,i3)*tzs2(i1,i2,i3)+tz(i1,i2,i3)*tzt2(i1,i2,i3)
        uxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*ur2(i1,i2,i3,kd)
        uyy21(i1,i2,i3,kd)=0
        uxy21(i1,i2,i3,kd)=0
        uxz21(i1,i2,i3,kd)=0
        uyz21(i1,i2,i3,kd)=0
        uzz21(i1,i2,i3,kd)=0
        ulaplacian21(i1,i2,i3,kd)=uxx21(i1,i2,i3,kd)
        uxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx22(i1,i2,i3))*us2(i1,i2,i3,kd)
        uyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(syy22(i1,i2,i3))*us2(i1,i2,i3,kd)
        uxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss2(i1,i2,i3,kd)+rxy22(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*us2(i1,i2,i3,kd)
        uxz22(i1,i2,i3,kd)=0
        uyz22(i1,i2,i3,kd)=0
        uzz22(i1,i2,i3,kd)=0
        ulaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*us2(i1,i2,i3,kd)
        uxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sx(i1,i2,i3)**2*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*us2(i1,i2,i3,kd)+txx23(i1,i2,i3)*ut2(i1,i2,i3,kd)
        uyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sy(i1,i2,i3)**2*uss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*ur2(i1,i2,i3,kd)+syy23(i1,i2,i3)*us2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*ut2(i1,i2,i3,kd)
        uzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr2(i1,i2,i3,kd)+sz(i1,i2,i3)**2*uss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*urs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*urt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+szz23(i1,i2,i3)*us2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
        uxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*utt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust2(i1,i2,i3,kd)+rxy23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*us2(i1,i2,i3,kd)+txy23(i1,i2,i3)*ut2(i1,i2,i3,kd)
        uxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*uss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*utt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*urs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust2(i1,i2,i3,kd)+rxz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*us2(i1,i2,i3,kd)+txz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
        uyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr2(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*uss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*utt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*urs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*urt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust2(i1,i2,i3,kd)+ryz23(i1,i2,i3)*ur2(i1,i2,i3,kd)+syz23(i1,i2,i3)*us2(i1,i2,i3,kd)+tyz23(i1,i2,i3)*ut2(i1,i2,i3,kd)
        ulaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*urr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*uss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*utt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*urs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*urt2(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*ust2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*ur2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*us2(i1,i2,i3,kd)+(txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*ut2(i1,i2,i3,kd)
    !============================================================================================
    ! Define derivatives for a rectangular grid
    !
    !============================================================================================
        h12(kd) = 1./(2.*dx(kd))
        h22(kd) = 1./(dx(kd)**2)
        ux23r(i1,i2,i3,kd)=(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))*h12(0)
        uy23r(i1,i2,i3,kd)=(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))*h12(1)
        uz23r(i1,i2,i3,kd)=(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))*h12(2)
        uxx23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)) )*h22(0)
        uyy23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd)) )*h22(1)
        uxy23r(i1,i2,i3,kd)=(ux23r(i1,i2+1,i3,kd)-ux23r(i1,i2-1,i3,kd))*h12(1)
        uzz23r(i1,i2,i3,kd)=(-2.*u(i1,i2,i3,kd)+(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd)) )*h22(2)
        uxz23r(i1,i2,i3,kd)=(ux23r(i1,i2,i3+1,kd)-ux23r(i1,i2,i3-1,kd))*h12(2)
        uyz23r(i1,i2,i3,kd)=(uy23r(i1,i2,i3+1,kd)-uy23r(i1,i2,i3-1,kd))*h12(2)
        ux21r(i1,i2,i3,kd)= ux23r(i1,i2,i3,kd)
        uy21r(i1,i2,i3,kd)= uy23r(i1,i2,i3,kd)
        uz21r(i1,i2,i3,kd)= uz23r(i1,i2,i3,kd)
        uxx21r(i1,i2,i3,kd)= uxx23r(i1,i2,i3,kd)
        uyy21r(i1,i2,i3,kd)= uyy23r(i1,i2,i3,kd)
        uzz21r(i1,i2,i3,kd)= uzz23r(i1,i2,i3,kd)
        uxy21r(i1,i2,i3,kd)= uxy23r(i1,i2,i3,kd)
        uxz21r(i1,i2,i3,kd)= uxz23r(i1,i2,i3,kd)
        uyz21r(i1,i2,i3,kd)= uyz23r(i1,i2,i3,kd)
        ulaplacian21r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)
        ux22r(i1,i2,i3,kd)= ux23r(i1,i2,i3,kd)
        uy22r(i1,i2,i3,kd)= uy23r(i1,i2,i3,kd)
        uz22r(i1,i2,i3,kd)= uz23r(i1,i2,i3,kd)
        uxx22r(i1,i2,i3,kd)= uxx23r(i1,i2,i3,kd)
        uyy22r(i1,i2,i3,kd)= uyy23r(i1,i2,i3,kd)
        uzz22r(i1,i2,i3,kd)= uzz23r(i1,i2,i3,kd)
        uxy22r(i1,i2,i3,kd)= uxy23r(i1,i2,i3,kd)
        uxz22r(i1,i2,i3,kd)= uxz23r(i1,i2,i3,kd)
        uyz22r(i1,i2,i3,kd)= uyz23r(i1,i2,i3,kd)
        ulaplacian22r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)+uyy23r(i1,i2,i3,kd)
        ulaplacian23r(i1,i2,i3,kd)=uxx23r(i1,i2,i3,kd)+uyy23r(i1,i2,i3,kd)+uzz23r(i1,i2,i3,kd)
        uxxx22r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
        uyyy22r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
        uxxy22r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
        uxyy22r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
        uxxxx22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4)
        uyyyy22r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)
        uxxyy22r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
    ! 2D laplacian squared = u.xxxx + 2 u.xxyy + u.yyyy
        uLapSq22r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +2.*(u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
        uxxx23r(i1,i2,i3,kd)=(-2.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
        uyyy23r(i1,i2,i3,kd)=(-2.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
        uzzz23r(i1,i2,i3,kd)=(-2.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
        uxxy23r(i1,i2,i3,kd)=( uxx22r(i1,i2+1,i3,kd)-uxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
        uxyy23r(i1,i2,i3,kd)=( uyy22r(i1+1,i2,i3,kd)-uyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
        uxxz23r(i1,i2,i3,kd)=( uxx22r(i1,i2,i3+1,kd)-uxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
        uyyz23r(i1,i2,i3,kd)=( uyy22r(i1,i2,i3+1,kd)-uyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
        uxzz23r(i1,i2,i3,kd)=( uzz22r(i1+1,i2,i3,kd)-uzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
        uyzz23r(i1,i2,i3,kd)=( uzz22r(i1,i2+1,i3,kd)-uzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
        uxxxx23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))+(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4)
        uyyyy23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))+(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)
        uzzzz23r(i1,i2,i3,kd)=(6.*u(i1,i2,i3,kd)-4.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))+(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)**4)
        uxxyy23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))   +   (u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
        uxxzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd)+u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))   +   (u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
        uyyzz23r(i1,i2,i3,kd)=( 4.*u(i1,i2,i3,kd)     -2.*(u(i1,i2+1,i3,kd)  +u(i1,i2-1,i3,kd)+  u(i1,i2  ,i3+1,kd)+u(i1,i2  ,i3-1,kd))   +   (u(i1,i2+1,i3+1,kd)+u(i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
    ! 3D laplacian squared = u.xxxx + u.yyyy + u.zzzz + 2 (u.xxyy + u.xxzz + u.yyzz )
        uLapSq23r(i1,i2,i3,kd)= ( 6.*u(i1,i2,i3,kd)   - 4.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))    +(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))    +(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 6.*u(i1,i2,i3,kd)    -4.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))    +(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)  +u(i1-1,i2,i3,kd)  +u(i1  ,i2+1,i3,kd)+u(i1  ,i2-1,i3,kd))   +2.*(u(i1+1,i2+1,i3,kd)+u(i1-1,i2+1,i3,kd)+u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*u(i1,i2,i3,kd)     -4.*(u(i1+1,i2,i3,kd)  +u(i1-1,i2,i3,kd)  +u(i1  ,i2,i3+1,kd)+u(i1  ,i2,i3-1,kd))   +2.*(u(i1+1,i2,i3+1,kd)+u(i1-1,i2,i3+1,kd)+u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*u(i1,i2,i3,kd)     -4.*(u(i1,i2+1,i3,kd)  +u(i1,i2-1,i3,kd)  +u(i1,i2  ,i3+1,kd)+u(i1,i2  ,i3-1,kd))   +2.*(u(i1,i2+1,i3+1,kd)+u(i1,i2-1,i3+1,kd)+u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
        d14(kd) = 1./(12.*dr(kd))
        d24(kd) = 1./(12.*dr(kd)**2)
        ur4(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*d14(0)
        us4(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*d14(1)
        ut4(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*d14(2)
        urr4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*d24(0)
        uss4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*d24(1)
        utt4(i1,i2,i3,kd)=(-30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*d24(2)
        urs4(i1,i2,i3,kd)=(8.*(ur4(i1,i2+1,i3,kd)-ur4(i1,i2-1,i3,kd))-(ur4(i1,i2+2,i3,kd)-ur4(i1,i2-2,i3,kd)))*d14(1)
        urt4(i1,i2,i3,kd)=(8.*(ur4(i1,i2,i3+1,kd)-ur4(i1,i2,i3-1,kd))-(ur4(i1,i2,i3+2,kd)-ur4(i1,i2,i3-2,kd)))*d14(2)
        ust4(i1,i2,i3,kd)=(8.*(us4(i1,i2,i3+1,kd)-us4(i1,i2,i3-1,kd))-(us4(i1,i2,i3+2,kd)-us4(i1,i2,i3-2,kd)))*d14(2)
        rxr4(i1,i2,i3)=(8.*(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))-(rx(i1+2,i2,i3)-rx(i1-2,i2,i3)))*d14(0)
        rxs4(i1,i2,i3)=(8.*(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))-(rx(i1,i2+2,i3)-rx(i1,i2-2,i3)))*d14(1)
        rxt4(i1,i2,i3)=(8.*(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))-(rx(i1,i2,i3+2)-rx(i1,i2,i3-2)))*d14(2)
        ryr4(i1,i2,i3)=(8.*(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))-(ry(i1+2,i2,i3)-ry(i1-2,i2,i3)))*d14(0)
        rys4(i1,i2,i3)=(8.*(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))-(ry(i1,i2+2,i3)-ry(i1,i2-2,i3)))*d14(1)
        ryt4(i1,i2,i3)=(8.*(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))-(ry(i1,i2,i3+2)-ry(i1,i2,i3-2)))*d14(2)
        rzr4(i1,i2,i3)=(8.*(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))-(rz(i1+2,i2,i3)-rz(i1-2,i2,i3)))*d14(0)
        rzs4(i1,i2,i3)=(8.*(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))-(rz(i1,i2+2,i3)-rz(i1,i2-2,i3)))*d14(1)
        rzt4(i1,i2,i3)=(8.*(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))-(rz(i1,i2,i3+2)-rz(i1,i2,i3-2)))*d14(2)
        sxr4(i1,i2,i3)=(8.*(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))-(sx(i1+2,i2,i3)-sx(i1-2,i2,i3)))*d14(0)
        sxs4(i1,i2,i3)=(8.*(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))-(sx(i1,i2+2,i3)-sx(i1,i2-2,i3)))*d14(1)
        sxt4(i1,i2,i3)=(8.*(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))-(sx(i1,i2,i3+2)-sx(i1,i2,i3-2)))*d14(2)
        syr4(i1,i2,i3)=(8.*(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))-(sy(i1+2,i2,i3)-sy(i1-2,i2,i3)))*d14(0)
        sys4(i1,i2,i3)=(8.*(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))-(sy(i1,i2+2,i3)-sy(i1,i2-2,i3)))*d14(1)
        syt4(i1,i2,i3)=(8.*(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))-(sy(i1,i2,i3+2)-sy(i1,i2,i3-2)))*d14(2)
        szr4(i1,i2,i3)=(8.*(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))-(sz(i1+2,i2,i3)-sz(i1-2,i2,i3)))*d14(0)
        szs4(i1,i2,i3)=(8.*(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))-(sz(i1,i2+2,i3)-sz(i1,i2-2,i3)))*d14(1)
        szt4(i1,i2,i3)=(8.*(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))-(sz(i1,i2,i3+2)-sz(i1,i2,i3-2)))*d14(2)
        txr4(i1,i2,i3)=(8.*(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))-(tx(i1+2,i2,i3)-tx(i1-2,i2,i3)))*d14(0)
        txs4(i1,i2,i3)=(8.*(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))-(tx(i1,i2+2,i3)-tx(i1,i2-2,i3)))*d14(1)
        txt4(i1,i2,i3)=(8.*(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))-(tx(i1,i2,i3+2)-tx(i1,i2,i3-2)))*d14(2)
        tyr4(i1,i2,i3)=(8.*(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))-(ty(i1+2,i2,i3)-ty(i1-2,i2,i3)))*d14(0)
        tys4(i1,i2,i3)=(8.*(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))-(ty(i1,i2+2,i3)-ty(i1,i2-2,i3)))*d14(1)
        tyt4(i1,i2,i3)=(8.*(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))-(ty(i1,i2,i3+2)-ty(i1,i2,i3-2)))*d14(2)
        tzr4(i1,i2,i3)=(8.*(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))-(tz(i1+2,i2,i3)-tz(i1-2,i2,i3)))*d14(0)
        tzs4(i1,i2,i3)=(8.*(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))-(tz(i1,i2+2,i3)-tz(i1,i2-2,i3)))*d14(1)
        tzt4(i1,i2,i3)=(8.*(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))-(tz(i1,i2,i3+2)-tz(i1,i2,i3-2)))*d14(2)
        ux41(i1,i2,i3,kd)= rx(i1,i2,i3)*ur4(i1,i2,i3,kd)
        uy41(i1,i2,i3,kd)=0
        uz41(i1,i2,i3,kd)=0
        ux42(i1,i2,i3,kd)= rx(i1,i2,i3)*ur4(i1,i2,i3,kd)+sx(i1,i2,i3)*us4(i1,i2,i3,kd)
        uy42(i1,i2,i3,kd)= ry(i1,i2,i3)*ur4(i1,i2,i3,kd)+sy(i1,i2,i3)*us4(i1,i2,i3,kd)
        uz42(i1,i2,i3,kd)=0
        ux43(i1,i2,i3,kd)=rx(i1,i2,i3)*ur4(i1,i2,i3,kd)+sx(i1,i2,i3)*us4(i1,i2,i3,kd)+tx(i1,i2,i3)*ut4(i1,i2,i3,kd)
        uy43(i1,i2,i3,kd)=ry(i1,i2,i3)*ur4(i1,i2,i3,kd)+sy(i1,i2,i3)*us4(i1,i2,i3,kd)+ty(i1,i2,i3)*ut4(i1,i2,i3,kd)
        uz43(i1,i2,i3,kd)=rz(i1,i2,i3)*ur4(i1,i2,i3,kd)+sz(i1,i2,i3)*us4(i1,i2,i3,kd)+tz(i1,i2,i3)*ut4(i1,i2,i3,kd)
        rxx41(i1,i2,i3)= rx(i1,i2,i3)*rxr4(i1,i2,i3)
        rxx42(i1,i2,i3)= rx(i1,i2,i3)*rxr4(i1,i2,i3)+sx(i1,i2,i3)*rxs4(i1,i2,i3)
        rxy42(i1,i2,i3)= ry(i1,i2,i3)*rxr4(i1,i2,i3)+sy(i1,i2,i3)*rxs4(i1,i2,i3)
        rxx43(i1,i2,i3)=rx(i1,i2,i3)*rxr4(i1,i2,i3)+sx(i1,i2,i3)*rxs4(i1,i2,i3)+tx(i1,i2,i3)*rxt4(i1,i2,i3)
        rxy43(i1,i2,i3)=ry(i1,i2,i3)*rxr4(i1,i2,i3)+sy(i1,i2,i3)*rxs4(i1,i2,i3)+ty(i1,i2,i3)*rxt4(i1,i2,i3)
        rxz43(i1,i2,i3)=rz(i1,i2,i3)*rxr4(i1,i2,i3)+sz(i1,i2,i3)*rxs4(i1,i2,i3)+tz(i1,i2,i3)*rxt4(i1,i2,i3)
        ryx42(i1,i2,i3)= rx(i1,i2,i3)*ryr4(i1,i2,i3)+sx(i1,i2,i3)*rys4(i1,i2,i3)
        ryy42(i1,i2,i3)= ry(i1,i2,i3)*ryr4(i1,i2,i3)+sy(i1,i2,i3)*rys4(i1,i2,i3)
        ryx43(i1,i2,i3)=rx(i1,i2,i3)*ryr4(i1,i2,i3)+sx(i1,i2,i3)*rys4(i1,i2,i3)+tx(i1,i2,i3)*ryt4(i1,i2,i3)
        ryy43(i1,i2,i3)=ry(i1,i2,i3)*ryr4(i1,i2,i3)+sy(i1,i2,i3)*rys4(i1,i2,i3)+ty(i1,i2,i3)*ryt4(i1,i2,i3)
        ryz43(i1,i2,i3)=rz(i1,i2,i3)*ryr4(i1,i2,i3)+sz(i1,i2,i3)*rys4(i1,i2,i3)+tz(i1,i2,i3)*ryt4(i1,i2,i3)
        rzx42(i1,i2,i3)= rx(i1,i2,i3)*rzr4(i1,i2,i3)+sx(i1,i2,i3)*rzs4(i1,i2,i3)
        rzy42(i1,i2,i3)= ry(i1,i2,i3)*rzr4(i1,i2,i3)+sy(i1,i2,i3)*rzs4(i1,i2,i3)
        rzx43(i1,i2,i3)=rx(i1,i2,i3)*rzr4(i1,i2,i3)+sx(i1,i2,i3)*rzs4(i1,i2,i3)+tx(i1,i2,i3)*rzt4(i1,i2,i3)
        rzy43(i1,i2,i3)=ry(i1,i2,i3)*rzr4(i1,i2,i3)+sy(i1,i2,i3)*rzs4(i1,i2,i3)+ty(i1,i2,i3)*rzt4(i1,i2,i3)
        rzz43(i1,i2,i3)=rz(i1,i2,i3)*rzr4(i1,i2,i3)+sz(i1,i2,i3)*rzs4(i1,i2,i3)+tz(i1,i2,i3)*rzt4(i1,i2,i3)
        sxx42(i1,i2,i3)= rx(i1,i2,i3)*sxr4(i1,i2,i3)+sx(i1,i2,i3)*sxs4(i1,i2,i3)
        sxy42(i1,i2,i3)= ry(i1,i2,i3)*sxr4(i1,i2,i3)+sy(i1,i2,i3)*sxs4(i1,i2,i3)
        sxx43(i1,i2,i3)=rx(i1,i2,i3)*sxr4(i1,i2,i3)+sx(i1,i2,i3)*sxs4(i1,i2,i3)+tx(i1,i2,i3)*sxt4(i1,i2,i3)
        sxy43(i1,i2,i3)=ry(i1,i2,i3)*sxr4(i1,i2,i3)+sy(i1,i2,i3)*sxs4(i1,i2,i3)+ty(i1,i2,i3)*sxt4(i1,i2,i3)
        sxz43(i1,i2,i3)=rz(i1,i2,i3)*sxr4(i1,i2,i3)+sz(i1,i2,i3)*sxs4(i1,i2,i3)+tz(i1,i2,i3)*sxt4(i1,i2,i3)
        syx42(i1,i2,i3)= rx(i1,i2,i3)*syr4(i1,i2,i3)+sx(i1,i2,i3)*sys4(i1,i2,i3)
        syy42(i1,i2,i3)= ry(i1,i2,i3)*syr4(i1,i2,i3)+sy(i1,i2,i3)*sys4(i1,i2,i3)
        syx43(i1,i2,i3)=rx(i1,i2,i3)*syr4(i1,i2,i3)+sx(i1,i2,i3)*sys4(i1,i2,i3)+tx(i1,i2,i3)*syt4(i1,i2,i3)
        syy43(i1,i2,i3)=ry(i1,i2,i3)*syr4(i1,i2,i3)+sy(i1,i2,i3)*sys4(i1,i2,i3)+ty(i1,i2,i3)*syt4(i1,i2,i3)
        syz43(i1,i2,i3)=rz(i1,i2,i3)*syr4(i1,i2,i3)+sz(i1,i2,i3)*sys4(i1,i2,i3)+tz(i1,i2,i3)*syt4(i1,i2,i3)
        szx42(i1,i2,i3)= rx(i1,i2,i3)*szr4(i1,i2,i3)+sx(i1,i2,i3)*szs4(i1,i2,i3)
        szy42(i1,i2,i3)= ry(i1,i2,i3)*szr4(i1,i2,i3)+sy(i1,i2,i3)*szs4(i1,i2,i3)
        szx43(i1,i2,i3)=rx(i1,i2,i3)*szr4(i1,i2,i3)+sx(i1,i2,i3)*szs4(i1,i2,i3)+tx(i1,i2,i3)*szt4(i1,i2,i3)
        szy43(i1,i2,i3)=ry(i1,i2,i3)*szr4(i1,i2,i3)+sy(i1,i2,i3)*szs4(i1,i2,i3)+ty(i1,i2,i3)*szt4(i1,i2,i3)
        szz43(i1,i2,i3)=rz(i1,i2,i3)*szr4(i1,i2,i3)+sz(i1,i2,i3)*szs4(i1,i2,i3)+tz(i1,i2,i3)*szt4(i1,i2,i3)
        txx42(i1,i2,i3)= rx(i1,i2,i3)*txr4(i1,i2,i3)+sx(i1,i2,i3)*txs4(i1,i2,i3)
        txy42(i1,i2,i3)= ry(i1,i2,i3)*txr4(i1,i2,i3)+sy(i1,i2,i3)*txs4(i1,i2,i3)
        txx43(i1,i2,i3)=rx(i1,i2,i3)*txr4(i1,i2,i3)+sx(i1,i2,i3)*txs4(i1,i2,i3)+tx(i1,i2,i3)*txt4(i1,i2,i3)
        txy43(i1,i2,i3)=ry(i1,i2,i3)*txr4(i1,i2,i3)+sy(i1,i2,i3)*txs4(i1,i2,i3)+ty(i1,i2,i3)*txt4(i1,i2,i3)
        txz43(i1,i2,i3)=rz(i1,i2,i3)*txr4(i1,i2,i3)+sz(i1,i2,i3)*txs4(i1,i2,i3)+tz(i1,i2,i3)*txt4(i1,i2,i3)
        tyx42(i1,i2,i3)= rx(i1,i2,i3)*tyr4(i1,i2,i3)+sx(i1,i2,i3)*tys4(i1,i2,i3)
        tyy42(i1,i2,i3)= ry(i1,i2,i3)*tyr4(i1,i2,i3)+sy(i1,i2,i3)*tys4(i1,i2,i3)
        tyx43(i1,i2,i3)=rx(i1,i2,i3)*tyr4(i1,i2,i3)+sx(i1,i2,i3)*tys4(i1,i2,i3)+tx(i1,i2,i3)*tyt4(i1,i2,i3)
        tyy43(i1,i2,i3)=ry(i1,i2,i3)*tyr4(i1,i2,i3)+sy(i1,i2,i3)*tys4(i1,i2,i3)+ty(i1,i2,i3)*tyt4(i1,i2,i3)
        tyz43(i1,i2,i3)=rz(i1,i2,i3)*tyr4(i1,i2,i3)+sz(i1,i2,i3)*tys4(i1,i2,i3)+tz(i1,i2,i3)*tyt4(i1,i2,i3)
        tzx42(i1,i2,i3)= rx(i1,i2,i3)*tzr4(i1,i2,i3)+sx(i1,i2,i3)*tzs4(i1,i2,i3)
        tzy42(i1,i2,i3)= ry(i1,i2,i3)*tzr4(i1,i2,i3)+sy(i1,i2,i3)*tzs4(i1,i2,i3)
        tzx43(i1,i2,i3)=rx(i1,i2,i3)*tzr4(i1,i2,i3)+sx(i1,i2,i3)*tzs4(i1,i2,i3)+tx(i1,i2,i3)*tzt4(i1,i2,i3)
        tzy43(i1,i2,i3)=ry(i1,i2,i3)*tzr4(i1,i2,i3)+sy(i1,i2,i3)*tzs4(i1,i2,i3)+ty(i1,i2,i3)*tzt4(i1,i2,i3)
        tzz43(i1,i2,i3)=rz(i1,i2,i3)*tzr4(i1,i2,i3)+sz(i1,i2,i3)*tzs4(i1,i2,i3)+tz(i1,i2,i3)*tzt4(i1,i2,i3)
        uxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*ur4(i1,i2,i3,kd)
        uyy41(i1,i2,i3,kd)=0
        uxy41(i1,i2,i3,kd)=0
        uxz41(i1,i2,i3,kd)=0
        uyz41(i1,i2,i3,kd)=0
        uzz41(i1,i2,i3,kd)=0
        ulaplacian41(i1,i2,i3,kd)=uxx41(i1,i2,i3,kd)
        uxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx42(i1,i2,i3))*us4(i1,i2,i3,kd)
        uyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(syy42(i1,i2,i3))*us4(i1,i2,i3,kd)
        uxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss4(i1,i2,i3,kd)+rxy42(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*us4(i1,i2,i3,kd)
        uxz42(i1,i2,i3,kd)=0
        uyz42(i1,i2,i3,kd)=0
        uzz42(i1,i2,i3,kd)=0
        ulaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*us4(i1,i2,i3,kd)
        uxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sx(i1,i2,i3)**2*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust4(i1,i2,i3,kd)+rxx43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxx43(i1,i2,i3)*us4(i1,i2,i3,kd)+txx43(i1,i2,i3)*ut4(i1,i2,i3,kd)
        uyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sy(i1,i2,i3)**2*uss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust4(i1,i2,i3,kd)+ryy43(i1,i2,i3)*ur4(i1,i2,i3,kd)+syy43(i1,i2,i3)*us4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*ut4(i1,i2,i3,kd)
        uzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr4(i1,i2,i3,kd)+sz(i1,i2,i3)**2*uss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*urs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*urt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust4(i1,i2,i3,kd)+rzz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+szz43(i1,i2,i3)*us4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
        uxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*utt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust4(i1,i2,i3,kd)+rxy43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*us4(i1,i2,i3,kd)+txy43(i1,i2,i3)*ut4(i1,i2,i3,kd)
        uxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*uss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*utt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*urs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust4(i1,i2,i3,kd)+rxz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*us4(i1,i2,i3,kd)+txz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
        uyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr4(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*uss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*utt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*urs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*urt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust4(i1,i2,i3,kd)+ryz43(i1,i2,i3)*ur4(i1,i2,i3,kd)+syz43(i1,i2,i3)*us4(i1,i2,i3,kd)+tyz43(i1,i2,i3)*ut4(i1,i2,i3,kd)
        ulaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*urr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*uss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*utt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*urs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*urt4(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*ust4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+ryy43(i1,i2,i3)+rzz43(i1,i2,i3))*ur4(i1,i2,i3,kd)+(sxx43(i1,i2,i3)+syy43(i1,i2,i3)+szz43(i1,i2,i3))*us4(i1,i2,i3,kd)+(txx43(i1,i2,i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*ut4(i1,i2,i3,kd)
    !============================================================================================
    ! Define derivatives for a rectangular grid
    !
    !============================================================================================
        h41(kd) = 1./(12.*dx(kd))
        h42(kd) = 1./(12.*dx(kd)**2)
        ux43r(i1,i2,i3,kd)=(8.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd)))*h41(0)
        uy43r(i1,i2,i3,kd)=(8.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd)))*h41(1)
        uz43r(i1,i2,i3,kd)=(8.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd)))*h41(2)
        uxx43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd)) )*h42(0) 
        uyy43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd)) )*h42(1) 
        uzz43r(i1,i2,i3,kd)=( -30.*u(i1,i2,i3,kd)+16.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd)) )*h42(2)
        uxy43r(i1,i2,i3,kd)=( (u(i1+2,i2+2,i3,kd)-u(i1-2,i2+2,i3,kd)- u(i1+2,i2-2,i3,kd)+u(i1-2,i2-2,i3,kd)) +8.*(u(i1-1,i2+2,i3,kd)-u(i1-1,i2-2,i3,kd)-u(i1+1,i2+2,i3,kd)+u(i1+1,i2-2,i3,kd) +u(i1+2,i2-1,i3,kd)-u(i1-2,i2-1,i3,kd)-u(i1+2,i2+1,i3,kd)+u(i1-2,i2+1,i3,kd))+64.*(u(i1+1,i2+1,i3,kd)-u(i1-1,i2+1,i3,kd)- u(i1+1,i2-1,i3,kd)+u(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
        uxz43r(i1,i2,i3,kd)=( (u(i1+2,i2,i3+2,kd)-u(i1-2,i2,i3+2,kd)-u(i1+2,i2,i3-2,kd)+u(i1-2,i2,i3-2,kd)) +8.*(u(i1-1,i2,i3+2,kd)-u(i1-1,i2,i3-2,kd)-u(i1+1,i2,i3+2,kd)+u(i1+1,i2,i3-2,kd) +u(i1+2,i2,i3-1,kd)-u(i1-2,i2,i3-1,kd)- u(i1+2,i2,i3+1,kd)+u(i1-2,i2,i3+1,kd)) +64.*(u(i1+1,i2,i3+1,kd)-u(i1-1,i2,i3+1,kd)-u(i1+1,i2,i3-1,kd)+u(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
        uyz43r(i1,i2,i3,kd)=( (u(i1,i2+2,i3+2,kd)-u(i1,i2-2,i3+2,kd)-u(i1,i2+2,i3-2,kd)+u(i1,i2-2,i3-2,kd)) +8.*(u(i1,i2-1,i3+2,kd)-u(i1,i2-1,i3-2,kd)-u(i1,i2+1,i3+2,kd)+u(i1,i2+1,i3-2,kd) +u(i1,i2+2,i3-1,kd)-u(i1,i2-2,i3-1,kd)-u(i1,i2+2,i3+1,kd)+u(i1,i2-2,i3+1,kd)) +64.*(u(i1,i2+1,i3+1,kd)-u(i1,i2-1,i3+1,kd)-u(i1,i2+1,i3-1,kd)+u(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
        ux41r(i1,i2,i3,kd)= ux43r(i1,i2,i3,kd)
        uy41r(i1,i2,i3,kd)= uy43r(i1,i2,i3,kd)
        uz41r(i1,i2,i3,kd)= uz43r(i1,i2,i3,kd)
        uxx41r(i1,i2,i3,kd)= uxx43r(i1,i2,i3,kd)
        uyy41r(i1,i2,i3,kd)= uyy43r(i1,i2,i3,kd)
        uzz41r(i1,i2,i3,kd)= uzz43r(i1,i2,i3,kd)
        uxy41r(i1,i2,i3,kd)= uxy43r(i1,i2,i3,kd)
        uxz41r(i1,i2,i3,kd)= uxz43r(i1,i2,i3,kd)
        uyz41r(i1,i2,i3,kd)= uyz43r(i1,i2,i3,kd)
        ulaplacian41r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)
        ux42r(i1,i2,i3,kd)= ux43r(i1,i2,i3,kd)
        uy42r(i1,i2,i3,kd)= uy43r(i1,i2,i3,kd)
        uz42r(i1,i2,i3,kd)= uz43r(i1,i2,i3,kd)
        uxx42r(i1,i2,i3,kd)= uxx43r(i1,i2,i3,kd)
        uyy42r(i1,i2,i3,kd)= uyy43r(i1,i2,i3,kd)
        uzz42r(i1,i2,i3,kd)= uzz43r(i1,i2,i3,kd)
        uxy42r(i1,i2,i3,kd)= uxy43r(i1,i2,i3,kd)
        uxz42r(i1,i2,i3,kd)= uxz43r(i1,i2,i3,kd)
        uyz42r(i1,i2,i3,kd)= uyz43r(i1,i2,i3,kd)
        ulaplacian42r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)+uyy43r(i1,i2,i3,kd)
        ulaplacian43r(i1,i2,i3,kd)=uxx43r(i1,i2,i3,kd)+uyy43r(i1,i2,i3,kd)+uzz43r(i1,i2,i3,kd)
        d16(kd) = 1./(60.*dr(kd))
        d26(kd) = 1./(180.*dr(kd)**2)
        ur6(i1,i2,i3,kd)=(45.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-9.*(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd))+(u(i1+3,i2,i3,kd)-u(i1-3,i2,i3,kd)))*d16(0)
        us6(i1,i2,i3,kd)=(45.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-9.*(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd))+(u(i1,i2+3,i3,kd)-u(i1,i2-3,i3,kd)))*d16(1)
        ut6(i1,i2,i3,kd)=(45.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-9.*(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd))+(u(i1,i2,i3+3,kd)-u(i1,i2,i3-3,kd)))*d16(2)
        urr6(i1,i2,i3,kd)=(-490.*u(i1,i2,i3,kd)+270.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-27.*(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd))+2.*(u(i1+3,i2,i3,kd)+u(i1-3,i2,i3,kd)) )*d26(0)
        uss6(i1,i2,i3,kd)=(-490.*u(i1,i2,i3,kd)+270.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-27.*(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd))+2.*(u(i1,i2+3,i3,kd)+u(i1,i2-3,i3,kd)) )*d26(1)
        utt6(i1,i2,i3,kd)=(-490.*u(i1,i2,i3,kd)+270.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-27.*(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd))+2.*(u(i1,i2,i3+3,kd)+u(i1,i2,i3-3,kd)) )*d26(2)
        urs6(i1,i2,i3,kd)=(45.*(ur6(i1,i2+1,i3,kd)-ur6(i1,i2-1,i3,kd))-9.*(ur6(i1,i2+2,i3,kd)-ur6(i1,i2-2,i3,kd))+(ur6(i1,i2+3,i3,kd)-ur6(i1,i2-3,i3,kd)))*d16(1)
        urt6(i1,i2,i3,kd)=(45.*(ur6(i1,i2,i3+1,kd)-ur6(i1,i2,i3-1,kd))-9.*(ur6(i1,i2,i3+2,kd)-ur6(i1,i2,i3-2,kd))+(ur6(i1,i2,i3+3,kd)-ur6(i1,i2,i3-3,kd)))*d16(2)
        ust6(i1,i2,i3,kd)=(45.*(us6(i1,i2,i3+1,kd)-us6(i1,i2,i3-1,kd))-9.*(us6(i1,i2,i3+2,kd)-us6(i1,i2,i3-2,kd))+(us6(i1,i2,i3+3,kd)-us6(i1,i2,i3-3,kd)))*d16(2)
        rxr6(i1,i2,i3)=(45.*(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))-9.*(rx(i1+2,i2,i3)-rx(i1-2,i2,i3))+(rx(i1+3,i2,i3)-rx(i1-3,i2,i3)))*d16(0)
        rxs6(i1,i2,i3)=(45.*(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))-9.*(rx(i1,i2+2,i3)-rx(i1,i2-2,i3))+(rx(i1,i2+3,i3)-rx(i1,i2-3,i3)))*d16(1)
        rxt6(i1,i2,i3)=(45.*(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))-9.*(rx(i1,i2,i3+2)-rx(i1,i2,i3-2))+(rx(i1,i2,i3+3)-rx(i1,i2,i3-3)))*d16(2)
        ryr6(i1,i2,i3)=(45.*(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))-9.*(ry(i1+2,i2,i3)-ry(i1-2,i2,i3))+(ry(i1+3,i2,i3)-ry(i1-3,i2,i3)))*d16(0)
        rys6(i1,i2,i3)=(45.*(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))-9.*(ry(i1,i2+2,i3)-ry(i1,i2-2,i3))+(ry(i1,i2+3,i3)-ry(i1,i2-3,i3)))*d16(1)
        ryt6(i1,i2,i3)=(45.*(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))-9.*(ry(i1,i2,i3+2)-ry(i1,i2,i3-2))+(ry(i1,i2,i3+3)-ry(i1,i2,i3-3)))*d16(2)
        rzr6(i1,i2,i3)=(45.*(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))-9.*(rz(i1+2,i2,i3)-rz(i1-2,i2,i3))+(rz(i1+3,i2,i3)-rz(i1-3,i2,i3)))*d16(0)
        rzs6(i1,i2,i3)=(45.*(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))-9.*(rz(i1,i2+2,i3)-rz(i1,i2-2,i3))+(rz(i1,i2+3,i3)-rz(i1,i2-3,i3)))*d16(1)
        rzt6(i1,i2,i3)=(45.*(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))-9.*(rz(i1,i2,i3+2)-rz(i1,i2,i3-2))+(rz(i1,i2,i3+3)-rz(i1,i2,i3-3)))*d16(2)
        sxr6(i1,i2,i3)=(45.*(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))-9.*(sx(i1+2,i2,i3)-sx(i1-2,i2,i3))+(sx(i1+3,i2,i3)-sx(i1-3,i2,i3)))*d16(0)
        sxs6(i1,i2,i3)=(45.*(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))-9.*(sx(i1,i2+2,i3)-sx(i1,i2-2,i3))+(sx(i1,i2+3,i3)-sx(i1,i2-3,i3)))*d16(1)
        sxt6(i1,i2,i3)=(45.*(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))-9.*(sx(i1,i2,i3+2)-sx(i1,i2,i3-2))+(sx(i1,i2,i3+3)-sx(i1,i2,i3-3)))*d16(2)
        syr6(i1,i2,i3)=(45.*(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))-9.*(sy(i1+2,i2,i3)-sy(i1-2,i2,i3))+(sy(i1+3,i2,i3)-sy(i1-3,i2,i3)))*d16(0)
        sys6(i1,i2,i3)=(45.*(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))-9.*(sy(i1,i2+2,i3)-sy(i1,i2-2,i3))+(sy(i1,i2+3,i3)-sy(i1,i2-3,i3)))*d16(1)
        syt6(i1,i2,i3)=(45.*(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))-9.*(sy(i1,i2,i3+2)-sy(i1,i2,i3-2))+(sy(i1,i2,i3+3)-sy(i1,i2,i3-3)))*d16(2)
        szr6(i1,i2,i3)=(45.*(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))-9.*(sz(i1+2,i2,i3)-sz(i1-2,i2,i3))+(sz(i1+3,i2,i3)-sz(i1-3,i2,i3)))*d16(0)
        szs6(i1,i2,i3)=(45.*(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))-9.*(sz(i1,i2+2,i3)-sz(i1,i2-2,i3))+(sz(i1,i2+3,i3)-sz(i1,i2-3,i3)))*d16(1)
        szt6(i1,i2,i3)=(45.*(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))-9.*(sz(i1,i2,i3+2)-sz(i1,i2,i3-2))+(sz(i1,i2,i3+3)-sz(i1,i2,i3-3)))*d16(2)
        txr6(i1,i2,i3)=(45.*(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))-9.*(tx(i1+2,i2,i3)-tx(i1-2,i2,i3))+(tx(i1+3,i2,i3)-tx(i1-3,i2,i3)))*d16(0)
        txs6(i1,i2,i3)=(45.*(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))-9.*(tx(i1,i2+2,i3)-tx(i1,i2-2,i3))+(tx(i1,i2+3,i3)-tx(i1,i2-3,i3)))*d16(1)
        txt6(i1,i2,i3)=(45.*(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))-9.*(tx(i1,i2,i3+2)-tx(i1,i2,i3-2))+(tx(i1,i2,i3+3)-tx(i1,i2,i3-3)))*d16(2)
        tyr6(i1,i2,i3)=(45.*(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))-9.*(ty(i1+2,i2,i3)-ty(i1-2,i2,i3))+(ty(i1+3,i2,i3)-ty(i1-3,i2,i3)))*d16(0)
        tys6(i1,i2,i3)=(45.*(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))-9.*(ty(i1,i2+2,i3)-ty(i1,i2-2,i3))+(ty(i1,i2+3,i3)-ty(i1,i2-3,i3)))*d16(1)
        tyt6(i1,i2,i3)=(45.*(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))-9.*(ty(i1,i2,i3+2)-ty(i1,i2,i3-2))+(ty(i1,i2,i3+3)-ty(i1,i2,i3-3)))*d16(2)
        tzr6(i1,i2,i3)=(45.*(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))-9.*(tz(i1+2,i2,i3)-tz(i1-2,i2,i3))+(tz(i1+3,i2,i3)-tz(i1-3,i2,i3)))*d16(0)
        tzs6(i1,i2,i3)=(45.*(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))-9.*(tz(i1,i2+2,i3)-tz(i1,i2-2,i3))+(tz(i1,i2+3,i3)-tz(i1,i2-3,i3)))*d16(1)
        tzt6(i1,i2,i3)=(45.*(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))-9.*(tz(i1,i2,i3+2)-tz(i1,i2,i3-2))+(tz(i1,i2,i3+3)-tz(i1,i2,i3-3)))*d16(2)
        ux61(i1,i2,i3,kd)= rx(i1,i2,i3)*ur6(i1,i2,i3,kd)
        uy61(i1,i2,i3,kd)=0
        uz61(i1,i2,i3,kd)=0
        ux62(i1,i2,i3,kd)= rx(i1,i2,i3)*ur6(i1,i2,i3,kd)+sx(i1,i2,i3)*us6(i1,i2,i3,kd)
        uy62(i1,i2,i3,kd)= ry(i1,i2,i3)*ur6(i1,i2,i3,kd)+sy(i1,i2,i3)*us6(i1,i2,i3,kd)
        uz62(i1,i2,i3,kd)=0
        ux63(i1,i2,i3,kd)=rx(i1,i2,i3)*ur6(i1,i2,i3,kd)+sx(i1,i2,i3)*us6(i1,i2,i3,kd)+tx(i1,i2,i3)*ut6(i1,i2,i3,kd)
        uy63(i1,i2,i3,kd)=ry(i1,i2,i3)*ur6(i1,i2,i3,kd)+sy(i1,i2,i3)*us6(i1,i2,i3,kd)+ty(i1,i2,i3)*ut6(i1,i2,i3,kd)
        uz63(i1,i2,i3,kd)=rz(i1,i2,i3)*ur6(i1,i2,i3,kd)+sz(i1,i2,i3)*us6(i1,i2,i3,kd)+tz(i1,i2,i3)*ut6(i1,i2,i3,kd)
        rxx61(i1,i2,i3)= rx(i1,i2,i3)*rxr6(i1,i2,i3)
        rxx62(i1,i2,i3)= rx(i1,i2,i3)*rxr6(i1,i2,i3)+sx(i1,i2,i3)*rxs6(i1,i2,i3)
        rxy62(i1,i2,i3)= ry(i1,i2,i3)*rxr6(i1,i2,i3)+sy(i1,i2,i3)*rxs6(i1,i2,i3)
        rxx63(i1,i2,i3)=rx(i1,i2,i3)*rxr6(i1,i2,i3)+sx(i1,i2,i3)*rxs6(i1,i2,i3)+tx(i1,i2,i3)*rxt6(i1,i2,i3)
        rxy63(i1,i2,i3)=ry(i1,i2,i3)*rxr6(i1,i2,i3)+sy(i1,i2,i3)*rxs6(i1,i2,i3)+ty(i1,i2,i3)*rxt6(i1,i2,i3)
        rxz63(i1,i2,i3)=rz(i1,i2,i3)*rxr6(i1,i2,i3)+sz(i1,i2,i3)*rxs6(i1,i2,i3)+tz(i1,i2,i3)*rxt6(i1,i2,i3)
        ryx62(i1,i2,i3)= rx(i1,i2,i3)*ryr6(i1,i2,i3)+sx(i1,i2,i3)*rys6(i1,i2,i3)
        ryy62(i1,i2,i3)= ry(i1,i2,i3)*ryr6(i1,i2,i3)+sy(i1,i2,i3)*rys6(i1,i2,i3)
        ryx63(i1,i2,i3)=rx(i1,i2,i3)*ryr6(i1,i2,i3)+sx(i1,i2,i3)*rys6(i1,i2,i3)+tx(i1,i2,i3)*ryt6(i1,i2,i3)
        ryy63(i1,i2,i3)=ry(i1,i2,i3)*ryr6(i1,i2,i3)+sy(i1,i2,i3)*rys6(i1,i2,i3)+ty(i1,i2,i3)*ryt6(i1,i2,i3)
        ryz63(i1,i2,i3)=rz(i1,i2,i3)*ryr6(i1,i2,i3)+sz(i1,i2,i3)*rys6(i1,i2,i3)+tz(i1,i2,i3)*ryt6(i1,i2,i3)
        rzx62(i1,i2,i3)= rx(i1,i2,i3)*rzr6(i1,i2,i3)+sx(i1,i2,i3)*rzs6(i1,i2,i3)
        rzy62(i1,i2,i3)= ry(i1,i2,i3)*rzr6(i1,i2,i3)+sy(i1,i2,i3)*rzs6(i1,i2,i3)
        rzx63(i1,i2,i3)=rx(i1,i2,i3)*rzr6(i1,i2,i3)+sx(i1,i2,i3)*rzs6(i1,i2,i3)+tx(i1,i2,i3)*rzt6(i1,i2,i3)
        rzy63(i1,i2,i3)=ry(i1,i2,i3)*rzr6(i1,i2,i3)+sy(i1,i2,i3)*rzs6(i1,i2,i3)+ty(i1,i2,i3)*rzt6(i1,i2,i3)
        rzz63(i1,i2,i3)=rz(i1,i2,i3)*rzr6(i1,i2,i3)+sz(i1,i2,i3)*rzs6(i1,i2,i3)+tz(i1,i2,i3)*rzt6(i1,i2,i3)
        sxx62(i1,i2,i3)= rx(i1,i2,i3)*sxr6(i1,i2,i3)+sx(i1,i2,i3)*sxs6(i1,i2,i3)
        sxy62(i1,i2,i3)= ry(i1,i2,i3)*sxr6(i1,i2,i3)+sy(i1,i2,i3)*sxs6(i1,i2,i3)
        sxx63(i1,i2,i3)=rx(i1,i2,i3)*sxr6(i1,i2,i3)+sx(i1,i2,i3)*sxs6(i1,i2,i3)+tx(i1,i2,i3)*sxt6(i1,i2,i3)
        sxy63(i1,i2,i3)=ry(i1,i2,i3)*sxr6(i1,i2,i3)+sy(i1,i2,i3)*sxs6(i1,i2,i3)+ty(i1,i2,i3)*sxt6(i1,i2,i3)
        sxz63(i1,i2,i3)=rz(i1,i2,i3)*sxr6(i1,i2,i3)+sz(i1,i2,i3)*sxs6(i1,i2,i3)+tz(i1,i2,i3)*sxt6(i1,i2,i3)
        syx62(i1,i2,i3)= rx(i1,i2,i3)*syr6(i1,i2,i3)+sx(i1,i2,i3)*sys6(i1,i2,i3)
        syy62(i1,i2,i3)= ry(i1,i2,i3)*syr6(i1,i2,i3)+sy(i1,i2,i3)*sys6(i1,i2,i3)
        syx63(i1,i2,i3)=rx(i1,i2,i3)*syr6(i1,i2,i3)+sx(i1,i2,i3)*sys6(i1,i2,i3)+tx(i1,i2,i3)*syt6(i1,i2,i3)
        syy63(i1,i2,i3)=ry(i1,i2,i3)*syr6(i1,i2,i3)+sy(i1,i2,i3)*sys6(i1,i2,i3)+ty(i1,i2,i3)*syt6(i1,i2,i3)
        syz63(i1,i2,i3)=rz(i1,i2,i3)*syr6(i1,i2,i3)+sz(i1,i2,i3)*sys6(i1,i2,i3)+tz(i1,i2,i3)*syt6(i1,i2,i3)
        szx62(i1,i2,i3)= rx(i1,i2,i3)*szr6(i1,i2,i3)+sx(i1,i2,i3)*szs6(i1,i2,i3)
        szy62(i1,i2,i3)= ry(i1,i2,i3)*szr6(i1,i2,i3)+sy(i1,i2,i3)*szs6(i1,i2,i3)
        szx63(i1,i2,i3)=rx(i1,i2,i3)*szr6(i1,i2,i3)+sx(i1,i2,i3)*szs6(i1,i2,i3)+tx(i1,i2,i3)*szt6(i1,i2,i3)
        szy63(i1,i2,i3)=ry(i1,i2,i3)*szr6(i1,i2,i3)+sy(i1,i2,i3)*szs6(i1,i2,i3)+ty(i1,i2,i3)*szt6(i1,i2,i3)
        szz63(i1,i2,i3)=rz(i1,i2,i3)*szr6(i1,i2,i3)+sz(i1,i2,i3)*szs6(i1,i2,i3)+tz(i1,i2,i3)*szt6(i1,i2,i3)
        txx62(i1,i2,i3)= rx(i1,i2,i3)*txr6(i1,i2,i3)+sx(i1,i2,i3)*txs6(i1,i2,i3)
        txy62(i1,i2,i3)= ry(i1,i2,i3)*txr6(i1,i2,i3)+sy(i1,i2,i3)*txs6(i1,i2,i3)
        txx63(i1,i2,i3)=rx(i1,i2,i3)*txr6(i1,i2,i3)+sx(i1,i2,i3)*txs6(i1,i2,i3)+tx(i1,i2,i3)*txt6(i1,i2,i3)
        txy63(i1,i2,i3)=ry(i1,i2,i3)*txr6(i1,i2,i3)+sy(i1,i2,i3)*txs6(i1,i2,i3)+ty(i1,i2,i3)*txt6(i1,i2,i3)
        txz63(i1,i2,i3)=rz(i1,i2,i3)*txr6(i1,i2,i3)+sz(i1,i2,i3)*txs6(i1,i2,i3)+tz(i1,i2,i3)*txt6(i1,i2,i3)
        tyx62(i1,i2,i3)= rx(i1,i2,i3)*tyr6(i1,i2,i3)+sx(i1,i2,i3)*tys6(i1,i2,i3)
        tyy62(i1,i2,i3)= ry(i1,i2,i3)*tyr6(i1,i2,i3)+sy(i1,i2,i3)*tys6(i1,i2,i3)
        tyx63(i1,i2,i3)=rx(i1,i2,i3)*tyr6(i1,i2,i3)+sx(i1,i2,i3)*tys6(i1,i2,i3)+tx(i1,i2,i3)*tyt6(i1,i2,i3)
        tyy63(i1,i2,i3)=ry(i1,i2,i3)*tyr6(i1,i2,i3)+sy(i1,i2,i3)*tys6(i1,i2,i3)+ty(i1,i2,i3)*tyt6(i1,i2,i3)
        tyz63(i1,i2,i3)=rz(i1,i2,i3)*tyr6(i1,i2,i3)+sz(i1,i2,i3)*tys6(i1,i2,i3)+tz(i1,i2,i3)*tyt6(i1,i2,i3)
        tzx62(i1,i2,i3)= rx(i1,i2,i3)*tzr6(i1,i2,i3)+sx(i1,i2,i3)*tzs6(i1,i2,i3)
        tzy62(i1,i2,i3)= ry(i1,i2,i3)*tzr6(i1,i2,i3)+sy(i1,i2,i3)*tzs6(i1,i2,i3)
        tzx63(i1,i2,i3)=rx(i1,i2,i3)*tzr6(i1,i2,i3)+sx(i1,i2,i3)*tzs6(i1,i2,i3)+tx(i1,i2,i3)*tzt6(i1,i2,i3)
        tzy63(i1,i2,i3)=ry(i1,i2,i3)*tzr6(i1,i2,i3)+sy(i1,i2,i3)*tzs6(i1,i2,i3)+ty(i1,i2,i3)*tzt6(i1,i2,i3)
        tzz63(i1,i2,i3)=rz(i1,i2,i3)*tzr6(i1,i2,i3)+sz(i1,i2,i3)*tzs6(i1,i2,i3)+tz(i1,i2,i3)*tzt6(i1,i2,i3)
        uxx61(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr6(i1,i2,i3,kd)+(rxx62(i1,i2,i3))*ur6(i1,i2,i3,kd)
        uyy61(i1,i2,i3,kd)=0
        uxy61(i1,i2,i3,kd)=0
        uxz61(i1,i2,i3,kd)=0
        uyz61(i1,i2,i3,kd)=0
        uzz61(i1,i2,i3,kd)=0
        ulaplacian61(i1,i2,i3,kd)=uxx61(i1,i2,i3,kd)
        uxx62(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr6(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*urs6(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*uss6(i1,i2,i3,kd)+(rxx62(i1,i2,i3))*ur6(i1,i2,i3,kd)+(sxx62(i1,i2,i3))*us6(i1,i2,i3,kd)
        uyy62(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr6(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*urs6(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*uss6(i1,i2,i3,kd)+(ryy62(i1,i2,i3))*ur6(i1,i2,i3,kd)+(syy62(i1,i2,i3))*us6(i1,i2,i3,kd)
        uxy62(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr6(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs6(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss6(i1,i2,i3,kd)+rxy62(i1,i2,i3)*ur6(i1,i2,i3,kd)+sxy62(i1,i2,i3)*us6(i1,i2,i3,kd)
        uxz62(i1,i2,i3,kd)=0
        uyz62(i1,i2,i3,kd)=0
        uzz62(i1,i2,i3,kd)=0
        ulaplacian62(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*urr6(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*urs6(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*uss6(i1,i2,i3,kd)+(rxx62(i1,i2,i3)+ryy62(i1,i2,i3))*ur6(i1,i2,i3,kd)+(sxx62(i1,i2,i3)+syy62(i1,i2,i3))*us6(i1,i2,i3,kd)
        uxx63(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr6(i1,i2,i3,kd)+sx(i1,i2,i3)**2*uss6(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt6(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*urs6(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*urt6(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust6(i1,i2,i3,kd)+rxx63(i1,i2,i3)*ur6(i1,i2,i3,kd)+sxx63(i1,i2,i3)*us6(i1,i2,i3,kd)+txx63(i1,i2,i3)*ut6(i1,i2,i3,kd)
        uyy63(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr6(i1,i2,i3,kd)+sy(i1,i2,i3)**2*uss6(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt6(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*urs6(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*urt6(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust6(i1,i2,i3,kd)+ryy63(i1,i2,i3)*ur6(i1,i2,i3,kd)+syy63(i1,i2,i3)*us6(i1,i2,i3,kd)+tyy63(i1,i2,i3)*ut6(i1,i2,i3,kd)
        uzz63(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr6(i1,i2,i3,kd)+sz(i1,i2,i3)**2*uss6(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt6(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*urs6(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*urt6(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust6(i1,i2,i3,kd)+rzz63(i1,i2,i3)*ur6(i1,i2,i3,kd)+szz63(i1,i2,i3)*us6(i1,i2,i3,kd)+tzz63(i1,i2,i3)*ut6(i1,i2,i3,kd)
        uxy63(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr6(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss6(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*utt6(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs6(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*urt6(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust6(i1,i2,i3,kd)+rxy63(i1,i2,i3)*ur6(i1,i2,i3,kd)+sxy63(i1,i2,i3)*us6(i1,i2,i3,kd)+txy63(i1,i2,i3)*ut6(i1,i2,i3,kd)
        uxz63(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr6(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*uss6(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*utt6(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*urs6(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*urt6(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust6(i1,i2,i3,kd)+rxz63(i1,i2,i3)*ur6(i1,i2,i3,kd)+sxz63(i1,i2,i3)*us6(i1,i2,i3,kd)+txz63(i1,i2,i3)*ut6(i1,i2,i3,kd)
        uyz63(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr6(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*uss6(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*utt6(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*urs6(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*urt6(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust6(i1,i2,i3,kd)+ryz63(i1,i2,i3)*ur6(i1,i2,i3,kd)+syz63(i1,i2,i3)*us6(i1,i2,i3,kd)+tyz63(i1,i2,i3)*ut6(i1,i2,i3,kd)
        ulaplacian63(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*urr6(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*uss6(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*utt6(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*urs6(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*urt6(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*ust6(i1,i2,i3,kd)+(rxx63(i1,i2,i3)+ryy63(i1,i2,i3)+rzz63(i1,i2,i3))*ur6(i1,i2,i3,kd)+(sxx63(i1,i2,i3)+syy63(i1,i2,i3)+szz63(i1,i2,i3))*us6(i1,i2,i3,kd)+(txx63(i1,i2,i3)+tyy63(i1,i2,i3)+tzz63(i1,i2,i3))*ut6(i1,i2,i3,kd)
    !============================================================================================
    ! Define derivatives for a rectangular grid
    !
    !============================================================================================
        h16(kd) = 1./(60.*dx(kd))
        h26(kd) = 1./(180.*dx(kd)**2)
        ux63r(i1,i2,i3,kd)=(45.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-9.*(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd))+(u(i1+3,i2,i3,kd)-u(i1-3,i2,i3,kd)))*h16(0)
        uy63r(i1,i2,i3,kd)=(45.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-9.*(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd))+(u(i1,i2+3,i3,kd)-u(i1,i2-3,i3,kd)))*h16(1)
        uz63r(i1,i2,i3,kd)=(45.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-9.*(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd))+(u(i1,i2,i3+3,kd)-u(i1,i2,i3-3,kd)))*h16(2)
        uxx63r(i1,i2,i3,kd)=(-490.*u(i1,i2,i3,kd)+270.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-27.*(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd))+2.*(u(i1+3,i2,i3,kd)+u(i1-3,i2,i3,kd)) )*h26(0)
        uyy63r(i1,i2,i3,kd)=(-490.*u(i1,i2,i3,kd)+270.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-27.*(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd))+2.*(u(i1,i2+3,i3,kd)+u(i1,i2-3,i3,kd)) )*h26(1)
        uzz63r(i1,i2,i3,kd)=(-490.*u(i1,i2,i3,kd)+270.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-27.*(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd))+2.*(u(i1,i2,i3+3,kd)+u(i1,i2,i3-3,kd)) )*h26(2)
        uxy63r(i1,i2,i3,kd)=(45.*(ux63r(i1,i2+1,i3,kd)-ux63r(i1,i2-1,i3,kd))-9.*(ux63r(i1,i2+2,i3,kd)-ux63r(i1,i2-2,i3,kd))+(ux63r(i1,i2+3,i3,kd)-ux63r(i1,i2-3,i3,kd)))*h16(1)
        uxz63r(i1,i2,i3,kd)=(45.*(ux63r(i1,i2,i3+1,kd)-ux63r(i1,i2,i3-1,kd))-9.*(ux63r(i1,i2,i3+2,kd)-ux63r(i1,i2,i3-2,kd))+(ux63r(i1,i2,i3+3,kd)-ux63r(i1,i2,i3-3,kd)))*h16(2)
        uyz63r(i1,i2,i3,kd)=(45.*(uy63r(i1,i2,i3+1,kd)-uy63r(i1,i2,i3-1,kd))-9.*(uy63r(i1,i2,i3+2,kd)-uy63r(i1,i2,i3-2,kd))+(uy63r(i1,i2,i3+3,kd)-uy63r(i1,i2,i3-3,kd)))*h16(2)
        ux61r(i1,i2,i3,kd)= ux63r(i1,i2,i3,kd)
        uy61r(i1,i2,i3,kd)= uy63r(i1,i2,i3,kd)
        uz61r(i1,i2,i3,kd)= uz63r(i1,i2,i3,kd)
        uxx61r(i1,i2,i3,kd)= uxx63r(i1,i2,i3,kd)
        uyy61r(i1,i2,i3,kd)= uyy63r(i1,i2,i3,kd)
        uzz61r(i1,i2,i3,kd)= uzz63r(i1,i2,i3,kd)
        uxy61r(i1,i2,i3,kd)= uxy63r(i1,i2,i3,kd)
        uxz61r(i1,i2,i3,kd)= uxz63r(i1,i2,i3,kd)
        uyz61r(i1,i2,i3,kd)= uyz63r(i1,i2,i3,kd)
        ulaplacian61r(i1,i2,i3,kd)=uxx63r(i1,i2,i3,kd)
        ux62r(i1,i2,i3,kd)= ux63r(i1,i2,i3,kd)
        uy62r(i1,i2,i3,kd)= uy63r(i1,i2,i3,kd)
        uz62r(i1,i2,i3,kd)= uz63r(i1,i2,i3,kd)
        uxx62r(i1,i2,i3,kd)= uxx63r(i1,i2,i3,kd)
        uyy62r(i1,i2,i3,kd)= uyy63r(i1,i2,i3,kd)
        uzz62r(i1,i2,i3,kd)= uzz63r(i1,i2,i3,kd)
        uxy62r(i1,i2,i3,kd)= uxy63r(i1,i2,i3,kd)
        uxz62r(i1,i2,i3,kd)= uxz63r(i1,i2,i3,kd)
        uyz62r(i1,i2,i3,kd)= uyz63r(i1,i2,i3,kd)
        ulaplacian62r(i1,i2,i3,kd)=uxx63r(i1,i2,i3,kd)+uyy63r(i1,i2,i3,kd)
        ulaplacian63r(i1,i2,i3,kd)=uxx63r(i1,i2,i3,kd)+uyy63r(i1,i2,i3,kd)+uzz63r(i1,i2,i3,kd)
        d18(kd) = 1./(840.*dr(kd))
        d28(kd) = 1./(5040.*dr(kd)**2)
        ur8(i1,i2,i3,kd)=(672.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-168.*(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd))+32.*(u(i1+3,i2,i3,kd)-u(i1-3,i2,i3,kd))-3.*(u(i1+4,i2,i3,kd)-u(i1-4,i2,i3,kd)))*d18(0)
        us8(i1,i2,i3,kd)=(672.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-168.*(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd))+32.*(u(i1,i2+3,i3,kd)-u(i1,i2-3,i3,kd))-3.*(u(i1,i2+4,i3,kd)-u(i1,i2-4,i3,kd)))*d18(1)
        ut8(i1,i2,i3,kd)=(672.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-168.*(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd))+32.*(u(i1,i2,i3+3,kd)-u(i1,i2,i3-3,kd))-3.*(u(i1,i2,i3+4,kd)-u(i1,i2,i3-4,kd)))*d18(2)
        urr8(i1,i2,i3,kd)=(-14350.*u(i1,i2,i3,kd)+8064.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-1008.*(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd))+128.*(u(i1+3,i2,i3,kd)+u(i1-3,i2,i3,kd))-9.*(u(i1+4,i2,i3,kd)+u(i1-4,i2,i3,kd)) )*d28(0)
        uss8(i1,i2,i3,kd)=(-14350.*u(i1,i2,i3,kd)+8064.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-1008.*(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd))+128.*(u(i1,i2+3,i3,kd)+u(i1,i2-3,i3,kd))-9.*(u(i1,i2+4,i3,kd)+u(i1,i2-4,i3,kd)) )*d28(1)
        utt8(i1,i2,i3,kd)=(-14350.*u(i1,i2,i3,kd)+8064.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-1008.*(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd))+128.*(u(i1,i2,i3+3,kd)+u(i1,i2,i3-3,kd))-9.*(u(i1,i2,i3+4,kd)+u(i1,i2,i3-4,kd)) )*d28(2)
        urs8(i1,i2,i3,kd)=(672.*(ur8(i1,i2+1,i3,kd)-ur8(i1,i2-1,i3,kd))-168.*(ur8(i1,i2+2,i3,kd)-ur8(i1,i2-2,i3,kd))+32.*(ur8(i1,i2+3,i3,kd)-ur8(i1,i2-3,i3,kd))-3.*(ur8(i1,i2+4,i3,kd)-ur8(i1,i2-4,i3,kd)))*d18(1)
        urt8(i1,i2,i3,kd)=(672.*(ur8(i1,i2,i3+1,kd)-ur8(i1,i2,i3-1,kd))-168.*(ur8(i1,i2,i3+2,kd)-ur8(i1,i2,i3-2,kd))+32.*(ur8(i1,i2,i3+3,kd)-ur8(i1,i2,i3-3,kd))-3.*(ur8(i1,i2,i3+4,kd)-ur8(i1,i2,i3-4,kd)))*d18(2)
        ust8(i1,i2,i3,kd)=(672.*(us8(i1,i2,i3+1,kd)-us8(i1,i2,i3-1,kd))-168.*(us8(i1,i2,i3+2,kd)-us8(i1,i2,i3-2,kd))+32.*(us8(i1,i2,i3+3,kd)-us8(i1,i2,i3-3,kd))-3.*(us8(i1,i2,i3+4,kd)-us8(i1,i2,i3-4,kd)))*d18(2)
        rxr8(i1,i2,i3)=(672.*(rx(i1+1,i2,i3)-rx(i1-1,i2,i3))-168.*(rx(i1+2,i2,i3)-rx(i1-2,i2,i3))+32.*(rx(i1+3,i2,i3)-rx(i1-3,i2,i3))-3.*(rx(i1+4,i2,i3)-rx(i1-4,i2,i3)))*d18(0)
        rxs8(i1,i2,i3)=(672.*(rx(i1,i2+1,i3)-rx(i1,i2-1,i3))-168.*(rx(i1,i2+2,i3)-rx(i1,i2-2,i3))+32.*(rx(i1,i2+3,i3)-rx(i1,i2-3,i3))-3.*(rx(i1,i2+4,i3)-rx(i1,i2-4,i3)))*d18(1)
        rxt8(i1,i2,i3)=(672.*(rx(i1,i2,i3+1)-rx(i1,i2,i3-1))-168.*(rx(i1,i2,i3+2)-rx(i1,i2,i3-2))+32.*(rx(i1,i2,i3+3)-rx(i1,i2,i3-3))-3.*(rx(i1,i2,i3+4)-rx(i1,i2,i3-4)))*d18(2)
        ryr8(i1,i2,i3)=(672.*(ry(i1+1,i2,i3)-ry(i1-1,i2,i3))-168.*(ry(i1+2,i2,i3)-ry(i1-2,i2,i3))+32.*(ry(i1+3,i2,i3)-ry(i1-3,i2,i3))-3.*(ry(i1+4,i2,i3)-ry(i1-4,i2,i3)))*d18(0)
        rys8(i1,i2,i3)=(672.*(ry(i1,i2+1,i3)-ry(i1,i2-1,i3))-168.*(ry(i1,i2+2,i3)-ry(i1,i2-2,i3))+32.*(ry(i1,i2+3,i3)-ry(i1,i2-3,i3))-3.*(ry(i1,i2+4,i3)-ry(i1,i2-4,i3)))*d18(1)
        ryt8(i1,i2,i3)=(672.*(ry(i1,i2,i3+1)-ry(i1,i2,i3-1))-168.*(ry(i1,i2,i3+2)-ry(i1,i2,i3-2))+32.*(ry(i1,i2,i3+3)-ry(i1,i2,i3-3))-3.*(ry(i1,i2,i3+4)-ry(i1,i2,i3-4)))*d18(2)
        rzr8(i1,i2,i3)=(672.*(rz(i1+1,i2,i3)-rz(i1-1,i2,i3))-168.*(rz(i1+2,i2,i3)-rz(i1-2,i2,i3))+32.*(rz(i1+3,i2,i3)-rz(i1-3,i2,i3))-3.*(rz(i1+4,i2,i3)-rz(i1-4,i2,i3)))*d18(0)
        rzs8(i1,i2,i3)=(672.*(rz(i1,i2+1,i3)-rz(i1,i2-1,i3))-168.*(rz(i1,i2+2,i3)-rz(i1,i2-2,i3))+32.*(rz(i1,i2+3,i3)-rz(i1,i2-3,i3))-3.*(rz(i1,i2+4,i3)-rz(i1,i2-4,i3)))*d18(1)
        rzt8(i1,i2,i3)=(672.*(rz(i1,i2,i3+1)-rz(i1,i2,i3-1))-168.*(rz(i1,i2,i3+2)-rz(i1,i2,i3-2))+32.*(rz(i1,i2,i3+3)-rz(i1,i2,i3-3))-3.*(rz(i1,i2,i3+4)-rz(i1,i2,i3-4)))*d18(2)
        sxr8(i1,i2,i3)=(672.*(sx(i1+1,i2,i3)-sx(i1-1,i2,i3))-168.*(sx(i1+2,i2,i3)-sx(i1-2,i2,i3))+32.*(sx(i1+3,i2,i3)-sx(i1-3,i2,i3))-3.*(sx(i1+4,i2,i3)-sx(i1-4,i2,i3)))*d18(0)
        sxs8(i1,i2,i3)=(672.*(sx(i1,i2+1,i3)-sx(i1,i2-1,i3))-168.*(sx(i1,i2+2,i3)-sx(i1,i2-2,i3))+32.*(sx(i1,i2+3,i3)-sx(i1,i2-3,i3))-3.*(sx(i1,i2+4,i3)-sx(i1,i2-4,i3)))*d18(1)
        sxt8(i1,i2,i3)=(672.*(sx(i1,i2,i3+1)-sx(i1,i2,i3-1))-168.*(sx(i1,i2,i3+2)-sx(i1,i2,i3-2))+32.*(sx(i1,i2,i3+3)-sx(i1,i2,i3-3))-3.*(sx(i1,i2,i3+4)-sx(i1,i2,i3-4)))*d18(2)
        syr8(i1,i2,i3)=(672.*(sy(i1+1,i2,i3)-sy(i1-1,i2,i3))-168.*(sy(i1+2,i2,i3)-sy(i1-2,i2,i3))+32.*(sy(i1+3,i2,i3)-sy(i1-3,i2,i3))-3.*(sy(i1+4,i2,i3)-sy(i1-4,i2,i3)))*d18(0)
        sys8(i1,i2,i3)=(672.*(sy(i1,i2+1,i3)-sy(i1,i2-1,i3))-168.*(sy(i1,i2+2,i3)-sy(i1,i2-2,i3))+32.*(sy(i1,i2+3,i3)-sy(i1,i2-3,i3))-3.*(sy(i1,i2+4,i3)-sy(i1,i2-4,i3)))*d18(1)
        syt8(i1,i2,i3)=(672.*(sy(i1,i2,i3+1)-sy(i1,i2,i3-1))-168.*(sy(i1,i2,i3+2)-sy(i1,i2,i3-2))+32.*(sy(i1,i2,i3+3)-sy(i1,i2,i3-3))-3.*(sy(i1,i2,i3+4)-sy(i1,i2,i3-4)))*d18(2)
        szr8(i1,i2,i3)=(672.*(sz(i1+1,i2,i3)-sz(i1-1,i2,i3))-168.*(sz(i1+2,i2,i3)-sz(i1-2,i2,i3))+32.*(sz(i1+3,i2,i3)-sz(i1-3,i2,i3))-3.*(sz(i1+4,i2,i3)-sz(i1-4,i2,i3)))*d18(0)
        szs8(i1,i2,i3)=(672.*(sz(i1,i2+1,i3)-sz(i1,i2-1,i3))-168.*(sz(i1,i2+2,i3)-sz(i1,i2-2,i3))+32.*(sz(i1,i2+3,i3)-sz(i1,i2-3,i3))-3.*(sz(i1,i2+4,i3)-sz(i1,i2-4,i3)))*d18(1)
        szt8(i1,i2,i3)=(672.*(sz(i1,i2,i3+1)-sz(i1,i2,i3-1))-168.*(sz(i1,i2,i3+2)-sz(i1,i2,i3-2))+32.*(sz(i1,i2,i3+3)-sz(i1,i2,i3-3))-3.*(sz(i1,i2,i3+4)-sz(i1,i2,i3-4)))*d18(2)
        txr8(i1,i2,i3)=(672.*(tx(i1+1,i2,i3)-tx(i1-1,i2,i3))-168.*(tx(i1+2,i2,i3)-tx(i1-2,i2,i3))+32.*(tx(i1+3,i2,i3)-tx(i1-3,i2,i3))-3.*(tx(i1+4,i2,i3)-tx(i1-4,i2,i3)))*d18(0)
        txs8(i1,i2,i3)=(672.*(tx(i1,i2+1,i3)-tx(i1,i2-1,i3))-168.*(tx(i1,i2+2,i3)-tx(i1,i2-2,i3))+32.*(tx(i1,i2+3,i3)-tx(i1,i2-3,i3))-3.*(tx(i1,i2+4,i3)-tx(i1,i2-4,i3)))*d18(1)
        txt8(i1,i2,i3)=(672.*(tx(i1,i2,i3+1)-tx(i1,i2,i3-1))-168.*(tx(i1,i2,i3+2)-tx(i1,i2,i3-2))+32.*(tx(i1,i2,i3+3)-tx(i1,i2,i3-3))-3.*(tx(i1,i2,i3+4)-tx(i1,i2,i3-4)))*d18(2)
        tyr8(i1,i2,i3)=(672.*(ty(i1+1,i2,i3)-ty(i1-1,i2,i3))-168.*(ty(i1+2,i2,i3)-ty(i1-2,i2,i3))+32.*(ty(i1+3,i2,i3)-ty(i1-3,i2,i3))-3.*(ty(i1+4,i2,i3)-ty(i1-4,i2,i3)))*d18(0)
        tys8(i1,i2,i3)=(672.*(ty(i1,i2+1,i3)-ty(i1,i2-1,i3))-168.*(ty(i1,i2+2,i3)-ty(i1,i2-2,i3))+32.*(ty(i1,i2+3,i3)-ty(i1,i2-3,i3))-3.*(ty(i1,i2+4,i3)-ty(i1,i2-4,i3)))*d18(1)
        tyt8(i1,i2,i3)=(672.*(ty(i1,i2,i3+1)-ty(i1,i2,i3-1))-168.*(ty(i1,i2,i3+2)-ty(i1,i2,i3-2))+32.*(ty(i1,i2,i3+3)-ty(i1,i2,i3-3))-3.*(ty(i1,i2,i3+4)-ty(i1,i2,i3-4)))*d18(2)
        tzr8(i1,i2,i3)=(672.*(tz(i1+1,i2,i3)-tz(i1-1,i2,i3))-168.*(tz(i1+2,i2,i3)-tz(i1-2,i2,i3))+32.*(tz(i1+3,i2,i3)-tz(i1-3,i2,i3))-3.*(tz(i1+4,i2,i3)-tz(i1-4,i2,i3)))*d18(0)
        tzs8(i1,i2,i3)=(672.*(tz(i1,i2+1,i3)-tz(i1,i2-1,i3))-168.*(tz(i1,i2+2,i3)-tz(i1,i2-2,i3))+32.*(tz(i1,i2+3,i3)-tz(i1,i2-3,i3))-3.*(tz(i1,i2+4,i3)-tz(i1,i2-4,i3)))*d18(1)
        tzt8(i1,i2,i3)=(672.*(tz(i1,i2,i3+1)-tz(i1,i2,i3-1))-168.*(tz(i1,i2,i3+2)-tz(i1,i2,i3-2))+32.*(tz(i1,i2,i3+3)-tz(i1,i2,i3-3))-3.*(tz(i1,i2,i3+4)-tz(i1,i2,i3-4)))*d18(2)
        ux81(i1,i2,i3,kd)= rx(i1,i2,i3)*ur8(i1,i2,i3,kd)
        uy81(i1,i2,i3,kd)=0
        uz81(i1,i2,i3,kd)=0
        ux82(i1,i2,i3,kd)= rx(i1,i2,i3)*ur8(i1,i2,i3,kd)+sx(i1,i2,i3)*us8(i1,i2,i3,kd)
        uy82(i1,i2,i3,kd)= ry(i1,i2,i3)*ur8(i1,i2,i3,kd)+sy(i1,i2,i3)*us8(i1,i2,i3,kd)
        uz82(i1,i2,i3,kd)=0
        ux83(i1,i2,i3,kd)=rx(i1,i2,i3)*ur8(i1,i2,i3,kd)+sx(i1,i2,i3)*us8(i1,i2,i3,kd)+tx(i1,i2,i3)*ut8(i1,i2,i3,kd)
        uy83(i1,i2,i3,kd)=ry(i1,i2,i3)*ur8(i1,i2,i3,kd)+sy(i1,i2,i3)*us8(i1,i2,i3,kd)+ty(i1,i2,i3)*ut8(i1,i2,i3,kd)
        uz83(i1,i2,i3,kd)=rz(i1,i2,i3)*ur8(i1,i2,i3,kd)+sz(i1,i2,i3)*us8(i1,i2,i3,kd)+tz(i1,i2,i3)*ut8(i1,i2,i3,kd)
        rxx81(i1,i2,i3)= rx(i1,i2,i3)*rxr8(i1,i2,i3)
        rxx82(i1,i2,i3)= rx(i1,i2,i3)*rxr8(i1,i2,i3)+sx(i1,i2,i3)*rxs8(i1,i2,i3)
        rxy82(i1,i2,i3)= ry(i1,i2,i3)*rxr8(i1,i2,i3)+sy(i1,i2,i3)*rxs8(i1,i2,i3)
        rxx83(i1,i2,i3)=rx(i1,i2,i3)*rxr8(i1,i2,i3)+sx(i1,i2,i3)*rxs8(i1,i2,i3)+tx(i1,i2,i3)*rxt8(i1,i2,i3)
        rxy83(i1,i2,i3)=ry(i1,i2,i3)*rxr8(i1,i2,i3)+sy(i1,i2,i3)*rxs8(i1,i2,i3)+ty(i1,i2,i3)*rxt8(i1,i2,i3)
        rxz83(i1,i2,i3)=rz(i1,i2,i3)*rxr8(i1,i2,i3)+sz(i1,i2,i3)*rxs8(i1,i2,i3)+tz(i1,i2,i3)*rxt8(i1,i2,i3)
        ryx82(i1,i2,i3)= rx(i1,i2,i3)*ryr8(i1,i2,i3)+sx(i1,i2,i3)*rys8(i1,i2,i3)
        ryy82(i1,i2,i3)= ry(i1,i2,i3)*ryr8(i1,i2,i3)+sy(i1,i2,i3)*rys8(i1,i2,i3)
        ryx83(i1,i2,i3)=rx(i1,i2,i3)*ryr8(i1,i2,i3)+sx(i1,i2,i3)*rys8(i1,i2,i3)+tx(i1,i2,i3)*ryt8(i1,i2,i3)
        ryy83(i1,i2,i3)=ry(i1,i2,i3)*ryr8(i1,i2,i3)+sy(i1,i2,i3)*rys8(i1,i2,i3)+ty(i1,i2,i3)*ryt8(i1,i2,i3)
        ryz83(i1,i2,i3)=rz(i1,i2,i3)*ryr8(i1,i2,i3)+sz(i1,i2,i3)*rys8(i1,i2,i3)+tz(i1,i2,i3)*ryt8(i1,i2,i3)
        rzx82(i1,i2,i3)= rx(i1,i2,i3)*rzr8(i1,i2,i3)+sx(i1,i2,i3)*rzs8(i1,i2,i3)
        rzy82(i1,i2,i3)= ry(i1,i2,i3)*rzr8(i1,i2,i3)+sy(i1,i2,i3)*rzs8(i1,i2,i3)
        rzx83(i1,i2,i3)=rx(i1,i2,i3)*rzr8(i1,i2,i3)+sx(i1,i2,i3)*rzs8(i1,i2,i3)+tx(i1,i2,i3)*rzt8(i1,i2,i3)
        rzy83(i1,i2,i3)=ry(i1,i2,i3)*rzr8(i1,i2,i3)+sy(i1,i2,i3)*rzs8(i1,i2,i3)+ty(i1,i2,i3)*rzt8(i1,i2,i3)
        rzz83(i1,i2,i3)=rz(i1,i2,i3)*rzr8(i1,i2,i3)+sz(i1,i2,i3)*rzs8(i1,i2,i3)+tz(i1,i2,i3)*rzt8(i1,i2,i3)
        sxx82(i1,i2,i3)= rx(i1,i2,i3)*sxr8(i1,i2,i3)+sx(i1,i2,i3)*sxs8(i1,i2,i3)
        sxy82(i1,i2,i3)= ry(i1,i2,i3)*sxr8(i1,i2,i3)+sy(i1,i2,i3)*sxs8(i1,i2,i3)
        sxx83(i1,i2,i3)=rx(i1,i2,i3)*sxr8(i1,i2,i3)+sx(i1,i2,i3)*sxs8(i1,i2,i3)+tx(i1,i2,i3)*sxt8(i1,i2,i3)
        sxy83(i1,i2,i3)=ry(i1,i2,i3)*sxr8(i1,i2,i3)+sy(i1,i2,i3)*sxs8(i1,i2,i3)+ty(i1,i2,i3)*sxt8(i1,i2,i3)
        sxz83(i1,i2,i3)=rz(i1,i2,i3)*sxr8(i1,i2,i3)+sz(i1,i2,i3)*sxs8(i1,i2,i3)+tz(i1,i2,i3)*sxt8(i1,i2,i3)
        syx82(i1,i2,i3)= rx(i1,i2,i3)*syr8(i1,i2,i3)+sx(i1,i2,i3)*sys8(i1,i2,i3)
        syy82(i1,i2,i3)= ry(i1,i2,i3)*syr8(i1,i2,i3)+sy(i1,i2,i3)*sys8(i1,i2,i3)
        syx83(i1,i2,i3)=rx(i1,i2,i3)*syr8(i1,i2,i3)+sx(i1,i2,i3)*sys8(i1,i2,i3)+tx(i1,i2,i3)*syt8(i1,i2,i3)
        syy83(i1,i2,i3)=ry(i1,i2,i3)*syr8(i1,i2,i3)+sy(i1,i2,i3)*sys8(i1,i2,i3)+ty(i1,i2,i3)*syt8(i1,i2,i3)
        syz83(i1,i2,i3)=rz(i1,i2,i3)*syr8(i1,i2,i3)+sz(i1,i2,i3)*sys8(i1,i2,i3)+tz(i1,i2,i3)*syt8(i1,i2,i3)
        szx82(i1,i2,i3)= rx(i1,i2,i3)*szr8(i1,i2,i3)+sx(i1,i2,i3)*szs8(i1,i2,i3)
        szy82(i1,i2,i3)= ry(i1,i2,i3)*szr8(i1,i2,i3)+sy(i1,i2,i3)*szs8(i1,i2,i3)
        szx83(i1,i2,i3)=rx(i1,i2,i3)*szr8(i1,i2,i3)+sx(i1,i2,i3)*szs8(i1,i2,i3)+tx(i1,i2,i3)*szt8(i1,i2,i3)
        szy83(i1,i2,i3)=ry(i1,i2,i3)*szr8(i1,i2,i3)+sy(i1,i2,i3)*szs8(i1,i2,i3)+ty(i1,i2,i3)*szt8(i1,i2,i3)
        szz83(i1,i2,i3)=rz(i1,i2,i3)*szr8(i1,i2,i3)+sz(i1,i2,i3)*szs8(i1,i2,i3)+tz(i1,i2,i3)*szt8(i1,i2,i3)
        txx82(i1,i2,i3)= rx(i1,i2,i3)*txr8(i1,i2,i3)+sx(i1,i2,i3)*txs8(i1,i2,i3)
        txy82(i1,i2,i3)= ry(i1,i2,i3)*txr8(i1,i2,i3)+sy(i1,i2,i3)*txs8(i1,i2,i3)
        txx83(i1,i2,i3)=rx(i1,i2,i3)*txr8(i1,i2,i3)+sx(i1,i2,i3)*txs8(i1,i2,i3)+tx(i1,i2,i3)*txt8(i1,i2,i3)
        txy83(i1,i2,i3)=ry(i1,i2,i3)*txr8(i1,i2,i3)+sy(i1,i2,i3)*txs8(i1,i2,i3)+ty(i1,i2,i3)*txt8(i1,i2,i3)
        txz83(i1,i2,i3)=rz(i1,i2,i3)*txr8(i1,i2,i3)+sz(i1,i2,i3)*txs8(i1,i2,i3)+tz(i1,i2,i3)*txt8(i1,i2,i3)
        tyx82(i1,i2,i3)= rx(i1,i2,i3)*tyr8(i1,i2,i3)+sx(i1,i2,i3)*tys8(i1,i2,i3)
        tyy82(i1,i2,i3)= ry(i1,i2,i3)*tyr8(i1,i2,i3)+sy(i1,i2,i3)*tys8(i1,i2,i3)
        tyx83(i1,i2,i3)=rx(i1,i2,i3)*tyr8(i1,i2,i3)+sx(i1,i2,i3)*tys8(i1,i2,i3)+tx(i1,i2,i3)*tyt8(i1,i2,i3)
        tyy83(i1,i2,i3)=ry(i1,i2,i3)*tyr8(i1,i2,i3)+sy(i1,i2,i3)*tys8(i1,i2,i3)+ty(i1,i2,i3)*tyt8(i1,i2,i3)
        tyz83(i1,i2,i3)=rz(i1,i2,i3)*tyr8(i1,i2,i3)+sz(i1,i2,i3)*tys8(i1,i2,i3)+tz(i1,i2,i3)*tyt8(i1,i2,i3)
        tzx82(i1,i2,i3)= rx(i1,i2,i3)*tzr8(i1,i2,i3)+sx(i1,i2,i3)*tzs8(i1,i2,i3)
        tzy82(i1,i2,i3)= ry(i1,i2,i3)*tzr8(i1,i2,i3)+sy(i1,i2,i3)*tzs8(i1,i2,i3)
        tzx83(i1,i2,i3)=rx(i1,i2,i3)*tzr8(i1,i2,i3)+sx(i1,i2,i3)*tzs8(i1,i2,i3)+tx(i1,i2,i3)*tzt8(i1,i2,i3)
        tzy83(i1,i2,i3)=ry(i1,i2,i3)*tzr8(i1,i2,i3)+sy(i1,i2,i3)*tzs8(i1,i2,i3)+ty(i1,i2,i3)*tzt8(i1,i2,i3)
        tzz83(i1,i2,i3)=rz(i1,i2,i3)*tzr8(i1,i2,i3)+sz(i1,i2,i3)*tzs8(i1,i2,i3)+tz(i1,i2,i3)*tzt8(i1,i2,i3)
        uxx81(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr8(i1,i2,i3,kd)+(rxx82(i1,i2,i3))*ur8(i1,i2,i3,kd)
        uyy81(i1,i2,i3,kd)=0
        uxy81(i1,i2,i3,kd)=0
        uxz81(i1,i2,i3,kd)=0
        uyz81(i1,i2,i3,kd)=0
        uzz81(i1,i2,i3,kd)=0
        ulaplacian81(i1,i2,i3,kd)=uxx81(i1,i2,i3,kd)
        uxx82(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*urr8(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*urs8(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*uss8(i1,i2,i3,kd)+(rxx82(i1,i2,i3))*ur8(i1,i2,i3,kd)+(sxx82(i1,i2,i3))*us8(i1,i2,i3,kd)
        uyy82(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*urr8(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*urs8(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*uss8(i1,i2,i3,kd)+(ryy82(i1,i2,i3))*ur8(i1,i2,i3,kd)+(syy82(i1,i2,i3))*us8(i1,i2,i3,kd)
        uxy82(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr8(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs8(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss8(i1,i2,i3,kd)+rxy82(i1,i2,i3)*ur8(i1,i2,i3,kd)+sxy82(i1,i2,i3)*us8(i1,i2,i3,kd)
        uxz82(i1,i2,i3,kd)=0
        uyz82(i1,i2,i3,kd)=0
        uzz82(i1,i2,i3,kd)=0
        ulaplacian82(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*urr8(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*urs8(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*uss8(i1,i2,i3,kd)+(rxx82(i1,i2,i3)+ryy82(i1,i2,i3))*ur8(i1,i2,i3,kd)+(sxx82(i1,i2,i3)+syy82(i1,i2,i3))*us8(i1,i2,i3,kd)
        uxx83(i1,i2,i3,kd)=rx(i1,i2,i3)**2*urr8(i1,i2,i3,kd)+sx(i1,i2,i3)**2*uss8(i1,i2,i3,kd)+tx(i1,i2,i3)**2*utt8(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*urs8(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*urt8(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*ust8(i1,i2,i3,kd)+rxx83(i1,i2,i3)*ur8(i1,i2,i3,kd)+sxx83(i1,i2,i3)*us8(i1,i2,i3,kd)+txx83(i1,i2,i3)*ut8(i1,i2,i3,kd)
        uyy83(i1,i2,i3,kd)=ry(i1,i2,i3)**2*urr8(i1,i2,i3,kd)+sy(i1,i2,i3)**2*uss8(i1,i2,i3,kd)+ty(i1,i2,i3)**2*utt8(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*urs8(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*urt8(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*ust8(i1,i2,i3,kd)+ryy83(i1,i2,i3)*ur8(i1,i2,i3,kd)+syy83(i1,i2,i3)*us8(i1,i2,i3,kd)+tyy83(i1,i2,i3)*ut8(i1,i2,i3,kd)
        uzz83(i1,i2,i3,kd)=rz(i1,i2,i3)**2*urr8(i1,i2,i3,kd)+sz(i1,i2,i3)**2*uss8(i1,i2,i3,kd)+tz(i1,i2,i3)**2*utt8(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*urs8(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*urt8(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*ust8(i1,i2,i3,kd)+rzz83(i1,i2,i3)*ur8(i1,i2,i3,kd)+szz83(i1,i2,i3)*us8(i1,i2,i3,kd)+tzz83(i1,i2,i3)*ut8(i1,i2,i3,kd)
        uxy83(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*urr8(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*uss8(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*utt8(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*urs8(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*urt8(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*ust8(i1,i2,i3,kd)+rxy83(i1,i2,i3)*ur8(i1,i2,i3,kd)+sxy83(i1,i2,i3)*us8(i1,i2,i3,kd)+txy83(i1,i2,i3)*ut8(i1,i2,i3,kd)
        uxz83(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*urr8(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*uss8(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*utt8(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*urs8(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*urt8(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*ust8(i1,i2,i3,kd)+rxz83(i1,i2,i3)*ur8(i1,i2,i3,kd)+sxz83(i1,i2,i3)*us8(i1,i2,i3,kd)+txz83(i1,i2,i3)*ut8(i1,i2,i3,kd)
        uyz83(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*urr8(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*uss8(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*utt8(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*urs8(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*urt8(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*ust8(i1,i2,i3,kd)+ryz83(i1,i2,i3)*ur8(i1,i2,i3,kd)+syz83(i1,i2,i3)*us8(i1,i2,i3,kd)+tyz83(i1,i2,i3)*ut8(i1,i2,i3,kd)
        ulaplacian83(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*urr8(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*uss8(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*utt8(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*urs8(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*urt8(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*ust8(i1,i2,i3,kd)+(rxx83(i1,i2,i3)+ryy83(i1,i2,i3)+rzz83(i1,i2,i3))*ur8(i1,i2,i3,kd)+(sxx83(i1,i2,i3)+syy83(i1,i2,i3)+szz83(i1,i2,i3))*us8(i1,i2,i3,kd)+(txx83(i1,i2,i3)+tyy83(i1,i2,i3)+tzz83(i1,i2,i3))*ut8(i1,i2,i3,kd)
    !============================================================================================
    ! Define derivatives for a rectangular grid
    !
    !============================================================================================
        h18(kd) = 1./(840.*dx(kd))
        h28(kd) = 1./(5040.*dx(kd)**2)
        ux83r(i1,i2,i3,kd)=(672.*(u(i1+1,i2,i3,kd)-u(i1-1,i2,i3,kd))-168.*(u(i1+2,i2,i3,kd)-u(i1-2,i2,i3,kd))+32.*(u(i1+3,i2,i3,kd)-u(i1-3,i2,i3,kd))-3.*(u(i1+4,i2,i3,kd)-u(i1-4,i2,i3,kd)))*h18(0)
        uy83r(i1,i2,i3,kd)=(672.*(u(i1,i2+1,i3,kd)-u(i1,i2-1,i3,kd))-168.*(u(i1,i2+2,i3,kd)-u(i1,i2-2,i3,kd))+32.*(u(i1,i2+3,i3,kd)-u(i1,i2-3,i3,kd))-3.*(u(i1,i2+4,i3,kd)-u(i1,i2-4,i3,kd)))*h18(1)
        uz83r(i1,i2,i3,kd)=(672.*(u(i1,i2,i3+1,kd)-u(i1,i2,i3-1,kd))-168.*(u(i1,i2,i3+2,kd)-u(i1,i2,i3-2,kd))+32.*(u(i1,i2,i3+3,kd)-u(i1,i2,i3-3,kd))-3.*(u(i1,i2,i3+4,kd)-u(i1,i2,i3-4,kd)))*h18(2)
        uxx83r(i1,i2,i3,kd)=(-14350.*u(i1,i2,i3,kd)+8064.*(u(i1+1,i2,i3,kd)+u(i1-1,i2,i3,kd))-1008.*(u(i1+2,i2,i3,kd)+u(i1-2,i2,i3,kd))+128.*(u(i1+3,i2,i3,kd)+u(i1-3,i2,i3,kd))-9.*(u(i1+4,i2,i3,kd)+u(i1-4,i2,i3,kd)) )*h28(0)
        uyy83r(i1,i2,i3,kd)=(-14350.*u(i1,i2,i3,kd)+8064.*(u(i1,i2+1,i3,kd)+u(i1,i2-1,i3,kd))-1008.*(u(i1,i2+2,i3,kd)+u(i1,i2-2,i3,kd))+128.*(u(i1,i2+3,i3,kd)+u(i1,i2-3,i3,kd))-9.*(u(i1,i2+4,i3,kd)+u(i1,i2-4,i3,kd)) )*h28(1)
        uzz83r(i1,i2,i3,kd)=(-14350.*u(i1,i2,i3,kd)+8064.*(u(i1,i2,i3+1,kd)+u(i1,i2,i3-1,kd))-1008.*(u(i1,i2,i3+2,kd)+u(i1,i2,i3-2,kd))+128.*(u(i1,i2,i3+3,kd)+u(i1,i2,i3-3,kd))-9.*(u(i1,i2,i3+4,kd)+u(i1,i2,i3-4,kd)) )*h28(2)
        uxy83r(i1,i2,i3,kd)=(672.*(ux83r(i1,i2+1,i3,kd)-ux83r(i1,i2-1,i3,kd))-168.*(ux83r(i1,i2+2,i3,kd)-ux83r(i1,i2-2,i3,kd))+32.*(ux83r(i1,i2+3,i3,kd)-ux83r(i1,i2-3,i3,kd))-3.*(ux83r(i1,i2+4,i3,kd)-ux83r(i1,i2-4,i3,kd)))*h18(1)
        uxz83r(i1,i2,i3,kd)=(672.*(ux83r(i1,i2,i3+1,kd)-ux83r(i1,i2,i3-1,kd))-168.*(ux83r(i1,i2,i3+2,kd)-ux83r(i1,i2,i3-2,kd))+32.*(ux83r(i1,i2,i3+3,kd)-ux83r(i1,i2,i3-3,kd))-3.*(ux83r(i1,i2,i3+4,kd)-ux83r(i1,i2,i3-4,kd)))*h18(2)
        uyz83r(i1,i2,i3,kd)=(672.*(uy83r(i1,i2,i3+1,kd)-uy83r(i1,i2,i3-1,kd))-168.*(uy83r(i1,i2,i3+2,kd)-uy83r(i1,i2,i3-2,kd))+32.*(uy83r(i1,i2,i3+3,kd)-uy83r(i1,i2,i3-3,kd))-3.*(uy83r(i1,i2,i3+4,kd)-uy83r(i1,i2,i3-4,kd)))*h18(2)
        ux81r(i1,i2,i3,kd)= ux83r(i1,i2,i3,kd)
        uy81r(i1,i2,i3,kd)= uy83r(i1,i2,i3,kd)
        uz81r(i1,i2,i3,kd)= uz83r(i1,i2,i3,kd)
        uxx81r(i1,i2,i3,kd)= uxx83r(i1,i2,i3,kd)
        uyy81r(i1,i2,i3,kd)= uyy83r(i1,i2,i3,kd)
        uzz81r(i1,i2,i3,kd)= uzz83r(i1,i2,i3,kd)
        uxy81r(i1,i2,i3,kd)= uxy83r(i1,i2,i3,kd)
        uxz81r(i1,i2,i3,kd)= uxz83r(i1,i2,i3,kd)
        uyz81r(i1,i2,i3,kd)= uyz83r(i1,i2,i3,kd)
        ulaplacian81r(i1,i2,i3,kd)=uxx83r(i1,i2,i3,kd)
        ux82r(i1,i2,i3,kd)= ux83r(i1,i2,i3,kd)
        uy82r(i1,i2,i3,kd)= uy83r(i1,i2,i3,kd)
        uz82r(i1,i2,i3,kd)= uz83r(i1,i2,i3,kd)
        uxx82r(i1,i2,i3,kd)= uxx83r(i1,i2,i3,kd)
        uyy82r(i1,i2,i3,kd)= uyy83r(i1,i2,i3,kd)
        uzz82r(i1,i2,i3,kd)= uzz83r(i1,i2,i3,kd)
        uxy82r(i1,i2,i3,kd)= uxy83r(i1,i2,i3,kd)
        uxz82r(i1,i2,i3,kd)= uxz83r(i1,i2,i3,kd)
        uyz82r(i1,i2,i3,kd)= uyz83r(i1,i2,i3,kd)
        ulaplacian82r(i1,i2,i3,kd)=uxx83r(i1,i2,i3,kd)+uyy83r(i1,i2,i3,kd)
        ulaplacian83r(i1,i2,i3,kd)=uxx83r(i1,i2,i3,kd)+uyy83r(i1,i2,i3,kd)+uzz83r(i1,i2,i3,kd)
    ! Difference approximations um (old time)
        umr2(i1,i2,i3,kd)=(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))*d12(0)
        ums2(i1,i2,i3,kd)=(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))*d12(1)
        umt2(i1,i2,i3,kd)=(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))*d12(2)
        umrr2(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)) )*d22(0)
        umss2(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd)) )*d22(1)
        umrs2(i1,i2,i3,kd)=(umr2(i1,i2+1,i3,kd)-umr2(i1,i2-1,i3,kd))*d12(1)
        umtt2(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd)) )*d22(2)
        umrt2(i1,i2,i3,kd)=(umr2(i1,i2,i3+1,kd)-umr2(i1,i2,i3-1,kd))*d12(2)
        umst2(i1,i2,i3,kd)=(ums2(i1,i2,i3+1,kd)-ums2(i1,i2,i3-1,kd))*d12(2)
        umrrr2(i1,i2,i3,kd)=(-2.*(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))+(um(i1+2,i2,i3,kd)-um(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
        umsss2(i1,i2,i3,kd)=(-2.*(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))+(um(i1,i2+2,i3,kd)-um(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
        umttt2(i1,i2,i3,kd)=(-2.*(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))+(um(i1,i2,i3+2,kd)-um(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
        umx21(i1,i2,i3,kd)= rx(i1,i2,i3)*umr2(i1,i2,i3,kd)
        umy21(i1,i2,i3,kd)=0
        umz21(i1,i2,i3,kd)=0
        umx22(i1,i2,i3,kd)= rx(i1,i2,i3)*umr2(i1,i2,i3,kd)+sx(i1,i2,i3)*ums2(i1,i2,i3,kd)
        umy22(i1,i2,i3,kd)= ry(i1,i2,i3)*umr2(i1,i2,i3,kd)+sy(i1,i2,i3)*ums2(i1,i2,i3,kd)
        umz22(i1,i2,i3,kd)=0
        umx23(i1,i2,i3,kd)=rx(i1,i2,i3)*umr2(i1,i2,i3,kd)+sx(i1,i2,i3)*ums2(i1,i2,i3,kd)+tx(i1,i2,i3)*umt2(i1,i2,i3,kd)
        umy23(i1,i2,i3,kd)=ry(i1,i2,i3)*umr2(i1,i2,i3,kd)+sy(i1,i2,i3)*ums2(i1,i2,i3,kd)+ty(i1,i2,i3)*umt2(i1,i2,i3,kd)
        umz23(i1,i2,i3,kd)=rz(i1,i2,i3)*umr2(i1,i2,i3,kd)+sz(i1,i2,i3)*ums2(i1,i2,i3,kd)+tz(i1,i2,i3)*umt2(i1,i2,i3,kd)
        umxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*umrr2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*umr2(i1,i2,i3,kd)
        umyy21(i1,i2,i3,kd)=0
        umxy21(i1,i2,i3,kd)=0
        umxz21(i1,i2,i3,kd)=0
        umyz21(i1,i2,i3,kd)=0
        umzz21(i1,i2,i3,kd)=0
        umlaplacian21(i1,i2,i3,kd)=umxx21(i1,i2,i3,kd)
        umxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*umrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*umss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*umr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3))*ums2(i1,i2,i3,kd)
        umyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*umrr2(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*umss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*umr2(i1,i2,i3,kd)+(syy22(i1,i2,i3))*ums2(i1,i2,i3,kd)
        umxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*umrr2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*umrs2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*umss2(i1,i2,i3,kd)+rxy22(i1,i2,i3)*umr2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*ums2(i1,i2,i3,kd)
        umxz22(i1,i2,i3,kd)=0
        umyz22(i1,i2,i3,kd)=0
        umzz22(i1,i2,i3,kd)=0
        umlaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*umrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*umss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*umr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*ums2(i1,i2,i3,kd)
        umxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*umrr2(i1,i2,i3,kd)+sx(i1,i2,i3)**2*umss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*umtt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*umrs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*umrt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*umst2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*umr2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*ums2(i1,i2,i3,kd)+txx23(i1,i2,i3)*umt2(i1,i2,i3,kd)
        umyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*umrr2(i1,i2,i3,kd)+sy(i1,i2,i3)**2*umss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*umtt2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*umrs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*umrt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*umst2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*umr2(i1,i2,i3,kd)+syy23(i1,i2,i3)*ums2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*umt2(i1,i2,i3,kd)
        umzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*umrr2(i1,i2,i3,kd)+sz(i1,i2,i3)**2*umss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*umtt2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*umrs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*umrt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*umst2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*umr2(i1,i2,i3,kd)+szz23(i1,i2,i3)*ums2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*umt2(i1,i2,i3,kd)
        umxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*umrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*umss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*umtt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*umrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*umst2(i1,i2,i3,kd)+rxy23(i1,i2,i3)*umr2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*ums2(i1,i2,i3,kd)+txy23(i1,i2,i3)*umt2(i1,i2,i3,kd)
        umxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*umrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*umss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*umtt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*umrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*umst2(i1,i2,i3,kd)+rxz23(i1,i2,i3)*umr2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*ums2(i1,i2,i3,kd)+txz23(i1,i2,i3)*umt2(i1,i2,i3,kd)
        umyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*umrr2(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*umss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*umtt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*umrs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*umrt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*umst2(i1,i2,i3,kd)+ryz23(i1,i2,i3)*umr2(i1,i2,i3,kd)+syz23(i1,i2,i3)*ums2(i1,i2,i3,kd)+tyz23(i1,i2,i3)*umt2(i1,i2,i3,kd)
        umlaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*umrr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*umss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*umtt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*umrs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*umrt2(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*umst2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*umr2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*ums2(i1,i2,i3,kd)+(txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*umt2(i1,i2,i3,kd)
    !============================================================================================
    ! Define derivatives for a rectangular grid
    !
    !============================================================================================
        umx23r(i1,i2,i3,kd)=(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))*h12(0)
        umy23r(i1,i2,i3,kd)=(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))*h12(1)
        umz23r(i1,i2,i3,kd)=(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))*h12(2)
        umxx23r(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)) )*h22(0)
        umyy23r(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd)) )*h22(1)
        umxy23r(i1,i2,i3,kd)=(umx23r(i1,i2+1,i3,kd)-umx23r(i1,i2-1,i3,kd))*h12(1)
        umzz23r(i1,i2,i3,kd)=(-2.*um(i1,i2,i3,kd)+(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd)) )*h22(2)
        umxz23r(i1,i2,i3,kd)=(umx23r(i1,i2,i3+1,kd)-umx23r(i1,i2,i3-1,kd))*h12(2)
        umyz23r(i1,i2,i3,kd)=(umy23r(i1,i2,i3+1,kd)-umy23r(i1,i2,i3-1,kd))*h12(2)
        umx21r(i1,i2,i3,kd)= umx23r(i1,i2,i3,kd)
        umy21r(i1,i2,i3,kd)= umy23r(i1,i2,i3,kd)
        umz21r(i1,i2,i3,kd)= umz23r(i1,i2,i3,kd)
        umxx21r(i1,i2,i3,kd)= umxx23r(i1,i2,i3,kd)
        umyy21r(i1,i2,i3,kd)= umyy23r(i1,i2,i3,kd)
        umzz21r(i1,i2,i3,kd)= umzz23r(i1,i2,i3,kd)
        umxy21r(i1,i2,i3,kd)= umxy23r(i1,i2,i3,kd)
        umxz21r(i1,i2,i3,kd)= umxz23r(i1,i2,i3,kd)
        umyz21r(i1,i2,i3,kd)= umyz23r(i1,i2,i3,kd)
        umlaplacian21r(i1,i2,i3,kd)=umxx23r(i1,i2,i3,kd)
        umx22r(i1,i2,i3,kd)= umx23r(i1,i2,i3,kd)
        umy22r(i1,i2,i3,kd)= umy23r(i1,i2,i3,kd)
        umz22r(i1,i2,i3,kd)= umz23r(i1,i2,i3,kd)
        umxx22r(i1,i2,i3,kd)= umxx23r(i1,i2,i3,kd)
        umyy22r(i1,i2,i3,kd)= umyy23r(i1,i2,i3,kd)
        umzz22r(i1,i2,i3,kd)= umzz23r(i1,i2,i3,kd)
        umxy22r(i1,i2,i3,kd)= umxy23r(i1,i2,i3,kd)
        umxz22r(i1,i2,i3,kd)= umxz23r(i1,i2,i3,kd)
        umyz22r(i1,i2,i3,kd)= umyz23r(i1,i2,i3,kd)
        umlaplacian22r(i1,i2,i3,kd)=umxx23r(i1,i2,i3,kd)+umyy23r(i1,i2,i3,kd)
        umlaplacian23r(i1,i2,i3,kd)=umxx23r(i1,i2,i3,kd)+umyy23r(i1,i2,i3,kd)+umzz23r(i1,i2,i3,kd)
        umxxx22r(i1,i2,i3,kd)=(-2.*(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))+(um(i1+2,i2,i3,kd)-um(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
        umyyy22r(i1,i2,i3,kd)=(-2.*(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))+(um(i1,i2+2,i3,kd)-um(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
        umxxy22r(i1,i2,i3,kd)=( umxx22r(i1,i2+1,i3,kd)-umxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
        umxyy22r(i1,i2,i3,kd)=( umyy22r(i1+1,i2,i3,kd)-umyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
        umxxxx22r(i1,i2,i3,kd)=(6.*um(i1,i2,i3,kd)-4.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))+(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )/(dx(0)**4)
        umyyyy22r(i1,i2,i3,kd)=(6.*um(i1,i2,i3,kd)-4.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))+(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )/(dx(1)**4)
        umxxyy22r(i1,i2,i3,kd)=( 4.*um(i1,i2,i3,kd)     -2.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)+um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))   +   (um(i1+1,i2+1,i3,kd)+um(i1-1,i2+1,i3,kd)+um(i1+1,i2-1,i3,kd)+um(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
    ! 2D laplacian squared = um.xxxx + 2 um.xxyy + um.yyyy
        umLapSq22r(i1,i2,i3,kd)= ( 6.*um(i1,i2,i3,kd)   - 4.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))    +(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*um(i1,i2,i3,kd)    -4.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))    +(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 8.*um(i1,i2,i3,kd)     -4.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)+um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))   +2.*(um(i1+1,i2+1,i3,kd)+um(i1-1,i2+1,i3,kd)+um(i1+1,i2-1,i3,kd)+um(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
        umxxx23r(i1,i2,i3,kd)=(-2.*(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))+(um(i1+2,i2,i3,kd)-um(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
        umyyy23r(i1,i2,i3,kd)=(-2.*(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))+(um(i1,i2+2,i3,kd)-um(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
        umzzz23r(i1,i2,i3,kd)=(-2.*(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))+(um(i1,i2,i3+2,kd)-um(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
        umxxy23r(i1,i2,i3,kd)=( umxx22r(i1,i2+1,i3,kd)-umxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
        umxyy23r(i1,i2,i3,kd)=( umyy22r(i1+1,i2,i3,kd)-umyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
        umxxz23r(i1,i2,i3,kd)=( umxx22r(i1,i2,i3+1,kd)-umxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
        umyyz23r(i1,i2,i3,kd)=( umyy22r(i1,i2,i3+1,kd)-umyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
        umxzz23r(i1,i2,i3,kd)=( umzz22r(i1+1,i2,i3,kd)-umzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
        umyzz23r(i1,i2,i3,kd)=( umzz22r(i1,i2+1,i3,kd)-umzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
        umxxxx23r(i1,i2,i3,kd)=(6.*um(i1,i2,i3,kd)-4.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))+(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )/(dx(0)**4)
        umyyyy23r(i1,i2,i3,kd)=(6.*um(i1,i2,i3,kd)-4.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))+(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )/(dx(1)**4)
        umzzzz23r(i1,i2,i3,kd)=(6.*um(i1,i2,i3,kd)-4.*(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd))+(um(i1,i2,i3+2,kd)+um(i1,i2,i3-2,kd)) )/(dx(2)**4)
        umxxyy23r(i1,i2,i3,kd)=( 4.*um(i1,i2,i3,kd)     -2.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)+um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))   +   (um(i1+1,i2+1,i3,kd)+um(i1-1,i2+1,i3,kd)+um(i1+1,i2-1,i3,kd)+um(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
        umxxzz23r(i1,i2,i3,kd)=( 4.*um(i1,i2,i3,kd)     -2.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd)+um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd))   +   (um(i1+1,i2,i3+1,kd)+um(i1-1,i2,i3+1,kd)+um(i1+1,i2,i3-1,kd)+um(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
        umyyzz23r(i1,i2,i3,kd)=( 4.*um(i1,i2,i3,kd)     -2.*(um(i1,i2+1,i3,kd)  +um(i1,i2-1,i3,kd)+  um(i1,i2  ,i3+1,kd)+um(i1,i2  ,i3-1,kd))   +   (um(i1,i2+1,i3+1,kd)+um(i1,i2-1,i3+1,kd)+um(i1,i2+1,i3-1,kd)+um(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
    ! 3D laplacian squared = um.xxxx + um.yyyy + um.zzzz + 2 (um.xxyy + um.xxzz + um.yyzz )
        umLapSq23r(i1,i2,i3,kd)= ( 6.*um(i1,i2,i3,kd)   - 4.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))    +(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*um(i1,i2,i3,kd)    -4.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))    +(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 6.*um(i1,i2,i3,kd)    -4.*(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd))    +(um(i1,i2,i3+2,kd)+um(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 8.*um(i1,i2,i3,kd)     -4.*(um(i1+1,i2,i3,kd)  +um(i1-1,i2,i3,kd)  +um(i1  ,i2+1,i3,kd)+um(i1  ,i2-1,i3,kd))   +2.*(um(i1+1,i2+1,i3,kd)+um(i1-1,i2+1,i3,kd)+um(i1+1,i2-1,i3,kd)+um(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*um(i1,i2,i3,kd)     -4.*(um(i1+1,i2,i3,kd)  +um(i1-1,i2,i3,kd)  +um(i1  ,i2,i3+1,kd)+um(i1  ,i2,i3-1,kd))   +2.*(um(i1+1,i2,i3+1,kd)+um(i1-1,i2,i3+1,kd)+um(i1+1,i2,i3-1,kd)+um(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*um(i1,i2,i3,kd)     -4.*(um(i1,i2+1,i3,kd)  +um(i1,i2-1,i3,kd)  +um(i1,i2  ,i3+1,kd)+um(i1,i2  ,i3-1,kd))   +2.*(um(i1,i2+1,i3+1,kd)+um(i1,i2-1,i3+1,kd)+um(i1,i2+1,i3-1,kd)+um(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
        umr4(i1,i2,i3,kd)=(8.*(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))-(um(i1+2,i2,i3,kd)-um(i1-2,i2,i3,kd)))*d14(0)
        ums4(i1,i2,i3,kd)=(8.*(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))-(um(i1,i2+2,i3,kd)-um(i1,i2-2,i3,kd)))*d14(1)
        umt4(i1,i2,i3,kd)=(8.*(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))-(um(i1,i2,i3+2,kd)-um(i1,i2,i3-2,kd)))*d14(2)
        umrr4(i1,i2,i3,kd)=(-30.*um(i1,i2,i3,kd)+16.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))-(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )*d24(0)
        umss4(i1,i2,i3,kd)=(-30.*um(i1,i2,i3,kd)+16.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))-(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )*d24(1)
        umtt4(i1,i2,i3,kd)=(-30.*um(i1,i2,i3,kd)+16.*(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd))-(um(i1,i2,i3+2,kd)+um(i1,i2,i3-2,kd)) )*d24(2)
        umrs4(i1,i2,i3,kd)=(8.*(umr4(i1,i2+1,i3,kd)-umr4(i1,i2-1,i3,kd))-(umr4(i1,i2+2,i3,kd)-umr4(i1,i2-2,i3,kd)))*d14(1)
        umrt4(i1,i2,i3,kd)=(8.*(umr4(i1,i2,i3+1,kd)-umr4(i1,i2,i3-1,kd))-(umr4(i1,i2,i3+2,kd)-umr4(i1,i2,i3-2,kd)))*d14(2)
        umst4(i1,i2,i3,kd)=(8.*(ums4(i1,i2,i3+1,kd)-ums4(i1,i2,i3-1,kd))-(ums4(i1,i2,i3+2,kd)-ums4(i1,i2,i3-2,kd)))*d14(2)
        umx41(i1,i2,i3,kd)= rx(i1,i2,i3)*umr4(i1,i2,i3,kd)
        umy41(i1,i2,i3,kd)=0
        umz41(i1,i2,i3,kd)=0
        umx42(i1,i2,i3,kd)= rx(i1,i2,i3)*umr4(i1,i2,i3,kd)+sx(i1,i2,i3)*ums4(i1,i2,i3,kd)
        umy42(i1,i2,i3,kd)= ry(i1,i2,i3)*umr4(i1,i2,i3,kd)+sy(i1,i2,i3)*ums4(i1,i2,i3,kd)
        umz42(i1,i2,i3,kd)=0
        umx43(i1,i2,i3,kd)=rx(i1,i2,i3)*umr4(i1,i2,i3,kd)+sx(i1,i2,i3)*ums4(i1,i2,i3,kd)+tx(i1,i2,i3)*umt4(i1,i2,i3,kd)
        umy43(i1,i2,i3,kd)=ry(i1,i2,i3)*umr4(i1,i2,i3,kd)+sy(i1,i2,i3)*ums4(i1,i2,i3,kd)+ty(i1,i2,i3)*umt4(i1,i2,i3,kd)
        umz43(i1,i2,i3,kd)=rz(i1,i2,i3)*umr4(i1,i2,i3,kd)+sz(i1,i2,i3)*ums4(i1,i2,i3,kd)+tz(i1,i2,i3)*umt4(i1,i2,i3,kd)
        umxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*umrr4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*umr4(i1,i2,i3,kd)
        umyy41(i1,i2,i3,kd)=0
        umxy41(i1,i2,i3,kd)=0
        umxz41(i1,i2,i3,kd)=0
        umyz41(i1,i2,i3,kd)=0
        umzz41(i1,i2,i3,kd)=0
        umlaplacian41(i1,i2,i3,kd)=umxx41(i1,i2,i3,kd)
        umxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*umrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*umss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*umr4(i1,i2,i3,kd)+(sxx42(i1,i2,i3))*ums4(i1,i2,i3,kd)
        umyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*umrr4(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*umss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*umr4(i1,i2,i3,kd)+(syy42(i1,i2,i3))*ums4(i1,i2,i3,kd)
        umxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*umrr4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*umrs4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*umss4(i1,i2,i3,kd)+rxy42(i1,i2,i3)*umr4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*ums4(i1,i2,i3,kd)
        umxz42(i1,i2,i3,kd)=0
        umyz42(i1,i2,i3,kd)=0
        umzz42(i1,i2,i3,kd)=0
        umlaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*umrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*umss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*umr4(i1,i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*ums4(i1,i2,i3,kd)
        umxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*umrr4(i1,i2,i3,kd)+sx(i1,i2,i3)**2*umss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*umtt4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*umrs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*umrt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*umst4(i1,i2,i3,kd)+rxx43(i1,i2,i3)*umr4(i1,i2,i3,kd)+sxx43(i1,i2,i3)*ums4(i1,i2,i3,kd)+txx43(i1,i2,i3)*umt4(i1,i2,i3,kd)
        umyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*umrr4(i1,i2,i3,kd)+sy(i1,i2,i3)**2*umss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*umtt4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*umrs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*umrt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*umst4(i1,i2,i3,kd)+ryy43(i1,i2,i3)*umr4(i1,i2,i3,kd)+syy43(i1,i2,i3)*ums4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*umt4(i1,i2,i3,kd)
        umzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*umrr4(i1,i2,i3,kd)+sz(i1,i2,i3)**2*umss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*umtt4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*umrs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*umrt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*umst4(i1,i2,i3,kd)+rzz43(i1,i2,i3)*umr4(i1,i2,i3,kd)+szz43(i1,i2,i3)*ums4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*umt4(i1,i2,i3,kd)
        umxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*umrr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*umss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*umtt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*umrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*umst4(i1,i2,i3,kd)+rxy43(i1,i2,i3)*umr4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*ums4(i1,i2,i3,kd)+txy43(i1,i2,i3)*umt4(i1,i2,i3,kd)
        umxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*umrr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*umss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*umtt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*umrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*umst4(i1,i2,i3,kd)+rxz43(i1,i2,i3)*umr4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*ums4(i1,i2,i3,kd)+txz43(i1,i2,i3)*umt4(i1,i2,i3,kd)
        umyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*umrr4(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*umss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*umtt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*umrs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*umrt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*umst4(i1,i2,i3,kd)+ryz43(i1,i2,i3)*umr4(i1,i2,i3,kd)+syz43(i1,i2,i3)*ums4(i1,i2,i3,kd)+tyz43(i1,i2,i3)*umt4(i1,i2,i3,kd)
        umlaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*umrr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*umss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*umtt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*umrs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*umrt4(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*umst4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+ryy43(i1,i2,i3)+rzz43(i1,i2,i3))*umr4(i1,i2,i3,kd)+(sxx43(i1,i2,i3)+syy43(i1,i2,i3)+szz43(i1,i2,i3))*ums4(i1,i2,i3,kd)+(txx43(i1,i2,i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*umt4(i1,i2,i3,kd)
    !============================================================================================
    ! Define derivatives for a rectangular grid
    !
    !============================================================================================
        umx43r(i1,i2,i3,kd)=(8.*(um(i1+1,i2,i3,kd)-um(i1-1,i2,i3,kd))-(um(i1+2,i2,i3,kd)-um(i1-2,i2,i3,kd)))*h41(0)
        umy43r(i1,i2,i3,kd)=(8.*(um(i1,i2+1,i3,kd)-um(i1,i2-1,i3,kd))-(um(i1,i2+2,i3,kd)-um(i1,i2-2,i3,kd)))*h41(1)
        umz43r(i1,i2,i3,kd)=(8.*(um(i1,i2,i3+1,kd)-um(i1,i2,i3-1,kd))-(um(i1,i2,i3+2,kd)-um(i1,i2,i3-2,kd)))*h41(2)
        umxx43r(i1,i2,i3,kd)=( -30.*um(i1,i2,i3,kd)+16.*(um(i1+1,i2,i3,kd)+um(i1-1,i2,i3,kd))-(um(i1+2,i2,i3,kd)+um(i1-2,i2,i3,kd)) )*h42(0) 
        umyy43r(i1,i2,i3,kd)=( -30.*um(i1,i2,i3,kd)+16.*(um(i1,i2+1,i3,kd)+um(i1,i2-1,i3,kd))-(um(i1,i2+2,i3,kd)+um(i1,i2-2,i3,kd)) )*h42(1) 
        umzz43r(i1,i2,i3,kd)=( -30.*um(i1,i2,i3,kd)+16.*(um(i1,i2,i3+1,kd)+um(i1,i2,i3-1,kd))-(um(i1,i2,i3+2,kd)+um(i1,i2,i3-2,kd)) )*h42(2)
        umxy43r(i1,i2,i3,kd)=( (um(i1+2,i2+2,i3,kd)-um(i1-2,i2+2,i3,kd)- um(i1+2,i2-2,i3,kd)+um(i1-2,i2-2,i3,kd)) +8.*(um(i1-1,i2+2,i3,kd)-um(i1-1,i2-2,i3,kd)-um(i1+1,i2+2,i3,kd)+um(i1+1,i2-2,i3,kd) +um(i1+2,i2-1,i3,kd)-um(i1-2,i2-1,i3,kd)-um(i1+2,i2+1,i3,kd)+um(i1-2,i2+1,i3,kd))+64.*(um(i1+1,i2+1,i3,kd)-um(i1-1,i2+1,i3,kd)- um(i1+1,i2-1,i3,kd)+um(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
        umxz43r(i1,i2,i3,kd)=( (um(i1+2,i2,i3+2,kd)-um(i1-2,i2,i3+2,kd)-um(i1+2,i2,i3-2,kd)+um(i1-2,i2,i3-2,kd)) +8.*(um(i1-1,i2,i3+2,kd)-um(i1-1,i2,i3-2,kd)-um(i1+1,i2,i3+2,kd)+um(i1+1,i2,i3-2,kd) +um(i1+2,i2,i3-1,kd)-um(i1-2,i2,i3-1,kd)- um(i1+2,i2,i3+1,kd)+um(i1-2,i2,i3+1,kd)) +64.*(um(i1+1,i2,i3+1,kd)-um(i1-1,i2,i3+1,kd)-um(i1+1,i2,i3-1,kd)+um(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
        umyz43r(i1,i2,i3,kd)=( (um(i1,i2+2,i3+2,kd)-um(i1,i2-2,i3+2,kd)-um(i1,i2+2,i3-2,kd)+um(i1,i2-2,i3-2,kd)) +8.*(um(i1,i2-1,i3+2,kd)-um(i1,i2-1,i3-2,kd)-um(i1,i2+1,i3+2,kd)+um(i1,i2+1,i3-2,kd) +um(i1,i2+2,i3-1,kd)-um(i1,i2-2,i3-1,kd)-um(i1,i2+2,i3+1,kd)+um(i1,i2-2,i3+1,kd)) +64.*(um(i1,i2+1,i3+1,kd)-um(i1,i2-1,i3+1,kd)-um(i1,i2+1,i3-1,kd)+um(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
        umx41r(i1,i2,i3,kd)= umx43r(i1,i2,i3,kd)
        umy41r(i1,i2,i3,kd)= umy43r(i1,i2,i3,kd)
        umz41r(i1,i2,i3,kd)= umz43r(i1,i2,i3,kd)
        umxx41r(i1,i2,i3,kd)= umxx43r(i1,i2,i3,kd)
        umyy41r(i1,i2,i3,kd)= umyy43r(i1,i2,i3,kd)
        umzz41r(i1,i2,i3,kd)= umzz43r(i1,i2,i3,kd)
        umxy41r(i1,i2,i3,kd)= umxy43r(i1,i2,i3,kd)
        umxz41r(i1,i2,i3,kd)= umxz43r(i1,i2,i3,kd)
        umyz41r(i1,i2,i3,kd)= umyz43r(i1,i2,i3,kd)
        umlaplacian41r(i1,i2,i3,kd)=umxx43r(i1,i2,i3,kd)
        umx42r(i1,i2,i3,kd)= umx43r(i1,i2,i3,kd)
        umy42r(i1,i2,i3,kd)= umy43r(i1,i2,i3,kd)
        umz42r(i1,i2,i3,kd)= umz43r(i1,i2,i3,kd)
        umxx42r(i1,i2,i3,kd)= umxx43r(i1,i2,i3,kd)
        umyy42r(i1,i2,i3,kd)= umyy43r(i1,i2,i3,kd)
        umzz42r(i1,i2,i3,kd)= umzz43r(i1,i2,i3,kd)
        umxy42r(i1,i2,i3,kd)= umxy43r(i1,i2,i3,kd)
        umxz42r(i1,i2,i3,kd)= umxz43r(i1,i2,i3,kd)
        umyz42r(i1,i2,i3,kd)= umyz43r(i1,i2,i3,kd)
        umlaplacian42r(i1,i2,i3,kd)=umxx43r(i1,i2,i3,kd)+umyy43r(i1,i2,i3,kd)
        umlaplacian43r(i1,i2,i3,kd)=umxx43r(i1,i2,i3,kd)+umyy43r(i1,i2,i3,kd)+umzz43r(i1,i2,i3,kd)
    ! Define difference approximations for v
        vr2(i1,i2,i3,kd)=(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))*d12(0)
        vs2(i1,i2,i3,kd)=(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))*d12(1)
        vt2(i1,i2,i3,kd)=(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))*d12(2)
        vrr2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)) )*d22(0)
        vss2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd)) )*d22(1)
        vrs2(i1,i2,i3,kd)=(vr2(i1,i2+1,i3,kd)-vr2(i1,i2-1,i3,kd))*d12(1)
        vtt2(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd)) )*d22(2)
        vrt2(i1,i2,i3,kd)=(vr2(i1,i2,i3+1,kd)-vr2(i1,i2,i3-1,kd))*d12(2)
        vst2(i1,i2,i3,kd)=(vs2(i1,i2,i3+1,kd)-vs2(i1,i2,i3-1,kd))*d12(2)
        vrrr2(i1,i2,i3,kd)=(-2.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
        vsss2(i1,i2,i3,kd)=(-2.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
        vttt2(i1,i2,i3,kd)=(-2.*(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))+(v(i1,i2,i3+2,kd)-v(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
        vx21(i1,i2,i3,kd)= rx(i1,i2,i3)*vr2(i1,i2,i3,kd)
        vy21(i1,i2,i3,kd)=0
        vz21(i1,i2,i3,kd)=0
        vx22(i1,i2,i3,kd)= rx(i1,i2,i3)*vr2(i1,i2,i3,kd)+sx(i1,i2,i3)*vs2(i1,i2,i3,kd)
        vy22(i1,i2,i3,kd)= ry(i1,i2,i3)*vr2(i1,i2,i3,kd)+sy(i1,i2,i3)*vs2(i1,i2,i3,kd)
        vz22(i1,i2,i3,kd)=0
        vx23(i1,i2,i3,kd)=rx(i1,i2,i3)*vr2(i1,i2,i3,kd)+sx(i1,i2,i3)*vs2(i1,i2,i3,kd)+tx(i1,i2,i3)*vt2(i1,i2,i3,kd)
        vy23(i1,i2,i3,kd)=ry(i1,i2,i3)*vr2(i1,i2,i3,kd)+sy(i1,i2,i3)*vs2(i1,i2,i3,kd)+ty(i1,i2,i3)*vt2(i1,i2,i3,kd)
        vz23(i1,i2,i3,kd)=rz(i1,i2,i3)*vr2(i1,i2,i3,kd)+sz(i1,i2,i3)*vs2(i1,i2,i3,kd)+tz(i1,i2,i3)*vt2(i1,i2,i3,kd)
        vxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*vr2(i1,i2,i3,kd)
        vyy21(i1,i2,i3,kd)=0
        vxy21(i1,i2,i3,kd)=0
        vxz21(i1,i2,i3,kd)=0
        vyz21(i1,i2,i3,kd)=0
        vzz21(i1,i2,i3,kd)=0
        vlaplacian21(i1,i2,i3,kd)=vxx21(i1,i2,i3,kd)
        vxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*vss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*vr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3))*vs2(i1,i2,i3,kd)
        vyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*vss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*vr2(i1,i2,i3,kd)+(syy22(i1,i2,i3))*vs2(i1,i2,i3,kd)
        vxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*vrr2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*vss2(i1,i2,i3,kd)+rxy22(i1,i2,i3)*vr2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*vs2(i1,i2,i3,kd)
        vxz22(i1,i2,i3,kd)=0
        vyz22(i1,i2,i3,kd)=0
        vzz22(i1,i2,i3,kd)=0
        vlaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*vss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*vr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*vs2(i1,i2,i3,kd)
        vxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*vrr2(i1,i2,i3,kd)+sx(i1,i2,i3)**2*vss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*vtt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*vrs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*vrt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*vst2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*vr2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*vs2(i1,i2,i3,kd)+txx23(i1,i2,i3)*vt2(i1,i2,i3,kd)
        vyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*vrr2(i1,i2,i3,kd)+sy(i1,i2,i3)**2*vss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*vtt2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*vrs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*vrt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*vst2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*vr2(i1,i2,i3,kd)+syy23(i1,i2,i3)*vs2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*vt2(i1,i2,i3,kd)
        vzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*vrr2(i1,i2,i3,kd)+sz(i1,i2,i3)**2*vss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*vtt2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*vrs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*vrt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*vst2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*vr2(i1,i2,i3,kd)+szz23(i1,i2,i3)*vs2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*vt2(i1,i2,i3,kd)
        vxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*vrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*vss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*vtt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*vrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*vst2(i1,i2,i3,kd)+rxy23(i1,i2,i3)*vr2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*vs2(i1,i2,i3,kd)+txy23(i1,i2,i3)*vt2(i1,i2,i3,kd)
        vxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*vrr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*vss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*vtt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*vrt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*vst2(i1,i2,i3,kd)+rxz23(i1,i2,i3)*vr2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*vs2(i1,i2,i3,kd)+txz23(i1,i2,i3)*vt2(i1,i2,i3,kd)
        vyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*vrr2(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*vss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*vtt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*vrs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*vrt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*vst2(i1,i2,i3,kd)+ryz23(i1,i2,i3)*vr2(i1,i2,i3,kd)+syz23(i1,i2,i3)*vs2(i1,i2,i3,kd)+tyz23(i1,i2,i3)*vt2(i1,i2,i3,kd)
        vlaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*vrr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*vss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*vtt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*vrs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*vrt2(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*vst2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*vr2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*vs2(i1,i2,i3,kd)+(txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*vt2(i1,i2,i3,kd)
    !============================================================================================
    ! Define derivatives for a rectangular grid
    !
    !============================================================================================
        vx23r(i1,i2,i3,kd)=(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))*h12(0)
        vy23r(i1,i2,i3,kd)=(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))*h12(1)
        vz23r(i1,i2,i3,kd)=(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))*h12(2)
        vxx23r(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)) )*h22(0)
        vyy23r(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd)) )*h22(1)
        vxy23r(i1,i2,i3,kd)=(vx23r(i1,i2+1,i3,kd)-vx23r(i1,i2-1,i3,kd))*h12(1)
        vzz23r(i1,i2,i3,kd)=(-2.*v(i1,i2,i3,kd)+(v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd)) )*h22(2)
        vxz23r(i1,i2,i3,kd)=(vx23r(i1,i2,i3+1,kd)-vx23r(i1,i2,i3-1,kd))*h12(2)
        vyz23r(i1,i2,i3,kd)=(vy23r(i1,i2,i3+1,kd)-vy23r(i1,i2,i3-1,kd))*h12(2)
        vx21r(i1,i2,i3,kd)= vx23r(i1,i2,i3,kd)
        vy21r(i1,i2,i3,kd)= vy23r(i1,i2,i3,kd)
        vz21r(i1,i2,i3,kd)= vz23r(i1,i2,i3,kd)
        vxx21r(i1,i2,i3,kd)= vxx23r(i1,i2,i3,kd)
        vyy21r(i1,i2,i3,kd)= vyy23r(i1,i2,i3,kd)
        vzz21r(i1,i2,i3,kd)= vzz23r(i1,i2,i3,kd)
        vxy21r(i1,i2,i3,kd)= vxy23r(i1,i2,i3,kd)
        vxz21r(i1,i2,i3,kd)= vxz23r(i1,i2,i3,kd)
        vyz21r(i1,i2,i3,kd)= vyz23r(i1,i2,i3,kd)
        vlaplacian21r(i1,i2,i3,kd)=vxx23r(i1,i2,i3,kd)
        vx22r(i1,i2,i3,kd)= vx23r(i1,i2,i3,kd)
        vy22r(i1,i2,i3,kd)= vy23r(i1,i2,i3,kd)
        vz22r(i1,i2,i3,kd)= vz23r(i1,i2,i3,kd)
        vxx22r(i1,i2,i3,kd)= vxx23r(i1,i2,i3,kd)
        vyy22r(i1,i2,i3,kd)= vyy23r(i1,i2,i3,kd)
        vzz22r(i1,i2,i3,kd)= vzz23r(i1,i2,i3,kd)
        vxy22r(i1,i2,i3,kd)= vxy23r(i1,i2,i3,kd)
        vxz22r(i1,i2,i3,kd)= vxz23r(i1,i2,i3,kd)
        vyz22r(i1,i2,i3,kd)= vyz23r(i1,i2,i3,kd)
        vlaplacian22r(i1,i2,i3,kd)=vxx23r(i1,i2,i3,kd)+vyy23r(i1,i2,i3,kd)
        vlaplacian23r(i1,i2,i3,kd)=vxx23r(i1,i2,i3,kd)+vyy23r(i1,i2,i3,kd)+vzz23r(i1,i2,i3,kd)
        vxxx22r(i1,i2,i3,kd)=(-2.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
        vyyy22r(i1,i2,i3,kd)=(-2.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
        vxxy22r(i1,i2,i3,kd)=( vxx22r(i1,i2+1,i3,kd)-vxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
        vxyy22r(i1,i2,i3,kd)=( vyy22r(i1+1,i2,i3,kd)-vyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
        vxxxx22r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )/(dx(0)**4)
        vyyyy22r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(1)**4)
        vxxyy22r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)+v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))   +   (v(i1+1,i2+1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
    ! 2D laplacian squared = v.xxxx + 2 v.xxyy + v.yyyy
        vLapSq22r(i1,i2,i3,kd)= ( 6.*v(i1,i2,i3,kd)   - 4.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd))    +(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*v(i1,i2,i3,kd)    -4.*(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))    +(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 8.*v(i1,i2,i3,kd)     -4.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)+v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))   +2.*(v(i1+1,i2+1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
        vxxx23r(i1,i2,i3,kd)=(-2.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
        vyyy23r(i1,i2,i3,kd)=(-2.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
        vzzz23r(i1,i2,i3,kd)=(-2.*(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))+(v(i1,i2,i3+2,kd)-v(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
        vxxy23r(i1,i2,i3,kd)=( vxx22r(i1,i2+1,i3,kd)-vxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
        vxyy23r(i1,i2,i3,kd)=( vyy22r(i1+1,i2,i3,kd)-vyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
        vxxz23r(i1,i2,i3,kd)=( vxx22r(i1,i2,i3+1,kd)-vxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
        vyyz23r(i1,i2,i3,kd)=( vyy22r(i1,i2,i3+1,kd)-vyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
        vxzz23r(i1,i2,i3,kd)=( vzz22r(i1+1,i2,i3,kd)-vzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
        vyzz23r(i1,i2,i3,kd)=( vzz22r(i1,i2+1,i3,kd)-vzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
        vxxxx23r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd))+(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )/(dx(0)**4)
        vyyyy23r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))+(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(1)**4)
        vzzzz23r(i1,i2,i3,kd)=(6.*v(i1,i2,i3,kd)-4.*(v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd))+(v(i1,i2,i3+2,kd)+v(i1,i2,i3-2,kd)) )/(dx(2)**4)
        vxxyy23r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)+v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))   +   (v(i1+1,i2+1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
        vxxzz23r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd)+v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd))   +   (v(i1+1,i2,i3+1,kd)+v(i1-1,i2,i3+1,kd)+v(i1+1,i2,i3-1,kd)+v(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
        vyyzz23r(i1,i2,i3,kd)=( 4.*v(i1,i2,i3,kd)     -2.*(v(i1,i2+1,i3,kd)  +v(i1,i2-1,i3,kd)+  v(i1,i2  ,i3+1,kd)+v(i1,i2  ,i3-1,kd))   +   (v(i1,i2+1,i3+1,kd)+v(i1,i2-1,i3+1,kd)+v(i1,i2+1,i3-1,kd)+v(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
    ! 3D laplacian squared = v.xxxx + v.yyyy + v.zzzz + 2 (v.xxyy + v.xxzz + v.yyzz )
        vLapSq23r(i1,i2,i3,kd)= ( 6.*v(i1,i2,i3,kd)   - 4.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd))    +(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*v(i1,i2,i3,kd)    -4.*(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))    +(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 6.*v(i1,i2,i3,kd)    -4.*(v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd))    +(v(i1,i2,i3+2,kd)+v(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 8.*v(i1,i2,i3,kd)     -4.*(v(i1+1,i2,i3,kd)  +v(i1-1,i2,i3,kd)  +v(i1  ,i2+1,i3,kd)+v(i1  ,i2-1,i3,kd))   +2.*(v(i1+1,i2+1,i3,kd)+v(i1-1,i2+1,i3,kd)+v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*v(i1,i2,i3,kd)     -4.*(v(i1+1,i2,i3,kd)  +v(i1-1,i2,i3,kd)  +v(i1  ,i2,i3+1,kd)+v(i1  ,i2,i3-1,kd))   +2.*(v(i1+1,i2,i3+1,kd)+v(i1-1,i2,i3+1,kd)+v(i1+1,i2,i3-1,kd)+v(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*v(i1,i2,i3,kd)     -4.*(v(i1,i2+1,i3,kd)  +v(i1,i2-1,i3,kd)  +v(i1,i2  ,i3+1,kd)+v(i1,i2  ,i3-1,kd))   +2.*(v(i1,i2+1,i3+1,kd)+v(i1,i2-1,i3+1,kd)+v(i1,i2+1,i3-1,kd)+v(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
        vr4(i1,i2,i3,kd)=(8.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))-(v(i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)))*d14(0)
        vs4(i1,i2,i3,kd)=(8.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))-(v(i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)))*d14(1)
        vt4(i1,i2,i3,kd)=(8.*(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))-(v(i1,i2,i3+2,kd)-v(i1,i2,i3-2,kd)))*d14(2)
        vrr4(i1,i2,i3,kd)=(-30.*v(i1,i2,i3,kd)+16.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd))-(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )*d24(0)
        vss4(i1,i2,i3,kd)=(-30.*v(i1,i2,i3,kd)+16.*(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))-(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )*d24(1)
        vtt4(i1,i2,i3,kd)=(-30.*v(i1,i2,i3,kd)+16.*(v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd))-(v(i1,i2,i3+2,kd)+v(i1,i2,i3-2,kd)) )*d24(2)
        vrs4(i1,i2,i3,kd)=(8.*(vr4(i1,i2+1,i3,kd)-vr4(i1,i2-1,i3,kd))-(vr4(i1,i2+2,i3,kd)-vr4(i1,i2-2,i3,kd)))*d14(1)
        vrt4(i1,i2,i3,kd)=(8.*(vr4(i1,i2,i3+1,kd)-vr4(i1,i2,i3-1,kd))-(vr4(i1,i2,i3+2,kd)-vr4(i1,i2,i3-2,kd)))*d14(2)
        vst4(i1,i2,i3,kd)=(8.*(vs4(i1,i2,i3+1,kd)-vs4(i1,i2,i3-1,kd))-(vs4(i1,i2,i3+2,kd)-vs4(i1,i2,i3-2,kd)))*d14(2)
        vx41(i1,i2,i3,kd)= rx(i1,i2,i3)*vr4(i1,i2,i3,kd)
        vy41(i1,i2,i3,kd)=0
        vz41(i1,i2,i3,kd)=0
        vx42(i1,i2,i3,kd)= rx(i1,i2,i3)*vr4(i1,i2,i3,kd)+sx(i1,i2,i3)*vs4(i1,i2,i3,kd)
        vy42(i1,i2,i3,kd)= ry(i1,i2,i3)*vr4(i1,i2,i3,kd)+sy(i1,i2,i3)*vs4(i1,i2,i3,kd)
        vz42(i1,i2,i3,kd)=0
        vx43(i1,i2,i3,kd)=rx(i1,i2,i3)*vr4(i1,i2,i3,kd)+sx(i1,i2,i3)*vs4(i1,i2,i3,kd)+tx(i1,i2,i3)*vt4(i1,i2,i3,kd)
        vy43(i1,i2,i3,kd)=ry(i1,i2,i3)*vr4(i1,i2,i3,kd)+sy(i1,i2,i3)*vs4(i1,i2,i3,kd)+ty(i1,i2,i3)*vt4(i1,i2,i3,kd)
        vz43(i1,i2,i3,kd)=rz(i1,i2,i3)*vr4(i1,i2,i3,kd)+sz(i1,i2,i3)*vs4(i1,i2,i3,kd)+tz(i1,i2,i3)*vt4(i1,i2,i3,kd)
        vxx41(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*vrr4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*vr4(i1,i2,i3,kd)
        vyy41(i1,i2,i3,kd)=0
        vxy41(i1,i2,i3,kd)=0
        vxz41(i1,i2,i3,kd)=0
        vyz41(i1,i2,i3,kd)=0
        vzz41(i1,i2,i3,kd)=0
        vlaplacian41(i1,i2,i3,kd)=vxx41(i1,i2,i3,kd)
        vxx42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*vrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*vss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3))*vr4(i1,i2,i3,kd)+(sxx42(i1,i2,i3))*vs4(i1,i2,i3,kd)
        vyy42(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*vrr4(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*vss4(i1,i2,i3,kd)+(ryy42(i1,i2,i3))*vr4(i1,i2,i3,kd)+(syy42(i1,i2,i3))*vs4(i1,i2,i3,kd)
        vxy42(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*vrr4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*vrs4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*vss4(i1,i2,i3,kd)+rxy42(i1,i2,i3)*vr4(i1,i2,i3,kd)+sxy42(i1,i2,i3)*vs4(i1,i2,i3,kd)
        vxz42(i1,i2,i3,kd)=0
        vyz42(i1,i2,i3,kd)=0
        vzz42(i1,i2,i3,kd)=0
        vlaplacian42(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*vrr4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*vss4(i1,i2,i3,kd)+(rxx42(i1,i2,i3)+ryy42(i1,i2,i3))*vr4(i1,i2,i3,kd)+(sxx42(i1,i2,i3)+syy42(i1,i2,i3))*vs4(i1,i2,i3,kd)
        vxx43(i1,i2,i3,kd)=rx(i1,i2,i3)**2*vrr4(i1,i2,i3,kd)+sx(i1,i2,i3)**2*vss4(i1,i2,i3,kd)+tx(i1,i2,i3)**2*vtt4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*vrs4(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*vrt4(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*vst4(i1,i2,i3,kd)+rxx43(i1,i2,i3)*vr4(i1,i2,i3,kd)+sxx43(i1,i2,i3)*vs4(i1,i2,i3,kd)+txx43(i1,i2,i3)*vt4(i1,i2,i3,kd)
        vyy43(i1,i2,i3,kd)=ry(i1,i2,i3)**2*vrr4(i1,i2,i3,kd)+sy(i1,i2,i3)**2*vss4(i1,i2,i3,kd)+ty(i1,i2,i3)**2*vtt4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*vrs4(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*vrt4(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*vst4(i1,i2,i3,kd)+ryy43(i1,i2,i3)*vr4(i1,i2,i3,kd)+syy43(i1,i2,i3)*vs4(i1,i2,i3,kd)+tyy43(i1,i2,i3)*vt4(i1,i2,i3,kd)
        vzz43(i1,i2,i3,kd)=rz(i1,i2,i3)**2*vrr4(i1,i2,i3,kd)+sz(i1,i2,i3)**2*vss4(i1,i2,i3,kd)+tz(i1,i2,i3)**2*vtt4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*vrs4(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*vrt4(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*vst4(i1,i2,i3,kd)+rzz43(i1,i2,i3)*vr4(i1,i2,i3,kd)+szz43(i1,i2,i3)*vs4(i1,i2,i3,kd)+tzz43(i1,i2,i3)*vt4(i1,i2,i3,kd)
        vxy43(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*vrr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*vss4(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*vtt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*vrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*vst4(i1,i2,i3,kd)+rxy43(i1,i2,i3)*vr4(i1,i2,i3,kd)+sxy43(i1,i2,i3)*vs4(i1,i2,i3,kd)+txy43(i1,i2,i3)*vt4(i1,i2,i3,kd)
        vxz43(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*vrr4(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*vss4(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*vtt4(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*vrt4(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*vst4(i1,i2,i3,kd)+rxz43(i1,i2,i3)*vr4(i1,i2,i3,kd)+sxz43(i1,i2,i3)*vs4(i1,i2,i3,kd)+txz43(i1,i2,i3)*vt4(i1,i2,i3,kd)
        vyz43(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*vrr4(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*vss4(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*vtt4(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*vrs4(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*vrt4(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*vst4(i1,i2,i3,kd)+ryz43(i1,i2,i3)*vr4(i1,i2,i3,kd)+syz43(i1,i2,i3)*vs4(i1,i2,i3,kd)+tyz43(i1,i2,i3)*vt4(i1,i2,i3,kd)
        vlaplacian43(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*vrr4(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*vss4(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*vtt4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*vrs4(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*vrt4(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*vst4(i1,i2,i3,kd)+(rxx43(i1,i2,i3)+ryy43(i1,i2,i3)+rzz43(i1,i2,i3))*vr4(i1,i2,i3,kd)+(sxx43(i1,i2,i3)+syy43(i1,i2,i3)+szz43(i1,i2,i3))*vs4(i1,i2,i3,kd)+(txx43(i1,i2,i3)+tyy43(i1,i2,i3)+tzz43(i1,i2,i3))*vt4(i1,i2,i3,kd)
    !============================================================================================
    ! Define derivatives for a rectangular grid
    !
    !============================================================================================
        vx43r(i1,i2,i3,kd)=(8.*(v(i1+1,i2,i3,kd)-v(i1-1,i2,i3,kd))-(v(i1+2,i2,i3,kd)-v(i1-2,i2,i3,kd)))*h41(0)
        vy43r(i1,i2,i3,kd)=(8.*(v(i1,i2+1,i3,kd)-v(i1,i2-1,i3,kd))-(v(i1,i2+2,i3,kd)-v(i1,i2-2,i3,kd)))*h41(1)
        vz43r(i1,i2,i3,kd)=(8.*(v(i1,i2,i3+1,kd)-v(i1,i2,i3-1,kd))-(v(i1,i2,i3+2,kd)-v(i1,i2,i3-2,kd)))*h41(2)
        vxx43r(i1,i2,i3,kd)=( -30.*v(i1,i2,i3,kd)+16.*(v(i1+1,i2,i3,kd)+v(i1-1,i2,i3,kd))-(v(i1+2,i2,i3,kd)+v(i1-2,i2,i3,kd)) )*h42(0) 
        vyy43r(i1,i2,i3,kd)=( -30.*v(i1,i2,i3,kd)+16.*(v(i1,i2+1,i3,kd)+v(i1,i2-1,i3,kd))-(v(i1,i2+2,i3,kd)+v(i1,i2-2,i3,kd)) )*h42(1) 
        vzz43r(i1,i2,i3,kd)=( -30.*v(i1,i2,i3,kd)+16.*(v(i1,i2,i3+1,kd)+v(i1,i2,i3-1,kd))-(v(i1,i2,i3+2,kd)+v(i1,i2,i3-2,kd)) )*h42(2)
        vxy43r(i1,i2,i3,kd)=( (v(i1+2,i2+2,i3,kd)-v(i1-2,i2+2,i3,kd)- v(i1+2,i2-2,i3,kd)+v(i1-2,i2-2,i3,kd)) +8.*(v(i1-1,i2+2,i3,kd)-v(i1-1,i2-2,i3,kd)-v(i1+1,i2+2,i3,kd)+v(i1+1,i2-2,i3,kd) +v(i1+2,i2-1,i3,kd)-v(i1-2,i2-1,i3,kd)-v(i1+2,i2+1,i3,kd)+v(i1-2,i2+1,i3,kd))+64.*(v(i1+1,i2+1,i3,kd)-v(i1-1,i2+1,i3,kd)- v(i1+1,i2-1,i3,kd)+v(i1-1,i2-1,i3,kd)))*(h41(0)*h41(1))
        vxz43r(i1,i2,i3,kd)=( (v(i1+2,i2,i3+2,kd)-v(i1-2,i2,i3+2,kd)-v(i1+2,i2,i3-2,kd)+v(i1-2,i2,i3-2,kd)) +8.*(v(i1-1,i2,i3+2,kd)-v(i1-1,i2,i3-2,kd)-v(i1+1,i2,i3+2,kd)+v(i1+1,i2,i3-2,kd) +v(i1+2,i2,i3-1,kd)-v(i1-2,i2,i3-1,kd)- v(i1+2,i2,i3+1,kd)+v(i1-2,i2,i3+1,kd)) +64.*(v(i1+1,i2,i3+1,kd)-v(i1-1,i2,i3+1,kd)-v(i1+1,i2,i3-1,kd)+v(i1-1,i2,i3-1,kd)) )*(h41(0)*h41(2))
        vyz43r(i1,i2,i3,kd)=( (v(i1,i2+2,i3+2,kd)-v(i1,i2-2,i3+2,kd)-v(i1,i2+2,i3-2,kd)+v(i1,i2-2,i3-2,kd)) +8.*(v(i1,i2-1,i3+2,kd)-v(i1,i2-1,i3-2,kd)-v(i1,i2+1,i3+2,kd)+v(i1,i2+1,i3-2,kd) +v(i1,i2+2,i3-1,kd)-v(i1,i2-2,i3-1,kd)-v(i1,i2+2,i3+1,kd)+v(i1,i2-2,i3+1,kd)) +64.*(v(i1,i2+1,i3+1,kd)-v(i1,i2-1,i3+1,kd)-v(i1,i2+1,i3-1,kd)+v(i1,i2-1,i3-1,kd)) )*(h41(1)*h41(2))
        vx41r(i1,i2,i3,kd)= vx43r(i1,i2,i3,kd)
        vy41r(i1,i2,i3,kd)= vy43r(i1,i2,i3,kd)
        vz41r(i1,i2,i3,kd)= vz43r(i1,i2,i3,kd)
        vxx41r(i1,i2,i3,kd)= vxx43r(i1,i2,i3,kd)
        vyy41r(i1,i2,i3,kd)= vyy43r(i1,i2,i3,kd)
        vzz41r(i1,i2,i3,kd)= vzz43r(i1,i2,i3,kd)
        vxy41r(i1,i2,i3,kd)= vxy43r(i1,i2,i3,kd)
        vxz41r(i1,i2,i3,kd)= vxz43r(i1,i2,i3,kd)
        vyz41r(i1,i2,i3,kd)= vyz43r(i1,i2,i3,kd)
        vlaplacian41r(i1,i2,i3,kd)=vxx43r(i1,i2,i3,kd)
        vx42r(i1,i2,i3,kd)= vx43r(i1,i2,i3,kd)
        vy42r(i1,i2,i3,kd)= vy43r(i1,i2,i3,kd)
        vz42r(i1,i2,i3,kd)= vz43r(i1,i2,i3,kd)
        vxx42r(i1,i2,i3,kd)= vxx43r(i1,i2,i3,kd)
        vyy42r(i1,i2,i3,kd)= vyy43r(i1,i2,i3,kd)
        vzz42r(i1,i2,i3,kd)= vzz43r(i1,i2,i3,kd)
        vxy42r(i1,i2,i3,kd)= vxy43r(i1,i2,i3,kd)
        vxz42r(i1,i2,i3,kd)= vxz43r(i1,i2,i3,kd)
        vyz42r(i1,i2,i3,kd)= vyz43r(i1,i2,i3,kd)
        vlaplacian42r(i1,i2,i3,kd)=vxx43r(i1,i2,i3,kd)+vyy43r(i1,i2,i3,kd)
        vlaplacian43r(i1,i2,i3,kd)=vxx43r(i1,i2,i3,kd)+vyy43r(i1,i2,i3,kd)+vzz43r(i1,i2,i3,kd)
        fr2(i1,i2,i3,kd)=(f(i1+1,i2,i3,kd)-f(i1-1,i2,i3,kd))*d12(0)
        fs2(i1,i2,i3,kd)=(f(i1,i2+1,i3,kd)-f(i1,i2-1,i3,kd))*d12(1)
        ft2(i1,i2,i3,kd)=(f(i1,i2,i3+1,kd)-f(i1,i2,i3-1,kd))*d12(2)
        frr2(i1,i2,i3,kd)=(-2.*f(i1,i2,i3,kd)+(f(i1+1,i2,i3,kd)+f(i1-1,i2,i3,kd)) )*d22(0)
        fss2(i1,i2,i3,kd)=(-2.*f(i1,i2,i3,kd)+(f(i1,i2+1,i3,kd)+f(i1,i2-1,i3,kd)) )*d22(1)
        frs2(i1,i2,i3,kd)=(fr2(i1,i2+1,i3,kd)-fr2(i1,i2-1,i3,kd))*d12(1)
        ftt2(i1,i2,i3,kd)=(-2.*f(i1,i2,i3,kd)+(f(i1,i2,i3+1,kd)+f(i1,i2,i3-1,kd)) )*d22(2)
        frt2(i1,i2,i3,kd)=(fr2(i1,i2,i3+1,kd)-fr2(i1,i2,i3-1,kd))*d12(2)
        fst2(i1,i2,i3,kd)=(fs2(i1,i2,i3+1,kd)-fs2(i1,i2,i3-1,kd))*d12(2)
        frrr2(i1,i2,i3,kd)=(-2.*(f(i1+1,i2,i3,kd)-f(i1-1,i2,i3,kd))+(f(i1+2,i2,i3,kd)-f(i1-2,i2,i3,kd)) )*d22(0)*d12(0)
        fsss2(i1,i2,i3,kd)=(-2.*(f(i1,i2+1,i3,kd)-f(i1,i2-1,i3,kd))+(f(i1,i2+2,i3,kd)-f(i1,i2-2,i3,kd)) )*d22(1)*d12(1)
        fttt2(i1,i2,i3,kd)=(-2.*(f(i1,i2,i3+1,kd)-f(i1,i2,i3-1,kd))+(f(i1,i2,i3+2,kd)-f(i1,i2,i3-2,kd)) )*d22(2)*d12(2)
        fx21(i1,i2,i3,kd)= rx(i1,i2,i3)*fr2(i1,i2,i3,kd)
        fy21(i1,i2,i3,kd)=0
        fz21(i1,i2,i3,kd)=0
        fx22(i1,i2,i3,kd)= rx(i1,i2,i3)*fr2(i1,i2,i3,kd)+sx(i1,i2,i3)*fs2(i1,i2,i3,kd)
        fy22(i1,i2,i3,kd)= ry(i1,i2,i3)*fr2(i1,i2,i3,kd)+sy(i1,i2,i3)*fs2(i1,i2,i3,kd)
        fz22(i1,i2,i3,kd)=0
        fx23(i1,i2,i3,kd)=rx(i1,i2,i3)*fr2(i1,i2,i3,kd)+sx(i1,i2,i3)*fs2(i1,i2,i3,kd)+tx(i1,i2,i3)*ft2(i1,i2,i3,kd)
        fy23(i1,i2,i3,kd)=ry(i1,i2,i3)*fr2(i1,i2,i3,kd)+sy(i1,i2,i3)*fs2(i1,i2,i3,kd)+ty(i1,i2,i3)*ft2(i1,i2,i3,kd)
        fz23(i1,i2,i3,kd)=rz(i1,i2,i3)*fr2(i1,i2,i3,kd)+sz(i1,i2,i3)*fs2(i1,i2,i3,kd)+tz(i1,i2,i3)*ft2(i1,i2,i3,kd)
        fxx21(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*frr2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*fr2(i1,i2,i3,kd)
        fyy21(i1,i2,i3,kd)=0
        fxy21(i1,i2,i3,kd)=0
        fxz21(i1,i2,i3,kd)=0
        fyz21(i1,i2,i3,kd)=0
        fzz21(i1,i2,i3,kd)=0
        flaplacian21(i1,i2,i3,kd)=fxx21(i1,i2,i3,kd)
        fxx22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2)*frr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3))*frs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2)*fss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3))*fr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3))*fs2(i1,i2,i3,kd)
        fyy22(i1,i2,i3,kd)=(ry(i1,i2,i3)**2)*frr2(i1,i2,i3,kd)+2.*(ry(i1,i2,i3)*sy(i1,i2,i3))*frs2(i1,i2,i3,kd)+(sy(i1,i2,i3)**2)*fss2(i1,i2,i3,kd)+(ryy22(i1,i2,i3))*fr2(i1,i2,i3,kd)+(syy22(i1,i2,i3))*fs2(i1,i2,i3,kd)
        fxy22(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*frr2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*frs2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*fss2(i1,i2,i3,kd)+rxy22(i1,i2,i3)*fr2(i1,i2,i3,kd)+sxy22(i1,i2,i3)*fs2(i1,i2,i3,kd)
        fxz22(i1,i2,i3,kd)=0
        fyz22(i1,i2,i3,kd)=0
        fzz22(i1,i2,i3,kd)=0
        flaplacian22(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2)*frr2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3))*frs2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2)*fss2(i1,i2,i3,kd)+(rxx22(i1,i2,i3)+ryy22(i1,i2,i3))*fr2(i1,i2,i3,kd)+(sxx22(i1,i2,i3)+syy22(i1,i2,i3))*fs2(i1,i2,i3,kd)
        fxx23(i1,i2,i3,kd)=rx(i1,i2,i3)**2*frr2(i1,i2,i3,kd)+sx(i1,i2,i3)**2*fss2(i1,i2,i3,kd)+tx(i1,i2,i3)**2*ftt2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*sx(i1,i2,i3)*frs2(i1,i2,i3,kd)+2.*rx(i1,i2,i3)*tx(i1,i2,i3)*frt2(i1,i2,i3,kd)+2.*sx(i1,i2,i3)*tx(i1,i2,i3)*fst2(i1,i2,i3,kd)+rxx23(i1,i2,i3)*fr2(i1,i2,i3,kd)+sxx23(i1,i2,i3)*fs2(i1,i2,i3,kd)+txx23(i1,i2,i3)*ft2(i1,i2,i3,kd)
        fyy23(i1,i2,i3,kd)=ry(i1,i2,i3)**2*frr2(i1,i2,i3,kd)+sy(i1,i2,i3)**2*fss2(i1,i2,i3,kd)+ty(i1,i2,i3)**2*ftt2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*sy(i1,i2,i3)*frs2(i1,i2,i3,kd)+2.*ry(i1,i2,i3)*ty(i1,i2,i3)*frt2(i1,i2,i3,kd)+2.*sy(i1,i2,i3)*ty(i1,i2,i3)*fst2(i1,i2,i3,kd)+ryy23(i1,i2,i3)*fr2(i1,i2,i3,kd)+syy23(i1,i2,i3)*fs2(i1,i2,i3,kd)+tyy23(i1,i2,i3)*ft2(i1,i2,i3,kd)
        fzz23(i1,i2,i3,kd)=rz(i1,i2,i3)**2*frr2(i1,i2,i3,kd)+sz(i1,i2,i3)**2*fss2(i1,i2,i3,kd)+tz(i1,i2,i3)**2*ftt2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*sz(i1,i2,i3)*frs2(i1,i2,i3,kd)+2.*rz(i1,i2,i3)*tz(i1,i2,i3)*frt2(i1,i2,i3,kd)+2.*sz(i1,i2,i3)*tz(i1,i2,i3)*fst2(i1,i2,i3,kd)+rzz23(i1,i2,i3)*fr2(i1,i2,i3,kd)+szz23(i1,i2,i3)*fs2(i1,i2,i3,kd)+tzz23(i1,i2,i3)*ft2(i1,i2,i3,kd)
        fxy23(i1,i2,i3,kd)=rx(i1,i2,i3)*ry(i1,i2,i3)*frr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sy(i1,i2,i3)*fss2(i1,i2,i3,kd)+tx(i1,i2,i3)*ty(i1,i2,i3)*ftt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sy(i1,i2,i3)+ry(i1,i2,i3)*sx(i1,i2,i3))*frs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*ty(i1,i2,i3)+ry(i1,i2,i3)*tx(i1,i2,i3))*frt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*ty(i1,i2,i3)+sy(i1,i2,i3)*tx(i1,i2,i3))*fst2(i1,i2,i3,kd)+rxy23(i1,i2,i3)*fr2(i1,i2,i3,kd)+sxy23(i1,i2,i3)*fs2(i1,i2,i3,kd)+txy23(i1,i2,i3)*ft2(i1,i2,i3,kd)
        fxz23(i1,i2,i3,kd)=rx(i1,i2,i3)*rz(i1,i2,i3)*frr2(i1,i2,i3,kd)+sx(i1,i2,i3)*sz(i1,i2,i3)*fss2(i1,i2,i3,kd)+tx(i1,i2,i3)*tz(i1,i2,i3)*ftt2(i1,i2,i3,kd)+(rx(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sx(i1,i2,i3))*frs2(i1,i2,i3,kd)+(rx(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*tx(i1,i2,i3))*frt2(i1,i2,i3,kd)+(sx(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*tx(i1,i2,i3))*fst2(i1,i2,i3,kd)+rxz23(i1,i2,i3)*fr2(i1,i2,i3,kd)+sxz23(i1,i2,i3)*fs2(i1,i2,i3,kd)+txz23(i1,i2,i3)*ft2(i1,i2,i3,kd)
        fyz23(i1,i2,i3,kd)=ry(i1,i2,i3)*rz(i1,i2,i3)*frr2(i1,i2,i3,kd)+sy(i1,i2,i3)*sz(i1,i2,i3)*fss2(i1,i2,i3,kd)+ty(i1,i2,i3)*tz(i1,i2,i3)*ftt2(i1,i2,i3,kd)+(ry(i1,i2,i3)*sz(i1,i2,i3)+rz(i1,i2,i3)*sy(i1,i2,i3))*frs2(i1,i2,i3,kd)+(ry(i1,i2,i3)*tz(i1,i2,i3)+rz(i1,i2,i3)*ty(i1,i2,i3))*frt2(i1,i2,i3,kd)+(sy(i1,i2,i3)*tz(i1,i2,i3)+sz(i1,i2,i3)*ty(i1,i2,i3))*fst2(i1,i2,i3,kd)+ryz23(i1,i2,i3)*fr2(i1,i2,i3,kd)+syz23(i1,i2,i3)*fs2(i1,i2,i3,kd)+tyz23(i1,i2,i3)*ft2(i1,i2,i3,kd)
        flaplacian23(i1,i2,i3,kd)=(rx(i1,i2,i3)**2+ry(i1,i2,i3)**2+rz(i1,i2,i3)**2)*frr2(i1,i2,i3,kd)+(sx(i1,i2,i3)**2+sy(i1,i2,i3)**2+sz(i1,i2,i3)**2)*fss2(i1,i2,i3,kd)+(tx(i1,i2,i3)**2+ty(i1,i2,i3)**2+tz(i1,i2,i3)**2)*ftt2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*sx(i1,i2,i3)+ ry(i1,i2,i3)*sy(i1,i2,i3)+rz(i1,i2,i3)*sz(i1,i2,i3))*frs2(i1,i2,i3,kd)+2.*(rx(i1,i2,i3)*tx(i1,i2,i3)+ ry(i1,i2,i3)*ty(i1,i2,i3)+rz(i1,i2,i3)*tz(i1,i2,i3))*frt2(i1,i2,i3,kd)+2.*(sx(i1,i2,i3)*tx(i1,i2,i3)+ sy(i1,i2,i3)*ty(i1,i2,i3)+sz(i1,i2,i3)*tz(i1,i2,i3))*fst2(i1,i2,i3,kd)+(rxx23(i1,i2,i3)+ryy23(i1,i2,i3)+rzz23(i1,i2,i3))*fr2(i1,i2,i3,kd)+(sxx23(i1,i2,i3)+syy23(i1,i2,i3)+szz23(i1,i2,i3))*fs2(i1,i2,i3,kd)+(txx23(i1,i2,i3)+tyy23(i1,i2,i3)+tzz23(i1,i2,i3))*ft2(i1,i2,i3,kd)
    !============================================================================================
    ! Define derivatives for a rectangular grid
    !
    !============================================================================================
        fx23r(i1,i2,i3,kd)=(f(i1+1,i2,i3,kd)-f(i1-1,i2,i3,kd))*h12(0)
        fy23r(i1,i2,i3,kd)=(f(i1,i2+1,i3,kd)-f(i1,i2-1,i3,kd))*h12(1)
        fz23r(i1,i2,i3,kd)=(f(i1,i2,i3+1,kd)-f(i1,i2,i3-1,kd))*h12(2)
        fxx23r(i1,i2,i3,kd)=(-2.*f(i1,i2,i3,kd)+(f(i1+1,i2,i3,kd)+f(i1-1,i2,i3,kd)) )*h22(0)
        fyy23r(i1,i2,i3,kd)=(-2.*f(i1,i2,i3,kd)+(f(i1,i2+1,i3,kd)+f(i1,i2-1,i3,kd)) )*h22(1)
        fxy23r(i1,i2,i3,kd)=(fx23r(i1,i2+1,i3,kd)-fx23r(i1,i2-1,i3,kd))*h12(1)
        fzz23r(i1,i2,i3,kd)=(-2.*f(i1,i2,i3,kd)+(f(i1,i2,i3+1,kd)+f(i1,i2,i3-1,kd)) )*h22(2)
        fxz23r(i1,i2,i3,kd)=(fx23r(i1,i2,i3+1,kd)-fx23r(i1,i2,i3-1,kd))*h12(2)
        fyz23r(i1,i2,i3,kd)=(fy23r(i1,i2,i3+1,kd)-fy23r(i1,i2,i3-1,kd))*h12(2)
        fx21r(i1,i2,i3,kd)= fx23r(i1,i2,i3,kd)
        fy21r(i1,i2,i3,kd)= fy23r(i1,i2,i3,kd)
        fz21r(i1,i2,i3,kd)= fz23r(i1,i2,i3,kd)
        fxx21r(i1,i2,i3,kd)= fxx23r(i1,i2,i3,kd)
        fyy21r(i1,i2,i3,kd)= fyy23r(i1,i2,i3,kd)
        fzz21r(i1,i2,i3,kd)= fzz23r(i1,i2,i3,kd)
        fxy21r(i1,i2,i3,kd)= fxy23r(i1,i2,i3,kd)
        fxz21r(i1,i2,i3,kd)= fxz23r(i1,i2,i3,kd)
        fyz21r(i1,i2,i3,kd)= fyz23r(i1,i2,i3,kd)
        flaplacian21r(i1,i2,i3,kd)=fxx23r(i1,i2,i3,kd)
        fx22r(i1,i2,i3,kd)= fx23r(i1,i2,i3,kd)
        fy22r(i1,i2,i3,kd)= fy23r(i1,i2,i3,kd)
        fz22r(i1,i2,i3,kd)= fz23r(i1,i2,i3,kd)
        fxx22r(i1,i2,i3,kd)= fxx23r(i1,i2,i3,kd)
        fyy22r(i1,i2,i3,kd)= fyy23r(i1,i2,i3,kd)
        fzz22r(i1,i2,i3,kd)= fzz23r(i1,i2,i3,kd)
        fxy22r(i1,i2,i3,kd)= fxy23r(i1,i2,i3,kd)
        fxz22r(i1,i2,i3,kd)= fxz23r(i1,i2,i3,kd)
        fyz22r(i1,i2,i3,kd)= fyz23r(i1,i2,i3,kd)
        flaplacian22r(i1,i2,i3,kd)=fxx23r(i1,i2,i3,kd)+fyy23r(i1,i2,i3,kd)
        flaplacian23r(i1,i2,i3,kd)=fxx23r(i1,i2,i3,kd)+fyy23r(i1,i2,i3,kd)+fzz23r(i1,i2,i3,kd)
        fxxx22r(i1,i2,i3,kd)=(-2.*(f(i1+1,i2,i3,kd)-f(i1-1,i2,i3,kd))+(f(i1+2,i2,i3,kd)-f(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
        fyyy22r(i1,i2,i3,kd)=(-2.*(f(i1,i2+1,i3,kd)-f(i1,i2-1,i3,kd))+(f(i1,i2+2,i3,kd)-f(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
        fxxy22r(i1,i2,i3,kd)=( fxx22r(i1,i2+1,i3,kd)-fxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
        fxyy22r(i1,i2,i3,kd)=( fyy22r(i1+1,i2,i3,kd)-fyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
        fxxxx22r(i1,i2,i3,kd)=(6.*f(i1,i2,i3,kd)-4.*(f(i1+1,i2,i3,kd)+f(i1-1,i2,i3,kd))+(f(i1+2,i2,i3,kd)+f(i1-2,i2,i3,kd)) )/(dx(0)**4)
        fyyyy22r(i1,i2,i3,kd)=(6.*f(i1,i2,i3,kd)-4.*(f(i1,i2+1,i3,kd)+f(i1,i2-1,i3,kd))+(f(i1,i2+2,i3,kd)+f(i1,i2-2,i3,kd)) )/(dx(1)**4)
        fxxyy22r(i1,i2,i3,kd)=( 4.*f(i1,i2,i3,kd)     -2.*(f(i1+1,i2,i3,kd)+f(i1-1,i2,i3,kd)+f(i1,i2+1,i3,kd)+f(i1,i2-1,i3,kd))   +   (f(i1+1,i2+1,i3,kd)+f(i1-1,i2+1,i3,kd)+f(i1+1,i2-1,i3,kd)+f(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
    ! 2D laplacian squared = f.xxxx + 2 f.xxyy + f.yyyy
        fLapSq22r(i1,i2,i3,kd)= ( 6.*f(i1,i2,i3,kd)   - 4.*(f(i1+1,i2,i3,kd)+f(i1-1,i2,i3,kd))    +(f(i1+2,i2,i3,kd)+f(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*f(i1,i2,i3,kd)    -4.*(f(i1,i2+1,i3,kd)+f(i1,i2-1,i3,kd))    +(f(i1,i2+2,i3,kd)+f(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 8.*f(i1,i2,i3,kd)     -4.*(f(i1+1,i2,i3,kd)+f(i1-1,i2,i3,kd)+f(i1,i2+1,i3,kd)+f(i1,i2-1,i3,kd))   +2.*(f(i1+1,i2+1,i3,kd)+f(i1-1,i2+1,i3,kd)+f(i1+1,i2-1,i3,kd)+f(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
        fxxx23r(i1,i2,i3,kd)=(-2.*(f(i1+1,i2,i3,kd)-f(i1-1,i2,i3,kd))+(f(i1+2,i2,i3,kd)-f(i1-2,i2,i3,kd)) )*h22(0)*h12(0)
        fyyy23r(i1,i2,i3,kd)=(-2.*(f(i1,i2+1,i3,kd)-f(i1,i2-1,i3,kd))+(f(i1,i2+2,i3,kd)-f(i1,i2-2,i3,kd)) )*h22(1)*h12(1)
        fzzz23r(i1,i2,i3,kd)=(-2.*(f(i1,i2,i3+1,kd)-f(i1,i2,i3-1,kd))+(f(i1,i2,i3+2,kd)-f(i1,i2,i3-2,kd)) )*h22(1)*h12(2)
        fxxy23r(i1,i2,i3,kd)=( fxx22r(i1,i2+1,i3,kd)-fxx22r(i1,i2-1,i3,kd))/(2.*dx(1))
        fxyy23r(i1,i2,i3,kd)=( fyy22r(i1+1,i2,i3,kd)-fyy22r(i1-1,i2,i3,kd))/(2.*dx(0))
        fxxz23r(i1,i2,i3,kd)=( fxx22r(i1,i2,i3+1,kd)-fxx22r(i1,i2,i3-1,kd))/(2.*dx(2))
        fyyz23r(i1,i2,i3,kd)=( fyy22r(i1,i2,i3+1,kd)-fyy22r(i1,i2,i3-1,kd))/(2.*dx(2))
        fxzz23r(i1,i2,i3,kd)=( fzz22r(i1+1,i2,i3,kd)-fzz22r(i1-1,i2,i3,kd))/(2.*dx(0))
        fyzz23r(i1,i2,i3,kd)=( fzz22r(i1,i2+1,i3,kd)-fzz22r(i1,i2-1,i3,kd))/(2.*dx(1))
        fxxxx23r(i1,i2,i3,kd)=(6.*f(i1,i2,i3,kd)-4.*(f(i1+1,i2,i3,kd)+f(i1-1,i2,i3,kd))+(f(i1+2,i2,i3,kd)+f(i1-2,i2,i3,kd)) )/(dx(0)**4)
        fyyyy23r(i1,i2,i3,kd)=(6.*f(i1,i2,i3,kd)-4.*(f(i1,i2+1,i3,kd)+f(i1,i2-1,i3,kd))+(f(i1,i2+2,i3,kd)+f(i1,i2-2,i3,kd)) )/(dx(1)**4)
        fzzzz23r(i1,i2,i3,kd)=(6.*f(i1,i2,i3,kd)-4.*(f(i1,i2,i3+1,kd)+f(i1,i2,i3-1,kd))+(f(i1,i2,i3+2,kd)+f(i1,i2,i3-2,kd)) )/(dx(2)**4)
        fxxyy23r(i1,i2,i3,kd)=( 4.*f(i1,i2,i3,kd)     -2.*(f(i1+1,i2,i3,kd)+f(i1-1,i2,i3,kd)+f(i1,i2+1,i3,kd)+f(i1,i2-1,i3,kd))   +   (f(i1+1,i2+1,i3,kd)+f(i1-1,i2+1,i3,kd)+f(i1+1,i2-1,i3,kd)+f(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)
        fxxzz23r(i1,i2,i3,kd)=( 4.*f(i1,i2,i3,kd)     -2.*(f(i1+1,i2,i3,kd)+f(i1-1,i2,i3,kd)+f(i1,i2,i3+1,kd)+f(i1,i2,i3-1,kd))   +   (f(i1+1,i2,i3+1,kd)+f(i1-1,i2,i3+1,kd)+f(i1+1,i2,i3-1,kd)+f(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)
        fyyzz23r(i1,i2,i3,kd)=( 4.*f(i1,i2,i3,kd)     -2.*(f(i1,i2+1,i3,kd)  +f(i1,i2-1,i3,kd)+  f(i1,i2  ,i3+1,kd)+f(i1,i2  ,i3-1,kd))   +   (f(i1,i2+1,i3+1,kd)+f(i1,i2-1,i3+1,kd)+f(i1,i2+1,i3-1,kd)+f(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
    ! 3D laplacian squared = f.xxxx + f.yyyy + f.zzzz + 2 (f.xxyy + f.xxzz + f.yyzz )
        fLapSq23r(i1,i2,i3,kd)= ( 6.*f(i1,i2,i3,kd)   - 4.*(f(i1+1,i2,i3,kd)+f(i1-1,i2,i3,kd))    +(f(i1+2,i2,i3,kd)+f(i1-2,i2,i3,kd)) )/(dx(0)**4) +( 6.*f(i1,i2,i3,kd)    -4.*(f(i1,i2+1,i3,kd)+f(i1,i2-1,i3,kd))    +(f(i1,i2+2,i3,kd)+f(i1,i2-2,i3,kd)) )/(dx(1)**4)  +( 6.*f(i1,i2,i3,kd)    -4.*(f(i1,i2,i3+1,kd)+f(i1,i2,i3-1,kd))    +(f(i1,i2,i3+2,kd)+f(i1,i2,i3-2,kd)) )/(dx(2)**4)  +( 8.*f(i1,i2,i3,kd)     -4.*(f(i1+1,i2,i3,kd)  +f(i1-1,i2,i3,kd)  +f(i1  ,i2+1,i3,kd)+f(i1  ,i2-1,i3,kd))   +2.*(f(i1+1,i2+1,i3,kd)+f(i1-1,i2+1,i3,kd)+f(i1+1,i2-1,i3,kd)+f(i1-1,i2-1,i3,kd)) )/(dx(0)**2*dx(1)**2)+( 8.*f(i1,i2,i3,kd)     -4.*(f(i1+1,i2,i3,kd)  +f(i1-1,i2,i3,kd)  +f(i1  ,i2,i3+1,kd)+f(i1  ,i2,i3-1,kd))   +2.*(f(i1+1,i2,i3+1,kd)+f(i1-1,i2,i3+1,kd)+f(i1+1,i2,i3-1,kd)+f(i1-1,i2,i3-1,kd)) )/(dx(0)**2*dx(2)**2)+( 8.*f(i1,i2,i3,kd)     -4.*(f(i1,i2+1,i3,kd)  +f(i1,i2-1,i3,kd)  +f(i1,i2  ,i3+1,kd)+f(i1,i2  ,i3-1,kd))   +2.*(f(i1,i2+1,i3+1,kd)+f(i1,i2-1,i3+1,kd)+f(i1,i2+1,i3-1,kd)+f(i1,i2-1,i3-1,kd)) )/(dx(1)**2*dx(2)**2)
     ! 2D laplacian squared = u.xxxx + 2 u.xxyy + u.yyyy
          lap2d2Pow2(i1,i2,i3,c)= ( 6.*u(i1,i2,i3,c)   - 4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))    +(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*dxi4 +( 6.*u(i1,i2,i3,c)    -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))    +(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dyi4  +( 8.*u(i1,i2,i3,c)     -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))   +2.*(u(i1+1,i2+1,i3,c)+u(i1-1,i2+1,i3,c)+u(i1+1,i2-1,i3,c)+u(i1-1,i2-1,i3,c)) )*dxdyi2
     ! 3D laplacian squared = u.xxxx + u.yyyy + u.zzzz + 2 (u.xxyy + u.xxzz + u.yyzz )
          lap3d2Pow2(i1,i2,i3,c)= ( 6.*u(i1,i2,i3,c)   - 4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))    +(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*dxi4 +(  +6.*u(i1,i2,i3,c)    -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))    +(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dyi4+(  +6.*u(i1,i2,i3,c)    -4.*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))    +(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) )*dzi4+(8.*u(i1,i2,i3,c)     -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))   +2.*(u(i1+1,i2+1,i3,c)+u(i1-1,i2+1,i3,c)+u(i1+1,i2-1,i3,c)+u(i1-1,i2-1,i3,c)) )*dxdyi2 +(8.*u(i1,i2,i3,c)     -4.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c)+u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))   +2.*(u(i1+1,i2,i3+1,c)+u(i1-1,i2,i3+1,c)+u(i1+1,i2,i3-1,c)+u(i1-1,i2,i3-1,c)) )*dxdzi2 +(8.*u(i1,i2,i3,c)     -4.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c)+u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))   +2.*(u(i1,i2+1,i3+1,c)+u(i1,i2-1,i3+1,c)+u(i1,i2+1,i3-1,c)+u(i1,i2-1,i3-1,c)) )*dydzi2 
   !    ** 4th order ****
          lap2d4(i1,i2,i3,c)=( -30.*u(i1,i2,i3,c)     +16.*(u(i1+1,i2,i3,c)+u(i1-1,i2,i3,c))     -(u(i1+2,i2,i3,c)+u(i1-2,i2,i3,c)) )*dxsq12i + ( -30.*u(i1,i2,i3,c)     +16.*(u(i1,i2+1,i3,c)+u(i1,i2-1,i3,c))     -(u(i1,i2+2,i3,c)+u(i1,i2-2,i3,c)) )*dysq12i 
          lap3d4(i1,i2,i3,c)=lap2d4(i1,i2,i3,c)+ ( -30.*u(i1,i2,i3,c)      +16.*(u(i1,i2,i3+1,c)+u(i1,i2,i3-1,c))      -(u(i1,i2,i3+2,c)+u(i1,i2,i3-2,c)) )*dzsq12i 
          lap2d4Pow2(i1,i2,i3,c)=(-30.*lap2d4(i1,i2,i3,c)+16.*(lap2d4(i1+1,i2,i3,c)+lap2d4(i1-1,i2,i3,c))-(lap2d4(i1+2,i2,i3,c)+lap2d4(i1-2,i2,i3,c)))*dxsq12i+(-30.*lap2d4(i1,i2,i3,c)+16.*(lap2d4(i1,i2+1,i3,c)+lap2d4(i1,i2-1,i3,c))-(lap2d4(i1,i2+2,i3,c)+lap2d4(i1,i2-2,i3,c)))*dysq12i
          lap3d4Pow2(i1,i2,i3,c)=(-30.*lap3d4(i1,i2,i3,c)+16.*(lap3d4(i1+1,i2,i3,c)+lap3d4(i1-1,i2,i3,c))-(lap3d4(i1+2,i2,i3,c)+lap3d4(i1-2,i2,i3,c)))*dxsq12i+(-30.*lap3d4(i1,i2,i3,c)+16.*(lap3d4(i1,i2+1,i3,c)+lap3d4(i1,i2-1,i3,c))-(lap3d4(i1,i2+2,i3,c)+lap3d4(i1,i2-2,i3,c)))*dysq12i+(-30.*lap3d4(i1,i2,i3,c)+16.*(lap3d4(i1,i2,i3+1,c)+lap3d4(i1,i2,i3-1,c))-(lap3d4(i1,i2,i3+2,c)+lap3d4(i1,i2,i3-2,c)))*dzsq12i
          lap2d4Pow3(i1,i2,i3,c)=(-30.*lap2d4Pow2(i1,i2,i3,c)+16.*(lap2d4Pow2(i1+1,i2,i3,c)+lap2d4Pow2(i1-1,i2,i3,c))-(lap2d4Pow2(i1+2,i2,i3,c)+lap2d4Pow2(i1-2,i2,i3,c)))*dxsq12i+(-30.*lap2d4Pow2(i1,i2,i3,c)+16.*(lap2d4Pow2(i1,i2+1,i3,c)+lap2d4Pow2(i1,i2-1,i3,c))-(lap2d4Pow2(i1,i2+2,i3,c)+lap2d4Pow2(i1,i2-2,i3,c)))*dysq12i
          lap3d4Pow3(i1,i2,i3,c)=(-30.*lap3d4Pow2(i1,i2,i3,c)+16.*(lap3d4Pow2(i1+1,i2,i3,c)+lap3d4Pow2(i1-1,i2,i3,c))-(lap3d4Pow2(i1+2,i2,i3,c)+lap3d4Pow2(i1-2,i2,i3,c)))*dxsq12i+(-30.*lap3d4Pow2(i1,i2,i3,c)+16.*(lap3d4Pow2(i1,i2+1,i3,c)+lap3d4Pow2(i1,i2-1,i3,c))-(lap3d4Pow2(i1,i2+2,i3,c)+lap3d4Pow2(i1,i2-2,i3,c)))*dysq12i+(-30.*lap3d4Pow2(i1,i2,i3,c)+16.*(lap3d4Pow2(i1,i2,i3+1,c)+lap3d4Pow2(i1,i2,i3-1,c))-(lap3d4Pow2(i1,i2,i3+2,c)+lap3d4Pow2(i1,i2,i3-2,c)))*dzsq12i
     ! D-zero in time (really undivided)
          DztU(i1,i2,i3,n) = (un(i1,i2,i3,n)-um(i1,i2,i3,n))
     !...........end   statement functions
     ! write(*,*) 'Inside advWave...'
          cc             = rpar( 0)  ! this is c
          dt             = rpar( 1)
          dx(0)          = rpar( 2)
          dx(1)          = rpar( 3)
          dx(2)          = rpar( 4)
          dr(0)          = rpar( 5)
          dr(1)          = rpar( 6)
          dr(2)          = rpar( 7)
          t              = rpar( 8)
          ep             = rpar( 9)
          sosupParameter = rpar(10)
          omega          = rpar(11) ! for helmholtz 
          bImp( 0)       = rpar(12) ! beta2 : coefficient for implicit time-stepping
          bImp( 1)       = rpar(13) ! beta4 : coefficient for implicit time-stepping
          bImp( 2)       = rpar(14) ! beta6 (for future)
          bImp( 3)       = rpar(15) ! beta8 (for future)
          gridCFL        = rpar(16)
          dy=dx(1)  ! Are these needed?
          dz=dx(2)
          option                       = ipar( 0)
          grid                         = ipar( 1)
          gridType                     = ipar( 2)
          orderOfAccuracy              = ipar( 3)
          orderInTime                  = ipar( 4)
          addForcing                   = ipar( 5)
          forcingOption                = ipar( 6)
          numberOfForcingFunctions     = ipar( 7)
          fcur                         = ipar( 8) 
          debug                        = ipar( 9)
          gridIsImplicit               = ipar(10)
          useUpwindDissipation         = ipar(11)  ! explicit upwind dissipation
          useImplicitUpwindDissipation = ipar(12)  ! true if upwind-dissipation is on for implicit time-stepping
          preComputeUpwindUt           = ipar(13)
          numberOfFrequencies          = ipar(14)
          adjustOmega                  = ipar(15)
          solveHelmholtz               = ipar(16)
          adjustHelmholtzForUpwinding  = ipar(17)
          modifiedEquationApproach     = ipar(18)
          fprev = mod(fcur-1+numberOfForcingFunctions,max(1,numberOfForcingFunctions))
          fnext = mod(fcur+1                         ,max(1,numberOfForcingFunctions))
     ! ** fix me ***
          timeSteppingMethod=modifiedEquationTimeStepping
     ! Set dr(:) = dx(:) for 6th-order derivatives
          if( gridType.eq.rectangular )then
              do axis=0,2
                  dr(axis)=dx(axis)
              end do
          end if
     ! ---- Compute the coefficients in the implicit time-stepping scheme ----
          beta2=bImp(0)
          beta4=bImp(1)
          alpha2 = (1.-beta2)/2.
          alpha4 = (alpha2-beta4-1./12.)/2. 
          cImp(-1,0)=alpha2
          cImp( 0,0)= beta2
          cImp( 1,0)=alpha2
          cImp(-1,1)=alpha4
          cImp( 0,1)= beta4
          cImp( 1,1)=alpha4  
     ! ! addDissipation=.true. if we add the dissipation in the dis(i1,i2,i3,c) array
     ! !  if combineDissipationWithAdvance.ne.0 we compute the dissipation on the fly in the time step
     ! !  rather than pre-computing it in diss(i1,i2,i3,c)
     ! addDissipation = adc.gt.0.
     ! adcdt=adc*dt
          csq=cc**2
          dtsq=dt**2
          cdtsq=(cc**2)*(dt**2)
          cdt=cc*dt
     ! new: 
          cdtPow4By12  = cdt**4/12.
          cdtPow6By360 = cdt**6/360 
          cdtsq12=cdtsq*cdtsq/12.  ! c^4 dt^4 /12 
     ! cdt4by360=(cdt)**4/360.  ! (c*dt)^4/360 
     ! cdt6by20160=cdt**6/(8.*7.*6.*5.*4.*3.)
          cdtSqBy12= cdtsq/12.   ! c^2*dt*2/12
          dt4by12=dtsq*dtsq/12.
          cdtdx = (cc*dt/dx(0))**2
          cdtdy = (cc*dt/dy)**2
          cdtdz = (cc*dt/dz)**2
          dxsqi=1./(dx(0)**2)
          dysqi=1./(dy**2)
          dzsqi=1./(dz**2)
          dxsq12i=1./(12.*dx(0)**2)
          dysq12i=1./(12.*dy**2)
          dzsq12i=1./(12.*dz**2)
          dxi4=1./(dx(0)**4)
          dyi4=1./(dy**4)
          dxdyi2=1./(dx(0)*dx(0)*dy*dy)
          dzi4=1./(dz**4)
          dxdzi2=1./(dx(0)*dx(0)*dz*dz)
          dydzi2=1./(dy*dy*dz*dz)
          if( option.eq.1 )then 
            useSosupDissipation = 1
          else
            useSosupDissipation = 0
          end if
            if( useImplicitUpwindDissipation.eq.1 )then
       ! Upwind dissipation for implicit time time-stepping
       ! from formImplicitTimeSteppingMatrix
       ! // fourth-order dissipation for 2nd-order scheme:
       ! Real upwindCoeff4[4][5] = { 1.,-4.,6.,-4.,1.,
       !                             1.,-3.,3.,-1.,0.,   // extrap right-most point D-^3 u(2)
       !                             0.,-1.,3.,-3.,1.,   // extrap left -most point D+^3 u(-2)
       !                             0.,-1.,2.,-1.,0.
       !                           };
       ! // sixth-order dissipation for 4th-order scheme
       ! Real upwindCoeff6[4][7] = {1.,-6.,15.,-20.,15.,-6.,1.,
       !                            1.,-5.,10.,-10., 5.,-1.,0.,  // extrap right-most point D-^5 u(3)
       !                            0.,-1., 5.,-10.,10.,-5.,1.,  // extrap left -most point D+^5 u(-3)
       !                            0.,-1., 4., -6., 4.,-1.,0.
       !                           };
              if( orderOfAccuracy.eq.2 )then
         ! UPW stencils for the 4 cases (depends on whether all points exist in the wider stencil )
                  upwindCoeff(-2,0)=1.; upwindCoeff(-1,0)=-4.; upwindCoeff(0,0)=6.; upwindCoeff(1,0)=-4.; upwindCoeff(2,0)=1.   ! centred
                  upwindCoeff(-2,1)=1.; upwindCoeff(-1,1)=-3.; upwindCoeff(0,1)=3.; upwindCoeff(1,1)=-1.; upwindCoeff(2,1)=0.   ! left biased
                  upwindCoeff(-2,2)=0.; upwindCoeff(-1,2)=-1.; upwindCoeff(0,2)=3.; upwindCoeff(1,2)=-3.; upwindCoeff(2,2)=1.   ! right biased
                  upwindCoeff(-2,3)=0.; upwindCoeff(-1,3)=-1.; upwindCoeff(0,3)=2.; upwindCoeff(1,3)=-1.; upwindCoeff(2,3)=0.   ! centred, lower-order
              else if( orderOfAccuracy.eq.4 )then
                  upwindCoeff(-3,0)=1.; upwindCoeff(-2,0)=-6.; upwindCoeff(-1,0)=15.; upwindCoeff(0,0)=-20.; upwindCoeff(1,0)=15.; upwindCoeff(2,0)=-6.; upwindCoeff(3,0)=1.
                  upwindCoeff(-3,1)=1.; upwindCoeff(-2,1)=-5.; upwindCoeff(-1,1)=10.; upwindCoeff(0,1)=-10.; upwindCoeff(1,1)= 5.; upwindCoeff(2,1)=-1.; upwindCoeff(3,1)=0.
                  upwindCoeff(-3,2)=0.; upwindCoeff(-2,2)=-1.; upwindCoeff(-1,2)= 5.; upwindCoeff(0,2)=-10.; upwindCoeff(1,2)=10.; upwindCoeff(2,2)=-5.; upwindCoeff(3,2)=1.
                  upwindCoeff(-3,3)=0.; upwindCoeff(-2,3)=-1.; upwindCoeff(-1,3)= 4.; upwindCoeff(0,3)= -6.; upwindCoeff(1,3)= 4.; upwindCoeff(2,3)=-1.; upwindCoeff(3,3)=0.
              else
                  write(*,'("advWave: upwind coefficients are not set yet for orderOfAccuracy=",i4)') orderOfAccuracy
                  stop 2222
              end if
            end if
          if( (.false. .or. debug.gt.1) .and. t.le.3*dt )then
              write(*,'("advWave: option=",i4," grid=",i4)') option,grid
              write(*,'("advWave: orderOfAccuracy=",i2," orderInTime=",i2  )') orderOfAccuracy,orderInTime
              write(*,'("advWave: addForcing=",i2," forcingOption=",i2)') addForcing,forcingOption
              write(*,'("advWave: useUpwindDissipation=",i2," (explicit), useImplicitUpwindDissipation=",i2," (implicit)")') useUpwindDissipation,useImplicitUpwindDissipation
              write(*,'("advWave: useSosupDissipation=",i2," (1= add upwind dissipation in this stage)")') useSosupDissipation
              write(*,'("advWave: t,dt,c,omega,gridCFL=",5(1pe10.2,1x))') t,dt,cc,omega,gridCFL
              write(*,'("advWave: gridIsImplicit=",i2," adjustOmega=",i2," solveHelmholtz=",i2)') gridIsImplicit,adjustOmega,solveHelmholtz
              if( forcingOption.eq.helmholtzForcing )then
                  write(*,'("advWave: numberOfFrequencies=",i2)') numberOfFrequencies
                  write(*,'("advWave: frequencyArray=",(1pe12.4,1x))') (frequencyArray(freq),freq=0,numberOfFrequencies-1)
              end if
              if( gridIsImplicit.eq.1 )then
                  write(*,'("  Implicit coeff: cImp(-1:1,0) = ",3(1pe10.2,1x), "(for 2nd-order)")') cImp(-1,0),cImp(0,0),cImp(1,0)
                  write(*,'("  Implicit coeff: cImp(-1:1,1) = ",3(1pe10.2,1x), "(for 4th-order)")') cImp(-1,1),cImp(0,1),cImp(1,1)
              end if
          end if
          if( forcingOption.eq.helmholtzForcing )then
       ! --- solving the Helmholtz problem ---
              if( t.le.dt .and. debug.gt.1 )then
                  write(*,'("advWave: numberOfFrequencies=",i6," omega=",1pe12.4," frequencyArray(0)=",1pe12.4)') numberOfFrequencies,omega,frequencyArray(0)
              end if
              if( numberOfFrequencies.le.0 )then
                  write(*,'("advWave: ERROR: numberOfFrequencies=",i6," is <= 0")') numberOfFrequencies
                  stop 0123
              end if
              if( numberOfFrequencies.eq.1  .and. frequencyArray(0) .ne. omega )then
                  write(*,'("advWave: ERROR: frequencyArray(0)=",1pe12.4," is not equal to omega=",1pe12.4)') frequencyArray(0),omega
                  stop 1234
              end if
              if( numberOfFrequencies.gt.maxFreq )then
                  write(*,'("advWave: ERROR: numberOfFrequencies > maxFreq=",i6," .. FIX ME")') maxFreq
                  stop 2345
              end if
       ! if( numberOfFrequencies.gt.1 .and. gridIsImplicit.eq.1 )then
       !   write(*,'("advWave: ERROR: numberOfFrequencies > 1 and implicit time-stepping : FINISH ME")') 
       !   stop 3456  
       ! end if
              do freq=0,numberOfFrequencies-1
                  cosFreqt(freq) = cos(frequencyArray(freq)*t)
              end do
          end if
          if( useSosupDissipation.ne.0 .or. useImplicitUpwindDissipation.eq.1 )then
       ! ---- coefficients of upwind dissipation ---
       !   **NOTE**: These must match the values in implicit.bC 
       ! if( .true. )then
       ! *new* way to define coefficients: (see cgWave.pdf)
              adSosup = cc*dt/( sqrt(1.*nd) * 2**(orderOfAccuracy+1) ) 
       ! if( orderOfAccuracy.eq.6 )then
       !   adSosup = adSosup*.5     ! *** TEMP : ME66 seems to be unstable at cfl=.9 and original adSosup
       ! end if
         ! upwinding always takes a positive coefficient now *wdh* Feb 5, 2022
         ! if( mod(orderOfAccuracy/2,2).eq.1 )then
         !   adSosup = - adSosup  ! negative for orders 2,6,10, ...
         ! end if
       ! else 
       !   ! **** OLD WAY ****
       !   ! Coefficients in the sosup dissipation from Jordan Angel
       !   if( orderOfAccuracy.eq.2 )then
       !    adSosup=-cc*dt*1./8.
       !    if( preComputeUpwindUt.eq.1 )then
       !      ! We need to reduce the upwind coefficient for stability if we pre-compute uDot: (see CgWave documentation)
       !      adSosup=adSosup/sqrt(1.*nd)
       !    end if 
       !   else if( orderOfAccuracy.eq.4 )then 
       !     adSosup=cc*dt*5./288.
       !     if( .false. )then 
       !        adSosup = adSosup*.5 ! ****TEST****
       !     end if
       !   else if( orderOfAccuracy.eq.6 )then 
       !     adSosup=-cc*dt*31./8640.
       !   else
       !     stop 1005
       !   end if
       ! end if
              uDotFactor=.5  ! By default uDot is D-zero and so we scale (un-um) by .5 --> .5*(un-um)/(dt)
       ! sosupParameter=gamma in sosup scheme  0<= gamma <=1   0=centered scheme
              adSosup=sosupParameter*adSosup
              if( gridIsImplicit.ne.0 )then
                  write(*,'("advWave: gridIsImplicit: REDUCE UPWIND DISS COEFF by gridCFL=",e10.2)') gridCFL
                  adSosup = adSosup/gridCFL 
              end if
              if( (.false. .or. debug.gt.1) .and. t.le.2*dt )then
                  write(*,'("advMxWave: grid=",i3," gridType=",i2," orderOfAccuracy=",i2," useImplicitUpwindDissipation=",i2)') grid,gridType,orderOfAccuracy,useImplicitUpwindDissipation
                  write(*,'("         : t,dt,adSosup=",3e10.2," adSosup/(c*dt)=",e12.4)')t,dt,adSosup,adSosup/(cc*dt)
                  write(*,'("         : useSosupDissipation=",i2," sosupParameter=",1pe10.2," preComputeUpwindUt=",i2)') useSosupDissipation,sosupParameter,preComputeUpwindUt
         ! write(*,'("advMxUp: updateDissipation=",i2)') updateDissipation
         ! write(*,'("advMxUp: useNewForcingMethod=",i2)') useNewForcingMethod
              end if
       ! ! Coefficients of the sosup dissipation with Cartesian grids:
       ! cdSosupx= adSosup/dx(0)
       ! cdSosupy= adSosup/dx(1)
       ! cdSosupz= adSosup/dx(2)
       ! Note: these next values are only used for rectangular grids. (curvilinear grid values are computed in the loops)
              adxSosup(0)=  uDotFactor*adSosup/dx(0)
              adxSosup(1)=  uDotFactor*adSosup/dx(1)
              adxSosup(2)=  uDotFactor*adSosup/dx(2)
          end if
          if( useSosupDissipation.eq.1 .and. preComputeUpwindUt.eq.1 )then
              if( option.ne.1 )then
                  write(*,'("advWaveOpt:ERROR: useSosupDissipation.eq.1 BUT option.ne.1")')
                  stop 6663
              end if
         ! precompute "uDot" = dt*du/dt used in the dissipation and store in v 
         ! we need uDot at enough ghost points for the dissipation operator 
                  if( debug.gt.3 .and. t.le.3.*dt )then
                    write(*,'(" advWave: add UPWIND DISSIPATION: Evaluate v= uDot")') 
                  end if
                  numGhost=orderOfAccuracy/2
                  numGhost=numGhost+1
                  m1a=n1a-numGhost
                  m1b=n1b+numGhost
                  m2a=n2a-numGhost
                  m2b=n2b+numGhost
                  if( nd.eq.2 )then
                  m3a=n3a
                  m3b=n3b
                  else
                    m3a=n3a-numGhost
                    m3b=n3b+numGhost
                  end if
         ! write(*,'(" numGhost=",i2," m1a,m1b,m2a,m2b=",4i4)') numGhost,m1a,m1b,m2a,m2b
         ! We need v at ghost outside interpolation points -- do not use mask here
                  if( (adjustHelmholtzForUpwinding.eq.0) .or. adjustOmega.eq.0 .or. solveHelmholtz.eq.0 )then 
                          do i3=m3a,m3b
                          do i2=m2a,m2b
                          do i1=m1a,m1b
                          v(i1,i2,i3,0)=un(i1,i2,i3,0)-um(i1,i2,i3,0)
                          end do
                          end do
                          end do
                  else 
           ! ***THIS OPTION IS UNDER DEVELOPMENT***
           ! subtract off upwinding applied to current Helmholtz solution to cancel the effect of upwinding.
           ! Periodic solution:
           !    up = vk * cos( omega * t )
           !   D0t( up  ) = vk * D0t( cos( omega * t ) )
           ! NOTE: upwinding is called at the start of the step to previous times (t) (t-dt) and (t-2*dt)
           ! write(*,'("advOpt: adjust upwinding: t=",e10.3," freq=",e14.6)') t,frequencyArray(0)
           ! cosineFactor = cos(frequencyArray(0)*(t+0.*dt)) - cos(frequencyArray(0)*(t-2.*dt))
           ! Compute some things needed in the loops below
                      if( adjustHelmholtzForUpwinding.eq.1 )then
                          do freq=0,numberOfFrequencies-1
               ! NOTE: upwinding is called at the start of the step to previous times (t) (t-dt) and (t-2*dt)
                              cosineFactor(freq) = cos(frequencyArray(freq)*(t+0.*dt)) - cos(frequencyArray(freq)*(t-2.*dt))
                          end do
                      end if
                      if( .true. )then
             ! ----- Adjustment for upwinding in Helmholtz solves ----
                              do i3=m3a,m3b
                              do i2=m2a,m2b
                              do i1=m1a,m1b
                              v(i1,i2,i3,0)= (un(i1,i2,i3,0)-um(i1,i2,i3,0)) 
                              do freq=0,numberOfFrequencies-1
                                  v(i1,i2,i3,0) = v(i1,i2,i3,0) - vh(i1,i2,i3,freq)*cosineFactor(freq)
                              end do
                              end do
                              end do
                              end do
                      else if( .false. )then
             ! **TESTING : 
             !   CHECK: u = vk * cos( omega * t )
             ! BUT DO NOT ADD UPWINDING 
                          maxDiff=0.
             ! beginLoops(i1,i2,i3,m1a,m1b,m2a,m2b,m3a,m3b)
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                              v(i1,i2,i3,0)= um(i1,i2,i3,0) - vh(i1,i2,i3,0)*cos(frequencyArray(0)*(t-2.*dt))
               ! write(*,'("(i1,i2)=",2i4," um=",e12.4," vh*cos=",e12.4," diff=",e8.2)') i1,i2,um(i1,i2,i3,0), vh(i1,i2,i3,0)*cos(frequencyArray(0)*(t-2.*dt)),v(i1,i2,i3,0)
                              maxDiff = max(maxDiff,v(i1,i2,i3,0));
                              end do
                              end do
                              end do
                          write(*,'("advOpt: adjust upwinding: check: maxDiff(um-vh*cos)=",e10.3," (no ghost)")') maxDiff
                          maxDiff=0.
                              do i3=m3a,m3b
                              do i2=m2a,m2b
                              do i1=m1a,m1b
                              v(i1,i2,i3,0)= um(i1,i2,i3,0) - vh(i1,i2,i3,0)*cos(frequencyArray(0)*(t-2.*dt))
               ! write(*,'("(i1,i2)=",2i4," um=",e12.4," vh*cos=",e12.4," diff=",e8.2)') i1,i2,um(i1,i2,i3,0), vh(i1,i2,i3,0)*cos(frequencyArray(0)*(t-2.*dt)),v(i1,i2,i3,0)
                              maxDiff = max(maxDiff,v(i1,i2,i3,0));
                              end do
                              end do
                              end do
                          write(*,'("advOpt: adjust upwinding: check: maxDiff(um-vh*cos)=",e10.3," (with ghost)")') maxDiff      
                          maxDiff=0.
             ! beginLoops(i1,i2,i3,m1a,m1b,m2a,m2b,m3a,m3b)
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                              v(i1,i2,i3,0)= un(i1,i2,i3,0) - vh(i1,i2,i3,0)*cos(frequencyArray(0)*(t+0.*dt))
                              maxDiff = max(maxDiff,v(i1,i2,i3,0));
                              end do
                              end do
                              end do
                          write(*,'("advOpt: adjust upwinding: check: maxDiff(un-vh*cos)=",e10.3," (no ghost)")') maxDiff 
                          maxDiff=0.
                              do i3=m3a,m3b
                              do i2=m2a,m2b
                              do i1=m1a,m1b
                              v(i1,i2,i3,0)= un(i1,i2,i3,0) - vh(i1,i2,i3,0)*cos(frequencyArray(0)*(t+0.*dt))
                              maxDiff = max(maxDiff,v(i1,i2,i3,0));
                              end do
                              end do
                              end do
                          write(*,'("advOpt: adjust upwinding: check: maxDiff(un-vh*cos)=",e10.3, " (with ghost)")') maxDiff            
                              do i3=m3a,m3b
                              do i2=m2a,m2b
                              do i1=m1a,m1b
               ! v(i1,i2,i3,0)= (un(i1,i2,i3,0)-um(i1,i2,i3,0)) - vh(i1,i2,i3,0)*cosineFactor
                              v(i1,i2,i3,0)= 0.
                              end do
                              end do
                              end do
                      end if  
                  end if
          end if 
     ! write(*,'(" advWave: timeSteppingMethod=",i2)') timeSteppingMethod
          if( timeSteppingMethod.eq.defaultTimeStepping )then
            write(*,'(" advWave:ERROR: timeSteppingMethod=defaultTimeStepping -- this should be set")')
        ! '
            stop 83322
          end if
          if( useSosupDissipation.eq.0 )then
              if( gridIsImplicit.eq.0 )then 
         ! ------- EXPLIICT update the solution ---------
                          if( (orderOfAccuracy.eq.6 .or. debug.gt.3) .and. t.lt.2*dt )then
                              write(*,'("advWave: ADVANCE dim=2 order=2 orderInTime=2, grid=curvilinear... t=",e10.2)') t
                          end if
             ! --- TAYLOR TIME-STEPPING --- 
                          m=0 ! component number 
                          ec = 0 ! component number
                          if( forcingOption.eq.helmholtzForcing )then
                              coswt = cos(omega*t)
                          end if 
                          fv(m)=0.
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                  if( mask(i1,i2,i3).gt.0 )then
                                  if( forcingOption.eq.twilightZoneForcing )then
                                      if( nd.eq.2 )then
                                              call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,ev(m) )
                                              call ogDeriv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evtt(m) )
                                              call ogDeriv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evxx(m) )
                                              call ogDeriv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evyy(m) )
                                          fv(m) = evtt(m) - csq*( evxx(m) + evyy(m) )
                                      else
                                              call ogDeriv(ep, 0,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,ev(m) )
                                              call ogDeriv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evtt(m) )
                                              call ogDeriv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evxx(m) )
                                              call ogDeriv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evyy(m) )
                                              call ogDeriv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evzz(m) )
                                          fv(m) = evtt(m) - csq*( evxx(m) + evyy(m)  + evzz(m) )
                                    end if
                                else if( forcingOption.eq.helmholtzForcing )then
                   ! forcing for solving the Helmholtz equation   
                   ! NOTE: change sign of forcing since for Helholtz we want to solve
                   !      ( omega^2 I + c^2 Delta) w = f 
                   ! fv(m) = -f(i1,i2,i3,0)*coswt  
                                      fv(m)=0.
                                      do freq=0,numberOfFrequencies-1 
                                          omega = frequencyArray(freq)
                                          coswt = cosFreqt(freq)    
                      ! if( i1.eq.2 .and. i2.eq.2 )then 
                      !   write(*,'(" adv: forcing f(i1,i2,i3)=",1pe12.4," coswt=",1pe12.4," t=",1pe12.4," omega=",1pe12.4)') f(i1,i2,i3,0),coswt,t,omega
                      ! end if
                      ! fv(m) = -f(i1,i2,i3,0)*coswt  
                                            fv(m) = fv(m) - f(i1,i2,i3,freq)*coswt
                                      end do ! do freq  
                                else if( addForcing.ne.0 )then  
                                      fv(m) = f(i1,i2,i3,0)
                                end if
                ! --- SECOND 2 ---
                  ! --- TWO DIMENSIONS ---
                                        un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( uxx22(i1,i2,i3,0)  + uyy22(i1,i2,i3,0) ) + dtSq*fv(m)
              ! write(*,'("i1,i2=",2i3," u-ue=",e10.2)') i1,i2,u(i1,i2,i3,m)-ev(m)
              ! write(*,'(" uxx-uxxe =",e10.2)') uxx22r(i1,i2,i3,0)-evxx(m)
              ! OGDERIV2D( 0,0,0,0,i1,i2,i3,t+dt, ec, ev(m)  )
              ! write(*,'(" un-ue=",e10.2)') un(i1,i2,i3,m)-ev(m)
                                  end if
                              end do
                              end do
                              end do
              else
         ! --- IMPLICIT: Fill in RHS to implicit time-stepping -----
                          if( (.true. .or. debug.gt.3) .and. t.lt.2*dt )then
                              write(*,'("advWave: ADVANCE IMPLICIT dim=2 order=2 orderInTime=2, grid=curvilinear... t=",e10.2)') t
                          end if
             ! --- IMPLICIT TAYLOR TIME-STEPPING --- 
             !  D+t D-t u = c^2 Delta( cImp(1) *u^{n+1} + cImp(0) *u^n + cImp(-1)* u^{n-1} )
                          m=0 ! component number 
                          ec = 0 ! component number
             ! #If "curvilinear" eq "curvilinear"
             !   #If "2" eq "4" && "2" eq "4"
             !     computeLaplacianOrder2(2)
             !   #End
             ! #End
                          if( forcingOption.eq.helmholtzForcing )then
                              do freq=0,numberOfFrequencies-1
                                  coswt = cos(frequencyArray(freq)*t) ! is this used?
                                  omega = frequencyArray(freq)
                                  coswtAve(freq) = cImp(-1,0)*cos(omega*(t-dt)) + cImp(0,0)*cos(omega*t) + cImp(1,0)*cos(omega*(t+dt))
                              end do
               ! coswt = cos(omega*t)
               ! coswtAve = cImp(-1,0)*cos(omega*(t-dt)) + cImp(0,0)*cos(omega*t) + cImp(1,0)*cos(omega*(t+dt))    
                          end if 
                          fv(m)=0.
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                  if( mask(i1,i2,i3).gt.0 )then
                                  if( forcingOption.eq.twilightZoneForcing )then
                                      if( nd.eq.2 )then
                     ! OGDERIV2D( 0,0,0,0,i1,i2,i3,t, ec, ev(m)  )
                                              call ogDeriv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,t, ec,evtt(m) )
                                          fv(m) = evtt(m)
                     ! We weight the TZ forcing with the implicit weights to make the solution exact for polynomials
                                          do mt=-1,1
                                              tm = t + dt*mt
                                                  call ogDeriv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,tm, ec,evxx(m) )
                                                  call ogDeriv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),0.,tm, ec,evyy(m) )
                                              fv(m) = fv(m) -csq*( cImp(mt,0)*( evxx(m) + evyy(m) )  )
                                          end do
                                      else
                     ! OGDERIV3D( 0,0,0,0,i1,i2,i3,t, ec, ev(m)  )
                                              call ogDeriv(ep, 2,0,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),t, ec,evtt(m) )
                                          fv(m) = evtt(m)
                     ! We weight the TZ forcing with the implicit weights to make the solution exact for polynomials
                                          do mt=-1,1
                                              tm = t + dt*mt       
                                                  call ogDeriv(ep, 0,2,0,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),tm, ec,evxx(m) )
                                                  call ogDeriv(ep, 0,0,2,0, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),tm, ec,evyy(m) )
                                                  call ogDeriv(ep, 0,0,0,2, xy(i1,i2,i3,0),xy(i1,i2,i3,1),xy(i1,i2,i3,2),tm, ec,evzz(m) )
                                              fv(m) = fv(m) -csq*( cImp(mt,0)*( evxx(m) + evyy(m) + evzz(m) )  ) 
                                          end do
                                    end if
                                else if( forcingOption.eq.helmholtzForcing )then
                  ! forcing for solving the Helmholtz equation   
                  ! NOTE: change sign of forcing since for Helholtz we want to solve
                  !      ( omega^2 I + c^2 Delta) w = f    
                  ! Define
                  !   coswtAve = cImp(-1)*cos(omega*(t-dt)) + cImp(0)*cos(omega*t) + cImp(1)*cos(omega*(t+dt))
                                      fv(m)=0. 
                                      do freq=0,numberOfFrequencies-1
                                          omega = frequencyArray(freq)
                                              fv(m) = fv(m) -f(i1,i2,i3,freq)*coswtAve(freq)    
                                      end do
                                else if( addForcing.ne.0 )then  
                                      fv(m) = f(i1,i2,i3,0)
                                end if
                ! --- SECOND 2 ---
                  ! --- TWO DIMENSIONS ---
                                        un(i1,i2,i3,m)= 2.*u(i1,i2,i3,m) - um(i1,i2,i3,m) + (cdtSq)*( cImp( 0,0) *(  uxx22(i1,i2,i3,0) +  uyy22(i1,i2,i3,0) ) +  cImp(-1,0) *( umxx22(i1,i2,i3,0) + umyy22(i1,i2,i3,0) ) )  + dtSq*fv(m)
              ! write(*,'("i1,i2=",2i3," u-ue=",e10.2)') i1,i2,u(i1,i2,i3,m)-ev(m)
              ! write(*,'(" uxx-uxxe =",e10.2)') uxx22r(i1,i2,i3,0)-evxx(m)
              ! OGDERIV2D( 0,0,0,0,i1,i2,i3,t+dt, ec, ev(m)  )
              ! write(*,'(" un-ue=",e10.2)') un(i1,i2,i3,m)-ev(m)
                                  end if
                              end do
                              end do
                              end do
         ! --- Add contributions from upwind dissipation ---
                  if( useImplicitUpwindDissipation.eq.1 )then
                          if( (.true. .or. debug.gt.3) .and. t.lt.2*dt )then
                              write(*,'("addUpwindDissImplicit: UPWIND DISS dim=2 order=2 grid=curvilinear... t=",e10.2)') t
                              write(*,'(" adxSosup=",3e12.4)') adxSosup(0), adxSosup(1),adxSosup(2)
                          end if
             ! from formImplicitTimeSteppingMatrix
               ! // fourth-order dissipation for 2nd-order scheme:
               ! Real upwindCoeff4[4][5] = { 1.,-4.,6.,-4.,1.,
               !                             1.,-3.,3.,-1.,0.,   // extrap right-most point D-^3 u(2)
               !                             0.,-1.,3.,-3.,1.,   // extrap left -most point D+^3 u(-2)
               !                             0.,-1.,2.,-1.,0.
               !                           };
               ! // sixth-order dissipation for 4th-order scheme
               ! Real upwindCoeff6[4][7] = {1.,-6.,15.,-20.,15.,-6.,1.,
               !                            1.,-5.,10.,-10., 5.,-1.,0.,  // extrap right-most point D-^5 u(3)
               !                            0.,-1., 5.,-10.,10., 5.,1.,  // extrap left -most point D+^5 u(-3)
               !                            0.,-1., 4., -6., 4.,-1.,0.
               !                           };
             ! -- Note: Could adjust loop bounds to avoid Dirichlet boundaries
                          m=0 ! component number 
                          ec = 0 ! component number
             ! Compute some things needed in the loops below
                          if( adjustHelmholtzForUpwinding.eq.1 )then
                              do freq=0,numberOfFrequencies-1
                                  cosineFactor(freq) = cos(frequencyArray(freq)*(t+dt)) - cos(frequencyArray(freq)*(t-dt))
                              end do
                          end if 
             ! stencil width = order + 1
             ! upwind stencil = stencilWidth + 2 = order+1 + 2
                          upwindHalfStencilWidth = (orderOfAccuracy+2)/2 
                              do i3=n3a,n3b
                              do i2=n2a,n2b
                              do i1=n1a,n1b
                                  if( mask(i1,i2,i3).gt.0 )then
                 ! --- Compute UPW coefficient for curvilinear grids ---
                                        do dir=0,1
                      ! diss-coeff ~= 1/(change in x along direction r(dir) )
                      ! Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
                                            adxSosup(dir) = adSosup*uDotFactor*sqrt( rsxy(i1,i2,i3,dir,0)**2 + rsxy(i1,i2,i3,dir,1)**2 )/dr(dir) 
                                        end do
               ! --- loop over directions ---
                              do dir=0,nd-1
                                  idv(0)=0; idv(1)=0; idv(2)=0
                                  idv(dir)=1 ! active direction
                 ! check if left-most and right-most entries in the upwind stencil are valid 
                                  i1l = i1-upwindHalfStencilWidth*idv(0); i1r = i1+upwindHalfStencilWidth*idv(0);
                                  i2l = i2-upwindHalfStencilWidth*idv(1); i2r = i2+upwindHalfStencilWidth*idv(1);
                                  i3l = i3-upwindHalfStencilWidth*idv(2); i3r = i3+upwindHalfStencilWidth*idv(2);
                 !  Note: there are at most four cases at any order, since we have order/2 layers of interpolation points 
                 !   Example, order=2, upwind-order=4
                 !     X---X---C---X---X           C = center point = valid discretization point
                 !     X---X---C---X               missing right-most 
                 !         X---C---X---X           missing left-most
                 !         X---C---X               missing left and right-most
                                  upwCase=0 
                                  if( mask(i1l,i2l,i3l) .ne. 0 .and. mask(i1r,i2r,i3r) .ne. 0 ) then
                                      upwCase=0  ! centred, full-width stencil
                                  else if( mask(i1l,i2l,i3l) .ne. 0 ) then
                                      upwCase=1  ! left biased stencil
                                  else if( mask(i1r,i2r,i3r) .ne. 0 ) then
                                      upwCase=2  ! right biased stencil   
                                  else  
                                      upwCase=3  ! centred smaller stencil     
                                  end if
                                  upw = 0. 
                                  do iStencil=-upwindHalfStencilWidth,upwindHalfStencilWidth
                                      j1 = i1 + iStencil*idv(0);  j2 = i2 + iStencil*idv(1);  j3 = i3 + iStencil*idv(2)
                                      umj = um(j1,j2,j3,0)
                                      if( adjustHelmholtzForUpwinding.eq.1 )then
                                          do freq=0,numberOfFrequencies-1
                                              umj = umj + cosineFactor(freq)*vh(j1,j2,j3,freq)
                                          end do
                                      end if
                                      upw = upw + upwindCoeff(iStencil,upwCase)*umj
                   ! upw = upw + upwindCoeff(iStencil,upwCase)*um(j1,j2,j3,ec)  ! *** CHECK ME 
                   ! write(*,'("upw-rhs: i1,i2=",2i4," j1,j2=",2i4," upwindCoeff=",1pe9.2, " um=",1pe9.2," upw=",1pe9.2)') i1,i2,j1,j2,upwindCoeff(iStencil,upwCase),um(j1,j2,j3,ec),upw
                                  end do 
                 ! if( abs(upw).gt.1e-10 )then
                 !   write(*,'(">>upw-rhs: i1,i2=",2i4," upw=",1pe9.2)') i1,i2,upw
                 ! end if 
                 ! This is the coeff of um in
                 !         + adxSosup(dir)*(UpwindStencil)( un - um )
                                  un(i1,i2,i3,ec) = un(i1,i2,i3,ec) - adxSosup(dir)*upw 
                 ! TEST un(i1,i2,i3,ec) = un(i1,i2,i3,ec) - adxSosup(dir)*upw 
                              end do
                                  end if
                              end do
                              end do
                              end do
                  end if 
              end if
          else
       ! ---- add upwind dissipation -----
       ! preComputeUpwindUt : true=precompute Ut in upwind dissipation,  (uses v=uDot computed above)
       !                      false=compute Ut inline in Gauss-Seidel fashion 
              if( preComputeUpwindUt.eq.1 )then
         ! precompute Ut in upwind dissipation,  (uses v=uDot computed above)
                      if( debug.gt.3 .and. t.lt.4*dt )then
                          write(*,'("addUpwindDiss: UPWIND DISS using u.t=v dim=2 order=2 grid=curvilinear... t=",e10.2)') t
                          write(*,'(" adxSosup=",3e12.4)') adxSosup(0), adxSosup(1),adxSosup(2)
                      end if
                      m=0 ! component number 
                      ec = 0 ! component number
                          do i3=n3a,n3b
                          do i2=n2a,n2b
                          do i1=n1a,n1b
                              if( mask(i1,i2,i3).gt.0 )then
              ! --- SECOND 2 ---
                ! --- TWO DIMENSIONS ---
                                      do dir=0,1
                     ! diss-coeff ~= 1/(change in x along direction r(dir) )
                     ! Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
                                          adxSosup(dir) = adSosup*uDotFactor*sqrt( rsxy(i1,i2,i3,dir,0)**2 + rsxy(i1,i2,i3,dir,1)**2 )/dr(dir) 
                                      end do
                                    un(i1,i2,i3,ec)=un(i1,i2,i3,ec)+(-6.*v(i1,i2,i3,ec)+4.*(v(i1+1,i2,i3,ec)+v(i1-1,i2,i3,ec))-(v(i1+2,i2,i3,ec)+v(i1-2,i2,i3,ec)))*adxSosup(0)+(-6.*v(i1,i2,i3,ec)+4.*(v(i1,i2+1,i3,ec)+v(i1,i2-1,i3,ec))-(v(i1,i2+2,i3,ec)+v(i1,i2-2,i3,ec)))*adxSosup(1)
                              end if
                          end do
                          end do
                          end do
              else
         ! compute Ut inline in Gauss-Seidel fashion (this is more stable)
                      if( debug.gt.3 .and. t.lt.4*dt )then
                          write(*,'("addUpwindDiss: UPWIND DISS using u.t=Dztu dim=2 order=2 grid=curvilinear... t=",e10.2)') t
                          write(*,'(" adxSosup=",3e12.4)') adxSosup(0), adxSosup(1),adxSosup(2)
                      end if
                      m=0 ! component number 
                      ec = 0 ! component number
                          do i3=n3a,n3b
                          do i2=n2a,n2b
                          do i1=n1a,n1b
                              if( mask(i1,i2,i3).gt.0 )then
              ! --- SECOND 2 ---
                ! --- TWO DIMENSIONS ---
                                      do dir=0,1
                     ! diss-coeff ~= 1/(change in x along direction r(dir) )
                     ! Assuming a nearly orthogonal grid gives ||dx|| = || grad_x(r_i) || / dr_i 
                                          adxSosup(dir) = adSosup*uDotFactor*sqrt( rsxy(i1,i2,i3,dir,0)**2 + rsxy(i1,i2,i3,dir,1)**2 )/dr(dir) 
                                      end do
                                    un(i1,i2,i3,ec)=un(i1,i2,i3,ec)+(-6.*Dztu(i1,i2,i3,ec)+4.*(Dztu(i1+1,i2,i3,ec)+Dztu(i1-1,i2,i3,ec))-(Dztu(i1+2,i2,i3,ec)+Dztu(i1-2,i2,i3,ec)))*adxSosup(0)+(-6.*Dztu(i1,i2,i3,ec)+4.*(Dztu(i1,i2+1,i3,ec)+Dztu(i1,i2-1,i3,ec))-(Dztu(i1,i2+2,i3,ec)+Dztu(i1,i2-2,i3,ec)))*adxSosup(1)
                              end if
                          end do
                          end do
                          end do
              end if
          end if
          return
          end
