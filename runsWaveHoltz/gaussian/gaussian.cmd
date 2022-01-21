#
#  cgWaveHoltz command script: 
# 
#     cgwh waveHoltz.cmd -g=<grid-name> -x0=<f> -y0=<f> -omega=<f> -solver=[none|fixedPoint|krylov] -tol=<f> -tp=<f> ...
#                         -kx=<f> -ky=<f> -kz=<f> -forcing=[gaussian|sine] -adjustOmega=[0|1] -maxIterations=<>
#                          -upwind=[0|1] -imode=[0|1] -bcApproach=[cbc|lcbc|oneSided] -go=[go|og|halt]
#
#   -solver=[fixedPoint|krylov] : fiexed-point or Krylov 
#   -imode=1 : do not wait in cgWave
#   -maxIterations=<>
#
$go="go"; $forcing="gaussian"; 
$omega=30.1; 
# $beta=50.; # $x0=0.5; $y0=0.5; $z0=0.5; 
@beta= (); @amp = (); @x0 = (); @y0 = (); @z0 = (); # this must be null for GetOptions to work, defaults are given below
$t0=0.; 
# $amp=1.; 
$numPeriods=1; $omegaSOR=1; $tol=1.e-3; 
$numberOfFrequencies=1; 
@freq = ();  # this must be null for GetOptions to work, defaults are given below
# beta2=0 = TRAP, beta2=.5 = Full weighting 
$beta2=.0; $beta4=0.; $beta6=0.; $beta8=0.; # weights in implicit time-stepping 
$upwind=0; # new way
$tp=.5; $imode=0; 
$solver="fixedPoint";  $kx=1; $ky=1; $kz=1; $maxIterations=100; $adjustOmega=0; 
$matlab="cgWaveHoltz"; $show="gaussian.show"; 
$cfl=.9; $bc="d"; $ts="explicit"; $dtMax=1; 
$bcApproach="oneSided"; # bc Approach : cbc, lcbc, oneSided
$orderInTime=-1;  # -1 = use default
GetOptions( "omega=f"=>\$omega,"x0=f{1,}"=>\@x0,"y0=f{1,}"=>\@y0,"z0=f{1,}"=>@$z0,"beta=f{1,}"=>\@beta,"numPeriods=i"=>\$numPeriods,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"cfl=f"=>\$cfl,"tp=f"=>\$tp,"iMode=i"=>\$imode,\
            "solver=s"=>\$solver,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"maxIterations=i"=>\$maxIterations,"matlab=s"=>\$matlab,\
            "go=s"=>\$go,"forcing=s"=>\$forcing,"bc=s"=>\$bc,"ts=s"=>\$ts,"orderInTime=i"=>\$orderInTime,\
            "dtMax=f"=>\$dtMax,"adjustOmega=i"=>\$adjustOmega,"amp=f{1,}"=>\@amp,"show=s"=>\$show,\
            "bcApproach=s"=>\$bcApproach,"upwind=i"=>\$upwind,"beta2=f"=>\$beta2,"beta4=f"=>\$beta4,"beta6=f"=>\$beta6,\
            "numberOfFrequencies=i"=>\$numberOfFrequencies,"nf=i"=>\$numberOfFrequencies,"freq=f{1,}"=>\@freq );
# 
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
if( $amp[0] eq "" ){ @amp=(1,1,1,1,1,1,1,1,1); }
if( $beta[0] eq "" ){ @beta=(50,50,50,50,50,50,50,50); }
if( $x0[0] eq ""  ){ @x0=(0,0,0,0,0,0,0,0); }
if( $y0[0] eq ""  ){ @y0=(0,0,0,0,0,0,0,0); }
if( $z0[0] eq ""  ){ @z0=(0,0,0,0,0,0,0,0); }
if( $freq[0] eq "" ){ @freq=(1,2,3,4,5,6,7,8,9,10); }
# 
$omega = $freq[0];
# 
# pause
## Gaussian params $beta $x0 $y0 0 (beta,x0,y0,z0)
omega $omega
maximum number of iterations $maxIterations 
tol $tol 
number of periods $numPeriods
adjust omega $adjustOmega
number of frequencies $numberOfFrequencies
frequencies $freq[0] $freq[1] $freq[2] $freq[3] $freq[4] $freq[5] $freq[6] $freq[7] $freq[8]
matlab filename: $matlab 
# pause
exit
# ------ Start cgWave setup ------
# time-stepping : explicit/implicit
$ts
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
$cmd
dtMax $dtMax
implicit weights $beta2 $beta4 $beta6 $beta8
cfl $cfl
interactiveMode $imode 
tPlot $tp 
# -- Here is input for cgWave 
#
$cmd="#";
if( $bcApproach eq "oneSided" ){ $cmd="useOneSidedBCs"; }
if( $bcApproach eq "cbc"      ){ $cmd="useCompatibilityBCs"; }
if( $bcApproach eq "lcbc"     ){ $cmd="useLocalCompatibilityBCs"; }
$cmd
#
bc=$bc
#
# if( $ad4>0. ){ $upwind=1; }# for backward compatibility
upwind dissipation $upwind
# artificial dissipation $ad4
# 
helmholtzForcing
solve Helmholtz 1
#
# ------ Gaussian sources ----
# amp, beta, omega, p, 
# $amp = 10.*$omega*$omega; 
#
$numTerms=1; # each source has 1 term 
$p=1;        # power in exponent
$cmd =  "user defined forcing...\n" \
      . " gaussian sources\n";
for( $i=0; $i<$numberOfFrequencies; $i++ ){\
   $cmd .= "$numTerms\n $amp[$i] $beta[$i] $omega $p $x0[$i] $y0[$i] $z0[$i] $t0\n";\
}
$cmd .= "exit";
# 
#   --- eigenmode-like exact solution : sin(kx*x)*sin(ky*y)*cos(omega*t)
if( $forcing eq "sine" ){ \
 $cmd = " user defined known solution...\n" \
      . " box helmholtz\n" \
      . "  $omega $kx $ky $kz\n" \
      . " done\n" \
      . " user defined forcing...\n" \
      . "   box Helmholtz\n" \
      . " exit\n"; }
$cmd
#
exit
# --- end cgWave setup  ---
show file $show
contour
exit
if( $solver eq "fixedPoint" ){ $cmd="compute with fixed-point"; }elsif( $solver eq "krylov" ){ $cmd="compute with petsc"; }else{ $cmd="#" }; 
if( $go eq "go" && $cmd ne "#" ){ $cmd .= "\n exit"; }
$cmd


# For disk--
plot forcing
DISPLAY AXES:0 0
plot boundaries (toggle)
set view:0 -0.0848943 0.253021 0 1.04747 1 0 0 0 0.5 -0.866025 0 0.866025 0.5
hardcopy file name:0 diskHelmholtzGaussianForcing0.ps
hardcopy save:0
plot:f1
hardcopy file name:0 diskHelmholtzGaussianForcing1.ps
hardcopy save:0
exit