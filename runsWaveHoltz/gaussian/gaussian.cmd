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
$omega=30.1; $beta=50.; $x0=0.5; $y0=0.5; $z0=0.5; $t0=0.; $amp=1.; 
$numPeriods=1; $omegaSOR=1; $tol=1.e-3; 
$upwind=0; # new way
$tp=.5; $imode=0; 
$solver="fixedPoint";  $kx=1; $ky=1; $kz=1; $maxIterations=100; $adjustOmega=0; 
$matlab="cgWaveHoltz"; $show="gaussian.show"; 
$cfl=.9; $bc="d"; $ts="explicit"; $dtMax=1; 
$bcApproach="oneSided"; # bc Approach : cbc, lcbc, oneSided
$orderInTime=-1;  # -1 = use default
GetOptions( "omega=f"=>\$omega,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"numPeriods=i"=>\$numPeriods,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"cfl=f"=>\$cfl,"tp=f"=>\$tp,"iMode=i"=>\$imode,\
            "solver=s"=>\$solver,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"maxIterations=i"=>\$maxIterations,"matlab=s"=>\$matlab,\
            "go=s"=>\$go,"forcing=s"=>\$forcing,"bc=s"=>\$bc,"ts=s"=>\$ts,"orderInTime=i"=>\$orderInTime,\
            "dtMax=f"=>\$dtMax,"adjustOmega=i"=>\$adjustOmega,"amp=f"=>\$amp,"show=s"=>\$show,\
            "bcApproach=s"=>\$bcApproach,"upwind=i"=>\$upwind );
# 
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
# 
# pause
Gaussian params $beta $x0 $y0 0 (beta,x0,y0,z0)
omega $omega
maximum number of iterations $maxIterations 
tol $tol 
number of periods $numPeriods
adjust omega $adjustOmega
matlab filename: $matlab 
# pause
exit
# ------ Start cgWave setup ------
# time-stepping : explicit/implicit
$ts
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
$cmd
dtMax $dtMax
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
# ------ Gaussian source ----
# amp, beta, omega, p, 
# $amp = 10.*$omega*$omega; 
$p=1; 
if( $forcing eq "gaussian"){ \
 $cmd = "user defined forcing...\n" \
      . " gaussian sources\n" \
      . "   1\n" \
      . " $amp $beta $omega $p $x0 $y0 $z0 $t0\n" \
      . " exit"; }
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