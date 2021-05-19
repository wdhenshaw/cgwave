#
#  cgWave command script for source term 
# 
#     cgwave source.cmd -g=<grid-name> -x0=<f> -y0=<f> -omega=<f> -solver=[none|fixedPoint|krylov] -tol=<f> -tp=<f> ...
#                         -kx=<f> -ky=<f> -kz=<f> -forcing=[gaussian|sine] -adjustOmega=[0|1]
#                         -ad4=<f> -imode=[0|1] -go=[go|og|halt]
#
#   -solver=[fixedPoint|krylov] : fiexed-point or Krylov 
#   -imode=1 : do not wait in cgWave
#   -maxIt=<>
#
$go="go"; $forcing="gaussian"; 
$omega=30.1; $beta=50.; $x0=0.5; $y0=0.5; $z0=0.5; $t0=0.
$beta=400; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; $ad4=1; $tp=.5; $imode=0; 
$solver="fixedPoint";  $kx=1; $ky=1; $kz=1; $maxIt=500; $adjustOmega=0; 
$matlab="cgWaveHoltz";
$cfl=.7; $bc="d"; $ts="explicit"; $dtMax=1; 
$orderInTime=-1;  # -1 = use default
GetOptions( "omega=f"=>\$omega,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"numPeriods=i"=>\$numPeriods,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"ad4=f"=>\$ad4,"cfl=f"=>\$cfl,"tp=f"=>\$tp,"iMode=i"=>\$imode,\
            "solver=s"=>\$solver,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"maxIt=i"=>\$maxIt,"matlab=s"=>\$matlab,\
            "go=s"=>\$go,"forcing=s"=>\$forcing,"bc=s"=>\$bc,"ts=s"=>\$ts,"orderInTime=i"=>\$orderInTime,\
            "dtMax=f"=>\$dtMax,"adjustOmega=i"=>\$adjustOmega );
# 
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
# 
# ------ Start cgWave setup ------
# time-stepping : explicit/implicit
$ts
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
$cmd
dtMax $dtMax
interactiveMode $imode 
tPlot $tp 
# -- Here is input for cgWave 
bc=$bc
#
artificial dissipation $ad4
#
zero initial condition
# 
helmholtzForcing
# solve Helmholtz 1
#
# ------ Gaussian source ----
# amp, beta, omega, p, 
$amp = 10.*$omega*$omega; $p=1; 
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
# 
solve




# --- end cgWave setup  ---
# contour
# exit
if( $solver eq "fixedPoint" ){ $cmd="compute with fixed-point"; }elsif( $solver eq "krylov" ){ $cmd="compute with petsc"; }else{ $cmd="#" }; 
if( $go eq "go" && $cmd ne "#" ){ $cmd .= "\n exit"; }
$cmd