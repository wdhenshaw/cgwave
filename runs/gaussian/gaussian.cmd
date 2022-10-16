#
#  cgWave command script: 
# 
#     cgwave gaussian.cmd -g=<grid-name> -x0=<f> -y0=<f> -omega=<f> -tol=<f> -tp=<f> ...
#                         -kx=<f> -ky=<f> -kz=<f> -forcing=[gaussian|sine] -upwind=[0|1] -imode=[0|1] 
#                         -bcApproach=[cbc|lcbc|oneSided] -meApproach=[std,ha] -rectangular=[implicit|explicit] -go=[go|og|halt]
#
#   -imode=1 : do not wait in cgWave
#
$go="go"; $forcing="gaussian"; 
$omega=30.1; $beta=50.; $x0=0.5; $y0=0.5; $z0=0.5; $t0=0.; $amp=1.; $debug=0; 
$ad4=0;    # old way
$upwind=1; # new way
$tf=5.; $tp=.5; $imode=0; 
$kx=1; $ky=1; $kz=1; 
$bcApproach="oneSided"; # bc Approach : cbc, lcbc, oneSided
$meApproach="std"; # or "ha"
$matlab="cgWave"; $show="gaussian.show"; 
$cfl=.9; $bc="d"; $ts="explicit"; $dtMax=1; 
$rectangular="implicit"; # for ts=implicit, set rectangular=explicit to treat rectangular grids explicitly
$orderInTime=-1;  # -1 = use default
$chooseImplicitTimeStepFromCFL=1; 
GetOptions( "omega=f"=>\$omega,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"numPeriods=i"=>\$numPeriods,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"ad4=f"=>\$ad4,"cfl=f"=>\$cfl,"tp=f"=>\$tp,"tf=f"=>\$tf,"iMode=i"=>\$imode,\
            "solver=s"=>\$solver,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"maxIterations=i"=>\$maxIterations,"matlab=s"=>\$matlab,\
            "go=s"=>\$go,"forcing=s"=>\$forcing,"bc=s"=>\$bc,"ts=s"=>\$ts,"orderInTime=i"=>\$orderInTime,\
            "dtMax=f"=>\$dtMax,"adjustOmega=i"=>\$adjustOmega,"amp=f"=>\$amp,"show=s"=>\$show,"upwind=i"=>\$upwind,\
            "debug=i"=>\$debug,"bcApproach=s"=>\$bcApproach,"meApproach=s"=>\$meApproach,"rectangular=s"=>\$rectangular );
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
if( $rectangular eq "explicit" ){ $chooseImplicitTimeStepFromCFL=0; }
choose implicit dt from cfl $chooseImplicitTimeStepFromCFL 
cfl $cfl
interactiveMode $imode 
debug $debug
tPlot $tp 
tFinal $tf
# 
$cmd="#";
if( $bcApproach eq "oneSided" ){ $cmd="useOneSidedBCs"; }
if( $bcApproach eq "cbc"      ){ $cmd="useCompatibilityBCs"; }
if( $bcApproach eq "lcbc"     ){ $cmd="useLocalCompatibilityBCs"; }
$cmd
# -- 
bc=$bc
$cmd="#";
if( $meApproach eq "std" ){ $cmd="standard modified equation"; }
if( $meApproach eq "ha" ){ $cmd="hierarchical modified equation"; }
$cmd
#
if( $ts eq "implicit" ){ $cmd="choose grids for implicit\n  rectangular=$rectangular\n done"; }else{ $cmd="#"; }
$cmd
#
#
if( $ad4>0. ){ $upwind=1; }# for backward compatibility
upwind dissipation $upwind
# artificial dissipation $ad4
# 
helmholtzForcing
solve Helmholtz 0
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
solve
contour
exit
if( $go eq "go" ){ $cmd .= "movie mode\n exit"; }else{ $cmd="#"; }
$cmd



# --- end cgWave setup  ---
show file $show
contour
exit
if( $solver eq "fixedPoint" ){ $cmd="compute with fixed-point"; }elsif( $solver eq "krylov" ){ $cmd="compute with petsc"; }else{ $cmd="#" }; 
if( $go eq "go" && $cmd ne "#" ){ $cmd .= "\n exit"; }
$cmd