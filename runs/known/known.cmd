#
#  cgWave: Compute to some "known" solutions
#
#   cgWave [-noplot] known.cmd -g=<grid-name> -known=[pw|gpw|boxHelmholtz|polyPeriodic|diskEig|annulusEig] 
#           -upwind=[0|1] -computeErrors=[0|1]
#           -setKnownOnBoundaries=[0|1] -bcApproach=[cbc|lcbc|oneSided]
#
#   pw = plane wave 
#   gpw = Gaussian plane wave 
$known="planeWave";
$amp=1; $kx=1.0; $ky=0; $kz=0; $omega=3.; 
$beta=20; $x0=.5; $y0=.0; $z0=.0; $k0=0.; # for Gaussian plane wave 
# $ad4=0;   # old
$upwind=0; 
$debug=3;  $go="halt"; $bc="d"; $dissFreq=1; 
$bcApproach="oneSided"; # bc Approach : cbc, lcbc, oneSided
$useKnownFirstStep=0; 
$computeErrors=1; $setKnownOnBoundaries=1; 
$tf=5.; $tp=.05; $cfl=.9; 
$ts="explicit"; $dtMax=1e10; 
$orderInTime=-1;  # -1 = use default
$degreeInSpace=2; $degreeInTime=2; 
$nBessel=1; $mTheta=1; 
$show=""; $flushFrequency=10; 
GetOptions( "cfl=f"=>\$cfl,"amp=f"=>\$amp,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"debug=i"=>\$debug,\
            "tf=f"=>\$tf,"tp=f"=>\$tp,"bc=s"=>\$bc,"dissFreq=i"=>\$dissFreq,"omega=f"=>\$omega,\
            "known=s"=>\$known,"orderInTime=i"=>\$orderInTime,"ts=s"=>\$ts,"dtMax=f"=>\$dtMax,"upwind=i"=>\$upwind,\
            "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"k0=f"=>\$k0,"beta=f"=>\$beta,"computeErrors=i"=>\$computeErrors,\
            "setKnownOnBoundaries=s"=>\$setKnownOnBoundaries,"show=s"=>\$show,"useKnownFirstStep=i"=>\$useKnownFirstStep,\
            "flushFrequency=i"=>\$flushFrequency,"bcApproach=s"=>\$bcApproach,"nBessel=i"=>\$nBessel,"mTheta=i"=>\$mTheta,\
            "go=s"=>\$go );
# 
#
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
if( $known eq "diskEig" ){ $setKnownOnBoundaries=0; }
if( $known eq "annulusEig" ){ $setKnownOnBoundaries=0; }
# time-stepping: (explicit or implicit)
$ts
# pause
tFinal $tf
tPlot $tp
cfl $cfl 
dtMax $dtMax
debug $debug
omega $omega
set known on boundaries $setKnownOnBoundaries
use known for first step $useKnownFirstStep
compute errors $computeErrors
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
$cmd
if( $show ne "" ){ $cmd="show file name $show\n save show file 1\n flush frequency $flushFrequency"; }else{ $cmd="#"; }
$cmd
# 
$cmd="#";
if( $bcApproach eq "oneSided" ){ $cmd="useOneSidedBCs"; }
if( $bcApproach eq "cbc"      ){ $cmd="useCompatibilityBCs"; }
if( $bcApproach eq "lcbc"     ){ $cmd="useLocalCompatibilityBCs"; }
$cmd
#
bc=$bc
#
if( $known eq "boxHelmholtz" ){ $cmd="helmholtzForcing\n user defined forcing...\n box Helmholtz\n exit"; }else{ $cmd="#"; }
if( $known eq "polyPeriodic" ){ $cmd="helmholtzForcing\n user defined forcing...\n poly periodic\n exit"; }
$cmd
#
user defined known solution...
  if( $known eq "gpw" ){ $cmd="gaussian plane wave\n wave numbers: $kx $ky $kz\n beta: $beta\n k0: $k0\n offset: $x0 $y0 $z0"; }
  if( $known eq "planeWave" || $known eq "pw" ){ $cmd="plane wave\n $amp $kx $ky $kz"; }
  if( $known eq "boxHelmholtz"){ $cmd="box helmholtz\n $omega $kx $ky $kz"; }
  $degreeInSpaceForPolyPeriodic=1; 
  if( $known eq "polyPeriodic"){ $cmd="poly periodic\n $omega $degreeInSpaceForPolyPeriodic"; }
  # -- disk eigenfunction: 
  # n,m,a,amp,bcOpt 
  $rad=1; 
  if( $bc eq "dirichlet" ){ $bcOpt=0; }else{ $bcOpt=1; } 
  if( $known eq "diskEig" ){ $cmd="disk eigenfunction\n $nBessel $mTheta $rad $amp $bcOpt"; }  
  if( $known eq "annulusEig" ){ $cmd="annulus eigenfunction\n $nBessel $mTheta $amp $bcOpt"; }  
  $cmd
done
#
# if( $ad4>0. ){ $upwind=1; }# for backward compatibility
upwind dissipation $upwind
# artificial dissipation $ad4
dissipation frequency $dissFreq
exit
#
solve
#
contour
  #   plot contour lines (toggle)
exit
#
if( $go eq "go" ){ $cmd="exit"; }else{ $cmd="#"; }
$cmd
