#
#  cgWave: Compute to some "known" solutions
#
#   cgWave [-noplot] known.cmd -g=<grid-name> -known=[pw|gpw|boxHelmholtz|polyPeriodic|squareEig|diskEig|annulusEig|sphereEig] 
#           -upwind=[0|1] -computeErrors=[0|1]
#           -setKnownOnBoundaries=[0|1] -bcApproach=[cbc|lcbc|oneSided] -assignInterpNeighbours=[extrap|interp]
#           -meApproach=[std,ha|stencil]
#           -ts=[explicit|implicit] -rectangular=[implicit|explicit]
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
$useKnownFirstStep=0; $checkKnown=0; 
$computeErrors=1; $plotScatteredField=0; 
$setKnownOnBoundaries=-1; #-1 : use default
$tf=5.; $tp=.05; $cfl=.9; 
$ts="explicit"; $dtMax=1e10; $implicitUpwind=0; 
$rectangular="implicit"; # for ts=implicit, set rectangular=explicit to treat rectangular grids explicitly
$beta2=.5; $beta4=0.; $beta6=0.; $beta8=0.; # weights in implicit time-stepping 
$meApproach="std"; # or "ha"
$orderInTime=-1;  # -1 = use default
$degreeInSpace=2; $degreeInTime=2; 
$assignInterpNeighbours="extrap"; # by default extrap interp neighbours 
$nBessel=1; $mTheta=1; 
$mPhi=1; $mr=1; 
$show=""; $flushFrequency=10; 
$solveri="yale"; $maxiti=2000; $rtoli=1.0e-10; $atoli=1.0e-10; # parameters for implicit time-stepping solver
GetOptions( "cfl=f"=>\$cfl,"amp=f"=>\$amp,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"debug=i"=>\$debug,\
            "tf=f"=>\$tf,"tp=f"=>\$tp,"bc=s"=>\$bc,"dissFreq=i"=>\$dissFreq,"omega=f"=>\$omega,"beta2=f"=>\$beta2,"beta4=f"=>\$beta4,\
            "known=s"=>\$known,"orderInTime=i"=>\$orderInTime,"ts=s"=>\$ts,"dtMax=f"=>\$dtMax,"upwind=i"=>\$upwind,\
            "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"k0=f"=>\$k0,"beta=f"=>\$beta,"computeErrors=i"=>\$computeErrors,\
            "setKnownOnBoundaries=s"=>\$setKnownOnBoundaries,"show=s"=>\$show,"useKnownFirstStep=i"=>\$useKnownFirstStep,\
            "flushFrequency=i"=>\$flushFrequency,"bcApproach=s"=>\$bcApproach,"meApproach=s"=>\$meApproach,\
            "nBessel=i"=>\$nBessel,"mPhi=i"=>\$mPhi,"mTheta=i"=>\$mTheta,"mr=i"=>\$mr,"rectangular=s"=>\$rectangular,\
            "solveri=s"=>\$solveri,"rtoli=f"=>\$rtoli,"atoli=f"=>\$atoli,"maxiti=i"=>\$maxiti,"plotScatteredField=i"=>\$plotScatteredField,\
            "assignInterpNeighbours=s"=>\$assignInterpNeighbours,"checkKnown=s"=>\$checkKnown,"implicitUpwind=i"=>\$implicitUpwind,"go=s"=>\$go );
# 
#
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="exact"; }
# if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
# by default we do NOT set known on boundaries for these solutions:
if( $setKnownOnBoundaries eq "-1" && ($known eq "diskEig" || $known eq "squareEig" || $known eq "annulusEig" || $known eq "sphereEig") ){ $setKnownOnBoundaries=0; }
if( $setKnownOnBoundaries eq "-1" ){ $setKnownOnBoundaries=1; }
# if( $known eq "annulusEig" ){ $setKnownOnBoundaries=0; }
# if( $known eq "sphereEig" ){ $setKnownOnBoundaries=0; }
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
plot scattered field $plotScatteredField
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
$cmd="#";
if( $meApproach eq "std" ){ $cmd="standard modified equation"; }
if( $meApproach eq "ha" ){ $cmd="hierarchical modified equation"; }
if( $meApproach eq "stencil" ){ $cmd="stencil modified equation"; }
$cmd
#
if( $ts eq "implicit" ){ $cmd="choose grids for implicit\n  rectangular=$rectangular\n done"; }else{ $cmd="#"; }
$cmd
# implicitUpwind = 1 : include upwinding in implicit matrix
implicit upwind $implicitUpwind
# beta2=0 : trap, beta2=.5 = FW
implicit weights $beta2 $beta4
#
if( $assignInterpNeighbours eq "interp" ){ $cmd="interpolateInterpNeighbours"; }else{ $cmd="#"; }
$cmd
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
  if( $bc eq "exact" ){ $bcOpt=0; }
  if( $known eq "diskEig" ){ $cmd="disk eigenfunction\n $nBessel $mTheta $rad $amp $bcOpt"; }  
  if( $known eq "annulusEig" ){ $cmd="annulus eigenfunction\n $nBessel $mTheta $amp $bcOpt"; }  
  if( $known eq "sphereEig" ){ $cmd="sphere eigenfunction\n $mPhi $mTheta $mr $rad $amp $bcOpt"; }  
  if( $known eq "squareEig" ){ $cmd="square eigenfunction\n $kx $ky $kz"; }  
  $cmd
done
#
# if( $ad4>0. ){ $upwind=1; }# for backward compatibility
upwind dissipation $upwind
#
if( $ts eq "implicit" ){ $cmd="include $ENV{CGWAVE}/runs/include/implicitOptions.h"; }else{ $cmd="#"; }
$cmd
#
# artificial dissipation $ad4
dissipation frequency $dissFreq
exit
if( $checkKnown eq "1" ){ $cmd="check known solution\n "; }else{ $cmd="#"; }
$cmd
#
solve
#
contour
  #   plot contour lines (toggle)
exit
#
if( $go eq "go" ){ $cmd="exit"; }else{ $cmd="#"; }
$cmd
