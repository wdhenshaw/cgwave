#
#  cgWaveHoltz script: 
#   known = [boxHelmholtz|polyPeriodic|computedHelmholtz]
#   solver = [fixedPoint|krylov] : fixed-point or Krylov 
#   imode =1 : do not wait in cgWave
#
$known="boxHelmholtz"; $setKnownOnBoundaries=0;  $degreeInSpaceForPolyPeriodic=1; 
$omega=30.1; $beta=50.; $x0=0.5; $y0=0.5; $z0=0.5; $t0=0.; $go="halt"; $matlab="cgWaveHoltz"; 
$beta=400; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; 
$numberOfFrequencies=1; 
@freq = ();  # this must be null for GetOptions to work, defaults are given below
$ts="explicit"; $dtMax=1e10; 
# beta2=0 = TRAP, beta2=.5 = Full weighting 
$beta2=.0; $beta4=0.; $beta6=0.; $beta8=0.; # weights in implicit time-stepping 
# $ad4=0; # old way
$upwind=0; # new way
$tp=.5; $imode=0; $adjustOmega=0; $show="trigHelmholtz.show"; 
$solver="fixedPoint";  $kx=1; $ky=1; $kz=1; $maxIterations=100; $orderInTime=-1;  # -1 = use default
$cfl=.9; $bc="d"; 
$bcApproach="oneSided"; # bc Approach : cbc, lcbc, oneSided
GetOptions( "omega=f"=>\$omega,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"numPeriods=i"=>\$numPeriods,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"cfl=f"=>\$cfl,"tp=f"=>\$tp,"iMode=i"=>\$imode,\
            "solver=s"=>\$solver,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"adjustOmega=i"=>\$adjustOmega,"known=s"=>\$known,\
            "matlab=s"=>\$matlab,"maxIterations=i"=>\$maxIterations,"upwind=i"=>\$upwind,"bcApproach=s"=>\$bcApproach,\
            "degreeInSpaceForPolyPeriodic=i"=>\$degreeInSpaceForPolyPeriodic,"orderInTime=i"=>\$orderInTime,\
            "beta2=f"=>\$beta2,"beta4=f"=>\$beta4,"beta6=f"=>\$beta6,"ts=s"=>\$ts,"dtMax=f"=>\$dtMax,\
            "numberOfFrequencies=i"=>\$numberOfFrequencies,"nf=i"=>\$numberOfFrequencies,"freq=f{1,}"=>\@freq,\
            "show=s"=>\$show,"go=s"=>\$go );
# 
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
if( $freq[0] eq "" ){ @freq=(1,2,3,4,5,6,7,8,9,10); }
# 
$omega = $freq[0]; 
# pause
# -------- Start CgWaveHoltz Options ------
# debug $debug
Gaussian params $beta $x0 $y0 0 (beta,x0,y0,z0)
omega $omega
omegaSOR $omegaSOR
tol $tol 
number of frequencies $numberOfFrequencies
frequencies $freq[0] $freq[1] $freq[2] $freq[3] $freq[4] $freq[5] $freq[6] $freq[7] $freq[8]
number of periods $numPeriods
adjust omega $adjustOmega
matlab filename: $matlab
exit
# ------ Start cgWave Options ------
interactiveMode $imode 
# time-stepping: (explicit or implicit)
$ts
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
$cmd
cfl $cfl 
tPlot $tp 
dtMax $dtMax
implicit weights $beta2 $beta4 $beta6 $beta8
if( $known eq "polyPeriodic" ){ $setKnownOnBoundaries=1; }
set known on boundaries $setKnownOnBoundaries
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
#
# helmholtzForcing
solve Helmholtz 1
#
adjust omega $adjustOmega
#   --- eigenmode-like solution 
user defined known solution...
  if( $known eq "boxHelmholtz"){ $cmd="box Helmholtz\n $omega $kx $ky $kz"; }
  if( $known eq "computedHelmholtz"){ $cmd="computed Helmholtz\n $omega $kx $ky $kz"; }
  if( $known eq "polyPeriodic"){ $cmd="poly periodic\n $omega $degreeInSpaceForPolyPeriodic"; }
  $cmd
  # box helmholtz
  # $omega $kx $ky $kz
done
$cmd="#"; 
if( ($known eq "boxHelmholtz") || ($known eq "computedHelmholtz") ){ $cmd="helmholtzForcing\n user defined forcing...\n box Helmholtz\n exit"; }
if( $known eq "polyPeriodic" ){ $cmd="helmholtzForcing\n user defined forcing...\n poly periodic\n exit"; }
$cmd
# user defined forcing...
#   box Helmholtz
# exit
exit
# --- end cgWave ---
show file $show
max iterations $maxIterations
if( $solver eq "fixedPoint" ){ $cmd="compute with fixed-point"; }else{ $cmd="#"; }
if( $solver eq "krylov" ){ $cmd="compute with petsc"; }
$cmd 
contour
exit
if( $go eq "go" ){ $cmd="exit"; }else{ $cmd="#"; }
$cmd 




contour
exit
exit


# Non-interactive mode -- just run and exit 
if( $imode eq 1 ){ $cmd="movie mode\n exit"; }else{ $cmd="#"; }
$cmd
