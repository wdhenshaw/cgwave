#
#  cgWave: test a solution to the Helmholtz equation
# 
#        cgWave gridName -cmd=helmholtz.cmd -kx= -ky= -kz = -bc=[d|n|e|r] -bcApproach=[cbc|lcbc|oneSided]
#
$ad4=0; $debug=3;
$tf=-1.; # adjusted below if not set
$tp=.1; $bc="d"; $go="halt"; 
$omega=2;  $kx=1.0; $ky=1; $kz=1; $solveHelmholtz=0; $computeTimeIntegral=1; 
$cfl=.9; 
$pi = atan2(1.,1.)*4.; 
$degreeInSpace=2; $degreeInTime=2;
$bcApproach="oneSided"; # bc Approach : cbc, lcbc, oneSided
$ts="explicit";
$dtMax=1e10; 
$orderInTime=-1;  # -1 = use default 
#
GetOptions( "cfl=f"=>\$cfl,"omega=f"=>\$omega,"ad4=f"=>\$ad4,"debug=i"=>\$debug,"solveHelmholtz=i"=>\$solveHelmholtz,\
            "bc=s"=>\$bc,"tf=f"=>\$tf,"tp=f"=>\$tp,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz, \
            "ts=s"=>\$ts,"dtMax=f"=>\$dtMax,"orderInTime=i"=>\$orderInTime,\
            "computeTimeIntegral=i"=>\$computeTimeIntegral,"bcApproach=s"=>\$bcApproach,"go=s"=>\$go );
# 
if( $tf < 0. ){ $tf = 2.*$pi/$omega; }
#
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
# time-stepping: (explicit or implicit)
$ts
# pause
tFinal $tf
tPlot $tp
cfl $cfl 
debug $debug
dtMax $dtMax
#
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
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
user defined known solution...
box Helmholtz
  $omega $kx $ky $kz
done
helmholtzForcing
solve Helmholtz $solveHelmholtz
compute time integral $computeTimeIntegral 
#
user defined forcing...
  box Helmholtz
exit
set known on boundaries 1
# 
artificial dissipation $ad4
exit
#
solve
# contour
# exit
if( $go eq "go" ){ $cmd="exit"; }else{ $cmd="#"; }
$cmd



