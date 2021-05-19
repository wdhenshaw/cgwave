#
#  cgWave: Compute to some "known" solutions
#
#   cgWave [-noplot] known.cmd -g=<grid-name> -known=[planeAve|boxHelmholtz]
#
$known="planeWave";
$amp=1; $kx=1.0; $ky=0; $kz=0; $omega=3.; $ad4=0; $debug=3;  $go="halt"; $bc="d"; $dissFreq=1; 
$tf=5.; $tp=.05; $cfl=.9; 
$ts="explicit"; $dtMax=1e10; 
$orderInTime=-1;  # -1 = use default
$degreeInSpace=2; $degreeInTime=2; 
GetOptions( "cfl=f"=>\$cfl,"ad4=f"=>\$ad4,"amp=f"=>\$amp,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"debug=i"=>\$debug,\
            "tf=f"=>\$tf,"tp=f"=>\$tp,"bc=s"=>\$bc,"dissFreq=i"=>\$dissFreq,"omega=f"=>\$omega,\
            "known=s"=>\$known,"orderInTime=i"=>\$orderInTime,"ts=s"=>\$ts,"dtMax=f"=>\$dtMax, "go=s"=>\$go );
# 
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
dtMax $dtMax
debug $debug
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
$cmd
#
bc=$bc
#
if( $known eq "boxHelmholtz" ){ $cmd="helmholtzForcing\n user defined forcing...\n box Helmholtz\n exit"; }else{ $cmd="#"; }
$cmd
#
user defined known solution...
  if( $known eq "planeWave"){ $cmd="plane wave\n $amp $kx $ky $kz"; }
  if( $known eq "boxHelmholtz"){ $cmd="box helmholtz\n $omega $kx $ky $kz"; }
  $cmd
done
#
artificial dissipation $ad4
dissipation frequency $dissFreq
exit
#
solve
#
contour
exit
#
if( $go eq "go" ){ $cmd="exit"; }else{ $cmd="#"; }
$cmd
