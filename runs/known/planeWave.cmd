#
#  cgWave: test a plane wave
#
$amp=1; $kx=1.0; $ky=0; $kz=0; $ad4=0; $debug=3;  $go="halt"; $bc="d"; $dissFreq=1; 
$tf=5.; $tp=.05; $cfl=.9; 
$degreeInSpace=2; $degreeInTime=2; 
$setKnownOnBoundaries=1; # set known solution on boundaries
GetOptions( "cfl=f"=>\$cfl,"ad4=f"=>\$ad4,"amp=f"=>\$amp,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"debug=i"=>\$debug,\
            "tf=f"=>\$tf,"tp=f"=>\$tp,"bc=s"=>\$bc,"dissFreq=i"=>\$dissFreq,\
            "setKnownOnBoundaries=i"=>\$setKnownOnBoundaries,"go=s"=>\$go);
# 
#
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
# pause
tFinal $tf
tPlot $tp
cfl $cfl 
debug $debug
#
bc=$bc
#
user defined known solution...
  plane wave
    $amp $kx $ky $kz
done
set known on boundaries $setKnownOnBoundaries
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


exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit

exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
exit
