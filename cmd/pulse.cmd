#
#  cgWave: test a Gaussian pulse
#        cgWave gridName -cmd=pulse.cmd -bc=[d|n|e|r]
#
$omega=30.1; $x0=0; $y0=0; $z0=0; $beta=100; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; $ad4=0; $debug=3; 
$tz="polynomial"; $tf=5.; $tp=.1; $bc="d"; $go="halt"; 
$cfl=.9; 
$degreeInSpace=2; $degreeInTime=2; 
GetOptions( "tz=s"=>\$tz,"degreeInSpace=i"=>\$degreeInSpace, "degreeInTime=i"=>\$degreeInTime,"cfl=f"=>\$cfl,"ad4=f"=>\$ad4,\
            "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"debug=i"=>\$debug,\
            "bc=s"=>\$bc,"tf=f"=>\$tf,"tp=f"=>\$tp, "go=s"=>\$go );
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
Gaussian params $beta $x0 $y0 0 (beta,x0,y0,z0)
artificial dissipation $ad4
exit
#
solve
contour
exit
if( $go eq "go" ){ $cmd="movie mode\n exit\n exit"; }else{ $cmd="#"; }
$cmd



