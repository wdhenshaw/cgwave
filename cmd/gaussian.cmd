#
#  cgWaveHoltz script: 
#
$omega=30.1; $x0=0; $y0=0; $z0=0; $beta=400; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; 
$ad4=1; # old
$upwind=1; 
$cfl=.7; 
GetOptions( "omega=f"=>\$omega,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"numPeriods=i"=>\$numPeriods,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"ad4=f"=>\$ad4,"cfl=f"=>\$cfl,"upwind=i"=>\$upwind );
# 
# pause
cfl $cfl 
Gaussian params $beta $x0 $y0 0 (beta,x0,y0,z0)
omega $omega
omegaSOR $omegaSOR
if( $ad4>0. ){ $upwind=1; }# for backward compatibility
upwind dissipation $upwind
# artificial dissipation $ad4
tol $tol 
number of periods $numPeriods
exit
#


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
