#
#  cgWave: test twilightZone
#     cgWave [-noplot] tz.cmd -g=<grid-name> -bc=[d|n] -cfl=<f> -tz=[polyl|trig] -degreeInSpace=<i> -degreeInTime=<i> ...
#                              -upwind=[0|1] -fx= -fy= -fz= -ft= -ts=[explicit|implicit] -rectangular=[implicit|explicit] -debug=<i>
#                              -bcApproach=[cbc|lcbc|oneSided]
#
# $aa=1./15.;
# printf("aa=%24.20e\n",$aa);
$omega=30.1; $x0=0; $y0=0; $z0=0; $beta=400; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; 
# $ad4=0; # OLD
$upwind=0; # new
$debug=3; 
$ts="explicit"; 
$rectangular="implicit"; # for ts=implicit, set rectangular=explicit to treat rectangular grids explicitly
$beta2=.5; $beta4=0.; $beta6=0.; $beta8=0.; # weights in implicit time-stepping 
$dtMax=1e10; 
$bc="d"; 
$bcApproach="oneSided"; # bc Approach : cbc, lcbc, oneSided
$useKnownFirstStep=0; 
$orderInTime=-1;  # -1 = use default
$tz="polynomial"; 
$tf=1.; $tp=.1; $cfl=.9; $go="halt; "
$degreeInSpace=2; $degreeInTime=2; $fx=2.; $fy=2.; $fz=2.; $ft=2.; 
GetOptions( "tz=s"=>\$tz,"degreeInSpace=i"=>\$degreeInSpace, "degreeInTime=i"=>\$degreeInTime,"cfl=f"=>\$cfl,\
            "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"debug=i"=>\$debug,"orderInTime=i"=>\$orderInTime,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"bc=s"=>\$bc,"tf=f"=>\$tf,"tp=f"=>\$tp,"ts=s"=>\$ts,"dtMax=f"=>\$dtMax,\
            "fx=f"=>\$fx,"fy=f"=>\$fy,"fz=f"=>\$fz,"ft=f"=>\$ft,"rectangular=s"=>\$rectangular,\
            "beta2=f"=>\$beta2,"beta4=f"=>\$beta4,"beta6=f"=>\$beta6,"upwind=i"=>\$upwind,"bcApproach=s"=>\$bcApproach,\
            "useKnownFirstStep=i"=>\$useKnownFirstStep,"go=s"=>\$go );
#
if( $tz eq "trig" ){ $tz="trigonometric"; }
if( $tz eq "poly" ){ $tz="polynomial"; }
# 
# time-stepping: (explicit or implicit)
$ts
#
cfl $cfl 
tPlot $tp 
tFinal $tf
dtMax $dtMax
#
$cmd="#";
if( $bcApproach eq "oneSided" ){ $cmd="useOneSidedBCs"; }
if( $bcApproach eq "cbc"      ){ $cmd="useCompatibilityBCs"; }
if( $bcApproach eq "lcbc"     ){ $cmd="useLocalCompatibilityBCs"; }
$cmd
#
if( $ts eq "implicit" ){ $cmd="choose grids for implicit\n  rectangular=$rectangular\n done"; }else{ $cmd="#"; }
$cmd
#
implicit weights $beta2 $beta4 $beta6 $beta8
#
debug $debug 
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
$cmd
#
if( $bc eq "d" ){ $cmd="bc=dirichlet"; }else{ $cmd="bc=neumann"; }
$cmd
turn on forcing 1
twilightZoneForcing
$tz 
#
use known for first step $useKnownFirstStep
#
degreeInSpace $degreeInSpace
degreeInTime $degreeInTime
#
trig frequencies $fx $fy $fz $ft 
#
#Gaussian params $beta $x0 $y0 0 (beta,x0,y0,z0)
# omega $omega
# omegaSOR $omegaSOR
# if( $ad4>0. ){ $upwind=1; }# for backward compatibility
upwind dissipation $upwind
# artificial dissipation $ad4
# tol $tol 
# number of periods $numPeriods
exit
#
solve
contour
exit
if( $go eq "go" ){ $cmd="exit"; }else{ $cmd="#"; }
$cmd



