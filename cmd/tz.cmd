#
#  cgWave: test twilightZone
#     cgWave [-noplot] tz.cmd -g=<grid-name> -bc=[d|n] -cfl=<f> -tz=[polyl|trig] -degreeInSpace=<i> -degreeInTime=<i> ...
#                             =fx= -fy= -fz= -ft= -ts=[explicit|implicit] -rectangular=[implicit|explicit] -debug=<i>
#
$omega=30.1; $x0=0; $y0=0; $z0=0; $beta=400; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; $ad4=0; $debug=3; 
$ts="explicit"; 
$rectangular="implicit"; # for ts=implicit, set rectangular=explicit to treat rectangular grids explicitly
$cImp1=.25; $cImp2=0.5; $cImp3=.25;  # weights in implicit time-stepping 
$dtMax=1e10; 
$bc="d"; 
$orderInTime=-1;  # -1 = use default
$tz="polynomial"; 
$tf=1.; $tp=.1; $cfl=.9; $go="halt; "
$degreeInSpace=2; $degreeInTime=2; $fx=2.; $fy=2.; $fz=2.; $ft=2.; 
GetOptions( "tz=s"=>\$tz,"degreeInSpace=i"=>\$degreeInSpace, "degreeInTime=i"=>\$degreeInTime,"cfl=f"=>\$cfl,"ad4=f"=>\$ad4,\
            "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"debug=i"=>\$debug,"orderInTime=i"=>\$orderInTime,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"bc=s"=>\$bc,"tf=f"=>\$tf,"tp=f"=>\$tp,"ts=s"=>\$ts,"dtMax=f"=>\$dtMax,\
            "fx=f"=>\$fx,"fy=f"=>\$fy,"fz=f"=>\$fz,"ft=f"=>\$ft,"rectangular=s"=>\$rectangular,\
            "cImp1=f"=>\$cImp1,"cImp2=f"=>\$cImp2,"cImp3=f"=>\$cImp3,"go=s"=>\$go );
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
if( $ts eq "implicit" ){ $cmd="choose grids for implicit\n  rectangular=$rectangular\n done"; }else{ $cmd="#"; }
$cmd
#
implicit weights $cImp1 $cImp2 $cImp3
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
degreeInSpace $degreeInSpace
degreeInTime $degreeInTime
#
trig frequencies $fx $fy $fz $ft 
#
#Gaussian params $beta $x0 $y0 0 (beta,x0,y0,z0)
# omega $omega
# omegaSOR $omegaSOR
artificial dissipation $ad4
# tol $tol 
# number of periods $numPeriods
exit
#
solve
contour
exit
if( $go eq "go" ){ $cmd="exit"; }else{ $cmd="#"; }
$cmd



