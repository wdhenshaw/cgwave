#
#  cgWave: test twilightZone
#     cgWave [-noplot] tz.cmd -g=<grid-name> -bc=[d|n|e] -cfl=<f> -tz=[polyl|trig] -degreeInSpace=<i> -degreeInTime=<i> ...
#                              -upwind=[0|1] -fx= -fy= -fz= -ft= -ts=[explicit|implicit] -rectangular=[implicit|explicit] -debug=<i>
#                              -bcApproach=[cbc|lcbc|oneSided] -meApproach=[std|ha|stencil] -implicitUpwind=[0|1]
#
# $aa=1./15.;
# printf("aa=%24.20e\n",$aa);
$omega=30.1; $x0=0; $y0=0; $z0=0; $beta=400; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; 
$upwind=0; $nuc=1; # new
$debug=3; $debugmg=1; $debugOges=0; 
$ts="explicit"; $implicitUpwind=0;
$chooseTimeStepFromExplicitGrids=1; # 1=choose dt from explicit grids and cfl unless all grids are implicit 
$rectangular="implicit"; # for ts=implicit, set rectangular=explicit to treat rectangular grids explicitly
$beta2=.5; $beta4=0.; $beta6=0.; $beta8=0.; # weights in implicit time-stepping 
$dtMax=1e10; $damp=0; 
$bc="d"; $bc1=""; $bc2=""; $bc3=""; $bc4=""; $bc5=""; $bc6=""; 
$bcApproach="oneSided"; # bc Approach : cbc, lcbc, oneSided
$meApproach="std"; # or "ha"
$useKnownFirstStep=0; $takeImplicitFirstStep=0; 
$orderInTime=-1;  # -1 = use default
$tz="polynomial"; 
$tf=1.; $tp=.1; $cfl=.9; $go="halt; "
$degreeInSpace=2; $degreeInTime=2; 
$fx=2.; $fy=-1; $fz=-1; $ft=-1; # -1 : set equal to $fx
$solveri="yale"; $maxiti=2000; $rtoli=1.0e-10; $atoli=1.0e-10; # parameters for implicit time-stepping solver
#
GetOptions( "tz=s"=>\$tz,"degreeInSpace=i"=>\$degreeInSpace, "degreeInTime=i"=>\$degreeInTime,"cfl=f"=>\$cfl,\
            "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"debug=i"=>\$debug,"orderInTime=i"=>\$orderInTime,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"bc=s"=>\$bc,"tf=f"=>\$tf,"tp=f"=>\$tp,"ts=s"=>\$ts,"dtMax=f"=>\$dtMax,\
            "fx=f"=>\$fx,"fy=f"=>\$fy,"fz=f"=>\$fz,"ft=f"=>\$ft,"rectangular=s"=>\$rectangular,"nuc=i"=>\$nuc,\
            "beta2=f"=>\$beta2,"beta4=f"=>\$beta4,"beta6=f"=>\$beta6,"upwind=i"=>\$upwind,"bcApproach=s"=>\$bcApproach,\
            "useKnownFirstStep=i"=>\$useKnownFirstStep,"meApproach=s"=>\$meApproach,"implicitUpwind=i"=>\$implicitUpwind,\
            "bc1=s"=>\$bc1,"bc2=s"=>\$bc2,"bc3=s"=>\$bc3,"bc4=s"=>\$bc4,"bc5=s"=>\$bc5,"bc6=s"=>\$bc6,\
            "solveri=s"=>\$solveri,"rtoli=f"=>\$rtoli,"atoli=f"=>\$atoli,"maxiti=i"=>\$maxiti,"debugmg=i"=>\$debugmg,"debugOges=i"=>\$debugOges,\
            "takeImplicitFirstStep=i"=>\$takeImplicitFirstStep,"damp=f"=>\$damp,"chooseTimeStepFromExplicitGrids=i"=>\$chooseTimeStepFromExplicitGrids,\
            "go=s"=>\$go );
#
if( $tz eq "trig" ){ $tz="trigonometric"; }
if( $tz eq "poly" ){ $tz="polynomial"; }
if( $fy eq -1 ){ $fy=$fx; }
if( $fz eq -1 ){ $fz=$fx; }
if( $ft eq -1 ){ $ft=$fx; }
# 
# time-stepping: (explicit or implicit)
$ts
#
choose time step from explicit grids $chooseTimeStepFromExplicitGrids
#
cfl $cfl 
tPlot $tp 
tFinal $tf
dtMax $dtMax
damp $damp
numUpwindCorrections $nuc
#
$cmd="#";
if( $bcApproach eq "oneSided" ){ $cmd="useOneSidedBCs"; }
if( $bcApproach eq "cbc"      ){ $cmd="useCompatibilityBCs"; }
if( $bcApproach eq "lcbc"     ){ $cmd="useLocalCompatibilityBCs"; }
$cmd
#
$cmd="#";
if( $meApproach eq "std" ){ $cmd="standard modified equation"; }
if( $meApproach eq "ha" ){ $cmd="hierarchical modified equation"; }
if( $meApproach eq "stencil" ){ $cmd="stencil modified equation"; }
# printf("meApproach=$meApproach\n");
# printf("cmd=$cmd\n");
# pause
$cmd
#
#  Set options for implicit time-stepping: 
if( $ts eq "implicit" ){ $cmd="include $ENV{CGWAVE}/runs/include/implicitOptions.h"; }else{ $cmd="#"; }
$cmd
take implicit first step $takeImplicitFirstStep
#
debug $debug 
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
$cmd
#
if( $bc eq "d" ){ $cmd="bc=dirichlet"; }elsif( $bc eq "n" ){ $cmd="bc=neumann"; }elsif( $bc eq "e" ){ $cmd="bc=exact"; }elsif( $bc eq "a" ){ $cmd="bc=absorbing"; }else{ $cmd="bc=dirichlet"; }
$cmd
$cmd="#"; 
if( $bc1 ne "" ){ $cmd .="\n bcNumber1=$bc1"; }
if( $bc2 ne "" ){ $cmd .="\n bcNumber2=$bc2"; }
if( $bc3 ne "" ){ $cmd .="\n bcNumber3=$bc3"; }
if( $bc4 ne "" ){ $cmd .="\n bcNumber4=$bc4"; }
if( $bc5 ne "" ){ $cmd .="\n bcNumber5=$bc5"; }
if( $bc6 ne "" ){ $cmd .="\n bcNumber6=$bc6"; }
$cmd
#
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
implicit upwind $implicitUpwind
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




