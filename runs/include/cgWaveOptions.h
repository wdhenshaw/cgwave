#
#  GENERIC CGWAVE COMMAND FILE -- INCLUDED IN OTHER SPECIFIC COMMAND FILES
#
# -bcApproach=[cbc|lcbc|oneSided]
# -bc=[d|n|e]                                : set all boundaries to d=Dirichlet, n=Neumann, e=exact
# -bc[1|2|3|4|...] = [d|n|e]                 : -bc1=d : set boundaries with boundary condition 1 to dirichlet
# -cfl=<f>
# -chooseTimeStepFromExplicitGrids=[0|1]     : 1=choose dt from explicit grids and cfl or from implicit grids if all are implicit
# -damp=<f>                                  : coefficient of linear damping
# -debug=[0|1|3|7|15|...]                    : debug bit flag, debug = 1 + 2 + 4 + ...
# -degreeInSpace=<i>
# -degreeInTime=<i> 
# -dtMax=<f>                                 : maximum time step
# -go=[halt|go|og]                           : 
# -implicitUpwind=[0|1]
# -meApproach=[std|ha|stencil]
# -orderInTime=[-1|2|4|6|8]                   :-1 = use default
# -rectangular=[explicit|implicit]            : for ts=implicit, set rectangular=explicit to treat rectangular grids explicitly 
# -takeImplicitFirstStep=[0|1]                :  
# -tf=<f>                                     : final time
# -tp=<f>                                     : times to plot/print 
# -ts=[explicit|implicit]                     : time-stepping method
# -tz=[polyl|trig]
# -upwind=[0|1]
# -useKnownFirstStep=[0|1]
# 
# implicit time-stepping options: (used in implicitOptions.h)
# 
#

if( $bcApproach eq ""                      ){ $bcApproach="oneSided"; }
if( $bc eq ""                              ){ $bc="d"; }
if( $cfl eq ""                             ){ $cfl=0.95; }
if( $chooseTimeStepFromExplicitGrids eq "" ){ $chooseTimeStepFromExplicitGrids=0; }
if( $damp eq ""                            ){ $damp=0.; }
if( $debug eq ""                           ){ $debug=0; }
if( $degreeInSpace eq ""                   ){ $degreeInSpace=2; }
if( $degreeInTime eq ""                    ){ $degreeInTime=2; }
if( $dtMax eq ""                           ){ $dtMax=1e10; }
if( $go eq ""                              ){ $go="halt"; }
if( $meApproach eq ""                      ){ $meApproach="std"; }
if( $orderInTime eq ""                     ){ $orderInTime=-1; }
if( $rectangular eq ""                     ){ $rectangular="implicit"; }
if( $takeImplicitFirstStep eq ""           ){ $takeImplicitFirstStep=0; }
if( $tf eq ""                              ){ $tf=1.; }
if( $tp eq ""                              ){ $tp=.1; }
if( $ts eq ""                              ){ $ts="explicit"; }
if( $tz eq ""                              ){ $tz="poly"; }
if( $upwind eq ""                          ){ $upwind=0; }
if( $useKnownFirstStep eq ""               ){ $useKnownFirstStep=0; }


GetOptions( "tz=s"=>\$tz,"degreeInSpace=i"=>\$degreeInSpace, "degreeInTime=i"=>\$degreeInTime,"cfl=f"=>\$cfl,\
            "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"debug=i"=>\$debug,"orderInTime=i"=>\$orderInTime,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"bc=s"=>\$bc,"tf=f"=>\$tf,"tp=f"=>\$tp,"ts=s"=>\$ts,"dtMax=f"=>\$dtMax,\
            "fx=f"=>\$fx,"fy=f"=>\$fy,"fz=f"=>\$fz,"ft=f"=>\$ft,"rectangular=s"=>\$rectangular,\
            "beta2=f"=>\$beta2,"beta4=f"=>\$beta4,"beta6=f"=>\$beta6,"upwind=i"=>\$upwind,"bcApproach=s"=>\$bcApproach,\
            "useKnownFirstStep=i"=>\$useKnownFirstStep,"meApproach=s"=>\$meApproach,"implicitUpwind=i"=>\$implicitUpwind,\
            "bc1=s"=>\$bc1,"bc2=s"=>\$bc2,"bc3=s"=>\$bc3,"bc4=s"=>\$bc4,"bc5=s"=>\$bc5,"bc6=s"=>\$bc6,\
            "solveri=s"=>\$solveri,"rtoli=f"=>\$rtoli,"atoli=f"=>\$atoli,"maxiti=i"=>\$maxiti,"debugmg=i"=>\$debugmg,"debugOges=i"=>\$debugOges,\
            "takeImplicitFirstStep=i"=>\$takeImplicitFirstStep,"damp=f"=>\$damp,"chooseTimeStepFromExplicitGrids=i"=>\$chooseTimeStepFromExplicitGrids,\
            "go=s"=>\$go );
#


$debug=3; $debugmg=1; $debugOges=0; 
$beta2=.5; $beta4=0.; $beta6=0.; $beta8=0.; # weights in implicit time-stepping 
$dtMax=1e10; $damp=0; 
$bc="d"; $bc1=""; $bc2=""; $bc3=""; $bc4=""; $bc5=""; $bc6=""; 


$fx=2.; $fy=-1; $fz=-1; $ft=-1; # -1 : set equal to $fx
$solveri="yale"; $maxiti=2000; $rtoli=1.0e-10; $atoli=1.0e-10; # parameters for implicit time-stepping solver




# $aa=1./15.;
# printf("aa=%24.20e\n",$aa);
$omega=30.1; $x0=0; $y0=0; $z0=0; $beta=400; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; 
# $ad4=0; # OLD
$upwind=0; # new
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




