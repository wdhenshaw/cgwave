#
#  cgWave: Compute solutions with a random initial condition
#
#     cgwave random.cmd -g=<grid-name> -x0=<f> -y0=<f> -omega=<f> -tol=<f> -tp=<f> ...
#                         -kx=<f> -ky=<f> -kz=<f> -ic=[gaussian] -upwind=[0|1] -imode=[0|1] 
#                         -bcApproach=[cbc|lcbc|oneSided] -assignInterpNeighbours=[extrap|interp]
#                         -go=[go|og|halt]
#
#   -imode=1 : do not wait in cgWave
#
$go="go"; $ic="gaussian"; 
$omega=30.1; $beta=50.; $x0=0.5; $y0=0.5; $z0=0.5; $t0=0.; $amp=1.; $debug=0; 
$ad4=0;    # old way
$upwind=1; # new way
$tf=5.; $tp=.5; $imode=1; 
$kx=1; $ky=1; $kz=1; 
$bcApproach="oneSided"; # bc Approach : cbc, lcbc, oneSided
$matlab="cgWave"; $show="gaussian.show"; $flushFrequency=1000; $computeEnergy=0;
$cfl=.9; $bc="d"; $ts="explicit"; $dtMax=1;
$rectangular="implicit"; # for ts=implicit, set rectangular=explicit to treat rectangular grids explicitly 
$takeImplicitFirstStep=0;
$orderInTime=-1;  # -1 = use default
$assignInterpNeighbours="extrap"; # by default extrap interp neighbours 
$chooseTimeStepFromExplicitGrids=1; # 1=choose dt from explicit grids and cfl unless all grids are implicit
$solveri="yale"; $maxiti=2000; $rtoli=1.0e-10; $atoli=1.0e-10; # parameters for implicit time-stepping solver
#
GetOptions( "omega=f"=>\$omega,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"numPeriods=i"=>\$numPeriods,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"ad4=f"=>\$ad4,"cfl=f"=>\$cfl,"tp=f"=>\$tp,"tf=f"=>\$tf,"iMode=i"=>\$imode,\
            "solver=s"=>\$solver,"k0=f"=>\$kx,"k0=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"maxIterations=i"=>\$maxIterations,"matlab=s"=>\$matlab,\
            "go=s"=>\$go,"ic=s"=>\$ic,"bc=s"=>\$bc,"ts=s"=>\$ts,"orderInTime=i"=>\$orderInTime,"rectangular=s"=>\$rectangular,\
            "dtMax=f"=>\$dtMax,"adjustOmega=i"=>\$adjustOmega,"amp=f"=>\$amp,"show=s"=>\$show,"upwind=i"=>\$upwind,\
            "debug=i"=>\$debug,"bcApproach=s"=>\$bcApproach,"assignInterpNeighbours=s"=>\$assignInterpNeighbours,\
            "takeImplicitFirstStep=i"=>\$takeImplicitFirstStep,"chooseTimeStepFromExplicitGrids=i"=>\$chooseTimeStepFromExplicitGrids,\
            "solveri=s"=>\$solveri,"rtoli=f"=>\$rtoli,"atoli=f"=>\$atoli,"maxiti=i"=>\$maxiti,\
            "computeEnergy=i"=>\$computeEnergy );
# 
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
# 
# ------ Start cgWave setup ------
# time-stepping : explicit/implicit
$ts
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
$cmd
#
choose time step from explicit grids $chooseTimeStepFromExplicitGrids
#
if( $show ne "" ){ $cmd="show file name $show\n save show file 1\n flush frequency $flushFrequency"; }else{ $cmd="#"; }
$cmd
dtMax $dtMax
cfl $cfl
# interactiveMode $imode 
debug $debug
tPlot $tp 
tFinal $tf
#
compute energy $computeEnergy
# 
$cmd="#";
if( $bcApproach eq "oneSided" ){ $cmd="useOneSidedBCs"; }
if( $bcApproach eq "cbc"      ){ $cmd="useCompatibilityBCs"; }
if( $bcApproach eq "lcbc"     ){ $cmd="useLocalCompatibilityBCs"; }
$cmd
# -- 
bc=$bc
if( $assignInterpNeighbours eq "interp" ){ $cmd="interpolateInterpNeighbours"; }else{ $cmd="#"; }
$cmd
#  Set options for implicit time-stepping: 
if( $ts eq "implicit" ){ $cmd="include $ENV{CGWAVE}/runs/include/implicitOptions.h"; }else{ $cmd="#"; }
$cmd
take implicit first step $takeImplicitFirstStep
#
solve Helmholtz 0
compute errors 0
#
if( $ad4>0. ){ $upwind=1; }# for backward compatibility
upwind dissipation $upwind
#
random initial condition
#
# user defined known solution...
#   offset: $x0, $y0, $z0 (x0,y0,z0)
#   beta: $beta
#   k0: $k0
#   modulated Gaussian
# exit
#
exit
solve
contour
exit
if( $go eq "go" ){ $cmd = "exit\n exit"; }else{ $cmd="#"; }
$cmd



open graphics

contour

if( $go eq "go" ){ $cmd .= "movie mode\n exit"; }else{ $cmd="#"; }
$cmd



# --- end cgWave setup  ---
show file $show
contour
exit
if( $solver eq "fixedPoint" ){ $cmd="compute with fixed-point"; }elsif( $solver eq "krylov" ){ $cmd="compute with petsc"; }else{ $cmd="#" }; 
if( $go eq "go" && $cmd ne "#" ){ $cmd .= "\n exit"; }
$cmd
