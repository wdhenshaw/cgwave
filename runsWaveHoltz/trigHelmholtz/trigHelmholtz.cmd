#
#  cgWaveHoltz script: 
#   known = [boxHelmholtz|polyPeriodic|computedHelmholtz]
#   solver = [fixedPoint|krylov] : fixed-point or Krylov 
#   imode =1 : do not wait in cgWave
#   go = [halt|go|direct|krylov|both]
#      go=dk  : direct and krylov
#      go=dfk : solve for direct, fixed-point, krylov
#
$known="boxHelmholtz"; $setKnownOnBoundaries=0;  $degreeInSpaceForPolyPeriodic=1; 
$omega=30.1; $beta=50.; $x0=0.5; $y0=0.5; $z0=0.5; $t0=0.; $go="halt"; $matlab="cgWaveHoltz"; 
$beta=400; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; 
$numberOfFrequencies=1; 
@freq = ();  # this must be null for GetOptions to work, defaults are given below
$ts="explicit"; $dtMax=1e10; 
$minStepsPerPeriod=10;   # min is 5
# beta2=0 = TRAP, beta2=.5 = Full weighting 
$beta2=.0; $beta4=0.; $beta6=0.; $beta8=0.; # weights in implicit time-stepping 
# $ad4=0; # old way
$upwind=0; # new way
$tp=.5; $imode=0; $adjustOmega=0; $show="trigHelmholtz.show"; 
$solver="fixedPoint";  $kx=1; $ky=1; $kz=1; $maxIterations=100; $orderInTime=-1;  # -1 = use default
$cfl=.9; $bc="d"; 
# parameter for the modulated Gaussian: beta, [x0,y0,z0], k0 
$k0=1; # parameter for the modulated Gaussian, k0<0 : k0 -> -k0/(2*pi)
$bcApproach="oneSided"; # bc Approach : cbc, lcbc, oneSided
$deflateWaveHoltz=0; $numToDeflate=1; $eigenVectorFile="eigenVectors.hdf"; 
$solverh="yale"; $maxith=2000; $rtolh=1.e-6; $atolh=1.e-5; $restart=50; $iluh=5; # parameters for direct Helmholtz solver
$solveri="yale"; $maxiti=2000; $rtoli=1.e-6; $atoli=1.e-5; # parameters for implicit time-stepping solver
GetOptions( "omega=f"=>\$omega,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"k0=f"=>\$k0,"beta=f"=>\$beta,"numPeriods=i"=>\$numPeriods,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"cfl=f"=>\$cfl,"tp=f"=>\$tp,"iMode=i"=>\$imode,\
            "solver=s"=>\$solver,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"adjustOmega=i"=>\$adjustOmega,"known=s"=>\$known,\
            "matlab=s"=>\$matlab,"maxIterations=i"=>\$maxIterations,"upwind=i"=>\$upwind,"bcApproach=s"=>\$bcApproach,\
            "degreeInSpaceForPolyPeriodic=i"=>\$degreeInSpaceForPolyPeriodic,"orderInTime=i"=>\$orderInTime,\
            "beta2=f"=>\$beta2,"beta4=f"=>\$beta4,"beta6=f"=>\$beta6,"ts=s"=>\$ts,"dtMax=f"=>\$dtMax,\
            "numberOfFrequencies=i"=>\$numberOfFrequencies,"nf=i"=>\$numberOfFrequencies,"freq=f{1,}"=>\@freq,\
            "solverh=s"=>\$solverh,"rtolh=f"=>\$rtolh,"atolh=f"=>\$atolh,"maxith=i"=>\$maxith,"restart=i"=>\$restart,"iluh=i"=>\$iluh,\
            "solveri=s"=>\$solveri,"rtoli=f"=>\$rtoli,"atoli=f"=>\$atoli,"maxiti=i"=>\$maxiti,\
            "deflateWaveHoltz=i"=>\$deflateWaveHoltz,"numToDeflate=i"=>\$numToDeflate,\
            "eigenVectorFile=s"=>\$eigenVectorFile,"show=s"=>\$show,"minStepsPerPeriod=i"=>\$minStepsPerPeriod,"go=s"=>\$go );
# 
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
if( $freq[0] eq "" ){ @freq=(1,2,3,4,5,6,7,8,9,10); }
#
if( $solveri eq "best" ){ $solveri="choose best iterative solver"; }
# 
$omega = $freq[0]; 
# pause
# -------- Start CgWaveHoltz Options ------
# debug $debug
Gaussian params $beta $x0 $y0 0 (beta,x0,y0,z0)
omega $omega
omegaSOR $omegaSOR
tol $tol 
number of frequencies $numberOfFrequencies
frequencies $freq[0] $freq[1] $freq[2] $freq[3] $freq[4] $freq[5] $freq[6] $freq[7] $freq[8]
number of periods $numPeriods
adjust omega $adjustOmega
matlab filename: $matlab
#
# choose parameters for the direct Helmholtz solver
direct solver parameters
  if( $solverh eq "best" ){ $solverh="choose best iterative solver"; }
  $solverh
  number of incomplete LU levels
    $iluh
  number of GMRES vectors
    $restart
  maximum number of iterations
    #
    $maxith
  relative tolerance
    $rtolh
  absolute tolerance
    $atolh    
  # gmres restart length
  #   30
  # pause
 #  PETSc
exit
# pause
exit
# ------ Start cgWave Options ------
interactiveMode $imode 
# time-stepping: (explicit or implicit)
$ts
if( $orderInTime > 0 ){ $cmd="orderInTime $orderInTime"; }else{ $cmd="#"; }
$cmd
cfl $cfl 
tPlot $tp 
dtMax $dtMax
implicit weights $beta2 $beta4 $beta6 $beta8
if( $known eq "polyPeriodic" ){ $setKnownOnBoundaries=1; }
set known on boundaries $setKnownOnBoundaries
# -- Here is input for cgWave 
#
$cmd="#";
if( $bcApproach eq "oneSided" ){ $cmd="useOneSidedBCs"; }
if( $bcApproach eq "cbc"      ){ $cmd="useCompatibilityBCs"; }
if( $bcApproach eq "lcbc"     ){ $cmd="useLocalCompatibilityBCs"; }
$cmd
#
bc=$bc
#
# if( $ad4>0. ){ $upwind=1; }# for backward compatibility
upwind dissipation $upwind
#
deflate WaveHoltz $deflateWaveHoltz
number to deflate $numToDeflate
eigenVectorFile $eigenVectorFile
min steps per period $minStepsPerPeriod
# pause
#
implicit solver parameters
  # NOTE: bcgs = bi-CG stab
  if( $solveri ne "yale" ){ $cmd="choose best iterative solver\n $solveri"; }else{ $cmd="choose best direct solver"; }
  $cmd
  number of incomplete LU levels
    3
  number of GMRES vectors
    20
  maximum number of iterations
    #
    $maxiti
  relative tolerance
    $rtoli
  absolute tolerance
    $atoli
 # 
  multigrid parameters
    choose good parameters: 1
    residual tolerance $rtoli
    absolute tolerance $atoli
    # debug
    #   1
    # show smoothing rates
    # Coarse level solver:  
    Oges parameters
      choose best direct solver
      # choose best iterative solver
      relative tolerance
        1.e-10
      absolute tolerance
        1.e-10
      number of incomplete LU levels
      3
    exit     
  exit
  #pause
exit
#
# helmholtzForcing
solve Helmholtz 1
#
adjust omega $adjustOmega
#   --- eigenmode-like solution 
user defined known solution...
  $pi = atan2(1.,1.)*4.;
  if( $k0 < 0 ){ $k0=-$k0/(2.*$pi); }
  if( $kx < 0 ){ $kx=-$kx/(2.*$pi); }
  if( $ky < 0 ){ $ky=-$ky/(2.*$pi); }
  if( $kz < 0 ){ $kz=-$kz/(2.*$pi); }
  if( $known eq "boxHelmholtz"){ $cmd="box Helmholtz\n $omega $kx $ky $kz"; }
  if( $known eq "computedHelmholtz"){ $cmd="computed Helmholtz\n $omega $kx $ky $kz"; }
  if( $known eq "polyPeriodic"){ $cmd="poly periodic\n $omega $degreeInSpaceForPolyPeriodic"; }
  if( $known eq "modulatedGaussian"){ $cmd="offset: $x0 $y0 $z0\n beta: $beta\n k0: $k0\n modulated Gaussian"; }
  $cmd
  # box helmholtz
  # $omega $kx $ky $kz
done
$cmd="#"; 
if( ($known eq "boxHelmholtz") || ($known eq "computedHelmholtz") ){ $cmd="helmholtzForcing\n user defined forcing...\n box Helmholtz\n exit"; }
if( $known eq "polyPeriodic" ){ $cmd="helmholtzForcing\n user defined forcing...\n poly periodic\n exit"; }
if( $known eq "modulatedGaussian" ){ $cmd="helmholtzForcing\n user defined forcing...\n modulated Gaussian\n exit"; }
$cmd
# 
exit
# --- end cgWave ---
show file $show
max iterations $maxIterations
if( $solver eq "fixedPoint" ){ $cmd="compute with fixed-point"; }else{ $cmd="#"; }
if( $solver eq "krylov" ){ $cmd="compute with petsc"; }
$cmd 
contour
exit
$cmd="#"; 
if( $go eq "go" ){ $cmd="compute with krylov\nexit"; }
if( $go eq "og" ){ $cmd="open graphics"; }
if( $go eq "direct" ){ $cmd="solve Helmholtz directly\n exit"; }
if( $go eq "d" ){ $cmd="solve Helmholtz directly\n exit"; }
if( $go eq "ds" ){ $cmd="solve Helmholtz directly\n save show file\n exit"; }
if( $go eq "krylov" ){ $cmd="compute with krylov\n exit"; }
if( $go eq "dk" ){ $cmd="solve Helmholtz directly\n zero initial condition\n compute with krylov\nexit"; }
if( $go eq "df" ){ $cmd="solve Helmholtz directly\n zero initial condition\ncompute with fixed-point\nexit"; }
if( $go eq "dfk" ){ $cmd="solve Helmholtz directly\n zero initial condition\ncompute with fixed-point\n zero initial condition\n compute with krylov\nexit"; }
$cmd 




contour
exit
exit


# Non-interactive mode -- just run and exit 
if( $imode eq 1 ){ $cmd="movie mode\n exit"; }else{ $cmd="#"; }
$cmd
