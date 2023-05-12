#
#  cgWaveHoltz script: 
#   solver = [fixedPoint|krylov] : fiexed-point or Krylov 
#   imode =1 : do not wait in cgWave
#
# old $omega=30.1; 
$numberOfFrequencies=1; 
@freq = ();  # this must be null for GetOptions to work, defaults are given below
$beta=50.; $x0=0.5; $y0=0.5; $z0=0.5; $t0=0.; $go="halt"; 
$beta=400; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; $upwind=1; $tp=.5; $imode=0; $adjustOmega=0; 
$solver="fixedPoint";  $kx=1; $ky=1; $kz=1;
$cfl=.7; $bc="d"; 
GetOptions( "x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"numPeriods=i"=>\$numPeriods,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"upwind=f"=>\$upwind,"cfl=f"=>\$cfl,"tp=f"=>\$tp,"iMode=i"=>\$imode,\
            "solver=s"=>\$solver,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"adjustOmega=i"=>\$adjustOmega,"freq=f{1,}"=>\@freq,\
            ,"go=s"=>\$go );
# 
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
if( $freq[0] eq "" ){ @freq=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15); }
# 
$omega = $freq[0];
# 
# pause
# -------- Start CgWaveHoltz Options ------
cfl $cfl 
Gaussian params $beta $x0 $y0 0 (beta,x0,y0,z0)
omega $omega
omegaSOR $omegaSOR
tol $tol 
number of periods $numPeriods
adjust omega $adjustOmega
number of frequencies $numberOfFrequencies
frequencies $freq[0] $freq[1] $freq[2] $freq[3] $freq[4] $freq[5] $freq[6] $freq[7] $freq[8] $freq[9] $freq[10] $freq[11]] $freq[12]
exit
# ------ Start cgWave Option ------
interactiveMode $imode 
tPlot $tp 
# -- Here is input for cgWave 
bc=$bc
#
upwind dissipation $upwind
## artificial dissipation $ad4
#
helmholtzForcing
solve Helmholtz 1
adjust omega $adjustOmega
#   --- eigenmode-like solution 
user defined known solution...
box Helmholtz
 $omega $kx $ky $kz
done
user defined forcing...
  box Helmholtz
exit
exit
# --- end cgWave ---
if( $solver eq "fixedPoint" ){ $cmd="compute with fixed-point"; }else{ $cmd="#"; }
if( $solver eq "krylov" ){ $cmd="compute with petsc"; }
$cmd 
contour
exit
if( $go eq "go" ){ $cmd="exit"; }else{ $cmd="#"; }
$cmd 




contour
exit
exit


# Non-interactive mode -- just run and exit 
if( $imode eq 1 ){ $cmd="movie mode\n exit"; }else{ $cmd="#"; }
$cmd
