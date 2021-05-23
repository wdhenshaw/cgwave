#
#  cgWaveHoltz script: 
#   solver = [fixedPoint|krylov] : fiexed-point or Krylov 
#   imode =1 : do not wait in cgWave
#
$omega=30.1; $beta=50.; $x0=0.5; $y0=0.5; $z0=0.5; $t0=0.; $go="halt"; 
$beta=400; $numPeriods=1; $omegaSOR=1; $tol=1.e-3; $ad4=1; $tp=.5; $imode=0; $adjustOmega=0; 
$solver="fixedPoint";  $kx=1; $ky=1; $kz=1;
$cfl=.7; $bc="d"; 
GetOptions( "omega=f"=>\$omega,"x0=f"=>\$x0,"y0=f"=>\$y0,"z0=f"=>\$z0,"beta=f"=>\$beta,"numPeriods=i"=>\$numPeriods,\
            "omegaSOR=f"=>\$omegaSOR,"tol=f"=>\$tol,"ad4=f"=>\$ad4,"cfl=f"=>\$cfl,"tp=f"=>\$tp,"iMode=i"=>\$imode,\
            "solver=s"=>\$solver,"kx=f"=>\$kx,"ky=f"=>\$ky,"kz=f"=>\$kz,"adjustOmega=i"=>\$adjustOmega,"go=s"=>\$go );
# 
if( $bc eq "d" ){ $bc="dirichlet"; }
if( $bc eq "n" ){ $bc="neumann"; }
if( $bc eq "e" ){ $bc="evenSymmetry"; }
if( $bc eq "r" ){ $bc="radiation"; }
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
exit
# ------ Start cgWave Option ------
interactiveMode $imode 
tPlot $tp 
# -- Here is input for cgWave 
bc=$bc
#
artificial dissipation $ad4
#
helmholtzForcing
solve Helmholtz 1
adjust omega $adjustOmega
#   --- eigenmode-like solution 
user defined known solution...
box helmholtz
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
