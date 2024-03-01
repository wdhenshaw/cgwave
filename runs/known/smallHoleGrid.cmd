###### small hole grid, made with scattering problem in mind #####
#
# ogen -noplot smallHoleGrid.cmd -order=2 -rad=.1 -factor=2 
# ogen -noplot smallHoleGrid.cmd -order=2 -rad=.05 -factor=2
# ogen -noplot smallHoleGrid.cmd -order=2 -rad=.01 -factor=2
#
# Use these:
# ogen -noplot smallHoleGrid.cmd -order=2 -numGhost=2 -rad=.01 -factor=4
# ogen -noplot smallHoleGrid.cmd -order=4 -numGhost=3 -rad=.01 -factor=4
#
# FIXED radial distance: use these for convergence studies:
#  ogen -noplot smallHoleGrid.cmd -prefix=holeGridRad0p1 -order=2 -rgd=fixed -rad=.1 -factor=2
#
# Initial version from Alli Carson
#
#some things you can pick
$prefix = "smallHoleGrid"; 
$rgd="var"; # set to "fixed" for fixed radial grid
$factor = 1; # factor of refinement, 1,2...
$order = 2; # 2 or 4
$numGhost=-1;  # set to >0 to explicitly set the number of ghost lines
$ml=0;
$xa=-2; $xb=2.; $ya=-2; $yb=2;
$rad=.1; # inner radius
$deltaRad = .25; # for fixed radial distance
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"rad=f"=> \$rad,"numGhost=i"=> \$numGhost,"prefix=s"=> \$prefix,\
            "rgd=s"=> \$rgd,"deltaRad=f"=> \$deltaRad );
# 
if( $numGhost eq -1  ){ $numGhost = $order/2 + 1; } # add one for upwinding
#
$ds=.1/$factor;
$pi = 4.*atan2(1.,1.);
#
#
if( $order eq 2 ){$orderString="second";}elsif($order eq 4){$orderString="fourth";}
#
#
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
# now the ogen code
create mappings 
  2D Mappings... 
  rectangle 
    mapping parameters... 
    name: backgroundGrid 
    set corners: $xa $xb $ya $yb
    lines
      $nx = intmg( ($xb-$xa)/$ds +1.5 ); 
      $ny = intmg( ($yb-$ya)/$ds +1.5 ); 
      $nx $ny    
    boundary conditions: 1 2 3 4 
    exit 
#
  # this nr only chooses the radial size of the hole grid , stretching will choose a new $dr
  $nr = intmg( 8+$numGhost ); 
  $outerRad = $rad + $nr*$ds;
  if( $rgd eq "fixed" ){ $outerRad = $rad + $deltaRad; $nr=intmg( $deltaRad/$ds + 1.5 ); }
  annulus 
    center: 0 0
    radii: $rad, $outerRad 
    lines 
      $nTheta= intmg(2*$pi*$outerRad/$ds + 1 );
      $nTheta $nr
    # open graphics
    boundary conditions
       -1 -1 5 0
    mappingName
      unstretchedAnnulus
    close mapping parameters 
    exit
# 
  stretch coordinates 
    STRT:stretch resolution factor 1.0
    Stretch r2:exp to linear
    $dTheta = 2*$pi*$rad/$nTheta; # theta spacing on inner radius
    $drMax=$ds; 
    $drMin = $dTheta;
    STP:stretch r2 expl: min dx, max dx $drMin $drMax
    mappingName
       annulus
    # Stretch r2:exp 
    # STP:stretch r2 exp: parameters 0 1 1 5.5 0.1 (a0,ar,a,b,c)
    # close r2 stretching parameters
    exit
  exit
generate an overlapping grid
  backgroundGrid
  annulus
  done
  change parameters
    interpolation type
      explicit for all grids
    order of accuracy
      $orderString order
    ghost points
      all
      $numGhost $numGhost $numGhost $numGhost
    exit
    # open graphics
    compute overlap
  exit
save a grid
if( $rgd eq "fixed" ){ $prefix .= "Fixed"; }
$gridFileName = $prefix . "e$factor.order$order.hdf";
$gridFileName
grid
exit
