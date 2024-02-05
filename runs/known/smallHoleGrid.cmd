###### small hole grid, made with scattering problem in mind #####
#
# ogen -noplot smallHoleGrid.cmd -order=2 -numGhost=2
# ogen -noplot smallHoleGrid.cmd -order=4 -numGhost=3
#
# Initial version from Alli Carson
#
#some things you can pick
$factor = 0; # factor of refinement, 0,1,2...
$order = 4; # 2 or 4
$numGhost=2;
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"ng=f"=> \$numGhost );
# 
if( $order eq 4 ){$numGhost=3;}
#
#
$backgroundlines = 150*2**$factor;
$radiallines = 85*2**$factor;
$angularlines = 90*2**factor;
#
#
if( $ft eq -1 ){ $ft=$fx; }
if( $order eq 2 ){$orderString="second";}elsif($order eq 4){$orderString="fourth";}
#
#
# now the ogen code
create mappings 
  2D Mappings... 
  rectangle 
    mapping parameters... 
    name: backgroundGrid 
    lines: $backgroundlines $backgroundlines
    set corners: -5,5 -5,5
    boundary conditions: 1 1 1 1
    close mapping parameters 
    exit 
  annulus 
    center: 0 0
    radii: 0.005, 1 
    mapping parameters... 
    lines: $angularlines $radiallines
    boundary conditions: -1 -1 5 0 (left,right,bot,top) 
    close mapping parameters 
    exit
  transform Mappings... 
  stretch coordinates 
    Stretch r2:exp 
    STP:stretch r2 exp: parameters 0 1 1 5.5 0.1 (a0,ar,a,b,c)
    close r2 stretching parameters
    exit
  exit
generate an overlapping grid
  backgroundGrid
  stretched-Annulus
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
    compute overlap
  exit
save a grid
smallHoleGride$factor.order$order.hdf
grid
exit
