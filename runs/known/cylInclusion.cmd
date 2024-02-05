###### small hole grid/ cylindrical inclusion, made with cylindrical scattering problem in mind #####
#
#
#some things you can pick
$factor = 0; # factor of refinement, 0,1,2...
$order = 2; # 2 or 4
$raidusi = .005; # radius of inclusion, <.1 suggested, you may have to change other variables to make overlap decent
#
#
$backgroundlines = 300*2**$factor;
$radiallines = 85*2**$factor;
$angularlines = 100*2**factor;
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
    set corners: -10,10 -10,10
    boundary conditions: 1 1 1 1
    close mapping parameters 
    exit 
  annulus 
    center: 0 0
    radii: 0.005, 1 
    mapping parameters... 
    lines: $angularlines $radiallines
    boundary conditions: -1 -1 2 0 (left,right,bot,top) 
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
  compute overlap
  change parameters
    interpolation type
      explicit for all grids
    order of accuracy
    $orderString order
    exit
  exit
save a grid
cylInclusione$factor.order$order.hdf
grid
exit
