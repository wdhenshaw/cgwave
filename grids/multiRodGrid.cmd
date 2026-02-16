#
#   Multiple rods in domain 
#
#
# usage: ogen [noplot] multiRodGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i]
# 
# examples:
#     ogen -noplot multiRodGrid -interp=e -order=2 -factor=2 
#     ogen -noplot multiRodGrid -interp=e -order=2 -factor=4 
#     ogen -noplot multiRodGrid -interp=e -order=4 -factor=2 
#     ogen -noplot multiRodGrid -interp=e -order=4 -factor=4
#
# Sharper corners:
#   ogen -noplot multiRodGrid -prefix=multiRodTwoLayerGrid -configFile=rodConfig1.h -interp=e -numPerSegment=6 -order=4 -factor=8
#   ogen -noplot multiRodGrid -prefix=multiRodTwoLayerGrid -configFile=rodConfig1.h -interp=e -numPerSegment=6 -order=4 -factor=16
#   ogen -noplot multiRodGrid -prefix=multiRodTwoLayerGrid -configFile=rodConfig1.h -interp=e -numPerSegment=6 -order=4 -factor=32
#
# Rotating rod: 
#   ogen -noplot multiRodGrid -prefix=multiRodRotating -configFile=rodConfig2.h -interp=e -numPerSegment=6 -order=4 -factor=8
#   ogen -noplot multiRodGrid -prefix=multiRodRotating -configFile=rodConfig2.h -interp=e -numPerSegment=6 -order=4 -factor=16
#   ogen -noplot multiRodGrid -prefix=multiRodRotating -configFile=rodConfig2.h -interp=e -numPerSegment=6 -order=4 -factor=32
#
# Thinner: 
#   ogen -noplot multiRodGrid -prefix=multiThinRodRotating -configFile=rodConfig2.h -interp=e -numPerSegment=6 -rodWidth=0.075 -order=4 -factor=8
#   ogen -noplot multiRodGrid -prefix=multiThinRodRotating -configFile=rodConfig2.h -interp=e -numPerSegment=6 -rodWidth=0.075 -order=4 -factor=16
#   ogen -noplot multiRodGrid -prefix=multiThinRodRotating -configFile=rodConfig2.h -interp=e -numPerSegment=6 -rodWidth=0.075 -order=4 -factor=32
#
$prefix="multipleRodsGrid"; 
$configFile = "rodConfig1.h"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-1.5; $xb=1.5; $ya=-1.5; $yb=1.5; $ml=0; 
$numGhost=-1;     # if this value is set, then use this number of ghost points
$numPerSegment=2; # number of points per segment in the definition of the rod curve, increase to make sharper coners
$rodHeight=1.0; $rodWidth=.125;
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"prefix=s"=> \$prefix, "numGhost=i"=>\$numGhost,\
            "ml=i"=>\$ml,"numPerSegment=i"=>\$numPerSegment,"configFile=s"=> \$configFile,"rodWidth=f"=> \$rodWidth );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.05/$factor;
#
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
create mappings
#
$gridNames="#"; 
$count=0;        # counts grids 
#
# --- build the ROD grid 
#
#   ---- define the rod ---
include $ENV{CGWAVE}/grids/rodGrid.h
#
# $delta = $rodHeight/sqrt(2); 
#
$mapName="rodGridBase";
#
# convert to Nurbs for parallel
#
# Define a subroutine to convert a Mapping to a Nurbs Mapping
sub convertToNurbs\
{ local($old,$new,$angle,$rotationAxis,$xShift,$yShift,$zShift)=@_; \
  $cmds = "nurbs \n" . \
   "interpolate from mapping with options\n" . \
   " $old \n" . \
   " parameterize by index (uniform)\n" . \
   " number of ghost points to include\n $numGhost\n" . \
   " choose degree\n" . \
   "  3 \n" . \
   " # number of points to interpolate\n" . \
   " #  11 21 5 \n" . \
   "done\n" . \
   "rotate \n" . \
   " $angle $rotationAxis \n" . \
   "   0. 0. 0.\n" . \
   "shift\n" . \
   " $xShift $yShift $zShift\n" . \
   "mappingName\n" . \
   " $new\n" . \
   "exit"; \
}
#
$numGhost=$ng+1; # N.B. to avoid negative volumes in the ghost points interpolate ghost too in Nurbs.
#
# 
# Include the configuration file
#
include $ENV{CGWAVE}/grids/$configFile
printf("numRods=$numRods\n");
#
$cmd=""; $zShift=0; $rotationAxis=2;
#
for( $jj=0; $jj<$numRods; $jj++){ \
 $cmd .= "include $ENV{CGWAVE}/grids/buildRodGrid.h\n"; \
 }
$cmd .="#";  
# -- execute commands --
$cmd
#
#
  rectangle
    mappingName
      backGround
   set corners
    $xa $xb $ya $yb
   lines
    $nx = intmg( ($xb-$xa)/$ds +1.5 ); 
    $ny = intmg( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
    boundary conditions
      1 2 3 4
  exit
exit 
#
generate an overlapping grid
  backGround
  # innerBackGround
  $gridNames
  done
  change parameters
 # choose implicit or explicit interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
    exit
   # pause
   # open graphics
  compute overlap
exit
#
save an overlapping grid
$name
triangle
exit
