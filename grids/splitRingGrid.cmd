#
#   Spit Ring resonator 
#
# usage: ogen [noplot] splitRingGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i]
# 
# examples:
#
# Gap on top: 
#     ogen -noplot splitRingGrid -prefix=splitRingGridTopGap -angle=90 -interp=e -order=2 -factor=4 
#     ogen -noplot splitRingGrid -prefix=splitRingGridTopGap -angle=90 -interp=e -order=2 -factor=8
# Order=4:
#     ogen -noplot splitRingGrid -prefix=splitRingGridTopGap -angle=90 -interp=e -order=4 -factor=4
#     ogen -noplot splitRingGrid -prefix=splitRingGridTopGap -angle=90 -interp=e -order=4 -factor=8
#
# Gap on left: 
#     ogen -noplot splitRingGrid -prefix=splitRingGridLeftGap -angle=180 -interp=e -order=2 -factor=4 
#     ogen -noplot splitRingGrid -prefix=splitRingGridLeftGap -angle=180 -interp=e -order=2 -factor=8
#     ogen -noplot splitRingGrid -prefix=splitRingGridLeftGap -angle=180 -interp=e -order=4 -factor=4 
#     ogen -noplot splitRingGrid -prefix=splitRingGridLeftGap -angle=180 -interp=e -order=4 -factor=8
#
#     ogen -noplot splitRingGrid -prefix=splitRingGridLeftGap -angle=180 -interp=e -order=4 -numGhost=3 -factor=4 
#     ogen -noplot splitRingGrid -prefix=splitRingGridLeftGap -angle=180 -interp=e -order=4 -numGhost=3 -factor=8
#     ogen -noplot splitRingGrid -prefix=splitRingGridLeftGap -angle=180 -interp=e -order=4 -numGhost=3 -factor=16
#     ogen -noplot splitRingGrid -prefix=splitRingGridLeftGap -angle=180 -interp=e -order=4 -numGhost=3 -factor=32
#     ogen -noplot splitRingGrid -prefix=splitRingGridLeftGap -angle=180 -interp=e -order=4 -numGhost=3 -factor=64
$prefix="splitRingGrid"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-1.25; $xb=1.25; $ya=-1.25; $yb=1.25; $ml=0; 
$addBackGround=0; 
$angle=0; # rotation angle
$numGhost=-1;  # if this value is set, then use this number of ghost points
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,"angle=f"=> \$angle,\
            "interp=s"=> \$interp,"name=s"=> \$name,"prefix=s"=> \$prefix, "addBackGround=i"=>\$addBackGround,"numGhost=i"=>\$numGhost);
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
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
#  Coarse background 
#
  rectangle
    mappingName
      backGround
   set corners
    $xa $xb $ya $yb
   lines
    $nx = int( ($xb-$xa)/$ds + 1.5 ); 
    $ny = int( ($yb-$ya)/$ds + 1.5 ); 
    $nx $ny
  exit
#
$gridNames="#"; 
#
$count=0;
include $ENV{CGWAVE}/grids/splitRing.h
$mapName="splitRingGridBase";
#
# convert to Nurbs for parallel
#
# Define a subroutine to convert a Mapping to a Nurbs Mapping
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
   " 0. 0. 0.\n" . \
   "shift\n" . \
   " $xShift $yShift $zShift\n" . \
   "mappingName\n" . \
   " $new\n" . \
   "exit"; \
}
#
$numGhost=$ng+1; # N.B. to avoid negative volumes in the ghost points interpolate ghost too in Nurbs.
#
# $xScale=1; $yScale=1;
$xShift=0; $yShift=0; $zShift=0; $rotationAxis=2; 
# $bc = "-1 -1 $splitRingBC 0 0 0";
# $share = "0 0 100 0 0 0";
$count=$count+1; $gridName="splitRing$count"; $gridNames .= "\n" . $gridName;
convertToNurbs($mapName,$gridName,$angle,$rotationAxis,$xShift,$yShift,$zShift);
$cmds
#
# Old: 
# 
# $xr=.0; $yr=0; # center for rotation
# $xShift=0.0; $yShift=0; $xScale=1; $yScale=1; 
# $count=$count+1; $gridName="splitRing$count"; $gridNames .= "\n" . $gridName;
# include $ENV{CGWAVE}/grids/transform.h
#
exit 
#
generate an overlapping grid
  # if( $addBackGround eq 1 ){ $cmd="backGround\n innerBackGround"; }else{ $cmd="innerBackGround"; }
  # $cmd 
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
splitRingGrid
exit
