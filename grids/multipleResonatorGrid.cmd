#
#   Spit Ring resonator 
#
# usage: ogen [noplot] multipleResonatorGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i]
# 
# examples:
#     ogen -noplot multipleResonatorGrid -angle=90 -interp=e -order=2 -factor=4 
#     ogen -noplot multipleResonatorGrid -angle=90 -interp=e -order=2 -factor=8
# 
#     ogen -noplot multipleResonatorGrid -interp=e -order=4 -factor=8
#
#   ogen -noplot multipleResonatorGrid -prefix=splitRingGridAngle90 -angle=90 -interp=e -order=2 -factor=4 
#
$prefix="multipleResonatorGrid"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-1.25; $xb=2.5; $ya=-1.25; $yb=1.25; $ml=0; 
$addBackGround=0; 
$angle=0; # rotation angle
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,"angle=f"=> \$angle,\
            "interp=s"=> \$interp,"name=s"=> \$name,"prefix=s"=> \$prefix, "addBackGround=i"=>\$addBackGround );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
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
$xr=0.0; $yr=0.0; # center for rotation
$xShift=0.0; $yShift=0; $xScale=1; $yScale=1;
$count=$count+1; $gridName="splitRing$count"; $gridNames .= "\n" . $gridName;
include $ENV{CGWAVE}/grids/transform.h
#
# 
$xr=0.0; $yr=0; # center for rotation
$xShift=1.35; $yShift=0; $xScale=1.2; $yScale=1.2;
$xShift=1.35; $yShift=0; $xScale=1.0; $yScale=1.0;
$count=$count+1; $gridName="splitRing$count"; $gridNames .= "\n" . $gridName;
include $ENV{CGWAVE}/grids/transform.h
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
multipleResonatorGrid
exit
