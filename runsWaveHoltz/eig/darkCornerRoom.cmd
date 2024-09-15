#
# usage: ogen [noplot] darkCornerRoom -factor=<num> -order=[2/4/6/8] -interp=[e/i]
# 
# examples:
#     ogen -noplot darkCornerRoom -interp=e -order=2 -factor=2 
#     ogen -noplot darkCornerRoom -interp=e -order=2 -factor=4
#     ogen -noplot darkCornerRoom -interp=e -order=2 -factor=6
# 
#     ogen -noplot darkCornerRoom -interp=e -order=4 -factor=2
#     ogen -noplot darkCornerRoom -interp=e -order=4 -factor=4
# Multigrid:
#     ogen -noplot darkCornerRoom -interp=e -order=2 -factor=4 -ml=2
#     ogen -noplot darkCornerRoom -interp=e -order=2 -factor=8 -ml=2
#     ogen -noplot darkCornerRoom -interp=e -order=2 -factor=16 -ml=3
#     ogen -noplot darkCornerRoom -interp=e -order=2 -factor=32 -ml=4
#
#     ogen -noplot darkCornerRoom -interp=e -order=4 -factor=4 -ml=2
#     ogen -noplot darkCornerRoom -interp=e -order=4 -factor=8 -ml=2
#     ogen -noplot darkCornerRoom -interp=e -order=4 -factor=16 -ml=3
#     ogen -noplot darkCornerRoom -interp=e -order=4 -factor=32 -ml=4
#  ==============================================================
#   Geometry for un-illumination problem (from Daniel Appelo)
#
#
# Topology and BCs
#     bc1,bc2,bc3,bc4 are also share flag numbers
#                       | w  |
#          --------+         +---------
#         /        |bc1   bc2| h       \
#       |        +-+         +-+        | <---- foci = sqrt( b2^2 - a2^2 )
#       |        |             |        |
#       |         \-----------/         |
#       |              bc5              |
#       |                               |
#       | bc8          O                |bc6
#       |                               |
#       |                               |
#       |                               |
#       |             bc7               |
#       |         /-----------\         |
#       |        |             |        |
#       |        +-+         +-+        |
#        \         |bc4   bc3|         / 
#         ---------+         +---------
#
#
# ==================================================================
#
$prefix="darkCornerRoomGrid"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=0; $xb=8.5; $ya=-2.; $yb=2.; $ml=0; 
#
$w=.5; 
$a1=2; $b1=1.;              # ellipse 1 semi-axes
$a2=3; $b2=6.;              # ellipse 2 semi-axes
$h=$b2 - sqrt($b2*$b2-$a2*$a2); # ellipse 1 is centered at the foci of ellipse 2, foci=sqrt( b2^2-a2^2)
#
$x1=0; $y1=$b2-$h;  # center of ellipse 1
$x2=$a1; $y2=0;     # centre of ellipse 2
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"ml=i"=>\$ml,"name=s"=> \$name);
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.1/$factor;
#
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
# 
#$shape="cross";
# $ng=1; # num ghost 
# $order=2; 
$pi = 4.*atan2(1.,1.);
# $ds=1./20.;
#
#
create mappings
#  ==============================================================
# Half ellipse insert:
#
#       -a1     -w            +w       $a1 
#                |             |  h 
#        3---2---1      X      N--------      
#        |           (x1,y1)           |
#         \                           /
#          \                         /
#           ------------------------
#                ellipse : a1, b1
#
$degree=3;  # degree of Nurbs 
#
$ns=0;
$xc[$ns]=-$w;             $yc[$ns]=$y1; $ns=$ns+1;
$xc[$ns]=-$w*.5 -$a1*.5; $yc[$ns]=$y1; $ns=$ns+1;
$xc[$ns]=-$w*.25 -$a1*.75; $yc[$ns]=$y1; $ns=$ns+1;
# points on ellipse
$numEllipse=11; 
for( $i=0; $i<$numEllipse; $i++ ){ $s=($i)/($numEllipse-1); $theta=$pi*$s; \
    $xc[$ns]=$x1-$a1*cos($theta); $yc[$ns]=$y1-$b1*sin($theta); $ns=$ns+1; \
}
$xc[$ns]= $w*.25 + $a1*.75; $yc[$ns]=$y1; $ns=$ns+1;
$xc[$ns]= $w*.5 + $a1*.5; $yc[$ns]=$y1; $ns=$ns+1;
$xc[$ns]= $w;             $yc[$ns]=$y1; $ns=$ns+1;
#
$numpt=$ns; 
#  ==============================================================
#
# Create a Nurbs
#
$arcLength=0.;
$cmd="#";
for( $i=0; $i<$numpt; $i++ ){ $x=$xc[$i]; $y=$yc[$i]; \
   if( $i > 0 ){ $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 );}\
   $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; }\
$knots="#"; for( $i=$degree-1; $i<$numpt-($degree-1); $i++ ){ $s=$i/($numpt-2); $knots .= "\n $s"; } 
printf("numpt=$numpt, arcLength=$arcLength\n");
printf("knots=[$knots]\n");
printf("cmd=[$cmd]\n");
nurbs (curve)
  enter control points
    $degree
    $numpt
    $knots
    $cmd 
 parameterize by chord length
 #
 lines
  $lines=intmg($arcLength/$ds + 1.5 );
  $lines
 mappingName
   curveBoundaryInitial1
exit
# -- interpolate the initial NURBS so that we have an arclength parameterization --
nurbs (curve)
  interpolate from a mapping
    curveBoundaryInitial1
  mappingName
   curveBoundary1
exit 
# Note: buildVolume grid uses the $arcLength variable
$startCurve="curveBoundary1";
$gridName = "innerEllipse";
$direction="forward";
$bc = "1 2 5 0 0 0";
$share =  "1 2 0 0 0 0"; 
include $ENV{CGWAVE}/runsWaveHoltz/eig/buildVolumeGrid.h
#
#  ==============================================================
#  ==============================================================
#  Outer Half ellipse
#
#
#       X      N--------.      
#                        \
#                         \
#                          |
#        (0,0)     (x2,y2) | ellipse : a2, b2
#                          |
#                         /
#                        /
#          +    1---2---3
#          <-w->        (x2,-b2)
#
$degree=3;  # degree of Nurbs 
#
#
$ns=0;
$xc[$ns]=$w;              $yc[$ns]=$y2-$b2; $ns=$ns+1;
$xc[$ns]=$w*.5 + $x2*.5;  $yc[$ns]=$y2-$b2; $ns=$ns+1;
$xc[$ns]=$w*.25+ $x2*.75; $yc[$ns]=$y2-$b2; $ns=$ns+1;
# points on ellipse
$numEllipse=31; 
for( $i=0; $i<$numEllipse; $i++ ){ $s=($i)/($numEllipse-1); $theta=$pi*$s; \
    $xc[$ns]=$x2+$a2*sin($theta); $yc[$ns]=$y2-$b2*cos($theta); $ns=$ns+1; \
}
$xc[$ns]=$w*.25 +$x2*.75; $yc[$ns]=$y2+$b2; $ns=$ns+1;
$xc[$ns]=$w*.5  +$x2*.5;  $yc[$ns]=$y2+$b2; $ns=$ns+1;
$xc[$ns]=$w;              $yc[$ns]=$y2+$b2; $ns=$ns+1;
#
$numpt=$ns; 
#
#
# Create a Nurbs
#
$arcLength=0.;
$cmd="#";
for( $i=0; $i<$numpt; $i++ ){ $x=$xc[$i]; $y=$yc[$i]; \
   if( $i > 0 ){ $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 );}\
   $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; }\
$knots="#"; for( $i=$degree-1; $i<$numpt-($degree-1); $i++ ){ $s=$i/($numpt-2); $knots .= "\n $s"; } 
printf("numpt=$numpt, arcLength=$arcLength\n");
printf("knots=[$knots]\n");
printf("cmd=[$cmd]\n");
nurbs (curve)
  enter control points
    $degree
    $numpt
    $knots
    $cmd 
 parameterize by chord length
 #
 lines
  $lines=intmg($arcLength/$ds + 1.5 );
  $lines
 mappingName
   curveBoundaryInitial2
exit
# -- interpolate the initial NURBS so that we have an arclength parameterization --
nurbs (curve)
  interpolate from a mapping
    curveBoundaryInitial2
  mappingName
   curveBoundary2
exit 
#
#
$startCurve="curveBoundary2";
$gridName = "outerEllipse";
$direction="backward";
$bc    = "3 2 6 0 0 0";
$share = "3 2 0 0 0 0"; 
include $ENV{CGWAVE}/runsWaveHoltz/eig/buildVolumeGrid.h
#
# Define a subroutine to convert a Mapping to a Nurbs Mapping
sub convertToNurbs\
{ local($old,$new,$xScale,$yScale,$xShift,$yShift,$bc,$share)=@_; \
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
   "scale \n" . \
   " $xScale $yScale \n" . \
   "shift\n" . \
   " $xShift $yShift $zShift\n" . \
   "boundary conditions\n" . \
   "  $bc\n" . \
   "share\n" . \
   "  $share\n" . \
      "mappingName\n" . \
   " $new\n" . \
   "exit"; \
}
#
$numGhost=$ng+1; # N.B. to avoid negative volumes in the ghost points interpolate ghost too in Nurbs.
#
$xScale=1; $yScale=1;
$xShift=0; $yShift=0; $zShift=0; 
$bc    = "1 2 5 0 0 0";
$share = "1 2 0 0 0 0";
convertToNurbs(innerEllipse,innerEllipse1,$xScale,$yScale,$xShift,$yShift,$bc,$share);
$cmds
#
$xScale=1; $yScale=-1;
$xShift=0; $yShift=0; $zShift=0; 
$bc = "4 3 7 0 0 0";
$share =  "4 3 0 0 0 0";
convertToNurbs(innerEllipse,innerEllipse2,$xScale,$yScale,$xShift,$yShift,$bc,$share);
$cmds
#
#
$xScale=1; $yScale=1;
$bc    = "3 2 6 0 0 0";
$share = "3 2 0 0 0 0"; 
convertToNurbs(outerEllipse,outerEllipse1,$xScale,$yScale,$xShift,$yShift,$bc,$share);
$cmds
# 
$xScale=-1; $yScale=1;
$bc    = "4 1 8 0 0 0";
$share = "4 1 0 0 0 0"; 
convertToNurbs(outerEllipse,outerEllipse2,$xScale,$yScale,$xShift,$yShift,$bc,$share);
$cmds
#
#
#  Background 
#
  $xa=-($x2+$a2); $xb=-$xa; 
  $ya=-($y2+$b2); $yb=-$ya; 
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
      0 0 0 0 
  exit
#
# End grids
#
$capWidth = .5; # width of cap grid 
#
# ----- top right
  $xa=$w;  $xb=$xa+ $capWidth;
  $ya=$y1; $yb=$y2+$b2; 
  rectangle
    mappingName
      topRightCap
   set corners
    $xa $xb $ya $yb
   lines
    $nx = intmg( ($xb-$xa)/$ds +1.5 ); 
    $ny = intmg( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
    boundary conditions
      2 0 0 0 
    share
      2 0 0 0
  exit
#
#  top left 
  $xa=-$w-$capWidth;  $xb=-$w;
  $ya=$y1; $yb=$y2+$b2; 
  rectangle
    mappingName
      topLeftCap
   set corners
    $xa $xb $ya $yb
   lines
    $nx = intmg( ($xb-$xa)/$ds +1.5 ); 
    $ny = intmg( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
    boundary conditions
      0 1 0 0 
    share
      0 1 0 0
  exit  
# -- bottom right
  $xa=$w;  $xb=$xa+ $capWidth;
  $yb=-$y1; $ya=-($y2+$b2); 
  rectangle
    mappingName
      bottomRightCap
   set corners
    $xa $xb $ya $yb
   lines
    $nx = intmg( ($xb-$xa)/$ds +1.5 ); 
    $ny = intmg( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
    boundary conditions
      3 0 0 0 
    share
      3 0 0 0
  exit
#
# --- bottom left 
  $xa=-$w-$capWidth;  $xb=-$w;
  $yb=-$y1; $ya=-($y2+$b2); 
  rectangle
    mappingName
      bottomLeftCap
   set corners
    $xa $xb $ya $yb
   lines
    $nx = intmg( ($xb-$xa)/$ds +1.5 ); 
    $ny = intmg( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
    boundary conditions
      0 4 0 0 
    share
      0 4 0 0
  exit  
  exit
generate an overlapping grid
  backGround
  topRightCap
  topLeftCap
  bottomRightCap
  bottomLeftCap
  innerEllipse1
  outerEllipse1
  innerEllipse2
  outerEllipse2
  done choosing mappings
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
darkCornerRoom
exit    
