#
#  Grid for mulitple knife edges (started from tipGrid.cmd)
#
#    Specify any number of upward and downward knife edges
#
#   ogen -noplot multiKnivesGrid -prefix=eightKnivesGrid -nke=4 -xke=-0.25 -0.15 0.15 0.25 -nker=4 -xker=-0.2 -0.1 0.1 0.2 -interp=e -order=2 -factor=16
#   ogen -noplot multiKnivesGrid -prefix=eightKnivesGrid -nke=4 -xke=-0.35 -0.275 0.275 0.35 -nker=4 -xker=-0.2 -0.125 0.125 0.2 -interp=e -order=2 -factor=32
# 
#   ogen -noplot multiKnivesGrid -prefix=twelveKnivesGrid -nke=6 -xke=-0.35 -0.25 -0.15 0.15 0.25 0.35 -nker=6 -xker=-0.3 -0.2 -0.1 0.1 0.2 0.3 -interp=e -order=2 -factor=16
#   ogen -noplot multiKnivesGrid -prefix=twelveKnivesGrid -nke=6 -xke=-0.35 -0.25 -0.15 0.15 0.25 0.35 -nker=6 -xker=-0.3 -0.2 -0.1 0.1 0.2 0.3 -interp=e -order=4 -factor=16
#
# SIXTEEN: 
# ogen -noplot multiKnivesGrid -prefix=sixteenKnivesGrid -nke=8 -xke=-.45 -0.35 -0.25 -0.15 0.15 0.25 0.35 .45 -nker=8 -xker=-0.4 -0.3 -0.2 -0.1 0.1 0.2 0.3 0.4 -interp=e -order=2 -factor=256
# ogen -noplot multiKnivesGrid -prefix=sixteenKnivesGrid -nke=8 -xke=-.45 -0.35 -0.25 -0.15 0.15 0.25 0.35 .45 -nker=8 -xker=-0.4 -0.3 -0.2 -0.1 0.1 0.2 0.3 0.4 -interp=e -order=2 -factor=16
# ogen -noplot multiKnivesGrid -prefix=sixteenKnivesGrid -nke=8 -xke=-.45 -0.35 -0.25 -0.15 0.15 0.25 0.35 .45 -nker=8 -xker=-0.4 -0.3 -0.2 -0.1 0.1 0.2 0.3 0.4 -interp=e -order=4 -factor=32
#
# GEOMETRY OF THE TIP OF KNIFE: 
#
#
#             + V4       + y3
#           /   \        |  
#      V5  +     + V3 y2 | t   
#         /       \      |
#     V6 +         + V2  - y1 
#        |         |     | |
#        |         |     | |
#        |         |     | h
#        |         |     | | 
#        |         |     | |
#       V7         V1    | y0 
#        <--- w --->
#
$w=.02; # tip width 
$t=.05;
$h=.025; 
$y3 = .5;
$y2 = $y3 - .5*$t;  
$y1 = $y3 - $t;
$y0 = $y1 - $h;
# 
$prefix="multiKnivesGrid"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-0.5; $xb=0.5; $ya=0; $yb=0.55; $ml=0;
$periodic=""; 
$numGhost=-1; # if this value is set, then use this number of ghost points
$nke=1;       # number of upward knife edges
$nker=1;      # number of downward (reversed) knife edges
@xke = (); @xker=();  # these must be null for GetOptions to work, defaults are given below
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"prefix=s"=> \$prefix,"periodic=s"=>\$periodic,\
            "numGhost=i"=>\$numGhost,"w=f"=> \$w,"nke=i"=>\$nke,"nker=i"=>\$nker,"xke=f{1,}"=>\@xke,"xker=f{1,}"=>\@xker );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix=""; 
if( $periodic eq "p" ){ $suffix = "p"; }
if( $periodic eq "np" ){ $suffix = "np"; }
if( $periodic eq "pn" ){ $suffix = "pn"; }
$suffix .= ".order$order";
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
#
if( $xke[0]  eq "" ){ @xke=(0,0,0,0); }
if( $xker[0] eq "" ){ @xker=(.2,-.2,0,0); }
# 
$ds=.1/$factor;
#
create mappings
  2D Mappings...
  #  ------------------------------------------------------
  #  ------------- MASTER UPWARD KNIFE TIP ----------------
  #  ------------------------------------------------------
  smoothedPolygon
    $wh=$w*.5;
    $wq=$w*.25;
    $v3y = $h+$t*.5; 
    $v5y = $h+$t*.5; 
    $v4y = $h+$t; 
    vertices
      7
      $wh  $y0
      $wh  $y1
      $wq  $y2
        0  $y3
      -$wq $y2
      -$wh $y1
      -$wh $y0
    # 
    # $nr=20; $nDist = ($nr-12)*$ds; 
    # $nr=14+$ng; $nDist = ($nr-10)*$ds; 
     $nr=18+$ng; $nDist = ($nr-10)*$ds; 
    n-dist
      fixed normal distance
      -$nDist
    n-stretch
      1 20
    sharpness
      15
      15
      15
      15
      15
      15
      15      
    t-stretch
      0 10
    .10 10
    .10 10
    .30 30
    .10 10
    .10 10
      0 10 
    lines
      # $stretchFactor=1.0; # add extra points to account for stretching
      $stretchFactor=1.75; # add extra points to account for stretching
      $arcLength = $stretchFactor*( $y3-$y0 ) + .5*$w; 
      $ns = int( 2*($arcLength)/$ds + 1.5 ); 
      $ns $nr       
    mappingName
      tipGrid
    boundary conditions
      0 0 5 0 
    share
      0 0 5 0 
    correct corners
    corners
    specify positions of corners
     $xc=$wh+$nDist; 
     $wh  $y0
    -$wh  $y0
     $xc  $y0
    -$xc  $y0    
    # open graphics
  exit
  #
  # side grids
  #
  rectangle
    mappingName
     leftSide
   set corners
    $nrs=3 + $ng; 
    #
    $xas=-$wh-$nrs*$ds; $xbs=-$wh; 
    $yas=0; $ybs=$y0 + ($ng+$order)*$ds; 
    $xas $xbs $yas $ybs
   lines
    $nx = int( ($xbs-$xas)/$ds + 1.5 ); 
    $ny = int( ($ybs-$yas)/$ds + 1.5 ); 
    $nx $ny
   boundary conditions
     0 5 3 0 
   share
    0 5 3 0
  exit 
  #
  rectangle
    mappingName
     rightSide
   set corners
    $xas=$wh; $xbs=$wh+$nrs*$ds; 
    $yas=0; $ybs=$y0 + ($ng+$order)*$ds; 
    $xas $xbs $yas $ybs
   lines
    $nx = int( ($xbs-$xas)/$ds +1.5 ); 
    $ny = int( ($ybs-$yas)/$ds +1.5 ); 
    $nx $ny
   boundary conditions
     5 0 3 0 
   share
    5 0 3 0 
  exit     
#
# open graphics
  # Background
  rectangle
    mappingName
      backGround
   set corners
    $xa $xb $ya $yb
   lines
    $nx = int( ($xb-$xa)/$ds +1.5 ); 
    $ny = int( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
   boundary conditions
    1 2 3 4
    share
    0 0 3 4
  exit
  # Upside down knife: 
  #  --------------------------------------------------------
  #  ------------- MASTER DOWNWARD KNIFE TIP ----------------
  #  --------------------------------------------------------
  $xShift=0;
  #
  $wr=$w; # tip width 
  $tr=.05;
  $hr=.025; 
    $whr=$wr*.5;
    $wqr=$wr*.25;  
  $y3r = .1;
  $y2r = $y3r + .5*$tr;  
  $y1r = $y3r + $tr;
  $y0r = $y1r + $hr;  
  $v1x = -$whr + $xShift;
  $v2x = -$whr + $xShift;
  $v3x = -$wqr + $xShift;
  $v4x =    0  + $xShift;
  $v5x =  $wqr + $xShift;
  $v6x =  $whr + $xShift;
  $v7x =  $whr + $xShift;
  smoothedPolygon
    # $v3y = $h+$t*.5; 
    # $v5y = $h+$t*.5; 
    # $v4y = $h+$t; 
    vertices
      7
      $v1x  $y0r
      $v2x  $y1r
      $v3x  $y2r
      $v4x  $y3r
      $v5x  $y2r
      $v6x  $y1r
      $v7x  $y0r
    # 
    # $nr=20; $nDist = ($nr-12)*$ds; 
    #  $nr=14+$ng; $nDist = ($nr-10)*$ds; 
    n-dist
      fixed normal distance
      -$nDist
    n-stretch
      1 20
    sharpness
      15
      15
      15
      15
      15
      15
      15      
    t-stretch
      0 10
    .10 10
    .10 10
    .30 30
    .10 10
    .10 10
      0 10 
    lines
      # $stretchFactor=1.0; # add extra points to account for stretching
      # $arcLength = $stretchFactor*( $y3r-$y0r ) + .5*$w; 
      # $ns = int( 2*($arcLength)/$ds + 1.5 ); 
      $ns $nr       
    mappingName
      tipGridReversed
    boundary conditions
      0 0 5 0 
    share
      0 0 6 0 
    correct corners
    corners
    specify positions of corners
     $cx1=$v1x-$nDist; 
     $cx7=$v7x+$nDist; 
     $v1x $y0r 
     $v7x $y0r
     $cx1 $y0r
     $cx7 $y0r
     # open graphics
  exit
  #
  # side grids (reversed)
  #
  rectangle
    mappingName
     leftSideReversed
   set corners
    $nrs=3 + $ng; 
    #
    $xas=$v1x-$nrs*$ds; $xbs=$v1x; 
    $yas=$y0r - ($ng+$order)*$ds;  $ybs=$yb;
    $xas $xbs $yas $ybs
   lines
    $nx = int( ($xbs-$xas)/$ds + 1.5 ); 
    $ny = int( ($ybs-$yas)/$ds + 1.5 ); 
    $nx $ny
   boundary conditions
     0 5 0 4  
   share
     0 6 0 4 
  exit 
  #
  rectangle
    mappingName
     rightSideReversed
   set corners
    $xas=$v7x; $xbs=$v7x+$nrs*$ds; 
    $yas=$y0r - ($ng+$order)*$ds; $ybs=$yb; 
    $xas $xbs $yas $ybs
   lines
    $nx = int( ($xbs-$xas)/$ds +1.5 ); 
    $ny = int( ($ybs-$yas)/$ds +1.5 ); 
    $nx $ny
   boundary conditions
     5 0 0 4  
   share
     6 0 0 4 
  exit    
#
#   CREATE ACTUAL UPWARD AND DOWNWARD KNIFE EDGES BY TRANSLATIONS
$knifeMappingNames="";
$count=0; 
#
# ======================================================
# Define a function to build an inner-annulus, outer-annulus and inner square.
# usage:
#   makeDisk(radius,xCenter,yCenter)
# =====================================================
sub makeKnifeEdge \
{ local($knifeName,$leftSideName,$rightSideName,$xShift)=@_; \
  $count = $count + 1; \
  $knifeMappingNames = $knifeMappingNames . "   $knifeName$count\n". "   $leftSideName$count\n" . "   $rightSideName$count\n"; \
  $commands = \
    "rotate/scale/shift\n" . \
    "  transform which mapping?\n" . \
    "  $knifeName\n" . \
    "  shift\n" . \
    "    $xShift 0\n" . \
    "  mappingName\n" . \
    "  $knifeName$count\n" . \
    "  exit\n" . \
    "rotate/scale/shift\n" . \
    "  transform which mapping?\n" . \
    "  $leftSideName\n" . \
    "  shift\n" . \
    "    $xShift 0\n" . \
    "  mappingName\n" . \
    "  $leftSideName$count\n" . \
    "  exit\n" . \
    "rotate/scale/shift\n" . \
    "  transform which mapping?\n" . \
    "  $rightSideName\n" . \
    "  shift\n" . \
    "    $xShift 0\n" . \
    "  mappingName\n" . \
    "  $rightSideName$count\n" . \
    "  exit\n"; \
}
#
sub makeKnives \
{ local($nke,$xke,$nker,$xker)=@_; \
  local $cmds; $cmds=""; \
  for( $i=0; $i<$nke; $i++ ){ \
    makeKnifeEdge("tipGrid","leftSide","rightSide",$xke[$i]); \
    $cmds = $cmds . $commands; \
  }\
  for( $i=0; $i<$nker; $i++ ){ \
    makeKnifeEdge("tipGridReversed","leftSideReversed","rightSideReversed",$xker[$i]); \
    $cmds = $cmds . $commands; \
  }\
  $commands=$cmds; \
}
makeKnives($nke,$xke,$nker,$xker);
$commands
# makeKnifeEdge("tipGrid","leftSide","rightSide",0.);
# $commands
# makeKnifeEdge("tipGridReversed","leftSideReversed","rightSideReversed",.25);
# $commands
#
#
#
#  $xshiftr2=-.5; 
#  rotate/scale/shift
#    transform which mapping?
#    tipGridReversed
#    shift
#      $xshiftr2 0
#    mappingName
#    tipGridReversed2
#    exit
#  rotate/scale/shift
#    transform which mapping?
#    leftSideReversed
#    shift
#      $xshiftr2 0
#    mappingName
#    leftSideReversed2
#    exit
#  rotate/scale/shift
#    transform which mapping?
#    rightSideReversed
#    shift
#      $xshiftr2 0
#    mappingName
#    rightSideReversed2
#    exit
# open graphics
#
#   view Mappings...
#     tipGrid
#     tipGridReversed
#     leftSideReversed
#     rightSideReverse
#     backGround    
# open graphics
exit
#
generate an overlapping grid
  backGround
  $knifeMappingNames
#  leftSide
#  rightSide
#  tipGrid
#  # 
#  leftSideReversed
#  rightSideReversed
#  tipGridReversed 
#  # 
#  leftSideReversed2
#  rightSideReversed2
#  tipGridReversed2     
  done
# 
  change parameters
    ###$solidDomains
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
 multiKnivesGrid
exit


  exit this menu
