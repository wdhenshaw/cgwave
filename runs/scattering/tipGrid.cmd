#
#  ogen -noplot tipGrid -factor=4
#  ogen -noplot tipGrid -interp=e -factor=8
#  ogen -noplot tipGrid -interp=e -factor=16
#  ogen -noplot tipGrid -interp=e -factor=32
#  ogen -noplot tipGrid -interp=e -order=4 -numGhost=3 -factor=32
# 
# THINNER: 
#  ogen -noplot tipGrid -prefix=tipGridThin -w=.02 -interp=e -factor=16
#  ogen -noplot tipGrid -prefix=tipGridThin -w=.02 -interp=e -factor=32
#
#
#  ogen -noplot tipGrid -prefix=tipGridThin -w=.02 -interp=e -order=4 -numGhost=3 -factor=32
#  ogen -noplot tipGrid -prefix=tipGridThin -w=.02 -interp=e -order=4 -numGhost=3 -factor=64
#
#             + V4       + y3
#           /   \        |  
#      V5  +     + V3    | t - y2  
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
$w=.05; # tip width 
$t=.05;
$h=.025; 
$y3 = .5;
$y2 = $y3 - .5*$t;  
$y1 = $y3 - $t;
$y0 = $y1 - $h;
# 
$prefix="tipGrid"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-1.25; $xb=1; $ya=0; $yb=1.; $ml=0;
$periodic=""; 
$numGhost=-1; # if this value is set, then use this number of ghost points
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"prefix=s"=> \$prefix,"periodic=s"=>\$periodic,\
            "numGhost=i"=>\$numGhost,"w=f"=> \$w );
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
$ds=.1/$factor;
#
create mappings
  2D Mappings...
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
     $nr=14+$ng; $nDist = ($nr-10)*$ds; 
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
    .25 40
    .10 10
    .10 10
      0 10 
    lines
      $stretchFactor=1.0; # add extra points to account for stretching
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
  # side grids
  rectangle
    mappingName
     leftSide
   set corners
    $nrs=3 + $ng; 
    #
    $xas=-$wh-$nrs*$ds; $xbs=-$wh; 
    $yas=0; $ybs=$y0 + ($ng+2)*$ds; 
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
    $yas=0; $ybs=$y0 + ($ng+2)*$ds; 
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
    0 0 3 0
  exit
exit
#
generate an overlapping grid
  backGround
  leftSide
  rightSide
  tipGrid
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
 tipGrid
exit


  exit this menu
