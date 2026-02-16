# ----- ROD GRID 
#    Start on a flat edge since the BC is clamped. 
#    NOTE: go in a counter clockwise direction
#
# 
#                 P3     P2    
#             ^    +-----+    
#             |    |     |    
#             |    |     |    
#             |    |     |      P3      CH = rodHeight
#             |    |     |              CW = rodWidth 
#          CH |  P4|  C  |P1, P7    
#             |    |     |           
#             |    |     |       
#             |    |     |    
#             \/   |     |    
#               P5 +-----+  P6  
#                ->  CW   <- 
#                                   
#
# ---These next parameters may be optionally set --
if( $rodHeight eq "" ){ $rodHeight=1.0; }
if( $rodWidth eq "" ){ $rodWidth=0.125; }
if( $rodBC eq "" ){ $rodBC=7; }
if( $numPerSegment eq "" ){ $numPerSegment=2; } # number of control points per segment -- increase to make corners sharper
#
$degree=3;  # degree of Nurbs 
$cx=0.; $cy=0.;  # center for the rod
#
$nc=7;  # number of control points 
#
    $x1 = $cx +.5*$rodWidth; $y1=$cy;                   # center-right point on I-beam
    $x2 = $x1;               $y2=$y1+.5*$rodHeight;
    $x3 = $cx -.5*$rodWidth; $y3=$y2;
    $x4 = $x3;               $y4=$y1; 
    $x5 = $x3;               $y5=$cy-.5*$rodHeight;
    $x6 = $x1;               $y6=$y5;
    $x7 = $x1;               $y7=$y1; 
#
@xc=($x1,$x2,$x3,$x4,$x5,$x6,$x7);
@yc=($y1,$y2,$y3,$y4,$y5,$y6,$y7);
$radX=0; $radY=0.; # for inner background grid  **FIX ME**
$cmd="#";  $ns=0;  $arcLength=0.; $x0=$xc[0]; $y0=$yc[0];
for( $ic=0; $ic<$nc-1; $ic++ ){ $numpt=$numPerSegment; if( $ic == $nc-2 ){ $numpt=$numPerSegment+1;} \
for( $i=0; $i<$numpt; $i++ ){ $s=($i)/($numPerSegment); $ns=$ns+1;  \
   $x=(1.-$s)*$xc[$ic]+$s*$xc[$ic+1]; \
   $y=(1.-$s)*$yc[$ic]+$s*$yc[$ic+1];   \
   $radX = max($radX,abs($x)); $radY = max($radY,abs($y));  \
   $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 ); $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; } }\
  $knots="#"; 
  for( $i=$degree-1; $i<$ns-($degree-1); $i++ ){ $s=$i/($ns-2); $knots .= "\n $s"; } \
#
nurbs (curve)
  periodicity
   2
  enter control points
    $degree
    $ns
    $knots
    $cmd 
 parameterize by chord length
 #
 lines
  $lines=intmg($arcLength/$ds + 1.5 );
  $lines
 mappingName
   rodBoundaryInitial
exit
# -- interpolate the initial NURBS so that we have an arclength parameterization --
nurbs (curve)
  interpolate from a mapping
    rodBoundaryInitial
  mappingName
   rodBoundary
exit 
# 
# -- Make a hyperbolic grid --
#
  $nr = intmg( 6 + $order/2 );
  hyperbolic
    forward
    $nDist=($nr-2)*$ds;
    distance to march $nDist
    $nrm=$nr-1; 
    lines to march $nrm
    $nTheta = int($arcLength/$ds+1.5);
    points on initial curve $nTheta
    uniform dissipation 0.05
    volume smooths $numberOfVolumeSmooths
    equidistribution 0 (in [0,1])
    #
    # spacing: geometric
    # geometric stretch factor 1.05 
    #
    generate
    # open graphics
    boundary conditions
      -1 -1 $rodBC 0 0 0
    share 
       0 0 100 0 0 0
    name rodGridBase
  exit

