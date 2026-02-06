#
#  ----- Build a split annulus resonator ------
#
#     
#     ra = inner radius of annulus
#     rb = outer radius 
#     theta1 = start curve here on rb (go counter clockwise)
#     phi = half open angle of split 
#     theta2 = angle of bottom split
#
#
$degree=3;  # degree of Nurbs 
if( $splitAnnulusBC eq "" ){ $splitAnnulusBC=8; }
$pi = 4.*atan2(1.,1.);
# $ra =.6; $rb =1.;
$ra =.3; $rb =.5;
$phi = $pi/16.; # half gap angle
$theta1 = $pi;
$theta2 = 2*$pi-$phi;
#
$cmd="#";  $ns=0;  $arcLength=0.;
#
# Segment I : outer annulus, lower half
#
$NI = 15; $theta1=$pi; $theta2=2.*$pi-$phi; 
$dTheta = ($theta2-$theta1)/($NI); 
$theta=$theta1; 
$x0=$rb*cos($theta); $y0=$rb*sin($theta);
for( $i=0; $i<=$NI; $i++ ){ $theta= $theta1 + ($i)*$dTheta; $ns=$ns+1;  \
   $x=$rb*cos($theta); $y=$rb*sin($theta);   \
   $radX = max($radX,abs($x)); $radY = max($radY,abs($y));  \
   $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 ); $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; }
# 
# Segment II -- lower radial segement
#
$NII=3;
$theta=-$phi;
$dr = ($rb-$ra)/($NII);
for( $i=1; $i<=$NII; $i++ ){ $r= $rb - ($i)*$dr; $ns=$ns+1;  \
   $x=$r*cos($theta); $y=$r*sin($theta);   \
   $radX = max($radX,abs($x)); $radY = max($radY,abs($y));  \
   $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 ); $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; }
#
# Segment III : inner annulus
$NIII = 20; $theta1=2*$pi-$phi; $theta2=$phi; 
$dTheta = ($theta2-$theta1)/($NIII); 
for( $i=1; $i<=$NIII; $i++ ){ $theta= $theta1 + ($i)*$dTheta; $ns=$ns+1;  \
   $x=$ra*cos($theta); $y=$ra*sin($theta);   \
   $radX = max($radX,abs($x)); $radY = max($radY,abs($y));  \
   $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 ); $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; }
# 
# Segment IV -- upper radial segement
#
$NIV=3;
$theta=$phi;
$dr = ($rb-$ra)/($NIV);
for( $i=1; $i<=$NIV; $i++ ){ $r= $ra + ($i)*$dr; $ns=$ns+1;  \
   $x=$r*cos($theta); $y=$r*sin($theta);   \
   $radX = max($radX,abs($x)); $radY = max($radY,abs($y));  \
   $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 ); $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; }   
#
# Segment V -- upper annlus
$NV = $NI; $theta1=$phi; $theta2=$pi; 
$dTheta = ($theta2-$theta1)/($NV); 
for( $i=1; $i<=$NV; $i++ ){ $theta= $theta1 + ($i)*$dTheta; $ns=$ns+1;  \
   $x=$rb*cos($theta); $y=$rb*sin($theta);   \
   $radX = max($radX,abs($x)); $radY = max($radY,abs($y));  \
   $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 ); $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; }
#
$knots="#"; for( $i=$degree-1; $i<$ns-($degree-1); $i++ ){ $s=$i/($ns-2); $knots .= "\n $s"; }
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
  # open graphics
   splitAnnulusBoundaryInitial
exit
# -- interpolate the initial NURBS so that we have an arclength parameterization --
nurbs (curve)
  interpolate from a mapping
    splitAnnulusBoundaryInitial
 lines
  $lines=intmg($arcLength/$ds + 1.5 );
  $lines
  mappingName
   splitAnnulusBoundary
exit 
# 
# -- Make a hyperbolic grid --
#
  $nr = intmg( 6 + $order/2 );
  hyperbolic
    forward
    $nDist=($nr-2.5)*$ds;
    distance to march $nDist
    $nrm=$nr-1; 
    lines to march $nrm
    $nTheta = int($arcLength/$ds+1.5);
    points on initial curve $nTheta
    uniform dissipation 0.075
    volume smooths $numberOfVolumeSmooths
    equidistribution 0 (in [0,1])
    #
    spacing: geometric
    geometric stretch factor 1.05 
    #
    generate
    # open graphics
    boundary conditions
      -1 -1 $splitAnnulusBC 0 0 0
    share 
       0 0 100 0 0 0
    name splitAnnulusGridBase
    # open graphics
  exit