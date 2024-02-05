#
# Random holes in a square
#
#
# usage: ogen [-noplot] holesGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] -nCylX=<i> -nCyly=<i> -deltaX=<f> -deltaY=<f> ...
#                             -xa=<> -xb=<> -ya=<> -yb=<> -blf=<num> -ml=<> -rgd=[fixed|var] -offsetColumns=[0|1] -numGhost=<i> -periodic=[|p|np|pn]
# 
#  -blf : boundary-layer-factor : blf>1 : make grid lines near boundary this many times smaller
#  -ml = number of (extra) multigrid levels to support
#  -rgd : var=variable : decrease radial grid distance as grids are refined. fixed=fix radial grid distance
#  -xa, -xb, -ya, -yb : bounds on the back ground grid
#  -cx, -cy : center for the annulus
# 
# examples:
#
#     ogen -noplot holesGrid -order=2 -interp=e -blf=2 -innerRad=.05 -periodic=np -xa=-3 -nCylx=3 -nCyly=5 -deltaX=0.5 -deltaY=0.5 -offsetColumns=1 -factor=4 
# 
#     ogen -noplot holesGrid -order=2 -interp=e -blf=2 -innerRad=.05 -periodic=np -xa=-3 -nCylx=7 -nCyly=15 -deltaX=.25 -deltaY=.25 -factor=8 
#     ogen -noplot holesGrid -order=4 -interp=e -blf=2 -innerRad=.05 -periodic=np -xa=-3 -nCylx=7 -nCyly=15 -deltaX=.25 -deltaY=.25 -numGhost=3 -factor=16
# 
#     ogen -noplot holesGrid -order=4 -interp=e -blf=2 -innerRad=.02 -periodic=np -xa=-3 -nCylx=7 -nCyly=17 -deltaX=.2 -deltaY=.2 -numGhost=3 -factor=16
#     ogen -noplot holesGrid -order=4 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01 -innerRad=.01 -periodic=np -xa=-3 -nCylx=7 -nCyly=18 -deltaX=.2 -deltaY=.2 -numGhost=3 -factor=16
#     ogen -noplot holesGrid -order=4 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01 -innerRad=.01 -periodic=np -xa=-3 -nCylx=7 -nCyly=26 -deltaX=.15 -deltaY=.15 -numGhost=3 -factor=16
#     ogen -noplot holesGrid -order=4 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01 -innerRad=.01 -periodic=np -xa=-3 -nCylx=11 -nCyly=25 -deltaX=.15 -deltaY=.15 -numGhost=3 -factor=16
#
#     1x26 
#     ogen -noplot holesGrid -order=4 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01 -innerRad=.01 -periodic=np -xa=-3 -nCylx=1 -nCyly=26 -deltaX=.15 -deltaY=.15 -numGhost=3 -factor=16
#     3x26 
#     ogen -noplot holesGrid -order=4 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01 -innerRad=.01 -periodic=np -xa=-3 -nCylx=3 -nCyly=26 -deltaX=.15 -deltaY=.15 -numGhost=3 -factor=16
#     7x26 
#     ogen -noplot holesGrid -order=4 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01 -innerRad=.01 -periodic=np -xa=-3 -nCylx=7 -nCyly=26 -deltaX=.15 -deltaY=.15 -numGhost=3 -factor=16
#     ogen -noplot holesGrid -order=4 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01Long -innerRad=.01 -periodic=np -xa=-3.5 -xb=3 -nCylx=7 -nCyly=26 -deltaX=.15 -deltaY=.15 -numGhost=3 -factor=16
#     ogen -noplot holesGrid -order=2 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01Long -innerRad=.01 -periodic=np -xa=-3.5 -xb=3 -nCylx=7 -nCyly=26 -deltaX=.15 -deltaY=.15 -factor=16
#     7x26 offset
#     ogen -noplot holesGrid -order=4 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01Offset -offsetColumns=1 -innerRad=.01 -periodic=np -xa=-3 -nCylx=7 -nCyly=26 -deltaX=.15 -deltaY=.15 -numGhost=3 -factor=16
#     ogen -noplot holesGrid -order=4 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01LongOffset -offsetColumns=1 -innerRad=.01 -periodic=np -xa=-3.5 -xb=3 -nCylx=7 -nCyly=26 -deltaX=.15 -deltaY=.15 -numGhost=3 -factor=16
#     ogen -noplot holesGrid -order=2 -interp=e -blf=4 -stretchRes=2 -prefix=holesGridRad01LongOffset -offsetColumns=1 -innerRad=.01 -periodic=np -xa=-3.5 -xb=3 -nCylx=7 -nCyly=26 -deltaX=.15 -deltaY=.15 -factor=16
# 
# 
# TROUBLE SOLVING WITH HIGH STRETCHING RATIOS
#  ogen -noplot holesGrid -order=4 -interp=e -blf=6 -stretchRes=1 -innerRad=.01 -periodic=np -xa=-3 -nCylx=7 -nCyly=19 -deltaX=.2 -deltaY=.2 -numGhost=3 -factor=8
# 
#  ogen -noplot holesGrid -order=4 -interp=e -blf=6 -stretchRes=1 -innerRad=.01 -periodic=np -xa=-3 -nCylx=11 -nCyly=27 -deltaX=.15 -deltaY=.15 -numGhost=3 -factor=16
#
#
$prefix="holesGrid";  $rgd="var";
$bcSquare="d"; # old way
$periodic="";  # new way 
$dw=""; $iw=""; 
$order=2; $factor=1; $interp="i"; $ml=0; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-2.; $xb=2.; $ya=-2.; $yb=2.; 
$cx=0.; $cy=0.;  # center for the annulus
$blf=1;  # this means no stretching
$deltaRadius0=.3; # radius for rgd fixed
$numGhost=-1;  # if this value is set, then use this number of ghost points
$innerRad=.5; 
$stretchRes=1.5; # stretch resolution factor to increase points in stretching direction
$offsetColumns=0; # 1=offset every second column by .5*$deltaY
#
$nCylx=2; $nCyly=3; $deltaX=1; $deltaY=1;
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"ml=i"=>\$ml,"blf=f"=> \$blf, "prefix=s"=> \$prefix,\
            "cx=f"=>\$cx,"cy=f"=>\$cy,"rgd=s"=> \$rgd,"bcSquare=s"=>\$bcSquare,"numGhost=i"=>\$numGhost,\
            "iw=i"=>\$iw,"dw=i"=>\$dw,"periodic=s"=>\$periodic,"innerRad=f"=>\$innerRad,"offsetColumns=i"=>\$offsetColumns,\
            "nCylx=i"=>\$nCylx,"nCyly=i"=>\$nCyly,"deltaX=f"=>\$deltaX,"deltaY=f"=>\$deltaY,"stretchRes=f"=>\$stretchRes );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=3; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=4; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
if( $rgd eq "fixed" ){ $prefix = $prefix . "Fixed"; }
# if( $bcSquare eq "p" ){ $prefix = $prefix . "p"; }
$suffix=""; 
if( $periodic eq "p" ){ $suffix = "p"; }
if( $periodic eq "np" ){ $suffix = "np"; }
if( $periodic eq "pn" ){ $suffix = "pn"; }
if( $iw eq "" ){ $suffix .= ".order$order"; }else{ $suffix .= ".Iw$iw" . "Dw$dw" . "$bc"; }
# $suffix = ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
# if( $blf ne 1 ){ $suffix .= ".s$blf"; }
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $name eq "" ){$name = $prefix . "N$nCylx" . "M$nCyly" . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.1/$factor;
$pi = 4.*atan2(1.,1.);
# 
if( $dw eq "" ){ $dw = $order+1; $iw=$order+1; }
# parallel ghost lines: for ogen we need at least:
#       .5*( iw -1 )   : implicit interpolation 
#       .5*( iw+dw-2 ) : explicit interpolation
$parallelGhost=($iw-1)/2;
if( $interp eq "e" ){  $parallelGhost=($iw+$dw-2)/2; }
if( $parallelGhost<1 ){ $parallelGhost=1; } 
minimum number of distributed ghost lines
  $parallelGhost
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
create mappings
#
rectangle
  set corners
    $xa $xb $ya $yb
  lines
    $nx = intmg( ($xb-$xa)/$ds +1.5 ); 
    $ny = intmg( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
  boundary conditions
    # $sbc="1 2 3 4";
    # if( $bcSquare eq "p" ){ $sbc = "-1 -1 3 4"; }
    # $sbc
    if( $periodic eq "p" ){ $bc ="-1 -1 -1 -1"; }\
    elsif( $periodic eq "np" ){ $bc ="1 2 -1 -1"; }\
    elsif( $periodic eq "pn" ){ $bc ="-1 -1 3 4"; }else{ $bc="1 2 3 4"; }   
    $bc 
  mappingName
  square
exit
#
# if 
if( $blf>1 ){ $annulusName="AnnulusUnStretched"; $stretchAnnulusName="Annulus"; }else{ $annulusName="Annulus"; $stretchAnnulusName="AnnulusStretched"; }
Annulus
  # Make sure there are at least 4 points on the coarsest MG level
  # $nr = max( 5+ $ng + 2*($order-2), 2**($ml+2) );
  $nr = max( 5+ $ng, 2**($ml+2) );
  $nr = intmg( $nr );
  # $innerRad=.5; 
  $outerRad = $innerRad + ($nr-1)*$ds;
  if( $rgd eq "fixed" ){ $outerRad = $innerRad + $deltaRadius0; $nr=intmg( $deltaRadius0/$ds + 2.5 ); }
  center: $cx $cy
  inner and outer radii
    $innerRad $outerRad
  lines
    if( $blf>1 ){ $nr = $nr + 4; } # extra grid lines to account for stretching
    # $nTheta = intmg( 2.*$pi*($innerRad+$outerRad)*.5/$ds + 1.5 );
    $nTheta = intmg( 2.*$pi*( .1*$innerRad+ .9*$outerRad )/$ds + .5 );
    $nTheta $nr
  boundary conditions
    -1 -1 5 0
  share
     0  0 5 0
  mappingName
   $annulusName
exit
#
# optionally stretch the grid lines next to the cylinder
# 
 stretch coordinates 
  transform which mapping? 
    $annulusName 
  multigrid levels $ml
  # add extra resolution in the stretching direction: 
  #
  stretch resolution factor $stretchRes
  # exponential to linear stretching: 
   Stretch r2:exp to linear
   STP:stretch r2 expl: position 0
   $dxMin = $ds/$blf; 
   STP:stretch r2 expl: min dx, max dx $dxMin $ds
  #Stretch r2:itanh
  #STP:stretch r2 itanh: position and min dx 0 $dxMin
  #stretch grid
  STRT:name $stretchAnnulusName
 exit
#
# --- Now make a collection of holes ------------
#
$count=0; $mappingsList="#"; 
sub makeDisk\
{ local($xc,$yc)=@_; \
  $count = $count + 1; $layer=1; \
  $mappingName="annulus$count";     \
  $mappingsList .= "\n$mappingName";     \
  $commands = \
  "rotate/scale/shift \n" . \
  "  transform which mapping? \n" . \
  "  Annulus \n" . \
  "  shift \n" . \
  "    $xc $yc\n" . \
  "  mappingName \n" . \
  "  $mappingName \n" . \
  "  exit \n"; \
}
# ===================================================
#   Make an array of cylinders
#     makeDiskArray(radius,nCylx,nCyly,x0,dx0,y0,dy0)
#     
# Make cylinders at centers
#      (x0+i*dx0,y0+j*dy0)  i=0,..,nx0, j=0,..,ny0
# 
# Result: $commands
# =====================================================
sub makeDiskArray \
{ local($nCylx,$nCyly,$x0,$dx0,$y0,$dy0)=@_; \
  local $cmds; $cmds=""; \
  if( $offsetColumns ){ $mod2=2; }else{ $mod2=1; } \
  for( $i=0; $i<$nCylx; $i++ ){ \
  for( $j=-($i%$mod2); $j<$nCyly; $j++ ){ \
    makeDisk($x0+$i*$dx0,$y0+$j*$dy0 + ($i%$mod2)*$dy0*.5 ); \
    $cmds = $cmds . $commands; \
  }}\
  $commands=$cmds; \
}
#
# 
$x0=-($nCylx-1)*$deltaX*.5; $y0=-($nCyly-1)*$deltaY*.5;  
makeDiskArray($nCylx,$nCyly,$x0,$deltaX,$y0,$deltaY);
#
$commands
# printf("command=$commands\n");
# printf("mappingsList=$mappingsList\n");
#
exit
generate an overlapping grid
    square
    $mappingsList
  done
  change parameters
    # choose implicit or explicit interpolation
    interpolation type
      $interpType
    #
    $cmd =" order of accuracy\n $orderOfAccuracy";
    $cmd
    #
    ghost points
      all
      # $ngp = $ng+1;
      $ngp = $ng;
      $ng $ng $ng $ngp $ng $ng
  exit
#   display intermediate results
# open graphics 
#
  compute overlap
#*  display computed geometry
  exit
#
# save an overlapping grid
save a grid (compressed)
# printf(" name=$name\n");
$name
cic
exit

