#
# Two adjacent squares
#
#
# usage: ogen [-noplot] twoSquaresGrid -nx1=<i> -nx2=<i> -xa1-<f> -xb1=<f> -xa2=<f> -xb2=<f> -ya=<f> -yb=<f> -ny=<i> ...
#              -order=[2/4/6/8] -interp=[e/i] -per=[0|1]
# 
# examples:
#     ogen -noplot twoSquaresGrid -nx1=54 -nx2=16 -xb1=.325 -xa2=.25 -order=2 -per=1 -interp=e 
#     ogen -noplot twoSquaresGrid -nx1=54 -nx2=16 -xb1=.325 -xa2=.25 -order=4 -per=1 -interp=e 
#
$prefix="twoSquaresGrid"; 
$xa1=-1.; $xb1=.1; $nx1=11; # left grid 
$xa2=-.1; $xb2=1.; $nx2=11; # right grid 
$ya=0.; $yb=1.; $ny=11; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$per=0; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=i"=> \$factor,"interp=s"=> \$interp,"per=i"=>\$per,"prefix=s"=>\$prefix,\
             "xa1=f"=>\$xa1,"xb1=f"=>\$xb1,"xa2=f"=>\$xa2,"xb2=f"=>\$xb2,"ya=f"=>\$ya,"yb=f"=>\$yb,\
             "nx1=i"=>\$nx1,"nx2=i"=>\$nx2,"ny=i"=>\$ny);
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=5; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $per eq 1 ){ $suffix .= "p"; }
$name = $prefix . "$interp" . $suffix . ".hdf";
# 
# $ds=.1/$factor;
# $width = ($order-2)/2;
# if( $interp eq "e" ){ $width=$width+1.; }
# $overlap = $ds*$width + $ds*.125;
# 
create mappings
  rectangle
    set corners
     $xa1 $xb1 $ya $yb 
    lines
      # $ny=int( ($yb-$ya)/$ds+1.5 );
      $nx1 $ny
    boundary conditions
      if( $per eq 0 ){ $cmd="1 0 2 2"; }else{ $cmd="1 0 -1 -1"; }
      $cmd
    share
      0 0 1 2
    mappingName
      leftSquare
    exit
# 
  rectangle
    set corners
     $xa2 $xb2 $ya $yb 
    lines
      $nx2 $ny
    boundary conditions
      if( $per eq 0 ){ $cmd="0 1 2 2"; }else{ $cmd="0 1 -1 -1"; }
      $cmd
      # 0 1 2 2 
    share
      0 0 1 2
    mappingName
      rightSquare
    exit
#
  exit
#
generate an overlapping grid
  leftSquare
  rightSquare
  done
  change parameters
    ghost points
      all
       $ng $ng $ng $ng $ng $ng 
    order of accuracy
      $orderOfAccuracy
    interpolation type
      $interpType
  exit
  compute overlap
exit
save a grid (compressed)
$name
sis
exit
