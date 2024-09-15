#
# plotStuff plotTipGrid.cmd -grid=/home/henshw/grids/tipGridSmalle64.order2.hdf -name=tipGridSmall64
#
$grid="sice1.order2.hdf"; $name="sicGrid"; 
GetOptions( "grid=s"=>\$grid,"name=s"=>\$name );
#
$grid
#
#
  # grid colour 3 BLACK
  # colour grid lines from chosen name
#
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 4
  bigger
  $plotName = $name . ".ps";
  hardcopy file name:0 $plotName
  hardcopy save:0  
  pause 
#
  set view:0 -0.00230064 -0.263913 0 6.40532 1 0 0 0 1 0 0 0 1
#   set view:0 0.000529337 -0.323342 0 26.1748 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom1.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0 
  pause
#
  set view:0 0.00016565 -0.327571 0 277.647 1 0 0 0 1 0 0 0 1
  plot
  $plotName = $name . "Zoom2.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0
  pause
