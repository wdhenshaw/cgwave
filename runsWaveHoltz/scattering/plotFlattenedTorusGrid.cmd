#
# plotStuff plotFlattenedTorusGrid.cmd -grid=/home/henshw/grids/flatTorusGridALongXe4.order2.hdf -name=flatTorusG4
#
$grid=""; $name=""; 
GetOptions( "grid=s"=>\$grid,"name=s"=>\$name );
#
$grid
#
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
#
  toggle grid 0 0
  grid colour 1 BRASS
  # grid colour 2 BRASS
  # plot block boundaries 0
  # plot grid lines 0
  set view:0 0.00725076 -0.00725076 0 1.3011 0.866025 0.17101 -0.469846 0 0.939693 0.34202 0.5 -0.296198 0.813798
  smaller
pause
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 4
  plot
  $plotName = $name . ".ps";
  hardcopy file name:0 $plotName
  hardcopy save:0  