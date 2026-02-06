#
# plotStuff plotMultiObjectsGrid.cmd -grid=/home/henshw/grids/multiObjectsGride3.order2 -name=multiObjectsGrid
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
  colour boundaries black
#
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 2
  set view:0 0.00966767 -0.00483384 0 1.28894 1 0 0 0 1 0 0 0 1
  plot
  $plotName = $name . ".ps";
  hardcopy file name:0 $plotName
  hardcopy save:0 
  pause 
  #
  set view:0 0.239748 -0.182816 0 6.25581 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0 