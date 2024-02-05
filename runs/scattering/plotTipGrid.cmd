#
# plotStuff plotTipGrid.cmd -grid=tipGridThine64.order2.hdf -name=tipGridThine64
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
#
  set view:0 -0.0955491 0.156643 0 5.66411 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom1.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0 
#
  set view:0 -0.0844205 0.0373909 0 20.4668 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom2.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0
  #
  line width scale factor:0 8
  plot
  set view:0 -0.0833007 0.00459446 0 454.429 1 0 0 0 1 0 0 0 1   
  $plotName = $name . "Zoom3.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0   