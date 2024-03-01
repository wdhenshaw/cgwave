#
# plotStuff plotSmallHoleGrid.cmd -grid=holeGridRad0p01Fixede2.order2 -name=holeGridRad0p01G2
#
$grid="sice1.order2.hdf"; $name="sicGrid"; 
GetOptions( "grid=s"=>\$grid,"name=s"=>\$name );
#
$grid
#
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
   coarsening factor 1
  bigger:0
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 3
#  
  plot interpolation points 1
  $plotName = $name . ".ps";
  hardcopy file name:0 $plotName
  hardcopy save:0
  #
  set view:0 0.00136505 0.00401252 0 6.38852 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0
#
  set view:0 -1.38707e-06 -0.000304041 0 71.6055 1 0 0 0 1 0 0 0 1
  hardcopy file name:0 holeGridRad0p01G2Zoom2.ps
  hardcopy save:0



  plot interpolation points 1
  colour interpolation points 1
  hardcopy file name:0 $name.ps
  hardcopy save:0
