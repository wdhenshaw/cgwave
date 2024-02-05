#
# plotStuff plotGrid.cmd -grid=holesGridRad01LongN7M26e16np.order2 -name=holesGridG16
# plotStuff plotGrid.cmd -grid=holesGridRad01LongOffsetN7M26e16np.order2 -name=offsetHolesGridG16
#
$grid="sice1.order2.hdf"; $name="sicGrid"; 
GetOptions( "grid=s"=>\$grid,"name=s"=>\$name );
#
$grid
  grid colour 15 RED
  pick to colour grids
  colour grid lines from chosen name
  grid colour 16 RED
  grid colour 42 RED
  grid colour 41 RED
  grid colour 14 RED
  grid colour 40 RED
  grid colour 68 RED
  grid colour 67 RED
  grid colour 66 RED
  grid colour 69 RED
  grid colour 43 RED  
  set view:0 0.0190566 -0.052984 0 20.7345 1 0 0 0 1 0 0 0 1
#
pause
#
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 3
#
  plot interpolation points 1
  # colour interpolation points 1
  hardcopy file name:0 $name.ps
  hardcopy save:0
