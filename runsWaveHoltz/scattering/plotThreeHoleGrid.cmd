#
# plotStuff plotThreeHoleGrid.cmd -grid=/home/henshw/grids/threeHolesGridN2M1e2.order2.hdf -name=threeHolesGrid
#
$grid=""; $name=""; 
GetOptions( "grid=s"=>\$grid,"name=s"=>\$name );
#
$grid
#
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 4
  bigger
  plot
  $plotName = $name . ".ps";
  hardcopy file name:0 $plotName
  hardcopy save:0  
