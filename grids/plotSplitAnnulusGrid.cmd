#
# plotStuff plotSplitAnnulusGrid.cmd -grid=/home/henshw/runs/cgwave/grids/splitAnnulusGridLeftGape4.order2.hdf -name=splitAnnulusGrid
# plotStuff plotSplitAnnulusGrid.cmd -grid=/home/henshw/runs/cgwave/grids/splitAnnulusGridLeftGape8.order2.hdf -name=splitAnnulusGrid
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
  set view:0 0.0120955 -0.0193463 0 1.25741 1 0 0 0 1 0 0 0 1
  plot
  $plotName = $name . ".ps";
  hardcopy file name:0 $plotName
  hardcopy save:0 
  pause 
  #
  set view:0 0.218608 -0.0132167 0 6.92377 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0 