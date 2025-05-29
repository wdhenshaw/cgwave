#
# plotStuff plotSixteenKnivesGrid.cmd -grid=/home/henshw/grids/sixteenKnivesGride64.order2.hdf -name=sixteenKnivesG64
#
$grid="sixteenKnivesGride16.order2.hdf"; $name="sixteenKnivesGride16"; 
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
  plot
  $plotName = $name . ".ps";
  hardcopy file name:0 $plotName
  hardcopy save:0  
  pause 
#
  set view:0 0.603792 -0.271477 0 8.29226 1 0 0 0 1 0 0 0 1
#   set view:0 0.000529337 -0.323342 0 26.1748 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom1.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0 
  pause
#
  # set view:0 0.675519 -0.322781 0 37.0314 1 0 0 0 1 0 0 0 1
  set view:0 0.675008 -0.326838 0 205.674 1 0 0 0 1 0 0 0 1
  plot
  $plotName = $name . "Zoom2.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0
  pause