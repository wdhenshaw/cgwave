#
# plotStuff plotGrid3d.cmd -grid=pipeClosedze2.order2.hdf -name=pipeGridG2
#
$grid="pipeze2.order2.hdf"; $name="pipeGridG2"; 
GetOptions( "grid=s"=>\$grid,"name=s"=>\$name );
#
$grid
#
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  colour boundaries by chosen name
  pick colour...
  PIC:sky blue
  pick closest 1
  grid colour 1 SKYBLUE
  PIC:aquamarine
  grid colour 0 AQUAMARINE
  smaller 1.05
  set view:0 0.00951662 -0.00634441 0 0.985119 0.939693 0.0593912 -0.336824 0 0.984808 0.173648 0.34202 -0.163176 0.925417
#
  hardcopy file name:0 $name.ps
  hardcopy save:0
