# 
#   plotStuff plotDarkCornerRoom3dGrid.cmd -show=darkCornerRoom3dGride1.order2.hdf -name=darkCornerRoom3dGridG1 -cf=2
#
$show="darkCornerRoom3dGride1.order2.hdf"; 
$cf=2; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name,"cf=i"=>\$cf );
#
$show
#
  toggle shaded surfaces 3 0
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  plot block boundaries 0
  plot block boundaries 1
  coarsening factor $cf
  pick closest 1
  colour boundaries by grid number
  set view:0 0 -0.00483384 0 1.27308 0.866025 0.0868241 -0.492404 0 0.984808 0.173648 0.5 -0.150384 0.852869
# 
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  hardcopy file name:0 darkCornerRoom3dGridG1.ps
  hardcopy save:0
pause  
#
set view:0 0.00486011 -0.232272 0 2.56394 0.866025 -0.17101 -0.469846 0 0.939693 -0.34202 0.5 0.296198 0.813798
#
hardcopy file name:0 darkCornerRoom3dGridG1Zoom.ps
hardcopy save:0