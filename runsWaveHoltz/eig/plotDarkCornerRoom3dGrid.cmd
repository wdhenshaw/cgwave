# 
#   plotStuff plotDarkCornerRoom3dGrid.cmd -show=darkCornerRoom3dGride1.order2.hdf -name=darkCornerRoom3dGridG1 
#
$show="darkCornerRoom3dGride1.order2.hdf";  $name="";
$cf=2; $res=1024; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name,"cf=i"=>\$cf,"res=i"=>\$res );
#
$show
#
  plot block boundaries 0
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  coarsening factor 4
  CLIP:0 0 on
  coarsening factor 2
  CLIP:0 0 distance 2.500000e-01
  set view:0 0 0 0 1.27308 1 0 0 0 0.984808 0.173648 0 -0.173648 0.984808
  colour boundaries by grid number
  grid colour 3 BRASS
# 
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 3
  plot
#
  hardcopy file name:0 darkCornerRoom3dGridG1.ps
  hardcopy save:0
  pause
#
  set view:0 -0.0154277 0.57266 0 2.87581 1 0 0 0 0.965926 0.258819 0 -0.258819 0.965926
#  -- zoom at bottom
  hardcopy file name:0 darkCornerRoom3dGridG1Zoom.ps
  hardcopy save:0






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


# --- for TOP HALF --
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