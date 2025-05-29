# 
#   plotStuff plotSolution3d.cmd -show=sphereOneG4Freq18.show -name=sphereG4O2Omega18 -solution=1
#
$show="gaussianSquare.show"; $solution="-1"; $name="plot"; $field="Ey"; $emin=0; $emax=-1; $numFreq=1; $clines=0; $cf=2; $comp="v"; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
$res=1024; # hardcopy resolution
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"numFreq=i"=>\$numFreq,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"cf=i"=>\$cf,"res=i"=>\$res, "comp=s"=>\$comp );
#
$show
#
if( $comp eq "absv" ){ $cmd="derived types\n absoluteValue\n v0  (off)\n done\n exit"; }else{ $cmd="#" }
$cmd
contour
  $compi = $comp . "0"; 
  plot:$compi
  contour lines 0
  plot the grid
    plot shaded surfaces (3D) 0
    coarsening factor $cf
    plot block boundaries 0
  exit this menu
  # delete contour plane 1
  # delete contour plane 0
  DISPLAY AXES:0 0
  if( $res ne 1024 ){ $cmd="hardcopy vertical resolution:0 $res\n hardcopy horizontal resolution:0 $res\n line width scale factor:0 3\nplot"; }else{ $cmd="#"; }
$cmd 
# darkCornerRoom3d:
 set view:0 -0.0879154 -0.00906344 0 1.05079 0.939693 0.116978 -0.321394 0 0.939693 0.34202 0.34202 -0.321394 0.883022
  hardcopy file name:0 darkCornerRoom3dG2O2Freq3p89.ps
  hardcopy save:0

  $plotName = "$name.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
