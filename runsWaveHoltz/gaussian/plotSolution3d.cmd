# 
#   plotStuff plotSolution.cmd -show=sphereOneG4Freq18.show -name=sphereG4O2Omega18 -solution=1
#
$show="gaussianSquare.show"; $solution="-1"; $name="plot"; $field="Ey"; $emin=0; $emax=-1; $numFreq=1; $clines=0; $cf=2; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"numFreq=i"=>\$numFreq,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"cf=i"=>\$cf );
#
$show
#
contour
  contour lines 0
  plot the grid
    plot shaded surfaces (3D) 0
    coarsening factor $cf
    plot block boundaries 0
  exit this menu
  # delete contour plane 1
  # delete contour plane 0
  DISPLAY AXES:0 0


  $plotName = "$name.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
