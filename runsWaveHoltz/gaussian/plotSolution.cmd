#
# 
#   plotStuff plotSolution.cmd -show=gaussianSquare.show -name=gaussianSquareOmega -solution=1
#   plotStuff plotSolution.cmd -show=gaussianShapesOmega9p777.show -name=gaussianShapesOmega9p777 -solution=1
#   plotStuff plotSolution.cmd -show=gaussianDiskOmega8p1.show -name=gaussianDiskOmega8p1 -solution=1
#   plotStuff plotSolution.cmd -show=annulusOmega9p777.show -name=annulusOmega9p777 -solution=1
#
#   plotStuff plotSolution.cmd -show=shapesOmega66.show -name=shapesOmega66 -solution=1
#   plotStuff plotSolution.cmd -show=shapesOmega67.show -name=shapesOmega67 -solution=1
#   plotStuff plotSolution.cmd -show=shapesOmega68.show -name=shapesOmega68 -solution=1
#   plotStuff plotSolution.cmd -show=shapesOmega69.show -name=shapesOmega69 -solution=1
#   plotStuff plotSolution.cmd -show=shapesOmega70.show -name=shapesOmega70 -solution=1
#
# plotStuff plotSolution.cmd -show=shapesG4O4freq45Np80 -name=shapesG4O4freq45Np80 -solution=1
# plotStuff plotSolution.cmd -show=shapesG4O4freq44Np80 -name=shapesG4O4freq44Np80 -solution=1
# plotStuff plotSolution.cmd -show=shapesG4O4freq43Np80 -name=shapesG4O4freq43Np80 -solution=1
#
$show="gaussianSquare.show"; $solution="-1"; $name="plot"; $field="Ey"; $emin=0; $emax=-1; $numFreq=1; $clines=0; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"numFreq=i"=>\$numFreq,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax );
#
$show
#
contour
  plot:v0
  if( $clines ==0 ){ $cmd="plot contour lines (toggle)"; }else{ $cmd="#"; }
  $cmd 
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  if( $emax > $emin ){ $cmd="min max $emin $emax"; }else{ $cmd="#"; }
  $cmd
exit
solution: $solution
pause
#
DISPLAY AXES:0 0
# DISPLAY LABELS:0 0
# DISPLAY COLOUR BAR:0 0
# 
$cmd="#"; 
for( $i=0; $i<$numFreq; $i++ ){ $cmd .= "\n plot:v$i\n \$plotName = $name . \"v$i.ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd



plot
$plotName = $name . "WaveHoltz.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0


#
# -- forcing:
contour
  plot:f0
  plot contour lines (toggle)
  plot boundaries (toggle)
  vertical scale factor 0.7
exit
x-r:0 65
$plotName = $name . "WaveHoltzForcing.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0




# -- OLD STUFF
Foreground colour:0 white
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
# 
plot
$plotName = $name . "EfieldNorm.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0