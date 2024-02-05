# 
#   plotStuff plotSolution.cmd -show=shapes.show -name=shapesMGPWt1p0 -vmin=-1 -vmax=1 -solution=11
# 
#   plotStuff plotSolution.cmd -show=diskG4O4.show -name=diskEigsG4O4Imp -solution=11
#
#   plotStuff plotSolution.cmd -show=cylScatG8O4SPIE.show -name=cylScatG8O4SPIE -solution=11
#
$show="gaussianSquare.show"; $solution="-1"; $name="plot"; $field="Ey"; $vmin=0; $vmax=-1; $clines=1; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave, "clines=i"=>\$clines,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,"vmin=f"=>\$vmin,"vmax=f"=>\$vmax );
#
$show
#
contour
  plot:v
  if( $clines eq "0" ){ $cmd="plot contour lines (toggle)"; }else{ $cmd="#"; }
  $cmd
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  if( $vmax > $vmin ){ $cmd="min max $vmin $vmax"; }else{ $cmd="#"; }
  $cmd
exit
solution: $solution
pause
# for disk: 
set view:0 -0.1 -0.0533736 0 1.06284 1 0 0 0 1 0 0 0 1
#
DISPLAY AXES:0 0
# DISPLAY LABELS:0 0
#
$plotName = $name . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
plot:err
$plotName = $name . "Err.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0



DISPLAY COLOUR BAR:0 0
bigger
# 
plot
$plotName = $name . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0