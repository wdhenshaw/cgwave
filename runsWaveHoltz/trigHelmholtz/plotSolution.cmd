#
# 
#   plotStuff plotSolution.cmd -show=shapesHelmholtzFreq11p5 -solution=1
#
$show="gaussianSquare"; $solution="-1"; $name=""; $field="Ey"; $emin=0; $emax=-1; $res=1024; $cl=0; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"res=i"=>\$res,"cl=i"=>\$cl );
#
if( $name eq "" ){ $name=$show; }
$show
#
contour
  plot:v
  if( $cl eq "0" ){ $cmd="plot contour lines (toggle)"; } else{ $cmd="#"; }
  $cmd
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  if( $emax > $emin ){ $cmd="min max $emin $emax"; }else{ $cmd="#"; }
  $cmd
exit
x-
solution: $solution
pause
#
DISPLAY AXES:0 0
if( $res ne 1024 ){ $cmd="hardcopy vertical resolution:0 $res\n hardcopy horizontal resolution:0 $res\n line width scale factor:0 3"; }else{ $cmd="#"; }
$cmd 
# DISPLAY LABELS:0 0
# DISPLAY COLOUR BAR:0 0
# 
plot
$plotName = $name . "WaveHoltz.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
plot:err0
$plotName = $name . "WaveHoltzErr.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0