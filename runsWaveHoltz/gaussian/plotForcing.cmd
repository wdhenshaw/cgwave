#
# plotStuff plotForcing.cmd -show=shapesG4O4ThreeFreqNp10 -name=shapesThreeFreqForcing -clines=1 -solution=1 -numFreq=3
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
  plot:f0
  plot boundaries (toggle)
  if( $clines ==0 ){ $cmd="plot contour lines (toggle)"; }else{ $cmd="#"; }
  $cmd 
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.5
  if( $emax > $emin ){ $cmd="min max $emin $emax"; }else{ $cmd="#"; }
  $cmd
exit
solution: $solution
#
DISPLAY AXES:0 0
x-r:0 60
pause
# DISPLAY LABELS:0 0
# DISPLAY COLOUR BAR:0 0
# 
$cmd="#"; 
for( $i=0; $i<$numFreq; $i++ ){ $cmd .= "\n plot:f$i\n \$plotName = $name . \"f$i.ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd


# For disk--
plot forcing
DISPLAY AXES:0 0
plot boundaries (toggle)
set view:0 -0.0848943 0.253021 0 1.04747 1 0 0 0 0.5 -0.866025 0 0.866025 0.5
hardcopy file name:0 diskHelmholtzGaussianForcing0.ps
hardcopy save:0
plot:f1
hardcopy file name:0 diskHelmholtzGaussianForcing1.ps
hardcopy save:0
plot:f2
hardcopy file name:0 diskHelmholtzGaussianForcing2.ps
hardcopy save:0
exit
#
reset:0
x-:0
x-:0
contour
plot:v0
hardcopy file name:0 diskHelmholtzGaussianSolution0.ps
hardcopy save:0
plot:v1
hardcopy file name:0 diskHelmholtzGaussianSolution1.ps
hardcopy save:0
plot:v2
hardcopy file name:0 diskHelmholtzGaussianSolution2.ps
hardcopy save:0