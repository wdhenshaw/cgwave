#
# 
#   plotStuff plotSolution3d.cmd -show=pipeG4O2ImpEig.show -name=pipeG2O2Eigs -start=99
#   plotStuff plotSolution3d.cmd -show=pipeG4O2ImpEig.show -name=pipeG2O2Eigs -start=98 -cs=.1
#   plotStuff plotSolution3d.cmd -show=pipeG4O2ImpEig.show -name=pipeG2O2Eigs -start=97
#   plotStuff plotSolution3d.cmd -show=pipeG4O2ImpEig.show -name=pipeG2O2Eigs -start=96 -cs=.1
#   plotStuff plotSolution3d.cmd -show=pipeG4O2ImpEig.show -name=pipeG2O2Eigs -start=9 -cs=.1
#   plotStuff plotSolution3d.cmd -show=pipeG4O2ImpEig.show -name=pipeG2O2Eigs -start=45
#
# plotStuff plotSolution3d.cmd -show=sphereG2O2ImpEig.show -name=sphereG2O2Eigs -start=0 
# plotStuff plotSolution3d.cmd -show=sphereG2O2ImpEig.show -name=sphereG2O2Eigs -start=9
#
$show="pipeG2O2ImpEig.show"; $solution="-1"; $name="plot"; $field="phi"; $cf=2; 
$emin=0; $emax=-1; $numFreq=1; $clines=0; 
$stride=1; $start=0; $cs=0; 
$numToSave=1; # save this many components
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"cs=f"=>\$cs,\
      "stride=i"=>\$stride, "numToSave=i"=>\$numToSave,"start=i"=>\$start,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"cf=i"=>\$cf );
#
$show
#
contour
  pick to delete contour planes
  # pipe: 
  # delete contour plane 2
  # add contour plane  0.00000e+00  0.00000e+00  1.00000e+00 -1.45649e-02  2.20647e-01  4.45946e-01 
  # 
  if( $emax > $emin ){ $cmd="min max $emin $emax"; }else{ $cmd="#"; }
  $cmd 
  # min max -1 1
  plot the grid
    plot shaded surfaces (3D) 0
    plot block boundaries 0
    coarsening factor $cf
  exit this menu
  DISPLAY AXES:0 0
  # pipe:
  set view:0 -0.117447 -0.0113293 0 0.814769 0.866025 0.17101 -0.469846 0 0.939693 0.34202 0.5 -0.296198 0.813798
  # sphere 
  set view:0 -0.0996487 -0.0231948 0 1.05019 0.866025 0.17101 -0.469846 0 0.939693 0.34202 0.5 -0.296198 0.813798
#  
  contour shift $cs 
  +shift contour planes
#
exit
plot:phi$start
pause 
#
  line width scale factor:0 3
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
$cmd="#"; 
for( $i=$start; $i<$start+$numToSave; $i=$i+$stride ){ $cmd .= "\n plot:$field$i\n \$plotName = $name . \"$field$i.ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd


next
next
bigger:0
smaller:0
smaller:0
x-:0
x-:0
DISPLAY AXES:0 0