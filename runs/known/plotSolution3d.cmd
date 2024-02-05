#
# 
#   plotStuff plotSolution3d.cmd -show=pipeG4O4.show -name=pipeG4O4Eig -solution=4
#   plotStuff plotSolution3d.cmd -show=sphereG4O4.show -name=sphereG4O4Eig -solution=6 -cs=-.2 -cf=3
#
$show="pipeG4O4.show"; $solution="-1"; $name="plot"; $field="phi"; $cf=2; 
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
solution: $solution
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
  $cmd="#";
  if( $show =~ "pipe*"   ){ $cmd = "set view:0 -0.117447 -0.0113293 0 0.814769 0.866025 0.17101 -0.469846 0 0.939693 0.34202 0.5 -0.296198 0.813798"; }
  if( $show =~ "sphere*" ){ $cmd = "set view:0 -0.0966823 -0.0231948 0 1.0601 0.866025 0.17101 -0.469846 0 0.939693 0.34202 0.5 -0.296198 0.813798"; }
  $cmd
  printf("cmd=$cmd\n");
  #  
  contour shift $cs 
  +shift contour planes
#
exit
pause
#
  line width scale factor:0 3
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
# 
plot:u 
$plotName = $name . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
plot:err
$plotName = $name . "Err.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0




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