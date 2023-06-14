#
# 
#   plotStuff plotDarkCornerRoom3d.cmd -show=darkCornerRoom3dG1O2EigsFreq8.show -name=darkCornerRoom3dG1O2EigsFreq8 -sol=2
#   plotStuff plotDarkCornerRoom3d.cmd -show=darkCornerRoom3dG1O2Eigs.show -name=darkCornerRoom3dG1O2Eigs -sol=144 149
#   plotStuff plotDarkCornerRoom3d.cmd -show=darkCornerRoom3dG1O2Eigs.show -name=darkCornerRoom3dG1O2Eigs -sol=148 150
#   plotStuff plotDarkCornerRoom3d.cmd -show=darkCornerRoom3dG1O2Eigs.show -name=darkCornerRoom3dG1O2Eigs -sol=74 78
#   plotStuff plotDarkCornerRoom3d.cmd -show=darkCornerRoom3dG1O2EigsEv32.show 
#
$show="pipeG2O2ImpEig.show"; $name="plot"; $field="phi"; $cf=2; 
$emin=0; $emax=-1; $numFreq=1; $clines=0; 
$cs=0; 
$numToSave=1; # save this many components
@sol= ();  # this must be null for GetOptions to work, defaults are given below
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name,"cs=f"=>\$cs,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"cf=i"=>\$cf,"sol=i{1,}"=>\@sol );
#
$show
#
DISPLAY AXES:0 0
contour
  delete contour plane 1
  # Add plane near "center-line" (bottom)
  add contour plane  0.00000e+00  1.00000e+00  0.00000e+00  .5  .01  0 
  # contour plane near middle
  add contour plane  0.00000e+00  1.00000e+00  0.00000e+00 0 2.5 0.
  plot the grid
    plot block boundaries 0
    toggle grid 3 0
    plot grid lines 0
    toggle boundary 1 1 4 0
    toggle boundary 0 1 0 0
    grid colour 4 BRASS
    grid colour 1 BRASS
    grid colour 2 BRASS
    set view:0 -0.0845921 0.108761 0 1.15331 0.939693 0.116978 -0.321394 0 0.939693 0.34202 0.34202 -0.321394 0.883022
  exit
exit
pause
# 
$cmd="#"; 
for( $i=0; $i<=$#sol; $i++ ){ $cmd .= "\n solution:$sol[$i]\n \$plotName = $name . \"$field$sol[$i].ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd

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