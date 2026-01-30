#
# 
#   plotStuff plotFlatTorus.cmd -show=flatTorusG8Freq20Ev256 -cf=6 -field=abs -sol=2 
#
#  NOTE: if $cf >0 then plot a wire frame, other-wise plot a brass surface 
#
$show="pipeG2O2ImpEig.show"; $name=""; $field="phi"; $cf=-1; 
$emin=0; $emax=-1; $numFreq=1; $clines=1; 
$cs=0; 
$numToSave=1; # save this many components
@sol= ();  # this must be null for GetOptions to work, defaults are given below
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name,"cs=f"=>\$cs,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"cf=i"=>\$cf,"sol=i{1,}"=>\@sol );
#
if( $name eq "" ){ $name=$show; }
#
$show
#
if( $field eq "abs" ){ $cmd="derived types\n absoluteValue\n phi  (off)\n done\n exit\n plot:absphi"; }else{ $cmd="#" }
$cmd
DISPLAY AXES:0 0
solution $sol[0]
#
contour 
  contour lines 0
  delete contour plane 0
  plot the grid
    toggle grid 0 0
    plot block boundaries 0
    # if $cf >0 then plot a wire frame, other-wise plot a brass surface 
    if( $cf ne "0" ){ $cmd="plot shaded surfaces (3D) 0\n coarsening factor $cf"; }else{ $cmd="plot grid lines 0\n grid colour 1 BRASS"; }
    $cmd
  exit this menu
exit
#
set view:0 -0.121224 0 0 0.82236 0.766044 0.219846 -0.604023 0 0.939693 0.34202 0.642788 -0.262003 0.719846
pause
# 
$cmd="#"; 
for( $i=0; $i<=$#sol; $i++ ){ $cmd .= "\n solution:$sol[$i]\n \$plotName = $name . \"$field$sol[$i].ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd


contour

  delete contour plane 1
  # Add plane near "center-line" (bottom)
  add contour plane  0.00000e+00  1.00000e+00  0.00000e+00  .5  .01  0 
  # contour plane near middle
  add contour plane  0.00000e+00  1.00000e+00  0.00000e+00 0 2.5 0.
  contour lines $clines
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