#
# 
#   plotStuff plotSolution.cmd -show=diskG4O2Eigs.show -name=diskG4O2Eigs -sol=1 2 
#
$show="gaussianSquare.show"; $solution="-1"; $name="plot"; $field="phi"; $emin=0; $emax=-1; $numFreq=1; $clines=0; 
# $tSave=1; $stride=1; $start=0; 
$view=0; 
# $numToSave=5; # save solution at these time intervals
# List of solutions to plot and save 
@sol= ();  # this must be null for GetOptions to work, defaults are given below
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "stride=i"=>\$stride, "numToSave=i"=>\$numToSave,"start=i"=>\$start,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"view=i"=>\$view,"sol=i{1,}"=>\@sol );
#
if( $sol[0] eq "" ){ @sol=($solution); }
$show
#
#  *new way* Eigenvectors stored as "time-steps" May 1, 2023
solution: $solution
contour
  if( $clines ==0 ){ $cmd="plot contour lines (toggle)"; }else{ $cmd="#"; }
  $cmd 
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  if( $emax > $emin ){ $cmd="min max $emin $emax"; }else{ $cmd="#"; }
  $cmd
exit
$cmd="#";
if( $view eq 1 ){ $cmd="set view:0 -0.083521 -0.0027465 0 1.16326 1 0 0 0 1 0 0 0 1"; } # drakCornerRoom
$cmd
pause
#
DISPLAY AXES:0 0
# DISPLAY LABELS:0 0
# DISPLAY COLOUR BAR:0 0
# 
$cmd="#"; 
for( $i=0; $i<=$#sol; $i++ ){ $cmd .= "\n solution:$sol[$i]\n \$plotName = $name . \"$field$sol[$i].ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd


$cmd="#"; 
for( $i=$start; $i<$start+$numToSave; $i=$i+$stride ){ $cmd .= "\n plot:$field$i\n \$plotName = $name . \"$field$i.ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd


# *********** OLD WAY : eigenvectors stored as components *******
contour
  plot:$field$start
  if( $clines ==0 ){ $cmd="plot contour lines (toggle)"; }else{ $cmd="#"; }
  $cmd 
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  if( $emax > $emin ){ $cmd="min max $emin $emax"; }else{ $cmd="#"; }
  $cmd
exit
$cmd="#";
if( $view eq 1 ){ $cmd="set view:0 -0.083521 -0.0027465 0 1.16326 1 0 0 0 1 0 0 0 1"; } # drakCornerRoom
$cmd
solution: $solution
pause
#
DISPLAY AXES:0 0
# DISPLAY LABELS:0 0
# DISPLAY COLOUR BAR:0 0
# 
$cmd="#"; 
for( $i=$start; $i<$start+$numToSave; $i=$i+$stride ){ $cmd .= "\n plot:$field$i\n \$plotName = $name . \"$field$i.ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd