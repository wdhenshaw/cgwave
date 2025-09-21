#
#   plotStuff plotSolution.cmd -show=diskG4O2Eigs.show -name=diskG4O2Eigs -sol=1 2 
# RPI: 
#   plotStuff plotSolution.cmd -show=rpiG4O2Freq15.show -name=rpiG4O2 -field=abs -sol=10 11
#   plotStuff plotSolution.cmd -show=rpiG16O2Freq40.show -name=rpiG16O2 -field=abs -sol=10 11
# Movie:
#   plotStuff plotSolution.cmd -show=sixteenKnivesGridG128O4Freq300Ev1024.show -name=sixteenKnivesEigsHighFreq -field=abs -plotMovie=1 -numMovie=20 -solution=1200
#   ppm2mpeg sixteenKnivesEigsHighFreq 1 100
#
$show="gaussianSquare.show"; $solution="-1"; $name="plot"; $field="phi"; $emin=0; $emax=-1; $numFreq=1; $clines=0; 
$plotMovie=0; $numMovie=100; 
$numRepeat=20; # repeat each frame this many times 
# $tSave=1; $stride=1; $start=0; 
$view=0; 
# $numToSave=5; # save solution at these time intervals
# List of solutions to plot and save 
@sol= ();  # this must be null for GetOptions to work, defaults are given below
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "stride=i"=>\$stride, "numToSave=i"=>\$numToSave,"start=i"=>\$start,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"view=i"=>\$view,"sol=i{1,}"=>\@sol,\
      "plotMovie=i"=>\$plotMovie, "numMovie=i"=>\$numMovie );
#
if( $sol[0] eq "" ){ @sol=($solution); }
$show
#
if( $field eq "abs" ){ $cmd="derived types\n absoluteValue\n phi  (off)\n done\n exit\n plot:absphi"; }else{ $cmd="#" }
$cmd
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
# Movie of eigenvectors:
# hardcopy format:0 ppm
# DISPLAY COLOUR BAR:0 0
# set view:0 0.00241692 0 0 1.28494 1 0 0 0 1 0 0 0 1
# $cmd="#"; 
# #  -- plot movie --
# #
# for( $i=1; $i<=$numMovie*$numRepeat; $i++ ){ $cmd .= "\n \$plotName = \"$name$i.ppm\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; if( ($i % $numRepeat) == $numRepeat-1 ){ $cmd .= "\n next";}  }
# if( $plotMovie eq 0 ) { $cmd="#"; }
# $cmd
# 
# -- save indicated solutions
$cmd="#"; 
for( $i=0; $i<=$#sol; $i++ ){ $cmd .= "\n solution:$sol[$i]\n \$plotName = $name . \"$field$sol[$i].ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
if( $plotMovie eq 1 ) { $cmd="#"; }
$cmd




hardcopy file name:0 knife1.ppm
hardcopy save:0
next
hardcopy file name:0 knife2.ppm
hardcopy save:0


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