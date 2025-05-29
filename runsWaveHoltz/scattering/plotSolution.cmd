#
# 
#   plotStuff plotSolution.cmd -show=scatCylG4kx3.show -name=scatCylG4kx3 
#   plotStuff plotSolution.cmd -show=scatCylG4SuperGridkx3.show -name=scatCylG4SUperGridkx3 
#   plotStuff plotSolution.cmd -show=shapesG4kx4.show -name=shapesG4kx4
#
# plot abs(urTotal), abs(uiTotal)
#  plotStuff plotSolution.cmd -show=diskG32kx16.show -name=diskG32kx16 -field=absurTotal absuiTotal -res=2048
# -- split ring resonator:
# plotStuff plotSolution.cmd -show=splitRingScatO4G4kx4.show -name=splitRingScatO4G4kx4 -field=absurTotal absuiTotal 
#
$show="scatCylG4kx3.show"; $solution="-1"; $name="plot"; 
@field = (); 
$emin=0; $emax=-1; $numFreq=1; $clines=0; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
$res=1024; # hardcopy resolution
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"numFreq=i"=>\$numFreq,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"field=s{1,}"=>\@field,"res=i"=>\$res );
#
if( $field[0] eq "" ){ @field=( "ur", "ui", "urTotal", "uiTotal" ); }
# 
$show
#
derived types
#
  absoluteValue
    ur  (off)
    ui  (off)  
    urTotal 
    uiTotal 
  done
  specify velocity components
    2 3 0
  speed  
exit
#
plot:$field[0]
contour
  # plot:v0
  if( $clines ==0 ){ $cmd="plot contour lines (toggle)"; }else{ $cmd="#"; }
  $cmd 
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  if( $emax > $emin ){ $cmd="min max $emin $emax"; }else{ $cmd="#"; }
  # $cmd
  # plot:$field[1]
  $cmd
  plot:$field[0]
exit
DISPLAY AXES:0 0
x-:0
pause
if( $res ne 1024 ){ $cmd="hardcopy vertical resolution:0 $res\n hardcopy horizontal resolution:0 $res\n line width scale factor:0 3"; }else{ $cmd="#"; }
$cmd
$cmd="#"; 
for( $i=0; $i<=$#field; $i++ ){ $cmd .= "\n plot:$field[$i]\n \$plotName = $name . \"$field[$i].ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd



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
