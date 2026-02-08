#
#  plotStuff plotComplexSolution.cmd -show=splitRingG8O4Freq50 -solution=1 -field=absur absui speed
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
if( $name eq "plot" ){ $name=$show; }
if( $field[0] eq "" ){ @field=( "ur", "ui", "urTotal", "uiTotal" ); }
# 
$show
#
derived types
#
  absoluteValue
    all
  #
  specify velocity components
    0 1 1
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






# 
#
#
$show="gaussianSquare.show"; $solution="-1"; $name="plot"; $field="Ey"; $emin=0; $emax=-1; $numFreq=1; $clines=0; $comp="v"; 
$res=1024; # hardcopy resolution
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave, "comp=s"=>\$comp,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"numFreq=i"=>\$numFreq,"clines=i"=>\$clines,\
       "emin=f"=>\$emin,"emax=f"=>\$emax,"numFreq=i"=>\$numFreq,"res=i"=>\$res,"field=s{1,}"=>\@field );
#
$show
if( $name eq "plot" ){ $name=$show; }
#
if( $comp eq "absv" ){ $cmd="derived types\n absoluteValue\n all\n   specify velocity components\n 0 1 0\n speed\n exit"; }else{ $cmd="#" }
$cmd
#
contour
  #plot:v0
  $compi = $comp . "0"; 
  plot:$compi
  if( $clines ==0 ){ $cmd="plot contour lines (toggle)"; }else{ $cmd="#"; }
  $cmd 
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  if( $emax > $emin ){ $cmd="min max $emin $emax"; }else{ $cmd="#"; }
  $cmd
exit
solution: $solution
pause
#
DISPLAY AXES:0 0
if( $res ne 1024 ){ $cmd="hardcopy vertical resolution:0 $res\n hardcopy horizontal resolution:0 $res\n line width scale factor:0 3"; }else{ $cmd="#"; }
$cmd 
# DISPLAY LABELS:0 0
# DISPLAY COLOUR BAR:0 0
# 
$cmd="#"; 
for( $i=0; $i<$numFreq; $i++ ){ $cmd .= "\n plot:$comp$i\n \$plotName = $name . \"$comp$i.ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd
exit


plot
$plotName = $name . "WaveHoltz.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0


#
# -- forcing:
contour
  plot:f0
  plot contour lines (toggle)
  plot boundaries (toggle)
  vertical scale factor 0.7
exit
x-r:0 65
$plotName = $name . "WaveHoltzForcing.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0




# -- OLD STUFF
Foreground colour:0 white
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 4
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
# 
plot
$plotName = $name . "EfieldNorm.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0