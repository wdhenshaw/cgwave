#
#  plotStuff plotFlattenedTorusSolution.cmd -show=flattenedTorusG4O4kx3 -field=uTotalNorm
#  plotStuff plotFlattenedTorusSolution.cmd -show=flattenedTorusG8O4kx5 -field=absurTotal absuiTotal uTotalNorm
#  plotStuff plotFlattenedTorusSolution.cmd -show=flattenedTorusG8O4kx5 -name=flattenedTorusG8O4kx5View2 -view=2 -field=absurTotal absuiTotal uTotalNorm
#
$show="flattenedTorusG4O4kx3"; $solution="-1"; $name=""; 
@field = (); 
$emin=0; $emax=-1; $numFreq=1; $clines=0; 
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
$res=1024; # hardcopy resolution
$view=1;   # default view 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"numFreq=i"=>\$numFreq,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"field=s{1,}"=>\@field,"res=i"=>\$res,"view=i"=>\$view );
#
if( $name eq "" ){ $name=$show; }
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
  # 2-norm of [urTotal,uiTotal]
  twoNorm
   uTotalNorm
    2
    2 3 
exit
# -- plot the pill box grid 
grid
  toggle grid 0 0
  grid colour 1 BRASS
  plot block boundaries 0
  plot grid lines 0
 exit this menu
#
plot:$field[0]
contour
  # plot:v0
  delete contour plane 0
  if( $clines ==0 ){ $cmd="plot contour lines (toggle)"; }else{ $cmd="#"; }
  $cmd 
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  # coarsening factor 1 (<0 : adaptive)
  # vertical scale factor 0.
  if( $emax > $emin ){ $cmd="min max $emin $emax"; }else{ $cmd="#"; }
  # $cmd
  # plot:$field[1]
  $cmd
  plot:$field[0]
exit
DISPLAY AXES:0 0
$cmd="set view:0 -0.1 -0.0120846 0 1.02477 0.939693 0.116978 -0.321394 0 0.939693 0.34202 0.34202 -0.321394 0.883022";
if( $view eq 2 ){ $cmd = "set view:0 -0.121567 -0.00686625 0 0.904596 0.34202 0.163176 -0.925417 0 0.984808 0.173648 0.939693 -0.0593912 0.336824"; }
$cmd
pause
if( $res ne 1024 ){ $cmd="hardcopy vertical resolution:0 $res\n hardcopy horizontal resolution:0 $res\n line width scale factor:0 3"; }else{ $cmd="#"; }
$cmd
$cmd="#"; 
for( $i=0; $i<=$#field; $i++ ){ $cmd .= "\n plot:$field[$i]\n \$plotName = $name . \"$field[$i].ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd


