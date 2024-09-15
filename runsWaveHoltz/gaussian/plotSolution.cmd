#
# 
#   plotStuff plotSolution.cmd -show=gaussianSquare.show -name=gaussianSquareOmega -solution=1
#   plotStuff plotSolution.cmd -show=gaussianShapesOmega9p777.show -name=gaussianShapesOmega9p777 -solution=1
#   plotStuff plotSolution.cmd -show=gaussianDiskOmega8p1.show -name=gaussianDiskOmega8p1 -solution=1
#   plotStuff plotSolution.cmd -show=annulusOmega9p777.show -name=annulusOmega9p777 -solution=1
#
#   plotStuff plotSolution.cmd -show=shapesOmega66.show -name=shapesOmega66 -solution=1
#   plotStuff plotSolution.cmd -show=shapesOmega67.show -name=shapesOmega67 -solution=1
#   plotStuff plotSolution.cmd -show=shapesOmega68.show -name=shapesOmega68 -solution=1
#   plotStuff plotSolution.cmd -show=shapesOmega69.show -name=shapesOmega69 -solution=1
#   plotStuff plotSolution.cmd -show=shapesOmega70.show -name=shapesOmega70 -solution=1
#
# plotStuff plotSolution.cmd -show=shapesG4O4freq45Np80 -name=shapesG4O4freq45Np80 -solution=1
# plotStuff plotSolution.cmd -show=shapesG4O4freq44Np80 -name=shapesG4O4freq44Np80 -solution=1
# plotStuff plotSolution.cmd -show=shapesG4O4freq43Np80 -name=shapesG4O4freq43Np80 -solution=1
#
#
# plotStuff plotSolution.cmd -show=multiObjectsFreq10Np8G4Order2Deflate64 -name=multiObjectsFreq10Np8G4Order2Deflate64 -solution=1
# plotStuff plotSolution.cmd -show=multiObjectsFreq10Np8G8Order2Deflate64 -name=multiObjectsFreq10Np8G8Order2Deflate64 -solution=1
# plotStuff plotSolution.cmd -show=multiObjectsFreq10Np8G4Order4Deflate128.show -name=multiObjectsFreq10Np8G4Order4Deflate128 -solution=1
# Plot DHS: 
# plotStuff plotSolution.cmd -show=multiObjectsFreq10Np2G8Order4DHS.show -name=multiObjectsFreq10G8Order4DHS -solution=1
#
# Order=2:
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq5Np8G4Order2Deflate64 -name=darkCornerRoomFreq5Np8G4Order2Deflate64 -solution=1
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq10Np8Nit10NG4Order2Deflate64 -name=darkCornerRoomFreq10Np8Nit10NG4Order2Deflate64 -solution=1
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq10Np8G8Order2Deflate64 -name=darkCornerRoomFreq10Np8G8Order2Deflate64 -solution=1
# Plot DHS: 
#   plotStuff plotSolution.cmd -show=darkCornerRoomFreq10Np8G16Order2 -name=darkCornerRoomFreq10Np8G16Order2 -solution=1 -comp=uDHS
# Order=4 
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq10Np8G4Order4Deflate64 -name=darkCornerRoomFreq10Np8G4Order4Deflate64 -solution=1
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq10Np8G8Order4Deflate64 -name=darkCornerRoomFreq10Np8G8Order4Deflate64 -solution=1
# 
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq15Np8Nit10NG4Order4Deflate100 -name=darkCornerRoomFreq15Np8Nit10NG4Order4Deflate100 -solution=1
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq10Np8Nit10NG4Order4Deflate256 -name=darkCornerRoomFreq10Np8Nit10NG4Order4Deflate256 -solution=1
#
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq15Np8Nit10NG4Order4Deflate256 -name=darkCornerRoomFreq15Np8Nit10NG4Order4Deflate256 -solution=1
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq5Np8Nit10NG4Order4Deflat646 -name=darkCornerRoomFreq5Np8Nit10NG4Order4Deflate46 -solution=1
# 
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq30G16O4 -name=darkCornerRoomFreq30G16O4 -solution=1 -comp=absv -emin=0 -emax=1
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq30G8O4X3Y3 -name=darkCornerRoomFreq30G8O4X3Y3 -solution=1 -comp=absv -emin=0 -emax=1
# plotStuff plotSolution.cmd -show=darkCornerRoomFreq40G8O4X3Y3 -name=darkCornerRoomFreq40G8O4X3Y3 -solution=1 -comp=absv -emin=0 -emax=0.3
#
$show="gaussianSquare.show"; $solution="-1"; $name="plot"; $field="Ey"; $emin=0; $emax=-1; $numFreq=1; $clines=0; $comp="v"; 
$res=1024; # hardcopy resolution
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave, "comp=s"=>\$comp,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"numFreq=i"=>\$numFreq,"clines=i"=>\$clines,\
      "field=s"=>\$field,"emin=f"=>\$emin,"emax=f"=>\$emax,"numFreq=i"=>\$numFreq,"res=i"=>\$res );
#
$show
if( $name eq "plot" ){ $name=$show; }
#
if( $comp eq "absv" ){ $cmd="derived types\n absoluteValue\n v0  (off)\n done\n exit"; }else{ $cmd="#" }
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