# 
#   plotStuff plotRandom -show=randomCICOrder2G4.show -name=randomCICO2
#   plotStuff plotRandom -show=randomCICOrder4G4.show -name=randomCICO4
# 
#   plotStuff plotRandom -show=randomDiskOrder2G4.show -name=randomDiskO2
#   plotStuff plotRandom -show=randomDiskOrder4G4.show -name=randomDiskO4
# 
#   plotStuff plotRandom -show=randomPipeOrder2G4.show -name=randomPipeO2
#   plotStuff plotRandom -show=randomPipeOrder4G4.show -name=randomPipeO4
# 
#   plotStuff plotRandom -show=randomSphereOrder2G4.show -name=randomSphereO2
#   plotStuff plotRandom -show=randomSphereOrder4G4.show -name=randomSphereO4
#
#
$show="diffract16.show"; $min=0.; $max=-1.; $grids=1;  $umin=0.; $umax=-1.; $solution=1; $nd=2; 
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "min=f"=>\$min,"max=f"=>\$max, "umin=f"=>\$umin,"umax=f"=>\$umax,\
            "solution=i"=>\$solution,"nd=i"=>\$nd );
#
$show
# 
# --- save sequence info ----
plot sequence:solutionNorms
  uNorm
  add energy
  save results to a matlab file
    $matlab = $name . ".m"; 
    $matlab
exit
exit