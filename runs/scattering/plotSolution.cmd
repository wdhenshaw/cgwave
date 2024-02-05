# 
# 1x26 SPIE
#  plotStuff plotSolution.cmd -show=holesRad01G16N1M26SPIE.show -name=holesRad01G16N1M26SPIE -vmin=-.85 -vMax=.85 -start=11 -stride=5 -last=36
# 3x26 SPIE
#  plotStuff plotSolution.cmd -show=holesRad01G16N3M26SPIE.show -name=holesRad01G16N3M26SPIE -vmin=-.85 -vMax=.85 -start=11 -stride=5 -last=36
# 7x26 SPIE
#  plotStuff plotSolution.cmd -show=holesRad01G16N7M26SPIE.show -name=holesRad01G16N7M26SPIE -solution=36 -vmin=-.85 -vMax=.85 -start=11 -stride=5 -last=36
#  plotStuff plotSolution.cmd -show=holesRad01G16N7M26SPIE.show -name=holesRad01G16N7M26SPIE -field=absu -plotBoundary=0 -vmin=0 -vMax=1 -solution=36 -start=11 -stride=5 -last=36
# 7x26 EME
#  plotStuff plotSolution.cmd -show=holesRad01G16N7M26EME.show -name=holesRad01G16N7M26EME -solution=36 -vmin=-.85 -vMax=.85 -start=11 -stride=5 -last=36
#  plotStuff plotSolution.cmd -show=holesRad01G16N7M26EME.show -name=holesRad01G16N7M26EME -field=absu -plotBoundary=0 -vmin=0 -vMax=1 -solution=36 -start=11 -stride=5 -last=36
# 7x26 Offset SPIE
#  plotStuff plotSolution.cmd -show=holesRad01OffsetG16N7M26SPIE.show -name=holesRad01OffsetG16N7M26SPIE -solution=36 -vmin=-.85 -vMax=.85 -start=11 -stride=5 -last=36
#  plotStuff plotSolution.cmd -show=holesRad01OffsetG16N7M26SPIE.show -name=holesRad01OffsetG16N7M26SPIE -field=absu -plotBoundary=0 -vmin=0 -vMax=1 -solution=36 -start=11 -stride=5 -last=36
#
# 7x16 SPIE: 
#  plotStuff plotSolution.cmd -show=holesRad01G16N7M26SPIE.show -name=holesRad01G16N7M26SPIE -solution=36
# 7x16 EME
#  plotStuff plotSolution.cmd -show=holesRad01G16N7M26EME.show -name=holesRad01G16N7M26EME -solution=36
#
# MOVIE
#  plotStuff plotSolution.cmd -show=holesRad01G16N7M26SPIEMovie.show -movie=1 -movieName=holesRad01N7M26 -vmin=-.95 -vMax=.8
#  plotStuff plotSolution.cmd -show=holesRad01OffsetG16N7M26SPIEMovie.show -movie=1 -movieName=holesRad01OffsetN7M26 -vmin=-.95 -vMax=.8
#
# Thin tip grid
# plotStuff plotSolution.cmd -show=tipGridG64O2SPIEk80.show -name=tipGridG64O2SPIEk80 -view=1 -solution=21 -start=21 -last=21
# plotStuff plotSolution.cmd -show=tipGridG64O2SPIEk40.show -name=tipGridG64O2SPIEk40 -view=1 -solution=21 -start=21 -last=21
# plotStuff plotSolution.cmd -show=tipGridG64O2SPIEk20.show -name=tipGridG64O2SPIEk20 -view=1 -solution=21 -start=21 -last=21
# plotStuff plotSolution.cmd -show=tipGridG64O2SPIEk10.show -name=tipGridG64O2SPIEk10 -view=1 -solution=21 -start=21 -last=21
#
# plot abs(u)
# plotStuff plotSolution.cmd -show=tipGridG64O2SPIEk20.show -name=tipGridG64O2SPIEk20 -view=1 -field=absu -solution=21 -start=21 -last=21
#
# plotStuff plotSolution.cmd -show=tipGridG64O5SPIEk20.show -name=tipGridG64O5SPIEk20.show -view=1 -field=absu -solution=21 -start=13 -last=21 -stride=2
#
$show="gaussianSquare.show"; $solution="-1"; $name="plot"; $field="u"; $vmin=0; $vmax=-1; $clines=0; $movie=0; $movieName="movie";
$start=1; $stride=5; $last=5; # save solution at these time intervals
$view=0; $plotBoundary=1; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave, "clines=i"=>\$clines,\
      "start=i"=>\$start,"last=i"=>\$last,"stride=i"=>\$stride,"vmin=f"=>\$vmin,"vmax=f"=>\$vmax,"movie=i"=>\$movie,\
       "movieName=s"=>\$movieName,"view=i"=>\$view,"plotBoundary=i"=>\$plotBoundary, "field=s"=>\$field );
#
$show
#
#
derived types
absoluteValue
u  (off)
done
exit
plot:$field
contour
  if( $clines eq "0" ){ $cmd="plot contour lines (toggle)"; }else{ $cmd="#"; }
  $cmd
  # set view:0 0.0694864 -0.0362538 0 2.34381 1 0 0 0 1 0 0 0 1
  coarsening factor 1 (<0 : adaptive)
  vertical scale factor 0.
  if( $vmax > $vmin ){ $cmd="min max $vmin $vmax"; }else{ $cmd="#"; }
  $cmd
  if( $plotBoundary eq "0" ){ $cmd="plot boundaries (toggle)"; }else{ $cmd="#"; }
  $cmd 
exit
#
DISPLAY AXES:0 0
# Movie: 
$cmd="#"; 
if( $movie eq "1"){ \
  $cmd =   "bigger:0\n" \
         . "movie file name: $movieName\n" \
         . "solution: 101\n" \
         . "movie frames: 800\n" \
         . "DISPLAY COLOUR BAR:0 0\n" \
         . "save movie files 1\n" \
         . "pause\n" \
         . "show movie\n" \
         . "  "; \
}
$cmd
#
$cmd="#";
if( $view eq "1" ){ $cmd="set view:0 -0.113215 -0.0129243 0 2.39824 1 0 0 0 1 0 0 0 1"; }
$cmd
solution: $solution
#
#
pause
$cmd="#"; 
for( $i=$start; $i<=$last; $i=$i+$stride ){ $cmd .= "\n solution: $i\n \$plotName = \"$name$field$i.ps\"; \n hardcopy file name:0 \$plotName\n hardcopy save:0"; }
$cmd


pause
# for disk: 
set view:0 -0.1 -0.0533736 0 1.06284 1 0 0 0 1 0 0 0 1
#
DISPLAY AXES:0 0
# DISPLAY LABELS:0 0
#
$plotName = $name . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
plot:err
$plotName = $name . "Err.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0



DISPLAY COLOUR BAR:0 0
bigger
# 
plot
$plotName = $name . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0