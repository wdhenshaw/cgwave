#
#  Compute self-convergence rates using the comp code:
# Usage:
#     comp compRates -caseName=[diskO4Freq60|diskO2Freq30] 
#
# Examples:
#   comp compRates -caseName=diskO4Freq60 
#   comp compRates -caseName=diskO2Freq30
#
#   comp compRates -caseName=rectangleO2Freq29
#   comp compRates -caseName=rectangleO4Freq29
#   comp compRates -caseName=rectangleO4Freq29SpacingRatio2
#   comp compRates -caseName=rectangleO4Freq29SpacingRatio1p5
#
#   comp compRates -caseName=rectangle4By1O2Freq25
#   comp compRates -caseName=rectangle4By1O2SpacingRatio1p5Freq25
#
#   comp compRates -caseName=rectangle4By1O4Freq25
#   comp compRates -caseName=rectangle4By1O4SpacingRatio1p5Freq25
# 
#   comp compRates -caseName=rectangle8By1O2Freq25
#   comp compRates -caseName=rectangle8By1O2SpacingRatio1p5Freq25
#
#   comp compRates -caseName=rectangle8By1O4Freq25
#   comp compRates -caseName=rectangle8By1O4SpacingRatio1p5Freq25
#   comp compRates -caseName=rectangle8By1O4SpacingRatio2p0Freq25
# 
$caseName="diskO4Freq60"; $iw=3; 
# get command line arguments
GetOptions( "caseName=s"=>\$caseName,"iw=i"=>\$iw );
#
if( $caseName eq "diskO4Freq60" ){ $iw=5; $showFiles = "diskG16O4Freq60.show\n diskG32O4Freq60.show \n diskG48O4Freq60.show \n diskG64O4Freq60.show ";  $comp="0\n 0\n 0\n 0"; }
# if( $caseName eq "diskO2Freq30" ){ $iw=3; $showFiles = "diskG32O2Freq30.show \n diskG64O2Freq30.show \n diskG128O2Freq30.show "; $comp="0\n 0\n 0"; }
if( $caseName eq "diskO2Freq30" ){ $iw=3; $showFiles = "diskG32O2Freq30.show \n diskG64O2Freq30.show \n diskG96O2Freq30.show \n diskG128O2Freq30.show ";  $comp="0\n 0\n 0\n 0"; }
if( $caseName eq "rectangleO2Freq29" ){ $iw=3; $showFiles = "rectangle16By1G8O2Freq29.show \n rectangle16By1G16O2Freq29.show \n rectangle16By1G32O2Freq29.show\n rectangle16By1G64O2Freq29.show ";  $comp="0\n 0\n 0\n 0"; }
if( $caseName eq "rectangleO2Freq29" ){ $iw=3; $showFiles = "rectangle16By1G8O2Freq29.show \n rectangle16By1G16O2Freq29.show \n rectangle16By1G32O2Freq29.show";  $comp="0\n 0\n 0\n 0"; }
#
if( $caseName eq "rectangleO4Freq29" ){ $iw=5; $showFiles = "rectangle16By1G4O4Freq29.show \nrectangle16By1G8O4Freq29.show \n rectangle16By1G16O4Freq29.show \n rectangle16By1G32O4Freq29.show";  $comp="0\n 0\n 0\n 0"; }
if( $caseName eq "rectangleO4Freq29SpacingRatio2" ){ $iw=5; $showFiles = "rectangle16By1G4O4SpacingRatio2Freq29.show \nrectangle16By1G8O4SpacingRatio2Freq29.show \n rectangle16By1G16O4SpacingRatio2Freq29.show \n rectangle16By1G32O4SpacingRatio2Freq29.show";  $comp="0\n 0\n 0\n 0"; }
if( $caseName eq "rectangleO4Freq29SpacingRatio1p5" ){ $iw=5; $showFiles = "rectangle16By1G4O4SpacingRatio1p5Freq29.show \nrectangle16By1G8O4SpacingRatio1p5Freq29.show \n rectangle16By1G16O4SpacingRatio1p5Freq29.show \n rectangle16By1G32O4SpacingRatio1p5Freq29.show";  $comp="0\n 0\n 0\n 0"; }
#
if( $caseName eq "rectangle4By1O2Freq25" ){ $iw=3; $showFiles = "rectangle4By1G16O2Freq25.show \nrectangle4By1G32O2Freq25.show \n rectangle4By1G64O2Freq25.show\n rectangle4By1G128O2Freq25.show";  $comp="0\n 0\n 0\n 0"; }
if( $caseName eq "rectangle4By1O2SpacingRatio1p5Freq25" ){ $iw=3; $showFiles = "rectangle4By1G16O2SpacingRatio1p5Freq25.show \nrectangle4By1G32O2SpacingRatio1p5Freq25.show \n rectangle4By1G64O2SpacingRatio1p5Freq25.show\n rectangle4By1G128O2SpacingRatio1p5Freq25.show";  $comp="0\n 0\n 0\n 0"; }
#
if( $caseName eq "rectangle4By1O4Freq25" ){ $iw=5; $showFiles = "rectangle4By1G4O4Freq25.show\n rectangle4By1G8O4Freq25.show \nrectangle4By1G16O4Freq25.show \nrectangle4By1G32O4Freq25.show \n rectangle4By1G64O4Freq25.show";  $comp="0\n 0\n 0\n 0\n 0"; }
if( $caseName eq "rectangle4By1O4SpacingRatio1p5Freq25" ){ $iw=5; $showFiles = "rectangle4By1G4O4SpacingRatio1p5Freq25.show \nrectangle4By1G8O4SpacingRatio1p5Freq25.show \nrectangle4By1G16O4SpacingRatio1p5Freq25.show \nrectangle4By1G32O4SpacingRatio1p5Freq25.show \n rectangle4By1G64O4SpacingRatio1p5Freq25.show";  $comp="0\n 0\n 0\n 0\n 0"; }
#
if( $caseName eq "rectangle8By1O2Freq25" ){ $iw=5; $showFiles = "rectangle8By1G4O2Freq25.show\n rectangle8By1G8O2Freq25.show \nrectangle8By1G16O2Freq25.show \nrectangle8By1G32O2Freq25.show \n rectangle8By1G64O2Freq25.show";  $comp="0\n 0\n 0\n 0\n 0"; }
# O4: 
if( $caseName eq "rectangle8By1O4Freq25" ){ $iw=5; $showFiles = "rectangle8By1G4O4Freq25.show\n rectangle8By1G8O4Freq25.show \nrectangle8By1G16O4Freq25.show \nrectangle8By1G32O4Freq25.show";  $comp="0\n 0\n 0\n 0"; }
if( $caseName eq "rectangle8By1O4SpacingRatio1p5Freq25" ){ $iw=5; $showFiles = "rectangle8By1G4O4SpacingRatio1p5Freq25.show\n rectangle8By1G8O4SpacingRatio1p5Freq25.show \nrectangle8By1G16O4SpacingRatio1p5Freq25.show \nrectangle8By1G32O4SpacingRatio1p5Freq25.show";  $comp="0\n 0\n 0\n 0"; }
if( $caseName eq "rectangle8By1O4SpacingRatio2p0Freq25" ){ $iw=5; $showFiles = "rectangle8By1G4O4SpacingRatio2p0Freq25.show\n rectangle8By1G8O4SpacingRatio2p0Freq25.show \nrectangle8By1G16O4SpacingRatio2p0Freq25.show \nrectangle8By1G32O4SpacingRatio2p0Freq25.show";  $comp="0\n 0\n 0\n 0"; }
#
interpolation width: $iw
specify files (coarse to fine)
$showFiles
exit
enter components to use per file
$comp
choose a solution
 1
#  3 
compute errors

