#
#  Compute self-convergence rates using the comp code:
# Usage:
#     comp compRates -caseName=[diskO4Freq60|diskO2Freq30] 
#
# Examples:
#   comp compRates -caseName=diskO4Freq60 
#   comp compRates -caseName=diskO2Freq30
# 
$caseName="diskO4Freq60"; $iw=3; 
# get command line arguments
GetOptions( "caseName=s"=>\$caseName,"iw=i"=>\$iw );
#
if( $caseName eq "diskO4Freq60" ){ $iw=5; $showFiles = "diskG16O4Freq60.show\n diskG32O4Freq60.show \n diskG48O4Freq60.show \n diskG64O4Freq60.show ";  $comp="0\n 0\n 0\n 0"; }
# if( $caseName eq "diskO2Freq30" ){ $iw=3; $showFiles = "diskG32O2Freq30.show \n diskG64O2Freq30.show \n diskG128O2Freq30.show "; $comp="0\n 0\n 0"; }
if( $caseName eq "diskO2Freq30" ){ $iw=3; $showFiles = "diskG32O2Freq30.show \n diskG64O2Freq30.show \n diskG96O2Freq30.show \n diskG128O2Freq30.show ";  $comp="0\n 0\n 0\n 0"; }
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

