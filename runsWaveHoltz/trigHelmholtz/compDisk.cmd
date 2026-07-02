#
#  Compute self-convergence rates:
# Usage:
#   comp compDisk
#
# -- order 4 -- 
interpolation width: 5
specify files (coarse to fine)
diskG2O4Freq30.show 
diskG4O4Freq30.show 
diskG8O4Freq30.show 
diskG16O4Freq30.show 
diskG32O4Freq30.show 
diskG64O4Freq30.show 
exit
choose a solution
 1
#  3 
compute errors

# -- order 2 -- 
interpolation width: 3
specify files (coarse to fine)
diskG4O2Freq30.show 
diskG8O2Freq30.show 
diskG16O2Freq30.show 
diskG32O2Freq30.show 
diskG64O2Freq30.show 
diskG128O2Freq30.show 
exit
choose a solution
 1
#  3 
compute errors



