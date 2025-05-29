#
#  Compute self-convergence rates:
# Usage:
#     comp compKnives
#
# 
# -- eight knivess order 4 -- 
interpolation width: 5
specify files (coarse to fine)
diskG8O4Freq40.show 
diskG16O4Freq40.show 
diskG32O4Freq40.show 
# eightKnivesGrid256O4Freq60.show
exit
choose a solution
 1
#  3 
compute errors


# -- eight knivess order 4 -- 
interpolation width: 5
specify files (coarse to fine)
eightKnivesGrid32O4Freq60.show
eightKnivesGrid64O4Freq60.show
eightKnivesGrid128O4Freq60.show
# eightKnivesGrid256O4Freq60.show
exit
choose a solution
 1
#  3 
compute errors


# -- eight knvies order 2 -- not bad
interpolation width: 3
specify files (coarse to fine)
eightKnivesGrid64O2Freq40.show
eightKnivesGrid128O2Freq40.show
eightKnivesGrid256O2Freq40.show
exit
choose a solution
 1
#  3 
compute errors



interpolation width: 5
specify files (coarse to fine)
sixteenKnivesGrid64O4Freq40.show
sixteenKnivesGrid128O4Freq40.show
sixteenKnivesGrid256O4Freq40.show
exit
choose a solution
 1
#  3 
compute errors



interpolation width: 3
specify files (coarse to fine)
# sixteenKnivesGrid64O2Freq40.show
sixteenKnivesGrid128O2Freq40.show
sixteenKnivesGrid256O2Freq40.show
sixteenKnivesGrid400O2Freq40.show
exit
choose a solution
 1
#  3 
compute errors



interpolation width: 5
specify files (coarse to fine)
sixteenKnivesGrid32O4Freq175.show
sixteenKnivesGride64O4Freq175.show
sixteenKnivesGride128O4Freq175.show
exit
choose a solution
 1
#  3 
compute errors

# -----------------double ellipse: ORDER=2 omega=80  source at (0,.3)------------------------------------
specify files (coarse to fine)
# sixteenKnivesGrid64O2Freq80.show
sixteenKnivesGrid128O2Freq80.show
sixteenKnivesGrid256O2Freq80.show
sixteenKnivesGrid400O2Freq80.show
exit
choose a solution
 1
#  3 
compute errors

