#
#  Compute self-convergence rates:
# Usage:
#     comp compDoubleEllipse
#
#  double ellipse: omega=10
# -----------------double ellipse: ORDER=2 omega=10  source at (-1.4,-5.6)------------------------------------
specify files (coarse to fine)
darkCornerRoomFreq10Np8G4Order2Deflate64.show
darkCornerRoomFreq10Np8G8Order2Deflate64.show
darkCornerRoomFreq10Np8G16Order2Deflate64.show
exit
choose a solution
 1
#  3 
compute errors

# -----------------double ellipse: ORDER=4 omega=10  source at (-1.4,-5.6)------------------------------------
specify files (coarse to fine)
darkCornerRoomFreq10Np8G2Order4Deflate64.show
darkCornerRoomFreq10Np8G4Order4Deflate64.show
darkCornerRoomFreq10Np8G8Order4Deflate64.show
exit
choose a solution
 1
#  3 
compute errors

# -------- double ellipse: omega=10  source at (1,1) -----------
specify files (coarse to fine)
darkCornerRoomFreq10Np8G2Order4X1Y1.show 
darkCornerRoomFreq10Np8G4Order4X1Y1.show 
darkCornerRoomFreq10Np8G8Order4X1Y1.show 
exit
choose a solution
 1
#  3 
compute errors

