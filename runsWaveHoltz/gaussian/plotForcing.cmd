# For disk--
plot forcing
DISPLAY AXES:0 0
plot boundaries (toggle)
set view:0 -0.0848943 0.253021 0 1.04747 1 0 0 0 0.5 -0.866025 0 0.866025 0.5
hardcopy file name:0 diskHelmholtzGaussianForcing0.ps
hardcopy save:0
plot:f1
hardcopy file name:0 diskHelmholtzGaussianForcing1.ps
hardcopy save:0
plot:f2
hardcopy file name:0 diskHelmholtzGaussianForcing2.ps
hardcopy save:0
exit
#
reset:0
x-:0
x-:0
contour
plot:v0
hardcopy file name:0 diskHelmholtzGaussianSolution0.ps
hardcopy save:0
plot:v1
hardcopy file name:0 diskHelmholtzGaussianSolution1.ps
hardcopy save:0
plot:v2
hardcopy file name:0 diskHelmholtzGaussianSolution2.ps
hardcopy save:0