#
#   Compute CONVergence rates for CgWave
#
# FINISH ME  -- THIS VERSION FROM CgMx
#
common = $(CG)/common

wimp = ~/Dropbox/AMP/wimp

# The convergence rate script is called conv.p
all:  
	ln -sf $(common)/src/conv.p .

cgWave:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgWave

# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++++++++++++++ RUNS FOR WIMP PAPER ++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++

# ======== Sphere eigenmodes ======
sphereEigen: sphereEigenIME sphereEigenSPIE

# these take a while to run
sphereEigenIME:
	conv.p tests.conv -option=sphereEig -cg=sphereFixed -order=2 -tf=.4 -tp=.4 -upwind=1 -bc=d -cfl=100 -dtMax=2e-2 -amp=1 -mPhi=2 -mTheta=1 -mr=1 -ts=implicit -solveri=cgs -takeImplicitFirstStep=1 -ng=3
	cp sphereEigstdoneSidedsphereFixedImplicitUWCOrder2maxNorm.m $(wimp)/matlab/
	conv.p tests.conv -option=sphereEig -cg=sphereFixed -order=4 -tf=.4 -tp=.4 -upwind=1 -bc=d -bcApproach=cbc -cfl=100 -dtMax=2e-2 -amp=1 -mPhi=2 -mTheta=1 -mr=1 -ts=implicit -solveri=cgs -takeImplicitFirstStep=1 -ng=3
	cp sphereEigstdcbcsphereFixedImplicitUWCcbcOrder4maxNorm.m $(wimp)/matlab/	

sphereEigenSPIE:
	conv.p tests.conv -option=sphereEig -cg=sphereFixed -order=2 -tf=.4 -tp=.4 -upwind=1 -bc=d -cfl=.9 -amp=1 -mPhi=2 -mTheta=1 -mr=1 -ts=implicit -solveri=cgs -takeImplicitFirstStep=1 -rectangular=explicit -ng=3
	cp sphereEigstdoneSidedsphereFixedSPIEUWCOrder2maxNorm.m $(wimp)/matlab/
	conv.p tests.conv -option=sphereEig -cg=sphereFixed -order=4 -tf=.4 -tp=.4 -upwind=1 -bc=d -bcApproach=cbc -cfl=.9 -amp=1 -mPhi=2 -mTheta=1 -mr=1 -ts=implicit -solveri=cgs -takeImplicitFirstStep=1 -rectangular=explicit -ng=3
	cp sphereEigstdcbcsphereFixedSPIEUWCcbcOrder4maxNorm.m $(wimp)/matlab/	

# ======== 3D cylinder eigenmodes ======
cylEigen: cylEigenIME cylEigenSPIE

cylEigenIME:
	conv.p tests.conv -option=diskEig -cg=pipezFixed -order=2 -tf=.6 -tp=.6 -upwind=1 -bc=d -cfl=10 -dtMax=.02 -nTheta=1 -mr=1 -mz=1 -rad=.5 -ts=implicit -solveri=cgs -takeImplicitFirstStep=1 -ng=3
	cp diskEigstdoneSidedpipezFixedImplicitUWCOrder2maxNorm.m $(wimp)/matlab/
	conv.p tests.conv -option=diskEig -cg=pipezFixed -order=4 -tf=.6 -tp=.6 -upwind=1 -bc=d -bcApproach=cbc -cfl=10 -dtMax=.02 -nTheta=1 -mr=1 -mz=1 -rad=.5 -ts=implicit -solveri=cgs -takeImplicitFirstStep=1 -ng=3
	cp diskEigstdcbcpipezFixedImplicitUWCcbcOrder4maxNorm.m $(wimp)/matlab/	

cylEigenSPIE:
	conv.p tests.conv -option=diskEig -cg=pipezFixed -order=2 -tf=.6 -tp=.6 -upwind=1 -bc=d -cfl=10 -dtMax=.02 -nTheta=1 -mr=1 -mz=1 -rad=.5 -ts=implicit -solveri=cgs -takeImplicitFirstStep=1 -rectangular=explicit -ng=3
	cp diskEigstdoneSidedpipezFixedSPIEUWCOrder2maxNorm.m $(wimp)/matlab/	
	conv.p tests.conv -option=diskEig -cg=pipezFixed -order=4 -tf=.6 -tp=.6 -upwind=1 -bc=d -bcApproach=cbc -cfl=10 -dtMax=.02 -nTheta=1 -mr=1 -mz=1 -rad=.5 -ts=implicit -solveri=cgs -takeImplicitFirstStep=1 -rectangular=explicit -ng=3
	cp diskEigstdcbcpipezFixedSPIEUWCcbcOrder4maxNorm.m $(wimp)/matlab/	

# ======= disk eigenmode =========
diskEigen: diskEigenIME diskEigenSPIE

# disk IME-UW-PC Dirichlet and Neumann
diskEigenIME:
	conv.p tests.conv -option=diskEig -cg=sicFixed -order=2 -tf=.7 -tp=.7 -upwind=1 -bc=d -cfl=10 -dtMax=.02 -nTheta=2 -mr=1  -ts=implicit -takeImplicitFirstStep=1 -ng=4
	cp diskEigstdoneSidedsicFixedImplicitUWCOrder2maxNorm.m $(wimp)/matlab/
	conv.p tests.conv -option=diskEig -cg=sic -order=4 -tf=.7 -tp=.7 -upwind=1 -bc=d -bcApproach=cbc -cfl=10 -dtMax=.02 -nTheta=2 -mr=1 -ts=implicit -takeImplicitFirstStep=1 -ng=4
	cp diskEigstdcbcsicImplicitUWCcbcOrder4maxNorm.m $(wimp)/matlab/
	conv.p tests.conv -option=diskEig -cg=sicFixed -order=2 -tf=.7 -tp=.7 -upwind=1 -bc=n -cfl=10 -dtMax=.02 -nTheta=2 -mr=1  -ts=implicit -takeImplicitFirstStep=1 -ng=4
	cp diskEigstdoneSidedsicFixedImplicitUWCNeumannOrder2maxNorm.m $(wimp)/matlab
	conv.p tests.conv -option=diskEig -cg=sic -order=4 -tf=.7 -tp=.7 -upwind=1 -bc=n -bcApproach=cbc -cfl=10 -dtMax=.02 -nTheta=2 -mr=1  -ts=implicit -takeImplicitFirstStep=1 -ng=4
	cp diskEigstdcbcsicImplicitUWCcbcNeumannOrder4maxNorm.m $(wimp)/matlab/

diskEigenSPIE:
	conv.p tests.conv -option=diskEig -cg=sicFixed -order=2 -tf=.7 -tp=.7 -upwind=1 -bc=d -cfl=.9 -nTheta=2 -mr=1  -ts=implicit -takeImplicitFirstStep=1 -rectangular=explicit -ng=4
	cp diskEigstdoneSidedsicFixedSPIEUWCOrder2maxNorm.m $(wimp)/matlab/
	conv.p tests.conv -option=diskEig -cg=sicFixed -order=4 -tf=.7 -tp=.7 -upwind=1 -bc=d -bcApproach=cbc -cfl=.9 -nTheta=2 -mr=1 -ts=implicit -takeImplicitFirstStep=1 -rectangular=explicit -ng=4
	cp diskEigstdcbcsicFixedSPIEUWCcbcOrder4maxNorm.m $(wimp)/matlab/	


# ======= scattering from a cylinder =========
cylScat: cylScatIME cylScatSPIE

cylScatIME:
	conv.p tests.conv -option=cylScat -cg=cicFixed -order=2 -kx=2 -tf=.4 -tp=.4 -upwind=1 -bc=d -cfl=100 -dtMax=2e-2 -ts=implicit -ng=4
	cp cylScatstdoneSidedcicFixedImplicitUWCOrder2maxNorm.m $(wimp)/matlab
	conv.p tests.conv -option=cylScat -cg=cicFixed -order=4 -kx=2 -tf=.4 -tp=.4 -upwind=1 -bc=d -cfl=100 -dtMax=2e-2 -ts=implicit -ng=4
	cp cylScatstdoneSidedcicFixedImplicitUWCOrder4maxNorm.m $(wimp)/matlab/	

cylScatSPIE:
	conv.p tests.conv -option=cylScat -cg=cicFixed -order=2 -kx=2 -tf=.4 -tp=.4 -upwind=1 -bc=d -cfl=.9 -ts=implicit -rectangular=explicit -ng=4
	cp cylScatstdoneSidedcicFixedSPIEUWCOrder2maxNorm.m $(wimp)/matlab/	
	conv.p tests.conv -option=cylScat -cg=cicFixed -order=4 -kx=2 -tf=.4 -tp=.4 -upwind=1 -bc=d -cfl=.9 -ts=implicit -rectangular=explicit -ng=4
	cp cylScatstdoneSidedcicFixedSPIEUWCOrder4maxNorm.m $(wimp)/matlab/



