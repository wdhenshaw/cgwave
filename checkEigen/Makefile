#
#  Makefile to run regression tests for EigenWave
#
check:  buildGrids
	@echo "\n ------ RUN EIGENWAVE CHECKS ------\n"
	./checkEigen.p -solver=cgwh -testsFile=testsEigenWave

check-quiet:  
	./checkEigen.p -verbose=0 -solver=cgwh -testsFile=testsEigenWave

replace:  
	@echo "\n ------ REPLACE EIGENWAVE CHECKS ------\n"
	./checkEigen.p -solver=cgwh -testsFile=testsEigenWave	-replace
	
# build grids that are used by the regression tests (and not already in Overture/sampleGrids)
buildGrids:  grids/square64.order2.hdf grids/square32.order4.hdf grids/sice4.order2.hdf grids/box1.order2.hdf
	touch buildGrids
	rm -f ogen.cmd ogen.check ogen.log

grids/square64.order2.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/squareArg -order=2 -nx=64
	mv square64.order2.hdf grids/

grids/square32.order4.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/squareArg -order=4 -nx=32 
	mv square32.order4.hdf grids/

grids/box1.order2.hdf:
	$(Overture)/bin/ogen -noplot boxArg -order=2 -factor=1
	mv box1.order2.hdf grids/

grids/sice4.order2.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/sicArg -factor=4 -interp=e -order=2
	mv sice4.order2.hdf grids/