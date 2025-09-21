#
#  Makefile to run regression tests for WaveHoltz
#
check: buildGrids
	@echo "\n ------ RUN WAVEHOLTZ CHECKS ------\n"
	./checkWaveHoltz.p -solver=cgwh -testsFile=testsWaveHoltz

check-quiet:  
	./checkWaveHoltz.p -verbose=0 -solver=cgwh -testsFile=testsWaveHoltz

replace:  
	@echo "\n ------ REPLACE WAVEHOLTZ CHECKS ------\n"
	./checkWaveHoltz.p -solver=cgwh -testsFile=testsWaveHoltz	-replace
	

# build grids that are used by the regression tests (and not already in Overture/sampleGrids)
buildGrids: grids/square64.order2.hdf grids/square32.order4.ng3.hdf grids/square32.order2.hdf grids/cice2.order2.hdf
	touch buildGrids
	rm -f ogen.cmd ogen.check ogen.log	

grids/square64.order2.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/squareArg -order=2 -nx=64
	mv square64.order2.hdf grids/

grids/square32.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/squareArg -order=4 -nx=32 -numGhost=3
	mv square32.order4.ng3.hdf grids/

grids/square32.order2.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/squareArg -order=2 -nx=32 
	mv square32.order2.hdf grids/


grids/cice2.order2.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/cicArg -factor=2 -interp=e -order=2
	mv cice2.order2.hdf grids/	