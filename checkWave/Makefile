# Type "make check" to run the regression tests
#
# Type "make" in this folder to run all checks.
#
# Type checkWave.p <checkName.check> to run one check.
#    e.g. "checkWave.p diskEig6LCBC.check"
#
# To replace a check file use
#   checkWave.p -replace <checkName.check> 

check:  buildGrids buildCheckCheckFiles
	./checkWave.p -solver=cgWave
	
# 	@echo "\n ------NOW RUN CGWAVEHOLTZ CHECKS ------\n"
# 	./checkWave.p -solver=cgwh -testsFile=testsWaveHoltz

check-quiet: buildGrids buildCheckCheckFiles
	./checkWave.p -solver=cgWave -verbose=0

# build grids that are used by the regression tests (and not already in Overture/sampleGrids)
buildGrids: grids/square16.order2.hdf grids/square16.order4.ng3.hdf grids/square16.order6.ng4.hdf grids/square16.order8.ng5.hdf \
            grids/nonSquare16.order2.hdf grids/nonSquare16.order4.hdf grids/nonSquare16.order4.ng3.hdf grids/nonSquare16.order6.ng4.hdf grids/nonSquare16.order8.ng5.hdf \
            grids/cice2.order2.hdf grids/cice2.order4.ng3.hdf \
            grids/box1.order2.hdf grids/box1.order4.hdf grids/nonBox1.order2.hdf grids/nonBox1.order4.hdf \
            grids/sice2.order2.hdf grids/sice2.order4.ng3.hdf grids/sice2.order6.ng4.hdf grids/sice2.order8.ng5.hdf
	touch buildGrids
	rm -f ogen.cmd ogen.check ogen.log


# Build the program that processes the check files:
buildCheckCheckFiles:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) checkCheckFiles


grids/square16.order2.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=2 -nx=16
	mv square16.order2.hdf grids/

grids/square16.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=4 -nx=16 -numGhost=3
	mv square16.order4.ng3.hdf grids/

grids/square16.order6.ng4.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=6 -nx=16 -numGhost=4
	mv square16.order6.ng4.hdf grids/

grids/square16.order8.ng5.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=8 -nx=16 -numGhost=5	
	mv square16.order8.ng5.hdf grids/

grids/nonSquare16.order2.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=2 -n=16 
	mv nonSquare16.order2.hdf grids/

grids/nonSquare16.order4.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=4 -n=16 
	mv nonSquare16.order4.hdf grids/

grids/nonSquare16.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=4 -n=16 -numGhost=3	
	mv nonSquare16.order4.ng3.hdf grids/

grids/nonSquare16.order6.ng4.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=6 -n=16 -numGhost=4
	mv nonSquare16.order6.ng4.hdf grids/	

grids/nonSquare16.order8.ng5.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=8 -n=16 -numGhost=5
	mv nonSquare16.order8.ng5.hdf grids/	


grids/cice2.order2.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/cicArg -factor=2 -interp=e -order=2
	mv cice2.order2.hdf grids/

grids/cice2.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/cicArg -factor=2 -interp=e -order=4 -numGhost=3
	mv cice2.order4.ng3.hdf grids/

grids/box1.order2.hdf:
	$(Overture)/bin/ogen -noplot boxArg -order=2 -factor=1
	mv box1.order2.hdf grids/

grids/box1.order4.hdf:
	$(Overture)/bin/ogen -noplot boxArg -order=4 -factor=1
	mv box1.order4.hdf grids/


grids/nonBox1.order2.hdf:
	$(Overture)/bin/ogen -noplot boxArg -prefix=nonBox -order=2 -factor=1
	mv nonBox1.order2.hdf grids/

grids/nonBox1.order4.hdf:
	$(Overture)/bin/ogen -noplot boxArg -prefix=nonBox -order=4 -factor=1
	mv nonBox1.order4.hdf grids/

grids/sice2.order2.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/sicArg -factor=2 -interp=e -order=2
	mv sice2.order2.hdf grids/

grids/sice2.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/sicArg -factor=2 -interp=e -order=4 -numGhost=3	
	mv sice2.order4.ng3.hdf grids/

grids/sice2.order6.ng4.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/sicArg -factor=2 -interp=e -order=6 -numGhost=4
	mv sice2.order6.ng4.hdf grids/			

grids/sice2.order8.ng5.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/sicArg -factor=2 -interp=e -order=8 -numGhost=5
	mv sice2.order8.ng5.hdf grids/

#tube1.order4.hdf:
#	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/tube1.order4.cmd
#
#square16.order4.hdf:
#	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/square16.order4.cmd
#
#cici1.order4.hdf:
#	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/cicArg -factor=1 interp=i -order=4
#
#box16.order4.hdf:
#	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/box16.order4.cmd
#
#box1.order2.hdf:
#	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/boxArg.cmd -order=2 -factor=1 
#
#nonBox1.order2.hdf:
#	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/boxArg.cmd -prefix=nonBox -order=2 -factor=1
#
#box1.order4.ng3.hdf:
#	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/boxArg.cmd  -order=4 -factor=1 -numGhost=3
#
#nonBox1.order4.ng3.hdf:
#	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/boxArg.cmd -prefix=nonBox -order=4 -factor=1 -numGhost=3
#
#innerOuter4.order4.hdf:
#	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/innerOuter -factor=4 -order=4 -deltaRad=.5 -interp=e -name="innerOuter4.order4.hdf"
#
#twoBoxesInterfacee1.order4.hdf:
#	$(Overture)/bin/ogen noplot $(Overture)/sampleGrids/twoBoxesInterface -order=4 -interp=e -factor=1
#
#square20.order2.hdf:
#	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=2 -nx=20 
#
#nonSquare20.order2.hdf:
#	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=2 -n=20 
#
## nore ghost points for sosup: 
#square20.order4.ng3.hdf:
#	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=4 -nx=20 -numGhost=3
#
#nonSquare20.order4.ng3.hdf:
#	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=4 -n=20 -numGhost=3
#
#square20.order6.ng4.hdf:
#	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=6 -nx=20 -numGhost=4
#
#solidSphereInABoxe2.order2.hdf:
#	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/solidSphereInABox -order=2 -interp=e -factor=2
#
#solidSphereInABoxe2.order4.ng3.hdf:
#	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/solidSphereInABox -order=4 -interp=e -numGhost=3 -factor=2



checkCheckFiles:	
	@cd $(common); $(MAKE) -$(MAKEFLAGS) $(BINDIR)checkCheckFiles

# cgmx:
# 	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgmx

clean:
	rm -f *.hdf



