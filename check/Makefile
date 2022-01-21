# Type "make check" to run the regression tests
# MASTER VERSION of check.p and checkCheckFiles.C is in cg/common/src --> move to Overture! )

check:  buildGrids buildCheckCheckFiles
	./check.p -solver=cgWave
	@echo "\n ------NOW RUN CGWAVEHOLTZ CHECKS ------\n"
	./check.p -solver=cgwh -testsFile=testsWaveHoltz

check-quiet: buildGrids buildCheckCheckFiles
	./check.p -solver=cgWave -verbose=0
	./check.p -solver=cgwh -testsFile=testsWaveHoltz -verbose=0

# build grids that are used by the regression tests (and not already in Overture/sampleGrids)
buildGrids: square16.order2.hdf square16.order4.ng3.hdf nonSquare16.order2.hdf nonSquare16.order4.hdf cice2.order2.hdf cice2.order4.ng3.hdf \
            box1.order2.hdf nonBox1.order2.hdf nonBox1.order4.hdf sice2.order2.hdf sice4.order4.ng3.hdf
	touch buildGrids


# Build the program that processes the check files:
buildCheckCheckFiles:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) checkCheckFiles


square16.order2.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=2 -nx=16

square16.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/squareArg -order=4 -nx=16 -numGhost=3

nonSquare16.order2.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=2 -n=16 

nonSquare16.order4.hdf:
	$(Overture)/bin/ogen -noplot  $(Overture)/sampleGrids/nonSquare -order=4 -n=16 

cice2.order2.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/cicArg -factor=2 -interp=e -order=2

cice2.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/cicArg -factor=2 -interp=e -order=4 -numGhost=3

box1.order2.hdf:
	$(Overture)/bin/ogen -noplot boxArg -order=2 -factor=1

nonBox1.order2.hdf:
	$(Overture)/bin/ogen -noplot boxArg -prefix=nonBox -order=2 -factor=1

nonBox1.order4.hdf:
	$(Overture)/bin/ogen -noplot boxArg -prefix=nonBox -order=4 -factor=1

sice2.order2.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/sicArg -factor=2 -interp=e -order=2

sice4.order4.ng3.hdf:
	$(Overture)/bin/ogen -noplot $(Overture)/sampleGrids/sicArg -factor=2 -interp=e -order=4 -numGhost=3			

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

cgmx:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgmx

clean:
	rm -f *.hdf



