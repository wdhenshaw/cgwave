
cgwh= $(HOME)/research/cgwave/bin/cgwh
eigenWave= $(HOME)/Dropbox/eigenWave


# CASES FOR EigenWave PAPER

# SQUARE
# cgwh -noplot eig.cmd -g=square128.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=30 -numPeriods=1  -matlab=eigSquareO2ImpEig12 -eigenVectorFile=/home/henshw/runs/eig/square128O2EigsEv512.show  -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=12 -numEigs=24 -eigenSolver=default -go=k


# ---- DISK EIGS IN PAPER ---
diskFreq10:
	$(cgwh) -noplot eig.cmd -g=sice4.order2  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskG4O2Freq10 -eigenVectorFile=/home/henshw/runs/eig/sice4O2EigsEv256.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -freq=10 -numEigs=42 -go=k > diskG4O2Freq10.out
	$(cgwh) -noplot eig.cmd -g=sice4.order4  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskG4O4Freq10 -eigenVectorFile=/home/henshw/runs/eig/sice4O4EigsEv256.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -freq=10 -numEigs=42 -go=k > diskG4O4Freq10.out
	cp diskG4O2Freq10Krylov.m $(eigenWave)/results
	cp diskG4O2Freq10Krylov.tex $(eigenWave)/tables
	cp diskG4O4Freq10Krylov.m $(eigenWave)/results
	cp diskG4O4Freq10Krylov.tex $(eigenWave)/tables




# ----------- DOUBLE ELLIPSE -----
doubleEllipse: doubleEllipseG4O2 doubleEllipseG4O4
doubleEllipseG4O4: 
	setsid $(cgwh) -noplot eig.cmd -g=darkCornerRoomGride4.order4.hdf -nf=1 -adjustOmega=0 -ts=implicit -cfl=100000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=6 -eigenVectorFile=darkCornerRoomG4O4EigsEv1500.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=11 -numEigs=256 -eigenSolver=default -matlab=darkCornerRoomG4O4Eigs -show=darkCornerRoomG4O4Eig11.show -go=ks > darkCornerRoomG4O4Eigs11.out
	cp darkCornerRoomG4O4EigsKrylov.m $(eigenWave)/results/
	cp darkCornerRoomG4O4EigsKrylov.tex $(eigenWave)/tables

doubleEllipseG4O2: 
	setsid $(cgwh) -noplot eig.cmd -g=darkCornerRoomGride4.order2.hdf -nf=1 -adjustOmega=0 -ts=implicit -cfl=100000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=6 -eigenVectorFile=darkCornerRoomG4O2EigsEv1500.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=11 -numEigs=256 -eigenSolver=default -matlab=darkCornerRoomG4O2Eigs -show=darkCornerRoomG4O2Eig11.show -go=ks > darkCornerRoomG4O2Eigs11.out
	cp darkCornerRoomG4O2EigsKrylov.m $(eigenWave)/results/
	cp darkCornerRoomG4O2EigsKrylov.tex $(eigenWave)/tables

#  --------  VARY MIN STEPS PER PERIOD -----
# cgwh -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=/home/henshw/runs/eig/sice8O4EigsEv256.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits10 -minStepsPerPeriod=10 -go=og 


varyNitsDisk:
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits5  -minStepsPerPeriod=5  -go=k > disk8Freq15Nits5.out 
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits6  -minStepsPerPeriod=6  -go=k > disk8Freq15Nits6.out 
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits7  -minStepsPerPeriod=7  -go=k > disk8Freq15Nits7.out 
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits8  -minStepsPerPeriod=8  -go=k > disk8Freq15Nits8.out 
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits9  -minStepsPerPeriod=9  -go=k > disk8Freq15Nits9.out 
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits10 -minStepsPerPeriod=10 -go=k > disk8Freq15Nits10.out
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits11 -minStepsPerPeriod=11 -go=k > disk8Freq15Nits11.out
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits12 -minStepsPerPeriod=12 -go=k > disk8Freq15Nits12.out
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits13 -minStepsPerPeriod=13 -go=k > disk8Freq15Nits13.out
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits14 -minStepsPerPeriod=14 -go=k > disk8Freq15Nits14.out
	$(cgwh) -noplot eig.cmd -g=sice8.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=2 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=15 -numEigs=50 -eigenSolver=default -matlab=disk8Freq15Nits15 -minStepsPerPeriod=15 -go=k > disk8Freq15Nits15.out


varyNitsSquare:
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits5  -minStepsPerPeriod=5  -go=k 
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits6  -minStepsPerPeriod=6  -go=k 
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits7  -minStepsPerPeriod=7  -go=k 
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits8  -minStepsPerPeriod=8  -go=k 
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits9  -minStepsPerPeriod=9  -go=k 
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits10 -minStepsPerPeriod=10 -go=k 
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits11 -minStepsPerPeriod=11 -go=k 
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits12 -minStepsPerPeriod=12 -go=k 
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits13 -minStepsPerPeriod=13 -go=k 
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits14 -minStepsPerPeriod=14 -go=k 
	$(cgwh) -noplot eig.cmd -g=square64.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-14 -upwind=0 -imode=1 -maxIterations=60 -numPeriods=1 -eigenVectorFile=/home/henshw/runs/eig/squareG64O2Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=10 -numEigs=12 -eigenSolver=default -matlab=square64Freq10Nits15 -minStepsPerPeriod=15 -go=k 


# ------------- MULTIGRID FOR A SQUARE -- runs to show O(N) scaling
squareOrder2mg:
	$(cgwh) -noplot eig.cmd -g=square32.order2   -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO2Freq12G32   -go=k > squareO2Freq12G32out 
	$(cgwh) -noplot eig.cmd -g=square64.order2   -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO2Freq12G64   -go=k > squareO2Freq12G64.out 
	$(cgwh) -noplot eig.cmd -g=square128.order2  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO2Freq12G128  -go=k > squareO2Freq12G128.out 
	$(cgwh) -noplot eig.cmd -g=square256.order2  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO2Freq12G256  -go=k > squareO2Freq12G256.out 
	$(cgwh) -noplot eig.cmd -g=square512.order2  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO2Freq12G512  -go=k > squareO2Freq12G512.out 
	$(cgwh) -noplot eig.cmd -g=square1024.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO2Freq12G1024 -go=k > squareO2Freq12G1024.out 
	$(cgwh) -noplot eig.cmd -g=square2048.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO2Freq12G2048 -go=k > squareO2Freq12G2048.out 

squareOrder4mg:
	$(cgwh) -noplot eig.cmd -g=square32.order4   -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO4Freq12G32   -go=k > squareO4Freq12G32out 
	$(cgwh) -noplot eig.cmd -g=square64.order4   -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO4Freq12G64   -go=k > squareO4Freq12G64.out 
	$(cgwh) -noplot eig.cmd -g=square128.order4  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO4Freq12G128  -go=k > squareO4Freq12G128.out 
	$(cgwh) -noplot eig.cmd -g=square256.order4  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO4Freq12G256  -go=k > squareO4Freq12G256.out 
	$(cgwh) -noplot eig.cmd -g=square512.order4  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO4Freq12G512  -go=k > squareO4Freq12G512.out 
	$(cgwh) -noplot eig.cmd -g=square1024.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO4Freq12G1024 -go=k > squareO4Freq12G1024.out 
	$(cgwh) -noplot eig.cmd -g=square2048.order4 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=200 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenVectorFile=none -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -orderOfExtrapolation=-1 -freq=12 -numEigs=16  -matlab=squareO4Freq12G2048 -go=k > squareO4Freq12G2048.out 

# ------------- MULTIGRID FOR A DISK -- runs to show O(N) scaling



# -------------- DISK FIXED ORDER 2 ----------
diskFixedOrder2mg: diskFixed4 diskFixed8 diskFixed16 diskFixed32 diskFixed64 diskFixed128

diskFixed128:
	setsid $(cgwh) -noplot eig.cmd -g=sicFixede128.order2.ml6 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskFixedO2Freq6G128 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskFixedO2Freq6G128.out
	cp diskFixedO2Freq6G128Krylov.m $(eigenWave)/results

diskFixed64:
	setsid $(cgwh) -noplot eig.cmd -g=sicFixede64.order2.ml5  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskFixedO2Freq6G64  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskFixedO2Freq6G64.out
	cp diskFixedO2Freq6G64Krylov.m $(eigenWave)/results

diskFixed32:
	setsid $(cgwh) -noplot eig.cmd -g=sicFixede32.order2.ml5  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskFixedO2Freq6G32  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskFixedO2Freq6G32.out
	cp diskFixedO2Freq6G32Krylov.m $(eigenWave)/results

diskFixed16:
	setsid $(cgwh) -noplot eig.cmd -g=sicFixede16.order2.ml4  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskFixedO2Freq6G16  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskFixedO2Freq6G16.out
	cp diskFixedO2Freq6G16Krylov.m $(eigenWave)/results

diskFixed8:
	setsid $(cgwh) -noplot eig.cmd -g=sicFixede8.order2.ml4   -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskFixedO2Freq6G8   -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskFixedO2Freq6G8.out
	cp diskFixedO2Freq6G8Krylov.m $(eigenWave)/results

diskFixed4:
	setsid $(cgwh) -noplot eig.cmd -g=sicFixede4.order2.ml3   -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskFixedO2Freq6G4   -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskFixedO2Freq6G4.out
	cp diskFixedO2Freq6G4Krylov.m $(eigenWave)/results


# -------------- DISK ORDER 2 ----------
diskOrder2mg: disk4 disk8 disk16 disk32 disk64 disk128

# G128: -- more levels
disk128ml6:
	setsid $(cgwh) -noplot eig.cmd -g=sice128.order2.ml6 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO2Freq6G128 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO2Freq6G128ml6.out

# G128: 6.7M points	
disk128:
	setsid $(cgwh) -noplot eig.cmd -g=sice128.order2.ml5 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO2Freq6G128 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO2Freq6G128.out

# G64: 1.7M points 
disk64ml6:
	setsid $(cgwh) -noplot eig.cmd -g=sice64.order2.ml6  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO2Freq6G64  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO2Freq6G64ml6.out
disk64:
	setsid $(cgwh) -noplot eig.cmd -g=sice64.order2.ml5  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO2Freq6G64  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO2Freq6G64.out

# G32: 424K points
# ml=6 
disk32ml6:
	setsid $(cgwh) -noplot eig.cmd -g=sice32.order2.ml6  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO2Freq6G32  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO2Freq6G32ml6.out

#ml=5  
disk32ml5:
	setsid $(cgwh) -noplot eig.cmd -g=sice32.order2.ml5  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO2Freq6G32  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO2Freq6G32ml5.out

# ml=4
disk32:
	setsid $(cgwh) -noplot eig.cmd -g=sice32.order2.ml4  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO2Freq6G32  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO2Freq6G32.out

# G16 : 110K points 
disk16:
	setsid $(cgwh) -noplot eig.cmd -g=sice16.order2.ml4  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO2Freq6G16  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO2Freq6G16.out

# G8 
disk8:
	setsid $(cgwh) -noplot eig.cmd -g=sice8.order2.ml3   -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO2Freq6G8   -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO2Freq6G8.out

# G4
disk4:
	setsid $(cgwh) -noplot eig.cmd -g=sice4.order2.ml3   -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO2Freq6G4   -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO2Freq6G4.out

# --------- DISK ORDER 4 -------------
diskOrder4mg: diskO4G4 diskO4G8 diskO4G16 diskO4G32 diskO4G64 diskO4G128

# G128: 6.7M points 
diskO4G128:
	setsid $(cgwh) -noplot eig.cmd -g=sice128.order4.ml5 -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO4Freq6G128 -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO4Freq6G128.out

# G64: 1.7M points 
diskO4G64:
	setsid $(cgwh) -noplot eig.cmd -g=sice64.order4.ml5  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO4Freq6G64  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO4Freq6G64.out

# G32: 424K points 
diskO4G32:
	setsid $(cgwh) -noplot eig.cmd -g=sice32.order4.ml4  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO4Freq6G32  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO4Freq6G32.out

# G16 : 110K points 
diskO4G16:
	setsid $(cgwh) -noplot eig.cmd -g=sice16.order4.ml4  -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO4Freq6G16  -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO4Freq6G16.out

# G8 
diskO4G8:
	setsid $(cgwh) -noplot eig.cmd -g=sice8.order4.ml3   -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO4Freq6G8   -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO4Freq6G8.out

# G4
diskO4G4:
	setsid $(cgwh) -noplot eig.cmd -g=sice4.order4.ml2   -nf=1 -adjustOmega=0 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskO4Freq6G4   -eigenVectorFile=none -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -eigenSolver=default -show=junk.show -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -freq=6 -numEigs=16 -go=k > diskO4Freq6G4.out

