cgwh= $(HOME)/research/cgwave/bin/cgwh
WAHO= $(HOME)/Dropbox/AMP/waho


# ============= DOUBLE ELLIPSE GRID CONVERGENCE  =============

doubleEllipseGridConvergence: doubleEllipseG2Freq10 doubleEllipseG4Freq10 doubleEllipseG8Freq10

doubleEllipseG2Freq10:	
	setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride2.order4.hdf -ts=implicit -cfl=1000 -x0=1.0 -y0=1.0 -nf=1 -amp=-50 -beta=3 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG2O4EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G2Order4X1Y1 -show=darkCornerRoomFreq10Np8G2Order4X1Y1.show -go=dfks > darkCornerRoomFreq10Np8G2O4X1Y1.out
doubleEllipseG4Freq10:	
	setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride4.order4.hdf -ts=implicit -cfl=1000 -x0=1.0 -y0=1.0 -nf=1 -amp=-50 -beta=3 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG4O4EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G4Order4X1Y1 -show=darkCornerRoomFreq10Np8G4Order4X1Y1.show -go=dfks > darkCornerRoomFreq10Np8G4O4X1Y1.out
doubleEllipseG8Freq10:	
	setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride8.order4.hdf -ts=implicit -cfl=1000 -x0=1.0 -y0=1.0 -nf=1 -amp=-50 -beta=3 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG8O4EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G8Order4X1Y1 -show=darkCornerRoomFreq10Np8G8Order4X1Y1.show -go=dfks > darkCornerRoomFreq10Np8G8O4X1Y1.out
doubleEllipseG16Freq10:	
	setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride16.order4.hdf -ts=implicit -cfl=1000 -x0=1.0 -y0=1.0 -nf=1 -amp=-50 -beta=3 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG16O4EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G8Order4X1Y1 -show=darkCornerRoomFreq10Np8G16Order4X1Y1.show -go=dfks > darkCornerRoomFreq10Np8G16O4X1Y1.out


# ============= DISK WITH MULTIGRID =============

diskmg: disk4O4mg disk8O4mg disk16O4mg disk32O4mg disk64O4mg

disk64O4mg:
	setsid $(cgwh) -noplot gaussian.cmd -g=sice64.order4.ml5 -x0=-0.2 -y0=-0.3 -nf=1 -freq=20 -amp=-400 -beta=10 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=2 -minStepsPerPeriod=10 -cfl=10000 -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -solverh=yale -restart=100 -iluh=150 -maxith=10000 -rtolh=1.0e-10 -atolh=1.0e-10 -eigenVectorFile=/home/henshw/runs/eig/sice64Order4ml5EigsEv256.show -deflateWaveHoltz=1 -numToDeflate=32 -eigTol=1.0e-6 -matlab=sice64Freq20Deflate32 -show=diskFreq20G64O4.show -go=fks > sice64Freq20Deflate32.out


disk32O4mg:
	setsid $(cgwh) -noplot gaussian.cmd -g=sice32.order4.ml4 -x0=-0.2 -y0=-0.3 -nf=1 -freq=20 -amp=-400 -beta=10 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=2 -minStepsPerPeriod=10 -cfl=10000 -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -solverh=yale -restart=100 -iluh=150 -maxith=10000 -rtolh=1.0e-10 -atolh=1.0e-10 -eigenVectorFile=/home/henshw/runs/eig/sice32Order4ml4EigsEv256.show -deflateWaveHoltz=1 -numToDeflate=32 -eigTol=1.0e-6 -matlab=sice32Freq20Deflate32 -show=diskFreq20G32O4.show -go=dfks > sice32Freq20Deflate32.out


disk16O4mg:
	setsid $(cgwh) -noplot gaussian.cmd -g=sice16.order4.ml4 -x0=-0.2 -y0=-0.3 -nf=1 -freq=20 -amp=-400 -beta=10 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=2 -minStepsPerPeriod=10 -cfl=10000 -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -solverh=yale -restart=100 -iluh=150 -maxith=10000 -rtolh=1.0e-10 -atolh=1.0e-10 -eigenVectorFile=/home/henshw/runs/eig/sice16Order4ml4EigsEv256.show -deflateWaveHoltz=1 -numToDeflate=32 -eigTol=1.0e-6 -matlab=sice16Freq20Deflate32 -show=diskFreq20G16O4.show -go=dfks > sice16Freq20Deflate32.out

disk8O4mg:
	setsid $(cgwh) -noplot gaussian.cmd -g=sice8.order4.ml3  -x0=-0.2 -y0=-0.3 -nf=1 -freq=20 -amp=-400 -beta=10 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=2 -minStepsPerPeriod=10 -cfl=10000 -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -solverh=yale -restart=100 -iluh=150 -maxith=10000 -rtolh=1.0e-10 -atolh=1.0e-10 -eigenVectorFile=/home/henshw/runs/eig/sice8Order4ml3EigsEv256.show  -deflateWaveHoltz=1 -numToDeflate=32 -eigTol=1.0e-6 -matlab=sice8Freq20Deflate32  -show=diskFreq20G8O4.show -go=dfks > sice8Freq20Deflate32.out

disk4O4mg:
	setsid $(cgwh) -noplot gaussian.cmd -g=sice4.order4.ml2  -x0=-0.2 -y0=-0.3 -nf=1 -freq=20 -amp=-400 -beta=10 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=2 -minStepsPerPeriod=10 -cfl=10000 -solveri=multigrid -rtoli=1e-10 -atoli=1e-10 -solverh=yale -restart=100 -iluh=150 -maxith=10000 -rtolh=1.0e-10 -atolh=1.0e-10 -eigenVectorFile=/home/henshw/runs/eig/sice4Order4ml2EigsEv256.show  -deflateWaveHoltz=1 -numToDeflate=32 -eigTol=1.0e-6 -matlab=sice4Freq20Deflate32  -show=diskFreq20G4O4.show -go=dfks > sice4Freq20Deflate32.out

# ============ SQUARE WITH MULTIGRID ===============

# NITS=10 NP=2 CPU=71 

squaremg256: square256mg square256dhs
squaremg512: square512mg square512dhs
squaremg1024: square1024mg square1024dhs
squaremg2048: square2048mg square2048dhs

square256mg:
	setsid $(cgwh) -noplot gaussian.cmd -g=square256.order2 -x0=0.4 -y0=0.4 -nf=1 -freq=11 -amp=100 -beta=20 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-8 -upwind=0 -imode=1 -maxIterations=100 -numPeriods=2 -cfl=10000 -solverh=gmres -rtolh=1.0e-9 -atolh=1.0e-9 -maxith=50000 -solveri=multigrid -rtoli=1e-8 -atoli=1e-8 -go=k > square256MultiGridKrylov.out

# Direct Helmholtz solve: 
square256dhs:
	setsid $(cgwh) -noplot gaussian.cmd -g=square256.order2 -x0=0.4 -y0=0.4 -nf=1 -freq=11 -amp=100 -beta=20 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-8 -upwind=0 -imode=1 -maxIterations=100 -numPeriods=2 -cfl=10000 -solverh=gmres -rtolh=1.0e-9 -atolh=1.0e-9 -maxith=50000 -solveri=multigrid -rtoli=1e-8 -atoli=1e-8 -go=d > square256DHS.out


square512mg:
	setsid $(cgwh) -noplot gaussian.cmd -g=square512.order2 -x0=0.4 -y0=0.4 -nf=1 -freq=11 -amp=100 -beta=20 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-8 -upwind=0 -imode=1 -maxIterations=100 -numPeriods=2 -cfl=10000 -solverh=gmres -rtolh=1.0e-9 -atolh=1.0e-9 -maxith=50000 -solveri=multigrid -rtoli=1e-8 -atoli=1e-8 -go=k > square5124MultiGridKrylov.out

square512dhs:
	setsid $(cgwh) -noplot gaussian.cmd -g=square512.order2 -x0=0.4 -y0=0.4 -nf=1 -freq=11 -amp=100 -beta=20 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-8 -upwind=0 -imode=1 -maxIterations=100 -numPeriods=2 -cfl=10000 -solverh=gmres -rtolh=1.0e-9 -atolh=1.0e-9 -maxith=50000 -solveri=multigrid -rtoli=1e-8 -atoli=1e-8 -go=d > square5124DHS.out


square1024mg:
	setsid $(cgwh) -noplot gaussian.cmd -g=square1024.order2 -x0=0.4 -y0=0.4 -nf=1 -freq=11 -amp=100 -beta=20 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-8 -upwind=0 -imode=1 -maxIterations=100 -numPeriods=2 -cfl=10000 -solverh=gmres -rtolh=1.0e-9 -atolh=1.0e-9 -maxith=50000 -solveri=multigrid -rtoli=1e-8 -atoli=1e-8 -go=k > square1024MultiGridKrylov.out

square1024dhs:
	setsid $(cgwh) -noplot gaussian.cmd -g=square1024.order2 -x0=0.4 -y0=0.4 -nf=1 -freq=11 -amp=100 -beta=20 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-8 -upwind=0 -imode=1 -maxIterations=100 -numPeriods=2 -cfl=10000 -solverh=gmres -rtolh=1.0e-9 -atolh=1.0e-9 -maxith=50000 -solveri=multigrid -rtoli=1e-8 -atoli=1e-8 -go=d > square1024DHS.out

square2048mg:
	setsid $(cgwh) -noplot gaussian.cmd -g=square2048.order2 -x0=0.4 -y0=0.4 -nf=1 -freq=11 -amp=100 -beta=20 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-8 -upwind=0 -imode=1 -maxIterations=100 -numPeriods=2 -cfl=10000 -solverh=gmres -rtolh=1.0e-9 -atolh=1.0e-9 -maxith=50000 -solveri=multigrid -rtoli=1e-8 -atoli=1e-8 -go=k > square2048MultiGridKrylov.out

square2048dhs:
	setsid $(cgwh) -noplot gaussian.cmd -g=square2048.order2 -x0=0.4 -y0=0.4 -nf=1 -freq=11 -amp=100 -beta=20 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-8 -upwind=0 -imode=1 -maxIterations=100 -numPeriods=2 -cfl=10000 -solverh=gmres -rtolh=1.0e-9 -atolh=1.0e-9 -maxith=5000000 -solveri=multigrid -rtoli=1e-8 -atoli=1e-8 -go=d > square2048DHS.out


# =========== double-ellipse =======

# G2 O2
doubleEllipseG2O2:
	setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride2.order2.hdf -ts=implicit -cfl=1000 -x0=-1.4 -y0=-5.6 -nf=1 -amp=-400 -beta=10 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG2O2EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G2Order2Deflate64 -show=darkCornerRoomFreq10Np8G2Order2Deflate64.show -go=dfks > darkCornerRoomFreq10Np8G2Order2Deflate64.out
# G4 O2
doubleEllipseG4O2:
	setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride4.order2.hdf -ts=implicit -cfl=1000 -x0=-1.4 -y0=-5.6 -nf=1 -amp=-400 -beta=10 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG4O2EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G4Order2Deflate64 -show=darkCornerRoomFreq10Np8G4Order2Deflate64.show -go=dfks > darkCornerRoomFreq10Np8G4Order2Deflate64.out

# G8 O2
doubleEllipseG8O2:
	setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride8.order2.hdf -ts=implicit -cfl=1000 -x0=-1.4 -y0=-5.6 -nf=1 -amp=-400 -beta=10 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG8O2EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G8Order2Deflate64 -show=darkCornerRoomFreq10Np8G8Order2Deflate64.show -go=dfks > darkCornerRoomFreq10Np8G8Order2Deflate64.out

# G2 O4
doubleEllipseG2O4:
	setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride2.order4.hdf -ts=implicit -cfl=1000 -x0=-1.4 -y0=-5.6 -nf=1 -amp=-400 -beta=10 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG2O4EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G2Order4Deflate64 -show=darkCornerRoomFreq10Np8G2Order4Deflate64.show -go=dfks > darkCornerRoomFreq10Np8G2Order4Deflate64.out
# G4 O4
doubleEllipseG4O4:
	setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride4.order4.hdf -ts=implicit -cfl=1000 -x0=-1.4 -y0=-5.6 -nf=1 -amp=-400 -beta=10 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG4O4EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G4Order4Deflate64 -show=darkCornerRoomFreq10Np8G4Order4Deflate64.show -go=dfks > darkCornerRoomFreq10Np8G4Order4Deflate64.out

#G8 O4
doubleEllipseG8O4:
	setsid $(cgwh) -noplot gaussian.cmd -g=darkCornerRoomGride8.order4.hdf -ts=implicit -cfl=1000 -x0=-1.4 -y0=-5.6 -nf=1 -amp=-400 -beta=10 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=8 -freq=10 -bc=d -implicitUpwind=0 -minStepsPerPeriod=10 -eigenVectorFile=/home/henshw/runs/eig/darkCornerRoomG8O4EigsEv1024.show -deflateWaveHoltz=1 -numToDeflate=64 -eigTol=1.0e-5 -matlab=darkCornerRoomFreq10Np8G8Order4Deflate64 -show=darkCornerRoomFreq10Np8G8Order4Deflate64.show -go=dfks > darkCornerRoomFreq10Np8G8Order4Deflate64.out

# ================= PIPE ==========

pipe: pipeG2O4e pipeG2O4i

pipeG2O4i:
	setsid $(cgwh) -noplot gaussian.cmd -g=pipeze2.order4p.ng3 -x0=0.1 -y0=0.1 -z0=0.5 -nf=1 -freq=8.5 -amp=150 -beta=10 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=4 -cfl=10000 -solveri=cgs -rtoli=1e-10 -atoli=1e-10 -solverh=yale -restart=100 -iluh=50 -maxith=10000 -rtolh=1.0e-12 -atolh=1.0e-12 -matlab=pipeFreq8p5Np4G2O4Imp -go=fk > pipeFreq8p5Np4G2O4Imp.out

pipeG2O4e:
	setsid $(cgwh) -noplot gaussian.cmd -g=pipeze2.order4p.ng3 -x0=0.1 -y0=0.1 -z0=0.5 -nf=1 -freq=8.5 -amp=150 -beta=10 -adjustOmega=1 -ts=explicit -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -imode=1 -maxIterations=50 -numPeriods=4 -cfl=.95 -solveri=cgs -rtoli=1e-10 -atoli=1e-10 -solverh=yale -restart=100 -iluh=50 -maxith=10000 -rtolh=1.0e-12 -atolh=1.0e-12 -matlab=pipeFreq8p5Np4G2O4Exp -go=fk > pipeFreq8p5Np4G2O4Exp.out

# ================= disk ==================
# EXPLICIT NF=1, NP=1, O4, G8
diskOneFreqNp1Exp:
	$(cgwh) -noplot gaussian.cmd -g=sice8.order4.ng3 -x0=.25  -y0=.25 -nf=1 -freq=8.1 -amp=150 -beta=150 -adjustOmega=1 -ts=explicit -orderInTime=2 -solver=none -tol=1e-12 -upwind=1 -imode=1 -maxIterations=100 -numPeriods=1 -cfl=0.95 -matlab=diskOneFreqNp1G8Order4Exp -go=dfk > diskOneFreqNp1G8Order4Exp.out

# EXPLICIT NF=1, NP=4
diskOneFreqNp4Exp:
	$(cgwh) -noplot gaussian.cmd -g=sice4.order2 -x0=.25  -y0=.25 -nf=1 -freq=8.1 -amp=150 -beta=150 -adjustOmega=1 -ts=explicit -orderInTime=2 -solver=none -tol=1e-12 -upwind=1 -imode=1 -maxIterations=100 -numPeriods=4 -cfl=0.95 -matlab=diskOneFreqNp4G4Order2Exp -go=dfk > diskOneFreqNp4G8Order4Exp.out

diskOneFreqNp1Imp:
	$(cgwh) -noplot gaussian.cmd -g=sice4.order4.ng3 -x0=.25  -y0=.25 -nf=1 -freq=8.1 -amp=150 -beta=150 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-12 -upwind=1 -imode=1 -maxIterations=100 -numPeriods=1 -cfl=10000 -matlab=diskOneFreqNp1G4Order4Imp -go=dfk > diskOneFreqNp1G4Order4Imp.out

diskOneFreqNp4Imp:
	$(cgwh) -noplot gaussian.cmd -g=sice4.order4.ng3 -x0=.25  -y0=.25 -nf=1 -freq=8.1 -amp=150 -beta=150 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-12 -upwind=1 -imode=1 -maxIterations=100 -numPeriods=4 -cfl=10000 -matlab=diskOneFreqNp4G4Order4Imp -go=dfk > diskOneFreqNp4G4Order4Imp.out

# BOX IMP:
# cgwh -noplot gaussian.cmd -g=box2.order4.ng3 -x0=0.4 0.3 -0.2 0.25 -y0=0.6 0.25 -0.3 0.15 -nf=1 -freq=8.1 3 9 18 22 -amp=150 150 150 150 -beta=150 150 150 150 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-12 -upwind=1 -imode=1 -maxIterations=100 -numPeriods=4 -cfl=10000 -matlab=boxOneFreqNp4G2Order4Imp -go=dfk > boxOneFreqNp4G2Order4Imp.out

# BOX EXP 
# cgwh gaussian.cmd -g=box2.order4.ng3 -x0=0.4 0.3 -0.2 0.25 -y0=0.6 0.25 -0.3 0.15 -nf=1 -freq=8.1 3 9 18 22 -amp=150 150 150 150 -beta=150 150 150 150 -adjustOmega=1 -ts=explicit -orderInTime=2 -solver=none -tol=1e-12 -upwind=1 -imode=1 -maxIterations=100 -numPeriods=4 -cfl=.95 -matlab=boxOneFreqNp4G2Order4Exp -go=dfk

# SPHERE2 Order 4 NP=4 - IMP  INCLUDE DIRECT SOLVE: 
# setsid $cgwh -noplot gaussian.cmd -g=spheree2.order4.ng3 -x0=0.1 0.3 -0.2 0.25 -y0=0.1 0.25 -0.3 0.15 -nf=1 -freq=8.5 3 9 18 22 -amp=150 150 150 150 -beta=150 150 150 150 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-12 -upwind=1 -imode=1 -maxIterations=50 -numPeriods=4 -cfl=10000 -solveri=cgs -rtoli=1e-8 -atoli=1e-8 -solverh=gmres -restart=100 -iluh=100 -maxith=10000 -rtolh=1.0e-10 -atolh=1.0e-10 -matlab=sphereOneFreqNp4G2Order4Imp -go=dfk > sphereOneFreqNp4G2Order4Imp.out

# SHAPES HIGH FREQ: 
# cgwh -noplot gaussian.cmd -g=shapese4.order4.ng3 -x0=0.25 -0.55 0.5 -y0=-0.6 0.5 0.5 -nf=1 -freq=45 -amp=100 50 50 -beta=400 400 400 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=1 -imode=1 -maxIterations=50 -numPeriods=80 -cfl=10000 -matlab=shapesG4O4freq45Np80 -show=shapesG4O4freq45Np80.show -go=dks > shapesG4O4freq45Np80.out

# cgwh -noplot gaussian.cmd -g=shapese4.order4.ng3 -x0=0.25 -0.55 0.5 -y0=-0.6 0.5 0.5 -nf=1 -freq=44 -amp=100 50 50 -beta=400 400 400 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=1 -imode=1 -maxIterations=50 -numPeriods=80 -cfl=10000 -matlab=shapesG4O4freq44Np80 -show=shapesG4O4freq44Np80.show -go=dks > shapesG4O4freq44Np80.out

# vmax=.2
# cgwh -noplot gaussian.cmd -g=shapese4.order4.ng3 -x0=0.25 -0.55 0.5 -y0=-0.6 0.5 0.5 -nf=1 -freq=43 -amp=100 50 50 -beta=400 400 400 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-10 -upwind=1 -imode=1 -maxIterations=50 -numPeriods=80 -cfl=10000 -matlab=shapesG4O4freq43Np80 -show=shapesG4O4freq43Np80.show -go=dks > shapesG4O4freq43Np80.out

COPY_DIR = $(WAHO)/codes
COPY_FILES = \
   diskOneFreqNp1G4Order4ImpFP.m diskOneFreqNp1G4Order4ImpKrylov.m \
   diskOneFreqNp4G4Order4ImpFP.m diskOneFreqNp4G4Order4ImpKrylov.m

copyFiles: $(COPY_FILES)
	for u in $(COPY_FILES); do echo $$u; cp -f $$u $(COPY_DIR); done	