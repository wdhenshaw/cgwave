
cgwh= $(HOME)/research/cgwave/bin/cgwh
WAHO= $(HOME)/Dropbox/AMP/waho


#SQUARE POLLUTION ERRORS

diskPollutionO4:
	$(cgwh) -noplot trigHelmholtz.cmd -g=square8.order4   -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG8O4Freq30.show   -go=ds >squareG8O4Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=square16.order4  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG16O4Freq30.show  -go=ds >squareG16O4Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=square32.order4  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG32O4Freq30.show  -go=ds >squareG32O4Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=square64.order4  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG64O4Freq30.show  -go=ds >squareG64O4Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=square128.order4 -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG128O4Freq30.show -go=ds >squareG128O4Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=square256.order4 -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG256O4Freq30.show -go=ds >squareG256O4Freq30.out

diskPollutionO2:
	$(cgwh) -noplot trigHelmholtz.cmd -g=square8.order2   -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG8O2Freq30.show   -go=ds >squareG8O2Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=square16.order2  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG16O2Freq30.show  -go=ds >squareG16O2Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=square32.order2  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG32O2Freq30.show  -go=ds >squareG32O2Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=square64.order2  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG64O2req30.show   -go=ds >squareG64O2Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=square128.order2 -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG128O2Freq30.show -go=ds >squareG128O2Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=square256.order2 -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=squareG256O2Freq30.show -go=ds >squareG256O2Freq30.out



# DISK POLLUTION ERRORS:

diskPollutionO4:
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice2.order4  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG2O4Freq30.show  -go=ds >diskG2O4Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice4.order4  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG4O4Freq30.show  -go=ds >diskG4O4Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice8.order4  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG8O4Freq30.show  -go=ds >diskG8O4Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice16.order4 -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG16O4Freq30.show -go=ds >diskG16O4Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice32.order4 -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG32O4Freq30.show -go=ds >diskG32O4Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice64.order4 -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG64O4Freq30.show -go=ds >diskG64O4Freq30.out

diskPollutionO2:
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice4.order2   -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG4O2Freq30.show   -go=ds >diskG4O2Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice8.order2   -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG8O2Freq30.show   -go=ds >diskG8O2Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice16.order2  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG16O2Freq30.show  -go=ds >diskG16O2Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice32.order2  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG32O2Freq30.show  -go=ds >diskG32O2Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice64.order2  -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG64O2Freq30.show  -go=ds >diskG64O2Freq30.out
	$(cgwh) -noplot trigHelmholtz.cmd -g=sice128.order2 -nf=1 -freq=30 -kx=3.4 -ky=3.4 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=0 -show=diskG128O2Freq30.show -go=ds >diskG128O2Freq30.out


# Confirmation of convergence analysis:
# Np=4 
# cgwh -noplot trigHelmholtz.cmd -g=square64.order4 -numPeriods=4 -nf=3 -freq=5.1 12 15.5 -omega=5.1 -kx=1.25 -ky=1.25 -ts=explicit -dtMax=.1 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-12 -tp=1000 -maxIterations=100 -matlab=square64Freq3Np4 -go=dfk > square64Freq3Np4.out
# Np=2 
# cgwh -noplot trigHelmholtz.cmd -g=square64.order4 -numPeriods=2 -nf=3 -freq=5.1 12 15.5 -omega=5.1 -kx=1.25 -ky=1.25 -ts=explicit -dtMax=.1 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-12 -tp=1000 -maxIterations=100 -matlab=square64Freq3Np2 -go=dfk > square64Freq3Np2.out
# Np=1 
# cgwh -noplot trigHelmholtz.cmd -g=square64.order4 -numPeriods=1 -nf=3 -freq=5.1 12 15.5 -omega=5.1 -kx=1.25 -ky=1.25 -ts=explicit -dtMax=.1 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-12 -tp=1000 -maxIterations=100 -matlab=square64Freq3Np1 -go=dfk > square64Freq3Np1.out


#
# ================= Three frequencies case study ==================
threeFreqOrder2: threeFreqImplicitSquare64 threeFreqExplicitSquare64 threeFreqImplicitSquare128 threeFreqExplicitSquare128
threeFreqImplicitSquare64:
	$(cgwh) -noplot trigHelmholtz.cmd -g=square64.order2 -numPeriods=4 -nf=3 -freq=5.1 12 15.5 -kx=1.75 -ky=1.75 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=1  -matlab=square64o2ImpThreeFreqNp4 -go=dfk

threeFreqExplicitSquare64:
	$(cgwh) -noplot trigHelmholtz.cmd -g=square64.order2 -numPeriods=4 -nf=3 -freq=5.1 12 15.5 -kx=1.75 -ky=1.75 -ts=explicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=.95 -upwind=1  -matlab=square64o2ExpThreeFreqNp4 -go=dfk

threeFreqImplicitSquare128:
	$(cgwh) -noplot trigHelmholtz.cmd -g=square128.order2 -numPeriods=4 -nf=3 -freq=5.1 12 15.5 -kx=1.75 -ky=1.75 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=1  -matlab=square128o2ImpThreeFreqNp4 -go=dfk

threeFreqExplicitSquare128:
	$(cgwh) -noplot trigHelmholtz.cmd -g=square128.order2 -numPeriods=4 -nf=3 -freq=5.1 12 15.5 -kx=1.75 -ky=1.75 -ts=explicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=.95 -upwind=1  -matlab=square128o2ExpThreeFreqNp4 -go=dfk

# -- order=4 
threeFreqOrder4: threeFreqImplicitSquare64Order4 threeFreqExplicitSquare64Order4 threeFreqImplicitSquare128Order4 threeFreqExplicitSquare128Order4
threeFreqImplicitSquare64Order4:
	$(cgwh) -noplot trigHelmholtz.cmd -g=square64.order4.ng3 -numPeriods=4 -nf=3 -freq=5.1 12 15.5 -kx=1.75 -ky=1.75 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=1  -matlab=square64o4ImpThreeFreqNp4 -go=dfk

threeFreqExplicitSquare64Order4:
	$(cgwh) -noplot trigHelmholtz.cmd -g=square64.order4.ng3 -numPeriods=4 -nf=3 -freq=5.1 12 15.5 -kx=1.75 -ky=1.75 -ts=explicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=.95 -upwind=1  -matlab=square64o4ExpThreeFreqNp4 -go=dfk

threeFreqImplicitSquare128Order4:
	$(cgwh) -noplot trigHelmholtz.cmd -g=square128.order4.ng3 -numPeriods=4 -nf=3 -freq=5.1 12 15.5 -kx=1.75 -ky=1.75 -ts=implicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=1000 -upwind=1  -matlab=square128o4ImpThreeFreqNp4 -go=dfk

threeFreqExplicitSquare128Order4:
	$(cgwh) -noplot trigHelmholtz.cmd -g=square128.order4.ng3 -numPeriods=4 -nf=3 -freq=5.1 12 15.5 -kx=1.75 -ky=1.75 -ts=explicit -dtMax=100 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-14 -tp=1000 -maxIterations=70 -cfl=.95 -upwind=1  -matlab=square128o4ExpThreeFreqNp4 -go=dfk



# copy files to paper directory:
COPY_DIR = $(WAHO)/codes
COPY_FILES = \
   square64o2ImpThreeFreqNp4FP.m square64o2ImpThreeFreqNp4Krylov.m square64o2ExpThreeFreqNp4FP.m square64o2ExpThreeFreqNp4Krylov.m \
   square128o2ImpThreeFreqNp4FP.m square128o2ImpThreeFreqNp4Krylov.m square128o2ExpThreeFreqNp4FP.m square128o2ExpThreeFreqNp4Krylov.m \
   square64o4ImpThreeFreqNp4FP.m square64o4ImpThreeFreqNp4Krylov.m square64o4ExpThreeFreqNp4FP.m square64o4ExpThreeFreqNp4Krylov.m \
   square128o4ImpThreeFreqNp4FP.m square128o4ImpThreeFreqNp4Krylov.m square128o4ExpThreeFreqNp4FP.m square128o4ExpThreeFreqNp4Krylov.m 

copyFiles: $(COPY_FILES)
	for u in $(COPY_FILES); do echo $$u; cp -f $$u $(COPY_DIR); done	

