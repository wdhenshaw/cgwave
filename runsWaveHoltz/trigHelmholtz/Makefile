
cgwh= $(HOME)/research/cgwave/bin/cgwh
WAHO= $(HOME)/Dropbox/AMP/waho


# Confirmationf of convergence analysis:
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

