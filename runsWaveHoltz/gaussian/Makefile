cgwh= $(HOME)/research/cgwave/bin/cgwh
WAHO= $(HOME)/Dropbox/AMP/waho

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