

# CASES FOR EigenWave PAPER


# DISK ORDER 2 -- save show for plotting
# cgwh eig.cmd -g=sice4.order2 -nf=1 -adjustOmega=0 -ts=implicit -cfl=1000 -orderInTime=2 -solver=none -tol=1e-16 -upwind=0 -imode=1 -maxIterations=30 -matlab=diskG4O2Eig10 -eigenVectorFile=diskG4O2Ev128.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -numPeriods=1 -freq=10 -numEigs=42 -go=halt -eigenSolver=default -show=diskG4O2Eigs.show

# ---- SQUARE ---
cgwh eig.cmd -g=square128.order4.ng3 -nf=1 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-14 -upwind=1 -imode=1 -maxIterations=30 -numPeriods=1 -cfl=1000 -matlab=eigSquareO4ImpEig4 -eigenVectorFile=square128O4Ev64.show -deflateWaveHoltz=0 -numToDeflate=1 -computeEigs=1 -freq=4.4 -go=halt

cgwh eig.cmd -g=square128.order4.ng3 -nf=1 -adjustOmega=1 -ts=explicit -orderInTime=2 -solver=none -tol=1e-14 -upwind=1 -imode=1 -maxIterations=30 -numPeriods=1 -cfl=0.95 -matlab=eigSquareO4ExpEig4 -eigenVectorFile=square128O4Ev64.show -deflateWaveHoltz=0 -numToDeflate=1 -computeEigs=1 -freq=4.4 -go=halt


cgwh eig.cmd -g=square128.order2 -nf=1 -adjustOmega=1 -ts=implicit -orderInTime=2 -solver=none -tol=1e-14 -upwind=1 -imode=1 -maxIterations=30 -numPeriods=1 -cfl=1000 -matlab=eigSquareO2ImpEig4 -eigenVectorFile=square128Ev64.show -deflateWaveHoltz=0 -numToDeflate=1 -computeEigs=1 -freq=4.4 -go=halt

cgwh eig.cmd -g=square128.order2 -nf=1 -adjustOmega=1 -ts=explicit -orderInTime=2 -solver=none -tol=1e-14 -upwind=1 -imode=1 -maxIterations=30 -numPeriods=1 -cfl=0.95 -matlab=eigSquareO2ExpEig4 -eigenVectorFile=square128Ev64.show -deflateWaveHoltz=0 -numToDeflate=1 -computeEigs=1 -freq=4.4 -go=halt



# --- DISK --- 
cgwh eig.cmd -g=sice4.order4.ng3 -x0=0.25 -nf=1 -amp=100 -beta=400 -adjustOmega=1 -ts=implicit -cfl=10000 -orderInTime=2 -solver=none -tol=1e-12 -upwind=1 -imode=1 -maxIterations=200 -matlab=diskG4O4ImpEig11p8 -eigenVectorFile=diskG4O4Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=11.8 -numPeriods=8 -numRitz=20 -show=diskG4O4ImpEig11p8.show -go=halt

cgwh eig.cmd -g=sice4.order4.ng3 -nf=1 -adjustOmega=1 -ts=explicit -cfl=.95 -orderInTime=2 -solver=none -tol=1e-12 -upwind=1 -imode=1 -maxIterations=200 -matlab=diskG4O4ExpEig11p8 -eigenVectorFile=diskG4O4Ev64.show -deflateWaveHoltz=0 -numToDeflate=0 -computeEigs=1 -freq=11.8 -numPeriods=8 -numRitz=20 -go=halt