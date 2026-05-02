cgwh= $(HOME)/research/cgwave/bin/cgwh
WAHO= $(HOME)/Dropbox/AMP/waho
complexWaveHoltz= $(HOME)/Dropbox/complexWaveHoltz

CWH = $(HOME)/Dropbox/complexWaveHoltz


# ============ COMPLEX WAVEHOLTZ PAPER ========

diskO4CWH:
	$(cgwh) -noplot scat.cmd -g=cice4.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=4 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -maxIterations=70 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=scatCylO4G4kx4.show -takeImplicitFirstStep=1 -matlab=cicG4O4kx4Brrrrd -go=dfk
	cp cicG4O4kx4Brrrrd*.m $(CWH)/results

diskO4FWCWH:
	$(cgwh) -noplot scat.cmd -g=cice4.order4 -ts=implicit -beta2=0.5 -cfl=1000 -nf=1 -amp=1 -kx=4 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -maxIterations=70 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=scatCylO4G4kx4.show -takeImplicitFirstStep=1 -matlab=cicG4O4kx4BrrrrdFW -go=dfk
	cp cicG4O4kx4BrrrrdFW*.m $(CWH)/results

# O2
diskO2CWH:
	$(cgwh) -noplot scat.cmd -g=cice4.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=4 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -maxIterations=70 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=scatCylO2G4kx4.show -takeImplicitFirstStep=1 -matlab=cicG4O2kx4Brrrrd -go=dfk
	cp cicG4O2kx4Brrrrd*.m $(CWH)/results


# O2 full weighting
diskO2FWCWH:
	$(cgwh) -noplot scat.cmd -g=cice4.order2 -ts=implicit -beta2=0.5 -cfl=1000 -nf=1 -amp=1 -kx=4 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-12 -upwind=0 -maxIterations=70 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=scatCylO2G4kx4.show -takeImplicitFirstStep=1 -matlab=cicG4O2kx4BrrrrdFW -go=dfk
	cp cicG4O2kx4BrrrrdFW*.m $(CWH)/results


# ============= SCATTERING FROM 3 DISKS  =============

threeHoles: threeHolesG4kx1 threeHolesG8kx2 threeHolesG8kx3 threeHolesG16kx7 threeHolesG32kx12


threeHolesG4kx1:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e4.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=1 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=70 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG8O4Scatkx4.show -takeImplicitFirstStep=1 -matlab=threeHolesG4O4Scatkx1 -go=fk
	cp threeHolesG4O4Scatkx1FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG4O4Scatkx1FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results


threeHolesG8kx2:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e8.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=2 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=70 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG8O4Scatkx4.show -takeImplicitFirstStep=1 -matlab=threeHolesG8O4Scatkx2 -go=fk
	cp threeHolesG8O4Scatkx2FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG8O4Scatkx2FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results

threeHolesG8kx3:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e8.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=3 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=700 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG8O4Scatkx4.show -takeImplicitFirstStep=1 -matlab=threeHolesG8O4Scatkx3 -go=fks
	cp threeHolesG8O4Scatkx3FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG8O4Scatkx3FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results


threeHolesG8kx4:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e8.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=4 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=700 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG8O4Scatkx4.show -takeImplicitFirstStep=1 -matlab=threeHolesG8O4Scatkx4 -go=fks
	cp threeHolesG8O4Scatkx4FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG8O4Scatkx4FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results

# G16 kx=5 
threeHolesG16kx5:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e16.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=5 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=7000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG16O4Scatkx5.show -takeImplicitFirstStep=1 -matlab=threeHolesG16O4Scatkx5 -go=fks
	cp threeHolesG16O4Scatkx5FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG16O4Scatkx5FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results


# G16 kx=6 
threeHolesG16kx6:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e16.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=6 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=7000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG16O4Scatkx6.show -takeImplicitFirstStep=1 -matlab=threeHolesG16O4Scatkx6 -go=fks
	cp threeHolesG16O4Scatkx6FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG16O4Scatkx6FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results


# G16 kx=7
threeHolesG16kx7:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e16.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=7 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=7000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG16O4Scatkx7.show -takeImplicitFirstStep=1 -matlab=threeHolesG16O4Scatkx7 -go=fks
	cp threeHolesG16O4Scatkx7FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG16O4Scatkx7FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results


# G32 kx=9
threeHolesG32kx9:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e32.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=9 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=7000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG32O4Scatkx9.show -takeImplicitFirstStep=1 -matlab=threeHolesG32O4Scatkx9 -go=fks
	cp threeHolesG32O4Scatkx9FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG32O4Scatkx9FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results


# G32 kx=10
threeHolesG32kx10:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e32.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=10 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=7000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG32O4Scatkx10.show -takeImplicitFirstStep=1 -matlab=threeHolesG32O4Scatkx10 -go=fks
	cp threeHolesG32O4Scatkx10FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG32O4Scatkx10FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results

# G32 kx=11
threeHolesG32kx11:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e32.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=11 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=7000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG32O4Scatkx11.show -takeImplicitFirstStep=1 -matlab=threeHolesG32O4Scatkx11 -go=fks
	cp threeHolesG32O4Scatkx11FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG32O4Scatkx11FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results


# G32 kx=12 
threeHolesG32kx12:
	$(cgwh) -noplot scat.cmd -g=threeHolesGridN2M1e32.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=12 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=7000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -debug=0 -eigenVectorFile=none -show=threeHolesG32O4Scatkx12.show -takeImplicitFirstStep=1 -matlab=threeHolesG32O4Scatkx12 -go=fks
	cp threeHolesG32O4Scatkx12FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp threeHolesG32O4Scatkx12FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results






# ============= SCATTERING FROM A CYLINDER + EM BCS  =============

# disk + EM BC's : CHECK FOR NUMBER OF ITERATIONS VERSUS OMEGA 
disk: diskG4 diskG8 diskG16 diskG32 diskG64

diskG4:
	$(cgwh) -noplot scat.cmd -g=cice4.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=2 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG4kx2 -show=diskG4kx2.show -go=dfks > diskG4kx2.out
	cp diskG4kx2FD22TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp diskG4kx2FD22TSINp4FP.m $(complexWaveHoltz)/matlab/results


diskG8:
	$(cgwh) -noplot scat.cmd -g=cice8.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=4 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG8kx4 -show=diskG8kx4.show -go=dfks > diskG8kx4.out
	cp diskG8kx4FD22TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp diskG8kx4FD22TSINp4FP.m $(complexWaveHoltz)/matlab/results

diskG16:
	setsid $(cgwh) -noplot scat.cmd -g=cice16.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=8 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG16kx8 -show=diskG16kx8.show -go=dfks > diskG16kx8.out
	cp diskG16kx8FD22TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp diskG16kx8FD22TSINp4FP.m $(complexWaveHoltz)/matlab/results

diskG32:
	setsid $(cgwh) -noplot scat.cmd -g=cice32.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=16 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG32kx16 -show=diskG32kx16.show -go=fks > diskG32kx16.out
	cp diskG32kx16FD22TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp diskG32kx16FD22TSINp4FP.m $(complexWaveHoltz)/matlab/results

diskG64:
	setsid $(cgwh) -noplot scat.cmd -g=cice64.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=32 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG64kx32 -show=diskG64kx32.show -go=fks > diskG64kx32.out
	cp diskG64kx32FD22TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp diskG64kx32FD22TSINp4FP.m $(complexWaveHoltz)/matlab/results

# -- Order=4
diskO4: diskG2O4 diskG4O4 diskG8O4 diskG16O4 diskG32O4 

diskG2O4:
	$(cgwh) -noplot scat.cmd -g=cice2.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=2 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG2kx2 -show=diskG2O4kx2.show -go=dfks > diskG2O4kx2.out
	cp diskG2kx2FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp diskG2kx2FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results


diskG4O4:
	$(cgwh) -noplot scat.cmd -g=cice4.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=4 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG4kx4 -show=diskG4O4kx4.show -go=dfks > diskG4O4kx4.out
	cp diskG4kx4FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp diskG4kx4FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results

diskG8O4:
	$(cgwh) -noplot scat.cmd -g=cice8.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=8 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG8kx8 -show=diskG8kx8.show -go=dfks > diskG8kx8.out
	cp diskG8kx8FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp diskG8kx8FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results


diskG16O4:
	setsid $(cgwh) -noplot scat.cmd -g=cice16.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=16 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG16kx16 -show=diskG16kx16.show -go=fks > diskG16kx16.out
	cp diskG16kx16FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp diskG16kx16FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results

diskG32O4:
	setsid $(cgwh) -noplot scat.cmd -g=cice32.order4 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=32 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG32kx32 -show=diskG32kx32.show -go=fks > diskG32kx32.out
	cp diskG32kx32FD24TSINp4gmres.m $(complexWaveHoltz)/matlab/results
	cp diskG32kx32FD24TSINp4FP.m $(complexWaveHoltz)/matlab/results


# ============= SCATTERING FROM A CYLINDER  + SUPERGRID  =============
# disk + Supergrid: 
disksg: disk4sg disk8sg disk16sg disk32sg

disk4sg:
	setsid $(cgwh) -noplot scat.cmd -g=cice4.order2  -ts=implicit -cfl=10000 -nf=1 -amp=1 -kx=2  -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=1  -maxIterations=200 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=1 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=1 -adjustErrorsForSuperGrid=1 -debug=0 -matlab=diskG4kx2sg -show=diskG4kx2sg.show -go=ks > diskG4kx2sg.out
	cp diskG4kx2sgKrylov.m $(WAHO)/codes

disk8sg:
	setsid $(cgwh) -noplot scat.cmd -g=cice8.order2  -ts=implicit -cfl=10000 -nf=1 -amp=1 -kx=4  -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=1  -maxIterations=200 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=1 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=1 -adjustErrorsForSuperGrid=1 -debug=0 -matlab=diskG8kx4sg  -show=diskG8kx4sg.show -go=ks > diskG8kx4sg.out
	cp diskG8kx4sgKrylov.m $(WAHO)/codes

disk16sg:
	setsid $(cgwh) -noplot scat.cmd -g=cice16.order2 -ts=implicit -cfl=10000 -nf=1 -amp=1 -kx=8  -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=1  -maxIterations=200 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=1 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=1 -adjustErrorsForSuperGrid=1 -debug=0 -matlab=diskG16kx8sg  -show=diskG16kx8sg.show -go=ks > diskG16kx8sg.out
	cp diskG16kx8sgKrylov.m $(WAHO)/codes

disk32sg:
	setsid $(cgwh) -noplot scat.cmd -g=cice32.order2 -ts=implicit -cfl=10000 -nf=1 -amp=1 -kx=16 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=1  -maxIterations=200 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=1 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=1 -adjustErrorsForSuperGrid=1 -debug=0 -matlab=diskG32kx16sg -show=diskG32kx16sg.show -go=ks > diskG32kx16sg.out
	cp diskG32kx16sgKrylov.m $(WAHO)/codes



# ============= SCATTERING FROM THREE SHAPES   =============

shapes: shapesG4 shapesG8 shapesG16 shapesG32

shapesG4:
	setsid $(cgwh) -noplot scat.cmd -g=shapesBiggere4.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=4 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0  -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=1 -adjustErrorsForSuperGrid=1 -matlab=shapesG4kx4 -show=shapesG4kx4.show -debug=0 -go=dks > shapesG4kx4.out 
	cp shapesG4kx4Krylov.m $(WAHO)/codes
shapesG8:
	setsid $(cgwh) -noplot scat.cmd -g=shapesBiggere8.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=8 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0  -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=1 -adjustErrorsForSuperGrid=1 -matlab=shapesG84kx -show=shapesG8kx8.show -debug=0 -go=dks > shapesG84kx.out 
	cp shapesG8kx8Krylov.m $(WAHO)/codes

shapesG16:
	setsid $(cgwh) -noplot scat.cmd -g=shapesBiggere16.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=16 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0  -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=1 -adjustErrorsForSuperGrid=1 -matlab=shapesG16kx16 -show=shapesG16kx16.show -debug=0 -go=dks > shapesG16kx16.out
	cp shapesG16kx16Krylov.m $(WAHO)/codes

shapesG32:
	setsid $(cgwh) -noplot scat.cmd -g=shapesBiggere32.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=32 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0  -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=1 -adjustErrorsForSuperGrid=1 -matlab=shapesG32kx32 -show=shapesG32kx32.show -debug=0 -go=ks > shapesG32kx32.out 
	cp shapesG32kx32Krylov.m $(WAHO)/codes

# ============= SCATTERING FROM MULTI OBJECTS   =============

multiObjects: multiObjectsG4 multiObjectsG8

multiObjectsG4:	
	setsid $(cgwh) -noplot scat.cmd -g=multiObjectsGride4.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=0 -ky=4 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-8 -upwind=0 -maxIterations=10000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=multiObjectsG4ky4 -show=multiObjectsG4ky4.show -go=dks > multiObjectsG4ky4.out

multiObjectsG8:
	setsid $(cgwh) -noplot scat.cmd -g=multiObjectsGride8.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=0 -ky=8 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-8 -upwind=0 -maxIterations=10000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=multiObjectsG8ky8 -show=multiObjectsG8ky8.show -go=dks > multiObjectsG8ky8.out