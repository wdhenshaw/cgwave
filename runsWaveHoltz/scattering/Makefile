cgwh= $(HOME)/research/cgwave/bin/cgwh
WAHO= $(HOME)/Dropbox/AMP/waho


# ============= SCATTERING FROM A CYLINDER   =============


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

disk: diskG8 diskG16 diskG32 diskG64

diskG8:
	$(cgwh) -noplot scat.cmd -g=cice8.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=4 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG8kx4 -show=diskG8kx4.show -go=dks > diskG8kx4.out
	cp diskG8kx4Krylov.m $(WAHO)/codes

diskG16:
	setsid $(cgwh) -noplot scat.cmd -g=cice16.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=8 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG16kx8 -show=diskG16kx8.show -go=dks > diskG16kx8.out
	cp diskG16kx8Krylov.m $(WAHO)/codes

diskG32:
	setsid $(cgwh) -noplot scat.cmd -g=cice32.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=16 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG32kx16 -show=diskG32kx16.show -go=ks > diskG32kx16.out
	cp diskG32kx16Krylov.m $(WAHO)/codes

diskG64:
	setsid $(cgwh) -noplot scat.cmd -g=cice64.order2 -ts=implicit -cfl=1000 -nf=1 -amp=1 -kx=32 -ky=0 -kz=0 -adjustOmega=1 -orderInTime=2 -solver=none -tol=1e-10 -upwind=0 -maxIterations=1000 -numPeriods=4 -filterTimeDerivative=1 -damp=0 -adjustHelmholtzForUpwinding=0 -useSuperGrid=0 -superGridWidth=0.15 -bc=d -bc1=a -bc2=a -bc3=a -bc4=a -implicitUpwind=0 -adjustPlotsForSuperGrid=0 -adjustErrorsForSuperGrid=0 -debug=0 -matlab=diskG64kx32 -show=diskG64kx32.show -go=ks > diskG64kx32.out
	cp diskG64kx32Krylov.m $(WAHO)/codes

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