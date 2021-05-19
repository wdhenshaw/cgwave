#
#   Compute CONVergence rates for CgWave
#
# FINISH ME  -- THIS VERSION FROM CgMx
#
common = $(CG)/common

# The convergence rate script is called conv.p
all:  
	ln -sf $(common)/src/conv.p .

cgWave:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgWave


#- 
#- 
#- # directory to store interface tables:
#- interfaceTables = ./interfaceTables
#- 
#- #  -- GDM interface convergence runs
#- #   Order = 2, 4
#- #   modeGDM = 
#- #   Np = 
#- 
#- 
#- # ---- GDM DISK in a SQUARE ---
#- #   Orders = 2,4
#- #   Modes = -1 (default NR-mode) and 1 (R-mode)
#- gdmDisk: gdmDiskOrder2  gdmDiskOrder4
#- 
#- gdmDiskOrder2:
#- 	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.5 -tp=.1 -order=2 -name="innerOuterFixede" -dm=gdm -modeGDM=-1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -diss=0 -useSosupDissipation=1 -ng=4
#- 	mv dielectricCylGDMNp1Np2Order2maxNormConvTable.tex $(interfaceTables)
#- 	mv dielectricCylNFDTDOrder2GDMNp1Np2maxNorm.m  $(interfaceTables)
#- 	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.5 -tp=.1 -order=2 -name="innerOuterFixede" -dm=gdm -modeGDM=1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -diss=0 -useSosupDissipation=1 -ng=4
#- 	mv dielectricCylGDMNp1Np2Mode1Order2maxNormConvTable.tex $(interfaceTables)
#- 	mv dielectricCylNFDTDOrder2GDMNp1Np2GdmMode1maxNorm.m $(interfaceTables)
#- 
#- gdmDiskOrder4:
#- 	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.5 -tp=.1 -order=4 -name="innerOuterFixede" -dm=gdm -modeGDM=-1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -diss=0 -useSosupDissipation=1 -ng=4
#- 	mv dielectricCylGDMNp1Np2Order4maxNormConvTable.tex $(interfaceTables)
#- 	mv dielectricCylNFDTDOrder4GDMNp1Np2maxNorm.m  $(interfaceTables)
#- 	conv.p dielectric.conv -cyl=1 -kx=3 -tf=.5 -tp=.1 -order=4 -name="innerOuterFixede" -dm=gdm -modeGDM=1 -npv=1 2  -a01=1 1.5 2. -a11=.1 .07 .05 -b01=1 1.5 2.  -b11=.2 .1 .05 -a02=1.2 1.4 -a12=.1 .15 -b02=.8 1.2 -b12=.05 .09  -diss=0 -useSosupDissipation=1 -ng=4
#- 	mv dielectricCylGDMNp1Np2Mode1Order4maxNormConvTable.tex $(interfaceTables)
#- 	mv dielectricCylNFDTDOrder4GDMNp1Np2GdmMode1maxNorm.m $(interfaceTables)


