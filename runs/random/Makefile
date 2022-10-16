
all: disk cic pipe

cgWave = /home/henshw/Dropbox/research/cgwave/bin/cgWave
tFinal = 100000

# ------------------- SPHERE ------
sphere:
	$(cgWave) -noplot random.cmd -g=spheree4.order2.hdf -tf=10000 -tp=100 -bcApproach=cbc -assignInterpNeighbours=interp -show=randomSphereOrder2G4.show -go=go > randomSphereOrder2G4.out
	$(cgWave) -noplot random.cmd -g=spheree4.order4.ng3.hdf -tf=10000 -tp=100 -bcApproach=cbc -assignInterpNeighbours=interp -show=randomSphereOrder4G4.show -go=go > randomSphereeOrder4G4.out

sphere2:
	setsid $(cgWave) -noplot random.cmd -g=spheree4.order2.hdf -tf=10000 -tp=100 -bcApproach=cbc -assignInterpNeighbours=interp -show=randomSphereOrder2G4.show -go=go > randomSphereOrder2G4.out &

sphere4:
	setsid $(cgWave) -noplot random.cmd -g=spheree4.order4.ng3.hdf -tf=10000 -tp=100 -bcApproach=cbc -assignInterpNeighbours=interp -show=randomSphereOrder4G4.show -go=go > randomSphereeOrder4G4.out &


# ------------------- PIPE ------
pipe:
	$(cgWave) -noplot random.cmd -g=pipeze4.order2.hdf -tf=10000 -tp=100 -bcApproach=cbc -assignInterpNeighbours=interp -show=randomPipeOrder2G4.show -go=go > randomPipeOrder2G4.out
	$(cgWave) -noplot random.cmd -g=pipeze4.order4.ng3.hdf -tf=10000 -tp=100 -bcApproach=cbc -assignInterpNeighbours=interp -show=randomPipeOrder4G4.show -go=go > randomPipeOrder4G4.out


# ------------------- DISK ------
disk:
	$(cgWave) -noplot random.cmd -g=sice4.order2 -tf=$(tFinal) -tp=100 -bcApproach=cbc -assignInterpNeighbours=interp -show=randomDiskOrder2G4.show -go=go > randomDiskOrder2G4.out
	$(cgWave) -noplot random.cmd -g=sice4.order4.ng3 -tf=$(tFinal) -tp=100 -bcApproach=cbc -assignInterpNeighbours=interp -show=randomDiskOrder4G4.show -go=go > randomDiskOrder4G4.out

# ------------------ CIRCLE IN CHANNEL ---
cic:
	$(cgWave) -noplot random.cmd -g=cice4.order2 -tf=$(tFinal) -tp=100 -bcApproach=cbc -assignInterpNeighbours=interp -show=randomCICOrder2G4.show -go=go > randomCICOrder2G4.out
	$(cgWave) -noplot random.cmd -g=cice4.order4.ng3 -tf=$(tFinal) -tp=100 -bcApproach=cbc -assignInterpNeighbours=interp -show=randomCICOrder4G4.show -go=go > randomCICOrder4G4.out
