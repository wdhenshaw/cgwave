# Run cgwave and extract timing and memory usage for cgwave/doc/performance.tex
# Examples: 
# make -s square
# make -s nonSquare

cgWave = $(CGWAVE)/bin/cgWave

meApproach = ha
# meApproach = std

all: square nonSquare

# double back-slash 
DBS = '\\\\'

# square: 
square:
	echo "%     Grid      & ord  & scheme & pts   & total    &  adv    & adv-r-c  &   up  &  bc  & interp"  $(DBS)
	$(cgWave) -noplot ../known/known.cmd -g=square1024.order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=square1024.order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=square1024.order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	echo "\hline"
	
# square: measure cycles: NOTE -bc=e (missing order=8 CBC)
squareGrid = square1024
# squareGrid = square512
squareCycles:
squareCycles:
	echo "%     Grid      & ord  & scheme & pts   & total    &  adv    & adv-r-c  &   up  &  bc  & interp" $(DBS)
	$(cgWave) -noplot ../known/known.cmd -g=$(squareGrid).order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(squareGrid).order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(squareGrid).order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(squareGrid).order8.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	echo "\hline"

# nonSquare: 
nonSquare:
	echo "%     Grid      & ord  & scheme & pts   & total    &  adv    & adv-r-c  &   up  &  bc  & interp" $(DBS)
	$(cgWave) -noplot ../known/known.cmd -g=nonSquare1024.order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=nonSquare1024.order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=nonSquare1024.order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	echo "\hline"

# nonSquare: measure cycles
nonSquareGrid = nonSquare1024
# nonSquareGrid = nonSquare512
nonSquareCycles:
	echo "%     Grid      & ord  & scheme & pts   & total    &  adv    & adv-r-c  &   up  &  bc  & interp"  $(DBS)
	$(cgWave) -noplot ../known/known.cmd -g=$(nonSquareGrid).order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(nonSquareGrid).order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(nonSquareGrid).order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(nonSquareGrid).order8.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	echo "\hline"

# box: Use bc=e for now
box:  
	echo "%     Grid      & ord  & scheme & pts   & total    &  adv    & adv-r-c  &   up  &  bc  & interp"  $(DBS)
	$(cgWave) -noplot ../known/known.cmd -g=box16.order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=box16.order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=box16.order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p cgWave.out
	echo "\hline"

boxGrid = box16
# boxGrid = box8
boxCycles:  
	echo "%     Grid      & ord  & scheme & pts   & total    &  adv    & adv-r-c  &   up  &  bc  & interp"  $(DBS)
	$(cgWave) -noplot ../known/known.cmd -g=$(boxGrid).order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(boxGrid).order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(boxGrid).order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(boxGrid).order8.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	echo "\hline"


# nonBox: Use bc=e for now
nonBox:
	echo "%     Grid      & ord  & scheme & pts   & total    &  adv    & adv-r-c  &   up  &  bc  & interp"  $(DBS)
	$(cgWave) -noplot ../known/known.cmd -g=nonBox16.order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=nonBox16.order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=nonBox16.order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=nonBox16.order8.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p cgWave.out
	echo "\hline"	

# nonBox: Use bc=e for now
# nonBoxGrid = nonBox16
nonBoxGrid = nonBox8
nonBoxCycles:
	echo "%     Grid      & ord  & scheme & pts   & total    &  adv    & adv-r-c  &   up  &  bc  & interp"  $(DBS)
	$(cgWave) -noplot ../known/known.cmd -g=$(nonBoxGrid).order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(nonBoxGrid).order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(nonBoxGrid).order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=$(nonBoxGrid).order8.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=-1 -debug=0 -meApproach=$(meApproach) -bc=e -go=go > junk.out; getPerfResults.p -cycles cgWave.out
	echo "\hline"	
