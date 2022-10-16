# Run cgwave and extract timing and memory usage for cgwave/doc/performance.tex
# Examples: 
# make -s square
# make -s nonSquare

cgWave = $(CGWAVE)/bin/cgWave

all: square nonSquare

# square: 
square:
	$(cgWave) -noplot ../known/known.cmd -g=square1024.order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=2 -debug=0 -meApproach=ha -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=square1024.order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=4 -debug=0 -meApproach=ha -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=square1024.order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=6 -debug=0 -meApproach=ha -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	echo "\hline"
	

# nonSquare: 
nonSquare:
	$(cgWave) -noplot ../known/known.cmd -g=nonSquare1024.order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=2 -debug=0 -meApproach=ha -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=nonSquare1024.order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=4 -debug=0 -meApproach=ha -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=nonSquare1024.order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=6 -debug=0 -meApproach=ha -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	echo "\hline"

# box: 
box:
	$(cgWave) -noplot ../known/known.cmd -g=box16.order2.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=2 -debug=0 -meApproach=ha -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=box16.order4.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=4 -debug=0 -meApproach=ha -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	$(cgWave) -noplot ../known/known.cmd -g=box16.order6.hdf -kx=1 -ky=1 -kz=1 -tf=.1 -tp=.1 -known=squareEig -upwind=0 -orderInTime=6 -debug=0 -meApproach=ha -bcApproach=cbc -go=go > junk.out; getPerfResults.p cgWave.out
	echo "\hline"