#
#  Makefile to run regression tests for EigenWave
#
check:  
	@echo "\n ------ RUN EIGENWAVE CHECKS ------\n"
	./checkEigen.p -solver=cgwh -testsFile=testsEigenWave

check-quiet:  
	./checkEigen.p -verbose=0 -solver=cgwh -testsFile=testsEigenWave

replace:  
	@echo "\n ------ REPLACE EIGENWAVE CHECKS ------\n"
	./checkEigen.p -solver=cgwh -testsFile=testsEigenWave	-replace
	