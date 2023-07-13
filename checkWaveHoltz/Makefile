#
#  Makefile to run regression tests for WaveHoltz
#
check:  
	@echo "\n ------ RUN WAVEHOLTZ CHECKS ------\n"
	./checkWaveHoltz.p -solver=cgwh -testsFile=testsWaveHoltz

check-quiet:  
	./checkWaveHoltz.p -verbose=0 -solver=cgwh -testsFile=testsWaveHoltz

replace:  
	@echo "\n ------ REPLACE WAVEHOLTZ CHECKS ------\n"
	./checkWaveHoltz.p -solver=cgwh -testsFile=testsWaveHoltz	-replace
	