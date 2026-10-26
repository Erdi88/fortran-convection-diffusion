FC     ?= gfortran
FFLAGS ?= -O2 -Wall

conv_diff: src/conv_diff.f90
	$(FC) $(FFLAGS) -o $@ $<

run: conv_diff
	cd run && ../conv_diff

plot: run
	python plot.py

clean:
	rm -f conv_diff run/file02 run/file_step50

.PHONY: run plot clean
