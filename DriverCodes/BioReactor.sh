#!/bin/bash

NAME='BioReactor'

for MIN_LEVEL in 7; do
    for MAX_LEVEL in 9; do
	
	rm -rf Data_all Data_specific Fig_vor Fig_vol Fig_tr Fig_oxy
	mkdir Data_all Data_specific Fig_vor Fig_vol Fig_tr Fig_oxy

	qcc -I -O2 -w -fopenmp -Wall BioReactor.c -o $NAME -L$BASILISK/gl -lglutils -lfb_tiny -lm
	export OMP_NUM_THREADS=2

	# Bioreactor parameters
	# 1. Width
	# 2. Rocking angle (degree)
	# 3. Rocking rpm
	
	./$NAME 0.25 7 37.5
    done
done
