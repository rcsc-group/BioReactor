#!/bin/bash

for MIN_LEVEL in 7; do
    for MAX_LEVEL in 9; do
	rm -rf Bioreactor-MIN$MIN_LEVEL-MAX$MAX_LEVEL
	mkdir Bioreactor-MIN$MIN_LEVEL-MAX$MAX_LEVEL

	qcc -I -O2 -w -fopenmp -Wall Bioreactor.c -o Bioreactor-MIN$MIN_LEVEL-MAX$MAX_LEVEL/Bioreactor -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
	export OMP_NUM_THREADS=2

	cd Bioreactor-MIN$MIN_LEVEL-MAX$MAX_LEVEL

	mkdir Slices
	mkdir Interfaces
	mkdir Animations
	mkdir Data_all

	# Bioreactor
	# 1. Width
	# 2. Rocking angle (degree)
	# 3. Rocking rpm
	# 4. MINLEVEL
	# 5. MAXLEVEL	
	
	./Bioreactor 0.25 7 40 $MIN_LEVEL $MAX_LEVEL
    done
done
