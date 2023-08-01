#!/bin/bash

NAME='Drop103'

# We're using a uniform grid; MIN_LEVEL and MAX_LEVEL aren't used in the simulation.
for MIN_LEVEL in 9; do
    for MAX_LEVEL in 9; do
	rm -rf Data_all Data_specific Fig_vor Fig_vol Fig_tr Fig_oxy
	
	qcc -I -O2 -w -fopenmp -Wall Bioreactor.c -o $NAME -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
	export OMP_NUM_THREADS=2
	# CC99='mpicc -std=c99' qcc -I /gpfs/scratch/mkim79/basilisk -Wall -O2 -D_MPI=1 230327_bio.c -o $NAME -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

	mkdir Data_all Data_specific Fig_vor Fig_vol Fig_tr Fig_oxy

	# Bioreactor
	# 1. Width
	# 2. Rocking angle (degree)
	# 3. Rocking rpm
	# 4. MINLEVEL
	# 5. MAXLEVEL	
	
	./$NAME 0.25 7 32.5 $MIN_LEVEL $MAX_LEVEL
    done
done
