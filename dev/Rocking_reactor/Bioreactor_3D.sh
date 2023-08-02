#!/bin/bash

NAME='Bioreactor3D'

# We're using a uniform grid; MIN_LEVEL and MAX_LEVEL aren't used in the simulation.
for MIN_LEVEL in 7; do # 9 as default 
    for MAX_LEVEL in 7; do # 9 as default
	rm -rf Data_all Data_specific Fig_vor Fig_vol Fig_tr Fig_oxy Slices Interfaces
	
	qcc -I -O2 -w -fopenmp -Wall Bioreactor_3D.c -o $NAME -lm # -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
	export OMP_NUM_THREADS=8
	# CC99='mpicc -std=c99' qcc -I /gpfs/scratch/mkim79/basilisk -Wall -O2 -D_MPI=1 230327_bio.c -o $NAME -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

	mkdir Data_all Data_specific Fig_vor Fig_vol Fig_tr Fig_oxy Slices Interfaces

	# Bioreactor
	# 1. Width
	# 2. Rocking angle (degree)
	# 3. Rocking rpm
	# 4. MINLEVEL
	# 5. MAXLEVEL	
	
	./$NAME 0.25 7 20.0 $MIN_LEVEL $MAX_LEVEL
    done
done
