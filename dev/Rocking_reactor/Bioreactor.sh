#!/bin/bash

for MIN_LEVEL in 7; do
    for MAX_LEVEL in 9; do
	# rm -rf Bioreactor-MIN$MIN_LEVEL-MAX$MAX_LEVEL
	# mkdir Bioreactor-MIN$MIN_LEVEL-MAX$MAX_LEVEL
	rm -rf New14
	mkdir New14

	# qcc -O2 -w -fopenmp -Wall Bioreactor.c -lm -o Bioreactor-MIN$MIN_LEVEL-MAX$MAX_LEVEL/Bioreactor
	# export OMP_NUM_THREADS=1
	CC99='mpicc -std=c99' qcc -I /gpfs/scratch/mkim79/basilisk -Wall -O2 -D_MPI=1 Bioreactor.c -o New14/New14 -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm 
	# CC99='mpicc -std=c99' qcc -I /gpfs/scratch/mkim79/basilisk -Wall -O2 -D_MPI=1 Bioreactor.c -o New13/New13 -lm 


	# cd Bioreactor-MIN$MIN_LEVEL-MAX$MAX_LEVEL
	
	#mkdir Slices
	#mkdir Interfaces
	#mkdir Animations
	#mkdir Data_all

	# Bioreactor
	# 1. Width
	# 2. Rocking angle (degree)
	# 3. Rocking rpm
	# 4. MINLEVEL
	# 5. MAXLEVEL	
	
	#./Bioreactor 0.25 7 40 $MIN_LEVEL $MAX_LEVEL
	sbatch Bioreactor.slurm
	cd ..	
    done
done
