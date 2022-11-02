#!/bin/bash

for MIN_LEVEL in 9 10; do
    for MAX_LEVEL in 10 11 12; do
	rm -rf BIO-MIN$MIN_LEVEL-MAX$MAX_LEVEL
	mkdir BIO-MIN$MIN_LEVEL-MAX$MAX_LEVEL

	qcc -O2 -w -fopenmp -Wall bio_rec_new5.c -lm -o BIO-MIN$MIN_LEVEL-MAX$MAX_LEVEL/Bioreactor
	export OMP_NUM_THREADS=4

	cd BIO-MIN$MIN_LEVEL-MAX$MAX_LEVEL

	# Bioreactor
	# 1. Width
	# 2. Rocking angle (degree)
	# 3. Rocking rpm
	# 4. MINLEVEL
	# 5. MAXLEVEL	
	
	./Bioreactor 0.4 7 40 $MIN_LEVEL $MAX_LEVEL
	cd ..	
    done
done


#rm -rf test1/
#mkdir test1
# qcc -O2 -w -fopenmp bio_rot_sq_oxy.c -o bio_rot_sq_oxy/bio_rot_sq_oxy -lm

# CC99='mpicc -std=c99' qcc -I /gpfs/scratch/mkim79/basilisk -Wall -O2 -D_MPI=1 bio_rec_new5.c -o oxy282/oxy282 -lm 

# CC99='mpicc -std=c99' qcc -I /gpfs/scratch/mkim79/basilisk -Wall -O2 -D_MPI=1 bio_rec_new5.c -o oxy34/oxy34 -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

# CC99='mpicc -std=c99' qcc -I /gpfs/scratch/mkim79/basilisk -Wall -O2 -D_MPI=1 bio_rec_new5.c -o oxy36/oxy36 -lm

#qcc -O2 -w -fopenmp -Wall bio_rec_new5.c -lm -o test1/test1
#export OMP_NUM_THREADS=2

#cd test1
#./test1 10

# CC99='mpicc -std=c99' qcc -I /gpfs/scratch/mkim79/basilisk -Wall -O2 -D_MPI=1 bio_sq_new_AMR.c -o new19/new19 -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

# everything from -L$BASILISK/gl onward is just for visualisation (but needs a correct installation, see .c)
# qcc -O2 -w -fopenmp RTI_AMR_diff_post.c -o RTI_AMR_diff_post_out/RTI_AMR_diff_post -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

# sbatch bio_rec_new2.slurm
