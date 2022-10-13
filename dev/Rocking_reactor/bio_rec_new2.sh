#!/bin/bash

rm -rf rec38/
mkdir rec38

qcc -O2 -w -fopenmp bio_rec_new5.c -o rec38/rec38 -lm

# CC99='mpicc -std=c99' qcc -I /gpfs/scratch/mkim79/basilisk -Wall -O2 -D_MPI=1 bio_rec_new5.c -o rec38/rec38 -lm 

# CC99='mpicc -std=c99' qcc -I /gpfs/scratch/mkim79/basilisk -Wall -O2 -D_MPI=1 bio_sq_new_AMR.c -o new19/new19 -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

# everything from -L$BASILISK/gl onward is just for visualisation (but needs a correct installation, see .c)
# qcc -O2 -w -fopenmp RTI_AMR_diff_post.c -o RTI_AMR_diff_post_out/RTI_AMR_diff_post -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm

# sbatch bio_rec_new2.slurm

export OMP_NUM_THREADS=4

cd rec38
./rec38
