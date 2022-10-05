#!/bin/bash

# Compile code
qcc -O2 -w -fopenmp -Wall bio_rec_new4.c -lm -o BioReactor -L$BASILISK/gl -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11

# Specify parallelisation features		
export OMP_NUM_THREADS=4
		
# Run executable
./BioReactor

