#!/bin/bash

#SBATCH -J oxy282         # Job name
#SBATCH -o oxy282_out-%j.o         # Name of output file
#SBATCH -e oxy282_err-%j.e          # Name of error file
#SBATCH -p batch           # Name of queue or partition
#SBATCH -N 1                 # Total # of nodes (24-48 cores per each node / see https://docs.ccv.brown.edu/oscar/system-overview)
#SBATCH -n 10                # Total # of mpi tasks (not parallelized: choose 1)
#SBATCH -t 10:00:00         # Run time (hh:mm:ss)
#SBATCH --mail-user=weizh348@gmail.com
#SBATCH --mail-type=all     # Send email at the beginning and end of your job

cd oxy282
mkdir Data_all

srun --mpi=pmix ./oxy282

#########################################################################
# Notes:
# module avail: see all programs available
# module load: load a module to make a program available
# sbatch <scriptname>: submit a job to the queue
# squeue -u <username>: check on your jobs in the queue
# scancel <jobid>: cancel a running job
# myjobinfo -j <jobid>: see memory usage of completed jobs


