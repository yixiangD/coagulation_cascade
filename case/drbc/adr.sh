#!/bin/bash
#SBATCH -p normal
# Request an hour of runtime:
#SBATCH --time=48:00:00

# Use 2 nodes with 8 tasks each, for 16 MPI tasks:
#SBATCH -N 10
#SBATCH -n 640

# Specify a job name:
#SBATCH -J drbc

# Specify an output file
#SBATCH -o myMPIJob-%j.out
#SBATCH -e myMPIJob-%j.out


# Run a command


srun --mpi=pmi2 ../src_cascade/lmp_stam -in chan.adr.run
