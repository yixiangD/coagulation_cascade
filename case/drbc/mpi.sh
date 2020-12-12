#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=13:00:00

# Use 2 nodes with 8 tasks each, for 16 MPI tasks:
#SBATCH --nodes=2
#SBATCH --tasks-per-node=24

# Specify a job name:
#SBATCH -J MyMPIJob

# Specify an output file
#SBATCH -o MyMPIJob-%j.out


# Run a command


#mpirun -n 12 ./FCM_DIPOLE step.conf
#srun --mpi=pmi2 ~/lmp_src_group/src_cascade/lmp_stam -in chan.adr.run
srun --mpi=pmi2 ~/lmp_src_group/src_cascade/lmp_stam -in chan.rheo.run
#srun --mpi=pmi2 ~/lmp_src_group/src_cascade/lmp_stam -in chan.equil.run
