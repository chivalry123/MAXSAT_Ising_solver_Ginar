#!/bin/bash -l
#SBATCH -J test_vasp 
#SBATCH -p debug
#SBATCH -N 20
#SBATCH -t 0:30:00

srun -n 480 ising
