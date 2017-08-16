#!/bin/bash -l
#SBATCH -J test_vasp 
#SBATCH --qos=premium
#SBATCH -p regular
#SBATCH -N 10
#SBATCH -t 5:30:00

srun -n 240 ising
