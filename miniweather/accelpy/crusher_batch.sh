#!/bin/bash
#SBATCH -A cli133
#SBATCH -J ysjob
#SBATCH -o %x-%j.out
#SBATCH -t 00:01:00
#SBATCH -p ecp
#SBATCH -N 1

srun -n3 --ntasks-per-node=3 python ./miniweather_accelpy.py
