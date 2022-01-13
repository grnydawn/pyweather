#!/bin/bash


if [[ "$SLURM_PROCID" == "0" ]]
then
	python -m cProfile -s tottime ./miniweather_mpi.py > miniweather.prof 

else
	python ./miniweather_mpi.py
fi
