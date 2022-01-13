#!/bin/bash


if [[ "$SLURM_PROCID" == "0" ]]
then
	python -m cProfile -s tottime ./miniweather_accelpy.py > miniweather_accelpy.prof 

else
	python ./miniweather_accelpy.py
fi
