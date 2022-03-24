# NOTE: to install mpi4py manually, in login node, run "env MPICC=cc python -m pip install mpi4py --force --no-binary :all:"

PYTHON = $(shell which python)
PROJECT ?= ecp
WALLTIME ?= "00:01:00"
#NPROCS ?= 64
NPROCS ?= 32
#NPROCS ?= 1

spock: clean
	srun -n${NPROCS} python ./miniweather_accelpy.py

crusher: clean
	#module load craype-accel-amd-gfx90a rocm; srun -N1 -n${NPROCS} -Acli133 -t${WALLTIME} python ./miniweather_accelpy.py
	#module load craype-accel-amd-gfx90a rocm; srun -n${NPROCS} python ./miniweather_accelpy.py -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork
	module load craype-accel-amd-gfx90a rocm; srun -n${NPROCS} python ./miniweather_accelpy.py

summit: clean
	#jsrun -n ${NPROCS} ${PYTHON} ./miniweather_accelpy.py -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork
	jsrun -n ${NPROCS} ${PYTHON} ./miniweather_accelpy.py -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork

debug: clean
	python ./miniweather_accelpy.py --debug

crusher_debug: clean
	module load craype-accel-amd-gfx90a rocm; srun -n1 python ./miniweather_accelpy.py --debug

summit_debug: clean
	jsrun -n 1 ${PYTHON} ./miniweather_accelpy.py -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork

profile: clean
	srun -n ${NPROCS} -A cli133 -t ${WALLTIME} -p ${PROJECT} -o miniweather_accelpy.log  ./profile.sh


crusher_profile: clean
	module load craype-accel-amd-gfx90a rocm; srun -n${NPROCS} profile.sh -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork

valgrind_spock: clean
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes python ./miniweather_accelpy.py --debug &> valgrind.out


clean:
	rm -rf *.slab *.log *.prof *.o *.s core _accelpy_debug_ miniweather_accel._tmpwork_
