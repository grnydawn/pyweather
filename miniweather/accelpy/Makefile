PYTHON = $(shell which python)
PROJECT ?= ecp
WALLTIME ?= "00:01:00"
NPROCS ?= 3

spock: clean
	srun -n ${NPROCS} -A cli133 -t ${WALLTIME} -p ${PROJECT} -o miniweather_accelpy.log  python ./miniweather_accelpy.py

summit: clean
	#jsrun -n ${NPROCS} ${PYTHON} ./miniweather_accelpy.py -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork
	jsrun -n ${NPROCS} ${PYTHON} ./miniweather_accelpy.py -o ${MEMBERWORK}/cli115/miniweather_accel.slab

debug: clean
	python ./miniweather_accelpy.py

summit_debug: clean
	jsrun -n 1 ${PYTHON} ./miniweather_accelpy.py -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork

profile: clean
	srun -n ${NPROCS} -A cli133 -t ${WALLTIME} -p ${PROJECT} -o miniweather_accelpy.log  ./profile.sh


valgrind_spock: clean
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes python ./miniweather_accelpy.py --debug &> valgrind.out


clean:
	rm -rf _accel_debug_ *.slab *.log *.prof
