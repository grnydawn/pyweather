
ACCOUNT ?= CLI133
PROJECT ?= ecp
WALLTIME ?= "00:01:00"
NPROCS ?= 3

spock: clean
	srun -n ${NPROCS} -A ${ACCOUNT} -t ${WALLTIME} -p ${PROJECT} -o miniweather_accelpy.log  python ./miniweather_accelpy.py


debug: clean
	python ./miniweather_accelpy.py

profile: clean
	srun -n ${NPROCS} -A ${ACCOUNT} -t ${WALLTIME} -p ${PROJECT} -o miniweather_accelpy.log  ./profile.sh


valgrind_spock: clean
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes python ./miniweather_accelpy.py --debug &> valgrind.out


clean:
	rm -rf _accel_debug_ *.slab *.log *.prof
