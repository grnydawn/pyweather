
ACCOUNT ?= CLI133
PROJECT ?= ecp
WALLTIME ?= "00:01:00"
NPROCS ?= 1

spock: clean
	srun -n ${NPROCS} -A ${ACCOUNT} -t ${WALLTIME} -p ${PROJECT} -o miniweather_accelpy.log  python ./miniweather_accelpy.py


debug: clean
	python ./miniweather_accelpy.py

profile: clean
	srun -n ${NPROCS} -A ${ACCOUNT} -t ${WALLTIME} -p ${PROJECT} -o miniweather_accelpy.log  ./profile.sh


clean:
	rm -f *.slab *.log *.prof
