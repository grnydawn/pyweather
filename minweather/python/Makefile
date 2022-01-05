
ACCOUNT ?= CLI133
PROJECT ?= ecp

spock:
	srun -n 3 -A ${ACCOUNT} -t 00:00:30 -p ${PROJECT} -o miniweather.log  python ./miniweather_mpi.py



clean:
	*.log *.slab
