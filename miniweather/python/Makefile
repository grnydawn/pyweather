
PYTHON = $(shell which python)
NRANKS = 1
WALLTIME = "00:00:30"

ACCOUNT_SPOCK ?= CLI133
PROJECT_SPOCK ?= ecp

PROJECT_SUMMIT ?= cli115

spock: clean
	srun -n ${NRANKS} -A ${ACCOUNT_SPOCK} -t ${WALLTIME} -p ${PROJECT_SPOCK} -o miniweather.log  python ./miniweather_mpi.py

summit: clean
	jsrun -n ${NRANKS} ${PYTHON} ./miniweather_mpi.py -w ${MEMBERWORK}/${LSB_PROJECT_NAME}/pyweatherwork

profile_spock: clean
	srun -n ${NRANKS} -A ${ACCOUNT_SPOCK} -t ${WALLTIME} -p ${PROJECT_SPOCK} -o miniweather.log  ./profile.sh

profile_summit: clean
	jsrun -n ${NRANKS} -p ${PROJECT_SUMMIT} -o miniweather.log  ./profile.sh

valgrind_spock: clean
	srun -n ${NRANKS} -A ${ACCOUNT_SPOCK} -t ${WALLTIME} -p ${PROJECT_SPOCK} -o miniweather.log  valgrind --leak-check=yes python ./miniweather_mpi.py

clean:
	rm -f *.slab *.log *.prof
