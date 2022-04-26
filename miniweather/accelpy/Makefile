PYTHON = $(shell which python)
PROJECT ?= ecp
WALLTIME ?= "00:01:00"
NPROCS ?= 4

LANG    ?= "fortran" # "c"
ACCEL   ?= "fortran" # "c" "openmp" "openacc" "omptarger"
VENDOR  ?= "*" # "pgi" "gnu" "cray" "amd" "intel" "ibm"
NX      ?= 100 # 400 # 2000 # number of local grid cells in the x-dimension
NZ      ?= 50  # 200 # 1000 # number of local grid cells in the z-dimension
TIME    ?= 10  # 5          # total simulation time in seconds
FREQ    ?= 1   # 5          # frequency to perform output in seconds
OUTFILE ?= "miniweather_mpi.slab" # output data file in pyslabs format

spock: clean
	srun -n${NPROCS} python ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -x ${NX} -z ${NZ} -s ${TIME} -f ${FREQ} -o ${OUTFILE}

crusher: clean
	#module load craype-accel-amd-gfx90a rocm; srun -N1 -n${NPROCS} -Acli133 -t${WALLTIME} python ./miniweather_accelpy.py
	#module load craype-accel-amd-gfx90a rocm; srun -n${NPROCS} python ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -x ${NX} -z ${NZ} -s ${TIME} -f ${FREQ} -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork
	module load craype-accel-amd-gfx90a rocm; srun -n${NPROCS} python ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -x ${NX} -z ${NZ} -s ${TIME} -f ${FREQ} -o ${OUTFILE}

summit: clean
	#jsrun -n ${NPROCS} -g 1 ${PYTHON} ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork
	jsrun -n ${NPROCS} -g 1 ${PYTHON} ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -x ${NX} -z ${NZ} -s ${TIME} -f ${FREQ} -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork

debug: clean
	python ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -x ${NX} -z ${NZ} -s ${TIME} -f ${FREQ} -o ${OUTFILE}

crusher_debug: clean
	#module load craype-accel-amd-gfx90a rocm; srun -n1 python ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -x ${NX} -z ${NZ} -s ${TIME} -f ${FREQ} -o ${OUTFILE} --debug
	module load craype-accel-amd-gfx90a rocm; srun -n1 -c64 python ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -x ${NX} -z ${NZ} -s ${TIME} -f ${FREQ} -o ${OUTFILE}

summit_debug: clean
	#jsrun -n 21 -c 2 -a 1 -g 1 ${PYTHON} ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -x ${NX} -z ${NZ} -s ${TIME} -f ${FREQ} -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork
	jsrun -n1 -c4 -g 1 ${PYTHON} ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -x ${NX} -z ${NZ} -s ${TIME} -f ${FREQ} -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork

profile: clean
	srun -n ${NPROCS} -A cli133 -t ${WALLTIME} -p ${PROJECT} -o miniweather_accelpy.log  ./profile.sh

summit_nvprof: clean
	jsrun -n1 -c1 -g 1  nvprof --profile-child-processes ${PYTHON} ./miniweather_accelpy.py -l ${LANG} -a ${ACCEL} -v ${VENDOR} -x ${NX} -z ${NZ} -s ${TIME} -f ${FREQ} -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork


crusher_profile: clean
	#module load craype-accel-amd-gfx90a rocm; srun -n${NPROCS} profile.sh -o ${MEMBERWORK}/cli115/miniweather_accel.slab -w ${MEMBERWORK}/cli115/pyweatherwork
	module load craype-accel-amd-gfx90a rocm; srun -n1 -c64 ./profile.sh

valgrind_spock: clean
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes python ./miniweather_accelpy.py --debug &> valgrind.out


clean:
	rm -rf _accel_debug_ *.slab *.log *.prof *.so *.mod
