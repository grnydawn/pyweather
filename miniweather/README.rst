================================
miniWeather in Python
================================


This folder contains miniWeather version that are ported in Python.

Prerequisites
===================

To run Python versions of miniWeather, please install "numpy", "mpi4py", "accelpy" and "pyslabs" packages as shown below before running the versions::

	$ pip install numpy mpi4py accelpy pyslabs

Notes::

	To install mpi4py, a working version of MPI should be available.
	Before installing mpi4py, please make sure that "mpicc --version" commmand
	generates a compiler version information.


Notes to the users of Summit system in OLCF::

	You may need to get an interactive node to install mpi4py.
	Please try following command to install mpi4py on a summit interactive node

	$ module load gcc
	$ pip install mpi4py


Pure Python version
===================

To drive the execution of the version, a Makefile is created (https://github.com/grnydawn/pyweather/blob/master/miniweather/python/Makefile). In the file, the Makefile targets are developed for Summit and Spock systems at ORNL.

To run the Makefile targets for Summit, you need to get an interactive node.

You may need to modify some arguments of shell command in the targets to accommodate your own environment. For example, the following is the command to run the version on Summit.

>>> jsrun -n ${NRANKS} ${PYTHON} ./miniweather_mpi.py -o ${MEMBERWORK}/cli115/miniweather_mpi.slab -w ${MEMBERWORK}/cli115/pyweatherwork

“-o” is an optional argument of miniweather_mpi.py to set the file path of the output data file. Because, on Summit, a program that runs on a computing node is not allowed to write a file in the user's home directory, a path to a scratch file is specified with “-o” option.

Similarly, “-w” option specifies a path to a working directory in a scratch file space. Working directory is required to create temporary files for a parallel output data file.

Modify the paths to match with your environment.

Following is a screen output from running the command.

nx_glob=100, nz_glob=50
dx=200.000000, dz=200.000000, dt=0.666667
Elapsed Time:      0.000 /     10.000
Elapsed Time:      0.667 /     10.000
Elapsed Time:      1.333 /     10.000
Elapsed Time:      2.000 /     10.000
Elapsed Time:      2.667 /     10.000
Elapsed Time:      3.333 /     10.000
Elapsed Time:      4.000 /     10.000
Elapsed Time:      4.667 /     10.000
Elapsed Time:      5.333 /     10.000
Elapsed Time:      6.000 /     10.000
Elapsed Time:      6.667 /     10.000
Elapsed Time:      7.333 /     10.000
Elapsed Time:      8.000 /     10.000
Elapsed Time:      8.667 /     10.000
Elapsed Time:      9.333 /     10.000
Elapsed Time:     10.000 /     10.000
CPU Time: 7.334666
d_mass: -0.000111
d_te: -0.000188

Nx_glob is the total grid size in x dimension.
Nz_glob is the total grid size in z dimension.
Dx and dz are the length of each grid in x and z dimensions each.
Dt is a time increment at each time step.

In this simulation, the total simulation time is 10 seconds.

CPU Time is the elapsed time for time stepping.
D_mass is the difference in doman mass between the initial value and the final value.
D_te is the difference in doman energy between the initial value and the final value.


Python+Accelerator version
=============================


With AccelPy, Python programmers can use multiple compiler-based programming models including GPU programming in their Python script.

At the time of this writing, Fortran code is used in this version.

Python miniWeather - AccelPy - Fortran
----------------------------------------

To support the execution of this version, a Makefile is created (https://github.com/grnydawn/pyweather/blob/master/miniweather/accelpy/Makefile)

As explained in the previous section, you may need to modify command arguments to accommodate your environment. Especially, you may need to get an interactive node first if you are on the Summit system of ORNL.


The version offloaded “compute_tendencies_z” to CPU using Fortran-backend.

You can find the following code in “compute_tendencies_z” function  of “miniweather_accelpy.py”.


# The following is a list of input variables.
       inputs = [
            self.hs, self.nx, self.nz, NUM_VARS, state, self.hv_beta,
            self.dz, self.dt, self.sten_size, self.ID_DENS, self.ID_UMOM,
            self.ID_WMOM, self.ID_RHOT, self.hy_dens_int, self.c0,
            self.gamma, self.hy_pressure_int, self.grav, self.hy_dens_theta_int
        ]

# The following is a list of output variables
        outputs = [
            self.flux, self.tend
        ]

# an order file is read. Further explanation in later section.
       with open(ord_tend_z) as fp:
            order = accelpy.Order(fp.read())

# Create Accel object with the following arguments
#  inputs : a list of inputs (either numpy arrays or scalar variables)
#  order : an Order object created right before
#  outputs : a list of inputs (either numpy arrays or scalar variables)
#  kind : specify the kind of backend (fortran in this case)
#  debug : True or False, Turn on/off debugging mode

        accel = accelpy.Accel(*inputs, order, *outputs, kind=kind, debug=self.debug)

# specifies how to create threads. In this case, nteams=1, nworkers_per_team=1
        accel.run(nteams, nworkers_per_team)

# finish offloading
        accel.stop()


“Ord_tend_z” order file

To specify the content of the offloading, the user needs to create an order file that contains the offload code.

The following is the content of the order file (https://github.com/grnydawn/pyweather/blob/master/miniweather/accelpy/tend_z.ord)

# list of input argument names

inputs = ("hs", "nx", "nz", "NUM_VARS", "state", "hv_beta", "dz", "dt",
          "sten_size", "ID_DENS", "ID_UMOM", "ID_WMOM", "ID_RHOT",
          "hy_dens_int", "c0", "gamma", "hy_pressure_int", "grav",
          "hy_dens_theta_int")

#   list of output argument names
outputs = ("flux", "tend")

# Python function to set the names
set_argnames(inputs, outputs)


# fortran section
[fortran]

! fortran code copied from miniWeather_mpi.F90

    integer , parameter :: rp = selected_real_kind(15)

    integer :: i,k,ll,s
    real(rp) :: r,u,w,t,p, stencil(4), d3_vals(NUM_VARS), vals(NUM_VARS), hv_coef
    !Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * dz / (16*dt)

   …

       if (ll == ID_WMOM+1) then
            tend(i,k,ID_WMOM+1) = tend(i,k,ID_WMOM+1) - state(i,k,ID_DENS+1)*grav
          endif
        enddo
      enddo
    enddo


PySlabs Parallel I/O
========================

To generate data in parallel, PySlabs parallel I/O Python package is used (https://github.com/grnydawn/pyslabs)

To plot the output in PySlabs data format, a Python plotting script is created (https://github.com/grnydawn/pyweather/blob/master/miniweather/slabplot.py) with a Makefile to run the script(https://github.com/grnydawn/pyweather/blob/master/miniweather/Makefile)

Once you install the pyslabs package, “slabs” command is also installed together. Try to run following “slabs” command:

>>> slabs info <slab data file>
