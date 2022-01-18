================================
Miniweather portings in Python
================================


This folder contains miniWeather ported to Python.

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

