================================
miniWeather in Python
================================

miniWeather is a mini app simulating weather-like flows for training in parallelizing accelerated HPC architectures developed by Dr. Matthew Norman (https://github.com/mrnorman/miniWeather)

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

In "python" folder, a Fortran miniweather MPI version is ported in Python. Please see the readme file in the folder for details.

Python+Accelerator version
=============================

In "accelpy" folder, a Fortran miniweather MPI version is ported in Python with various accelerator extensions. Please see the readme file in the folder for details.

With AccelPy, Python programmers can use multiple compiler-based programming models including GPU programming in their Python script.

At the time of this writing, Fortran code is used in this version.


AccelPy - Scalable Accelerator Interface in Python
======================================================

To extend Python programs with multiple compiler-based programming models, `AccelPy package <https://github.com/grnydawn/accelpy>`_ is used.


PySlabs Parallel I/O
========================

To generate data in parallel, `PySlabs parallel I/O Python package <https://github.com/grnydawn/pyslabs>`_ is used.

Once you install the pyslabs package, “slabs” command is also installed together. Try to run following “slabs” command once you created a PySlabs data file::

        >>> slabs info <slab data file>


Plotting the output
--------------------

"slabplot.py" demonstrate how to plot output data in PySlabs data file. See Makefile for running the plotting script under "plot" Makefile target.
