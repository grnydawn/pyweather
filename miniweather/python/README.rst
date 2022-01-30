================================
Pure Python miniWeather
================================

To drive the execution of the version, "Makefile" is created. In the file, the Makefile targets are developed for Summit and Spock systems at ORNL.

To run the Makefile targets for Summit, you need to get an interactive node.

You may need to modify some arguments of shell command in the targets to accommodate your own environment. For example, the following is the command to run the version on Summit::

        >>> jsrun -n ${NRANKS} ${PYTHON} ./miniweather_mpi.py \
                  -o ${MEMBERWORK}/cli115/miniweather_mpi.slab \
                  -w ${MEMBERWORK}/cli115/pyweatherwork

“-o” is an optional argument of miniweather_mpi.py to set the file path of the output data file. Because, on Summit, a program that runs on a computing node is not allowed to write a file in the user's home directory, a path to a scratch file is specified with “-o” option.

Similarly, “-w” option specifies a path to a working directory in a scratch file space. Working directory is required to create temporary files for a parallel output data file.

Modify the paths to match with your environment.

Following is a screen output from running the command::

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

