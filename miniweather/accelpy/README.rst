
================================
Python+Accelerator miniWeather
================================


With AccelPy, Python programmers can use multiple compiler-based programming models including GPU programming in their Python script.

At the time of this writing, Fortran code is used in this version.

Python miniWeather - AccelPy - Fortran
----------------------------------------

To support the execution of this version, "Makefile" is created.

As explained in the previous section, you may need to modify command arguments to accommodate your environment. Especially, you may need to get an interactive node first if you are on the Summit system of ORNL.


The version offloaded “compute_tendencies_z” to CPU using Fortran-backend.

You can find the following code in “compute_tendencies_z” function  of “miniweather_accelpy.py”::

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


“ord_tend_z” order file
------------------------

To specify the content of the offloading, the user needs to create an order file that contains the offload code.

The following is the content of the order file::

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

