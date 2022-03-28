
MODULE modompkernel

public runkernel_0

CONTAINS

INTEGER (C_INT64_T) FUNCTION runkernel_0(hs, nx, nz, NUM_VARS, hv_beta, dz, dt, sten_size, ID_DENS, ID_UMOM, ID_WMOM, ID_RHOT, c0, gamma, grav) BIND(C, name="runkernel_0")
    USE, INTRINSIC :: ISO_C_BINDING
    USE MODOMPDATA, ONLY : state
USE MODOMPDATA, ONLY : hy_dens_int
USE MODOMPDATA, ONLY : hy_pressure_int
USE MODOMPDATA, ONLY : hy_dens_theta_int
USE MODOMPDATA, ONLY : flux
USE MODOMPDATA, ONLY : tend

    INTEGER (C_INT64_T), INTENT(IN) :: hs
INTEGER (C_INT64_T), INTENT(IN) :: nx
INTEGER (C_INT64_T), INTENT(IN) :: nz
INTEGER (C_INT64_T), INTENT(IN) :: NUM_VARS
REAL (C_DOUBLE), INTENT(IN) :: hv_beta
REAL (C_DOUBLE), INTENT(IN) :: dz
REAL (C_DOUBLE), INTENT(IN) :: dt
INTEGER (C_INT64_T), INTENT(IN) :: sten_size
INTEGER (C_INT64_T), INTENT(IN) :: ID_DENS
INTEGER (C_INT64_T), INTENT(IN) :: ID_UMOM
INTEGER (C_INT64_T), INTENT(IN) :: ID_WMOM
INTEGER (C_INT64_T), INTENT(IN) :: ID_RHOT
REAL (C_DOUBLE), INTENT(IN) :: c0
REAL (C_DOUBLE), INTENT(IN) :: gamma
REAL (C_DOUBLE), INTENT(IN) :: grav

    
    integer , parameter :: rp = selected_real_kind(15)

    integer :: i,k,ll,s
    real(rp) :: r,u,w,t,p, stencil(4), d3_vals(NUM_VARS), vals(NUM_VARS), hv_coef

    !Compute the hyperviscosity coeficient
    hv_coef = -hv_beta * dz / (16*dt)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! TODO: THREAD ME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !Compute fluxes in the x-direction for each cell
!    !!$omp target teams distribute parallel do collapse(2) map(to: state) map(alloc: stencil,vals,d3_vals)
!    !!$omp target teams num_teams(nz+1) map(to: state) private(stencil,vals,d3_vals)
!    !!$omp target teams distribute parallel do simd collapse(2) private(stencil,vals,d3_vals)
!    !!$omp target teams distribute parallel do simd collapse(2) map(to: state)  map(alloc: stencil,vals,d3_vals)
    !$omp target teams num_teams(nz+1) map(to: state, hy_dens_int) map(alloc: stencil,vals,d3_vals)
    !$omp distribute
    do k = 1 , nz+1
      !$omp parallel do private(stencil,vals,d3_vals)
      do i = 1 , nx
        !Use fourth-order interpolation from four cell averages to compute the value at the interface in question
        do ll = 1 , NUM_VARS
          do s = 1 , sten_size
            stencil(s) = state(i,k-hs-1+s,ll)
          enddo
          !Fourth-order-accurate interpolation of the state
          vals(ll) = -stencil(1)/12 + 7*stencil(2)/12 + 7*stencil(3)/12 - stencil(4)/12
          !First-order-accurate interpolation of the third spatial derivative of the state
          d3_vals(ll) = -stencil(1) + 3*stencil(2) - 3*stencil(3) + stencil(4)
        enddo

        !Compute density, u-wind, w-wind, potential temperature, and pressure (r,u,w,t,p respectively)
        ! the values of ID_* are incremented by 1 because Fortran index starts from 1, not 0
        r = vals(ID_DENS) + hy_dens_int(k)
        u = vals(ID_UMOM) / r
        w = vals(ID_WMOM) / r
        t = ( vals(ID_RHOT) + hy_dens_theta_int(k) ) / r
        p = C0*(r*t)**gamma - hy_pressure_int(k)
        !Enforce vertical boundary condition and exact mass conservation
        if (k == 1 .or. k == nz+1) then
          w                = 0
          d3_vals(ID_DENS) = 0
        endif

        !Compute the flux vector with hyperviscosity
        flux(i,k,ID_DENS) = r*w     - hv_coef*d3_vals(ID_DENS)
        flux(i,k,ID_UMOM) = r*w*u   - hv_coef*d3_vals(ID_UMOM)
        flux(i,k,ID_WMOM) = r*w*w+p - hv_coef*d3_vals(ID_WMOM)
        flux(i,k,ID_RHOT) = r*w*t   - hv_coef*d3_vals(ID_RHOT)
      enddo
      !$omp end parallel do
    enddo
    !$omp end target teams
!
!    !$omp target update from(flux)
!
!    !Use the fluxes to compute tendencies for each cell
    !$omp target teams distribute parallel do collapse(3)
    do ll = 1 , NUM_VARS
      do k = 1 , nz
        do i = 1 , nx
          tend(i,k,ll) = -( flux(i,k+1,ll) - flux(i,k,ll) ) / dz
          if (ll == ID_WMOM) then
            tend(i,k,ID_WMOM) = tend(i,k,ID_WMOM) - state(i,k,ID_DENS)*grav
          endif
        enddo
      enddo
    enddo
!
    !$omp target update from(tend)

    runkernel_0 = 0

END FUNCTION

END MODULE
