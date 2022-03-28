
module modompdata
USE, INTRINSIC :: ISO_C_BINDING

REAL (C_DOUBLE), DIMENSION(:, :, :), POINTER :: state
REAL (C_DOUBLE), DIMENSION(:, :, :), POINTER :: state_tmp
REAL (C_DOUBLE), DIMENSION(:), POINTER :: hy_dens_cell
REAL (C_DOUBLE), DIMENSION(:), POINTER :: hy_dens_theta_cell
REAL (C_DOUBLE), DIMENSION(:), POINTER :: hy_dens_int
REAL (C_DOUBLE), DIMENSION(:), POINTER :: hy_dens_theta_int
REAL (C_DOUBLE), DIMENSION(:), POINTER :: hy_pressure_int
REAL (C_DOUBLE), DIMENSION(:, :, :), POINTER :: flux
REAL (C_DOUBLE), DIMENSION(:, :, :), POINTER :: tend

public dataenter_0, dataexit_0
public state, state_tmp, hy_dens_cell, hy_dens_theta_cell, hy_dens_int, hy_dens_theta_int, hy_pressure_int, flux, tend

contains

INTEGER (C_INT64_T) FUNCTION dataenter_0(lcio00, lci00, lci01, lci02, lci03, lci04, lci05, lal00, lal01, lal02, lal03, lal04, lal05) BIND(C, name="dataenter_0")
    USE, INTRINSIC :: ISO_C_BINDING

    REAL (C_DOUBLE), DIMENSION(-1:102, -1:52, 4), INTENT(INOUT), TARGET :: lcio00
REAL (C_DOUBLE), DIMENSION(-1:102, -1:52, 4), INTENT(INOUT), TARGET :: lci00
REAL (C_DOUBLE), DIMENSION(-1:52), INTENT(INOUT), TARGET :: lci01
REAL (C_DOUBLE), DIMENSION(-1:52), INTENT(INOUT), TARGET :: lci02
REAL (C_DOUBLE), DIMENSION(51), INTENT(INOUT), TARGET :: lci03
REAL (C_DOUBLE), DIMENSION(51), INTENT(INOUT), TARGET :: lci04
REAL (C_DOUBLE), DIMENSION(51), INTENT(INOUT), TARGET :: lci05
REAL (C_DOUBLE), DIMENSION(101, 51, 4), INTENT(INOUT), TARGET :: lal00
REAL (C_DOUBLE), DIMENSION(100, 50, 4), INTENT(INOUT), TARGET :: lal01

    state => lcio00
state_tmp => lci00
hy_dens_cell => lci01
hy_dens_theta_cell => lci02
hy_dens_int => lci03
hy_dens_theta_int => lci04
hy_pressure_int => lci05
flux => lal00
tend => lal01

    !$omp target enter data map(to:state)
!$omp target enter data map(to:state_tmp, hy_dens_cell, hy_dens_theta_cell, hy_dens_int, hy_dens_theta_int, hy_pressure_int)
!$omp target enter data map(alloc:flux, tend)

    dataenter = 0

END FUNCTION

INTEGER (C_INT64_T) FUNCTION dataexit_0() BIND(C, name="dataexit_0")
    USE, INTRINSIC :: ISO_C_BINDING

    !$omp target exit data map(from:state)

    dataexit = 0

END FUNCTION

end module
