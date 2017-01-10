SUBROUTINE run_doeclim(ns, time_out, forcing_in, t2co_in, kappa_in, temp_out, &
                        heatflux_mixed_out, heatflux_interior_out)

    USE global
    USE doeclim

    implicit none

    integer(i4b), intent(IN) :: ns
    real(DP), intent(IN) :: t2co_in
    real(DP), intent(IN) :: kappa_in
    real(DP), dimension(ns), intent(IN) :: forcing_in
    real(DP), dimension(ns), intent(OUT) :: time_out
    real(DP), dimension(ns), intent(OUT) :: temp_out
    real(DP), dimension(ns), intent(OUT) :: heatflux_mixed_out
    real(DP), dimension(ns), intent(OUT) :: heatflux_interior_out

    integer(i4b) :: i
    integer(i4b) :: start_year = 1850

! Assign global variables.
    nsteps = ns
    deltat = 1.0d0

    call init_doeclim_arrays()

    call init_doeclim_parameters(t2co_in, kappa_in)

    do i = 1,nsteps
        call doeclimtimestep_simple(i, forcing_in(i), temp_out(i))

        time_out(i) = start_year + (i-1)*deltat
    end do

    heatflux_mixed_out = heatflux_mixed
    heatflux_interior_out = heatflux_interior

    call dealloc_doeclim()

    RETURN

END SUBROUTINE run_doeclim
