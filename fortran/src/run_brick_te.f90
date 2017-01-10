!=================================================================================
! BRICK-TE: ! Simple model to simulate contribution of thermosteric expansion (TE)
! to global sea-level rise. (Grinsted, Moore and Jevrejeva 2010, equations 1/2
! that were originally applied to global sea-level)
!=================================================================================

!---------------------------------------------------------------------------------
subroutine run_brick_te(ns, tstep, brick_te_a, brick_te_b, brick_te_tau, &
                               brick_te_TE_0, Gl_Temp, TE_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     Gl_Temp  	Global mean surface temperature [degC]
! |     ns        	Number of timesteps
! |
! |    Parameters:
! |     tstep     	time step
! |     a         	sensitivity of equilibrium TE [m/degC]
! |     b         	equilibrium TE [m] for temperature Tg = 0
! |     tau       	timescale (efolding time)
! |
! |    Initial conditions:
! |     TE_0        initial thermosteric expansion
! |
! | Outputs:
! |     TE_out	Thermosteric expansion
!  =========================================================================

    USE global
    USE brick_te

    implicit none

    integer(i4b), intent(IN) :: ns ! time series length

! parameters
    real(DP),     intent(IN) :: tstep

    real(DP),     intent(IN) :: brick_te_a
    real(DP),     intent(IN) :: brick_te_b
    real(DP),     intent(IN) :: brick_te_tau

! intial conditions    
    real(DP),     intent(IN) :: brick_te_TE_0

! input variables
    real(DP), dimension(ns), intent(IN)  :: Gl_Temp

! output variables    
    real(DP), dimension(ns), intent(OUT) :: TE_out

    integer(i4b) :: i ! time step counter
    

! Initialize brick_te (parameters and initial variable values)
    i = 1
    call init_brick_te(tstep, brick_te_a, brick_te_b, brick_te_tau, &
                          brick_te_TE_0, TE_out(i) )
                                  
! estimate outputs
    do i=2,ns
    
        ! global sea level rise
        call brick_te_step(Gl_Temp(i-1), TE_out(i) )
                        
    end do

    RETURN

end subroutine run_brick_te