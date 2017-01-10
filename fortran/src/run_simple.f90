!=================================================================================
! SIMPLE: simple model for Greenland ice-sheet volume [m sle] (Bakker et al 2014)
!=================================================================================

!---------------------------------------------------------------------------------
subroutine run_simple(ns, tstep, simple_a, simple_b, simple_alpha, simple_beta, &
                               simple_V0, Grl_Temp, GIS_Volume_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     Grl_Temp  	Greenland (or global) mean surface temperature [degC]
! |     ns        	Number of timesteps
! |
! |    Parameters:
! |     tstep     	time step
! |     a         	sensitivity of equilibrium volume Veq [m sle/degC]
! |     b         	equilibrium volume Veq [m sle] for temperature Tg = 0
! |     alpha     	sensitivity of exponential decay rate (1/tau)
! |     beta      	exponential decay rate [1 / K] at Tg = 0
! |
! |    Initial conditions:
! |     V0        	initial ice-sheet volume [m sle]
! |
! | Outputs:
! |     GIS_Volume	Volume Greenland ice sheet [m sle]
!  =========================================================================

    USE global
    USE simple

    implicit none

    integer(i4b), intent(IN) :: ns ! time series length

! parameters
    real(DP),     intent(IN) :: tstep

    real(DP),     intent(IN) :: simple_a
    real(DP),     intent(IN) :: simple_b
    real(DP),     intent(IN) :: simple_alpha
    real(DP),     intent(IN) :: simple_beta

! intial conditions    
    real(DP),     intent(IN) :: simple_V0

! input variables
    real(DP), dimension(ns), intent(IN)  :: Grl_Temp

! output variables    
    real(DP), dimension(ns), intent(OUT) :: GIS_Volume_out

    integer(i4b) :: i ! time step counter
    

! Initialize simple (parameters and initial variable values)
    i = 1
    call init_simple(tstep, simple_a, simple_b, simple_alpha, &
                          simple_beta, simple_V0, GIS_Volume_out(i) )
                                  
! estimate outputs
    do i=2,ns
    
        ! global sea level rise
        call simple_step(Grl_Temp(i-1), GIS_Volume_out(i) )
                        
    end do

    RETURN

end subroutine run_simple