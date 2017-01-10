!=================================================================================
! GSIC_magicc: simple model for glacier and small icecap contribution to 
! SLR [m] (Wigley 2005)
!=================================================================================

!---------------------------------------------------------------------------------
subroutine run_gsic_magicc(ns, tstep, gsic_magicc_beta0, gsic_magicc_V0, &
                           gsic_magicc_n, gsic_magicc_Gs0,  gsic_magicc_Teq, &
                           Gl_Temp, SL_contribution_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     Gl_Temp  	Global mean surface temperature [degC]
! |     ns        	Number of timesteps
! |
! |    Parameters:
! |     tstep     	time step
! |     beta0     initial mass balance sensitivity (how long it takes GSIC to
! |               respond to increasing temps) [m/yr/C]
! |     V0        initial volume = max(Gs) [meter sle]
! |     n         exponent for area-volume scaling [-]
! |     Gs0       Gs[1]: the corrected corresponding sea-level rise in 1961 [m]
! |     Teq       equilibrium temperature (at which there is no change) [deg C]
! |
! | Outputs:
! |     SL_contribution_out
!  =========================================================================

    USE global
    USE gsic_magicc

    implicit none

    integer(i4b), intent(IN) :: ns ! time series length

! parameters
    real(DP),     intent(IN) :: tstep

    real(DP),     intent(IN) :: gsic_magicc_beta0
    real(DP),     intent(IN) :: gsic_magicc_V0
    real(DP),     intent(IN) :: gsic_magicc_n
    real(DP),     intent(IN) :: gsic_magicc_Teq

! intial conditions    
    real(DP),     intent(IN) :: gsic_magicc_Gs0

! input variables
    real(DP), dimension(ns), intent(IN)  :: Gl_Temp

! output variables    
    real(DP), dimension(ns), intent(OUT) :: SL_contribution_out

    integer(i4b) :: i ! time step counter
    

! Initialize gsic_magicc (parameters and initial variable values)
    i = 1
    call init_gsic_magicc(tstep, gsic_magicc_beta0, gsic_magicc_V0, &
                                 gsic_magicc_n,     gsic_magicc_Teq, &
                                 gsic_magicc_Gs0,   SL_contribution_out(i) )
                                  
! estimate outputs
    do i=2,ns
    
        ! global sea level rise
        call gsic_magicc_step(Gl_Temp(i-1), SL_contribution_out(i) )
                        
    end do

    RETURN

end subroutine run_gsic_magicc