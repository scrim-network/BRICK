!=================================================================================
! BRICK (TE + GSIC): 
!=================================================================================

!---------------------------------------------------------------------------------
subroutine run_brick_te_gsic(ns, tstep, &
                             te_a, te_b, te_tau, te_TE_0, &
                             gsic_beta0, gsic_V0, gsic_n, gsic_Gs0, gsic_Teq, &
                             Gl_Temp, TE_out, GSIC_out, SLR_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     Gl_Temp  	Global mean surface temperature [degC]
! |     ns        	Number of timesteps
! |
! |     tstep     	time step
! |
! |
! |    Parameters / initial conditions (TE):
! |     a         sensitivity of equilibrium TE [m/degC]
! |     b         equilibrium TE [m] for temperature Tg = 0
! |     tau       imescale (efolding time)
! |     TE_0      initial thermosteric expansion
! |
! |
! |    Parameters / initial conditions (GSIC):
! |     beta0     initial mass balance sensitivity (how long it takes GSIC to
! |               respond to increasing temps) [m/yr/C]
! |     V0        initial volume = max(Gs) [meter sle]
! |     n         exponent for area-volume scaling [-]
! |     Gs0       Gs[1]: the corrected corresponding sea-level rise in 1961 [m]
! |     Teq       equilibrium temperature (at which there is no change) [deg C]
! |   
! |
! | Outputs:
! |     TE_out	  Thermosteric expansion
! |     GSIC_out  Sea-level contribution from glaciers and small ice caps
! |
! |     SLR_out   Sea-level rise (total all components)
!  ===============================================================================

    USE global
    USE brick_te
    USE gsic_magicc


    implicit none

    integer(i4b), intent(IN) :: ns ! time series length


! parameters / initial conditions
    real(DP),     intent(IN) :: tstep

  ! Thermosteric expansion
    real(DP),     intent(IN) :: te_a
    real(DP),     intent(IN) :: te_b
    real(DP),     intent(IN) :: te_tau
    real(DP),     intent(IN) :: te_TE_0

  ! GSIC contribution
    real(DP),     intent(IN) :: gsic_beta0
    real(DP),     intent(IN) :: gsic_V0
    real(DP),     intent(IN) :: gsic_n
    real(DP),     intent(IN) :: gsic_Teq
    real(DP),     intent(IN) :: gsic_Gs0

! input variables
    real(DP), dimension(ns), intent(IN)  :: Gl_Temp

! output variables    
    real(DP), dimension(ns), intent(OUT) :: TE_out
    real(DP), dimension(ns), intent(OUT) :: GSIC_out
    real(DP), dimension(ns), intent(OUT) :: SLR_out


    integer(i4b) :: i ! time step counter
    

! Initialize brick_te (parameters and initial variable values)
    i = 1
    call init_brick_te(tstep, te_a, te_b, te_tau, te_TE_0, TE_out(i) )
    
    call init_gsic_magicc(tstep, gsic_beta0, gsic_V0,  gsic_n, &
                                 gsic_Teq,   gsic_Gs0, GSIC_out(i) )
    
    SLR_out(i) = TE_out(i) + GSIC_out(i)
                                 
                                  
! estimate outputs
    do i=2,ns
    
        call brick_te_step(   Gl_Temp(i-1), TE_out(i) )
        call gsic_magicc_step(Gl_Temp(i-1), GSIC_out(i) )
        
        SLR_out(i) = TE_out(i) + GSIC_out(i)
                        
    end do

    RETURN

end subroutine run_brick_te_gsic