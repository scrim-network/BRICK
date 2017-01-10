!=================================================================================
!  Subroutines of GSIC_magicc:
!  Simple model to simulate global contribution of Glaciers and Small Ice Caps
! (GSIC) to sea-level rise (Wigley and Raper 2005, application of equation 4/5)
!=================================================================================
!
! Private parameters/variables 'globally' used within module
!
!   tstep     time step
!
!   beta0     initial mass balance sensitivity (how long it takes GSIC to respond to
!             increasing temps) [m/yr/C]
!   V0        initial volume = max(Gs) [meter sle]
!   n         exponent for area-volume scaling [-]
!   Gs0       Gs[1]: the corrected corresponding sea-level rise in 1961 [m]
!   Teq       equilibrium temperature (at which there is no change) [deg C]
!
!   Gs        cumulative sea-level contribution since t0 (i=1)a [m]
!             by definition Gs(1) = 0
!=================================================================================


module gsic_magicc

    USE global
    implicit none
    private
    
! parameters:
    real(DP) :: tstep
    real(DP) :: beta0, V0, n, Gs0, Teq

! variables
    real(DP) :: Gs       ! Cumulative sea-level contribution
        
! public subroutines
    public :: gsic_magicc_step, init_gsic_magicc


contains


!------------------------------------------------------------------------------
subroutine init_gsic_magicc(time_step, SMB_sensitivity, initial_volume, &
                       area_volume_scaling_exponent, equil_temp, &
                       initial_contribution, SeaLevel)
!  =========================================================================
! |  Initialize the SIMPLE parameters and initial variables.                                   |
!  =========================================================================

    real(DP), intent(IN) :: time_step
    real(DP), intent(IN) :: SMB_sensitivity
    real(DP), intent(IN) :: initial_volume
    real(DP), intent(IN) :: area_volume_scaling_exponent
    real(DP), intent(IN) :: initial_contribution
    real(DP), intent(IN) :: equil_temp
    
    real(DP), intent(OUT) :: SeaLevel
    
  
! Assign values to model parameters
    tstep = time_step
    beta0 = SMB_sensitivity
    V0    = initial_volume
    n     = area_volume_scaling_exponent
    Gs0   = initial_contribution
    Teq   = equil_temp

! Initial values
    Gs       = initial_contribution
    SeaLevel = Gs
    
end subroutine init_gsic_magicc
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine gsic_magicc_step(Tg, SeaLevel)
!  ==========================================================================
! | Calculate current state from previous state
! |
! | Input:
! |       Tg:     Greenland/Global mean surface temperature (degC)
! |
! | Output:
! |       SeaLevel:  Sea level contribution [m]
!  ==========================================================================

    implicit none

    real(DP), intent(IN)  :: Tg

    real(DP), intent(OUT) :: SeaLevel
    
! Start model
    SeaLevel = Gs + tstep * (beta0 * (Tg - Teq) * (1.-(Gs/V0))**n)
    Gs       = SeaLevel

end subroutine gsic_magicc_step
!------------------------------------------------------------------------------

END MODULE gsic_magicc