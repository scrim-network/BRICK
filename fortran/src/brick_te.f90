!=================================================================================
!  Subroutines to run BRICK-TE:
! Simple model to simulate contribution of thermosteric expansion (TE) to global sea-
! level rise. (Grinsted, Moore and Jevrejeva 2010, equations 1/2 that were originally
! applied to global sea-level)
!=================================================================================
!
! Private parameters/variables 'globally' used within module
!
!   tstep     time step
!
!   a         sensitivity of equilibrium volume Veq [m sle/degC]
!   b         equilibrium volume Veq [m sle] for temperature Tg = 0
!   tau       time-scale of exponential decay (e-folding time) [years]
!   TE_0      initial thermosteric expansion
!
!   TE         current thermosteric expansion
!=================================================================================

module brick_te

    USE global
    implicit none
    private
    
! parameters:
    real(DP) :: tstep
    real(DP) :: a, b, tau, TE_0

! variables
    real(DP) :: TE       ! thermosteric expansion
        
! public subroutines
    public :: brick_te_step, init_brick_te


contains


!------------------------------------------------------------------------------
subroutine init_brick_te(time_step, equil_sensitivity, equil_T0, &
                       timescale, Initial_TE, thermosteric)
!  =========================================================================
! |  Initialize the BRICK-TE parameters and initial variables.                                   |
!  =========================================================================

    real(DP), intent(IN) :: time_step
    real(DP), intent(IN) :: equil_sensitivity
    real(DP), intent(IN) :: equil_T0
    real(DP), intent(IN) :: timescale
    real(DP), intent(IN) :: Initial_TE
    
    real(DP), intent(OUT) :: thermosteric
    
  
! Assign values to model parameters
    tstep = time_step
    a     = equil_sensitivity
    b     = equil_T0
    tau   = timescale
    TE_0  = Initial_TE

! Initial values
    thermosteric = TE_0
    TE           = thermosteric
    
end subroutine init_brick_te
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine brick_te_step(Tg, thermosteric)
!  ==========================================================================
! | Calculate current state from previous state
! |
! | Input:
! |       Tg:     Greenland/Global mean surface temperature (degC)
! |
! | Output:
! |       thermosteric:    thermosteric expansion [m]
!  ==========================================================================

    implicit none

    real(DP), intent(IN)  :: Tg

    real(DP), intent(OUT) :: thermosteric
    
    real(DP) :: TEeq
    
! Start model
    TEeq  = a * Tg + b                   ! equilibrium TE

    thermosteric = TE + tstep * ((TEeq - TE) / tau)
    TE           = thermosteric


end subroutine brick_te_step
!------------------------------------------------------------------------------

END MODULE brick_te