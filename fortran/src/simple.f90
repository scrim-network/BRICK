!=================================================================================
!  Subroutines to run SIMPLE:
!  Simple, mechanistically motivated model of the Greenland ice sheet volume in
!  response to temperature variations (Bakker, Applegate and Keller 2016, eq 1-3)
!=================================================================================
!
! Private parameters/variables 'globally' used within module
!
!   tstep     time step
!
!   a         sensitivity of equilibrium volume Veq [m sle/degC]
!   b         equilibrium volume Veq [m sle] for temperature Tg = 0
!   alpha     sensitivity of exponential decay rate (1/tau)
!   beta      exponential decay rate [1 / K] at Tg = 0
!   V0        initial ice-sheet volume [m sle]
!
!   V         current volume
!=================================================================================

module simple

    USE global
    implicit none
    private
    
! parameters:
    real(DP) :: tstep
    real(DP) :: a, b, alpha, beta, V0

! variables
    real(DP) :: V       ! Volume ice sheet
        
! public subroutines
    public :: simple_step, init_simple


contains


!------------------------------------------------------------------------------
subroutine init_simple(time_step, equil_sensitivity, equil_T0, &
                       rate_sensitivity, rate_T0, Initial_Vol, Vol)
!  =========================================================================
! |  Initialize the SIMPLE parameters and initial variables.                                   |
!  =========================================================================

    real(DP), intent(IN) :: time_step
    real(DP), intent(IN) :: equil_sensitivity
    real(DP), intent(IN) :: equil_T0
    real(DP), intent(IN) :: rate_sensitivity
    real(DP), intent(IN) :: rate_T0
    real(DP), intent(IN) :: Initial_Vol
    
    real(DP), intent(OUT) :: Vol
    
  
! Assign values to model parameters
    tstep = time_step
    a     = equil_sensitivity
    b     = equil_T0
    alpha = rate_sensitivity
    beta  = rate_T0
    V0    = Initial_Vol

! Initial values
    Vol   = V0
    V     = Vol
    
end subroutine init_simple
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine simple_step(Tg, Vol)
!  ==========================================================================
! | Calculate current state from previous state
! |
! | Input:
! |       Tg:     Greenland/Global mean surface temperature (degC)
! |
! | Output:
! |       Vol:    Ice sheet's volume [m3]
!  ==========================================================================

    implicit none

    real(DP), intent(IN)  :: Tg

    real(DP), intent(OUT) :: Vol
    
    real(DP) :: Veq, tau
    
! Start model
    Veq  = a * Tg + b                   ! 'virtual' equilibrium volume (equation 2)
    tau  = 1. / (alpha * Tg + beta)      ! timescale                    (equation 3)

    Vol  = max(0., V + tstep * ((Veq - V) / tau)  )
    V    = Vol


end subroutine simple_step
!------------------------------------------------------------------------------

END MODULE simple