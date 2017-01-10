!=================================================================================
!  Subroutines to run BRICK-TE:
! Simple model to simulate contribution of thermal expansion (TE) to global sea-
! level rise. (Grinsted, Moore and Jevrejeva 2010, equations 1/2 that were originally
! applied to global sea-level)
!
!  Modified to use first-order explicit integration forwards in time from the
!  initial condition given (_step_forward), and first-order implicit integration
!  backwards (_step_backward) in time to the beginning of the climate simulation
!=================================================================================
!
! Private parameters/variables 'globally' used within module
!
!   tstep     time step
!
!   a         sensitivity of equilibrium volume Veq [m sle/degC]
!   b         equilibrium volume Veq [m sle] for temperature Tg = 0
!   tau       time-scale of exponential decay (e-folding time) [years]
!   TE_0      initial thermal expansion
!
!   TE         current thermal expansion
!================================================================================
! Copyright 2016 Tony Wong, Alexander Bakker
! This file is part of BRICK (Building blocks for Relevant Ice and Climate
! Knowledge). BRICK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! BRICK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with BRICK.  If not, see <http://www.gnu.org/licenses/>.
!================================================================================

module brick_te

    USE global
    implicit none
    private

! parameters:
    real(DP) :: tstep
    real(DP) :: a, b, tau, TE_0

! variables
    real(DP) :: TE       ! thermal expansion

! public subroutines
    public :: brick_te_step_forward, brick_te_step_backward, init_brick_te


contains


!------------------------------------------------------------------------------
subroutine init_brick_te(time_step, equil_sensitivity, equil_T0, &
                       timescale, Initial_TE, thermal)
!  =========================================================================
! |  Initialize the BRICK-TE parameters and initial variables.                                   |
!  =========================================================================

    real(DP), intent(IN) :: time_step
    real(DP), intent(IN) :: equil_sensitivity
    real(DP), intent(IN) :: equil_T0
    real(DP), intent(IN) :: timescale
    real(DP), intent(IN) :: Initial_TE

    real(DP), intent(OUT) :: thermal


! Assign values to model parameters
    tstep = time_step
    a     = equil_sensitivity
    b     = equil_T0
    tau   = timescale
    TE_0  = Initial_TE

! Initial values
    thermal = TE_0
    TE           = thermal

end subroutine init_brick_te
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine brick_te_step_forward(Tg, thermal)
!  ==========================================================================
! | Calculate current state from previous state
! |
! | Input:
! |       Tg:     Greenland/Global mean surface temperature (degC)
! |
! | Output:
! |       thermal:    thermal expansion [m]
!  ==========================================================================

    implicit none

    real(DP), intent(IN)  :: Tg

    real(DP), intent(OUT) :: thermal

    real(DP) :: TEeq

! Start model
    TEeq  = a * Tg + b                   ! equilibrium TE

    thermal = TE + tstep * ((TEeq - TE) / tau)
    TE           = thermal


end subroutine brick_te_step_forward
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine brick_te_step_backward(Tg, thermal_current, thermal_previous)
!  ==========================================================================
! | Calculate current state from previous state
! |
! | Input:
! |       Tg:               Greenland/Global mean surface temperature (degC)
! |       thermal_current:  current thermal expansion sea-level rise [m]
! |
! | Output:
! |       thermal_previous: previous iteration's thermal expansion [m]
!  ==========================================================================

    implicit none

    real(DP), intent(IN)  :: Tg
    real(DP), intent(IN)  :: thermal_current

    real(DP), intent(OUT) :: thermal_previous

    real(DP) :: TEeq

! Start model
    TEeq  = a * Tg + b                   ! equilibrium TE

    thermal_previous = thermal_current - tstep * ((TEeq - thermal_current) / tau)

end subroutine brick_te_step_backward
!------------------------------------------------------------------------------

END MODULE brick_te
