!=================================================================================
!  Subroutines to run SIMPLE:
!  Simple, mechanistically motivated model of the Greenland ice sheet volume in
!  response to temperature variations (Bakker, Applegate and Keller 2016, eq 1-3)
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
!   alpha     sensitivity of exponential decay rate (1/tau)
!   beta      exponential decay rate [1 / K] at Tg = 0
!   V0        initial ice-sheet volume [m sle]
!
!   V         current volume
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
    public :: simple_step_forward, simple_step_backward, init_simple


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
subroutine simple_step_forward(Tg, Vol)
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


end subroutine simple_step_forward
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine simple_step_backward(Tg, Vol_Current, Vol_Previous)
!  ==========================================================================
! | Calculate previous state from current state
! |
! | Input:
! |       Tg:           Greenland/Global mean surface temperature (degC)
! |       Vol_Current:  GIS volume at same time as Tg
! |
! | Output:
! |       Vol_Previous: Ice sheet's volume at previous time  [m3]
!  ==========================================================================

    implicit none

    real(DP), intent(IN)  :: Tg
    real(DP), intent(IN)  :: Vol_Current

    real(DP), intent(OUT) :: Vol_Previous

    real(DP) :: Veq, tau

! Start model
    Veq  = a * Tg + b                   ! 'virtual' equilibrium volume (equation 2)
    tau  = 1. / (alpha * Tg + beta)      ! timescale                    (equation 3)

    Vol_Previous = max(0., Vol_Current - tstep * ((Veq - Vol_Current) / tau)  )

end subroutine simple_step_backward
!------------------------------------------------------------------------------

END MODULE simple
