!=================================================================================
!  Subroutines of GSIC_magicc:
!  Simple model to simulate global contribution of Glaciers and Small Ice Caps
! (GSIC) to sea-level rise (Wigley and Raper 2005, application of equation 4/5)
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
!   beta0     initial mass balance sensitivity (how long it takes GSIC to respond to
!             increasing temps) [m/yr/C]
!   V0        initial volume = max(Gs) [meter sle]
!   n         exponent for area-volume scaling [-]
!   Gs0       Gs[i0]: the corrected corresponding sea-level rise in year given by i0 [m]
!   Teq       equilibrium temperature (at which there is no change) [deg C]
!
!   Gs        cumulative sea-level contribution since t0 (i=1)a [m]
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
    public :: gsic_magicc_step_forward, gsic_magicc_step_backward, init_gsic_magicc


contains


!------------------------------------------------------------------------------
subroutine init_gsic_magicc(time_step, SMB_sensitivity, initial_volume, &
                       area_volume_scaling_exponent, equil_temp, &
                       initial_contribution, SeaLevel)
!  =========================================================================
! |  Initialize the GSIC parameters and initial variables.                                   |
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
subroutine gsic_magicc_step_forward(Tg, SeaLevel)
!  ==========================================================================
! | Calculate current state from previous state
! | This is standard "forward Euler", first order explicit integration.
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

end subroutine gsic_magicc_step_forward
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine gsic_magicc_step_backward(Tg, SeaLevel_Current, SeaLevel_Previous)
!  ==========================================================================
! | Calculate previous state from current state
! | This is standard "backwards differentiation", first order implicit.
! |
! | Input:
! |       Tg:               Greenland/Global mean surface temperature (degC)
! |       SeaLevel_Current: GSIC sea level contribution in time step i (m SLE)
! |
! | Output:
! |       SeaLevel_Previous: sea level contribution in time step i-1 [m SLE]
!  ==========================================================================

    implicit none

    real(DP), intent(IN)  :: Tg
    real(DP), intent(IN)  :: SeaLevel_Current

    real(DP), intent(OUT) :: SeaLevel_Previous

! Start model
    SeaLevel_Previous = SeaLevel_Current - tstep * (beta0 * (Tg - Teq) * (1.-(SeaLevel_Current/V0))**n)

end subroutine gsic_magicc_step_backward
!------------------------------------------------------------------------------

END MODULE gsic_magicc
