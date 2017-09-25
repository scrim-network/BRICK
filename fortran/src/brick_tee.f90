!=================================================================================
!  Subroutines to run BRICK-TEE:
! Explicitly estimate contribution of thermal expansion (TE) to global sea-
! level rise, given an ocean heat change.
! 
! Ben Vega-Westhoff, July 2017
!=================================================================================
!
! Private parameters/variables 'globally' used within module
!
!   c         heat capacity of conservative temperature [J/kg/K]
!   a         global ocean-averaged thermal expansion coefficient [kg/m3/K]
!   rho_te    approximate density of global ocean [kg/m3]
!   sa        global ocean surface area [m2]
!   TE_0      initial thermal expansion
!
!   TE        current thermal expansion
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

module brick_tee

    USE global
    implicit none
    private

! parameters:
    real(DP) :: c, a, rho_te, sa, TE_0

! variables
    real(DP) :: TE       ! thermal expansion

! public subroutines
    public :: brick_tee_step_forward, brick_tee_step_backward, init_brick_tee


contains


!------------------------------------------------------------------------------
subroutine init_brick_tee(h_capacity, expansion_coeff, & 
                                density0, ocsa, Initial_TE, thermal)
!  =========================================================================
! |  Initialize the BRICK-TEE parameters and initial variables.                                   |
!  =========================================================================

    real(DP), intent(IN) :: h_capacity
    real(DP), intent(IN) :: expansion_coeff
    real(DP), intent(IN) :: density0
    real(DP), intent(IN) :: ocsa
    real(DP), intent(IN) :: Initial_TE

    real(DP), intent(OUT) :: thermal
    
! Assign values to model parameters
    c       = h_capacity
    a       = expansion_coeff
    rho_te  = density0
    sa      = ocsa
    TE_0    = Initial_TE

! Initial values
    thermal = TE_0
    TE      = thermal

end subroutine init_brick_tee
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine brick_tee_step_forward(deltaH, thermal_previous, thermal_current)
!  ==========================================================================
! | Calculate current state from previous state
! |
! | Input:
! |       deltaH           :  change in ocean heat [J]     
! |       thermal_previous :  previous time step thermal expansion [m]
! |
! | Output:
! |       thermal_current  :  current time step thermal expansion [m]
!  ==========================================================================

    implicit none

    real(DP), intent(IN) :: deltaH
    real(DP), intent(IN) :: thermal_previous

    real(DP), intent(OUT) :: thermal_current


! Start model
    thermal_current = thermal_previous + deltaH * a / ( c * rho_te * rho_te * sa ) 
    TE              = thermal_current

end subroutine brick_tee_step_forward
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine brick_tee_step_backward(deltaH, thermal_current, thermal_previous)
!  ==========================================================================
! | Calculate previous state from current state
! |
! | Input:
! |       deltaH           :  change in ocean heat [J]     
! |       thermal_current :  current time step thermal expansion [m]
! |
! | Output:
! |       thermal_previous  :  previous time step thermal expansion [m]
!  ==========================================================================

    implicit none

    real(DP), intent(IN) :: deltaH
    real(DP), intent(IN) :: thermal_current

    real(DP), intent(OUT) :: thermal_previous


! Start model
    thermal_previous = thermal_current - deltaH * a / ( c * rho_te * rho_te * sa )

end subroutine brick_tee_step_backward
!------------------------------------------------------------------------------



END MODULE brick_tee
