!=================================================================================
! BRICK-TEE: ! Explicitly estimate contribution of thermal expansion
! (TE) to global sea-level rise, given an ocean heat change.
! 
! Ben Vega-Westhoff, July 2017
!================================================================================
! Copyright 2016 Tony Wong, Alexander Bakker
! This file is part of BRICK (Building blocks for Relevant Ice and
! Climate
! Knowledge). BRICK is free software: you can redistribute it and/or
! modify
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

!---------------------------------------------------------------------------------
subroutine run_brick_tee(ns, brick_tee_c, brick_tee_a, & 
                               brick_tee_rho, ocsa, brick_tee_TE_0, &
                               deltaH, brick_tee_i0, TE_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     deltaH          change in ocean heat [J]
! |     ns              Number of timesteps
! |
! |    Parameters:
! |     c               heat capacity of conservative temp [J/kg/K]
! |     a               coefficient of therm. expansion [kg/m3/degC]
! |     rho             ocean-average density [kg/m3]
! |     ocsa            ocean surface area [m2]
! |
! |    Initial conditions:
! |     TE_0        initial thermal expansion [m]
! |     i0          index of the "initial" year (TE_0 initial expansion)
! |
! | Outputs:
! |     TE_out  Thermosteric expansion [m]
!  =========================================================================

    USE global
    USE brick_tee

    implicit none

    integer(i4b), intent(IN) :: ns ! time series length

! parameters
    real(DP),     intent(IN) :: brick_tee_c
    real(DP),     intent(IN) :: brick_tee_a
    real(DP),     intent(IN) :: brick_tee_rho
    real(DP),     intent(IN) :: ocsa

! intial conditions
    real(DP),     intent(IN) :: brick_tee_TE_0
    real(DP),     intent(IN) :: brick_tee_i0

! input variables
    real(DP), dimension(ns), intent(IN)  :: deltaH

! output variables
    real(DP), dimension(ns), intent(OUT) :: TE_out

    integer(i4b) :: i   ! time step counter
    integer(i4b) :: i0  ! initial condition index

! error check
    if (brick_tee_i0 < 1) print *, 'ERROR - i0 TE < 1'

! Initialize brick_te (parameters and initial variable values)
    i0 = brick_tee_i0
    call init_brick_tee(brick_tee_c, brick_tee_a, brick_tee_rho, &
                          ocsa, brick_tee_TE_0, TE_out(i0) )
! estimate outputs

! forward integration, from i0 to end of simulation
    do i=(i0+1),ns

        ! global sea level rise
        call brick_tee_step_forward(deltaH(i-1), TE_out(i-1), TE_out(i) )
    end do

! backward integration, from i0 to beginning of simulation
    if (i0 > 1) then
        do i=i0,2,-1

            ! global sea level rise
            call brick_tee_step_backward(deltaH(i-1), TE_out(i), TE_out(i-1) )

        end do
    end if

    RETURN

end subroutine run_brick_tee

