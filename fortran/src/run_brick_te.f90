!=================================================================================
! BRICK-TE: ! Simple model to simulate contribution of thermal expansion (TE)
! to global sea-level rise. (Grinsted, Moore and Jevrejeva 2010, equations 1/2
! that were originally applied to global sea-level)
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

!---------------------------------------------------------------------------------
subroutine run_brick_te(ns, tstep, brick_te_a, brick_te_b, brick_te_tau, &
                               brick_te_TE_0, Gl_Temp, brick_te_i0, TE_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     Gl_Temp  	Global mean surface temperature [degC]
! |     ns        	Number of timesteps
! |
! |    Parameters:
! |     tstep     	time step
! |     a         	sensitivity of equilibrium TE [m/degC]
! |     b         	equilibrium TE [m] for temperature Tg = 0
! |     tau       	timescale (efolding time)
! |
! |    Initial conditions:
! |     TE_0        initial thermal expansion
! |     i0          index of the "initial" year (TE_0 initial expansion)
! |
! | Outputs:
! |     TE_out	Thermosteric expansion
!  =========================================================================

    USE global
    USE brick_te

    implicit none

    integer(i4b), intent(IN) :: ns ! time series length

! parameters
    real(DP),     intent(IN) :: tstep

    real(DP),     intent(IN) :: brick_te_a
    real(DP),     intent(IN) :: brick_te_b
    real(DP),     intent(IN) :: brick_te_tau

! intial conditions
    real(DP),     intent(IN) :: brick_te_TE_0
    real(DP),     intent(IN) :: brick_te_i0

! input variables
    real(DP), dimension(ns), intent(IN)  :: Gl_Temp

! output variables
    real(DP), dimension(ns), intent(OUT) :: TE_out

    integer(i4b) :: i   ! time step counter
    integer(i4b) :: i0  ! initial condition index

! error check
    if (brick_te_i0 < 1) print *, 'ERROR - i0 TE < 1'

! Initialize brick_te (parameters and initial variable values)
    i0 = brick_te_i0
    call init_brick_te(tstep, brick_te_a, brick_te_b, brick_te_tau, &
                          brick_te_TE_0, TE_out(i0) )

! estimate outputs

! forward integration, from i0 to end of simulation
    do i=(i0+1),ns

        ! global sea level rise
        call brick_te_step_forward(Gl_Temp(i-1), TE_out(i) )

    end do

! backward integration, from i0 to beginning of simulation
    if (i0 > 1) then
        do i=i0,2,-1

            ! global sea level rise
            call brick_te_step_backward(Gl_Temp(i), TE_out(i), TE_out(i-1) )

        end do
    end if

    RETURN

end subroutine run_brick_te
