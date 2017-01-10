!================================================================================
! SIMPLE: simple model for Greenland ice-sheet volume [m sle] (Bakker et al 2014)
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
subroutine run_simple(ns, tstep, simple_a, simple_b, simple_alpha, simple_beta, &
                               simple_V0, Grl_Temp, simple_i0, GIS_Volume_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     Grl_Temp  	Greenland (or global) mean surface temperature [degC]
! |     ns        	Number of timesteps
! |
! |    Parameters:
! |     tstep     	time step
! |     a         	sensitivity of equilibrium volume Veq [m sle/degC]
! |     b         	equilibrium volume Veq [m sle] for temperature Tg = 0
! |     alpha     	sensitivity of exponential decay rate (1/tau)
! |     beta      	exponential decay rate [1 / K] at Tg = 0
! |
! |    Initial conditions:
! |     V0        	initial ice-sheet volume [m sle]
! |     i0          index of the "initial" year (GIS_Volume[i0]=V0)
! |
! | Outputs:
! |     GIS_Volume	Volume Greenland ice sheet [m sle]
!  =========================================================================

    USE global
    USE simple

    implicit none

    integer(i4b), intent(IN) :: ns ! time series length

! parameters
    real(DP),     intent(IN) :: tstep

    real(DP),     intent(IN) :: simple_a
    real(DP),     intent(IN) :: simple_b
    real(DP),     intent(IN) :: simple_alpha
    real(DP),     intent(IN) :: simple_beta
    real(DP),     intent(IN) :: simple_i0

! intial conditions
    real(DP),     intent(IN) :: simple_V0

! input variables
    real(DP), dimension(ns), intent(IN)  :: Grl_Temp

! output variables
    real(DP), dimension(ns), intent(OUT) :: GIS_Volume_out

    integer(i4b) :: i   ! time step counter
    integer(i4b) :: i0  ! initial condition index

! error check
    if (simple_i0 < 1) print *, 'ERROR - i0 SIMPLE < 1: ', simple_i0

! Initialize simple (parameters and initial variable values)
    i0 = simple_i0
    call init_simple(tstep, simple_a, simple_b, simple_alpha, &
                          simple_beta, simple_V0, GIS_Volume_out(i0) )

! estimate outputs

! forward integration, from i0 to end of simulation
    do i=(i0+1),ns

        ! global sea level rise
        call simple_step_forward(Grl_Temp(i-1), GIS_Volume_out(i) )

    end do

! backward integration, from i0 to beginning of simulation
    if (i0 > 1) then
        do i=i0,2,-1

            ! global sea level rise
            call simple_step_backward(Grl_Temp(i), GIS_Volume_out(i), GIS_Volume_out(i-1) )

        end do
    end if

    RETURN

end subroutine run_simple
