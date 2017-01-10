!=================================================================================
! GSIC_magicc: simple model for glacier and small icecap contribution to
! SLR [m] (Wigley and Raper 2005)
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
subroutine run_gsic_magicc(ns, tstep, gsic_magicc_beta0, gsic_magicc_V0, &
                           gsic_magicc_n, gsic_magicc_Gs0,  gsic_magicc_Teq, &
                           Gl_Temp, gsic_magicc_i0, SL_contribution_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     Gl_Temp  	Global mean surface temperature [degC]
! |     ns        	Number of timesteps
! |
! |    Parameters:
! |     tstep     	time step
! |     beta0     initial mass balance sensitivity (how long it takes GSIC to
! |               respond to increasing temps) [m/yr/C]
! |     V0        initial volume = max(Gs) [meter sle]
! |     n         exponent for area-volume scaling [-]
! |     Gs0       Gs[1]: the corrected corresponding sea-level rise in 1961 [m]
! |     Teq       equilibrium temperature (at which there is no change) [deg C]
! |     i0        index of the "initial" year (Gs[i0]=Gs0, with V0 GSIC volume left)
! |
! | Outputs:
! |     SL_contribution_out
!  =========================================================================

    USE global
    USE gsic_magicc

    implicit none

    integer(i4b), intent(IN) :: ns ! time series length

! parameters
    real(DP),     intent(IN) :: tstep

    real(DP),     intent(IN) :: gsic_magicc_beta0
    real(DP),     intent(IN) :: gsic_magicc_V0
    real(DP),     intent(IN) :: gsic_magicc_n
    real(DP),     intent(IN) :: gsic_magicc_Teq
    real(DP),     intent(IN) :: gsic_magicc_i0

! intial conditions
    real(DP),     intent(IN) :: gsic_magicc_Gs0

! input variables
    real(DP), dimension(ns), intent(IN)  :: Gl_Temp

! output variables
    real(DP), dimension(ns), intent(OUT) :: SL_contribution_out

    integer(i4b) :: i   ! time step counter
    integer(i4b) :: i0  ! initial condition index

! error check
    if (gsic_magicc_i0 < 1) print *, 'ERROR - i0 GSIC < 1'

! Initialize gsic_magicc (parameters and initial variable values)
    i0 = gsic_magicc_i0
    call init_gsic_magicc(tstep, gsic_magicc_beta0, gsic_magicc_V0, &
                                 gsic_magicc_n,     gsic_magicc_Teq, &
                                 gsic_magicc_Gs0,   SL_contribution_out(i0) )

! estimate outputs

! forward integration, from i0 to end of simulation
    do i=(i0+1),ns

        ! global sea level rise
        call gsic_magicc_step_forward(Gl_Temp(i-1), SL_contribution_out(i) )

    end do

! backward integration, from i0 to beginning of simulation
    if (i0 > 1) then
        do i=i0,2,-1

            ! global sea level rise
            call gsic_magicc_step_backward(Gl_Temp(i), SL_contribution_out(i), SL_contribution_out(i-1) )

        end do
    end if

    RETURN

end subroutine run_gsic_magicc
