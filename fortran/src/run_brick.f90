!===============================================================================
! Fully coupled (tightly, integrated together) BRICK model in Fortran 90.
!
! The input/output of model parameters and simulated sea-level rise, temperature
! and ocean heat uptake is not as elegant as it could be. This is quick and
! dirty testing of tight coupling of model components, for eventual integration
! of BRICK into more comprehensive assessment models. 
!
!===============================================================================
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
!===============================================================================

!===============================================================================

subroutine run_brick(ns, tstep, &
                     forcing_in, doeclim_t2co, doeclim_kappa, &                                             ! doeclim inputs
                     gsic_magicc_beta0, gsic_magicc_V0, gsic_magicc_n, gsic_magicc_Gs0,  gsic_magicc_Teq, & ! gsic-magicc inputs
                     simple_a, simple_b, simple_alpha, simple_beta, simple_V0, &                            ! gis-simple inputs
                     brick_te_a, brick_te_b, brick_te_invtau, brick_te_V0, &                               ! te inputs
                     dais_parameters, &                                                                     ! ais-dais inputs
                     time_out, temp_out, heatflux_mixed_out, heatflux_interior_out, &                       ! doeclim outputs
                     sl_te_out)                                                                             ! te outputs
                     sl_gsic_out, &                                                                         ! gsic-magicc outputs
                     sl_gis_out, GIS_Volume_out, &                                                          ! gis-simple outputs
                     sl_ais_out, AIS_Radius_out, AIS_Volume_out, &                                          ! ais-dais outputs

! these need to be incorporated as couplings within the model structure
!                               Ant_Temp,           Ant_Sea_Level, &
!                               Ant_Sur_Ocean_Temp, Ant_SL_rate,   &

!===============================================================================
! Inputs:
! Variables:
!   ns        	        Number of timesteps
!   forcing_in          radiative forcing [W/m2]
! 
! Parameters:
!   tstep            	time step [y]
! 
!   doeclim_t2co        cliamte sensitivity to doubling CO2 [deg C]
!   doeclim_kappa       ocean vertical diffusivity [cm^2 s^-1]
!   gsic_magicc_beta0   initial mass balance sensitivity (how long it takes GSIC to respond to increasing temps) [m/yr/C]
!   gsic_magicc_V0      initial volume = max(Gs) [meter sle]
!   gsic_magicc_n       exponent for area-volume scaling [-]
!   gsic_magicc_Gs0     Gs[t=1], the corrected corresponding sea-level rise in first year [m]
!   gsic_magicc_Teq     equilibrium temperature (at which there is no change) [deg C]
!   simple_a            sensitivity of equilibrium volume Veq [m sle/degC]
!   simple_b            equilibrium volume Veq [m sle] for temperature Tg = 0
!   simple_alpha        sensitivity of exponential decay rate (1/tau)
!   simple_beta         exponential decay rate [1 / K] at Tg = 0
!   simple_V0           initial ice-sheet volume [m sle]
!   brick_te_a         	sensitivity of equilibrium TE [m/degC]
!   brick_te_b         	equilibrium TE [m] for temperature Tg = 0
!   brick_te_invtau     1/timescale (efolding time) [1/y]
!   brick_te_V0        initial thermal expansion [m SLE]
!   dais_parameters     13 parameters for DAIS-ANTO model (details within the DAIS model structure)
! 
! Outputs:
!   temp_out  	        Global mean surface temperature [degC]
!   ocheat_out  	    ocean heat uptake [10^22 J]
!   sl_te_out           sea-level rise relative to 1850 (or beg. of run) from thermal expansion [m SLE]
!   sl_gis_out          sea-level rise relative to 1850 (or beg. of run) from Greenland ice sheet [m SLE]
!   sl_ais_out          sea-level rise relative to 1850 (or beg. of run) from Antarctic ice sheet [m SLE]
!   sl_gsic_out         sea-level rise relative to 1850 (or beg. of run) from glaciers and ice caps [m SLE]
!===============================================================================

    USE global

    implicit none

    integer(i4b), intent(IN) :: ns ! time series length
    real(DP), dimension(ns), intent(IN) :: forcing_in

! parameters
    real(DP),     intent(IN) :: tstep

    real(DP),     intent(IN) :: doeclim_t2co
    real(DP),     intent(IN) :: doeclim_kappa
    real(DP),     intent(IN) :: gsic_magicc_beta0
    real(DP),     intent(IN) :: gsic_magicc_n
    real(DP),     intent(IN) :: gsic_magicc_Teq
    real(DP),     intent(IN) :: simple_a
    real(DP),     intent(IN) :: simple_b
    real(DP),     intent(IN) :: simple_alpha
    real(DP),     intent(IN) :: simple_beta
    real(DP),     intent(IN) :: brick_te_a
    real(DP),     intent(IN) :: brick_te_b
    real(DP),     intent(IN) :: brick_te_invtau
    real(DP), dimension(20), intent(IN) :: dais_parameters

! intial conditions (may be viewed as parameters)
    real(DP),     intent(IN) :: gsic_magicc_V0
    real(DP),     intent(IN) :: gsic_magicc_Gs0
    real(DP),     intent(IN) :: simple_V0
    real(DP),     intent(IN) :: brick_te_V0

! output variables
    real(DP), dimension(ns), intent(OUT) :: temp_out
    real(DP), dimension(ns), intent(OUT) :: ocheat_out
    real(DP), dimension(ns), intent(OUT) :: heatflux_mixed_out
    real(DP), dimension(ns), intent(OUT) :: heatflux_interior_out
    real(DP), dimension(ns), intent(OUT) :: sl_te_out
    real(DP), dimension(ns), intent(OUT) :: sl_gis_out
    real(DP), dimension(ns), intent(OUT) :: sl_ais_out
    real(DP), dimension(ns), intent(OUT) :: sl_gsic_out

    integer(i4b) :: i   ! time step counter
!===============================================================================

! Initialize brick (parameters and initial variable values)
    i=1
    call init_brick(tstep, doeclim_t2co, doeclim_kappa, &
                    gsic_magicc_beta0, gsic_magicc_n, gsic_magicc_Teq, gsic_magicc_V0, gsic_magicc_Gs0, &
                    brick_te_a, brick_te_b, brick_te_invtau, brick_te_V0, &
                    simple_a, simple_b, simple_alpha, simple_beta, simple_V0, &
                    dais_parameters, &
                    temp_out(i), ocheat_out(i), &
                    sl_gsic_out(i), sl_te_out(i), sl_gis_out(i), sl_ais_out(i) )

! estimate outputs

! forward integration, from beginning to end of simulation
    do i=1,ns

        ! step the model forward
        call brick_step_forward(Gl_Temp(i-1), TE_out(i) )

    end do

    RETURN

!===============================================================================

end subroutine run_brick
