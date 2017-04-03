!===============================================================================
! Fully coupled (tightly, integrated together) BRICK model in Fortran 90.
! Version without DOECLIM, sea-level rise components only. This version takes
! temperature as forcing.
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

subroutine run_brick(ns               , tstep          ,                  &
                     temp_forcing_in  ,                                   &    ! proxy DOECLIM input
                     gsic_magicc_beta0, gsic_magicc_V0 , gsic_magicc_n  , &    ! GSIC input
                     gsic_magicc_Gs0  , gsic_magicc_Teq,                  &
                     brick_te_a       , brick_te_b     , brick_te_invtau, &    ! TE input
                     brick_te_V0      ,                                   &
                     simple_a         , simple_b       , simple_alpha   , &    ! SIMPLE input
                     simple_beta      , simple_V0      ,                  &
                     anto_a           , anto_b         , slope_Ta2Tg    , &    ! DAIS input
                     intercept_Ta2Tg  , dais_b0        , dais_slope     , &
                     dais_mu          , dais_h0        , dais_c         , &
                     dais_chr         , dais_P0        , dais_kappa     , &
                     dais_nu          , dais_f0        , dais_gamma     , &
                     dais_alpha       , dais_Tfrz      , dais_rho_w     , &
                     dais_rho_i       , dais_rho_m     , dais_Toc0      , &
                     dais_Rad0        , dais_Aoc       , dais_lf        , &
                     sl_gsic_out      , sl_te_out      , sl_gis_out     , &    ! sea level output
                     sl_ais_out       , sl_out)

!===============================================================================
! Inputs:
! Variables:
!   ns        	        Number of timesteps
!   time_in             model time points [y] (removed)
!   temp_forcing_in     temperature forcing, relative to 1850 [deg C]
! 
! Parameters:
!   tstep            	time step [y]
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
!   anto_a              sensitivity of Toc (Antarctic ocean subsurface) on Tg (global)
!   anto_b              Toc when Tg=0
!   slope_Ta2Tg         slope of regression of Tg (global) as linear function of Ta (antarctic)
!   intercept_Ta2Tg     intercept of regression of Tg as a linear function of Ta
!   dais_b0             undisturbed bed height at the continent center, m
!   dais_slope          slope of ice sheet bed before loading, -
!   dais_mu             profile parameter for parabolic ice surface (related to ice stress), m0.5
!   dais_h0             height of runoff line at AIS surface temperaure of 0 deg C, m
!   dais_c              temperaure sensitivity of runoff line height, m/deg C
!   dais_chr            runoff line height uncertainty parameter, unitless
!   dais_P0             annual precipitation for AIS surf temp Ta=0, m (ice equivalent)
!   dais_kappa          coefficient for exponential dependency of precip on Ta, 1/degC
!   dais_nu             propportionality constant relating runoff to precip, 1/m0.5 1/yr0.5
!   dais_f0             proportionality constant for ice flow at groudning line, m/yr
!   dais_gamma          power for relation of ice flow speed to water depth, -
!   dais_alpha          partition parameter for efect of ocean subsurf temp on ice flux, -
!   dais_Tfrz           freezing temperature of sea water, deg C
!   dais_rho_w          sea water density, kg/m3
!   dais_rho_i          ice density, kg/m3
!   dais_rho_m          rock density, kg/m3
!   dais_Toc0           Antarctic ocean subsurface reference temperature, deg C
!   dais_Rad0           Antarctic ice sheet reference radius, m
!   dais_Aoc            ocean surface area, m2
!   dais_lf             AIS shore mean AIS sea level fingerprint, -
! 
! Outputs:
!   sl_te_out           sea-level rise relative to 1850 (or beg. of run) from thermal expansion [m SLE]
!   sl_gis_out          sea-level rise relative to 1850 (or beg. of run) from Greenland ice sheet [m SLE]
!   sl_ais_out          sea-level rise relative to 1850 (or beg. of run) from Antarctic ice sheet [m SLE]
!   sl_gsic_out         sea-level rise relative to 1850 (or beg. of run) from glaciers and ice caps [m SLE]
!   sl_out              sea-level rise relative to 1850 (or beg. of run), total [m]
!===============================================================================

    USE global
    USE brick

    implicit none

    integer(i4b), intent(IN) :: ns ! time series length
    real(DP), dimension(ns), intent(IN) :: temp_forcing_in

! parameters
    real(DP),     intent(IN) :: tstep
    real(DP),     intent(IN) :: gsic_magicc_beta0
    real(DP),     intent(IN) :: gsic_magicc_n
    real(DP),     intent(IN) :: gsic_magicc_V0
    real(DP),     intent(IN) :: gsic_magicc_Gs0
    real(DP),     intent(IN) :: gsic_magicc_Teq
    real(DP),     intent(IN) :: simple_a
    real(DP),     intent(IN) :: simple_b
    real(DP),     intent(IN) :: simple_alpha
    real(DP),     intent(IN) :: simple_beta
    real(DP),     intent(IN) :: simple_V0
    real(DP),     intent(IN) :: brick_te_a
    real(DP),     intent(IN) :: brick_te_b
    real(DP),     intent(IN) :: brick_te_invtau
    real(DP),     intent(IN) :: brick_te_V0
    real(DP),     intent(IN) :: anto_a
    real(DP),     intent(IN) :: anto_b
    real(DP),     intent(IN) :: slope_Ta2Tg
    real(DP),     intent(IN) :: intercept_Ta2Tg
    real(DP),     intent(IN) :: dais_b0
    real(DP),     intent(IN) :: dais_slope
    real(DP),     intent(IN) :: dais_mu
    real(DP),     intent(IN) :: dais_h0
    real(DP),     intent(IN) :: dais_c
    real(DP),     intent(IN) :: dais_chr
    real(DP),     intent(IN) :: dais_P0
    real(DP),     intent(IN) :: dais_kappa
    real(DP),     intent(IN) :: dais_nu
    real(DP),     intent(IN) :: dais_f0
    real(DP),     intent(IN) :: dais_gamma
    real(DP),     intent(IN) :: dais_alpha
    real(DP),     intent(IN) :: dais_Tfrz
    real(DP),     intent(IN) :: dais_rho_w
    real(DP),     intent(IN) :: dais_rho_i
    real(DP),     intent(IN) :: dais_rho_m
    real(DP),     intent(IN) :: dais_Toc0
    real(DP),     intent(IN) :: dais_Rad0
    real(DP),     intent(IN) :: dais_Aoc
    real(DP),     intent(IN) :: dais_lf

! output variables
    real(DP), dimension(ns), intent(OUT) :: sl_gsic_out
    real(DP), dimension(ns), intent(OUT) :: sl_gis_out
    real(DP), dimension(ns), intent(OUT) :: sl_te_out
    real(DP), dimension(ns), intent(OUT) :: sl_ais_out
    real(DP), dimension(ns), intent(OUT) :: sl_out

! internal variables
    real(DP), dimension(ns) :: vol_gis_out
    real(DP), dimension(ns) :: rad_ais_out
    real(DP), dimension(ns) :: vol_ais_out
    integer(i4b) :: i   ! time step counter
!===============================================================================

! Initialize brick (parameters and initial variable values)

    ! assign global variables
    deltat = 1.0d0
    nsteps = ns
    Tfrz = dais_Tfrz

    i=1
    call init_brick(i, tstep, &
                    gsic_magicc_beta0, gsic_magicc_V0, gsic_magicc_n, &
                    gsic_magicc_Teq, gsic_magicc_Gs0, sl_gsic_out(i), &
                    brick_te_a, brick_te_b, brick_te_invtau, &
                    brick_te_V0, sl_te_out(i), &
                    simple_a, simple_b, simple_alpha, simple_beta, &
                    simple_V0, sl_gis_out(i), vol_gis_out(i), &
                    dais_b0, dais_slope, dais_mu, dais_h0, dais_c, &
                    dais_chr, dais_P0, dais_kappa, dais_nu, dais_f0, &
                    dais_gamma, dais_alpha, dais_Tfrz, dais_rho_w, &
                    dais_rho_i, dais_rho_m, dais_Toc0, dais_Rad0, &
                    dais_Aoc, dais_lf, &
                    sl_ais_out(i), rad_ais_out(i), vol_ais_out(i), &
                    sl_out(i))

! estimate outputs

! forward integration, from beginning to end of simulation
    do i=2,ns

        ! step the model forward
        call brick_step_forward(i, temp_forcing_in(i-1), &
                                sl_gsic_out(i-1), sl_gsic_out(i), &
                                sl_te_out(i-1), sl_te_out(i), &
                                sl_gis_out(i-1), vol_gis_out(i-1), sl_gis_out(i), vol_gis_out(i), &
                                anto_a, anto_b, slope_Ta2Tg, intercept_Ta2Tg, &
                                sl_ais_out(i-1), rad_ais_out(i-1), vol_ais_out(i-1), &
                                sl_ais_out(i), rad_ais_out(i), vol_ais_out(i), &
                                sl_out(i-1), sl_out(i))

    end do

    RETURN

!===============================================================================

end subroutine run_brick
