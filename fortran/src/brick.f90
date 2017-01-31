!===============================================================================
!  Subroutines to run BRICK
!===============================================================================
!   Description to come...
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

module brick

    USE global
    USE brick_te
    USE dais
    USE doeclim
    USE gsic_magicc
    USE simple

    implicit none
    private

!!TODO
!!TODO - use the parameters as defined in the submodel modules?
!!          they are all "private" so should not have problems stepping on toes
!!TODO
! parameters:
!    real(DP) :: tstep
!    real(DP) :: S_doeclim, kappa_doeclim
!    real(DP) :: beta0_gsic, n_gsic, Teq_gsic, V0_gsic, Gs0_gsic
!    real(DP) :: a_te, b_te, invtau_te, V0_te
!    real(DP) :: a_simple, b_simple, alpha_simple, beta_simple, V0_simple
!    real(DP) :: b0_dais, slope_dais, mu_dais, h0_dais, c_dais, &
!                P0_dais, kappa_dais, nu_dais, f0_dais, gamma_dais, &
!                alpha_dais, Tf_dais, rho_w_dais, rho_i_dais, rho_m_dais, &
!                Toc0_dais, Rad0_dais, Aoc_dais, lf_dais

! variables
    real(DP) :: V_te        ! thermal expansion
    real(DP) :: V_simple        ! Greenland ice sheet volume
    real(DP) :: V_gsic        ! glacier and ice cap contribution to sea level
    real(DP) :: R_dais        ! Antarctic ice sheet radius (m)
    real(DP) :: V_dais        ! Antarctic ice sheet volume (m3)

! public subroutines
    public :: brick_step_forward, init_brick


contains


!===============================================================================
subroutine init_brick(  tstep_in,
                        S_doeclim_in, kappa_doeclim_in, &
                        beta0_gsic_in, n_gsic_in, Teq_gsic_in, V0_gsic_in, Gs0_gsic_in, &
                        a_te_in, b_te_in, invtau_te_in, V0_te_in, &
                        a_simple_in, b_simple_in, alpha_simple_in, beta_simple_in, V0_simple_in, &
                        parameters_dais_in, &
                        temp_init_out, heatflux_mixed_init_out, heatflux_interior_init_out, &
                        sl_gsic_init_out, sl_te_init_out, sl_gis_init_out, sl_ais_init_out, &
                        vol_gis_init_out, &
                        rad_ais_init_out, vol_ais_init_out, sl_init_out)
!  =========================================================================
!   Initialize the BRICK parameters and initial variables
!  =========================================================================

    real(DP), intent(IN) :: time_step
    real(DP), intent(IN) :: S_doeclim_in
    real(DP), intent(IN) :: kappa_doeclim_in
    real(DP), intent(IN) :: beta0_gsic_in
    real(DP), intent(IN) :: n_gsic_in
    real(DP), intent(IN) :: Teq_gsic_in
    real(DP), intent(IN) :: V0_gsic_in
    real(DP), intent(IN) :: Gs0_gsic_in
    real(DP), intent(IN) :: a_te_in
    real(DP), intent(IN) :: b_te_in
    real(DP), intent(IN) :: invtau_te_in
    real(DP), intent(IN) :: V0_te_in
    real(DP), intent(IN) :: a_simple_in
    real(DP), intent(IN) :: b_simple_in
    real(DP), intent(IN) :: alpha_simple_in
    real(DP), intent(IN) :: beta_simple_in
    real(DP), intent(IN) :: V0_simple_in
    real(DP), dimension(20), intent(IN) :: parameters_dais_in

    real(DP), intent(OUT) :: temp_init_out
    real(DP), intent(OUT) :: heatflux_mixed_init_out
    real(DP), intent(OUT) :: heatflux_interior_init_out
    real(DP), intent(OUT) :: sl_gsic_init_out
    real(DP), intent(OUT) :: sl_te_init_out
    real(DP), intent(OUT) :: sl_gis_init_out
    real(DP), intent(OUT) :: sl_ais_init_out
    real(DP), intent(OUT) :: vol_gis_init_out
    real(DP), intent(OUT) :: rad_ais_init_out
    real(DP), intent(OUT) :: vol_ais_init_out
    real(DP), intent(OUT) :: sl_init_out

    real(DP) :: tau_te
    real(DP) :: sea_level_noAIS

! Assign values to model parameters, and initialize values
    tstep = time_step

! DOECLIM
    call init_doeclim_arrays()
    call init_doeclim_parameters(S_doeclim_in, kappa_doeclim_in)
    temp_init_out   = 0.
    heatflux_mixed_init_out = 0.
    heatflux_interior_init_out = 0.

! GSIC-MAGICC
    call init_gsic_magicc(  tstep, beta0_gsic_in, V0_gsic_in, n_gsic_in, &
                            Teq_gsic_in, Gs0_gsic_in, sl_gsic_init_out)

! TE
    tau_te = 1/invtau_te_in
    call init_brick_te( tstep, a_te_in, b_te_in, tau_te, V0_te_in, &
                        sl_te_init_out)

! GIS-SIMPLE
    call init_simple(   tstep, a_simple_in, b_simple_in, alpha_simple_in, &
                        beta_simple_in, V0_simple_in, vol_gis_init_out)
    sl_gis_init_out = 0.

! AIS-DAIS
    sea_level_noAIS = sl_gsic_init_out + sl_te_init_out + sl_gis_init_out
    call init_dais( tstep, parameters_dais_in, sea_level_noAIS, &
                    rad_ais_init_out, vol_ais_init_out)
    sl_ais_init_out = 0.

! Add up total inital sea level
    sl_init_out = sl_gsic_init_out + sl_te_init_out + &
                  sl_gis_init_out  + sl_ais_init_out

    RETURN

end subroutine init_brick
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine brick_step_forward(  nstep, forcing_current, Tg_current, heatflux_mixed_current, heatflux_interior_current, &
                                sl_gsic_previous, sl_gsic_current, &
                                sl_te_previous, sl_te_current, &
                                vol_gis_previous, vol_gis_current, &
                                sl_gis_previous, sl_gis_current, &
                                a_anto, b_anto, Tfrz, &
                                slope_Ta2Tg, intercept_Ta2Tg, &
                                rad_ais_current, vol_ais_previous, vol_ais_current, sl_ais_current, &
                                sl_previous, sl_current)
!------------------------------------------------------------------------------
! Calculate current state from previous state
! 
! Input:
!  nstep                current time step
!  forcing_current      current time step radiative forcing
!  sl_gsic_previous     previous time step cumulative GSIC contribution to SL [m]
!  sl_te_previous       previous time step cumulative TE contribution to SL [m]
!  vol_gis_previous     previous time step Greenland ice sheet volume [m SLE]
!  sl_gis_previous      previous time step cumulative GIS contribution to SL [m]
!  a_anto               sensitivity of Toc to Tg
!  b_anto               Toc for Tg=0 deg C
!  Tfrz                 freezing temperature of ocean water (deg C)
!  slope_Ta2Tg          slope of the regression of Tg as linear function of Ta
!  intercept_Ta2Tg      intercept of the regression of Tg as linear function of Ta
!  vol_ais_previous     previous time step Antarctice ice sheet volume [m^3]
!  sl_ais_previous      previous time step cumulative AIS contribution to SL [m]
!
! Output:
!  Tg_current           current time step surface temperature
!  heatflux_mixed_current       current time step heat flux into ocean mixed layer
!  heatflux_interior_current    current time step heat flux into ocean interior
!  sl_gsic_current      current time step cumulative GSIC contribution to SL [m]
!  sl_te_current        current time step cumulative TE contribution to SL [m]
!  vol_gis_current      current time step Greenland ice sheet volume [m SLE]
!  sl_gis_current       current time step cumulative GIS contribution to SL [m]
!  rad_ais_current      current time step Antarctic ice sheet radius [m]
!  vol_ais_current      current time step Antarctic ice sheet volume [m^3]
!  sl_ais_current       current time step cumulative AIS contribution to SL [m]
!  sl_current           current time step total sea level [m]
!
!------------------------------------------------------------------------------

    implicit none

    integer(i4b), intent(IN)  :: nstep
    real(DP), intent(IN)  :: forcing_current
    real(DP), intent(IN)  :: sl_gsic_previous
    real(DP), intent(IN)  :: sl_te_previous
    real(DP), intent(IN)  :: vol_gis_previous
    real(DP), intent(IN)  :: sl_gis_previous
    real(DP), intent(IN)  :: vol_ais_previous
    real(DP), intent(IN)  :: sl_ais_previous

    real(DP), intent(OUT) :: Tg_current
    real(DP), intent(OUT) :: sl_gsic_current
    real(DP), intent(OUT) :: sl_te_current
    real(DP), intent(OUT) :: vol_gis_current
    real(DP), intent(OUT) :: sl_gis_current
    real(DP), intent(OUT) :: rad_ais_current
    real(DP), intent(OUT) :: vol_ais_current
    real(DP), intent(OUT) :: sl_ais_current
    real(DP), intent(OUT) :: sl_current

    real(DP) :: sea_level_noAIS, change_sea_level_noAIS
    real(DP) :: Ta_current, Toc_current, ctmp, Tfrz

    ! Start the show.

    ! DOECLIM
    call doeclimtimestep_simple(nstep, forcing_current, Tg_current)

    ! GSIC-MAGICC
    call gsic_magicc_step_forward(Tg_current, sl_gsic_previous, sl_gsic_current)

    ! TE
    call brick_te_step_forward(Tg_current, sl_te_previous, sl_te_current)

    ! GIS-SIMPLE
    call simple_step_forward(Tg_current, vol_gis_previous, vol_gis_current)
    sl_gis_current = sl_gis_previous + (vol_gis_previous - vol_gis_current)

    ! AIS-DAIS
    change_sea_level_noAIS =    (sl_gsic_current - sl_gsic_previous) + &
                                (sl_te_current   - sl_te_previous)   + &
                                (sl_gis_current  - sl_gis_previous)
    sea_level_noAIS = sl_previous + change_sea_level_noAIS

    ! scale temperatures, accounting for relative to 1850 (or whenever iniital point is)
    ctmp = (Tfrz-b_anto)/a_anto
    Toc_current = Tfrz + ((a_anto*Tg_current + b_anto - Tfrz)/(1. + exp(-Tg_current + ctmp)))
    Ta_current = (Tg_current - intercept_Ta2Tg)/slope_Ta2Tg

    ! run DAIS
    call dais_step( Ta_current, sea_level_noAIS, Toc_current, &
                    change_sea_level_noAIS, rad_ais_current, vol_ais_current)
    sl_ais_current = sl_ais_previous + &
                     (57. - sl_ais_previous)*(1. - (vol_ais_current/vol_ais_previous))

    ! add up total sea level
    sl_current = sl_gsic_current + sl_te_current + sl_gis_current + sl_ais_current
    

end subroutine brick_step_forward
!------------------------------------------------------------------------------

END MODULE brick
