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

module brick2

    USE global
    USE doeclim
    USE gsic_magicc
    USE simple

    implicit none
    private

! variables
    real(DP) :: tstep

! public subroutines
    public :: brick_step_forward, init_brick

contains

!===============================================================================
subroutine init_brick(nstep, tstep_in, forcing_current, &
                      S_doeclim_in, kappa_doeclim_in, &
                      temp_init_out, heatflux_mixed_init_out, heatflux_interior_init_out, &
                      beta0_gsic_magicc_in, V0_gsic_magicc_in, n_gsic_magicc_in, &
                      Teq_gsic_magicc_in, Gs0_gsic_magicc_in, sl_gsic_init_out, &
                      a_simple_in, b_simple_in, alpha_simple_in, beta_simple_in, &
                      V0_simple_in, sl_gis_init_out, vol_gis_init_out)
!  =========================================================================
!   Initialize the BRICK parameters and initial variables
!  =========================================================================

    integer(i4b), intent(IN)  :: nstep
    real(DP), intent(IN) :: tstep_in
    real(DP), intent(IN) :: forcing_current
    real(DP), intent(IN) :: S_doeclim_in
    real(DP), intent(IN) :: kappa_doeclim_in
    real(DP), intent(IN) :: beta0_gsic_magicc_in
    real(DP), intent(IN) :: V0_gsic_magicc_in
    real(DP), intent(IN) :: n_gsic_magicc_in
    real(DP), intent(IN) :: Teq_gsic_magicc_in
    real(DP), intent(IN) :: Gs0_gsic_magicc_in
    real(DP), intent(IN) :: a_simple_in
    real(DP), intent(IN) :: b_simple_in
    real(DP), intent(IN) :: alpha_simple_in
    real(DP), intent(IN) :: beta_simple_in
    real(DP), intent(IN) :: V0_simple_in

    real(DP), intent(OUT) :: temp_init_out
    real(DP), intent(OUT) :: heatflux_mixed_init_out
    real(DP), intent(OUT) :: heatflux_interior_init_out
    real(DP), intent(OUT) :: sl_gsic_init_out
    real(DP), intent(OUT) :: sl_gis_init_out
    real(DP), intent(OUT) :: vol_gis_init_out

! Assign values to model parameters, and initialize values
    tstep = tstep_in

! DOECLIM
    call init_doeclim_arrays()
    call init_doeclim_parameters(S_doeclim_in, kappa_doeclim_in)
    call doeclimtimestep_simple(nstep, forcing_current, temp_init_out)
    heatflux_mixed_init_out = heatflux_mixed(nstep)
    heatflux_interior_init_out = heatflux_interior(nstep)

! GSIC-MAGICC
    call init_gsic_magicc(tstep, beta0_gsic_magicc_in, V0_gsic_magicc_in, n_gsic_magicc_in, &
                          Teq_gsic_magicc_in, Gs0_gsic_magicc_in, sl_gsic_init_out)

! GIS-SIMPLE
    call init_simple(tstep, a_simple_in, b_simple_in, alpha_simple_in, &
                     beta_simple_in, V0_simple_in, vol_gis_init_out)
    sl_gis_init_out = V0_simple_in - vol_gis_init_out

end subroutine init_brick
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine brick_step_forward(nstep, forcing_current, &
                              Tg_previous, Tg_current, heatflux_mixed_current, heatflux_interior_current, &
                              sl_gsic_previous, sl_gsic_current, &
                              sl_gis_previous, vol_gis_previous, sl_gis_current, vol_gis_current)
!------------------------------------------------------------------------------
! Calculate current state from previous state
! 
! Input:
!  nstep                current time step
!  forcing_current      current time step radiative forcing
!  Tg_previous          previous time step surface temperature
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
    real(DP), intent(IN)  :: Tg_previous
    real(DP), intent(IN)  :: sl_gsic_previous
    real(DP), intent(IN)  :: sl_gis_previous
    real(DP), intent(IN)  :: vol_gis_previous

    real(DP), intent(OUT) :: Tg_current
    real(DP), intent(OUT) :: heatflux_mixed_current
    real(DP), intent(OUT) :: heatflux_interior_current
    real(DP), intent(OUT) :: sl_gsic_current
    real(DP), intent(OUT) :: sl_gis_current
    real(DP), intent(OUT) :: vol_gis_current

    ! Start the show.

    ! DOECLIM
    call doeclimtimestep_simple(nstep, forcing_current, Tg_current)
    heatflux_mixed_current = heatflux_mixed(nstep)
    heatflux_interior_current = heatflux_interior(nstep)

    ! GSIC-MAGICC
    call gsic_magicc_step_forward(Tg_previous, sl_gsic_previous, sl_gsic_current)

    ! GIS-SIMPLE
    call simple_step_forward(Tg_previous, vol_gis_previous, vol_gis_current)
    sl_gis_current = sl_gis_previous + (vol_gis_previous - vol_gis_current)

end subroutine brick_step_forward
!------------------------------------------------------------------------------

END MODULE brick2
