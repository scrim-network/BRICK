!    SNEASY:  Simple Nonlinear EArth SYstem model, a composite of 
!               DOECLIM, Carbon Cycle, and MOC Boxmodel
!
!    Copyright (C) 2009  Klaus Keller, Nathan Urban, and Brian Tuttle
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!    Klaus Keller, klaus@psu.edu
!    Nathan Urban, nurban@princeton.edu
!    Brian Tuttle, btuttle@psu.edu
!
!===========================================================================

SUBROUTINE init_sneasy(tab_ni, tab_nj, oceanomtab)
    
    USE global
    USE CCM

    implicit none

    integer(i4b), intent(IN) :: tab_ni
    integer(i4b), intent(IN) :: tab_nj
    real(DP), dimension(tab_ni,tab_nj), intent(IN) :: oceanomtab

    ATsize = (/ tab_ni, tab_nj /)   ! Ocean anomaly table size

    call alloc_anomtab()

    anomtable = oceanomtab

    RETURN

END SUBROUTINE init_sneasy
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE fin_sneasy()

    USE global
    USE CCM

    implicit none

    call dealloc_anomtab()

    RETURN

END SUBROUTINE fin_sneasy
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE run_sneasy_model(ns, timestep, &
                    Clim_sens, Oc_vert_diff, Aer_ampl, Soil_resp, CO2_fert, &
                    TC_diff, HS_NA_surf, init_MOC, init_CO2, carbon_emis, &
                    aer_forc, other_forc, moc, forcing, co2, ocflux, temp, oheat)
!  =========================================================================
! | Inputs:
! |     ns              = number of timesteps
! |     timestep        = size of timesteps
!!|     tab_ni          = number of rows in ocean anomaly table
!!|     tab_nj          = number of columns in ocean anomaly table
!!|     oceanomtab      = ocean anomaly table
! |     Clim_sens       = climate sensitivity
! |     Oc_vert_diff    = ocean vertical diffusivity
! |     Aer_ampl        = aerosol amplification factor
! |     Soil_resp       = soil respiration
! |     CO2_fert        = carbon fertilization
! |     TC_diff         = thermocline diffusivity
! |     HS_NA_surf      = hydrological sensitivity, surface North Atlantic
! |     init_MOC        = initial MOC strength
! |     init_CO2        = initial atmospheric CO2
! |     carbon_emis     = carbon emissions timeseries
! |     aer_forc        = aerosol forcing time series
! |     other_forc      = timeseries of non-CO2, non-aerosol forcing
! |
! | Outputs:
! |     moc             = MOC strength [Sv]
! |     forcing         = radiative forcing [W/m^2]
! |     co2             = atmospheric CO2 [ppm]
! |     ocflux          = atmosphere-ocean heat flux
! |     temp            = global mean temperature anomaly [K re. preindust.]
! |     oheat           = interior ocean heat anomaly [10^22 J]
! |
! | Internal:
! |     nonCO2_forc     = nonCO2 radiative forcing time series
! |
!  =========================================================================

    USE global
    USE sneasy
    USE doeclim
    USE CCM

    implicit none

    integer(i4b), intent(IN) :: ns
    real(DP), intent(IN) :: timestep
    real(DP), intent(IN) :: Clim_sens
    real(DP), intent(IN) :: Oc_vert_diff
    real(DP), intent(IN) :: Aer_ampl
    real(DP), intent(IN) :: Soil_resp
    real(DP), intent(IN) :: CO2_fert
    real(DP), intent(IN) :: TC_diff
    real(DP), intent(IN) :: HS_NA_surf
    real(DP), intent(IN) :: init_MOC
    real(DP), intent(IN) :: init_CO2
    real(DP), dimension(ns), intent(IN) :: carbon_emis
    real(DP), dimension(ns), intent(IN) :: aer_forc
    real(DP), dimension(ns), intent(IN) :: other_forc
    real(DP) :: nonCO2_forc
    real(DP), dimension(ns), intent(OUT) :: moc
    real(DP), dimension(ns), intent(OUT) :: forcing
    real(DP), dimension(ns), intent(OUT) :: co2
    real(DP), dimension(ns), intent(OUT) :: ocflux
    real(DP), dimension(ns), intent(OUT) :: temp
    real(DP), dimension(ns), intent(OUT) :: oheat
    integer(i4b) :: tx

! Assign global variables.
    nsteps = ns
    deltat = timestep

! Allocate SNES arrays
    call init_SNES_arrays()

    call init_SNES_params(Clim_sens, Oc_vert_diff, Soil_resp, CO2_fert, &
                          TC_diff, HS_NA_surf, init_MOC, init_CO2)

    do tx=1,ns
    
        nonCO2_forc = Aer_ampl * aer_forc(tx) + other_forc(tx) ! Estimate non-CO2 radiative forcing

        call SNES_step(tx, carbon_emis(tx), nonCO2_forc, moc(tx), &
                        forcing(tx), temp(tx))
    end do

! Assign values to output parameters from the component modules global vars.
!  CCM:
    co2 = atmco2
    ocflux = atm_oc_flux

!  DOECLIM:
    oheat = heat_interior

! Clean up.
    call dealloc_SNES()

    RETURN

END SUBROUTINE run_sneasy_model
