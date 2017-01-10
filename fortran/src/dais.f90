!=================================================================================
!  Subroutines to run DAIS:
!  simple model for Antarctic ice-sheet radius/volume [m sle] (Schaffer 2014)
!======================================================================================
!
! Private parameters/variables globally used within module
!
!   tstep     time step
!
!   b0        Undisturbed bed height at the continent center [m]
!   slope     Slope of ice sheet bed before loading [-]
!   mu        Profile parameter for parabolic ice surface (related to ice stress) [m0.5]
!   h0        hr(Ta=0): Height of runoff line at Ta = 0 [m]
!   c         Sensitivity of Height of runoff line (hr) [m/degC]
!   P0        P(Ta=0): Annual precipitation for Ta = 0 [m (ice equivalent)]
!   kappa     Coefficient for the exponential dependency of precipitation on Ta [degC-1]
!   nu        Proportionality constant relating runoff decrease with height to
!             precipitation [m^(-1/2) yr^(-1/2)]
!   f0        Proportionality constant for ice flow at grounding line [m/yr]
!   gamma     Power for the relation of ice flow speed to water depth [-]
!   alpha     Partition parameter for effect of ocean subsurface temperature on ice flux [-]
!   Toc_0     Present-day, high latitude ocean subsurface temperature [degC]
!   Rad0      Reference ice sheet radius [m]
!   dSLais    logical that tells if < dSL > represents contribution of
!              - all components (including AIS)        -> dSLais = 1
!              - all otherm components (excluding AIS) -> dSLais = 0
!   lf        Mean AIS fingerprint at AIS shore
!
!   Tf        Freecing temperature sea water [degC]
!   ro_w      (Sea) water density [kg/m3]
!   ro_i      Ice density [kg/m3]
!   ro_m      Rock density [kg/m3]
!   Aoc       Surface of the ocean [m2]
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

module dais

    USE global
    implicit none
    private

! parameters:
    real(DP) :: tstep

    real(DP) :: b0
    real(DP) :: slope
    real(DP) :: mu
    real(DP) :: h0
    real(DP) :: c
    real(DP) :: P0
    real(DP) :: kappa
    real(DP) :: nu
    real(DP) :: f0
    real(DP) :: gamma
    real(DP) :: alpha
    real(DP) :: Tf

    real(DP) :: Toc_0
    real(DP) :: Rad0
    real(DP) :: Aoc
    real(DP) :: includes_dSLais
    real(DP) :: lf

    real(DP) :: del
    real(DP) :: eps1
    real(DP) :: eps2

! variables
    real(DP) :: R       ! Radius ice sheet
    real(DP) :: V       ! Volume ice sheet

! public subroutines
    public :: dais_step, init_dais


contains


!------------------------------------------------------------------------------
subroutine init_dais(time_step, parameters, SL, Rad, Vol)
!  =========================================================================
! |  Initialize the DAIS parameters and initial variables.                                   |
!  =========================================================================

    real(DP), intent(IN) :: time_step
    real(DP), dimension(20), intent(IN) :: parameters
    real(DP), intent(IN)  :: SL
    real(DP), intent(OUT) :: Rad
    real(DP), intent(OUT) :: Vol

    real(DP) :: rc, rho_w, rho_i, rho_m


! Assign values to model parameters
    tstep  = time_step
    b0     = parameters(1)
    slope  = parameters(2)
    mu     = parameters(3)
    h0     = parameters(4)
    c      = parameters(5)
    P0     = parameters(6)
    kappa  = parameters(7)
    nu     = parameters(8)
    f0     = parameters(9)
    gamma  = parameters(10)
    alpha  = parameters(11)
    Tf     = parameters(12)
    rho_w  = parameters(13)
    rho_i  = parameters(14)
    rho_m  = parameters(15)
    Toc_0  = parameters(16)
    Rad0   = parameters(17)
    Aoc    = parameters(18)
    lf     = parameters(19)
    includes_dSLais = parameters(20)

! Initialize intermediate parameters
    del  = rho_w / rho_i
    eps1 = rho_i /(rho_m - rho_i)
    eps2 = rho_w /(rho_m - rho_i)

! Initial values
    R      = Rad0
    rc     = (b0 - SL)/slope
    V      = Pi * (1+eps1) * ( (8./15.) * mu**0.5 * R**2.5 - (1./3.)*slope*R**3)
    if(R > rc) then
       V   = V - Pi*eps2 * ( (2./3.)  * slope*(R**3-rc**3)-b0*(R**2-rc**2) )
    end if

    Rad = R
    Vol = V

end subroutine init_dais
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine dais_step(Ta, SL, Toc, dSL, Rad, Vol)
!  ==========================================================================
! | Calculate current state from previous state
! |
! | Input:
! |   (from timestep = i-1)
! |       Ta:     Antarctic mean surface temperature (degC)
! |       SL:     Sea level (m)
! |       Toc:    High latitude ocean subsurface temperatures [degC]
! |   (from timestep = i)
! |       dSL:    Sea level rate (m/yr)
! |
! | Output:
! |       Rad:    Ice sheet's radius [m]
! |       Vol:    Ice sheet's volume [m3]
!  ==========================================================================

    implicit none

    real(DP), intent(IN)  :: Ta
    real(DP), intent(IN)  :: SL
    real(DP), intent(IN)  :: Toc
    real(DP), intent(IN)  :: dSL

    real(DP), intent(OUT) :: Rad
    real(DP), intent(OUT) :: Vol

    real(DP) :: hr, rc, P, beta
    real(DP) :: rR, Btot
    real(DP) :: mit, F, ISO
    real(DP) :: Hw, Speed
    real(DP) :: fac
    real(DP) :: c_iso

! Start model
    hr   = h0 + c * Ta        ! equation 5
    rc   = (b0 - SL)/slope    ! application of equation 1 (paragraph after eq3)
    P    = P0 * exp(kappa*Ta) ! equation 6
    beta = nu * P**(0.5)      ! equation 7 (corrected with respect to text)

! Total mass accumulation on ice sheet (equation 8)
    if(hr > 0) then
      rR   = R - ((hr - b0 + slope*R)**2) / mu

      Btot = P * Pi * R**2 - &
        Pi * beta * (hr - b0 + slope*R) * (R*R - rR*rR) - &
        (4. * Pi * beta * mu**0.5 *   (R-rR)**2.5) / 5.  + &
        (4. * Pi * beta * mu**0.5 * R*(R-rR)**1.5) / 3.
    else
      Btot = P * Pi*R**2
    end if

! In case there is no marine ice sheet / grounding line
    F   = 0.   ! no ice flux
    ISO = 0.   ! (third term equation 14) NAME?
    fac = Pi * (1.+eps1) * (4./3. * mu**0.5 * R**1.5 - slope*R**2) ! ratio dV/dR (eq 14)

! In case there is a marine ice sheet / grounding line
    if(R > rc) then
      fac   = fac - ( 2.*pi*eps2 * (slope*R**2 - b0*R) ) ! correction fac (eq 14)

      Hw = slope*R - b0 + SL  ! equation 10

  ! Ice speed at grounding line (equation 11)
      Speed = f0 * &
        ((1.-alpha) + alpha * ((Toc - Tf)/(Toc_0 - Tf))**2) * &
        (Hw**gamma) / ( (slope*Rad0 - b0)**(gamma-1.) )
      F     = 2.*Pi*R * del * Hw * Speed   ! equation 9

      ! ISO term depends on dSL_tot (third term equation 14 !! NAME)
      c_iso = 2.*Pi*eps2* (slope*rc**2 - (b0/slope)*rc)

      ! first term is zero if dSL represents only non-AIS components (dSLais=0)
      ! second term is zero if dSL represents all components (dSLais=1)
      ISO =    includes_dSLais  *        c_iso         *  dSL +                        &!dSL = dSL_tot
         (1.-includes_dSLais) * ((1.-c_iso)/c_iso) * (dSL - lf * (Btot - F) / Aoc)  !dSL = dSL_nonAIS
    end if

    ! Ice sheet volume (equation 13)
    R      = R + tstep*(Btot-F+ISO)/fac
    V      = V + tstep*(Btot-F+ISO)

    Rad = R
    Vol = V


end subroutine dais_step
!------------------------------------------------------------------------------

END MODULE dais
