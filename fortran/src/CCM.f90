!  CCM:   Carbon Cycle Model
!
!  Copyright (C) 2009 D. Ricciuto, B. Tuttle, K. Keller
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!------------------------------------------------------------------------------
! Nonlinear impulse response carbon/climate model
! based on the NICCS model from Hooss et al. (2001)
! original version by DMRicciuto 7/16/2004
!
! Notes from model-code-draft-3/src/model.f90 (D. McInerney):
!   Terrestrial model:  4 box model based on Meyer et al. (1999)
!   Model to be calibrated based on flux tower data synthesis
!   Oceanic model currently using linear IRFs
!   Climate model calibrated to Hamburg AOGCM (T only)
!
! See model details in:
!   Ricciuto, D. M., K. J. Davis, and K. Keller (2008), A Bayesian calibration 
!       of a simple carbon cycle model: The role of observations in estimating 
!       and reducing uncertainty, Global Biogeochem. Cycles, 22, GB2030, 
!       doi:10.1029/2006GB002908.
!
!------------------------------------------------------------------------------
!  13 Mar 2009  Brian Tuttle <btuttle@psu.edu> received 
!               globalinversion_fortran/model.f from Dan Ricciuto.
!   May 2009    Rewrote model() subroutine as CCM.f90 module, including 
!               initialization, (de)allocation, and CC_model subroutines.
!   Aug 2009    Incorporated CCM.f90 into EarthSystem module.
!               Removed usetemp logical switch as well as the simple impulse
!               response carbon/climate model in lieue of externally 
!               computed temperature forcing.
!------------------------------------------------------------------------------

MODULE CCM

    USE global

    implicit none

    private

! Input parameters:
    real(DP) :: Cs          ! Climate sensitivity
    real(DP) :: Q10         ! Respiration Temperature sens.
    real(DP) :: Beta        ! Carbon Fertilization param.
    real(DP) :: Eta         ! Thermocline transfer velocity

    real(DP), dimension(2) :: trm
    real(DP), dimension(2) :: err
    real(DP), dimension(:,:), allocatable, public :: anomtable
    integer(i4b), dimension(2), public :: ATsize 

! Model parameters:
    real(DP) :: n2, n3, n4, hs, h1, h2, h3, h4, npp0
    real(DP), dimension(:,:), allocatable :: tpools
    real(DP), dimension(:,:), allocatable :: ocanom
    real(DP), dimension(4) :: Ftp, Goc

! Model factors:
    real(DP) :: r3f, tp1f, tp2f, tp3f
!    real(DP) :: deltat

! Climate model variables:
!    real(DP) :: a1, a2, tao1, tao2
    real(DP), dimension(:), allocatable, public :: atmco2
    real(DP), dimension(:), allocatable, public :: landflux
    real(DP), dimension(:), allocatable, public :: atm_oc_flux
    real(DP), dimension(:), allocatable, public :: emissions    ![GtC/yr]

    public :: init_CCM_arrays, init_CCM_parameters, CC_model, dealloc_CCM
    public :: alloc_anomtab, dealloc_anomtab

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE init_CCM_arrays()
!  =========================================================================
! |  Allocate and initialize Carbon Cycle Model arrays.  Set first array    |
! |  elements to preindustrial values.                                      |
!  =========================================================================

    implicit none

! Allocate global arrays.
    call alloc_CCM()

! Define factors used in the calculation.
    r3f = 45.0d0 / 120.0d0
    tp1f = 35.0d0 / 60.0d0
    tp2f = 25.0d0 / 60.0d0
    tp3f = 1.0d0 / 12.0d0

! Initialize terrestrial pools (equilibrium preindustrial values)
    tpools = 0.0d0
    Ftp = 0.0d0

    tpools(1,1) = 100.0d0     ! Non-woody vegetation [GtC]
    tpools(1,2) = 500.0d0     ! Woody vegetation [GtC]
    tpools(1,3) = 120.0d0     ! Detritus [GtC]
    tpools(1,4) = 1500.0d0    ! Soil carbon [GtC]

    landflux = 0.0d0
      
    ocanom = 0.0d0
    Goc = 0.0d0

    atmco2(:) = 0.0d0
    atmco2(1) = 285.2d0        ! [ppm]

    RETURN

END SUBROUTINE init_CCM_arrays
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE init_CCM_parameters(Clim_sens, Soil_resp, Carb_fert, Therm_diff)
!  =========================================================================
! |  Initialize the Carbon Cycle Model parameters.                          |
!  =========================================================================

    implicit none

    real(DP), intent(IN) :: Clim_sens
    real(DP), intent(IN) :: Soil_resp
    real(DP), intent(IN) :: Carb_fert
    real(DP), intent(IN) :: Therm_diff

! Assign inputs to global variables.
    Cs = Clim_sens
    Q10 = Soil_resp
    Beta = Carb_fert
    Eta = Therm_diff

! Define model parameters.
    hs = 64.0d0    !layer depths (m)
    h1 = 672.0d0
    h2 = 419.0d0
    h3 = 1136.0d0
    h4 = 2382.0d0
    n2 = Eta             !default 16.88,   diffusion coeffs [m/yr]
    n3 = 9.04d0 
    n4 = 6.32d0

    npp0 = 60.0d0             ! [GtC/yr]

! Define climate parameters.
!    a1 = 0.290d0
!    a2 = 0.710d0
!    tao1 = 400.0d0
!    tao2 = 12.0d0

! Initialize control parameters.
!    usetemp = .true.

    RETURN

END SUBROUTINE init_CCM_parameters
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE CC_model(t, temp, CO2_emissions, diagnostic)
!  =========================================================================
! |  Carbon Cycle Model, single time step
! |
! |  Input parameters:
! |     t:      time index
! |     temp:   temperature forcing [K]
! |     CO2_emissions:  [GtC/yr]
! |
! |  Optional output:
! |     diagnostic: array of model variables
!  =========================================================================
      
      implicit none

      integer(i4b), intent(IN) :: t
      real(DP), intent(IN) :: temp
      real(DP), intent(IN) :: CO2_emissions
      real(DP), dimension(:,:), intent(OUT), OPTIONAL :: diagnostic
      real(DP) :: resp3, resp4
      real(DP) :: fracinoc, netemissions
      real(DP) :: npp, resp_h
      real(DP) :: Q10temp
!      real(DP) :: forc
!      integer(i4b) :: i
     
!        if (tempobs(t+1) .lt. -900 .and. usetemp .eq. .true.) then 
!            usetemp = .false.
!        end if

    Q10temp = Q10**(temp/10.0d0)

    ! Calculate Net Primary Productivity.   (eq2, Ricciuto 2008)
        npp = npp0 * (1.0d0+Beta*log(atmco2(t)/atmco2(1)))

    ! Calculate Heterotrophic respiration rate.     (eq3, Ricciuto 2008)
        resp3 = tpools(t,3) * r3f * Q10temp
        resp4 = tpools(t,4) * 0.01d0 * Q10temp
        resp_h = resp3+resp4

        landflux(t) = resp_h - npp
	
    ! Set terrestrial pool sizes for next timestep	
    ! Prepart F(tpools,t) = d(tpools(t))/dt
        Ftp(1) = npp*tp1f - 0.35d0*tpools(t,1)
        Ftp(2) = npp*tp2f - 0.05d0*tpools(t,2)
        Ftp(3) = 0.35d0*tpools(t,1) + 0.04d0*tpools(t,2) - tp3f*tpools(t,3) - &
                    resp3
        Ftp(4) = 0.01d0*tpools(t,2) + tp3f*tpools(t,3) - resp4
      
        tpools(t+1,:) = tpools(t,:) + deltat*Ftp(:)


	    netemissions = CO2_emissions + landflux(t)
	   
    ! Find fracinoc using the ocean anomaly table.
        call anom_interp(fracinoc, temp, ocanom(t,1)+netemissions)

    
    ! Compute carbon anomaly in ocean mixed layer/atmosphere
        Goc(1) = netemissions-(n2/hs)*ocanom(t,1)* &
                            fracinoc+(n2/h2)*ocanom(t,2)
    ! Compute carbon anomaly in ocean layers 2-4
        Goc(2) = (n2/hs)*ocanom(t,1)*fracinoc- &
                            ((n2+n3)/h2)*ocanom(t,2)+(n3/h3)*ocanom(t,3)
        Goc(3) = (n3/h2)*ocanom(t,2)-((n3+n4)/h3)*ocanom(t,3)+ &
                            (n4/h4)*ocanom(t,4)
        Goc(4) = (n4/h3)*ocanom(t,3) - (n4/h4)* ocanom(t,4)

        ocanom(t+1,:) = ocanom(t,:) + deltat*Goc(:)
    
    ! Compute flux into ocean
    !    atm_oc_flux(t) = -((ocanom(t+1,1)-ocanom(t,1))*fracinoc + &
    !             ocanom(t+1,2)-ocanom(t,2)+ocanom(t+1,3)-ocanom(t,3) + &
    !                    ocanom(t+1,4)-ocanom(t,4)) / deltat
                        
    ! Compute flux into ocean 
        atm_oc_flux(t) = ((ocanom(t+1,2) + ocanom(t+1,3) + ocanom(t+1,4)) + & 
                ocanom(t+1,1)*fracinoc - ((ocanom(t,2) + ocanom(t,3) + ocanom(t,4)) + &
                        ocanom(t,1)*fracinoc)) / deltat
                        
        atmco2(t+1) = atmco2(1)+(ocanom(t+1,1)*(1-fracinoc))/2.13d0

        emissions(t) = CO2_emissions

!    ! Simple climate model
!      
!        if (usetemp .eq. .false.) then 
!            do i=1,t
!                forc = log(atmco2(i+1)/atmco2(1))-log(atmco2(i)/atmco2(1))
!                temp(t+1) = temp(t+1) + 1.0/log(2.0)*Cs*(a1*(1-exp(-(t-i) / &
!                            tao1))+a2*(1-exp(-(t-i)/tao2)))*(forc)
!            end do
!        else
!            temp(t+1)=tempobs(t+1)
!        end if

    if (present(diagnostic)) then 
        diagnostic(t+1,1) = fracinoc
        diagnostic(t+1,2) = netemissions
        diagnostic(t+1,3) = ocanom(t+1,1)
        diagnostic(t+1,4) = ocanom(t+1,2)
        diagnostic(t+1,5) = ocanom(t+1,3)
        diagnostic(t+1,6) = ocanom(t+1,4)
        diagnostic(t+1,7) = resp_h
        diagnostic(t+1,8) = npp
        diagnostic(t+1,9) = resp3
        diagnostic(t+1,10) = resp4
        diagnostic(t+1,11) = tpools(t+1,1)
        diagnostic(t+1,12) = tpools(t+1,2)
        diagnostic(t+1,13) = tpools(t+1,3)
        diagnostic(t+1,14) = tpools(t+1,4)
        diagnostic(t+1,15) = Q10temp
        diagnostic(t+1,16) = atmco2(t+1)

    end if

    RETURN

END SUBROUTINE CC_model
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE anom_interp(frac_in_ocean, ref_temp, ref_emis)
!  =========================================================================
! |  A simple linear 2D interpolation over the ocean anomaly table to find  |
! |  the fracinoc variable.                                                 |
!  =========================================================================

    implicit none

    real(DP), intent(OUT) :: frac_in_ocean
    real(DP), intent(IN) :: ref_temp, ref_emis
    integer(i4b) :: templox, temphix, emislox, emishix 
    real(DP) :: tempx, emisx, FIO_tlo, FIO_thi

! Upper and lower bound indices (integers).
     templox = min(max(floor(ref_temp*10.0d0+10.0d0)+1,1),ATsize(1))
     temphix = max(min(ceiling(ref_temp*10.0d0+10.0d0)+1,ATsize(1)),1)
     emislox = min(max(floor((ref_emis)/2.0d0)+1,1),ATsize(2))
     emishix = max(min(ceiling((ref_emis)/2.0d0)+1,ATsize(2)),1)
! Target indices (reals).
     tempx = (ref_temp*10.0d0+10.0d0)+1.0d0
     emisx = ((ref_emis)/2.0d0)+1.0d0

! First interpolate anomtable in the emission direction.
    if (emislox == emishix) then
        FIO_tlo = anomtable(templox,emislox)
        FIO_thi = anomtable(temphix,emislox)
    else
        FIO_tlo = (anomtable(templox,emishix)-anomtable(templox,emislox)) * &
         (emisx-emislox) / real(emishix-emislox) + anomtable(templox,emislox)
        FIO_thi = (anomtable(temphix,emishix)-anomtable(temphix,emislox)) * &
         (emisx-emislox) / real(emishix-emislox) + anomtable(temphix,emislox)
    end if
! Then interpolate the result in the temperature direction.
    if (templox == temphix) then
        frac_in_ocean = FIO_tlo
    else
        frac_in_ocean = (FIO_thi - FIO_tlo) * (tempx-templox) / &
            real(temphix-templox) + FIO_tlo
    end if

    RETURN

END SUBROUTINE anom_interp
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE alloc_CCM()

    implicit none

    integer(i4b) :: astat

    astat = 0
    allocate(tpools(nsteps+1,4), STAT=astat)
    allocate(ocanom(nsteps+1,4), STAT=astat)
    allocate(atmco2(nsteps+1), STAT=astat)
    allocate(landflux(nsteps), STAT=astat)
    allocate(atm_oc_flux(nsteps), STAT=astat)
    allocate(emissions(nsteps), STAT=astat)

    if (astat > 0) print *, "Problem allocating CCM arrays."

    RETURN

END SUBROUTINE alloc_CCM
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE dealloc_CCM()

    implicit none

    integer(i4b) :: astat

    astat = 0
    deallocate(tpools, STAT=astat)
    deallocate(ocanom, STAT=astat)
    deallocate(atmco2, STAT=astat)
    deallocate(landflux, STAT=astat)
    deallocate(atm_oc_flux, STAT=astat)
    deallocate(emissions, STAT=astat)

    if (astat > 0) print *, "Problem deallocating CCM arrays."

    RETURN

END SUBROUTINE dealloc_CCM
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE alloc_anomtab()

    implicit none

    integer(i4b) :: astat

    astat = 0
    allocate(anomtable(ATsize(1),ATsize(2)), STAT=astat)
    if (astat > 0) print *, "Problem allocating ocean anomaly table."

    RETURN

END SUBROUTINE alloc_anomtab
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE dealloc_anomtab()

    implicit none

    integer(i4b) :: astat

    astat = 0
    deallocate(anomtable, STAT=astat)
    if (astat > 0) print *, "Problem deallocating ocean anomaly table."

    RETURN

END SUBROUTINE dealloc_anomtab
!------------------------------------------------------------------------------

END MODULE
