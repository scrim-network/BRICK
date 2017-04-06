! MOC_boxmodel:  time-dependent solution of the Zickfeld/Rahmstorf box model
!
!  Copyright (C)
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
!----------------------------------------------------------------------------
!  original author: Kirsten Zickfeld
!  obtained from Klaus Keller Oct 2004
!  Klaus received from Kirsten April 30, 2003
!  Modified 6/22/05 by K. Brennan to delay warming by 50 years and
!     add noise to overturning (=m*fact).  
!  Modified 1/5/07 by N. Urban to add noise to surface temperatures and 
!     salinity instead of to overturning.
!  Modified 1/31/07 by N. Urban to accept temperature forcing as an argument
!     rather than hardcoding a linear->constant trend.
!  Translated 8/27/08 by Brian Tuttle <btuttle@psu.edu> from Matlab to
!     Fortran 90 module.
!----------------------------------------------------------------------------
!
!
!  The model equations are (see Ocean Dyn. 54, 8 (2004), eqs. 2-9,1):
!  (only valid for positive m!):
!       dT1/step = lam1*(Tr1-T1) + m/V1*(T4-T1)
!       dT2/step = lam2*(Tr2-T2) + m/V2*(T3-T2)
!       dT3/step = lam3*(Tr3-T3) + m/V3*(T1-T3)
!       dT4/step =                 m/V4*(T2-T4)
!       dS1/step = SO*f1/V1      + m/V1*(S4-S1)
!       dS2/step = -SO*f2/V2     + m/V2*(S3-S2)
!       dS3/step = SO*(f2-f1)/V3 + m/V3*(S1-S3)
!       dS4/step =                 m/V4*(S2-S4)
!
!       m = k * (beta*(S2-S1) - alpha*(T2-T1)) ; m = flow rate
!
!   Box model parameters:
!       f1,f2 = freshwater fluxes
!       k = empirical flow constant
!       alpha, beta = expansion coefficients
!       Tri = restoring temperatures
!       lami = thermal relaxation constants
!       Vi = volume ratio of box i
!
!============================================================================
MODULE moc_boxmodel

    USE global
#ifdef MOC_NOISE
    USE rndseed
#endif

    implicit none

    private

    real(DP) :: fact, f1, f2
! Physical constants
    real(DP), parameter :: alpha = 0.00017d0       ! 1/K
    real(DP), parameter :: beta  = 0.0008d0        ! 1/psu
    real(DP), parameter :: S0    = 35.0d0          ! reference salinity, psu

! k is derived from a fit to the CLIMBER hysteresis
    real(DP), parameter :: k = 25.4d0              ! hydraulic constant, 1/yr

! Climate parameters:

!  hydrological sensitivities derived from CLIMBER-2 (Sv/K)
!  hysens = DF_i/DTreg_i according to Eq. 26-27
    real(DP), parameter :: hysens1 = -0.005d0   ! DF1/DT3
    real(DP)            :: hysens2 !=  0.013     ! DF2/DT2 (input parameter)
    real(DP), parameter :: hysens3 = 0.0d0      ! DF3/DT2 (F3 = meltwater) 0.0
    real(DP), parameter :: hysens4 = 0.0d0      ! 0.0

! regional patterns of temperature increase derived from CLIMBER-2
! pati = DTreg_i/DTglobal according to Eq. 25
    real(DP), parameter :: pat1 = 0.86d0        ! South Atlantic
    real(DP), parameter :: pat2 = 1.07d0        ! North Atlantic
    real(DP), parameter :: pat3 = 0.79d0        ! Tropical Atlantic
    real(DP), parameter :: patSH = 0.93d0       ! Southern Hemisphere
    real(DP), parameter :: patNH = 1.07d0       ! Northern Hemisphere


! model parameters:
    real(DP) :: lam1, lam2, lam3
    real(SP), parameter :: V1 = 1.1       ! volume ratio of box i w.r.t. boxvol
    real(SP), parameter :: V2 = 0.4
    real(SP), parameter :: V3 = 0.68
    real(SP), parameter :: V4 = 0.0544    ! 0.04*1.36

!  Tri are derived from a fit to the CLIMBER hysteresis
    real(DP), parameter :: Tr1 = 6.6d0          ! degC
    real(DP), parameter :: Tr2 = 2.7d0          ! degC
    real(DP), parameter :: Tr3 = 11.7d0         ! degC

    real(DP) :: sigma_T, sigma_S

! Solution vectors, previous and current time step:
    real(DP), dimension(:), allocatable :: m, T1, T2, T3, T4
    real(DP), dimension(:), allocatable :: S1, S2, S3, S4
!    real(DP) :: step

    public :: init_moc_boxmodel, init_moc_parameters, moc_mod_var_forced
    public :: dealloc_moc


CONTAINS
!-----------------------------------------------------------------------------
!SUBROUTINE init_moc_boxmodel(n_steps, t_step)
SUBROUTINE init_moc_boxmodel()
!  ===========================================================================
! | This routine allocates and initializes the arrays used in the MOC 
! | boxmodel.  It is separate from init_moc_parameters so that the parameters
! | can be changed without reinitializing the arrays.
! |
! | This subroutine calls init_moc_parameters.  To call them separately,
! | comment out that line as indicated.
! |
! | Input parameters:
! |   n_steps      :  number of time steps
! |   t_step       :  length of time step
! |
!  ===========================================================================

    implicit none

!    integer(i4b), intent(IN) :: n_steps
!    real(DP), intent(IN) :: t_step

! physical constants:
    integer(i4b), parameter :: secy = 31536000 ! seconds per year
    real(DP), parameter :: csp = 4000.0d0      ! specific heat capacity, J/K/Kg
    real(DP), parameter :: ro  = 1025.0d0      ! density, Kg/m/m/m
! Model parameters:
    real(DP), parameter :: F2_Sv = 0.065d0     ! freshwater transport, Sv
    real(DP), parameter :: gamma = 23.1d0      ! thermal coupling const. W/m^2/K
    real(SP), parameter :: boxvol = 1e17       ! reference box volume, m^3
    real(SP), parameter :: A1 = 4.6e13         ! box surface area, m^2
    real(SP), parameter :: A2 = 1e13
    real(SP), parameter :: A3 = 6.8e13
    real(SP) :: Dx1, Dx2, Dx3

!    step = t_step

    call alloc_moc(nsteps + 1)

! conversion factor of m and F from 1/yr into Sv
    fact = boxvol/1e6/real(secy,8)

! model parameters
    Dx1 = 3000             ! V1*boxvol/A1; ! depth of box i, m
    Dx2 = 3000             ! V2*boxvol/A2;
    Dx3 = V3*boxvol/A3

    lam1 = gamma*secy/csp/ro/Dx1 ! thermal relaxation constants
    lam2 = gamma*secy/csp/ro/Dx2
    lam3 = gamma*secy/csp/ro/Dx3

    f2 = F2_Sv/fact

    RETURN

END SUBROUTINE init_moc_boxmodel
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
SUBROUTINE init_moc_parameters(hysens, m0, sigma_noise_T, sigma_noise_S, moc1)
!  ===========================================================================
! |  Initialize MOC global variables.
! |
! | Input parameters:
! |   hysens: hydrological sensitivity in box 2 (surface North Atlantic)
! |            [in Sv/degrees C]
! |   m0: equilibrium MOC strength [Sv]
! |   sigma_noise_T:  standard deviation of surface box temperature
! |                   perturbations [in degrees C]
! |   sigma_noise_S:  standard deviation of surface box salinity perturbations
! |
! | Output value:
! |   moc1: initial MOC strength  [Sv]
!  ===========================================================================

    implicit none

    real(DP), intent(IN) :: hysens, m0
    real(DP), intent(IN) :: sigma_noise_T
    real(DP), intent(IN) :: sigma_noise_S
    real(DP), intent(OUT) :: moc1
! Initial condition values:
    real(DP) :: deltam, mi, DTi, fi

! physical constants:
! model parameters:
!    real(DP) :: tau1, tau2, tau3
! equilibrium values:
    real(DP) :: me
    real(DP) :: DT
!    real(DP) :: T1e, T2e, T3e, T4e
!    real(DP) :: S1e, S2e, S3e, S4e
!    real(DP) :: DS, DSN, DSS

! Define global constants.
    sigma_T = sigma_noise_T
    sigma_S = sigma_noise_S
    hysens2 = hysens

! model parameters
!    tau1 = 1.0d0/lam1      ! restoring times, yr
!    tau2 = 1.0d0/lam2
!    tau3 = 1.0d0/lam3

! time stepping
!_______________________________________________
!   delay: duration of zero temperature forcing (in years)
!ft = (length(temp_forcing(1,:))-1)/step; % forcing time in years
!tt = delay+ft;          % total time in years
!nt = tt/step;             % total time in steps
!ns = nt+1;              % total number of time steps
!wst=delay/step;  % last time step before warming begins (warming start time)
!wstart_year = temp_forcing(1,1); % year at which warming begins
!start_year = wstart_year - delay; % year at which model run begins
!t = start_year + [0:nt]*step; % times at which model is evaluated

! Equilibrium conditions:
!   Me = m0      ! default: 22.6 Sv
!   me = Me/fact
    me = m0/fact

! Compute the corresponding steady state freshwater flux f1:
!DT=(T2-T1); % accoring to Eq. 17

    DT = (me*lam1*(lam3*(Tr3-Tr1)/V2+lam2*(Tr2-Tr1)/V3)+lam1*lam2*lam3* &
    (Tr2-Tr1))/((me/V1+lam1)*(me/V2+lam2)*(me/V3+lam3)-me*me*me/V1/V2/V3)

    f1 = -me*(me+k*alpha*DT)/k/beta/S0
!    F1 = f1*fact

! Compute the other equilibrium values.

!    T1e = (me*DT/V1+lam1*Tr1)/lam1
!    T3e = (lam3*Tr3+me/V3*T1e)/(lam3+me/V3)
!    T2e = (lam2*Tr2+me/V2*T3e)/(lam2+me/V2)
!    T4e = T2e

!    S1e = S0*f1/me
!    S2e = 0.0d0
!    S3e = S0*f2/me
!    S4e = 0.0d0

!    DS = S2e - S1e
!    DSN = S3e - S2e
!    DSS = S3e - S1e

! Initial conditions
!_______________________________________________________
! Now add a perturbation deltam to the equilibrium flow
! to initialize the model away from equilibrium.

    deltam = 0.0d0     ! initial deviation of flow from equilibrium in Sv
    mi = (m0 + deltam)/fact
    DTi = (mi*lam1*(lam3*(Tr3-Tr1)/V2+lam2*(Tr2-Tr1)/V3)+lam1*lam2*lam3* &
          (Tr2-Tr1))/((mi/V1+lam1)*(mi/V2+lam2)*(mi/V3+lam3)-mi*mi*mi/V1/V2/V3)
    fi = -mi*(mi+k*alpha*DTi)/k/beta/S0

! Initialize solution vectors:
    m  = 0.0d0
    T1 = 0.0d0
    T2 = 0.0d0
    T3 = 0.0d0
    T4 = 0.0d0
    S1 = 0.0d0
    S2 = 0.0d0
    S3 = 0.0d0
    S4 = 0.0d0

! Compute initial values of T,S (assuming initial state is a steady state).
    m(1)  = mi
    T1(1) = (mi*DTi/V1+lam1*Tr1)/lam1               ! Eq. 8
    T3(1) = (lam3*Tr3+mi/V3*T1(1))/(lam3+mi/V3)     ! Eq. 9
    T2(1) = (lam2*Tr2+mi/V2*T3(1))/(lam2+mi/V2)     ! Eq. 10
    T4(1) = T2(1)
    S1(1) = S0*fi/mi
    S2(1) = 0.0d0
    S3(1) = S0*f2/mi
    S4(1) = 0.0d0

#ifdef MOC_NOISE
! Add white noise to T,S in surface boxes (NMU)
    T1(1) = T1(1) + randn()*sigma_noise_T
    T2(1) = T2(1) + randn()*sigma_noise_T
    T3(1) = T3(1) + randn()*sigma_noise_T
    S1(1) = S1(1) + randn()*sigma_noise_S
    S2(1) = S2(1) + randn()*sigma_noise_S
    S3(1) = S3(1) + randn()*sigma_noise_S
#endif

    moc1 = m(1)*fact

    RETURN

END SUBROUTINE init_moc_parameters
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
SUBROUTINE moc_mod_var_forced(tx, temp_forcing, moc)
!  ===========================================================================
! |  Single time step of the MOC box model.
! |
! |  Input:
! |    temp_forcing: temperature forcing (DegC)
! |  Output:
! |    moc: MOC strength (Sv)
! |
!  ===========================================================================

    implicit none

    integer(i4b), intent(IN) :: tx 
    real(DP), intent(IN) :: temp_forcing
    real(DP), intent(OUT) :: moc

! Solution vectors:
    real(DP) :: W1, W2, W3, W1p, W2p
    real(DP) :: h1, h2, h3, h4
    integer(i4b) :: n

!________________________________________________________
! * * * compute THC response * * *
!________________________________________________________

! define forcing time series
!________________________________________________________
 
! ... then forced as in 'temp_forcing', w/ regional pattern scaling
    W1 = pat1*temp_forcing
    W2 = pat2*temp_forcing
    W3 = pat3*temp_forcing
    W1p = patSH*temp_forcing
    W2p = patNH*temp_forcing
! KK --
!    temp_global = W1p/patSH       ! W1p*1/patSH
! KK --

!    h1 = -f1*S0     ! salt flux
!    h2 = -f2*S0     ! salt flux

! add freshwater perturbation to f1
    h1 = -S0*(hysens1*W1p/fact + f1)

! add freshwater perturbation to f2
    h2 = -S0*(hysens2*W2p/fact + f2)

! add meltwater perturbation (f3)
    h3 = -hysens3*W2p/fact*S0

! "Latif" effect
    h4 = -hysens4*W3/fact*S0

! ____________________________________________________
! compute m(t)
! ____________________________________________________

! KB changes --
! add red noise to overturning with white noise component eta, alpha=0.8
!alph=0.8;
!eta=(randn(ns,1)*sigma_noise)/fact;
!global ar_noise; 
!ar_noise=zeros(ns,1);
! KB changes --

! solve DEs:
! add white noise perturbations to surface box temperatures and salinities (NMU)
!do n=2,nt
    T1(tx+1) = T1(tx) + (lam1*(Tr1+W1-T1(tx)) + &
                m(tx)/V1*(T4(tx)-T1(tx))) * deltat
    T2(tx+1) = T2(tx) + (lam2*(Tr2+W2-T2(tx)) + &
                m(tx)/V2*(T3(tx)-T2(tx))) * deltat
    T3(tx+1) = T3(tx) + (lam3*(Tr3+W3-T3(tx)) + &
                m(tx)/V3*(T1(tx)-T3(tx))) * deltat
    T4(tx+1) = T4(tx) + (                         &
                m(tx)/V4*(T2(tx)-T4(tx))) * deltat
#ifdef MOC_NOISE    
    T1(tx+1) = T1(tx+1) + randn()*sigma_T
    T2(tx+1) = T2(tx+1) + randn()*sigma_T
    T3(tx+1) = T3(tx+1) + randn()*sigma_T
#endif

    S1(tx+1) = S1(tx) + (-h1/V1                 + &
                m(tx)/V1*(S4(tx)-S1(tx))) * deltat
    S2(tx+1) = S2(tx) + ((h2+h3)/V2             + &
                m(tx)/V2*(S3(tx)-S2(tx))) * deltat
    S3(tx+1) = S3(tx) + ((h1-h2-h4)/V3          + &
                m(tx)/V3*(S1(tx)-S3(tx))) * deltat
    S4(tx+1) = S4(tx) + (                       + &
                m(tx)/V4*(S2(tx)-S4(tx))) * deltat
#ifdef MOC_NOISE    
    S1(tx+1) = S1(tx+1) + randn()*sigma_S
    S2(tx+1) = S2(tx+1) + randn()*sigma_S
    S3(tx+1) = S3(tx+1) + randn()*sigma_S
#endif
    
    m(tx+1) = k*(beta*(S2(tx+1)-S1(tx+1)) - alpha*(T2(tx+1)-T1(tx+1)))
! KB --
!    ar_noise(n)=0+alph*(ar_noise(n-1)-0)+eta(n); % AR(1) noise zero mean, var 1 
!    m(n)=m(n)+ar_noise(n); % add ar_noise(n) to m(n)
! KB --    
    
    if (m(tx+1)<0) m(tx+1)=0

!end do

    moc = m(tx+1)*fact

! Switch indices for the next iteration.
!    cur = switch(cur)
!    prev = switch(prev)

    RETURN

END SUBROUTINE moc_mod_var_forced
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
SUBROUTINE alloc_moc(n_steps)

    implicit none

    integer(i4b), intent(IN) :: n_steps
    integer(i4b) :: astat

    astat = 0
    allocate(m(n_steps), STAT=astat)
    allocate(T1(n_steps), STAT=astat)
    allocate(T2(n_steps), STAT=astat)
    allocate(T3(n_steps), STAT=astat)
    allocate(T4(n_steps), STAT=astat)
    allocate(S1(n_steps), STAT=astat)
    allocate(S2(n_steps), STAT=astat)
    allocate(S3(n_steps), STAT=astat)
    allocate(S4(n_steps), STAT=astat)
    if (astat > 0) print *, "Problem with alloc_moc."

    RETURN

END SUBROUTINE alloc_moc
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
SUBROUTINE dealloc_moc()

    implicit none

    integer(i4b) :: astat

    astat = 0
    deallocate(m, STAT=astat)
    deallocate(T1, STAT=astat)
    deallocate(T2, STAT=astat)
    deallocate(T3, STAT=astat)
    deallocate(T4, STAT=astat)
    deallocate(S1, STAT=astat)
    deallocate(S2, STAT=astat)
    deallocate(S3, STAT=astat)
    deallocate(S4, STAT=astat)
    if (astat > 0) print *, "Problem with dealloc_moc."

    RETURN

END SUBROUTINE dealloc_moc
!-----------------------------------------------------------------------------
END MODULE moc_boxmodel
