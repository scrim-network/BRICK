!  DOECLIM:  Diffusion Ocean Energy balance CLIMate model
!
!  Copyright (C) 2007 E. Kriegler
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
! Simple climate model DOECLIM
!
! calculates sea surface and land air temperature response to radiative forcing
! based on an energy balance model with 1-D diffusion ocean
!
! Constructed by Elmar Kriegler (EK),
! Potsdam Institute for Climate Impact Research
! Date: 06.02.2005
!
! References for the historical forcing values can be found in (Reference EK05):
! Kriegler, E (2005) Imprecise probability analysis for integrated assessment
! of climate change. Ph.D. thesis. University of Potsdam, 256 pp.
! opus.kobv.de/ubp/volltexte/2005/561/
!
! Model equations are described in EK05 and (Reference TK07):
! Tanaka, K, Kriegler, E, Bruckner, T, Hooss, G, Knorr, W, Raddatz, T (2007)
! Aggregated carbon cycle, atmospheric chemistry, and climate model (ACC2):
! Description of the forward and inverse modes, Reports on Earth System Science ! 40/2007,
! Max Planck Institute for Meteorology, Hamburg, 199 pp.
! www.mpimet.mpg.de/fileadmin/publikationen/Reports/BzE_40.pdf
!
!==============================================================================
!
! Updates:
! 22.05.2007 Hammer-Hollingsworth numerical correction included (EK)
! 23.05.2007 Ocean heat uptake added (EK)
! 12.02.2008 Translated to Fortran90 (Marlos Goes <mpg14@psu.edu>)
! 15.08.2009 Written as Fortran90 module (Brian Tuttle <btuttle@psu.edu>)
!  
!==============================================================================
!
! Global Parameters:
!   ak      slope coeff. for land-sea heat exchange
!   bk      inters. coeff. for land-sea heat exch.
!   bsi     marine air warming enhancement
!   cal     heat cap. of land-troposph. system
!   cas     heat cap. of ocean ML-troposph.
!   csw     specific heat capacity of 1m^3 seawater [Wa/m^3/K]
!   deltat  time step size [years]
!   flnd    land fraction
!   fso     ocean frac. area below 60m
!   kcon    conversion factor [cm2/s->m2/a]
!   q2co    2xCo2 forcing increase [W/m^2]
!   rlam    clim sens. over land enhancement
!   zbot    depth of interior ocean
!   
!   temp_landair:       land air temperature anomaly (K)
!   temp_sst:           sea surface temperature anomaly (K)
!   heat_mixed:         mixed layer heat anomaly (10^22 J)
!   heat_interior:      interior ocean heat anomaly (10^22 J)
!   heatflux_mixed:     heat uptake of the mixed layer (W/m^2)
!   heatflux_interior:  heat uptake of the interior ocean (W/m^2)
!
!==============================================================================
MODULE doeclim

    USE global

    implicit none

    private

    real(DP), parameter :: ak   = 0.31d0
    real(DP), parameter :: bk   = 1.59d0
    real(DP), parameter :: bsi  = 1.3d0
    real(DP), parameter :: cal  = 0.52d0
    real(DP), parameter :: cas  = 7.80d0
    real(DP), parameter :: csw  = 0.13d0
    real(DP), parameter :: flnd = 0.29d0
    real(DP), parameter :: fso  = 0.95d0
    real(DP), parameter :: kcon = 3155.8d0
    real(DP), parameter :: q2co = 3.7d0
    real(DP), parameter :: rlam = 1.43d0
    real(DP), parameter :: zbot = 4000d0
    real(DP), parameter :: earth_area = 5100656.D8      ! [m^2]
    real(DP), parameter :: secs_per_Year = 31556926d0

    real(DP), dimension(:), allocatable, public :: temp_landair
    real(DP), dimension(:), allocatable, public :: temp_sst
    real(DP), dimension(:), allocatable, public :: heat_mixed
    real(DP), dimension(:), allocatable, public :: heat_interior
    real(DP), dimension(:), allocatable, public :: heatflux_mixed
    real(DP), dimension(:), allocatable, public :: heatflux_interior
    real(DP), dimension(:), allocatable, public :: QL, Q0

    real(DP), dimension(:), allocatable :: Ker
    real(DP), dimension(2,2) :: IB, Adoe
    real(DP) :: taucfl, taukls, taucfs, tauksl, taudif, taubot
    real(DP) :: powtoheat

    public :: init_doeclim_arrays, doeclimtimestep_simple
    public :: init_doeclim_parameters, dealloc_doeclim

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE init_doeclim_arrays()
!  ==========================================================================
! |  This routine allocates and initializes global arrays for DOECLIM.  It   |
! |  is separate from init_doeclim_parameters so that the parameters can be  |
! |  changed without reinitializing the arrays.                              |
!  ==========================================================================

    implicit none

! Allocate global arrays.
    call alloc_doeclim()

! Initialize global arrays to zero.
    temp_landair  = 0.0d0
    temp_sst = 0.0d0
    heat_mixed = 0.0d0
    heat_interior = 0.0d0
    heatflux_mixed = 0.0d0
    heatflux_interior = 0.0d0
    QL = 0.0d0
    Q0 = 0.0d0

    RETURN

END SUBROUTINE init_doeclim_arrays
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE init_doeclim_parameters(t2co, kappa)
!  =========================================================================
! |  Initialize variables for DOECLIM.                                      |
! |                                                                         |
! |  Input parameters:                                                      |
! |     t2co:   climate sensitivity to 2xCO2 (K); default = 3               |
! |     kappa:  vertical ocean diffusivity (cm^2 s^-1); default = 0.55      |
! |     nsteps: number of steps (length of forcing and response vectors)    |
! |                                                                         |
!  =========================================================================

    implicit none

    real(DP), intent(IN)  :: t2co, kappa
    real(DP) :: cfl, cfs, kls  !, tmean
    real(DP), dimension(nsteps) :: KT0,KTA1,KTB1,KTA2,KTB2,KTA3,KTB3
!    real(DP) :: KTB3(1:nsteps)
    real(DP) :: ocean_area, keff, cnum, cden
    real(DP), dimension(2,2) :: IBaux, Baux, Cdoe
    integer(i4b) :: i,j

! DEPENDENT MODEL PARAMETERS
      ocean_area = (1.-flnd)*earth_area
!      print *,'e area',earth_area
      powtoheat = ocean_area*secs_per_Year / 1.D22

      cnum= rlam*flnd + bsi * (1.-flnd)

      cden = rlam * flnd - ak *(rlam-bsi)


!!!!! vertical diffusivity in [m^2/a]

      keff   = kcon * kappa

!!!!! climate feedback strength over land

      cfl = flnd *cnum/cden*q2co/t2co-bk*(rlam-bsi)/cden

!!!!! climate feedback strength over ocean

      cfs = (rlam * flnd - ak / (1.-flnd) * (rlam-bsi))                  &
      * cnum / cden * q2co / t2co + rlam * flnd / (1.-flnd) * bk *       &
      (rlam - bsi) / cden

!!!!! land-sea heat exchange coefficient

      kls = bk * rlam * flnd / cden - ak * flnd * cnum                  &
      / cden * q2co / t2co

!!!!! interior ocean warming time scale

      taubot = zbot**2 / keff

!!!!! ocean heat diff. time scale

      taudif = cas**2 / csw**2 * Pi / keff

!!!!! ocean response time scale

      taucfs = cas / cfs

!!!!! land response time scale

      taucfl = cal / cfl

!!!!! sea-land heat exchange time scale

      tauksl  = (1.-flnd) * cas / kls

!!!!! land-sea heat exchange time scale

      taukls  = flnd * cal / kls


!      enddo
!   Zeroth Order

      KT0(nsteps) = 4d0-2d0*SQRT(2.)

!   First Order

      KTA1(nsteps) = -8d0*exp(-taubot/deltat) + &
                    4d0*sqrt(2.)*exp(-0.5d0*taubot/deltat)

      KTB1(nsteps) = 4d0*sqrt(Pi*taubot/deltat) * &
            (1d0+erf(sqrt(0.5d0*taubot/deltat)) - 2d0*erf(sqrt(taubot/deltat)))

!   Second order

      KTA2(nsteps) =  8d0*exp(-4.*taubot/deltat) - &
                4d0*sqrt(2.)*exp(-2.*taubot/deltat)

      KTB2(nsteps) = -8d0*sqrt(Pi*taubot/deltat) *  &
            (1.+ erf(sqrt(2.*taubot/deltat)) - 2.*erf(2.*sqrt(taubot/deltat)) )

!    Third Order

      KTA3(nsteps) = -8.*exp(-9.*taubot/deltat) + &
            4d0*sqrt(2.)*exp(-4.5*taubot/deltat)

      KTB3(nsteps) = 12.*sqrt(Pi*taubot/deltat) * &
          (1d0 +erf(sqrt(4.5*taubot/deltat)) - 2.*erf(3.*sqrt(taubot/deltat)) )

!%Hammer and Hollingsworth correction (Equation 2.3.27, TK07):
!%Switched on (To switch off, comment out lines below)
!     do i=1,N_samples_lambda_star
      Cdoe(1,1) = 1./taucfl**2+1./taukls**2                      &
          +2./taucfl/taukls+bsi/taukls/tauksl
      Cdoe(1,2) = -bsi/taukls**2-bsi/taucfl/taukls            &
      -bsi/taucfs/taukls-bsi**2/taukls/tauksl
      Cdoe(2,1) = -bsi/tauksl**2-1./taucfs/tauksl             &
       -1./taucfl/tauksl-1./taukls/tauksl
      Cdoe(2,2) =  1./taucfs**2+bsi**2/tauksl**2                &
       +2.*bsi/taucfs/tauksl+bsi/taukls/tauksl
      Cdoe=Cdoe*(deltat**2/12.)

!%------------------------------------------------------------------
!% Matrices of difference equation system B*T(i+1) = Q(i) + A*T(i)
!% T = (TL,TO)
!% (Equation A.27, EK05, or Equations 2.3.24 and 2.3.27, TK07)
      Baux(1,1) = 1. + deltat/(2.*taucfl) + deltat/(2.*taukls)
      Baux(1,2) = -deltat/(2.*taukls)*bsi
      Baux(2,1) = -deltat/(2.*tauksl)
      Baux(2,2) = 1. + deltat/(2.*taucfs) + deltat/(2.*tauksl)*bsi +    &
      2.*fso*sqrt(deltat/taudif)
      Baux=Baux+Cdoe

! Calculate inverse of B
      call migs(Baux,2,IBaux)!,indx)
!      Bdoe(:,:,i)=Baux
      IB(:,:)=IBaux(:,:)
!      Bdoe=Bdoe+Cdoe



       do i=1,nsteps-1

!  Zeroth Order

      KT0(i) = 4d0*sqrt(dble(nsteps+1-i)) - 2.*sqrt(dble(nsteps+2-i))     &
            - 2d0*sqrt(dble(nsteps-i))

!  First Order

      KTA1(i) = -8d0*sqrt(dble(nsteps+1-i)) *                              &
      exp(-taubot/deltat/(nsteps+1-i)) +                                &
      4d0*sqrt(dble(nsteps+2-i)) *exp(-taubot/deltat/(nsteps+2-i)) +           &
      4d0*sqrt(dble(nsteps-i)) *exp(-taubot/deltat/(nsteps-i))

      KTB1(i) =  4d0*sqrt(Pi*taubot/deltat) * (                        &
      erf(sqrt(taubot/deltat/(nsteps-i))) +                             &
      erf(sqrt(taubot/deltat/(nsteps+2-i))) -                           &
      2d0*erf(sqrt(taubot/deltat/(nsteps+1-i))) )

!% Second Order

      KTA2(i) =  8.*sqrt(dble(nsteps+1-i)) *                               &
      exp(-4.*taubot/deltat/(nsteps+1-i))- 4.*sqrt(dble(nsteps+2-i))*           &
      exp(-4.*taubot/deltat/(nsteps+2-i))- 4.*sqrt(dble(nsteps-i)) *            &
      exp(-4.*taubot/deltat/(nsteps-i))

      KTB2(i) = -8.*sqrt(Pi*taubot/deltat) * (                         &
      erf(2.*sqrt(taubot/deltat/(dble(nsteps-i)))) +                          &
      erf(2.*sqrt(taubot/deltat/dble(nsteps+2-i))) -                        &
       2.*erf(2.*sqrt(taubot/deltat/dble(nsteps+1-i))) )

!% Third Order

      KTA3(i) = -8.*sqrt(dble(nsteps+1-i)) *                              &
      exp(-9.*taubot/deltat/(nsteps+1.-i)) + 4.*sqrt(dble(nsteps+2-i))*        &
      exp(-9.*taubot/deltat/(nsteps+2.-i)) + 4.*sqrt(dble(nsteps-i))*          &
      exp(-9.*taubot/deltat/(nsteps-i))

      KTB3(i) = 12.*sqrt(Pi*taubot/deltat) * (                         &
      erf(3.*sqrt(taubot/deltat/(nsteps-i))) +                          &
      erf(3.*sqrt(taubot/deltat/(nsteps+2-i))) -                        &
      2.*erf(3.*sqrt(taubot/deltat/(nsteps+1-i))) )

       enddo

      Ker = KT0+KTA1+KTB1+KTA2+KTB2+KTA3+KTB3

      Adoe(1,1) = 1d0 - deltat/(2.*taucfl) - deltat/(2.*taukls)
      Adoe(1,2) =  deltat/(2.*taukls)*bsi
      Adoe(2,1) =  deltat/(2.*tauksl)
      Adoe(2,2) = 1d0 - deltat/(2.*taucfs) - deltat/(2.*tauksl)*bsi +   &
      Ker(nsteps)*fso*sqrt(deltat/taudif)
     Adoe=Adoe+Cdoe

!    print *,'Adoe=',Adoe(:,:,1)
!    Cdoe=Cdoe*(deltat**2/12.)
!    Adoe=Adoe+Cdoe
!    Bdoe=Bdoe+Cdoe

    RETURN

END SUBROUTINE init_doeclim_parameters
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE doeclimtimestep_simple(n,forcing,temp)
!  ==========================================================================
! | Simple climate model DOECLIM
! |
! | calculates sea surface and land air temperature response to radiative 
! | forcing based on an energy balance model with 1-D diffusion ocean
! |
! | *** computes single time step ***
! | *** initialize with init_doeclim ***
! | *** then iterate this function ***
! |
! | Input:
! |       n:        current time step
! |       forcing:  global radiative forcing (top of atmosphere) (W/m^2)
! |
! | Output:
! |       temp: global mean temperature anomaly (K), relative to preindustrial
! |
! | Assumptions: 
! |       land surface temperature = land air temperature
! |       mixed layer temperature  = sea surface temperatures 
! |                                = marine air temperature divided by bsi
!  ==========================================================================

    implicit none

    integer(i4b),intent(IN) :: n
    real(DP),intent(OUT) :: temp
    real(DP), dimension(2) :: DQ, DPAST, QC, DTEAUX
    real(DP), dimension(2,nsteps) :: DTE
    real(DP) :: DelQL, DelQ0
    real(DP) :: forcing
    integer(i4b) :: i

     DTE(1,:) = temp_landair
     DTE(2,:) = temp_sst
!% assume land and ocean forcings are equal to global forcing
     QL(n) = forcing !-forcing(1)!-1.4D-1!-FORC(1)    !forcing
     Q0(n) = forcing !-forcing(1)!-1.4D-1!-FORC(1)    !forcing

   if (n.gt.1) then
     DelQL = QL(n) - QL(n-1)
     DelQ0 = Q0(n) - Q0(n-1)

!    % Assumption: linear forcing change between n and n+1
!     do i=1,nsteps
     QC(1) = (DelQL/cal*(1./taucfl+1./taukls)-bsi*DelQ0/cas/taukls)
     QC(2) = (DelQ0/cas*(1./taucfs+bsi/tauksl)-DelQL/cal/tauksl)
!     enddo
     QC = QC* deltat**2/12.
!    % -------------------------- INITIAL CONDITIONS ------------------------
!    % Initialization of temperature and forcing vector:
!    %Factor 1/2 in front of Q in Equation A.27, EK05, and Equation 2.3.27, TK07 is a typo!
!    %Assumption: linear forcing change between n and n+1
     DQ(1) = 0.5d0*deltat/cal*(QL(n)+QL(n-1))
     DQ(2) = 0.5d0*deltat/cas*(Q0(n)+Q0(n-1))
     DQ = DQ + QC

!    % -------------- SOLVE MODEL ------------------------------------
!    % Calculate temperatures
!     DPAST = zeros(2,1)
     DPAST = 0.0d0
     do i=1,n-1
        DPAST(2) = DPAST(2)+DTE(2,i)*Ker(nsteps-n+i)
     enddo
     DPAST(2) = DPAST(2)*fso * sqrt(deltat/taudif)
!     DPAST(2) = fso * sqrt(deltat/taudif) * DTE(2,(1:n-1))*Ker(nsteps-n+1:nsteps-1)  !was transposed

!     DTE(:,n) = IB * ( DQ + DPAST + A*DTE(:,n-1) )

     DTEAUX(1) = Adoe(1,1)*DTE(1,n-1)+Adoe(1,2)*DTE(2,n-1)
     DTEAUX(2) = Adoe(2,1)*DTE(1,n-1)+Adoe(2,2)*DTE(2,n-1)

     DTE(1,n) = IB(1,1)*(DQ(1)+DPAST(1)+DTEAUX(1))+                  &
                IB(1,2)*(DQ(2)+DPAST(2)+DTEAUX(2))
     DTE(2,n) = IB(2,1)*(DQ(1)+DPAST(1)+DTEAUX(1))+                  &
                IB(2,2)*(DQ(2)+DPAST(2)+DTEAUX(2))

     temp_landair(n) = DTE(1,n)
     temp_sst(n) = DTE(2,n)

!% Calculate ocean heat uptake [W/m^2]
!% heatflux(n) captures in the heat flux in the period between n-1 and n
!% Numerical implementation of Equation 2.7, EK05, or Equation 2.3.13, TK07)
!% ------------------------------------------------------------------------

     heatflux_mixed(n) = cas*( DTE(2,n)-DTE(2,n-1) )

!     heatflux_interior(n) = cas*fso/sqrt(taudif*deltat)*            &
!     ( 2.*DTE(2,n) - DTE(2,(1:n-1))*Ker(nsteps-n+2:nsteps) )  !was transposed
!     heatflux_interior(1)=0.
     do i=1,n-1
        heatflux_interior(n) = heatflux_interior(n)+DTE(2,i)*Ker(nsteps-n+1+i)
     enddo
     heatflux_interior(n) = cas*fso/sqrt(taudif*deltat)*(2.*DTE(2,n) -       &
                            heatflux_interior(n))

     heat_mixed(n) = heat_mixed(n-1) +heatflux_mixed(n) *(powtoheat*deltat)

     heat_interior(n) = heat_interior(n-1) + heatflux_interior(n) *      &
                        (fso*powtoheat*deltat)

   endif
   
   temp = flnd*temp_landair(n) + (1.-flnd)*bsi*temp_sst(n)


END SUBROUTINE doeclimtimestep_simple
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE MIGS(FV,N,X)!,INDX)
!  ==========================================================================
! |  Subroutine to invert matrix A(N,N) with the inverse stored
! |  in X(N,N) in the output.  Copyright (c) Tao Pang 2001.
! |
! |  Updated 10/24/2001.
! |
! |  ----------------------   Program 4.4   -----------------------------
! |
! |*************************************************************************
! |*                                                                       *
! |* Please Note:                                                          *
! |*                                                                       *
! |* (1) This computer program is written by Tao Pang in conjunction with  *
! |*     his book, "An Introduction to Computational Physics," published   *
! |*     by Cambridge University Press in 1997.                            *
! |*                                                                       *
! |* (2) No warranties, express or implied, are made for this program.     *
! |*                                                                       *
! |*************************************************************************
! |
!  ==========================================================================

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    INTEGER :: I,J,K
    INTEGER, DIMENSION(N) :: INDX
    real(DP), INTENT(IN), DIMENSION(N,N):: FV
    real(DP), INTENT(OUT), DIMENSION(N,N):: X
    real(DP), DIMENSION(N,N) :: VF
    real(DP), DIMENSION(N,N) :: B
!
    VF=FV
    DO I = 1, N
      DO J = 1, N
        B(I,J) = 0.
      END DO
    END DO
    DO I = 1, N
      B(I,I) = 1.
    END DO
!
!    print *,VF!,N,INDX
    CALL ELGS(VF,N,INDX)
!
    DO I = 1, N-1
      DO J = I+1, N
        DO K = 1, N
          B(INDX(J),K) = B(INDX(J),K)-VF(INDX(J),I)*B(INDX(I),K)
        END DO
      END DO
    END DO
!
    DO I = 1, N
      X(N,I) = B(INDX(N),I)/VF(INDX(N),N)
      DO J = N-1, 1, -1
        X(J,I) = B(INDX(J),I)
        DO K = J+1, N
          X(J,I) = X(J,I)-VF(INDX(J),K)*X(K,I)
        END DO
        X(J,I) =  X(J,I)/VF(INDX(J),J)
      END DO
    END DO
!    print *, X
END SUBROUTINE MIGS
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE ELGS(VF,N,INDX)
!  ==========================================================================
! |
! |  Subroutine to perform the partial-pivoting Gaussian elimination.
! |  A(N,N) is the original matrix in the input and transformed matrix
! |  plus the pivoting element ratios below the diagonal in the output.
! |  INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
! |
!  ==========================================================================

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    INTEGER :: I,J,K,ITMP
    INTEGER, INTENT(OUT), DIMENSION(N) :: INDX
    real(DP) :: C1,PI,PI1,PJ
    real(DP), INTENT(INOUT), DIMENSION(N,N) :: VF
    real(DP), DIMENSION (N) :: C
!
! Initialize the index
!
!    print *,VF
    DO I = 1, N
      INDX(I) = I
    END DO
!
! Find the rescaling factors, one from each row
!
    DO I = 1, N
      C1= 0.
      DO J = 1, N
        C1 = MAX(C1,ABS(VF(I,J)))
      END DO
      C(I) = C1
    END DO
!
! Search the pivoting (largest) element from each column
! Search the pivoting (largest) element from each column
!
    DO J = 1, N-1
      PI1 = 0.
      DO I = J, N
        PI = ABS(VF(INDX(I),J))/C(INDX(I))
        IF (PI.GT.PI1) THEN
          PI1 = PI
          K   = I
        ENDIF
      END DO
!
! Interchange the rows via INDX(N) to record pivoting order
!
      ITMP    = INDX(J)
      INDX(J) = INDX(K)
      INDX(K) = ITMP
      DO I = J+1, N
        PJ  = VF(INDX(I),J)/VF(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
        VF(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
        DO K = J+1, N
          VF(INDX(I),K) = VF(INDX(I),K)-PJ*VF(INDX(J),K)
        END DO
      END DO
    END DO
!
END SUBROUTINE ELGS
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE alloc_doeclim()

    implicit none
    
    integer(i4b) :: astat = 0

    allocate(temp_landair(nsteps), STAT=astat)
    allocate(temp_sst(nsteps), STAT=astat)
    allocate(heat_mixed(nsteps), STAT=astat)
    allocate(heat_interior(nsteps), STAT=astat)
    allocate(heatflux_mixed(nsteps), STAT=astat)
    allocate(heatflux_interior(nsteps), STAT=astat)
    allocate(QL(nsteps), STAT=astat)
    allocate(Q0(nsteps), STAT=astat)
    allocate(Ker(nsteps), STAT=astat)

    if (astat > 0) print *, "Problem allocating doeclim arrays."

    RETURN

END SUBROUTINE alloc_doeclim
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE dealloc_doeclim()

    implicit none
    
    integer(i4b) :: astat = 0

    deallocate(temp_landair, STAT=astat)
    deallocate(temp_sst, STAT=astat)
    deallocate(heat_mixed, STAT=astat)
    deallocate(heat_interior, STAT=astat)
    deallocate(heatflux_mixed, STAT=astat)
    deallocate(heatflux_interior, STAT=astat)
    deallocate(QL, STAT=astat)
    deallocate(Q0, STAT=astat)

    deallocate(Ker, STAT=astat)

    if (astat > 0) print *, "Problem deallocating doeclim arrays."

    RETURN

END SUBROUTINE dealloc_doeclim
!------------------------------------------------------------------------------

END MODULE doeclim
