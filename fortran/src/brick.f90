!=================================================================================
!  Subroutines to run BRICK
!=================================================================================
!   Description to come...
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
    real(DP) :: tstep
    real(DP) :: S_doeclim, kappa_doeclim
    real(DP) :: beta0_gsic, n_gsic, Teq_gsic, V0_gsic, Gs0_gsic
    real(DP) :: a_te, b_te, invtau_te, V0_te
    real(DP) :: a_simple, b_simple, alpha_simple, beta_simple, V0_simple
    real(DP) :: b0_dais, slope_dais, mu_dais, h0_dais, c_dais, &
                P0_dais, kappa_dais, nu_dais, f0_dais, gamma_dais, &
                alpha_dais, Tf_dais, rho_w_dais, rho_i_dais, rho_m_dais, &
                Toc0_dais, Rad0_dais, Aoc_dais, lf_dais

! variables
    real(DP) :: V_te        ! thermal expansion
    real(DP) :: V_simple        ! Greenland ice sheet volume
    real(DP) :: V_gsic        ! glacier and ice cap contribution to sea level
    real(DP) :: R_dais        ! Antarctic ice sheet radius (m)
    real(DP) :: V_dais        ! Antarctic ice sheet volume (m3)

! public subroutines
    public :: brick_step_forward, init_brick


contains


!------------------------------------------------------------------------------
subroutine init_brick(  tstep_in,
                        S_doeclim_in, kappa_doeclim_in, &
                        beta0_gsic_in, n_gsic_in, Teq_gsic_in, V0_gsic_in, Gs0_gsic_in, &
                        a_te_in, b_te_in, invtau_te_in, V0_te_in, &
                        a_simple_in, b_simple_in, alpha_simple_in, beta_simple_in, V0_simple_in, &
                        parameters_dais_in, &
                        temp_init_out, ocheat_init_out, &
                        sl_gsic_init_out, sl_te_init_out, sl_gis_init_out, sl_ais_init_out, &
                        rad_ais_init_out, vol_ais_init_out)
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
    real(DP), dimension(20), intent(IN) :: parameters

    real(DP), intent(OUT) :: temp_init_out
    real(DP), intent(OUT) :: ocheat_init_out
    real(DP), intent(OUT) :: sl_gsic_init_out
    real(DP), intent(OUT) :: sl_te_init_out
    real(DP), intent(OUT) :: sl_gis_init_out
    real(DP), intent(OUT) :: sl_ais_init_out
    real(DP), intent(OUT) :: rad_ais_init_out
    real(DP), intent(OUT) :: vol_ais_init_out


! Assign values to model parameters
    tstep = time_step

    S_doeclim       = S_doeclim_in
    kappa_doeclim   = kappa_doeclim_in

    beta0_gsic  = beta0_gsic_in
    n_gsic      = n_gsic_in
    Teq_gsic    = Teq_gsic_in
    V0_gsic     = V0_gsic_in
    Gs0_gsic    = Gs0_gsic_in

    a_te        = a_te_in
    b_te        = b_te_in
    invtau_te   = invtau_te_in
    V0_te       = V0_te_in

    a_simple    = a_simple_in
    b_simple    = b_simple_in
    alpha_simple= alpha_simple_in
    beta_simple = beta_simple_in
    V0_simple   = V0_simple_in

    b0_dais     = parameters(1)
    slope_dais  = parameters(2)
    mu_dais     = parameters(3)
    h0_dais     = parameters(4)
    c_dais      = parameters(5)
    P0_dais     = parameters(6)
    kappa_dais  = parameters(7)
    nu_dais     = parameters(8)
    f0_dais     = parameters(9)
    gamma_dais  = parameters(10)
    alpha_dais  = parameters(11)
    Tf_dais     = parameters(12)
    rho_w_dais  = parameters(13)
    rho_i_dais  = parameters(14)
    rho_m_dais  = parameters(15)
    Toc0_dais   = parameters(16)
    Rad0_dais   = parameters(17)
    Aoc_dais    = parameters(18)
    lf_dais     = parameters(19)

! Initial values
    temp_init_out   = 0.
    ocheat_init_out = 0.
    sl_gsic_init_out= Gs0_gsic
    sl_te_init_out  = V0_te
    sl_gis_init_out =
    sl_ais_init_out =

   ! TODO
   !TODO
    !TODO here now 

end subroutine init_brick
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine brick_te_step_forward(Tg, thermal)
!  ==========================================================================
! | Calculate current state from previous state
! |
! | Input:
! |       Tg:     Greenland/Global mean surface temperature (degC)
! |
! | Output:
! |       thermal:    thermal expansion [m]
!  ==========================================================================

    implicit none

    real(DP), intent(IN)  :: Tg

    real(DP), intent(OUT) :: thermal

    real(DP) :: TEeq

! Start model
    TEeq  = a * Tg + b                   ! equilibrium TE

    thermal = TE + tstep * ((TEeq - TE) / tau)
    TE           = thermal


end subroutine brick_te_step_forward
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
subroutine brick_te_step_backward(Tg, thermal_current, thermal_previous)
!  ==========================================================================
! | Calculate current state from previous state
! |
! | Input:
! |       Tg:               Greenland/Global mean surface temperature (degC)
! |       thermal_current:  current thermal expansion sea-level rise [m]
! |
! | Output:
! |       thermal_previous: previous iteration's thermal expansion [m]
!  ==========================================================================

    implicit none

    real(DP), intent(IN)  :: Tg
    real(DP), intent(IN)  :: thermal_current

    real(DP), intent(OUT) :: thermal_previous

    real(DP) :: TEeq

! Start model
    TEeq  = a * Tg + b                   ! equilibrium TE

    thermal_previous = thermal_current - tstep * ((TEeq - thermal_current) / tau)

end subroutine brick_te_step_backward
!------------------------------------------------------------------------------

END MODULE brick_te
