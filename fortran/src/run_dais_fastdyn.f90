!=================================================================================
! DAIS: simple model for Antarctic ice-sheet radius/volume [m sle] (Schaffer 2014)
! Includes the fast dynamics emulator of Wong et al. (2017)
!=================================================================================

!---------------------------------------------------------------------------------
subroutine run_dais(ns, tstep, dais_parameters,                   &
                               Ant_Temp,           Ant_Sea_Level, &
                               Ant_Sur_Ocean_Temp, Ant_SL_rate,   &
                               L_Fast_Dynamics,    AIS_Radius_out,&
                               AIS_Volume_out,     AIS_disint_out)
!  ===============================================================================
! | Inputs:
! |    Variables:
! |     Ta        Antarctic mean surface temperature [degC]
! |     SL       (Global mean) sea level [m]
! |     Toc       High latitude ocean subsurface temperatures [degC]
! |     ns        Number of timesteps
! |
! |    Parameters:
! |     tstep     time step
! |
! |     b0        Undisturbed bed height at the continent center [m]
! |     slope     Slope of ice sheet bed before loading [-]
! |     mu        Profile parameter for parabolic ice surface (related to
! |               ice stress) [m0.5]
! |     h0        hr(Ta=0): Height of runoff line at Ta = 0 [m]
! |     c         Sensitivity of Height of runoff line (hr) [m/degC]
! |     P0        P(Ta=0): Annual precipitation for Ta = 0 [m (ice equivalent)]
! |     kappa     Coefficient for the exponential dependency of precipitation on Ta
! |               [degC-1]
! |     nu        Proportionality constant relating runoff decrease with height to
! |               precipitation [m^(-1/2) yr^(-1/2)]
! |     f0        Proportionality constant for ice flow at grounding line [m/yr]
! |     gamma     Power for the relation of ice flow speed to water depth [-]
! |     alpha     Partition parameter for effect of ocean subsurface temperature on
! |               ice flux [-]
! |     Tcrit     disintegration temperature [deg C]
! |     lambda    disintegration rate [m/y]
! |     L_Fast_Dynamics  logical saying whether the fast dynamics emulator is to be used
! |     Toc_0     Present-day, high latitude ocean subsurface temperature [degC]
! |
! |    Initial conditions:
! |     Rad0      Reference ice sheet radius [m]
! |     dSL0      Sea-level rate [m/yr]
! |
! |    Constants:
! |     Tf        Freecing temperature sea water [degC]
! |     rho_w     (Sea) water density [kg/m3]
! |     rho_i     Ice density [kg/m3]
! |     rho_m     Rock density [kg/m3]
! |
! | Outputs:
! |     Rad       Radius Antarctic ice sheet [m]
! |     Vais      Volume Antarctic ice sheet [m sle]
! |     disint    Volume of disintegrated ice, due to fast dynamics [m sle]
!  =========================================================================

    USE global
    USE dais

    implicit none

    integer(i4b), intent(IN) :: ns ! time series length

! parameters
    real(DP),     intent(IN) :: tstep
    real(DP), dimension(23), intent(IN) :: dais_parameters

! input variables
    logical, intent(IN) :: L_Fast_Dynamics
    real(DP), dimension(ns), intent(IN) :: Ant_Temp
    real(DP), dimension(ns), intent(IN) :: Ant_Sea_Level
    real(DP), dimension(ns), intent(IN) :: Ant_Sur_Ocean_Temp
    real(DP), dimension(ns), intent(IN) :: Ant_SL_rate

! output variables
    real(DP), dimension(ns), intent(OUT) :: AIS_Radius_out
    real(DP), dimension(ns), intent(OUT) :: AIS_Volume_out
    real(DP), dimension(ns), intent(OUT) :: AIS_disint_out

    integer(i4b) :: i ! time step counter

! Initialize dais (parameters and initial variable values)
    i = 1
    call init_dais(tstep, dais_parameters, Ant_Sea_Level(i), L_Fast_Dynamics, &
                   AIS_Radius_out(i), AIS_Volume_out(i), AIS_disint_out(i) )

! estimate outputs
    do i=2,ns

        ! global sea level rise (Note, timestep sl_rate should be i rather than i-1)
        call dais_step(Ant_Temp(i-1), Ant_Sea_Level(i-1), &
                       Ant_Sur_Ocean_Temp(i-1), Ant_SL_rate(i), &
                       AIS_Radius_out(i), AIS_Volume_out(i), AIS_disint_out(i) )
        AIS_disint_out(i) = AIS_disint_out(i)+AIS_disint_out(i-1)

    end do

    RETURN

end subroutine run_dais
