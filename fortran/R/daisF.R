# =======================================================================================
# DAIS-fortran90 (# estimation by calling fortran routine)
# DAIS: Simple model for Antarctic ice-sheet volume [m sle] (Schaffer 2014)
# =======================================================================================
#
#  Requires (input variables):
#  - Ta        Antarctic mean surface temperature [degC]
#  - SL        Sea level [m]
#  - Toc       High latitude ocean subsurface temperatures [degC]
#  - dSL       Rate of sea-level change [m/yr] of either
#              - all components (including AIS)        -> includes_dSLais = 1
#              - all other components (excluding AIS)  -> includes_dSLais = 0
#
#  Simulates (output variables):
#  - Rad       Ice sheet radius [m]
#  - Vais      Volume of Antarctic ice sheet [m3]

#  Internal variables:
#  - b         Undisturbed bed profile [m]
#  - h         Ice sheet surface height [m]
#  - Btot      Total mass accumulation rate on the ice-sheet [ ??? ]
#  - F         Total ice flux across the grounding line [m3/yr]
#  - rc        Distance from the continent center to where the ice sheets enters the sea [m]
#  - hr        Height of runoff line above which precipitation accumulates as snow [m]
#  - P         Precipitation [m (ice equivalent)]
#  - beta      Mass balance gradient [ m-0.5 ]
#  - rR        Distance from the continent center to where the runoff line intersects the ice sheet surface [m]
#  - Hw        Water depth at grounding line
#
#  Coordinates:
#  - r         Radial coordinate
#
#  Parameters:
#  - b0        Undisturbed bed height at the continent center [m]
#  - slope     Slope of ice sheet bed before loading [-]
#  - mu        Profile parameter for parabolic ice sheet surface (related to ice stress) [m0.5]
#  - h0        hr(Ta=0): Height of runoff line at Ta = 0 [m]
#  - c         Sensitivity of Height of runoff line (hr) [m/degC]
#  - P0        P(Ta=0): Annual precipitation for Ta = 0 [m (ice equivalent)]
#  - kappa     Coefficient for the exponential dependency of precipitation on Ta [degC-1]
#  - nu        Proportionality constant relating runoff decrease with height to precipitation [m^(-1/2) yr^(-1/2)]
#  - f0        Proportionality constant for ice flow at grounding line [m/yr]
#  - gamma     Power for the relation of ice flow speed to water depth [-]
#  - alpha     Partition parameter for effect of ocean subsurface temperature on ice flux [-]
#  - Toc_0     Present-day, high latitude ocean subsurface temperature [degC]
#  - Rad0      Reference ice sheet radius [m]
#  - includes_dSLais    logical that tells if < dSL > represents contribution of
#              - all components (including AIS)        -> includes_dSLais = 1
#              - all otherm components (excluding AIS) -> includes_dSLais = 0
#  - tstep     time step
#
#  Constants:
#  - Tf        Freecing temperature sea water [degC]
#  - rho_w      (Sea) water density [kg/m3]
#  - rho_i      Ice density [kg/m3]
#  - rho_m      Rock density [kg/m3]
#  - Aoc        Surface of the ocean [m2]
#  - lf         Mean AIS fingerprint at AIS shore
##==============================================================================
## Copyright 2016 Tony Wong, Alexander Bakker
## This file is part of BRICK (Building blocks for Relevant Ice and Climate
## Knowledge). BRICK is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## BRICK is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with BRICK.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

# =======================================================================================
# load fortran subroutine (# to check if library is loaded is.loaded("run_dais") )
if(.Platform$OS.type == "unix") {
    dyn.load("../fortran/dais.so")
} else {
    dyn.load("../fortran/dais")
}

daisF <- function(
  tstep = 1,
  b0    = 775,
  slope = 6 * 10^(-4),
  mu    = 8.7,
  h0    = 1471,
  c     = 95,
  P0    = 0.35,
  kappa = 4 * 10^(-2),
  nu    = 1.2 * 10^(-2),
  f0    = 1.2,
  gamma = 2.5,
  alpha = 0.5,
  Tf    = -1.8,
  rho_w  = 1030,
  rho_i  = 917,
  rho_m  = 4000,
  Toc_0 = 0.72,
  Rad0  = 1.864 * 10^6,
  Aoc   = 3.619e14,
  lf    = -1.18,
  includes_dSLais= 1,
  Ta,        # Antarctic mean surface temperature [degC]
  SL,        # (global mean) sea level [m]
  Toc,       # High latitude ocean subsurface temperatures [degC]
  dSL        #
) {

  # determine series length
  ns <- length(Ta)

  parameters <- c(b0,
                  slope,
                  mu,
                  h0,
                  c,
                  P0,
                  kappa,
                  nu,
                  f0,
                  gamma,
                  alpha,
                  Tf,
                  rho_w,
                  rho_i,
                  rho_m,
                  Toc_0,
                  Rad0,
                  Aoc,
                  lf,
                  includes_dSLais)

  # call fortran
  f.output <- .Fortran("run_dais",
                  ns                 = ns,
                  tstep              = as.double(tstep),
                  dais_parameters    = as.double(parameters),
                  Ant_Temp           = as.double(Ta),
                  Ant_Sea_Level      = as.double(SL),
                  Ant_Sur_Ocean_Temp = as.double(Toc),
                  Ant_SL_rate        = as.double(dSL),
                  AIS_Radius_out     = as.double(rep(-999.99,ns)),
                  AIS_Volume_out     = as.double(rep(-999.99,ns))
  )
  Vsle = 57*(1-f.output$AIS_Volume_out/f.output$AIS_Volume_out[1]) #Takes steady state present day volume to correspond to 57m SLE

    return(Vsle)

}
# =================================================================================
