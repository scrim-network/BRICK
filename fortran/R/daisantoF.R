# =======================================================================================
# daisanto: couples
# DAIS: Simple model for Antarctic ice-sheet volume [m sle] (Schaffer 2014)
# ANTO: Simple scaling to estimate Antarctic ocean-surface temperatures from atmospheric
#       Tg
# =======================================================================================
#
#  Requires (input variables):
#  - Tg        global mean surface temperature anomaly [deg C]
#  - SL        Sea level [m]
#  - slope.Ta2Tg       slope of linear regression between Antarctic temperature and global surface temperature anomalies
#  - intercept.Ta2Tg   intercept of this linear regression
#                       (if slope and intercept are not provided, assumed to be Antarctic temperature already)
#
#  Simulates (output variables):
#  - Rad       Ice sheet radius [m]
#  - Vais      Volume of Antarctic ice sheet [m3]
#  - Ta        Antarctic mean surface temperature [degC] (or global temperature)
#  - Toc       High latitude ocean subsurface temperatures [degC] (calculated by ANTO)

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
#  - anto.a    Sensitivity Ta [degC/degC]
#  - anto.b    Ta(Tg=0)
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
#  - dSL0      Initial sea level rate
#  - tstep     time step
#
#  Constants:
#  - Tf        Freecing temperature sea water [degC]
#  - rho_w      (Sea) water density [kg/m3]
#  - rho_i      Ice density [kg/m3]
#  - rho_m      Rock density [kg/m3]
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

source("../fortran/R/daisF.R")
source("../R/anto.R")

daisantoF <- function(
  tstep = 1,
  anto.a = 0.26,
  anto.b = 0.62,
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
  slope.Ta2Tg = 1,
  intercept.Ta2Tg = 0,
  Tg,        # global mean surface temperature anomaly [deg C]
  SL,
  dSL=0){

  Toc <- anto(a=anto.a, b=anto.b, Tf=Tf, Tg=Tg)

  # Reconstruct Antarctic surface temperature from global temperature anomalies.
  # Assumes Ta is global temperature anomaly. If it is not, then do not supply
  # slope.Ta2Tg and intercept.Ta2Tg to daisantoF(...) and it will do nothing.

  Ta.recon = (Tg-intercept.Ta2Tg)/slope.Ta2Tg   # undo linear regression Tg ~ Ta

  Contribution_AIS <- daisF(
    tstep  = tstep,
    b0     = b0,
    slope  = slope,
    mu     = mu,
    h0     = h0,
    c      =  c,
    P0     = P0,
    kappa  = kappa,
    nu     = nu,
    f0     = f0,
    gamma  = gamma,
    alpha  = alpha,
    Tf     = Tf,
    rho_w  = rho_w,
    rho_i  = rho_i,
    rho_m  = rho_m,
    Toc_0  = Toc_0,
    Rad0   = Rad0,
    Aoc    = Aoc,
    lf     = lf,
    includes_dSLais = includes_dSLais,
    Ta     = Ta.recon,
    SL     = SL,
    Toc    = Toc,
    dSL    = dSL)

return(Contribution_AIS)

}
