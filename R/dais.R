# =======================================================================================
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
#  - P         Precipitation [m/y (ice equivalent)]
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
#  - P0        P(Ta=0): Annual precipitation for Ta = 0 [m/y (ice equivalent)]
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
#  - lf         Mean AIS fingerprint at AIS shore
#
#  - tstep     time step
#
#  Constants:
#  - Tf        Freecing temperature sea water [degC]
#  - rho_w      (Sea) water density [kg/m3]
#  - rho_i      Ice density [kg/m3]
#  - rho_m      Rock density [kg/m3]
#  - Aoc        Surface of the ocean [m2]
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

dais <- function(
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
  Toc,        # High latitude ocean subsurface temperatures [degC]
  dSL) {

  # Initialize intermediate parameters
  del  <- rho_w/rho_i                # Ratio sea water and ice density [-]
  eps1 <- rho_i/(rho_m - rho_i)       # Ratio ice density and density difference between rock and ice [-]
  eps2 <- rho_w/(rho_m - rho_i)       # Ratio sea water density and density difference between rock and ice [-]

  # Define vectors with state variables
  np     <- length(Ta)
  Rad    <- rep(NA,np)               # Radius of ice sheet
  Vais   <- rep(NA,np)               # Ice volume

  # just some testing variables (to be discarded frmo the script later on)
  #isodepiso <- rep(NA,np)
    #dISO <- rep(NA,np)
  #dAIS <- rep(NA,np)


  # Initial conditions
  R       <- Rad0                  # gets updated at end of loop
  rc      <- (b0 - SL[1])/slope    # application of equation 1 (paragraph after eq3)
  V       <- pi * (1+eps1) * ( (8/15) *  mu^0.5 * R^2.5 - (1/3)*slope*R^3 ) # first estimate volume
  if(R>rc) { # in case of marine ice sheet correction term
    V <- V - pi*eps2 * ( (2/3) * slope*(R^3-rc^3)-b0*(R^2-rc^2) )
  }

  Vais[1] <- V
  Rad[1]  <- R
  #isodepiso[1] <- 0
  #dISO[1] =0
  #dAIS[1]=0

  # Run model
  for(i in 2:np) {

    hr   <- h0 + c * Ta[i-1]        # equation 5
    rc   <- (b0 - SL[i-1])/slope    # application of equation 1 (paragraph after eq3)
    P    <- P0 * exp(kappa*Ta[i-1]) # equation 6
    beta <- nu * P^(0.5)            # equation 7 (corrected with respect to text)

    # Total mass accumulation on ice sheet (equation 8)
    if(hr > 0) {
      rR   <- R - ((hr - b0 + slope*R)^2) / mu

      Btot <- P * pi * R^2 -
        pi * beta * (hr - b0 + slope*R) * (R*R - rR*rR) -
        (4 * pi * beta * mu^0.5 *   (R-rR)^2.5) / 5  +
        (4 * pi * beta * mu^0.5 * R*(R-rR)^1.5) / 3
    } else {
      Btot <- P * pi*R^2
    }


    # In case there is no marine ice sheet / grounding line
    F   <- 0                                                    # no ice flux
    ISO <- 0                                                    # (third term equation 14) NAME?
    fac <- pi * (1+eps1) * (4/3 * mu^0.5 * R^1.5 - slope*R^2)   # ratio dV/dR

    # In case there is a marine ice sheet / grounding line
    if (R > rc) {

      fac <- fac - 2*pi*eps2 * (slope*R^2 - b0*R)     # correction fac (=ratio dV/dR)

      Hw    <- slope*R - b0 + SL[i-1]                           # (equation 10)
      Speed <- f0 *                                             # Ice speed at grounding line (eq 11)
        ((1-alpha) + alpha* ((Toc[i-1] - Tf)/(Toc_0 - Tf))^2) *
        (Hw^gamma) / ( (slope*Rad0 - b0)^(gamma-1) )
      F     <- 2*pi*R * del * Hw * Speed                        # Ice flux (eq 9)

      # ISO term depends on dSL_tot (third term equation 14 !! NAME)
      c_iso <- 2*pi*eps2* (slope*rc^2 - (b0/slope)*rc)          # ratio ISO / dSL (all components)

      # first term is zero if dSL represents only non-AIS components (includes_dSLais=0)
      # second term is zero if dSL represents all components (includes_dSLais=1)
      ISO   <-    includes_dSLais  * c_iso             *  dSL[i] +                         #dSL = dSL_tot
               (1-includes_dSLais) * ((1-c_iso)/c_iso) * (dSL[i] - lf * (Btot - F) / Aoc)  #dSL = dSL_nonAIS
    }

    # Ice sheet volume (equation 13)
    R       <- R + tstep*(Btot-F+ISO)/fac
    V       <- V + tstep*(Btot-F+ISO)

    Rad[i]  <- R
    Vais[i] <- V

    # just some testing features (to be discarded later)
    #isodepiso[i] <- isodepiso[i-1] + ISO
    #dISO[i] <- ISO
    #dAIS[i] <- (Btot-F+ISO)
  }

  Vsle = 57*(1-Vais/Vais[1]) #Takes steady state present day volume to correspond to 57m SLE

  # just some testing features (to be discarded later)
  #return(data.frame(Vsle=Vsle, ISO=isodepiso/Aoc, dISO=dISO/Aoc, dAIS=dAIS/Aoc, dSL=dSL))
  return(Vsle)
}
