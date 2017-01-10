# =======================================================================================
#  Simple, mechanistically motivated model of the Greenland ice sheet volume [m sle] in
#  response to temperature variations (Bakker, Applegate and Keller 2016, equations 1-3)
#
# Modified to use first-order explicit integration forwards in time from the initial
# condition given, and first-order implicit integration backwards in time to the
# beginning of the climate simulation
# =======================================================================================
#
#  Requires (input variables):
#  - Tg        global temperature [degC]
#
#  Simulates (output variables):
#  - Vgrl      Volume of Greenland ice sheet [m sea-level equivalent (sle)]
#  - sle.gis   sea level rise (m) from Greenland ice sheet
#
#  Parameters:
#  - a         sensitivity of equilibrium volume Veq [m sle/degC]
#  - b         equilibrium volume Veq [m sle] for temperature Tg = 0
#  - alpha     sensitivity of exponential decay rate (1/tau)
#  - beta      exponential decay rate [1 / K] at Tg = 0
#  - V0        initial ice-sheet volume [m sle]
#  - tstep     time step
#  - i0        index of the "initial" year (V[i0]=V0, with V0 GIS volume left)
#
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

simple <- function(a      = -0.827,
                   b      =  7.242,
                   alpha  =  1.630e-4,
                   beta   =  2.845e-05,
                   V0     =  7.242,
                   tstep  = 1,
                   i0     = 1,
                   Tg)
{

  np      <- length(Tg)

  Vgrl    <- rep(NA, np)
  sle.gis <- rep(NA, np)
  Vgrl[i0] <- V0

  # Error check
  if (i0 < 1) {print('ERROR - i0 GIS < 1')}

  # first-order explicit time step from initial condition to end of simulation
  for(i in (i0+1):np){
    Veq     <- a * Tg[i-1] + b              # 'virtual' equilibrium volume (equation 2)
    tauinv  <- alpha * Tg[i-1] + beta       # 1/timescale (tauinv = 1/tau) (equation 3)

    Vgrl[i] <- max(0,                                                   #  (equations 1 and 4)
                   Vgrl[i-1]  +  tstep * ( ( Veq - Vgrl[i-1]) * tauinv)  )
  }

  # first order implicit "backward differentiation" from initial condition to
  # beginning of simulation
  if(i0>1) {
    for(i in seq(i0,2,by=-1)){
      Veq     <- a * Tg[i] + b              # 'virtual' equilibrium volume (equation 2)
      tauinv  <- alpha * Tg[i] + beta       # 1/timescale (tauinv = 1/tau) (equation 3)

      Vgrl[i-1] <- max(0,                                                   #  (equations 1 and 4)
                     Vgrl[i]  -  tstep * ( ( Veq - Vgrl[i]) * tauinv)  )
    }
  }

  # This is sea-level rise due to GIS *relative* to time given by i0, with
  # volume given by V0
  sle.gis <- V0 - Vgrl

  return(list(Vgrl=Vgrl,sle.gis=sle.gis))
}
