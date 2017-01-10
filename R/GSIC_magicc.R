# =======================================================================================
# Simple model to simulate global contribution of Glaciers and Small Ice Caps (GSIC) to
# sea-level rise (Wigley and Raper 2005, application of equation 4/5)
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
#  - Gs        cumulative sea-level contribution since t0 (i=1)a [m]
#              by definition Gs(1) = 0
#
#  Parameters:
#  - beta0     initial mass balance sensitivity (how long it takes GSIC to respond to
#              increasing temps) [m/yr/C]
#  - V0        initial volume = max(Gs) [meter sle]
#  - n         exponent for area-volume scaling [-]
#  - Gs0       Gs[1]: the corrected corresponding sea-level rise in 1961 [m]
#  - Teq       equilibrium temperature (at which there is no change) [deg C]
#  - tstep     time step [years]
#  - i0        index of the "initial" year (Gs[i0]=Gs0, with V0 GSIC volume left)
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


gsic_magicc = function(beta0 = 0.000577,
                       V0    = 0.4,
                       n     = 0.82,
                       Gs0   = 0,
                       Teq   = -0.15,
                       tstep = 1,
                       i0    = 1,
                       Tg) {

  ns    <- length(Tg)

  Gs    <- rep(NA, ns)
  Gs[i0] <- Gs0

  # Error check
  if (i0 < 1) {print('ERROR - i0 GSIC < 1')}

  # first-order explicit time step from initial condition to end of simulation
  for(i in (i0+1):ns){
    Gs[i] <- Gs[i-1] + tstep * (beta0 * (Tg[i-1] - Teq) * (1-(Gs[i-1]/V0))^n)
  }

  # first order implicit "backward differentiation" from initial condition to
  # beginning of simulation
  if(i0>1) {
    for(i in seq(i0,2,by=-1)){
      Gs[i-1] <- Gs[i] - tstep * (beta0 * (Tg[i] - Teq) * (1-(Gs[i]/V0))^n)
    }
  }

  return(Gs)
}
