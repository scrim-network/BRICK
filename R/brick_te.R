# =======================================================================================
# BRICK_TE:
# Simple model to simulate contribution of thermosteric expansion (TE) to global sea-
# level rise. (Grinsted, Moore and Jevrejeva 2010, equations 1/2 that were originally
# applied to global sea-level)
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
#  - TE        cumulative sea-level contribution from thermosteric expansion [m]
#              since t0 (i=1)
#              by definition Gs(1) = 0
#
#  Internal variables:
#  - TEeq      Equilibrium TE [m]
#
#  Parameters:
#  - a         sensitivity of TE equilibrium [m/K]
#  - b         equilibrium TE [m] for temperature anomaly Tg = 0
#  - invtau    1/time-scale of exponential decay (e-folding time) [years^-1]
#  - TE_0      initial sea level
#  - tstep     time step
#  - i0        index of the "initial" year (TE[i0]=TE0)
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

brick_te = function(a     = 0.5,
                    b     = 0.0,
                    invtau= 0.005,
                    TE_0  = 0.0,
                    tstep = 1,
                    i0    = 1,
                    Tg) {

  np    <- length(Tg)

  TE    <- rep(NA, np)
  TE[i0] <- TE_0

  # Error check
  if (i0 < 1) {print('ERROR - i0 TE < 1')}

  # first-order explicit time step from initial condition to end of simulation
  for(i in (i0+1):np){
    TEeq  <- a*Tg[i-1] + b
    TE[i] <- TE[i-1]  +  tstep * ( ( TEeq - TE[i-1]) * invtau) # sample 1/tau

  }

  # first order implicit "backward differentiation" from initial condition to
  # beginning of simulation
  if(i0>1) {
    for(i in seq(i0,2,by=-1)){
      TEeq  <- a*Tg[i] + b
      TE[i-1] <- TE[i] - tstep * ( ( TEeq - TE[i]) * invtau) # sample 1/tau
    }
  }

  return(TE)
}
