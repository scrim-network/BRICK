##==============================================================================
## brick_lws:
## Simple model to simulate contribution of land water storage changes to global
## sea level changes
##
## Modified to use first-order explicit integration forwards in time from the initial
## condition given, and first-order implicit integration backwards in time to the
## beginning of the climate simulation
##==============================================================================
##
##  Requires (input variables):
##  - Tg        global temperature [degC] (LWS does not depend on Tg, just used
##              for the time steps)
##
##  Simulates (output variables):
##  - LWS       cumulative sea-level contribution from land water storage
##              changes [m] since t0 (i=1); by definition LWS(1) = 0
##
##  Internal variables:
##  -
##
##  Parameters:
##  - lws0      initial sea level contribution from LWS
##    Using Dieng et al 2015 (doi:10.1088/1748-9326/10/12/124010):
##  - lws.mean  mean of lws contribution [m/year] (default 0.30 mm/y)
##  - lws.sd    stdev of lws contribution [m/year] (default 0.18 mm/y)
##  - tstep     time step
##  - i0        index of the "initial" year (LWS[i0]=LWS_0)
##
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

brick_lws = function(lws.mean = 0.00030,
                     lws.sd   = 0.00018,
                     lws0     = 0.0,
                     tstep    = 1,
                     i0       = 1,
                     Tg) {

  np <- length(Tg)
  LWS <- rep(NA, np)

  # Error check
  if (i0 < 1) {print('ERROR - i0 LWS < 1')}

  # first-order explicit time step from initial condition to end of simulation
  # well... that is basically what extrapolation of LWS trend is, with 1 year
  # time step.

  # scale lws.mean and lws.sd to represent time step (if not 1)
  lws.mean <- lws.mean*tstep
  lws.sd   <- lws.sd*tstep

  # cumulative sea level changes from land water storage (trends must be m/year)
  LWS <- cumsum(rnorm(n=np, mean=lws.mean, sd=lws.sd))

  # impose the initial condition
  LWS <- LWS + lws0 - LWS[i0]

  return(LWS)
}
