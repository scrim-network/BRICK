# =======================================================================================
#  Simple, mechanistically motivated model of global mean sea leave (GMSL) [m sle]
#  Based on Rahmstorf 2007, Equation 1
# =======================================================================================
#
#  Requires (input variables):
#  - Tg        global temperature [degC]
#
#  Simulates (output variables):
#  - gmsl      global mean sea level (m sle)]
#
#  Parameters:
#  - a         sensitivity of equilibrium sea level [m sle/year/degC]
#  - Teq       equilibrium temperature (at which there is no change in sea
#               level) [deg C]
#  - SL0       initial sea-level (take equal to 0 and normalize to match obs)
#               (consistent with normalizing the forcing temperatures to 0 at
#                 the beginning of the simulation, which is done in the driver)
#  - tstep     time step
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

gmsl_r07 <- function( a      = 0.0034,
                      Teq    = 0,
                      SL0    = 0,
                      tstep  = 1,
                      Tg)
{

  np      <- length(Tg)

  gmsl    <- rep(NA, np)
  gmsl[1] <- SL0

  for (i in 2:np) {

    gmsl[i] <- gmsl[i-1]  +  tstep * a * ( Tg[i-1]-Teq )

  }

  return(gmsl)
}
