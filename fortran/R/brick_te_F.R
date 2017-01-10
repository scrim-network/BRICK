# =======================================================================================
# BRICK_TE fortran90 (# estimation by calling fortran routine)
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

# =======================================================================================
# load fortran subroutine (# to check if library is loaded is.loaded("run_brick_te") )
# dyn.load("../fortran/brick_te.so")
if(.Platform$OS.type == "unix") {
    dyn.load("../fortran/brick_te.so")
} else {
    dyn.load("../fortran/brick_te")
}

brick_te_F <- function(a = 0.5,
                       b = 0.0,
                       invtau = 0.005,
                       TE_0 = 0.0,
                       tstep = 1,
                       i0 = 1,
                       Tg) {

  # determine series length
  ns <- length(Tg)

  tau = 1/invtau

  # call fortran
  f.output <- .Fortran("run_brick_te",
                  ns            = ns,
                  tstep         = as.double(tstep),
                  brick_te_a    = as.double(a),
                  brick_te_b    = as.double(b),
                  brick_te_tau  = as.double(tau),
                  brick_te_TE_0 = as.double(TE_0),
                  Gl_Temp       = as.double(Tg),
                  brick_te_i0   = as.double(i0),
                  TE_out        = as.double(rep(-999.99,ns))
  )
  return(f.output$TE_out)

}
# =================================================================================
