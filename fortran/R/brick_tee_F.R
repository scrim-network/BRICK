# =======================================================================================
# BRICK_TEE fortran90 (# estimation by calling fortran routine)
# Explicitly estimate contribution of thermosteric expansion (TE) to global sea-
# level rise, given an ocean heat change.
# 
# Ben Vega-Westhoff, July 2017
# =======================================================================================
#
#  Requires (input variables):
#  - deltaH	change in ocean heat [J]
#
#  Simulates (output variables):
#  - TE        cumulative sea-level contribution from thermosteric expansion [m]
#              since t0 (i=1)
#              by definition Gs(1) = 0
#
#  Parameters:
#  - c	       heat capacity of conservative temperature [J/kg/K]
#  - a	       global ocean-averaged thermal expansion coefficient [kg/m3/K]
#  - rho       approximate density of global ocean [kg/m3]
#  - sa        global ocean surface area [m2]
#  - TE_0      initial sea level
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
    dyn.load("../fortran/brick_tee.so")
} else {
    dyn.load("../fortran/brick_tee")
}

brick_tee_F <- function(c = 3991.86795711963, #heat capacity of conservative temp (Graham and McDougall, 2013)
                       a = 0.16, #thermal exp. coeff, global estimate from Roquet et al., 2015
                       rho = 1027., #ocean-avg density
		       sa = 3.619e14, #ocean surface area, Eakins and Sharmin 2010
                       TE_0 = 0.0,
                       i0 = 1,
                       deltaH) {

  # determine series length
  ns <- length(deltaH) + 1
  
  # call fortran
  f.output <- .Fortran("run_brick_tee",
                  ns             = as.integer(ns),
                  brick_tee_c    = as.double(c),
                  brick_tee_a    = as.double(a),
                  brick_tee_rho  = as.double(rho),
                  ocsa 		 = as.double(sa),
                  brick_tee_TE_0 = as.double(TE_0),
                  deltaH         = as.double(deltaH),
		  brick_tee_i0   = as.double(i0),
                  TE_out         = as.double(rep(-999.99,ns))
  )
  return(f.output$TE_out)

}
# =================================================================================
