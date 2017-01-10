# =======================================================================================
# SIMPLE-fortran90 (# estimation by calling fortran routine)
# SIMPLE: Simple model for Greenland ice-sheet volume [m sle] (Bakker et al 2014)
#
# Modified to use first-order explicit integration forwards in time from the initial
# condition given, and first-order implicit integration backwards in time to the
# beginning of the climate simulation
# =======================================================================================
#
# #  Requires (input variables):
#  - Tg        global temperature [degC]
#
#  Simulates (output variables):
#  - Vgrl      Volume of Greenland ice sheet [m sea-level equivalent (sle)]
#
#  Internal variables:
#  - Veq       Equilibrium Volume [m sle]
#  - tau       Decay rate [1/ degC]
#
#  Parameters:
#  - a         sensitivity of equilibrium volume Veq [m sle/degC]
#  - b         equilibrium volume Veq [m sle] for temperature Tg = 0
#  - alpha     sensitivity of exponential decay rate (1/tau)
#  - beta      exponential decay rate [1 / degC] at Tg = 0
#  - Vgrl_0    initial ice-sheet volume [m sle]
#  - tstep     time step
#  - i0        index of the "initial" year (V[i0]=V0, with V0 GIS volume left)
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
# load fortran subroutine (# to check if library is loaded is.loaded("run_simple") )
# dyn.load("../fortran/simple.so")
if(.Platform$OS.type == "unix") {
    dyn.load("../fortran/simple.so")
} else {
    dyn.load("../fortran/simple")
}

simpleF <- function(a      = -0.827,
                    b      =  7.242,
                    alpha  =  1.630e-4,
                    beta   =  2.845e-05,
                    V0     =  7.242,
                    tstep  = 1,
                    i0     = 1,
                    Tg)
{

  # determine series length
  ns <- length(Tg)

  # call fortran
  f.output <- .Fortran("run_simple",
                  ns            = ns,
                  tstep         = as.double(tstep),
                  simple_a      = as.double(a),
                  simple_b      = as.double(b),
                  simple_alpha  = as.double(alpha),
                  simple_beta   = as.double(beta),
                  simple_V0     = as.double(V0),
                  Grl_Temp      = as.double(Tg),
                  simple_i0     = as.double(i0),
                  GIS_Volume_out = as.double(rep(-999.99,ns))
  )
  sle.gis <- V0 - f.output$GIS_Volume_out

  output <- list(Vgrl=f.output$GIS_Volume_out, sle.gis=sle.gis)

  return(output)

}
# =================================================================================
