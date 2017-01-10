# =======================================================================================
# Script to call GSIC_magicc (Fortran90):
# Simple model to simulate global contribution of Glaciers and Small Ice Caps (GSIC) to
#  sea-level rise (Wigley and Raper 2005, application of equation 4/5)
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

# =======================================================================================
# load fortran subroutine (# to check if library is loaded is.loaded("run_simple") )
# dyn.load("../fortran/gsic_magicc.so")
if(.Platform$OS.type == "unix") {
    dyn.load("../fortran/gsic_magicc.so")
} else {
    dyn.load("../fortran/gsic_magicc")
}

gsic_magiccF = function(beta0 = 0.000577,
                       V0    = 0.4,
                       n     = 0.82,
                       Gs0   = 0,
                       Teq   = -0.15,
                       tstep = 1,
                       i0    = 1,
                       Tg) {

  # determine series length
  ns <- length(Tg)

  # call fortran
  f.output <- .Fortran("run_gsic_magicc",
                  ns                = ns,
                  tstep             = as.double(tstep),
                  gsic_magicc_beta0 = as.double(beta0),
                  gsic_magicc_V0    = as.double(V0),
                  gsic_magicc_n     = as.double(n),
                  gsic_magicc_Gs0   = as.double(Gs0),
                  gsic_magicc_Teq   = as.double(Teq),
                  Gl_Temp           = as.double(Tg),
                  gsic_magicc_i0    = as.double(i0),
                  SL_contribution_out = as.double(rep(-999.99,ns))
  )
  return(f.output$SL_contribution_out)

}
# =================================================================================
