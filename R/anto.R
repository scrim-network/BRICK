# =======================================================================================
# anto:
# simple scaling to estimate Antarctic ocean temperature (Toc) from global temperature (Tg)
# =======================================================================================
#
#  Requires (input variables):
#  - Tg        Global surface temperature [anomaly degC]
#
#  Simulates (output variables):
#  - Toc       Temperature of ocean surface at Antarctica [degC]
#
#  Parameters:
#  - a        Sensitivity Ta [degC/degC]
#  - b        Ta(Tg=0)
#  - Tf       Freezing temperature sea-water (lower bound)
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


anto <- function(a  =  0.26,
                  b  =  0.62,
                  Tf = -1.8,
                  Tg) {

  c  <- (Tf-b)/a
  Toc <-  Tf + (a*Tg + b-Tf) / (1 + exp(-Tg+c))

  return(Toc)
}
