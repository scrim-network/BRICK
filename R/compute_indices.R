##==============================================================================
## Function to compute the indices within a vector of observational data and a
## vector of model output, that are of the same time points, for comparison.
##
## Input:
##  obs.time      times of observations
##  mod.time      times of model output
##
## Output:
##  oidx          indices within observational data vector that are same as
##  midx          these indices within the model output vector
##
## Original code is from Kelsey Ruckert (klr324@psu.edu).
## Modified for BRICK by Tony Wong (twong@psu.edu).
##
## Questions? Tony Wong (twong@psu.edu)
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

compute_indices <- function(
                              obs.time,
                              mod.time
                              ){

  oidx = NULL; midx = NULL
  endyear = mod.time[(length(mod.time))]
  if(min(obs.time) <= endyear) {
    oidx = which(obs.time <= endyear)
    midx = obs.time[oidx] - mod.time[1] + 1
    oidx = oidx[which(midx>0)] # only positive indices are within the simulation
    midx = midx[which(midx>0)]
  }
  return(list(oidx=oidx,midx=midx))
}

##==============================================================================
## End
##==============================================================================
