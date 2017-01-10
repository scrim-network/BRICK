##==============================================================================
## Script to grab the Greenland ice sheet mass balance data to calibrate SIMPLE
## model.
##
## Read in the Greenland ice sheet observations in SLE, observation errors,
## years, and historic and emission temps.
## Read in the Greenland Mass balance data from 1958 - 2009
## The Mass balance is estimated as MB = surface mass balance + Discharge (Discharge is negative)
## Original data is from Sasgen et al (2012)
##
## Questions? Tony Wong <twong@psu.edu>
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


dat <- read.csv('../data/Greenland_OBS_MAR_InSAR_updated.csv')
obs.gis.time <- dat[1:52,9]
obs.gis <- dat[1:52,12] # Annual Mass Balance in meters sea level equivalence
obs.gis.err <- dat[1,15] # The error is +/- 30 Gt

idx = compute_indices(obs.time=obs.gis.time, mod.time=mod.time)
oidx.gis = idx$oidx; midx.gis = idx$midx

##==============================================================================
## End
##==============================================================================
