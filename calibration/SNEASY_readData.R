## SNEASY_readData.R
##
##==============================================================================
## Copyright 2017 Tony Wong, Alexander Bakker
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
#
# Nathan M. Urban (nurban@princeton.edu)
# Woodrow Wilson School, Princeton
#
# Modified for BRICK by Tony Wong (twong@psu.edu), 4 April 2017
# originally loadobs.R
#
# Loads observational data.
#
# Also defines index variables to align data and model output
# vectors.

# HADCRUT3 annual global mean surface temperature
# http://www.cru.uea.ac.uk/cru/data/temperature/hadcrut3gl.txt
dat = read.table("../data/sneasy_tmp/temp.txt",header=TRUE)
obs.temp = dat[,2]-mean(dat[1:20,2])
obs.temp.time = dat[,1]

# annual global ocean heat content (0-3000 m)
# Gouretski & Koltermann, "How much is the ocean really warming?",
# Geophysical Research Letters 34, L01610 (2007).
dat = read.table("../data/sneasy_tmp/gouretski_ocean_heat_3000m.txt",skip=1)
obs.ocheat = dat[,2]
obs.ocheat.time = dat[,1]

# Mauna Loa instrumental CO2
dat = read.table("../data/sneasy_tmp/co2instobs.txt",header=TRUE)
obs.co2inst = dat[,2]
obs.co2inst.time = dat[,1]

# Law Dome ice core CO2
dat = read.table("../data/sneasy_tmp/co2iceobs.txt")
obs.co2ice = dat[,2]
obs.co2ice.err = rep(4, length(obs.co2ice))
obs.co2ice.time = dat[,1]

# decadal ocean carbon fluxes, McNeil et al. (2003)
obs.ocflux = c(1.6, 2.0)
obs.ocflux.err = c(0.4, 0.4)
obs.ocflux.time = c(1985, 1995)

# MOC strength, Bryden et al. (2005)
dat = read.table("../data/sneasy_tmp/bryden_moc_strength.txt")
obs.moc.bryden = dat[,2]
obs.moc.bryden.err = rep(6, rep(length(obs.moc.bryden)))
obs.moc.bryden.time = dat[,1]
# MOC strength, Lumpkin & Speer (2007)
dat = read.table("../data/sneasy_tmp/lumpkin_speer_moc_strength.txt")
obs.moc.ls = dat[,2]
obs.moc.ls.err = 2.5
obs.moc.ls.time = dat[,1]
# MOC strength, Cunningham et al. (2007)
dat = read.table("../data/sneasy_tmp/cunningham_moc_strength.txt")
obs.moc.cunn = dat[,2]
obs.moc.cunn.err = 2
obs.moc.cunn.time = dat[,1]

obs.moc = c(obs.moc.bryden, obs.moc.ls, obs.moc.cunn)
obs.moc.err = c(obs.moc.bryden.err, obs.moc.ls.err, obs.moc.cunn.err)
obs.moc.time = c(obs.moc.bryden.time, obs.moc.ls.time, obs.moc.cunn.time)

# index variables to align model output vectors with observation vectors
compute.indices = function(obs, obs.time)
{
	oidx = NULL; midx = NULL
	if(min(obs.time) <= endyear) {
		oidx = which(obs.time <= endyear)
		midx = obs.time[oidx] - mod.time[1] + 1
	}
	return(list(oidx=oidx,midx=midx))
}

idx = compute.indices(obs.temp, obs.temp.time)
oidx.temp = idx$oidx; midx.temp = idx$midx

idx = compute.indices(obs.ocheat, obs.ocheat.time)
oidx.ocheat = idx$oidx; midx.ocheat = idx$midx

idx = compute.indices(obs.co2inst, obs.co2inst.time)
oidx.co2inst = idx$oidx; midx.co2inst = idx$midx
idx = compute.indices(obs.co2ice, obs.co2ice.time)
oidx.co2ice = idx$oidx; midx.co2ice = idx$midx

idx = compute.indices(obs.ocflux, obs.ocflux.time)
oidx.ocflux = idx$oidx; midx.ocflux = idx$midx

idx = compute.indices(obs.moc, obs.moc.time)
oidx.moc = idx$oidx; midx.moc = idx$midx

##==============================================================================
## End
##==============================================================================
