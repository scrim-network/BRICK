# =======================================================================================
# Create projections of local sea level and storm surge for Annapolis.
# Use storm surge projections a la Grinsted et al 2013, with T=temperature
# anomaly relative to 1980-2000 average.
# And temperature projections from DOECLIM (within BRICK v0.1) under RCP2.6,
# 4.5 and 8.5.
# =======================================================================================
#
#   Required settings (define below):
# - lat.proj    latitude of point at which you want local sea level (currently only support 1 pt)
# - lon.proj    longitude of point at which you want local sea level
# - dat.dir     directory with tide gauge data (m) in the form [date | time | sl]
#               where date is YYYYMMDD, time is hh:mm, sl is meters, and these
#               are the column names (for referencing easily)
# - l.nonstat   which gev parameters to be non-stationary? these will covary
#               with temperature anomaly relative to 1980-2000 mean.
#               excepts [location | shape | scale]. Note that Coles et al. 2001
#               (book: An Introduction to Statistical Modeling of Extreme Values)
#               suggest to leave the scale parameter stationary.
# - filename.brick  BRICK output file to base projections on
# - niter.gev   number of MCMC iterations to use to estimate (non-)stationary
#               GEV parameters
#
#   Simulates (output variables):
# - lsl.out     local sea level rise (m sle)
#
#   Parameters:
# - filename.fingerprints = "../fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"
#               sea level fingerprinting data (Slangen et al 2014)
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

# set things here to get different site projections

dat.dir <- '../data/tidegauge_Annapolis/'
filetype <- 'txt'
septype <- '\t'

lat.proj <- 38+(59/60)          # Annapolis tide gauge at 38° 59' N
lon.proj <- -(76+(28.9/60))     # 76° 28.9' W (NOAA Tides and Currents website)

l.nonstat <- vector('list',3); names(l.nonstat) <- c('location','shape','scale')
l.nonstat$location <- TRUE
l.nonstat$shape <- TRUE
l.nonstat$scale <- FALSE
niter.gev <- 1e4
burnin <- 0.5

filename.brick <- '../output_model/BRICK-fastdyn_physical_gamma_31Jan2017.nc'

# unless other data sets are quite different, below here likely does not need
# to be modified

##==============================================================================

## Read tide gauge data
files.tg <- list.files(path=dat.dir,pattern=filetype)

data <- read.table(paste(dat.dir,files.tg[1],sep=''), header = TRUE, sep=septype)
if(length(files.tg) > 1) {
    for (ff in 2:length(files.tg)) {
        data <- rbind(data, read.table(paste(dat.dir,files.tg[ff],sep=''), header = TRUE, sep=septype))
    }
}

##==============================================================================


##==============================================================================

## Make projections
source('../R/BRICK_project_LSL_surge.R')

proj.out <- brick_project_lsl_surge(lat.in = lat.proj,
                                    lon.in = lon.proj,
                                    filename.brick = filename.brick,
                                    tg.data = data,
                                    l.nonstat = l.nonstat,
                                    niter.gev = niter.gev,
                                    burnin = burnin
                                    )


# TODO -- some sort of output file?

##==============================================================================
## End
##==============================================================================
