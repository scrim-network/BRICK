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
# - gev.out     ...?
# - output.file ...?
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

# Note: you should interactively run the tide gauge/storm surge projections.
# This should NOT be routinized. It is not a "one-size-fits-all" type of
# calculation.

dat.dir <- '../data/tidegauge_Annapolis/'
filetype <- 'txt'
septype <- '\t'

lat.proj <- 38+(59/60)          # Annapolis tide gauge at 38° 59' N
lon.proj <- -(76+(28.9/60))     # 76° 28.9' W (NOAA Tides and Currents website)
scen.rcp <- c('rcp26','rcp45','rcp85')

l.nonstat <- vector('list',3); names(l.nonstat) <- c('location','shape','scale')
l.nonstat$location <- FALSE
l.nonstat$shape <- FALSE
l.nonstat$scale <- FALSE
niter.gev <- 1e4
burnin <- 0.5

filename.brick <- '../output_model/BRICK-fastdyn_physical_gamma_31Jan2017.nc'

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.projout <- paste('../output_model/BRICK_project-lsl-surge_Annapolis_',today,'.nc',sep='')
filename.lslout  <- paste('../output_model/BRICK_project-lsl_Annapolis_',today,'.csv', sep="")

# unless other data sets are quite different, below here likely does not need
# to be modified
##==============================================================================


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
## Make local sea level projections

#install.packages('ncdf4')
library(ncdf4)

# Initialize lists for projections
temp <- vector('list', length(scen.rcp)); names(temp) <- scen.rcp
slr_gis <- temp
slr_ais <- temp
slr_gsic <- temp
slr_te <- temp
gmsl <- temp
lsl_proj <- temp

# Read sea level output
ncdata <- nc_open(filename.brick)
  for (rcp in scen.rcp) {
    if(rcp=='rcp26') {
        temp[[rcp]]     <- ncvar_get(ncdata, 'temp_RCP26')
        slr_gis[[rcp]]  <- ncvar_get(ncdata, 'GIS_RCP26')
        slr_ais[[rcp]]  <- ncvar_get(ncdata, 'AIS_RCP26')
        slr_gsic[[rcp]] <- ncvar_get(ncdata, 'GSIC_RCP26')
        slr_te[[rcp]]   <- ncvar_get(ncdata, 'TE_RCP26')
        gmsl[[rcp]]     <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
    } else if(rcp=='rcp45') {
        temp[[rcp]]     <- ncvar_get(ncdata, 'temp_RCP45')
        slr_gis[[rcp]]  <- ncvar_get(ncdata, 'GIS_RCP45')
        slr_ais[[rcp]]  <- ncvar_get(ncdata, 'AIS_RCP45')
        slr_gsic[[rcp]] <- ncvar_get(ncdata, 'GSIC_RCP45')
        slr_te[[rcp]]   <- ncvar_get(ncdata, 'TE_RCP45')
        gmsl[[rcp]]     <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
    } else if(rcp=='rcp60') {
        temp[[rcp]]     <- ncvar_get(ncdata, 'temp_RCP60')
        slr_gis[[rcp]]  <- ncvar_get(ncdata, 'GIS_RCP60')
        slr_ais[[rcp]]  <- ncvar_get(ncdata, 'AIS_RCP60')
        slr_gsic[[rcp]] <- ncvar_get(ncdata, 'GSIC_RCP60')
        slr_te[[rcp]]   <- ncvar_get(ncdata, 'TE_RCP60')
        gmsl[[rcp]]     <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP60')
    } else if(rcp=='rcp85') {
        temp[[rcp]]     <- ncvar_get(ncdata, 'temp_RCP85')
        slr_gis[[rcp]]  <- ncvar_get(ncdata, 'GIS_RCP85')
        slr_ais[[rcp]]  <- ncvar_get(ncdata, 'AIS_RCP85')
        slr_gsic[[rcp]] <- ncvar_get(ncdata, 'GSIC_RCP85')
        slr_te[[rcp]]   <- ncvar_get(ncdata, 'TE_RCP85')
        gmsl[[rcp]]     <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
    } else {
        print(paste('Error - unrecognized RCP scenario: ',rcp,sep=''))
    }
  }
  t.proj     <- ncvar_get(ncdata, 'time_proj')
  n.ensemble <- ncol(temp[[1]])
nc_close(ncdata)

# Routine for projecting local sea level
source('../R/BRICK_LSL.R')

# Get LSL projections for each RPC scenario
for (rcp in scen.rcp) {
    lsl_proj[[rcp]] = brick_lsl(lat.in   = lat.proj,
                                lon.in   = lon.proj,
                                n.time   = length(t.proj),
                                slr_gis  = slr_gis[[rcp]],
                                slr_gsic = slr_gsic[[rcp]],
						        slr_ais  = slr_ais[[rcp]],
						        slr_te   = slr_te[[rcp]])
}
##==============================================================================


##==============================================================================
## Write CSV with monthly projections

# normalize realtive to 1983-2001
inorm <- which(t.proj==1983):which(t.proj==2001)
lsl26 <- lsl_proj$rcp26 - t(replicate(nrow(lsl_proj$rcp26),apply(lsl_proj$rcp26[inorm,],2,mean)))
lsl45 <- lsl_proj$rcp45 - t(replicate(nrow(lsl_proj$rcp45),apply(lsl_proj$rcp45[inorm,],2,mean)))
lsl85 <- lsl_proj$rcp85 - t(replicate(nrow(lsl_proj$rcp85),apply(lsl_proj$rcp85[inorm,],2,mean)))

lsl.out <- cbind(t.proj,
                 apply(lsl26,1,median),
                 apply(lsl45,1,median),
                 apply(lsl85,1,median))
colnames(lsl.out) <- c('Year','RCP26','RCP45','RCP85')
write.table(lsl.out, file=filename.lslout, sep=",", qmethod="double", row.names=FALSE)

##==============================================================================


##==============================================================================
# get an ensemble of projections of GEV parameters

# routien for projecting storm surge (Grinsted et al 2013)
source('../R/BRICK_project_surge.R')

gev_proj <- brick_surge_grinsted(temperature = temp,
                                 time.proj   = t.proj,
                                 tidegauge   = data,
                                 l.nonstat   = l.nonstat,
                                 niter.gev   = niter.gev,
                                 burnin      = burnin)
##==============================================================================


##==============================================================================
## Write a netCDF ensemble output file including each of the RCP scenarios:
## (1) global total sea level, (2) local sea level. Also will want each
## contribution to global sea level rise, for the hindcast plots

dim.tproj <- ncdim_def('time_proj', 'years', as.double(t.proj))
dim.lat <- ncdim_def('lat', 'deg N', as.double(length(lat.proj)))
dim.lon <- ncdim_def('lon', 'deg E', as.double(length(lon.proj)))

dim.ensemble <- ncdim_def('ens', 'ensemble member', as.double(1:ncol(lsl_proj$rcp26)), unlim=TRUE)

lat.out <- ncvar_def('lat.lsl', 'deg N', list(dim.lat), -999,
                  longname = 'latitude of local sea level point')
lon.out <- ncvar_def('lon.lsl', 'deg N', list(dim.lon), -999,
                  longname = 'longitude of local sea level point')

lsl26 <- ncvar_def('LocalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level (RCP26)')
gsl26 <- ncvar_def('GlobalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level (RCP26)')
temp26 <- ncvar_def('temp_RCP26', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP26)')
loc26 <- ncvar_def('gev_location_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'location parameter (mu) of storm surge GEV distribution (RCP26)')
sha26 <- ncvar_def('gev_shape_RCP26', '', list(dim.tproj, dim.ensemble), -999,
                  longname = 'shape parameter (xi) of storm surge GEV distribution (RCP26)')
sca26 <- ncvar_def('gev_scale_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'scale parameter (sigma) of storm surge GEV distribution (RCP26)')

lsl45 <- ncvar_def('LocalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level (RCP45)')
gsl45 <- ncvar_def('GlobalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level (RCP45)')
temp45 <- ncvar_def('temp_RCP45', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP45)')
loc45 <- ncvar_def('gev_location_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'location parameter (mu) of storm surge GEV distribution (RCP45)')
sha45 <- ncvar_def('gev_shape_RCP45', '', list(dim.tproj, dim.ensemble), -999,
                  longname = 'shape parameter (xi) of storm surge GEV distribution (RCP45)')
sca45 <- ncvar_def('gev_scale_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'scale parameter (sigma) of storm surge GEV distribution (RCP45)')

lsl85 <- ncvar_def('LocalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level (RCP85)')
gsl85 <- ncvar_def('GlobalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level (RCP85)')
temp85 <- ncvar_def('temp_RCP85', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP85)')
loc85 <- ncvar_def('gev_location_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'location parameter (mu) of storm surge GEV distribution (RCP85)')
sha85 <- ncvar_def('gev_shape_RCP85', '', list(dim.tproj, dim.ensemble), -999,
                  longname = 'shape parameter (xi) of storm surge GEV distribution (RCP85)')
sca85 <- ncvar_def('gev_scale_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'scale parameter (sigma) of storm surge GEV distribution (RCP85)')

outnc <- nc_create(filename.projout,
                list( gsl26, gsl45, gsl85, lsl26, lsl45, lsl85, temp26, temp45, temp85,
                loc26, loc45, loc85, sha26, sha45, sha85, sca26, sca45, sca85, lat.out, lon.out),
                force_v4 = TRUE)

ncvar_put(outnc, lat.out, lat.proj)
ncvar_put(outnc, lon.out, lon.proj)

ncvar_put(outnc, lsl26, lsl_proj$rcp26)
ncvar_put(outnc, gsl26, gmsl$rcp26)
ncvar_put(outnc, temp26, temperature$rcp26)
ncvar_put(outnc, loc26, gev_proj$gev.proj$location$rcp26)
ncvar_put(outnc, sha26, gev_proj$gev.proj$shape$rcp26)
ncvar_put(outnc, sca26, gev_proj$gev.proj$scale$rcp26)

ncvar_put(outnc, lsl45, lsl_proj$rcp45)
ncvar_put(outnc, gsl45, gmsl$rcp45)
ncvar_put(outnc, temp45, temperature$rcp45)
ncvar_put(outnc, loc45, gev_proj$gev.proj$location$rcp45)
ncvar_put(outnc, sha45, gev_proj$gev.proj$shape$rcp45)
ncvar_put(outnc, sca45, gev_proj$gev.proj$scale$rcp45)

ncvar_put(outnc, lsl85, lsl_proj$rcp85)
ncvar_put(outnc, gsl85, gmsl$rcp85)
ncvar_put(outnc, temp85, temperature$rcp85)
ncvar_put(outnc, loc85, gev_proj$gev.proj$location$rcp85)
ncvar_put(outnc, sha85, gev_proj$gev.proj$shape$rcp85)
ncvar_put(outnc, sca85, gev_proj$gev.proj$scale$rcp85)

nc_close(outnc)
##==============================================================================


##==============================================================================
## Some plots

## Note -- later (if desired) can remove this part and put into a separate routine

#install.packages('ncdf4')
library(ncdf4)

filename.projin <- "../output_model/BRICK_project-lsl-surge_Annapolis_20Feb2017.nc"

# Initialize lists for projections
scen.rcp <- c('rcp26','rcp45','rcp85')
temp.proj <- vector('list', length(scen.rcp)); names(temp.proj) <- scen.rcp
lsl.proj <- temp.proj
gmsl.proj <- temp.proj
gev.proj <- vector('list',3); names(gev.proj) <- c('location','shape','scale')
gev.proj$location <- temp.proj
gev.proj$shape <- temp.proj
gev.proj$scale <- temp.proj

# Read sea level output
ncdata <- nc_open(filename.projin)
  for (rcp in scen.rcp) {
    if(rcp=='rcp26') {
        temp.proj[[rcp]]         <- ncvar_get(ncdata, 'temp_RCP26')
        gmsl.proj[[rcp]]         <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
        lsl.proj[[rcp]]          <- ncvar_get(ncdata, 'LocalSeaLevel_RCP26')
        gev.proj$location[[rcp]] <- ncvar_get(ncdata, 'gev_location_RCP26')
        gev.proj$shape[[rcp]]    <- ncvar_get(ncdata, 'gev_shape_RCP26')
        gev.proj$scale[[rcp]]    <- ncvar_get(ncdata, 'gev_scale_RCP26')
    } else if(rcp=='rcp45') {
        temp.proj[[rcp]]         <- ncvar_get(ncdata, 'temp_RCP45')
        gmsl.proj[[rcp]]         <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
        lsl.proj[[rcp]]          <- ncvar_get(ncdata, 'LocalSeaLevel_RCP45')
        gev.proj$location[[rcp]] <- ncvar_get(ncdata, 'gev_location_RCP45')
        gev.proj$shape[[rcp]]    <- ncvar_get(ncdata, 'gev_shape_RCP45')
        gev.proj$scale[[rcp]]    <- ncvar_get(ncdata, 'gev_scale_RCP45')
    } else if(rcp=='rcp60') {
        temp.proj[[rcp]]         <- ncvar_get(ncdata, 'temp_RCP60')
        gmsl.proj[[rcp]]         <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP60')
        lsl.proj[[rcp]]          <- ncvar_get(ncdata, 'LocalSeaLevel_RCP60')
        gev.proj$location[[rcp]] <- ncvar_get(ncdata, 'gev_location_RCP60')
        gev.proj$shape[[rcp]]    <- ncvar_get(ncdata, 'gev_shape_RCP60')
        gev.proj$scale[[rcp]]    <- ncvar_get(ncdata, 'gev_scale_RCP60')
    } else if(rcp=='rcp85') {
        temp.proj[[rcp]]         <- ncvar_get(ncdata, 'temp_RCP85')
        gmsl.proj[[rcp]]         <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
        lsl.proj[[rcp]]          <- ncvar_get(ncdata, 'LocalSeaLevel_RCP85')
        gev.proj$location[[rcp]] <- ncvar_get(ncdata, 'gev_location_RCP85')
        gev.proj$shape[[rcp]]    <- ncvar_get(ncdata, 'gev_shape_RCP85')
        gev.proj$scale[[rcp]]    <- ncvar_get(ncdata, 'gev_scale_RCP85')
    } else {
        print(paste('Error - unrecognized RCP scenario: ',rcp,sep=''))
    }
  }
  t.proj     <- ncvar_get(ncdata, 'time_proj')
  n.ensemble <- ncdata$dim$ens$len
nc_close(ncdata)

## Get statistical libraries needed for GEV calculations
library(extRemes)
library(fExtremes)
library(ismev)
library(lubridate)
library(zoo)
library(Bolstad)        # integrate (normalize) via Simpson's rule

## Calculate pdfs and sfs for sea-level rise and storm surge
iproj <- which(t.proj==2100)
inorm <- which(t.proj==1985):which(t.proj==2005)

## Initialize some list arrays
init <- vector('list',length(scen.rcp)); names(init) <- scen.rcp
lsl.norm <- init
f.lsl <- init
f.surge <- init
cdf.lsl <- init
cdf.surge <- init
sf.lsl <- init
sf.surge <- init

## distribution of sea-level, subsidence, and storm surges
lsl.lower <- 0
lsl.upper <- 20
lsl.n <- 2^11
lsl.x <- seq(lsl.lower, lsl.upper, length.out=lsl.n)
lsl.dx <- median(diff(lsl.x))

for (rcp in scen.rcp) {
    lsl.norm[[rcp]] <- lsl.proj[[rcp]] - t(replicate(nrow(lsl.proj[[rcp]]),apply(lsl.proj[[rcp]][inorm,],2,mean)))
    f.lsl[[rcp]]    <- density(x=lsl.norm[[rcp]][iproj,], from=lsl.lower, to=lsl.upper, n=lsl.n, kernel='gaussian')
    f.lsl[[rcp]]    <- f.lsl[[rcp]]$y/sintegral(x=f.lsl[[rcp]]$x, fx=f.lsl[[rcp]]$y)$value
    cdf.lsl[[rcp]]  <- cumsum(f.lsl[[rcp]]*lsl.dx)
    cdf.lsl[[rcp]]  <- cdf.lsl[[rcp]] - cdf.lsl[[rcp]][1]     # normalize (begin at 0)
    cdf.lsl[[rcp]]  <- cdf.lsl[[rcp]]/cdf.lsl[[rcp]][lsl.n]   # normalize (end at 1)
    sf.lsl[[rcp]]   <- 1-cdf.lsl[[rcp]]

## TODO
## TODO -- herenow -- add distributions for storm surge
## TODO

  }
}


##=====================================
## Figure 1 -- pd and survival functions of sea-level and storm surge in


##=====================================


##==============================================================================


##==============================================================================
## End
##==============================================================================
