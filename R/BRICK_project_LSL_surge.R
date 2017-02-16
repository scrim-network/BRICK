# =======================================================================================
# Create projections of local sea level and storm surge.
# Use storm surge projections a la Grinsted et al 2013, with T=temperature
# anomaly relative to 1980-2000 average.
# And temperature projections from DOECLIM (within BRICK v0.1) under RCP2.6,
# 4.5 and 8.5.
# =======================================================================================
#
#   Requires (input variables):
# - lat.in      latitude of point at which you want local sea level (currently only support 1 pt)
# - lon.in      longitude of point at which you want local sea level
# - filename.brick  BRICK model output file to fingerprint LSL projections off of
# - tg.data     tide gauge data in the form [date | time | sl] where date is
#               YYYYMMDD, time is hh:mm, sl is meters, and these are the column
#               names (for referencing easily)
# - l.nonstat   which gev parameters to be non-stationary? these will covary
#               with temperature anomaly relative to 1980-2000 mean.
#               Expects list with names $location, $shape, $scale. Note that
#               Coles et al. 2001 (book: An Introduction to Statistical Modeling
#               of Extreme Values) suggest to leave the scale parameter stationary.
# - filename.brick  BRICK output file to base projections on
# - niter.gev   number of MCMC iterations to use to estimate (non-)stationary
#               GEV parameters
# - burnin      between 0 and 1, what fraction of the MCMC chain should be
#               discarded for burn-in? (default is 0.5)
#
#   Simulates (output variables):
# - lsl.out     local sea level rise (m sle) (Nens x Ntime)
# - gev.location    time series of location parameter (mu) (Nens x Ntime)
# - gev.shape       time series of shape parameter (xi) (Nens x Ntime)
# - gev.scale       time series of location parameter (sigma) (Nens x Ntime)
#
#TODO --
#TODO --
#TODO -- do the GEV parameters (non-stationary) for each of rcp26, 45 and 85.
#TODO -- --> actually, make general, with rcp as a list item for outputs
#TODO --        (like scenarios pipeline)
#TODO --
#TODO --

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

brick_project_lsl_surge <- function(
    lat.in = 0,
    lon.in = 0,
    filename.brick = NULL,
    tg.data = NULL,
    l.nonstat,
    niter.gev,
    burnin=0.5
    )
{

#===============================================================================
# Local sea level
#===============================================================================

  #install.packages('ncdf4')
  library(ncdf4)

  # Read sea level output
  ncdata <- nc_open(filename.brick)
    temp_rcp26     <- ncvar_get(ncdata, 'temp_RCP26')
    slr_gis_rcp26  <- ncvar_get(ncdata, 'GIS_RCP26')
    slr_gsic_rcp26 <- ncvar_get(ncdata, 'GSIC_RCP26')
    slr_ais_rcp26  <- ncvar_get(ncdata, 'AIS_RCP26')
    slr_te_rcp26   <- ncvar_get(ncdata, 'TE_RCP26')
    temp_rcp45     <- ncvar_get(ncdata, 'temp_RCP45')
    slr_gis_rcp45  <- ncvar_get(ncdata, 'GIS_RCP45')
    slr_gsic_rcp45 <- ncvar_get(ncdata, 'GSIC_RCP45')
    slr_ais_rcp45  <- ncvar_get(ncdata, 'AIS_RCP45')
    slr_te_rcp45   <- ncvar_get(ncdata, 'TE_RCP45')
    temp_rcp85     <- ncvar_get(ncdata, 'temp_RCP85')
    slr_gis_rcp85  <- ncvar_get(ncdata, 'GIS_RCP85')
    slr_gsic_rcp85 <- ncvar_get(ncdata, 'GSIC_RCP85')
    slr_ais_rcp85  <- ncvar_get(ncdata, 'AIS_RCP85')
    slr_te_rcp85   <- ncvar_get(ncdata, 'TE_RCP85')
    mod.time       <- ncvar_get(ncdata, 'time_proj')
    n.time         <- length(mod.time)
    n.ensemble     <- ncol(slr_gis_rcp26)
  nc_close(ncdata)

#TODO --
#TODO --
#TODO -- do all of these for each of rcp26, 45 and 85.
#TODO -- --> actually, make general, with rcp as a list item for outputs
#TODO --        (like scenarios pipeline)
#TODO --
#TODO --
  # testing for now only with RCP85
  slr_gis <- slr_gis_rcp85
  slr_gsic <- slr_gsic_rcp85
  slr_ais <- slr_ais_rcp85
  slr_te <- slr_te_rcp85

  # which of the two dimensions is time? if only one time series, turn into
  # a matrix for generality
  dims = dim(slr_gis)
  if(is.null(dims)) {
    slr_gis = as.matrix(slr_gis)
    slr_ais = as.matrix(slr_ais)
    slr_gsic = as.matrix(slr_gsic)
    slr_te = as.matrix(slr_te)
    itime = which(dim(slr_gis)==n.time)
  } else {
    itime = which(dims==n.time)
  }
  n.ensemble = dims[which(dims!=n.time)]
  if(is.null(n.ensemble)) {n.ensemble=1}

  # set up output with number of time points as number of rows, then re-shape the
  # output to match what was entered as input (slr_gis)
  lsl.out <- mat.or.vec(n.time,n.ensemble)

  # reshape the input to match this convention
  if(itime!=1) {
    slr_gis_reshape = t(slr_gis)
    slr_ais_reshape = t(slr_ais)
    slr_gsic_reshape = t(slr_gsic)
    slr_te_reshape = t(slr_te)
  } else {
    slr_gis_reshape = slr_gis
    slr_ais_reshape = slr_ais
    slr_gsic_reshape = slr_gsic
    slr_te_reshape = slr_te
  }

  # read the sea-level fingerprints in only one
  filename.fingerprints = "../fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"

  ncdata <- nc_open(filename.fingerprints)
    lat = ncvar_get(ncdata, 'lat')
    lon = ncvar_get(ncdata, 'lon')
    fp.gsic = ncvar_get(ncdata, 'GLAC')
    fp.gis = ncvar_get(ncdata, 'GIS')
    fp.ais = ncvar_get(ncdata, 'AIS')
  nc_close(ncdata)

  # convert longitude to degrees East, and find the fingerprinting data location
  # closest to the local sea level lat/lon given
  if(lon.in < 0) {lon.in=lon.in+360}	# convert longitude to degrees East
  ilat = which( abs(lat-lat.in)==min(abs(lat-lat.in)) )
  ilon = which( abs(lon-lon.in)==min(abs(lon-lon.in)) )

  # it is possible there were multiple lats/lons 'closest' to your given point
  # take the average of the non-NA of these
  fp.ais.loc = mean(fp.ais[ilon,ilat],na.rm=TRUE)
  fp.gsic.loc = mean(fp.gsic[ilon,ilat],na.rm=TRUE)
  fp.gis.loc = mean(fp.gis[ilon,ilat],na.rm=TRUE)
  fp.te.loc = 1.0		# TE response is to global mean temperature, so global mean sea level response is same everywhere

  # check if the nearest spot ended up on land.
  # If it did, take the average everywhere around the location.
  if(is.na(fp.ais.loc) | is.na(fp.gsic.loc) | is.na(fp.gis.loc) | is.na(fp.te.loc)) {
    fp.ais.loc = mean(fp.ais[(ilon-1):(ilon+1),(ilat-1):(ilat+1)], na.rm=TRUE)
    fp.gsic.loc = mean(fp.gsic[(ilon-1):(ilon+1),(ilat-1):(ilat+1)], na.rm=TRUE)
    fp.gis.loc = mean(fp.gis[(ilon-1):(ilon+1),(ilat-1):(ilat+1)], na.rm=TRUE)
    fp.te.loc = 1.0
  }

  # error message if something is still wrong
  if(is.na(fp.ais.loc) | is.na(fp.gsic.loc) | is.na(fp.gis.loc) | is.na(fp.te.loc)) {
    print('WARNING -- local sea level fingerprints are NaN')
  }

  # for each ensemble member, fingerprint the time series
  for (i in 1:n.ensemble) {
    lsl.out[,i] = fp.gis.loc * slr_gis_reshape[,i] +
                  fp.ais.loc * slr_ais_reshape[,i] +
                  fp.gsic.loc * slr_gsic_reshape[,i] +
                  fp.te.loc * slr_te_reshape[,i]
  }

  # reshape output to match the input
  if (itime!=1) {
    lsl.out = t(lsl.out)
  }

#===============================================================================
# Storm surge
#===============================================================================

library(extRemes)
library(fExtremes)
library(ismev)
library(lubridate)
library(zoo)

# Read temperature data for non-stationary GEV co-variate
# HADCRUT4 annual global mean surface temperature
# Note: all ensemble member files have same time series of ucnertainties, so just
# grabbing the first one.
data.temperature <- read.table("../data/HadCRUT.4.4.0.0.annual_ns_avg.txt")
obs.temp <- data.temperature[,2]
obs.temp.time <- data.temperature[,1]
obs.temp.err <- read.table("../data/HadCRUT.4.4.0.0.annual_ns_avg_realisations/HadCRUT.4.4.0.0.annual_ns_avg.1.txt")[,3]
ibeg <- which(obs.temp.time==years.unique[1])
iend <- which(obs.temp.time==max(years.unique))
obs.temp <- obs.temp[ibeg:iend]
obs.temp.time <- obs.temp.time[ibeg:iend]
obs.temp.err <- obs.temp.err[ibeg:iend]

# Normalize temperature relative to 1980-2000 mean, and create data 'obs'
ibeg <- which(obs.temp.time==1980)
iend <- which(obs.temp.time==2000)
obs.temp <- obs.temp - mean(obs.temp[ibeg:iend])
obs <- data.frame(cbind(years.unique,obs.temp))
colnames(obs) <- c('year','temp')

# Normalize the temperature projections too
ibeg <- which(mod.time==1980)
iend <- which(mod.time==2000)
temp_rcp26 <- temp_rcp26 - t(replicate(nrow(temp_rcp26),apply(temp_rcp26[ibeg:iend,],2,mean)))
temp_rcp45 <- temp_rcp45 - t(replicate(nrow(temp_rcp45),apply(temp_rcp45[ibeg:iend,],2,mean)))
temp_rcp85 <- temp_rcp85 - t(replicate(nrow(temp_rcp85),apply(temp_rcp85[ibeg:iend,],2,mean)))

# Get tide gauge data and prepare to analyze.
tg.data$date2 <- as.Date(as.character(tg.data$date), format="%Y%m%d")
years         <- floor(as.numeric(as.yearmon(tg.data$date2)))
years.unique  <- unique(years)
n.years       <- length(years.unique)
lsl.mean      <- rep(0,length(n.years))
lsl.max       <- rep(0,length(n.years))
tg.data$sl.norm <- rep(NA,length(years))

# 1. get unique years of data
# 2. subtract annual block means
# 3. get annual block maxima <-- these are what we fit GEV distribution to
for (tt in 1:n.years) {
    ind.thisyear <- which(years==years.unique[tt])
    lsl.mean[tt] <- mean(tg.data$sl[ind.thisyear])
    tg.data$sl.norm[ind.thisyear] <- tg.data$sl[ind.thisyear] - lsl.mean[tt]
    lsl.max[tt] <- max( tg.data$sl.norm[ind.thisyear] )
}

# cases: which GEV parameters are to covary with temperature?
if(!all(unlist(l.nonstat))) {
    # stationary GEV case -- just return a time series for each that is constant

#TODO

} else if(l.nonstat$location & !l.nonstat$shape & !l.nonstat$scale) {
    # non-stationary location only

#TODO

} else if(l.nonstat$location & !l.nonstat$shape & !l.nonstat$scale) {
    # non-stationary location and shape only

    # fit a preliminary maximum likelihood estimate
    gev.mle <- fevd(coredata(lsl.max), type='GEV',
                    data=obs, location.fun=~temp, shape.fun=~temp) # extRemes
    gev.mle2 <- gev.fit(coredata(lsl.max), show = FALSE, ydat=obs, mul=2, shl=2)   # ismev package

    # initial parameters estimates
    init[[1]] <- gev.mle$results$par[1]
    init[[2]] <- 0
    init[[3]] <- gev.mle$results$par[2]
    init[[4]] <- gev.mle$results$par[3]
    init[[5]] <- 0
    names(init) <- names(gev.mle$results$par)
    std.err <- gev.mle3$se
    names(std.err) <- names(gev.mle$results$par)

    # proposal sd estimates from a preliminary longer chain, with the default sd
    params.prop <- vector('list',1) # list order: mu0, m1, log.scale, xi0, xi1
    params.prop[[1]] <- as.numeric(1*std.err)
    params.prop[[1]][3] <- as.numeric(2*std.err[3])
    names(params.prop) <- 'sd'

    params.pri <- vector('list',1) # list order: mu0, m1, log.scale, xi0, xi1
    params.pri[[1]] <- as.numeric(4*std.err)
    names(params.pri) <- 'v'

    # Fit non-stationary GEV
    gev.bayes <- fevd(coredata(lsl.max), type="GEV", method="Bayesian",
                      data=obs, location.fun=~temp, shape.fun=~temp,
                      initial=init, iter=niter.gev, proposalParams=params.prop, priorParams=params.pri)

    gev.est <- gev.bayes$results[(round(burnin*niter.gev)+1):niter.gev, 1:length(std.err)]
    i.scale <- match('log.scale',names(gev.bayes$results[1,]))
    gev.est[,i.scale] <- exp(gev.est[,i.scale])   # account for log(scale) from MCMC
    colnames(gev.est) <- names(gev.bayes$results[1,1:ncol(gev.est)])

    # check acceptance rates - these should be close to 0.44 (for Metropolis/Gibbs
    # hybrid method, essentially single-variable)
    paccept <- apply(gev.bayes$chain.info[2:nrow(gev.bayes$chain.info),],2,sum)/nrow(gev.bayes$chain.info)
    paccept <- paccept[1:length(std.err)]
    print(paste('Acceptance rates should be around 0.44:',paccept))
    print('If they are not, go into BRICK_project_LSL_surge.R and tune the proposal distributions using params.prop')

    # draw parameter sets
    ind.ensemble <- sample( seq(1,nrow(gev.est)), size=n.ensemble, replace=FALSE)
    gev.sample <- gev.est[ind.ensemble,]
    location.rcp85 <- mat.or.vec(n.time,n.ensemble)
    shape.rcp85 <- mat.or.vec(n.time,n.ensemble)
    scale.rcp85 <- mat.or.vec(n.time,n.ensemble)
    for (i in 1:n.ensemble) {
        location.rcp85[,i] <- gev.sample[i,'mu0'] + gev.sample[i,'mu1']*temp_rcp85[,i]
        shape.rcp85[,i] <- gev.sample[i,'xi0'] + gev.sample[i,'xi1']*temp_rcp85[,i]
        scale.rcp85[,i] <- rep(exp(gev.sample[i,'log.scale']), n.time)
    }

} else {
    print('Sorry, that combination of non-stationary GEV parameters is not supported.')
#TODO -- add support for other combinations
}



#===============================================================================

    output <- list(lsl.out, location.rcp85, shape.rcp85, scale.rcp85)
    names(output) <- c('lsl','location','shape','scale')
    return(output)
}

#===============================================================================
# end
#===============================================================================
