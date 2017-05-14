##==============================================================================
##  processingPipeline_VanDantzig_allScenarios.R
##
## Script to read and process the physical model output for all three fast
## dyanmics scenarios (and the 3 RCP scenarios within each), and both stationary
## and non-stationary storm surge. This is so you do not need to re-run the full
## processingPipeline_BRICKscenarios.R script 3 times. Instead, once you have the
## BRICK_physical_fd-XXX_[datestamp].nc files, you can feed them in here
## and just run this one.
##
## Why, you may ask? Well, the Van Dantzig takes a while. It is not coded in
## Fortran because R has the nice statistical utilities (calculating the CDF of
## the GEV functions), so we leave it in R. Plus, transparency is furthered in a
## more gentle programming language - nice for decision analysis.
##
## Questions? Tony Wong (twong@psu.edu)
##
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

rm(list=ls())

##==============================================================================
## Input file names from previous parameter sampling/calibration/simulation
## If you are running this for the first time (after
## 'processingPipeline_BRICKscenarios.R', )

## Sea-level rise projections
#filename.allslr <- "../output_model/BRICK_physical_allslr_16Apr2017.nc"
filename.gamma <- "../output_model/BRICK_physical_fd-gamma_08May2017.nc"
filename.uniform <- "../output_model/BRICK_physical_fd-uniform_08May2017.nc"

## GEV parameters, fit from tide gauge data
#filename.gevstat <- '../output_calibration/BRICK_GEVsample-AnnMean_16Apr2017.nc'
filename.gevmcmc <- '../output_calibration/BRICK_estimateGEV-AnnMean_12Apr2017.nc'

## Surge level increase factors (USACE)
#filename.surgefactor <- '../output_calibration/BRICK_surgefactor_16Apr2017.nc'

## Van Dantzig model parameters
#filename.vdparams <- '../output_calibration/BRICK-VanDantzig_parameters_16Apr2017.nc'

##==============================================================================
## What time horizon to consider?
currentyear <- 2015
endyear <- 2065

##==============================================================================
library(ncdf4)

today=Sys.Date(); today=format(today,format="%d%b%Y")

# no fast dynamics case uses the same projections as gamma/uniform, but with the
# fast dynamic AIS disintegration contribution neglected (subtracted)
scen.ais <- c("none","gamma","uniform")
scen.rcp <- c("rcp26","rcp45","rcp60","rcp85")
slr <- vector('list',length(scen.rcp)); names(slr) <- scen.rcp
slr.nofd <- vector('list',length(scen.rcp)); names(slr) <- scen.rcp
for (rcp in scen.rcp) {
  slr[[rcp]] <- vector('list',length(scen.ais)); names(slr[[rcp]]) <- scen.ais
  slr.nofd[[rcp]] <- vector('list',length(scen.ais)); names(slr.nofd[[rcp]]) <- scen.ais
}

##==============================================================================
if(exists('filename.allslr')) {
  ncdata <- nc_open(filename.allslr)
    slr$rcp26$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP26')
    slr$rcp45$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP45')
    slr$rcp60$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP60')
    slr$rcp85$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP85')
    slr$rcp26$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP26')
    slr$rcp45$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP45')
    slr$rcp60$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP60')
    slr$rcp85$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP85')
    slr$rcp26$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP26')
    slr$rcp45$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP45')
    slr$rcp60$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP60')
    slr$rcp85$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP85')
    mod.time <- ncvar_get(ncdata, 'time_proj')
    n.ens <- ncol(slr$rcp26$none)
  nc_close(ncdata)
} else {
  # get all ensembles, clip to make sure they are all the same size, then write
  # the 'allslr' file
  ncdata <- nc_open(filename.gamma)
    n.gamma <- ncol(ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP26'))
    slr$rcp26$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_RCP26')
    slr$rcp45$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_RCP45')
    slr$rcp60$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_RCP60')
    slr$rcp85$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_RCP85')
    slr.nofd$rcp26$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP26')
    slr.nofd$rcp45$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP45')
    slr.nofd$rcp60$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP60')
    slr.nofd$rcp85$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP85')
    mod.time <- ncvar_get(ncdata, 'time_proj')
  nc_close(ncdata)

  ncdata <- nc_open(filename.uniform)
    n.uniform <- ncol(ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP26'))
    slr$rcp26$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_RCP26')
    slr$rcp45$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_RCP45')
    slr$rcp60$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_RCP60')
    slr$rcp85$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_RCP85')
    slr.nofd$rcp26$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP26')
    slr.nofd$rcp45$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP45')
    slr.nofd$rcp60$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP60')
    slr.nofd$rcp85$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP85')
  nc_close(ncdata)

  n.ens <- min(n.uniform, n.gamma)
  for (rcp in scen.rcp) {
    for (ais in c('gamma','uniform')) {
      slr[[rcp]][[ais]] <- slr[[rcp]][[ais]][,1:n.ens]
      slr.nofd[[rcp]][[ais]] <- slr.nofd[[rcp]][[ais]][,1:n.ens]
    }
    slr[[rcp]]$none <- slr.nofd[[rcp]]$gamma
    irnd <- sample(1:n.ens, floor(0.5*n.ens), replace=FALSE)
    slr[[rcp]]$none[,irnd] <- slr.nofd[[rcp]]$uniform[,irnd]
  }

  # write the ensembles to a file so you don't have ot do this again
  today=Sys.Date(); today=format(today,format="%d%b%Y")
  filename.slrout = paste('../output_model/BRICK_physical_allslr_',today,'.nc',sep="")
  dim.time <- ncdim_def('n.time', '', 1:nrow(slr$rcp26$none), unlim=FALSE)
  dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:ncol(slr$rcp26$none), unlim=TRUE)
  slr.rcp26.none <- ncvar_def('LocalSeaLevel_nofd_RCP26', '', list(dim.time,dim.ensemble), -999)
  slr.rcp45.none <- ncvar_def('LocalSeaLevel_nofd_RCP45', '', list(dim.time,dim.ensemble), -999)
  slr.rcp60.none <- ncvar_def('LocalSeaLevel_nofd_RCP60', '', list(dim.time,dim.ensemble), -999)
  slr.rcp85.none <- ncvar_def('LocalSeaLevel_nofd_RCP85', '', list(dim.time,dim.ensemble), -999)
  slr.rcp26.gamma <- ncvar_def('LocalSeaLevel_gamma_RCP26', '', list(dim.time,dim.ensemble), -999)
  slr.rcp45.gamma <- ncvar_def('LocalSeaLevel_gamma_RCP45', '', list(dim.time,dim.ensemble), -999)
  slr.rcp60.gamma <- ncvar_def('LocalSeaLevel_gamma_RCP60', '', list(dim.time,dim.ensemble), -999)
  slr.rcp85.gamma <- ncvar_def('LocalSeaLevel_gamma_RCP85', '', list(dim.time,dim.ensemble), -999)
  slr.rcp26.gamma.nofd <- ncvar_def('LocalSeaLevel_gamma_nofd_RCP26', '', list(dim.time,dim.ensemble), -999)
  slr.rcp45.gamma.nofd <- ncvar_def('LocalSeaLevel_gamma_nofd_RCP45', '', list(dim.time,dim.ensemble), -999)
  slr.rcp60.gamma.nofd <- ncvar_def('LocalSeaLevel_gamma_nofd_RCP60', '', list(dim.time,dim.ensemble), -999)
  slr.rcp85.gamma.nofd <- ncvar_def('LocalSeaLevel_gamma_nofd_RCP85', '', list(dim.time,dim.ensemble), -999)
  slr.rcp26.uniform <- ncvar_def('LocalSeaLevel_uniform_RCP26', '', list(dim.time,dim.ensemble), -999)
  slr.rcp45.uniform <- ncvar_def('LocalSeaLevel_uniform_RCP45', '', list(dim.time,dim.ensemble), -999)
  slr.rcp60.uniform <- ncvar_def('LocalSeaLevel_uniform_RCP60', '', list(dim.time,dim.ensemble), -999)
  slr.rcp85.uniform <- ncvar_def('LocalSeaLevel_uniform_RCP85', '', list(dim.time,dim.ensemble), -999)
  slr.rcp26.uniform.nofd <- ncvar_def('LocalSeaLevel_uniform_nofd_RCP26', '', list(dim.time,dim.ensemble), -999)
  slr.rcp45.uniform.nofd <- ncvar_def('LocalSeaLevel_uniform_nofd_RCP45', '', list(dim.time,dim.ensemble), -999)
  slr.rcp60.uniform.nofd <- ncvar_def('LocalSeaLevel_uniform_nofd_RCP60', '', list(dim.time,dim.ensemble), -999)
  slr.rcp85.uniform.nofd <- ncvar_def('LocalSeaLevel_uniform_nofd_RCP85', '', list(dim.time,dim.ensemble), -999)
  time.proj <- ncvar_def('time_proj', '', list(dim.time), -999)
  outnc <- nc_create(filename.slrout,
                     list(slr.rcp26.none, slr.rcp45.none, slr.rcp60.none, slr.rcp85.none,
                          slr.rcp26.gamma, slr.rcp45.gamma, slr.rcp60.gamma, slr.rcp85.gamma,
                          slr.rcp26.gamma.nofd, slr.rcp45.gamma.nofd, slr.rcp60.gamma.nofd, slr.rcp85.gamma.nofd,
                          slr.rcp26.uniform, slr.rcp45.uniform, slr.rcp60.uniform, slr.rcp85.uniform,
                          slr.rcp26.uniform.nofd, slr.rcp45.uniform.nofd, slr.rcp60.uniform.nofd, slr.rcp85.uniform.nofd,
                          time.proj))
  ncvar_put(outnc, slr.rcp26.none, slr$rcp26$none)
  ncvar_put(outnc, slr.rcp45.none, slr$rcp45$none)
  ncvar_put(outnc, slr.rcp60.none, slr$rcp60$none)
  ncvar_put(outnc, slr.rcp85.none, slr$rcp85$none)
  ncvar_put(outnc, slr.rcp26.gamma, slr$rcp26$gamma)
  ncvar_put(outnc, slr.rcp45.gamma, slr$rcp45$gamma)
  ncvar_put(outnc, slr.rcp60.gamma, slr$rcp60$gamma)
  ncvar_put(outnc, slr.rcp85.gamma, slr$rcp85$gamma)
  ncvar_put(outnc, slr.rcp26.gamma.nofd, slr.nofd$rcp26$gamma)
  ncvar_put(outnc, slr.rcp45.gamma.nofd, slr.nofd$rcp45$gamma)
  ncvar_put(outnc, slr.rcp60.gamma.nofd, slr.nofd$rcp60$gamma)
  ncvar_put(outnc, slr.rcp85.gamma.nofd, slr.nofd$rcp85$gamma)
  ncvar_put(outnc, slr.rcp26.uniform, slr$rcp26$uniform)
  ncvar_put(outnc, slr.rcp45.uniform, slr$rcp45$uniform)
  ncvar_put(outnc, slr.rcp60.uniform, slr$rcp60$uniform)
  ncvar_put(outnc, slr.rcp85.uniform, slr$rcp85$uniform)
  ncvar_put(outnc, slr.rcp26.uniform.nofd, slr.nofd$rcp26$uniform)
  ncvar_put(outnc, slr.rcp45.uniform.nofd, slr.nofd$rcp45$uniform)
  ncvar_put(outnc, slr.rcp60.uniform.nofd, slr.nofd$rcp60$uniform)
  ncvar_put(outnc, slr.rcp85.uniform.nofd, slr.nofd$rcp85$uniform)
  ncvar_put(outnc, time.proj, mod.time)
  nc_close(outnc)
}

##==============================================================================
## Draw/read (stationary) GEV parameters for each of the n.ens SOW
source('../calibration/stormsurge_NOLA.R')
library(extRemes)
library(zoo)
library(coda)

if(exists('filename.gevstat')){
	ncdata <- nc_open(filename.gevstat)
	gev.names <- ncvar_get(ncdata, 'GEV_names')
	gev.params <- t(ncvar_get(ncdata, 'GEV_parameters'))
	nc_close(ncdata)
	colnames(gev.params) <- gev.names
} else {
    if(exists('filename.gevmcmc')){
        # sample from the previous MCMC results
        ncdata <- nc_open(filename.gevmcmc)
    	gev.names <- ncvar_get(ncdata, 'GEV_names')
    	gev.mcmc <- t(ncvar_get(ncdata, 'GEV_parameters'))
    	nc_close(ncdata)
    	colnames(gev.mcmc) <- gev.names
        isamp <- sample(x=1:nrow(gev.mcmc), size=n.ens, replace=TRUE)
        gev.params <- gev.mcmc[isamp,]
        # and write the results to a file so you can reproduce this later
        lmax=0
        for (i in 1:ncol(gev.params)){lmax=max(lmax,nchar(colnames(gev.params)[i]))}
        today=Sys.Date(); today=format(today,format="%d%b%Y")
        filename.gevshort = paste('../output_calibration/BRICK_GEVsample-AnnMean_',today,'.nc',sep="")
        dim.parameters <- ncdim_def('n.parameters', '', 1:ncol(gev.params), unlim=FALSE)
        dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
        dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(gev.params), unlim=TRUE)
        parameters.var <- ncvar_def('GEV_parameters', '', list(dim.parameters,dim.ensemble), -999)
        parnames.var <- ncvar_def('GEV_names', '', list(dim.name,dim.parameters), prec='char')
        outnc <- nc_create(filename.gevshort, list(parameters.var,parnames.var))
        ncvar_put(outnc, parameters.var, t(gev.params))
        ncvar_put(outnc, parnames.var, colnames(gev.params))
        nc_close(outnc)
    } else {
        # sd.prop from a preliminary MLE experiment which found 90% CI for the GEV
        # parameters that gives stdev estimates of 13.5 (mu), .2 (log(sigma)), .25 (xi)
        # if interpreted as a 4-sigma window.
        # sd.prop = [22,.15,.23] is tuned to get about .44 acceptance rate (Rosenthal 2011)
        # Note that sd.prop and other prior function/parameters are set within
        # this script. The calibration is not a one-size-fits-all type deal; you
        # *will* need to make detailed modifications for different data or
        # different statistical models for fitting.
        gev.est <- BRICK_estimateGEV_NOLA(niter=5e5, burnin=0.5, N=n.ens)
        gev.params <- gev.est[[2]]
    }
}

##==============================================================================
## Draw/read surge.factor parameters for each of the n.ens SOW

## Ref: Table 1.2, p.1-22, of USACE "Hurricane and Storm Damage Risk Reduction
## System Design Guidelines" (with revisions through June 2012). Accessible
## here (as of 30 Nov 2016):
## http://www.mvn.usace.army.mil/Portals/56/docs/engineering/HurrGuide/EntireDocument.pdf

## set min/max surge factors
sf.min <- 1.5
sf.max <- 2

## sample the surge factors for each model run. min/max (or pick a different
## distribution to sample from) are set above
if(exists('filename.surgefactor')){
  ncdata <- nc_open(filename.surgefactor)
  surge.factor <- as.vector(ncvar_get(ncdata, 'surge_factor'))
  nc_close(ncdata)
} else {
  if(is.null(n.ens)) {n.ens <- 1}
  surge.factor <- runif(n=n.ens, min=sf.min, max=sf.max)
  today=Sys.Date(); today=format(today,format="%d%b%Y")
  filename.surgefactor = paste('../output_calibration/BRICK_surgefactor_',today,'.nc',sep="")
  dim.parameters <- ncdim_def('n.parameters', '', 1, unlim=FALSE)
  dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:n.ens, unlim=TRUE)
  parameters.var <- ncvar_def('surge_factor', '', list(dim.parameters,dim.ensemble), -999)
  outnc <- nc_create(filename.surgefactor, list(parameters.var))
  ncvar_put(outnc, parameters.var, surge.factor)
  nc_close(outnc)
}

##==============================================================================
## Draw/read Van Dantzig parameters for each of the n.ens SOW
source("../R/VD_NOLA.R")					# contains parameter sampling routine (called below)

if(exists('filename.vdparams')) {
  # read pre-existing parameters file for the Van Dantzig model
  ncdata <- nc_open(filename.vdparams)
  params.vd <- t(ncvar_get(ncdata, 'VanDantzig_parameters'))
  parnames.vd <- ncvar_get(ncdata, 'VanDantzig_parnames')
  nc_close(ncdata)
  colnames(params.vd) <- parnames.vd
  if(nrow(params.vd) < n.ens) {print('ERROR - number of Van Dantzig parameters is less than number of sea-level rise realizations')}
  if(nrow(params.vd) > n.ens){params.vd <- params.vd[sample(1:nrow(params.vd), size=n.ens, replace=FALSE) , ]}
  params.vd <- as.data.frame(params.vd)
} else {
  # if you do not have a Van Dantzig parameters file to read, sample and write
  # a new one.
  p_zero_p = 0.0038                 # Initial flood frequency (1/yr) with zero height increase (Van Dantzig (1956))
  alpha_p = 2.6                     # Exponential flood frequency constant (Van Dantzig (1956))
  V_p = c(7.5e+9, 3e+10)            # Value of goods protected by dike (based on estimates in Jonkman et al. (2009)) (US$)
  delta_prime_p = c(0.02, 0.06)     # Discount rate (percent/year) (based on estimates in Jonkman et al. (2009))
  I_range = c(-0.5,1.0)             # Investment cost uncertainty range (as fraction of investment cost) (Jonkman and Dutch Perspective use -50%, +100%)
  subs_rate = 0.0056                # Rate of land subsidence (meter/year) (Dixon et al. (2006))

  ## Sample parameters (using Dutch Perspective as a guide). Uses Latin Hypercube,
  ## assigns each ensemble member a set of parameters.
  params.vd = parameter_sampling_DP(ncol(slr$rcp26$none), p_zero_p, alpha_p, V_p, delta_prime_p, I_range, subs_rate)

  ## Write file with Van Dantzig parameters
  today=Sys.Date(); today=format(today,format="%d%b%Y")
  filename.vdparams <- paste('../output_calibration/BRICK-VanDantzig_parameters_',today,'.nc',sep='')

  ## Get maximum length of parameter name, for width of array to write to netcdf
  ## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
  ## to get back into the shape BRICK expects, just transpose it
  lmax=0
  for (i in 1:ncol(params.vd)){lmax=max(lmax,nchar(colnames(params.vd)[i]))}

  library(ncdf4)
  dim.vdparams <- ncdim_def('params', 'parameters', (1:ncol(params.vd)) , unlim=FALSE)
  dim.ensemble <- ncdim_def('ens', 'ensemble member', (1:nrow(params.vd)) , unlim=TRUE)
  dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
  parameters.var <- ncvar_def('VanDantzig_parameters', '', list(dim.vdparams,dim.ensemble), -999)
  parnames.var <- ncvar_def('VanDantzig_parnames', '', list(dim.name,dim.vdparams), prec='char')
  outnc <- nc_create(filename.vdparams, list(parameters.var,parnames.var))
  ncvar_put(outnc, parameters.var, t(params.vd))
  ncvar_put(outnc, parnames.var, colnames(params.vd))
  nc_close(outnc)
}

##==============================================================================
## For each of the n.ens SOW, we now have:
## params.vd[sow,], surge.factor[sow], gev.params[sow,], slr[[rcp]][[ais]][,sow]
##==============================================================================

## Comment/modify in case you only want to (re-)run a selection of the scenarios
scen.ais <- c("none","gamma","uniform")
#scen.ais <- c("uniform","gamma")
#scen.ais <- 'none'

for (ais in scen.ais) {

  ## Set up the file names for output
  filename.vdout = paste("../output_model/VanDantzig_fd-",ais,"_",endyear,"_",today,".nc",sep="")

  ## Trim the sea-level rise projections to currentyear:endyear
  sea_level_rcp26 <- slr$rcp26[[ais]][which(mod.time==currentyear):which(mod.time==endyear),]
  sea_level_rcp45 <- slr$rcp45[[ais]][which(mod.time==currentyear):which(mod.time==endyear),]
  sea_level_rcp60 <- slr$rcp60[[ais]][which(mod.time==currentyear):which(mod.time==endyear),]
  sea_level_rcp85 <- slr$rcp85[[ais]][which(mod.time==currentyear):which(mod.time==endyear),]
  time.proj <- currentyear:endyear

  ## Yields ss$year, ss$surge.rise
  ss.rcp26 <- BRICK_estimateStormSurge_NOLA_usace(time.proj=time.proj, sea_level=sea_level_rcp26, surge.factor=surge.factor)
  ss.rcp45 <- BRICK_estimateStormSurge_NOLA_usace(time.proj=time.proj, sea_level=sea_level_rcp45, surge.factor=surge.factor)
  ss.rcp60 <- BRICK_estimateStormSurge_NOLA_usace(time.proj=time.proj, sea_level=sea_level_rcp60, surge.factor=surge.factor)
  ss.rcp85 <- BRICK_estimateStormSurge_NOLA_usace(time.proj=time.proj, sea_level=sea_level_rcp85, surge.factor=surge.factor)

  ## Evaluate flood risk analysis model for each of these realizations

  ## Set up the Van Dantzig simulations (only need to do this once, so separate)
  #source('../R/BRICK_VanDantzig.R')
  source('../R/BRICK_VanDantzig_setup.R')

  ## Set up the model for the actual simulations
  source('../R/BRICK_VanDantzig_run.R')

  # stationary storm surge GEV
  vandantzig.ensemble.rcp26.ssst <- brick_vandantzig(currentyear, endyear, sea_level=sea_level_rcp26, time=time.proj, ss.gev=gev.params, params.vd=params.vd)
  vandantzig.ensemble.rcp45.ssst <- brick_vandantzig(currentyear, endyear, sea_level=sea_level_rcp45, time=time.proj, ss.gev=gev.params, params.vd=params.vd)
  vandantzig.ensemble.rcp60.ssst <- brick_vandantzig(currentyear, endyear, sea_level=sea_level_rcp60, time=time.proj, ss.gev=gev.params, params.vd=params.vd)
  vandantzig.ensemble.rcp85.ssst <- brick_vandantzig(currentyear, endyear, sea_level=sea_level_rcp85, time=time.proj, ss.gev=gev.params, params.vd=params.vd)

  # non-stationary storm surge GEV
  vandantzig.ensemble.rcp26.ssns <- brick_vandantzig(currentyear, endyear, sea_level=sea_level_rcp26, time=time.proj, ss.gev=gev.params, surge.rise=ss.rcp26$surge.rise, params.vd=params.vd)
  vandantzig.ensemble.rcp45.ssns <- brick_vandantzig(currentyear, endyear, sea_level=sea_level_rcp45, time=time.proj, ss.gev=gev.params, surge.rise=ss.rcp45$surge.rise, params.vd=params.vd)
  vandantzig.ensemble.rcp60.ssns <- brick_vandantzig(currentyear, endyear, sea_level=sea_level_rcp60, time=time.proj, ss.gev=gev.params, surge.rise=ss.rcp60$surge.rise, params.vd=params.vd)
  vandantzig.ensemble.rcp85.ssns <- brick_vandantzig(currentyear, endyear, sea_level=sea_level_rcp85, time=time.proj, ss.gev=gev.params, surge.rise=ss.rcp85$surge.rise, params.vd=params.vd)

  ## Write output file
  dim.heightening <- ncdim_def('H', 'meters', vandantzig.ensemble.rcp26.ssst$Heightening[,1])
  dim.ensemble <- ncdim_def('ens', 'ensemble member', (1:ncol(vandantzig.ensemble.rcp26.ssst$Heightening)))
  dim.gev <- ncdim_def('param.gev', '[-]', (1:3))
  dim.vdparams <- ncdim_def('param.vd', '[-]', (1:6))

  surge.factor.params <- ncvar_def('surge_factor', '[-]', list(dim.ensemble), -999,
                  longname = 'storm surge level increase = surge_factor * sea level increase')
  vd.params <- ncvar_def('VD_params', '[-]', list(dim.ensemble, dim.vdparams), -999,
                  longname = 'Van Dantzig model parameters')
  gev.stat <- ncvar_def('gev_stat', '[-]', list(dim.ensemble, dim.gev), -999,
                  longname = 'stationary storm surge GEV parameters [loc, scale, shape]')

  cost.rcp26.nonstat <- ncvar_def('ExpectedCost_RCP26_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost (RCP26, non-stationary storm surge)')
  loss.rcp26.nonstat <- ncvar_def('ExpectedLoss_RCP26_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss (RCP26, non-stationary storm surge)')
  investment.rcp26.nonstat <- ncvar_def('ExpectedInvestment_RCP26_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment (RCP26, non-stationary storm surge)')
  preturn.rcp26.nonstat <- ncvar_def('ExpectedPreturn_RCP26_nonstat', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period (RCP26, non-stationary storm surge)')

  cost.rcp45.nonstat <- ncvar_def('ExpectedCost_RCP45_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost (RCP45, non-stationary storm surge)')
  loss.rcp45.nonstat <- ncvar_def('ExpectedLoss_RCP45_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss (RCP45, non-stationary storm surge)')
  investment.rcp45.nonstat <- ncvar_def('ExpectedInvestment_RCP45_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment (RCP45, non-stationary storm surge)')
  preturn.rcp45.nonstat <- ncvar_def('ExpectedPreturn_RCP45_nonstat', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period (RCP45, non-stationary storm surge)')

  cost.rcp60.nonstat <- ncvar_def('ExpectedCost_RCP60_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost (RCP60, non-stationary storm surge)')
  loss.rcp60.nonstat <- ncvar_def('ExpectedLoss_RCP60_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss (RCP60, non-stationary storm surge)')
  investment.rcp60.nonstat <- ncvar_def('ExpectedInvestment_RCP60_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment (RCP60, non-stationary storm surge)')
  preturn.rcp60.nonstat <- ncvar_def('ExpectedPreturn_RCP60_nonstat', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period (RCP60, non-stationary storm surge)')

  cost.rcp85.nonstat <- ncvar_def('ExpectedCost_RCP85_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost (RCP85, non-stationary storm surge)')
  loss.rcp85.nonstat <- ncvar_def('ExpectedLoss_RCP85_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss (RCP85, non-stationary storm surge)')
  investment.rcp85.nonstat <- ncvar_def('ExpectedInvestment_RCP85_nonstat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment (RCP85, non-stationary storm surge)')
  preturn.rcp85.nonstat <- ncvar_def('ExpectedPreturn_RCP85_nonstat', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period (RCP85, non-stationary storm surge)')

  cost.rcp26.stat <- ncvar_def('ExpectedCost_RCP26_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost (RCP26, stationary storm surge)')
  loss.rcp26.stat <- ncvar_def('ExpectedLoss_RCP26_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss (RCP26, stationary storm surge)')
  investment.rcp26.stat <- ncvar_def('ExpectedInvestment_RCP26_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment (RCP26, stationary storm surge)')
  preturn.rcp26.stat <- ncvar_def('ExpectedPreturn_RCP26_stat', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period (RCP26, stationary storm surge)')

  cost.rcp45.stat <- ncvar_def('ExpectedCost_RCP45_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost (RCP45, stationary storm surge)')
  loss.rcp45.stat <- ncvar_def('ExpectedLoss_RCP45_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss (RCP45, stationary storm surge)')
  investment.rcp45.stat <- ncvar_def('ExpectedInvestment_RCP45_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment (RCP45, stationary storm surge)')
  preturn.rcp45.stat <- ncvar_def('ExpectedPreturn_RCP45_stat', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period (RCP45, stationary storm surge)')

  cost.rcp60.stat <- ncvar_def('ExpectedCost_RCP60_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost (RCP60, stationary storm surge)')
  loss.rcp60.stat <- ncvar_def('ExpectedLoss_RCP60_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss (RCP60, stationary storm surge)')
  investment.rcp60.stat <- ncvar_def('ExpectedInvestment_RCP60_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment (RCP60, stationary storm surge)')
  preturn.rcp60.stat <- ncvar_def('ExpectedPreturn_RCP60_stat', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period (RCP60, stationary storm surge)')

  cost.rcp85.stat <- ncvar_def('ExpectedCost_RCP85_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost (RCP85, stationary storm surge)')
  loss.rcp85.stat <- ncvar_def('ExpectedLoss_RCP85_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss (RCP85, stationary storm surge)')
  investment.rcp85.stat <- ncvar_def('ExpectedInvestment_RCP85_stat', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment (RCP85, stationary storm surge)')
  preturn.rcp85.stat <- ncvar_def('ExpectedPreturn_RCP85_stat', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period (RCP85, stationary storm surge)')

  outnc <- nc_create(filename.vdout,
										list(surge.factor.params, vd.params, gev.stat,
										cost.rcp26.nonstat, loss.rcp26.nonstat, investment.rcp26.nonstat, preturn.rcp26.nonstat,
										cost.rcp45.nonstat, loss.rcp45.nonstat, investment.rcp45.nonstat, preturn.rcp45.nonstat,
										cost.rcp60.nonstat, loss.rcp60.nonstat, investment.rcp60.nonstat, preturn.rcp60.nonstat,
										cost.rcp85.nonstat, loss.rcp85.nonstat, investment.rcp85.nonstat, preturn.rcp85.nonstat,
										cost.rcp26.stat, loss.rcp26.stat, investment.rcp26.stat, preturn.rcp26.stat,
										cost.rcp45.stat, loss.rcp45.stat, investment.rcp45.stat, preturn.rcp45.stat,
										cost.rcp60.stat, loss.rcp60.stat, investment.rcp60.stat, preturn.rcp60.stat,
										cost.rcp85.stat, loss.rcp85.stat, investment.rcp85.stat, preturn.rcp85.stat),
										force_v4 = TRUE)

  ncvar_put(outnc, surge.factor.params, ss.rcp26$surge.factor)
  ncvar_put(outnc, vd.params, as.matrix(params.vd))
  ncvar_put(outnc, gev.stat, gev.params)

  # must check for possible 0% chance of failure, which gives infinite return
  # period and netCDF chokes. (these are for gigantic build strategies and do not
  # affect the scientific results, since the investment costs are enormous)
  pfail.tmp <- vandantzig.ensemble.rcp26.ssns$Average_p_fail
  itmp <- which(pfail.tmp==0); if(length(itmp)>0) {pfail.tmp[itmp]=min(pfail.tmp[-itmp])}
  ncvar_put(outnc, cost.rcp26.nonstat, vandantzig.ensemble.rcp26.ssns$Expected_costs)
  ncvar_put(outnc, loss.rcp26.nonstat, vandantzig.ensemble.rcp26.ssns$Expected_loss)
  ncvar_put(outnc, investment.rcp26.nonstat, vandantzig.ensemble.rcp26.ssns$Investment)
  ncvar_put(outnc, preturn.rcp26.nonstat, 1/pfail.tmp)

  pfail.tmp <- vandantzig.ensemble.rcp45.ssns$Average_p_fail
  itmp <- which(pfail.tmp==0); if(length(itmp)>0) {pfail.tmp[itmp]=min(pfail.tmp[-itmp])}
  ncvar_put(outnc, cost.rcp45.nonstat, vandantzig.ensemble.rcp45.ssns$Expected_costs)
  ncvar_put(outnc, loss.rcp45.nonstat, vandantzig.ensemble.rcp45.ssns$Expected_loss)
  ncvar_put(outnc, investment.rcp45.nonstat, vandantzig.ensemble.rcp45.ssns$Investment)
  ncvar_put(outnc, preturn.rcp45.nonstat, 1/pfail.tmp)

  pfail.tmp <- vandantzig.ensemble.rcp60.ssns$Average_p_fail
  itmp <- which(pfail.tmp==0); if(length(itmp)>0) {pfail.tmp[itmp]=min(pfail.tmp[-itmp])}
  ncvar_put(outnc, cost.rcp60.nonstat, vandantzig.ensemble.rcp60.ssns$Expected_costs)
  ncvar_put(outnc, loss.rcp60.nonstat, vandantzig.ensemble.rcp60.ssns$Expected_loss)
  ncvar_put(outnc, investment.rcp60.nonstat, vandantzig.ensemble.rcp60.ssns$Investment)
  ncvar_put(outnc, preturn.rcp60.nonstat, 1/pfail.tmp)

  pfail.tmp <- vandantzig.ensemble.rcp85.ssns$Average_p_fail
  itmp <- which(pfail.tmp==0); if(length(itmp)>0) {pfail.tmp[itmp]=min(pfail.tmp[-itmp])}
  ncvar_put(outnc, cost.rcp85.nonstat, vandantzig.ensemble.rcp85.ssns$Expected_costs)
  ncvar_put(outnc, loss.rcp85.nonstat, vandantzig.ensemble.rcp85.ssns$Expected_loss)
  ncvar_put(outnc, investment.rcp85.nonstat, vandantzig.ensemble.rcp85.ssns$Investment)
  ncvar_put(outnc, preturn.rcp85.nonstat, 1/pfail.tmp)

  pfail.tmp <- vandantzig.ensemble.rcp26.ssst$Average_p_fail
  itmp <- which(pfail.tmp==0); if(length(itmp)>0) {pfail.tmp[itmp]=min(pfail.tmp[-itmp])}
  ncvar_put(outnc, cost.rcp26.stat, vandantzig.ensemble.rcp26.ssst$Expected_costs)
  ncvar_put(outnc, loss.rcp26.stat, vandantzig.ensemble.rcp26.ssst$Expected_loss)
  ncvar_put(outnc, investment.rcp26.stat, vandantzig.ensemble.rcp26.ssst$Investment)
  ncvar_put(outnc, preturn.rcp26.stat, 1/pfail.tmp)

  pfail.tmp <- vandantzig.ensemble.rcp45.ssst$Average_p_fail
  itmp <- which(pfail.tmp==0); if(length(itmp)>0) {pfail.tmp[itmp]=min(pfail.tmp[-itmp])}
  ncvar_put(outnc, cost.rcp45.stat, vandantzig.ensemble.rcp45.ssst$Expected_costs)
  ncvar_put(outnc, loss.rcp45.stat, vandantzig.ensemble.rcp45.ssst$Expected_loss)
  ncvar_put(outnc, investment.rcp45.stat, vandantzig.ensemble.rcp45.ssst$Investment)
  ncvar_put(outnc, preturn.rcp45.stat, 1/pfail.tmp)

  pfail.tmp <- vandantzig.ensemble.rcp60.ssst$Average_p_fail
  itmp <- which(pfail.tmp==0); if(length(itmp)>0) {pfail.tmp[itmp]=min(pfail.tmp[-itmp])}
  ncvar_put(outnc, cost.rcp60.stat, vandantzig.ensemble.rcp60.ssst$Expected_costs)
  ncvar_put(outnc, loss.rcp60.stat, vandantzig.ensemble.rcp60.ssst$Expected_loss)
  ncvar_put(outnc, investment.rcp60.stat, vandantzig.ensemble.rcp60.ssst$Investment)
  ncvar_put(outnc, preturn.rcp60.stat, 1/pfail.tmp)

  pfail.tmp <- vandantzig.ensemble.rcp85.ssst$Average_p_fail
  itmp <- which(pfail.tmp==0); if(length(itmp)>0) {pfail.tmp[itmp]=min(pfail.tmp[-itmp])}
  ncvar_put(outnc, cost.rcp85.stat, vandantzig.ensemble.rcp85.ssst$Expected_costs)
  ncvar_put(outnc, loss.rcp85.stat, vandantzig.ensemble.rcp85.ssst$Expected_loss)
  ncvar_put(outnc, investment.rcp85.stat, vandantzig.ensemble.rcp85.ssst$Investment)
  ncvar_put(outnc, preturn.rcp85.stat, 1/pfail.tmp)

  nc_close(outnc)

}


##==============================================================================
## End
##==============================================================================
