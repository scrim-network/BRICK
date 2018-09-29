##==============================================================================
## processingPipeline_BRICKexperiments.R
##
## Pipeline for processing DAIS calibration results and BRICK rest-of-model
## calibration results.
## This is done for the standard BRICK experiment (control), as well as
## - the GSIC-SIMPLE experiment
## - the BRICK-GMSL experiment (this one only to create ensembles of hindcasts
##    and projections)
##
## 1. Data are read, hindcasts are set up. Parameters for the DAIS paleoclimatic
##    calibration and the DOECLIM+GSIC+GIS+TE modern calibration are drawn. The
##    number of parameter combinations (initial ensemble members) is specified
##    in the section of the script for the user to modify. (n.ensemble)
## 2. These full BRICK model parameter sets run the full model to obtain hindcasts
##    of global mean sea level. The parameters are calibrated to global mean sea
##    level data (Church and White, 2011). This is done using rejection sampling.
## 3. These fully calibrated parameter sets are written to a netCDF file whose
##    name is given by [filename.parameters], specified by the user.
## 4. Projections of sea level and its components are made to 2100. These, as
##    well as the hindcasts, are written to the file [filename.brickout]. This
##    file contains the physical model output.
## 5. If running the control experiment, projections of the components of sea-
##    level rise are fingerprinted to New Orleans, Louisiana. This yields local
##    sea-level rise for "NOLA".
## 6. If running the control experiment, local NOLA sea level is used in a
##
##  Required input: (set below)
##    filename.rho_simple_fixed  csv file with the fixed value for rho.simple
##    filename.DAIScalibration  DAIS calibration parameter posterior draws
##    filename.BRICKcalibration BRICK (non-DAIS) calibration parameter posterior draws
##    filename.parameters       [output] file name for post-calibrated parameters
##    filename.brickout         [output] file name for BRICK physical model results
##    filename.vdout            [output] file name for the Van Dantzig model results
##
##  Output:
##    BRICK-model_postcalibratedParameters_[datestamp].nc
##                            post-calibrated parameters file
##    BRICK-model_postcalibratedParameters_SIMPLE-GSIC_[datestamp].nc
##                            post-calibrated parameters file for SIMPLE-GSIC experiment
##    BRICK-model_drawcalibratedParameters_R07_[datestamp].nc
##                            calibrated parameter draws file for GMSL experiment
##                            (no post-calibration, but output the ensemble members)
##    BRICK-model_physical_[datestamp].nc
##                            netCDF4 file with the BRICK physical model output
##                            (should be as many runs as post-calibrated parameters).
##                            Includes SIMPLE-GSIC experimental results.
##                            Includes GMSL-R07 experimental results.
##    BRICK-model_VanDantzig_[datestamp].nc
##                            netCDF4 file with the BRICK physical model output
##                            (should be as many runs as post-calibrated parameters
##                            on the postcalibrated parameters file)
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

## Initial set-up
library(ncdf4)

t.beg = proc.time()

##==============================================================================
##==============================================================================

## Define some directories and settings

# BRICK working directory (the calibration one)
#setwd('~/codes/BRICK/calibration')
setwd('/home/scrim/axw322/codes/BRICK/calibration')

# Do the Van Dantzig with RCP8.5?
l.dovandantzig <- FALSE

# How many DAIS paleo simulations in your sample? (If you do a huge BRICK ensemble,
# it may be infeasible to do that number of 240,000 time step long simulations for
# DAIS paleo runs.)
n.samplepaleo <- 1e2

## Define the files you want to process/read/create

# Set up a filename for saving RData images along the way
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.saveprogress <- paste('brick_processing_',today,'.RData',sep='')

experiment='c'  ## Which model set-up? (c = control, with MAGICC-GSIC; e = experiment,
                ## with SIMPLE-GSIC; g = experiment, with BRICK-GMSL)
appen=''        ## Append file name? In case you process multiple files in one day
today=Sys.Date(); today=format(today,format="%d%b%Y")

l.aisfastdy = TRUE        # including AIS fast dynamics in the DAIS version used? (must be consistent with how DAIS_calib_driver.R was run)
n.ensemble = 1000      # total proposed ensemble before rejection sampling
n.ensemble.gmsl = 500    # pick n.ensemble for BRICK-GMSL to match control
n.ensemble.report = n.ensemble

if(experiment=='c'){
#  filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_12Aug2017.csv"
  filename.BRICKcalibration = "../output_calibration/BRICK-model_calibratedParameters_16Aug2017.nc"
  filename.DAIScalibration  = "../output_calibration/DAIS_calibratedParameters_17Aug2017.nc"
  filename.parameters       = paste('../output_calibration/BRICK-model_postcalibratedParameters_control_',today,appen,'.nc', sep="")
  filename.brickout         = paste('../output_model/BRICK-model_physical_control_',today,appen,'.nc',sep="")
  filename.vdout            = paste('../output_model/VanDantzig_RCP85_control_',today,appen,'.nc',sep="")
}
if(experiment=='e'){
  filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_01Nov2016.csv"
  filename.BRICKcalibration = "../output_calibration/BRICK-model_calibratedParameters_gsic-simple_01Nov2016.nc"
  filename.DAIScalibration  = "../output_calibration/DAIS_calibratedParameters_11Aug2016.nc"
  filename.parameters       = paste('../output_calibration/BRICK-model_postcalibratedParameters_SIMPLE-GSIC_',today,appen,'.nc', sep="")
  filename.brickout         = paste('../output_model/BRICK-model_physical_SIMPLE-GSIC_',today,appen,'.nc',sep="")
  # note: do not produce a Van Dantzig analysis file for the SIMPLE-GSIC experiment
}
if(experiment=='g'){
  filename.rho_simple_fixed = NULL
  filename.BRICKcalibration = "../output_calibration/BRICK-model_calibratedParameters_R07_01Nov2016.nc"
  filename.DAIScalibration  = NULL
  filename.parameters       = paste('../output_calibration/BRICK-model_drawcalibratedParameters_R07_',today,appen,'.nc', sep="")
  filename.brickout         = paste('../output_model/BRICK-model_physical_R07_',today,appen,'.nc',sep="")
  # note: do not produce a Van Dantzig analysis file for the GMSL experiment
  n.ensemble = n.ensemble.gmsl
}

## Add heteroscedastic (observational) error into the hindcasts?
l.ar1.hetero = TRUE

## Fingerprints of sea-level rise sources on local sea-level rise
lat.fp = 29.95      # latitude of location to fingerprint local sea level rise (>0 is North, <0 is South)
lon.fp = -90.07      # longitude of location ... (>0 is East, <0 is West)

## Mean and standard deviation for sampling land water storage (LWS)
## contributions to GMSL
## -- Using IPCC AR5 (Church et al. 2013) --
#lws.mean <- 0.38           # mm/y
#lws.sd   <- (0.49-.26)/4   # mm/y (take the IPCC 5-95% range as +/-2sigma)
## -- Using Dieng et al 2015 (doi:10.1088/1748-9326/10/12/124010)--
lws.mean <- 0.30           # mm/y
lws.sd   <- 0.18           # mm/y

##==============================================================================
##==============================================================================
##
##  NOTHING BELOW HERE SHOULD NEED MODIFIED
##
##==============================================================================
##==============================================================================




##==============================================================================
##==============================================================================
# simulate stationary AR(1) process (approximate - faster, better convergence, and
# results not sensitive to use of this as opposed to exact AR1)
ar1.sim = function(N,rho1,sigma) {
  x = rep(NA,N)
  for(i in 2:N)
    if(length(sigma)>1) {
      x[1] = sigma[1]/sqrt(1-rho1^2)
      x[i] = rho1*x[i-1] + rnorm(1,sd=sigma[i])
    } else {
      x[1] = sigma/sqrt(1-rho1^2)
      x[i] = rho1*x[i-1] + rnorm(1,sd=sigma)
    }
  return(x)
}
##==============================================================================
##==============================================================================




##==============================================================================
##==============================================================================
## Combine calibrated parameters from DAIS-fast dynamics and BRICK-RoM

## Make posterior parameter draws to run an ensemble of simulations and make
## projections of SLR (include DAIS, from calibrated parameters)
## Note: DAIS is part of BRICK. that the first parameter set is BRICK "rest-of-
## model", just that part that is not DAIS. DAIS is calibrated based on paleo
## data, so dealt with separately.

## First, read in the calibrated model parameters
print(paste('Reading calibrated non-DAIS model parameters...'))
ncdata <- nc_open(filename.BRICKcalibration)
parameters.brick = ncvar_get(ncdata, 'BRICK_parameters')
parnames.brick = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters.brick = t(parameters.brick)
colnames(parameters.brick) = parnames.brick
print(paste('  ... read ',nrow(parameters.brick),' model parameters'))

## Check to make sure you are not asking for more ensemble members than there
## are parameters available.
n.ensemble = min( n.ensemble, nrow(parameters.brick))

if(experiment!='g') {
  ## Second, read in the DAIS calibrated parameters
  print(paste('Reading calibrated DAIS model parameters...'))
  ncdata <- nc_open(filename.DAIScalibration)
  parameters.dais = ncvar_get(ncdata, 'DAIS_parameters')
  parnames.dais = ncvar_get(ncdata, 'parnames')
  nc_close(ncdata)
  parameters.dais = t(parameters.dais)
  colnames(parameters.dais) = parnames.dais
  print(paste('  ... read ',nrow(parameters.dais),' DAIS model parameters'))

  ## Check to make sure you are not asking for more ensemble members than there
  ## are parameters available.
  n.ensemble = min( n.ensemble, nrow(parameters.dais))

  ## Draw parameters for calibrated DAIS model
  parameters.ensemble.dais = mat.or.vec(n.ensemble , ncol(parameters.dais) )
  ind.dais=sample( seq(1,nrow(parameters.dais)), size=n.ensemble, replace=FALSE)
  for (p in 1:ncol(parameters.dais)){
    for (i in 1:n.ensemble){
        #parameters.ensemble.dais[i,p] = rnorm( 1, mean=parameters.dais[ind.dais[i],p], sd=bandwidths.dais[1,p])
      parameters.ensemble.dais[i,p] = parameters.dais[ind.dais[i],p]
    }
  }
}

## Draw parameters for calibrated models
print(paste('Creating possible parameter combinations for calibrated models...'))

## Use the "rnorm" versions if you want to draw from the multivariate kernel
## density estimates. This is unnecessary when using ~100,000s of posterior
## parameters from MCMC, however.
parameters.ensemble.brick = mat.or.vec(n.ensemble , ncol(parameters.brick) )
ind.ensemble = sample( seq(1,nrow(parameters.brick)), size=n.ensemble, replace=FALSE)
for (p in 1:ncol(parameters.brick)){
  for (i in 1:n.ensemble){
    #parameters.ensemble.brick[i,p] = rnorm( 1, mean=parameters.brick[ind.ensemble[i],p], sd=bandwidths.brick[1,p])
    parameters.ensemble.brick[i,p] = parameters.brick[ind.ensemble[i],p]
  }
}
print(paste(' ... done creating parameter sets'))

if(experiment=='g') {
  ## Set up the parameter indices within a new "parnames"; original calibrated
  ## model parnames was saved (parnames.brick/dais) in case you need it later
  parameters = parameters.ensemble.brick
  parnames = parnames.brick
  rownames(parameters)=NULL
} else {
  ## Set up the parameter indices within a new "parnames", which includes the
  ## DAIS parameters; original calibrated model parnames was saved (above) in case
  ## you need it later
  parameters = cbind(parameters.ensemble.brick, parameters.ensemble.dais)
  parnames = c(parnames.brick, parnames.dais)
  rownames(parameters)=NULL
}

## Set up the model for hindcasts
l.project = FALSE
begyear = 1850
endyear = 2009
mod.time= begyear:endyear
begyear.norm = 1986
endyear.norm = 2005
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time=length(mod.time)

## Get the forcing data
forcing = read.csv( '../data/forcing_hindcast.csv', header=TRUE )

## Source some useful functions for manipulating data
source('../R/forcing_total.R')          # function to add up the total forcing
source('../R/compute_indices.R')    # function to determine the model and
                    # data indices for comparisons

## Source the model(s) and data
source('../fortran/R/doeclimF.R')      # the DOECLIM model - resets the mod.time
source('../fortran/R/daisanto_fastdynF.R')  # DAIS (Antarctic Ice Sheet) model, with 'anto'
                    # for translating Tg (surface temperature anomaly)
                    # to Ta (Antarctic temperature reduced to sea level)
                    # and Toc (Antarctic/high-latitude ocean temperature)
source('../fortran/R/GSIC_magiccF.R')   # the GSIC model
source('../fortran/R/brick_te_F.R')     # TE (thermal expansion) model
source('../fortran/R/brick_tee_F.R')    # TEE (explicit thermal expansion) model
source('../fortran/R/simpleF.R')        # GIS (Greenland Ice Sheet) model
source('../R/gmsl_r07.R')        # the GMSL model

## Read the data sets for hindcast comparisons
source('../calibration/DOECLIM_readData.R')
source('../calibration/GSIC_readData.R')
source('../calibration/SIMPLE_readData.R')
source('../calibration/DAIS_readData.R')
source('../calibration/TE_readData.R')

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.temp,midx.ocheat,midx.gis,midx.gsic,midx.sl)
names(midx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )
oidx.all        = list(oidx.temp,oidx.ocheat,oidx.gis,oidx.gsic,oidx.sl)
names(oidx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.temp, obs.ocheat, obs.gis, obs.gsic, obs.sl)
names(obs.all) = c(    "temp"  , "ocheat"  , "gis"  , "gsic"  , "sl" )
obs.err.all        = list( obs.temp.err, obs.ocheat.err, obs.gis.err, obs.gsic.err, obs.sl.err)
names(obs.err.all) = c(    "temp"      , "ocheat"      , "gis"      , "gsic"      , "sl"      )

## Set the indices for normalization that are consistent with each data set
ind.norm.data = data.frame(
    c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
    c(which(mod.time==1850),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
    c(which(mod.time==1870),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

## Set the indices of the initial condition for each sub-model
i0 = vector("list",nrow(ind.norm.data)); names(i0)=as.character(ind.norm.data[,1])

## GSIC initial conditions are actually relative to 1990 (Wigley and Raper 2005)
## Re-set these. The simulation is relative to 1990, but results and comparison
## to data is relative to 1960.
i0$gsic = which(mod.time==1990)

## GIS initial conditions are relative to 1961-1990
i0$gis = which(mod.time==1961)

luse.sneasy   = FALSE    # Simple Nonlinear EArth SYstem model (DOECLIM+CCM)
luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # thermal expansion contribution to SLR
luse.tee      = FALSE   # explicit thermal expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = TRUE    # Antarctic ice sheet model
luse.lws      = FALSE   # for now, use LWS sampler hard-coded
luse.brick = cbind(luse.sneasy,luse.doeclim, luse.gsic, luse.te, luse.tee, luse.simple, luse.dais, luse.lws)

## Source the appropriate BRICK model for your purposes
if(experiment=='c') {source('../R/BRICK_coupledModel.R')}
if(experiment=='e') {source('../R/BRICK_coupledModel_SIMPLE-GSIC.R')}
if(experiment=='g') {
  luse.sneasy   = FALSE    # Simple Nonlinear EArth SYstem model (DOECLIM+CCM)
  luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
  luse.gsic     = FALSE   # glaciers and small ice caps contribution to SLR
  luse.te       = FALSE   # thermal expansion contribution to SLR
  luse.tee	= FALSE   # explicit thermal expansion contribution to SLR
  luse.simple   = FALSE   # Greenland ice sheet model
  luse.dais     = FALSE   # Antarctic ice sheet model
  luse.gmsl     = TRUE    # Example of adding your own model component - GMSL
  luse.lws      = FALSE   # land water storage
  luse.brick = cbind(luse.sneasy,luse.doeclim, luse.gsic, luse.te, luse.tee, luse.simple, luse.dais, luse.gmsl, luse.lws)
  source('../R/BRICK_coupledModel_R07.R')
}

## Initialize matrix to store model ensemble output
brick.out = vector("list", n.ensemble)

## Initialize flag for possibly bad runs (DAIS+BRICK parameters could go wrong,
## because the other model components were calibrated without DAIS, and vice
## versa)
badruns = rep(0, n.ensemble)

## Run the sample, and enjoy a nice progress bar
print(paste('Starting ',n.ensemble,' model hindcasts...',sep=''))
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {
  brick.out[[i]] = brick_model(parameters.in     = as.numeric(parameters[i,]),
                               parnames.in       = parnames,
                               forcing.in        = forcing,
                               l.project         = l.project,
                               slope.Ta2Tg.in    = slope.Ta2Tg,
                               intercept.Ta2Tg.in= intercept.Ta2Tg,
                               mod.time          = mod.time,
                               ind.norm.data     = ind.norm.data,
                               ind.norm.sl       = ind.norm,
                               luse.brick        = luse.brick,
                               i0                = i0,
                               l.aisfastdy       = l.aisfastdy)
  setTxtProgressBar(pb, i)
}
close(pb)
print(paste(' ... done running model hindcasts'))

## Save progress so far via workspace image
save.image(file=filename.saveprogress)

## Before post-calibration, need to add the modeled statistical noise back in.
## Only using sea-level rise data, so only need to modify GSIC, GIS.
## using the statistical parameters for AR1, AR1 and Gaussian noise, respectively
## Do not do for AIS, because var.dais was fit to paleo data-model mismatch, not
## representative of the current era.

## Read rho.simple from file
rho.simple.fixed=NA
if(exists('filename.rho_simple_fixed') & experiment!='g') {rho.simple.fixed = as.numeric(read.csv(filename.rho_simple_fixed))}

## Gather the fields for each simulation (easy referencing for plotting and
## analysis)
slr.out    = mat.or.vec(n.ensemble,length(mod.time))
temp.out   = mat.or.vec(n.ensemble,length(mod.time))
ocheat.out = mat.or.vec(n.ensemble,length(mod.time))
gsic.out   = mat.or.vec(n.ensemble,length(mod.time))
gis.out    = mat.or.vec(n.ensemble,length(mod.time))
ais.out    = mat.or.vec(n.ensemble,length(mod.time))
ais.disint.out = mat.or.vec(n.ensemble,length(mod.time))
te.out     = mat.or.vec(n.ensemble,length(mod.time))

## Will also normalize the output to "ind.norm" (1961-1990? 1986-2005? (Mengel, IPCC))
slr.out.norm    = slr.out
temp.out.norm   = temp.out
ocheat.out.norm = ocheat.out
gsic.out.norm   = gsic.out
gis.out.norm    = gis.out
ais.out.norm    = ais.out
te.out.norm      = te.out

## And add statistical noise
slr.norm.stat    = slr.out
temp.norm.stat   = temp.out
ocheat.norm.stat = ocheat.out
gsic.norm.stat   = gsic.out
gis.norm.stat    = gis.out
ais.norm.stat    = ais.out
te.norm.stat     = te.out

## Go through each simulation and collect, normalize and add modeled error to
## the fields

print(paste('Starting to add up total sea level rise from model hindcasts...',sep=''))

pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {

  T0=parameters[i,match("T0",parnames)]
  H0=parameters[i,match("H0",parnames)]

  # set the results (note that we already have slr.out)
  temp.out[i,]   = brick.out[[i]]$doeclim.out$temp + T0
  ocheat.out[i,] = brick.out[[i]]$doeclim.out$ocheat + H0

  # Normalize the output to "ind.norm.data"
  temp.out.norm[i,]   = temp.out[i,]  -mean(temp.out[i,1:20])
  ocheat.out.norm[i,] = ocheat.out[i,]#-mean(ocheat.out[i,ind.norm.data[which(ind.norm.data[,1]=='ocheat'),2]:ind.norm.data[which(ind.norm.data[,1]=='ocheat'),3]])

  # Add the statistcal model for AR1 (or otherwise) noise
  sigma.T     =parameters[i,match("sigma.T"     ,parnames)]
  rho.T       =parameters[i,match("rho.T"       ,parnames)]
  sigma.H     =parameters[i,match("sigma.H"     ,parnames)]
  rho.H       =parameters[i,match("rho.H"       ,parnames)]
  err.temp   = rep(sigma.T,n.time); if(l.ar1.hetero) {err.temp[midx.temp]=sqrt(sigma.T^2 + obs.temp.err[oidx.temp]^2)}
  err.ocheat = rep(sigma.H,n.time); if(l.ar1.hetero) {err.ocheat[midx.ocheat]=sqrt(sigma.H^2+obs.ocheat.err[oidx.ocheat]^2)}
  temp.norm.stat[i,]   = temp.out.norm[i,]   + ar1.sim(n.time, rho.T, err.temp)
  ocheat.norm.stat[i,] = ocheat.out.norm[i,] + ar1.sim(n.time, rho.H, err.ocheat)

  if(experiment=='g') {
    slr.out[i,] = brick.out[[i]]$gmsl.out

    # Normalize the output to "ind.norm.data"
    slr.out.norm[i,] = slr.out[i,]  -mean(slr.out[i,ind.norm.data[which(ind.norm.data[,1]=='sl'),2]:ind.norm.data[which(ind.norm.data[,1]=='sl'),3]])

    # Add the statistcal model for AR1 (or otherwise) noise
    sigma.gmsl=parameters[i,match("sigma.gmsl",parnames)]
    rho.gmsl  =parameters[i,match("rho.gmsl"  ,parnames)]
    err.gmsl  =rep(sigma.gmsl,n.time); if(l.ar1.hetero) {err.gmsl[midx.sl]=sqrt(sigma.gmsl^2 + obs.sl.err[oidx.sl]^2)}
    slr.norm.stat[i,] = slr.out.norm[i,] + ar1.sim(n.time, rho.gmsl, err.gmsl)

  } else {
    if(experiment=='c') {gsic.out[i,]   = brick.out[[i]]$gsic.out}
    if(experiment=='e') {gsic.out[i,]   = brick.out[[i]]$gsic.out$sle.gis}
    gis.out[i,]    = brick.out[[i]]$simple.out$sle.gis
    ais.out[i,]    = brick.out[[i]]$dais.out$Vais
    ais.disint.out[i,] = brick.out[[i]]$dais.out$Vdisint
    te.out[i,]     = brick.out[[i]]$te.out

    # Normalize the output to "ind.norm.data"
    gsic.out.norm[i,]   = gsic.out[i,]  -mean(gsic.out[i,ind.norm.data[which(ind.norm.data[,1]=='gsic'),2]:ind.norm.data[which(ind.norm.data[,1]=='gsic'),3]])
    gis.out.norm[i,]    = gis.out[i,]   -mean(gis.out[i,ind.norm.data[which(ind.norm.data[,1]=='gis'),2]:ind.norm.data[which(ind.norm.data[,1]=='gis'),3]])
    ais.out.norm[i,]    = ais.out[i,]   -mean(ais.out[i,ind.norm.data[which(ind.norm.data[,1]=='ais'),2]:ind.norm.data[which(ind.norm.data[,1]=='ais'),3]])
    te.out.norm[i,]     = te.out[i,]    -mean(te.out[i,ind.norm.data[which(ind.norm.data[,1]=='te'),2]:ind.norm.data[which(ind.norm.data[,1]=='te'),3]])

    # Add the statistcal model for AR1 (or otherwise) noise
    sigma.gsic  =parameters[i,match("sigma.gsic"  ,parnames)]
    rho.gsic    =parameters[i,match("rho.gsic"    ,parnames)]
    sigma.simple=parameters[i,match("sigma.simple",parnames)]
    rho.simple  =parameters[i,match("rho.simple"  ,parnames)]
    var.dais    =parameters[i,match("var.dais"    ,parnames)]
    if(is.null(rho.simple) | is.na(rho.simple)) rho.simple=rho.simple.fixed

    err.gsic = rep(sigma.gsic,n.time); if(l.ar1.hetero) {err.gsic[midx.gsic]=sqrt(sigma.gsic^2+obs.gsic.err[oidx.gsic]^2)}
    err.gis = rep(sigma.simple,n.time); if(l.ar1.hetero) {err.gis[midx.gis]=sqrt(sigma.simple^2+obs.gis.err^2)}

    gsic.norm.stat[i,]   = gsic.out.norm[i,]   + ar1.sim(n.time, rho.gsic, err.gsic)
    gis.norm.stat[i,]    = gis.out.norm[i,]    + ar1.sim(n.time, rho.simple, err.gis)
    ais.norm.stat[i,]    = ais.out.norm[i,]    #+ rnorm(  n.time, mean=0,sd=sqrt(var.dais))
    te.norm.stat[i,]     = te.out.norm[i,]

    slr.norm.stat[i,] = gsic.norm.stat[i,] +
                        gis.norm.stat[i,]  +
                        ais.norm.stat[i,]  +
                        te.norm.stat[i,]
  }

  slr.norm.stat[i,] = slr.norm.stat[i,] - mean(slr.norm.stat[i,ind.norm.data[which(ind.norm.data[,1]=='sl'),2]:ind.norm.data[which(ind.norm.data[,1]=='sl'),3]])

  setTxtProgressBar(pb, i)
}
close(pb)
print(paste(' ... done adding up model hindcast sea level rise and contributions'))

## Save progress so far via workspace image
save.image(file=filename.saveprogress)

if(experiment=='g') {
  ind.survive = 1:nrow(slr.norm.stat)
} else {
    # Rejection sampling, with target distribution as the likelihood of the sea-
  # level rise data, proposing uniformly across the ensemble members.
  survive = rep(0, n.ensemble)

  ## Make sure SLR data are also normalized
  ibeg=which(obs.sl.time==mod.time[ind.norm.data[which(ind.norm.data[,1]=='sl'),2]])
  iend=which(obs.sl.time==mod.time[ind.norm.data[which(ind.norm.data[,1]=='sl'),3]])
  obs.sl = obs.sl - mean(obs.sl[ibeg:iend])

  # calibrate to the Church and White data with land water subtracted out
  # and uncertainties added in quadrature
  # assumed budget: TE+AIS+GIS+GSIC+LWS = GMSL
  # 1901-1990: –0.11 [–0.16 to –0.06] (5-95% range)
  lw.time.1900 <- 1900:1989
  i1900 <- which(mod.time==lw.time.1900[1]):which(mod.time==lw.time.1900[length(lw.time.1900)])
  lw.1900 <- (-0.11/1000)*(lw.time.1900 - 1900)
  lw.err.1900 <- (0.25*(-0.06--0.16)/1000)*sqrt(lw.time.1900 - lw.time.1900[1])
  # 1971-2010: 0.12 [0.03 to 0.22]
  lw.time.1970 <- 1970:2009
  i1970 <- which(mod.time==lw.time.1970[1]):which(mod.time==lw.time.1970[length(lw.time.1970)])
  lw.1970 <- (0.12/1000)*(lw.time.1970 - lw.time.1970[1])
  lw.err.1970 <- (0.25*(0.2-0.03)/1000)*sqrt(lw.time.1970 - lw.time.1970[1])
  # 1993-2010: 0.38 [0.26 to 0.49]
  lw.time.1992 <- 1992:2009
  i1992 <- which(mod.time==lw.time.1992[1]):which(mod.time==lw.time.1992[length(lw.time.1992)])
  lw.1992 <- (0.38/1000)*(lw.time.1992 - lw.time.1992[1])
  lw.err.1992 <- (0.25*(0.49-0.26)/1000)*sqrt(lw.time.1992 - lw.time.1992[1])

  # normalize, subtract and add error in quadrature
  obs.sl.lw.1900 <- obs.sl[which(obs.sl.time==lw.time.1900[1]):which(obs.sl.time==lw.time.1900[length(lw.time.1900)])] - obs.sl[which(obs.sl.time==lw.time.1900[1])]
  obs.sl.lw.1970 <- obs.sl[which(obs.sl.time==lw.time.1970[1]):which(obs.sl.time==lw.time.1970[length(lw.time.1970)])] - obs.sl[which(obs.sl.time==lw.time.1970[1])]
  obs.sl.lw.1992 <- obs.sl[which(obs.sl.time==lw.time.1992[1]):which(obs.sl.time==lw.time.1992[length(lw.time.1992)])] - obs.sl[which(obs.sl.time==lw.time.1992[1])]

  obs.sl.lw.1900 <- obs.sl.lw.1900 - lw.1900
  obs.sl.lw.1970 <- obs.sl.lw.1970 - lw.1970
  obs.sl.lw.1992 <- obs.sl.lw.1992 - lw.1992

  obs.sl.lw.err.1900 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1900[1]):which(obs.sl.time==lw.time.1900[length(lw.time.1900)])]^2 + lw.err.1900^2)
  obs.sl.lw.err.1970 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1970[1]):which(obs.sl.time==lw.time.1970[length(lw.time.1970)])]^2 + lw.err.1970^2)
  obs.sl.lw.err.1992 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1992[1]):which(obs.sl.time==lw.time.1992[length(lw.time.1992)])]^2 + lw.err.1992^2)

  # calculate likelihood as the product of the three independent likelihoods
  resid.1900 <- obs.sl.lw.1900 - obs.sl.lw.1900
  llik.1900 <- sum(dnorm(resid.1900, sd=obs.sl.lw.err.1900, log=TRUE))
  resid.1970 <- obs.sl.lw.1970 - obs.sl.lw.1970
  llik.1970 <- sum(dnorm(resid.1970, sd=obs.sl.lw.err.1970, log=TRUE))
  resid.1992 <- obs.sl.lw.1992 - obs.sl.lw.1992
  llik.1992 <- sum(dnorm(resid.1992, sd=obs.sl.lw.err.1992, log=TRUE))
  lik.max <- (llik.1900 + llik.1970 + llik.1992)/n.ensemble

  imod.1900 <- which(mod.time==lw.time.1900[1]):which(mod.time==lw.time.1900[length(lw.time.1900)])
  imod.1970 <- which(mod.time==lw.time.1970[1]):which(mod.time==lw.time.1970[length(lw.time.1970)])
  imod.1992 <- which(mod.time==lw.time.1992[1]):which(mod.time==lw.time.1992[length(lw.time.1992)])

  uni.rnd = log(runif(n.ensemble))
  for (i in 1:n.ensemble) {
    resid.1900 <- obs.sl.lw.1900 - (slr.norm.stat[i,imod.1900]-slr.norm.stat[i,imod.1900[1]])
    resid.1970 <- obs.sl.lw.1970 - (slr.norm.stat[i,imod.1970]-slr.norm.stat[i,imod.1970[1]])
    resid.1992 <- obs.sl.lw.1992 - (slr.norm.stat[i,imod.1992]-slr.norm.stat[i,imod.1992[1]])
    llik.1900 <- sum(dnorm(resid.1900, sd=obs.sl.lw.err.1900, log=TRUE))
    llik.1970 <- sum(dnorm(resid.1970, sd=obs.sl.lw.err.1970, log=TRUE))
    llik.1992 <- sum(dnorm(resid.1992, sd=obs.sl.lw.err.1992, log=TRUE))
    lik.mem <- llik.1900 + llik.1970 + llik.1992
    if( uni.rnd[i] <= lik.mem-lik.max) {survive[i]=1}
  }
  ind.survive = which( as.logical(survive))
  print(paste('Calibration to sea level data by rejection sampling leaves ',length(ind.survive),' full calibrated ensemble members',sep=''))
}

slr.out.good = slr.norm.stat[ind.survive,]
parameters.good = parameters[ind.survive,]
colnames(parameters.good) = parnames

## Save progress so far via workspace image
save.image(file=filename.saveprogress)

##==============================================================================
##==============================================================================

if(experiment!='g') {

## DAIS paleo runs with the post-calibrated parameters
## Already have dSL, SL, Tg.recon.

parameters=parameters.good

## How many members do you want in your ensemble?
n.ensemble = nrow(parameters)
n.parameters = ncol(parameters)
n.sample <- min(n.samplepaleo, n.ensemble)
ind.sample <- sample( 1:n.ensemble, size=n.sample, replace=FALSE)
parameters.sample <- parameters[ind.sample,]

n.paleo = length(SL)
dais.paleo = mat.or.vec(n.sample, n.paleo)
date = seq(-239999,16,1) #240 Kyr BP to 2100AD at one year intervals of the forcings
norm.period=c(1961,1990)
ibeg=which(date==(norm.period[1]-2000))
iend=which(date==(norm.period[2]-2000))
ind.norm.paleo=ibeg:iend
t.paleo = date

# only do the paleo hindcast for the control model
if(experiment=='c') {

print(paste('Starting ',n.sample,' DAIS paleo hindcast simulations, using the post-calibrated parameters ...',sep=''))
pb <- txtProgressBar(min=0,max=n.sample,initial=0,style=3);
for (i in 1:n.sample) {

  anto.a=parameters.sample[i,match("anto.a",parnames)]
  anto.b=parameters.sample[i,match("anto.b",parnames)]
  gamma =parameters.sample[i,match("gamma" ,parnames)]
  alpha =parameters.sample[i,match("alpha.dais" ,parnames)]
  mu =parameters.sample[i,match("mu" ,parnames)]
  nu =parameters.sample[i,match("nu" ,parnames)]
  P0 =parameters.sample[i,match("P0" ,parnames)]
  kappa =parameters.sample[i,match("kappa.dais" ,parnames)]
  f0 =parameters.sample[i,match("f0" ,parnames)]
  h0 =parameters.sample[i,match("h0" ,parnames)]
  c =parameters.sample[i,match("c" ,parnames)]
  b0 =parameters.sample[i,match("b0" ,parnames)]
  slope =parameters.sample[i,match("slope" ,parnames)]
  var.dais =parameters.sample[i,match("var.dais" ,parnames)]
  if(l.aisfastdy) {
    Tcrit =parameters.sample[i,match("Tcrit" ,parnames)]
    lambda =parameters.sample[i,match("lambda" ,parnames)]
  } else {
    Tcrit = NULL
    lambda = NULL
  }

  dais.tmp = daisanto_fastdynF(
                       anto.a=anto.a, anto.b=anto.b,
                       gamma=gamma  , alpha=alpha  ,
                       mu=mu        , nu=nu        ,
                       P0=P0        , kappa=kappa  ,
                       f0=f0        , h0=h0        ,
                       c=c          , b0=b0        ,
                       slope=slope  , l.aisfastdy=l.aisfastdy,
                       Tcrit=Tcrit  , lambda=lambda,
                       slope.Ta2Tg=slope.Ta2Tg, intercept.Ta2Tg=intercept.Ta2Tg,
                       Tg=Tg.recon  , SL=SL ,
                       dSL=dSL      , includes_dSLais=1
                       )

  # Subtract off the 1961-1990 normalization period
  dais.norm = dais.tmp$Vais - mean(dais.tmp$Vais[ind.norm.paleo])

  # Add the modeled error back in
  dais.paleo[i,] = dais.norm + rnorm(n.paleo, mean=0,sd=sqrt(var.dais))

  setTxtProgressBar(pb, i)
}
close(pb)

}

parameters.good = parameters
colnames(parameters.good) = parnames

## Save progress so far via workspace image
save.image(file=filename.saveprogress)

## Get 5-95% CI for hindcasts
## Initialize arrays for the output
dais.paleo.05 = rep(NA,length(date)); dais.paleo.50 = rep(NA,length(date)); dais.paleo.95 = rep(NA,length(date))
dais.paleo.max= rep(NA,length(date)); dais.paleo.min= rep(NA,length(date));

## Source a useful script, to allow for assigning multiple outputs at once
source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
if(experiment=='c'){
pb <- txtProgressBar(min=0,max=length(date),initial=0,style=3);
for (t in 1:length(date)){
  c(dais.paleo.05[t] , dais.paleo.50[t] , dais.paleo.95[t] , dais.paleo.max[t], dais.paleo.min[t]) := quantile(dais.paleo[,t],c(0.05,.50,.95,1,0), na.rm=TRUE)
  setTxtProgressBar(pb, t)
}
close(pb)
print(paste(' ... done with the paleo simulations!',sep=''))
}

## Make smoothed version of the AIS paleo results, so the plots are not massive
n.avg=100  # number of years in averaging period
n.time.avg=ceiling(length(date)/n.avg)
dais.paleo.05.avg = rep(NA,n.time.avg)
dais.paleo.50.avg = rep(NA,n.time.avg)
dais.paleo.95.avg = rep(NA,n.time.avg)
dais.paleo.max.avg = rep(NA,n.time.avg)
dais.paleo.min.avg = rep(NA,n.time.avg)
date.avg=seq(date[1],date[length(date)],by=n.avg)

if(experiment=='c'){
print(paste('Smoothing the paleo simulations with ',n.avg,'-year averages ...',sep=''))
pb <- txtProgressBar(min=0,max=n.time.avg,initial=0,style=3);
for (t in 1:(n.time.avg-1)){
  dais.paleo.05.avg[t] = mean(dais.paleo.05[((t-1)*n.avg+1) : (t*n.avg)])
  dais.paleo.50.avg[t] = mean(dais.paleo.50[((t-1)*n.avg+1) : (t*n.avg)])
  dais.paleo.95.avg[t] = mean(dais.paleo.95[((t-1)*n.avg+1) : (t*n.avg)])
  dais.paleo.max.avg[t] = mean(dais.paleo.max[((t-1)*n.avg+1) : (t*n.avg)])
  dais.paleo.min.avg[t] = mean(dais.paleo.min[((t-1)*n.avg+1) : (t*n.avg)])
  setTxtProgressBar(pb, t)
}
dais.paleo.05.avg[n.time.avg] = mean(dais.paleo.05[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.50.avg[n.time.avg] = mean(dais.paleo.50[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.95.avg[n.time.avg] = mean(dais.paleo.95[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.max.avg[n.time.avg] = mean(dais.paleo.max[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.min.avg[n.time.avg] = mean(dais.paleo.min[((n.time.avg-1)*n.avg+1) : length(date)])
close(pb)
print(paste(' ... done smoothing the paleo simulations!',sep=''))
}

}  # end check if(experiment!='g') (that one does not use DAIS)

##==============================================================================
##==============================================================================

## Save the hindcasts, trimming down first for ind.survive (Church and White, 2011)
## and then for AIS Vmin
gsic.hind = t(gsic.norm.stat[ind.survive,])
te.hind = t(te.norm.stat[ind.survive,])
gis.hind = t(gis.norm.stat[ind.survive,])
ais.hind = t(ais.norm.stat[ind.survive,])
temp.hind = t(temp.norm.stat[ind.survive,])
ocheat.hind = t(ocheat.norm.stat[ind.survive,])
gsl.hind = t(slr.out.good)
t.hind = mod.time

## Save progress so far via workspace image
save.image(file=filename.saveprogress)

##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## Write post-calibrated parameter netCDF output file
## Note: benefit of netCDF is that multiple post-calibrated parameters files can
## be concatenated along the "unlimited" dimension of n.ensemble, to yield larger
## ensembles and more robust statistics.

print(paste("Writing post-calibrated parameters to file ",filename.parameters,sep=""))

lmax=0
for (i in 1:length(parnames)){lmax=max(lmax,nchar(parnames[i]))}

library(ncdf4)
dim.parameters <- ncdim_def('n.parameters', '', 1:ncol(parameters.good), unlim=FALSE)
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(parameters.good), unlim=TRUE)
parameters.var <- ncvar_def('BRICK_parameters', '', list(dim.parameters,dim.ensemble), -999)
parnames.var <- ncvar_def('parnames', '', list(dim.name,dim.parameters), prec='char')
outnc <- nc_create(filename.parameters, list(parameters.var,parnames.var))
ncvar_put(outnc, parameters.var, t(parameters.good))
ncvar_put(outnc, parnames.var, parnames)
nc_close(outnc)

##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================

## Make RCP2.6, 4.5 and 8.5 projections to 2100. Need:
## (1) global total sea level, (2) global sea level without fast dynamics,
## (3) local (NOLA) sea level, (4) local sea level without fast dynamics

parameters=parameters.good
parnames = colnames(parameters)

n.ensemble = nrow(parameters)      # ... or take all of them
n.parameters = ncol(parameters)    # number of distinct model parameters

## Display output to let the user know what is going on, and how many ensemble
## members will be processed into projections.
print(paste('Starting projections for ensemble of ',n.ensemble,' members...',sep=''))

## Set up the model projections
l.project = TRUE
begyear = 1850
endyear = 2100; #if(!l.project & endyear>2009) print('l.project and endyear not compatible')
mod.time= begyear:endyear
t.proj = mod.time # save for output

## NOTE: IPCC generally are relative to 1986-2005.
## Therefore use that period to be commensurate with their results
begyear.norm = 1986
endyear.norm = 2005
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time = length(mod.time)

## Want model projections under RCP2.6, 4.5 and 8.5. Initialize a list, each
## will hold the ensemble model output for a different forcing scenario.
n.scen = 3
proj.out = vector("list", n.scen)
forc.scen = vector("list", n.scen)
badruns.scen = vector("list", n.scen)
names.scen = c('rcp26','rcp45','rcp85')
names(proj.out) = names.scen

## Get the forcing data
forc.scen[[1]] = read.csv( '../data/forcing_rcp26.csv', header=TRUE )
forc.scen[[2]] = read.csv( '../data/forcing_rcp45.csv', header=TRUE )
forc.scen[[3]] = read.csv( '../data/forcing_rcp85.csv', header=TRUE )

## Loop over forcing scenarios
for (ff in 1:n.scen) {

  forcing = forc.scen[[ff]]

  ## Initialize matrix to store model ensemble output and flag for bad runs
  ## (could go bad because of different forcing from hindcasts; different
  ## temperatures lead to different stability requirements for SIMPLE)
  brick.out = vector("list", n.ensemble)
  badruns   = rep(0, n.ensemble)

  ## Run the sample, and enjoy a nice progress bar
  pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3);
  for (i in 1:n.ensemble) {

    brick.out[[i]] = brick_model(parameters.in     = parameters[i,],
                                 parnames.in       = parnames,
                                 forcing.in        = forcing,
                                 l.project         = l.project,
                                 slope.Ta2Tg.in    = slope.Ta2Tg,
                                 intercept.Ta2Tg.in= intercept.Ta2Tg,
                                 mod.time          = mod.time,
                                 ind.norm.data     = ind.norm.data,
                                 ind.norm.sl       = ind.norm,
                                 luse.brick        = luse.brick,
                                 i0                = i0,
                                 l.aisfastdy       = l.aisfastdy)

    # check if the run turned out bad
    if( is.na(brick.out[[i]]$slr.out[length(mod.time)]) ) {badruns[i]=1}

    setTxtProgressBar(pb, i)
  }
  close(pb)

  badruns.scen[[ff]] = badruns
  proj.out[[ff]] = brick.out
}

## Save progress so far via workspace image
save.image(file=filename.saveprogress)

## Filter out any bad runs
ind.good = which(badruns.scen[[1]]==0 & badruns.scen[[2]]==0 & badruns.scen[[3]]==0)
n.ensemble = length(ind.good)
parameters = parameters[ind.good,]

## Trim down the hindcasts to match (some runs might go off the rails with RCP
## forcings but didn't with the hindcast forcings)
gsl.hind = gsl.hind[,ind.good]
gsic.hind = gsic.hind[,ind.good]
te.hind = te.hind[,ind.good]
gis.hind = gis.hind[,ind.good]
ais.hind = ais.hind[,ind.good]
temp.hind = temp.hind[,ind.good]
ocheat.hind = ocheat.hind[,ind.good]

brick.rcp26 = proj.out[[1]][ind.good]
brick.rcp45 = proj.out[[2]][ind.good]
brick.rcp85 = proj.out[[3]][ind.good]

## Gather the fields for each simulation (easy referencing for plotting and
## analysis)
proj.rcp26 = vector("list", n.scen)
proj.rcp45 = vector("list", n.scen)
proj.rcp85 = vector("list", n.scen)

## Make each projections list a list of the SLR contributions and other fields
## Get the 90% CI at each time step, of all model output fields

if(experiment=='g') {
  names.output = c('slr','temp','ocheat')
} else {
  names.output = c('slr','gsic','gis','ais','disint','te','lws','temp','ocheat','slr.nola')
}
n.output = length(names.output)

proj.rcp26 = vector("list", n.output)
proj.rcp45 = vector("list", n.output)
proj.rcp85 = vector("list", n.output)

names(proj.rcp26) = names.output
names(proj.rcp45) = names.output
names(proj.rcp85) = names.output

for (j in 1:n.output) {
  proj.rcp26[[j]] = mat.or.vec(n.ensemble, n.time)
  proj.rcp45[[j]] = mat.or.vec(n.ensemble, n.time)
  proj.rcp85[[j]] = mat.or.vec(n.ensemble, n.time)
}

## Go through each simulation and collect, normalize and add modeled error to
## the fields

pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)

for (i in 1:n.ensemble) {

  T0=parameters[i,match("T0",parnames)]
  H0=parameters[i,match("H0",parnames)]

  # set the results (note that we already have slr.out)
  proj.rcp26$temp[i,] = brick.rcp26[[i]]$doeclim.out$temp + T0
  proj.rcp26$ocheat[i,] = brick.rcp26[[i]]$doeclim.out$ocheat + H0

  proj.rcp45$temp[i,] = brick.rcp45[[i]]$doeclim.out$temp + T0
  proj.rcp45$ocheat[i,] = brick.rcp45[[i]]$doeclim.out$ocheat + H0

  proj.rcp85$temp[i,] = brick.rcp85[[i]]$doeclim.out$temp + T0
  proj.rcp85$ocheat[i,] = brick.rcp85[[i]]$doeclim.out$ocheat + H0

  # Normalize the output to "ind.norm" (1961-1990? 1986-2005 (Mengel)?).
  # Normalize ocean heat uptake too, for sake of plotting (it will be plotted as
  # the amount of heat taken up by ocean since 1986-2005 period)
  proj.rcp26$temp[i,] = proj.rcp26$temp[i,] - mean( proj.rcp26$temp[i,ind.norm])
  proj.rcp26$ocheat[i,] = proj.rcp26$ocheat[i,] - mean( proj.rcp26$ocheat[i,ind.norm])

  proj.rcp45$temp[i,] = proj.rcp45$temp[i,] - mean( proj.rcp45$temp[i,ind.norm])
  proj.rcp45$ocheat[i,] = proj.rcp45$ocheat[i,] - mean( proj.rcp45$ocheat[i,ind.norm])

  proj.rcp85$temp[i,] = proj.rcp85$temp[i,] - mean( proj.rcp85$temp[i,ind.norm])
  proj.rcp85$ocheat[i,] = proj.rcp85$ocheat[i,] - mean( proj.rcp85$ocheat[i,ind.norm])

  # Add the statistcal model for AR1 (or otherwise) noise?
  # Edit: not for DAIS, because that noise matches paleo data, which has
  # different uncertainties than modern data.
  sigma.T  =parameters[i,match("sigma.T"  ,parnames)]
  rho.T    =parameters[i,match("rho.T"    ,parnames)]
  sigma.H  =parameters[i,match("sigma.H"  ,parnames)]
  rho.H    =parameters[i,match("rho.H"    ,parnames)]

  proj.rcp26$temp[i,] = proj.rcp26$temp[i,] + ar1.sim(n.time, rho.T, sigma.T)
  proj.rcp26$ocheat[i,] = proj.rcp26$ocheat[i,] + ar1.sim(n.time, rho.H, sigma.H)

  proj.rcp45$temp[i,] = proj.rcp45$temp[i,] + ar1.sim(n.time, rho.T, sigma.T)
  proj.rcp45$ocheat[i,] = proj.rcp45$ocheat[i,] + ar1.sim(n.time, rho.H, sigma.H)

  proj.rcp85$temp[i,] = proj.rcp85$temp[i,] + ar1.sim(n.time, rho.T, sigma.T)
  proj.rcp85$ocheat[i,] = proj.rcp85$ocheat[i,] + ar1.sim(n.time, rho.H, sigma.H)

  # set projected sea-level if GMSL experiment
  if(experiment=='g') {
    proj.rcp26$slr[i,] = brick.rcp26[[i]]$gmsl.out
    proj.rcp45$slr[i,] = brick.rcp45[[i]]$gmsl.out
    proj.rcp85$slr[i,] = brick.rcp85[[i]]$gmsl.out
  }

  if(experiment!='g') {

    # set the results (note that we already have slr.out)
    if(experiment=='c') {proj.rcp26$gsic[i,] = brick.rcp26[[i]]$gsic.out}
    if(experiment=='e') {proj.rcp26$gsic[i,] = brick.rcp26[[i]]$gsic.out$sle.gis}
    proj.rcp26$gis[i,] = brick.rcp26[[i]]$simple.out$sle.gis
    proj.rcp26$ais[i,] = brick.rcp26[[i]]$dais.out$Vais
    proj.rcp26$disint[i,] = brick.rcp26[[i]]$dais.out$Vdisint
    proj.rcp26$te[i,] = brick.rcp26[[i]]$te.out

    if(experiment=='c') {proj.rcp45$gsic[i,] = brick.rcp45[[i]]$gsic.out}
    if(experiment=='e') {proj.rcp45$gsic[i,] = brick.rcp45[[i]]$gsic.out$sle.gis}
    proj.rcp45$gis[i,] = brick.rcp45[[i]]$simple.out$sle.gis
    proj.rcp45$ais[i,] = brick.rcp45[[i]]$dais.out$Vais
    proj.rcp45$disint[i,] = brick.rcp45[[i]]$dais.out$Vdisint
    proj.rcp45$te[i,] = brick.rcp45[[i]]$te.out

    if(experiment=='c') {proj.rcp85$gsic[i,] = brick.rcp85[[i]]$gsic.out}
    if(experiment=='e') {proj.rcp85$gsic[i,] = brick.rcp85[[i]]$gsic.out$sle.gis}
    proj.rcp85$gis[i,] = brick.rcp85[[i]]$simple.out$sle.gis
    proj.rcp85$ais[i,] = brick.rcp85[[i]]$dais.out$Vais
    proj.rcp85$disint[i,] = brick.rcp85[[i]]$dais.out$Vdisint
    proj.rcp85$te[i,] = brick.rcp85[[i]]$te.out

    # Normalize the output to "ind.norm" (1961-1990? 1986-2005 (Mengel)?).
    # Normalize ocean heat uptake too, for sake of plotting (it will be plotted as
    # the amount of heat taken up by ocean since 1986-2005 period)
    proj.rcp26$gsic[i,] = proj.rcp26$gsic[i,] - mean( proj.rcp26$gsic[i,ind.norm])
    proj.rcp26$gis[i,] = proj.rcp26$gis[i,] - mean( proj.rcp26$gis[i,ind.norm])
    proj.rcp26$ais[i,] = proj.rcp26$ais[i,] - mean( proj.rcp26$ais[i,ind.norm])
    proj.rcp26$te[i,] = proj.rcp26$te[i,] - mean( proj.rcp26$te[i,ind.norm])

    proj.rcp45$gsic[i,] = proj.rcp45$gsic[i,] - mean( proj.rcp45$gsic[i,ind.norm])
    proj.rcp45$gis[i,] = proj.rcp45$gis[i,] - mean( proj.rcp45$gis[i,ind.norm])
    proj.rcp45$ais[i,] = proj.rcp45$ais[i,] - mean( proj.rcp45$ais[i,ind.norm])
    proj.rcp45$te[i,] = proj.rcp45$te[i,] - mean( proj.rcp45$te[i,ind.norm])

    proj.rcp85$gsic[i,] = proj.rcp85$gsic[i,] - mean( proj.rcp85$gsic[i,ind.norm])
    proj.rcp85$gis[i,] = proj.rcp85$gis[i,] - mean( proj.rcp85$gis[i,ind.norm])
    proj.rcp85$ais[i,] = proj.rcp85$ais[i,] - mean( proj.rcp85$ais[i,ind.norm])
    proj.rcp85$te[i,] = proj.rcp85$te[i,] - mean( proj.rcp85$te[i,ind.norm])

    # Add the statistcal model for AR1 (or otherwise) noise?
    # Edit: not for DAIS, because that noise matches paleo data, which has
    # different uncertainties than modern data.
    sigma.gsic  =parameters[i,match("sigma.gsic"  ,parnames)]
    rho.gsic    =parameters[i,match("rho.gsic"    ,parnames)]
    sigma.simple=parameters[i,match("sigma.simple",parnames)]
    rho.simple  =parameters[i,match("rho.simple"  ,parnames)]
    var.dais    =parameters[i,match("var.dais"    ,parnames)]
    if(is.na(rho.simple)) rho.simple=rho.simple.fixed

    proj.rcp26$gsic[i,] = proj.rcp26$gsic[i,] + ar1.sim(n.time, rho.gsic, sigma.gsic)
    proj.rcp26$gis[i,] = proj.rcp26$gis[i,] + ar1.sim(n.time, rho.simple, sigma.simple)

    proj.rcp45$gsic[i,] = proj.rcp45$gsic[i,] + ar1.sim(n.time, rho.gsic, sigma.gsic)
    proj.rcp45$gis[i,] = proj.rcp45$gis[i,] + ar1.sim(n.time, rho.simple, sigma.simple)

    proj.rcp85$gsic[i,] = proj.rcp85$gsic[i,] + ar1.sim(n.time, rho.gsic, sigma.gsic)
    proj.rcp85$gis[i,] = proj.rcp85$gis[i,] + ar1.sim(n.time, rho.simple, sigma.simple)

    # Add contributions to land water storage. /1000 to convert to meters.
    # This is done for each ensemble member and each RCP.
    # All other contributions are normalized to 1986-2005 (some have the
    # estimated noise added, so mean(1986-2005) not nec. equall to 0), so
    # normalize the lws.est contribution to this period.

    # RCP2.6
    proj.rcp26$lws[i,] <- cumsum(rnorm(n=n.time, mean=lws.mean, sd=lws.sd)) /1000
    proj.rcp26$lws[i,] <- proj.rcp26$lws[i,] - mean(proj.rcp26$lws[i,ind.norm])

    # RCP4.5
    proj.rcp45$lws[i,] <- cumsum(rnorm(n=n.time, mean=lws.mean, sd=lws.sd)) /1000
    proj.rcp45$lws[i,] <- proj.rcp45$lws[i,] - mean(proj.rcp45$lws[i,ind.norm])

    # RCP8.5
    proj.rcp85$lws[i,] <- cumsum(rnorm(n=n.time, mean=lws.mean, sd=lws.sd)) /1000
    proj.rcp85$lws[i,] <- proj.rcp85$lws[i,] - mean(proj.rcp85$lws[i,ind.norm])

    # Add up to total sea-level rise
    proj.rcp26$slr[i,] = proj.rcp26$gsic[i,] +
                         proj.rcp26$gis[i,] +
                         proj.rcp26$ais[i,] +
                         proj.rcp26$te[i,] +
                         proj.rcp26$lws[i,]

    proj.rcp45$slr[i,] = proj.rcp45$gsic[i,] +
                         proj.rcp45$gis[i,] +
                         proj.rcp45$ais[i,] +
                         proj.rcp45$te[i,]
                         proj.rcp26$lws[i,]

    proj.rcp85$slr[i,] = proj.rcp85$gsic[i,] +
                         proj.rcp85$gis[i,] +
                         proj.rcp85$ais[i,] +
                         proj.rcp85$te[i,]
                         proj.rcp26$lws[i,]

  } # end if(experiment!='g')

  # And normalize sea-level rise
  proj.rcp26$slr[i,] = proj.rcp26$slr[i,] - mean(proj.rcp26$slr[i,ind.norm])
  proj.rcp45$slr[i,] = proj.rcp45$slr[i,] - mean(proj.rcp45$slr[i,ind.norm])
  proj.rcp85$slr[i,] = proj.rcp85$slr[i,] - mean(proj.rcp85$slr[i,ind.norm])

  setTxtProgressBar(pb, i)
}
close(pb)

## Save progress so far via workspace image
save.image(file=filename.saveprogress)

print(paste('... finished projections!',sep=''))

##==============================================================================
##==============================================================================





if(experiment!='g') {
##==============================================================================
##==============================================================================
## Local sea level rise
## Note: leaving out land water storage. Minor contribution, and no good local
## estimate.

print(paste('Beginning fingerprinting to local sea level rise...',sep=''))

source('../R/BRICK_LSL.R')

proj.rcp26$slr.nola = brick_lsl(lat.in=lat.fp,
                                lon.in=lon.fp,
                n.time=length(t.proj),
                slr_gis = proj.rcp26$gis,
                slr_gsic = proj.rcp26$gsic,
                slr_ais = proj.rcp26$ais,
                slr_te = proj.rcp26$te)
proj.rcp45$slr.nola = brick_lsl(lat.in=lat.fp,
                                lon.in=lon.fp,
                n.time=length(t.proj),
                slr_gis = proj.rcp45$gis,
                slr_gsic = proj.rcp45$gsic,
                slr_ais = proj.rcp45$ais,
                slr_te = proj.rcp45$te)
proj.rcp85$slr.nola = brick_lsl(lat.in=lat.fp,
                                lon.in=lon.fp,
                n.time=length(t.proj),
                slr_gis = proj.rcp85$gis,
                slr_gsic = proj.rcp85$gsic,
                slr_ais = proj.rcp85$ais,
                slr_te = proj.rcp85$te)

# And normalize sea-level rise
for (i in 1:n.ensemble) {
  proj.rcp26$slr.nola[i,] = proj.rcp26$slr.nola[i,] - mean(proj.rcp26$slr.nola[i,ind.norm])
  proj.rcp45$slr.nola[i,] = proj.rcp45$slr.nola[i,] - mean(proj.rcp45$slr.nola[i,ind.norm])
  proj.rcp85$slr.nola[i,] = proj.rcp85$slr.nola[i,] - mean(proj.rcp85$slr.nola[i,ind.norm])
}

## Save progress so far via workspace image
save.image(file=filename.saveprogress)

print(paste('... finished local sea level rise',sep=''))

##==============================================================================
##==============================================================================
}






##==============================================================================
##==============================================================================
## Write a netCDF ensemble output file including each of the RCP scenarios:
## (1) global total sea level, (2) local (NOLA) sea level. Also will want each
## contribution to global sea level rise, for the hindcast plots

library(ncdf4)

dim.tproj <- ncdim_def('time_proj', 'years', as.double(t.proj))
#dim.ensemble <- ncdim_def('ens', 'ensemble member', as.double(1:nrow(proj.rcp26$slr)))
dim.ensemble <- ncdim_def('ens', 'ensemble member', as.double(1:nrow(proj.rcp26$slr)), unlim=TRUE)
dim.thind <- ncdim_def('time_hind', 'years', as.double(t.hind))

if(experiment!='g') {
  dim.tpaleo <- ncdim_def('time_paleo', 'year paleo', as.double(t.paleo))
  dim.tpaleo.avg <- ncdim_def('time_paleo_avg', 'year avg paleo', as.double(date.avg))

  ais.paleo.05 <- ncvar_def('AIS_paleo_q05', 'meters', list(dim.tpaleo), -999,
                  longname = 'AIS paleo contribution to sea level (5th quantile)')
  ais.paleo.50 <- ncvar_def('AIS_paleo_q50', 'meters', list(dim.tpaleo), -999,
                  longname = 'AIS paleo contribution to sea level (median)')
  ais.paleo.95 <- ncvar_def('AIS_paleo_q95', 'meters', list(dim.tpaleo), -999,
                  longname = 'AIS paleo contribution to sea level (95th quantile)')
  ais.paleo.max <- ncvar_def('AIS_paleo_max', 'meters', list(dim.tpaleo), -999,
                  longname = 'AIS paleo contribution to sea level (maximum)')
  ais.paleo.min <- ncvar_def('AIS_paleo_min', 'meters', list(dim.tpaleo), -999,
                  longname = 'AIS paleo contribution to sea level (minimum)')

  ais.paleo.05.avg <- ncvar_def('AIS_paleo_avg_q05', 'meters', list(dim.tpaleo.avg), -999,
                  longname = 'AIS paleo contribution to sea level (smoothed, 5th quantile)')
  ais.paleo.50.avg <- ncvar_def('AIS_paleo_avg_q50', 'meters', list(dim.tpaleo.avg), -999,
                  longname = 'AIS paleo contribution to sea level (smoothed, median)')
  ais.paleo.95.avg <- ncvar_def('AIS_paleo_avg_q95', 'meters', list(dim.tpaleo.avg), -999,
                  longname = 'AIS paleo contribution to sea level (smoothed, 95th quantile)')
  ais.paleo.max.avg <- ncvar_def('AIS_paleo_avg_max', 'meters', list(dim.tpaleo.avg), -999,
                  longname = 'AIS paleo contribution to sea level (smoothed, maximum)')
  ais.paleo.min.avg <- ncvar_def('AIS_paleo_avg_min', 'meters', list(dim.tpaleo.avg), -999,
                  longname = 'AIS paleo contribution to sea level (smoothed, minimum)')

  gsic.hindcast <- ncvar_def('GSIC_hind', 'meters', list(dim.thind, dim.ensemble), -999,
                      longname = 'Glaciers and small ice caps contribution to GSL (hindcast)')
  te.hindcast <- ncvar_def('TE_hind', 'meters', list(dim.thind, dim.ensemble), -999,
                      longname = 'Thermal expansion contribution to GSL (hindcast)')
  gis.hindcast <- ncvar_def('GIS_hind', 'meters', list(dim.thind, dim.ensemble), -999,
                      longname = 'Greenland ice sheet contribution to GSL (hindcast)')
  ais.hindcast <- ncvar_def('AIS_hind', 'meters', list(dim.thind, dim.ensemble), -999,
                      longname = 'Antarctic ice sheet contribution to GSL (hindcast)')

  lsl.rcp26 <- ncvar_def('LocalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level (RCP26)')
  gsic.rcp26 <- ncvar_def('GSIC_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GSIC contribution to sea level (RCP26)')
  te.rcp26 <- ncvar_def('TE_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'TE contribution to sea level (RCP26)')
  gis.rcp26 <- ncvar_def('GIS_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GIS contribution to sea level (RCP26)')
  ais.rcp26 <- ncvar_def('AIS_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'AIS contribution to sea level (RCP26)')
  disint.rcp26 <- ncvar_def('AIS_DISINT_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'AIS fast dynamics contribution to sea level (RCP26)')
  lws.rcp26 <- ncvar_def('LWS_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'estimated LWS contribution to sea level (RCP26)')

  lsl.rcp45 <- ncvar_def('LocalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level (RCP45)')
  gsic.rcp45 <- ncvar_def('GSIC_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GSIC contribution to sea level (RCP45)')
  te.rcp45 <- ncvar_def('TE_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'TE contribution to sea level (RCP45)')
  gis.rcp45 <- ncvar_def('GIS_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GIS contribution to sea level (RCP45)')
  ais.rcp45 <- ncvar_def('AIS_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'AIS contribution to sea level (RCP45)')
  disint.rcp45 <- ncvar_def('AIS_DISINT_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'AIS fast dynamics contribution to sea level (RCP45)')
  lws.rcp45 <- ncvar_def('LWS_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'estimated LWS contribution to sea level (RCP45)')

  lsl.rcp85 <- ncvar_def('LocalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level (RCP85)')
  gsic.rcp85 <- ncvar_def('GSIC_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GSIC contribution to sea level (RCP85)')
  te.rcp85 <- ncvar_def('TE_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'TE contribution to sea level (RCP85)')
  gis.rcp85 <- ncvar_def('GIS_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GIS contribution to sea level (RCP85)')
  ais.rcp85 <- ncvar_def('AIS_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'AIS contribution to sea level (RCP85)')
  disint.rcp85 <- ncvar_def('AIS_DISINT_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'AIS fast dynamics contribution to sea level (RCP85)')
  lws.rcp85 <- ncvar_def('LWS_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'estimated LWS contribution to sea level (RCP85)')
}

gsl.hindcast <- ncvar_def('GlobalSeaLevel_hind', 'meters', list(dim.thind, dim.ensemble), -999,
                      longname = 'Global sea level (hindcast)')
temp.hindcast <- ncvar_def('temp_hind', 'deg C', list(dim.thind, dim.ensemble), -999,
                      longname = 'Global average surface temperature anomaly (hindcast)')
ocheat.hindcast <- ncvar_def('ocheat_hind', '10^22 J', list(dim.thind, dim.ensemble), -999,
                      longname = 'Ocean heat uptake (hindcast)')

gsl.rcp26 <- ncvar_def('GlobalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level (RCP26)')
temp.rcp26 <- ncvar_def('temp_RCP26', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP26)')
ocheat.rcp26 <- ncvar_def('ocheat_RCP26', '10^22 J', list(dim.tproj, dim.ensemble), -999,
                  longname = 'ocean heat uptake (RCP26)')

gsl.rcp45 <- ncvar_def('GlobalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level (RCP45)')
temp.rcp45 <- ncvar_def('temp_RCP45', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP45)')
ocheat.rcp45 <- ncvar_def('ocheat_RCP45', '10^22 J', list(dim.tproj, dim.ensemble), -999,
                  longname = 'ocean heat uptake (RCP45)')

gsl.rcp85 <- ncvar_def('GlobalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level (RCP85)')
temp.rcp85 <- ncvar_def('temp_RCP85', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP85)')
ocheat.rcp85 <- ncvar_def('ocheat_RCP85', '10^22 J', list(dim.tproj, dim.ensemble), -999,
                  longname = 'ocean heat uptake (RCP85)')

today=Sys.Date(); today=format(today,format="%d%b%Y")

if(experiment=='g') {
  outnc <- nc_create(filename.brickout,
                    list( gsl.rcp26, gsl.rcp45, gsl.rcp85,
                           temp.rcp26, ocheat.rcp26,
                           temp.rcp45, ocheat.rcp45,
                           temp.rcp85, ocheat.rcp85,
                          gsl.hindcast, temp.hindcast, ocheat.hindcast),
                    force_v4 = TRUE)
} else {
  outnc <- nc_create(filename.brickout,
                    list( gsl.rcp26, gsl.rcp45, gsl.rcp85, lsl.rcp26, lsl.rcp45, lsl.rcp85,
                          gsic.rcp26, te.rcp26, gis.rcp26, ais.rcp26, temp.rcp26, ocheat.rcp26, lws.rcp26, disint.rcp26,
                          gsic.rcp45, te.rcp45, gis.rcp45, ais.rcp45, temp.rcp45, ocheat.rcp45, lws.rcp45, disint.rcp45,
                          gsic.rcp85, te.rcp85, gis.rcp85, ais.rcp85, temp.rcp85, ocheat.rcp85, lws.rcp85, disint.rcp85,
                          gsl.hindcast, gsic.hindcast, te.hindcast, gis.hindcast, ais.hindcast, temp.hindcast, ocheat.hindcast,
                          ais.paleo.05, ais.paleo.50, ais.paleo.95, ais.paleo.max, ais.paleo.min,
                          ais.paleo.05.avg, ais.paleo.50.avg, ais.paleo.95.avg, ais.paleo.max.avg, ais.paleo.min.avg),
                    force_v4 = TRUE)
}

if(experiment!='g') {
  ncvar_put(outnc, lsl.rcp26, t(proj.rcp26$slr.nola))
  ncvar_put(outnc, gsic.rcp26, t(proj.rcp26$gsic))
  ncvar_put(outnc, te.rcp26, t(proj.rcp26$te))
  ncvar_put(outnc, gis.rcp26, t(proj.rcp26$gis))
  ncvar_put(outnc, ais.rcp26, t(proj.rcp26$ais))
  ncvar_put(outnc, disint.rcp26, t(proj.rcp26$disint))
  ncvar_put(outnc, lws.rcp26, t(proj.rcp26$lws))

  ncvar_put(outnc, lsl.rcp45, t(proj.rcp45$slr.nola))
  ncvar_put(outnc, gsic.rcp45, t(proj.rcp45$gsic))
  ncvar_put(outnc, te.rcp45, t(proj.rcp45$te))
  ncvar_put(outnc, gis.rcp45, t(proj.rcp45$gis))
  ncvar_put(outnc, ais.rcp45, t(proj.rcp45$ais))
  ncvar_put(outnc, disint.rcp45, t(proj.rcp45$disint))
  ncvar_put(outnc, lws.rcp45, t(proj.rcp45$lws))

  ncvar_put(outnc, lsl.rcp85, t(proj.rcp85$slr.nola))
  ncvar_put(outnc, gsic.rcp85, t(proj.rcp85$gsic))
  ncvar_put(outnc, te.rcp85, t(proj.rcp85$te))
  ncvar_put(outnc, gis.rcp85, t(proj.rcp85$gis))
  ncvar_put(outnc, ais.rcp85, t(proj.rcp85$ais))
  ncvar_put(outnc, disint.rcp85, t(proj.rcp85$disint))
  ncvar_put(outnc, lws.rcp85, t(proj.rcp85$lws))

  ncvar_put(outnc, gsic.hindcast, gsic.hind)
  ncvar_put(outnc, te.hindcast, te.hind)
  ncvar_put(outnc, gis.hindcast, gis.hind)
  ncvar_put(outnc, ais.hindcast, ais.hind)

  ncvar_put(outnc, ais.paleo.05, dais.paleo.05)
  ncvar_put(outnc, ais.paleo.50, dais.paleo.50)
  ncvar_put(outnc, ais.paleo.95, dais.paleo.95)
  ncvar_put(outnc, ais.paleo.max, dais.paleo.max)
  ncvar_put(outnc, ais.paleo.min, dais.paleo.min)
  ncvar_put(outnc, ais.paleo.05.avg, dais.paleo.05.avg)
  ncvar_put(outnc, ais.paleo.50.avg, dais.paleo.50.avg)
  ncvar_put(outnc, ais.paleo.95.avg, dais.paleo.95.avg)
  ncvar_put(outnc, ais.paleo.max.avg, dais.paleo.max.avg)
  ncvar_put(outnc, ais.paleo.min.avg, dais.paleo.min.avg)
}

ncvar_put(outnc, gsl.rcp26, t(proj.rcp26$slr))
ncvar_put(outnc, temp.rcp26, t(proj.rcp26$temp))
ncvar_put(outnc, ocheat.rcp26, t(proj.rcp26$ocheat))

ncvar_put(outnc, gsl.rcp45, t(proj.rcp45$slr))
ncvar_put(outnc, temp.rcp45, t(proj.rcp45$temp))
ncvar_put(outnc, ocheat.rcp45, t(proj.rcp45$ocheat))

ncvar_put(outnc, gsl.rcp85, t(proj.rcp85$slr))
ncvar_put(outnc, temp.rcp85, t(proj.rcp85$temp))
ncvar_put(outnc, ocheat.rcp85, t(proj.rcp85$ocheat))

ncvar_put(outnc, gsl.hindcast, gsl.hind)
ncvar_put(outnc, temp.hindcast, temp.hind)
ncvar_put(outnc, ocheat.hindcast, ocheat.hind)

nc_close(outnc)

##==============================================================================
##==============================================================================






if(l.dovandantzig) {
if(experiment=='c') {  # only do Van Dantzig analysis if control experiment
##==============================================================================
##==============================================================================
## Beginning here from a previous experiment?
if(FALSE){
  filename.in = "../output_model/BRICK-model_physical_control_01Nov2016.nc"
  ncdata <- nc_open(filename.in)
  sea_level = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
  mod.time =ncvar_get(ncdata, 'time_proj')
  nc_close(ncdata)
} else {
  sea_level=t(proj.rcp85$slr.nola)
}

## Evaluate flood risk analysis model for each of these realizations
source('../R/BRICK_VanDantzig.R')

vandantzig.ensemble = brick_vandantzig(sea_level = sea_level, time = mod.time)

## Save progress so far via workspace image
save.image(file=filename.saveprogress)

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## Output results of van Dantzig to netCDF file

library(ncdf4)

dim.heightening <- ncdim_def('H', 'meters', vandantzig.ensemble$Heightening[,1])
dim.ensemble <- ncdim_def('ens', 'ensemble member', (1:ncol(vandantzig.ensemble$Heightening)))

cost <- ncvar_def('ExpectedCost', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost')
loss <- ncvar_def('ExpectedLoss', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss')
investment <- ncvar_def('ExpectedInvestment', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment')
preturn <- ncvar_def('ExpectedPreturn', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period')

outnc <- nc_create(filename.vdout,
                    list(cost, loss, investment, preturn),
                    force_v4 = TRUE)

ncvar_put(outnc, cost, vandantzig.ensemble$Expected_costs)
ncvar_put(outnc, loss, vandantzig.ensemble$Expected_loss)
ncvar_put(outnc, investment, vandantzig.ensemble$Investment)
ncvar_put(outnc, preturn, 1/vandantzig.ensemble$Average_p_fail)

nc_close(outnc)

##==============================================================================
##==============================================================================
}  # end check if experiment=='c'
}  # end check if l.dovandantzig

t.end = proc.time()
time.minutes = (t.end-t.beg)[3]/60

print(paste('it took',time.minutes,'minutes to process an initial ensemble of',n.ensemble.report,'simulations'))



##==============================================================================
## End
##==============================================================================
