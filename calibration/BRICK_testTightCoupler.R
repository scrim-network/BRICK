##==============================================================================
##	This script is for testing running ensembles of fully-forward coupled
##	simulations by using the results for the 1850 "initial condition" as the
##	uncertain initial condition in:
##		DOECLIM	--	T0, H0
##		GSIC		--	Gs0 (may need to adjust V0 higher)
##		SIMPLE	--	V0
##		TE			--	TE0
##		DAIS		--	n/a (initial ice sheet volume calculated from other parameters)
##
##	Questions? -- Tony Wong <twong@psu.edu
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

setwd('~/codes/BRICK/calibration')

library(ncdf4)
n.ensemble <- 10

# read calibrated model simulations
filename.brick <- '../output_model/BRICK-model_physical_control_01Nov2016.nc'
ncdata <- nc_open(filename.brick)
  t.hind <- ncvar_get(ncdata, 'time_hind')
  t.proj <- ncvar_get(ncdata, 'time_proj')
  gmsl.ctrl.hind <- ncvar_get(ncdata, 'GlobalSeaLevel_hind')
  temp.ctrl.hind <- ncvar_get(ncdata, 'temp_hind')
  ocheat.ctrl.hind <- ncvar_get(ncdata, 'ocheat_hind')
  ais.ctrl.hind <- ncvar_get(ncdata, 'AIS_hind')
  gis.ctrl.hind <- ncvar_get(ncdata, 'GIS_hind')
  te.ctrl.hind <- ncvar_get(ncdata, 'TE_hind')
	gsic.ctrl.hind <- ncvar_get(ncdata, 'GSIC_hind')
  gmsl.ctrl.rcp26 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
  gmsl.ctrl.rcp45 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
  gmsl.ctrl.rcp85 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
nc_close(ncdata)

#TODO

# draw ensemble parameters
filename.parameters <- '../output_calibration/BRICK-model_postcalibratedParameters_control_01Nov2016.nc'
ncdata <- nc_open(filename.parameters)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters <- t(parameters)
colnames(parameters) <- parnames

#TODO

# use for initial condition the 1850 hindcast (set up index i0 for this)
#TODO

# call the tightly coupled model
mod.time <- t.hind
l.project = FALSE

begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)

source('../R/compute_indices.R')				# function to determine the model and
source('../calibration/DOECLIM_readData.R')		# read DOECLIM calibration data
source('../calibration/GSIC_readData.R')			# read GSIC calibration data
source('../calibration/TE_readData.R')        # read TE data
source('../calibration/SIMPLE_readData.R')		# GIS data, and trends in mass balance
source('../calibration/DAIS_readData.R')			# DAIS forcing data (if at all uncoupled)

# initialize matrix to store model ensemble output
brick.out <- vector("list", n.ensemble)
source('../R/forcing_total.R')					# function to add up the total forcing
forcing = read.csv( '../data/forcing_hindcast.csv', header=TRUE )

# source the models
source('../R/BRICK_coupledModel.R')
source('../fortran/R/doeclimF.R')       # the DOECLIM model
source('../fortran/R/GSIC_magiccF.R')   # the GSIC model
source('../fortran/R/brick_te_F.R')     # TE (thermosteric expansion) model
source('../fortran/R/simpleF.R')        # GIS (Greenland Ice Sheet) model
source('../fortran/R/daisantoF.R')			# DAIS (Antarctic Ice Sheet) model

# set up the model
luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # thermosteric expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = TRUE    # Antarctic ice sheet model
luse.mymodel  = FALSE   # Example of adding your own model component
luse.brick = cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais)
source('../calibration/BRICK_parameterSetup.R')

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

## Run the sample, and enjoy a nice progress bar
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {

	brick.out[[i]] <- brick_model(parameters.in			= as.numeric(parameters[i,]),
																parnames.in				= parnames,
																forcing.in				= forcing,
																l.project					= l.project,
																slope.Ta2Tg.in		= slope.Ta2Tg,
																intercept.Ta2Tg.in= intercept.Ta2Tg,
																mod.time					= mod.time,
																ind.norm.data 		= ind.norm.data,
																ind.norm.sl 			= ind.norm,
																luse.brick 				= luse.brick,
																i0=i0
																)

  setTxtProgressBar(pb, i)
}
close(pb)

iH0 = match("H0",parnames)
iT0 = match("T0",parnames)
iGs0= match("Gs0",parnames)

par(mfrow=c(2,2))

plot(mod.time, brick.out[[1]]$doeclim.out$temp+parameters[1,iT0], type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$doeclim.out$temp+parameters[i,iT0], type='l', col='black')}
points(obs.temp.time, obs.temp, pch=16, col='red')

plot(mod.time, brick.out[[1]]$doeclim.out$ocheat+parameters[1,iH0], type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$doeclim.out$ocheat+parameters[i,iH0], type='l', col='black')}
points(obs.ocheat.time, obs.ocheat, pch=16, col='red')

plot(mod.time, brick.out[[1]]$gsic.out+parameters[1,iGs0], type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$gsic.out+parameters[i,iGs0], type='l', col='black')}
points(obs.gsic.time, obs.gsic, pch=16, col='red')

plot(mod.time, brick.out[[1]]$simple.out$sle.gis, type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$simple.out$sle.gis, type='l', col='black')}
points(obs.gis.time, obs.gis, pch=16, col='red')

##==============================================================================
## Okay, it all works. Now screw around with the "fully-forward" BRICK model.
## Want to code this as a two-step (run_BRICK_forward and BRICK_step_forward)
## R-calling-into-Fortran model, similar to the other BRICK sub-models.
##==============================================================================

##==============================================================================
## Test running just hte first ensemble member in "fully-forward" mode versus
## the "old" mode. Trouble-shooting time.
##==============================================================================






##==============================================================================
## End
##==============================================================================
