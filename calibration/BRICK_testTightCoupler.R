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

TODO

# draw ensemble parameters
filename.parameters <- '../output_calibration/BRICK-model_postcalibratedParameters_control_01Nov2016.nc'
ncdata <- nc_open(filename.parameters)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters <- t(parameters)
colnames(parameters) <- parnames

TODO

# use for initial condition the 1850 hindcast (set up index i0 for this)
TODO

# call the tightly coupled model
mod.time <- t.hind



# initialize matrix to store model ensemble output
brick.out <- vector("list", n.ensemble)

## Run the sample, and enjoy a nice progress bar
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {

	brick.out[[i]] <- brick_model(parameters.in			= as.numeric(parameters.ensemble[i,]),
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

plot(mod.time, brick.out[[1]]$doeclim.out$temp+parameters.ensemble[1,iT0], type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$doeclim.out$temp+parameters.ensemble[i,iT0], type='l', col='black')}
points(obs.temp.time, obs.temp, pch=16, col='red')

plot(mod.time, brick.out[[1]]$doeclim.out$ocheat+parameters.ensemble[1,iH0], type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$doeclim.out$ocheat+parameters.ensemble[i,iH0], type='l', col='black')}
points(obs.ocheat.time, obs.ocheat, pch=16, col='red')

plot(mod.time, brick.out[[1]]$gsic.out+parameters.ensemble[1,iGs0], type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$gsic.out+parameters.ensemble[i,iGs0], type='l', col='black')}
points(obs.gsic.time, obs.gsic, pch=16, col='red')

plot(mod.time, brick.out[[1]]$simple.out$sle.gis, type='l', col='black')
for (i in 1:n.ensemble) {lines(mod.time, brick.out[[i]]$simple.out$sle.gis, type='l', col='black')}
points(obs.gis.time, obs.gis, pch=16, col='red')

##==============================================================================
## End
##==============================================================================
