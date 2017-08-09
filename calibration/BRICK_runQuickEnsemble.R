##==============================================================================
##	This script is for running a quick ensemble of hindcast simulations drawn
##	from the second half of an MCMC calibration chain.
##
##	400 is the default ensemble size. This is so the ensemble is large enough to
##	resolve the 0.5% and 99.5% quantiles. Obviously, more members is better, but
##	the purpose here is a sanity check, not science results.
##
##	This script is set up to be run interactively within the BRICK_calib_driver.R
##	script.
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

n.ensemble = 400

## Draw ensemble from the second have of "chain1" (which we assume we have from
## the BRICK_calib_driver.R script)

parameters.ensemble = mat.or.vec(n.ensemble , ncol(chain1) )
ind = sample( seq( round(0.5*nrow(chain1)),nrow(chain1) ), size=n.ensemble, replace=FALSE)
for (p in 1:ncol(chain1)){
	for (i in 1:n.ensemble){
		parameters.ensemble[i,p] = chain1[ind[i],p]
	}
}


## Initialize matrix to store model ensemble output
brick.out = vector("list", n.ensemble)

## Run the sample, and enjoy a nice progress bar
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {

	brick.out[[i]] = brick_model(	parameters.in			= as.numeric(parameters.ensemble[i,]),
																parnames.in				= parnames,
																forcing.in				= forcing,
																l.project					= l.project,
																slope.Ta2Tg.in		= slope.Ta2Tg,
																intercept.Ta2Tg.in= intercept.Ta2Tg,
																mod.time					= mod.time,
																ind.norm.data 		= ind.norm.data,
																ind.norm.sl 			= ind.norm,
																luse.brick 				= luse.brick,
																i0                = i0,
                                l.aisfastdy       = l.aisfastdy)

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
