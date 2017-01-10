##==============================================================================
## Define the coupled BRICK model and its inputs/outputs
## The inputs/outputs will vary based on which components you use, but an overview
## is given below.
## This version steps all models forward in time together
##
## Input:
##	parameters.in					input vector of model parameters
##	parnames.in						vector of parameter names
##	forcing.in						matrix of radiative forcing input
##	l.project							making projections or hindcasts?
##	slope.Ta2Tg.in				slope of Antarctic vs global mean temperature regression
##	intercept.Ta2Tg.in		intercept of Antarctic vs global mean temperature regression
##	mod.time							time (in years) of the model simulation
##	obs.temp							mean global surface temperature anomalies, used if no climate module
##	ind.norm.data					indices within the model output for setting zero anomaly of calibration data fields
##	ind.norm.sl						indices within model output for setting zero sea level
##	timestep							model timestep [years]
##	i0										index of reference year, within mod.time. For initial conditions to sub-models.
##
## Requires:
##  luse.brick, includes: luse.doeclim, luse.gsic, luse.te, luse.simple,
##                        luse.dais, and luse.XXX, where XXX
##                        may be replaced with your favorite model component
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
source('../R/BRICK_coupledModel_stepForward.R')

brick_model = function(
												parameters.in,
												parnames.in,
												forcing.in,
												l.project = FALSE,
												slope.Ta2Tg.in = 1,
												intercept.Ta2Tg.in = 0,
												mod.time,
												obs.temp = NULL,
												ind.norm.data = NULL,
												ind.norm.sl = NULL,
												luse.brick,
												timestep = 1,
												i0
												){

	# Initialize the list of output (do NOT grow lists/arrays in R)
	# The +1 is to have total global mean sea level (relative to ind.relative) in
	# the output.
	brick.out = vector('list',sum(luse.brick)+1)
	SL.out = rep(0,length(mod.time))
	outcnt=1

	#=============================================================================

	# initialize output arrays. all relative to 1850 (or wherever initial condition is)
	temp.out <- rep(NA,length(mod.time)); temp.out[1] <- 0

	l.fprint=TRUE

	# main time-stepping loop
	for (t in 1:(length(mod.time)-1)) {
		brick.new <- brick_model_stepforward(	parameters = parameters,
																					parnames = parnames,
																					forcing = forcing,
																					l.project = FALSE,
																					slope.Ta2Tg = slope.Ta2Tg,
																					intercept.Ta2Tg = intercept.Ta2Tg,
																					mod.time = mod.time[t:(t+1)],
																					luse.brick = luse.brick,
																					timestep = 1,
																					SL.old = SL.out[t],
																					l.fprint = l.fprint
																					)
		temp.out[t+1] <- brick.new$doeclim.out$temp[2]
	}

	# normalize?
	SL.couple = SL.couple - mean(SL.couple[ind.norm.sl])

	itmp = ind.norm.data[match("sl",ind.norm.data[,1]),2]:ind.norm.data[match("sl",ind.norm.data[,1]),3]
	SL.couple = slr.out
  SL.couple = SL.couple - mean(SL.couple[itmp])
	dSL.couple = c(-999,diff(slr.out))
		include_dSLais = 0		# in coupled model, feeding AIS dSL without AIS contribution

		## Check to make sure output from other models was reasonable
		if(any(is.na(SL.couple))) {
			slr.out <- rep(NA,length(mod.time))
			brick.out[[outcnt]] <- slr.out; names(brick.out)[outcnt] <- "dais.out"; outcnt <- outcnt+1;
		}	else {
			# run dais
		}

		## Set up the radiative forcing
		forcing.total <- forcing_total(forcing=forcing.in,
																	alpha.doeclim=alpha.doeclim,
																	l.project=l.project,
																	begyear=mod.time[1],
																	endyear=mod.time[length(mod.time)]
																	)

		## Run BRICK
#		brick.out <-


		## Normalize outputs to observations

	# Normalize temperature and ocean heat to match the observations
	itmp <- ind.norm.data[match("temp",ind.norm.data[,1]),2]:ind.norm.data[match("temp",ind.norm.data[,1]),3]
	doeclim.out$temp <- doeclim.out$temp - mean(doeclim.out$temp[itmp])

	itmp <- ind.norm.data[match("gsic",ind.norm.data[,1]),2]:ind.norm.data[match("gsic",ind.norm.data[,1]),3]
	gsic.out.norm <- gsic.out - mean(gsic.out[itmp])

	itmp <- ind.norm.data[match("te",ind.norm.data[,1]),2]:ind.norm.data[match("te",ind.norm.data[,1]),3]
	te.out.norm <- te.out - mean(te.out[itmp])

	itmp = ind.norm.data[match("gis",ind.norm.data[,1]),2]:ind.norm.data[match("gis",ind.norm.data[,1]),3]
	simple.out$sle.gis = simple.out$sle.gis - mean(simple.out$sle.gis[itmp])

	#=============================================================================
	# Total sea-level rise

	## Add the SLR to the output
	brick.out[[sum(luse.brick)+1]] = slr.out - mean(slr.out[ind.norm.sl])
  names(brick.out)[sum(luse.brick)+1]="slr.out"

	## Check to make sure all the output made it
	if(outcnt!=sum(luse.brick)+1) print('ERROR - missing model output!')

	#=============================================================================

	return(brick.out)
}



##==============================================================================
## End
##==============================================================================
