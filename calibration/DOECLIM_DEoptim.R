##==============================================================================
##  File = "DOECLIM_DEoptim.R"
##
##  Author: Tony Wong <twong@psu.edu>
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

minimize_residuals_doeclim = function(	parameters.in,
																				parnames.in,
																				forcing.in,
																				mod.time,
																				l.project=FALSE,
																				oidx,
																				midx,
																				obs,
																				obs.err,
																				ind.norm.data
																				){

  S            =parameters.in[match("S"            ,parnames.in)]
	kappa.doeclim=parameters.in[match("kappa.doeclim",parnames.in)]
	alpha.doeclim=parameters.in[match("alpha.doeclim",parnames.in)]
	T0           =parameters.in[match("T0"           ,parnames.in)]
	H0           =parameters.in[match("H0"           ,parnames.in)]

	## Set up the radiative forcing
	forcing.total = forcing_total(  forcing  =forcing.in , alpha.doeclim =alpha.doeclim,
																	l.project=l.project  , begyear       =mod.time[1]  ,
																	endyear  =mod.time[length(mod.time)])

	model.out = doeclimF( S=S,
												kappa=kappa.doeclim,
												forcing.total=forcing.total,
												mod.time=mod.time
												)

	## Normalize temperature and ocean heat to match the observations
	itmp = ind.norm.data[match("temp",ind.norm.data[,1]),2]:ind.norm.data[match("temp",ind.norm.data[,1]),3]
	model.out$temp = model.out$temp - mean(model.out$temp[itmp])

	#itmp = ind.norm.data[match("ocheat",ind.norm.data[,1]),2]:ind.norm.data[match("ocheat",ind.norm.data[,1]),3]
	#model.out$ocheat = model.out$ocheat - mean(model.out$ocheat[itmp])

  err.sum = sum(abs(obs$temp[oidx$temp]     - (model.out$temp[midx$temp]    +T0))/
    		obs.err$temp[oidx$temp]) +
						sum(abs(obs$ocheat[oidx$ocheat] - (model.out$ocheat[midx$ocheat]+H0))/
				obs.err$ocheat[oidx$ocheat])

	if(is.nan(err.sum)) err.sum=Inf
	return(err.sum)

}

##==============================================================================
## End
##==============================================================================
