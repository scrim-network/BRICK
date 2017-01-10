##==============================================================================
##  File = "BRICK_DEoptim.R"
##
## Input:
##	parameters.in     parameters over which to optimize the coupled BRICK model
##  parnames.in    		standard names of the parameters (see README file)
##  l.project         making projections using RCP8.5 forcing (TRUE) or historical
##                    forcing data (FALSE)?
##
##
## Notes:
##  -- Minimizing on the mean absolute residuals, normalized by the observational
##     standard errors. If you minimize over the sum of the absolute residuals,
##     you are biased towards the temperature data, because there is much more of
##     it.
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

minimize_residuals_brick = function(
																		parameters.in,
																		parnames.in,
																		forcing.in,
																		l.project=FALSE,
																		slope.Ta2Tg.in=1,
																		intercept.Ta2Tg.in=0,
																		mod.time,
																		obs.temp=NULL,
																		ind.norm.data=NULL,
																		ind.norm.sl=NULL,
																		midx,
																		oidx,
																		obs,
																		obs.err,
																		trends.te,
																		luse.brick,
																		i0
																		){

	T0 = parameters.in[match("T0",parnames.in)]
	H0 = parameters.in[match("H0",parnames.in)]

	brick.out = brick_model(parameters.in,
													parnames.in       =parnames.in,
													l.project         =l.project,
													forcing.in        =forcing.in,
													slope.Ta2Tg.in    =slope.Ta2Tg.in,
													intercept.Ta2Tg.in=intercept.Ta2Tg.in,
													mod.time          =mod.time,
													ind.norm.data     =ind.norm.data,
													ind.norm.sl       =ind.norm.sl,
													luse.brick        =luse.brick,
													i0								=i0
													)

	doeclim.norm.resid=0
	gsic.norm.resid=0
	simple.norm.resid=0
	te.norm.resid=0
	sl.norm.resid=0

	if(luse.brick[,"luse.doeclim"]) {
  	if(!is.null(oidx$temp)) {
			doeclim.norm.resid = doeclim.norm.resid +
						mean(abs( (obs$temp[oidx$temp]-(brick.out$doeclim.out$temp[midx$temp]+T0)        )/
    											obs.err$temp[oidx$temp]     ))
		}
		if(!is.null(oidx$ocheat)) {
			doeclim.norm.resid = doeclim.norm.resid +
						mean(abs( (obs$ocheat[oidx$ocheat]-(brick.out$doeclim.out$ocheat[midx$ocheat]+H0))/
											  	obs.err$ocheat[oidx$ocheat] ))
		}
	}

  if(!is.null(oidx$gsic) & luse.brick[,"luse.gsic"]) {
		gsic.norm.resid = mean(abs( (obs$gsic[oidx$gsic]-brick.out$gsic.out[midx$gsic])/obs.err$gsic[oidx$gsic] ))
	}

  if(luse.brick[,"luse.te"]) {
    # Note 1: the trends from IPCC are in mm/year, and model output is m
    # Note 2: these calculate the least squares regression slope coefficients. It
    # is more than twice as fast to calcualte by hand like this than to use R's
    # "lm(...)" function.
    # Note 3: Need 1000*trends.mod because they're in meters, but trends.te is mm
		trends.mod = rep(0, nrow(trends.te))
		for (i in 1:nrow(trends.te)) {
			x = seq(trends.te[i,6],trends.te[i,7]);              barx = mean(x);
			y = brick.out$te.out[trends.te[i,6]:trends.te[i,7]]; bary = mean(y);
	  	trends.mod[i] = sum( (x-rep(barx,length(x)))*(y-rep(bary,length(y))))/sum( (x-rep(barx,length(x)))^2 )
		}
	  resid.trends = 1000*trends.mod - trends.te[,1]
		err.trends   = 0.5*(trends.te[,3]-trends.te[,2])
		te.norm.resid = mean(abs(resid.trends/ err.trends ))
	}

  if(!is.null(oidx$gis) & luse.brick[,"luse.simple"]) {
	  simple.resid  = (obs$gis[oidx$gis]-brick.out$simple.out$sle.gis[midx$gis])/obs.err$gis
  	if(!all(is.finite(simple.resid))) {
      simple.norm.resid = Inf
    } else {
			simple.norm.resid = mean(abs( simple.resid))
		}
	}

	#sl.norm.resid = sl.norm.resid + mean(abs( (brick.out$slr.out[midx$sl] - obs$sl[oidx$sl])/obs.err$sl[oidx$sl] ))

	err.sum = doeclim.norm.resid + gsic.norm.resid + simple.norm.resid +
						te.norm.resid      + sl.norm.resid

	if(is.nan(err.sum)) err.sum=Inf
	return(err.sum)

}

##==============================================================================
## End
##==============================================================================
