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
##	modified for BRICK v0.2 (fully-forward version) Tony Wong, 15 Feb 2017
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
	forcing.raw,
	l.project=FALSE,
	tstep=1,
	slope.Ta2Tg=1,
	intercept.Ta2Tg=0,
	mod.time,
	obs.temp=NULL,
	ind.norm.data=NULL,
	ind.norm.sl=NULL,
	midx,
	oidx,
	obs,
	obs.err,
	trends.te,
	trends.ais,
	luse.brick
){

    brick.out <- brickF(tstep=tstep,
                        mod.time=mod.time,
                        forcing.raw = forcing.raw,
						l.project = l.project,
                        S.doeclim = parameters.in[match("S.doeclim",parnames.in)],
                        kappa.doeclim = parameters.in[match("kappa.doeclim",parnames.in)],
						alpha.doeclim = parameters.in[match("alpha.doeclim",parnames.in)],
                        T0.doeclim = parameters.in[match("T0.doeclim",parnames.in)],
                        H0.doeclim = parameters.in[match("H0.doeclim",parnames.in)],
                        beta0.gsic = parameters.in[match("beta0.gsic",parnames.in)],
						V0.gsic = parameters.in[match("V0.gsic",parnames.in)],
                        n.gsic = parameters.in[match("n.gsic",parnames.in)],
                        Gs0.gsic = parameters.in[match("Gs0.gsic",parnames.in)],
                        a.simple = parameters.in[match("a.simple",parnames.in)],
                        b.simple = parameters.in[match("b.simple",parnames.in)],
                        alpha.simple = parameters.in[match("alpha.simple",parnames.in)],
                        beta.simple = parameters.in[match("beta.simple",parnames.in)],
                        V0.simple = parameters.in[match("V0.simple",parnames.in)],
                        a.te = parameters.in[match("a.te",parnames.in)],
                        b.te = parameters.in[match("b.te",parnames.in)],
                        invtau.te = parameters.in[match("invtau.te",parnames.in)],
                        V0.te = parameters.in[match("V0.te",parnames.in)],
                        a.anto = parameters.in[match("anto.a",parnames.in)],
                        b.anto = parameters.in[match("anto.b",parnames.in)],
                        slope.Ta2Tg = slope.Ta2Tg,
                        intercept.Ta2Tg = intercept.Ta2Tg,
                        b0.dais = parameters.in[match("b0",parnames.in)],
                        slope.dais = parameters.in[match("slope",parnames.in)],
                        mu.dais = parameters.in[match("mu",parnames.in)],
                        h0.dais = parameters.in[match("h0",parnames.in)],
                        c.dais = parameters.in[match("c",parnames.in)],
                        P0.dais = parameters.in[match("P0",parnames.in)],
                        kappa.dais = parameters.in[match("kappa.dais",parnames.in)],
                        nu.dais = parameters.in[match("nu",parnames.in)],
                        f0.dais = parameters.in[match("f0",parnames.in)],
                        gamma.dais = parameters.in[match("gamma",parnames.in)],
                        alpha.dais = parameters.in[match("alpha.dais",parnames.in)]
                        )

	doeclim.norm.resid <- 0
	gsic.norm.resid <- 0
	simple.norm.resid <- 0
	te.norm.resid <- 0
	ais.norm.resid <- 0
	sl.norm.resid <- 0

	# Temperature/ocean heat uptake contributions
	if(luse.brick[,"luse.doeclim"]) {
  		if(!is.null(oidx$temp)) {
			itmp <- ind.norm.data[which(ind.norm.data[,1]=='temp'),2]:ind.norm.data[which(ind.norm.data[,1]=='temp'),3]
			temperature.model <- brick.out$temp_out - mean(brick.out$temp_out[itmp])
			doeclim.norm.resid <- doeclim.norm.resid +
						mean(abs( (obs$temp[oidx$temp]-(temperature.model[midx$temp])        )/
    											obs.err$temp[oidx$temp]     ))
		}
		if(!is.null(oidx$ocheat)) {
			# ocean heat uptake does not need to be normalized
			doeclim.norm.resid <- doeclim.norm.resid +
						mean(abs( (obs$ocheat[oidx$ocheat]-(brick.out$ocheat[midx$ocheat]))/
											  	obs.err$ocheat[oidx$ocheat] ))
		}
	}

	# GSIC contribution
	if(!is.null(oidx$gsic) & luse.brick[,"luse.gsic"]) {
		itmp <- ind.norm.data[which(ind.norm.data[,1]=='gsic'),2]:ind.norm.data[which(ind.norm.data[,1]=='gsic'),3]
		gsic.model <- brick.out$sl_gsic_out - mean(brick.out$sl_gsic_out[itmp])
		gsic.norm.resid = mean(abs( (obs$gsic[oidx$gsic]-gsic.model[midx$gsic])/obs.err$gsic[oidx$gsic] ))
	}

	# TE contribution
	if(luse.brick[,"luse.te"]) {
    	# Note 1: the trends from IPCC are in mm/year, and model output is m
    	# Note 2: these calculate the least squares regression slope coefficients. It
    	# is more than twice as fast to calcualte by hand like this than to use R's
    	# "lm(...)" function.
    	# Note 3: Need 1000*trends.mod because they're in meters, but trends.te is mm
		# Note 4: No normalization needed, because these are trends (slopes)

		trends.mod = rep(0, nrow(trends.te))
		for (i in 1:nrow(trends.te)) {
			x = seq(trends.te[i,6],trends.te[i,7]);                 barx = mean(x);
			y = brick.out$sl_te_out[trends.te[i,6]:trends.te[i,7]]; bary = mean(y);
	  		trends.mod[i] = sum( (x-rep(barx,length(x)))*(y-rep(bary,length(y))))/sum( (x-rep(barx,length(x)))^2 )
		}
	  	resid.trends = 1000*trends.mod - trends.te[,1]
		err.trends   = 0.5*(trends.te[,3]-trends.te[,2])
		te.norm.resid = mean(abs(resid.trends)/ err.trends )
	}

	# GIS contribution
  	if(!is.null(oidx$gis) & luse.brick[,"luse.simple"]) {
		# observations are relative to 1960-1990 average (SIMPLE_readData.R)
		itmp <- ind.norm.data[which(ind.norm.data[,1]=='gis'),2]:ind.norm.data[which(ind.norm.data[,1]=='gis'),3]
		gis.model <- brick.out$sl_gis_out - mean(brick.out$sl_gis_out[itmp])
		simple.resid  = (obs$gis[oidx$gis] - gis.model[midx$gis])/obs.err$gis
		if(!all(is.finite(simple.resid))) {
      		simple.norm.resid = Inf
    	} else {
			simple.norm.resid = mean(abs( simple.resid))
		}

	}

	# AIS contribution
	if(!is.null(oidx$ais) & luse.brick[,"luse.dais"]) {
		# First part is from Shepherd et al 2012 instrumental point (from
		# Ruckert et al 2017, or Wong et al 2017)
		itmp <- ind.norm.data[which(ind.norm.data[,1]=='ais'),2]:ind.norm.data[which(ind.norm.data[,1]=='ais'),3]
		ais.model <- brick.out$sl_ais_out - mean(brick.out$sl_ais_out[itmp])
		ais.instr.resid <- abs(ais.model[midx$ais] - obs$ais[oidx$ais])/obs.err$ais[oidx$ais]
		ais.norm.resid <- ais.norm.resid + ais.instr.resid

		# Second part is from IPCC trends
		trends.mod = rep(0, nrow(trends.ais))
		for (i in 1:nrow(trends.ais)) {
			x = seq(trends.ais[i,6],trends.ais[i,7]);                  barx = mean(x);
			y = brick.out$sl_ais_out[trends.ais[i,6]:trends.ais[i,7]]; bary = mean(y);
	  		trends.mod[i] = sum( (x-rep(barx,length(x)))*(y-rep(bary,length(y))))/sum( (x-rep(barx,length(x)))^2 )
		}
	  	resid.trends = 1000*trends.mod - trends.ais[,1]
		err.trends   = 0.5*(trends.ais[,3]-trends.ais[,2])
		ais.trend.resid = mean(abs(resid.trends)/ err.trends )

		ais.norm.resid <- ais.norm.resid + ais.trend.resid
	}

	# GMSL contribution
	sl.model <- brick.out$sl_out - mean(brick.out$sl_out[ind.norm.sl])
	resid.sl_lw.1900 <- mean( abs(obs$sl_lw$r1900$sl - sl.model[obs$sl_lw$r1900$midx]) / obs$sl_lw$r1900$err )
	resid.sl_lw.1970 <- mean( abs(obs$sl_lw$r1970$sl - sl.model[obs$sl_lw$r1970$midx]) / obs$sl_lw$r1970$err )
	resid.sl_lw.1992 <- mean( abs(obs$sl_lw$r1992$sl - sl.model[obs$sl_lw$r1992$midx]) / obs$sl_lw$r1992$err )

	sl.norm.resid = sl.norm.resid + resid.sl_lw.1900 + resid.sl_lw.1970 + resid.sl_lw.1992

	# Tally up the total normalized residual.
	err.sum = doeclim.norm.resid + gsic.norm.resid + simple.norm.resid +
						te.norm.resid      + sl.norm.resid

	if(is.nan(err.sum)) err.sum=Inf
	return(err.sum)

}

##==============================================================================
## End
##==============================================================================
