##==============================================================================
##
##  -file = "BRICK_assimLikelihood_AR1homo.R"   Origional code written July 2014
##  - Author: Yawen Guan (yig5031@psu.edu)
##  - Edited to run SLR model by: Kelsey Ruckert (klr324@psu.edu)
##  - Edited to run DOECLIM model by: Tony Wong (twong@psu.edu)
##  - Edited for BRICK by: Tony Wong (twong@psu.edu)
##
##  -This function computes the log likelihood for a zero-mean AR1 process from
##       observations as described in  Ruckert et al. (2016).
##       For further description and references, please read the paper
##       and the appendix.
##
##   -NOTE: Descriptions of how to use this for other observation and models
##       can be found in the R package in review "VAR1"
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

##==============================================================================
## AR1 model for errors (centered at 0 -- X(t)=rho*X(t-1)+eps(t), eps(t)~N(0,sigma1)
##		For clarity -- sigma.proc is for AR(1) process;
##					-- sigma1 (sampled) is the whitened innovation sigma
## Estimate the log likelihood of the AR1 process

if(TRUE){
## APPROX AR1?
logl.ar1 = function(r,sigma1,rho1,eps1=0) # default obs error is 0
{
  n = length(r) # r is the residuals
  if(length(eps1)==1) eps1 = rep(eps1,n)

	logl=0
	if(n>1) {
  	w = r[2:n] - rho1*r[1:(n-1)] # this process whitens the residuals
  	logl = logl + sum(dnorm(w,sd=sqrt((sigma1)^2+(eps1[c(-1)])^2),log=TRUE)) # add in the sum of
		    # density of the whitened residuals with a standard deviation of the
		    # variance and the obs. errors
  }
  return(logl)
}
##
} else {
## EXACT AR1?
library(mvtnorm)
logl.ar1 <-
  function(r,sigma1,rho1,eps1) # default obs error is 0. sigma1 is standard error.
  {
  library(mvtnorm)

    n = length(r)
    sigma.proc = sigma1/sqrt(1-rho1^2) # stationary process variance sigma.proc^2
    if(all(eps1==0)){
      logl = dnorm(r[1],sd=sigma.proc,log=TRUE)
      if(n>1) {
        w = r[2:n] - rho1*r[1:(n-1)] # whitened residuals
        # logl = logl + sum(dnorm(w,sd=sigma1,log=TRUE))
        # This is what we had before to make the computation faster.
        # This approximation should not change the result, but it is worth trying
        logl = logl + sum(dnorm(w,sd=sqrt(sigma1^2+eps1[-1]^2),log=TRUE))
      }
    }else{
      H <- abs(outer(1:n, 1:n, "-"))
      v = sigma.proc^2*rho1^H
      if(length(eps1)>1) {v = v+diag(eps1^2)
      } else {v = v+diag(rep(eps1^2,n))}
      # Need R package "mvtnorm"
      logl = dmvnorm(r,sigma=v,log=TRUE)
    }
    return(logl)
  }
##
}
##==============================================================================
## rest of the statistical model
##==============================================================================
log.lik = function( parameters.in,
                    parnames.in,
                    forcing.in,
                    l.project=FALSE,
                    rho.simple.in=NULL,
                    sigma.simple.in=NULL,
                    slope.Ta2Tg.in=1,
                    intercept.Ta2Tg.in=0,
                    mod.time,
                    ind.norm.data,
                    ind.norm.sl,
                    midx,
                    oidx,
                    obs,
                    obs.err,
                    trends.te,
                    luse.brick,
                    i0
                   ){

	## Run the coupled BRICK model
	brick.out = brick_model(  parameters.in=parameters.in,
                            parnames.in=parnames.in,
                            forcing.in=forcing.in,
                            l.project=l.project,
                            slope.Ta2Tg.in=slope.Ta2Tg.in,
                            intercept.Ta2Tg.in=intercept.Ta2Tg.in,
                            mod.time=mod.time,
                            ind.norm.data=ind.norm.data,
                            ind.norm.sl=ind.norm.sl,
                            luse.brick=luse.brick,
                            i0=i0
                          )

  ## Calculate contribution from DOECLIM temperature
	llik.temp = 0
  if(!is.null(oidx$temp) & luse.brick[,"luse.doeclim"]) {

    # Grab the DOECLIM statistical parameters
    T0        =parameters.in[match("T0"     ,parnames.in)]
  	sigma.T   =parameters.in[match("sigma.T",parnames.in)]
  	rho.T     =parameters.in[match("rho.T"  ,parnames.in)]

    # Calculate the DOECLIM temperature residuals; apply AR1 error model
    resid.temp= obs$temp[oidx$temp] - (brick.out$doeclim.out$temp[midx$temp]+T0)
    llik.temp = logl.ar1(resid.temp, sigma.T, rho.T, obs.err$temp[oidx$temp]) # AR(1)

  }

  ## Calculate contribution from DOECLIM ocean heat
  llik.ocheat = 0
  if(!is.null(oidx$ocheat) & luse.brick[,"luse.doeclim"]) {

    # Grab the DOECLIM statistical parameters
  	H0        =parameters.in[match("H0"     ,parnames.in)]
  	sigma.H   =parameters.in[match("sigma.H",parnames.in)]
  	rho.H     =parameters.in[match("rho.H"  ,parnames.in)]

    # Calculate the DOECLIM ocean heat residuals; apply AR1 error model
    resid.ocheat= obs$ocheat[oidx$ocheat] - (brick.out$doeclim.out$ocheat[midx$ocheat]+H0)
    llik.ocheat = logl.ar1(resid.ocheat, sigma.H, rho.H, obs.err$ocheat[oidx$ocheat]) # AR(1)

  }

  ## Calculate contribution from GSIC SLR
  llik.gsic = 0
  if(!is.null(oidx$gsic) & luse.brick[,"luse.gsic"]) {

    # Grab the GSIC statistical parameters
		sigma.gsic=parameters.in[match("sigma.gsic",parnames.in)]
		rho.gsic  =parameters.in[match("rho.gsic"  ,parnames.in)]

    # Calculate the GSIC residuals; apply AR1 error model
    resid.gsic= obs$gsic[oidx$gsic] - brick.out$gsic.out[midx$gsic]
    llik.gsic = logl.ar1(resid.gsic, sigma.gsic, rho.gsic, obs.err$gsic[oidx$gsic]) # AR(1)
  }

  ## Calculate contribution from thermosteric expansion
  llik.te = 0
  if(luse.brick[,"luse.te"]) {

    # Calculate the SLR residuals - only proceed if all TE SLR < total SLR
    # (all after the first 20 years, that is, because they are the 0 point)

    resid.sl.te = (obs$sl[oidx$sl]-mean(obs$sl[oidx$sl[1:20]])) -
                   (brick.out$te.out[midx$sl]-mean(brick.out$te.out[midx$sl[1:20]]))

    if(all(resid.sl.te[20:length(resid.sl.te)]>0)){

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
      llik.te = sum (dnorm(resid.trends, mean=rep(0,length(resid.trends)), sd = sqrt(err.trends^2), log=TRUE))

    } else {
      llik.te = -Inf
    }
  }

  ## Calculate contribution from SIMPLE (Greenland Ice Sheet)
  llik.simple = 0
  if(!is.null(oidx$gis) & luse.brick[,"luse.simple"]) {

    # Grab the SIMPLE statistical parameters
    sigma.simple=parameters.in[match("sigma.simple",parnames.in)]
    rho.simple  =parameters.in[match("rho.simple"  ,parnames.in)]

    # Overwrite the SIMPLE statistical parameters if values were fed into MCMC
    if(!is.null(rho.simple.in  )) rho.simple  =rho.simple.in
    if(!is.null(sigma.simple.in)) sigma.simple=sigma.simple.in

    # Calibrate SIMPLE based on GIS data alone?
    resid.simple = obs$gis[oidx$gis] - brick.out$simple.out$sle.gis[midx$gis] #Calculate the residuals

    if(!all(is.finite(resid.simple))) {
      llik.simple = -Inf
    } else {
      llik.simple  = logl.ar1(r=resid.simple, sigma1=sigma.simple,
                              rho1=rho.simple, eps1=obs.err$gis) # AR(1) #Set up the likelihood function
    }
  }

  ## Calculate contribution from Antarctic ice sheet
  llik.dais = 0
  if(luse.brick[,"luse.dais"]) {

    # Calculate the SLR residuals - only proceed if all AIS SLR < total SLR
    # (all after the first 20 years, that is, because they are the 0 point)

    resid.sl.ais = (obs$sl[oidx$sl]-mean(obs$sl[oidx$sl[1:20]])) -
                   (brick.out$dais.out[midx$sl]-mean(brick.out$dais.out[midx$sl[1:20]]))

    if(all(resid.sl.ais[20:length(resid.sl.ais)]>0)){
      lik.dais = 0
    } else {
      llik.dais = -Inf
    }
  }

  # Calculate contribution from total sea level rise
  llik.sl = 0
if(FALSE) {
  if(!is.null(oidx$sl)) {

    resid.sl = brick.out$slr.out[midx$sl] - obs$sl[oidx$sl]
    llik.sl  = sum (dnorm(resid.sl, mean=rep(0,length(resid.sl)), sd = sqrt(obs.err$sl[oidx$sl]^2), log=TRUE))

  }
}
	# Assume residual time series are independent
  llik = llik.temp + llik.ocheat + llik.gsic + llik.te + llik.simple + llik.dais + llik.sl

	return(llik)
}
##==============================================================================
## (log of the) prior probability
log.pri = function(parameters.in , parnames.in, bound.lower.in, bound.upper.in,
                    shape.in, scale.in )
{

	# Pluck off the model and statistical parameters (only ones without uniform prior)
  ind.S      = match("S",parnames.in)
  ind.invtau = match("invtau.te",parnames.in)
  ind.rho.gsic = match("rho.gsic",parnames.in)
  lpri.S      = 0
  lpri.invtau = 0

	in.range = all(parameters.in > bound.lower.in) & all(parameters.in < bound.upper.in)

	if(in.range){
		lpri.uni = 0									# Sum of all uniform priors (log(1)=0)
#    lpri.S = log(dcauchy(parameters.in[ind.S],location=3,scale=2) / 	# S has truncated Cauchy(3,2) prior
#					(pcauchy(bound.upper[ind.S],location=3,scale=2)-pcauchy(bound.lower[ind.S],location=3,scale=2)))
    if(!is.na(ind.invtau)) {lpri.invtau = dgamma( parameters.in[ind.invtau], shape=shape.in, scale=scale.in, log=TRUE)}
		lpri = lpri.uni + lpri.S + lpri.invtau
#    lpri = lpri.uni
	} else {
		lpri = -Inf
	}

	return(lpri)
}
##==============================================================================
## (log of the) posterior distribution:  posterior ~ likelihood * prior
log.post = function(  parameters.in,
                      parnames.in,
                      forcing.in,
                      bound.lower.in,
                      bound.upper.in,
                      l.project=FALSE,
                      rho.simple.in=NULL,
                      sigma.simple.in=NULL,
                      shape.in,
                      scale.in,
                      slope.Ta2Tg.in=1,
                      intercept.Ta2Tg.in=0,
                      mod.time,
                      ind.norm.data,
                      ind.norm.sl,
                      midx,
                      oidx,
                      obs,
                      obs.err,
                      trends.te,
                      luse.brick,
                      i0
                      ){

  llik = 0
	lpri = log.pri( parameters.in=parameters.in,
                  parnames.in=parnames.in,
                  bound.lower.in=bound.lower.in,
                  bound.upper.in=bound.upper.in,
                  shape.in=shape.in,
                  scale.in=scale.in )
  if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
  	llik = log.lik(   parameters.in=parameters.in,
                      parnames.in=parnames.in,
                      forcing.in=forcing.in,
                      l.project=l.project,
                      rho.simple.in=rho.simple.in,
                      sigma.simple.in=sigma.simple.in,
                      slope.Ta2Tg.in=slope.Ta2Tg.in,
                      intercept.Ta2Tg.in=intercept.Ta2Tg.in,
                      mod.time=mod.time,
                      ind.norm.data=ind.norm.data,
                      ind.norm.sl=ind.norm.sl,
                      midx=midx,
                      oidx=oidx,
                      obs=obs,
										  obs.err=obs.err,
											trends.te=trends.te,
                      luse.brick=luse.brick,
                      i0=i0
                      )
    lpost = llik + lpri
  } else {
  	lpost = -Inf
  }
  return(lpost)
}
##==============================================================================
## End
##==============================================================================
