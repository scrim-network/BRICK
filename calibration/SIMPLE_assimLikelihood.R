##==============================================================================
# R Function:
# compute (log) likelihood for observations
# Inputs: Model parameter P (vector), choose method AR/VAR
# Outputs
#
#  -original file = "RaR.R"   Code written July 2014
#  - Author: Yawen Guan (yig5031@psu.edu)
#  - Modified: Tony Wong (twong@psu.edu)
#
#  -This function finds the log likelihood for a zero-mean AR1 process and the
#       simulated lag-1 autocorrelation coefficient for the model. Both outputs are
#       used in the MCMC likelihood programs as described in Ruckert et al. (GRL 2015). For
#       further description and references, please read the paper and the appendix.
#
#   -NOTE: The innovation variance = sigma^2 and the lag-1 autocorrelation
#       coefficient = rho1. In this function the initial values are ignored to
#       reduce bias and to make consistent with a VAR code. Description of a VAR
#       code can be found in the R package in review "VAR1"
#
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
## Estimate the log likelihood of the AR1 process
logl.ar1 = function(resid.in, sigma.in, rho.in, eps.in=0)
{
  n = length(resid.in)
  if(length(eps.in)==1) eps.in = rep(eps.in,n)

  #sigma.proc = sigma.in/sqrt(1-rho.in^2) # stationary process variance sigma.proc^2
  #logl = dnorm(resid.in[1],sd=sigma.proc+eps.in[1],log=TRUE)

  logl=0
  if(n>1) {
    resid.white = resid.in[2:n] - rho.in*resid.in[1:(n-1)] # this process whitens the residuals
    logl = logl + sum(dnorm(resid.white, sd=sqrt(sigma.in^2+eps.in[c(-1)]^2), log=TRUE)) # add in the sum of
  	          # density of the whitened residuals with a standard deviation of the
              # variance and the obs. errors
  }

  return(logl)
}
##==============================================================================
## log-likelihood of model output/data comparison, given the parameter values
log.lik = function( parameters.in,
                    parnames.in,
                    forcing.temp.in,
                    oidx,
                    midx,
                    obs,
                    obs.err,
                    ind.norm.data
                    )
{
	a.simple     = parameters.in[match("a.simple"    ,parnames.in)]
  b.simple     = parameters.in[match("b.simple"    ,parnames.in)]
  alpha.simple = parameters.in[match("alpha.simple",parnames.in)]
  beta.simple  = parameters.in[match("beta.simple" ,parnames.in)]
  V0           = parameters.in[match("V0"          ,parnames.in)]
  sigma.simple = parameters.in[match("sigma.simple",parnames.in)]
  rho.simple   = parameters.in[match("rho.simple"  ,parnames.in)]

  simple.out = simple(a=a.simple      , b=b.simple, alpha=alpha.simple ,
                       beta=beta.simple, V0=V0     , Tg=forcing.temp.in )

	## Subtract off normalization period
	itmp = ind.norm.data[match("gis",ind.norm.data[,1]),2]:ind.norm.data[match("gis",ind.norm.data[,1]),3]
	simple.out$sle.gis = simple.out$sle.gis - mean(simple.out$sle.gis[itmp])

  llik.simple  = 0

  resid.simple = obs$gis[oidx$gis] - simple.out$sle.gis[midx$gis]  #Calculate the residuals

  if(!all(is.finite(resid.simple))) {
    llik.simple = -Inf
  } else {
    llik.simple  = logl.ar1(resid.in=resid.simple, sigma.in=sigma.simple,
                            rho.in=rho.simple, eps.in=obs.err$gis) # AR(1) #Set up the likelihood function
  }

  llik = llik.simple # assume residuals are independent

  return(llik)
}
##==============================================================================

##==============================================================================
## (log) prior distributions of model parameters (all assumed to be uniform)
log.pri = function( parameters.in,
                    bound.lower.in,
                    bound.upper.in
                    )
{
  in.range = all(parameters.in > bound.lower.in) & all(parameters.in < bound.upper.in)

  if(in.range) {
    lpri=0
  } else {
    lpri = -Inf
  }

  return(lpri)
}
##==============================================================================

##==============================================================================
## (log) posterior distribution:  posterior ~ likelihood * prior
log.post = function(  parameters.in,
                      parnames.in,
                      bound.lower.in,
                      bound.upper.in,
                      forcing.temp.in,
                      midx,
                      oidx,
                      obs,
                      obs.err,
                      ind.norm.data
                      )
{
  lpri = log.pri( parameters.in=parameters.in,
                  bound.lower.in=bound.lower.in,
                  bound.upper.in=bound.upper.in
                  )
  if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
    lpost = log.lik(  parameters.in=parameters.in,
                      parnames.in=parnames.in,
                      forcing.temp.in=forcing.temp.in,
                      midx=midx,
                      oidx=oidx,
                      obs=obs,
                      obs.err=obs.err,
                      ind.norm.data=ind.norm.data
                      ) + lpri
  } else {
    lpost = -Inf
  }

  return(lpost)
}
##==============================================================================

##==============================================================================
## End
##==============================================================================
