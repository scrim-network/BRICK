##==============================================================================
#
#  -original file = "RaR.R"   Code written July 2014
#  - Author: Yawen Guan (yig5031@psu.edu)
#
#  -This function finds the log likelihood for a zero-mean AR1 process and the
#       simulated lag-1 autocorrelation coefficient for the model. Both outputs are
#       used in the MCMC likelihood programs as described in Ruckert et al. (2016). For
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

# Estimate the log likelihood of the AR1 process
logl.ar1 = function(r,sigma1,rho1,eps1=0) # default obs error is 0
{
  n = length(r) # r is the residuals

  #sigma.proc = sigma1/sqrt(1-rho1^2) # stationary process variance sigma.proc^2
  #logl = dnorm(r[1],sd=sigma.proc+eps1[1],log=TRUE)

  logl=0
  if(n>1) {
    w = r[2:n] - rho1*r[1:(n-1)] # this process whitens the residuals
    logl = logl + sum(dnorm(w,sd=sqrt((sigma1)^2+(eps1[c(-1)])^2),log=TRUE)) # add in the sum of
           # density of the whitened residuals with a standard deviation of the
           # variance and the obs. errors
  }
  return(logl)
}

##==============================================================================
#
#  -original file = "WRobs_likelihood_AR.R"   Code written July 2014
#  - Author: Yawen Guan (yig5031@psu.edu)
#  - Edited to run SLR model by: Kelsey Ruckert (klr324@psu.edu)
#
#  -This function computes the log likelihood for a zero-mean AR1 process from
#       observations as described in  Ruckert et al. (2016).
#       For further description and references, please read the paper
#       and the appendix.
#
#   -NOTE: Descriptions of how to use this for other observation and models
#       can be found in the R package in review "VAR1"
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
##==============================================================================
log.lik = function(	parameters.in,
										parnames.in,
										forcing.temp.in,
										oidx,
										midx,
										obs,
										obs.err,
                    ind.norm.data
									  ){

	beta0  = parameters.in[match("beta0"  ,parnames.in)]
	V0.gsic= parameters.in[match("V0.gsic",parnames.in)]
	n      = parameters.in[match("n"      ,parnames.in)]
	Gs0    = parameters.in[match("Gs0"    ,parnames.in)]
	sigma.gsic=parameters.in[match("sigma.gsic",parnames.in)]
	rho.gsic  =parameters.in[match("rho.gsic"  ,parnames.in)]

  #Calculate a model simulation based on a randomly sampled parameter set
  gsic.out = gsic_magiccF(beta0=beta0, V0=V0.gsic, n=n, Gs0=Gs0, Tg=forcing.temp.in)

	## Subtract off normalization period model GSIC output as the zero point
	itmp = ind.norm.data[match("gsic",ind.norm.data[,1]),2]:ind.norm.data[match("gsic",ind.norm.data[,1]),3]
	gsic.out.norm = gsic.out - mean(gsic.out[itmp])

  llik.gsic  = 0
  resid.gsic = obs$gsic[oidx$gsic] - gsic.out.norm[midx$gsic] #Calculate the residuals
  llik.y  = logl.ar1(resid.gsic, sigma.gsic, rho.gsic, obs.err$gsic) # AR(1) #Set up the likelihood function

  llik = llik.y # assume residuals are independent

  return(llik)
}
##==============================================================================
log.pri = function( parameters.in,
                    bound.lower.in,
                    bound.upper.in
                    )
{

  in.range = all(parameters.in >= bound.lower.in) & all(parameters.in <= bound.upper.in)

  if(in.range) {
    lpri=0
  } else {
    lpri = -Inf
  }

  return(lpri)
}
##==============================================================================
# (log) posterior distribution:  posterior ~ likelihood * prior
log.post = function(  parameters.in,
                      parnames.in,
                      bound.lower.in,
                      bound.upper.in,
                      forcing.temp.in,
                      oidx,
                      midx,
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
                      oidx=oidx,
                      midx=midx,
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
## End
##==============================================================================
