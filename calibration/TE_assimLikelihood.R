##==============================================================================
## Likelihood functions for (pre)calibrating thermosteric expansion model in
## BRICK.
## The data used are:
##    (1) The IPCC AR5 Ch13 1971-2010 trend in TE sea level rise: 0.8+/-0.3 mm/year
##    (2) The IPCC AR5 Ch13 1993-2010 trend in TE sea level rise: 1.1+/-0.3 mm/year
##    (3) Step function likelihood, for whether or not modeled TE SLR exceeds
##        the total SLR (Church and White (2011))
##
## Questions? Tony Wong (twong@psu.edu)
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
## Log-likelihood
log.lik = function( parameters.in,
                    parnames.in,
                    forcing.temp.in,
                    oidx,
                    midx,
                    obs,
                    obs.err,
                    trends.te,
                    ind.norm
                    )
{

  a.te     =parameters.in[match("a.te"     ,parnames.in)]
  b.te     =parameters.in[match("b.te"     ,parnames.in)]
  invtau.te=parameters.in[match("invtau.te",parnames.in)]
  TE0      =parameters.in[match("TE0"      ,parnames.in)]

  te.out = brick_te_F(a=a.te , b=b.te , invtau=invtau.te, TE_0=TE0, Tg=forcing.temp.in)

  # Normalize to the 1961-1990 mean? Or first 20 years? (1880-1900)
  # Use first 20 years and only check residuals after that for exceeding total SLR
  # because there will be noise around the trend in the first 20 years; don't want
  # to throw out the run because of a little noise.

  te.out = te.out - mean(te.out[ind.norm])

  # Calculate the SLR residuals - only proceed if all TE SLR < total SLR
  # (all after the first 20 years, that is, because they are the 0 point)
  resid.sl = (obs$sl[oidx$sl]-mean(obs$sl[oidx$sl[1:20]])) -
             (te.out[midx$sl]-mean(te.out[midx$sl[1:20]]))

  if(all(resid.sl[20:length(resid.sl)]>0)){

    # Note 1: the trends from IPCC are in mm/year, and model output is m
    # Note 2: these calculate the least squares regression slope coefficients. It
    # is more than twice as fast to calcualte by hand like this than to use R's
    # "lm(...)" function.

		trends.mod = rep(0, nrow(trends.te))
		for (i in 1:nrow(trends.te)) {
			x = seq(trends.te[i,6],trends.te[i,7]);    barx = mean(x);
			y = te.out[trends.te[i,6]:trends.te[i,7]]; bary = mean(y);
	  	trends.mod[i] = sum( (x-rep(barx,length(x)))*(y-rep(bary,length(y))))/sum( (x-rep(barx,length(x)))^2 )
		}
	  resid.trends = 1000*trends.mod - trends.te[,1]
		err.trends   = 0.5*(trends.te[,3]-trends.te[,2])
    llik = sum (dnorm(resid.trends, mean=rep(0,length(resid.trends)), sd = sqrt(err.trends^2), log=TRUE))

  } else {
    llik = -Inf
  }

  llik
}
##==============================================================================
## Log-prior
log.pri = function( parameters.in,
                    parnames.in,
                    bound.lower.in,
                    bound.upper.in,
                    shape.in,
                    scale.in
                    )
{

  invtau = parameters.in[match("invtau.te",parnames.in)]

  in.range = all(parameters.in > bound.lower.in) & all(parameters.in < bound.upper.in)
  lpri.invtau = 0

  if(in.range) {
    # Assign a gamma prior to 1/tau
    lpri.invtau = dgamma( invtau, shape=shape.in, scale=scale.in, log=TRUE)
    lpri = lpri.invtau
  } else {
    lpri = -Inf
  }

  lpri
}
##==============================================================================
## Log-posterior distribution:  posterior ~ likelihood * prior
log.post = function(  parameters.in,
                      parnames.in,
                      forcing.temp.in,
                      bound.lower.in,
                      bound.upper.in,
                      oidx,
                      midx,
                      obs,
                      obs.err,
                      trends.te,
                      shape.in,
                      scale.in,
                      ind.norm
                      )
{
  lpri = log.pri( parameters.in,
                  parnames.in=parnames.in,
                  bound.lower.in=bound.lower.in,
                  bound.upper.in=bound.upper.in,
                  shape.in=shape.in,
                  scale.in=scale.in
                  )
  if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
    llik = log.lik(  parameters.in=parameters.in,
                     parnames.in=parnames.in,
                     forcing.temp.in=forcing.temp.in,
                     oidx=oidx,
                     midx=midx,
                     obs=obs,
                     obs.err=obs.err,
                     trends.te=trends.te,
                     ind.norm
                     )
    lpost = llik + lpri
  } else {
    lpost = -Inf
  }
  lpost
}
##==============================================================================
## End
##==============================================================================
