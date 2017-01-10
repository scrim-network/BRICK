##==============================================================================
## R Function: GMSL_assimLikelihood.R
## -Calibration of GMSL model (Rahmstorf 2007)
##
## compute (log) likelihood for observations
## The observations are independent and identically distributed (IID)
##
## -Original authors: Yawen Guan and Kelsey Ruckert (klr324@psu.edu)
## - Modified by Tony Wong (twong@psu.edu) for running within BRICK model framework
## - Modified by Tony Wong (twong@psu.edu) for running GMSL model
##==============================================================================
## -June 10 2015 (original)
## - 6 July 2016 (modified)
## - 19 October 2016 (modified)
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
                    obs.in,
                    obs.err.in,
                    ind.norm.in,
                    Tg.in=NULL
                    )
{

  a.gmsl=parameters.in[match("a.gmsl",parnames.in)]
  Teq.gmsl=parameters.in[match("Teq.gmsl",parnames.in)]

  gmsl.out = gmsl_r07(
                       a = a.gmsl,
                       Teq = Teq.gmsl,
                       Tg = Tg.in
                       )

  gmsl.out.norm = gmsl.out - mean(gmsl.out[ind.norm.in])

  llik0  = 0

  # Calculate the residuals
  resid.sl = obs.in - gmsl.out.norm

  # Calculate the likelihood. The observations are not correlated. They are independent
  llik0 = sum (dnorm(resid.sl, mean=rep(0,length(resid.sl)), sd = obs.err.in, log=TRUE))

  llik = llik0 # assume residuals are independent

  llik
}
##==============================================================================
## Log-prior
log.pri = function(parameters.in, parnames.in, bound.lower.in, bound.upper.in)
{

  in.range = all(parameters.in < bound.upper.in) &
             all(parameters.in > bound.lower.in)

  if(in.range) {
    lpri=0
  } else {
    lpri = -Inf
  }

  lpri
}
##==============================================================================
## Log-posterior distribution:  posterior ~ likelihood * prior
log.post = function(parameters.in,
                    parnames.in,
                    bound.lower.in,
                    bound.upper.in,
                    obs.in, obs.err.in, obs.step.in,
                    ind.norm.in,
                    Tg.in=NULL
                    )
{
  lpri = log.pri(parameters.in , parnames.in=parnames.in,
                 bound.lower.in=bound.lower.in, bound.upper.in=bound.upper.in )
  if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
    lpost = log.lik( parameters.in=parameters.in  , parnames.in=parnames.in ,
                     obs.in=obs.in                , obs.err.in=obs.err.in   ,
                     Tg.in=Tg.in                  , ind.norm.in=ind.norm.in )
              +  lpri
  } else {
    lpost = -Inf
  }
  lpost
}
