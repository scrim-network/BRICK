##==============================================================================
#
#  -original file = "DOECLIM_assimLikelihood_AR1hetero.R"   Origional code written July 2014
#  - Author: Yawen Guan (yig5031@psu.edu)
#  - Edited to run SLR model by: Kelsey Ruckert (klr324@psu.edu)
#  - Edited to run DOECLIM model by: Tony Wong (twong@psu.edu)
#  - Edited to run in BRICK model framework July 2016 by: Tony Wong (twong@psu.edu)
#
#  -This function computes the log likelihood for a zero-mean AR1 process from
#       observations as described in  Ruckert et al. (2016).
#       For further description and references, please read the paper
#       and the appendix.
#
#   -NOTE: Descriptions of how to use this for other observation and models
#       can be found in the R package in review "VAR1"
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
## AR1 model (centered at 0 -- X(t)=rho*X(t-1)+eps(t), eps(t)~N(0,sigma1)
## For clarity -- sigma.proc is for AR(1) process;
##             -- sigma1 (sampled) is the whitened innovation sigma
##
## Estimate the log likelihood of the AR1 process
logl.ar1 = function(r,sigma1,rho1,eps1=0) # default obs error is 0
{
  n = length(r) # r is the residuals

	logl=0
	if(n>1) {
    w = r[2:n] - rho1*r[1:(n-1)] # this process whitens the residuals
    logl = logl + sum(dnorm(w,sd=sqrt((sigma1)^2+(eps1[c(-1)])^2),log=TRUE)) # add in the sum of
			# density of the whitened residuals with a standard deviation of the
			# variance and the obs. errors
      # eps1 using errors estimates from temp and ocean heat which change with year -> heteroscedastic
  }
  return(logl)
}
##==============================================================================
## rest of the statistical model
##==============================================================================
log.lik = function(	parameters.in,
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

	# Pluck off the model and statistical parameters
  S            =parameters.in[match("S"            ,parnames.in)]
	kappa.doeclim=parameters.in[match("kappa.doeclim",parnames.in)]
	alpha.doeclim=parameters.in[match("alpha.doeclim",parnames.in)]
	T0           =parameters.in[match("T0"           ,parnames.in)]
	H0           =parameters.in[match("H0"           ,parnames.in)]
	sigma.T      =parameters.in[match("sigma.T"      ,parnames.in)]
	rho.T        =parameters.in[match("rho.H"        ,parnames.in)]
	sigma.H      =parameters.in[match("sigma.H"      ,parnames.in)]
	rho.H        =parameters.in[match("rho.H"        ,parnames.in)]

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

	llik.temp = 0

  if(!is.null(oidx$temp)) {
      resid.temp = obs$temp[oidx$temp] - (model.out$temp[midx$temp]+T0)
      llik.temp = logl.ar1(resid.temp, sigma.T, rho.T, obs.err$temp[oidx$temp]) # AR(1)
  }

  llik.ocheat = 0
  if(!is.null(oidx$ocheat)) {
      resid.ocheat = obs$ocheat[oidx$ocheat] - (model.out$ocheat[midx$ocheat]+H0)
      llik.ocheat = logl.ar1(resid.ocheat, sigma.H, rho.H, obs.err$ocheat[oidx$ocheat]) # AR(1)
  }

	# Assume residual time series are independent
    llik = llik.temp + llik.ocheat

	return(llik)
}
##==============================================================================
## (log of the) prior probability
log.pri = function( parameters.in,
                    parnames.in,
                    bound.lower.in,
                    bound.upper.in,
                    l.informed.prior.S
                    ){

	# Pluck off the model and statistical parameters
  ind.S = match("S"            ,parnames.in)
  S            =parameters.in[ind.S]
	kappa.doeclim=parameters.in[match("kappa.doeclim",parnames.in)]
	alpha.doeclim=parameters.in[match("alpha.doeclim",parnames.in)]
	T0           =parameters.in[match("T0"           ,parnames.in)]
	H0           =parameters.in[match("H0"           ,parnames.in)]
	sigma.T      =parameters.in[match("sigma.T"      ,parnames.in)]
	rho.T        =parameters.in[match("rho.H"        ,parnames.in)]
	sigma.H      =parameters.in[match("sigma.H"      ,parnames.in)]
	rho.H        =parameters.in[match("rho.H"        ,parnames.in)]

	
	in.range.vec = ((parameters.in >= bound.lower.in) & (parameters.in <= bound.upper.in))
	if(l.informed.prior.S) {
	  in.range = all(in.range.vec[-ind.S])
	  lpri.S = dlnorm(S, meanlog=1.10704, sdlog=0.264, log=TRUE)
	} else {
	  in.range = all(in.range.vec)
	  lpri.S = 0
	}

	if(in.range){
		lpri.uni = 0									# Sum of all uniform priors (log(1)=0)
		#lpri.S = log(dcauchy(S,location=3,scale=2) / 	# S has truncated Cauchy(3,2) prior
		#			(pcauchy(bound.upper.in[1],location=3,scale=2)-pcauchy(bound.lower.in[1],location=3,scale=2)))
		lpri = lpri.uni + lpri.S
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
                      l.informed.prior.S=FALSE,
                      mod.time,
                      midx,
                      oidx,
                      obs,
                      obs.err,
                      ind.norm.data
                      ){

	lpri = log.pri( parameters.in=parameters.in,
                  parnames.in=parnames.in,
                  bound.lower.in=bound.lower.in,
                  bound.upper.in=bound.upper.in,
	                l.informed.prior.S=l.informed.prior.S
                  )
  	if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
    	lpost = log.lik(  parameters.in=parameters.in,
                        parnames.in=parnames.in,
                        forcing.in=forcing.in,
                        l.project=l.project,
                        mod.time=mod.time,
                        midx=midx,
                        oidx=oidx,
                        obs=obs,
                        obs.err=obs.err,
                        ind.norm.data
                        ) + lpri
  	} else {
    	lpost = -Inf
  	}
  	return(lpost)
}
##==============================================================================
## End
##==============================================================================
