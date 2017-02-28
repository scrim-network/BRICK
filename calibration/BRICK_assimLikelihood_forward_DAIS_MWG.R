##==============================================================================
##
##  -file = "BRICK_assimLikelihood_forward.R"
##  - Origional code written July 2014
##  - Author: Yawen Guan (yig5031@psu.edu)
##  - Edited to run SLR model by: Kelsey Ruckert (klr324@psu.edu)
##  - Edited to run DOECLIM model by: Tony Wong (twong@psu.edu)
##  - Edited for BRICK by: Tony Wong (twong@psu.edu)
##  - Edited for BRICK-forward by: Tony Wong (twong@psu.edu) (14 Feb 2017)
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
log.lik_dais = function( parameters.in,
                    parnames.in,
                    parameters.fixed,
                    parnames.fixed,
                    forcing.raw,
                    l.project=FALSE,
                    rho.simple.fixed=NULL,
                    sigma.simple.fixed=NULL,
                    slope.Ta2Tg=1,
                    intercept.Ta2Tg=0,
                    mod.time,
                    tstep,
                    ind.norm.data,
                    ind.norm.sl,
                    midx,
                    oidx,
                    obs,
                    obs.err,
                    trends.te,
                    trends.ais,
                    luse.brick
                   ){

	## Run the coupled BRICK model
    brick.out <- brickF(tstep=tstep,
                        mod.time=mod.time,
                        forcing.raw = forcing.raw,
						l.project = l.project,
                        S.doeclim = parameters.fixed[match("S.doeclim",parnames.fixed)],
                        kappa.doeclim = parameters.fixed[match("kappa.doeclim",parnames.fixed)],
						alpha.doeclim = parameters.fixed[match("alpha.doeclim",parnames.fixed)],
                        T0.doeclim = parameters.fixed[match("T0.doeclim",parnames.fixed)],
                        H0.doeclim = parameters.fixed[match("H0.doeclim",parnames.fixed)],
                        beta0.gsic = parameters.fixed[match("beta0.gsic",parnames.fixed)],
						V0.gsic = parameters.fixed[match("V0.gsic",parnames.fixed)],
                        n.gsic = parameters.fixed[match("n.gsic",parnames.fixed)],
                        Gs0.gsic = parameters.fixed[match("Gs0.gsic",parnames.fixed)],
                        a.simple = parameters.fixed[match("a.simple",parnames.fixed)],
                        b.simple = parameters.fixed[match("b.simple",parnames.fixed)],
                        alpha.simple = parameters.fixed[match("alpha.simple",parnames.fixed)],
                        beta.simple = parameters.fixed[match("beta.simple",parnames.fixed)],
                        V0.simple = parameters.fixed[match("V0.simple",parnames.fixed)],
                        a.te = parameters.fixed[match("a.te",parnames.fixed)],
                        b.te = parameters.fixed[match("b.te",parnames.fixed)],
                        invtau.te = parameters.fixed[match("invtau.te",parnames.fixed)],
                        V0.te = parameters.fixed[match("V0.te",parnames.fixed)],
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

  ## Calculate contribution from DOECLIM temperature
  llik.temp = 0
  if(!is.null(oidx$temp) & luse.brick[,"luse.doeclim"]) {

    # Get the DOECLIM statistical parameters
  	sigma.T   =parameters.fixed[match("sigma.T",parnames.fixed)]
  	rho.T     =parameters.fixed[match("rho.H"  ,parnames.fixed)]

    # Normalize temperature
	itmp <- ind.norm.data[which(ind.norm.data[,1]=='temp'),2]:ind.norm.data[which(ind.norm.data[,1]=='temp'),3]
	temperature.model <- brick.out$temp_out - mean(brick.out$temp_out[itmp])

    # Calculate the DOECLIM temperature residuals; apply AR1 error model
    resid.temp= obs$temp[oidx$temp] - temperature.model[midx$temp]
    llik.temp = logl.ar1(resid.temp, sigma.T, rho.T, obs.err$temp[oidx$temp]) # AR(1)

  }

  ## Calculate contribution from DOECLIM ocean heat
  llik.ocheat = 0
  if(!is.null(oidx$ocheat) & luse.brick[,"luse.doeclim"]) {

    # Grab the DOECLIM statistical parameters
  	sigma.H   =parameters.fixed[match("sigma.H",parnames.fixed)]
  	rho.H     =parameters.fixed[match("rho.H"  ,parnames.fixed)]

    # Normalize ocean heat uptake
    # not needed; H0 in model offsets this.

    # Calculate the DOECLIM ocean heat residuals; apply AR1 error model
    resid.ocheat= obs$ocheat[oidx$ocheat] - brick.out$ocheat[midx$ocheat]
    llik.ocheat = logl.ar1(resid.ocheat, sigma.H, rho.H, obs.err$ocheat[oidx$ocheat]) # AR(1)
  }

  ## Calculate contribution from GSIC SLR
  llik.gsic = 0
  if(!is.null(oidx$gsic) & luse.brick[,"luse.gsic"]) {

    # Grab the GSIC statistical parameters
	sigma.gsic=parameters.fixed[match("sigma.gsic",parnames.fixed)]
    rho.gsic  =parameters.fixed[match("rho.gsic"  ,parnames.fixed)]
#DEBUG
#rho.gsic <- 0.85

    # Normalize GSIC
	itmp <- ind.norm.data[which(ind.norm.data[,1]=='gsic'),2]:ind.norm.data[which(ind.norm.data[,1]=='gsic'),3]
	gsic.model <- brick.out$sl_gsic_out - mean(brick.out$sl_gsic_out[itmp])

    # Calculate the GSIC residuals; apply AR1 error model
    resid.gsic= obs$gsic[oidx$gsic] - gsic.model[midx$gsic]
    llik.gsic = logl.ar1(resid.gsic, sigma.gsic, rho.gsic, obs.err$gsic[oidx$gsic]) # AR(1)
  }

  ## Calculate contribution from thermal expansion
  llik.te = 0
  if(luse.brick[,"luse.te"]) {

    # Calculate the SLR residuals - only proceed if all TE SLR < total SLR
    # (all after the first 20 years, that is, because they are the 0 point)

#    resid.sl.te = (obs$sl[oidx$sl]-mean(obs$sl[oidx$sl[1:20]])) -
#                   (brick.out$sl_te_out[midx$sl]-mean(brick.out$sl_te_out[midx$sl[1:20]]))

#    if(all(resid.sl.te[20:length(resid.sl.te)]>0)){

      # Note 1: the trends from IPCC are in mm/year, and model output is m
      # Note 2: these calculate the least squares regression slope coefficients. It
      # is more than twice as fast to calcualte by hand like this than to use R's
      # "lm(...)" function.
      # Note 3: Need 1000*trends.mod because they're in meters, but trends.te is mm

if(FALSE){
	  trends.mod = rep(0, nrow(trends.te))
	  for (i in 1:nrow(trends.te)) {
		  x = seq(trends.te[i,6],trends.te[i,7]);                 barx = mean(x);
		  y = brick.out$sl_te_out[trends.te[i,6]:trends.te[i,7]]; bary = mean(y);
	  	  trends.mod[i] = sum( (x-rep(barx,length(x)))*(y-rep(bary,length(y))))/sum( (x-rep(barx,length(x)))^2 )
	  }
      resid.trends = 1000*trends.mod - trends.te[,1]
	  err.trends   = 0.5*(trends.te[,3]-trends.te[,2])
      llik.te = sum (dnorm(resid.trends, mean=rep(0,length(resid.trends)), sd = sqrt(err.trends^2), log=TRUE))
} else {
      #obs.err.te1 <- sqrt( 0.5*diff(trends.te[1,6:7]) ) * 0.5*(trends.te[1,3]-trends.te[1,2])/1000
      #obs.te1 <- trends.te[1,1]*0.5*diff(trends.te[1,6:7])/1000
      mod.te1 <- brick.out$sl_te_out[trends.te[1,6]+0.5*diff(trends.te[1,6:7])] - brick.out$sl_te_out[trends.te[1,6]]

      #obs.err.te2 <- sqrt( 0.5*diff(trends.te[2,6:7]) ) * 0.5*(trends.te[2,3]-trends.te[2,2])/1000
      #obs.te2 <- trends.te[2,1]*0.5*diff(trends.te[2,6:7])/1000
      mod.te2 <- brick.out$sl_te_out[trends.te[2,6]+0.5*diff(trends.te[2,6:7])] - brick.out$sl_te_out[trends.te[2,6]]

      llik.te = sum (dnorm(c(mod.te1, mod.te2), mean=c(obs.te1,obs.te2), sd=c(obs.err.te1,obs.err.te2), log=TRUE))
}

#    } else {
#      llik.te = -Inf
#    }
  }

  ## Calculate contribution from SIMPLE (Greenland Ice Sheet)
  llik.simple = 0
  if(!is.null(oidx$gis) & luse.brick[,"luse.simple"]) {

    # Get the SIMPLE statistical parameters
    sigma.simple=parameters.fixed[match("sigma.simple",parnames.fixed)]
    rho.simple  =parameters.fixed[match("rho.simple"  ,parnames.fixed)]

    # Overwrite the SIMPLE statistical parameters if values were fed into MCMC
    if(!is.null(rho.simple.fixed  )) rho.simple  =rho.simple.fixed
    if(!is.null(sigma.simple.fixed)) sigma.simple=sigma.simple.fixed

    # Normalize GIS; observations are relative to 1960-1990 average (SIMPLE_readData.R)
	itmp <- ind.norm.data[which(ind.norm.data[,1]=='gis'),2]:ind.norm.data[which(ind.norm.data[,1]=='gis'),3]
	gis.model <- brick.out$sl_gis_out - mean(brick.out$sl_gis_out[itmp])

    # Calibrate SIMPLE based on GIS data alone?
    resid.simple = obs$gis[oidx$gis] - gis.model[midx$gis] #Calculate the residuals

#    if(!all(is.finite(resid.simple))) {
#      llik.simple = -Inf
#    } else {
      llik.simple  = logl.ar1(r=resid.simple, sigma1=sigma.simple,
                              rho1=rho.simple, eps1=obs.err$gis) # AR(1) #Set up the likelihood function
#    }
  }

  ## Calculate contribution from Antarctic ice sheet
  llik.dais = 0
  if(luse.brick[,"luse.dais"]) {

    # Calculate the SLR residuals - only proceed if all AIS SLR < total SLR
    # (all after the first 20 years, that is, because they are the 0 point)

#    resid.sl.ais <- (obs$sl[oidx$sl]-mean(obs$sl[oidx$sl[1:20]])) -
#                    (brick.out$sl_ais_out[midx$sl]-mean(brick.out$sl_ais_out[midx$sl[1:20]]))

#    if(all(resid.sl.ais[20:length(resid.sl.ais)]>0)){

        var.dais <- parameters.in[match("var.dais",parnames.in)]

        # First part is from Shepherd et al 2012 instrumental point (from
        # Ruckert et al 2017, or Wong et al 2017)
        itmp            <- ind.norm.data[which(ind.norm.data[,1]=='ais'),2]:ind.norm.data[which(ind.norm.data[,1]=='ais'),3]
        ais.model       <- brick.out$sl_ais_out - mean(brick.out$sl_ais_out[itmp])
        ais.instr.resid <- obs$ais[oidx$ais] - ais.model[midx$ais]

#        llik.instr <- sum(dnorm(ais.instr.resid, mean=0, sd=sqrt(var.dais + obs.err$ais[oidx$ais]^2), log=TRUE))
#        llik.instr <- sum(dnorm(ais.instr.resid, mean=0, sd=obs.err$ais[oidx$ais], log=TRUE))
        llik.instr <- sum(dnorm(ais.instr.resid, mean=0, sd=sqrt(var.dais), log=TRUE))

if(FALSE){
        # Second part is from IPCC trends
        trends.mod = rep(0, nrow(trends.ais))
        for (i in 1:nrow(trends.ais)) {
            x = seq(trends.ais[i,6],trends.ais[i,7]);                  barx = mean(x);
            y = brick.out$sl_ais_out[trends.ais[i,6]:trends.ais[i,7]]; bary = mean(y);
            trends.mod[i] = sum( (x-rep(barx,length(x)))*(y-rep(bary,length(y))))/sum( (x-rep(barx,length(x)))^2 )
        }
        resid.trends <- 1000*trends.mod - trends.ais[,1]
        err.trends   <- 0.5*(trends.ais[,3]-trends.ais[,2])
        llik.trends  <- sum(dnorm(resid.trends, mean=rep(0,length(resid.trends)), sd=sqrt(err.trends^2), log=TRUE))
} else {
      #obs.err.ais1 <- sqrt( 0.5*diff(trends.ais[1,6:7]) ) * 0.5*(trends.ais[1,3]-trends.ais[1,2])/1000
      #obs.ais1 <- trends.ais[1,1]*0.5*diff(trends.ais[1,6:7])/1000
      mod.ais1 <- brick.out$sl_te_out[trends.ais[1,6]+0.5*diff(trends.ais[1,6:7])] - brick.out$sl_te_out[trends.ais[1,6]]

      #obs.err.ais2 <- sqrt( 0.5*diff(trends.ais[2,6:7]) ) * 0.5*(trends.ais[2,3]-trends.ais[2,2])/1000
      #obs.ais2 <- trends.ais[2,1]*0.5*diff(trends.ais[2,6:7])/1000
      mod.ais2 <- brick.out$sl_te_out[trends.ais[2,6]+0.5*diff(trends.ais[2,6:7])] - brick.out$sl_te_out[trends.ais[2,6]]

      llik.trends = sum (dnorm(c(mod.ais1, mod.ais2), mean=c(obs.ais1,obs.ais2), sd=c(obs.err.ais1,obs.err.ais2), log=TRUE))
}
        llik.dais    <- llik.instr + llik.trends
#    } else {
#        llik.dais    <- -Inf
#    }
  }

  # Calculate contribution from total sea level rise
  # (subtracting out land water storage trends and adding errors in quadrature
  # to errors from GMSL; Church et al 2013 (IPCC AR5), Church and White 2011)
  var.gmsl        <- parameters.fixed[match("var.gmsl",parnames.fixed)]
  rho.gmsl.1900   <- parameters.fixed[match("rho.gmsl.1900",parnames.fixed)]
  sigma.gmsl.1900 <- parameters.fixed[match("sigma.gmsl.1900",parnames.fixed)]

  sl.model <- brick.out$sl_out - mean(brick.out$sl_out[ind.norm.sl])
#  resid.sl_lw.1900 <- obs$sl_lw$r1900$sl - sl.model[obs$sl_lw$r1900$midx]
#  resid.sl_lw.1970 <- obs$sl_lw$r1970$sl - sl.model[obs$sl_lw$r1970$midx]
#  resid.sl_lw.1992 <- obs$sl_lw$r1992$sl - sl.model[obs$sl_lw$r1992$midx]
  resid.sl_lw.1900 <- obs$sl_lw$r1900$sl - (sl.model[obs$sl_lw$r1900$midx]-sl.model[max(obs$sl_lw$r1900$midx)])
  resid.sl_lw.1970 <- obs$sl_lw$r1970$sl - (sl.model[obs$sl_lw$r1970$midx]-sl.model[max(obs$sl_lw$r1970$midx)])
  resid.sl_lw.1992 <- obs$sl_lw$r1992$sl - (sl.model[obs$sl_lw$r1992$midx]-sl.model[max(obs$sl_lw$r1992$midx)])

  # AR1 model for llik.sl.1900
  llik.sl.1900  = logl.ar1(r=resid.sl_lw.1900, sigma1=sigma.gmsl.1900,
                           rho1=rho.gmsl.1900, eps1=obs$sl_lw$r1900$err) # AR(1) #Set up the likelihood function
  # Others are normally distributed (preliminary check of ACF is okay)
  llik.sl.1970 <- sum(dnorm(resid.sl_lw.1970, mean=rep(0,length(resid.sl_lw.1970)), sd=sqrt(var.gmsl+obs$sl_lw$r1970$err^2), log=TRUE))
  llik.sl.1992 <- sum(dnorm(resid.sl_lw.1992, mean=rep(0,length(resid.sl_lw.1992)), sd=sqrt(var.gmsl+obs$sl_lw$r1992$err^2), log=TRUE))
#  llik.sl      <- llik.sl.1900 + llik.sl.1970 + llik.sl.1992
  llik.sl      <- llik.sl.1900 + llik.sl.1992   # cut out 1970-2009 because it will be correalted with 1900-1990 and 1992-2009
#  llik.sl      <- llik.sl.1900

  # Assume residual time series are independent
  llik = llik.temp + llik.ocheat + llik.gsic + llik.te + llik.simple + llik.dais + llik.sl
#  llik = llik.dais

  return(llik)
}
##==============================================================================
## (log of the) prior probability
log.pri_dais = function(parameters.in,
                   parnames.in,
                   parameters.fixed,
                   parnames.fixed,
                   bound.lower.fixed,
                   bound.upper.fixed,
                   bound.lower,
                   bound.upper,
                   shape.invtau,
                   scale.invtau,
                   dais.prior.fit,
                   gsic.prior.fit,
                   te.prior.fit,
                   simple.prior.fit,
                   anto.prior.fit,
                   doeclim.prior.fit)
{

    # Get the model and statistical parameters (only ones without uniform prior)
    ind.S          <- match("S.doeclim",parnames.fixed)
    ind.invtau     <- match("invtau.te",parnames.fixed)
    ind.Gs0        <- match("Gs0.gsic",parnames.fixed)
    ind.V0.gsic    <- match("V0.gsic",parnames.fixed)
    ind.V0.te      <- match("V0.te",parnames.fixed)
    ind.V0.simple  <- match("V0.simple",parnames.fixed)
    ind.vdais      <- match("var.dais",parnames.in)
    ind.vgmsl      <- match("var.gmsl",parnames.fixed)
    ind.rgsic      <- match("rho.gsic",parnames.fixed)
    ind.gsic       <- c(match("beta0.gsic",parnames.fixed), match("n.gsic",parnames.fixed), match("rho.gsic",parnames.fixed))
    ind.te         <- c(match('a.te',parnames.fixed), match('b.te',parnames.fixed), match('invtau.te',parnames.fixed))
    ind.doeclim    <- c(match('S.doeclim',parnames.fixed)    , match('kappa.doeclim',parnames.fixed),
                        match('alpha.doeclim',parnames.fixed), match('T0.doeclim',parnames.fixed)   ,
                        match('H0.doeclim',parnames.fixed))
    ind.simple     <- c(match('a.simple',parnames.fixed)    , match('b.simple',parnames.fixed),
                        match('alpha.simple',parnames.fixed), match('beta.simple',parnames.fixed))
    ind.anto       <- c(match('anto.a',parnames.in), match('anto.b',parnames.in))
    ind.dais       <- match('gamma',parnames.in):match('slope',parnames.in)
    lpri.S         <- 0
    lpri.invtau    <- 0
    lpri.dais      <- 0
    lpri.gsic      <- 0
    lpri.doeclim   <- 0
    lpri.te        <- 0
    lpri.simple    <- 0
    lpri.anto      <- 0
    lpri.Gs0       <- 0
    lpri.V0.gsic   <- 0
    lpri.V0.te     <- 0
    lpri.V0.simple <- 0
    lpri.vdais     <- 0
    lpri.vgmsl     <- 0
    lpri.rgsic     <- 0

    # some of the parameters have infinite support in one or both directions
    ind.inf.up <- c(ind.vdais)
    ind.inf.dn <- c()
    #in.range <- all(parameters.in > bound.lower) &
    #            all(parameters.in < bound.upper)
	in.range <- all(parameters.in > bound.lower) &
                all(parameters.in[-ind.inf.up] < bound.upper[-ind.inf.up])

	if(in.range){
        lpri.uni = 0									# Sum of all uniform priors (log(1)=0)
#        lpri.S         <- log(dcauchy(parameters.in[ind.S],location=3,scale=2) / 	# S has truncated Cauchy(3,2) prior
#					         (pcauchy(bound.upper[ind.S],location=3,scale=2)-pcauchy(bound.lower[ind.S],location=3,scale=2)))
#        lpri.invtau    <- dgamma( parameters.in[ind.invtau], shape=shape.invtau, scale=scale.invtau, log=TRUE )
        lpri.Gs0       <- dbeta(range.to.beta(parameters.fixed[ind.Gs0], bound.lower.fixed[ind.Gs0], bound.upper.fixed[ind.Gs0]),
                          shape1=Gs0.gsic.prior[1], shape2=Gs0.gsic.prior[2], log=TRUE)
        lpri.V0.gsic   <- dbeta(range.to.beta(parameters.fixed[ind.V0.gsic], bound.lower.fixed[ind.V0.gsic], bound.upper.fixed[ind.V0.gsic]),
                          shape1=V0.gsic.prior[1], shape2=V0.gsic.prior[2], log=TRUE)
        lpri.V0.te     <- dbeta(range.to.beta(parameters.fixed[ind.V0.te], bound.lower.fixed[ind.V0.te], bound.upper.fixed[ind.V0.te]),
                          shape1=V0.te.prior[1], shape2=V0.te.prior[2], log=TRUE)
        lpri.V0.simple <- dbeta(range.to.beta(parameters.fixed[ind.V0.simple], bound.lower.fixed[ind.V0.simple], bound.upper.fixed[ind.V0.simple]),
                          shape1=V0.simple.prior[1], shape2=V0.simple.prior[2], log=TRUE)
#        lpri.rgsic     <- dbeta(range.to.beta(parameters.fixed[ind.rgsic], bound.lower.fixed[ind.rgsic], bound.upper.fixed[ind.rgsic]),
#                          shape1=rho.gsic.prior[1], shape2=rho.gsic.prior[2], log=TRUE)
        lpri.vdais     <- log(densigamma(x=parameters.in[ind.vdais], alpha=2.2, beta=3e-5) )
        lpri.vgmsl     <- log(densigamma(x=parameters.fixed[ind.vgmsl], alpha=6.26, beta=3.27e-4) )
        lpri.gsic      <- log(dmsn( parameters.fixed[ind.gsic], gsic.prior.fit$beta, gsic.prior.fit$Omega, gsic.prior.fit$alpha)/gsic.prior.fit$cnorm)
        lpri.te        <- log(dmsn( parameters.fixed[ind.te], te.prior.fit$beta, te.prior.fit$Omega, te.prior.fit$alpha)/te.prior.fit$cnorm)
        lpri.doeclim   <- log(dmsn( parameters.fixed[ind.doeclim], doeclim.prior.fit$beta, doeclim.prior.fit$Omega, doeclim.prior.fit$alpha)/doeclim.prior.fit$cnorm)
        lpri.simple    <- log(dmsn( parameters.fixed[ind.simple], simple.prior.fit$beta, simple.prior.fit$Omega, simple.prior.fit$alpha)/simple.prior.fit$cnorm)
        lpri.anto      <- log(dmsn( parameters.in[ind.anto], anto.prior.fit$beta, anto.prior.fit$Omega, anto.prior.fit$alpha)/anto.prior.fit$cnorm)
        lpri.dais      <- log(dmsn( parameters.in[ind.dais], dais.prior.fit$beta, dais.prior.fit$Omega, dais.prior.fit$alpha)/dais.prior.fit$cnorm)
        lpri = lpri.uni +
               lpri.doeclim + lpri.gsic + lpri.dais + lpri.te + lpri.simple + lpri.anto +
               lpri.vdais + lpri.vgmsl + lpri.Gs0 + lpri.V0.gsic + lpri.V0.te + lpri.V0.simple
	} else {
        lpri = -Inf
	}

	return(lpri)
}
##==============================================================================
## (log of the) posterior distribution:  posterior ~ likelihood * prior
log.post_dais = function(  parameters.in,
                      parnames.in,
                      parameters.fixed,
                      parnames.fixed,
                      bound.lower.fixed,
                      bound.upper.fixed,
                      forcing.raw,
                      bound.lower,
                      bound.upper,
                      l.project=FALSE,
                      rho.simple.fixed=NULL,
                      sigma.simple.fixed=NULL,
                      shape.invtau,
                      scale.invtau,
                      slope.Ta2Tg=1,
                      intercept.Ta2Tg=0,
                      mod.time,
                      tstep=1,
                      ind.norm.data,
                      ind.norm.sl,
                      midx,
                      oidx,
                      obs,
                      obs.err,
                      trends.te,
                      trends.ais,
                      dais.prior.fit=NULL,
                      gsic.prior.fit=NULL,
                      te.prior.fit=NULL,
                      simple.prior.fit=NULL,
                      anto.prior.fit=NULL,
                      doeclim.prior.fit=NULL,
                      luse.brick
                      ){

  llik = 0
  lpri = log.pri_dais( parameters.in      = parameters.in,
                  parnames.in        = parnames.in,
                  parameters.fixed   = parameters.fixed,
                  parnames.fixed     = parnames.fixed,
                  bound.lower.fixed  = bound.lower.fixed,
                  bound.upper.fixed  = bound.upper.fixed,
                  bound.lower        = bound.lower,
                  bound.upper        = bound.upper,
                  shape.invtau       = shape.invtau,
                  scale.invtau       = scale.invtau,
                  dais.prior.fit     = dais.prior.fit,
                  gsic.prior.fit     = gsic.prior.fit,
                  te.prior.fit       = te.prior.fit,
                  simple.prior.fit   = simple.prior.fit,
                  anto.prior.fit     = anto.prior.fit,
                  doeclim.prior.fit  = doeclim.prior.fit)
  if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
    llik = log.lik_dais( parameters.in=parameters.in,
                    parnames.in=parnames.in,
                    parameters.fixed=parameters.fixed,
                    parnames.fixed=parnames.fixed,
                    forcing.raw=forcing.raw,
                    l.project=l.project,
                    rho.simple.fixed=rho.simple.fixed,
                    sigma.simple.fixed=sigma.simple.fixed,
                    slope.Ta2Tg=slope.Ta2Tg,
                    intercept.Ta2Tg=intercept.Ta2Tg,
                    mod.time,
                    tstep=tstep,
                    ind.norm.data=ind.norm.data,
                    ind.norm.sl=ind.norm.sl,
                    midx=midx,
                    oidx=oidx,
                    obs=obs,
                    obs.err=obs.err,
                    trends.te=trends.te,
                    trends.ais=trends.ais,
                    luse.brick=luse.brick)
    lpost = llik + lpri
  } else {
  	lpost = -Inf
  }
  return(lpost)
}
##==============================================================================
## End
##==============================================================================
