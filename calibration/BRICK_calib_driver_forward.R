##==============================================================================
## Script to define, run and calibrate the coupled BRICK model
##
## To modify the model parameters, their prior ranges, or suggested initial values,
## edit the script 'BRICK_parameterSetup.R' before sourcing it.
##
## To modify which model components you use, edit the settings below under
## 'Which model(s) will you use?'
##
##	Questions? -- Tony Wong <twong@psu.edu
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

#setwd('~/brick/calibration')

rm(list=ls())												# Clear all previous variables

## Install packages/get libraries
library(DEoptim)

## Set the seed (for reproducibility)
#set.seed(1234)
set.seed(as.double(Sys.time())) # should yield same distributions...

## Do you want to use RCP8.5 to make projections? (l.project=TRUE)
## Or do you want to use historical data to make hindcasts? (l.project=FALSE)
## Note -- l.project=FALSE => Kriegler (2005) data, same as Urban and Keller (2010)
l.project = FALSE
begyear = 1850
endyear = 2009; #if(!l.project & endyear>2009) print('l.project and endyear not compatible')
mod.time= begyear:endyear
begyear.norm = 1961
endyear.norm = 1990
ind.norm.sl = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
tstep <- 1	# model time step (years)

## Source the models
source('../fortran/R/brickF.R')		# the full BRICK model

## Source some useful functions for manipulating data
source('../R/forcing_total.R')		# function to add up the total forcing
source('../R/compute_indices.R')	# function to determine the model and
									# data indices for comparisons
##==============================================================================

## Read the model calibration data sets
source('../calibration/DOECLIM_readData.R')		# read DOECLIM calibration data
source('../calibration/GSIC_readData.R')		# read GSIC calibration data
source('../calibration/SIMPLE_readData.R')		# GIS data, and trends in mass balance
source('../calibration/DAIS_readData.R')		# DAIS forcing data (if at all uncoupled)
source('../calibration/TE_readData.R')			# read TE data

##==============================================================================

## Estimate land water storage accounting in global mean sea level budget.
## Subtract these contributions off of the GMSL data, and add errors in quadr.
source('BRICK_estimateLandWater_IPCC.R')

##==============================================================================

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.temp,midx.ocheat,midx.gis,midx.gsic,midx.ais,midx.sl)
names(midx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"ais"   ,"sl"   )
oidx.all        = list(oidx.temp,oidx.ocheat,oidx.gis,oidx.gsic,oidx.ais,oidx.sl)
names(oidx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"ais"   ,"sl"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.temp, obs.ocheat, obs.gis, obs.gsic, obs.ais, obs.sl, obs.sl_lw)
names(obs.all) = c(    "temp"  , "ocheat"  , "gis"  , "gsic"  , "ais"  , "sl"  , "sl_lw"  )
obs.err.all        = list( obs.temp.err, obs.ocheat.err, obs.gis.err, obs.gsic.err, obs.ais.err, obs.sl.err)
names(obs.err.all) = c(    "temp"      , "ocheat"      , "gis"      , "gsic"      , "ais"	   , "sl"      )

## Set the indices for normalization that are consistent with each data set
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1850),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1992),which(mod.time==1961)) ,
		c(which(mod.time==1870),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1992),which(mod.time==1990)) )

##==============================================================================

## Get the forcing data
if(l.project) {
  forcing = read.csv( '../data/forcing_rcp85.csv', header=TRUE )
} else {
  forcing = read.csv( '../data/forcing_hindcast.csv', header=TRUE )
}

##==============================================================================
## Which model(s) will you use?
## If you want to plug in your own model, insert the relevant "luse.XXX" line
## below, as well as into the "luse.brick = ..." command.
luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # thermosteric expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = TRUE    # Antarctic ice sheet model
luse.brick = cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais)

##==============================================================================
## Define parameters and their prior ranges
## -> Note: 'parnames' is defined here, which establishes how the parameters
##    are passed around into DEoptim, MCMC, likelihood functions, and the models
source('../calibration/BRICK_parameterSetup.R')

##==============================================================================
## Use differential optimization (DEoptim) algorithm to find suitable initial
## parameters for the MCMC chains
require('DEoptim')
library(DEoptim)
source('../calibration/BRICK_DEoptim.R')
p0.deoptim=p0								# initialize optimized initial parameters
niter.deoptim=200					  		# number of iterations for DE optimization
NP.deoptim=11*length(index.model)			# population size for DEoptim (do at least 10*[N parameters])
F.deoptim=0.8								# as suggested by Storn et al (2006)
CR.deoptim=0.9								# as suggested by Storn et al (2006)
outDEoptim <- DEoptim(minimize_residuals_brick, bound.lower[index.model], bound.upper[index.model],
				DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
				parnames.in=parnames[index.model], forcing.raw=forcing       , l.project=l.project      ,
				tstep=tstep						 , slope.Ta2Tg=slope.Ta2Tg	 , intercept.Ta2Tg=intercept.Ta2Tg,
				ind.norm.data=ind.norm.data      , ind.norm.sl=ind.norm.sl   , mod.time=mod.time        ,
				oidx = oidx.all                  , midx = midx.all           , obs=obs.all              ,
				obs.err = obs.err.all            , trends.te = trends.te     , trends.ais = trends.ais	,
				luse.brick = luse.brick
				)
p0.deoptim[index.model] = outDEoptim$optim$bestmem

## Run the model and examine output at these parameter values
brick.out <- brickF(tstep=tstep,
                    mod.time=mod.time,
                    forcing.raw = forcing,
					l.project = l.project,
                    S.doeclim = p0.deoptim[match("S.doeclim",parnames)],
                    kappa.doeclim = p0.deoptim[match("kappa.doeclim",parnames)],
					alpha.doeclim = p0.deoptim[match("alpha.doeclim",parnames)],
                    T0.doeclim = p0.deoptim[match("T0.doeclim",parnames)],
                    H0.doeclim = p0.deoptim[match("H0.doeclim",parnames)],
                    beta0.gsic = p0.deoptim[match("beta0.gsic",parnames)],
					V0.gsic = p0.deoptim[match("V0.gsic",parnames)],
                    n.gsic = p0.deoptim[match("n.gsic",parnames)],
                    Gs0.gsic = p0.deoptim[match("Gs0.gsic",parnames)],
                    a.simple = p0.deoptim[match("a.simple",parnames)],
                    b.simple = p0.deoptim[match("b.simple",parnames)],
                    alpha.simple = p0.deoptim[match("alpha.simple",parnames)],
                    beta.simple = p0.deoptim[match("beta.simple",parnames)],
                    V0.simple = p0.deoptim[match("V0.simple",parnames)],
                    a.te = p0.deoptim[match("a.te",parnames)],
                    b.te = p0.deoptim[match("b.te",parnames)],
                    invtau.te = p0.deoptim[match("invtau.te",parnames)],
                    V0.te = p0.deoptim[match("V0.te",parnames)],
                    a.anto = p0.deoptim[match("anto.a",parnames)],
                    b.anto = p0.deoptim[match("anto.b",parnames)],
                    slope.Ta2Tg = slope.Ta2Tg,
                    intercept.Ta2Tg = intercept.Ta2Tg,
                    b0.dais = p0.deoptim[match("b0",parnames)],
                    slope.dais = p0.deoptim[match("slope",parnames)],
                    mu.dais = p0.deoptim[match("mu",parnames)],
                    h0.dais = p0.deoptim[match("h0",parnames)],
                    c.dais = p0.deoptim[match("c",parnames)],
                    P0.dais = p0.deoptim[match("P0",parnames)],
                    kappa.dais = p0.deoptim[match("kappa.dais",parnames)],
                    nu.dais = p0.deoptim[match("nu",parnames)],
                    f0.dais = p0.deoptim[match("f0",parnames)],
                    gamma.dais = p0.deoptim[match("gamma",parnames)],
                    alpha.dais = p0.deoptim[match("alpha.dais",parnames)]
                    )

par(mfrow=c(4,2))
	# plot 1 -- DOECLIM, temperature match
plot(obs.temp.time[oidx.temp], obs.temp[oidx.temp], pch=20, ylab='surface temperature anomaly [deg C]', xlab='year')
lines(mod.time, brick.out$temp_out, col='red', lwd=2)
	# plot 2 -- DOECLIM, ocean heat match
plot(obs.ocheat.time[oidx.ocheat], obs.ocheat[oidx.ocheat], pch=20, ylab='ocean heat uptake [10^22 J]', xlab='year')
lines(mod.time[midx.ocheat], brick.out$ocheat[midx.ocheat], col='red', lwd=2)
	# plot 3 -- GSIC match
plot(obs.gsic.time, obs.gsic, pch=20, ylab = 'sea-level equivalence [m]', xlab='year')
lines(mod.time[midx.gsic], brick.out$sl_gsic_out[midx.gsic]-brick.out$sl_gsic_out[midx.gsic[1]], col="blue", lwd=2)
	# plot 4 -- GIS match
plot(obs.gis.time[oidx.gis], obs.gis[oidx.gis], pch=20, ylab = 'GIS mass balance change (SLE [m])', xlab='year')
lines(mod.time[midx.gis], brick.out$sl_gis_out[midx.gis]-mean(brick.out$sl_gis_out[ind.norm.sl]), col="blue", lwd=2)
	# plot 5 -- TE
plot(mod.time, brick.out$sl_te_out - mean(brick.out$sl_te_out[ind.norm.sl]), type='l', col='purple', lwd=2, ylab = 'TE SLR [m]', xlab='year')
	# plot 6 -- AIS
plot(mod.time, brick.out$sl_ais_out - mean(brick.out$sl_ais_out[ind.norm.sl]), type='l', col='purple', lwd=2, ylab = 'AIS SLR [m]', xlab='year')
	# plot 7 -- SLR match
plot(obs.sl.time[oidx.sl], obs.sl[oidx.sl], pch=20, ylab = 'sea level [m]', xlab='year')
lines(mod.time, brick.out$sl_out-mean(brick.out$sl_out[ind.norm.sl]), col="purple", lwd=2)
##==============================================================================

## Fix the SIMPLE (GIS) rho.simple statistical parameter at value from the
## optimized coupled model. There are problems of convergence when it is
## free-running in the calibration, caused by the high degree of autocorrelation
## in the GIS residuals. The sigma.simple parameter also was trouble, but not
## anymore (if convergence issues arise, this might be why).

sigma.simple.fixed = NULL
rho.simple.fixed   = NULL

## Read an old rho.simple.fixed?
if(TRUE){
	rho.simple.fixed = as.numeric(read.csv('../output_calibration/rho_simple_fixed_15Feb2017.csv'))
} else {
	## If rho/sigma.simple.fixed = NULL, then will be calibrated
	resid = brick.out$sl_gis_out[midx.gis] - obs.gis[oidx.gis]
	ac = acf(resid, lag.max=5, plot=FALSE, main="")

	## If not in the declared parameter names, then fix the estimate of the
	## AR(1) process sqrt(variance) and autocorrelation
	if(is.na(match("sigma.simple",parnames))) sigma.simple.fixed = sd(resid)
	if(is.na(match("rho.simple",parnames)))   rho.simple.fixed = ac$acf[2]

	## rho.simple.fixed should be somewhere around 0.85-0.90
	## Fix it at a value you decide. Model results are insensitive to this choice,
	## provided it is realistic (between 0.85 and 1). Write to a file.

	today=Sys.Date(); today=format(today,format="%d%b%Y")
	filename=paste('../output_calibration/rho_simple_fixed_',today,'.csv', sep="")
	write.table(rho.simple.fixed, file=filename, sep=",", qmethod="double", row.names=FALSE)
}

print(paste('rho.simple.fixed=',rho.simple.fixed))

##==============================================================================
## Establish a gamma prior for drawing 1/tau.
## gamma distribution is the conjugate prior for the uncertain parameter beta
## in exponentially distributed random variable V(t):
##    dV/dt = -beta*V => V(t) = exp(-beta*t)
## For the TE model, this is not quite accurate, but similar enough (and
## previous testing with uniform priors on 1/tau confirm) that the gamma is an
## excellent approximation of the distribution of 1/tau (or tau)

invtau.hat = 1/200   # preliminary guess at what the expected value of 1/tau should be
                     # (The timescale and extent of thermal expansion of the oceans due to climate change, Marcelja, 2009)
q05 = 1/1290         # 82-1290 y are bounds given by Mengel et al (2015) for tau
q95 = 1/82           # use them as the 5-95% bounds for our 1/tau
shape.invtau = NULL	 # initialize
scale.invtau = NULL	 # initialize

## Gamma distribution has only 2 degrees of freedom. So can only fit 2 of these
## 3 requirements.

## Fit the quantiles -- it naturally arises that the mean/median of the resulting
## gamma distribution for 1/tau is around tau=200 y anyhow
rmse.quantiles <- function(parameters,q05.in,q95.in){
  shape.in=parameters[1]
  scale.in=parameters[2]
  q05.hat = qgamma(0.05, shape=shape.in, scale=scale.in, lower.tail=TRUE)
  q95.hat = qgamma(0.95, shape=shape.in, scale=scale.in, lower.tail=TRUE)
  rmse = sqrt(0.5*((q05.hat-q05.in)^2 + (q95.hat-q95.in)^2))
  return(rmse)
}

niter.deoptim=1000		   	# number of iterations for DE optimization
NP.deoptim=50	      		# population size for DEoptim (do at least 10*[N parameters])
F.deoptim=0.8						# as suggested by Storn et al (2006)
CR.deoptim=0.9					# as suggested by Storn et al (2006)
outDEoptim <- DEoptim(rmse.quantiles, c(0,0), c(100,100),
				DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
				q05.in=q05, q95.in=q95)
shape.invtau = outDEoptim$optim$bestmem[1]
scale.invtau = outDEoptim$optim$bestmem[2]
print(paste('shape.invtau=',shape.invtau,' (1.81?)'))
print(paste('scale.invtau=',scale.invtau,' (0.00275?)'))

## This should yield results somewhere in the ballpark of
##
##     shape.invtau = 1.81    and     scale.invtau = 0.00275
##
## If yours does not, check whether the gamma distribution defined by your
## shape and scale fit the invtau.hat=1/200, q05=1/1290 and q95=1/82
## requirements above. If not, try re-running, possibly with more population
## members (NP.deoptim) or more iterations (niter.deoptim).
##==============================================================================

##==============================================================================
## Now set up the coupled model calibration
##	-- will need to whip up a coupled model likelihood function file
## 	-- will calibrate with log.post(...) (in the above LL file) calling 'brick_model'
## Notes:
##	-- With DOECLIM+GSIC+TE+SIMPLE, all R models, requires ~1900s for 1e6 iterations.
##==============================================================================
## Set up and run the MCMC calibration



# 	xx   xx  xxxxx  xxxxx    xxxxx       xx    xx  xxxxxxx  xx          xx  xx
#	xx   xx  xx     xx  xx   xx          xxxx  xx  xx   xx  xx          xx  xx
#	xxxxxxx  xxxx   xxxxx    xxxx        xx xx xx  xx   xx   xx   xx   xx   xx
#	xx   xx  xx     xx  xx   xx          xx  xxxx  xx   xx    xx xxxx xx
#	xx   xx  xxxxx  xx   xx  xxxxx       xx    xx  xxxxxxx     xxx  xxx     xx



## Source the statistical models
source('../calibration/BRICK_assimLikelihood_forward.R')

## MCMC calibration
require('adaptMCMC')
library(adaptMCMC)						# use robust adaptive Metropolis
accept.mcmc = 0.234						# Optimal as # parameters->infinity
										# (Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = 1e6						# number of iterations for MCMC
gamma.mcmc = 0.5						# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)			# remove first ?? of chains for burn-in
stopadapt.mcmc = round(niter.mcmc*1.0)	# stop adapting after ?? iterations? (niter*1 => don't stop)

##==============================================================================
## Actually run the calibration
if(FALSE){
t.beg=proc.time()											# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out1 = MCMC(log.post, niter.mcmc, p0.deoptim, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
					gamma=gamma.mcmc               , list=TRUE                  , n.start=round(0.01*niter.mcmc)    ,
					parnames.in=parnames           , forcing.raw=forcing        , l.project=l.project          		,
					slope.Ta2Tg=slope.Ta2Tg		   , intercept.Ta2Tg=intercept.Ta2Tg, tstep=tstep					,
					ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm.sl    , mod.time=mod.time                 ,
					oidx = oidx.all                , midx = midx.all            , obs=obs.all                       ,
					obs.err = obs.err.all          , trends.te = trends.te      , trends.ais = trends.ais			,
					bound.lower.in=bound.lower     , bound.upper.in=bound.upper , shape.in=shape.invtau      		,
					scale.in=scale.invtau          , rho.simple.in=rho.simple.fixed, sigma.simple.in=sigma.simple.fixed,
					luse.brick=luse.brick)
t.end=proc.time()											# save timing
chain1 = amcmc.out1$samples
}

## Extend and run more MCMC samples?
if(FALSE){
t.beg=proc.time()
amcmc.extend1 = MCMC.add.samples(amcmc.out1, niter.mcmc,
					parnames.in=parnames           , forcing.raw=forcing        , l.project=l.project          		,
					slope.Ta2Tg=slope.Ta2Tg		   , intercept.Ta2Tg=intercept.Ta2Tg, tstep=tstep					,
					ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm.sl    , mod.time=mod.time                 ,
					oidx = oidx.all                , midx = midx.all            , obs=obs.all                       ,
					obs.err = obs.err.all          , trends.te = trends.te      , trends.ais = trends.ais			,
					bound.lower.in=bound.lower     , bound.upper.in=bound.upper , shape.in=shape.invtau      		,
					scale.in=scale.invtau          , rho.simple.in=rho.simple.fixed, sigma.simple.in=sigma.simple.fixed,
					luse.brick=luse.brick)
t.end=proc.time()
chain1 = amcmc.extend1$samples
}

## If you want to run 2 (or more) chains in parallel (save time, more sampling)
if(TRUE){
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.par1 = MCMC.parallel(log.post, niter.mcmc, p0.deoptim, n.chain=4, n.cpu=4,
					dyn.libs=c('../fortran/brick.so','../fortran/doeclim.so','../fortran/brick_te.so','../fortran/gsic_magicc.so','../fortran/simple.so'),
					scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
					gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
					parnames.in=parnames           , forcing.raw=forcing        , l.project=l.project          		,
					slope.Ta2Tg=slope.Ta2Tg		   , intercept.Ta2Tg=intercept.Ta2Tg, tstep=tstep					,
					ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm.sl    , mod.time=mod.time                 ,
					oidx = oidx.all                , midx = midx.all            , obs=obs.all                       ,
					obs.err = obs.err.all          , trends.te = trends.te      , trends.ais = trends.ais			,
					bound.lower.in=bound.lower     , bound.upper.in=bound.upper , shape.in=shape.invtau      		,
					scale.in=scale.invtau          , rho.simple.in=rho.simple.fixed, sigma.simple.in=sigma.simple.fixed,
					luse.brick=luse.brick)
t.end=proc.time()											# save timing
chain1=amcmc.par1[[1]]$samples
chain2=amcmc.par1[[2]]$samples
chain3=amcmc.par1[[3]]$samples
chain4=amcmc.par1[[4]]$samples
}

save.image(file = "BRICK_calib_MCMC.RData")
##==============================================================================

## Diagnostic plots
if(FALSE) {
## Check #1: History plots
par(mfrow=c(5,5))
for (pp in 1:length(parnames)) {
	plot(chain1[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="")
}
}

##==============================================================================

## Determine when (in increments of 10,000 iterations, using Gelman and Rubin
## diagnostic) the two parallel chains, which are assumed to be
## chain1=amcmc.par1[[1]]$samples and chain2=amcmc.par1[[2]]$samples
## as defined above in the MCMC.parallel(...) command.

## Initialize the testing of the Gelman and Rubin diagnostics
niter.test = seq(from=100000, to=nrow(chain1), by=50000)
gr.stat = rep(NA,length(niter.test))

## Calculate the statistic at a few spots throughout the chain. Once it is
## close to 1 (people often use GR<1.1 or 1.05), the between-chain variability
## is indistinguishable from the within-chain variability, and they are
## converged. It is only after the chains are converged that you should use the
## parameter values as posterior draws, for analysis.
for (i in 1:length(niter.test)){
  mcmc1 = as.mcmc(chain1[1:niter.test[i],])
  mcmc2 = as.mcmc(chain2[1:niter.test[i],])
	mcmc3 = as.mcmc(chain3[1:niter.test[i],])
	mcmc4 = as.mcmc(chain4[1:niter.test[i],])
  mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2, mcmc3, mcmc4))
  gr.stat[i] = gelman.diag(mcmc_chain_list)[2]
}

## Plot GR statistics as a function of iterations, decide where to cut off
## chains and use the tails of both for analysis
plot(niter.test,gr.stat)

## GR statistic gets and stays below 1.1 by iteration ...? (set this to n.burnin)
##	Wong et al 2016 -- n.min = 2e5, but discard entire first half of chain
n.burnin = 5e5
n.sample = nrow(chain1)-n.burnin
parameters1=chain1[(n.burnin+1):nrow(chain1),]
parameters2=chain2[(n.burnin+1):nrow(chain1),]
parameters3=chain3[(n.burnin+1):nrow(chain1),]
parameters4=chain4[(n.burnin+1):nrow(chain1),]
parameters.posterior = rbind(parameters1,parameters2, parameters3, parameters4)
n.parameters = ncol(parameters.posterior)

## Histograms
par(mfrow=c(5,5))
for (pp in 1:length(parnames)) {
	hist(parameters.posterior[,pp], xlab=parnames[pp], main='')
}

## Write the calibrated parameters file (netCDF version)

## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape BRICK expects, just transpose it
lmax=0
for (i in 1:length(parnames)){lmax=max(lmax,nchar(parnames[i]))}

## Name the output file
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.parameters = paste('../output_calibration/BRICK-model_calibratedParameters_control_',today,'.nc',sep="")

library(ncdf4)
dim.parameters <- ncdim_def('n.parameters', '', 1:ncol(parameters.posterior), unlim=FALSE)
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(parameters.posterior), unlim=TRUE)
parameters.var <- ncvar_def('BRICK_parameters', '', list(dim.parameters,dim.ensemble), -999)
parnames.var <- ncvar_def('parnames', '', list(dim.name,dim.parameters), prec='char')
outnc <- nc_create(filename.parameters, list(parameters.var,parnames.var))
ncvar_put(outnc, parameters.var, t(parameters.posterior))
ncvar_put(outnc, parnames.var, parnames)
nc_close(outnc)

##==============================================================================

##==============================================================================
## End
##==============================================================================
