##==============================================================================
## Script to define, run and calibrate the coupled BRICK model
##	=> this version is DOECLIM+GMSL from Rahmstorf 2007
##
## To modify the model parameters, their prior ranges, or suggested initial values,
## edit the script 'BRICK_parameterSetup_R07.R' before sourcing it.
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
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)

## Source the models
source('../fortran/R/doeclimF.R')       # the DOECLIM model
source('../R/gmsl_r07.R')								# the GMSL model

## Source some useful functions for manipulating data
source('../R/forcing_total.R')					# function to add up the total forcing
source('../R/compute_indices.R')				# function to determine the model and
																				# data indices for comparisons
##==============================================================================

## Read the model calibration data sets
source('../calibration/DOECLIM_readData.R')		# read DOECLIM calibration data
source('../calibration/GSIC_readData.R')			# read GSIC calibration data
source('../calibration/TE_readData.R')        # read TE data
source('../calibration/SIMPLE_readData.R')		# GIS data, and trends in mass balance
#source('../calibration/DAIS_readData.R')			# DAIS forcing data (if at all uncoupled)

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.temp,midx.ocheat,midx.gis,midx.gsic,midx.sl)
names(midx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )
oidx.all        = list(oidx.temp,oidx.ocheat,oidx.gis,oidx.gsic,oidx.sl)
names(oidx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.temp, obs.ocheat, obs.gis, obs.gsic, obs.sl)
names(obs.all) = c(    "temp"  , "ocheat"  , "gis"  , "gsic"  , "sl" )
obs.err.all        = list( obs.temp.err, obs.ocheat.err, obs.gis.err, obs.gsic.err, obs.sl.err)
names(obs.err.all) = c(    "temp"      , "ocheat"      , "gis"      , "gsic"      , "sl"      )

## Set the indices for normalization that are consistent with each data set
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1850),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
		c(which(mod.time==1870),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

## Set the indices of the initial condition for each sub-model
i0 = vector("list",nrow(ind.norm.data)); names(i0)=as.character(ind.norm.data[,1])

## GSIC initial conditions are actually relative to 1990 (Wigley and Raper 2005)
## Re-set these. The simulation is relative to 1990, but results and comparison
## to data is relative to 1960.
i0$gsic = which(mod.time==1990)

## GIS initial conditions are relative to 1961-1990
i0$gis = which(mod.time==1961)

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
luse.gsic     = FALSE    # glaciers and small ice caps contribution to SLR
luse.te       = FALSE    # thermosteric expansion contribution to SLR
luse.simple   = FALSE    # Greenland ice sheet model
luse.dais     = FALSE    # Antarctic ice sheet model
luse.gmsl  		= TRUE   # Example of adding your own model component - GMSL
luse.brick = cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais, luse.gmsl)

##==============================================================================
## Define parameters and their prior ranges
## -> Note: 'parnames' is defined here, which establishes how the parameters
##    are passed around into DEoptim, MCMC, likelihood functions, and the models
source('../calibration/BRICK_parameterSetup_R07.R')

##==============================================================================
## Define the coupled model
## -> need it defined before DEoptim, so we can calculate the objective function
source('../R/BRICK_coupledModel_R07.R')

##==============================================================================
## Use differential optimization (DEoptim) algorithm to find suitable initial
## parameters for the MCMC chains
require('DEoptim')
library(DEoptim)
source('../calibration/BRICK_DEoptim_R07.R')
p0.deoptim=p0													# initialize optimized initial parameters
niter.deoptim=30								  		# number of iterations for DE optimization
NP.deoptim=11*length(index.model)			# population size for DEoptim (do at least 10*[N parameters])
F.deoptim=0.8													# as suggested by Storn et al (2006)
CR.deoptim=0.9												# as suggested by Storn et al (2006)
outDEoptim <- DEoptim(minimize_residuals_brick, bound.lower[index.model], bound.upper[index.model],
				DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
				parnames.in=parnames[index.model], forcing.in=forcing        , l.project=l.project      ,
				ind.norm.data=ind.norm.data      , ind.norm.sl=ind.norm      , mod.time=mod.time        ,
				oidx = oidx.all                  , midx = midx.all           , obs=obs.all              ,
				obs.err = obs.err.all            , luse.brick = luse.brick	 , i0=i0
				)
p0.deoptim[index.model] = outDEoptim$optim$bestmem

## Run the model and examine output at these parameter values
brick.out = brick_model(parameters.in=p0.deoptim,
												parnames.in=parnames,
												forcing.in=forcing,
												l.project=l.project,
												#slope.Ta2Tg.in=slope.Ta2Tg,
												#intercept.Ta2Tg.in=intercept.Ta2Tg,
												mod.time=mod.time,
												ind.norm.data = ind.norm.data,
												ind.norm.sl = ind.norm,
												luse.brick = luse.brick,
												i0 = i0
												)

par(mfrow=c(3,2))
	# plot 1 -- DOECLIM, temperature match
plot(obs.temp.time[oidx.temp], obs.temp[oidx.temp], pch=20, ylab='surface temperature anomaly [deg C]', xlab='year')
lines(brick.out$doeclim.out$time[midx.temp], brick.out$doeclim.out$temp[midx.temp]+p0.deoptim[4], col='red', lwd=2)
	# plot 2 -- DOECLIM, ocean heat match
plot(obs.ocheat.time[oidx.ocheat], obs.ocheat[oidx.ocheat], pch=20, ylab='ocean heat uptake [10^22 J]', xlab='year')
lines(brick.out$doeclim.out$time[midx.ocheat], brick.out$doeclim.out$ocheat[midx.ocheat]+p0.deoptim[5], col='red', lwd=2)
	# plot 3 -- GSIC match
plot(obs.gsic.time, obs.gsic, pch=20, ylab = 'sea-level equivalence [m]', xlab='year')
lines(obs.gsic.time, brick.out$gsic.out[midx.gsic], col="blue", lwd=2)
	# plot 4 -- GIS match
plot(obs.gis.time[oidx.gis], obs.gis[oidx.gis], pch=20, ylab = 'GIS mass balance change (SLE [m])', xlab='year')
lines(obs.gis.time[oidx.gis], brick.out$simple.out$sle.gis[midx.gis], col="blue", lwd=2)
	# plot 5 -- SLR match
plot(obs.sl.time[oidx.sl], obs.sl[oidx.sl], pch=20, ylab = 'sea level [m]', xlab='year')
lines(brick.out$doeclim.out$time[midx.sl], brick.out$slr.out[midx.sl], col="purple", lwd=2)
##==============================================================================


##==============================================================================
## Now set up the coupled model calibration
##	-- will need to whip up a coupled model likelihood function file
## 	-- will calibrate with log.post(...) (in the above LL file) calling 'brick_model'
## Note: with DOECLIM+GMSL_R07, takes about 5 mins for 1e5 iterations, in
## parallel with 4 chains. Run on a 2012 Macbook with 2 cores x 2 threads/core.
##==============================================================================
## Set up and run the MCMC calibration

## Source the statistical models
source('../calibration/BRICK_assimLikelihood_R07.R')

## MCMC calibration
require('adaptMCMC')
library(adaptMCMC)										# use robust adaptive Metropolis
accept.mcmc = 0.234										# Optimal as # parameters->infinity
																			#	(Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = 1e5											# number of iterations for MCMC
gamma.mcmc = 0.5											# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)				# remove first ?? of chains for burn-in
stopadapt.mcmc = round(niter.mcmc*1.0)# stop adapting after ?? iterations? (niter*1 => don't stop)

##==============================================================================
## Actually run the calibration
if(FALSE){
t.beg=proc.time()											# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out1 = MCMC(log.post, niter.mcmc, p0.deoptim, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc               , list=TRUE                  , n.start=round(0.01*niter.mcmc)    ,
									parnames.in=parnames           , forcing.in=forcing         , l.project=l.project            		,
									ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time                 ,
									oidx = oidx.all                , midx = midx.all            , obs=obs.all                       ,
									obs.err = obs.err.all          , bound.lower.in=bound.lower ,	bound.upper.in=bound.upper				,
									luse.brick=luse.brick					 , i0=i0
									)
t.end=proc.time()											# save timing
chain1 = amcmc.out1$samples
}

## Extend and run more MCMC samples?
if(FALSE){
t.beg=proc.time()
amcmc.extend1 = MCMC.add.samples(amcmc.out1, niter.mcmc,
									parnames.in=parnames           , forcing.in=forcing         , l.project=l.project            		,
									ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time                 ,
									oidx = oidx.all                , midx = midx.all            , obs=obs.all                       ,
									obs.err = obs.err.all          , bound.lower.in=bound.lower ,	bound.upper.in=bound.upper				,
									luse.brick=luse.brick					 , i0=i0
									)
t.end=proc.time()
chain1 = amcmc.extend1$samples
}

## If you want to run 2 (or more) chains in parallel (save time, more sampling)
if(TRUE){
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.par1 = MCMC.parallel(log.post, niter.mcmc, p0.deoptim, n.chain=4, n.cpu=4,
									dyn.libs=c('../fortran/doeclim.so'),
									scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames           , forcing.in=forcing         , l.project=l.project            		,
									ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time                 ,
									oidx = oidx.all                , midx = midx.all            , obs=obs.all                       ,
									obs.err = obs.err.all          , bound.lower.in=bound.lower ,	bound.upper.in=bound.upper				,
									luse.brick=luse.brick					 , i0=i0
									)
t.end=proc.time()											# save timing
chain1=amcmc.par1[[1]]$samples
chain2=amcmc.par1[[2]]$samples
chain3=amcmc.par1[[3]]$samples
chain4=amcmc.par1[[4]]$samples
}

##==============================================================================

## Diagnostic plots
if(FALSE) {
## Check #1: History plots
par(mfrow=c(4,4))
for (pp in 1:length(parnames)) {
	plot(chain1[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="")
}
}

##==============================================================================

## Determine when (in increments of xx,000 iterations, using Gelman and Rubin
## diagnostic) the two parallel chains, which are assumed to be
## chain1=amcmc.par1[[1]]$samples and chain2=amcmc.par1[[2]]$samples
## as defined above in the MCMC.parallel(...) command.

## Initialize the testing of the Gelman and Rubin diagnostics
niter.test = seq(from=10000, to=nrow(chain1), by=10000)
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
##	Wong et al 2016 -- n.min = 50,000, discard entire first half of chain
n.burnin = 5e4
n.sample = nrow(chain1)-n.burnin
parameters1=chain1[(n.burnin+1):nrow(chain1),]
parameters2=chain2[(n.burnin+1):nrow(chain1),]
parameters3=chain3[(n.burnin+1):nrow(chain1),]
parameters4=chain4[(n.burnin+1):nrow(chain1),]
parameters.posterior = rbind(parameters1,parameters2, parameters3, parameters4)
n.parameters = ncol(parameters.posterior)

## Histograms
par(mfrow=c(4,4))
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
filename.parameters = paste('../output_calibration/BRICK-model_calibratedParameters_R07_',today,'.nc',sep="")

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
