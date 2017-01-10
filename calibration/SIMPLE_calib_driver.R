##==============================================================================
## Driver script for calibration of SIMPLE model for the Greenland Ice Sheet.
## Requires:
##    simple.R               the physical model
##    SIMPLE_readData.R      reads data to calibrate the model against
##    SIMPLE_optim.R         simplified model version to use "optim" on and find initial parameters
##    SIMPLE_likelihood.R    likelihood functions, posteriors, priors, for MCMC
##
## Parameters:
## [1] a.simple        sensitivity of equilibrium volume Veq
## [2] b.simple        equilibrium volume Veq for temperature Tg=0
## [3] alpha.simple    sensitivity of exponential decay rate
## [4] beta.simple     exponential decay rate [1/K] at Tg=0
## [5] Vgrl0           initial ice sheet volume [m SLE]
## [6] sigma.simple    AR1 error model innovation variance (whitened noise) (statistical parameter)
## [7] rho.simple      AR1 SLR residual lag-1 autocorrelation coefficient (statistical parameter)
##
## Questions? Tony Wong <twong@psu.edu>
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

## Preliminary bookkeeping
## -- Clear environment and graphical settings
rm(list =ls())
#graphics.off()

## -- Libraries
library(compiler)
enableJIT(3) # speeds up the compiling / see compiler code
library(DEoptim)

## -- Set the seed
set.seed(1234)

## -- Source the physical model
source('../fortran/R/simpleF.R')
source('../R/simple.R')
source('../R/compute_indices.R')

##==============================================================================
## Set model parameters and their prior ranges

## No rho.simple? (fixed)
#parnames   =c("a.simple" ,"b.simple" ,"alpha.simple","beta.simple","V0"      ,"sigma.simple") # parameters names
#p0         =c(-0.825     , 7.36      , 1.63e-4      , 2.85e-5     , 7.36     , 5e-4         ) # initial parameter guesses
#bound.lower=c( -4        , 5.888     , 0            , 0           , 5.888    , 0            ) # prior range lower bounds
#bound.upper=c( -1e-3     , 8.832     , 1e-3         , 1e-3        , 8.832    , 0.002        ) # prior range upper bounds
#step.mcmc  =c( 0.2       , 0.05      , 1e-5         , 1e-5        , 0.05     , 0.0001       ) # step sizes for initial MCMC
#index.model=c(1,2,3,4,5)			# which are model parameters? (index within parnames.simple)

## With rho.simple? (statistically modeled)
parnames   =c("a.simple" ,"b.simple" ,"alpha.simple","beta.simple","V0"      ,"sigma.simple","rho.simple") # parameters names
p0         =c(-0.825     , 7.36      , 1.63e-4      , 2.85e-5     , 7.36     , 5e-4         , 0.5        ) # initial parameter guesses
bound.lower=c( -4        , 5.888     , 0            , 0           , 5.888    , 0            , -.999      ) # prior range lower bounds
bound.upper=c( -1e-3     , 8.832     , 1e-3         , 1e-3        , 8.832    , 0.002        ,  .999      ) # prior range upper bounds
step.mcmc  =c( 0.2       , 0.05      , 1e-5         , 1e-5        , 0.05     , 0.0001       ,  .01       ) # step sizes for initial MCMC
index.model=c(1,2,3,4,5)			# which are model parameters? (index within parnames.simple)

##==============================================================================

##==============================================================================
## Set up model run details (length, start date, end date,...)
## -- First, specify beginning and end of model run
begyear = 1958
endyear = 2009
mod.time= begyear:endyear

## -- Next, read the data, set oidx.gis, midx.gis, indices for obs/model comparison
source('../calibration/SIMPLE_readData.R')
source('../calibration/DOECLIM_readData.R')

ibeg = which(obs.temp.time==begyear); iend = which(obs.temp.time==endyear)
if(length(ibeg)*length(iend)==0) print('ERROR - begyear/endyear outside forcing data')

forcing.time = obs.temp.time[ibeg:iend]
forcing.temp = obs.temp[ibeg:iend]

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.gis)
names(midx.all) = c(   "gis"   )
oidx.all        = list(oidx.gis)
names(oidx.all) = c(   "gis"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.gis)
names(obs.all) = c(    "gis"  )
obs.err.all        = list( obs.gis.err)
names(obs.err.all) = c(    "gis"      )

## What indices should be used to normalize in same way as data?
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1961),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
		c(which(mod.time==1990),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

##==============================================================================

##==============================================================================
## Use differential optimization to find good initial parameters for the MCMC
source('../calibration/SIMPLE_DEoptim.R')

niter.deoptim=500										# number of iterations for DE optimization
NP.deoptim=10*length(index.model)			# population size for DEoptim (do at least 10*[N parameters])
F.deoptim=0.8													# as suggested by Storn et al (2006)
CR.deoptim=0.9												# as suggested by Storn et al (2006)
outDEoptim <- DEoptim(minimize_residuals_simple, bound.lower[index.model], bound.upper[index.model],
				DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
				parnames.in=parnames,	forcing.temp.in=forcing.temp,
				oidx=oidx.all, midx=midx.all, obs=obs.all, obs.err=obs.err.all, ind.norm.data=ind.norm.data
				)
p0.deoptim = p0
p0.deoptim[index.model] = outDEoptim$optim$bestmem

## Get a preliminary model run, for better estimates of sigma.simple and rho.simple,
## the statistical parameters
simple.out = simple(a=p0.deoptim[1]   , b=p0.deoptim[2]     , alpha=p0.deoptim[3] ,
                    beta=p0.deoptim[4], V0=p0.deoptim[5], Tg=forcing.temp )
## Subtract off normalization period
itmp = ind.norm.data[match("gis",ind.norm.data[,1]),2]:ind.norm.data[match("gis",ind.norm.data[,1]),3]
simple.out$sle.gis = simple.out$sle.gis - mean(simple.out$sle.gis[itmp])

resid = simple.out$sle.gis[midx.gis] - obs.gis[oidx.gis]
rho = rep(NA,3)
ac = acf(resid, lag.max=5, plot=FALSE, main="")
rho[1]=ac$acf[1]
rho[2]=ac$acf[2]
rho[3]=ac$acf[3]
rho[4]=ac$acf[4]
rho[5]=ac$acf[5]

p0.deoptim[6:7] = c(sd(resid),rho[2])

## Plot model against obs
plot(forcing.time,simple.out$sle.gis,type='l'); points(obs.gis.time,obs.gis,pch=20);
##==============================================================================

##==============================================================================
## Set up the MCMC calibration
## -- Source the statistical models
source('../calibration/SIMPLE_assimLikelihood.R')

## -- Grab the adaptive MCMC library (use library(mcmc) if you want vanilla MCMC)
## Install relevant packages, if necessary
#install.packages('adaptMCMC')
library(adaptMCMC)

## Set up MCMC calibration settings
accept.mcmc = 0.234										# Optimal as # parameters->infinity
																			#	(Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = 1e5											# number of iterations for MCMC
gamma.mcmc = 0.5											# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)				# remove first ?? of chains for burn-in
stopadapt.mcmc = round(niter.mcmc*1.0)# stop adapting after ?? iterations? (niter*1 => don't stop)
##==============================================================================

##==============================================================================
## Actually run the MCMC calibration
t.beg=proc.time()											# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out1 = MCMC(log.post, niter.mcmc, p0.deoptim, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames, bound.lower.in=bound.lower, bound.upper.in=bound.upper,
									forcing.temp.in=forcing.temp, obs=obs.all, obs.err=obs.err.all,
									oidx=oidx.all, midx=midx.all, ind.norm.data=ind.norm.data
									)
t.end=proc.time()											# save timing
chain1 = amcmc.out1$samples

## Extend and run more MCMC samples?
if(FALSE){
t.beg=proc.time()
amcmc.extend1 = MCMC.add.samples(amcmc.out1, niter.mcmc,
								parnames.in=parnames, bound.lower.in=bound.lower, bound.upper.in=bound.upper,
								forcing.temp.in=forcing.temp, obs=obs.all, obs.err=obs.err.all,
								oidx=oidx.all, midx=midx.all, ind.norm.data=ind.norm.data
								)
t.end=proc.time()
chain1 = amcmc.extend1$samples
}

## If you want to run 2 (or more) chains in parallel (save time, more sampling)
if(FALSE){
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.par1 = MCMC.parallel(log.post, niter.mcmc, p0.deoptim, n.chain=4, n.cpu=4, dyn.libs='../fortran/simple.so',
									scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames, bound.lower.in=bound.lower, bound.upper.in=bound.upper,
									forcing.temp.in=forcing.temp, obs=obs.all, obs.err=obs.err.all,
									oidx=oidx.all, midx=midx.all, ind.norm.data=ind.norm.data
									)
t.end=proc.time()											# save timing
chain1=amcmc.par1[[1]]$samples
chain2=amcmc.par1[[2]]$samples
chain3=amcmc.par1[[3]]$samples
chain4=amcmc.par1[[4]]$samples
}

## Extend and run more MCMC samples from parallel? (note: extension is not in parallel)
if(FALSE){
t.beg=proc.time()
amcmc.extend1 = MCMC.add.samples(amcmc.par1, niter.mcmc,
								parnames.in=parnames, bound.lower.in=bound.lower, bound.upper.in=bound.upper,
								forcing.temp.in=forcing.temp, obs=obs.all, obs.err=obs.err.all,
								oidx=oidx.all, midx=midx.all, ind.norm.data=ind.norm.data
								)
t.end=proc.time()
chain1=amcmc.extend1[[1]]$samples
chain2=amcmc.extend1[[2]]$samples
}

##==============================================================================


##==============================================================================
## Calibration output diagnostics
## Check #1: History plots
par(mfrow=c(3,3))
for (pp in 1:length(parnames)) {
	plot(chain1[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="")
}

## Check #2: Heidelberger and Welch's convergence diagnostic:
heidel.diag(chain1, eps=0.1, pvalue=0.05)

## Check #3: Gelman and Rubin's convergence diagnostic:
# Converged when the potental scale reduction factor is less than 1.1
set.seed(222)
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out2 = MCMC(log.post, niter.mcmc, p0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc, gamma=0.5,
					list=TRUE, n.start=round(0.01*niter.mcmc))	# use Robust Adaptive MCMC (RAM)
chain2 = amcmc.out2$samples
t.end=proc.time()											# save timing

mcmc1 = as.mcmc(chain1)
mcmc2 = as.mcmc(chain2)

mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2))
gelman.diag(mcmc_chain_list)

# Reset seed back to original
set.seed(1234)

## Histograms
par(mfrow=c(3,3))
for (pp in 1:length(parnames)) {
	hist(chain1[burnin:niter.mcmc,pp], xlab=parnames[pp],xlim=c(bound.lower[pp],bound.upper[pp]))
}
##==============================================================================

## Determine when (in increments of ??,000 iterations, using Gelman and Rubin
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
## SIMPLE is converged within 10,000s-100,000 of iterations.
n.burnin = 100000
n.sample = nrow(chain1)-n.burnin
parameters1=chain1[(n.burnin+1):nrow(chain1),]
parameters2=chain2[(n.burnin+1):nrow(chain1),]
parameters3=chain3[(n.burnin+1):nrow(chain1),]
parameters4=chain4[(n.burnin+1):nrow(chain1),]
parameters.posterior = rbind(parameters1,parameters2, parameters3, parameters4)
n.parameters = ncol(parameters.posterior)

## Histograms
par(mfrow=c(3,3))
for (pp in 1:length(parnames)) {
	hist(parameters.posterior[,pp], xlab=parnames[pp], main='')
}

## Fit PDFs to the parameter distributions
pdf.all=vector('list',n.parameters)
n.node=200
for (pp in 1:n.parameters){
  tmp = density(parameters.posterior[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.all[[pp]] = tmp; names(pdf.all)[pp]=parnames[pp]
}

##==============================================================================
## Plot the PDFs?
par(mfrow=c(3,3))
for (pp in 1:n.parameters){
  plot(pdf.all[[pp]]$x,pdf.all[[pp]]$y,type='l',xlab=parnames[pp]);
}

##==============================================================================

## Write a CSV file with the successful parameter combinations and bandwidths
## Structure of the CSV file is as follows. 5 columns, 2*n.sample rows (2 chains
## with n.samples each).
##    First row: Parameter names.
##    Rows 2-??: The calibrated parameter values.
##    Last row:  The bandwidths. These are the standard
##               deviations of the normal distributions one should sample from.
##               The idea is that if you pick a row out of this CSV file, and
##               draw random-normally with these bandwidths (stdevs) around each
##               parameter value, you are sampling from the joint distribution.
bandwidths=rep(NA,n.parameters)
for (i in 1:n.parameters){
  bandwidths[i]=pdf.all[[i]]$bw
}

## Write calibrated parameters to a csv file
to.file = rbind(parameters.posterior,bandwidths) # bandwidths are in the last row
rownames(to.file)=NULL
colnames(to.file)=parnames
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename=paste('../output_calibration/SIMPLE_precalibrationMCMC_parameters_',today,'.csv', sep="")
write.table(to.file, file=filename, sep=",", qmethod="double", row.names=FALSE)

##==============================================================================
## End
##==============================================================================
