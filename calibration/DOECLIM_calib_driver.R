##==============================================================================
##
##	Implementation of MCMC (RAM) calibration of DOECLIM energy balance model
##	(Diffusion Ocean Energy CLIMate model, Kriegler, 2005; Tanaka et al, 2007)
##	+ robust adaptive MCMC (RAM) calibration for model parameters (Vihola, 2001)
##	+ differential evolution (DE) optimization to find suitable initial parameter guesses
##		(Storn and Price, 1997; Price et al, 2005)
##
##	Questions? -- Tony Wong <twong@psu.edu>
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

rm(list=ls())												# Clear all previous variables

## Install packages/get libraries
library(DEoptim)

## Set the seed (for reproducibility)
set.seed(1234)

## Make projections or hindcasts? When woudl you like to start/end your model?
l.project=FALSE
begyear = 1850
endyear = 2009
mod.time= begyear:endyear

## Get the forcing data
if(l.project) {
  forcing = read.csv( '../data/forcing_rcp85.csv', header=TRUE )
} else {
  forcing = read.csv( '../data/forcing_hindcast.csv', header=TRUE )
}

## Read temperature and ocean heat data for assimilation

source('../fortran/R/doeclimF.R')
source('../R/compute_indices.R')
source('../R/forcing_total.R')
source('../calibration/DOECLIM_readData.R')

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.temp,midx.ocheat)
names(midx.all) = c(   "temp"   ,"ocheat"   )
oidx.all        = list(oidx.temp,oidx.ocheat)
names(oidx.all) = c(   "temp"   ,"ocheat"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.temp, obs.ocheat)
names(obs.all) = c(    "temp"  , "ocheat"  )
obs.err.all        = list( obs.temp.err, obs.ocheat.err)
names(obs.err.all) = c(    "temp"      , "ocheat"      )

## What indices should be used to normalize in same way as data?
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1961),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
		c(which(mod.time==1990),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

## (can revise to be better later, for now use Urban and Keller, 2010, approximate values)
parnames   =c("S" ,"kappa.doeclim","alpha.doeclim","T0"  ,"H0" ,"sigma.T","sigma.H","rho.T","rho.H")	# parameters names
p0		     =c(3.1 , 3.5           , 1.1           , -0.06, -33 , 0.1     , 2       , 0.55  , 0.9   )	# initial parameter guesses
bound.lower=c(0.1 , 0.1           , 0             , -0.3 , -50 , 0.05    , 0.1     , 0     , 0     )	# prior range lower bounds
bound.upper=c(10  , 4             , 2             ,  0.3 ,   0 , 5       , 10      , 0.99  , 0.99  )	# prior range upper bounds
step.mcmc	 =c(0.16, 0.17          , 0.025         , 0.003, 0.9 , 5e-4    , 0.025   , 0.007 , 0.006 )	# step size for parameters in MCMC (proposals)
index.model=c(1,2,3,4,5)		# which are model parameters
index.stat = (length(index.model)+1):length(p0)

##==============================================================================

## Use 'DEoptim' (differential evolution optimization) to find better initial parameters
## + p0 initial parameter guesses based on Urban and Keller (2010)
## + These are okay, and work, but can improve using differential evolution optimization
##   (as long as you use a large enough vector population (at least 10*[# parameters]))

source('../calibration/DOECLIM_DEoptim.R')
niter.deoptim=200											# number of iterations for DE optimization
NP.deoptim=10*length(p0)							# population size for DEoptim (do at least 10*[N parameters])
F.deoptim=0.8													# as suggested by Storn et al (2006)
CR.deoptim=0.9												# as suggested by Storn et al (2006)
outDEoptim <- DEoptim(minimize_residuals_doeclim, bound.lower[1:5], bound.upper[1:5],
				DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
				parnames.in=parnames,	forcing.in=forcing, l.project=l.project, mod.time=mod.time,
				oidx=oidx.all, midx=midx.all, obs=obs.all, obs.err=obs.err.all, ind.norm.data=ind.norm.data
				)
p0[index.model] = outDEoptim$optim$bestmem

##==============================================================================

## MCMC calibration
## + Cite Metropolis et al (1953) and Hasting (1973) for any Metropolis-Hastings
## + Also cite Vihola (2012) if you use "adaptMCMC" library's "MCMC" with adapt=TRUE
## + log.post is in the 'DOECLIM_assimLikelihood.R' module, and defines
##   the statistical model for calibration

## Source the statistical model
source('../calibration/DOECLIM_assimLikelihood.R')

## Install (if needed) relevant packages
#install.packages('adaptMCMC')
library(adaptMCMC)										# use robust adaptive Metropolis

accept.mcmc = 0.234										# Optimal acceptance rate as # parameters->infinity
																			#	(Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = 5e4											# number of iterations for MCMC
gamma.mcmc = 0.5											# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)				# how much to remove for burn-in
stopadapt.mcmc = round(niter.mcmc*1)	# stop adapting after how long (if ever)?
t.beg=proc.time()											# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out1 = MCMC(log.post, niter.mcmc, p0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames,	forcing.in=forcing, l.project=l.project, mod.time=mod.time,
									oidx=oidx.all, midx=midx.all, obs=obs.all, obs.err=obs.err.all,
									bound.lower.in=bound.lower, bound.upper.in=bound.upper,
									ind.norm.data=ind.norm.data
									)
t.end=proc.time()											# save timing
chain1 = amcmc.out1$samples

## Extend and run more MCMC samples?
if(FALSE){
t.beg=proc.time()
amcmc.extend1 = MCMC.add.samples(amcmc.out1, niter.mcmc,
									parnames.in=parnames,	forcing.in=forcing, l.project=l.project, mod.time=mod.time,
									oidx=oidx.all, midx=midx.all, obs=obs.all, obs.err=obs.err.all,
									bound.lower.in=bound.lower, bound.upper.in=bound.upper,
									ind.norm.data=ind.norm.data
									)
t.end=proc.time()
chain1 = amcmc.extend1$samples
}

## If you want to run 2 (or more) chains in parallel (save time, more sampling)
if(FALSE){
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.par1 = MCMC.parallel(log.post, niter.mcmc, p0, n.chain=4, n.cpu=4, dyn.libs='../fortran/doeclim.so',
									scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames,	forcing.in=forcing, l.project=l.project, mod.time=mod.time,
									oidx=oidx.all, midx=midx.all, obs=obs.all, obs.err=obs.err.all,
									bound.lower.in=bound.lower, bound.upper.in=bound.upper,
									ind.norm.data=ind.norm.data
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
									parnames.in=parnames,	forcing.in=forcing, l.project=l.project, mod.time=mod.time,
									oidx=oidx.all, midx=midx.all, obs=obs.all, obs.err=obs.err.all,
									bound.lower.in=bound.lower, bound.upper.in=bound.upper,
									ind.norm.data=ind.norm.data
									)
t.end=proc.time()
chain1=amcmc.extend1[[1]]$samples
chain2=amcmc.extend1[[2]]$samples
}

##==============================================================================

## Check for convergence
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
amcmc.out2 = MCMC(log.post, niter.mcmc, p0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc, gamma=gamma.mcmc,
					list=TRUE, n.start=round(0.01*niter.mcmc), l.project=l.project)
chain2 = amcmc.out2$samples

mcmc1 = as.mcmc(chain1)
mcmc2 = as.mcmc(chain2)

mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2))
gelman.diag(mcmc_chain_list)

# Reset seed back to original
set.seed(1234)

## Histograms
par(mfrow=c(3,3))
for (pp in 1:length(parnames)) {
	hist(chain1[burnin:niter.mcmc,pp], xlab=parnames[pp])
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
## DOECLIM is converged within 10,000s of iterations.
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
filename=paste('../output_calibration/DOECLIM_precalibrationMCMC_parameters_',today,'.csv', sep="")
write.table(to.file, file=filename, sep=",", qmethod="double", row.names=FALSE)

##==============================================================================
## End
##==============================================================================
