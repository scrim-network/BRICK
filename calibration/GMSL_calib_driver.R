##==============================================================================
## Script to calibrate GMSL model (Rahmstorf 2007), coupled to DOECLIM,
## for comparison with BRICK
##
## Parameters to calibrate, and standard values:
## [1-9] standard DOECLIM (see Urban et al, 2010, for example, or DOECLIM
## calibration codes within BRICK)
## [10] a.gmsl = 3.4 mm/year/degC  # sensitivity of sea-level rise to temperature
## [11] T0 = 0                     # equilibrium temperature (at which there is no rise)
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

rm(list =ls()) #Clear global environment

## Setup packages and libraries
#install.packages('compiler')
#install.packages('pscl')
library(compiler)
enableJIT(3)

## Set the seed (for reproducibility)
set.seed(1234)

## Read the data forcing for hindcasts and projections. Yields:
##  obs.sl, obs.sl.time, obs.sl.err (Church and White, 2011, GMSL data)
##  obs.temp, obs.temp.time, obs.temp.err (Morice et al, 2012, HadCRUT4 temperature data)

begyear.norm = 1961
endyear.norm = 1990
begyear = 1880
endyear = 2013
mod.time=begyear:endyear
ind.norm.sl = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
source('../R/compute_indices.R')				# function to determine the model and
																				# data indices for comparisons
source('../calibration/TE_readData.R')
source('../calibration/DOECLIM_readData.R')
temp.couple = obs.temp[which(obs.temp.time==begyear):which(obs.temp.time==endyear)]
temp.couple = temp.couple - mean(temp.couple[1:20])

## Set model parameters names and prior ranges
parnames    = c('a.gmsl','Teq.gmsl')
parameters0 = c( 0.0034 , -0.3     )
bound.lower = c( 0.0    , -1.5     )
bound.upper = c( 0.0035 , 1.5      )

## Source the GMSL model
source('../R/gmsl_r07.R')

##==============================================================================

gmsl0 = gmsl_r07(
                  a=parameters0[match("a.gmsl",parnames)],
                  Teq=parameters0[match("Teq.gmsl",parnames)],
                  Tg=temp.couple
                  )
gmsl0 = gmsl0-mean(gmsl0[ind.norm.sl])

plot(obs.sl.time,obs.sl); lines(mod.time,gmsl0,type='l',col='red')

##==============================================================================





##==============================================================================
## MCMC simulation
##==============================================================================

## Source the statistical models
source('../calibration/GMSL_assimLikelihood.R')

#install.packages('adaptMCMC')
library(adaptMCMC)
accept.mcmc = 0.234										# Optimal as # parameters->infinity
																			#	(Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = 1e5											# number of iterations for MCMC
gamma.mcmc = 0.5											# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)				# remove first ?? of chains for burn-in
stopadapt.mcmc = round(niter.mcmc*1.0)# stop adapting after ?? iterations? (niter*1 => don't stop)
step.mcmc = (bound.upper-bound.lower)*.05
##==============================================================================
## Actually run the calibration

t.beg=proc.time()											# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out1 = MCMC(log.post, niter.mcmc, parameters0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames    , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  obs.in=obs.sl           , obs.err.in=obs.sl.err     , ind.norm=ind.norm.sl      ,
                  Tg.in=temp.couple       )
t.end=proc.time()											# save timing
chain1 = amcmc.out1$samples

## If you want to run 2 (or more) chains in parallel (save time, more sampling)
## Quick note about using MCMC.parallel and MCMC.add.samples:
##    You can add samples to parallel chains by simply replacing "amcmc.out1"
##    with "amcmc.par1" in the MCMC.add.samples command above BUT these R
##    functions do not have the capability to add samples in parallel. It does
##    the task, but does it in serial after the initial chains were in parallel.
if(FALSE){
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.par1 = MCMC.parallel(log.post, niter.mcmc, parameters0, n.chain=4, n.cpu=4,
                  scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames    , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  obs.in=obs.sl           , obs.err.in=obs.sl.err     , ind.norm=ind.norm.sl      ,
                  Tg.in=temp.couple       )
t.end=proc.time()											# save timing
chain1=amcmc.par1[[1]]$samples
chain2=amcmc.par1[[2]]$samples
chain3=amcmc.par1[[3]]$samples
chain4=amcmc.par1[[4]]$samples
}

##==============================================================================
## Check history plots
par(mfrow=c(4,4))
for (pp in 1:length(parnames)) {
	plot(chain1[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="")
}

## Histograms
par(mfrow=c(4,4))
for (pp in 1:length(parnames)) {
	hist(chain1[burnin:niter.mcmc,pp], xlab=parnames[pp])
}
##==============================================================================

## Determine when (in increments of xx,000 iterations, using Gelman and Rubin
## diagnostic) the two (or more) parallel chains, which are assumed to be
## chain1=amcmc.par1[[1]]$samples and chain2=amcmc.par1[[2]]$samples
## as defined above in the MCMC.parallel(...) command.

## Initialize the testing of the Gelman and Rubin diagnostics
niter.test = seq(from=50000, to=niter.mcmc, by=10000)
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

##==============================================================================

## GR statistic gets and stays below 1.1 by iteration ...? (set this to n.burnin)
n.burnin = 120000
n.sample = nrow(chain1)-n.burnin
parameters1=chain1[(n.burnin+1):nrow(chain1),]
parameters2=chain2[(n.burnin+1):nrow(chain1),]
parameters3=chain3[(n.burnin+1):nrow(chain1),]
parameters4=chain4[(n.burnin+1):nrow(chain1),]
parameters.posterior = rbind(parameters1,parameters2, parameters3, parameters4)
n.parameters = ncol(parameters.posterior)

## Fit PDFs to the parameter distributions
pdf.all=vector('list',n.parameters)
n.node=100
for (pp in 1:n.parameters){
  tmp = density(parameters.posterior[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.all[[pp]] = tmp; names(pdf.all)[pp]=parnames[pp]
}

##==============================================================================
## Plot the PDFs?
par(mfrow=c(4,4))
for (pp in 1:n.parameters){
  plot(pdf.all[[pp]]$x,pdf.all[[pp]]$y,type='l',xlab=parnames[pp],ylab='density');
}

## Plot all relative to older estimates?
dat.parameters = read.csv("../output_calibration/DAIS_calibratedParameters_29Jul2016.csv")
parameters.old = dat.parameters[1:(nrow(dat.parameters)-1), ]
n.parameters.old=ncol(parameters.old)
pdf.old=vector('list',n.parameters.old)
n.node=200
for (pp in 1:n.parameters.old){
  tmp = density(parameters.old[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.old[[pp]] = tmp; names(pdf.old)[pp]=parnames[pp]
}
par(mfrow=c(4,4))
for (pp in 1:n.parameters){
	plot(pdf.all[[pp]]$x,pdf.all[[pp]]$y,type='l',xlab=parnames[pp],col='red');
  	lines(pdf.old[[pp]]$x,pdf.old[[pp]]$y,type='l',col='black');
}
##==============================================================================

## Note on bandwidths and number of nodes for Kernel density estimate fitting:
##    With hundreds of thousands of parameters in the distributions we seek to
##    fit KDEs to, the fits and bandwidths are fairly insensitive to the number
##    of KDE nodes you plop down (n.node, above). With 300,000 = n.parameters,
##    Tony tested n.node = 20, 50, 100, 200 and 1000. The bandwidths for
##    n.node >= 50 are all the same to within 7 significant figures.

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

## Plot the distributions
par(mfrow=c(4,4))
for (p in 1:n.parameters) {
  plot(pdf.all[[p]]$x,pdf.all[[p]]$y,type='l',xlab=parnames[p],ylab='density',main="")
}

## Write the calibrated parameters file (csv version - antiquated)
#to.file = rbind(parameters.posterior,bandwidths)  # bandwidths are in the last row
#rownames(to.file)=NULL
#colnames(to.file)=parnames
#today=Sys.Date(); today=format(today,format="%d%b%Y")
#filename=paste('../output_calibration/DAIS_calibratedParameters_',today,'.csv', sep="")
#write.table(to.file, file=filename, sep=",", qmethod="double", row.names=FALSE)

# Write the calibrated parameters file (netCDF version - better)

## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape BRICK expects, just transpose it
lmax=0
for (i in 1:length(parnames)){lmax=max(lmax,nchar(parnames[i]))}

library(ncdf4)
dim.parameters <- ncdim_def('n.parameters', '', 1:ncol(parameters.posterior), unlim=FALSE)
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(parameters.posterior), unlim=TRUE)
parameters.var <- ncvar_def('DAIS_parameters', '', list(dim.parameters,dim.ensemble), -999)
parnames.var <- ncvar_def('parnames', '', list(dim.name,dim.parameters), prec='char')
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.daisparameters = paste('../output_calibration/DAIS_calibratedParameters_',today,'.nc',sep="")
outnc <- nc_create(filename.daisparameters, list(parameters.var,parnames.var))
ncvar_put(outnc, parameters.var, t(parameters.posterior))
ncvar_put(outnc, parnames.var, parnames)
nc_close(outnc)

##==============================================================================
## End
##==============================================================================
