##==============================================================================
##	Implementation of Markov chain Monte Carlo calibration of
##	GSIC (Glaciers and Small Ice Caps model, Wigley and Raper, 2005)
##	+ robust adaptive MCMC (RAM) calibration for model parameters (Vihola, 2001)
##	+ differential evolution (DE) optimization to find suitable initial parameter guesses
##		(Storn and Price, 1997; Price et al, 2005)
##
##	Original code by Kelsey Ruckert <klr324@psu.edu>
##	Modifications for BRICK framework by Tony Wong <twong@psu.edu> and
##										 Alexander Bakker <bakker@psu.edu>
##
##	Questions? -- Tony Wong <twong@psu.edu>
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
## Preliminary bookkeeping
## -- Clear environment and graphical settings
rm(list =ls())

## -- Libraries
library(compiler)
library(DEoptim)
enableJIT(3) # speeds up the compiling / see compiler code

## -- Set the seed
set.seed(1234)

## -- Source the physical model
source('../fortran/R/GSIC_magiccF.R')
##==============================================================================

##==============================================================================
## Set up model run details
## -- First, specify the beginning and end of the simulation, and timestep (years)
begyear = 1890
endyear = 2003
mod.time= begyear:endyear

## -- Next, read the data and set indices for model/data comparison.
##		use the HadCRUT4 temperature data (which is used to calibrate DOECLIM)
source('../R/compute_indices.R')
source("../calibration/DOECLIM_readData.R")
source("../calibration/GSIC_readData.R")

## -- Trim/extend forcing data, surface temperature anomaly, to match the desired
##		begyear and endyear.

ibeg = which(obs.temp.time==begyear); iend = which(obs.temp.time==endyear)
if(length(ibeg)*length(iend)==0) print('ERROR - begyear/endyear outside forcing data')

forcing.time = obs.temp.time[ibeg:iend]
forcing.temp = obs.temp[ibeg:iend]

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.gsic)
names(midx.all) = c(   "gsic"   )
oidx.all        = list(oidx.gsic)
names(oidx.all) = c(   "gsic"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.gsic)
names(obs.all) = c(    "gsic"  )
obs.err.all        = list( obs.gsic.err)
names(obs.err.all) = c(    "gsic"      )

## What indices should be used to normalize in same way as data?
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1961),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
		c(which(mod.time==1990),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

##==============================================================================

##==============================================================================
## Define parameters and their prior ranges and step sizes for MCMC

parnames   =c("beta0","V0.gsic","n"  ,"Gs0"    , "sigma.gsic", "rho.gsic")	# parameters names
p0		     =c(0.00058, 0.41    , 0.82, 0.0     , 0.00045     , 0.5       )	# initial parameter guesses
bound.lower=c(0      , 0.3     , 0.55, -0.0041 , 0           , -0.999    )	# prior range lower bounds
bound.upper=c(0.041  , 0.5     , 1.0 ,  0.0041 , 0.00150     ,  0.999    )	# prior range upper bounds
step.mcmc	 =c(0.01   , 0.01    , 0.1 , 0.01    , 0.0001      , 0.01      )	# step size for parameters in MCMC (proposals)
index.model=c(1,2,3,4)			# which are model parameters? (index within parnames.gsic)

##==============================================================================

##==============================================================================
## Find Initial Parameter & Initial Hindcast
#Initial parameters:
# [1] beta0 = 0.000577 [m/yr/C]   # initial mass balance sensitivity (how long it takes GSIC to respond to increasing temps) [m/yr/C]
# [2] V0.gsic=0.41               	# (meter SLE) total value the GSIC can contribute to sea level
# [3] n = 0.82 [dimensonless]    	# (exponent) relates to the speed of melting (larger the n the less rapid the melt) [dimensionless]
# [4] Gs0 or Gs[1] = 0.01945 [m]	# Intial value of GSIC contribution to sea-level

## Run DEoptim in R to find good initial parameters
source("../calibration/GSIC_DEoptim.R")					# source the DE optimization model, includes residual calculation
niter.deoptim=1000								# number of iterations for DE optimization
NP.deoptim=10*length(index.model)	# population size for DEoptim
																	#   (do at least 10*[N parameters])
F.deoptim=0.8											# as suggested by Storn et al (2006)
CR.deoptim=0.9										# as suggested by Storn et al (2006)
outDEoptim <- DEoptim(minimize_residuals_gsic, bound.lower[index.model], bound.upper[index.model],
			DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
			parnames.in=parnames,	forcing.temp.in=forcing.temp,
			oidx=oidx.all, midx=midx.all, obs=obs.all, obs.err=obs.err.all, ind.norm.data=ind.norm.data
			)
p0.deoptim = p0
p0.deoptim[index.model] = outDEoptim$optim$bestmem

#Run the model with the initial parameters to create a simulation of the observations
gsic.out = gsic_magiccF(beta0=p0.deoptim[1], V0=p0.deoptim[2], n=p0.deoptim[3],
						Gs0=p0.deoptim[4], Tg=forcing.temp)
# Must normalize by the melt at year before the first time of observations, if
# begyear != 1960, since this is commensurate with what the data are
# (../GSIC-MAGICC/GSICobservations_UPDATED.csv)
gsic.out.norm = gsic.out - gsic.out[which(forcing.time==obs.gsic.time[oidx.gsic[1]])-1]

# Plot the fit?
plot(obs.gsic.time, obs.gsic, pch=20, ylab = "Sea-level equivalence [m]", xlab="Year")
lines(forcing.time, gsic.out.norm, col="blue", lwd=2)

##==============================================================================
## Calculate the Residuals  & AR(1) Coefficient
#Calculate Residuals
gsic.residuals = obs.gsic[oidx.gsic] - gsic.out.norm[midx.gsic]

## Estimate and save the lag-1 autocorrelation coefficient (rho[2])
rho=rep(NA,3)
ac=acf(gsic.residuals, lag.max=5, plot=TRUE, main="")	# apply auto-correlation to determine
														#	correlation coefficients
rho[1]=ac$acf[1]
rho[2]=ac$acf[2]
rho[3]=ac$acf[3]
rho[4]=ac$acf[4]
rho[5]=ac$acf[5]

##==============================================================================
## Run MCMC Calibration

# Step 3: source the statistical model
source("../calibration/GSIC_assimLikelihood.R")

# Step 4: Set up the initial parameters from the DEoptim best guess parameters
p0.deoptim[5:6] = c(sd(gsic.residuals), rho[2])


# Install libraries, if needed
#install.packages('adaptMCMC')
library(adaptMCMC)										# use robust adaptive Metropolis-Hastings

# Step 5: Set up the step size, burnin, and number of iterations to run
niter.mcmc = 3e5										# number of iterations for MCMC
gamma.mcmc = 0.5											# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = seq(1,0.01*niter.mcmc,1)						# 1% burnin
accept.mcmc = 0.234										# Optimal as # parameters->infinity
														#	(Gelman et al, 1996; Roberts et al, 1997)

# Step 6: Run the MCMC chain
t.beg=proc.time()
amcmc.out1 = MCMC(log.post, niter.mcmc, p0.deoptim, scale=step.mcmc, adapt=TRUE,
					acc.rate=accept.mcmc, gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
					parnames.in=parnames, bound.lower.in=bound.lower, bound.upper.in=bound.upper,
					forcing.temp.in=forcing.temp, obs=obs.all, obs.err=obs.err.all,
					oidx=oidx.all, midx=midx.all, ind.norm.data=ind.norm.data
					)
t.end=proc.time()
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
amcmc.par1 = MCMC.parallel(log.post, niter.mcmc, p0, n.chain=4, n.cpu=4, dyn.libs='../fortran/gsic_magicc.so',
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
## Test parameter convergence
library(coda)

## Check #1: History plots
par(mfrow=c(3,2))
for (pp in 1:length(parnames)) {
	plot(chain1[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="",
					  ylim=c(bound.lower[pp],bound.upper[pp]))
}

## Check #2: Heidelberger and Welch's convergence diagnostic:
heidel.diag(chain1, eps=0.1, pvalue=0.05)

## Check #3: Gelman and Rubin's convergence diagnostic:
# Converged when the potental scale reduction factor is less than 1.1
set.seed(222)
amcmc.out2 = MCMC(log.post, niter.mcmc, p0, scale=step, adapt=TRUE, acc.rate=accept.mcmc, gamma=0.5,
					list=TRUE, n.start=round(0.01*niter.mcmc))
chain2 = amcmc.out2$samples

mcmc1 = as.mcmc(chain1)
mcmc2 = as.mcmc(chain2)

mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2))
gelman.diag(mcmc_chain_list)

# Reset seed back to original
set.seed(1234)

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
## GSIC is converged within 10,000s of iterations.
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
filename=paste('../output_calibration/GSIC_precalibrationMCMC_parameters_',today,'.csv', sep="")
write.table(to.file, file=filename, sep=",", qmethod="double", row.names=FALSE)

##==============================================================================
## End
##==============================================================================
