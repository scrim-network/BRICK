##==============================================================================
## Script to precalibrate TE (thermosteric expansion) model for use within the
## BRICK coupled modeling framework.
##
## Parameters to calibrate:
## [1] a.te       # sensitivity of TE equilibrium [m/K]
## [2] b.te       # equilibrium TE [m] for temperature anomaly Tg = 0
## [3] invtau.te  # 1/time-scale of exponential decay (e-folding time) [1/years]
## [4] TE_0       # initial sea level
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

niter.deoptim=100		   	# number of iterations for DE optimization
NP.deoptim=20	      		# population size for DEoptim (do at least 10*[N parameters])
F.deoptim=0.8						# as suggested by Storn et al (2006)
CR.deoptim=0.9					# as suggested by Storn et al (2006)
outDEoptim <- DEoptim(rmse.quantiles, c(0,0), c(100,100),
				DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
				q05.in=q05, q95.in=q95)
shape.invtau = outDEoptim$optim$bestmem[1]
scale.invtau = outDEoptim$optim$bestmem[2]

## This should yield results somewhere in the ballpark of
##
##     shape.invtau = 1.81    and     scale.invtau = 0.00275
##
## If yours does not, check whether the gamma distribution defined by your
## shape and scale fit the invtau.hat=1/200, q05=1/1290 and q95=1/82
## requirements above.

##==============================================================================

## Setup packages and libraries
#install.packages('lhs')
#install.packages('compiler')
#install.packages('pscl')
require(lhs)
library(pscl) # install inverse gamma distribution
library(compiler)
enableJIT(3)
enableJIT(3)

## Set the seed (for reproducibility)
set.seed(1234)

## Set up the model simulation and read the data forcing for hindcasts and precalibration.

begyear = 1850
endyear = 2009; #if(!l.project & endyear>2009) print('l.project and endyear not compatible')
mod.time= begyear:endyear
begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)

source('../R/compute_indices.R')
source('../calibration/TE_readData.R')
source('../calibration/DOECLIM_readData.R')

ibeg = which(obs.temp.time==begyear); iend = which(obs.temp.time==endyear)
if(length(ibeg)*length(iend)==0) print('ERROR - begyear/endyear outside forcing data')

forcing.time = obs.temp.time[ibeg:iend]
forcing.temp = obs.temp[ibeg:iend]

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.sl)
names(midx.all) = c(   "sl"   )
oidx.all        = list(oidx.sl)
names(oidx.all) = c(   "sl"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.sl)
names(obs.all) = c(    "sl"  )
obs.err.all        = list( obs.sl.err)
names(obs.err.all) = c(    "sl"      )

## What indices should be used to normalize in same way as data?
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1961),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
		c(which(mod.time==1990),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )


## Set parameter names
parnames    = c('a.te','b.te','invtau.te','TE0')

## Set the upper and lower bounds
## 24.2 mm SLR is uncertainty in first annual data point (1880) from Church and
## White (2011) SLR data, so is appropriate 1-sigma for TE0. Also happens to be
## the largest uncertainty associated with those data, so it is appropriate to
## use as 1-sigma for the b parameter.
    ## >>>>>>>> SAMPLING 1/TAU <<<<<<<<<
p0		      =c(0.3616, 0.5   , 1/200      , 0.0   )        # initial parameter guesses
bound.lower=c(0.0   , 0.0   , 0          ,-0.0484)        # prior range lower bounds
bound.upper=c(0.8595, 2.193 , 1          , 0.0484)        # prior range upper bounds
    ## >>>>>>>> SAMPLING TAU <<<<<<<<<
#p0		     =c( 0.3616  , 0.186 , 200       ,  0.0        )  # initial parameter guesses
#bound.lower=c( 0.0     , 0     , 14.89     , -2*24.2/1000)  # prior range lower bounds
#bound.upper=c( 0.8595  , 2.193 , 1357.11   ,  2*24.2/1000)  # prior range upper bounds

## effectively get rid of b and TE0?
#bound.lower[c(2,4)] = bound.lower[c(2,4)]*1e-10
#bound.upper[c(2,4)] = bound.upper[c(2,4)]*1e-10

## effectively get rid of TE0?
#bound.lower[c(4)] = bound.lower[c(4)]*1e-10
#bound.upper[c(4)] = bound.upper[c(4)]*1e-10

step.mcmc	 =0.05*(bound.upper-bound.lower) # set size for parameters in MCMC (proposals)
index.model=c(1,2,3,4)			# which are model parameters? (index within parnames)

## Grab the TE model
source('../fortran/R/brick_te_F.R')

##==============================================================================
## Get a standard hindcast and future projections (a "control" model)
te.out = brick_te_F(a   =p0[match("a.te"  ,parnames)],
                  b   =p0[match("b.te"  ,parnames)],
                  invtau =p0[match("invtau.te",parnames)],
                  TE_0=p0[match("TE0"   ,parnames)],
                  Tg  =obs.temp                   )

## Relative to 1961-1990 mean, just like temperature and sea level data
#te.out = te.out - mean(te.out[i1961_1990])
te.out = te.out - mean(te.out[ind.norm])
##==============================================================================

##==============================================================================
## Set up (pre-)calibration windows around the data
## The point: We run the model at many parameter values and see which ones
##            send the simulation through a window (width determined by the
##            observational errors) around the observational data.
## The data we precalibrate against are the total changes in global average sea
## level (Church and White (2011)).
## If sea level rises, then the contribution from thermosteric expansion should
## not exceed the total sea level rise. If it does, then this parameter
## combination is implausible.
## If sea level drops, then cannot say anything. This is because the other
## factors (ice sheets, glaciers, etc) may be raising sea level, but unsure
## about TE contribution.

##==============================================================================
## MCMC simulation
##==============================================================================

## Source the statistical models
source('../calibration/TE_assimLikelihood.R')

##library(mcmc)												# use vanilla Metropolis
library(adaptMCMC)										# use robust adaptive Metropolis
accept.mcmc = 0.234										# Optimal as # parameters->infinity
																			#	(Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = 3e5											# number of iterations for MCMC
gamma.mcmc = 0.5											# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)				# remove first ?? of chains for burn-in
stopadapt.mcmc = round(niter.mcmc*1.0)# stop adapting after ?? iterations? (niter*1 => don't stop)
step.mcmc = (bound.upper-bound.lower)*.05
##==============================================================================

## Actually run the calibration
t.beg=proc.time()											# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out1 = MCMC(log.post, niter.mcmc, p0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames        , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  forcing.temp.in=forcing.temp, obs=obs.all               , obs.err=obs.err.all       ,
                  trends.te=trends.te         , oidx=oidx.all             , midx=midx.all             ,
                  shape.in=shape.invtau       , scale.in=scale.invtau     , ind.norm=ind.norm
                  )
t.end=proc.time()											# save timing
chain1 = amcmc.out1$samples

## If you want to run 2 (or more) chains in parallel (save time, more sampling)
if(FALSE){
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.par1 = MCMC.parallel(log.post, niter.mcmc, p0, n.chain=4, n.cpu=4, dyn.libs='../fortran/brick_te.so',
									scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames        , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  forcing.temp.in=forcing.temp, obs=obs.all               , obs.err=obs.err.all       ,
                  trends.te=trends.te         , oidx=oidx.all             , midx=midx.all             ,
                  shape.in=shape.invtau       , scale.in=scale.invtau     , ind.norm=ind.norm
                  )
t.end=proc.time()											# save timing
chain1=amcmc.par1[[1]]$samples
chain2=amcmc.par1[[2]]$samples
chain3=amcmc.par1[[3]]$samples
chain4=amcmc.par1[[4]]$samples
}

## Extend and run more MCMC samples?
if(FALSE){
t.beg=proc.time()
amcmc.extend1 = MCMC.add.samples(amcmc.out1, niter.mcmc,
									parnames.in=parnames        , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  forcing.temp.in=forcing.temp, obs=obs.all               , obs.err=obs.err.all       ,
                  trends.te=trends.te         , oidx=oidx.all             , midx=midx.all             ,
                  shape.in=shape.invtau       , scale.in=scale.invtau     , ind.norm=ind.norm
                  )
t.end=proc.time()											# save timing
chain1 = amcmc.extend1$samples
}

##==============================================================================

## Check #1: History plots
par(mfrow=c(2,2))
for (pp in 1:length(parnames)) {
	plot(chain1[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="")
}

## Check #2: Heidelberger and Welch's convergence diagnostic:
heidel.diag(chain1, eps=0.1, pvalue=0.05)

## Check #3: Gelman and Rubin's convergence diagnostic:
# Converged when the potental scale reduction factor is less than 1.1
set.seed(222)
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out2 = MCMC(log.post, niter.mcmc, p0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames    , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  sl.in=obs.sl            , sl.err.in=obs.sl.std      , trends.in=trends          ,
                  ind.trends.in=ind.trends, err.trends.in=err.trends  , shape.in=shape.invtau     ,
                  scale.in=scale.invtau   , temp.in=obs.temp )
chain2 = amcmc.out2$samples
t.end=proc.time()											# save timing

mcmc1 = as.mcmc(chain1)
mcmc2 = as.mcmc(chain2)

mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2))
gelman.diag(mcmc_chain_list)

# Reset seed back to original
set.seed(1234)

## Histograms
par(mfrow=c(2,2))
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
## TE is converged within 10,000s-100,000 of iterations.
n.burnin = 100000
n.sample = nrow(chain1)-n.burnin
parameters1=chain1[(n.burnin+1):nrow(chain1),]
parameters2=chain2[(n.burnin+1):nrow(chain1),]
parameters3=chain3[(n.burnin+1):nrow(chain1),]
parameters4=chain4[(n.burnin+1):nrow(chain1),]
parameters.posterior = rbind(parameters1,parameters2, parameters3, parameters4)
n.parameters = ncol(parameters.posterior)

## Histograms
par(mfrow=c(2,2))
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
par(mfrow=c(2,2))
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
filename=paste('../output_calibration/TE_precalibrationMCMC_parameters_',today,'.csv', sep="")
write.table(to.file, file=filename, sep=",", qmethod="double", row.names=FALSE)

##==============================================================================
## End
##==============================================================================






##==============================================================================
##==============================================================================
## Further testing of sampling 1/tau (invtau) versus tau itself

## get TE simulation results (cumulative SLR relative to 1880 (??)) by sampling 1/tau
te_invtau = mat.or.vec(dim(p_invtau)[1],length(te.out))
for(i in 1:dim(p_invtau)[1]) {
  # Run the TE model at each iteration of latin hypercube sample parameters
  te.out = brick_te(a=p_invtau[i,ind.a],     b=p_invtau[i,ind.b],
                    tau=p_invtau[i,ind.tau], TE_0=p_invtau[i,ind.TE0],
                    Tg=obs.temp )
  # Normalize to the 1961-1990 mean
#  te.out = te.out - mean(te.out[i1961_1990])
  te.out = te.out - mean(te.out[1:20])
  # Save
  te_invtau[i,] = te.out
}

## get TE simulation results (cumulative SLR relative to 1880 (??)) by sampling tau
te_tau = mat.or.vec(dim(p_tau)[1],length(te.out))
for(i in 1:dim(p_tau)[1]) {
  # Run the TE model at each iteration of latin hypercube sample parameters
  te.out = brick_te(a=p_tau[i,ind.a],     b=p_tau[i,ind.b],
                    tau=p_tau[i,ind.tau], TE_0=p_tau[i,ind.TE0],
                    Tg=obs.temp )
  # Normalize to the 1961-1990 mean
#  te.out = te.out - mean(te.out[i1961_1990])
  te.out = te.out - mean(te.out[1:20])
  # Save
  te_tau[i,] = te.out
}

if(FALSE){ # commands to save because you can re-use them
  te_noTE0nob = te_tau
  te_noTE0 = te_tau
}

## Plot the trajectories and histograms of the last time slice
par(mfrow=c(2,2))
plot(obs.time,te_invtau[1,],type='l'); for (i in 2:dim(p_invtau)[1]){lines(obs.time,te_invtau[i,],type='l')}
plot(obs.time,te_tau[1,],type='l'); for (i in 2:dim(p_tau)[1]){lines(obs.time,te_tau[i,],type='l')}
hist(te_invtau[,length(te.out)])
hist(te_tau[,length(te.out)])

## Plot of with all parameters, no TE0, no TE0+no b
par(mfrow=c(2,2))
plot(obs.time,te_invtau[1,],type='l',ylab="TE SLR (m)",xlab="year",main="all parameters",ylim=c(-0.05,0.05))
  for (i in 2:dim(p_invtau)[1]){lines(obs.time,te_invtau[i,],type='l')}
plot(obs.time,te_noTE0[1,],type='l',ylab="TE SLR (m)",xlab="year",main="no TE0",ylim=c(-0.05,0.05))
  for (i in 2:dim(p_noTE0)[1]){lines(obs.time,te_noTE0[i,],type='l')}
plot(obs.time,te_noTE0nob[1,],type='l',ylab="TE SLR (m)",xlab="year",main="no b, no TE0",ylim=c(-0.05,0.05))
  for (i in 2:dim(p_noTE0nob)[1]){lines(obs.time,te_noTE0nob[i,],type='l')}

##==============================================================================
##==============================================================================
