##==============================================================================
## Script to precalibrate TE (thermosteric expansion) model for use within the
## BRICK coupled modeling framework.
##
## Parameters to calibrate:
## [1] a.te       # sensitivity of ice flow to sea level
## [2] b.te       # sensitivity of ice flow to ocean subsurface temperature
## [3] tau.te     # Profile parameter related to ice stress [m^(1/2)]
## [4] TE_0       # Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
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

## Read the data forcing for hindcasts and precalibration. Yields:
## obs.time
## obs.temp
## obs.sl
## obs.sl.std
source('../brick_te/TE_readData.R')

## What are the observed sea level rises? What indices do they occur? This will
## serve as the precalibration data (described below)
obs.slr = diff(obs.sl)
irise = which(obs.slr>0)

## Set parameter names
parnames    = c('a.te','b.te','tau.te','TE0')

## Set the upper and lower bounds
## 24.2 mm SLR is uncertainty in first annual data point (1880) from Church and
## White (2011) SLR data, so is appropriate 1-sigma for TE0. Also happens to be
## the largest uncertainty associated with those data, so it is appropriate to
## use as 1-sigma for the b parameter.
    ## >>>>>>>> SAMPLING 1/TAU <<<<<<<<<
p0		     =c(0.3616 ,    0.186       , invtau.hat  ,  0.0        )  # initial parameter guesses
bound.lower=c(0.0      , 0 , 0 , -2*24.2/1000)  # prior range lower bounds
bound.upper=c(0.8595   ,  2.193 , qgamma(.999,shape=shape.invtau,scale=scale.invtau) ,  2*24.2/1000)  # prior range upper bounds
    ## >>>>>>>> SAMPLING TAU <<<<<<<<<
#p0		     =c(0.045878 ,    0.0       , 65  ,  0.0        )  # initial parameter guesses
#bound.lower=c(0.0045878, -2*24.2/1000 , 30 , -2*24.2/1000)  # prior range lower bounds
#bound.upper=c(0.45878  ,  2*24.2/1000 , 650 ,  2*24.2/1000)  # prior range upper bounds

## effectively get rid of b and TE0?
bound.lower[c(2,4)] = bound.lower[c(2,4)]*1e-10
bound.upper[c(2,4)] = bound.upper[c(2,4)]*1e-10

## effectively get rid of TE0?
bound.lower[c(4)] = bound.lower[c(4)]*1e-10
bound.upper[c(4)] = bound.upper[c(4)]*1e-10

step.mcmc	 =0.05*(bound.upper-bound.lower) # set size for parameters in MCMC (proposals)
index.model=c(1,2,3,4)			# which are model parameters? (index within parnames)

## Grab the TE model
source('../brick_te/brick_te.R')

##==============================================================================
## Get a standard hindcast and future projections (a "control" model)
te.out = brick_te(a   =p0[match("a.te"  ,parnames)],
                  b   =p0[match("b.te"  ,parnames)],
                  tau =p0[match("tau.te",parnames)],
                  TE_0=p0[match("TE0"   ,parnames)],
                  Tg  =obs.temp                   )

## Relative to 1961-1990 mean, just like temperature and sea level data
#te.out = te.out - mean(te.out[i1961_1990])
te.out = te.out - mean(te.out[1:20])
slr.te = diff(te.out)
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
## Latin hypercube sampling (LHS) for pre-calibration (precal)
## The point: LHS disperses a given (n_samples) number of points throughout the
##            parameter-space. "maximinLHS" is a good function to use because it
##            maximizes the minimum distance between sample points, and we don't
##            want clumping. We will run the model with all these different
##            parameter combinations. The ones which produce model output which
##            match the sets of observations at the dates of observations are
##            considered successes.
n_samples = 2000                 # define the number of samples to take
n_parameters = length(parnames) # number of parameters to optimize over

## Define a LHS sampling function
## Note 1: The calls to qunif are needed because LHS in R returns samples
##         distributed on a unit hypercube.
fn.LHS <- function(n_samples, x, y, bound.upper, bound.lower) {
  x <- maximinLHS(n_samples, n_parameters)
  y <- x
  y[,1] <- qunif(x[,1], bound.lower[1], bound.upper[1])
  y[,2] <- qunif(x[,2], bound.lower[2], bound.upper[2])
#  y[,3] <- qunif(x[,3], bound.lower[3], bound.upper[3])
  y[,3] <- qgamma(x[,3], shape=shape.invtau, scale=scale.invtau)  # sample 1/tau from a gamma
  y[,4] <- qunif(x[,4], bound.lower[4], bound.upper[4])
  return(as.data.frame(y))
}

## Generate a LHS sample from the box parameter-space defined by
## bound.upper and bound.lower
t0=proc.time()
parameters = fn.LHS(n_samples, x, y, bound.upper, bound.lower)
t1=proc.time()
colnames(parameters)=c(parnames)
print(paste(n_samples,' samples required ',round((t1-t0)[3]*10)/10,'s to set up LHS'))

## Vector of 0 or 1 if the simulation produces plausible TE SLR results
survive.all = rep(NULL,n_samples)
survive.slr = rep(NULL,n_samples)
survive.1971 = rep(NULL,n_samples)
survive.1993 = rep(NULL,n_samples)
ind.good= NULL

## Set up the indices. Do it this way to keep the parameters and their names
## most general
ind.a  =match("a.te",parnames);   ind.b  =match("b.te",parnames);
ind.tau=match("tau.te",parnames); ind.TE0=match("TE0",parnames);

## Actually run the LHS. And set up a progress bar
pb <- txtProgressBar(min=0,max=n_samples,initial=0,style=3)
t0=proc.time()                   # save the timing information
for(i in 1:n_samples) {

  # Run the TE model at each iteration of latin hypercube sample parameters
  te.out = brick_te(a=parameters[i,ind.a],     b=parameters[i,ind.b],
                    tau=parameters[i,ind.tau], TE_0=parameters[i,ind.TE0],
                    Tg=obs.temp )

  # Normalize to the 1961-1990 mean
#  te.out = te.out - mean(te.out[i1961_1990])
  te.out = te.out - mean(te.out[1:20])

  # Calculate the sea level rise
  mod.slr = diff(te.out)

  # Check if the modeled TE SLR < observed total SLR in all years in which sea
  # level increased
  survive.slr[i] = all(te.out[20:length(te.out)]<=(obs.sl[20:length(obs.sl)]+0*obs.sl.std[20:length(obs.sl)]))

  # Check if matches IPCC AR5, Ch 13, observed trends:
  # 1971-2010: 0.8 +/- 0.3 mm/year (note that te.out is in meters)
  # 1993-2010: 1.1 +/- 0.3 mm/year
  trend.1971 = lm(te.out[i1971_2010]~i1971_2010)$coefficients[2]
  trend.1993 = lm(te.out[i1993_2010]~i1993_2010)$coefficients[2]
  survive.1971[i]=abs(trend.1971*1000-0.8) <= 0.3*2
  survive.1993[i]=abs(trend.1993*1000-1.1) <= 0.3*2

  # Update the progress bar
  setTxtProgressBar(pb, i)

}
t1=proc.time()                   # save the timing information
close(pb)
print(paste(n_samples,' samples required ',round((t1-t0)[3]*10)/10,'s to run LHS'))

## What are the indices of the good parameter combinations?
ind.good = which(survive.slr & survive.1971 & survive.1993)

#p_tau = parameters[ind.good,]
#p_invtau = parameters[ind.good,]

## Plot some histograms
par(mfrow=c(2,2))
for (pp in 1:length(parnames)){
  hist(parameters[ind.good,pp] , xlim=c(bound.lower[pp],bound.upper[pp]) ,
                                 xlab=parnames[pp] )
}

##==============================================================================
## End
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
## End
##==============================================================================



##
## (the rest of this stuff is leftover from copyign DAIS precalibration; it's
## code left there because one might want to use some of it)
##



## Save workspace image - you do not want to re-simulate all those!
save.image(file = "DAIS_precalibration_LHS.RData")

## Make sure this is in SLE volume (m), not m^3 volume
if(dais.models.lhs[1,1]>1e5) {
  for (i in 1:n_samples){
    dais.models.lhs[i,] = 57*(1-dais.models.lhs[i,]/Volo)
  }
}
save.image(file = "DAIS_precalibration_LHS.RData")

## Apply the biases (remember, you sampled the variance, so square-root it!)
##    True world = model + bias   (for a good read, check out Higdon et al, 2004)
## Model bias is assumed constant over entire model simulation.
## While you are at it... Apply scaling to get modeled Antarctic volume (DAIS
## output) in SLE (for BRICK). The 57 is from Shaffer (2014)
## Subtract off the 1961-1990 mean
dais.models.lhs.anom = mat.or.vec(n_samples, length(AIS_melt))
for (i in 1:n_samples){
  dais.models.lhs.anom[i,] = dais.models.lhs[i,] - mean(dais.models.lhs[i,SL.1961_1990])
}

## Memory efficiency with large sample
rm(dais.models.lhs) # done with this one, since we have xxx.wbais version
save.image(file = "DAIS_precalibration_LHS.RData")

## Finally, add back in the noise
bias = sqrt(parameters[,12])
dais.precal = mat.or.vec(n_samples,length(AIS_melt))
for(i in 1:n_samples){
  dais.precal[i,] = dais.models.lhs.anom[i,] + rnorm(1,mean=0,sd=bias[i])
}

## Memory efficiency with large sample
rm(dais.models.lhs.wbias) # done with this one, since we have xxx.precal version
save.image(file = "DAIS_precalibration_LHS.RData")

## Write csv of SLE values for the targeted years
surv.targ = matrix(c(dais.precal[1:n_samples,120000],
                     dais.precal[1:n_samples,220000],
                     dais.precal[1:n_samples,234000],
                     dais.precal[1:n_samples,240002]), nrow=n_samples, ncol=4)
colnames(surv.targ, do.NULL = FALSE)
colnames(surv.targ) <- c("Last Interglacial", "Last Glacial Max", "Holocene", "93-2011 Trend")
write.csv(surv.targ, file="surviving_targets_LHS.csv")

## Define function to check if model run passes through the precalibration
## windows around the data (original code from Kelsey Ruckert <klr324@psu.edu>)
surviveTarget = function(window, surtarget){
  lbw = window[1]
  upw = window[2]
  sv = surtarget
  which((sv >= lbw & sv <= upw))
}

## Check which simulations fit through each precalibration window around the data
surLIG = surviveTarget(windows[1,], surv.targ[,1])
surLGM = surviveTarget(windows[2,], surv.targ[,2])
surMH = surviveTarget(windows[3,], surv.targ[,3])
sur9311trend = surviveTarget(windows[4,], surv.targ[,4])

## Find the runs that pass through all constraints
total.sur = c(surLIG,surLGM,surMH,sur9311trend)
template <- table(as.vector(total.sur))
all = names(template)[template == max(template)]
sur.all = as.numeric(all)

## How many simulations fit each constraint? Which one is the toughest hoop to jump through?
percent.include = c((n_samples/n_samples)*100, (length(surLIG)/n_samples)*100,
                    (length(surLGM)/n_samples)*100, (length(surMH)/n_samples)*100,
                    (length(sur9311trend)/n_samples)*100)
constraints = c("No constraints","Last integlacial","Last glacial maximum","Mid-Holocene","Instrumental period")
table.parts = matrix(c(constraints, percent.include), nrow=5, ncol=2)
write.csv(table.parts, file="constraint_trend_percent.csv")

## Fit PDFs to the parameter distributions
## Just marginals, but turns out an 11-dimensional joint KDE is a real doozie...
pdf.all=vector('list',12)
for (pp in 1:12){
  if(pp!=12) tmp = density(parameters[sur.all,pp],kernel='gaussian',
                          n=100,from=bound.lower[pp],to=bound.upper[pp])
  if(pp==12) tmp = density(parameters[sur.all,pp],kernel='gaussian', n=100)
  pdf.all[[pp]] = tmp; names(pdf.all)[pp]=colnames(parameters)[pp]
}

## Write a CSV file with the successful parameter combinations and bandwidths
bandwidths=rep(NA,length(colnames(parameters)))
for (i in 1:length(colnames(parameters))){
  bandwidths[i]=pdf.all[[i]]$bw
}
write.csv(rbind(parameters[sur.all,],bandwidths), file="successful_parameters_2.csv")

## Save workspace image - you do not want to re-process all those!
save.image(file = "DAIS_precalibration_LHS.RData")

## Plot the PDFs
par(mfrow=c(4,3))
for (pp in 1:length(parnames)) {
	if(pp!=2) plot(pdf.all[[pp]]$x,pdf.all[[pp]]$y, type="l", ylab="density", xlab=parnames[pp],
                 main="", xlim=c(bound.lower[pp],bound.upper[pp]))
  if(pp==2) plot(pdf.all[[pp]]$x,pdf.all[[pp]]$y, type="l", ylab="density", xlab=parnames[pp],
                 main="Gaussian kernels", xlim=c(bound.lower[pp],bound.upper[pp]))
}

##==============================================================================


##==============================================================================
## End
##==============================================================================
