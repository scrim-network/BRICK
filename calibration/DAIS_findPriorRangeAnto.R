##==============================================================================
## Scratch work -- finding suitable prior ranges for anto.a and anto.b, the
## parameters for fitting functional form for Toc, given Tg
##
## Toc            Antarctic ocean surface temperature (deg C)
## Tg (Tg.recon)  (reconstructed and regressed) mean global surface temperature anomalies (deg C)
## Ta             (reconstructed) Antarctic atmospheric temperature, reduced to sea level (deg C)
## Tf             freezing temperature of sea water (-1.8 deg C)
## anto fit:
##     c  = (Tf-b)/a
##    Toc = Tf + (a*Tg + b-Tf) / (1+exp(-Tg+c))
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

library(DEoptim)

source('../dais/DAIS_readData.R')
Ta = TA     # Antarctic temperature
Toc = TO    # ocean temperature
SL = SL     # sea level

##==============================================================================
## Set up (pre-)calibration windows around the data
## The point: We run the model at many parameter values and see which ones
##            send the simulation through a window (width determined by the
##            observational errors) around the observational data.
## These windows are presented in Shaffer (2014) and Shepherd et al. (2012)
## 1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
## We want the cumulative sea-level equivalent in meters for the year 2002
## Note the conversion: 360Gt = 1mm SLE
estimate.SLE.rate = abs(-71/360)/1000
time.years = 2002-1992      # using the midpoint of the 19-year interval
mid.cum.SLE_2002 = estimate.SLE.rate*time.years
i1992 = which(date==-8)

estimate.SLE.rate.error = abs(-53/360)/1000     #1-sigma error
estimate.SLE.error = sqrt(time.years)*abs(-53/360)/1000 #1-sigma error
        # (*sqrt(10) because 10 years of potentially accumulated error:
        #  total error^2 = year 1 error^2 + year 2 error^2 + ... year 10 error^2
        #                = 10*year X error^2)
SE2_2002 = estimate.SLE.error*2 #2-sigma error

positive_2SE = mid.cum.SLE_2002 + SE2_2002 # Add the 2 standard error to the mean value
negative_2SE = mid.cum.SLE_2002 - SE2_2002 # Subtract the 2 standard error to the mean value

l.usetrend=FALSE
if(l.usetrend) {
  upper.wind = c(6.0, -6.9, -1.25, estimate.SLE.rate+2*estimate.SLE.rate.error) # Kelsey's, modified to match
  lower.wind = c(1.8, -15.8, -4.0, estimate.SLE.rate-2*estimate.SLE.rate.error) # rate in 1992-2012
} else {
  upper.wind = c(6.0, -6.9, -1.25, positive_2SE) # Original from Kelsey
  lower.wind = c(1.8, -15.8, -4.0, negative_2SE)
  #upper.wind = c(5.5, -8, -2, positive_2SE) # From Shaffer 2014, p 1809
  #lower.wind = c(2.5, -17, -4, negative_2SE)
}

windows = matrix(c(lower.wind, upper.wind), nrow = 4, ncol=2)
obs.targets = (windows[,2]+windows[,1])*.5    # middle of window = obs to compare model to
obs.err = (windows[,2]-windows[,1])*.5       # half-width of window = uncertainty
obs.err=0.5*obs.err                           # assume all windows are 2*stdErr
                                              # (last one actually is)
## Create a vector with each observation year
## 120kyr, 20Kyr, 6kyr (before present), and 2002
obs.years = c(120000, 220000, 234000, 240002)
##==============================================================================

##==============================================================================
## Linear regression between reconstructed Antarctic surface temperatures (Ta)
## and global surface temperature anomalies (obs.temp). This will yield a
## lengthy (240,000 years before present - present) time series of global average
## surface temperature anomalies, akin to how DAIS will be forced in the coupled
## BRICK model.
dat = read.table("../doeclim/HadCRUT.4.4.0.0.annual_ns_avg.txt")
Tg = dat[,2]-mean(dat[1:20,2])
Tg.time = dat[,1]
i1997 = which(Tg.time==1997) # end fit at 1997 because Ta is obs after that
ibeg=which(date==-150); iend=ibeg+length(Tg.time[1:i1997])-1;
fit = lm(Tg[1:i1997] ~ Ta[ibeg:iend])
Tg.recon = fit$coefficients[1] + fit$coefficients[2]*Ta
      # Tg.recon is now an appropriately scaled global surface temperature anomaly
      # akin to output from DOECLIM.

## Set up parameter names, ranges, ...
## Set Best Case (Case #4) from Shaffer (2014)
parnames    = c('anto.a','anto.b','gamma','alpha','mu' ,'nu'  ,'P0' ,'kappa','f0' ,'h0'  ,'c' , 'b0','slope')
parameters0 = c(  0.26  ,  -18   ,  2    , 0.35  , 8.7 , 0.012, 0.35, 0.04  , 1.2 , 1471 , 95 , 775 , 0.0006)
##==============================================================================

##==============================================================================
## Set the upper and lower bounds, and revise from Shaffer, 2014
bound.lower=parameters0-(parameters0*0.5);  bound.upper=parameters0 + (parameters0*0.5);
bound.lower[match('gamma',parnames)]=0.5;   bound.upper[match('gamma',parnames)]=17/4;
bound.lower[match('alpha',parnames)]=0;     bound.upper[match('alpha',parnames)]=1;
bound.lower[match('b0',parnames)]=725;      bound.upper[match('b0',parnames)]=825;
bound.lower[match('slope',parnames)]=4.5e-4;bound.upper[match('slope',parnames)]=7.5e-4;
bound.lower[match('anto.a',parnames)]=0;    bound.upper[match('anto.a',parnames)]=1;
bound.lower[match('anto.b',parnames)]=-100;    bound.upper[match('anto.b',parnames)]=50;
## Try doubling the range for priors, not for gamma and alpha, since those are
## from Shaffer 2014
widths=bound.upper-bound.lower
bound.lower[c(1,2,5:length(bound.lower))] = bound.lower[c(1,2,5:length(bound.lower))] -
                                        0.5*widths[c(1,2,5:length(bound.lower))]
bound.upper[c(1,2,5:length(bound.lower))] = bound.upper[c(1,2,5:length(bound.lower))] +
                                        0.5*widths[c(1,2,5:length(bound.lower))]
##==============================================================================

## Fetch the DAIS model
source('../dais/daisantoF.R')

##==============================================================================
## Define a function to optimize
rmse.dais <- function( parameters.in, parnames.in,
                        obs.in, obs.step.in, obs.err.in, Tg.in)
{

  anto.a=parameters.in[match("anto.a",parnames.in)]
  anto.b=parameters.in[match("anto.b",parnames.in)]
  gamma =parameters.in[match("gamma" ,parnames.in)]
  alpha =parameters.in[match("alpha" ,parnames.in)]
  mu =parameters.in[match("mu" ,parnames.in)]
  nu =parameters.in[match("nu" ,parnames.in)]
  P0 =parameters.in[match("P0" ,parnames.in)]
  kappa =parameters.in[match("kappa" ,parnames.in)]
  f0 =parameters.in[match("f0" ,parnames.in)]
  h0 =parameters.in[match("h0" ,parnames.in)]
  c =parameters.in[match("c" ,parnames.in)]
  b0 =parameters.in[match("b0" ,parnames.in)]
  slope =parameters.in[match("slope" ,parnames.in)]
  sigma =parameters.in[match("sigma" ,parnames.in)]

  dais.out = daisantoF(anto.a=anto.a, anto.b=anto.b,
                       gamma=gamma  , alpha=alpha  ,
                       mu=mu        , nu=nu        ,
                       P0=P0        , kappa=kappa  ,
                       f0=f0        , h0=h0        ,
                       c=c          , b0=b0        ,
                       slope=slope  ,
                       Ta=Tg.in,
                       SL=SL )
  dais.out.norm = dais.out - mean(dais.out[SL.1961_1990])
  rmse=0
  #get the residuals
  resid = c( obs.in[1] - dais.out.norm[obs.step.in[1]] ,
             obs.in[2] - dais.out.norm[obs.step.in[2]] ,
             obs.in[3] - dais.out.norm[obs.step.in[3]] ,
             obs.in[4] - (dais.out.norm[obs.step.in[4]]-dais.out.norm[i1992]) )
  rmse = sqrt( mean(resid^2) )
  return(rmse)
}
##==============================================================================
## Define another function to optimize, only anto parameters
rmse.anto <- function( parameters.in, parnames.in,
                        obs.in, obs.step.in, obs.err.in, Tg.in)
{

  anto.a=parameters.in[match("anto.a",parnames.in)]
  anto.b=parameters.in[match("anto.b",parnames.in)]

  rmse=0

  #get the residuals
if(FALSE){  #model performance?
  dais.out = daisantoF(anto.a=anto.a, anto.b=anto.b,
                       Ta=Tg.in,
                       SL=SL )
  dais.out.norm = dais.out - mean(dais.out[SL.1961_1990])
  resid = c( obs.in[1] - dais.out.norm[obs.step.in[1]] ,
             obs.in[2] - dais.out.norm[obs.step.in[2]] ,
             obs.in[3] - dais.out.norm[obs.step.in[3]] ,
             obs.in[4] - (dais.out.norm[obs.step.in[4]]-dais.out.norm[i1992]) )
} else {  # matching actual Toc
  Toc.fit = anto(anto.a, anto.b, Tg=Tg.in)
  resid = Toc.fit-Toc
}
  rmse = sqrt( mean(resid^2) )
  return(rmse)
}
##==============================================================================

##==============================================================================
## Use DEoptim to obtain a minimum RMSE over the parameters
niter.deoptim=10										# number of iterations for DE optimization
NP.deoptim=10*length(parnames)			# population size for DEoptim (do at least 10*[N parameters])
F.deoptim=0.8													# as suggested by Storn et al (2006)
CR.deoptim=0.9												# as suggested by Storn et al (2006)
outDEoptim <- DEoptim(rmse.dais, bound.lower, bound.upper,
				DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
				parnames.in=parnames, obs.in=obs.targets, obs.step.in=obs.years, obs.err.in=obs.errs, Tg.in=Ta)
p0 = outDEoptim$optim$bestmem


outDEoptim <- DEoptim(rmse.anto, bound.lower[1:2], bound.upper[1:2],
				DEoptim.control(NP=10*2,itermax=50,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
				parnames.in=parnames[1:2], obs.in=obs.targets, obs.step.in=obs.years, obs.err.in=obs.errs, Tg.in=Ta)
p0 = outDEoptim$optim$bestmem
##==============================================================================

##==============================================================================
## Do an MCMC calibration of just anto.a and anto.b to get an idea of prior ranges
##==============
library(adaptMCMC)										# use robust adaptive Metropolis

## Define a function to *maximize* in the MCMC
rmse.anto.neg <- function( parameters.in, parnames.in,
                        obs.in, obs.step.in, obs.err.in, Tg.in,
                        bound.lower.in, bound.upper.in)
{

if(all(parameters.in > bound.lower.in & parameters.in < bound.upper.in)) {

  anto.a=parameters.in[match("anto.a",parnames.in)]
  anto.b=parameters.in[match("anto.b",parnames.in)]

  rmse=0

  #get the residuals
  if(FALSE){  #model performance?
    dais.out = daisantoF(anto.a=anto.a, anto.b=anto.b,
                          Ta=Tg.in, SL=SL )
    dais.out.norm = dais.out - mean(dais.out[SL.1961_1990])
    resid = c( obs.in[1] - dais.out.norm[obs.step.in[1]] ,
               obs.in[2] - dais.out.norm[obs.step.in[2]] ,
               obs.in[3] - dais.out.norm[obs.step.in[3]] ,
               obs.in[4] - (dais.out.norm[obs.step.in[4]]-dais.out.norm[i1992]) )
  } else {  # matching actual Toc
    Toc.fit = anto(anto.a, anto.b, Tg=Tg.in)
    resid = Toc.fit-Toc
  }
#  rmse.neg = -1*sqrt( mean(resid^2) )
  rmse.neg = sum (dnorm(resid, mean=rep(0,length(Toc.fit)), sd = .5, log=TRUE))
} else {
  rmse.neg=-Inf
}
  return(rmse.neg)
}
##==============
niter.mcmc = 1e4
step.mcmc = c(0.05, 0.05)
accept.mcmc = 0.234
gamma.mcmc = 0.5
burnin = round(niter.mcmc*.5)

t.beg=proc.time()											# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out1 = MCMC(rmse.anto.neg, niter.mcmc, p0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
				parnames.in=parnames[1:2], obs.in=obs.targets, obs.step.in=obs.years, obs.err.in=obs.errs, Tg.in=Ta,
        bound.lower.in=bound.lower[1:2], bound.upper.in=bound.upper[1:2])
t.end=proc.time()											# save timing
chain1 = amcmc.out1$samples

## Check #1: History plots
par(mfrow=c(2,1))
for (pp in 1:length(parnames)) {
	plot(chain1[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="")
}

par(mfrow=c(2,1))
for (pp in 1:length(parnames)) {
	hist(chain1[burnin:niter.mcmc,pp], xlab=parnames[pp])
}
##==============================================================================
## Should find:
p0.Ta = c(0.1317,3.097)   # calibrated anto using Ta
p0.Tg = c(0.1574,0.6677)  # calibrated anto using Tg.recon

dais.out = daisF( Ta=Ta, Toc=Toc, SL=SL)
daisanto.out.Ta = daisantoF( anto.a=p0.Ta[1], anto.b=p0.Ta[2],
                             Ta=Ta, SL=SL )
daisanto.out.Tg = daisantoF( anto.a=p0.Tg[1], anto.b=p0.Tg[2],
                             Ta=Tg.recon, SL=SL )

## daisanto.out.Tg gives us an idea of how wrong we might be if we calibrate based
## on incorrect Antarctic surface temperature.

## Based on these parameters (p0.Ta, p0.Tg), see that we want to let
##   anto.a range from 0 to 1
##   anto.b range from 0 (no negative?) to at least 4

##==============================================================================
## End
##==============================================================================
