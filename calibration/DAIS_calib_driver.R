##==============================================================================
## Script to calibrate DAIS model for use within BRICK coupled modeling
## framework.
##
## Parameters to calibrate, and standard values (Shaffer, 2014):
## [1] gamma = 1/2-17/4  # sensitivity of ice flow to sea level
## [2] alpha = 0-1       # sensitivity of ice flow to ocean subsurface temperature
## [3] mu = 8.7          # Profile parameter related to ice stress [m^(1/2)]
## [4] nu = 1.2e-2       # Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
## [5] P0 = 0.35         # Precipitation at 0C [m of ice/yr]
## [6] kappa = 4e-2      # Relates precipitation to temperature [K^-1]
## [7] f0 = 1.2          # Constant of proportionality for ice speed [m/yr]
## [8] h0 = 1471         # Initial value for runoff line calculation [m]
## [9] c = 95            # Second value for runoff line calculation [m]
## [10] b0 = 775         # Height of bed at the center of the continent [m]
## [11] slope = 6e-4     # Slope of the bed
##
## Much of the original code is thanks to Kelsey Ruckert <klr324@psu.edu>
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

## Switch to your BRICK calibration directory
setwd('/home/scrim/axw322/codes/BRICK/calibration')
#setwd('/Users/tony/codes/BRICK/calibration')

## Set up MCMC stuff here so that it can be automated for HPC
nnode_mcmc000 <- 8
niter_mcmc000 <- 1e6

## Set up a filename for saving RData images along the way
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.saveprogress <- paste('DAIS_calib_MCMC_',today,'.RData',sep='')

## Use the AIS fast dynamics emulator (Wong et al., 2017)
l.aisfastdy <- TRUE
fd.priors <- 'gamma'			# which prior distirbutions to use?

## Set the seed (for reproducibility)
set.seed(1234)

## Setup packages and libraries
#install.packages('compiler')
#install.packages('pscl')
library(pscl) # install inverse gamma distribution
library(compiler)
enableJIT(3)
enableJIT(3)

## Read the data forcing for hindcasts and projections. Yields:
##  Ta (Antarctic temperature reduced to sea-level)
##  Toc (High latitude subsurface ocean temperature)
##  SL, obs.sl (Reconstructed sea level, Church and White (2011) modern-era sea level)
##  dSL (time rate of change of sea level)
##  date (240 Kyr before present to 2100 AD at 1-year intervals of forcings)
l.project=FALSE
source('../calibration/DAIS_readData.R')

## Set Best Case (Case #4) from Shaffer (2014), and parameter ranges
    # >>> with anto? <<<
if(l.aisfastdy) {
  parnames    = c('anto.a','anto.b','gamma','alpha.dais','mu'  ,'nu'  ,'P0' ,'kappa.dais','f0' ,'h0'  ,'c'  , 'b0','slope' ,'lambda','Tcrit','var.dais')
  parameters0 = c(  0.1574, 0.6677 ,  2    , 0.35       , 8.7  , 0.012, 0.35, 0.04       , 1.2 , 1471 , 95  , 775 , 0.0006 , 0.01   , -15   , 0.0004656)
  bound.lower = c( 0.0    , 0      ,  0.5  , 0          , 7.05 , 0.003,0.026, 0.025      , 0.6 , 735.5, 47.5, 740 , 0.00045, 0.005  , -20   , 0        )
  bound.upper = c( 1.0    , 2      ,  4.25 , 1          , 13.65, 0.015, 1.5 , 0.085      , 1.8 ,2206.5,142.5, 820 , 0.00075, 0.015  , -10   , 2        )
} else {
  parnames    = c('anto.a','anto.b','gamma','alpha.dais','mu'  ,'nu'  ,'P0' ,'kappa.dais','f0' ,'h0'  ,'c'  , 'b0','slope' ,'var.dais')
  parameters0 = c(  0.1574, 0.6677 ,  2    , 0.35       , 8.7  , 0.012, 0.35, 0.04       , 1.2 , 1471 , 95  , 775 , 0.0006 , 0.0004656)
  bound.lower = c( 0.0    , 0      ,  0.5  , 0          , 7.05 , 0.003,0.026, 0.025      , 0.6 , 735.5, 47.5, 740 , 0.00045, 0        )
  bound.upper = c( 1.0    , 2      ,  4.25 , 1          , 13.65, 0.015, 1.5 , 0.085      , 1.8 ,2206.5,142.5, 820 , 0.00075, 2        )
}
    # >>> with anto? and same prior ranges as Ruckert et al (2016, in prep) <<<
#parnames    = c('anto.a','anto.b','gamma','alpha.dais','mu'  ,'nu'  ,'P0' ,'kappa.dais','f0' ,'h0'  ,'c'  , 'b0','slope' ,'var.dais')
#parameters0 = c( 0.1574 , 0.6677 , 2     , 0.35       , 8.7  , 0.012, 0.35, 0.04       , 1.2 , 1471 , 95  , 775 , 0.0006 , 0.0004656)
#bound.upper = c( 1.0    , 2      ,  4.25 , 1          , 13.05, 0.018,0.525, 0.06       , 1.8 ,2206.5,142.5, 825 , 0.00075, 2        )
#bound.lower = c( 0.0    , 0      ,  0.5  , 0          , 4.35 , 0.006,0.175, 0.02       , 0.6 , 735.5, 47.5, 725 , 0.00045, 0        )
    # >>> without anto? and same prior ranges as Ruckert et al (2016, in prep) <<<
#parnames    = c('gamma','alpha.dais','mu'  ,'nu'  ,'P0' ,'kappa.dais','f0' ,'h0'  ,'c'  , 'b0','slope' ,'var.dais')
#parameters0 = c(  2    , 0.35       , 8.7  , 0.012, 0.35, 0.04       , 1.2 , 1471 , 95  , 775 , 0.0006 , 0.0004656)
#bound.upper = c(  4.25 , 1          , 13.05, 0.018,0.525, 0.06       , 1.8 ,2206.5,142.5, 825 , 0.00075, 2        )
#bound.lower = c(  0.5  , 0          , 4.35 , 0.006,0.175, 0.02       , 0.6 , 735.5, 47.5, 725 , 0.00045, 0        )

## For prior distributions:
alpha.var = 2     # alpha parameter for inverse gamma for var (E[x]=beta/(alpha+1))
beta.var = 1      # beta parameter for inverse gamma for var (uncertainty parameter)
                  # note that the upper bound on var.dais is not actually imposed; just for plotting
shape.lambda = 8.1              # gives 5% quantile at lambda=0.005 and
rate.lambda = 100*shape.lambda  # gives mean at 0.01 m/yr, DeConto and Pollard (2016)
rate.Tcrit = 1.37               # gives 5% quantile at Tcrit = -10 deg C
shape.Tcrit = 15*rate.Tcrit     # gives mean at -15 deg C (negative requires multiplication of Tcrit by -1)

## Set the f0, h0 and c bounds, +/- a factor of "fac" from the defaults (Shaffer (2014))
## (bounds are from Ruckert et al (2016))
#fac = 0.5   # factor by which f0, h0 and c ranges are set
#bound.upper[match('f0',parnames)]=(1+fac)*parameters0[match('f0',parnames)]
#bound.upper[match('h0',parnames)]=(1+fac)*parameters0[match('h0',parnames)]
#bound.upper[match('c',parnames)]=(1+fac)*parameters0[match('c',parnames)]
#bound.lower[match('f0',parnames)]=(1-fac)*parameters0[match('f0',parnames)]
#bound.lower[match('h0',parnames)]=(1-fac)*parameters0[match('h0',parnames)]
#bound.lower[match('c',parnames)]=(1-fac)*parameters0[match('c',parnames)]

## Source the DAIS model
## (note that the _fastdyn version is used either way; if you selected
##  l.aisfastdy=FALSE, it is neglected in the Fortran codes called)
source('../fortran/R/daisanto_fastdynF.R')
##==============================================================================

## Get a standard hindcast and future projections (a "control" model)
if(l.aisfastdy) {
  Tcrit0 <- parameters0[match('Tcrit',parnames)]
  lambda0 <- parameters0[match('lambda',parnames)]
} else {
  Tcrit0 <- NULL
  lambda0 <- NULL
}
if(TRUE){
  AIS_melt = daisanto_fastdynF(
                anto.a = parameters0[match("anto.a",parnames)],
                anto.b = parameters0[match("anto.b",parnames)],
                gamma  = parameters0[match("gamma",parnames)],
                alpha  = parameters0[match("alpha.dais",parnames)],
                mu     = parameters0[match("mu",parnames)],
                nu     = parameters0[match("nu",parnames)],
                P0     = parameters0[match("P0",parnames)],
                kappa  = parameters0[match("kappa.dais",parnames)],
                f0     = parameters0[match("f0",parnames)],
                h0     = parameters0[match("h0",parnames)],
                c      = parameters0[match("c",parnames)],
                b0     = parameters0[match("b0",parnames)],
                slope  = parameters0[match("slope",parnames)],
                Tcrit  = Tcrit0,
                lambda = lambda0,
                l.aisfastdy = l.aisfastdy,
                Tg     = Tg.recon,
                slope.Ta2Tg = slope.Ta2Tg,
                intercept.Ta2Tg = intercept.Ta2Tg,
                SL     = SL,
                dSL    = dSL)
}
if(FALSE){
  AIS_melt = dais_fastdynF(
                gamma = parameters0[match("gamma",parnames)],
                alpha = parameters0[match("alpha.dais",parnames)],
                mu    = parameters0[match("mu",parnames)],
                nu    = parameters0[match("nu",parnames)],
                P0    = parameters0[match("P0",parnames)],
                kappa = parameters0[match("kappa.dais",parnames)],
                f0    = parameters0[match("f0",parnames)],
                h0    = parameters0[match("h0",parnames)],
                c     = parameters0[match("c",parnames)],
                b0    = parameters0[match("b0",parnames)],
                slope = parameters0[match("slope",parnames)],
                Tcrit = Tcrit0,
                lambda= lambda0,
                l.aisfastdy = l.aisfastdy,
                Ta    = Ta,
                Toc   = Toc,
                SL    = SL,
                dSL   = dSL,
                includes_dSLais=1   # yes, these data include AIS contribution to dSL
                )
}

##==============================================================================
## Set up (pre-)calibration windows around the data
## The point: We run the model at many parameter values and see which ones
##            send the simulation through a window (width determined by the
##            observational errors) around the observational data.
## These windows are presented in Shaffer (2014) and Shepherd et al. (2012)
## 1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
## We want the cumulative sea-level equivalent in meters for the year 2002
## Note the conversion: 360Gt = 1mm SLE
## A fifth window is added to match IPCC AR5 Ch13 (page 1151) AIS SLR trend:
## 0.27 +/- 0.11 mm/year (convert to m/year here)

## Precal windows 5-?:
## Last "precalibration window" is 1993-2010 mean trend, from the IPCC AR5 Ch13
## (Page 1151), for AIS SLR contribution: 0.27 +- 0.11 mm/year
## Note that model output is in meters SLE and these trends are mm, so a
## conversion is necessary.

trends.ais = c(0.27 , 0.08 , 0.40 )/1000   # m/year (after the /1000)
trends.err = c(0.11 , 0.185, 0.205)/1000   # m/year (after the /1000)
trends.2up = trends.ais+2*trends.err
trends.2dn = trends.ais-2*trends.err
ind.trends = mat.or.vec( length(trends.ais), 2)
ind.trends[1,] = c(which(date==-7) , which(date==10)) # 1993-2010
ind.trends[2,] = c(which(date==-8) , which(date== 1)) # 1992-2001
ind.trends[3,] = c(which(date== 2) , which(date==11)) # 2002-2011

## Precal window 4:
## Adding observational constraint
estimate.SLE.rate = abs(-71/360)/1000
time.years = 2002-1992      # using the midpoint of the 19-year interval
mid.cum.SLE_2002 = estimate.SLE.rate*time.years
i1992 = which(date==-8)

estimate.SLE.rate.error = abs(-53/360)/1000     #1-sigma error
estimate.SLE.error = sqrt(time.years)*estimate.SLE.rate.error #1-sigma error
        # (*sqrt(10) because 10 years of potentially accumulated error:
        #  total error^2 = year 1 error^2 + year 2 error^2 + ... year 10 error^2
        #                = 10*year X error^2)
SE2_2002 = estimate.SLE.error*2 #2-sigma error

positive_2SE = mid.cum.SLE_2002 + SE2_2002 # Add the 2 standard error to the mean value
negative_2SE = mid.cum.SLE_2002 - SE2_2002 # Subtract the 2 standard error to the mean value

## Precal windows 1-3:
## from Shaffer (2014). modified by Kelsey
upper.wind = c(7.4, -6.9, -1.25, positive_2SE) # Windows 2-3 from Kelsey, Window 1 from DeConto and Pollard 2016
lower.wind = c(3.6, -15.8, -4.0, negative_2SE)
#upper.wind = c(6.0, -6.9, -1.25, positive_2SE) # Windows 1-3 from Kelsey
#lower.wind = c(1.8, -15.8, -4.0, negative_2SE)
#upper.wind = c(5.5, -8 , -2, positive_2SE) # Windows 1-3 fFrom Shaffer 2014, p 1809
#lower.wind = c(2.5, -17, -4, negative_2SE)

windows = matrix(c(lower.wind, upper.wind), nrow = length(upper.wind), ncol=2)
obs.targets = (windows[,2]+windows[,1])*.5    # middle of window = obs to compare model to
obs.err = (windows[,2]-windows[,1])*.5       # half-width of window = uncertainty
obs.err=0.5*obs.err                           # assume all windows are 2*stdErr
                                              # (last two actually are)
## Create a vector with each observation year
## 120kyr, 20Kyr, 6kyr (before present), 2002, and 1993 (first year of the IPCC
## trend)
obs.years = c(120000, 220000, 234000, 240002)

##==============================================================================




##==============================================================================
## Preliminary Latin Hypercube to find decent starting parameter values

require(lhs)

# Use the same statistical model - calculate log-likelihood, pick starting values
# to maximize this.
source('../calibration/DAIS_assimLikelihood.R')

t0=proc.time()
# Draw LHS sample
n.lhs = 1000
parameters.lhs <- randomLHS(n.lhs, length(parnames))

# Transform unif(0,1) to the parameter bounds
for (i in 1:length(parnames)) {
  parameters.lhs[,i] <- qunif(parameters.lhs[,i], bound.lower[i], bound.upper[i])
}
colnames(parameters.lhs)=parnames

lpost.lhs = rep(NA,n.lhs)
t1=proc.time()
pb <- txtProgressBar(min=0,max=n.lhs,initial=0,style=3)
for (j in 1:n.lhs) {
  lpost.lhs[j] = log.post(  parameters.in = as.numeric(parameters.lhs[j,]),
                            parnames.in = parnames,
                            bound.lower.in = bound.lower,
                            bound.upper.in = bound.upper,
                            obs.in = obs.targets,
                            obs.err.in = obs.err,
                            obs.step.in = obs.years,
                            trends.ais.in = trends.ais,
                            trends.err.in = trends.err,
                            ind.trends.in = ind.trends,
                            ind.norm.in = ind.relative,
                            alpha.var = alpha.var,
                            beta.var = beta.var,
                            shape.lambda = shape.lambda,
                            rate.lambda = rate.lambda,
                            shape.Tcrit = shape.Tcrit,
                            rate.Tcrit = rate.Tcrit,
                            l.aisfastdy = l.aisfastdy,
                            slope.Ta2Tg.in = slope.Ta2Tg,
                            intercept.Ta2Tg.in = intercept.Ta2Tg,
                            Tg.in = Tg.recon,
                            SL.in = SL,
                            dSL.in = dSL
                            )
  setTxtProgressBar(pb, j)
}
close(pb)
t2=proc.time()

parameters0.lhs = parameters.lhs[which(lpost.lhs==max(lpost.lhs)),]

## Get an updated model to compare against control
if(l.aisfastdy) {
  Tcrit0 <- parameters0.lhs[match('Tcrit',parnames)]
  lambda0 <- parameters0.lhs[match('lambda',parnames)]
} else {
  Tcrit0 <- NULL
  lambda0 <- NULL
}
AIS_melt0.lhs = daisanto_fastdynF(anto.a=parameters0.lhs[match("anto.a",parnames)],
                                  anto.b=parameters0.lhs[match("anto.b",parnames)],
                                  gamma = parameters0.lhs[match("gamma",parnames)],
                                  alpha = parameters0.lhs[match("alpha.dais",parnames)],
                                  mu    = parameters0.lhs[match("mu",parnames)],
                                  nu    = parameters0.lhs[match("nu",parnames)],
                                  P0    = parameters0.lhs[match("P0",parnames)],
                                  kappa = parameters0.lhs[match("kappa.dais",parnames)],
                                  f0    = parameters0.lhs[match("f0",parnames)],
                                  h0    = parameters0.lhs[match("h0",parnames)],
                                  c     = parameters0.lhs[match("c",parnames)],
                                  b0    = parameters0.lhs[match("b0",parnames)],
                                  slope = parameters0.lhs[match("slope",parnames)],
                                  Tcrit = Tcrit0,
                                  lambda= lambda0,
                                  l.aisfastdy = l.aisfastdy,
#                                  Ta    = Ta,
#                                  Toc   = Toc,
                                  Tg    = Tg.recon,
                                  slope.Ta2Tg = slope.Ta2Tg,
                                  intercept.Ta2Tg = intercept.Ta2Tg,
                                  SL    = SL,
                                  dSL   = dSL,
                                  includes_dSLais=1
                                  )
##==============================================================================


##==============================================================================
## MCMC simulation
##==============================================================================

## Source the statistical models
source('../calibration/DAIS_assimLikelihood.R')

#install.packages('adaptMCMC')
library(adaptMCMC)
accept.mcmc = 0.234										# Optimal as # parameters->infinity
																			#	(Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = niter_mcmc000											# number of iterations for MCMC
nnode.mcmc = nnode_mcmc000                        # number of logical CPUs for parallel MCMC, if you will use it
gamma.mcmc = 0.5											# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)				# remove first ?? of chains for burn-in
stopadapt.mcmc = round(niter.mcmc*1.0)# stop adapting after ?? iterations? (niter*1 => don't stop)
step.mcmc = (bound.upper-bound.lower)*.05
##==============================================================================
## Actually run the calibration

if(FALSE){
t.beg=proc.time()											# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out1 = MCMC(log.post, niter.mcmc, parameters0.lhs, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames    , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  obs.in=obs.targets      , obs.err.in=obs.err        , obs.step.in=obs.years     ,
                  trends.ais.in=trends.ais, trends.err.in=trends.err  , ind.trends.in=ind.trends  ,
                  ind.norm=ind.relative   , alpha.var=alpha.var       , beta.var=beta.var         ,
                  slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , Tg.in=Tg.recon,
                  shape.lambda=shape.lambda  , rate.lambda=rate.lambda, shape.Tcrit=shape.Tcrit   ,
                  rate.Tcrit=rate.Tcrit      , l.aisfastdy=l.aisfastdy,
                  #Ta.in=Ta  , Toc.in=Toc ,
                  SL.in=SL  , dSL.in=dSL)
t.end=proc.time()											# save timing
chain1 = amcmc.out1$samples
}

## Extend and run more MCMC samples?
if(FALSE){
t.beg=proc.time()
amcmc.extend1 = MCMC.add.samples(amcmc.out1, niter.mcmc,
									parnames.in=parnames    , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  obs.in=obs.targets      , obs.err.in=obs.err        , obs.step.in=obs.years     ,
                  trends.ais.in=trends.ais, trends.err.in=trends.err  , ind.trends.in=ind.trends  ,
                  ind.norm=ind.relative   , alpha.var=alpha.var       , beta.var=beta.var         ,
                  slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , Tg.in=Tg.recon,
                  shape.lambda=shape.lambda  , rate.lambda=rate.lambda, shape.Tcrit=shape.Tcrit   ,
                  rate.Tcrit=rate.Tcrit      , l.aisfastdy=l.aisfastdy,
                  #Ta.in=Ta  , Toc.in=Toc ,
                  SL.in=SL  , dSL.in=dSL)
t.end=proc.time()											# save timing
chain1 = amcmc.extend1$samples
}

## If you want to run 2 (or more) chains in parallel (save time, more sampling)
## Quick note about using MCMC.parallel and MCMC.add.samples:
##    You can add samples to parallel chains by simply replacing "amcmc.out1"
##    with "amcmc.par1" in the MCMC.add.samples command above BUT these R
##    functions do not have the capability to add samples in parallel. It does
##    the task, but does it in serial after the initial chains were in parallel.
if(TRUE){
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.par1 = MCMC.parallel(log.post, niter.mcmc, parameters0.lhs, n.chain=nnode.mcmc, n.cpu=nnode.mcmc,
                  dyn.libs='../fortran/dais_fastdyn.so', scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames    , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  obs.in=obs.targets      , obs.err.in=obs.err        , obs.step.in=obs.years     ,
                  trends.ais.in=trends.ais, trends.err.in=trends.err  , ind.trends.in=ind.trends  ,
                  ind.norm=ind.relative   , alpha.var=alpha.var       , beta.var=beta.var         ,
                  slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , Tg.in=Tg.recon,
                  shape.lambda=shape.lambda  , rate.lambda=rate.lambda, shape.Tcrit=shape.Tcrit   ,
                  rate.Tcrit=rate.Tcrit      , l.aisfastdy=l.aisfastdy,
                  #Ta.in=Ta  , Toc.in=Toc ,
                  SL.in=SL  , dSL.in=dSL)
t.end=proc.time()											# save timing
#chain1=amcmc.par1[[1]]$samples
#chain2=amcmc.par1[[2]]$samples
#chain3=amcmc.par1[[3]]$samples
#chain4=amcmc.par1[[4]]$samples
}

## Save workspace image - you do not want to re-simulate all those!
save.image(file=filename.saveprogress)

##==============================================================================
if(FALSE){
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
}
##==============================================================================

## Determine when (in increments of 10,000 iterations, using Gelman and Rubin
## diagnostic) the two parallel chains, which are assumed to be
## chain1=amcmc.par1[[1]]$samples and chain2=amcmc.par1[[2]]$samples
## as defined above in the MCMC.parallel(...) command.

## Initialize the testing of the Gelman and Rubin diagnostics
niter.test = seq(from=50000, to=niter.mcmc, by=10000)
gr.test = rep(NA,length(niter.test))

## Calculate the statistic at a few spots throughout the chain. Once it is
## close to 1 (people often use GR<1.1 or 1.05), the between-chain variability
## is indistinguishable from the within-chain variability, and they are
## converged. It is only after the chains are converged that you should use the
## parameter values as posterior draws, for analysis.
if(nnode.mcmc == 1) {
  # don't do GR stats, just cut off first half of chains
  print('only one chain; will lop off first half for burn-in instead of doing GR diagnostics')
} else if(nnode.mcmc > 1) {
  # this case is FAR more fun
  # accumulate the names of the soon-to-be mcmc objects
  string.mcmc.list <- 'mcmc1'
  for (m in 2:nnode.mcmc) {
    string.mcmc.list <- paste(string.mcmc.list, ', mcmc', m, sep='')
  }
  for (i in 1:length(niter.test)) {
    for (m in 1:nnode.mcmc) {
      # convert each of the chains into mcmc object
      eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc.par1[[m]]$samples[1:niter.test[i],])', sep='')))
    }
    eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))
    gr.test[i] <- as.numeric(gelman.diag(mcmc_chain_list)[2])
  }
} else {print('error - nnode.mcmc < 1 makes no sense')}

## Save workspace image
save.image(file=filename.saveprogress)

## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO
## <<<<<<< automate the rest of this as in the MESS codes
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< TODO

## Plot GR statistics as a function of iterations, decide where to cut off
## chains and use the tails of both for analysis
#plot(niter.test,gr.test)

#===============================================================================
# Chop off burn-in
#===============================================================================
#

# Note: here, we are only using the Gelman and Rubin diagnostic. But this is
# only after looking at the quantile stability as iterations increase, as well
# as the Heidelberger and Welch diagnostics, which suggest the chains are okay.
# 'ifirst' is the first spot where the GR stat gets to and stays below gr.max
# for all of the models.
# save a separate ifirst for each experiment
ifirst <- NULL
if(nnode.mcmc==1) {
  ifirst <- round(0.5*niter.mcmc)
} else {
  gr.max <- 1.1
  lgr <- rep(NA, length(niter.test))
  for (i in 1:length(niter.test)) {lgr[i] <- gr.test[i] < gr.max}
  for (i in seq(from=length(niter.test), to=1, by=-1)) {
    if( all(lgr[i:length(lgr)]) ) {ifirst <- niter.test[i]}
  }
}

if(nnode.mcmc > 1) {
  chains_burned <- vector('list', nnode.mcmc)
  for (m in 1:nnode.mcmc) {
   chains_burned[[m]] <- amcmc.par1[[m]]$samples[(ifirst+1):niter.mcmc,]
  }
} else {
  chains_burned <- amcmc.out1$samples[(ifirst+1):niter.mcmc,]
}

# Combine all of the chains from 'ifirst' to 'niter_mcmc' into a potpourri of
# [alleged] samples from the posterior. Only saving the transition covariance
# matrix for one of the chains (if in parallel).
if(nnode.mcmc==1) {
  parameters.posterior <- chains_burned
  covjump.posterior <- amcmc.out1$cov.jump
} else {
  parameters.posterior <- chains_burned[[1]]
  covjump.posterior <- amcmc.par1[[1]]$cov.jump
  for (m in 2:nnode.mcmc) {
    parameters.posterior <- rbind(parameters.posterior, chains_burned[[m]])
  }
}
n.parameters = ncol(parameters.posterior)

# save results in case you need to revisit later
save.image(file=filename.saveprogress)

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
