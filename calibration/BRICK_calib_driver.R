##==============================================================================
## Script to define, run and calibrate the coupled BRICK model
##
## To modify the model parameters, their prior ranges, or suggested initial values,
## edit the script 'BRICK_parameterSetup.R' before sourcing it.
##
## To modify which model components you use, edit the settings below under
## 'Which model(s) will you use?'
##
##  Questions? -- Tony Wong <twong@psu.edu
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

rm(list=ls())                        # Clear all previous variables

## Switch to your BRICK calibration directory
setwd('/home/scrim/axw322/codes/BRICK/calibration')
#setwd('/Users/tony/codes/BRICK/calibration')

## Set up MCMC stuff here so that it can be automated for HPC
nnode_mcmc000 <- 8
niter_mcmc000 <- 2e6

## Show plots? (probably want FALSE on HPC, non-interactive)
l.doplots <- FALSE

## Set up a filename for saving RData images along the way
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.saveprogress <- paste('BRICK_calib_MCMC_',today,'.RData',sep='')

## Set the seed (for reproducibility)
set.seed(1234)
#set.seed(as.double(Sys.time())) # should yield same distributions... (a good test!)

## Do you want to use RCP8.5 to make projections? (l.project=TRUE)
## Or do you want to use historical data to make hindcasts? (l.project=FALSE)
## Note -- l.project=FALSE => Kriegler (2005) data, same as Urban and Keller (2010)
l.project = FALSE
#begyear = 1765  # SNEASY start date
begyear = 1850  # DOECLIM start date
endyear = 2009
tstep   = 1
mod.time= seq(from=begyear, to=endyear, by=tstep)
begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)

## Source the models
source('../fortran/R/sneasyF.R')        # the SNEASY model (includes DOECLIM, and carbon cycle)
source('../fortran/R/doeclimF.R')       # the DOECLIM model
source('../fortran/R/GSIC_magiccF.R')   # the GSIC model
source('../fortran/R/brick_te_F.R')     # TE (thermosteric expansion) model
source('../fortran/R/brick_tee_F.R')    # TEE (explicit thermosteric expansion) model
source('../fortran/R/simpleF.R')        # GIS (Greenland Ice Sheet) model
source('../fortran/R/daisanto_fastdynF.R') # DAIS (Antarctic Ice Sheet) model
source('../R/brick_lws.R')              # LWS (land water storage)

## Source some useful functions for manipulating data
source('../R/forcing_total.R')          # function to add up the total forcing
source('../R/compute_indices.R')        # function to determine the model and
                                        # data indices for comparisons
##==============================================================================

## Read the model calibration data sets

  ## TODO
  ## TODO -- revise SNEASY_readData.R to match other components
  ## TODO

#source('../calibration/SNEASY_readData.R')    # read SNEASY calibration data
source('../calibration/DOECLIM_readData.R')   # read DOECLIM calibration data
source('../calibration/GSIC_readData.R')      # read GSIC calibration data
source('../calibration/TE_readData.R')        # read TE data
source('../calibration/SIMPLE_readData.R')    # GIS data, and trends in mass balance
#source('../calibration/DAIS_readData.R')     # DAIS forcing data (if at all uncoupled)


  ## TODO
  ## TODO -- add SNEASY calibration data to midx, oidx, obs, obs.err,
  ## TODO -- ind.norm.data, and i0
  ## TODO

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
## Which model(s) will you use?
## If you want to plug in your own model, insert the relevant "luse.XXX" line
## below, as well as into the "luse.brick = ..." command.
luse.sneasy   = FALSE    # Simple Nonlinear EArth SYstem model (DOECLIM+CCM)
luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # thermosteric expansion contribution to SLR
luse.tee      = FALSE   # explicit thermosteric expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = FALSE    # Antarctic ice sheet model
luse.lws      = FALSE    # land water storage
luse.brick = cbind(luse.sneasy, luse.doeclim, luse.gsic, luse.te, luse.tee, 
		   luse.simple, luse.dais, luse.lws)

## If you are using DAIS, include the fast dynamics emulator?
l.aisfastdy = FALSE
if(!luse.dais) {l.aisfastdy = FALSE} # force FALSE if not using DAIS

if(luse.te & luse.tee) {
  luse.tee = FALSE 
  print('Only use 1 thermosteric expansion model; switching off explicit model.')
}
##==============================================================================
## Define parameters and their prior ranges
## -> Note: 'parnames' is defined here, which establishes how the parameters
##    are passed around into DEoptim, MCMC, likelihood functions, and the models
source('../calibration/BRICK_parameterSetup.R')

##==============================================================================
## Get the forcing data

if(luse.sneasy) {

  ## SNEASY
  setup.sneasy() # call this to initialize SNEASY

## TODO
## TODO -- (YG, TW) likely will need modified, to make sure forcing is as we want it
## TODO

  forcing <- vector('list',3)
  names(forcing) <- c('co2','aero','other')

  # load emissions data time series.
  rcp8.5.emis = read.csv("../data/sneasy_tmp/RCP85_EMISSIONS.csv")
  ibeg = which(rcp8.5.emis[,1]==begyear)
  iend = which(rcp8.5.emis[,1]==endyear)
  if(length(iend)*length(ibeg)==0) {print('ERROR - emissions data does not span mod.time')}
  emis = rcp8.5.emis[ibeg:iend,2] + rcp8.5.emis[ibeg:iend,3] # fossil + land use
  emisdata = data.frame(cbind(mod.time, emis))
  colnames(emisdata) = c("year","co2")
  forcing$co2 = emisdata$co2

  forcing.dat = read.table("../data/sneasy_tmp/forcing_rcp85.txt", header=TRUE)
  ibeg = which(forcing.dat[,1]==begyear)
  iend = which(forcing.dat[,1]==endyear)
  if(length(iend)*length(ibeg)==0) {print('ERROR - other radiative forcing data does not span mod.time')}
  aero.tmp  <- forcing.dat$aerosol.direct + forcing.dat$aerosol.indirect # aerosol radiative forcing
  other.tmp <- forcing.dat$ghg.nonco2 + forcing.dat$solar + forcing.dat$volcanic + forcing.dat$other
  forcing$aero  <- aero.tmp[ibeg:iend]
  forcing$other <- other.tmp[ibeg:iend]

} else {

  if(l.project) {
    forcing = read.csv( '../data/forcing_rcp85.csv', header=TRUE )
  } else {
    forcing = read.csv( '../data/forcing_hindcast.csv', header=TRUE )
  }

}

##==============================================================================
## Define the coupled model
## -> need it defined before DEoptim, so we can calculate the objective function
source('../R/BRICK_coupledModel.R')

##==============================================================================
## Use differential optimization (DEoptim) algorithm to find suitable initial
## parameters for the MCMC chains

##TODO
##TODO -- TW, YG -- need to incorporate SNEASY CO2 (and MOC data?) into optimization
##TODO -- Also note that some of the SNEASY parameter bounds are infinite, which
##TODO -- causes problems in optimization.
##TODO

library(DEoptim)
source('../calibration/BRICK_DEoptim.R')
p0.deoptim=p0                          # initialize optimized initial parameters
niter.deoptim=200                      # number of iterations for DE optimization
NP.deoptim=11*length(index.model)      # population size for DEoptim (do at least 10*[N parameters])
F.deoptim=0.8                          # as suggested by Storn et al (2006)
CR.deoptim=0.9                        # as suggested by Storn et al (2006)
outDEoptim <- DEoptim(minimize_residuals_brick, bound.lower[index.model], bound.upper[index.model],
        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,CR=CR.deoptim,trace=FALSE),
        parnames.in=parnames[index.model], forcing.in=forcing        , l.project=l.project      ,
        #slope.Ta2Tg.in=slope.Ta2Tg       , intercept.Ta2Tg.in=intercept.Ta2Tg,
        ind.norm.data=ind.norm.data      , ind.norm.sl=ind.norm      , mod.time=mod.time        ,
        tstep=tstep                      , oidx = oidx.all           , midx = midx.all          ,
        obs=obs.all                      , obs.err = obs.err.all     , trends.te = trends.te    ,
        luse.brick = luse.brick           , i0 = i0                   , l.aisfastdy = l.aisfastdy )
p0.deoptim[index.model] = outDEoptim$optim$bestmem

## Run the model and examine output at these parameter values
brick.out = brick_model(parameters.in=p0.deoptim,
                        parnames.in=parnames,
                        forcing.in=forcing,
                        l.project=l.project,
                        #slope.Ta2Tg.in=slope.Ta2Tg,
                        #intercept.Ta2Tg.in=intercept.Ta2Tg,
                        mod.time=mod.time,
                        tstep=tstep,
                        ind.norm.data = ind.norm.data,
                        ind.norm.sl = ind.norm,
                        luse.brick = luse.brick,
                        i0 = i0,
                        l.aisfastdy = l.aisfastdy)

##TODO -- TW -- modify plotting for only the components that are used (incl SNEASY)
if(l.doplots) {
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
}

##==============================================================================
## Establish a gamma prior for drawing 1/tau.
## gamma distribution is the conjugate prior for the uncertain parameter beta
## in exponentially distributed random variable V(t):
##    dV/dt = -beta*V => V(t) = exp(-beta*t)
## For the TE model, this is not quite accurate, but similar enough (and
## previous testing with uniform priors on 1/tau confirm) that the gamma is an
## excellent approximation of the distribution of 1/tau (or tau)

shape.invtau = NULL    # initialize
scale.invtau = NULL    # initialize

if(luse.te) {

  invtau.hat = 1/200     # preliminary guess at what the expected value of 1/tau should be
                         # (The timescale and extent of thermal expansion of the oceans due to climate change, Marcelja, 2009)
  q05 = 1/1290           # 82-1290 y are bounds given by Mengel et al (2015) for tau
  q95 = 1/82             # use them as the 5-95% bounds for our 1/tau

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

  niter.deoptim=1000        # number of iterations for DE optimization
  NP.deoptim=50             # population size for DEoptim (do at least 10*[N parameters])
  F.deoptim=0.8             # as suggested by Storn et al (2006)
  CR.deoptim=0.9            # as suggested by Storn et al (2006)
  outDEoptim <- DEoptim(rmse.quantiles, c(0,0), c(100,100),
                        DEoptim.control(NP=NP.deoptim,itermax=niter.deoptim,F=F.deoptim,
            CR=CR.deoptim,trace=FALSE),  q05.in=q05, q95.in=q95)
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

}
##==============================================================================

##==============================================================================
## Now set up the coupled model calibration
##  -- will need to whip up a coupled model likelihood function file
##   -- will calibrate with log.post(...) (in the above LL file) calling 'brick_model'
## Notes:
##  -- With DOECLIM+GSIC+TE+SIMPLE, all R models, requires ~1900s for 1e6 iterations.
##==============================================================================
## Set up and run the MCMC calibration

##TODO
##TODO -- TW, YG -- need to incorporate SNEASY CO2 (and MOC data?) into calibration
##TODO

## Source the statistical models
source('../calibration/BRICK_assimLikelihood.R')

## MCMC calibration
require('adaptMCMC')
library(adaptMCMC)                # use robust adaptive Metropolis
accept.mcmc = 0.234               # Optimal as # parameters->infinity
                                  # (Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = niter_mcmc000        # number of iterations for MCMC
nnode.mcmc = nnode_mcmc000        # number of nodes for parallel MCMC
gamma.mcmc = 0.5                  # rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)    # remove first ?? of chains for burn-in (not used)
stopadapt.mcmc = round(niter.mcmc*1.0)# stop adapting after ?? iterations? (niter*1 => don't stop)

##==============================================================================
## Actually run the calibration
if(FALSE){
t.beg=proc.time()                      # save timing (running millions of iterations so best to have SOME idea...)
amcmc.out1 = MCMC(log.post, niter.mcmc, p0.deoptim, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
                  gamma=gamma.mcmc               , list=TRUE                  , n.start=round(0.01*niter.mcmc),
                  parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
                  #slope.Ta2Tg.in=slope.Ta2Tg    , intercept.Ta2Tg.in=intercept.Ta2Tg,
                  ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
                  oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
                  obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
                  bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
                  luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       )
t.end=proc.time()                      # save timing
chain1 = amcmc.out1$samples
}

## Extend and run more MCMC samples?
if(FALSE){
t.beg=proc.time()
amcmc.extend1 = MCMC.add.samples(amcmc.out1, niter.mcmc,
                    parnames.in=parnames           , forcing.in=forcing         , l.project=l.project            ,
                    #slope.Ta2Tg.in=slope.Ta2Tg    , intercept.Ta2Tg.in=intercept.Ta2Tg,
                    ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time              ,
                    oidx = oidx.all                , midx = midx.all            , obs=obs.all                    ,
                    obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower     ,
                    bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau          ,
                    luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy        )
t.end=proc.time()
chain1 = amcmc.extend1$samples
}

## If you want to run 2 (or more) chains in parallel (save time, more sampling)
if(TRUE){
t.beg=proc.time()                    # save timing (running millions of iterations so best to have SOME idea...)
amcmc.par1 = MCMC.parallel(log.post, niter.mcmc, p0.deoptim, n.chain=nnode.mcmc, n.cpu=nnode.mcmc,
                  dyn.libs=c('../fortran/doeclim.so','../fortran/brick_te.so','../fortran/brick_tee.so','../fortran/gsic_magicc.so','../fortran/simple.so'),
                  scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
                  gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
                  parnames.in=parnames           , forcing.in=forcing         , l.project=l.project           ,
                  #slope.Ta2Tg.in=slope.Ta2Tg    , intercept.Ta2Tg.in=intercept.Ta2Tg,
                  ind.norm.data=ind.norm.data    , ind.norm.sl=ind.norm       , mod.time=mod.time             ,
                  oidx = oidx.all                , midx = midx.all            , obs=obs.all                   ,
                  obs.err = obs.err.all          , trends.te = trends.te      , bound.lower.in=bound.lower    ,
                  bound.upper.in=bound.upper     , shape.in=shape.invtau      , scale.in=scale.invtau         ,
                  luse.brick=luse.brick          , i0=i0                      , l.aisfastdy=l.aisfastdy       )
t.end=proc.time()                      # save timing
}

if(luse.sneasy) {cleanup.sneasy()}  # deallocates memory after SNEASY is done

## Save workspace image - you do not want to re-simulate all those!
save.image(file=filename.saveprogress)

##==============================================================================

## Diagnostic plots
if(l.doplots) {
## Check #1: History plots
par(mfrow=c(5,5))
for (pp in 1:length(parnames)) {
  plot(chain1[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="")
}
}

##==============================================================================

## Determine when (in increments of 50,000 iterations, using Gelman and Rubin
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

## Plot GR statistics as a function of iterations, decide where to cut off
## chains and use the tails of both for analysis
if(l.doplots) {
plot(niter.test,gr.test)
}

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

## Histograms
if(FALSE) {
par(mfrow=c(5,5))
for (pp in 1:length(parnames)) {
  hist(parameters.posterior[,pp], xlab=parnames[pp], main='')
}
}

## Write the calibrated parameters file (netCDF version)

## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape BRICK expects, just transpose it
lmax=0
for (i in 1:length(parnames)){lmax=max(lmax,nchar(parnames[i]))}

## Name the output file
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.parameters = paste('../output_calibration/BRICK-model_calibratedParameters_',today,'.nc',sep="")

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
