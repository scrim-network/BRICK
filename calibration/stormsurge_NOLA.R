##==============================================================================
##  stormsurge_NOLA.R
##
##  Two routines for estimating storm surge. One based on tide gauge data, and
##  one based on the US Army Corps of Engineers Coastal Louisiana flood defense
##  manual (circa 2008).
##
## Questions? Tony Wong (twong@psu.edu)
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

##==============================================================================
##
##  BRICK_estimateGEV_NOLA
##
##  Estimation of storm surge GEV parameters for New Orleans, based off of
##  Grand Isle, LA tide gauge data
##
##  Summary:
##  1. use maximum likelihood estimator to get initial guess for GEV parameters
##      some may consider this to be poor practice. If you want, you can run
##      several chains from multiple initial conditions. (Did this, looks fine
##      even with only 5000 iterations, using the default stepsizes and dispersed
##      initial conditions. (TW, 3 Jan 2016))
##  2. run a longer MCMC for the GEV parameters
##  3. chop off a burn-in (first half of chain) and leave the rest of the chains
##      to draw parameters from.
##  4. write to netcdf file
##
##  Required input:
##  niter               number of MCMC iterations to do (total, with burnin included)
##  sd.prop             standard deviations for the Gaussian random walk
##  sd.pri              standard deviations for the normal priors
##                      elsements are: [loc, log(scale), shape]
##  N                   number of ensemble members (draws) to write to a different
##                      file for later reading
##
##  Returns:
##  filename.gevstat    netcdf file with MCMC results for GEV parameters (creates this)
##  parameters.posterior   the gev estimates (chain of them niter*burnin long)
##==============================================================================

BRICK_estimateGEV_NOLA <- function(
  niter=500000,
  N=NULL
  ){

# estimate first year storm surge GEV, using tide gauge data from Grand Isle.
print('Reading tide gauge data...')

dat.dir <- '../data/tideGauge_GrandIsle/'
files.tg <- list.files(path=dat.dir,pattern='csv')

# all dates/times are relative to first observation (1980-11-11 16:00 GMT)
data <- read.csv(paste(dat.dir,'CO-OPS__8761724__hr.csv',sep=''))
origin <- as.POSIXlt(data$Date.Time, tz='GMT')[1]

# first file to start the data array
data.new <- read.csv(paste(dat.dir,files.tg[1],sep=''))
time.new <- as.POSIXlt(data.new$Date.Time, tz='GMT')
time.new.rel <- difftime(time.new, origin, unit='days')
data.new[,1] <- time.new.rel
data.all <- data.new

for (ff in 2:length(files.tg)) {
  data.new <- read.csv(paste(dat.dir,files.tg[ff],sep=''))
  time.new <- as.POSIXlt(data.new$Date.Time, tz='GMT')
  time.new.rel <- difftime(time.new, origin, unit='days')
  data.new[,1] <- time.new.rel
  data.all <- rbind(data.all,data.new)
}

# sort the tide gauge data array by the time.relative column
# (not necessarily in order, and possibly some overlap between adjacent year files)
data.all.sort <- data.all[order(data.all$Date.Time),]

# Get water levels in mm
data.all.sort$Water.Level <- data.all.sort$Water.Level*1000

print('...done!')

print('Estimating stationary GEV parameters...')

# bin by year, or otherwise make the timestamp actually useful
# Note: to get actual datestamp back, just do data.all.sort$Date.TIme + origin
days.year <- 365.25 # number of days in a year, on average
nyear <- floor(as.numeric(max(data.all.sort$Date.Time)/days.year))

# initialize these list arrays
data.all.sort$Water.Level.avg <- data.all.sort$Water.Level
data.all.sort$lsl <- data.all.sort$Water.Level
data.all.sort$lsl.avg <- data.all.sort$Water.Level
data.all.sort$lsl.avg.rel <- data.all.sort$Water.Level
lsl.max <- rep(NA,nyear)  # annual block maxima

# One method:
# - Calculate annual means
# - Subtract these off
# - Calculate annual block maxima
# - Fit GEV
# Alternative method:
# - Use NOAA/USACE SLR trend of 9.24 mm/year (https://tidesandcurrents.noaa.gov/est/curves.shtml?stnid=8761724)
# - Subtract this SLR trend from the tidge gauge data.
# - Normalize the result to have mean of zero.
# - Calculate annual block maxima, fit GEV to these, as usual.
# - Yields result which is lower than the annual means method (both are nice ways
# to account for sea-level rise.)
slr.rate <- 9.24/365  # sea level rise per day
slr <- as.numeric(slr.rate*data.all.sort$Date.Time)

data.all.sort$lsl <- data.all.sort$Water.Level - slr
data.all.sort$lsl <- data.all.sort$lsl - mean(data.all.sort$lsl)
for (tt in 1:nyear) {
  iyear <- which(data.all.sort$Date.Time <  days.year*tt &
                 data.all.sort$Date.Time >= days.year*(tt-1))
  data.all.sort$lsl.avg[iyear] <- rep(mean(data.all.sort$Water.Level[iyear]),length(iyear))
  data.all.sort$lsl.avg.rel[iyear] <- data.all.sort$Water.Level[iyear]-data.all.sort$lsl.avg[iyear]
  # Method 1:
  lsl.max[tt] <- max(data.all.sort$lsl.avg.rel[iyear])
  method.name <- '-AnnMean'
  # Method 2:
  #lsl.max[tt] <- max(data.all.sort$lsl[iyear])
  #method.name <- '-NOAAtrend'
}

# First, estimate GEV parameters using maximum likelihood.
gev.mle <- fevd(coredata(lsl.max)) # extremes
init <- vector("list",3)
init[[1]]=gev.mle$results$par[1]
init[[2]]=gev.mle$results$par[2]
init[[3]]=gev.mle$results$par[3]
names(init) <- names(gev.mle$results$par)

# Define prior distribution function for the GEV parameters
fprior <- function( theta, mu, sigma, lb, ub){
    pri.loc <- dnorm(x=theta[1], mean=mu[1], sd=sigma[1], log=TRUE)
    pri.sca <- dnorm(x=exp(theta[2]), mean=mu[2], sd=sigma[2], log=TRUE)
    if(theta[3] > lb[3] & theta[3] < ub[3]) pri.sha <- 0
    else pri.sha <- -Inf
    pri <- pri.loc + pri.sca + pri.sha
    return(pri)
}

pri.mu <- c(init$location, init$scale, init$shape)
pri.sd <- c(54, 30, 1)
#pri.lb <- c(-Inf, -Inf, 0 )
pri.lb <- c(-Inf, -Inf, -5)
pri.ub <- c( Inf,  Inf, 5 )
params.pri <- vector('list',4)
params.pri[[1]] <- pri.mu
params.pri[[2]] <- pri.sd
params.pri[[3]] <- pri.lb
params.pri[[4]] <- pri.ub
names(params.pri) <- c('mu','sigma','lb','ub')

# These were tuned to obtain acceptance rates ~44% (Rosenthal et al, eg)
params.prop <- vector('list',1)
params.prop[[1]] <- c(22, .15, .23)
names(params.prop) <- 'sd'

# first chain
# use these as the initial guesses for the Bayesian GEV parameter estimation
set.seed(111)
gev.bayes1 <- fevd(coredata(lsl.max), type="GEV", method="Bayesian",
                  initial=init, iter=niter, proposalParams=params.prop,
                  priorFun='fprior', priorParams=params.pri)
gev.est1 <- gev.bayes1$results[ ,1:3]
gev.est1[,2] <- exp(gev.est1[,2])   # account for log(scale) from MCMC
colnames(gev.est1) <- c('location','scale','shape')
apply(gev.bayes1$chain.info[2:niter,1:3],2,sum)/niter

# another chain with different seed and initial condition
set.seed(222)
init.new <- init
init.new$location <- rnorm(n=1, mean=init[[1]], sd=params.pri$sigma[1])
init.new$scale    <- rnorm(n=1, mean=init[[2]], sd=params.pri$sigma[2])
init.new$shape    <- runif(n=1, min=pri.lb[3] , max=pri.ub[3])
gev.bayes2 <- fevd(coredata(lsl.max), type="GEV", method="Bayesian",
                  initial=init.new, iter=niter, proposalParams=params.prop,
                  priorFun='fprior', priorParams=params.pri)
gev.est2 <- gev.bayes2$results[ ,1:3]
gev.est2[,2] <- exp(gev.est2[,2])   # account for log(scale) from MCMC
colnames(gev.est2) <- c('location','scale','shape')
apply(gev.bayes2$chain.info[2:niter,1:3],2,sum)/niter

## Convergence checks
heidel.diag(gev.est1, eps=0.1, pvalue=0.05)
heidel.diag(gev.est2, eps=0.1, pvalue=0.05)

niter.test = seq(from=0.2*niter, to=niter, by=50000)
gr.stat = rep(NA,length(niter.test))
for (i in 1:length(niter.test)){
  mcmc1 = as.mcmc(gev.est1[1:niter.test[i],])
  mcmc2 = as.mcmc(gev.est2[1:niter.test[i],])
  mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2))
  gr.stat[i] = gelman.diag(mcmc_chain_list)[2]
}

plot(niter.test,gr.stat, xlab='iteration', ylab='GR statistic')

for (i in 1:length(niter.test)) {
    lconv <- all(gr.stat[i:length(niter.test)] < 1.05)
    if(lconv) {
        n.burnin <- niter.test[i]
        break
    }
}

n.sample <- niter-n.burnin
parameters1 <- gev.est1[(n.burnin+1):niter,]
parameters2 <- gev.est2[(n.burnin+1):niter,]
parameters.posterior = rbind(parameters1, parameters2)
n.parameters = ncol(parameters.posterior)



## prepare and write output file
## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape BRICK expects, just transpose it
lmax=0
for (i in 1:ncol(parameters.posterior)){lmax=max(lmax,nchar(colnames(parameters.posterior)[i]))}

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.gevstat = paste('../output_calibration/BRICK_estimateGEV',method.name,'_',today,'.nc',sep="")

library(ncdf4)
dim.parameters <- ncdim_def('n.parameters', '', 1:ncol(parameters.posterior), unlim=FALSE)
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(parameters.posterior), unlim=TRUE)
parameters.var <- ncvar_def('GEV_parameters', '', list(dim.parameters,dim.ensemble), -999)
parnames.var <- ncvar_def('GEV_names', '', list(dim.name,dim.parameters), prec='char')
outnc <- nc_create(filename.gevstat, list(parameters.var,parnames.var))
ncvar_put(outnc, parameters.var, t(parameters.posterior))
ncvar_put(outnc, parnames.var, colnames(parameters.posterior))
nc_close(outnc)

## Write shorter file, with only the number (N) for the given ensemble
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.gevshort = paste('../output_calibration/BRICK_GEVsample',method.name,'_',today,'.nc',sep="")
gev.params <- parameters.posterior[sample(1:nrow(parameters.posterior), size=N, replace=FALSE) , ]

dim.parameters <- ncdim_def('n.parameters', '', 1:ncol(gev.params), unlim=FALSE)
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(gev.params), unlim=TRUE)
parameters.var <- ncvar_def('GEV_parameters', '', list(dim.parameters,dim.ensemble), -999)
parnames.var <- ncvar_def('GEV_names', '', list(dim.name,dim.parameters), prec='char')
outnc <- nc_create(filename.gevshort, list(parameters.var,parnames.var))
ncvar_put(outnc, parameters.var, t(gev.params))
ncvar_put(outnc, parnames.var, colnames(gev.params))
nc_close(outnc)

print('...done!')

return(list(parameters.posterior,gev.params))

}

##==============================================================================


##==============================================================================
##
##  BRICK_estimateStormSurge_NOLA_usace
##
##  Estimation of sea level + storm surge heights for New Orleans, based off of
##  Table 1.2, p.1-22, of USACE "Hurricane and Storm Damage Risk Reduction
##  System Design Guidelines" (with revisions through June 2012). Accessible
##  here (as of 30 Nov 2016):
##
##  http://www.mvn.usace.army.mil/Portals/56/docs/engineering/HurrGuide/EntireDocument.pdf
##
##  Summary:
##  * surge increase is 1.5 to 2 times the sea level increase
##    (so sea level needs to be relative to the first time slice)
##
##  Required input:
##  time.proj     years for which storm surge GEV parameters are needed for projections
##  sea_level     sea level rise (m)
##  surge.factor  storm surge increase relative to sea level increase (unitless)
##
##  Returns:
##  surge.rise  height (m) above original (stationary, from gev.stat) storm surge
##              the nonstationary storm surge is.
##==============================================================================

BRICK_estimateStormSurge_NOLA_usace <- function(
  time.proj,
  sea_level,
  surge.factor
  ){

# check that the number of rows of sea level matches time.proj; reshape if not.
dims = dim(sea_level)
if(dims[1]!=length(time.proj)) {sea_level <- t(sea_level)}

# NB - currently no check for if number of time slices in time.proj does not match
# the number of time slices in sea_level.

# normalize sea level rise to the first year, because the increase in surge is
# relative to the increase in sea level (above first year).
sea_level <- sea_level - t(replicate(nrow(sea_level),sea_level[1,]))

# (surge rise) = (surge factor) * (sea level rise)
# surge.factor = 0 indicates the surge has not worsened, so the gev.stat (above)
# is appropriate. surge.rise is the height above the original surge by which the
# situation has worsened. The surge.rise levels are relative to the risen sea
# levels. This will be subtracted from the effective dike height, then the
# original stationary GEV can be applied to get return levels.

surge.rise <- t(replicate(nrow(sea_level), surge.factor)) * sea_level

output = list(time.proj, surge.rise, surge.factor)
names(output) = c('year','surge.rise','surge.factor')

  return(output)
}
##==============================================================================


##==============================================================================
##
##  BRICK_estimateStormSurge_NOLA_tg
##
##  Estimation of storm surge parameters for New Orleans, Louisiana.
##  Using tide gauge data from Grand Isle, LA; Pensacola, FL; Galveston, TX.
##  GEV function - stationary (gev.stat) or non-stationary (gev.nonstat).
##
##  Required input:
##  time.proj   years for which storm surge GEV parameters are needed for projections
##
##  Returns:
##  gev.year    years for which storm surge GEV parameters are given
##  gev.stat    stationary GEV parameters, each row a different year, column 1
##              is location, column 2 is shape, column 3 is scale
##              (Note: for gev.stat, all rows should be the same)
##  gev.nonstat non-stationary GEV parameters, each row a different year, column
##              1 is location, column 2 is shape, column 3 is scale
##
##==============================================================================

BRICK_estimateStormSurge_NOLA_tg <- function(
  time.proj
  ){

#install.packages("extRemes")
#install.packages("fExtremes")
#install.packages('ismev')
#install.packages('zoo')
library(extRemes)
library(fExtremes)
library(ismev)
library(zoo)

## Initialize output
gev.stat = mat.or.vec(length(time.proj),3)
gev.nonstat = mat.or.vec(length(time.proj),3)

print('Reading tide gauge data...')

# these data are in millimeters
dat.la=read.csv('/Users/axw322/Box Sync/Wong-Projects/BRICK_scenarios/tideGauge_hourly_GrandIsle.csv')
dat.fl=read.csv('/Users/axw322/Box Sync/Wong-Projects/BRICK_scenarios/tideGauge_hourly_Pensacola.csv')
dat.tx=read.csv('/Users/axw322/Box Sync/Wong-Projects/BRICK_scenarios/tideGauge_hourly_Galveston.csv')

print('...done!')

print('Estimating GEV parameters...')

## Estimate stationary GEV parameters using Grand Isle data
dat=dat.la
lsl=dat[,5]
year=unique(dat[,1])
nyear=length(year)
lsl.res=rep(NA,length(lsl))
lsl.avg=rep(NA,nyear)
lsl.max=rep(NA,nyear)
for (t in 1:nyear) {
  itmp = which(dat[,1]==year[t])
  lsl.avg[t] = mean(lsl[itmp])
  lsl.res[itmp] = lsl[itmp]-lsl.avg[t]
  lsl.max[t] = max(lsl.res[itmp])
}
#gev.hind <- gev.fit(lsl.max,show=FALSE) # ismev
gev.hind <- fevd(coredata(lsl.max)) # extremes

# extrapolate to time.proj (easy for stationary case...)
for (i in 1:3) {gev.stat[,i] = gev.hind$results$par[i]}

colnames(gev.stat) = c('location','scale','shape')
gev.stat = data.frame(gev.stat)

## Estimate non-stationary GEV parameters using Pensacola and/or Galveston data
# non-stationary fits
# -> using the Grand Isle data, can do 2 of the GEV parameters non-st, but not all 3
# -> and cannot do 'tmp3' for just the scale and shape parameters

dat=dat.fl
lsl=dat[,5]
year=unique(dat[,1])
nyear=length(year)
lsl.res=rep(NA,length(lsl))
lsl.avg=rep(NA,nyear)
lsl.max=rep(NA,nyear)
for (t in 1:nyear) {
  itmp = which(dat[,1]==year[t])
  lsl.avg[t] = mean(lsl[itmp])
  lsl.res[itmp] = lsl[itmp]-lsl.avg[t]
  lsl.max[t] = max(lsl.res[itmp])
}

yr=matrix(cbind(1:length(year),1:length(year),1:length(year)),ncol=3)

# need to use the XXinit inputs...
iblock1 = 1:41
iblock2 = (length(year)-41):length(year)
time.hind = c(mean(year[iblock1]) , mean(year[iblock2]))

gev.hind1 <- fevd(coredata(lsl.max[iblock1])) # extremes
gev.hind2 <- fevd(coredata(lsl.max[iblock2])) # extremes

gev.hind = rbind(gev.hind1$results$par,gev.hind2$results$par)

#tmp=gev.fit(lsl.max, ydat = yr, show=FALSE)
#tmp0=gev.fit(lsl.max, ydat = yr, mul=1, show=FALSE)
#tmp1=gev.fit(lsl.max, ydat = yr, mul=1, sigl=2, show=FALSE)
#tmp2=gev.fit(lsl.max, ydat = yr, mul=1, shl=3, show=FALSE)
#tmp3=gev.fit(lsl.max, ydat = yr, sigl=2, shl=3, show=FALSE)
#tmp4=gev.fit(lsl.max, ydat = yr, mul=1, sigl=2, shl=3, show=FALSE)
#gev.ns1 = tmp1$vals
#gev.ns2 = tmp2$vals

if(FALSE){
ss.level=rep(NA,nrow(gev.hind))
for (t in 1:nrow(gev.hind)){ss.level[t]=qgev(xi=gev.hind[t,'shape'],beta=gev.hind[t,'scale'],mu=gev.hind[t,'location'],p=0.99)}
plot(time.hind,ss.level)
}

# i=2 is linear fit in log(sigma)
if(FALSE){
i=1; fit <- lm(gev.hind[,i] ~ time.hind)$coefficients
gev.nonstat[,i] = fit[1] + fit[2]*time.proj
i=2; fit <- lm(log(gev.hind[,i]) ~ time.hind)$coefficients
gev.nonstat[,i] = exp(fit[1] + fit[2]*time.proj)
i=3; fit <- lm(gev.hind[,i] ~ time.hind)$coefficients
gev.nonstat[,i] = fit[1] + fit[2]*time.proj
}
# instead, start in 2015 from the stationary GEV parameters
i=1; fit <- lm(gev.hind[,i] ~ time.hind)$coefficients
gev.nonstat[,i] = gev.stat[1,i] + fit[2]*(time.proj-2015)
i=2; fit <- lm(log(gev.hind[,i]) ~ time.hind)$coefficients
gev.nonstat[,i] = exp(log(gev.stat[1,i]) + fit[2]*(time.proj-2015))
i=3; fit <- lm(gev.hind[,i] ~ time.hind)$coefficients
gev.nonstat[,i] = gev.stat[1,i] + fit[2]*(time.proj-2015)

colnames(gev.nonstat) = c('location','scale','shape')
gev.nonstat = data.frame(gev.nonstat)

gev.year = time.proj

if(FALSE){
ss.level=rep(NA,nrow(gev.nonstat))
for (t in 1:nrow(gev.nonstat)){ss.level[t]=qgev(xi=gev.nonstat[t,'shape'],beta=gev.nonstat[t,'scale'],mu=gev.nonstat[t,'location'],p=0.99)}
plot(gev.year,ss.level)
}

# for test case, only assume location parameter is non-stationary
#gev.nonstat[,'scale'] = gev.stat[1,'scale']
#gev.nonstat[,'shape'] = gev.stat[1,'shape']

print('...done!')

output = list(gev.year,gev.stat,gev.nonstat)
names(output) = c('year','stat','nonstat')

  return(output)
}
##==============================================================================


##==============================================================================
## End
##==============================================================================
