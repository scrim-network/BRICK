# =======================================================================================
# Create projections of local sea level and storm surge.
# Use storm surge projections a la Grinsted et al 2013, with T=temperature
# anomaly relative to 1980-2000 average.
# And temperature projections from DOECLIM (within BRICK v0.1) under RCP2.6,
# 4.5 and 8.5.
# =======================================================================================
#
#   Requires (input variables):
# - temperature  projections of global mean surface temperature anomaly
# - time.proj    years associated with temperature projections
# - tidegauge    tide gauge data in the form [date | time | sl] where date is
#                YYYYMMDD, time is hh:mm, sl is meters, and these are the column
#                names (for referencing easily)
# - l.nonstat    which gev parameters to be non-stationary? these will covary
#                with temperature anomaly relative to 1980-2000 mean.
#                Expects list with names $location, $shape, $scale. Note that
#                Coles et al. 2001 (book: An Introduction to Statistical Modeling
#                of Extreme Values) suggest to leave the scale parameter stationary.
# - niter.gev    number of MCMC iterations to use to estimate (non-)stationary
#                GEV parameters
# - burnin       between 0 and 1, what fraction of the MCMC chain should be
#                discarded for burn-in? (default is 0.5)
#
#   Simulates (output variables):
# - gev.location    time series of location parameter (mu) (Nens x Ntime)
# - gev.shape       time series of shape parameter (xi) (Nens x Ntime)
# - gev.scale       time series of location parameter (sigma) (Nens x Ntime)
#                   does each of these for the RCP scenarios that are the "names"
#                   of the input temperature list
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

brick_surge_grinsted <- function(
    temperature = NULL,
    time.proj = NULL,
    tidegauge = NULL,
    l.nonstat,
    niter.gev,
    burnin=0.5)
{

#===============================================================================
# Storm surge
#===============================================================================

library(extRemes)
library(fExtremes)
library(ismev)
library(lubridate)
library(zoo)

# Get tide gauge data and prepare to analyze.
tidegauge$date2 <- as.Date(as.character(tidegauge$date), format="%Y%m%d")
years         <- floor(as.numeric(as.yearmon(tidegauge$date2)))
years.unique  <- unique(years)
n.years       <- length(years.unique)
lsl.mean      <- rep(0,length(n.years))
lsl.max       <- rep(0,length(n.years))
tidegauge$sl.norm <- rep(NA,length(years))

# 1. get unique years of data
# 2. subtract annual block means
# 3. get annual block maxima <-- these are what we fit GEV distribution to
for (tt in 1:n.years) {
    ind.thisyear <- which(years==years.unique[tt])
    lsl.mean[tt] <- mean(tidegauge$sl[ind.thisyear])
    tidegauge$sl.norm[ind.thisyear] <- tidegauge$sl[ind.thisyear] - lsl.mean[tt]
    lsl.max[tt] <- max( tidegauge$sl.norm[ind.thisyear] )
}

# Read temperature data for non-stationary GEV co-variate
# HADCRUT4 annual global mean surface temperature
# Note: all ensemble member files have same time series of ucnertainties, so just
# grabbing the first one.
data.temperature <- read.table("../data/HadCRUT.4.4.0.0.annual_ns_avg.txt")
obs.temp <- data.temperature[,2]
obs.temp.time <- data.temperature[,1]
obs.temp.err <- read.table("../data/HadCRUT.4.4.0.0.annual_ns_avg_realisations/HadCRUT.4.4.0.0.annual_ns_avg.1.txt")[,3]
ibeg <- which(obs.temp.time==years.unique[1])
iend <- which(obs.temp.time==max(years.unique))
obs.temp <- obs.temp[ibeg:iend]
obs.temp.time <- obs.temp.time[ibeg:iend]
obs.temp.err <- obs.temp.err[ibeg:iend]

# Normalize temperature relative to 1980-2000 mean, and create data 'obs'
ibeg <- which(obs.temp.time==1980)
iend <- which(obs.temp.time==2000)
obs.temp <- obs.temp - mean(obs.temp[ibeg:iend])
obs <- data.frame(cbind(years.unique,obs.temp))
colnames(obs) <- c('year','temp')

# Normalize the temperature projections too
ibeg <- which(time.proj==1980)
iend <- which(time.proj==2000)
n.time <- nrow(temperature[[1]])
for (rcp in names(temperature)) {
    temperature[[rcp]] <- temperature[[rcp]] - t(replicate(nrow(temperature[[rcp]]),apply(temperature[[rcp]][ibeg:iend,],2,mean)))
}

# Initialize the projections of the three GEV parameters; use "temperature" to
# assign the proper names and storage of double prec variables
gev_proj <- vector('list',3); names(gev_proj) <- c('location','shape','scale')
for (i in 1:3) {gev_proj[[i]] <- temperature}

# cases: which GEV parameters are to covary with temperature?
if(all(!unlist(l.nonstat))) {

    # stationary GEV case -- just return a time series for each that is constant

    # fit a preliminary maximum likelihood estimate
    gev.mle <- fevd(coredata(lsl.max), type='GEV') # extRemes
    gev.mle2 <- gev.fit(coredata(lsl.max), show = FALSE)   # ismev package

    # initial parameters estimates
    init <- vector('list',3)
    init[[1]] <- gev.mle$results$par[1]
    init[[2]] <- gev.mle$results$par[2]
    init[[3]] <- gev.mle$results$par[3]
    names(init) <- names(gev.mle$results$par)
    std.err <- gev.mle2$se
    names(std.err) <- names(gev.mle$results$par)

    # proposal sd estimates from a preliminary longer chain, with the default sd
    params.prop <- vector('list',1) # list order: mu0, m1, log.scale, xi0, xi1
    params.prop[[1]] <- as.numeric(1*std.err)
    params.prop[[1]][2] <- as.numeric(2*std.err[2])
    names(params.prop) <- 'sd'

    params.pri <- vector('list',1) # list order: mu0, m1, log.scale, xi0, xi1
    params.pri[[1]] <- as.numeric(4*std.err)
    names(params.pri) <- 'v'

    # Fit non-stationary GEV
    gev.bayes <- fevd(coredata(lsl.max), type="GEV", method="Bayesian",
                      initial=init, iter=niter.gev, proposalParams=params.prop, priorParams=params.pri)

    gev.est <- gev.bayes$results[(round(burnin*niter.gev)+1):niter.gev, 1:length(std.err)]
    i.scale <- match('log.scale',names(gev.bayes$results[1,]))
    gev.est[,i.scale] <- exp(gev.est[,i.scale])   # account for log(scale) from MCMC
    colnames(gev.est) <- names(gev.bayes$results[1,1:ncol(gev.est)])
    colnames(gev.est)[i.scale] <- 'scale'

    # check acceptance rates - these should be close to 0.44 (for Metropolis/Gibbs
    # hybrid method, essentially single-variable)
    paccept <- apply(gev.bayes$chain.info[2:nrow(gev.bayes$chain.info),],2,sum)/nrow(gev.bayes$chain.info)
    paccept <- paccept[1:length(std.err)]
    print(paste('Acceptance rates should be around 0.44:',paccept))
    print('If they are not, go into BRICK_project_surge.R and tune the proposal distributions using params.prop')

    # draw parameter sets
    ind.ensemble <- sample( seq(1,nrow(gev.est)), size=n.ensemble, replace=FALSE)
    gev.sample <- gev.est[ind.ensemble,]
    for (i in 1:n.ensemble) {
        for (rcp in names(temperature)) {
            gev_proj$location[[rcp]][,i] <- rep(gev.sample[i,'location'], n.time)
            gev_proj$shape[[rcp]][,i] <- rep(gev.sample[i,'shape'], n.time)
            gev_proj$scale[[rcp]][,i] <- rep(gev.sample[i,'scale'], n.time)
        }
    }

} else if(l.nonstat$location & !l.nonstat$shape & !l.nonstat$scale) {

    # non-stationary location only

    # fit a preliminary maximum likelihood estimate
    gev.mle <- fevd(coredata(lsl.max), type='GEV',
                    data=obs, location.fun=~temp) # extRemes
    gev.mle2 <- gev.fit(coredata(lsl.max), show = FALSE, ydat=obs, mul=2)   # ismev package

    # initial parameters estimates
    init <- vector('list',4)
    init[[1]] <- gev.mle$results$par[1]
    init[[2]] <- 0
    init[[3]] <- gev.mle$results$par[2]
    init[[4]] <- gev.mle$results$par[3]
    names(init) <- names(gev.mle$results$par)
    std.err <- gev.mle2$se
    names(std.err) <- names(gev.mle$results$par)

    # proposal sd estimates from a preliminary longer chain, with the default sd
    params.prop <- vector('list',1) # list order: mu0, m1, log.scale, xi0, xi1
    params.prop[[1]] <- as.numeric(1*std.err)
    params.prop[[1]][3] <- as.numeric(2*std.err[3])
    names(params.prop) <- 'sd'

    params.pri <- vector('list',1) # list order: mu0, m1, log.scale, xi0, xi1
    params.pri[[1]] <- as.numeric(4*std.err)
    names(params.pri) <- 'v'

    # Fit non-stationary GEV
    gev.bayes <- fevd(coredata(lsl.max), type="GEV", method="Bayesian",
                      data=obs, location.fun=~temp,
                      initial=init, iter=niter.gev, proposalParams=params.prop, priorParams=params.pri)

if(FALSE) {
    # fit multiple chains to get GR diagnostics - otherwise, one long chain
    gev.bayes1 <- fevd(coredata(lsl.max), type="GEV", method="Bayesian",
                      data=obs, location.fun=~temp,
                      initial=init, iter=niter.gev, proposalParams=params.prop, priorParams=params.pri)
    gev.bayes2 <- fevd(coredata(lsl.max), type="GEV", method="Bayesian",
                      data=obs, location.fun=~temp,
                      initial=init, iter=niter.gev, proposalParams=params.prop, priorParams=params.pri)
    gev.bayes3 <- fevd(coredata(lsl.max), type="GEV", method="Bayesian",
                      data=obs, location.fun=~temp,
                      initial=init, iter=niter.gev, proposalParams=params.prop, priorParams=params.pri)
    gev.bayes4 <- fevd(coredata(lsl.max), type="GEV", method="Bayesian",
                      data=obs, location.fun=~temp,
                      initial=init, iter=niter.gev, proposalParams=params.prop, priorParams=params.pri)
    chain1 <- gev.bayes1$results[,1:4]
    chain2 <- gev.bayes2$results[,1:4]
    chain3 <- gev.bayes3$results[,1:4]
    chain4 <- gev.bayes4$results[,1:4]
    niter.test = seq(from=1000, to=nrow(chain1), by=1000)
    gr.stat = rep(NA,length(niter.test))
    for (i in 1:length(niter.test)){
        mcmc1 = as.mcmc(chain1[1:niter.test[i],])
        mcmc2 = as.mcmc(chain2[1:niter.test[i],])
        mcmc3 = as.mcmc(chain3[1:niter.test[i],])
        mcmc4 = as.mcmc(chain4[1:niter.test[i],])
        mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2, mcmc3, mcmc4))
        gr.stat[i] = gelman.diag(mcmc_chain_list)[2]
    }
    i1 <- niter.test[which(gr.stat<1.05)[1]]+1
    iend <- nrow(chain1)
    gev.posterior <- rbind(chain1[i1:iend,], chain2[i1:iend,], chain3[i1:iend,], chain4[i1:iend,])

    gev.est <- gev.posterior
    i.scale <- match('log.scale',names(gev.bayes$results[1,]))
    gev.est[,i.scale] <- exp(gev.est[,i.scale])   # account for log(scale) from MCMC
    colnames(gev.est) <- names(gev.bayes$results[1,1:ncol(gev.est)])
    colnames(gev.est)[i.scale] <- 'scale'

    # check acceptance rates - these should be close to 0.44 (for Metropolis/Gibbs
    # hybrid method, essentially single-variable)
    paccept <- apply(gev.bayes$chain.info[2:nrow(gev.bayes$chain.info),],2,sum)/nrow(gev.bayes$chain.info)
    paccept <- paccept[1:length(std.err)]
    print(paste('Acceptance rates should be around 0.44:',paccept))
    print('If they are not, go into BRICK_project_surge.R and tune the proposal distributions using params.prop')

    # draw parameter sets
    ind.ensemble <- sample( seq(1,nrow(gev.est)), size=n.ensemble, replace=FALSE)
    gev.sample <- gev.est[ind.ensemble,]
    for (i in 1:n.ensemble) {
        for (rcp in names(temperature)) {
            gev_proj$location[[rcp]][,i] <- gev.sample[i,'mu0'] + gev.sample[i,'mu1']*temperature[[rcp]][,i]
            gev_proj$shape[[rcp]][,i] <- rep(gev.sample[i,'shape'], n.time)
            gev_proj$scale[[rcp]][,i] <- rep(gev.sample[i,'scale'], n.time)
        }
    }

} else {
    gev.est <- gev.bayes$results[(round(burnin*niter.gev)+1):niter.gev, 1:length(std.err)]
    i.scale <- match('log.scale',names(gev.bayes$results[1,]))
    gev.est[,i.scale] <- exp(gev.est[,i.scale])   # account for log(scale) from MCMC
    colnames(gev.est) <- names(gev.bayes$results[1,1:ncol(gev.est)])
    colnames(gev.est)[i.scale] <- 'scale'

    # check acceptance rates - these should be close to 0.44 (for Metropolis/Gibbs
    # hybrid method, essentially single-variable)
    paccept <- apply(gev.bayes$chain.info[2:nrow(gev.bayes$chain.info),],2,sum)/nrow(gev.bayes$chain.info)
    paccept <- paccept[1:length(std.err)]
    print(paste('Acceptance rates should be around 0.44:',paccept))
    print('If they are not, go into BRICK_project_surge.R and tune the proposal distributions using params.prop')

    # draw parameter sets
    ind.ensemble <- sample( seq(1,nrow(gev.est)), size=n.ensemble, replace=FALSE)
    gev.sample <- gev.est[ind.ensemble,]
    for (i in 1:n.ensemble) {
        for (rcp in names(temperature)) {
            gev_proj$location[[rcp]][,i] <- gev.sample[i,'mu0'] + gev.sample[i,'mu1']*temperature[[rcp]][,i]
            gev_proj$shape[[rcp]][,i] <- rep(gev.sample[i,'shape'], n.time)
            gev_proj$scale[[rcp]][,i] <- rep(gev.sample[i,'scale'], n.time)
        }
    }
}

} else if(l.nonstat$location & l.nonstat$shape & !l.nonstat$scale) {

    # non-stationary location and shape only

    # fit a preliminary maximum likelihood estimate
    gev.mle <- fevd(coredata(lsl.max), type='GEV',
                    data=obs, location.fun=~temp, shape.fun=~temp) # extRemes
    gev.mle2 <- gev.fit(coredata(lsl.max), show = FALSE, ydat=obs, mul=2, shl=2)   # ismev package

    # initial parameters estimates
    init <- vector('list',5)
    init[[1]] <- gev.mle$results$par[1]
    init[[2]] <- 0
    init[[3]] <- gev.mle$results$par[2]
    init[[4]] <- gev.mle$results$par[3]
    init[[5]] <- 0
    names(init) <- names(gev.mle$results$par)
    std.err <- gev.mle2$se
    names(std.err) <- names(gev.mle$results$par)

    # proposal sd estimates from a preliminary longer chain, with the default sd
    params.prop <- vector('list',1) # list order: mu0, m1, log.scale, xi0, xi1
    params.prop[[1]] <- as.numeric(1*std.err)
    params.prop[[1]][3] <- as.numeric(2*std.err[3])
    names(params.prop) <- 'sd'

    params.pri <- vector('list',1) # list order: mu0, m1, log.scale, xi0, xi1
    params.pri[[1]] <- as.numeric(4*std.err)
    names(params.pri) <- 'v'

    # Fit non-stationary GEV
    gev.bayes <- fevd(coredata(lsl.max), type="GEV", method="Bayesian",
                      data=obs, location.fun=~temp, shape.fun=~temp,
                      initial=init, iter=niter.gev, proposalParams=params.prop, priorParams=params.pri)

    gev.est <- gev.bayes$results[(round(burnin*niter.gev)+1):niter.gev, 1:length(std.err)]
    i.scale <- match('log.scale',names(gev.bayes$results[1,]))
    gev.est[,i.scale] <- exp(gev.est[,i.scale])   # account for log(scale) from MCMC
    colnames(gev.est) <- names(gev.bayes$results[1,1:ncol(gev.est)])
    colnames(gev.est)[i.scale] <- 'scale'

    # check acceptance rates - these should be close to 0.44 (for Metropolis/Gibbs
    # hybrid method, essentially single-variable)
    paccept <- apply(gev.bayes$chain.info[2:nrow(gev.bayes$chain.info),],2,sum)/nrow(gev.bayes$chain.info)
    paccept <- paccept[1:length(std.err)]
    print(paste('Acceptance rates should be around 0.44:',paccept))
    print('If they are not, go into BRICK_project_surge.R and tune the proposal distributions using params.prop')

    # draw parameter sets
    ind.ensemble <- sample( seq(1,nrow(gev.est)), size=n.ensemble, replace=FALSE)
    gev.sample <- gev.est[ind.ensemble,]
    for (i in 1:n.ensemble) {
        for (rcp in names(temperature)) {
            gev_proj$location[[rcp]][,i] <- gev.sample[i,'mu0'] + gev.sample[i,'mu1']*temperature[[rcp]][,i]
            gev_proj$shape[[rcp]][,i] <- gev.sample[i,'xi0'] + gev.sample[i,'xi1']*temperature[[rcp]][,i]
            gev_proj$scale[[rcp]][,i] <- rep(gev.sample[i,'scale'], n.time)
        }
    }

} else {
    print('Sorry, that combination of non-stationary GEV parameters is not supported.')
    #TODO -- add support for other combinations
}


# get projections of the 99% quantile (1:100 level) of 
if(FALSE) {
    tmp2000 <- rep(0,n.ensemble); i2000 <- which(t.proj==2000)
    tmp2050 <- rep(0,n.ensemble); i2050 <- which(t.proj==2050)
    tmp2100 <- rep(0,n.ensemble); i2100 <- which(t.proj==2100)
    for (i in 1:n.ensemble) {
        tmp2000[i] <- qevd(p=0.99,loc=gev_proj$location$rcp85[i2000,i],scale=gev_proj$scale$rcp85[i2000,i],shape=gev_proj$shape$rcp85[i2000,i])
        tmp2050[i] <- qevd(p=0.99,loc=gev_proj$location$rcp85[i2050,i],scale=gev_proj$scale$rcp85[i2050,i],shape=gev_proj$shape$rcp85[i2050,i])
        tmp2100[i] <- qevd(p=0.99,loc=gev_proj$location$rcp85[i2100,i],scale=gev_proj$scale$rcp85[i2100,i],shape=gev_proj$shape$rcp85[i2100,i])
    }

}
#===============================================================================

    output <- list(gev_proj, time.proj)
    names(output) <- c('gev.proj','time.proj')
    return(output)
}

#===============================================================================
# end
#===============================================================================
