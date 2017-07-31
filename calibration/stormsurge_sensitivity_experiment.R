##==============================================================================
##  stormsurge_sensitivity_experiment.R
##
##  Estimate storm surge GEV parameters for Galveston, Texas and Pensacola,
##  Florida using 35 year blocks and assess how sensitive the parameters are to
##  data availability/length.
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

setwd('~/codes/BRICK/calibration')

rm(list=ls())

nnode_mcmc000 <- 4
niter_mcmc000 <- 5e4
gamma_mcmc000 <- 0.5

library(adaptMCMC)
library(extRemes)
library(date)
library(Hmisc)
library(zoo)

# Name the saved progress RData workspace image file
today <- Sys.Date(); today=format(today,format="%d%b%Y")
filename.saveprogress <- paste('../output_calibration/stormsurge_sensitivity_mcmc_',today,'.RData', sep='')

##==============================================================================

data.tg <- vector('list', 2); names(data.tg) <- c('galveston','pensacola')
for (dd in 1:length(data.tg)) {data.tg[[dd]] <- vector('list', 5); names(data.tg[[dd]]) <- c('year','month','day','hour','sl')}

data.tmp <- read.table('../data/rqh0775_galveston.csv', header=FALSE, sep=',')
for (col in 1:5) {data.tg$galveston[[col]] <- data.tmp[[col]]}

data.tmp <- read.table('../data/rqh0762a_pensacola.csv', header=FALSE, sep=',')
for (col in 1:5) {data.tg$pensacola[[col]] <- data.tmp[[col]]}

for (dd in 1:length(data.tg)) {
  names(data.tg[[dd]]) <- c('year','month','day','hour','sl')
  data.tg[[dd]]$time.days <- as.numeric(mdy.date(month=data.tg[[dd]]$month, day=data.tg[[dd]]$day,
                                                  year=data.tg[[dd]]$year)) + data.tg[[dd]]$hour/24
}

##==============================================================================
## Take annual block maxima for ALL years (will get rid of the ones with too
## much missing data later)
##==============================================================================

ind.years.to.remove <- vector('list', length(data.tg)); names(ind.years.to.remove) <- names(data.tg)
for (dd in 1:length(data.tg)) {
  years.with.data <- unique(data.tg[[dd]]$year)
  data.tg[[dd]]$lsl.max <- rep(NA, length(years.with.data))
  data.tg[[dd]]$lsl.mean <- rep(NA, length(years.with.data))
  data.tg[[dd]]$lsl.norm <- rep(NA, length(data.tg[[dd]]$sl))
  data.tg[[dd]]$year.max <- years.with.data
  for (year in years.with.data) {
    ind.this.year <- which(data.tg[[dd]]$year == year)
    data.tg[[dd]]$lsl.mean[match(year,years.with.data)] <- mean(data.tg[[dd]]$sl[ind.this.year])
    data.tg[[dd]]$lsl.norm[ind.this.year] <- data.tg[[dd]]$sl[ind.this.year] - data.tg[[dd]]$lsl.mean[match(year,years.with.data)]
    data.tg[[dd]]$lsl.max[match(year,years.with.data)] <- max(data.tg[[dd]]$lsl.norm[ind.this.year])
    # Find any years with less than 90% of the data and mark them for removal
    ndata.fullyear <- yearDays(paste(year,'-1-1', sep='')) * 24
    perc.data.this.year <- length(ind.this.year) / ndata.fullyear
    if(perc.data.this.year < 0.9) {ind.years.to.remove[[dd]] <- c(ind.years.to.remove[[dd]], match(year, years.with.data))}
  }
  data.tg[[dd]]$lsl.max <- data.tg[[dd]]$lsl.max[-ind.years.to.remove[[dd]]]
  data.tg[[dd]]$lsl.mean <- data.tg[[dd]]$lsl.mean[-ind.years.to.remove[[dd]]]
  data.tg[[dd]]$year.max <- data.tg[[dd]]$year.max[-ind.years.to.remove[[dd]]]
}

##==============================================================================
## Decompose into sets of 35-year blocks.
## Moving window, by 16 years each time. Gives roughly half-overlap between
## adjacent blocks, and good coverage of the 101 years (goes to the 99th year)
## from Galveston and complete coverage of the 83 years from Pensacola (once we
## remove the years with too much missing data).
##==============================================================================

block.size   <- 35 # how many years in each block?
block.offset <- 8 # how many years are the blocks shifted relative to neighboring blocks?

ind.block.right.endpt <- vector('list', length(data.tg)); names(ind.block.right.endpt) <- names(data.tg)
ind.block.left.endpt <- vector('list', length(data.tg)); names(ind.block.left.endpt) <- names(data.tg)
nblocks <- rep(NA, length(data.tg)); names(nblocks) <- names(data.tg)
block.names <- vector('list', length(data.tg)); names(block.names) <- names(data.tg)

for (dd in 1:length(data.tg)) {
  ind.block.right.endpt[[dd]] <- seq(from=length(data.tg[[dd]]$year.max), to=block.size, by=-block.offset)
  nblocks[[dd]] <- length(ind.block.right.endpt[[dd]])
  ind.block.left.endpt[[dd]] <- ind.block.right.endpt[[dd]]-block.size+1

  # create some other data objects for calibration for each experiment block
  for (bb in 1:nblocks[[dd]]) {block.names[[dd]] <- c(block.names[[dd]], paste('block',bb,sep=''))}
  data.tg[[dd]]$blocks <- vector('list', nblocks[[dd]]); names(data.tg[[dd]]$blocks) <- block.names[[dd]]
  for (bb in 1:nblocks[[dd]]) {
    data.tg[[dd]]$blocks[[bb]] <- vector('list', 2); names(data.tg[[dd]]$blocks[[bb]]) <- c('year.max','lsl.max')
    data.tg[[dd]]$blocks[[bb]]$year.max <- data.tg[[dd]]$year.max[ind.block.left.endpt[[dd]][bb]:ind.block.right.endpt[[dd]][bb]]
    data.tg[[dd]]$blocks[[bb]]$lsl.max <- data.tg[[dd]]$lsl.max[ind.block.left.endpt[[dd]][bb]:ind.block.right.endpt[[dd]][bb]]
  }
}

##==============================================================================
## Preliminary MCMC steps
##==============================================================================

## Set up prior distributions
parnames <- c('mu','sigma','xi')
priors <- vector('list', 3); names(priors) <- parnames

# uniform priors
for (par in 1:3) {priors[[par]] <- vector('list', 3); names(priors[[par]]) <- c('type','lower','upper'); priors[[par]]$type <- 'uniform'}
priors$mu$lower <- 0
priors$mu$upper <- 5000
priors$sigma$lower <- 0
priors$sigma$upper <- 1000
priors$xi$lower <- -3
priors$xi$upper <- 3

## Load likelihood functions from auxiliary file
source('stormsurge_sensitivity_experiment_aux-likelihood.R')

nnode_mcmc <- nnode_mcmc000
niter_mcmc <- niter_mcmc000
gamma_mcmc <- gamma_mcmc000
startadapt_mcmc <- max(500,round(0.05*niter_mcmc))
stopadapt_mcmc <- round(niter_mcmc*1.0)
accept_mcmc_few <- 0.44         # optimal for only one parameter
accept_mcmc_many <- 0.234       # optimal for many parameters
accept_mcmc <- accept_mcmc_many + (accept_mcmc_few - accept_mcmc_many)/length(parnames)

amcmc_out <- vector('list', length(data.tg)); names(amcmc_out) <- names(data.tg)
for (dd in 1:length(data.tg)) {amcmc_out[[dd]] <- vector('list', nblocks[[dd]]); names(amcmc_out[[dd]]) <- block.names[[dd]]}

##==============================================================================
## Actually run the MCMC - fit GEV parameters (stationary model)
##==============================================================================

for (dd in 1:length(data.tg)) {
  print(paste('Starting calibrations for site ',names(data.tg)[dd],' (',nnode_mcmc,' cores x ',niter_mcmc,' iterations)...', sep=''))
  for (bb in 1:nblocks[[dd]]) {
    print(paste('-- block ',bb,' / ',nblocks[[dd]], sep=''))
    tbeg <- proc.time()
    # initial parameter estimates from maximum likelihood
    gev.mle <- fevd(coredata(data.tg[[dd]]$blocks[[bb]]$lsl.max))
    initial_parameters <- as.numeric(gev.mle$results$par)
    step_mcmc <- sqrt(diag(1/gev.mle$results$hessian)) # hessian is like curvative, so radius of curvature ~ 1/hessian, a step estimate
    if(nnode_mcmc==1) {
      amcmc_out[[dd]][[bb]] <- MCMC(log_post_gev, niter_mcmc, initial_parameters,
                                   adapt=TRUE, acc.rate=accept_mcmc, scale=step_mcmc,
                                   gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                                   parnames=parnames, data_calib=data.tg[[dd]]$blocks[[bb]]$lsl.max,
                                   priors=priors)
    } else if (nnode_mcmc > 1) {
      amcmc_out[[dd]][[bb]] <- MCMC.parallel(log_post_gev, niter_mcmc, initial_parameters,
                                   n.chain=nnode_mcmc, n.cpu=nnode_mcmc, packages='extRemes',
                                   scale=step_mcmc, adapt=TRUE, acc.rate=accept_mcmc,
                                   gamma=gamma_mcmc, list=TRUE, n.start=startadapt_mcmc,
                                   parnames=parnames, data_calib=data.tg[[dd]]$blocks[[bb]]$lsl.max,
                                   priors=priors)

    } else {print('nnode_mcmc must be positive integer')}
    tend <- proc.time()
    print(paste('... done. Took ',round(as.numeric(tend-tbeg)[3]/60,2),' minutes', sep=''))
    print(paste('... Saving MCMC workspace results as .RData file (',filename.saveprogress,') to read and use later...',sep=''))
    save.image(file=filename.saveprogress)
    print('...done.')
  }
}

##==============================================================================
## Convergence diagnostics and burn-in removal
##==============================================================================

## Gelman and Rubin diagnostics - determine and chop off for burn-in
niter.test <- seq(from=round(0.1*niter_mcmc), to=niter_mcmc, by=round(0.05*niter_mcmc))
gr.test <- vector('list', length(data.tg)); names(gr.test) <- names(data.tg)
for (dd in 1:length(data.tg)) {
  gr.test[[dd]] <- mat.or.vec(length(niter.test), nblocks[[dd]])
  colnames(gr.test[[dd]]) <- block.names[[dd]]
}
gr.tmp <- rep(NA, length(niter.test))

for (dd in 1:length(data.tg)) {
  for (bb in 1:nblocks[[dd]]) {
    if(nnode_mcmc == 1) {
      # don't do GR stats, just cut off first half of chains
      print('only one chain; will lop off first half for burn-in instead of doing GR diagnostics')
    } else if(nnode_mcmc > 1) {
      # this case is FAR more fun
      # accumulate the names of the soon-to-be mcmc objects
      string.mcmc.list <- 'mcmc1'
      for (m in 2:nnode_mcmc) {
        string.mcmc.list <- paste(string.mcmc.list, ', mcmc', m, sep='')
      }
      for (i in 1:length(niter.test)) {
        for (m in 1:nnode_mcmc) {
          # convert each of the chains into mcmc object
          eval(parse(text=paste('mcmc',m,' <- as.mcmc(amcmc_out[[dd]][[bb]][[m]]$samples[1:niter.test[i],])', sep='')))
        }
        eval(parse(text=paste('mcmc_chain_list = mcmc.list(list(', string.mcmc.list , '))', sep='')))

        gr.test[[dd]][i,bb] <- as.numeric(gelman.diag(mcmc_chain_list)[2])
      }
    } else {print('error - nnode_mcmc < 1 makes no sense')}
  }
}

# 'ifirst' is the first spot where the GR stat gets to and stays below gr.max
# for all of the models. save a separate ifirst for each site
ifirst <- rep(NA, length(data.tg)); names(ifirst) <- names(data.tg)
gr.max <- 1.1
for (dd in 1:length(data.tg)) {
  if(nnode_mcmc==1) {
    ifirst[[dd]] <- round(0.5*niter_mcmc)
  } else {
    lgr <- rep(NA, length(niter.test))
    for (i in 1:length(niter.test)) {lgr[i] <- all(gr.test[[dd]][i,] < gr.max)}
    for (i in seq(from=length(niter.test), to=1, by=-1)) {
      if( all(lgr[i:length(lgr)]) ) {ifirst[[dd]] <- niter.test[i]}
    }
  }
}

chains_burned <- vector('list', length(data.tg)); names(chains_burned) <- names(data.tg)
for (dd in 1:length(data.tg)) {
  chains_burned[[dd]] <- vector('list', nblocks[[dd]]); names(chains_burned[[dd]]) <- block.names[[dd]]
  for (bb in 1:nblocks[[dd]]) {
    if(nnode_mcmc > 1) {
      chains_burned[[dd]][[bb]] <- vector('list', nnode_mcmc)
      for (m in 1:nnode_mcmc) {
        chains_burned[[dd]][[bb]][[m]] <- amcmc_out[[dd]][[bb]][[m]]$samples[(ifirst[[dd]]+1):niter_mcmc,]
      }
    } else {
      chains_burned[[dd]][[bb]] <- amcmc_out[[dd]][[bb]]$samples[(ifirst[[dd]]+1):niter_mcmc,]
    }
  }
}

##==============================================================================
## Combine the parallel chains and thin to a manageable size?
##==============================================================================

# If no thinning, then this initialization will remain
chains_burned_thinned <- chains_burned

n.sample <- 12586

for (dd in 1:length(data.tg)) {
  for (bb in 1:nblocks[[dd]]) {
    if(nnode_mcmc == 1) {
      ind.sample <- sample(x=1:nrow(chains_burned[[dd]][[bb]]), size=n.sample, replace=FALSE)
      chains_burned_thinned[[dd]][[bb]] <- chains_burned[[dd]][[bb]][ind.sample,]
    } else {
      n.sample.sub <- rep(NA, nnode_mcmc)
      # the case where desired sample size is divisible by the number of chains
      if(round(n.sample/nnode_mcmc) == n.sample/nnode_mcmc) {
        n.sample.sub[1:nnode_mcmc] <- n.sample/nnode_mcmc
      } else {
      # the case where it is not
        n.sample.sub[2:nnode_mcmc] <- round(n.sample/nnode_mcmc)
        n.sample.sub[1] <- n.sample - sum(n.sample.sub[2:nnode_mcmc])
      }
      for (m in 1:nnode_mcmc) {
        ind.sample <- sample(x=1:nrow(chains_burned[[dd]][[bb]][[m]]), size=n.sample.sub[m], replace=FALSE)
        chains_burned_thinned[[dd]][[bb]][[m]] <- chains_burned[[dd]][[bb]][[m]][ind.sample,]
      }
    }
  }
}

# Combine all of the chains from 'ifirst' to 'niter_mcmc' into a potpourri of
# [alleged] samples from the posterior. Only saving the transition covariance
# matrix for one of the chains (if in parallel).
parameters.posterior <- vector('list', length(data.tg)); names(parameters.posterior) <- names(data.tg)
for (dd in 1:length(data.tg)) {
  parameters.posterior[[dd]] <- vector('list', nblocks[[dd]]); names(parameters.posterior[[dd]]) <- block.names[[dd]]
  for (bb in 1:nblocks[[dd]]) {
    if(nnode_mcmc==1) {
      parameters.posterior[[dd]][[bb]] <- chains_burned_thinned[[dd]][[bb]]
    } else {
      parameters.posterior[[dd]][[bb]] <- chains_burned_thinned[[dd]][[bb]][[1]]
      for (m in 2:nnode_mcmc) {
        parameters.posterior[[dd]][[bb]] <- rbind(parameters.posterior[[dd]][[bb]], chains_burned_thinned[[dd]][[bb]][[m]])
      }
    }
  }
}

##==============================================================================
## Calculate 100-year return level for each of the blocks, for each site.
##==============================================================================

returnperiod.of.interest <- 100 # in years

returnlevel <- vector('list', length(data.tg)); names(returnlevel) <- names(data.tg)
for (dd in 1:length(data.tg)) {
  returnlevel[[dd]] <- vector('list', nblocks[[dd]]); names(returnlevel[[dd]]) <- block.names[[dd]]
  for (bb in 1:nblocks[[dd]]) {
    returnlevel[[dd]][[bb]] <- rep(NA, nrow(parameters.posterior[[dd]][[bb]]))
    returnlevel[[dd]][[bb]] <- sapply(1:nrow(parameters.posterior[[dd]][[bb]]), function(i) {
                                      qevd(p=0.01, loc=parameters.posterior[[dd]][[bb]][i,1],
                                           scale=parameters.posterior[[dd]][[bb]][i,2],
                                           shape=parameters.posterior[[dd]][[bb]][i,3],
                                           type='GEV', lower.tail=FALSE)})
  }
}

##==============================================================================
## Get the NOLA GEV results to superimpose
##==============================================================================

filename.vandantzig.nofd <- '../output_model/VanDantzig_fd-none_2065_08May2017.nc'
ncdata <- nc_open(filename.vandantzig.nofd)
  heightening <- ncvar_get(ncdata, 'H')
  parameters.nola <- ncvar_get(ncdata, 'gev_stat')
  colnames(parameters.nola) <- c('mu','sigma','xi')
nc_close(ncdata)

returnlevel.nola <- rep(NA, nrow(parameters.nola))
returnlevel.nola <- sapply(1:nrow(parameters.nola), function(i) {
                                  qevd(p=0.01, loc=parameters.nola[i,1],
                                       scale=parameters.nola[i,2],
                                       shape=parameters.nola[i,3],
                                       type='GEV', lower.tail=FALSE)})

## Save progress to revisit later
print(paste('... Saving MCMC workspace results as .RData file (',filename.saveprogress,') to read and use later...',sep=''))
save.image(file=filename.saveprogress)
print('...done.')

##==============================================================================
## Get kernel density estimates for each of the distributions and plot
##==============================================================================

returnlevel.kde <- vector('list', length(data.tg)); names(returnlevel.kde) <- names(data.tg)
for (dd in 1:length(data.tg)) {
  returnlevel.kde[[dd]] <- vector('list', nblocks[[dd]]); names(returnlevel.kde[[dd]]) <- block.names[[dd]]
  for (bb in 1:nblocks[[dd]]) {
    returnlevel.kde[[dd]][[bb]] <- density(returnlevel[[dd]][[bb]], from=0, to=15000, n=512)
    returnlevel.kde[[dd]][[bb]]$x <- returnlevel.kde[[dd]][[bb]]$x/1000 # convert to m from mm
  }
}
returnlevel.kde.nola <- density(returnlevel.nola, from=0, to=15000, n=512)
returnlevel.kde.nola$x <- returnlevel.kde.nola$x/1000 # convert from mm to m

# note - can quantify the uncertainty in the PP-year return level by checking out
# the variance in the maximum likelihood posterior estimates?

# block years
block.years <- cbind(data.tg[[1]]$year.max[ind.block.left.endpt[[1]]], data.tg[[1]]$year.max[ind.block.right.endpt[[1]]])
names.block.years <- rep(NA, max(nblocks))
for (bb in 1:length(names.block.years)) {names.block.years[bb] <- paste(block.years[bb,1],block.years[bb,2],sep='-')}

##==============================================================================
## Make the actual figures
##==============================================================================

##
## Figure -- two panels (Galveston and Pensacola); each panel has the distributions
##           of 100-year return level as estimated using each block of 35 years
##

plotdir='~/Box\ Sync/Wong-Projects/BRICK_scenarios/figures/'
pdf(paste(plotdir,'stormsurge_sensitivity_experiments.pdf',sep=''),width=7,height=3.5,colormodel='cmyk')
par(mfrow=c(1,2), mai=c(1,.5,.15,.3))
dd=1 # galveston
plot(returnlevel.kde[[dd]]$block1$x, returnlevel.kde[[dd]]$block1$y,
     type='l', lwd=2, col='darkblue', xlim=c(0,10), ylim=c(0,6.1e-4),
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(returnlevel.kde[[dd]]$block2$x, returnlevel.kde[[dd]]$block2$y, lwd=2, col='darkcyan')
lines(returnlevel.kde[[dd]]$block3$x, returnlevel.kde[[dd]]$block3$y, lwd=2, col='cornflowerblue')
lines(returnlevel.kde[[dd]]$block4$x, returnlevel.kde[[dd]]$block4$y, lwd=2, col='aquamarine')
lines(returnlevel.kde[[dd]]$block5$x, returnlevel.kde[[dd]]$block5$y, lwd=2, col='cadetblue1')
lines(returnlevel.kde.nola$x, returnlevel.kde.nola$y, lwd=2, lty=2, col='black')
axis(1,seq(0,15,2),cex.axis=1.2)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=1);
mtext('100-year return level [m]\nGalveston, Texas', side=1, line=3.5, cex=1);
mtext(side=3, text=expression(bold('   a')), line=-1, cex=.9, adj=0);
text(4.2,5.8e-4, 'Years of data:', pos=4)
legend(4, 5.8e-4, c(names.block.years, 'New Orleans', '(1980-2016)'), lty=c(1,1,1,1,1,2, NA), lwd=2, cex=1.0, bty='n',
       col=c('darkblue','darkcyan','cornflowerblue','aquamarine','cadetblue1','black'))
dd=2 # pensacola
plot(returnlevel.kde[[dd]]$block1$x, returnlevel.kde[[dd]]$block1$y,
     type='l', lwd=2, col='darkblue', xlim=c(0,10), ylim=c(0,1.1e-3),
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(returnlevel.kde[[dd]]$block2$x, returnlevel.kde[[dd]]$block2$y, lwd=2, col='darkcyan')
lines(returnlevel.kde[[dd]]$block3$x, returnlevel.kde[[dd]]$block3$y, lwd=2, col='cornflowerblue')
lines(returnlevel.kde[[dd]]$block4$x, returnlevel.kde[[dd]]$block4$y, lwd=2, col='aquamarine')
lines(returnlevel.kde.nola$x, returnlevel.kde.nola$y, lwd=2, lty=2, col='black')
axis(1,seq(0,15,2),cex.axis=1.2)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=1);
mtext('100-year return level [m]\nPensacola, Florida', side=1, line=3.5, cex=1);
mtext(side=3, text=expression(bold('   b')), line=-1, cex=.9, adj=0);
dev.off()

##==============================================================================
## Box-whisker plot for diagnosing the need for non-stationary approach
##==============================================================================

for (dd in 1:length(data.tg)) {
  for (bb in 1:nblocks[[dd]]) {
    ind.thisblock <- which(data.tg[[dd]]$year <= data.tg[[dd]]$year.max[ind.block.right.endpt[[dd]][bb]] &
                           data.tg[[dd]]$year >= data.tg[[dd]]$year.max[ind.block.left.endpt[[dd]][bb]])
    data.tg[[dd]]$blocks[[bb]]$lsl.norm <- data.tg[[dd]]$lsl.norm[ind.thisblock]
  }
}


##
## Figure -- two panels (top/bottom, Galveston and Pensacola); each panel has
##           the box-whisker plot within each block of the detrended/processed
##           hourly tide gauge data from that block; to show potential trends
##

# TODO -- use data.tg[[dd]]$blocks[[bb]]$lsl.norm quantiles to define boxes

# TODO -- to get better resolution here for trends, may want to stagger the blocks by 10 years?
# Maybe shoot for twice as many blocks, but only present a subset (for aesthetics)
# (so stagger 8 years, and we can still present these)

##
##==============================================================================
## End
##==============================================================================
##
