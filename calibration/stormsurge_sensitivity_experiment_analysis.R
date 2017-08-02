##==============================================================================
##  stormsurge_sensitivity_experiment_analysis.R
##
##  Estimate storm surge GEV parameters for Galveston, Texas and Pensacola,
##  Florida using 35 year blocks and assess how sensitive the parameters are to
##  data availability/length.
##
##  This script does the analysis for the storm surge result sensitviity to data
##  period, based on (subsets of the) Galveston and Pensacola tide gauge records.
##  This assumes you have run either 'stormsurge_sensitivity_experiment.R'
##  (interactively) or 'stormsurge_sensitivity_experiment_driver.R' (thrown to
##  HPC using the 'data_sensitivity_run.pbs' batch submit script), either of
##  which will result in the mcmc RData object read in below (and you will of
##  course need to modify the dates and local directory structure).
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

setwd('~/codes/BRICK/calibration')

rm(list=ls())

library(adaptMCMC)
library(extRemes)
library(date)
library(Hmisc)
library(zoo)

## Name the saved progress RData workspace image file that your
## 'stormsurge_sensitivity_experiment[_driver, possible].R' results were saved to.
filename.sensitivityexperiment <- '../output_calibration/stormsurge_sensitivity_mcmc_31Jul2017.RData'

## Directory where you would like to store the figures this script creates
plotdir='~/Box\ Sync/Wong-Projects/BRICK_scenarios/figures/'


##==============================================================================
## Read the MCMC results for Pensacola and Galveston block experiments
##==============================================================================

load(filename.sensitivityexperiment)

##==============================================================================
## Define colors to use for each block
##==============================================================================

block.colors <- colorRampPalette(c("darkslateblue","royalblue","turquoise1"),space="Lab")(max(nblocks))
block.colors.lighter <- paste(block.colors, "70", sep="")

##==============================================================================
## Figure -- two panels (Galveston and Pensacola); each panel has the distributions
##           of 100-year return level as estimated using each block of 35 years
##==============================================================================

pdf(paste(plotdir,'stormsurge_sensitivity_distributions.pdf',sep=''),width=7,height=3.5,colormodel='cmyk')
par(mfrow=c(1,2), mai=c(1,.5,.15,.3))
dd=1 # galveston
plot(returnlevel.kde[[dd]]$block1$x, returnlevel.kde[[dd]]$block1$y,
     type='l', lwd=2, col='darkblue', xlim=c(0,10), ylim=c(0,6.1e-4),
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(returnlevel.kde[[dd]]$block3$x, returnlevel.kde[[dd]]$block3$y, lwd=2, col='darkcyan')
lines(returnlevel.kde[[dd]]$block5$x, returnlevel.kde[[dd]]$block5$y, lwd=2, col='cornflowerblue')
lines(returnlevel.kde[[dd]]$block7$x, returnlevel.kde[[dd]]$block7$y, lwd=2, col='aquamarine')
lines(returnlevel.kde[[dd]]$block9$x, returnlevel.kde[[dd]]$block9$y, lwd=2, col='cadetblue1')
lines(returnlevel.kde.nola$x, returnlevel.kde.nola$y, lwd=2, lty=2, col='black')
axis(1,seq(0,15,2),cex.axis=1.2)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=1);
mtext('100-year return level [m]\nGalveston, Texas', side=1, line=3.5, cex=1);
mtext(side=3, text=expression(bold('   a')), line=-1, cex=.9, adj=0);
text(4.2,5.8e-4, 'Years of data:', pos=4)
legend(4, 5.8e-4, c(names.block.years[c(1,3,5,7,9)], 'New Orleans', '(1980-2016)'), lty=c(1,1,1,1,1,2, NA), lwd=2, cex=1.0, bty='n',
       col=c('darkblue','darkcyan','cornflowerblue','aquamarine','cadetblue1','black'))
dd=2 # pensacola
plot(returnlevel.kde[[dd]]$block1$x, returnlevel.kde[[dd]]$block1$y,
     type='l', lwd=2, col='darkblue', xlim=c(0,10), ylim=c(0,1.1e-3),
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(returnlevel.kde[[dd]]$block3$x, returnlevel.kde[[dd]]$block3$y, lwd=2, col='darkcyan')
lines(returnlevel.kde[[dd]]$block5$x, returnlevel.kde[[dd]]$block5$y, lwd=2, col='cornflowerblue')
lines(returnlevel.kde[[dd]]$block7$x, returnlevel.kde[[dd]]$block7$y, lwd=2, col='aquamarine')
lines(returnlevel.kde.nola$x, returnlevel.kde.nola$y, lwd=2, lty=2, col='black')
axis(1,seq(0,15,2),cex.axis=1.2)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1, cex=1);
mtext('100-year return level [m]\nPensacola, Florida', side=1, line=3.5, cex=1);
mtext(side=3, text=expression(bold('   b')), line=-1, cex=.9, adj=0);
dev.off()


##==============================================================================
## Figure -- two panels (top/bottom, Galveston and Pensacola); each panel has
##           the box-whisker plot within each block of the detrended/processed
##           hourly tide gauge data from that block; to show potential trends
##==============================================================================

## Want the box/whisker plot of the distribution of return levels that generated
## the KDEs from the prvious plot, and to show these marching forward through
## time.

## Calculate the quantiles to plot
quantiles.to.grab <- c(.05, .25, .5, .75, .95)
quantile.names <- rep(NULL, length(quantiles.to.grab))
for (qq in 1:length(quantiles.to.grab)) {
  if(quantiles.to.grab[qq] >= .10) {
    quantile.names[qq] <- paste('q',100*quantiles.to.grab[qq], sep='')
  } else if(quantiles.to.grab[qq] < .10 & quantiles.to.grab[qq] >= 0) {
    quantile.names[qq] <- paste('q0',100*quantiles.to.grab[qq], sep='')
  }
}
returnlevel.quantiles <- vector('list', length(data.tg)); names(returnlevel.quantiles) <- names(data.tg)
for (dd in 1:length(data.tg)) {
  returnlevel.quantiles[[dd]] <- mat.or.vec(nblocks[[dd]], length(quantiles.to.grab))
  colnames(returnlevel.quantiles[[dd]]) <- quantile.names
  for (bb in 1:nblocks[[dd]]) {
    # the /1000 is to convert to m from mm
    returnlevel.quantiles[[dd]][bb,] <- quantile(returnlevel[[dd]][[bb]], quantiles.to.grab)/1000
  }
}

## Useful for plotting - centers of the time blocks used in the experiments
block.years.center <- apply(X=block.years, MARGIN=1, FUN=median)

## And useful for plotting the NOLA ensemble alongside these results
returnlevel.quantiles.nola <- quantile(returnlevel.nola, quantiles.to.grab)/1000; names(returnlevel.quantiles.nola) <- quantile.names
block.years.center.nola <- max(block.years.center) + median(abs(diff(block.years.center)))


## The actual figure

pdf(paste(plotdir,'stormsurge_sensitivity_boxwhisker.pdf',sep=''),width=4,height=6,colormodel='cmyk')
par(mfrow=c(2,1), mai=c(.8,.7,.15,.2))
halfwidth <- 2 # half the width of the boxes, in years
dd=1 # galveston
# put the first median bar down, to get hte plot started
plot(c(block.years.center[1]-halfwidth, block.years.center[1]+halfwidth), rep(returnlevel.quantiles[[dd]][1,'q50'],2),
     type='l', lwd=3, col='black', xlim=c(1920,2000), ylim=c(0,21), xlab='', ylab='', las=1)
# now add the darker 25-75% range polygon before the median bars, ...
for (bb in 1:nblocks[[dd]]) {
    times.beg.end <- c(block.years.center[bb]-halfwidth, block.years.center[bb]+halfwidth)
    polygon(c(times.beg.end,rev(times.beg.end)), c(returnlevel.quantiles[[dd]][bb,c('q25','q25')],rev(returnlevel.quantiles[[dd]][bb,c('q75','q75')])),
            col=block.colors[bb], border=NA)
}
# ... and add the lighter 5-95% range polygon before the median bars too...
for (bb in 1:nblocks[[dd]]) {
    times.beg.end <- c(block.years.center[bb]-halfwidth, block.years.center[bb]+halfwidth)
    polygon(c(times.beg.end,rev(times.beg.end)), c(returnlevel.quantiles[[dd]][bb,c('q05','q05')],rev(returnlevel.quantiles[[dd]][bb,c('q95','q95')])),
            col=block.colors.lighter[bb], border=NA)
}
# ... so the bars are on top
for (bb in 1:nblocks[[dd]]) {lines(c(block.years.center[bb]-halfwidth, block.years.center[bb]+halfwidth),
                                   rep(returnlevel.quantiles[[dd]][bb,'q50'],2), lwd=3, col='black')}
# finally, add the NOLA data on the far right
#lines(c(2000,2000), c(-1000,1000), lty=3, lwd=2, col='black')
#text(2002, 7, 'New Orleans', srt=90, pos=4)
#times.beg.end <- c(block.years.center.nola-halfwidth, block.years.center.nola+halfwidth)
#polygon(c(times.beg.end,rev(times.beg.end)), c(returnlevel.quantiles.nola[c('q05','q05')],rev(returnlevel.quantiles.nola[c('q95','q95')])),
#        col='gray65', border=NA);
#polygon(c(times.beg.end,rev(times.beg.end)), c(returnlevel.quantiles.nola[c('q25','q25')],rev(returnlevel.quantiles.nola[c('q75','q75')])),
#        col='gray25', border=NA);
#lines(c(block.years.center.nola-halfwidth, block.years.center.nola+halfwidth),
#      rep(returnlevel.quantiles.nola['q50'],2), lwd=3, col='black')
text(1960, 20, 'Galveston, Texas', pos=4)
mtext('Year', side=1, line=2.4, cex=1);
mtext('100-year return level [m]', side=2, line=2.2, cex=1);
mtext(side=3, text=expression(bold('   a')), line=-1, cex=1, adj=0);

dd=2 # pensacola
# first median to get it started
plot(c(block.years.center[1]-halfwidth, block.years.center[1]+halfwidth), rep(returnlevel.quantiles[[dd]][1,'q50'],2),
     type='l', lwd=3, col='black', xlim=c(1920,2000), ylim=c(0,21), xlab='', ylab='', las=1)
# now add the darker 25-75% range polygon before the median bars, ...
for (bb in 1:nblocks[[dd]]) {
    times.beg.end <- c(block.years.center[bb]-halfwidth, block.years.center[bb]+halfwidth)
    polygon(c(times.beg.end,rev(times.beg.end)), c(returnlevel.quantiles[[dd]][bb,c('q25','q25')],rev(returnlevel.quantiles[[dd]][bb,c('q75','q75')])),
            col=block.colors[bb], border=NA);
}
# ... and add the lighter 5-95% range polygon before the median bars too...
for (bb in 1:nblocks[[dd]]) {
    times.beg.end <- c(block.years.center[bb]-halfwidth, block.years.center[bb]+halfwidth)
    polygon(c(times.beg.end,rev(times.beg.end)), c(returnlevel.quantiles[[dd]][bb,c('q05','q05')],rev(returnlevel.quantiles[[dd]][bb,c('q95','q95')])),
            col=block.colors.lighter[bb], border=NA);
}
# ... so the bars are on top
for (bb in 1:nblocks[[dd]]) {lines(c(block.years.center[bb]-halfwidth, block.years.center[bb]+halfwidth),
                                   rep(returnlevel.quantiles[[dd]][bb,'q50'],2), lwd=3, col='black')}
# finally, add the NOLA data on the far right
#lines(c(2000,2000), c(-1000,1000), lty=3, lwd=2, col='black')
#text(2002, 7, 'New Orleans', srt=90, pos=4)
#times.beg.end <- c(block.years.center.nola-halfwidth, block.years.center.nola+halfwidth)
#polygon(c(times.beg.end,rev(times.beg.end)), c(returnlevel.quantiles.nola[c('q05','q05')],rev(returnlevel.quantiles.nola[c('q95','q95')])),
#       col='gray65', border=NA);
#polygon(c(times.beg.end,rev(times.beg.end)), c(returnlevel.quantiles.nola[c('q25','q25')],rev(returnlevel.quantiles.nola[c('q75','q75')])),
#       col='gray25', border=NA);
#lines(c(block.years.center.nola-halfwidth, block.years.center.nola+halfwidth),
#     rep(returnlevel.quantiles.nola['q50'],2), lwd=3, col='black')
text(1960, 20, 'Pensacola, Florida', pos=4)
mtext('Year', side=1, line=2.4, cex=1);
mtext('100-year return level [m]', side=2, line=2.2, cex=1);
mtext(side=3, text=expression(bold('   b')), line=-1, cex=1, adj=0);
dev.off()



##
##==============================================================================
## End
##==============================================================================
##
