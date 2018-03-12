##==============================================================================
## Plots and tables for BRICK scenarios paper (...)
##
## These routines also are provided as guidance for users to develop their own
## plotting and analysis routines.
##
## Questions? Tony Wong (twong@psu.edu)
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

setwd('~/codes/BRICK/R')

rm(list=ls())

## Initial set-up
library(ncdf4)

## File name for the BRICK physical model output (netCDF4)
filename.brick.allslr = '../output_model/BRICK_physical_allslr_08May2017.nc'

## File name for the Van Dantzig model output (netCDF4)
## Each of these also has x3 RCP scenarios, x2 storm surge scenarios
filename.vandantzig.nofd = '../output_model/VanDantzig_fd-none_2065_08May2017.nc'
filename.vandantzig.uniform = '../output_model/VanDantzig_fd-uniform_2065_08May2017.nc'
filename.vandantzig.gamma = '../output_model/VanDantzig_fd-gamma_2065_08May2017.nc'

## File name for the BRICK post-calibrated parameters (netcdf) (the BRICK output came from these guys)
filename.parameters.uniform = '../output_calibration/BRICK_postcalibratedParameters_fd-uniform_08May2017.nc'
filename.parameters.gamma = '../output_calibration/BRICK_postcalibratedParameters_fd-gamma_08May2017.nc'

## Other files
filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_07May2017.csv"
filename.fingerprints = "../fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"

## Other useful scripts
source('../Useful/colorblindPalette.R') # Get nice plotting colors: mycol array
source('../Useful/MultipleOutput.R')    # defines the useful ":=" operator

## And set the IPCC RCP colors
col26 <- c(0, 0, 255)/255
col45 <- c(121, 188, 255)/255
col60 <- c(255, 130, 45)/255
col85 <- c(255, 0, 0)/255

## Where would you like to save the plots?
plotdir='~/Box\ Sync/Wong-Projects/BRICK_scenarios/figures/'

##==============================================================================






##==============================================================================
##==============================================================================
## Grab the van Dantzig output and store, for all scenarios
## Also store the local sea level, and the surge.factor parameters for each AIS
## scenario.

# first level are RCP scenarios
gev.names <- c('location','scale','shape')
scen.rcp <- c('rcp26','rcp45','rcp60','rcp85'); n.rcp <- length(scen.rcp)
preturn = vector("list",n.rcp); names(preturn)=scen.rcp
lsl = vector("list",n.rcp); names(lsl)=scen.rcp
sf.level.all = vector("list",n.rcp); names(sf.level.all)=scen.rcp
sf.surge.all = vector("list",n.rcp); names(sf.surge.all)=scen.rcp
init = vector("list",n.rcp); names(init)=scen.rcp

# second level are AIS fast dynamics scenarios
scen.ais = c('none','gamma','uniform'); n.ais <- length(scen.ais)
scen.ss = c('st','ns'); n.ss <- length(scen.ss)
surge.factor <- vector('list',n.ss); names(surge.factor) <- scen.ss
for (rcp in scen.rcp) {
  preturn[[rcp]] = vector("list",n.ais); names(preturn[[rcp]])=scen.ais
  lsl[[rcp]] = vector("list",n.ais); names(lsl[[rcp]])=scen.ais
  sf.level.all[[rcp]] = vector("list",n.ais); names(sf.level.all[[rcp]])=scen.ais
  sf.surge.all[[rcp]] = vector("list",n.ais); names(sf.surge.all[[rcp]])=scen.ais
  init[[rcp]] = vector("list",n.ais); names(init[[rcp]])=scen.ais

# third level are storm surge scenarios
  for (ais in scen.ais) {
    preturn[[rcp]][[ais]] = vector("list",n.ss); names(preturn[[rcp]][[ais]])=scen.ss
    sf.level.all[[rcp]][[ais]] = vector("list",n.ss); names(sf.level.all[[rcp]][[ais]])=scen.ss
    sf.surge.all[[rcp]][[ais]] = vector("list",n.ss); names(sf.surge.all[[rcp]][[ais]])=scen.ss
  }
}

# read no fast dyanmics results
ncdata <- nc_open(filename.vandantzig.nofd)
  heightening <- ncvar_get(ncdata, 'H')
  surge.factor$ns <- ncvar_get(ncdata, 'surge_factor') # all RCPs have same surge factor
  surge.factor$st <- rep(0,length(surge.factor$ns))
  gev.stat <- ncvar_get(ncdata, 'gev_stat')
  VD.params <- ncvar_get(ncdata, 'VD_params')
  colnames(gev.stat) <- gev.names

  preturn$rcp26$none$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat')
  preturn$rcp45$none$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat')
  preturn$rcp60$none$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP60_stat')
  preturn$rcp85$none$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat')

  preturn$rcp26$none$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat')
  preturn$rcp45$none$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat')
  preturn$rcp60$none$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP60_nonstat')
  preturn$rcp85$none$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat')
nc_close(ncdata)

# read gamma fast dynamics prior results
ncdata <- nc_open(filename.vandantzig.gamma)
  #surge.factor <- ncvar_get(ncdata, 'surge_factor') # don't reread, same
  gev.stat <- ncvar_get(ncdata, 'gev_stat')
  colnames(gev.stat) <- gev.names

  preturn$rcp26$gamma$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat')
  preturn$rcp45$gamma$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat')
  preturn$rcp60$gamma$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP60_stat')
  preturn$rcp85$gamma$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat')

  preturn$rcp26$gamma$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat')
  preturn$rcp45$gamma$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat')
  preturn$rcp60$gamma$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP60_nonstat')
  preturn$rcp85$gamma$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat')
nc_close(ncdata)

# read uniform fast dynamics prior results
ncdata <- nc_open(filename.vandantzig.uniform)
  #surge.factor <- ncvar_get(ncdata, 'surge_factor') # don't reread
  gev.stat <- ncvar_get(ncdata, 'gev_stat')
  colnames(gev.stat) <- gev.names

  preturn$rcp26$uniform$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat')
  preturn$rcp45$uniform$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat')
  preturn$rcp60$uniform$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP60_stat')
  preturn$rcp85$uniform$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat')

  preturn$rcp26$uniform$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat')
  preturn$rcp45$uniform$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat')
  preturn$rcp60$uniform$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP60_nonstat')
  preturn$rcp85$uniform$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat')
nc_close(ncdata)

n.ensemble <- c(ncol(preturn$rcp85$none$ns) , ncol(preturn$rcp85$gamma$ns) , ncol(preturn$rcp85$uniform$ns) )
names(n.ensemble) <- scen.ais
n.height <- length(heightening)

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================

#install.packages("extRemes")
#install.packages("fExtremes")
#install.packages('ismev')
#install.packages('zoo')
library(extRemes)
library(fExtremes)
library(ismev)
library(zoo)
library(fields)

ncdata <- nc_open(filename.brick.allslr)
  mod.time <- ncvar_get(ncdata, 'time_proj')
  lsl$rcp26$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP26')
  lsl$rcp45$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP45')
  lsl$rcp60$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP60')
  lsl$rcp85$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP85')
  lsl$rcp26$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP26')
  lsl$rcp45$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP45')
  lsl$rcp60$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP60')
  lsl$rcp85$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP85')
  lsl$rcp26$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP26')
  lsl$rcp45$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP45')
  lsl$rcp60$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP60')
  lsl$rcp85$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP85')
nc_close(ncdata)

## distribution of subsidence
iproj <- which(mod.time==2065)
inorm <- which(mod.time==2015)
subs <- VD.params[,6] # m/year
slr.subs <- -subs*(mod.time[iproj]-mod.time[inorm])


# Use Simpson's rule for all integrals
# example: sintegral(x=pdf.x, fx=pdf.fit[1,])$value
library(Bolstad)

## Get the return levels for the stationary GEV storm surge case
surge.level <- n.ensemble['none'] # desired storm surge level

# empirical CDF
ecdf.vals  <- seq(from=0, to=1, length.out=(surge.level))

# empirical survival function
esf.vals <- 1-ecdf.vals

#n.lev <- 2^13  # largest power of two less than the ensemble size
n.lev <- length(ecdf.vals)  # determined by what our empirical survival functions can resolve

# sf.level will give the ensemble median level (sea level + surge + subs) for
# each of the empircal CDF levels (that we can resolve using our ensemble size).
# sealev is empirical CDF/SF and PDF for the scenario sea-level rise only
# surlev.avg is same as sf.level, but only surge (no sea level or subsidence)
sf.level <- init
sf.surge <- init
sf.sealev <- init
pdf.sealev <- init
pdf.surge  <- init
pdf.level  <- init

print('Warning: this step takes a while. More than ten minutes, less than an hour...')
for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    sf.level[[rcp]][[ais]] <- vector('list',n.ss); names(sf.level[[rcp]][[ais]]) <- scen.ss
    sf.surge[[rcp]][[ais]] <- vector('list',n.ss); names(sf.surge[[rcp]][[ais]]) <- scen.ss
    pdf.surge[[rcp]][[ais]] <- vector('list',n.ss); names(pdf.surge[[rcp]][[ais]]) <- scen.ss
    pdf.level[[rcp]][[ais]] <- vector('list',n.ss); names(pdf.level[[rcp]][[ais]]) <- scen.ss
    for (ss in scen.ss) {
      sf.level.all[[rcp]][[ais]][[ss]] <- mat.or.vec(n.ensemble[[ais]], n.lev)
      sf.surge.all[[rcp]][[ais]][[ss]] <- mat.or.vec(n.ensemble[[ais]], n.lev)
      # note that the GEV returns levels in mm, so /1000
      sf.level.all[[rcp]][[ais]][[ss]] <- t(sapply(1:n.ensemble[[ais]], function(sow) {(1+surge.factor[[ss]][sow])*(lsl[[rcp]][[ais]][iproj,sow]-lsl[[rcp]][[ais]][inorm,sow]) +
                                            slr.subs[sow] +
                                            0.001*qgev(ecdf.vals, xi=gev.stat[sow,'shape'], mu=gev.stat[sow,'location'], beta=gev.stat[sow,'scale'])}))
      sf.surge.all[[rcp]][[ais]][[ss]] <- t(sapply(1:n.ensemble[[ais]], function(sow) {surge.factor[[ss]][sow]*(lsl[[rcp]][[ais]][iproj,sow]-lsl[[rcp]][[ais]][inorm,sow]) +
                                            0.001*qgev(ecdf.vals, xi=gev.stat[sow,'shape'], mu=gev.stat[sow,'location'], beta=gev.stat[sow,'scale'])}))
      # SOW is the first dimension (row) of level, so apply median to each column
      sf.level[[rcp]][[ais]][[ss]] <- apply(sf.level.all[[rcp]][[ais]][[ss]], 2, median)
      sf.surge[[rcp]][[ais]][[ss]] <- apply(sf.surge.all[[rcp]][[ais]][[ss]], 2, median)

      # on a 16-GB RAM machine, you will probably kill your memory if you do not
      # clear out the old results for all of the end members (especially as the
      # ensembles become larger than 1000-10000...). you do not need them anyway
      sf.level.all[[rcp]][[ais]][[ss]] <- NULL
      sf.surge.all[[rcp]][[ais]][[ss]] <- NULL
    }
    sl.tmp <- lsl[[rcp]][[ais]][iproj,]-lsl[[rcp]][[ais]][inorm,]
    sf.sealev[[rcp]][[ais]] <- sl.tmp[order(sl.tmp)]
  }
}
# Save RData image file to resume later?
save.image(file = "BRICK_scenarios_analysis.RData")


lsl.lower <- 0
lsl.upper <- 5
lsl.n <- 2^9
lsl.breaks <- seq(from=lsl.lower, to=lsl.upper, length.out=lsl.n)
kern <- 'gaussian'

# fit kernel density estimate (default to Gaussian with 512 nodes) to the chunks
# of probability from the empirical survival function
# normalize with /sintegral(x=pdf.x, fx=tmp$y)$value
tmp <- density(x=sf.sealev[[rcp]][[ais]], from=lsl.lower, to=lsl.upper, n=lsl.n, kernel=kern)
pdf.x <- tmp$x
for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    tmp <- density(x=sf.sealev[[rcp]][[ais]], from=lsl.lower, to=lsl.upper, n=lsl.n, kernel=kern)
    pdf.sealev[[rcp]][[ais]] <- tmp$y/sintegral(x=pdf.x, fx=tmp$y)$value
    for (ss in scen.ss) {
      tmp <- density(x=sf.surge[[rcp]][[ais]][[ss]], from=lsl.lower, to=lsl.upper, n=lsl.n, kernel=kern)
      pdf.surge[[rcp]][[ais]][[ss]] <- tmp$y/sintegral(x=pdf.x, fx=tmp$y)$value
      tmp <- density(x=sf.level[[rcp]][[ais]][[ss]], from=lsl.lower, to=lsl.upper, n=lsl.n, kernel=kern)
      pdf.level[[rcp]][[ais]][[ss]] <- tmp$y/sintegral(x=pdf.x, fx=tmp$y)$value
    }
  }
}

## Return periods, from Van Dantzig output
# if overtopping is XX fraction of total failure probability, then
# p_fail_total = (1/XX)*p_fail_overtopping
# return_period = 1/p_fail_total = XX*return_period_overtopping
do.rcp <- c(1,2,4) # do only 2.6, 4.5 and 8.5? or all of them?
n.total <- length(do.rcp) * n.ais * n.ss
return.period <- mat.or.vec(n.total*n.ensemble['none'], 2); cnt <- 1
return.period.build <- mat.or.vec(n.total*n.ensemble['none'], 2); cnt <- 1
return.period.fragile60 <- mat.or.vec(n.total*n.ensemble['none'], 2); cnt <- 1
return.period.fragile80 <- mat.or.vec(n.total*n.ensemble['none'], 2); cnt <- 1
scen.names <- rep(NA,n.total)
rcp.tmp <- c('RCP2.6','RCP4.5','RCP6.0','RCP8.5'); rcp.tmp <- rcp.tmp[do.rcp]
names(rcp.tmp) <- scen.rcp[do.rcp]
for (rcp in scen.rcp[do.rcp]) {for (ais in scen.ais) {for (ss in scen.ss) {
    scen.names[cnt] <- paste(rcp.tmp[rcp],ais,ss)
    return.period[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 1] <- cnt
    return.period.build[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 1] <- cnt
    return.period.fragile60[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 1] <- cnt
    return.period.fragile80[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 1] <- cnt
    return.period[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- preturn[[rcp]][[ais]][[ss]][1,] # first height is 0 m build
    return.period.build[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- preturn[[rcp]][[ais]][[ss]][2,] # second height is 0.91 m build
    return.period.fragile60[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- (6/10)*preturn[[rcp]][[ais]][[ss]][1,]
    return.period.fragile80[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- (8/10)*preturn[[rcp]][[ais]][[ss]][1,]
    #return.period[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- rp[[rcp]][[ais]][[ss]]
    cnt <- cnt+1
}}}
colnames(return.period) <- c('Scenario','ReturnPeriod')
colnames(return.period.build) <- c('Scenario','ReturnPeriod')
colnames(return.period.fragile60) <- c('Scenario','ReturnPeriod')
colnames(return.period.fragile80) <- c('Scenario','ReturnPeriod')

##==============================================================================
##==============================================================================



save.image('../output_calibration/BRICK_analysis.RData')



##==============================================================================
## FIGURES
##==============================================================================






##==============================================================================

##=========     (column 1) pdfs of sea level rise and storm surge, all scenarios
## FIGURE 1     (column 2) corresponding survival functions
##=========

pdf(paste(plotdir,'pdfs_sf_slr_surge.pdf',sep=''),width=8,height=7,colormodel='cmyk')

par(mfrow=c(2,2), mai=c(.6,.63,.2,.26))
# (a) pdfs of local sea-level rise (2065)
plot(pdf.x, pdf.sealev$rcp26$none, type='l', xlim=c(0,1.1), ylim=c(0,12), lwd=1.5,
     col=rgb(col26[1],col26[2],col26[3]), lty=1,
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(pdf.x, pdf.sealev$rcp26$gamma, type='l', col=rgb(col26[1],col26[2],col26[3]), lty=2)
lines(pdf.x, pdf.sealev$rcp26$uniform, type='l', col=rgb(col26[1],col26[2],col26[3]), lty=3)
lines(pdf.x, pdf.sealev$rcp45$none, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=1)
lines(pdf.x, pdf.sealev$rcp45$gamma, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=2)
lines(pdf.x, pdf.sealev$rcp45$uniform, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=3)
lines(pdf.x, pdf.sealev$rcp85$none, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=1)
lines(pdf.x, pdf.sealev$rcp85$gamma, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=2)
lines(pdf.x, pdf.sealev$rcp85$uniform, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=3)

axis(1,seq(0,1,0.2),lab=c("0","0.2","0.4","0.6","0.8","1"), cex.axis=1.2)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1.2, cex=1);
mtext('Projected sea level in 2065 [m]', side=1, line=2.3, cex=1);
mtext(side=3, text=expression(bold('   a')), line=-1, cex=.9, adj=0);

legend(0.5,11.5,c("RCP2.6","RCP4.5","RCP8.5","no FD","FD, gamma","FD, uniform"),
       lty=c(1,1,1,1,2,3), lwd=2, cex=1.2,
       col=c(rgb(col26[1],col26[2],col26[3]),rgb(col45[1],col45[2],col45[3]),rgb(col85[1],col85[2],col85[3]),'black','black','black'),
       bty='n')

# (b) pdfs of storm surge
par(mai=c(.6,.63,.2,.26))
plot(pdf.x, pdf.surge$rcp26$none$st, type='l', xlim=c(0,4), ylim=c(0,3.5), lwd=1.5,
     col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lty=1,
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(pdf.x, pdf.surge$rcp26$none$ns, type='l', col=rgb(col26[1],col26[2],col26[3]), lty=1)
lines(pdf.x, pdf.surge$rcp26$gamma$ns, type='l', col=rgb(col26[1],col26[2],col26[3]), lty=2)
lines(pdf.x, pdf.surge$rcp26$uniform$ns, type='l', col=rgb(col26[1],col26[2],col26[3]), lty=3)
lines(pdf.x, pdf.surge$rcp45$none$ns, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=1)
lines(pdf.x, pdf.surge$rcp45$gamma$ns, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=2)
lines(pdf.x, pdf.surge$rcp45$uniform$ns, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=3)
lines(pdf.x, pdf.surge$rcp85$none$ns, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=1)
lines(pdf.x, pdf.surge$rcp85$gamma$ns, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=2)
lines(pdf.x, pdf.surge$rcp85$uniform$ns, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=3)

axis(1,seq(0,4,0.5),lab=c("0","","1","","2","","3","","4"), cex.axis=1.2)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1.2, cex=1);
mtext('Projected surge level in 2065 [m]', side=1, line=2.3, cex=1);
mtext(side=3, text=expression(bold('   b')), line=-1, cex=.9, adj=0);
text(.12,3,"stationary", pos=4, cex=1.3)
text(1.2,2.5,"non-stationary", pos=4, cex=1.3)

# (c) survival functions of sea-level rise
par(mai=c(.6,.63,.2,.26))
# rcp2.6
plot(sf.sealev$rcp26$none, log10(esf.vals), type='l', xlim=c(0,1), ylim=c(-3.2,0),
     lty=1, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i')
lines(sf.sealev$rcp26$gamma, log10(esf.vals), type='l',
      lty=2, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(sf.sealev$rcp26$uniform, log10(esf.vals), type='l',
      lty=3, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
# rcp4.5
lines(sf.sealev$rcp45$none, log10(esf.vals), type='l',
      lty=1, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(sf.sealev$rcp45$gamma, log10(esf.vals), type='l',
      lty=2, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(sf.sealev$rcp45$uniform, log10(esf.vals), type='l',
      lty=3, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
# rcp8.5
lines(sf.sealev$rcp85$none, log10(esf.vals), type='l',
      lty=1, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(sf.sealev$rcp85$gamma, log10(esf.vals), type='l',
      lty=2, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(sf.sealev$rcp85$uniform, log10(esf.vals), type='l',
      lty=3, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
# x-axis, ticks and text
axis(1,seq(0,1,0.2),lab=c("0","0.2","0.4","0.6","0.8","1"), cex.axis=1.2)
mtext('Projected sea level in 2065 [m]', side=1, line=2.3, cex=1);
# y-axis, ticks and text
axis(2, at=seq(-3,0), label=parse(text=paste("10^", seq(-3,0), sep="")), las=1, cex.axis=1.2, mgp=c(3,.7,0))
mtext('Survival function [1-CDF]', side = 2, line=2.6);
# add a panel label
mtext(side=3, text=expression(bold('   c')), line=-1.1, cex=.9, adj=0);
# add the 1/100, 1/500, 1/1000 levels
lines(c(-4,4),c(-2,-2),lty=2,col='black'); text(0.16,-1.85,"1/100 level", cex=1.2)
lines(c(-4,4),c(-log10(500),-log10(500)),lty=2,col='black'); text(0.16,-2.55,"1/500 level", cex=1.2);
#lines(c(-4,4),c(-3,-3),lty=2,col='black'); text(0.16,-2.85,"1/1000 level", cex=1.2)

# (d) survival functions of storm surge
par(mai=c(.6,.63,.2,.26))
# rcp2.6
plot(sf.surge$rcp26$none$st, log10(esf.vals), type='l', xlim=c(0,4), ylim=c(-2.2,0),
     lty=1, col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=1.5, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i')
lines(sf.surge$rcp26$none$ns, log10(esf.vals), type='l',
      lty=1, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(sf.surge$rcp26$gamma$ns, log10(esf.vals), type='l',
      lty=2, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(sf.surge$rcp26$uniform$ns, log10(esf.vals), type='l',
      lty=3, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
# rcp4.5
lines(sf.surge$rcp45$none$ns, log10(esf.vals), type='l',
      lty=1, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(sf.surge$rcp45$gamma$ns, log10(esf.vals), type='l',
      lty=2, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(sf.surge$rcp45$uniform$ns, log10(esf.vals), type='l',
      lty=3, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
# rcp8.5
lines(sf.surge$rcp85$none$ns, log10(esf.vals), type='l',
      lty=1, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(sf.surge$rcp85$gamma$ns, log10(esf.vals), type='l',
      lty=2, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(sf.surge$rcp85$uniform$ns, log10(esf.vals), type='l',
      lty=3, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
# x-axis, ticks and text
axis(1,seq(0,4,0.5),lab=c("0","","1","","2","","3","","4"), cex.axis=1.2)
mtext('Projected surge level in 2065 [m]', side=1, line=2.3, cex=1);
# y-axis, ticks and text
axis(2, at=seq(-3,0), label=parse(text=paste("10^", seq(-3,0), sep="")), las=1, cex.axis=1.2, mgp=c(3,.7,0))
mtext('Survival function [1-CDF]', side = 2, line=2.6);
# add a panel label
mtext(side=3, text=expression(bold('   d')), line=-1.1, cex=.9, adj=0);
# add the 1/100, 1/500, 1/1000 levels
lines(c(-4,4),c(-2,-2),lty=2,col='black'); text(0.6,-1.88,"1/100 level", cex=1.2)
#lines(c(-4,4),c(-log10(500),-log10(500)),lty=2,col='black'); text(0.6,-1.55,"1/500 level", cex=1.2);
#lines(c(-4,4),c(-3,-3),lty=2,col='black'); text(0.1,-2.85,"1/1,000 level", cex=1.2)

dev.off()

##==============================================================================








##==============================================================================

##=========   (row 1) pdfs of total sea+storm level, all scenarios
## FIGURE 2   (row 2) survival functions
##=========   (row 3) a legend

pdf(paste(plotdir,'pdfs_sf_sea+surge.pdf',sep=''),width=3.5,height=7,colormodel='cmyk')

par(mfrow=c(3,1), mai=c(.25,.63,.1,.2))
# (a) pdfs of sea-level rise + storm surge
# stationary
plot(pdf.x, pdf.level$rcp26$none$st, type='l', xlim=c(0,4), ylim=c(0,3), lwd=1.5,
     col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lty=1,
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(pdf.x, pdf.level$rcp26$gamma$st, type='l', col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lty=2)
lines(pdf.x, pdf.level$rcp26$uniform$st, type='l', col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lty=3)
lines(pdf.x, pdf.level$rcp45$none$st, type='l', col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lty=1)
lines(pdf.x, pdf.level$rcp45$gamma$st, type='l', col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lty=2)
lines(pdf.x, pdf.level$rcp45$uniform$st, type='l', col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lty=3)
lines(pdf.x, pdf.level$rcp85$none$st, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lty=1)
lines(pdf.x, pdf.level$rcp85$gamma$st, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lty=2)
lines(pdf.x, pdf.level$rcp85$uniform$st, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lty=3)
# non-stationary
lines(pdf.x, pdf.level$rcp26$none$ns, type='l', col=rgb(col26[1],col26[2],col26[3]), lty=1)
lines(pdf.x, pdf.level$rcp26$gamma$ns, type='l', col=rgb(col26[1],col26[2],col26[3]), lty=2)
lines(pdf.x, pdf.level$rcp26$uniform$ns, type='l', col=rgb(col26[1],col26[2],col26[3]), lty=3)
lines(pdf.x, pdf.level$rcp45$none$ns, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=1)
lines(pdf.x, pdf.level$rcp45$gamma$ns, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=2)
lines(pdf.x, pdf.level$rcp45$uniform$ns, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=3)
lines(pdf.x, pdf.level$rcp85$none$ns, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=1)
lines(pdf.x, pdf.level$rcp85$gamma$ns, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=2)
lines(pdf.x, pdf.level$rcp85$uniform$ns, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=3)

axis(1,seq(0,5,0.5),lab=c("0","","1","","2","","3","","4","","5"), cex.axis=1.3)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1.2, cex=.8);
mtext('Projected sea+surge level in 2065 [m]', side=1, line=2.4, cex=.8);
mtext(side=3, text=expression(bold('   a')), line=-1, cex=.9, adj=0);


# (b) survival functions of sea + storm surge level
par(mai=c(.05,.63,.3,.2))
# stationary
plot(sf.level$rcp26$none$st, log10(esf.vals), type='l', xlim=c(0,4), ylim=c(-2.2,0),
      lty=1, col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lwd=1.5, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i')
lines(sf.level$rcp26$gamma$st, log10(esf.vals), type='l',
      lty=2, col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lwd=1.5)
lines(sf.level$rcp26$uniform$st, log10(esf.vals), type='l',
      lty=3, col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lwd=1.5)
lines(sf.level$rcp45$none$st, log10(esf.vals), type='l',
      lty=1, col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lwd=1.5)
lines(sf.level$rcp45$gamma$st, log10(esf.vals), type='l',
      lty=2, col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lwd=1.5)
lines(sf.level$rcp45$uniform$st, log10(esf.vals), type='l',
      lty=3, col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lwd=1.5)
lines(sf.level$rcp85$none$st, log10(esf.vals), type='l',
      lty=1, col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=1.5)
lines(sf.level$rcp85$gamma$st, log10(esf.vals), type='l',
      lty=2, col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=1.5)
lines(sf.level$rcp85$uniform$st, log10(esf.vals), type='l',
      lty=3, col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=1.5)
# non-stationary
lines(sf.level$rcp26$none$ns, log10(esf.vals), type='l',
      lty=1, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(sf.level$rcp26$gamma$ns, log10(esf.vals), type='l',
      lty=2, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(sf.level$rcp26$uniform$ns, log10(esf.vals), type='l',
      lty=3, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(sf.level$rcp45$none$ns, log10(esf.vals), type='l',
      lty=1, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(sf.level$rcp45$gamma$ns, log10(esf.vals), type='l',
      lty=2, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(sf.level$rcp45$uniform$ns, log10(esf.vals), type='l',
      lty=3, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(sf.level$rcp85$none$ns, log10(esf.vals), type='l',
      lty=1, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(sf.level$rcp85$gamma$ns, log10(esf.vals), type='l',
      lty=2, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(sf.level$rcp85$uniform$ns, log10(esf.vals), type='l',
      lty=3, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)

# x-axis, ticks and text
axis(1,seq(0,5,0.5),lab=c("0","","1","","2","","3","","4","","5"), cex.axis=1.3)
mtext('Projected sea+surge level in 2065 [m]', side=1, line=2.6, cex=.8);
# y-axis, ticks and text
axis(2, at=seq(-3,0), label=parse(text=paste("10^", seq(-3,0), sep="")), las=1, cex.axis=1.3, mgp=c(3,.7,0))
mtext('Survival function [1-CDF]', side = 2, line=2.8, cex=0.8);
# add a panel label
mtext(side=3, text=expression(bold('   b')), line=-1.2, cex=.9, adj=0);
# add the 1/100, 1/500, 1/1000 levels
lines(c(-4,6),c(-2,-2),lty=2,col='black'); text(0.64,-1.88,"1/100 level", cex=1.3)
lines(c(-4,6),c(-log10(500),-log10(500)),lty=2,col='black'); text(0.48,-2.55,"1/500 level", cex=1.3);
#lines(c(-4,6),c(-3,-3),lty=2,col='black'); text(0.1,-2.85,"1/1000 level", cex=1.3)

par(mai=c(.05,.05,.05,.05))
# stationary
plot( 9, 9, xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE, xlim=c(0,1),ylim=c(0,1), xlab='', ylab='')

legend(0.15,0.75,c("RCP2.6, stationary surge","RCP4.5, stationary surge","RCP8.5, stationary surge","RCP2.6, non-stationary surge","RCP4.5, non-stationary surge","RCP8.5, non-stationary surge","no FD","FD, gamma","FD, uniform"),
       lty=c(1,1,1,1,1,1,1,2,3), lwd=2, cex=1.2,
       col=c(rgb(mycol[10,1],mycol[10,2],mycol[10,3]),rgb(mycol[8,1],mycol[8,2],mycol[8,3]),rgb(mycol[6,1],mycol[6,2],mycol[6,3]),rgb(col26[1],col26[2],col26[3]),rgb(col45[1],col45[2],col45[3]),rgb(col85[1],col85[2],col85[3]),'black','black','black'),
       bty='n')

dev.off()

##==============================================================================









##==============================================================================

##=========   box-whisker plots of return periods with no building.
## FIGURE 3   center all around 100-year, use log10-scale.
##=========   Note 1: post-processing done in Google Slides. A copy of this
##            document (with edit access, because I kept a different copy for
##            myself :-) ) is available at the URL below:
##            https://docs.google.com/presentation/d/1OdEt-FhmUC_xo0QtgWHdODRfvwAnR7p7EUG0NbHxMCg/edit?usp=sharing
##            Note 2: to be nice to others who may find this useful, make yourself
##            a copy instead of modifying the slides at the link above.

pdf(paste(plotdir,'returnperiods.pdf',sep=''),width=5,height=5.5,colormodel='cmyk')

b1 <- boxplot(ReturnPeriod~Scenario,data=return.period, horizontal=TRUE, log='x',
              xlab="", ylab="", outline=FALSE, ylim=c(8,10000), xaxt='n', yaxt='n', plot=FALSE)
tmp <- data.frame(cbind(1:n.total, b1$stats[3,])); colnames(tmp) <- c('scen','retp')
tmp2 <- tmp[order(tmp$retp),]
rp.order <- tmp2[,1]
return.period2 <- return.period
for (i in 1:n.total) {
    #return.period2[((i-1)*n.ensemble['none']+1):(i*n.ensemble['none']),1] <- rp.order[i]
    return.period2[((i-1)*n.ensemble['none']+1):(i*n.ensemble['none']),2] <- return.period[((rev(rp.order)[i]-1)*n.ensemble['none']+1):((rev(rp.order)[i])*n.ensemble['none']),2]
}
scen.names2 <- scen.names[rev(rp.order)]

par(mfrow=c(1,1), mai=c(.6,1.66,.1,.05))
boxplot(ReturnPeriod~Scenario,data=return.period2, horizontal=TRUE, log='x',
        xlab="", ylab="", outline=FALSE, ylim=c(8,75000), xaxt='n', yaxt='n', xaxs='i',xlim=c(1,21))
axis(1,c(1,10,1e2,1e3,1e4),lab=c("1","10","100","1,000","10,000"), cex.axis=1.1, mgp=c(3,.6,.0))
mtext('Return period [years]', side=1, line=1.9, cex=1.1)
axis(2,seq(1,n.total),lab=scen.names2, cex.axis=1.1, las=1, mgp=c(3,.6,0))
mtext('Scenario', side=3, line=1.3, cex=1.1, adj=-.3)
lines(c(100,100),c(-100,100),lty=1, lwd=1.5, col=rgb(col85[1],col85[2],col85[3]))   # Master Plan general safety
lines(c(500,500),c(-100,20.5),lty=5, lwd=1.5, col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]))   # Master Plan critical infrastructure
lines(c(5000,5000),c(-100,19),lty=6, lwd=1.5, col=rgb(mycol[3,1],mycol[3,2],mycol[3,3])) # Dutch Perspective economically-efficient

dev.off()



##==============================================================================
## Analysis of the high-risk upper tails of the probabilistic return period
## estimates:
thing1 <- preturn$rcp85$none$ns[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP8.5, no AIS FD, non-stationary storm surge that miss 500-year protection: ',ecdf.vals[which.min(abs(thing2-500))],sep=''))

thing1 <- preturn$rcp85$gamma$ns[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP8.5, gamma FD priors, non-stationary storm surge that miss 500-year protection: ',ecdf.vals[which.min(abs(thing2-500))],sep=''))

thing1 <- preturn$rcp85$uniform$ns[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP8.5, uniform FD priors, non-stationary storm surge that miss 500-year protection: ',ecdf.vals[which.min(abs(thing2-500))],sep=''))

## ... best/worst case SOW that miss the 100-year protection standard:
thing1 <- preturn$rcp26$none$st[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP2.6, no AIS FD, stationary storm surge that miss 500-year protection: ',ecdf.vals[which.min(abs(thing2-500))],sep=''))

thing1 <- preturn$rcp85$gamma$ns[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP8.5, gamma FD priors, non-stationary storm surge that miss 500-year protection: ',ecdf.vals[which.min(abs(thing2-500))],sep=''))

##==============================================================================
## ... best/worst case SOW that miss the 100-year protection standard:
thing1 <- preturn$rcp26$uniform$st[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP2.6, no AIS FD, stationary storm surge that miss 100-year protection: ',ecdf.vals[which.min(abs(thing2-100))],sep=''))

thing1 <- preturn$rcp85$uniform$ns[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP8.5, gamma FD priors, non-stationary storm surge that miss 100-year protection: ',ecdf.vals[which.min(abs(thing2-100))],sep=''))

## ... best/worst case SOW that miss the 500-year protection standard:
thing1 <- preturn$rcp26$none$st[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP2.6, no AIS FD, stationary storm surge that miss 500-year protection: ',ecdf.vals[which.min(abs(thing2-500))],sep=''))

thing1 <- preturn$rcp85$uniform$ns[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP8.5, gamma FD priors, non-stationary storm surge that miss 500-year protection: ',ecdf.vals[which.min(abs(thing2-500))],sep=''))

##==============================================================================
## ... best/worst case SOW that fulfill the 100-year protection standard:
thing1 <- preturn$rcp26$uniform$st[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP2.6, no AIS FD, stationary storm surge that fulfill 100-year protection: ',esf.vals[which.min(abs(thing2-100))],sep=''))

thing1 <- preturn$rcp85$uniform$ns[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP8.5, gamma FD priors, non-stationary storm surge that fulfill 100-year protection: ',esf.vals[which.min(abs(thing2-100))],sep=''))

## ... best/worst case SOW that fulfill the 500-year protection standard:
thing1 <- preturn$rcp26$none$st[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP2.6, no AIS FD, stationary storm surge that fulfill 500-year protection: ',esf.vals[which.min(abs(thing2-500))],sep=''))

thing1 <- preturn$rcp85$uniform$ns[1,]
thing2 <- thing1[order(thing1)]
print(paste('Fraction of SOW in RCP8.5, uniform FD priors, non-stationary storm surge that fulfill 500-year protection: ',esf.vals[which.min(abs(thing2-500))],sep=''))
##==============================================================================








##==============================================================================

##==========   same as figure 3, but with additional 3 feet of heightening
## FIGURE S1
##==========

pdf(paste(plotdir,'returnperiods_build.pdf',sep=''),width=5,height=5.5,colormodel='cmyk')

b1 <- boxplot(ReturnPeriod~Scenario,data=return.period.build, horizontal=TRUE, log='x',
              xlab="", ylab="", outline=FALSE, ylim=c(8,10000), xaxt='n', yaxt='n', plot=FALSE)
tmp <- data.frame(cbind(1:n.total, b1$stats[3,])); colnames(tmp) <- c('scen','retp')
tmp2 <- tmp[order(tmp$retp),]
rp.order <- tmp2[,1]
return.period.build2 <- return.period.build
for (i in 1:n.total) {
    return.period.build2[((i-1)*n.ensemble['none']+1):(i*n.ensemble['none']),2] <- return.period.build[((rev(rp.order)[i]-1)*n.ensemble['none']+1):((rev(rp.order)[i])*n.ensemble['none']),2]
}
scen.names2 <- scen.names[rev(rp.order)]

par(mfrow=c(1,1), mai=c(.6,1.66,.1,.05))
boxplot(ReturnPeriod~Scenario,data=return.period.build2, horizontal=TRUE, log='x',
        xlab="", ylab="", outline=FALSE, ylim=c(8,75000), xaxt='n', yaxt='n', xaxs='i',xlim=c(1,21))
axis(1,c(1,10,1e2,1e3,1e4),lab=c("1","10","100","1,000","10,000"), cex.axis=1.1, mgp=c(3,.6,.0))
mtext('Return period [years]', side=1, line=1.9, cex=1.1)
axis(2,seq(1,n.total),lab=scen.names2, cex.axis=1.1, las=1, mgp=c(3,.6,0))
mtext('Scenario', side=3, line=1.3, cex=1.1, adj=-.3)
lines(c(100,100),c(-100,100),lty=1, lwd=1.5, col=rgb(col85[1],col85[2],col85[3]))   # Master Plan general safety
lines(c(500,500),c(-100,20.5),lty=5, lwd=1.5, col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]))   # Master Plan critical infrastructure
lines(c(5000,5000),c(-100,19),lty=6, lwd=1.5, col=rgb(mycol[3,1],mycol[3,2],mycol[3,3])) # Dutch Perspective economically-efficient

dev.off()

##==============================================================================

##==========   same as figure 3, but with structural fragility at 60% overtopping
## FIGURE S2   (set above) instead of assuming 100% overtopping
##==========

pdf(paste(plotdir,'returnperiods_fragile60.pdf',sep=''),width=5,height=5.5,colormodel='cmyk')

b1 <- boxplot(ReturnPeriod~Scenario,data=return.period.fragile60, horizontal=TRUE, log='x',
              xlab="", ylab="", outline=FALSE, ylim=c(8,10000), xaxt='n', yaxt='n', plot=FALSE)
tmp <- data.frame(cbind(1:n.total, b1$stats[3,])); colnames(tmp) <- c('scen','retp')
tmp2 <- tmp[order(tmp$retp),]
rp.order <- tmp2[,1]
return.period.fragile602 <- return.period.fragile60
for (i in 1:n.total) {
    return.period.fragile602[((i-1)*n.ensemble['none']+1):(i*n.ensemble['none']),2] <- return.period.fragile60[((rev(rp.order)[i]-1)*n.ensemble['none']+1):((rev(rp.order)[i])*n.ensemble['none']),2]
}
scen.names2 <- scen.names[rev(rp.order)]

par(mfrow=c(1,1), mai=c(.6,1.66,.1,.05))
boxplot(ReturnPeriod~Scenario,data=return.period.fragile602, horizontal=TRUE, log='x',
        xlab="", ylab="", outline=FALSE, ylim=c(8,75000), xaxt='n', yaxt='n', xaxs='i',xlim=c(1,21))
axis(1,c(1,10,1e2,1e3,1e4),lab=c("1","10","100","1,000","10,000"), cex.axis=1.1, mgp=c(3,.6,.0))
mtext('Return period [years]', side=1, line=1.9, cex=1.1)
axis(2,seq(1,n.total),lab=scen.names2, cex.axis=1.1, las=1, mgp=c(3,.6,0))
mtext('Scenario', side=3, line=1.3, cex=1.1, adj=-.3)
lines(c(100,100),c(-100,100),lty=1, lwd=1.5, col=rgb(col85[1],col85[2],col85[3]))   # Master Plan general safety
lines(c(500,500),c(-100,20.5),lty=5, lwd=1.5, col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]))   # Master Plan critical infrastructure
lines(c(5000,5000),c(-100,19),lty=6, lwd=1.5, col=rgb(mycol[3,1],mycol[3,2],mycol[3,3])) # Dutch Perspective economically-efficient

dev.off()

##==============================================================================

##==========   same as figure 3, but with structural fragility at 80% overtopping
## FIGURE S3   (set above) instead of assuming 100% overtopping
##==========

pdf(paste(plotdir,'returnperiods_fragile80.pdf',sep=''),width=5,height=5.5,colormodel='cmyk')

b1 <- boxplot(ReturnPeriod~Scenario,data=return.period.fragile80, horizontal=TRUE, log='x',
              xlab="", ylab="", outline=FALSE, ylim=c(8,10000), xaxt='n', yaxt='n', plot=FALSE)
tmp <- data.frame(cbind(1:n.total, b1$stats[3,])); colnames(tmp) <- c('scen','retp')
tmp2 <- tmp[order(tmp$retp),]
rp.order <- tmp2[,1]
return.period.fragile802 <- return.period.fragile80
for (i in 1:n.total) {
    return.period.fragile802[((i-1)*n.ensemble['none']+1):(i*n.ensemble['none']),2] <- return.period.fragile80[((rev(rp.order)[i]-1)*n.ensemble['none']+1):((rev(rp.order)[i])*n.ensemble['none']),2]
}
scen.names2 <- scen.names[rev(rp.order)]

par(mfrow=c(1,1), mai=c(.6,1.66,.1,.05))
boxplot(ReturnPeriod~Scenario,data=return.period.fragile802, horizontal=TRUE, log='x',
        xlab="", ylab="", outline=FALSE, ylim=c(8,75000), xaxt='n', yaxt='n', xaxs='i',xlim=c(1,21))
axis(1,c(1,10,1e2,1e3,1e4),lab=c("1","10","100","1,000","10,000"), cex.axis=1.1, mgp=c(3,.6,.0))
mtext('Return period [years]', side=1, line=1.9, cex=1.1)
axis(2,seq(1,n.total),lab=scen.names2, cex.axis=1.1, las=1, mgp=c(3,.6,0))
mtext('Scenario', side=3, line=1.3, cex=1.1, adj=-.3)
lines(c(100,100),c(-100,100),lty=1, lwd=1.5, col=rgb(col85[1],col85[2],col85[3]))   # Master Plan general safety
lines(c(500,500),c(-100,20.5),lty=5, lwd=1.5, col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]))   # Master Plan critical infrastructure
lines(c(5000,5000),c(-100,19),lty=6, lwd=1.5, col=rgb(mycol[3,1],mycol[3,2],mycol[3,3])) # Dutch Perspective economically-efficient

dev.off()

##==============================================================================

##==========
## FIGURE 4   control Sobol sensitivity analysis, radial convergence plot
## FIGURE S4  Sobol plot with h0 and c removed (AIS runoff line height parameters)
## FIGURE S5  Sobol plot with GEV parameters removed
##==========

## Plots for Sobol sensitivity analysis for drivers of flood risk, plus SOM

## These are generated separately in BRICK_Sobol_plotting.R

##==============================================================================

##==========
## FIGURE S6  Generated using stormsurge_sensitivity_experiment_analysis.R
## FIGURE S7  Generated using stormsurge_sensitivity_experiment_analysis.R
##==========

##==============================================================================




##==============================================================================
## End
##==============================================================================
