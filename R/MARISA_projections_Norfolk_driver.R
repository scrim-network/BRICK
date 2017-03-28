# =======================================================================================
# Create projections of local sea level and storm surge for Sewall Point.
# Use storm surge projections a la Grinsted et al 2013, with T=temperature
# anomaly relative to 1980-2000 average.
# And temperature projections from DOECLIM (within BRICK v0.1) under RCP2.6,
# 4.5 and 8.5.
#
# NOTE -- NOT CURRENTLY SUPPORTING TIDE GAUGE DATA AND STORM SURGE PROJECTION
#           FOR SEWALL PT --
#
# =======================================================================================
#
#   Required settings (define below):
# - lat.proj    latitude of point at which you want local sea level (currently only support 1 pt)
# - lon.proj    longitude of point at which you want local sea level
# - dat.dir     directory with tide gauge data (m) in the form [date | time | sl]
#               where date is YYYYMMDD, time is hh:mm, sl is meters, and these
#               are the column names (for referencing easily)
# - l.nonstat   which gev parameters to be non-stationary? these will covary
#               with temperature anomaly relative to 1980-2000 mean.
#               excepts [location | shape | scale]. Note that Coles et al. 2001
#               (book: An Introduction to Statistical Modeling of Extreme Values)
#               suggest to leave the scale parameter stationary.
# - filename.brick  BRICK output file to base projections on
# - niter.gev   number of MCMC iterations to use to estimate (non-)stationary
#               GEV parameters
#
#   Simulates (output variables):
# - lsl.out     local sea level rise (m sle)
# - gev.out     ...?
# - output.file ...?
#
#   Parameters:
# - filename.fingerprints = "../fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"
#               sea level fingerprinting data (Slangen et al 2014)
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

# set things here to get different site projections

# Note: you should interactively run the tide gauge/storm surge projections.
# This should NOT be routinized. It is not a "one-size-fits-all" type of
# calculation.

dat.dir <- '../data/tidegauge_SewellsPoint/'
plotdir <- '~/Box\ Sync/Wong-Projects/MARISSA/figures/'
filetype <- 'txt'
septype <- '\t'

lat.proj <- 36+(56.8/60)          # Sewall Point tide gauge at 36° 56.8' N
lon.proj <- -(76+(19.8/60))       # 76° 19.8' W (NOAA Tides and Currents website)
scen.rcp <- c('rcp26','rcp45','rcp85')

l.nonstat <- vector('list',3); names(l.nonstat) <- c('location','shape','scale')
l.nonstat$location <- TRUE
l.nonstat$shape <- FALSE
l.nonstat$scale <- FALSE
niter.gev <- 1e4
burnin <- 0.5

filename.brick <- '../output_model/BRICK-fastdyn_physical_gamma_31Jan2017.nc'

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.projout <- paste('../output_model/BRICK_project-lsl-surge_Norfolk_',today,'.nc',sep='')
filename.lslout  <- paste('../output_model/BRICK_project-lsl_Norfolk_',today,'.csv', sep="")

# unless other data sets are quite different, below here likely does not need
# to be modified
##==============================================================================


##==============================================================================
## Read tide gauge data
files.tg <- list.files(path=dat.dir,pattern=filetype)

data <- read.table(paste(dat.dir,files.tg[1],sep=''), header = TRUE, sep=septype)
if(length(files.tg) > 1) {
    for (ff in 2:length(files.tg)) {
        data <- rbind(data, read.table(paste(dat.dir,files.tg[ff],sep=''), header = TRUE, sep=septype))
    }
}
##==============================================================================


##==============================================================================
## Make local sea level projections

#install.packages('ncdf4')
library(ncdf4)

# Initialize lists for projections
temp <- vector('list', length(scen.rcp)); names(temp) <- scen.rcp
slr_gis <- temp
slr_ais <- temp
slr_gsic <- temp
slr_te <- temp
gmsl <- temp
lsl_proj <- temp

# Read sea level output
ncdata <- nc_open(filename.brick)
  for (rcp in scen.rcp) {
    if(rcp=='rcp26') {
        temp[[rcp]]     <- ncvar_get(ncdata, 'temp_RCP26')
        slr_gis[[rcp]]  <- ncvar_get(ncdata, 'GIS_RCP26')
        slr_ais[[rcp]]  <- ncvar_get(ncdata, 'AIS_RCP26')
        slr_gsic[[rcp]] <- ncvar_get(ncdata, 'GSIC_RCP26')
        slr_te[[rcp]]   <- ncvar_get(ncdata, 'TE_RCP26')
        gmsl[[rcp]]     <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
    } else if(rcp=='rcp45') {
        temp[[rcp]]     <- ncvar_get(ncdata, 'temp_RCP45')
        slr_gis[[rcp]]  <- ncvar_get(ncdata, 'GIS_RCP45')
        slr_ais[[rcp]]  <- ncvar_get(ncdata, 'AIS_RCP45')
        slr_gsic[[rcp]] <- ncvar_get(ncdata, 'GSIC_RCP45')
        slr_te[[rcp]]   <- ncvar_get(ncdata, 'TE_RCP45')
        gmsl[[rcp]]     <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
    } else if(rcp=='rcp60') {
        temp[[rcp]]     <- ncvar_get(ncdata, 'temp_RCP60')
        slr_gis[[rcp]]  <- ncvar_get(ncdata, 'GIS_RCP60')
        slr_ais[[rcp]]  <- ncvar_get(ncdata, 'AIS_RCP60')
        slr_gsic[[rcp]] <- ncvar_get(ncdata, 'GSIC_RCP60')
        slr_te[[rcp]]   <- ncvar_get(ncdata, 'TE_RCP60')
        gmsl[[rcp]]     <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP60')
    } else if(rcp=='rcp85') {
        temp[[rcp]]     <- ncvar_get(ncdata, 'temp_RCP85')
        slr_gis[[rcp]]  <- ncvar_get(ncdata, 'GIS_RCP85')
        slr_ais[[rcp]]  <- ncvar_get(ncdata, 'AIS_RCP85')
        slr_gsic[[rcp]] <- ncvar_get(ncdata, 'GSIC_RCP85')
        slr_te[[rcp]]   <- ncvar_get(ncdata, 'TE_RCP85')
        gmsl[[rcp]]     <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
    } else {
        print(paste('Error - unrecognized RCP scenario: ',rcp,sep=''))
    }
  }
  t.proj     <- ncvar_get(ncdata, 'time_proj')
  n.ensemble <- ncol(temp[[1]])
nc_close(ncdata)

# Routine for projecting local sea level
source('../R/BRICK_LSL.R')

# Get LSL projections for each RPC scenario
for (rcp in scen.rcp) {
    lsl_proj[[rcp]] = brick_lsl(lat.in   = lat.proj,
                                lon.in   = lon.proj,
                                n.time   = length(t.proj),
                                slr_gis  = slr_gis[[rcp]],
                                slr_gsic = slr_gsic[[rcp]],
						        slr_ais  = slr_ais[[rcp]],
						        slr_te   = slr_te[[rcp]])
}
##==============================================================================


##==============================================================================
## Write CSV with monthly projections

# normalize realtive to 1983-2001
inorm <- which(t.proj==1983):which(t.proj==2001)
lsl26 <- lsl_proj$rcp26 - t(replicate(nrow(lsl_proj$rcp26),apply(lsl_proj$rcp26[inorm,],2,mean)))
lsl45 <- lsl_proj$rcp45 - t(replicate(nrow(lsl_proj$rcp45),apply(lsl_proj$rcp45[inorm,],2,mean)))
lsl85 <- lsl_proj$rcp85 - t(replicate(nrow(lsl_proj$rcp85),apply(lsl_proj$rcp85[inorm,],2,mean)))

lsl.out <- cbind(t.proj,
                 apply(lsl26,1,median),
                 apply(lsl45,1,median),
                 apply(lsl85,1,median))
colnames(lsl.out) <- c('Year','RCP26','RCP45','RCP85')
write.table(lsl.out, file=filename.lslout, sep=",", qmethod="double", row.names=FALSE)

##==============================================================================





##==============================================================================
## Plot --  akin to NOAA 2017 GMSL/US coastal report

## round up all of the GMSL projections, and LSL (already fingerprinted, so
## no need to screw around with the slr_xxx)
ns <- length(gmsl)
nt <- nrow(gmsl[[1]])
ne <- ncol(gmsl[[1]])
gmsl.all <- mat.or.vec(nt,ns*ne)
lsl.all  <- mat.or.vec(nt,ns*ne)
for (i in 1:ns) {
    gmsl.all[,((i-1)*ne+1):(i*ne)] <- gmsl[[i]]
    lsl.all[,((i-1)*ne+1):(i*ne)]  <- lsl_proj[[i]]
}

## normalize relative to 1991-2009 mean
ind.norm <- which(t.proj==1991):which(t.proj==2009)
gmsl.all <- gmsl.all - t(replicate(nrow(gmsl.all),apply(gmsl.all[ind.norm,],2,mean)))
lsl.all  <- lsl.all  - t(replicate(nrow(lsl.all) ,apply(lsl.all[ind.norm,] ,2,mean)))

## bin into those that hit the targets from NOAA 2017 report for the 6 scenarios
## 1. low                 28-32 cm GMSL rise in 2100
## 2. intermediate-low    48-52
## 3. intermediate        98-102
## 4. intermediate-high   145-155
## 5. high                195-205
## 6. extreme             235-265

scen.slr <- c('Low','Intermediate-low','Intermediate','Intermediate-high','High','Extreme')
windows <- data.frame( cbind( c(.28, .48,  .98, 1.45, 1.95, 2.35)   # lower bounds
                             ,c(.32, .52, 1.02, 1.55, 2.05, 2.65))) # upper bounds
windows.low <- data.frame( cbind( c(2030, 2050, 2100),
                                  c( .08,  .13,  .28),
                                  c( .10,  .17,  .32) ))
i2030 <- which(t.proj==2030)
i2050 <- which(t.proj==2050)
i2100 <- which(t.proj==2100)
ind.bin <- vector('list',length(scen.slr))
for (i in 1:length(ind.bin)) {
#    if(scen.slr[i]=='Low') {
#        ind.bin[[i]] <- which( gmsl.all[i2030] >= windows.low[1,2] &
#                               gmsl.all[i2030] <= windows.low[1,3] &
#                               gmsl.all[i2050] >= windows.low[2,2] &
#                               gmsl.all[i2050] <= windows.low[2,3] &
#                               gmsl.all[i2100] >= windows.low[3,2] &
#                               gmsl.all[i2100] <= windows.low[3,3] )
#    } else {
        ind.bin[[i]] <- which( gmsl.all[i2100,] >= windows[i,1] &
                               gmsl.all[i2100,] <= windows[i,2] )
#    }
}

## pick out the LSL projections that fit each scenario
lsl.scen <- vector('list',length(scen.slr)); names(lsl.scen) <- scen.slr
for (i in 1:length(scen.slr)) {
    lsl.scen[[i]] <- lsl.all[,ind.bin[[i]]]
}

## plot
source('../Useful/colorblindPalette.R')

plot(t.proj, apply(lsl.scen[[6]], 1, mean), type='l', lwd=2,
     xlim=c(1960, 2060), ylim=c(-.4,1.2), col=mycol.rgb[11],
     xlab='Year', ylab='Regional sea level [meters]')
lines(t.proj, apply(lsl.scen[[5]], 1, mean), lwd=2, col=mycol.rgb[13])
lines(t.proj, apply(lsl.scen[[4]], 1, mean), lwd=2, col=mycol.rgb[15])
lines(t.proj, apply(lsl.scen[[3]], 1, mean), lwd=2, col=mycol.rgb[9])
lines(t.proj, apply(lsl.scen[[2]], 1, mean), lwd=2, col=mycol.rgb[2])
lines(t.proj, apply(lsl.scen[[1]], 1, mean), lwd=2, col=mycol.rgb[6])

## Medians in each NOAA report scenario throughout the next century
## (Table 5/6, NOAA report)
time.noaa <- seq(from=2000, to=2100, by=10)
gmsl.noaa <- vector('list',6)
names(gmsl.noaa) <- c('low','intlow','int','inthigh','high','extreme')
gmsl.noaa$low     <- c(0, 0.03, 0.06, 0.09, 0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.30)
gmsl.noaa$intlow  <- c(0, 0.04, 0.08, 0.13, 0.18, 0.24, 0.29, 0.35, 0.4, 0.45, 0.50)
gmsl.noaa$int     <- c(0, 0.04, 0.10, 0.16, 0.25, 0.34, 0.45, 0.57, 0.71, 0.85, 1.0)
gmsl.noaa$inthigh <- c(0, 0.05, 0.10, 0.19, 0.30, 0.44, 0.60, 0.79, 1.0, 1.2, 1.5)
gmsl.noaa$high    <- c(0, 0.05, 0.11, 0.21, 0.36, 0.54, 0.77, 1.0, 1.3, 1.7, 2.0)
gmsl.noaa$extreme <- c(0, 0.04, 0.11, 0.24, 0.41, 0.63, 0.90, 1.2, 1.6, 2.0, 2.5)
##==============================================================================

##==============================================================================
## Plot -- pdf and survival functions

## And set the IPCC RCP colors
col26 <- c(0, 0, 255)/255
col45 <- c(121, 188, 255)/255
col60 <- c(255, 130, 45)/255
col85 <- c(255, 0, 0)/255

i2100 <- which(t.proj==2100)

## Normalize relative to 1991-2009, as in NOAA 2017 report projections
ind.norm <- which(t.proj==1991):which(t.proj==2009)
lsl_proj$rcp26 <- lsl_proj$rcp26 - t(replicate(nrow(lsl_proj$rcp26) ,apply(lsl_proj$rcp26[ind.norm,] ,2,mean)))
lsl_proj$rcp45 <- lsl_proj$rcp45 - t(replicate(nrow(lsl_proj$rcp45) ,apply(lsl_proj$rcp45[ind.norm,] ,2,mean)))
lsl_proj$rcp85 <- lsl_proj$rcp85 - t(replicate(nrow(lsl_proj$rcp85) ,apply(lsl_proj$rcp85[ind.norm,] ,2,mean)))

slr2100.rcp26 <- lsl_proj$rcp26[i2100,]
slr2100.rcp45 <- lsl_proj$rcp45[i2100,]
slr2100.rcp85 <- lsl_proj$rcp85[i2100,]

n.pdf <- 2^11       # 2^11 is largest power of 2 less than # ensemble members
pdf.slr2100.rcp26 <- density(slr2100.rcp26,from=0,to=4,na.rm=TRUE, n=n.pdf)
pdf.slr2100.rcp45 <- density(slr2100.rcp45,from=0,to=4,na.rm=TRUE, n=n.pdf)
pdf.slr2100.rcp85 <- density(slr2100.rcp85,from=0,to=4,na.rm=TRUE, n=n.pdf)

cdf.slr2100.rcp26 <- rep(0,length(pdf.slr2100.rcp26$x))
cdf.slr2100.rcp45 <- rep(0,length(pdf.slr2100.rcp45$x))
cdf.slr2100.rcp85 <- rep(0,length(pdf.slr2100.rcp85$x))

# empirical CDF
ecdf.vals  <- seq(from=0, to=1, length.out=length(slr2100.rcp26))
ecdf.rcp26 <- slr2100.rcp26[order(slr2100.rcp26)]
ecdf.rcp45 <- slr2100.rcp45[order(slr2100.rcp45)]
ecdf.rcp85 <- slr2100.rcp85[order(slr2100.rcp85)]

# empirical survival function
esf.vals <- 1-ecdf.vals

pdf(paste(plotdir,'distributions_SLR2100_pdf+sf_Norfolk.pdf',sep=''),width=5,height=7,colormodel='cmyk')
par(mfrow=c(2,1), mai=c(.85,.74,.1,.15))

plot(pdf.slr2100.rcp85$x,pdf.slr2100.rcp85$y, type='l', xlim=c(0,3), ylim=c(0,5.2), lty=1,
     col=rgb(col85[1],col85[2],col85[3]), lwd=3, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i',axes=FALSE);
axis(1,seq(0,3,0.5),lab=c("0","0.5","1","1.5","2","2.5","3"))
u <- par("usr")
arrows(0, u[3],0, u[4], code = 2, xpd = TRUE)
mtext('Probability density', side=2, line=1.3);
mtext('Projected local mean sea level in 2100\nrelative to 1991-2009 average [m]', side=1, line=3);
mtext(side=3, text=expression(bold('   a')), line=-1, cex=.9, adj=0);

lines(pdf.slr2100.rcp45$x,pdf.slr2100.rcp45$y, type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=3);
lines(pdf.slr2100.rcp26$x,pdf.slr2100.rcp26$y, type='l', col=rgb(col26[1],col26[2],col26[3]), lwd=3);

# add NOAA scenarios as shaded regions/bars
lines(c(gmsl.noaa$low[11],gmsl.noaa$low[11]),         c(-10,10), lty=2, lwd=2, col=mycol.rgb[11])
lines(c(gmsl.noaa$intlow[11],gmsl.noaa$intlow[11]),   c(-10,10), lty=2, lwd=2, col=mycol.rgb[13])
lines(c(gmsl.noaa$int[11],gmsl.noaa$int[11]),         c(-10,10), lty=2, lwd=2, col=mycol.rgb[15])
lines(c(gmsl.noaa$inthigh[11],gmsl.noaa$inthigh[11]), c(-10,10), lty=2, lwd=2, col=mycol.rgb[9] )
lines(c(gmsl.noaa$high[11],gmsl.noaa$high[11]),       c(-1,1.1), lty=2, lwd=2, col=mycol.rgb[2] )
lines(c(gmsl.noaa$extreme[11],gmsl.noaa$extreme[11]), c(-1,1.1), lty=2, lwd=2, col=mycol.rgb[6] )

legend(1.63,5.2,c("RCP2.6","RCP4.5","RCP8.5",
                  "NOAA, Low","NOAA, Int-Low", "NOAA, Int", "NOAA, Int-High", "NOAA, High","NOAA, Extreme"),
       lty=c(1,1,1,2,2,2,2,2,2), lwd=2,
       col=c(rgb(col26[1],col26[2],col26[3]),rgb(col45[1],col45[2],col45[3]),rgb(col85[1],col85[2],col85[3]),
             mycol.rgb[11], mycol.rgb[13], mycol.rgb[15], mycol.rgb[9], mycol.rgb[2], mycol.rgb[6]))

plot(ecdf.rcp85,log10(esf.vals),type='l', xlim=c(0,3), ylim=c(-3.45,0), lty=1,
     col=rgb(col85[1],col85[2],col85[3]), lwd=2, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i');
axis(1,seq(0,3,0.5),lab=c("0","0.5","1","1.5","2","2.5","3"))
mtext('Survival function [1-CDF]', side = 2, line=2.6);
mtext('Projected local mean sea level in 2100\nrelative to 1991-2009 average [m]', side=1, line=3);
mtext(side=3, text=expression(bold('   b')), line=-1.1, cex=.9, adj=0);

axis(2, at=seq(-4,0), label=parse(text=paste("10^", seq(-4,0), sep="")), las=1)

lines(ecdf.rcp45,log10(esf.vals),type='l',lty=1, col=rgb(col45[1],col45[2],col45[3]), lwd=2);
lines(ecdf.rcp26,log10(esf.vals),type='l',lty=1, col=rgb(col26[1],col26[2],col26[3]), lwd=2);

lines(c(-4,4),c(-2,-2),lty=2,col='black'); text(0.35,-1.85,"1:100 level");
lines(c(-4,4),c(-3,-3),lty=2,col='black'); text(0.35,-2.85,"1:1000 level");

# add NOAA scenarios as shaded regions/bars
lines(c(gmsl.noaa$low[11],gmsl.noaa$low[11]),         c(-1.6,10), lty=2, lwd=2, col=mycol.rgb[11])
lines(c(gmsl.noaa$intlow[11],gmsl.noaa$intlow[11]),   c(-1.6,10), lty=2, lwd=2, col=mycol.rgb[13])
lines(c(gmsl.noaa$int[11],gmsl.noaa$int[11]),         c(-10,10), lty=2, lwd=2, col=mycol.rgb[15])
lines(c(gmsl.noaa$inthigh[11],gmsl.noaa$inthigh[11]), c(-10,10), lty=2, lwd=2, col=mycol.rgb[9] )
lines(c(gmsl.noaa$high[11],gmsl.noaa$high[11]),       c(-10,10), lty=2, lwd=2, col=mycol.rgb[2] )
lines(c(gmsl.noaa$extreme[11],gmsl.noaa$extreme[11]), c(-10,10), lty=2, lwd=2, col=mycol.rgb[6] )

dev.off()

##==============================================================================


##==============================================================================
## Plot -- projections

## Get quantiles to plot
iproj <- which(t.proj==2000):which(t.proj==2100)
iend <- iproj[length(iproj)]
lsl_proj$q05$rcp26 <- rep(0,length(t.proj))
lsl_proj$q50$rcp26 <- rep(0,length(t.proj))
lsl_proj$q95$rcp26 <- rep(0,length(t.proj))
lsl_proj$q05$rcp45 <- rep(0,length(t.proj))
lsl_proj$q50$rcp45 <- rep(0,length(t.proj))
lsl_proj$q95$rcp45 <- rep(0,length(t.proj))
lsl_proj$q05$rcp85 <- rep(0,length(t.proj))
lsl_proj$q50$rcp85 <- rep(0,length(t.proj))
lsl_proj$q95$rcp85 <- rep(0,length(t.proj))
lsl_proj$max <- rep(0,length(t.proj))
lsl_proj$min <- rep(0,length(t.proj))
for (t in 1:length(t.proj)) {
    lsl_proj$q50$rcp26[t] <- quantile(lsl_proj$rcp26[t,], 0.5)
    lsl_proj$q05$rcp26[t] <- quantile(lsl_proj$rcp26[t,], .05)
    lsl_proj$q95$rcp26[t] <- quantile(lsl_proj$rcp26[t,], .95)
    lsl_proj$q50$rcp45[t] <- quantile(lsl_proj$rcp45[t,], 0.5)
    lsl_proj$q05$rcp45[t] <- quantile(lsl_proj$rcp45[t,], .05)
    lsl_proj$q95$rcp45[t] <- quantile(lsl_proj$rcp45[t,], .95)
    lsl_proj$q50$rcp85[t] <- quantile(lsl_proj$rcp85[t,], 0.5)
    lsl_proj$q05$rcp85[t] <- quantile(lsl_proj$rcp85[t,], .05)
    lsl_proj$q95$rcp85[t] <- quantile(lsl_proj$rcp85[t,], .95)
    lsl_proj$max[t]       <- max(lsl_proj$rcp85[t,])
    lsl_proj$min[t]       <- min(lsl_proj$rcp26[t,])
}

pdf(paste(plotdir,'projections_SLR2100_Norfolk.pdf',sep=''),width=4,height=6,colormodel='cmyk')

par(mfrow=c(2,1), mai=c(.5,.65,.20,.3))
plot(t.proj[iproj],lsl_proj$q50$rcp85[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
     xlim=c(2000,2100), ylim=c(0,3), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
axis(1, seq(2000,2100,by=20));
axis(2, seq(0,3,by=.5), lab=c('0','','1','','2','','3'), las=1);
mtext(side=1, text='Year', line=2.2, cex=1);
mtext(side=2, text='Local mean sea level [m]', line=2.2, cex=1);
#mtext(side=3, text=expression(bold(' a')), line=-1, cex=1, adj=0);
polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(lsl_proj$q95$rcp85[iproj],rev(lsl_proj$q05$rcp85[iproj])),
        col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
lines(t.proj[iproj],lsl_proj$q50$rcp45[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2.5);
polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(lsl_proj$q95$rcp45[iproj],rev(lsl_proj$q05$rcp45[iproj])),
        col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
lines(t.proj[iproj],lsl_proj$q50$rcp26[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(lsl_proj$q95$rcp26[iproj],rev(lsl_proj$q05$rcp26[iproj])),
        col=rgb(col26[1],col26[2],col26[3],.5), border=NA);
# max/min of our ensemble
polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(lsl_proj$max[iproj],rev(lsl_proj$min[iproj])),
        col=rgb(.5,.5,.5,.23), border=NA);

# superimpose the NOAA scenarios
lines(time.noaa, gmsl.noaa$low,     lty=2, lwd=2, col=mycol.rgb[6] )
lines(time.noaa, gmsl.noaa$intlow,  lty=2, lwd=2, col=mycol.rgb[2] )
lines(time.noaa, gmsl.noaa$int,     lty=2, lwd=2, col=mycol.rgb[9] )
lines(time.noaa, gmsl.noaa$inthigh, lty=2, lwd=2, col=mycol.rgb[15])
lines(time.noaa, gmsl.noaa$high,    lty=2, lwd=2, col=mycol.rgb[13])
lines(time.noaa, gmsl.noaa$extreme, lty=2, lwd=2, col=mycol.rgb[11])

# legend in second panel
par(mai=c(.05,.05,.05,.05))
# stationary
plot( 9, 9, xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE, xlim=c(0,1),ylim=c(0,1), xlab='', ylab='')

legend(0.15, 0.75,
    c("RCP2.6, 5-95% range", "RCP4.5, 5-95% range", "RCP8.5, 5-95% range", "Ensemble max/min range",
        "NOAA, Low", "NOAA, Intermediate-Low", "NOAA, Intermediate", "NOAA, Intermediate-High", "NOAA, High", "NOAA, Extreme"),
    lty=c(1,1,1,1,2,2,2,2,2,2), lwd=c(9,9,9,9,3,3,3,3,3,3), , bty='n', cex=1,
    col=c(rgb(col26[1],col26[2],col26[3]),rgb(col45[1],col45[2],col45[3]),rgb(col85[1],col85[2],col85[3]),rgb(.5,.5,.5),
        mycol.rgb[6],mycol.rgb[2],mycol.rgb[9],mycol.rgb[15],mycol.rgb[13],mycol.rgb[11]))

dev.off()
##==============================================================================


plot(t.proj, apply(lsl.scen[[6]], 1, mean), type='l', lwd=2,
     xlim=c(1960, 2060), ylim=c(-.4,1.2), col=mycol.rgb[11],
     xlab='Year', ylab='Regional sea level [meters]')
lines(t.proj, apply(lsl.scen[[5]], 1, mean), lwd=2, col=mycol.rgb[13])
lines(t.proj, apply(lsl.scen[[4]], 1, mean), lwd=2, col=mycol.rgb[15])
lines(t.proj, apply(lsl.scen[[3]], 1, mean), lwd=2, col=mycol.rgb[9])
lines(t.proj, apply(lsl.scen[[2]], 1, mean), lwd=2, col=mycol.rgb[2])
lines(t.proj, apply(lsl.scen[[1]], 1, mean), lwd=2, col=mycol.rgb[6])




##==============================================================================
# get an ensemble of projections of GEV parameters

# routine for projecting storm surge (Grinsted et al 2013)
source('../R/BRICK_project_surge.R')

gev_proj <- brick_surge_grinsted(temperature = temp,
                                 time.proj   = t.proj,
                                 tidegauge   = data,
                                 l.nonstat   = l.nonstat,
                                 niter.gev   = niter.gev,
                                 burnin      = burnin)
##==============================================================================


##==============================================================================
## Write a netCDF ensemble output file including each of the RCP scenarios:
## (1) global total sea level, (2) local sea level. Also will want each
## contribution to global sea level rise, for the hindcast plots

dim.tproj <- ncdim_def('time_proj', 'years', as.double(t.proj))
dim.lat <- ncdim_def('lat', 'deg N', as.double(length(lat.proj)))
dim.lon <- ncdim_def('lon', 'deg E', as.double(length(lon.proj)))

dim.ensemble <- ncdim_def('ens', 'ensemble member', as.double(1:ncol(lsl_proj$rcp26)), unlim=TRUE)

lat.out <- ncvar_def('lat.lsl', 'deg N', list(dim.lat), -999,
                  longname = 'latitude of local sea level point')
lon.out <- ncvar_def('lon.lsl', 'deg N', list(dim.lon), -999,
                  longname = 'longitude of local sea level point')

lsl26 <- ncvar_def('LocalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level (RCP26)')
gsl26 <- ncvar_def('GlobalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level (RCP26)')
temp26 <- ncvar_def('temp_RCP26', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP26)')
loc26 <- ncvar_def('gev_location_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'location parameter (mu) of storm surge GEV distribution (RCP26)')
sha26 <- ncvar_def('gev_shape_RCP26', '', list(dim.tproj, dim.ensemble), -999,
                  longname = 'shape parameter (xi) of storm surge GEV distribution (RCP26)')
sca26 <- ncvar_def('gev_scale_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'scale parameter (sigma) of storm surge GEV distribution (RCP26)')

lsl45 <- ncvar_def('LocalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level (RCP45)')
gsl45 <- ncvar_def('GlobalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level (RCP45)')
temp45 <- ncvar_def('temp_RCP45', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP45)')
loc45 <- ncvar_def('gev_location_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'location parameter (mu) of storm surge GEV distribution (RCP45)')
sha45 <- ncvar_def('gev_shape_RCP45', '', list(dim.tproj, dim.ensemble), -999,
                  longname = 'shape parameter (xi) of storm surge GEV distribution (RCP45)')
sca45 <- ncvar_def('gev_scale_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'scale parameter (sigma) of storm surge GEV distribution (RCP45)')

lsl85 <- ncvar_def('LocalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level (RCP85)')
gsl85 <- ncvar_def('GlobalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level (RCP85)')
temp85 <- ncvar_def('temp_RCP85', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP85)')
loc85 <- ncvar_def('gev_location_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'location parameter (mu) of storm surge GEV distribution (RCP85)')
sha85 <- ncvar_def('gev_shape_RCP85', '', list(dim.tproj, dim.ensemble), -999,
                  longname = 'shape parameter (xi) of storm surge GEV distribution (RCP85)')
sca85 <- ncvar_def('gev_scale_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'scale parameter (sigma) of storm surge GEV distribution (RCP85)')

outnc <- nc_create(filename.projout,
                list( gsl26, gsl45, gsl85, lsl26, lsl45, lsl85, temp26, temp45, temp85,
                loc26, loc45, loc85, sha26, sha45, sha85, sca26, sca45, sca85, lat.out, lon.out),
                force_v4 = TRUE)

ncvar_put(outnc, lat.out, lat.proj)
ncvar_put(outnc, lon.out, lon.proj)

ncvar_put(outnc, lsl26, lsl_proj$rcp26)
ncvar_put(outnc, gsl26, gmsl$rcp26)
ncvar_put(outnc, temp26, temperature$rcp26)
ncvar_put(outnc, loc26, gev_proj$gev.proj$location$rcp26)
ncvar_put(outnc, sha26, gev_proj$gev.proj$shape$rcp26)
ncvar_put(outnc, sca26, gev_proj$gev.proj$scale$rcp26)

ncvar_put(outnc, lsl45, lsl_proj$rcp45)
ncvar_put(outnc, gsl45, gmsl$rcp45)
ncvar_put(outnc, temp45, temperature$rcp45)
ncvar_put(outnc, loc45, gev_proj$gev.proj$location$rcp45)
ncvar_put(outnc, sha45, gev_proj$gev.proj$shape$rcp45)
ncvar_put(outnc, sca45, gev_proj$gev.proj$scale$rcp45)

ncvar_put(outnc, lsl85, lsl_proj$rcp85)
ncvar_put(outnc, gsl85, gmsl$rcp85)
ncvar_put(outnc, temp85, temperature$rcp85)
ncvar_put(outnc, loc85, gev_proj$gev.proj$location$rcp85)
ncvar_put(outnc, sha85, gev_proj$gev.proj$shape$rcp85)
ncvar_put(outnc, sca85, gev_proj$gev.proj$scale$rcp85)

nc_close(outnc)
##==============================================================================


##==============================================================================
## Some plots

## Note -- later (if desired) can remove this part and put into a separate routine

#install.packages('ncdf4')
library(ncdf4)

filename.projin <- "../output_model/BRICK_project-lsl-surge_Annapolis_20Feb2017.nc"

# Initialize lists for projections
scen.rcp <- c('rcp26','rcp45','rcp85')
temp.proj <- vector('list', length(scen.rcp)); names(temp.proj) <- scen.rcp
lsl.proj <- temp.proj
gmsl.proj <- temp.proj
gev.proj <- vector('list',3); names(gev.proj) <- c('location','shape','scale')
gev.proj$location <- temp.proj
gev.proj$shape <- temp.proj
gev.proj$scale <- temp.proj

# Read sea level output
ncdata <- nc_open(filename.projin)
  for (rcp in scen.rcp) {
    if(rcp=='rcp26') {
        temp.proj[[rcp]]         <- ncvar_get(ncdata, 'temp_RCP26')
        gmsl.proj[[rcp]]         <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
        lsl.proj[[rcp]]          <- ncvar_get(ncdata, 'LocalSeaLevel_RCP26')
        gev.proj$location[[rcp]] <- ncvar_get(ncdata, 'gev_location_RCP26')
        gev.proj$shape[[rcp]]    <- ncvar_get(ncdata, 'gev_shape_RCP26')
        gev.proj$scale[[rcp]]    <- ncvar_get(ncdata, 'gev_scale_RCP26')
    } else if(rcp=='rcp45') {
        temp.proj[[rcp]]         <- ncvar_get(ncdata, 'temp_RCP45')
        gmsl.proj[[rcp]]         <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
        lsl.proj[[rcp]]          <- ncvar_get(ncdata, 'LocalSeaLevel_RCP45')
        gev.proj$location[[rcp]] <- ncvar_get(ncdata, 'gev_location_RCP45')
        gev.proj$shape[[rcp]]    <- ncvar_get(ncdata, 'gev_shape_RCP45')
        gev.proj$scale[[rcp]]    <- ncvar_get(ncdata, 'gev_scale_RCP45')
    } else if(rcp=='rcp60') {
        temp.proj[[rcp]]         <- ncvar_get(ncdata, 'temp_RCP60')
        gmsl.proj[[rcp]]         <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP60')
        lsl.proj[[rcp]]          <- ncvar_get(ncdata, 'LocalSeaLevel_RCP60')
        gev.proj$location[[rcp]] <- ncvar_get(ncdata, 'gev_location_RCP60')
        gev.proj$shape[[rcp]]    <- ncvar_get(ncdata, 'gev_shape_RCP60')
        gev.proj$scale[[rcp]]    <- ncvar_get(ncdata, 'gev_scale_RCP60')
    } else if(rcp=='rcp85') {
        temp.proj[[rcp]]         <- ncvar_get(ncdata, 'temp_RCP85')
        gmsl.proj[[rcp]]         <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
        lsl.proj[[rcp]]          <- ncvar_get(ncdata, 'LocalSeaLevel_RCP85')
        gev.proj$location[[rcp]] <- ncvar_get(ncdata, 'gev_location_RCP85')
        gev.proj$shape[[rcp]]    <- ncvar_get(ncdata, 'gev_shape_RCP85')
        gev.proj$scale[[rcp]]    <- ncvar_get(ncdata, 'gev_scale_RCP85')
    } else {
        print(paste('Error - unrecognized RCP scenario: ',rcp,sep=''))
    }
  }
  t.proj     <- ncvar_get(ncdata, 'time_proj')
  n.ensemble <- ncdata$dim$ens$len
nc_close(ncdata)

## Get statistical libraries needed for GEV calculations
library(extRemes)
library(fExtremes)
library(ismev)
library(lubridate)
library(zoo)
library(Bolstad)        # integrate (normalize) via Simpson's rule

## Calculate pdfs and sfs for sea-level rise and storm surge
iproj <- which(t.proj==2100)
inorm <- which(t.proj==1985):which(t.proj==2005)

## Initialize some list arrays
init <- vector('list',length(scen.rcp)); names(init) <- scen.rcp
lsl.norm <- init
f.lsl <- init
f.surge <- init
cdf.lsl <- init
cdf.surge <- init
sf.lsl <- init
sf.surge <- init

## distribution of sea-level, subsidence, and storm surges
lsl.lower <- 0
lsl.upper <- 20
lsl.n <- 2^11
lsl.x <- seq(lsl.lower, lsl.upper, length.out=lsl.n)
lsl.dx <- median(diff(lsl.x))

for (rcp in scen.rcp) {
    lsl.norm[[rcp]] <- lsl.proj[[rcp]] - t(replicate(nrow(lsl.proj[[rcp]]),apply(lsl.proj[[rcp]][inorm,],2,mean)))
    f.lsl[[rcp]]    <- density(x=lsl.norm[[rcp]][iproj,], from=lsl.lower, to=lsl.upper, n=lsl.n, kernel='gaussian')
    f.lsl[[rcp]]    <- f.lsl[[rcp]]$y/sintegral(x=f.lsl[[rcp]]$x, fx=f.lsl[[rcp]]$y)$value
    cdf.lsl[[rcp]]  <- cumsum(f.lsl[[rcp]]*lsl.dx)
    cdf.lsl[[rcp]]  <- cdf.lsl[[rcp]] - cdf.lsl[[rcp]][1]     # normalize (begin at 0)
    cdf.lsl[[rcp]]  <- cdf.lsl[[rcp]]/cdf.lsl[[rcp]][lsl.n]   # normalize (end at 1)
    sf.lsl[[rcp]]   <- 1-cdf.lsl[[rcp]]

## TODO
## TODO -- herenow -- add distributions for storm surge
## TODO

  }
}


##=====================================
## Figure 1 -- pd and survival functions of sea-level and storm surge in


##=====================================


##==============================================================================


##==============================================================================
## End
##==============================================================================
