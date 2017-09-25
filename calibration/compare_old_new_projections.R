##==============================================================================
## Compare a set of old projections (or perhaps one set of model structures)
## with a set of new projections (or an alternative formulation).
##
## Required:
##   - two file names, containing BRICK projections to compare
##     these are assumed to be the results from the processingPipeline script.
##   - project.years, the years in which you would like pdfs and survival
##     functions of global mean sea level
##   - ref.period, the beginning and end of the reference period, to which the
##     global mean sea level projections will be relative
##
## Yields:
##   - some numbers (quantiles) of projections at user-defined years relative to
##     a user-defined reference period
##   - pdfs and survival functions for each of the elements of project.years
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================
## Copyright 2017 Tony Wong, Alexander Bakker
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

rm(list=ls())

## navigate to your local "BRICK/output_model" directory

setwd('~/codes/BRICK/output_model')

## define files that result from the rejection sampling calibration to GMSL
## (these file names are relative to the "BRICK/output_model" directory)

filename.old <- 'BRICK-model_physical_fd-gamma_13Aug2017.nc'
filename.new <- 'BRICK-model_physical_fd-gamma_17Aug2017.nc'

## define reference period (beginning and end years)
## if this is NULL, projections will be relative to 1986-2005 global mean sea
## level

ref.period <- c(1986,2005)

## define *single* years in which you would like the sea-level projections
## (have the option later of asking for projections as mean over a period)
## If this is NULL, projections will be returned for 2100

project.years <- c(2065,2100)

## which RCP scenarios do you want? subset of {rcp26, rcp45, rcp85}
scen.rcp <- c('rcp26','rcp45','rcp85')

## which quantiles do you want?
## If this is NULL, default will be min, 1%, 5%, 50%, 95% 99% and max.
quantiles.to.get <- c(0, .01, .05, .5, .95, .99, 1)

##==============================================================================
## end user-defined stuff
##==============================================================================

## required libraries
library(ncdf4)

## And set the IPCC RCP colors
col26 <- c(0, 0, 255)/255
col45 <- c(121, 188, 255)/255
col60 <- c(255, 130, 45)/255
col85 <- c(255, 0, 0)/255

## check the inputs, replace with defaults if needed
if(is.null(ref.period)) {ref.period <- c(1986,2005)}
if(is.null(project.years)) {project.years <- 2100}

##==============================================================================

## read model results

gmsl.new        <- gmsl.old        <- vector('list', length(scen.rcp))
names(gmsl.new) <- names(gmsl.old) <- scen.rcp

ncdata <- nc_open(filename.old)
  tproj.old <- ncvar_get(ncdata, 'time_proj')
  if(is.element('rcp26', scen.rcp)) {gmsl.old$rcp26 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')}
  if(is.element('rcp45', scen.rcp)) {gmsl.old$rcp45 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')}
  if(is.element('rcp85', scen.rcp)) {gmsl.old$rcp85 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')}
nc_close(ncdata)

ncdata <- nc_open(filename.new)
  tproj.new <- ncvar_get(ncdata, 'time_proj')
  if(is.element('rcp26', scen.rcp)) {gmsl.new$rcp26 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')}
  if(is.element('rcp45', scen.rcp)) {gmsl.new$rcp45 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')}
  if(is.element('rcp85', scen.rcp)) {gmsl.new$rcp85 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')}
nc_close(ncdata)

## throw an error if the reference period is not contained within one set of
## projections or the other

if(!(min(tproj.new) <= ref.period[1] & max(tproj.new) >= ref.period[2] &
     min(tproj.old) <= ref.period[1] & max(tproj.new) >= ref.period[2])) {
  print('ERROR: reference period not contained within one or both sets of projections')
}

## check that tproj is the first (rows) dimension of gmsl$rcpXX

for (rcp in scen.rcp) {
  # check gmsl.old
  dims <- dim(gmsl.old[[rcp]])
  if (dims[1] != length(tproj.old)) {
    if (dims[2] == length(tproj.old)) {
      # dimensions are just swapped, so transpose
      gmsl.old[[rcp]] <- t(gmsl.old[[rcp]])
    } else {
      # neither dimension has the proper length, so throw an error
      print('ERROR: neither dimension of gmsl.old matches the tproj length')
    }
  }
  # check gmsl.new
  dims <- dim(gmsl.new[[rcp]])
  if (dims[1] != length(tproj.new)) {
    if (dims[2] == length(tproj.new)) {
      # dimensions are just swapped, so transpose
      gmsl.new[[rcp]] <- t(gmsl.new[[rcp]])
    } else {
      # neither dimension has the proper length, so throw an error
      print('ERROR: neither dimension of gmsl.new matches the tproj length')
    }
  }
}

## trim projections so they cover the same time period, and get a common tproj

tproj.beg <- min(c(tproj.old, tproj.new))
tproj.end <- max(c(tproj.old, tproj.new))
tproj <- seq(from=tproj.beg, to=tproj.end, by=1)

## trim the GMSL projections
for (rcp in scen.rcp) {
  # trim the gmsl.old arrays
  tbeg <- which(tproj.old==tproj[1])
  tend <- which(tproj.old==max(tproj))
  gmsl.old[[rcp]] <- gmsl.old[[rcp]][tbeg:tend,]
  # trim the gmsl.new arrays
  tbeg <- which(tproj.new==tproj[1])
  tend <- which(tproj.new==max(tproj))
  gmsl.new[[rcp]] <- gmsl.new[[rcp]][tbeg:tend,]
}

## make sure both sets are normalized to the reference period, ref.period

# what are the indices within tproj?
ind.norm <- which(tproj==ref.period[1]):which(tproj==ref.period[2])

for (rcp in scen.rcp) {
  gmsl.old[[rcp]] <- sapply(1:ncol(gmsl.old[[rcp]]), function(i) {gmsl.old[[rcp]][,i] - mean(gmsl.old[[rcp]][ind.norm, i])})
  gmsl.new[[rcp]] <- sapply(1:ncol(gmsl.new[[rcp]]), function(i) {gmsl.new[[rcp]][,i] - mean(gmsl.new[[rcp]][ind.norm, i])})
}

##==============================================================================
## create a table of quantiles
##==============================================================================

# create names for the years in which you want projections
# start with a 'y' (for 'year') because R will not let you reference list
# objects using numbers
project.years.names <- rep(NA, length(project.years))
for (ind.year in 1:length(project.years)) {
  project.years.names[ind.year] <- paste('y',project.years[ind.year], sep='')
}
names(project.years) <- project.years.names

# initialize lists to hold the quantiles for each of the projection years
quantiles.new        <- quantiles.old        <- vector('list', length(project.years.names))
names(quantiles.new) <- names(quantiles.old) <- project.years.names

# and for each of the projection years, need each of the RCP scenarios
for (year in project.years.names) {
  quantiles.new[[year]]           <- quantiles.old[[year]]           <- mat.or.vec(nr=length(scen.rcp), nc=length(quantiles.to.get))
  rownames(quantiles.new[[year]]) <- rownames(quantiles.old[[year]]) <- scen.rcp
  colnames(quantiles.new[[year]]) <- colnames(quantiles.old[[year]]) <- quantiles.to.get
  for (rcp in scen.rcp) {
    iproj <- which(tproj==project.years[[year]])
    quantiles.old[[year]][rcp,] <- quantile(gmsl.old[[rcp]][iproj,], quantiles.to.get)
    quantiles.new[[year]][rcp,] <- quantile(gmsl.new[[rcp]][iproj,], quantiles.to.get)
  }
}

# TODO
# TODO -- make this output into a CSV table
# TODO

##==============================================================================
## pdf and survival function figure
##==============================================================================

# fit kernel density estimates to each of the sets of projections
#tmp <- density(x=sf.sealev[[rcp]][[ais]], from=lsl.lower, to=lsl.upper, n=lsl.n, kernel=kern)

# TODO
# TODO!
# TODO

# calculate empirical survival function (by first calculating the empirical cdf)
# for each set of projections

ecdf.vals.old <- seq(from=0, to=1, length.out=ncol(gmsl.old[[1]]))
ecdf.vals.new <- seq(from=0, to=1, length.out=ncol(gmsl.new[[1]]))
esf.vals.old <- 1-ecdf.vals.old
esf.vals.new <- 1-ecdf.vals.new

# get a sample that is normally-distributed with mean and standard deviation
# matching that of the 'new' and 'old' ensembles. then you can superimpose the
# normal samples and see how much heavier the tails are
normal.sample.old <- rnorm(n=length(esf.vals.old), mean=mean(gmsl.old$rcp85[251,]), sd=sd(gmsl.old$rcp85[251,]))
normal.sample.new <- rnorm(n=length(esf.vals.new), mean=mean(gmsl.new$rcp85[251,]), sd=sd(gmsl.new$rcp85[251,]))

plot(sort(gmsl.old$rcp85[251,]), log10(esf.vals.old), type='l', ylim=c(-5,0), lwd=2)
lines(sort(gmsl.new$rcp85[251,]), log10(esf.vals.new), col='red', lwd=2)
lines(sort(normal.sample.old), log10(esf.vals.old), col='black', lty=2, lwd=2)
lines(sort(normal.sample.new), log10(esf.vals.new), col='red', lty=2, lwd=2)

# TODO
# TODO finish this
# TODO

##==============================================================================
## End
##==============================================================================
