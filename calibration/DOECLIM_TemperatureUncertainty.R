##==============================================================================
## Script to estimate the total ucnertainty in the HadCRUT4 temperature anomaly
## data.
## Read in all of the monthly files, which include these uncertainty estimates,
## and aggregate them into yearly uncertainty. The yearly files provided only
## have the ensemble 5-95% range, which is (likely) an overestimate of the
## actual uncertainty.
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

## Each of the HadCRUT.4.4.0.0.annual_ns_avg.XX.txt files is an individual
## ensemble member.
## Column 1 = date (year)
## Column 2 = median
## Column 3 = 1-sigma uncertainty (includes Columns 4, 5, and 6 uncertainties)
dat.dir = "../data/HadCRUT.4.4.0.0.annual_ns_avg_realisations/"
files = list.files(dat.dir)

## Number of ensemble members in HadCRUT4
n.ensemble = length(files)

## Read a single file to get the sizes of the data fields
dat = read.table(paste(dat.dir,files[1],sep=""))
n.year = nrow(dat)

## Initialize things
tg.medians = mat.or.vec( n.year , n.ensemble )
tg.errs = mat.or.vec( n.year , n.ensemble )

## Go through the files and grab the ensemble member medians and errors by year
for (i in 1:n.ensemble) {
  dat = read.table(paste(dat.dir,files[i],sep=""))
  tg.medians[,i] = dat[,2]
  tg.errs[,i] = dat[,3]
}

## Note: all of the ensemble members have the same uncertainties assigned to them.
## So it doesn't really matter whether you take mean/median/whatever of their
## individual uncertainties.
tg.err = tg.errs[,1]

## Check the ensemble median of medians against the global ensemble average file
tg.median = rep(NA,n.year)
tg.mean = rep(NA,n.year)
for (t in 1:n.year) {
  tg.median[t] = median(tg.medians[t,])
  tg.mean[t] = mean(tg.medians[t,])
}


##==============================================================================
## End
##==============================================================================
