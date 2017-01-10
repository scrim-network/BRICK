##==============================================================================
## Script to read global average sea level data (Church and White (2011)) for
## precalibration of thermosteric expansion (TE) model in BRICK.
## Also grabs the annual global mean surface temperature anomaly forcing data
## (HadCRUT4 data, same as is used to calibrate DOECLIM).
##
## Returns:
##  obs.sl          Church and White (2011) observations of global mean sea level
##  obs.sl.time     times (years AD) of these sea level observations
##  obs.sl.err      uncertainties associated with these sea level data
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

## Get the sea level data; new version, Church and White updated to include through 2013
dat = read.table("../data/GMSL_ChurchWhite2011_yr_2015.txt")
obs.sl      = dat[,2]/1000              # data are in mm
obs.sl.time = dat[,1]-0.5 # they give half-year so it is unambiguous as to which year they mean
obs.sl.err  = dat[,3]/1000

## Useful: set the indices of 1961-1990 (everything is anomaly relative to this
## time period, so make sure the sea levels are too)
ibeg=which(obs.sl.time==begyear.norm)
iend=which(obs.sl.time==endyear.norm)
obs.sl = obs.sl - mean(obs.sl[ibeg:iend])

idx = compute_indices(obs.time=obs.sl.time, mod.time=mod.time)
oidx.sl = idx$oidx; midx.sl = idx$midx

## Set trends from IPCC AR5 Ch 13 (observational) to match:
## trends.te = [Column 1=trend ; Column 2-3=90% window ; Column 4-5=beginning/ending year ;
##              Column 6-7=beginning/ending model indices for trend period]
## Note: IPCC trends are through 2010, but for the hindcast calibration, the
## forcing data go through 2009.
trends.te = mat.or.vec( 2 , 7)

trends.te[1,1:5] = c( 0.8 , 0.5 , 1.1 , 1971 , 2009 )
trends.te[2,1:5] = c( 1.1 , 0.8 , 1.4 , 1993 , 2009 )

## Compute the indices for the model to compare trends
for (i in 1:nrow(trends.te)) {
  idx = compute_indices(obs.time=(trends.te[i,4]:trends.te[i,5]) , mod.time=mod.time)
  trends.te[i,6:7] = c(idx$midx[1],(idx$midx[length(idx$midx)]))
}

##==============================================================================
## End
##==============================================================================
