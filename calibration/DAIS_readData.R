#############################################################
## -original file = DAIS_data.R
## -Antarctic Ice Sheet (AIS) model
#############################################################
## -Model can be found in Shaffer_GMDD_2014
## -Matlab codes and forcings can be found at:
## -www.dcess.dk under "DCESS models"
## -Sea-level values are reconstructed from:
## -Waelbroeck et al 2002(21-240kyr BP), Clark et al 2012(7000-21000BP),
## -Lambeck et al 2010(6000BP-1869AD), and Church & White 2011(1870-2010AD)
## -Future Sea level values and and rates are calculated using the Rahmstorf (2007) model, DEoptim, and RCP8.5
## -Future air and ocean temperatures were calculated by Robert Nicholas and Varada Vaidya using the CNRM-CMIP5 model
## -Author: Kelsey Ruckert (klr324@psu.edu)
##############################################################
## -June 17, 2014 #Updates June 10 2015
## - Modified 6 July 2016 by Tony Wong (twong@psu.edu) for use in the BRICK
## model framework
##############################################################
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

## read in forcing data

date = seq(-239999,300,1) #240 Kyr BP to 2100AD at one year intervals of the forcings
dSL = scan("../data/future_GSL.txt")   #Time rate of change of sea-level
Ta  = scan("../data/future_TA.txt")     #Antarctic temperature reduced to sea-level
Toc = scan("../data/future_TO.txt")     #High latitude subsurface ocean temperature
SL  = scan("../data/future_SL.txt")     #Reconstructed sea-level

## Trim to only to 2016
iend=which(date==16)
date=date[1:iend]
dSL=dSL[1:iend]
SL=SL[1:iend]
Ta=Ta[1:iend]
Toc=Toc[1:iend]

## Set the normalization period, and find indices for the DAIS runs
## Only set the normalization period here, otherwise you might screw things up
norm.period=c(1961,1990)
ibeg=which(date==(norm.period[1]-2000))
iend=which(date==(norm.period[2]-2000))
ind.relative=ibeg:iend

## Read Church and White (2011) updated sea level data for pre-calibration
sl.dat.new = read.table("../data/GMSL_ChurchWhite2011_yr_2015.txt")     #Reconstructed sea-level
SL.time= sl.dat.new[,1]-0.5     # times are in half-year
SL.new = sl.dat.new[,2]/1000    # data are in mm
SL.err = sl.dat.new[,3]/1000
ibeg=which(SL.time==norm.period[1]); iend=which(SL.time==norm.period[2]);
SL.new = SL.new - mean(SL.new[ibeg:iend])               # make sure SL data are relative to 1961-1990

## beginning and ending indices (within "date") of the Church and White sea level
## data, to impose an upper bound on what the AIS can contribute to SLR
midx.sl = which(date==(SL.time[1]-2000)):which(date==(SL.time[length(SL.time)]-2000))

##==============================================================================

## Linear regression between reconstructed Antarctic surface temperatures (Ta)
## and global surface temperature anomalies (obs.temp). This will yield a
## lengthy (240,000 years before present - present) time series of global average
## surface temperature anomalies, akin to how DAIS will be forced in the coupled
## BRICK model.

dat = read.table("../data/HadCRUT.4.4.0.0.annual_ns_avg.txt")
Tg = dat[,2]-mean(dat[1:20,2])
Tg.time = dat[,1]
i1997 = which(Tg.time==1997) # end fit at 1997 because Ta is obs after that
ibeg=which(date==-150); iend=ibeg+length(Tg.time[1:i1997])-1;
fit = lm(Tg[1:i1997] ~ Ta[ibeg:iend])
intercept.Ta2Tg = fit$coefficients[1]
slope.Ta2Tg = fit$coefficients[2]
Tg.recon = fit$coefficients[1] + fit$coefficients[2]*Ta
      # Tg.recon is now an appropriately scaled global surface temperature anomaly
      # akin to output from DOECLIM.

##==============================================================================
## End
##==============================================================================
