##==============================================================================
## Script to calibrate Antarctic Ocean temperature emulator for use within
## BRICK coupled modeling framework.
##
## Parameters to calibrate, and standard values (Shaffer, 2014):
## [1] a.anto
## [2] b.anto
##
## Questions? Tony Wong <twong@psu.edu>
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

#rm(list =ls()) #Clear global environment

## Setup packages and libraries
#install.packages('compiler')
#install.packages('pscl')
#install.packages('sn')
#install.packages('fMultivar')
library(sn)
library(fMultivar)
library(pscl) # install inverse gamma distribution
library(compiler)
enableJIT(3)

## Read the data forcing for hindcasts and projections. Yields:
##  Ta (Antarctic temperature reduced to sea-level)
##  Toc (High latitude subsurface ocean temperature)
##  SL, obs.sl (Reconstructed sea level, Church and White (2011) modern-era sea level)
##  dSL (time rate of change of sea level)
##  date (240 Kyr before present to 2100 AD at 1-year intervals of forcings)
## Note: already have from previous call in BRICK_calib_driver.R.
#source('../calibration/DAIS_readData.R')

## To run BRICK in a full-forward capacity, want to calibrate so that Tg relative
## to 1850 can give us Toc relative to Toc(1961-1990 mean)=0.72
date.ce <- date+2000
i1961 <- which(date.ce==1961):which(date.ce==1990)
# verify that it is normalize to 0.72 in 1961-1990 mean
print(paste('should be about 0.72: ',mean(Toc[i1961]),sep=''))
Toc.norm <- Toc[which(date.ce==Tg.time[1]):which(date.ce==Tg.time[length(Tg.time)])]
# normalize Tg relative to 1850, since that's where temperature forcing can start
Tg.norm <- Tg - Tg[1]

## Source the DAIS/ANTO model
source('../fortran/R/daisantoF.R')

## Set up the parameters to calibrate
parnames.anto    = c('anto.a','anto.b')
parameters0.anto = c(  0.1574, 0.6677 )
bound.lower.anto = c( 0.0    , 0      )
bound.upper.anto = c( 1.0    , 2      )
##==============================================================================
rmse = function(parameters.in,
                parnames.in,
                Toc.obs.in,
                Tg.in
) {
    a.anto <- parameters.in[match("anto.a",parnames.in)]
    b.anto <- parameters.in[match("anto.b",parnames.in)]
    mod <- anto(a=a.anto, b=b.anto, Tf=-1.8, Tg=Tg.in)
    rmse <- sqrt( mean( (mod-Toc.obs.in)^2 ))

    rmse
}
##==============================================================================

#Tg.recon.norm <- Tg.recon - mean(Tg.recon[ind.relative])
Tg.recon.norm <- Tg.recon - Tg.recon[which(date.ce==1850)]
Toc.recon.norm <- Toc # no subtraction here - want with 1961-1990 mean = 0.72degC
Toc.anto.norm <- anto(a=.25, b=.05, Tf=-1.8, Tg=Tg.recon.norm)
Toc.anto.norm <- Toc.anto.norm - (mean(Toc.anto.norm[ind.relative])-0.72)
plot(Tg.recon.norm, type='l')
lines(Toc.recon.norm, type='l', col='blue')
lines(Toc.anto.norm, type='l', col='red')

##==============================================================================
## Preliminary Latin Hypercube to find decent parameter ranges/values

require(lhs)

t0=proc.time()
# Draw LHS sample
n.lhs = 10000
parameters.lhs <- randomLHS(n.lhs, length(parnames.anto))

# Transform unif(0,1) to the parameter bounds
for (i in 1:length(parnames.anto)) {
  parameters.lhs[,i] <- qunif(parameters.lhs[,i], bound.lower.anto[i], bound.upper.anto[i])
}
colnames(parameters.lhs)=parnames.anto
t1=proc.time()

rmse.lhs = rep(NA,n.lhs)

pb <- txtProgressBar(min=0,max=n.lhs,initial=0,style=3)
for (j in 1:n.lhs) {
  rmse.lhs[j] = rmse(   parameters.in = as.numeric(parameters.lhs[j,]),
                        parnames.in = parnames.anto,
                        Toc.obs.in = Toc.recon.norm,
                        Tg.in = Tg.recon.norm
                        )
  setTxtProgressBar(pb, j)
}
close(pb)
t2=proc.time()

tmp <- data.frame( cbind(rmse.lhs,parameters.lhs))
tmp.sort <- tmp[order(rmse.lhs),]

# pick a nice cut-off for "low enough" RMSE
quantile(tmp.sort[,1],c(.01,.05,.1))

par(mfrow=c(2,2))
# these indicate that restricting RMSE < 0.5 provides a nice constraint on
# anto.a and anto.b
itmp <- which(tmp.sort[,1]<0.5)
plot(tmp.sort[,2],tmp.sort[,1],xlab='anto.a',ylab='RMSE')
plot(tmp.sort[,3],tmp.sort[,1],xlab='anto.b',ylab='RMSE')
plot(tmp.sort[itmp,2],tmp.sort[itmp,3],xlab='anto.a',ylab='anto.b')

# get the min/max range of these low-RMSE parameter sets
print(paste('low-RMSE anto.a range=',quantile(tmp.sort[itmp,2],0),' to ',quantile(tmp.sort[itmp,2],1)))
print(paste('low-RMSE anto.b range=',quantile(tmp.sort[itmp,3],0),' to ',quantile(tmp.sort[itmp,3],1)))
print(paste('min-RMSE (anto.a, anto.b)=(',tmp.sort[1,2],',',tmp.sort[1,3],')'))

# joint prior for anto.a and anto.b?

# first, run precalibration LHS in scratch_findANTObounds_precalibration.R
# will yield tmp.sort, where columns 2 and 3 are anto.a and anto.b, respectively


anto.prior <- msnFit(tmp.sort[itmp,2:3])
anto.prior.fit <- anto.prior@fit$estimated
# account for truncated prior
itmp <- c(match('anto.a',parnames),match('anto.b',parnames))
anto.prior.fit$cnorm <- as.numeric(pmsn( bound.upper[itmp], anto.prior.fit$beta, anto.prior.fit$Omega, anto.prior.fit$alpha) -
                         pmsn( bound.lower[itmp], anto.prior.fit$beta, anto.prior.fit$Omega, anto.prior.fit$alpha))

##==============================================================================
## End
##==============================================================================
