##==============================================================================
##	DAIS_priors.R
##
## Define DAIS parameter prior ranges, based on BRICK v0.1 experiments with
## calibrating standalone DAIS model using paleo data.
##
## Question? Tony Wong <twong@psu.edu>
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

#install.packages('fitdistrplus')
library(fitdistrplus)
library(ncdf4)

# read calibrated BRICK parameters file
filename.parameters <- '../output_calibration/BRICK-model_postcalibratedParameters_control_01Nov2016.nc'
ncdata <- nc_open(filename.parameters)
parameters.fit = ncvar_get(ncdata, 'BRICK_parameters')
parnames.fit = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters.fit = t(parameters.fit)
colnames(parameters.fit) = parnames.fit

# function to scale a random variable x on [a,b] to [0,1]
range.to.beta <- function(x, a, b){
    return((x-a)/(b-a))
}
# function to scale a random variable y on [0,1] to [a,b]
beta.to.range <- function(y, a, b){
    return(y*(b-a)+a)
}

# fit beta distribution to each of the DAIS parameters, except var.dais and ANTO
parnames.tmp <- parnames.dais[3:(length(parnames.dais)-1)]
dais.priors <- mat.or.vec(2,length(parnames.tmp))
colnames(dais.priors) <- parnames.tmp

for (p in 1:length(parnames.tmp)) {
    p.fit <- match(parnames.tmp[p],parnames.fit)
    ptmp <- range.to.beta(parameters.fit[,p.fit],
                          bound.lower[match(parnames.tmp[p],parnames)],
                          bound.upper[match(parnames.tmp[p],parnames)])
    dais.priors[,p] <- fitdist(ptmp, 'beta')$estimate
}

# overwrite a couple highly correlated one with a 2D skew-normal prior fit

# install.packages('sn')
# install.packages('fMultivar')
library(sn)
library(fMultivar)

# gamma-slope-P0
itmp <- c(match('gamma',parnames.fit), match('slope',parnames.fit), match('P0',parnames.fit))
gamma_slope_P0.prior <- msnFit(parameters.fit[,itmp])
gamma_slope_P0.prior.fit <- gamma_slope_P0.prior@fit$estimated
dais.priors <- dais.priors[,-match('gamma',colnames(dais.priors))]
dais.priors <- dais.priors[,-match('slope',colnames(dais.priors))]
dais.priors <- dais.priors[,-match('P0',colnames(dais.priors))]

# c-h0
itmp <- c(match('c',parnames.fit), match('h0',parnames.fit))
c_h0.prior <- msnFit(parameters.fit[,itmp])
c_h0.prior.fit <- c_h0.prior@fit$estimated
dais.priors <- dais.priors[,-match('h0',colnames(dais.priors))]
dais.priors <- dais.priors[,-match('c',colnames(dais.priors))]


if(FALSE){  # if you want to visualize these priors
X <- parameters.fit[,itmp]
ans <- c_h0.prior                                           # <<< for example
plot(hexBinning(X[,1], X[, 2], bins = 30), main="Skew Normal")
N <- 101
x <- seq(min(X[, 1]), max(X[, 1]), l=N)
y <- seq(min(X[, 2]), max(X[, 2]), l=N)
u <- grid2d(x, y)$x
v <- grid2d(x, y)$y
XY <- cbind(u, v)
param <- ans@fit$estimate
Z <- matrix(dmsn(XY, param[[1]][1,], param[[2]], param[[3]]), ncol=N)
contour(x, y, Z, add=TRUE, col="green", lwd=2)
grid(col="brown", lty=3)
}

##==============================================================================
## End
##==============================================================================
