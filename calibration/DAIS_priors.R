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

#install.packages('sn')
#install.packages('fMultivar')
#install.packages('fitdistrplus')
#install.packages('ncdf4')
library(fitdistrplus)
library(ncdf4)
library(sn)
library(fMultivar)

# read calibrated BRICK parameters file
filename.parameters <- '../output_calibration/BRICK-model_postcalibratedParameters_control_01Nov2016.nc'
ncdata <- nc_open(filename.parameters)
parameters.fit <- ncvar_get(ncdata, 'BRICK_parameters')
parnames.fit <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters.fit <- t(parameters.fit)
colnames(parameters.fit) <- parnames.fit

# all of them? (not ANTO)
# 3rd on is gamma, length-1 is slope
parnames.tmp <- parnames.dais[3:(length(parnames.dais)-1)]

# if this is the chr version, get rid of h0 and c from the fit parameters and
# create a vector of the indices of the rights ones within parnames that
# do correspond to the fit
ichr <- match('chr',parnames.tmp)
if(!is.null(ichr)) {parnames.tmp <- parnames.tmp[-ichr]}
ind.dais.prior    <- rep(NA, length(parnames.tmp))
bounds.upper.dais <- rep(NA, length(parnames.tmp))
bounds.lower.dais <- rep(NA, length(parnames.tmp))
itmp <- NULL
for (p in 1:length(parnames.tmp)) {
    itmp <- c(itmp,match(parnames.tmp[p],parnames.fit))
    ind.dais.prior[p]    <- match(parnames.tmp[p],parnames)
    bounds.upper.dais[p] <- bound.upper[ind.dais.prior[p]]
    bounds.lower.dais[p] <- bound.lower[ind.dais.prior[p]]
}
dais.all.prior <- msnFit(parameters.fit[,itmp])
dais.prior.fit <- dais.all.prior@fit$estimated

# account for truncation; these are normal (skewed and multivariate, but still
# normally distributed), so farthest from origin has highest CDF and closest to
# origin has lowest.
dais.prior.fit$cnorm <- as.numeric(pmsn( bounds.upper.dais, dais.prior.fit$beta, dais.prior.fit$Omega, dais.prior.fit$alpha) -
                         pmsn( bounds.lower.dais, dais.prior.fit$beta, dais.prior.fit$Omega, dais.prior.fit$alpha))

# and store the indices too
dais.prior.fit$ind <- ind.dais.prior

##==============================================================================
## End
##==============================================================================
