##==============================================================================
##	This script is for testing running ensembles of fully-forward coupled
##	simulations by using the results for the 1850 "initial condition" as the
##	uncertain initial condition in:
##		DOECLIM	--	T0, H0
##		GSIC		--	Gs0 (may need to adjust V0 higher)
##		SIMPLE	--	V0
##		TE			--	TE0
##		DAIS		--	n/a (initial ice sheet volume calculated from other parameters)
##
##	Questions? -- Tony Wong <twong@psu.edu
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

rm(list=ls())

setwd('~/codes/BRICK/calibration')

library(ncdf4)

# read BRICK v0.1 calibrated model simulations
filename.brick <- '../output_model/BRICK-model_physical_control_01Nov2016.nc'
ncdata <- nc_open(filename.brick)
  t.hind <- ncvar_get(ncdata, 'time_hind')
  t.proj <- ncvar_get(ncdata, 'time_proj')
  gmsl.ctrl.hind <- ncvar_get(ncdata, 'GlobalSeaLevel_hind')
  temp.ctrl.hind <- ncvar_get(ncdata, 'temp_hind')
  ocheat.ctrl.hind <- ncvar_get(ncdata, 'ocheat_hind')
  ais.ctrl.hind <- ncvar_get(ncdata, 'AIS_hind')
  gis.ctrl.hind <- ncvar_get(ncdata, 'GIS_hind')
  te.ctrl.hind <- ncvar_get(ncdata, 'TE_hind')
  gsic.ctrl.hind <- ncvar_get(ncdata, 'GSIC_hind')
  gmsl.ctrl.rcp26 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
  gmsl.ctrl.rcp45 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
  gmsl.ctrl.rcp85 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
nc_close(ncdata)

# read corresponding ensemble parameters
filename.parameters <- '../output_calibration/BRICK-model_postcalibratedParameters_control_01Nov2016.nc'
ncdata <- nc_open(filename.parameters)
parameters.fit <- ncvar_get(ncdata, 'BRICK_parameters')
parnames.fit <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters.fit <- t(parameters.fit)
colnames(parameters.fit) <- parnames.fit


# Determine suitable bounds for the initial condtions in BRICK v0.2 (fully-
# forward integrating, single fortran function call)


# Gs0 is the initial sea-level contribution from GSIC. It is defined to be 0,
# but there is uncertainty. The previous calibration runs give an idea as to
# what this uncertainty should be.
Gs0.gsic <- gsic.ctrl.hind[1,] - mean(gsic.ctrl.hind[1,])

# V0.gsic is the initial GSIC volume. In the previous calibration runs, this was
# defined in 1960 (Wigley and Raper, 2005), but we need it in 1850 (when the
# BRICK-[other model] runs would start). The previous calibration (as described
# in Wong et al 2017, GMDD) starts GSIC in 1961 and integrates forward
# to 2009, and backward to 1850 to obtain the estimates below for V0.gsic.
V0.gsic <- parameters.fit[,'V0.gsic'] +
            (parameters.fit[,'Gs0'] - gsic.ctrl.hind[1,])

# V0.te is the initial TE sea-level contribution. It is defined to be 0, but
# (like Gs0) there is uncertainty. The previous calibratoin runs give an idea as
# to what this uncertainty is.
V0.te <- te.ctrl.hind[1,] - mean(te.ctrl.hind[1,])

# V0.simple is the initial GIS volume (not sea-level contribution). Previous
# calibration has this defined to be relative to 1961-1990.
i1961 <- which(t.hind==1961)
V0.simple <- parameters.fit[,'V0'] +
            (gis.ctrl.hind[i1961,] - gis.ctrl.hind[1,])

print(paste('Gs0.gsic bounds: ',min(Gs0.gsic),' to ',max(Gs0.gsic),sep=''))
print(paste('V0.gsic bounds: ',min(V0.gsic),' to ',max(V0.gsic),sep=''))
print(paste('V0.te bounds: ',min(V0.te),' to ',max(V0.te),sep=''))
print(paste('V0.simple bounds: ',min(V0.simple),' to ',max(V0.simple),sep=''))


##==============================================================================

# fit some prior distributions

#install.packages('fitdistrplus')
library(fitdistrplus)

# function to scale a random variable x on [a,b] to [0,1]
range.to.beta <- function(x, a, b){
    return((x-a)/(b-a))
}
# function to scale a random variable y on [0,1] to [a,b]
beta.to.range <- function(y, a, b){
    return(y*(b-a)+a)
}

# fit beta distribution as prior

ptmp <- range.to.beta(Gs0.gsic, -0.037, 0.024)
Gs0.gsic.prior <- fitdist(ptmp, 'beta')$estimate

ptmp <- range.to.beta(V0.gsic, 0.31, 0.53)
V0.gsic.prior <- fitdist(ptmp, 'beta')$estimate

itmp <- which(V0.te >= -0.046 & V0.te <= 0.074)
ptmp <- range.to.beta(V0.te[itmp], -0.046, 0.074)
V0.te.prior <- fitdist(ptmp, 'beta')$estimate

itmp <- which(V0.simple >= 7.15 & V0.simple <= 7.58)
ptmp <- range.to.beta(V0.simple[itmp], 7.15, 7.58)
V0.simple.prior <- fitdist(ptmp, 'beta')$estimate

itmp <- which(parameters.fit[,'rho.gsic'] >= 0 & parameters.fit[,'rho.gsic'] <= 1)
ptmp <- range.to.beta(parameters.fit[itmp,'rho.gsic'], 0, 1)
rho.gsic.prior <- fitdist(ptmp, 'beta')$estimate

# sigmaT and sigmaH, fit as variance with inverse gamma
# to match mean and variance of the sample, pick shape1=2+mean/variance and
# shape2=mean*(shape1-1)
xtmp <- parameters.fit[,'sigma.T']^2
sigmaT.prior <- matrix(c(2+((mean(xtmp)^2)/var(xtmp)),
                         mean(xtmp)*(1+((mean(xtmp)^2)/var(xtmp)))), 1,2)
colnames(sigmaT.prior) <- c('shape','rate')

xtmp <- parameters.fit[,'sigma.H']^2
sigmaH.prior <- matrix(c(2+((mean(xtmp)^2)/var(xtmp)),
                         mean(xtmp)*(1+((mean(xtmp)^2)/var(xtmp)))), 1,2)
colnames(sigmaH.prior) <- c('shape','rate')

##==============================================================================

# joint prior for a.te, b.te and invtau.te?

# install.packages('sn')
# install.packages('fMultivar')
library(sn)
library(fMultivar)

itmp <- c(match('a.te',parnames.fit), match('b.te',parnames.fit), match('invtau.te',parnames.fit))
te.prior <- msnFit(parameters.fit[,itmp])
te.prior.fit <- te.prior@fit$estimated
# account for truncated prior
te.prior.fit$cnorm <- as.numeric(pmsn( bound.upper[match(parnames.fit[itmp],parnames)], te.prior.fit$beta, te.prior.fit$Omega, te.prior.fit$alpha) -
                       pmsn( bound.lower[match(parnames.fit[itmp],parnames)], te.prior.fit$beta, te.prior.fit$Omega, te.prior.fit$alpha))

# formulate bivariate gamma?
a.bivg <- 1.4638434
rate.invtau <- 1/scale.invtau
bivg <- function(x,y,shape1,rate1,shape2){
    tmp <- (rate1^shape1)*(x^(shape1+shape2-1))*(y^(shape1-1))*exp(-x*(shape1+y))
    tmp <- tmp/gamma(shape1)
    tmp <- tmp/gamma(shape2)
    return(tmp)
}

##==============================================================================

# joint prior for n.gsic, beta0.gsic, rho.gsic?

itmp <- c(match('beta0',parnames.fit), match('n',parnames.fit), match('rho.gsic',parnames.fit))
gsic.prior <- msnFit(parameters.fit[,itmp])
gsic.prior.fit <- gsic.prior@fit$estimated
# account for truncated prior
# CAREFUL - hardcoded in the change in parnames between old and new
gsic.prior.fit$cnorm <- as.numeric(pmsn( bound.upper[c(6,8,37)], gsic.prior.fit$beta, gsic.prior.fit$Omega, gsic.prior.fit$alpha) -
                       pmsn( bound.lower[c(6,8,37)], gsic.prior.fit$beta, gsic.prior.fit$Omega, gsic.prior.fit$alpha))

##==============================================================================

# joint prior for DOECLIM parameters?

# install.packages('sn')
# install.packages('fMultivar')
library(sn)
library(fMultivar)

itmp <- c(match('S',parnames.fit), match('kappa.doeclim',parnames.fit),
          match('alpha.doeclim',parnames.fit), match('T0',parnames.fit), match('H0',parnames.fit))
doeclim.prior <- msnFit(parameters.fit[,itmp])
doeclim.prior.fit <- doeclim.prior@fit$estimated
# account for truncated prior
# CAREFUL - hardcoded in the change in parnames between old and new
doeclim.prior.fit$cnorm <- as.numeric(pmsn( bound.upper[1:5], doeclim.prior.fit$beta, doeclim.prior.fit$Omega, doeclim.prior.fit$alpha) -
                            pmsn( bound.lower[1:5], doeclim.prior.fit$beta, doeclim.prior.fit$Omega, doeclim.prior.fit$alpha))

##==============================================================================+

# joint prior for a.simple, b.simple, alpha.simple, beta.simple?

itmp <- c(match('a.simple',parnames.fit), match('b.simple',parnames.fit),
          match('alpha.simple',parnames.fit), match('beta.simple',parnames.fit))
simple.prior <- msnFit(parameters.fit[,itmp])
simple.prior.fit <- simple.prior@fit$estimated
# account for truncated prior
simple.prior.fit$cnorm <- as.numeric(pmsn( bound.upper[match(parnames.fit[itmp],parnames)], simple.prior.fit$beta, simple.prior.fit$Omega, simple.prior.fit$alpha) -
                           pmsn( bound.lower[match(parnames.fit[itmp],parnames)], simple.prior.fit$beta, simple.prior.fit$Omega, simple.prior.fit$alpha))

##==============================================================================

# joint prior for anto.a and anto.b?

# first, run precalibration LHS in scratch_findANTObounds_precalibration.R
# will yield tmp.sort, where columns 2 and 3 are anto.a and anto.b, respectively

source('scratch_findANTObounds_precalibration.R')

anto.prior <- msnFit(tmp.sort[itmp,2:3])
anto.prior.fit <- anto.prior@fit$estimated
# account for truncated prior
itmp <- c(match('anto.a',parnames),match('anto.b',parnames))
anto.prior.fit$cnorm <- as.numeric(pmsn( bound.upper[itmp], anto.prior.fit$beta, anto.prior.fit$Omega, anto.prior.fit$alpha) -
                         pmsn( bound.lower[itmp], anto.prior.fit$beta, anto.prior.fit$Omega, anto.prior.fit$alpha))

##==============================================================================
# if you want to visualize any of these priors, fit to the data
X <- parameters.fit[,itmp]
ans <- gsic.prior                                           # <<< for example
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
##==============================================================================

# principal component analysis for dimension reduction in DAIS?

file='/Users/axw322/codes/BRICK_Experiments_backup/output_model/BRICK-model_physical_control_01Nov2016.nc'
ncdata <- nc_open(file)
  ais.hind = ncvar_get(ncdata, 'AIS_hind')
nc_close(ncdata)
thing1 <- cbind(parameters.fit, ais.hind[160,])
colnames(thing1)[40] <- 'ais'
#aisdat <- thing1[,c(26:38,40)]
aisdat <- thing1[,26:38]
ais.pca <- prcomp(log(aisdat), center = TRUE, scale. = TRUE)
##==============================================================================
## End
##==============================================================================
