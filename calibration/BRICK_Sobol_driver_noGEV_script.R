##==============================================================================
## Sobol sensitvity analysis for drivers of flood risk.
## Need to set up the design matrix with each of the three levers of RCP
## scenario, AIS fast dynamics contribution to sea level (in 2065), and surge
## factor parameter (0 or between 1.5 and 2), assocaited with the model response
## of average AEP between 2015 and 2065 (or return period, inversely).
##
## Submit with the following to quiet the output files (except the .log, below):
##  qsub -e /dev/null -o /dev/null rsobol_run.sh
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

rm(list =ls()) #Clear global environment

##########################################
# setting things here, not to be touched #
N.sample000 <- 122000 #122000
N.use000    <- 122000 #122000
N.boot000   <- 50000 #50000
N.core000   <- 15     # number of cores for parallel Sobol analysis
lbuild      <- TRUE
lfullAIS    <- TRUE   # full prior range from Wong et al 2017 (fast dynamics) for h0 and c? (fixed if FALSE)
lfullGEV    <- FALSE   # full uncertainty in GEV parameters? (fixed if FALSE) (not used)
begyear     <- 1850   # this is just beginning of the projections; assessment begins in 2015
endyear     <- 2065   # must be equal to end year of flood risk assessment
appen       <- '-Build-AIS-2065'
name.output.rdata <- 'BRICK_Sobol_2065noGEV.RData'
setwd('/home/scrim/axw322/codes/BRICK/calibration')
##########################################

##==============================================================================
## Need preliminary function to map ranges from [0,1] and back
map.range <- function(X, lbin, ubin, lbout, ubout){
    Y <- lbout + (ubout-lbout)*( (X-lbin)/(ubin-lbin) )
    return(Y)
}
##==============================================================================



##==============================================================================
## Preliminary setup

## Setup packages and libraries
#install.packages('compiler')
#install.packages('pscl')
#library(pscl) # install inverse gamma distribution
#library(compiler)
#enableJIT(3)

## Set the seed (for reproducibility)
set.seed(1234)
##==============================================================================



##==============================================================================
## Model setup and parameters

## Mean and standard deviation for sampling land water storage (LWS)
## contributions to GMSL
## -- Using IPCC AR5 (Church et al. 2013) --
#lws.mean <- 0.38           # mm/y
#lws.sd   <- (0.49-.26)/4   # mm/y (take the IPCC 5-95% range as +/-2sigma)
## -- Using Dieng et al 2015 (doi:10.1088/1748-9326/10/12/124010)--
lws.mean000 <- 0.30           # mm/y
lws.sd000   <- 0.18           # mm/y

## Read the calibrated parameter sets (after rejection sampling to GMSL data)
library(ncdf4)

filename.parameters <- '~/codes/BRICK/output_calibration/BRICK_postcalibratedParameters_fd-gamma_08May2017.nc'
ncdata <- nc_open(filename.parameters)
parnames.brick   <-   ncvar_get(ncdata, 'parnames')
parameters.brick <- t(ncvar_get(ncdata, 'BRICK_parameters'))
nc_close(ncdata)
colnames(parameters.brick) <- parnames.brick

## Remove the statistical parameters
irem <- c(match('sigma.T',parnames.brick),
          match('sigma.H',parnames.brick),
          match('rho.T',parnames.brick),
          match('rho.H',parnames.brick),
          match('sigma.gsic',parnames.brick),
          match('rho.gsic',parnames.brick),
          match('sigma.simple',parnames.brick),
          match('var.dais',parnames.brick))
parnames.brick <- parnames.brick[-irem]
parameters.brick <- parameters.brick[,-irem]

## Source the BRICK models
source('../fortran/R/doeclimF.R')       # the DOECLIM model
source('../fortran/R/GSIC_magiccF.R')   # the GSIC model
source('../fortran/R/brick_te_F.R')     # TE (thermosteric expansion) model
source('../fortran/R/simpleF.R')        # GIS (Greenland Ice Sheet) model
source('../fortran/R/daisantoF.R')		# DAIS (Antarctic Ice Sheet) model
source('../fortran/R/daisanto_fastdynF.R')
source('../R/BRICK_coupledModel_fastdyn.R')

luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # thermal expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = TRUE    # Antarctic ice sheet model
luse.brick = cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais)

## Source some useful functions for manipulating data
source('../R/forcing_total.R')			# function to add up the total forcing
source('../R/compute_indices.R')		# function to determine the model and
						# data indices for comparison

## Set up the models
l.project = TRUE
#begyear = 1850
#endyear = 2065; #if(!l.project & endyear>2009) print('l.project and endyear not compatible')
mod.time= begyear:endyear
begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)

## Set the indices for normalization that are consistent with each data set
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1850),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
		c(which(mod.time==1870),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

## Get the sloep and intercept for global mean surface temperature-Antarctic
## surface temperature relatinship
source('../calibration/DOECLIM_readData.R')
source('../calibration/DAIS_readData.R')

## Set the indices of the initial condition for each sub-model
i0 = vector("list",nrow(ind.norm.data)); names(i0)=as.character(ind.norm.data[,1])

## GSIC initial conditions are actually relative to 1990 (Wigley and Raper 2005)
## Re-set these. The simulation is relative to 1990, but results and comparison
## to data is relative to 1960.
i0$gsic = which(mod.time==1990)

## GIS initial conditions are relative to 1961-1990
i0$gis = which(mod.time==1961)
##==============================================================================



##==============================================================================
## Read the four RCP forcings

forcing.rcp26 = read.csv( '../data/forcing_rcp26.csv', header=TRUE )
forcing.rcp45 = read.csv( '../data/forcing_rcp45.csv', header=TRUE )
forcing.rcp60 = read.csv( '../data/forcing_rcp6.csv', header=TRUE )
forcing.rcp85 = read.csv( '../data/forcing_rcp85.csv', header=TRUE )
##==============================================================================



##==============================================================================
## Fit kernel density estimates to each

source('../calibration/BRICK_parameterSetup_fastdyn.R')

# re-order bound.lower and bound.upper (associated with parnames) to get
# bound.lower.brick and bound.upper.brick to match parnames.brick from file,
# which we will sample from the KDE of
bound.lower.brick <- rep(-999,length(parnames.brick))
bound.upper.brick <- rep(-999,length(parnames.brick))
for (pp in 1:length(parnames.brick)) {
    itmp <- match(parnames.brick[pp], parnames)
    bound.lower.brick[pp] <- bound.lower[itmp]
    bound.upper.brick[pp] <- bound.upper[itmp]
}

# remove h0 and c to test sensitivity
if (!lfullAIS) {
  bound.lower.brick <- bound.lower.brick[-match('h0',parnames.brick)]
  bound.upper.brick <- bound.upper.brick[-match('h0',parnames.brick)]
  parameters.brick  <- parameters.brick[,-match('h0',parnames.brick)]
  parnames.brick    <- parnames.brick[-match('h0',parnames.brick)]
  bound.lower.brick <- bound.lower.brick[-match('c',parnames.brick)]
  bound.upper.brick <- bound.upper.brick[-match('c',parnames.brick)]
  parameters.brick  <- parameters.brick[,-match('c',parnames.brick)]
  parnames.brick    <- parnames.brick[-match('c',parnames.brick)]
}

# gamma priors might be outside of the prior range. cut them
itmp <- match('Tcrit',parnames.brick)
irem <- which(parameters.brick[,'Tcrit'] < bound.lower.brick[itmp] |
              parameters.brick[,'Tcrit'] > bound.upper.brick[itmp] )
if(length(irem)>0) {parameters.brick <- parameters.brick[-irem,]}
itmp <- match('lambda',parnames.brick)
irem <- which(parameters.brick[,'lambda'] < bound.lower.brick[itmp] |
              parameters.brick[,'lambda'] > bound.upper.brick[itmp] )
if(length(irem)>0) {parameters.brick <- parameters.brick[-irem,]}

# fit KDEs and store in list.
# 1. Box kernels using 'hist' -- using freq=false gives a box kernel density
# estimate, which integrates to 1 against dx, which is calculated and stored.
# 2. Gaussian (or other) kernels using 'density' -- use Gaussian kernels and
# restrict to the bounds of the posterior parameter estimates. Permits values
# outside the posterior range observed by drawing from the Gaussian kernels, but
# will restrict to the prior ranges still (to avoid unreasonable model output).
# 3. Note that 'nnode' is only for visualization purposes - the actual KDE
# object is fit to the *data* not the nodes.
kerntype <- 'gaussian'

kde.brick <- vector('list', length(parnames.brick))
kde.bw <- rep(-999, length(parnames.brick))
names(kde.brick) <- parnames.brick
nnode <- 1024

for (pp in 1:length(parnames.brick)) {
  if(kerntype=='box'){
#    kde.brick[[pp]] <- hist(parameters.brick[,pp],
#                            breaks=seq(from=bound.lower.brick[pp],
#                                       to=bound.upper.brick[pp],
#                                       length.out=nnode),
#                            freq=FALSE,
#                            plot=FALSE, warn.unused=FALSE)
#    kde.brick[[pp]]$dx <- median(diff(kde.brick[[pp]]$mids))
  } else if(kerntype=='gaussian') {
    kde.brick[[pp]] <- density(x=parameters.brick[,pp],
                               from=min(parameters.brick[,pp]),
                               to=max(parameters.brick[,pp]),
#                               from=bound.lower.brick[pp],
#                               to=bound.upper.brick[pp],
                               n=nnode,
                               kernel=kerntype)
    kde.bw[pp] <- kde.brick[[pp]]$bw
  } else {
    print('ERROR - unrecognized kerntype')
  }
}
##==============================================================================



##==============================================================================
## Sample the parameters to create an ensemble of BRICK model parameters
## n.ensemble gives the size of each of the two Sobol' samples

#n.sample <- 500
n.sample <- N.sample000

## Sample BRICK parameters
## And note that the kernel here is the box of width 'dx', so to draw from the
## multivariate KDE, pick an index i from {1,2,...,N.ensemble} and then draw
## uniform random from each of [lower.bound[pp]_i , upper.bound[pp]_i]
## And using box kernel, so should not find any parameter values outside of the
## bounds.
n.ensemble <- nrow(parameters.brick)
parameters.sample1.brick <- mat.or.vec(n.sample, ncol(parameters.brick))
parameters.sample2.brick <- mat.or.vec(n.sample, ncol(parameters.brick))
isample1 <- sample(x=seq(1,n.ensemble), size=n.sample, replace=TRUE)
isample2 <- sample(x=seq(1,n.ensemble), size=n.sample, replace=TRUE)

print(paste('Sampling ',n.sample,' BRICK parameters from ',kerntype,' kernels...',sep=''))
pb <- txtProgressBar(min=0,max=n.sample,initial=0,style=3)
if(kerntype=='box') {
    for (i in 1:n.sample) {
        p1.brick <- as.numeric(parameters.brick[isample1[i],])
        p2.brick <- as.numeric(parameters.brick[isample2[i],])
        for (pp in 1:length(p1.brick)) {
            # find which bin each parameter of this set is in
            ibin1 <- which.min(abs(as.numeric(p1.brick[pp]) - kde.brick[[pp]]$mids))
            ibin2 <- which.min(abs(as.numeric(p2.brick[pp]) - kde.brick[[pp]]$mids))
            # draw uniform randomly from this bin
            parameters.sample1.brick[i,pp] <- runif(n=1,
                                                min=kde.brick[[pp]]$breaks[ibin1],
                                                max=kde.brick[[pp]]$breaks[ibin1+1])
            parameters.sample2.brick[i,pp] <- runif(n=1,
                                                min=kde.brick[[pp]]$breaks[ibin2],
                                                max=kde.brick[[pp]]$breaks[ibin2+1])
            # map these to [0,1] for the Sobol
            parameters.sample1.brick[i,pp] <- map.range(parameters.sample1.brick[i,pp],
                                                    lbin =bound.lower.brick[pp],
                                                    ubin =bound.upper.brick[pp],
                                                    lbout=0,
                                                    ubout=1)
            parameters.sample2.brick[i,pp] <- map.range(parameters.sample2.brick[i,pp],
                                                    lbin =bound.lower.brick[pp],
                                                    ubin =bound.upper.brick[pp],
                                                    lbout=0,
                                                    ubout=1)
        }
        setTxtProgressBar(pb, i)
    }
} else if(kerntype=='gaussian') {
    for (i in 1:n.sample) {
        p1.brick <- as.numeric(parameters.brick[isample1[i],])
        p2.brick <- as.numeric(parameters.brick[isample2[i],])
        for (pp in 1:length(p1.brick)) {
            # draw sample distributed normally with standard deviation = bandwidth
            # from density estimation (above), and centered at the sampled parameter
            # values.
            parameters.sample1.brick[i,pp] <- rnorm(n=1, mean=p1.brick[pp], sd=kde.bw[pp])
            parameters.sample2.brick[i,pp] <- rnorm(n=1, mean=p2.brick[pp], sd=kde.bw[pp])

            # re-sample for parameters.sample1.brick if outside of prior ranges
            lresample <- parameters.sample1.brick[i,pp] > bound.upper.brick[pp] |
                         parameters.sample1.brick[i,pp] < bound.lower.brick[pp]
            while(lresample) {
                parameters.sample1.brick[i,pp] <- rnorm(n=1, mean=p1.brick[pp], sd=kde.bw[pp])
                lresample <- parameters.sample1.brick[i,pp] > bound.upper.brick[pp] |
                             parameters.sample1.brick[i,pp] < bound.lower.brick[pp]
            }
            # re-sample for parameters.sample2.brick if outside of prior ranges
            lresample <- parameters.sample2.brick[i,pp] > bound.upper.brick[pp] |
                         parameters.sample2.brick[i,pp] < bound.lower.brick[pp]
            while(lresample) {
                parameters.sample2.brick[i,pp] <- rnorm(n=1, mean=p2.brick[pp], sd=kde.bw[pp])
                lresample <- parameters.sample2.brick[i,pp] > bound.upper.brick[pp] |
                             parameters.sample2.brick[i,pp] < bound.lower.brick[pp]
            }

            # map to [0,1] for Sobol
            parameters.sample1.brick[i,pp] <- map.range(parameters.sample1.brick[i,pp],
                                                        lbin =bound.lower.brick[pp],
                                                        ubin =bound.upper.brick[pp],
                                                        lbout=0,
                                                        ubout=1)
            parameters.sample2.brick[i,pp] <- map.range(parameters.sample2.brick[i,pp],
                                                        lbin =bound.lower.brick[pp],
                                                        ubin =bound.upper.brick[pp],
                                                        lbout=0,
                                                        ubout=1)
        }
        setTxtProgressBar(pb, i)
    }
}
close(pb)
print('... done.')

## Add to these the storm surge rise parameters. Draw from [0,1] because that is
## what the Sobol expects, then map to [bound.lower/upper.surge] during function
## call.
bound.lower.surge <- 0
bound.upper.surge <- 2
# Draw the parameters as part of Latin hypercube with RCP later
#parameters.sample1.surge <- runif(n=n.sample, 0, 1)
#parameters.sample2.surge <- runif(n=n.sample, 0, 1)

## Add to these the Van Dantzig parameters. Draw from [0,1] and map to assumed
## distributions during function call.
## Not actually using the first 5 (related to economic optimization), so only
## subsidence needed.
p_zero_p = 0.0038              # Initial flood frequency (1/yr) with zero height increase (Van Dantzig (1956))
alpha_p = 2.6                  # Exponential flood frequency constant (Van Dantzig (1956))
V_p = c(7.5e+9, 3e+10)         # Value of goods protected by dike (based on estimates in Jonkman et al. (2009)) (US$)
delta_prime = c(0.02, 0.06)    # Discount rate (percent/year) (based on estimates in Jonkman et al. (2009))
I_range = c(-0.5,1.0)          # Investment cost uncertainty range (as fraction of investment cost) (Jonkman and Dutch Perspective use -50%, +100%)
sub_rate = 0.0056              # Rate of land subsidence (meter/year) (Dixon et al. (2006))
#bound.lower.vd <- c( 0 , -Inf, 7.5e9, 0.02, -0.5, 0  )
#bound.upper.vd <- c(Inf,  Inf, 3e10 , 0.06,  1.0, Inf)
bound.lower.subs <- 0
bound.upper.subs <- Inf
#parameters.sample1.vd <- randomLHS(n.sample, 6)
#parameters.sample2.vd <- randomLHS(n.sample, 6)
#parnames.vd <- c("p0.vd", "alpha.vd", "V0.vd", "delta.vd", "Iunc.vd", "sub_rate.vd")
# Draw the parameters as part of Latin hypercube with RCP later
#parameters.sample1.subs <- runif(n=n.sample, 0, 1)
#parameters.sample2.subs <- runif(n=n.sample, 0, 1)

## Add GEV parameter fits for storm surge (stationary)
## draw quantiles from U[0,1], 3D Latin hypercube (to get good coverage) and map
## to the parameter estimates corresponding to those quantiles within the
## BRICK_sobol call

filename.gevstat <- '../output_calibration/BRICK_estimateGEV-AnnMean_12Apr2017.nc'
ncdata <- nc_open(filename.gevstat)
parnames.gev <- ncvar_get(ncdata, 'GEV_names')
gev.mcmc <- t(ncvar_get(ncdata, 'GEV_parameters'))
nc_close(ncdata)
colnames(gev.mcmc) <- parnames.gev

## Don't allow shape (xi) < 0, because that changes the fundamental behavior
irem <- which(gev.mcmc[,'shape'] <= 0)
if(length(irem) > 0) gev.mcmc <- gev.mcmc[-irem,]
bound.lower.gev <- apply(gev.mcmc, 2, min)
bound.upper.gev <- apply(gev.mcmc, 2, max)
iloc <- match('location',parnames.gev)
isha <- match('shape',parnames.gev)
isca <- match('scale',parnames.gev)

# Do not draw as LHS - this misses the correlation structure that was obtained
# by calibrating using extRemes package MCMC.
isample1 <- sample(seq(1, nrow(gev.mcmc)), size=n.sample, replace=TRUE)
isample2 <- sample(seq(1, nrow(gev.mcmc)), size=n.sample, replace=TRUE)
parameters.sample1.gev <- gev.mcmc[isample1,]
parameters.sample2.gev <- gev.mcmc[isample2,]
# Map parameters to [0,1] for Sobol
for (pp in 1:3) {
    parameters.sample1.gev[,pp] <- map.range(parameters.sample1.gev[,pp],
                                             lbin =bound.lower.gev[pp],
                                             ubin =bound.upper.gev[pp],
                                             lbout=0,
                                             ubout=1)
    parameters.sample2.gev[,pp] <- map.range(parameters.sample2.gev[,pp],
                                             lbin =bound.lower.gev[pp],
                                             ubin =bound.upper.gev[pp],
                                             lbout=0,
                                             ubout=1)
}

## Add land water storage to parameter bounds. Sampling from normal distribution
## so these don't really matter
bound.lower.lws <- -Inf
bound.upper.lws <- Inf

## Add build scenario to parameters and bounds. Sample from [-3,3] feet
bound.lower.build <- -3*.3048
bound.upper.build <- 3*.3048

## Finally, add RCP scenario to parameters and bounds. Sample from [0,1], and
## each 1/4 gives different scenario (2.6, 4.5, 6.0, 8.5).
bound.lower.rcp <- 0
bound.upper.rcp <- 1

# Draw the surge and subsidence parameters as part of Latin hypercube with RCP
# Also draw build heights from [-3 ft to +3 ft]. Note that obviously we will not
# build -3 ft, but this is to assess the impact of varying levee height on the
# flood risk.
library(lhs)

lhs.sample1 <- randomLHS(n.sample, 5)
lhs.sample2 <- randomLHS(n.sample, 5)
parameters.sample1.rcp <- lhs.sample1[,1]
parameters.sample2.rcp <- lhs.sample2[,1]
parameters.sample1.surge <- lhs.sample1[,2]
parameters.sample2.surge <- lhs.sample2[,2]
parameters.sample1.subs <- lhs.sample1[,3]
parameters.sample2.subs <- lhs.sample2[,3]
parameters.sample1.build <- lhs.sample1[,4]
parameters.sample2.build <- lhs.sample2[,4]
parameters.sample1.lws <- lhs.sample1[,5]
parameters.sample2.lws <- lhs.sample2[,5]

## add Van Dantzig, surge, and RCP to parnames and bounds (and optionally build)
if (lbuild) {
  parameters.sobol1 <- cbind(parameters.sample1.brick, parameters.sample1.lws  ,
                             parameters.sample1.surge, parameters.sample1.subs ,
                             parameters.sample1.rcp  , parameters.sample1.build)
  parameters.sobol2 <- cbind(parameters.sample2.brick, parameters.sample2.lws  ,
                             parameters.sample2.surge, parameters.sample2.subs ,
                             parameters.sample2.rcp  , parameters.sample2.build)
  bound.lower.sobol <- c(bound.lower.brick, bound.lower.lws, bound.lower.surge ,
                         bound.lower.subs , bound.lower.rcp, bound.lower.build )
  bound.upper.sobol <- c(bound.upper.brick, bound.upper.lws, bound.upper.surge ,
                         bound.upper.subs , bound.upper.rcp, bound.upper.build )
  parnames.sobol <- c(parnames.brick, 'lws.mean', 'surge.factor', 'subs.rate', 'rcp','build')
  ind.brick <- 1:length(parnames.brick)
  ind.lws   <- length(parnames.brick)+1
  ind.surge <- match('surge.factor',parnames.sobol)
  ind.subs  <- match('subs.rate'   ,parnames.sobol)
  ind.rcp   <- length(parnames.sobol)-1
  ind.build <- length(parnames.sobol)
} else {
  parameters.sobol1 <- cbind(parameters.sample1.brick, parameters.sample1.lws  ,
                             parameters.sample1.surge, parameters.sample1.subs ,
                             parameters.sample1.rcp  )#, parameters.sample1.build)
  parameters.sobol2 <- cbind(parameters.sample2.brick, parameters.sample2.lws  ,
                             parameters.sample2.surge, parameters.sample2.subs ,
                             parameters.sample2.rcp  )#, parameters.sample2.build)
  bound.lower.sobol <- c(bound.lower.brick, bound.lower.lws, bound.lower.surge ,
                         bound.lower.subs , bound.lower.rcp)#, bound.lower.build )
  bound.upper.sobol <- c(bound.upper.brick, bound.upper.lws, bound.upper.surge ,
                         bound.upper.subs , bound.upper.rcp)#, bound.upper.build )
  parnames.sobol <- c(parnames.brick, 'lws.mean', 'surge.factor', 'subs.rate', 'rcp')#,'build')
  ind.brick <- 1:length(parnames.brick)
  ind.lws   <- length(parnames.brick)+1
  ind.surge <- match('surge.factor',parnames.sobol)
  ind.subs  <- match('subs.rate'   ,parnames.sobol)
  ind.rcp   <- length(parnames.sobol)
}

parameters.sobol1 <- data.frame(parameters.sobol1)
parameters.sobol2 <- data.frame(parameters.sobol2)
colnames(parameters.sobol1) <- parnames.sobol
colnames(parameters.sobol2) <- parnames.sobol
##==============================================================================



##==============================================================================
## Define model for Sobol
## Needs to take in a data frame, with each column as a different model
## parameter, and return an output data frame with same number of rows (one
## column response vector).

## Define some global variables to use within the BRICK_sobol function calls
# response over 2015-2065
i2015 <- which(mod.time==2015)
imodend <- which(mod.time==endyear)
nt    <- 0:(imodend-i2015)

# Fingerprints of sea-level rise sources on local sea-level rise
lat.fp = 29.95			# latitude of location to fingerprint local sea level rise (>0 is North, <0 is South)
lon.fp = -90.07			# longitude of location ... (>0 is East, <0 is West)
filename.fingerprints = "../fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"
ncdata <- nc_open(filename.fingerprints)
  lat = ncvar_get(ncdata, 'lat')
  lon = ncvar_get(ncdata, 'lon')
  fp.gsic = ncvar_get(ncdata, 'GLAC')
  fp.gis = ncvar_get(ncdata, 'GIS')
  fp.ais = ncvar_get(ncdata, 'AIS')
nc_close(ncdata)
# convert longitude to degrees East, and find the fingerprinting data location
# closest to the local sea level lat/lon given
if(lon.fp < 0) {lon.fp=lon.fp+360}	# convert longitude to degrees East
ilat = which( abs(lat-lat.fp)==min(abs(lat-lat.fp)) )
ilon = which( abs(lon-lon.fp)==min(abs(lon-lon.fp)) )
# it is possible there were multiple lats/lons 'closest' to your given point
# take the average of the non-NA of these
fp.loc      <- vector('list',5); names(fp.loc) <- c('ais','gis','gsic','te','lws')
fp.loc$ais  <- mean(fp.ais[ilon,ilat],na.rm=TRUE)
fp.loc$gsic <- mean(fp.gsic[ilon,ilat],na.rm=TRUE)
fp.loc$gis  <- mean(fp.gis[ilon,ilat],na.rm=TRUE)
fp.loc$te   <- 1.0		# TE response is to global mean temperature, so global mean sea level response is same everywhere
fp.loc$lws  <- 1.0		# assume LWS is uniform

# levee height
H0 <- 16*0.3048 # initial levee height (m) (change to 16ft = 4.877m
                            # rough average of NOLA central levee ring)
H0 <- H0-(4*0.3048)        # initial subsidence of 4 ft

library(fExtremes)

# using mapply is generally faster than a 'for' loop, but no progress bar, which
# is comforting

# fixed at the MLE GEV parameters
parameters.gev.fixed <- c(577.051961, 141.104958,   0.360476) #location, scale, shape

export.names <- c('brick_model', 'doeclimF', 'gsic_magiccF', 'simpleF',
                  'brick_te_F', 'daisanto_fastdynF', 'anto', 'forcing_total',
                  'flux.to.heat', 'parnames.brick', 'slope.Ta2Tg',
                  'intercept.Ta2Tg','mod.time','ind.norm.data', 'i2015',
                  'imodend', 'fp.loc', 'nt', 'H0', 'isca', 'iloc', 'isha', 'parameters.gev.fixed',
                  'ind.norm','luse.brick','i0','lbuild',
                  'lws.mean000', 'lws.sd000',
                  'forcing.rcp26', 'forcing.rcp45',
                  'forcing.rcp60', 'forcing.rcp85')



##==================##
## parallel version ##
##==================##

brick_sobol_par <- function(dataframe.in){

    # initialize output
    nr            <- nrow(dataframe.in)
    lsl.proj      <- rep(0,length(nt))
    surge.rise    <- rep(0,length(nt))
    lsl.norm.subs <- rep(0,length(nt))
    H_eff         <- rep(0,length(nt))
    output        <- rep(0,nr)

    # map input from [0,1] range to actual parameters
    parameters.brick  <- dataframe.in[,ind.brick]
    for (pp in 1:length(ind.brick)) {
        parameters.brick[,pp] <- map.range(parameters.brick[,pp], 0, 1, bound.lower.brick[pp], bound.upper.brick[pp])
    }

    lws.mean <- qnorm(dataframe.in[,ind.lws], mean=lws.mean000, sd=lws.sd000)

    surge.factor <- map.range(dataframe.in[,ind.surge], 0, 1, bound.lower.surge, bound.upper.surge)

    subs.rate <- qlnorm(dataframe.in[,ind.subs], log(sub_rate), 0.4) # sdlog to yield about stdev from Dixon et al. (2006), 2.5mm/y

    parameters.rcp <- dataframe.in[,ind.rcp]

    if (lbuild) {
      parameters.build <- map.range(dataframe.in[,ind.build], 0, 1, bound.lower.build, bound.upper.build)
    } else {
      parameters.build <- rep(0,nr)
    }

    # testing in parallel
    #install.packages('foreach')
    #install.packages('doParallel')
    library(foreach)
    library(doParallel)
    cores=detectCores()
#    cl <- makeCluster(cores[1]-1) #not to overload your computer
    cl <- makeCluster(N.core000)
    print(paste('Starting cluster with ',N.core000,' cores', sep=''))
    registerDoParallel(cl)

    # make sea-level rise projections
    print(paste('Starting ',nr,' model projections (inluding flood risk)...',sep=''))
    finalOutput <- foreach(i=1:nr, .combine=c,
                                   .packages='fExtremes',
                                   .export=export.names,
                                   .inorder=FALSE) %dopar% {

        dyn.load("../fortran/doeclim.so")
        dyn.load("../fortran/gsic_magicc.so")
        dyn.load("../fortran/brick_te.so")
        dyn.load("../fortran/dais_fastdyn.so")
        dyn.load("../fortran/simple.so")

        # map RCP
        if(parameters.rcp[i] >= 0 & parameters.rcp[i] <= 0.25) {
            forcing <- forcing.rcp26
        } else if(parameters.rcp[i] > 0.25 & parameters.rcp[i] <= 0.5) {
            forcing <- forcing.rcp45
        } else if(parameters.rcp[i] > 0.5 & parameters.rcp[i] <= 0.75) {
            forcing <- forcing.rcp60
        } else if(parameters.rcp[i] > 0.75 & parameters.rcp[i] <= 1) {
            forcing <- forcing.rcp85
        }
        brick.out <- brick_model(parameters.in      = as.numeric(parameters.brick[i,]),
                                 parnames.in        = parnames.brick,
				 forcing.in         = forcing,
				 l.project          = TRUE,
				 slope.Ta2Tg.in     = slope.Ta2Tg,
				 intercept.Ta2Tg.in = intercept.Ta2Tg,
				 mod.time           = mod.time,
				 ind.norm.data      = ind.norm.data,
				 ind.norm.sl        = ind.norm,
				 luse.brick         = luse.brick,
				 i0                 = i0)
        slr.gsic <- brick.out$gsic.out[i2015:imodend]           - brick.out$gsic.out[i2015]
        slr.te   <- brick.out$te.out[i2015:imodend]             - brick.out$te.out[i2015]
        slr.gis  <- brick.out$simple.out$sle.gis[i2015:imodend] - brick.out$simple.out$sle.gis[i2015]
        slr.ais  <- brick.out$dais.out$Vais[i2015:imodend]      - brick.out$dais.out$Vais[i2015]
        slr.lws  <- cumsum(rnorm(n=length(mod.time), mean=lws.mean[i], sd=lws.sd000)) /1000
        slr.lws  <- slr.lws[i2015:imodend] - slr.lws[i2015]
        lsl.proj <- fp.loc$gsic*slr.gsic + fp.loc$te  *slr.te  +
                    fp.loc$gis *slr.gis  + fp.loc$ais *slr.ais + fp.loc$lws*slr.lws
        surge.rise <- surge.factor[i] * lsl.proj
        lsl.norm.subs <- subs.rate[i]*nt + lsl.proj + surge.rise
        H_eff <- H0 - lsl.norm.subs + parameters.build[i]

        # failure probability is the survival function of the storm surge GEV at effective dike height
        # "by hand" calculation faster; assumes shape parameter /= 0 (which is
        # okay, given that we use 0 as lower bound on prior for it)
        p_fail <- as.numeric(1-pgev(q=1000*H_eff, xi=parameters.gev.fixed[isha], mu=parameters.gev.fixed[iloc], beta=parameters.gev.fixed[isca]))
        #ttmp <- (1+parameters.gev.fixed[i,isha]*((1000*H_eff-parameters.gev.fixed[i,iloc])/parameters.gev.fixed[i,isca]))^(-1/parameters.gev.fixed[i,isha])
        #p_fail <- 1-exp(-ttmp)

        # using average annual exceedance probability as the response
        output[i] <- mean(p_fail)
    }
    print(paste(' ... done.'))

    stopCluster(cl)

    # for Sobol, output must be centered at 0
    output.avg <- mean(finalOutput)
    finalOutput <- finalOutput - output.avg
    return(finalOutput)
}

##==============================================================================








##==============================================================================
## Actually run the Sobol
## First-order indices give size of effects on the output response variance
## that are directly (first-order) related to the variable in question.
## Total-order indices give the size of the effects that are related to the
## variable and all interactions between that variable and others, regardless
## of the order.

library(sensitivity)

#n.use <- 100
n.use <- N.use000
if(n.use > nrow(parameters.sobol1)) {print('ERROR - asking to use more samples than you drew')}

# get first-, second- and total-order indices
print(paste('Starting Sobol analysis using samples of size ',n.use,sep=''))
print(paste('  will require ',n.use*(2*ncol(parameters.sobol1)+2),' simulations',sep=''))

t.out <- system.time(s.out <- sobolSalt(model=brick_sobol_par,
                             parameters.sobol1[1:n.use,],
                             parameters.sobol2[1:n.use,],
                             scheme='B',
                             nboot=N.boot000))

print(paste('Sobol analysis took ',signif(as.numeric(t.out[3]/60),4),' minutes to complete.',sep=''))

#system.time(s.out <- sobol(model=brick_sobol_ser,
#                             parameters.sobol1[1:1000,],
#                             parameters.sobol2[1:1000,],
#                             order=2,
#                             nboot=100))

#s.total$T[rev(order(s.total$T[,1])),]
#s.total$S[rev(order(s.total$S[,1])),]

#plot(s.out, choice=1)

if(is.null(name.output.rdata)) {name.output.rdata <- 'BRICK_Sobol_tmp.RData'}
save.image(file = name.output.rdata)
##==============================================================================



##==============================================================================
## Write an output file like the modified code from Perry will expect

##TODO
##TODO -- check how 'Sx_conf' is used; make sure you feed it in correctly
##TODO

today=Sys.Date(); today=format(today,format="%d%b%Y")
file.sobolout1 <- paste('../output_calibration/BRICK_Sobol-1-tot_',today,appen,'.txt',sep='')
file.sobolout2 <- paste('../output_calibration/BRICK_Sobol-2_',today,appen,'.txt',sep='')

headers.1st.tot <- matrix(c('Parameter', 'S1', 'S1_conf_low', 'S1_conf_high',
                            'ST', 'ST_conf_low', 'ST_conf_high'), nrow=1)
output.1st.tot  <- data.frame(cbind( parnames.sobol,
                                     s.out$S[,1],
                                     s.out$S[,4],
                                     s.out$S[,5],
                                     s.out$T[,1],
                                     s.out$T[,4],
                                     s.out$T[,5]))
write.table(headers.1st.tot, file=file.sobolout1, append=FALSE, sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.1st.tot , file=file.sobolout1, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)

headers.2nd     <- matrix(c('Parameter_1', 'Parameter_2', 'S2', 'S2_conf_low',
                            'S2_conf_high'), nrow=1)
output2.indices <- s.out$S2[,1]
output2.conf1   <- s.out$S2[,4]
output2.conf2   <- s.out$S2[,5]

# 2nd order index names ordered as: (assuming 39 parameters)
# 1. parnames.sobol[1]-parnames.sobol[2]
# 2. parnames.sobol[1]-parnames.sobol[3]
# 3. parnames.sobol[1]-parnames.sobol[4]
# ... etc ...
# 38. parnames.sobol[1]-parnames.sobol[39] << N=2:39 => p1-p[N]
# 39. parnames.sobol[2]-parnames.sobol[3]
# 40. parnames.sobol[2]-parnames.sobol[4]
# 38+37. parnames.sobol[2]-parnames.sobol[39] << N=3:39 => p2-p[N]
# ... etc ...
names2  <- rownames(s.out$S2)
names2a <- rep(NA, length(names2))
names2b <- rep(NA, length(names2))
cnt <- 1
for (i in seq(from=1, to=(length(parnames.sobol)-1), by=1)) {           # i = index of first name
    for (j in seq(from=(i+1), to=(length(parnames.sobol)), by=1)) {   # j = index of second name
        names2a[cnt] <- parnames.sobol[i]
        names2b[cnt] <- parnames.sobol[j]
        cnt <- cnt+1
    }
}

output.2nd <- data.frame(cbind( names2a,
                                names2b,
                                output2.indices,
                                output2.conf1,
                                output2.conf2 ))
write.table(headers.2nd    , file=file.sobolout2, append=FALSE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.2nd     , file=file.sobolout2, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
##==============================================================================


##==============================================================================
## End
##==============================================================================
