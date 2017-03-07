##==============================================================================
## Sobol sensitvity analysis for drivers of flood risk.
## Need to set up the design matrix with each of the three levers of RCP
## scenario, AIS fast dynamics contribution to sea level (in 2065), and surge
## factor parameter (0 or between 1.5 and 2), assocaited with the model response
## of average AEP between 2015 and 2065 (or return period, inversely).
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
library(pscl) # install inverse gamma distribution
library(compiler)
enableJIT(3)

## Set the seed (for reproducibility)
set.seed(1234)
##==============================================================================



##==============================================================================
## Model setup and parameters

## Read the calibrated parameter sets (after rejection sampling to GMSL data)
filename.parameters <- '~/codes/BRICK/output_calibration/BRICK-fastdyn_postcalibratedParameters_gamma_31Jan2017.csv'
parameters.brick <- read.csv(filename.parameters)
parnames.brick <- colnames(parameters.brick)

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
begyear = 1850
endyear = 2065; #if(!l.project & endyear>2009) print('l.project and endyear not compatible')
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

n.sample <- 10000

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
for (i in 1:n.sample) {
    p1.brick <- as.numeric(parameters.brick[isample1[i],])
    p2.brick <- as.numeric(parameters.brick[isample2[i],])
    if(kerntype=='box') {
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
    } else if(kerntype=='gaussian') {
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
    }
    setTxtProgressBar(pb, i)
}
close(pb)
print('... done.')

## Add to these the storm surge rise parameters. Draw from [0,1] because that is
## what the Sobol expects, then map to [bound.lower/upper.surge] during function
## call.
bound.lower.surge <- 0
bound.upper.surge <- 2
parameters.sample1.surge <- runif(n=n.sample, 0, 1)
parameters.sample2.surge <- runif(n=n.sample, 0, 1)

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
parameters.sample1.subs <- runif(n=n.sample, 0, 1)
parameters.sample2.subs <- runif(n=n.sample, 0, 1)

## Add GEV parameter fits for storm surge (stationary)
## draw quantiles from U[0,1], 3D Latin hypercube (to get good coverage) and map
## to the parameter estimates corresponding to those quantiles within the
## BRICK_sobol call
filename.gevstat <- '../output_calibration/BRICK_estimateGEV-AnnMean_07Feb2017.nc'
ncdata <- nc_open(filename.gevstat)
parnames.gev <- ncvar_get(ncdata, 'GEV_names')
gev.mcmc <- t(ncvar_get(ncdata, 'GEV_parameters'))
nc_close(ncdata)
colnames(gev.mcmc) <- parnames.gev

## Don't allow shape (xi) < 0, because that changes the fundamental behavior
irem <- which(gev.mcmc[,'shape'] <= 0)
gev.mcmc <- gev.mcmc[-irem,]

parameters.sample1.gev <- randomLHS(n.sample, 3)
parameters.sample2.gev <- randomLHS(n.sample, 3)
bound.lower.gev <- c(-Inf,  0 , -Inf)
bound.upper.gev <- c( Inf, Inf, Inf )

## Finally, add RCP scenario to parameters and bounds. Sample from [0,1], and
## each 1/4 gives different scenario (2.6, 4.5, 6.0, 8.5).
bound.lower.rcp <- 0
bound.upper.rcp <- 1
parameters.sample1.rcp <- runif(n=n.sample, 0, 1)
parameters.sample2.rcp <- runif(n=n.sample, 0, 1)

## add Van Dantzig, surge, and RCP to parnames and bounds
parameters.sobol1 <- cbind(parameters.sample1.brick, parameters.sample1.surge,
                           parameters.sample1.subs , parameters.sample1.gev  ,
                           parameters.sample1.rcp)
parameters.sobol2 <- cbind(parameters.sample2.brick, parameters.sample2.surge,
                           parameters.sample2.subs , parameters.sample2.gev  ,
                           parameters.sample2.rcp)
bound.lower.sobol <- c(bound.lower.brick, bound.lower.surge, bound.lower.subs, bound.lower.gev, bound.lower.rcp)
bound.upper.sobol <- c(bound.upper.brick, bound.upper.surge, bound.upper.subs, bound.upper.gev, bound.upper.rcp)
parnames.sobol <- c(parnames.brick, 'surge.factor', 'subs.rate', parnames.gev, 'rcp')
ind.brick <- 1:length(parnames.brick)
ind.surge <- match('surge.factor',parnames.sobol)
ind.subs  <- match('subs.rate'   ,parnames.sobol)
ind.gev   <- match(parnames.gev[1], parnames.sobol):match(parnames.gev[3], parnames.sobol)
ind.rcp   <- length(parnames.sobol)

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
i2065 <- which(mod.time==2065)
nt    <- 0:(i2065-i2015)

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
fp.ais.loc = mean(fp.ais[ilon,ilat],na.rm=TRUE)
fp.gsic.loc = mean(fp.gsic[ilon,ilat],na.rm=TRUE)
fp.gis.loc = mean(fp.gis[ilon,ilat],na.rm=TRUE)
fp.te.loc = 1.0		# TE response is to global mean temperature, so global mean sea level response is same everywhere

# levee height
H0 <- 16*0.3048 # initial levee height (m) (change to 16ft = 4.877m
                            # rough average of NOLA central levee ring)
H0 <- H0-(4*0.3048)        # initial subsidence of 4 ft

library(fExtremes)

# using mapply is generally faster than a 'for' loop
brick_sobol <- function(dataframe.in){

    # initialize output
    nr         <- nrow(dataframe.in)
    brick.out  <- vector("list", nr)
    lsl.proj   <- mat.or.vec(length(nt), nr)
    surge.rise <- mat.or.vec(length(nt), nr)
    output     <- rep(0,nr)

    # map input from [0,1] range to actual parameters
    parameters.brick  <- dataframe.in[,ind.brick]
    for (pp in 1:length(ind.brick)) {
        parameters.brick[,pp] <- map.range(parameters.brick[,pp], 0, 1, bound.lower.brick[pp], bound.upper.brick[pp])
    }

    surge.factor <- map.range(dataframe.in[,ind.surge], 0, 1, bound.lower.surge, bound.upper.surge)

    subs.rate <- qlnorm(dataframe.in[,ind.subs], log(sub_rate), 0.4) # sdlog to yield about stdev from Dixon et al. (2006), 2.5mm/y

    parameters.rcp <- dataframe.in[,ind.rcp]

    parameters.gev <- dataframe.in[,ind.gev]
    parameters.gev[,match('location',parnames.gev)] <- quantile(gev.mcmc[,'location'], parameters.gev[,match('location',parnames.gev)])
    parameters.gev[,match('shape',parnames.gev)]    <- quantile(gev.mcmc[,'shape']   , parameters.gev[,match('shape',parnames.gev)])
    parameters.gev[,match('scale',parnames.gev)]    <- quantile(gev.mcmc[,'scale']   , parameters.gev[,match('scale',parnames.gev)])

    # make sea-level rise projections
    print(paste('Starting ',nr,' model projections...',sep=''))
    pb <- txtProgressBar(min=0,max=nr,initial=0,style=3)
    for (i in 1:nr) {
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
        brick.out[[i]] <- brick_model(parameters.in      = as.numeric(parameters.brick[i,]),
                                     parnames.in        = parnames.brick,
									 forcing.in         = forcing,
									 l.project		    = TRUE,
									 slope.Ta2Tg.in	    = slope.Ta2Tg,
									 intercept.Ta2Tg.in = intercept.Ta2Tg,
									 mod.time			= mod.time,
									 ind.norm.data 		= ind.norm.data,
									 ind.norm.sl 		= ind.norm,
									 luse.brick 	    = luse.brick,
									 i0					= i0)
        slr.gsic <- brick.out[[i]]$gsic.out[i2015:i2065] - brick.out[[i]]$gsic.out[i2015]
        slr.te   <- brick.out[[i]]$te.out[i2015:i2065] - brick.out[[i]]$te.out[i2015]
        slr.gis  <- brick.out[[i]]$simple.out$sle.gis[i2015:i2065] - brick.out[[i]]$simple.out$sle.gis[i2015]
        slr.ais  <- brick.out[[i]]$dais.out$Vais[i2015:i2065] - brick.out[[i]]$dais.out$Vais[i2015]
        lsl.proj[,i] <- fp.gsic.loc*slr.gsic + fp.te.loc  *slr.te  +
                        fp.gis.loc *slr.gis  + fp.ais.loc *slr.ais
        surge.rise[,i] <- surge.factor[i] * lsl.proj[,i]
        setTxtProgressBar(pb, i)
    }
    close(pb)
    print(paste(' ... done.'))

    # make flood risk assessment
    print(paste('Starting ',nr,' flood risk assessments...',sep=''))
    pb <- txtProgressBar(min=0,max=nr,initial=0,style=3)
    for (i in 1:nr) {

        lsl.norm.subs <- subs.rate[i] * nt + lsl.proj[,i] + surge.rise[,i]
        H_eff <- H0 - lsl.norm.subs

        # failure probability is the survival function of the storm surge GEV at effective dike height
        p_fail <- as.numeric(1-pgev(q=1000*H_eff, xi=parameters.gev[i,'shape'], mu=parameters.gev[i,'location'], beta=parameters.gev[i,'scale']))

        # using average annual exceedance probability as the response
        output[i] <- mean(p_fail)

        setTxtProgressBar(pb, i)
    }
    close(pb)
    print(paste(' ... done.'))

    # for Sobol, output must be centered at 0
    output <- output - mean(output)
    return(output)
}
##==============================================================================



##==============================================================================
## Actually run the Sobol
## First-order indices give size of effects on the output response variance
## that are directly (first-order) related to the variable in question.
## Total-order indices give the size of the effects that are related to the
## variable and all interactions between that variable and others, regardless
## of the order.

# get first-order and total-order indices, MCMC-generated
system.time(s.total <- sobolmartinez(model=brick_sobol,
                             parameters.sobol1[1:100,],
                             parameters.sobol2[1:100,],
                             nboot=0))

s.total$T[rev(order(s.total$T[,1])),]

s.total$S[rev(order(s.total$S[,1])),]

# get first- and second-order indices, MCMC-generated
system.time(s.inter <- sobol(model=brick_sobol,
                            parameters.sobol1[1:100,],
                            parameters.sobol2[1:100,],
                            order=2,
                            nboot=0))

print(s.out)
plot(s.out)
##==============================================================================



##==============================================================================
## Write an output file like the modified code from Perry will expect

today=Sys.Date(); today=format(today,format="%d%b%Y")
file.sobolout1 <- paste('../output_calibration/BRICK_Sobol-1-tot_',today,'.txt',sep='')
file.sobolout2 <- paste('../output_calibration/BRICK_Sobol-2_',today,'.txt',sep='')

headers.1st.tot <- matrix(c('Parameter', 'S1', 'S1_conf', 'ST', 'ST_conf'), nrow=1)
output.1st.tot  <- data.frame(cbind( parnames.sobol             ,
                                     s.total$S[,1]              ,
                                     s.total$S[,5]-s.total$S[,1],
                                     s.total$T[,1]              ,
                                     s.total$T[,5]-s.total$T[,1] ))
write.table(headers.1st.tot, file=file.sobolout1, append=FALSE, sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.1st.tot , file=file.sobolout1, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)

headers.2nd     <- matrix(c('Parameter_1', 'Parameter_2', 'S2', 'S2_conf'), nrow=1)
output2.indices <- s.inter$S[(length(parnames.sobol)+1):nrow(s.inter$S),1]
output2.conf    <- s.inter$S[(length(parnames.sobol)+1):nrow(s.inter$S),1] -
                   s.inter$S[(length(parnames.sobol)+1):nrow(s.inter$S),1]

names2  <- rownames(s.inter$S)[(length(parnames.sobol)+1):nrow(s.inter$S)]
names2a <- rep(NA, length(names2))
names2b <- rep(NA, length(names2))
for (i in 1:length(names2)) {
    iast       <- regexpr('[*]', names2[i])[1]
    names2a[i] <- substr(names2[i], start=1, stop=(iast-1))
    names2b[i] <- substr(names2[i], start=(iast+1), stop=nchar(names2[i]))
}

output.2nd <- data.frame(cbind( names2a,
                                names2b,
                                output2.indices,
                                output2.conf ))
write.table(headers.2nd    , file=file.sobolout2, append=FALSE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
write.table(output.2nd     , file=file.sobolout2, append=TRUE , sep = " ",
            quote=FALSE    , row.names = FALSE , col.names=FALSE)
##==============================================================================



##==============================================================================
## Radial plots

## Code history:
# this on by Tony Wong, 6 March 2017, Penn State
###################################
# Adapted from 'radialPlot_vanDantzig.R'
# Originally authored by: Perry Oddo
# Pennsylvania State University
# poddo@psu.edu
###################################
# Adapted from 'radialConvergeTest.R'
# Originally authored by: Calvin Whealton
# Cornell University
# caw324@cornell.edu
####################################
# Code for radial Sobol Analysis plot
# Original code available at:
# https://github.com/calvinwhealton/SensitivityAnalysisPlots
####################################

# Libraries----
library(RColorBrewer) # good color palettes
library(graphics)     # used when plotting polygons
library(plotrix)      # used when plotting circles

# Functions in other files
source('../calibration/BRICK_Sobol_functions.R')

# Set number of parameters being analyzed
n_params <- 39

# Set Sobol indices file name
Sobol_file_1 <- "../output_calibration/BRICK_Sobol-1-tot_06Mar2017.txt"
Sobol_file_2 <- "../output_calibration/BRICK_Sobol-2_06Mar2017.txt"

####################################
# Import data from sensitivity analysis

# First- and total-order indices
s1st <- read.csv(Sobol_file_1,
                  sep=' ',
                  header=TRUE,
                  nrows = n_params,
                  as.is=c(TRUE,rep(FALSE,4)))

# Import second-order indices
s2_table <- read.csv(Sobol_file_2,
               sep=' ',
               header=TRUE,
               nrows = n_params*(n_params-1)/2,
               as.is=c(TRUE,rep(FALSE,3)))

# Convert second-order to upper-triangular matrix
s2 <- matrix(nrow=n_params, ncol = n_params, byrow = FALSE)
s2[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2)
s2 <- as.data.frame(s2)
colnames(s2) <- rownames(s2) <- s1st$Parameter

# Convert confidence intervals to upper-triangular matrix
s2_conf <- matrix(nrow=n_params, ncol = n_params, byrow = FALSE)
s2_conf[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf)
s2_conf <- as.data.frame(s2_conf)
colnames(s2_conf) <- rownames(s2_conf) <- s1st$Parameter

####################################
# Determine which indices are statistically significant
# S1 & ST: using the confidence intervals
#s1st1 <- stat_sig_s1st(s1st
#                      ,method="sig"
#                      ,sigCri='either')

# S1 & ST: using greater than a given value
s1st1 <- stat_sig_s1st(s1st
                      ,method="gtr"
                      ,greater=0.01
                      ,sigCri='either')

# S2: using the confidence intervals
#s2_sig1 <- stat_sig_s2(s2
#                       ,s2_conf
#                       ,method='sig')

# S2: using greater than a given value
s2_sig1 <- stat_sig_s2(s2
                       ,s2_conf
                       ,greater=0.01
                       ,method='gtr')

####################################
# Define groups for the variables and the color schemes
# Defining lists of the variables for each group
#name_list1 <- list('Sea Level' = parnames.sobol[c(ind.brick,ind.subs)],
#                   'Storm Surge' = parnames.sobol[c(ind.surge, ind.gev)],
#                   'Emissions' = parnames.sobol[c(ind.rcp)])
name_list1 <- list('Temperature' = parnames.sobol[c(ind.brick[1:5])],
                   'Sea Level: GSIC' = parnames.sobol[6:9],
                   'Sea Level: TE' = parnames.sobol[10:13],
                   'Sea Level: GIS' = parnames.sobol[14:18],
                   'Sea Level: AIS' = parnames.sobol[19:33],
                   'Storm Surge' = parnames.sobol[c(ind.surge, ind.gev)],
                   'Emissions' = parnames.sobol[c(ind.rcp)])

# add Parameter symbols to plot
name_symbols <- c('S', expression(kappa[D]), expression(alpha[D]),
                  expression('T'[0]), expression('H'[0]), expression(beta[0]),
                  expression('V'['0,GSIC']), 'n', expression('G'['s,0']),
                  expression('a'['TE']), expression('b'['TE']),
                  expression(1/tau['TE']), expression('V'['0,TE']),
                  expression('a'['GIS']), expression('b'['GIS']),
                  expression(alpha['GIS']), expression(beta['GIS']),
                  expression('V'['0,GIS']), expression('a'['ANTO']),
                  expression('b'['ANTO']), expression(gamma), expression(alpha['AIS']),
                  expression(mu), expression(nu), expression('P'[0]),
                  expression(kappa['AIS']), expression('f'[0]),
                  expression('h'[0]), 'c', expression('b'[0]), 'slope',
                  expression(lambda), expression('T'['crit']),
                  expression('C'['surge']), 'subs', expression(mu),
                  expression(sigma), expression(xi), 'RCP')

# Parameter descriptions
#param_desc <- c("Value of goods", "Discount rate", "Construction cost",
#                "Subsidence rate", "Location", "Scale", "Shape",
#                "SLR anomaly", "SLR rate", "Acceleration", "Abrupt SLR rate", "Abrupt SLR year")

source('../Useful/colorblindPalette.R')

colScheme1 <- c(myblue, "red3", mygreen)
colScheme2 <- c(myblue, "red3", mygreen)

# defining list of colors for each group
#col_list1 <- list("Sea Level" = colScheme1[1],
#                  "Storm Surge" = colScheme1[2],
#                  "Emissions"=colScheme1[3])
col_list1 <- list("Temperature"     = rgb(mycol[11,1],mycol[11,2],mycol[11,3]),
                  "Sea Level: GSIC" = rgb(mycol[3,1],mycol[3,2],mycol[3,3]),
                  "Sea Level: TE"   = rgb(mycol[9,1],mycol[9,2],mycol[9,3]),
                  "Sea Level: GIS"  = rgb(mycol[2,1],mycol[2,2],mycol[2,3]),
                  "Sea Level: AIS"  = rgb(mycol[7,1],mycol[7,2],mycol[7,3]),
                  "Storm Surge"     = rgb(mycol[6,1],mycol[6,2],mycol[6,3]),
                  "Emissions"       = rgb(mycol[13,1],mycol[13,2],mycol[13,3]))

# using function to assign variables and colors based on group
s1st1 <- gp_name_col(name_list1
                     ,col_list1
                     ,s1st1)

s1st1$symbols <- name_symbols

#s1st1$desc <- param_desc

#s1st2 <- gp_name_col(name_list2
#                     , col_list2
#                     ,s1st2)

# plotting results
#pdf("Figures/test.pdf")
plotRadCon(df=s1st1
           ,s2=s2
           ,scaling = .45
           ,s2_sig=s2_sig1
           ,filename = './sobol_fig_test'
           ,plotType = 'EPS'
           ,legLoc = "bottomcenter",cex = .76
)

#dev.off()
##==============================================================================

#Not yet modified below here


##==============================================================================

## Note on bandwidths and number of nodes for Kernel density estimate fitting:
##    With hundreds of thousands of parameters in the distributions we seek to
##    fit KDEs to, the fits and bandwidths are fairly insensitive to the number
##    of KDE nodes you plop down (n.node, above). With 300,000 = n.parameters,
##    Tony tested n.node = 20, 50, 100, 200 and 1000. The bandwidths for
##    n.node >= 50 are all the same to within 7 significant figures.

## Write a CSV file with the successful parameter combinations and bandwidths
## Structure of the CSV file is as follows. 5 columns, 2*n.sample rows (2 chains
## with n.samples each).
##    First row: Parameter names.
##    Rows 2-??: The calibrated parameter values.
##    Last row:  The bandwidths. These are the standard
##               deviations of the normal distributions one should sample from.
##               The idea is that if you pick a row out of this CSV file, and
##               draw random-normally with these bandwidths (stdevs) around each
##               parameter value, you are sampling from the joint distribution.
bandwidths=rep(NA,n.parameters)
for (i in 1:n.parameters){
  bandwidths[i]=pdf.all[[i]]$bw
}

## Plot the distributions
par(mfrow=c(4,4))
for (p in 1:n.parameters) {
  plot(pdf.all[[p]]$x,pdf.all[[p]]$y,type='l',xlab=parnames[p],ylab='density',main="")
}

## Write the calibrated parameters file (csv version - antiquated)
#to.file = rbind(parameters.posterior,bandwidths)  # bandwidths are in the last row
#rownames(to.file)=NULL
#colnames(to.file)=parnames
#today=Sys.Date(); today=format(today,format="%d%b%Y")
#filename=paste('../output_calibration/DAIS_calibratedParameters_',today,'.csv', sep="")
#write.table(to.file, file=filename, sep=",", qmethod="double", row.names=FALSE)

# Write the calibrated parameters file (netCDF version - better)

## Get maximum length of parameter name, for width of array to write to netcdf
## this code will write an n.parameters (rows) x n.ensemble (columns) netcdf file
## to get back into the shape BRICK expects, just transpose it
lmax=0
for (i in 1:length(parnames)){lmax=max(lmax,nchar(parnames[i]))}

library(ncdf4)
dim.parameters <- ncdim_def('n.parameters', '', 1:ncol(parameters.posterior), unlim=FALSE)
dim.name <- ncdim_def('name.len', '', 1:lmax, unlim=FALSE)
dim.ensemble <- ncdim_def('n.ensemble', 'ensemble member', 1:nrow(parameters.posterior), unlim=TRUE)
parameters.var <- ncvar_def('DAIS_parameters', '', list(dim.parameters,dim.ensemble), -999)
parnames.var <- ncvar_def('parnames', '', list(dim.name,dim.parameters), prec='char')
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.daisparameters = paste('../output_calibration/DAIS_calibratedParameters_',today,'.nc',sep="")
outnc <- nc_create(filename.daisparameters, list(parameters.var,parnames.var))
ncvar_put(outnc, parameters.var, t(parameters.posterior))
ncvar_put(outnc, parnames.var, parnames)
nc_close(outnc)

##==============================================================================
## End
##==============================================================================


## Extra code snippets


##==============================================================================
## Latin Hypercube to get parameters to send into Sobol
## Could fit to distributions from Wong et al 2017 instead of the uniforms.

require(lhs)

# Draw LHS samples (need two)
# Okay to leave on U[0,1] because brick_sobol function will scale up to the
# parameter ranges.
n.lhs = 1e2
parameters.lhs1 <- randomLHS(n.lhs, length(parnames))
parameters.lhs2 <- randomLHS(n.lhs, length(parnames))
colnames(parameters.lhs1) <- parnames
colnames(parameters.lhs2) <- parnames
parameters.lhs1 <- data.frame(parameters.lhs1)
parameters.lhs2 <- data.frame(parameters.lhs2)
##==============================================================================



##==============================================================================
## Map to [0,1]
parameters.sample1 <- parameters.sample.orig1
parameters.sample2 <- parameters.sample.orig2
for (j in 1:ncol(parameters.sample)) {
    parameters.sample1[,j] <- map.range(parameters.sample.orig1[,j], bound.lower[j], bound.upper[j], 0, 1)
    parameters.sample2[,j] <- map.range(parameters.sample.orig2[,j], bound.lower[j], bound.upper[j], 0, 1)
}
##==============================================================================



##==============================================================================
# Some of the parameter combinations will yield complete melt, and NaNs.
# Filter these out.
dais.lhs1 <- dais_sobol(parameters.lhs01.1)
dais.lhs2 <- dais_sobol(parameters.lhs01.2)

ibad1 <- which(is.na(dais.lhs1))
ibad2 <- which(is.na(dais.lhs2))

parameters.lhs01.1 <- parameters.lhs01.1[-ibad1,]
parameters.lhs01.2 <- parameters.lhs01.2[-ibad2,]

nmax <- min(nrow(parameters.lhs01.1), nrow(parameters.lhs01.2))
parameters.lhs01.1 <- parameters.lhs01.1[1:nmax,]
parameters.lhs01.2 <- parameters.lhs01.2[1:nmax,]
##==============================================================================



## Check the samples
par(mfrow=c(2,2))
pp <- 12
hist(parameters.brick[,pp], xlab=parnames.brick[pp], xlim=quantile(parameters.brick[,pp], c(0,1)))
plot(kde.brick[[pp]]$x, kde.brick[[pp]]$y, type='l', xlim=quantile(parameters.brick[,pp], c(0,1)))
hist(map.range(parameters.sample1.brick[,pp], lbin=0, ubin=1, lbout=bound.lower.brick[pp], ubout=bound.upper.brick[pp]), xlim=quantile(parameters.brick[,pp], c(0,1)))
hist(map.range(parameters.sample1.brick[,pp], lbin=0, ubin=1, lbout=min(parameters.brick[pp]), ubout=max(parameters.brick[pp])), xlim=quantile(parameters.brick[,pp], c(0,1)))
