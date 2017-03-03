##==============================================================================
## Script to calibrate DAIS model for use within BRICK coupled modeling
## framework.
##
## Parameters to calibrate, and standard values (Shaffer, 2014):
## [1] gamma = 1/2-17/4  # sensitivity of ice flow to sea level
## [2] alpha = 0-1       # sensitivity of ice flow to ocean subsurface temperature
## [3] mu = 8.7          # Profile parameter related to ice stress [m^(1/2)]
## [4] nu = 1.2e-2       # Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
## [5] P0 = 0.35         # Precipitation at 0C [m of ice/yr]
## [6] kappa = 4e-2      # Relates precipitation to temperature [K^-1]
## [7] f0 = 1.2          # Constant of proportionality for ice speed [m/yr]
## [8] h0 = 1471         # Initial value for runoff line calculation [m]
## [9] c = 95            # Second value for runoff line calculation [m]
## [10] b0 = 775         # Height of bed at the center of the continent [m]
## [11] slope = 6e-4     # Slope of the bed
##
## Much of the original code is thanks to Kelsey Ruckert <klr324@psu.edu>
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

rm(list =ls()) #Clear global environment

## Setup packages and libraries
#install.packages('compiler')
#install.packages('pscl')
library(pscl) # install inverse gamma distribution
library(compiler)
enableJIT(3)
enableJIT(3)

## Set the seed (for reproducibility)
set.seed(1234)

## Read the data forcing for hindcasts and projections. Yields:
##  Ta (Antarctic temperature reduced to sea-level)
##  Toc (High latitude subsurface ocean temperature)
##  SL, obs.sl (Reconstructed sea level, Church and White (2011) modern-era sea level)
##  dSL (time rate of change of sea level)
##  date (240 Kyr before present to 2100 AD at 1-year intervals of forcings)
l.project=FALSE
date = seq(-239999,300,1) #240 Kyr BP to 2100AD at one year intervals of the forcings
dSL = scan("../data/future_GSL.txt")   #Time rate of change of sea-level
Ta  = scan("../data/future_TA.txt")     #Antarctic temperature reduced to sea-level
Toc = scan("../data/future_TO.txt")     #High latitude subsurface ocean temperature
SL  = scan("../data/future_SL.txt")     #Reconstructed sea-level

## Set Best Case (Case #4) from Shaffer (2014), and parameter ranges
    # >>> with anto? <<<
parnames.dais    = c('anto.a','anto.b','gamma','alpha.dais','mu'  ,'nu'  ,'P0' ,'kappa.dais','f0' ,'h0'  ,'c'  , 'b0','slope' )
parameters0.dais = c(  0.1574, 0.6677 ,  2    , 0.35       , 8.7  , 0.012, 0.35, 0.04       , 1.2 , 1471 , 95  , 775 , 0.0006 )
bound.lower.dais = c( 0.0    , 0      ,  0.5  , 0          , 7.05 , 0.003,0.026, 0.025      , 0.6 , 735.5, 47.5, 740 , 0.00045)
bound.upper.dais = c( 1.0    , 2      ,  4.25 , 1          , 13.65, 0.015, 1.5 , 0.085      , 1.8 ,2206.5,142.5, 820 , 0.00075)
# Test narrower widths for h0 and c
bound.lower.dais = c( 0.0    , 0      ,  0.5  , 0          , 7.05 , 0.003,0.026, 0.025      , 0.6 , 1271 , 75  , 740 , 0.00045)
bound.upper.dais = c( 1.0    , 2      ,  4.25 , 1          , 13.65, 0.015, 1.5 , 0.085      , 1.8 , 1671 , 115 , 820 , 0.00075)

alpha.var = 2     # alpha parameter for inverse gamma for var (E[x]=beta/(alpha+1))
beta.var = 1      # beta parameter for inverse gamma for var (uncertainty parameter)
                  # note that the upper bound on var.dais is not actually imposed; just for plotting

## Source the DAIS model
source('../fortran/R/daisantoF.R')
##==============================================================================



##==============================================================================
## Need preliminary function to map ranges from [0,1] and back
map.range <- function(X, lbin, ubin, lbout, ubout){
    Y <- lbout + (ubout-lbout)*( (X-lbin)/(ubin-lbin) )
    return(Y)
}
##==============================================================================



##==============================================================================
## Latin Hypercube to get parameters to send into Sobol
## Could fit to distributions from Ruckert et al 2017 or Wong et al 2017
## instead of the uniforms...

require(lhs)

# Draw LHS samples (need two)
# Okay to leave on U[0,1] because dais_sobol function will scale up to the
# parameter ranges.
n.lhs = 1e5
parameters.lhs01.1 <- randomLHS(n.lhs, length(parnames.dais))
parameters.lhs01.2 <- randomLHS(n.lhs, length(parnames.dais))
colnames(parameters.lhs01.1) <- parnames.dais
colnames(parameters.lhs01.2) <- parnames.dais
parameters.lhs01.1 <- data.frame(parameters.lhs01.1)
parameters.lhs01.2 <- data.frame(parameters.lhs01.2)
##==============================================================================



##==============================================================================
## Or get parameters from Ruckert et al 2017?
load('~/Downloads/DAIS_calib_MCMC_C1234_relative_8e5.RData')
dais.mcmc <- DAIS_chains[ceiling(nrow(DAIS_chains)*0.5):nrow(DAIS_chains),]

## Or get parameters from the DAIS standalone calibration in Wong et al 2017?
library(ncdf4)
filename.DAIScalibration = "~/codes/BRICK_Experiments_backup/output_calibration/DAIS_calibratedParameters_11Aug2016.nc"
ncdata <- nc_open(filename.DAIScalibration)
dais.mcmc <- t(ncvar_get(ncdata, 'DAIS_parameters'))
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)

## Or get parameters from the DAIS within BRICK calibration?
## Overwrite from the above, since you need to get "parnames"
## Also, need to cut off a burn-in from these.
library(ncdf4)
filename.DAIScalibration = "~/codes/BRICK/output_calibration/BRICK_calibDAISonly_01Mar2017.nc"
ncdata <- nc_open(filename.DAIScalibration)
dais.mcmc <- t(ncvar_get(ncdata, 'DAIS_parameters'))
nc_close(ncdata)
dais.mcmc <- dais.mcmc[7e5:nrow(dais.mcmc),]

colnames(dais.mcmc) <- parnames
ind.sample <- sample(x=1:nrow(dais.mcmc), size=1e5, replace=FALSE)
parameters.sample <- data.frame(dais.mcmc[ind.sample,])
#parameters.sample <- parameters.sample[,-c(12,13)]      # remove satatistical parmetesr
parameters.sample <- parameters.sample[,-c(14)]      # remove satatistical parmetesr
bound.lower.dais <- bound.lower.dais[-c(14)]
bound.upper.dais <- bound.upper.dais[-c(14)]

## Map to [0,1]
parameters.sample01 <- parameters.sample
for (j in 1:ncol(parameters.sample)) {
    parameters.sample01[,j] <- map.range(parameters.sample[,j], bound.lower.dais[j], bound.upper.dais[j], 0, 1)
}

##==============================================================================




##==============================================================================
## Define model for Sobol
## Needs to take in a data frame, with each column as a different model
## parameter, and return an output data frame with same number of rows (one
## column response vector).

# reponse in 2000, and only run 1850-2000, because that's what we will run
# BRICK for
i2000 <- which(date==0)
i1850 <- which(date==-150)

# using mapply is generally faster than a 'for' loop
dais_sobol <- function(dataframe.in){
    nr <- nrow(dataframe.in)
    output <- rep(0,nr)
    input  <- dataframe.in
    input[,1] <- map.range(dataframe.in[,1], 0, 1, bound.lower.dais[1], bound.upper.dais[1])
    input[,2] <- map.range(dataframe.in[,2], 0, 1, bound.lower.dais[2], bound.upper.dais[2])
    input[,3] <- map.range(dataframe.in[,3], 0, 1, bound.lower.dais[3], bound.upper.dais[3])
    input[,4] <- map.range(dataframe.in[,4], 0, 1, bound.lower.dais[4], bound.upper.dais[4])
    input[,5] <- map.range(dataframe.in[,5], 0, 1, bound.lower.dais[5], bound.upper.dais[5])
    input[,6] <- map.range(dataframe.in[,6], 0, 1, bound.lower.dais[6], bound.upper.dais[6])
    input[,7] <- map.range(dataframe.in[,7], 0, 1, bound.lower.dais[7], bound.upper.dais[7])
    input[,8] <- map.range(dataframe.in[,8], 0, 1, bound.lower.dais[8], bound.upper.dais[8])
    input[,9] <- map.range(dataframe.in[,9], 0, 1, bound.lower.dais[9], bound.upper.dais[9])
    input[,10] <- map.range(dataframe.in[,10], 0, 1, bound.lower.dais[10], bound.upper.dais[10])
    input[,11] <- map.range(dataframe.in[,11], 0, 1, bound.lower.dais[11], bound.upper.dais[11])
    input[,12] <- map.range(dataframe.in[,12], 0, 1, bound.lower.dais[12], bound.upper.dais[12])
    input[,13] <- map.range(dataframe.in[,13], 0, 1, bound.lower.dais[13], bound.upper.dais[13])
    dataframe.in <- input
    AIS_melt <- mapply(daisantoF,
                        anto.a=dataframe.in[,1],
                        anto.b=dataframe.in[,2],
                        gamma=dataframe.in[,3],
                        alpha=dataframe.in[,4],
                        mu=dataframe.in[,5],
                        nu    = dataframe.in[,6],
                        P0    = dataframe.in[,7],
                        kappa = dataframe.in[,8],
                        f0    = dataframe.in[,9],
                        h0    = dataframe.in[,10],
                        c     = dataframe.in[,11],
                        b0    = dataframe.in[,12],
                        slope = dataframe.in[,13],
                        MoreArgs=list(Tg    = Tg.recon[i1850:i2000],
                        slope.Ta2Tg = slope.Ta2Tg,
                        intercept.Ta2Tg = intercept.Ta2Tg,
                        SL    = SL[i1850:i2000],
                        dSL   = dSL[i1850:i2000],
                        includes_dSLais=1))
    output <- AIS_melt[nrow(AIS_melt),]
    output <- output - mean(output)
    return(output)
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




##==============================================================================
## Actually run the Sobol
## First-order indices give size of effects on the output response variance
## that are directly (first-order) related to the variable in question.
## Total-order indices give the size of the effects that are related to the
## variable and all interactions between that variable and others, regardless
## of the order.

# get first-order and total-order indices, MCMC-generated
system.time(s.total.mcmc <- sobolmartinez(model=dais_sobol,
                             parameters.sample01[1:10000,],
                             parameters.sample01[10001:20000,],
                             nboot=100))
# compare to LHS-generated
system.time(s.total.lhs <- sobolmartinez(model=dais_sobol,
                                X1=parameters.lhs01.1[1:1000,],
                                X2=parameters.lhs01.2[1:1000,],
                                nboot=100))

par(mfrow=c(2,1))
plot(s.total.lhs); lines(c(-100,100),c(0,0), col='red')
plot(s.total.mcmc); lines(c(-100,100),c(0,0), col='red')

# get first- and second-order indices, MCMC-generated
system.time(s.inter.mcmc <- sobol(model=dais_sobol,
                            parameters.sample01[1:50000,],
                            parameters.sample01[50001:100000,],
                            order=2,
                            nboot=1000))
# compare to LHS-generated
system.time(s.inter.lhs <- sobol(model=dais_sobol,
                            parameters.lhs01.1[1:50000,],
                            parameters.lhs01.2[1:50000,],
                            order=2,
                            nboot=1000))


print(s.out)
plot(s.out)

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
