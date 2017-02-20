##==============================================================================
## Plots and tables for BRICK scenarios paper (...)
##
## These routines also are provided as guidance for users to develop their own
## plotting and analysis routines.
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

rm(list=ls())

## Initial set-up
library(ncdf4)

## File name for the BRICK physical model output (netCDF4)
#filename.brick.nofd = '../output_model/BRICK-model_physical_fd-gamma_08Dec2016.nc'#no-FD case uses gamma, without disintegration
#filename.brick.uniform = '../output_model/BRICK-model_physical_fd-uniform_08Dec2016.nc'
#filename.brick.gamma = '../output_model/BRICK-model_physical_fd-gamma_08Dec2016.nc'
filename.brick.allslr = '../output_model/BRICK-model_physical_allslr_20Feb2017.nc'

## File name for the Van Dantzig model output (netCDF4)
## Each of these also has x3 RCP scenarios, x2 storm surge scenarios
filename.vandantzig.nofd = '../output_model/VanDantzig_fd-none_2065_07Feb2017.nc'
filename.vandantzig.uniform = '../output_model/VanDantzig_fd-uniform_2065_07Feb2017.nc'
filename.vandantzig.gamma = '../output_model/VanDantzig_fd-gamma_2065_07Feb2017.nc'

## File name for the BRICK post-calibrated parameters (netcdf) (the BRICK output came from these guys)
filename.parameters.nofd = '../output_calibration/BRICK-model_postcalibratedParameters_fd-none_08Dec2016.nc'
filename.parameters.uniform = '../output_calibration/BRICK-model_postcalibratedParameters_fd-uniform_08Dec2016.nc'
filename.parameters.gamma = '../output_calibration/BRICK-model_drawcalibratedParameters_fd-gamma_08Dec2016.nc'

## Other files
filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_01Nov2016.csv"
filename.fingerprints = "../fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"

## Other useful scripts
source('../Useful/colorblindPalette.R') # Get nice plotting colors: mycol array
source('../Useful/MultipleOutput.R')    # defines the useful ":=" operator

## set up the RCP forcings' default colors within mycol
c85=6;
c45=4;
c26=2;

## And set the IPCC RCP colors
col26 <- c(0, 0, 255)/255
col45 <- c(121, 188, 255)/255
col60 <- c(255, 130, 45)/255
col85 <- c(255, 0, 0)/255

## Where would you like to save the plots?
plotdir='~/Box\ Sync/Wong-Projects/BRICK_scenarios/figures/'

## Get nice plotting colors: mycol array
source('../Useful/colorblindPalette.R')

##==============================================================================






##==============================================================================
##==============================================================================
## Grab the van Dantzig output and store, for all scenarios
## Also store the local sea level, and the surge.factor parameters for each AIS
## scenario.

# first level are RCP scenarios
gev.names <- c('location','scale','shape')
scen.rcp = c('rcp26','rcp45','rcp85')
cost = vector("list",3); names(cost)=scen.rcp
loss = vector("list",3); names(loss)=scen.rcp
investment = vector("list",3); names(investment)=scen.rcp
preturn = vector("list",3); names(preturn)=scen.rcp
lsl = vector("list",3); names(lsl)=scen.rcp
init = vector("list",3); names(init)=scen.rcp

# second level are AIS fast dynamics scenarios
scen.ais = c('none','gamma','uniform')
scen.ss = c('st','ns')
for (rcp in scen.rcp) {
  cost[[rcp]] = vector("list",3); names(cost[[rcp]])=scen.ais
  loss[[rcp]] = vector("list",3); names(loss[[rcp]])=scen.ais
  investment[[rcp]] = vector("list",3); names(investment[[rcp]])=scen.ais
  preturn[[rcp]] = vector("list",3); names(preturn[[rcp]])=scen.ais
  lsl[[rcp]] = vector("list",3); names(lsl[[rcp]])=scen.ais
  init[[rcp]] = vector("list",3); names(init[[rcp]])=scen.ais

# third level are storm surge scenarios
  for (ais in scen.ais) {
    cost[[rcp]][[ais]] = vector("list",2); names(cost[[rcp]][[ais]])=scen.ss
    loss[[rcp]][[ais]] = vector("list",2); names(loss[[rcp]][[ais]])=scen.ss
    investment[[rcp]][[ais]] = vector("list",2); names(investment[[rcp]][[ais]])=scen.ss
    preturn[[rcp]][[ais]] = vector("list",2); names(preturn[[rcp]][[ais]])=scen.ss
  }
}

## What is the economically efficient heightening for each ensemble member?
## (and the index, so we can grab the return period)

iopt = cost     # optimal heightening index
ipre = cost     # prescribed return period index
Ropt = cost     # optimal return period
Hopt = cost     # optimal heightening

# read no fast dyanmics results
ncdata <- nc_open(filename.vandantzig.nofd)
  heightening <- ncvar_get(ncdata, 'H')
  surge.factor <- ncvar_get(ncdata, 'surge_factor') # all RCPs have same surge factor
  gev.stat <- ncvar_get(ncdata, 'gev_stat')
  VD.params <- ncvar_get(ncdata, 'VD_params')
  colnames(gev.stat) <- gev.names

  cost$rcp26$none$st<- ncvar_get(ncdata, 'ExpectedCost_RCP26_stat')
  loss$rcp26$none$st<- ncvar_get(ncdata, 'ExpectedLoss_RCP26_stat')
  investment$rcp26$none$st<- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_stat')
  preturn$rcp26$none$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat')

  cost$rcp45$none$st<- ncvar_get(ncdata, 'ExpectedCost_RCP45_stat')
  loss$rcp45$none$st<- ncvar_get(ncdata, 'ExpectedLoss_RCP45_stat')
  investment$rcp45$none$st<- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_stat')
  preturn$rcp45$none$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat')

  cost$rcp85$none$st<- ncvar_get(ncdata, 'ExpectedCost_RCP85_stat')
  loss$rcp85$none$st<- ncvar_get(ncdata, 'ExpectedLoss_RCP85_stat')
  investment$rcp85$none$st<- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_stat')
  preturn$rcp85$none$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat')

  cost$rcp26$none$ns<- ncvar_get(ncdata, 'ExpectedCost_RCP26_nonstat')
  loss$rcp26$none$ns<- ncvar_get(ncdata, 'ExpectedLoss_RCP26_nonstat')
  investment$rcp26$none$ns<- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_nonstat')
  preturn$rcp26$none$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat')

  cost$rcp45$none$ns<- ncvar_get(ncdata, 'ExpectedCost_RCP45_nonstat')
  loss$rcp45$none$ns<- ncvar_get(ncdata, 'ExpectedLoss_RCP45_nonstat')
  investment$rcp45$none$ns<- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_nonstat')
  preturn$rcp45$none$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat')

  cost$rcp85$none$ns<- ncvar_get(ncdata, 'ExpectedCost_RCP85_nonstat')
  loss$rcp85$none$ns<- ncvar_get(ncdata, 'ExpectedLoss_RCP85_nonstat')
  investment$rcp85$none$ns<- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_nonstat')
  preturn$rcp85$none$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat')
nc_close(ncdata)

# read gamma fast dynamics prior results
ncdata <- nc_open(filename.vandantzig.gamma)
  surge.factor <- ncvar_get(ncdata, 'surge_factor') # all RCPs have same surge factor
  gev.stat <- ncvar_get(ncdata, 'gev_stat')
  colnames(gev.stat) <- gev.names

  cost$rcp26$gamma$st<- ncvar_get(ncdata, 'ExpectedCost_RCP26_stat')
  loss$rcp26$gamma$st<- ncvar_get(ncdata, 'ExpectedLoss_RCP26_stat')
  investment$rcp26$gamma$st<- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_stat')
  preturn$rcp26$gamma$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat')

  cost$rcp45$gamma$st<- ncvar_get(ncdata, 'ExpectedCost_RCP45_stat')
  loss$rcp45$gamma$st<- ncvar_get(ncdata, 'ExpectedLoss_RCP45_stat')
  investment$rcp45$gamma$st<- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_stat')
  preturn$rcp45$gamma$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat')

  cost$rcp85$gamma$st<- ncvar_get(ncdata, 'ExpectedCost_RCP85_stat')
  loss$rcp85$gamma$st<- ncvar_get(ncdata, 'ExpectedLoss_RCP85_stat')
  investment$rcp85$gamma$st<- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_stat')
  preturn$rcp85$gamma$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat')

  cost$rcp26$gamma$ns<- ncvar_get(ncdata, 'ExpectedCost_RCP26_nonstat')
  loss$rcp26$gamma$ns<- ncvar_get(ncdata, 'ExpectedLoss_RCP26_nonstat')
  investment$rcp26$gamma$ns<- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_nonstat')
  preturn$rcp26$gamma$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat')

  cost$rcp45$gamma$ns<- ncvar_get(ncdata, 'ExpectedCost_RCP45_nonstat')
  loss$rcp45$gamma$ns<- ncvar_get(ncdata, 'ExpectedLoss_RCP45_nonstat')
  investment$rcp45$gamma$ns<- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_nonstat')
  preturn$rcp45$gamma$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat')

  cost$rcp85$gamma$ns<- ncvar_get(ncdata, 'ExpectedCost_RCP85_nonstat')
  loss$rcp85$gamma$ns<- ncvar_get(ncdata, 'ExpectedLoss_RCP85_nonstat')
  investment$rcp85$gamma$ns<- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_nonstat')
  preturn$rcp85$gamma$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat')
nc_close(ncdata)

# read uniform fast dynamics prior results
ncdata <- nc_open(filename.vandantzig.uniform)
  surge.factor <- ncvar_get(ncdata, 'surge_factor') # all RCPs have same surge factor
  gev.stat <- ncvar_get(ncdata, 'gev_stat')
  colnames(gev.stat) <- gev.names

  cost$rcp26$uniform$st<- ncvar_get(ncdata, 'ExpectedCost_RCP26_stat')
  loss$rcp26$uniform$st<- ncvar_get(ncdata, 'ExpectedLoss_RCP26_stat')
  investment$rcp26$uniform$st<- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_stat')
  preturn$rcp26$uniform$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat')

  cost$rcp45$uniform$st<- ncvar_get(ncdata, 'ExpectedCost_RCP45_stat')
  loss$rcp45$uniform$st<- ncvar_get(ncdata, 'ExpectedLoss_RCP45_stat')
  investment$rcp45$uniform$st<- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_stat')
  preturn$rcp45$uniform$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat')

  cost$rcp85$uniform$st<- ncvar_get(ncdata, 'ExpectedCost_RCP85_stat')
  loss$rcp85$uniform$st<- ncvar_get(ncdata, 'ExpectedLoss_RCP85_stat')
  investment$rcp85$uniform$st<- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_stat')
  preturn$rcp85$uniform$st<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat')

  cost$rcp26$uniform$ns<- ncvar_get(ncdata, 'ExpectedCost_RCP26_nonstat')
  loss$rcp26$uniform$ns<- ncvar_get(ncdata, 'ExpectedLoss_RCP26_nonstat')
  investment$rcp26$uniform$ns<- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_nonstat')
  preturn$rcp26$uniform$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat')

  cost$rcp45$uniform$ns<- ncvar_get(ncdata, 'ExpectedCost_RCP45_nonstat')
  loss$rcp45$uniform$ns<- ncvar_get(ncdata, 'ExpectedLoss_RCP45_nonstat')
  investment$rcp45$uniform$ns<- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_nonstat')
  preturn$rcp45$uniform$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat')

  cost$rcp85$uniform$ns<- ncvar_get(ncdata, 'ExpectedCost_RCP85_nonstat')
  loss$rcp85$uniform$ns<- ncvar_get(ncdata, 'ExpectedLoss_RCP85_nonstat')
  investment$rcp85$uniform$ns<- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_nonstat')
  preturn$rcp85$uniform$ns<- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat')
nc_close(ncdata)

n.ensemble <- c(ncol(cost$rcp85$none$ns) , ncol(cost$rcp85$gamma$ns) , ncol(cost$rcp85$uniform$ns) )
names(n.ensemble) <- scen.ais
n.height <- length(heightening)

## Initialize arrays for average cost, regret, return period among states of world
## within a scenario
Ctmp <- rep(NA,n.height)
Rtmp <- rep(NA,n.height)
Vtmp <- rep(NA,n.height)

# And initialize for Laplace-averaged (all scenarios) Cost
Cavg <- rep(0,n.height)    # initialize with 0, otherwise recurions isn't general
Ravg <- rep(0,n.height)
Vavg <- rep(0,n.height)
CVARavg <- rep(0,n.height)
navg <- 0          # store the total number of points in the average so far, to use "recursion"

## What is the economically efficient heightening for each scenario?
## (and the index, so we can grab the return period)
ieum <- vector("list",3); names(ieum)<-scen.rcp     # EUM-optimal heightening index for scenario
Ceum <- vector("list",3); names(Ceum)<-scen.rcp     # EUM-optimal cost for scenario
Reum <- vector("list",3); names(Reum)<-scen.rcp     # EUM-optimal return period scenario
Heum <- vector("list",3); names(Heum)<-scen.rcp     # EUM-optimal heightening for scenario
Hpre <- vector("list",3); names(Hpre)<-scen.rcp     # heightening for 1:100 reliability

ireg <- vector("list",3); names(ireg)<-scen.rcp     # min regret heightening index for scenario
Creg <- vector("list",3); names(Creg)<-scen.rcp     # min regret cost for scenario
Rreg <- vector("list",3); names(Rreg)<-scen.rcp     # min regret return period scenario
Hreg <- vector("list",3); names(Hreg)<-scen.rcp     # min regret heightening for scenario

icva <- vector("list",3); names(icva)<-scen.rcp     # min cvar heightening index for scenario
Ccva <- vector("list",3); names(Ccva)<-scen.rcp     # min cvar cost for scenario
Rcva <- vector("list",3); names(Rcva)<-scen.rcp     # min cvar return period scenario
Hcva <- vector("list",3); names(Hcva)<-scen.rcp     # min cvar heightening for scenario
alpha <- 0.05 # tail for cvar

Cens.eum <- vector("list",3); names(Cens.eum)<-scen.rcp     # each SOW expected cost, at EUM-optimal heightening
Rens.eum <- vector("list",3); names(Rens.eum)<-scen.rcp     # each SOW return period, at EUM-optimal heightening

iopt <- vector("list",3); names(iopt)<-scen.rcp     # optimal heightening index for SOW
Copt <- vector("list",3); names(Copt)<-scen.rcp     # optimal cost for SOW
Ropt <- vector("list",3); names(Ropt)<-scen.rcp     # optimal return period SOW
Hopt <- vector("list",3); names(Hopt)<-scen.rcp     # optimal heightening for SOW

C.h <- vector("list",3); names(C.h)<-scen.rcp     # average cost at h heightening (each scenario)
V.h <- vector("list",3); names(V.h)<-scen.rcp     # average regret at h heightening (each scenario)
R.h <- vector("list",3); names(R.h)<-scen.rcp     # average return period at h heightening (each scenario)

for (rcp in scen.rcp) {
  ieum[[rcp]] <- vector("list",3); names(ieum[[rcp]])<-scen.ais
  Ceum[[rcp]] <- vector("list",3); names(Ceum[[rcp]])<-scen.ais
  Reum[[rcp]] <- vector("list",3); names(Reum[[rcp]])<-scen.ais
  Heum[[rcp]] <- vector("list",3); names(Heum[[rcp]])<-scen.ais
  Hpre[[rcp]] <- vector("list",3); names(Hpre[[rcp]])<-scen.ais
  ireg[[rcp]] <- vector("list",3); names(ireg[[rcp]])<-scen.ais
  Creg[[rcp]] <- vector("list",3); names(Creg[[rcp]])<-scen.ais
  Rreg[[rcp]] <- vector("list",3); names(Rreg[[rcp]])<-scen.ais
  Hreg[[rcp]] <- vector("list",3); names(Hreg[[rcp]])<-scen.ais
  icva[[rcp]] <- vector("list",3); names(icva[[rcp]])<-scen.ais
  Ccva[[rcp]] <- vector("list",3); names(Ccva[[rcp]])<-scen.ais
  Rcva[[rcp]] <- vector("list",3); names(Rcva[[rcp]])<-scen.ais
  Hcva[[rcp]] <- vector("list",3); names(Hcva[[rcp]])<-scen.ais
  Cens.eum[[rcp]] <- vector("list",3); names(Cens.eum[[rcp]])<-scen.ais
  Rens.eum[[rcp]] <- vector("list",3); names(Rens.eum[[rcp]])<-scen.ais
  iopt[[rcp]] <- vector("list",3); names(iopt[[rcp]])<-scen.ais
  Copt[[rcp]] <- vector("list",3); names(Copt[[rcp]])<-scen.ais
  Ropt[[rcp]] <- vector("list",3); names(Ropt[[rcp]])<-scen.ais
  Hopt[[rcp]] <- vector("list",3); names(Hopt[[rcp]])<-scen.ais
  V.h[[rcp]] <- vector("list",3); names(V.h[[rcp]])<-scen.ais
  C.h[[rcp]] <- vector("list",3); names(C.h[[rcp]])<-scen.ais
  R.h[[rcp]] <- vector("list",3); names(R.h[[rcp]])<-scen.ais
  for (ais in scen.ais) {
    ieum[[rcp]][[ais]] <- vector("list",2); names(ieum[[rcp]][[ais]])<-scen.ss
    Ceum[[rcp]][[ais]] <- vector("list",2); names(Ceum[[rcp]][[ais]])<-scen.ss
    Reum[[rcp]][[ais]] <- vector("list",2); names(Reum[[rcp]][[ais]])<-scen.ss
    Heum[[rcp]][[ais]] <- vector("list",2); names(Heum[[rcp]][[ais]])<-scen.ss
    Hpre[[rcp]][[ais]] <- vector("list",2); names(Hpre[[rcp]][[ais]])<-scen.ss
    ireg[[rcp]][[ais]] <- vector("list",2); names(ireg[[rcp]][[ais]])<-scen.ss
    Creg[[rcp]][[ais]] <- vector("list",2); names(Creg[[rcp]][[ais]])<-scen.ss
    Rreg[[rcp]][[ais]] <- vector("list",2); names(Rreg[[rcp]][[ais]])<-scen.ss
    Hreg[[rcp]][[ais]] <- vector("list",2); names(Hreg[[rcp]][[ais]])<-scen.ss
    icva[[rcp]][[ais]] <- vector("list",2); names(icva[[rcp]][[ais]])<-scen.ss
    Ccva[[rcp]][[ais]] <- vector("list",2); names(Ccva[[rcp]][[ais]])<-scen.ss
    Rcva[[rcp]][[ais]] <- vector("list",2); names(Rcva[[rcp]][[ais]])<-scen.ss
    Hcva[[rcp]][[ais]] <- vector("list",2); names(Hcva[[rcp]][[ais]])<-scen.ss
    Cens.eum[[rcp]][[ais]] <- vector("list",2); names(Cens.eum[[rcp]][[ais]])<-scen.ss
    Rens.eum[[rcp]][[ais]] <- vector("list",2); names(Rens.eum[[rcp]][[ais]])<-scen.ss
    iopt[[rcp]][[ais]] <- vector("list",2); names(iopt[[rcp]][[ais]])<-scen.ss
    Copt[[rcp]][[ais]] <- vector("list",2); names(Copt[[rcp]][[ais]])<-scen.ss
    Ropt[[rcp]][[ais]] <- vector("list",2); names(Ropt[[rcp]][[ais]])<-scen.ss
    Hopt[[rcp]][[ais]] <- vector("list",2); names(Hopt[[rcp]][[ais]])<-scen.ss
    V.h[[rcp]][[ais]] <- vector("list",2); names(V.h[[rcp]][[ais]])<-scen.ss
    C.h[[rcp]][[ais]] <- vector("list",2); names(C.h[[rcp]][[ais]])<-scen.ss
    R.h[[rcp]][[ais]] <- vector("list",2); names(R.h[[rcp]][[ais]])<-scen.ss
    for (ss in scen.ss) {
      ieum[[rcp]][[ais]][[ss]] <- NA
      Ceum[[rcp]][[ais]][[ss]] <- NA
      Reum[[rcp]][[ais]][[ss]] <- NA
      Heum[[rcp]][[ais]][[ss]] <- NA
      Hpre[[rcp]][[ais]][[ss]] <- NA
      ireg[[rcp]][[ais]][[ss]] <- NA
      Creg[[rcp]][[ais]][[ss]] <- NA
      Rreg[[rcp]][[ais]][[ss]] <- NA
      Hreg[[rcp]][[ais]][[ss]] <- NA
      icva[[rcp]][[ais]][[ss]] <- NA
      Ccva[[rcp]][[ais]][[ss]] <- NA
      Rcva[[rcp]][[ais]][[ss]] <- NA
      Hcva[[rcp]][[ais]][[ss]] <- NA
      iopt[[rcp]][[ais]][[ss]] <- rep(NA,n.ensemble[[ais]])
      Copt[[rcp]][[ais]][[ss]] <- rep(NA,n.ensemble[[ais]])
      Ropt[[rcp]][[ais]][[ss]] <- rep(NA,n.ensemble[[ais]])
      Hopt[[rcp]][[ais]][[ss]] <- rep(NA,n.ensemble[[ais]])
      V.h[[rcp]][[ais]][[ss]] <- rep(NA,n.height)
      C.h[[rcp]][[ais]][[ss]] <- rep(NA,n.height)
      R.h[[rcp]][[ais]][[ss]] <- rep(NA,n.height)
    }
  }
}

# calculate the optimal strategy (EUM, min-Regret, ...?) for each SOW/ensemble
# member, in each scenario

for (ais in scen.ais) {
  for (rcp in scen.rcp) {
    for (ss in scen.ss) {

      # calculate the economically efficient solution for each SOW, within each scenario
      # (needed for the regret calculation)
      for (sow in 1:n.ensemble[[ais]]) {
        Copt[[rcp]][[ais]][[ss]][sow] <- min(cost[[rcp]][[ais]][[ss]][,sow])
        itmp <- which(cost[[rcp]][[ais]][[ss]][,sow]==Copt[[rcp]][[ais]][[ss]][sow])
        if(length(itmp)>1) {itmp <- median(itmp)}
        iopt[[rcp]][[ais]][[ss]][sow] <- itmp
        Hopt[[rcp]][[ais]][[ss]][sow] <- heightening[itmp]
        Ropt[[rcp]][[ais]][[ss]][sow] <- preturn[[rcp]][[ais]][[ss]][itmp,sow]
      }
      for (h in 1:n.height) {
        C.h[[rcp]][[ais]][[ss]][h] <- mean( cost[[rcp]][[ais]][[ss]][h,] )
        R.h[[rcp]][[ais]][[ss]][h] <- mean( preturn[[rcp]][[ais]][[ss]][h,] )
        V.h[[rcp]][[ais]][[ss]][h] <- mean( cost[[rcp]][[ais]][[ss]][h,] - Copt[[rcp]][[ais]][[ss]] )

        # tally up for Laplace averaging
        Cavg[h] <- (navg/(navg+n.ensemble[[ais]]))*Cavg[h] + (n.ensemble[[ais]]/(navg+n.ensemble[[ais]]))*C.h[[rcp]][[ais]][[ss]][h]
        Ravg[h] <- (navg/(navg+n.ensemble[[ais]]))*Ravg[h] + (n.ensemble[[ais]]/(navg+n.ensemble[[ais]]))*R.h[[rcp]][[ais]][[ss]][h]
        Vavg[h] <- (navg/(navg+n.ensemble[[ais]]))*Vavg[h] + (n.ensemble[[ais]]/(navg+n.ensemble[[ais]]))*V.h[[rcp]][[ais]][[ss]][h]
      }

      # Expected utility maximization
      Ceum[[rcp]][[ais]][[ss]] <- min(C.h[[rcp]][[ais]][[ss]])
      itmp <- which(C.h[[rcp]][[ais]][[ss]]==min(C.h[[rcp]][[ais]][[ss]]))
      if(length(itmp)>1) {itmp <- median(itmp)}
      ieum[[rcp]][[ais]][[ss]] <- itmp
      Heum[[rcp]][[ais]][[ss]] <- heightening[ieum[[rcp]][[ais]][[ss]]]
      Reum[[rcp]][[ais]][[ss]] <- R.h[[rcp]][[ais]][[ss]][itmp]
      Cens.eum[[rcp]][[ais]][[ss]] <- cost[[rcp]][[ais]][[ss]][itmp,]
      Rens.eum[[rcp]][[ais]][[ss]] <- preturn[[rcp]][[ais]][[ss]][itmp,]

      # Expected regret minimization
      Creg[[rcp]][[ais]][[ss]] <- min(V.h[[rcp]][[ais]][[ss]])
      itmp <- which(V.h[[rcp]][[ais]][[ss]]==min(V.h[[rcp]][[ais]][[ss]]))
      if(length(itmp)>1) {itmp <- median(itmp)}
      ireg[[rcp]][[ais]][[ss]] <- itmp
      Hreg[[rcp]][[ais]][[ss]] <- heightening[ieum[[rcp]][[ais]][[ss]]]
      Rreg[[rcp]][[ais]][[ss]] <- R.h[[rcp]][[ais]][[ss]][itmp]

      # CVaR minimization; use Ctmp as temporary storage.
      for (h in 1:n.height) {
        Ctmp[h] <- mean( rev(sort(cost[[rcp]][[ais]][[ss]][h,]))[1:ceiling(alpha*n.ensemble[[ais]])] )
      }
      Ccva[[rcp]][[ais]][[ss]] <- min(Ctmp)
      itmp <- which(Ctmp==min(Ctmp))
      if(length(itmp)>1) {itmp <- median(itmp)}
      icva[[rcp]][[ais]][[ss]] <- itmp
      Hcva[[rcp]][[ais]][[ss]] <- heightening[icva[[rcp]][[ais]][[ss]]]
      Rcva[[rcp]][[ais]][[ss]] <- R.h[[rcp]][[ais]][[ss]][itmp]

    }
  }
}

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================

#install.packages("extRemes")
#install.packages("fExtremes")
#install.packages('ismev')
#install.packages('zoo')
library(extRemes)
library(fExtremes)
library(ismev)
library(zoo)

# nodes for pdfs of sea-level rise, storm surge, and flood risk, and list
# vectors for storing the density estimates for each scenario
f.sealev <- init
f.surlev <- init
cdf.sealev <- init
cdf.surlev <- init
sf.sealev <- init
sf.surlev <- init
f.surge.rise <- init
f.seasublev <- init
f.seasurlev <- init
cdf.seasurlev <- init
sf.seasurlev <- init
lsl.norm <- init

lsl.lower <- 0
lsl.upper <- 20
lsl.n <- 2^11 # powers of 2 are good for density estimation in R
              # 2^10 = 1024 ~ number of ensemble members

flood.n <- lsl.n

ncdata <- nc_open(filename.brick.allslr)
  mod.time <- ncvar_get(ncdata, 'time_proj')
  lsl$rcp26$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP26')
  lsl$rcp45$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP45')
  lsl$rcp85$none <- ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP85')
  lsl$rcp26$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP26')
  lsl$rcp45$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP45')
  lsl$rcp85$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_gamma_RCP85')
  lsl$rcp26$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP26')
  lsl$rcp45$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP45')
  lsl$rcp85$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_uniform_RCP85')
nc_close(ncdata)

## initialize
f.surge.stat <- vector("list",length(scen.rcp)); names(f.surge.stat) <- scen.rcp
cdf.surge.stat <- vector("list",length(scen.rcp)); names(cdf.surge.stat) <- scen.rcp
sf.surge.stat <- vector("list",length(scen.rcp)); names(sf.surge.stat) <- scen.rcp
f.surge.avg <- vector("list",length(scen.rcp)); names(f.surge.avg) <- scen.rcp
for (rcp in scen.rcp) {
  f.surge.stat[[rcp]] <- vector("list",length(scen.ais)); names(f.surge.stat[[rcp]]) <- scen.ais
  cdf.surge.stat[[rcp]] <- vector("list",length(scen.ais)); names(cdf.surge.stat[[rcp]]) <- scen.ais
  sf.surge.stat[[rcp]] <- vector("list",length(scen.ais)); names(sf.surge.stat[[rcp]]) <- scen.ais
  f.surge.avg[[rcp]] <- vector("list",length(scen.ais)); names(f.surge.avg[[rcp]]) <- scen.ais
  for (ais in scen.ais) {
    f.surlev[[rcp]][[ais]] <- vector("list",2); names(f.surlev[[rcp]][[ais]])=scen.ss
    f.seasurlev[[rcp]][[ais]] <- vector("list",2); names(f.seasurlev[[rcp]][[ais]])=scen.ss
    f.surge.stat[[rcp]][[ais]] <- mat.or.vec(n.ensemble[[ais]],flood.n)
    cdf.surge.stat[[rcp]][[ais]] <- mat.or.vec(n.ensemble[[ais]],flood.n)
    sf.surge.stat[[rcp]][[ais]] <- mat.or.vec(n.ensemble[[ais]],flood.n)
  }
}

tmp <- density(x=lsl$rcp26$none, from=lsl.lower, to=lsl.upper, n=lsl.n)
x.lsl <- tmp$x
dx.lsl <- mean(diff(x.lsl))
kern <- 'gaussian'



##TODO
# change all integrals to simpson's rule
library(Bolstad)
#sintegral(x=x.lsl, fx=pdf.fit[1,])$value
##TODO




## Get the return levels for the stationary GEV storm surge case
surge.level <- n.ensemble['none'] # desired storm surge level;
  # determine this resolution by size of ensemble (can't resolve the 1% tail if
  # your ensemble is only 99 members, for example)
q <- seq(0,1,length.out= surge.level +1)  # quantile array

# Find fit to the stationary surge distirbutions
cdf.fit <- mat.or.vec(n.ensemble['none'],lsl.n)
pdf.fit <- mat.or.vec(n.ensemble['none'],lsl.n)
sf.fit <- mat.or.vec(n.ensemble['none'],lsl.n)
fit.surge <- mat.or.vec(n.ensemble['none'],length(q))

# These cdf and survival function fits are not
for (i in 1:n.ensemble['none']) {
    cdf.fit[i,] <- pevd(1000*x.lsl, shape=gev.stat[i,'shape'], loc=gev.stat[i,'location'], scale=gev.stat[i,'scale'])
    pdf.fit[i,] <- c(0, diff(cdf.fit[i,]))
    sf.fit[i,] <- 1-cdf.fit[i,]
    fit.surge[i,] <- qgev(q, xi=gev.stat[i,'shape'], mu=gev.stat[i,'location'], beta=gev.stat[i,'scale'])
}
fit.surge <- fit.surge[,-n.ensemble['none']]

# Plot them all
#plot(fit.surge[1,]/1000, log10(1-q[2:1294]), type='l', xlim=c(0,10), ylim=c(-4,0)); for (i in 1:n.ensemble['none']) {lines(fit.surge[i,]/1000, log10(1-q[2:1294]))}

if(FALSE) {

    # use the old surge-level fitting way (kind of a riff of Kelsey's SF codes)

# Get the pdfs for each SOW stationary storm surge
f.tmp <- mat.or.vec(n.ensemble['none'], n.ensemble['none']+1)
for (sow in 1:n.ensemble['none']) {
  for (node in 2:(n.ensemble['none']+1)) {
    f.tmp[sow,node-1] <- (q[node]-q[node-1])/(fit.surge[sow,node]-fit.surge[sow,node-1])
  }
}
# the above bit leaves:
# f.tmp = [N.ens x N.node], N.node=N.ens+1 (base maximum surge level we can get off ensemble size)
#                           f.tmp is the probability density of the stationary surge for each SOW
# fit.surge = [N.ens x N.node], gives the locations (along sea level axis) of the
#                               nodes for the f.tmp fit
# Now interpolate to the x.lsl grid
f.tmp.lsl <- mat.or.vec(n.ensemble['none'], lsl.n)
for (sow in 1:n.ensemble['none']) {
  f.tmp.lsl[sow,] <- approx(x=fit.surge[sow,]/1000, y=f.tmp[sow,], xout=x.lsl)$y
  f.tmp.lsl[sow,] <- f.tmp.lsl[sow,]/sum(dx.lsl*f.tmp.lsl[sow,], na.rm=TRUE)
}
f.tmp.lsl[which(is.nan(f.tmp.lsl))] <- 0

# Cut off the storm surge stationary distributions beyond what we can resolve
# with our model
max.surge <- fit.surge[,surge.level]/1000

} else {

    # use the pevd (fit CDF) based distributions
    f.tmp.lsl <- mat.or.vec(n.ensemble['none'], lsl.n)
    for (sow in 1:n.ensemble['none']) {
        f.tmp.lsl[sow,] <- pdf.fit[sow,]/sintegral(x=x.lsl, fx=pdf.fit[sow,])$value
    }

}


## distribution of subsidence
iproj <- which(mod.time==2065)
inorm <- which(mod.time==2015)

subs <- VD.params[,6] # m/year
f.subs <- density(x=subs*(mod.time[iproj]-mod.time[inorm]), from=lsl.lower, to=lsl.upper, n=lsl.n, kernel=kern)
f.subs <- f.subs$y/sum(f.subs$y*dx.lsl) # normalize

## distribution of sea-level, subsidence, and storm surges
for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    lsl.norm[[rcp]][[ais]] <- lsl[[rcp]][[ais]][iproj,] - lsl[[rcp]][[ais]][inorm,]
    f.sealev[[rcp]][[ais]] <- density(x=lsl.norm[[rcp]][[ais]], from=lsl.lower, to=lsl.upper, n=lsl.n, kernel=kern)
    f.sealev[[rcp]][[ais]] <- f.sealev[[rcp]][[ais]]$y/sum(f.sealev[[rcp]][[ais]]$y*dx.lsl) # normalize
    cdf.sealev[[rcp]][[ais]] <- cumsum(f.sealev[[rcp]][[ais]]*dx.lsl)
    cdf.sealev[[rcp]][[ais]] <- cdf.sealev[[rcp]][[ais]]-cdf.sealev[[rcp]][[ais]][1] # normalize (beginning=0)
    cdf.sealev[[rcp]][[ais]] <- cdf.sealev[[rcp]][[ais]]/cdf.sealev[[rcp]][[ais]][lsl.n] # normalize (end=1)
    sf.sealev[[rcp]][[ais]] <- 1-cdf.sealev[[rcp]][[ais]]
    f.surge.rise[[rcp]][[ais]] <- density(x=lsl.norm[[rcp]][[ais]]*surge.factor, from=lsl.lower, to=lsl.upper, n=lsl.n, kernel=kern)
    f.surge.rise[[rcp]][[ais]] <- f.surge.rise[[rcp]][[ais]]$y/sum(f.surge.rise[[rcp]][[ais]]$y*dx.lsl) # normalize
    # and initialize storm surge
    for (ss in scen.ss) {
      f.surlev[[rcp]][[ais]][[ss]] <- rep(NA, lsl.n)
    }
    # and sea level, accountign for subsidence as effectively more sea-level rise
    f.seasublev[[rcp]][[ais]] <- convolve(x=f.sealev[[rcp]][[ais]], y=rev(f.subs))
    f.seasublev[[rcp]][[ais]] <- f.seasublev[[rcp]][[ais]]/sum(f.seasublev[[rcp]][[ais]]*dx.lsl) # normalize
  }
}

## distribution of storm surge level
## convolution of surge.rise distribution + gev.stat distribution

for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    for (sow in 1:n.ensemble[[ais]]) {
      f.surge.stat[[rcp]][[ais]][sow,] <- as.numeric(dgev(x=x.lsl*1000,
                                                        xi=gev.stat[sow,'shape'],
                                                        mu=gev.stat[sow,'location'],
                                                        beta=gev.stat[sow,'scale'],
                                                        log=FALSE))
      f.surge.stat[[rcp]][[ais]][sow,which(x.lsl > max.surge[sow])] <- 0 # cut off above maximum surge
      f.surge.stat[[rcp]][[ais]][sow,] <- f.surge.stat[[rcp]][[ais]][sow,]/sum( dx.lsl*f.surge.stat[[rcp]][[ais]][sow,])
      cdf.surge.stat[[rcp]][[ais]][sow,] <- cumsum(f.surge.stat[[rcp]][[ais]][sow,]*dx.lsl)
      sf.surge.stat[[rcp]][[ais]][sow,] <- 1-cdf.surge.stat[[rcp]][[ais]][sow,]
    }

    # one option: average the probabilities across all ensemble members, and re-normalize
    if(FALSE){
    f.surge.avg[[rcp]][[ais]] <- apply(f.surge.stat[[rcp]][[ais]], 2, mean)
    f.surge.avg[[rcp]][[ais]] <- f.surge.avg[[rcp]][[ais]]/sum(dx.lsl*f.surge.avg[[rcp]][[ais]])
    for (ss in scen.ss) {
      if(ss=='st') {
        f.surlev[[rcp]][[ais]][[ss]] <- f.surge.avg[[rcp]][[ais]]
        f.surlev[[rcp]][[ais]][[ss]] <- f.surlev[[rcp]][[ais]][[ss]]/sum(f.surlev[[rcp]][[ais]][[ss]]*dx.lsl) # normalize
        f.seasurlev[[rcp]][[ais]][[ss]] <- convolve(x=f.surlev[[rcp]][[ais]][[ss]], y=rev(f.sealev[[rcp]][[ais]]))
        f.seasurlev[[rcp]][[ais]][[ss]] <- f.seasurlev[[rcp]][[ais]][[ss]]/sum(dx.lsl*f.seasurlev[[rcp]][[ais]][[ss]]) # normalize
      } else {
        f.surlev[[rcp]][[ais]][[ss]] <- convolve(x=f.surge.rise[[rcp]][[ais]], y=rev(f.surge.avg[[rcp]][[ais]]))
        f.surlev[[rcp]][[ais]][[ss]] <- f.surlev[[rcp]][[ais]][[ss]]/sum(dx.lsl*f.surlev[[rcp]][[ais]][[ss]]) # normalize
        f.seasurlev[[rcp]][[ais]][[ss]] <- convolve(x=f.surlev[[rcp]][[ais]][[ss]], y=rev(f.sealev[[rcp]][[ais]]))
        f.seasurlev[[rcp]][[ais]][[ss]] <- f.seasurlev[[rcp]][[ais]][[ss]]/sum(dx.lsl*f.seasurlev[[rcp]][[ais]][[ss]]) # normalize
      }
      cdf.seasurlev[[rcp]][[ais]][[ss]] <- cumsum(f.seasurlev[[rcp]][[ais]][[ss]]*dx.lsl)
      cdf.seasurlev[[rcp]][[ais]][[ss]] <- cdf.seasurlev[[rcp]][[ais]][[ss]]-cdf.seasurlev[[rcp]][[ais]][[ss]][1] # normalize (beginning=0)
      cdf.seasurlev[[rcp]][[ais]][[ss]] <- cdf.seasurlev[[rcp]][[ais]][[ss]]/cdf.seasurlev[[rcp]][[ais]][[ss]][lsl.n] # normalize (end=1)
      sf.seasurlev[[rcp]][[ais]][[ss]] <- 1-cdf.seasurlev[[rcp]][[ais]][[ss]]
      cdf.surlev[[rcp]][[ais]][[ss]] <- cumsum(f.surlev[[rcp]][[ais]][[ss]]*dx.lsl)
      cdf.surlev[[rcp]][[ais]][[ss]] <- cdf.surlev[[rcp]][[ais]][[ss]]-cdf.surlev[[rcp]][[ais]][[ss]][1] # normalize (beginning=0)
      cdf.surlev[[rcp]][[ais]][[ss]] <- cdf.surlev[[rcp]][[ais]][[ss]]/cdf.surlev[[rcp]][[ais]][[ss]][lsl.n] # normalize (end=1)
      sf.surlev[[rcp]][[ais]][[ss]] <- 1-cdf.surlev[[rcp]][[ais]][[ss]]
    }
    }
    # second option: get seasurlev for each SOW - this one is good, gets distribution of return periods as it should be
    f.surlev[[rcp]][[ais]] <- vector('list', length(scen.ss)); names(f.surlev[[rcp]][[ais]]) <- scen.ss
    f.seasurlev[[rcp]][[ais]] <- vector('list', length(scen.ss)); names(f.seasurlev[[rcp]][[ais]]) <- scen.ss
    for (ss in scen.ss) {
      f.surlev[[rcp]][[ais]][[ss]] <- mat.or.vec(n.ensemble[[ais]], lsl.n)
      cdf.surlev[[rcp]][[ais]][[ss]] <- mat.or.vec(n.ensemble[[ais]], lsl.n)
      sf.surlev[[rcp]][[ais]][[ss]] <- mat.or.vec(n.ensemble[[ais]], lsl.n)
      f.seasurlev[[rcp]][[ais]][[ss]] <- mat.or.vec(n.ensemble[[ais]], lsl.n)
      cdf.seasurlev[[rcp]][[ais]][[ss]] <- mat.or.vec(n.ensemble[[ais]], lsl.n)
      sf.seasurlev[[rcp]][[ais]][[ss]] <- mat.or.vec(n.ensemble[[ais]], lsl.n)
      for (sow in 1:n.ensemble[[ais]]) {
        if(ss=='st') {
          f.surlev[[rcp]][[ais]][[ss]][sow,] <- f.surge.stat[[rcp]][[ais]][sow,]
          #f.surlev[[rcp]][[ais]][[ss]][sow,] <- f.tmp.lsl[sow,]
          f.surlev[[rcp]][[ais]][[ss]][sow,] <- f.surlev[[rcp]][[ais]][[ss]][sow,]/sum(f.surlev[[rcp]][[ais]][[ss]][sow,]*dx.lsl) # normalize
          f.seasurlev[[rcp]][[ais]][[ss]][sow,] <- convolve(x=f.surlev[[rcp]][[ais]][[ss]][sow,], y=rev(f.seasublev[[rcp]][[ais]]))
          f.seasurlev[[rcp]][[ais]][[ss]][sow,] <- f.seasurlev[[rcp]][[ais]][[ss]][sow,]/sum(dx.lsl*f.seasurlev[[rcp]][[ais]][[ss]][sow,]) # normalize
        } else {
          f.surlev[[rcp]][[ais]][[ss]][sow,] <- convolve(x=f.surge.rise[[rcp]][[ais]], y=rev(f.surge.stat[[rcp]][[ais]][sow,]))
          #f.surlev[[rcp]][[ais]][[ss]][sow,] <- convolve(x=f.surge.rise[[rcp]][[ais]], y=rev(f.tmp.lsl[sow,]))
          f.surlev[[rcp]][[ais]][[ss]][sow,] <- f.surlev[[rcp]][[ais]][[ss]][sow,]/sum(f.surlev[[rcp]][[ais]][[ss]][sow,]*dx.lsl) # normalize
          f.seasurlev[[rcp]][[ais]][[ss]][sow,] <- convolve(x=f.surlev[[rcp]][[ais]][[ss]][sow,], y=rev(f.seasublev[[rcp]][[ais]]))
          f.seasurlev[[rcp]][[ais]][[ss]][sow,] <- f.seasurlev[[rcp]][[ais]][[ss]][sow,]/sum(dx.lsl*f.seasurlev[[rcp]][[ais]][[ss]][sow,]) # normalize
        }
        cdf.seasurlev[[rcp]][[ais]][[ss]][sow,] <- cumsum(f.seasurlev[[rcp]][[ais]][[ss]][sow,]*dx.lsl)
        cdf.seasurlev[[rcp]][[ais]][[ss]][sow,] <- cdf.seasurlev[[rcp]][[ais]][[ss]][sow,]-cdf.seasurlev[[rcp]][[ais]][[ss]][sow,1] # normalize (beginning=0)
        cdf.seasurlev[[rcp]][[ais]][[ss]][sow,] <- cdf.seasurlev[[rcp]][[ais]][[ss]][sow,]/cdf.seasurlev[[rcp]][[ais]][[ss]][sow,lsl.n] # normalize (end=1)
        sf.seasurlev[[rcp]][[ais]][[ss]][sow,] <- 1-cdf.seasurlev[[rcp]][[ais]][[ss]][sow,]
        cdf.surlev[[rcp]][[ais]][[ss]][sow,] <- cumsum(f.surlev[[rcp]][[ais]][[ss]][sow,]*dx.lsl)
        cdf.surlev[[rcp]][[ais]][[ss]][sow,] <- cdf.surlev[[rcp]][[ais]][[ss]][sow,]-cdf.surlev[[rcp]][[ais]][[ss]][sow,1] # normalize (beginning=0)
        cdf.surlev[[rcp]][[ais]][[ss]][sow,] <- cdf.surlev[[rcp]][[ais]][[ss]][sow,]/cdf.surlev[[rcp]][[ais]][[ss]][sow,lsl.n] # normalize (end=1)
        sf.surlev[[rcp]][[ais]][[ss]][sow,] <- 1-cdf.surlev[[rcp]][[ais]][[ss]][sow,]
      }
    }
  }
}

## Return periods.
## Calculate as the tail probability above the height of the levee system (h0=16ft)
subs.rate <- 0.0056
H0 <- 16*0.3048 # initial levee height
h0 <- H0 #- subs.rate*(mod.time[iproj]-mod.time[inorm])  # ... reduced by subsidence (already accounted for in f.subs)
h0 <- h0 - (4*0.3048) # reduced by elevation already subsided below sea level (p. 134 of USACE manual)
iflood <- which(x.lsl > h0)

# if overtopping is XX fraction of total failure probability, then
# p_fail_total = (1/XX)*p_fail_overtopping
# return_period = 1/p_fail_total = XX*return_period_overtopping
return.period <- mat.or.vec(18*n.ensemble['none'], 2); cnt <- 1
return.period.build <- mat.or.vec(18*n.ensemble['none'], 2); cnt <- 1
return.period.fragile60 <- mat.or.vec(18*n.ensemble['none'], 2); cnt <- 1
return.period.fragile80 <- mat.or.vec(18*n.ensemble['none'], 2); cnt <- 1
scen.names <- rep(NA,18)
rcp.tmp <- c('RCP2.6','RCP4.5','RCP8.5'); names(rcp.tmp) <- scen.rcp
for (rcp in scen.rcp) {for (ais in scen.ais) {for (ss in scen.ss) {
    scen.names[cnt] <- paste(rcp.tmp[rcp],ais,ss)
    return.period[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 1] <- cnt
    return.period.build[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 1] <- cnt
    return.period.fragile60[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 1] <- cnt
    return.period.fragile80[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 1] <- cnt
    return.period[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- preturn[[rcp]][[ais]][[ss]][1,] # first height is 0 m build
    return.period.build[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- preturn[[rcp]][[ais]][[ss]][2,] # second height is 0.91 m build
    return.period.fragile60[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- (6/10)*preturn[[rcp]][[ais]][[ss]][1,]
    return.period.fragile80[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- (8/10)*preturn[[rcp]][[ais]][[ss]][1,]
    #return.period[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- rp[[rcp]][[ais]][[ss]]
    #return.period[(n.ensemble['none']*(cnt-1)+1):(n.ensemble['none']*cnt), 2] <- 1/apply(dx.lsl*f.seasurlev[[rcp]][[ais]][[ss]][,iflood], 1, sum)
    cnt <- cnt+1
}}}
colnames(return.period) <- c('Scenario','ReturnPeriod')
colnames(return.period.build) <- c('Scenario','ReturnPeriod')
colnames(return.period.fragile60) <- c('Scenario','ReturnPeriod')
colnames(return.period.fragile80) <- c('Scenario','ReturnPeriod')

##==============================================================================
##==============================================================================
## Sort out Sobol sensitivity analysis for drivers of flood risk

source('BRICK_nola_scenarios_sobol.R')

##==============================================================================
##==============================================================================





##==============================================================================
## FIGURES
##==============================================================================






##==============================================================================

##=========     (column 1) pdfs of sea level rise and storm surge, all scenarios
## FIGURE 1     (column 2) corresponding survival functions
##=========

pdf(paste(plotdir,'pdfs_sf_slr_surge.pdf',sep=''),width=8,height=7,colormodel='cmyk')

par(mfrow=c(2,2), mai=c(.6,.63,.2,.26))
# (a) pdfs of local sea-level rise (2065)
plot(x.lsl, f.sealev$rcp26$none, type='l', xlim=c(0,1.1), ylim=c(0,10), lwd=1.5,
     col=rgb(col26[1],col26[2],col26[3]), lty=1,
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(x.lsl, f.sealev$rcp26$gamma, type='l', col=rgb(col26[1],col26[2],col26[3]), lty=2)
lines(x.lsl, f.sealev$rcp26$uniform, type='l', col=rgb(col26[1],col26[2],col26[3]), lty=3)
lines(x.lsl, f.sealev$rcp45$none, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=1)
lines(x.lsl, f.sealev$rcp45$gamma, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=2)
lines(x.lsl, f.sealev$rcp45$uniform, type='l', col=rgb(col45[1],col45[2],col45[3]), lty=3)
lines(x.lsl, f.sealev$rcp85$none, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=1)
lines(x.lsl, f.sealev$rcp85$gamma, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=2)
lines(x.lsl, f.sealev$rcp85$uniform, type='l', col=rgb(col85[1],col85[2],col85[3]), lty=3)

axis(1,seq(0,1,0.2),lab=c("0","0.2","0.4","0.6","0.8","1"), cex.axis=1.2)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1.2, cex=1);
mtext('Projected sea level in 2065 [m]', side=1, line=2.3, cex=1);
mtext(side=3, text=expression(bold('   a')), line=-1, cex=.9, adj=0);

legend(0.5,9.5,c("RCP2.6","RCP4.5","RCP8.5","no FD","FD, gamma","FD, uniform"),
       lty=c(1,1,1,1,2,3), lwd=2, cex=1.2,
       col=c(rgb(col26[1],col26[2],col26[3]),rgb(col45[1],col45[2],col45[3]),rgb(col85[1],col85[2],col85[3]),'black','black','black'),
       bty='n')

# (b) pdfs of storm surge
par(mai=c(.6,.63,.2,.26))
plot(x.lsl, apply(f.surlev$rcp26$none$st,2,mean), type='l', xlim=c(0,4), ylim=c(0,3), lwd=1.5,
     col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lty=1,
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(x.lsl, apply(f.surlev$rcp26$none$ns,2,mean), type='l', col=rgb(col26[1],col26[2],col26[3]), lty=1)
lines(x.lsl, apply(f.surlev$rcp26$gamma$ns,2,mean), type='l', col=rgb(col26[1],col26[2],col26[3]), lty=2)
lines(x.lsl, apply(f.surlev$rcp26$uniform$ns,2,mean), type='l', col=rgb(col26[1],col26[2],col26[3]), lty=3)
lines(x.lsl, apply(f.surlev$rcp45$none$ns,2,mean), type='l', col=rgb(col45[1],col45[2],col45[3]), lty=1)
lines(x.lsl, apply(f.surlev$rcp45$gamma$ns,2,mean), type='l', col=rgb(col45[1],col45[2],col45[3]), lty=2)
lines(x.lsl, apply(f.surlev$rcp45$uniform$ns,2,mean), type='l', col=rgb(col45[1],col45[2],col45[3]), lty=3)
lines(x.lsl, apply(f.surlev$rcp85$none$ns,2,mean), type='l', col=rgb(col85[1],col85[2],col85[3]), lty=1)
lines(x.lsl, apply(f.surlev$rcp85$gamma$ns,2,mean), type='l', col=rgb(col85[1],col85[2],col85[3]), lty=2)
lines(x.lsl, apply(f.surlev$rcp85$uniform$ns,2,mean), type='l', col=rgb(col85[1],col85[2],col85[3]), lty=3)

axis(1,seq(0,4,0.5),lab=c("0","","1","","2","","3","","4"), cex.axis=1.2)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1.2, cex=1);
mtext('Projected surge level in 2065 [m]', side=1, line=2.3, cex=1);
mtext(side=3, text=expression(bold('   b')), line=-1, cex=.9, adj=0);
text(.5,2.6,"stationary", pos=4, cex=1.3)
text(1.1,2,"non-stationary", pos=4, cex=1.3)

# (c) survival functions of sea-level rise
par(mai=c(.6,.63,.2,.26))
# rcp2.6
plot( x.lsl, log10(sf.sealev$rcp26$none), type='l', xlim=c(0,1), ylim=c(-3.2,0),
      lty=1, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i')
lines(x.lsl, log10(sf.sealev$rcp26$gamma), type='l',
      lty=2, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(x.lsl, log10(sf.sealev$rcp26$uniform), type='l',
      lty=3, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
# rcp4.5
lines(x.lsl, log10(sf.sealev$rcp45$none), type='l',
      lty=1, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(x.lsl, log10(sf.sealev$rcp45$gamma), type='l',
      lty=2, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(x.lsl, log10(sf.sealev$rcp45$uniform), type='l',
      lty=3, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
# rcp8.5
lines(x.lsl, log10(sf.sealev$rcp85$none), type='l',
      lty=1, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(x.lsl, log10(sf.sealev$rcp85$gamma), type='l',
      lty=2, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(x.lsl, log10(sf.sealev$rcp85$uniform), type='l',
      lty=3, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
# x-axis, ticks and text
axis(1,seq(0,1,0.2),lab=c("0","0.2","0.4","0.6","0.8","1"), cex.axis=1.2)
mtext('Projected sea level in 2065 [m]', side=1, line=2.3, cex=1);
# y-axis, ticks and text
axis(2, at=seq(-3,0), label=parse(text=paste("10^", seq(-3,0), sep="")), las=1, cex.axis=1.2, mgp=c(3,.7,0))
mtext('Survival function [1-CDF]', side = 2, line=2.6);
# add a panel label
mtext(side=3, text=expression(bold('   c')), line=-1.1, cex=.9, adj=0);
# add the 1/100, 1/500, 1/1000 levels
lines(c(-4,4),c(-2,-2),lty=2,col='black'); text(0.16,-1.85,"1/100 level", cex=1.2)
lines(c(-4,4),c(-log10(500),-log10(500)),lty=2,col='black'); text(0.16,-2.55,"1/500 level", cex=1.2);
#lines(c(-4,4),c(-3,-3),lty=2,col='black'); text(0.16,-2.85,"1/1000 level", cex=1.2)

# (d) survival functions of storm surge
par(mai=c(.6,.63,.2,.26))
# rcp2.6
plot( x.lsl, log10(apply(sf.surlev$rcp26$none$st,2,mean)), type='l', xlim=c(0,4), ylim=c(-2.2,0),
      lty=1, col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=1.5, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i')
lines(x.lsl, log10(apply(sf.surlev$rcp26$none$ns,2,mean)), type='l',
      lty=1, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.surlev$rcp26$gamma$ns,2,mean)), type='l',
      lty=2, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.surlev$rcp26$uniform$ns,2,mean)), type='l',
      lty=3, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
# rcp4.5
lines(x.lsl, log10(apply(sf.surlev$rcp45$none$ns,2,mean)), type='l',
      lty=1, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.surlev$rcp45$gamma$ns,2,mean)), type='l',
      lty=2, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.surlev$rcp45$uniform$ns,2,mean)), type='l',
      lty=3, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
# rcp8.5
lines(x.lsl, log10(apply(sf.surlev$rcp85$none$ns,2,mean)), type='l',
      lty=1, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.surlev$rcp85$gamma$ns,2,mean)), type='l',
      lty=2, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.surlev$rcp85$uniform$ns,2,mean)), type='l',
      lty=3, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
# x-axis, ticks and text
axis(1,seq(0,4,0.5),lab=c("0","","1","","2","","3","","4"), cex.axis=1.2)
mtext('Projected surge level in 2065 [m]', side=1, line=2.3, cex=1);
# y-axis, ticks and text
axis(2, at=seq(-3,0), label=parse(text=paste("10^", seq(-3,0), sep="")), las=1, cex.axis=1.2, mgp=c(3,.7,0))
mtext('Survival function [1-CDF]', side = 2, line=2.6);
# add a panel label
mtext(side=3, text=expression(bold('   d')), line=-1.1, cex=.9, adj=0);
# add the 1/100, 1/500, 1/1000 levels
lines(c(-4,4),c(-2,-2),lty=2,col='black'); text(0.6,-1.88,"1/100 level", cex=1.2)
#lines(c(-4,4),c(-log10(500),-log10(500)),lty=2,col='black'); text(0.6,-1.55,"1/500 level", cex=1.2);
#lines(c(-4,4),c(-3,-3),lty=2,col='black'); text(0.1,-2.85,"1/1,000 level", cex=1.2)

dev.off()

##==============================================================================








##==============================================================================

##=========   (row 1) pdfs of total sea+storm level, all scenarios
## FIGURE 2   (row 2) survival functions
##=========   (row 3) (if needed) corresponding survival functions?

pdf(paste(plotdir,'pdfs_sf_sea+surge.pdf',sep=''),width=3.5,height=7,colormodel='cmyk')

par(mfrow=c(3,1), mai=c(.25,.63,.1,.2))
# (a) pdfs of sea-level rise + storm surge
# stationary
plot(x.lsl, apply(f.seasurlev$rcp26$none$st,2,mean), type='l', xlim=c(0,5), ylim=c(0,2), lwd=1.5,
     col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lty=1,
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(x.lsl, apply(f.seasurlev$rcp26$gamma$st,2,mean), type='l', col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lty=2)
lines(x.lsl, apply(f.seasurlev$rcp26$uniform$st,2,mean), type='l', col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lty=3)
lines(x.lsl, apply(f.seasurlev$rcp45$none$st,2,mean), type='l', col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lty=1)
lines(x.lsl, apply(f.seasurlev$rcp45$gamma$st,2,mean), type='l', col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lty=2)
lines(x.lsl, apply(f.seasurlev$rcp45$uniform$st,2,mean), type='l', col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lty=3)
lines(x.lsl, apply(f.seasurlev$rcp85$none$st,2,mean), type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lty=1)
lines(x.lsl, apply(f.seasurlev$rcp85$gamma$st,2,mean), type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lty=2)
lines(x.lsl, apply(f.seasurlev$rcp85$uniform$st,2,mean), type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lty=3)
# non-stationary
lines(x.lsl, apply(f.seasurlev$rcp26$none$ns,2,mean), type='l', col=rgb(col26[1],col26[2],col26[3]), lty=1)
lines(x.lsl, apply(f.seasurlev$rcp26$gamma$ns,2,mean), type='l', col=rgb(col26[1],col26[2],col26[3]), lty=2)
lines(x.lsl, apply(f.seasurlev$rcp26$uniform$ns,2,mean), type='l', col=rgb(col26[1],col26[2],col26[3]), lty=3)
lines(x.lsl, apply(f.seasurlev$rcp45$none$ns,2,mean), type='l', col=rgb(col45[1],col45[2],col45[3]), lty=1)
lines(x.lsl, apply(f.seasurlev$rcp45$gamma$ns,2,mean), type='l', col=rgb(col45[1],col45[2],col45[3]), lty=2)
lines(x.lsl, apply(f.seasurlev$rcp45$uniform$ns,2,mean), type='l', col=rgb(col45[1],col45[2],col45[3]), lty=3)
lines(x.lsl, apply(f.seasurlev$rcp85$none$ns,2,mean), type='l', col=rgb(col85[1],col85[2],col85[3]), lty=1)
lines(x.lsl, apply(f.seasurlev$rcp85$gamma$ns,2,mean), type='l', col=rgb(col85[1],col85[2],col85[3]), lty=2)
lines(x.lsl, apply(f.seasurlev$rcp85$uniform$ns,2,mean), type='l', col=rgb(col85[1],col85[2],col85[3]), lty=3)

axis(1,seq(0,5,0.5),lab=c("0","","1","","2","","3","","4","","5"), cex.axis=1.3)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, length=.15, xpd = TRUE)
mtext('Probability density', side=2, line=1.2, cex=.8);
mtext('Projected sea+surge level in 2065 [m]', side=1, line=2.4, cex=.8);
mtext(side=3, text=expression(bold('   a')), line=-1, cex=.9, adj=0);


# (b) survival functions of sea + storm surge level
par(mai=c(.05,.63,.3,.2))
# stationary
plot( x.lsl, log10(apply(sf.seasurlev$rcp26$none$st,2,mean)), type='l', xlim=c(0,5), ylim=c(-2.2,0),
      lty=1, col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lwd=1.5, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i')
lines(x.lsl, log10(apply(sf.seasurlev$rcp26$gamma$st,2,mean)), type='l',
      lty=2, col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp26$uniform$st,2,mean)), type='l',
      lty=3, col=rgb(mycol[10,1],mycol[10,2],mycol[10,3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp45$none$st,2,mean)), type='l',
      lty=1, col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp45$gamma$st,2,mean)), type='l',
      lty=2, col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp45$uniform$st,2,mean)), type='l',
      lty=3, col=rgb(mycol[8,1],mycol[8,2],mycol[8,3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp85$none$st,2,mean)), type='l',
      lty=1, col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp85$gamma$st,2,mean)), type='l',
      lty=2, col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp85$uniform$st,2,mean)), type='l',
      lty=3, col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=1.5)
# non-stationary
lines(x.lsl, log10(apply(sf.seasurlev$rcp26$none$ns,2,mean)), type='l',
      lty=1, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp26$gamma$ns,2,mean)), type='l',
      lty=2, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp26$uniform$ns,2,mean)), type='l',
      lty=3, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp45$none$ns,2,mean)), type='l',
      lty=1, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp45$gamma$ns,2,mean)), type='l',
      lty=2, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp45$uniform$ns,2,mean)), type='l',
      lty=3, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp85$none$ns,2,mean)), type='l',
      lty=1, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp85$gamma$ns,2,mean)), type='l',
      lty=2, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)
lines(x.lsl, log10(apply(sf.seasurlev$rcp85$uniform$ns,2,mean)), type='l',
      lty=3, col=rgb(col85[1],col85[2],col85[3]), lwd=1.5)

# x-axis, ticks and text
axis(1,seq(0,5,0.5),lab=c("0","","1","","2","","3","","4","","5"), cex.axis=1.3)
mtext('Projected sea+surge level in 2065 [m]', side=1, line=2.6, cex=.8);
# y-axis, ticks and text
axis(2, at=seq(-3,0), label=parse(text=paste("10^", seq(-3,0), sep="")), las=1, cex.axis=1.3, mgp=c(3,.7,0))
mtext('Survival function [1-CDF]', side = 2, line=2.8, cex=0.8);
# add a panel label
mtext(side=3, text=expression(bold('   b')), line=-1.2, cex=.9, adj=0);
# add the 1/100, 1/500, 1/1000 levels
lines(c(-4,6),c(-2,-2),lty=2,col='black'); text(0.8,-1.88,"1/100 level", cex=1.3)
lines(c(-4,6),c(-log10(500),-log10(500)),lty=2,col='black'); text(0.6,-2.55,"1/500 level", cex=1.3);
#lines(c(-4,6),c(-3,-3),lty=2,col='black'); text(0.1,-2.85,"1/1000 level", cex=1.3)


# TODO
# TODO -- add legend to 3rd panel
# TODO
par(mai=c(.05,.05,.05,.05))
# stationary
plot( 9, 9, xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE, xlim=c(0,1),ylim=c(0,1), xlab='', ylab='')

legend(0.15,0.75,c("RCP2.6, stationary surge","RCP4.5, stationary surge","RCP8.5, stationary surge","RCP2.6, non-stationary surge","RCP4.5, non-stationary surge","RCP8.5, non-stationary surge","no FD","FD, gamma","FD, uniform"),
       lty=c(1,1,1,1,1,1,1,2,3), lwd=2, cex=1.2,
       col=c(rgb(mycol[10,1],mycol[10,2],mycol[10,3]),rgb(mycol[8,1],mycol[8,2],mycol[8,3]),rgb(mycol[6,1],mycol[6,2],mycol[6,3]),rgb(col26[1],col26[2],col26[3]),rgb(col45[1],col45[2],col45[3]),rgb(col85[1],col85[2],col85[3]),'black','black','black'),
       bty='n')

dev.off()

##==============================================================================









##==============================================================================

##=========   box-whisker plots of return periods with no building.
## FIGURE 3   center all around 100-year, use log10-scale.
##=========   Note: post-processing done in Google Slides.b1.

pdf(paste(plotdir,'returnperiods.pdf',sep=''),width=5,height=5.5,colormodel='cmyk')

b1 <- boxplot(ReturnPeriod~Scenario,data=return.period, horizontal=TRUE, log='x',
              xlab="", ylab="", outline=FALSE, ylim=c(8,10000), xaxt='n', yaxt='n', plot=FALSE)
tmp <- data.frame(cbind(1:18, b1$stats[3,])); colnames(tmp) <- c('scen','retp')
tmp2 <- tmp[order(tmp$retp),]
rp.order <- tmp2[,1]
return.period2 <- return.period
for (i in 1:18) {
    #return.period2[((i-1)*n.ensemble['none']+1):(i*n.ensemble['none']),1] <- rp.order[i]
    return.period2[((i-1)*n.ensemble['none']+1):(i*n.ensemble['none']),2] <- return.period[((rev(rp.order)[i]-1)*n.ensemble['none']+1):((rev(rp.order)[i])*n.ensemble['none']),2]
}
scen.names2 <- scen.names[rev(rp.order)]

par(mfrow=c(1,1), mai=c(.6,1.66,.1,.05))
boxplot(ReturnPeriod~Scenario,data=return.period2, horizontal=TRUE, log='x',
        xlab="", ylab="", outline=FALSE, ylim=c(8,75000), xaxt='n', yaxt='n', xaxs='i',xlim=c(1,21))
axis(1,c(1,10,1e2,1e3,1e4),lab=c("1","10","100","1,000","10,000"), cex.axis=1.1, mgp=c(3,.6,.0))
mtext('Return period [years]', side=1, line=1.9, cex=1.1)
axis(2,seq(1,18),lab=scen.names2, cex.axis=1.1, las=1, mgp=c(3,.6,0))
mtext('Scenario', side=3, line=1.3, cex=1.1, adj=-.3)
lines(c(100,100),c(-100,100),lty=1, lwd=1.5, col=rgb(col85[1],col85[2],col85[3]))   # Master Plan general safety
lines(c(500,500),c(-100,20.5),lty=5, lwd=1.5, col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]))   # Master Plan critical infrastructure
lines(c(5000,5000),c(-100,19),lty=6, lwd=1.5, col=rgb(mycol[3,1],mycol[3,2],mycol[3,3])) # Dutch Perspective economically-efficient

dev.off()

##==============================================================================








##==============================================================================

##==========   same as figure 3, but with additional 3 feet of heightening
## FIGURE S1
##==========

pdf(paste(plotdir,'returnperiods_build.pdf',sep=''),width=5,height=5.5,colormodel='cmyk')

b1 <- boxplot(ReturnPeriod~Scenario,data=return.period.build, horizontal=TRUE, log='x',
              xlab="", ylab="", outline=FALSE, ylim=c(8,10000), xaxt='n', yaxt='n', plot=FALSE)
tmp <- data.frame(cbind(1:18, b1$stats[3,])); colnames(tmp) <- c('scen','retp')
tmp2 <- tmp[order(tmp$retp),]
rp.order <- tmp2[,1]
return.period.build2 <- return.period.build
for (i in 1:18) {
    return.period.build2[((i-1)*n.ensemble['none']+1):(i*n.ensemble['none']),2] <- return.period.build[((rev(rp.order)[i]-1)*n.ensemble['none']+1):((rev(rp.order)[i])*n.ensemble['none']),2]
}
scen.names2 <- scen.names[rev(rp.order)]

par(mfrow=c(1,1), mai=c(.6,1.66,.1,.05))
boxplot(ReturnPeriod~Scenario,data=return.period.build2, horizontal=TRUE, log='x',
        xlab="", ylab="", outline=FALSE, ylim=c(8,75000), xaxt='n', yaxt='n', xaxs='i',xlim=c(1,21))
axis(1,c(1,10,1e2,1e3,1e4),lab=c("1","10","100","1,000","10,000"), cex.axis=1.1, mgp=c(3,.6,.0))
mtext('Return period [years]', side=1, line=1.9, cex=1.1)
axis(2,seq(1,18),lab=scen.names2, cex.axis=1.1, las=1, mgp=c(3,.6,0))
mtext('Scenario', side=3, line=1.3, cex=1.1, adj=-.3)
lines(c(100,100),c(-100,100),lty=1, lwd=1.5, col=rgb(col85[1],col85[2],col85[3]))   # Master Plan general safety
lines(c(500,500),c(-100,20.5),lty=5, lwd=1.5, col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]))   # Master Plan critical infrastructure
lines(c(5000,5000),c(-100,19),lty=6, lwd=1.5, col=rgb(mycol[3,1],mycol[3,2],mycol[3,3])) # Dutch Perspective economically-efficient

dev.off()

##==============================================================================

##==========   same as figure 3, but with structural fragility at 60% overtopping
## FIGURE S2   (set above) instead of assuming 100% overtopping
##==========

pdf(paste(plotdir,'returnperiods_fragile60.pdf',sep=''),width=5,height=5.5,colormodel='cmyk')

b1 <- boxplot(ReturnPeriod~Scenario,data=return.period.fragile60, horizontal=TRUE, log='x',
              xlab="", ylab="", outline=FALSE, ylim=c(8,10000), xaxt='n', yaxt='n', plot=FALSE)
tmp <- data.frame(cbind(1:18, b1$stats[3,])); colnames(tmp) <- c('scen','retp')
tmp2 <- tmp[order(tmp$retp),]
rp.order <- tmp2[,1]
return.period.fragile602 <- return.period.fragile60
for (i in 1:18) {
    return.period.fragile602[((i-1)*n.ensemble['none']+1):(i*n.ensemble['none']),2] <- return.period.fragile60[((rev(rp.order)[i]-1)*n.ensemble['none']+1):((rev(rp.order)[i])*n.ensemble['none']),2]
}
scen.names2 <- scen.names[rev(rp.order)]

par(mfrow=c(1,1), mai=c(.6,1.66,.1,.05))
boxplot(ReturnPeriod~Scenario,data=return.period.fragile602, horizontal=TRUE, log='x',
        xlab="", ylab="", outline=FALSE, ylim=c(8,75000), xaxt='n', yaxt='n', xaxs='i',xlim=c(1,21))
axis(1,c(1,10,1e2,1e3,1e4),lab=c("1","10","100","1,000","10,000"), cex.axis=1.1, mgp=c(3,.6,.0))
mtext('Return period [years]', side=1, line=1.9, cex=1.1)
axis(2,seq(1,18),lab=scen.names2, cex.axis=1.1, las=1, mgp=c(3,.6,0))
mtext('Scenario', side=3, line=1.3, cex=1.1, adj=-.3)
lines(c(100,100),c(-100,100),lty=1, lwd=1.5, col=rgb(col85[1],col85[2],col85[3]))   # Master Plan general safety
lines(c(500,500),c(-100,20.5),lty=5, lwd=1.5, col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]))   # Master Plan critical infrastructure
lines(c(5000,5000),c(-100,19),lty=6, lwd=1.5, col=rgb(mycol[3,1],mycol[3,2],mycol[3,3])) # Dutch Perspective economically-efficient

dev.off()

##==============================================================================

##==========   same as figure 3, but with structural fragility at 80% overtopping
## FIGURE S3   (set above) instead of assuming 100% overtopping
##==========

pdf(paste(plotdir,'returnperiods_fragile80.pdf',sep=''),width=5,height=5.5,colormodel='cmyk')

b1 <- boxplot(ReturnPeriod~Scenario,data=return.period.fragile80, horizontal=TRUE, log='x',
              xlab="", ylab="", outline=FALSE, ylim=c(8,10000), xaxt='n', yaxt='n', plot=FALSE)
tmp <- data.frame(cbind(1:18, b1$stats[3,])); colnames(tmp) <- c('scen','retp')
tmp2 <- tmp[order(tmp$retp),]
rp.order <- tmp2[,1]
return.period.fragile802 <- return.period.fragile80
for (i in 1:18) {
    return.period.fragile802[((i-1)*n.ensemble['none']+1):(i*n.ensemble['none']),2] <- return.period.fragile80[((rev(rp.order)[i]-1)*n.ensemble['none']+1):((rev(rp.order)[i])*n.ensemble['none']),2]
}
scen.names2 <- scen.names[rev(rp.order)]

par(mfrow=c(1,1), mai=c(.6,1.66,.1,.05))
boxplot(ReturnPeriod~Scenario,data=return.period.fragile802, horizontal=TRUE, log='x',
        xlab="", ylab="", outline=FALSE, ylim=c(8,75000), xaxt='n', yaxt='n', xaxs='i',xlim=c(1,21))
axis(1,c(1,10,1e2,1e3,1e4),lab=c("1","10","100","1,000","10,000"), cex.axis=1.1, mgp=c(3,.6,.0))
mtext('Return period [years]', side=1, line=1.9, cex=1.1)
axis(2,seq(1,18),lab=scen.names2, cex.axis=1.1, las=1, mgp=c(3,.6,0))
mtext('Scenario', side=3, line=1.3, cex=1.1, adj=-.3)
lines(c(100,100),c(-100,100),lty=1, lwd=1.5, col=rgb(col85[1],col85[2],col85[3]))   # Master Plan general safety
lines(c(500,500),c(-100,20.5),lty=5, lwd=1.5, col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]))   # Master Plan critical infrastructure
lines(c(5000,5000),c(-100,19),lty=6, lwd=1.5, col=rgb(mycol[3,1],mycol[3,2],mycol[3,3])) # Dutch Perspective economically-efficient

dev.off()

##==============================================================================









##==============================================================================
##==============================================================================
## end of main paper stuff. so far, nothing below this line is included AT ALL
##==============================================================================
##==============================================================================







##==============================================================================
##==============================================================================
## Scrap/scratch work:


rcp.best <- 'rcp26'; ais.best <- 'none'; ss.best <- 'st'
rcp.worst <- 'rcp85'; ais.worst <- 'uniform'; ss.worst <- 'ns'

pfail.best.prot <- 1/Ropt[[rcp.best]][[ais.best]][[ss.best]]
pfail.worst.prot <- 1/Ropt[[rcp.worst]][[ais.worst]][[ss.worst]]
pfail.best.noac <- 1/preturn[[rcp.best]][[ais.best]][[ss.best]][1,]
pfail.worst.noac <- 1/preturn[[rcp.worst]][[ais.worst]][[ss.worst]][1,]

pfail.low <- 0
pfail.high <- 0.1
pfail.n <- 2^8

f.pfail.best.noac <- density(x=pfail.best.noac, from=pfail.low, to=pfail.high, n=pfail.n)
f.pfail.worst.noac <- density(x=pfail.worst.noac, from=pfail.low, to=pfail.high, n=pfail.n)
f.pfail.best.prot <- density(x=pfail.best.prot, from=pfail.low, to=pfail.high, n=pfail.n)
f.pfail.worst.prot <- density(x=pfail.worst.prot, from=pfail.low, to=pfail.high, n=pfail.n)

# (c) pdfs of flood risk (average annual exceedance probability)
par(mai=c(.5,.35,.1,.5))
plot(f.pfail.best.noac$x, f.pfail.best.noac$y, type='l', col=mycol.rgb[2], lty=2,
     xlim=c(0,0.04), ylim=c(0,155), xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(f.pfail.worst.noac$x, f.pfail.worst.noac$y, type='l', col=mycol.rgb[11], lty=2)
lines(f.pfail.best.prot$x, f.pfail.best.prot$y, type='l', col=mycol.rgb[2], lty=1)
lines(f.pfail.worst.prot$x, f.pfail.worst.prot$y, type='l', col=mycol.rgb[11], lty=1)
xlabels <- seq(0,0.04,0.01)
axis(1,xlabels, cex.axis=1.3)
u <- par("usr")
mtext('Probability density', side=2, line=1.2, cex=.9);
mtext('Projected annual flood probability', side=1, line=2.3, cex=.9);
mtext(side=3, text=expression(bold('   c')), line=-1.4, cex=.9, adj=0);
arrows(0, u[3], 0, .97*u[4], code = 2, xpd = TRUE)

legend(0.018,150,c("best case","worst case","with protection","no action"),
       lty=c(1,1,1,2), lwd=2, cex=1.3,
       col=c(mycol.rgb[2],mycol.rgb[11],'black','black'),
       bty='n')

dev.off()
##=========

## Numbers to accompany Figure 1c (tail probabilities, above 1% flood probability)
rcp.best <- 'rcp26'; ais.best <- 'none'; ss.best <- 'st'
rcp.worst <- 'rcp85'; ais.worst <- 'uniform'; ss.worst <- 'ns'

tmp <- 1/Ropt[[rcp.best]][[ais.best]][[ss.best]]; itmp <- which(tmp>0.01); tail.best.prot <- length(itmp)/length(tmp)
tmp <- 1/Ropt[[rcp.worst]][[ais.worst]][[ss.worst]]; itmp <- which(tmp>0.01); tail.worst.prot <- length(itmp)/length(tmp)
tmp <- 1/preturn[[rcp.best]][[ais.best]][[ss.best]][1,]; itmp <- which(tmp>0.01); tail.best.noac <- length(itmp)/length(tmp)
tmp <- 1/preturn[[rcp.worst]][[ais.worst]][[ss.worst]][1,]; itmp <- which(tmp>0.01); tail.worst.noac <- length(itmp)/length(tmp)
tails <- rbind(c(tail.best.prot,tail.worst.prot), c(tail.best.noac,tail.worst.noac))
colnames(tails) <- c('best','worst')
rownames(tails) <- c('prot','noac')


rcp.tmp <- 'rcp26'; ais.tmp <- 'none'; ss.tmp <- 'st'
tmp <- 1/preturn[[rcp.tmp]][[ais.tmp]][[ss.tmp]][1,]; itmp <- which(tmp>0.01); tail.tmp <- length(itmp)/length(tmp)



##==============================================================================

pdf(paste(plotdir,'oneScen_EUM_best_vandantzig.pdf',sep=''),width=3.5,height=3.5,colormodel='cmyk')

rcp <- 'rcp26'
ais <- 'none'
ss <- 'st'
itmp <- NA

par(mfrow=c(1,1), mai=c(.61,.6,.7,.2))

plot(heightening, C.h[[rcp]][[ais]][[ss]]/1e9, type='l', xlim=c(0,5), ylim=c(7,9), lwd=2, col='blue',
     xaxs='i', yaxs='i', xlab='', ylab='', xaxt='n');
u <- par("usr")
points(Heum[[rcp]][[ais]][[ss]], Ceum[[rcp]][[ais]][[ss]]/1e9, pch=15, col='blue', cex=2)
lines(c(0,Heum[[rcp]][[ais]][[ss]]), c(Ceum[[rcp]][[ais]][[ss]],Ceum[[rcp]][[ais]][[ss]])/1e9, lty=2);
lines(c(Heum[[rcp]][[ais]][[ss]],Heum[[rcp]][[ais]][[ss]]), c(Ceum[[rcp]][[ais]][[ss]],0)/1e9, lty=2);
lines(c(Heum[[rcp]][[ais]][[ss]],Heum[[rcp]][[ais]][[ss]]), c(Ceum[[rcp]][[ais]][[ss]]/1e9, u[4]), lty=2);
mtext('Dike heightening [m]', side=1, line=2);
mtext('Expected costs [US$ billion]', side=2, line=2);

xlabels <- seq(u[1],u[2],1)
axis(1,xlabels, cex.axis=1.)
#for (i in 1:length(xlabels)) {itmp[i] <- which(heightening==xlabels[i])}
#Rlabels <- signif(R.h[[rcp]][[ais]][[ss]][itmp],2)
#axis(3,at=xlabels,labels=Rlabels, cex.axis=1.)
#mtext('Return period [years]', side=3, line=2.3, cex=1);


dev.off()

##==============================================================================

##=========
## FIGURE 3 (GOOD version, for rcp85.uniform.ns)
##=========

print(paste('the scenario with the highest expected damages (with protection) is:',names(unlist(Ceum))[which(unlist(Ceum)==max(unlist(Ceum)))]))

rcp.best <- 'rcp26'; ais.best <- 'none'; ss.best <- 'st'
rcp.worst <- 'rcp85'; ais.worst <- 'uniform'; ss.worst <- 'ns'

allcosts.noac <- NULL
allcosts.prot <- NULL
for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    for (ss in scen.ss) {
      allcosts.noac <- c(allcosts.noac, cost[[rcp]][[ais]][[ss]][1,])
      allcosts.prot <- c(allcosts.prot, Copt[[rcp]][[ais]][[ss]])
    }
  }
}
bestcosts.noac <- cost[[rcp.best]][[ais.best]][[ss.best]][1,]
worstcosts.noac <- cost[[rcp.worst]][[ais.worst]][[ss.worst]][1,]

range.noac <- signif(c(C.h[[rcp.best]][[ais.best]][[ss.best]][1],C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1])/1e9,2)
range.prot <- signif(c(Ceum[[rcp.best]][[ais.best]][[ss.best]],Ceum[[rcp.worst]][[ais.worst]][[ss.worst]])/1e9,2)

pdf(paste(plotdir,'protection_vs_noaction_ranges_vandantzig.pdf',sep=''),width=3.5,height=3.5,colormodel='cmyk')

itmp <- NA
xlabels <- seq(0,5,1)
ylabels <- seq(7,9,.5)
iend <- which(heightening==xlabels[length(xlabels)])

par(mfrow=c(1,1), mai=c(.65,.6,.1,.1))

plot(heightening[1:iend], C.h[[rcp.best]][[ais.best]][[ss.best]][1:iend]/1e9, type='l', xlim=c(0,5.3), ylim=c(7.2,9.1), lwd=2,
     col=mycol.rgb[2], xaxs='i', yaxs='i', xlab='', ylab='', xaxt='n', axes=FALSE);
lines(heightening[1:iend], C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1:iend]/1e9,
      col=mycol.rgb[11], lwd=2);
u <- par("usr")

mtext('Dike heightening [m]', side=1, line=2);
mtext('Expected costs [US$ billion]', side=2, line=2);

axis(1,xlabels, cex.axis=1.)
axis(2,ylabels, cex.axis=1.)

points(Heum[[rcp.best]][[ais.best]][[ss.best]], Ceum[[rcp.best]][[ais.best]][[ss.best]]/1e9, pch=15,
       col=mycol.rgb[2], cex=2);
lines(c(0,xlabels[length(xlabels)]),
      c(C.h[[rcp.best]][[ais.best]][[ss.best]][1],C.h[[rcp.best]][[ais.best]][[ss.best]][1])/1e9, lty=2);
lines(c(Heum[[rcp.best]][[ais.best]][[ss.best]],xlabels[length(xlabels)]),
      c(Ceum[[rcp.best]][[ais.best]][[ss.best]],Ceum[[rcp.best]][[ais.best]][[ss.best]])/1e9, lty=2);
lines(c(Heum[[rcp.best]][[ais.best]][[ss.best]],Heum[[rcp.best]][[ais.best]][[ss.best]]),
      c(Ceum[[rcp.best]][[ais.best]][[ss.best]],0)/1e9, lty=2);
polygon(c(1*xlabels[length(xlabels)],1.05*xlabels[length(xlabels)],1.05*xlabels[length(xlabels)],1*xlabels[length(xlabels)]),
        c(C.h[[rcp.best]][[ais.best]][[ss.best]][1],C.h[[rcp.best]][[ais.best]][[ss.best]][1],C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1],C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1])/1e9,
        col=mycol.rgb[10], border=NA)
text(3.5,8.49,'no-action damages:', pos=1, cex=1.)
text(3.9,8.36,paste('$',range.noac[1],'-',range.noac[2],' billion',sep=''), pos=1, cex=1.)

points(Heum[[rcp.worst]][[ais.worst]][[ss.worst]], Ceum[[rcp.worst]][[ais.worst]][[ss.worst]]/1e9, pch=15,
       col=mycol.rgb[11], cex=2);
lines(c(0,xlabels[length(xlabels)]),
      c(C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1],C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1])/1e9, lty=2);
lines(c(Heum[[rcp.worst]][[ais.worst]][[ss.worst]],xlabels[length(xlabels)]),
      c(Ceum[[rcp.worst]][[ais.worst]][[ss.worst]],Ceum[[rcp.worst]][[ais.worst]][[ss.worst]])/1e9, lty=2);
lines(c(Heum[[rcp.worst]][[ais.worst]][[ss.worst]],Heum[[rcp.worst]][[ais.worst]][[ss.worst]]),
      c(Ceum[[rcp.worst]][[ais.worst]][[ss.worst]],0)/1e9, lty=2);
polygon(c(1*xlabels[length(xlabels)],1.05*xlabels[length(xlabels)],1.05*xlabels[length(xlabels)],1*xlabels[length(xlabels)]),
        c(Ceum[[rcp.best]][[ais.best]][[ss.best]],Ceum[[rcp.best]][[ais.best]][[ss.best]],Ceum[[rcp.worst]][[ais.worst]][[ss.worst]],Ceum[[rcp.worst]][[ais.worst]][[ss.worst]])/1e9,
        col=mycol.rgb[10], border=NA)
text(4.1,7.73,paste('$',range.prot[1],'-',range.prot[2],' billion',sep=''), pos=1, cex=1.)
text(4.2,7.6,'damages with', pos=1, cex=1.)
text(4.4,7.47,'  protection', pos=1, cex=1.)

legend(-.05,9.2,c("best case","worst case"),
       lty=c(1,1), lwd=2, cex=1,
       col=c(mycol.rgb[2],mycol.rgb[11]),
       bty='n')

dev.off()

##==============================================================================

##=========
## FIGURE 3 (version for rcp85.gamma.ns)
##=========

print(paste('the scenario with the highest expected damages (with protection) is:',names(unlist(Ceum))[which(unlist(Ceum)==max(unlist(Ceum)))]))

rcp.best <- 'rcp26'; ais.best <- 'none'; ss.best <- 'st'
rcp.worst <- 'rcp85'; ais.worst <- 'gamma'; ss.worst <- 'ns'

allcosts.noac <- NULL
allcosts.prot <- NULL
for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    for (ss in scen.ss) {
      allcosts.noac <- c(allcosts.noac, cost[[rcp]][[ais]][[ss]][1,])
      allcosts.prot <- c(allcosts.prot, Copt[[rcp]][[ais]][[ss]])
    }
  }
}
bestcosts.noac <- cost[[rcp.best]][[ais.best]][[ss.best]][1,]
worstcosts.noac <- cost[[rcp.worst]][[ais.worst]][[ss.worst]][1,]

range.noac <- signif(c(C.h[[rcp.best]][[ais.best]][[ss.best]][1],C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1])/1e9,2)
range.prot <- signif(c(Ceum[[rcp.best]][[ais.best]][[ss.best]],Ceum[[rcp.worst]][[ais.worst]][[ss.worst]])/1e9,2)

#pdf(paste(plotdir,'protection_vs_noaction_ranges_vandantzig.pdf',sep=''),width=3.5,height=3.5,colormodel='cmyk')

itmp <- NA
xlabels <- seq(0,5,1)
ylabels <- seq(7,9,.5)
iend <- which(heightening==xlabels[length(xlabels)])

par(mfrow=c(1,1), mai=c(.62,.6,.1,.1))

plot(heightening[1:iend], C.h[[rcp.best]][[ais.best]][[ss.best]][1:iend]/1e9, type='l', xlim=c(0,5.3), ylim=c(7.2,9.1), lwd=2,
     col=mycol.rgb[2], xaxs='i', yaxs='i', xlab='', ylab='', xaxt='n', axes=FALSE);
lines(heightening[1:iend], C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1:iend]/1e9,
      col=mycol.rgb[11], lwd=2);
u <- par("usr")

mtext('Dike heightening [m]', side=1, line=2);
mtext('Expected costs [US$ billion]', side=2, line=2);

axis(1,xlabels, cex.axis=1.)
axis(2,ylabels, cex.axis=1.)

points(Heum[[rcp.best]][[ais.best]][[ss.best]], Ceum[[rcp.best]][[ais.best]][[ss.best]]/1e9, pch=15,
       col=mycol.rgb[2], cex=2);
lines(c(0,xlabels[length(xlabels)]),
      c(C.h[[rcp.best]][[ais.best]][[ss.best]][1],C.h[[rcp.best]][[ais.best]][[ss.best]][1])/1e9, lty=2);
lines(c(Heum[[rcp.best]][[ais.best]][[ss.best]],xlabels[length(xlabels)]),
      c(Ceum[[rcp.best]][[ais.best]][[ss.best]],Ceum[[rcp.best]][[ais.best]][[ss.best]])/1e9, lty=2);
lines(c(Heum[[rcp.best]][[ais.best]][[ss.best]],Heum[[rcp.best]][[ais.best]][[ss.best]]),
      c(Ceum[[rcp.best]][[ais.best]][[ss.best]],0)/1e9, lty=2);
polygon(c(1*xlabels[length(xlabels)],1.05*xlabels[length(xlabels)],1.05*xlabels[length(xlabels)],1*xlabels[length(xlabels)]),
				c(C.h[[rcp.best]][[ais.best]][[ss.best]][1],C.h[[rcp.best]][[ais.best]][[ss.best]][1],C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1],C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1])/1e9,
        col=mycol.rgb[10], border=NA)
text(3.45,8.68,'no-action damages:', pos=1, cex=1.)
text(3.85,8.55,paste('$',range.noac[1],'-',range.noac[2],' billion',sep=''), pos=1, cex=1.)

points(Heum[[rcp.worst]][[ais.worst]][[ss.worst]], Ceum[[rcp.worst]][[ais.worst]][[ss.worst]]/1e9, pch=15,
       col=mycol.rgb[11], cex=2);
lines(c(0,xlabels[length(xlabels)]),
      c(C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1],C.h[[rcp.worst]][[ais.worst]][[ss.worst]][1])/1e9, lty=2);
lines(c(Heum[[rcp.worst]][[ais.worst]][[ss.worst]],xlabels[length(xlabels)]),
      c(Ceum[[rcp.worst]][[ais.worst]][[ss.worst]],Ceum[[rcp.worst]][[ais.worst]][[ss.worst]])/1e9, lty=2);
lines(c(Heum[[rcp.worst]][[ais.worst]][[ss.worst]],Heum[[rcp.worst]][[ais.worst]][[ss.worst]]),
      c(Ceum[[rcp.worst]][[ais.worst]][[ss.worst]],0)/1e9, lty=2);
polygon(c(1*xlabels[length(xlabels)],1.05*xlabels[length(xlabels)],1.05*xlabels[length(xlabels)],1*xlabels[length(xlabels)]),
				c(Ceum[[rcp.best]][[ais.best]][[ss.best]],Ceum[[rcp.best]][[ais.best]][[ss.best]],Ceum[[rcp.worst]][[ais.worst]][[ss.worst]],Ceum[[rcp.worst]][[ais.worst]][[ss.worst]])/1e9,
        col=mycol.rgb[10], border=NA)
text(4.2,7.67,paste('$',range.prot[1],'-',range.prot[2],' billion',sep=''), pos=1, cex=1.)
text(4.3,7.54,'damages with', pos=1, cex=1.)
text(4.5,7.41,'  protection', pos=1, cex=1.)

legend(.1,9.2,c("best case","worst case"),
       lty=c(1,1), lwd=2, cex=1,
       col=c(mycol.rgb[2],mycol.rgb[11]),
       bty='n')

dev.off()

##==============================================================================

##=========
## FIGURE 3 (version with return period)
##=========

#pdf(paste(plotdir,'oneScen_EUM_worst_vandantzig.pdf',sep=''),width=3.5,height=3.5,colormodel='cmyk')

rcp <- 'rcp85'
ais <- 'gamma'
ss <- 'ns'
itmp <- NA

par(mfrow=c(1,1), mai=c(.61,.6,.7,.2))

plot(heightening, C.h[[rcp]][[ais]][[ss]]/1e9, type='l', xlim=c(0,5), ylim=c(7,9), lwd=2, col='red',
     xaxs='i', yaxs='i', xlab='', ylab='', xaxt='n');
u <- par("usr")

points(Heum[[rcp]][[ais]][[ss]], Ceum[[rcp]][[ais]][[ss]]/1e9, pch=15, col='red', cex=2);
lines(c(0,Heum[[rcp]][[ais]][[ss]]), c(Ceum[[rcp]][[ais]][[ss]],Ceum[[rcp]][[ais]][[ss]])/1e9, lty=2);
lines(c(Heum[[rcp]][[ais]][[ss]],Heum[[rcp]][[ais]][[ss]]), c(Ceum[[rcp]][[ais]][[ss]],0)/1e9, lty=2);
lines(c(Heum[[rcp]][[ais]][[ss]],Heum[[rcp]][[ais]][[ss]]), c(Ceum[[rcp]][[ais]][[ss]]/1e9, u[4]), lty=2);

mtext('Dike heightening [m]', side=1, line=2);
mtext('Expected costs [US$ billion]', side=2, line=2);

#xlabels <- seq(u[1],u[2],1)
#for (i in 1:length(xlabels)) {itmp[i] <- which(heightening==xlabels[i])}
#Rlabels <- signif(R.h[[rcp]][[ais]][[ss]][itmp],2)
axis(1,xlabels, cex.axis=1.)
axis(3,at=xlabels,labels=Rlabels, cex.axis=1.)
mtext('Return period [years]', side=3, line=2.3, cex=1);

dev.off()

##==============================================================================

##=========
## FIGURE 3 (version with all scenarios)
##=========

#pdf(paste(plotdir,'allScen_EUM_vandantzig.pdf',sep=''),width=3.5,height=3.5,colormodel='cmyk')

itmp <- NA
col.rcp <- t(data.frame(c('green','blue','red'))); colnames(col.rcp) <- scen.rcp; rownames(col.rcp) <- NULL

par(mfrow=c(1,1), mai=c(.61,.6,.7,.2))

plot(heightening, C.h$rcp26$none$st/1e9, type='l', xlim=c(0,5), ylim=c(7,9), lwd=2, col='black',
     xaxs='i', yaxs='i', xlab='', ylab='', xaxt='n');
for (rcp in scen.rcp) {
    for (ais in scen.ais) {
        for (ss in scen.ss) {
            lines(heightening, C.h[[rcp]][[ais]][[ss]]/1e9, col=col.rcp[1,rcp], lwd=2);
            points(Heum[[rcp]][[ais]][[ss]], Ceum[[rcp]][[ais]][[ss]]/1e9, pch=15, col=col.rcp[1,rcp], cex=1.6);
            lines(c(0,Heum[[rcp]][[ais]][[ss]]), c(Ceum[[rcp]][[ais]][[ss]],Ceum[[rcp]][[ais]][[ss]])/1e9, lty=2);
            lines(c(Heum[[rcp]][[ais]][[ss]],Heum[[rcp]][[ais]][[ss]]), c(0, u[4]), lty=2);
        }
    }
}
u <- par("usr")
mtext('Dike heightening [m]', side=1, line=2);
mtext('Expected costs [US$ billion]', side=2, line=2);
xlabels <- seq(u[1],u[2],1)
for (i in 1:length(xlabels)) {itmp[i] <- which(heightening==xlabels[i])}
Rlabels <- signif(R.h[[rcp]][[ais]][[ss]][itmp],2)
axis(1,xlabels, cex.axis=1.)
axis(3,at=xlabels,labels=Rlabels, cex.axis=1.)
mtext('Return period [years]', side=3, line=2.3, cex=1);

legend(2,6,c("RCP2.6","RCP4.5","RCP8.5"),
       lty=c(1,1,1), lwd=2, cex=1,
       col=c('green','blue','red'),
       bty='n')

dev.off()

##==============================================================================
##==============================================================================

## TABLE -- EUM-optimal, in the mean sense

H.rcp26 <- data.frame(rbind(c(Heum$rcp26$none$st Heum$rcp26$gamma$st Heum$rcp26$uniform$st),
                            c(Heum$rcp26$none$ns Heum$rcp26$gamma$ns Heum$rcp26$uniform$ns)))
H.rcp45 <- data.frame(rbind(c(Heum$rcp45$none$st Heum$rcp45$gamma$st Heum$rcp45$uniform$st),
                            c(Heum$rcp45$none$ns Heum$rcp45$gamma$ns Heum$rcp45$uniform$ns)))
H.rcp85 <- data.frame(rbind(c(Heum$rcp85$none$st Heum$rcp85$gamma$st Heum$rcp85$uniform$st),
                            c(Heum$rcp85$none$ns Heum$rcp85$gamma$ns Heum$rcp85$uniform$ns)))

C.rcp26 <- data.frame(rbind(c(Ceum$rcp26$none$st Ceum$rcp26$gamma$st Ceum$rcp26$uniform$st),
                            c(Ceum$rcp26$none$ns Ceum$rcp26$gamma$ns Ceum$rcp26$uniform$ns)))/1e9
C.rcp45 <- data.frame(rbind(c(Ceum$rcp45$none$st Ceum$rcp45$gamma$st Ceum$rcp45$uniform$st),
                            c(Ceum$rcp45$none$ns Ceum$rcp45$gamma$ns Ceum$rcp45$uniform$ns)))/1e9
C.rcp85 <- data.frame(rbind(c(Ceum$rcp85$none$st Ceum$rcp85$gamma$st Ceum$rcp85$uniform$st),
                            c(Ceum$rcp85$none$ns Ceum$rcp85$gamma$ns Ceum$rcp85$uniform$ns)))/1e9

C.min <- min(unlist(Ceum))

## Cost table

table.rcp26.no <- 100*(C.rcp26.no - C.min.no)/C.min.no
colnames(table.rcp26.no) <- scen.ais
rownames(table.rcp26.no) <- c("RCP26, stat, no VOSL","RCP26, nonstat, no VOSL")

table.rcp45.no <- 100*(C.rcp45.no - C.min.no)/C.min.no
colnames(table.rcp45.no) <- scen.ais
rownames(table.rcp45.no) <- c("RCP45, stat, no VOSL","RCP45, nonstat, no VOSL")

table.rcp85.no <- 100*(C.rcp85.no - C.min.no)/C.min.no
colnames(table.rcp85.no) <- scen.ais
rownames(table.rcp85.no) <- c("RCP85, stat, no VOSL","RCP85, nonstat, no VOSL")

table.rcp26.yes <- 100*(C.rcp26.yes - C.min.yes)/C.min.yes
colnames(table.rcp26.yes) <- scen.ais
rownames(table.rcp26.yes) <- c("RCP26, stat, VOSL","RCP26, nonstat, VOSL")

table.rcp45.yes <- 100*(C.rcp45.yes - C.min.yes)/C.min.yes
colnames(table.rcp45.yes) <- scen.ais
rownames(table.rcp45.yes) <- c("RCP45, stat, VOSL","RCP45, nonstat, VOSL")

table.rcp85.yes <- 100*(C.rcp85.yes - C.min.yes)/C.min.yes
colnames(table.rcp85.yes) <- scen.ais
rownames(table.rcp85.yes) <- c("RCP85, stat, VOSL","RCP85, nonstat, VOSL")

table.cost <- rbind(c("ensemble mean total cost percent increase from best-case","",""),table.rcp26.no, c("","",""), table.rcp45.no, c("","",""), table.rcp85.no,c("","",""),table.rcp26.yes, c("","",""), table.rcp45.yes, c("","",""), table.rcp85.yes)


## Return period table

table.rcp26.no <- R.rcp26.no
colnames(table.rcp26.no) <- scen.ais
rownames(table.rcp26.no) <- c("RCP26, stat, no VOSL","RCP26, nonstat, no VOSL")

table.rcp45.no <- R.rcp45.no
colnames(table.rcp45.no) <- scen.ais
rownames(table.rcp45.no) <- c("RCP45, stat, no VOSL","RCP45, nonstat, no VOSL")

table.rcp85.no <- R.rcp85.no
colnames(table.rcp85.no) <- scen.ais
rownames(table.rcp85.no) <- c("RCP85, stat, no VOSL","RCP85, nonstat, no VOSL")

table.rcp26.yes <- R.rcp26.yes
colnames(table.rcp26.yes) <- scen.ais
rownames(table.rcp26.yes) <- c("RCP26, stat, VOSL","RCP26, nonstat, VOSL")

table.rcp45.yes <- R.rcp45.yes
colnames(table.rcp45.yes) <- scen.ais
rownames(table.rcp45.yes) <- c("RCP45, stat, VOSL","RCP45, nonstat, VOSL")

table.rcp85.yes <- R.rcp85.yes
colnames(table.rcp85.yes) <- scen.ais
rownames(table.rcp85.yes) <- c("RCP85, stat, VOSL","RCP85, nonstat, VOSL")

table.returnperiod <- rbind(c("ensemble mean return period","",""),table.rcp26.no, c("","",""), table.rcp45.no, c("","",""), table.rcp85.no,c("","",""),table.rcp26.yes, c("","",""), table.rcp45.yes, c("","",""), table.rcp85.yes)


## Heightening table

table.rcp26.no <- H.rcp26.no
colnames(table.rcp26.no) <- scen.ais
rownames(table.rcp26.no) <- c("RCP26, stat, no VOSL","RCP26, nonstat, no VOSL")

table.rcp45.no <- H.rcp45.no
colnames(table.rcp45.no) <- scen.ais
rownames(table.rcp45.no) <- c("RCP45, stat, no VOSL","RCP45, nonstat, no VOSL")

table.rcp85.no <- H.rcp85.no
colnames(table.rcp85.no) <- scen.ais
rownames(table.rcp85.no) <- c("RCP85, stat, no VOSL","RCP85, nonstat, no VOSL")

table.rcp26.yes <- H.rcp26.yes
colnames(table.rcp26.yes) <- scen.ais
rownames(table.rcp26.yes) <- c("RCP26, stat, VOSL","RCP26, nonstat, VOSL")

table.rcp45.yes <- H.rcp45.yes
colnames(table.rcp45.yes) <- scen.ais
rownames(table.rcp45.yes) <- c("RCP45, stat, VOSL","RCP45, nonstat, VOSL")

table.rcp85.yes <- H.rcp85.yes
colnames(table.rcp85.yes) <- scen.ais
rownames(table.rcp85.yes) <- c("RCP85, stat, VOSL","RCP85, nonstat, VOSL")

table.heightening <- rbind(c("ensemble mean heightening","",""),table.rcp26.no, c("","",""), table.rcp45.no, c("","",""), table.rcp85.no,c("","",""),table.rcp26.yes, c("","",""), table.rcp45.yes, c("","",""), table.rcp85.yes)

table.mean <- cbind(table.cost,rep("",18),table.heightening,rep("",18),table.returnperiod)
write.csv(table.mean, "../output_model/scenarios_mean.csv")

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================

cost.avg.st <- apply(cost$rcp26$none$st, 1, mean )
cost.avg.ns <- apply(cost$rcp26$none$ns, 1, mean )
cost.avg.fd.st <- apply(cost$rcp26$gamma$st, 1, mean )
cost.avg.fd.ns <- apply(cost$rcp26$gamma$ns, 1, mean )

plot(heightening, cost.avg.ns, type='l', xlim=c(0,3), ylim=c(3.5e9, 4.5e9))
  lines(heightening, cost.avg.st, col='red')
  lines(heightening, cost.avg.fd.ns, col='blue')
  lines(heightening, cost.avg.fd.st, col='purple')

# at the moment, looks like the formulation using the stationary GEV is the
# main driver of costs.

# in rcp2.6, note that purple and red (stationary cases) are about same trade-offs
# but green and blue are higher cost - surge is main driver



cost.avg.st <- apply(cost$rcp85$none$st, 1, mean )
cost.avg.ns <- apply(cost$rcp85$none$ns, 1, mean )
cost.avg.fd.st <- apply(cost$rcp85$gamma$st, 1, mean )
cost.avg.fd.ns <- apply(cost$rcp85$gamma$ns, 1, mean )

plot(heightening, cost.avg.ns, type='l', xlim=c(0,3), ylim=c(3.5e9, 4.5e9))
  lines(heightening, cost.avg.st, col='red')
  lines(heightening, cost.avg.fd.ns, col='blue')
  lines(heightening, cost.avg.fd.st, col='purple')




# SCRATCH
if(FALSE) {
# vvv
rcp <- 'rcp85'
ais <- 'gamma'
ss <- 'ns'
f.seasub.sur <- mat.or.vec(n.ensemble['none'], 1292)
f.surrise.sur <- mat.or.vec(n.ensemble['none'], 1292)

for (sow in 1:n.ensemble['none']) {
  f.seasub.sur[sow,] <- approx(x=x.lsl, y=f.seasublev[[rcp]][[ais]], xout=fit.surge[sow,2:1293]/1000)$y
  f.seasub.sur[sow,which(is.na(f.seasub.sur[sow,]))] <- 0
  f.seasub.sur[sow,] <- f.seasub.sur[sow,]/sintegral(x=fit.surge[sow,2:1293]/1000, fx=f.seasub.sur[sow,])$value
  f.surrise.sur[sow,] <- approx(x=x.lsl, y=f.surge.rise[[rcp]][[ais]], xout=fit.surge[sow,2:1293]/1000)$y
  f.surrise.sur[sow,which(is.na(f.surrise.sur[sow,]))] <- 0
  f.surrise.sur[sow,] <- f.surrise.sur[sow,]/sintegral(x=fit.surge[sow,2:1293]/1000, fx=f.surrise.sur[sow,])$value
}

x.tmp <- seq(from=fit.surge[sow,2], to=fit.surge[sow,1293], length.out=1292)/1000
f.seasub.tmp <- approx(x=fit.surge[sow,2:1293]/1000, y=f.seasub.sur[sow,], xout=x.tmp)$y
f.surrise.tmp <- approx(x=fit.surge[sow,2:1293]/1000, y=f.surrise.sur[sow,], xout=x.tmp)$y
f.tmp.tmp <- approx(x=fit.surge[sow,2:1293]/1000, y=f.tmp[sow,2:1293], xout=x.tmp)$y

tmp <- convolve(x=f.surrise.tmp, y=rev(f.tmp.tmp), type='open')
dx <- median(diff(x.tmp))
x.tmp2 <- seq(from=x.tmp[1], by=dx, length.out=length(tmp))
tmp2 <- tmp/sintegral(x=x.tmp2, fx=tmp)$value

x.tmp <- seq(from=fit.surge[sow,2], to=fit.surge[sow,1293], length.out=1292)/1000

tmp3 <- convolve(x=f.seasub.tmp, y=rev(tmp2), type='open')
x.tmp3 <- seq(from=x.tmp[1], by=dx, length.out=length(tmp3))
tmp3 <- tmp3/sintegral(x=x.tmp3, fx=tmp3)$value

# TODO  -



for (sow in 1:n.ensemble['none']) {
  for (node in 2:(n.ensemble['none']+1)) {
    f.seasub.sur[sow,node-1] <- (q[node]-q[node-1])/(fit.surge[sow,node]-fit.surge[sow,node-1])
  }
}
# ^^^
}
# SCRATCH




# check the 1:1293 (N.ensemble) level surge for all the stationary GEV sets
surgelevs <- rep(NA,n.ensemble['none'])
for (i in 1:n.ensemble['none']) {
  surgelevs[i] <- 0.001*qgev(1-1/n.ensemble['none'], xi=gev.stat[i,'shape'], mu=gev.stat[i,'location'], beta=gev.stat[i,'scale'])
}


rp <- init

for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    rp[[rcp]][[ais]] <- vector("list",2); names(rp[[rcp]][[ais]])=scen.ss
    for (ss in scen.ss) {
      rp[[rcp]][[ais]][[ss]] <- rep(NA, n.ensemble[[ais]])
      for (sow in 1:n.ensemble[[ais]]) {
        if(ss=='st') {h.eff <- h0 - lsl.norm[[rcp]][[ais]][sow]}
        if(ss=='ns') {h.eff <- h0 - (1+surge.factor[sow])*lsl.norm[[rcp]][[ais]][sow]}
        rp[[rcp]][[ais]][[ss]][sow] <-1/(1-pgev(q=1000*h.eff, xi=gev.stat[i,'shape'], mu=gev.stat[i,'location'], beta=gev.stat[i,'scale']))
      }
    }
  }
}

##==============================================================================

##==============================================================================
## End
##==============================================================================
