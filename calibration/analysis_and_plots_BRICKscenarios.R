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
filename.brick.nofd = '../output_model/BRICK-model_physical_fd-gamma_08Dec2016.nc'#no-FD case uses gamma, without disintegration
filename.brick.uniform = '../output_model/BRICK-model_physical_fd-uniform_08Dec2016.nc'
filename.brick.gamma = '../output_model/BRICK-model_physical_fd-gamma_08Dec2016.nc'
filename.brick.allslr = '../output_model/BRICK-model_physical_allslr_09Jan2017.nc'

## File name for the Van Dantzig model output (netCDF4)
## Each of these also has x3 RCP scenarios, x2 storm surge scenarios
filename.vandantzig.nofd = '../output_model/VanDantzig_fd-none_2065_09Jan2017.nc'
filename.vandantzig.uniform = '../output_model/VanDantzig_fd-uniform_2065_09Jan2017.nc'
filename.vandantzig.gamma = '../output_model/VanDantzig_fd-gamma_2065_09Jan2017.nc'

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

# second level are AIS fast dynamics scenarios
scen.ais = c('none','gamma','uniform')
scen.ss = c('st','ns')
for (rcp in scen.rcp) {
  cost[[rcp]] = vector("list",3); names(cost[[rcp]])=scen.ais
  loss[[rcp]] = vector("list",3); names(loss[[rcp]])=scen.ais
  investment[[rcp]] = vector("list",3); names(investment[[rcp]])=scen.ais
  preturn[[rcp]] = vector("list",3); names(preturn[[rcp]])=scen.ais
  lsl[[rcp]] = vector("list",3); names(lsl[[rcp]])=scen.ais

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

## Laplace-averaging minimization of expected costs, expected regret (all
## scenarios, equal weight)

#todo?


# experiment with 50y vs 100y time horizon?
#costs <- as.matrix(cbind(signif(unlist(Ceum)/1e9,6)))
#p.inc.cost.ss <- 100*(costs[c(2,4,6,8,10,12),]-costs[c(1,3,5,7,9,11),])/costs[c(1,3,5,7,9,11),]
#p.inc.cost.fd <- 100*(costs[c(3,4,7,8,11,12),]-costs[c(1,2,5,6,9,10),])/costs[c(1,2,5,6,9,10),]
#p.inc.cost.rcp <- 100*(costs[c(5,6,7,8,9,10,11,12),]-costs[c(1,2,3,4,1,2,3,4),])/costs[c(1,2,3,4,1,2,3,4),]
#inc.cost.ss <- (costs[c(2,4,6,8,10,12),]-costs[c(1,3,5,7,9,11),])
#inc.cost.fd <- (costs[c(3,4,7,8,11,12),]-costs[c(1,2,5,6,9,10),])
#inc.cost.rcp <- (costs[c(5,6,7,8,9,10,11,12),]-costs[c(1,2,3,4,1,2,3,4),])



##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## FIGURE -- pdfs of (a) local sea-level rise, (b) storm surge, (c) flood risk
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
f.sealev <- lsl
f.surlev <- lsl
f.surge.rise <- lsl
f.seasurlev <- lsl
f.flood <- lsl

lsl.lower <- 0
lsl.upper <- 10
lsl.n <- 2^10 # powers of 2 are good for density estimation in R
              # 2^10 = 1024 is first one greater than number of ensemble members

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
f.surge.avg <- vector("list",length(scen.rcp)); names(f.surge.avg) <- scen.rcp
for (rcp in scen.rcp) {
  f.surge.stat[[rcp]] <- vector("list",length(scen.ais)); names(f.surge.stat[[rcp]]) <- scen.ais
  f.surge.avg[[rcp]] <- vector("list",length(scen.ais)); names(f.surge.avg[[rcp]]) <- scen.ais
  for (ais in scen.ais) {
    f.surlev[[rcp]][[ais]] <- vector("list",2); names(f.surlev[[rcp]][[ais]])=scen.ss
    f.seasurlev[[rcp]][[ais]] <- vector("list",2); names(f.seasurlev[[rcp]][[ais]])=scen.ss
    f.flood[[rcp]][[ais]] <- vector("list",2); names(f.flood[[rcp]][[ais]])=scen.ss
    f.surge.stat[[rcp]][[ais]] <- mat.or.vec(n.ensemble[[ais]],flood.n)
  }
}

## distribution of sea-level and storm surges
iproj <- which(mod.time==2065)
inorm <- which(mod.time==2015)
for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    lsl[[rcp]][[ais]] <- lsl[[rcp]][[ais]][iproj,] - lsl[[rcp]][[ais]][inorm,]
    f.sealev[[rcp]][[ais]] <- density(x=lsl[[rcp]][[ais]], from=lsl.lower, to=lsl.upper, n=lsl.n)
    f.surge.rise[[rcp]][[ais]] <- density(x=lsl[[rcp]][[ais]]*surge.factor, from=lsl.lower, to=lsl.upper, n=lsl.n)
    # and initialize storm surge
    for (ss in scen.ss) {
      f.surlev[[rcp]][[ais]][[ss]] <- rep(NA, lsl.n)
    }
  }
}

x.lsl <- f.sealev$rcp85$none$x
dx.lsl <- mean(diff(x.lsl))

## distribution of storm surge level
## convolution of surge.rise distribution + gev.stat distribution

## stationary case (x*1000, because GEV parameters assume levels are in mm, not m)
#f.surge.stat <- as.numeric(dgev(x=x.lsl*1000, xi=gev.stat[1,'shape'], mu=gev.stat[1,'location'], beta=gev.stat[1,'scale'], log=FALSE))
#f.surge.stat <- f.surge.stat/sum( dx.lsl*f.surge.stat)

for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    for (i in 1:n.ensemble[[ais]]) {
      f.surge.stat[[rcp]][[ais]][i,] <- as.numeric(dgev(x=x.lsl*1000,
                                                        xi=gev.stat[i,'shape'],
                                                        mu=gev.stat[i,'location'],
                                                        beta=gev.stat[i,'scale'],
                                                        log=FALSE))
      f.surge.stat[[rcp]][[ais]][i,] <- f.surge.stat[[rcp]][[ais]][i,]/sum( dx.lsl*f.surge.stat[[rcp]][[ais]][i,])
    }

# need to add an index for each ensemble member

    # one option: average the probabilities across all ensemble members, and re-normalize
    f.surge.avg[[rcp]][[ais]] <- apply(f.surge.stat[[rcp]][[ais]], 2, mean)
    f.surge.avg[[rcp]][[ais]] <- f.surge.avg[[rcp]][[ais]]/sum(dx.lsl*f.surge.avg[[rcp]][[ais]])

    for (ss in scen.ss) {
      if(ss=='st') {
        f.surlev[[rcp]][[ais]][[ss]] <- f.surge.avg[[rcp]][[ais]]
        f.seasurlev[[rcp]][[ais]][[ss]] <- convolve(x=f.surlev[[rcp]][[ais]][[ss]], y=rev(f.sealev[[rcp]][[ais]]$y))
        f.seasurlev[[rcp]][[ais]][[ss]] <- f.seasurlev[[rcp]][[ais]][[ss]]/sum(dx.lsl*f.seasurlev[[rcp]][[ais]][[ss]])
      } else {
        f.surlev[[rcp]][[ais]][[ss]] <- convolve(x=f.surge.rise[[rcp]][[ais]]$y, y=rev(f.surge.avg[[rcp]][[ais]]))
        f.surlev[[rcp]][[ais]][[ss]] <- f.surlev[[rcp]][[ais]][[ss]]/sum(dx.lsl*f.surlev[[rcp]][[ais]][[ss]])
        f.seasurlev[[rcp]][[ais]][[ss]] <- convolve(x=f.surlev[[rcp]][[ais]][[ss]], y=rev(f.sealev[[rcp]][[ais]]$y))
        f.seasurlev[[rcp]][[ais]][[ss]] <- f.seasurlev[[rcp]][[ais]][[ss]]/sum(dx.lsl*f.seasurlev[[rcp]][[ais]][[ss]])
      }
    }
  }
}



##=========
## FIGURE 1
##=========
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





pdf(paste(plotdir,'pdfs_slr-ss-risk.pdf',sep=''),width=3.5,height=8,colormodel='cmyk')

par(mfrow=c(3,1), mai=c(.5,.35,.05,.5))
# (a) pdfs of local sea-level rise (2065)
plot(x.lsl, f.sealev$rcp26$none$y, type='l', xlim=c(0,1.1), ylim=c(0,10), lwd=1.5,
     col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=1,
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(x.lsl, f.sealev$rcp26$gamma$y, type='l', col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=2)
lines(x.lsl, f.sealev$rcp26$uniform$y, type='l', col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=3)
lines(x.lsl, f.sealev$rcp45$none$y, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=1)
lines(x.lsl, f.sealev$rcp45$gamma$y, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=2)
lines(x.lsl, f.sealev$rcp45$uniform$y, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=3)
lines(x.lsl, f.sealev$rcp85$none$y, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=1)
lines(x.lsl, f.sealev$rcp85$gamma$y, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=2)
lines(x.lsl, f.sealev$rcp85$uniform$y, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=3)

axis(1,seq(0,1,0.2),lab=c("0","0.2","0.4","0.6","0.8","1"), cex.axis=1.3)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, xpd = TRUE)
mtext('Probability density', side=2, line=1.2, cex=.9);
mtext('Projected sea level in 2065 relative to 2015 [m]', side=1, line=2.3, cex=.9);
mtext(side=3, text=expression(bold('   a')), line=-1, cex=.9, adj=0);

legend(0.5,9.5,c("RCP2.6","RCP4.5","RCP8.5","no FD","FD, gamma","FD, uniform"),
       lty=c(1,1,1,1,2,3), lwd=2, cex=1.2,
       col=c(rgb(mycol[13,1],mycol[13,2],mycol[13,3]),rgb(mycol[2,1],mycol[2,2],mycol[2,3]),rgb(mycol[11,1],mycol[11,2],mycol[11,3]),'black','black','black'),
       bty='n')

# (b) pdfs of storm surge
par(mai=c(.45,.35,.15,.5))
plot(x.lsl, f.surlev$rcp26$none$st, type='l', xlim=c(0,3), ylim=c(0,4), lwd=1.5,
     col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lty=1,
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(x.lsl, f.surlev$rcp26$none$ns, type='l', col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=2)
lines(x.lsl, f.surlev$rcp26$gamma$ns, type='l', col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=2)
lines(x.lsl, f.surlev$rcp26$uniform$ns, type='l', col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=3)
lines(x.lsl, f.surlev$rcp45$none$ns, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=1)
lines(x.lsl, f.surlev$rcp45$gamma$ns, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=2)
lines(x.lsl, f.surlev$rcp45$uniform$ns, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=3)
lines(x.lsl, f.surlev$rcp85$none$ns, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=1)
lines(x.lsl, f.surlev$rcp85$gamma$ns, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=2)
lines(x.lsl, f.surlev$rcp85$uniform$ns, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=3)

axis(1,seq(0,3,0.5),lab=c("0","0.5","1","1.5","2","2.5","3"), cex.axis=1.3)
u <- par("usr")
arrows(0, u[3],0, .95*u[4], code = 2, xpd = TRUE)
mtext('Probability density', side=2, line=1.2, cex=.9);
mtext('Projected surge level in 2065 relative to 2015 [m]', side=1, line=2.3, cex=.9);
mtext(side=3, text=expression(bold('   b')), line=-1, cex=.9, adj=0);
text(.53,3.7,"stationary", pos=4, cex=1.3)
text(1.01,2.5,"non-stationary", pos=4, cex=1.3)

#legend(0.7,4,c("stationary","all others: non-stationary\n (see panel (a) legend)"),
#       lty=c(1,NA), lwd=2, cex=1.2,
#       col=c(rgb(mycol[6,1],mycol[6,2],mycol[6,3]), NA),
#       bty='n')

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









##=========
## FIGURE 2
##=========

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

##=========













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

##=========








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

##=========







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

##=========








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

##=========










##==============================================================================
##==============================================================================



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
## end of main paper stuff. so far, nothing below this line is included AT ALL
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







##==============================================================================
## End
##==============================================================================
