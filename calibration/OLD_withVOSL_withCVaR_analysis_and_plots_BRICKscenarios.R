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
filename.brick.nofd = '../output_model/BRICK-model_physical_fd-none_01Dec2016.nc'
filename.brick.uniform = '../output_model/BRICK-model_physical_fd-uniform_01Dec2016.nc'
filename.brick.gamma = '../output_model/BRICK-model_physical_fd-gamma_01Dec2016.nc'

## File name for the Van Dantzig model output (netCDF4)
## Each of these also has x3 RCP scenarios, x2 storm surge scenarios
filename.vandantzig.nofd = '../output_model/VanDantzig_fd-none_05Dec2016.nc'
filename.vandantzig.uniform = '../output_model/VanDantzig_fd-uniform_05Dec2016.nc'
filename.vandantzig.gamma = '../output_model/VanDantzig_fd-gamma_05Dec2016.nc'

## File name for the BRICK post-calibrated parameters (netcdf) (the BRICK output came from these guys)
filename.parameters.nofd = '../output_calibration/BRICK-model_postcalibratedParameters_fd-none_01Dec2016.nc'
filename.parameters.uniform = '../output_calibration/BRICK-model_postcalibratedParameters_fd-uniform_01Dec2016.nc'
filename.parameters.gamma = '../output_calibration/BRICK-model_drawcalibratedParameters_fd-gamma_01Dec2016.nc'

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
scen.rcp = c('rcp26','rcp45','rcp85')
cost = vector("list",3); names(cost)=scen.rcp
loss = vector("list",3); names(loss)=scen.rcp
investment = vector("list",3); names(investment)=scen.rcp
preturn = vector("list",3); names(preturn)=scen.rcp
lsl = vector("list",3); names(lsl)=scen.rcp

# second level are AIS fast dynamics scenarios
scen.ais = c('none','gamma','uniform'); surge.factor = vector("list",length(scen.ais)); names(surge.factor)=scen.ais
scen.ss = c('st','ns')
scen.vosl = c('no','yes')
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

# fourth level is whether or not value of lives lost is taken into account
    for (ss in scen.ss) {
      cost[[rcp]][[ais]][[ss]] = vector("list",2); names(cost[[rcp]][[ais]][[ss]])=scen.vosl
      loss[[rcp]][[ais]][[ss]] = vector("list",2); names(loss[[rcp]][[ais]][[ss]])=scen.vosl
      investment[[rcp]][[ais]][[ss]] = vector("list",2); names(investment[[rcp]][[ais]][[ss]])=scen.vosl
      preturn[[rcp]][[ais]][[ss]] = vector("list",2); names(preturn[[rcp]][[ais]][[ss]])=scen.vosl
    }
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
  # no value of life
  heightening <- ncvar_get(ncdata, 'H')
  gev.stat <- ncvar_get(ncdata, 'gev_stat'); gev.stat <- t(as.data.frame(gev.stat)); colnames(gev.stat) <- c('location','scale','shape')
  surge.factor$none <- ncvar_get(ncdata, 'surge_factor_RCP26') # all RCPs have same surge factor

  cost$rcp26$none$st$no <- ncvar_get(ncdata, 'ExpectedCost_RCP26_stat')
  loss$rcp26$none$st$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_stat')
  investment$rcp26$none$st$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_stat')
  preturn$rcp26$none$st$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat')

  cost$rcp45$none$st$no <- ncvar_get(ncdata, 'ExpectedCost_RCP45_stat')
  loss$rcp45$none$st$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_stat')
  investment$rcp45$none$st$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_stat')
  preturn$rcp45$none$st$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat')

  cost$rcp85$none$st$no <- ncvar_get(ncdata, 'ExpectedCost_RCP85_stat')
  loss$rcp85$none$st$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_stat')
  investment$rcp85$none$st$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_stat')
  preturn$rcp85$none$st$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat')

  cost$rcp26$none$ns$no <- ncvar_get(ncdata, 'ExpectedCost_RCP26_nonstat')
  loss$rcp26$none$ns$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_nonstat')
  investment$rcp26$none$ns$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_nonstat')
  preturn$rcp26$none$ns$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat')

  cost$rcp45$none$ns$no <- ncvar_get(ncdata, 'ExpectedCost_RCP45_nonstat')
  loss$rcp45$none$ns$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_nonstat')
  investment$rcp45$none$ns$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_nonstat')
  preturn$rcp45$none$ns$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat')

  cost$rcp85$none$ns$no <- ncvar_get(ncdata, 'ExpectedCost_RCP85_nonstat')
  loss$rcp85$none$ns$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_nonstat')
  investment$rcp85$none$ns$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_nonstat')
  preturn$rcp85$none$ns$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat')

  # account for some value to life
  cost$rcp26$none$st$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP26_stat_vosl')
  loss$rcp26$none$st$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_stat_vosl')
  investment$rcp26$none$st$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_stat_vosl')
  preturn$rcp26$none$st$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat_vosl')

  cost$rcp45$none$st$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP45_stat_vosl')
  loss$rcp45$none$st$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_stat_vosl')
  investment$rcp45$none$st$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_stat_vosl')
  preturn$rcp45$none$st$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat_vosl')

  cost$rcp85$none$st$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP85_stat_vosl')
  loss$rcp85$none$st$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_stat_vosl')
  investment$rcp85$none$st$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_stat_vosl')
  preturn$rcp85$none$st$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat_vosl')

  cost$rcp26$none$ns$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP26_nonstat_vosl')
  loss$rcp26$none$ns$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_nonstat_vosl')
  investment$rcp26$none$ns$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_nonstat_vosl')
  preturn$rcp26$none$ns$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat_vosl')

  cost$rcp45$none$ns$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP45_nonstat_vosl')
  loss$rcp45$none$ns$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_nonstat_vosl')
  investment$rcp45$none$ns$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_nonstat_vosl')
  preturn$rcp45$none$ns$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat_vosl')

  cost$rcp85$none$ns$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP85_nonstat_vosl')
  loss$rcp85$none$ns$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_nonstat_vosl')
  investment$rcp85$none$ns$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_nonstat_vosl')
  preturn$rcp85$none$ns$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat_vosl')
nc_close(ncdata)

# read gamma fast dynamics prior results
ncdata <- nc_open(filename.vandantzig.gamma)
  surge.factor$gamma <- ncvar_get(ncdata, 'surge_factor_RCP26') # all RCPs have same surge factor

  # no value of life
  cost$rcp26$gamma$st$no <- ncvar_get(ncdata, 'ExpectedCost_RCP26_stat')
  loss$rcp26$gamma$st$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_stat')
  investment$rcp26$gamma$st$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_stat')
  preturn$rcp26$gamma$st$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat')

  cost$rcp45$gamma$st$no <- ncvar_get(ncdata, 'ExpectedCost_RCP45_stat')
  loss$rcp45$gamma$st$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_stat')
  investment$rcp45$gamma$st$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_stat')
  preturn$rcp45$gamma$st$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat')

  cost$rcp85$gamma$st$no <- ncvar_get(ncdata, 'ExpectedCost_RCP85_stat')
  loss$rcp85$gamma$st$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_stat')
  investment$rcp85$gamma$st$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_stat')
  preturn$rcp85$gamma$st$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat')

  cost$rcp26$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedCost_RCP26_nonstat')
  loss$rcp26$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_nonstat')
  investment$rcp26$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_nonstat')
  preturn$rcp26$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat')

  cost$rcp45$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedCost_RCP45_nonstat')
  loss$rcp45$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_nonstat')
  investment$rcp45$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_nonstat')
  preturn$rcp45$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat')

  cost$rcp85$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedCost_RCP85_nonstat')
  loss$rcp85$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_nonstat')
  investment$rcp85$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_nonstat')
  preturn$rcp85$gamma$ns$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat')

  # account for some value to life
  cost$rcp26$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP26_stat_vosl')
  loss$rcp26$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_stat_vosl')
  investment$rcp26$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_stat_vosl')
  preturn$rcp26$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat_vosl')

  cost$rcp45$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP45_stat_vosl')
  loss$rcp45$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_stat_vosl')
  investment$rcp45$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_stat_vosl')
  preturn$rcp45$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat_vosl')

  cost$rcp85$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP85_stat_vosl')
  loss$rcp85$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_stat_vosl')
  investment$rcp85$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_stat_vosl')
  preturn$rcp85$gamma$st$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat_vosl')

  cost$rcp26$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP26_nonstat_vosl')
  loss$rcp26$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_nonstat_vosl')
  investment$rcp26$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_nonstat_vosl')
  preturn$rcp26$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat_vosl')

  cost$rcp45$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP45_nonstat_vosl')
  loss$rcp45$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_nonstat_vosl')
  investment$rcp45$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_nonstat_vosl')
  preturn$rcp45$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat_vosl')

  cost$rcp85$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP85_nonstat_vosl')
  loss$rcp85$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_nonstat_vosl')
  investment$rcp85$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_nonstat_vosl')
  preturn$rcp85$gamma$ns$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat_vosl')
nc_close(ncdata)

# read uniform fast dynamics prior results
ncdata <- nc_open(filename.vandantzig.uniform)
  surge.factor$uniform <- ncvar_get(ncdata, 'surge_factor_RCP26') # all RCPs have same surge factor

  # no value of life
  cost$rcp26$uniform$st$no <- ncvar_get(ncdata, 'ExpectedCost_RCP26_stat')
  loss$rcp26$uniform$st$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_stat')
  investment$rcp26$uniform$st$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_stat')
  preturn$rcp26$uniform$st$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat')

  cost$rcp45$uniform$st$no <- ncvar_get(ncdata, 'ExpectedCost_RCP45_stat')
  loss$rcp45$uniform$st$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_stat')
  investment$rcp45$uniform$st$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_stat')
  preturn$rcp45$uniform$st$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat')

  cost$rcp85$uniform$st$no <- ncvar_get(ncdata, 'ExpectedCost_RCP85_stat')
  loss$rcp85$uniform$st$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_stat')
  investment$rcp85$uniform$st$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_stat')
  preturn$rcp85$uniform$st$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat')

  cost$rcp26$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedCost_RCP26_nonstat')
  loss$rcp26$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_nonstat')
  investment$rcp26$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_nonstat')
  preturn$rcp26$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat')

  cost$rcp45$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedCost_RCP45_nonstat')
  loss$rcp45$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_nonstat')
  investment$rcp45$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_nonstat')
  preturn$rcp45$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat')

  cost$rcp85$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedCost_RCP85_nonstat')
  loss$rcp85$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_nonstat')
  investment$rcp85$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_nonstat')
  preturn$rcp85$uniform$ns$no <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat')

  # account for some value to life
  cost$rcp26$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP26_stat_vosl')
  loss$rcp26$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_stat_vosl')
  investment$rcp26$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_stat_vosl')
  preturn$rcp26$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_stat_vosl')

  cost$rcp45$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP45_stat_vosl')
  loss$rcp45$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_stat_vosl')
  investment$rcp45$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_stat_vosl')
  preturn$rcp45$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_stat_vosl')

  cost$rcp85$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP85_stat_vosl')
  loss$rcp85$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_stat_vosl')
  investment$rcp85$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_stat_vosl')
  preturn$rcp85$uniform$st$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_stat_vosl')

  cost$rcp26$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP26_nonstat_vosl')
  loss$rcp26$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP26_nonstat_vosl')
  investment$rcp26$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP26_nonstat_vosl')
  preturn$rcp26$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP26_nonstat_vosl')

  cost$rcp45$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP45_nonstat_vosl')
  loss$rcp45$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP45_nonstat_vosl')
  investment$rcp45$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP45_nonstat_vosl')
  preturn$rcp45$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP45_nonstat_vosl')

  cost$rcp85$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedCost_RCP85_nonstat_vosl')
  loss$rcp85$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedLoss_RCP85_nonstat_vosl')
  investment$rcp85$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedInvestment_RCP85_nonstat_vosl')
  preturn$rcp85$uniform$ns$yes <- ncvar_get(ncdata, 'ExpectedPreturn_RCP85_nonstat_vosl')
nc_close(ncdata)

n.ensemble <- c(ncol(cost$rcp85$none$ns$no) , ncol(cost$rcp85$gamma$ns$no) , ncol(cost$rcp85$uniform$ns$no) )
names(n.ensemble) <- scen.ais
n.height <- length(heightening)

## Initialize arrays for average cost and return period among states of world
## within a scenario
Cavg <- rep(NA,n.height)
Ravg <- rep(NA,n.height)

## What is the economically efficient heightening for each ensemble member?
## (and the index, so we can grab the return period)
ieum <- vector("list",3); names(ieum)<-scen.rcp     # EUM-optimal heightening index for scenario
Ceum <- vector("list",3); names(Ceum)<-scen.rcp     # EUM-optimal cost for scenario
Reum <- vector("list",3); names(Reum)<-scen.rcp     # EUM-optimal return period scenario
Heum <- vector("list",3); names(Heum)<-scen.rcp     # EUM-optimal heightening for scenario

Cens.eum <- vector("list",3); names(Cens)<-scen.rcp     # each SOW expected cost, at EUM-optimal heightening
Rens.eum <- vector("list",3); names(Rens)<-scen.rcp     # each SOW return period, at EUM-optimal heightening

Crob <- vector("list",3); names(Crob)<-scen.rcp     # robust cost (CVaR)
irob <- vector("list",3); names(irob)<-scen.rcp     # robust index
Rrob <- vector("list",3); names(Rrob)<-scen.rcp     # robust return period
Hrob <- vector("list",3); names(Hrob)<-scen.rcp     # robust heightening

for (rcp in scen.rcp) {
  ieum[[rcp]] <- vector("list",3); names(ieum[[rcp]])<-scen.ais
  Ceum[[rcp]] <- vector("list",3); names(Ceum[[rcp]])<-scen.ais
  Reum[[rcp]] <- vector("list",3); names(Reum[[rcp]])<-scen.ais
  Heum[[rcp]] <- vector("list",3); names(Heum[[rcp]])<-scen.ais
  Cens[[rcp]] <- vector("list",3); names(Cens[[rcp]])<-scen.ais
  Rens[[rcp]] <- vector("list",3); names(Rens[[rcp]])<-scen.ais
  irob[[rcp]] <- vector("list",3); names(irob[[rcp]])<-scen.ais
  Crob[[rcp]] <- vector("list",3); names(Crob[[rcp]])<-scen.ais
  Rrob[[rcp]] <- vector("list",3); names(Rrob[[rcp]])<-scen.ais
  Hrob[[rcp]] <- vector("list",3); names(Hrob[[rcp]])<-scen.ais
  for (ais in scen.ais) {
    ieum[[rcp]][[ais]] <- vector("list",2); names(ieum[[rcp]][[ais]])<-scen.ss
    Ceum[[rcp]][[ais]] <- vector("list",2); names(Ceum[[rcp]][[ais]])<-scen.ss
    Reum[[rcp]][[ais]] <- vector("list",2); names(Reum[[rcp]][[ais]])<-scen.ss
    Heum[[rcp]][[ais]] <- vector("list",2); names(Heum[[rcp]][[ais]])<-scen.ss
    Cens[[rcp]][[ais]] <- vector("list",2); names(Cens[[rcp]][[ais]])<-scen.ss
    Rens[[rcp]][[ais]] <- vector("list",2); names(Rens[[rcp]][[ais]])<-scen.ss
    irob[[rcp]][[ais]] <- vector("list",2); names(irob[[rcp]][[ais]])<-scen.ss
    Crob[[rcp]][[ais]] <- vector("list",2); names(Crob[[rcp]][[ais]])<-scen.ss
    Rrob[[rcp]][[ais]] <- vector("list",2); names(Rrob[[rcp]][[ais]])<-scen.ss
    Hrob[[rcp]][[ais]] <- vector("list",2); names(Hrob[[rcp]][[ais]])<-scen.ss
    for (ss in scen.ss) {
      ieum[[rcp]][[ais]][[ss]] <- vector("list",2); names(ieum[[rcp]][[ais]][[ss]])<-scen.vosl
      Ceum[[rcp]][[ais]][[ss]] <- vector("list",2); names(Ceum[[rcp]][[ais]][[ss]])<-scen.vosl
      Reum[[rcp]][[ais]][[ss]] <- vector("list",2); names(Reum[[rcp]][[ais]][[ss]])<-scen.vosl
      Heum[[rcp]][[ais]][[ss]] <- vector("list",2); names(Heum[[rcp]][[ais]][[ss]])<-scen.vosl
      Cens[[rcp]][[ais]][[ss]] <- vector("list",2); names(Cens[[rcp]][[ais]][[ss]])<-scen.vosl
      Rens[[rcp]][[ais]][[ss]] <- vector("list",2); names(Rens[[rcp]][[ais]][[ss]])<-scen.vosl
      irob[[rcp]][[ais]][[ss]] <- vector("list",2); names(irob[[rcp]][[ais]][[ss]])<-scen.vosl
      Crob[[rcp]][[ais]][[ss]] <- vector("list",2); names(Crob[[rcp]][[ais]][[ss]])<-scen.vosl
      Rrob[[rcp]][[ais]][[ss]] <- vector("list",2); names(Rrob[[rcp]][[ais]][[ss]])<-scen.vosl
      Hrob[[rcp]][[ais]][[ss]] <- vector("list",2); names(Hrob[[rcp]][[ais]][[ss]])<-scen.vosl
      for (vosl in scen.vosl) {
        ieum[[rcp]][[ais]][[ss]][[vosl]] <- NA
        Ceum[[rcp]][[ais]][[ss]][[vosl]] <- NA
        Reum[[rcp]][[ais]][[ss]][[vosl]] <- NA
        Heum[[rcp]][[ais]][[ss]][[vosl]] <- NA
        irob[[rcp]][[ais]][[ss]][[vosl]] <- NA
        Crob[[rcp]][[ais]][[ss]][[vosl]] <- NA
        Rrob[[rcp]][[ais]][[ss]][[vosl]] <- NA
        Hrob[[rcp]][[ais]][[ss]][[vosl]] <- NA
      }
    }
  }
}

# calculate the optimal strategy (EUM, min-Regret, ...?) for each ensemble
# member, in each scenario

for (ais in scen.ais) {
  for (rcp in scen.rcp) {
    for (ss in scen.ss) {
      for (vosl in scen.vosl) {
        for (dh in 1:n.height) {
          Cavg[dh] <- mean( cost[[rcp]][[ais]][[ss]][[vosl]][dh,] )
          Ravg[dh] <- mean( preturn[[rcp]][[ais]][[ss]][[vosl]][dh,] )
        }
        Ceum[[rcp]][[ais]][[ss]][[vosl]] <- min(Cavg)
        itmp <- which(Cavg==min(Cavg))
        if(length(itmp)>1) {itmp <- median(itmp)}
        ieum[[rcp]][[ais]][[ss]][[vosl]] <- itmp
        Heum[[rcp]][[ais]][[ss]][[vosl]] <- heightening[ieum[[rcp]][[ais]][[ss]][[vosl]]]
        Reum[[rcp]][[ais]][[ss]][[vosl]] <- Ravg[itmp]
        Cens[[rcp]][[ais]][[ss]][[vosl]] <- cost[[rcp]][[ais]][[ss]][[vosl]][itmp,]
        Rens[[rcp]][[ais]][[ss]][[vosl]] <- preturn[[rcp]][[ais]][[ss]][[vosl]][itmp,]
      }
    }
  }
}

##==============================================================================
##==============================================================================


##==============================================================================
## eumimal in the mean sense
##==============================================================================

R.rcp26.no <- data.frame(rbind( c(Reum$rcp26$none$st$no, Reum$rcp26$gamma$st$no, Reum$rcp26$uniform$st$no),
                                c(Reum$rcp26$none$ns$no, Reum$rcp26$gamma$ns$no, Reum$rcp26$uniform$ns$no)))
R.rcp45.no <- data.frame(rbind( c(Reum$rcp45$none$st$no, Reum$rcp45$gamma$st$no, Reum$rcp45$uniform$st$no),
                                c(Reum$rcp45$none$ns$no, Reum$rcp45$gamma$ns$no, Reum$rcp45$uniform$ns$no)))
R.rcp85.no <- data.frame(rbind( c(Reum$rcp85$none$st$no, Reum$rcp85$gamma$st$no, Reum$rcp85$uniform$st$no),
                                c(Reum$rcp85$none$ns$no, Reum$rcp85$gamma$ns$no, Reum$rcp85$uniform$ns$no)))

H.rcp26.no <- data.frame(rbind( c(Heum$rcp26$none$st$no, Heum$rcp26$gamma$st$no, Heum$rcp26$uniform$st$no),
                                c(Heum$rcp26$none$ns$no, Heum$rcp26$gamma$ns$no, Heum$rcp26$uniform$ns$no)))
H.rcp45.no <- data.frame(rbind( c(Heum$rcp45$none$st$no, Heum$rcp45$gamma$st$no, Heum$rcp45$uniform$st$no),
                                c(Heum$rcp45$none$ns$no, Heum$rcp45$gamma$ns$no, Heum$rcp45$uniform$ns$no)))
H.rcp85.no <- data.frame(rbind( c(Heum$rcp85$none$st$no, Heum$rcp85$gamma$st$no, Heum$rcp85$uniform$st$no),
                                c(Heum$rcp85$none$ns$no, Heum$rcp85$gamma$ns$no, Heum$rcp85$uniform$ns$no)))

C.rcp26.no <- data.frame(rbind( c(Ceum$rcp26$none$st$no, Ceum$rcp26$gamma$st$no, Ceum$rcp26$uniform$st$no),
                                c(Ceum$rcp26$none$ns$no, Ceum$rcp26$gamma$ns$no, Ceum$rcp26$uniform$ns$no)))/1e9
C.rcp45.no <- data.frame(rbind( c(Ceum$rcp45$none$st$no, Ceum$rcp45$gamma$st$no, Ceum$rcp45$uniform$st$no),
                                c(Ceum$rcp45$none$ns$no, Ceum$rcp45$gamma$ns$no, Ceum$rcp45$uniform$ns$no)))/1e9
C.rcp85.no <- data.frame(rbind( c(Ceum$rcp85$none$st$no, Ceum$rcp85$gamma$st$no, Ceum$rcp85$uniform$st$no),
                                c(Ceum$rcp85$none$ns$no, Ceum$rcp85$gamma$ns$no, Ceum$rcp85$uniform$ns$no)))/1e9

# means, with VOSL
R.rcp26.yes <- data.frame(rbind( c(Reum$rcp26$none$st$yes, Reum$rcp26$gamma$st$yes, Reum$rcp26$uniform$st$yes),
                                 c(Reum$rcp26$none$ns$yes, Reum$rcp26$gamma$ns$yes, Reum$rcp26$uniform$ns$yes)))
R.rcp45.yes <- data.frame(rbind( c(Reum$rcp45$none$st$yes, Reum$rcp45$gamma$st$yes, Reum$rcp45$uniform$st$yes),
                                 c(Reum$rcp45$none$ns$yes, Reum$rcp45$gamma$ns$yes, Reum$rcp45$uniform$ns$yes)))
R.rcp85.yes <- data.frame(rbind( c(Reum$rcp85$none$st$yes, Reum$rcp85$gamma$st$yes, Reum$rcp85$uniform$st$yes),
                                 c(Reum$rcp85$none$ns$yes, Reum$rcp85$gamma$ns$yes, Reum$rcp85$uniform$ns$yes)))

H.rcp26.yes <- data.frame(rbind( c(Heum$rcp26$none$st$yes, Heum$rcp26$gamma$st$yes, Heum$rcp26$uniform$st$yes),
                                 c(Heum$rcp26$none$ns$yes, Heum$rcp26$gamma$ns$yes, Heum$rcp26$uniform$ns$yes)))
H.rcp45.yes <- data.frame(rbind( c(Heum$rcp45$none$st$yes, Heum$rcp45$gamma$st$yes, Heum$rcp45$uniform$st$yes),
                                 c(Heum$rcp45$none$ns$yes, Heum$rcp45$gamma$ns$yes, Heum$rcp45$uniform$ns$yes)))
H.rcp85.yes <- data.frame(rbind( c(Heum$rcp85$none$st$yes, Heum$rcp85$gamma$st$yes, Heum$rcp85$uniform$st$yes),
                                 c(Heum$rcp85$none$ns$yes, Heum$rcp85$gamma$ns$yes, Heum$rcp85$uniform$ns$yes)))

C.rcp26.yes <- data.frame(rbind( c(Ceum$rcp26$none$st$yes, Ceum$rcp26$gamma$st$yes, Ceum$rcp26$uniform$st$yes),
                                 c(Ceum$rcp26$none$ns$yes, Ceum$rcp26$gamma$ns$yes, Ceum$rcp26$uniform$ns$yes)))/1e9
C.rcp45.yes <- data.frame(rbind( c(Ceum$rcp45$none$st$yes, Ceum$rcp45$gamma$st$yes, Ceum$rcp45$uniform$st$yes),
                                 c(Ceum$rcp45$none$ns$yes, Ceum$rcp45$gamma$ns$yes, Ceum$rcp45$uniform$ns$yes)))/1e9
C.rcp85.yes <- data.frame(rbind( c(Ceum$rcp85$none$st$yes, Ceum$rcp85$gamma$st$yes, Ceum$rcp85$uniform$st$yes),
                                 c(Ceum$rcp85$none$ns$yes, Ceum$rcp85$gamma$ns$yes, Ceum$rcp85$uniform$ns$yes)))/1e9

C.min.no <- min(c(as.vector(as.matrix(C.rcp26.no)) ,as.vector(as.matrix(C.rcp45.no)) ,as.vector(as.matrix(C.rcp85.no)) ))
C.min.yes<- min(c(as.vector(as.matrix(C.rcp26.yes)),as.vector(as.matrix(C.rcp45.yes)),as.vector(as.matrix(C.rcp85.yes))))
C.min.no <- min(C.min.no, C.min.yes)
C.min.yes<- min(C.min.no, C.min.yes)

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

flood.lower <- 1/1000
flood.upper <- 1/100
flood.n <- 2^10

ncdata <- nc_open(filename.brick.nofd)
  mod.time <- ncvar_get(ncdata, 'time_proj')
  lsl$rcp26$none <- ncvar_get(ncdata, 'LocalSeaLevel_RCP26')
  lsl$rcp45$none <- ncvar_get(ncdata, 'LocalSeaLevel_RCP45')
  lsl$rcp85$none <- ncvar_get(ncdata, 'LocalSeaLevel_RCP85')
nc_close(ncdata)

ncdata <- nc_open(filename.brick.gamma)
  lsl$rcp26$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_RCP26')
  lsl$rcp45$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_RCP45')
  lsl$rcp85$gamma <- ncvar_get(ncdata, 'LocalSeaLevel_RCP85')
nc_close(ncdata)

ncdata <- nc_open(filename.brick.uniform)
  lsl$rcp26$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_RCP26')
  lsl$rcp45$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_RCP45')
  lsl$rcp85$uniform <- ncvar_get(ncdata, 'LocalSeaLevel_RCP85')
nc_close(ncdata)

## initialize
for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    f.surlev[[rcp]][[ais]] <- vector("list",2); names(f.surlev[[rcp]][[ais]])=scen.ss
    f.seasurlev[[rcp]][[ais]] <- vector("list",2); names(f.seasurlev[[rcp]][[ais]])=scen.ss
    f.flood[[rcp]][[ais]] <- vector("list",2); names(f.flood[[rcp]][[ais]])=scen.ss
  }
}

## distribution of sea-level
iproj <- which(mod.time==2065)
inorm <- which(mod.time==2015)
for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    lsl[[rcp]][[ais]] <- lsl[[rcp]][[ais]][iproj,] - lsl[[rcp]][[ais]][inorm,]
    f.sealev[[rcp]][[ais]] <- density(x=lsl[[rcp]][[ais]], from=lsl.lower, to=lsl.upper, n=lsl.n)
    f.surge.rise[[rcp]][[ais]] <- density(x=lsl[[rcp]][[ais]]*surge.factor[[ais]], from=lsl.lower, to=lsl.upper, n=lsl.n)
    # and initialize storm surge
    for (ss in scen.ss) {
      f.surlev[[rcp]][[ais]][[ss]] <- rep(NA, length(f.sealev[[rcp]][[ais]]$y))
    }
  }
}

x.lsl <- f.sealev$rcp85$none$x
dx.lsl <- mean(diff(x.lsl))

## distribution of storm surge level
## --> convolution of surge.rise distribution + gev.stat distribution

## stationary case (x*1000, because GEV parameters assume levels are in mm, not m)
f.surge.stat <- as.numeric(dgev(x=x.lsl*1000, xi=gev.stat[1,'shape'], mu=gev.stat[1,'location'], beta=gev.stat[1,'scale'], log=FALSE))
f.surge.stat <- f.surge.stat/sum( dx.lsl*f.surge.stat)

iproj <- which(mod.time==2065)
inorm <- which(mod.time==2015)
for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    for (ss in scen.ss) {
      if(ss=='st') {
        f.surlev[[rcp]][[ais]][[ss]] <- f.surge.stat
        f.seasurlev[[rcp]][[ais]][[ss]] <- convolve(x=f.surlev[[rcp]][[ais]][[ss]], y=rev(f.sealev[[rcp]][[ais]]$y))
        f.seasurlev[[rcp]][[ais]][[ss]] <- f.seasurlev[[rcp]][[ais]][[ss]]/sum(dx.lsl*f.seasurlev[[rcp]][[ais]][[ss]])
      } else {
        f.surlev[[rcp]][[ais]][[ss]] <- convolve(x=f.surge.rise[[rcp]][[ais]]$y, y=rev(f.surge.stat))
        f.surlev[[rcp]][[ais]][[ss]] <- f.surlev[[rcp]][[ais]][[ss]]/sum(dx.lsl*f.surlev[[rcp]][[ais]][[ss]])
        f.seasurlev[[rcp]][[ais]][[ss]] <- convolve(x=f.surlev[[rcp]][[ais]][[ss]], y=rev(f.sealev[[rcp]][[ais]]$y))
        f.seasurlev[[rcp]][[ais]][[ss]] <- f.seasurlev[[rcp]][[ais]][[ss]]/sum(dx.lsl*f.seasurlev[[rcp]][[ais]][[ss]])
      }
    }
  }
}

## distribution of flood risk
## (each SOW within scenario, using the EUM-optimal heightening for that scenario)
for (rcp in scen.rcp) {
  for (ais in scen.ais) {
    for (ss in scen.ss) {
      f.flood[[rcp]][[ais]][[ss]] <- density(x=1/Rens[[rcp]][[ais]][[ss]]$no, from=flood.lower, to=flood.upper, n=flood.n)
    }
  }
}

x.flood <- f.flood$rcp85$none$st$x
dx.flood <- mean(diff(x.flood))



## THE ACTUAL FIGURE

pdf(paste(plotdir,'pdfs_slr-ss-risk.pdf',sep=''),width=3.5,height=8,colormodel='cmyk')

par(mfrow=c(3,1), mai=c(.6,.35,.2,.5))
# (a) pdfs of local sea-level rise (2065)
plot(x.lsl, f.sealev$rcp26$none$st$y, type='l', xlim=c(0,1.1), ylim=c(0,10), lwd=1.5,
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
par(mai=c(.65,.35,.15,.5))
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
text(1.01,2.4,"non-stationary", pos=4, cex=1.3)

#legend(0.7,4,c("stationary","all others: non-stationary\n (see panel (a) legend)"),
#       lty=c(1,NA), lwd=2, cex=1.2,
#       col=c(rgb(mycol[6,1],mycol[6,2],mycol[6,3]), NA),
#       bty='n')

# (c) pdfs of flood risk (average annual exceedance probability)
par(mai=c(.5,.35,.5,.5))
plot(x.flood, f.flood$rcp26$none$st$y, type='l', xlim=c(0.0034,0.0046), ylim=c(0,6000), lwd=1.5,
     col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=1,
     xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i', axes=FALSE)
lines(x.flood, f.flood$rcp26$gamma$st$y, type='l', col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=2)
lines(x.flood, f.flood$rcp26$uniform$st$y, type='l', col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=3)
lines(x.flood, f.flood$rcp45$none$st$y, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=1)
lines(x.flood, f.flood$rcp45$gamma$st$y, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=2)
lines(x.flood, f.flood$rcp45$uniform$st$y, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=3)
lines(x.flood, f.flood$rcp85$none$st$y, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=1)
lines(x.flood, f.flood$rcp85$gamma$st$y, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=2)
lines(x.flood, f.flood$rcp85$uniform$st$y, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=3)

lines(x.flood, f.flood$rcp26$none$ns$y, type='l', col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=2)
lines(x.flood, f.flood$rcp26$gamma$ns$y, type='l', col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=2)
lines(x.flood, f.flood$rcp26$uniform$ns$y, type='l', col=rgb(mycol[13,1],mycol[13,2],mycol[13,3]), lty=3)
lines(x.flood, f.flood$rcp45$none$ns$y, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=1)
lines(x.flood, f.flood$rcp45$gamma$ns$y, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=2)
lines(x.flood, f.flood$rcp45$uniform$ns$y, type='l', col=rgb(mycol[2,1],mycol[2,2],mycol[2,3]), lty=3)
lines(x.flood, f.flood$rcp85$none$ns$y, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=1)
lines(x.flood, f.flood$rcp85$gamma$ns$y, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=2)
lines(x.flood, f.flood$rcp85$uniform$ns$y, type='l', col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lty=3)

xlabels <- seq(0.0034,0.0046,0.0002)
axis(1,xlabels, cex.axis=1.3)
axis(3,at=xlabels,labels=round(1/xlabels), cex.axis=1.3)
u <- par("usr")
mtext('Probability density', side=2, line=1.2, cex=.9);
mtext('Projected annual flood probability', side=1, line=2.3, cex=.9);
mtext('Projected return period [years]', side=3, line=2.3, cex=.9);
mtext(side=3, text=expression(bold('   c')), line=-1.4, cex=.9, adj=0);
arrows(u[1], u[3], u[1], .94*u[4], code = 2, xpd = TRUE)
text(.00383,5300,"stationary", pos=4, cex=1.3)
text(.00390,3800,"non-stationary", pos=4, cex=1.3)


dev.off()




##==============================================================================
##==============================================================================





##==============================================================================
## CVaR(alpha-tail)
##==============================================================================

# calculate CVaR for optimal cost distributions (Crob), and for heightening and
# return period (Hrob, Rrob)
alpha <- 0.05

for (ais in scen.ais) {
  for (rcp in scen.rcp) {
    for (ss in scen.ss) {
      for (vosl in scen.vosl) {
        q.cvar <- rev(sort(Copt[[rcp]][[ais]][[ss]][[vosl]]))[floor(alpha*length(Copt[[rcp]][[ais]][[ss]][[vosl]]))]
        i.cvar <- which(Copt[[rcp]][[ais]][[ss]][[vosl]] >= q.cvar)
        Crob[[rcp]][[ais]][[ss]][[vosl]] <- mean( Copt[[rcp]][[ais]][[ss]][[vosl]][i.cvar] )
        Hrob[[rcp]][[ais]][[ss]][[vosl]] <- mean( Hopt[[rcp]][[ais]][[ss]][[vosl]][i.cvar] )
        Rrob[[rcp]][[ais]][[ss]][[vosl]] <- mean( Ropt[[rcp]][[ais]][[ss]][[vosl]][i.cvar] )
      }
    }
  }
}

# no VOSL
R.rcp26.no <- data.frame(rbind( c(Rrob$rcp26$none$st$no, Rrob$rcp26$gamma$st$no, Rrob$rcp26$uniform$st$no),
                                c(Rrob$rcp26$none$ns$no, Rrob$rcp26$gamma$ns$no, Rrob$rcp26$uniform$ns$no)))
R.rcp45.no <- data.frame(rbind( c(Rrob$rcp45$none$st$no, Rrob$rcp45$gamma$st$no, Rrob$rcp45$uniform$st$no),
                                c(Rrob$rcp45$none$ns$no, Rrob$rcp45$gamma$ns$no, Rrob$rcp45$uniform$ns$no)))
R.rcp85.no <- data.frame(rbind( c(Rrob$rcp85$none$st$no, Rrob$rcp85$gamma$st$no, Rrob$rcp85$uniform$st$no),
                                c(Rrob$rcp85$none$ns$no, Rrob$rcp85$gamma$ns$no, Rrob$rcp85$uniform$ns$no)))

H.rcp26.no <- data.frame(rbind( c(Hrob$rcp26$none$st$no, Hrob$rcp26$gamma$st$no, Hrob$rcp26$uniform$st$no),
                                c(Hrob$rcp26$none$ns$no, Hrob$rcp26$gamma$ns$no, Hrob$rcp26$uniform$ns$no)))
H.rcp45.no <- data.frame(rbind( c(Hrob$rcp45$none$st$no, Hrob$rcp45$gamma$st$no, Hrob$rcp45$uniform$st$no),
                                c(Hrob$rcp45$none$ns$no, Hrob$rcp45$gamma$ns$no, Hrob$rcp45$uniform$ns$no)))
H.rcp85.no <- data.frame(rbind( c(Hrob$rcp85$none$st$no, Hrob$rcp85$gamma$st$no, Hrob$rcp85$uniform$st$no),
                                c(Hrob$rcp85$none$ns$no, Hrob$rcp85$gamma$ns$no, Hrob$rcp85$uniform$ns$no)))

C.rcp26.no <- data.frame(rbind( c(Crob$rcp26$none$st$no, Crob$rcp26$gamma$st$no, Crob$rcp26$uniform$st$no),
                                c(Crob$rcp26$none$ns$no, Crob$rcp26$gamma$ns$no, Crob$rcp26$uniform$ns$no)))/1e9
C.rcp45.no <- data.frame(rbind( c(Crob$rcp45$none$st$no, Crob$rcp45$gamma$st$no, Crob$rcp45$uniform$st$no),
                                c(Crob$rcp45$none$ns$no, Crob$rcp45$gamma$ns$no, Crob$rcp45$uniform$ns$no)))/1e9
C.rcp85.no <- data.frame(rbind( c(Crob$rcp85$none$st$no, Crob$rcp85$gamma$st$no, Crob$rcp85$uniform$st$no),
                                c(Crob$rcp85$none$ns$no, Crob$rcp85$gamma$ns$no, Crob$rcp85$uniform$ns$no)))/1e9

# with VOSL
R.rcp26.yes <- data.frame(rbind( c(Rrob$rcp26$none$st$yes, Rrob$rcp26$gamma$st$yes, Rrob$rcp26$uniform$st$yes),
                                 c(Rrob$rcp26$none$ns$yes, Rrob$rcp26$gamma$ns$yes, Rrob$rcp26$uniform$ns$yes)))
R.rcp45.yes <- data.frame(rbind( c(Rrob$rcp45$none$st$yes, Rrob$rcp45$gamma$st$yes, Rrob$rcp45$uniform$st$yes),
                                 c(Rrob$rcp45$none$ns$yes, Rrob$rcp45$gamma$ns$yes, Rrob$rcp45$uniform$ns$yes)))
R.rcp85.yes <- data.frame(rbind( c(Rrob$rcp85$none$st$yes, Rrob$rcp85$gamma$st$yes, Rrob$rcp85$uniform$st$yes),
                                 c(Rrob$rcp85$none$ns$yes, Rrob$rcp85$gamma$ns$yes, Rrob$rcp85$uniform$ns$yes)))

H.rcp26.yes <- data.frame(rbind( c(Hrob$rcp26$none$st$yes, Hrob$rcp26$gamma$st$yes, Hrob$rcp26$uniform$st$yes),
                                 c(Hrob$rcp26$none$ns$yes, Hrob$rcp26$gamma$ns$yes, Hrob$rcp26$uniform$ns$yes)))
H.rcp45.yes <- data.frame(rbind( c(Hrob$rcp45$none$st$yes, Hrob$rcp45$gamma$st$yes, Hrob$rcp45$uniform$st$yes),
                                 c(Hrob$rcp45$none$ns$yes, Hrob$rcp45$gamma$ns$yes, Hrob$rcp45$uniform$ns$yes)))
H.rcp85.yes <- data.frame(rbind( c(Hrob$rcp85$none$st$yes, Hrob$rcp85$gamma$st$yes, Hrob$rcp85$uniform$st$yes),
                                 c(Hrob$rcp85$none$ns$yes, Hrob$rcp85$gamma$ns$yes, Hrob$rcp85$uniform$ns$yes)))

C.rcp26.yes <- data.frame(rbind( c(Crob$rcp26$none$st$yes, Crob$rcp26$gamma$st$yes, Crob$rcp26$uniform$st$yes),
                                 c(Crob$rcp26$none$ns$yes, Crob$rcp26$gamma$ns$yes, Crob$rcp26$uniform$ns$yes)))/1e9
C.rcp45.yes <- data.frame(rbind( c(Crob$rcp45$none$st$yes, Crob$rcp45$gamma$st$yes, Crob$rcp45$uniform$st$yes),
                                 c(Crob$rcp45$none$ns$yes, Crob$rcp45$gamma$ns$yes, Crob$rcp45$uniform$ns$yes)))/1e9
C.rcp85.yes <- data.frame(rbind( c(Crob$rcp85$none$st$yes, Crob$rcp85$gamma$st$yes, Crob$rcp85$uniform$st$yes),
                                 c(Crob$rcp85$none$ns$yes, Crob$rcp85$gamma$ns$yes, Crob$rcp85$uniform$ns$yes)))/1e9

C.min.no <- min(c(as.vector(as.matrix(C.rcp26.no)) ,as.vector(as.matrix(C.rcp45.no)) ,as.vector(as.matrix(C.rcp85.no)) ))
C.min.yes<- min(c(as.vector(as.matrix(C.rcp26.yes)),as.vector(as.matrix(C.rcp45.yes)),as.vector(as.matrix(C.rcp85.yes))))
C.min.no <- min(C.min.no, C.min.yes)
C.min.yes<- min(C.min.no, C.min.yes)

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

table.cost <- rbind(c(paste("ensemble CVaR(",alpha,") total cost percent increase from best case",sep=""),"",""),table.rcp26.no, c("","",""), table.rcp45.no, c("","",""), table.rcp85.no,c("","",""),table.rcp26.yes, c("","",""), table.rcp45.yes, c("","",""), table.rcp85.yes)


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

table.returnperiod <- rbind(c(paste("ensemble CVaR(",alpha,") return period",sep=""),"",""),table.rcp26.no, c("","",""), table.rcp45.no, c("","",""), table.rcp85.no,c("","",""),table.rcp26.yes, c("","",""), table.rcp45.yes, c("","",""), table.rcp85.yes)


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

table.heightening <- rbind(c(paste("ensemble CVaR(",alpha,") heightening",sep=""),"",""),table.rcp26.no, c("","",""), table.rcp45.no, c("","",""), table.rcp85.no,c("","",""),table.rcp26.yes, c("","",""), table.rcp45.yes, c("","",""), table.rcp85.yes)

table.cvar <- cbind(table.cost,rep("",18),table.heightening,rep("",18),table.returnperiod)
write.csv(table.cvar, "../output_model/scenarios_cvar.csv")

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
##==============================================================================
## end of main paper stuff. so far, nothing below this line is included AT ALL
##==============================================================================
##==============================================================================



##==============================================================================
##==============================================================================







##==============================================================================
## End
##==============================================================================
