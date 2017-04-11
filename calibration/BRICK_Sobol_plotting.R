##==============================================================================
## Sobol sensitvity analysis for drivers of flood risk.
## --> plotting routines
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


##==============================================================================
## Radial plots

## Code history:
# this one by Tony Wong, 6 March 2017, Penn State. Modified a bit realtive to
# the previous versions listed below.
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
n_params <- 38

# Set Sobol indices file name
Sobol_file_1 <- "../output_calibration/BRICK_Sobol-1-tot_31Mar2017-Build-GEV-2065.txt"
Sobol_file_2 <- "../output_calibration/BRICK_Sobol-2_31Mar2017-Build-GEV-2065.txt"

noGEV <- FALSE
noHR <- TRUE

####################################
# Import data from sensitivity analysis

# First- and total-order indices
s1st <- read.csv(Sobol_file_1,
                  sep=' ',
                  header=TRUE,
                  nrows = n_params,
                  as.is=c(TRUE,rep(FALSE,5)))

parnames.sobol <- s1st[,1]

# Import second-order indices
s2_table <- read.csv(Sobol_file_2,
               sep=' ',
               header=TRUE,
               nrows = n_params*(n_params-1)/2,
               as.is=c(TRUE,rep(FALSE,4)))

# Convert second-order to upper-triangular matrix
s2 <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2)
s2 <- as.data.frame(s2)
colnames(s2) <- rownames(s2) <- s1st$Parameter

# Convert confidence intervals to upper-triangular matrix
s2_conf_low <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2_conf_high <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2_conf_low[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_low)
s2_conf_high[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_high)

s2_conf_low <- as.data.frame(s2_conf_low)
s2_conf_high <- as.data.frame(s2_conf_high)
colnames(s2_conf_low) <- rownames(s2_conf_low) <- s1st$Parameter
colnames(s2_conf_high) <- rownames(s2_conf_high) <- s1st$Parameter

####################################
# Determine which indices are statistically significant

# S1 & ST: using the confidence intervals
s1st1 <- stat_sig_s1st(s1st
                      ,method="congtr"
                      ,greater=0.01
                      ,sigCri='either')

# S1 & ST: using greater than a given value
#s1st1 <- stat_sig_s1st(s1st
#                      ,method="gtr"
#                      ,greater=0.01
#                      ,sigCri='either')

# S2: using the confidence intervals
s2_sig1 <- stat_sig_s2(s2
                       ,s2_conf_low
                       ,s2_conf_high
                       ,method='congtr'
                       ,greater=0.01
                       )

# S2: using greater than a given value
#s2_sig1 <- stat_sig_s2(s2
#                       ,s2_conf
#                       ,greater=0.02
#                       ,method='gtr')

####################################
# Define groups for the variables and the color schemes
# Defining lists of the variables for each group

name_list1 <- list('Temperature' = parnames.sobol[1:5]
                   ,'Sea Level:\n   Glaciers & Ice Caps' = parnames.sobol[6:9]
                   ,'Sea Level:\nThermal Expansion' = parnames.sobol[10:13]
                   ,'Sea Level:\nGreenland Ice Sheet' = parnames.sobol[14:18]
                   ,'Sea Level:\nAntarctic Ice Sheet' = parnames.sobol[19:33]
                   ,'Land\nSubsidence' = parnames.sobol[35]
                   #,'Storm Surge' = parnames.sobol[c(34,36:38)]
                   ,'Storm Surge' = parnames.sobol[c(34)]
                   ,'Emissions' = parnames.sobol[36]
                   ,'Protection' = parnames.sobol[37]
                   )

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
                  expression('C'['surge']), 'subs'
                  #, expression(mu), expression(sigma), expression(xi)
                  , 'RCP'
                  , 'build'
                  )

# modify the names and symbols if we are plotting one of the SOM figures
if(noGEV) {
    name_list1 <- list('Temperature' = parnames.sobol[1:5]
                       ,'Sea Level:\n   Glaciers & Ice Caps' = parnames.sobol[6:9]
                       ,'Sea Level:\nThermal Expansion' = parnames.sobol[10:13]
                       ,'Sea Level:\nGreenland Ice Sheet' = parnames.sobol[14:18]
                       ,'Sea Level:\nAntarctic Ice Sheet' = parnames.sobol[19:33]
                       ,'Land\nSubsidence' = parnames.sobol[35]
                       #,'Storm Surge' = parnames.sobol[c(34,36:38)]
                       ,'Storm Surge' = parnames.sobol[c(34)]
                       ,'Emissions' = parnames.sobol[36]
                       ,'Protection' = parnames.sobol[37]
                       )

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
                      expression('C'['surge']), 'subs'
                      #, expression(mu), expression(sigma), expression(xi)
                      , 'RCP'
                      , 'build'
                      )
} else if(noHR) {
    name_list1 <- list('Temperature' = parnames.sobol[1:5]
                       ,'Sea Level:\n   Glaciers & Ice Caps' = parnames.sobol[6:9]
                       ,'Sea Level:\nThermal Expansion' = parnames.sobol[10:13]
                       ,'Sea Level:\nGreenland Ice Sheet' = parnames.sobol[14:18]
                       ,'Sea Level:\nAntarctic Ice Sheet' = parnames.sobol[19:31]
                       ,'Land\nSubsidence' = parnames.sobol[33]
                       ,'Storm Surge' = parnames.sobol[c(32,34:36)]
                       ,'Emissions' = parnames.sobol[37]
                       ,'Protection' = parnames.sobol[38]
                       )

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
                      expression(kappa['AIS']), expression('f'[0])
                      , expression('b'[0]), 'slope'
                      ,expression(lambda), expression('T'['crit'])
                      ,expression('C'['surge']), 'subs'
                      , expression(mu), expression(sigma), expression(xi)
                      , 'RCP'
                      , 'build'
                      )
}

source('../Useful/colorblindPalette.R')

# defining list of colors for each group
col_list1 <- list("Temperature"     = rgb(mycol[11,1],mycol[11,2],mycol[11,3])
                  ,'Sea Level:\n   Glaciers & Ice Caps' = rgb(mycol[3,1],mycol[3,2],mycol[3,3])
                  ,'Sea Level:\nThermal Expansion'   = rgb(mycol[9,1],mycol[9,2],mycol[9,3])
                  ,'Sea Level:\nGreenland Ice Sheet'  = rgb(mycol[2,1],mycol[2,2],mycol[2,3])
                  ,'Sea Level:\nAntarctic Ice Sheet'  = rgb(mycol[7,1],mycol[7,2],mycol[7,3])
                  ,'Land\nSubsidence' = rgb(mycol[1,1],mycol[1,2],mycol[1,3])
                  ,"Storm Surge"     = rgb(mycol[6,1],mycol[6,2],mycol[6,3])
                  ,"Emissions"       = rgb(mycol[13,1],mycol[13,2],mycol[13,3])
                  ,"Protection"      = rgb(mycol[12,1],mycol[12,2],mycol[12,3])
                  )

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

# swap surge.factor and subs.rate to get all the surge parameters together
s1st1.swap <- s1st1
s1st1.swap[match('subs.rate',parnames.sobol),] <- s1st1[match('surge.factor',parnames.sobol),]
s1st1.swap[match('surge.factor',parnames.sobol),] <- s1st1[match('subs.rate',parnames.sobol),]

s2.swap <- s2
s2.swap['subs.rate',] <- s2['surge.factor',]
s2.swap['surge.factor',] <- s2['subs.rate',]
s2.swap[,'subs.rate'] <- s2[,'surge.factor']
s2.swap[,'surge.factor'] <- s2[,'subs.rate']

s2_sig1.swap <- s2_sig1
s2_sig1.swap[match('subs.rate',parnames.sobol),] <- s2_sig1[match('surge.factor',parnames.sobol),]
s2_sig1.swap[match('surge.factor',parnames.sobol),] <- s2_sig1[match('subs.rate',parnames.sobol),]
s2_sig1.swap[,match('subs.rate',parnames.sobol)] <- s2_sig1[,match('surge.factor',parnames.sobol)]
s2_sig1.swap[,match('surge.factor',parnames.sobol)] <- s2_sig1[,match('subs.rate',parnames.sobol)]

plotRadCon(df=s1st1.swap
           ,s2=s2.swap
           ,scaling = .33
           ,s2_sig=s2_sig1.swap
           #,filename = '~/Box\ Sync/Wong-Projects/BRICK_scenarios/figures/sobol_spider_2065Build'
           ,filename = '~/Box\ Sync/Wong-Projects/BRICK_scenarios/figures/sobol_spider_2065Build-noHR'
           #,filename = '~/Box\ Sync/Wong-Projects/BRICK_scenarios/figures/sobol_spider_2065Build-noGEV'
           #,filename = './sobol_fig_test3'
           ,plotType = 'EPS'
           ,gpNameMult=1.5
           ,RingThick=0.1
           ,legLoc = "bottomcenter",cex = .76
           ,s1_col = rgb(mycol[3,1],mycol[3,2],mycol[3,3])
           ,st_col = rgb(mycol[6,1],mycol[6,2],mycol[6,3])
           ,line_col = rgb(mycol[10,1],mycol[10,2],mycol[10,3])
           ,STthick = 0.5
           ,legFirLabs=c(.05,.85), legTotLabs=c(.10,.90), legSecLabs=c(.02,.1)
           ,lBuildRCPhoriz=FALSE
           ,lnoGEVhoriz=FALSE
           ,lnoHRhoriz=TRUE
)

##==============================================================================






##==============================================================================
##==============================================================================
## Scratch work below here
##==============================================================================
##==============================================================================






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
