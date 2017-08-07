##==============================================================================
## Sobol sensitvity analysis for drivers of flood risk.
## --> plotting routines
##
## Code history:
## this one by Tony Wong, 6 March 2017, Penn State. Modified a bit realtive to
## the previous versions listed below.
###################################
## Adapted from 'radialPlot_vanDantzig.R'
## Originally authored by: Perry Oddo
## Pennsylvania State University
## poddo@psu.edu
###################################
## Adapted from 'radialConvergeTest.R'
## Originally authored by: Calvin Whealton
## Cornell University
## caw324@cornell.edu
####################################
## Code for radial Sobol Analysis plot
## Original code available at:
## https://github.com/calvinwhealton/SensitivityAnalysisPlots
####################################
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

control <- FALSE
noGEV <- FALSE
noHR <- TRUE
noRCP <- FALSE

plotdir <- '~/Box\ Sync/Wong-Projects/BRICK_scenarios/figures/'

if(control) {
    # Set number of parameters being analyzed
    n_params <- 41
    # Set Sobol indices file names
    Sobol_file_1 <- "../output_calibration/BRICK_Sobol-1-tot_04Aug2017-Build-AIS-GEV-2065.txt"
    Sobol_file_2 <- "../output_calibration/BRICK_Sobol-2_04Aug2017-Build-AIS-GEV-2065.txt"
} else if(noGEV) {
    n_params <- 38
    Sobol_file_1 <- "../output_calibration/BRICK_Sobol-1-tot_04Aug2017-Build-AIS-2065.txt"
    Sobol_file_2 <- "../output_calibration/BRICK_Sobol-2_04Aug2017-Build-AIS-2065.txt"
} else if(noHR) {
    n_params <- 39
    Sobol_file_1 <- "../output_calibration/BRICK_Sobol-1-tot_06Aug2017-Build-GEV-2065.txt"
    Sobol_file_2 <- "../output_calibration/BRICK_Sobol-2_06Aug2017-Build-GEV-2065.txt"
} else if(noRCP) {
    n_params <- 40
    Sobol_file_1 <- "../output_calibration/BRICK_Sobol-1-tot_26Apr2017-Build-AIS-GEV-RCP85-2065.txt"
    Sobol_file_2 <- "../output_calibration/BRICK_Sobol-2_26Apr2017-Build-AIS-GEV-RCP85-2065.txt"
}

##==============================================================================
## Radial plots


## Libraries----
library(RColorBrewer) # good color palettes
library(graphics)     # used when plotting polygons
library(plotrix)      # used when plotting circles

## Functions in other files
source('../calibration/BRICK_Sobol_functions.R')

## Import data from sensitivity analysis
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

sig.cutoff <- 0.01

# S1 & ST: using the confidence intervals
s1st1 <- stat_sig_s1st(s1st
                      ,method="congtr"
                      ,greater=sig.cutoff
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
                       ,greater=sig.cutoff
                       )

# S2: using greater than a given value
#s2_sig1 <- stat_sig_s2(s2
#                       ,s2_conf
#                       ,greater=0.02
#                       ,method='gtr')

####################################
# Define groups for the variables and the color schemes
# Defining lists of the variables for each group

if(control) {
    name_list1 <- list('Temperature' = parnames.sobol[1:5]
                   ,'Sea Level:\n   Glaciers & Ice Caps' = parnames.sobol[6:9]
                   ,'Sea Level:\nThermal Expansion' = parnames.sobol[10:13]
                   ,'Sea Level:\nGreenland Ice Sheet' = parnames.sobol[14:18]
                   ,'Sea Level:\nAntarctic Ice Sheet' = parnames.sobol[19:33]
                   ,'Sea Level:\nLand Water Storage' = parnames.sobol[34]
                   ,'Land Subsidence' = parnames.sobol[36]
                   ,'Storm Surge' = parnames.sobol[c(35,37:39)]
                   ,'Emissions' = parnames.sobol[40]
                   ,'Protection' = parnames.sobol[41]
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
                  expression(mu['AIS']), expression(nu), expression('P'[0]),
                  expression(kappa['AIS']), expression('f'[0]),
                  expression('h'[0]), 'c', expression('b'[0]), 'slope',
                  expression(lambda), expression('T'['crit']), expression('S'['LWS']),
                  expression('C'['surge']), 'subs'
                  , expression(mu), expression(sigma), expression(xi)
                  , 'RCP'
                  , 'build'
                  )

# modify the names and symbols if we are plotting one of the SOM figures
} else if(noRCP) {
    name_list1 <- list('Temperature' = parnames.sobol[1:5]
                   ,'Sea Level:\n   Glaciers & Ice Caps' = parnames.sobol[6:9]
                   ,'Sea Level:\nThermal Expansion' = parnames.sobol[10:13]
                   ,'Sea Level:\nGreenland Ice Sheet' = parnames.sobol[14:18]
                   ,'Sea Level:\nAntarctic Ice Sheet' = parnames.sobol[19:33]
                   ,'Sea Level:\nLand Water Storage' = parnames.sobol[34]
                   ,'Land Subsidence' = parnames.sobol[36]
                   ,'Storm Surge' = parnames.sobol[c(35,37:39)]
                   #,'Emissions' = parnames.sobol[40]
                   ,'Protection' = parnames.sobol[40]
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
                  expression(mu['AIS']), expression(nu), expression('P'[0]),
                  expression(kappa['AIS']), expression('f'[0]),
                  expression('h'[0]), 'c', expression('b'[0]), 'slope',
                  expression(lambda), expression('T'['crit']), expression('S'['LWS']),
                  expression('C'['surge']), 'subs'
                  , expression(mu), expression(sigma), expression(xi)
                  #, 'RCP'
                  , 'build'
                  )

} else if(noGEV) {
    name_list1 <- list('Temperature' = parnames.sobol[1:5]
                   ,'Sea Level:\n   Glaciers & Ice Caps' = parnames.sobol[6:9]
                   ,'Sea Level:\nThermal Expansion' = parnames.sobol[10:13]
                   ,'Sea Level:\nGreenland Ice Sheet' = parnames.sobol[14:18]
                   ,'Sea Level:\nAntarctic Ice Sheet' = parnames.sobol[19:33]
                   ,'Sea Level:\nLand Water Storage' = parnames.sobol[34]
                   ,'Land Subsidence' = parnames.sobol[36]
                   ,'Storm Surge' = parnames.sobol[35]
                   ,'Emissions' = parnames.sobol[37]
                   ,'Protection' = parnames.sobol[38]
                   )
    name_symbols <- c('S', expression(kappa[D]), expression(alpha[D]),
                  expression('T'[0]), expression('H'[0]), expression(beta[0]),
                  expression('V'['0,GSIC']), 'n', expression('G'['s,0']),
                  expression('a'['TE']), expression('b'['TE']),
                  expression(1/tau['TE']), expression('V'['0,TE']),
                  expression('a'['GIS']), expression('b'['GIS']),
                  expression(alpha['GIS']), expression(beta['GIS']),
                  expression('V'['0,GIS']), expression('a'['ANTO']),
                  expression('b'['ANTO']), expression(gamma), expression(alpha['AIS']),
                  expression(mu['AIS']), expression(nu), expression('P'[0]),
                  expression(kappa['AIS']), expression('f'[0]),
                  expression('h'[0]), 'c', expression('b'[0]), 'slope',
                  expression(lambda), expression('T'['crit']), expression('S'['LWS']),
                  expression('C'['surge']), 'subs'
                  , 'RCP'
                  , 'build'
                  )
} else if(noHR) {
    name_list1 <- list('Temperature' = parnames.sobol[1:5]
                   ,'Sea Level:\n   Glaciers & Ice Caps' = parnames.sobol[6:9]
                   ,'Sea Level:\nThermal Expansion' = parnames.sobol[10:13]
                   ,'Sea Level:\nGreenland Ice Sheet' = parnames.sobol[14:18]
                   ,'Sea Level:\nAntarctic Ice Sheet' = parnames.sobol[19:31]
                   ,'Sea Level:\nLand Water Storage' = parnames.sobol[32]
                   ,'Land Subsidence' = parnames.sobol[34]
                   ,'Storm Surge' = parnames.sobol[c(33,35:37)]
                   ,'Emissions' = parnames.sobol[38]
                   ,'Protection' = parnames.sobol[39]
                   )
    name_symbols <- c('S', expression(kappa[D]), expression(alpha[D]),
                  expression('T'[0]), expression('H'[0]), expression(beta[0]),
                  expression('V'['0,GSIC']), 'n', expression('G'['s,0']),
                  expression('a'['TE']), expression('b'['TE']),
                  expression(1/tau['TE']), expression('V'['0,TE']),
                  expression('a'['GIS']), expression('b'['GIS']),
                  expression(alpha['GIS']), expression(beta['GIS']),
                  expression('V'['0,GIS']), expression('a'['ANTO']),
                  expression('b'['ANTO']), expression(gamma), expression(alpha['AIS']),
                  expression(mu['AIS']), expression(nu), expression('P'[0]),
                  expression(kappa['AIS']), expression('f'[0]),
                  expression('b'[0]), 'slope',
                  expression(lambda), expression('T'['crit']), expression('S'['LWS']),
                  expression('C'['surge']), 'subs'
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
                  ,'Sea Level:\nGreenland Ice Sheet' = rgb(mycol[2,1],mycol[2,2],mycol[2,3])
                  ,'Sea Level:\nAntarctic Ice Sheet' = rgb(mycol[7,1],mycol[7,2],mycol[7,3])
                  ,'Sea Level:\nLand Water Storage'  = rgb(mycol[2,1],mycol[2,2],mycol[2,3])
                  ,'Land Subsidence' = rgb(mycol[1,1],mycol[1,2],mycol[1,3])
                  ,"Storm Surge"     = rgb(mycol[6,1],mycol[6,2],mycol[6,3])
                  ,"Emissions"       = rgb(mycol[13,1],mycol[13,2],mycol[13,3])
                  ,"Protection"      = rgb(mycol[12,1],mycol[12,2],mycol[12,3])
                  )

# using function to assign variables and colors based on group
s1st1 <- gp_name_col(name_list1
                     ,col_list1
                     ,s1st1)

s1st1$symbols <- name_symbols

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

if(control) {
    plot.filename <- paste(plotdir,'sobol_spider_2065Build',sep='')
    #plot.horiz <- c(TRUE,FALSE,FALSE)
    plot.horiz <- c(FALSE,FALSE,FALSE,TRUE)
} else if(noRCP) {
    plot.filename <- paste(plotdir,'sobol_spider_2065Build-RCP85',sep='')
    #plot.horiz <- c(FALSE,FALSE,TRUE)
    plot.horiz <- c(FALSE,FALSE,FALSE,TRUE)
} else if(noHR) {
    plot.filename <- paste(plotdir,'sobol_spider_2065Build-noHR',sep='')
    #plot.horiz <- c(FALSE,FALSE,TRUE)
    plot.horiz <- c(FALSE,FALSE,FALSE,TRUE)
} else if(noGEV) {
    plot.filename <- paste(plotdir,'sobol_spider_2065Build-noGEV',sep='')
    #plot.horiz <- c(FALSE,TRUE,FALSE)
    plot.horiz <- c(FALSE,FALSE,FALSE,TRUE)
}

plotRadCon(df=s1st1.swap
           ,s2=s2.swap
           ,scaling = .4
           ,s2_sig=s2_sig1.swap
           ,filename = plot.filename
           ,plotType = 'EPS'
           ,gpNameMult=1.5
           ,RingThick=0.1
           ,legLoc = "bottomcenter",cex = .76
           ,s1_col = rgb(mycol[3,1],mycol[3,2],mycol[3,3])
           ,st_col = rgb(mycol[6,1],mycol[6,2],mycol[6,3])
           ,line_col = rgb(mycol[10,1],mycol[10,2],mycol[10,3])
           ,STthick = 0.5
           ,legFirLabs=c(.05,.77), legTotLabs=c(.05,.83), legSecLabs=c(.02,.05)
           ,lBuildRCPhoriz=plot.horiz[1]
           ,lnoHRhoriz=plot.horiz[2]
           ,lnoGEVhoriz=plot.horiz[3]
           ,lsetback=plot.horiz[4]
)

##
## Further analysis for the text:
##

# what are the highest first-order indices?
s1.sort <- s1st[rev(order(s1st[,'S1'])),1:4]
itmp <- which(s1.sort[,'S1'] > sig.cutoff & s1.sort[,'S1_conf_low']*s1.sort[,'S1_conf_high'] > 0)
s1.sort <- s1.sort[itmp,]
print('********************************')
print('significant first-order indices:')
print(s1.sort)
print('********************************')

# what are the highest total-order indices?
st.sort <- s1st[rev(order(s1st[,'ST'])),c(1,5:7)]
itmp <- which(st.sort[,'ST'] > sig.cutoff & st.sort[,'ST_conf_low']*st.sort[,'ST_conf_high'] > 0)
st.sort <- st.sort[itmp,]
print('********************************')
print('significant total-order indices:')
print(st.sort)
print('********************************')

# what are the highest second-order interaction indices?
s2.sort <- s2_table[rev(order(s2_table[,3])),]
itmp <- which(s2.sort[,'S2'] > sig.cutoff & s2.sort[,'S2_conf_low']*s2.sort[,'S2_conf_high'] > 0)
s2.sort <- s2.sort[itmp,]
print('********************************')
print('significant second-order indices:')
print(s2.sort)
print('********************************')

##==============================================================================



##==============================================================================
## End
##==============================================================================
