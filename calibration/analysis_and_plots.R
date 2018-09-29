##==============================================================================
## Plots and tables for BRICK model paper (Wong et al 2016)
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
filename.brick.magicc = '../output_model/BRICK-model_physical_control_02Apr2017.nc'
filename.brick.simple = '../output_model/BRICK-model_physical_SIMPLE-GSIC_02Apr2017.nc'
filename.brick.gmsl   = '../output_model/BRICK-model_physical_R07_03Apr2017.nc'

## File name for the Van Dantzig model output (netCDF4)
filename.vandantzig = '../output_model/VanDantzig_RCP85_control_02Apr2017.nc'

## File name for the BRICK post-calibrated parameters (csv) (the BRICK output came from these guys)
filename.parameters.magicc  = '../output_calibration/BRICK-model_postcalibratedParameters_control_02Apr2017.nc'
filename.parameters.simple  = '../output_calibration/BRICK-model_postcalibratedParameters_SIMPLE-GSIC_02Apr2017.nc'
filename.parameters.gmsl    = '../output_calibration/BRICK-model_drawcalibratedParameters_R07_03Apr2017.nc'

## Other files
filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_01Nov2016.csv"
filename.fingerprints = "../fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"

## Other useful scripts
source('../Useful/colorblindPalette.R') # Get nice plotting colors: mycol array
source('../Useful/MultipleOutput.R')    # defines the useful ":=" operator

## And set the IPCC RCP colors
col26 <- c(0, 0, 255)/255
col45 <- c(121, 188, 255)/255
col60 <- c(255, 130, 45)/255
col85 <- c(255, 0, 0)/255

## Set colors to use for control model, and observations/experimental model,
## in all plots. This is indexed within "mycol", from "colorblindPalette.R".
colmod <- 2
colobs <- 11

## Where would you like to save the plots?
plotdir='~/Box\ Sync/Wong-Projects/BRICK_model/figures/test/'

##==============================================================================
logl.ar1 = function(r,sigma1,rho1,eps1=0) # default obs error is 0
{
  n = length(r) # r is the residuals
  if(length(eps1)==1) eps1 = rep(eps1,n)

	logl=0
	if(n>1) {
  	w = r[2:n] - rho1*r[1:(n-1)] # this process whitens the residuals
  	logl = logl + sum(dnorm(w,sd=sqrt((sigma1)^2+(eps1[c(-1)])^2),log=TRUE)) # add in the sum of
		    # density of the whitened residuals with a standard deviation of the
		    # variance and the obs. errors
  }
  return(logl)
}
##==============================================================================





##==============================================================================
##==============================================================================
## FIGURE 2 -- R07 MODEL FOR GMSL RISE VS BRICK
##=========

ncdata <- nc_open(filename.brick.magicc)
  t.hind = ncvar_get(ncdata, 'time_hind')
  t.proj = ncvar_get(ncdata, 'time_proj')
  gmsl.ctrl.hind = ncvar_get(ncdata, 'GlobalSeaLevel_hind')
  gmsl.ctrl.rcp26 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
  gmsl.ctrl.rcp45 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
  gmsl.ctrl.rcp85 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
nc_close(ncdata)

ncdata <- nc_open(filename.brick.gmsl)
  gmsl.r07.hind = ncvar_get(ncdata, 'GlobalSeaLevel_hind')
  gmsl.r07.rcp26 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
  gmsl.r07.rcp45 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
  gmsl.r07.rcp85 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
nc_close(ncdata)

## Initialize arrays for the output
gmsl.ctrl.hind.05 = rep(NA,length(t.hind)); gmsl.ctrl.hind.50 = rep(NA,length(t.hind)); gmsl.ctrl.hind.95 = rep(NA,length(t.hind))
gmsl.ctrl.rcp26.05 = rep(NA,length(t.proj)); gmsl.ctrl.rcp26.50 = rep(NA,length(t.proj)); gmsl.ctrl.rcp26.95 = rep(NA,length(t.proj))
gmsl.ctrl.rcp45.05 = rep(NA,length(t.proj)); gmsl.ctrl.rcp45.50 = rep(NA,length(t.proj)); gmsl.ctrl.rcp45.95 = rep(NA,length(t.proj))
gmsl.ctrl.rcp85.05 = rep(NA,length(t.proj)); gmsl.ctrl.rcp85.50 = rep(NA,length(t.proj)); gmsl.ctrl.rcp85.95 = rep(NA,length(t.proj))
gmsl.r07.hind.05 = rep(NA,length(t.hind)); gmsl.r07.hind.50 = rep(NA,length(t.hind)); gmsl.r07.hind.95 = rep(NA,length(t.hind))
gmsl.r07.rcp26.05 = rep(NA,length(t.proj)); gmsl.r07.rcp26.50 = rep(NA,length(t.proj)); gmsl.r07.rcp26.95 = rep(NA,length(t.proj))
gmsl.r07.rcp45.05 = rep(NA,length(t.proj)); gmsl.r07.rcp45.50 = rep(NA,length(t.proj)); gmsl.r07.rcp45.95 = rep(NA,length(t.proj))
gmsl.r07.rcp85.05 = rep(NA,length(t.proj)); gmsl.r07.rcp85.50 = rep(NA,length(t.proj)); gmsl.r07.rcp85.95 = rep(NA,length(t.proj))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.hind)){
	c(gmsl.ctrl.hind.05[t], gmsl.ctrl.hind.50[t], gmsl.ctrl.hind.95[t]) := quantile(gmsl.ctrl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gmsl.r07.hind.05[t], gmsl.r07.hind.50[t], gmsl.r07.hind.95[t]) := quantile(gmsl.r07.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
}
for (t in 1:length(t.proj)) {
  c(gmsl.ctrl.rcp26.05[t], gmsl.ctrl.rcp26.50[t], gmsl.ctrl.rcp26.95[t]) := quantile(gmsl.ctrl.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gmsl.ctrl.rcp45.05[t], gmsl.ctrl.rcp45.50[t], gmsl.ctrl.rcp45.95[t]) := quantile(gmsl.ctrl.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gmsl.ctrl.rcp85.05[t], gmsl.ctrl.rcp85.50[t], gmsl.ctrl.rcp85.95[t]) := quantile(gmsl.ctrl.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gmsl.r07.rcp26.05[t], gmsl.r07.rcp26.50[t], gmsl.r07.rcp26.95[t]) := quantile(gmsl.r07.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gmsl.r07.rcp45.05[t], gmsl.r07.rcp45.50[t], gmsl.r07.rcp45.95[t]) := quantile(gmsl.r07.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gmsl.r07.rcp85.05[t], gmsl.r07.rcp85.50[t], gmsl.r07.rcp85.95[t]) := quantile(gmsl.r07.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
}

begyear = t.hind[1]
endyear = t.hind[length(t.hind)]
mod.time= begyear:endyear
begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time = length(mod.time)

source('../R/compute_indices.R')				# function to determine the model and
																				# data indices for comparisons

## Source the data for hindcast comparisons
source('../calibration/DOECLIM_readData.R')
source('../calibration/GSIC_readData.R')
source('../calibration/SIMPLE_readData.R')
source('../calibration/DAIS_readData.R')
source('../calibration/TE_readData.R')

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.temp,midx.ocheat,midx.gis,midx.gsic,midx.sl)
names(midx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )
oidx.all        = list(oidx.temp,oidx.ocheat,oidx.gis,oidx.gsic,oidx.sl)
names(oidx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.temp, obs.ocheat, obs.gis, obs.gsic, obs.sl)
names(obs.all) = c(    "temp"  , "ocheat"  , "gis"  , "gsic"  , "sl" )
obs.err.all        = list( obs.temp.err, obs.ocheat.err, obs.gis.err, obs.gsic.err, obs.sl.err)
names(obs.err.all) = c(    "temp"      , "ocheat"      , "gis"      , "gsic"      , "sl"      )

## Set the indices for normalization that are consistent with each data set
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1850),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
		c(which(mod.time==1870),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

## Also normalize the available data to this time period.
ibeg=which(obs.sl.time==mod.time[ind.norm[1]])
iend=which(obs.sl.time==mod.time[ind.norm[length(ind.norm)]])
obs.sl.norm = obs.sl - mean(obs.sl[ibeg:iend])

## Get the parameters for GSIC-SIMPLE and GSIC-MAGICC, so the hindcasts can
## account for the AR1 error in the likelihood.
ncdata <- nc_open(filename.parameters.magicc)
parameters = ncvar_get(ncdata, 'BRICK_parameters')
parnames = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters.magicc = t(parameters)
colnames(parameters.magicc) = parnames

ncdata <- nc_open(filename.parameters.gmsl)
parameters = ncvar_get(ncdata, 'BRICK_parameters')
parnames = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters.gmsl = t(parameters)
colnames(parameters.gmsl) = parnames

##
## calculation of RMSE, mean bias, AIC, BIC, ... (other model comparison metrics)
##

brick.resid = mat.or.vec(length(midx.sl),ncol(gmsl.ctrl.hind))
r07.resid = mat.or.vec(length(midx.sl),ncol(gmsl.r07.hind))
llik.brick = rep(NA,ncol(gmsl.ctrl.hind))
llik.r07 = rep(NA,ncol(gmsl.r07.hind))
sigma.gmsl = parameters.gmsl[,match('sigma.gmsl',colnames(parameters.gmsl))]
rho.gmsl = parameters.gmsl[,match('rho.gmsl',colnames(parameters.gmsl))]

for (i in 1:ncol(gmsl.ctrl.hind)) {
  brick.resid[,i] = obs.all$sl[oidx.sl]-gmsl.ctrl.hind[midx.sl,i]
  llik.brick[i] = sum(dnorm(brick.resid[,i],sd=obs.sl.err[oidx.sl],log=TRUE))
}
for (i in 1:ncol(gmsl.r07.hind)) {
  r07.resid[,i] = obs.all$sl[oidx.sl]-gmsl.r07.hind[midx.sl,i]
  llik.r07[i] = logl.ar1(r07.resid[,i], sigma.gmsl[i], rho.gmsl[i], obs.err.all$sl[oidx.all$sl])
}

bic.brick = -2*max(llik.brick)+ncol(parameters.magicc)*log(length(oidx.sl))
bic.r07 = -2*max(llik.r07)+ncol(parameters.gmsl)*log(length(oidx.sl))

aic.brick = 2*ncol(parameters.magicc)-2*max(llik.brick)
aic.r07 = 2*ncol(parameters.gmsl)-2*max(llik.r07)

rmse.brick = sqrt(mean(brick.resid[,which(llik.brick==max(llik.brick))]^2))
rmse.r07 = sqrt(mean(r07.resid[,which(llik.r07==max(llik.r07))]^2))

gmsl.ctrl.hind.mle = gmsl.ctrl.hind[,which(llik.brick==max(llik.brick))]
gmsl.r07.hind.mle = gmsl.r07.hind[,which(llik.r07==max(llik.r07))]

##
## 5-95% CI of full BRICK and BRICK-R07 hindcasts, with obs
##

pdf(paste(plotdir,'gmsl_comparison.pdf',sep=''),width=3.5,height=5.25,colormodel='cmyk')

n.sig = 2         # how many sigma to plot around the obs?
itmp=midx.sl[1]:midx.sl[length(midx.sl)]

# >>> FULL BRICK <<<
par(mfrow=c(2,1), mai=c(.45,.7,.25,.08))
plot(mod.time[itmp], gmsl.ctrl.hind.mle[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(1880,2010), ylim=c(-.2,.12), cex.lab=1.2, cex.axis=1.0, xaxs='i', yaxs='i');
  mtext(side=2, text='GMSL, BRICK-full [m SLE]', line=2.3, cex=1.0);
  mtext(side=3, text=expression(bold('  a')), line=-1.5, cex=1.0, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gmsl.ctrl.hind.95[itmp],rev(gmsl.ctrl.hind.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
	text(1930, -.125, paste("AIC = ",signif(aic.brick,4),sep=""), pos = 4, cex = 1)
  text(1930, -.150, paste("BIC = ",signif(bic.brick,4),sep=""), pos = 4, cex = 1)
  text(1930, -.175, paste("RMSE = ",signif(rmse.brick,2)," m",sep=""), pos = 4, cex = 1)

  legend(1892,.13,c("2-sigma range, obs.","5-95% range, model"),
         col=c(rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3])),
         lwd=2, bty='n', cex=1.0)

# >>> BRICK-GMSL <<<
par(mai=c(.65,.7,.05,.08))
plot(mod.time[itmp], gmsl.r07.hind.mle[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(1880,2010), ylim=c(-.2,.12), cex.lab=1.2, cex.axis=1.0, xaxs='i', yaxs='i');
  mtext(side=1, text='Year', line=2.0, cex=1.0);
  mtext(side=2, text='GMSL, BRICK-GMSL [m SLE]', line=2.3, cex=1.0);
  mtext(side=3, text=expression(bold('  b')), line=-1.5, cex=1.0, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gmsl.r07.hind.95[itmp],rev(gmsl.r07.hind.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
	text(1930, -.125, paste("AIC = ",signif(aic.r07,4),sep=""), pos = 4, cex = 1)
  text(1930, -.150, paste("BIC = ",signif(bic.r07,4),sep=""), pos = 4, cex = 1)
  text(1930, -.175, paste("RMSE = ",signif(rmse.r07,2)," m",sep=""), pos = 4, cex = 1)

dev.off()

##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## FIGURE 3 -- GSIC HINDCASTS, MAGICC VS SIMPLE
##=========

ncdata <- nc_open(filename.brick.magicc)
  t.hind      = ncvar_get(ncdata, 'time_hind')
  gsic.magicc= ncvar_get(ncdata, 'GSIC_hind')
nc_close(ncdata)

ncdata <- nc_open(filename.brick.simple)
  gsic.simple = ncvar_get(ncdata, 'GSIC_hind')
nc_close(ncdata)

## Initialize arrays for the output
gsic.magicc.05 = rep(NA,length(t.hind));		gsic.magicc.50 = rep(NA,length(t.hind));		gsic.magicc.95 = rep(NA,length(t.hind))
gsic.simple.05 = rep(NA,length(t.hind));		gsic.simple.50 = rep(NA,length(t.hind));		gsic.simple.95 = rep(NA,length(t.hind))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.hind)){
	c(gsic.magicc.05[t], gsic.magicc.50[t], gsic.magicc.95[t]) := quantile(gsic.magicc[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gsic.simple.05[t], gsic.simple.50[t], gsic.simple.95[t]) := quantile(gsic.simple[t,], c(0.05,.50,.95), na.rm=TRUE)
}

begyear = t.hind[1]
endyear = t.hind[length(t.hind)]
mod.time= begyear:endyear
begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time = length(mod.time)

source('../R/compute_indices.R')				# function to determine the model and
																				# data indices for comparisons

## Source the data for hindcast comparisons
source('../calibration/DOECLIM_readData.R')
source('../calibration/GSIC_readData.R')
source('../calibration/SIMPLE_readData.R')
source('../calibration/DAIS_readData.R')
source('../calibration/TE_readData.R')

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.temp,midx.ocheat,midx.gis,midx.gsic,midx.sl)
names(midx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )
oidx.all        = list(oidx.temp,oidx.ocheat,oidx.gis,oidx.gsic,oidx.sl)
names(oidx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.temp, obs.ocheat, obs.gis, obs.gsic, obs.sl)
names(obs.all) = c(    "temp"  , "ocheat"  , "gis"  , "gsic"  , "sl" )
obs.err.all        = list( obs.temp.err, obs.ocheat.err, obs.gis.err, obs.gsic.err, obs.sl.err)
names(obs.err.all) = c(    "temp"      , "ocheat"      , "gis"      , "gsic"      , "sl"      )

## Set the indices for normalization that are consistent with each data set
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1850),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
		c(which(mod.time==1870),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

## Also normalize the available data to this time period.
ibeg=which(obs.gsic.time==mod.time[ind.norm[1]])
iend=which(obs.gsic.time==mod.time[ind.norm[length(ind.norm)]])
obs.gsic.norm = obs.gsic #- mean(obs.gsic[ibeg:iend]) # GSIC does not need normalized - already is normalized to 1960

## Get the statistical model parameters (so the hindcasts can account
## for the AR1 error in the likelihood), because this was done in the
## processingPipeline script.
ncdata <- nc_open(filename.parameters.magicc)
parameters = ncvar_get(ncdata, 'BRICK_parameters')
parnames = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters.magicc = t(parameters)
colnames(parameters.magicc) = parnames

ncdata <- nc_open(filename.parameters.simple)
parameters = ncvar_get(ncdata, 'BRICK_parameters')
parnames = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters.simple = t(parameters)
colnames(parameters.simple) = parnames

##
## calculation of RMSE, mean bias, AIC, BIC, ... (other model comparison metrics)
##

magicc.resid = mat.or.vec(length(midx.gsic),ncol(gsic.magicc))
simple.resid = mat.or.vec(length(midx.gsic),ncol(gsic.simple))
llik.magicc = rep(NA,ncol(gsic.magicc))
llik.simple = rep(NA,ncol(gsic.simple))
sigma.magicc = parameters.magicc[,match('sigma.gsic',colnames(parameters.magicc))]
sigma.simple = parameters.simple[,match('sigma.gsic',colnames(parameters.simple))]
rho.magicc = parameters.magicc[,match('rho.gsic',colnames(parameters.magicc))]
rho.simple = parameters.simple[,match('rho.gsic',colnames(parameters.simple))]

for (i in 1:ncol(gsic.magicc)) {
  magicc.resid[,i] = obs.all$gsic[oidx.gsic]-gsic.magicc[midx.gsic,i]
  llik.magicc[i] = logl.ar1(magicc.resid[,i], sigma.magicc[i], rho.magicc[i], obs.err.all$gsic[oidx.all$gsic])
}
for (i in 1:ncol(gsic.simple)) {
  simple.resid[,i] = obs.all$gsic[oidx.gsic]-gsic.simple[midx.gsic,i]
  llik.simple[i] = logl.ar1(simple.resid[,i], sigma.simple[i], rho.simple[i], obs.err.all$gsic[oidx.all$gsic])
}

bic.magicc = -2*max(llik.magicc)+6*log(length(oidx.gsic))
bic.simple = -2*max(llik.simple)+7*log(length(oidx.gsic))

aic.magicc = 2*6-2*max(llik.magicc)
aic.simple = 2*7-2*max(llik.simple)

rmse.magicc = sqrt(mean(magicc.resid[,which(llik.magicc==max(llik.magicc))]^2))
rmse.simple = sqrt(mean(simple.resid[,which(llik.simple==max(llik.simple))]^2))

gsic.magicc.mle = gsic.magicc[,which(llik.magicc==max(llik.magicc))]
gsic.simple.mle = gsic.simple[,which(llik.simple==max(llik.simple))]

##
## 5-95% CI of GSIC-MAGICC and -SIMPLE hindcasts, with obs
##

pdf(paste(plotdir,'gsic_comparison.pdf',sep=''),width=3.5,height=5.25,colormodel='cmyk')

n.sig = 2         # how many sigma to plot around the obs?
itmp=midx.gsic[1]:midx.gsic[length(midx.gsic)]

# >>> GSIC-MAGICC <<<
par(mfrow=c(2,1), mai=c(.45,.7,.25,.08))
plot(mod.time[itmp], gsic.magicc.mle[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(1960,2003), ylim=c(-.02,.05), cex.lab=1.2, cex.axis=1.0, xaxs='i', yaxs='i');
  mtext(side=2, text='GIC, MAGICC [m SLE]', line=2.3, cex=1.0);
  mtext(side=3, text=expression(bold('  a')), line=-1.5, cex=1.0, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gsic.magicc.95[itmp],rev(gsic.magicc.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  text(1976, -0.006, paste("AIC = ",signif(aic.magicc,4),sep=""), pos = 4, cex = 1)
  text(1976, -0.011, paste("BIC = ",signif(bic.magicc,4),sep=""), pos = 4, cex = 1)
	text(1976, -0.016, paste("RMSE = ",signif(rmse.magicc,2)," m",sep=""), pos = 4, cex = 1)

  legend(1964,0.051,c("2-sigma range, obs.","5-95% range, model"),
         col=c(rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3])),
         lwd=2, bty='n', cex=1.0)

# >>> GSIC-SIMPLE <<<
par(mai=c(.65,.7,.05,.08))
plot(mod.time[itmp], gsic.simple.mle[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(1960,2003), ylim=c(-.02,.05), cex.lab=1.2, cex.axis=1.0, xaxs='i', yaxs='i');
  mtext(side=1, text='Year', line=2.2, cex=1.0);
  mtext(side=2, text='GIC, SIMPLE [m SLE]', line=2.3, cex=1.0);
  mtext(side=3, text=expression(bold('  b')), line=-1.5, cex=1.0, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gsic.simple.95[itmp],rev(gsic.simple.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  text(1976, -0.007, paste("AIC = ",signif(aic.simple,4),sep=""), pos = 4, cex = 1)
  text(1976, -0.012, paste("BIC = ",signif(bic.simple,4),sep=""), pos = 4, cex = 1)
  text(1976, -0.017, paste("RMSE = ",signif(rmse.simple,2)," m",sep=""), pos = 4, cex = 1)

dev.off()

##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## FIGURE 4 -- PROJECTIONS OF LOCAL SLR, MAP FROM FINGERPRINTING
##=========

library(gplots)
library(fields)
library(maps)

ncdata <- nc_open(filename.fingerprints)
  lat = ncvar_get(ncdata, 'lat')
  lon = ncvar_get(ncdata, 'lon')
  fp.gsic = ncvar_get(ncdata, 'GLAC')
  fp.gis = ncvar_get(ncdata, 'GIS')
  fp.ais = ncvar_get(ncdata, 'AIS')
nc_close(ncdata)

## Rearrange longitudes so centered on Europe/Africa
lon = c((lon[181:360]-360),lon[1:180])
fp.gsic = rbind(fp.gsic[181:360,],fp.gsic[1:180,])
fp.gis = rbind(fp.gis[181:360,],fp.gis[1:180,])
fp.ais = rbind(fp.ais[181:360,],fp.ais[1:180,])

ncdata <- nc_open(filename.brick.magicc)
  t.proj = ncvar_get(ncdata, 'time_proj')
  gsic.rcp26 = ncvar_get(ncdata, 'GSIC_RCP26')
  gis.rcp26 = ncvar_get(ncdata, 'GIS_RCP26')
  ais.rcp26 = ncvar_get(ncdata, 'AIS_RCP26')
  te.rcp26 = ncvar_get(ncdata, 'TE_RCP26')
  gsic.rcp45 = ncvar_get(ncdata, 'GSIC_RCP45')
  gis.rcp45 = ncvar_get(ncdata, 'GIS_RCP45')
  ais.rcp45 = ncvar_get(ncdata, 'AIS_RCP45')
  te.rcp45 = ncvar_get(ncdata, 'TE_RCP45')
  gsic.rcp85 = ncvar_get(ncdata, 'GSIC_RCP85')
  gis.rcp85 = ncvar_get(ncdata, 'GIS_RCP85')
  ais.rcp85 = ncvar_get(ncdata, 'AIS_RCP85')
  te.rcp85 = ncvar_get(ncdata, 'TE_RCP85')
  gmsl.rcp26 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
  gmsl.rcp45 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
  gmsl.rcp85 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
nc_close(ncdata)

i2100 = which(t.proj==2100)
slr.rcp26 = fp.gsic*apply(gsic.rcp26,1,median)[i2100] +
            fp.gis*apply(gis.rcp26,1,median)[i2100] +
            fp.ais*apply(ais.rcp26,1,median)[i2100] +
            apply(te.rcp26,1,median)[i2100]
slr.rcp45 = fp.gsic*apply(gsic.rcp45,1,median)[i2100] +
            fp.gis*apply(gis.rcp45,1,median)[i2100] +
            fp.ais*apply(ais.rcp45,1,median)[i2100] +
            apply(te.rcp45,1,median)[i2100]
slr.rcp85 = fp.gsic*apply(gsic.rcp85,1,median)[i2100] +
            fp.gis*apply(gis.rcp85,1,median)[i2100] +
            fp.ais*apply(ais.rcp85,1,median)[i2100] +
            apply(te.rcp85,1,median)[i2100]

itmp=which(!is.na(slr.rcp26)); lims26=round(c(quantile(slr.rcp26[itmp],.025), quantile(slr.rcp26[itmp],.975))*100)/100
itmp=which(!is.na(slr.rcp45)); lims45=round(c(quantile(slr.rcp45[itmp],.025), quantile(slr.rcp45[itmp],.975))*100)/100
itmp=which(!is.na(slr.rcp85)); lims85=round(c(quantile(slr.rcp85[itmp],.025), quantile(slr.rcp85[itmp],.975))*100)/100


##================================================
## Results for mid-latitudes, tropics and arctic
##================================================
i.arctic = which(lat>60)
i.midlat = which( (lat< -30 & lat>-60) | (lat>30 & lat< 60) )
i.tropic = which( lat<30 & lat> -30)

tmp = slr.rcp26[,i.arctic]; itmp=which(!is.na(tmp)); arctic.26=median(tmp[itmp])
tmp = slr.rcp26[,i.midlat]; itmp=which(!is.na(tmp)); midlat.26=median(tmp[itmp])
tmp = slr.rcp26[,i.tropic]; itmp=which(!is.na(tmp)); tropic.26=median(tmp[itmp])

tmp = slr.rcp45[,i.arctic]; itmp=which(!is.na(tmp)); arctic.45=median(tmp[itmp])
tmp = slr.rcp45[,i.midlat]; itmp=which(!is.na(tmp)); midlat.45=median(tmp[itmp])
tmp = slr.rcp45[,i.tropic]; itmp=which(!is.na(tmp)); tropic.45=median(tmp[itmp])

tmp = slr.rcp85[,i.arctic]; itmp=which(!is.na(tmp)); arctic.85=median(tmp[itmp])
tmp = slr.rcp85[,i.midlat]; itmp=which(!is.na(tmp)); midlat.85=median(tmp[itmp])
tmp = slr.rcp85[,i.tropic]; itmp=which(!is.na(tmp)); tropic.85=median(tmp[itmp])

gmsl.26 = median(gmsl.rcp26[i2100,])
gmsl.45 = median(gmsl.rcp45[i2100,])
gmsl.85 = median(gmsl.rcp85[i2100,])

## Report results to the screen
print(paste('>>>>>>>>>> 2100 medians under RCP2.6 <<<<<<<<<<<<<'))
print(paste('GMSL = ',gmsl.26,' m',sep=""))
print(paste('tropic SLR = ',tropic.26,' m',sep=""))
print(paste('midlat SLR = ',midlat.26,' m',sep=""))
print(paste('arctic SLR = ',arctic.26,' m',sep=""))
print(paste('>>>>>>>>>> 2100 medians under RCP4.5 <<<<<<<<<<<<<'))
print(paste('GMSL = ',gmsl.45,' m',sep=""))
print(paste('tropic SLR = ',tropic.45,' m',sep=""))
print(paste('midlat SLR = ',midlat.45,' m',sep=""))
print(paste('arctic SLR = ',arctic.45,' m',sep=""))
print(paste('>>>>>>>>>> 2100 medians under RCP8.5 <<<<<<<<<<<<<'))
print(paste('GMSL = ',gmsl.85,' m',sep=""))
print(paste('tropic SLR = ',tropic.85,' m',sep=""))
print(paste('midlat SLR = ',midlat.85,' m',sep=""))
print(paste('arctic SLR = ',arctic.85,' m',sep=""))

##================================================


plotme = slr.rcp85
slr26.filt = slr.rcp26
slr45.filt = slr.rcp45
slr85.filt = slr.rcp85
lims = c(-1.5,1.5)

## Filtering the SLR values outside of the plotting range
ifilt.hi = which(slr26.filt>lims[2]); slr26.filt[ifilt.hi] = lims[2]; ifilt.lo = which(slr26.filt<lims[1]); slr26.filt[ifilt.lo] = lims[1];
ifilt.hi = which(slr45.filt>lims[2]); slr45.filt[ifilt.hi] = lims[2]; ifilt.lo = which(slr45.filt<lims[1]); slr45.filt[ifilt.lo] = lims[1];
ifilt.hi = which(slr85.filt>lims[2]); slr85.filt[ifilt.hi] = lims[2]; ifilt.lo = which(slr85.filt<lims[1]); slr85.filt[ifilt.lo] = lims[1];

frac.neg = abs(lims[1])/diff(lims)
frac.pos = abs(lims[2])/diff(lims)
ncols = 300
n.pos = floor(ncols*frac.pos)
n.neg = ncols-n.pos
#cols.pos = colorRampPalette(c("#FF2222","white"),space="Lab")(n.pos)
cols.pos = colorRampPalette(c("red","orange","yellow","white"),space="Lab")(n.pos)
cols.neg = colorRampPalette(c("white","blue"),space="Lab")(n.neg)
cols.b2r = rev(c(cols.pos,cols.neg))

#png(paste(plotdir,'slr_projections_map.tif',sep=''), width=330,height=650,res=300,units='px')
png(paste(plotdir,'slr_projections_map.tif',sep=''), width=3.5,height=6.89,res=300,units='in')
#pdf(paste(plotdir,'slr_projections_map.pdf',sep=''), width=3.5,height=8,colormodel='cmyk')

par(mfrow=c(3,1), mai=c(.2,.5,.7,.9))
image(x=lon,y=lat,z=slr26.filt, zlim=lims, col=cols.b2r,
      xlab='',ylab='', xaxt='n', yaxt='n')
map("world", add = TRUE, interior=FALSE)
axis(1, at=seq(-180,180,by=45), labels=c('-180','','-90','','0','','90','','180'), cex.axis=1.4)
mtext('Longitude',side = 1,line = 2.2)
axis(2, at=seq(-90,90,by=30), labels=c('-90','-60','-30','0','30','60','90'), cex.axis=1.4)
mtext('Latitude',side=2,line=2.2)
mtext('Global mean sea level\nrelative to 1986-2005 average\nRCP2.6',side=3,line=0.4)
mtext('[meters]', side=3, adj=1.5, line=-1.5)
mtext(side=3, text=expression(bold(' a')), line=-1.2, cex=1.0, adj=0);
par(mai=c(.35,.5,.55,.9))
image(x=lon,y=lat,z=slr45.filt, zlim=lims, col=cols.b2r,
      xlab='',ylab='', xaxt='n', yaxt='n')
map("world", add = TRUE, interior=FALSE)
axis(1, at=seq(-180,180,by=45), labels=c('-180','','-90','','0','','90','','180'), cex.axis=1.4)
mtext('Longitude',side = 1,line = 2.2)
axis(2, at=seq(-90,90,by=30), labels=c('-90','-60','-30','0','30','60','90'), cex.axis=1.4)
mtext('Latitude',side=2,line=2.2)
mtext('RCP4.5',side=3,line=0.4)
mtext(side=3, text=expression(bold(' b')), line=-1.2, cex=1.0, adj=0);
par(mai=c(.5,.5,.4,.9))
image(x=lon,y=lat,z=slr85.filt, zlim=lims, col=cols.b2r,
      xlab='',ylab='', xaxt='n', yaxt='n')
map("world", add = TRUE, interior=FALSE)
axis(1, at=seq(-180,180,by=45), labels=c('-180','','-90','','0','','90','','180'), cex.axis=1.4)
mtext('Longitude',side = 1,line = 2.2)
axis(2, at=seq(-90,90,by=30), labels=c('-90','-60','-30','0','30','60','90'), cex.axis=1.4)
mtext('Latitude',side=2,line=2.2)
mtext('RCP8.5',side=3,line=0.4)
mtext(side=3, text=expression(bold(' c')), line=-1.2, cex=1.0, adj=0);

par(fig=c(.2,1,0,1), mai=c(.3,2.5,.3,.1))
image.plot(zlim=lims,legend.only=TRUE, col=cols.b2r, cex=1.4, legend.shrink = 0.75,
           axis.args=list(cex.axis=1.4,
                          at=seq(-1.5,1.5,by=0.1),
                          labels=c("-1.5","","","","","-1","","","","","-0.5","","","","","0",
                                   "","","","","0.5","","","","","1","","","","","1.5")))

dev.off()

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## FIGURE 5 -- VAN DANTZIG ANALYSIS
##=========

## Grab the van Dantzig output
ncdata <- nc_open(filename.vandantzig)
  heightening     = ncvar_get(ncdata, 'H')
  cost            = ncvar_get(ncdata, 'ExpectedCost')
  loss            = ncvar_get(ncdata, 'ExpectedLoss')
  investment      = ncvar_get(ncdata, 'ExpectedInvestment')
  preturn         = ncvar_get(ncdata, 'ExpectedPreturn')
nc_close(ncdata)

n.ensemble = ncol(cost)
n.height = nrow(cost)

## What is the optimal heightening (minimum expected costs) for each ensemble
## member? (and the index, so we can grab the return period)
iopt=rep(NA,n.ensemble)						# index of optimal heightening, return period
rp.opt=rep(NA,n.ensemble)		      # return period of optimal heightening (including fast dynamics)
Hopt=rep(NA,n.ensemble)           # optimal heightening, accounting for fast dynamics

for (i in 1:n.ensemble) {
	iopt[i] = which(cost[,i]==min(cost[,i]))
  Hopt[i] = heightening[iopt[i]]
  rp.opt[i] = preturn[iopt[i],i]
}

Hopt.quan = c(quantile(Hopt, 0.05), quantile(Hopt, 0.50), quantile(Hopt, 0.95))
rp.opt.quan = c(quantile(rp.opt, 0.05), quantile(rp.opt, 0.50), quantile(rp.opt, 0.95))

H.text = round(Hopt.quan[c(1,3)],2)
rp.text = round(rp.opt.quan[c(1,3)],0)

preturn.avg = apply(preturn, 1, mean)
cost.avg    = apply(cost, 1, mean)
loss.avg    = apply(loss, 1, mean)
investment.avg    = apply(investment, 1, mean)
iopt.avg = which(cost.avg==min(cost.avg))

preturn.pre = 100
#preturn.pre = round(preturn.nofd.avg[iopt.nofd.avg])

ipre = rep(NA,n.ensemble)
for (i in 1:n.ensemble) {
  ipre[i] = which(preturn[,i]>=preturn.pre)[1]
}
ipre.avg = which(preturn.avg >= preturn.pre)[1]

tmp=lm( investment.avg ~ log10(preturn.avg))
inter = tmp$coefficients[[1]]
slope = tmp$coefficients[[2]]

## Calculate ensemble median and 5-95% range for costs as a function of
## (a) heightening, and (b) return period
cost.h = mat.or.vec(length(heightening),3)

for (t in 1:length(heightening)) {
  cost.h[t,1] = quantile(cost[t,],.05); cost.h[t,2] = quantile(cost[t,],.50); cost.h[t,3] = quantile(cost[t,],.95);
}

## Return period is trickier. Need to bin, because not linear.
rp.bin = seq(from=log10(min(preturn)), to=log10(max(preturn)), length.out=100)
cost.rp = mat.or.vec(length(rp.bin),3)
for (t in 1:length(rp.bin)) {
  itmp=which(log10(preturn)>rp.bin[t] & log10(preturn)<rp.bin[t+1])
  cost.rp[t,1] = quantile(cost[itmp],0.05, na.rm=TRUE); cost.rp[t,2] = quantile(cost[itmp],0.50, na.rm=TRUE); cost.rp[t,3] = quantile(cost[itmp],0.95, na.rm=TRUE);
}
itrim=which(!is.na(cost.rp[,1]) & !is.na(cost.rp[,2]) & !is.na(cost.rp[,3]))
rp.bin=rp.bin[itrim]
cost.rp = cost.rp[itrim,]

rp.plot = c(2,3,4)
H.rp.plot = rep(0,length(rp.plot))
for (i in 1:length(H.rp.plot)) {
  tmp = abs(preturn.avg-(10^rp.plot[i]))
  itmp = which(tmp==min(tmp))
  H.rp.plot[i] = heightening[itmp]
}

print('====================================================================')
print(paste('rec. heightening =',heightening[iopt.avg]/.3048,'ft'))
print('====================================================================')
print(paste('rec. heightening =',heightening[iopt.avg],'m'))
print('====================================================================')


## BOTH RETURN PERIOD AND HEIGHTENING, SAME COST AXIS

pdf(paste(plotdir,'vandantzig_RP+H.pdf',sep=''),width=3.5,height=3.5,colormodel='cmyk')

conv=1e9  # convert from $ to billions or millions of $? (for nicer looking axes)

## RETURN PERIOD ALONG TOP, HEIGHTENING ALONG BOTTOM

par(mfrow=c(1,1), mai=c(.65,.65,.65,.1))

plot(heightening,cost.avg/conv,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.2),type='l',xlim=c(0,3),ylim=c(0,8e9)/conv,
  xlab='', ylab='', xaxt='n', xaxs='i');
  mtext('Heightening [m]',side = 1,line = 2.2)
  mtext('Expected costs [billion US $]',side=2,line=2.2);
  axis(1, at=seq(0,3,by=0.5), label=c('0','','1','','2','','3'))
  axis(2, at=seq(0,8,by=2), label=rep("",5), las=1)
  polygon(c(heightening,rev(heightening)), c(cost.h[,3],rev(cost.h[,1]))/conv, col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.3), border=NA);
  lines(heightening,cost.avg/conv,type='l',col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]),lwd=3);
  points(heightening[iopt.avg],cost.avg[iopt.avg]/conv,type='p',pch=16,cex=1.3)
  arrows(Hopt.quan[1], 0.1, Hopt.quan[3], 0.1, length = 0.05, angle = 90, code = 3)
  text(0.48, 0.7, paste('90% range = ',H.text[1],'-',H.text[2],' m',sep=""), pos = 4, cex = 1)

  axis(3, at=H.rp.plot, label=parse(text=paste("10^", rp.plot, sep="")))
  mtext('Return period [years]',side = 3,line = 2.2)
  arrows(Hopt.quan[1], 7.9, Hopt.quan[3], 7.9, length = 0.05, angle = 90, code = 3)
  text(0.48, 7.1, paste('90% range = ',rp.text[1],'-',rp.text[2],' years',sep=""), pos = 4, cex = 1)

dev.off()

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## APPENDIX 1 -- PRIOR AND POSTERIOR PARAMETER RANGES
##===========

source('../Useful/MultipleOutput.R') # defines the ":=" operator

##==============================================================================
## The GMSL emulator experiment

## First, get prior ranges from each BRICK_parameterSetup.R version (for the
## three model experiments). Do the control one last because it should not be
## overwritten. These comments are not repeated for the SIMPLE-GSIC or Control
## experiments.

luse.sneasy <- FALSE
luse.doeclim <- TRUE
luse.gsic <- FALSE
luse.te <- FALSE
luse.tee <- FALSE
luse.simple <- FALSE
luse.dais <- FALSE
luse.gmsl	<- TRUE
luse.brick <- cbind(luse.sneasy,luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais, luse.gmsl)
source('BRICK_parameterSetup_R07.R')
ind.gmsl <- rep(NA,length(parnames.gmsl))
for (i in 1:length(parnames.gmsl)) {
  ind.gmsl[i] <- which(parnames==parnames.gmsl[i])
}
bounds.gmsl <- cbind(bound.lower.gmsl, bound.upper.gmsl)

## Now ind.gmsl, parnames.gmsl and bounds.gmsl hold the indices within parameters.gmsl
## and the prior ranges for the GMSL model parameters.

## Second, get the posterior parameters.
ncdata <- nc_open(filename.parameters.gmsl)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
nc_close(ncdata)
parameters.gmsl <- parameters[ind.gmsl,]
rownames(parameters.gmsl) <- parnames.gmsl

post.gmsl <- mat.or.vec(length(parnames.gmsl),3)
for (i in 1:length(parnames.gmsl)) {
  post.gmsl[i,] <- signif( quantile(parameters.gmsl[i,],c(.05,.5,.95)) , 2)
}

rownames(post.gmsl) <- parnames.gmsl
colnames(post.gmsl) <- c('5%','Median','95%')

table.gmsl <- cbind(bounds.gmsl, post.gmsl)
write.csv(table.gmsl, "../output_calibration/parameters.gmsl.csv")
##==============================================================================

##==============================================================================
## The SIMPLE-GIC experiment

luse.sneasy <- FALSE
luse.doeclim <- TRUE
luse.gsic <- TRUE
luse.te <- TRUE
luse.tee <- FALSE
luse.simple <- TRUE
luse.dais <- FALSE
luse.brick <- cbind(luse.sneasy,luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais)
source('BRICK_parameterSetup_SIMPLE-GSIC.R')
ind.gsic <- rep(NA,length(parnames.gsic))
for (i in 1:length(parnames.gsic)) {ind.gsic[i] <- which(parnames==parnames.gsic[i])}
bounds.gsic <- cbind(bound.lower.gsic, bound.upper.gsic)
ncdata <- nc_open(filename.parameters.simple)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
nc_close(ncdata)
parameters.gsic <- parameters[ind.gsic,]
rownames(parameters.gsic) <- parnames.gsic
post.gsic <- mat.or.vec(length(parnames.gsic),3)
for (i in 1:length(parnames.gsic)) {post.gsic[i,] <- signif( quantile(parameters.gsic[i,],c(.05,.5,.95)) , 2)}
rownames(post.gsic) <- parnames.gsic
colnames(post.gsic) <- c('5%','Median','95%')
table.gsic <- cbind(bounds.gsic, post.gsic)
write.csv(table.gsic, "../output_calibration/parameters.gic-simple.csv")
##==============================================================================

##==============================================================================
## The Control experiment

## For this one, need *all* of the models' parameters

luse.sneasy <- FALSE
luse.doeclim <- TRUE
luse.gsic <- TRUE
luse.te <- TRUE
luse.tee <- FALSE
luse.simple <- TRUE
luse.dais <- TRUE
luse.brick <- cbind(luse.sneasy, luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais)
source('BRICK_parameterSetup.R')

## DOECLIM

bounds.doeclim <- cbind(bound.lower.doeclim, bound.upper.doeclim)
ncdata <- nc_open(filename.parameters.magicc)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
ind.doeclim <- rep(NA,length(parnames.doeclim))
for (i in 1:length(parnames.doeclim)) {ind.doeclim[i] <- which(parnames==parnames.doeclim[i])}
parameters.doeclim <- parameters[ind.doeclim,]
rownames(parameters.doeclim) <- parnames.doeclim
post.doeclim <- mat.or.vec(length(parnames.doeclim),3)
for (i in 1:length(parnames.doeclim)) {post.doeclim[i,] <- signif( quantile(parameters.doeclim[i,],c(.05,.5,.95)) , 2)}
rownames(post.doeclim) <- parnames.doeclim
colnames(post.doeclim) <- c('5%','Median','95%')
table.doeclim <- cbind(bounds.doeclim, post.doeclim)
write.csv(table.doeclim, "../output_calibration/parameters.doeclim.csv")

## GIC-MAGICC

bounds.magicc <- cbind(bound.lower.gsic, bound.upper.gsic)
ncdata <- nc_open(filename.parameters.magicc)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
ind.magicc <- rep(NA,length(parnames.gsic))
for (i in 1:length(parnames.gsic)) {ind.magicc[i] <- which(parnames==parnames.gsic[i])}
parameters.magicc <- parameters[ind.magicc,]
rownames(parameters.magicc) <- parnames.gsic
post.magicc <- mat.or.vec(length(parnames.gsic),3)
for (i in 1:length(parnames.gsic)) {post.magicc[i,] <- signif( quantile(parameters.magicc[i,],c(.05,.5,.95)) , 2)}
rownames(post.magicc) <- parnames.gsic
colnames(post.magicc) <- c('5%','Median','95%')
table.magicc <- cbind(bounds.magicc, post.magicc)
write.csv(table.magicc, "../output_calibration/parameters.gic-magicc.csv")

## GIS

bounds.simple <- cbind(bound.lower.simple, bound.upper.simple)
ncdata <- nc_open(filename.parameters.magicc)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
ind.simple <- rep(NA,length(parnames.simple))
for (i in 1:length(parnames.simple)) {ind.simple[i] <- which(parnames==parnames.simple[i])}
parameters.simple <- parameters[ind.simple,]
rownames(parameters.simple) <- parnames.simple
post.simple <- mat.or.vec(length(parnames.simple),3)
for (i in 1:length(parnames.simple)) {post.simple[i,] <- signif( quantile(parameters.simple[i,],c(.05,.5,.95)) , 2)}
rownames(post.simple) <- parnames.simple
colnames(post.simple) <- c('5%','Median','95%')
table.simple <- cbind(bounds.simple, post.simple)
write.csv(table.simple, "../output_calibration/parameters.gis.csv")

## TE

bounds.te <- cbind(bound.lower.te, bound.upper.te)
ncdata <- nc_open(filename.parameters.magicc)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
ind.te <- rep(NA,length(parnames.te))
for (i in 1:length(parnames.te)) {ind.te[i] <- which(parnames==parnames.te[i])}
parameters.te <- parameters[ind.te,]
rownames(parameters.te) <- parnames.te
post.te <- mat.or.vec(length(parnames.te),3)
for (i in 1:length(parnames.te)) {post.te[i,] <- signif( quantile(parameters.te[i,],c(.05,.5,.95)) , 2)}
rownames(post.te) <- parnames.te
colnames(post.te) <- c('5%','Median','95%')
table.te <- cbind(bounds.te, post.te)
write.csv(table.te, "../output_calibration/parameters.te.csv")

## DAIS

bounds.dais <- cbind(bound.lower.dais, bound.upper.dais)
ncdata <- nc_open(filename.parameters.magicc)
parameters <- ncvar_get(ncdata, 'BRICK_parameters')
parnames <- ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
ind.dais <- rep(NA,length(parnames.dais))
for (i in 1:length(parnames.dais)) {ind.dais[i] <- which(parnames==parnames.dais[i])}
parameters.dais <- parameters[ind.dais,]
rownames(parameters.dais) <- parnames.dais
post.dais <- mat.or.vec(length(parnames.dais),3)
for (i in 1:length(parnames.dais)) {post.dais[i,] <- signif( quantile(parameters.dais[i,],c(.05,.5,.95)) , 2)}
rownames(post.dais) <- parnames.dais
colnames(post.dais) <- c('5%','Median','95%')
table.dais <- cbind(bounds.dais, post.dais)
write.csv(table.dais, "../output_calibration/parameters.dais.csv")
##==============================================================================


##==============================================================================
##==============================================================================
## Supplemental Figure
## -- Hindcast, model vs observational data

ncdata <- nc_open(filename.brick.magicc)
  t.hind      = ncvar_get(ncdata, 'time_hind')
  gsl.hind    = ncvar_get(ncdata, 'GlobalSeaLevel_hind')
  temp.hind   = ncvar_get(ncdata, 'temp_hind')
  ocheat.hind = ncvar_get(ncdata, 'ocheat_hind')
  gis.hind    = ncvar_get(ncdata, 'GIS_hind')
  te.hind     = ncvar_get(ncdata, 'TE_hind')
  gsic.hind   = ncvar_get(ncdata, 'GSIC_hind')
  t.paleo     = ncvar_get(ncdata, 'time_paleo_avg')
  ais.paleo.05= ncvar_get(ncdata, 'AIS_paleo_avg_q05')
  ais.paleo.50= ncvar_get(ncdata, 'AIS_paleo_avg_q50')
  ais.paleo.95= ncvar_get(ncdata, 'AIS_paleo_avg_q95')
nc_close(ncdata)

## Initialize arrays for the output
slr.05 = rep(NA,length(t.hind));		slr.50 = rep(NA,length(t.hind));		slr.95 = rep(NA,length(t.hind))
gsic.05 = rep(NA,length(t.hind));		gsic.50 = rep(NA,length(t.hind));		gsic.95 = rep(NA,length(t.hind))
gis.05 = rep(NA,length(t.hind));		gis.50 = rep(NA,length(t.hind));		gis.95 = rep(NA,length(t.hind))
te.05 = rep(NA,length(t.hind));			te.50 = rep(NA,length(t.hind));			te.95 = rep(NA,length(t.hind))
temp.05 = rep(NA,length(t.hind));		temp.50 = rep(NA,length(t.hind));		temp.95 = rep(NA,length(t.hind))
ocheat.05 = rep(NA,length(t.hind)); ocheat.50 = rep(NA,length(t.hind));	ocheat.95 = rep(NA,length(t.hind))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.hind)){
	c(slr.05[t], slr.50[t], slr.95[t])					:= quantile(gsl.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
	c(gsic.05[t], gsic.50[t], gsic.95[t])				:= quantile(gsic.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
	c(gis.05[t], gis.50[t], gis.95[t])					:= quantile(gis.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
	c(te.05[t], te.50[t], te.95[t])							:= quantile(te.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
	c(temp.05[t], temp.50[t], temp.95[t])				:= quantile(temp.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
	c(ocheat.05[t], ocheat.50[t], ocheat.95[t])	:= quantile(ocheat.hind[t,], c(0.05,.50,.95), na.rm=TRUE)
}

begyear = t.hind[1]
endyear = t.hind[length(t.hind)]
mod.time= begyear:endyear
begyear.norm = 1961
endyear.norm = 1990
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time = length(mod.time)

source('../R/compute_indices.R')  # function to determine the model and
                                  # data indices for comparisons

## Source the data for hindcast comparisons
source('../calibration/DOECLIM_readData.R')
source('../calibration/GSIC_readData.R')
source('../calibration/SIMPLE_readData.R')
source('../calibration/DAIS_readData.R')
source('../calibration/TE_readData.R')

## Gather up all the data/model indices for comparisons. use lists to avoid
## enormous amounts of input to the MCMC functions
midx.all        = list(midx.temp,midx.ocheat,midx.gis,midx.gsic,midx.sl)
names(midx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )
oidx.all        = list(oidx.temp,oidx.ocheat,oidx.gis,oidx.gsic,oidx.sl)
names(oidx.all) = c(   "temp"   ,"ocheat"   ,"gis"   ,"gsic"   ,"sl"   )

## Gather up all the observations for comparisons
obs.all        = list( obs.temp, obs.ocheat, obs.gis, obs.gsic, obs.sl)
names(obs.all) = c(    "temp"  , "ocheat"  , "gis"  , "gsic"  , "sl" )
obs.err.all        = list( obs.temp.err, obs.ocheat.err, obs.gis.err, obs.gsic.err, obs.sl.err)
names(obs.err.all) = c(    "temp"      , "ocheat"      , "gis"      , "gsic"      , "sl"      )

## Set the indices for normalization that are consistent with each data set
ind.norm.data = data.frame(
		c( "temp"              , "ocheat"            , "gsic"             , "gis"               , "te"                 , "ais"               , "sl"                ) ,
		c(which(mod.time==1850),which(mod.time==1960),which(mod.time==1960),which(mod.time==1960),which(mod.time==1961),which(mod.time==1961),which(mod.time==1961)) ,
		c(which(mod.time==1870),which(mod.time==1990),which(mod.time==1960),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990),which(mod.time==1990)) )

## Also normalize the available data to this time period.
ibeg=which(obs.temp.time==mod.time[1])
iend=which(obs.temp.time==mod.time[20])
obs.temp.norm = obs.temp - mean(obs.temp[ibeg:iend])

ibeg=which(obs.ocheat.time==mod.time[ind.norm[1]])
iend=which(obs.ocheat.time==mod.time[ind.norm[length(ind.norm)]])
obs.ocheat.norm = obs.ocheat - mean(obs.ocheat[ibeg:iend])

ibeg=which(obs.gsic.time==mod.time[ind.norm[1]])
iend=which(obs.gsic.time==mod.time[ind.norm[length(ind.norm)]])
obs.gsic.norm = obs.gsic #- mean(obs.gsic[ibeg:iend]) # GSIC does not need normalized - already is normalized to 1960

ibeg=which(obs.gis.time==mod.time[ind.norm[1]])
iend=which(obs.gis.time==mod.time[ind.norm[length(ind.norm)]])
obs.gis.norm = obs.gis - mean(obs.gis[ibeg:iend])

ibeg=which(obs.sl.time==mod.time[ind.norm[1]])
iend=which(obs.sl.time==mod.time[ind.norm[length(ind.norm)]])
obs.sl.norm = obs.sl - mean(obs.sl[ibeg:iend])

## TE trends
# read in TE_readData.R -- "trends.te"

## DAIS calibration windows

trends.ais = c(0.27 , 0.08 , 0.40 )/1000   # m/year (after the /1000)
trends.err = c(0.11 , 0.185, 0.205)/1000   # m/year (after the /1000)
trends.2up = trends.ais+2*trends.err
trends.2dn = trends.ais-2*trends.err
ind.trends = mat.or.vec( length(trends.ais), 2)
ind.trends[1,] = c(which(date==-7) , which(date==10)) # 1993-2010
ind.trends[2,] = c(which(date==-8) , which(date== 1)) # 1992-2001
ind.trends[3,] = c(which(date== 2) , which(date==11)) # 2002-2011

## Precal window 4:
## Adding observational constraint
estimate.SLE.rate = abs(-71/360)/1000
time.years = 2002-1992      # using the midpoint of the 19-year interval
mid.cum.SLE_2002 = estimate.SLE.rate*time.years
i1992 = which(date==-8)

estimate.SLE.rate.error = abs(-53/360)/1000     #1-sigma error
estimate.SLE.error = sqrt(time.years)*estimate.SLE.rate.error #1-sigma error
        # (*sqrt(10) because 10 years of potentially accumulated error:
        #  total error^2 = year 1 error^2 + year 2 error^2 + ... year 10 error^2
        #                = 10*year X error^2)
SE2_2002 = estimate.SLE.error*2 #2-sigma error

positive_2SE = mid.cum.SLE_2002 + SE2_2002 # Add the 2 standard error to the mean value
negative_2SE = mid.cum.SLE_2002 - SE2_2002 # Subtract the 2 standard error to the mean value

## Precal windows 1-3:
## from Shaffer (2014). modified by Kelsey
#upper.wind = c(7.4, -6.9, -1.25, positive_2SE) # Windows 1-3 from Kelsey
#lower.wind = c(3.6, -15.8, -4.0, negative_2SE)
upper.wind = c(6.0, -6.9, -1.25, positive_2SE) # Windows 1-3 from Kelsey
lower.wind = c(1.8, -15.8, -4.0, negative_2SE)
#upper.wind = c(5.5, -8 , -2, positive_2SE) # Windows 1-3 fFrom Shaffer 2014, p 1809
#lower.wind = c(2.5, -17, -4, negative_2SE)

windows = matrix(c(lower.wind, upper.wind), nrow = length(upper.wind), ncol=2)
obs.targets = (windows[,2]+windows[,1])*.5    # middle of window = obs to compare model to
obs.err = (windows[,2]-windows[,1])*.5       # half-width of window = uncertainty
obs.err=0.5*obs.err                           # assume all windows are 2*stdErr
                                              # (last two actually are)
## Create a vector with each observation year
## 120kyr, 20Kyr, 6kyr (before present), 2002, and 1993 (first year of the IPCC
## trend)
obs.years = c(120000, 220000, 234000, 240002)

##
## 5-95% CI of hindcasts, with obs in there
##

pdf(paste(plotdir,'hindcasts.pdf',sep=''),width=7,height=6,colormodel='cmyk')
n.sig = 2         # how many sigma to plot around the obs?
layout(cbind(c(1,4,7),c(2,5,7),c(3,6,7)))
par(mai=c(.5,.7,.2,.08))

# >>> SURFACE TEMPERATURE <<<
plot(mod.time[midx.temp], temp.50[midx.temp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(1850,2016), ylim=c(-.3,1.5), cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Surface temperature\n[deg C]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' a')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time[midx.temp],rev(mod.time[midx.temp])), c(temp.95,rev(temp.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(obs.temp.time[oidx.temp], obs.temp.norm[oidx.temp], type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.temp.time[oidx.temp],rev(obs.temp.time[oidx.temp])),
					c(obs.temp.norm[oidx.temp]+n.sig*obs.temp.err[oidx.temp],rev(obs.temp.norm[oidx.temp]-n.sig*obs.temp.err[oidx.temp])),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
  #legend(1839,1.65,c("5-95% range, modeled","2-sigma range, observations"),
  #       col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]),rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])), lwd=2, bty='n', cex=1.2)

# >>> OCEAN HEAT <<<
itmp=midx.ocheat[1]:nrow(ocheat.hind)
plot(mod.time[itmp], ocheat.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-50,50), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Ocean heat uptake\n[10^22 J]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(ocheat.95[itmp],rev(ocheat.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(obs.ocheat.time, obs.ocheat.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.ocheat.time,rev(obs.ocheat.time)), c(obs.ocheat.norm+n.sig*obs.ocheat.err,rev(obs.ocheat.norm-n.sig*obs.ocheat.err)),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);

# >>> GSIC <<<
itmp=midx.gsic[1]:nrow(gsic.hind)
plot(mod.time[itmp], gsic.50[itmp], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.01,.04), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Glaciers and\nsmall ice caps [m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' c')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gsic.95[itmp],rev(gsic.05[itmp])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);

# >>> GIS <<<
plot(mod.time, gis.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(gis.95,rev(gis.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);

# >>> TE <<<
x1971=seq(trends.te[1,4],trends.te[1,5])
c1971=mean(x1971); yc1971=te.50[which(mod.time==mean(x1971))]
lo1971 = yc1971+(trends.te[1,2]/1000)*(x1971-c1971)
hi1971 = yc1971+(trends.te[1,3]/1000)*(x1971-c1971)
y1971 = yc1971+(trends.te[1,1]/1000)*(x1971-c1971)
x1993=seq(trends.te[2,4],trends.te[2,5])
c1993=mean(x1993); yc1993=te.50[which(mod.time==mean(x1993))]
lo1993 = yc1993+(trends.te[2,2]/1000)*(x1993-c1993)
hi1993 = yc1993+(trends.te[2,3]/1000)*(x1993-c1993)
y1993 = yc1993+(trends.te[2,1]/1000)*(x1993-c1993)

plot(mod.time, te.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.04,.04), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Thermal expansion\n[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(te.95,rev(te.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
	lines(x1971,y1971, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
  polygon(c(x1971,rev(x1971)), c(lo1971,rev(hi1971)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);
	lines(x1993,y1993, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
	polygon(c(x1993,rev(x1993)), c(lo1993,rev(hi1993)), col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);

# >>> TOTAL SLR <<<
plot(mod.time, slr.50, type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' f')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(slr.95,rev(slr.05)), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
  lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
          col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.5), border=NA);

# >>> AIS PALEO, SMOOTHED <<<
ipaleo=which(t.paleo==-149999):which(t.paleo==1)
plot(t.paleo[ipaleo], ais.paleo.50[ipaleo], type='l', col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]), lwd=2, xlab='',
     ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]), ylim=c(-20,10),
		 cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' g')), line=.25, cex=.9, adj=0);
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(ais.paleo.95[ipaleo],rev(ais.paleo.05[ipaleo])), col=rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3],.5), border=NA);
	for (i in 1:3) {
  	polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
	 			  	c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
          	col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3],.7), border=NA);
	}
	i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3]))
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  legend(-92000,12, c("5-95% range, model" , "2-sigma range, observations"), lwd=2, bty='n', cex=1.2,
         col=c(rgb(mycol[colmod,1],mycol[colmod,2],mycol[colmod,3]) , rgb(mycol[colobs,1],mycol[colobs,2],mycol[colobs,3])) )

dev.off()
##==============================================================================
##==============================================================================


##==============================================================================
##==============================================================================
## Supplemental Figure
## -- Sea-level projections to 2100

ncdata <- nc_open(filename.brick.magicc)
  t.proj = ncvar_get(ncdata, 'time_proj')
  slr.rcp26 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
  slr.rcp45 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
  slr.rcp85 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
  te.rcp26 = ncvar_get(ncdata, 'TE_RCP26')
  te.rcp45 = ncvar_get(ncdata, 'TE_RCP45')
  te.rcp85 = ncvar_get(ncdata, 'TE_RCP85')
  gis.rcp26 = ncvar_get(ncdata, 'GIS_RCP26')
  gis.rcp45 = ncvar_get(ncdata, 'GIS_RCP45')
  gis.rcp85 = ncvar_get(ncdata, 'GIS_RCP85')
  gsic.rcp26 = ncvar_get(ncdata, 'GSIC_RCP26')
  gsic.rcp45 = ncvar_get(ncdata, 'GSIC_RCP45')
  gsic.rcp85 = ncvar_get(ncdata, 'GSIC_RCP85')
  ais.rcp26 = ncvar_get(ncdata, 'AIS_RCP26')
  ais.rcp45 = ncvar_get(ncdata, 'AIS_RCP45')
  ais.rcp85 = ncvar_get(ncdata, 'AIS_RCP85')
  temp.rcp26 = ncvar_get(ncdata, 'temp_RCP26')
  temp.rcp45 = ncvar_get(ncdata, 'temp_RCP45')
  temp.rcp85 = ncvar_get(ncdata, 'temp_RCP85')
  ocheat.rcp26 = ncvar_get(ncdata, 'ocheat_RCP26')
  ocheat.rcp45 = ncvar_get(ncdata, 'ocheat_RCP45')
  ocheat.rcp85 = ncvar_get(ncdata, 'ocheat_RCP85')
nc_close(ncdata)


## Initialize arrays for the output
slr.rcp26.05 = rep(NA,length(t.proj)); slr.rcp26.50 = rep(NA,length(t.proj)); slr.rcp26.95 = rep(NA,length(t.proj))
slr.rcp45.05 = rep(NA,length(t.proj)); slr.rcp45.50 = rep(NA,length(t.proj)); slr.rcp45.95 = rep(NA,length(t.proj))
slr.rcp85.05 = rep(NA,length(t.proj)); slr.rcp85.50 = rep(NA,length(t.proj)); slr.rcp85.95 = rep(NA,length(t.proj))
ais.rcp26.05 = rep(NA,length(t.proj)); ais.rcp26.50 = rep(NA,length(t.proj)); ais.rcp26.95 = rep(NA,length(t.proj))
ais.rcp45.05 = rep(NA,length(t.proj)); ais.rcp45.50 = rep(NA,length(t.proj)); ais.rcp45.95 = rep(NA,length(t.proj))
ais.rcp85.05 = rep(NA,length(t.proj)); ais.rcp85.50 = rep(NA,length(t.proj)); ais.rcp85.95 = rep(NA,length(t.proj))
gis.rcp26.05 = rep(NA,length(t.proj)); gis.rcp26.50 = rep(NA,length(t.proj)); gis.rcp26.95 = rep(NA,length(t.proj))
gis.rcp45.05 = rep(NA,length(t.proj)); gis.rcp45.50 = rep(NA,length(t.proj)); gis.rcp45.95 = rep(NA,length(t.proj))
gis.rcp85.05 = rep(NA,length(t.proj)); gis.rcp85.50 = rep(NA,length(t.proj)); gis.rcp85.95 = rep(NA,length(t.proj))
gsic.rcp26.05 = rep(NA,length(t.proj)); gsic.rcp26.50 = rep(NA,length(t.proj)); gsic.rcp26.95 = rep(NA,length(t.proj))
gsic.rcp45.05 = rep(NA,length(t.proj)); gsic.rcp45.50 = rep(NA,length(t.proj)); gsic.rcp45.95 = rep(NA,length(t.proj))
gsic.rcp85.05 = rep(NA,length(t.proj)); gsic.rcp85.50 = rep(NA,length(t.proj)); gsic.rcp85.95 = rep(NA,length(t.proj))
te.rcp26.05 = rep(NA,length(t.proj)); te.rcp26.50 = rep(NA,length(t.proj)); te.rcp26.95 = rep(NA,length(t.proj))
te.rcp45.05 = rep(NA,length(t.proj)); te.rcp45.50 = rep(NA,length(t.proj)); te.rcp45.95 = rep(NA,length(t.proj))
te.rcp85.05 = rep(NA,length(t.proj)); te.rcp85.50 = rep(NA,length(t.proj)); te.rcp85.95 = rep(NA,length(t.proj))
temp.rcp26.05 = rep(NA,length(t.proj)); temp.rcp26.50 = rep(NA,length(t.proj)); temp.rcp26.95 = rep(NA,length(t.proj))
temp.rcp45.05 = rep(NA,length(t.proj)); temp.rcp45.50 = rep(NA,length(t.proj)); temp.rcp45.95 = rep(NA,length(t.proj))
temp.rcp85.05 = rep(NA,length(t.proj)); temp.rcp85.50 = rep(NA,length(t.proj)); temp.rcp85.95 = rep(NA,length(t.proj))
ocheat.rcp26.05 = rep(NA,length(t.proj)); ocheat.rcp26.50 = rep(NA,length(t.proj)); ocheat.rcp26.95 = rep(NA,length(t.proj))
ocheat.rcp45.05 = rep(NA,length(t.proj)); ocheat.rcp45.50 = rep(NA,length(t.proj)); ocheat.rcp45.95 = rep(NA,length(t.proj))
ocheat.rcp85.05 = rep(NA,length(t.proj)); ocheat.rcp85.50 = rep(NA,length(t.proj)); ocheat.rcp85.95 = rep(NA,length(t.proj))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.proj)){
  c(slr.rcp26.05[t], slr.rcp26.50[t], slr.rcp26.95[t]) := quantile(slr.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(slr.rcp45.05[t], slr.rcp45.50[t], slr.rcp45.95[t]) := quantile(slr.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(slr.rcp85.05[t], slr.rcp85.50[t], slr.rcp85.95[t]) := quantile(slr.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(ais.rcp26.05[t], ais.rcp26.50[t], ais.rcp26.95[t]) := quantile(ais.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(ais.rcp45.05[t], ais.rcp45.50[t], ais.rcp45.95[t]) := quantile(ais.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(ais.rcp85.05[t], ais.rcp85.50[t], ais.rcp85.95[t]) := quantile(ais.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gis.rcp26.05[t], gis.rcp26.50[t], gis.rcp26.95[t]) := quantile(gis.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gis.rcp45.05[t], gis.rcp45.50[t], gis.rcp45.95[t]) := quantile(gis.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gis.rcp85.05[t], gis.rcp85.50[t], gis.rcp85.95[t]) := quantile(gis.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gsic.rcp26.05[t], gsic.rcp26.50[t], gsic.rcp26.95[t]) := quantile(gsic.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gsic.rcp45.05[t], gsic.rcp45.50[t], gsic.rcp45.95[t]) := quantile(gsic.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(gsic.rcp85.05[t], gsic.rcp85.50[t], gsic.rcp85.95[t]) := quantile(gsic.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(te.rcp26.05[t], te.rcp26.50[t], te.rcp26.95[t]) := quantile(te.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(te.rcp45.05[t], te.rcp45.50[t], te.rcp45.95[t]) := quantile(te.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(te.rcp85.05[t], te.rcp85.50[t], te.rcp85.95[t]) := quantile(te.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(temp.rcp26.05[t], temp.rcp26.50[t], temp.rcp26.95[t]) := quantile(temp.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(temp.rcp45.05[t], temp.rcp45.50[t], temp.rcp45.95[t]) := quantile(temp.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(temp.rcp85.05[t], temp.rcp85.50[t], temp.rcp85.95[t]) := quantile(temp.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(ocheat.rcp26.05[t], ocheat.rcp26.50[t], ocheat.rcp26.95[t]) := quantile(ocheat.rcp26[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(ocheat.rcp45.05[t], ocheat.rcp45.50[t], ocheat.rcp45.95[t]) := quantile(ocheat.rcp45[t,], c(0.05,.50,.95), na.rm=TRUE)
  c(ocheat.rcp85.05[t], ocheat.rcp85.50[t], ocheat.rcp85.95[t]) := quantile(ocheat.rcp85[t,], c(0.05,.50,.95), na.rm=TRUE)
}

iproj = which(t.proj==2000):which(t.proj==2100)
i2100 = which(t.proj==2100)

print('==============================================================')
print('min/5%/50%/95%/max of 2100 sea level relative to 1986-2005:')
print(paste('RCP2.6: ',quantile(slr.rcp26[251,],c(0,.05,.50,.95,1))))
print(paste('RCP4.5: ',quantile(slr.rcp45[251,],c(0,.05,.50,.95,1))))
print(paste('RCP8.5: ',quantile(slr.rcp85[251,],c(0,.05,.50,.95,1))))
print('==============================================================')

i2050 <- which(t.proj==2050)
print('==============================================================')
print('min/5%/50%/95%/max of 2050 sea level relative to 1986-2005:')
print(paste('RCP2.6: ',quantile(slr.rcp26[i2050,],c(0,.05,.50,.95,1))))
print(paste('RCP4.5: ',quantile(slr.rcp45[i2050,],c(0,.05,.50,.95,1))))
print(paste('RCP8.5: ',quantile(slr.rcp85[i2050,],c(0,.05,.50,.95,1))))
print('==============================================================')

##==============================================================================

##
## 5-95% CI and median of 2100 SLR and all contributions under the RCP scenarios
##

slr.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
 													c(slr.rcp26.50[i2100],slr.rcp45.50[i2100],slr.rcp85.50[i2100]),
													c(slr.rcp26.05[i2100],slr.rcp45.05[i2100],slr.rcp85.05[i2100]),
													c(slr.rcp26.95[i2100],slr.rcp45.95[i2100],slr.rcp85.95[i2100])
													))
gsic.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
 													c(gsic.rcp26.50[i2100],gsic.rcp45.50[i2100],gsic.rcp85.50[i2100]),
													c(gsic.rcp26.05[i2100],gsic.rcp45.05[i2100],gsic.rcp85.05[i2100]),
													c(gsic.rcp26.95[i2100],gsic.rcp45.95[i2100],gsic.rcp85.95[i2100])
													))
gis.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
 													c(gis.rcp26.50[i2100],gis.rcp45.50[i2100],gis.rcp85.50[i2100]),
													c(gis.rcp26.05[i2100],gis.rcp45.05[i2100],gis.rcp85.05[i2100]),
													c(gis.rcp26.95[i2100],gis.rcp45.95[i2100],gis.rcp85.95[i2100])
													))
te.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
 													c(te.rcp26.50[i2100],te.rcp45.50[i2100],te.rcp85.50[i2100]),
													c(te.rcp26.05[i2100],te.rcp45.05[i2100],te.rcp85.05[i2100]),
													c(te.rcp26.95[i2100],te.rcp45.95[i2100],te.rcp85.95[i2100])
													))
ais.ci90 = data.frame( cbind( c('RCP2.6','RCP4.5','RCP8.5'),
 													c(ais.rcp26.50[i2100],ais.rcp45.50[i2100],ais.rcp85.50[i2100]),
													c(ais.rcp26.05[i2100],ais.rcp45.05[i2100],ais.rcp85.05[i2100]),
													c(ais.rcp26.95[i2100],ais.rcp45.95[i2100],ais.rcp85.95[i2100])
													))
## get into format for Latex table
row.gsic = paste('Glaciers and small ice caps &',
									1000*signif(gsic.rcp26.50[i2100],4),' (',1000*signif(gsic.rcp26.05[i2100],4),'-',1000*signif(gsic.rcp26.95[i2100],4),') &',
									1000*signif(gsic.rcp45.50[i2100],4),' (',1000*signif(gsic.rcp45.05[i2100],4),'-',1000*signif(gsic.rcp45.95[i2100],4),') &',
									1000*signif(gsic.rcp85.50[i2100],4),' (',1000*signif(gsic.rcp85.05[i2100],4),'-',1000*signif(gsic.rcp85.95[i2100],4),') \\',
									sep='')
row.te = paste('Thermal expansion &',
									1000*signif(te.rcp26.50[i2100],4),' (',1000*signif(te.rcp26.05[i2100],4),'-',1000*signif(te.rcp26.95[i2100],4),') &',
									1000*signif(te.rcp45.50[i2100],4),' (',1000*signif(te.rcp45.05[i2100],4),'-',1000*signif(te.rcp45.95[i2100],4),') &',
									1000*signif(te.rcp85.50[i2100],4),' (',1000*signif(te.rcp85.05[i2100],4),'-',1000*signif(te.rcp85.95[i2100],4),') \\',
									sep='')
row.gis = paste('Greenland Ice Sheet &',
									1000*signif(gis.rcp26.50[i2100],4),' (',1000*signif(gis.rcp26.05[i2100],4),'-',1000*signif(gis.rcp26.95[i2100],4),') &',
									1000*signif(gis.rcp45.50[i2100],4),' (',1000*signif(gis.rcp45.05[i2100],4),'-',1000*signif(gis.rcp45.95[i2100],4),') &',
									1000*signif(gis.rcp85.50[i2100],4),' (',1000*signif(gis.rcp85.05[i2100],4),'-',1000*signif(gis.rcp85.95[i2100],4),') \\',
									sep='')
row.ais = paste('Antarctic Ice Sheet &',
									1000*signif(ais.rcp26.50[i2100],4),' (',1000*signif(ais.rcp26.05[i2100],4),'-',1000*signif(ais.rcp26.95[i2100],4),') &',
									1000*signif(ais.rcp45.50[i2100],4),' (',1000*signif(ais.rcp45.05[i2100],4),'-',1000*signif(ais.rcp45.95[i2100],4),') &',
									1000*signif(ais.rcp85.50[i2100],4),' (',1000*signif(ais.rcp85.05[i2100],4),'-',1000*signif(ais.rcp85.95[i2100],4),') \\',
									sep='')
row.slr = paste('Total sea level &',
									1000*signif(slr.rcp26.50[i2100],4),' (',1000*signif(slr.rcp26.05[i2100],4),'-',1000*signif(slr.rcp26.95[i2100],4),') &',
									1000*signif(slr.rcp45.50[i2100],4),' (',1000*signif(slr.rcp45.05[i2100],4),'-',1000*signif(slr.rcp45.95[i2100],4),') &',
									1000*signif(slr.rcp85.50[i2100],4),' (',1000*signif(slr.rcp85.05[i2100],4),'-',1000*signif(slr.rcp85.95[i2100],4),') \\',
									sep='')

pdf(paste(plotdir,'projections_SLR_total.pdf',sep=''),width=3.5,height=2.45,colormodel='cmyk')
par(mfrow=c(1,1))
# RCP85
par(mai=c(.65,.65,.20,.2))
plot(t.proj[iproj],slr.rcp85.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2100), ylim=c(0,2), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2100,by=20)); axis(2, seq(0,2,by=.25), lab=c('0','','0.5','','1','','1.5','','2'));
		 mtext(side=2, text='Global mean sea level [m]', line=2.2, cex=1);
   mtext(side=1, text='Year', line=2.2, cex=1);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp85.95[iproj],rev(slr.rcp85.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + UNIFORM RCP45
  lines(t.proj[iproj],slr.rcp45.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp45.95[iproj],rev(slr.rcp45.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
# + UNIFORM RCP26
  lines(t.proj[iproj],slr.rcp26.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp26.95[iproj],rev(slr.rcp26.05[iproj])),
          col=rgb(col26[1],col26[2],col26[3],.5), border=NA);
# + legend
  legend(t.proj[iproj[1]],2,c("5-95% range,","RCP2.6","RCP4.5","RCP8.5"),
         lty=c(NA,1,1,1), lwd=3, bty='n', cex=1,
         col=c(NA , rgb(col26[1],col26[2],col26[3]) ,
               rgb(col45[1],col45[2],col45[3]) , rgb(col85[1],col85[2],col85[3])))

dev.off()
##==============================================================================
##==============================================================================


##==============================================================================
## End
##==============================================================================
