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
filename.brick.magicc = '../output_model/BRICK-model_physical_control_01Nov2016.nc'
filename.brick.simple = '../output_model/BRICK-model_physical_SIMPLE-GSIC_01Nov2016.nc'
filename.brick.gmsl   = '../output_model/BRICK-model_physical_R07_01Nov2016.nc'

## File name for the Van Dantzig model output (netCDF4)
filename.vandantzig = '../output_model/VanDantzig_RCP85_control_01Nov2016.nc'

## File name for the BRICK post-calibrated parameters (csv) (the BRICK output came from these guys)
filename.parameters.magicc  = '../output_calibration/BRICK-model_postcalibratedParameters_control_01Nov2016.nc'
filename.parameters.simple  = '../output_calibration/BRICK-model_postcalibratedParameters_SIMPLE-GSIC_01Nov2016.nc'
filename.parameters.gmsl    = '../output_calibration/BRICK-model_drawcalibratedParameters_R07_01Nov2016.nc'

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
plotdir='~/Box\ Sync/Wong-Projects/BRICK_model/figures/'

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
## FIGURE 1 -- R07 MODEL FOR GMSL RISE VS BRICK
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

pdf(paste(plotdir,'gmsl_comparison.pdf',sep=''),width=4,height=6,colormodel='cmyk')

n.sig = 2         # how many sigma to plot around the obs?
itmp=midx.sl[1]:midx.sl[length(midx.sl)]

# >>> FULL BRICK <<<
par(mfrow=c(2,1), mai=c(.45,.7,.25,.08))
plot(mod.time[itmp], gmsl.ctrl.hind.mle[itmp], type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(1880,2010), ylim=c(-.2,.12), cex.lab=1.2, cex.axis=1.0, xaxs='i', yaxs='i');
  mtext(side=2, text='GMSL, BRICK-full [m SLE]', line=2.3, cex=1.0);
  mtext(side=3, text=expression(bold('  a')), line=-1.5, cex=1.0, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gmsl.ctrl.hind.95[itmp],rev(gmsl.ctrl.hind.05[itmp])), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
  lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);
	text(1930, -.125, paste("AIC = ",signif(aic.brick,4),sep=""), pos = 4, cex = 1)
  text(1930, -.150, paste("BIC = ",signif(bic.brick,4),sep=""), pos = 4, cex = 1)
  text(1930, -.175, paste("RMSE = ",signif(rmse.brick,2)," m SLE",sep=""), pos = 4, cex = 1)

  legend(1892,.13,c("2-sigma range, observations","5-95% range, model"),
         col=c(rgb(mycol[6,1],mycol[6,2],mycol[6,3]), rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3])),
         lwd=2, bty='n', cex=1.0)

# >>> BRICK-R07 <<<
par(mai=c(.65,.7,.05,.08))
plot(mod.time[itmp], gmsl.r07.hind.mle[itmp], type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(1880,2010), ylim=c(-.2,.12), cex.lab=1.2, cex.axis=1.0, xaxs='i', yaxs='i');
  mtext(side=1, text='Year', line=2.0, cex=1.0);
  mtext(side=2, text='GMSL, BRICK-R07 [m SLE]', line=2.3, cex=1.0);
  mtext(side=3, text=expression(bold('  b')), line=-1.5, cex=1.0, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gmsl.r07.hind.95[itmp],rev(gmsl.r07.hind.05[itmp])), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
  lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);
	text(1930, -.125, paste("AIC = ",signif(aic.r07,4),sep=""), pos = 4, cex = 1)
  text(1930, -.150, paste("BIC = ",signif(bic.r07,4),sep=""), pos = 4, cex = 1)
  text(1930, -.175, paste("RMSE = ",signif(rmse.r07,2)," m SLE",sep=""), pos = 4, cex = 1)

dev.off()

##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## FIGURE 2 -- GSIC HINDCASTS, MAGICC VS SIMPLE
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

pdf(paste(plotdir,'gsic_comparison.pdf',sep=''),width=4,height=6,colormodel='cmyk')

n.sig = 2         # how many sigma to plot around the obs?
itmp=midx.gsic[1]:midx.gsic[length(midx.gsic)]

# >>> GSIC-MAGICC <<<
par(mfrow=c(2,1), mai=c(.45,.7,.25,.08))
plot(mod.time[itmp], gsic.magicc.mle[itmp], type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(1960,2003), ylim=c(-.02,.04), cex.lab=1.2, cex.axis=1.0, xaxs='i', yaxs='i');
  mtext(side=2, text='GIC, MAGICC [m SLE]', line=2.3, cex=1.0);
  mtext(side=3, text=expression(bold('  a')), line=-1.5, cex=1.0, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gsic.magicc.95[itmp],rev(gsic.magicc.05[itmp])), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
  lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);
  text(1976, -0.006, paste("AIC = ",signif(aic.magicc,4),sep=""), pos = 4, cex = 1)
  text(1976, -0.011, paste("BIC = ",signif(bic.magicc,4),sep=""), pos = 4, cex = 1)
	text(1976, -0.016, paste("RMSE = ",signif(rmse.magicc,2)," m SLE",sep=""), pos = 4, cex = 1)

  legend(1964,0.0413,c("2-sigma range, observations","5-95% range, model"),
         col=c(rgb(mycol[6,1],mycol[6,2],mycol[6,3]), rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3])),
         lwd=2, bty='n', cex=1.0)

# >>> GSIC-SIMPLE <<<
par(mai=c(.65,.7,.05,.08))
plot(mod.time[itmp], gsic.simple.mle[itmp], type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(1960,2003), ylim=c(-.02,.04), cex.lab=1.2, cex.axis=1.0, xaxs='i', yaxs='i');
  mtext(side=1, text='Year', line=2.2, cex=1.0);
  mtext(side=2, text='GIC, SIMPLE [m SLE]', line=2.3, cex=1.0);
  mtext(side=3, text=expression(bold('  b')), line=-1.5, cex=1.0, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gsic.simple.95[itmp],rev(gsic.simple.05[itmp])), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
  lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);
  text(1976, -0.006, paste("AIC = ",signif(aic.simple,4),sep=""), pos = 4, cex = 1)
  text(1976, -0.011, paste("BIC = ",signif(bic.simple,4),sep=""), pos = 4, cex = 1)
	text(1976, -0.016, paste("RMSE = ",signif(rmse.simple,2)," m SLE",sep=""), pos = 4, cex = 1)

dev.off()

##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## FIGURE 3 -- PROJECTIONS OF LOCAL SLR, MAP FROM FINGERPRINTING
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

png(paste(plotdir,'slr_projections_map.tif',sep=''), width=330,height=650,units='px')

par(mfrow=c(3,1), mai=c(.5,.5,.5,1))
image(x=lon,y=lat,z=slr26.filt, zlim=lims, col=cols.b2r,
      xlab='',ylab='', xaxt='n', yaxt='n')
  map("world", add = TRUE, interior=FALSE)
  axis(1, at=seq(-180,180,by=30), labels=c('-180','','-120','','-60','','0','','60','','120','','180'), cex.axis=1.4)
  mtext('Longitude',side = 1,line = 2.2)
  axis(2, at=seq(-90,90,by=30), labels=c('-90','-60','-30','0','30','60','90'), cex.axis=1.4)
  mtext('Latitude',side=2,line=2.2)
  mtext('Global mean sea level relative to 1986-2005 average\nRCP2.6',side=3,line=0.4)
  mtext('[meters]', side=3, adj=1.3, line=-3.3)
  mtext(side=3, text=expression(bold(' a')), line=-1.5, cex=1.0, adj=0);
par(mai=c(.5,.5,.5,1))
image(x=lon,y=lat,z=slr45.filt, zlim=lims, col=cols.b2r,
      xlab='',ylab='', xaxt='n', yaxt='n')
  map("world", add = TRUE, interior=FALSE)
  axis(1, at=seq(-180,180,by=30), labels=c('-180','','-120','','-60','','0','','60','','120','','180'), cex.axis=1.4)
  mtext('Longitude',side = 1,line = 2.2)
  axis(2, at=seq(-90,90,by=30), labels=c('-90','-60','-30','0','30','60','90'), cex.axis=1.4)
  mtext('Latitude',side=2,line=2.2)
  mtext('RCP4.5',side=3,line=0.4)
  mtext(side=3, text=expression(bold(' b')), line=-1.5, cex=1.0, adj=0);
par(mai=c(.5,.5,.5,1))
image(x=lon,y=lat,z=slr85.filt, zlim=lims, col=cols.b2r,
      xlab='',ylab='', xaxt='n', yaxt='n')
  map("world", add = TRUE, interior=FALSE)
  axis(1, at=seq(-180,180,by=30), labels=c('-180','','-120','','-60','','0','','60','','120','','180'), cex.axis=1.4)
  mtext('Longitude',side = 1,line = 2.2)
  axis(2, at=seq(-90,90,by=30), labels=c('-90','-60','-30','0','30','60','90'), cex.axis=1.4)
  mtext('Latitude',side=2,line=2.2)
  mtext('RCP8.5',side=3,line=0.4)
  mtext(side=3, text=expression(bold(' c')), line=-1.5, cex=1.0, adj=0);

par(fig=c(.2,1,0,1))
image.plot(zlim=lims,legend.only=TRUE, col=cols.b2r, cex=1.4, legend.shrink = 0.85,
           axis.args=list(cex.axis=1.4,
           at=seq(-1.5,1.5,by=0.1),
           labels=c("-1.5","","","","","-1","","","","","-0.5","","","","","0",
           "","","","","0.5","","","","","1","","","","","1.5")))

dev.off()

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## FIGURE 4 -- VAN DANTZIG ANALYSIS
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

plot(heightening,cost.avg/conv,col=rgb(mycol[11,1],mycol[11,2],mycol[11,3],.2),type='l',xlim=c(0,3),ylim=c(0,8e9)/conv,
  xlab='', ylab='', xaxt='n', xaxs='i');
  mtext('Heightening [m]',side = 1,line = 2.2)
  mtext('Expected costs [billion US $]',side=2,line=2.2);
  axis(1, at=seq(0,3,by=0.5), label=c('0','','1','','2','','3'))
  axis(2, at=seq(0,8,by=2), label=rep("",5), las=1)
  polygon(c(heightening,rev(heightening)), c(cost.h[,3],rev(cost.h[,1]))/conv, col=rgb(mycol[11,1],mycol[11,2],mycol[11,3],.3), border=NA);
  lines(heightening,cost.avg/conv,type='l',col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]),lwd=3);
  points(heightening[iopt.avg],cost.avg[iopt.avg]/conv,type='p',pch=16,cex=1.3)
  arrows(Hopt.quan[1], 0.1, Hopt.quan[3], 0.1, length = 0.05, angle = 90, code = 3)
  text(0.45, 0.7, paste('90% range = ',H.text[1],'-',H.text[2],' m',sep=""), pos = 4, cex = 1)

  axis(3, at=H.rp.plot, label=parse(text=paste("10^", rp.plot, sep="")))
  mtext('Return period [years]',side = 3,line = 2.2)
  arrows(Hopt.quan[1], 7.9, Hopt.quan[3], 7.9, length = 0.05, angle = 90, code = 3)
  text(0.45, 7.1, paste('90% range = ',rp.text[1],'-',rp.text[2],' years',sep=""), pos = 4, cex = 1)

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

luse.doeclim <- TRUE
luse.gsic <- FALSE
luse.te <- FALSE
luse.simple <- FALSE
luse.dais <- FALSE
luse.gmsl	<- TRUE
luse.brick <- cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais, luse.gmsl)
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

luse.doeclim <- TRUE
luse.gsic <- TRUE
luse.te <- TRUE
luse.simple <- TRUE
luse.dais <- FALSE
luse.brick <- cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais)
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

luse.doeclim <- TRUE
luse.gsic <- TRUE
luse.te <- TRUE
luse.simple <- TRUE
luse.dais <- TRUE
luse.brick <- cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais)
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
## End
##==============================================================================
