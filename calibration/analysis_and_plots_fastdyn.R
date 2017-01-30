##==============================================================================
## plots and tables for fast dynamics deep uncertainty paper
##
## Directory references are assuming code is run in AISfastdynamics/calibration/
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

## Initial set-up
library(ncdf4)

## File name for the BRICK physical model output (netCDF4)
filename.brick.uniform  = '../output_model/BRICK-fastdyn_physical_uniform_29Jan2017.nc'
filename.brick.gamma    = '../output_model/BRICK-fastdyn_physical_gamma_29Jan2017.nc'

## File name for the Van Dantzig model output (netCDF4)
filename.vandantzig.uniform = '../output_model/vanDantzig_RCP85_uniform_29Jan2017.nc'
filename.vandantzig.gamma   = '../output_model/vanDantzig_RCP85_gamma_29Jan2017.nc'

## File name for the BRICK post-calibrated parameters (csv) (the BRICK and van Dantzig output came from these guys)
filename.parameters.uniform = '../output_calibration/BRICK-fastdyn_postcalibratedParameters_uniform_29Jan2017.csv'
filename.parameters.gamma   = '../output_calibration/BRICK-fastdyn_postcalibratedParameters_gamma_29Jan2017.csv'

## Other files
filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_06Sep2016.csv"
filename.daisnofd = '../output_model/DAISfastdyn_noFD-paleo-ensemble_08Sep2016.nc'
filename.DAIScalibration = "../output_calibration/DAISfastdyn_precalibrationMCMC_parameters_gamma_21Aug2016.csv"

## Get nice plotting colors: mycol array
source('../Useful/colorblindPalette.R')

## And set the IPCC RCP colors
col26 <- c(0, 0, 255)/255
col45 <- c(121, 188, 255)/255
col60 <- c(255, 130, 45)/255
col85 <- c(255, 0, 0)/255

## Where would you like to save the plots?
plotdir='~/Box\ Sync/Wong-Projects/AIS_fast_dynamics/ToSubmit_CC/figures/'
#plotdir='../'

##==============================================================================
##==============================================================================
## EXTENDED DATA FIGURE 1 -- HINDCASTS VS OBSERVATIONS
##=======================

ncdata <- nc_open(filename.brick.gamma)
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

ncdata <- nc_open(filename.daisnofd)
  ais.nofd.05= ncvar_get(ncdata, 'AIS_paleo_avg_q05')
  ais.nofd.50= ncvar_get(ncdata, 'AIS_paleo_avg_q50')
  ais.nofd.95= ncvar_get(ncdata, 'AIS_paleo_avg_q95')
  ais.nofd.min= ncvar_get(ncdata, 'AIS_paleo_avg_min')
  ais.nofd.max= ncvar_get(ncdata, 'AIS_paleo_avg_max')
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
## from Shaffer (2014). modified by Ruckert et al (2016)
upper.wind = c(7.4, -6.9, -1.25, positive_2SE) # Windows 1-3 from Kelsey
lower.wind = c(3.6, -15.8, -4.0, negative_2SE)

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

pdf(paste(plotdir,'hindcasts_tmp.pdf',sep=''),width=7,height=6,colormodel='cmyk')
n.sig = 2         # how many sigma to plot around the obs?
layout(cbind(c(1,4,7),c(2,5,7),c(3,6,7)))
par(mai=c(.5,.7,.2,.08))

# >>> SURFACE TEMPERATURE <<<
plot(mod.time[midx.temp], temp.50[midx.temp], type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, xlab='',
     ylab='', xlim=c(1850,2016), ylim=c(-.3,1.5), cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Surface temperature\n[deg C]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' a')), line=-1.5, cex=.9, adj=0);
  polygon(c(mod.time[midx.temp],rev(mod.time[midx.temp])), c(temp.95,rev(temp.05)), col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
  lines(obs.temp.time[oidx.temp], obs.temp.norm[oidx.temp], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.temp.time[oidx.temp],rev(obs.temp.time[oidx.temp])),
					c(obs.temp.norm[oidx.temp]+n.sig*obs.temp.err[oidx.temp],rev(obs.temp.norm[oidx.temp]-n.sig*obs.temp.err[oidx.temp])),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);
  #legend(1839,1.65,c("5-95% range, modeled","2-sigma range, observations"),
  #       col=c(rgb(col45[1],col45[2],col45[3]),rgb(mycol[6,1],mycol[6,2],mycol[6,3])), lwd=2, bty='n', cex=1.2)

# >>> OCEAN HEAT <<<
plot(mod.time, ocheat.50, type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-50,50), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Ocean heat uptake\n[10^22 J]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b')), line=-1.5, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(ocheat.95,rev(ocheat.05)), col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
  lines(obs.ocheat.time, obs.ocheat.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.ocheat.time,rev(obs.ocheat.time)), c(obs.ocheat.norm+n.sig*obs.ocheat.err,rev(obs.ocheat.norm-n.sig*obs.ocheat.err)),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);

# >>> GSIC <<<
plot(mod.time, gsic.50, type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.01,.04), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Glaciers and\nsmall ice caps [m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' c')), line=-1.5, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(gsic.95,rev(gsic.05)), col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
  lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);

# >>> GIS <<<
plot(mod.time, gis.50, type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d')), line=-1.5, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(gis.95,rev(gis.05)), col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
  lines(obs.gis.time, obs.gis.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gis.time,rev(obs.gis.time)), c(obs.gis.norm+n.sig*sqrt(obs.gis.err^2),rev(obs.gis.norm-n.sig*sqrt(obs.gis.err^2))),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);

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

plot(mod.time, te.50, type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.04,.04), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Thermal expansion\n[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=-1.5, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(te.95,rev(te.05)), col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
	lines(x1971,y1971, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
  polygon(c(x1971,rev(x1971)), c(lo1971,rev(hi1971)), col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);
	lines(x1993,y1993, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
	polygon(c(x1993,rev(x1993)), c(lo1993,rev(hi1993)), col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);

# >>> TOTAL SLR <<<
plot(mod.time, slr.50, type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, xlab='',
     ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' f')), line=-1.5, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(slr.95,rev(slr.05)), col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
  lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);

# >>> AIS PALEO, SMOOTHED <<<
ipaleo=which(t.paleo==-149999):which(t.paleo==1)
plot(t.paleo[ipaleo], ais.paleo.50[ipaleo], type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, xlab='',
     ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]), ylim=c(-25,25),
		 cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' g')), line=-1.5, cex=.9, adj=0);
  lines(t.paleo[ipaleo], ais.nofd.50[ipaleo], type='l', col=rgb(mycol[12,1],mycol[12,2],mycol[12,3]), lwd=2);
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(ais.nofd.max[ipaleo],rev(ais.nofd.min[ipaleo])), col=rgb(mycol[13,1],mycol[13,2],mycol[13,3],.5), border=NA);
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(ais.paleo.95[ipaleo],rev(ais.paleo.05[ipaleo])), col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
	for (i in 1:3) {
  	polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
	 			  	c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
          	col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.7), border=NA);
	}
	i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]))
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  legend(-90000,29,c("5-95% range, model, with fast dynamics","2-sigma range, observations","Max. likelihood ensemble range, no fast dynamics"),
         col=c(rgb(col45[1],col45[2],col45[3]),rgb(mycol[6,1],mycol[6,2],mycol[6,3]),rgb(col26[1],col26[2],col26[3])), lwd=2, bty='n', cex=1.2)

dev.off()

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## FIGURE 2 -- FAST DYNAMICS PROJECTIONS (DIFFERENT RCP SCENARIOS)
##=========

ncdata <- nc_open(filename.brick.uniform)
  t.proj = ncvar_get(ncdata, 'time_proj')
  gsl.rcp26.uni = ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
  gsl.rcp45.uni = ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
  gsl.rcp85.uni = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
  gsl.nofd.rcp26.uni = ncvar_get(ncdata, 'GlobalSeaLevel_nofd_RCP26')
  gsl.nofd.rcp45.uni = ncvar_get(ncdata, 'GlobalSeaLevel_nofd_RCP45')
  gsl.nofd.rcp85.uni = ncvar_get(ncdata, 'GlobalSeaLevel_nofd_RCP85')
nc_close(ncdata)
ncdata <- nc_open(filename.brick.gamma)
  gsl.rcp26.gam = ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
  gsl.rcp45.gam = ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
  gsl.rcp85.gam = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
  gsl.nofd.rcp26.gam = ncvar_get(ncdata, 'GlobalSeaLevel_nofd_RCP26')
  gsl.nofd.rcp45.gam = ncvar_get(ncdata, 'GlobalSeaLevel_nofd_RCP45')
  gsl.nofd.rcp85.gam = ncvar_get(ncdata, 'GlobalSeaLevel_nofd_RCP85')
nc_close(ncdata)


## Initialize arrays for the output
gsl.rcp26.uni.05 = rep(NA,length(t.proj)); gsl.rcp26.uni.50 = rep(NA,length(t.proj)); gsl.rcp26.uni.95 = rep(NA,length(t.proj))
gsl.rcp45.uni.05 = rep(NA,length(t.proj)); gsl.rcp45.uni.50 = rep(NA,length(t.proj)); gsl.rcp45.uni.95 = rep(NA,length(t.proj))
gsl.rcp85.uni.05 = rep(NA,length(t.proj)); gsl.rcp85.uni.50 = rep(NA,length(t.proj)); gsl.rcp85.uni.95 = rep(NA,length(t.proj))
gsl.nofd.rcp26.uni.05 = rep(NA,length(t.proj)); gsl.nofd.rcp26.uni.50 = rep(NA,length(t.proj)); gsl.nofd.rcp26.uni.95 = rep(NA,length(t.proj))
gsl.nofd.rcp45.uni.05 = rep(NA,length(t.proj)); gsl.nofd.rcp45.uni.50 = rep(NA,length(t.proj)); gsl.nofd.rcp45.uni.95 = rep(NA,length(t.proj))
gsl.nofd.rcp85.uni.05 = rep(NA,length(t.proj)); gsl.nofd.rcp85.uni.50 = rep(NA,length(t.proj)); gsl.nofd.rcp85.uni.95 = rep(NA,length(t.proj))
gsl.rcp26.gam.05 = rep(NA,length(t.proj)); gsl.rcp26.gam.50 = rep(NA,length(t.proj)); gsl.rcp26.gam.95 = rep(NA,length(t.proj))
gsl.rcp45.gam.05 = rep(NA,length(t.proj)); gsl.rcp45.gam.50 = rep(NA,length(t.proj)); gsl.rcp45.gam.95 = rep(NA,length(t.proj))
gsl.rcp85.gam.05 = rep(NA,length(t.proj)); gsl.rcp85.gam.50 = rep(NA,length(t.proj)); gsl.rcp85.gam.95 = rep(NA,length(t.proj))
gsl.nofd.rcp26.gam.05 = rep(NA,length(t.proj)); gsl.nofd.rcp26.gam.50 = rep(NA,length(t.proj)); gsl.nofd.rcp26.gam.95 = rep(NA,length(t.proj))
gsl.nofd.rcp45.gam.05 = rep(NA,length(t.proj)); gsl.nofd.rcp45.gam.50 = rep(NA,length(t.proj)); gsl.nofd.rcp45.gam.95 = rep(NA,length(t.proj))
gsl.nofd.rcp85.gam.05 = rep(NA,length(t.proj)); gsl.nofd.rcp85.gam.50 = rep(NA,length(t.proj)); gsl.nofd.rcp85.gam.95 = rep(NA,length(t.proj))
fastdyn.rcp26.uni.05 = rep(NA,length(t.proj)); fastdyn.rcp26.uni.50 = rep(NA,length(t.proj)); fastdyn.rcp26.uni.95 = rep(NA,length(t.proj))
fastdyn.rcp45.uni.05 = rep(NA,length(t.proj)); fastdyn.rcp45.uni.50 = rep(NA,length(t.proj)); fastdyn.rcp45.uni.95 = rep(NA,length(t.proj))
fastdyn.rcp85.uni.05 = rep(NA,length(t.proj)); fastdyn.rcp85.uni.50 = rep(NA,length(t.proj)); fastdyn.rcp85.uni.95 = rep(NA,length(t.proj))
fastdyn.rcp26.gam.05 = rep(NA,length(t.proj)); fastdyn.rcp26.gam.50 = rep(NA,length(t.proj)); fastdyn.rcp26.gam.95 = rep(NA,length(t.proj))
fastdyn.rcp45.gam.05 = rep(NA,length(t.proj)); fastdyn.rcp45.gam.50 = rep(NA,length(t.proj)); fastdyn.rcp45.gam.95 = rep(NA,length(t.proj))
fastdyn.rcp85.gam.05 = rep(NA,length(t.proj)); fastdyn.rcp85.gam.50 = rep(NA,length(t.proj)); fastdyn.rcp85.gam.95 = rep(NA,length(t.proj))

source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
for (t in 1:length(t.proj)){
    c(gsl.rcp26.uni.05[t], gsl.rcp26.uni.50[t], gsl.rcp26.uni.95[t]) := quantile(gsl.rcp26.uni[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(gsl.rcp45.uni.05[t], gsl.rcp45.uni.50[t], gsl.rcp45.uni.95[t]) := quantile(gsl.rcp45.uni[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(gsl.rcp85.uni.05[t], gsl.rcp85.uni.50[t], gsl.rcp85.uni.95[t]) := quantile(gsl.rcp85.uni[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(gsl.nofd.rcp26.uni.05[t], gsl.nofd.rcp26.uni.50[t], gsl.nofd.rcp26.uni.95[t]) := quantile(gsl.nofd.rcp26.uni[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(gsl.nofd.rcp45.uni.05[t], gsl.nofd.rcp45.uni.50[t], gsl.nofd.rcp45.uni.95[t]) := quantile(gsl.nofd.rcp45.uni[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(gsl.nofd.rcp85.uni.05[t], gsl.nofd.rcp85.uni.50[t], gsl.nofd.rcp85.uni.95[t]) := quantile(gsl.nofd.rcp85.uni[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(gsl.rcp26.gam.05[t], gsl.rcp26.gam.50[t], gsl.rcp26.gam.95[t]) := quantile(gsl.rcp26.gam[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(gsl.rcp45.gam.05[t], gsl.rcp45.gam.50[t], gsl.rcp45.gam.95[t]) := quantile(gsl.rcp45.gam[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(gsl.rcp85.gam.05[t], gsl.rcp85.gam.50[t], gsl.rcp85.gam.95[t]) := quantile(gsl.rcp85.gam[t,], c(0.05,.50,.95), na.rm=TRUE)
	c(gsl.nofd.rcp26.gam.05[t], gsl.nofd.rcp26.gam.50[t], gsl.nofd.rcp26.gam.95[t]) := quantile(gsl.nofd.rcp26.gam[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(gsl.nofd.rcp45.gam.05[t], gsl.nofd.rcp45.gam.50[t], gsl.nofd.rcp45.gam.95[t]) := quantile(gsl.nofd.rcp45.gam[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(gsl.nofd.rcp85.gam.05[t], gsl.nofd.rcp85.gam.50[t], gsl.nofd.rcp85.gam.95[t]) := quantile(gsl.nofd.rcp85.gam[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(fastdyn.rcp26.uni.05[t], fastdyn.rcp26.uni.50[t], fastdyn.rcp26.uni.95[t]) := quantile(gsl.rcp26.uni[t,]-gsl.nofd.rcp26.uni[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(fastdyn.rcp45.uni.05[t], fastdyn.rcp45.uni.50[t], fastdyn.rcp45.uni.95[t]) := quantile(gsl.rcp45.uni[t,]-gsl.nofd.rcp45.uni[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(fastdyn.rcp85.uni.05[t], fastdyn.rcp85.uni.50[t], fastdyn.rcp85.uni.95[t]) := quantile(gsl.rcp85.uni[t,]-gsl.nofd.rcp85.uni[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(fastdyn.rcp26.gam.05[t], fastdyn.rcp26.gam.50[t], fastdyn.rcp26.gam.95[t]) := quantile(gsl.rcp26.gam[t,]-gsl.nofd.rcp26.gam[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(fastdyn.rcp45.gam.05[t], fastdyn.rcp45.gam.50[t], fastdyn.rcp45.gam.95[t]) := quantile(gsl.rcp45.gam[t,]-gsl.nofd.rcp45.gam[t,], c(0.05,.50,.95), na.rm=TRUE)
    c(fastdyn.rcp85.gam.05[t], fastdyn.rcp85.gam.50[t], fastdyn.rcp85.gam.95[t]) := quantile(gsl.rcp85.gam[t,]-gsl.nofd.rcp85.gam[t,], c(0.05,.50,.95), na.rm=TRUE)
}

iproj = which(t.proj==2000):which(t.proj==2100)
ifirst.rcp85 = rep(0,ncol(gsl.rcp26.gam))
for (i in 1:length(ifirst.rcp85)) {ifirst.rcp85[i]=which(gsl.rcp85.gam[,i]-gsl.nofd.rcp85.gam[,i]>0)[1]}
ifirst.rcp45 = rep(0,ncol(gsl.rcp26.gam))
for (i in 1:length(ifirst.rcp45)) {ifirst.rcp45[i]=which(gsl.rcp45.gam[,i]-gsl.nofd.rcp45.gam[,i]>0)[1]}
disint.rcp45=quantile(t.proj[ifirst.rcp45],c(.05,.50,.95),na.rm=TRUE)
disint.rcp85=quantile(t.proj[ifirst.rcp85],c(.05,.50,.95),na.rm=TRUE)

print('==============================================================')
print('min/5%/50%/95%/max of 2100 sea level relative to 1986-2005:')
print(paste('RCP2.6: ',quantile(gsl.rcp26.gam[251,],c(0,.05,.50,.95,1))))
print(paste('RCP4.5: ',quantile(gsl.rcp45.gam[251,],c(0,.05,.50,.95,1))))
print(paste('RCP8.5: ',quantile(gsl.rcp85.gam[251,],c(0,.05,.50,.95,1))))
print('==============================================================')

print('==============================================================')
print('min/5%/50%/95%/max of 2100 fast dynamics contribution:')
print(paste('RCP2.6: ',quantile(gsl.rcp26.gam[251,]-gsl.nofd.rcp26.gam[251,],c(0,.05,.50,.95,1))))
print(paste('RCP4.5: ',quantile(gsl.rcp45.gam[251,]-gsl.nofd.rcp45.gam[251,],c(0,.05,.50,.95,1))))
print(paste('RCP8.5: ',quantile(gsl.rcp85.gam[251,]-gsl.nofd.rcp85.gam[251,],c(0,.05,.50,.95,1))))
print('==============================================================')

print('==============================================================')
print('5%/50%/95% of fast dynamics onset:')
print(paste('RCP4.5: ',disint.rcp45))
print(paste('RCP8.5: ',disint.rcp85))
print('==============================================================')


pdf(paste(plotdir,'projections_SLR_total.pdf',sep=''),width=4,height=5.5,colormodel='cmyk')
par(mfrow=c(2,1))
# GAMMA RCP85
par(mai=c(.3,.83,.20,.2))
plot(t.proj[iproj],gsl.rcp85.gam.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2100), ylim=c(0,2), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2100,by=20)); axis(2, seq(0,2,by=.25), lab=c('0','','0.5','','1','','1.5','','2'));
		 mtext(side=2, text='Total sea level [m]', line=2.2, cex=1);
     mtext(side=3, text=expression(bold(' a')), line=-1, cex=1, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp85.gam.95[iproj],rev(gsl.rcp85.gam.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + UNIFORM RCP45
	lines(t.proj[iproj],gsl.rcp45.gam.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp45.gam.95[iproj],rev(gsl.rcp45.gam.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
# + UNIFORM RCP26
	lines(t.proj[iproj],gsl.rcp26.gam.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp26.gam.95[iproj],rev(gsl.rcp26.gam.05[iproj])),
          col=rgb(col26[1],col26[2],col26[3],.5), border=NA);
# + legend
  legend(t.proj[iproj[1]]+10,2,c("5-95% range,",
                                "RCP2.6",
																"RCP4.5",
																"RCP8.5"),
         lty=c(NA,1,1,1), lwd=3, col=c(NA,rgb(col85[1],col85[2],col85[3]),rgb(col45[1],col45[2],col45[3]),rgb(col26[1],col26[2],col26[3])), bty='n', cex=1)
# GAMMA RCP85, FD contribution
par(mai=c(.65,.83,.2,.2))
plot(t.proj[iproj],fastdyn.rcp85.gam.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2100), ylim=c(0,0.75), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2100,by=20)); axis(2, seq(0,2,by=.2)); #axis(3, seq(2000,2100,by=20)); axis(4, seq(0,2,by=.2));
		 mtext(side=1, text='Year', line=2.1, cex=1);
		 mtext(side=2, text='Fast dynamics sea \nlevel contribution [m]', line=2.2, cex=1);
     mtext(side=3, text=expression(bold(' b')), line=-1, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp85.gam.95[iproj],rev(fastdyn.rcp85.gam.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + UNIFORM RCP45, FD contribution
	lines(t.proj[iproj],fastdyn.rcp45.gam.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp45.gam.95[iproj],rev(fastdyn.rcp45.gam.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);

dev.off()

##==============================================================================
##==============================================================================









##==============================================================================
##==============================================================================
## EXTENDED DATA FIGURE 3 -- FAST DYNAMICS PROJECTIONS TO 2200
##=======================

iproj = which(t.proj==2000):which(t.proj==2200)
iend=iproj[length(iproj)]

print('==============================================================')
print('min/5%/50%/95%/max of 2100 sea level relative to 1986-2005:')
print(paste('RCP2.6: ',quantile(gsl.rcp26.gam[iend,],c(0,.05,.50,.95,1))))
print(paste('RCP4.5: ',quantile(gsl.rcp45.gam[iend,],c(0,.05,.50,.95,1))))
print(paste('RCP8.5: ',quantile(gsl.rcp85.gam[iend,],c(0,.05,.50,.95,1))))
print('==============================================================')

print('==============================================================')
print('min/5%/50%/95%/max of 2100 fast dynamics contribution:')
print(paste('RCP2.6: ',quantile(gsl.rcp26.gam[iend,]-gsl.nofd.rcp26.gam[iend,],c(0,.05,.50,.95,1))))
print(paste('RCP4.5: ',quantile(gsl.rcp45.gam[iend,]-gsl.nofd.rcp45.gam[iend,],c(0,.05,.50,.95,1))))
print(paste('RCP8.5: ',quantile(gsl.rcp85.gam[iend,]-gsl.nofd.rcp85.gam[iend,],c(0,.05,.50,.95,1))))
print('==============================================================')


pdf(paste(plotdir,'projections_SLR_total_2200.pdf',sep=''),width=4,height=5.5,colormodel='cmyk')
par(mfrow=c(2,1))
# GAMMA RCP85
par(mai=c(.3,.83,.20,.2))
plot(t.proj[iproj],gsl.rcp85.gam.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2200), ylim=c(0,8.5), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 		 axis(1, seq(2000,2200,by=25),label=c('2000','','2050','','2100','','2150','','2200')); axis(2, seq(0,8,by=1), lab=c('0','','2','','4','','6','','8'));
		 mtext(side=2, text='Total sea level [m]', line=2.2, cex=1);
     mtext(side=3, text=expression(bold(' a')), line=-1, cex=1, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp85.gam.95[iproj],rev(gsl.rcp85.gam.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + UNIFORM RCP45
	lines(t.proj[iproj],gsl.rcp45.gam.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp45.gam.95[iproj],rev(gsl.rcp45.gam.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
# + UNIFORM RCP26
	lines(t.proj[iproj],gsl.rcp26.gam.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp26.gam.95[iproj],rev(gsl.rcp26.gam.05[iproj])),
          col=rgb(col26[1],col26[2],col26[3],.5), border=NA);
# + legend
  legend(t.proj[iproj[1]]+10,8,c("5-95% range,",
                                "RCP2.6",
                                "RCP4.5",
								"RCP8.5"),
         lty=c(NA,1,1,1), lwd=3, col=c(NA,rgb(col85[1],col85[2],col85[3]),rgb(col45[1],col45[2],col45[3]),rgb(col26[1],col26[2],col26[3])), bty='n', cex=1)
# GAMMA RCP85, FD contribution
par(mai=c(.65,.83,.2,.2))
plot(t.proj[iproj],fastdyn.rcp85.gam.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2200), ylim=c(0,2.5), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2200,by=25),label=c('2000','','2050','','2100','','2150','','2200')); axis(2, seq(0,2.5,by=0.5), lab=c('0','','1','','2',''));
		 mtext(side=1, text='Year', line=2.1, cex=1);
		 mtext(side=2, text='Fast dynamics sea \nlevel contribution [m]', line=2.2, cex=1);
     mtext(side=3, text=expression(bold(' b')), line=-1, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp85.gam.95[iproj],rev(fastdyn.rcp85.gam.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + UNIFORM RCP45, FD contribution
	lines(t.proj[iproj],fastdyn.rcp45.gam.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp45.gam.95[iproj],rev(fastdyn.rcp45.gam.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);

dev.off()

##
## OR DO WE WANT THE ONE WITH BOTH PRIORS?
##


pdf(paste(plotdir,'projections_SLR_total_2200_bothPriors.pdf',sep=''),width=6.2,height=5.5,colormodel='cmyk')
par(mfrow=c(2,2))
# UNIFORM RCP85
par(mai=c(.3,.75,.25,0))
plot(t.proj[iproj],gsl.rcp85.uni.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2200), ylim=c(0,8.7), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 		 axis(1, seq(2000,2200,by=25),label=c('2000','','2050','','2100','','2150','','2200')); axis(2, seq(0,8,by=1), lab=c('0','','2','','4','','6','','8'));
		 mtext(side=2, text='Total sea level [m]', line=2.2, cex=.9);
     mtext(side=3, text='Uniform priors', line=0.3, cex=.9);
     mtext(side=3, text=expression(bold(' a')), line=-1, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp85.uni.95[iproj],rev(gsl.rcp85.uni.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + UNIFORM RCP26
	lines(t.proj[iproj],gsl.rcp26.uni.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp26.uni.95[iproj],rev(gsl.rcp26.uni.05[iproj])),
          col=rgb(col26[1],col26[2],col26[3],.5), border=NA);
# + UNIFORM RCP45
	lines(t.proj[iproj],gsl.rcp45.uni.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp45.uni.95[iproj],rev(gsl.rcp45.uni.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
# + legend
  legend(t.proj[iproj[1]]+10,8,c("5-95% range,",
                                "RCP2.6",
																"RCP4.5",
																"RCP8.5"),
         lty=c(NA,1,1,1), lwd=3, col=c(NA,rgb(col26[1],col26[2],col26[3]),rgb(col45[1],col45[2],col45[3]),rgb(col85[1],col85[2],col85[3])), bty='n', cex=.9)
# GAMMA RCP85
par(mai=c(.3,.6,.25,.2))
plot(t.proj[iproj],gsl.rcp85.gam.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2200), ylim=c(0,8.7), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 		 axis(1, seq(2000,2200,by=25),label=c('2000','','2050','','2100','','2150','','2200')); axis(2, seq(0,8,by=1), lab=c('0','','2','','4','','6','','8'));
     mtext(side=3, text='Gamma priors', line=0.3, cex=.9);
     mtext(side=3, text=expression(bold(' b')), line=-1, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp85.gam.95[iproj],rev(gsl.rcp85.gam.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + GAMMA RCP26
	lines(t.proj[iproj],gsl.rcp26.gam.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp26.gam.95[iproj],rev(gsl.rcp26.gam.05[iproj])),
          col=rgb(col26[1],col26[2],col26[3],.5), border=NA);
# + GAMMA RCP45
	lines(t.proj[iproj],gsl.rcp45.gam.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp45.gam.95[iproj],rev(gsl.rcp45.gam.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
# UNIFORM RCP85, FD contribution
par(mai=c(.6,.75,.1,0))
plot(t.proj[iproj],fastdyn.rcp85.uni.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2200), ylim=c(0,2.5), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2200,by=25),label=c('2000','','2050','','2100','','2150','','2200')); axis(2, seq(0,2.5,by=0.5), lab=c('0','','1','','2',''));
		 mtext(side=1, text='Year', line=2.1, cex=.9);
		 mtext(side=2, text='Fast dynamics sea \nlevel contribution [m]', line=2.2, cex=.9);
     mtext(side=3, text=expression(bold(' c')), line=-1, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp85.uni.95[iproj],rev(fastdyn.rcp85.uni.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + UNIFORM RCP45, FD contribution
	lines(t.proj[iproj],fastdyn.rcp45.uni.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp45.uni.95[iproj],rev(fastdyn.rcp45.uni.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
# GAMMA RCP85, FD contribution
par(mai=c(.6,.6,.1,.2))
plot(t.proj[iproj],fastdyn.rcp85.gam.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2200), ylim=c(0,2.5), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2200,by=25),label=c('2000','','2050','','2100','','2150','','2200')); axis(2, seq(0,2.5,by=0.5), lab=c('0','','1','','2',''));
		 mtext(side=1, text='Year', line=2.1, cex=.9);
     mtext(side=3, text=expression(bold(' d')), line=-1, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp85.gam.95[iproj],rev(fastdyn.rcp85.gam.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + GAMMA RCP45, FD contribution
	lines(t.proj[iproj],fastdyn.rcp45.gam.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp45.gam.95[iproj],rev(fastdyn.rcp45.gam.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);

dev.off()

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## FIGURE 3 -- PROJECTIONS, PROBABILITY DISTRIBUTIONS, SURVIVAL FUNCTIONS
##=========

i2100 = which(t.proj==2100)

slr2100.rcp26.uni = gsl.rcp26.uni[i2100,]
slr2100.rcp45.uni = gsl.rcp45.uni[i2100,]
slr2100.rcp85.uni = gsl.rcp85.uni[i2100,]
slr2100.nofd.rcp26.uni = gsl.nofd.rcp26.uni[i2100,]
slr2100.nofd.rcp45.uni = gsl.nofd.rcp45.uni[i2100,]
slr2100.nofd.rcp85.uni = gsl.nofd.rcp85.uni[i2100,]

slr2100.rcp26.gam = gsl.rcp26.gam[i2100,]
slr2100.rcp45.gam = gsl.rcp45.gam[i2100,]
slr2100.rcp85.gam = gsl.rcp85.gam[i2100,]
slr2100.nofd.rcp26.gam = gsl.nofd.rcp26.gam[i2100,]
slr2100.nofd.rcp45.gam = gsl.nofd.rcp45.gam[i2100,]
slr2100.nofd.rcp85.gam = gsl.nofd.rcp85.gam[i2100,]

pdf.slr2100.rcp26.uni = density(slr2100.rcp26.uni,from=0,to=4,na.rm=TRUE)
pdf.slr2100.rcp45.uni = density(slr2100.rcp45.uni,from=0,to=4,na.rm=TRUE)
pdf.slr2100.rcp85.uni = density(slr2100.rcp85.uni,from=0,to=4,na.rm=TRUE)
pdf.slr2100.nofd.rcp26.uni = density(slr2100.nofd.rcp26.uni,from=0,to=4,na.rm=TRUE)
pdf.slr2100.nofd.rcp45.uni = density(slr2100.nofd.rcp45.uni,from=0,to=4,na.rm=TRUE)
pdf.slr2100.nofd.rcp85.uni = density(slr2100.nofd.rcp85.uni,from=0,to=4,na.rm=TRUE)

pdf.slr2100.rcp26.gam = density(slr2100.rcp26.gam,from=0,to=4,na.rm=TRUE)
pdf.slr2100.rcp45.gam = density(slr2100.rcp45.gam,from=0,to=4,na.rm=TRUE)
pdf.slr2100.rcp85.gam = density(slr2100.rcp85.gam,from=0,to=4,na.rm=TRUE)
pdf.slr2100.nofd.rcp26.gam = density(slr2100.nofd.rcp26.gam,from=0,to=4,na.rm=TRUE)
pdf.slr2100.nofd.rcp45.gam = density(slr2100.nofd.rcp45.gam,from=0,to=4,na.rm=TRUE)
pdf.slr2100.nofd.rcp85.gam = density(slr2100.nofd.rcp85.gam,from=0,to=4,na.rm=TRUE)

cdf.slr2100.rcp26.uni = rep(0,length(pdf.slr2100.rcp26.uni$x))
cdf.slr2100.rcp45.uni = rep(0,length(pdf.slr2100.rcp45.uni$x))
cdf.slr2100.rcp85.uni = rep(0,length(pdf.slr2100.rcp85.uni$x))
cdf.slr2100.nofd.rcp26.uni = rep(0,length(pdf.slr2100.nofd.rcp26.uni$x))
cdf.slr2100.nofd.rcp45.uni = rep(0,length(pdf.slr2100.nofd.rcp45.uni$x))
cdf.slr2100.nofd.rcp85.uni = rep(0,length(pdf.slr2100.nofd.rcp85.uni$x))

cdf.slr2100.rcp26.gam = rep(0,length(pdf.slr2100.rcp26.gam$x))
cdf.slr2100.rcp45.gam = rep(0,length(pdf.slr2100.rcp45.gam$x))
cdf.slr2100.rcp85.gam = rep(0,length(pdf.slr2100.rcp85.gam$x))
cdf.slr2100.nofd.rcp26.gam = rep(0,length(pdf.slr2100.nofd.rcp26.gam$x))
cdf.slr2100.nofd.rcp45.gam = rep(0,length(pdf.slr2100.nofd.rcp45.gam$x))
cdf.slr2100.nofd.rcp85.gam = rep(0,length(pdf.slr2100.nofd.rcp85.gam$x))

dx=median(diff(pdf.slr2100.rcp26.uni$x))
x=pdf.slr2100.rcp26.uni$x

for (i in 2:length(x)){
  cdf.slr2100.rcp26.uni[i] = sum(c(cdf.slr2100.rcp26.uni[i-1],dx*pdf.slr2100.rcp26.uni$y[i]))
  cdf.slr2100.rcp45.uni[i] = sum(c(cdf.slr2100.rcp45.uni[i-1],dx*pdf.slr2100.rcp45.uni$y[i]))
  cdf.slr2100.rcp85.uni[i] = sum(c(cdf.slr2100.rcp85.uni[i-1],dx*pdf.slr2100.rcp85.uni$y[i]))
  cdf.slr2100.nofd.rcp26.uni[i] = sum(c(cdf.slr2100.nofd.rcp26.uni[i-1],dx*pdf.slr2100.nofd.rcp26.uni$y[i]))
  cdf.slr2100.nofd.rcp45.uni[i] = sum(c(cdf.slr2100.nofd.rcp45.uni[i-1],dx*pdf.slr2100.nofd.rcp45.uni$y[i]))
  cdf.slr2100.nofd.rcp85.uni[i] = sum(c(cdf.slr2100.nofd.rcp85.uni[i-1],dx*pdf.slr2100.nofd.rcp85.uni$y[i]))
  cdf.slr2100.rcp26.gam[i] = sum(c(cdf.slr2100.rcp26.gam[i-1],dx*pdf.slr2100.rcp26.gam$y[i]))
  cdf.slr2100.rcp45.gam[i] = sum(c(cdf.slr2100.rcp45.gam[i-1],dx*pdf.slr2100.rcp45.gam$y[i]))
  cdf.slr2100.rcp85.gam[i] = sum(c(cdf.slr2100.rcp85.gam[i-1],dx*pdf.slr2100.rcp85.gam$y[i]))
  cdf.slr2100.nofd.rcp26.gam[i] = sum(c(cdf.slr2100.nofd.rcp26.gam[i-1],dx*pdf.slr2100.nofd.rcp26.gam$y[i]))
  cdf.slr2100.nofd.rcp45.gam[i] = sum(c(cdf.slr2100.nofd.rcp45.gam[i-1],dx*pdf.slr2100.nofd.rcp45.gam$y[i]))
  cdf.slr2100.nofd.rcp85.gam[i] = sum(c(cdf.slr2100.nofd.rcp85.gam[i-1],dx*pdf.slr2100.nofd.rcp85.gam$y[i]))
}

cdf.slr2100.rcp26.uni = cdf.slr2100.rcp26.uni/cdf.slr2100.rcp26.uni[length(cdf.slr2100.rcp26.uni)]
cdf.slr2100.rcp45.uni = cdf.slr2100.rcp45.uni/cdf.slr2100.rcp45.uni[length(cdf.slr2100.rcp45.uni)]
cdf.slr2100.rcp85.uni = cdf.slr2100.rcp85.uni/cdf.slr2100.rcp85.uni[length(cdf.slr2100.rcp85.uni)]
cdf.slr2100.nofd.rcp26.uni = cdf.slr2100.nofd.rcp26.uni/cdf.slr2100.nofd.rcp26.uni[length(cdf.slr2100.nofd.rcp26.uni)]
cdf.slr2100.nofd.rcp45.uni = cdf.slr2100.nofd.rcp45.uni/cdf.slr2100.nofd.rcp45.uni[length(cdf.slr2100.nofd.rcp45.uni)]
cdf.slr2100.nofd.rcp85.uni = cdf.slr2100.nofd.rcp85.uni/cdf.slr2100.nofd.rcp85.uni[length(cdf.slr2100.nofd.rcp85.uni)]

cdf.slr2100.rcp26.gam = cdf.slr2100.rcp26.gam/cdf.slr2100.rcp26.gam[length(cdf.slr2100.rcp26.gam)]
cdf.slr2100.rcp45.gam = cdf.slr2100.rcp45.gam/cdf.slr2100.rcp45.gam[length(cdf.slr2100.rcp45.gam)]
cdf.slr2100.rcp85.gam = cdf.slr2100.rcp85.gam/cdf.slr2100.rcp85.gam[length(cdf.slr2100.rcp85.gam)]
cdf.slr2100.nofd.rcp26.gam = cdf.slr2100.nofd.rcp26.gam/cdf.slr2100.nofd.rcp26.gam[length(cdf.slr2100.nofd.rcp26.gam)]
cdf.slr2100.nofd.rcp45.gam = cdf.slr2100.nofd.rcp45.gam/cdf.slr2100.nofd.rcp45.gam[length(cdf.slr2100.nofd.rcp45.gam)]
cdf.slr2100.nofd.rcp85.gam = cdf.slr2100.nofd.rcp85.gam/cdf.slr2100.nofd.rcp85.gam[length(cdf.slr2100.nofd.rcp85.gam)]

sur.slr2100.rcp26.uni = 1-cdf.slr2100.rcp26.uni
sur.slr2100.rcp45.uni = 1-cdf.slr2100.rcp45.uni
sur.slr2100.rcp85.uni = 1-cdf.slr2100.rcp85.uni
sur.slr2100.nofd.rcp26.uni = 1-cdf.slr2100.nofd.rcp26.uni
sur.slr2100.nofd.rcp45.uni = 1-cdf.slr2100.nofd.rcp45.uni
sur.slr2100.nofd.rcp85.uni = 1-cdf.slr2100.nofd.rcp85.uni

sur.slr2100.rcp26.gam = 1-cdf.slr2100.rcp26.gam
sur.slr2100.rcp45.gam = 1-cdf.slr2100.rcp45.gam
sur.slr2100.rcp85.gam = 1-cdf.slr2100.rcp85.gam
sur.slr2100.nofd.rcp26.gam = 1-cdf.slr2100.nofd.rcp26.gam
sur.slr2100.nofd.rcp45.gam = 1-cdf.slr2100.nofd.rcp45.gam
sur.slr2100.nofd.rcp85.gam = 1-cdf.slr2100.nofd.rcp85.gam


pdf(paste(plotdir,'distributions_SLR2100_pdf+sf.pdf',sep=''),width=5,height=7,colormodel='cmyk')
par(mfrow=c(2,1), mai=c(.85,.74,.1,.15))

plot(x,pdf.slr2100.rcp85.uni$y, type='l', xlim=c(0,3), ylim=c(0,4.7), lty=1,
     col=rgb(col85[1],col85[2],col85[3]), lwd=1.5, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i',axes=FALSE);
  axis(1,seq(0,3,0.5),lab=c("0","0.5","1","1.5","2","2.5","3"))
  u <- par("usr")
  arrows(0, u[3],0, u[4], code = 2, xpd = TRUE)
  mtext('Probability density', side=2, line=1.3);
  mtext('Projected sea level in 2100\nrelative to 1986-2005 average [m]', side=1, line=3);
  mtext(side=3, text=expression(bold('   a')), line=-1, cex=.9, adj=0);

  lines(x,pdf.slr2100.rcp45.uni$y, type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=1.5);
  lines(x,pdf.slr2100.rcp26.uni$y, type='l', col=rgb(col26[1],col26[2],col26[3]), lwd=1.5);

  lines(x,pdf.slr2100.nofd.rcp85.uni$y, type='l', col=rgb(col85[1],col85[2],col85[3]), lwd=1.5, lty=2);
  lines(x,pdf.slr2100.nofd.rcp45.uni$y, type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=1.5, lty=2);
  lines(x,pdf.slr2100.nofd.rcp26.uni$y, type='l', col=rgb(col26[1],col26[2],col26[3]), lwd=1.5, lty=2);

  legend(1.1,4.7,c("RCP2.6","RCP4.5","RCP8.5","including fast dynamics","neglecting fast dynamics"),
        lty=c(1,1,1,1,2), lwd=2, col=c(rgb(col26[1],col26[2],col26[3]),rgb(col45[1],col45[2],col45[3]),rgb(col85[1],col85[2],col85[3]),'black','black'),
        bty='n')

plot(x,log10(sur.slr2100.rcp85.uni),type='l', xlim=c(0,3), ylim=c(-3.3,0), lty=1,
     col=rgb(col85[1],col85[2],col85[3]), lwd=1.5, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i');
  axis(1,seq(0,3,0.5),lab=c("0","0.5","1","1.5","2","2.5","3"))
  mtext('Survival function [1-CDF]', side = 2, line=2.6);
  mtext('Projected sea level in 2100\nrelative to 1986-2005 average [m]', side=1, line=3);
  mtext(side=3, text=expression(bold('   b')), line=-1.1, cex=.9, adj=0);

	axis(2, at=seq(-4,0), label=parse(text=paste("10^", seq(-4,0), sep="")), las=1)

  lines(x,log10(sur.slr2100.nofd.rcp85.uni),type='l',col=rgb(col85[1],col85[2],col85[3]), lty=2, lwd=1.5);
  lines(x,log10(sur.slr2100.rcp45.uni),type='l',lty=1, col=rgb(col45[1],col45[2],col45[3]), lwd=1.5);    lines(x,log10(sur.slr2100.nofd.rcp45.uni),type='l',col=rgb(col45[1],col45[2],col45[3]), lty=2, lwd=1.5);
  lines(x,log10(sur.slr2100.rcp26.uni),type='l',lty=1, col=rgb(col26[1],col26[2],col26[3]), lwd=1.5); lines(x,log10(sur.slr2100.nofd.rcp26.uni),type='l',col=rgb(col26[1],col26[2],col26[3]), lty=2, lwd=1.5);

	lines(c(-4,4),c(-2,-2),lty=2,col='black'); text(0.35,-1.85,"1:100 level");
	lines(c(-4,4),c(-3,-3),lty=2,col='black'); text(0.35,-2.85,"1:1000 level");

dev.off()

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## FIGURE 4 -- VAN DANTZIG
##=========

## Grab the van Dantzig output
## (switch to RCP26 or 45 if you want to evaluate those, but not the focus
## of the main paper)
#filename.vandantzig.gamma   = '../output_model/vanDantzig_RCP26_gamma_29Jan2017.nc'
#filename.vandantzig.gamma   = '../output_model/vanDantzig_RCP45_gamma_29Jan2017.nc'

ncdata <- nc_open(filename.vandantzig.gamma)
  heightening     = ncvar_get(ncdata, 'H')
  cost            = ncvar_get(ncdata, 'ExpectedCost')
  cost.nofd       = ncvar_get(ncdata, 'ExpectedCost_nofd')
  loss            = ncvar_get(ncdata, 'ExpectedLoss')
  loss.nofd       = ncvar_get(ncdata, 'ExpectedLoss_nofd')
  investment      = ncvar_get(ncdata, 'ExpectedInvestment')
  investment.nofd = ncvar_get(ncdata, 'ExpectedInvestment_nofd')
  preturn         = ncvar_get(ncdata, 'ExpectedPreturn')
  preturn.nofd    = ncvar_get(ncdata, 'ExpectedPreturn_nofd')
nc_close(ncdata)

n.ensemble = ncol(cost)
n.height = nrow(cost)

## What is the optimal heightening (minimum expected costs) for each ensemble
## member? (and the index, so we can grab the return period)
iopt=rep(NA,n.ensemble)			# index of optimal heightening, return period
iopt.nofd=rep(NA,n.ensemble)	# index of optimal heightening, return period
rp.opt=rep(NA,n.ensemble)		# return period of optimal heightening (including fast dynamics)
rp.nofd.opt=rep(NA,n.ensemble)	# return period of optimal heightening (without considering fast dynamics)
rp.nofd.act=rep(NA,n.ensemble)	# actual return period (when you consider fast dynamics)
Hopt=rep(NA,n.ensemble)         # optimal heightening, accounting for fast dynamics
Hopt.nofd=rep(NA,n.ensemble)    # what you *think* is the optimal heightening, neglecting fast dynamics
waste = rep(NA,n.ensemble)      # difference between expected costs if you had built optimally (accounting
                                # for fast dynamics) and if you built ignorant of the fast dynamics

for (i in 1:n.ensemble) {
    iopt[i] = which(cost[,i]==min(cost[,i]))
	iopt.nofd[i] = which(cost.nofd[,i]==min(cost.nofd[,i]))
    Hopt[i] = heightening[iopt[i]]
	Hopt.nofd[i] = heightening[iopt.nofd[i]]
    rp.nofd.opt[i] = preturn.nofd[iopt.nofd[i],i]
	rp.nofd.act[i] = preturn[iopt.nofd[i],i]
    rp.opt[i] = preturn[iopt[i],i]
    waste[i] = cost[iopt.nofd[i],i] - cost[iopt[i],i]
}

quan.nofd.opt = c(quantile(rp.nofd.opt, 0.05) , quantile(rp.nofd.opt, 0.50) , quantile(rp.nofd.opt, 0.95) )
quan.nofd.act = c(quantile(rp.nofd.act, 0.05) , quantile(rp.nofd.act, 0.50) , quantile(rp.nofd.act, 0.95) )
quan.opt      = c(quantile(rp.opt, 0.05)      , quantile(rp.opt, 0.50)      , quantile(rp.opt, 0.95)      )

preturn.avg = apply(preturn, 1, mean)
cost.avg    = apply(cost, 1, mean)
loss.avg    = apply(loss, 1, mean)
investment.avg    = apply(investment, 1, mean)
preturn.nofd.avg  = apply(preturn.nofd, 1, mean)
cost.nofd.avg     = apply(cost.nofd, 1, mean)
loss.nofd.avg     = apply(loss.nofd, 1, mean)
investment.nofd.avg = apply(investment.nofd, 1, mean)
iopt.avg = which(cost.avg==min(cost.avg))
iopt.nofd.avg = which(cost.nofd.avg==min(cost.nofd.avg))

overspending.avg = signif((cost.avg[iopt.nofd.avg]-cost.nofd.avg[iopt.nofd.avg])/1e6,2)
waste.avg = signif((cost.avg[iopt.nofd.avg]-cost.avg[iopt.avg])/1e6,2)

#preturn.pre = 10000                                  # prescribed return period? (years)
preturn.pre = round(preturn.nofd.avg[iopt.nofd.avg])  # use optimal return period? (assuming no-FD)

ipre = rep(NA,n.ensemble)
ipre.nofd = rep(NA,n.ensemble)
for (i in 1:n.ensemble) {
  ipre[i] = which(preturn[,i]>=preturn.pre)[1]
  ipre.nofd[i] = which(preturn.nofd[,i]>=preturn.pre)[1]
}
ipre.avg = which(preturn.avg >= preturn.pre)[1]
ipre.nofd.avg = which(preturn.nofd.avg >= preturn.pre)[1]

tmp=lm( investment.avg ~ log10(preturn.avg))
inter = tmp$coefficients[[1]]
slope = tmp$coefficients[[2]]
tmp=lm( investment.nofd.avg ~ log10(preturn.nofd.avg))
inter.nofd = tmp$coefficients[[1]]
slope.nofd = tmp$coefficients[[2]]

## Calculate ensemble median and 5-95% range for costs as a function of
## (a) heightening, and (b) return period
cost.h = mat.or.vec(length(heightening),3)
cost.nofd.h = mat.or.vec(length(heightening),3)

for (t in 1:length(heightening)) {
  cost.h[t,1] = quantile(cost[t,],.05); cost.nofd.h[t,1] = quantile(cost.nofd[t,],.05)
  cost.h[t,2] = quantile(cost[t,],.50); cost.nofd.h[t,2] = quantile(cost.nofd[t,],.50)
  cost.h[t,3] = quantile(cost[t,],.95); cost.nofd.h[t,3] = quantile(cost.nofd[t,],.95)
}

## Return period is trickier. Need to bin, because not linear.
rp.bin = seq(from=log10(min(c(preturn,preturn.nofd))), to=log10(max(c(preturn,preturn.nofd))), length.out=100)
cost.rp = mat.or.vec(length(rp.bin),3)
cost.nofd.rp = mat.or.vec(length(rp.bin),3)
for (t in 1:length(rp.bin)) {
  itmp=which(log10(preturn)>rp.bin[t] & log10(preturn)<rp.bin[t+1])
  cost.rp[t,1] = quantile(cost[itmp],0.05, na.rm=TRUE); cost.rp[t,2] = quantile(cost[itmp],0.50, na.rm=TRUE); cost.rp[t,3] = quantile(cost[itmp],0.95, na.rm=TRUE);
  itmp=which(log10(preturn.nofd)>rp.bin[t] & log10(preturn.nofd)<rp.bin[t+1])
  cost.nofd.rp[t,1] = quantile(cost.nofd[itmp],0.05, na.rm=TRUE); cost.nofd.rp[t,2] = quantile(cost.nofd[itmp],0.50, na.rm=TRUE); cost.nofd.rp[t,3] = quantile(cost.nofd[itmp],0.95, na.rm=TRUE);
}
itrim=which(!is.na(cost.nofd.rp[,1]) & !is.na(cost.nofd.rp[,2]) & !is.na(cost.nofd.rp[,3]) &
            !is.na(cost.rp[,1]) & !is.na(cost.rp[,2]) & !is.na(cost.rp[,3]) )
rp.bin=rp.bin[itrim]
cost.rp = cost.rp[itrim,]
cost.nofd.rp = cost.nofd.rp[itrim,]

#incheight.avg = signif((heightening[iopt.avg]-heightening[iopt.nofd.avg])/.3048,2)#feet
incheight.avg = signif((heightening[iopt.avg]-heightening[iopt.nofd.avg]),2)      #meters

print('====================================================================')
print(paste('no-FD rec. heightening =',heightening[iopt.nofd.avg]/.3048,'ft'))
print(paste('w/-FD rec. heightening =',heightening[iopt.avg]/.3048,'ft'))
print('====================================================================')
print(paste('no-FD rec. heightening =',heightening[iopt.nofd.avg],'m'))
print(paste('w/-FD rec. heightening =',heightening[iopt.avg],'m'))
print('====================================================================')


# this version does not have the shaded uncertainty regions on it

# BOTH RETURN PERIOD AND HEIGHTENING, SAME COST AXIS
## Note: the inset lines might require some hard-coded finessing.

## RETURN PERIOD

pdf(paste(plotdir,'vandantzig_RP+H_noUnc.pdf',sep=''),width=7.25,height=4,colormodel='cmyk')

inset.x = c(2.63,3.4)
inset.y = c(1.85e9,2.32e9)
conv=1e9  # convert from $ to billions or millions of $? (for nicer looking axes)
pfig=0.55

par(mfrow=c(1,2), fig=c(0,pfig,0,1), mai=c(.7,.7,.06,.2))
plot(log10(preturn.avg),cost.avg/conv, col=rgb(col85[1],col85[2],col85[3]),
     type='l',xlim=c(1.5,7),ylim=c(0,8e9)/conv,xaxt='n', xlab='', ylab='');
mtext('Return period [years]',side = 1,line = 2.1)
mtext('Expected costs [billion US $]',side=2,line=2.2);
axis(1, at=seq(1,7), label=parse(text=paste("10^", seq(1,7), sep="")), las=1)
lines(log10(preturn.avg)      , cost.avg/conv      , type='l', col=rgb(col85[1],col85[2],col85[3]), lwd=2);
lines(log10(preturn.nofd.avg) , cost.nofd.avg/conv , type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2);
lines(log10(preturn.avg)      , investment.avg/conv      , type='l', col=rgb(col85[1],col85[2],col85[3]), lwd=2, lty=6);
lines(log10(preturn.nofd.avg) , investment.nofd.avg/conv , type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, lty=6);
mtext(side=3, text=expression(bold(' a')), line=-2, cex=1, adj=0);
points(log10(preturn.pre),cost.nofd.avg[iopt.nofd.avg]/conv,type='p',pch=16,cex=1.4)

polygon(c(inset.x,rev(inset.x)),c(inset.y[1],inset.y[1],inset.y[2],inset.y[2])/conv,col=rgb(0,0,0,0),border=TRUE,lwd=1.5)
legend(4.15,1.2,c("Total cost","Investment"), col=c('black','black'), lty=c(1,6), lwd=2, bty='n', cex=1.1)

u <- par("usr")
v <- c(
    grconvertX(u[1:2], "user", "ndc"),
    grconvertY(u[3:4], "user", "ndc")
)
py=.58
px=.28
v <- c( (v[1]+v[2])*px, v[2], (v[4]+v[3])*py, v[4] )

lines( c(inset.x[2],u[2]*.995), c(inset.y[1]/conv,v[3]*u[4]+2*u[3]), lty=1, lwd=1);
lines( c(inset.x[1]-.05,v[1]*(u[2]-u[1])+u[1])+.05, c(inset.y[1]/conv,v[3]*u[4]+2*u[3]), lty=1, lwd=1);
lines( c(log10(preturn.pre),log10(preturn.pre)), c(-9,investment.avg[ipre.nofd.avg]/conv), lty=2, lwd=2);
lines( c(log10(preturn.avg[ipre.nofd.avg]), log10(preturn.avg[ipre.nofd.avg])),
       c(investment.avg[ipre.nofd.avg]/conv, -9), lty=2, lwd=2);

rect(u[2], u[4], (u[1]+u[2])*px*.99, (u[4]-u[3])*py, col="white")
par( fig=v, new=TRUE, mar=c(0,0,0,0) )

plot(log10(preturn.avg),cost.avg/conv, col=rgb(col85[1],col85[2],col85[3]),
     type='l',xlim=inset.x,ylim=inset.y/conv, xlab='', ylab='', xaxt='n', yaxt='n');
axis(1, at=seq(3.8,4.0,by=.2), label=parse(text=paste("10^", seq(3.8,4.0,by=.2), sep="")), las=1)
lines(log10(preturn.avg)      , cost.avg/conv      , type='l', col=rgb(col85[1],col85[2],col85[3]), lwd=2);
lines(log10(preturn.nofd.avg) , cost.nofd.avg/conv , type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2);
lines(log10(preturn.avg)      , investment.avg/conv      , type='l', col=rgb(col85[1],col85[2],col85[3]), lwd=2, lty=6);
lines(log10(preturn.nofd.avg) , investment.nofd.avg/conv , type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, lty=6);

#lines( c(log10(preturn.pre),log10(preturn.pre)), c(0,(log10(preturn.pre)*slope.nofd+inter.nofd))/conv, lty=2, lwd=2);
lines( c(log10(preturn.pre),log10(preturn.pre)), c(0,cost.nofd.avg[iopt.nofd.avg])/conv, lty=2, lwd=2);

xtmp = (log10(preturn.pre)*slope.nofd+inter.nofd-inter)/slope

arrows( log10(preturn.pre), (log10(preturn.pre)*slope.nofd+inter.nofd)/conv,
        xtmp, (log10(preturn.pre)*slope.nofd+inter.nofd)/conv,
        code=2, length=.11, angle=45, lty=1, lwd=2);
lines( c(xtmp, xtmp), c((log10(preturn.pre)*slope.nofd+inter.nofd)/conv, 0), lty=2, lwd=2);
#points(log10(preturn.pre), (log10(preturn.pre)*slope.nofd+inter.nofd)/conv,type='p',pch=18,cex=2)
points(log10(preturn.pre),cost.nofd.avg[iopt.nofd.avg]/conv,type='p',pch=16,cex=1.4)

text(0.998*log10(preturn.pre), 0.985*mean(c(inset.y[1],investment.avg[ipre.nofd.avg]))/conv,
     bquote(atop("Goal:",paste("1:",.(prettyNum(preturn.pre,big.mark=',')),"y"))), pos = 4, cex = 1.1, srt = 0)
text(0.95*log10(preturn.avg[ipre.nofd.avg]), .97*investment.avg[ipre.nofd.avg]/conv,
     bquote(atop("  Actual:",paste("1:",.(prettyNum(round(preturn.avg[ipre.nofd.avg]),big.mark=',')),"y"))), pos = 3, cex = 1.1, srt = 0)
#  text(0.997*xtmp, 0.968*mean(c(inset.y[1],investment.avg[ipre.nofd.avg]))/conv,
#        paste("Actual: 1:",round(preturn.avg[ipre.nofd.avg]),"y",sep=""), pos=4, cex=1.1, srt=0)

## HEIGHTENING

inset.x = c(1.34,1.61)
inset.y = c(2.22e9,2.38e9)
conv=1e9  # convert from $ to billions or millions of $? (for nicer looking axes)

par(fig=c(pfig,1,0,1), new=TRUE, mai=c(.7,.1,.06,.06))
plot(heightening,cost.avg/conv,col=rgb(col85[1],col85[2],col85[3],.2),type='l',xlim=c(0,3),ylim=c(0,8e9)/conv,
     xlab='', ylab='', xaxt='n', yaxt='n');
mtext('Heightening [m]',side = 1,line = 2.1)
axis(1, at=seq(0,3,by=0.5), label=c('0','','1','','2','','3'))
axis(2, at=seq(0,8,by=2), label=rep("",5), las=1)
lines(heightening,cost.avg/conv,type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=3);
lines(heightening,cost.nofd.avg/conv,type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=3);
points(heightening[iopt.nofd.avg],cost.nofd.avg[iopt.nofd.avg]/conv,type='p',pch=16,cex=1.)
#  points(heightening[iopt.nofd.avg],cost.avg[iopt.nofd.avg]/conv,type='p',pch=17,cex=1.)
points(heightening[iopt.avg],cost.avg[iopt.avg]/conv,type='p',pch=15,cex=1.)
polygon(c(inset.x,rev(inset.x)),c(inset.y[1]-.15e9,inset.y[1]-.15e9,inset.y[2]+.15e9,inset.y[2]+.15e9)/conv,col=rgb(0,0,0,0),border=TRUE,lwd=1.5)
legend(.35,1.2,c("fast dynamics included", "fast dynamics neglected"),
       lty=c(1,1), lwd=2, col=c(rgb(col85[1],col85[2],col85[3]),rgb(col45[1],col45[2],col45[3])), cex=1.1, bty='n')
mtext(side=3, text=expression(bold(' b')), line=-2, cex=1, adj=0);

u <- par("usr")
v <- c(
    grconvertX(u[1:2], "user", "ndc"),
    grconvertY(u[3:4], "user", "ndc")
)
py=1.1/2
px=.9/2
v <- c( (v[1]+v[2])*px, v[2], (v[4]+v[3])*py, v[4] )

lines( c(inset.x[2],u[2]*1.01), c(inset.y[1]/conv+.5*u[3],v[3]*u[4]+2*u[3]), lty=1, lwd=1);
lines( c(inset.x[1],(u[2]-u[1])*v[1]*(1-pfig)+u[1]), c(inset.y[1]/conv+.5*u[3],v[3]*u[4]+2*u[3]), lty=1, lwd=1);

rect(u[2], u[4], (u[1]+u[2])*px-px*.99, (u[4]-u[3])*py+u[3]*(1-py), col="white")
par( fig=v, new=TRUE, mar=c(0,0,0,0) )

plot(heightening,cost.avg/conv,type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=3, xlim=inset.x, ylim=inset.y/conv,
     xlab='', ylab='', yaxt='n', xaxt='n');
lines(heightening,cost.nofd.avg/conv,type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=3);
lines(heightening,cost.avg/conv,type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=3);

lines(c(heightening[iopt.nofd.avg],heightening[iopt.nofd.avg]),c(0,cost.avg[iopt.nofd.avg])/conv,type='l',lty=2,lwd=2,col='black')
lines(c(heightening[iopt.nofd.avg],heightening[iopt.avg]),c(cost.avg[iopt.nofd.avg],cost.avg[iopt.nofd.avg])/conv,type='l',lty=2,lwd=2,col='black')
arrows(heightening[iopt.nofd.avg], cost.avg[iopt.avg]/conv, .995*heightening[iopt.avg], cost.avg[iopt.avg]/conv,
       length=.12, angle=40, lty=1, lwd=1.5, code=2, col='black')
arrows(heightening[iopt.avg], cost.avg[iopt.nofd.avg]/conv, heightening[iopt.avg], cost.avg[iopt.avg]/conv,
       length=.11, angle=40, lty=1, lwd=1.5, code=1, col='black')

points(heightening[iopt.nofd.avg],cost.nofd.avg[iopt.nofd.avg]/conv,type='p',pch=16,cex=1.4)
points(heightening[iopt.nofd.avg],cost.avg[iopt.nofd.avg]/conv,type='p',pch=17,cex=1.4)
points(heightening[iopt.avg],cost.avg[iopt.avg]/conv,type='p',pch=15,cex=1.4)

text(0.91*heightening[iopt.avg], 1.006*cost.avg[iopt.nofd.avg]/conv,
     paste("+$",waste.avg," million",sep=""), pos = 4, cex = 1.1, srt = 0)
text(heightening[iopt.nofd.avg], .992*cost.avg[iopt.avg]/conv,
     paste("+",incheight.avg," m",sep=""), pos = 4, cex = 1.1, srt = 0)

#  text(1.0*heightening[iopt.nofd.avg], .99*mean(c(cost.avg[iopt.nofd.avg],cost.nofd.avg[iopt.nofd.avg]))/conv,
#       paste("+$",overspending.avg," million",sep=""), pos = 4, cex = 1.1, srt = 0)

dev.off()

##==============================================================================
##==============================================================================





# this version has the uncertainty regions on it

## BOTH RETURN PERIOD AND HEIGHTENING, SAME COST AXIS
## Note: the inset lines might require some hard-coded finessing.

## RETURN PERIOD

pdf(paste(plotdir,'vandantzig_RP+H.pdf',sep=''),width=7.25,height=4,colormodel='cmyk')

inset.x = c(2.63,3.4)
inset.y = c(1.85e9,2.32e9)
conv=1e9  # convert from $ to billions or millions of $? (for nicer looking axes)
pfig=0.55

par(mfrow=c(1,2), fig=c(0,pfig,0,1), mai=c(.7,.7,.06,.2))
plot(log10(preturn.avg),cost.avg/conv, col=rgb(col85[1],col85[2],col85[3]),
			type='l',xlim=c(1.5,7),ylim=c(0,8e9)/conv,xaxt='n', xlab='', ylab='');
      mtext('Return period [years]',side = 1,line = 2.1)
      mtext('Expected costs [billion US $]',side=2,line=2.2);
			axis(1, at=seq(1,7), label=parse(text=paste("10^", seq(1,7), sep="")), las=1)
  polygon(c(rp.bin,rev(rp.bin)), c(cost.rp[,3],rev(cost.rp[,1]))/conv, col=rgb(col85[1],col85[2],col85[3],.3), border=NA);
  polygon(c(rp.bin,rev(rp.bin)), c(cost.nofd.rp[,3],rev(cost.nofd.rp[,1]))/conv, col=rgb(col45[1],col45[2],col45[3],.3), border=NA);
  lines(log10(preturn.avg)      , cost.avg/conv      , type='l', col=rgb(col85[1],col85[2],col85[3]), lwd=2);
  lines(log10(preturn.nofd.avg) , cost.nofd.avg/conv , type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2);
  lines(log10(preturn.avg)      , investment.avg/conv      , type='l', col=rgb(col85[1],col85[2],col85[3]), lwd=2, lty=6);
  lines(log10(preturn.nofd.avg) , investment.nofd.avg/conv , type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, lty=6);
  mtext(side=3, text=expression(bold(' a')), line=-2, cex=1, adj=0);
  points(log10(preturn.pre),cost.nofd.avg[iopt.nofd.avg]/conv,type='p',pch=16,cex=1.4)

  polygon(c(inset.x,rev(inset.x)),c(inset.y[1],inset.y[1],inset.y[2],inset.y[2])/conv,col=rgb(0,0,0,0),border=TRUE,lwd=1.5)
  legend(4.15,1.2,c("Total cost","Investment"), col=c('black','black'), lty=c(1,6), lwd=2, bty='n', cex=1.1)

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
py=.58
px=.28
v <- c( (v[1]+v[2])*px, v[2], (v[4]+v[3])*py, v[4] )

  lines( c(inset.x[2],u[2]), c(inset.y[1]/conv,v[3]*u[4]+2*u[3]), lty=1, lwd=1);
  lines( c(inset.x[1]-.05,v[1]*(u[2]-u[1])+u[1])+.05, c(inset.y[1]/conv,v[3]*u[4]+2*u[3]), lty=1, lwd=1);
  lines( c(log10(preturn.pre),log10(preturn.pre)), c(-9,investment.avg[ipre.nofd.avg]/conv), lty=2, lwd=2);
  lines( c(log10(preturn.avg[ipre.nofd.avg]), log10(preturn.avg[ipre.nofd.avg])),
         c(investment.avg[ipre.nofd.avg]/conv, -9), lty=2, lwd=2);

rect(u[2], u[4], (u[1]+u[2])*px, (u[4]-u[3])*py, col="white")
par( fig=v, new=TRUE, mar=c(0,0,0,0) )

plot(log10(preturn.avg),cost.avg/conv, col=rgb(col85[1],col85[2],col85[3]),
			type='l',xlim=inset.x,ylim=inset.y/conv, xlab='', ylab='', xaxt='n', yaxt='n');
      axis(1, at=seq(3.8,4.0,by=.2), label=parse(text=paste("10^", seq(3.8,4.0,by=.2), sep="")), las=1)
  polygon(c(rp.bin,rev(rp.bin)), c(cost.rp[,3],rev(cost.rp[,1]))/conv, col=rgb(col85[1],col85[2],col85[3],.3), border=NA);
  polygon(c(rp.bin,rev(rp.bin)), c(cost.nofd.rp[,3],rev(cost.nofd.rp[,1]))/conv, col=rgb(col45[1],col45[2],col45[3],.3), border=NA);
  lines(log10(preturn.avg)      , cost.avg/conv      , type='l', col=rgb(col85[1],col85[2],col85[3]), lwd=2);
  lines(log10(preturn.nofd.avg) , cost.nofd.avg/conv , type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2);
  lines(log10(preturn.avg)      , investment.avg/conv      , type='l', col=rgb(col85[1],col85[2],col85[3]), lwd=2, lty=6);
  lines(log10(preturn.nofd.avg) , investment.nofd.avg/conv , type='l', col=rgb(col45[1],col45[2],col45[3]), lwd=2, lty=6);

  #lines( c(log10(preturn.pre),log10(preturn.pre)), c(0,(log10(preturn.pre)*slope.nofd+inter.nofd))/conv, lty=2, lwd=2);
  lines( c(log10(preturn.pre),log10(preturn.pre)), c(0,cost.nofd.avg[iopt.nofd.avg])/conv, lty=2, lwd=2);

  xtmp = (log10(preturn.pre)*slope.nofd+inter.nofd-inter)/slope

  arrows( log10(preturn.pre), (log10(preturn.pre)*slope.nofd+inter.nofd)/conv,
          xtmp, (log10(preturn.pre)*slope.nofd+inter.nofd)/conv,
          code=2, length=.11, angle=45, lty=1, lwd=2);
  lines( c(xtmp, xtmp), c((log10(preturn.pre)*slope.nofd+inter.nofd)/conv, 0), lty=2, lwd=2);
  #points(log10(preturn.pre), (log10(preturn.pre)*slope.nofd+inter.nofd)/conv,type='p',pch=18,cex=2)
  points(log10(preturn.pre),cost.nofd.avg[iopt.nofd.avg]/conv,type='p',pch=16,cex=1.4)

  text(0.998*log10(preturn.pre), 0.985*mean(c(inset.y[1],investment.avg[ipre.nofd.avg]))/conv,
       bquote(atop("Goal:",paste("1:",.(prettyNum(preturn.pre,big.mark=',')),"y"))), pos = 4, cex = 1.1, srt = 0)
  text(0.95*log10(preturn.avg[ipre.nofd.avg]), .97*investment.avg[ipre.nofd.avg]/conv,
       bquote(atop("  Actual:",paste("1:",.(prettyNum(round(preturn.avg[ipre.nofd.avg]),big.mark=',')),"y"))), pos = 3, cex = 1.1, srt = 0)
#  text(0.997*xtmp, 0.968*mean(c(inset.y[1],investment.avg[ipre.nofd.avg]))/conv,
#        paste("Actual: 1:",round(preturn.avg[ipre.nofd.avg]),"y",sep=""), pos=4, cex=1.1, srt=0)

## HEIGHTENING

inset.x = c(1.34,1.61)
inset.y = c(2.22e9,2.38e9)
conv=1e9  # convert from $ to billions or millions of $? (for nicer looking axes)

par(fig=c(pfig,1,0,1), new=TRUE, mai=c(.7,.1,.06,.06))
plot(heightening,cost.avg/conv,col=rgb(col85[1],col85[2],col85[3],.2),type='l',xlim=c(0,3),ylim=c(0,8e9)/conv,
  xlab='', ylab='', xaxt='n', yaxt='n');
  mtext('Heightening [m]',side = 1,line = 2.1)
  axis(1, at=seq(0,3,by=0.5), label=c('0','','1','','2','','3'))
  axis(2, at=seq(0,8,by=2), label=rep("",5), las=1)
  polygon(c(heightening,rev(heightening)), c(cost.h[,3],rev(cost.h[,1]))/conv, col=rgb(col85[1],col85[2],col85[3],.3), border=NA);
  polygon(c(heightening,rev(heightening)), c(cost.nofd.h[,3],rev(cost.nofd.h[,1]))/conv, col=rgb(col45[1],col45[2],col45[3],.3), border=NA);
  lines(heightening,cost.avg/conv,type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=3);
  lines(heightening,cost.nofd.avg/conv,type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=3);
  points(heightening[iopt.nofd.avg],cost.nofd.avg[iopt.nofd.avg]/conv,type='p',pch=16,cex=1.)
#  points(heightening[iopt.nofd.avg],cost.avg[iopt.nofd.avg]/conv,type='p',pch=17,cex=1.)
  points(heightening[iopt.avg],cost.avg[iopt.avg]/conv,type='p',pch=15,cex=1.)
  polygon(c(inset.x,rev(inset.x)),c(inset.y[1]-.15e9,inset.y[1]-.15e9,inset.y[2]+.15e9,inset.y[2]+.15e9)/conv,col=rgb(0,0,0,0),border=TRUE,lwd=1.5)
  legend(.35,1.2,c("fast dynamics included", "fast dynamics neglected"),
         lty=c(1,1), lwd=2, col=c(rgb(col85[1],col85[2],col85[3]),rgb(col45[1],col45[2],col45[3])), cex=1.1, bty='n')
  mtext(side=3, text=expression(bold(' b')), line=-2, cex=1, adj=0);

u <- par("usr")
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
py=1.1/2
px=.9/2
v <- c( (v[1]+v[2])*px, v[2], (v[4]+v[3])*py, v[4] )

  lines( c(inset.x[2],u[2]), c(inset.y[1]/conv+.5*u[3],v[3]*u[4]+2*u[3]), lty=1, lwd=1);
  lines( c(inset.x[1],(u[2]-u[1])*v[1]*(1-pfig)+u[1]), c(inset.y[1]/conv+.5*u[3],v[3]*u[4]+2*u[3]), lty=1, lwd=1);

rect(u[2], u[4], (u[1]+u[2])*px-px, (u[4]-u[3])*py+u[3]*(1-py), col="white")
par( fig=v, new=TRUE, mar=c(0,0,0,0) )

plot(heightening,cost.avg/conv,type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=3, xlim=inset.x, ylim=inset.y/conv,
  xlab='', ylab='', yaxt='n', xaxt='n');
  polygon(c(heightening,rev(heightening)), c(cost.h[,3],rev(cost.h[,1]))/conv, col=rgb(col85[1],col85[2],col85[3],.3), border=NA);
  polygon(c(heightening,rev(heightening)), c(cost.nofd.h[,3],rev(cost.nofd.h[,1]))/conv, col=rgb(col45[1],col45[2],col45[3],.3), border=NA);
  lines(heightening,cost.nofd.avg/conv,type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=3);
  lines(heightening,cost.avg/conv,type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=3);

  lines(c(heightening[iopt.nofd.avg],heightening[iopt.nofd.avg]),c(0,cost.avg[iopt.nofd.avg])/conv,type='l',lty=2,lwd=2,col='black')
  lines(c(heightening[iopt.nofd.avg],heightening[iopt.avg]),c(cost.avg[iopt.nofd.avg],cost.avg[iopt.nofd.avg])/conv,type='l',lty=2,lwd=2,col='black')
  arrows(heightening[iopt.nofd.avg], cost.avg[iopt.avg]/conv, .995*heightening[iopt.avg], cost.avg[iopt.avg]/conv,
         length=.12, angle=40, lty=1, lwd=1.5, code=2, col='black')
  arrows(heightening[iopt.avg], cost.avg[iopt.nofd.avg]/conv, heightening[iopt.avg], cost.avg[iopt.avg]/conv,
         length=.11, angle=40, lty=1, lwd=1.5, code=1, col='black')

  points(heightening[iopt.nofd.avg],cost.nofd.avg[iopt.nofd.avg]/conv,type='p',pch=16,cex=1.4)
#  points(heightening[iopt.nofd.avg],cost.avg[iopt.nofd.avg]/conv,type='p',pch=17,cex=1.4)
  points(heightening[iopt.avg],cost.avg[iopt.avg]/conv,type='p',pch=15,cex=1.4)

  text(0.91*heightening[iopt.avg], 1.006*cost.avg[iopt.nofd.avg]/conv,
       paste("+$",waste.avg," million",sep=""), pos = 4, cex = 1.1, srt = 0)
  text(heightening[iopt.nofd.avg], .992*cost.avg[iopt.avg]/conv,
       paste("+",incheight.avg," m",sep=""), pos = 4, cex = 1.1, srt = 0)

#  text(1.0*heightening[iopt.nofd.avg], .99*mean(c(cost.avg[iopt.nofd.avg],cost.nofd.avg[iopt.nofd.avg]))/conv,
#       paste("+$",overspending.avg," million",sep=""), pos = 4, cex = 1.1, srt = 0)

dev.off()

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## SI 1 -- FIGURE 5 -- CALIBRATED PARAMETER MARGINAL DISTRIBUTIONS
##=================


dat.uni = read.csv(filename.parameters.uniform)
parameters.uni = dat.uni[1:(nrow(dat.uni)-1), ]
parnames = colnames(parameters.uni)

dat.gam = read.csv(filename.parameters.gamma)
parameters.gam = dat.gam[1:(nrow(dat.gam)-1), ]

## Fit PDFs to the parameter distributions
np = ncol(parameters.uni)
bound.lower = rep(0, np)
bound.upper = rep(0, np)
for (pp in 1:np) {
	bound.lower[pp] = min(c(parameters.uni[,pp],parameters.gam[,pp])) - 0.0*(max(c(parameters.uni[,pp],parameters.gam[,pp]))-min(c(parameters.uni[,pp],parameters.gam[,pp])))
  bound.upper[pp] = max(c(parameters.uni[,pp],parameters.gam[,pp])) + 0.0*(max(c(parameters.uni[,pp],parameters.gam[,pp]))-min(c(parameters.uni[,pp],parameters.gam[,pp])))
}

pdf.uni=vector('list',np)
pdf.gam=vector('list',np)
n.node=100
for (pp in 1:np){
  tmp = density(parameters.uni[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.uni[[pp]] = tmp; names(pdf.uni)[pp]=colnames(parameters.uni)[pp]
  tmp = density(parameters.gam[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.gam[[pp]] = tmp; names(pdf.gam)[pp]=colnames(parameters.gam)[pp]
}

pdf(paste(plotdir,'parameters_marginal_distributions.pdf',sep=''),width=6,height=9)#,colormodel='cmyk')

par(mfrow=c(7,6), mai=c(.6,.3,.1,.1))
for (pp in 1:np) {
	if(pp!=19) {plot(pdf.uni[[pp]]$x, pdf.uni[[pp]]$y, type='l', main="",
									xlab=parnames[pp], ylab='', cex.axis=1.2, cex.lab=1.2, yaxt="n",
                  ylim=c(0,max(c(pdf.uni[[pp]]$y,pdf.gam[[pp]]$y))))
              lines(pdf.gam[[pp]]$x, pdf.gam[[pp]]$y, type='l', col='red', lty=2)
	} else {plot(pdf.uni[[pp]]$x, pdf.uni[[pp]]$y, type='l', main="",
									xlab=parnames[pp], ylab='', cex.axis=1.2, cex.lab=1.2, yaxt="n",
                  ylim=c(0,max(c(pdf.uni[[pp]]$y,pdf.gam[[pp]]$y))))
          lines(pdf.gam[[pp]]$x, pdf.gam[[pp]]$y, type='l', col='red', lty=2)
					mtext('Density', side = 2, line = -1.3, outer = TRUE, at = NA, cex = 1.2)

	}
}
par(mai=c(.1,.1,.1,.1)); plot(1, type="n", axes=F, xlab="", ylab="");
  legend(0.35,1.4,c("Uniform",'  priors',"Gamma",'  priors'),
         col=c('black',NA,'red',NA), lty=c(1,NA,2,NA), bty='n', cex=1.2, lwd=1.3)
dev.off()

##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## SI 2, TABLE 1 -- PRIOR RANGES AND 90% CI FOR PARAMETERS
##==============

luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # Thermal expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = TRUE    # Antarctic ice sheet model
luse.brick = cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais)
source('../calibration/BRICK_parameterSetup_fastdyn.R') # get bound.lower/upper
parnames.bounds=parnames

dat.uni = read.csv(filename.parameters.uniform)
parameters.uni = dat.uni[1:(nrow(dat.uni)-1), ]
parnames = colnames(parameters.uni)

rho.simple.fixed = signif(as.numeric(read.csv(filename.rho_simple_fixed)),digits=4)

parameters=parameters.gam

source('../Useful/MultipleOutput.R') # defines the ":=" operator

p05=rep(0, ncol(parameters))
p50=rep(0, ncol(parameters))
p95=rep(0, ncol(parameters))
for (i in 1:ncol(parameters)) {
	c(p05[i],p50[i],p95[i]) := quantile( parameters[,i], c(.05,.50,.95))
}

p05=signif(p05,digits=4)
p50=signif(p50,digits=4)
p95=signif(p95,digits=4)

ci90=cbind(parnames,p50,p05,p95)

tmp = rbind( ci90[c(1:5,19:22,6:9,23:24,10:13,14:18,25,26:41),])
lb = rep(0,nrow(tmp)); ub = lb
for (i in 1:nrow(tmp)){
  itmp = match(parnames.bounds[i],tmp[,1])
  lb[itmp] = bound.lower[i]; ub[itmp] = bound.upper[i]
}
tmp = cbind(tmp,lb,ub)
ci90.rearr = rbind( tmp[1:match('sigma.simple',tmp[,1]),],
                    c("rho.simple",rho.simple.fixed,rho.simple.fixed,rho.simple.fixed,rho.simple.fixed,rho.simple.fixed),
                    tmp[(match('sigma.simple',tmp[,1])+1):41,])

models = c('DOECLIM','DOECLIM','DOECLIM','DOECLIM','DOECLIM','DOECLIM','DOECLIM','DOECLIM','DOECLIM',
           'GSIC','GSIC','GSIC','GSIC','GSIC','GSIC',
           'TE','TE','TE','TE',
           'SIMPLE-GIS','SIMPLE-GIS','SIMPLE-GIS','SIMPLE-GIS','SIMPLE-GIS','SIMPLE-GIS','SIMPLE-GIS',
           'AntOc','AntOc',
           'DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS-FD','DAIS-FD','DAIS')
units = c('deg C','cm^2 s^{-1}','-','deg C','10^{22} J','deg C','10^{22} J','-','-',
          'm yr^{-1} deg C^{-1}','m SLE','-','m SLE','m SLE','-',
          'm deg C^{-1}','m SLE','yr^{-1}','m SLE',
          'm SLE deg C^{-1}','m SLE','yr^{-1} deg C^{-1}','yr^{-1}','m SLE','m SLE','-',
          'deg C deg C^{-1}','deg C',
          '-','-','m^{1/2}','m^{-1/2} yr^{-1/2}','m yr^{-1}','deg C^{-1}','m yr^{-1}','m','m deg C^{-1}','m','-','m yr^{-1}','deg C','m^2 SLE')
output.table = data.frame( cbind(ci90.rearr[,1], models, ci90.rearr[,2:6], units))
colnames(output.table) = c('Parameter','Model','Median','5% quantile','95% quantile','Lower bound','Upper bound','Units')

## Write calibrated parameters to a csv file
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename=paste(plotdir,'BRICK-fastdyn_ExtendedDataTable1_',today,'.csv', sep="")
write.table(output.table, file=filename, sep=",", qmethod="double", row.names=FALSE)

## Will need to save the file as .xlsx for Supplementary Information

## Latex version: (too large for a single page)
if(FALSE){
n.dig=4

row.S = paste('$S$ & DOECLIM & ',signif(p50[match("S",parnames)],n.dig),'&',signif(p05[match("S",parnames)],n.dig),'&',signif(p95[match("S",parnames)],n.dig),'&',signif(bound.lower[match("S",parnames.bounds)],n.dig),'&',signif(bound.upper[match("S",parnames.bounds)],n.dig),' ($\\degree$C)\\\\',sep='')
row.kappa.doeclim = paste('$\\kappa_{DOECLIM}$ & DOECLIM & ',signif(p50[match("kappa.doeclim",parnames)],n.dig),'&',signif(p05[match("kappa.doeclim",parnames)],n.dig),'&',signif(p95[match("kappa.doeclim",parnames)],n.dig),'&',signif(bound.lower[match("kappa.doeclim",parnames.bounds)],n.dig),'&',signif(bound.upper[match("kappa.doeclim",parnames.bounds)],n.dig),' (cm$^2$ s$^{-1}$)\\\\',sep='')
row.alpha.doeclim = paste('$\\alpha_{DOECLIM}$ & DOECLIM & ',signif(p50[match("alpha.doeclim",parnames)],n.dig),'&',signif(p05[match("alpha.doeclim",parnames)],n.dig),'&',signif(p95[match("alpha.doeclim",parnames)],n.dig),'&',signif(bound.lower[match("alpha.doeclim",parnames.bounds)],n.dig),'&',signif(bound.upper[match("alpha.doeclim",parnames.bounds)],n.dig),' (-)\\\\',sep='')
row.T0 = paste('$T_0$ & DOECLIM & ',signif(p50[match("T0",parnames)],n.dig),'&',signif(p05[match("T0",parnames)],n.dig),'&',signif(p95[match("T0",parnames)],n.dig),'&',signif(bound.lower[match("T0",parnames.bounds)],n.dig),'&',signif(bound.upper[match("T0",parnames.bounds)],n.dig),' ($\\degree$C)\\\\',sep='')
row.H0 = paste('$H_0$ & DOECLIM & ',signif(p50[match("H0",parnames)],n.dig),'&',signif(p05[match("H0",parnames)],n.dig),'&',signif(p95[match("H0",parnames)],n.dig),'&',signif(bound.lower[match("H0",parnames.bounds)],n.dig),'&',signif(bound.upper[match("H0",parnames.bounds)],n.dig),' ($10^{22}$ J)\\\\',sep='')
row.sigma.T = paste('$\\sigma_T$ & DOECLIM & ',signif(p50[match("sigma.T",parnames)],n.dig),'&',signif(p05[match("sigma.T",parnames)],n.dig),'&',signif(p95[match("sigma.T",parnames)],n.dig),'&',signif(bound.lower[match("sigma.T",parnames.bounds)],n.dig),'&',signif(bound.upper[match("sigma.T",parnames.bounds)],n.dig),' ($\\degree$C)\\\\',sep='')
row.sigma.H = paste('$\\sigma_H$ & DOECLIM & ',signif(p50[match("sigma.H",parnames)],n.dig),'&',signif(p05[match("sigma.H",parnames)],n.dig),'&',signif(p95[match("sigma.H",parnames)],n.dig),'&',signif(bound.lower[match("sigma.H",parnames.bounds)],n.dig),'&',signif(bound.upper[match("sigma.H",parnames.bounds)],n.dig),' ($10^{22}$ J)\\\\',sep='')
row.rho.T = paste('$\\rho_T$ & DOECLIM & ',signif(p50[match("rho.T",parnames)],n.dig),'&',signif(p05[match("rho.T",parnames)],n.dig),'&',signif(p95[match("rho.T",parnames)],n.dig),'&',signif(bound.lower[match("rho.T",parnames.bounds)],n.dig),'&',signif(bound.upper[match("rho.T",parnames.bounds)],n.dig),' (-)\\\\',sep='')
row.rho.H = paste('$\\rho_H$ & DOECLIM & ',signif(p50[match("rho.H",parnames)],n.dig),'&',signif(p05[match("rho.H",parnames)],n.dig),'&',signif(p95[match("rho.H",parnames)],n.dig),'&',signif(bound.lower[match("rho.H",parnames.bounds)],n.dig),'&',signif(bound.upper[match("rho.H",parnames.bounds)],n.dig),' (-)\\\\',sep='')
row.beta0 = paste('$\\beta_0$ & GSIC & ',signif(p50[match("beta0",parnames)],n.dig),'&',signif(p05[match("beta0",parnames)],n.dig),'&',signif(p95[match("beta0",parnames)],n.dig),'&',signif(bound.lower[match("beta0",parnames.bounds)],n.dig),'&',signif(bound.upper[match("beta0",parnames.bounds)],n.dig),' (m yr$^{-1}$ $\\degree$C$^{-1}$)\\\\',sep='')
row.V0.gsic = paste('$V_{0,GSIC}$ & GSIC & ',signif(p50[match("V0.gsic",parnames)],n.dig),'&',signif(p05[match("V0.gsic",parnames)],n.dig),'&',signif(p95[match("V0.gsic",parnames)],n.dig),'&',signif(bound.lower[match("V0.gsic",parnames.bounds)],n.dig),'&',signif(bound.upper[match("V0.gsic",parnames.bounds)],n.dig),' (m SLE)\\\\',sep='')
row.n = paste('$n$ & GSIC & ',signif(p50[match("n",parnames)],n.dig),'&',signif(p05[match("n",parnames)],n.dig),'&',signif(p95[match("n",parnames)],n.dig),'&',signif(bound.lower[match("n",parnames.bounds)],n.dig),'&',signif(bound.upper[match("n",parnames.bounds)],n.dig),' (-)\\\\',sep='')
row.Gs0 = paste('$G_{s,0}$ & GSIC & ',signif(p50[match("Gs0",parnames)],n.dig),'&',signif(p05[match("Gs0",parnames)],n.dig),'&',signif(p95[match("Gs0",parnames)],n.dig),'&',signif(bound.lower[match("Gs0",parnames.bounds)],n.dig),'&',signif(bound.upper[match("Gs0",parnames.bounds)],n.dig),' (m SLE)\\\\',sep='')
row.sigma.gsic = paste('$\\sigma_{GSIC}$ & GSIC & ',signif(p50[match("sigma.gsic",parnames)],n.dig),'&',signif(p05[match("sigma.gsic",parnames)],n.dig),'&',signif(p95[match("sigma.gsic",parnames)],n.dig),'&',signif(bound.lower[match("sigma.gsic",parnames.bounds)],n.dig),'&',signif(bound.upper[match("sigma.gsic",parnames.bounds)],n.dig),' (m SLE)\\\\',sep='')
row.rho.gsic = paste('$\\rho_{GSIC}$ & GSIC & ',signif(p50[match("rho.gsic",parnames)],n.dig),'&',signif(p05[match("rho.gsic",parnames)],n.dig),'&',signif(p95[match("rho.gsic",parnames)],n.dig),'&',signif(bound.lower[match("rho.gsic",parnames.bounds)],n.dig),'&',signif(bound.upper[match("rho.gsic",parnames.bounds)],n.dig),' (-)\\\\',sep='')
row.a.te = paste('$a_{TE}$ & TE & ',signif(p50[match("a.te",parnames)],n.dig),'&',signif(p05[match("a.te",parnames)],n.dig),'&',signif(p95[match("a.te",parnames)],n.dig),'&',signif(bound.lower[match("a.te",parnames.bounds)],n.dig),'&',signif(bound.upper[match("a.te",parnames.bounds)],n.dig),' (m $\\degree$C$^{-1}$)\\\\',sep='')
row.b.te = paste('$b_{TE}$ & TE & ',signif(p50[match("b.te",parnames)],n.dig),'&',signif(p05[match("b.te",parnames)],n.dig),'&',signif(p95[match("b.te",parnames)],n.dig),'&',signif(bound.lower[match("b.te",parnames.bounds)],n.dig),'&',signif(bound.upper[match("b.te",parnames.bounds)],n.dig),' (m SLE)\\\\',sep='')
row.invtau.te = paste('$1/\\tau_{TE}$ & TE & ',signif(p50[match("invtau.te",parnames)],n.dig),'&',signif(p05[match("invtau.te",parnames)],n.dig),'&',signif(p95[match("invtau.te",parnames)],n.dig),'&',signif(bound.lower[match("invtau.te",parnames.bounds)],n.dig),'&',signif(bound.upper[match("invtau.te",parnames.bounds)],n.dig),' (yr$^{-1}$)\\\\',sep='')
row.TE0 = paste('$TE_0$ & TE & ',signif(p50[match("TE0",parnames)],n.dig),'&',signif(p05[match("TE0",parnames)],n.dig),'&',signif(p95[match("TE0",parnames)],n.dig),'&',signif(bound.lower[match("TE0",parnames.bounds)],n.dig),'&',signif(bound.upper[match("TE0",parnames.bounds)],n.dig),' (m SLE)\\\\',sep='')
row.a.simple = paste('$a_{GIS}$ & SIMPLE-GIS & ',signif(p50[match("a.simple",parnames)],n.dig),'&',signif(p05[match("a.simple",parnames)],n.dig),'&',signif(p95[match("a.simple",parnames)],n.dig),'&',signif(bound.lower[match("a.simple",parnames.bounds)],n.dig),'&',signif(bound.upper[match("a.simple",parnames.bounds)],n.dig),' (m SLE $\\degree$C$^{-1}$)\\\\',sep='')
row.b.simple = paste('$b_{GIS}$ & SIMPLE-GIS & ',signif(p50[match("b.simple",parnames)],n.dig),'&',signif(p05[match("b.simple",parnames)],n.dig),'&',signif(p95[match("b.simple",parnames)],n.dig),'&',signif(bound.lower[match("b.simple",parnames.bounds)],n.dig),'&',signif(bound.upper[match("b.simple",parnames.bounds)],n.dig),' (m SLE)\\\\',sep='')
row.alpha.simple = paste('$\\alpha_{GIS}$ & SIMPLE-GIS & ',signif(p50[match("alpha.simple",parnames)],n.dig),'&',signif(p05[match("alpha.simple",parnames)],n.dig),'&',signif(p95[match("alpha.simple",parnames)],n.dig),'&',signif(bound.lower[match("alpha.simple",parnames.bounds)],n.dig),'&',signif(bound.upper[match("alpha.simple",parnames.bounds)],n.dig),' (yr$^{-1}$ \\degree C$^{-1}$)\\\\',sep='')
row.beta.simple = paste('$\\beta_{GIS}$ & SIMPLE-GIS & ',signif(p50[match("beta.simple",parnames)],n.dig),'&',signif(p05[match("beta.simple",parnames)],n.dig),'&',signif(p95[match("beta.simple",parnames)],n.dig),'&',signif(bound.lower[match("beta.simple",parnames.bounds)],n.dig),'&',signif(bound.upper[match("beta.simple",parnames.bounds)],n.dig),' (yr$^{-1}$)\\\\',sep='')
row.V0 = paste('$V_{0,GIS}$ & SIMPLE-GIS & ',signif(p50[match("V0",parnames)],n.dig),'&',signif(p05[match("V0",parnames)],n.dig),'&',signif(p95[match("V0",parnames)],n.dig),'&',signif(bound.lower[match("V0",parnames.bounds)],n.dig),'&',signif(bound.upper[match("V0",parnames.bounds)],n.dig),' (m SLE)\\\\',sep='')
row.sigma.simple = paste('$\\sigma_{GIS}$ & SIMPLE-GIS & ',signif(p50[match("sigma.simple",parnames)],n.dig),'&',signif(p05[match("sigma.simple",parnames)],n.dig),'&',signif(p95[match("sigma.simple",parnames)],n.dig),'&',signif(bound.lower[match("sigma.simple",parnames.bounds)],n.dig),'&',signif(bound.upper[match("sigma.simple",parnames.bounds)],n.dig),' (m SLE)\\\\',sep='')
row.rho.simple = paste('$\\rho_{GIS}$ & SIMPLE-GIS & ',signif(rho.simple.fixed,n.dig),'& * & * & * & * (-) \\\\',sep='')
row.anto.a = paste('$a_{AntOc}$ & AntOc & ',signif(p50[match("anto.a",parnames)],n.dig),'&',signif(p05[match("anto.a",parnames)],n.dig),'&',signif(p95[match("anto.a",parnames)],n.dig),'&',signif(bound.lower[match("anto.a",parnames.bounds)],n.dig),'&',signif(bound.upper[match("anto.a",parnames.bounds)],n.dig),' ($\\degree$C $T_{AntOc}$ $\\degree$C$^{-1}$ $T_{g}$)\\\\',sep='')
row.anto.b = paste('$b_{AntOc}$ & AntOc & ',signif(p50[match("anto.b",parnames)],n.dig),'&',signif(p05[match("anto.b",parnames)],n.dig),'&',signif(p95[match("anto.b",parnames)],n.dig),'&',signif(bound.lower[match("anto.b",parnames.bounds)],n.dig),'&',signif(bound.upper[match("anto.b",parnames.bounds)],n.dig),' ($\\degree$C $T_{AntOc}$)\\\\',sep='')
row.gamma = paste('$\\gamma$ & DAIS & ',signif(p50[match("gamma",parnames)],n.dig),'&',signif(p05[match("gamma",parnames)],n.dig),'&',signif(p95[match("gamma",parnames)],n.dig),'&',signif(bound.lower[match("gamma",parnames.bounds)],n.dig),'&',signif(bound.upper[match("gamma",parnames.bounds)],n.dig),' (-)\\\\',sep='')
row.alpha.dais = paste('$\\alpha_{DAIS}$ & DAIS & ',signif(p50[match("alpha.dais",parnames)],n.dig),'&',signif(p05[match("alpha.dais",parnames)],n.dig),'&',signif(p95[match("alpha.dais",parnames)],n.dig),'&',signif(bound.lower[match("alpha.dais",parnames.bounds)],n.dig),'&',signif(bound.upper[match("alpha.dais",parnames.bounds)],n.dig),' (-)\\\\',sep='')
row.mu = paste('$\\mu$ & DAIS & ',signif(p50[match("mu",parnames)],n.dig),'&',signif(p05[match("mu",parnames)],n.dig),'&',signif(p95[match("mu",parnames)],n.dig),'&',signif(bound.lower[match("mu",parnames.bounds)],n.dig),'&',signif(bound.upper[match("mu",parnames.bounds)],n.dig),' (m$^{1/2}$)\\\\',sep='')
row.nu = paste('$\\nu$ & DAIS & ',signif(p50[match("nu",parnames)],n.dig),'&',signif(p05[match("nu",parnames)],n.dig),'&',signif(p95[match("nu",parnames)],n.dig),'&',signif(bound.lower[match("nu",parnames.bounds)],n.dig),'&',signif(bound.upper[match("nu",parnames.bounds)],n.dig),' (m$^{-1/2}$ yr$^{-1/2}$)\\\\',sep='')
row.P0 = paste('$P_0$ & DAIS & ',signif(p50[match("P0",parnames)],n.dig),'&',signif(p05[match("P0",parnames)],n.dig),'&',signif(p95[match("P0",parnames)],n.dig),'&',signif(bound.lower[match("P0",parnames.bounds)],n.dig),'&',signif(bound.upper[match("P0",parnames.bounds)],n.dig),' (m yr$^{-1}$)\\\\',sep='')
row.kappa.dais = paste('$\\kappa_{DAIS} $ & DAIS & ',signif(p50[match("kappa.dais",parnames)],n.dig),'&',signif(p05[match("kappa.dais",parnames)],n.dig),'&',signif(p95[match("kappa.dais",parnames)],n.dig),'&',signif(bound.lower[match("kappa.dais",parnames.bounds)],n.dig),'&',signif(bound.upper[match("kappa.dais",parnames.bounds)],n.dig),' ($\\degree$C$^{-1}$)\\\\',sep='')
row.f0 = paste('$f_0$ & DAIS & ',signif(p50[match("f0",parnames)],n.dig),'&',signif(p05[match("f0",parnames)],n.dig),'&',signif(p95[match("f0",parnames)],n.dig),'&',signif(bound.lower[match("f0",parnames.bounds)],n.dig),'&',signif(bound.upper[match("f0",parnames.bounds)],n.dig),' (m yr$^{-1}$)\\\\',sep='')
row.h0 = paste('$h_0$ & DAIS & ',signif(p50[match("h0",parnames)],n.dig),'&',signif(p05[match("h0",parnames)],n.dig),'&',signif(p95[match("h0",parnames)],n.dig),'&',signif(bound.lower[match("h0",parnames.bounds)],n.dig),'&',signif(bound.upper[match("h0",parnames.bounds)],n.dig),' (m)\\\\',sep='')
row.c = paste('$c$ & DAIS & ',signif(p50[match("c",parnames)],n.dig),'&',signif(p05[match("c",parnames)],n.dig),'&',signif(p95[match("c",parnames)],n.dig),'&',signif(bound.lower[match("c",parnames.bounds)],n.dig),'&',signif(bound.upper[match("c",parnames.bounds)],n.dig),' (m $\\degree$C$^{-1}$)\\\\',sep='')
row.b0 = paste('$b_0$ & DAIS & ',signif(p50[match("b0",parnames)],n.dig),'&',signif(p05[match("b0",parnames)],n.dig),'&',signif(p95[match("b0",parnames)],n.dig),'&',signif(bound.lower[match("b0",parnames.bounds)],n.dig),'&',signif(bound.upper[match("b0",parnames.bounds)],n.dig),' (m)\\\\',sep='')
row.slope = paste('$s$ & DAIS & ',signif(p50[match("slope",parnames)],n.dig),'&',signif(p05[match("slope",parnames)],n.dig),'&',signif(p95[match("slope",parnames)],n.dig),'&',signif(bound.lower[match("slope",parnames.bounds)],n.dig),'&',signif(bound.upper[match("slope",parnames.bounds)],n.dig),' (-)\\\\',sep='')
row.lambda = paste('$\\lambda$ & DAIS-FD & ',signif(p50[match("lambda",parnames)],n.dig),'&',signif(p05[match("lambda",parnames)],n.dig),'&',signif(p95[match("lambda",parnames)],n.dig),'&',signif(bound.lower[match("lambda",parnames.bounds)],n.dig),'&',signif(bound.upper[match("lambda",parnames.bounds)],n.dig),' (m yr$^{-1}$)\\\\',sep='')
row.Tcrit = paste('$T_{crit}$ & DAIS-FD & ',signif(p50[match("Tcrit",parnames)],n.dig),'&',signif(p05[match("Tcrit",parnames)],n.dig),'&',signif(p95[match("Tcrit",parnames)],n.dig),'&',signif(bound.lower[match("Tcrit",parnames.bounds)],n.dig),'&',signif(bound.upper[match("Tcrit",parnames.bounds)],n.dig),' ($\\degree$C)\\\\',sep='')
row.var.dais = paste('$\\sigma^2_{DAIS}$ & DAIS & ',signif(p50[match("var.dais",parnames)],n.dig),'&',signif(p05[match("var.dais",parnames)],n.dig),'&',signif(p95[match("var.dais",parnames)],n.dig),'&',signif(bound.lower[match("var.dais",parnames.bounds)],n.dig),'&',signif(bound.upper[match("var.dais",parnames.bounds)],n.dig),' (m$^2$ SLE)\\\\',sep='')

table=data.frame(rbind(cat(row.S),cat(row.kappa.doeclim),cat(row.alpha.doeclim),cat(row.T0),cat(row.H0),cat(row.sigma.T),cat(row.sigma.H),cat(row.rho.T),cat(row.rho.H)),
											cat(row.beta0),cat(row.V0.gsic),cat(row.n),cat(row.Gs0),cat(row.sigma.gsic),cat(row.rho.gsic),
											cat(row.a.te),cat(row.b.te),cat(row.invtau.te),cat(row.TE0),
											cat(row.a.simple),cat(row.b.simple),cat(row.alpha.simple),cat(row.beta.simple),cat(row.V0),cat(row.sigma.simple),cat(row.rho.simple),
											cat(row.anto.a),cat(row.anto.b),
											cat(row.gamma),cat(row.alpha.dais),cat(row.mu),cat(row.nu),cat(row.P0),cat(row.kappa.dais),cat(row.f0),cat(row.h0),cat(row.c),cat(row.b0),cat(row.slope),cat(row.lambda),cat(row.Tcrit),cat(row.var.dais),
											row.names=NULL)
}
##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## SI 3 -- FIGURE 6 -- Like Figure 2, but with both uniform and gamma priors
##=================

pdf(paste(plotdir,'projections_SLR_total_bothPriors.pdf',sep=''),width=6.2,height=5.5)#,colormodel='cmyk')
par(mfrow=c(2,2))
# UNIFORM RCP85
par(mai=c(.3,.75,.25,0))
plot(t.proj[iproj],gsl.rcp85.uni.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2100), ylim=c(0,2), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2100,by=20)); axis(2, seq(0,2,by=.25)); #axis(3, seq(2000,2100,by=20)); axis(4, seq(0,2,by=.25));
		 mtext(side=2, text='Total sea level [m]', line=2.2, cex=.9);
     mtext(side=3, text='Uniform priors', line=0.3, cex=.9);
     mtext(side=3, text=expression(bold(' a')), line=-1, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp85.uni.95[iproj],rev(gsl.rcp85.uni.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + UNIFORM RCP26
	lines(t.proj[iproj],gsl.rcp26.uni.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp26.uni.95[iproj],rev(gsl.rcp26.uni.05[iproj])),
          col=rgb(col26[1],col26[2],col26[3],.5), border=NA);
# + UNIFORM RCP45
	lines(t.proj[iproj],gsl.rcp45.uni.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp45.uni.95[iproj],rev(gsl.rcp45.uni.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
# + legend
  legend(t.proj[iproj[1]]+10,2,c("5-95% range,",
                                "RCP2.6",
																"RCP4.5",
																"RCP8.5"),
         lty=c(NA,1,1,1), lwd=3, col=c(NA,rgb(col26[1],col26[2],col26[3]),rgb(col45[1],col45[2],col45[3]),rgb(col85[1],col85[2],col85[3])), bty='n', cex=.9)
# GAMMA RCP85
par(mai=c(.3,.6,.25,.2))
plot(t.proj[iproj],gsl.rcp85.gam.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2100), ylim=c(0,2), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2100,by=20)); axis(2, seq(0,2,by=.25)); #axis(3, seq(2000,2100,by=20)); axis(4, seq(0,2,by=.25));
     mtext(side=3, text='Gamma priors', line=0.3, cex=.9);
     mtext(side=3, text=expression(bold(' b')), line=-1, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp85.gam.95[iproj],rev(gsl.rcp85.gam.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],.5), border=NA);
# + GAMMA RCP26
	lines(t.proj[iproj],gsl.rcp26.gam.50[iproj],type='l',col=rgb(col26[1],col26[2],col26[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp26.gam.95[iproj],rev(gsl.rcp26.gam.05[iproj])),
          col=rgb(col26[1],col26[2],col26[3],0.5), border=NA);
# + GAMMA RCP45
	lines(t.proj[iproj],gsl.rcp45.gam.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsl.rcp45.gam.95[iproj],rev(gsl.rcp45.gam.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);
# UNIFORM RCP85, FD contribution
par(mai=c(.6,.75,.1,0))
plot(t.proj[iproj],fastdyn.rcp85.uni.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2100), ylim=c(0,0.75), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2100,by=20)); axis(2, seq(0,2,by=.2)); #axis(3, seq(2000,2100,by=20)); axis(4, seq(0,2,by=.2));
		 mtext(side=1, text='Year', line=2.1, cex=.9);
		 mtext(side=2, text='Fast dynamics sea \nlevel contribution [m]', line=2.2, cex=.9);
     mtext(side=3, text=expression(bold(' c')), line=-1, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp85.uni.95[iproj],rev(fastdyn.rcp85.uni.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],0.5), border=NA);
# + UNIFORM RCP45, FD contribution
	lines(t.proj[iproj],fastdyn.rcp45.uni.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp45.uni.95[iproj],rev(fastdyn.rcp45.uni.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],0.5), border=NA);
# GAMMA RCP85, FD contribution
par(mai=c(.6,.6,.1,.2))
plot(t.proj[iproj],fastdyn.rcp85.gam.50[iproj],type='l',col=rgb(col85[1],col85[2],col85[3]),lwd=2, ann='',
		 xlim=c(2000,2100), ylim=c(0,0.75), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2100,by=20)); axis(2, seq(0,2,by=.2)); #axis(3, seq(2000,2100,by=20)); axis(4, seq(0,2,by=.2));
		 mtext(side=1, text='Year', line=2.1, cex=.9);
     mtext(side=3, text=expression(bold(' d')), line=-1, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp85.gam.95[iproj],rev(fastdyn.rcp85.gam.05[iproj])),
          col=rgb(col85[1],col85[2],col85[3],0.5), border=NA);
# + GAMMA RCP45, FD contribution
	lines(t.proj[iproj],fastdyn.rcp45.gam.50[iproj],type='l',col=rgb(col45[1],col45[2],col45[3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(fastdyn.rcp45.gam.95[iproj],rev(fastdyn.rcp45.gam.05[iproj])),
          col=rgb(col45[1],col45[2],col45[3],.5), border=NA);

dev.off()

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## FIGURE 1 -- Distribution of Tcrit, scaled to global mean surface temperature
##=========

dat.uni = read.csv(filename.parameters.uniform)
parameters.uni = dat.uni[1:(nrow(dat.uni)-1), ]
parnames = colnames(parameters.uni)

dat.gam = read.csv(filename.parameters.gamma)
parameters.gam = dat.gam[1:(nrow(dat.gam)-1), ]

Tacrit.uni = parameters.uni[,match("Tcrit",parnames)]
Tacrit.gam = parameters.gam[,match("Tcrit",parnames)]

Tgcrit.uni = slope.Ta2Tg*Tacrit.uni + intercept.Ta2Tg
Tgcrit.gam = slope.Ta2Tg*Tacrit.gam + intercept.Ta2Tg

quan.Tgcrit.uni = quantile(Tgcrit.uni,c(0,.05,.50,.95,.99,1))
quan.Tgcrit.gam = quantile(Tgcrit.gam,c(0,.05,.50,.95,.99,1))

pdf.Tgcrit.uni = density(Tgcrit.uni,from=0,to=max(c(Tgcrit.uni,Tgcrit.gam)),na.rm=TRUE)
pdf.Tgcrit.gam = density(Tgcrit.gam,from=0,to=max(c(Tgcrit.uni,Tgcrit.gam)),na.rm=TRUE)
#pdf.Tgcrit.uni = density(Tgcrit.uni,from=min(c(Tgcrit.uni,Tgcrit.gam)),to=max(c(Tgcrit.uni,Tgcrit.gam)),na.rm=TRUE)
#pdf.Tgcrit.gam = density(Tgcrit.gam,from=min(c(Tgcrit.uni,Tgcrit.gam)),to=max(c(Tgcrit.uni,Tgcrit.gam)),na.rm=TRUE)

cdf.Tgcrit.uni = rep(0,length(pdf.Tgcrit.uni$x))
cdf.Tgcrit.gam = rep(0,length(pdf.Tgcrit.gam$x))

dx=median(diff(pdf.Tgcrit.uni$x))
x=pdf.Tgcrit.uni$x

for (i in 2:length(x)){
  cdf.Tgcrit.uni[i] = sum(c(cdf.Tgcrit.uni[i-1],dx*pdf.Tgcrit.uni$y[i]))
  cdf.Tgcrit.gam[i] = sum(c(cdf.Tgcrit.gam[i-1],dx*pdf.Tgcrit.gam$y[i]))
}

cdf.Tgcrit.uni = cdf.Tgcrit.uni/cdf.Tgcrit.uni[length(cdf.Tgcrit.uni)]
cdf.Tgcrit.gam = cdf.Tgcrit.gam/cdf.Tgcrit.gam[length(cdf.Tgcrit.gam)]

sur.Tgcrit.uni = 1-cdf.Tgcrit.uni
sur.Tgcrit.gam = 1-cdf.Tgcrit.gam

## Calculate the quantile at the COP21 2-degree warming
tmp=abs(x-2);
print('================================================================')
print(paste('At 2-degree warming, have CDF=',cdf.Tgcrit.gam[which(tmp==min(tmp))]))
print(paste('median (5-95% CI) =',quantile(Tgcrit.gam,c(.5,.05,.95))))
print('================================================================')

## Make sure obs.temp is normalized relative to 1850-1870
obs.temp = obs.temp - mean(obs.temp[1:20])
temp.2015 = obs.temp[which(obs.temp.time==2015)]

pdf(paste(plotdir,'distributions_Tgcrit.pdf',sep=''),width=5,height=3.5,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.9,.74,.05,.35))

plot(x,pdf.Tgcrit.gam$y, type='l', xlim=c(-.2,4.1), ylim=c(0,1.5),lty=1,
     col=rgb(col45[1],col45[2],col45[3]), lwd=1.5, xlab='', ylab='', yaxt='n', axes=FALSE,yaxs='i');
  axis(1,seq(-0.5,4.5,0.5),lab=c("","0","","1","","2","","3","","4",""))
  u <- par("usr")
  arrows(u[1], u[3], u[1], u[4], code = 2, xpd = TRUE)
  mtext('Probability density', side=2, line=1.3);
  mtext('Global surface temperature anomaly (relative to 1850-1870\naverage) triggering fast dynamics disintegration [deg C]', side=1, line=3);

  lines(c(0,0),c(-10,1.38), type='l', col='black', lwd=1.5);
  text(0.5,1.35,'Pre-industrial', pos=3, cex=1.1, srt=0);

  lines(c(temp.2015,temp.2015),c(-10,1.25), type='l', col='black', lwd=1.5);
  text(.8,1.3,'Current (2015)', pos=4, cex=1.1, srt=0);

  lines(c(2,2),c(-10,1.1), type='l', col='black', lwd=1.5);
  text(2.18,1.07,'COP21', pos=3, cex=1.1, srt=0);

  itmp=which(pdf.Tgcrit.gam$x<2); xtmp=pdf.Tgcrit.gam$x[itmp]; ytmp=pdf.Tgcrit.gam$y[itmp]
  polygon(c(xtmp,rev(xtmp)), c(ytmp,rep(0,length(ytmp))),
          col=rgb(col85[1],col85[2],col85[3],0.5), border=NA);

dev.off()

##==============================================================================
##==============================================================================








##==============================================================================
##==============================================================================
## EXTENDED DATA FIGURE 1 -- JOINT MARGINAL DISTRIBUTION OF TCRIT AND LAMBDA
##                            AFTER PALEO CALIBRATION ONLY
##=======================

dat.dais = read.csv(filename.DAIScalibration)
parameters.dais = dat.dais[1:(nrow(dat.dais)-1),]
bandwidths.dais = dat.dais[nrow(dat.dais)      ,]
parnames.dais   = colnames(parameters.dais)

Tcrit = slope.Ta2Tg*parameters.dais[,match("Tcrit",parnames.dais)]+intercept.Ta2Tg
lambda= parameters.dais[,match("lambda",parnames.dais)]

#itmp=seq(1,nrow(Tcrit),by=100) # thin so plotting is faster/smaller file sizes
fdpar=data.frame(Tcrit,lambda)

## Also draw an empirical joint marginal prior distribution for the two FD parameters
## (from gamma)
shape.lambda = 8.1              # gives 5% quantile at lambda=0.005 and
rate.lambda = 100*shape.lambda  # gives mean at 0.01 m/yr, DeConto and Pollard (2016)
rate.Tcrit = 1.37               # gives 5% quantile at Tcrit = -10 deg C
shape.Tcrit = 15*rate.Tcrit     # gives mean at -15 deg C (negative requires multiplication of Tcrit by -1)

lambda.emp = rgamma(n=length(Tcrit), shape=shape.lambda, rate=rate.lambda)
Tcrit.emp  = slope.Ta2Tg*(-rgamma(n=length(Tcrit), shape=shape.Tcrit,  rate=rate.Tcrit))+intercept.Ta2Tg

## preliminary plot to get the natural bounds
h2 <- hist2d(fdpar, nbins=70, FUN=function(x) log(length(x)))
lims = par("usr")
#itmp=which(Tcrit.emp>lims[1] & Tcrit.emp<lims[2] & lambda.emp>lims[3] & lambda.emp<lims[4])
itmp=which((Tcrit.emp-intercept.Ta2Tg)/slope.Ta2Tg> -25 & (Tcrit.emp-intercept.Ta2Tg)/slope.Ta2Tg< -5 & lambda.emp>lims[3] & lambda.emp<lims[4]) # Tcrit too tightly constrained!
fdpar.emp = data.frame(Tcrit.emp[itmp],lambda.emp[itmp])

library(gplots)
library(fields)

pdf(paste(plotdir,'paleoConstraint_FDparameters.pdf',sep=''),width=4,height=4)#,colormodel='cmyk')

par(mfrow=c(1,1))
par(fig=c(0,1,.5,1),mai=c(.4,.65,.2,.1))
  h1 <- hist2d(fdpar.emp, nbins=400, FUN=function(x) log(length(x)), xlim=c(0,7), col=topo.colors(5))
  mtext(expression(lambda (m/yr)), side=2, line=2);
  mtext(expression(bold(' a')), side=3, line=-1, adj=0);
par(fig=c(0,.66,0,.5),mai=c(.65,.65,.01,.01),new=TRUE)
  h2 <- hist2d(fdpar, nbins=200, FUN=function(x) log(length(x)), xlim=c(0,4.2), col=topo.colors(5))
  mtext(expression(T[crit]*(degC)), side=1, line=2.2);
  mtext(expression(lambda (m/yr)), side=2, line=2);
  mtext(expression(bold(' b')), side=3, line=-1, adj=0);
par(fig=c(0,1,0,.5))
  image.plot(zlim=c(0,5), legend.only=TRUE,col=topo.colors(1000), legend.shrink = .6, legend.lab=NULL,
            legend.width = 1.5, horizontal = FALSE,
            axis.args=list(at=c(0,5),labels=c("Less\nlikely", "More\nlikely")))

dev.off()


##==============================================================================
##==============================================================================








##==============================================================================
## End
##==============================================================================
