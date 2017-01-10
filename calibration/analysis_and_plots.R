##==============================================================================
## plots and tables for robust SLR paper
##
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

rm(list=ls())

## Initial set-up
library(ncdf4)

## File name for the BRICK physical model output (netCDF4)
filename.brick  = '../output_model/BRICK-robustSLR_physical_14Sep2016.nc'

## File name for the BRICK post-calibrated parameters (csv) (the BRICK output came from these guys)
filename.parameters = '../output_calibration/BRICK-robustSLR_postcalibratedParameters_14Sep2016.csv'

## Other files
filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_06Sep2016.csv"

## Get nice plotting colors: mycol array
source('../Useful/colorblindPalette.R')

## set up the RCP forcings' default colors within mycol
c85=6;
c45=4;
c26=2;

## Where would you like to save the plots?
plotdir='~/Box\ Sync/Robust_SLR_Projections/Figures/'

##==============================================================================
##==============================================================================
## FIGURE 1 -- HINDCASTS VS OBSERVATIONS
##=========

ncdata <- nc_open(filename.brick)
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
plot(mod.time[midx.temp], temp.50[midx.temp], type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(1850,2016), ylim=c(-.3,1.5), cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Surface temperature\n[deg C]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' a')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time[midx.temp],rev(mod.time[midx.temp])), c(temp.95,rev(temp.05)), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
  lines(obs.temp.time[oidx.temp], obs.temp.norm[oidx.temp], type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.temp.time[oidx.temp],rev(obs.temp.time[oidx.temp])),
					c(obs.temp.norm[oidx.temp]+n.sig*obs.temp.err[oidx.temp],rev(obs.temp.norm[oidx.temp]-n.sig*obs.temp.err[oidx.temp])),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);
  #legend(1839,1.65,c("5-95% range, modeled","2-sigma range, observations"),
  #       col=c(rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),rgb(mycol[6,1],mycol[6,2],mycol[6,3])), lwd=2, bty='n', cex=1.2)

# >>> OCEAN HEAT <<<
itmp=midx.ocheat[1]:nrow(ocheat.hind)
plot(mod.time[itmp], ocheat.50[itmp], type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-50,50), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Ocean heat uptake\n[10^22 J]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' b')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(ocheat.95[itmp],rev(ocheat.05[itmp])), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
  lines(obs.ocheat.time, obs.ocheat.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.ocheat.time,rev(obs.ocheat.time)), c(obs.ocheat.norm+n.sig*obs.ocheat.err,rev(obs.ocheat.norm-n.sig*obs.ocheat.err)),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);

# >>> GSIC <<<
itmp=midx.gsic[1]:nrow(gsic.hind)
plot(mod.time[itmp], gsic.50[itmp], type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.01,.04), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Glaciers and\nsmall ice caps [m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' c')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time[itmp],rev(mod.time[itmp])), c(gsic.95[itmp],rev(gsic.05[itmp])), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
  lines(obs.gsic.time, obs.gsic.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.gsic.time,rev(obs.gsic.time)), c(obs.gsic.norm+n.sig*obs.gsic.err,rev(obs.gsic.norm-n.sig*obs.gsic.err)),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);

# >>> GIS <<<
plot(mod.time, gis.50, type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.003,.01), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Greenland Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' d')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(gis.95,rev(gis.05)), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
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

plot(mod.time, te.50, type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(1950,2016), ylim=c(-.04,.04), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Thermal expansion\n[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' e')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(te.95,rev(te.05)), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
	lines(x1971,y1971, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
  polygon(c(x1971,rev(x1971)), c(lo1971,rev(hi1971)), col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);
	lines(x1993,y1993, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
	polygon(c(x1993,rev(x1993)), c(lo1993,rev(hi1993)), col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);

# >>> TOTAL SLR <<<
plot(mod.time, slr.50, type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(1850,2016), ylim=c(-.3,.2), cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year', line=2.3, cex=.9);
  mtext(side=2, text='Total sea level [m]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' f')), line=.25, cex=.9, adj=0);
  polygon(c(mod.time,rev(mod.time)), c(slr.95,rev(slr.05)), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
  lines(obs.sl.time, obs.sl.norm, type='l', col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]), lwd=2);
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  polygon(c(obs.sl.time,rev(obs.sl.time)), c(obs.sl.norm+n.sig*obs.sl.err,rev(obs.sl.norm-n.sig*obs.sl.err)),
          col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.5), border=NA);

# >>> AIS PALEO, SMOOTHED <<<
ipaleo=which(t.paleo==-149999):which(t.paleo==1)
plot(t.paleo[ipaleo], ais.paleo.50[ipaleo], type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=2, xlab='',
     ylab='', xlim=c(t.paleo[ipaleo[1]],t.paleo[ipaleo[length(ipaleo)]]), ylim=c(-20,10),
		 cex.lab=1.2, cex.axis=1.2);
  mtext(side=1, text='Year [before present]', line=2.3, cex=.9);
  mtext(side=2, text='Antarctic Ice Sheet\n[m SLE]', line=2.3, cex=.9);
  mtext(side=3, text=expression(bold(' g')), line=.25, cex=.9, adj=0);
  polygon(c(t.paleo[ipaleo],rev(t.paleo[ipaleo])), c(ais.paleo.95[ipaleo],rev(ais.paleo.05[ipaleo])), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
	for (i in 1:3) {
  	polygon(c( date[c(obs.years[i]-1000,obs.years[i]+1000)], rev(date[c(obs.years[i]-1000,obs.years[i]+1000)]) ),
	 			  	c( c(windows[i,2],windows[i,2])                , rev(c(windows[i,1],windows[i,1]))                 ),
          	col=rgb(mycol[6,1],mycol[6,2],mycol[6,3],.7), border=NA);
	}
	i=4; points(date[obs.years[i]],mean(windows[i,]),pch=15,col=rgb(mycol[6,1],mycol[6,2],mycol[6,3]))
	lines(c(-1e6,1e6),c(0,0),type='l',lty=2,col='black');
  legend(-90000,29,c("5-95% range, model","2-sigma range, observations"),
         col=c(rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),rgb(mycol[6,1],mycol[6,2],mycol[6,3]),rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3])), lwd=2, bty='n', cex=1.2)

dev.off()

##==============================================================================
##==============================================================================






##==============================================================================
##==============================================================================
## FIGURE 2 -- SLR PROJECTIONS (DIFFERENT RCP SCENARIOS), VS MENGEL
##=========

ncdata <- nc_open(filename.brick)
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

##==============================================================================

## Get Mengel et al (2016) 5-95% interval for comparison
## NOTE: these are relative to 1986-2005 period!
ind.m16 = which(mod.time==1986):which(mod.time==2005)

slr.rcp26.m16 = c(.2799,.3938,0.5555)
slr.rcp45.m16 = c(.3708,.5290,0.7727)
slr.rcp85.m16 = c(.5741,.8455,1.3120)

te.rcp26.m16 = c(.0662,.1490,.2280)
te.rcp45.m16 = c(.0858,.1940,.3030)
te.rcp85.m16 = c(.1200,.2910,.4540)

gsic.rcp26.m16 = c(.0624,.0790,.1030)
gsic.rcp45.m16 = c(.0725,.0932,.1220)
gsic.rcp85.m16 = c(.0848,.1090,.1470)

gis.sid.rcp26.m16 = c(.0351,.0474,.0872)
gis.sid.rcp45.m16 = c(.0415,.0557,.1090)
gis.sid.rcp85.m16 = c(.0508,.0741,.1470)

gis.smb.rcp26.m16 = c(.0401,.0697,.1160)
gis.smb.rcp45.m16 = c(.0693,.1170,.2140)
gis.smb.rcp85.m16 = c(.1520,.2660,.5180)

gis.rcp26.m16 = c(gis.sid.rcp26.m16[2]+gis.smb.rcp26.m16[2]-sqrt(diff(gis.sid.rcp26.m16)^2 + diff(gis.smb.rcp26.m16)^2)[1],
									gis.sid.rcp26.m16[2]+gis.smb.rcp26.m16[2],
									gis.sid.rcp26.m16[2]+gis.smb.rcp26.m16[2]+sqrt(diff(gis.sid.rcp26.m16)^2 + diff(gis.smb.rcp26.m16)^2)[2])
gis.rcp45.m16 = c(gis.sid.rcp45.m16[2]+gis.smb.rcp45.m16[2]-sqrt(diff(gis.sid.rcp45.m16)^2 + diff(gis.smb.rcp45.m16)^2)[1],
									gis.sid.rcp45.m16[2]+gis.smb.rcp45.m16[2],
									gis.sid.rcp45.m16[2]+gis.smb.rcp45.m16[2]+sqrt(diff(gis.sid.rcp45.m16)^2 + diff(gis.smb.rcp45.m16)^2)[2])
gis.rcp85.m16 = c(gis.sid.rcp85.m16[2]+gis.smb.rcp85.m16[2]-sqrt(diff(gis.sid.rcp85.m16)^2 + diff(gis.smb.rcp85.m16)^2)[1],
									gis.sid.rcp85.m16[2]+gis.smb.rcp85.m16[2],
									gis.sid.rcp85.m16[2]+gis.smb.rcp85.m16[2]+sqrt(diff(gis.sid.rcp85.m16)^2 + diff(gis.smb.rcp85.m16)^2)[2])

ais.sid.rcp26.m16 = c(.0404,.0644,.0910)
ais.sid.rcp45.m16 = c(.0559,.0854,.1240)
ais.sid.rcp85.m16 = c(.0888,.1280,.1890)

ais.smb.rcp26.m16 = c(-.0263,-.0160,-.0079)
ais.smb.rcp45.m16 = c(-.0337,-.0203,-.00996)
ais.smb.rcp85.m16 = c(-.0483,-.0286,-.0138)

ais.rcp26.m16 = c(ais.sid.rcp26.m16[2]+ais.smb.rcp26.m16[2]-sqrt(diff(ais.sid.rcp26.m16)^2 + diff(ais.smb.rcp26.m16)^2)[1],
									ais.sid.rcp26.m16[2]+ais.smb.rcp26.m16[2],
									ais.sid.rcp26.m16[2]+ais.smb.rcp26.m16[2]+sqrt(diff(ais.sid.rcp26.m16)^2 + diff(ais.smb.rcp26.m16)^2)[2])
ais.rcp45.m16 = c(ais.sid.rcp45.m16[2]+ais.smb.rcp45.m16[2]-sqrt(diff(ais.sid.rcp45.m16)^2 + diff(ais.smb.rcp45.m16)^2)[1],
									ais.sid.rcp45.m16[2]+ais.smb.rcp45.m16[2],
									ais.sid.rcp45.m16[2]+ais.smb.rcp45.m16[2]+sqrt(diff(ais.sid.rcp45.m16)^2 + diff(ais.smb.rcp45.m16)^2)[2])
ais.rcp85.m16 = c(ais.sid.rcp85.m16[2]+ais.smb.rcp85.m16[2]-sqrt(diff(ais.sid.rcp85.m16)^2 + diff(ais.smb.rcp85.m16)^2)[1],
									ais.sid.rcp85.m16[2]+ais.smb.rcp85.m16[2],
									ais.sid.rcp85.m16[2]+ais.smb.rcp85.m16[2]+sqrt(diff(ais.sid.rcp85.m16)^2 + diff(ais.smb.rcp85.m16)^2)[2])

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


pdf(paste(plotdir,'projections_SLR_total_vs_Mengel.pdf',sep=''),width=5,height=3.5,colormodel='cmyk')
par(mfrow=c(1,1))
# UNIFORM RCP85
par(mai=c(.65,.65,.20,.4))
plot(t.proj[iproj],slr.rcp85.50[iproj],type='l',col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]),lwd=2, ann='',
		 xlim=c(2000,2111), ylim=c(0,1.7), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(2000,2100,by=20)); axis(2, seq(0,2,by=.25), lab=c('0','','0.5','','1','','1.5','','2'));
		 mtext(side=2, text='Total sea level [m]', line=2.2, cex=1);
     mtext(side=1, text='Year', line=2.2, cex=1);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp85.95[iproj],rev(slr.rcp85.05[iproj])),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA);
# + UNIFORM RCP45
	lines(t.proj[iproj],slr.rcp45.50[iproj],type='l',col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp45.95[iproj],rev(slr.rcp45.05[iproj])),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
# + UNIFORM RCP26
	lines(t.proj[iproj],slr.rcp26.50[iproj],type='l',col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp26.95[iproj],rev(slr.rcp26.05[iproj])),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA);
# + legend
  legend(t.proj[iproj[1]]+10,1.5,c("5-95% range,",
                                "RCP2.6",
																"RCP4.5",
																"RCP8.5"),
         lty=c(NA,1,1,1), lwd=3, col=c(NA,rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]),rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3])), bty='n', cex=1)

# + "Mengel" label
#	text(2110, 1.8, "Mengel et al, 2016", pos = 2, cex = 1, srt = 90)
	 mtext(side = 4, text = "Mengel et al, 2016", line = 0.5, cex=1.2)
# + Shaded box around Mengel
	polygon(c(t.proj[i2100]+0.5,t.proj[i2100]+11,t.proj[i2100]+11,t.proj[i2100]+.5),
					c(10,10,-10,-10), col=rgb(.5,.5,.5,.15), border=NA)

# + Mengel RCP26
  polygon(c(t.proj[i2100]+1,t.proj[i2100]+4,t.proj[i2100]+4,t.proj[i2100]+1),
					c(slr.rcp26.m16[3],slr.rcp26.m16[3],slr.rcp26.m16[1],slr.rcp26.m16[1]),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA)
  lines(c(t.proj[i2100]+1,t.proj[i2100]+4),c(slr.rcp26.m16[2],slr.rcp26.m16[2]),type='l',lwd=2,col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]));
# + Mengel RCP45
  polygon(c(t.proj[i2100]+4,t.proj[i2100]+7,t.proj[i2100]+7,t.proj[i2100]+4),
					c(slr.rcp45.m16[3],slr.rcp45.m16[3],slr.rcp45.m16[1],slr.rcp45.m16[1]),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA)
  lines(c(t.proj[i2100]+4,t.proj[i2100]+7),c(slr.rcp45.m16[2],slr.rcp45.m16[2]),type='l',lwd=2,col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]));
# + Mengel RCP85
  polygon(c(t.proj[i2100]+7,t.proj[i2100]+10,t.proj[i2100]+10,t.proj[i2100]+7),
					c(slr.rcp85.m16[3],slr.rcp85.m16[3],slr.rcp85.m16[1],slr.rcp85.m16[1]),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],0.5), border=NA)
  lines(c(t.proj[i2100]+7,t.proj[i2100]+10),c(slr.rcp85.m16[2],slr.rcp85.m16[2]),type='l',lwd=2,col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]));

dev.off()

##==============================================================================
##==============================================================================








##==============================================================================
##==============================================================================
## FIGURE 3 -- PROJECTIONS OF SLR WITH WAIS
##=========

## Get a WAIS simulation to find a worst-case scenario
source('../R/WAIS.R')
wais.out = WAIS(start.year=2016, t=t.proj)
wais.worst = wais.out[1]-wais.out
wais.out = WAIS(dSmx = 0.04731183,start.year = 2040, t=t.proj)
wais.dp16 = wais.out[1]-wais.out

##
## 5-95% CI of SLR projections -- including potential WAIS contribution
##

i2095=which(t.proj==2095)

pdf(paste(plotdir,'projections_SLR_total_DP16_wais.pdf',sep=''),width=5,height=5,colormodel='cmyk')
par(mfrow=c(1,1), mai=c(.7,.7,.2,.25))
# RCP85
plot(t.proj[iproj],slr.rcp85.50[iproj],type='l',col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]),lwd=2, xlab='', ylab='',
		 xlim=c(2000,2100), ylim=c(0,5), xaxt='n', yaxt='n', xaxs='i', yaxs='i', lty=1);
		 axis(1, seq(1960,2100,by=20)); axis(2, seq(0,5,by=.50), lab=c('0','','1','','2','','3','','4','','5'))
		 mtext(side=1, text='Year', line=2.1, cex=1.2);
		 mtext(side=2, text='Sea level [m]', line=2.1, cex=1.2);
	polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(slr.rcp85.95[iproj],rev(slr.rcp85.05[iproj])),
	        col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],0.5), border=NA);

lines(t.proj[iproj],slr.rcp85.50[iproj]+wais.worst[iproj], col=rgb(mycol[11,1],mycol[11,2],mycol[11,3]), lwd=2)
	polygon(c(t.proj[iproj],rev(t.proj[iproj])),
					c(slr.rcp85.95[iproj]+wais.worst[iproj],rev(slr.rcp85.05[iproj]+wais.worst[iproj])),
	        col=rgb(mycol[11,1],mycol[11,2],mycol[11,3],0.5), border=NA);

lines(t.proj[iproj],slr.rcp85.50[iproj]+wais.dp16[iproj], col=rgb(mycol[12,1],mycol[12,2],mycol[12,3]), lwd=2)
	polygon(c(t.proj[iproj],rev(t.proj[iproj])),
					c(slr.rcp85.95[iproj]+wais.dp16[iproj],rev(slr.rcp85.05[iproj]+wais.dp16[iproj])),
	        col=rgb(mycol[13,1],mycol[13,2],mycol[13,3],0.5), border=NA);

# + WAIS deep uncertainty block
  polygon(c(t.proj[iproj],rev(t.proj[iproj])),
					c(slr.rcp85.95[iproj]+wais.dp16[iproj],rev(slr.rcp85.05[iproj]+wais.worst[iproj])),
          col=rgb(.5,.5,.5,0.3), border=NA);
# + Deep uncertainty arrow
  arrows(2095, slr.rcp85.95[i2095]+wais.dp16[i2095]+.05, 2095, slr.rcp85.05[i2095]+wais.worst[i2095]-.05,
					length=0.15, angle=45, lty=1, lwd=6, code=3)
	text(2080.5,2.6,"deep \nuncertainty",cex=1.2,srt=0)
# + legend
  legend(t.proj[iproj[1]],4.9,c("5-95% range, RCP8.5",
																"+ WAIS, DP16",
																"+ WAIS, worst-case"),
    lty=c(1,1,1), lwd=3, col=c(rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]),rgb(mycol[12,1],mycol[12,2],mycol[12,3]),rgb(mycol[11,1],mycol[11,2],mycol[11,3])), bty='n')
dev.off()

##==============================================================================
##==============================================================================








##==============================================================================
##==============================================================================
## FIGURE 3 -- PROJECTIONS OF EACH COMPONENT OF SLR, VS MENGEL
##=========

##
## 5-95% CI of projections for each component of SLR
##

pdf(paste(plotdir,'projections_SLR_contributions-temp-ocheat.pdf',sep=''),width=6,height=7,colormodel='cmyk')
par(mfrow=c(3,2))
## >>> SURFACE TEMPERATURE <<<
# RCP85
par(mai=c(.4,.5,.1,.05))
plot(t.proj[iproj],temp.rcp85.50[iproj],type='l',col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]),lwd=2, xlab='', ylab='',
		 xlim=c(2000,2111), ylim=c(0,5), xaxt='n', yaxt='n', xaxs='i', yaxs='i', cex.lab=1.1);
		 axis(1, seq(1960,2100,by=20), cex.axis=1.1); axis(2, seq(0,5,by=.50), cex.axis=1.1)
		 mtext(side = 2, text = "Surface temperature [deg C]", line = 2.2, cex=0.8)
     mtext(side=3, text=expression(bold(' a')), line=-1.2, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(temp.rcp85.95[iproj],rev(temp.rcp85.05[iproj])),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA);
# + RCP26
lines(t.proj[iproj],temp.rcp26.50[iproj],type='l',col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(temp.rcp26.95[iproj],rev(temp.rcp26.05[iproj])),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA);
# + RCP45
lines(t.proj[iproj],temp.rcp45.50[iproj],type='l',col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(temp.rcp45.95[iproj],rev(temp.rcp45.05[iproj])),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
# + Legend
legend(t.proj[iproj[1]]-2,4.7,c("5-95% range, RCP2.6",
																	"5-95% range, RCP4.5",
																	"5-95% range, RCP8.5"),
         lty=c(1,1,1), lwd=3,
         col=c(rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]),rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3])),
         bty='n', cex=1.1)

## >>> OCEAN HEAT UPTAKE <<<
# RCP85
par(mai=c(.4,.5,.1,.05))
plot(t.proj[iproj],ocheat.rcp85.50[iproj],type='l',col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]),lwd=2, xlab='', ylab='',
		 xlim=c(2000,2111), ylim=c(0,400), xaxt='n', yaxt='n', xaxs='i', yaxs='i', cex.lab=1.1);
		 axis(1, seq(1960,2100,by=20),cex.axis=1.1); axis(2, seq(0,400,by=50),cex.axis=1.1)
		 mtext(side = 2, text = "Ocean heat uptake [10^22 J]", line = 2.2, cex=0.8)
     mtext(side=3, text=expression(bold(' b')), line=-1.2, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(ocheat.rcp85.95[iproj],rev(ocheat.rcp85.05[iproj])),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA);
# + RCP26
lines(t.proj[iproj],ocheat.rcp26.50[iproj],type='l',col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(ocheat.rcp26.95[iproj],rev(ocheat.rcp26.05[iproj])),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA);
# + RCP45
lines(t.proj[iproj],ocheat.rcp45.50[iproj],type='l',col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(ocheat.rcp45.95[iproj],rev(ocheat.rcp45.05[iproj])),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);

## >>> AIS SLR <<<
# RCP85
par(mai=c(.4,.5,.1,.05))
plot(t.proj[iproj],ais.rcp85.50[iproj],type='l',col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]),lwd=2, xlab='', ylab='',
		 xlim=c(2000,2111), ylim=c(0,1), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(1960,2100,by=20), cex.axis=1.1); axis(2, seq(0,1,by=.10), cex.axis=1.1)
		 mtext(side = 2, text = "Antarctic Ice Sheet [m SLE]", line = 2.2, cex=0.8)
     mtext(side=3, text=expression(bold(' c')), line=-1.2, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(ais.rcp85.95[iproj],rev(ais.rcp85.05[iproj])),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA);
# + RCP26
lines(t.proj[iproj],ais.rcp26.50[iproj],type='l',col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(ais.rcp26.95[iproj],rev(ais.rcp26.05[iproj])),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA);
# + RCP45
lines(t.proj[iproj],ais.rcp45.50[iproj],type='l',col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(ais.rcp45.95[iproj],rev(ais.rcp45.05[iproj])),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
# + "Mengel" label
	text(2109, .9, "Mengel et al, 2016", pos = 2, cex = 1.1, srt = 90)
# + Mengel RCP26
  polygon(c(t.proj[i2100]+1,t.proj[i2100]+4,t.proj[i2100]+4,t.proj[i2100]+1),
					c(ais.rcp26.m16[3],ais.rcp26.m16[3],ais.rcp26.m16[1],ais.rcp26.m16[1]),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA)
  lines(c(t.proj[i2100]+1,t.proj[i2100]+4),c(ais.rcp26.m16[2],ais.rcp26.m16[2]),type='l',lwd=2,col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]));
# + Mengel RCP45
  polygon(c(t.proj[i2100]+4,t.proj[i2100]+7,t.proj[i2100]+7,t.proj[i2100]+4),
					c(ais.rcp45.m16[3],ais.rcp45.m16[3],ais.rcp45.m16[1],ais.rcp45.m16[1]),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA)
  lines(c(t.proj[i2100]+4,t.proj[i2100]+7),c(ais.rcp45.m16[2],ais.rcp45.m16[2]),type='l',lwd=2,col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]));
# + Mengel RCP85
  polygon(c(t.proj[i2100]+7,t.proj[i2100]+10,t.proj[i2100]+10,t.proj[i2100]+7),
					c(ais.rcp85.m16[3],ais.rcp85.m16[3],ais.rcp85.m16[1],ais.rcp85.m16[1]),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA)
  lines(c(t.proj[i2100]+7,t.proj[i2100]+10),c(ais.rcp85.m16[2],ais.rcp85.m16[2]),type='l',lwd=2,col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]));
# + Shaded box around Mengel
	polygon(c(t.proj[i2100]+0.5,t.proj[i2100]+11,t.proj[i2100]+11,t.proj[i2100]+.5),
					c(10,10,-10,-10), col=rgb(.5,.5,.5,.15), border=NA)

## >>> GIS SLR <<<
# RCP85
par(mai=c(.4,.5,.1,.05))
plot(t.proj[iproj],gis.rcp85.50[iproj],type='l',col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]),lwd=2, xlab='', ylab='',
		 xlim=c(2000,2111), ylim=c(0,1), xaxt='n', yaxt='n', xaxs='i', yaxs='i', cex.lab=1.1);
		 axis(1, seq(1960,2100,by=20), cex.axis=1.1); axis(2, seq(0,1,by=.10), cex.axis=1.1)
		 mtext(side = 2, text = "Greenland Ice Sheet [m SLE]", line = 2.2, cex=0.8)
     mtext(side=3, text=expression(bold(' d')), line=-1.2, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gis.rcp85.95[iproj],rev(gis.rcp85.05[iproj])),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA);
# + RCP26
lines(t.proj[iproj],gis.rcp26.50[iproj],type='l',col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gis.rcp26.95[iproj],rev(gis.rcp26.05[iproj])),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA);
# + RCP45
lines(t.proj[iproj],gis.rcp45.50[iproj],type='l',col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gis.rcp45.95[iproj],rev(gis.rcp45.05[iproj])),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA);
# + Mengel RCP26
  polygon(c(t.proj[i2100]+1,t.proj[i2100]+4,t.proj[i2100]+4,t.proj[i2100]+1),
					c(gis.rcp26.m16[3],gis.rcp26.m16[3],gis.rcp26.m16[1],gis.rcp26.m16[1]),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA)
  lines(c(t.proj[i2100]+1,t.proj[i2100]+4),c(gis.rcp26.m16[2],gis.rcp26.m16[2]),type='l',lwd=2,col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]));
# + Mengel RCP45
  polygon(c(t.proj[i2100]+4,t.proj[i2100]+7,t.proj[i2100]+7,t.proj[i2100]+4),
					c(gis.rcp45.m16[3],gis.rcp45.m16[3],gis.rcp45.m16[1],gis.rcp45.m16[1]),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA)
  lines(c(t.proj[i2100]+4,t.proj[i2100]+7),c(gis.rcp45.m16[2],gis.rcp45.m16[2]),type='l',lwd=2,col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]));
# + Mengel RCP85
  polygon(c(t.proj[i2100]+7,t.proj[i2100]+10,t.proj[i2100]+10,t.proj[i2100]+7),
					c(gis.rcp85.m16[3],gis.rcp85.m16[3],gis.rcp85.m16[1],gis.rcp85.m16[1]),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA)
  lines(c(t.proj[i2100]+7,t.proj[i2100]+10),c(gis.rcp85.m16[2],gis.rcp85.m16[2]),type='l',lwd=2,col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]));
# + Shaded box around Mengel
	polygon(c(t.proj[i2100]+0.5,t.proj[i2100]+11,t.proj[i2100]+11,t.proj[i2100]+.5),
					c(10,10,-10,-10), col=rgb(.5,.5,.5,.15), border=NA)

## >>> GSIC SLR <<<
# RCP85
par(mai=c(.5,.5,.1,.05))
plot(t.proj[iproj],gsic.rcp85.50[iproj],type='l',col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]),lwd=2, xlab='', ylab='',
		 xlim=c(2000,2111), ylim=c(0,1), xaxt='n', yaxt='n', xaxs='i', yaxs='i', cex.lab=1.1);
		 axis(1, seq(1960,2100,by=20), cex.axis=1.1); axis(2, seq(0,1,by=.10), cex.axis=1.1)
     mtext(side = 1, text = "Year", line = 2.1, cex=0.8)
		 mtext(side = 2, text = "Glaciers, small ice caps [m SLE]", line = 2.2, cex=0.8)
     mtext(side=3, text=expression(bold(' e')), line=-1.2, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsic.rcp85.95[iproj],rev(gsic.rcp85.05[iproj])),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA);
# + RCP26
lines(t.proj[iproj],gsic.rcp26.50[iproj],type='l',col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsic.rcp26.95[iproj],rev(gsic.rcp26.05[iproj])),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA);
# + RCP45
lines(t.proj[iproj],gsic.rcp45.50[iproj],type='l',col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(gsic.rcp45.95[iproj],rev(gsic.rcp45.05[iproj])),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.35), border=NA);
# + Shaded box around Mengel
	polygon(c(t.proj[i2100]+0.5,t.proj[i2100]+11,t.proj[i2100]+11,t.proj[i2100]+.5),
					c(10,10,-10,-10), col=rgb(.5,.5,.5,.15), border=NA)
# + Mengel RCP26
  polygon(c(t.proj[i2100]+1,t.proj[i2100]+4,t.proj[i2100]+4,t.proj[i2100]+1),
					c(gsic.rcp26.m16[3],gsic.rcp26.m16[3],gsic.rcp26.m16[1],gsic.rcp26.m16[1]),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA)
  lines(c(t.proj[i2100]+1,t.proj[i2100]+4),c(gsic.rcp26.m16[2],gsic.rcp26.m16[2]),type='l',lwd=2,col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]));
# + Mengel RCP45
  polygon(c(t.proj[i2100]+4,t.proj[i2100]+7,t.proj[i2100]+7,t.proj[i2100]+4),
					c(gsic.rcp45.m16[3],gsic.rcp45.m16[3],gsic.rcp45.m16[1],gsic.rcp45.m16[1]),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA)
  lines(c(t.proj[i2100]+4,t.proj[i2100]+7),c(gsic.rcp45.m16[2],gsic.rcp45.m16[2]),type='l',lwd=2,col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]));
# + Mengel RCP85
  polygon(c(t.proj[i2100]+7,t.proj[i2100]+10,t.proj[i2100]+10,t.proj[i2100]+7),
					c(gsic.rcp85.m16[3],gsic.rcp85.m16[3],gsic.rcp85.m16[1],gsic.rcp85.m16[1]),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA)
  lines(c(t.proj[i2100]+7,t.proj[i2100]+10),c(gsic.rcp85.m16[2],gsic.rcp85.m16[2]),type='l',lwd=2,col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]));

## >>> TE SLR <<<
# RCP85
par(mai=c(.5,.5,.1,.05))
plot(t.proj[iproj],te.rcp85.50[iproj],type='l',col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]),lwd=2, xlab='', ylab='',
		 xlim=c(2000,2111), ylim=c(0,1), xaxt='n', yaxt='n', xaxs='i', yaxs='i');
		 axis(1, seq(1960,2100,by=20), cex.axis=1.1); axis(2, seq(0,1,by=.10), cex.axis=1.1)
     mtext(side = 1, text = "Year", line = 2.1, cex=0.8)
		 mtext(side = 2, text = "Thermal expansion [m SLE]", line = 2.2, cex=0.8)
     mtext(side=3, text=expression(bold(' f')), line=-1.2, cex=.9, adj=0);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(te.rcp85.95[iproj],rev(te.rcp85.05[iproj])),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA);
# + RCP26
lines(t.proj[iproj],te.rcp26.50[iproj],type='l',col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(te.rcp26.95[iproj],rev(te.rcp26.05[iproj])),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA);
# + RCP45
lines(t.proj[iproj],te.rcp45.50[iproj],type='l',col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),lwd=2);
  polygon(c(t.proj[iproj],rev(t.proj[iproj])), c(te.rcp45.95[iproj],rev(te.rcp45.05[iproj])),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.35), border=NA);
# + Mengel RCP26
  polygon(c(t.proj[i2100]+1,t.proj[i2100]+4,t.proj[i2100]+4,t.proj[i2100]+1),
					c(te.rcp26.m16[3],te.rcp26.m16[3],te.rcp26.m16[1],te.rcp26.m16[1]),
          col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3],.5), border=NA)
  lines(c(t.proj[i2100]+1,t.proj[i2100]+4),c(te.rcp26.m16[2],te.rcp26.m16[2]),type='l',lwd=2,col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]));
# + Mengel RCP45
  polygon(c(t.proj[i2100]+4,t.proj[i2100]+7,t.proj[i2100]+7,t.proj[i2100]+4),
					c(te.rcp45.m16[3],te.rcp45.m16[3],te.rcp45.m16[1],te.rcp45.m16[1]),
          col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3],.5), border=NA)
  lines(c(t.proj[i2100]+4,t.proj[i2100]+7),c(te.rcp45.m16[2],te.rcp45.m16[2]),type='l',lwd=2,col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]));
# + Mengel RCP85
  polygon(c(t.proj[i2100]+7,t.proj[i2100]+10,t.proj[i2100]+10,t.proj[i2100]+7),
					c(te.rcp85.m16[3],te.rcp85.m16[3],te.rcp85.m16[1],te.rcp85.m16[1]),
          col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3],.5), border=NA)
  lines(c(t.proj[i2100]+7,t.proj[i2100]+10),c(te.rcp85.m16[2],te.rcp85.m16[2]),type='l',lwd=2,col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]));
# + Shaded box around Mengel
	polygon(c(t.proj[i2100]+0.5,t.proj[i2100]+11,t.proj[i2100]+11,t.proj[i2100]+.5),
					c(10,10,-10,-10), col=rgb(.5,.5,.5,.15), border=NA)

dev.off()

##==============================================================================
##==============================================================================







##==============================================================================
##==============================================================================
## FIGURE 4 -- PROJECTIONS, PROBABILITY DISTRIBUTIONS, SURVIVAL FUNCTIONS
##=========

i2100 = which(t.proj==2100)

slr2100.rcp26 = slr.rcp26[i2100,]
slr2100.rcp45 = slr.rcp45[i2100,]
slr2100.rcp85 = slr.rcp85[i2100,]

pdf.slr2100.rcp26 = density(slr2100.rcp26,from=0,to=4,na.rm=TRUE)
pdf.slr2100.rcp45 = density(slr2100.rcp45,from=0,to=4,na.rm=TRUE)
pdf.slr2100.rcp85 = density(slr2100.rcp85,from=0,to=4,na.rm=TRUE)

cdf.slr2100.rcp26 = rep(0,length(pdf.slr2100.rcp26$x))
cdf.slr2100.rcp45 = rep(0,length(pdf.slr2100.rcp45$x))
cdf.slr2100.rcp85 = rep(0,length(pdf.slr2100.rcp85$x))

dx=median(diff(pdf.slr2100.rcp26$x))
x=pdf.slr2100.rcp26$x

for (i in 2:length(x)){
  cdf.slr2100.rcp26[i] = sum(c(cdf.slr2100.rcp26[i-1],dx*pdf.slr2100.rcp26$y[i]))
  cdf.slr2100.rcp45[i] = sum(c(cdf.slr2100.rcp45[i-1],dx*pdf.slr2100.rcp45$y[i]))
  cdf.slr2100.rcp85[i] = sum(c(cdf.slr2100.rcp85[i-1],dx*pdf.slr2100.rcp85$y[i]))
}

cdf.slr2100.rcp26 = cdf.slr2100.rcp26/cdf.slr2100.rcp26[length(cdf.slr2100.rcp26)]
cdf.slr2100.rcp45 = cdf.slr2100.rcp45/cdf.slr2100.rcp45[length(cdf.slr2100.rcp45)]
cdf.slr2100.rcp85 = cdf.slr2100.rcp85/cdf.slr2100.rcp85[length(cdf.slr2100.rcp85)]

sur.slr2100.rcp26 = 1-cdf.slr2100.rcp26
sur.slr2100.rcp45 = 1-cdf.slr2100.rcp45
sur.slr2100.rcp85 = 1-cdf.slr2100.rcp85

c(q005.rcp26,q05.rcp26,q95.rcp26,q995.rcp26) := quantile(slr2100.rcp26,c(.005,.05,.95,.995),na.rm=TRUE)
c(q005.rcp45,q05.rcp45,q95.rcp45,q995.rcp45) := quantile(slr2100.rcp45,c(.005,.05,.95,.995),na.rm=TRUE)
c(q005.rcp85,q05.rcp85,q95.rcp85,q995.rcp85) := quantile(slr2100.rcp85,c(.005,.05,.95,.995),na.rm=TRUE)


pdf(paste(plotdir,'distributions_SLR2100_pdf+sf.pdf',sep=''),width=5,height=7,colormodel='cmyk')
par(mfrow=c(2,1), mai=c(.85,.74,.1,.15))

plot(x,pdf.slr2100.rcp85$y, type='l', xlim=c(0,2.5), ylim=c(-1.25,4.7), lty=1,
     col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]), lwd=1.5, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i',axes=FALSE);
  axis(1,seq(0,3,0.5),lab=c("0","0.5","1","1.5","2","2.5","3"))
  u <- par("usr")
  arrows(0, u[3],0, u[4], code = 2, xpd = TRUE)
  mtext('Probability density', side=2, line=1.3);
  mtext('Projected sea level in 2100\nrelative to 1986-2005 average [m]', side=1, line=3);
  mtext(side=3, text=expression(bold('   a')), line=-1, cex=.9, adj=0);
  lines(x,pdf.slr2100.rcp45$y, type='l', col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=1.5);
  lines(x,pdf.slr2100.rcp26$y, type='l', col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]), lwd=1.5);

	polygon(c(q95.rcp26,q95.rcp26,q05.rcp26,q05.rcp26), c(-0.07,-0.15,-0.15,-0.07), col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]), border=NA)
	polygon(c(q95.rcp45,q95.rcp45,q05.rcp45,q05.rcp45), c(-0.20,-0.28,-0.28,-0.20), col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), border=NA)
	polygon(c(q95.rcp85,q95.rcp85,q05.rcp85,q05.rcp85), c(-0.33,-0.41,-0.41,-0.33), col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]), border=NA)

# + Mengel RCP26 5-95% CI
  arrows(slr.rcp26.m16[1], -0.7, slr.rcp26.m16[3], -0.7, length=0.05, angle=90, lty=1, lwd=2, code=3, col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]))
# + Mengel RCP45 5-95% CI
  arrows(slr.rcp45.m16[1], -0.85, slr.rcp45.m16[3], -0.85, length=0.05, angle=90, lty=1, lwd=2, code=3, col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]))
# + Mengel RCP26 5-95% CI
  arrows(slr.rcp85.m16[1], -1, slr.rcp85.m16[3], -1, length=0.05, angle=90, lty=1, lwd=2, code=3, col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]))

# + Box around Mengel's estimates
	polygon(c(-100,-100,100,100), c(-.5,-100,-100,-.5), col=rgb(.5,.5,.5,.15), border=NA)
# + "Mengel" label
	text(2.45, -1, "Mengel et al, 2016", pos = 2, cex = 1, srt = 0)
# + "This study" label
	text(2.45, -0.3, "This study", pos = 2, cex = 1, srt = 0)

  legend(1.1,4.7,c("RCP2.6","RCP4.5","RCP8.5"),
        lty=c(1,1,1), lwd=2, col=c(rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]),rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]),rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3])),
        bty='n')

plot(x,log10(sur.slr2100.rcp85),type='l', xlim=c(0,2.5), ylim=c(-3.5,0), lty=1,
     col=rgb(mycol[c85,1],mycol[c85,2],mycol[c85,3]), lwd=1.5, xlab='', ylab='', xaxt='n', yaxt='n', xaxs='i', yaxs='i');
  axis(1,seq(0,3,0.5),lab=c("0","0.5","1","1.5","2","2.5","3"))
  mtext('Survival function [1-CDF]', side = 2, line=2.6);
  mtext('Projected sea level in 2100\nrelative to 1986-2005 average [m]', side=1, line=3);
  mtext(side=3, text=expression(bold('   b')), line=-1.1, cex=.9, adj=0);

	axis(2, at=seq(-4,0), label=parse(text=paste("10^", seq(-4,0), sep="")), las=1)

  lines(x,log10(sur.slr2100.rcp45),type='l',lty=1, col=rgb(mycol[c45,1],mycol[c45,2],mycol[c45,3]), lwd=1.5);
  lines(x,log10(sur.slr2100.rcp26),type='l',lty=1, col=rgb(mycol[c26,1],mycol[c26,2],mycol[c26,3]), lwd=1.5);

	lines(c(-4,4),c(-2,-2),lty=2,col='black'); text(0.35,-1.85,"1:100 level");
	lines(c(-4,4),c(-3,-3),lty=2,col='black'); text(0.35,-2.85,"1:1000 level");

# + Pfeffer (2008) upper bound of 2 m
  lines(c(2,2),c(-4,-0.6),lty=2,lwd=2,col='black'); text(1.9,-0.3,"Pfeffer et al. (2008):\n 2 m rise upper bound")

dev.off()

##==============================================================================
##==============================================================================








##==============================================================================
##==============================================================================
## SI 1 -- FIGURE 5 -- CALIBRATED PARAMETER MARGINAL DISTRIBUTIONS
##=================

dat.uni = read.csv(filename.parameters)
parameters.uni = dat.uni[1:(nrow(dat.uni)-1), ]
parnames = colnames(parameters.uni)

## Fit PDFs to the parameter distributions
np = ncol(parameters.uni)
bound.lower = rep(0, np)
bound.upper = rep(0, np)
for (pp in 1:np) {
	bound.lower[pp] = min(parameters.uni[,pp]) - 0.0*(max(parameters.uni[,pp])-min(parameters.uni[,pp]))
  bound.upper[pp] = max(parameters.uni[,pp]) + 0.0*(max(parameters.uni[,pp])-min(parameters.uni[,pp]))
}

pdf.uni=vector('list',np)
n.node=100
for (pp in 1:np){
  tmp = density(parameters.uni[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.uni[[pp]] = tmp; names(pdf.uni)[pp]=colnames(parameters.uni)[pp]
}

## Swap sigma.simple and (rho.gsic, sigma.gsic), so the label doesn't extend
## beyond the frame
plotorder=1:np
ind.sigma.simple = match('sigma.simple',parnames)
ind.sigma.gsic = match('sigma.gsic',parnames)
tmp=c(ind.sigma.gsic,ind.sigma.gsic+1)
plotorder[ind.sigma.gsic]=ind.sigma.simple
plotorder[(ind.sigma.gsic+1):ind.sigma.simple]=tmp

## Make the parameter names prettier, and use Greek letters/subscripts
parnames.fancy=parnames
parnames.fancy[match('kappa.doeclim',parnames)]=expression(kappa[DOECLIM])
parnames.fancy[match('alpha.doeclim',parnames)]=expression(alpha[DOECLIM])
parnames.fancy[match('T0',parnames)]=expression(T[0])
parnames.fancy[match('H0',parnames)]=expression(H[0])
parnames.fancy[match('beta0',parnames)]=expression(beta[0])
parnames.fancy[match('V0.gsic',parnames)]=expression(V[0,GSIC])
parnames.fancy[match('Gs0',parnames)]=expression(G[s,0])
parnames.fancy[match('a.te',parnames)]=expression(a[TE])
parnames.fancy[match('b.te',parnames)]=expression(b[TE])
parnames.fancy[match('invtau.te',parnames)]=expression(1/tau[TE])
parnames.fancy[match('TE0',parnames)]=expression(TE[0])
parnames.fancy[match('a.simple',parnames)]=expression(a[SIMPLE])
parnames.fancy[match('b.simple',parnames)]=expression(b[SIMPLE])
parnames.fancy[match('alpha.simple',parnames)]=expression(alpha[SIMPLE])
parnames.fancy[match('beta.simple',parnames)]=expression(beta[SIMPLE])
parnames.fancy[match('V0',parnames)]=expression(V[0,SIMPLE])
parnames.fancy[match('sigma.T',parnames)]=expression(sigma[T])
parnames.fancy[match('sigma.H',parnames)]=expression(sigma[H])
parnames.fancy[match('rho.T',parnames)]=expression(rho[T])
parnames.fancy[match('rho.H',parnames)]=expression(rho[H])
parnames.fancy[match('sigma.simple',parnames)]=expression(sigma[SIMPLE])
parnames.fancy[match('sigma.gsic',parnames)]=expression(sigma[GSIC])
parnames.fancy[match('rho.gsic',parnames)]=expression(rho[GSIC])
parnames.fancy[match('anto.a',parnames)]=expression(a[ANTO])
parnames.fancy[match('anto.b',parnames)]=expression(b[ANTO])
parnames.fancy[match('gamma',parnames)]=expression(gamma)
parnames.fancy[match('alpha.dais',parnames)]=expression(alpha[DAIS])
parnames.fancy[match('mu',parnames)]=expression(mu)
parnames.fancy[match('nu',parnames)]=expression(nu)
parnames.fancy[match('P0',parnames)]=expression(P[0])
parnames.fancy[match('kappa.dais',parnames)]=expression(kappa[DAIS])
parnames.fancy[match('f0',parnames)]=expression(f[0])
parnames.fancy[match('h0',parnames)]=expression(h[0])
parnames.fancy[match('b0',parnames)]=expression(b[0])
parnames.fancy[match('var.dais',parnames)]=expression(Var[DAIS])


pdf(paste(plotdir,'parameters_marginal_distributions.pdf',sep=''),width=6,height=9,colormodel='cmyk')

par(mfrow=c(8,5), mai=c(.45,.4,.1,.1))
for (pp in plotorder) {
	if(pp!=19) {plot(pdf.uni[[pp]]$x, pdf.uni[[pp]]$y, type='l', main="",
									xlab='', ylab='', cex.axis=1.2, cex.lab=1.2, yaxt="n")
                  mtext(parnames.fancy[pp],side=1,line=2.4)
	} else {plot(pdf.uni[[pp]]$x, pdf.uni[[pp]]$y, type='l', main="",
									xlab='', ylab='', cex.axis=1.2, cex.lab=1.2, yaxt="n")
                  mtext(parnames.fancy[pp],side=1,line=2.4)
					mtext('Probability density', side = 2, line = -1.5, outer = TRUE, at = NA, cex = 1.2)
	}
}
par(mai=c(.1,.1,.1,.1)); plot(1, type="n", axes=F, xlab="", ylab="");
#  legend(0.35,1.4,c("Uniform",'  priors',"Gamma",'  priors'),
#         col=c('black',NA,'red',NA), lty=c(1,NA,2,NA), bty='n', cex=1.2, lwd=1.3)
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
luse.waisdi   = FALSE		# West Antartic ice sheet emulator
luse.brick = cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais, luse.waisdi)
source('../calibration/BRICK_parameterSetup.R') # get bound.lower/upper
parnames.bounds=parnames

dat.uni = read.csv(filename.parameters)
parameters.uni = dat.uni[1:(nrow(dat.uni)-1), ]
parnames = colnames(parameters.uni)

rho.simple.fixed = signif(as.numeric(read.csv(filename.rho_simple_fixed)),digits=4)

parameters=parameters.uni

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

tmp = rbind( ci90[c(1:5,19:22,6:9,23:24,10:13,14:18,25,26:39),])
lb = rep(0,nrow(tmp)); ub = lb
for (i in 1:nrow(tmp)){
  itmp = match(parnames.bounds[i],tmp[,1])
  lb[itmp] = bound.lower[i]; ub[itmp] = bound.upper[i]
}
tmp = cbind(tmp,lb,ub)
ci90.rearr = rbind( tmp[1:match('sigma.simple',tmp[,1]),],
                    c("rho.simple",rho.simple.fixed,rho.simple.fixed,rho.simple.fixed,rho.simple.fixed,rho.simple.fixed),
                    tmp[(match('sigma.simple',tmp[,1])+1):39,])

models = c('DOECLIM','DOECLIM','DOECLIM','DOECLIM','DOECLIM','DOECLIM','DOECLIM','DOECLIM','DOECLIM',
           'GSIC','GSIC','GSIC','GSIC','GSIC','GSIC',
           'TE','TE','TE','TE',
           'SIMPLE-GIS','SIMPLE-GIS','SIMPLE-GIS','SIMPLE-GIS','SIMPLE-GIS','SIMPLE-GIS','SIMPLE-GIS',
           'AntOc','AntOc',
           'DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS','DAIS')
units = c('deg C','cm^2 s^{-1}','-','deg C','10^{22} J','deg C','10^{22} J','-','-',
          'm yr^{-1} deg C^{-1}','m SLE','-','m SLE','m SLE','-',
          'm deg C^{-1}','m SLE','yr^{-1}','m SLE',
          'm SLE deg C^{-1}','m SLE','yr^{-1} deg C^{-1}','yr^{-1}','m SLE','m SLE','-',
          'deg C deg C^{-1}','deg C',
          '-','-','m^{1/2}','m^{-1/2} yr^{-1/2}','m yr^{-1}','deg C^{-1}','m yr^{-1}','m','m deg C^{-1}','m','-','m^2 SLE')
output.table = data.frame( cbind(ci90.rearr[,1], models, ci90.rearr[,2:6], units))
colnames(output.table) = c('Parameter','Model','Median','5% quantile','95% quantile','Lower bound','Upper bound','Units')

## Write calibrated parameters to a csv file
today=Sys.Date(); today=format(today,format="%d%b%Y")
filename=paste(plotdir,'BRICK-robustSLR_ExtendedDataTable1_',today,'.csv', sep="")
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
row.var.dais = paste('$\\sigma^2_{DAIS}$ & DAIS & ',signif(p50[match("var.dais",parnames)],n.dig),'&',signif(p05[match("var.dais",parnames)],n.dig),'&',signif(p95[match("var.dais",parnames)],n.dig),'&',signif(bound.lower[match("var.dais",parnames.bounds)],n.dig),'&',signif(bound.upper[match("var.dais",parnames.bounds)],n.dig),' (m$^2$ SLE)\\\\',sep='')

table=data.frame(rbind(cat(row.S),cat(row.kappa.doeclim),cat(row.alpha.doeclim),cat(row.T0),cat(row.H0),cat(row.sigma.T),cat(row.sigma.H),cat(row.rho.T),cat(row.rho.H)),
											cat(row.beta0),cat(row.V0.gsic),cat(row.n),cat(row.Gs0),cat(row.sigma.gsic),cat(row.rho.gsic),
											cat(row.a.te),cat(row.b.te),cat(row.invtau.te),cat(row.TE0),
											cat(row.a.simple),cat(row.b.simple),cat(row.alpha.simple),cat(row.beta.simple),cat(row.V0),cat(row.sigma.simple),cat(row.rho.simple),
											cat(row.anto.a),cat(row.anto.b),
											cat(row.gamma),cat(row.alpha.dais),cat(row.mu),cat(row.nu),cat(row.P0),cat(row.kappa.dais),cat(row.f0),cat(row.h0),cat(row.c),cat(row.b0),cat(row.slope),cat(row.var.dais),
											row.names=NULL)
}
##==============================================================================
##==============================================================================







##==============================================================================
## End
##==============================================================================
