##==============================================================================
##
## Pipeline for processing DAIS calibration results and BRICK rest-of-model
## calibration results.
##
##  Required input:
##    filename.DAIScalibration  DAIS calibration parameter posterior draws
##    filename.BRICKcalibration BRICK rest-of-model calibration parameter posterior draws
##
##  Output:
##    BRICK-robustSLR_postcalibratedParameters_[datestamp].csv
##                            post-calibrated parameters file
##		BRICK-robustSLR_physical_[datestamp].nc
##														netCDF4 file with the BRICK physical model output
##														(should be as many runs as post-calibrated parameters
##														on the CSV file)
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

rm(list=ls())


t.beg = proc.time()

##==============================================================================
##==============================================================================
## Define the files you want to process

filename.BRICKcalibration = "../output_calibration/BRICK_calibratedParameters_12Aug2016.nc"
#filename.BRICKcalibration = "../output_calibration/BRICK_calibratedParameters_12Aug2016.csv"
filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_06Sep2016.csv"

appen=''		## Append file name? In case you process multiple files in one day
today=Sys.Date(); today=format(today,format="%d%b%Y")

## Define the files you want to read/create
filename.DAIScalibration = "../output_calibration/DAIS_calibratedParameters_11Aug2016.nc"
#filename.DAIScalibration = "../output_calibration/DAIS_calibratedParameters_11Aug2016.csv"
filename.parameters = paste('../output_calibration/BRICK-robustSLR_postcalibratedParameters_',today,appen,'.csv', sep="")
filename.brickout = paste('../output_model/BRICK-robustSLR_physical_',today,appen,'.nc',sep="")

n.ensemble = 50000

## Post-calibrate sea-level rise to n.sigma-wide window (n.sigma>0), prc.close% of model
## within an abs(n.sigma) window (n.sigma<0), or nothing (n.sigma==0)?
n.sigma = 4
prc.close = .90

##==============================================================================
##==============================================================================
# simulate stationary AR(1) process
ar1.sim = function(N,rho1,sigma) {
	x = rep(NA,N)
	for(i in 2:N)
		if(length(sigma)>1) {
			x[1] = sigma[1]/sqrt(1-rho1^2)
			x[i] = rho1*x[i-1] + rnorm(1,sd=sigma[i])
		} else {
			x[1] = sigma/sqrt(1-rho1^2)
			x[i] = rho1*x[i-1] + rnorm(1,sd=sigma)
		}
	return(x)
}
##==============================================================================
##==============================================================================




##==============================================================================
##==============================================================================
## Combine calibrated parameters from DAIS-fast dynamics and BRICK-RoM

## Make posterior parameter draws to run an ensemble of simulations and make
## projections of SLR (include DAIS, from calibrated parameters)
## Note: DAIS is part of BRICK. that the first parameter set is BRICK "rest-of-
## model", just that part that is not DAIS. DAIS is calibrated based on paleo
## data, so dealt with separately.

## First, read in the calibrated model parameters
print(paste('Reading calibrated non-DAIS model parameters...'))
	# csv version - antiquated
#dat.brick = read.csv(filename.BRICKcalibration)
#parameters.brick = dat.brick[1:(nrow(dat.brick)-1), ]
#bandwidths.brick = dat.brick[nrow(dat.brick)      , ]
#parnames.brick = colnames(parameters.brick)
	# netcdf version - good!
ncdata <- nc_open(filename.BRICKcalibration)
parameters.brick = ncvar_get(ncdata, 'BRICK_parameters')
parnames.brick = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters.brick = t(parameters.brick)
colnames(parameters.brick) = parnames.brick
print(paste('  ... read ',nrow(parameters.brick),' model parameters'))

#t.beg = proc.time()

## Second, read in the DAIS calibrated parameters
print(paste('Reading calibrated DAIS model parameters...'))
	# csv version - antiquated
#dat.dais = read.csv(filename.DAIScalibration)
#parameters.dais = dat.dais[1:(nrow(dat.dais)-1),]
#bandwidths.dais = dat.dais[nrow(dat.dais)      ,]
#parnames.dais   = colnames(parameters.dais)
	# netcdf version - good!
ncdata <- nc_open(filename.DAIScalibration)
parameters.dais = ncvar_get(ncdata, 'DAIS_parameters')
parnames.dais = ncvar_get(ncdata, 'parnames')
nc_close(ncdata)
parameters.dais = t(parameters.dais)
colnames(parameters.dais) = parnames.dais
print(paste('  ... read ',nrow(parameters.dais),' DAIS model parameters'))

## Define number of ensemble members. Check to make sure you are not asking for
## more ensemble members than there are parameters available

#n.ensemble=4000 ## set earlier, with the file names
n.ensemble = min( n.ensemble, nrow(parameters.brick), nrow(parameters.dais))

## Draw parameters for calibrated models

## Use the "rnorm" versions if you want to draw from the multivariate kernel
## density estimates. This is unnecessary when using ~100,000s of posterior
## parameters from MCMC, however.
parameters.ensemble.brick = mat.or.vec(n.ensemble , ncol(parameters.brick) )
ind.ensemble = sample( seq(1,nrow(parameters.brick)), size=n.ensemble, replace=FALSE)
for (p in 1:ncol(parameters.brick)){
	for (i in 1:n.ensemble){
	  #parameters.ensemble.brick[i,p] = rnorm( 1, mean=parameters.brick[ind.ensemble[i],p], sd=bandwidths.brick[1,p])
		parameters.ensemble.brick[i,p] = parameters.brick[ind.ensemble[i],p]
	}
}

## Draw parameters for calibrated DAIS model
parameters.ensemble.dais = mat.or.vec(n.ensemble , ncol(parameters.dais) )
ind.dais=sample( seq(1,nrow(parameters.dais)), size=n.ensemble, replace=FALSE)
for (p in 1:ncol(parameters.dais)){
	for (i in 1:n.ensemble){
	  #parameters.ensemble.dais[i,p] = rnorm( 1, mean=parameters.dais[ind.dais[i],p], sd=bandwidths.dais[1,p])
		parameters.ensemble.dais[i,p] = parameters.dais[ind.dais[i],p]
	}
}

## Set up the parameter indices within a new "parnames", which includes the
## DAIS parameters; original calibrated model parnames was saved (above) in case
## you need it later
parameters = cbind(parameters.ensemble.brick, parameters.ensemble.dais)
parnames = c(parnames.brick, parnames.dais)
rownames(parameters)=NULL

## Set up the model for hindcasts
l.project = FALSE
begyear = 1850
endyear = 2009
mod.time= begyear:endyear
begyear.norm = 1986
endyear.norm = 2005
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time=length(mod.time)

## Get the forcing data
forcing = read.csv( '../data/forcing_hindcast.csv', header=TRUE )

## Source some useful functions for manipulating data
source('../R/forcing_total.R')					# function to add up the total forcing
source('../R/compute_indices.R')				# function to determine the model and
																				# data indices for comparisons

## Source the model(s) and data
source('../fortran/R/doeclimF.R')  			# the DOECLIM model - resets the mod.time
source('../fortran/R/daisantoF.R')   		# DAIS (Antarctic Ice Sheet) model, with 'anto'
																				# for translating Tg (surface temperature anomaly)
																				# to Ta (Antarctic temperature reduced to sea level)
																				# and Toc (Antarctic/high-latitude ocean temperature)
source('../fortran/R/GSIC_magiccF.R') 	# the GSIC model
source('../fortran/R/brick_te_F.R')   	# TE (thermal expansion) model
source('../fortran/R/simpleF.R')      	# GIS (Greenland Ice Sheet) model

## Read the data sets for hindcast comparisons
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

luse.doeclim  = TRUE    # diffusion-ocean-energy balance climate model
luse.gsic     = TRUE    # glaciers and small ice caps contribution to SLR
luse.te       = TRUE    # thermal expansion contribution to SLR
luse.simple   = TRUE    # Greenland ice sheet model
luse.dais     = TRUE    # Antarctic ice sheet model
luse.brick = cbind(luse.doeclim, luse.gsic, luse.te, luse.simple, luse.dais)

source('../R/BRICK_coupledModel.R')

## Initialize matrix to store model ensemble output
brick.out = vector("list", n.ensemble)

## Initialize flag for bad runs (possible crashes include instability in SIMPLE
## or perturbed parameters being outside of plausible ranges)
badruns = rep(0, n.ensemble)

## Run the sample, and enjoy a nice progress bar
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {

	brick.out[[i]] = brick_model(	parameters.in			= as.numeric(parameters[i,]),
																parnames.in				= parnames,
																forcing.in				= forcing,
																l.project					= l.project,
																slope.Ta2Tg.in		= slope.Ta2Tg,
																intercept.Ta2Tg.in= intercept.Ta2Tg,
																mod.time					= mod.time,
																ind.norm.data 		= ind.norm.data,
																ind.norm.sl 			= ind.norm,
																luse.brick 				= luse.brick
																)

  setTxtProgressBar(pb, i)
}
close(pb)

## Before post-calibration, need to add the modeled statistical noise back in.
## Only using sea-level rise data, so only need to modify GSIC, GIS.
## using the statistical parameters for AR1, AR1 and Gaussian noise, respecively
## Do not do for AIS, because var.dais was fit to paleo data-model mismatch, not
## representative of the current era.

## Read rho.simple from file
rho.simple.fixed = as.numeric(read.csv(filename.rho_simple_fixed))

## Gather the fields for each simulation (easy referencing for plotting and
## analysis)
slr.out    = mat.or.vec(n.ensemble,length(mod.time))
temp.out   = mat.or.vec(n.ensemble,length(mod.time))
ocheat.out = mat.or.vec(n.ensemble,length(mod.time))
gsic.out   = mat.or.vec(n.ensemble,length(mod.time))
gis.out    = mat.or.vec(n.ensemble,length(mod.time))
ais.out    = mat.or.vec(n.ensemble,length(mod.time))
te.out     = mat.or.vec(n.ensemble,length(mod.time))

## Will also normalize the output to "ind.norm" (1961-1990? 1986-2005? (Mengel, IPCC))
slr.out.norm    = slr.out
temp.out.norm   = temp.out
ocheat.out.norm = ocheat.out
gsic.out.norm   = gsic.out
gis.out.norm    = gis.out
ais.out.norm		= ais.out
te.out.norm		  = te.out

## And add statistical noise
slr.norm.stat    = slr.out
temp.norm.stat   = temp.out
ocheat.norm.stat = ocheat.out
gsic.norm.stat   = gsic.out
gis.norm.stat    = gis.out
ais.norm.stat    = ais.out
te.norm.stat	   = te.out

## Go through each simulation and collect, normalize and add modeled error to
## the fields

pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {

	T0=parameters[i,match("T0",parnames)]
	H0=parameters[i,match("H0",parnames)]

  # set the results (note that we already have slr.out)
  temp.out[i,]   = brick.out[[i]]$doeclim.out$temp + T0
  ocheat.out[i,] = brick.out[[i]]$doeclim.out$ocheat + H0
  gsic.out[i,]   = brick.out[[i]]$gsic.out
  gis.out[i,]    = brick.out[[i]]$simple.out$sle.gis
	ais.out[i,]    = brick.out[[i]]$dais.out
	te.out[i,]     = brick.out[[i]]$te.out

  # Normalize the output to "ind.norm.data"
  temp.out.norm[i,]   = temp.out[i,]  -mean(temp.out[i,1:20])
  ocheat.out.norm[i,] = ocheat.out[i,]#-mean(ocheat.out[i,ind.norm.data[which(ind.norm.data[,1]=='ocheat'),2]:ind.norm.data[which(ind.norm.data[,1]=='ocheat'),3]])
  gsic.out.norm[i,]   = gsic.out[i,]  -mean(gsic.out[i,ind.norm.data[which(ind.norm.data[,1]=='gsic'),2]:ind.norm.data[which(ind.norm.data[,1]=='gsic'),3]])
  gis.out.norm[i,]    = gis.out[i,]   -mean(gis.out[i,ind.norm.data[which(ind.norm.data[,1]=='gis'),2]:ind.norm.data[which(ind.norm.data[,1]=='gis'),3]])
  ais.out.norm[i,]    = ais.out[i,]   -mean(ais.out[i,ind.norm.data[which(ind.norm.data[,1]=='ais'),2]:ind.norm.data[which(ind.norm.data[,1]=='ais'),3]])
	te.out.norm[i,]     = te.out[i,]    -mean(te.out[i,ind.norm.data[which(ind.norm.data[,1]=='te'),2]:ind.norm.data[which(ind.norm.data[,1]=='te'),3]])

	# Add the statistcal model for AR1 (or otherwise) noise
	sigma.T     =parameters[i,match("sigma.T"     ,parnames)]
	rho.T       =parameters[i,match("rho.T"       ,parnames)]
	sigma.H     =parameters[i,match("sigma.H"     ,parnames)]
	rho.H       =parameters[i,match("rho.H"       ,parnames)]
	sigma.gsic  =parameters[i,match("sigma.gsic"  ,parnames)]
	rho.gsic    =parameters[i,match("rho.gsic"    ,parnames)]
  sigma.simple=parameters[i,match("sigma.simple",parnames)]
  rho.simple  =parameters[i,match("rho.simple"  ,parnames)]
	var.dais    =parameters[i,match("var.dais"    ,parnames)]
	if(is.null(rho.simple) | is.na(rho.simple)) rho.simple=rho.simple.fixed

	# Edit: do NOT include the DAIS noise, because it was modeled using paleo
	# data, which have very different uncertainties than modern instrumental
	# period data.
	err.temp = rep(sigma.T,n.time); err.temp[midx.temp]=sqrt(sigma.T^2 + obs.temp.err[oidx.temp]^2)
	err.ocheat = rep(sigma.H,n.time); err.ocheat[midx.ocheat]=sqrt(sigma.H^2+obs.ocheat.err[oidx.ocheat]^2)
	err.gsic = rep(sigma.gsic,n.time); err.gsic[midx.gsic]=sqrt(sigma.gsic^2+obs.gsic.err[oidx.gsic]^2)
	err.gis = rep(sigma.simple,n.time); err.gis[midx.gis]=sqrt(sigma.simple^2+obs.gis.err^2)

	temp.norm.stat[i,]   = temp.out.norm[i,]   + ar1.sim(n.time, rho.T, err.temp)
	ocheat.norm.stat[i,] = ocheat.out.norm[i,] + ar1.sim(n.time, rho.H, err.ocheat)
	gsic.norm.stat[i,]   = gsic.out.norm[i,]   + ar1.sim(n.time, rho.gsic, err.gsic)
	gis.norm.stat[i,]    = gis.out.norm[i,]    + ar1.sim(n.time, rho.simple, err.gis)
	ais.norm.stat[i,]    = ais.out.norm[i,]    #+ rnorm(  n.time, mean=0,sd=sqrt(var.dais))
	te.norm.stat[i,]     = te.out.norm[i,]

	slr.norm.stat[i,] = gsic.norm.stat[i,] +
										  gis.norm.stat[i,]  +
										  ais.norm.stat[i,]  +
										  te.norm.stat[i,]

  slr.norm.stat[i,] = slr.norm.stat[i,] - mean(slr.norm.stat[i,ind.norm.data[which(ind.norm.data[,1]=='sl'),2]:ind.norm.data[which(ind.norm.data[,1]=='sl'),3]])

  setTxtProgressBar(pb, i)
}
close(pb)

## Post-calibrate using Church and White (2011) data
## filter through a tunnel around the SLR data

## Make sure SLR data are also normalized
ibeg=which(obs.sl.time==mod.time[ind.norm.data[which(ind.norm.data[,1]=='sl'),2]])
iend=which(obs.sl.time==mod.time[ind.norm.data[which(ind.norm.data[,1]=='sl'),3]])
obs.sl = obs.sl - mean(obs.sl[ibeg:iend])

ind.window.mod = midx.sl
ind.window.obs = oidx.sl

#n.sigma=4	# use 4 to reproduce the spread in Church and White data well (nice envelope in hindcast)

if(n.sigma>0) {
	window.obs = mat.or.vec(2, length(ind.window.obs))
	window.obs[1,] = obs.sl[ind.window.obs]-n.sigma*obs.sl.err[ind.window.obs]
	window.obs[2,] = obs.sl[ind.window.obs]+n.sigma*obs.sl.err[ind.window.obs]
	window.mod = slr.norm.stat[,ind.window.mod]
	survive = rep(0, nrow(window.mod))
	for (i in 1:nrow(window.mod)) {
  	survive[i] = all( (window.mod[i,] > window.obs[1,]) & (window.mod[i,] < window.obs[2,]) )
  	if(is.na(survive[i])){survive[i]=0}
	}
	ind.survive = which( as.logical(survive))
} else if(n.sigma<0) {
	window.obs = mat.or.vec(2, length(ind.window.obs))
	window.obs[1,] = obs.sl[ind.window.obs]-abs(n.sigma)*obs.sl.err[ind.window.obs]
	window.obs[2,] = obs.sl[ind.window.obs]+abs(n.sigma)*obs.sl.err[ind.window.obs]
	window.mod = slr.norm.stat[,ind.window.mod]
	survive = rep(0, nrow(window.mod))
	for (i in 1:nrow(window.mod)) {
		ind.in = which(window.mod[i,]>window.obs[1,] & window.mod[i,]<window.obs[2,])
		prc.in = length(ind.in)/length(oidx.sl)
  	survive[i] = prc.in>prc.close
  	if(is.na(survive[i])){survive[i]=0}
	}
	ind.survive = which( as.logical(survive))
} else if(n.sigma==0) {
	ind.survive = 1:nrow(slr.norm.stat)
}

slr.out.good = slr.norm.stat[ind.survive,]
parameters.good = parameters[ind.survive,]
colnames(parameters.good) = parnames



##==============================================================================
##==============================================================================
## DAIS paleo runs with the post-calibrated parameters
## Already have dSL, SL, Tg.recon.

parameters=parameters.good

## How many members do you want in your ensemble?
n.ensemble = nrow(parameters)
n.parameters = ncol(parameters)

## Run only a sample of the DAIS paleo, because they are huge
#n.sample = 500
#ind.sample = sample( seq(1,nrow(parameters)), size=n.sample, replace=FALSE)
#parameters.sample = parameters[ind.sample,]

## Run all of them?
n.sample = n.ensemble
parameters.sample = parameters

source('../fortran/R/daisantoF.R')
n.paleo = length(SL)
dais.paleo = mat.or.vec(n.sample, n.paleo)
date = seq(-239999,16,1) #240 Kyr BP to 2100AD at one year intervals of the forcings
norm.period=c(1961,1990)
ibeg=which(date==(norm.period[1]-2000))
iend=which(date==(norm.period[2]-2000))
ind.norm.paleo=ibeg:iend
t.paleo = date

pb <- txtProgressBar(min=0,max=n.sample,initial=0,style=3);
for (i in 1:n.sample) {

	anto.a=parameters.sample[i,match("anto.a",parnames)]
  anto.b=parameters.sample[i,match("anto.b",parnames)]
  gamma =parameters.sample[i,match("gamma" ,parnames)]
  alpha =parameters.sample[i,match("alpha.dais" ,parnames)]
  mu =parameters.sample[i,match("mu" ,parnames)]
  nu =parameters.sample[i,match("nu" ,parnames)]
  P0 =parameters.sample[i,match("P0" ,parnames)]
  kappa =parameters.sample[i,match("kappa.dais" ,parnames)]
  f0 =parameters.sample[i,match("f0" ,parnames)]
  h0 =parameters.sample[i,match("h0" ,parnames)]
  c =parameters.sample[i,match("c" ,parnames)]
  b0 =parameters.sample[i,match("b0" ,parnames)]
  slope =parameters.sample[i,match("slope" ,parnames)]
  var.dais =parameters.sample[i,match("var.dais" ,parnames)]

	dais.tmp = daisantoF(
                       anto.a=anto.a, anto.b=anto.b,
                       slope.Ta2Tg=slope.Ta2Tg, intercept.Ta2Tg=intercept.Ta2Tg,
                       gamma=gamma  , alpha=alpha  ,
                       mu=mu        , nu=nu        ,
                       P0=P0        , kappa=kappa  ,
                       f0=f0        , h0=h0        ,
                       c=c          , b0=b0        ,
                       slope=slope  ,
                       Tg=Tg.recon  , SL=SL , dSL=dSL, includes_dSLais=1
											 )

	# Subtract off the 1961-1990 normalization period
	dais.norm = dais.tmp - mean(dais.tmp[ind.norm.paleo])

	# Add the modeled error back in
	dais.paleo[i,] = dais.norm + rnorm(n.paleo, mean=0,sd=sqrt(var.dais))

	setTxtProgressBar(pb, i)
}
close(pb)

parameters.good = parameters
colnames(parameters.good) = parnames

## Get 5-95% CI for hindcasts
## Initialize arrays for the output
dais.paleo.05 = rep(NA,length(date)); dais.paleo.50 = rep(NA,length(date)); dais.paleo.95 = rep(NA,length(date))
dais.paleo.max= rep(NA,length(date)); dais.paleo.min= rep(NA,length(date));

## Source a useful script, to allow for assigning multiple outputs at once
source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
pb <- txtProgressBar(min=0,max=length(date),initial=0,style=3);
for (t in 1:length(date)){
	c(dais.paleo.05[t] , dais.paleo.50[t] , dais.paleo.95[t] , dais.paleo.max[t], dais.paleo.min[t]) := quantile(dais.paleo[,t],c(0.05,.50,.95,1,0), na.rm=TRUE)
  setTxtProgressBar(pb, t)
}
close(pb)

## Make smoothed version of the AIS paleo results, so the plots are not massive
n.avg=100	# number of years in averaging period
n.time.avg=ceiling(length(date)/n.avg)
dais.paleo.05.avg = rep(NA,n.time.avg)
dais.paleo.50.avg = rep(NA,n.time.avg)
dais.paleo.95.avg = rep(NA,n.time.avg)
dais.paleo.max.avg = rep(NA,n.time.avg)
dais.paleo.min.avg = rep(NA,n.time.avg)
date.avg=seq(date[1],date[length(date)],by=n.avg)

pb <- txtProgressBar(min=0,max=n.time.avg,initial=0,style=3);
for (t in 1:(n.time.avg-1)){
	dais.paleo.05.avg[t] = mean(dais.paleo.05[((t-1)*n.avg+1) : (t*n.avg)])
	dais.paleo.50.avg[t] = mean(dais.paleo.50[((t-1)*n.avg+1) : (t*n.avg)])
	dais.paleo.95.avg[t] = mean(dais.paleo.95[((t-1)*n.avg+1) : (t*n.avg)])
	dais.paleo.max.avg[t] = mean(dais.paleo.max[((t-1)*n.avg+1) : (t*n.avg)])
	dais.paleo.min.avg[t] = mean(dais.paleo.min[((t-1)*n.avg+1) : (t*n.avg)])
  setTxtProgressBar(pb, t)
}
dais.paleo.05.avg[n.time.avg] = mean(dais.paleo.05[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.50.avg[n.time.avg] = mean(dais.paleo.50[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.95.avg[n.time.avg] = mean(dais.paleo.95[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.max.avg[n.time.avg] = mean(dais.paleo.max[((n.time.avg-1)*n.avg+1) : length(date)])
dais.paleo.min.avg[n.time.avg] = mean(dais.paleo.min[((n.time.avg-1)*n.avg+1) : length(date)])
close(pb)

##==============================================================================
##==============================================================================

## Save the hindcasts, trimming down first for ind.survive (Church and White, 2011)
## and then for AIS Vmin
gsic.hind = t(gsic.norm.stat[ind.survive,])
te.hind = t(te.norm.stat[ind.survive,])
gis.hind = t(gis.norm.stat[ind.survive,])
ais.hind = t(ais.norm.stat[ind.survive,])
temp.hind = t(temp.norm.stat[ind.survive,])
ocheat.hind = t(ocheat.norm.stat[ind.survive,])
gsl.hind = t(slr.out.good)
t.hind = mod.time

##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## Write post-calibrated parameter output file

## Fit multivariate kernel density estimate to the distributions of these
## post-calibrated parameter sets. THIS is the multivariate distribution of
## "good" parameters. We take 50-150% of the min/max range so as to ensure the
## post-calibrated parameter draws capture the tail behavior.
n.parameters=ncol(parameters.good)
pdf.all=vector('list',n.parameters)
n.node=100
for (pp in 1:n.parameters){
	lower.bound = 0.5*min(parameters.good[,pp])
	upper.bound = 1.5*max(parameters.good[,pp])
  tmp = density(parameters.good[,pp],kernel='gaussian',
                n=n.node, from=lower.bound, to=upper.bound)
  pdf.all[[pp]] = tmp; names(pdf.all)[pp]=parnames[pp]
}
bandwidths=rep(NA,n.parameters)
for (i in 1:n.parameters){
  bandwidths[i]=pdf.all[[i]]$bw
}

## Write these post-calibrated parameter sets to a CSV file for easy grabbing
## later? filename.parameters set at beginning of this script
to.file = rbind(parameters.good, bandwidths) # bandwidths are the last row
rownames(to.file)=NULL
colnames(to.file)=parnames
write.table(to.file, file=filename.parameters, sep=",", qmethod="double", row.names=FALSE)
print(paste("Writing post-calibrated parameters to file ",filename.parameters,sep=""))
##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## Pick up here?
if(FALSE){
	filename.postcalibration = "../output_calibration/BRICK-robustSLR_postcalibratedParameters_22Aug2016.csv"
	dat.brick = read.csv(filename.postcalibration)
	parameters.good = dat.brick[1:(nrow(dat.brick)-1), ]
	parnames.brick = colnames(parameters.good)
}

## Make RCP2.6, 4.5 and 8.5 projections to 2100. Need:
## (1) global total sea level, (2) global sea level without fast dynamics,
## (3) local (NOLA) sea level, (4) local sea level without fast dynamics

parameters=parameters.good
parnames = colnames(parameters)

## How many members do you want in your ensemble?
n.ensemble = nrow(parameters)
n.parameters = ncol(parameters)

## If asking for too many ensemble members, reset to just the parameters
n.ensemble = min(n.ensemble, nrow(parameters))

## Draw parameters for the ensemble members
parameters.ensemble = mat.or.vec(n.ensemble , n.parameters )
ind.ensemble = sample( seq(1,nrow(parameters)), size=n.ensemble, replace=FALSE)
for (p in 1:n.parameters){
	for (i in 1:n.ensemble){
#	  parameters.ensemble[i,p] = rnorm( 1, mean=parameters[ind.ensemble[i],p], sd=bandwidths[1,p])
	  parameters.ensemble[i,p] = parameters[ind.ensemble[i],p]
	}
}

## Set up the model projections
l.project = TRUE
begyear = 1850
endyear = 2100; #if(!l.project & endyear>2009) print('l.project and endyear not compatible')
mod.time= begyear:endyear
t.proj = mod.time # save for output

## NOTE: IPCC generally are relative to 1986-2005.
## Therefore use that period to be commensurate with their results
begyear.norm = 1986
endyear.norm = 2005
ind.norm = which(mod.time==begyear.norm):which(mod.time==endyear.norm)
n.time = length(mod.time)

## Want model projections under RCP2.6, 4.5 and 8.5. Initialize a list, each
## will hold the ensemble model output for a different forcing scenario.
n.scen = 3
proj.out = vector("list", n.scen)
forc.scen = vector("list", n.scen)
badruns.scen = vector("list", n.scen)
names.scen = c('rcp26','rcp45','rcp85')
names(proj.out) = names.scen

## Get the forcing data
forc.scen[[1]] = read.csv( '../data/forcing_rcp26.csv', header=TRUE )
forc.scen[[2]] = read.csv( '../data/forcing_rcp45.csv', header=TRUE )
forc.scen[[3]] = read.csv( '../data/forcing_rcp85.csv', header=TRUE )

## Loop over forcing scenarios
for (ff in 1:n.scen) {

	forcing = forc.scen[[ff]]

	## Initialize matrix to store model ensemble output and flag for bad runs
	## (could go bad because of different forcing from hindcasts; different
	## temperatures lead to different stability requirements for SIMPLE)
	brick.out = vector("list", n.ensemble)
	badruns   = rep(0, n.ensemble)

	## Run the sample, and enjoy a nice progress bar
	pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3);
	for (i in 1:n.ensemble) {

		brick.out[[i]] = brick_model(	parameters.in			= parameters.ensemble[i,],
																	parnames.in				= parnames,
																	forcing.in				= forcing,
																	l.project					= l.project,
																	slope.Ta2Tg.in		= slope.Ta2Tg,
																	intercept.Ta2Tg.in= intercept.Ta2Tg,
																	mod.time					= mod.time,
																	ind.norm.data 		= ind.norm.data,
																	ind.norm.sl 			= ind.norm,
																	luse.brick 				= luse.brick
																	)

		# check if the run turned out bad
		if( is.na(brick.out[[i]]$slr.out[length(mod.time)]) ) {badruns[i]=1}

  	setTxtProgressBar(pb, i)
	}
	close(pb)

	badruns.scen[[ff]] = badruns
	proj.out[[ff]] = brick.out
}

## Filter out any bad runs
ind.good = which(badruns.scen[[1]]==0 & badruns.scen[[2]]==0 & badruns.scen[[3]]==0)
n.ensemble = length(ind.good)
parameters.ensemble = parameters.ensemble[ind.good,]

## Trim down the hindcasts to match (some runs might go off the rails with RCP
## forcings but didn't with the hindcast forcings)
gsl.hind = gsl.hind[,ind.good]
gsic.hind = gsic.hind[,ind.good]
te.hind = te.hind[,ind.good]
gis.hind = gis.hind[,ind.good]
ais.hind = ais.hind[,ind.good]
temp.hind = temp.hind[,ind.good]
ocheat.hind = ocheat.hind[,ind.good]

brick.rcp26 = proj.out[[1]][ind.good]
brick.rcp45 = proj.out[[2]][ind.good]
brick.rcp85 = proj.out[[3]][ind.good]

## Gather the fields for each simulation (easy referencing for plotting and
## analysis)
proj.rcp26 = vector("list", n.scen)
proj.rcp45 = vector("list", n.scen)
proj.rcp85 = vector("list", n.scen)

## Make each projections list a list of the SLR contributions and other fields
## Get the 90% CI at each time step, of all model output fields

names.output = c('slr','gsic','gis','ais','te','temp','ocheat')
n.output = length(names.output)

proj.rcp26 = vector("list", n.output)
proj.rcp45 = vector("list", n.output)
proj.rcp85 = vector("list", n.output)

names(proj.rcp26) = names.output
names(proj.rcp45) = names.output
names(proj.rcp85) = names.output

for (j in 1:n.output) {
	proj.rcp26[[j]] = mat.or.vec(n.ensemble, n.time)
	proj.rcp45[[j]] = mat.or.vec(n.ensemble, n.time)
	proj.rcp85[[j]] = mat.or.vec(n.ensemble, n.time)
}

## Go through each simulation and collect, normalize and add modeled error to
## the fields
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {

	T0=parameters.ensemble[i,match("T0",parnames)]
	H0=parameters.ensemble[i,match("H0",parnames)]

	# set the results (note that we already have slr.out)
	proj.rcp26$gsic[i,] = brick.rcp26[[i]]$gsic.out
	proj.rcp26$gis[i,] = brick.rcp26[[i]]$simple.out$sle.gis
	proj.rcp26$ais[i,] = brick.rcp26[[i]]$dais.out
	proj.rcp26$te[i,] = brick.rcp26[[i]]$te.out
	proj.rcp26$temp[i,] = brick.rcp26[[i]]$doeclim.out$temp + T0
	proj.rcp26$ocheat[i,] = brick.rcp26[[i]]$doeclim.out$ocheat + H0

	proj.rcp45$gsic[i,] = brick.rcp45[[i]]$gsic.out
	proj.rcp45$gis[i,] = brick.rcp45[[i]]$simple.out$sle.gis
	proj.rcp45$ais[i,] = brick.rcp45[[i]]$dais.out
	proj.rcp45$te[i,] = brick.rcp45[[i]]$te.out
	proj.rcp45$temp[i,] = brick.rcp45[[i]]$doeclim.out$temp + T0
	proj.rcp45$ocheat[i,] = brick.rcp45[[i]]$doeclim.out$ocheat + H0

	proj.rcp85$gsic[i,] = brick.rcp85[[i]]$gsic.out
	proj.rcp85$gis[i,] = brick.rcp85[[i]]$simple.out$sle.gis
	proj.rcp85$ais[i,] = brick.rcp85[[i]]$dais.out
	proj.rcp85$te[i,] = brick.rcp85[[i]]$te.out
	proj.rcp85$temp[i,] = brick.rcp85[[i]]$doeclim.out$temp + T0
	proj.rcp85$ocheat[i,] = brick.rcp85[[i]]$doeclim.out$ocheat + H0

	# Normalize the output to "ind.norm" (1961-1990? 1986-2005 (Mengel)?).
	# Normalize ocean heat uptake too, for sake of plotting (it will be plotted as
	# the amount of heat taken up by ocean since 1986-2005 period)
	proj.rcp26$gsic[i,] = proj.rcp26$gsic[i,] - mean( proj.rcp26$gsic[i,ind.norm])
	proj.rcp26$gis[i,] = proj.rcp26$gis[i,] - mean( proj.rcp26$gis[i,ind.norm])
	proj.rcp26$ais[i,] = proj.rcp26$ais[i,] - mean( proj.rcp26$ais[i,ind.norm])
	proj.rcp26$te[i,] = proj.rcp26$te[i,] - mean( proj.rcp26$te[i,ind.norm])
	proj.rcp26$temp[i,] = proj.rcp26$temp[i,] - mean( proj.rcp26$temp[i,ind.norm])
	proj.rcp26$ocheat[i,] = proj.rcp26$ocheat[i,] - mean( proj.rcp26$ocheat[i,ind.norm])

	proj.rcp45$gsic[i,] = proj.rcp45$gsic[i,] - mean( proj.rcp45$gsic[i,ind.norm])
	proj.rcp45$gis[i,] = proj.rcp45$gis[i,] - mean( proj.rcp45$gis[i,ind.norm])
	proj.rcp45$ais[i,] = proj.rcp45$ais[i,] - mean( proj.rcp45$ais[i,ind.norm])
	proj.rcp45$te[i,] = proj.rcp45$te[i,] - mean( proj.rcp45$te[i,ind.norm])
	proj.rcp45$temp[i,] = proj.rcp45$temp[i,] - mean( proj.rcp45$temp[i,ind.norm])
	proj.rcp45$ocheat[i,] = proj.rcp45$ocheat[i,] - mean( proj.rcp45$ocheat[i,ind.norm])

	proj.rcp85$gsic[i,] = proj.rcp85$gsic[i,] - mean( proj.rcp85$gsic[i,ind.norm])
	proj.rcp85$gis[i,] = proj.rcp85$gis[i,] - mean( proj.rcp85$gis[i,ind.norm])
	proj.rcp85$ais[i,] = proj.rcp85$ais[i,] - mean( proj.rcp85$ais[i,ind.norm])
	proj.rcp85$te[i,] = proj.rcp85$te[i,] - mean( proj.rcp85$te[i,ind.norm])
	proj.rcp85$temp[i,] = proj.rcp85$temp[i,] - mean( proj.rcp85$temp[i,ind.norm])
	proj.rcp85$ocheat[i,] = proj.rcp85$ocheat[i,] - mean( proj.rcp85$ocheat[i,ind.norm])

	# Add the statistcal model for AR1 (or otherwise) noise?
	# Edit: not for DAIS, because that noise matches paleo data, which has
	# different uncertainties than modern data.
	sigma.T  =parameters.ensemble[i,match("sigma.T"  ,parnames)]
	rho.T    =parameters.ensemble[i,match("rho.T"    ,parnames)]
	sigma.H  =parameters.ensemble[i,match("sigma.H"  ,parnames)]
	rho.H    =parameters.ensemble[i,match("rho.H"    ,parnames)]
	sigma.gsic  =parameters.ensemble[i,match("sigma.gsic"  ,parnames)]
	rho.gsic    =parameters.ensemble[i,match("rho.gsic"    ,parnames)]
	sigma.simple=parameters.ensemble[i,match("sigma.simple",parnames)]
	rho.simple  =parameters.ensemble[i,match("rho.simple"  ,parnames)]
	var.dais    =parameters.ensemble[i,match("var.dais"    ,parnames)]
	if(is.na(rho.simple)) rho.simple=rho.simple.fixed

	proj.rcp26$temp[i,] = proj.rcp26$temp[i,] + ar1.sim(n.time, rho.T, sigma.T)
	proj.rcp26$ocheat[i,] = proj.rcp26$ocheat[i,] + ar1.sim(n.time, rho.H, sigma.H)
	proj.rcp26$gsic[i,] = proj.rcp26$gsic[i,] + ar1.sim(n.time, rho.gsic, sigma.gsic)
	proj.rcp26$gis[i,] = proj.rcp26$gis[i,] + ar1.sim(n.time, rho.simple, sigma.simple)

	proj.rcp45$temp[i,] = proj.rcp45$temp[i,] + ar1.sim(n.time, rho.T, sigma.T)
	proj.rcp45$ocheat[i,] = proj.rcp45$ocheat[i,] + ar1.sim(n.time, rho.H, sigma.H)
	proj.rcp45$gsic[i,] = proj.rcp45$gsic[i,] + ar1.sim(n.time, rho.gsic, sigma.gsic)
	proj.rcp45$gis[i,] = proj.rcp45$gis[i,] + ar1.sim(n.time, rho.simple, sigma.simple)

	proj.rcp85$temp[i,] = proj.rcp85$temp[i,] + ar1.sim(n.time, rho.T, sigma.T)
	proj.rcp85$ocheat[i,] = proj.rcp85$ocheat[i,] + ar1.sim(n.time, rho.H, sigma.H)
	proj.rcp85$gsic[i,] = proj.rcp85$gsic[i,] + ar1.sim(n.time, rho.gsic, sigma.gsic)
	proj.rcp85$gis[i,] = proj.rcp85$gis[i,] + ar1.sim(n.time, rho.simple, sigma.simple)

	# Add up to total sea-level rise
	proj.rcp26$slr[i,] = proj.rcp26$gsic[i,] +
												proj.rcp26$gis[i,] +
												proj.rcp26$ais[i,] +
												proj.rcp26$te[i,]

	proj.rcp45$slr[i,] = proj.rcp45$gsic[i,] +
												proj.rcp45$gis[i,] +
												proj.rcp45$ais[i,] +
												proj.rcp45$te[i,]

	proj.rcp85$slr[i,] = proj.rcp85$gsic[i,] +
												proj.rcp85$gis[i,] +
												proj.rcp85$ais[i,] +
												proj.rcp85$te[i,]


	# And normalize sea-level rise
	proj.rcp26$slr[i,] = proj.rcp26$slr[i,] - mean(proj.rcp26$slr[i,ind.norm])
	proj.rcp45$slr[i,] = proj.rcp45$slr[i,] - mean(proj.rcp45$slr[i,ind.norm])
	proj.rcp85$slr[i,] = proj.rcp85$slr[i,] - mean(proj.rcp85$slr[i,ind.norm])

  setTxtProgressBar(pb, i)
}
close(pb)
##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## Write a netCDF ensemble output file including each of the RCP scenarios for
## global total sea level and each contribution to global sea level rise, for
## the hindcast plots and projections.

library(ncdf4)

dim.tproj <- ncdim_def('time_proj', 'years', as.double(t.proj))
dim.ensemble <- ncdim_def('ens', 'ensemble member', as.double(1:nrow(proj.rcp26$slr)))
dim.thind <- ncdim_def('time_hind', 'years', as.double(t.hind))
dim.tpaleo <- ncdim_def('time_paleo', 'year paleo', as.double(t.paleo))
dim.tpaleo.avg <- ncdim_def('time_paleo_avg', 'year avg paleo', as.double(date.avg))

ais.paleo.05 <- ncvar_def('AIS_paleo_q05', 'meters', list(dim.tpaleo), -999,
                  longname = 'AIS paleo contribution to sea level (5th quantile)')
ais.paleo.50 <- ncvar_def('AIS_paleo_q50', 'meters', list(dim.tpaleo), -999,
                  longname = 'AIS paleo contribution to sea level (median)')
ais.paleo.95 <- ncvar_def('AIS_paleo_q95', 'meters', list(dim.tpaleo), -999,
                  longname = 'AIS paleo contribution to sea level (95th quantile)')
ais.paleo.max <- ncvar_def('AIS_paleo_max', 'meters', list(dim.tpaleo), -999,
                  longname = 'AIS paleo contribution to sea level (maximum)')
ais.paleo.min <- ncvar_def('AIS_paleo_min', 'meters', list(dim.tpaleo), -999,
                  longname = 'AIS paleo contribution to sea level (minimum)')

ais.paleo.05.avg <- ncvar_def('AIS_paleo_avg_q05', 'meters', list(dim.tpaleo.avg), -999,
                  longname = 'AIS paleo contribution to sea level (smoothed, 5th quantile)')
ais.paleo.50.avg <- ncvar_def('AIS_paleo_avg_q50', 'meters', list(dim.tpaleo.avg), -999,
                  longname = 'AIS paleo contribution to sea level (smoothed, median)')
ais.paleo.95.avg <- ncvar_def('AIS_paleo_avg_q95', 'meters', list(dim.tpaleo.avg), -999,
                  longname = 'AIS paleo contribution to sea level (smoothed, 95th quantile)')
ais.paleo.max.avg <- ncvar_def('AIS_paleo_avg_max', 'meters', list(dim.tpaleo.avg), -999,
                  longname = 'AIS paleo contribution to sea level (smoothed, maximum)')
ais.paleo.min.avg <- ncvar_def('AIS_paleo_avg_min', 'meters', list(dim.tpaleo.avg), -999,
                  longname = 'AIS paleo contribution to sea level (smoothed, minimum)')

gsl.hindcast <- ncvar_def('GlobalSeaLevel_hind', 'meters', list(dim.thind, dim.ensemble), -999,
                  		longname = 'Global sea level (hindcast)')
gsic.hindcast <- ncvar_def('GSIC_hind', 'meters', list(dim.thind, dim.ensemble), -999,
                  		longname = 'Glaciers and small ice caps contribution to GSL (hindcast)')
te.hindcast <- ncvar_def('TE_hind', 'meters', list(dim.thind, dim.ensemble), -999,
                  		longname = 'Thermal expansion contribution to GSL (hindcast)')
gis.hindcast <- ncvar_def('GIS_hind', 'meters', list(dim.thind, dim.ensemble), -999,
                  		longname = 'Greenland ice sheet contribution to GSL (hindcast)')
ais.hindcast <- ncvar_def('AIS_hind', 'meters', list(dim.thind, dim.ensemble), -999,
                  		longname = 'Antarctic ice sheet contribution to GSL (hindcast)')
temp.hindcast <- ncvar_def('temp_hind', 'deg C', list(dim.thind, dim.ensemble), -999,
                  		longname = 'Global average surface temperature anomaly (hindcast)')
ocheat.hindcast <- ncvar_def('ocheat_hind', '10^22 J', list(dim.thind, dim.ensemble), -999,
                  		longname = 'Ocean heat uptake (hindcast)')

gsic.rcp26 <- ncvar_def('GSIC_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GSIC contribution to sea level (RCP26)')
te.rcp26 <- ncvar_def('TE_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'TE contribution to sea level (RCP26)')
gis.rcp26 <- ncvar_def('GIS_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GIS contribution to sea level (RCP26)')
ais.rcp26 <- ncvar_def('AIS_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'AIS contribution to sea level (RCP26)')
temp.rcp26 <- ncvar_def('temp_RCP26', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP26)')
ocheat.rcp26 <- ncvar_def('ocheat_RCP26', '10^22 J', list(dim.tproj, dim.ensemble), -999,
                  longname = 'ocean heat uptake (RCP26)')

gsl.rcp26 <- ncvar_def('GlobalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level accounting for fast dynamics (RCP26)')
gsic.rcp26 <- ncvar_def('GSIC_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GSIC contribution to sea level (RCP26)')
te.rcp26 <- ncvar_def('TE_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'TE contribution to sea level (RCP26)')
gis.rcp26 <- ncvar_def('GIS_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GIS contribution to sea level (RCP26)')
ais.rcp26 <- ncvar_def('AIS_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'AIS contribution to sea level (RCP26)')
temp.rcp26 <- ncvar_def('temp_RCP26', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP26)')
ocheat.rcp26 <- ncvar_def('ocheat_RCP26', '10^22 J', list(dim.tproj, dim.ensemble), -999,
                  longname = 'ocean heat uptake (RCP26)')

gsl.rcp45 <- ncvar_def('GlobalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level accounting for fast dynamics (RCP45)')
gsic.rcp45 <- ncvar_def('GSIC_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GSIC contribution to sea level (RCP45)')
te.rcp45 <- ncvar_def('TE_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'TE contribution to sea level (RCP45)')
gis.rcp45 <- ncvar_def('GIS_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GIS contribution to sea level (RCP45)')
ais.rcp45 <- ncvar_def('AIS_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'AIS contribution to sea level (RCP45)')
temp.rcp45 <- ncvar_def('temp_RCP45', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP45)')
ocheat.rcp45 <- ncvar_def('ocheat_RCP45', '10^22 J', list(dim.tproj, dim.ensemble), -999,
                  longname = 'ocean heat uptake (RCP45)')

gsl.rcp85 <- ncvar_def('GlobalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Global sea level accounting for fast dynamics (RCP85)')
gsic.rcp85 <- ncvar_def('GSIC_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GSIC contribution to sea level (RCP85)')
te.rcp85 <- ncvar_def('TE_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'TE contribution to sea level (RCP85)')
gis.rcp85 <- ncvar_def('GIS_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'GIS contribution to sea level (RCP85)')
ais.rcp85 <- ncvar_def('AIS_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'AIS contribution to sea level (RCP85)')
temp.rcp85 <- ncvar_def('temp_RCP85', 'deg C', list(dim.tproj, dim.ensemble), -999,
                  longname = 'global mean surface temperature anomaly (RCP85)')
ocheat.rcp85 <- ncvar_def('ocheat_RCP85', '10^22 J', list(dim.tproj, dim.ensemble), -999,
                  longname = 'ocean heat uptake (RCP85)')

today=Sys.Date(); today=format(today,format="%d%b%Y")

outnc <- nc_create(filename.brickout,
										list( gsl.rcp26, gsl.rcp45, gsl.rcp85,
													gsic.rcp26, te.rcp26, gis.rcp26, ais.rcp26, temp.rcp26, ocheat.rcp26,
													gsic.rcp45, te.rcp45, gis.rcp45, ais.rcp45, temp.rcp45, ocheat.rcp45,
													gsic.rcp85, te.rcp85, gis.rcp85, ais.rcp85, temp.rcp85, ocheat.rcp85,
													gsl.hindcast, gsic.hindcast, te.hindcast, gis.hindcast, ais.hindcast, temp.hindcast, ocheat.hindcast,
													ais.paleo.05, ais.paleo.50, ais.paleo.95, ais.paleo.max, ais.paleo.min,
													ais.paleo.05.avg, ais.paleo.50.avg, ais.paleo.95.avg, ais.paleo.max.avg, ais.paleo.min.avg),
										force_v4 = TRUE)

ncvar_put(outnc, gsl.rcp26, t(proj.rcp26$slr))
ncvar_put(outnc, gsic.rcp26, t(proj.rcp26$gsic))
ncvar_put(outnc, te.rcp26, t(proj.rcp26$te))
ncvar_put(outnc, gis.rcp26, t(proj.rcp26$gis))
ncvar_put(outnc, ais.rcp26, t(proj.rcp26$ais))
ncvar_put(outnc, temp.rcp26, t(proj.rcp26$temp))
ncvar_put(outnc, ocheat.rcp26, t(proj.rcp26$ocheat))

ncvar_put(outnc, gsl.rcp45, t(proj.rcp45$slr))
ncvar_put(outnc, gsic.rcp45, t(proj.rcp45$gsic))
ncvar_put(outnc, te.rcp45, t(proj.rcp45$te))
ncvar_put(outnc, gis.rcp45, t(proj.rcp45$gis))
ncvar_put(outnc, ais.rcp45, t(proj.rcp45$ais))
ncvar_put(outnc, temp.rcp45, t(proj.rcp45$temp))
ncvar_put(outnc, ocheat.rcp45, t(proj.rcp45$ocheat))

ncvar_put(outnc, gsl.rcp85, t(proj.rcp85$slr))
ncvar_put(outnc, gsic.rcp85, t(proj.rcp85$gsic))
ncvar_put(outnc, te.rcp85, t(proj.rcp85$te))
ncvar_put(outnc, gis.rcp85, t(proj.rcp85$gis))
ncvar_put(outnc, ais.rcp85, t(proj.rcp85$ais))
ncvar_put(outnc, temp.rcp85, t(proj.rcp85$temp))
ncvar_put(outnc, ocheat.rcp85, t(proj.rcp85$ocheat))

ncvar_put(outnc, gsl.hindcast, gsl.hind)
ncvar_put(outnc, gsic.hindcast, gsic.hind)
ncvar_put(outnc, te.hindcast, te.hind)
ncvar_put(outnc, gis.hindcast, gis.hind)
ncvar_put(outnc, ais.hindcast, ais.hind)
ncvar_put(outnc, temp.hindcast, temp.hind)
ncvar_put(outnc, ocheat.hindcast, ocheat.hind)

ncvar_put(outnc, ais.paleo.05, dais.paleo.05)
ncvar_put(outnc, ais.paleo.50, dais.paleo.50)
ncvar_put(outnc, ais.paleo.95, dais.paleo.95)
ncvar_put(outnc, ais.paleo.max, dais.paleo.max)
ncvar_put(outnc, ais.paleo.min, dais.paleo.min)
ncvar_put(outnc, ais.paleo.05.avg, dais.paleo.05.avg)
ncvar_put(outnc, ais.paleo.50.avg, dais.paleo.50.avg)
ncvar_put(outnc, ais.paleo.95.avg, dais.paleo.95.avg)
ncvar_put(outnc, ais.paleo.max.avg, dais.paleo.max.avg)
ncvar_put(outnc, ais.paleo.min.avg, dais.paleo.min.avg)

nc_close(outnc)

##==============================================================================
##==============================================================================

t.end = proc.time()




##==============================================================================
## End
##==============================================================================
