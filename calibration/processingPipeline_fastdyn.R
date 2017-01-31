##==============================================================================
##
## Pipeline for processing DAIS-fast dynamics calibration results and BRICK
## rest-of-model (ROM) calibration results.
##
##  Required input:
##    filename.DAIScalibration  DAIS-fast dynamics calibration parameter posterior draws
##    filename.BRICKcalibration BRICK rest-of-model calibration parameter posterior draws
##
##  Output:
##    BRICK-fastdyn_postcalibratedParameters_[datestamp].csv
##                            post-calibrated parameters file (CSV format)
##		BRICK-fastdyn_physical_[uniform/gamma]_[datestamp].nc
##														BRICK physical model output (netcdf format)
##		vanDantzig_RCP85_[uniform/gamma]_[datestamp].nc
##														Van Dantzig optimization output (netcdf format)
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

rm(list=ls())

library(ncdf4)

t.beg = proc.time()

##==============================================================================
##==============================================================================
## Define the files you want to process

filename.BRICKcalibration = "../output_calibration/BRICK_calibratedParameters_12Aug2016.nc"
#filename.BRICKcalibration = "../output_calibration/BRICK_calibratedParameters_12Aug2016.csv"
filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_06Sep2016.csv"

priors='u'	## Which priors? u=uniform, g=gamma
appen=''		## Append file name? In case you process multiple files in one day
today=Sys.Date(); today=format(today,format="%d%b%Y")

## Define the files you want to read/create
if(priors=='u'){
	filename.DAIScalibration = "../output_calibration/DAISfastdyn_calibratedParameters_uniform_29Jan2017.nc"
	filename.parameters = paste('../output_calibration/BRICK-fastdyn_postcalibratedParameters_uniform_',today,appen,'.csv', sep="")
	filename.brickout = paste('../output_model/BRICK-fastdyn_physical_uniform_',today,appen,'.nc',sep="")
	filename.vdout = paste('../output_model/vanDantzig_RCP85_uniform_',today,appen,'.nc',sep="")
}
if(priors=='g'){
	filename.DAIScalibration = "../output_calibration/DAISfastdyn_calibratedParameters_gamma_29Jan2017.nc"
	filename.parameters = paste('../output_calibration/BRICK-fastdyn_postcalibratedParameters_gamma_',today,appen,'.csv', sep="")
	filename.brickout = paste('../output_model/BRICK-fastdyn_physical_gamma_',today,appen,'.nc',sep="")
	filename.vdout = paste('../output_model/vanDantzig_RCP85_gamma_',today,appen,'.nc',sep="")
}

n.ensemble = 30000
n.ensemble.report = n.ensemble

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
## Combine calibrated parameters from DAIS-fast dynamics and BRICK

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
source('../fortran/R/daisanto_fastdynF.R')   		# DAIS (Antarctic Ice Sheet) model, with 'anto'
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

source('../R/BRICK_coupledModel_fastdyn.R')

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
disint.out = mat.or.vec(n.ensemble,length(mod.time))

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
	ais.out[i,]    = brick.out[[i]]$dais.out$Vais
	disint.out[i,] = brick.out[[i]]$dais.out$Vdisint
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

# Rejection sampling, with target distribution as the likelihood of the sea-
# level rise data, proposing uniformly across the ensemble members.
survive = rep(0, n.ensemble)

## Make sure SLR data are also normalized
ibeg=which(obs.sl.time==mod.time[ind.norm.data[which(ind.norm.data[,1]=='sl'),2]])
iend=which(obs.sl.time==mod.time[ind.norm.data[which(ind.norm.data[,1]=='sl'),3]])
obs.sl = obs.sl - mean(obs.sl[ibeg:iend])


if(FALSE) {

	# calibrate to the Church and White data with an assumed closed sea level
	# budget (TE+AIS+GIS+GSIC = GMSL)

# maximum likelihood of sea-level rise data would occur if you saw the exact
# same data out in the wild
resid = obs.sl-obs.sl
lik.max = sum(dnorm(resid,sd=obs.sl.err,log=TRUE))/n.ensemble

# suppose each of the concomitant parameter combinations yields a sea-level rise
# simulation that is uniformly distributed with probability lik.max.
# this is certainly an upper bound on the likelihood function given the actual
# sea level rise data, evaluated at the simulated sea level rise (lik.mem).
# the ratio of lik.mem to lik.max is the acceptance ratio for rejection sampling
uni.rnd = log(runif(n.ensemble))
for (i in 1:n.ensemble) {
	resid = obs.sl[oidx.sl] - slr.norm.stat[i,midx.sl]
	lik.mem = sum(dnorm(resid, sd=obs.sl.err, log=TRUE))
	if( uni.rnd[i] <= lik.mem-lik.max) {survive[i]=1}
}
ind.survive = which( as.logical(survive))
print(paste('Calibration to sea level data by rejection sampling leaves ',length(ind.survive),' full calibrated ensemble members',sep=''))

slr.out.good = slr.norm.stat[ind.survive,]
parameters.good = parameters[ind.survive,]
colnames(parameters.good) = parnames

} else {

	# calibrate to the Church and White data with land water subtracted out
	# and uncertainties added in quadrature
	# assumed budget: TE+AIS+GIS+GSIC+LWS = GMSL
	#1901-1990: –0.11 [–0.16 to –0.06] (5-95% range)
	lw.time.1900 <- 1900:1989
	i1900 <- which(mod.time==lw.time.1900[1]):which(mod.time==lw.time.1900[length(lw.time.1900)])
	lw.1900 <- (-0.11/1000)*(lw.time.1900 - 1900)
	lw.err.1900 <- (0.25*(-0.06--0.16)/1000)*sqrt(lw.time.1900 - lw.time.1900[1])
	#1971-2010: 0.12 [0.03 to 0.22]
	lw.time.1970 <- 1970:2009
	i1970 <- which(mod.time==lw.time.1970[1]):which(mod.time==lw.time.1970[length(lw.time.1970)])
	lw.1970 <- (0.12/1000)*(lw.time.1970 - lw.time.1970[1])
	lw.err.1970 <- (0.25*(0.2-0.03)/1000)*sqrt(lw.time.1970 - lw.time.1970[1])
	#1993-2010: 0.38 [0.26 to 0.49]
	lw.time.1992 <- 1992:2009
	i1992 <- which(mod.time==lw.time.1992[1]):which(mod.time==lw.time.1992[length(lw.time.1992)])
	lw.1992 <- (0.38/1000)*(lw.time.1992 - lw.time.1992[1])
	lw.err.1992 <- (0.25*(0.49-0.26)/1000)*sqrt(lw.time.1992 - lw.time.1992[1])

	# normalize, subtract and add error in quadrature
	obs.sl.lw.1900 <- obs.sl[which(obs.sl.time==lw.time.1900[1]):which(obs.sl.time==lw.time.1900[length(lw.time.1900)])] - obs.sl[which(obs.sl.time==lw.time.1900[1])]
	obs.sl.lw.1970 <- obs.sl[which(obs.sl.time==lw.time.1970[1]):which(obs.sl.time==lw.time.1970[length(lw.time.1970)])] - obs.sl[which(obs.sl.time==lw.time.1970[1])]
	obs.sl.lw.1992 <- obs.sl[which(obs.sl.time==lw.time.1992[1]):which(obs.sl.time==lw.time.1992[length(lw.time.1992)])] - obs.sl[which(obs.sl.time==lw.time.1992[1])]

	obs.sl.lw.1900 <- obs.sl.lw.1900 - lw.1900
	obs.sl.lw.1970 <- obs.sl.lw.1970 - lw.1970
	obs.sl.lw.1992 <- obs.sl.lw.1992 - lw.1992

	obs.sl.lw.err.1900 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1900[1]):which(obs.sl.time==lw.time.1900[length(lw.time.1900)])]^2 + lw.err.1900^2)
	obs.sl.lw.err.1970 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1970[1]):which(obs.sl.time==lw.time.1970[length(lw.time.1970)])]^2 + lw.err.1970^2)
	obs.sl.lw.err.1992 <- sqrt(obs.sl.err[which(obs.sl.time==lw.time.1992[1]):which(obs.sl.time==lw.time.1992[length(lw.time.1992)])]^2 + lw.err.1992^2)

	# calculate likelihood as the product of the three independent likelihoods
	resid.1900 <- obs.sl.lw.1900 - obs.sl.lw.1900
	llik.1900 <- sum(dnorm(resid.1900, sd=obs.sl.lw.err.1900, log=TRUE))
	resid.1970 <- obs.sl.lw.1970 - obs.sl.lw.1970
	llik.1970 <- sum(dnorm(resid.1970, sd=obs.sl.lw.err.1970, log=TRUE))
	resid.1992 <- obs.sl.lw.1992 - obs.sl.lw.1992
	llik.1992 <- sum(dnorm(resid.1992, sd=obs.sl.lw.err.1992, log=TRUE))
	lik.max <- (llik.1900 + llik.1970 + llik.1992)/n.ensemble

	imod.1900 <- which(mod.time==lw.time.1900[1]):which(mod.time==lw.time.1900[length(lw.time.1900)])
	imod.1970 <- which(mod.time==lw.time.1970[1]):which(mod.time==lw.time.1970[length(lw.time.1970)])
	imod.1992 <- which(mod.time==lw.time.1992[1]):which(mod.time==lw.time.1992[length(lw.time.1992)])

	uni.rnd = log(runif(n.ensemble))
	for (i in 1:n.ensemble) {
		resid.1900 <- obs.sl.lw.1900 - (slr.norm.stat[i,imod.1900]-slr.norm.stat[i,imod.1900[1]])
		resid.1970 <- obs.sl.lw.1970 - (slr.norm.stat[i,imod.1970]-slr.norm.stat[i,imod.1970[1]])
		resid.1992 <- obs.sl.lw.1992 - (slr.norm.stat[i,imod.1992]-slr.norm.stat[i,imod.1992[1]])
		llik.1900 <- sum(dnorm(resid.1900, sd=obs.sl.lw.err.1900, log=TRUE))
		llik.1970 <- sum(dnorm(resid.1970, sd=obs.sl.lw.err.1970, log=TRUE))
		llik.1992 <- sum(dnorm(resid.1992, sd=obs.sl.lw.err.1992, log=TRUE))
		lik.mem <- llik.1900 + llik.1970 + llik.1992
		if( uni.rnd[i] <= lik.mem-lik.max) {survive[i]=1}
	}
	ind.survive = which( as.logical(survive))
	print(paste('Calibration to sea level data by rejection sampling leaves ',length(ind.survive),' full calibrated ensemble members',sep=''))

	slr.out.good = slr.norm.stat[ind.survive,]
	parameters.good = parameters[ind.survive,]
	colnames(parameters.good) = parnames

}


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

source('../fortran/R/daisanto_fastdynF.R')
n.paleo = length(SL)
dais.paleo = mat.or.vec(n.sample, n.paleo)
date = seq(-239999,16,1) #240 Kyr BP to 2100AD at one year intervals of the forcings
norm.period=c(1961,1990)
ibeg=which(date==(norm.period[1]-2000))
iend=which(date==(norm.period[2]-2000))
ind.norm.paleo=ibeg:iend
t.paleo = date

## Post-calibrate also on the minimum LIG AIS volume?
ind.vmin = rep(0,n.ensemble)
#Vmin = 18e15	# minimum AIS volume in LIG (m^3)
Vmin = 0	# minimum AIS volume in LIG (m^3)

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
  lambda =parameters.sample[i,match("lambda" ,parnames)]
  Tcrit =parameters.sample[i,match("Tcrit" ,parnames)]
  var.dais =parameters.sample[i,match("var.dais" ,parnames)]

	dais.tmp = daisanto_fastdynF(
                       anto.a=anto.a, anto.b=anto.b,
                       slope.Ta2Tg=slope.Ta2Tg, intercept.Ta2Tg=intercept.Ta2Tg,
                       gamma=gamma  , alpha=alpha  ,
                       mu=mu        , nu=nu        ,
                       P0=P0        , kappa=kappa  ,
                       f0=f0        , h0=h0        ,
                       c=c          , b0=b0        ,
                       slope=slope  ,
                       Tg=Tg.recon  , SL=SL , dSL=dSL, includes_dSLais=1,
											 Tcrit=Tcrit	, lambda=lambda
											 )

	# Subtract off the 1961-1990 normalization period
	dais.norm = dais.tmp$Vais - mean(dais.tmp$Vais[ind.norm.paleo])

	# Add the modeled error back in
	dais.paleo[i,] = dais.norm + rnorm(n.paleo, mean=0,sd=sqrt(var.dais))

	if(min(dais.tmp$Vm3)<Vmin) {ind.vmin[i]=1}

  setTxtProgressBar(pb, i)
}
close(pb)

## Post-calibrate out too low of Vais as well.
ind.survive.vmin = which(ind.vmin==0)
slr.out.good = slr.out.good[ind.survive.vmin,]
parameters.good = parameters[ind.survive.vmin,]
colnames(parameters.good) = parnames

## Trim down the model output
n.ensemble = length(ind.survive.vmin)
dais.paleo=dais.paleo[ind.survive.vmin,]

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

gsic.hind = gsic.hind[,ind.survive.vmin]
te.hind = te.hind[,ind.survive.vmin]
gis.hind = gis.hind[,ind.survive.vmin]
ais.hind = ais.hind[,ind.survive.vmin]
temp.hind = temp.hind[,ind.survive.vmin]
ocheat.hind = ocheat.hind[,ind.survive.vmin]

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
	filename.postcalibration = "../output_calibration/BRICK-fastdyn_postcalibratedParameters_22Aug2016.csv"
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
endyear = 2200; #if(!l.project & endyear>2009) print('l.project and endyear not compatible')
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
parameters.good = parameters.ensemble[ind.good,]

gsl.hind = gsl.hind[,ind.good]
gsic.hind = gsic.hind[,ind.good]
te.hind = te.hind[,ind.good]
gis.hind = gis.hind[,ind.good]
ais.hind = ais.hind[,ind.good]
temp.hind = temp.hind[,ind.good]
ocheat.hind = ocheat.hind[,ind.good]

## Fingerprints of sea-level rise sources on New Orleans local sea-level rise
fp.ais = 1.1
fp.gsic = 0.89
fp.gis = 0.81
fp.te = 1.0

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

names.output = c('slr','gsic','gis','ais','disint','te','temp','ocheat','slr.nofd','slr.nola','slr.nola.nofd')
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
	proj.rcp26$ais[i,] = brick.rcp26[[i]]$dais.out$Vais
	proj.rcp26$disint[i,] = brick.rcp26[[i]]$dais.out$Vdisint
	proj.rcp26$te[i,] = brick.rcp26[[i]]$te.out
	proj.rcp26$temp[i,] = brick.rcp26[[i]]$doeclim.out$temp + T0
	proj.rcp26$ocheat[i,] = brick.rcp26[[i]]$doeclim.out$ocheat + H0

	proj.rcp45$gsic[i,] = brick.rcp45[[i]]$gsic.out
	proj.rcp45$gis[i,] = brick.rcp45[[i]]$simple.out$sle.gis
	proj.rcp45$ais[i,] = brick.rcp45[[i]]$dais.out$Vais
	proj.rcp45$disint[i,] = brick.rcp45[[i]]$dais.out$Vdisint
	proj.rcp45$te[i,] = brick.rcp45[[i]]$te.out
	proj.rcp45$temp[i,] = brick.rcp45[[i]]$doeclim.out$temp + T0
	proj.rcp45$ocheat[i,] = brick.rcp45[[i]]$doeclim.out$ocheat + H0

	proj.rcp85$gsic[i,] = brick.rcp85[[i]]$gsic.out
	proj.rcp85$gis[i,] = brick.rcp85[[i]]$simple.out$sle.gis
	proj.rcp85$ais[i,] = brick.rcp85[[i]]$dais.out$Vais
	proj.rcp85$disint[i,] = brick.rcp85[[i]]$dais.out$Vdisint
	proj.rcp85$te[i,] = brick.rcp85[[i]]$te.out
	proj.rcp85$temp[i,] = brick.rcp85[[i]]$doeclim.out$temp + T0
	proj.rcp85$ocheat[i,] = brick.rcp85[[i]]$doeclim.out$ocheat + H0

	# Normalize the output to "ind.norm" (1961-1990? 1986-2005 (Mengel)?).
	# Normalize ocean heat uptake too, for sake of plotting (it will be plotted as
	# the amount of heat taken up by ocean since 1986-2005 period)
	proj.rcp26$gsic[i,] = proj.rcp26$gsic[i,] - mean( proj.rcp26$gsic[i,ind.norm])
	proj.rcp26$gis[i,] = proj.rcp26$gis[i,] - mean( proj.rcp26$gis[i,ind.norm])
	proj.rcp26$ais[i,] = proj.rcp26$ais[i,] - mean( proj.rcp26$ais[i,ind.norm])
	proj.rcp26$disint[i,] = proj.rcp26$disint[i,] - mean( proj.rcp26$disint[i,ind.norm])
	proj.rcp26$te[i,] = proj.rcp26$te[i,] - mean( proj.rcp26$te[i,ind.norm])
	proj.rcp26$temp[i,] = proj.rcp26$temp[i,] - mean( proj.rcp26$temp[i,ind.norm])
	proj.rcp26$ocheat[i,] = proj.rcp26$ocheat[i,] - mean( proj.rcp26$ocheat[i,ind.norm])

	proj.rcp45$gsic[i,] = proj.rcp45$gsic[i,] - mean( proj.rcp45$gsic[i,ind.norm])
	proj.rcp45$gis[i,] = proj.rcp45$gis[i,] - mean( proj.rcp45$gis[i,ind.norm])
	proj.rcp45$ais[i,] = proj.rcp45$ais[i,] - mean( proj.rcp45$ais[i,ind.norm])
	proj.rcp45$disint[i,] = proj.rcp45$disint[i,] - mean( proj.rcp45$disint[i,ind.norm])
	proj.rcp45$te[i,] = proj.rcp45$te[i,] - mean( proj.rcp45$te[i,ind.norm])
	proj.rcp45$temp[i,] = proj.rcp45$temp[i,] - mean( proj.rcp45$temp[i,ind.norm])
	proj.rcp45$ocheat[i,] = proj.rcp45$ocheat[i,] - mean( proj.rcp45$ocheat[i,ind.norm])

	proj.rcp85$gsic[i,] = proj.rcp85$gsic[i,] - mean( proj.rcp85$gsic[i,ind.norm])
	proj.rcp85$gis[i,] = proj.rcp85$gis[i,] - mean( proj.rcp85$gis[i,ind.norm])
	proj.rcp85$ais[i,] = proj.rcp85$ais[i,] - mean( proj.rcp85$ais[i,ind.norm])
	proj.rcp85$disint[i,] = proj.rcp85$disint[i,] - mean( proj.rcp85$disint[i,ind.norm])
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
	proj.rcp26$slr.nola[i,] = fp.gsic*proj.rcp26$gsic[i,] +
												    fp.gis*proj.rcp26$gis[i,] +
												    fp.ais*proj.rcp26$ais[i,] +
												    fp.te*proj.rcp26$te[i,]

	proj.rcp45$slr[i,] = proj.rcp45$gsic[i,] +
												proj.rcp45$gis[i,] +
												proj.rcp45$ais[i,] +
												proj.rcp45$te[i,]
	proj.rcp45$slr.nola[i,] = fp.gsic*proj.rcp45$gsic[i,] +
												    fp.gis*proj.rcp45$gis[i,] +
												    fp.ais*proj.rcp45$ais[i,] +
												    fp.te*proj.rcp45$te[i,]

	proj.rcp85$slr[i,] = proj.rcp85$gsic[i,] +
												proj.rcp85$gis[i,] +
												proj.rcp85$ais[i,] +
												proj.rcp85$te[i,]
	proj.rcp85$slr.nola[i,] = fp.gsic*proj.rcp85$gsic[i,] +
												    fp.gis*proj.rcp85$gis[i,] +
												    fp.ais*proj.rcp85$ais[i,] +
												    fp.te*proj.rcp85$te[i,]

	proj.rcp26$slr.nofd[i,] = proj.rcp26$slr[i,] - fp.ais*proj.rcp26$disint[i,]
	proj.rcp45$slr.nofd[i,] = proj.rcp45$slr[i,] - fp.ais*proj.rcp45$disint[i,]
	proj.rcp85$slr.nofd[i,] = proj.rcp85$slr[i,] - fp.ais*proj.rcp85$disint[i,]
	proj.rcp26$slr.nola.nofd[i,] = proj.rcp26$slr.nola[i,] - fp.ais*proj.rcp26$disint[i,]
	proj.rcp45$slr.nola.nofd[i,] = proj.rcp45$slr.nola[i,] - fp.ais*proj.rcp45$disint[i,]
	proj.rcp85$slr.nola.nofd[i,] = proj.rcp85$slr.nola[i,] - fp.ais*proj.rcp85$disint[i,]

	# And normalize sea-level rise
	proj.rcp26$slr[i,] = proj.rcp26$slr[i,] - mean(proj.rcp26$slr[i,ind.norm])
	proj.rcp45$slr[i,] = proj.rcp45$slr[i,] - mean(proj.rcp45$slr[i,ind.norm])
	proj.rcp85$slr[i,] = proj.rcp85$slr[i,] - mean(proj.rcp85$slr[i,ind.norm])
	proj.rcp26$slr.nola[i,] = proj.rcp26$slr.nola[i,] - mean(proj.rcp26$slr.nola[i,ind.norm])
	proj.rcp45$slr.nola[i,] = proj.rcp45$slr.nola[i,] - mean(proj.rcp45$slr.nola[i,ind.norm])
	proj.rcp85$slr.nola[i,] = proj.rcp85$slr.nola[i,] - mean(proj.rcp85$slr.nola[i,ind.norm])

	proj.rcp26$slr.nofd[i,] = proj.rcp26$slr.nofd[i,] - mean(proj.rcp26$slr.nofd[i,ind.norm])
	proj.rcp45$slr.nofd[i,] = proj.rcp45$slr.nofd[i,] - mean(proj.rcp45$slr.nofd[i,ind.norm])
	proj.rcp85$slr.nofd[i,] = proj.rcp85$slr.nofd[i,] - mean(proj.rcp85$slr.nofd[i,ind.norm])
  proj.rcp26$slr.nola.nofd[i,] = proj.rcp26$slr.nola.nofd[i,] - mean(proj.rcp26$slr.nola.nofd[i,ind.norm])
  proj.rcp45$slr.nola.nofd[i,] = proj.rcp45$slr.nola.nofd[i,] - mean(proj.rcp45$slr.nola.nofd[i,ind.norm])
  proj.rcp85$slr.nola.nofd[i,] = proj.rcp85$slr.nola.nofd[i,] - mean(proj.rcp85$slr.nola.nofd[i,ind.norm])

  setTxtProgressBar(pb, i)
}
close(pb)
##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## Write a netCDF ensemble output file including each of the RCP scenarios:
## (1) global total sea level, (2) global sea level without fast dynamics,
## (3) local (NOLA) sea level, (4) local sea level without fast dynamics
## Also will want each contribution to global sea level rise, for the hindcast
## plots

library(ncdf4)

dim.tproj <- ncdim_def('time_proj', 'years', as.double(t.proj))
dim.ensemble <- ncdim_def('ens', 'ensemble member', (1:n.ensemble))
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
gsl.nofd.rcp26 <- ncvar_def('GlobalSeaLevel_nofd_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  		longname = 'Global sea level without accounting for fast dynamics (RCP26)')
lsl.rcp26 <- ncvar_def('LocalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level accounting for fast dynamics (RCP26)')
lsl.nofd.rcp26 <- ncvar_def('LocalSeaLevel_nofd_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  		longname = 'Local sea level without accounting for fast dynamics (RCP26)')
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
gsl.nofd.rcp45 <- ncvar_def('GlobalSeaLevel_nofd_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  		longname = 'Global sea level without accounting for fast dynamics (RCP45)')
lsl.rcp45 <- ncvar_def('LocalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level accounting for fast dynamics (RCP45)')
lsl.nofd.rcp45 <- ncvar_def('LocalSeaLevel_nofd_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  		longname = 'Local sea level without accounting for fast dynamics (RCP45)')
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
gsl.nofd.rcp85 <- ncvar_def('GlobalSeaLevel_nofd_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  		longname = 'Global sea level without accounting for fast dynamics (RCP85)')
lsl.rcp85 <- ncvar_def('LocalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'Local sea level accounting for fast dynamics (RCP85)')
lsl.nofd.rcp85 <- ncvar_def('LocalSeaLevel_nofd_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  		longname = 'Local sea level without accounting for fast dynamics (RCP85)')
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
										list( gsl.rcp26, gsl.nofd.rcp26, lsl.rcp26, lsl.nofd.rcp26,
													gsl.rcp45, gsl.nofd.rcp45, lsl.rcp45, lsl.nofd.rcp45,
													gsl.rcp85, gsl.nofd.rcp85, lsl.rcp85, lsl.nofd.rcp85,
													gsic.rcp26, te.rcp26, gis.rcp26, ais.rcp26, temp.rcp26, ocheat.rcp26,
													gsic.rcp45, te.rcp45, gis.rcp45, ais.rcp45, temp.rcp45, ocheat.rcp45,
													gsic.rcp85, te.rcp85, gis.rcp85, ais.rcp85, temp.rcp85, ocheat.rcp85,
													gsl.hindcast, gsic.hindcast, te.hindcast, gis.hindcast, ais.hindcast, temp.hindcast, ocheat.hindcast,
													ais.paleo.05, ais.paleo.50, ais.paleo.95, ais.paleo.max, ais.paleo.min,
													ais.paleo.05.avg, ais.paleo.50.avg, ais.paleo.95.avg, ais.paleo.max.avg, ais.paleo.min.avg),
										force_v4 = TRUE)

ncvar_put(outnc, gsl.rcp26, t(proj.rcp26$slr))
ncvar_put(outnc, gsl.nofd.rcp26, t(proj.rcp26$slr.nofd))
ncvar_put(outnc, lsl.rcp26, t(proj.rcp26$slr.nola))
ncvar_put(outnc, lsl.nofd.rcp26, t(proj.rcp26$slr.nola.nofd))
ncvar_put(outnc, gsic.rcp26, t(proj.rcp26$gsic))
ncvar_put(outnc, te.rcp26, t(proj.rcp26$te))
ncvar_put(outnc, gis.rcp26, t(proj.rcp26$gis))
ncvar_put(outnc, ais.rcp26, t(proj.rcp26$ais))
ncvar_put(outnc, temp.rcp26, t(proj.rcp26$temp))
ncvar_put(outnc, ocheat.rcp26, t(proj.rcp26$ocheat))

ncvar_put(outnc, gsl.rcp45, t(proj.rcp45$slr))
ncvar_put(outnc, gsl.nofd.rcp45, t(proj.rcp45$slr.nofd))
ncvar_put(outnc, lsl.rcp45, t(proj.rcp45$slr.nola))
ncvar_put(outnc, lsl.nofd.rcp45, t(proj.rcp45$slr.nola.nofd))
ncvar_put(outnc, gsic.rcp45, t(proj.rcp45$gsic))
ncvar_put(outnc, te.rcp45, t(proj.rcp45$te))
ncvar_put(outnc, gis.rcp45, t(proj.rcp45$gis))
ncvar_put(outnc, ais.rcp45, t(proj.rcp45$ais))
ncvar_put(outnc, temp.rcp45, t(proj.rcp45$temp))
ncvar_put(outnc, ocheat.rcp45, t(proj.rcp45$ocheat))

ncvar_put(outnc, gsl.rcp85, t(proj.rcp85$slr))
ncvar_put(outnc, gsl.nofd.rcp85, t(proj.rcp85$slr.nofd))
ncvar_put(outnc, lsl.rcp85, t(proj.rcp85$slr.nola))
ncvar_put(outnc, lsl.nofd.rcp85, t(proj.rcp85$slr.nola.nofd))
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





##==============================================================================
##==============================================================================
## Pick up here?
## To run RCP2.6 or 4.5 ensembles through Van Dantzig analysis, change "RCP85"
## in the "sea_level" and "sea_level_nofd" lines below to your RCP of choice.
if(FALSE){
	setwd('~/codes/BRICK/calibration')
	library(ncdf4)
	filename.in = "../output_model/BRICK-fastdyn_physical_gamma_29Jan2017.nc"
	ncdata <- nc_open(filename.in)
	sea_level = ncvar_get(ncdata, 'LocalSeaLevel_RCP45')
	sea_level_nofd = ncvar_get(ncdata, 'LocalSeaLevel_nofd_RCP45')
	mod.time =ncvar_get(ncdata, 'time_proj')
	nc_close(ncdata)
	today=Sys.Date(); today=format(today,format="%d%b%Y")
	filename.vdout = paste('../output_model/vanDantzig_RCP45_gamma_',today,'.nc',sep="")
} else {
	sea_level=t(proj.rcp85$slr.nola)
	sea_level_nofd=t(proj.rcp85$slr.nola.nofd)
}

## Evaluate van Dantzig model for each of these realizations

## Set up
currentyear		= 2015		# start year
endyear			= 2100		# reevaluation year
lowleveeheight	= 0			# lowest levee heightening considered (m)
highleveeheight	= 10		# largest levee heightening considered (m)
increment		= 0.05		# increment of levee heightening range (m)

## Time horizon until next evaluation of levees (years)
T = endyear - currentyear
time_frame=mod.time

## Range of considered levee heights (meters)
X = seq(from = lowleveeheight, to = highleveeheight, by = increment)

## Create a matrix of sea-level data starting from the current year for each model simulation.
## Note that if you included the NOLA fingerprinting above (in by default), then
## what is called "global_sea_level" below is actually local, but does not include
## subsidence.
global_sea_level = data.matrix(sea_level[match(currentyear, time_frame):match(endyear, time_frame),])
global_sea_level_nofd = data.matrix(sea_level_nofd[match(currentyear, time_frame):match(endyear, time_frame),])

## Set up basic information for the van Dantzig model.

## Number of observations (ensemble members):
n_obs = length(global_sea_level[1, ])

## Conversion from feet to meters:
feet_to_meter = function(ft) { ft * 0.3048 }

## Investment costs for levee heightening. Table # 2 in Jonkman et al. (2009)
design_surge_level_ft = c(9, 11, 13, 15, 17, 21, 25)
design_surge_level_M = feet_to_meter(design_surge_level_ft)

investments_US = c(2.2E09, 2.4E09, 2.6E09, 2.9E09, 3.1E09, 3.6E09, 4.1E09)

I_table = matrix(c(design_surge_level_M, investments_US), ncol=2)
colnames(I_table) = c("Levee heightening", "Associated cost")

p_zero_p = 0.0038                 # Initial flood frequency (1/yr) with zero height increase (Van Dantzig (1956))
alpha_p = 2.6                     # Exponential flood frequency constant (Van Dantzig (1956))
V_p = c(5e+9, 3e+10)              # Value of goods protected by dike (based on estimates in Jonkman et al. (2009)) (US$)
delta_prime_p = c(0.02, 0.06)     # Discount rate (percent/year) (based on estimates in Jonkman et al. (2009))
I_range = c(-0.5,1.0)             # Investment cost uncertainty range (as fraction of investment cost) (Jonkman and Dutch Perspective use -50%, +100%)
subs_rate = 0.0056                # Rate of land subsidence (meter/year) (Dixon et al. (2006))

investment.fit <- lm(I_table[,2]~I_table[,1])
intercept.h2i = investment.fit$coefficients[[1]]
slope.h2i = investment.fit$coefficients[[2]]

source("../R/VD_NOLA.R")

## Sample parameters (using Dutch Perspective as a guide)
## Uses Latin Hypercube to sample Van Dantzig parameters, assigns each ensemble
## member of SLR a set of parameters.
params.vd = parameter_sampling_DP(n_obs, p_zero_p, alpha_p, V_p, delta_prime_p, I_range, subs_rate)

## Test simulation to get the names
vanDantzig.out = VD_NOLA_R(params.vd[1,], I_table, T, X, global_sea_level[,1])

n.ensemble = nrow(params.vd)
vanDantzig.ensemble = vector("list", ncol(vanDantzig.out))
vanDantzig.nofd.ensemble = vector("list", ncol(vanDantzig.out))
names(vanDantzig.ensemble) = names(vanDantzig.out)
names(vanDantzig.nofd.ensemble) = names(vanDantzig.out)
for (i in 1:ncol(vanDantzig.out)) {
  vanDantzig.ensemble[[i]] = mat.or.vec(nrow(vanDantzig.out), n.ensemble)
  vanDantzig.nofd.ensemble[[i]] = mat.or.vec(nrow(vanDantzig.out), n.ensemble)
}

## Run the simulations
pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3)
for (i in 1:n.ensemble) {
  vanDantzig.out = VD_NOLA_R(params.vd[i,], I_table, T, X, global_sea_level[,i])
  for (j in 1:ncol(vanDantzig.out)) {
		vanDantzig.ensemble[[j]][,i] = vanDantzig.out[,j]
	}
  vanDantzig.out = VD_NOLA_R(params.vd[i,], I_table, T, X, global_sea_level_nofd[,i])
  for (j in 1:ncol(vanDantzig.out)) {
		vanDantzig.nofd.ensemble[[j]][,i] = vanDantzig.out[,j]
	}
  setTxtProgressBar(pb, i)
}
close(pb)

##==============================================================================
##==============================================================================





##==============================================================================
##==============================================================================
## Output results of van Dantzig to netCDF file

library(ncdf4)

dim.heightening <- ncdim_def('H', 'meters', as.double(X))
dim.ensemble <- ncdim_def('ens', 'ensemble member', (1:n.ensemble))

cost <- ncvar_def('ExpectedCost', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected cost accounting for fast dynamics')
cost.nofd <- ncvar_def('ExpectedCost_nofd', 'US$', list(dim.heightening, dim.ensemble), -999,
                  		longname = 'Expected cost without accounting for fast dynamics')
loss <- ncvar_def('ExpectedLoss', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected loss accounting for fast dynamics')
loss.nofd <- ncvar_def('ExpectedLoss_nofd', 'US$', list(dim.heightening, dim.ensemble), -999,
                  		longname = 'Expected loss without accounting for fast dynamics')
investment <- ncvar_def('ExpectedInvestment', 'US$', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected investment accounting for fast dynamics')
investment.nofd <- ncvar_def('ExpectedInvestment_nofd', 'US$', list(dim.heightening, dim.ensemble), -999,
                  		longname = 'Expected investment without accounting for fast dynamics')
preturn <- ncvar_def('ExpectedPreturn', 'years', list(dim.heightening, dim.ensemble), -999,
                  longname = 'Expected return period accounting for fast dynamics')
preturn.nofd <- ncvar_def('ExpectedPreturn_nofd', 'years', list(dim.heightening, dim.ensemble), -999,
                  		longname = 'Expected return period without accounting for fast dynamics')

today=Sys.Date(); today=format(today,format="%d%b%Y")

outnc <- nc_create(filename.vdout,
										list(cost, cost.nofd, loss, loss.nofd, investment, investment.nofd, preturn, preturn.nofd),
										force_v4 = TRUE)

ncvar_put(outnc, cost, vanDantzig.ensemble$Expected_costs)
ncvar_put(outnc, cost.nofd, vanDantzig.nofd.ensemble$Expected_costs)
ncvar_put(outnc, loss, vanDantzig.ensemble$Expected_loss)
ncvar_put(outnc, loss.nofd, vanDantzig.nofd.ensemble$Expected_loss)
ncvar_put(outnc, investment, vanDantzig.ensemble$Investment)
ncvar_put(outnc, investment.nofd, vanDantzig.nofd.ensemble$Investment)
ncvar_put(outnc, preturn, 1/vanDantzig.ensemble$Average_p_fail)
ncvar_put(outnc, preturn.nofd, 1/vanDantzig.nofd.ensemble$Average_p_fail)

nc_close(outnc)

##==============================================================================
##==============================================================================

t.end = proc.time()

print(paste('it took ',(t.end[3]-t.beg[3])/60,' minutes to process an initial ensemble of ',n.ensemble.report,' to a fully calibrated ensemble of ',n.ensemble,sep=''))


##==============================================================================
## End
##==============================================================================
