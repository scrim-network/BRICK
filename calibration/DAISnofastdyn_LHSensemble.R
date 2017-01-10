##==============================================================================
## Script to do a Latin hypercube sample of DAIS (no fast dynamics) parameters
## and pre-calibrate against same likelihood function as the fast dynamics DAIS
## version was MCMC-calibrated against. LHS precalibration is used because the
## MCMC fails to converge without the fast dynamics mechanism.
##
## Set the BRICK calibrated parameters file below, as well as the rho.simple.fixed
## option. Each LHS sample of DAIS parameters that fits the precalibration
## likelihood (step function for LIG window and heaviside function for total
## sea-level rise versus AIS component of sea-level rise) is combined with a set
## of calibrated rest-of-model (BRICK) parameters.
##
## Questions? Tony Wong <twong@psu.edu>
##==============================================================================

rm(list =ls()) #Clear global environment

filename.BRICKcalibration = "../output_calibration/BRICK_calibratedParameters_12Aug2016.csv"
filename.rho_simple_fixed = "../output_calibration/rho_simple_fixed_06Sep2016.csv"

today=Sys.Date(); today=format(today,format="%d%b%Y")
filename.daisout = paste('../output_model/DAISfastdyn_noFD-paleo-ensemble_',today,'.nc',sep="")

n.lhs = 200

##==============================================================================
t.beg=proc.time()

## Setup packages and libraries
#install.packages('compiler')
#install.packages('pscl')
library(pscl) # install inverse gamma distribution
library(compiler)
enableJIT(3)
enableJIT(3)

## Set the seed (for reproducibility)
set.seed(1234)

## Read the data forcing for hindcasts and projections. Yields:
##  Ta (Antarctic temperature reduced to sea-level)
##  Toc (High latitude subsurface ocean temperature)
##  SL, obs.sl (Reconstructed sea level, Church and White (2011) modern-era sea level)
##  dSL (time rate of change of sea level)
##  date (240 Kyr before present to 2100 AD at 1-year intervals of forcings)
l.project=FALSE
source('../calibration/DAIS_readData.R')

## Set Best Case (Case #4) from Shaffer (2014), and parameter ranges
parnames    = c('anto.a','anto.b','gamma','alpha.dais','mu'  ,'nu'  ,'P0' ,'kappa.dais','f0' ,'h0'  ,'c'  , 'b0','slope' ,'var.dais')
parameters0 = c(  0.1574, 0.6677 ,  2    , 0.35       , 8.7  , 0.012, 0.35, 0.04       , 1.2 , 1471 , 95  , 775 , 0.0006 , 0.0004656)
bound.lower = c( 0.0    , 0      ,  0.5  , 0          , 7.05 , 0.003,0.026, 0.025      , 0.6 , 735.5, 47.5, 740 , 0.00045, 0        )
bound.upper = c( 1.0    , 2      ,  4.25 , 1          , 13.65, 0.015, 1.5 , 0.085      , 1.8 ,2206.5,142.5, 820 , 0.00075, 2        )

alpha.var = 2     # alpha parameter for inverse gamma for var (E[x]=beta/(alpha+1))
beta.var = 1      # beta parameter for inverse gamma for var (uncertainty parameter)
                  # note that the upper bound on var.dais is not actually imposed; just for plotting

## Source the DAIS model
source('../fortran/R/daisantoF.R')
##==============================================================================

## Get a standard hindcast and future projections (a "control" model)
AIS_melt = daisantoF(
                anto.a=parameters0[match("anto.a",parnames)],
                anto.b=parameters0[match("anto.b",parnames)],
                gamma=parameters0[match("gamma",parnames)],
                alpha=parameters0[match("alpha.dais",parnames)],
                mu=parameters0[match("mu",parnames)],
                nu=parameters0[match("nu",parnames)],
                P0=parameters0[match("P0",parnames)],
                kappa=parameters0[match("kappa.dais",parnames)],
                f0=parameters0[match("f0",parnames)],
                h0=parameters0[match("h0",parnames)],
                c=parameters0[match("c",parnames)],
                b0=parameters0[match("b0",parnames)],
                slope=parameters0[match("slope",parnames)],
                Tg=Tg.recon,
                slope.Ta2Tg=slope.Ta2Tg,
                intercept.Ta2Tg=intercept.Ta2Tg,
                SL=SL,
                dSL=dSL)

##==============================================================================
## Set up (pre-)calibration windows around the data
## The point: We run the model at many parameter values and see which ones
##            send the simulation through a window (width determined by the
##            observational errors) around the observational data.
## These windows are presented in Shaffer (2014) and Shepherd et al. (2012)
## 1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
## We want the cumulative sea-level equivalent in meters for the year 2002
## Note the conversion: 360Gt = 1mm SLE
## A fifth window is added to match IPCC AR5 Ch13 (page 1151) AIS SLR trend:
## 0.27 +/- 0.11 mm/year (convert to m/year here)

## Precal windows 5-?:
## Last "precalibration window" is 1993-2010 mean trend, from the IPCC AR5 Ch13
## (Page 1151), for AIS SLR contribution: 0.27 +- 0.11 mm/year
## Note that model output is in meters SLE and these trends are mm, so a
## conversion is necessary.

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
upper.wind = c(7.4, -6.9, -1.25, positive_2SE) # Windows 2-3 from Kelsey, Window 1 from DeConto and Pollard 2016
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

##==============================================================================




##==============================================================================
## Preliminary Latin Hypercube to find decent starting parameter values

require(lhs)

# Use the same statistical model - calculate log-likelihood, pick starting values
# to maximize this.
source('../calibration/DAIS_assimLikelihood.R')

t0=proc.time()
# Draw LHS sample
#n.lhs = 2000 #(set above, with the file names)
parameters.lhs <- randomLHS(n.lhs, length(parnames))

# Transform unif(0,1) to the parameter bounds
for (i in 1:length(parnames)) {
  parameters.lhs[,i] <- qunif(parameters.lhs[,i], bound.lower[i], bound.upper[i])
}
colnames(parameters.lhs)=parnames
t1=proc.time()

lpost.lhs = rep(NA,n.lhs)

pb <- txtProgressBar(min=0,max=n.lhs,initial=0,style=3)
for (j in 1:n.lhs) {
  lpost.lhs[j] = log.post(  parameters.in = as.numeric(parameters.lhs[j,]),
                            parnames.in = parnames,
                            bound.lower.in = bound.lower,
                            bound.upper.in = bound.upper,
                            obs.in = obs.targets,
                            obs.err.in = obs.err,
                            obs.step.in = obs.years,
                            trends.ais.in = trends.ais,
                            trends.err.in = trends.err,
                            ind.trends.in = ind.trends,
                            ind.norm.in = ind.relative,
                            alpha.var = alpha.var,
                            beta.var = beta.var,
                            #shape.lambda = shape.lambda,
                            #rate.lambda = rate.lambda,
                            #shape.Tcrit = shape.Tcrit,
                            #rate.Tcrit = rate.Tcrit,
                            slope.Ta2Tg.in = slope.Ta2Tg,
                            intercept.Ta2Tg.in = intercept.Ta2Tg,
                            Tg.in = Tg.recon,
                            SL.in = SL,
                            dSL.in = dSL
                            )
  setTxtProgressBar(pb, j)
}
close(pb)
t2=proc.time()

##==============================================================================

igood=which(!is.infinite(lpost.lhs))
parameters.dais = parameters.lhs[igood,]
parnames.dais   = colnames(parameters.dais)

lpost=lpost.lhs[igood]

## Get the indices WITHIN PARAMETERS.DAIS of the
##  n.best highest likelihood simulations
n.best = round(.1*length(igood))
ibest = rev(order(lpost))[1:n.best]
parameters = parameters.dais[ibest,]

##==============================================================================
##==============================================================================

n.ensemble = nrow(parameters)
n.parameters = ncol(parameters)

source('../fortran/R/daisantoF.R')
source('../calibration/DAIS_readData.R')

n.paleo = length(SL)
date = seq(-239999,16,1) #240 Kyr BP to 2100AD at one year intervals of the forcings
norm.period=c(1961,1990)
ibeg=which(date==(norm.period[1]-2000))
iend=which(date==(norm.period[2]-2000))
ind.norm.paleo=ibeg:iend
t.paleo = date

## Make smoothed version of the AIS paleo results, so the plots are not massive
n.avg=100	# number of years in averaging period
n.time.avg=ceiling(length(date)/n.avg)
dais.paleo.05.avg = rep(NA,n.time.avg)
dais.paleo.50.avg = rep(NA,n.time.avg)
dais.paleo.95.avg = rep(NA,n.time.avg)
dais.paleo.max.avg = rep(NA,n.time.avg)
dais.paleo.min.avg = rep(NA,n.time.avg)
date.avg=seq(date[1],date[length(date)],by=n.avg)

## Only save the smoothed ones
dais.paleo.avg = mat.or.vec(n.ensemble, n.time.avg)
dais.smooth = rep(0,n.time.avg)

pb <- txtProgressBar(min=0,max=n.ensemble,initial=0,style=3);
for (i in 1:n.ensemble) {

	anto.a=parameters[i,match("anto.a",parnames)]
  anto.b=parameters[i,match("anto.b",parnames)]
  gamma =parameters[i,match("gamma" ,parnames)]
  alpha =parameters[i,match("alpha.dais" ,parnames)]
  mu =parameters[i,match("mu" ,parnames)]
  nu =parameters[i,match("nu" ,parnames)]
  P0 =parameters[i,match("P0" ,parnames)]
  kappa =parameters[i,match("kappa.dais" ,parnames)]
  f0 =parameters[i,match("f0" ,parnames)]
  h0 =parameters[i,match("h0" ,parnames)]
  c =parameters[i,match("c" ,parnames)]
  b0 =parameters[i,match("b0" ,parnames)]
  slope =parameters[i,match("slope" ,parnames)]
  var.dais =parameters[i,match("var.dais" ,parnames)]

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
	dais.norm = dais.norm + rnorm(n.paleo, mean=0,sd=sqrt(var.dais))

  # smooth
  dais.smooth[1:(n.time.avg-1)] = apply( matrix(dais.norm[1:(floor(length(dais.norm)/n.avg)*n.avg)], nrow=n.avg) ,2,'mean')
  dais.smooth[n.time.avg] = mean(dais.norm[(floor(length(dais.norm)/n.avg)*n.avg+1):n.paleo])
  dais.paleo.avg[i,] = dais.smooth

  setTxtProgressBar(pb, i)
}
close(pb)

## Get 5-95% CI for hindcasts

## Source a useful script, to allow for assigning multiple outputs at once
source('../Useful/MultipleOutput.R') # defines the ":=" operator

## Actually tally up the data
pb <- txtProgressBar(min=0,max=length(date.avg),initial=0,style=3);
for (t in 1:length(date.avg)){
	c(dais.paleo.05.avg[t] , dais.paleo.50.avg[t] , dais.paleo.95.avg[t] , dais.paleo.max.avg[t], dais.paleo.min.avg[t]) := quantile(dais.paleo.avg[,t],c(0.05,.50,.95,1,0), na.rm=TRUE)
  setTxtProgressBar(pb, t)
}
close(pb)


library(ncdf4)

dim.ensemble <- ncdim_def('ens', 'ensemble member', (1:n.ensemble))
dim.tpaleo.avg <- ncdim_def('time_paleo_avg', 'year avg paleo', as.double(date.avg))

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

outnc <- nc_create(filename.daisout,
										list( ais.paleo.05.avg, ais.paleo.50.avg, ais.paleo.95.avg, ais.paleo.max.avg, ais.paleo.min.avg),
										force_v4 = TRUE)

ncvar_put(outnc, ais.paleo.05.avg, dais.paleo.05.avg)
ncvar_put(outnc, ais.paleo.50.avg, dais.paleo.50.avg)
ncvar_put(outnc, ais.paleo.95.avg, dais.paleo.95.avg)
ncvar_put(outnc, ais.paleo.max.avg, dais.paleo.max.avg)
ncvar_put(outnc, ais.paleo.min.avg, dais.paleo.min.avg)

nc_close(outnc)

##==============================================================================

t.end=proc.time()

##==============================================================================
##==============================================================================



##==============================================================================
## End
##==============================================================================
