##==============================================================================
## Script to calibrate Antarctic Ocean temperature emulator for use within
## BRICK coupled modeling framework.
##
## Parameters to calibrate, and standard values (Shaffer, 2014):
## [1] a.anto
## [2] b.anto
##
## Questions? Tony Wong <twong@psu.edu>
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

rm(list =ls()) #Clear global environment

## Setup packages and libraries
#install.packages('compiler')
#install.packages('pscl')
library(pscl) # install inverse gamma distribution
library(compiler)
enableJIT(3)

## Set the seed (for reproducibility)
set.seed(1234)

## Read the data forcing for hindcasts and projections. Yields:
##  Ta (Antarctic temperature reduced to sea-level)
##  Toc (High latitude subsurface ocean temperature)
##  SL, obs.sl (Reconstructed sea level, Church and White (2011) modern-era sea level)
##  dSL (time rate of change of sea level)
##  date (240 Kyr before present to 2100 AD at 1-year intervals of forcings)
source('../calibration/DAIS_readData.R')

## To run BRICK in a full-forward capacity, need everything relative to 1850
## (the starting point)
date.ce <- date+2000
Toc.norm <- Toc[which(date.ce==Tg.time[1]):which(date.ce==Tg.time[length(Tg.time)])]
Toc.norm <- Toc.norm - Toc.norm[1]
Tg.norm <- Tg - Tg[1]


## Source the DAIS/ANTO model
source('../fortran/R/daisantoF.R')

## Set up the parameters to calibrate
parnames    = c('anto.a','anto.b')
parameters0 = c(  0.1574, 0.6677 )
bound.lower = c( 0.0    , 0      )
bound.upper = c( 1.0    , 2      )
##==============================================================================
rmse = function(parameters.in,
                parnames.in,
                Toc.obs.in,
                Tg.in
) {
    a.anto <- parameters.in[match("anto.a",parnames.in)]
    b.anto <- parameters.in[match("anto.b",parnames.in)]
    mod <- anto(a=a.anto, b=b.anto, Tf=-1.8, Tg=Tg.in)
    rmse <- sqrt( mean( (mod-Toc.obs.in)^2 ))

    rmse
}
##==============================================================================

Tg.recon.norm <- Tg.recon - mean(Tg.recon[ind.relative])
Toc.recon.norm <- Toc - mean(Toc[ind.relative])
Toc.anto.norm <- anto(a=.25, b=.05, Tf=-1.8, Tg=Tg.recon.norm)
Toc.anto.norm <- Toc.anto.norm - mean(Toc.anto.norm[ind.relative])
plot(Tg.recon.norm, type='l')
lines(Toc.recon.norm, type='l', col='blue')
lines(Toc.anto.norm, type='l', col='red')

##==============================================================================
## Preliminary Latin Hypercube to find decent starting parameter values

require(lhs)

t0=proc.time()
# Draw LHS sample
n.lhs = 1000
parameters.lhs <- randomLHS(n.lhs, length(parnames))

# Transform unif(0,1) to the parameter bounds
for (i in 1:length(parnames)) {
  parameters.lhs[,i] <- qunif(parameters.lhs[,i], bound.lower[i], bound.upper[i])
}
colnames(parameters.lhs)=parnames
t1=proc.time()

rmse.lhs = rep(NA,n.lhs)

pb <- txtProgressBar(min=0,max=n.lhs,initial=0,style=3)
for (j in 1:n.lhs) {
  rmse.lhs[j] = rmse(   parameters.in = as.numeric(parameters.lhs[j,]),
                        parnames.in = parnames,
                        Toc.obs.in = Toc.recon.norm,
                        Tg.in = Tg.recon.norm
                        )
  setTxtProgressBar(pb, j)
}
close(pb)
t2=proc.time()

tmp <- data.frame( cbind(rmse.lhs,parameters.lhs))
tmp.sort <- tmp[order(rmse.lhs),]






## Get an updated model to compare against control
AIS_melt0.lhs = daisantoF(
                anto.a=parameters0.lhs[match("anto.a",parnames)],
                anto.b=parameters0.lhs[match("anto.b",parnames)],
                gamma = parameters0.lhs[match("gamma",parnames)],
                alpha = parameters0.lhs[match("alpha.dais",parnames)],
                mu    = parameters0.lhs[match("mu",parnames)],
                nu    = parameters0.lhs[match("nu",parnames)],
                P0    = parameters0.lhs[match("P0",parnames)],
                kappa = parameters0.lhs[match("kappa.dais",parnames)],
                f0    = parameters0.lhs[match("f0",parnames)],
                h0    = parameters0.lhs[match("h0",parnames)],
                c     = parameters0.lhs[match("c",parnames)],
                b0    = parameters0.lhs[match("b0",parnames)],
                slope = parameters0.lhs[match("slope",parnames)],
#                Ta    = Ta,
#                Toc   = Toc,
                Tg    = Tg.recon,
                slope.Ta2Tg = slope.Ta2Tg,
                intercept.Ta2Tg = intercept.Ta2Tg,
                SL    = SL,
                dSL   = dSL,
                includes_dSLais=1
                )
##==============================================================================



##==============================================================================
## MCMC simulation
##==============================================================================

## Source the statistical models
source('../calibration/DAIS_assimLikelihood.R')

#install.packages('adaptMCMC')
library(adaptMCMC)
accept.mcmc = 0.234										# Optimal as # parameters->infinity
																			#	(Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = 5e5											# number of iterations for MCMC
gamma.mcmc = 0.5											# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)				# remove first ?? of chains for burn-in
stopadapt.mcmc = round(niter.mcmc*1.0)# stop adapting after ?? iterations? (niter*1 => don't stop)
step.mcmc = (bound.upper-bound.lower)*.05
##==============================================================================
## Actually run the calibration

if(FALSE){
t.beg=proc.time()											# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out1 = MCMC(log.post, niter.mcmc, parameters0.lhs, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames    , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  obs.in=obs.targets      , obs.err.in=obs.err        , obs.step.in=obs.years     ,
                  trends.ais.in=trends.ais, trends.err.in=trends.err  , ind.trends.in=ind.trends  ,
                  ind.norm=ind.relative   , alpha.var=alpha.var       , beta.var=beta.var         ,
                  slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , Tg.in=Tg.recon,
                  #Ta.in=Ta  , Toc.in=Toc ,
                  SL.in=SL  , dSL.in=dSL)
t.end=proc.time()											# save timing
chain1 = amcmc.out1$samples
}

## Extend and run more MCMC samples?
if(FALSE){
t.beg=proc.time()
amcmc.extend1 = MCMC.add.samples(amcmc.out1, niter.mcmc,
									parnames.in=parnames    , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  obs.in=obs.targets      , obs.err.in=obs.err        , obs.step.in=obs.years     ,
                  trends.ais.in=trends.ais, trends.err.in=trends.err  , ind.trends.in=ind.trends  ,
                  ind.norm=ind.relative   , alpha.var=alpha.var       , beta.var=beta.var         ,
                  slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , Tg.in=Tg.recon,
                  #Ta.in=Ta  , Toc.in=Toc ,
                  SL.in=SL  , dSL.in=dSL)
t.end=proc.time()											# save timing
chain1 = amcmc.extend1$samples
}

## If you want to run 2 (or more) chains in parallel (save time, more sampling)
## Quick note about using MCMC.parallel and MCMC.add.samples:
##    You can add samples to parallel chains by simply replacing "amcmc.out1"
##    with "amcmc.par1" in the MCMC.add.samples command above BUT these R
##    functions do not have the capability to add samples in parallel. It does
##    the task, but does it in serial after the initial chains were in parallel.
if(TRUE){
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.par1 = MCMC.parallel(log.post, niter.mcmc, parameters0.lhs, n.chain=4, n.cpu=4,
                  dyn.libs='../fortran/dais.so', scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames    , bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  obs.in=obs.targets      , obs.err.in=obs.err        , obs.step.in=obs.years     ,
                  trends.ais.in=trends.ais, trends.err.in=trends.err  , ind.trends.in=ind.trends  ,
                  ind.norm=ind.relative   , alpha.var=alpha.var       , beta.var=beta.var         ,
                  slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg , Tg.in=Tg.recon,
                  #Ta.in=Ta  , Toc.in=Toc ,
                  SL.in=SL  , dSL.in=dSL)
t.end=proc.time()											# save timing
chain1=amcmc.par1[[1]]$samples
chain2=amcmc.par1[[2]]$samples
chain3=amcmc.par1[[3]]$samples
chain4=amcmc.par1[[4]]$samples
}

## Save workspace image - you do not want to re-simulate all those!
save.image(file = "DAIS-noFD_calib_MCMC.RData")

##==============================================================================
if(FALSE){
## Check #1: History plots
par(mfrow=c(4,4))
for (pp in 1:length(parnames)) {
	plot(chain1[,pp], type="l", ylab=parnames[pp], xlab="Number of Runs", main="")
}

## Check #2: Heidelberger and Welch's convergence diagnostic:
heidel.diag(chain1, eps=0.1, pvalue=0.05)

## Check #3: Gelman and Rubin's convergence diagnostic:
# Converged when the potental scale reduction factor is less than 1.1
set.seed(222)
t.beg=proc.time()										# save timing (running millions of iterations so best to have SOME idea...)
amcmc.out2 = MCMC(log.post, niter.mcmc, parameters0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,
									gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc),
									parnames.in=parnames, bound.lower.in=bound.lower, bound.upper.in=bound.upper,
                  obs.in=obs.targets  , obs.err.in=obs.err        , obs.step.in=obs.years     ,
                  alpha.var=alpha.var , beta.var=beta.var         ,
                  slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg ,
                  Ta.in=Tg.recon      , SL.in=SL )
chain2 = amcmc.out2$samples
t.end=proc.time()											# save timing

mcmc1 = as.mcmc(chain1)
mcmc2 = as.mcmc(chain2)

mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2))
gelman.diag(mcmc_chain_list)

# Reset seed back to original
set.seed(1234)

## Histograms
par(mfrow=c(4,4))
for (pp in 1:length(parnames)) {
	hist(chain1[burnin:niter.mcmc,pp], xlab=parnames[pp])
}
}
##==============================================================================

## Determine when (in increments of 50,000 iterations, using Gelman and Rubin
## diagnostic) the two parallel chains, which are assumed to be
## chain1=amcmc.par1[[1]]$samples and chain2=amcmc.par1[[2]]$samples
## as defined above in the MCMC.parallel(...) command.

## Initialize the testing of the Gelman and Rubin diagnostics
niter.test = seq(from=50000, to=niter.mcmc, by=10000)
gr.stat = rep(NA,length(niter.test))

## Calculate the statistic at a few spots throughout the chain. Once it is
## close to 1 (people often use GR<1.1 or 1.05), the between-chain variability
## is indistinguishable from the within-chain variability, and they are
## converged. It is only after the chains are converged that you should use the
## parameter values as posterior draws, for analysis.
for (i in 1:length(niter.test)){
  mcmc1 = as.mcmc(chain1[1:niter.test[i],])
  mcmc2 = as.mcmc(chain2[1:niter.test[i],])
  mcmc_chain_list = mcmc.list(list(mcmc1, mcmc2))
  gr.stat[i] = gelman.diag(mcmc_chain_list)[2]
}

## Plot GR statistics as a function of iterations, decide where to cut off
## chains and use the tails of both for analysis
plot(niter.test,gr.stat)

##==============================================================================

## GR statistic gets and stays below 1.1 by iteration ...? (set this to n.burnin)
n.burnin = 120000
n.sample = nrow(chain1)-n.burnin
parameters1=chain1[(n.burnin+1):nrow(chain1),]
parameters2=chain2[(n.burnin+1):nrow(chain1),]
parameters3=chain3[(n.burnin+1):nrow(chain1),]
parameters4=chain4[(n.burnin+1):nrow(chain1),]
parameters.posterior = rbind(parameters1,parameters2, parameters3, parameters4)
n.parameters = ncol(parameters.posterior)

## Fit PDFs to the parameter distributions
pdf.all=vector('list',n.parameters)
n.node=100
for (pp in 1:n.parameters){
  tmp = density(parameters.posterior[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.all[[pp]] = tmp; names(pdf.all)[pp]=parnames[pp]
}

##==============================================================================
## Plot the PDFs?
par(mfrow=c(4,4))
for (pp in 1:n.parameters){
  plot(pdf.all[[pp]]$x,pdf.all[[pp]]$y,type='l',xlab=parnames[pp],ylab='density');
}

if(FALSE){
## Plot all relative to older estimates?
dat.parameters = read.csv("../output_calibration/DAIS_calibratedParameters_29Jul2016.csv")
parameters.old = dat.parameters[1:(nrow(dat.parameters)-1), ]
n.parameters.old=ncol(parameters.old)
pdf.old=vector('list',n.parameters.old)
n.node=200
for (pp in 1:n.parameters.old){
  tmp = density(parameters.old[,pp],kernel='gaussian',
                n=n.node,from=bound.lower[pp],to=bound.upper[pp])
  pdf.old[[pp]] = tmp; names(pdf.old)[pp]=parnames[pp]
}
par(mfrow=c(4,4))
for (pp in 1:n.parameters){
	plot(pdf.all[[pp]]$x,pdf.all[[pp]]$y,type='l',xlab=parnames[pp],col='red');
  	lines(pdf.old[[pp]]$x,pdf.old[[pp]]$y,type='l',col='black');
}
}
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
