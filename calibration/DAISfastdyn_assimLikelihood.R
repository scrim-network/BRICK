##==============================================================================
## R Function: DAISfastdyn_assimLikelihood.R
## -Antarctic Ice Sheet (AIS) model
##
## compute (log) likelihood for observations
## The observations are independent and identically distributed (IID)
##
## -Original authors: Yawen Guan and Kelsey Ruckert (klr324@psu.edu)
## - Modified by Tony Wong (twong@psu.edu) for running within BRICK model framework
##==============================================================================
## -June 10 2015 (original)
## - 6 July 2016 (modified)
##==============================================================================

##==============================================================================
## Log-likelihood
log.lik = function( parameters.in,
                    parnames.in,
                    obs.in,
                    obs.err.in,
                    obs.step.in,
                    ind.norm.in,
                    SL.in,
                    dSL.in,
                    trends.ais.in,
                    trends.err.in,
                    ind.trends.in,
                    slope.Ta2Tg.in=1,
                    intercept.Ta2Tg.in=0,
                    Tg.in=NULL,
                    Ta.in=NULL,
                    Toc.in=NULL
                    )
{

  anto.a=parameters.in[match("anto.a",parnames.in)]
  anto.b=parameters.in[match("anto.b",parnames.in)]
  gamma =parameters.in[match("gamma" ,parnames.in)]
  alpha =parameters.in[match("alpha.dais" ,parnames.in)]
  mu =parameters.in[match("mu" ,parnames.in)]
  nu =parameters.in[match("nu" ,parnames.in)]
  P0 =parameters.in[match("P0" ,parnames.in)]
  kappa =parameters.in[match("kappa.dais" ,parnames.in)]
  f0 =parameters.in[match("f0" ,parnames.in)]
  h0 =parameters.in[match("h0" ,parnames.in)]
  c =parameters.in[match("c" ,parnames.in)]
  b0 =parameters.in[match("b0" ,parnames.in)]
  slope =parameters.in[match("slope" ,parnames.in)]
  var.dais =parameters.in[match("var.dais" ,parnames.in)]
  Tcrit =parameters.in[match("Tcrit" ,parnames.in)]
  lambda =parameters.in[match("lambda" ,parnames.in)]

#  dais.out = dais_fastdynF(
  dais.out = daisanto_fastdynF(
                    anto.a=anto.a, anto.b=anto.b,
                    gamma=gamma  , alpha=alpha  ,
                    mu=mu        , nu=nu        ,
                    P0=P0        , kappa=kappa  ,
                    f0=f0        , h0=h0        ,
                    c=c          , b0=b0        ,
                    slope=slope  , SL=SL.in     ,
                    dSL=dSL.in   , includes_dSLais=1,
                    #Ta=Ta.in     , Toc=Toc.in   ,
                    Tg=Tg.in      ,
                    slope.Ta2Tg = slope.Ta2Tg,
                    intercept.Ta2Tg = intercept.Ta2Tg,
                    Tcrit=Tcrit   , lambda=lambda
                    )

  dais.out.norm = dais.out$Vais - mean(dais.out$Vais[ind.norm.in])

  llik0  = 0

  # Calculate the residuals

  # First check if the modeled AIS-only SLR exceeds total
  # Use first 20 years and only check residuals after that for exceeding total SLR
  # because there will be noise around the trend in the first 20 years; don't want
  # to throw out the run because of a little noise.

  resid.sl = ( SL.new-mean(SL.new[1:20]) ) - ( dais.out$Vais[midx.sl]-mean(dais.out$Vais[midx.sl[1:20]]) )

  # Only admit runs which fit the Pliocene window, to 1-sigma
  resid.lig = abs(obs.in[1]+0-dais.out.norm[obs.step.in[1]])

if(!is.na(resid.lig)){
  if(all(resid.sl[20:length(resid.sl)]>0) & (resid.lig < 2.0*obs.err.in[1])){

    # Calculate the modeled trend for the IPCC periods
    # 1993-2010
    x=seq(ind.trends[1,1],ind.trends[1,2]); y=dais.out.norm[ind.trends[1,1]:ind.trends[1,2]]
    barx=mean(x); bary=mean(y)
    trend.1993 = sum( (x-barx)*(y-bary) )/sum( (x-barx)^2 )
    # 1992-2001
    x=seq(ind.trends[2,1],ind.trends[2,2]); y=dais.out.norm[ind.trends[2,1]:ind.trends[2,2]]
    barx=mean(x); bary=mean(y)
    trend.1992 = sum( (x-barx)*(y-bary) )/sum( (x-barx)^2 )
    # 2002-2011
    x=seq(ind.trends[3,1],ind.trends[3,2]); y=dais.out.norm[ind.trends[3,1]:ind.trends[3,2]]
    barx=mean(x); bary=mean(y)
    trend.2002 = sum( (x-barx)*(y-bary) )/sum( (x-barx)^2 )

    resid = c( obs.in[1] - dais.out.norm[obs.step.in[1]] ,
               obs.in[2] - dais.out.norm[obs.step.in[2]] ,
               obs.in[3] - dais.out.norm[obs.step.in[3]] ,
               obs.in[4] - (dais.out.norm[obs.step.in[4]]-dais.out.norm[i1992]) ,
               trends.ais.in[1] - trend.1993 ,
               trends.ais.in[2] - trend.1992 ,
               trends.ais.in[3] - trend.2002
               )

    # Calculate the likelihood. The observations are not correlated. They are independent
    llik0 = sum (dnorm(resid, mean=rep(0,length(resid)), sd = c(sqrt(var.dais + obs.err.in[1:4]^2), trends.err.in), log=TRUE))
#    llik0 = sum (dnorm(resid, mean=rep(0,length(resid)), sd = c(sqrt(var.dais + obs.err.in[1:4]^2)), log=TRUE))

    llik = llik0 # assume residuals are independent
  } else {
    llik = -Inf
  }
} else {
llik=-Inf
}
  llik
}
##==============================================================================
## Log-prior
log.pri = function( parameters.in,
                    parnames.in,
                    alpha.var.in,
                    beta.var.in,
                    shape.lambda=NULL,
                    rate.lambda=NULL,
                    shape.Tcrit=NULL,
                    rate.Tcrit=NULL,
                    bound.lower.in,
                    bound.upper.in
                    )
{

  ind.lambda = match("lambda",parnames.in)
  ind.var = match("var.dais",parnames.in)
  var.dais =parameters.in[ind.var]
  lambda =parameters.in[match("lambda",parnames.in)]
  Tcrit =parameters.in[match("Tcrit",parnames.in)]

  # var.dais has inverse gamma prior, so there is a lower bound at 0 but no
  # upper bound
  if(fd.priors=='u') {
    in.range = all(parameters.in[1:(ind.var-1)] < bound.upper.in[1:(ind.var-1)]) &
               all(parameters.in > bound.lower.in)
  }
  if(fd.priors=='g') {
    in.range = all(parameters.in[1:(ind.lambda-1)] < bound.upper.in[1:(ind.lambda-1)]) &
               all(parameters.in[1:(ind.lambda-1)] > bound.lower.in[1:(ind.lambda-1)]) &
               all(parameters.in[ind.var] > bound.lower.in[ind.var]) &
               all(parameters.in[ind.lambda] > 0) &
               all(parameters.in[ind.lambda+1] < 0)
  }

  var_pri = 0
  lambda_pri = 0
  Tcrit_pri = 0

  if(in.range) {
    var_pri = (-alpha.var.in - 1)*log(var.dais) + (-beta.var.in/var.dais)
    if(fd.priors=='g') {
      lambda_pri = dgamma(x=lambda, shape=shape.lambda, rate=rate.lambda, log=TRUE)
      Tcrit_pri = dgamma(x=-Tcrit, shape=shape.Tcrit, rate=rate.Tcrit, log=TRUE)
    }
    lpri=0 + var_pri + lambda_pri + Tcrit_pri
  } else {
    lpri = -Inf
  }

  lpri
}
##==============================================================================
## Log-posterior distribution:  posterior ~ likelihood * prior
log.post = function(parameters.in,
                    parnames.in,
                    bound.lower.in,
                    bound.upper.in,
                    obs.in, obs.err.in, obs.step.in,
                    trends.ais.in , trends.err.in , ind.trends.in ,
                    ind.norm.in,
                    alpha.var,
                    beta.var,
                    shape.lambda=NULL,
                    rate.lambda=NULL,
                    shape.Tcrit=NULL,
                    rate.Tcrit=NULL,
                    slope.Ta2Tg.in=1,
                    intercept.Ta2Tg.in=0,
                    Tg.in=NULL,
                    Toc.in=NULL,
                    Ta.in=NULL,
                    SL.in,
                    dSL.in )
{
  lpri = log.pri(parameters.in , parnames.in=parnames.in,
                 alpha.var.in=alpha.var, beta.var.in=beta.var,
                 shape.lambda=shape.lambda, rate.lambda=rate.lambda,
                 shape.Tcrit=shape.Tcrit, rate.Tcrit=rate.Tcrit,
                 bound.lower.in=bound.lower.in, bound.upper.in=bound.upper.in )
  if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
    lpost = log.lik( parameters.in=parameters.in  , parnames.in=parnames.in ,
                     obs.in=obs.in                , obs.err.in=obs.err.in   ,
                     obs.step.in=obs.step.in      , ind.trends.in=ind.trends.in ,
                     trends.ais.in=trends.ais.in  , trends.err.in = trends.err.in ,
                     slope.Ta2Tg.in=slope.Ta2Tg.in, intercept.Ta2Tg.in=intercept.Ta2Tg.in, Tg.in=Tg.in,
                     #Ta.in=Ta.in                  , Toc.in=Toc.in     ,
                     SL.in=SL.in                  , dSL.in=dSL.in     ,
                     ind.norm.in=ind.norm.in )  +  lpri
  } else {
    lpost = -Inf
  }
  lpost
}
