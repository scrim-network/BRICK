#===============================================================================
# likelihood function(s), priors, posterior
# for GEV distribution
#
# Questions? Tony Wong (twong@psu.edu)
#===============================================================================


#
#===============================================================================
# log(prior) for gev model
#===============================================================================
#
log_prior_gev <- function(parameters,
                          parnames,
                          priors,
                          auxiliary=NULL
){
  lpri <- 0

  for (par in parnames) {
    parameter.value <- as.numeric(parameters[match(par,parnames)])
    if(priors[[par]]$type=='normal') {
      lpri <- lpri + dnorm(x=parameter.value, mean=priors[[par]]$mean, sd=priors[[par]]$sd, log=TRUE)
    } else if(priors[[par]]$type=='gamma') {
      lpri <- lpri + dgamma(x=parameter.value, shape=priors[[par]]$shape, rate=priors[[par]]$rate, log=TRUE)
    } else if(priors[[par]]$type=='uniform') {
      lpri <- lpri + dunif(x=parameter.value, min=priors[[par]]$lower, max=priors[[par]]$upper, log=TRUE)
    }
  }

  return(lpri)
}
#===============================================================================


#
#===============================================================================
# log(likelihood) for gev model
#===============================================================================
#
log_like_gev <- function(parameters,
                         parnames,
                         data_calib,
                         auxiliary=NULL
){
  llik <- 0
  #print(parameters)
  n.param <- length(parnames)
  if(n.param==3) {
    # fit a standard stationary GEV
    mu <- parameters[match('mu',parnames)]
    sigma <- parameters[match('sigma',parnames)]
    xi <- parameters[match('xi',parnames)]
  } else if(n.param==4) {
    # location parameter nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    sigma <- parameters[match('sigma',parnames)]
    xi <- parameters[match('xi',parnames)]
    mu <- mu0 + mu1*auxiliary
  } else if(n.param==5) {
    # location and scale parameters nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    xi <- parameters[match('xi',parnames)]
    mu <- mu0 + mu1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
  } else if(n.param==6) {
    # location, scale and shape all nonstationary
    mu0 <- parameters[match('mu0',parnames)]
    mu1 <- parameters[match('mu1',parnames)]
    sigma0 <- parameters[match('sigma0',parnames)]
    sigma1 <- parameters[match('sigma1',parnames)]
    xi0 <- parameters[match('xi0',parnames)]
    xi1 <- parameters[match('xi1',parnames)]
    mu <- mu0 + mu1*auxiliary
    sigma <- exp(sigma0 + sigma1*auxiliary)
    xi <- xi0 + xi1*auxiliary
  } else {print('ERROR - invalid number of parameters for GEV')}
  llik <- sum(devd(data_calib, loc=mu, scale=sigma, shape=xi, log=TRUE, type='GEV'))
  return(llik)
}
#===============================================================================


#
#===============================================================================
# log(post) for gev model
#===============================================================================
#
log_post_gev <- function(parameters,
                         parnames,
                         data_calib,
                         priors,
                         auxiliary
){
  lpost <- 0
  llik <- 0
  lpri <- 0

  # calculate prior
  lpri <- log_prior_gev(parameters=parameters,
                        parnames=parnames,
                        priors=priors,
                        auxiliary=auxiliary)

  if(is.finite(lpri)){
    # calculate likelihood (only if parameters pass the prior test)
    llik <- log_like_gev(parameters=parameters,
                         parnames=parnames,
                         data_calib=data_calib,
                         auxiliary=auxiliary)
  }

  lpost <- lpri + llik
  return(lpost)
}
#===============================================================================


#===============================================================================
# End
#===============================================================================
