# ====================================================================================
# Tuning of Van Dantzig model to NOLA
#
# by Alexander Bakker
# ====================================================================================

source("../R/VanDantzig.R")

# ====================================================================================
# Load R-libs
# ====================================================================================
require(lhs)

# ====================================================================================
# VD_NOLA_R:
#    VanDantzig tuned to NOLA (New Orleans, Louisiana)
# ====================================================================================
VD_NOLA_R <-
  function(params, I_table, T, X, global_sea_level){

    nx <- length(X)    # number of dike heightenings to be evaluated
    np <- nrow(params) # number of parameter sets
    nt <- T+1          # length of time horizon (including initial year t0)
    ts <- 0:T          # considered years

    outcome <- as.data.frame(
      matrix(c(X, rep(NA,6*nx)), nx, 7, byrow=F,
             dimnames=list(NULL,c("Heightening","Expected_costs","Expected_loss","Investment","Total_p_fail","Average_p_fail","Maximum_p_fail"))))

    if(np == 1) {

      # failure loss for all parameter samples and timesteps
      failure_loss <- as.matrix(failure_loss_VD(params$V0, params$delta, ts), nrow=nt,ncol=1)

      # global sea level change
      GSL <- global_sea_level
      # normalized to t0
      GSL = GSL-GSL[1]

      # relative local sea level rise (subsidence + GSL)
      LSL <- t(params$sub_rate * t(ts)) + GSL

      for(i in 1:nx) {
        investment <- (intercept.h2i + slope.h2i*X[i])*(1 + params$I_unc)

        p_fail     <- as.matrix(failure_probability_VD(params$p0, params$alpha, LSL, X[i] ), nrow=nt,ncol=1)

        outcome[i,2:7] <- VanDantzig_outcomes(p_fail, failure_loss, investment)
      }
    } else {

      # failure loss for all parameter samples and timesteps
      failure_loss <- t(sapply(1:length(ts), function(j) {
        failure_loss_VD(params$V0, params$delta, ts[j]) }))

      # global sea level rise
      #GSL <- t(sapply(1:length(ts), function(j) {
      # pola_R(a=params$init_gsl_anomaly, b=params$init_gsl_rate,
      #        c=params$gsl_acceleration, c_star=params$abrupt_rate_increase,
      #        t_star=params$timing_abrupt_increase,
      #        t=ts[j])
      #}))

      GSL <- global_sea_level

      # relative local sea level rise (subsidence + GSL)
      LSL <- t(params$sub_rate %*% t(ts)) + GSL

      for(i in 1:nx) {
#        investment <-
#          approx(x=I_table[,1], y=I_table[,2], X[i])$y #* (1 + params$I_unc)
          investment <- intercept.h2i + slope.h2i*X[i]

        p_fail     <- t(sapply(1:length(ts), function(j) {
          failure_probability_VD(
            params$p0, params$alpha, LSL[j,], X[i]) }))

        outcome[i,2:7] <- VanDantzig_outcomes(t(p_fail), failure_loss, investment)
      }
    }

    return(outcome)
  }
# ====================================================================================




# ====================================================================================
# parameter_sampling_DP:
#    Samples parameters based on the 'Dutch Perspectives' appendix E applying a
#    Latin Hypercube Sampling.
#
# (subroutine adjusted from R-code Perry Oddo (2015) )
# ====================================================================================
parameter_sampling_DP <- function(n_obs, p0, alpha, V0, delta, I_range, sub_rate) {
  z <- randomLHS(n_obs, 6)

  z[,1] <- qlnorm(z[,1],       log(p0),       0.25)
  z[,2] <-  qnorm(z[,2],         alpha,        0.1)
  z[,3] <-  qunif(z[,3],         V0[1],      V0[2])
  z[,4] <-  qunif(z[,4],      delta[1],   delta[2])
  z[,5] <-  qunif(z[,5],    I_range[1], I_range[2])
  z[,6] <- qlnorm(z[,6], log(sub_rate),        0.4) # sdlog to yield about stdev from Dixon et al. (2006)

  parameters         <- as.data.frame(z)
  names(parameters)  <- c("p0", "alpha", "V0", "delta", "I_unc", "sub_rate")

  return(parameters)
}
# ====================================================================================


# ====================================================================================
#
# ====================================================================================
