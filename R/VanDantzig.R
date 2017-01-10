# ====================================================================================
# Extended Van Dantzig model for the assessment of flood risk and optimizing
# strategies
#
# file contains several subroutines to apply the analysis
#
# by Alexander Bakker
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

require(lhs)

# ====================================================================================
# VanDantzig_R
#    applies classic VanDantzig
# ====================================================================================
VanDantzig_R <-
  function(params, T, X, local_sea_level){

    nx <- length(X)    # number of dike heightenings to be evaluated
    np <- nrow(params) # number of parameter sets
    nt <- T+1          # length of time horizon (including initial year t0)
    ts <- 0:T          # considered years

    outcome <- as.data.frame(
      matrix(c(X, rep(NA,6*nx)), nx, 7, byrow=F,
             dimnames=list(NULL,c("Heightening","Expected_costs","Expected_loss","Investment","Total_p_fail","Average_p_fail","Maximum_p_fail"))))

    if(np == 1) {

      # failure loss for all parameter samples and timesteps
      failure_loss <- as.matrix(failure_loss_VD(params$V0, params$delta, ts),
                                nrow=nt,ncol=1)

      # LSL <- pola_R(a=params$init_gsl_anomaly[1], b=params$init_gsl_rate[1],
      #               c=params$gsl_acceleration[1], c_star=params$abrupt_rate_increase[1],
      #               t_star=params$timing_abrupt_increase[1],
      #               t=ts)

      LSL <- local_sea_level

      for(i in 1:nx) {
        investment <- params$k * X[i] + params$I0
        p_fail     <- as.matrix(failure_probability_VD(
                                params$p0, params$alpha, LSL, X[i] ),
                                nrow=nt,ncol=1)

        outcome[i,2:7] <- VanDantzig_outcomes(p_fail, failure_loss, investment)
      }
    } else {

      # failure loss for all parameter samples and timesteps
      failure_loss <- t(sapply(1:length(ts), function(j) {
        failure_loss_VD(params$V0, params$delta, ts[j]) }))

      # LSL <- t(sapply(1:length(ts), function(j) {
      #   pola_R(a=params$init_gsl_anomaly, b=params$init_gsl_rate,
      #          c=params$gsl_acceleration, c_star=params$abrupt_rate_increase,
      #          t_star=params$timing_abrupt_increase,
      #          t=ts[j])
      # }))
      LSL <- local_sea_level

      for(i in 1:nx) {
        investment <- params$k * X[i] + params$I0
        p_fail     <- t(sapply(1:length(ts), function(j) {
          failure_probability_VD(
            params$p0, params$alpha, LSL[j,], X[i]) }))

        outcome[i,2:7] <- VanDantzig_outcomes(p_fail, failure_loss, investment)
      }
    }

    return(outcome)
  }
# ====================================================================================


# ====================================================================================
# VanDantzig_outcomes:
#    returns the classic Van Dantzig outcomes, given a transient projection of
#    failure_probability (p_fail), losses in case of failure (failure_loss) and
#    investment costs (investment).
# ====================================================================================
VanDantzig_outcomes <- function(p_fail, failure_loss, investment) {

  # expected loss and costs
  expected_investment <- mean(investment)
  expected_loss       <- sum(apply(p_fail * failure_loss,1,mean))
  expected_costs      <- expected_loss + expected_investment

  # total, average and maximum failure probability during planning horizon
  p_fail_agg     <- apply(p_fail,1,mean)   # mean failure probability of all parameter sets
  total_p_fail   <- 1 - (prod(1 - p_fail_agg))
  ave_p_fail     <- mean(p_fail_agg)
  max_p_fail     <- max(p_fail_agg)

  return(c(expected_costs,
           expected_loss,
           expected_investment,
           total_p_fail,
           ave_p_fail,
           max_p_fail))
}
# ====================================================================================


# ====================================================================================
# failure_loss_VD:
#    Estimates economic loss in case of failure as a function of time in a classic
#    Van Dantzig fashion.
#    - discounting equation -
# ====================================================================================
failure_loss_VD <- function(V0   = 2e10, # Value at t0
                            delta = 0.02, # combines discount rate and value increase
                            ts) {
  loss <- V0 * (1 + delta)^(-ts)
  return(loss)
}
# ====================================================================================



# ====================================================================================
# failure_probability_VD:
#    Estimates the probability of dike failure (due to overtopping) as a function of
#    dike heightening (X) and LSL
#
# p0      failure probability at H0=4.25 m
# alpha   overtopping exponential constant
# LSL     local sea level anomaly (m) compared to t0
# H0      initial dike height (m)
# x       dike heightening (m)
# ss.gev  storm surge GEV parameters, time series; using mm, so use H_effective*1000
# ====================================================================================
failure_probability_VD <- function(
  p0 = 0.0038,
  alpha = 2.6,
  LSL,
  H0 = 4.25,
  x,
  ss.gev = NULL
  ) {

  p_fail = rep(NA,length(LSL))

# Note -- The time (t) loop is SLOW! Currently done in-line in VD_NOLA.R

  if(is.null(ss.gev)) {
    p_fail <- p0 * exp(- alpha * (x - LSL))
  } else {
    #p_fail <- p0 * exp(- alpha * (x - LSL))
    H_effective = H0+x-LSL
    p_fail <- 1-sapply(1:nrow(ss.gev), function(t) {pgev(q=1000*H_effective, xi=ss.gev[t,'shape'], mu=ss.gev[t,'location'], beta=ss.gev[t,'scale']) })[,1]
#    p_fail <- 1-sapply(1:nrow(ss.gev), function(t) {pgev(q=1000*H_effective, xi=ss.gev[t,'shape'], mu=ss.gev[t,'location'], beta=ss.gev[t,'scale']) })[,1]
#    for (t in 1:length(LSL)){
#      p_fail[t] = 1-pgev(q=1000*H_effective[t], xi=ss.gev[t,'shape'], mu=ss.gev[t,'location'], beta=ss.gev[t,'scale'])
#    }
  }

  return(p_fail)
}
# ====================================================================================


# ====================================================================================
# sample_parameters_VD:
#    Applies a Latin Hypercube Sampling and transformed according to the assumed
#    distributions of the parameters.
#
# (subroutine adjusted from R-code Perry Oddo (2015) )
# ====================================================================================
sample_parameters_VD <- function(n_obs, p0, alpha, V0, delta, k, I0) {

  z <- randomLHS(n_obs, 6)

  z[,1] <- qlnorm(z[,1], log(p0),       0.25)
  z[,2] <-  qnorm(z[,2],     alpha,     0.1)
  z[,3] <-  qnorm(z[,3],     V0,        1e3)
  z[,4] <- qlnorm(z[,4], log(delta),    0.1)
  z[,5] <-  qnorm(z[,5],     k,         4.0)
  z[,6] <- 0.0

  parameters         <- as.data.frame(z)
  names(parameters)  <- c("p0", "alpha", "V0", "delta", "k", "I0")

  return(parameters)
}
# ====================================================================================
