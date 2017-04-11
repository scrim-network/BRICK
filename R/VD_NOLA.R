# ====================================================================================
# Tuning of Van Dantzig model to New Orleans, Louisiana (NOLA)
#
# Original code by Alexander Bakker (July 2016)
# Modified code by Tony Wong (November 2016-present)
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

source("../R/VanDantzig.R") # this script contains the relevant Van Dantzig analysis functions

# ====================================================================================
# Load R-libs
# ====================================================================================
require(lhs)

# ====================================================================================
# VD_NOLA_R:
#    VanDantzig tuned to NOLA (New Orleans, Louisiana)
# ====================================================================================
VD_NOLA_R <- function(
  params,
  T,
  X,
  local_sea_level,
  intercept.h2i,
  slope.h2i,
  ss.gev=NULL,
  surge.rise=NULL,
  H0=4.25
  ){


    nx <- length(X)    # number of dike heightenings to be evaluated
    np <- nrow(params) # number of parameter sets
    nt <- T+1          # length of time horizon (including initial year t0)
    ts <- 0:T          # considered years

    outcome <- as.data.frame(
      matrix(c(X, rep(NA,6*nx)), nx, 7, byrow=F,
             dimnames=list(NULL,c("Heightening","Expected_costs","Expected_loss","Investment","Total_p_fail","Average_p_fail","Maximum_p_fail"))))

    # failure loss for all timesteps, accounting for V0 (the ucnertain value of
    # goods protected by the dike ring) and delta (the uncertain net discount
    # rate)
    failure_loss <- as.matrix(failure_loss_VD(params$V0, params$delta, ts), nrow=nt,ncol=1)

    # normalized local sea level change (taking 0 sea level as local_sea_level[1])
    LSL.norm <- local_sea_level - local_sea_level[1]

    # relative local sea level rise (subsidence + LSL.norm + surge rise (if non-stat))
    if(is.null(surge.rise)) {
      LSL.norm.subs <- t(params$sub_rate * t(ts)) + LSL.norm
    } else {
      LSL.norm.subs <- t(params$sub_rate * t(ts)) + LSL.norm + surge.rise
    }

    # correction from previous version: slope and intercept are fitted to
    # Jonkman et al (2009) Table 2, which is "surge level" (total height) vs
    # investment needed to build and maintain that structure
    investment <- (intercept.h2i + slope.h2i*(H0+X))*(1 + params$I_unc)

    # effective dike height is initial height + heightening - local sea level (accounting for subsidence and sea level rise)
    H_eff <- H0 + matrix(data=rep(X,nt),ncol=nt) - matrix(rep(LSL.norm.subs,nx),nrow=nx,byrow=TRUE)

    # Failure probability is the survival function of the storm surge GEV at effective dike height
    # NOTE: It is about 10x faster to calculate by hand than to use the function pgev
    # Note that this assumes shape parameter does not equal 0, and will return
    # NAN if you are outside the support (as determined by the three parameters)
    ###p_fail <- 1-sapply(1:nt, function(j) {pgev(q=1000*H_eff[,j], xi=ss.gev['shape'], mu=ss.gev['location'], beta=ss.gev['scale']) })
    ttmp <- (1+ss.gev['shape']*((1000*H_eff-ss.gev['location'])/ss.gev['scale']))^(-1/ss.gev['shape'])
    p_fail <- 1-exp(-ttmp)

    for(i in 1:nx) {
      outcome[i,2:7] <- VanDantzig_outcomes(matrix(p_fail[i,]), failure_loss, investment[i])
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
  z[,6] <- qlnorm(z[,6], log(sub_rate),        0.4) # sdlog to yield about stdev from Dixon et al. (2006), 2.5mm/y

  parameters         <- as.data.frame(z)
  names(parameters)  <- c("p0", "alpha", "V0", "delta", "I_unc", "sub_rate")

  return(parameters)
}
# ====================================================================================


# ====================================================================================
# End
# ====================================================================================
