# =======================================================================================
# BRICK_TE:
# Simple model to simulate contribution of thermosteric expansion (TE) to global sea-
# level rise. (Grinsted, Moore and Jevrejeva 2010, equations 1/2 that were originally
# applied to global sea-level)
# =======================================================================================
#
#  Requires (input variables):
#  - Tg        global temperature [degC]
#
#  Simulates (output variables):
#  - TE        cumulative sea-level contribution from thermosteric expansion [m]
#              since t0 (i=1)
#              by definition Gs(1) = 0
#
#  Internal variables:
#  - TEeq      Equilibrium TE [m]
#
#  Parameters:
#  - a         sensitivity of TE equilibrium [m/K]
#  - b         equilibrium TE [m] for temperature anomaly Tg = 0
#  - invtau    1/time-scale of exponential decay (e-folding time) [years^-1]
#  - TE_0      initial sea level
#  - tstep     time step
# =======================================================================================

brick_te = function(a     = 0.5,
                    b     = 0.0,
                    invtau= 0.005,
                    TE_0  = 0.0,
                    tstep = 1,
                    Tg) {

  np    <- length(Tg)

  TE    <- rep(NA, np)
  TE[1] <- TE_0

  for(i in 2:np){
#    TEeq  <- a*Tg[i-1] + b
    TE[i] <- TE[i-1]  +  tstep * ( ( a*Tg[i-1] + b - TE[i-1]) * invtau) # sample 1/tau

  }

  return(TE)
}
