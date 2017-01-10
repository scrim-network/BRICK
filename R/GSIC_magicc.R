# =======================================================================================
# Simple model to simulate global contribution of Glaciers and Small Ice Caps (GSIC) to
# sea-level rise (Wigley and Raper 2005, application of equation 4/5)
# =======================================================================================
#
#  Requires (input variables):
#  - Tg        global temperature [degC]
#
#  Simulates (output variables):
#  - Gs        cumulative sea-level contribution since t0 (i=1)a [m]
#              by definition Gs(1) = 0
#
#  Parameters:
#  - beta0     initial mass balance sensitivity (how long it takes GSIC to respond to
#              increasing temps) [m/yr/C]
#  - V0        initial volume = max(Gs) [meter sle]
#  - n         exponent for area-volume scaling [-]
#  - Gs0       Gs[1]: the corrected corresponding sea-level rise in 1961 [m]
#  - Teq       equilibrium temperature (at which there is no change) [deg C]
#  - tstep     time step [years]
#
# =======================================================================================

gsic_magicc = function(beta0 = 0.000577,
                       V0    = 0.4,
                       n     = 0.82,
                       Gs0   = 0,
                       Teq   = -0.15,
                       tstep = 1,
                       Tg) {

  ns    <- length(Tg)

  Gs    <- rep(NA, ns)
  Gs[1] <- Gs0

  for(i in 2:ns){
    Gs[i] <- Gs[i-1] + tstep * (beta0 * (Tg[i-1] - Teq) * (1-(Gs[i-1]/V0))^n)
  }

  return(Gs)
}
