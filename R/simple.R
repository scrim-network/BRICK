# =======================================================================================
#  Simple, mechanistically motivated model of the Greenland ice sheet volume [m sle] in
#  response to temperature variations (Bakker, Applegate and Keller 2016, equations 1-3)
# =======================================================================================
#
#  Requires (input variables):
#  - Tg        global temperature [degC]
#
#  Simulates (output variables):
#  - Vgrl      Volume of Greenland ice sheet [m sea-level equivalent (sle)]
#  - sle.gis   sea level rise (m) from Greenland ice sheet
#
#  Parameters:
#  - a         sensitivity of equilibrium volume Veq [m sle/degC]
#  - b         equilibrium volume Veq [m sle] for temperature Tg = 0
#  - alpha     sensitivity of exponential decay rate (1/tau)
#  - beta      exponential decay rate [1 / K] at Tg = 0
#  - V0        initial ice-sheet volume [m sle]
#  - tstep     time step
#
# =======================================================================================

simple <- function(a      = -0.827,
                   b      =  7.242,
                   alpha  =  1.630e-4,
                   beta   =  2.845e-05,
                   V0     =  7.242,
                   tstep  = 1,
                   Tg)
{

  np      <- length(Tg)

  Vgrl    <- rep(NA, np)
  sle.gis <- rep(NA, np)
  Vgrl[1] <- V0

  for (i in 2:np) {
    Veq     <- a * Tg[i-1] + b              # 'virtual' equilibrium volume (equation 2)
    tauinv  <- alpha * Tg[i-1] + beta       # 1/timescale (tauinv = 1/tau) (equation 3)

    Vgrl[i] <- max(0,                                                   #  (equations 1 and 4)
                   Vgrl[i-1]  +  tstep * ( ( Veq - Vgrl[i-1]) * tauinv)  )

    # check tau to make sure it does not violate stability requirements
    #if(tauinv > (2/tstep) | tauinv<0) {Vgrl=rep(Inf,np); break;}

  }

  sle.gis <- V0 - Vgrl

  return(list(Vgrl=Vgrl,sle.gis=sle.gis))
}
