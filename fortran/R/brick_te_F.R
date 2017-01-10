# =======================================================================================
# BRICK_TE fortran90 (# estimation by calling fortran routine)
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

# =======================================================================================
# load fortran subroutine (# to check if library is loaded is.loaded("run_brick_te") )
dyn.load("../fortran/brick_te.so")

brick_te_F <- function(a     = 0.5,
                       b     = 0.0,
                       invtau= 0.005,
                       TE_0    = 0.0,
                       tstep = 1,
                       Tg) {

  # determine series length
  ns <- length(Tg)

  tau = 1/invtau

  # call fortran
  f.output <- .Fortran("run_brick_te",
                  ns            = ns,
                  tstep         = as.double(tstep),
                  brick_te_a    = as.double(a),
                  brick_te_b    = as.double(b),
                  brick_te_tau  = as.double(tau),
                  brick_te_TE_0 = as.double(TE_0),
                  Gl_Temp       = as.double(Tg),
                  TE_out        = as.double(rep(-999.99,ns))
  )
  return(f.output$TE_out)

}
# =================================================================================
