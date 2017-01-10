# =======================================================================================
# SIMPLE-fortran90 (# estimation by calling fortran routine)
# SIMPLE: Simple model for Greenland ice-sheet volume [m sle] (Bakker et al 2014)
# =======================================================================================
#
# #  Requires (input variables):
#  - Tg        global temperature [degC]
#
#  Simulates (output variables):
#  - Vgrl      Volume of Greenland ice sheet [m sea-level equivalent (sle)]
#
#  Internal variables:
#  - Veq       Equilibrium Volume [m sle]
#  - tau       Decay rate [1/ degC]
#
#  Parameters:
#  - a         sensitivity of equilibrium volume Veq [m sle/degC]
#  - b         equilibrium volume Veq [m sle] for temperature Tg = 0
#  - alpha     sensitivity of exponential decay rate (1/tau)
#  - beta      exponential decay rate [1 / degC] at Tg = 0
#  - Vgrl_0    initial ice-sheet volume [m sle]
#  - tstep     time step
# =======================================================================================

# =======================================================================================
# load fortran subroutine (# to check if library is loaded is.loaded("run_simple") )
dyn.load("../fortran/simple.so")

simpleF <- function(a      = -0.827,
                    b      =  7.242,
                    alpha  =  1.630e-4,
                    beta   =  2.845e-05,
                    V0     =  7.242,
                    tstep  = 1,
                    Tg)
{

  # determine series length
  ns <- length(Tg)

  # call fortran
  f.output <- .Fortran("run_simple",
                  ns            = ns,
                  tstep         = as.double(tstep),
                  simple_a      = as.double(a),
                  simple_b      = as.double(b),
                  simple_alpha  = as.double(alpha),
                  simple_beta   = as.double(beta),
                  simple_V0     = as.double(V0),
                  Grl_Temp      = as.double(Tg),
                  GIS_Volume_out = as.double(rep(-999.99,ns)),
                  lunstable_out = as.logical(FALSE)
  )
  if(f.output$lunstable_out) {f.output$GIS_Volume_out = rep(Inf,ns)}
  sle.gis <- V0 - f.output$GIS_Volume_out

  output <- list(Vgrl=f.output$GIS_Volume_out, sle.gis=sle.gis)

  return(output)

}
# =================================================================================
