# =======================================================================================
# Script to call GSIC_magicc (Fortran90):
# Simple model to simulate global contribution of Glaciers and Small Ice Caps (GSIC) to 
#  sea-level rise (Wigley and Raper 2005, application of equation 4/5)
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
# =======================================================================================

# =======================================================================================
# load fortran subroutine (# to check if library is loaded is.loaded("run_simple") )
dyn.load("../fortran/gsic_magicc.so")

gsic_magiccF = function(beta0 = 0.000577,
                       V0    = 0.4,
                       n     = 0.82,
                       Gs0   = 0,
                       Teq   = -0.15,
                       tstep = 1,
                       Tg) {

  # determine series length
  ns <- length(Tg)
  
  # call fortran
  f.output <- .Fortran("run_gsic_magicc",
                  ns                = ns,
                  tstep             = as.double(tstep),
                  gsic_magicc_beta0 = as.double(beta0),
                  gsic_magicc_V0    = as.double(V0),
                  gsic_magicc_n     = as.double(n),
                  gsic_magicc_Gs0   = as.double(Gs0),
                  gsic_magicc_Teq   = as.double(Teq),
                  Gl_Temp           = as.double(Tg),
                  SL_contribution_out = as.double(rep(-999.99,ns))
  )
  return(f.output$SL_contribution_out)
  
}
# =================================================================================