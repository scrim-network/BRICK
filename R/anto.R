# =======================================================================================
# anto:
# simple scaling to estimate Antarctic ocean temperature (Toc) from global temperature (Tg)
# =======================================================================================
#
#  Requires (input variables):
#  - Tg        Global surface temperature [anomaly degC]
#
#  Simulates (output variables):
#  - Toc       Temperature of ocean surface at Antarctica [degC]
#
#  Parameters:
#  - a        Sensitivity Ta [degC/degC]
#  - b        Ta(Tg=0)
#  - Tf       Freezing temperature sea-water (lower bound)
# =======================================================================================

anto <- function(a  =  0.26,
                  b  =  0.62,
                  Tf = -1.8,
                  Tg) {

  c  <- (Tf-b)/a
  Toc <-  Tf + (a*Tg + b-Tf) / (1 + exp(-Tg+c))

  return(Toc)
}
