# =======================================================================================
# WAIS:
# Representation of deeply uncertain WAIS collapse as a function of 'time'
# (due to combined effect cliff failure and hydrofracturing)
# =======================================================================================
#
#  Requires (input variables):
#  - t        time [years from anomaly year]
#
#  Simulates (output variables):
#  - Vol      Ice volume left [m sle]
#
#  Parameters:
#  - dSmx       max disintegration rate [m sle/yr]
#  - start.year Starting year disintegration (defined as year where st.yr.pc % is disintegrated)
#  - st.yr.pc   Percentage of disintegration where disintegration is defined to start
#  - V0         Original volume [m sle]
# =======================================================================================

WAIS <- function(dSmx       = 3.3/40,
                 start.year = 0,
                 st.yr.pc   = 0.01,
                 V0         = 3.3,
                 t) {

  a <- dSmx/V0 / 0.25
  b <- start.year + log((1-st.yr.pc)/st.yr.pc) / a

  Vol <- V0 * (1 - 1/(1+exp(-a*(t - b))) )

  return(Vol)
}
