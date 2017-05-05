# =======================================================================================
# Local sea level mapping.
# Will return a time series of local sea level for the given lat/lon point, with
# same length as the slr_xxx input sea-level rise time series.
#
# If slr_xxx is fed in as a matrix of multiple time series, then return local
# sea level for each of them.
#
# Currently no support to output global map (must be done post-processing).
#
# Currently no support to output global map and an LSL time series for a specific
# location (this routine must be called twice, once for each).
# =======================================================================================
#
#   Requires (input variables):
# - lat.in      latitude of point at which you want local sea level (currently only support 1 pt)
# - lon.in      longitude of point at which you want local sea level
# - global      (logical) whether you want a global map or not
# - n.time      number of points in each time series
# - slr_gis     Greenland ice sheet component of sea-level rise (m sle)
# - slr_gsic    glaciers and small ice caps component of sea-level rise (m sle)
# - slr_ais     Antarctic ice sheet component of sea-level rise (m sle)
# - slr_te      thermal expansion component of sea-level rise (m sle)
# - slr_lws     land water storage component of sea-level rise (m sle)
#
#   Simulates (output variables):
# - lsl.out     local sea level rise (m sle)
#
#   Parameters:
# - filename.fingerprints = "../fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"
#               sea level fingerprinting data (Slangen et al 2014)
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

brick_lsl <- function(  lat.in = 0,
                        lon.in = 0,
                        global = FALSE,
                        n.time = 0,
                        slr_gis = 0,
                        slr_gsic = 0,
                        slr_ais = 0,
                        slr_te = 0,
                        slr_lws = 0
                        )
{

  # which of the two dimensions is time? if only one time series, turn into
  # a matrix for generality
  dims = dim(slr_gis)
  if(is.null(dims)) {
    slr_gis = as.matrix(slr_gis)
    slr_ais = as.matrix(slr_ais)
    slr_gsic = as.matrix(slr_gsic)
    slr_te = as.matrix(slr_te)
    slr_lws = as.matrix(slr_lws)
    itime = which(dim(slr_gis)==n.time)
  } else {
    itime = which(dims==n.time)
  }
  n.ensemble = dims[which(dims!=n.time)]
  if(is.null(n.ensemble)) {n.ensemble=1}

  # should choke if gis, ais, te, or gsic is missing, but LWS shoudl be optional
  # if missing, resize to fit other fields
  if(is.null(dim(slr_lws))) {
    slr_lws = mat.or.vec(nr=nrow(slr_gis), nc=ncol(slr_gis))
  }

  # set up output with number of time points as number of rows, then re-shape the
  # output to match what was entered as input (slr_gis)
  lsl.out <- mat.or.vec(n.time,n.ensemble)

  # reshape the input to match this convention
  if(itime!=1) {
    slr_gis_reshape = t(slr_gis)
    slr_ais_reshape = t(slr_ais)
    slr_gsic_reshape = t(slr_gsic)
    slr_te_reshape = t(slr_te)
    slr_lws_reshape = t(slr_lws)
  } else {
    slr_gis_reshape = slr_gis
    slr_ais_reshape = slr_ais
    slr_gsic_reshape = slr_gsic
    slr_te_reshape = slr_te
    slr_lws_reshape = slr_lws
  }

  # read the sea-level fingerprints in only one
  filename.fingerprints = "../fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"

  ncdata <- nc_open(filename.fingerprints)
    lat = ncvar_get(ncdata, 'lat')
    lon = ncvar_get(ncdata, 'lon')
    fp.gsic = ncvar_get(ncdata, 'GLAC')
    fp.gis = ncvar_get(ncdata, 'GIS')
    fp.ais = ncvar_get(ncdata, 'AIS')
  nc_close(ncdata)

  # convert longitude to degrees East, and find the fingerprinting data location
  # closest to the local sea level lat/lon given
  if(lon.in < 0) {lon.in=lon.in+360}	# convert longitude to degrees East
  ilat = which( abs(lat-lat.in)==min(abs(lat-lat.in)) )
  ilon = which( abs(lon-lon.in)==min(abs(lon-lon.in)) )

  # it is possible there were multiple lats/lons 'closest' to your given point
  # take the average of the non-NA of these
  fp.ais.loc = mean(fp.ais[ilon,ilat],na.rm=TRUE)
  fp.gsic.loc = mean(fp.gsic[ilon,ilat],na.rm=TRUE)
  fp.gis.loc = mean(fp.gis[ilon,ilat],na.rm=TRUE)
  fp.te.loc = 1.0       # TE response is to global mean temperature, so global mean sea level response is same everywhere
  fp.lws.loc = 1.0      # assume LWS changes uniformly (likely not quite true,
                        # but a small contribution anyway)

  # check if the nearest spot ended up on land.
  # If it did, take the average everywhere around the location.
  if(is.na(fp.ais.loc) | is.na(fp.gsic.loc) | is.na(fp.gis.loc) | is.na(fp.te.loc)) {
    fp.ais.loc = mean(fp.ais[(ilon-1):(ilon+1),(ilat-1):(ilat+1)], na.rm=TRUE)
    fp.gsic.loc = mean(fp.gsic[(ilon-1):(ilon+1),(ilat-1):(ilat+1)], na.rm=TRUE)
    fp.gis.loc = mean(fp.gis[(ilon-1):(ilon+1),(ilat-1):(ilat+1)], na.rm=TRUE)
    fp.te.loc = 1.0
    fp.lws.loc = 1.0
  }

  # error message if something is still wrong
  if(is.na(fp.ais.loc) | is.na(fp.gsic.loc) | is.na(fp.gis.loc) | is.na(fp.te.loc)) {
    print('WARNING -- local sea level fingerprints are NaN')
  }

  # for each ensemble member, fingerprint the time series
  for (i in 1:n.ensemble) {
    lsl.out[,i] = fp.gis.loc * slr_gis_reshape[,i] +
                  fp.ais.loc * slr_ais_reshape[,i] +
                  fp.gsic.loc * slr_gsic_reshape[,i] +
                  fp.te.loc * slr_te_reshape[,i] +
                  fp.lws.loc * slr_lws_reshape[,i]
  }

  # reshape output to match the input
  if (itime!=1) {
    lsl.out = t(lsl.out)
  }

  return(lsl.out)
}
