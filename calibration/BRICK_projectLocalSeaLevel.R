##==============================================================================
##
## Function to project local sea level from a given set of global mean sea level
## contributions and a given lat/lon coordinate.
##
##  Required input: (set below)
##		lat.in										[input] latitude (deg N) you want projected local sea level
##		lon.in										[input] longitude (deg E) you want projected local sea level
##		filename.brickin					[input] file name for BRICK physical model results for NOLA
##		filename.brickout					[output] file name for BRICK physical model results for Vietnam
##
##  Output:
##		[filename.brickout]			netCDF4 file with the BRICK physical model output
##														(should be as many runs as post-calibrated parameters).
##		rc											return code; 0 if no problems, 1 if problems
##
## Questions? Tony Wong (twong@psu.edu)
##
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

BRICK_projectLocalSeaLevel <- function(
	lat.in,
	lon.in,
	filename.brickin,
	filename.brickout
	){

	rc=0

	##==============================================================================
	##==============================================================================
	## Read from post-calibrated physical model results netCDF output file

	print(paste("Reading post-calibrated model output from file ",filename.brickin,sep=""))

		ncdata <- nc_open(filename.brickin)

		gmsl.rcp26 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
		gmsl.rcp45 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
		gmsl.rcp85 = ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')

		ais.rcp26 = ncvar_get(ncdata, 'AIS_RCP26')
		ais.rcp45 = ncvar_get(ncdata, 'AIS_RCP45')
		ais.rcp85 = ncvar_get(ncdata, 'AIS_RCP85')

		gis.rcp26 = ncvar_get(ncdata, 'GIS_RCP26')
		gis.rcp45 = ncvar_get(ncdata, 'GIS_RCP45')
		gis.rcp85 = ncvar_get(ncdata, 'GIS_RCP85')

		gsic.rcp26 = ncvar_get(ncdata, 'GSIC_RCP26')
		gsic.rcp45 = ncvar_get(ncdata, 'GSIC_RCP45')
		gsic.rcp85 = ncvar_get(ncdata, 'GSIC_RCP85')

		te.rcp26 = ncvar_get(ncdata, 'TE_RCP26')
		te.rcp45 = ncvar_get(ncdata, 'TE_RCP45')
		te.rcp85 = ncvar_get(ncdata, 'TE_RCP85')

		temp.rcp26 = ncvar_get(ncdata, 'temp_RCP26')
		temp.rcp45 = ncvar_get(ncdata, 'temp_RCP45')
		temp.rcp85 = ncvar_get(ncdata, 'temp_RCP85')

		ocheat.rcp26 = ncvar_get(ncdata, 'ocheat_RCP26')
		ocheat.rcp45 = ncvar_get(ncdata, 'ocheat_RCP45')
		ocheat.rcp85 = ncvar_get(ncdata, 'ocheat_RCP85')

		t.proj =ncvar_get(ncdata, 'time_proj')
		ens =ncvar_get(ncdata, 'ens')

		nc_close(ncdata)

	## Other set up
	## NOTE: IPCC generally are relative to 1986-2005.
	## Therefore use that period to be commensurate with their results
	begyear.norm = 1986
	endyear.norm = 2005
	ind.norm = which(t.proj==begyear.norm):which(t.proj==endyear.norm)
	n.ensemble = length(ens)

	##==============================================================================
	##==============================================================================







	##==============================================================================
	##==============================================================================
	## Local sea level rise

	print(paste('Beginning fingerprinting to local sea level rise...',sep=''))
	print(paste('  lat=',lat.in,', lon=',lon.in,sep=''))

	source('../R/BRICK_LSL.R')

	lsl.rcp26 = brick_lsl(lat.in=lat.in,
																	lon.in=lon.in,
																	n.time=length(t.proj),
																	slr_gis = gis.rcp26,
																	slr_gsic = gsic.rcp26,
																	slr_ais = ais.rcp26,
																	slr_te = te.rcp26
																	)
	lsl.rcp45 = brick_lsl(lat.in=lat.in,
																	lon.in=lon.in,
																	n.time=length(t.proj),
																	slr_gis = gis.rcp45,
																	slr_gsic = gsic.rcp45,
																	slr_ais = ais.rcp45,
																	slr_te = te.rcp45
																	)
	lsl.rcp85 = brick_lsl(lat.in=lat.in,
																	lon.in=lon.in,
																	n.time=length(t.proj),
																	slr_gis = gis.rcp85,
																	slr_gsic = gsic.rcp85,
																	slr_ais = ais.rcp85,
																	slr_te = te.rcp85
																	)

	# And normalize sea-level rise
	for (i in 1:n.ensemble) {
		lsl.rcp26[,i] = lsl.rcp26[,i] - mean(lsl.rcp26[ind.norm,i])
		lsl.rcp45[,i] = lsl.rcp45[,i] - mean(lsl.rcp45[ind.norm,i])
		lsl.rcp85[,i] = lsl.rcp85[,i] - mean(lsl.rcp85[ind.norm,i])
	}

	if(any(is.na(lsl.rcp26))) {rc=1; print('LSL RCP26 has NA');}
	if(any(is.na(lsl.rcp45))) {rc=1; print('LSL RCP45 has NA');}
	if(any(is.na(lsl.rcp85))) {rc=1; print('LSL RCP85 has NA');}

	print(paste('... finished local sea level rise',sep=''))

	##==============================================================================
	##==============================================================================







	##==============================================================================
	##==============================================================================
	## Write a netCDF ensemble output file including each of the RCP scenarios:
	## (1) global total sea level, (2) local sea level. Also will want each
	## contribution to global sea level rise, for the hindcast plots

	library(ncdf4)

	dim.tproj <- ncdim_def('time_proj', 'years', as.double(t.proj))
	dim.lat <- ncdim_def('lat', 'deg N', as.double(length(lat.in)))
	dim.lon <- ncdim_def('lon', 'deg E', as.double(length(lon.in)))

	dim.ensemble <- ncdim_def('ens', 'ensemble member', as.double(1:ncol(lsl.rcp26)), unlim=TRUE)

	lat.out <- ncvar_def('lat.lsl', 'deg N', list(dim.lat), -999,
	                  longname = 'latitude of local sea level point')
	lon.out <- ncvar_def('lon.lsl', 'deg N', list(dim.lon), -999,
	                  longname = 'longitude of local sea level point')

	lsl26 <- ncvar_def('LocalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'Local sea level (RCP26)')
	gsic26 <- ncvar_def('GSIC_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'GSIC contribution to sea level (RCP26)')
	te26 <- ncvar_def('TE_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'TE contribution to sea level (RCP26)')
	gis26 <- ncvar_def('GIS_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'GIS contribution to sea level (RCP26)')
	ais26 <- ncvar_def('AIS_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'AIS contribution to sea level (RCP26)')
	gsl26 <- ncvar_def('GlobalSeaLevel_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'Global sea level (RCP26)')
	temp26 <- ncvar_def('temp_RCP26', 'deg C', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'global mean surface temperature anomaly (RCP26)')
	ocheat26 <- ncvar_def('ocheat_RCP26', '10^22 J', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'ocean heat uptake (RCP26)')

	lsl45 <- ncvar_def('LocalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'Local sea level (RCP45)')
	gsic45 <- ncvar_def('GSIC_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'GSIC contribution to sea level (RCP45)')
	te45 <- ncvar_def('TE_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'TE contribution to sea level (RCP45)')
	gis45 <- ncvar_def('GIS_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'GIS contribution to sea level (RCP45)')
	ais45 <- ncvar_def('AIS_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'AIS contribution to sea level (RCP45)')
	gsl45 <- ncvar_def('GlobalSeaLevel_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'Global sea level (RCP45)')
	temp45 <- ncvar_def('temp_RCP45', 'deg C', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'global mean surface temperature anomaly (RCP45)')
	ocheat45 <- ncvar_def('ocheat_RCP45', '10^22 J', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'ocean heat uptake (RCP45)')

	lsl85 <- ncvar_def('LocalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'Local sea level (RCP85)')
	gsic85 <- ncvar_def('GSIC_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'GSIC contribution to sea level (RCP85)')
	te85 <- ncvar_def('TE_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'TE contribution to sea level (RCP85)')
	gis85 <- ncvar_def('GIS_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'GIS contribution to sea level (RCP85)')
	ais85 <- ncvar_def('AIS_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'AIS contribution to sea level (RCP85)')
	gsl85 <- ncvar_def('GlobalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'Global sea level (RCP85)')
	temp85 <- ncvar_def('temp_RCP85', 'deg C', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'global mean surface temperature anomaly (RCP85)')
	ocheat85 <- ncvar_def('ocheat_RCP85', '10^22 J', list(dim.tproj, dim.ensemble), -999,
	                  longname = 'ocean heat uptake (RCP85)')

	outnc <- nc_create(filename.brickout,
											list( gsl26, gsl45, gsl85, lsl26, lsl45, lsl85,
														gsic26, te26, gis26, ais26, temp26, ocheat26,
														gsic45, te45, gis45, ais45, temp45, ocheat45,
														gsic85, te85, gis85, ais85, temp85, ocheat85,
														lat.out, lon.out),
											force_v4 = TRUE)

	ncvar_put(outnc, lat.out, lat.in)
	ncvar_put(outnc, lon.out, lon.in)

	ncvar_put(outnc, lsl26, lsl.rcp26)
	ncvar_put(outnc, gsl26, gmsl.rcp26)
	ncvar_put(outnc, temp26, temp.rcp26)
	ncvar_put(outnc, ocheat26, ocheat.rcp26)
	ncvar_put(outnc, gsic26, gsic.rcp26)
	ncvar_put(outnc, te26, te.rcp26)
	ncvar_put(outnc, gis26, gis.rcp26)
	ncvar_put(outnc, ais26, ais.rcp26)

	ncvar_put(outnc, lsl45, lsl.rcp45)
	ncvar_put(outnc, gsl45, gmsl.rcp45)
	ncvar_put(outnc, temp45, temp.rcp45)
	ncvar_put(outnc, ocheat45, ocheat.rcp45)
	ncvar_put(outnc, gsic45, gsic.rcp45)
	ncvar_put(outnc, te45, te.rcp45)
	ncvar_put(outnc, gis45, gis.rcp45)
	ncvar_put(outnc, ais45, ais.rcp45)

	ncvar_put(outnc, lsl85, lsl.rcp85)
	ncvar_put(outnc, gsl85, gmsl.rcp85)
	ncvar_put(outnc, temp85, temp.rcp85)
	ncvar_put(outnc, ocheat85, ocheat.rcp85)
	ncvar_put(outnc, gsic85, gsic.rcp85)
	ncvar_put(outnc, te85, te.rcp85)
	ncvar_put(outnc, gis85, gis.rcp85)
	ncvar_put(outnc, ais85, ais.rcp85)

	nc_close(outnc)

	print(paste('... finished writing output file ',filename.brickout,sep=''))

##==============================================================================
##==============================================================================

	return(rc)
}

##==============================================================================
## End
##==============================================================================
