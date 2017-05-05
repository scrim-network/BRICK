##==============================================================================
## script to add the inferred land water storage contributions onto the "supported"
## projections/hindcasts from BRICK v0.2 accompanying the revised ms submitted
## to GMD.
## Why? Well, because Tony forgot to put them on the netCDF results file. but
## SLR = SLR_GSIC + SLR_GIS + SLR_AIS + SLR_TE + SLR_LWS, so we can infer the
## missing SLR_LWS and put it on a revised netCDF file.
## And all are relative to 1986-2005 (the projections)
##
## Questions? Tony Wong (twong@psu.edu)
##==============================================================================

library(ncdf4)

# initialize
slr <- vector('list',3); names(slr) <- c('rcp26','rcp45','rcp85')
gsic <- slr; gis <- slr; ais <- slr; te <- slr; lws <- slr

# file name from the original simulations
filename.projections <- '../output_model/BRICK-model_physical_control_02Apr2017.nc'

# read the data
ncdata <- nc_open(filename.projections)

slr$rcp85 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
gsic$rcp85 <- ncvar_get(ncdata, 'GSIC_RCP85')
gis$rcp85 <- ncvar_get(ncdata, 'GIS_RCP85')
ais$rcp85 <- ncvar_get(ncdata, 'AIS_RCP85')
te$rcp85 <- ncvar_get(ncdata, 'TE_RCP85')

slr$rcp45 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
gsic$rcp45 <- ncvar_get(ncdata, 'GSIC_RCP45')
gis$rcp45 <- ncvar_get(ncdata, 'GIS_RCP45')
ais$rcp45 <- ncvar_get(ncdata, 'AIS_RCP45')
te$rcp45 <- ncvar_get(ncdata, 'TE_RCP45')

slr$rcp26 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
gsic$rcp26 <- ncvar_get(ncdata, 'GSIC_RCP26')
gis$rcp26 <- ncvar_get(ncdata, 'GIS_RCP26')
ais$rcp26 <- ncvar_get(ncdata, 'AIS_RCP26')
te$rcp26 <- ncvar_get(ncdata, 'TE_RCP26')

t.proj <- ncvar_get(ncdata, 'time_proj')
n.ens <- length(ncvar_get(ncdata, 'ens'))
nc_close(ncdata)

# set inferred land water storage contributions
lws$rcp26 <- slr$rcp26 - gsic$rcp26 - gis$rcp26 - ais$rcp26 - te$rcp26
lws$rcp45 <- slr$rcp45 - gsic$rcp45 - gis$rcp45 - ais$rcp45 - te$rcp45
lws$rcp85 <- slr$rcp85 - gsic$rcp85 - gis$rcp85 - ais$rcp85 - te$rcp85

# add these to output file

#ncdata <- nc_open(filename.projections, write=TRUE)
ncdata <- nc_open('../output_model/output.nc', write=TRUE)

dim.tproj <- ncdim_def('time_proj', 'years', as.double(t.proj))
dim.ensemble <- ncdim_def('ens', 'ensemble member', as.double(1:n.ens), unlim=TRUE)

lws.rcp26 <- ncvar_def('LWS_RCP26', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'estimated LWS contribution to sea level (RCP26)')
lws.rcp45 <- ncvar_def('LWS_RCP45', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'estimated LWS contribution to sea level (RCP45)')
lws.rcp85 <- ncvar_def('LWS_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999,
                  longname = 'estimated LWS contribution to sea level (RCP85)')

ncdata <- ncvar_add( ncdata, lws.rcp26, verbose=FALSE, indefine=FALSE )
ncdata <- ncvar_add( ncdata, lws.rcp45, verbose=FALSE, indefine=FALSE )
ncdata <- ncvar_add( ncdata, lws.rcp85, verbose=FALSE, indefine=FALSE )

ncvar_put(ncdata, lws.rcp26, lws$rcp26)
ncvar_put(ncdata, lws.rcp45, lws$rcp45)
ncvar_put(ncdata, lws.rcp85, lws$rcp85)

nc_close(ncdata)


# check sea-level budget on modified history file
ncdata <- nc_open('../output_model/output.nc')
slr$rcp85 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP85')
gsic$rcp85 <- ncvar_get(ncdata, 'GSIC_RCP85')
gis$rcp85 <- ncvar_get(ncdata, 'GIS_RCP85')
ais$rcp85 <- ncvar_get(ncdata, 'AIS_RCP85')
te$rcp85 <- ncvar_get(ncdata, 'TE_RCP85')
slr$rcp45 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP45')
gsic$rcp45 <- ncvar_get(ncdata, 'GSIC_RCP45')
gis$rcp45 <- ncvar_get(ncdata, 'GIS_RCP45')
ais$rcp45 <- ncvar_get(ncdata, 'AIS_RCP45')
te$rcp45 <- ncvar_get(ncdata, 'TE_RCP45')
slr$rcp26 <- ncvar_get(ncdata, 'GlobalSeaLevel_RCP26')
gsic$rcp26 <- ncvar_get(ncdata, 'GSIC_RCP26')
gis$rcp26 <- ncvar_get(ncdata, 'GIS_RCP26')
ais$rcp26 <- ncvar_get(ncdata, 'AIS_RCP26')
te$rcp26 <- ncvar_get(ncdata, 'TE_RCP26')
lws$rcp26 <- ncvar_get(ncdata, 'LWS_RCP26')
lws$rcp45 <- ncvar_get(ncdata, 'LWS_RCP45')
lws$rcp85 <- ncvar_get(ncdata, 'LWS_RCP85')
t.proj <- ncvar_get(ncdata, 'time_proj')
n.ens <- length(ncvar_get(ncdata, 'ens'))
nc_close(ncdata)

imbal26 <- slr$rcp26 - gsic$rcp26 - gis$rcp26 - ais$rcp26 - te$rcp26 - lws$rcp26
max(abs(imbal26))
imbal45 <- slr$rcp45 - gsic$rcp45 - gis$rcp45 - ais$rcp45 - te$rcp45 - lws$rcp45
max(abs(imbal45))
imbal85 <- slr$rcp85 - gsic$rcp85 - gis$rcp85 - ais$rcp85 - te$rcp85 - lws$rcp85
max(abs(imbal85))

# returns about 10^-9, which is certainly less than what we care about for
# sea-level rise; netcdf double should be 15 or 16 digits, but lose some on
# subtraction

##==============================================================================
## End
##==============================================================================
