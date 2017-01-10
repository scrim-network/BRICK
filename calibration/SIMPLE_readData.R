##==============================================================================
## Script to grab the Greenland ice sheet mass balance data to calibrate SIMPLE
## model.
##
## Read in the Greenland ice sheet observations in SLE, observation errors,
## years, and historic and emission temps.
## Read in the Greenland Mass balance data from 1958 - 2009
## The Mass balance is estimated as MB = surface mass balance + Discharge (Discharge is negative)
## Original data is from Sasgen et al (2012)
##
## Questions? Tony Wong <twong@psu.edu>
##==============================================================================

dat <- read.csv('../data/Greenland_OBS_MAR_InSAR_updated.csv')
obs.gis.time <- dat[1:52,9]
obs.gis <- dat[1:52,12] # Annual Mass Balance in meters sea level equivalence
obs.gis.err <- dat[1,15] # The error is +/- 30 Gt

idx = compute_indices(obs.time=obs.gis.time, mod.time=mod.time)
oidx.gis = idx$oidx; midx.gis = idx$midx

##==============================================================================
## End
##==============================================================================
