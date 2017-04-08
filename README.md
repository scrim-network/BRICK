![alt text](https://github.com/scrim-network/BRICK/blob/master/brick_logo.png "This is a brick!")
# BRICK v0.2

## Synopsis

BRICK (**B**uilding blocks for **R**elevant **I**ce and **C**limate **K**nowledge) is a modular semi-empirical modeling framework to simulate global temperature and sea-level rise. In the default model configuration, first, global mean surface temperature and ocean heat uptake are simulated by DOECLIM. Changes in global mean surface temperature drive changes in global mean sea level (GMSL). The contributions to GMSL from the Greenland and Antarctic ice sheets, thermal expansion, and glaciers and ice caps are simulated.

This repository contains the BRICK model source code and analysis scripts. Analysis codes and the main user interface are written in R and the main physics codes are in Fortran 90.

It contains the sub-models
* DOECLIM (Diffusion-Ocean-Energy CLIMate model)
* SNEASY (Simple Nonlinear EArth SYstem model)
* DAIS (Antarctic ice sheet)
* SIMPLE (Greenland ice sheet)
* GSIC-MAGICC (Glaciers and small ice caps)
* TE (thermal expansion)
* GMSL (Rahmstorf 2007, global mean sea level emulator)
* Van Dantzig flood risk assessment (Van Dantzig, 1956)

## Directory structure

./
   * BRICK "home" directory

./calibration/
   * R scripts related to the (pre-/post-)calibration of the physical models, including reading data, likelihood functions, and sub-model stand-alone calibration drivers (some of these may be antiquated)

./data/
   * data for radiative forcing (DOECLIM) and calibration

./fortran/
   * Fortran versions of all physical models are available, and are wrapped in R calling functions in ./fortran/R/

./output_model/
   * physical model output (i.e., temperature, sea-level rise)

./output_calibration/
   * statistical model output (i.e., posterior parameter values)

./R/
   * physical models in R

## Motivation

The motivation for the BRICK model is detailed in the [GMDD model description paper](http://www.geosci-model-dev-discuss.net/gmd-2016-303/):

> Here we introduce BRICK v0.1 ("Building blocks for Relevant Ice and Climate Knowledge"), a new model framework that focuses on accessibility, transparency, and flexibility while maintaining, as much as possible, the computational efficiency that make simple models so appealing. There is a wide range of potential applications for such a model. A simple framework enables uncertainty quantification via statistical calibration approaches (Higdon et al., 2004; Kennedy and O’Hagan, 2001), which would be infeasible with more computationally expensive models. A transparent modeling framework enables communication between scientists as well as communication with stakeholders. This leads to potential application of the model framework in decision support and education (Weaver et al., 2013). The present work expands on previous studies by (1) providing a platform of simple, but mechanistically motivated sea-level process models that resolve more processes, (2) providing a model framework that can facilitate model comparisons (for example, between our models and those of Nauels et al. (2016)), (3) exploring combined effects of key structural and parametric uncertainties, (4) explicitly demonstrating the flexibility of our framework for interchanging model components, and (5) explicitly demonstrating the utility of our model framework for informing decision analyses.

## Installation

### For the impatient

Model codes forked from the main Github repository, consistent with the [GMDD model description paper](http://www.geosci-model-dev-discuss.net/gmd-2016-303/), are available as a simple tarball from [this download server](https://download.scrim.psu.edu/Wong_etal_BRICK/). This code tarball includes everything needed to reproduce that work, as well as netCDF files containing the projections of sea-level rise from that work.

### Longer-term support

To obtain the model codes:
~~~~
git clone https://github.com/scrim-network/BRICK.git
~~~~

The calibrated parameter files are larger than we prefer to move around with the Github repository codes. The parameter files that correspond to the GMDD description paper (doi:10.5194/gmd-2016-303) are available from [this download server](https://download.scrim.psu.edu/Wong_etal_BRICK/). Of course, the ambitious user is invited to run his/her own calibrations, as detailed in the workflow below.

The following R packages are required. Within the workflow detailed below, there is a script to install all of these, so no need to copy-paste. Note that each command is issued separately so that if one throws an error, it is immediately obvious where the problem lies.
~~~~
install.packages(adaptMCMC)
install.packages(compiler)
install.packages(DEoptim)
install.packages(doParallel)
install.packages(fExtremes)
install.packages(fields)
install.packages(fMultivar)
install.packages(foreach)
install.packages(gplots)
install.packages(graphics)
install.packages(lhs)
install.packages(maps)
install.packages(methods)
install.packages(ncdf4)
install.packages(plotrix)
install.packages(pscl)
install.packages(RColorBrewer)
install.packages(sensitivity)
install.packages(sn)
install.packages(stats)
~~~~

## Workflow

### To reproduce the work of Wong, Bakker et al. (GMDD description paper)

1. Checkout the model codes, and enter
~~~~
git clone https://github.com/scrim-network/BRICK.git
cd BRICK
~~~~

2. Create the dynamic libraries necessary to run the model in Fortran. You might need to modify the `Makefile` to use your preferred Fortran compiler. Further help can be found at `BRICK/fortran/README`.
~~~~
cd fortran
make
~~~~

3. Open R and install the relevant R packages.
~~~~
cd ../calibration
R
source('BRICK_install_packages.R')
~~~~

4. Calibrate the default BRICK model configuration (DOECLIM+SIMPLE+GSIC+TE) parameters using modern data. This should not take longer than an hour or two on a modern multi-core computer.
~~~~
source('BRICK_calib_driver.R')
~~~~

5. Calibrate the DAIS parameters using paleoclimate data. This will take about 12 hours with a modern laptop.
~~~~
source('DAIS_calib_driver.R')
~~~~

6. Calibrate global mean sea-level rise model parameters using modern data and the Rahmstorf 2007 model.
~~~~
source('BRICK_calib_driver_R07.R')
~~~~

7. Calibrate the SIMPLE-GSIC parameters using modern data.
~~~~
source('BRICK_calib_driver_SIMPLE-GSIC.R')
~~~~

8. Combine modern and paleo calibration parameters, and calibrate the joint set to sea level data using rejection sampling; make hindcasts and projections of sea level; project local sea level and assess flood risks for the control model configuration. Note: if you run your own calibrations, files names must be edited to point to the correct files. By default, they point to the file names as used in the GMDD model description paper. To fully reproduce all of the experiments from that work, you must run the following script three times, where `experiment` is set to (i) c (control), (ii) e (SIMPLE-GSIC), and (iii) g (BRICK-GMSL).
~~~~
source('processingPipeline_BRICKexperiments.R')
~~~~

9. Create plots and analysis, as seen in the GMDD description paper. Note: file names must be edited to point to the correct files, if you have run your own calibrations. Also note: the directory in which you save the plots must be changed to match somewhere on your own machine.
~~~~
source('analysis_and_plots_BRICKexperiments.R')
~~~~

Note that in each of these scripts, some edits will be necessary. These will include pointing at the proper file names. The BRICK and DAIS calibration driver scripts produce calibrated parameter files with date-stamps in their names. You will need to make sure the processing pipeline script points at the current calibrated parameters files. The processing pipeline script, in turn, produces several netCDF output files from your fully calibrated BRICK parameters. These file names also include date-stamps. You will need to make sure the analysis and plotting script points at the current files. You will also need to modify the directory in which your plots will be saved.

## Code Example

Suppose you are a researcher who wishes to use the sea level projections from the BRICK GMDD model description paper in your own work. The following example will demonstrate how to use these projections and fingerprint the global sea level contributions to local mean sea level rise.

1. Checkout the model codes.
~~~~
git clone https://github.com/scrim-network/BRICK.git
~~~~

2. Open R, and navigate to the BRICK directory containing the `BRICK_LSL.R` script.
~~~~
R
setwd('~/codes/BRICK/R')
~~~~

3. Set a latitude and longitude at which you want to project local mean sea level. Here, we demonstrate using Key West, Florida (24.5551 deg N, 81.7800 deg W).
~~~~
lat <- 24.5551
lon <- -81.7800
~~~~

4. Read in the projections of the major contributions to global mean sea level for the control ensemble presented in the [BRICK model description paper](http://www.geosci-model-dev-discuss.net/gmd-2016-303/). Here, we demonstrate using RCP8.5.
~~~~
library(ncdf4)
filename.projections <- '../output_model/BRICK-model_physical_control_02Apr2017.nc'
ncdata <- nc_open(filename.projections)
slr.gsic <- ncvar_get(ncdata, 'GSIC_RCP85')
slr.gis <- ncvar_get(ncdata, 'GIS_RCP85')
slr.ais <- ncvar_get(ncdata, 'AIS_RCP85')
slr.te <- ncvar_get(ncdata, 'TE_RCP85')
t.proj <- ncvar_get(ncdata, 'time_proj')
n.ens <- length(ncvar_get(ncdata, 'ens'))
nc_close(ncdata)
~~~~

5. Use the `BRICK_LSL.R` script to project local mean sea level.
~~~~
source('BRICK_LSL.R')
lsl.proj <- brick_lsl(lat.in=lat, lon.in=lon, n.time=length(t.proj), slr_gis=slr.gis, slr_gsic=slr.gsic, slr_ais=slr.ais, slr_te=slr.te)
~~~~

6. Write output to a netCDF file
~~~~
filename.output <- '../output_model/BRICK_LSL_KeyWest_RCP85.nc'
dim.tproj <- ncdim_def('time_proj', 'years', as.double(t.proj))
dim.ensemble <- ncdim_def('ens', 'ensemble member', as.double(1:n.ens), unlim=TRUE)
lsl.rcp85 <- ncvar_def('LocalSeaLevel_RCP85', 'meters', list(dim.tproj, dim.ensemble), -999, longname = 'Local sea level (RCP85)')
outnc <- nc_create(filename.output, list(lsl.rcp85), force_v4 = TRUE)
ncvar_put(outnc, lsl.rcp85, t(lsl.proj))
nc_close(outnc)
~~~~

## Reference documentation to get new users started

/fortran/README
   * details for "make"-ing and using the dynamic libraries necessary to run the models in Fortran, with R wrapper functions

/data/README
   * description of each source of calibration/forcing data, including citations and download sources

/calibration/README_calibration_DAIS
   * describes the MCMC approach to calibrate DAIS

/calibration/README_calibration
   * describes the MCMC approach to calibrate the rest of the physical models in a coupled setting

/calibration/README_projections
   * describes the calibration approach to combine the DAIS and rest-of-model posterior parameter estimates, and use these fully calibrated parameters to make sea-level rise projections

## Contributors

Please enjoy the code and offer us any suggestions. It is our aim to make the model accessible and usable by all. We are always interested to hear about potential improvements to the model, both in the statistical calibration framework as well as the physical sub-models for climate and contributions to sea-level rise.

Questions? Tony Wong (twong@psu.edu)

## License

Copyright 2016 Tony Wong, Alexander Bakker

This file is part of BRICK (Building blocks for Relevant Ice and Climate Knowledge). BRICK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

BRICK is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with BRICK.  If not, see <http://www.gnu.org/licenses/>.
