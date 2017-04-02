# BRICK v0.1.1

## Synopsis

BRICK (**B**uilding blocks for **R**elevant **I**ce and **C**limate **K**nowledge) is a modular semi-empirical modeling framework to simulate global temperature and sea-level rise. In the default model configuration, first, global mean surface temperature and ocean heat uptake are simulated by DOECLIM. Changes in global mean surface temperature drive changes in global mean sea level (GMSL). The contributions to GMSL from the Greenland and Antarctic ice sheets, thermal expansion, and glaciers and ice caps are simulated. 

This repository contains the BRICK model source code and analysis scripts.

It contains the sub-models
* DOECLIM (climate)
* DAIS (Antarctic ice sheet)
* SIMPLE (Greenland ice sheet)
* GSIC-MAGICC (Glaciers and small ice caps)
* TE (thermal expansion)
* GMSL (Rahmstorf 2007, global mean sea level emulator)
* Van Dantzig flood risk assessment (Van Dantzig, 1956)

## Directory structure

.
:  BRICK "home" directory

./calibration/
: R scripts related to the (pre-/post-)calibration of the physical models, including reading data, likelihood functions, and sub-model stand-alone calibration drivers (some of these may be antiquated)

./data/
: data for radiative forcing (DOECLIM) and calibration

./fortran/
: Fortran versions of all physical models are available, and are wrapped in R calling functions in ./fortran/R/

./output_model/
: physical model output (i.e., temperature, sea-level rise)

./output_calibration/
: statistical model output (i.e., posterior parameter values)

./R/
: physical models in R

## Motivation

The motivation for the BRICK model is detailed in the [GMDD model description paper](http://www.geosci-model-dev-discuss.net/gmd-2016-303/):

> Here we introduce BRICK v0.1 ("Building blocks for Relevant Ice and Climate Knowledge"), a new model framework that focuses on accessibility, transparency, and flexibility while maintaining, as much as possible, the computational efficiency that make simple models so appealing. There is a wide range of potential applications for such a model. A simple framework enables uncertainty quantification via statistical calibration approaches (Higdon et al., 2004; Kennedy and O’Hagan, 2001), which would be infeasible with more computationally expensive models. A transparent modeling framework enables communication between scientists as well as communication with stakeholders. This leads to potential application of the model framework in decision support and education (Weaver et al., 2013). The present work expands on previous studies by (1) providing a platform of simple, but mechanistically motivated sea-level process models that resolve more processes, (2) providing a model framework that can facilitate model comparisons (for example, between our models and those of Nauels et al. (2016)), (3) exploring combined effects of key structural and parametric uncertainties, (4) explicitly demonstrating the flexibility of our framework for interchanging model components, and (5) explicitly demonstrating the utility of our model framework for informing decision analyses.

## Installation

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

**TONY TODO**

1. Create the dynamic libraries necessary to run the model in Fortran.
2. Directions on "make"-ing the file is found at /fortran/README
3. ./calibration/BRICK_calib_driver.R
    * calibrate the DOECLIM+SIMPLE+GSIC+TE parameters using modern data
4. ./calibration/DAIS_calib_driver.R
    * calibrate the DAIS parameters using paleoclimate data
5. ./calibration/BRICK_calib_driver_R07.R
    * calibrate Global sea-level rise parameters using modern data and the Rahmstorf 2007 model
6. ./calibration/BRICK_calib_driver_SIMPLE-GSIC.R
    * calibrate the SIMPLE-GSIC parameters using modern data
7. ./calibration/processingPipeline_BRICKexperiments.R
    * combine modern and paleo calibration parameters, and calibrate the joint set to sea level data using rejection sampling,
    * make hindcasts and projections of sea level
    * project local sea level and assess flood risks
    * NOTE: files names must be edited to point to the correct files.
      * Run control: set experiments to c
      * Run w/ SIMPLE-GSIC: set experiments to e
      * Run R07: set experiments to g
8. ./calibration/analysis_and_plots_BRICKexperiments.R
    * create plots and analysis for the Wong, Bakker et al (2016) figures
    * NOTE: files names must be edited to point to the correct files.

Note that in each of these scripts, some edits will be necessary. These will include pointing at the proper file names. The BRICK and DAIS calibration driver scripts produce calibrated parameter files with date-stamps in their names. You will need to make sure the processing pipeline script points at the current calibrated parameters files. The processing pipeline script, in turn, produces several netCDF output files from your fully calibrated BRICK parameters. These file names also include date-stamps. You will need to make sure the analysis and plotting script points at the current files. You will also need to modify the directory in which your plots will be saved.

## Code Example

**TONY TODO**

Show what the library does as concisely as possible, developers should be able to figure out **how** your project solves their problem by looking at the code example. Make sure the API you are showing off is obvious, and that your code is short and concise.

## Reference documentation to get new users started

/fortran/README
: details for "make"-ing and using the dynamic libraries necessary to run the models in Fortran, with R wrapper functions

/data/README
: description of each source of calibration/forcing data, including citations and download sources

/calibration/README_calibration_DAIS
: describes the MCMC approach to calibrate DAIS

/calibration/README_calibration
: describes the MCMC approach to calibrate the rest of the physical models in a coupled setting

/calibration/README_projections
: describes the calibration approach to combine the DAIS and rest-of-model posterior parameter estimates, and use these fully calibrated parameters to make sea-level rise projections

## Tests

**TONY TODO**

Describe and show how to run the tests with code examples.

## Contributors

Please enjoy the code and offer us any suggestions. It is our aim to make the model accessible and usable by all. We are always interested to hear about potential improvements to the model, both in the statistical calibration framework as well as the physical sub-models for climate and contributions to sea-level rise. 

Questions? Tony Wong (twong@psu.edu)

## License

Copyright 2016 Tony Wong, Alexander Bakker
This file is part of BRICK (Building blocks for Relevant Ice and Climate Knowledge). BRICK is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

BRICK is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with BRICK.  If not, see <http://www.gnu.org/licenses/>.
