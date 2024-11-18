#License info ----
#Copyright 2019 Province of British Columbia
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

#libraries
library(terra)
library(dplyr)
library(data.table)
library(climr)
library(sf)
library(raster)
library(tidyverse)


#pull in training data set (lat/longs)
trainpts<-read.csv("spatialdata/WNA_v13_50-200filtpts_15Nov.csv")

##pull in climate data for training pts
my_grid<-trainpts
my_grid<-dplyr::select(my_grid, -BGC)
colnames(my_grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects

#which variables do we want? 
varsl = c("CMD_sm", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
          "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", "SHM", "Tave_sm", 
          "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", "Tmin_at", "Tmin_sm", 
          "Tmin_sp", "Tmin_wt", "CMI")  #  "DD_0_sp","PPT_MJ", "PPT_JAS", "CMD.total")

## climr call- This will return the observed 1961-1990 climates for the raster grid points.
cache_clear()
gc()

#only use if need to re-run climr 
climlayer <- downscale(
  xyz = my_grid,  which_refmap = "refmap_climr",
  obs_periods = "2001_2020", 
  vars = varsl)

addVars <- function(dat) {
  dat[, PPT_MJ := PPT_05 + PPT_06]
  dat[, PPT_JAS := PPT_07 + PPT_08 + PPT_09]
  #dat[, PPT.dormant := PPT_at + PPT_wt]
  #dat[, CMD.def := pmax(0, 500 - PPT.dormant)]
  #dat[, CMDMax := CMD_07]   ## TODO: THIS IS NOT NECESSARILY CMD MAX
  dat[, CMD.total := CMD.def + CMD]
}

save(climlayer, file="trainingpts_w_clim_Final.Rdata")


