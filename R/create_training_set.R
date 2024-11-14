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

## CREATE TRAINING SET----
#libraries
library(terra)
library(reproducible)
library(dplyr)
library(data.table)
library(climr)
library(sf)
library(raster)

## this is for testing, otherwise make NULL
studyArea <- vect(ext(c(-125, -112, 43, 55)), crs = "EPSG:4326")

#pull in polygons
bgcs <- vect("spatialdata/WNA_BGC_v13_13Nov2024_2.gpkg")
gc()
studyArea <- project(studyArea, y = crs(bgcs, proj = TRUE))
#look at data
#bgcsdf <- as_tibble(bgcs, xy = TRUE) %>% group_by(BGC) %>% distinct(.)

#spbal for balanced acceptance sampling----
#NOT working 11/14/24
library(spbal)
#create sf object
bgcs_sf <- sf::st_as_sf(bgcs) #this needs to be the dissolved version

#create named number vector for each bgc
n_samples <- select(bgcsdf, BGC)
n_samples <- as.vector(t(n_samples))
numpts <- c(rep(5, 395)) #start with 5 each
#custom_numpts <- c(10, 20, 30, 40) #update based on size/climate of each BGC
n_samples <- setNames(numpts, n_samples)

#run BAS- 
set.seed(511)
result <- spbal::BAS(shapefile = bgcs_sf,
                     n = n_samples,
                     boundingbox = studyArea,
                     stratum = "BGC")

#pull in training set from Qgis-----
#200 pts per BGC- Will M. 11/14/24
trainpts<- read.csv("spatialdata/WNA_v13_200_rndpts.csv")
  
##assess climate variability within BGCs---- 

#Use DEM to call in downscaled climr data
dem <- rast("spatialdata/WNA_DEM_4326_clipped.tif")

## convert the DEM to a data.frame
my_grid <- as.data.frame(dem, cells = TRUE, xy = TRUE)
colnames(my_grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects
gc()

#which variables do we want? 
varsl = c("Tmax","Tmin", "PPT")

## climr call- This will return the observed 1961-1990 climates for the raster grid points.
cache_clear()
climlayer <- downscale(
  xyz = my_grid,  which_refmap = "refmap_climr",
  #obs_periods = "2001_2020", 
  vars = varsl)

#extract climate values 
trainpts0<-select(trainpts,xcoord, ycoord)
climdat<-extract(climlayer, trainpts0)
#climdat<-extract(dem, trainpts0) #try first with DEM 

#merge back with BGC info 
climdat<-cbind(trainpts0, climdat)
rm(trainpts0)
trainpts<-left_join(trainpts, climdat)
#check that elev matches and join worked 
#plot(trainpts$Elev1, trainpts$WNA_DEM_3005_clipped)

