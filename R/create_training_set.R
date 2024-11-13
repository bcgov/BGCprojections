## CREATE TRAINING SET 
#source('R/utils.R') #custom functions

library(terra)
library(reproducible)
library(tidyverse)
library(data.table)

## this is for testing, otherwise make NULL
#studyArea <- vect(ext(c(-125, -112, 43, 55)), crs = "EPSG:4326")

#pull in bgc polygons from OS 
bgcs <- vect("//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")  ## for poly validation if need be
#saveRDS(bgcs, file="spatialdata/bgc_poly.RData")

#read in bgc polygons -only on CGC local 
#bgcs<-readRDS(file="spatialdata/bgc_poly.RData")
#bgcs<-vect(bgcs)
gc()
studyArea <- project(studyArea, y = crs(bgcs, proj = TRUE))

bgcsdf<-as_tibble(bgcs, xy=T)%>%group_by(BGC)%>%distinct(.)

#spbal for balanced acceptance sampling
library(spbal)
#create sf object
bgcs_sf <- sf::st_as_sf(bgcs) #this needs to be the dissolved version 

#create named number vector for each bgc 
n_samples<-as.vector(t(bgcsdf))
numpts<-c(rep(5, 378))
#custom_numpts <- c(10, 20, 30, 40....) #update based on size/climate of each BGC- TO DO
n_samples <- setNames(numpts, n_samples)

#run BAS
set.seed(511)
result <- spbal::BAS(shapefile = bgcs_sf,
                     n = n_samples,
                     boundingbox = studyArea, 
                     stratum= "BGC" )


