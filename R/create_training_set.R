## CREATE TRAINING SET 
#source('R/utils.R') #custom functions

library(terra)
library(reproducible)
library(dplyr)
library(data.table)

## this is for testing, otherwise make NULL
studyArea <- vect(ext(c(-125, -112, 43, 55)), crs = "EPSG:4326")

#pull in bgc polygons from OS 
#bgcs <- vect("//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")  ## for poly validation if need be

#CGC locally
bgcs<-vect("spatialdata/WNA_BGC_v13_13Nov2024_2.gpkg")
gc()
studyArea <- project(studyArea, y = crs(bgcs, proj = TRUE))

bgcsdf<-as_tibble(bgcs, xy=T)%>%group_by(BGC)%>%distinct(.)

#spbal for balanced acceptance sampling
library(spbal)
#create sf object
bgcs_sf <- sf::st_as_sf(bgcs) #this needs to be the dissolved version 

#create named number vector for each bgc 
n_samples<-select(bgcsdf, BGC)
n_samples<-as.vector(t(n_samples))
numpts<-c(rep(5, 395)) #start with 5 each 
#custom_numpts <- c(10, 20, 30, 40....) #update based on size/climate of each BGC- TO DO
n_samples <- setNames(numpts, n_samples)

#run BAS
set.seed(511)
result <- spbal::BAS(shapefile = bgcs_sf,
                     n = n_samples,
                     boundingbox = studyArea, 
                     stratum= "BGC" )


