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
#only use if need to re-run climr---- 
##pull in climate data for training pts
my_grid<-trainpts
my_grid<-dplyr::select(my_grid, -BGC)
colnames(my_grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects

#which variables do we want? 
varsl = c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
          "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", "SHM", "Tave_sm", 
          "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", "Tmin_at", "Tmin_sm", 
          "Tmin_sp", "Tmin_wt", "CMI")  #  ""PPT_MJ", "PPT_JAS", "CMD.total")

## climr call- This will return the observed 1961-1990 climates for the raster grid points.
cache_clear()
gc()

#only use if need to re-run climr 
#climlayer <- downscale(
#  xyz = my_grid,  which_refmap = "refmap_climr",
#  obs_periods = "2001_2020", 
#  vars = varsl)

addVars <- function(dat) {
  dat[, PPT_MJ := PPT_05 + PPT_06]
  dat[, PPT_JAS := PPT_07 + PPT_08 + PPT_09]
  #dat[, PPT.dormant := PPT_at + PPT_wt]
  #dat[, CMD.def := pmax(0, 500 - PPT.dormant)]
  #dat[, CMDMax := CMD_07]   ## TODO: THIS IS NOT NECESSARILY CMD MAX
  dat[, CMD.total := CMD.def + CMD]
}

save(climlayer, file="trainingpts_w_clim_Final.Rdata")

#fit mods----
#load climate dataset
load(file="trainingpts_w_clim_Final.Rdata")
model_data<-climlayer
rm(climlayer)
model_data<-subset(model_data, PERIOD=="1961_1990")
#merge with BGC info 
trainpts<-rename(trainpts, id=X)
model_data<-left_join(model_data, trainpts)
model_data$BGC<-as.factor(model_data$BGC)
model_data<-dplyr::select(model_data, -id, -lat, -lon, -elev, -PERIOD)
model_data<-na.omit(model_data) #removes 47 obs with NAs in CMI
#make sure still >50 obs for all BGCs
BGSpts<-group_by(model_data, BGC)%>%summarise(ct=n())
min(BGCpts$ct)

# Fit RF models- response and probability  
BGC_RFresp<- ranger::ranger(BGC~ .,data = model_data, mtry= 5, classification = T, probability = F, 
                            num.trees = 501, splitrule =  "extratrees", min.node.size = 2,
                            importance = "permutation",write.forest = TRUE) 
save(BGC_RFresp, file="model_output/BGC_RFresp.Rdata")# CGC local only 

#they are in here 
#F:/OneDrive - Government of BC/WNA_BGC/Trained_Models/
  
BGC_RFprob<- ranger::ranger(BGC~ .,data = model_data, mtry= 5, classification = T, probability = F, 
                             num.trees = 501, splitrule =  "extratrees", min.node.size = 2,
                             importance = "permutation",write.forest = TRUE) 
save(BGC_RFresp, file="model_output/BGC_RFprob.Rdata")# CGC local only 
#they are in here 
#F:/OneDrive - Government of BC/WNA_BGC/Trained_Models/


#error estimates 
#oob estimates 
print(BGC_RFresp)#0.24
print(BGC_RFprob) #0.24

#confusion matrices on fitted data 
cf<-BGC_RFresp$confusion.matrix
accuracy <- sum(diag(cf)) / sum(cf) #acc = 0.76 

cf2<-BGC_RFprob$confusion.matrix
accuracy <- sum(diag(cf2)) / sum(cf2) #acc = 0.76 


