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
library(tidyverse)
library(terra)
library(data.table)
library(climr)
library(sf)
library(raster)#masks dplyr select!!
library(ranger)

#load trained mods----
#multi edatope model
load(file="model_output/BGC_RFresp.Rdata")

#predict----
#predict models over the full climate space of WNA  
WNADEM<-rast("spatialdata/WNA_DEM_4326_clipped.tif")
WNADEMlow<-aggregate(WNADEM, fact=13, funct=mean) #go from 30m to ~400m resolution
rm(WNADEM)
gc()
my_grid <- as.data.frame(WNADEMlow, cells = TRUE, xy = TRUE)
colnames(my_grid) <- c("id", "lon", "lat", "elev") # rename column names to what climr expects

#select climate variables from BGC model
varsl = c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", "EXT", 
          "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", "SHM", "Tave_sm", 
          "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", "Tmin_at", "Tmin_sm", 
          "Tmin_sp", "Tmin_wt", "CMI")  #  ""PPT_MJ", "PPT_JAS", "CMD.total")

## climr call- This will return the observed 1961-1990 climates for the raster grid points.
cache_clear()
gc()

#only use if need to re-run climr 
climlayer <- downscale(
  xyz = my_grid,  which_refmap = "refmap_climr",
#  obs_periods = "2001_2020", 
  vars = varsl)

addVars <- function(dat) {
  dat[, PPT_MJ := PPT_05 + PPT_06]
  dat[, PPT_JAS := PPT_07 + PPT_08 + PPT_09]
  #dat[, PPT.dormant := PPT_at + PPT_wt]
  #dat[, CMD.def := pmax(0, 500 - PPT.dormant)]
  #dat[, CMDMax := CMD_07]   ## TODO: THIS IS NOT NECESSARILY CMD MAX
  dat[, CMD.total := CMD.def + CMD]
}

save(climlayer, file="spatialdata/400m_climlayer_WNA.Rdata")

#
#make predictions
#load trained model 
rm(WNADEMlow)
gc()
load(file="model_output/BGC_RFresp.Rdata")

predsBGC<-predict(object = BGC_RFresp,data =climlayer, type='response') 
save(predsBGC, file = "model_output/BGCmodel_preds_400m_1961-1990.Rdata")

process_predictions <- function(preds) {
  # Convert predictions to a data frame
  preds_df <- as.data.frame(preds$predictions)
  
  # Create an ID column based on row names, converting to numeric
  preds_df$id <- as.numeric(as.character(row.names(preds_df)))
  
  # Rename columns to "ypred" and "id"
  colnames(preds_df) <- c("ypred", "id")
  
  return(preds_df)
}
predsBGC.df<-process_predictions(predsBGC)

#combine with lat/lon
my_grid$id<-row_number(my_grid)
predsBGC.df<-left_join(predsBGC.df, my_grid)
predsBGC_tidy<-as_tibble(predsBGC.df)#create tidy df 

library(tidyterra)

predsBGC_tidy$BGCnum<-as.numeric(predsBGC_tidy$ypred)
predsBGC_tidy<-rename(predsBGC_tidy, BGC=ypred)

BGCplot<-ggplot()+ 
  geom_raster(data = predsBGC_tidy, aes(x = lon, y = lat, fill = BGCnum)) +
  #geom_point(data = Fd, aes(x = Longitude, y = Latitude, fill= TotalA_class), shape=21, size=3) +
  scale_fill_grass_c(palette = 'corine', limits = c(1,366))+
  #labs(title = "Fd plot data")  + xlab(" ") + ylab(" ")+ 
  # geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(legend.position='none', axis.line=element_blank(),axis.text=element_blank(), axis.ticks = element_blank())

ggpubr::ggarrange(BGCplot, legend = 'top')

#combine back with BGCs so less colors to plot 
BGC_list<-read.csv("tables/WNA_BGCs_Info_v13_1.csv")
BGC_list<-select(BGC_list, Zone, Subzone, BGC)
predsBGC_tidy<-left_join(predsBGC_tidy, BGC_list)

cols <- c(
  "dodgerblue2", "blue1", "steelblue4", "skyblue2", 'cornflowerblue', 'darkblue',
  "green4","palegreen2","green1",'darkgreen','darkolivegreen', 
  "darkturquoise",'aquamarine2' ,
  "#E31A1C", 'darkred',
  "orchid1", "deeppink1",  "#FB9A99",'deeppink4',  "maroon",
  "coral1", "#FF7F00","#FDBF6F",'orange',
  "gold1", "yellow4", "yellow3","khaki2",
  "#6A3D9A", "blueviolet", "#CAB2D6",'darkslateblue',
  "black", "gray70",  'azure3','bisque2',
  "darkorange4", 'chocolate3'
)
predsBGC_tidy<-na.omit(predsBGC_tidy)

BGCplot<-ggplot()+ 
  geom_raster(data = predsBGC_tidy, aes(x = lon, y = lat, fill = Zone)) +
  #geom_point(data = Fd, aes(x = Longitude, y = Latitude, fill= TotalA_class), shape=21, size=3) +
  #scale_fill_grass_d(palette = 'corine')+
  scale_fill_manual(values=cols)+
  #labs(title = "Fd plot data")  +
  xlab(" ") + ylab(" ")+ 
  # geom_sf(data = bcboundlow.reproj, fill = NA, color='grey')+ 
  theme_classic() + theme(legend.position='none', axis.line=element_blank(),
                          axis.text=element_blank(), axis.ticks = element_blank())

ggpubr::ggarrange(BGCplot, legend = 'left')
unique(predsBGC_tidy$Zone)
