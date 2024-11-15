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
library(dplyr)
library(data.table)
library(climr)
library(sf)
library(raster)
library(tidyverse)

#pull in training pts from Qgis-----
#200 pts per BGC- Will M. 11/14/24
trainpts<- read.csv("spatialdata/WNA_v13_200_rndpts_15Nov.csv")
  
##pull in climate data for training pts---- 
my_grid<-trainpts
my_grid$id<-rownames(my_grid)
my_grid<-dplyr::select(my_grid, -BGC)
colnames(my_grid) <- c( "lon", "lat", "elev", "id") # rename column names to what climr expects

#which variables do we want? 
varsl = c("Tmax","Tmin", "PPT")
## climr call- This will return the observed 1961-1990 climates for the raster grid points.
cache_clear()
gc()
climlayer <- downscale(
  xyz = my_grid,  which_refmap = "refmap_climr",
  #obs_periods = "2001_2020", 
  vars = varsl)

save(climlayer, file="trainingpts_w_clim.Rdata")

#assess climate variability within BGCs
#merge back with BGC info 
load(file="trainingpts_w_clim.Rdata")
climlayer<-left_join(climlayer, my_grid)
climlayer<-left_join(climlayer, trainpts)
#check that elev matches and join worked 
plot(climlayer$Elev1, climlayer$elev) #good 
climlayer<-dplyr::select(climlayer, -Elev1, -xcoord, -ycoord)

save(climlayer, file="trainingpts_w_climALL.Rdata")

rm(my_grid)
rm(trainpts)
gc()

#look at data
library(ggplot2)

ggplot(climlayer, aes(x=BGC, y=Tmin))+ geom_boxplot()
ggplot(climlayer, aes(x=BGC, y=Tmax))+ geom_boxplot()
ggplot(climlayer, aes(x=BGC, y=PPT))+ geom_boxplot()


#generate estimates by BGC
min(climlayer$Tmin)
CVs<-group_by(climlayer, BGC)%>%summarise(avgTmin=mean(Tmin+15), avgTmax=mean(Tmax+15), avgPpt=mean(PPT),
                                          sdTmin=sd(Tmin+15), sdTmax=sd(Tmax+15), sdPpt=sd(PPT))%>%mutate(cvTmin=sdTmin/(avgTmin), 
                                                                                                    cvTmax=sdTmax/(avgTmax), 
                                                                                                     cvPpt=sdPpt/avgPpt)
CVlong<-pivot_longer(CVs, cols = c(cvTmin, cvTmax, cvPpt), names_to = 'climvar', values_to = 'cv')%>%distinct(.)

ggplot(CVlong, aes(x=BGC, y=cv))+ geom_point()+ facet_wrap(~climvar,scales = 'free')+ geom_hline(yintercept = 0.3)


#exclusion criteria
BGC_exclude<-subset(CVs, cvPpt>0.3)
BGC_exclude<-BGC_exclude$BGC

climlayer<-mutate(climlayer, exclude=if_else(BGC %in% BGC_exclude, 'Y', 'N'))

ggplot(subset(climlayer,exclude=='Y'), aes(x=BGC, y=Tmin, fill=cvPpt, colour = cvPpt))+ geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(subset(climlayer,exclude=='Y'), aes(x=BGC, y=Tmax))+ geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(subset(climlayer,exclude=='Y'|BGC=="CWHxs"|BGC=="IDFww"), aes(x=BGC, y=PPT,fill=cvPpt, colour = cvPpt ))+ geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#add in IDFww and CWHxs for comparison
#ggplot(subset(climlayer,exclude=='Y'|BGC=="CWHxs"|BGC=="IDFww"), aes(x=BGC, y=Tmin,fill=cvTmin, colour = cvTmin ))+ geom_boxplot()+
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggplot(subset(climlayer,exclude=='Y'|BGC=="CWHxs"|BGC=="IDFww"), aes(x=BGC, y=Tmax,fill=cvTmax, colour = cvTmax ))+ geom_boxplot()+
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggplot(subset(climlayer,BGC=="IDFww"|BGC=="CWHxs"), aes(x=BGC, y=PPT,fill=cvPpt, colour = cvPpt ))+ geom_boxplot()+
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#remove actual outliers by BGC
climlayer<-mutate(climlayer, remove= case_when(BGC=="BAFAun"& PPT>=3000~'Y',
                                               BGC=="CCHun_CA"& PPT>=1000~'Y',
                                               BGC=="CDFxm_CA"& PPT>=1500~'Y',
                                               BGC=="CMAun"& PPT>=5700~'Y',
                                               BGC=="CMXdm_OR"& PPT>=1800~'Y',
                                               BGC=="CMXmm_OR"& PPT>=1800~'Y',
                                               
                                               ))

#spbal for balanced acceptance sampling----
#NOT working 11/14/24
## this is for testing, otherwise make NULL
studyArea <- vect(ext(c(-125, -112, 43, 55)), crs = "EPSG:4326")

#pull in polygons
bgcs <- vect("spatialdata/WNA_BGC_v13_13Nov2024_2.gpkg")
gc()
studyArea <- project(studyArea, y = crs(bgcs, proj = TRUE))
#look at data
#bgcsdf <- as_tibble(bgcs, xy = TRUE) %>% group_by(BGC) %>% distinct(.)

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
