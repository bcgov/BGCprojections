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

#Generate climate attributed training data-----
#pull in training pts from Qgis- 200 pts per BGC
trainpts<- read.csv("spatialdata/WNA_v13_200_rndpts_15Nov.csv")
trainpts$id<-rownames(trainpts)

##pull in climate data for training pts
my_grid<-trainpts
my_grid<-dplyr::select(my_grid, -BGC)
colnames(my_grid) <- c( "lon", "lat", "elev", "id") # rename column names to what climr expects

#which variables do we want? 
varsl = c("Tmax","Tmin", "PPT")
## climr call- This will return the observed 1961-1990 climates for the raster grid points.
cache_clear()
gc()

#only use if need to re-run climr 
#climlayer <- downscale(
#  xyz = my_grid,  which_refmap = "refmap_climr",
#  obs_periods = "2001_2020", 
#  vars = varsl)

#save(climlayer, file="trainingpts_w_clim.Rdata")

#assess climate variability within BGCs----
#merge back with BGC info 
#load(file="trainingpts_w_clim.Rdata")
#climlayer<-left_join(climlayer, my_grid)
#climlayer<-left_join(climlayer, trainpts)
#check that elev matches and join worked 
#plot(climlayer$Elev1, climlayer$elev) #good 
#climlayer<-dplyr::select(climlayer, -Elev1, -xcoord, -ycoord)

#save again
#save(climlayer, file="trainingpts_w_clim.Rdata")

rm(my_grid)
rm(trainpts)
gc()

#look at data
library(ggplot2)
load(file="trainingpts_w_clim.Rdata")

ggplot(climlayer, aes(x=BGC, y=Tmin))+ geom_boxplot()
ggplot(climlayer, aes(x=BGC, y=Tmax))+ geom_boxplot()
ggplot(climlayer, aes(x=BGC, y=PPT))+ geom_boxplot()


#generate estimates of climate variation by BGC
min(climlayer$Tmin)#add 15 to all temps so non negative, non zero CVs 
CVs<-group_by(climlayer, BGC)%>%summarise(avgTmin=mean(Tmin+15), avgTmax=mean(Tmax+15), avgPpt=mean(PPT),
                                          sdTmin=sd(Tmin+15), sdTmax=sd(Tmax+15), sdPpt=sd(PPT))%>%mutate(cvTmin=sdTmin/(avgTmin), 
                                                                                                    cvTmax=sdTmax/(avgTmax), 
                                                                                                     cvPpt=sdPpt/avgPpt)
CVlong<-pivot_longer(CVs, cols = c(cvTmin, cvTmax, cvPpt), names_to = 'climvar', values_to = 'cv')%>%distinct(.)

#plot all CVs
cv1<-ggplot(CVlong, aes(x=BGC, y=cv))+ geom_point()+ facet_wrap(~climvar,scales = 'free')+ geom_hline(yintercept = 0.3)

#pull in BGC areas to scale
BGCareas<-read.csv('spatialdata/WNA_v13_BGC_area.csv')
min(BGCareas$area)
#rescale areas by smallest 
BGCareas$area2<-BGCareas$area/min(BGCareas$area)
min(BGCareas$area2)#now =1
max(BGCareas$area2)

CVs<-left_join(CVs, BGCareas)
CVs<-mutate(CVs, cvTmin_scaled=(cvTmin/area2)*100, 
            cvTmax_scaled=(cvTmax/area2)*100, 
            cvPpt_scaled=(cvPpt/area2)*100)

CVlong2<-pivot_longer(CVs, cols = c(cvTmin_scaled, cvTmax_scaled, cvPpt_scaled), names_to = 'climvar', values_to = 'cv_scaled')%>%distinct(.)

cv2<-ggplot(CVlong2, aes(x=BGC, y=cv_scaled))+ geom_point()+ facet_wrap(~climvar,scales = 'free')+ geom_hline(yintercept = 0.002)+
ylim(0,0.02)

gridExtra::grid.arrange(cv1, cv2)

#exclusion criteria- scaled CV> 0.003
BGC_exclude<-subset(CVlong2, cv_scaled>=0.003)
BGC_exclude<-unique(BGC_exclude$BGC) #25 BGCs - all relatively small/ transitional units  - most common ESSF, MH & IDF 
 
climlayer<-left_join(climlayer, CVs)
climlayer<-mutate(climlayer, exclude=if_else(BGC %in% BGC_exclude, 'Y', 'N'))

ggplot(subset(climlayer,exclude=='Y'), aes(x=BGC, y=Tmin, fill=cvPpt_scaled, colour = cvPpt))+ geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(subset(climlayer,exclude=='Y'), aes(x=BGC, y=Tmax, fill=cvPpt_scaled, colour = cvPpt))+ geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(subset(climlayer,exclude=='Y'), aes(x=BGC, y=log(PPT),fill=cvPpt_scaled, colour = cvPpt ))+ geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

climlayer_filt<-subset(climlayer, exclude=="N")#72224 obs
unique(climlayer_filt$BGC) #370 BGCs- removed 25

#remove actual outlier training points from remaining BGCs > +/- 3 SDevs from mean 
ggplot(climlayer_filt, aes(x=BGC, y=Tmin, fill=cvPpt_scaled, colour = cvPpt))+ geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

climlayer_filt<-group_by(climlayer_filt, BGC)%>%mutate(avgTmin=mean(Tmin), avgTmax=mean(Tmax), avgPpt=mean(PPT),
                                          sdTmin=sd(Tmin), sdTmax=sd(Tmax), sdPpt=sd(PPT))%>%
  mutate(sdevsTmin=(Tmin-avgTmin)/sdTmin, 
         sdevsTmax=(Tmax-avgTmax)/sdTmax, 
         sdevsPPT=(PPT-avgPpt)/sdPpt) %>%subset(.,  sdevsPPT>-3 &  sdevsPPT<3) %>%subset(., sdevsTmin>-3 & sdevsTmin<3)%>%
         subset(., sdevsTmax>-3 & sdevsTmax<3)

#71262 obs - 962 outliers removed 

ggplot(climlayer_filt, aes(x=BGC, y=Tmin, fill=cvPpt_scaled, colour = cvPpt))+ geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#make sure no BGCs are highly under sampled- min n=50
BGCpts<-group_by(climlayer_filt, BGC)%>%summarise(ct=n())
BGCkeep<-filter(BGCpts, ct>50)#removed 3 BGCs- IDFdk_WA,MSxx_NV,MSxk_WA
BGCkeep<-unique(BGCkeep$BGC)

climlayer_filt<-subset(climlayer_filt, BGC %in% BGCkeep)

#clean up and save
climlayer_filt<- dplyr::select(climlayer_filt, -id, -exclude, -area2)
save(climlayer_filt, file="trainingpts_w_clim_FILTERED.Rdata")
training_points_filt<-dplyr::select(climlayer_filt, lon, lat, elev, BGC)%>%relocate(BGC, .before = lon) 
write.csv(training_points_filt, "spatialdata/WNA_v13_50-200filtpts_15Nov.csv" )

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
