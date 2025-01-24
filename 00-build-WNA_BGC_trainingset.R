# Random Forest Model of Biogeoclimatic Units for Western North America
# Original script: Build_WNA_BGC_trainingset.Rmd by William H MacKenzie & Kiri Daust

# Updated by Deb Obrist (January 2025)

# Load packages: 
library(tidyverse)
library(terra)
library(climr)
library(reproducible) # For Cache function
library(data.table)
library(sf)
library(foreach) # for outlier removal function
library(tidymodels) # for prep() function from recipes package. 
library(themis) # for step_downsample() function
library(ranger) # For RF
library(caret) # For confusionMatrix()

# Source some functions: 
source("R/utils.R")

# Set default cache directories for the reproducible and climr packages (where intermediate results/downloaded data will be stored): 
options(reproducible.cachePath = "reproducible.cache/",
        climr.cache.path = "climr.cache/")

#### Create training points: ####
# Load in BGC polygons: 
# QUESTION 1: What version of the BGC data should I use? 
bgcs <- st_read("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/BGCProjections/Proxy_ObjectStorage/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")
# bgcs <- st_read("//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")


# And the DEM:
# QUESTION 2: Which DEM should I use?  
elev <- rast("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/BGCProjections/Proxy_ObjectStorage/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")
# elev <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")

# Reproject elev to Albers: 
elev <- project(elev, crs(bgcs))

# Define smaller study area for script building purposes. These are the extents in lat/long but we want them in Albers instead.

# Remove this later and run script on entire training area instead. 
# trainingarea <- ext(c(-125, -112, 43, 55))
studyarea <- ext(c(-123, -117, 49, 52.5))

# Create a SpatRaster to represent the extents in lat/long
dummy_raster <- rast(ext = studyarea, crs = "EPSG:4326", res = 0.1)  

# Reproject the dummy raster to Albers (EPSG:3005)
dummy_raster_albers <- project(dummy_raster, "EPSG:3005")

# Extract the reprojected extents
studyarea_albers <- ext(dummy_raster_albers)

# Crop the elev DEM and bgcs to just the smaller study area for faster execution: 
elev <- crop(elev, studyarea_albers)
bgcs <- st_crop(bgcs, studyarea_albers)

# From example in climr documentation: "Since climr is meant to be used to downscale climate variables in land, we will “clip” (set values outside the polygon to NAs) the raster using a land-only polygonRemove areas with water: 

# QUESTION 3: Do I need to do this? If so, what is a good land-only polygon to use? If not, what is the justification? (e.g., Colin mentioned that we want some overlap with ocean to make sure we hit islands etc)
# elev <- mask(elev, bgcs)

# This function makes a grid over the extent of bgcs, and fills it with a dummy variable (1L), converts to a spatial vector. It extracts elevation at the grid points using bilinear interpolation. 

# QUESTION 4: Gridsize = 2 km. Are we still happy with this? Justification from RMD is:  
# "A 2km grid seems to provide enough training points for most BGCs. Large non-vegetation land areas are excluded (lakes and glaciers primarily)."

# QUESTION 5: Do we want to do this? Or would it be simpler to just extract grid point every 800 m from DEM?
coords <- makePointCoords(bgcs, elev, gridSize = 2000) |>
  Cache()

# Need to rename lon and lat to x and y to work with subsetByExtent(), and because they are in Albers, not lat/long: 
setnames(coords, old = c("lon", "lat"), new = c("x", "y"))

# Extract BGC data from bgcs polygons and append as column to coords data:
points_sf <- st_as_sf(coords, coords = c("x", "y"), crs = 3005)
# points_sf <- st_transform(points_sf,3005) # Don't need because already 3005
bgc_att <- st_join(points_sf, bgcs)
bgc_att <- data.table(st_drop_geometry(bgc_att))

# bgc_att has 48108 unique IDs but 48112 rows. 
length(unique(bgc_att$id))
nrow(bgc_att)
bgc_att[duplicated(bgc_att$id), ] # ids 4979, 4980, 10485, and 46026 are duplicated. 

bgc_duplicates <- bgc_att[bgc_att$id %in% c(4979, 4980, 10485, 46026), ]
bgc_duplicates <- merge(bgc_duplicates, coords, by = c("id", "elev"))

# QUESTION 6: Figure out what's going on with these duplicates and what to do with them. 
# Check what's going on with these duplicates: 
checkarea <- ext(c(1651994 - 1000, 1651994 + 1000, 779989.6 - 10000, 779989.6 + 10000)) # xmin, xmax, ymin, ymax
bgcs_check <- st_crop(bgcs, checkarea)
plot(bgcs_check, xlim = c(min(checkarea[1]), max(checkarea[2])), 
     ylim = c(min(checkarea[3]), max(checkarea[4]))) 
points(x = bgc_duplicates$x[c(5, 6)], y = bgc_duplicates$y[c(5, 6)], col = "black", pch = 16)

# Remove duplicates for now: 
bgc_att <- unique(bgc_att, by = "id")

# Also remove rows where BGC is NA: 
bgc_att <- bgc_att[!is.na(bgc_att$BGC), ]

# Summarize how many points in each  BGC.   
# Need to check what the numbers are like when using the full extent. 
# NOTE - this step might not be necessary as there is a utils function to do this later on. 

# QUESTION 7: Does it make more sense to do this here and now or later, after we have climr vars? 

BGC_counts <- bgc_att[, .(Num = .N), by = .(BGC)] 
# BGC_counts <- BGC_counts[Num >= 10] # For now, retain only those where N > 10. 
bgc_att_filtered <- bgc_att[BGC %in% BGC_counts$BGC]

# Merge coords and BGC data from bgc_att: 
coords2 <- merge(coords, bgc_att_filtered, by = c("id", "elev"))

# Crops the coordinates to the study area defined above: 
# Note: coords_train will be identical to coords for now since I already cropped to the size of the smaller study area earlier but when I rerun with entire training area, it will be different. 
coords_train <- subsetByExtent(coords2, studyarea_albers)

# Make rectangular gap extents within the bounding box of the study area. 5L is the default number of gaps to create. 

# QUESTION 8: Is the idea that I will eventually make these gaps over the entire training area (i.e., BC, US, AB) to train a "final" model? 
gapextents <- makeGapExtents(studyarea_albers, 5L)

# Convert list of spatial extents into to polygons: 
gap_poly <- lapply(gapextents, vect, crs = "EPSG:3005")

# Combine all individual polygons into one spatial object. 
gap_poly <- do.call(rbind, gap_poly)

# Filters points in coords that fall within the gap polygons.
coords_gaps <- subsetByExtent(coords_train, gap_poly)

# Removes the points in coords_gaps from coords_train to produce a dataset of points that do not fall within the gaps. It does this by keeping only points in coords_train that do not match the id values in coords_gaps. 
coords_trainWgaps <- coords_train[!coords_gaps, on = "id"]

# Visualize: 
plot(elev, alpha = 0.8)

for(i in 1:5){
  plot(gapextents[[i]], add=T)
}

points(x = coords_train$x, y = coords_train$y, col = "grey50", cex = 0.001) # All points
points(x = coords_gaps$x, y = coords_gaps$y, col = "black", cex = 0.001) # Just the gaps
points(x = coords_trainWgaps$x, y = coords_trainWgaps$y, col = "white", cex = 0.001) # Everything but the gaps.

# Recombine coords_gaps and coords_trainingWgaps but with an extra column for "gap = 0 or 1" for holdouts. 
coords_gaps[, gap := "0"]
coords_trainWgaps[, gap:= "1"]

coords_all <- rbind(coords_gaps, coords_trainWgaps)

#### Get climate variables: ####
# coords_all must be in lat/long to work with climr. First, make it into a SpatVector: 
coords_spat <- vect(coords_all, geom = c("x", "y"), crs = "EPSG:3005")

# Reproject to lat/long: epgs 4326:
coords_spat_latlong <- project(coords_spat, "EPSG:4326")

# Extract the transformed coordinates (longitude and latitude)
coords_latlong <- as.data.table(geom(coords_spat_latlong))

# Add the transformed lon and lat columns back to coords_all: 
coords_all[, c("lon", "lat") := .(coords_latlong$x, coords_latlong$y)]

# coords_all must have the following column names for climr: id, lon, lat, elev: 
coords_all <- coords_all %>% 
  #  rename(lon = x, lat = y) %>% 
  select(id, lon, lat, elev, BGC, gap, x, y)

# Define variables needed: 

# QUESTION 9: Does it make sense to pull all variables we might want at once, and subset for model later? If so, what are all the relevant variables? 

# First, just simply PPT, Tmax, Tmin: 
vars_simple <- c("PPT", "Tmax", "Tmin")

# QUESTION 10: These are the variables previously selected. Colin and Kiri mentioned the possibility of testing three subsets of climate variables: seasonal or monthly PPT/Tmax/Tmin, a pairwise selection, and an "ecologically relevant" subset. Is this the ecologically relevant subset? If not, what should be included in that instead? (Note: check AddVars() function to add more.)
# Define a more complex set: 
vars_more <- c("DD5", "DD_0_at", "DD_0_wt", "PPT05", "PPT06", "PPT07", "PPT08",
               "PPT09", "CMD", "PPT_at", "PPT_wt", "CMD07", "SHM", "AHM", "NFFD", "PAS", "CMI")

# QUESTION 10: Is there existing code for the "pairwise" variable selection? 01-23-25 - https://github.com/bcgov/Build_WNA_BGC_model/blob/development/R/LocalFeatureSelection.R. 

# QUESTION 11: Lots of questions with this prompt - obs_periods, gcm_periods, gcms, ssps, max_run? 01-23-25 - Discussed with Kiri - This is all that is required for now. 

# Pull data from climr: 
clim_vars <- downscale(
  xyz = coords_all,
  which_refmap = "refmap_climr", 
  return_refperiod = TRUE, # Also return the 1961-1990 normals period. 
  vars = vars_simple,
  cache = TRUE)|>
  Cache()

# Subset coords_all to include only rows where the id column matches an id in clim_vars: 
# QUESTION 12: I'm not sure exactly why we need to do this - I guess in case there are some cases where climr didn't have data for all coordinates? Is that possible? 
coords_all <- coords_all[clim_vars[, .(id)], on = "id", nomatch = 0L] 

# Remove duplicates:
coords_all  <- coords_all %>% 
  distinct()

#### Assess climate variability within BGCs: ####
# First, add long, lat, elevation, BGC, x, y, and gap back in: 
trainData <- left_join(clim_vars, coords_all, relationship = "many-to-many") %>%
  distinct()

# For now, just look at the reference period: 
trainData <- trainData[PERIOD == "1961_1990"]

# QUESTION 13: Which BGCs should be included or not? Should they be removed based on number of points per BGC, or some metric from the climr data (e.g., CV for parameters within BGCs?) 

# Define bad BGCs and remove them: (These were selected in the RMarkdown script but I'm not sure why.) 
# badbgcs <- c("BWBSvk", "ICHmc1a", "MHun", "SBSun", "ESSFun", "SWBvk","MSdm3","ESSFdc3", "IDFdxx_WY", "MSabS", "FGff", "JPWmk_WY" )#, "ESSFab""CWHws2", "CWHwm", "CWHms1" , 
# trainData_bad <- trainData[BGC %in% badbgcs,]

# Set alpha for removal of outliers (2.5% = 3SD): 
# Question 14: Is this (inc. alpha of 0.025) standard practice? 
trainData <- removeOutlier(as.data.frame(trainData), alpha = .025, vars = vars_simple) |>
  Cache()

# QUESTION 14: Ok here it looks like we're removing BGCs where the number of points is less than 30. So what made the "bad" ones bad above? 

# Remove very small sample BGC units:
trainData <- rmLowSampleBGCs(trainData) |>
  Cache()

# QUESTION 15: How exactly does this work? Just ensures that at most, larger BGCs have at most 90x as many rows as smallest BGC? And randomly selects rows of that to keep? Should we do a sensitivity analysis here? 
# Subsample "oversampled" BGCs: 
dataBalance_recipe <- recipe(BGC ~ ., data =  trainData) |>
  step_downsample(BGC, under_ratio = 90) |>  ## subsamples "oversampled" BGCs
  prep()

## extract data.table
trainData_balanced <- dataBalance_recipe |>
  juice() |>
  as.data.table()

# Check numbers of BGCs: 
BGC_Nums <- trainData_balanced[,.(Num = .N), by = BGC]   

# Train ranger random forest model: 
trainData_balanced[, BGC := as.factor(BGC)]

cols <- c("BGC", vars_simple)

# QUESTION 16: How to decide on parameters here? Also need to figure out why it works when I manually separate holdout(gap) vs non-gap data but when I run within the function it crashes R. 
BGCmodel_full <- ranger(
  BGC ~ .,
  data = trainData_balanced[, ..cols],
  num.trees = 501,
  splitrule =  "extratrees",
  mtry = 2,
  min.node.size = 2,
  importance = "permutation",
 # case.weights = trainData_balanced$gap,
 # holdout = TRUE,
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
  
) |>
  Cache()

beepr::beep()

trainData_balanced_Wgaps <- trainData_balanced[gap == 1]

BGCmodel_Wgaps <- ranger(
  BGC ~ .,
  data = trainData_balanced_Wgaps[, ..cols],
  num.trees = 501,
  splitrule =  "extratrees",
  mtry = 2,
  min.node.size = 2,
  importance = "permutation",
  # case.weights = trainData_balanced$gap,
  # holdout = TRUE,
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
  
) |>
  Cache()
beepr::beep()

# QUESTION 17: How do we diagnose the models? What are we "happy" with? Also, what output specifically is required for CCISS? Accuracy/Precision/Recall/F1 score/AUC PR/ROC AUC/etc. 
# Check the models: 
conf_matrix_full <- caret::confusionMatrix(data = predictions(BGCmodel_full),
                reference = trainData_balanced$BGC)

conf_matrix_Wgaps <- caret::confusionMatrix(data = predictions(BGCmodel_Wgaps),
                       reference = trainData_balanced_Wgaps$BGC)

# Metrics (have not run)
precision <- conf_matrix_full$byClass["Precision"]
recall <- conf_matrix_full$byClass["Recall"]
f1 <- conf_matrix_full$byClass["F1"]
balanced_accuracy <- conf_matrix_full$byClass["Balanced Accuracy"]

print(BGCmodel_full) # OOB prediction error: 52.91%
print(BGCmodel_Wgaps) # OOB prediction error: 50.39%

# Save predictions (also need these for leaflet script): 
predictions_full <- predict(BGCmodel_full, data = trainData_balanced)$predictions
predictions_Wgaps <- predict(BGCmodel_Wgaps, data = trainData_balanced_Wgaps)$predictions
