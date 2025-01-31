# Random Forest Model of Biogeoclimatic Units for Western North America
# Original script: Build_WNA_BGC_trainingset.Rmd by William H MacKenzie & Kiri Daust

# Updated by Deb Obrist (January 2025)
rm(list = ls())

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
library(beepr)

# Source some functions: 
source("R/utils.R")

# Set default cache directories for the reproducible and climr packages (where intermediate results/downloaded data will be stored): 
options(reproducible.cachePath = "reproducible.cache/",
        climr.cache.path = "climr.cache/")

#### Create training points: ####
# Load in BGC polygons: 
# TO DO: Update this with v13 once available. 
bgcs <- st_read("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_15Nov2024.gpkg")
# bgcs <- st_read("//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")

# And the DEM 
# TO DO: Update to 30 m when finalized:
elev <- rast("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/BGCProjections/Proxy_ObjectStorage/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")
# elev <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")

# Reproject elev to Albers: 
elev <- project(elev, crs(bgcs))

# Define smaller test study area. These are the extents in lat/long but we want them in Albers instead.

# TO DO: Remove this later and run script on entire training area instead. 
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

# TO DO: Find out if this is necessary. If so, what is a good land-only polygon to use? If not, what is the justification? (e.g., Colin mentioned that we want some overlap with ocean to make sure we hit islands etc)

# From example in climr documentation: "Since climr is meant to be used to downscale climate variables in land, we will “clip” (set values outside the polygon to NAs) the raster using a land-only polygonRemove areas with water: 

# elev <- mask(elev, bgcs)


# TO DO: Decide/test whether it makes a difference if we randomly sample a set number of points (if so, how many) per BGC, or if we do a balanced sampling approach. If balanced, what kind of approach? 

# This function makes a grid over the extent of bgcs, and fills it with a dummy variable (1L), converts to a spatial vector. It extracts elevation at the grid points using bilinear interpolation. 
coords <- makePointCoords(bgcs, elev, gridSize = 2000) |>
  Cache()

# Need to rename lon and lat to x and y to work with subsetByExtent(), and because they are in Albers, not lat/long: 
setnames(coords, old = c("lon", "lat"), new = c("x", "y"))

# Extract BGC data from bgcs polygons and append as column to coords data:
points_sf <- st_as_sf(coords, coords = c("x", "y"), crs = 3005)
# points_sf <- st_transform(points_sf,3005) # Don't need because already 3005
bgc_att <- st_join(points_sf, bgcs)
bgc_att <- data.table(st_drop_geometry(bgc_att))

# bgc_att has 48108 unique IDs but 48109 rows. 
length(unique(bgc_att$id))
nrow(bgc_att)
bgc_att[duplicated(bgc_att$id), ] # ids 4979 is duplicated in v13

bgc_duplicates <- bgc_att[bgc_att$id %in% c(4979), ]
bgc_duplicates <- merge(bgc_duplicates, coords, by = c("id", "elev"))
bgc_duplicates_sf <- st_as_sf(bgc_duplicates, coords = c("x", "y"), crs = 3005)

# TO DO: Figure out what's going on with these duplicates and what to do about them. 
checkarea <- ext(c(1583996 - 2500, 1583996 + 2500, 827812.8 - 2500, 827812.8 + 2500)) # xmin, xmax, ymin, ymax
bgcs_check <- st_crop(bgcs, checkarea)
plot(st_geometry(bgcs_check), col = as.factor(bgcs_check$BGC))
plot(st_geometry(bgc_duplicates_sf), col = "black", pch = 16, add = TRUE)
text(st_coordinates(bgc_duplicates_sf), labels = bgc_duplicates_sf$id, cex = 0.7, pos = 3, col = "black")

# Remove this duplicate for now: 
bgc_att <- unique(bgc_att, by = "id")

# Also remove rows where BGC is NA: 
nrow(bgc_att[is.na(bgc_att$BGC)]) # 70
bgc_att <- bgc_att[!is.na(bgc_att$BGC), ]

# Merge coords and BGC data from bgc_att: 
coords <- merge.data.table(coords, bgc_att, by = c("id", "elev"))

# Crops the coordinates to the study area defined above: 
# Note: coords_train will be identical to coords for now since I already cropped to the size of the smaller study area earlier but when I rerun with entire training area, it will be different. 
coords_train <- subsetByExtent(coords, studyarea_albers)

# Make rectangular gap extents within the bounding box of the study area. 5L is the default number of gaps to create. 
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

# points(x = coords_train$x, y = coords_train$y, col = "grey50", cex = 0.001) # All points
points(x = coords_gaps$x, y = coords_gaps$y, col = "black", cex = 0.001) # Just the gaps
points(x = coords_trainWgaps$x, y = coords_trainWgaps$y, col = "white", cex = 0.001) # Everything but the gaps.

# Recombine coords_gaps and coords_trainingWgaps but with an extra column for "gap = yes or no" for holdouts. 
# coords_gaps[, gap := "yes"]
# coords_trainWgaps[, gap:= "no"]
set(coords_gaps, j = "gap", value = "yes")
set(coords_trainWgaps, j = "gap", value = "no")

# coords_all <- rbind(coords_gaps, coords_trainWgaps)
coords_all <- rbindlist(list(coords_gaps, coords_trainWgaps))

#### Local feature selection: ####
# Copy code on local feature selection here

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
setcolorder(coords_all, c("id", "lon", "lat", "elev", "BGC", "gap", "x", "y"))

# Pull data from climr:
clim_vars <- downscale(
  xyz = coords_all,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars(),
  cache = TRUE)|>
  Cache()

# Make an explicit copy of clim_vars to avoid modifying the original object by reference
clim_vars_all <- copy(clim_vars)

# Add the additional variables not available in climr. Use this for "kitchen-sink" scenario of all climate variables.
ccissr::addVars(clim_vars_all)

# Subset coords_all to include only rows where the id column matches an id in clim_vars: 
# QUESTION: I'm not sure exactly why we need to do this - I guess in case there are some cases where climr didn't have data for all coordinates? Is that possible? 
coords_all <- coords_all[clim_vars_all[, .(id)], on = "id", nomatch = 0L] 

# Remove duplicates:
coords_all  <- coords_all %>% 
  distinct()

# Bring coordinate data together with clim_vars_all data:
trainData <- left_join(clim_vars_all, coords_all) %>%
  distinct()

#### Remove outliers and filter BGCs: ####
# Set alpha for removal of outliers (2.5% = 3SD): 
# TO DO: Determine how sensitive the results are to this alpha value. 
trainData <- removeOutlier(as.data.frame(trainData), alpha = .025, vars = vars_all) |>
  Cache()

# To do: Decide how and if to filter BGCs. Start with not filtering them, document reasons why to or not to do it. 

# Previous round of bad BGCs: 
# Define bad BGCs and remove them: (These were selected in the RMarkdown script but I'm not sure why.) 
# badbgcs <- c("BWBSvk", "ICHmc1a", "MHun", "SBSun", "ESSFun", "SWBvk","MSdm3","ESSFdc3", "IDFdxx_WY", "MSabS", "FGff", "JPWmk_WY" )#, "ESSFab""CWHws2", "CWHwm", "CWHms1" , 
# trainData_bad <- trainData[BGC %in% badbgcs,]

# Remove very small sample BGC units (default cutoff = 30): SKIP FOR NOW. 
# trainData <- rmLowSampleBGCs(trainData) |>
#   Cache()

# # TO DO: Figure out if this is necessary/the best way to do it. It ensures that at most, larger BGCs have at most 90x as many rows as smallest BGC.
# BGCs_pre_downsample <- trainData %>% 
#   dplyr::group_by(BGC) %>% 
#   dplyr::summarize (n = n()) %>% 
#   dplyr::arrange(n)
# 
# # Subsample "oversampled" BGCs. This cuts them off at 90, if all BGCs are left in the sample, including those with only 1 point. 
# dataBalance_recipe <- recipe(BGC ~ ., data =  trainData) |>
#   step_downsample(BGC, under_ratio = 90) |>  ## subsamples "oversampled" BGCs
#   prep()

# Extract the data.table of balanced points:
# trainData_balanced <- dataBalance_recipe |>
#   juice() |>
#   as.data.table()

# BGCs_post_downsample <- trainData_balanced %>% 
#   dplyr::group_by(BGC) %>% 
#   dplyr::summarize (n = n()) %>% 
#   dplyr::arrange(n)

#### Assess climate variability within BGCs: ####
# Figure out if more BGCs need to be removed due to high variability in combination with small sample sizes. 

#### Set up climate variable combinations: ####
# First, just seasonal PPT, Tmax, Tmin: 
vars_simple <- c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", 
                 "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", 
                 "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt")

# Set selected through local feature selection, with features selected objectively. Every BGC is evaluated based on all BGCs that touch it. The final set of climate variables is the total list of each variable that was most important in each BGC. 
# TO DO: Set this up. 
vars_LFS <- c("")

# Expert set selected by Will MacKenzie in previous iteration of CCISS (also based on local feature selection): 
vars_expert <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", 
                 "EXT", "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", 
                 "SHM", "Tave_sm", "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", 
                 "Tmin_at", "Tmin_sm", "Tmin_sp", "Tmin_wt", "CMI", "PPT_MJ", 
                 "PPT_JAS", "CMD.total")

# Kitchen sink scenario: All vars.  
vars_all <- c(list_vars(), "PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total", "DD_delayed")

# Subset trainData to include the variables in different combinations of climr data: 
trainData <- as.data.table(trainData)
trainData[, BGC := as.factor(BGC)]

trainData_simple <- trainData[, c("id", "PERIOD", "gap", "x", "y", "lat", "lon", "elev", "BGC", ..vars_simple)]
trainData_expert <- trainData[, c("id", "PERIOD", "gap", "x", "y", "lat", "lon", "elev", "BGC", ..vars_expert)]
trainData_all <- trainData[, c("id", "PERIOD", "gap", "x", "y", "lat", "lon", "elev", "BGC", ..vars_all)]

#### Train ranger random forest model: ####
cols_simple <- c("BGC", vars_simple)
cols_expert <- c("BGC", vars_expert)
cols_all <- c("BGC", vars_all)

# Train model with simple variables, on points from the entire study area: 
BGCmodel_full_simple <- ranger(
  BGC ~ .,
  data = trainData_simple[, ..cols_simple],
  num.trees = 501,
  splitrule =  "extratrees",
  min.node.size = 2,
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
  
) |>
  Cache()
beepr::beep()

# Train model with simple variables, on points from area outside of the gaps: 
trainData_simple_Wgaps <- trainData_simple[gap == "no"]

BGCmodel_Wgaps_simple <- ranger(
  BGC ~ .,
  data = trainData_simple_Wgaps[, ..cols_simple],
  num.trees = 501,
  splitrule =  "extratrees",
  min.node.size = 2,
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
) |>
  Cache()
beepr::beep()

# Check the models: 
conf_matrix_full <- caret::confusionMatrix(data = predictions(BGCmodel_full_simple),
                                           reference = trainData_simple$BGC)
conf_matrix_Wgaps <- caret::confusionMatrix(data = predictions(BGCmodel_Wgaps_simple),
                       reference = trainData_simple_Wgaps$BGC)

print(BGCmodel_full_simple) # OOB prediction error: 28.98%
print(BGCmodel_Wgaps_simple) # OOB prediction error: 29.13%

#### Predictions: ####
# Generate predictions based on the entire study area raster: 
# First, reproject the raster to lat/long: 
elev_latlong <- project(elev, "EPSG:4326")

# Reproject the polygon to match the raster CRS
gap_poly_latlong <- terra::project(gap_poly, "EPSG:4326")

# Crop and mask in the reprojected space
elev_gaps <- crop(elev_latlong, ext(gap_poly_latlong))
elev_gaps <- mask(elev_gaps, gap_poly_latlong)

# Extract the transformed coordinates (longitude and latitude)
elev_latlong <- as.data.frame(elev_latlong, cells = TRUE, xy = TRUE)
elev_gaps <- as.data.frame(elev_gaps, cells = TRUE, xy = TRUE)

# Add the transformed lon and lat columns back to coords_all: 
# coords_all[, c("lon", "lat") := .(coords_latlong$x, coords_latlong$y)]

# Convert the raster to a dataframe: 
colnames(elev_latlong) <- c("id", "lon", "lat", "elev")
colnames(elev_gaps) <- c("id", "lon", "lat", "elev")

# Get climr data for these raster coordinates:
clim_vars_preds <- downscale(
  xyz = elev_latlong,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars(),
  cache = TRUE)|>
  Cache()

# Predictions based on full model and gap removed models: 
preds_full_simple <- predict(BGCmodel_full_simple, data = clim_vars_preds)
preds_Wgaps_simple <- predict(BGCmodel_Wgaps_simple, data = clim_vars_preds)







# Save into trainData_balanced_gaps testing dataframe:
trainData_balanced_gaps$predictions_full <- as.character(predictions_full$predictions)
trainData_balanced_gaps$predictions_Wgaps <- as.character(predictions_Wgaps$predictions)

# Check how well each model performed at predicting gaps:  
conf_matrix_testgaps_full <- caret::confusionMatrix(
  data = factor(trainData_balanced_gaps$predictions_full, levels = levels(trainData_balanced_gaps$BGC)),
  reference = factor(trainData_balanced_gaps$BGC, levels = levels(trainData_balanced_gaps$BGC))
)

conf_matrix_testgaps_Wgaps <- caret::confusionMatrix(
  data = factor(trainData_balanced_gaps$predictions_Wgaps, levels = levels(trainData_balanced_gaps$BGC)),
  reference = factor(trainData_balanced_gaps$BGC, levels = levels(trainData_balanced_gaps$BGC))
)

# Print accuracy of each on new data: 
accuracy_full_on_gaps <- mean(trainData_balanced_gaps$predictions_full == trainData_balanced_gaps$BGC)
print(paste("Accuracy on new data:", accuracy_full_on_gaps))

accuracy_Wgaps_on_gaps <- mean(trainData_balanced_gaps$predictions_Wgaps == trainData_balanced_gaps$BGC)
print(paste("Accuracy on new data:", accuracy_Wgaps_on_gaps))

# Look at feature importance: 
importance_scores_Wgaps <- BGCmodel_Wgaps$variable.importance
importance_scores_Wgaps <- sort(importance_scores_Wgaps, decreasing = TRUE)

# Create a data frame for plotting
importance_df <- data.frame(
  Feature = names(importance_scores_Wgaps),
  Importance = importance_scores_Wgaps
)

# Plot
ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Feature Importance",
       x = "Features",
       y = "Importance") +
  theme_minimal()
