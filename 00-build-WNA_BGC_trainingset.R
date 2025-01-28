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
# TO DO: Update this with v13 once available. 
bgcs <- st_read("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_15Nov2024.gpkg")
# bgcs <- st_read("//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")

# TO DO: Find 30 m DEM. (Will says best to use 30 m DEM because climr vars depend on elevation). 
# And the DEM:
elev <- rast("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_DEM_4326_clipped.tif")
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
coords <- merge(coords, bgc_att, by = c("id", "elev"))

# Crops the coordinates to the study area defined above: 
# Note: coords_train will be identical to coords for now since I already cropped to the size of the smaller study area earlier but when I rerun with entire training area, it will be different. 
coords_train <- subsetByExtent(coords2, studyarea_albers)

# Make rectangular gap extents within the bounding box of the study area. 5L is the default number of gaps to create. 

# QUESTION 8: Is the idea that I will eventually make these gaps over the entire training area (i.e., BC, US, AB) to train and check a "final" model? 
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
coords_gaps[, gap := "yes"]
coords_trainWgaps[, gap:= "no"]

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
list_vars()

# First, just seasonal PPT, Tmax, Tmin: 
vars_simple <- c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", 
                 "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", 
                 "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt")





# QUESTION 10: These are the variables previously selected. Colin and Kiri mentioned the possibility of testing three subsets of climate variables: seasonal or monthly PPT/Tmax/Tmin, a pairwise selection, and an "ecologically relevant" subset. Is this the ecologically relevant subset? If not, what should be included in that instead? (Note: check AddVars() function to add more.)
# Define a more complex set: 
vars_more <- c("DD5", "DD_0_at", "DD_0_wt", "PPT05", "PPT06", "PPT07", "PPT08",
               "PPT09", "CMD", "PPT_at", "PPT_wt", "CMD07", "SHM", "AHM", "NFFD", "PAS", "CMI")

# QUESTION 10: Is there existing code for the "pairwise" variable selection? 01-23-25 - https://github.com/bcgov/Build_WNA_BGC_model/blob/development/R/LocalFeatureSelection.R. 

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

# Run the RF model with the entire area: 
BGCmodel_full <- ranger(
  BGC ~ .,
  data = trainData_balanced[, ..cols],
  num.trees = 501,
  splitrule =  "extratrees",
  # mtry = 2, # Allow for default.
  min.node.size = 2,
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
  
) |>
  Cache()
beepr::beep()

# Run the model on area outside of the gaps: 
trainData_balanced_Wgaps <- trainData_balanced[gap == 1]

BGCmodel_Wgaps <- ranger(
  BGC ~ .,
  data = trainData_balanced_Wgaps[, ..cols],
  num.trees = 501,
  splitrule =  "extratrees",
 # mtry = 2, # Allow for default.
  min.node.size = 2,
  importance = "permutation",
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

print(BGCmodel_full) # OOB prediction error: 23.37%
print(BGCmodel_Wgaps) # OOB prediction error: 23.13%

# Generate predictions based on the coordinates in the gaps: 
# Make testing data with gap data (gap = 0): 
trainData_balanced_gaps <- trainData_balanced[gap == 0]

# Predictions based on full model and gap removed model: 
predictions_full <- predict(BGCmodel_full, data = trainData_balanced_gaps)
predictions_Wgaps <- predict(BGCmodel_Wgaps, data = trainData_balanced_gaps)

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
