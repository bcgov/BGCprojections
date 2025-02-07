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
library(beepr)
library(leaflet)

# Source functions: 
source("utils.R")

# Set default cache directories for the reproducible and climr packages (where intermediate results/downloaded data will be stored): 
options(reproducible.cachePath = "reproducible.cache/",
        climr.cache.path = "climr.cache/")

# Read in colours for reference of factors: 
subzones_colours_ref <- fread("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNAv13_Subzone_colours_2.csv") %>% 
  dplyr::select(!c(fid, MAP_LABEL, NSRNAME, ZONE)) %>% 
  dplyr::mutate(BGC_num = as.numeric(as.factor(BGC)))

#### Create training points: ####
# Load in BGC polygons: 
# TO DO: 
# Update this with v13 once available, source from object storage. 
bgcs <- st_read("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_15Nov2024.gpkg")
# bgcs <- st_read("//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")

# And the DEM 
# TO DO: 
# Update to 30 m when finalized, also source from object storage. 
elev <- rast("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/BGCProjections/Proxy_ObjectStorage/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")
# elev <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")

# First, reproject elev to Albers: 
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

# TO DO: 
# From example in climr documentation: "Since climr is meant to be used to downscale climate variables in land, we will “clip” (set values outside the polygon to NAs) the raster using a land-only polygonRemove areas with water: 

# elev <- mask(elev, bgcs)

# Find out if this is necessary. If so, what is a good land-only polygon to use? If not, what is the justification? (e.g., Colin mentioned that we want some overlap with ocean to make sure we hit islands etc)


# TO DO: 
# Decide/test whether it makes a difference if we randomly sample a set number of points (if so, how many) per BGC, or if we do a balanced sampling approach. If balanced, what kind of approach? 

# This function makes a grid over the extent of bgcs, and fills it with a dummy variable (1L), converts to a spatial vector. It extracts elevation at the grid points using bilinear interpolation. 
# TO DO: 
# Test to see if changing grid size makes a difference? Running now with 1000 because when I use 2000 there are some BGCs with only 1 point. 
coords <- makePointCoords(bgcs, elev, gridSize = 1000) |>
  Cache()

# Need to rename lon and lat to x and y to work with subsetByExtent(), and because they are in Albers, not lat/long:
setnames(coords, old = c("lon", "lat"), new = c("x", "y"))

# Extract BGC data from bgcs polygons and append as column to coords data:
points_sf <- st_as_sf(coords, coords = c("x", "y"), crs = 3005)
# points_sf <- st_transform(points_sf,3005) # Don't need because already 3005
bgc_att <- st_join(points_sf, bgcs)
bgc_att <- data.table(st_drop_geometry(bgc_att))

# bgc_att has 192010 unique IDs but 192012 rows. 
length(unique(bgc_att$id))
nrow(bgc_att)
bgc_att[duplicated(bgc_att$id), ] # ids 21321 and 34557 are duplicated in v13

bgc_duplicates <- bgc_att[bgc_att$id %in% c(21321, 34557), ]
bgc_duplicates <- merge(bgc_duplicates, coords, by = c("id", "elev"))
bgc_duplicates_sf <- st_as_sf(bgc_duplicates, coords = c("x", "y"), crs = 3005)

# TO DO: 
# Document why I'm removing these duplicates (likely overlapping polygons across single grid cell) 
checkarea <- ext(c(1593710  - 2500, 1593710 + 2500, 824076.6 - 2500, 824076.6 + 2500)) # xmin, xmax, ymin, ymax
bgcs_check <- st_crop(bgcs, checkarea)
plot(st_geometry(bgcs_check), col = as.factor(bgcs_check$BGC))
plot(st_geometry(bgc_duplicates_sf), col = "black", pch = 16, add = TRUE)
text(st_coordinates(bgc_duplicates_sf), labels = bgc_duplicates_sf$id, cex = 0.7, pos = 3, col = "black")

# Remove these duplicates:  
bgc_att <- unique(bgc_att, by = "id")

# Also remove rows where BGC is NA: 
nrow(bgc_att[is.na(bgc_att$BGC)]) # 294 rows with 1000 grid cells
bgc_att <- bgc_att[!is.na(bgc_att$BGC), ]

# Merge coords and BGC data from bgc_att: 
coords <- merge.data.table(coords, bgc_att, by = c("id", "elev"))

# Crop the coordinates to the study area defined above: 
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
# points(x = coords_gaps$x, y = coords_gaps$y, col = "black", cex = 0.001) # Just the gaps
# points(x = coords_trainWgaps$x, y = coords_trainWgaps$y, col = "white", cex = 0.001) # Everything but the gaps.

# Recombine coords_gaps and coords_trainingWgaps but with an extra column for "gap = yes or no" for holdouts. 
set(coords_gaps, j = "gap", value = "yes")
set(coords_trainWgaps, j = "gap", value = "no")

# coords_all <- rbind(coords_gaps, coords_trainWgaps)
coords_all <- rbindlist(list(coords_gaps, coords_trainWgaps))

#### Local feature selection: ####
# Copy code on local feature selection here

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
# TO DO: 
# Determine how sensitive the results are to this alpha value. 
trainData <- removeOutlier(as.data.frame(trainData), alpha = .025, vars = vars_all) |>
  Cache()

# TO DO: 
# Decide how and if to filter BGCs. Start with not filtering them, document reasons why to or not to do it. 

# How many points per BGC? Smallest number is 6 with gridsize = 1000, 10 < 30, 15 < 50. Largest: 10056. 
# (Smallest number is 1 with gridsize = 2000). 
BGCs_pre_downsample <- trainData %>%
  dplyr::group_by(BGC) %>%
  dplyr::summarize (n = n()) %>%
  dplyr::arrange(n)

# We do want to remove BAFA and "un" subzones (unvegetated and odd to predict): 
trainData <- trainData[!grepl("un|BAFA", trainData$BGC), ]

# Define bad BGCs and remove them: (These were selected in the RMarkdown script but I'm not sure why.) 
# badbgcs <- c("BWBSvk", "ICHmc1a", "MHun", "SBSun", "ESSFun", "SWBvk","MSdm3","ESSFdc3", "IDFdxx_WY", "MSabS", "FGff", "JPWmk_WY" )#, "ESSFab""CWHws2", "CWHwm", "CWHms1" , 
# trainData_bad <- trainData[BGC %in% badbgcs,]

# Remove very small sample BGC units (default cutoff = 30): 
# SKIP FOR NOW. Note - Courtney's most recent version removed BGCs with < 50.
# trainData <- rmLowSampleBGCs(trainData) |>
#   Cache()

# TO DO: Figure out if this is necessary/the best way to do it. It ensures that at most, larger BGCs have at most 90x as many rows as smallest BGC but it sets a ceiling for number of points. 

# # Subsample "oversampled" BGCs. This cuts them off at 90, if all BGCs are left in the sample, including those with only 1 point. 
# dataBalance_recipe <- recipe(BGC ~ ., data =  trainData) |>
#   step_downsample(BGC, under_ratio = 90) |>  ## subsamples "oversampled" BGCs
#   prep()

# Extract the data.table of balanced points:
# trainData_balanced <- dataBalance_recipe |>
#   juice() |>
#   as.data.table()

# See how many BGCs per with the downsample: 
# BGCs_post_downsample <- trainData_balanced %>% 
#   dplyr::group_by(BGC) %>% 
#   dplyr::summarize (n = n()) %>% 
#   dplyr::arrange(n)

#### Assess climate variability within BGCs: ####
# TO DO: 
# Figure out if more BGCs need to be removed due to high variability in combination with small sample sizes. 

#### Train ranger random forest model: ####
trainData <- as.data.table(trainData)
trainData[, BGC := as.factor(BGC)]
trainData_Wgaps <- trainData[gap == "no"]

cols_simple <- c("BGC", vars_simple)
cols_expert <- c("BGC", vars_expert)
cols_all <- c("BGC", vars_all)

# Train model with simple variables, on points from the entire study area: 
BGCmodel_full_simple <- ranger(
  BGC ~ .,
  data = trainData[, ..cols_simple],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
  
) |>
  Cache()

# Train model with simple variables, on points from area outside of the gaps: 
BGCmodel_Wgaps_simple <- ranger(
  BGC ~ .,
  data = trainData_Wgaps[, ..cols_simple],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, 
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
) |>
  Cache()

# Train model with expert selection of variables, on points from the entire study area: 
BGCmodel_full_expert <- ranger(
  BGC ~ .,
  data = trainData[, ..cols_expert],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, 
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
  
) |>
  Cache()

BGCmodel_Wgaps_expert <- ranger(
  BGC ~ .,
  data = trainData_Wgaps[, ..cols_expert],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, 
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
) |>
  Cache()

# Train model on all climate variables (kitchen sink scenario): 
BGCmodel_full_all <- ranger(
  BGC ~ .,
  data = trainData[, ..cols_all],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, 
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
) |>
  Cache()


BGCmodel_Wgaps_all <- ranger(
  BGC ~ .,
  data = trainData_Wgaps[, ..cols_all],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, 
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
) |>
  Cache()
beepr::beep()


# Check the models: 
conf_matrix_full_simple <- caret::confusionMatrix(data = predictions(BGCmodel_full_simple),
                                                  reference = trainData$BGC)
conf_matrix_Wgaps_simple <- caret::confusionMatrix(data = predictions(BGCmodel_Wgaps_simple),
                                                   reference = trainData_Wgaps$BGC)

conf_matrix_full_expert <- caret::confusionMatrix(data = predictions(BGCmodel_full_expert),
                                                  reference = trainData$BGC)
conf_matrix_Wgaps_expert <- caret::confusionMatrix(data = predictions(BGCmodel_Wgaps_expert),
                                                   reference = trainData_Wgaps$BGC)

conf_matrix_full_all <- caret::confusionMatrix(data = predictions(BGCmodel_full_all),
                                               reference = trainData$BGC)
conf_matrix_Wgaps_all <- caret::confusionMatrix(data = predictions(BGCmodel_Wgaps_all),
                                                reference = trainData_Wgaps$BGC)

print(BGCmodel_full_simple) # OOB prediction error: 21.87%
print(BGCmodel_Wgaps_simple) # OOB prediction error: 21.91%
print(BGCmodel_full_expert) # OOB prediction error: 22.99%
print(BGCmodel_Wgaps_expert) # OOB prediction error: 23.02%
print(BGCmodel_full_all) # OOB prediction error: 20.44%
print(BGCmodel_Wgaps_all) # OOB prediction error: 20.51%

#### Predictions: ####
# Rasterize the BGC polygons, assigning the values from the BGC column. 
bgcs_rast <- rasterize(bgcs, elev, field = "BGC")

# Align bgcs_rast and elev DEM to ensure resolution, extent, and CRS match: 
bgcs_rast <- resample(bgcs_rast, elev, method = "near")

# Merge bgcs_rast and elev into one multi-layer raster: 
bgcs_elev <- c(elev, bgcs_rast)
names(bgcs_elev) <- c("elev", "BGC")

# Reproject to lat/long to work with climr: 
bgcs_elev <- project(bgcs_elev, "EPSG:4326", method = "near")

# Make into data table: 
bgcs_elev_dt <- as.data.table(bgcs_elev, cells = TRUE, xy = TRUE) 
colnames(bgcs_elev_dt) <- c("id", "lon", "lat", "elev", "BGC")

# Get climr data for these raster coordinates:
clim_vars_preds2 <- downscale(
  xyz = bgcs_elev_dt,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars(),
  cache = TRUE)|>
  Cache()

clim_vars_preds <- copy(clim_vars_preds2)
ccissr::addVars(clim_vars_preds)

# Merge important info (lat, lon, elev, and BGC) back in: 
clim_vars_preds <- merge(clim_vars_preds, bgcs_elev_dt, by = "id") 

# Predictions based on full model and gap removed models: 
preds_full_simple <- predict(BGCmodel_full_simple, data = clim_vars_preds)
preds_Wgaps_simple <- predict(BGCmodel_Wgaps_simple, data = clim_vars_preds)
preds_full_expert <- predict(BGCmodel_full_expert, data = clim_vars_preds)
preds_Wgaps_expert <- predict(BGCmodel_Wgaps_expert, data = clim_vars_preds)
preds_full_all <- predict(BGCmodel_full_all, data = clim_vars_preds)
preds_Wgaps_all <- predict(BGCmodel_Wgaps_all, data = clim_vars_preds)

# Save into clim_vars_preds dataframe:
clim_vars_preds$preds_full_simple <- as.character(preds_full_simple$predictions)
clim_vars_preds$preds_Wgaps_simple <- as.character(preds_Wgaps_simple$predictions)
clim_vars_preds$preds_full_expert <- as.character(preds_full_expert$predictions)
clim_vars_preds$preds_Wgaps_expert <- as.character(preds_Wgaps_expert$predictions)
clim_vars_preds$preds_full_all <- as.character(preds_full_all$predictions)
clim_vars_preds$preds_Wgaps_all <- as.character(preds_Wgaps_all$predictions)

#### Check quality of predictions: ####
conf_matrix_preds_full_simple <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_full_simple, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

conf_matrix_preds_Wgaps_simple <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_Wgaps_simple, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

conf_matrix_preds_full_expert <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_full_expert, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

conf_matrix_preds_Wgaps_expert <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_Wgaps_expert, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

conf_matrix_preds_full_all <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_full_all, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

conf_matrix_preds_Wgaps_all <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_Wgaps_all, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

#### Look at feature importance: ####
# importance_scores_full_simple <- BGCmodel_full_simple$variable.importance
# importance_scores_full_simple<- sort(importance_scores_full_simple, decreasing = TRUE)
# 
# importance_scores_Wgaps_simple <- BGCmodel_Wgaps_simple$variable.importance
# importance_scores_Wgaps_simple<- sort(importance_scores_Wgaps_simple, decreasing = TRUE)
# 
# # Create a data frame for plotting
# importance_scores_full_simple_df <- data.frame(
#   Feature = names(importance_scores_full_simple),
#   Importance = importance_scores_full_simple
# )
# 
# importance_scores_Wgaps_simple_df <- data.frame(
#   Feature = names(importance_scores_Wgaps_simple),
#   Importance = importance_scores_Wgaps_simple
# )
# 
# # Plots
# ggplot(importance_scores_full_simple_df, aes(x = reorder(Feature, Importance), y = Importance)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   coord_flip() +
#   labs(title = "Feature Importance",
#        x = "Features",
#        y = "Importance") +
#   theme_minimal()
# 
# ggplot(importance_scores_Wgaps_simple_df, aes(x = reorder(Feature, Importance), y = Importance)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   coord_flip() +
#   labs(title = "Feature Importance",
#        x = "Features",
#        y = "Importance") +
#   theme_minimal()

#### Leaflet script: ####
# TO DO: 
# Try to simplify this code: 
# Merge subzone colors refs to match BGCs with BGC factors properly:
clim_vars_preds$preds_full_simple <- factor(clim_vars_preds$preds_full_simple, levels = subzones_colours_ref$BGC)
clim_vars_preds <- merge(clim_vars_preds, subzones_colours_ref, by.x = "preds_full_simple", by.y = "BGC") %>% 
  dplyr::rename(preds_full_simple_fct = BGC_num, RGB_full_simple = RGB)

clim_vars_preds$preds_Wgaps_simple <- factor(clim_vars_preds$preds_Wgaps_simple, levels = subzones_colours_ref$BGC)
clim_vars_preds <- merge(clim_vars_preds, subzones_colours_ref, by.x = "preds_Wgaps_simple", by.y = "BGC") %>% 
  dplyr::rename(preds_Wgaps_simple_fct = BGC_num, RGB_Wgaps_simple = RGB)

clim_vars_preds$preds_full_expert <- factor(clim_vars_preds$preds_full_expert, levels = subzones_colours_ref$BGC)
clim_vars_preds <- merge(clim_vars_preds, subzones_colours_ref, by.x = "preds_full_expert", by.y = "BGC") %>% 
  dplyr::rename(preds_full_expert_fct = BGC_num, RGB_full_expert = RGB)

clim_vars_preds$preds_Wgaps_expert <- factor(clim_vars_preds$preds_Wgaps_expert, levels = subzones_colours_ref$BGC)
clim_vars_preds <- merge(clim_vars_preds, subzones_colours_ref, by.x = "preds_Wgaps_expert", by.y = "BGC") %>% 
  dplyr::rename(preds_Wgaps_expert_fct = BGC_num, RGB_Wgaps_expert = RGB)

clim_vars_preds$preds_full_all <- factor(clim_vars_preds$preds_full_all, levels = subzones_colours_ref$BGC)
clim_vars_preds <- merge(clim_vars_preds, subzones_colours_ref, by.x = "preds_full_all", by.y = "BGC") %>% 
  dplyr::rename(preds_full_all_fct = BGC_num, RGB_full_all = RGB)

clim_vars_preds$preds_Wgaps_all <- factor(clim_vars_preds$preds_Wgaps_all, levels = subzones_colours_ref$BGC)
clim_vars_preds <- merge(clim_vars_preds, subzones_colours_ref, by.x = "preds_Wgaps_all", by.y = "BGC") %>% 
  dplyr::rename(preds_Wgaps_all_fct = BGC_num, RGB_Wgaps_all = RGB)

# Rasterize the predictions: 
template_raster <- bgcs_elev[[2]] 
n_cells <- ncell(template_raster)

preds_full_simple_rast <- rast(template_raster)
values(preds_full_simple_rast) <- NA    

preds_Wgaps_simple_rast <- rast(template_raster)
values(preds_Wgaps_simple_rast) <- NA         

preds_full_expert_rast <- rast(template_raster)
values(preds_full_expert_rast) <- NA    

preds_Wgaps_expert_rast <- rast(template_raster)
values(preds_Wgaps_expert_rast) <- NA         

preds_full_all_rast <- rast(template_raster)
values(preds_full_all_rast) <- NA    

preds_Wgaps_all_rast <- rast(template_raster)
values(preds_Wgaps_all_rast) <- NA         

# Map predictions to the template raster
# Assume 'id' matches the cell numbers in the raster
preds_full_simple_rast[clim_vars_preds$id] <- clim_vars_preds$preds_full_simple_fct
names(preds_full_simple_rast) <- "preds_full_simple_fct"

preds_Wgaps_simple_rast[clim_vars_preds$id] <- clim_vars_preds$preds_Wgaps_simple_fct
names(preds_Wgaps_simple_rast) <- "preds_Wgaps_simple_fct"

preds_full_expert_rast[clim_vars_preds$id] <- clim_vars_preds$preds_full_expert_fct
names(preds_full_expert_rast) <- "preds_full_expert_fct"

preds_Wgaps_expert_rast[clim_vars_preds$id] <- clim_vars_preds$preds_Wgaps_expert_fct
names(preds_Wgaps_expert_rast) <- "preds_Wgaps_expert_fct"

preds_full_all_rast[clim_vars_preds$id] <- clim_vars_preds$preds_full_all_fct
names(preds_full_all_rast) <- "preds_full_all_fct"

preds_Wgaps_all_rast[clim_vars_preds$id] <- clim_vars_preds$preds_Wgaps_all_fct
names(preds_Wgaps_all_rast) <- "preds_Wgaps_all_fct"

# Merge with bgcs data, reproject to lat/long as required by leaflet: 
bgcs2 <- merge(bgcs, subzones_colours_ref, by = "BGC")
bgcs2 <- st_transform(bgcs2, crs = 4326) # This is really fast... why can't I do this from the beginning? 

# Add gaps: 
# First, reproject to lat/long: 
gap_poly2 <- project(gap_poly, "EPSG:4326")

# Create a color factor mapping BGC to RGB colors
color_pal <- colorNumeric(
  palette = clim_vars_preds$RGB_full_simple,
  domain = clim_vars_preds$preds_full_simple_fct
)

# # Rename values in the predictions raster to BGC_num so that they match up with the reference: 
# names(preds_full_simple_rast) <- "BGC_num"

# Leaflet:
leaflet() %>%
  addTiles()  %>%
  addPolygons(
    data = gap_poly2,
    fillColor = "transparent",
    color = "black",
    weight = 2,
    opacity = 1,
    fillOpacity = 0,
    popup = ~paste("Gap Polygon"),
    group = "Gap Extents"
  ) %>%
  addPolygons(
    data = bgcs2,
    fillColor = ~RGB,
    color = ~RGB,
    weight = 1,
    opacity = 1,
    fillOpacity = 1,
    popup = ~paste("Zone:", BGC),
    group = "BGC Zones"
  ) %>%
  addRasterImage(
    preds_full_simple_rast,
    colors = color_pal,
    opacity = 1,
    group = "Preds: full, simple"
  ) %>%
  addRasterImage(
    preds_Wgaps_simple_rast,
    colors = color_pal,
    opacity = 1,
    group = "Preds: Wgaps, simple"
  ) %>%
  addRasterImage(
    preds_full_expert_rast,
    colors = color_pal,
    opacity = 1,
    group = "Preds: full, expert"
  ) %>%
  addRasterImage(
    preds_Wgaps_expert_rast,
    colors = color_pal,
    opacity = 1,
    group = "Preds: Wgaps, expert"
  ) %>%
  addRasterImage(
    preds_full_all_rast,
    colors = color_pal,
    opacity = 1,
    group = "Preds: full, all"
  ) %>%
  addRasterImage(
    preds_Wgaps_all_rast,
    colors = color_pal,
    opacity = 1,
    group = "Preds: Wgaps, all"
  ) %>%
  addLayersControl(
    overlayGroups = c("Gap Extents", 
                      "BGC Zones", 
                      "Preds: full, simple", 
                      "Preds: Wgaps, simple",
                      "Preds: full, expert",
                      "Preds: Wgaps, expert",
                      "Preds: full, all",
                      "Preds: Wgaps, all"),
    options = layersControlOptions(collapsed = FALSE)
  )

beepr::beep()
