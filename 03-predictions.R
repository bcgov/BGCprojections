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


# To do: 
# Read in models: 
BGCmodel_full_simple <- readRDS("RF_models/BGCmodel_full_simple_V1.rds")
BGCmodel_Wgaps_simple <- readRDS("RF_models/BGCmodel_Wgaps_simple_V1.rds")
BGCmodel_full_expert <- readRDS("RF_models/BGCmodel_full_simple_V1.rds")
BGCmodel_Wgaps_expert <- readRDS("RF_models/BGCmodel_Wgaps_expert_V1.rds")
BGCmodel_full_all <- readRDS("RF_models/BGCmodel_full_all_V1.rds")
BGCmodel_Wgaps_all <- readRDS("RF_models/BGCmodel_Wgaps_all_V1.rds")

# Read in and prepare DEM and BGCs again: 
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

#### Look at how they vary by size of BGC, number of points: ####
# This is where I was trying to figure out if there is an obvious "minimum" number of sample points required to make "good" predictions, which I was evaluating based on sensitivity. Not sure if this is the best metric. Also not sure what a "good" sensitivity is. 

# How many points per BGC? 
BGCs_pre_downsample <- trainData %>%
  dplyr::group_by(BGC) %>%
  dplyr::summarize (n = n()) %>%
  dplyr::arrange(n)

# Calculate area of each BGC: 
bgcs_areas <- st_area(bgcs)
bgcs_dt <- as.data.table(bgcs)
bgcs_dt[, area_sq_km := as.numeric(bgcs_areas) / 1e6]
bgcs_dt <- bgcs_dt[, .(BGC, area_sq_km)]

bgcs_dt <- merge(bgcs_dt, BGCs_pre_downsample, by = "BGC")

ggplot(bgcs_dt, aes(x = area_sq_km, y = n)) +
  geom_point() +
  theme_bw() +
  theme(panel.grid = element_blank())

sensitivity_full_simple <- conf_matrix_preds_full_simple$byClass[, "Sensitivity"]
sensitivity_Wgaps_simple <- conf_matrix_preds_Wgaps_simple$byClass[, "Sensitivity"]
sensitivity_full_expert <- conf_matrix_preds_full_expert$byClass[, "Sensitivity"]
sensitivity_Wgaps_expert <- conf_matrix_preds_Wgaps_expert$byClass[, "Sensitivity"]
sensitivity_full_all <- conf_matrix_preds_full_all$byClass[, "Sensitivity"]
sensitivity_Wgaps_all <- conf_matrix_preds_Wgaps_all$byClass[, "Sensitivity"]

sensitivity_dt <- data.table(
  BGC = rownames(conf_matrix_preds_full_simple$byClass),
  Sensitivity_full_simple = sensitivity_full_simple,
  Sensitivity_Wgaps_simple = sensitivity_Wgaps_simple,
  Sensitivity_full_expert = sensitivity_full_expert,
  Sensitivity_Wgaps_expert = sensitivity_Wgaps_expert,
  Sensitivity_full_all = sensitivity_full_all,
  Sensitivity_Wgaps_all = sensitivity_Wgaps_all
)

sensitivity_dt[, BGC := gsub("^Class: ", "", BGC)]

check_preds <- merge(bgcs_dt, sensitivity_dt, by = "BGC")

ggplot(check_preds, aes(x = log10(n), y = Sensitivity_full_simple)) +
  geom_point() +
  geom_smooth()

setdiff(sensitivity_dt$BGC, bgcs_dt$BGC) # BAFAun, CMAun, CMAun_WA, IMAun, IMAun_WA are not in here.

#### Look at feature importance: ####
# This is one of the things I accidentally deleted when I had my GitHub mishap earlier today and now this is a FAR earlier version. I don't know if it's relevant right now anyways! 


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


#### Export csvs: ####
write.csv(clim_vars_preds, "data-generated/clim_vars_preds.csv", row.names = FALSE)