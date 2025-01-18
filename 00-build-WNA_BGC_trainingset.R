# Random Forest Model of Biogeoclimatic Units for Western North America
# Original script: Build_WNA_BGC_trainingset.Rmd by William H MacKenzie & Kiri Daust

# Updated by Deb Obrist (January 2025)

# Load packages: 
library(tidyverse)
library(terra)
library(climr)
library(reproducible) # For Cache function

# Source some functions: 
source("R/utils.R")

#### Get training points: ####

# Create new file paths for Deb's temporary data location: 
bgcs <- vect("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/BGCProjections/Proxy_ObjectStorage/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")

elev <- rast("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/BGCProjections/Proxy_ObjectStorage/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")


# Reproject both to lat/long: 
bgcs <- project(bgcs, elev)

# Define smaller study area for script building purposes: 
# Remove this later and run script on training area instead. 
trainingarea <- ext(c(-125, -112, 43, 55))
studyarea <- ext(c(-123, -117, 49, 52.5))

# Crop the elev DEM and bgcs to just the smaller study area for faster execution: 
elev <- crop(elev, studyarea)
bgcs <- crop(bgcs, studyarea)

# From climr documentation: "Since climr is meant to be used to downscale climate variables in land, we will “clip” (set values outside the polygon to NAs) the raster using a land-only polygonRemove areas with water: * Confirm that I need to do this!* 
# elev <- mask(elev, bgcs)

# This function makes a grid over the extent of bgcs, and fills it with a dummy variable (1L), converts to a spatial vector. It extracts elevation at the grid points using bilinear interpolation. 
# Using gridSize = 0.018 for now because that's roughly 2 km latitude (but only 1.18 km longitude): 
coords <- makePointCoords(bgcs, elev, gridSize = 0.018) |>
  Cache()

# Need to rename lon and lat to x and y to work with subsetByExtent(): 
setnames(coords, old = c("lon", "lat"), new = c("x", "y"))

# This crops the coordinates to the study area defined above: 
# Note: coords_train will be identical to coords for now since I already cropped to the size of the smaller study area earlier but when I rerun with entire training area, it will be different. 
coords_train <- subsetByExtent(coords, studyarea)

# Make rectangular gap extents within the bounding box of the study area. 5L is the default number of gaps to create. 
gapextents <- makeGapExtents(studyarea, 5L)

# Converts list of spatial extents into to polygons: 
gap_poly <- lapply(gapextents, vect, crs = "EPSG:4326")

# Combines all individual polygons into one spatial object. 
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
points(x = coords_train$x, y = coords_train$y, col = "black", cex = 0.001)
points(x = coords_gaps$x, y = coords_gaps$y, col = "grey50", cex = 0.001)

#### Get climate variables: ####
# Define variables needed: 
# First, just simply PPT, Tmax, Tmin: 
vars_simple <- c("PPT", "Tmax", "Tmin")

# Define a more complex set (ecologically relevant?): 
vars_more <- c("DD5", "DD_0_at", "DD_0_wt", "PPT05", "PPT06", "PPT07", "PPT08",
               "PPT09", "CMD", "PPT_at", "PPT_wt", "CMD07", "SHM", "AHM", "NFFD", "PAS", "CMI")

# All for pairwise variable selection? *Come back to this.*

# coords_train must have the following column names for climr: id, lon, lat, elev: 
coords_train <- coords_train %>% 
  rename(lon = x, lat = y) %>% 
  select(id, lon, lat, elev)

# Pull from climr: 
clim_vars <- downscale(
  xyz = coords_train,
  which_refmap = "refmap_climr", 
  obs_periods = "2001_2020", # Courtney's code has this. 
  gcm_periods = "2021_2040", # I suppose I should do all actually? Come back here.  
  # gcms = list_gcms()[c(1, 4, 5, 6, 7, 10, 11, 12)], # 8 GCMs recommended in Mahony et al. 2022
  gcms = "CanESM5", # Just picking one for now to make the left join run faster: 
  ssps = "ssp245",
  max_run = 2, # How many should I run? 
  return_refperiod = TRUE, # Also return the 1961-1990 normals period. 
  vars = vars_simple,
  cache = TRUE)|>
  Cache()

# Subset coords_train and coords_trainWgaps to include only rows where the id column matches an id in clim_vars:
coords_train <- coords_train[clim_vars[, .(id)], on = "id", nomatch = 0L]
coords_trainWgaps <- coords_trainWgaps[clim_vars[, .(id)], on = "id", nomatch = 0L]

# Assess climate variability within BGCs:
# First, add long, lat, and elevation back in: 
clim_vars <- left_join(clim_vars, coords_train, relationship = "many-to-many") %>% 
  distinct()

# Figure out which BGC each point is in. 
# Turn data.table object into a SpatVector:
coords_train_vect <- vect(coords_train, geom = c("lon", "lat"), crs = crs(bgcs))
points(coords_train_vect, col = "red")
# Extract values (BGC) at point locations: 
# S4 method for class 'SpatVector,SpatVector'
# extract(x, y)

test <- extract(bgcs, coords_train_vect)

