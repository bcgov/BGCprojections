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

# Source some functions: 
source("R/utils.R")

# Set default cache directories for the reproducible and climr packages (where intermediate results/downloaded data will be stored): 
options(reproducible.cachePath = "reproducible.cache/",
        climr.cache.path = "climr.cache/")

#### Create training points: ####
# Load in BGC polygons: 
bgcs <- st_read("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/BGCProjections/Proxy_ObjectStorage/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")
# bgcs <- st_read("//objectstore2.nrs.bcgov/ffec/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")


# And the DEM:
elev <- rast("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/BGCProjections/Proxy_ObjectStorage/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")
# elev <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")

# Reproject elev to Albers: 
elev <- project(elev, crs(bgcs))

# Define smaller study area for script building purposes. These are the extents in lat/long. Want them in Albers instead.

# Remove this later and run script on training area instead. 
# trainingarea <- ext(c(-125, -112, 43, 55))
studyarea <- ext(c(-123, -117, 49, 52.5))

# Create a SpatRaster to represent the extents in lat/long
# dummy_raster <- rast(ext = trainingarea, crs = "EPSG:4326", res = 0.1)
dummy_raster <- rast(ext = studyarea, crs = "EPSG:4326", res = 0.1)  

# Reproject the dummy raster to Albers (EPSG:3005)
dummy_raster_albers <- project(dummy_raster, "EPSG:3005")

# Extract the reprojected extents
studyarea_albers <- ext(dummy_raster_albers)

# Crop the elev DEM and bgcs to just the smaller study area for faster execution: 
elev <- crop(elev, studyarea_albers)
bgcs <- st_crop(bgcs, studyarea_albers)

# From climr documentation: "Since climr is meant to be used to downscale climate variables in land, we will “clip” (set values outside the polygon to NAs) the raster using a land-only polygonRemove areas with water: * Confirm that I need to do this!* 
# elev <- mask(elev, bgcs)

# This function makes a grid over the extent of bgcs, and fills it with a dummy variable (1L), converts to a spatial vector. It extracts elevation at the grid points using bilinear interpolation. 
# Using gridSize = 0.018 for now because that's roughly 2 km latitude (but only 1.18 km longitude): 
coords <- makePointCoords(bgcs, elev, gridSize = 2000) |>
  Cache()

# Need to rename lon and lat to x and y to work with subsetByExtent(), and because they are in Albers, not lat/long: 
setnames(coords, old = c("lon", "lat"), new = c("x", "y"))

# Extract BGC data from bgcs polygons and append as column to coords data:
points_sf <- st_as_sf(coords, coords = c("x", "y"), crs = 3005)
# points_sf <- st_transform(points_sf,3005) # Don't need because already 3005
bgc_att <- st_join(points_sf, bgcs)
bgc_att <- data.table(st_drop_geometry(bgc_att))

# bgc_att has 48108 unique IDs but 48108 rows. 
length(unique(bgc_att$id))
nrow(bgc_att)
bgc_att[duplicated(bgc_att$id), ] # 4979, 4980, 10485, and 46026 are duplicated. 

# REMOVE DUPLICATES FOR NOW: 
bgc_att[bgc_att$id == 4979, ]
bgc_att[bgc_att$id == 4980, ]
bgc_att[bgc_att$id == 10485, ]
bgc_att[bgc_att$id == 46026, ]

# Remove duplicates for now: 


# Summarize how many points in each  BGC, for now, retain only those where N > 10 for now:  
BGC_counts <- bgc_att[, .(Num = .N), by = .(BGC)] 
BGC_counts <- BGC_counts[Num >= 10]
bgc_att_filtered <- bgc_att[BGC %in% BGC_counts$BGC]

# Merge coords and BGC data from bgc_att: 
coords2 <- merge(coords, bgc_att_filtered, by = c("id", "elev"))

# This crops the coordinates to the study area defined above: 
# Note: coords_train will be identical to coords for now since I already cropped to the size of the smaller study area earlier but when I rerun with entire training area, it will be different. 
coords_train <- subsetByExtent(coords2, studyarea_albers)

# Make rectangular gap extents within the bounding box of the study area. 5L is the default number of gaps to create. 
gapextents <- makeGapExtents(studyarea_albers, 5L)

# Converts list of spatial extents into to polygons: 
gap_poly <- lapply(gapextents, vect, crs = "EPSG:3005")

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

points(x = coords_train$x, y = coords_train$y, col = "grey50", cex = 0.001) # All points
points(x = coords_gaps$x, y = coords_gaps$y, col = "black", cex = 0.001) # Just the gaps
points(x = coords_trainWgaps$x, y = coords_trainWgaps$y, col = "white", cex = 0.001) # Everything but the gaps.

# Recombine coords_gaps and coords_trainingWgaps but with an extra column for "Gap = Yes or No". 
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
# First, just simply PPT, Tmax, Tmin: 
vars_simple <- c("PPT", "Tmax", "Tmin")

# Define a more complex set (ecologically relevant?): 
vars_more <- c("DD5", "DD_0_at", "DD_0_wt", "PPT05", "PPT06", "PPT07", "PPT08",
               "PPT09", "CMD", "PPT_at", "PPT_wt", "CMD07", "SHM", "AHM", "NFFD", "PAS", "CMI")

# All for pairwise variable selection? *Come back to this.*

# Pull from climr (just gaps for now? Check w Colin): 
clim_vars <- downscale(
  xyz = coords_all,
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
# coords_train <- coords_train[clim_vars[, .(id)], on = "id", nomatch = 0L]

coords_gaps <- coords_gaps[clim_vars[, .(id)], on = "id", nomatch = 0L] 
# Remove duplicates:
coords_gaps  <- coords_gaps %>% 
  distinct()

# coords_trainWgaps <- coords_trainWgaps[clim_vars[, .(id)], on = "id", nomatch = 0L]

# Assess climate variability within BGCs:
# First, add long, lat, and elevation back in: 
# clim_vars_all <- left_join(clim_vars, coords_train, relationship = "many-to-many") %>% 
#  distinct()

# clim_vars_Wgaps <- left_join(clim_vars, coords_trainWgaps, relationship = "many-to-many") %>% 
#  distinct()

clim_vars_gaps <- left_join(clim_vars, coords_gaps, relationship = "many-to-many") %>% 
  distinct() 

# Figure out which BGC each point is in. 
# Turn data.table object into a SpatVector:
# coords_train_vect <- vect(coords_train, geom = c("lon", "lat"), crs = crs(bgcs))
# coords_trainWgaps_vect <- vect(coords_trainWgaps, geom = c("lon", "lat"), crs = crs(bgcs))
coords_gaps_vect <- vect(coords_gaps, geom = c("lon", "lat"), crs = crs(bgcs))

# Extract values (BGC) at point locations: 
s <- sample(1:dim(coords_gaps)[1], 20)
system.time({
  test <- terra::extract(bgcs, coords_gaps[s, c(2,3)])
})



# Once I have the BGCs for each point, I can see which BGCs are bad and filter those out: 
# BGC_counts <- clim_vars_gaps[, .(Num = .N), by = .(BGC)]   ## (Not sure what the criteria used here is)

# Define bad BGCs and remove them: 
badbgcs <- c("BWBSvk", "ICHmc1a", "MHun", "SBSun", "ESSFun", "SWBvk","MSdm3","ESSFdc3", "IDFdxx_WY", "MSabS", "FGff", "JPWmk_WY" )#, "ESSFab""CWHws2", "CWHwm", "CWHms1" , 
trainData <- trainData[!BGC %in% badbgcs,]

## set alpha for removal of outliers (2.5% = 3SD)
trainData <- removeOutlier(as.data.frame(trainData), alpha = .025, vars = vs_final) |>
  Cache()

# Remove very small sample units: 
trainData <- rmLowSampleBGCs(trainData) |>
  Cache()

# Subsample "oversampled" BGCs: 
dataBalance_recipe <- recipe(BGC ~ ., data =  trainData) |>
  step_downsample(BGC, under_ratio = 90) |>  ## subsamples "oversampled" BGCs
  prep()

## extract data.table
trainData_balanced <- dataBalance_recipe |>
  juice() |>
  as.data.table()

# Train ranger random forest model: 

# BGC_Nums <- trainData_balanced[,.(Num = .N), by = BGC]   ## for inspection