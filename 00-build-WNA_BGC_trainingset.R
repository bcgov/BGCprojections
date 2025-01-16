# Random Forest Model of Biogeoclimatic Units for Western North America
# Original script: Build_WNA_BGC_trainingset.Rmd by William H MacKenzie & Kiri Daust

# Load packages: 
library(terra)

# Source some functions: 
source("R/utils.R")


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

# Visualize: 
plot(elev)

for(i in 1:5){
  plot(gapextents[[i]], add=T)
}

# Converts list of spatial extents into to polygons: 
gap_poly <- lapply(gapextents, vect, crs = "EPSG:4326")

# Combines all individual polygons into one spatial object. 
gap_poly <- do.call(rbind, gap_poly)

# Filters points in coords that fall within the gap polygons.
coords_gaps <- subsetByExtent(coords_train, gap_poly)

# Removes the points in coords_gaps from coords_train to produce a dataset of points that do not fall within the gaps. It does this by keeping only points in coords_train that do not match the id values in coords_gaps. 
coords_trainWgaps <- coords_train[!coords_gaps, on = "id"]
