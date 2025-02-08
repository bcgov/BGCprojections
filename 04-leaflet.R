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
library(leaflet)

# Read in predictions data: 
# NOTE - this is quite large so I'm not going to push to Github but it can be made locally by running the 03 script. 
clim_vars_preds <- read.csv("data-generated/clim_vars_preds.csv")

# Maybe make a separate 00 script with leaflet set up so that I can just run leaflet here. 
# Also try to simplify this code. The set up is in part copied and pasted from other scripts, especially the part setting up the raster template because I wasn't sure if it was a good idea to save copies of those large spatial objects...

#### Set up: ####
# Read in colours for reference of factors: 
subzones_colours_ref <- fread("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNAv13_Subzone_colours_2.csv") %>% 
  dplyr::select(!c(fid, MAP_LABEL, NSRNAME, ZONE)) %>% 
  dplyr::mutate(BGC_num = as.numeric(as.factor(BGC)))


clim_vars_preds <- read.csv("data-generated/clim_vars_preds.csv")

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

#### Leaflet: ####
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
