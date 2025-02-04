library(leaflet)

# Define color palettes for trainData_balanced_gaps2# Define color palettes for the models
subzones_colours_ref <- fread("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNAv13_Subzone_colours_2.csv") %>% 
  dplyr::select(!c(fid, MAP_LABEL, NSRNAME, ZONE)) %>% 
  dplyr::mutate(BGC_num = as.numeric(as.factor(BGC)))

# Reproject elev raster to WGS84 (EPSG:4326)
preds_full_simple_DEM <- project(elev, "EPSG:4326")
# elev_wgs <- project(elev, "EPSG:4326")

# Convert the factor levels to numeric values
clim_vars_preds$preds_full_simple_num <- as.numeric(clim_vars_preds$preds_full_simple)

# Initialize the raster values with NA
values(preds_full_simple_DEM) <- NA

# Assign predictions to the corresponding raster cells by valid ID
preds_full_simple_DEM[clim_vars_preds$id] <- clim_vars_preds$preds_full_simple_num

# Also merge with bgcs data for comparison. Reproject to lat/long as required by leaflet: 
bgcs2 <- merge(bgcs, subzones_colours_ref, by = "BGC")
bgcs2 <- st_transform(bgcs2, crs = 4326) 


plot(preds_full_simple_DEM)

# Add gaps: 
# First, reproject to lat/long: 
gap_poly2 <- project(gap_poly, "EPSG:4326")

# Create a color factor mapping BGC to RGB colors
color_pal <- colorNumeric(
  palette = subzones_colours_ref$RGB,  
  domain = values(preds_full_simple_DEM)   
)

color_pal_elev <- colorNumeric(
  palette = viridis::viridis(256), 
  domain = values(elev_wgs), 
  na.color = "transparent"  
)

?colorNumeric
# Check if the assignment worked (optional)
head(values(preds_full_simple_DEM))

# Leaflet visualization
leaflet() %>%
  addTiles() %>%
  addRasterImage(
    preds_full_simple_DEM, 
    colors = color_pal, 
 #   opacity = 0.5, 
    group = "Preds: full, simple"
  ) %>%
  addLayersControl(
    overlayGroups = c("Preds: full, simple"),
    options = layersControlOptions(collapsed = FALSE)
  )

# Leaflet:
leaflet(clim_vars_preds) %>%
  addTiles() %>%
  addPolygons(
    data = bgcs2,            
    fillColor = ~RGB,        
    color = ~RGB,            
    weight = 1,              
    opacity = 1,             
    fillOpacity = 0.5,       
    popup = ~paste("Zone:", BGC),  
    group = "BGC Zones"      
  ) %>%
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
  addRasterImage(
    preds_full_simple_DEM, 
    # colors = color_pal, 
    opacity = 0.5,    
    group = "Preds: full, simple"
  ) %>%
  # addRasterImage(
  #   elev_wgs,  # Raster object (must be in EPSG:4326)
  #   colors = color_pal_elev,  # Apply the color palette
  #   opacity = 0.7,  # Set opacity for the raster
  #   group = "Elevation"
  # ) %>% 
  addLayersControl(
    overlayGroups = c("BGC Zones", "Gap Extents", "Preds: full, simple", "Elevation"),  
    options = layersControlOptions(collapsed = FALSE)  
  )

beepr::beep()
