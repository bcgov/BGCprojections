library(leaflet)

# Define color palettes for trainData_balanced_gaps2# Define color palettes for the models
subzones_colours_ref <- fread("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNAv13_Subzone_colours_2.csv") %>% 
  dplyr::select(!c(fid, MAP_LABEL, NSRNAME, ZONE)) %>% 
  dplyr::mutate(BGC_num = as.numeric(as.factor(BGC)))

# Copy DEM: 
preds_full_simple_DEM <- elev

# Initialize the raster values with NA
values(preds_full_simple_DEM) <- NA

# Assign predictions to raster cells by ID
preds_full_simple_DEM[clim_vars_preds$id] <- clim_vars_preds$preds_full_simple_num

# Also merge with bgcs data for comparison. Reproject to lat/long as required by leaflet: 
bgcs2 <- merge(bgcs, subzones_colours_ref, by = "BGC")
bgcs2 <- st_transform(bgcs2, crs = 4326) 

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

plot(preds_full_simple_DEM)

preds_full_simple_DEM <- project(preds_full_simple_DEM, "EPSG:4326")

# Leaflet:
leaflet(clim_vars_preds) %>%
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
    preds_full_simple_DEM, 
    colors = color_pal, 
    opacity = 1,    
    group = "Preds: full, simple"
  ) %>% 
  addLayersControl(
    overlayGroups = c("Gap Extents", "BGC Zones", "Preds: full, simple"),  
    options = layersControlOptions(collapsed = FALSE)  
  )

beepr::beep()
