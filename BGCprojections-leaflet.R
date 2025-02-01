library(leaflet)

# Define color palettes for trainData_balanced_gaps2# Define color palettes for the models
subzones_colours_ref <- fread("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNAv13_Subzone_colours_2.csv") %>% 
  dplyr::select(!c(fid, MAP_LABEL, NSRNAME, ZONE)) 

# Reproject elev raster to WGS84 (EPSG:4326)
elev <- project(elev, "EPSG:4326")

# Assign predictions to raster cells
elev[clim_vars_preds[, id]] <- clim_vars_preds[, preds_full_simple]

# Also merge with bgcs data for comparison. Reproject to lat/long as required by leaflet: 
bgcs2 <- merge(bgcs, subzones_colours_ref, by = "BGC")
bgcs2 <- st_transform(bgcs2, crs = 4326) 

# Add gaps: 
# First, reproject to lat/long: 
gap_poly2 <- project(gap_poly, "EPSG:4326")

# Create a color factor mapping BGC to RGB colors
color_pal <- colorFactor(
  palette = subzones_colours_ref$RGB,  # The RGB values from your table
  domain = subzones_colours_ref$BGC   # The BGC categories from your table
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
    elev, # This is the reprojected raster
    colors = color_pal,  # Apply the color palette
    opacity = 0.5,     # Set the opacity of the raster layer
    group = "Preds: full, simple"
  ) %>%
  addLayersControl(
    overlayGroups = c("BGC Zones", "Gap Extents", "Preds: full, simple"),  
    options = layersControlOptions(collapsed = FALSE)  
  )

beepr::beep()

mismatched_predictions <- trainData_balanced_gaps[predictions_full != BGC]
mismatched_predictions2 <- trainData_balanced_gaps[predictions_Wgaps != BGC]
mismatched_predictions3 <- trainData_balanced_gaps[predictions_Wgaps != predictions_full]

nrow(mismatched_predictions)
nrow(mismatched_predictions2)
nrow(mismatched_predictions3)

# # Match up colors and zones: 
# trainData_balanced_gaps[, predictions_full := as.character(predictions_full)]
# subzones_colours_ref[, BGC := as.character(BGC)]
# 
# trainData_balanced_gaps <- trainData_balanced_gaps[subzones_colours_ref, 
#                                                    on = .(predictions_full = BGC), 
#                                                    nomatch = 0][
#                                                      , color_full := RGB
#                                                    ][, RGB := NULL]
# 
# trainData_balanced_gaps <- trainData_balanced_gaps[subzones_colours_ref, 
#                                                    on = .(predictions_Wgaps = BGC), 
#                                                    nomatch = 0][
#                                                      , color_Wgaps := RGB
#                                                    ][, RGB := NULL]

# Ensure preds_full_simple is a factor with correct levels
clim_vars_preds$preds_full_simple <- factor(clim_vars_preds$preds_full_simple, levels = subzones_colours_ref$BGC)

leaflet() %>%
  addTiles() %>%
  addPolygons(
    data = bgcs2,            
    fillColor = ~color_pal(BGC),
    color = ~color_pal(BGC),
    weight = 1,              
    opacity = 1,             
    fillOpacity = 0.5,       
    popup = ~paste("Zone:", BGC),
    group = "BGC Zones"
  ) %>%
  addRasterImage(
    elev,  
    colors = color_pal,  
    opacity = 0.5,     
    group = "Preds: full, simple"
  ) %>%
  addLayersControl(
    overlayGroups = c("BGC Zones", "Gap Extents", "Preds: full, simple"),
    options = layersControlOptions(collapsed = FALSE)
  )
