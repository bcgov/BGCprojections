library(leaflet)

# Define color palettes for trainData_balanced_gaps2# Define color palettes for the models
subzones_colours_ref <- fread("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNAv13_Subzone_colours_2.csv") %>% 
  dplyr::select(!c(fid, MAP_LABEL, NSRNAME, ZONE)) 

# Create template raster: 
gap_DEM <- crop(elev, ext(gap_poly))
gap_DEM <- mask(gap_DEM, gap_poly)
plot(gap_DEM)

# Convert predictions_full to a factor
trainData_balanced_gaps[, predictions_full := as.factor(predictions_full)]

# Assign predictions to raster cells
gap_DEM[trainData_balanced_gaps[, id]] <- trainData_balanced_gaps[, predictions_full]

# Plot the resulting raster
plot(gap_DEM, main = "Predictions Full")



# Convert trainData_balanced_gaps coordinates into a spatial object
train_coords <- trainData_balanced_gaps[, .(lon, lat)]  # Extract lon/lat
train_sp <- vect(train_coords, crs = crs(gap_DEM))  # Make it a spatial object with the same CRS as the raster

# Get the corresponding cell indices in the raster based on the coordinates
cell_indices <- cellFromXY(gap_DEM, train_sp)

# Assign the values to the raster based on those indices
values(gap_DEM)[cell_indices] <- as.character(trainData_balanced_gaps$predictions_full)

# Check the result
summary(values(gap_DEM))






# Match up colors and zones: 
trainData_balanced_gaps[, predictions_full := as.character(predictions_full)]
subzones_colours_ref[, BGC := as.character(BGC)]

trainData_balanced_gaps <- trainData_balanced_gaps[subzones_colours_ref, 
                                                   on = .(predictions_full = BGC), 
                                                   nomatch = 0][
                                                     , color_full := RGB
                                                   ][, RGB := NULL]

trainData_balanced_gaps <- trainData_balanced_gaps[subzones_colours_ref, 
                                                   on = .(predictions_Wgaps = BGC), 
                                                   nomatch = 0][
                                                     , color_Wgaps := RGB
                                                   ][, RGB := NULL]

# Also merge with bgcs data for comparison. Reproject to lat/long as required by leaflet: 
bgcs2 <- merge(bgcs, subzones_colours_ref, by = "BGC")
bgcs2 <- st_transform(bgcs2, crs = 4326) 

# Add gaps: 
# First, reproject to lat/long: 
gap_poly2 <- project(gap_poly, "EPSG:4326")

# Leaflet:
leaflet(trainData_balanced_gaps) %>%
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
  addCircleMarkers(
    lat = ~lat, lng = ~lon,  
    color = ~color_full,           
    radius = 5, 
    fillOpacity = 0.7,
    popup = ~paste("Prediction:", predictions_full),
    group = "Predictions - All"  
  ) %>%
  addCircleMarkers(
    lat = ~lat, lng = ~lon,  
    color = ~color_Wgaps,            
    radius = 5,
    fillOpacity = 0.7,
    popup = ~paste("Prediction with Gaps:", predictions_Wgaps),
    group = "Predictions with Gaps"  
  ) %>% 
  addLayersControl(
    overlayGroups = c("Predictions - All", "Predictions with Gaps", "BGC Zones", "Gap Extents"),  
    options = layersControlOptions(collapsed = FALSE)  
  )

beepr::beep()

mismatched_predictions <- trainData_balanced_gaps[predictions_full != BGC]
mismatched_predictions2 <- trainData_balanced_gaps[predictions_Wgaps != BGC]
mismatched_predictions3 <- trainData_balanced_gaps[predictions_Wgaps != predictions_full]

nrow(mismatched_predictions)
nrow(mismatched_predictions2)
nrow(mismatched_predictions3)