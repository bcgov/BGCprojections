library(leaflet)

# Combine data into a single data frame for plotting
trainData_balanced_TEST <- trainData_balanced

# Add in the predictions for the model version with no gaps/holdouts. 
trainData_balanced_TEST$predictions_full <- predictions_full

# Deal with gaps/holdouts: 
trainData_balanced_TEST$predictions_Wgaps <- NA
trainData_balanced_TEST$predictions_Wgaps[trainData_balanced_TEST$gap == 1] <- as.character(predictions_Wgaps)

mismatched_predictions <- trainData_balanced_TEST[predictions_Wgaps != BGC]


# Define color palettes for the models
subzones_colours_ref <- fread("C:/Users/dobrist/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNAv13_Subzone_colours_2.csv") %>% 
  dplyr::select(!c(fid, MAP_LABEL, NSRNAME, ZONE)) 

# Match up colors and zones: 
trainData_balanced_TEST <- merge(
  trainData_balanced_TEST,
  subzones_colours_ref,  
  by.x = "predictions_full",  
  by.y = "BGC",
  all.x = TRUE
)

# Rename the color column for clarity
setnames(trainData_balanced_TEST, "RGB", "color_full")

# Add colors for predictions_Wgaps
trainData_balanced_TEST <- merge(
  trainData_balanced_TEST,
  subzones_colours_ref,
  by.x = "predictions_Wgaps", # Match predictions_Wgaps with the BGC in bgc_colors
  by.y = "BGC",
  all.x = TRUE
)

# Rename the second color column
setnames(trainData_balanced_TEST, "RGB", "color_Wgaps")


# Some colors missing: 
setdiff(trainData_balanced_TEST$BGC, trainData_balanced_TEST2$BGC) # "CWHdm"  "CWHxm1" "IMAunp"
nrow(trainData_balanced_TEST[trainData_balanced$BGC %in% c("CWHdm", "CWHxm1", "IMAunp")])


# There is not a colour assigned to "CWHdm"  "CWHxm1" "IMAunp". Ignore for now. 

# Also merge with bgcs data for comparison. Reproject to lat/long as required by leaflet: 
bgcs2 <- merge(bgcs, subzones_colours_ref, by = "BGC")
bgcs2 <- st_transform(bgcs2, crs = 4326) 

# Add gaps: 
# First, reproject to lat/long: 
gap_poly2 <- project(gap_poly, "EPSG:4326")

# Leaflet:
leaflet(trainData_balanced_TEST) %>%
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

mismatched_predictions <- trainData_balanced_TEST[predictions_full != predictions_Wgaps]
mismatched_predictions <- trainData_balanced_TEST[predictions_full != predictions_Wgaps]
