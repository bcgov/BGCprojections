library(leaflet)

# Combine data into a single data frame for plotting
trainData_balanced_TEST <- trainData_balanced

# Add in the predictions for the model version with no gaps/holdouts. 
trainData_balanced_TEST$predictions_full <- predictions_full

# Deal with gaps/holdouts: 
trainData_balanced_TEST$predictions_Wgaps <- NA
trainData_balanced_TEST$predictions_Wgaps[trainData_balanced_TEST$gaps == 1] <- predictions_Wgaps


# Create Leaflet map
leaflet(trainData_balanced) %>%
  addTiles() %>%
  # Add circles for predictions with holdout
  addCircleMarkers(
    lat = ~lat, lon = ~lon,
    color = ~colorNumeric("YlOrRd", predictions_withholdout)(predictions_withholdout),
    radius = 5, fillOpacity = 0.7,
    popup = ~paste("With Holdout:", predictions_withholdout)
  ) %>%
  addLegend(
    position = "bottomright",
    pal = colorNumeric("YlOrRd", trainData_balanced$predictions_withholdout),
    values = ~predictions_withholdout,
    title = "Predictions With Holdout"
  ) %>%
  addLayersControl(
    overlayGroups = c("Model with Holdout", "Model without Holdout"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  # Add circles for predictions without holdout
  addCircleMarkers(
    lat = ~lat, lon = ~lon,
    color = ~colorNumeric("YlGnBu", predictions_withoutholdout)(predictions_withoutholdout),
    radius = 5, fillOpacity = 0.7,
    popup = ~paste("Without Holdout:", predictions_withoutholdout),
    group = "Model without Holdout"
  ) %>%
  addLegend(
    position = "bottomright",
    pal = colorNumeric("YlGnBu", trainData_balanced$predictions_withoutholdout),
    values = ~predictions_withoutholdout,
    title = "Predictions Without Holdout"
  )
