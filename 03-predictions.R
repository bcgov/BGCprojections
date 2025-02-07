# To do: 
# Read in models: 

#### Predictions: ####
# Rasterize the BGC polygons, assigning the values from the BGC column. 
bgcs_rast <- rasterize(bgcs, elev, field = "BGC")

# Align bgcs_rast and elev DEM to ensure resolution, extent, and CRS match: 
bgcs_rast <- resample(bgcs_rast, elev, method = "near")

# Merge bgcs_rast and elev into one multi-layer raster: 
bgcs_elev <- c(elev, bgcs_rast)
names(bgcs_elev) <- c("elev", "BGC")

# Reproject to lat/long to work with climr: 
bgcs_elev <- project(bgcs_elev, "EPSG:4326", method = "near")

# Make into data table: 
bgcs_elev_dt <- as.data.table(bgcs_elev, cells = TRUE, xy = TRUE) 
colnames(bgcs_elev_dt) <- c("id", "lon", "lat", "elev", "BGC")

# Get climr data for these raster coordinates:
clim_vars_preds2 <- downscale(
  xyz = bgcs_elev_dt,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars(),
  cache = TRUE)|>
  Cache()

clim_vars_preds <- copy(clim_vars_preds2)
ccissr::addVars(clim_vars_preds)

# Merge important info (lat, lon, elev, and BGC) back in: 
clim_vars_preds <- merge(clim_vars_preds, bgcs_elev_dt, by = "id") 

# Predictions based on full model and gap removed models: 
preds_full_simple <- predict(BGCmodel_full_simple, data = clim_vars_preds)
preds_Wgaps_simple <- predict(BGCmodel_Wgaps_simple, data = clim_vars_preds)
preds_full_expert <- predict(BGCmodel_full_expert, data = clim_vars_preds)
preds_Wgaps_expert <- predict(BGCmodel_Wgaps_expert, data = clim_vars_preds)
preds_full_all <- predict(BGCmodel_full_all, data = clim_vars_preds)
preds_Wgaps_all <- predict(BGCmodel_Wgaps_all, data = clim_vars_preds)

# Save into clim_vars_preds dataframe:
clim_vars_preds$preds_full_simple <- as.character(preds_full_simple$predictions)
clim_vars_preds$preds_Wgaps_simple <- as.character(preds_Wgaps_simple$predictions)
clim_vars_preds$preds_full_expert <- as.character(preds_full_expert$predictions)
clim_vars_preds$preds_Wgaps_expert <- as.character(preds_Wgaps_expert$predictions)
clim_vars_preds$preds_full_all <- as.character(preds_full_all$predictions)
clim_vars_preds$preds_Wgaps_all <- as.character(preds_Wgaps_all$predictions)

#### Check quality of predictions: ####
conf_matrix_preds_full_simple <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_full_simple, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

conf_matrix_preds_Wgaps_simple <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_Wgaps_simple, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

conf_matrix_preds_full_expert <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_full_expert, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

conf_matrix_preds_Wgaps_expert <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_Wgaps_expert, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

conf_matrix_preds_full_all <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_full_all, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

conf_matrix_preds_Wgaps_all <- caret::confusionMatrix(
  data = factor(clim_vars_preds$preds_Wgaps_all, levels = levels(clim_vars_preds$BGC)),
  reference = factor(clim_vars_preds$BGC, levels = levels(clim_vars_preds$BGC))
)

sensitivity_values <- cm$byClass[, "Sensitivity"]
class_labels <- rownames(cm$byClass)

#### Look at feature importance: ####
# importance_scores_full_simple <- BGCmodel_full_simple$variable.importance
# importance_scores_full_simple<- sort(importance_scores_full_simple, decreasing = TRUE)
# 
# importance_scores_Wgaps_simple <- BGCmodel_Wgaps_simple$variable.importance
# importance_scores_Wgaps_simple<- sort(importance_scores_Wgaps_simple, decreasing = TRUE)
# 
# # Create a data frame for plotting
# importance_scores_full_simple_df <- data.frame(
#   Feature = names(importance_scores_full_simple),
#   Importance = importance_scores_full_simple
# )
# 
# importance_scores_Wgaps_simple_df <- data.frame(
#   Feature = names(importance_scores_Wgaps_simple),
#   Importance = importance_scores_Wgaps_simple
# )
# 
# # Plots
# ggplot(importance_scores_full_simple_df, aes(x = reorder(Feature, Importance), y = Importance)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   coord_flip() +
#   labs(title = "Feature Importance",
#        x = "Features",
#        y = "Importance") +
#   theme_minimal()
# 
# ggplot(importance_scores_Wgaps_simple_df, aes(x = reorder(Feature, Importance), y = Importance)) +
#   geom_bar(stat = "identity", fill = "steelblue") +
#   coord_flip() +
#   labs(title = "Feature Importance",
#        x = "Features",
#        y = "Importance") +
#   theme_minimal()
