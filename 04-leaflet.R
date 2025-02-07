#### Leaflet script: ####
# TO DO: 
# Read in predictions data: 
# Maybe make a separate 00 script with leaflet set up so that I can just run leaflet here. 
# Also try to simplify this code. 


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

# Leaflet:
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
