## update WNA BGC with the new BEC13. 
## Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(sf)

# read in previous version of WNA BGC and reclassify minor units
bgcs <- st_read("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_15Nov2024.gpkg")
bgcs <- st_transform(bgcs,4326)
bgcs <- rasterize(bgcs, dem, field = "BGC")
bgcs <- subst(bgcs, "CMAwh", "CMAun")
bgcs <- subst(bgcs, "BWBScm", "BWBScmE")

# new BEC13
bgcs.v13 <- st_read("//objectstore2.nrs.bcgov/ffec/WNA_BGC/BEC13Draft_BEC13v1wLabelUpdates.shp")
bgcs.v13 <- st_transform(bgcs.v13,4326)
bgcs.v13 <- rasterize(bgcs.v13, dem, field = "BGC")
bgcs.updated <- bgcs.v13
bgcs.updated <- subst(bgcs.updated, "ICHmc1a", "ICHmc1")

# Extract category tables
bgcs_cats <- levels(bgcs)[[1]]
bgcs_updated_cats <- levels(bgcs.updated)[[1]]

# Merge unique categories by BGC names
all_levels <- unique(c(bgcs_cats$BGC, bgcs_updated_cats$BGC))

# Sort alphabetically by BGC name
all_levels <- sort(all_levels)

# Create a new factor level mapping with sorted BGC names
lookup <- data.frame(value = seq_along(all_levels), BGC = all_levels)

# Create a mapping from old values to new sorted values
bgcs_map <- merge(bgcs_cats, lookup, by = "BGC", all.x = TRUE, sort = FALSE)[, c("value.x", "value.y")]
bgcs_updated_map <- merge(bgcs_updated_cats, lookup, by = "BGC", all.x = TRUE, sort = FALSE)[, c("value.x", "value.y")]

# Rename columns for clarity
colnames(bgcs_map) <- colnames(bgcs_updated_map) <- c("old", "new")

# Convert NA to identity mapping (for any unmatched values)
bgcs_map$new[is.na(bgcs_map$new)] <- bgcs_map$old[is.na(bgcs_map$new)]
bgcs_updated_map$new[is.na(bgcs_updated_map$new)] <- bgcs_updated_map$old[is.na(bgcs_updated_map$new)]

# Reclassify rasters to use the sorted numeric values
bgcs <- classify(bgcs, bgcs_map, others = NA)
bgcs.updated <- classify(bgcs.updated, bgcs_updated_map, others = NA)

# Merge rasters
bgcs.merged <- cover(bgcs.updated, bgcs)

# Drop unused levels
used_values <- unique(na.omit(values(bgcs.merged, mat = FALSE)))  # Get used numeric codes
final_lookup <- lookup[lookup$value %in% used_values, ]           # Keep only relevant levels

# Assign cleaned and sorted category table
levels(bgcs.merged) <- list(final_lookup)

# Verify final categories
print(levels(bgcs.merged))

# mask out ocean areas
bgcs.merged <- mask(bgcs.merged, bgcs)
plot(bgcs.merged)

writeRaster(bgcs.merged, "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_22Feb2025_raster.tif", overwrite=TRUE)

