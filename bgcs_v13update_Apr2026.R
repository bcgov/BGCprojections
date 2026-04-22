## create rasterized WNA BGCs 
## Colin Mahony colin.mahony@gov.bc.ca

library(terra)
library(sf)
library(data.table)

# import DEM
dem <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling//WNA_DEM_4326_clipped.tif")

# rasterize
bgcs <- st_read("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_13Apr2026.gpkg")
bgcs <- st_transform(bgcs,4326)
bgcs <- rasterize(bgcs, dem, field = "BGC")

# ---------------------------------
# substitute labels for very small or mislabeled units

# mapping
ct <- cats(bgcs)[[1]]
rcl <- cbind(
  ct$value[match(c("CMAwh","BWBScm","ICHmc1a"), ct$BGC)],
  ct$value[match(c("CMAun","BWBScmE","ICHmc1"), ct$BGC)]
)

# reclassify
bgcs2 <- classify(bgcs, rcl,
                  filename="bgcs_reclass.tif",
                  overwrite=TRUE)

# restore categorical structure
ct2 <- cats(bgcs)[[1]]
for (i in seq_len(nrow(lookup))) {
  ct2$BGC[ct2$BGC == lookup$from[i]] <- lookup$to[i]
}
ct2 <- ct2[!duplicated(ct2$BGC), ]
bgcs2 <- as.factor(bgcs2)
levels(bgcs2) <- ct2









# 1. mapping
ct <- cats(bgcs)[[1]]

rcl <- cbind(
  ct$value[match(c("CMAwh","BWBScm","ICHmc1a"), ct$BGC)],
  ct$value[match(c("CMAun","BWBScmE","ICHmc1"), ct$BGC)]
)

# 2. reclassify
bgcs2 <- classify(bgcs, rcl,
                  filename="bgcs_reclass.tif",
                  overwrite=TRUE)

# 3. restore categorical structure (IMPORTANT: use ORIGINAL cats)
vals <- unique(values(bgcs2))
ct <- cats(bgcs)[[1]]
ct2$BGC[ct2$BGC == "CMAwh"]  <- "CMAun"
ct2$BGC[ct2$BGC == "BWBScm"] <- "BWBScmE"
ct2$BGC[ct2$BGC == "ICHmc1a"]<- "ICHmc1"
ct <- ct[ct$value %in% vals, ]
levels(bgcs2) <- ct


# ---------------------------------
# write raster and clean up disk

writeRaster(bgcs2, "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_13Apr2026_raster.tif", overwrite=TRUE)

rm(bgcs2)
gc()
file.remove("bgcs_reclass.tif")
