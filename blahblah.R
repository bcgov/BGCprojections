library(terra)
source("R/utils.R")


# Create new file paths for Deb's temporary data location: 
bgcs <- vect("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/BGCProjections/Proxy_ObjectStorage/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022/WNA_BGC_v12_5Apr2022.gpkg")

# elev <- rast("//objectstore2.nrs.bcgov/ffec/DEM/DEM_NorAm/NA_Elevation/data/northamerica/northamerica_elevation_cec_2023.tif")
# Use 800 m DEM instead: 
elev <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/BGCProjections/Proxy_ObjectStorage/DEM/DEM_Composite_WNA_800m/composite_WNA_dem.tif")

# Reproject both to lat/long: 
bgcs <- project(bgcs, elev)

## make subsets of the study area for hold-outs (gaps):
# This code makes rectangular gap extents within the bounding box of the study area. 5L is the default number of gaps to create. 
studyarea <- ext(c(-123, -117, 49, 52.5))
elev <- crop(elev, studyarea)
bgcs <- crop(bgcs, studyarea)
gapextents <- makeGapExtents(studyarea, 5L)
plot(temp)
for(i in 1:5){
  plot(gapextents[[i]], add=T)
}
