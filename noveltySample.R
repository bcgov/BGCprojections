## generate a sample of points for each BGC to use in novel climates measurement
## Colin Mahony colin.mahony@gov.bc.ca
## February 2025

# remotes::install_github("bcgov/ccissr@development")

library(ccissr)
library(terra)
library(data.table)

# import DEM
dem <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling//WNA_DEM_4326_clipped.tif")

# import BGCs
bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_22Feb2025_raster.tif")

points_novelty <- bgc_trainingSample(dem, bgcs,
                                        scheme = "asymptotic", asymptote = 200, shape=1,
                                        plotDiagnostics = FALSE
)
dim(points_novelty)
table(points_novelty$BGC)
write.csv(points_novelty, "//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_simple200.csv", row.names = FALSE)
