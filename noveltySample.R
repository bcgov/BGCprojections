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

# ---------------------------------
# testing
library(climr)
clim.pts <- downscale(xyz = points_subsample,
                      vars = list_vars())
clim.pts <- points_subsample[clim.pts, on = "id"]

# Calculate the centroid climate for the training points
clim.pts.mean <- clim.pts[, lapply(.SD, mean), by = BGC, .SDcols = -c("id", "PERIOD")]

par(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.75, 0.25, 0), tck=-0.01)
x <- clim.pts$MAT
y <- log2(clim.pts$MAP)
x1 <- clim.pts.mean$MAT
y1 <- log2(clim.pts.mean$MAP)
x2 <- clim.pts[BGC=="CWHxm_WA", MAT]
y2 <- log2(clim.pts[BGC=="CWHxm_WA", MAP])
plot(x, y, pch=16, col="grey", cex=0.3)
text(x1,y1, clim.pts.mean$BGC, cex=0.5)
points(x2,y2, pch=16, cex=0.5, col="blue")

