## BGC model for WNA
## Colin Mahony colin.mahony@gov.bc.ca
## February 2025

library(climr)
library(data.table)
library(terra)
library(sf)
library(foreach) # for outlier removal function

# Source functions TODO: move these into ccissr as independent functions
source("utils.R")
source("bgc_trainingSample.R")

studyarea <- ext(c(-123, -122, 50.5, 51)) 

# import DEM
dem <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling//WNA_DEM_4326_clipped.tif")
# dem <- crop(dem, studyarea)
X <- dem # for plotting

# import BGCs
bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_22Feb2025_raster.tif")
# bgcs <- crop(bgcs, studyarea)

bgcs_info <- fread("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNA_BGCs_Info_v13_1.csv")

trainingSample_V1 <- bgc_trainingSample(dem, bgcs, bgcs_info = bgcs_info,
                           scheme = "asymptotic", asymptote = 2000, 
                           removeOutliers = TRUE, alpha.BC = .0027, alpha.nonBC = .05, 
                           removeOutliers.vars = c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt"), 
                           climaticVariance = TRUE, climaticVariance.var = "MAT",
                           bgcs.remove = "MSSDun_NV",
                           plotDiagnostics = TRUE, 
                           plot.dir = "//objectstore2.nrs.bcgov/ffec/BGC_models", 
                           plot.name = "diagnostics_v1"
                           )
dim(trainingSample_V1)
write.csv(trainingSample_V1, "//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v1.csv", row.names = FALSE)

trainingSample_V2 <- bgc_trainingSample(dem, bgcs, bgcs_info = bgcs_info,
                                        scheme = "squareRoot", 
                                        squareRoot.multiplier = 5, 
                                        removeOutliers = TRUE, alpha.BC = .0027, alpha.nonBC = .05, 
                                        removeOutliers.vars = c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt"), 
                                        climaticVariance = TRUE, climaticVariance.var = "MAT",
                                        bgcs.remove = "MSSDun_NV",
                                        plotDiagnostics = TRUE, 
                                        plot.dir = "//objectstore2.nrs.bcgov/ffec/BGC_models", 
                                        plot.name = "diagnostics_v2"
)
dim(trainingSample_V2)
write.csv(trainingSample_V2, "//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v2.csv", row.names = FALSE)

trainingSample_V3 <- bgc_trainingSample(dem, bgcs, bgcs_info = bgcs_info,
                                        scheme = "asymptotic", asymptote = 2000, 
                                        removeOutliers = TRUE, alpha.BC = .0027, alpha.nonBC = .05, 
                                        removeOutliers.vars = c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt"), 
                                        climaticVariance = FALSE,
                                        bgcs.remove = "MSSDun_NV",
                                        plotDiagnostics = TRUE, 
                                        plot.dir = "//objectstore2.nrs.bcgov/ffec/BGC_models", 
                                        plot.name = "diagnostics_v3"
)
dim(trainingSample_V3)
write.csv(trainingSample_V3, "//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v3.csv", row.names = FALSE)

trainingSample_V4 <- bgc_trainingSample(dem, bgcs, bgcs_info = bgcs_info,
                                        scheme = "squareRoot", 
                                        squareRoot.multiplier = 10, 
                                        removeOutliers = TRUE, alpha.BC = .0027, alpha.nonBC = .05, 
                                        removeOutliers.vars = c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt"), 
                                        climaticVariance = TRUE, climaticVariance.var = "MAT",
                                        bgcs.remove = "MSSDun_NV",
                                        plotDiagnostics = TRUE, 
                                        plot.dir = "//objectstore2.nrs.bcgov/ffec/BGC_models", 
                                        plot.name = "diagnostics_v4"
)
dim(trainingSample_V4)
write.csv(trainingSample_V4, "//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v4.csv", row.names = FALSE)
trainingSample_V4 <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v4.csv")


## -------------------------------------------------
## -------------------------------------------------
## STEP 2 - Train RF model
## -------------------------------------------------
## -------------------------------------------------

library(climr)
library(data.table)
library(terra)
library(ranger) # For RF
library(caret) # For confusionMatrix()
library(leaflet)

## -------------------------------------------------
## define variable set
vars_expert <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", 
                 "EXT", "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", 
                 "SHM", "Tave_sm", "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", 
                 "Tmin_at", "Tmin_sm", "Tmin_sp", "Tmin_wt", "CMI", "PPT_MJ", 
                 "PPT_JAS", "CMD.total")

## -------------------------------------------------
## V1 RF model

#read in training sample generated in the last step
points <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v1.csv")

## climate data for all points
clim <- downscale(
  xyz = points,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim)

trainData <- merge(points, clim, by="id")

trainData[, BGC := as.factor(BGC)]

# Train model 
BGCmodel <- ranger(
  BGC ~ .,
  data = trainData[, c("BGC", vars_expert), with = FALSE],
  num.trees = 250,
  splitrule =  "extratrees",
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 14,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  keep.inbag = FALSE, 
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, "//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V1.rds") 

## -------------------------------------------------
## V2 RF model

#read in training sample generated in the last step
points <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v2.csv")

## climate data for all points
clim <- downscale(
  xyz = points,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim)

trainData <- merge(points, clim, by="id")

trainData[, BGC := as.factor(BGC)]

# Train model 
BGCmodel <- ranger(
  BGC ~ .,
  data = trainData[, c("BGC", vars_expert), with = FALSE],
  num.trees = 250,
  splitrule =  "extratrees",
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 14,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  keep.inbag = FALSE, 
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, "//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V2.rds") 

# Train model 
BGCmodel <- ranger(
  BGC ~ .,
  data = trainData[, c("BGC", vars_expert), with = FALSE],
  num.trees = 500,
  splitrule =  "extratrees",
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 14,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  keep.inbag = FALSE, 
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, "//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V2.1.rds") 

## -------------------------------------------------
## V3 RF model

#read in training sample generated in the last step
points <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v3.csv")

## climate data for all points
clim <- downscale(
  xyz = points,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim)

trainData <- merge(points, clim, by="id")

trainData[, BGC := as.factor(BGC)]

# Train model 
BGCmodel <- ranger(
  BGC ~ .,
  data = trainData[, c("BGC", vars_expert), with = FALSE],
  num.trees = 250,
  splitrule =  "extratrees",
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 14,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  keep.inbag = FALSE, 
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, "//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V3.rds") 

## -------------------------------------------------
## V4 RF model

#read in training sample generated in the last step
points <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v4.csv")

## climate data for all points
clim <- downscale(
  xyz = points,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim)

trainData <- merge(points, clim, by="id")

trainData[, BGC := as.factor(BGC)]

# Train model 
BGCmodel <- ranger(
  BGC ~ .,
  data = trainData[, c("BGC", vars_expert), with = FALSE],
  num.trees = 250,
  splitrule =  "extratrees",
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 14,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  keep.inbag = FALSE, 
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, "//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V4.rds") 

## -------------------------------------------------
## V4.1 RF model

#read in training sample generated in the last step
points <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/points_WNA_v4.csv")

## climate data for all points
clim <- downscale(
  xyz = points,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim)

trainData <- merge(points, clim, by="id")

trainData[, BGC := as.factor(BGC)]

# Train model 
BGCmodel <- ranger(
  BGC ~ .,
  data = trainData[, c("BGC", vars_expert), with = FALSE],
  num.trees = 500,
  splitrule =  "extratrees",
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 7,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  keep.inbag = FALSE, 
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, "//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V4.1.rds") 

## -------------------------------------------------
## -------------------------------------------------
## STEP 3 - Results
## -------------------------------------------------
## -------------------------------------------------


library(climr)
library(data.table)
library(terra)
library(leaflet)
library(ranger) # For RF

## -------------------------------------------------
## color scheme for bgc units
subzones_colours_ref <- fread("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNAv13_Subzone_colours_2.csv") %>% 
  dplyr::select(!c(fid, MAP_LABEL, NSRNAME, ZONE)) %>% 
  dplyr::mutate(BGC_num = as.numeric(as.factor(BGC)))
color_pal <- colorFactor(
  palette = subzones_colours_ref$RGB,
  domain = subzones_colours_ref$BGC_num
)

# data.table to store results of sensitivity analyses
results.error <- data.table(
  model       = c("v1", "v2", "v2.1", "v3", "v4", "nov2024"),
  Bamfield    = NA_real_,
  Kamloops    = NA_real_,
  Pemberton   = NA_real_,
  Smithers   = NA_real_,
  BC   = NA_real_
)

## -------------------------------------------------
## alternative RF models

## Courtney's model
loaded_objects <- load("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/Trained_Models/BGC_RFresp.Rdata")
print(loaded_objects)

BGCmodel_v1 <- readRDS("//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V1.rds")
BGCmodel_v2 <- readRDS("//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V2.rds")
BGCmodel_v2.1 <- readRDS("//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V2.1.rds")
BGCmodel_v3 <- readRDS("//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V3.rds")
BGCmodel_v4 <- readRDS("//objectstore2.nrs.bcgov/ffec/BGC_models/BGCmodel_WNA_V4.rds")


## -------------------------------------------------
## case study areas
studynames <- c("Bamfield", "Kamloops", "Pemberton", "Smithers", "BC")
studyname <- "BC"
for(studyname in studynames){
  studyarea <- if(studyname=="Bamfield") ext(c(-125.25, -124, 48.5, 49.125)) else 
    if(studyname=="Kamloops") ext(c(-121, -120, 50.5, 51)) else 
      if(studyname=="Pemberton") ext(c(-124, -122, 50, 51)) else 
        if(studyname=="Smithers") ext(c(-129, -126, 54.25, 55.25)) else 
          NULL
  
if(studyname == "BC"){
  bdy <- vect("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/bdy.BC.shp")
  bdy <- project(bdy, "EPSG:4326")
  dem <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_DEM_4326_clipped.tif")
  dem <- aggregate(dem, fact=9)
  dem <- mask(dem, bdy) 
  dem <- crop(dem, ext(bdy) )
  plot(dem)
  X <- dem # template raster for testing
  values(X) <- NA
  
  bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_22Feb2025_raster.tif")
  bgcs <- crop(bgcs, ext(bdy))
  bgcs <- project(bgcs, dem, method = "near")
  bgcs <- mask(bgcs, dem) 
  bgc_levels <- cats(bgcs)[[1]]  # Extract the category mapping
  vals_char <- bgc_levels$BGC[match(values(bgcs), bgc_levels$value)] # Convert numeric values to character labels
  values(bgcs) <- factor(vals_char, levels = subzones_colours_ref$BGC) # refactor with full set of levels
  bgcs_3857 <- project(bgcs, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
} else {
  # study area DEM
  dem <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_DEM_4326_clipped.tif")
  dem <- crop(dem, studyarea)
  plot(dem)
  X <- dem # template raster for testing
  values(X) <- NA
  
  # study area bgcs
  bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_22Feb2025_raster.tif")
  bgcs <- crop(bgcs, studyarea)
  bgc_levels <- cats(bgcs)[[1]]  # Extract the category mapping
  vals_char <- bgc_levels$BGC[match(values(bgcs), bgc_levels$value)] # Convert numeric values to character labels
  values(bgcs) <- factor(vals_char, levels = subzones_colours_ref$BGC) # refactor with full set of levels
  bgcs_3857 <- project(bgcs, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
  
}

## study area points
points <- as.data.table(dem, cells=T, xy=T)
colnames(points) <- c("id", "lon", "lat", "elev")
points <- points[,c(2,3,4,1)] #restructure for climr input
values(X) <- NA; values(X)[points$id] <- points$el ; plot(X)

## -------------------------------------------------
## climate data for all points

clim_ref <- downscale(
  xyz = points,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim_ref)

clim_proj <- downscale(
  xyz = points,
  gcms = list_gcms()[4],
  ssps = list_ssps()[2],
  gcm_periods = list_gcm_periods()[3],
  run_nm = list_runs_ssp(list_gcms()[4], list_ssps()[2])[3],
  which_refmap = "refmap_climr",
  return_refperiod = FALSE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim_proj)

## -------------------------------------------------
## predictions

## v1 model predictions for reference period
preds_ref_v1_vec <- predict(BGCmodel_v1, data = clim_ref)$prediction
preds_ref_v1 <- X
preds_ref_v1[points$id] <- factor(preds_ref_v1_vec, levels = subzones_colours_ref$BGC)
preds_ref_v1 <- project(preds_ref_v1, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v1 model predictions for future period
preds_proj_v1_vec <- predict(BGCmodel_v1, data = clim_proj)$prediction
preds_proj_v1 <- X
preds_proj_v1[points$id] <- factor(preds_proj_v1_vec, levels = subzones_colours_ref$BGC)
preds_proj_v1 <- project(preds_proj_v1, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v2 model predictions for reference period
preds_ref_v2_vec <- predict(BGCmodel_v2, data = clim_ref)$prediction
preds_ref_v2 <- X
preds_ref_v2[points$id] <- factor(preds_ref_v2_vec, levels = subzones_colours_ref$BGC)
preds_ref_v2 <- project(preds_ref_v2, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v2 model predictions for future period
preds_proj_v2_vec <- predict(BGCmodel_v2, data = clim_proj)$prediction
preds_proj_v2 <- X
preds_proj_v2[points$id] <- factor(preds_proj_v2_vec, levels = subzones_colours_ref$BGC)
preds_proj_v2 <- project(preds_proj_v2, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v2 model predictions for reference period
preds_ref_v2.1_vec <- predict(BGCmodel_v2.1, data = clim_ref)$prediction
preds_ref_v2.1 <- X
preds_ref_v2.1[points$id] <- factor(preds_ref_v2.1_vec, levels = subzones_colours_ref$BGC)
preds_ref_v2.1 <- project(preds_ref_v2.1, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v2.1 model predictions for future period
preds_proj_v2.1_vec <- predict(BGCmodel_v2.1, data = clim_proj)$prediction
preds_proj_v2.1 <- X
preds_proj_v2.1[points$id] <- factor(preds_proj_v2.1_vec, levels = subzones_colours_ref$BGC)
preds_proj_v2.1 <- project(preds_proj_v2.1, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v3 model predictions for reference period
preds_ref_v3_vec <- predict(BGCmodel_v3, data = clim_ref)$prediction
preds_ref_v3 <- X
preds_ref_v3[points$id] <- factor(preds_ref_v3_vec, levels = subzones_colours_ref$BGC)
preds_ref_v3 <- project(preds_ref_v3, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v3 model predictions for future period
preds_proj_v3_vec <- predict(BGCmodel_v3, data = clim_proj)$prediction
preds_proj_v3 <- X
preds_proj_v3[points$id] <- factor(preds_proj_v3_vec, levels = subzones_colours_ref$BGC)
preds_proj_v3 <- project(preds_proj_v3, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v4 model predictions for reference period
preds_ref_v4_vec <- predict(BGCmodel_v4, data = clim_ref)$prediction
preds_ref_v4 <- X
preds_ref_v4[points$id] <- factor(preds_ref_v4_vec, levels = subzones_colours_ref$BGC)
preds_ref_v4 <- project(preds_ref_v4, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v4 model predictions for future period
preds_proj_v4_vec <- predict(BGCmodel_v4, data = clim_proj)$prediction
preds_proj_v4 <- X
preds_proj_v4[points$id] <- factor(preds_proj_v4_vec, levels = subzones_colours_ref$BGC)
preds_proj_v4 <- project(preds_proj_v4, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## nov2024 model predictions for reference period
preds_ref_nov2024_vec <- predict(BGC_RFresp, data = clim_ref)$prediction
preds_ref_nov2024 <- X
preds_ref_nov2024[points$id] <- factor(preds_ref_nov2024_vec, levels = subzones_colours_ref$BGC)
preds_ref_nov2024 <- project(preds_ref_nov2024, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries.

## nov2024 model predictions for future period
preds_proj_nov2024_vec <- predict(BGC_RFresp, data = clim_proj)$prediction
preds_proj_nov2024 <- X
preds_proj_nov2024[points$id] <- factor(preds_proj_nov2024_vec, levels = subzones_colours_ref$BGC)
preds_proj_nov2024 <- project(preds_proj_nov2024, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries.

## -------------------------------------------------
## leaflet map

map <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%  # Add satellite imagery
  addTiles(group = "OSM Basemap") %>%  # Keep OpenStreetMap as an option
  addRasterImage(preds_ref_v1, colors = color_pal, opacity = 1, group = "v1: baseline") %>%
  addRasterImage(preds_proj_v1, colors = color_pal, opacity = 1, group = "v1: future") %>%
  addRasterImage(preds_ref_v2, colors = color_pal, opacity = 1, group = "v2: baseline") %>%
  addRasterImage(preds_proj_v2, colors = color_pal, opacity = 1, group = "v2: future") %>%
  addRasterImage(preds_ref_v2.1, colors = color_pal, opacity = 1, group = "v2.1: baseline") %>%
  addRasterImage(preds_proj_v2.1, colors = color_pal, opacity = 1, group = "v2.1: future") %>%
  addRasterImage(preds_ref_v3, colors = color_pal, opacity = 1, group = "v3: baseline") %>%
  addRasterImage(preds_proj_v3, colors = color_pal, opacity = 1, group = "v3: future") %>%
  addRasterImage(preds_ref_v4, colors = color_pal, opacity = 1, group = "v4: baseline") %>%
  addRasterImage(preds_proj_v4, colors = color_pal, opacity = 1, group = "v4: future") %>%
  # addRasterImage(preds_ref_nov2024, colors = color_pal, opacity = 1, group = "nov2024: baseline") %>%
  # addRasterImage(preds_proj_nov2024, colors = color_pal, opacity = 1, group = "nov2024: future") %>%
  addRasterImage(bgcs_3857, colors = color_pal, opacity = 1, group = "BGC") %>%
  addLayersControl(
    baseGroups = c("Satellite", "OSM Basemap"),  # Base layer switcher
    overlayGroups = c("v1: baseline", "v1: future", 
                      "v2: baseline", "v2: future", 
                      "v2.1: baseline", "v2.1: future", 
                      "v3: baseline", "v3: future", 
                      "v4: baseline", "v4: future", 
                      # "nov2024: baseline", "nov2024: future", 
                      "BGC"),
    options = layersControlOptions(collapsed = FALSE)
  )
print(map)

## -------------------------------------------------
## Compute classification error

results.error[1, which(names(results.error)==studyname)] <- mean(values(bgcs_3857, mat = FALSE) != values(preds_ref_v1, mat = FALSE), na.rm=T)
results.error[2, which(names(results.error)==studyname)] <- mean(values(bgcs_3857, mat = FALSE) != values(preds_ref_v2, mat = FALSE), na.rm=T)
results.error[3, which(names(results.error)==studyname)] <- mean(values(bgcs_3857, mat = FALSE) != values(preds_ref_v2.1, mat = FALSE), na.rm=T)
results.error[4, which(names(results.error)==studyname)] <- mean(values(bgcs_3857, mat = FALSE) != values(preds_ref_v3, mat = FALSE), na.rm=T)
results.error[5, which(names(results.error)==studyname)] <- mean(values(bgcs_3857, mat = FALSE) != values(preds_ref_v4, mat = FALSE), na.rm=T)
results.error[6, which(names(results.error)==studyname)] <- mean(values(bgcs_3857, mat = FALSE) != values(preds_ref_nov2024, mat = FALSE), na.rm=T)

print(studyname)
}
write.csv(results.error, "//objectstore2.nrs.bcgov/ffec/BGC_models/results.error.csv", row.names = FALSE)
results.error <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/results.error.csv")
results.error[, BC := NA_real_]

# Plot of error results
plot.dir = "//objectstore2.nrs.bcgov/ffec/BGC_models"
plot.name = "SamplingTrials_error"
png(filename=paste0(plot.dir, "/", plot.name, ".png",sep="."), type="cairo", units="in", width=6.5, height=4, pointsize=10, res=300)
error_matrix <- as.matrix(results.error[, -1, with = FALSE])
rownames(error_matrix) <- results.error$model
barplot(error_matrix, beside = TRUE, col = rainbow(nrow(error_matrix)), 
        legend.text = rownames(error_matrix), args.legend = list(x = "topleft", bty="n"),
        main = "Prediction error of alternative BGC models",
        xlab = "Case Study", ylab = "Error")
box()
dev.off()
