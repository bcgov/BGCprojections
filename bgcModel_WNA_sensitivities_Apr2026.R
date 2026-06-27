## BGC model for WNA
## Colin Mahony colin.mahony@gov.bc.ca
## April 2026
## This is an update for the refined US BEC (completed march 2026). the approach is to use the model selected based on the Feb2025 sensitivity analysis. 

library(climr)
library(ccissr)
library(data.table)
library(terra)
library(sf)
library(foreach) # for outlier removal function

# dir <- "//objectstore2.nrs.bcgov/ffec/BGC_models/" 
dir <- "C:/Users/CMAHONY/Data/BGC_models/" #local copy, for speed

# Source functions TODO: move these into ccissr as independent functions
source("utils.R")

studyarea <- ext(c(-123, -122, 50.5, 51)) 

# import DEM
dem <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling//WNA_DEM_4326_clipped.tif")
# dem <- crop(dem, studyarea)
X <- dem # for plotting

# import BGCs
bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_13Apr2026_raster.tif")
# bgcs <- crop(bgcs, studyarea)

# BGC metadata table. TODO: modify this chunk when WNA_BGCs.csv is made into a ccissr data object. 
wna_bgcs <- fread("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNA_BGCs.csv")
ct <- cats(bgcs)[[1]]
ct <- ct[!is.na(ct$value), ] # keep only rows with real values (exclude the NA class)
ct <- ct[ct$BGC != "(None)", ]
bgcs_info <- data.table(Source   = NA_character_, BGC = ct$BGC)
bgcs_info[, Source := fifelse(grepl("_(WA|OR|CA|ID|NV|UT|CO|WY|MT|WC|OC)$", BGC), "USA_", NA_character_)]
bgcs_info[wna_bgcs[DataSet == "AB", .(BGC, DataSet)],  on = .(BGC = BGC), Source := "AB_"]
bgcs_info[is.na(Source), Source := "BC_"]
setorder(bgcs_info, Source, BGC)
f <- freq(bgcs)
setorder(f, count)
f
tail(f,50)
setorder(f, value)
f

# # aggregate because the raster is too big for the bgc_trainingSample() memory requirements
# bgcs <- aggregate(bgcs, fact=2)
# levels(bgcs) <- ct
# dem <- aggregate(dem, fact=2)

trainingSample_V2 <- bgc_trainingSample(dem, bgcs, bgcs_info = bgcs_info,
                                        scheme = "squareRoot", 
                                        squareRoot.multiplier = 5, 
                                        removeOutliers = TRUE, alpha.BC = .0027, alpha.nonBC = .05, 
                                        removeOutliers.vars = c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt"), 
                                        climaticVariance = TRUE, climaticVariance.var = "MAT",
                                        plotDiagnostics = TRUE, 
                                        plot.dir = "//objectstore2.nrs.bcgov/ffec/BGC_models", 
                                        plot.name = "diagnostics_v2"
)
dim(trainingSample_V2)

write.csv(trainingSample_V2, paste0(dir, "points_WNA_v2.csv"), row.names = FALSE)

trainingSample_V4 <- bgc_trainingSample(dem, bgcs, bgcs_info = bgcs_info,
                                        scheme = "squareRoot", 
                                        squareRoot.multiplier = 10, 
                                        removeOutliers = TRUE, alpha.BC = .0027, alpha.nonBC = .05, 
                                        removeOutliers.vars = c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt"), 
                                        climaticVariance = TRUE, climaticVariance.var = "MAT",
                                        plotDiagnostics = TRUE, 
                                        plot.dir = "//objectstore2.nrs.bcgov/ffec/BGC_models", 
                                        plot.name = "diagnostics_v4"
)
dim(trainingSample_V4)

write.csv(trainingSample_V4, paste0(dir, "points_WNA_v4.csv"), row.names = FALSE)

# downsample non-BC units (using the default n=2000 asymptote)
trainingSample_V4 <- fread(paste0(dir, "points_WNA_v4.csv"))
bgcs_nonBC <- bgcs_info[grep("USA_|AB_", bgcs_info$Source), BGC]
# #simple downsampling
# trainingSample_V4a <- trainingSample_V4[
#   , if (BGC[1] %in% bgcs_nonBC) .SD[sample(.N, ceiling(.N * 0.2))]
#   else .SD,
#   by = BGC
# ]
#asymptotic downsampling
trainingSample_V4a <- trainingSample_V4[
  , {
    if (BGC[1] %in% bgcs_nonBC) {
      n_keep <- subsample_asymptotic(.N)
      if (n_keep <= 0) return(NULL)
      .SD[sample(.N, n_keep)]
    } else {
      .SD
    }
  },
  by = BGC
]

write.csv(trainingSample_V4a, paste0(dir, "points_WNA_v4a.csv"), row.names = FALSE)

dim(trainingSample_V4a)
sort(table(trainingSample_V4a$BGC))
sort(table(trainingSample_V4$BGC))

trainingSample_V5 <- bgc_trainingSample(dem, bgcs, bgcs_info = bgcs_info,
                                        scheme = "squareRoot", 
                                        squareRoot.multiplier = 20, 
                                        removeOutliers = TRUE, alpha.BC = .0027, alpha.nonBC = .05, 
                                        removeOutliers.vars = c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt"), 
                                        climaticVariance = TRUE, climaticVariance.var = "MAT",
                                        bgcs.remove = "MSSDun_NV",
                                        plotDiagnostics = TRUE, 
                                        plot.dir = "//objectstore2.nrs.bcgov/ffec/BGC_models", 
                                        plot.name = "diagnostics_v5"
)
dim(trainingSample_V5)

write.csv(trainingSample_V5, paste0(dir, "points_WNA_v5.csv"), row.names = FALSE)

## -------------------------------------------------
## technical report figure of training sample

points <- fread(paste0(dir, "points_WNA_v4a.csv"))
popn <- setDT(freq(bgcs))
head(popn)
samp <- as.data.table(table(points$BGC))
samp

setnames(samp, "V1", "BGC")
samp[, count := popn[.SD, on = .(value = BGC), x.count]]
samp[, count := popn[.SD, on = .(value = BGC), x.count]]

bgcs_nonBC <- bgcs_info[grep("USA_|AB_", bgcs_info$Source), BGC]


figdir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(figdir, "sampleSize_v4a.png",sep="."), type="cairo", units="in", width=6.5, height=6.5, pointsize=12, res=300)
par(mfrow=c(1,1), mar = c(3,4,1,1), mgp=c(1.75,0.25,0), tck=-0.005)

# plot setup
x <- log10(samp[, count])
y <- log2(samp[, N])
z <- samp[, BGC]
plot(x, y, col = "white",
     xaxt = "n", yaxt = "n",
     ylab = "",
     xlab = paste0("BGC unit area (number of cells)"),
)
axis(1, at = seq(1,20), labels = round(10^seq(1,20)))
axis(2, at = 1:99, labels = 2^(1:99), las=2)
par(mgp=c(2.75,0.25,0))
title(ylab = "Sample size of BGC Unit")

# line of square-root proportionality
z <- 10^seq(0, 7, 0.01)
x <- log10(z)
y <- log2(z^0.5*10)
lines(x,y, col="gray")

# plot of BGC subsample size vs BGC area - BC UNITS
x <- log10(samp[-which(BGC %in% bgcs_nonBC), count])
y <- log2(samp[-which(BGC %in% bgcs_nonBC), N])
z <- samp[-which(BGC %in% bgcs_nonBC), BGC]
text(x, y, labels = z, cex = 0.5)

# plot of BGC subsample size vs BGC area - NON-BC UNITS
x <- log10(samp[which(BGC %in% bgcs_nonBC), count])
y <- log2(samp[which(BGC %in% bgcs_nonBC), N])
z <- samp[which(BGC %in% bgcs_nonBC), BGC]
text(x, y, labels = z, cex = 0.5, col="blue")

legend("topleft", legend=c("BC units", "Non-BC units"), fill = c("black", "blue"), bty="n")

dev.off()
par(mfrow=c(1,1))


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
                 "EXT", "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", 
                 "SHM", "Tave_sm", "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", 
                 "Tmin_at", "Tmin_sm", "Tmin_sp", "Tmin_wt", "CMI_an", "PPT_MJ", 
                 "PPT_JAS", "CMD.total")
## -------------------------------------------------
## V2.2 RF model

#read in training sample generated in the last step
points <- fread(paste0(dir, "points_WNA_v2.csv"))

## climate data for all points
clim <- downscale(
  xyz = points,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim)

trainData <- merge(points, clim, by="id")

# remove points that climr does not return full data for 
num_cols <- names(trainData)[sapply(trainData, is.numeric)]
bad_rows <- trainData[, !Reduce(`&`, lapply(.SD, is.finite)), .SDcols = num_cols]
sum(bad_rows)
trainData <- trainData[!bad_rows]

trainData[, BGC := as.factor(BGC)]

# Train model 
BGCmodel <- ranger(
  BGC ~ .,
  data = trainData[, c("BGC", vars_expert), with = FALSE],
  num.trees = 500,
  splitrule =  "extratrees", # way faster than gini, not likely a big performance difference, but we should test this. 
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 14,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  keep.inbag = FALSE, 
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, paste0(dir, "BGCmodel_WNA_V2.2.rds"))

## -------------------------------------------------
## V4.2 RF model - v4 sample with asymptotic reduction in non-BC units

#read in training sample generated in the last step
points <- fread(paste0(dir, "points_WNA_v4a.csv"))

# count of training points by region
points_counts <- points[bgcs_info, on = "BGC", nomatch = 0][, .N, by = Source]
points_counts[Source %in% c("AB_", "USA_", "BC_")]

## climate data for all points
clim <- downscale(
  xyz = points,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim)

trainData <- merge(points, clim, by="id")

# remove points that climr does not return full data for 
num_cols <- names(trainData)[sapply(trainData, is.numeric)]
bad_rows <- trainData[, !Reduce(`&`, lapply(.SD, is.finite)), .SDcols = num_cols]
sum(bad_rows)
trainData <- trainData[!bad_rows]

trainData[, BGC := as.factor(BGC)]

# Train model 
BGCmodel <- ranger(
  BGC ~ .,
  data = trainData[, c("BGC", vars_expert), with = FALSE],
  num.trees = 500,
  splitrule =  "extratrees", # way faster than gini, not likely a big performance difference, but we should test this. 
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 14,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  keep.inbag = FALSE, 
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, paste0(dir, "BGCmodel_WNA_V4.2.rds"))

## -------------------------------------------------
## V4.3 RF model - [v4 sample with asymptotic reduction in non-BC units] + [Tuning]

#read in training sample generated in the last step
points <- fread(paste0(dir, "points_WNA_v4a.csv"))

## climate data for all points
clim <- downscale(
  xyz = points,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim)

trainData <- merge(points, clim, by="id")

# remove points that climr does not return full data for 
num_cols <- names(trainData)[sapply(trainData, is.numeric)]
bad_rows <- trainData[, !Reduce(`&`, lapply(.SD, is.finite)), .SDcols = num_cols]
sum(bad_rows)
trainData <- trainData[!bad_rows]

trainData[, BGC := as.factor(BGC)]

# Train model 
BGCmodel <- ranger(
  BGC ~ .,
  data = trainData[, c("BGC", vars_expert), with = FALSE],
  num.trees = 1000,
  splitrule =  "extratrees", # way faster than gini, not likely a big performance difference, but we should test this. 
  min.node.size = 1, # Default is 1. 
  num.random.splits = 1,
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 14,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  replace = FALSE,
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, paste0(dir, "BGCmodel_WNA_V4.3.rds"))

## -------------------------------------------------
## V4.4 RF model - [v4 sample with asymptotic reduction in non-BC units] + [Tuning based on sensitivity analyses]

#read in training sample generated in the last step
points <- fread(paste0(dir, "points_WNA_v4a.csv"))

## climate data for all points
clim <- downscale(
  xyz = points,
  which_refmap = "refmap_climr",
  return_refperiod = TRUE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim)

trainData <- merge(points, clim, by="id")

# remove points that climr does not return full data for 
num_cols <- names(trainData)[sapply(trainData, is.numeric)]
bad_rows <- trainData[, !Reduce(`&`, lapply(.SD, is.finite)), .SDcols = num_cols]
sum(bad_rows)
trainData <- trainData[!bad_rows]

trainData[, BGC := as.factor(BGC)]

# Train model 
BGCmodel <- ranger(
  BGC ~ .,
  data = trainData[, c("BGC", vars_expert), with = FALSE],
  num.trees = 1000,
  splitrule =  "extratrees", # way faster than gini, not likely a big performance difference, but we should test this. 
  min.node.size = 1, # Default is 1. 
  num.random.splits = 1,
  # importance = 'permutation',
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 14,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  replace = FALSE,
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, paste0(dir, "BGCmodel_WNA_V4.4.rds"))


## -------------------------------------------------
## -------------------------------------------------
## STEP 3 - Results
## -------------------------------------------------
## -------------------------------------------------


library(climr)
library(ccissr)
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

## -------------------------------------------------
## alternative RF models

# dir <- "//objectstore2.nrs.bcgov/ffec/BGC_models/" 
dir <- "C:/Users/CMAHONY/Data/BGC_models/" #local copy, for speed

BGCmodel_v4.2 <- readRDS(paste0(dir, "BGCmodel_WNA_v4.2.rds"))
BGCmodel_V4.2gini <- readRDS(paste0(dir, "BGCmodel_WNA_V4.2gini.rds")) #Kiri trained this on Thufir using the Gini split rule. 
BGCmodel_v4.3 <- readRDS(paste0(dir, "BGCmodel_WNA_v4.3.rds")) # v4a sample with alternative model tunings for potentially improved prediction performance (less overfitting)
BGCmodel_v4.4 <- readRDS(paste0(dir, "BGCmodel_WNA_v4.4.rds")) # v4a sample with model tunings based on sensitivity analysis


# data.table to store results of sensitivity analyses
results.error <- data.table(
  model       = c("v4.2", "V4.2gini", "V4.3", "V4.4"),
  Bamfield    = NA_real_,
  Kamloops    = NA_real_,
  Pemberton   = NA_real_,
  Smithers   = NA_real_,
  BC   = NA_real_
)


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
    dem <- aggregate(dem, fact=5)
    dem <- mask(dem, bdy) 
    dem <- crop(dem, ext(bdy) )
    # plot(dem)
    X <- dem # template raster for testing
    values(X) <- NA
    
    bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_13Apr2026_raster.tif")
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
    # plot(dem)
    X <- dem # template raster for testing
    values(X) <- NA
    
    # study area bgcs
    bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_13Apr2026_raster.tif")
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
  clim_ref[!is.finite(CMD.total), CMD.total := 0]
  
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
  clim_proj[!is.finite(CMD.total), CMD.total := 0]
  
  ## -------------------------------------------------
  ## predictions
  
  ## v4.2 model predictions for reference period
  preds_ref_v4.2_vec <- predict(BGCmodel_v4.2, data = clim_ref)$prediction
  preds_ref_v4.2 <- X
  preds_ref_v4.2[points$id] <- factor(preds_ref_v4.2_vec, levels = subzones_colours_ref$BGC)
  preds_ref_v4.2 <- project(preds_ref_v4.2, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
  
  ## v4.2 model predictions for future period
  preds_proj_v4.2_vec <- predict(BGCmodel_v4.2, data = clim_proj)$prediction
  preds_proj_v4.2 <- X
  preds_proj_v4.2[points$id] <- factor(preds_proj_v4.2_vec, levels = subzones_colours_ref$BGC)
  preds_proj_v4.2 <- project(preds_proj_v4.2, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries.

  ## V4.2gini model predictions for reference period
  preds_ref_V4.2gini_vec <- predict(BGCmodel_V4.2gini, data = clim_ref)$prediction
  preds_ref_V4.2gini <- X
  preds_ref_V4.2gini[points$id] <- factor(preds_ref_V4.2gini_vec, levels = subzones_colours_ref$BGC)
  preds_ref_V4.2gini <- project(preds_ref_V4.2gini, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries.

  ## V4.2gini model predictions for future period
  preds_proj_V4.2gini_vec <- predict(BGCmodel_V4.2gini, data = clim_proj)$prediction
  preds_proj_V4.2gini <- X
  preds_proj_V4.2gini[points$id] <- factor(preds_proj_V4.2gini_vec, levels = subzones_colours_ref$BGC)
  preds_proj_V4.2gini <- project(preds_proj_V4.2gini, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries.

  ## v4.3 model predictions for reference period
  preds_ref_v4.3_vec <- predict(BGCmodel_v4.3, data = clim_ref)$prediction
  preds_ref_v4.3 <- X
  preds_ref_v4.3[points$id] <- factor(preds_ref_v4.3_vec, levels = subzones_colours_ref$BGC)
  preds_ref_v4.3 <- project(preds_ref_v4.3, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
  
  ## v4.3 model predictions for future period
  preds_proj_v4.3_vec <- predict(BGCmodel_v4.3, data = clim_proj)$prediction
  preds_proj_v4.3 <- X
  preds_proj_v4.3[points$id] <- factor(preds_proj_v4.3_vec, levels = subzones_colours_ref$BGC)
  preds_proj_v4.3 <- project(preds_proj_v4.3, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries.
  
  ## v4.4 model predictions for reference period
  preds_ref_v4.4_vec <- predict(BGCmodel_v4.4, data = clim_ref)$prediction
  preds_ref_v4.4 <- X
  preds_ref_v4.4[points$id] <- factor(preds_ref_v4.4_vec, levels = subzones_colours_ref$BGC)
  preds_ref_v4.4 <- project(preds_ref_v4.4, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
  
  ## v4.4 model predictions for future period
  preds_proj_v4.4_vec <- predict(BGCmodel_v4.4, data = clim_proj)$prediction
  preds_proj_v4.4 <- X
  preds_proj_v4.4[points$id] <- factor(preds_proj_v4.4_vec, levels = subzones_colours_ref$BGC)
  preds_proj_v4.4 <- project(preds_proj_v4.4, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries.
  
  
  ## -------------------------------------------------
  ## leaflet map
  
  map <- leaflet() %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%  # Add satellite imagery
    addTiles(group = "OSM Basemap") %>%  # Keep OpenStreetMap as an option
    addRasterImage(preds_ref_v4.2, colors = color_pal, opacity = 1, group = "v4.2: baseline") %>%
    addRasterImage(preds_proj_v4.2, colors = color_pal, opacity = 1, group = "v4.2: future") %>%
    addRasterImage(preds_ref_V4.2gini, colors = color_pal, opacity = 1, group = "V4.2gini: baseline") %>%
    addRasterImage(preds_proj_V4.2gini, colors = color_pal, opacity = 1, group = "V4.2gini: future") %>%
    # addRasterImage(preds_ref_v4.3, colors = color_pal, opacity = 1, group = "v4.3: baseline") %>%
    # addRasterImage(preds_proj_v4.3, colors = color_pal, opacity = 1, group = "v4.3: future") %>%
    # addRasterImage(preds_ref_v4.4, colors = color_pal, opacity = 1, group = "v4.4: baseline") %>%
    # addRasterImage(preds_proj_v4.4, colors = color_pal, opacity = 1, group = "v4.4: future") %>%
    addRasterImage(bgcs_3857, colors = color_pal, opacity = 1, group = "BGC") %>%
    addLayersControl(
      baseGroups = c("Satellite", "OSM Basemap"),  # Base layer switcher
      overlayGroups = c("v4.2: baseline", "v4.2: future",
                        "V4.2gini: baseline", "V4.2gini: future",
                        # "v4.3: baseline", "v4.3: future",
                        # "v4.4: baseline", "v4.4: future",
                        "BGC"),
      options = layersControlOptions(collapsed = FALSE)
    )
  print(map)

  ## -------------------------------------------------
  ## Compute classification error
  
  results.error[1, which(names(results.error)==studyname)] <- mean(values(bgcs_3857, mat = FALSE) != values(preds_ref_v4.2, mat = FALSE), na.rm=T)
  results.error[2, which(names(results.error)==studyname)] <- mean(values(bgcs_3857, mat = FALSE) != values(preds_ref_V4.2gini, mat = FALSE), na.rm=T)
  results.error[3, which(names(results.error)==studyname)] <- mean(values(bgcs_3857, mat = FALSE) != values(preds_ref_v4.3, mat = FALSE), na.rm=T)
  results.error[4, which(names(results.error)==studyname)] <- mean(values(bgcs_3857, mat = FALSE) != values(preds_ref_v4.4, mat = FALSE), na.rm=T)
  
  print(studyname)
}
write.csv(results.error, paste0(dir, "results.error.csv"), row.names = FALSE)


# results.error[, BC := NA_real_]

# Plot of error results
results.error <- fread(paste0(dir, "results.error.csv"))
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
