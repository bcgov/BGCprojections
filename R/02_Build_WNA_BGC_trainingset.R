## This script must be sourced

## FIT AND EVALUATE A BGC MODEL

dPath <- unlist(options("reproducible.destinationPath"))

## this is for testing, otherwise make NULL
studyArea <- vect(ext(c(-125, -112, 43, 55)), crs = "EPSG:4326")

dwnldFromObjSto(prefix = "~/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022",
                bucket = "gmrtde",
                path = dPath)
bgcs <- vect(file.path(dPath, "WNA_BGC_v12_5Apr2022.gpkg"))  ## for poly validation if need be
.gc()
studyArea <- project(studyArea, y = crs(bgcs, proj = TRUE))

cacheExtra <- list(studyArea, table(bgcs$BGC))
bgcs <- Cache(postProcessTerra,
              from = bgcs,
              cropTo = if (is.null(studyArea)) NA else studyArea,
              projectTo = NA,
              maskTo = NA,
              .cacheExtra = cacheExtra,
              userTags = "bgcs", 
              omitArgs = c("from", "userTags", "cropTo", "projectTo", "maskTo"))
.gc()


dwnldFromObjSto(prefix = "~/DEM/DEM_NorAm/NA_Elevation/data/northamerica",
                bucket = "gmrtde",
                path = dPath)
elev <- rast(file.path(dPath, "northamerica_elevation_cec_2023.tif"))

cacheExtra <- list(summary(elev), table(bgcs$BGC))
elev <- Cache(postProcessTerra,
              from = elev,
              cropTo = bgcs,
              projectTo = bgcs,
              maskTo = NA,
              .cacheExtra = cacheExtra,
              userTags = "elev", 
              omitArgs = c("from", "userTags", "cropTo", "projectTo", "maskTo"))
.gc()

trainCoords <- Cache(makePointCoords,
                      bgc_poly = bgcs,
                      elev = elev,
                      .cacheExtra = cacheExtra,
                      userTags = "coords",
                      omitArgs = c("userTags", "bgc_poly", "elev"))
trainCoords <- trainCoords[!is.na(elev),]

climrVars <- c("DD5", "DD_0_at", "DD_0_wt", "PPT05", "PPT06", "PPT07", "PPT08",
                 "PPT09", "CMD", "PPT_at", "PPT_wt", "CMD07", "SHM", "AHM", "NFFD", "PAS", "CMI")
## add BGCs
trainCoords <- Cache(terra::intersect,  ## much faster than extract.
                      x = vect(trainCoords, geom = c("lon", "lat"), crs = crs(bgcs)), 
                      y = bgcs,
                      .cacheExtra = list(summary(trainCoords), cacheExtra[[2]]),
                      userTags = "coords_bgc",
                      omitArgs = c("userTags", "x", "y"))

## project to climr-compatible projection "EPSG:4326"
trainCoords <- Cache(project,
                      x = trainCoords, 
                      y = crs("EPSG:4326", proj = TRUE),
                      .cacheExtra = summary(trainCoords),
                      userTags = "coords_proj",
                      omitArgs = c("userTags", "x"))

## back to DT
trainCoords <- as.data.table(trainCoords, geom = "XY") |>
  setnames(old = c("x", "y"), new = c("lon", "lat"))
trainCoords[, id := seq_along(id)]  ## re-do to have contiguous IDs

climData <- Cache(
  getClimate,
  coords = trainCoords, 
  which_normal = "normal_composite", 
  return_normal = TRUE, 
  vars = climrVars, 
  cache = TRUE,
  .cacheExtra = list(cacheExtra[[2]], summary(trainCoords)),
  userTags = "climData",
  omitArgs = c("userTags", "coords", "bgcs"))

## rm no-longer necessary objects
rm(bgcs)
.gc()

setDT(climData)   ## this shouldn't be necessary, submit issue/reprex to reproducible.
addVars(climData)

climPredictors <- c("DD5", "DD_delayed", "PPT_MJ", "PPT_JAS", 
              "CMD.total", "CMI", "CMDMax", "SHM", "AHM", "NFFD", "PAS")

trainData <- climData[, .SD, .SDcols = c("BGC", "id", climPredictors)]
trainData <- trainData[complete.cases(trainData)]

## subset coords objects to ids with data
trainCoords <- trainCoords[trainData[, .(id)], on = "id", nomatch = 0L]

# BGC_counts <- trainData[, .(Num = .N), by = .(BGC)]   ## for inspection

## clean trainind dataset - remove bad BGCs, outliers and BGCs with low sample sizes.
badBGCs <- c("BWBSvk", "ICHmc1a", "MHun", "SBSun", "ESSFun", "SWBvk", "MSdm3",
             "ESSFdc3", "IDFdxx_WY", "MSabS", "FGff", "JPWmk_WY" ) #, "ESSFab""CWHws2", "CWHwm", "CWHms1" , 
if (!is.na(badBGCs)) {
  trainData <- trainData[!BGC %in% badBGCs,]
}

trainData <- removeOutlier(as.data.frame(trainData), alpha = .025, vars = climPredictors) |>
  Cache()

trainData <- rmLowSampleBGCs(trainData) |>
  Cache()

## subsampling
dataBalance_recipe <- recipe(BGC ~ ., data =  trainData) |>
  step_downsample(BGC, under_ratio = 90) |>  ## subsamples "oversampled" BGCs
  prep()

## extract data.table
trainData <- dataBalance_recipe |>
  juice() |>
  as.data.table()

# trainData[,.(Num = .N), by = BGC]   ## for inspection

trainData[, BGC := as.factor(BGC)]

## FIT MODEL WITH FULL DATA SET OR CROSS-VALIDATION
# cols <- c("BGC", climPredictors) ## no longer used?

## convert to sf for spatial cv below
trainData_SF <- trainCoords[, .(id, lon, lat)][trainData, on = "id"]
trainData_SF <- sf::st_as_sf(trainData_SF, coords = c("lon", "lat"))
sf::st_crs(trainData_SF) <- crs(elev, proj = TRUE)

rm(elev);
.gc()

## make a modelling task
trainData_SF$BGC <- as.factor(trainData_SF$BGC)
tsk_bgc <- as_task_classif_st(trainData_SF, target = "BGC")

## choose model type and eval metric
# lrn_rf_resp <- lrn("classif.ranger", 
#                    predict_type = "response",
#                    num.trees = 501,
#                    splitrule =  "extratrees",
#                    mtry = 4,
#                    min.node.size = 2,
#                    importance = "permutation",
#                    write.forest = TRUE)

## we'll also fit a probability type model for current/normal period predictions
lrn_rf_prob <- lrn("classif.ranger", 
                   predict_type = "prob",
                   num.trees = 501,
                   splitrule =  "extratrees",
                   mtry = 4,
                   min.node.size = 2,
                   importance = "permutation",
                   write.forest = TRUE)

measure_acc <- msrs(c("classif.acc", "classif.ce", "oob_error"))

## fit full models --------------------------------------------
# .gc()
# lrn_rf_resp <- trainModel(lrn_rf_resp, tsk_bgc) |>
#   Cache(userTags = "lrn_rf_resp",
#         .cacheExtra = summary(trainData_SF),
#         omitArgs = c("task"))  ## for some reason task changes when it should be the same
# .gc()
# 
# lrn_rf_resp$model

## running out of mem when using full WNA
.gc()
# lrn_rf_prob <- trainModel(lrn_rf_prob, tsk_bgc) |>
#   Cache(userTags = "lrn_rf_prob",
#         .cacheExtra = list(mlr3misc::calculate_hash(lrn_rf_prob),  ## these hashes are not stable. need to find another way of caching.
#                            mlr3misc::calculate_hash(trainData_SF)),
#         omitArgs = c("learner", "task"), showSimilar = TRUE)  ## for some reason task changrm(es when it should be the same
lrn_rf_prob <- loadFromCache(cacheId = "7c1eea2779fa66f0")
.gc()
lrn_rf_prob$model
## this plot is useless, but code is kept here for fut reference
# autoplot(lrn_rf_resp$predict(tsk_bgc)) +
#   theme(legend.position = "none")

## evaluate model with CV ------------------------------------
## define and apply cv strategy
# folds <- 10
# cv_strategy <- rsmp("repeated_spcv_coords", folds = folds, repeats = 1)
# 
# future::plan("multisession", 
#              workers = ifelse(folds <= future::availableCores(), 
#                               folds, 
#                               future::availableCores() - 2))
# RF_cv_resp <- mlr3::resample(tsk_bgc, lrn_rf_resp, cv_strategy, store_models = TRUE) |>
#   Cache()
# 
# RF_cv_prob <- mlr3::resample(tsk_bgc, lrn_rf_prob, cv_strategy, store_models = TRUE) |>
#   Cache()
# future:::ClusterRegistry("stop")
# 
# ## eval metrics
# ## TODO: find a pretty way to export these
# RF_cv_prob$score(measure_acc)  ## eval of each fold
# RF_cv_prob$aggregate(measure_acc)  ## aggregated scores
# RF_cv_prob$prediction()$score(msrs(c("classif.acc", "classif.ce"),
#                                    average = "micro"))  ## pool predictions across resampling iterations into one Prediction object and then compute the measure on this directly
# ## confusion matrix
# RF_cv_prob$prediction()$confusion
# 
# ## TODO: save/export models and evaluation metrics