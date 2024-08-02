defineModule(sim, list(
  name = "BGCWNAmodel_trainAndProj",
  description = paste("This module integrates thedata preparation, model fitting and prediction/projection",
                      "steps of the biogeoclimatic (BGC) zone projections worflow. By default, it prepares the training",
                      "data as a table of 'current' BGC zone, elevation and climate data obtained for the desired",
                      "study area (must be in Western Canada + Western US), at points sampled by BGC using spatially and",
                      "climatically balanced sampling. Trainigng cliamte data refers to the reference climate period",
                      "1961-1990. The 'prediction' data uses climate values obtained from a selection of General Circulation Models",
                      " (GCMs), Shared Socioeconomic Pathway scenarios (SSPs), model runs, historic and future periods",
                      "and extracted at the centre points of an hexagonal grid. The historic period refers to 2001-2020 and",
                      "is only used to ... <TO COMPLETE>"),
  keywords = c("biogeoclimatic zones", "climate", "data preparation", "CCISS"),
  authors = structure(list(
    list(given = "Ceres", family = "Barros", role = c("aut"), email = "ceres.barros@nrcan-rncan.gc.ca", comment = NULL),
    list(given = "Kiri", family = "Daust", role = c("aut"), email = "kiri.daust@gov.bc.ca", comment = NULL),
    list(given = "Will", family = "MacKenzie", role = c("aut"), email = "will.mackenzie@gov.bc.ca", comment = NULL),
    list(given = "Colin", family = "Mahony", role = c("aut"), email = "colin.mahony@gov.bc.ca", comment = NULL)
  ), class = "person"),
  childModules = character(0),
  version = list(BGCWNAmodel_trainAndProj = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "BGCWNAmodel_trainAndProj.Rmd"),
  loadOrder = list(),
  reqdPkgs = list(
    ## SpaDES/caching-related:
    "reproducible",
    "SpaDES.core (>= 2.1.5)", 
    ## other:
    "aws.s3",
    "bcgov/climr@devl (HEAD)",
    "bcgov/ccissr@development (HEAD)",
    "crayon",
    "data.table",
    "foreach",
    "gdalUtilities",
    "ggplot2",
    "mlr3learners",
    "mlr3spatial",
    "mlr3spatiotempcv",
    "mlr3viz",
    "ranger",
    "sf",
    "smotefamily", 
    "terra",
    "themis",
    "tidymodels"
  ),
  parameters = bindrows(
    defineParameter("accuracyMetrics", "character", c("classif.acc", "classif.ce", "oob_error"),
                    NA, NA, 
                    paste("A vector of accuracy metrics to use when evaluating fitted models.",
                          "Passed to `mlr3::msrs()`. See `as.data.table(mlr_measures)` for a",
                          "list of measures and `?mlr3::mlr_measures` for details. Note that",
                          "the accuracy metrics chosen must be appropriate for the models listed in",
                          "`sim$modelsToTrain`")),
    defineParameter("assertions", "logical", TRUE, NA, NA, 
                    paste("Turn on/off assertions -- i.e. checks ",
                          "for object conformity and expected code behaviour.")),
    defineParameter("badBGCs", "character", 
                    c("BWBSvk", "ICHmc1a", "MHun", "SBSun", "ESSFun", "SWBvk", "MSdm3",
                      "ESSFdc3", "IDFdxx_WY", "MSabS", "FGff", "JPWmk_WY"), 
                    NA, NA,
                    paste("OPTIONAL. BCG zones that to remove from the model fitting",
                          "procedure and subsequent analyses. By default, these are zones with very low",
                          "sample sizes or that were shown to behave strangely in a predictive",
                          "setting (e.g., overly predicted due to a broad climate space) during",
                          "exploratory analyses. Set to `NA_character` if no BGCs are to be excluded.")),
    defineParameter("balance", "logical", TRUE, NA, NA,
                    paste("Should `sim$trainData` be balanced? See `sim$trainData` for details.")),
    defineParameter("climrVars", "character", 
                    c("DD5", "DDsub0_at", "DDsub0_wt", "PPT_05", "PPT_06", "PPT_07", 
                      "PPT_08", "PPT_09", "CMD", "PPT_at", "PPT_wt", "CMD_07", 
                      "SHM", "AHM", "NFFD", "PAS", "CMI"), 
                    NA, NA,
                    paste("Climate variables to obtain from `climr::downscale`.",
                          "May include variables listed in `P(sim)$climPredictors` or variables",
                          "necessary to produce them (using `addVars()`). See `climr::list_vars()`." )),
    defineParameter("climPredictors", "character", 
                    c("DD5", "DD_delayed", "PPT_MJ", "PPT_JAS", "CMD.total", 
                      "CMI", "CMDMax", "SHM", "AHM", "NFFD", "PAS"), 
                    NA, NA,
                    paste("Climate covariares used to fit BGC zone ~ climate model. May include variables",
                          "and others derived from them using `addVars()`. See `climr::list_vars()`",
                          "and `?addVars()` for options.")),
    defineParameter("CVfolds", "numeric", 10, 1, Inf,
                    paste("Number of folds to use for model k-fold cross-validation.",
                          "Passed to `mlr3::rsmp('repeated_spcv_coords', folds)`")),
    defineParameter("GCMs", "character", climr::list_gcms(), NA, NA,
                    paste("List of General Circulation Models (GCMs) to obtain climate projections for",
                          "`sim$predData`.")),
    defineParameter("GCMperiods", "character", climr::list_gcm_periods(), NA, NA,
                    paste("List of future periods to obtain climate projections for `sim$predData`.")),
    defineParameter("GCMruns", "numeric", 3, 1, Inf,
                    paste("Number of General Circulation Models (GCMs) runs to obtain climate projections for",
                          "`sim$predData`. Note that the ensemble mean across runs will not be used.")),
    defineParameter("gridSize", "numeric", 200L, 1, Inf,
                    paste("Distance in m between points of the regular grid used to extract",
                          "BGC, elevation (used for climate downscaling) and climate data to",
                          "fit the BGC zones predictive model.")),
    defineParameter("SSPs", "character", climr::list_ssps(), NA, NA,
                    paste("List of Shared Socioeconomic Pathway scenarios (SSPs) to obtain climate projections for",
                          "`sim$predData`.")),
    defineParameter("nthread", "numeric", 1, 1, Inf,
                    paste("Passed to `climr::downscale(..., nthread)` to parallelise climate downscaling operations")),
    defineParameter("nthread_mlr3", "list", list("lrn_rf_resp" = 3, "lrn_rf_prob" = 1), NA, NA,
                    paste("Named list of number of cores to use to parallellise model training",
                          "(each element is passed to `mlr3::set_threads(Learner, n)` before running `<Learner>$train()`)",
                          "and model evaluation (passed to `future` before running `mlr3::resample`).",
                          "`names(P(sim)$nthread_mlr3)` must be the same as `names(sim$modelsToTrain)`",
                          "so that the number of cores used can differ by model algorithm. By default, 3 cores are",
                          "used for `sim$modelsToTrain$lrn_rf_resp` and 1 for `sim$modelsToTrain$lrn_rf_prob`.")),
    ## default global params.
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", c(".inputObjects"), NA, NA,
                    "Should caching of events or module be used? We advise not caching",
                    "events other than '.inputObjects' as this module produces a large 'simList'",
                    "that can take long to digest. Internal caching still occurs in each event.")
  ),
  inputObjects = bindrows(
    expectsInput("bgcs", "SpatVector", 
                 paste("A polygon SpatVector of 'current' biogeoclimatic",
                       "zones across Western Canada and the US. This map",
                       "was produced ...TODO. Hosted privately.", 
                       "Please contact developers to obtain access."), 
                 sourceURL = NA),
    expectsInput("districts", "SpatVector", 
                 paste("A polygon SpatVector of forest districts in British",
                       "Columbia. Obtained from ...TODO. Hosted privately.", 
                       "Please contact developers to obtain access."), 
                 sourceURL = NA),
    expectsInput("elev", "SpatRaster", 
                 paste("A DEM SpatRaster. Elevation values are used to dowscale",
                       "climate data using `climr`. Defaults to Western Canada and",
                       "US raster at 250m resolution from ...TODO. Hosted privately.", 
                       "Please contact developers to obtain access."), 
                 sourceURL = NA),
    expectsInput("hexgrid", "SpatVector", 
                 paste("A hexagonal grid covering BC at 400m 'resolution'.", 
                       "Hexagon center points will be used as locations to predict",
                       "BGC zones under novel climate conditions. Hosted privately.", 
                       "Please contact developers to obtain access."), 
                 sourceURL = NA),
    expectsInput("modelsToTrain", "list", 
                 paste("A named list of `mlr3` model learner calls (`mlr3::lrn`).",
                       "List names will be used as output model names. By default,",
                       "two Random Forest models will be fit: one predicting the response",
                       "variable (BGCs), another predicting probabilities of each response",
                       "category. To see these calls, run `simInit` and consult `sim$modelsToTrain`",
                       "or navigate to the `.inputObjects` section of the module `.R` script."), 
                 sourceURL = NA),
    expectsInput("studyArea", "SpatVector", 
                 paste("OPTIONAL. A polygon of a study area to subset `sim$bgcs`",
                       "and, consequently, the training data `sim$trainData` to."), 
                 sourceURL = NA),
    expectsInput("studyArea_predict", "SpatVector", 
                 paste("OPTIONAL. A polygon of a study area in to subset `sim$hexgrid`",
                       "and, consequently, the prediction data `sim$predData` to."), 
                 sourceURL = NA),
    expectsInput("validStrategies", "list",
                 paste("Named list of validation strategies passed to models defined in `sim$modelsToTrain`.",
                       "Each entry should que a quoted call to `mlr3::rsmp`, that specifies the strategy",
                       "and parameters used (e.g. `quote(mlr3::rsmp('cv', folds = 10))`). List `names()` must", 
                       "match `names(sim$modelsToTrain)`. By default, repeated coordinate-based k-means clustering",
                       "used for all models (see `?mlr3spatiotempcv::ResamplingRepeatedSpCVCoords`). For other",
                       "strategies and more detail see `as.data.table(mlr_resamplings)` and `?mlr_resamplings`.",
                       sourceURL = NA))
  ),
  outputObjects = bindrows(
    createsOutput("fittedBGCmodels", "list",
                  paste("A named list containing the fitted model (`$BGCmodel`) and validation",
                        "(`$BGCmodel_val`) objects computed by `mlr3` for each learner and validation",
                        "strategy specified by `sim$modelsToTrain` and `sim$validStrategies`, respectively.",
                        "`names(fittedBGCmodels)` will follow `names(sim$modelsToTrain)`.")),
    createsOutput("trainData", "data.table", 
                  paste("A training dataset to fit predictive models of",
                        "biogeoclimatic zones (BGCs) ~ climate relationships. Includes",
                        "a column of 'current' BGCs (response), columns of longitude ('lon') and latitude ('lat'),", 
                        "columns of climate covariates defined by `P(sim)$climPredictors`",
                        "and a 'CRS' attribute pertaining to the projection used in 'lon', 'lat' columns.",
                        "Climate covariate values correspond to the 1961-1990 reference period.",
                        "Rows correspond to point locations sampled",
                        "from a regular grid of `P(sim)$gridSize` spacing.",
                        "If `P(sim)$balance` is `TRUE`, `sim$trainData` is subsampled",
                        "so that 1) outliers within each BGC zone are removed; 2) BGC zones with",
                        "low sample sizes are removed; and 3) oversampled BGC zones are subsampled.",
                        "For 1), outliers are defined as observations whose climatic distance from",
                        "the mean climate (calculated using Mahalanobis distance) is > 0.975 quantile",
                        "of a Chi-square distribution. For 2), see `P(sim)$badBGCs`. For 3), data is",
                        "balanced so that the ratio of the minority-to-majority frequencies is 0.9",
                        "(i.e. majority levels, here BGCs, are downsampled so that they have at most 90% of the",
                        "points of the least occurring level)."))
  )
))


doEvent.BGCWNAmodel_trainAndProj = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- trainingDataEvent(sim)
      
      if (P(sim)$balance) {
        sim <- balanceDataEvent(sim)
      }
      ## schedule future events
      sim <- scheduleEvent(sim, start(sim), "BGCWNAmodel_trainAndProj", "fitModel")
      sim <- scheduleEvent(sim, start(sim), "BGCWNAmodel_trainAndProj", "predict", .normal()+1)
    },
    fitModel = {
      ## TODO: check memopry requirements to fit model to full data
      sim <- fitBGCmodelEvent(sim)
    },
    predictHistorical = {
      ## we may need a separate event when using the probability-type model to project 
      ## 2001-2020 data -- we will also need some sort of parameter that tells us which model to use for that
    },
    predictFuture = {
      ## this needs to change. we will have one event function that 
      ## uses a chain of functions that will do the following steps by chunks of data
      ## 1. break up the prediction point locations into chunks
      ## 2. for each chunk, get/prep data for prediction using climr - many GCMs, runs, periods, ssps
      ## 3. for each chunk predict BGCs for all climate projections
      ## 4. for each chunk calculate ensemble BGCs at the zonal and subzonal levels
      ## 5. push each chunk to PostGRS
      
      sim <- predictionDataEvent(sim)  ## rename to   predictBGCmodelEvent
      
      ## TODO: bring code from R/03_RunPredictHex.R
      
      ## if predicting again - not necessayr, here for pedagocic reasons :)
      # sim <- scheduleEvent(sim, time(sim) + 1, "BGCWNAmodel_trainAndProj", "predict", .normal()+1)
    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

trainingDataEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  cacheTags <- c(currentModule(sim), "function:trainingDataEvent")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  
  ## PREP SPATIAL DATA -------------------------------
  ## check SA and crop if need be
  if (!is.null(sim$studyArea)) {
    CRS <- crs(sim$bgcs, proj = TRUE)
    
    ## July 2024 - NRCan network blocking HTTPS, which causes terra to fail to reproject.
    ## see https://stackoverflow.com/questions/76973035/certgetcertificatechain-trust-error-cert-trust-is-untrusted-root-gdal-error-1
    studyArea <- st_as_sf(sim$studyArea)
    studyArea <- st_transform(studyArea, y = CRS)
    
    sim$studyArea <- vect(studyArea)
    
    ## July 2024 - NRCan network blocking HTTPS, which causes terra to fail to reproject.
    ## see https://stackoverflow.com/questions/76973035/certgetcertificatechain-trust-error-cert-trust-is-untrusted-root-gdal-error-1
    cacheExtra <- list(sim$studyArea, table(sim$bgcs$BGC))
    bgcs <- Cache(postProcessTo,
                  from = st_as_sf(sim$bgcs),
                  cropTo = if (is.null(sim$studyArea)) NA else st_as_sf(sim$studyArea),
                  projectTo = NA,
                  maskTo = NA,
                  .cacheExtra = cacheExtra,
                  userTags = c(cacheTags, "bgcs"),
                  omitArgs = c("from", "userTags", "cropTo", "projectTo", "maskTo"))
    sim$bgcs <- vect(bgcs)
    .gc()
  } else {
    message(magenta("studyArea not provided. No cropping of `sim$bgcs` will be done."))
  }
  
  message(magenta("Preparing training data..."))
  
  cacheExtra <- list(summary(sim$elev), table(sim$bgcs$BGC))
  elev <- Cache(postProcessTerra,
                from = sim$elev,
                cropTo = sim$bgcs,
                projectTo = sim$bgcs,
                maskTo = NA,
                .cacheExtra = cacheExtra,
                userTags = c(cacheTags, "elev_train"), 
                omitArgs = c("from", "userTags", "cropTo", "projectTo", "maskTo"))
  .gc()
  
  ## MAKE TRAINING DATA -------------------------------
  ## generate points at regularly spaced grid and extract elevation using interpolation
  trainCoords <- Cache(makePointCoords,
                       bgc_poly = sim$bgcs,
                       elev = elev,
                       gridSize = P(sim)$gridSize,
                       .cacheExtra = cacheExtra,
                       userTags = c(cacheTags, "trainCoords"),
                       omitArgs = c("userTags", "bgc_poly", "elev"))
  trainCoords <- trainCoords[!is.na(elev),]
  
  rm(elev)
  .gc()
  
  ## checks
  if (P(sim)$assertions) {
    test <- vect(trainCoords, geom = c("lon", "lat"), crs = crs(sim$bgcs)) |>
      geom()
    if (nrow(test[!complete.cases(test),])) {
      stop("NAs found in lat/long values for elevation. Please debug 'trainingDataEvent' function")
    }
  }
  
  ## add BGCs
  trainCoords <- Cache(terra::intersect,  ## much faster than extract.
                       x = vect(trainCoords, geom = c("lon", "lat"), crs = crs(sim$bgcs)), 
                       y = sim$bgcs,
                       .cacheExtra = list(summary(trainCoords), cacheExtra[[2]]),
                       userTags = c(cacheTags, "trainCoords_bgc"),
                       omitArgs = c("userTags", "x", "y"))
  
  if (P(sim)$assertions) {
    test <- geom(trainCoords)
    
    if (nrow(test[!complete.cases(test),])) {
      stop("NAs found in lat/long or BGC zones. Please debug 'trainingDataEvent' function")
    }
  }
  
  
  ## project to climr-compatible projection "EPSG:4326"
  trainCoords <- Cache(projectTo,
                       from = trainCoords, 
                       projectTo = crs("EPSG:4326", proj = TRUE),
                       .cacheExtra = summary(trainCoords),
                       userTags = c(cacheTags, "trainCoords_proj"),
                       omitArgs = c("userTags", "x"))
  
  if (P(sim)$assertions) {
    test <- geom(trainCoords)
    
    if (nrow(test[!complete.cases(test),])) {
      stop("NAs generated after reprojecting. Please debug 'trainingDataEvent' function")
    }
  }
  
  ## back to DT
  trainCoords <- as.data.table(trainCoords, geom = "XY") |>
    setnames(old = c("x", "y"), new = c("lon", "lat"))
  trainCoords[, id := seq_along(id)]  ## re-do to have contiguous IDs
  
  browser() ## overwrite/clear cache and remove loadCache
  ## get climate data -- takes almost 1h
  # climData <- Cache(
  #   getClimate,
  #   coords = trainCoords, 
  #   which_refmap = "auto", 
  #   return_refperiod = TRUE, ## fitting on 1961-1990 period 
  #   vars = P(sim)$climrVars, 
  #   byCombo = TRUE,
  #   cache = TRUE,
  #   useCache = TRUE,
  #   .cacheExtra = list(summary(trainCoords)),
  #   userTags = c(cacheTags, "climData"),
  #   omitArgs = c("userTags", "trainCoords", "bgcs"))
  climData <- loadFromCache(cacheId = "3842396fdc799084")
  setDT(climData)  ## something with the cache is screwing up the data.table
  
  ## add more climate variables
  addVars(climData)
  
  ## subset data to variables of interest
  trainData <- climData[, .SD, .SDcols = c("BGC", "id", P(sim)$climPredictors)]
  trainData <- trainData[complete.cases(trainData)]
  
  ## add back original lon, lat, elev and a CRS attribute
  trainData <- trainCoords[, .(id, elev, lon, lat)][trainData, on = "id", nomatch = 0L]
  attr(trainData, "CRS") <- crs("EPSG:4326", proj = TRUE)
  
  ## export to sim
  sim$trainData <- trainData
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

balanceDataEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  cacheTags <- c(currentModule(sim), "function:balanceDataEvent")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  
  message(magenta("Balancing training data..."))
  trainDataCRS <- attr(sim$trainData, "CRS")
  
  message(magenta("Number of data points at:"))
  message(nrow(sim$trainData))
  
  if (!all(is.na(P(sim)$badBGCs))) {
    message(magenta("Excluding the following BGCs from further analyses:"))
    message(magenta("  ", paste(P(sim)$badBGCs, collapse = ", ")))
    sim$trainData <- sim$trainData[!BGC %in% P(sim)$badBGCs,]
  }
  
  sim$trainData <- removeOutlier(as.data.frame(sim$trainData), alpha = .025, vars = climPredictors) |>
    Cache(userTags = c(cacheTags, "outliers"))
  
  sim$trainData <- rmLowSampleBGCs(sim$trainData) |>
    Cache(userTags = c(cacheTags, "rmLowSampleBGCs"))
  
  dataBalance_recipe <- recipe(BGC ~ ., data =  sim$trainData) |>
    step_downsample(BGC, under_ratio = 90) |>  ## subsamples "oversampled" BGCs
    prep()
  
  ## extract data.table
  sim$trainData <- dataBalance_recipe |>
    juice() |>
    as.data.table()
  
  # sim$trainData[,.(Num = .N), by = BGC]   ## for inspection
  
  sim$trainData[, BGC := as.factor(BGC)]
  
  ## add back CRS attribute
  attr(sim$trainData, "CRS") <- trainDataCRS
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

fitBGCmodelEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  cacheTags <- c(currentModule(sim), "function:fitBGCmodelEvent")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  
  ## checks
  if (P(sim)$assertions) {
    assertNames(sim$modelsToTrain, P(sim)$nthread_mlr3)
    assertNames(sim$modelsToTrain, sim$validStrategies)
  }
  
  ## make sf object
  trainData_SF <- sf::st_as_sf(sim$trainData, coords = c("lon", "lat"))
  sf::st_crs(trainData_SF) <- attr(sim$trainData, "CRS")
  
  ## make a modelling task
  trainData_SF$BGC <- as.factor(trainData_SF$BGC)
  tsk_bgc <- as_task_classif_st(trainData_SF, target = "BGC")
  
  ## create model "learners" -- choose model type 
  ## a response-type model
  modelList <- sapply(sim$modelsToTrain, eval, USE.NAMES = TRUE, simplify = FALSE)
  
  ## eval metrics
  ## TODO: PASS A LIST
  measure_acc <- msrs(P(sim)$accuracyMetrics)
  
  ## set number of cores to parellelise over
  nthreads <- P(sim)$nthread_mlr3[names(modelList)] ## make sure order is the same
  Map(f = set_threads,
      x = modelList,
      n = nthreads)
  
  ## train models
  modelList <- Map(
    lrner = modelList,
    modname = names(modelList),
    f = trainLearners,
    MoreArgs = list("task" = tsk_bgc, "cacheTags" = cacheTags)
  )
  
  ## evaluate model with CV ------------------------------------
  ## define and apply cv strategy
  valStrategiesList <- sapply(sim$validStrategies, eval)
  
  fittedBGCmodels <- Map(
    lrner = modelList,
    strategy = valStrategiesList,
    modname = names(modelList),
    nthreads = nthreads,
    f = evalLearners,
    MoreArgs = list(
      "task" = tsk_bgc,
      "cacheTags" = cacheTags
    )
  )
  
  ## inspect:
  ## TODO: find a pretty way to export these
  # cvList[[1]]$score(measure_acc)  ## eval of each fold
  # cvList[[1]]$aggregate(measure_acc)  ## aggregated scores
  # cvList[[1]]$prediction()$score(msrs(c("classif.acc", "classif.ce"),
  #                                    average = "micro"))  ## pool predictions across resampling iterations into one Prediction object and then compute the measure on this directly
  # ## confusion matrix
  # cvList[[1]]$prediction()$confusion
  
  ## export to sim
  ## export a list of lists. for each model name put model and validation object.
  sim$fittedBGCmodels <- fittedBGCmodels
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

predictionDataEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  ## TODO: move into predictBGCmodel
  cacheTags <- c(currentModule(sim), "function:predictionDataEvent")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  
  if (!is.null(sim$studyArea_predict)) {
    CRS <- crs(sim$bgcs, proj = TRUE)
    
    ## July 2024 - NRCan network blocking HTTPS, which causes terra to fail to reproject.
    ## see https://stackoverflow.com/questions/76973035/certgetcertificatechain-trust-error-cert-trust-is-untrusted-root-gdal-error-1
    studyArea_predict <- st_as_sf(sim$studyArea_predict)
    studyArea_predict <- st_transform(studyArea_predict, y = CRS)
    
    sim$studyArea_predict <- vect(studyArea_predict)
    sim$hexgrid <- GDALcrop("HexGrid400m_Sept2021.gpkg", "HexGrid400m_Sept2021_crop.gpkg",
                            dPath, studyArea = sim$studyArea_predict) |>
      Cache(userTags = c(cacheTags, "hexgrid"))
    .gc()
  } else {
    message(magenta("studyArea_predict not provided. No cropping of hexgrid will be done."))
  }
  
  ## get centroids
  cacheExtra <- c(summary(sim$hexgrid), crs(sim$hexgrid, proj = TRUE), ext(sim$hexgrid))
  hexcentroids <- Cache(centroids,
                        x = sim$hexgrid[1:1000],
                        .cacheExtra = cacheExtra,
                        userTags = c(cacheTags, "hexcentroids"), 
                        omitArgs = c("x", "userTags"))
  .gc()
  
  ## add districts
  cacheExtra <- c(summary(hexcentroids), crs(hexcentroids, proj = TRUE), ext(hexcentroids),
                  summary(sim$districts), crs(sim$districts, proj = TRUE), ext(sim$districts))
  hexcentroids <- Cache(terra::intersect,
                        x = hexcentroids, 
                        y = sim$districts, 
                        .cacheExtra = cacheExtra,
                        userTags = c(cacheTags, "hexcentroids_districts"),
                        omitArgs = c("userTags", "x", "y")) 
  
  ## extract elevation with interpolation - crop first
  cacheExtra <- c(summary(sim$elev), crs(sim$elev, proj = TRUE), ext(sim$elev),
                  summary(hexcentroids), crs(hexcentroids, proj = TRUE), ext(hexcentroids))
  elev <- Cache(postProcessTerra,
                from = sim$elev,
                cropTo = hexcentroids,
                projectTo = hexcentroids,
                maskTo = NA,
                .cacheExtra = cacheExtra,
                userTags = c(cacheTags, "elev"),
                omitArgs = c("from", "userTags", "cropTo", "projectTo", "maskTo"))
  .gc()
  
  predCoords <- Cache(terra::extract,
                      x = elev, 
                      y = hexcentroids, 
                      method = "bilinear",
                      cells = TRUE,
                      xy =  TRUE,
                      .cacheExtra = cacheExtra,
                      userTags = c(cacheTags, "hexcoords"),
                      omitArgs = c("userTags", "x", "y")
  ) |>
    as.data.table()
  
  .gc()
  predCoords[, cell := NULL]
  predCoords[, dist_code := hexcentroids$dist_code]
  
  ## change column names to climr standard
  setnames(predCoords, c("id", "elev", "lon", "lat", "dist_code"))
  
  ## exclude NAs
  predCoords <- predCoords[complete.cases(predCoords)]
  
  ## project to climr-compatible projection EPSG:4326
  ## project to EPSG:4326
  points <- vect(predCoords[, .(lon, lat)], crs = crs(elev, proj = TRUE))
  points <- project(points, crs("EPSG:4326", proj = TRUE))
  
  predCoords[, `:=`(lon = geom(points)[, "x"],
                    lat = geom(points)[, "y"])]
  
  ## free RAM
  rm(elev, hexcentroids)
  .gc()
  
  ## get climate data from climr
  ## extract by subsets of points, then write to csv with append = TRUE
  cacheExtra <- c(summary(predCoords), crs(points, proj = TRUE), ext(points))
  projClimData <- Cache(
    getClimate,
    coords = predCoords, 
    which_refmap = "auto", 
    return_refperiod = FALSE, 
    gcms = P(sim)$GCMs,
    ssps = P(sim)$SSPs,
    gcm_periods = P(sim)$GCMperiods,
    max_run = P(sim)$GCMruns,
    vars = P(sim)$climrVars, 
    nthread = P(sim)$nthread,
    cache = TRUE,
    byCombo = TRUE,
    .cacheExtra = cacheExtra,
    userTags = c(cacheTags, "projClimData"),
    omitArgs = c("userTags", "coords"))
  
  setDT(projClimData) ## something with the cache is screwing up the data.table
  
  ## add more climate variables
  addVars(projClimData)
  
  ## subset data to variables of interest
  predData <- projClimData[, .SD, .SDcols = c("id", P(sim)$climPredictors)]
  predData <- predData[complete.cases(predData)]
  
  ## add back original lon lat, and a CRS attribute
  predData <- predCoords[predData, on = "id", nomatch = 0L]
  attr(predData, "CRS") <- crs("EPSG:4326", proj = TRUE)
  
  # sim$predData <- predData   ## not exporting anymore
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

predictBGCmodelEvent <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:predictBGCmodelEvent")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")
  
  # ! ----- EDIT BELOW ----- ! #
  ## BGCs
  if (!suppliedElsewhere("bgcs", sim)) {
    dwnldFromObjSto(prefix = "~/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022",
                    bucket = "gmrtde",
                    path = dPath)
    sim$bgcs <- vect(file.path(dPath, "WNA_BGC_v12_5Apr2022.gpkg"))
    .gc()
  }
  
  ## if available crop now to save MEM
  if (!is.null(sim$studyArea)) {
    CRS <- crs(sim$bgcs, proj = TRUE)
    
    ## July 2024 - NRCan network blocking HTTPS, which causes terra to fail to reproject.
    ## see https://stackoverflow.com/questions/76973035/certgetcertificatechain-trust-error-cert-trust-is-untrusted-root-gdal-error-1
    studyArea <- st_as_sf(sim$studyArea)
    studyArea <- st_transform(studyArea, y = CRS)
    
    sim$studyArea <- vect(studyArea)
    
    ## July 2024 - NRCan network blocking HTTPS, which causes terra to fail to reproject.
    ## see https://stackoverflow.com/questions/76973035/certgetcertificatechain-trust-error-cert-trust-is-untrusted-root-gdal-error-1
    cacheExtra <- list(sim$studyArea, table(sim$bgcs$BGC))
    bgcs <- Cache(postProcessTo,
                  from = st_as_sf(sim$bgcs),
                  cropTo = if (is.null(sim$studyArea)) NA else st_as_sf(sim$studyArea),
                  projectTo = NA,
                  maskTo = NA,
                  .cacheExtra = cacheExtra,
                  userTags = c(cacheTags, "bgcs"),
                  omitArgs = c("from", "userTags", "cropTo", "projectTo", "maskTo"))
    sim$bgcs <- vect(bgcs)
    .gc()
  }
  
  ## DEM
  if (!suppliedElsewhere("elev", sim)) {
    dwnldFromObjSto(prefix = "~/DEM/DEM_NorAm/NA_Elevation/data/northamerica",
                    bucket = "gmrtde",
                    path = dPath)
    sim$elev <- rast(file.path(dPath, "northamerica_elevation_cec_2023.tif"))
  }
  
  ## HEX grid
  if (!suppliedElsewhere("hexgrid", sim)) {
    dwnldFromObjSto(prefix = "~/CCISS_Working/WNA_BGC/HexGrid400m_Sept2021",
                    bucket = "gmrtde",
                    path = dPath)
    
    sim$hexgrid <- vect(file.path(dPath, "HexGrid400m_Sept2021.gpkg"))
  }
  
  ## if available crop now to save MEM
  if (!is.null(sim$studyArea_predict)) {
    CRS <- crs(sim$bgcs, proj = TRUE)
    
    ## July 2024 - NRCan network blocking HTTPS, which causes terra to fail to reproject.
    ## see https://stackoverflow.com/questions/76973035/certgetcertificatechain-trust-error-cert-trust-is-untrusted-root-gdal-error-1
    studyArea_predict <- st_as_sf(sim$studyArea_predict)
    studyArea_predict <- st_transform(studyArea_predict, y = CRS)
    
    sim$studyArea_predict <- vect(studyArea_predict)
    sim$hexgrid <- GDALcrop("HexGrid400m_Sept2021.gpkg", "HexGrid400m_Sept2021_crop.gpkg",
                            dPath, studyArea = sim$studyArea_predict) |>
      Cache(userTags = c(cacheTags, "hexgrid"))
    .gc()
  }
  
  ## Forest districts
  if (!suppliedElsewhere("districts", sim)) {
    dwnldFromObjSto(prefix = "~/CommonTables/Forest_Region_District",
                    bucket = "gmrtde",
                    path = dPath)
    
    ## crop first to save memory  during intersect
    sim$districts <- GDALcrop("Forest_Region_District.gpkg", "Forest_Region_District_crop.gpkg",
                              dPath, studyArea = sim$studyArea_predict) |>
      Cache(userTags = c(cacheTags, "districts"))
    
    sim$districts <- sim$districts["ORG_UNIT"]
    sim$districts$ORG_UNIT <- as.character(sim$districts$ORG_UNIT)
    sim$districts$ORG_UNIT[sim$districts$ORG_UNIT == "DSS"] <- "CAS"
    names(sim$districts) <- "dist_code"
    .gc()
  }
  
  ## Model-related inputs
  if (!suppliedElsewhere("modelsToTrain", sim)) {
    sim$modelsToTrain <- list(
      "lrn_rf_resp" = quote(lrn("classif.ranger",
                                predict_type = "response",
                                num.trees = 501,
                                splitrule =  "extratrees",
                                mtry = 4,
                                min.node.size = 2,
                                importance = "permutation",
                                write.forest = TRUE)),
      
      ## and probability type model for current/normal period predictions
      "lrn_rf_prob" = quote(lrn("classif.ranger", 
                                predict_type = "prob",
                                num.trees = 501,
                                splitrule =  "extratrees",
                                mtry = 4,
                                min.node.size = 2,
                                importance = "permutation",
                                write.forest = TRUE))
    )
  }
  
  if (!suppliedElsewhere("validStrategies", sim)) {
    sim$validStrategies <- list(
      "lrn_rf_resp" = quote(rsmp("repeated_spcv_coords", folds = 10, repeats = 1)),
      "lrn_rf_prob" = quote(rsmp("repeated_spcv_coords", folds = 10, repeats = 1))
    )
  }
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
