defineModule(sim, list(
  name = "BGC_WNA_dataPrep",
  description = paste("Prepares training data to fit a biogeoclimatic (BGC) zones model,",
                      "relating observed BGCs and downscaled climate variables",
                      "(BGC_WNA_model module), as well as data to predict BGC zones under novel",
                      "climate conditions (BGC_WNA_predict module). The training",
                      "data set includes 'current' BGC zone, elevation and climate data obtained for Western Canada and",
                      "US, or a study area in this region, at points spaced in a regular grid (spacing",
                      "can be defined by the user). The 'prediction' dataset includes the same",
                      "climate covariates obtained from a selection of General Circulation Models (GCMs),",
                      "Shared Socioeconomic Pathway scenarios (SSPs), model runs, historic and future periods",
                      "and extracted at the centre points of an hexagonal grid."),
  keywords = c("biogeoclimatic zones", "climate", "data preparation", "CCISS"),
  authors = structure(list(
    list(given = "Ceres", family = "Barros", role = c("aut"), email = "ceres.barros@gov.bc.ca", comment = NULL),
    list(given = "Colin", family = "Mahony", role = c("aut"), email = "colin.mahony@gov.bc.ca", comment = NULL),
    list(given = "Kiri", family = "Daust", role = c("aut", "cre"), email = "kiri.daust@gov.bc.ca", comment = NULL),
    list(given = "Will", family = "MacKenzie", role = c("aut", "cre"), email = "will.mackenzie@gov.bc.ca", comment = NULL)
  ), class = "person"),
  childModules = character(0),
  version = list(BGC_WNA_dataPrep = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "BGC_WNA_dataPrep.Rmd"),
  loadOrder = list(before = c("BGC_WNA_model", "BGC_WNA_predict")),
  reqdPkgs = list(
    ## SpaDES/caching-related:
    "reproducible",
    "SpaDES.core (>= 2.1.5)", 
    ## other:
    "aws.s3",
    "bcgov/climr@main (HEAD)",
    "bcgov/ccissr@development (HEAD)",
    "crayon",
    "data.table",
    "foreach",
    "gdalUtilities",
    "ggplot2",
    "sf",
    "terra",
    "themis",
    "tidymodels"),
  parameters = bindrows(
    defineParameter("assertions", "logical", TRUE, NA, NA, paste("Turn on/off assertions -- i.e. checks ",
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
    defineParameter("GCMs", "character", climr::list_gcms(), NA, NA,
                    paste("List of General Circulation Models (GCMs) to obtain climate projections for",
                          "`sim$predData`.")),
    defineParameter("GCMperiods", "character", climr::list_gcm_periods(), NA, NA,
                    paste("List of future periods to obtain climate projections for `sim$predData`.")),
    defineParameter("GCMruns", "integer", 3L, 1L, Inf,
                    paste("Number of General Circulation Models (GCMs) runs to obtain climate projections for",
                          "`sim$predData`. Note that the ensemble mean across runs will not be used.")),
    defineParameter("gridSize", "integer", 2000L, 1, Inf,
                    paste("Distance in m between points of the regular grid used to extract",
                          "BGC, elevation (used for cliamte downscaling) and climate data to",
                          "fit the BGC zones predictive model.")),
    defineParameter("SSPs", "character", climr::list_ssps(), NA, NA,
                    paste("List of Shared Socioeconomic Pathway scenarios (SSPs) to obtain climate projections for",
                          "`sim$predData`.")),
    defineParameter("nthread", "integer", 1L, 1L, Inf,
                    paste("Passed to `climr::downscale(..., nthread)` to parallelise climate downscaling operations")),
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
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "bgcs", objectClass = "SpatVector", 
                 desc = paste0("A polygon SpatVector of 'current' biogeoclimatic",
                               "zones across Western Canada and the US. This map",
                               "was produced ...TODO. Hosted privately.", 
                               "Please contact developers to obtain access."), 
                 sourceURL = NA),
    expectsInput(objectName = "districts", objectClass = "SpatVector", 
                 desc = paste0("A polygon SpatVector of forest districts in British",
                               "Columbia. Obtained from ...TODO. Hosted privately.", 
                               "Please contact developers to obtain access."), 
                 sourceURL = NA),
    expectsInput(objectName = "elev", objectClass = "SpatRaster", 
                 desc = paste0("A DEM SpatRaster. Elevation values are used to dowscale",
                               "climate data using `climr`. Defaults to Western Canada and",
                               "US raster at 250m resolution from ...TODO. Hosted privately.", 
                               "Please contact developers to obtain access."), 
                 sourceURL = NA),
    expectsInput(objectName = "hexgrid", objectClass = "SpatVector", 
                 desc = paste0("A hexagonal grid covering BC at 400m 'resolution'.", 
                               "Hexagon center points will be used as locations to predict",
                               "BGC zones under novel climate conditions. Hosted privately.", 
                               "Please contact developers to obtain access."), 
                 sourceURL = NA),
    expectsInput(objectName = "studyArea", objectClass = "SpatVector", 
                 desc = paste0("OPTIONAL. A polygon of a study area to subset the training data",
                               "`sim$trainData` to."), 
                 sourceURL = NA),
    expectsInput(objectName = "studyArea_predict", objectClass = "SpatVector", 
                 desc = paste0("OPTIONAL. A polygon of a study area to subset the prediction data",
                               "`sim$predictData` to."), 
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "predData", objectClass = "data.table", 
                  desc = paste("Data used to project biogeoclimatic zones (BGCs)",
                               "under novel climate conditions. Should have the same format",
                               "as `sim$trainData`, but point locations (i.e. rows) correspond to",
                               "center points of hexagons in `sim$hexgrid`.")),
    createsOutput(objectName = "trainData", objectClass = "data.table", 
                  desc = paste("A training dataset to fit predictive models of",
                               "biogeoclimatic zones (BGCs) ~ climate relationships. Includes",
                               "a column of 'current' BGCs (response), columns of longitude ('lon') and latitude ('lat'),", 
                               "columns of climate covariates defined by `P(sim)$climPredictors`",
                               "and a 'CRS' attribute pertaining to the projection used in 'lon', 'lat' columns.",
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


doEvent.BGC_WNA_dataPrep = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- trainingDataEvent(sim)

      ## schedule future events
        sim <- scheduleEvent(sim, time(sim), currentModule(sim), "getPredictionData")
        
      ## schedule future events
      if (P(sim)$balance) {
        sim <- scheduleEvent(sim, time(sim), currentModule(sim), "balanceData")
      }
    },
    getPredictionData = {
      sim <- predictionDataEvent(sim)
    },
    balanceData = {
      sim <- balanceDataEvent(sim)
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
  
  message("Preparing training data...")
  
  ## FITTING DATA -------------------------------
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
  
  ## get climate data
  climData <- Cache(
    getClimate,
    coords = trainCoords, 
    which_refmap = "auto", 
    return_refperiod = TRUE, 
    vars = P(sim)$climrVars, 
    cache = TRUE,
    .cacheExtra = list(summary(trainCoords)),
    userTags = c(cacheTags, "climData"),
    omitArgs = c("userTags", "trainCoords", "bgcs"))
  
  setDT(climData)   ## this shouldn't be necessary, submit issue/reprex to reproducible.
  
  ## add more climate variables
  addVars(climData)
  
  ## subset data to variables of interest
  trainData <- climData[, .SD, .SDcols = c("BGC", "id", P(sim)$climPredictors)]
  trainData <- trainData[complete.cases(trainData)]
  
  ## add back lon lat, and a CRS attribute
  trainData <- trainCoords[trainData, on = "id", nomatch = 0L]
  attr(trainData, "CRS") <- crs("EPSG:4326", proj = TRUE)
  
  ## export to sim -- keep trainDataCRS
  sim$trainData <- trainData
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

balanceDataEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  cacheTags <- c(currentModule(sim), "function:balanceDataEvent")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)

  message("Balancing training data...")
  
  if (!is.na(P(sim)$badBGCs)) {
    message(orange("Excluding the following BGCs from further analyses:"))
    message(orange("  ", paste(badBGCs, collapse = ", ")))
    sim$trainData <- sim$trainData[!BGC %in% badBGCs,]
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
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

predictionDataEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  cacheTags <- c(currentModule(sim), "function:predictionDataEvent")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  
  ## get centroids
  cacheExtra <- summary(sim$hexgrid)
  hexcentroids <- Cache(centroids,
                        x = sim$hexgrid,
                        .cacheExtra = cacheExtra,
                        userTags = c(cacheTags, "hexcentroids"), 
                        omitArgs = c("x", "userTags"))
  .gc()
  
  ## add districts
  hexcentroids <- Cache(terra::intersect,
                        x = hexcentroids, 
                        y = sim$districts, 
                        .cacheExtra = list(summary(hexcentroids), summary(sim$districts)),
                        userTags = c(cacheTags, "hexcentroids_districts"),
                        omitArgs = c("userTags", "x", "y")) 
  
  ## extract elevation with interpolation - crop first
  cacheExtra <- list(summary(sim$elev), summary(hexcentroids)) 
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
                      omitArgs = c("userTags", "x", "y")) |>
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
  browser()  ## out of memory if using a original test SA
  ## get climate data from climr
  ## extract by subsets of points, then write to csv with append = TRUE
  projectedClimate <- getClimate(predCoords,
                                 ## climr args:
                                 which_normal = "normal_composite",
                                 gcm_models = P(sim)$GCMs,
                                 ssp = P(sim)$SSPs, 
                                 gcm_period = P(sim)$GCMperiods, 
                                 max_run = P(sim)$GCMruns,
                                 nthread = P(sim)$nthread,
                                 ## .getClimVars args
                                 byCombo = TRUE, 
                                 outFormat = "disk",
                                 filename = file.path(outputPath(sim), "projectedClimate.csv"))
  
  sim$trainData
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #
  ## BGCs
  dwnldFromObjSto(prefix = "~/CCISS_Working/WNA_BGC/WNA_BGC_v12_5Apr2022",
                  bucket = "gmrtde",
                  path = dPath)
  sim$bgcs <- vect(file.path(dPath, "WNA_BGC_v12_5Apr2022.gpkg"))
  .gc()

  if (!is.null(sim$studyArea)) {
    CRS <- crs(sim$bgcs, proj = TRUE)
    
    ## July 2024 - NRCan network blocking HTTPS, which causes terra to fail to reproject.
    ## see https://stackoverflow.com/questions/76973035/certgetcertificatechain-trust-error-cert-trust-is-untrusted-root-gdal-error-1
    studyArea <- st_as_sf(sim$studyArea)
    studyArea <- st_transform(studyArea, y = CRS)
    
    sim$studyArea <- vect(studyArea)
  } else {
    message(blue("studyArea not provided. No cropping will be done."))
  }
  
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
  
  ## DEM
  dwnldFromObjSto(prefix = "~/DEM/DEM_NorAm/NA_Elevation/data/northamerica",
                  bucket = "gmrtde",
                  path = dPath)
  sim$elev <- rast(file.path(dPath, "northamerica_elevation_cec_2023.tif"))
  
  ## HEX grid
  dwnldFromObjSto(prefix = "~/CCISS_Working/WNA_BGC/HexGrid400m_Sept2021",
                  bucket = "gmrtde",
                  path = dPath)
  sim$hexgrid <- GDALcrop("HexGrid400m_Sept2021.gpkg", "HexGrid400m_Sept2021_crop.gpkg",
                          dPath, studyArea = sim$studyArea_predict) |>
    Cache(userTags = c(cacheTags, "hexgrid"))
  .gc()
  
  ## Forest districts
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
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
