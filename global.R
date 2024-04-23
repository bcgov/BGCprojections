## check we're "above" the SpaDES project level
if (basename(getwd()) != "CCISS") {
  stop("Please make sure workign directory is at the GH repository/CCISS folder level")
} else {
  rver <- paste0(version[["major"]], ".", strsplit(version[["minor"]], "[.]")[[1]][1])
  pkgPath <- normalizePath(file.path("packages", rver), winslash = "/", mustWork = FALSE)
  dir.create(pkgPath, recursive = TRUE, showWarnings = FALSE)
  .libPaths(pkgPath)
  
  # if (!requireNamespace("Require", quietly = TRUE) |
  # tryCatch(packageVersion("Require", lib.loc = pkgPath) < "0.3.1.9041", error = function(e) TRUE)) {
  remotes::install_github("PredictiveEcology/Require@simplify2", update = FALSE)
  # }
  
  if (!requireNamespace("SpaDES.project", quietly = TRUE) | 
      tryCatch(packageVersion("SpaDES.project", lib.loc = pkgPath) < "0.0.8.9031", error = function(e) TRUE)) {
    remotes::install_github("PredictiveEcology/SpaDES.project@transition", update = FALSE)
  }
  library(SpaDES.project)
  
  out <- setupProject(
    name = "CCISS_BGCworkflow",
    paths = list("inputPath" = "data/",
                 "outputPath" = "outputs/",
                 "modulePath" = "modules/",
                 "packagePath" = pkgPath, 
                 "cachePath" = "reproducible.cache/"),
    options = list("reproducible.cachePath" = "reproducible.cache/",
                   "reproducible.destinationPath" = paths$inputPath,
                   "climr.cache.path" = "climr.cache/"),
    ## install, but don't load these:
    packages = c("future", 
                 "PredictiveEcology/reproducible@modsForLargeArchives",
                 "sf"), 
    modules = c("BGC_WNA_dataPrep"),
    Restart = FALSE,
    times = list(start = 1, end = 1), 
    studyArea =  {
      # terra::vect(terra::ext(c(-125, -112, 43, 55)), crs = "EPSG:4326")  ## not enough mem with 32Gb
      terra::vect(terra::ext(c(-122, -118, 46, 50)), crs = "EPSG:4326") ## south BC, w/ US
    },
    studyArea_predict = {
      studyArea
    }
  )
  
  BGCprojections_simInit <- do.call(SpaDES.core::simInit, out) |>
    reproducible::Cache(userTags = "simInit")
  # Mar15 11:09:53 simInit Global parameter(s) not used in any module: .studyAreaName.
  # Elpsed time for simInit: 1.321427 mins
  # Saving large object (fn: SpaDES.core::simInit, cacheId: a44bf279ada6880e) to Cache: 1.1 Gb Done!
  #   There were 26 warnings (use warnings() to see them)
  # > Sys.time()
  # [1] "2024-03-15 11:17:09 PDT"
  
  BGCprojections_simOut <- SpaDES.core::spades(BGCprojections_simInit)
  
  # ## install/load packages
  # Require::Require(c(
  #   "bcgov/climr@devl (HEAD)", 
  #   "data.table",
  #   "foreach",
  #   "ggplot2", 
  #   "reproducible",
  #   # "smotefamily", 
  #   "mlr3learners",
  #   "mlr3spatial", 
  #   "mlr3spatiotempcv",
  #   "mlr3viz",
  #   "SpaDES.core", 
  #   "terra", 
  #   "themis",
  #   "tidymodels"
  # ))
  
  ## install, but don't load these.
  # Require::Require(c(
  #   "future",
  #   "sf"), 
  #   require = FALSE)
  
  ## in 03_RunPredictHex -- not sure what we'll need
  # require(foreach)
  # require(dplyr)
  # require(reshape2)
  # library(doParallel)
  # library(tidyr)
  # require(sf)
  # require(RPostgreSQL)
  # library(disk.frame)
  # require(RPostgres)
  
  
  
  
  ## source scripts
  source("SpaDES/R/02_Build_WNA_BGC_trainingset.R")
  source("SpaDES/R/03_RunPredictHex.R")
  source("SpaDES/R/04_CreateBGCfutsMap.R")  
}
