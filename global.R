
repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
install.packages("SpaDES.project", repos = repos)

library(SpaDES.project)

out <- SpaDES.project::setupProject(
  name = "BGCprojections",
  paths = list("inputPath" = "data/",
               "outputPath" = "outputs/",
               "modulePath" = "modules/",
               "packagePath" = "packages/", 
               "cachePath" = "reproducible.cache/",
               "projectPath" = "."),
  options = list("reproducible.cachePath" = "reproducible.cache/",
                 "reproducible.destinationPath" = paths$inputPath,
                 "reproducible.useMemoise" = TRUE,
                 "climr.cache.path" = "climr.cache/",
                 "repos" = repos),
  ## install, but don't load these:
  packages = c("future", 
               "sf",
               "CeresBarros/climr@devl"), ## override
  modules = c("bcgov/BGCprojections@main/modules/BGCWNAmodel_trainAndProj"),
  overwrite = FALSE,
  Restart = TRUE,
  times = list(start = 1, end = 1), 
  ## reduce things for testing
  params = list(
    BGCWNAmodel_trainAndProj = list(
      "GCMperiods" = climr::list_gcm_periods()[1:2],
      "GCMruns" = 1,
      "GCMs" = climr::list_gcms()[c(1, 4, 5)],
      "SSPs" = climr::list_ssps()[2:4],
      "nthread_mlr3" = list("lrn_rf_resp" = 3)
    )
  ),
  studyArea = {
    # terra::vect(terra::ext(c(-125, -112, 43, 55)), crs = "EPSG:4326")  ## not enough mem with 32Gb
    terra::vect(terra::ext(c(-122, -118, 46, 50)), crs = "EPSG:4326") ## south BC, w/ US
  },
  studyArea_predict = {
    studyArea
  },
  modelsToTrain = {
    list(
      "lrn_rf_resp" = quote(mlr3::lrn("classif.ranger",
                                      predict_type = "response",
                                      num.trees = 501,
                                      splitrule =  "extratrees",
                                      mtry = 4,
                                      min.node.size = 2,
                                      importance = "permutation",
                                      write.forest = TRUE))
    )
  },
  validStrategies = {
    list(
      "lrn_rf_resp" = quote(mlr3::rsmp("repeated_spcv_coords", folds = 10, repeats = 1))
    ) 
  },
  
  BGCprojections_simInit <- SpaDES.core::simInit2(out)
  
  # unloadNamespace("ccissr")
  # unloadNamespace("climr")
  # devtools::load_all("C:/Users/cbarros/GitHub/BCGov/climr/")
  # devtools::load_all("C:/Users/cbarros/GitHub/BCGov/ccissr/")
  # 
  # SpaDES.core::restartSpades()
  
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
  
