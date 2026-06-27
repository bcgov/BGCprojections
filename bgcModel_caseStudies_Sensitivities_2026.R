## BGC projections sensitivity analyses
## Colin Mahony colin.mahony@gov.bc.ca
## February 2025

library(climr)
library(data.table)
library(terra)
library(sf)
library(foreach) # for outlier removal function

# Source functions: 
source("utils.R")

# dir <- "//objectstore2.nrs.bcgov/ffec/BGC_models/" 
dir <- "C:/Users/CMAHONY/Data/BGC_models/" #local copy, for speed

# Source functions TODO: move these into ccissr as independent functions
source("utils.R")

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

## -------------------------------------------------
## data.table to store results of sensitivity analyses

results <- data.table(
  sensitivity = character(),
  studyname   = character(),
  varSet      = character(),
  num.trees     = integer(),
  num.random.splits = integer(),
  mtry = character(),
  min.node.size = integer(),
  replace = logical(),
  sample.fraction = numeric(),
  oob_error   = numeric(),
  gap_error   = numeric()
)

#regional study areas
studynames <- c("Bamfield", "Kamloops", "Pemberton")
studyname <- "Pemberton"
#loop regional study areas
for(studyname in studynames){
  
  ## -------------------------------------------------
  ## -------------------------------------------------
  ## STEP 1 - training and prediction data - Attribute a table of points that correspond to the raster cells of a high-resolution study area
  ## -------------------------------------------------
  ## -------------------------------------------------
  
  ## -------------------------------------------------
  ## Study area DEM and points table
  studyarea <- if(studyname=="Bamfield") ext(c(-125.25, -124, 48.5, 49.125)) else 
    if(studyname=="Kamloops") ext(c(-121, -120, 50.5, 51)) else 
      if(studyname=="Pemberton") ext(c(-124, -122, 50, 51)) else 
        NULL
  
  # create a DEM
  dem.studyarea <- crop(dem, studyarea)
  X <- dem.studyarea # template raster for testing
  
  ## make the points table
  points <- as.data.table(dem.studyarea, cells=T, xy=T)
  colnames(points) <- c("id", "lon", "lat", "elev")
  points <- points[,c(2,3,4,1)] #restructure for climr input
  values(X)[points$id] <- points$el ; plot(X)
  
  ## -------------------------------------------------
  ## attribute the points table with BGC label
  
  bgcs.studyarea <- crop(bgcs, studyarea)
  points.bgc <- as.data.table(bgcs.studyarea, cells=T, xy=T)
  colnames(points.bgc) <- c("id", "lon", "lat", "BGC")
  points.bgc[, BGC := factor(BGC)]
  points[, BGC := points.bgc[.SD, on = "id", BGC]]
  points <- points[!is.na(points$BGC), ]
  values(X) <- NA; values(X)[points$id] <- points$BGC ; plot(X)
  
  ## -------------------------------------------------
  ## divide study area into training and testing areas using an internal checkerboard
  
  points_sf <- st_as_sf(points, coords = c("lon", "lat"), crs = 4326)
  
  # Make rectangular gap extents within the bounding box of the study area. 5L is the default number of gaps to create. 
  gapextents <- makeGapExtents(studyarea=studyarea, 5L)
  
  # Convert list of spatial extents into to polygons: 
  gap_poly <- lapply(gapextents, vect)
  
  # Combine all individual polygons into one spatial object (convert to sf for speed). 
  gap_poly <- st_as_sf(do.call(rbind, gap_poly))
  st_crs(gap_poly) <- st_crs(points_sf) # assign CRS
  plot(gap_poly, add=T)
  
  # identify points that fall within the gap polygons.
  points_gaps <- st_intersection(points_sf, gap_poly)
  plot(points_gaps, add=T)
  points_gaps_dt <- as.data.table(points_gaps)
  
  # add the gaps as a logical field in points
  points[, gap := points_gaps_dt[.SD, on = "id", elev]]
  points[, gap := !is.na(gap)]
  values(X) <- NA; values(X)[points$id] <- points$gap ; plot(X)
  
  ## -------------------------------------------------
  ## training sample (use the method developed in bgcModel_WNA_sensitivities_Apr2026.R)
  
  dem.training <- dem.studyarea
  values(dem.training)[points$id][points$gap] <- NA ; plot(dem.training)
  
  bgcs.training <- bgcs.studyarea
  bgcs.training[points$id[points$gap]] <- NA ; plot(bgcs.training)
  
  
  trainsample <- bgc_trainingSample(dem.training, bgcs.training, bgcs_info = bgcs_info,
                                    scheme = "squareRoot", 
                                    squareRoot.multiplier = 10, 
                                    removeOutliers = TRUE, alpha.BC = .0027, alpha.nonBC = .05, 
                                    removeOutliers.vars = c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt"), 
                                    climaticVariance = TRUE, climaticVariance.var = "MAT",
                                    plotDiagnostics = TRUE, 
                                    plot.dir = "C:/Users/CMAHONY/Data/BGC_models/", 
                                    plot.name = "diagnostics_v4_Pemberton"
  )
  dim(trainsample)
  
  # add the subsample as a new logical field in points
  points[, trainsample := trainsample[.SD, on = "id", BGC]]
  points[, trainsample := !is.na(trainsample)]
  values(X) <- NA; values(X)[points$id] <- points$trainsample ; plot(X)
  
  # write.csv(points, "//objectstore2.nrs.bcgov/ffec/BGC_models/points_studyArea.csv", row.names = FALSE)
  
  
  ## -------------------------------------------------
  ## -------------------------------------------------
  ## STEP 2 - Sensitivity Analyses
  ## -------------------------------------------------
  ## -------------------------------------------------
  
  library(climr)
  library(data.table)
  library(terra)
  library(ranger) # For RF
  library(caret) # For confusionMatrix()
  library(leaflet)
  
  
  ## -------------------------------------------------
  ## define variable sets
  vars_simple <- c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", 
                   "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", 
                   "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt")
  
  vars_expert <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", 
                   "EXT", "MWMT", "NFFD_sm", "NFFD_sp", "PAS_an", "PAS_sp", 
                   "SHM", "Tave_sm", "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin_an", 
                   "Tmin_at", "Tmin_sm", "Tmin_sp", "Tmin_wt", "CMI_an", "PPT_MJ", 
                   "PPT_JAS", "CMD.total")

  vars_seasonal <- c(list_vars(set = "Seasonal"), "PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total")
  
  vars_bioclimate <- c("DD5_an", "CMD.total", "PPT_an", "MCMT", "MWMT", "TD", "FFP", "MSP")
  
  
  
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
    gcms = list_gcms()[5],
    ssps = list_ssps()[2],
    gcm_periods = list_gcm_periods()[3],
    run_nm = list_runs_ssp(list_gcms()[5], list_ssps()[2])[3],
    which_refmap = "refmap_climr",
    return_refperiod = FALSE, # Also return the 1961-1990 normals period.
    vars = list_vars()
  )
  ccissr::addVars(clim_proj)
  
  ## -------------------------------------------------
  ## training data
  
  trainData <- merge(points[trainsample==TRUE], clim_ref, by="id")
  
  ## -------------------------------------------------
  ## Colors and Baseline BGC raster
  
  ## color scheme for bgc units
  subzones_colours_ref <- fread("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISSv13_latest_tool_materials/WNAv13_Subzone_colours_2.csv") %>% 
    dplyr::select(!c(fid, MAP_LABEL, NSRNAME, ZONE)) %>% 
    dplyr::mutate(BGC_num = as.numeric(as.factor(BGC)))
  color_pal <- colorFactor(
    palette = subzones_colours_ref$RGB,
    domain = subzones_colours_ref$BGC_num
  )
  
  ## baseline BGC raster for leaflet maps
  bgc_ref <- X
  bgc_ref[points$id] <- factor(points$BGC, levels = subzones_colours_ref$BGC)
  bgc_ref <- project(bgc_ref, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
  plot(bgc_ref)
  
  
  ## -------------------------------------------------
  ## define sensitivities and hyperparameter settings
  ## -------------------------------------------------
  
  sensitivities <- c("num.trees", "num.random.splits", "mtry", "min.node.size", "replace", "sample.fraction", "varSet")
  
  ## -------------------------------------------------
  ## Sensitivity Analysis #1 - hyperparameter num.trees
  
  ## train models in a loop of permutations of variable sets and gaps
  parameters.num.trees <- c(100, 200, 500, 1000)
  
  ## -------------------------------------------------
  ## Sensitivity Analysis #2 - hyperparameter num.random.splits
  # num.random.splits (most important for ExtraTrees)
  # Controls how many random thresholds are tried per variable
  # 
  # Effect:
  #   
  #   Low → more randomness → higher bias, lower variance, faster
  #   High → closer to Gini behavior → lower bias, slower
  # 
  # Guidelines:
  #   
  #   Start: 5–10
  #   If accuracy is lagging: increase toward 10–20
  # 
  # This is the main “dial” between Gini-like and ExtraTrees-like behavior.
  
  parameters.num.random.splits <- c(1,3,5,10)
  
  ## -------------------------------------------------
  ## Sensitivity Analysis #3 - hyperparameter mtry
  # Number of variables tried per split
  # 
  # Effect with extratrees:
  #   
  #   Still important, but randomness from splits reduces sensitivity
  # 
  # Guidelines:
  #   
  #   Regression: mtry ≈ p/3 (default is fine as a baseline)
  #   High collinearity (climate predictors): try smaller values (e.g., p/5, sqrt(p))
  # 
  # Smaller mtry + extratrees = strong decorrelation (often good for spatial transferability)
  
  parameters.mtry <- c(10,5,3,1)
  
  ## -------------------------------------------------
  ## Sensitivity Analysis #4 - hyperparameter min.node.size
  # Minimum samples per terminal node
  # 
  # Effect:
  #   
  #   Small → complex trees (risk overfitting spatial noise)
  #   Large → smoother predictions (often desirable for rasters)
  # 
  # Guidelines for spatial data:
  #   
  #   Start: 5–20
  #   Large N (millions): 20–100
  #   If maps look “speckled”: increase this
  
  parameters.min.node.size <- c(1,3,10,50)
  
  ## -------------------------------------------------
  ## Sensitivity Analysis #5 - hyperparameter replace
  # Whether sampling is with replacement
  # Recommendation:
  #   
  #   For extratrees: often replace = FALSE
  #   This mimics classical ExtraTrees behavior and:
  #     reduces variance
  #     improves speed
  
  parameters.replace <- c(TRUE,FALSE,TRUE,FALSE)
  
  ## -------------------------------------------------
  ## Sensitivity Analysis #6 - hyperparameter sample.fraction
  # sample.fraction
  # Fraction of data used per tree
  # 
  # Key interaction:
  #   
  #   ExtraTrees already adds randomness → you can often:
  #     avoid full bootstrap
  #     reduce sample size per tree
  # 
  # Guidelines:
  #   
  #   Try: 0.5–0.8
  #   For huge datasets: even 0.2–0.5 can work
  # 
  # Big win for speed without much accuracy loss
  
  parameters.sample.fraction <- c(0.2,0.4,0.6,0.8)
  
  ## -------------------------------------------------
  ## Sensitivity Analysis #7 - variable sets

  parameters.varSet <- c("simple", "expert", "seasonal", "bioclimate")

  ## -------------------------------------------------
  ## loop through the sensitivities
  ## -------------------------------------------------
  
  for(sensitivity in sensitivities){
    
    parameters <- get(paste("parameters", sensitivity, sep="."))
    
    for(parameter in parameters){
      
      # Train model with simple variables, on points from the entire study area:
      BGCmodel <- ranger(
        BGC ~ .,
        data = trainData[, c("BGC", get(paste0("vars_", if(sensitivity=="varSet") parameter else "expert"))), with = FALSE],
        splitrule =  "extratrees",
        num.trees = if(sensitivity=="num.trees") parameter else 1000,
        num.random.splits = if(sensitivity=="num.random.splits") parameter else 1,
        mtry = if(sensitivity=="mtry") parameter else floor(length(vars_expert)^0.5),
        min.node.size = if(sensitivity=="min.node.size") parameter else 1,
        replace = if(sensitivity=="replace") parameter else FALSE,
        sample.fraction = if(sensitivity=="sample.fraction") parameter else ifelse(sensitivity == "replace" & parameter == TRUE, 1, 0.632),
        importance = "none",
        write.forest = TRUE,
        classification = TRUE,
        probability = FALSE
      ) 
      print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
      
      ## predictions for reference period
      preds <- predict(BGCmodel, data = clim_ref)$prediction
      assign(paste("preds_ref", sensitivity, parameter, sep="_"), preds) # assign the vector
      preds_rast <- X
      preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
      preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
      assign(paste("preds_ref", sensitivity, parameter, "3857", sep="_"), preds_rast) # assign the raster object
      
      # Compute gap misclassification error
      test_points <- points[gap == TRUE]
      test_pred <- preds[test_points$id]  # Get predicted values for test set
      test_bgcs <- test_points$BGC  # Actual BGC
      gap_error <- mean(test_pred != test_bgcs)
      
      ## predictions for future period
      preds <- predict(BGCmodel, data = clim_proj)$prediction
      assign(paste("preds_proj", sensitivity, parameter, sep="_"), preds) # assign the vector
      preds_rast <- X
      preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
      preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
      assign(paste("preds_proj", sensitivity, parameter, "3857", sep="_"), preds_rast) # assign the raster object
      
      # Store results
      results <- rbind(results, data.table(
        sensitivity = sensitivity,  
        studyname = studyname,  
        varSet   = if(sensitivity=="varSet") parameter else "expert",
        num.trees = if(sensitivity=="num.trees") parameter else 1000,
        num.random.splits = if(sensitivity=="num.random.splits") parameter else 1,
        mtry = if(sensitivity=="mtry") parameter else floor(length(vars_expert)^0.5),
        min.node.size = if(sensitivity=="min.node.size") parameter else 1,
        replace = if(sensitivity=="replace") parameter else FALSE,
        sample.fraction = if(sensitivity=="sample.fraction") parameter else ifelse(sensitivity == "replace" & parameter == TRUE, 1, 0.632),
        oob_error = BGCmodel$prediction.error,
        gap_error = gap_error
      ), fill = TRUE)
      
      print(parameter) 
    }
    
    ## leaflet map for reference period
    leaflet() %>%
      addTiles() %>%
      addPolygons(data = gap_poly, fillColor = "transparent", color = "black", weight = 2, opacity = 1, fillOpacity = 0, popup = ~paste("Gap Polygon"), group = "Gap Extents") %>%
      addRasterImage(get(paste("preds_ref", sensitivity, parameters[1], "3857", sep="_")), colors = color_pal, opacity = 1, group = paste0(sensitivity, "=", parameters[1])) %>%
      addRasterImage(get(paste("preds_ref", sensitivity, parameters[2], "3857", sep="_")), colors = color_pal, opacity = 1, group = paste0(sensitivity, "=", parameters[2])) %>%
      addRasterImage(get(paste("preds_ref", sensitivity, parameters[3], "3857", sep="_")), colors = color_pal, opacity = 1, group = paste0(sensitivity, "=", parameters[3])) %>%
      addRasterImage(get(paste("preds_ref", sensitivity, parameters[4], "3857", sep="_")), colors = color_pal, opacity = 1, group = paste0(sensitivity, "=", parameters[4])) %>%
      addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
      addLayersControl(overlayGroups = c(paste0(sensitivity, "=", parameters[1]), paste0(sensitivity, "=", parameters[2]), paste0(sensitivity, "=", parameters[3]), paste0(sensitivity, "=", parameters[4]), "BGC"), options = layersControlOptions(collapsed = FALSE))
    
    ## leaflet map for future period
    leaflet() %>%
      addTiles() %>%
      addPolygons(data = gap_poly, fillColor = "transparent", color = "black", weight = 2, opacity = 1, fillOpacity = 0, popup = ~paste("Gap Polygon"), group = "Gap Extents") %>%
      addRasterImage(get(paste("preds_proj", sensitivity, parameters[1], "3857", sep="_")), colors = color_pal, opacity = 1, group = paste0(sensitivity, "=", parameters[1])) %>%
      addRasterImage(get(paste("preds_proj", sensitivity, parameters[2], "3857", sep="_")), colors = color_pal, opacity = 1, group = paste0(sensitivity, "=", parameters[2])) %>%
      addRasterImage(get(paste("preds_proj", sensitivity, parameters[3], "3857", sep="_")), colors = color_pal, opacity = 1, group = paste0(sensitivity, "=", parameters[3])) %>%
      addRasterImage(get(paste("preds_proj", sensitivity, parameters[4], "3857", sep="_")), colors = color_pal, opacity = 1, group = paste0(sensitivity, "=", parameters[4])) %>%
      addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
      addLayersControl(overlayGroups = c(paste0(sensitivity, "=", parameters[1]), paste0(sensitivity, "=", parameters[2]), paste0(sensitivity, "=", parameters[3]), paste0(sensitivity, "=", parameters[4]), "BGC"), options = layersControlOptions(collapsed = FALSE))
    
    print(sensitivity)
  } #end of biggest loop
  
  print(studyname)
} #end of biggest loop

write.csv(results, "//objectstore2.nrs.bcgov/ffec/BGC_models/sensitivity_results_2026_hyperparameters.csv", row.names = FALSE)

## -------------------------------------------------
## Conclusions: 

# num.trees: 
    # larger num.trees reduces obb error moderately and gap error marginally. 
    # larger num.trees reduces speckling marginally in baseline prediction, but not beyond 500 trees
    # larger num.trees reduces speckling noticeably in extrapolation (future), with 1000 better than 500
    # Conclusion: use num.trees=1000. 
      #this is because the dataset is large and unbalanced, and splitrule = "extratrees" needs more trees to stabilize than Gini

# num.random.splits: 
    # default is 1, but i tried higher values
    # 3 gave .0025 lower oob error but 0.002 higher gap error
    # very little effect on speckling 
    # noticeable effect on future prediction, with 1 giving more spatially cohesive results. 
    # stick with 1

# mtry: 
    # default is sqrt(number of predictors). smaller numbers handle collinearity better
    # smaller mtry reduces oob error slightly (0.004) but increases gap error marginally (0.001). 
    # substantial effect on future prediction (extrapolation), with mtry=10 showing more bec change in valleys
    # conclusion: stick with default sqrt(N). 

# min.node.size: 
    # default is 1. large values give less complex trees, meaning less risk of overfitting and higher risk of speckling
    # high node size (50) gave much less speckling, but a large increase in oob error (0.075) and gap error (0.025)
    # increased node size did not noticably improve extrapolation behaviour
    # stick with default

# replace: 
    # default is TRUE. 
    # very minor differences in baseline, with slightly lower errors for FALSE (0.001)
    # FALSE gives spatially more cohesive extrapolation
    # go with FALSE. apparently this is acceptable with extratrees rule. 

# sample.fraction: 
    # 0.6 is optimal for error, but 0.2 reduces speckling substantially in baseline
    # not much difference in extrapolation
    # conclusion: stick with default 0.632
