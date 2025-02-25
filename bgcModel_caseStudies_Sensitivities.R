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

# data.table to store results of sensitivity analyses
results <- data.table(
  sensitivity = character(),
  studyname   = character(),
  gapSet      = character(),
  varSet      = character(),
  subsample   = character(),
  maxSample   = integer,
  numTree     = integer(),
  outlierPct  = numeric(),
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
  dem <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling//WNA_DEM_4326_clipped.tif")
  dem <- crop(dem, studyarea)
  X <- dem # template raster for testing
  
  ## make the points table
  points <- as.data.table(dem, cells=T, xy=T)
  colnames(points) <- c("id", "lon", "lat", "elev")
  points <- points[,c(2,3,4,1)] #restructure for climr input
  values(X)[points$id] <- points$el ; plot(X)
  
  ## -------------------------------------------------
  ## attribute the points table with BGC label
  
  bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_22Feb2025_raster.tif")
  bgcs <- crop(bgcs, studyarea)
  points.bgc <- as.data.table(bgcs, cells=T, xy=T)
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
  ## Remove outliers 
  
  vars_simple <- c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt")
  
  # Pull climr data for the subsample:
  clim <- downscale(
    xyz = points,
    which_refmap = "refmap_climr",
    return_refperiod = TRUE, # Also return the 1961-1990 normals period.
    vars = vars_simple
  )
  
  ## Remove outliers (0.0027 excludes 3-sigma outliers)
  clim <- merge(points, clim, by = "id", all = FALSE)
  nonOutlier <- removeOutlier(clim, alpha = .0027, vars = vars_simple) # [colin] changed to vars_simple because vars_seasonal likely would violate multivariate normality. I also added a log-transformation to the remove_outliers function, otherwise ratio variables would violate normality. 
  
  # add the nonOutliers as a new logical field in points
  points[, nonOutlier := nonOutlier[.SD, on = "id", BGC]]
  points[, nonOutlier := !is.na(nonOutlier)]
  values(X) <- NA; values(X)[points$id] <- points$nonOutlier ; plot(X)
  
  ## ----------------------------------------------------------------------------
  ## get a subsample from each BGC unit's grid points based on population size 
  
  # function for subsample size
  subsample <- function(N, threshold = 200, asymptote = 2000, shape = 0.0002) {
    above_thresh <- N > threshold
    N[above_thresh] <- threshold + (asymptote - threshold) * (1 - exp(-shape * (N[above_thresh] - threshold)))
    return(round(N))
  }
  # visualization
  x <- seq(0,10000, 50)
  y <- subsample(x)
  plot(x, y, type = "l", col = "blue", lwd = 2, main = "Subsampling Function (Exponential)",
       xlab = "population size", ylab = "sample size")
  abline(a = 0, b = 1, col = "gray", lty = 2)
  
  # apply function to BGC units
  pointcount <-  table(points[(nonOutlier)]$BGC)
  samplesize1 <- subsample(pointcount)
  points_subsample1 <- points[(nonOutlier), .SD[sample(.N, min(ifelse(BGC %in% names(samplesize1), samplesize1[BGC], .N), .N))], by = BGC]
  dim(points_subsample1)
  
  # visualize
  par(mfrow=c(1,2))
  plot(sort((pointcount)), type = "l", ylab="Number of points")
  lines(sort((samplesize1)), type = "l", lty=2)
  legend("topleft", bty="n", legend = c("BGC grid point count", "sample size"), lty=c(1,2))
  plot(sort(log10(pointcount)), type = "l")
  lines(sort(log10(samplesize1)), type = "l", lty=2)
  par(mfrow=c(1,1))
  
  # add the subsample as a new logical field in points
  points[, subsample1 := points_subsample1[.SD, on = "id", BGC]]
  points[, subsample1 := !is.na(subsample1)]
  values(X) <- NA; values(X)[points$id] <- points$subsample1 ; plot(X)
  
  ## ----------------------------------------------------------------------------
  ## get a subsample from each BGC unit's grid points based on spatial dispersion (clumpiness) of unit
  
  # get a metric of percent of cells that have same adjacent class
  library(landscapemetrics) # for pladj
    pladj <- lsm_c_pladj(bgcs)  # Per-class adjacency metric
  class_mapping <- cats(bgcs)[[1]]  # Extract numeric-class mapping
  pladj <- merge(pladj, class_mapping, by.x = "class", by.y = "value", all.x = TRUE)
  colnames(pladj)[colnames(pladj) == "category"] <- "BGC"
  print(pladj)
   
  # # apply square root function to BGC units
  # pointcount <-  table(points$BGC)
  # samplesize <- pointcount^0.5*10 #square root multiplied by 10 (1:1 for N=100)
  # points_subsample2 <- points[(nonOutlier), .SD[sample(.N, min(ifelse(BGC %in% names(samplesize), samplesize[BGC], .N), .N))], by = BGC]
  # dim(points_subsample2)
  # 
  # # plot of pladj vs BGC subsample size
  # x <- log2(as.vector(table(points_subsample$BGC)))
  # y <- (pladj$value/100)^2
  # plot(x, y, col = "white",
  #      xaxt = "n",
  #      xlab = "Sample size of BGC Unit",
  #      ylab = "pladj")
  # axis(1, at = seq(1,20), labels = round(2^seq(1,20)))
  # text(x, y, labels = pladj$BGC, cex = 0.5)

  # modify square root subsampling with amplifying for pladj complexity
  pointcount <-  table(points[(nonOutlier)]$BGC)
  pladj_factor <- (pladj$value[match(names(pointcount), pladj$BGC)]/100)^2
  samplesize2 <- round(pointcount^0.5*10/pladj_factor) #square root multiplied by 10 (1:1 for N=100)
  points_subsample2 <- points[(nonOutlier), .SD[sample(.N, min(ifelse(BGC %in% names(samplesize2), samplesize2[BGC], .N), .N))], by = BGC]
  dim(points_subsample2)
  
  # add the subsample as a new logical field in points
  points[, subsample2 := points_subsample2[.SD, on = "id", BGC]]
  points[, subsample2 := !is.na(subsample2)]
  values(X) <- NA; values(X)[points$id] <- points$subsample2 ; plot(X)
  
  x <- as.vector(samplesize1)
  y <- as.vector(samplesize2)
  plot(x,y, col="white", xlab="subsample1", ylab="subsample2")
  text(x, y, labels = names(samplesize1), cex = 0.5)
  abline(a = 0, b = 1, col = "gray", lty = 2)
  
  ## -----------------------------------------------------------------------
  ## write points to file
  
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
  
  # # recreate the study area DEM
  # dem <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling//WNA_DEM_4326_clipped.tif")
  # dem <- crop(dem, studyarea)
  # X <- dem # template raster for testing
  # values(X) <- NA
  # 
  # #read in point attributes generated in the last step
  # points <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/points_studyArea.csv")
  # values(X) <- NA; values(X)[points$id] <- points$gap ; plot(X)
  
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
  ## define variable sets to use in sensitivity analyses
  
  vars_simple <- c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", 
                   "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", 
                   "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt")
  
  vars_expert <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", 
                   "EXT", "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", 
                   "SHM", "Tave_sm", "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", 
                   "Tmin_at", "Tmin_sm", "Tmin_sp", "Tmin_wt", "CMI", "PPT_MJ", 
                   "PPT_JAS", "CMD.total")
  
  vars_seasonal <- c(list_vars(set = "Seasonal"), "PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total", "DD_delayed")
  
  
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
  ## Sensitivity Analysis #0 - training subsample
  
  varSet <- "simple"
  
  ## train models in a loop of permutations of subsamples and gaps
  gapSets <- c("WithGaps", "NoGaps")
  subsamples <- c("subsample1", "subsample2")
  gapSet <- "WithGaps"
  for(gapSet in gapSets){
    subsample <- "subsample1"
    for(subsample in subsamples){
      
      trainData <- merge(points[get(subsample)], clim_ref, by="id")
      if(gapSet=="WithGaps") trainData <- trainData[gap==FALSE]
      
      trainData[, BGC := as.factor(BGC)]
      
      # Train model with simple variables, on points from the entire study area:
      BGCmodel <- ranger(
        BGC ~ .,
        data = trainData[, c("BGC", get(paste0("vars_", varSet))), with = FALSE],
        num.trees = 501,
        splitrule =  "extratrees",
        # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
        importance = "permutation",
        write.forest = TRUE,
        classification = TRUE,
        probability = FALSE
      ) 
      # assign(paste("BGCmodel", gapSet, subsample, sep="_"), BGCmodel) # optionally assign the model object
      print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
      
      ## predictions for reference period
      preds <- predict(BGCmodel, data = clim_ref)$prediction
      assign(paste("preds_ref", gapSet, subsample, sep="_"), preds) # assign the vector
      preds_rast <- X
      preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
      preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
      assign(paste("preds_ref", gapSet, subsample, "3857", sep="_"), preds_rast) # assign the raster object
      
      # Compute gap misclassification error
      test_points <- points[gap == TRUE]
      test_pred <- preds[test_points$id]  # Get predicted values for test set
      test_bgcs <- test_points$BGC  # Actual BGC
      gap_error <- mean(test_pred != test_bgcs)
      
      ## predictions for future period
      preds <- predict(BGCmodel, data = clim_proj)$prediction
      assign(paste("preds_proj", gapSet, subsample, sep="_"), preds) # assign the vector
      preds_rast <- X
      preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
      preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
      assign(paste("preds_proj", gapSet, subsample, "3857", sep="_"), preds_rast) # assign the raster object
      
      # Store results
      results <- rbind(results, data.table(
        sensitivity = "subsample",  
        studyname = studyname,  
        gapSet    = gapSet,
        varSet    = varSet,
        subsample = subsample,
        maxSample = 2000,
        numTree   = 501,
        outlierPct= 0.0027,
        oob_error = BGCmodel$prediction.error,
        gap_error = gap_error
      ), fill = TRUE)
      
      print(subsample) 
    }
    print(gapSet) 
  }
  
  ## leaflet map for reference period
  leaflet() %>%
    addTiles() %>%
    addPolygons(data = gap_poly, fillColor = "transparent", color = "black", weight = 2, opacity = 1, fillOpacity = 0, popup = ~paste("Gap Polygon"), group = "Gap Extents") %>%
    addRasterImage(preds_ref_NoGaps_subsample1_3857, colors = color_pal, opacity = 1, group = "Preds: full, subsample1") %>%
    addRasterImage(preds_ref_WithGaps_subsample1_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, subsample1") %>%
    addRasterImage(preds_ref_NoGaps_subsample2_3857, colors = color_pal, opacity = 1, group = "Preds: full, subsample2") %>%
    addRasterImage(preds_ref_WithGaps_subsample2_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, subsample2") %>%
    addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
    addLayersControl(overlayGroups = c("Gap Extents", "Preds: full, subsample1", "Preds: WithGaps, subsample1", "Preds: full, subsample2", "Preds: WithGaps, subsample2", "BGC"), options = layersControlOptions(collapsed = FALSE))
  
  ## leaflet map for future period
  leaflet() %>%
    addTiles() %>%
    addPolygons(data = gap_poly, fillColor = "transparent", color = "black", weight = 2, opacity = 1, fillOpacity = 0, popup = ~paste("Gap Polygon"), group = "Gap Extents") %>%
    addRasterImage(preds_proj_NoGaps_subsample1_3857, colors = color_pal, opacity = 1, group = "Preds: full, subsample1") %>%
    addRasterImage(preds_proj_WithGaps_subsample1_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, subsample1") %>%
    addRasterImage(preds_proj_NoGaps_subsample2_3857, colors = color_pal, opacity = 1, group = "Preds: full, subsample2") %>%
    addRasterImage(preds_proj_WithGaps_subsample2_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, subsample2") %>%
    addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
    addLayersControl(overlayGroups = c("Gap Extents", "Preds: full, subsample1", "Preds: WithGaps, subsample1", "Preds: full, subsample2", "Preds: WithGaps, subsample2", "BGC"), options = layersControlOptions(collapsed = FALSE))
  
  ## -------------------------------------------------
  ## Sensitivity Analysis #1 - variable set and gaps
  
  ## train models in a loop of permutations of variable sets and gaps
  subsample <- "subsample1"
  gapSets <- c("WithGaps", "NoGaps")
  varSets <- c("simple", "expert", "seasonal")
  gapSet <- "WithGaps"
  for(gapSet in gapSets){
    varSet <- "simple"
    for(varSet in varSets){
      
      trainData <- merge(points[(subsample)], clim_ref, by="id")
      if(gapSet=="WithGaps") trainData <- trainData[gap==FALSE]
      
      trainData[, BGC := as.factor(BGC)]
      
      # Train model with simple variables, on points from the entire study area:
      BGCmodel <- ranger(
        BGC ~ .,
        data = trainData[, c("BGC", get(paste0("vars_", varSet))), with = FALSE],
        num.trees = 501,
        splitrule =  "extratrees",
        # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
        importance = "permutation",
        write.forest = TRUE,
        classification = TRUE,
        probability = FALSE
      ) 
      # assign(paste("BGCmodel", gapSet, varSet, sep="_"), BGCmodel) # optionally assign the model object
      print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
      
      ## predictions for reference period
      preds <- predict(BGCmodel, data = clim_ref)$prediction
      assign(paste("preds_ref", gapSet, varSet, sep="_"), preds) # assign the vector
      preds_rast <- X
      preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
      preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
      assign(paste("preds_ref", gapSet, varSet, "3857", sep="_"), preds_rast) # assign the raster object
      
      # Compute gap misclassification error
      test_points <- points[gap == TRUE]
      test_pred <- preds[test_points$id]  # Get predicted values for test set
      test_bgcs <- test_points$BGC  # Actual BGC
      gap_error <- mean(test_pred != test_bgcs)
      
      ## predictions for future period
      preds <- predict(BGCmodel, data = clim_proj)$prediction
      assign(paste("preds_proj", gapSet, varSet, sep="_"), preds) # assign the vector
      preds_rast <- X
      preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
      preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
      assign(paste("preds_proj", gapSet, varSet, "3857", sep="_"), preds_rast) # assign the raster object
      
      # Store results
      results <- rbind(results, data.table(
        sensitivity = "varSet",  
        studyname = studyname,  
        gapSet    = gapSet,
        varSet    = varSet,
        subsample = subsample,
        maxSample = 2000,
        numTree   = 501,
        outlierPct= 0.0027,
        oob_error = BGCmodel$prediction.error,
        gap_error = gap_error
      ), fill = TRUE)
      
      print(varSet) 
    }
    print(gapSet) 
  }
  
  ## leaflet map for reference period
  leaflet() %>%
    addTiles() %>%
    addPolygons(data = gap_poly, fillColor = "transparent", color = "black", weight = 2, opacity = 1, fillOpacity = 0, popup = ~paste("Gap Polygon"), group = "Gap Extents") %>%
    addRasterImage(preds_ref_NoGaps_simple_3857, colors = color_pal, opacity = 1, group = "Preds: full, simple") %>%
    addRasterImage(preds_ref_WithGaps_simple_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, simple") %>%
    addRasterImage(preds_ref_NoGaps_expert_3857, colors = color_pal, opacity = 1, group = "Preds: full, expert") %>%
    addRasterImage(preds_ref_WithGaps_expert_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, expert") %>%
    addRasterImage(preds_ref_NoGaps_seasonal_3857, colors = color_pal, opacity = 1, group = "Preds: full, seasonal") %>%
    addRasterImage(preds_ref_WithGaps_seasonal_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, seasonal") %>%
    addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
    addLayersControl(overlayGroups = c("Gap Extents", "Preds: full, simple", "Preds: WithGaps, simple", "Preds: full, expert", "Preds: WithGaps, expert", "Preds: full, seasonal", "Preds: WithGaps, seasonal", "BGC"), options = layersControlOptions(collapsed = FALSE))
  
  ## leaflet map for future period
  leaflet() %>%
    addTiles() %>%
    addPolygons(data = gap_poly, fillColor = "transparent", color = "black", weight = 2, opacity = 1, fillOpacity = 0, popup = ~paste("Gap Polygon"), group = "Gap Extents") %>%
    addRasterImage(preds_proj_NoGaps_simple_3857, colors = color_pal, opacity = 1, group = "Preds: full, simple") %>%
    addRasterImage(preds_proj_WithGaps_simple_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, simple") %>%
    addRasterImage(preds_proj_NoGaps_expert_3857, colors = color_pal, opacity = 1, group = "Preds: full, expert") %>%
    addRasterImage(preds_proj_WithGaps_expert_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, expert") %>%
    addRasterImage(preds_proj_NoGaps_seasonal_3857, colors = color_pal, opacity = 1, group = "Preds: full, seasonal") %>%
    addRasterImage(preds_proj_WithGaps_seasonal_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, seasonal") %>%
    addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
    addLayersControl(overlayGroups = c("Gap Extents", "Preds: full, simple", "Preds: WithGaps, simple", "Preds: full, expert", "Preds: WithGaps, expert", "Preds: full, seasonal", "Preds: WithGaps, seasonal", "BGC"), options = layersControlOptions(collapsed = FALSE))
  
  
  for(varSet in varSets[-3]){
    
    ## -------------------------------------------------
    ## Sensitivity Analysis #2 - training sample size
    gapSet <- "WithGaps"
    for(gapSet in gapSets){
    maxSamples <- c(50, 200, 500, 2000, 8000, 999999)
    for(maxSample in maxSamples){
      
      trainData <- merge(points, clim_ref, by="id")
      if(gapSet=="WithGaps") trainData <- trainData[gap==FALSE]
      
      samplesize_test <- table(trainData$BGC)
      samplesize_test[samplesize_test>maxSample] <- maxSample
      trainData <- trainData[, .SD[sample(.N, min(ifelse(BGC %in% names(samplesize_test), samplesize_test[BGC], .N), .N))], by = BGC]
      
      trainData[, BGC := as.factor(BGC)]
      
      # Train model with simple variables, on points from the entire study area:
      BGCmodel <- ranger(
        BGC ~ .,
        data = trainData[, c("BGC", get(paste0("vars_", varSet))), with = FALSE],
        num.trees = 501,
        splitrule =  "extratrees",
        # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
        importance = "permutation",
        write.forest = TRUE,
        classification = TRUE,
        probability = FALSE
      ) 
      # assign(paste("BGCmodel", maxSample, sep="_"), BGCmodel) # optionally assign the model object
      print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
      
      ## predictions for reference period
      preds <- predict(BGCmodel, data = clim_ref)$prediction
      assign(paste("preds_ref", maxSample, sep="_"), preds) # assign the vector
      preds_rast <- X
      preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
      preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
      assign(paste("preds_ref", maxSample, "3857", sep="_"), preds_rast) # assign the raster object
      
      # Compute gap misclassification error
      test_points <- points[gap == TRUE]
      test_pred <- preds[test_points$id]  # Get predicted values for test set
      test_bgcs <- test_points$BGC  # Actual BGC
      gap_error <- mean(test_pred != test_bgcs)
      
      ## predictions for future period
      preds <- predict(BGCmodel, data = clim_proj)$prediction
      assign(paste("preds_proj", maxSample, sep="_"), preds) # assign the vector
      preds_rast <- X
      preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
      preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
      assign(paste("preds_proj", maxSample, "3857", sep="_"), preds_rast) # assign the raster object
      
      # Store results
      results <- rbind(results, data.table(
        sensitivity = "maxSample",  
        studyname = studyname,  
        gapSet    = gapSet,
        varSet    = varSet,
        subsample = subsample,
        maxSample = maxSample,
        numTree   = 501,
        outlierPct= 0.0027,
        oob_error = BGCmodel$prediction.error,
        gap_error = gap_error
      ), fill = TRUE)
      
      print(maxSample) 
    }
    print(gapSet) 
    }
    
    ## leaflet map for reference period
    leaflet() %>%
      addTiles() %>%
      addRasterImage(preds_ref_50_3857, colors = color_pal, opacity = 1, group = "Preds: N=50") %>%
      addRasterImage(preds_ref_200_3857, colors = color_pal, opacity = 1, group = "Preds: N=200") %>%
      addRasterImage(preds_ref_500_3857, colors = color_pal, opacity = 1, group = "Preds: N=500") %>%
      addRasterImage(preds_ref_2000_3857, colors = color_pal, opacity = 1, group = "Preds: N=2000") %>%
      addRasterImage(preds_ref_8000_3857, colors = color_pal, opacity = 1, group = "Preds: N=8000") %>%
      addRasterImage(preds_ref_999999_3857, colors = color_pal, opacity = 1, group = "Preds: N=all") %>%
      addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
      addLayersControl(overlayGroups = c("Preds: N=50", "Preds: N=200", "Preds: N=500", "Preds: N=2000", "Preds: N=8000", "Preds: N=all", "BGC"), options = layersControlOptions(collapsed = FALSE))
    
    ## leaflet map for future period
    leaflet() %>%
      addTiles() %>%
      addRasterImage(preds_proj_50_3857, colors = color_pal, opacity = 1, group = "Preds: N=50") %>%
      addRasterImage(preds_proj_200_3857, colors = color_pal, opacity = 1, group = "Preds: N=200") %>%
      addRasterImage(preds_proj_500_3857, colors = color_pal, opacity = 1, group = "Preds: N=500") %>%
      addRasterImage(preds_proj_2000_3857, colors = color_pal, opacity = 1, group = "Preds: N=2000") %>%
      addRasterImage(preds_proj_8000_3857, colors = color_pal, opacity = 1, group = "Preds: N=8000") %>%
      addRasterImage(preds_proj_999999_3857, colors = color_pal, opacity = 1, group = "Preds: N=all") %>%
      addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
      addLayersControl(overlayGroups = c("Preds: N=50", "Preds: N=200", "Preds: N=500", "Preds: N=2000", "Preds: N=8000", "Preds: N=all", "BGC"), options = layersControlOptions(collapsed = FALSE))
    
    ## -------------------------------------------------
    ## Sensitivity Analysis #3 - hyperparameter num.trees
    gapSet <- "NoGaps"
    maxSample <- 2000
    subsample <- "subsample1"
    
    ## train models in a loop of permutations of variable sets and gaps
    numTrees <- c(50, 100, 200, 500)
    for(numTree in numTrees){
      
      trainData <- merge(points[(subsample)], clim_ref, by="id")
      
      samplesize_reduced <- samplesize
      samplesize_reduced[samplesize_reduced>maxSample] <- maxSample
      trainData <- trainData[, .SD[sample(.N, min(ifelse(BGC %in% names(samplesize_reduced), samplesize_reduced[BGC], .N), .N))], by = BGC]
      
      trainData[, BGC := as.factor(BGC)]
      
      # Train model with simple variables, on points from the entire study area:
      BGCmodel <- ranger(
        BGC ~ .,
        data = trainData[, c("BGC", get(paste0("vars_", varSet))), with = FALSE],
        num.trees = numTree,
        splitrule =  "extratrees",
        # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
        importance = "permutation",
        write.forest = TRUE,
        classification = TRUE,
        probability = FALSE
      ) 
      # assign(paste("BGCmodel_numTree", numTree, sep="_"), BGCmodel) # optionally assign the model object
      print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
      
      ## predictions for reference period
      preds <- predict(BGCmodel, data = clim_ref)$prediction
      assign(paste("preds_ref_numTree", numTree, sep="_"), preds) # assign the vector
      preds_rast <- X
      preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
      preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
      assign(paste("preds_ref_numTree", numTree, "3857", sep="_"), preds_rast) # assign the raster object
      
      # Compute gap misclassification error
      test_points <- points[gap == TRUE]
      test_pred <- preds[test_points$id]  # Get predicted values for test set
      test_bgcs <- test_points$BGC  # Actual BGC
      gap_error <- mean(test_pred != test_bgcs)
      
      ## predictions for future period
      preds <- predict(BGCmodel, data = clim_proj)$prediction
      assign(paste("preds_proj_numTree", numTree, sep="_"), preds) # assign the vector
      preds_rast <- X
      preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
      preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
      assign(paste("preds_proj_numTree", numTree, "3857", sep="_"), preds_rast) # assign the raster object
      
      # Store results
      results <- rbind(results, data.table(
        sensitivity = "numTree",  
        studyname = studyname,  
        gapSet    = gapSet,
        varSet    = varSet,
        subsample = subsample,
        maxSample = maxSample,
        numTree   = numTree,
        outlierPct= 0.0027,
        oob_error = BGCmodel$prediction.error,
        gap_error = gap_error
      ), fill = TRUE)
      
      print(numTree) 
    }
    
    
    ## leaflet map for reference period
    leaflet() %>%
      addTiles() %>%
      addRasterImage(preds_ref_numTree_50_3857, colors = color_pal, opacity = 1, group = "Preds: numTree=50") %>%
      addRasterImage(preds_ref_numTree_100_3857, colors = color_pal, opacity = 1, group = "Preds: numTree=100") %>%
      addRasterImage(preds_ref_numTree_200_3857, colors = color_pal, opacity = 1, group = "Preds: numTree=200") %>%
      addRasterImage(preds_ref_numTree_500_3857, colors = color_pal, opacity = 1, group = "Preds: numTree=500") %>%
      addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
      addLayersControl(overlayGroups = c("Preds: numTree=50", "Preds: numTree=100", "Preds: numTree=200", "Preds: numTree=500", "BGC"), options = layersControlOptions(collapsed = FALSE))
    
    ## leaflet map for future period
    leaflet() %>%
      addTiles() %>%
      addRasterImage(preds_proj_numTree_50_3857, colors = color_pal, opacity = 1, group = "Preds: numTree=50") %>%
      addRasterImage(preds_proj_numTree_100_3857, colors = color_pal, opacity = 1, group = "Preds: numTree=100") %>%
      addRasterImage(preds_proj_numTree_200_3857, colors = color_pal, opacity = 1, group = "Preds: numTree=200") %>%
      addRasterImage(preds_proj_numTree_500_3857, colors = color_pal, opacity = 1, group = "Preds: numTree=500") %>%
      addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
      addLayersControl(overlayGroups = c("Preds: numTree=50", "Preds: numTree=100", "Preds: numTree=200", "Preds: numTree=500", "BGC"), options = layersControlOptions(collapsed = FALSE))
    
    
    
    
    ## -------------------------------------------------
    ## Sensitivity Analysis #4 - outlier removal

    ## train models in a loop of permutations of variable sets and gaps
    gapSets <- c("WithGaps", "NoGaps")
    outlierPcts <- c(0, 0.0027, 0.05, 0.32)
    gapSet <- "WithGaps"
    for(gapSet in gapSets){
      outlierPct <- 0.32
      for(outlierPct in outlierPcts){
        
        trainData <- merge(points, clim_ref, by="id", all = FALSE)
        
        ## Remove outliers from the full population
        trainData <- removeOutlier(trainData, alpha = outlierPct, vars = vars_simple) 
        
        ## get a subsample from each BGC unit's grid points based on population size 
        pointcount <-  table(trainData$BGC)
        samplesize <- subsample(pointcount, threshold = 200, asymptote = 2000, shape = 0.0001)
        trainData <- trainData[, .SD[sample(.N, min(ifelse(BGC %in% names(samplesize), samplesize[BGC], .N), .N))], by = BGC]
        
        if(gapSet=="WithGaps") trainData <- trainData[gap==FALSE]
        
        trainData[, BGC := as.factor(BGC)]
        
        # Train model with simple variables, on points from the entire study area:
        BGCmodel <- ranger(
          BGC ~ .,
          data = trainData[, c("BGC", get(paste0("vars_", varSet))), with = FALSE],
          num.trees = 501,
          splitrule =  "extratrees",
          # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
          importance = "permutation",
          write.forest = TRUE,
          classification = TRUE,
          probability = FALSE
        ) 
        # assign(paste("BGCmodel_outlierPct", outlierPct, gapSet, sep="_"), BGCmodel) # optionally assign the model object
        print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
        
        ## predictions for reference period
        preds <- predict(BGCmodel, data = clim_ref)$prediction
        assign(paste("preds_ref_outlierPct", outlierPct, gapSet, sep="_"), preds) # assign the vector
        preds_rast <- X
        preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
        preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
        assign(paste("preds_ref_outlierPct", outlierPct, gapSet, "3857", sep="_"), preds_rast) # assign the raster object
        
        # Compute gap misclassification error
        test_points <- points[gap == TRUE]
        test_pred <- preds[test_points$id]  # Get predicted values for test set
        test_bgcs <- test_points$BGC  # Actual BGC
        gap_error <- mean(test_pred != test_bgcs)
        
        ## predictions for future period
        preds <- predict(BGCmodel, data = clim_proj)$prediction
        assign(paste("preds_proj_outlierPct", outlierPct, gapSet, sep="_"), preds) # assign the vector
        preds_rast <- X
        preds_rast[points$id] <- factor(preds, levels = subzones_colours_ref$BGC)
        preds_rast <- project(preds_rast, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 
        assign(paste("preds_proj_outlierPct", outlierPct, gapSet, "3857", sep="_"), preds_rast) # assign the raster object
        
        # Store results
        results <- rbind(results, data.table(
          sensitivity = "outlierPcts",  
          studyname = studyname,  
          gapSet    = gapSet,
          varSet    = varSet,
          subsample = subsample,
          maxSample = maxSample,
          numTree   = 501,
          outlierPct= outlierPct,
          oob_error = BGCmodel$prediction.error,
          gap_error = gap_error
        ), fill = TRUE)
        
        print(outlierPct) 
      }
      print(gapSet) 
    }
    
    ## leaflet map for reference period
    leaflet() %>%
      addTiles() %>%
      addPolygons(data = gap_poly, fillColor = "transparent", color = "black", weight = 2, opacity = 1, fillOpacity = 0, popup = ~paste("Gap Polygon"), group = "Gap Extents") %>%
      addRasterImage(preds_ref_outlierPct_0_NoGaps_3857, colors = color_pal, opacity = 1, group = "Preds: Full, 0% outliers") %>%
      addRasterImage(preds_ref_outlierPct_0_WithGaps_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, 0% outliers") %>%
      addRasterImage(preds_ref_outlierPct_0.0027_NoGaps_3857, colors = color_pal, opacity = 1, group = "Preds: Full, 0.3% outliers") %>%
      addRasterImage(preds_ref_outlierPct_0.0027_WithGaps_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, 0.3% outliers") %>%
      addRasterImage(preds_ref_outlierPct_0.05_NoGaps_3857, colors = color_pal, opacity = 1, group = "Preds: Full, 5% outliers") %>%
      addRasterImage(preds_ref_outlierPct_0.05_WithGaps_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, 5% outliers") %>%
      addRasterImage(preds_ref_outlierPct_0.32_NoGaps_3857, colors = color_pal, opacity = 1, group = "Preds: Full, 32% outliers") %>%
      addRasterImage(preds_ref_outlierPct_0.32_WithGaps_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, 32% outliers") %>%
      addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
      addLayersControl(overlayGroups = c("Gap Extents", "Preds: Full, 0% outliers", "Preds: WithGaps, 0% outliers", "Preds: Full, 0.3% outliers", "Preds: WithGaps, 0.3% outliers", "Preds: Full, 5% outliers", "Preds: WithGaps, 5% outliers", "Preds: Full, 32% outliers", "Preds: WithGaps, 32% outliers", "BGC"), options = layersControlOptions(collapsed = FALSE))
    
    ## leaflet map for reference period
    leaflet() %>%
      addTiles() %>%
      addPolygons(data = gap_poly, fillColor = "transparent", color = "black", weight = 2, opacity = 1, fillOpacity = 0, popup = ~paste("Gap Polygon"), group = "Gap Extents") %>%
      addRasterImage(preds_proj_outlierPct_0_NoGaps_3857, colors = color_pal, opacity = 1, group = "Preds: Full, 0% outliers") %>%
      addRasterImage(preds_proj_outlierPct_0_WithGaps_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, 0% outliers") %>%
      addRasterImage(preds_proj_outlierPct_0.0027_NoGaps_3857, colors = color_pal, opacity = 1, group = "Preds: Full, 0.3% outliers") %>%
      addRasterImage(preds_proj_outlierPct_0.0027_WithGaps_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, 0.3% outliers") %>%
      addRasterImage(preds_proj_outlierPct_0.05_NoGaps_3857, colors = color_pal, opacity = 1, group = "Preds: Full, 5% outliers") %>%
      addRasterImage(preds_proj_outlierPct_0.05_WithGaps_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, 5% outliers") %>%
      addRasterImage(preds_proj_outlierPct_0.32_NoGaps_3857, colors = color_pal, opacity = 1, group = "Preds: Full, 32% outliers") %>%
      addRasterImage(preds_proj_outlierPct_0.32_WithGaps_3857, colors = color_pal, opacity = 1, group = "Preds: WithGaps, 32% outliers") %>%
      addRasterImage(bgc_ref, colors = color_pal, opacity = 1, group = "BGC") %>%
      addLayersControl(overlayGroups = c("Gap Extents", "Preds: Full, 0% outliers", "Preds: WithGaps, 0% outliers", "Preds: Full, 0.3% outliers", "Preds: WithGaps, 0.3% outliers", "Preds: Full, 5% outliers", "Preds: WithGaps, 5% outliers", "Preds: Full, 32% outliers", "Preds: WithGaps, 32% outliers", "BGC"), options = layersControlOptions(collapsed = FALSE))
    
    
    
    print(varSet)
  }
  
  
  print(studyname)
} #end of biggest loop

write.csv(results, "//objectstore2.nrs.bcgov/ffec/BGC_models/sensitivity_results_v2.csv", row.names = FALSE)


## -------------------------------------------------
## -------------------------------------------------
## STEP 3 - sensitivity analysis results
## -------------------------------------------------
## -------------------------------------------------

results <- fread("//objectstore2.nrs.bcgov/ffec/BGC_models/sensitivity_results.csv")

## -------------------------------------------------
## Plot 1 - gap error for the maxSample experiment

# Filter data for maxSample sensitivity
maxsample_results <- results[results$sensitivity == "maxSample" & gapSet =="NoGaps", ]

# Get unique studyname and varSet combinations
study_var_combinations <- unique(maxsample_results[, c("studyname", "varSet")])

# Set up colors and line types
colors <- rainbow(nrow(study_var_combinations))  # Unique colors for each combination
line_types <- 1:length(unique(maxsample_results$studyname))  # Different line types for studyname

# Open a blank plot
plot(
  NULL, NULL, 
  xlim = range(maxsample_results$maxSample), 
  ylim = range(maxsample_results$gap_error), 
  log = "x",  # Log scale for maxSample
  xlab = "Max Sample", 
  ylab = "Gap Error", 
  main = "Gap Error vs Max Sample"
)

# Loop through each combination and plot lines
for (i in seq_len(nrow(study_var_combinations))) {
  subset_data <- maxsample_results[
    maxsample_results$studyname == study_var_combinations$studyname[i] &
      maxsample_results$varSet == study_var_combinations$varSet[i], 
  ]
  
  lines(
    subset_data$maxSample, subset_data$gap_error, 
    col = colors[i], 
    lty = line_types[which(unique(maxsample_results$studyname) == study_var_combinations$studyname[i])],
    type = "o",  # Points and lines
    pch = 16  # Solid dots
  )
}

# Add legend
legend(
  "topright", 
  legend = paste(study_var_combinations$studyname, "-", study_var_combinations$varSet), 
  col = colors, 
  lty = rep(line_types, length.out = length(colors)), 
  pch = 16
)

## -------------------------------------------------
## Plot 2 - gap error for the varSet experiment

# Filter data for varSet sensitivity
varset_results <- results[results$sensitivity == "varSet", ]

# Create a combined label for studyname and gapSet
varset_results$study_gap <- paste(varset_results$studyname, varset_results$gapSet, sep = "\n")

# Compute mean gap error for each (varSet, study_gap) combination
agg_data <- aggregate(gap_error ~ varSet + study_gap, data = varset_results, FUN = mean, na.rm = TRUE)

# Convert to wide format
library(reshape2)
gap_error_wide <- dcast(agg_data, varSet ~ study_gap, value.var = "gap_error", fill = NA)

# Convert to matrix (exclude varSet column)
gap_error_matrix <- as.matrix(gap_error_wide[, -1])

# Set row names to varSet values
rownames(gap_error_matrix) <- gap_error_wide$varSet

# Define colors for bars
colors <- c("gray25", "gray50", "gray75")

# Create grouped bar plot
barplot(
  gap_error_matrix, beside = TRUE, col = colors,
  names.arg = colnames(gap_error_matrix), las = 1,
  ylab = "Gap Error",
  main = "Gap Error by varSet, studyname, and gapSet",
  legend.text = rownames(gap_error_matrix),
  args.legend = list(x = "topleft", bty = "n", title="Variable set")
)

