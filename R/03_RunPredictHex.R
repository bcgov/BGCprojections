## read in climateBC files, predict, and send to db
## Kiri Daust, Dec 2020

## This script makes all predictions to the grid of hex centroids (from 400m hex grid)
## and uploads them to the server

## Load hex grid from where we'll extract the centroids
## get hexgrid from object storage
dPath <- unlist(options("reproducible.destinationPath"))
## get file from object storage and save locally. then use prepInputs for the rest
dwnldFromObjSto(prefix = "~/CCISS_Working/WNA_BGC/HexGrid400m_Sept2021",
                bucket = "gmrtde",
                path = dPath)

studyArea_predict <- project(studyArea, y = crs(vect(file.path(dPath, "HexGrid400m_Sept2021.gpkg")), proj = TRUE))

hexgrid <- GDALcrop("HexGrid400m_Sept2021.gpkg", "HexGrid400m_Sept2021_crop.gpkg",
                   dPath, studyArea = studyArea_predict) |>
  Cache(userTags = "hexgrid")

## get centroids
cacheExtra <- list(summary(hexgrid))
hexcentroids <- Cache(centroids,
                      x = hexgrid,
                      .cacheExtra = cacheExtra,
                      userTags = "hexcentroids", 
                      omitArgs = c("x", "userTags"))
.gc()

## add forest regions
dwnldFromObjSto(prefix = "~/CommonTables/Forest_Region_District",
                bucket = "gmrtde",
                path = dPath)

## crop first to save memory  during intersect
districts <- GDALcrop("Forest_Region_District.gpkg", "Forest_Region_District_crop.gpkg",
                      dPath, studyArea = hexcentroids) |>
  Cache(userTags = "districts")

districts <- districts["ORG_UNIT"]
districts$ORG_UNIT <- as.character(districts$ORG_UNIT)
districts$ORG_UNIT[districts$ORG_UNIT == "DSS"] <- "CAS"
names(districts) <- "dist_code"

.gc()
hexcentroids <- Cache(terra::intersect,
                      x = hexcentroids, 
                      y = districts, 
                      .cacheExtra = list(summary(hexcentroids), summary(districts)),
                      userTags = c("hexcentroids_districts"),
                      omitArgs = c("userTags", "x", "y")) 

## add elevation
dwnldFromObjSto(prefix = "~/DEM/DEM_NorAm/NA_Elevation/data/northamerica",
                bucket = "gmrtde",
                path = dPath)
elev <- rast(file.path(dPath, "northamerica_elevation_cec_2023.tif"))
cacheExtra <- list(summary(elev), summary(hexcentroids)) 
elev <- Cache(postProcessTerra,
              from = elev,
              cropTo = hexcentroids,
              projectTo = hexcentroids,
              maskTo = NA,
              .cacheExtra = cacheExtra,
              userTags = "elev",
              omitArgs = c("from", "userTags", "cropTo", "projectTo", "maskTo"))
.gc()

## need extract here because interpolation is necessary
predCoords <- Cache(terra::extract,
                   x = elev, 
                   y = hexcentroids, 
                   method = "bilinear",
                   cells = TRUE,
                   xy =  TRUE,
                   .cacheExtra = list(summary(hexcentroids), summary(elev)),
                   userTags = "predCoords",
                   omitArgs = c("userTags", "x", "y")) |>
  as.data.table()
predCoords[, cell := NULL]
predCoords[, dist_code := hexcentroids$dist_code]

## change column names to climr standard
setnames(predCoords, c("id", "elev", "lon", "lat", "dist_code"))

## exclude NAs
predCoords <- predCoords[complete.cases(predCoords)]

## project to EPSG:4326
points <- vect(predCoords[, .(lon, lat)], crs = crs(elev, proj = TRUE))
points <- project(points, crs("EPSG:4326", proj = TRUE))

predCoords[, `:=`(lon = geom(points)[, "x"],
                  lat = geom(points)[, "y"])]
rm(points)
.gc()

## get climate data from climr
gcms <- list_gcm()
runs <- 3 ## exclude ensembleMean later and note that not all GCMs have 3 runs
ssps <- list_ssp()
periods <- list_gcm_period()

projectedClimate <- getClimate(predCoords,
                               which_normal = "normal_composite",
                               gcm_models = gcms,
                               ssp = ssps, 
                               gcm_period = periods, 
                               max_run = runs)

##helper predict function if the tiles are too big for prediction all at once
##doesn't return anything, adds in place

# might still need to predict (on hex grid centroids.) per tile, to avoid exhausting memory
# Y1 are all predictors
# note: may not have enough RAM to build full prediction table
# note: upload each predicted table to server -- we will have future tables N = 13 GCMs X 5 SSP x 3 runs x 5 periods + historic tables

tile_predict <- function(Y1, maxSize = 6000000){
  n = nrow(Y1)
  brks <- seq(1, n, by = maxSize)
  brks <- c(brks, n)
  Y1[, BGC.pred := NA_character_]
  ## TODO: lapply, maybe parallelise?
  for (j in 1:(length(brks)-1)) {
    browser() ## use mlr3 clas prediction model here
    Y1[brks[j]:brks[j+1], BGC.pred := predict(BGCmodel2, Y1[brks[j]:brks[j+1],-c(1:3)],num.threads = 14)[['predictions']]]
  }
  
  return(Y1)
}


##Current Period Probability
## for current period we want a probability type of prediction, so need to fit another full model

setwd("~/BC_HexGrid/")
load("./WNAv12_Subzone_11_Var_ranger_probability_24Aug21.Rdata")
varImport <- c("ID1","ID2","AHM", "CMD_sp", "DD5", "DD5_sm", 
               "DD5_sp", "MCMT", "MWMT", "NFFD", "PAS_sp", 
               "PAS_wt", "SHM", "PPT_at","PPT_wt","PPT05", "PPT06", "PPT07", "PPT08", 
               "PPT09", "CMD07","CMD","DD_0_at","DD_0_wt","CMI","PAS","EXT","bFFP")
files <- list.files("./ClimBC_Out/")  ## this is getting all the tile files.
setwd("~/BC_HexGrid/ClimBC_Out/")
for(input in filesUse){
  cat("Processing",input)
  dat <- fread(input,select = varImport)
  addVars(dat)
  vars <- BGCmodel2p[["forest"]][["independent.variable.names"]]
  varList = c("SiteNo", "BGC", vars)
  setnames(dat, old = c("ID1","ID2"), new = c("SiteNo", "BGC"))
  Y1 <- dat
  Y1=Y1[,..varList]
  Y1 <- na.omit(Y1)
  
  ## TODO: use mlr3 probability model here.
  grid.pred <- predict(BGCmodel2p, data = Y1[,-c(1:2)])
  predMat <- grid.pred$predictions
  predMat <- as.data.table(predMat) 
  predMat[,`:=`(SiteNo = Y1$SiteNo,BGC = Y1$BGC)]
  predMat <- melt(predMat, id.vars = c("SiteNo","BGC"), variable.name = "Pred",value.name = "Prob")
  predMat <- predMat[Prob > 0.1,]
  setorder(predMat,SiteNo)
  predMat[,ProbAdj := Prob/sum(Prob), by = SiteNo]
  predMat[,Prob := NULL]
  #predMat[bgcID, bgc := i.ID2, on = c(ID = "ID1")]
  predMat[,period := "1991"] ## change to 2001-2020 or "2001"
  # setnames(predMat,old = "SiteNo",new = "OldIdx")
  # 
  # newID <- fread("RCB_CrosswalkTable.csv")
  # predMat[newID, SiteNo := i.Index, on = "OldIdx"]
  # predMat[,OldIdx := NULL]
  setcolorder(predMat, c("SiteNo","period","BGC","Pred","ProbAdj"))
  setnames(predMat,c("siteno","period","bgc","bgc_pred","prob"))
  con <- dbConnect(drv, user = "postgres", host = "138.197.168.220",password = "PowerOfBEC", port = 5432, dbname = "cciss") ### for local use
  dbWriteTable(con, "cciss_prob12", predMat,row.names = F, append = T)
  dbDisconnect(con)
  rm(dat,Y1,grid.pred,predMat)
  gc()
}


## Normal Period - with probability type model
tableName <- "cciss_historic"
datDir <- "~/Desktop/BCHex_ClimateBC/Normal/"
varImport <- c("ID1","ID2","AHM", "CMD_sp", "DD5", "DD5_sm", 
               "DD5_sp", "MCMT", "MWMT", "NFFD", "PAS_sp", 
               "PAS_wt", "SHM", "PPT_at","PPT_wt","PPT05", "PPT06", "PPT07", "PPT08", 
               "PPT09", "CMD07","CMD","DD_0_at","DD_0_wt","CMI","PAS")

## do by chunk not tiles.
for(i in 0:13){
  # cat("Processing tile",i,"... \n")
  # IDName <- "NewID"
  # if(i == 0) IDName <- "ID1"
  # varImport[1] <- IDName
  dat <- fread(paste0(datDir,"Tile",i,"_Norm.csv"),select = varImport) ##point to climateBC data
  dat <- fread("~/BC_HexGrid/ClimBC_Out/RCB_ClimBC_Normal_1991_2020MSY.csv",select = varImport)
  addVars(dat)
  
  vars <- BGCmodel2[["forest"]][["independent.variable.names"]]
  varList = c("SiteNo", "BGC", vars)
  setnames(dat, old = c("ID1","ID2"), new = c("SiteNo", "BGC"))
  Y1 <- dat
  Y1=Y1[,..varList]
  Y1 <- na.omit(Y1)
  
  ##Predict normal period subzones######
  Y1[,BGC.pred := predict(BGCmodel2, Y1[,-c(1:2)])[['predictions']]]
  gc()
  Y1[,Period := "Normal61"]
  Y1 <- Y1[,.(Period,SiteNo,BGC,BGC.pred)]
  setnames(Y1, c("period","siteno","bgc","bgc_pred"))
  dbWriteTable(con, tableName, Y1,row.names = F, append = T)
  rm(Y1,dat)
  gc()
  
}

## Current Period 1991-2019
###This will use a probability model and weight the historic BGC slightly for conservatism.
### The probabilities will be used in the species selection apps this will Likely need to store the data in the future table as there will now be multiple plausible BGC states
### In addition however, a majority vote layer (incorporating the present BGC weighting) should also be created for use in the CreateBCMaps function - this could remained stored in the historic period tables
### For the CreateBCMaps
tableName <- "cciss_historic"
datDir <-"E:/BGC_Hex/BCHex_ClimateBC/Current/"
varImport <- c("ID1","ID2","AHM", "bFFP", "CMD_sp", "DD5", "DD5_sm", 
               "DD5_sp", "Eref_sm", "Eref_sp", "MCMT", "MWMT", "NFFD", "PAS_sp", 
               "PAS_wt", "SHM", "Tmax_sm","PPT_at","PPT_wt","PPT05", "PPT06", "PPT07", "PPT08", 
               "PPT09", "CMD07","CMD","DD_0_at","DD_0_wt")
for(i in 0:13){
  cat("Processing tile",i,"... \n")
  IDName <- "NewID"
  if(i == 0) IDName <- "ID1"
  varImport[1] <- IDName
  dat <- fread(paste0(datDir,"Tile",i,"_Curr.csv"),select = varImport)
  Y1 <- addVars(dat)
  vars <- BGCmodel2[["forest"]][["independent.variable.names"]]
  varList = c("SiteNo", "BGC", vars)
  setnames(Y1, old = c(IDName,"ID2"), new = c("SiteNo", "BGC"))
  Y1=Y1[,..varList]
  Y1 <- Y1[Tmax_sm > -100,]
  Y1 <- na.omit(Y1)
  
  ##Predict 1991-2019 subzones######
  Y1[,BGC.pred := predict(BGCmodel2, Y1[,-c(1:2)])[['predictions']]]
  gc()
  Y1[,Period := "Current91"]
  Y1 <- Y1[,.(Period,SiteNo,BGC,BGC.pred)]
  setnames(Y1, c("period","siteno","bgc","bgc_pred"))
  dbWriteTable(con, tableName, Y1,row.names = F, append = T)
  rm(Y1,dat)
  gc()
  
}

###now do joins and create indices
dbExecute(con, "create table historic_sf as select cciss_historic.*,grid_dist.dist_code,grid_dist.geom from cciss_historic,grid_dist where cciss_historic.siteno = grid_dist.siteno")
dbExecute(con, "create table future_sf as select cciss_future.*,grid_dist.dist_code,grid_dist.geom from cciss_future,grid_dist where cciss_future.siteno = grid_dist.siteno")
dbExecute(con, "create index fut_sf_idx on future_sf (dist_code,scenario,gcm,futureperiod)")
dbExecute(con, "create index hist_sf_idx on historic_sf (dist_code,period)")
dbExecute(con, "create index fut_idx on cciss_future (siteno)")
dbExecute(con, "create index hist_idx on cciss_historic (siteno)")
