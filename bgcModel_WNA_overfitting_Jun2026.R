## test if BGC model for WNA show evidence of overfitting due to large training sample
## Colin Mahony colin.mahony@gov.bc.ca
## June 2026

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

## -------------------------------------------------
## -------------------------------------------------
## STEP 1a - Training samples - remove points for cross-validation
## -------------------------------------------------
## -------------------------------------------------

# import training sample
trainingSample_V4a <- fread(paste0(dir, "points_WNA_v4a.csv"))
bgcs_nonBC <- bgcs_info[grep("USA_|AB_", bgcs_info$Source), BGC]
dim(trainingSample_V4a)

## ----------------------------
## 2. Define CV grid

bdy <- vect("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/bdy.BC.shp")
bdy <- project(bdy, "EPSG:4326")
studyarea <- ext(bdy)

gapwidth <- 0.025
spacing <- 2L

nx <- round(1 / gapwidth)
ny <- nx

xbreaks <- seq(studyarea[1], studyarea[2], length.out = nx + 1)
ybreaks <- seq(studyarea[3], studyarea[4], length.out = ny + 1)

trainingSample_V4a[, inside :=
                     lon >= min(xbreaks) & lon <= max(xbreaks) &
                     lat >= min(ybreaks) & lat <= max(ybreaks)
]

## ----------------------------
## GRID INDEXING (only valid inside)
trainingSample_V4a[, i := findInterval(lon, xbreaks, all.inside = TRUE)]
trainingSample_V4a[, j := findInterval(lat, ybreaks, all.inside = TRUE)]

## ----------------------------
## GAP MASK (ONLY INSIDE STUDY AREA)
trainingSample_V4a[, gap :=
                     inside &
                     ((i - 1) %% spacing == 0) &
                     ((j - 1) %% spacing == 0)
]

## drop gap points and clean up
trainingSample_V4a <- trainingSample_V4a[gap == FALSE]
trainingSample_V4a[, c("inside", "i", "j", "gap") := NULL]

## ----------------------------
## plot points
plot(trainingSample_V4a[, .(lon, lat)], pch=16, cex=0.1, xlim=c(-140, -114), ylim=c(47, 60))

write.csv(trainingSample_V4a, paste0(dir, "points_WNA_v4a_CV.csv"), row.names = FALSE)
dim(trainingSample_V4a)


## -------------------------------------------------
## -------------------------------------------------
## STEP 1b - Training samples - downsample asymptotically
## -------------------------------------------------
## -------------------------------------------------

#asymptotic downsampling of BC units with the default limit of N=2000
trainingSample_V4b <- trainingSample_V4a[
  , {
    if (!(BGC[1] %in% bgcs_nonBC)) {
      n_keep <- subsample_asymptotic(.N)
      if (n_keep <= 0) return(NULL)
      .SD[sample(.N, n_keep)]
    } else {
      .SD
    }
  },
  by = BGC
]
write.csv(trainingSample_V4b, paste0(dir, "points_WNA_v4b_CV.csv"), row.names = FALSE)
dim(trainingSample_V4b)

#asymptotic downsampling of BC units with a limit of N=1000
trainingSample_V4c <- trainingSample_V4a[
  , {
    n_keep <- subsample_asymptotic(.N, asymptote = 1000)
    if (n_keep <= 0) return(NULL)
    .SD[sample(.N, n_keep)]
  },
  by = BGC
]
write.csv(trainingSample_V4c, paste0(dir, "points_WNA_v4c_CV.csv"), row.names = FALSE)
dim(trainingSample_V4c)

#asymptotic downsampling of BC units with a limit of N=500
trainingSample_V4d <- trainingSample_V4a[
  , {
    n_keep <- subsample_asymptotic(.N, asymptote = 500)
    if (n_keep <= 0) return(NULL)
    .SD[sample(.N, n_keep)]
  },
  by = BGC
]
write.csv(trainingSample_V4d, paste0(dir, "points_WNA_v4d_CV.csv"), row.names = FALSE)
dim(trainingSample_V4d)
plot(trainingSample_V4d[, .(lon, lat)], pch=16, cex=0.1)



## -------------------------------------------------
## -------------------------------------------------
## STEP 1c - confirm training sample
## -------------------------------------------------
## -------------------------------------------------

points <- fread(paste0(dir, "points_WNA_v4d_CV.csv"))
popn <- setDT(freq(bgcs))
head(popn)
samp <- as.data.table(table(points$BGC))
samp

setnames(samp, "V1", "BGC")
samp[, count := popn[.SD, on = .(value = BGC), x.count]]
samp[, count := popn[.SD, on = .(value = BGC), x.count]]

bgcs_nonBC <- bgcs_info[grep("USA_|AB_", bgcs_info$Source), BGC]


# figdir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
# png(filename=paste0(figdir, "sampleSize_v4a.png",sep="."), type="cairo", units="in", width=6.5, height=6.5, pointsize=12, res=300)
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

# dev.off()
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
## v4.2a RF model 

#read in training sample generated in the last step
points <- fread(paste0(dir, "points_WNA_v4a_CV.csv"))

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
  splitrule =  "gini", # way faster than gini, not likely a big performance difference, but we should test this. 
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = 'none',
  write.forest = TRUE,
  classification = TRUE,
  num.threads = 14,  # Adjust to fewer threads to reduce memory usage
  probability = FALSE, 
  keep.inbag = FALSE, 
) 
print(paste0("OOB error: ", round(100*BGCmodel$prediction.error, 2), "%")) # OOB prediction error
saveRDS(BGCmodel, paste0(dir, "BGCmodel_WNA_v4a_CV.rds"))



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

BGCmodel_v4.2 <- readRDS(paste0(dir, "BGCmodel_WNA_V4.2gini.rds"))
BGCmodel_v4.2a <- readRDS(paste0(dir, "BGCmodel_WNA_v4a_CV.rds"))
BGCmodel_v4.2b <- readRDS(paste0(dir, "BGCmodel_WNA_v4b_CV.rds"))
BGCmodel_v4.2c <- readRDS(paste0(dir, "BGCmodel_WNA_v4c_CV.rds"))
BGCmodel_v4.2d <- readRDS(paste0(dir, "BGCmodel_WNA_V4d_CV.rds")) 


## -------------------------------------------------
## data

bdy <- vect("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/bdy.BC.shp")
bdy <- project(bdy, "EPSG:4326")
dem <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_DEM_4326_clipped.tif")

# reduce resolution selecting center cell rather than taking an average, to avoid elevation error
fact <- 5
center <- (fact^2 + 1) / 2
dem <- aggregate(dem, fact = fact, fun = function(x, ...) x[center])

dem <- mask(dem, bdy) 
dem <- crop(dem, ext(bdy) )
# plot(dem)
X <- dem # template raster for testing
values(X) <- NA

studyarea <- ext(bdy)

bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_13Apr2026_raster.tif")
bgcs <- crop(bgcs, ext(bdy))
bgcs <- project(bgcs, dem, method = "near")
bgcs <- mask(bgcs, dem) 
bgc_levels <- cats(bgcs)[[1]]  # Extract the category mapping
vals_char <- bgc_levels$BGC[match(values(bgcs), bgc_levels$value)] # Convert numeric values to character labels
values(bgcs) <- factor(vals_char, levels = subzones_colours_ref$BGC) # refactor with full set of levels
bgcs_3857 <- project(bgcs, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## study area points
points <- as.data.table(dem, cells=T, xy=T)
colnames(points) <- c("id", "lon", "lat", "elev")
points <- points[,c(2,3,4,1)] #restructure for climr input
values(X) <- NA; values(X)[points$id] <- points$el ; plot(X)

# extract BGC values and add to the points table
points.bgc <- as.data.table(bgcs, cells=T, xy=T)
colnames(points.bgc) <- c("id", "lon", "lat", "BGC")
points[, BGC := points.bgc[.SD, on = "id", BGC]]
points <- points[!is.na(points$BGC), ]
points[, BGC := factor(BGC)]

## identify holdouts (gaps) created in step1a
points[, inside :=
         lon >= studyarea[1] & lon <= studyarea[2] &
         lat >= studyarea[3] & lat <= studyarea[4]
]
points[, i := findInterval(lon, xbreaks, all.inside = TRUE)]
points[, j := findInterval(lat, ybreaks, all.inside = TRUE)]
points[, gap :=
         inside &
         ((i - 1) %% spacing == 0) &
         ((j - 1) %% spacing == 0)
]
points[, c("inside","i","j") := NULL]



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
clim_ref[!is.finite(Eref_sp), Eref_sp := 0]




## -------------------------------------------------
## predictions

## v4.2 model predictions for reference period
preds_ref_v4.2_vec <- predict(BGCmodel_v4.2, data = clim_ref)$prediction
preds_ref_v4.2 <- X
preds_ref_v4.2[points$id] <- factor(preds_ref_v4.2_vec, levels = subzones_colours_ref$BGC)
preds_ref_v4.2 <- project(preds_ref_v4.2, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v4.2 model predictions for reference period
preds_ref_v4.2a_vec <- predict(BGCmodel_v4.2a, data = clim_ref)$prediction
preds_ref_v4.2a <- X
preds_ref_v4.2a[points$id] <- factor(preds_ref_v4.2a_vec, levels = subzones_colours_ref$BGC)
preds_ref_v4.2a <- project(preds_ref_v4.2a, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 

## v4.2b model predictions for reference period
preds_ref_v4.2b_vec <- predict(BGCmodel_v4.2b, data = clim_ref)$prediction
preds_ref_v4.2b <- X
preds_ref_v4.2b[points$id] <- factor(preds_ref_v4.2b_vec, levels = subzones_colours_ref$BGC)
preds_ref_v4.2b <- project(preds_ref_v4.2b, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries.

## v4.2c model predictions for reference period
preds_ref_v4.2c_vec <- predict(BGCmodel_v4.2c, data = clim_ref)$prediction
preds_ref_v4.2c <- X
preds_ref_v4.2c[points$id] <- factor(preds_ref_v4.2c_vec, levels = subzones_colours_ref$BGC)
preds_ref_v4.2c <- project(preds_ref_v4.2c, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries.

## v4.2d model predictions for reference period
preds_ref_v4.2d_vec <- predict(BGCmodel_v4.2d, data = clim_ref)$prediction
preds_ref_v4.2d <- X
preds_ref_v4.2d[points$id] <- factor(preds_ref_v4.2d_vec, levels = subzones_colours_ref$BGC)
preds_ref_v4.2d <- project(preds_ref_v4.2d, "EPSG:3857", method = "near") #have to manually resample to web mercator with nearest neighbour sampling otherwise leaflet will do so using bilinear interpolation which corrupts the factor levels at polygon boundaries. 


## -------------------------------------------------
## analysis of holdout error
## -------------------------------------------------



trainingSample_V4a <- fread(paste0(dir, "points_WNA_v4a_CV.csv"))
bgc_levels <- unique(trainingSample_V4a$BGC)

# table for predictions
results <- data.table(
  truth = factor(points$BGC, levels = bgc_levels),
  pred_v4.2  = factor(preds_ref_v4.2_vec, levels = bgc_levels),
  pred_v4.2a  = factor(preds_ref_v4.2a_vec, levels = bgc_levels),
  pred_v4.2b  = factor(preds_ref_v4.2b_vec, levels = bgc_levels),
  pred_v4.2c  = factor(preds_ref_v4.2c_vec, levels = bgc_levels),
  pred_v4.2d  = factor(preds_ref_v4.2d_vec, levels = bgc_levels),
  gap   = points$gap
)

models <- c("v4.2d", "v4.2c", "v4.2b", "v4.2a")
asymptote <- c("n=500", "n=1000", "n=2000", expression(n == infinity))

error.fullModel <- mean(results$truth != results$pred_v4.2, na.rm=T)

error_total <- rbindlist(lapply(models, function(m) {
  pred <- results[[paste0("pred_", m)]]
  truth <- results$truth
  data.table(
    model = m,
    gap = results$gap,
    correct = pred == truth
  )[, .(error = 1 - mean(correct, na.rm = TRUE)), by = .(model, gap)]
}))

err <- dcast(error_total, model ~ gap, value.var = "error")
colnames(err) <- c("model", "error_traindata", "error_holdout")
err[, model := factor(model, levels = models)]
setorder(err, model)

figdir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(figdir, "HoldoutError_BC.png"), type="cairo", units="in", width=4.5, height=4, pointsize=10, res=300)

par(oma = c(0,0,0,0), mar=c(3,3,1,1), mgp=c(1.5, 0.25, 0), tck=-0.005, oma = c(0, 0, 0, 0))
ylim=c(0, max(err$error_holdout)*1.05)
plot(1:length(models), col="white", ylim=ylim, axes=F, ylab="", xlab="Maximum training sample per BGC unit", yaxs="i")
lines(err[,error_holdout])
lines(err[,error_traindata], lty=2)
points(err$error_holdout, pch=22, bg="gray")
points(err$error_traindata, pch=21, bg="white")
points(length(models), error.fullModel, pch=8)
axis(1, at=1:length(models), labels = asymptote)
axis(2, at=pretty(ylim), labels = pretty(ylim), las=2)
par(mgp=c(2, 0.25, 0)); title(ylab="Baseline prediction error")
legend("bottomright", c("Holdout regions", "Training regions", "Full model"), pch=c(22,21,8), pt.bg = c("gray", "white", NA), lty=c(1,2,NA), bty="n")
box()


#map of holdouts
par(mar=c(0,0,0,0), plt= c(0.15, 0.55, 0.125, 0.55), new=TRUE)
image(dem, col="white", axes=F)
plot(bdy, add=T, lwd=0.4)
plot(gap_clipped, add=T, col="gray", lwd=0.4)
# legend("bottomleft", legend="Training sample holdouts", fill = "gray", bty="n")
box()

dev.off()


## -------------------------------------------------
## Zone summary analysis of holdout error (holdout error is by subzone-variant, but summarized by zone)

error_bgc <- rbindlist(lapply(models, function(m) {
  pred <- results[[paste0("pred_", m)]]
  truth <- results$truth
  data.table(
    model = m,
    gap = results$gap,
    zone = wna_bgcs[match(results$truth, wna_bgcs$BGC), Zone],
    correct = pred == truth
  )[, .(error = 1 - mean(correct, na.rm = TRUE)), by = .(model, zone, gap)]
}))

figdir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(figdir, "HoldoutError_zone.png"), type="cairo", units="in", width=6.5, height=6.5, pointsize=10, res=300)

par(oma = c(4, 3, 0, 0), mfrow=c(4,4), mar = c(0, 0, .5, .5), mgp=c(1.5, 0.25, 0), cex=1, tck=-0.005, font=1)
ylim=c(0, max(error_bgc$error)*1.02)
zones <- rev(c("CDF", "CWH", "MH", "ESSF", "MS", "IDF", "PP", "BG", "ICH", "SBPS", "SBS", "BWBS", "SWB", "CMA", "BAFA", "IMA"))
asymptote <- c("500", "1000", "2000", "N/A")
for(zone_sel in zones){
  
  err <- dcast(error_bgc[zone==zone_sel], model ~ gap, value.var = "error")
  colnames(err) <- c("model", "error_traindata", "error_holdout")
  err[, model := factor(model, levels = models)]
  setorder(err, model)
  
  plot(1:length(models), col="white", ylim=ylim, axes=F, ylab="", xlab="", yaxs="i")
  lines(err[,error_holdout])
  lines(err[,error_traindata], lty=2)
  points(err$error_holdout, pch=22, bg="gray")
  points(err$error_traindata, pch=21, bg="white")
  axis(1, at=1:length(models), labels = asymptote, outer=TRUE, las=2)
  axis(2, at=pretty(ylim), labels = pretty(ylim), las=2, outer=TRUE)
  mtext(zone_sel, side = 3, line = -1.5, font=2, adj=0.95)
  box()
  if(zone_sel == zones[1]) legend("bottomleft", c("Holdout regions", "Training regions"), pch=c(22,21), pt.bg = c("gray", "white"), lty=c(1,2), bty="n", cex=0.9)
}

mtext("Maximum training sample per BGC unit", side = 1, outer = TRUE, line = 2.5)
mtext("Baseline prediction error", side = 2, outer = TRUE, line = 2)
box()

dev.off()

## -------------------------------------------------
## leaflet map

diff(xbreaks)[1]
diff(ybreaks)[1]
diff(xbreaks)[1]*111*cos(54.5 * pi / 180)
diff(ybreaks)[1]*111

cells <- CJ(i = 1:nx, j = 1:ny)

cells[, gap := (i - 1) %% spacing == 0 &
        (j - 1) %% spacing == 0]

gap_cells <- cells[gap == TRUE]
gap_sf <- st_sfc(lapply(seq_len(nrow(gap_cells)), function(k) {
  
  i <- gap_cells$i[k]
  j <- gap_cells$j[k]
  
  st_polygon(list(matrix(c(
    xbreaks[i],     ybreaks[j],
    xbreaks[i + 1], ybreaks[j],
    xbreaks[i + 1], ybreaks[j + 1],
    xbreaks[i],     ybreaks[j + 1],
    xbreaks[i],     ybreaks[j]
  ), ncol = 2, byrow = TRUE)))
  
}), crs = 4326)

gap_vect <- vect(gap_sf)
gap_clipped <- intersect(gap_vect, bdy)

map <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%  # Add satellite imagery
  addTiles(group = "OSM Basemap") %>%  # Keep OpenStreetMap as an option
  addRasterImage(preds_ref_v4.2, colors = color_pal, opacity = 1, group = "v4.2: baseline") %>%
  addRasterImage(preds_ref_v4.2a, colors = color_pal, opacity = 1, group = "v4.2a: baseline") %>%
  addRasterImage(preds_ref_v4.2b, colors = color_pal, opacity = 1, group = "v4.2b: baseline") %>%
  addRasterImage(preds_ref_v4.2c, colors = color_pal, opacity = 1, group = "v4.2c: baseline") %>%
  addRasterImage(preds_ref_v4.2d, colors = color_pal, opacity = 1, group = "v4.2d: baseline") %>%
  addRasterImage(bgcs_3857, colors = color_pal, opacity = 1, group = "BGC") %>%
  addPolygons(data = gap_clipped, color = "black", weight = 1, group = "Holdouts") %>%
  addLayersControl(
    baseGroups = c("Satellite", "OSM Basemap"),  # Base layer switcher
    overlayGroups = c("v4.2: baseline",
                      "v4.2a: baseline",
                      "v4.2b: baseline",
                      "v4.2c: baseline",
                      "v4.2d: baseline",
                      "BGC", 
                      "Holdouts"),
    options = layersControlOptions(collapsed = FALSE)
  )
print(map)

#check alignment of gaps and training points
plot(trainingSample_V4a[, .(lon, lat)], pch=16, cex=0.1, xlim=c(-129.1, -128.3), ylim=c(54.1, 54.5))
plot(gap_sf, add=T)
