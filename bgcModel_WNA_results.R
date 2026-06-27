## Results of the BGC modeling

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
library(MASS)
library(scales)



## -------------------------------------------------
## bgc metadata

bgcs_info <- fread("https://github.com/bcgov/ccissr/raw/refs/heads/feas_tables/tables/WNA_BGCs_Info.csv")[,-1]
names(bgcs_info) <- as.character(bgcs_info[1,])
bgcs_info <- bgcs_info[-1,]

## -------------------------------------------------
## RF model

# dir <- "//objectstore2.nrs.bcgov/ffec/BGC_models/" 
dir <- "C:/Users/CMAHONY/Data/BGC_models/" #local copy, for speed

BGCmodel <- readRDS(paste0(dir, "BGCmodel_WNA_V4.2gini.rds")) #Kiri trained this on Thufir using the Gini split rule. 

studyname <- "BC"
bdy <- vect("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/bdy.BC.shp")
bdy <- project(bdy, "EPSG:4326")
dem <- rast("C:/Users/CMAHONY/OneDrive - Government of BC/Data/DEM/dem_noram_5arcsec.tif")
dem <- crop(dem, ext(bdy) )

# reduce resolution selecting center cell rather than taking an average, to avoid elevation error
fact <- 5
center <- (fact^2 + 1) / 2
dem <- aggregate(dem, fact = fact, fun = function(x, ...) x[center])

dem <- mask(dem, bdy) 
# plot(dem)
X <- dem # template raster for testing
values(X) <- NA

bgcs <- rast("C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/ccissv13_workingfiles/BGC_modelling/WNA_BGC_v13_13Apr2026_raster.tif")
bgcs <- crop(bgcs, ext(bdy))
bgcs <- project(bgcs, dem, method = "near")
bgcs <- mask(bgcs, dem) 
bgc_levels <- cats(bgcs)[[1]]  # Extract the category mapping
vals_char <- bgc_levels$BGC[match(values(bgcs), bgc_levels$value)] # Convert numeric values to character labels
values(bgcs) <- factor(vals_char, levels = bgcs_info$BGC) # refactor with full set of levels

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
clim_ref[!is.finite(Eref_sp), Eref_sp := 0]

clim_proj <- downscale(
  xyz = points,
  gcms = list_gcms()[5],
  ssps = list_ssps()[2],
  gcm_periods = list_gcm_periods()[3],
  run_nm = list_runs_ssp(list_gcms()[5], list_ssps()[2])[4],
  which_refmap = "refmap_climr",
  return_refperiod = FALSE, # Also return the 1961-1990 normals period.
  vars = list_vars()
)
ccissr::addVars(clim_proj)
clim_proj[!is.finite(CMD.total), CMD.total := 0]
clim_proj[!is.finite(Eref_sp), Eref_sp := 0]

## -------------------------------------------------
## predictions

## V4.2gini model predictions for reference period
preds_ref_vec <- predict(BGCmodel, data = clim_ref)$prediction
preds_ref <- X
preds_ref[points$id] <- factor(preds_ref_vec, levels = bgcs_info$BGC)

## V4.2gini model predictions for future period
preds_proj_vec <- predict(BGCmodel, data = clim_proj)$prediction
preds_proj <- X
preds_proj[points$id] <- factor(preds_proj_vec, levels = subzones_colours_ref$BGC)


## -------------------------------------------------
## -------------------------------------------------
## bivariate kernel density distributions of baseline and future latitude and elevation for each BGC zone
## -------------------------------------------------
## -------------------------------------------------

zone.ref <- bgcs_info[match(preds_ref_vec, bgcs_info$BGC), Zone]

zone.proj <- bgcs_info[match(preds_proj_vec, bgcs_info$BGC), Zone]

zones <- unique(bgcs_info[DataSet=="BC", Zone ])

figdir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(figdir, "/ElevShift.zones.png"), type="cairo", units="in", width=6.5, height=6.5, pointsize=10, res=300)

mat <- matrix(c(17, 1:4, 17, 5:8, 17, 9:12, 17, 13:16, rep(18,5)),5, byrow=T)   #define the plotting order
layout(mat, widths=c(0.1,1,1,1,1), heights=c(1,1,1,1,0.1))   #set up the multipanel plot
par(mar=c(1.5,1.5,0.8,0.8), mgp=c(1.5,0.25,0), tck= -0.01)

for(zone in zones){
  elev.ref <- points[which(zone.ref==zone), elev]
  lat.ref <- points[which(zone.ref==zone), lat]
  elev.proj <- points[which(zone.proj==zone), elev]
  lat.proj <- points[which(zone.proj==zone), lat]
  
  x <- c(elev.ref, elev.proj)
  y <- c(lat.ref, lat.proj)
  
  # KDE parameters
  h.factor <- 0.1
  prob.threshold <- 0.95 #probability contour to draw as polygon
  
  # Kernel density estimation for ref period
  s <- sample(1:length(elev.ref), 10000, replace = TRUE) # sample to reduce computation
  k <- kde2d(elev.ref[s], lat.ref[s], n=500, h=c(h.factor*diff(range(x)), h.factor*diff(range(y))), lims = c(range(x)+c(-200, 200), range(y)+c(-2, 20)))
  dx <- diff(k$x[1:2]); dy <- diff(k$y[1:2]) # grid cell area
  z <- k$z / sum(k$z * dx * dy) # normalize density surface
  zvec <- as.vector(z); # flatten
  zsort <- sort(zvec, decreasing = TRUE) # sort
  cumprob <- cumsum(zsort * dx * dy) # cumulative probability mass
  lev <- zsort[min(which(cumprob >= prob.threshold))] # threshold for specified probability threshold HDR
  cl.ref <- contourLines(x = k$x, y = k$y, z = k$z, levels = lev)

  # Kernel density estimation for proj period
  s <- sample(1:length(elev.proj), 10000, replace = TRUE)
  k <- kde2d(elev.proj[s], lat.proj[s], n=500, h=c(h.factor*diff(range(x)), h.factor*diff(range(y))), lims = c(range(x)+c(-200, 200), range(y)+c(-2, 20)))
  dx <- diff(k$x[1:2]); dy <- diff(k$y[1:2]) # grid cell area
  z <- k$z / sum(k$z * dx * dy) # normalize density surface
  zvec <- as.vector(z); # flatten
  zsort <- sort(zvec, decreasing = TRUE) # sort
  cumprob <- cumsum(zsort * dx * dy) # cumulative probability mass
  lev <- zsort[min(which(cumprob >= prob.threshold))] # threshold for specified probability threshold HDR
  cl.proj <- contourLines(x = k$x, y = k$y, z = k$z, levels = lev)

  # all x coordinates
  xrange <- range(c(
    unlist(lapply(cl.ref, `[[`, "x")),
    unlist(lapply(cl.proj, `[[`, "x"))
  ))
  
  # all y coordinates
  yrange <- range(c(
    unlist(lapply(cl.ref, `[[`, "y")),
    unlist(lapply(cl.proj, `[[`, "y"))
  ))
  
  # plot(x,y, xaxs="i", yaxs="i", col="white")
  plot(1, type = "n", xlab = "", ylab="", xlim = range(xrange), ylim = range(yrange))
  
  # plot polygons
  for(i in 1:length(cl.ref)){polygon(cl.ref[[i]]$x, cl.ref[[i]]$y, border = "grey", lwd = 1, col = adjustcolor("grey", alpha.f = 0.2))}   # draw polygon
  for(i in 1:length(cl.proj)){polygon(cl.proj[[i]]$x, cl.proj[[i]]$y, border = "dodgerblue", lwd = 1, col = adjustcolor("dodgerblue", alpha.f = 0.2))}   # draw polygon
  # s <- sample(1:length(elev.ref), 10000)
  # points(elev.ref[s],lat.ref[s], col=alpha("dodgerblue", 0.25), pch=16, cex=0.5)
  
  mtext(paste0("(", letters[which(zones==zone)], ") ", zone), line=-1.5, side=3, adj=0.025, cex=0.8, font=1)
  
  if(zone==zones[1]){
    legend("bottomleft", legend=c("Baseline (1961-1990)", "Future (2041-2060)"), fill=alpha(c("grey", "dodgerblue"), 0.2), border=c("grey", "dodgerblue"), bty="n")
  }
  
  # plot(x,y, col="white", xlab = "Elevation (m)", ylab="Latitude (deg.)")
  # points(elev.ref,lat.ref, col=alpha("dodgerblue", 0.25), pch=16, cex=1.5)
  # points(elev.proj,lat.proj, col=alpha("grey40", 0.85), cex=1.5)
  
print(zone)
  
}
par(mar=c(0,0,0,0))

plot(1, type="n", axes=F, xlab="", ylab="")  
text(1,1, "Latitude (degrees)", srt=90, font=1,cex=1.3)  

par(mar=c(0,0,0,0))
plot(1, type="n", axes=F, xlab="", ylab="")  
text(1,1, "Elevation (m)", font=1,cex=1.3)  

dev.off()


## -------------------------------------------------
## -------------------------------------------------
## error analyses
## -------------------------------------------------
## -------------------------------------------------

# -------------------------------------------------
## total error

results.error <- mean(values(bgcs, mat = FALSE) != values(preds_ref, mat = FALSE), na.rm=T)

## -------------------------------------------------
## error summary (chatgpt)
## -------------------------------------------------

ref <- factor(
  values(bgcs, mat = FALSE),
  levels = seq_along(bgcs_info$BGC),
  labels = bgcs_info$BGC
)

pred <- factor(
  values(preds_ref, mat = FALSE),
  levels = seq_along(bgcs_info$BGC),
  labels = bgcs_info$BGC
)

# remove NA pairs
keep <- !is.na(ref) & !is.na(pred)

ref  <- ref[keep]
pred <- pred[keep]

# confusion matrix
cm <- table(
  Reference = ref,
  Predicted = pred
)

cm

# Category-wise metrics:

# total samples
n <- sum(cm)

# true positives
tp <- diag(cm)

# false negatives (omission error)
fn <- rowSums(cm) - tp

# false positives (commission error)
fp <- colSums(cm) - tp

# producer's accuracy (recall)
producer_acc <- tp / rowSums(cm)

# user's accuracy (precision)
user_acc <- tp / colSums(cm)

# omission error
omission_err <- fn/rowSums(cm)

# commission error
commission_err <- fp/colSums(cm)

# area bias ratio
area_bias <- colSums(cm) / rowSums(cm) 

# per-class accuracy summary
acc <- data.frame(
  class = rownames(cm),
  n_ref = rowSums(cm),
  n_pred = colSums(cm),
  correct = tp,
  producer_accuracy = producer_acc,
  user_accuracy = user_acc,
  omission_error = omission_err,
  commission_error = commission_err,
  area_bias = area_bias
)

acc <- acc[acc$class %in% bgcs_info[DataSet=="BC", BGC],]
acc <- acc[is.finite(acc$producer_acc),]
acc$zone <- bgcs_info[match(acc$class, bgcs_info$BGC), Zone]

setDT(acc)

# training sample
trainingSample <- fread(file = paste0(dir, "points_WNA_v4a.csv"))
sampSize <- as.data.table(table(trainingSample$BGC))
setnames(sampSize, "V1", "BGC")
acc[, training_sample := sampSize[.SD, on = .(BGC = class), N]]


## -------------------------------------------------
## diagnostic plots of baseline prediction error 
## -------------------------------------------------

figdir <- "C:/Users/CMAHONY/Government of BC/Future Forest Ecosystems Centre - CCISS - CCISS/CCISS_Manual/Figures/"
png(filename=paste0(figdir, "BaselinePrediction_ErrorDiagnostics.png",sep="."), type="cairo", units="in", width=6.5, height=6.5, pointsize=12, res=300)

par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.75,0.25, 0), tck=-0.005)
z <- acc$class

## -------------------------------------------------
# plot of baseline mismatch vs BGC unit area

cellsize <- res(dem)[1]*111.32*res(dem)[1]*111.32*cos(54.5*pi/180)
x <- log2(acc$n_ref*cellsize)
y <- 1-acc$producer_accuracy

plot(x,y, col="white", xlim=range(x)*c(0.96, 1.03), yaxt="n", xaxt="n", ylab="Baseline prediction error", xlab=expression("BGC unit area (" * km^2 * ")"))
# text(x,y,z, cex=0.5, font=2, col=bgcs_info[match(acc$class, bgcs_info$BGC), ZoneColour])
text(x,y,z, cex=0.5)
axis(1, at = seq(0,99,2), labels = 2^(seq(0,99,2)))
axis(2, at = pretty(y), labels = pretty(y), las=2)
mtext("(a)", side = 3, line=-1.25, adj=0.025, font=1)
mtext(paste0("r = ", round(cor(x,y), 2)), side = 3, line=-1.5, adj=0.975, font=1)

## -------------------------------------------------
# plot of omission error vs commission error

x <- acc$commission_error
y <- acc$omission_error

plot(x,y, col="white", yaxt="n", xlim=range(x)*c(0.96, 1.03), ylab="Omission error", xlab="Commission error")
lines(c(-99,99), c(-99,99), col="gray", lty=2)
text(x,y,z, cex=0.5)
axis(2, at = pretty(y), labels = pretty(y), las=2)
mtext("(b)", side = 3, line=-1.25, adj=0.025, font=1)
mtext(paste0("r = ", round(cor(x,y), 2)), side = 1, line=-1.5, adj=0.975, font=1)

## -------------------------------------------------
# plot of area bias vs BGC unit area

cellsize <- res(dem)[1]*111.32*res(dem)[1]*111.32*cos(54.5*pi/180)
x <- log2(acc$n_ref*cellsize)
y <- log2(acc$area_bias)

plot(x,y, col="white", xlim=range(x)*c(0.96, 1.03), yaxt="n", xaxt="n", ylab="Area bias ratio", xlab=expression("BGC unit area (" * km^2 * ")"))
# text(x,y,z, cex=0.5, font=2, col=bgcs_info[match(acc$class, bgcs_info$BGC), ZoneColour])
lines(c(1,99), c(0,0), col="gray", lty=2)
text(x,y,z, cex=0.5)
axis(1, at = seq(0,99,2), labels = 2^(seq(0,99,2)))
axis(2, at = log2(seq(0.2, 2, 0.1)), labels = seq(0.2, 2, 0.1), las=2)
mtext("(c)", side = 3, line=-1.25, adj=0.025, font=1)
mtext(paste0("r = ", round(cor(x,y), 2)), side = 1, line=-1.5, adj=0.975, font=1)

## -------------------------------------------------
# plot of area bias vs training sample

x <- log2(acc$training_sample)
y <- log2(acc$area_bias)

plot(x,y, col="white", xlim=range(x)*c(0.96, 1.03), yaxt="n", xaxt="n", ylab="Area bias ratio", xlab="Training sample size")
# text(x,y,z, cex=0.5, font=2, col=bgcs_info[match(acc$class, bgcs_info$BGC), ZoneColour])
lines(c(1,99), c(0,0), col="gray", lty=2)
text(x,y,z, cex=0.5)
axis(1, at = seq(1,99,2), labels = 2^(seq(1,99,2)))
axis(2, at = log2(seq(0.2, 2, 0.1)), labels = seq(0.2, 2, 0.1), las=2)
mtext("(d)", side = 3, line=-1.25, adj=0.025, font=1)
mtext(paste0("r = ", round(cor(x,y), 2)), side = 1, line=-1.5, adj=0.975, font=1)

dev.off()



## -------------------------------------------------
## zone-level error summary 
## -------------------------------------------------

# zone error 
zone_err <- acc[
  ,
  .(
    weighted_omission_error =
      weighted.mean(omission_error, w = n_ref, na.rm = TRUE),
    
    weighted_producer_accuracy =
      weighted.mean(producer_accuracy, w = n_ref, na.rm = TRUE),
    
    total_n_ref = sum(n_ref),
    mean_n_ref = mean(n_ref)
  ),
  by = zone
]

zone_err

ref.zone <- bgcs_info[match(ref, bgcs_info$BGC), Zone]
  
pred.zone <- bgcs_info[match(pred, bgcs_info$BGC), Zone]

# remove NA pairs
keep <- !is.na(ref.zone) & !is.na(pred.zone)

ref.zone  <- ref.zone[keep]
pred.zone <- pred.zone[keep]

# confusion matrix
cm <- table(
  Reference = ref.zone,
  Predicted = pred.zone
)

cm

# Category-wise metrics:

# total samples
n <- sum(cm)

# true positives
tp <- diag(cm)

# false negatives (omission error)
fn <- rowSums(cm) - tp

# false positives (commission error)
fp <- colSums(cm) - tp

# producer's accuracy (recall)
producer_acc <- tp / rowSums(cm)

# user's accuracy (precision)
user_acc <- tp / colSums(cm)

# omission error
omission_err <- 1 - producer_acc

# commission error
commission_err <- 1 - user_acc

# per-class accuracy summary
acc.zone <- data.frame(
  class = rownames(cm),
  n_ref = rowSums(cm),
  n_pred = colSums(cm),
  correct = tp,
  producer_accuracy = producer_acc,
  user_accuracy = user_acc,
  omission_error = omission_err,
  commission_error = commission_err
)

acc.zone <- acc.zone[!is.na(producer_acc),]

write.csv(acc.zone, "acc.zone.csv", row.names = FALSE)

