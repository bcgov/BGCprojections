# Load packages: 
library(tidyverse)
library(terra)
library(climr)
library(reproducible) # For Cache function
library(data.table)
library(sf)
library(foreach) # for outlier removal function
library(tidymodels) # for prep() function from recipes package. 
library(themis) # for step_downsample() function
library(ranger) # For RF
library(caret) # For confusionMatrix()


# Read in training data: 
trainData <- read.csv("data-generated/trainData-no-removals.csv")

# Assign variable sets (described more in 01-establish-training-points): 
vars_simple <- c("PPT_sp", "PPT_sm", "PPT_at", "PPT_wt", 
                 "Tmax_sp", "Tmax_sm", "Tmax_at", "Tmax_wt", 
                 "Tmin_sp", "Tmin_sm", "Tmin_at", "Tmin_wt")

vars_LFS <- c("")

vars_expert <- c("CMD_sm", "DDsub0_sp", "DD5_sp", "Eref_sm", "Eref_sp", 
                 "EXT", "MWMT", "NFFD_sm", "NFFD_sp", "PAS", "PAS_sp", 
                 "SHM", "Tave_sm", "Tave_sp", "Tmax_sm", "Tmax_sp", "Tmin", 
                 "Tmin_at", "Tmin_sm", "Tmin_sp", "Tmin_wt", "CMI", "PPT_MJ", 
                 "PPT_JAS", "CMD.total")

vars_all <- c(list_vars(), "PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total", "DD_delayed")

#### Train ranger random forest model: ####
trainData <- as.data.table(trainData)
trainData[, BGC := as.factor(BGC)]
trainData_Wgaps <- trainData[gap == "no"]

cols_simple <- c("BGC", vars_simple)
cols_expert <- c("BGC", vars_expert)
cols_all <- c("BGC", vars_all)

# Train model with simple variables, on points from the entire study area: 
# BGCmodel_full_simple <- ranger(
#   BGC ~ .,
#   data = trainData[, ..cols_simple],
#   num.trees = 501,
#   splitrule =  "extratrees",
#   # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
#   importance = "permutation",
#   write.forest = TRUE,
#   classification = TRUE,
#   probability = FALSE,
#   
# ) |>
#   Cache()
# 
# # Note: V1 = no BGCs removed. Saving these locally. 
# saveRDS(BGCmodel_full_simple, "RF_models/BGCmodel_full_simple_V1.rds")
BGCmodel_full_simple <- readRDS("RF_models/BGCmodel_full_simple_V1.rds")

# Train model with simple variables, on points from area outside of the gaps: 
BGCmodel_Wgaps_simple <- ranger(
  BGC ~ .,
  data = trainData_Wgaps[, ..cols_simple],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2,
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
) |>
  Cache()

# Note: V1 = no BGCs removed.
saveRDS(BGCmodel_Wgaps_simple, "RF_models/BGCmodel_Wgaps_simple_V1.rds")
BGCmodel_Wgaps_simple <- readRDS("RF_models/BGCmodel_Wgaps_simple_V1.rds")

# Train model with expert selection of variables, on points from the entire study area: 
# BGCmodel_full_expert <- ranger(
#   BGC ~ .,
#   data = trainData[, ..cols_expert],
#   num.trees = 501,
#   splitrule =  "extratrees",
#   # min.node.size = 2,
#   importance = "permutation",
#   write.forest = TRUE,
#   classification = TRUE,
#   probability = FALSE,
# 
# ) |>
#   Cache()
# 
# # Note: V1 = no BGCs removed.
# saveRDS(BGCmodel_full_expert, "RF_models/BGCmodel_full_expert_V1.rds")
BGCmodel_full_expert <- readRDS("RF_models/BGCmodel_full_simple_V1.rds")

BGCmodel_Wgaps_expert <- ranger(
  BGC ~ .,
  data = trainData_Wgaps[, ..cols_expert],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, 
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
) |>
  Cache()

# Note: V1 = no BGCs removed.
saveRDS(BGCmodel_Wgaps_expert, "RF_models/BGCmodel_Wgaps_expert_V1.rds")
BGCmodel_Wgaps_expert <- readRDS("RF_models/BGCmodel_Wgaps_expert_V1.rds")

# Train model on all climate variables (kitchen sink scenario): 
BGCmodel_full_all <- ranger(
  BGC ~ .,
  data = trainData[, ..cols_all],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, 
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
) |>
  Cache()

# Note: V1 = no BGCs removed.
saveRDS(BGCmodel_full_all, "RF_models/BGCmodel_full_all_V1.rds")
BGCmodel_full_all <- readRDS("RF_models/BGCmodel_full_all_V1.rds")

BGCmodel_Wgaps_all <- ranger(
  BGC ~ .,
  data = trainData_Wgaps[, ..cols_all],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, 
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
) |>
  Cache()
beepr::beep()

# Note: V1 = no BGCs removed.
saveRDS(BGCmodel_Wgaps_all, "RF_models/BGCmodel_Wgaps_all_V1.rds")
BGCmodel_Wgaps_all <- readRDS("RF_models/BGCmodel_Wgaps_all_V1.rds")

# Check the models: 
conf_matrix_full_simple <- caret::confusionMatrix(data = predictions(BGCmodel_full_simple),
                                                  reference = trainData$BGC)
conf_matrix_Wgaps_simple <- caret::confusionMatrix(data = predictions(BGCmodel_Wgaps_simple),
                                                   reference = trainData_Wgaps$BGC)

conf_matrix_full_expert <- caret::confusionMatrix(data = predictions(BGCmodel_full_expert),
                                                  reference = trainData$BGC)
conf_matrix_Wgaps_expert <- caret::confusionMatrix(data = predictions(BGCmodel_Wgaps_expert),
                                                   reference = trainData_Wgaps$BGC)

conf_matrix_full_all <- caret::confusionMatrix(data = predictions(BGCmodel_full_all),
                                               reference = trainData$BGC)
conf_matrix_Wgaps_all <- caret::confusionMatrix(data = predictions(BGCmodel_Wgaps_all),
                                                reference = trainData_Wgaps$BGC)

print(BGCmodel_full_simple) # OOB prediction error: 21.87%
print(BGCmodel_Wgaps_simple) # OOB prediction error: 21.91%
print(BGCmodel_full_expert) # OOB prediction error: 22.99%
print(BGCmodel_Wgaps_expert) # OOB prediction error: 23.02%
print(BGCmodel_full_all) # OOB prediction error: 20.44%
print(BGCmodel_Wgaps_all) # OOB prediction error: 20.51%