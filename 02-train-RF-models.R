# Read in training data sets: ####

#### Train ranger random forest model: ####
trainData <- as.data.table(trainData)
trainData[, BGC := as.factor(BGC)]
trainData_Wgaps <- trainData[gap == "no"]

cols_simple <- c("BGC", vars_simple)
cols_expert <- c("BGC", vars_expert)
cols_all <- c("BGC", vars_all)

# Train model with simple variables, on points from the entire study area: 
BGCmodel_full_simple <- ranger(
  BGC ~ .,
  data = trainData[, ..cols_simple],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, # Default is 1. Try with 1, see if model overfits. All other code used 2 so maybe that was why.
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
  
) |>
  Cache()

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

# Train model with expert selection of variables, on points from the entire study area: 
BGCmodel_full_expert <- ranger(
  BGC ~ .,
  data = trainData[, ..cols_expert],
  num.trees = 501,
  splitrule =  "extratrees",
  # min.node.size = 2, 
  importance = "permutation",
  write.forest = TRUE,
  classification = TRUE,
  probability = FALSE,
  
) |>
  Cache()

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