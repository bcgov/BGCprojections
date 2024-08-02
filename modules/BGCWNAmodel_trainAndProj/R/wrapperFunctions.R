#' Train `mlr3` models
#'
#' @param lrner an `[mlr3::Learner]` object. The model algorithm to fit.
#' @param task an `[mlr3::Task]` object. The data and response/predictor variables 
#'   to use.
#' @param modname character. Name of `lrner`. Passed to `reproducible::Cache(..., userTags)
#' @param cacheTags character. A vector of tags to pass to `reproducible::Cache(..., userTags)`
#'
#' @return `lrner` with the fitted models (which can be accessed with `lrner$model`)
#' @export
trainLearners <- function(lrner, task, modname = NULL, cacheTags = NULL) {
  for (i in 1:3) gc(reset = TRUE)
  
  modType <- capture.output(lrner$print())[1]
  
  cacheExtra <- list(
    modType,
    lrner$task_type, 
    lrner$param_set$values, 
    summary(task$data())
  )
  lrner$train(task) |>
    Cache(userTags = c(cacheTags, modname),
          .cacheExtra = cacheExtra,
          omitArgs = c("task", "userTags")) 
}

#' Evaluate trained `mlr3` models
#'
#' @inheritParams trainLearners 
#' @param strategy an `mlr3::Resampling` object. The model
#'  evaluation strategy passed to `mlr3::resample(..., resampling)`
#' @param nthreads numeric. Number of cores to parallelise over.
#'   Passed to `future::plan("multisession", workers = nthreads)`.
#'
#' @return
#' @export
evalLearners <- function(lrner, task, strategy, 
                         modname = NULL, nthreads = 1, cacheTags = NULL) {
  for (i in 1:3) gc(reset = TRUE)
  
  modType <- capture.output(lrner$print())[1]
  cvType <-  capture.output(strategy$print())[1]
  folds <- strategy$param_set$values$folds
  cacheExtra <- list(
      modType,
      lrner$task_type, 
      lrner$param_set$values,
      cvType,
      strategy$param_set$values,
      summary(task$data())
    )
  
  ## Note: Caching the outputs of `resample` will mean that we will not repeat
  ## a resampling if the inputs do not change. Meaning that stochasticity arising 
  ## from resampling during validation is not being accounted for. Presumably this is 
  ## fine because we don't want to run this multiple times.
  future::plan("multisession",
               workers = nthreads)
  out <- mlr3::resample(task, lrner, strategy, 
                        store_models = TRUE, store_backends = TRUE) |>
    Cache(userTags = c(modname, "cross-validation", cacheTags),
          .cacheExtra = cacheExtra,
          omitArgs = c("task", "learner", "resampling")) 
  future:::ClusterRegistry("stop")
  
  return(list("BGCmodel" = lrner, "BGCmodel_val" = out))
}