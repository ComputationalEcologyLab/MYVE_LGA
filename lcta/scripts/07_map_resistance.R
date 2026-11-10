
################################################################################
#
# Script description: Create resistance surface by predicting from fitted 
#                     gradient boosting machine
#
# Script author: Eric Palm
#
################################################################################

# Load required packages
sapply(c("dplyr", "stringr", "caret", "gbm", "terra"),
       require, character.only = T)

## ----------------------------------------------------------------------
## Best LCP iteration (1:5) by pooled RMSE
## ----------------------------------------------------------------------
source("lcta/scripts/cv_diagnostics.R")

## If running with a subset of variables (excluding karst vars and eco_niche)
## If not, comment this line out
source("lcta/scripts/run_config.R")


model_dir   <- paste0("lcta/models/", run_suffix)

find_best_iteration <- function(model_dir) {
  files <- list.files(model_dir, pattern = "^lcp_[1-9][0-9]*\\.rds$", full.names = TRUE)
  stopifnot(
    "no LCP iteration models found (lcp_1.rds..lcp_5.rds) -- run 04_lcp_iteration_cycle.R first" =
      length(files) > 0
  )
  
  results <- files %>%
    purrr::map_dfr(function(f) {
      fit  <- readRDS(f)
      rmse <- compute_fold_metrics(dplyr::semi_join(fit$pred, fit$bestTune))$global$rmse
      iter <- as.integer(stringr::str_match(basename(f), "^lcp_(\\d+)\\.rds$")[, 2])
      tibble::tibble(iteration = iter, rmse = rmse, file = f)
    })
  
  best <- results %>%
    dplyr::arrange(rmse, iteration) %>%
    dplyr::slice(1)
  
  stopifnot(
    "find_best_iteration() produced 0 rows -- every candidate file failed to score an RMSE" =
      nrow(best) == 1,
    "find_best_iteration()'s chosen file path is blank/NA" =
      is.character(best$file) && length(best$file) == 1 && !is.na(best$file) && nzchar(best$file),
    "find_best_iteration()'s chosen file doesn't exist on disk (stale/moved/deleted since list.files() ran)" =
      file.exists(best$file)
  )
  
  best
}

best  <- find_best_iteration(model_dir)
model <- readRDS(best$file)
message("Best LCP iteration: ", best$iteration, " (RMSE = ", round(best$rmse, 4), ")")


# Use the median geographic distance as a constant value for spatial predictions
env <- 
  list.files("data/SpatialLayers", full.names = T, pattern = ".tif$") %>% 
  terra::rast() %>% 
  `names<-`(c("crop_binary", "d_crop", "d_karst", "d_water",
              "elev", "eco_niche", "forest_prop_250", "forest_prop_gcam", "human_1", "human_2",
              "karst_prop", "lights", "slope", "tpi", "tri", "water_prop"
  )) %>% 
  {
    if (any(str_detect(colnames(model$trainingData), "\\beuc_geog\\b"))) {
      
      c(., .$slope %>% terra::setValues(., median(model$trainingData$euc_geog)) %>% 
       `names<-`("euc_geog"))
    } else {
      .
    }
  }

# Make spatial predictions of landscape resistance from the model using the covariate raster stack
pred_rast <- predict(env[[head(colnames(model$trainingData), -1)]], model) %>% 
  terra::mask(., env[["slope"]])

# Plot to make sure landscape resistance predictions look reasonable
terra::plot(pred_rast)

# Save raster to file
map_dir <- "lcta/outputs/maps"

dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

pred_rast %>% 
  terra::writeRaster(., stringr::str_c(map_dir, "/resistance_map", run_suffix, ".tif"), 
                     datatype = "FLT4S", overwrite = T)
