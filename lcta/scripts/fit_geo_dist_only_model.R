
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Fit geographic distance-only model
#
# Script author: Eric Palm
#
################################################################################

# Fit a straight-line geographic distance-only model to provide basis of 
# comparison for RMSEtest value from top model comparison

# Load required packages
sapply(c("doParallel", "dplyr", "tibble"),
       require, character.only = T)

# Load required scripts
source("lcta/scripts/fit_tools.R")     
source("lcta/scripts/cv_diagnostics.R")   

## ---- User inputs ----
run_suffix   <- "_subset"
model_dir    <- paste0("models", run_suffix, "/")
response_col <- "gd"
noCLS        <- 60
## ---------------------------------------------------------------

folds <- readRDS("lcta/results/bat_cv_folds.rds")
sub_data <- readRDS("lcta/results/straight_line_extraction.rds")
stopifnot(nrow(sub_data) == folds$n_pairs)

cl <- parallel::makeCluster(noCLS)
doParallel::registerDoParallel(cl)
fit_dist <- tryCatch(
  fit_gbm_cv(sub_data, "euc_geog", response_col, folds),
  finally = { parallel::stopCluster(cl); foreach::registerDoSEQ() }
)

# Calculate pooled RMSEtest
p <- dplyr::semi_join(fit_dist$pred, fit_dist$bestTune, by = names(fit_dist$bestTune)) %>%
  dplyr::filter(!is.na(obs), !is.na(pred))

list(rmse = compute_fold_metrics(dplyr::semi_join(fit_dist$pred, fit_dist$bestTune))$global$rmse)

