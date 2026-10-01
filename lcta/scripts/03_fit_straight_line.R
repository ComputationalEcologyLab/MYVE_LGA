
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Fit straight-line GBM
#
# Script author: Eric Palm
#
################################################################################

## ============================================================================
## Fit the iteration-0 (straight-line) GBM.
##
## Requires savePredictions = "final" in caret::trainControl()
##
## Reuses 02_select_variables.R's saved extraction
## (lcta/results/straight_line_extraction.rds) instead of re-extracting
## ============================================================================

# Load packages
sapply(c("dplyr", "parallel", "foreach", "doParallel"), 
       require, 
       character.only = T)

source("lcta/scripts/01_extract_covariates.R")
source("lcta/scripts/fit_tools.R")


## ---- User configuration -- 
response_col <- "gd"           
outdir       <- "lcta/models/"
noCLS        <- 60
## ---------------------------------------------------------------

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

folds             <- readRDS("lcta/results/bat_cv_folds.rds")
selection_results <- readRDS("lcta/results/bat_variable_selection.rds")
pairs <- readRDS("lcta/data/samples_pairwise.rds") %>%   # ADAPT: your master file path
  dplyr::mutate(index = dplyr::row_number())

predictor_cols <- selection_results$selected_vars

## iteration 0 = straight-line 
model_path <- function(iteration, dir = outdir) paste0(dir, "lcp_", iteration, ".rds")

cache_path <- "lcta/results/straight_line_extraction.rds"
sub_data <- if (file.exists(cache_path)) {
  readRDS(cache_path)
} else {
  sub_pairs <- pairs %>% dplyr::filter(index %in% folds$orig_index)
  lines <- lines_from_pairs(sub_pairs)
  cl <- parallel::makeCluster(noCLS)
  doParallel::registerDoParallel(cl)
  tryCatch(extract_covariates(lines, sub_pairs),
           finally = { parallel::stopCluster(cl); foreach::registerDoSEQ() })
}
stopifnot(nrow(sub_data) == folds$n_pairs)

fit <- fit_gbm_cv(sub_data, predictor_cols, response_col, folds)
saveRDS(fit, model_path(0))
