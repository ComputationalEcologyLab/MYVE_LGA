
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Variable selection for gradient boosting machines
#
# Script author: Eric Palm
#
################################################################################

## ============================================================================
## Variable selection -- run once, ONLY on the iteration-0 (straight-line) data
##
## SELECTION_MODE -- "standard" or "robust":
##   "standard" (default): forward selection at ONE fixed hyperparameter
##     setting 

##   "robust": runs select_variables_robustness_check() (already defined
##     in variable_selection_tools.R) -- forward selection independently under 3
##     fixed hyperparameter settings (modest/generous/regularized; see
##     default_gbm_param_settings in selection_tools.R). A candidate is
##     kept in the final selected_vars if AT LEAST
##     robust_majority_threshold (default 2) of the 3 settings selected it. 
##
## ============================================================================

# Load packages
sapply(c("dplyr", "parallel", "doParallel", "foreach"), 
       require, 
       character.only = T)

source("lcta/scripts/01_extract_covariates.R")
source("lcta/scripts/selection_helpers.R")

## ---- User configure -- 
response_col   <- "gd"                 
forced_vars    <- "euc_geog"  # forced variable(s); character(0) for none
candidate_vars <- var_names             # every raster-derived covariate
                                        # used in 01_extract_covariates.R

## Variable pairs/groups that are different representations of the same
## underlying thing -- forward selection only ever keeps ONE member of
## each group (see exclusive_groups in select_variables_forced(),
## selection_helpers.R): once one member is added, the rest of its group
## is permanently dropped from consideration for the remainder of the run.
exclusive_var_groups <- list(
  c("crop_binary", "d_crop"),
  c("karst_prop", "d_karst"),
  c("human_1", "human_2", "lights"),
  c("forest_prop_250", "forest_prop_gcam"),
  c("slope", "tri")
)

noCLS          <- 60   

SELECTION_MODE <- "robust"   # "standard" or "robust" 

robust_majority_threshold <- 2 # robust mode only: keep a candidate if
                                # selected by >= this many of the 3 fixed
                                # hyperparameter settings
## ---------------------------------------------------------------

folds <- readRDS("lcta/results/bat_cv_folds.rds")
pairs <- readRDS("lcta/data/samples_pairwise.rds") %>%   # ADAPT: your master file path
  dplyr::mutate(index = dplyr::row_number())

sub_pairs <- pairs %>% dplyr::filter(index %in% folds$orig_index)
lines     <- lines_from_pairs(sub_pairs)

cl <- parallel::makeCluster(noCLS)
doParallel::registerDoParallel(cl)

sub_data <- tryCatch(extract_covariates(lines, sub_pairs),
                      finally = { parallel::stopCluster(cl); foreach::registerDoSEQ() })

## sub_data must stay row-aligned with folds$index/indexOut positions --
stopifnot(nrow(sub_data) == folds$n_pairs)

dir.create("results", recursive = TRUE, showWarnings = FALSE)
saveRDS(sub_data, "results/straight_line_extraction.rds")

if (SELECTION_MODE == "robust") {

  robust <- select_variables_robustness_check(sub_data, folds, forced_vars, candidate_vars, response_col,
                                               exclusive_groups = exclusive_var_groups)

  message("Stability across ", length(default_gbm_param_settings), " hyperparameter settings:")
  print(robust$stability, n = Inf)

  keep <- robust$stability %>%
    dplyr::filter(n_settings_selected >= robust_majority_threshold) %>%
    dplyr::pull(variable)
  ## preserve candidate_vars' original order rather than stability's
  ## descending-support order
  selected_vars <- c(forced_vars, candidate_vars[candidate_vars %in% keep])

  message("Majority-vote selection (>= ", robust_majority_threshold, "/",
          length(default_gbm_param_settings), " settings): ", paste(selected_vars, collapse = ", "))

  selection_results <- list(selected_vars = selected_vars, mode = "robust",
                             history = robust$results_by_setting$modest$history,
                             robust_results = robust$results_by_setting, stability = robust$stability)

} else {
  res <- select_variables_forced(sub_data, folds, forced_vars, candidate_vars, response_col,
                                  exclusive_groups = exclusive_var_groups)
  print(res$history, n = Inf)
  message("Selected: ", paste(res$selected_vars, collapse = ", "))
  selection_results <- c(res, list(mode = "standard"))
}

saveRDS(selection_results, "lcta/results/bat_variable_selection.rds")

