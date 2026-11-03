
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Reselect variables using LCP iteration
#
# Script author: Eric Palm
#
################################################################################

## ============================================================================
## Reselect variables using LCP-refined covariates (from a completed LCP
## iteration you've already run) instead of straight-line ones, to check
## whether the straight-line-based selection in 02_select_variables.R
## (results/bat_variable_selection.rds) still holds using LCPs. 
## Set ITERATION below depending on which round you want to check against.
##
## Uses the LCP LINES ALREADY ON DISK from 04_lcp_iteration_cycle.R
## (outputs/lcp_lines/lcp_<ITERATION>.rds) -- no LCP computation happens
## here, this only re-extracts covariates along those existing lines and
## reruns the same forward search (select_variables_forced() from
## variable_selection_tools.R, unchanged) against that table instead of the
## straight-line one.
##
## Does NOT touch results/bat_variable_selection.rds -- saves to
## results/bat_variable_selection_lcp<ITERATION>.rds instead, and prints a
## comparison against the straight-line result so you can see exactly
## what would change before deciding whether to adopt it.
##
## SELECTION_MODE -- "standard" or "robust": same toggle as
## 02_select_variables.R. "standard" (default) runs
## select_variables_forced() at one fixed hyperparameter setting.
## "robust" runs select_variables_robustness_check() (3 fixed settings)
## and keeps a candidate if >= robust_majority_threshold of them selected it.
##
## COMPARE_TUNED_RMSE: select_variables_forced()'s RMSE comes from
## evaluate_variable_set()'s modest hyperparameter settings, which are NOT
## comparable to the fully-tuned RMSE saved on the existing iteration-
## <ITERATION> model on disk (from fit_gbm_cv()'s full tuneGrid search in
## 03/04). When the reselected set differs from what's currently locked
## in (results/bat_variable_selection.rds), this fits a FULLY TUNED model
## (fit_gbm_cv(), same tuneGrid as 03/04) using the reselected variables
## on the SAME LCP lines already used for reselection above, and compares
## its pooled RMSE against:
##   - the existing tuned iteration-<ITERATION> model (whatever variable
##     set is currently locked in), AND
##   - the existing tuned iteration-<ITERATION - 1> model (the previous
##     round).
## worth_adopting requires beating BOTH. Adopting the reselected set would 
## generate different LCP lines on the NEXT iteration going forward, 
## so a win here means "worth adopting and rerunning 03/04 to see if it holds 
## up," not "already proven better."
## ============================================================================

# Load packages
sapply(c("dplyr", "parallel", "foreach", "tibble"), 
       require, 
       character.only = T)

source("lcta/scripts/01_extract_covariates.R")
source("lcta/scripts/variable_selection_tools.R")
source("lcta/scripts/fit_tools.R")    # fit_gbm_cv(), for the tuned RMSE comparison
source("lcta/scripts/cv_diagnostics.R") # compute_fold_metrics() to compare tuned RMSE 

## If running with a subset of variables (excluding karst vars and eco_niche)
## If not, comment this line out
source("lcta/scripts/run_config.R")

## ---- User configure ----
response_col   <- "gd"                                
forced_vars    <- "euc_geog"                           
candidate_vars <- setdiff(var_names, drop_vars) # every raster-derived covariate

## Same mutually-exclusive variable groups as 02_select_variables.R -- see
## its header note and select_variables_forced()'s exclusive_groups
## argument (selection_helpers.R) for what this enforces.
exclusive_var_groups <- list(
  c("crop_binary", "d_crop"),
  c("karst_prop", "d_karst"),
  c("human_1", "human_2", "lights"),
  c("forest_prop_250", "forest_prop_gcam"),
  c("slope", "tri"),
  c("water_prop", "d_water")
)

exclusive_var_groups <-
  exclusive_var_groups %>% 
  lapply(intersect, candidate_vars) %>% 
  purrr::keep(~ length(.x) > 1)

ITERATION      <- 2     # which completed LCP iteration to reselect
model_dir     <- paste0("lcta/models/", run_suffix, "/")     # straight-line will be iteration 0
lines_dir     <- paste0("lcta/outputs/lcp_lines/", run_suffix, "/")
noCLS          <- 60

SELECTION_MODE <- "robust"   # "standard" or "robust" -- see the header
                                # note above for what changes
robust_majority_threshold <- 2 # robust mode only: keep a candidate if
                                # selected by >= this many of the 3 fixed
                                # hyperparameter settings

COMPARE_TUNED_RMSE <- TRUE     # see the COMPARE_TUNED_RMSE header note --
                                # fair tuned-vs-tuned RMSE check, only
                                # runs if the reselection actually changed
                                # the variable set
noCLS_tuned_fit <- 60          # separate cluster for the tuned-fit
                                # comparison -- created only after the
                                # extraction cluster (noCLS above) is
                                # fully stopped, never both alive at once
## ---------------------------------------------------------------

lines_path <- function(iteration, dir = lines_dir) paste0(dir, "lcp_", iteration, ".rds")

## Matches 03_fit_straight_line.R/04_lcp_iteration_cycle.R's model_path()
model_path <- function(iteration, dir = model_dir) paste0(dir, "lcp_", iteration, ".rds")

folds              <- readRDS("lcta/results/bat_cv_folds.rds")
selection_straight <- readRDS(paste0("lcta/results/bat_variable_selection", run_suffix, ".rds"))
pairs <- readRDS("lcta/data/samples_pairwise.rds") %>%  
  dplyr::mutate(index = dplyr::row_number())

sub_pairs <- pairs %>% dplyr::filter(index %in% folds$orig_index)

lines_file <- lines_path(ITERATION)
if (!file.exists(lines_file)) {
  stop("No saved lines at ", lines_file, " -- run 04_lcp_iteration_cycle.R for iteration ",
       ITERATION, " first")
}

lcp_lines <- readRDS(lines_file)

cl <- parallel::makeCluster(noCLS)
doParallel::registerDoParallel(cl)

sub_data <- tryCatch(
  extract_covariates(lcp_lines, sub_pairs, already_buffered = TRUE),
  finally = { parallel::stopCluster(cl); foreach::registerDoSEQ() }
)
stopifnot(nrow(sub_data) == folds$n_pairs)

if (SELECTION_MODE == "robust") {

  robust <- select_variables_robustness_check(sub_data, folds, forced_vars, candidate_vars, response_col,
                                               exclusive_groups = exclusive_var_groups)

  message("Stability across ", length(default_gbm_param_settings), " hyperparameter settings:")
  print(robust$stability, n = Inf)

  keep <- robust$stability %>%
    dplyr::filter(n_settings_selected >= robust_majority_threshold) %>%
    dplyr::pull(variable)
  selected_vars <- c(forced_vars, candidate_vars[candidate_vars %in% keep])

  message("Majority-vote selection (>= ", robust_majority_threshold, "/",
          length(default_gbm_param_settings), " settings): ", paste(selected_vars, collapse = ", "))

  selection_lcp <- list(selected_vars = selected_vars, mode = "robust",
                         history = robust$results_by_setting$modest$history,
                         robust_results = robust$results_by_setting, stability = robust$stability)

} else {
  selection_lcp <- select_variables_forced(sub_data, folds, forced_vars, candidate_vars, response_col,
                                            exclusive_groups = exclusive_var_groups)
  selection_lcp <- c(selection_lcp, list(mode = "standard"))
}

dir.create("results", recursive = TRUE, showWarnings = FALSE)
saveRDS(selection_lcp, paste0("lcta/results/bat_variable_selection_lcp", ITERATION, run_suffix, ".rds"))

## ----------------------------------------------------------------------
## Comparison -- straight-line selection vs. LCP-N reselection
## ----------------------------------------------------------------------
straight_vars <- selection_straight$selected_vars
lcp_vars      <- selection_lcp$selected_vars
changed       <- !setequal(straight_vars, lcp_vars)

message("Straight-line: ", paste(straight_vars, collapse = ", "))
message("LCP", ITERATION, " (", SELECTION_MODE, "): ", paste(lcp_vars, collapse = ", "))
message("Changed: ", changed)

## ----------------------------------------------------------------------
## COMPARE_TUNED_RMSE -- Only runs if changed == TRUE; 
## fits the reselected variables with the SAME tuning procedure (fit_gbm_cv()) 
## as the existing models, on the SAME LCP lines used for reselection above, and
## compares pooled RMSE against BOTH the current iteration's official model and
## the previous iteration's.
## ----------------------------------------------------------------------
if (COMPARE_TUNED_RMSE && !changed) {
  message("\nCOMPARE_TUNED_RMSE: variable set didn't change -- nothing to compare.")
} else if (COMPARE_TUNED_RMSE) {

  message("\nCOMPARE_TUNED_RMSE: fitting a fully-tuned model for the reselected variable set",
          " against the existing tuned iteration-", ITERATION, " and iteration-", ITERATION - 1,
          " models...")

  fit_cl <- parallel::makeCluster(noCLS_tuned_fit)
  doParallel::registerDoParallel(fit_cl)

  ## Reads whatever official model is saved at a given iteration and
  ## returns its pooled RMSE
  official_rmse <- function(it) {
    path <- model_path(it)
    if (!file.exists(path)) {
      message("  No existing iteration-", it, " model at ", path, " -- treating as NA.")
      return(NA_real_)
    }
    fit <- readRDS(path)
    compute_fold_metrics(dplyr::semi_join(fit$pred, fit$bestTune))$global$rmse
  }

  tuned_comparison <- tryCatch({

    ## New fit: reselected variables, fully tuned, on these SAME LCP lines.
    new_fit  <- fit_gbm_cv(sub_data, lcp_vars, response_col, folds)
    new_rmse <- compute_fold_metrics(dplyr::semi_join(new_fit$pred, new_fit$bestTune))$global$rmse

    old_rmse  <- official_rmse(ITERATION)
    prev_rmse <- official_rmse(ITERATION - 1)

    tibble::tibble(
      prev_iter = ITERATION - 1, prev_rmse = prev_rmse,
      old_iter  = ITERATION,     old_rmse  = old_rmse,
      new_rmse  = new_rmse,
      new_beats_old  = new_rmse < old_rmse,
      new_beats_prev = new_rmse < prev_rmse,
      worth_adopting = (new_rmse < old_rmse) & (new_rmse < prev_rmse)
    )
  }, finally = {
    parallel::stopCluster(fit_cl)
    foreach::registerDoSEQ()
  })

  print(tuned_comparison, n = Inf)
  message("worth_adopting requires beating BOTH prev and old",
          " if old is itself a regression from prev (RMSE doesn't always improve LCP-over-LCP).")
}

## If the variable set changed and you want to adopt the LCP-N result
## going forward: overwrite results/bat_variable_selection.rds with
## selection_lcp, then rerun 03_fit_straight_line.R and
## 04_lcp_iteration_cycle.R so the locked models reflect the updated
## variable set -- neither script does this automatically, on purpose.
