
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Variable selection tools
#
# Script author: Eric Palm
#
################################################################################

## ============================================================================
## Variable selection tools
## ============================================================================

# Load packages
sapply(c("dplyr", "purrr", "tibble", "gbm"),
       require, 
       character.only = T)

## ----------------------------------------------------------------------
## Fit + pool-evaluate one candidate variable set using spatial CV
## folds. Returns pooled RMSE across all held-out rows across all folds,
## as in  CAST::global_validation()
##
## Sets random seed to keep things reproducible.
## ----------------------------------------------------------------------
evaluate_variable_set <- function(sub_data, folds, vars, response_col,
                                   gbm_params = list(n.trees = 1000, interaction.depth = 1,
                                                      shrinkage = 0.01, n.minobsinnode = 10),
                                   seed = 1234) {

  set.seed(seed)

  fold_formula <- stats::reformulate(sprintf("`%s`", vars), response = sprintf("`%s`", response_col))

  fold_preds <- names(folds$index) %>%
    purrr::map_dfr(function(fold_name) {
      train_data <- sub_data[folds$index[[fold_name]], c(vars, response_col)]
      test_data  <- sub_data[folds$indexOut[[fold_name]], c(vars, response_col)]

      fold_fit <- gbm::gbm(
        formula = fold_formula, data = train_data, distribution = "gaussian",
        n.trees = gbm_params$n.trees, interaction.depth = gbm_params$interaction.depth,
        shrinkage = gbm_params$shrinkage, n.minobsinnode = gbm_params$n.minobsinnode
      )

      pred <- gbm::predict.gbm(fold_fit, newdata = test_data, n.trees = gbm_params$n.trees)
      tibble::tibble(obs = test_data[[response_col]], pred = pred)
    })

  sqrt(mean((fold_preds$obs - fold_preds$pred)^2))
}

## ----------------------------------------------------------------------
## Forward selection with a mandatory inclusion of pairwise 
## geographic distance ("forced vars").
## Approach based off CAST::ffs(globalval = TRUE, earlyStopping = TRUE): 
## Start from geographic distance-only model, then iteratively add whichever 
## remaining candidate most improves pooled RMSE stop if none do. 
## ----------------------------------------------------------------------
select_variables_forced <- function(sub_data, folds, forced_vars, candidate_vars, response_col,
                                     gbm_params = list(n.trees = 1000, interaction.depth = 1,
                                                       shrinkage = 0.01, n.minobsinnode = 10),
                                     seed = 1234, verbose = TRUE, label = NULL,
                                     exclusive_groups = list()) {

  ## label scenario types (modest, robust, standard)
  msg_prefix <- if (!is.null(label)) paste0("[", label, "] ") else ""

  ## exclusive_groups: sets of variables that are different representations of 
  ## the same underlying thing (e.g. several proxies for human disturbance).
  ## Once ANY member of a group is added to current_vars, every OTHER
  ## member of that group is permanently dropped from `remaining` -- so at
  ## most one variable per group can ever end up in the final selected set. 
  drop_groupmates <- function(added_var, remaining) {
    hit <- exclusive_groups %>% purrr::keep(~ added_var %in% .x)
    if (length(hit) == 0) return(remaining)
    excluded <- intersect(setdiff(unlist(hit), added_var), remaining)
    if (length(excluded) > 0 && verbose) {
      message(msg_prefix, "  ", added_var, " is in a mutually-exclusive group -- removing ",
              paste(excluded, collapse = ", "), " from further consideration")
    }
    setdiff(remaining, excluded)
  }

  current_vars <- forced_vars
  remaining    <- candidate_vars

  if (verbose) message(msg_prefix, "Base set: ", paste(current_vars, collapse = ", "))
  current_rmse <- evaluate_variable_set(sub_data, folds, current_vars, response_col, gbm_params, seed)
  if (verbose) message(msg_prefix, "  base RMSE = ", round(current_rmse, 5))

  history <- tibble::tibble(step = 0, added = NA_character_,
                             vars = paste(current_vars, collapse = ", "), rmse = current_rmse)

  step <- 1
  repeat {
    if (length(remaining) == 0) break

    if (verbose) message(msg_prefix, "Step ", step, ": testing ", length(remaining),
                          " candidate(s) -- ", paste(remaining, collapse = ", "))

    trial_rmse <- remaining %>%
      purrr::map_dbl(function(v) {
        rmse_v <- evaluate_variable_set(sub_data, folds, c(current_vars, v), response_col, gbm_params, seed)
        if (verbose) message(msg_prefix, "  testing ", v, " -> RMSE = ", round(rmse_v, 5))
        rmse_v
      })

    best_i <- which.min(trial_rmse)

    if (trial_rmse[best_i] >= current_rmse) {
      if (verbose) message(msg_prefix, "  no remaining candidate improves RMSE -- stopping")
      break   # early stopping: no remaining candidate improves RMSE
    }

    added        <- remaining[best_i]
    current_vars <- c(current_vars, added)
    current_rmse <- trial_rmse[best_i]
    if (verbose) message(msg_prefix, "  added ", added, " (RMSE = ", round(current_rmse, 5), ")")
    history <- history %>%
      dplyr::bind_rows(tibble::tibble(step = step, added = added,
                                       vars = paste(current_vars, collapse = ", "), rmse = current_rmse))

    remaining <- drop_groupmates(added, remaining[-best_i])
    step <- step + 1
  }

  if (verbose) message(msg_prefix, "Selected: ", paste(current_vars, collapse = ", "))

  list(selected_vars = current_vars, history = history)
}

## ----------------------------------------------------------------------
## Robustness check: run select_variables_forced() under several
## different FIXED hyperparameter settings and compare the resulting
## variable sets. Quick way to see whether deliberately modest
## hyperparameters (n.trees=1000, shrinkage=0.01) are excluding a
## variable that a more-boosted setting would keep.
## set robust = TRUE in 02_select_variables.R to run this 
## ----------------------------------------------------------------------
default_gbm_param_settings <- list(
  modest      = list(n.trees = 1000, interaction.depth = 1, shrinkage = 0.01,  n.minobsinnode = 10),
  generous    = list(n.trees = 3000, interaction.depth = 1, shrinkage = 0.005, n.minobsinnode = 10),
  regularized = list(n.trees = 2000, interaction.depth = 1, shrinkage = 0.005, n.minobsinnode = 25)
)

select_variables_robustness_check <- function(sub_data, folds, forced_vars, candidate_vars, response_col,
                                               gbm_param_settings = default_gbm_param_settings,
                                               verbose = TRUE, exclusive_groups = list()) {

  results_by_setting <- names(gbm_param_settings) %>%
    purrr::set_names() %>%
    purrr::map(function(setting_name) {
      select_variables_forced(sub_data, folds, forced_vars, candidate_vars, response_col,
                               gbm_params = gbm_param_settings[[setting_name]],
                               verbose = verbose, label = setting_name,
                               exclusive_groups = exclusive_groups)
    })

  stability <- candidate_vars %>%
    purrr::map_dfr(function(v) {
      selected_in <- names(results_by_setting) %>%
        purrr::keep(~ v %in% results_by_setting[[.x]]$selected_vars)
      tibble::tibble(variable = v,
                      n_settings_selected = length(selected_in),
                      n_settings_total    = length(results_by_setting),
                      settings_selected_in = paste(selected_in, collapse = ", "))
    }) %>%
    dplyr::arrange(dplyr::desc(n_settings_selected))

  list(results_by_setting = results_by_setting, stability = stability)
}
