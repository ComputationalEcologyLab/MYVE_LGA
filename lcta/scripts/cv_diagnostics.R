
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Spatial cross validation diagnostics
#
# Script author: Eric Palm
#
################################################################################

## ============================================================================
## Post-hoc validation and diagnostic tools -- run AFTER fitting a model
## ============================================================================

# Load packages
sapply(c("dplyr", "parallel", "foreach", "doParallel"), 
       require, 
       character.only = T)

## ============================================================================
## Fold-aware validation -- run AFTER fitting, against fit$pred (requires
## trainControl(savePredictions = "final")).
##
## RMSE is primary metric
## ============================================================================

rmse      <- function(obs, pred) sqrt(mean((obs - pred) ^ 2))
r_squared <- function(obs, pred) 1 - sum((obs - pred) ^ 2) / sum((obs - mean(obs)) ^ 2)

## Vectorized -- statistically identical to resampling n_folds folds with
## replacement and recomputing RMSE/R2 on the pooled rows each time

compute_fold_metrics <- function(pred_df, n_boot = 2000) {

  n_na <- sum(is.na(pred_df$obs) | is.na(pred_df$pred))
  if (n_na > 0) {
    pred_df <- pred_df %>% dplyr::filter(!is.na(obs), !is.na(pred))
    message("compute_fold_metrics(): dropping ", n_na, " row(s) with NA obs/pred",
            " (likely a fold that failed entirely under the chosen hyperparameters --",
            " see fit_gbm_cv()'s n.minobsinnode feasibility filter in fit_helpers.R)",
            " before computing metrics.")
  }

  fold_id <- pred_df$Resample %>% stringr::str_remove("^Fold") %>% as.integer()

  fold_stats <- tibble::tibble(fold = fold_id, obs = pred_df$obs, pred = pred_df$pred) %>%
    dplyr::group_by(fold) %>%
    dplyr::summarise(n = dplyr::n(), sse = sum((obs - pred)^2),
                      sum_obs = sum(obs), sum_obs2 = sum(obs^2), .groups = "drop") %>%
    dplyr::arrange(fold)

  ## Exclude folds with no rows so it isn't sampled as a real
  ## fold in the bootstrap below.
  n_empty_folds <- sum(fold_stats$n == 0)
  if (n_empty_folds > 0) {
    message("compute_fold_metrics(): ", n_empty_folds, " fold(s) had zero rows after",
            " dropping NAs -- excluding from the pooled estimate and the bootstrap CI.")
    fold_stats <- fold_stats %>% dplyr::filter(n > 0)
  }

  n_folds <- nrow(fold_stats)
  stopifnot("compute_fold_metrics(): every fold was empty after dropping NAs -- nothing to compute" =
              n_folds > 0)

  global_metrics <- tibble::tibble(
    rmse = sqrt(sum(fold_stats$sse) / sum(fold_stats$n)),
    r2   = 1 - sum(fold_stats$sse) / (sum(fold_stats$sum_obs2) - sum(fold_stats$sum_obs)^2 / sum(fold_stats$n))
  )

  ## n_folds x n_boot matrix: how many times each fold is drawn in each
  ## bootstrap replicate (equivalent to sampling n_folds folds with
  ## replacement, n_boot times)
  counts <- stats::rmultinom(n_boot, size = n_folds, prob = rep(1 / n_folds, n_folds))

  boot_n        <- as.numeric(t(counts) %*% fold_stats$n)
  boot_sse      <- as.numeric(t(counts) %*% fold_stats$sse)
  boot_sum_obs  <- as.numeric(t(counts) %*% fold_stats$sum_obs)
  boot_sum_obs2 <- as.numeric(t(counts) %*% fold_stats$sum_obs2)

  boot_metrics <- tibble::tibble(
    rmse = sqrt(boot_sse / boot_n),
    r2   = 1 - boot_sse / (boot_sum_obs2 - boot_sum_obs^2 / boot_n)
  )

  ci <- boot_metrics %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "metric", values_to = "value") %>%
    dplyr::group_by(metric) %>%
    dplyr::summarise(lo = stats::quantile(value, .025), median = stats::quantile(value, .5),
                      hi = stats::quantile(value, .975), .groups = "drop")

  list(global = global_metrics, ci = ci, boot_draws = boot_metrics)
}

## Use the pooled RMSE as in CAST::global_validation(fit) because
## knndm folds are commonly unbalanced in size.

##
## ============================================================================
## CV permutation importance -- fits per fold on TRAINING rows,
## but measures importance on each fold's HELD-OUT test rows: for each
## predictor, shuffle it in the test data, re-predict, and record how much
## test RMSE degrades relative to the unshuffled baseline. This is a held-out, 
## spatially-blocked measure of how much each variable contributes to explaining
## the response in data the model never trained on.
##
## delta_rmse is the importance scale: how much worse (in RMSE units) test
## predictions get when that variable's information is destroyed. Repeats
## the shuffle n_repeats times per variable per fold and averages, since a
## single permutation draw is noisy. 
## ============================================================================

cv_permutation_importance <- function(final_model, sub_data, predictor_cols, response_col, folds,
                                       n_repeats = 10, seed = 1) {

  set.seed(seed)

  ## accept either a caret train object (fit) or the raw gbm object
  ## (fit$finalModel) -- only the latter has $distribution/$n.trees/etc.
  ## at the top level, so unwrap automatically if a train object was passed
  ## (same as cv_variable_importance() above)
  if (inherits(final_model, "train")) final_model <- final_model$finalModel

  ## hyperparameters sourced from the saved final model; the model object
  ## itself is never used for prediction here -- each fold needs a fresh
  ## model fit on only that fold's training rows (final_model was fit on
  ## ALL of sub_data, so evaluating it directly against any fold's
  ## "held-out" rows wouldn't be held-out at all)
  n_trees           <- final_model$n.trees
  interaction_depth <- final_model$interaction.depth
  shrinkage         <- final_model$shrinkage
  n_minobsinnode    <- final_model$n.minobsinnode
  distribution      <- final_model$distribution$name

  stopifnot(
    "folds must be the fold list from 00_build_cv_folds.R (with named $index/$indexOut)" =
      !is.null(folds$index) && length(folds$index) > 0 && !is.null(names(folds$index)),
    "response_col not found in sub_data" =
      response_col %in% colnames(sub_data),
    "some predictor_cols not found in sub_data" =
      all(predictor_cols %in% colnames(sub_data))
  )

  fold_formula <- stats::reformulate(sprintf("`%s`", predictor_cols),
                                      response = sprintf("`%s`", response_col))

  imp_by_fold <- names(folds$index) %>%
    purrr::map_dfr(function(fold_name) {

      train_data <- sub_data[folds$index[[fold_name]], c(predictor_cols, response_col)]
      test_data  <- sub_data[folds$indexOut[[fold_name]], c(predictor_cols, response_col)]

      fold_fit <- gbm::gbm(
        formula = fold_formula, data = train_data, distribution = distribution,
        n.trees = n_trees, interaction.depth = interaction_depth,
        shrinkage = shrinkage, n.minobsinnode = n_minobsinnode
      )

      baseline_pred <- gbm::predict.gbm(fold_fit, newdata = test_data, n.trees = n_trees)
      baseline_rmse <- sqrt(mean((test_data[[response_col]] - baseline_pred) ^ 2))

      predictor_cols %>%
        purrr::map_dfr(function(var) {
          rmse_shuffled <- seq_len(n_repeats) %>%
            purrr::map_dbl(function(i) {
              shuffled <- test_data
              shuffled[[var]] <- sample(shuffled[[var]])
              shuffled_pred <- gbm::predict.gbm(fold_fit, newdata = shuffled, n.trees = n_trees)
              sqrt(mean((shuffled[[response_col]] - shuffled_pred) ^ 2))
            })
          tibble::tibble(variable = var, fold = fold_name,
                          baseline_rmse = baseline_rmse,
                          delta_rmse = mean(rmse_shuffled) - baseline_rmse)
        })
    })

  summary_tbl <- imp_by_fold %>%
    dplyr::summarise(
      mean_delta_rmse = mean(delta_rmse),
      sd_delta_rmse   = stats::sd(delta_rmse),
      lo = stats::quantile(delta_rmse, .025),
      hi = stats::quantile(delta_rmse, .975),
      .by = variable
    ) %>%
    dplyr::arrange(dplyr::desc(mean_delta_rmse))

  list(summary = summary_tbl, by_fold = imp_by_fold)
}
