
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Model fitting tools for gradient boosting machines
#
# Script author: Eric Palm
#
################################################################################

## ============================================================================
## Tools for GBMs used in 03_fit_straight_line.R and 04_lcp_iteration_cycle.R
## ============================================================================

# Load packages
sapply(c("caret", "dplyr", "gbm"), 
         require, 
         character.only = T)

# Create default hyperparameter grid 
default_gbm_grid <- expand.grid(
  interaction.depth = c(1, 2),
  n.trees            = seq(100, 3000, 10),
  shrinkage          = 0.005,
  n.minobsinnode     = c(5, 10, 15, 25, 50)
)

## Sub_data must be row-aligned with folds$index/indexOut
fit_gbm_cv <- function(sub_data, predictor_cols, response_col, folds,
                        t_grid = default_gbm_grid, seed = 1234, bag_fraction = 0.5) {

  ## kNNDM folds are unbalanced, so the smallest fold can be much
  ## smaller than the full training set. This drops tuneGrid rows that will be
  ## infeasible for this dataset's SMALLEST fold. gbm's
  ## default bag.fraction=0.5 means each tree is fit on a random half of a
  ## fold's training rows, and a tree of interaction.depth d needs roughly
  ## 2^d leaves, each needing >= n.minobsinnode observations. 
  ## Without this filter, gbm model will fail. 
  min_fold_train_n <- min(lengths(folds$index))
  feasible <- t_grid$n.minobsinnode <= floor(bag_fraction * min_fold_train_n / 2^t_grid$interaction.depth)
  if (!all(feasible)) {
    message("fit_gbm_cv(): dropping ", sum(!feasible), " of ", nrow(t_grid),
            " tuneGrid row(s) infeasible for this dataset's smallest fold (",
            min_fold_train_n, " training rows) -- see the comment above fit_gbm_cv().")
    t_grid <- t_grid[feasible, , drop = FALSE]
  }

  n_folds <- folds$n_folds
  seeds <- vector(mode = "list", length = n_folds + 1)
  set.seed(seed)
  for (i in seq_len(n_folds)) seeds[[i]] <- rep(seed, nrow(t_grid))
  seeds[[n_folds + 1]] <- seed

  ctrl <- caret::trainControl(
    method          = "cv",
    index           = folds$index,
    indexOut        = folds$indexOut,
    allowParallel   = TRUE,
    savePredictions = "final",
    returnResamp    = "none",
    seeds           = seeds
  )

  caret::train(
    x = sub_data[, predictor_cols],
    y = sub_data[[response_col]],
    method    = "gbm",
    trControl = ctrl,
    metric    = "RMSE",
    tuneGrid  = t_grid,
    verbose   = FALSE
  )
}
