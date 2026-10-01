
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Permutation variable importance figure
#
# Script author: Eric Palm
#
################################################################################

## ============================================================================
## Permutation variable importance figure
##
## Normalization and confidence intervals: 
## each variable's delta_rmse (and CI) expressed as a percent of the model's 
## RMSE, with a standard-error-based 95% CI of the mean across folds.
## ============================================================================

# Load packages
sapply(c("dplyr", "stringr", "tibble", "tibble", "ggplot2", "purrr"), 
       require, 
       character.only = T)

source("lcta/scripts/cv_diagnostics.R")   # compute_fold_metrics(), cv_permutation_importance()

## ---- User configure ----
model_dir    <- "lcta/models/"
fig_dir      <- "lcta/outputs/figures/"
response_col <- "gd"    
n_repeats    <- 100      
perm_seed    <- 1
## ---------------------------------------------------------------

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

folds <- readRDS("lcta/results/bat_cv_folds.rds")   

## ----------------------------------------------------------------------
## Best LCP iteration (1:5) by pooled RMSE
## ----------------------------------------------------------------------
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

best     <- find_best_iteration(model_dir)
best_fit <- readRDS(best$file)
message("Best LCP iteration: ", best$iteration, " (RMSE = ", round(best$rmse, 4), ")")

predictor_cols <- setdiff(colnames(best_fit$trainingData), ".outcome")
sub_data <- best_fit$trainingData %>% dplyr::rename(!!response_col := .outcome)

perm <- cv_permutation_importance(best_fit, sub_data, predictor_cols, response_col, folds,
                                   n_repeats = n_repeats, seed = perm_seed)

n_folds_this <- length(folds$index)

importance <- perm$summary %>%
  dplyr::transmute(
    variable = variable,
    se_delta_rmse = sd_delta_rmse / sqrt(n_folds_this),
    mean_pct = 100 * mean_delta_rmse / best$rmse,
    se_pct   = 100 * se_delta_rmse   / best$rmse,
    lo_pct   = mean_pct - 1.96 * se_pct,
    hi_pct   = mean_pct + 1.96 * se_pct
  ) %>%
  dplyr::arrange(mean_pct) %>% 
  dplyr::mutate(label = 
                  dplyr::case_when(variable == "euc_geog" ~ "Geographic distance",
                                   variable == "tpi" ~ "Terrain position index",
                                   variable == "crop_binary" ~ "Crop",
                                   variable == "human_1" ~ "Human footprint",
                                   variable == "elev" ~ "Elevation",
                                   variable == "d_water" ~ "Distance to water",
                                   T ~ "Distance to karst")) 

saveRDS(importance, "lcta/results/permutation_importance.rds")

## ----------------------------------------------------------------------
## Plot -- single panel, variables ordered ascending by importance 
## "effect_move" is hardcoded here but you could create the ALE plot
## here and determine this using the ALE data
## ----------------------------------------------------------------------
plot_data <- importance %>%
  dplyr::mutate(label = forcats::fct_reorder(label, mean_pct)) %>% 
  dplyr::mutate(effect_move = dplyr::if_else(variable %in% c("euc_geog", "crop", "elev"),
                                            "Inhibited movement", "Facilitated movement"))

p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = label, x = mean_pct, color = effect_move,
                                             xmin = lo_pct, xmax = hi_pct)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  ggplot2::geom_pointrange(orientation = "y", fatten = 3) +
  ggplot2::labs(
    x = "Permutation importance\n(mean % increase in RMSE when shuffled)",
    y = NULL,
  ) +
  ggplot2::theme_minimal(base_size = 13) + 
  scale_color_manual(values = ggokabeito::palette_okabe_ito(order = c(3,6))) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                 legend.position = "bottom",
                 axis.text = ggplot2::element_text(color = "black"),
                 panel.grid.major = ggplot2::element_line(linetype = "dotted",
                                                          color = "gray20",
                                                          linewidth = .2),
                 text = ggplot2::element_text(family = "Helvetica"),
                 axis.title.x = ggplot2::element_text(size = 11)) +
  ggplot2::labs(color = NULL) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, .01)))

ggplot2::ggsave(file.path(fig_dir, "permutation_importance.tiff"), p, width = 6, 
                height = 5, dpi = 600, bg = "white", compression = "lzw")

