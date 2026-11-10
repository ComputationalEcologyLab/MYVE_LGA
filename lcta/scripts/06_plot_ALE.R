
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Create accumulated local effects (ALE) plots from top model
#
# Script author: Eric Palm
#
################################################################################

## ============================================================================
## Accumulated local effects figure
##
## These allow you to plot the relationships between explanatory variables
## and the response variable while accounting for collinear predictors. 
## 
## See the `ale` package documentation here: https://tripartio.github.io/ale/
## ============================================================================

# Might be necessary to set for cluster management depending on your machine
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

# Load packages
sapply(c("dplyr", "gbm", "caret", "ale", "stringr", "parallel", "tidyr",
         "tibble", "purrr"), 
       require, 
       character.only = T)

## Functions

# Add random variable
call_rand <- function(x) {
  stringr::str_c("gbm::gbm(y ~ ", str_c(x$trainingData %>% colnames() %>% head(-1), collapse = " + "), " + random_variable",
                 ", data = rand_data", ", n.trees = ", x$bestTune$n.trees, 
                 ", interaction.depth = ", x$bestTune$interaction.depth,
                 ", shrinkage = ", x$bestTune$shrinkage, ", n.minobsinnode = ", x$bestTune$n.minobsinnode, ")")
}

# Set up boot data for model bootstrapping
call_boot <- function(x) {
  stringr::str_c("gbm::gbm(y ~ ", str_c(x$trainingData %>% colnames() %>% head(-1), collapse = " + "), ", data = boot_data",
                 ", n.trees = ", x$bestTune$n.trees, ", interaction.depth = ", x$bestTune$interaction.depth,
                 ", shrinkage = ", x$bestTune$shrinkage, ", n.minobsinnode = ", x$bestTune$n.minobsinnode, ")")
}

# Function to generate random variable distributions
aledist_it <- function(model, call) {
  ale::ALEpDist(model = model$finalModel, data = model$trainingData %>% dplyr::rename(y = .outcome),
                y_col = "y", pred_type = "response",
                random_model_call_string = call, seed = 1234,
                rand_it = 1000, parallel = parallel::detectCores(logical = FALSE) - 25)
}

# Model bootstrapping to get confidence intervals around ALE curves
boot_it <- function(model, call, dist) {
  ale::ModelBoot(model = model$finalModel, data = model$trainingData %>% dplyr::rename(y = .outcome),
                 y_col = "y", pred_type = "response",
                 model_call_string = call,
                 ale_options = list(max_num_bins = 20),
                 model_packages = c("caret", "gbm"),
                 ale_p = dist,
                 seed = 1234,
                 boot_centre = "median",
                 output_model_stats = FALSE,
                 output_model_coefs = FALSE,
                 boot_it = 500, parallel = parallel::detectCores(logical = FALSE) - 25)
}

## ----------------------------------------------------------------------
## Best LCP iteration (1:5) by pooled RMSE
## ----------------------------------------------------------------------
source("lcta/scripts/cv_diagnostics.R")
## If running with a subset of variables (excluding karst vars and eco_niche)
## If not, comment this line out
source("lcta/scripts/run_config.R")

model_dir <- paste0("lcta/models/", run_suffix)

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
model <- readRDS(best$file)
message("Best LCP iteration: ", best$iteration, " (RMSE = ", round(best$rmse, 4), ")")

####
call_r <- call_rand(model)

call_b <- call_boot(model)

distie <- aledist_it(model, call = call_r)

ale <- boot_it(model, call = call_b, dist = distie) 

saveRDS(ale, paste0("lcta/results/ale", run_suffix, ".rds"))

# Extract the ALE data and format for plot
for_plot <- 
  ale %>%
  get(., type = "boot") %>% 
  purrr::map(., ~dplyr::rename(., !!"x" := names(.)[1])) %>%
  purrr::map_dfr(., ~.x, .id = "variable") %>% 
  dplyr::rename(y = .y_median) %>% 
  dplyr::mutate(x = dplyr::if_else(stringr::str_detect(variable, "karst|water"), x/1000, x),
                label_ale = dplyr::case_when(variable == "euc_geog" ~ "Geographic distance (km)",
                                             variable == "tri" ~ "Terrain ruggedness index",
                                             variable == "tpi" ~ "Terrain position index",
                                             variable == "crop_binary" ~ "Crop (binary)",
                                             variable == "human_1" ~ "Human footprint",
                                             variable == "elev" ~ "Elevation (m)",
                                             variable == "d_water" ~ "Distance to water (km)",
                                             T ~ "Distance to karst (km)"))

# Get the ALER band
aler <- 
  ale@ale$single@params$y_summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., "stat") %>%
  tidyr::pivot_wider(names_from = stat, values_from = y) %>% 
  dplyr::select(aler_lo, aler_hi) 

# Variable importance to order panels in plot by their relative importance
imp <- readRDS(paste0("lcta/results/permutation_importance", run_suffix, ".rds"))

# Add in rugs for each variable to show where most data are in the distribution
rugs <- 
  model$trainingData %>% 
  dplyr::select(-.outcome) %>%
  tidyr::pivot_longer(., cols = everything(), names_to = "variable", values_to = "x") %>% 
  mutate(x = if_else(str_detect(variable, "karst|water"), x/1000, x),
         label_ale = dplyr::case_when(variable == "euc_geog" ~ "Geographic distance (km)",
                                      variable == "tri" ~ "Terrain ruggedness index",
                                      variable == "tpi" ~ "Terrain position index",
                                       variable == "crop_binary" ~ "Crop (binary)",
                                      variable == "human_1" ~ "Human footprint",
                                      variable == "elev" ~ "Elevation (m)",
                                      variable == "d_water" ~ "Distance to water (km)",
                                      variable == "d_karst" ~ "Distance to karst (km)"),
         y = min(for_plot$y))

# Restrict x-axis to the middle 99.8% of data
var_ranges_rug <- rugs %>% 
  dplyr::summarize(min = quantile(x, .001),
                   max = quantile(x, .999), .by = variable)

# Format the ALER band data
for_plot_aler <- for_plot %>% 
  dplyr::mutate(aler_lo = aler$aler_lo,
                aler_hi = aler$aler_hi) %>% 
  dplyr::inner_join(., var_ranges_rug) %>% 
  dplyr::filter(x > min,
                x < max) %>%
  dplyr::inner_join(., imp) %>% 
  dplyr::mutate(label_ale = forcats::fct_reorder(label_ale, mean_pct, .desc = T))

# Restrict x-axis to the middle 99.8% of data
var_ranges_ale <- for_plot_aler %>% 
  dplyr::summarize(min = quantile(x, .001),
                   max = quantile(x, .999), .by = variable)

# Prepare rug data for plotting
rugs_filtered <- rugs %>% 
  dplyr::inner_join(., var_ranges_ale) %>% 
  dplyr::inner_join(., imp) %>% 
  dplyr::filter(x > min,
                x < max) %>%
  dplyr::mutate(label_ale = forcats::fct_reorder(label_ale, mean_pct, .desc = T))

minmax <- c(min(for_plot_aler$.y_lo), c(max(for_plot_aler$.y_hi)))

for_plot_aler %>% 
  ggplot(.) +
  geom_ribbon(aes(x = x, ymin = .y_lo, ymax = .y_hi), alpha = .5) +
  geom_rect(data = aler, aes(xmin = -Inf, xmax = Inf, ymin = aler_lo, ymax = aler_hi),
            alpha = .25) +
  geom_hline(data = aler, aes(yintercept = aler_hi), linetype = "dashed", linewidth = .4) +
  geom_hline(data = aler, aes(yintercept = aler_lo), linetype = "dashed", linewidth = .4) +
  geom_rug(data = rugs_filtered, sides = "b", alpha = .05, length = unit(0.03, "npc"), aes(x = x, y = y)) +
  geom_line(aes(x = x, y = y), linewidth = .8) +
  theme_classic(base_size = 16) +
  facet_wrap(~label_ale, scales = "free_x", strip.position = "bottom") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text = element_text(color = "black"),
        axis.line = element_line(linewidth = .5),
        legend.position = "bottom",
        text = element_text(family = "Helvetica"),
        strip.text = element_text(vjust = 1.5, size = 14),
        plot.margin = margin(r = 15, l = 5, t = 5),
        panel.spacing.y = unit(-.3, "cm"),
        legend.key.width=unit(1.1, "cm")) +
  labs(x = NULL, y = "Predicted genetic distance") +
  scale_y_continuous(limits = minmax, expand = expansion(mult = c(.02, .01))) +
  scale_x_continuous(expand = expansion(mult = c(.005, 0))) 

ggsave(paste0("lcta/outputs/figures/ale_boot", run_suffix, ".tiff"), height = 9, width = 11, compression = "lzw")

