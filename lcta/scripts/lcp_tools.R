
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Least cost path tools
#
# Script author: Eric Palm
#
################################################################################

## ============================================================================
## LCP functions: prediction-surface construction, LCP line computation. 
## ============================================================================

# Load packages
sapply(c("dplyr", "purrr", "terra", "sp", "raster", "gdistance", "sf"), 
       require, 
       character.only = T)

source("lcta/scripts/01_extract_covariates.R")
source("lcta/scripts/fit_tools.R")

## Full raster objects for every candidate landscape variable. 
static_rasters <- static_paths %>% 
  purrr::map(terra::rast) %>% 
  purrr::set_names(var_names)

## ----------------------------------------------------------------------
## Build the prediction raster stack, restricted to whatever
## predictor_cols were selected. The training median of 
## geographic distance is used a single-pixel raster for predictions. 
## Set FORCED_VAR to NULL below if there's no forced variable.
## ----------------------------------------------------------------------
FORCED_VAR <- "euc_geog"   

build_prediction_stack <- function(predictor_cols, forced_value, domain) {
  layers <- list()
  for (v in intersect(predictor_cols, names(static_rasters))) layers[[v]] <- static_rasters[[v]]
  if (!is.null(FORCED_VAR) && FORCED_VAR %in% predictor_cols) {
    layers[[FORCED_VAR]] <- static_rasters[[1]] %>% terra::setValues(forced_value)
  }
  terra::rast(layers) %>% `names<-`(names(layers)) %>% terra::crop(terra::vect(domain), mask = TRUE)
}

# One resistance surface per LCP iteration
predict_resistance_surface <- function(fit, predictor_cols, forced_value, domain) {
  env <- build_prediction_stack(predictor_cols, forced_value, domain)
  terra::predict(env[[predictor_cols]], fit) %>%
    terra::crop(terra::vect(domain), mask = TRUE)
}

## ----------------------------------------------------------------------
## LCP line computation via gdistance.
## ----------------------------------------------------------------------
compute_lcp_lines <- function(resistance, sub_pairs, domain, buffer_m = BUFFER_M) {

  max_val <- terra::global(resistance, "max", na.rm = TRUE)$max
  min_val <- terra::global(resistance, "min", na.rm = TRUE)$min
  r_norm <- terra::app(resistance, function(x) (x - min_val) / (max_val - min_val) + .01) %>%
    terra::mask(resistance)   

  ## transition() needs a raster::RasterLayer, not a terra::SpatRaster 
  ## This snippet makes sure no transition matrix values are below 0
  tr <- gdistance::transition(raster::raster(r_norm), function(x) 1 / mean(x), directions = 8) %>%
    gdistance::geoCorrection(type = "c", multpl = FALSE)
  min_pos <- min(tr@transitionMatrix@x[tr@transitionMatrix@x > 0])
  tr@transitionMatrix@x[tr@transitionMatrix@x < 0] <- min_pos

  foreach::foreach(
    i = seq_len(nrow(sub_pairs)), .combine = rbind, .verbose = FALSE,
    .packages = c("sp", "sf", "raster", "dplyr", "gdistance"),
    .export = "PROJECT_CRS"
  ) %dopar% {
    idx <- sub_pairs$index[i]
    x_1 <- sub_pairs$x_1[i]; y_1 <- sub_pairs$y_1[i]
    x_2 <- sub_pairs$x_2[i]; y_2 <- sub_pairs$y_2[i]

    gdistance::shortestPath(tr,
                             sp::SpatialPoints(cbind(x_1, y_1)),
                             sp::SpatialPoints(cbind(x_2, y_2)),
                             output = "SpatialLines") %>%
      sf::st_as_sf() %>%
      sf::st_set_crs(PROJECT_CRS) %>%
      sf::st_buffer(buffer_m) %>%
      dplyr::mutate(index = idx)
  }
}
