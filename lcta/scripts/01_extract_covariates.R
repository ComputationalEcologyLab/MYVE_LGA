
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: Covariate extraction along straight lines or least cost paths
#
# Script author: Eric Palm
#
################################################################################

## ============================================================================
## Shared covariate extraction -- used identically for iteration 0
## (straight-line transects) and every LCP iteration. 
##
## Lines are buffered by BUFFER_M (default 500m) before extraction.
##
## ============================================================================

# Load packages
sapply(c("foreach", "parallel", "doParallel", "sf", "exactextractr", "terra",
         "dplyr"), 
       require, 
       character.only = T)

## ---- User configure ----
BUFFER_M <- 1000   # buffer applied to straight-line transects/LCP corridors
                   # before extracting covariates along them

covs_dir  <- "data/SpatialLayers/"

# Manually set names in the same order as the .tif files in the directory above
var_names <- c("crop_binary", "d_crop", "d_karst", "d_water", "elev", "eco_niche", 
               "forest_prop_250", "forest_prop_gcam", "human_1", "human_2",
               "karst_prop", "lights", "slope", "tpi", "tri", "water_prop")


static_paths <- list.files(covs_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(
  "number of .tif files in covs_dir doesn't match length(var_names) -- update whichever is out of date" =
    length(static_paths) == length(var_names)
)

## Project CRS -- read directly off an actual raster 

PROJECT_CRS <- terra::crs(terra::rast(static_paths[1]))
## ---------------------------------------------------------------

## ----------------------------------------------------------------------
## Build straight-line LINESTRINGs (iteration 0 only) from gen_table's
## y_1/x_1/y_2/x_2. LCP iterations skip this entirely -- they
## already have line geometry from the LCP-computation step, and feed
## straight into extract_covariates() below.
## ----------------------------------------------------------------------
lines_from_pairs <- function(gen_table, crs = PROJECT_CRS) {
  stopifnot(all(c("index", "y_1", "x_1", "y_2", "x_2") %in% names(gen_table)))

  purrr::map(seq_len(nrow(gen_table)), function(i) {
    sf::st_linestring(rbind(
      c(gen_table$x_1[i], gen_table$y_1[i]),
      c(gen_table$x_2[i], gen_table$y_2[i])
    ))
  }) %>%
    sf::st_sfc(crs = crs) %>%
    sf::st_sf(index = gen_table$index, geometry = .)
}

## ----------------------------------------------------------------------
## Extract every candidate covariate along a set of lines.
##
## `lines` : sf object with an `index` column matching gen_table. Either
## raw LINESTRINGs to be buffered here (straight-line, iteration 0 -- the
## default), OR an already-buffered corridor polygon (LCP iterations,
## whose lines come out of compute_lcp_lines() in 04_lcp_iteration_cycle.R
## already at BUFFER_M) -- use already_buffered = TRUE in that case, or
## this would buffer a second time on top of an already-buffered polygon.
## ----------------------------------------------------------------------
extract_covariates <- function(lines, gen_table, buffer_m = BUFFER_M, already_buffered = FALSE) {
  stopifnot("index" %in% names(lines))

  ## var_names/static_paths come from those defined manually above
  n_layers <- length(static_paths)

  buffered <- if (already_buffered) lines else sf::st_buffer(lines, buffer_m)

  ## Cluster setup
  out <- foreach::foreach(
    i = seq_len(nrow(buffered)), .combine = rbind,
    .packages = c("sf", "exactextractr", "terra"),
    .export = "static_paths"
  ) %dopar% {
    env_local <- terra::rast(static_paths)
    vals <- exactextractr::exact_extract(env_local, buffered[i, ])[[1]]
    vals <- vals[, seq_len(n_layers), drop = FALSE]
    c(index = buffered$index[i], colMeans(vals, na.rm = TRUE))
  }

  covariate_tbl <- as.data.frame(out) %>%
    `names<-`(c("index", var_names)) %>%
    tibble::as_tibble()

  ## join back onto gen_table for the response/grouping columns
  covariate_tbl %>% dplyr::inner_join(gen_table, by = "index")
}
