
################################################################################
# Analysis: Myotis velifer landscape genetics - least cost transect analysis
#
# Script title: 00_build_cv_folds
#
# Script author: Eric Palm
#
################################################################################


## ============================================================================
## Transect-based kNNDM spatial CV folds for fitting a GBM using transects.
##
## Builds one set of cross-validation folds for the whole dataset. Output
## is a single list formatted for caret::trainControl(index=, indexOut=).
##
## ============================================================================

# Load packages
sapply(c("dplyr", "tidyr", "purrr", "tibble", "stringr", "caret", "sf",
         "units", "twosamples"),
       require, 
       character.only = T)

source("lcta/scripts/01_extract_covariates.R")   # for PROJECT_CRS

set.seed(1)

## Set constants
MCP_BUFFER_M <- 10000   # buffer (m) around the sample-location MCP --
                         # defines the model domain, used both for
                         # predpoints here and reused as the raster-
                         # cropping domain in 04_lcp_iteration_cycle.R
N_PREDPOINTS <- 5000    # number of points sampled from the model domain
                         # for kNNDM's prediction-distance distribution
MIN_DIST     <- 0       # optional euc_geog filter on pairs 
## ---------------------------------------------------------------

## ----------------------------------------------------------------------
## Minimum convex polygon around a set of sf points -- used below
## as the model domain for predpoints. percentile = 100 keeps all points
## ----------------------------------------------------------------------
mcp <- function(x, percentile = 100) {
  centroid <- sf::st_centroid(sf::st_union(x))
  dist <- as.numeric(sf::st_distance(x, centroid))
  within_range <- dist <= stats::quantile(dist, percentile / 100)
  sf::st_convex_hull(sf::st_union(x[within_range, ]))
}

## ----------------------------------------------------------------------
## Transect-based kNNDM: uses point-to-nearest-TRANSECT
## distance rather than point-to-point distance, so sampling unit is a 
## pair of locations rather than a single location.
## ----------------------------------------------------------------------
knndm_transect <- function(tpoints, pt2transect, pred2transect, transect_ends,
                            k = 10, maxp = 0.5, linkf = "ward.D2", nk_len = 100) {

  tcoords <- sf::st_coordinates(tpoints)
  n_animals <- nrow(tcoords)
  gij <- apply(pred2transect, 1, min)

  cv_distances_transect <- function(clust_k) {
    vapply(seq_len(n_animals), function(i) {
      f <- clust_k[i]
      train_cols <- !(clust_k[transect_ends[, 1]] == f | clust_k[transect_ends[, 2]] == f)
      if (!any(train_cols)) return(NA_real_)
      min(pt2transect[i, train_cols])
    }, numeric(1))
  }

  hc <- stats::hclust(stats::dist(tcoords), method = linkf)

  clustgrid <- data.frame(
    nk = unique(as.integer(round(exp(seq(log(k), log(n_animals - 2), length.out = nk_len)))))
  )
  clustgrid$w_stat <- NA_real_
  clustgroups <- list()
  pcacoords <- stats::prcomp(tcoords, center = TRUE, scale. = FALSE, rank. = 1)

  for (nk in clustgrid$nk) {

    clust_nk <- stats::cutree(hc, k = nk)
    tabclust <- as.data.frame(table(clust_nk))
    tabclust$clust_k <- NA_integer_

    centr <- vapply(tabclust$clust_nk, function(x) {
      m <- matrix(colMeans(tcoords[clust_nk %in% x, , drop = FALSE]), nrow = 1)
      colnames(m) <- colnames(tcoords)
      stats::predict(pcacoords, m)
    }, numeric(1))
    tabclust$centrpca <- centr
    tabclust <- tabclust[order(tabclust$centrpca), ]

    clust_i <- 1
    for (i in seq_len(nrow(tabclust))) {
      if (tabclust$Freq[i] >= n_animals / k) {
        tabclust$clust_k[i] <- clust_i
        clust_i <- clust_i + 1
      }
    }
    remaining <- setdiff(seq_len(k), unique(tabclust$clust_k))
    if (length(remaining) > 0) {
      na_idx <- is.na(tabclust$clust_k)
      tabclust$clust_k[na_idx] <- rep(remaining, ceiling(nk / length(remaining)))[seq_len(sum(na_idx))]
    }

    id_map <- merge(data.frame(pos = seq_along(clust_nk), clust_nk = clust_nk), tabclust, by = "clust_nk")
    id_map <- id_map[order(id_map$pos), ]
    clust_k <- id_map$clust_k

    if (!any(table(clust_k) / length(clust_k) > maxp)) {
      gjstar_i <- cv_distances_transect(clust_k)
      clustgrid$w_stat[clustgrid$nk == nk] <- twosamples::wass_stat(stats::na.omit(gjstar_i), gij)
      clustgroups[[paste0("nk", nk)]] <- clust_k
    }
  }

  k_final <- clustgrid$nk[which.min(clustgrid$w_stat)]
  clust <- clustgroups[[paste0("nk", k_final)]]

  list(clusters = stats::setNames(clust, tpoints$id),
       w_stat = min(clustgrid$w_stat, na.rm = TRUE), q = k_final)
}

## ----------------------------------------------------------------------
## Try several candidate number of folds, k, and keep whichever gives the
## lowest W statistic as per CAST::knndm(). 
## ----------------------------------------------------------------------
choose_k <- function(tpoints, pt2transect, pred2transect, transect_ends,
                      k_candidates, maxp = 0.5, linkf = "ward.D2", nk_len = 100) {

  attempts <- k_candidates %>%
    purrr::map(function(k) {
      res <- tryCatch(
        knndm_transect(tpoints, pt2transect, pred2transect, transect_ends,
                        k = k, maxp = maxp, linkf = linkf, nk_len = nk_len),
        error = function(e) NULL
      )
      if (is.null(res)) return(NULL)
      tibble::tibble(k = k, w_stat = res$w_stat, q = res$q)
    }) %>%
    purrr::compact() %>%
    dplyr::bind_rows()

  list(table = attempts, best_k = attempts$k[which.min(attempts$w_stat)])
}

## ----------------------------------------------------------------------
## Build folds for the whole dataset. min_dist/max_dist optionally filter
## pairs by euc_geog before folds are built.
## ----------------------------------------------------------------------
build_folds <- function(pairs, min_dist = MIN_DIST, max_dist = MAX_DIST, crs = PROJECT_CRS,
                         mcp_buffer_m = MCP_BUFFER_M, n_predpoints = N_PREDPOINTS,
                         mcp_percentile = 100) {

  message("Building folds")

  orig_index <- pairs %>%
    dplyr::filter(.data$euc_geog > .env$min_dist, .data$euc_geog <= .env$max_dist) %>%
    dplyr::pull(index)

  if (length(orig_index) == 0) stop("No eligible pairs")

  pairs_sub <- pairs %>% dplyr::filter(index %in% orig_index)

  sample_locs <- dplyr::bind_rows(
    dplyr::transmute(pairs_sub, id = id_1, x = x_1, y = y_1),
    dplyr::transmute(pairs_sub, id = id_2, x = x_2, y = y_2)
  ) %>%
    dplyr::distinct()

  tpoints <- sample_locs %>%
    sf::st_as_sf(coords = c("x", "y"), crs = crs, remove = FALSE)

  ## Model domain: MCP around every actual sample location, 
  ## buffered by mcp_buffer_m. This is reused as the raster-cropping domain in 
  ## 04_lcp_iteration_cycle.R for consistency.
  model_domain <- sample_locs %>%
    sf::st_as_sf(coords = c("x", "y"), crs = crs) %>%
    mcp(percentile = mcp_percentile) %>%
    sf::st_buffer(mcp_buffer_m) %>%
    sf::st_as_sf()

  predpoints <- sf::st_sample(model_domain, size = n_predpoints, type = "regular")

  ## straight-line transect geometry from the endpoint coordinates
  transect_lines <- purrr::map(seq_len(nrow(pairs_sub)), function(i) {
    sf::st_linestring(rbind(c(pairs_sub$x_1[i], pairs_sub$y_1[i]),
                             c(pairs_sub$x_2[i], pairs_sub$y_2[i])))
  }) %>%
    sf::st_sfc(crs = crs)

  transect_ends <- cbind(match(pairs_sub$id_1, tpoints$id), match(pairs_sub$id_2, tpoints$id))

  pt2transect   <- sf::st_distance(tpoints, transect_lines)    %>% units::drop_units()
  pred2transect <- sf::st_distance(predpoints, transect_lines) %>% units::drop_units()

  n_animals <- nrow(tpoints)
  k_candidates <- c(4, 5, 6, 8, 10) %>% purrr::keep(~ .x <= max(2, floor(n_animals / 3)))
  if (length(k_candidates) == 0) k_candidates <- max(2, floor(n_animals / 3))

  k_search <- choose_k(tpoints, pt2transect, pred2transect, transect_ends, k_candidates)
  message("  k candidates tried: ", paste(k_search$table$k, collapse = ", "),
          " -- chose k = ", k_search$best_k, " (lowest W)")

  res <- knndm_transect(tpoints, pt2transect, pred2transect, transect_ends, k = k_search$best_k)
  clust <- res$clusters
  n_folds <- length(unique(clust))

  c1 <- clust[match(pairs_sub$id_1, names(clust))]
  c2 <- clust[match(pairs_sub$id_2, names(clust))]
  coin <- sample(c(1L, 2L), nrow(pairs_sub), replace = TRUE)
  test_cluster <- ifelse(c1 == c2, c1, ifelse(coin == 1L, c1, c2))

  fold_indices <- purrr::map(seq_len(n_folds), function(f) {
    touches_f <- (c1 == f) | (c2 == f)
    list(train = which(!touches_f), test = which(test_cluster == f))
  })

  list(
    index        = purrr::map(fold_indices, "train") %>% purrr::set_names(paste0("Fold", seq_len(n_folds))),
    indexOut     = purrr::map(fold_indices, "test")  %>% purrr::set_names(paste0("Fold", seq_len(n_folds))),
    orig_index   = orig_index,
    clusters     = clust,
    n_folds      = n_folds,
    w_stat       = res$w_stat,
    q            = res$q,
    k_search     = k_search$table,
    model_domain = model_domain,
    n_pairs      = nrow(pairs_sub),
    n_animals    = n_animals
  )
}

## ----------------------------------------------------------------------
## Build and save folds
## ----------------------------------------------------------------------

## Input data frame needs: index, id_1, id_2, x_1, y_1, x_2, y_2, euc_geog.
pairs <- readRDS("lcta/data/samples_pairwise.rds") %>%
  dplyr::mutate(index = dplyr::row_number())
## ----------------------------------------------------------------

folds <- build_folds(pairs)

dir.create("lcta/results", recursive = TRUE, showWarnings = FALSE)
saveRDS(folds, "lcta/results/bat_cv_folds.rds")
