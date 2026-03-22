# 00_functions.R
# Utility functions for landscape simulation and geometric analysis.
# Extracted from the geomUtil R package by Kyle W. Taylor <kwtaylor@ualberta.ca>
# Dependencies: terra, dplyr, tidyr, landscapeR

# Landscape generation ----------------------------------------------------

#' Create a blank raster landscape filled with forest cells
create_landscape <- function(dimension      = 316,
                             forest_code    = 10,
                             harvest_code   = 20,
                             retention_code = 30,
                             edge_code      = 40) {
  r <- terra::rast(matrix(forest_code, dimension, dimension))
  r <- terra::as.factor(r)
  lv <- data.frame(
    value = c(forest_code, harvest_code, retention_code, edge_code),
    class = c("Forest", "Harvest", "Retention", "Edge")
  )
  terra::set.cats(r, layer = 1, value = lv)
  r
}

#' Return cell indices in the interior of a raster, excluding a border buffer
get_interior_cells <- function(r, patch_size) {
  n_margin_cells <- ceiling(sqrt(patch_size))
  row_rng <- (1 + n_margin_cells):(nrow(r) - n_margin_cells)
  col_rng <- (1 + n_margin_cells):(ncol(r) - n_margin_cells)
  terra::cellFromRowColCombine(r, row_rng, col_rng)
}

#' Generate random seed cell indices within a raster
generate_seeds <- function(r,
                           n_seeds        = 10,
                           protect_border = FALSE,
                           patch_size     = 25 * 100,
                           rng_seed       = NULL) {
  stopifnot(inherits(r, "SpatRaster"))
  if (!is.null(rng_seed)) set.seed(as.integer(rng_seed))

  border_width <- round(sqrt(patch_size))

  if (protect_border && border_width > 0) {
    cells <- get_interior_cells(r, patch_size = border_width^2)
  } else {
    cells <- seq_len(terra::ncell(r))
  }

  stopifnot(n_seeds <= length(cells))
  selected <- sample(cells, n_seeds)
  selected <- selected[sample(seq_along(selected))]

  list(raster = r, seeds = selected)
}

#' Grow multiple harvest patches in parallel without overlap
add_patches_parallel <- function(x,
                                 seeds        = NULL,
                                 patch_size,
                                 forest_code  = 10,
                                 harvest_code = 20) {
  if (inherits(x, "list")) {
    r     <- x$raster
    seeds <- x$seeds
  } else {
    r <- x
  }
  stopifnot(inherits(r, "SpatRaster"))
  stopifnot(is.numeric(seeds), length(seeds) > 0)

  rasters <- lapply(seeds, function(s) {
    landscapeR::makeClass(
      context = r,
      size    = patch_size,
      npatch  = 1,
      pts     = s,
      bgr     = forest_code,
      val     = harvest_code,
      rast    = TRUE
    )
  })

  merged <- terra::app(terra::rast(rasters), fun = max)
  attr(merged, "patch_seeds") <- seeds
  merged
}

#' Add harvest patches and retention islands to a landscape
add_retention_harvest <- function(landscape_matrix,
                                  n_harvests        = 10,
                                  harvest_size      = 25 * 100,
                                  percent_retention = 0.05,
                                  forest_code       = 10,
                                  harvest_code      = 20,
                                  island_code       = 30,
                                  protect_border    = TRUE) {
  if (inherits(landscape_matrix, "list")) {
    r_in  <- landscape_matrix$raster
    seeds <- landscape_matrix$seeds
  } else {
    r_in  <- landscape_matrix
    seeds <- generate_seeds(
      r_in, n_seeds = n_harvests,
      protect_border = protect_border,
      patch_size = harvest_size
    )$seeds
  }

  patch_area_target <- round(n_harvests * harvest_size * (1 + percent_retention))

  r <- add_patches_parallel(
    x            = r_in,
    seeds        = seeds,
    patch_size   = round(harvest_size * (1 + percent_retention)),
    forest_code  = forest_code,
    harvest_code = harvest_code
  )

  patch_area_generated <- sum(terra::values(r) == harvest_code, na.rm = TRUE)
  if (patch_area_generated < patch_area_target) {
    shortfall <- patch_area_target - patch_area_generated
    r <- landscapeR::expandClass(
      context = r,
      class   = harvest_code,
      size    = shortfall,
      bgr     = forest_code
    )
  }

  seeds_used <- attr(r, "patch_seeds")

  island_size <- round(harvest_size * percent_retention)
  if (island_size > 0) {
    r <- landscapeR::makeClass(
      context = r,
      npatch  = n_harvests,
      size    = island_size,
      pts     = seeds_used,
      bgr     = harvest_code,
      val     = island_code,
      rast    = TRUE
    )
  }

  attr(r, "patch_seeds") <- seeds_used
  r
}

#' Apply an edge-effect buffer outward from harvest cells into forest cells
add_edge_effect <- function(r,
                            harvest_code   = 20,
                            edge_code      = 40,
                            edge_cell_dist = 1) {
  r_harvest <- r
  r_harvest[r_harvest != harvest_code] <- 9999
  r_dist  <- terra::distance(r_harvest, target = 9999)
  r_edge  <- r_dist |>
    terra::classify(rcl = matrix(c(0, edge_cell_dist, 0), ncol = 3), others = 1) |>
    (\(x) x * r)() |>
    terra::classify(rcl = matrix(c(0, edge_code), ncol = 2))
  r_edge
}

#' Count raster cells by class, returning a wide one-row tibble
count_classes <- function(harvest_raster,
                          forest_code  = 10,
                          harvest_code = 20,
                          island_code  = 30,
                          edge_code    = 40) {
  counts <- dplyr::tibble(value = terra::values(harvest_raster)) |>
    dplyr::mutate(
      class = dplyr::case_match(
        value,
        forest_code  ~ "Forest",
        harvest_code ~ "Harvest",
        island_code  ~ "Forest_Retention",
        edge_code    ~ "Edge"
      )
    ) |>
    dplyr::count(class, name = "n")

  dplyr::tibble(class = c("Forest", "Harvest", "Forest_Retention", "Edge")) |>
    dplyr::left_join(counts, by = "class") |>
    dplyr::mutate(n = dplyr::coalesce(n, 0L)) |>
    tidyr::pivot_wider(names_from = class, values_from = n)
}

#' Apply factor labels to numeric raster class codes
label_classes <- function(r,
                          codes  = c(10, 20, 30, 40),
                          labels = c("Forest", "Harvest", "Retention", "Edge")) {
  stopifnot(length(codes) == length(labels))
  r  <- terra::as.factor(r)
  lv <- data.frame(value = codes, class = labels)
  levels(r) <- lv
  r
}

# Geometric analysis -------------------------------------------------------

#' Exterior edge-effect loss area (Steiner formula approximation)
exterior_edge_loss <- function(area, edge_dist, shape_index = 2 * sqrt(pi)) {
  P <- sqrt(area) * shape_index
  P * edge_dist + pi * edge_dist^2
}

#' Interior edge-effect loss area for retention islands
interior_edge_loss <- function(area, edge_dist) {
  inradius <- sqrt(area / pi)
  out      <- numeric(length(area))
  idx      <- edge_dist < inradius
  out[idx]  <- 2 * edge_dist[idx] * sqrt(pi * area[idx]) - pi * edge_dist[idx]^2
  out[!idx] <- area[!idx]
  out
}

#' Generate theoretical edge-loss predictions across a retention gradient
make_example_data <- function(harvest_area  = 25 * 10^4,
                               edge_distance = c(10, 25, 50, 75, 100),
                               psi           = 2 * sqrt(pi),
                               xmin          = 0,
                               xmax          = 0.2) {
  p_retention <- seq(xmin, xmax, by = 0.001)

  dplyr::as_tibble(tidyr::expand_grid(
    harvest_area = harvest_area,
    p_retention  = p_retention,
    edge_dist    = edge_distance
  )) |>
    dplyr::mutate(
      harvest_footprint          = harvest_area * (1 + p_retention),
      retention_area             = harvest_area * p_retention,
      baseline_edge_loss         = exterior_edge_loss(harvest_area,     edge_dist, shape_index = psi),
      exterior_edge_loss         = exterior_edge_loss(harvest_footprint, edge_dist, shape_index = psi),
      interior_edge_loss         = interior_edge_loss(retention_area,    edge_dist),
      pct_agg_loss               = dplyr::if_else(retention_area != 0,
                                                  interior_edge_loss / retention_area,
                                                  NA_real_),
      total_edge_loss            = exterior_edge_loss + interior_edge_loss,
      additional_edge_loss       = total_edge_loss - baseline_edge_loss,
      additional_edge_loss_per_H = additional_edge_loss / harvest_area,
      additional_edge_loss_per_R = additional_edge_loss / retention_area
    )
}
