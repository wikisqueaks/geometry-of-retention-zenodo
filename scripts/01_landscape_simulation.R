# 01_landscape_simulation.R
#
# Runs the canonical landscape simulation for retention harvest edge effects.
# Produces a CSV of cell-class counts across a grid of retention percentages,
# edge-effect distances, and replication seeds.
#
# Each replication generates one random landscape seed, then applies all
# parameter combinations to that same seed -- isolating the effect of
# retention from landscape heterogeneity.
#
# Cell size: 10 m x 10 m (100 m2); 1 ha = 100 cells
# Landscape: 316 x 316 cells (~1000 ha)
# Harvest size: 25 ha (2500 cells)
#
# Outputs a timestamped CSV to data/

library(terra)
library(future)
library(future.apply)
library(furrr)
library(parallel)
library(tidyverse)
library(landscapeR)
library(glue)

source(here::here("scripts", "00_functions.R"))

retention_landscape_simulation <- function(
    replications   = 100,
    edge_distances = c(1, 2.5, 5, 7.5, 10),
    retention_min  = 0,
    retention_max  = 0.20,
    retention_step = 0.025,
    n_harvests     = 7,
    harvest_size   = 25 * 100,
    ls_dim         = 316,
    csv_output_dir = here::here("data"),
    tif_output_dir = NULL
) {
  percent_retention <- seq(
    from = retention_min,
    to   = retention_max,
    by   = retention_step
  )

  params <- expand.grid(
    edge_distance     = edge_distances,
    percent_retention = percent_retention,
    n_harvests        = n_harvests,
    harvest_size      = harvest_size
  )

  sim_name <- glue(
    "{runtime}_{harvests}_{retentions}_{edges}_N{replications}",
    runtime    = format(Sys.time(), "%Y%m%d%H%M%S"),
    harvests   = paste0("H", n_harvests),
    retentions = paste0("R", format(retention_max * 100), "b", format(retention_step * 1000)),
    edges      = paste0("E", paste0(edge_distances * 10, collapse = "."))
  )

  csv_output_path <- file.path(csv_output_dir, paste0(sim_name, ".csv"))
  dir.create(csv_output_dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.null(tif_output_dir)) {
    dir.create(tif_output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Parallel processing: replications are distributed across all available cores.
  # Each worker handles one replication group independently. The number of workers
  # is set automatically to the number of logical cores on the machine.
  # To limit core usage, replace detectCores() with a fixed integer (e.g., 4).
  plan(multisession, workers = parallel::detectCores())
  message("Running simulation: ", sim_name)

  results_df <- params |>
    slice(rep(1:n(), each = replications)) |>
    mutate(
      rep_id = rep(1:replications, times = nrow(params)),
      sim_id = paste0(
        "S", row_number(),
        "-E", edge_distance,
        "-R", 100 * percent_retention,
        "-REP", rep_id
      )
    ) |>
    group_split(rep_id) |>
    future_map_dfr(
      \(rep_group) {
        seed_data <- create_landscape(dimension = ls_dim) |>
          generate_seeds(
            n_seeds        = n_harvests,
            protect_border = TRUE,
            patch_size     = harvest_size,
            rng_seed       = rep_group$rep_id[1]
          )

        map_dfr(seq_len(nrow(rep_group)), \(i) {
          p_row <- rep_group[i, ]

          r <- add_retention_harvest(
            landscape_matrix  = seed_data,
            n_harvests        = p_row$n_harvests,
            harvest_size      = p_row$harvest_size,
            percent_retention = p_row$percent_retention
          ) |>
            add_edge_effect(edge_cell_dist = p_row$edge_distance)

          if (!is.null(tif_output_dir)) {
            out_path <- file.path(tif_output_dir, paste0(p_row$sim_id, ".tif"))
            terra::writeRaster(r, out_path, overwrite = TRUE)
          }

          metrics <- count_classes(r)
          cbind(
            sim_id            = p_row$sim_id,
            rep_id            = p_row$rep_id,
            edge_distance     = p_row$edge_distance,
            percent_retention = p_row$percent_retention,
            metrics
          )
        })
      },
      .options = furrr_options(seed = TRUE)
    )

  plan(sequential)

  write_csv(results_df, csv_output_path)
  message("Done. Results written to: ", csv_output_path)
}

# Run simulations -----------------------------------------------------------
# Single harvest (H1), 1000 replications
retention_landscape_simulation(
  replications = 1000,
  n_harvests   = 1
)

# Multi-harvest (H8), 1000 replications
retention_landscape_simulation(
  replications = 1000,
  n_harvests   = 8
)
