# 03_multi_patch_outputs.R
#
# Analyzes simulation results for the multi-harvest (H8) scenario and compares
# empirical simulation results against theoretical "donut model" predictions.
# Produces:
#   - A table of theoretical predictions from the geometric model
#   - A summary table of simulated mean additional edge (%) with 95% CIs

library(tidyverse)
library(kableExtra)
library(scales)
library(glue)

source(here::here("scripts", "00_functions.R"))

sim_data <- read_csv(here::here("data", "simulation_H8_N1000.csv"))

# Theoretical predictions from the geometric (donut) model ------------------
targets <- c(0.05, 0.10, 0.15, 0.20)

make_example_data() |>
  filter(p_retention %in% targets) |>
  select(edge_dist, p_retention, additional_edge_loss_per_H) |>
  pivot_wider(
    id_cols     = "p_retention",
    names_from  = "edge_dist",
    values_from = "additional_edge_loss_per_H",
    names_glue  = "{.name} m"
  ) |>
  mutate(
    Retention = scales::percent(p_retention),
    .before   = everything(),
    .keep     = "unused"
  ) |>
  mutate(across(where(is.numeric), \(x) round(x * 100, 1))) |>
  kbl(caption = "Theoretical additional edge-effect loss per unit harvest (%)")

# Compute additional edge relative to the 0%-retention baseline ---------------
edge_change_data <- sim_data |>
  mutate(
    edge_distance     = edge_distance * 10,
    Total_Forest      = Forest + Forest_Retention,
    # Round to avoid floating-point mismatches in downstream filters
    percent_retention = round(percent_retention, 4)
  ) |>
  group_by(rep_id, edge_distance) |>
  mutate(
    baseline_edge         = Edge[percent_retention == 0],
    additional_edge       = Edge - baseline_edge,
    additional_edge_per_H = additional_edge / Harvest
  ) |>
  ungroup()

# Summary table: mean and 95% CI of additional edge per unit harvest ----------
ci95 <- edge_change_data |>
  filter(percent_retention %in% c(0, targets)) |>
  group_by(edge_distance, rep_id) |>
  mutate(y0 = additional_edge_per_H[percent_retention == 0]) |>
  filter(percent_retention %in% targets) |>
  mutate(delta = additional_edge_per_H - y0) |>
  ungroup() |>
  group_by(edge_distance, percent_retention) |>
  summarise(
    mean = mean(delta),
    lwr  = quantile(delta, 0.025),
    upr  = quantile(delta, 0.975),
    .groups = "drop"
  )

ci95 |>
  mutate(
    percent_retention = scales::percent(percent_retention, suffix = ""),
    across(c(lwr, upr), \(x) format(round(x * 100, 1), nsmall = 1, trim = TRUE)),
    mean   = scales::percent(mean, accuracy = 0.1, suffix = ""),
    mean   = paste0("<b>", mean, "</b>"),
    conf95 = glue("({lwr}\u2013{upr})")
  ) |>
  pivot_wider(
    id_cols     = "percent_retention",
    names_from  = "edge_distance",
    values_from = c("mean", "conf95")
  ) |>
  select(1, 2, 7, 3, 8, 4, 9, 5, 10, 6, 11) |>
  kbl(
    align     = rep("r", 11),
    col.names = c("Retention (%)", rep(c("Mean", "95% CI"), 5)),
    escape    = FALSE
  ) |>
  add_header_above(c(" " = 1, "10 m" = 2, "25 m" = 2, "50 m" = 2, "75 m" = 2, "100 m" = 2)) |>
  add_header_above(c(" " = 1, "Distance of Edge Effect" = 10)) |>
  kable_material()

