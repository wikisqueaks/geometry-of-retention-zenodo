# Data and Code: The Geometry of Retention

This repository contains the simulation code and output data associated with:

> Taylor, K.W., Viliani, L., & Nielsen, S.E. (*in review*). The geometry of retention: paradoxical consequences of aggregated retention forestry.

**Authors:** Kyle W. Taylor, Leonardo Viliani, Scott E. Nielsen  
Department of Renewable Resources, University of Alberta, Edmonton, Alberta, Canada

**Funding:** Natural Sciences and Engineering Research Council of Canada Alliance Grant (NSERC, ALLRP 548285–19). Research conducted in partnership with the Boreal Ecosystem Recovery and Assessment (BERA) project (www.bera-project.org).

---

## Contents

```
.
├── scripts/
│   ├── 00_functions.R               # Shared functions (landscape generation & geometry)
│   ├── 01_landscape_simulation.R    # Runs the canonical landscape simulation
│   ├── 02_single_patch_outputs.R    # Analyses and figures for single-harvest scenario
│   └── 03_multi_patch_outputs.R     # Analyses and figures for multi-harvest scenario
└── data/
    ├── simulation_H1_N1000.csv      # Single-harvest simulation results (1000 replications)
    └── simulation_H8_N1000.csv      # Multi-harvest simulation results (1000 replications)
```

---

## Data Description

Both CSV files share the same schema:

| Column | Description |
|---|---|
| `sim_id` | Unique identifier encoding harvest count, edge distance, retention %, and replication |
| `rep_id` | Replication number (1–1000) |
| `edge_distance` | Edge-effect distance in raster cells (1 cell = 10 m; values: 1, 2.5, 5, 7.5, 10) |
| `percent_retention` | Retention as a proportion of harvest area (0–0.20) |
| `Forest` | Count of forest cells (class code 10) |
| `Harvest` | Count of harvested cells (class code 20) |
| `Forest_Retention` | Count of retention island cells (class code 30) |
| `Edge` | Count of edge-effect cells (class code 40) |

Landscape dimensions: 316 × 316 cells (each cell = 10 m × 10 m ≈ 1000 ha total).
Harvest size: 25 ha (2500 cells) per patch.

---

## Reproducing the Simulation

`01_landscape_simulation.R` re-runs the simulation from scratch. Results will be
written to `data/` with a timestamped filename. **Note:** the full simulation
(1000 replications × 9 retention levels × 5 edge distances) takes several hours
even with parallel processing. The pre-computed CSVs in `data/` are provided so
the analysis scripts can be run without re-running the simulation.

---

## Dependencies

### R packages

| Package | Role |
|---|---|
| `terra` | Raster landscape creation and manipulation |
| `landscapeR` | Patch generation within rasters |
| `tidyverse` | Data wrangling and plotting |
| `future` / `furrr` / `future.apply` | Parallel simulation |
| `kableExtra` | Summary tables |
| `scales` | Number formatting |
| `glue` | String interpolation |
| `here` | Project-relative file paths |

Install all at once:

```r
install.packages(c(
  "terra", "landscapeR", "tidyverse",
  "future", "furrr", "future.apply",
  "kableExtra", "scales", "glue", "here"
))
```

---

## Usage

Set the project root as your working directory (or open the `.Rproj` file), then
run the scripts in order:

```r
# Optional: re-run simulation (slow)
source("scripts/01_landscape_simulation.R")

# Analyse pre-computed results
source("scripts/02_single_patch_outputs.R")
source("scripts/03_multi_patch_outputs.R")
```

All scripts use `here::here()` for file paths and expect to be run from the
project root.

---

## Shared Functions (`scripts/00_functions.R`)

The shared functions were extracted from the `geomUtil` R package developed
alongside this project by K.W. Taylor (kwtaylor@ualberta.ca). They implement:

- **Landscape generation**: `create_landscape()`, `generate_seeds()`,
  `add_patches_parallel()`, `add_retention_harvest()`, `add_edge_effect()`,
  `count_classes()`, `label_classes()`
- **Geometric analysis**: `exterior_edge_loss()`, `interior_edge_loss()`,
  `make_example_data()`

---

## License

Code: [MIT License](LICENSE)  
Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
