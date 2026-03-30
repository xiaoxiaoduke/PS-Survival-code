# PS-Survival-Code

Demo code for the paper:
> "Implementing the principal stratum strategy for intercurrent events with survival outcomes: a tutorial"

## File Description

| File | Description |
|------|-------------|
| `util.R` | Helper functions used by the simulation |
| `generate_data.R` | Generates synthetic data for use in `example.R` |
| `example.R` | (1) Mixture modeling approach via R package `PStrata`; (2) Weighting approach via R package `mrPStrata` |

## Usage

Run scripts in this order:
1. `generate_data.R`
2. `example.R` (sources `util.R` automatically)

## Requirements

- R (≥ 4.0)
- Packages: `PStrata`, `mrPStrata`, `survival`, `tidyverse`

## Citation

