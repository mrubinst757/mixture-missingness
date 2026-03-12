# Sensitivity Analysis for Causal Inference with Informative Censoring

This repository contains R code for conducting sensitivity analysis when estimating causal effects in the presence of informative censoring. The methods allow researchers to bound treatment effects under various assumptions about the censoring mechanism and outcome distributions.

## Overview

When missing data is not missing at random (MNAR), standard causal inference methods produce biased estimates. This code implements:

- **Assumption-free bounds**: Sharp bounds on the Average Treatment Effect (ATE) without additional assumptions
- **Bounds under monotonicity**: Tighter bounds when outcomes are monotone in censoring status
- **Parameterized sensitivity analysis**: Bounds under user-specified sensitivity parameters:
  - `delta`: Fraction of censoring that is informative (vs. uninformative)
  - `tau`: Risk ratio between censored and uncensored subpopulations
- **Alternative estimands**: Z-ATE (composite outcome) and Separable Direct Effects (SDE)
- **Tipping point analysis**: Identifying sensitivity parameter values where conclusions would change

## Repository Structure

```
.
├── R/
│   ├── sim-fns.R        # Simulation data generation functions
│   ├── bound-fns.R      # Bound computation functions
│   ├── plot-fns.R       # Visualization functions
│   └── longdat.R        # Longitudinal data utilities
├── 00-simulation-specs.R    # Define simulation specifications
├── 01-run-sims.R            # Run Monte Carlo simulations
├── 02-make-plots.R          # Generate plots for paper
├── 03-calc-estimands.R      # Calculate estimands (worked examples)
├── 04-running-example.R     # Detailed running example with explanations
├── sim-inputs/              # Generated simulation inputs (created by script 00)
├── _output/                 # Simulation results (created by script 01)
└── plots/                   # Generated plots (created by script 02)
```

## Getting Started

### Dependencies

Install required R packages:

```r
install.packages(c("dplyr", "tidyr", "purrr", "ggplot2", "here",
                   "randplot", "ggrepel", "broom", "knitr", "kableExtra"))
```

For the full simulations with machine learning-based nuisance estimation:
```r
install.packages("SuperLearner")
install.packages("ranger")  # For random forest
```

### Quick Start: Running Example

The best place to start is [`04-running-example.R`](04-running-example.R), which provides a detailed worked example demonstrating:

- Data generation with informative censoring
- Computing true vs biased ATEs
- General (assumption-free) bounds
- Bounds under positive/negative monotonicity
- Bounds with bounded counterfactual outcome risk
- Tipping point analysis
- Sensitivity plot generation

Simply run:
```r
source("04-running-example.R")
```

## Running the Full Simulation Study

The simulation study evaluates finite-sample performance of the proposed bounds estimators. Run scripts in order:

### 1. Generate Simulation Specifications

```r
source("00-simulation-specs.R")
```

This creates population datasets under different data-generating processes and defines a grid of sensitivity parameter specifications. Outputs are saved to `sim-inputs/`.

### 2. Run Monte Carlo Simulations

```r
source("01-run-sims.R")
```

This runs two types of simulations:
- **Easy Simulation**: Known nuisance functions with added noise (tests at different convergence rates)
- **Hard Simulation**: Nuisance functions estimated via SuperLearner (more realistic)

Results are saved to `_output/simulations.rds` and `_output/simulations-sl.rds`.

**Note**: The full simulation study is computationally intensive and may take several hours to complete.

### 3. Generate Visualizations

```r
source("02-make-plots.R")
```

Creates publication-ready plots:
- **Plot 0**: Nuisance functions (propensity score, censoring, outcome models)
- **Plot 1**: True vs biased conditional treatment effects
- **Plot 2**: Assumption-free bounds
- **Plot 3**: Bounds under monotonicity assumptions
- **Plot 4**: Parameterized sensitivity analysis bounds

Plots are saved to `plots/`.

### 4. Calculate Estimands (Optional)

```r
source("03-calc-estimands.R")
```

Provides additional worked examples of calculating causal estimands and includes an interrupted time series illustration.

## Key Functions

### Bound Computation Functions (`R/bound-fns.R`)

- `ConstructPI()`: Compute plug-in estimates of bound components
- `CalculateTruth()`: Compute true values of estimands from population data
- `IterateSimParams()`: Compute bounds across multiple sensitivity parameter specifications

### Simulation Functions (`R/sim-fns.R`)

- `CreatePopulationData()`: Generate population data with informative censoring
- `EasySimulation()`: Run simulations with known nuisance functions
- `HardSimulation()`: Run simulations with estimated nuisance functions

### Plotting Functions (`R/plot-fns.R`)

- `Plot0()`: Visualize nuisance functions
- `Plot1()`: Compare true vs biased conditional treatment effects
- `Plot2()`: Show assumption-free bounds
- `Plot3()`: Show bounds under monotonicity
- `Plot4()`: Show parameterized sensitivity bounds

## Example: Computing Bounds

```r
# Load functions
source("R/bound-fns.R")

# Generate population data
set.seed(123)
data <- CreatePopulationData(
  population_size = 1e5,
  delta1 = 0.5,  # 50% of treatment arm censoring is informative
  delta0 = 0.5,  # 50% of control arm censoring is informative
  tau1 = 2,      # 2x risk ratio in treatment arm
  tau0 = 2       # 2x risk ratio in control arm
)

# Compute plug-in estimates
data_pi <- ConstructPI(data)

# Assumption-free bounds (no additional assumptions)
ate_lb_free <- with(data_pi, mean(phi.ate.pi - phi.11.pi - (phi.0.pi - phi.00.pi)))
ate_ub_free <- with(data_pi, mean(phi.ate.pi + (phi.1.pi - phi.11.pi) + phi.00.pi))
c(ate_lb_free, ate_ub_free)

# Bounds under positive monotonicity (delta = 0.8)
ate_lb_pos <- with(data_pi, mean(phi.ate.pi - 0.8 * (phi.0.pi - phi.00.pi)))
ate_ub_pos <- with(data_pi, mean(phi.ate.pi + 0.8 * (phi.1.pi - phi.11.pi)))
c(ate_lb_pos, ate_ub_pos)
```

## Sensitivity Parameters

### Delta (δ): Informative Censoring Fraction

- Range: [0, 1]
- Interpretation: Proportion of censoring that is informatively missing
- δ = 0: All censoring is uninformative (MAR holds)
- δ = 1: All censoring is informative (maximum bias)

### Tau (τ): Outcome Risk Ratio

- Range: (0, ∞), typically τ ≥ 1 for inequality bounds
- Interpretation: Ratio of outcome risk between censored and uncensored
- τ = 1: Outcomes are the same regardless of censoring type (MAR holds)
- τ > 1: Censored individuals have higher outcome risk (positive selection)
- τ < 1: Censored individuals have lower outcome risk (negative selection)

## Citation

This code generates the simulations and analysis for the accompanying research paper. If you use this code, please cite:

[Citation information to be added upon publication]

## License

[License information to be added]

## Contact

For questions or issues, please open an issue on this repository.
