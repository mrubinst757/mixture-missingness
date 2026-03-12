################################################################################
############################ RUN SIMULATIONS ###################################
################################################################################
#
# This script runs Monte Carlo simulations to evaluate the finite-sample
# performance of the proposed bounds estimators. Two simulation designs:
#
# 1. Easy Simulation: Known nuisance functions with added noise
#    - Tests at different convergence rates (alpha = 0.3 and 0.1)
#    - Faster to run, good for testing estimation procedure
#
# 2. Hard Simulation: Nuisance functions estimated via SuperLearner
#    - More realistic assessment of finite-sample performance
#    - Accounts for nuisance function estimation error
#
################################################################################

# Load required packages and functions
library(tidyverse); library(knitr); library(kableExtra)
source("R/sim-fns.R"); source("R/bound-fns.R")

# Load population datasets
data = readRDS("sim-inputs/data.rds")

################################################################################
##################### Simulation One: Simulated Eta ############################
################################################################################
# Run simulations with known nuisance functions plus simulated estimation error
# Using dataset 4 from the data list

set.seed(145)

# Simulation with alpha = 0.3 convergence rate (slower convergence)
test1 = EasySimulation(data$data[[4]], 1000, 1000,
                      mu_rate = 0.3, pi_rate = 0.3, gamma_rate = 0.3)

# Simulation with alpha = 0.1 convergence rate (faster convergence)
test2 = EasySimulation(data$data[[4]], 1000, 1000,
                       mu_rate = 0.1, pi_rate = 0.1, gamma_rate = 0.1)

# Save simulation results
saveRDS(list(test1, test2), "_output/simulations.rds")

# Process simulation results and create summary table
table_input = map(list(test1, test2), ~.x %>%
  # Combine all simulation replications
  invoke(rbind, .) %>%
  # Calculate performance metrics for each bound component
  group_by(Estimand = component) %>%
  summarize(Bias = mean(ests - truth),                    # Average bias
            RMSE = sqrt(mean((ests - truth)^2)),          # Root mean squared error
            Coverage = mean(covered))) %>%                # 95% CI coverage rate
  # Add convergence rate identifier
  map2(c(0.3, 0.1), ~mutate(.x, Alpha = .y)) %>%
  # Combine both convergence rates
  invoke(rbind, .) %>%
  # Reshape for table presentation
  nest(res = c(Bias, RMSE, Coverage)) %>%
  pivot_wider(names_from = Alpha, values_from = res) %>%
  unnest()


################################################################################
##################### Simulation Two: Estimated Eta ############################
################################################################################
# Run simulations with SuperLearner-estimated nuisance functions
# This is more computationally intensive but more realistic

set.seed(155)

test3 = HardSimulation(popdata = data$data[[4]],
                       sample_size = 2000,           # Larger sample for ML
                       number_sims = 1000,           # 1000 replications
                       sl_lib = c("SL.glm", "SL.ranger"))  # Use GLM and Random Forest

saveRDS(test3, "_output/simulations-sl.rds")

# Process SuperLearner simulation results
table_input1 = test3 %>%
  invoke(rbind, .) %>%
  # Filter out extreme variance estimates (likely due to convergence issues)
  filter(varests < 0.25) %>%
  group_by(Estimand = component) %>%
  summarize(Bias2 = mean(ests - truth),
            RMSE2 = sqrt(mean((ests - truth)^2)),
            Coverage2 = mean(covered))

# Combine all results into publication-ready table
table_input %>%
  left_join(table_input1) %>%
  # Round to 2 decimal places
  mutate_if(is.numeric, ~round(., 2)) %>%
  # Reorder rows (put psi_tilde first)
  slice(c(6, 1:5)) %>%
  # Generate LaTeX table
  knitr::kable(format = "latex", booktabs = TRUE, escape = FALSE, align = "c",
               caption = "Simulation results \\label{tab:1}",
               col.names = c("Estimand", rep(c("Bias", "RMSE", "Coverage"), 3))) %>%
  # Add header row for different simulation conditions
  add_header_above(header = c(" " = 1,
                              "$\\\\alpha = 0.3$" = 3,
                              "$\\\\alpha = 0.1$" = 3,
                              "SuperLearner" = 3), escape = FALSE)

