################################################################################
#################### SIMULATION SPECIFICATION SETUP ############################
################################################################################
#
# This script sets up the specifications for the simulation study by:
# 1. Generating population datasets with different data generating parameters
# 2. Creating a grid of sensitivity parameter specifications to evaluate
# 3. Combining datasets with specifications to create full simulation inputs
#
# The output is saved to sim-inputs/ directory for use by simulation scripts
#
################################################################################

# Load required functions
source("R/sim-fns.R")
source("R/bound-fns.R")

################################################################################
###################### STEP ONE: INITIALIZE DATASETS ###########################
################################################################################
# Generate population datasets under different data generating processes
# Each dataset has different true values of delta and tau parameters

data_list = expand.grid(
  delta1 = c(0.5, 0.2),  # Proportion of censoring that is informative under treatment
  delta0 = c(0.5, 0.2),  # Proportion of censoring that is informative under control
  tau1 = c(2, 0.5),      # Risk ratio mu_1(X,1)/mu_1(X,0) under treatment
  tau0 = c(2, 0.5)       # Risk ratio mu_0(X,1)/mu_0(X,0) under control
) %>%
  # Keep only scenarios where parameters match or differ across treatment arms
  # Avoids redundant scenarios where only one arm's parameters differ
  filter((delta1 == delta0 & tau1 == tau0) | (delta1 != delta0 & tau1 != tau0)) %>%
  as_tibble() %>%
  # Add "_true" suffix to distinguish true DGP values from assumed values
  set_names(paste0(names(.), "_true")) %>%
  # Generate population data for each parameter combination
  mutate(data = pmap(list(delta1 = delta1_true,
                          delta0 = delta0_true,
                          tau1 = tau1_true,
                          tau0 = tau0_true),
                     ~CreatePopulationData(population_size = 1e5, ...)))

# Save population datasets for use in simulations
saveRDS(data_list, "sim-inputs/data.rds")

################################################################################
################# STEP TWO: INITIALIZE ASSUMED SENS PARAMS #####################
################################################################################
# Create grid of sensitivity analysis specifications
# These represent different assumptions an analyst might make

# Base specifications: defines which estimand and assumptions to evaluate
base_specs = tibble(
  estimand = c(rep("ate", 5), rep("z-ate", 2), rep("sde", 2)),
  # Monotonicity assumptions (for ATE only)
  monotonicity = c("none", "positive", "negative", rep("none", 6)),
  # Outcome model assumptions (for ATE only)
  outcome_assumption = c(rep("none", 3), "equality", "inequality", rep("none", 4)),
  # Whether computing bounds (TRUE) or point-identified parameter (FALSE)
  bound = c(rep(TRUE, 3), FALSE, rep(TRUE, 2), FALSE, TRUE, FALSE)
)

# Create grid of assumed sensitivity parameter values to evaluate
param_spec = expand.grid(
  delta1 = c(0.8, 0.2),  # Assumed bounds on informative censoring fraction
  delta0 = c(0.8, 0.2),
  tau1 = c(0.5, 2)       # Assumed outcome risk ratios
) %>%
  mutate(tau0 = tau1) %>%  # Assume symmetric tau across treatment arms
  # Replicate base_specs for each parameter combination
  pmap(~mutate(base_specs, delta1 = ..1, delta0 = ..2, tau1 = ..3, tau0 = ..4)) %>%
  # For inequality bounds, ensure tau >= 1 (take reciprocal if needed)
  map(~mutate(.x, tau1 = if_else(tau1 < 1 & outcome_assumption != "equality", 1 / tau1, tau1))) %>%
  map(~mutate(.x, tau0 = if_else(tau0 < 1 & outcome_assumption != "equality", 1 / tau0, tau0))) %>%
  invoke(rbind, .) %>%
  # Set tau to NA for estimands that don't use it
  mutate_at(vars(matches("tau")), ~case_when(estimand != "ate" ~ NA_real_,
                                             estimand == "ate" ~ .)) %>%
  # Add special cases with asymmetric tau values
  bind_rows(
    tibble(
      estimand = c("ate", "ate"),
      monotonicity = c("none", "none"),
      outcome_assumption = rep("equality", 2),
      tau1 = c(0.5, 2), tau0 = c(2, 0.5), delta1 = 0.8, delta0 = 0.8
    )
  ) %>%
  distinct()  # Remove any duplicate specifications

################################################################################
##################### STEP THREE: COMBINE ALL SPECS  ###########################
################################################################################
# For each population dataset, compute bounds under all specifications
# and add metadata about whether assumptions match the truth

full_specs = map(data_list$data, ~IterateSimParams(param_spec, ConstructPI(.x))) %>%
  # Add true DGP parameters
  map2(
    map(1:nrow(data_list), ~slice(data_list %>% select(-data), .x)),
    bind_cols
  ) %>%
  # Add true values of estimands
  map2(
    map(data_list$data, CalculateTruth),
    left_join
  ) %>%
  # Unnest the computed bounds/estimates
  map(~unnest_wider(.x, estimands, names_sep = "")) %>%
  # Add indicators for whether assumed parameters match true parameters
  map(~.x %>%
        mutate(true_params = case_when(
          # All sensitivity parameters correctly specified
          delta1 == delta1_true & delta0 == delta0_true & tau1 == tau1_true & tau0 == tau0_true ~ "all",
          # Only censoring parameters correctly specified
          delta1 == delta1_true & delta0 == delta0_true ~ "pi_only",
          # Misspecified
          TRUE ~ "no"
        )) %>%
        # Indicator for whether true DGP satisfies monotonicity assumption
        mutate(motonicity_true = case_when(
          estimand == "ate" & tau1_true > 1 & tau0_true > 1 ~ "positive",
          estimand == "ate" & tau1_true < 1 & tau0_true < 1 ~ "negative",
          estimand == "ate" & tau1_true > 1 & tau0_true < 1 ~ "no",
          estimand == "ate" & tau1_true < 1 & tau0_true > 1 ~ "no",
          estimand != "ate" ~ NA_character_
        ))
      )

# Save full simulation specifications
saveRDS(full_specs, "sim-inputs/specs.rds")
