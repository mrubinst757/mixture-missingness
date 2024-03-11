source("R/sim-fns.R")
source("R/bound-fns.R")

################################################################################
###################### STEP ONE: INITIALIZE DATASETS ###########################
################################################################################
data_list = expand.grid(
  delta1 = c(0.5, 0.2),
  delta0 = c(0.5, 0.2),
  tau1 = c(2, 0.5),
  tau0 = c(2, 0.5)
) %>%
  filter((delta1 == delta0 & tau1 == tau0) | (delta1 != delta0 & tau1 != tau0)) %>%
  as_tibble() %>%
  set_names(paste0(names(.), "_true")) %>%
  mutate(data = pmap(list(delta1 = delta1_true,
                          delta0 = delta0_true,
                          tau1 = tau1_true,
                          tau0 = tau0_true),
                     ~CreatePopulationData(population_size = 1e5, ...)))

saveRDS(data_list, "sim-inputs/data.rds")

################################################################################
################# STEP TWO: INITIALIZE ASSUMED SENS PARAMS #####################
################################################################################

base_specs = tibble(
  estimand = c(rep("ate", 5), rep("z-ate", 2), rep("sde", 2)),
  monotonicity = c("none", "positive", "negative", rep("none", 6)),
  outcome_assumption = c(rep("none", 3), "equality", "inequality", rep("none", 4)),
  bound = c(rep(TRUE, 3), FALSE, rep(TRUE, 2), FALSE, TRUE, FALSE)
)

param_spec = expand.grid(
  delta1 = c(0.8, 0.2),
  delta0 = c(0.8, 0.2),
  tau1 = c(0.5, 2)
) %>%
  mutate(tau0 = tau1) %>%
  pmap(~mutate(base_specs, delta1 = ..1, delta0 = ..2, tau1 = ..3, tau0 = ..4)) %>%
  map(~mutate(.x, tau1 = if_else(tau1 < 1 & outcome_assumption != "equality", 1 / tau1, tau1))) %>%
  map(~mutate(.x, tau0 = if_else(tau0 < 1 & outcome_assumption != "equality", 1 / tau0, tau0))) %>%
  invoke(rbind, .) %>%
  mutate_at(vars(matches("tau")), ~case_when(estimand != "ate" ~ NA_real_,
                                             estimand == "ate" ~ .)) %>%
  bind_rows(
    tibble(
      estimand = c("ate", "ate"),
      monotonicity = c("none", "none"),
      outcome_assumption = rep("equality", 2),
      tau1 = c(0.5, 2), tau0 = c(2, 0.5), delta1 = 0.8, delta0 = 0.8
    )
  ) %>%
  distinct()

################################################################################
##################### STEP THREE: COMBINE ALL SPECS  ###########################
################################################################################
full_specs = map(data_list$data, ~IterateSimParams(param_spec, ConstructPI(.x))) %>%
  map2(
    map(1:nrow(data_list), ~slice(data_list %>% select(-data), .x)),
    bind_cols
  ) %>%
  map2(
    map(data_list$data, CalculateTruth),
    left_join
  ) %>%
  map(~unnest_wider(.x, estimands, names_sep = "")) %>%
  map(~.x %>%
        mutate(true_params = case_when(
          delta1 == delta1_true & delta0 == delta0_true & tau1 == tau1_true & tau0 == tau0_true ~ "all",
          delta1 == delta1_true & delta0 == delta0_true ~ "pi_only",
          TRUE ~ "no"
        )) %>%
        mutate(motonicity_true = case_when(
          estimand == "ate" & tau1_true > 1 & tau0_true > 1 ~ "positive",
          estimand == "ate" & tau1_true < 1 & tau0_true < 1 ~ "negative",
          estimand == "ate" & tau1_true > 1 & tau0_true < 1 ~ "no",
          estimand == "ate" & tau1_true < 1 & tau0_true > 1 ~ "no",
          estimand != "ate" ~ NA
        ))
      )

saveRDS(full_specs, "sim-inputs/specs.rds")
