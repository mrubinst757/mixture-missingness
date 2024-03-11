library(tidyverse); library(knitr); library(kableExtra)
source("R/sim-fns.R"); source("R/bound-fns.R")

data = readRDS("sim-inputs/data.rds")

################################################################################
##################### Simulation One: Simulated Eta ############################
################################################################################

set.seed(145)
test1 = EasySimulation(data$data[[4]], 1000, 1000,
                      mu_rate = 0.3, pi_rate = 0.3, gamma_rate = 0.3)

test2 = EasySimulation(data$data[[4]], 1000, 1000,
                       mu_rate = 0.1, pi_rate = 0.1, gamma_rate = 0.1)

table_input = map(list(test1, test2), ~.x %>%
  invoke(rbind, .) %>%
  group_by(Estimand = component) %>%
  summarize(Bias = mean(ests - truth),
            RMSE = sqrt(mean((ests - truth)^2)),
            Coverage = mean(covered))) %>%
  map2(c(0.3, 0.1), ~mutate(.x, Alpha = .y)) %>%
  invoke(rbind, .) %>%
  nest(res = c(Bias, RMSE, Coverage)) %>%
  pivot_wider(names_from = Alpha, values_from = res) %>%
  unnest()

table_input %>%
  mutate_if(is.numeric, ~round(., 2)) %>%
  slice(c(5, 2, 3, 1, 4)) %>%
  knitr::kable(format = "latex", booktabs = TRUE, escape = FALSE, align = "c",
               caption = "Simulation results \\label{tab:1}",
               col.names = c("Estimand", rep(c("Bias", "RMSE", "Coverage"), 2))) %>%
  add_header_above(header = c(" " = 1, "$\\\\alpha = 0.3$" = 3, "$\\\\alpha = 0.1$" = 3), escape = FALSE)

################################################################################
##################### Simulation One: Estimated Eta ############################
################################################################################
set.seed(155)
test3 = HardSimulation(popdata = data$data[[4]],
                       sample_size = 2000,
                       number_sims = 1000,
                       sl_lib = c("SL.glm", "SL.ranger"))

table_input1 = test3 %>%
  invoke(rbind, .) %>%
  filter(varests < 0.25) %>%
  group_by(Estimand = component) %>%
  summarize(Bias = mean(ests - truth),
            RMSE = sqrt(mean((ests - truth)^2)),
            Coverage = mean(covered)) %>%
  slice(c(5, 2, 3, 1, 4)) %>%
  mutate_if(is.numeric, ~round(., 2))

table_input1 %>%
  knitr::kable(format = "latex", booktabs = TRUE, escape = FALSE, align = "c",
               caption = "Simulation results: SuperLearner \\label{tab:2}")
