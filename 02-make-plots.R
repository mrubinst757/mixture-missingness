################################################################################
########################## GENERATE PLOTS ######################################
################################################################################
#
# This script generates visualizations for the paper:
# - Plot 0: Nuisance functions (propensity score, censoring, outcome models)
# - Plot 1: True vs biased conditional treatment effects
# - Plot 2: Assumption-free bounds
# - Plot 3: Bounds under monotonicity assumptions
# - Plot 4: Bounds under parameterized sensitivity analysis
#
# Each plot is generated for multiple datasets and saved to plots/ directory
#
################################################################################

# Load required functions and packages
source("R/plot-fns.R"); source("R/bound-fns.R"); library(ggrepel)

# Load population datasets
data = readRDS("sim-inputs/data.rds")

set.seed(123)

# Prepare data for plotting: compute plug-in estimates and sample for visualization
plot_data = data %>%
  mutate(data = map(data, ConstructPI)) %>%    # Add plug-in bound components
  mutate(data = map(data, PlotData))           # Sample and compute conditional estimands

################################################################################
### Plot 0: Visualize nuisance functions
################################################################################
p0s = map(plot_data$data, Plot0)
ggsave("plots/plot-0.png", p0s[[7]], height = 5, width = 9)

################################################################################
### Plot 1: Compare true vs biased conditional treatment effects
################################################################################
# Generate separate plots for ATE and for Z-ATE/SDE
p1s = map(plot_data$data, ~Plot1(.x, "ATE"))
p1sa = map(plot_data$data, ~Plot1(.x, c("ATE-Z", "SDE")))

# Save plot for dataset 7
ggsave("plots/plot-1.png", p1s[[7]], height = 5, width = 9)
ggsave("plots/plot-1a.png", p1sa[[7]], height = 5, width = 9)

################################################################################
### Plot 2: Assumption-free bounds
################################################################################
p2s = map(plot_data$data, ~Plot2(.x, "ATE"))
ggsave("plots/plot-2.png", p2s[[7]], height = 5, width = 9)

################################################################################
### Plot 3: Bounds under monotonicity assumptions
################################################################################
p3s = map(plot_data$data, Plot3)
ggsave("plots/plot-3.png", p3s[[7]], height = 5, width = 9)
################################################################################
### Numerical examples: compute specific bound values
################################################################################
dat = plot_data$data[[7]]

# Example 1: Bounds under positive monotonicity with delta = 1 (no censoring attenuation)
ate.lb.pos = mean(with(dat, phi.ate.pi - (phi.0.pi - phi.00.pi)))
ate.ub.pos = mean(with(dat, phi.ate.pi + (phi.1.pi - phi.11.pi)))
ate.lb.neg = mean(with(dat, phi.ate.pi - phi.11.pi))
ate.ub.neg = mean(with(dat, phi.ate.pi + phi.00.pi))

# Example 2: Bounds under monotonicity with delta = 0.7
ate.lb.pos7 = mean(with(dat, phi.ate.pi - 0.7 * (phi.0.pi - phi.00.pi)))
ate.ub.pos7 = mean(with(dat, phi.ate.pi + 0.7 * (phi.1.pi - phi.11.pi)))
ate.lb.neg7 = mean(with(dat, phi.ate.pi - 0.7 * phi.11.pi))
ate.ub.neg7 = mean(with(dat, phi.ate.pi + 0.7 * phi.00.pi))
c(ate.lb.neg7, ate.ub.neg7)
c(ate.lb.pos7, ate.ub.pos7)

# Example 3: Assumption-free bounds
ate.ub = mean(with(dat, phi.ate.pi + (phi.1.pi - phi.11.pi) + phi.00.pi))
ate.lb = mean(with(dat, phi.ate.pi - phi.11.pi - (phi.0.pi - phi.00.pi)))

################################################################################
### Plot 4: Parameterized sensitivity analysis bounds
################################################################################
# Compare assumption-free bounds to bounds under specified sensitivity parameters
p4s = map(plot_data$data, ~Plot4(.x, tau1 = 2, tau0 = 2,
                                 delta1 = 0.7, delta0 = 0.7,
                                 delta1l = 0.3, delta0l = 0.3,
                                 group = "ATE"))

p4sa = map(plot_data$data, ~Plot4(.x, tau1 = 2, tau0 = 2,
                                 delta1 = 0.7, delta0 = 0.7,
                                 delta1l = 0.3, delta0l = 0.3,
                                 group = c("ATE-Z", "SDE")))

p4s[[7]]
p4sa[[7]]

# Example 4: Bounds with asymmetric delta values and tau = 2
delta1 = 0.7; delta0 = 0.3; tau1 = tau0 = 2
ate.ub = mean(with(dat, phi.ate.pi + delta1 * phi.11.pi * (tau1 - 1) - delta0 * phi.00.pi * ((1 - tau0) / tau0)))
ate.lb = mean(with(dat, phi.ate.pi + delta1 * phi.11.pi * ((1 - tau1) / tau1) - delta0 * phi.00.pi * (tau0 - 1)))

# Save plots
ggsave("plots/plot-4.png", p4s[[7]], height = 5, width = 9)
ggsave("plots/plot-4a.png", p4sa[[7]], height = 5, width = 9)

