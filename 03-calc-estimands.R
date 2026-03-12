################################################################################
##################### CALCULATE ESTIMANDS EXAMPLE ##############################
################################################################################
#
# This script provides worked examples of calculating causal estimands
# and their biased counterparts in a simulated population. It demonstrates:
# - Computing true ATE vs observed (biased) ATE under MAR
# - Examining censoring rates and conditional outcome probabilities
# - Computing outcome risk ratios
#
################################################################################

# Load required functions
source("R/plot-fns.R"); source("R/bound-fns.R"); library(ggrepel)

# Load a single population dataset (dataset 7)
data = readRDS("sim-inputs/data.rds")$data[[7]]

################################################################################
### Calculate true and observed estimands
################################################################################

# True Average Treatment Effect (ground truth)
ate = with(data, mean(Y1 - Y0))

# Observed difference among uncensored (biased under MNAR)
obs = with(data, mean(mu.a1.u0 - mu.a0.u0))

################################################################################
### Examine censoring and outcomes
################################################################################

# Overall censoring rate
with(data, mean(C))

# Rate of informative censoring (U2)
with(data, mean(U2))

# Outcome probability among informationally censored
with(filter(data, U2 == 1), mean(mu.a1.u1))

# Outcome rate among uncensored
with(filter(data, U2 == 0), mean(Y))

# Outcome risk ratio: informative vs non-informative censoring
with(data, mean(mu.a1.u1/mu.a1.u0))

################################################################################
### Interrupted time series example (conceptual illustration)
################################################################################
# This section creates a simple visualization to illustrate interrupted time
# series with informative censoring

# Time grid and intervention points
time = seq(1, 10, 0.01)
int1 = 5  # First intervention time
int2 = 8  # Second intervention time

# Outcome trajectories: piecewise linear with slope changes at interventions
y00.0 = time * (time < int1) + (2*time - int1) * (time >= int1)  # Control trajectory
y00.1 = 2*time * (time < int1) + (3*time - int1) * (time >= int1)  # Treatment trajectory

# Create visualization of interrupted time series with censoring
tibble(
  time = time,
  y00.0 = y00.0,
  y00.1 = y00.1
) %>%
  gather(key, value, -time) %>%
  # Label treatment groups
  mutate(Treatment = case_when(
    grepl("0$", key) ~ "Control",
    grepl("1$", key) ~ "Treatment",
  )) %>%
  # Distinguish observed vs imputed periods
  mutate(Observed = case_when(
    time < int1 ~ "dashed",                                      # Pre-intervention: all observed
    time >= int2 ~ "round",                                      # Late post: all imputed
    time >= int1 & time < int2 & grepl("0$", key) ~ "dashed",  # Control: still observed
    time >= int1 & time < int2 & grepl("1$", key) ~  "round"   # Treatment: censored
  )) %>%
  # Add noise to simulate observed data points
  mutate(Points = round(value) + rnorm(nrow(.), 0, 4)) %>%
  mutate(tf = round(time, 0)) %>%
  group_by(tf, Treatment) %>%
  mutate(Points = mean(Points)) %>%
  # Only show points for observed periods
  mutate(Points = if_else(Observed == "dashed", Points, NA_real_)) %>%
  # Create plot
  ggplot(aes(x = time, y = value, color = Treatment)) +
  geom_line(aes(linetype = Observed)) +
  theme_minimal() +
  geom_point(aes(x = tf, y = Points)) +
  geom_vline(xintercept = 5, lty = "dashed", color = "firebrick") +  # Intervention 1
  geom_vline(xintercept = 8, lty = "dashed", color = "firebrick") +  # Intervention 2
  xlab("Time") +
  ylab("Outcome value") +
  scale_linetype_manual(labels = c("Observed", "Imputed"),
                        values = c(1, 8))

################################################################################
### Multiple intervention time series with random effects
################################################################################
# Demonstrates regression analysis with multiple interventions and time-varying
# treatment effects

# Intervention timepoints
T1 = 50  # First intervention
T2 = 75  # Second intervention

# Random treatment effects that kick in after interventions
fx1 = c(rep(0, 49), runif(51, 0, 1))  # Effect of treatment 1 (starts at T1)
fx2 = c(rep(0, 74), runif(26, 0, 1))  # Common effect (starts at T2)
fx3 = c(rep(0, 49), runif(51, 0, 1))  # Effect of treatment 2 (starts at T1)

# Generate data with piecewise linear trends and random effects
data = tibble(
  time = c(1:100),
  # Baseline outcome trajectories (before adding random effects)
  Y00.0 = time * (time < T1) + (2*(time - T1) + T1) * (time >= T1),
  Y00.1 = (1 + 2*time) * (time < T1) + (3*(time - T1) + 1 + 2*T1) * (time >= T1),
  Y00.2 = (2 + 3*time) * (time < T1) + (4*(time - T1) + 2 + 3*T1) * (time >= T1)
) %>%
  mutate(
    # Add random effects to observed outcomes
    Y.0 = Y00.0 + fx2,
    Y.1 = Y00.1 + fx1 + fx2,
    Y.2 = Y00.2 + fx3 + fx2
  )

# Reshape to long format for regression
ldata = data %>%
  select(time, Y.0, Y.1, Y.2) %>%
  gather(treatment, Y, -time) %>%
  mutate(treatment = as.numeric(gsub("Y\\.", "", treatment)),
         # Treatment indicators
         treat1 = if_else(treatment == 1, 1, 0),
         treat2 = if_else(treatment == 2, 1, 0),
         # Post-intervention indicators
         post1 = if_else(time >= T1, 1, 0),
         post2 = if_else(time >= T2, 1, 0),
         # Time since intervention (0 before intervention)
         time.1 = post1*time,
         time.2 = post2*time,
         # Factor versions for time-varying effects
         ftime.2 = factor(time.2),
         ftime.1 = factor(time.1))

# Fit regression model with treatment-specific time-varying effects
mod = lm(Y ~ time + treat1 + treat2 + time:treat1 + time:treat2 +
           post1 + time.1 + ftime.2 + treat1:ftime.1 + treat2:ftime.1, ldata)

# Extract and compare estimated treatment effects to true random effects
# Treatment 1 effects (should match fx1)
broom::tidy(mod) %>%
  filter(grepl("treat1\\:ftime\\.1", term)) %>%
  mutate_at('estimate', ~round(., 4)) %>%
  mutate(fx1 = fx1[50:length(fx1)]) %>%
  print.data.frame()

# Treatment 2 effects (should match fx3)
broom::tidy(mod) %>%
  filter(grepl("treat2\\:ftime\\.1", term)) %>%
  mutate_at('estimate', ~round(., 4)) %>%
  mutate(fx3 = fx3[50:length(fx3)]) %>%
  print.data.frame()

# Common time effects after T2 (should match fx2)
broom::tidy(mod) %>%
  filter(grepl("ftime\\.2", term)) %>%
  mutate_at('estimate', ~round(., 4)) %>%
  mutate(fx2 = fx2[75:length(fx2)]) %>%
  print.data.frame()


