################################################################################
######################## RUNNING EXAMPLE: BOUNDS ###############################
################################################################################
#
# This script provides a detailed worked example computing various bounds
# for a simple two-covariate scenario. Demonstrates:
# - Computing true ATE and biased ATE under MAR
# - General (assumption-free) bounds
# - Bounds under monotonicity assumptions
# - Bounds under bounded outcome assumptions
# - Tipping point analysis
# - Alternative estimands (Z-ATE, SDE)
#
################################################################################

################################################################################
### Data specification
################################################################################
# Generate a large population with binary covariate X

px = 0.7  # P(X=1)
X <- rbinom(5e7, 1, px)

# Outcome probabilities conditional on treatment, covariate, and censoring type
mu0 = 0.1 + 0.05*X      # mu_0(X, 0): P(Y=1 | A=0, X, U2=0)
mu1 = 2 * mu0           # mu_1(X, 0): P(Y=1 | A=1, X, U2=0) - twice as high
mu0s = 2 * mu0          # mu_0(X, 1): P(Y=1 | A=0, X, U2=1) - informatively censored, higher risk
mu1s = 2 * mu1          # mu_1(X, 1): P(Y=1 | A=1, X, U2=1)

# Censoring probabilities (informative censoring only)
pi0s = 0.07 + 0.05*X    # pi*_0(X): P(U2=1 | A=0, X)
pi1s = 0.14 + 0.05*X    # pi*_1(X): P(U2=1 | A=1, X) - higher under treatment

# Total censoring probabilities (observed)
pi0  = (3/2) * pi0s     # pi_0(X): P(C=1 | A=0, X) = pi*_0 + pi^u_0
pi1  = (3/2) * pi1s     # pi_1(X): P(C=1 | A=1, X)

overall_censoring = mean(0.5*(pi0 + pi1)) #P(C = 1)

################################################################################
### True parameters
################################################################################
# True ATE (accounting for informative censoring)
psi = mean(mu1 + pi1s*(mu1s - mu1) - mu0 - pi0s*(mu0s - mu0))

# Biased ATE (observable under MAR assumption)
psitilde = mean(mu1 - mu0)

# Proportion of bias
pctr = 1-psitilde/psi

c(psi, psitilde, pctr)

################################################################################
### General (assumption-free) bounds
################################################################################
# Sharp bounds with no sensitivity parameter assumptions (delta=1)
ell0 = psitilde - mean(pi1*mu1 + pi0*(1-mu0))
u0 = psitilde + mean(pi1*(1-mu1) + pi0*mu0)
c(ell0, u0)

################################################################################
### Bounds under positive monotonicity
################################################################################
# Assume mu_a(X,1) >= mu_a(X,0) for both a (informatively censored have higher risk)

# Monotonicity with delta=1 (all censoring is informative)
ell0_pos = psitilde - mean(pi0*(1-mu0))
u0_pos = psitilde + mean(pi1*(1-mu1))
c(ell0_pos, u0_pos)

# Monotonicity with delta=0.8 (80% of censoring is informative)
ell0_pos1 = psitilde - mean(0.8*pi0*(1-mu0))
u0_pos1 = psitilde + mean(0.8*pi1*(1-mu1))
c(ell0_pos1, u0_pos1)

################################################################################
### Bounds with bounded counterfactual outcome risk
################################################################################
# Assume (1/tau) <= mu_a(X,1)/mu_a(X,0) <= tau
# Here tau=3 means risk can be up to 3x higher or 3x lower

# Bounded risk with delta=1
ell0_tilde = psitilde + mean((1/3 - 1) * pi1*mu1 - (3-1)*pi0*mu0)
u0_tilde = psitilde + mean((3 - 1) * pi1*mu1 - (1/3-1)*pi0*mu0)
c(ell0_tilde, u0_tilde)

# Bounded risk with delta=0.8
ell0_tilde1 = psitilde + mean(0.8*(1/3 - 1) * pi1*mu1 - 0.8*(3-1)*pi0*mu0)
u0_tilde1 = psitilde + mean(0.8*(3 - 1) * pi1*mu1 - 0.8*(1/3-1)*pi0*mu0)
c(ell0_tilde1, u0_tilde1)

################################################################################
### Tipping point analysis
################################################################################
# Find sensitivity parameter values where conclusions would change

# Tipping point for positive monotonicity: value of delta where lower bound = 0
tipping_pos = psitilde / mean(pi0 * (1 - mu0))
tipping_pos

# Components of the bias
c(psitilde, mean(pi0*mu0), mean(pi1*mu1))

# Tipping point for equality assumption: value of tau where estimate = 0
tipping_eq = 1 + psitilde / mean(pi0*mu0)
tipping_eq

# Tipping point surface (delta vs tau)
delta1 = 1; tau = 10
slope = mean(mu1*pi1)/mean(mu0*pi0)
intercept = (psitilde / (tau - 1)) / mean(mu0*pi0)
c(intercept, slope)

################################################################################
### Alternative estimands
################################################################################

# Z-ATE: Treatment effect on composite outcome Z = I(U2=1 or Y=1)
psi1 = psitilde + mean(pi1s*(1-mu1) - pi0s*(1-mu0))

# SDE: Separable Direct Effect (excluding informative censoring pathway)
psi2 = psitilde - mean(pi0s*(mu1-mu0))

c(psi1, psi2)

# Bounds for Z-ATE under positive monotonicity
ell1 = psitilde - mean(pi0*(1-mu0))
u1   = psitilde + mean(pi1*(1-mu1))
c(ell1, u1)
psi1

# Bounds for SDE (using sign of treatment effect)
u2 = psitilde - mean(pi0*(mu1-mu0)*(mu1-mu0 <= 0))   # Upper bound
ell2 = psitilde - mean(pi0*(mu1-mu0)*(mu1-mu0 > 0))  # Lower bound
c(ell2, u2)

################################################################################
### Sensitivity plot example: Tipping point surface
################################################################################
# This section creates a contour plot showing combinations of (delta_0, delta_1)
# that yield ATE = 0 for different values of tau. These curves represent
# "tipping points" where conclusions would change from significant to null.

# Compute expected values of mu*pi (weighted outcome risk × censoring probability)
# These represent the contribution to bias from informative censoring
mugrave_0 = mean(mu0*pi0)  # E[mu_0(X) * pi_0(X)]: Control arm contribution
mugrave_1 = mean(mu1*pi1)  # E[mu_1(X) * pi_1(X)]: Treatment arm contribution

# Create grid of sensitivity parameter combinations where ATE = 0
# The relationship comes from the sensitivity analysis bound formula:
#   psi = psitilde + (tau-1)*(delta_1*mugrave_1 - delta_0*mugrave_0)
# Setting psi = 0 and solving for delta_0 gives:
#   delta_0 = delta_1*mugrave_1/mugrave_0 + psitilde/(tau-1)/mugrave_0
pd <- expand_grid(
  tau = seq(8, 15, by = 1),  # Risk ratio values (how much higher risk in censored)
  delta_1 = seq(0, mugrave_0/mugrave_1, length = 1000)  # Treatment informative fraction
) %>%
  # Compute corresponding delta_0 values that give ATE = 0
  mutate(delta_0 = delta_1*mugrave_1/mugrave_0 + psitilde/(tau-1)/mugrave_0) %>%
  # Keep only valid sensitivity parameter values (delta must be in [0,1])
  filter(delta_0 < 1)

# Create contour plot showing tipping point curves for each tau value
plot <- ggplot(pd, aes(x = delta_0, y = delta_1, group = tau, color = as.factor(tau))) +
  geom_line(linewidth = 1.5) +
  theme_rand() +
  labs(x = expression(delta[0]),  # Fraction informative censoring in control
       y = expression(delta[1]),  # Fraction informative censoring in treatment
       color = expression(tau)) + # Risk ratio for censored vs uncensored
  scale_color_manual(values = RandCatPal) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.text = element_text(size = 18),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18))

# Interpretation: Points below/left of a curve are consistent with ATE > 0,
# while points above/right suggest ATE < 0 for that tau value.
# Higher tau values (more extreme risk ratios) require larger delta values
# to overturn the positive finding.

ggsave('plots/sensitivity-plot.png', height = 5, width = 7)
