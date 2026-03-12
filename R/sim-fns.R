################################################################################
########################### SIMULATION FUNCTIONS ###############################
################################################################################
#
# This file contains functions for generating simulated data with censoring
# mechanisms and calculating causal estimands. These functions support
# sensitivity analysis for missing data under the Missing at Random (MAR)
# assumption.
#
################################################################################

# Helper functions for transformations
expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log(x / (1 - x))

#' GenerateData: simulate covariate, treatment, non-informative, and informative
#' censoring given sample size and known sensitivity parameters
#'
#' This function generates a synthetic population with:
#' - A continuous covariate X
#' - Binary treatment A with propensity based on X
#' - Two types of censoring indicators: U1 (non-informative) and U2 (informative)
#' - Censoring probabilities that depend on treatment and covariates
#'
#' @param population_size Sample size for the generated population
#' @param delta1 Sensitivity parameter relating pi*_1(X) to pi_1(X) under treatment
#'   This controls the strength of informative vs non-informative censoring when A=1
#' @param delta0 Sensitivity parameter relating pi*_0(X) to pi_0(X) under control
#'   This controls the strength of informative vs non-informative censoring when A=0
#'
#' @return dataframe containing data and data generating parameters for:
#'   - X: covariate (uniform on [-3, 3])
#'   - A: binary treatment assignment
#'   - ps: propensity score P(A=1|X)
#'   - U1.1, U2.1: censoring indicators under treatment
#'   - U1.0, U2.0: censoring indicators under control
#'   - C.1, C.0: overall censoring indicators for each treatment
#'   - gamma.c1, gamma.c0: censoring probabilities by treatment
#'   - gamma.u1.1, gamma.u2.1: component censoring probabilities under treatment
#'   - gamma.u1.0, gamma.u2.0: component censoring probabilities under control
GenerateData <- function(population_size,
                         delta1,
                         delta0) {
  N = population_size
  data = tibble(
    # Generate covariate uniformly distributed from -3 to 3
    X = runif(N, -3, 3),
    # Propensity score: probability of treatment given X
    ps = expit(X),
    # Binary treatment assignment based on propensity score
    A = rbinom(N, 1, ps),

    # Censoring probability under treatment (A=1): piecewise function of X
    # Different functional forms for X<0, 0<=X<1, and X>=1
    gamma.c1 = expit(X) * (X < 0) + expit(0.8*X) * (X >= 0 & X < 1) + expit(0.2 + 0.6*X)*(X >= 1),
    # Informative censoring probability: fraction delta1 of total censoring
    gamma.u2.1 = delta1 * gamma.c1,
    # Non-informative censoring probability: remaining fraction
    gamma.u1.1 = gamma.c1 * (1 - delta1),

    # Censoring probability under control (A=0): piecewise function of X
    gamma.c0 = expit(0.6*X) * (X < 0) + expit(0.5*X) * (X >= 0 & X < 1) + expit(0.1 + 0.4*X)*(X >= 1),
    # Informative censoring probability under control
    gamma.u2.0 = delta0 * gamma.c0,
    # Non-informative censoring probability under control
    gamma.u1.0 = gamma.c0 * (1 - delta0),

    # Generate censoring indicators under treatment (potential outcomes)
    U1.1 = rbinom(N, 1, gamma.u1.1),
    # U2 only occurs if U1 doesn't occur (mutually exclusive)
    U2.1 = ifelse(U1.1 == 1, 0, rbinom(N, 1, gamma.u2.1 / (1 - gamma.u1.1))),
    # Overall censoring under treatment
    C.1  = as.integer(U1.1 == 1 | U2.1 == 1),

    # Generate censoring indicators under control (potential outcomes)
    U1.0 = rbinom(N, 1, gamma.u1.0),
    U2.0 = ifelse(U1.0 == 1, 0, rbinom(N, 1, gamma.u2.0 / (1 - gamma.u1.0))),
    # Overall censoring under control
    C.0  = as.integer(U1.0 == 1 | U2.0 == 1),

    # Observed censoring indicators (based on actual treatment received)
    U1 = A*U1.1 + (1-A)*U1.0,
    U2 = A*U2.1 + (1-A)*U2.0,
    gamma = A * gamma.c1 + (1 - A) * gamma.c0,
    C  = ifelse(U1 == 1 | U2 == 1, 1, 0)
  )
  data
}

#' GenerateCEF: generate conditional expectation function for outcome given
#' informative censoring indicator, i.e. P(Y = 1 | X, A = a, U2 = u)
#'
#' This function generates outcome probabilities conditional on covariates,
#' treatment, and censoring type. The relationship between outcomes under
#' informative vs non-informative censoring is controlled by the tau parameter.
#'
#' @param data Input data from GenerateData
#' @param a_val Value of treatment assignment indicator (0 or 1)
#' @param tau Sensitivity parameter relating mu_a(X, 1) to mu_a(X, 0)
#'   If tau >= 1: outcome risk under informative censoring is tau times higher
#'   If tau < 1: outcome risk under informative censoring is tau times lower
#' @param coef Outcome coefficient on covariates X (default 1)
#'   Controls the strength of the X-outcome relationship
#'
#' @return Dataframe with added columns:
#'   - mu.a{a_val}.u0: P(Y=1 | X, A=a, U2=0) - outcome prob under non-informative censoring
#'   - mu.a{a_val}.u1: P(Y=1 | X, A=a, U2=1) - outcome prob under informative censoring
GenerateCEF <- function(data,
                        a_val,
                        tau,
                        coef = 1) {
  # Create variable names for the two outcome probabilities
  var_name0 = paste0("mu.a", a_val, ".u0")  # Non-informative censoring
  var_name1 = paste0("mu.a", a_val, ".u1")  # Informative censoring

  # If tau >= 1, informative censoring increases outcome risk
  if (tau >= 1) {
    data[[var_name1]] = expit(coef * data$X)
    data[[var_name0]] = data[[var_name1]] / tau
  }

  # If tau < 1, informative censoring decreases outcome risk
  if (tau < 1) {
    data[[var_name0]] = expit(coef * data$X)
    data[[var_name1]] = data[[var_name0]] * tau
  }
  data
}

#' GenerateOutcomes: Generate potential and observed outcomes
#'
#' This function generates all potential outcomes Y(a, u) for each combination
#' of treatment (a=0,1) and censoring type (u=0,1), then creates the observed
#' outcomes based on actual treatment assignment and censoring status.
#'
#' @param data Output from GenerateCEF containing outcome probabilities
#'
#' @return Dataframe with additional outcome columns:
#'   - Y11: potential outcome under treatment and informative censoring
#'   - Y10: potential outcome under treatment and non-informative censoring
#'   - Y01: potential outcome under control and informative censoring
#'   - Y00: potential outcome under control and non-informative censoring
#'   - Y1, Y0: potential outcomes under treatment/control (marginalized over censoring)
#'   - Y: observed outcome
#'   - mu1, mu0: average potential outcomes under treatment/control
#'   - Z1, Z0, Z: composite outcome (censoring OR outcome event)
GenerateOutcomes <- function(data) {
  N = nrow(data)
  data %>%
    mutate(
      # Generate all four potential outcomes Y(a, u2)
      Y11 = rbinom(N, 1, mu.a1.u1),  # Treatment + informative censoring
      Y10 = rbinom(N, 1, mu.a1.u0),  # Treatment + non-informative censoring
      Y01 = rbinom(N, 1, mu.a0.u1),  # Control + informative censoring
      Y00 = rbinom(N, 1, mu.a0.u0),  # Control + non-informative censoring

      # Potential outcomes marginalized over censoring type
      Y1 = U2.1 * Y11 + (1 - U2.1) * Y10,  # Under treatment
      Y0 = U2.0 * Y01 + (1 - U2.0) * Y00,  # Under control

      # Expected potential outcomes (probabilities)
      mu1 = gamma.u2.1 * mu.a1.u1 + (1 - gamma.u2.1) * mu.a1.u0,
      mu0 = gamma.u2.0 * mu.a0.u1 + (1 - gamma.u2.0) * mu.a0.u0,

      # Observed outcome (based on actual treatment received)
      Y  = A * Y1 + (1 - A) * Y0
    ) %>%
    mutate(
      # Composite outcome Z: equals 1 if either censored (U2=1) OR outcome occurred (Y=1)
      Z1 = if_else(U2.1 == 1 | Y1 == 1, 1, 0),  # Under treatment
      Z0 = if_else(U2.0 == 1 | Y0 == 1, 1, 0),  # Under control
      Z = A * Z1 + (1 - A) * Z0                  # Observed
    )
}

#' CreatePopulationData: Wrapper function to create a complete population dataset
#'
#' This is a convenience function that chains together GenerateData, GenerateCEF,
#' and GenerateOutcomes to create a complete synthetic population in one step.
#'
#' @param population_size Size of the population to generate
#' @param delta1 Sensitivity parameter for censoring under treatment
#'   Controls fraction of informative vs non-informative censoring when A=1
#' @param delta0 Sensitivity parameter for censoring under control
#'   Controls fraction of informative vs non-informative censoring when A=0
#' @param tau1 Sensitivity parameter relating outcome risk under treatment
#'   Ratio of mu_1(X, 1)/mu_1(X, 0) when tau1>=1, or inverse when tau1<1
#' @param tau0 Sensitivity parameter relating outcome risk under control
#'   Ratio of mu_0(X, 1)/mu_0(X, 0) when tau0>=1, or inverse when tau0<1
#'
#' @return Complete population dataframe with covariates, treatment, censoring,
#'   and outcome variables
CreatePopulationData <- function(population_size,
                                 delta1,
                                 delta0,
                                 tau1,
                                 tau0) {
  # Generate base data with censoring structure
  GenerateData(population_size, delta1 = delta1, delta0 = delta0) %>%
    # Add outcome model under treatment (coefficient = 1)
    GenerateCEF(a_val = 1, tau = tau1, coef = 1) %>%
    # Add outcome model under control (coefficient = 0.5)
    GenerateCEF(a_val = 0, tau = tau0, coef = 0.5) %>%
    # Generate realized outcomes
    GenerateOutcomes()
}

#' CalculateTruth: Calculate true causal estimands from population data
#'
#' Computes the ground truth values for three causal estimands:
#' 1. ATE (Average Treatment Effect): E[Y(1) - Y(0)]
#' 2. Z-ATE: E[Z(1) - Z(0)] where Z is composite outcome (censoring or outcome)
#' 3. SDE (Separable Direct Effect): Treatment effect excluding informative censoring pathway
#'
#' @param data Output from CreatePopulationData containing potential outcomes
#'
#' @return Tibble with two columns:
#'   - estimand: name of the causal parameter ("ate", "z-ate", or "sde")
#'   - truth: true population value of the estimand
CalculateTruth <- function(data) {
  # Average Treatment Effect on primary outcome Y
  ate  = with(data, mean(Y1 - Y0))

  # Average Treatment Effect on composite outcome Z (censoring or outcome)
  zate = with(data, mean(Z1 - Z0))

  # Separable Direct Effect: treatment effect among those who would not be
  # informationally censored under control
  sde  = with(data, mean((1 - gamma.u2.0) * (mu.a1.u0 - mu.a0.u0)))

  tibble(
    estimand = c("ate", "z-ate", "sde"),
    truth = c(ate, zate, sde)
  )
}

#' AddEstimationError: Add noise to nuisance functions to simulate estimation error
#'
#' This function mimics the estimation error that would occur in practice when
#' estimating nuisance functions (outcome models, censoring models, propensity scores)
#' from finite samples. The amount of error is controlled by convergence rate parameters.
#'
#' @param data Output from CreatePopulationData with true nuisance functions
#' @param mu.rate Convergence rate for outcome models mu_a(X, 0)
#'   Error magnitude scales as n^(-mu.rate), where n is sample size
#' @param gamma.rate Convergence rate for censoring probability models pi_a(X)
#'   Error magnitude scales as n^(-gamma.rate)
#' @param pi.rate Convergence rate for propensity score e(X)
#'   Error magnitude scales as n^(-pi.rate)
#'
#' @return Dataset with Gaussian noise added to nuisance functions on logit scale
#'   This preserves the (0,1) range of probabilities after back-transforming via expit
AddEstimationError <- function(data,
                     mu.rate,
                     gamma.rate,
                     pi.rate) {
  n = nrow(data)

  data = data %>%
    # Add noise to outcome models: mu.a0.u0, mu.a1.u0, etc.
    mutate_at(vars(matches("mu\\.")), ~expit(logit(.) + rnorm(n, n^(-mu.rate), n^(-mu.rate)))) %>%
    # Add noise to censoring probability models: gamma.c0, gamma.c1, etc.
    mutate_at(vars(matches("gamma\\.")), ~expit(logit(.) + rnorm(n, n^(-gamma.rate), n^(-gamma.rate)))) %>%
    # Add noise to propensity score
    mutate_at(vars(matches("ps")), ~expit(logit(.) + rnorm(n, n^(-pi.rate), n^(-pi.rate))))

  data
}

#' IterateSimParams: Calculate bounds for multiple parameter specifications
#'
#' This function takes a grid of sensitivity parameter specifications and computes
#' the corresponding bounds or point estimates for each specification. This is useful
#' for evaluating performance across different assumed sensitivity parameter values.
#'
#' @param param_spec Dataframe of simulation parameter specifications with columns:
#'   - estimand: target estimand ("ate", "z-ate", or "sde")
#'   - monotonicity: monotonicity assumption ("none", "positive", or "negative")
#'   - outcome_assumption: outcome model assumption ("none", "equality", "inequality")
#'   - delta1, delta0: censoring sensitivity parameters
#'   - tau1, tau0: outcome sensitivity parameters
#'   - bound: logical indicating whether to compute bounds (TRUE) or point estimate (FALSE)
#' @param eif_data Data containing nuisance components (from ConstructPI)
#'
#' @return Dataframe of parameter specifications with added column 'estimands'
#'   containing the computed bounds/estimates for each specification
IterateSimParams <- function(param_spec, eif_data) {
  param_spec %>%
    # For each row of specifications, call CensoredBounds with those parameters
    mutate(estimands = pmap(list(estimand = estimand,
                                 monotonicity = monotonicity,
                                 outcome_assumption = outcome_assumption,
                                 delta1 = delta1,
                                 delta0 = delta0,
                                 tau1 = tau1,
                                 tau0 = tau0,
                                 bound = bound), ~CensoredBounds(eif_data = eif_data,
                                                                 estimator = "pi",
                                                                 ...)))
}

#' EasySimulation: Monte Carlo simulation with known nuisance functions + noise
#'
#' This function runs a simplified simulation study where nuisance functions are
#' known from the population data generating process, and estimation error is added
#' via AddNoise. This is faster than HardSimulation and allows for testing the
#' influence function-based estimators without the complexity of machine learning.
#'
#' @param popdata Population data from CreatePopulationData
#' @param sample_size Size of each simulated sample
#' @param number_sims Number of Monte Carlo replications
#' @param mu_rate Convergence rate for outcome model estimation error (default 0.25)
#' @param pi_rate Convergence rate for propensity score estimation error (default 0.25)
#' @param gamma_rate Convergence rate for censoring model estimation error (default 0.25)
#'
#' @return List of length number_sims, where each element is a tibble containing:
#'   - component: name of the bound component (psi_tilde, Omega_1, ..., Omega_5)
#'   - ests: point estimate
#'   - varests: variance estimate
#'   - lci, uci: 95% confidence interval bounds
#'   - truth: true population value
#'   - covered: indicator for whether CI covers truth
EasySimulation <- function(popdata, sample_size, number_sims,
                           mu_rate = 0.25, pi_rate = 0.25, gamma_rate = 0.25) {

  # Construct plug-in estimates of bound components from population
  popdata = ConstructPI(popdata)

  # Calculate true values of all bound components
  truth = c(with(popdata, mean(phi.ate.pi)),           # psi_tilde: biased ATE
            with(popdata, mean(phi.1.pi)),             # Omega_1: censoring under trt
            with(popdata, mean(phi.0.pi)),             # Omega_2: censoring under control
            with(popdata, mean(phi.11.pi)),            # Omega_3: outcome*censoring (trt)
            with(popdata, mean(phi.00.pi)),            # Omega_4: outcome*censoring (ctrl)
            with(popdata, mean(phi.01.pos.pi - phi.00.pos.pi)))  # Omega_5: SDE component

  res = list()

  # Run Monte Carlo simulations
  for (i in 1:number_sims) {

    # Draw sample, add estimation error, construct influence functions
    samp = sample_n(popdata, sample_size) %>%
      AddNoise(mu_rate = mu_rate,
               pi_rate = pi_rate,
               gamma_rate = gamma_rate) %>%
      ConstructEIF()

    # Compute point estimates for all components
    ests = c(with(samp, mean(phi.ate)),
             with(samp, mean(phi.1)),
             with(samp, mean(phi.0)),
             with(samp, mean(phi.11)),
             with(samp, mean(phi.00)),
             with(samp, mean(phi.01.pos - phi.00.pos)))

    # Compute variance estimates (sample variance / n)
    varests = c(with(samp, var(phi.ate)),
                with(samp, var(phi.1)),
                with(samp, var(phi.0)),
                with(samp, var(phi.11)),
                with(samp, var(phi.00)),
                with(samp, var(phi.01.pos - phi.00.pos))) / sample_size

    # Store results for this simulation
    res[[i]] = tibble(
      component = c("$\\tilde{\\Psi}$", "$\\Omega_1$",
                    "$\\Omega_2$", "$\\Omega_3$", "$\\Omega_4$",
                    "$\\Omega_5$"),
      ests = ests,
      varests = varests,
      lci = ests - 1.96 * sqrt(varests),  # 95% CI lower bound
      uci = ests + 1.96 * sqrt(varests),  # 95% CI upper bound
      truth = truth,
      covered = ifelse(lci < truth & uci > truth, 1, 0)  # Coverage indicator
    )
  }
  res
}

#' HardSimulation: Monte Carlo simulation with SuperLearner-estimated nuisance functions
#'
#' This function runs a realistic simulation study where nuisance functions
#' (outcome models, censoring models, propensity scores) are estimated from
#' observed data using SuperLearner ensemble machine learning. This provides
#' a more realistic assessment of finite-sample performance.
#'
#' @param popdata Population data from CreatePopulationData
#' @param sample_size Size of each simulated sample
#' @param number_sims Number of Monte Carlo replications
#' @param sl_lib Character vector of SuperLearner algorithms to include in ensemble
#'   (e.g., c("SL.glm", "SL.ranger", "SL.earth"))
#'
#' @return List of length number_sims, where each element is a tibble containing:
#'   - component: name of the bound component
#'   - ests: point estimate
#'   - varests: variance estimate
#'   - lci, uci: 95% confidence interval bounds
#'   - truth: true population value
#'   - covered: indicator for whether CI covers truth
HardSimulation <- function(popdata, sample_size, number_sims, sl_lib) {

  # Construct plug-in estimates from population for truth calculations
  popdata = ConstructPI(popdata)

  # Calculate true values of all bound components
  truth = c(with(popdata, mean(phi.ate.pi)),
            with(popdata, mean(phi.1.pi)),
            with(popdata, mean(phi.0.pi)),
            with(popdata, mean(phi.11.pi)),
            with(popdata, mean(phi.00.pi)),
            with(popdata, mean(phi.01.pos.pi - phi.00.pos.pi)))

  res = list()

  # Run Monte Carlo simulations
  for (i in 1:number_sims) {
    # Draw sample with only observed variables (X, A, C, Y)
    samp = select(popdata, X, A, C, Y) %>%
      sample_n(sample_size) %>%
      mutate(id = 1:nrow(.),
             cluster_id = id)

    # Estimate nuisance functions using SuperLearner with cross-fitting
    samp = EstimateNuisanceFunctions(X_dat = tibble(id = samp$id,
                                                    cluster_id = samp$id,
                                                    X = samp$X),
                                     outcome_tx_dat = samp %>% select(-X),
                                     outer_folds = 2,    # 2-fold cross-fitting
                                     sl_folds = 2,       # 2-fold CV for SuperLearner
                                     sl_lib = sl_lib) %>%
      ConstructPI() %>%
      ConstructEIF()

    # Compute point estimates for all components
    ests = c(with(samp, mean(phi.ate)),
             with(samp, mean(phi.1)),
             with(samp, mean(phi.0)),
             with(samp, mean(phi.11)),
             with(samp, mean(phi.00)),
             with(samp, mean(phi.01.pos - phi.00.pos)))

    # Compute variance estimates
    varests = c(with(samp, var(phi.ate)),
             with(samp, var(phi.1)),
             with(samp, var(phi.0)),
             with(samp, var(phi.11)),
             with(samp, var(phi.00)),
             with(samp, var(phi.01.pos - phi.00.pos))) / sample_size

    # Store results for this simulation
    res[[i]] = tibble(
      component = c("$\\tilde{\\Psi}$", "$\\Omega_1$", "$\\Omega_2$",
                    "$\\Omega_3$", "$\\Omega_4$", "$\\Omega_5$"),
      ests = ests,
      varests = varests,
      lci = ests - 1.96 * sqrt(varests),
      uci = ests + 1.96 * sqrt(varests),
      truth = truth,
      covered = ifelse(lci < truth & uci > truth, 1, 0)
    )
  }
  res
}

#' Simulation: Generic simulation function for evaluating bounds
#'
#' This is a general-purpose simulation function that evaluates the performance
#' of bound estimates under specific assumptions. It draws repeated samples from
#' population data and computes bounds/estimates for each sample.
#'
#' @param population_data Population data with influence function components
#' @param number_simulations Number of Monte Carlo replications
#' @param sample_size Size of each simulated sample
#' @param estimand Target estimand ("ate", "z-ate", or "sde")
#' @param outcome_assumption Outcome model assumption ("none", "equality", "inequality")
#' @param monotonicity Monotonicity assumption ("none", "positive", "negative")
#' @param truth True value of the target parameter (for coverage evaluation)
#' @param delta1, delta0 Censoring sensitivity parameters
#' @param tau1, tau0 Outcome sensitivity parameters
#' @param bound Logical indicating whether to compute bounds (TRUE) or point estimate (FALSE)
#'
#' @return Dataframe with one row per simulation containing estimates and confidence intervals
Simulation <- function(population_data,
                       number_simulations,
                       sample_size,
                       estimand,
                       outcome_assumption,
                       monotonicity, truth,
                       delta1, delta0, tau1, tau0,
                       bound) {

  res = list()

  # Run repeated simulations
  for (i in 1:number_simulations) {
    # Draw sample from population
    sample <- sample_n(population_data, sample_size)

    # Compute bounds/estimates with specified assumptions
    res[[i]] <- CensoredBounds(estimand = estimand,
                               monotonicity = monotonicity,
                               delta1 = delta1, delta0 = delta0,
                               outcome_assumption = outcome_assumption,
                               tau1 = tau1, tau0 = tau0,
                               bound = bound,
                               simulation = TRUE,
                               eif_data = sample,
                               truth = truth)
  }
  # Combine all simulation results into single dataframe
  invoke(rbind, res)
}

#' AddNoise: Add estimation error with directional bias to nuisance functions
#'
#' Similar to AddEstimationError but adds noise with non-zero mean to create
#' directional bias in the nuisance function estimates. This can be used to
#' assess robustness to model misspecification.
#'
#' @param data Dataset with true nuisance functions
#' @param mu_rate Convergence rate for outcome model noise
#' @param pi_rate Convergence rate for propensity score noise
#' @param gamma_rate Convergence rate for censoring model noise
#'
#' @return Dataset with biased noise added to nuisance functions
#'   - Positive bias for mu.a0, gamma.c1, and ps
#'   - Negative bias for mu.a1 and gamma.c0
AddNoise <- function(data, mu_rate, pi_rate, gamma_rate) {
  n = nrow(data)
  data %>%
    # Add positively biased noise to control outcome models
    mutate_at(vars(matches("mu\\.a0")), ~expit(logit(.) + rnorm(n, 1/n^mu_rate, 1/n^mu_rate))) %>%
    # Add negatively biased noise to treatment outcome models
    mutate_at(vars(matches("mu\\.a1")), ~expit(logit(.) + rnorm(n, -2/n^mu_rate, 1/n^mu_rate))) %>%
    # Add positively biased noise to propensity score
    mutate_at(vars(matches("ps")), ~expit(logit(.) + rnorm(n, 2/n^pi_rate, 1/n^pi_rate))) %>%
    # Add positively biased noise to censoring model under treatment
    mutate_at(vars(matches("gamma\\.c1")), ~expit(logit(.) + rnorm(n, 1/n^gamma_rate, 1/n^gamma_rate))) %>%
    # Add negatively biased noise to censoring model under control
    mutate_at(vars(matches("gamma\\.c0")), ~expit(logit(.) + rnorm(n, -2/n^gamma_rate, 1/n^gamma_rate)))
}


