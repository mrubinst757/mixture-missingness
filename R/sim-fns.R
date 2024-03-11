#' GenerateData: simulate covariate, treatment, non-informative, and informative
#' censoring given sample size and known sensitivity parameters
#'
#' @param N population size
#' @param delta1 parameter relating \pi*_1(X) to \pi_1(X)
#' @param delta0 parameter relatin \pi_0*(X) to \pi_0(X)
#'
#' @return dataframe containing data and data generating parameters for
#' covariate X, treatment A, non-informative censoring indicator U1, and informative
#' censoring indicator U2.
GenerateData <- function(population_size,
                         delta1,
                         delta0) {
  N = population_size
  data = tibble(
    X = runif(N, -3, 3),
    ps = expit(X),
    A = rbinom(N, 1, ps),
    gamma.c1 = expit(X) * (X < 0) + expit(0.8*X) * (X >= 0 & X < 1) + expit(0.2 + 0.6*X)*(X >= 1),
    gamma.u2.1 = delta1 * gamma.c1,
    gamma.u1.1 = gamma.c1 * (1 - delta1) / (1 - delta1 * gamma.c1),
    gamma.c0 = expit(0.6*X) * (X < 0) + expit(0.5*X) * (X >= 0 & X < 1) + expit(0.1 + 0.4*X)*(X >= 1),
    gamma.u2.0 = delta0 * gamma.c0,
    gamma.u1.0 = gamma.c0 * (1 - delta0) / (1 - delta0 * gamma.c0),
    U1.1 = rbinom(N, 1, gamma.u1.1),
    U1.0 = rbinom(N, 1, gamma.u1.0),
    U2.1 = rbinom(N, 1, gamma.u2.1),
    U2.0 = rbinom(N, 1, gamma.u2.0),
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
#' @param data input data from GenerateData
#' @param a_val value of treatment assignment indicator
#' @param tau parameter relating \mu_a(X, 1) to \mu_a(X, 0)
#' @param coef outcome coefficient on covariates X
#'
#' @return dataframe containing \mu_a(X, 0) and \mu_a(X, 0)
GenerateCEF <- function(data,
                        a_val,
                        tau,
                        coef = 1) {
  var_name0 = paste0("mu.a", a_val, ".u0")
  var_name1 = paste0("mu.a", a_val, ".u1")

  if (tau >= 1) {
    data[[var_name1]] = expit(coef * data$X)
    data[[var_name0]] = data[[var_name1]] / tau
  }

  if (tau < 1) {
    data[[var_name0]] = expit(coef * data$X)
    data[[var_name1]] = data[[var_name0]] * tau
  }
  data
}

#' GenerateOutcomes
#'
#' @param data output from GenerateCEF
#'
#' @return generates potential outcomes for all individuals in the dataframe
GenerateOutcomes <- function(data) {
  N = nrow(data)
  data %>%
    mutate(
      Y11 = rbinom(N, 1, mu.a1.u1),
      Y10 = rbinom(N, 1, mu.a1.u0),
      Y01 = rbinom(N, 1, mu.a0.u1),
      Y00 = rbinom(N, 1, mu.a0.u0),
      Y1 = U2.1 * Y11 + (1 - U2.1) * Y10,
      Y0 = U2.0 * Y01 + (1 - U2.0) * Y00,
      mu1 = gamma.u2.1 * mu.a1.u1 + (1 - gamma.u2.1) * mu.a1.u0,
      mu0 = gamma.u2.0 * mu.a0.u1 + (1 - gamma.u2.0) * mu.a0.u0,
      Y  = A * Y1 + (1 - A) * Y0
    ) %>%
    mutate(
      Z1 = if_else(U2.1 == 1 | Y1 == 1, 1, 0),
      Z0 = if_else(U2.0 == 1 | Y0 == 1, 1, 0),
      Z = A * Z1 + (1 - A) * Z0
    )
}

#' CreatePopulationData
#'
#' @param population_size size of population
#' @param delta1 parameter relating \pi*_1(X) to \pi_1(X)
#' @param delta0 parameter relating \pi*_0(X) to \pi_0(X)
#' @param tau1 parameter relating \mu_1(X, 0) to \mu_1(X, 1)
#' @param tau0 parameter relating \mu_0(X, 0) to \mu_0(X, 1)
#'
#' @return population dataframe
CreatePopulationData <- function(population_size,
                                 delta1,
                                 delta0,
                                 tau1,
                                 tau0) {
  GenerateData(population_size, delta1 = delta1, delta0 = delta0) %>%
    GenerateCEF(a_val = 1, tau = tau1, coef = 1) %>%
    GenerateCEF(a_val = 0, tau = tau0, coef = 0.5) %>%
    GenerateOutcomes()
}

#' CalculateTruth: calculate true causal estimands of potential interest
#'
#' @param data output from CreatePopulationData
#'
#' @return dataframe containing true causal parameters
CalculateTruth <- function(data) {
  ate  = with(data, mean(Y1 - Y0))
  zate = with(data, mean(Z1 - Z0))
  sde  = with(data, mean((1 - gamma.u2.0) * (mu.a1.u0 - mu.a0.u0)))

  tibble(
    estimand = c("ate", "z-ate", "sde"),
    truth = c(ate, zate, sde)
  )
}

#' AddEstimationError: add noise at specified error rates to true data generating
#' functions to simulate estimation error
#'
#' @param data output from CreatePopulationData
#' @param mu.rate convergence rate on \mu_a(X, 0) (RMSE = O_p(n^-mu_rate))
#' @param gamma.rate convergence rate on \pi_a(X)
#' @param pi.rate convergence rate on e(X)
#'
#' @return dataset with error added to nuisance functions
AddEstimationError <- function(data,
                     mu.rate,
                     gamma.rate,
                     pi.rate) {
  n = nrow(data)

  data = data %>%
    mutate_at(vars(matches("mu\\.")), ~expit(logit(.) + rnorm(n, n^(-mu.rate), n^(-mu.rate)))) %>%
    mutate_at(vars(matches("gamma\\.")), ~expit(logit(.) + rnorm(n, n^(-gamma.rate), n^(-gamma.rate)))) %>%
    mutate_at(vars(matches("ps")), ~expit(logit(.) + rnorm(n, n^(-pi.rate), n^(-pi.rate))))

  data
}

#' IterateSimParams: add a column of true parameters to estimate given a dataframe
#' containing parameter specifications for simulations
#'
#' @param param_spec dataframe of simulation parameter specifications with columns
#' estimand, monotonicity, outcome assumption, delta1, delta0, tau1, tau0, and bound
#' @param eif_data data containing nuisance components
#'
#' @return dataframe of parameter specifications along with the true targeted
#' functionals
IterateSimParams <- function(param_spec, eif_data) {
  param_spec %>%
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

EasySimulation <- function(popdata, sample_size, number_sims,
                           mu_rate = 0.25, pi_rate = 0.25, gamma_rate = 0.25) {


  popdata = ConstructPI(popdata)
  truth = c(with(popdata, mean(phi.ate.pi)),
            with(popdata, mean(phi.1.pi - phi.0.pi)),
            with(popdata, mean(phi.11.pi - phi.00.pi)),
            with(popdata, mean(phi.01.pi - phi.00.pi)),
            with(popdata, mean(phi.01.pos.pi - phi.00.pos.pi)))

  res = list()

  for (i in 1:number_sims) {

    samp = sample_n(popdata, sample_size) %>%
      AddNoise(mu_rate = mu_rate,
               pi_rate = pi_rate,
               gamma_rate = gamma_rate) %>%
      ConstructEIF()

    ests = c(with(samp, mean(phi.ate)),
             with(samp, mean(phi.1 - phi.0)),
             with(samp, mean(phi.11 - phi.00)),
             with(samp, mean(phi.01 - phi.00)),
             with(samp, mean(phi.01.pos - phi.00.pos)))

    varests = c(with(samp, var(phi.ate)),
             with(samp, var(phi.1 - phi.0)),
             with(samp, var(phi.11 - phi.00)),
             with(samp, var(phi.01 - phi.00)),
             with(samp, var(phi.01.pos - phi.00.pos))) / sample_size

    res[[i]] = tibble(
      component = c("$\\tilde{\\Psi}$", "$\\Omega_{11}$",
                    "$\\Omega_{21}$", "$\\Omega_3$", "$\\Omega_{41}$"),
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

HardSimulation <- function(popdata, sample_size, number_sims, sl_lib) {


  popdata = ConstructPI(popdata)
  truth = c(with(popdata, mean(phi.ate.pi)),
            with(popdata, mean(phi.1.pi - phi.0.pi)),
            with(popdata, mean(phi.11.pi - phi.00.pi)),
            with(popdata, mean(phi.01.pi - phi.00.pi)),
            with(popdata, mean(phi.01.pos.pi - phi.00.pos.pi)))

  res = list()

  for (i in 1:number_sims) {
    samp = select(popdata, X, A, C, Y) %>%
      sample_n(sample_size) %>%
      mutate(id = 1:nrow(.),
             cluster_id = id)

    samp = EstimateNuisanceFunctions(X_dat = tibble(id = samp$id,
                                                    cluster_id = samp$id,
                                                    X = samp$X),
                                     outcome_tx_dat = samp %>% select(-X),
                                     outer_folds = 2,
                                     sl_folds = 2,
                                     sl_lib = sl_lib) %>%
      ConstructPI() %>%
      ConstructEIF()

    ests = c(with(samp, mean(phi.ate)),
             with(samp, mean(phi.1 - phi.0)),
             with(samp, mean(phi.11 - phi.00)),
             with(samp, mean(phi.01 - phi.00)),
             with(samp, mean(phi.01.pos - phi.00.pos)))

    varests = c(with(samp, var(phi.ate)),
                with(samp, var(phi.1 - phi.0)),
                with(samp, var(phi.11 - phi.00)),
                with(samp, var(phi.01 - phi.00)),
                with(samp, var(phi.01.pos - phi.00.pos))) / sample_size

    res[[i]] = tibble(
      component = c("$\\tilde{\\Psi}$", "$\\Omega_{11}$",
                    "$\\Omega_{21}$", "$\\Omega_3$", "$\\Omega_{41}$"),
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

Simulation <- function(population_data,
                       number_simulations,
                       sample_size,
                       estimand,
                       outcome_assumption,
                       monotonicity, truth,
                       delta1, delta0, tau1, tau0,
                       bound) {

  res = list()

  for (i in 1:number_simulations) {
    sample <- sample_n(population_data, sample_size)
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
  invoke(rbind, res)
}

AddNoise <- function(data, mu_rate, pi_rate, gamma_rate) {
  n = nrow(data)
  data %>%
    mutate_at(vars(matches("mu\\.a0")), ~expit(logit(.) + rnorm(n, 1/n^mu_rate, 1/n^mu_rate))) %>%
    mutate_at(vars(matches("mu\\.a1")), ~expit(logit(.) + rnorm(n, -2/n^mu_rate, 1/n^mu_rate))) %>%
    mutate_at(vars(matches("ps")), ~expit(logit(.) + rnorm(n, 2/n^pi_rate, 1/n^pi_rate))) %>%
    mutate_at(vars(matches("gamma\\.c1")), ~expit(logit(.) + rnorm(n, 1/n^gamma_rate, 1/n^gamma_rate))) %>%
    mutate_at(vars(matches("gamma\\.c0")), ~expit(logit(.) + rnorm(n, -2/n^gamma_rate, 1/n^gamma_rate))) #%>%
}


expit <- function(x) exp(x) / (1 + exp(x))

logit <- function(x) log(x / (1 - x))


