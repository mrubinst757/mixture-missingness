################################################################################
################### BOUNDS AND ESTIMATION FUNCTIONS ############################
################################################################################
#
# This file contains functions for:
# 1. Estimating nuisance functions (outcome models, censoring models, propensity scores)
#    using SuperLearner with cross-fitting
# 2. Constructing influence function-based (EIF) estimators of bound components
# 3. Constructing plug-in estimators of bound components
# 4. Computing sensitivity analysis bounds under various assumptions
#
# The main estimands supported are:
# - ATE: Average Treatment Effect
# - Z-ATE: Treatment effect on composite outcome (censoring or outcome)
# - SDE: Separable Direct Effect
#
################################################################################

library(tidyverse)
library(SuperLearner)
library(assertthat)

#' Estimate nuisance functions for double robust estimator of bounds
#'
#' @param outcome_tx_dat dataset with outcome and treatment variables.
#' @param X_dat covariate data
#' @param outer_folds number of folds for cross-validation and splitting
#' @param sl_lib list of SuperLearner libraries for SuperLearner
#' @param sl_folds number of folds for the SuperLearner
#'
#' @return dataset with estimated nuisance functions for bound components
EstimateNuisanceFunctions <- function(X_dat,
                                      outcome_tx_dat,
                                      outer_folds,
                                      sl_lib,
                                      sl_folds,
                                      folds_seed = NULL) {

  #######################
  #### Data Cleaning ####
  #######################

  # Sort data by cluster and individual ID for consistent processing
  X_dat = arrange(X_dat, cluster_id, id)
  outcome_tx_dat = arrange(outcome_tx_dat, cluster_id, id)
  clusters = unique(X_dat$cluster_id)
  number_per_cluster = as.numeric(table(X_dat$cluster_id))

  # Create cross-fitting folds at the cluster level (important for clustered data)
  # This ensures individuals from same cluster stay together
  set.seed(folds_seed)
  folds_agg <- sample(outer_folds, length(clusters), replace = TRUE)
  folds <- unlist(map2(folds_agg, number_per_cluster, ~rep(.x, .y)))
  set.seed(NULL)

  # Convert covariate data to model matrix (handles factors, interactions, etc.)
  X_dat <- X_dat %>%
    stats::model.matrix(~ . -1, data = .) %>%
    as.data.frame(.)

  # Rename columns to generic names for SuperLearner compatibility
  # Some SL algorithms have issues with special characters in variable names
  colnames(X_dat) <- paste0("V", seq(1, ncol(X_dat)))

  # Initialize storage for estimated nuisance functions
  # Will be filled in via cross-fitting below
  dat <- dplyr::select(outcome_tx_dat, id)
  dat$mu.a1.u0 <- NA; dat$mu.a0.u0 <- NA; dat$ps <- NA  # Outcome models and propensity score
  dat$gamma.c1 <- NA; dat$gamma.c0 <- NA                # Censoring models

  a_vals <- sort(unique(outcome_tx_dat$A))  # Unique treatment values (typically 0 and 1)

  ###################################
  ### Nuisance Function Estimation
  ###################################
  # Uses K-fold cross-fitting: train on K-1 folds, predict on held-out fold
  # This gives out-of-sample predictions for all observations

  for (FOLD in seq(1, outer_folds, 1)) {

    # Define train/test split for this fold
    TEST_ROWS <- folds == FOLD
    TRAIN_ROWS <- !TEST_ROWS

    # Special case: if only 1 fold, train and test on same data
    if (outer_folds == 1) {
      TRAIN_ROWS <- TEST_ROWS
    }

    ########################
    ### Pi: P(A = 1 | X) ###
    ########################
    # Estimate propensity score (probability of treatment given covariates)

    outcome <- outcome_tx_dat$A[TRAIN_ROWS]

    # Remove cluster_id and id columns before modeling
    Xdat.train = X_dat[TRAIN_ROWS, -c(1:2)]
    Xdat.test = X_dat[TEST_ROWS, -c(1:2)]

    # Handle edge case of single covariate (convert to dataframe)
    if (is.numeric(Xdat.train)) {
      Xdat.train = tibble(X = Xdat.train)
    }

    if (is.numeric(Xdat.test)) {
      Xdat.test = tibble(X = Xdat.test)
    }

    # Fit SuperLearner ensemble
    mod_pi <-
      SuperLearner::SuperLearner(
        Y = outcome,
        X = Xdat.train,
        SL.library = sl_lib,           # Ensemble of algorithms
        family = binomial(),            # Binary outcome
        cvControl = SuperLearner.CV.control(V = sl_folds)  # Internal CV for algorithm weighting
      )

    # Store out-of-sample predictions for test fold
    dat[["ps"]][TEST_ROWS] <- predict(mod_pi, newdata = Xdat.test)$pred[,1]

    ##################################
    ### Gamma: P(C = 1 | X, A = a) ###
    ##################################
    # Estimate censoring probability conditional on treatment and covariates
    # Fit separate model for each treatment level

    for (a in a_vals) {
      # Train on individuals who received treatment a, predict for all test observations
      outcome <- outcome_tx_dat$C[outcome_tx_dat$A == a & TRAIN_ROWS]

      Xdat.train = X_dat[outcome_tx_dat$A == a & TRAIN_ROWS, -c(1:2)]
      Xdat.test  = X_dat[TEST_ROWS, -c(1:2)]

      # Handle single covariate case
      if (is.numeric(Xdat.train)) {
        Xdat.train = tibble(X = Xdat.train)
      }

      if (is.numeric(Xdat.test)) {
        Xdat.test = tibble(X = Xdat.test)
      }

      # Fit SuperLearner for censoring model
      mod_gamma <-
        SuperLearner::SuperLearner(
          Y = outcome,
          X = Xdat.train,
          SL.library = sl_lib,
          family = binomial(),
          cvControl = SuperLearner.CV.control(V = sl_folds)
        )

      # Store predictions in column gamma.c0 or gamma.c1
      dat[[paste0("gamma.c", a)]][TEST_ROWS] <- predict(mod_gamma, newdata = Xdat.test)$pred[,1]
    }
    ###################################
    ### Mu: P(Y | X, A = a, C = 0) ####
    ###################################
    # Estimate outcome model among uncensored individuals (C=0)
    # This gives E[Y | X, A=a, C=0] = mu_a(X, 0) in the paper notation
    # Fit separate model for each treatment level

    for (a in a_vals) {
      # Train on uncensored individuals who received treatment a
      outcome <- outcome_tx_dat$Y[outcome_tx_dat$A == a & outcome_tx_dat$C == 0 & TRAIN_ROWS]

      # Detect outcome type (binary vs continuous)
      outcome_type = if_else(length(unique(outcome)) == 2, "Binary", "Continuous")

      # Use appropriate family for SuperLearner
      family_specification <- switch(outcome_type,
                                     "Binary" = binomial(),
                                     "Continuous" = gaussian())

      Xdat.train = X_dat[outcome_tx_dat$A == a & outcome_tx_dat$C == 0 & TRAIN_ROWS, -c(1:2)]
      Xdat.test  = X_dat[TEST_ROWS, -c(1:2)]

      # Handle single covariate case
      if (is.numeric(Xdat.train)) {
        Xdat.train = tibble(X = Xdat.train)
      }

      if (is.numeric(Xdat.test)) {
        Xdat.test = tibble(X = Xdat.test)
      }

      # Fit SuperLearner for outcome model
      mod_mu <-
        SuperLearner::SuperLearner(
          Y = outcome,
          X = Xdat.train,
          SL.library = sl_lib,
          family = family_specification,
          cvControl = SuperLearner.CV.control(V = sl_folds)
        )

      # Store predictions in column mu.a0.u0 or mu.a1.u0
      dat[[paste0("mu.a", a, ".u0")]][TEST_ROWS] <- predict(mod_mu, newdata = Xdat.test)$pred[,1]
    }
  }

  # Merge estimated nuisance functions back to original data
  outcome_tx_dat <- inner_join(outcome_tx_dat, dat, by = "id")
  assert_that(nrow(outcome_tx_dat) == nrow(X_dat), msg = "We lost people!?")

  outcome_tx_dat
}

#' ConstructEIF: given nuisance estimates, construct influence-function based estimators
#' of bound components
#'
#' @param data dataframe containing estimates of nuisance components with appropriate labelings
#' @param sd standard deviation for normal CDF approximation for indicator function
#'
#' @return dataframe with influence-function based estimates of bound components

ConstructEIF <- function(data, sd = 0.001) {

  ############################
  ### Basic ATE components ###
  ############################

  # Influence function for E[Y(1)]: treatment outcome
  # phi.y1.p1 is the empirical process term (centered at zero)
  data$phi.y1.p1 = with(data, ((A * (1 - C)) / (ps * (1 - gamma.c1))) * (Y - mu.a1.u0))
  # phi.y1 is the one-step estimator
  data$phi.y1 = with(data, phi.y1.p1 + mu.a1.u0)

  # Influence function for E[Y(0)]: control outcome
  data$phi.y0.p1 = with(data, (((1 - A) * (1 - C)) / ((1 - ps) * (1 - gamma.c0))) * (Y - mu.a0.u0))
  data$phi.y0 = with(data, phi.y0.p1 + mu.a0.u0)

  # Influence function for biased ATE (observable under MAR)
  data$phi.ate = with(data, phi.y1 - phi.y0)

  ############################
  ### Bound components     ###
  ############################

  # Omega_3: E[gamma_1(X) * mu_1(X, 0)] - outcome*censoring interaction under treatment
  # Three parts: empirical process for Y, empirical process for C, and plug-in
  data$phi.11.p1 <- with(data, (C == 0 & A == 1) / ((1 - gamma.c1) * ps) * (Y - mu.a1.u0) * gamma.c1)
  data$phi.11.p2 <- with(data, ((A / ps) * (C - gamma.c1) * mu.a1.u0))
  data$phi.11.p3 <- with(data, mu.a1.u0 * gamma.c1)
  data$phi.11 <- with(data, phi.11.p1 + phi.11.p2 + phi.11.p3)

  # E[gamma_0(X) * mu_1(X, 0)] - cross-treatment interaction (for SDE)
  data$phi.01.p1 <- with(data, ((C == 0 & A == 1) / (ps * (1 - gamma.c1))) * (Y - mu.a1.u0) * gamma.c0)
  data$phi.01.p2 <- with(data, (((1-A) / (1-ps)) * (C - gamma.c0) * mu.a1.u0))
  data$phi.01.p3 <- with(data, mu.a1.u0 * gamma.c0)
  data$phi.01 <- with(data, phi.01.p1 + phi.01.p2 + phi.01.p3)

  # Omega_4: E[gamma_0(X) * mu_0(X, 0)] - outcome*censoring interaction under control
  data$phi.00.p1 <- with(data, (C == 0 & A == 0) / ((1 - gamma.c0) * (1 - ps)) * (Y - mu.a0.u0) * gamma.c0)
  data$phi.00.p2 <- with(data, (((1-A) / (1-ps)) * (C - gamma.c0) * mu.a0.u0))
  data$phi.00.p3 <- with(data, mu.a0.u0 * gamma.c0)
  data$phi.00 <- with(data, phi.00.p1 + phi.00.p2 + phi.00.p3)

  # Omega_1: E[gamma_1(X)] - censoring probability under treatment
  data$phi.1.p1 <- with(data, (A / ps) * (C - gamma.c1))
  data$phi.1.p2 <- with(data, gamma.c1)
  data$phi.1 <- with(data, phi.1.p1 + phi.1.p2)

  # Omega_2: E[gamma_0(X)] - censoring probability under control
  data$phi.0.p1 <- with(data, ((1-A) / (1-ps)) * (C - gamma.c0))
  data$phi.0.p2 <- with(data, gamma.c0)
  data$phi.0 <- with(data, phi.0.p1 + phi.0.p2)

  ########################################
  ### SDE components with monotonicity ###
  ########################################
  # These components multiply by indicator functions (positive vs negative effect)
  # Approximated using normal CDF/PDF with small sd to avoid discontinuities

  # Omega_5 under positive monotonicity: E[gamma_0(X) * mu_1(X,0) * I{mu_1(X,0) > mu_0(X,0)}]
  data$phi.01.pos = with(data, phi.01 * pnorm(mu.a1.u0 - mu.a0.u0, sd = sd) +
                               phi.01.p3  * dnorm(mu.a1.u0 - mu.a0.u0, sd = sd) * (phi.01.p1 - phi.00.p1))

  # Omega_5 under negative monotonicity: E[gamma_0(X) * mu_1(X,0) * I{mu_1(X,0) < mu_0(X,0)}]
  data$phi.01.neg = with(data, phi.01 * pnorm(mu.a0.u0 - mu.a1.u0, sd = sd) +
                               phi.01.p3  * dnorm(mu.a0.u0 - mu.a1.u0, sd = sd) * (phi.00.p1 - phi.01.p1))

  # Corresponding terms for phi.00 (control outcome * censoring)
  data$phi.00.pos = with(data, phi.00 * pnorm(mu.a1.u0 - mu.a0.u0, sd = sd) +
                               phi.00.p3  * dnorm(mu.a1.u0 - mu.a0.u0, sd = sd) * (phi.01.p1 - phi.00.p1))

  data$phi.00.neg = with(data, phi.00 * pnorm(mu.a0.u0 - mu.a1.u0, sd = sd) +
                               phi.00.p3  * dnorm(mu.a0.u0 - mu.a1.u0, sd = sd) * (phi.00.p1 - phi.01.p1))

  data
}

#' ConstructPI: Constructs plug-in (non-influence function) estimates of bound components
#'
#' This function computes simpler plug-in estimators that just substitute estimated
#' nuisance functions into the target parameters. These are not influence function-based,
#' so they don't have the same asymptotic efficiency, but they are useful for:
#' 1. Quick point estimates when standard errors aren't needed
#' 2. Computing population truth values in simulations
#'
#' @param data Dataframe containing estimates of nuisance components with appropriate labelings
#' @param sd Standard deviation for normal CDF approximation for indicator function (default 0.001)
#'
#' @return Dataframe with plug-in estimates of bound components (columns with .pi suffix)
ConstructPI <- function(data, sd = 0.001) {
  # Biased ATE (observable under MAR): psi_tilde
  data$phi.ate.pi = with(data, mu.a1.u0 - mu.a0.u0)

  # Bound components (Omegas in paper notation)
  data$phi.11.pi = with(data, mu.a1.u0 * gamma.c1)  # Omega_3
  data$phi.01.pi = with(data, mu.a1.u0 * gamma.c0)  # Cross-term for SDE
  data$phi.00.pi = with(data, mu.a0.u0 * gamma.c0)  # Omega_4

  data$phi.1.pi = with(data, gamma.c1)  # Omega_1
  data$phi.0.pi = with(data, gamma.c0)  # Omega_2

  # SDE components with indicator approximation
  data$phi.01.pos.pi = with(data, phi.01.pi * pnorm(mu.a1.u0 - mu.a0.u0, sd = sd))

  data$phi.01.neg.pi = with(data, phi.01.pi * pnorm(mu.a0.u0 - mu.a1.u0, sd = sd))

  data$phi.00.pos.pi = with(data, phi.00.pi * pnorm(mu.a1.u0 - mu.a0.u0, sd = sd))

  data$phi.00.neg.pi = with(data, phi.00.pi * pnorm(mu.a0.u0 - mu.a1.u0, sd = sd))

  data
}

#' CensoredBounds: generates estimates of specified bounds given input sensitivity
#' parameters and assumptions
#'
#' @param eif_data dataframe containing influence functions of bound components
#' @param estimand target estimand: "ate" (treating Y as the outcome); "z-ate" (treating Y or U2 as the outcome), or "sde" (estimating a separable direct effect)
#' @param estimator influence-function based estimator ("eif") or plug-in estimator ("pi). note: no confidence
#' intervals are output for the plug-in estimates
#' @param monotonicity assumption on sign of relationship of \mu_a(X, 0) to \mu_a(X, 1): "none", "positive", or "negative" (for "ate" only)
#' @param delta1 sensitivity parameter upper bounding \pi*_1(X) / \pi_1(X)
#' @param delta0 sensitivity parameter upper bounding \pi*_0(X) / \pi_0(X)
#' @param delta1l sensitivity parameter lower bounding \pi*_1(X) / \pi_1(X)
#' @param delta0l sensitivity parameter lower bounding \pi*_0(X) / \pi_0(X)
#' @param outcome_assumption assumption relating \mu_a(X, 0) to \mu_a(X, 1) via sensitivity parameters tau_a ("ate" only)
#' @param tau1 sensitivity parameter relating \mu_1(X, 1) to \mu_1(X, 0) ("ate" only)
#' @param tau0 sensitivity parameter relating \mu_0(X, 1) to \mu_0(X, 0) ("ate" only)
#' @param bound bound or point-identified parameter desired
#' @param simulation TRUE/FALSE indicating whether inputs are from simulation
#' @param truth if simulation is TRUE, input giving the true targeted parameter
#'
#' @return dataframe with estimated parameters and confidence intervals
CensoredBounds <- function(eif_data,
                      estimand = "ate",
                      estimator = "eif",
                      monotonicity = "none",
                      delta1 = 1,
                      delta0 = delta1,
                      delta1l = 0,
                      delta0l = delta1l,
                      outcome_assumption = "none",
                      tau1 = NULL,
                      tau0 = tau1,
                      bound = TRUE,
                      simulation = FALSE,
                      truth = NULL) {

  # Input validation
  assertthat::assert_that(estimand %in% c("ate", "z-ate", "sde"),
                          msg = "Estimand must be either ate, z-ate, or sde")

  assertthat::assert_that(monotonicity %in% c("none", "positive", "negative"),
                          msg = "Monotonicity must be either none, positive, or negative")

  assertthat::assert_that(outcome_assumption %in% c("none", "equality", "inequality"),
                          msg = "outcome_assumptions must be either none, equality, or inequality")

  if (outcome_assumption != "none") {
    assertthat::assert_that(!is.null(tau1),
                            msg = "At minimum tau1 must be specified")
  }

  ##############################################################################
  ### Main computation: switch based on (estimand, outcome_assumption, monotonicity)
  ##############################################################################

  ######################
  ### ATE estimand   ###
  ######################

  if (estimand == "ate") {
    # ATE: Average Treatment Effect E[Y(1) - Y(0)]

    if (outcome_assumption == "none") {
      # No assumptions on relationship between mu_a(X,0) and mu_a(X,1)

      if (monotonicity == "none") {
        # Assumption-free bounds (sharpest without additional structure)
        if (estimator == "eif") {
          lb = with(eif_data, phi.ate - delta1 * phi.11 + delta0 * (phi.00 - phi.0))
          ub = with(eif_data, phi.ate + delta1 * (phi.1 - phi.11) + delta0 * phi.00)
          res = GenResTab(lb, ub)
        }
        if (estimator == "pi") {
          lb = with(eif_data, phi.ate.pi - delta1 * phi.11.pi + delta0 * (phi.00.pi - phi.0.pi))
          ub = with(eif_data, phi.ate.pi + delta1 * (phi.1.pi - phi.11.pi) + delta0 * phi.00.pi)
          res = c(mean(lb), mean(ub))
        }
      }
      if (monotonicity == "positive") {
        # Positive monotonicity: mu_a(X, 1) >= mu_a(X, 0) for both a
        # This tightens the bounds by removing some terms
        if (estimator == "eif") {
          lb = with(eif_data, phi.ate - delta0 * (phi.0 - phi.00))
          ub = with(eif_data, phi.ate + delta1 * (phi.1 - phi.11))
          res = GenResTab(lb, ub)
        }
        if (estimator == "pi") {
          lb = mean(with(eif_data, phi.ate.pi - delta0 * (phi.0.pi - phi.00.pi)))
          ub = mean(with(eif_data, phi.ate.pi + delta1 * (phi.1.pi - phi.11.pi)))
          res = c(lb, ub)
        }
      }
      if (monotonicity == "negative") {
        # Negative monotonicity: mu_a(X, 1) <= mu_a(X, 0) for both a
        if (estimator == "eif") {
          lb = with(eif_data, phi.ate - delta1 * phi.11)
          ub = with(eif_data, phi.ate + delta0 * phi.00)
          res = GenResTab(lb, ub)
        }
        if (estimator == "pi") {
          lb = mean(with(eif_data, phi.ate.pi - delta1 * phi.11.pi))
          ub = mean(with(eif_data, phi.ate.pi + delta0 * phi.00.pi))
          res = c(lb, ub)
        }
      }
    }
    if (outcome_assumption == "equality") {
      # Equality assumption: mu_a(X,1) = tau_a * mu_a(X,0)
      # This gives point identification (no bounds)
      if (estimator == "eif") {
        est = with(eif_data, phi.ate + delta1 * phi.11 * (tau1 - 1) - delta0 * phi.00 * (tau0 - 1))
        res = GenResTab(lb = NULL, ub = NULL, est = est, estimand = "estimand")
      }
      if (estimator == "pi") {
        res = mean(with(eif_data, phi.ate.pi + delta1 * phi.11.pi * (tau1 - 1) - delta0 * phi.00.pi * (tau0 - 1)))
      }
    }
    if (outcome_assumption == "inequality") {
      # Inequality assumption: (1/tau_a) <= mu_a(X,1)/mu_a(X,0) <= tau_a
      # This gives bounds (tighter than assumption-free)
      if (estimator == "eif") {
        lb = with(eif_data, phi.ate + delta1 * phi.11 * ((1 - tau1) / tau1) - delta0 * phi.00 * (tau0 - 1))
        ub = with(eif_data, phi.ate + delta1 * phi.11 * (tau1 - 1) - delta0 * phi.00 * ((1 - tau0) / tau0))
        res = GenResTab(lb, ub)
      }
      if (estimator == "pi") {
        lb = with(eif_data, phi.ate.pi + delta1 * phi.11.pi * ((1 - tau1) / tau1) - delta0 * phi.00.pi * (tau0 - 1))
        ub = with(eif_data, phi.ate.pi + delta1 * phi.11.pi * (tau1 - 1) - delta0 * phi.00.pi * ((1 - tau0) / tau0))
        res = c(mean(lb), mean(ub))
      }
      }
    }

  #######################
  ### Z-ATE estimand  ###
  #######################
  # Z-ATE: Treatment effect on composite outcome Z = I(U2=1 or Y=1)

  if (estimand == "z-ate") {
        if (bound) {
          if (estimator == "eif") {
            lb = with(eif_data, phi.ate + delta1l * (phi.1 - phi.11) - delta0 * (phi.0 - phi.00))
            ub = with(eif_data, phi.ate + delta1  * (phi.1 - phi.11) - deltal0 * (phi.0 - phi.00))
            res = GenResTab(lb, ub)
          }
          if (estimator == "pi") {
            lb = with(eif_data, phi.ate.pi + delta1l * (phi.1.pi - phi.11.pi) - delta0 * (phi.0.pi - phi.00.pi))
            ub = with(eif_data, phi.ate.pi + delta1 * (phi.1.pi - phi.11.pi)  - delta0l * (phi.0.pi - phi.00.pi))
            res = c(mean(lb), mean(ub))
          }
        }
        if(!bound) {
          if (estimator == "eif") {
            est = with(eif_data, phi.ate + delta1 * (phi.1 - phi.11) - delta0 * (phi.0 - phi.00))
            res = GenResTab(lb = NULL, ub = NULL, est = est, estimand = "estimand")
          }
          if (estimator == "pi") {
            res = mean(with(eif_data, phi.ate.pi + delta1 * (phi.1.pi - phi.11.pi) - delta0 * (phi.0.pi - phi.00.pi)))
          }
        }
  }

  ######################
  ### SDE estimand   ###
  ######################
  # SDE: Separable Direct Effect
  # Treatment effect excluding the pathway through informative censoring

  if (estimand == "sde") {
    if (bound) {
      if (estimator == "eif") {
        lb  = with(eif_data, phi.ate - delta0 * (phi.01.pos - phi.00.pos))
        ub  = with(eif_data, phi.ate - delta0 * (phi.01.neg - phi.00.neg))
        res = GenResTab(lb, ub)
      }
      if (estimator == "pi") {
        ub = mean(with(eif_data, phi.ate.pi - delta0 * (phi.01.pi - phi.00.pi) * (pnorm(-phi.ate.pi, sd = 0.001))))
        lb = mean(with(eif_data, phi.ate.pi - delta0 * (phi.01.pi - phi.00.pi) * (pnorm(phi.ate.pi, sd = 0.001))))
        res = c(lb, ub)
      }
    }
    if (!bound) {
      if (estimator == "eif") {
        est = with(eif_data, phi.ate - delta0 * (phi.01 - phi.00))
        res = GenResTab(lb = NULL, ub = NULL, est = est, estimand = "estimand")
      }
      if (estimator == "pi") {
        res = mean(with(eif_data, phi.ate.pi - delta0 * (phi.01.pi - phi.00.pi)))
      }
    }
  }
  if (simulation) {
    res$truth = truth
  }
  res
}

#' GenResTab: helper function for CensoredBounds that returns a dataframe of
#' influence-function based estimates desired estimand
#'
#' @param lb lower bound (for estimands = "bounds" only)
#' @param ub upper bound (for estimands = "bounds" only)
#' @param est estimated estimand (for estimands = "estimand" only)
#' @param estimand "bounds" or "estimand"
#' @param alpha type I error rate (default five percent)
#'
#' @return dataframe with influence-function based estimates and corresponding
#' confidence intervals
GenResTab <- function(lb,
                      ub,
                      est = NULL,
                      estimand = "bounds",
                      alpha = 0.05) {

  if (estimand == "bounds") {
    res = tibble(
      estimands = c("lb", "ub"),
      ests = c(mean(lb), mean(ub)),
      varests = c(var(lb)/length(lb), var(ub)/length(lb)),
      lci = ests - qnorm(1-alpha/2) * sqrt(varests),
      uci = ests + qnorm(1-alpha/2) * sqrt(varests)
    )
  }

  if (estimand == "estimand") {
    res = tibble(
      estimands = estimand,
      ests = mean(est),
      varests = var(est)/length(est),
      lci = ests - qnorm(1-alpha/2) * sqrt(varests),
      uci = ests + qnorm(1-alpha/2) * sqrt(varests)
    )
  }
  res
}

