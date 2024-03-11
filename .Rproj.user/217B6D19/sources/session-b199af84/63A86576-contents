library(tidyverse)

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

  ### Sort data
  X_dat = arrange(X_dat, cluster_id, id)
  outcome_tx_dat = arrange(outcome_tx_dat, cluster_id, id)
  clusters = unique(X_dat$cluster_id)
  number_per_cluster = as.numeric(table(X_dat$cluster_id))

  ### Create folds
  set.seed(folds_seed)
  folds_agg <- sample(outer_folds, length(clusters), replace = TRUE)
  folds <- unlist(map2(folds_agg, number_per_cluster, ~rep(.x, .y)))
  set.seed(NULL)

  ### Turn covariate data into model matrix
  X_dat <- X_dat %>%
    stats::model.matrix(~ . -1, data = .) %>%
    as.data.frame(.)

  ### Column naming hack to make models work with SuperLearner
  colnames(X_dat) <- paste0("V", seq(1, ncol(X_dat)))

  ### Initialize nuisance function columns

  dat <- dplyr::select(outcome_tx_dat, id)
  dat$mu.a1.u0 <- NA; dat$mu.a0.u0 <- NA; dat$ps <- NA
  dat$gamma.c1 <- NA; dat$gamma.c0 <- NA

  a_vals <- sort(unique(outcome_tx_dat$A))

  ###################################
  ### Nuisance Function Estimation
  ###################################

  for (FOLD in seq(1, outer_folds, 1)) {

    #cat("\n------------------------------------------------",
    #    "\n-- Estimating fold ", FOLD, " out of ", outer_folds,
    #    "\n------------------------------------------------")

    TEST_ROWS <- folds == FOLD
    TRAIN_ROWS <- !TEST_ROWS

    if (outer_folds == 1) {
      TRAIN_ROWS <- TEST_ROWS
    }

    ########################
    ### Pi: P(A = 1 | X) ###
    ########################

    #cat("\nEstimating Pi\n")

    outcome <- outcome_tx_dat$A[TRAIN_ROWS]

    Xdat.train = X_dat[TRAIN_ROWS, -c(1:2)]
    Xdat.test = X_dat[TEST_ROWS, -c(1:2)]

    if (is.numeric(Xdat.train)) {
      Xdat.train = tibble(X = Xdat.train)
    }

    if (is.numeric(Xdat.test)) {
      Xdat.test = tibble(X = Xdat.test)
    }

Y
    dat[["ps"]][TEST_ROWS] <- predict(mod_pi, newdata = Xdat.test)$pred[,1]

    #cat("\nEstimating Gamma\n")

    ##################################
    ### Gamma: P(C = 1 | X, A = a) ###
    ##################################

    ### Indicator for each treatment value
    for (a in a_vals) {
      outcome <- outcome_tx_dat$C[outcome_tx_dat$A == a & TRAIN_ROWS]

      Xdat.train = X_dat[outcome_tx_dat$A == a & TRAIN_ROWS, -c(1:2)]
      Xdat.test  = X_dat[TEST_ROWS, -c(1:2)]

      if (is.numeric(Xdat.train)) {
        Xdat.train = tibble(X = Xdat.train)
      }

      if (is.numeric(Xdat.test)) {
        Xdat.test = tibble(X = Xdat.test)
      }

      mod_gamma <-
        SuperLearner::SuperLearner(
          Y = outcome,
          X = Xdat.train,
          SL.library = sl_lib,
          family = binomial(),
          cvControl = SuperLearner.CV.control(V = sl_folds)
        )
      dat[[paste0("gamma.c", a)]][TEST_ROWS] <- predict(mod_gamma, newdata = Xdat.test)$pred[,1]
    }
    ###################################
    ### Mu: P(Y | X, A = a, C = 0) ####
    ###################################

    ### Indicator for each treatment value
    for (a in a_vals) {
      outcome <- outcome_tx_dat$Y[outcome_tx_dat$A == a & outcome_tx_dat$C == 0 & TRAIN_ROWS]

      outcome_type = if_else(length(unique(outcome)) == 2, "Binary", "Continuous")

      family_specification <- switch(outcome_type,
                                     "Binary" = binomial(),
                                     "Continuous" = gaussian())

      Xdat.train = X_dat[outcome_tx_dat$A == a & outcome_tx_dat$C == 0 & TRAIN_ROWS, -c(1:2)]
      Xdat.test  = X_dat[TEST_ROWS, -c(1:2)]

      if (is.numeric(Xdat.train)) {
        Xdat.train = tibble(X = Xdat.train)
      }

      if (is.numeric(Xdat.test)) {
        Xdat.test = tibble(X = Xdat.test)
      }

      mod_mu <-
        SuperLearner::SuperLearner(
          Y = outcome,
          X = Xdat.train,
          SL.library = sl_lib,
          family = family_specification,
          cvControl = SuperLearner.CV.control(V = sl_folds)
        )
      dat[[paste0("mu.a", a, ".u0")]][TEST_ROWS] <- predict(mod_mu, newdata = Xdat.test)$pred[,1]
    }
  }
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

  data$phi.y1.p1 = with(data, ((A * (1 - C)) / (ps * (1 - gamma.c1))) * (Y - mu.a1.u0))
  data$phi.y1 = with(data, phi.y1.p1 + mu.a1.u0)

  data$phi.y0.p1 = with(data, (((1 - A) * (1 - C)) / ((1 - ps) * (1 - gamma.c0))) * (Y - mu.a0.u0))
  data$phi.y0 = with(data, phi.y0.p1 + mu.a0.u0)

  data$phi.ate = with(data, phi.y1 - phi.y0)

  data$phi.11.p1 <- with(data, (C == 0 & A == 1) / ((1 - gamma.c1) * ps) * (Y - mu.a1.u0) * gamma.c1)
  data$phi.11.p2 <- with(data, ((A / ps) * (C - gamma.c1) * mu.a1.u0))
  data$phi.11.p3 <- with(data, mu.a1.u0 * gamma.c1)
  data$phi.11 <- with(data, phi.11.p1 + phi.11.p2 + phi.11.p3)

  data$phi.01.p1 <- with(data, ((C == 0 & A == 1) / (ps * (1 - gamma.c1))) * (Y - mu.a1.u0) * gamma.c0)
  data$phi.01.p2 <- with(data, (((1-A) / (1-ps)) * (C - gamma.c0) * mu.a1.u0))
  data$phi.01.p3 <- with(data, mu.a1.u0 * gamma.c0)
  data$phi.01 <- with(data, phi.01.p1 + phi.01.p2 + phi.01.p3)

  data$phi.00.p1 <- with(data, (C == 0 & A == 0) / ((1 - gamma.c0) * (1 - ps)) * (Y - mu.a0.u0) * gamma.c0)
  data$phi.00.p2 <- with(data, (((1-A) / (1-ps)) * (C - gamma.c0) * mu.a0.u0))
  data$phi.00.p3 <- with(data, mu.a0.u0 * gamma.c0)
  data$phi.00 <- with(data, phi.00.p1 + phi.00.p2 + phi.00.p3)

  data$phi.1.p1 <- with(data, (A / ps) * (C - gamma.c1))
  data$phi.1.p2 <- with(data, gamma.c1)
  data$phi.1 <- with(data, phi.1.p1 + phi.1.p2)

  data$phi.0.p1 <- with(data, ((1-A) / (1-ps)) * (C - gamma.c0))
  data$phi.0.p2 <- with(data, gamma.c0)
  data$phi.0 <- with(data, phi.0.p1 + phi.0.p2)

  data$phi.01.pos = with(data, phi.01 * pnorm(mu.a1.u0 - mu.a0.u0, sd = sd) +
                               phi.01.p3  * dnorm(mu.a1.u0 - mu.a0.u0, sd = sd) * (phi.01.p1 - phi.00.p1))

  data$phi.01.neg = with(data, phi.01 * pnorm(mu.a0.u0 - mu.a1.u0, sd = sd) +
                               phi.01.p3  * dnorm(mu.a0.u0 - mu.a1.u0, sd = sd) * (phi.00.p1 - phi.01.p1))

  data$phi.00.pos = with(data, phi.00 * pnorm(mu.a1.u0 - mu.a0.u0, sd = sd) +
                               phi.00.p3  * dnorm(mu.a1.u0 - mu.a0.u0, sd = sd) * (phi.01.p1 - phi.00.p1))

  data$phi.00.neg = with(data, phi.00 * pnorm(mu.a0.u0 - mu.a1.u0, sd = sd) +
                               phi.00.p3  * dnorm(mu.a0.u0 - mu.a1.u0, sd = sd) * (phi.00.p1 - phi.01.p1))

  data
}

#' ConstructPI: constructs plug-in estimates of bound components
#'
#' @param data dataframe containing estimates of nuisance components with appropriate labelings
#' @param sd standard deviation for normal CDF approximation for indicator function
#'
#' @return dataframe with plug-in estimates of bound components
ConstructPI <- function(data, sd = 0.001) {
  data$phi.ate.pi = with(data, mu.a1.u0 - mu.a0.u0)
  data$phi.11.pi = with(data, mu.a1.u0 * gamma.c1)
  data$phi.01.pi = with(data, mu.a1.u0 * gamma.c0)
  data$phi.00.pi = with(data, mu.a0.u0 * gamma.c0)

  data$phi.1.pi = with(data, gamma.c1)
  data$phi.0.pi = with(data, gamma.c0)

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

  if (estimand == "ate") {
    if (outcome_assumption == "none") {
      if (monotonicity == "none") {
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
      if (estimator == "eif") {
        est = with(eif_data, phi.ate + delta1 * phi.11 * (tau1 - 1) - delta0 * phi.00 * (tau0 - 1))
        res = GenResTab(lb = NULL, ub = NULL, est = est, estimand = "estimand")
      }
      if (estimator == "pi") {
        res = mean(with(eif_data, phi.ate.pi + delta1 * phi.11.pi * (tau1 - 1) - delta0 * phi.00.pi * (tau0 - 1)))
      }
    }
    if (outcome_assumption == "inequality") {
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

