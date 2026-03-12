################################################################################
############################# PLOTTING FUNCTIONS ###############################
################################################################################
#
# This file contains functions for creating visualizations of:
# - Nuisance functions (outcome models, censoring models, propensity scores)
# - Conditional treatment effects (CATE, CATE-Z, CSDE)
# - Bounds under various sensitivity analysis assumptions
#
################################################################################

#' PlotData: Prepare data for plotting by computing conditional estimands
#'
#' This function computes conditional (individual-level) versions of three estimands:
#' - CATE: conditional average treatment effect on outcome Y
#' - CATE-Z: conditional average treatment effect on composite outcome Z
#' - CSDE: conditional separable direct effect
#' Then samples 2000 observations for visualization.
#'
#' @param data Population data with nuisance functions (from ConstructPI)
#'
#' @return Sampled dataframe with added columns for conditional estimands
PlotData <- function(data) {
  data %>%
    mutate(
      # Conditional mean outcomes under treatment and control (marginalized over U2)
      mu1 = gamma.u2.1 * mu.a1.u1 + (1 - gamma.u2.1) * mu.a1.u0,
      mu0 = gamma.u2.0 * mu.a0.u1 + (1 - gamma.u2.0) * mu.a0.u0,

      # Conditional mean of composite outcome Z (censoring or outcome)
      mu1z = gamma.u2.1 + (1 - gamma.u2.1) * mu.a1.u0,
      mu0z = gamma.u2.0 + (1 - gamma.u2.0) * mu.a0.u0,

      # Conditional average treatment effects
      cate = mu1 - mu0,           # CATE: treatment effect on Y
      catez = mu1z - mu0z,        # CATE-Z: treatment effect on Z
      csde = (1 - gamma.u2.0) * (mu.a1.u0 - mu.a0.u0),    # CSDE: separable direct effect (true)
      csde.b = (1 - gamma.c0) * (mu.a1.u0 - mu.a0.u0)    # CSDE: biased version using gamma.c0
    ) %>%
    sample_n(2e3)  # Sample 2000 observations for plotting
}

#' Plot0: Visualize all nuisance functions
#'
#' Creates a faceted plot showing how all estimated nuisance functions vary
#' with the covariate X. Includes propensity score, censoring probabilities,
#' and outcome models.
#'
#' @param data Dataset from PlotData with nuisance functions
#'
#' @return ggplot object with three facets:
#'   - Propensity score model: P(A=1|X)
#'   - Censoring models: P(C=1|A,X) for both treatment levels
#'   - Outcome models: P(Y=1|A,C=0,X) and P(Y=1|X) for both treatment levels
Plot0 <- function(data) {
  data %>%
    select(X,
           `P(A = 1 | X)` = ps,
           `P(C = 1 | A = 1, X)` = gamma.c1,
           `P(C = 1 | A = 0, X)` = gamma.c0,
           `P(Y(1) = 1 | A = 1, C = 0, X)` = mu.a1.u0,
           `P(Y(0) = 1 | A = 0, C = 0, X)` = mu.a0.u0,
           `P(Y(1) = 1 | X)` = mu1,
           `P(Y(0) = 1 | X)` = mu0) %>%
    gather(Function, value, -X) %>%
    mutate(group = case_when(
      grepl("Y\\(1\\)|Y\\(0\\)", Function) ~ "Outcome models",
      grepl("C = ", Function) ~ "Censoring models",
      grepl("A =", Function) ~ "Propensity score model"
    )) %>%
    mutate(color = case_when(
      grepl("Y\\(1\\)", Function) ~ 1,
      grepl("Y\\(0\\)", Function) ~ 2,
      grepl("P\\(A = 1", Function) ~ 3,
      grepl("C = 1 \\| A = 1", Function) ~ 4,
      grepl("C = 1 \\| A = 0", Function) ~ 5
    )) %>%
    mutate(color = factor(color)) %>%
    group_by(Function) %>%
    mutate(
      label = if_else(value == max(value), Function, NA_character_)
    ) %>%
    ggplot(aes(x = X, y = value, group = Function, color = color, label = Function)) +
    geom_line(lwd = 1.5) +
    theme_minimal() +
    facet_wrap(~group) +
    ylab("") +
    ylim(c(0,1)) +
    xlim(c(-4,4)) +
    geom_label_repel(aes(label = label),
                     nudge_x = 1,
                     na.rm = TRUE,
                     xlim = c(-4,4),
                     ylim = c(0,1)) +
    guides(color="none") +
    theme(text = element_text(size = 14))
}

#' Plot1: Compare true vs biased conditional treatment effects
#'
#' Visualizes the true conditional estimands (CATE, CATE-Z, CSDE) alongside
#' their biased counterparts that would be estimated under MAR assumption.
#' Shows how bias varies across covariate values.
#'
#' @param data Dataset from PlotData
#' @param groups Character vector specifying which estimands to plot
#'   Options: "ATE", "ATE-Z", "SDE" (default: all three)
#'
#' @return ggplot object showing true vs biased estimands across X values
Plot1 <- function(data, groups = c("ATE", "ATE-Z", "SDE")) {
  data %>%
    select(X, CATE = cate, `Biased CATE` = phi.ate.pi,
           `CATE-Z` = catez, `Biased CATE-Z` = phi.ate.pi,
           `CSDE` = csde, `Biased CSDE` = csde.b) %>%
    gather(estimand, value, -X) %>%
    mutate(group = case_when(
      grepl("CATE$", estimand) ~ "ATE",
      grepl("CATE-Z$", estimand) ~ "ATE-Z",
      grepl("SDE$", estimand) ~ "SDE"
    )) %>%
    filter(group %in% groups) %>%
    mutate(`Inferential target` = case_when(
      grepl("Biased", estimand) ~ "Biased",
      TRUE ~ "True"
    )) %>%
    ggplot(aes(x = X, y = value, fill = `Inferential target`, color = `Inferential target`)) +
    geom_line(lwd = 1.1) +
    stat_smooth(method="lm", formula=y~1, se=FALSE, lty = "dashed", linewidth = 0.75)+
    theme_minimal() +
    facet_wrap(~group) +
    ylab("Conditional estimand \n (average estimand in dashed line)") +
    xlab("Covariate value") +
    scale_color_manual(values = c("firebrick", "skyblue"))
}

#' Plot2: Visualize assumption-free bounds
#'
#' Shows the sharp bounds (no sensitivity parameter assumptions) for three
#' estimands, plotted as functions of X. Compares bounds to true estimands.
#'
#' @param data Dataset from PlotData
#' @param groups Character vector specifying which estimands to plot
#'
#' @return ggplot object showing assumption-free bounds and true estimands
Plot2 <- function(data, groups = c("ATE", "ATE-Z", "SDE")) {

  # Compute assumption-free bounds for ATE (no monotonicity, outcome assumptions)
  data$ate.ub = with(data, phi.ate.pi + (phi.1.pi - phi.11.pi) + phi.00.pi)
  data$ate.lb = with(data, phi.ate.pi - phi.11.pi - (phi.0.pi - phi.00.pi))

  # Compute assumption-free bounds for Z-ATE
  data$zate.ub = with(data, phi.ate.pi - (phi.0.pi - phi.00.pi))
  data$zate.lb = with(data, phi.ate.pi + (phi.1.pi - phi.11.pi))

  # Compute assumption-free bounds for SDE (uses indicator for sign of treatment effect)
  data$sde.ub = with(data, phi.ate.pi - (phi.01.pi - phi.00.pi) * ((phi.01.pi < phi.00.pi)))
  data$sde.lb = with(data, phi.ate.pi - (phi.01.pi - phi.00.pi) * ((phi.01.pi > phi.00.pi)))

  gdata = tibble(
    estimand = c("ate.ub", "ate.lb", "zate.ub", "zate.lb", "sde.ub", "sde.lb", "ate", "zate", "sde"),
    yintercept = c(mean(data$ate.ub), mean(data$ate.lb), mean(data$zate.ub), mean(data$zate.lb),
                   mean(data$sde.ub), mean(data$sde.lb), mean(data$cate), mean(data$catez),
                   mean(data$csde)),
    group = c("ATE", "ATE", "ATE-Z", "ATE-Z", "SDE", "SDE", "ATE", "ATE-Z", "SDE"),
    Target = c(rep("Assumption-free bound", 6), rep("Estimand", 3))
  ) %>%
    filter(group %in% groups)

  data %>%
    select(X, CATE = cate, `CATE-UB` = ate.ub, `CATE-LB` = ate.lb,
           `ZCATE-LB` = zate.lb, `ZCATE-UB` = zate.ub, ZCATE = catez,
           SDE = csde, `SDE-LB` = sde.lb, `SDE-UB` = sde.ub) %>%
    gather(estimand, value, -X) %>%
    mutate(group = case_when(
      grepl("^CATE", estimand) ~ "ATE",
      grepl("ZCATE", estimand) ~ "ATE-Z",
      grepl("SDE", estimand) ~ "SDE"
    )) %>%
    mutate(Target = case_when(
      grepl("B$", estimand) ~ "Assumption-free bound",
      TRUE ~ "Estimand"
    )) %>%
    filter(group %in% groups) %>%
    ggplot(aes(x = X, y = value, group = estimand, color = Target)) +
    facet_wrap(~group) +
    geom_line(lwd = 1.1, alpha = 0.6) +
    theme_minimal() +
    geom_hline(data = gdata, aes(yintercept = yintercept,
                                 group = estimand, color = Target), lwd = 1.1,
               lty = "dashed", alpha = 0.6) +
    ylab("Conditional estimand") +
    xlab("Covariate value") +
    scale_color_manual(values = c("forestgreen", "skyblue"))
}

#' Plot3: Compare bounds under monotonicity assumptions
#'
#' Visualizes how positive and negative monotonicity assumptions tighten the
#' bounds for the ATE compared to the assumption-free case.
#'
#' @param data Dataset from PlotData
#'
#' @return ggplot object with two facets comparing bounds under:
#'   - Positive monotonicity: mu_a(X,1) >= mu_a(X,0)
#'   - Negative monotonicity: mu_a(X,1) <= mu_a(X,0)
Plot3 <- function(data) {
  # Assumption-free bounds (same as Plot2)
  data$ate.ub = with(data, phi.ate.pi + (phi.1.pi - phi.11.pi) + phi.00.pi)
  data$ate.lb = with(data, phi.ate.pi - phi.11.pi - (phi.0.pi - phi.00.pi))

  # Bounds under positive monotonicity
  data$ate.ub.pos = with(data, phi.ate.pi - (phi.0.pi - phi.00.pi))
  data$ate.lb.pos = with(data, phi.ate.pi + (phi.1.pi - phi.11.pi))

  # Bounds under negative monotonicity
  data$ate.ub.neg = with(data, phi.ate.pi + phi.00.pi)
  data$ate.lb.neg = with(data, phi.ate.pi - phi.11.pi)

  g1 = data %>%
    select(X, CATE = cate, `CATE-UB` = ate.ub, `CATE-LB` = ate.lb,
           `CATE-LB-POS` = ate.lb.pos, `CATE-UB-POS` = ate.ub.pos,
           `CATE-LB-NEG` = ate.lb.neg, `CATE-UB-NEG` = ate.ub.neg) %>%
    gather(estimand, value, -X) %>%
    mutate(group = case_when(
      grepl("POS", estimand) ~ "Monotonicity: Positive",
      grepl("NEG", estimand) ~ "Monotonicity: Negative",
      !grepl("POS|NEG", estimand) & grepl("UB|LB", estimand) ~ "Assumption-free",
      TRUE ~ "True target (CATE)"
    )) %>%
    filter(group %in% c("Monotonicity: Positive", "Assumption-free", "True target (CATE)")) %>%
    mutate(Group = "Monotonicity: Positive")

  g2 = data %>%
    select(X, CATE = cate, `CATE-UB` = ate.ub, `CATE-LB` = ate.lb,
           `CATE-LB-POS` = ate.lb.pos, `CATE-UB-POS` = ate.ub.pos,
           `CATE-LB-NEG` = ate.lb.neg, `CATE-UB-NEG` = ate.ub.neg) %>%
    gather(estimand, value, -X) %>%
    mutate(group = case_when(
      grepl("POS", estimand) ~ "Monotonicity: Positive",
      grepl("NEG", estimand) ~ "Monotonicity: Negative",
      !grepl("POS|NEG", estimand) & grepl("UB|LB", estimand) ~ "Assumption-free",
      TRUE ~ "True target (CATE)"
    )) %>%
    filter(group %in% c("Monotonicity: Negative", "Assumption-free", "True target (CATE)")) %>%
    mutate(Group = "Monotonicity: Negative")

  bind_rows(g1, g2) %>%
    ggplot(aes(x = X, y = value, group = estimand, fill = group, color = group)) +
    geom_line() +
    facet_wrap(~Group) +
    theme_minimal() +
    ylab("Estimand value")
}

#' Plot4: Compare parameterized bounds to assumption-free bounds
#'
#' Visualizes how specifying sensitivity parameters (delta, tau) tightens
#' the bounds compared to making no assumptions. Shows both sets of bounds
#' for comparison.
#'
#' @param data Dataset from PlotData
#' @param tau1, tau0 Outcome sensitivity parameters for treatment and control
#' @param delta1, delta0 Upper bounds on censoring sensitivity parameters
#' @param delta1l, delta0l Lower bounds on censoring sensitivity parameters
#' @param groups Character vector specifying which estimands to plot
#'
#' @return ggplot object comparing parameterized vs assumption-free bounds
Plot4 <- function(data, tau1, tau0, delta1, delta1l, delta0, delta0l,
                  groups = c("ATE", "ATE-Z", "SDE")) {

  # Parameterized bounds for ATE (using specified delta and tau values)
  data$ate.ub = with(data, phi.ate.pi + delta1 * phi.11.pi * (tau1 - 1) - delta0 * phi.00.pi * ((1 - tau0) / tau0))
  data$ate.lb = with(data, phi.ate.pi + delta1 * phi.11.pi * ((1 - tau1) / tau1) - delta0 * phi.00.pi * (tau0 - 1))

  # Parameterized bounds for Z-ATE (using delta range [delta_l, delta])
  data$zate.lb = with(data, phi.ate.pi + delta1l * (phi.1.pi - phi.11.pi) - delta0 * (phi.0.pi - phi.00.pi))
  data$zate.ub = with(data, phi.ate.pi + delta1 * (phi.1.pi - phi.11.pi)  - delta0l * (phi.0.pi - phi.00.pi))

  # Parameterized bounds for SDE
  data$sde.ub = with(data, phi.ate.pi - delta0 * (phi.01.pi - phi.00.pi) * ((phi.01.pi < phi.00.pi)))
  data$sde.lb = with(data, phi.ate.pi - delta0 * (phi.01.pi - phi.00.pi) * ((phi.01.pi > phi.00.pi)))

  # Assumption-free bounds for comparison (suffix ".f" for "free")
  data$ate.ub.f = with(data, phi.ate.pi + (phi.1.pi - phi.11.pi) + phi.00.pi)
  data$ate.lb.f = with(data, phi.ate.pi - phi.11.pi - (phi.0.pi - phi.00.pi))

  data$zate.ub.f = with(data, phi.ate.pi - (phi.0.pi - phi.00.pi))
  data$zate.lb.f = with(data, phi.ate.pi + (phi.1.pi - phi.11.pi))

  data$sde.ub.f = with(data, phi.ate.pi - (phi.01.pi - phi.00.pi) * ((phi.01.pi < phi.00.pi)))
  data$sde.lb.f = with(data, phi.ate.pi - (phi.01.pi - phi.00.pi) * ((phi.01.pi > phi.00.pi)))

  data %>%
    select(X, CATE = cate, ZCATE = catez, SDE = csde,
           `CATE-UB` = ate.ub, `CATE-LB` = ate.lb,
           `ZCATE-LB` = zate.lb, `ZCATE-UB` = zate.ub,
           `SDE-LB` = sde.lb, `SDE-UB` = sde.ub,
           `CATE-UB-F` = ate.ub.f, `CATE-LB-F` = ate.lb.f,
           `ZCATE-LB-F` = zate.lb.f, `ZCATE-UB-F` = zate.ub.f,
           `SDE-LB-F` = sde.lb.f, `SDE-UB-F` = sde.ub.f) %>%
    gather(estimand, value, -X) %>%
    mutate(group = case_when(
      grepl("^CATE", estimand) ~ "ATE",
      grepl("ZCATE", estimand) ~ "ATE-Z",
      grepl("SDE", estimand) ~ "SDE"
    )) %>%
    filter(group %in% groups) %>%
    mutate(Target = case_when(
      grepl("B$", estimand) ~ "Parameterized bound",
      grepl("F$", estimand) ~ "Assumption-free bound",
      TRUE ~ "Estimand"
    )) %>%
    ggplot(aes(x = X, y = value, group = estimand, color = Target)) +
    facet_wrap(~group) +
    geom_line(lwd = 1.1, alpha = 0.6) +
    theme_minimal() +
    ylab("Conditional estimand") +
    xlab("Covariate value") +
    scale_color_brewer(palette = "Dark2")
    #scale_color_manual(values = c("forestgreen", "", "skyblue"))
}
