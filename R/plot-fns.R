PlotData <- function(data) {
  data %>%
    mutate(mu1 = gamma.u2.1 * mu.a1.u1 + (1 - gamma.u2.1) * mu.a1.u0,
      mu0 = gamma.u2.0 * mu.a0.u1 + (1 - gamma.u2.0) * mu.a0.u0,
      mu1z = gamma.u2.1 + (1 - gamma.u2.1) * mu.a1.u0,
      mu0z = gamma.u2.0 + (1 - gamma.u2.0) * mu.a0.u0,
      cate = mu1 - mu0,
      catez = mu1z - mu0z,
      csde = (1 - gamma.u2.0) * (mu.a1.u0 - mu.a0.u0),
      csde.b = (1 - gamma.c0) * (mu.a1.u0 - mu.a0.u0)
    ) %>%
    sample_n(2e3)
}

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
    geom_line(lwd = 1) +
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
    guides(color="none")
}

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

Plot2 <- function(data, groups = c("ATE", "ATE-Z", "SDE")) {

  data$ate.ub = with(data, phi.ate.pi + (phi.1.pi - phi.11.pi) + phi.00.pi)
  data$ate.lb = with(data, phi.ate.pi - phi.11.pi - (phi.0.pi - phi.00.pi))

  data$zate.ub = with(data, phi.ate.pi - (phi.0.pi - phi.00.pi))
  data$zate.lb = with(data, phi.ate.pi + (phi.1.pi - phi.11.pi))

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

Plot3 <- function(data) {
  data$ate.ub = with(data, phi.ate.pi + (phi.1.pi - phi.11.pi) + phi.00.pi)
  data$ate.lb = with(data, phi.ate.pi - phi.11.pi - (phi.0.pi - phi.00.pi))

  data$ate.ub.pos = with(data, phi.ate.pi - (phi.0.pi - phi.00.pi))
  data$ate.lb.pos = with(data, phi.ate.pi + (phi.1.pi - phi.11.pi))

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

Plot4 <- function(data, tau1, tau0, delta1, delta1l, delta0, delta0l,
                  groups = c("ATE", "ATE-Z", "SDE")) {

  data$ate.ub = with(data, phi.ate.pi + delta1 * phi.11.pi * (tau1 - 1) - delta0 * phi.00.pi * ((1 - tau0) / tau0))
  data$ate.lb = with(data, phi.ate.pi + delta1 * phi.11.pi * ((1 - tau1) / tau1) - delta0 * phi.00.pi * (tau0 - 1))

  data$zate.lb = with(data, phi.ate.pi + delta1l * (phi.1.pi - phi.11.pi) - delta0 * (phi.0.pi - phi.00.pi))
  data$zate.ub = with(data, phi.ate.pi + delta1 * (phi.1.pi - phi.11.pi)  - delta0l * (phi.0.pi - phi.00.pi))

  data$sde.ub = with(data, phi.ate.pi - delta0 * (phi.01.pi - phi.00.pi) * ((phi.01.pi < phi.00.pi)))
  data$sde.lb = with(data, phi.ate.pi - delta0 * (phi.01.pi - phi.00.pi) * ((phi.01.pi > phi.00.pi)))

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
