library(tidyverse); expit <- function(x) exp(x) / (1 + exp(x))

GenerateLongData <- function(number_units,
                             number_time,
                             delta1,
                             delta0) {
  sample_size = number_units * number_time
  tibble(
    id = rep(1:number_units, each = number_time),
    time = rep(1:number_time, number_units)
  ) %>%
    group_by(id) %>%
    mutate(Z = rnorm(1, 0, 1)) %>%
    group_by(time) %>%
    mutate(dx = rnorm(1, 0, 1),
           dg = 0.025 * time + runif(1, 0, 0.25)) %>%
    ungroup() %>%
    mutate(X = dx + Z + rnorm(sample_size, 0, 1),
           gamma.c1 = expit(-1 + dg*X) * (X < 0) + expit(-1 + dg*0.5*X - 0.15*dg*X^2) * (X >= 0),
           gamma.u2.1 = delta1 * gamma.c1,
           gamma.u1.1 = gamma.c1 * (1 - delta1) / (1 - delta1 * gamma.c1),
           gamma.c0 = expit(-2 + dg*0.3*X) * (X < 0) + expit(-2 + dg*0.3*X - 0.05*dg*X^2) * (X >= 0),
           gamma.u2.0 = delta0 * gamma.c0,
           gamma.u1.0 = gamma.c0 * (1 - delta0) / (1 - delta0 * gamma.c0),
           U1.1 = rbinom(sample_size, 1, gamma.u1.1),
           U1.0 = rbinom(sample_size, 1, gamma.u1.0),
           U2.1 = rbinom(sample_size, 1, gamma.u2.1),
           U2.0 = rbinom(sample_size, 1, gamma.u2.0)) %>%
    group_by(id) %>%
    mutate(U1.0l = cumsum(U1.0), U2.0l = cumsum(U2.0),
           U1.1l = cumsum(U1.1), U2.1l = cumsum(U2.1),
           U1.0 = if_else(U1.0l > 0, 1, 0),
           U2.0 = if_else(U2.0l > 0, 1, 0),
           U1.1 = if_else(U1.1l > 0, 1, 0),
           U2.1 = if_else(U2.1l > 0, 1, 0),
           A = rbinom(1, 1, plogis(Z))) %>%
    select(-ends_with("l")) %>%
    mutate(C0 = if_else(U1.0 == 1 | U2.0 == 1, 1, 0),
           C1 = ifelse(U1.1 == 1 | U2.1 == 1, 1, 0),
           U1 = A * U1.1 + (1 - A) * U1.0,
           U2 = A * U2.1 + (1 - A) * U2.0,
           C  = A * C1 + (1 - A) * C0) %>%
    arrange(time, id)
}

GenerateLongCEF <- function(data,
                        a_val,
                        tau,
                        coef = 1) {
  if (!is.null(data[['time']])) {
    max_time = max(data[['time']])
  }
  if (is.null(data[['time']]))  {
    max_time = 1; data[['time']] = 1
  }

  res.u1 = c(); res.u0 = c()

  for (time in 1:max_time) {
    coef_t = coef * time / max_time

    if (tau >= 1) {
      new = expit(coef_t * data$X[data$time == time])
      res.u1 = c(res.u1, new)
      res.u0 = c(res.u0, new / tau)
    }
    if (tau < 1) {
      new = expit(coef_t * data$X[data$time == time])
      res.u1 = c(res.u1, new)
      res.u0 = c(res.u0, new * tau)
    }
  }
  data[[paste0("mu.a", a_val, ".u0")]] = res.u0
  data[[paste0("mu.a", a_val, ".u1")]] = res.u1

  data
}

GenerateLongOutcomes <- function(data) {
  N = nrow(data)
  data %>%
    ungroup() %>%
    mutate(
      Y11 = rbinom(N, 1, mu.a1.u1),
      Y10 = rbinom(N, 1, mu.a1.u0),
      Y01 = rbinom(N, 1, mu.a0.u1),
      Y00 = rbinom(N, 1, mu.a0.u0)) %>%
    arrange(id, time) %>%
    group_by(id) %>%
    mutate(Y11.l = cumsum(Y11), Y10.l = cumsum(Y10),
           Y01.l = cumsum(Y01), Y00.l = cumsum(Y00),
           Y11 = if_else(Y11.l > 0, 1, 0),
           Y10 = if_else(Y10.l > 0, 1, 0),
           Y01 = if_else(Y01.l > 0, 1, 0),
           Y00 = if_else(Y00.l > 0, 1, 0)) %>%
    select(-ends_with("l")) %>%
    mutate(
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

data = GenerateLongData(number_units = 100000,
                        number_time  = 3,
                        delta1 = 1,
                        delta0 = 1) %>%
  GenerateLongCEF(a_val = 1, tau = 1, coef = 1) %>%
  GenerateLongCEF(a_val = 0, tau = 1, coef = 0.5)

data = GenerateLongOutcomes(data)

data0 = GenerateData(population_size = 1000000,
                     delta1 = 1, delta0 = 1) %>%
  GenerateCEF(1, 1, 1) %>%
  GenerateCEF(0, 1, 0.5) %>%
  GenerateOutcomes()

with(filter(data, time == 1), mean(Y1 - Y0))
with(filter(data, time == 1), mean(mu1 - mu0))
#====

tdata = data %>%
  filter(time == 1) %>%
  group_by(id) %>%
  mutate(
    psi.1 = replace_na(lag(cumprod(1 - mu.a1.u0)), 1) * mu.a1.u0,
    psi.0 = replace_na(lag(cumprod(1 - mu.a0.u0)), 1) * mu.a0.u0,
    psi.1.lb = psi.1 - replace_na(lag(cumprod(1 - mu.a1.u0)), 1) * gamma.c1 * mu.a1.u0,
    psi.1.ub = psi.1 + replace_na(lag(cumprod(1 - mu.a1.u0)), 1) * gamma.c1 * (1 - mu.a1.u0),
    psi.0.lb = psi.0 - replace_na(lag(cumprod(1 - mu.a0.u0)), 1) * gamma.c0 * mu.a0.u0,
    psi.0.ub = psi.0 + replace_na(lag(cumprod(1 - mu.a0.u0)), 1) * gamma.c0 * (1 - mu.a0.u0)
  ) %>%
  group_by(id) %>%
  summarize(psi = mean(psi.1 - psi.0),
            psi.lb = mean(psi.1.lb - psi.0.ub),
            psi.ub = mean(psi.1.ub - psi.0.lb),
            X = X[time == 1])

tdata %>%
  sample_n(1000) %>%
  gather(key, value, psi, psi.lb, psi.ub) %>%
  ggplot(aes(x = X, y = value, fill = key, color = key)) +
  geom_line()

data0 %>%
  sample_n(1000) %>%
  ConstructPI() %>%
  mutate(psi.lb = phi.ate.pi - phi.11.pi + phi.0.pi - phi.00.pi) %%
  gather(key, value, phi.ate.pi, psi.lb) %>%
  ggplot(aes(x = X, y = value, fill = key, color = key)) +
  geom_line()

with(filter(data, time == 3), mean(Y1 - Y0))


#----
data %>%
  filter(time < 4) %>%
  group_by(id) %>%
  mutate(
    p1 = replace_na(lag(cumprod(1 - mu.a1.u0)), 1),
    p0 = replace_na(lag(cumprod(1 - mu.a0.u0)), 1),
    psi.1 = replace_na(lag(cumprod(1 - mu.a1.u0)), 1) * mu.a1.u0,
    psi.0 = replace_na(lag(cumprod(1 - mu.a0.u0)), 1) * mu.a0.u0
  ) %>%
  select(time, p1, p0, psi.1, psi.0) %>%
  group_by(time) %>%
  summarize(psi.1 = mean(psi.1),
            psi.0 = mean(psi.0),
            psi = psi.1 - psi.0) %>%
  summarize(psi = sum(psi))

with(filter(data, time == 3), mean(Y1 - Y0))

