source("R/plot-fns.R"); source("R/bound-fns.R"); library(ggrepel)

data = readRDS("sim-inputs/data.rds")
data
set.seed(123)

plot_data = data %>%
  mutate(data = map(data, ConstructPI)) %>%
  mutate(data = map(data, PlotData))

# Generate Plot 0
p0s = map(plot_data$data, Plot0)
ggsave("plots/plot-0.png", p0s[[7]], height = 5, width = 9)

# Generate Plot 1
p1s = map(plot_data$data, ~Plot1(.x, "ATE"))
p1sa = map(plot_data$data, ~Plot1(.x, c("ATE-Z", "SDE")))

ggsave("plots/plot-1.png", p1s[[7]], height = 5, width = 9)
ggsave("plots/plot-1a.png", p1sa[[7]], height = 5, width = 9)

# Generate Plot 2
p2s = map(plot_data$data, ~Plot2(.x, "ATE"))
ggsave("plots/plot-2.png", p2s[[7]], height = 5, width = 9)

# Generate Plot 3
p3s = map(plot_data$data, Plot3)
ggsave("plots/plot-3.png", p3s[[7]], height = 5, width = 9)
dat = plot_data$data[[7]]
# monotonicity, delta = 1
ate.lb.pos = mean(with(dat, phi.ate.pi - (phi.0.pi - phi.00.pi)))
ate.ub.pos = mean(with(dat, phi.ate.pi + (phi.1.pi - phi.11.pi)))
ate.lb.neg = mean(with(dat, phi.ate.pi - phi.11.pi))
ate.ub.neg = mean(with(dat, phi.ate.pi + phi.00.pi))

# monotonicity, delta = 0.7
ate.lb.pos7 = mean(with(dat, phi.ate.pi - 0.7 * (phi.0.pi - phi.00.pi)))
ate.ub.pos7 = mean(with(dat, phi.ate.pi + 0.7 * (phi.1.pi - phi.11.pi)))
ate.lb.neg7 = mean(with(dat, phi.ate.pi - 0.7 * phi.11.pi))
ate.ub.neg7 = mean(with(dat, phi.ate.pi + 0.7 * phi.00.pi))
c(ate.lb.neg7, ate.ub.neg7)
c(ate.lb.pos7, ate.ub.pos7)

ate.ub = mean(with(dat, phi.ate.pi + (phi.1.pi - phi.11.pi) + phi.00.pi))
ate.lb = mean(with(dat, phi.ate.pi - phi.11.pi - (phi.0.pi - phi.00.pi)))

# Generate Plot4
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

delta1 = 0.7; delta0 = 0.3; tau1 = tau0 =2
ate.ub = mean(with(dat, phi.ate.pi + delta1 * phi.11.pi * (tau1 - 1) - delta0 * phi.00.pi * ((1 - tau0) / tau0)))
ate.lb = mean(with(dat, phi.ate.pi + delta1 * phi.11.pi * ((1 - tau1) / tau1) - delta0 * phi.00.pi * (tau0 - 1)))


ggsave("plots/plot-4.png", p4s[[7]], height = 5, width = 9)
ggsave("plots/plot-4a.png", p4sa[[7]], height = 5, width = 9)

