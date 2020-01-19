## Simulation Study
rm(list = ls())
library(tidyverse)
library(mvtnorm)
setwd("~/Homework/SAMSI Bays Opt/Paper Modeling/Scripts/github_scripts/")

uniform_F <- function(size, l, min, max){
  out = matrix(runif(l * size, min, max) * sample(c(-1, 1), l * size, replace=TRUE), nrow=l, ncol = size)
  out
}


sim_theta <- function(size, level, amp, periods, l, k = 12){
  out = matrix(0, ncol=size, nrow=l)
  for(i in 1:size){
    s = rep(sinpi(seq(0, 2, length=ceiling(l / periods[i])))[-ceiling(l / periods[i])], periods[i]+1)
    out[,i] = level[i] + amp[i]*s[1:l]
  }
  out
}

plot_thetas <- function(theta){
  p = ncol(theta)
  plot_data = as.data.frame(theta) %>%
    mutate(row = row_number()) %>%
    gather("coef", "value", 1:p)
  ggplot(plot_data, aes(x = row, y = value, color=coef)) + geom_line()
}

l = 100 + 30 + 30 + 25 # analysis will be 100 time steps, 30 data points to define prior, 30 for training, forecast horizon of 25

size = 2
F = uniform_F(size, l, min = 1, max = 1)
hist(F)

amp = c(1, .01)
periods = c(2, 3)
levels = c(0, 3)
theta = sim_theta(size, level=levels, amp = amp, periods=periods, l, k = 12)

plot_thetas(theta)

y = apply(F*theta, 1, sum) + 0 #rnorm(l, mean = 0, sd = .1)
plot(y, type='l')

output = cbind(F, y)
setwd("~/Homework/SAMSI Bays Opt/Paper Modeling/Scripts/github_scripts/")
write_csv(as.data.frame(output), "./Initial Data/simulation_data.csv")

# Make professional plot of the simulated thetas
p = ncol(theta)
plot_data = as.data.frame(theta)
names(plot_data) = c("Theta1", "Theta2")
plot_data = plot_data %>%
  slice(-((l-30):l)) %>%
  mutate(row = row_number()) %>%
  gather("coef", "value", 1:p)
p = ggplot(plot_data, aes(x = row, y = value, color=coef)) + geom_line(size=1.5)
tticks = c(0, 30, 60, 90, 120, 150)
tdates = tticks - 30
p = standard_theme(p, ylab = "Coefficient Value") + scale_x_continuous("Time", breaks = tticks, labels = tdates)

make_jpg(p, "sim_thetas")
