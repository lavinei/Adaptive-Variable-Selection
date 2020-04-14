## Standard plots for an AVS/bma comparison

library(zoo)
library(ggplot2)
library(tidyr)
library(directlabels)
library(arrayhelpers)
library(gridExtra)
library(dplyr)
library(tidyverse)

basepath = "~/Homework/SAMSI Bays Opt/Paper Modeling/Scripts/Adaptive-Variable-Selection/"

setwd(basepath)
source("./Scripts/plots_comparison_functions.R")

training_period = 60

filename = "lpfds"

imgfolder = paste0(basepath, "Images/", filename)

if(!dir.exists(imgfolder)){
  dir.create(imgfolder)  
}

# Load the BMA data
simulation = paste("bma", filename, sep="_")
load(paste0(basepath, "Simulations/", simulation, ".RData"))
diff = FALSE
plot_bma = get_plot_data(cred_int_store, forecast_sample_store, diff, training_period) %>% mutate("model" = "bma")
kdens_bma = get_kdens_data(k_step_density_store, time, "bma", training_period)
bma_envs = environments

# Make the standard individual plots
setwd(imgfolder)
if(!dir.exists(simulation)){
  dir.create(simulation)  
}
setwd(paste("./", simulation, sep=''))
plot_data = get_plot_data(cred_int_store, forecast_sample_store, diff, training_period)
if(!intervention){
  interventions = NULL
}

# Make the pdfs
for(i in 1:num_series){
  series_number = i
  make_standard_plots(plot_data, environments[[series_number]]$model_store, environments[[series_number]]$predictor_inclusion_probabilities_store, series_number, interventions, k_step_density_store, k = forecast_length, training_period, paste(names[series_number], ".pdf", sep=""))
}

# Load the AVS data
setwd(basepath)
# setwd(paste("../Images/", filename, sep=""))
simulation = paste("AVS", filename, sep="_")
load(paste("./Simulations/", simulation, ".RData", sep=''))
diff = FALSE
plot_avs = get_plot_data(cred_int_store, forecast_sample_store, diff, training_period) %>% mutate("model" = "avs")
kdens_avs = get_kdens_data(k_step_density_store, time, "avs", training_period)
avs_envs = environments

# Make the standard individual plots
setwd(imgfolder)
if(!dir.exists(simulation)){
  dir.create(simulation)  
}
setwd(paste("./", simulation, sep=''))
plot_data = get_plot_data(cred_int_store, forecast_sample_store, diff, training_period)
if(!intervention){
  interventions = NULL
}

# Make the pdfs
for(i in 1:num_series){
  series_number = i
  make_standard_plots(plot_data, environments[[series_number]]$model_store, environments[[series_number]]$predictor_inclusion_probabilities_store, series_number, interventions, k_step_density_store, k = forecast_length, training_period, paste(names[series_number], ".pdf", sep=""))
}

# Make the comparison plots
setwd(imgfolder)


# Plot the marginal rMSE
plot_data = bind_rows(plot_bma, plot_avs)
p = rMSE_marg_combined(plot_data, names, scales="free_y") + scale_x_continuous(name = "Horizon", breaks = seq(4, k, by=4))
make_jpg(p, paste(filename, "marginal_rMSE_combined", sep="_"))

# Plot the marginal percent error
# plot_data = bind_rows(plot_bma, plot_avs)
# p = perc_error_marg_combined(plot_data, names) + scale_x_continuous(name = "Forecast Length")
# make_jpg(p, paste(filename, "marginal_perc_err_combined", sep="_"))

# Plot the objective fxn relative to AVS over time
kdens_data = get_kdens_comparison_data(list(kdens_avs, kdens_bma))

# p = LPDR_logitudinal_relative(kdens_data) + scale_x_continuous("Time", breaks = tticks, labels = tdates) + scale_y_continuous(name = "log density of k-step forecast relative to AVS", limits=c(-100, 0))
# make_jpg(p, paste(filename, "k_step_path_density_comparison_relative", sep="_"))

# Plot the objective fxn over time
p = LPDR_logitudinal(kdens_data) + 
  scale_x_continuous("Time", breaks = tticks, labels = tdates) + 
  scale_y_continuous(name = "log density of k-step forecast")#, limits=c(-100, 0))

make_jpg(p, paste(filename, "k_step_path_density_comparison", sep="_"))






# Plot the cumulative objective fxn over time
kdens_data_cumsum = kdens_data %>% spread(key = "model", value = "k_step_path_density") %>%
  mutate(bma = bma - avs) %>%
  transmute(time = time,
            bma = cumsum(bma)) %>%
  gather("model", "LPDR", 2)
p = LPDR(kdens_data_cumsum) + 
  scale_x_continuous("Time", breaks = tticks, labels = tdates) + 
  scale_y_continuous(name = "log density of k-step forecast \n relative to AVS")
make_jpg(p, paste(filename, "k_step_path_density_comparison_cumulative", sep="_"))

# Variable inclusion plots - avs
plot_list = lapply(avs_envs, function(e) variable_inclusion_plot(e$predictor_inclusion_probabilities_store, e$xnames, training_period))
walk2(plot_list, names, function(p, n) make_jpg(p, paste(filename, "AVS_variable_inclusion", n, sep="_"), width=12, height=8))

# Variable inclusion plots - bma
plot_list = lapply(bma_envs, function(e) variable_inclusion_plot(e$predictor_inclusion_probabilities_store, e$xnames, training_period))
walk2(plot_list, names, function(p, n) make_jpg(p, paste(filename, "BMA_variable_inclusion", n, sep="_"), width=12, height=8))

# Model plots - avs
plot_list = lapply(avs_envs, function(e) model_plot(e$model_store, e$xnames, e$Omega, training_period))
walk2(plot_list, names, function(p, n) make_jpg(p, paste(filename, "AVS_model", n, sep="_"), width=12*.7, height=8*.7, res=150))

# Model plots - bma
plot_list = lapply(bma_envs, function(e) model_plot(e$model_store, e$xnames, e$Omega, training_period))
walk2(plot_list, names, function(p, n) make_jpg(p, paste(filename, "BMA_model", n, sep="_"), width=12*.7, height=8*.7, res=150))


# Professional variable inclusion plots for the simulation study
tticks = seq(0, 100, length = 5)
tdates = seq(0, 100, length = 5)
# filename = "simulation_t_plus_k_density"

kdens_data$time = kdens_data$time - training_period


avs_envs[[1]]$xnames = c("Intercept", expression(x[1]), expression(x[2]))
bma_envs[[1]]$xnames = c("Intercept", expression(x[1]), expression(x[2]))

plot_list = lapply(avs_envs, function(e) variable_inclusion_plot(e$predictor_inclusion_probabilities_store, e$xnames, training_period, simulation=TRUE))
plot_list[[1]] = plot_list[[1]] + scale_x_continuous("Time", breaks = tticks, labels = tdates)
make_jpg(plot_list[[1]], paste(filename, "AVS_variable_inclusion", sep="_"))

plot_list = lapply(avs_envs, function(e) model_plot_sim(e$model_store, e$xnames, e$Omega, training_period))
plot_list[[1]] = plot_list[[1]] + scale_x_continuous("Time", breaks = tticks, labels = tdates)
make_jpg(plot_list[[1]], paste(filename, "AVS_model", sep="_"))

plot_list = lapply(bma_envs, function(e) variable_inclusion_plot(e$predictor_inclusion_probabilities_store, e$xnames, training_period, simulation=TRUE))
plot_list[[1]] = plot_list[[1]] + scale_x_continuous("Time", breaks = tticks, labels = tdates)
make_jpg(plot_list[[1]], paste(filename, "bma_variable_inclusion", sep="_"))

plot_list = lapply(bma_envs, function(e) model_plot_sim(e$model_store, e$xnames, e$Omega, training_period))
plot_list[[1]] = plot_list[[1]] + scale_x_continuous("Time", breaks = tticks, labels = tdates)
make_jpg(plot_list[[1]], paste(filename, "bma_model", sep="_"))


# Forecast plots
# Multiple lengths at the same time - avs
plot_list = lapply(1:num_series, function(i) forecast_multiple_lengths(plot_data, series_number = i, names, mod = "avs", forecast_lengths = c(1, 12)))
walk2(plot_list, names, function(p, n) make_jpg(p, paste(filename, "AVS_forecast_k1_k12", n, sep="_")))

# Multiple lengths at the same time - bma
plot_list = lapply(1:num_series, function(i) forecast_multiple_lengths(plot_data, series_number = i, names, mod = "bma", forecast_lengths = c(1, 12)))
walk2(plot_list, names, function(p, n) make_jpg(p, paste(filename, "BMA_forecast_k1_k12", n, sep="_")))

# Comparing AVS and BMA forecasts
plot_list = lapply(1:num_series, function(i) forecast_multiple_mods(plot_data, i, names, mods = c("avs", "bma"), forecast_length = 12))
walk2(plot_list, names, function(p, n) make_jpg(p, paste(filename, "AVS_BMA_forecast_k12", n, sep="_")))
