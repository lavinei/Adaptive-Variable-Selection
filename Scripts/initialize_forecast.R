## Initialize data, hyperparameters, and storage variables

library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
sourceCpp("./Scripts/utils_C.cpp",rebuild=TRUE)
library(jsonlite)
library(R.matlab)
library(stringr)
library(Matrix)
library(mvtnorm)
source("./Scripts/utils.R")
source("./Scripts/sss.R")
source("./Scripts/forecast_stationary.R")
source("./Scripts/score_functions.R")
library(compiler)
cmpfile(infile = "./Scripts/utils.R", outfile = "./Scripts/utils_cmp")
loadcmp("./Scripts/utils_cmp")
library(BAS)


if("OpenBlasThreads" %in% installed.packages()[,1]){
  library("OpenBlasThreads")
  set_num_threads(2)
}

library(directlabels)
library(tidyverse)
library(corrplot)
library(parallel)
library(doParallel)
library(hashmap)
library(abind)
library(purrr)

if(simulation){
  series_names = "sim"
  filename = "sim"
  names = series_names
  csv_name = "simulation_data"
  source("./Scripts/startup_simulation.R") 
  tau = 1
  max_lag = NULL
}else{
  # Define the series being used, the start date, and maximum lag
  series_names = c("Inflation", "Consumption", "Treasury10Yr")
  
  # Higher dimensional example:
  # series_names = c("Inflation", "Consumption", "Wage", "BAA", "Gold", "M2", "Treasury10Yr")
  
  start_date = "1/1/1993"
  max_lag = 12
  choose_lags = TRUE
  lags = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
  source("./Scripts/startup_macroecon.R")
  filename = "t_plus_k_density_ar_longtraining"
  
  # Set tau so that the AVS score is on the same scale as a standard likelihood
  
  # tau = 1 / 12 # Use tau = 1/k if the score function is lpfds_analytic
  tau = 1 # Use tau = 1 if score function is the t_plus_k_density_ar
} 

# Defining the score function, which generates a method to score models over all historical data
# This is wrapped in an objective function, which pulls old information from a running leaderboard, to save computation

# Main 3 options for sim data (not AR model): lpfds_analytic, t_plus_k_density, t_plus_kvec_density_generator(c(6, 12))
# Main 3 options for real data (AR models): lpfds_analytic, t_plus_k_density_ar_wrapper, t_plus_kvec_density_ar_generator_wrapper(c(6, 12))
score_fxn = t_plus_k_density_ar
score_fxn_alltime = score_alltime_generator(score_fxn)
#score_fxn_alltime = score_alltime_generator(score_fxn, nsamps = 50) #Default is to use 50 samples, can uncomment this line to change that
objective = objective_generator(score_fxn_alltime)


# Setting constants
sss_frequency = 1 # How often to run SSS to explore new models
niter = 2 # SSS iterations at each timestep
alpha = 0.98 # Model discount factor
delta = 0.98 # State vector discount factor
beta = 0.98 # Stochastic DLM variance discount factor
k = 12 # k referes to the number of predictive steps used in the model score
forecast_length = 12 # Typically this is equal to k... this is how far forward we want to simulate, and what level of score to keep track of
nsamps = 100 # Samples to draw when calculating the model LPFDS score
forecast_samps = 800 # Sample to draw when making predictions - much faster, so using more MC samples

maxmodel = 4 #Maximum model size
maxmodel_vec = rep(maxmodel, length(series_names))#Max model size for each series

maxstore = 50000 #Maximum number of models to store in the leaderboard
quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)
prior_type = "uniform"
#select_z = "predictive"
select_z_dist = "stationary" #either "stationary" or anything else
select_z = "postmode"
ncores = 8

#Defining storage variables
forecast_sample_store = NULL
mean_store = NULL
cred_int_store = NULL
runtime_store = NULL
sss_runtime_store = NULL
forecast_runtime_store = NULL
select_runtime_store = NULL
k_step_density_store = NULL
k_step_density_split_store = NULL

# Setting flags
store_predictor_matrix = FALSE

# Defining the simulation function (different based on a few factors, including the score_fxn)
simulate_path = ifelse(simulation,
                       simulate_path_generator(score_fxn, initialize_dlm_list, update_dlm_list, update_X_list),
                       simulate_path_generator(score_fxn, initialize_dlm_list_stationary, update_dlm_list_stationary, update_X_list_ar))


simulate_path_MC = ifelse(is_bma,
                          simulate_path_MC_BMA_generator(simulate_path),
                          simulate_path_MC_generator(simulate_path))

intervention=FALSE
if(intervention){
  intervention_dates = c("1/1/2001", "1/1/2008")
  interventions = as.numeric(sapply(intervention_dates, function(date) which(dates == date)))
}else{
  interventions = NULL
}

parallel = TRUE


if(parallel){
  cl = tryCatch({makeForkCluster(nnodes = ncores)}, error = function(e) {makeForkCluster(nnodes = 1)})
}


# Running a timestep
timesteps = nrow(environments[[1]]$X_full) - forecast_length
end_date = dates[timesteps]
