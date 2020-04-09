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
  ### For analysis of the univariate simulated data
  series_names = "sim"
  filename = "sim"
  names = series_names
  csv_name = "simulation_data"
  source("./Scripts/startup_simulation.R") 
  tau = 1
  max_lag = NULL
  
  # Main 3 options for sim data (not AR model): lpfds_analytic, t_plus_k_density, t_plus_kvec_density_generator(c(6, 12))
  k = 25 # k referes to the number of predictive steps used in the model score
  forecast_length = 25 # Typically this is equal to k... this is how far forward we want to simulate, and what level of score to keep track of
  score_fxn = t_plus_k_density
  

}else{
  ### For analysis of the multivariate macroeconomic data
  # Define the series being used, the start date, and maximum lag
  # Series availabe: InterestRate, Inflation, Unemployment, Wage, Consumption, M2, M1, Treasury10Yr, BAA, Oil, Gold
  
  lpfds = FALSE #If FALSE, then using the t_plus_k_density score function
  higher_dimension = FALSE # If TRUE, then doing a 7-dimensional example
  
  filename = ifelse(lpfds, "lpfds_newobj", "t_plus_k_density_ar_newobj")
  print(filename)

  if(higher_dimension){
    series_names = c("Inflation", "Consumption", "Wage", "BAA", "Gold", "M2", "Treasury10Yr")
    lags = c(1, 3, 6, 12)
    filename = paste(filename, "_higherdim", sep="")
    k = 24 # k referes to the number of predictive steps used in the model score
    forecast_length = 24 # Typically this is equal to k... this is how far forward we want to simulate, and what level of score to keep track of
    
  }else{
    
    series_names = c("Inflation", "Consumption", "Treasury10Yr")
    lags = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
    k = 12 # k referes to the number of predictive steps used in the model score
    forecast_length = 12 # Typically this is equal to k... this is how far forward we want to simulate, and what level of score to keep track of
  }
  
  start_date = "1/1/1991"
  max_lag = 12
  choose_lags = TRUE
  
  source("./Scripts/startup_macroecon.R")
  
  # Main 3 options for real data (AR models): lpfds_analytic, t_plus_k_density_ar_wrapper, t_plus_kvec_density_ar_generator_wrapper(c(6, 12))
  
  if(lpfds){
    score_fxn = lpfds_analytic
    # tau=1/k for the LPFD score puts AVS score on the same scale as a standard likelihood
    tau = 1 / k
    
  }else{
    score_fxn = t_plus_k_density_ar_wrapper
    # tau=1 for the t+k density score puts AVS score on the same scale as a standard likelihood
    tau = 1
    
  }

} 

# Defining the score function, which generates a method to score models over all historical data
# This is wrapped in an objective function, which pulls old information from a running leaderboard, to save computation
score_fxn_alltime = score_alltime_generator(score_fxn)
#score_fxn_alltime = score_alltime_generator(score_fxn, nsamps = 50) #Default is to use 50 monte carlo samples, can uncomment this line to change that
objective = objective_generator(score_fxn_alltime)

# Setting constants
sss_frequency = 1 # How often to run SSS to explore new models
niter = 2 # SSS iterations at each timestep
alpha = 0.98 # Model discount factor
delta = 0.98 # State vector discount factor
beta = 0.98 # Stochastic DLM variance discount factor
nsamps = 100 # Samples to draw when calculating the model LPFDS score
forecast_samps = 500 # Sample to draw when making predictions - much faster, so using more MC samples

maxmodel = 4 #Maximum model size
maxmodel_vec = rep(maxmodel, length(series_names))#Max model size for each series

maxstore = 50000 #Maximum number of models to store in the leaderboard
quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)
prior_type = "uniform"
# Choose the posterior modal model as our representative model.
# Other option is 'select_z = predictive', to choose model closest in predictive distribution to the mixture of models
select_z = "postmode" 

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
#mixture = TRUE
store_predictor_matrix = FALSE

# Defining the forecast simulation functions (different based on a few factors, including the score_fxn)
if(simulation){
  update_X_list = update_X_list
}else(
  if(choose_lags){
    update_X_list = partial(update_X_list_ar_chooselags, lags = lags)
  }else{
    update_X_list = update_X_list_ar
  }
)

simulate_path = ifelse(simulation,
                       simulate_path_generator(score_fxn, initialize_dlm_list, update_dlm_list, update_X_list),
                       simulate_path_generator(score_fxn, initialize_dlm_list_stationary, update_dlm_list_stationary, update_X_list))


simulate_path_MC = ifelse(is_bma,
                          simulate_path_MC_BMA_generator(simulate_path),
                          simulate_path_MC_generator(simulate_path))

# Add specific interventions, where more uncertainty is included? Generally FALSE
intervention=FALSE
if(intervention){
  intervention_dates = c("1/1/2001", "1/1/2008")
  interventions = as.numeric(sapply(intervention_dates, function(date) which(dates == date)))
}else{
  interventions = NULL
}

# Number of cores for parallelization
ncores = 8
parallel = TRUE

if(parallel){
  cl = tryCatch({makeForkCluster(nnodes = ncores)}, error = function(e) {makeForkCluster(nnodes = 1)})
}

# Number of timesteps in the analysis
timesteps = nrow(environments[[1]]$X_full) - forecast_length
end_date = dates[timesteps]
