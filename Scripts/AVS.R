# Running adaptive variable selection with a DDNM
setwd("~/Homework/SAMSI Bays Opt/Paper Modeling/Scripts/Adaptive-Variable-Selection/")

is_bma = FALSE
simulation = FALSE # If false, running the macroeconomic example

source("./Scripts/initialize_forecast.R")

for(i in 1:timesteps){

  start = Sys.time()
  
  source("./Scripts/one_timestep.R")
    
  # Remove old, unused models from each leaderboard so they don't grow in memory
  for(series_number in 1:num_series){
    num_mods = nrow(environments[[series_number]]$leaderboard)
    if(num_mods > maxstore){
      environments[[series_number]]$leaderboard = environments[[series_number]]$leaderboard[1:maxstore,]
    }
  }

  # Store variables
  forecast_sample_store = abind(forecast_sample_store, array(forecast_sample, dim=c(dim(forecast_sample), 1)))
  mean_store = abind(mean_store, array(mean, dim = c(dim(mean), 1)))
  cred_int_store = abind(cred_int_store, array(cred_int, dim = c(dim(cred_int), 1)))
  k_step_density_store = c(k_step_density_store, k_step_density) #These are log densities
  print(k_step_density)
  
  end = Sys.time()
  runtime_store = c(runtime_store, difftime(end, start, units="mins"))
  
  for(series_number in 1:num_series){
    # Saving the chosen model at each time step
    environments[[series_number]]$model_store = c(environments[[series_number]]$model_store, environments[[series_number]]$model)
    
    #Emptying out the model_list so it doesn't grow in the memory
    if(time %% sss_frequency == 0){
      rm(model_list, envir = environments[[series_number]]) 
    }
    
  
  }
  
  print(paste("Step", time, ":", dates[time]))
  
}

if(parallel){
  stopCluster(cl)
  gc()
}


print(sum(runtime_store)/60)
print(sum(sss_runtime_store)/60)
print(sum(forecast_runtime_store)/60)
print(sum(select_runtime_store)/60)

######## Save data ##########

quantile_names = c("lower", "cred_25", "prediction", "cred_75", "upper")
dimnames(cred_int_store) = list("quantile" = quantile_names, "k" = 1:forecast_length, "series" = 1:num_series, "time" = 1:time)
save(environments, cred_int_store, time, names, xnames, tdates, tticks, prior_type, select_z, k, tau, maxmodel, niter, alpha, beta, delta, num_series, runtime_store, sss_runtime_store, forecast_runtime_store, select_runtime_store, max_lag, X_full, y_full, nsamps, quantiles, ncores, k_step_density_store, intervention, forecast_length, leaderboard, interventions, maxstore, forecast_samps, forecast_sample_store, mean_store, maxmodel_vec, select_z_dist, diff, k_step_density_split_store, score_fxn, file = paste0("./Simulations/AVS_", filename, ".RData"))
