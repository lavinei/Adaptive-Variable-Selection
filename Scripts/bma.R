
setwd("~/Homework/SAMSI Bays Opt/Paper Modeling/Scripts/Adaptive-Variable-Selection/")

is_bma = TRUE
simulation = FALSE # If false, running the macroeconomic example

source("./Scripts/initialize_forecast.R")

bma = list()


for (series_number in 1:num_series){
  if(series_number == 1){
    environments[[series_number]]$bmadata = data.frame(environments[[series_number]]$y_full,
                                                       environments[[series_number]]$X_full[,-1])
    names(environments[[series_number]]$bmadata) = c("Y", environments[[series_number]]$xnames[-1])
  }else{
    environments[[series_number]]$bmadata = data.frame(environments[[series_number]]$y_full,
                                                       environments[[series_number]]$X_full[,-1])
    names(environments[[series_number]]$bmadata) = c("Y", environments[[series_number]]$xnames[-1])
  }
  
}

for (t in 1:timesteps){
  
  start = Sys.time()
  # For t in 1:T
  #   for each series:
  #     - Select n models with BMA. Discount / put less weight on older observations.
  #     - Forecast forwards k steps (select 1 model for each series, combine into a DDNM, and forecast)
  #     - Use forward filtering to get the DLM model coefficients

  # Select n models with BMA
  for (series_number in 1:num_series){
    if(t == 1){
      # Just use SSS for time = 0, assign all models a flat probability... there's no data to judge with!
      bma[[series_number]] = sss(objective, constraint, cset(environments[[series_number]]$dlm$model), environments[[series_number]]$Omega, 2)  
      models = lapply(names(bma[[series_number]]), extract_model_dlm)
      models = lapply(models, function(x) x - 1)
      probs = as.vector(bma[[series_number]], mode="numeric")
      bma[[series_number]] = list(which = models, postprobs = probs, probne0=rep(1/length(environments[[series_number]]$Omega), length(environments[[series_number]]$Omega)))
    }else{
      n = time
      bma[[series_number]] = bas.lm(Y ~ ., data=environments[[series_number]]$bmadata[1:n,], weights = beta^((n-1):0), n.models=100000, method="MCMC", renormalize = FALSE,
                                     modelprior=tr.beta.binomial(1,1,maxmodel_vec[series_number]-1), prior='g-prior', alpha=n)
    }
    
    # Save the variable inclusion probabilities
    environments[[series_number]]$predictor_inclusion_probabilities = bma[[series_number]]$probne0
    environments[[series_number]]$predictor_inclusion_probabilities_store = rbind(environments[[series_number]]$predictor_inclusion_probabilities_store,
                                                                                  environments[[series_number]]$predictor_inclusion_probabilities)
  }

  # Forecast ahead with DDNMs
  data_list = lapply(environments, function(x) x$X_full)
  observations_list = lapply(environments, function(x) x$y_full)
  ncol = lapply(data_list, ncol)
  
  model_list = lapply(bma, function(x) x$which)
  model_probs_list = lapply(bma, function(x) x$postprobs)
  
  for(i in 1:num_series){
    model_list[[i]] = lapply(model_list[[i]], function(x) x + 1)
  }
  
  if(parallel){
    clusterExport(cl, varlist=c("environments", "model_list", "model_probs_list", "data_list", "observations_list", "ncol", "time"))
    forecast_list = parLapply(cl, 1:forecast_samps, function(z) simulate_path_MC(model_list, model_probs_list, k = forecast_length, data_list = data_list, observations_list, ncol, store_predictor_matrix))
  }else if(!parallel){
    forecast_list = lapply(1:forecast_samps, function(z) simulate_path_MC(model_list, model_probs_list, k = forecast_length, data_list = data_list, observations_list, ncol, store_predictor_matrix))
  }
  forecast_sample = sapply(forecast_list, function(x) x[[1]], simplify="array")
  k_step_density = log(mean(exp(sapply(forecast_list, function(x) x[[2]], simplify="array"))))

  print(k_step_density)
  k_step_density_store = c(k_step_density_store, k_step_density)

  cred_int = apply(forecast_sample, c(1, 2), function(x) quantile(x, quantiles))  
  
  # Get the Monte Carlo mean
  mean = apply(forecast_sample, c(1, 2), function(x) mean(x))
  

  if(parallel){
    clusterExport(cl, varlist=c("environments"))
  }
  
  # Now advance one timestep, and add in data
  time = time + 1
  for(series_number in num_series:1){
    environments[[series_number]]$data = matrix(environments[[series_number]]$X_full[1:time,], nrow=time)
    if(series_number == num_series){
      environments[[series_number]]$data = environments[[series_number]]$X_full[1:(time+1),]
    }
    
    environments[[series_number]]$y = environments[[series_number]]$y_full[1:time]
  }
  
  # Save the mean and credible intervals
  mean_store = abind(mean_store, array(mean, dim = c(dim(mean), 1)))
  cred_int_store = abind(cred_int_store, array(cred_int, dim = c(dim(cred_int), 1)))
  
  # Save the forecast samples
  forecast_sample_store = abind(forecast_sample_store, array(forecast_sample, dim=c(dim(forecast_sample), 1)))
  
  # Save the "chosen" model as the posterior mode from BMA, for comparison with AVS
  for(series_number in 1:num_series){
    environments[[series_number]]$model_store = c(environments[[series_number]]$model_store, list(bma[[series_number]]$which[[which.max(bma[[series_number]]$postprobs)]] + 1))
  }
  
  end = Sys.time()
  runtime_store = c(runtime_store, difftime(end, start, units="mins"))
  
  print(time)
  print(paste("Step", time, ":", dates[time]))
}

if(parallel){
  stopCluster(cl)
  gc()
}

quantile_names = c("lower", "cred_25", "prediction", "cred_75", "upper")
dimnames(cred_int_store) = list("quantile" = quantile_names, "k" = 1:forecast_length, "series" = 1:num_series, "time" = 1:time)
save(environments, forecast_sample_store, k_step_density_store, diff, cred_int_store, time, names, xnames, tdates, tticks,  k, beta, delta, num_series, runtime_store, max_lag, nsamps, quantiles, ncores, intervention, forecast_length, interventions, mean_store, maxmodel_vec, k_step_density_split_store, score_fxn, file = paste0("./Simulations/bma_", filename, ".RData"))
