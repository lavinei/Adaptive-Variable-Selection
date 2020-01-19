### This is the main function for running the simulations
############################################## SET UP ENVIRONMENTAL VARIABLES ##############################################
walk(environments, function(e) e$simulated_array = array(NA, dim = c(time + k, ncol(e$X_full), nsamps)))
walk(environments, function(e) e$sim_count = 1)

data_list = lapply(environments, function(x) x$data)
observations_list = lapply(environments, function(x) x$y)
ncol = lapply(data_list, ncol)

############################################## RUN SSS ##############################################
sss_start = Sys.time()
if(time %% sss_frequency == 0){
  
  for(series_number in 1:num_series){
    
    # Setting individual series maxmodel lengths
    maxmodel = maxmodel_vec[series_number]
    
    if(parallel){
      clusterExport(cl, varlist=c("environments", "k", "series_number", "time", "maxmodel"), envir=globalenv())
    }
    
    environments[[series_number]]$sss_output = sss(objective, constraint, cset(environments[[series_number]]$dlm$model), environments[[series_number]]$Omega, niter)
    
    #   Extracting the list of models, and the vector of log-probabilities
    environments[[series_number]]$model_names = names(environments[[series_number]]$sss_output)
    environments[[series_number]]$models = lapply(names(environments[[series_number]]$sss_output), extract_model_dlm)
    #environments[[series_number]]$probs = as.vector(environments[[series_number]]$sss_output, mode="numeric") #Not actually the probs yet
    environments[[series_number]]$s_z = as.vector(environments[[series_number]]$sss_output, mode="numeric")
    
    # Now saving all of these models to the leaderboard and removing any models that are entered twice
    environments[[series_number]]$leaderboard = environments[[series_number]]$leaderboard %>%
      bind_rows(data.frame("name" = environments[[series_number]]$model_names, "score" = as.numeric(environments[[series_number]]$s_z), "time" = time) %>%
                  mutate(name = as.character(name))) %>%
      group_by(name) %>% filter(time == max(time)) %>% # removing models that are in there twice at different times
      distinct %>% # Removing duplicate rows
      arrange(desc(time), desc(score)) #Arranging in descending order by time and then model score
    
    for(model in environments[[series_number]]$models){
      environments[[series_number]]$model_list[[paste("Model:", toString(model))]] = dynamic_linear_model(time = time, F = environments[[series_number]]$data,
                                                                                                          y = environments[[series_number]]$y,
                                                                                                          series_number = series_number,
                                                                                                          model = model,
                                                                                                          dlm = environments[[series_number]]$model_list[[paste("Model:", toString(model))]],
                                                                                                          delta = delta,
                                                                                                          beta = beta)
    }
    
    # Re-updating the chosen model by pulling it from the model list
    environments[[series_number]]$dlm = environments[[series_number]]$model_list[[paste("Model:",toString(environments[[series_number]]$dlm$model))]]
  }
  
}

# Log the SSS runtime
sss_end = Sys.time()
sss_runtime_store = c(sss_runtime_store, difftime(sss_end, sss_start, units="mins"))


############################################## DEFINE THE PRIOR ####################################################
walk(environments, function(e) e$nmod = length(e$model_names))
if(prior_type == "uniform"){
  walk(environments, function(e) e$prior = rep(1, e$nmod))
}else if(prior_type == "seed"){
  # Giving the chosen model a prior probability of 0.5
  walk(environments, function(e) e$current_model = match(toString(e$dlm$model), sapply(e$models, function(x) toString(x))))
  walk(environments, function(e) e$prior = rep(1, length(e$model_names))/((e$nmod-1)*2))
  walk(environments, function(e) e$prior[e$current_model] = 0.5)
}

walk(environments, function(e) e$s_z0 = e$s_z[e$current_model])
walk(environments, function(e) e$probs = exp(tau*(rescale_probs(e$s_z))))
walk(environments, function(e) e$model_probs = (e$prior * e$probs) / sum(e$prior * e$probs))

############################################## SIMULATE FORWARD WITH BMA ##############################################
model_names_list = list()
model_probs_list = list()
for(s in 1:num_series){
  model_names_list = c(model_names_list, list(environments[[s]]$model_names))
  model_probs_list = c(model_probs_list, list(environments[[s]]$model_probs))
}

data_list = lapply(environments, function(x) x$X_full)
observations_list = lapply(environments, function(x) x$y_full)
ncol = lapply(data_list, ncol)

forecast_start = Sys.time()
if(parallel){
  clusterExport(cl, varlist=c("environments", "model_names_list", "model_probs_list", "data_list", "observations_list", "ncol", "time"))
  forecast_list = parLapply(cl, 1:forecast_samps, function(z) simulate_path_MC(model_names_list, model_probs_list, k = forecast_length, data_list = data_list, observations_list, ncol, store_predictor_matrix))
}else if(!parallel){
  forecast_list = lapply(1:forecast_samps, function(z) simulate_path_MC(model_names_list, model_probs_list, k = forecast_length, data_list = data_list, observations_list, ncol, store_predictor_matrix))
}

forecast_sample = sapply(forecast_list, function(x) x[[1]], simplify="array")

# Monte carlo average is taken on the real scale, not the log scale
k_step_density = log(mean(exp(sapply(forecast_list, function(x) x[[2]], simplify="array"))))

# Log the forecasting runtime
forecast_end = Sys.time()
forecast_runtime_store = c(forecast_runtime_store, difftime(forecast_end, forecast_start, units="mins"))

# Get the Monte Carlo credible intervals and means
cred_int = apply(forecast_sample, c(1, 2), function(x) quantile(x, quantiles))  
mean = apply(forecast_sample, c(1, 2), function(x) mean(x))  


############################################## OBSERVE TIME t RESPONSE, UPDATE MODEL PROBABILITIES ##############################################
# Add the next time point onto the series (data is stored in X_full) - need to tack it on to the end of the series
time = time + 1
for(series_number in num_series:1){
  environments[[series_number]]$data = matrix(environments[[series_number]]$X_full[1:time,], nrow=time)
  if(series_number == num_series){
    environments[[series_number]]$data = environments[[series_number]]$X_full[1:(time+1),]
  }
  
  environments[[series_number]]$y = environments[[series_number]]$y_full[1:time]
}

data_list = lapply(environments, function(x) x$data)
observations_list = lapply(environments, function(x) x$y)
ncol = lapply(data_list, ncol)

# Update all models with the new observation
for(series_number in 1:num_series){
  
  # Updating the model list
  environments[[series_number]]$p_y = NULL
  for(model_name in environments[[series_number]]$model_names){
    environments[[series_number]]$model_list[[model_name]] = dynamic_linear_model(time = time,
                                                                                  F = data_list[[series_number]],
                                                                                  y = observations_list[[series_number]],
                                                                                  dlm = environments[[series_number]]$model_list[[model_name]],
                                                                                  delta = delta,
                                                                                  beta = beta)
    
    
  }
  
  # Re-updating the chosen model by pulling it from the model list
  environments[[series_number]]$dlm = environments[[series_number]]$model_list[[paste("Model:",toString(environments[[series_number]]$dlm$model))]]
  
}

# Assign posterior probability to all models after new observation
for(series_number in 1:num_series){

  if(parallel){
    clusterExport(cl, varlist=c("environments", "k", "series_number", "time"), envir=globalenv())
    environments[[series_number]]$s_z = parSapply(cl, environments[[series_number]]$models, function(model) objective(model))
  }else{
    environments[[series_number]]$s_z = sapply(environments[[series_number]]$models, function(model) objective(model))
  }
    
  
  # Updating the leaderboard
  environments[[series_number]]$leaderboard = environments[[series_number]]$leaderboard %>%
    bind_rows(data.frame("name" = environments[[series_number]]$model_names, "score" = as.numeric(environments[[series_number]]$s_z), "time" = time) %>%
                mutate(name = as.character(name))) %>%
    group_by(name) %>% filter(time == max(time)) %>% # removing models that are in there twice at different times
    distinct %>% # Removing duplicate rows
    arrange(desc(time), desc(score)) #Arranging in descending order by time and then model score
  
  # Updating the posterior model probabilities
  current_model = match(toString(environments[[series_number]]$dlm$model), sapply(environments[[series_number]]$models, function(x) toString(x)))
  s_z0 = environments[[series_number]]$s_z[current_model]
  environments[[series_number]]$probs = exp(tau*(rescale_probs(environments[[series_number]]$s_z)))
  environments[[series_number]]$post_probs = environments[[series_number]]$prior * environments[[series_number]]$probs
  environments[[series_number]]$post_probs = environments[[series_number]]$post_probs / sum(environments[[series_number]]$post_probs)

  model_probs_list[[series_number]] = environments[[series_number]]$post_probs
  
  # Re-setting the simulated array to store the next simulation
  environments[[series_number]]$simulated_array = array(NA, dim = c(time + k, ncol(environments[[series_number]]$X_full), nsamps))
  environments[[series_number]]$sim_count = 1
  
  # Store variable inclusion probabilities
  environments[[series_number]]$tmp = t(sapply(environments[[series_number]]$models, function(model) environments[[series_number]]$Omega %in% model))
  environments[[series_number]]$predictor_inclusion_probabilities = apply(matrix(apply(environments[[series_number]]$tmp, 2, function(binary_model) binary_model*environments[[series_number]]$post_probs), ncol=ncol[[series_number]]), 2, function(x) sum(x))
  environments[[series_number]]$predictor_inclusion_probabilities_store = rbind(environments[[series_number]]$predictor_inclusion_probabilities_store,
                                                                                environments[[series_number]]$predictor_inclusion_probabilities)
  
}

############################################## SELECT NEXT CHOSEN MODEL Z ##############################################
select_start = Sys.time()
if(select_z == "postmode"){
  for(series_number in 1:num_series){
    chosen_model = which.max(environments[[series_number]]$post_probs)
    environments[[series_number]]$model = environments[[series_number]]$models[chosen_model]  
    environments[[series_number]]$dlm = environments[[series_number]]$model_list[[environments[[series_number]]$model_names[chosen_model]]]
  }
}else if(select_z == "predictive"){
  # Re-run the k-step ahead simulations, with the updated model probabilities
  
  # First forecast ahead from the model averaging, either with or without the stationary distribution
  if(select_z_dist == "stationary"){
    simulated_pts = replicate(nsamps, simulate_MC_path_stationary(model_names_list, model_probs_list, k = forecast_length, data_list = data_list, observations_list, ncol, store_loglik=FALSE), simplify="array")  
  }else{
    simulated_pts = replicate(nsamps, simulate_MC_path(model_names_list, model_probs_list, k = k, data_list = data_list, observations_list = observations_list, ncol = ncol))  
  }
  
  if(parallel){
    clusterExport(cl, varlist=c("environments", "simulated_pts"))
  }
  for(series_number in 1:num_series){
    
    # Calculate the KL-divergence between all individual predictive model densities and the BMA predictions
    
    # Right now this is done one series at a time
    # Since forecasts are made jointly, across all series, this is done by:
    # Looking at each of the simulated futures in TURN
    # Conditioning on the other series, and then getting the LPFD score for each model
    # And then averaging the LPFD scores for each model across all possible simulated futures
    # Since this is a lot of computation, only doing 10 MC samples per model & simulated future combination
    if(parallel){
      clusterExport(cl, varlist=c("series_number"))
      environments[[series_number]]$densities = parSapply(cl, environments[[series_number]]$model_names, function(model_name) KL_indep(nsamps, environments[[series_number]]$model_list[[model_name]]$m, environments[[series_number]]$model_list[[model_name]]$C,
                                                                                                                                       environments[[series_number]]$model_list[[model_name]]$n, environments[[series_number]]$model_list[[model_name]]$s, environments[[series_number]]$model_list[[model_name]]$model,
                                                                                                                                       forecast_length, series_number, time+forecast_length, environments[[series_number]]$simulated_array, environments[[series_number]]$y, simulated_pts[,series_number,],
                                                                                                                                       num_series, ncol(environments[[series_number]]$X_full), delta, beta))
    }else{
      environments[[series_number]]$densities = sapply(environments[[series_number]]$model_names, function(model_name) KL_indep(nsamps, environments[[series_number]]$model_list[[model_name]]$m, environments[[series_number]]$model_list[[model_name]]$C,
                                                                                                                      environments[[series_number]]$model_list[[model_name]]$n, environments[[series_number]]$model_list[[model_name]]$s, environments[[series_number]]$model_list[[model_name]]$model,
                                                                                                                      forecast_length, series_number, time+forecast_length, environments[[series_number]]$simulated_array, environments[[series_number]]$y, simulated_pts[,series_number,],
                                                                                                                      num_series, ncol(environments[[series_number]]$X_full), delta, beta))
    }
    
    # Choose the model with the minimum KL-divergence from the BMA mixture
    # I could go ahead and work out the density of each of these simulated paths under the BMA mixture
    # And the density of each path under the individual models
    # And then choose the model with the minimum KL-divergence
    # But that is equivalent to choosing the model with the highest likelihood for these points

    chosen_model = which.max(environments[[series_number]]$densities)
    environments[[series_number]]$model = environments[[series_number]]$models[chosen_model]
    environments[[series_number]]$dlm = environments[[series_number]]$model_list[[environments[[series_number]]$model_names[chosen_model]]]
    
  }
  
}

######################################

select_end = Sys.time()
select_runtime_store = c(select_runtime_store, difftime(select_end, select_start, units="mins"))
