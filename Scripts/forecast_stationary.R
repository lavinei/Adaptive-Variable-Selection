### Re-writing forecasting functions to include a stationarity constraint
### Or more specifically, a stability constrain

### Forecasting Simulation method (for each Monte Carlo draw):
### 1. Choose a DDNM, represented by 1 model for each univariate series
### 2. For each series, simulate the state vector k-steps ahead. At each step, jointly test for stationarity, and reject non-stationary samples.
### 3. Jointly simulate k-steps ahead, conditioned on the simulated state vectors

# model_names_list and model_probs_list should be a list of length = num_series
# Each entry in the list should be another list of models

# Function to construct the VAR representation of a DDNM
# Decompose that representation to test for stationarity

########################### FUNCTIONAL VERSIONS OF SIMULATE_PATH_STATIONARY ################################


is_stationary <- function(dlm_list, max_lag, num_series){
  # Initialize the dense G matrix
  G = matrix(0, nrow = max_lag*num_series, ncol=max_lag*num_series)
  
  # Initialize the coefficients & contemporaneous predictors from each model filled in with zeroes
  # Ignore the intercept
  coefs = matrix(0, nrow=num_series, ncol=(max_lag+1)*num_series)
  for(series_number in 1:num_series){
    coefs[series_number,dlm_list[[series_number]]$model[-1] + series_number - 1] = dlm_list[[series_number]]$theta[-1]
  }
  
  # Set up the contemporaneous predictor matrix and get I-Gam
  IGam = diag(num_series) - coefs[, 1:num_series]
    
  # Go from DDNM to VAR representation for coefficient matrices by multiplying coef =  (I-Gam)^-1 * coef for each coef matrix
  # Record the locations & values for a sparse matrix representation as well
  G[1:num_series,] = solve(IGam, coefs[,-(1:num_series)])
  
  # Add in the identity matrices to G
  G[(num_series+1):(num_series*max_lag), 1:(num_series*max_lag - num_series)] = diag(num_series*max_lag - num_series)
  
  # I only need the LARGEST eigenvalue - which luckily I can extract just that
  evs = eigen(G, symmetric=FALSE, only.values=TRUE)

  # Test for stationarity
  return(ifelse(abs(evs$values[1]) < 1, TRUE, FALSE))
}

sample_theta <- function(dlm){
  with(dlm, theta <- m + sqrt(v/s)*t(chol(C)) %*% as.matrix(rnorm(p)))
}

sample_v <- function(dlm){
  with(dlm, v <- 1/rgamma(1, n/2, (n*s)/2))
}

calculate_constants <- function(dlm){
  with(dlm, p <- length(m))
  with(dlm, W <- C*(1-delta)/delta)
  with(dlm, L <- t(chol(W)))
  dlm
}

sample_params <- function(dlm){
  sample_v(dlm)
  sample_theta(dlm)
  dlm
}

initialize_dlm_list <- function(dlm_list){
  dlm_list = lapply(dlm_list, function(dlm) calculate_constants(dlm))
  dlm_list = lapply(dlm_list, function(dlm) sample_params(dlm))
  dlm_list
}

initialize_dlm_list_stationary <- function(dlm_list){
  
  dlm_list = lapply(dlm_list, function(dlm) calculate_constants(dlm))
  
  # Try 10,000 times to initialize stationary parameters
  for(j in 1:10000){
    # Sample theta, v
    dlm_list = lapply(dlm_list, function(dlm) sample_params(dlm))
    
    # Check for stationarity
    theta_list = lapply(dlm_list, function(x) x$theta[-1])
    model_list = lapply(dlm_list, function(x) x$model[-1])
    p_list = sapply(dlm_list, function(x) x$p - 1)
    
    if(is_stationaryC(theta_list, model_list, p_list, max_lag, num_series)){
      break
    }
    
  }
  
  dlm_list
}

update_theta <- function(dlm){
  with(dlm, theta <- theta_old + sqrt(v/s) * L %*% as.matrix(rnorm(p)))
}

update_v <- function(dlm){
  with(dlm, v <- v_old * beta / rbeta(1, beta*n/2, (1-beta)*n/2))
}

update_constants <- function(dlm){
  with(dlm, v_old <- v)
  with(dlm, theta_old <- theta)
  dlm
}

update_params <- function(dlm){
  update_v(dlm)
  update_theta(dlm)
  dlm
}

update_dlm_list <- function(dlm_list){
  dlm_list = lapply(dlm_list, function(dlm) update_constants(dlm))
  dlm_list = lapply(dlm_list, function(dlm) update_params(dlm))
  dlm_list
}

update_dlm_list_stationary <- function(dlm_list){
  # Simulate the precisions and the states, while checking for stationarity
  dlm_list = lapply(dlm_list, function(dlm) update_constants(dlm))
  model_list = lapply(dlm_list, function(x) x$model[-1])
  p_list = sapply(dlm_list, function(x) x$p - 1)
  
  # Attempt to update 200 times, and then just move on
  for(j in 1:200){
    # Attempt to update
    dlm_list = lapply(dlm_list, function(dlm) update_params(dlm))
    theta_list = lapply(dlm_list, function(x) x$theta[-1])

    # Check for stationarity
    if(is_stationaryC(theta_list, model_list, p_list, max_lag, num_series)){
      break
    }
  }
  
  dlm_list
}

update_X_list_ar <- function(X_list, sim_pts, i, data_list, series_number, num_series){
  # Add the simulated point to all other series in the data list with a lower series number
  if(series_number > 1){
    X_list[[series_number-1]] = c(1, sim_pts[i, series_number], X_list[[series_number]][-1])
  }else if(series_number == 1){
    X_list[[num_series]] = c(1, sim_pts[i,], X_list[[num_series]][-c(1, (ncol[[num_series]] - num_series + 1):ncol[[num_series]])])
  }
  
  X_list
}

update_X_list_ar_chooselags <- function(X_list, sim_pts, i, data_list, series_number, num_series, lags=c()){

  # Start with the true future predictor values
  if(series_number > 1){
    X_list[[series_number-1]] = data_list[[series_number-1]][time+i,]

    # Fill in the contemporaneous simulated values
    X_list[[series_number-1]][2:(2 + num_series - series_number)] = sim_pts[i,series_number:num_series]

    # Fill in the other simulated values
    for(z in 1:length(lags)){
      l = lags[z]
      if(l < i){
        spot = 3 + num_series*z - series_number
        X_list[[series_number-1]][spot:(spot + num_series - 1)] = sim_pts[i - l,]
      }
    }

  }else if(series_number == 1){
    X_list[[num_series]] = data_list[[num_series]][time+i+1,]
    
    # Fill in the other simulated values
    for(z in 1:length(lags)){
      l = lags[z]
      if(l < (i+1)){
        spot = 2 + num_series*(z-1)
        X_list[[num_series]][spot:(spot + num_series - 1)] = sim_pts[i - l + 1,]
      }
    }
  }





  X_list
}


update_X_list <- function(X_list, sim_pts, i, data_list, series_number, num_series){
  # Get the next X values
  if(series_number > 1){
    X_list[[series_number - 1]] = data_list[[series_number - 1]][time + i,]
  }else if(series_number == 1){
    X_list[[num_series]] = data_list[[num_series]][time + i + 1,]
  }
  
  X_list
}

calculate_loglik <- function(s, time, dlm, X_full, y_full){
  # Calculate the predictive mean conditioned on observed values
  X = X_full[time, dlm$model]
  ft = crossprod(X, dlm$theta)
  # Loglikelihood of observed values
  dnorm(y_full[time], ft, sqrt(dlm$v), log=TRUE)
}

simulate_pt <- function(X, dlm){
  ft = crossprod(X[dlm$model], dlm$theta)
  ft + sqrt(dlm$v) * rnorm(1)
}

get_ft <- function(dlm, data, time){
  crossprod(data[time, dlm$model], dlm$theta)
}

simulate_path_generator <- function(score_fxn, initialize_dlm_list, update_dlm_list, update_X_list){
  function(dlm_list, k = 12, data_list, observations_list, ncol, store_predictor_matrix = TRUE){
    
    sim_pts = matrix(NA, nrow = k, ncol = num_series)
    X_list = vector("list", length = num_series)
    X_list[[num_series]] = data_list[[num_series]][time+1, ]
    
    # First, analyze the objective function applied into the future (must be analytic)
    

    # objective = vapply(1:num_series, function(i) score_fxn(dlm_list[[i]], k, i, time+k, environments[[i]]$X_full, environments[[i]]$y_full), numeric(1))
    # Joint k-step density objective:
    objective = c()
    
    # Second, forecast into the future through simulation
    
    # Initialize stationary samples for the precision and the state evolution variance
    dlm_list = initialize_dlm_list(dlm_list)
    
    for(i in 1:k){
      dlm_list = update_dlm_list(dlm_list)
      
      # Simulate the observations
      for(series_number in num_series:1){
        
        ### Joint t+k density, for horizon-specific forecasting
        if(i == k){
          # Create copies of the predictor matrix and simulated points
          # on this round, need to fill in the true values
          if(series_number == num_series){
            X_list_dens = X_list
            sim_pts_dens = sim_pts
          }
          
          objective = c(objective, log(model_dist(X_list_dens[[series_number]][dlm_list[[series_number]]$model], dlm_list[[series_number]], y = observations_list[[series_number]][time+k])))
          
          sim_pts_dens[i, series_number] = observations_list[[series_number]][time+k]
          X_list_dens = update_X_list(X_list_dens, sim_pts_dens, i, data_list, series_number, num_series)
        }
        
        # Simulate a new point in the path
        sim_pts[i, series_number] = simulate_pt(X_list[[series_number]], dlm_list[[series_number]])
        
        X_list = update_X_list(X_list, sim_pts, i, data_list, series_number, num_series)
        
        # Save the full simulated data list
        if(i == k & store_predictor_matrix){
          environments[[series_number]]$simulated_array[,,environments[[series_number]]$sim_count] = data_list[[series_number]]
          environments[[series_number]]$sim_count = environments[[series_number]]$sim_count + 1
        }
      }
    }
    
    # return(list(sim_pts, sum(objective), objective))  
    return(list(sim_pts, sum(objective)))  
    
  }
}

simulate_path_MC_generator <- function(simulate_path){
  function(model_names_list, model_probs_list, k = 12, data_list, observations_list, ncol, store_predictor_matrix = TRUE){
    if(length(model_names_list) != num_series | length(model_probs_list) != num_series){
      print("model_list and model_probs_list must be lists with length equal to num_series.")
      break
    }
    
    dlm_list = list()
    for(i in 1:num_series){
      dlm_list = c(dlm_list, environments[[i]]$model_list[[sample(x = model_names_list[[i]], size = 1, prob = model_probs_list[[i]])]])
    }
    
    # If the dlms haven't been calculated up to time t, then they need to be, using the real beta and delta
    # This part of the code is here just in case, but it should never actually be triggered, and a warning will print if it does
    for(series_number in num_series:1){
      if(dlm_list[[series_number]]$time < time){
        print("this dlm needed updating... seems wrong")
        dlm_list[[series_number]] = dynamic_linear_model(time = time, F = data_list[[series_number]],
                                                         y = observations_list[[series_number]], model = dlm$model,
                                                         series_number = series_number, dlm = dlm_list[[series_number]],
                                                         delta = delta, beta = beta)
        
      }
    }
    
  return(simulate_path(dlm_list, k, data_list, observations_list, ncol, store_predictor_matrix))
  }
  
}

simulate_path_MC_BMA_generator <- function(simulate_path){
  function(model_list, model_probs_list, k = 12, data_list, observations_list, ncol, store_predictor_matrix = TRUE){
    if(length(model_list) != num_series | length(model_probs_list) != num_series){
      print("model_list and model_probs_list must be lists with length equal to num_series.")
      break
    }
    
    models = list()
    for(i in 1:num_series){
      models = c(models, sample(x = model_list[[i]], size = 1, prob = model_probs_list[[i]]))
    }
    
    # Now calculating the DLMs
    dlm_list = list()
    for(series_number in 1:num_series){
      
      dlm_list = c(dlm_list, dynamic_linear_model(time = time, F = environments[[series_number]]$data,
                                                  y = environments[[series_number]]$y,
                                                  model = models[[series_number]],
                                                  series_number = series_number,
                                                  delta = delta,
                                                  beta = beta))
    }
    
    return(simulate_path(dlm_list, k, data_list, observations_list, ncol, store_predictor_matrix))
  }
}
