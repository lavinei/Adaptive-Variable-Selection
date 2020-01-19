################################## Functions #####################################

# Adaptive Variable Selection Functions
objective_generator <- function(score_fxn_alltime){
  function(set){
  
    # Turns out that R doesn't automatically copy environments... it's just passing a pointer here
    # So I make a new list of the chosen environments in 'one_timestep' before running SSS
    model = as.numeric(set)
    
    if(time == 0){
      score = 1
    }else if(time <= k){
      # Checking to see if this model is in the leaderboard
      old_score = (environments[[series_number]]$leaderboard %>% filter(name == paste("Model:",toString(model))) %>%
                     pull(score))[1]
      old_time = (environments[[series_number]]$leaderboard %>% filter(name == paste("Model:",toString(model))) %>%
                    pull(time))[1]
      
      # If it is in the leaderboard, use it!
      if(!is.na(old_score) & time == old_time){
        score = old_score
      }else{
        score = score_fxn_alltime(model = model, k=time, series_to_eval = series_number,
                                       start_time = time, end_time = time, data=environments[[series_number]]$data,
                                       observations = environments[[series_number]]$y)
      }
      
      
    }else if(time > k){
      # Checking to see if this model is in the leaderboard
      old_score = (environments[[series_number]]$leaderboard %>% filter(name == paste("Model:",toString(model))) %>%
                     pull(score))[1]
      old_time = (environments[[series_number]]$leaderboard %>% filter(name == paste("Model:",toString(model))) %>%
                    pull(time))[1]
      
      # If it is in the leaderboard, use it!
      if(is.na(old_score) || old_time <= k){
        score = score_fxn_alltime(model = model, k=k, series_to_eval = series_number,
                                       start_time = k, end_time = time, data=environments[[series_number]]$data,
                                       observations = environments[[series_number]]$y)
      }else if(old_time == time){
        score = old_score        
      }else{
        start_time = old_time+1
        end_time = time
        score = score_fxn_alltime(model = model, k=k, series_to_eval = series_number,
                                       start_time = start_time, end_time = end_time, data=environments[[series_number]]$data,
                                       observations = environments[[series_number]]$y)
        
        old_score = alpha^(time - old_time)*old_score # Assumes a constant alpha
        score = old_score + score
      }
      
    }
    return(score) #Switching to returning the score on the log-scale
  }
}

# Constraint function - all models must include a random walk intercept
constraint <- function(set){
  if(cset_contains_element(set,1) & length(set) <= maxmodel){
    TRUE
  }else{
    FALSE
  }
}

initialize_dlm <- function(model, series_number, delta=1, beta=1){
  
  p = length(model)
  initial_points = length(environments[[series_number]]$y_prior)
  
  # Initialize dlm parameters very loosely
  m = rep(0, p)
  C = diag(p)*10
  n = 1
  s = .1
  
  # Pull the prior data from global variables
  F = environments[[series_number]]$X_prior
  y = environments[[series_number]]$y_prior
  
  # Use forward filtering on prior data to set up dlm
  for(t in 1:initial_points){
    
    # Discount information
    n = delta * n
    C = C / beta
    
    # Find mean and variance at time t
    X = matrix(F[t, model])
    ft = crossprod(X, m)
    qt = crossprod(X, C) %*% X  + s
    qt = qt[1,1]

    # Posterior Update
    error = as.numeric(y[t] - ft)
    #At = C %*% as.matrix(F[t,model]) / qt
    At = C %*% X / qt
    rt = as.numeric((n + error^2 / qt) / (n + 1))
    
    n = n + 1
    s = s * rt
    m = m + At %*% error
    C = rt * (C - qt * tcrossprod(At, At))
  }
  
  return(list2env(list("model" = model,
                       "m" = m, 
                       "C" = C, 
                       "n" = n, 
                       "s" = s,
                       "time" = 0,
                       "likelihood" = NULL,
                       "preds" = NULL),
                  hash=TRUE))
  
}
# Function to evaluate the Kalman Filter for a given DLM
# The prior is theta0 ~ N(m0,  C0 v0/s0) and
# 1/v0 ~ Ga(n0/2, n0s0/2)
# delta = discount factor on sample size
# beta = discount factor on parameter variances
dynamic_linear_model <- function(time, F, y, model = NULL, series_number=NULL, save_preds=FALSE, dlm=NULL, delta = 1, beta = 1){
  
  if(is.null(dlm)){
    start_time = 1
    
    dlm = initialize_dlm(model, series_number, delta, beta)
    m = dlm$m
    C = dlm$C
    n = dlm$n
    s = dlm$s
    
    likelihood = NULL
    preds = NULL
    
  }else{
    start_time = dlm$time + 1
    
    m = dlm$m
    C = dlm$C
    n = dlm$n
    s = dlm$s
    
    likelihood = dlm$likelihood
    model = dlm$model
    preds = dlm$preds
  }
  
  if(time > 0 & time >= start_time){
    for(t in start_time:time){
      # Apply the discount factors
      
      # Option of a manual intervention to downweight models before a set time, to adapt quickly at the intervention point
      if(intervention == TRUE){
        if(t %in% interventions){
          n = 0.5 * n / delta
          C = C / 0.5 * beta  
        }
      }

      
      n = delta * n
      C = C / beta
      
      # likelihood of y
      X = matrix(F[t, model])
      # print(X)
      ft = crossprod(X, m)
      qt = crossprod(X, C) %*% X  + s
      qt = qt[1,1]
      likelihood = c(likelihood, dt((y[t] - ft)/sqrt(qt), n) / sqrt(qt))
      
      # Posterior Update
      error = as.numeric(y[t] - ft)
      #At = C %*% as.matrix(F[t,model]) / qt
      At = C %*% X / qt
      rt = as.numeric((n + error^2 / qt) / (n + 1))
      
      n = n + 1
      s = s * rt
      m = m + At %*% error
      C = rt * (C - qt * tcrossprod(At, At))
      
      if(save_preds){
        preds = c(preds, ft)
      }
      
    }
  }
  
  return(list2env(list("model" = model,
              "m" = m, 
              "C" = C, 
              "n" = n, 
              "s" = s,
              "time" = time,
              "likelihood" = likelihood,
              "preds" = preds),
              hash=TRUE))
}


### Extract the models into vector form from a string
extract_model_dlm <- function(model_string){
  start = str_locate(model_string, pattern = ": ")[2]
  mod = str_sub(model_string, start = start+1)
  mod = as.vector(sapply(str_split(mod, ", "), function(x) as.numeric(x)))
  return(sort(mod))
}

# naming system for dlms - model is a vector of predictor locations
dlm_name <- function(model){
  return(paste("Model:", toString(model)))  
}

#### Draw values from a given model
# Assume that theta, lambda are distributed Normal-Gamma
# Model = vector of column values pointing to the predictors used
# r, s = parameters so that lambda ~ Gamma(r/2, rs/2)
# a, R = parameters so that theta ~ N(a, R/(s*lambda)) 
# y = optional y value. Indicates that desired output is density, not random draws
# mean = flag to indicate that desired output is the mean of the t-distribution
# conditional = flag to indicate that we consider a model z to be conditioned on all other parameters being 0
model_dist_dlm <- function(time, F, dlm, y=NULL, mean=FALSE){
  X = F[time, dlm$model]
  model_dist(X = X, dlm = dlm, y = y, mean = mean)
}

model_dist <- function(X, dlm, y=NULL, mean=FALSE){
  ft = crossprod(X, dlm$m)
  qt = crossprod(X, dlm$C) %*% X  + dlm$s

  if(mean){
    return(ft)
  }
  
  if(is.null(y)){
    return(rt(1, dlm$n)*sqrt(qt) + ft)
  }else if(!(is.null(y))){
    return(dt((y - ft)/sqrt(qt), dlm$n)/sqrt(qt))
  }  
}

# Takes in a vector of probabilities on the log-scale. Renormalizes them by subtracting the max
# Would expect most/all of these log-probabilities to be highly negative, at bad times
rescale_probs <- function(probs){
  return(probs-max(probs))
}

## Given a vector of m KL-Divergences, want to use Newton's method to set p(z0) = 0.5
# Note that we assume the KL divergence of the true model is 0 (so exp(tau*KL(z0)) = 1)
# Which implies that the sum (p(z != z0)) = 1
# The true model should NOT be included in the list sent to this function
renormalize <- function(KL){
  # Initialize a value of tau
  tau = 1
  tolerance = 1E-3
  diff = 1
  # Newton-Raphson method
  while(diff > tolerance){
    func = 1 - sum(exp(-tau*KL))
    func_prime = sum(KL * exp(-tau*KL))
    tau = tau - func/func_prime
    diff = abs(func)
  }
  return(tau)
}

arrange_data <- function(series, max_lag = 12){
  lengths = sapply(series, length)
  lag_0 = sapply(series, function(x) length(x) > lengths[1])
  data = matrix(NA, nrow = lengths[1] - max_lag + 1, ncol = (max_lag * length(series) + sum(lag_0) + 1))
  
  # Now fill in the data frame
  data[,1] = 1 # set the intercept column
  
  # Fill in the lag-0 column, if any
  lag_0_cols = 0
  if(length(series) > 1){
    for(j in 2:length(series)){
      if(lag_0[j]){
        lag_0_cols = lag_0_cols + 1
        data[,1 + lag_0_cols] = series[[j]][(max_lag+1):lengths[j]]
        series[[j]] = series[[j]][-lengths[j]] #remove the same-time information now that it's been recorded
        lengths[j] = lengths[j] - 1
      }  
    }
  }
  
  
  # Fill in the regular columns
  for(j in 1:length(series)){
    for(i in 1:max_lag){
      data[ ,(i-1)*length(series) + j + lag_0_cols + 1] = series[[j]][(max_lag - i + 1):(lengths[j] - i + 1)]
    } 
  }
  
  return(data)
}

## Arrange_data_ddnm is a wrapper fxn
## This will automatically remove the most recent observations from all series with a lower number than the series we are arranging data for
## And then re-order the series into the appropriate ordering
## And then call the base arrange_data fxn

arrange_data_ddnm <- function(series, series_number, series_order, max_lag = 12){
  # Remove the most recent observations from all series with a lower number than the series we are arranging data for
  # so then we get the lag-0 information from all series with a higher number, but not the ones with the lower number
  for(i in 1:series_number){
    series[[i]] = series[[i]][-length(series[[i]])]
  }
  
  # Call arrange_data, with the proper re-ordering
  arrange_data(series, max_lag)
}

arrange_data_lags <- function(series, lags = c(1)){
  max_lag = max(lags)
  lengths = sapply(series, length)
  lag_0 = sapply(series, function(x) length(x) > lengths[1])
  data = matrix(NA, nrow = lengths[1] - max_lag + 1, ncol = (length(lags) * length(series) + sum(lag_0) + 1))
  
  # Now fill in the data frame
  data[,1] = 1 # set the intercept column
  c = 2
  
  # Fill in the lag-0 column, if any
  if(length(series) > 1){
    for(j in 2:length(series)){
      if(lag_0[j]){
        data[,c] = series[[j]][(max_lag+1):lengths[j]]
        series[[j]] = series[[j]][-lengths[j]] #remove the same-time information now that it's been recorded
        lengths[j] = lengths[j] - 1
        c = c + 1
      }  
    }
  }
  
  
  # Fill in the regular columns
  for(i in lags){
    for(j in 1:length(series)){
      data[ ,c] = series[[j]][(max_lag - i + 1):(lengths[j] - i + 1)]
      c = c + 1
    } 
  }
  
  return(data)
}

arrange_data_ddnm_lags <- function(series, series_number, lags = c(-1)){
  # Remove the most recent observations from all series with a lower number than the series we are arranging data for
  # so then we get the lag-0 information from all series with a higher number, but not the ones with the lower number
  for(i in 1:series_number){
    series[[i]] = series[[i]][-length(series[[i]])]
  }
  
  # Call arrange_data, with the proper re-ordering
  arrange_data_lags(series, lags)
}


## Quick function to define priors for my model using the first initial_points points of the series
# environment needs to have the basics already defined
encompassing_model_prior <- function(data, y, environment){
  #series_number = environment$series_number
  # data = environment$X_full[1:initial_points,]
  # y = environment$y_full[1:initial_points]
  
  model_data = data.frame("y" = y, data)
  # Hard coding in a discount factor here
  # environment$mod = lm(y ~ . - 1, data = model_data, weights = 0.98^((nrow(data)-1):0))
  environment$mod = lm(y ~ . - 1, data = model_data)
  environment$m0 = as.matrix(coef(environment$mod))
  environment$C0 = 2*summary(environment$mod)$cov.unscaled * summary(environment$mod)$sigma^2
  environment$s0 = 5*summary(environment$mod)$sigma^2
  environment$n0 = initial_points/2
  
  return(environment)
}
