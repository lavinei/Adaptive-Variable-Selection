## LPFDS score
lpfds_analytic <- function(dlm, k = 12, series_to_eval, time, data, y){
  cov = matrix(NA, nrow=k, ncol=k)
  mu = matrix(NA, nrow = k, ncol=1)
  X = matrix(data[(time-k+1):time, dlm$model], ncol=length(dlm$model))
  dlm$W = (1-delta)/delta * dlm$C

  for(i in 1:k){
    # R(i) is the variance of the coefficients at time t+i
    R_i = dlm$C + i * dlm$W
    
    # Mean and variance for observation y_(t+i)
    mu[i] = crossprod(X[i,], dlm$m)
    cov[i,i] = crossprod(X[i,], R_i) %*% X[i,]  + dlm$s
    
    # Covariance between y_(t+i) and y_(t+i+1)...y_(t+k)
    if(i < k){
      for(j in (i+1):k){
        cov[i, j] = cov[j, i] = crossprod(X[i,], R_i) %*% X[j,]
      }
    }
    
  }
  
  # Approximate degrees of freedom
  dof = beta^k * dlm$n
  
  # Turn covariance into the scale matrix
  sigma = cov * (dof-2)/dof
  
  # Log density of the k-step ahead path
  return(dmvt(y[(time-k+1):time], mu,  sigma, dof, log = TRUE))
}


t_plus_k_density <- function(dlm, k = 12, series_to_eval, time, data, y){
  
  # X = matrix(data[(time-k+1):time, dlm$model], ncol=length(dlm$model))
  dlm$W = (1-delta)/delta * dlm$C
  X = data[time, dlm$model]

  # Calculating the mean and variance of the forecast for time t+k
  mu_k = crossprod(X, dlm$m)
  R_k = dlm$C + k * dlm$W
  var_k = crossprod(X, R_k) %*% X + dlm$s

  # Approximate degrees of freedom
  dof = beta^k * dlm$n
  
  # Log density of the forecast for time t+k
  log(dt((y[time] - mu_k)/sqrt(var_k), dof)/sqrt(var_k))
}

t_plus_k_density_ar <- function(dlm, k = 12, series_to_eval, time, data, y){
  
  index = seq(num_series + 1, by = num_series, length = ifelse(choose_lags, length(lags), k))
  var = vector("numeric", length = k)
  mu = vector("numeric", length = k)
  dof = vector("numeric", length = k)
  sim_pts = vector("numeric", length = k)
  X = matrix(data[(time-k+1):time, dlm$model], ncol=length(dlm$model))
  dlm$W = (1-delta)/delta * dlm$C
  
  for(i in 1:k){
    # R(i) is the variance of the coefficients at time t+i
    R_i = dlm$C + i * dlm$W
    
    # Mean and variance for observation y_(t+i)
    X = data[time - k + i, dlm$model]
    mu[i] = crossprod(X, dlm$m)
    var[i] = crossprod(X, R_i) %*% X  + dlm$s
    dof[i] = beta^i * dlm$n
    
    sim_pts[i] = mu[i] + sqrt(var[i]) * rt(1, df = dof[i])
    if(i < k){
      
      if(choose_lags){
        ls = lags[lags<=i]
        data[time - k + i + 1, index[1:length(ls)]] = sim_pts[i - ls + 1]
      }else{
        data[time - k + i + 1, index[1:i]] = sim_pts[i:1]
      }
      
      y[time - k + i] = sim_pts[i]
      
    }
    
  }
  
  # Log density of the k-step ahead forecast
  return(log(dt((y[time] - mu[k])/sqrt(var[k]), dof[k])/sqrt(var[k])))
}

t_plus_k_density_ar_wrapper <- function(dlm, k = 12, series_to_eval, time, data, y, nsamps = 50){
  
  # Check if the model actually has an AR component, before resorting to the ar function, which is slower
  index = seq(num_series + 1, by = num_series, length = ifelse(choose_lags, length(lags), k))
  if(any(dlm$model %in% index)){
    mean(vapply(1:nsamps, function(z) t_plus_k_density_ar(dlm, k = k, series_to_eval, time, data, y), numeric(1)))
  }else{
    t_plus_k_density(dlm, k = k, series_to_eval, time, data, y)
  }
  
}

# This version is the rao-blackwellized version of the k-step forecast density, for speed
t_plus_kvec_density_ar_generator <- function(kvec){
  function(dlm, k = 12, series_to_eval, time, data, y){
    
    index = seq(num_series + 1, by = num_series, length = k)
    var = vector("numeric", length = k)
    mu = vector("numeric", length = k)
    dof = vector("numeric", length = k)
    sim_pts = vector("numeric", length = k)
    # X = matrix(data[(time-k+1):time, dlm$model], ncol=length(dlm$model))
    dlm$W = (1-delta)/delta * dlm$C
    
    vec = kvec[kvec <= k]
    if(length(vec) == 0){
      return(1)
    }
    
    for(i in 1:k){
      # R(i) is the variance of the coefficients at time t+i
      R_i = dlm$C + i * dlm$W
      
      # Mean and variance for observation y_(t+i)
      X = data[time - k + i, dlm$model]
      mu[i] = crossprod(X, dlm$m)
      var[i] = crossprod(X, R_i) %*% X  + dlm$s
      dof[i] = beta^i * dlm$n
      
      if(i %in% vec){
        # Plug in the true value into the design matrix
        sim_pts[i] = y[time - k + i]
      }else{
        # Otherwise plug in the simulated value
        sim_pts[i] = mu[i] + sqrt(var[i]) * rt(1, df = dof[i])
        y[time - k + i] = sim_pts[i]
      }
      
      if(i < k){
        data[time - k + i + 1, index[1:i]] = sim_pts[i:1]
      }
      
      
    }
    
    # Calculate the covariance between points in the kvec
    l = length(vec)
    cov = matrix(NA, nrow=l, ncol=l)
    diag(cov) = var[vec]
    if(l > 1){
      for(i in 1:(l-1)){
        for(j in (i+1):l){
          R_i = dlm$C + vec[i] * dlm$W
          cov[i, j] = cov[j, i] = crossprod(data[time - k + vec[i], dlm$model], R_i) %*% data[time - k + vec[j], dlm$model]
        }
      }
    }
    
    
    # Approximate degrees of freedom
    dof = dof[vec[l]]
    
    # Turn covariance into the scale matrix
    sigma = cov * (dof-2)/dof
    
    # Log density of the kvec-steps ahead path
    return(dmvt(y[time - k + vec], mu[vec],  as.matrix(sigma), dof, log = TRUE))
  }
}

t_plus_kvec_density_ar_generator_wrapper <- function(kvec){
  t_plus_kvec_density_ar = t_plus_kvec_density_ar_generator(kvec)
  
  function(dlm, k = 12, series_to_eval, time, data, y, nsamps = 50){
    mean(vapply(1:nsamps, function(z) t_plus_kvec_density_ar(dlm, k = k, series_to_eval, time, data, y), numeric(1)))
  }  
}


t_plus_kvec_density_generator <- function(kvec){
  function(dlm, k = 12, series_to_eval, time, data, y){
    cov = matrix(NA, nrow=k, ncol=k)
    mu = matrix(NA, nrow = k, ncol=1)
    X = matrix(data[(time-k+1):time, dlm$model], ncol=length(dlm$model))
    dlm$W = (1-delta)/delta * dlm$C
    
    for(i in 1:k){
      # R(i) is the variance of the coefficients at time t+i
      R_i = dlm$C + i * dlm$W
      
      # Mean and variance for observation y_(t+i)
      mu[i] = crossprod(X[i,], dlm$m)
      cov[i,i] = crossprod(X[i,], R_i) %*% X[i,]  + dlm$s
      
      # Covariance between y_(t+i) and y_(t+i+1)...y_(t+k)
      if(i < k){
        for(j in (i+1):k){
          cov[i, j] = cov[j, i] = crossprod(X[i,], R_i) %*% X[j,]
        }
      }
      
    }
    
    # Approximate degrees of freedom
    dof = beta^k * dlm$n
    
    # Turn covariance into the scale matrix
    sigma = cov * (dof-2)/dof
    
    # Log density of the k-step ahead path
    vec = kvec[kvec <= k]
    if(length(vec) == 0){
      return(1)
    }else{
      return(dmvt(y[time - k + vec], mu[vec],  as.matrix(sigma[vec, vec]), dof, log = TRUE))
    }
    
  }
}

score_alltime_generator <- function(score_fxn, nsamps = NA){
  if(is.na(nsamps)){
    function(model, k = 12, series_to_eval, start_time, end_time, data, observations){
      weights = alpha^((end_time - start_time):0)
      scores = NULL
      dlm = dynamic_linear_model(time = start_time - k, F = data, y = observations, model = model, series_number = series_to_eval, delta=delta, beta=beta)
      for(time in start_time:end_time){
        dlm = dynamic_linear_model(time = time-k, F = data, y = observations, model = model, series_number = series_to_eval, dlm=dlm, delta=delta, beta=beta)
        scores = c(scores, score_fxn(dlm, k, series_to_eval, time, data, observations))
      }
      
      return(sum(weights * scores))
    }
  }else{
    function(model, k = 12, series_to_eval, start_time, end_time, data, observations){
      weights = alpha^((end_time - start_time):0)
      scores = NULL
      dlm = dynamic_linear_model(time = start_time - k, F = data, y = observations, model = model, series_number = series_to_eval, delta=delta, beta=beta)
      for(time in start_time:end_time){
        dlm = dynamic_linear_model(time = time-k, F = data, y = observations, model = model, series_number = series_to_eval, dlm=dlm, delta=delta, beta=beta)
        scores = c(scores, score_fxn(dlm, k, series_to_eval, time, data, observations, round(nsamps * weights[time - start_time + 1], 0)))
      }
      
      return(sum(weights * scores))
    }
  }
  
}
