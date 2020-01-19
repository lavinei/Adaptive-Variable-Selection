#Startup Script

# Loading the data
full_dataset = read.csv("./Initial Data/Data 1965-2016.csv")

# Pulling the X, Y values we're using from the dataset,
first_observation = which(full_dataset$Date == start_date)
series_full = as.list(full_dataset %>% slice(first_observation:nrow(full_dataset)) %>% select(series_names))
series = as.list(full_dataset %>% slice((first_observation-max_lag):(first_observation-1)) %>% select(series_names))

if(choose_lags){
  X_full = arrange_data_ddnm_lags(as.list(full_dataset %>% slice((first_observation-max_lag):nrow(full_dataset)) %>% select(series_names)), 1, lags = lags)
  y_full = full_dataset %>% slice(first_observation:nrow(full_dataset)) %>% pull(series_names[1])
}else{
  X_full = arrange_data_ddnm(as.list(full_dataset %>% slice((first_observation-max_lag):nrow(full_dataset)) %>% select(series_names)), 1, series_order = 1:length(series_names), max_lag = max_lag)
  y_full = full_dataset %>% slice(first_observation:nrow(full_dataset)) %>% pull(series_names[1])
}

# Setting up all necessary variables
names = series_names
num_series = length(series_names)
environments = lapply(series_names, function(x) x = new.env(hash=TRUE))
names(environments) = series_names
leaderboard = data.frame()
dates = full_dataset %>% slice(first_observation:nrow(full_dataset)) %>% mutate(Date = as.character(Date)) %>% pull(Date)
tdates = sapply(dates, function(x) substr(x, nchar(x)-3, nchar(x)), USE.NAMES=FALSE)[seq(1, length(dates), by=24)]
tticks = 1:length(dates)
tticks = tticks[seq(1, length(dates), by=24)]
Rdates = as.Date(dates, "%m/%d/%Y")
xnames = NULL
xnames = c(xnames, "Intercept")
for(i in 2:num_series){
  xnames = c(xnames, paste("lag 0", series_names[i]))
}
if(choose_lags){
  for(i in lags){
    for(j in 1:num_series){
      xnames = c(xnames, paste("lag", i, series_names[j]))  
    }
  }
  
}else{
  for(i in 1:max_lag){
    for(j in 1:num_series){
      xnames = c(xnames, paste("lag", i, series_names[j]))  
    }
  }
}

## Define one environment / namespace per series in the DDNM

## Setting all the information and priors that are unique to each environment
initial_points = 60

# Select the prior data
if(choose_lags){
  X_prior = arrange_data_ddnm_lags(as.list(full_dataset %>% slice((first_observation-max_lag-initial_points):(first_observation-1)) %>% select(series_names)), 1, lags = lags)
  y_prior = full_dataset %>% slice((first_observation-initial_points):(first_observation-1)) %>% pull(series_names[1])
}else{
  X_prior = arrange_data_ddnm(as.list(full_dataset %>% slice((first_observation-max_lag-initial_points):(first_observation-1)) %>% select(series_names)), 1, series_order = 1:length(series_names), max_lag = max_lag)
  y_prior = full_dataset %>% slice((first_observation-initial_points):(first_observation-1)) %>% pull(series_names[1])
}

model = c(1, 2, 3, 4)
time = 0


for(i in 1:num_series){
  if(i == 1){
    environments[[i]]$X_full = X_full
    environments[[i]]$y_full = as.numeric(y_full)
    
    environments[[i]]$y = NULL
    environments[[i]]$series_order = 1:num_series
    environments[[i]]$Omega = cset(as.numeric(1:ncol(X_full)))
    environments[[i]]$series_number = i
    environments[[i]]$xnames = xnames
    environments[[i]]$data = NULL
    environments[[i]]$X_prior = X_prior
    environments[[i]]$y_prior = y_prior
    
    # Defining storage variables
    environments[[i]]$model_list = list()
    environments[[i]]$preds = NULL
    environments[[i]]$cred_lower = NULL
    environments[[i]]$cred_upper = NULL
    environments[[i]]$cred_25 = NULL
    environments[[i]]$cred_75 = NULL
    environments[[i]]$y_prob_store = NULL
    environments[[i]]$runtime_store = NULL
    environments[[i]]$model_store = list()
    environments[[i]]$predictor_inclusion_probabilities = NULL
    environments[[i]]$predictor_inclusion_probabilities_store = NULL
    
    # Defining running variables to track models
    environments[[i]]$model_list = new.env(hash=TRUE)
    environments[[i]]$leaderboard = data_frame()
    
    # Defining the inital model's dlm 
    environments[[i]]$model = model
    environments[[i]]$dlm = dynamic_linear_model(time, F = NULL, y = NULL, model = model, series_number = i)
    
  }else{
    environments[[i]]$X_full = X_full[,-(2:i)]
    environments[[i]]$y_full = X_full[,i]
    environments[[i]]$X_prior = X_prior[,-(2:i)]
    environments[[i]]$y_prior = X_prior[,i]
    
    environments[[i]]$y = NULL
    environments[[i]]$series_order = c(i:num_series, 1:(i-1))
    environments[[i]]$Omega = cset(as.numeric(1:ncol(environments[[i]]$X_full)))
    environments[[i]]$xnames = xnames[-(2:i)]
    environments[[i]]$series_number = i
    if(i == num_series){
      environments[[i]]$data = t(matrix(environments[[i]]$X_full[1,]))
    }else if(i < num_series){
      environments[[i]]$data = NULL
    }
    
    # Defining storage variables
    environments[[i]]$model_list = list()
    environments[[i]]$preds = NULL
    environments[[i]]$cred_lower = NULL
    environments[[i]]$cred_upper = NULL
    environments[[i]]$cred_25 = NULL
    environments[[i]]$cred_75 = NULL
    environments[[i]]$y_prob_store = NULL
    environments[[i]]$runtime_store = NULL
    environments[[i]]$model_store = list()
    environments[[i]]$predictor_inclusion_probabilities = NULL
    environments[[i]]$predictor_inclusion_probabilities_store = NULL
    
    # Defining running variables to track models
    environments[[i]]$model_list = new.env(hash=TRUE)
    environments[[i]]$leaderboard = data_frame()
    
    # Defining the inital model's dlm 
    environments[[i]]$model = model
    environments[[i]]$dlm = dynamic_linear_model(time, F = NULL, y = NULL, model = model, series_number = i)
    
  }
}
