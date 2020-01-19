#Startup Script


# Loading the data
full_dataset = read.csv(paste("./Initial Data/", csv_name, ".csv", sep=""))
nvars = ncol(full_dataset) - 1

# Pulling the X, Y values we're using from the dataset,
initial_points = 30
n = nrow(full_dataset)
X_full = cbind(1, as.matrix(full_dataset[-(1:initial_points), 1:nvars]))
y_full = full_dataset[-(1:initial_points), nvars+1]

# Adding in random dates, to give something to print out while being run
Date = read.csv("./Initial Data/Data 1965-2016.csv")$Date
full_dataset$Date = Date[(length(Date) - n + 1):length(Date)]

# Setting up all necessary variables
series_names = "y"
num_series = length(series_names)
environments = lapply(series_names, function(x) x = new.env(hash=TRUE))
names(environments) = series_names
leaderboard = data.frame()
dates = full_dataset %>% slice((initial_points + 1):n) %>% mutate(Date = as.character(Date)) %>% pull(Date)
tdates = sapply(dates, function(x) substr(x, nchar(x)-3, nchar(x)), USE.NAMES=FALSE)[seq(1, length(dates), by=24)]
tticks = 1:length(dates)
tticks = tticks[seq(1, length(dates), by=24)]
Rdates = as.Date(dates, "%m/%d/%Y") #Don't need this, but maybe at some point will be useful?
xnames = c("Intercept", names(full_dataset)[1:nvars])

l = nrow(X_full)
tticks = seq(1, l, by = 20)
tdates = seq(1, l, by = 20)

## Define one environment / namespace per series in the DDNM

## Setting all the information and priors that are unique to each environment

# Select the prior data
X_prior = cbind(1, as.matrix(full_dataset[1:initial_points, 1:nvars]))
y_prior = full_dataset[1:initial_points, nvars+1]


model = c(1, 2, 3)
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

