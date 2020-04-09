# Comparison plot functions

# Plot long term inclusion probability over time - actually would need to save this during the simulation, sadly
long_term_preds <- function(avs_envs, bma_envs, series_number){
  avs = avs_envs[[series_number]]$predictor_inclusion_probabilities_store
  bma = bma_envs[[series_number]]$predictor_inclusion_probabilities_store
  
  # Only get lags greater than 7 months - except can't get what I want from the marginal inclusion probabilities.
  
  plot_data = data.frame("time" = 1:time, avs, bma)
  p = ggplot()
}

LPDR_logitudinal <- function(kdens_data){
  p = ggplot(kdens_data, aes(x = time, y = k_step_path_density, color=model, linetype=model))
  p = p + geom_line(size=1.5)
  return(standard_theme(p, ylab = "log k-step Path Forecast Density"))
}

LPDR_logitudinal_relative <- function(kdens_data){
  models = unique(kdens_data$model)
  avs_idx = which(models == "avs")
  tmp = kdens_data %>%
    spread(key = model, value = k_step_path_density)
  avs_col = which(names(tmp) == "avs")
  model_cols = which(names(tmp) %in% models[-avs_idx])
  kdens_data = as.data.frame(lapply(model_cols, function(c) tmp[,c] - tmp[,avs_col]))
  kdens_data = cbind(time = tmp$time, kdens_data)

  kdens_data = kdens_data %>%
    gather("model", "k_step_path_density", seq(2, length(models)))
  p = ggplot(kdens_data, aes(x = time, y = k_step_path_density, color=model))
  p = p +
    geom_line(size=1) +
    geom_hline(yintercept = 0)
  return(standard_theme(p, ylab = "LPFD Score Relative to AVS"))
}

LPDR <- function(kdens_data){
  p = ggplot(kdens_data, aes(x = time, y = LPDR, color=model)) +
    geom_line() #+
    #geom_hline(yintercept = 0)
  standard_theme(p) +
    scale_y_continuous(name = "Cumulative LPFDS relative to AVS")
}

# This function isn't completely written yet
LPDR_window <-function(kdens_window_data, comps = 2, window = 48){
  kdens_window_data = kdens_window_data %>% mutate(
    tvvar = kdens_tvvar - kdens_avs,
    bma = kdens_bma - kdens_avs) %>%
    transmute(time = time)
  p = ggplot(kdens_data, aes(x = time, y = LPDR, color=model))
  p = p + geom_line() +
    theme_bw() +
    scale_x_continuous(name = "Time", breaks=tticks, labels=tdates) +
    ggtitle("Log k-step Path Density Ratio")
  
  return(p)
}

variable_inclusion_plot <- function(predictor_inclusion_probabilities, xnames, training_period, simulation=FALSE){
  pip = tbl_df(data.frame(round(predictor_inclusion_probabilities, 2)))
  vars = ncol(pip)
  names(pip) = 1:vars
  pip = pip %>% mutate("time" = row_number())
  pip_plot = pip %>% gather("variable", "inclusion_probability", 1:(ncol(pip)-1))
  pip_plot$variable = as.numeric(pip_plot$variable)
  pip_plot = pip_plot %>% filter(time > training_period)
  if(simulation){
    pip_plot = pip_plot %>% mutate(time = time - training_period)
  }
  
  p = ggplot(pip_plot, aes(x = time, y = variable)) + 
    geom_tile(aes(fill = inclusion_probability), color="black") 
  p = standard_theme(p)
  p + scale_y_continuous(name = "Covariates", breaks=1:length(xnames), labels = xnames) +
    scale_fill_continuous(name = "Inclusion Probability", limits= c(0, 1))
}

model_plot <- function(model_store, xnames, Omega, training_period, recession_bar = TRUE){
  tmp = sapply(model_store, function(model) Omega %in% model)
  plot_models = data.frame(tmp) %>% mutate("theta" = 1:nrow(tmp))
  names(plot_models)[1:(time+1)] = as.character(0:time)
  plot_models = gather(plot_models, "time", "in_model", 1:length(model_store)) %>%
    mutate(time = as.numeric(time)) %>%
    filter(time > training_period)
  names(plot_models)[1] = "theta"
  p = ggplot(plot_models, aes(x = as.numeric(time), y = theta))
  model_plot = p + geom_tile(aes(fill = in_model), color="black")
  model_plot = standard_theme(model_plot)
  
  if(recession_bar){
   model_plot = model_plot +
     geom_vline(xintercept=tticks[tdates == "2007"] + 11) + # December 2007 
     geom_vline(xintercept=tticks[tdates == "2009"] + 5) # end of June 2009
  }
  
  model_plot + 
    scale_y_continuous(name = "Covariates", breaks=1:length(xnames), labels = xnames) +
    scale_fill_discrete(name = "In Model") +
    theme(axis.title = element_text(size=16, face="bold"), #14
          axis.text.x = element_text(angle = 0, hjust = 1, size=13), #11
          axis.text.y = element_text(angle = 0, size=10), #11
          legend.position="none",
          plot.title = element_text(size=20),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(angle = 0, size = 14, face="bold"), #14
          strip.text.y = element_text(angle = -90, size = 14, face="bold"))
  

}

model_plot_sim <- function(model_store, xnames, Omega, training_period){
  tmp = sapply(model_store, function(model) Omega %in% model)
  plot_models = data.frame(tmp) %>% mutate("theta" = 1:nrow(tmp))
  names(plot_models)[1:(time+1)] = as.character(0:time)
  plot_models = gather(plot_models, "time", "in_model", 1:length(model_store)) %>%
    mutate(time = as.numeric(time)) %>%
    filter(time > training_period) %>%
    mutate(time = time - training_period)
  names(plot_models)[1] = "theta"
  p = ggplot(plot_models, aes(x = as.numeric(time), y = theta))
  model_plot = p + geom_tile(aes(fill = in_model), color="black")
  model_plot = standard_theme(model_plot)
  
  model_plot + 
    scale_y_continuous(name = "Covariates", breaks=1:length(xnames), labels = xnames) +
    theme(axis.title = element_text(size=16, face="bold"), #14
          axis.text.x = element_text(angle = 0, hjust = 1, size=13), #11
          axis.text.y = element_text(angle = 0, size=13), #11
          legend.text = element_text(size=14), #12
          legend.title = element_text(size=16),
          plot.title = element_text(size=20),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(angle = 0, size = 14, face="bold"), #14
          strip.text.y = element_text(angle = -90, size = 14, face="bold")) +
    scale_fill_discrete(name = "In Model")
}

perc_error_marg <- function(plot_data, series_number){
  dat = plot_data %>% filter(series == series_number & !(is.na(k))) %>% group_by(k, model) %>%
    transmute(percerror = 100*abs(prediction - true_value)/true_value) %>%
    summarize(Percent_Error = mean(percerror))

  # Comparison plot over time, for a specific number of steps ahead
  p = ggplot(dat, aes(x = k, y = Percent_Error, color=model))
  p = p + geom_line(size=1.5) + 
    geom_point(size=1.5) +
    theme_bw() +
    scale_x_continuous(breaks=1:k) +
    ggtitle(paste(names[series_number], ": Average Marginal % Error", sep="")) +
    theme(axis.title=element_text(size=14, face="bold"),
          axis.text.x = element_text(hjust = 1, size=10),
          axis.text.y = element_text(size=10),
          legend.text = element_text(size=10),
          plot.title = element_text(size=20))
  return(p)
}

perc_error_marg_combined <- function(plot_data, series_names){
  dat = plot_data %>% filter(!(is.na(k))) %>% group_by(k, model, series) %>%
    transmute(percerror = 100*abs(prediction - true_value)/true_value) %>%
    summarize(Percent_Error = mean(percerror))
  
  dat = dat %>%
    mutate(series = factor(series, labels = series_names))

  # Comparison plot over time, for a specific number of steps ahead
  p = ggplot(dat, aes(x = k, y = Percent_Error, color=model))
  p = p + geom_line(size=1.5) + 
    geom_point(size=1.5) +
    facet_wrap(~series, nrow=1)
  
  standard_theme(p, ylab="Percent Forecast Error") +
    scale_x_continuous(name = "Forecast Length (Months)", breaks=seq(1,k, by=2))
  
}

rMSE_marg <- function(plot_data, series_number){
  dat = plot_data %>% filter(series == series_number & !(is.na(k))) %>% group_by(k, model) %>%
    transmute(sqerror = (prediction - true_value)^2) %>%
    summarize(rMSFE = sqrt(mean(sqerror)))

  # Comparison plot over time, for a specific number of steps ahead
  p = ggplot(dat, aes(x = k, y = rMSFE, color=model))
  p = p + geom_line(size=1.5) + 
    geom_point(size=1.5)
  
  standard_theme(p, ylab = paste("rMSFE:", names[series_number], sep=" ")) +
    scale_x_continuous(name = "Forecast Length (Months)", breaks=seq(1,k, by=2))
}

rMSE_marg_combined <- function(plot_data, series_names, scales = "fixed"){
  dat = plot_data %>% filter(!(is.na(k))) %>% group_by(k, model, series) %>%
    transmute(sqerror = (prediction - true_value)^2) %>%
    summarize(rMSFE = sqrt(mean(sqerror)))
  
  dat = dat %>%
    mutate(series = factor(series, labels = series_names))

  # Comparison plot over time, for a specific number of steps ahead
  p = ggplot(dat, aes(x = k, y = rMSFE, color=model, linetype=model))
  p = p + geom_line(size=1.5) + 
     geom_point(size=2) +
    # scale_colour_discrete("")+
    # scale_linetype_manual("", guide=FALSE)+
    if(length(series_names) > 3){
      facet_wrap(~series, nrow=2, scales = scales)
    }else{
      facet_wrap(~series, nrow=1, scales = scales)
    }
     
  
  standard_theme(p, ylab = "rMSFE") +
    scale_x_continuous(name = "Forecast Length (Months)", breaks=seq(1,k, by=2)) +
    scale_color_brewer(name = "", type="qual", palette=6) #+
    # labs(color = "", linetype="")
}

forecast_combined <- function(plot_data, series_number, forecast_length){
  ylims = range(plot_data %>% filter(series == series_number) %>% pull(true_value))
  mean = mean(plot_data %>% filter(series == series_number) %>% pull(true_value))
  ylims[1] = ylims[1] - (mean - ylims[1])
  ylims[2] = ylims[2] - (mean - ylims[2])
  
  plot_data = plot_data %>%
    filter(series == series_number & k == forecast_length)
  
  p_preds = ggplot(plot_data, aes(x = time, y = prediction))
  predictive_plot = p_preds + 
    geom_ribbon(aes(ymax = upper, ymin=lower), fill = "grey90") +
    geom_ribbon(aes(ymax = cred_75, ymin=cred_25), fill="grey70") +
    geom_line(aes(y = true_value, color = "Data")) +
    geom_line(aes(color = "Predicted Value")) +
    #  geom_vline(xintercept = input$time) +
    ggtitle(paste(forecast_length, "-step Prediction vs Data"))  +
    ylab("Prediction") +
    theme_bw() +
    geom_vline(xintercept = (interventions - 1 + forecast_length), color="black") +
    theme(plot.title = element_text(size=22), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.text = element_text(size=11),
          axis.text.x = element_text(angle=90)) + 
    scale_color_brewer(type="qual", palette=6) +
    scale_x_continuous(name = "Time", breaks=tticks, labels=tdates) +
    coord_cartesian(ylim = ylims) +
    facet_wrap(~model, ncol=1, scales="free_y")
  
  return(predictive_plot)
}

forecast_combined_nocredint <- function(plot_data, forecast_length){
  ylims = range(plot_data %>% pull(true_value))
  mean = mean(plot_data %>% pull(true_value))
  ylims[1] = ylims[1] - (mean - ylims[1])
  ylims[2] = ylims[2] - (mean - ylims[2])
  
  plot_data = plot_data %>%
    filter(k == forecast_length)
  
  p_preds = ggplot(plot_data, aes(x = time, y = prediction, color=model))
  predictive_plot = p_preds + 
    geom_line() +
    #geom_ribbon(aes(ymax = upper, ymin=lower), fill = "grey90") +
    #geom_ribbon(aes(ymax = cred_75, ymin=cred_25), fill="grey70") +
    geom_line(aes(y = true_value, color = "Data"), color="black") +
    #  geom_vline(xintercept = input$time) +
    ggtitle(paste(forecast_length, "-step Prediction vs Data"))  +
    ylab("Prediction") +
    theme_bw() +
    geom_vline(xintercept = (interventions - 1 + forecast_length), color="black") +
    theme(plot.title = element_text(size=22), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.text = element_text(size=11),
          axis.text.x = element_text(angle=90)) + 
    scale_color_brewer(type="qual", palette=6) +
    scale_x_continuous(name = "Time", breaks=tticks, labels=tdates) +
    #coord_cartesian(ylim = ylims) +
    facet_wrap(~series, ncol=1, scales="free_y")
  
  return(predictive_plot)
}

forecast_multiple_lengths <- function(plot_data, series_number, series_names, mod, forecast_lengths){
  ylims = range(plot_data %>% filter(series == series_number) %>% pull(true_value))
  mean = mean(plot_data %>% filter(series == series_number) %>% pull(true_value))
  ylims[1] = ylims[1] - (mean - ylims[1])
  ylims[2] = ylims[2] - (mean - ylims[2])
  
  forecast_length_names = paste(forecast_lengths, "Months Ahead", sep=" ")
  
  plot_data = plot_data %>%
    filter(series == series_number & k %in% forecast_lengths & model == mod) %>%
    mutate(forecast_length_name = factor(k, levels = forecast_lengths, labels = forecast_length_names))
  
  
  p_preds = ggplot(plot_data, aes(x = time, y = prediction))
  predictive_plot = p_preds + 
    geom_ribbon(aes(ymax = upper, ymin=lower), fill = "grey90") +
    geom_ribbon(aes(ymax = cred_75, ymin=cred_25), fill="grey70") +
    geom_line(aes(y = prediction, color = "a", linetype='a'), size=1) +
    geom_line(aes(y = true_value, color = "b", linetype='b'), size=1) +
    geom_vline(xintercept = (interventions - 1 + forecast_lengths), color="black") +
    coord_cartesian(ylim = ylims) +
    facet_wrap(~forecast_length_name, ncol=1, strip.position = "left")

    return(standard_theme(predictive_plot, ylab = NULL))
  
    #ylab = paste(series_names[series_number], "Prediction", sep=" "))
           
    # theme_bw() +
    # theme(plot.title = element_text(size=22), panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(), axis.text = element_text(size=11),
    #       axis.text.x = element_text(angle=90)) + 
    # scale_color_brewer(type="qual", palette=6) +
    # scale_x_continuous(name = "Time", breaks=tticks, labels=tdates) +
    
    
  
  # return(predictive_plot)
}

forecast_multiple_mods <- function(plot_data, series_number, series_names, mods, forecast_length){
  ylims = range(plot_data %>% filter(series == series_number) %>% pull(true_value))
  mean = mean(plot_data %>% filter(series == series_number) %>% pull(true_value))
  ylims[1] = ylims[1] - (mean - ylims[1])
  ylims[2] = ylims[2] - (mean - ylims[2])
  
  # forecast_mod_names = paste(str_to_upper(mod), "Forecast", forecast_length, "Months Ahead", sep=" ")
  forecast_mod_names = map(mods, function(mod) paste(str_to_upper(mod), "Forecast", sep=" "))
  
  plot_data = plot_data %>%
    filter(series == series_number & k == forecast_length & model %in% mods) %>%
    mutate(forecast_mod_name = factor(model, levels = mods, labels = forecast_mod_names))
  
  
  p_preds = ggplot(plot_data, aes(x = time, y = prediction))
  predictive_plot = p_preds + 
    geom_ribbon(aes(ymax = upper, ymin=lower), fill = "grey90") +
    geom_ribbon(aes(ymax = cred_75, ymin=cred_25), fill="grey70") +
    geom_line(aes(y = true_value, color = "Data")) +
    geom_line(aes(color = "Predicted Value")) +
    coord_cartesian(ylim = ylims) +
    facet_wrap(~forecast_mod_name, ncol=1, strip.position = "left")
  
  return(standard_theme(predictive_plot, ylab = NULL))
}


error_time_smooth <- function(plot_data, series_number, months_ahead, smoothed_months){
  dat = plot_data %>% filter(series == series_number & !is.na(k) & k == months_ahead) %>%
    group_by(time, model) %>%
    transmute(rMSE = sqrt((prediction - true_value)^2)) %>%
    arrange(model, time) 
  names(dat) = c("Time", "model", "rMSE")
  dat = dat %>% group_by(model) %>% mutate(rollrMSE = rollmean(x=rMSE, k=smoothed_months, align="center", fill=NA, na.rm=T))
  
  # Command to add in a rolling window:
  # rollmean(dat$rMSE, 5, fill=NA, na.rm=T)
  
  p = ggplot(dat, aes(x = Time, y=rollrMSE, color=model))
  p = p +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_x_continuous(name = "Time", breaks=tticks, labels=tdates) +
    ggtitle(paste(months_ahead, "-Step Forecast rMSE:", smoothed_months, " Month Rolling Average", sep=""))
  
  return(p)
}

error_time <- function(plot_data, series_number, months_ahead){
  dat = plot_data %>% filter(series == series_number & !is.na(k) & k == months_ahead) %>%
    group_by(time, model) %>%
    transmute(rMSE = sqrt((prediction - true_value)^2)) %>%
    arrange(model, time) 
  names(dat) = c("Time", "model", "rMSE")
  
  # Command to add in a rolling window:
  # rollmean(dat$rMSE, 5, fill=NA, na.rm=T)
  
  p = ggplot(dat, aes(x = Time, y=rMSE, color=model))
  p = p +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_x_continuous(name = "Time", breaks=tticks, labels=tdates) +
    ggtitle(paste(months_ahead, "-Step Forecast rMSE", sep=""))
  
  return(p)
}

get_plot_data <- function(cred_int_store, forecast_sample_store, diff, training_period){
  k = forecast_length # hack, b/c the plotting uses k instead of forecast length; often these are the same
  
  if(diff){ # If the predictions were done on differences between time steps, instead of on the absolute values
    ##### DIFF POST-PROCESSING OF CRED_INT_STORE
    timesteps = time
    full_dataset = read.csv("~/Homework/SAMSI Bays Opt/Year2 Modeling/Initial Data/Data 1965-2016.csv")
    start_date = "1/1/1993"
    series_names = c("Inflation", "Unemployment", "Treasury10Yr")
    first_observation = which(full_dataset$Date == start_date)
    series_full = as.list(full_dataset %>% slice((first_observation-1):nrow(full_dataset)) %>% select(series_names))
    true_vals = data.frame(matrix(unlist(series_full), nrow=length(series_full[[1]])))
    true_values = true_vals[-1,]
    forecast_sample_cumsum = apply(forecast_sample_store, c(2,3,4), cumsum)
    forecast_sample_real = forecast_sample_cumsum
    for (t in 1:timesteps){
      for (s in 1:num_series){
        forecast_sample_real[,s,,t] = forecast_sample_cumsum[,s,,t] + true_vals[t,s]
      }
    }
    cred_int_store = apply(forecast_sample_real, c(1, 2, 4), function(x) quantile(x, quantiles))
    quantile_names = c("lower", "cred_25", "prediction", "cred_75", "upper")
    dimnames(cred_int_store) = list("quantile" = quantile_names, "k" = 1:forecast_length, "series" = 1:num_series, "time" = 1:time)
    
  }else{
    max_time = min(time + k, nrow(environments[[1]]$X_full))
    series_tmp = lapply(environments, function(x) x$y_full[1:max_time] )
    true_values = data.frame(matrix(unlist(series_tmp), nrow=length(series_tmp[[1]])))
  }
  
  # Replace the median prediction with the mean! See if it makes a difference
  cred_int_store[3,,,] = mean_store
  
  plot_data = array2df(cred_int_store, levels = list(TRUE, NA, TRUE, TRUE)) %>%
    rename(value = cred_int_store) %>%
    spread(key = quantile, value = value) %>%
    mutate(time = as.numeric(time), series = as.numeric(series)) %>%
    mutate(time_of_pred = time) %>%
    mutate(time = time + k - 1) #Note that now time is the time being predicted!
  
  
  names(true_values) = 1:num_series
  true_values = true_values %>% mutate(time = row_number()) %>% gather("series", "true_value", 1:num_series) %>%
    mutate(series = as.numeric(series))
  
  plot_data = full_join(plot_data, true_values, by = c("time", "series"))
  
  plot_data = plot_data %>% filter(time > training_period)
}

get_kdens_data <- function(k_step_density_store, time, model, training_period){
  out = data_frame(time = 1:time, kdens = k_step_density_store, model=model)
  out = out[-(1:training_period),]
  out
}

get_kdens_split_data <- function(k_step_density_split_store, time, model, training_period){
  out = as.data.frame(k_step_density_split_store)
  names(out) = names
  out$time = 1:time
  out$model = model
  
  out = out %>% 
    gather("series", "lpfds", 1:num_series) %>%
    filter(time > training_period)
  
  out
}

get_kdens_comparison_data <- function(kdens_list){
  models = lapply(kdens_list, function(x) x$model)
  kdens = lapply(kdens_list, function(x) x$kdens)
  time = lapply(kdens_list, function(x) x$time)
  
  kdens_data = data_frame(time = unlist(time), model = unlist(models), k_step_path_density = unlist(kdens))
  kdens_data
}

get_kdens_comparison_split_data <- function(kdens_list, series_name){
  kdens_list = lapply(kdens_list, function(d) d %>%
                        filter(series == series_name) %>%
                        rename(kdens = lpfds))
  get_kdens_comparison_data(kdens_list)
}

## Plot the predictions
make_standard_plots <- function(plot_data, model_store, predictor_inclusion_probabilities, series_number, interventions, k_step_density_store, k, training_period, pdf_name){
  
  xnames = environments[[series_number]]$xnames
  plot_list = list()
  
  ylims = range(plot_data %>% filter(series == series_number) %>% pull(true_value))
  mean = mean(plot_data %>% filter(series == series_number) %>% pull(true_value))
  ylims[1] = ylims[1] - (mean - ylims[1])
  ylims[2] = ylims[2] - (mean - ylims[2])
  for(i in 1:k){
    p_preds = ggplot(plot_data %>% filter(series == series_number & k == i), aes(x = time, y = prediction))
    predictive_plot = p_preds + 
      geom_ribbon(aes(ymax = upper, ymin=lower), fill = "grey90") +
      geom_ribbon(aes(ymax = cred_75, ymin=cred_25), fill="grey70") +
      geom_line(aes(y = true_value, color = "Data", linetype="Data")) +
      geom_line(aes(color = "Predicted Value", linetype="PredictedValue")) +
      #  geom_vline(xintercept = input$time) +
      ggtitle(paste(i, "-step Prediction vs Data"))  +
      ylab("Prediction") +
      theme_bw() +
      geom_vline(xintercept = (interventions-1+i), color="black") +
      theme(plot.title = element_text(size=22), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.text = element_text(size=11),
            axis.text.x = element_text(angle=90)) + 
      scale_color_brewer(type="qual", palette=6) +
      scale_x_continuous(name = "Time", breaks=tticks, labels=tdates) +
      coord_cartesian(ylim = ylims)
    plot_list = c(plot_list, list(predictive_plot))
  }
  
  
  ## Plot the models
  tmp = sapply(model_store, function(model) environments[[series_number]]$Omega %in% model)
  plot_models = data.frame(tmp) %>% mutate("theta" = 1:nrow(tmp))
  names(plot_models)[1:(time+1)] = as.character(0:time)
  plot_models = gather(plot_models, "time", "in_model", 1:length(model_store)) %>%
    mutate(time = as.numeric(time)) %>%
    filter(time > training_period)
  names(plot_models)[1] = "theta"
  p = ggplot(plot_models, aes(x = as.numeric(time), y = theta))
  model_plot = p + geom_tile(aes(fill = in_model), color="black") + 
    #ggtitle("Z-Vector over Time") + theme(plot.title = element_text(size=22)) + 
    scale_x_continuous(name = "Time", breaks=tticks, labels=tdates) + 
    scale_y_continuous(name = "Covariates", breaks=1:length(xnames), labels = xnames) +
    # scale_fill_brewer(palette = "Accent")+
    # scale_fill_manual(values=c("white", 'black'))+
    theme(axis.text.x = element_text(size=11, angle=90), axis.text.y = element_text(size=8, hjust=0))
  #ggsave("Z-Vector over Time.png", device="png", width=10, height=5)
  
  #ggsave("predictions with confidence bands.png", device="png", width=10, height=5)
  ## Plot the predictor inclusion probabilities across the models found by SSS
  pip = tbl_df(data.frame(round(predictor_inclusion_probabilities, 3)))
  vars = ncol(pip)
  names(pip) = 1:vars
  pip = pip %>% mutate("time" = row_number())
  pip_plot = pip %>% gather("variable", "inclusion_probability", 1:(ncol(pip)-1))
  pip_plot$variable = as.numeric(pip_plot$variable)
  pip_plot = pip_plot %>% filter(time > training_period)
  p = ggplot(pip_plot, aes(x = time, y = variable))
  variable_inclusion_plot = p + geom_tile(aes(fill = inclusion_probability), color="black") +
    scale_x_continuous(name = "Time", breaks=tticks, labels=tdates) + 
    scale_y_continuous(name = "Covariates", breaks=1:length(xnames), labels = xnames) +
    scale_fill_continuous(limits = c(0, 1)) +
    theme(axis.text.x = element_text(size=11, angle=90), axis.text.y = element_text(size=8, hjust=0))
  
  ## Plot the 12-step ahead density
  # dat = data.frame(k_step_density_store, row.names=NULL, check.names=FALSE) %>% mutate_all(funs(cumsum)) %>% mutate(Time = row_number())
  # names(dat)[1:length(names)] = names
  # 
  # dat = dat %>% gather("Series", "Cumulative_k_Step_Log_Predictive_Density", 1:3) %>% filter(Series == names[series_number])
  # p = ggplot(dat, aes(x = Time, y = Cumulative_k_Step_Log_Predictive_Density))
  # log_pred_density_plot = p + geom_line(aes(color = Series)) +
  #   scale_x_continuous(name = "Time of Prediction", breaks=tticks, labels=tdates) +
  #   scale_y_continuous(name = paste("Cumulative ", k, "-step Log Predictive Density", sep="")) +
  #   theme(axis.text.x = element_text(size=11, angle=90)) +
  #   ggtitle(paste("Cumulative ", k, "-step Log Predictive Density", sep="")) +
  #   theme_bw() +
  #   geom_vline(xintercept = (interventions-1+i), color="black")
  
  # Marginal Percent Error
  dat = plot_data %>% filter(series == series_number & !is.na(k)) %>% group_by(k) %>%
    transmute(percerror = 100*abs(prediction - true_value)/true_value) %>%
    summarize(mean(percerror))
  names(dat) = c("k", "Perc_Error")
  p = ggplot(dat, aes(x = k, y=Perc_Error))
  percerr_plot = p +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_x_continuous(breaks=1:k) +
    ggtitle(paste("Average Marginal % Error", sep=""))
  
  ## Plot the k-step ahead rMSE and MAD (i.e. the L2 and L1 norms of the k-step ahead path forecast error over time)
  dat = plot_data %>% filter(series == series_number & !is.na(k)) %>% group_by(time_of_pred) %>%
    transmute(error = abs(prediction - true_value), sqerror = (prediction - true_value)^2) %>%
    summarize(rMSE = sqrt(mean(sqerror)), mean(error))
  names(dat) = c("Time_of_Prediction", "rMSE", "MAD")
  dat = dat %>% gather("Metric", "Value", 2:3)
  
  #p = ggplot(dat, aes(x = Time_of_Prediction, y=Value, color=Metric))
  p = ggplot(dat %>% filter(Metric == 'rMSE'), aes(x = Time_of_Prediction, y=Value, color=Metric)) # or just the MSE
  error_plot = p +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_x_continuous(name = "Time of Prediction", breaks=tticks, labels=tdates) +
    ggtitle(paste(k, "-Step Path Root Mean Squared Error", sep=""))
  
  # This plots the marginal rMSE at 1:k steps ahead
  dat = plot_data %>% filter(series == series_number & !is.na(k)) %>% group_by(k) %>%
    transmute(error = abs(prediction - true_value), sqerror = (prediction - true_value)^2) %>%
    summarize(sqrt(mean(sqerror)), mean(error))
  names(dat) = c("k", "rMSE", "MAD")
  p = ggplot(dat, aes(x = k, y=rMSE))
  rMSE_plot = p +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_x_continuous(breaks=1:k) +
    ggtitle(paste("Marginal", k, "-Step Root Mean Squared Error", sep=""))
  
  ## Plot the k-step ahead MAD
  # This plots the marginal MAD at 1:k steps ahead
  p = ggplot(dat, aes(x = k, y=MAD))
  MAD_plot = p +
    geom_point() +
    geom_line() +
    theme_bw() +
    scale_x_continuous(breaks=1:k) +
    ggtitle(paste("Marginal", k, "-Step Root Mean Absolute Deviation", sep=""))
  
  
  plot_list = c(plot_list, list(model_plot, variable_inclusion_plot, #log_pred_density_plot,
                                percerr_plot, error_plot, rMSE_plot, MAD_plot))
  pdf(file = pdf_name, width=10, height=5)
  invisible(lapply(plot_list, function(plot) print(plot)))
  dev.off()
  #ggsave(filename = pdf_name, device = "pdf", marrangeGrob(grobs = plot_list, nrow = 3, ncol = 1))
}

standard_theme <- function(plot, xlab = "Time", ylab = "Value", extra_theme = NULL, legend_title = NULL){
  plot + theme_classic() +
    scale_x_continuous(name = xlab, breaks=tticks[seq(1, 12, by=2)], labels=tdates[seq(1, 12, by=2)]) +
    #ggtitle("k-step Path Forecast Density") + 
    scale_y_continuous(name = ylab) +
    theme(axis.title = element_text(size=16, face="bold"), #14
          axis.text.x = element_text(angle = 0, hjust = 1, size=13), #11
          axis.text.y = element_text(angle = 0, size=13), #11
          legend.text = element_text(size=14), #12
          legend.position="none",
          plot.title = element_text(size=20),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text = element_text(angle = 0, size = 14, face="bold"), #14
          strip.text.y = element_text(angle = -90, size = 14, face="bold")) +
    scale_color_brewer(name = "Color", type="qual", palette=6)
}

raw_data <- function(series, series_names, time, k){
  plot_data = as.data.frame(series) %>%
    mutate(time = -(k - 1):time) %>%
    filter(time >= 1) %>% 
    gather("series", "value", 1:num_series) %>%
    mutate(series = factor(series, levels = series_names, labels = series_names))
  
  p = ggplot(plot_data, aes(x = time, y = value)) +
    geom_line() + 
    facet_wrap(~series, dir="v", scales = "free_y", strip.position = "left")
  
  standard_theme(p, ylab = NULL)
}

make_jpg <- function(p, title, width = 8, height = 4, res = 400){
  jpeg(paste(title, ".jpg", sep=""), width=width, height=height, units="in", res=res)
  print(p)
  dev.off()
  return(NULL)
}
