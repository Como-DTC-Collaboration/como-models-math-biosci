library(grid)
library(gridExtra)
library(ggpmisc)

# Define the log-likelihood function
RealData_LogLikelihoodFn <- function(parameters, model=SEIRD(),
                                     inc_numbers=0, death_numbers=0) {
  
  # Set up parameters
  transmission_para <- list(beta=parameters[1],
                            kappa=parameters[2],
                            gamma=parameters[3],
                            mu=parameters[4])
  
  transmission_parameters(model) <- transmission_para
  initial_conditions(model) <- init_cond
  
  # Simulate model
  times <- seq(0, length(inc_numbers), by = 1)
  out_df <- run(model, times)
  
  # Organise data and scale it with size of population
  data_wide <- spread(out_df$changes, compartment, value)
  data_wide$Incidence <- abs(data_wide$Incidence) * population_size
  data_wide$Deaths <- abs(data_wide$Deaths) * population_size
  
  # Computation of log-likelihood function
  logIncidence <- log(data_wide$Incidence[-1])
  incidence_likelihood <- sum(inc_numbers * logIncidence - data_wide$Incidence[-1])
  logDeaths <- log(data_wide$Deaths[-1])
  death_likelihood <- sum(death_numbers * logDeaths - data_wide$Deaths[-1])
  
  (incidence_likelihood + death_likelihood)
}


# Define function to plot data and optimised trajectory
data_optimisation_plot <- function(parameters_df, data, time_length_scale=1){
  
    # Extract optimised parameters from data frame
    optimised_para <- subset(parameters_df, select=c(beta_opt, kappa_opt,
                                                   gamma_opt, mu_opt))
    colnames(optimised_para) <- c("beta", "kappa", "gamma", "mu")
  
    # Simulate the model for all optimised parameters
    for (repeats in 1:nrow(parameters_df)){
        my_model <- SEIRD()
    
        transmission_parameters(my_model) <- optimised_para[repeats, 1:4]
        initial_conditions(my_model) <- init_cond
    
        times <- seq(0, length(data$DailyCases) * time_length_scale, by = 1)
        out_df <- run(my_model, times)
    
      cases <- out_df$changes
      cases$value <- cases$value * population_size
      cases$repeat_id <- rep(repeats, length(cases$time))
    
      if (repeats == 1){
          cases_repeats <- cases
      } else{
          cases_repeats <- rbind(cases_repeats, cases)
      }
    }

    # Visualise simulated trajectory
    cases_repeats <- spread(cases_repeats, compartment, value)
    ymax <- ceiling(max(data$DailyCases, cases_repeats$Incidence)) 
    
    plots <- vector('list', 2)
    for (plot_ind in 1:2){
      plots[[plot_ind]] <- ggplot()}
      
    for (i in 1:nrow(parameters_df)){
        data_set <- cases_repeats[cases_repeats$repeat_id == i,]
        
  	  plots[[1]] <- plots[[1]] +
          geom_line(data=data_set, aes(x = time, y = Incidence,
                                       color="Projected Incidence"))
  	      
      plots[[2]] <- plots[[2]] +
          geom_line(data=data_set, aes(x = time, y = Deaths,
                                       color="Projected Deaths"))
    }
    plots[[1]] <- plots[[1]] +
      geom_line(data=data, aes(x = time, y = DailyCases,
                               color="Actual Incidence"))
    plots[[2]] <- plots[[2]] +
      geom_line(data=data, aes(x = time, y = DailyDeaths,
                               color="Actual Deaths"))
    plots[[1]] <- plots[[1]] +
        scale_color_manual(name = "",
                           values = c("Projected Incidence" = "#FF6699",
                                      "Actual Incidence" = "darkred")) +
        ggtitle("Incidence") +
        labs(x = "time (days)", y = "counts") +
        ylim(0, ymax) +
        theme(text = element_text(size = 12), legend.position = "bottom",
              legend.text = element_text(size = 10),
              legend.direction = "vertical",
              legend.box.spacing = unit(0.01, "cm"))
    plots[[2]] <- plots[[2]] +
        scale_color_manual(name = "",
                           values = c("Projected Deaths" = "#66CCFF",
                                      "Actual Deaths" = "blue")) +
        ggtitle("Deaths") +
        labs(x = "time (days)", y = "counts") +
        ylim(0, ymax) +
        theme(text = element_text(size = 12), legend.position = "bottom",
              legend.text = element_text(size = 10),
              legend.direction = "vertical",
              legend.box.spacing = unit(0.01, "cm"))
    grid.arrange(grobs = plots, nrow = 1, ncol = 2)
}


# Define function to plot data and optimised trajectory
optimisation_plot <- function(parameters_df, data,
							  time_length_scale=1) {
  
  # Extract optimised parameters from data frame
  optimised_para <- subset(parameters_df, select=c(beta_opt, kappa_opt,
                                                   gamma_opt, mu_opt))
  colnames(optimised_para) <- c("beta", "kappa", "gamma", "mu")
  
  # Simulate the model for all optimised parameters
  for (repeats in 1:nrow(parameters_df)){
    my_model <- SEIRD()
    
    transmission_parameters(my_model) <- optimised_para[repeats, 1:4]
    initial_conditions(my_model) <- init_cond
    
    times <- seq(0, length(data$DailyCases) * time_length_scale, by = 1)
    out_df <- run(my_model, times)
    
    cases <- out_df$changes
    cases$value <- cases$value * population_size
    cases$repeat_id <- rep(repeats, length(cases$time))
    
    if (repeats == 1){
      cases_repeats <- cases
    } else{
      cases_repeats <- rbind(cases_repeats, cases)
    }
  }

  # Visualise simulated trajectory
  cases_repeats <- spread(cases_repeats, compartment, value)
  ymax <- ceiling(max(cases_repeats$Incidence)) 
  
  plots <- ggplot()
  for (i in 1:nrow(parameters_df)){
    data_set <- cases_repeats[cases_repeats$repeat_id == i,]
    
    plots <- plots +
      geom_line(data=data_set, aes(x = time, y = Incidence,
                                   color="Projected Incidence")) +
      geom_line(data=data_set, aes(x = time, y = Deaths,
                                   color="Projected Deaths"))
  }
  plots <- plots +
    scale_color_manual(name = "",
                       values = c("Projected Incidence" = "#FF6699",
                                  "Projected Deaths" = "#66CCFF")) +
    labs(x = "time (days)", y = "counts") +
    ylim(0, ymax) +
    theme(text = element_text(size = 12), legend.position = "bottom",
          legend.text = element_text(size = 10))
  print(plots)
}


# Define log-likelihood function to allow fixing of parameters
LogLikelihoodFn <- function(parameters, model=SEIRD(), inc_numbers=0, death_numbers=0,
                            profile_parameters=profile_parameters, fixed_parameter=FALSE,
                            fixed_parameter_value=0) {
  
  # Set up parameters for the SEIRD model
  optimising_para <- simulating_para
  for (i in 1:length(profile_parameters)){
    optimising_para[profile_parameters[i]] <- parameters[i]
  }
  
  # Fixing parameter of interest
  if (fixed_parameter != FALSE){
    optimising_para[names(optimising_para) == fixed_parameter] <-
      fixed_parameter_value}
  transmission_parameters(model) <- optimising_para[c('beta', 'kappa', 'gamma', 'mu')]
  initial_conditions(model) <- optimising_para[c('S0', 'E0', 'I0', 'R0')]
  
  # Simulate model 
  times <- seq(0, length(inc_numbers), by = 1)
  out_df <- run(model, times)
  
  # Organise and scale data
  data_wide <- spread(out_df$changes, compartment, value)
  data_wide$Incidence <- abs(data_wide$Incidence) * population_size
  data_wide$Deaths <- abs(data_wide$Deaths) * population_size
  
  # Calculate log-likelihood value
  logIncidence <- log(data_wide$Incidence[-1])
  incidence_likelihood <- sum(inc_numbers * logIncidence - data_wide$Incidence[-1])
  logDeaths <- log(data_wide$Deaths[-1])
  death_likelihood <- sum(death_numbers * logDeaths - data_wide$Deaths[-1])

  (incidence_likelihood + death_likelihood)
}


profilelikelihood_opt <- function(profile_parameters, range_transmission,
                                  likelihoodfn, synthetic_data, reltol=1e-6){
  # Set up structure
  previous_param <- rep(0, length(profile_parameters))
  profile_likelihood <- data.frame(parameter=character(), fixed_value=double(),
                                   likelihood_value=double())
  profile_likelihood$parameter <- as.character(profile_likelihood$parameter)
  for (fixed_param in profile_parameters){
    profile_likelihood[paste('optim_', fixed_param, sep = "")] <- double()
  }
  
  # Set up constraints
  constraint_ui <- rbind(diag(length(profile_parameters)))
  constraint_ci <- c(rep(0,length(profile_parameters)))
  
  # Run optimisation for all values in the range defined above
  set.seed(0)
  for (param_name in profile_parameters){
    err_value <- 0
    fixed_values <- range_transmission %>% filter(parameter == param_name)
    fixed_values <- fixed_values$fixed_value
    for (i in 1:length(fixed_values)){
      if (param_name == 'beta' & fixed_values[i] == simulating_para$beta){
        init_guess <- simulating_para[profile_parameters]
      } else if (param_name == 'kappa' & fixed_values[i] == simulating_para$kappa){
        init_guess <- simulating_para[profile_parameters]
      } else if (param_name == 'gamma' & fixed_values[i] == simulating_para$gamma){
        init_guess <- simulating_para[profile_parameters]
      } else if (param_name == 'mu' & fixed_values[i] == simulating_para$mu){
        init_guess <- simulating_para[profile_parameters]
      } else {
        init_guess <- previous_param}
      init_guess <- as.numeric(init_guess)
      result <- constrOptim(init_guess, likelihoodfn, 'NULL', constraint_ui,
                            constraint_ci, method = "Nelder-Mead",
                            control=list(fnscale=-1, reltol=reltol),
                            model=model,
                            inc_numbers = synthetic_data$IncNoise,
                            death_numbers = synthetic_data$DeathNoise,
							profile_parameters = profile_parameters,
                            fixed_parameter = param_name,
                            fixed_parameter_value = fixed_values[i])
      
      # Save the log-likelihood value and value of optimised parameters 
      err_value <- append(err_value, result$value)
      previous_param <- result$par
      profile_likelihood[
        nrow(profile_likelihood) + 1,] <- c(
          param_name, fixed_values[i],
          result$value, result$par)
    }
  }
  ind <- seq(2, 3 + length(profile_parameters), by = 1)
  profile_likelihood[ ,ind] <- apply(profile_likelihood[ ,ind], 2,
                                     function(x) as.numeric(x))
  profile_likelihood
}

# Create function to plot the profile likelihood
profilelikelihood_plot <- function(profile_likelihood, profile_parameters){

  # Get minimum and maximum of all log-likelihood values
  ymin <- floor(min(profile_likelihood$likelihood_value))
  ymax <- ceiling(max(profile_likelihood$likelihood_value))

  # Plot the log-likelihood values for each parameters of interest
  profile_plot <- vector('list', length(profile_parameters))
  for (i in 1:length(profile_parameters)){
    plot_data <- profile_likelihood %>% filter(parameter == profile_parameters[i])
  	profile_plot[[i]] <- ggplot(plot_data, aes(x = fixed_value,
                                               y = likelihood_value)) +
      geom_vline(xintercept = as.numeric(simulating_para[profile_parameters[i]])) +
      geom_point() + 
      labs(x = paste("fixed value"), y = "log-likelihood value") +
      ggtitle(profile_parameters[i]) +
	  scale_y_continuous(limits = c(ymin, ymax), labels = function(x) format(x, scientific = TRUE)) +
      theme(text = element_text(size = 12))
  }
  grid.arrange(grobs=profile_plot, nrow = ceiling(length(profile_parameters)/2), ncol = 2)
}

