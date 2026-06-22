# simulate bh and ricker data with age structure, then fir the two models to the data
#estimate the parameters, have K and Smsy in generated quantities

# first simulate multiple datasets with the same parameters and keep track of estimated medians
# then simulate multiple datasets with different parameters and keep track of estimated medians

library(here)
library(ggplot2)
suppressPackageStartupMessages(library(rstan))
rstan_options("auto_write" = TRUE)
library(PNWColors)
library(tidyverse)
library(gsl)

bh_function_w_age <- function(mean_harvest, sd_harvest, K_max, alpha_mean, sigma_mean, ages, p_mean, burn_in = 50, variation = 0.9){
  
  if (sd_harvest^2 >= mean_harvest * (1 - mean_harvest)) {
    stop("sd_harvest is too large for the given mean_harvest. Decrease sd_harvest.")
  }
  
  # number of years - random between 20 and 50
  years <- sample(26:50, 1)
  
  total_years <- years + burn_in #simulate for total an then discard data for burn_in years
  
  # harvest rate varying with year, beta distribution with mean = 0.3, sd = 0.2
  
  harvest_rate <- rbeta(total_years, shape1 = mean_harvest * ((mean_harvest * (1 - mean_harvest) / sd_harvest^2) - 1), 
                        shape2 = (1 - mean_harvest) * ((mean_harvest * (1 - mean_harvest) / sd_harvest^2) - 1))
  
  
  prop <- generate_return_proportions(
    n_years = total_years, 
    p_mean = p_mean, 
    ages = ages, 
    u = variation # Using the variation value for coho from Peacock & Holt 2012 
  )
  
  
  
  #Random K
  K = sample(seq(0.8 * K_max, K_max), 1)
  
  alpha_sample <- rnorm(100, alpha_mean, sd = 1)
  
  alpha <- sample(alpha_sample[alpha_sample > 0 & alpha_sample < 10], 1)
  
  sigma_sample <- rnorm(100, sigma_mean, 1)
  
  sigma <- sample(sigma_sample[sigma_sample > 0 & sigma_sample < 2], 1)
  
  colnames_R <-paste0("R_",ages)
  
  R_S = data.frame(matrix(nrow = total_years, ncol = length(ages)+3))
  
  colnames(R_S) <- c("S", colnames_R, "Run", "R")
  
  R_S$year <- 1:total_years
  
  R_S$alpha <- alpha
  
  R_S$sigma <- sigma
  
  R_S$K <- K
  
  R_S$harvest_rate <- harvest_rate
  
  max_age <- max(ages)
  
  for(t in 1:total_years){
    epsilon <- rnorm(1, mean = 0, sd = sigma)
    if(t<=max_age){
      #random initial spawner abundance near carrying capacity for the first max_age years
      R_S$S[t] <- K * runif(1, 0.8, 1.2)
      
    } else{
      
      total_run_for_year <- 0
      for(a in ages) {
        # Calculate returns for specific age 'a' from brood year 't-a'
        age_col <- paste0("Age_", a)
        return_col <- paste0("R_", a)
        
        # Recruits generated 'a' years ago * proportion returning at age 'a'
        age_returns <- R_S$R[t - a] * prop[t - a, age_col]
        
        R_S[t, return_col] <- age_returns
        total_run_for_year <- total_run_for_year + age_returns
      }
      
      
      R_S$Run[t] <- total_run_for_year
      
      R_S$S[t] <- R_S$Run[t]*(1-harvest_rate[t])
      
      
    }
    
    if(R_S$S[t] < 2){
      
      R_S$R[t] <- 0
      
    } else{
      
      Rk <- exp(alpha)*K/(exp(alpha) -1)
      
      R_S$R[t] <- R_S$S[t]*(exp(alpha)/(1 + exp(alpha)*R_S$S[t]/Rk))*exp(epsilon)
      
    }
    
    
    
    
  }
  #discard burn in years
  
  R_S_final <- R_S[(burn_in + 1):total_years,]
  
  R_S_final$year <- 1:years
  
  R_S_final$ln_RS <- log(R_S_final$R / R_S_final$S)
  
  R_S_final$Smsy <- R_S_final$K/(exp(R_S_final$alpha/2)+1)
  
  return(R_S_final)
  
}


ric_function_w_age <- function(mean_harvest, sd_harvest, K_max, alpha_mean, sigma_mean, ages, p_mean, burn_in = 50, variation = 0.9){
  
  if (sd_harvest^2 >= mean_harvest * (1 - mean_harvest)) {
    stop("sd_harvest is too large for the given mean_harvest. Decrease sd_harvest.")
  }
  # number of years - random between 20 and 50
  years <- sample(26:50, 1)
  
  total_years <- years + burn_in #simulate for total an then discard data for burn_in years
  
  # harvest rate varying with year, beta distribution with mean = 0.3, sd = 0.2
  
  harvest_rate <- rbeta(total_years, shape1 = mean_harvest * ((mean_harvest * (1 - mean_harvest) / sd_harvest^2) - 1), 
                        shape2 = (1 - mean_harvest) * ((mean_harvest * (1 - mean_harvest) / sd_harvest^2) - 1))
  
  
  prop <- generate_return_proportions(
    n_years = total_years, 
    p_mean = p_mean, 
    ages = ages, 
    u = variation # 0.9 Using the variation value from Peacock & Holt 2012 
  )
  
  
  
  #Random K
  K = sample(seq(0.8 * K_max, K_max), 1)
  
  alpha_sample <- rnorm(100, alpha_mean, sd = 1)
  
  alpha <- sample(alpha_sample[alpha_sample>0 & alpha_sample < 10], 1)
  
  sigma_sample <- rnorm(100, sigma_mean, 1)
  
  sigma <- sample(sigma_sample[sigma_sample>0 & sigma_sample < 2], 1)
  
  colnames_R <-paste0("R_",ages)
  
  R_S = data.frame(matrix(nrow = total_years, ncol = length(ages)+3))
  
  colnames(R_S) <- c("S", colnames_R, "Run", "R")
  
  R_S$year <- 1:total_years
  
  R_S$alpha <- alpha
  
  R_S$sigma <- sigma
  
  R_S$K <- K
  
  R_S$harvest_rate <- harvest_rate
  
  max_age <- max(ages)
  
  for(t in 1:total_years){
    
    epsilon <- rnorm(1, mean = 0, sd = sigma)
    
    if(t <= max_age){
      #random initial spawner abundance near carrying capacity for the first max_age years
      R_S$S[t] <- K * runif(1, 0.8, 1.2)
      
    } else{
      
      total_run_for_year <- 0
      for(a in ages) {
        # Calculate returns for specific age 'a' from brood year 't-a'
        age_col <- paste0("Age_", a)
        return_col <- paste0("R_", a)
        
        # Recruits generated 'a' years ago * proportion returning at age 'a'
        age_returns <- R_S$R[t - a] * prop[t - a, age_col]
        
        R_S[t, return_col] <- age_returns
        total_run_for_year <- total_run_for_year + age_returns
      }
      
      
      R_S$Run[t] <- total_run_for_year
      
      R_S$S[t] <- R_S$Run[t]*(1-harvest_rate[t])
      
      
    }
    
    if(R_S$S[t] < 2){
      
      R_S$R[t] <- 0
      
    } else{
      
      Smax <- K/alpha
      R_S$R[t]  <- R_S$S[t]*(exp(alpha - R_S$S[t]/Smax))*exp(epsilon)
      
    }
    
    
    
    
  }
  #discard burn in years
  
  R_S_final <- R_S[(burn_in + 1):total_years,]
  
  R_S_final$year <- 1:years
  
  R_S_final$ln_RS <- log(R_S_final$R / R_S_final$S)
  
  R_S_final$Smsy <- (1-lambert_W0(exp(1-R_S_final$alpha)))*R_S_final$K/R_S_final$alpha
  
  return(R_S_final)
  
}


generate_return_proportions <- function(n_years, p_mean, ages, u = 0.9) {
  
  # Ensure probabilities sum to 1
  if (abs(sum(p_mean) - 1) > 1e-6) {
    stop("The mean proportions (p_mean) must sum to 1.")
  }
  
  num_ages <- length(ages)
  
  # Initialize the output matrix
  p_gt_matrix <- matrix(0, nrow = n_years, ncol = num_ages)
  colnames(p_gt_matrix) <- paste0("Age_", ages)
  
  for (t in 1:n_years) {
    
    # Generate standard normal deviates for each age in year t 
    epsilon_gt <- rnorm(n = num_ages, mean = 0, sd = 1)
    
    # Calculate the inner term: log(p_mean) + scaling factor depends on species * error
    inner_term <- log(p_mean) + (u * epsilon_gt)
    
    # Calculate the mean of the inner term across all ages 
    # This corresponds to the summation divided by (G - a1 + 1) 
    mean_inner_term <- mean(inner_term)
    
    # Calculate the dummy variable x_{g,t} 
    x_gt <- inner_term - mean_inner_term
    
    # Calculate the final proportions p_{g,t} using the multivariate logistic equation 
    p_gt <- exp(x_gt) / sum(exp(x_gt))
    
    # Store in the matrix
    p_gt_matrix[t, ] <- p_gt
  }
  
  return(p_gt_matrix)
}

chum_ages <- c(3, 4, 5, 6)
chum_p_mean <- c(0.15, 0.65, 0.18, 0.02)

ages <- chum_ages
p_mean <- chum_p_mean

bh_data <- bh_function_w_age(mean_harvest = 0.3, 
                             sd_harvest = 0.2, 
                             K_max = 10000, 
                             alpha_mean = 1.5, 
                             sigma_mean = 1, 
                             ages = chum_ages, 
                             p_mean = chum_p_mean)

bh_data$generating_model <- "Beverton-Holt"

ric_data <- ric_function_w_age(mean_harvest = 0.3, 
                               sd_harvest = 0.2, 
                               K_max = 10000, 
                               alpha_mean = 5.5, 
                               sigma_mean = 1, 
                               ages = chum_ages, 
                               p_mean = chum_p_mean)

ric_data$generating_model <- "Ricker"

# data <- ric_data %>% 
#   filter(!is.nan(ln_RS), !is.infinite(ln_RS)) %>% 
#   rbind(bh_data %>% filter(!is.nan(ln_RS), !is.infinite(ln_RS)))
# 
data <- NULL



for(i in 1:nsims){
  
  set.seed(12345+i)
  bh_data <- bh_function_w_age(mean_harvest = 0.3, 
                               sd_harvest = 0.2, 
                               K_max = 10000, 
                               alpha_mean = 1.5, 
                               sigma_mean = 1, 
                               ages = chum_ages, 
                               p_mean = chum_p_mean)
  
  bh_data$generating_model <- "Beverton-Holt"
  bh_data$sim <- i
  
  ric_data <- ric_function_w_age(mean_harvest = 0.3, 
                                 sd_harvest = 0.2, 
                                 K_max = 10000, 
                                 alpha_mean = 1.5, 
                                 sigma_mean = 1, 
                                 ages = chum_ages, 
                                 p_mean = chum_p_mean)
  
  ric_data$generating_model <- "Ricker"
  ric_data$sim <- i
  
  data <- data %>% 
    rbind(ric_data %>% filter(!is.nan(ln_RS), !is.infinite(ln_RS))) %>% 
    rbind(bh_data %>% filter(!is.nan(ln_RS), !is.infinite(ln_RS)))
  
  
}


generating_model <- c("Beverton-Holt", "Ricker")
fitting_model <- c("Beverton-Holt", "Ricker")

sim_ric_model <- stan_model(file = here("simulation",
                                        "stan_models",
                                        "code",
                                        "ric_simple_model_for_simulated_data.stan"))

sim_bh_model <- stan_model(file = here("simulation",
                                       "stan_models",
                                       "code",
                                       "bh_simple_model_for_simulated_data.stan"))

model_results_w_spawner_combined_df <- data.frame(simulation = numeric(),
                                                  generating_model = character(),
                                                  parameter = character(),
                                                  true_value = numeric(),
                                                  fitting_model = character(),
                                                  estimate_median = numeric(),
                                                  estimate_lower = numeric(),
                                                  estimate_upper = numeric(),
                                                  Rhat = numeric(),
                                                  error = numeric())

# sample in parallel

options(mc.cores = parallel::detectCores())

nsims <- 500

for(i in 1:nsims){
  
  
  model <- sample(generating_model, 1)
  
  if(model == "Beverton-Holt"){
    
    data <- bh_function_w_age(mean_harvest = 0.3, 
                              sd_harvest = 0.2, 
                              K = 10000, 
                              alpha_mean = 1.5, 
                              sigma_mean = 1, 
                              ages = chum_ages, 
                              p_mean = chum_p_mean)
    
    data$generating_model <- "Beverton-Holt"
    
  } else{
    
    data <- ric_function_w_age(mean_harvest = 0.3, 
                              sd_harvest = 0.2, 
                              K = 10000, 
                              alpha_mean = 1.5, 
                              sigma_mean = 1, 
                              ages = chum_ages, 
                              p_mean = chum_p_mean)
    
    data$generating_model <- "Ricker"
    
  }
  
  true_values <- data %>% 
    group_by(sigma, alpha, K, Smsy) %>% 
    summarize(sigma = mean(sigma), 
              # forestry_effect = mean(forestry_effect), 
              alpha = mean(alpha), 
              # Smax = mean(Smax),
              # Rk = mean(Rk),
              K = mean(K),
              Smsy = mean(Smsy),
              generating_model = first(generating_model)) %>% 
    pivot_longer(cols = c(sigma, alpha, K, Smsy), names_to = "parameter", values_to = "true_value")
    
  
  data_list <- list(
    N = nrow(data),
    year = data$year,
    spawners = data$S,
    ln_RS = data$ln_RS,
    # forestry = data$forestry,
    Rk_mean = max(data$R),
    Rk_sigma = max(data$R)*2,
    Smax_mean = data$S[which.max(data$R)],
    Smax_sigma = data$S[which.max(data$R)]*2,
    prior_alpha = 5
  )
  for(fit_model in fitting_model){
    
    
    set.seed(i)
    
    
    
    if(fit_model == "Beverton-Holt"){
      
      model_sampling <- rstan::sampling(sim_bh_model,
                                                   data = data_list,
                                                   iter = 2000,
                                                   chains = 6,
                                                   warmup = 1000,
                                                   verbose = FALSE)
      
      
      
      
      
    } else if(fit_model == "Ricker"){
      
      model_sampling <- rstan::sampling(sim_ric_model,
                                                   data = data_list,
                                                   iter = 2000,
                                                   chains = 6,
                                                   warmup = 1000,
                                                   verbose = FALSE)
      
      
      
      
      
      
    }
    
    Rhat_values <- data.frame(Rhat = round(summary(model_sampling)$summary[,"Rhat"],3)) %>% 
      rownames_to_column("parameter")
    
    
    model_results <- tidybayes::spread_draws(model_sampling, alpha, sigma, K, Smsy) %>%
      mutate(fitting_model = fit_model, simulation = i) %>%
      select(fitting_model, alpha, sigma, K, Smsy, simulation) %>% 
      pivot_longer(cols = c(alpha, K, sigma, Smsy), names_to = "parameter", values_to = "value") %>%
      group_by(fitting_model, parameter, simulation) %>%
      summarise(
        estimate_median = round(median(value),2),
        estimate_lower = round(quantile(value, 0.025),2),
        estimate_upper = round(quantile(value, 0.975),2)
      ) %>%
      ungroup() %>% 
      
      left_join(true_values, by = "parameter") %>% 
      left_join(Rhat_values, by = "parameter") %>% 
      # mutate(data_model = "Ricker") %>% 
      mutate(error = 100*(estimate_median - true_value)/true_value)
    
    
    
    # print(true_values)
    
    
    
    
    
    model_results_w_spawner_combined_df <- model_results_w_spawner_combined_df %>%
      bind_rows(model_results)
  }
  
  
}



pal <- PNWColors::pnw_palette("Starfish", 5)


model_results_w_spawner_combined_df_new <- model_results_w_spawner_combined_df %>% 
  group_by(simulation, generating_model, fitting_model) %>%
  mutate(alpha = true_value[parameter == "alpha"],
         sigma = true_value[parameter == "sigma"],
         # b_for = true_value[parameter == "b_for"],
         Smsy = true_value[parameter == "Smsy"],
         K = true_value[parameter == "K"]
         ) %>%
  mutate(alpha_estimate = estimate_median[parameter == "alpha"],
         sigma_estimate = estimate_median[parameter == "sigma"],
         # b_for_estimate = estimate_median[parameter == "b_for"],
         Smsy_estimate = estimate_median[parameter == "Smsy"],
         K_estimate = estimate_median[parameter == "K"]
         ) %>%
  ungroup() 


write_csv(model_results_w_spawner_combined_df_new, here("simulation", 
                                                        "stan_models",
                                                        "output", 
                                                        "simulation_fitting_results_w_age.csv"))



ggplot(model_results_w_spawner_combined_df_new %>% 
         filter(parameter == "alpha"), aes(x = alpha, y = estimate_median)) +
  # geom_point(aes(color = alpha),alpha = 0.5, size = 2) +
  geom_pointrange(aes(ymin= estimate_lower, ymax = estimate_upper, color = Smsy), size = 0.5, alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ paste("Data model: ",generating_model) + paste("Fitting model: ",fitting_model)) +
  labs(x = "True alpha", y = "Estimated alpha") +
  # scale_color_gradient2(name = 'alpha',
  #                       low = pal[2], mid = 'gray', high = pal[4], midpoint = 5) +
  scale_color_gradientn(name = 'Estimated S_msy',
                        colors = pal)+
  theme_classic()

ggplot(model_results_w_spawner_combined_df_new %>% 
         filter(parameter == "K"), aes(x = K, y = estimate_median)) +
  geom_point(alpha = 0.5, size = 2) +
  # geom_pointrange(aes(ymin= estimate_lower, ymax = estimate_upper, color = Smsy), size = 0.5, alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ paste("Data model: ",generating_model) + paste("Fitting model: ",fitting_model)) +
  labs(x = "True K", y = "Estimated K") +
  scale_x_log10()+
  scale_y_log10()+
  # scale_color_gradient2(name = 'alpha',
  #                       low = pal[2], mid = 'gray', high = pal[4], midpoint = 5) +
  scale_color_gradientn(name = 'Estimated S_msy',
                        colors = pal)+
  theme_classic()

ggplot(model_results_w_spawner_combined_df_new %>% 
         filter(parameter == "Smsy"), aes(x = Smsy, y = estimate_median)) +
  # geom_point(aes(color = alpha),alpha = 0.5, size = 2) +
  geom_pointrange(aes(ymin= estimate_lower, ymax = estimate_upper, color = alpha), size = 0.5, alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ paste("Data model: ",generating_model) + paste("Fitting model: ",fitting_model)) +
  labs(x = "True Smsy", y = "Estimated Smsy") +
  # scale_color_gradient2(name = 'alpha',
  #                       low = pal[2], mid = 'gray', high = pal[4], midpoint = 5) +
  scale_color_gradientn(name = 'Estimated alpha',
                        colors = pal)+
  theme_classic()













