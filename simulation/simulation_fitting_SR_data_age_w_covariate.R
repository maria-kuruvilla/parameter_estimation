# Goal - simulate different data sets with the different parameter 1000 times with covariate effect



if(Sys.info()[7] == "mariakur") {
  print("Running on local machine")
  library(cmdstanr)
  set_cmdstan_path("C:/Users/mariakur/.cmdstan/cmdstan-2.35.0")
} else {
  print("Running on server")
  .libPaths(new = "/home/mkuruvil/R_Packages")
  library(cmdstanr)
  set_cmdstan_path("/home/mkuruvil/R_Packages/cmdstan-2.35.0")
}

library(here)
library(ggplot2)
# suppressPackageStartupMessages(library(rstan))
# rstan_options("auto_write" = TRUE)
library(PNWColors)
library(tidyverse)
library(gsl)

#read data

carnation_data <- read_csv(here("data",
                 "carnation_data.csv"))

sst_data <- read_csv(here("data",
                 "sst_ersst",
                 "ersst_spring.csv"))

forestry <- scale(c(rep(0,56), carnation_data$disturbedarea_prct_cs))

sst <- scale(sst_data$spring_ersst)


bh_function_w_age_w_covariates <- function(mean_harvest, sd_harvest, K, alpha, sigma, b_for, b_sst, ages, p_mean, years, burn_in = 50, variation = 0.9){
  
  if (sd_harvest^2 >= mean_harvest * (1 - mean_harvest)) {
    stop("sd_harvest is too large for the given mean_harvest. Decrease sd_harvest.")
  }
  
  
  
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
  
  
  
  
  
  colnames_R <-paste0("R_",ages)
  
  R_S = data.frame(matrix(nrow = total_years, ncol = length(ages)+5))
  
  colnames(R_S) <- c("S", colnames_R, "Run", "R", "forestry", "sst")
  
  R_S$year <- 1:total_years
  
  R_S$alpha <- alpha
  
  R_S$sigma <- sigma
  
  R_S$K <- K
  
  R_S$harvest_rate <- harvest_rate
  
  max_age <- max(ages)
  
  for(t in 1:total_years){
    epsilon <- rnorm(1, mean = 0, sd = sigma)
    R_S$forestry[t] <- forestry[t]
    R_S$sst[t] <- sst[t]
    
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
      
      R_S$R[t] <- R_S$S[t]*(exp(alpha)/(1 + exp(alpha)*R_S$S[t]/Rk))*exp(epsilon)*exp(b_for*forestry[t])*exp(b_sst*sst[t])
      
    }
    
    
    
    
  }
  #discard burn in years
  
  R_S_final <- R_S[(burn_in + 1):total_years,]
  
  R_S_final$year <- 1:years
  
  R_S_final$ln_RS <- log(R_S_final$R / R_S_final$S)
  
  R_S_final$Smsy <- R_S_final$K/(exp(R_S_final$alpha/2)+1)
  
  return(R_S_final)
  
}


ric_function_w_age_w_covariates <- function(mean_harvest, sd_harvest, K, alpha, sigma, b_for, b_sst, ages, p_mean, years, burn_in = 50, variation = 0.9){
  
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
  
  
  

  
  colnames_R <-paste0("R_",ages)
  
  R_S = data.frame(matrix(nrow = total_years, ncol = length(ages)+5))
  
  colnames(R_S) <- c("S", colnames_R, "Run", "R", "forestry", "sst")
  
  R_S$year <- 1:total_years
  
  R_S$alpha <- alpha
  
  R_S$sigma <- sigma
  
  R_S$K <- K
  
  R_S$harvest_rate <- harvest_rate
  
  max_age <- max(ages)
  
  for(t in 1:total_years){
    
    epsilon <- rnorm(1, mean = 0, sd = sigma)
    
    R_S$forestry[t] <- forestry[t]
    
    R_S$sst[t] <- sst[t]
    
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
      R_S$R[t]  <- R_S$S[t]*(exp(alpha - R_S$S[t]/Smax))*exp(epsilon)*exp(b_for*forestry[t])*exp(b_sst*sst[t])
      
    }
    
    
    
    
  }
  #discard burn in years
  
  R_S_final <- R_S[(burn_in + 1):total_years,]
  
  R_S_final$year <- 1:years
  
  R_S_final$ln_RS <- log(R_S_final$R / R_S_final$S)
  
  R_S_final$Smsy <- (1-lambert_W0(exp(1-R_S_final$alpha)))*R_S_final$K/R_S_final$alpha
  
  return(R_S_final)
  
}


bh_function_w_age_w_covariate <- function(mean_harvest, sd_harvest, K, alpha, sigma, b_cov, ages, p_mean, years, burn_in = 50, variation = 0.9){
  
  if (sd_harvest^2 >= mean_harvest * (1 - mean_harvest)) {
    stop("sd_harvest is too large for the given mean_harvest. Decrease sd_harvest.")
  }
  
  
  
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
  
  
  
  
  
  colnames_R <-paste0("R_",ages)
  
  R_S = data.frame(matrix(nrow = total_years, ncol = length(ages)+4))
  
  colnames(R_S) <- c("S", colnames_R, "Run", "R", "covariate")
  
  R_S$year <- 1:total_years
  
  R_S$alpha <- alpha
  
  R_S$sigma <- sigma
  
  R_S$K <- K
  
  R_S$harvest_rate <- harvest_rate
  
  max_age <- max(ages)
  
  for(t in 1:total_years){
    epsilon <- rnorm(1, mean = 0, sd = sigma)
    
    R_S$covariate[t] <- covariate[t]
    
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
      
      R_S$R[t] <- R_S$S[t]*(exp(alpha)/(1 + exp(alpha)*R_S$S[t]/Rk))*exp(epsilon)*exp(b_cov*covariate[t])
      
    }
    
    
    
    
  }
  #discard burn in years
  
  R_S_final <- R_S[(burn_in + 1):total_years,]
  
  R_S_final$year <- 1:years
  
  R_S_final$ln_RS <- log(R_S_final$R / R_S_final$S)
  
  R_S_final$Smsy <- R_S_final$K/(exp(R_S_final$alpha/2)+1)
  
  R_S_final$b_cov <- b_cov
  
  
  return(R_S_final)
  
}


ric_function_w_age_w_covariate <- function(mean_harvest, sd_harvest, K, alpha, sigma, b_cov, ages, p_mean, years, burn_in = 50, variation = 0.9){
  
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
  
  
  
  
  
  colnames_R <-paste0("R_",ages)
  
  R_S = data.frame(matrix(nrow = total_years, ncol = length(ages)+4))
  
  colnames(R_S) <- c("S", colnames_R, "Run", "R", "covariate")
  
  R_S$year <- 1:total_years
  
  R_S$alpha <- alpha
  
  R_S$sigma <- sigma
  
  R_S$K <- K
  
  R_S$harvest_rate <- harvest_rate
  
  max_age <- max(ages)
  
  for(t in 1:total_years){
    
    epsilon <- rnorm(1, mean = 0, sd = sigma)
    
    R_S$covariate[t] <- covariate[t]
    
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
      R_S$R[t]  <- R_S$S[t]*(exp(alpha - R_S$S[t]/Smax))*exp(epsilon)*exp(b_cov*covariate[t])
      
    }
    
    
    
    
  }
  #discard burn in years
  
  R_S_final <- R_S[(burn_in + 1):total_years,]
  
  R_S_final$year <- 1:years
  
  R_S_final$ln_RS <- log(R_S_final$R / R_S_final$S)
  
  R_S_final$Smsy <- (1-lambert_W0(exp(1-R_S_final$alpha)))*R_S_final$K/R_S_final$alpha
  
  R_S_final$b_cov <- b_cov
  
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


generating_model <- c("Beverton-Holt", "Ricker")
fitting_model <- c("Beverton-Holt", "Ricker")

sim_ric_model <- cmdstanr::cmdstan_model(file.path(here("simulation",
                                        "stan_models",
                                        "code",
                                        "ric_simple_model_for_simulated_data.stan")))

sim_bh_model <- cmdstanr::cmdstan_model(file.path(here("simulation",
                                       "stan_models",
                                       "code",
                                       "bh_simple_model_for_simulated_data.stan")))

sim_ric_w_covariate_model <- cmdstanr::cmdstan_model(file.path(here("simulation",
                                        "stan_models",
                                        "code",
                                        "ric_w_covariate_simple_model_for_simulated_data.stan")))

sim_bh_w_covariate_model <- cmdstanr::cmdstan_model(file.path(here("simulation",
                                       "stan_models",
                                       "code",
                                       "bh_w_covariate_simple_model_for_simulated_data.stan")))


sim_ric_w_covariates_model <- cmdstanr::cmdstan_model(file.path(here("simulation",
                                                    "stan_models",
                                                    "code",
                                                    "ric_w_covariates_simple_model_for_simulated_data.stan")))

sim_bh_w_covariates_model <- cmdstanr::cmdstan_model(file.path(here("simulation",
                                                   "stan_models",
                                                   "code",
                                                   "bh_w_covariates_simple_model_for_simulated_data.stan")))

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

nsims <- 100

#first for one covariate

covariate <- forestry

for(i in 1:nsims){
  
  set.seed(12345+i)
  
  # number of years - random between 20 and 50
  years <- sample(26:45, 1)
  
  alpha_mean = 1.5
  
  sigma_mean = 1
  
  K_max = 10000
  
  #Random K
  K = sample(seq(0.8 * K_max, K_max), 1)
  
  alpha_sample <- rnorm(100, alpha_mean, sd = 1)
  
  alpha <- sample(alpha_sample[alpha_sample > 0 & alpha_sample < 10], 1)
  
  sigma_sample <- rnorm(100, sigma_mean, 1)
  
  sigma <- sample(sigma_sample[sigma_sample > 0 & sigma_sample < 2], 1)
  
  b_sst <- rnorm(1, 0, 1)
  
  b_for <- rnorm(1, 0, 1)
  
  model <- sample(generating_model, 1)
  
  if(model == "Beverton-Holt"){
    
    data <- bh_function_w_age_w_covariate(mean_harvest = 0.3, 
                              sd_harvest = 0.2, 
                              K = K, 
                              alpha = alpha, 
                              sigma = sigma, 
                              # b_for = b_for,
                              b_cov = b_sst,
                              ages = chum_ages, 
                              p_mean = chum_p_mean,
                              years = years)
    
    data$generating_model <- "Beverton-Holt"
    
  } else{
    
    data <- ric_function_w_age_w_covariate(mean_harvest = 0.3, 
                              sd_harvest = 0.2, 
                              K = K, 
                              alpha = alpha, 
                              sigma = sigma, 
                              # b_for = b_for,
                              b_cov = b_sst,
                              ages = chum_ages, 
                              p_mean = chum_p_mean,
                              years = years)
    
    data$generating_model <- "Ricker"
    
  }
  
  data <- data %>% 
        filter(!is.nan(ln_RS), !is.infinite(ln_RS))
  
  #if data has <2 rows, then go to next simulation
  if(nrow(data) < 2){
    next
  }
  
  
  true_values <- data %>% 
    # group_by(sigma, alpha, K, Smsy) %>% 
    summarize(sigma = mean(sigma), 
              # forestry_effect = mean(forestry_effect), 
              alpha = mean(alpha), 
              # Smax = mean(Smax),
              # Rk = mean(Rk),
              K = mean(K),
              Smsy = mean(Smsy),
              # b_for = mean(b_for),
              b_cov = mean(b_cov),
              min_S = min(S),
              generating_model = first(generating_model)) %>% 
    pivot_longer(cols = c(sigma, alpha, K, Smsy, b_cov), names_to = "parameter", values_to = "true_value")
    
  
  data_list <- list(
    N = nrow(data),
    year = data$year,
    spawners = data$S,
    ln_RS = data$ln_RS,
    # forestry = data$forestry,
    covariate = data$covariate,
    Rk_mean = max(data$R),
    Rk_sigma = max(data$R)*2,
    Smax_mean = data$S[which.max(data$R)],
    Smax_sigma = data$S[which.max(data$R)]*2,
    prior_alpha = 5
  )
  for(fit_model in fitting_model){
    
    
    set.seed(12345+i)
    
    
    
    if(fit_model == "Beverton-Holt"){
      
      model_sampling <- sim_bh_w_covariate_model$sample(data = data_list,
                                            iter_sampling  = 2000,
                                            chains = 6,
                                            iter_warmup = 1000)
      
      
      
      
      
    } else if(fit_model == "Ricker"){
      
      model_sampling <- sim_ric_w_covariate_model$sample(data = data_list,
                                             iter_sampling  = 2000,
                                             chains = 6,
                                             iter_warmup = 1000)
      
      
      
      
      
      
    }
    
    Rhat_values <- data.frame(Rhat = round(model_sampling$summary()$rhat,3)) %>% 
      mutate(parameter = model_sampling$summary()$variable)
    
    model_results <- data.frame(model_sampling$draws(variables=c("alpha", "sigma", "K", "Smsy", "b_cov"),format='draws_matrix')) %>%
      mutate(fitting_model = fit_model, simulation = i) %>%
      select(fitting_model, alpha, sigma, K, Smsy, b_cov, simulation) %>% 
      pivot_longer(cols = c(alpha, K, sigma, Smsy, b_cov), names_to = "parameter", values_to = "value") %>%
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



model_results_w_spawner_combined_df_new <- model_results_w_spawner_combined_df %>% 
  group_by(simulation, generating_model, fitting_model) %>%
  mutate(alpha = true_value[parameter == "alpha"],
         sigma = true_value[parameter == "sigma"],
         b_cov = true_value[parameter == "b_cov"],
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
                                                        paste0("simulation_fitting_results_w_covariate_forestry_",nsims,".csv")))



if(Sys.info()[1] == "Windows") {
  
  
  pal <- PNWColors::pnw_palette("Starfish", 5)
  
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
    coord_trans(y="log10") + 
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_y_continuous(labels = scales::comma, breaks = c(1000, 10000, 100000)) +
    facet_wrap(~ paste("Data model: ",generating_model) + paste("Fitting model: ",fitting_model), scale = "free_y") +
    labs(x = "True Smsy", y = "Estimated Smsy") +
    scale_color_gradientn(name = 'Estimated alpha',
                          colors = pal)+
    theme_classic()
  
  ggplot(model_results_w_spawner_combined_df_new %>% 
           filter(parameter == "b_cov"), aes(x = b_cov, y = estimate_median)) +
    # geom_point(aes(color = alpha),alpha = 0.5, size = 2) +
    geom_pointrange(aes(ymin= estimate_lower, ymax = estimate_upper, color = alpha), size = 0.5, alpha = 0.5)+
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_y_continuous(labels = scales::comma, breaks = c(1000, 10000, 100000)) +
    facet_wrap(~ paste("Data model: ",generating_model) + paste("Fitting model: ",fitting_model), scale = "free_y") +
    labs(x = "True Smsy", y = "Estimated Smsy") +
    scale_color_gradientn(name = 'Estimated alpha',
                          colors = pal)+
    theme_classic()
  
  
} else{
  print("running on server")
  
  
  
}





