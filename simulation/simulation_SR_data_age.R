# goal - simulate age structures Spawner-Recruit data
# similar to Peacock and Holt 2012, but for chum, and without straying

library(tidyverse)

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

simulated_chum_proportions <- generate_return_proportions(
  n_years = 20, 
  p_mean = chum_p_mean, 
  ages = chum_ages, 
  u = 0.9 # Using the variation value from Peacock & Holt 2012 
)

bh_data <- bh_function_w_age(mean_harvest = 0.3, 
                             sd_harvest = 0.2, 
                             K = 10000, 
                             alpha_mean = 1.5, 
                             sigma_mean = 1, 
                             ages = chum_ages, 
                             p_mean = chum_p_mean)

bh_data$generating_model <- "Beverton-Holt"

ric_data <- ric_function_w_age(mean_harvest = 0.3, 
                             sd_harvest = 0.2, 
                             K = 10000, 
                             alpha_mean = 5.5, 
                             sigma_mean = 1, 
                             ages = chum_ages, 
                             p_mean = chum_p_mean)

ric_data$generating_model <- "Ricker"

data <- rbind(ric_data,bh_data)

#plot as timeline

ggplot(data) + 
  geom_line(aes(x = year, y = S, color = generating_model, group = generating_model), 
            # color = "cadetblue"
  ) +
  labs(x = "Year",
       y = "Spawners") +
  theme_classic()


#plot data R_S vs S

ggplot(data) + 
  geom_point(aes(x = S, y = log(R/S), 
                 color = generating_model
                 ), 
             size = 2, alpha = 0.5) +
  facet_wrap(~ paste("alpha",alpha)+ paste("sigma",sigma)+ paste("Rk", Rk), scales = "free") + 
  # scale_color_gradient2(name = 'CPD std',
  #                       low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 0) +
  labs(#title = paste("alpha = ",mean(alpha), "sigma = ", mean(sigma), "forestry effect = ",mean(forestry_effect)),
    x = "Spawners (S)",
    y = "log(Recruits/Spawners) ") +
  theme_classic()


#Try generating data for pink salmon

pink_ages <- c(2)
pink_p_mean <- c(1)

ages <- pink_ages
p_mean <- pink_p_mean

simulated_pink_proportions <- generate_return_proportions(
  n_years = 20, 
  p_mean = pink_p_mean, 
  ages = pink_ages, 
  u = 0.9 # Using the variation value from Peacock & Holt 2012 
)

bh_data_pink <- bh_function_w_age(mean_harvest = 0.3, 
                             sd_harvest = 0.2, 
                             K = 10000, 
                             alpha_mean = 1.5, 
                             sigma_mean = 1, 
                             ages = pink_ages, 
                             p_mean = pink_p_mean)

bh_data_pink$generating_model <- "Beverton-Holt"

ric_data_pink <- ric_function_w_age(mean_harvest = 0.3, 
                               sd_harvest = 0.2, 
                               K = 10000, 
                               alpha_mean = 5.5, 
                               sigma_mean = 1, 
                               ages = pink_ages, 
                               p_mean = pink_p_mean)

ric_data_pink$generating_model <- "Ricker"

data <- rbind(ric_data_pink,bh_data_pink)



ggplot(data) + 
  geom_line(aes(x = year, y = S, color = generating_model, group = generating_model), 
            # color = "cadetblue"
  ) +
  labs(x = "Year",
       y = "Spawners") +
  theme_classic()

ggplot(data) + 
  geom_point(aes(x = S, y = log(R/S), 
                 color = generating_model
  ), 
  size = 2, alpha = 0.5) +
  facet_wrap(~ paste("alpha",alpha)+ paste("sigma",sigma)+ paste("Rk", Rk), scales = "free") + 
  # scale_color_gradient2(name = 'CPD std',
  #                       low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 0) +
  labs(#title = paste("alpha = ",mean(alpha), "sigma = ", mean(sigma), "forestry effect = ",mean(forestry_effect)),
    x = "Spawners (S)",
    y = "log(Recruits/Spawners) ") +
  theme_classic()

