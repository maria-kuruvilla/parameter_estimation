library(here);library(dplyr); library(stringr)
library(ggplot2)
library(tidyverse)
library(bayesplot)
library(patchwork)
library(hues)
library(GGally)
library(latex2exp)
library(ggforce)
library(PNWColors)


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
  
  R_S_final$ln_RS <- log(R_S_final$R/R_S_final$S)
  
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
  
  R_S_final$ln_RS <- log(R_S_final$R/R_S_final$S)
  
  R_S_final$year <- 1:years
  
  return(R_S_final)
  
}


set.seed(123)

K_max <- 1000

alpha_mean <- 1.2

sigma_mean <- 1

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


K = sample(seq(0.8 * K_max, K_max), 1)

alpha_sample <- rnorm(100, alpha_mean, sd = 1)

alpha <- sample(alpha_sample[alpha_sample > 0 & alpha_sample < 10], 1)

sigma_sample <- rnorm(100, sigma_mean, 1)

sigma <- sample(sigma_sample[sigma_sample > 0 & sigma_sample < 2], 1)

chum_ages <- c(3, 4, 5, 6)
chum_p_mean <- c(0.15, 0.65, 0.18, 0.02)

ages <- chum_ages
p_mean <- chum_p_mean

ricker_data <- ric_function_w_age(mean_harvest = 0.3, 
                              sd_harvest = 0.2, 
                              K = K, 
                              alpha = alpha, 
                              sigma = sigma, 
                              ages = chum_ages, 
                              p_mean = chum_p_mean)

bh_data <- bh_function_w_age(mean_harvest = 0.3, 
                            sd_harvest = 0.2, 
                            K = K, 
                            alpha = alpha, 
                            sigma = sigma, 
                            ages = chum_ages, 
                            p_mean = chum_p_mean)



bh_data$model <- "BH"

bh_data <- bh_data %>% 
  filter(!is.nan(ln_RS), !is.infinite(ln_RS))


ricker_data$model <- "Ricker"

ricker_data <- ricker_data %>% 
  filter(!is.nan(ln_RS), !is.infinite(ln_RS))

my_data <- rbind(bh_data, ricker_data)

ggplot(my_data) + 
  geom_point(aes(x = S, y = log(R/S), color = model), size = 2, alpha = 0.5) +
  # facet_wrap(~ paste("alpha",alpha)+ paste("sigma",sigma)+ paste("forestry effect",forestry_effect), scales = "free") + 
  # scale_color_gradient2(name = 'CPD std',
  #                       low = '#35978f', mid = 'gray', high = '#bf812d', midpoint = 0) +
  scale_color_manual(values = c("BH" = "#7F719F", "Ricker" = "#7CAE7A"))+
  labs(title = paste("Simulated Data with alpha =", round(alpha,2), "and sigma =", round(sigma,2)),
       x = "Spawners (S)",
       y = "log(Recruits/Spawners) ") +
  theme_classic()

# ggsave(filename = here("forestry_simulation", "output", "simulated_bh_ricker_data.png"), width = 7, height = 5)


# ===================================================================
# 1. RICKER OBJECTIVE FUNCTION (FULLY PROFILEABLE)
# ===================================================================
nll_ricker <- function(pars, data, alpha_fixed = NULL, K_fixed = NULL, sigma_fixed = NULL) {
  idx <- 1 # Index tracker for the 'pars' vector
  
  # 1. Extract or Fix Alpha
  if (!is.null(alpha_fixed)) {
    alpha <- alpha_fixed
  } else {
    alpha <- pars[idx]
    idx <- idx + 1
  }
  
  # 2. Extract or Fix K
  if (!is.null(K_fixed)) {
    K <- K_fixed
  } else {
    K <- pars[idx]
    idx <- idx + 1
  }
  
  # 3. Extract or Fix Sigma
  if (!is.null(sigma_fixed)) {
    sigma <- sigma_fixed
  } else {
    sigma <- pars[idx]
  }
  
  # Safety penalty for impossible values
  if (alpha <= 0 || K <= 0 || sigma <= 0) {
    return(1e9) 
  }
  
  # Calculate Negative Log-Likelihood
  Smax <- K / alpha
  predicted_ln_RS <- alpha - (data$S / Smax)
  nll <- -sum(dnorm(data$ln_RS, mean = predicted_ln_RS, sd = sigma, log = TRUE))
  
  return(nll)
}

# ===================================================================
# 2. BEVERTON-HOLT OBJECTIVE FUNCTION (FULLY PROFILEABLE)
# ===================================================================
nll_bh <- function(pars, data, alpha_fixed = NULL, K_fixed = NULL, sigma_fixed = NULL) {
  idx <- 1 
  
  if (!is.null(alpha_fixed)) {
    alpha <- alpha_fixed
  } else {
    alpha <- pars[idx]
    idx <- idx + 1
  }
  
  if (!is.null(K_fixed)) {
    K <- K_fixed
  } else {
    K <- pars[idx]
    idx <- idx + 1
  }
  
  if (!is.null(sigma_fixed)) {
    sigma <- sigma_fixed
  } else {
    sigma <- pars[idx]
  }
  
  if (alpha <= 0 || K <= 0 || sigma <= 0) {
    return(1e9) 
  }
  
  # Calculate Negative Log-Likelihood
  Rk <- (exp(alpha) * K) / (exp(alpha) - 1)
  predicted_ln_RS <- alpha - log(1 + (exp(alpha) * data$S) / Rk)
  nll <- -sum(dnorm(data$ln_RS, mean = predicted_ln_RS, sd = sigma, log = TRUE))
  
  return(nll)
}

# ===================================================================
# 3. PROFILE LIKELIHOOD
# ===================================================================


# profile likelihood for Carrying Capacity (K)

# 1. Define the range of K values you want to test
K_sequence <- seq(100, 2000, by = 100)


optimized_ricker <- optim(
  par = c(1.5, 1000, 0.5),
  fn = nll_ricker, 
  data = my_data, 
)

optimized_bh <- optim(
  par = c(1.5, 1000, 0.5),
  fn = nll_bh, 
  data = my_data, 
)

best_nll_ricker <- optimized_ricker$value

best_nll_bh <- optimized_bh$value

profile_df <- data.frame(K = K_sequence, nll_ricker = NA_real_, nll_bh = NA_real_)

# 3. Loop through the sequence
for(i in 1:length(K_sequence)) {
  
  # Because K is fixed, optim() only needs starting guesses for alpha and sigma
  initial_guesses <- c(1.5, 0.5) # c(alpha, sigma)
  
  # Run the optimizer
  fit_ricker <- optim(
    par = initial_guesses, 
    fn = nll_ricker, 
    data = my_data, 
    K_fixed = K_sequence[i] # <-- Lock K for this specific loop iteration!
  )
  
  fit_bh <- optim(
    par = initial_guesses, 
    fn = nll_bh, 
    data = my_data, 
    K_fixed = K_sequence[i] # <-- Lock K for this specific loop iteration!
  )
  
  # Store the best (lowest) NLL achieved for this specific value of K
  profile_df$nll_ricker[profile_df$K == K_sequence[i]] <- fit_ricker$value
  
  profile_df$nll_bh[profile_df$K == K_sequence[i]] <- fit_bh$value
}

# 4. Plot the profile
plot(K_sequence, profile_nll_scores, type = "l", lwd = 2,
     xlab = "Carrying Capacity (K)", 
     ylab = "Negative Log-Likelihood",
     main = "Profile Likelihood for K")

ggplot() +
  geom_line(data = profile_df, aes(x = K, y = nll_ricker - best_nll_ricker, color = "Ricker"), size = 1, alpha = 0.5) +
  geom_line(data = profile_df, aes(x = K, y = nll_bh - best_nll_bh, color = "BH"), size = 1, alpha = 0.5) +
  #add horizontal line at confidence interval
  geom_hline(yintercept = qchisq(0.95, df = 1) / 2, linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("BH" = "#7F719F", "Ricker" = "#7CAE7A"))+
  geom_vline(xintercept = K, linetype = "dashed", color = "cadetblue") +
  labs(title = "Profile Likelihood for K in Ricker model",
       x = "K",
       y = "Negative Log Likelihood") +
  theme_classic()

# profile likelihood for alpha

alpha_seq <- seq(0.1, 5, by = 0.1)


optimized_ricker <- optim(
  par = c(1.5, 1000, 0.5),
  fn = nll_ricker, 
  data = my_data, 
)

optimized_bh <- optim(
  par = c(1.5, 1000, 0.5),
  fn = nll_bh, 
  data = my_data, 
)

best_nll_ricker <- optimized_ricker$value

best_nll_bh <- optimized_bh$value

profile_df_alpha <- data.frame(alpha = alpha_seq, nll_ricker = NA_real_, nll_bh = NA_real_)


for(i in 1:length(alpha_seq)) {
  
  initial_guesses <- c(1000, 0.5) # c(K, sigma)
  
  fit_ricker <- optim(
    par = initial_guesses, 
    fn = nll_ricker, 
    data = my_data, 
    alpha_fixed = alpha_seq[i] # <-- Lock alpha for this specific loop iteration!
  )
  
  fit_bh <- optim(
    par = initial_guesses, 
    fn = nll_bh, 
    data = my_data, 
    alpha_fixed = alpha_seq[i] # <-- Lock alpha for this specific loop iteration!
  )
  
  # Store the best (lowest) NLL achieved for this specific value of K
  profile_df_alpha$nll_ricker[profile_df_alpha$alpha == alpha_seq[i]] <- fit_ricker$value
  
  profile_df_alpha$nll_bh[profile_df_alpha$alpha == alpha_seq[i]] <- fit_bh$value
}


ggplot() +
  geom_line(data = profile_df_alpha, aes(x = alpha, y = nll_ricker - min(nll_ricker), color = "Ricker"), size = 1, alpha = 0.5) +
  geom_line(data = profile_df_alpha, aes(x = alpha, y = nll_bh - min(nll_bh), color = "BH"), size = 1, alpha = 0.5) +
  #add horizontal line at confidence interval
  geom_hline(yintercept = qchisq(0.95, df = 1) / 2, linetype = "dashed", color = "gray", size = 0.75) +
  scale_color_manual(values = c("BH" = "#7F719F", "Ricker" = "#7CAE7A"), name = "Model")+
  geom_vline(xintercept = alpha, linetype = "dashed", color = "cadetblue", size = 0.75) +
  labs(title = "Profile likelihood for alpha",
       x = "Alpha",
       y = "Negative Log Likelihood") +
  annotate("text", x = alpha-0.2, y = 50, label = "True alpha", color = "cadetblue", angle = 90) +
  annotate("text", x = 4, y = qchisq(0.95, df = 1)/2 + 2 , label = "95% CI threshold", 
           color = "gray", angle = 0) +
  theme_classic()

#save

ggsave(filename = here("profile_likelihoods", "output", "profile_likelihood_alpha.png"), width = 6, height = 4)

#save profile df alpha

write.csv(profile_df_alpha, file = here("profile_likelihoods", "output", "profile_likelihood_alpha.csv"), row.names = FALSE)


# profile likelihood for alpha when K is also fixed

alpha_seq <- seq(0.1, 5, by = 0.1)

K_fixed <- 1000

profile_nll_scores_alpha <- numeric(length(alpha_seq))


for(i in 1:length(alpha_seq)) {
  
  initial_guesses <- c(0.5) # Only sigma is free now
  
  fit <- optim(
    par = initial_guesses, 
    fn = nll_ricker, 
    data = my_data, 
    alpha_fixed = alpha_seq[i], 
    K_fixed = K_fixed
  )
  
  profile_nll_scores_alpha[i] <- fit$value
}









