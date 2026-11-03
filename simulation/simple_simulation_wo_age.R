# Goal - simulate ricker and bh dynamics without harvest rate, age structure
# allow spawners to go below 2 

# write it in cmdstanr

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
library(tidyverse)
library(gsl)

bh_function <- function(alpha, n0, sigma, K, years){
  R_S = data.frame(S = rep(NA, years), R = rep(NA, years), year = 1:years)
  
  R_S$n0 <- n0
  
  R_S$alpha <- alpha
  
  R_S$sigma <- sigma
  
  R_S$K <- K
  
  # Beverton-Holt function
  for(t in 1:years){
    epsilon <- rnorm(1, mean = 0, sd = sigma)
    if(t==1){
      R_S$S[1] <- n0
      
    } else{
      R_S$S[t] <- R_S$R[t-1]
      
    }
    
    Rk <- K/(exp(alpha) - 1)
    R_S$R[t] <- R_S$S[t]*(exp(alpha)/(1 + exp(alpha)*R_S$S[t]/Rk))*exp(epsilon)
    
  }
  R_S$ln_RS <- log(R_S$R/R_S$S)
  R_S$Smsy <- R_S$K/(exp(R_S$alpha/2)+1)
  R_S$msy <- R_S$K*(exp(R_S$alpha/2) - 1)/(exp(R_S$alpha/2) + 1)
  R_S$Umsy <- 1 - exp(-alpha/2)
  
  return(R_S)
}



ric_function <- function(alpha, n0, sigma, K, years){
  R_S = data.frame(S = rep(NA, years), R = rep(NA, years), year = 1:years)
  
  R_S$n0 <- n0
  
  R_S$alpha <- alpha
  
  R_S$sigma <- sigma
  
  R_S$K <- K
  
  for(t in 1:years){
    epsilon <- rnorm(1, mean = 0, sd = sigma)
    if(t==1){
      R_S$S[1] <- n0
      
    } else{
      R_S$S[t] <- R_S$R[t-1]
      
    }
    
    Smax <- K/alpha
    
    R_S$R[t] <- R_S$S[t]*(exp(alpha - R_S$S[t]/Smax))*exp(epsilon)
    
  }
  R_S$ln_RS <- log(R_S$R/R_S$S)
  R_S$Smsy <- (1-lambert_W0(exp(1-R_S$alpha)))*R_S$K/R_S$alpha
  R_S$msy <- (R_S$Smsy* exp(R_S$alpha - R_S$Smsy/Smax)) - R_S$Smsy
  R_S$Umsy <- R_S$msy/(R_S$Smsy + R_S$msy)
  
  return(R_S)
}

alpha_mean <- 5.4
n0 <- 100
sigma_mean <- 1
K <- 10000

alpha_sample <- rnorm(100, alpha_mean, sd = 1)

alpha <- sample(alpha_sample[alpha_sample>0 & alpha_sample < 10], 1)

sigma_sample <- rnorm(100, sigma_mean, 1)

sigma <- sample(sigma_sample[sigma_sample>0 & sigma_sample < 2], 1)



# do simulation and fit many times


sim_ric_model <- cmdstanr::cmdstan_model(file.path(here("simulation",
                                                     "stan_models",
                                                     "code",
                                                     "ric_simple_model_for_simulated_data.stan")))

sim_bh_model <- cmdstanr::cmdstan_model(file.path(here("simulation",
                                                    "stan_models",
                                                    "code",
                                                    "bh_simple_model_for_simulated_data.stan")))


generating_model = c("Ricker", "Beverton-Holt")
fitting_model = c("Ricker", "Beverton-Holt")

nsims <- 200
model_results_df <- data.frame(simulation = numeric(),
                               data_model = character(),
                               parameter = character(),
                               true_value = numeric(),
                               fitting_model = character(),
                               estimate_mean = numeric(),
                               estimate_lower = numeric(),
                               estimate_upper = numeric(),
                               Rhat = numeric(),
                               error = numeric())


# sample in parallel

options(mc.cores = parallel::detectCores())

# nsims <- 500
start.time <- Sys.time()

for(i in 1:nsims){
  
  print(i)
  
  set.seed(12345+i)
  
  n0_mean = 200
  
  n0 = sample(seq(0.1 * n0_mean, 4 * n0_mean), 1)
  
  min_years <- 20
  
  years <- sample(seq(min_years, 62), 1)
  
  alpha_mean = 1.5
  
  sigma_mean = 1
  
  K_max = 10000
  
  #Random K
  K = sample(seq(0.8 * K_max, K_max), 1)
  
  alpha_sample <- rnorm(100, alpha_mean, sd = 1)
  
  alpha <- sample(alpha_sample[alpha_sample > 0 & alpha_sample < 10], 1)
  
  sigma_sample <- rnorm(100, sigma_mean, 1)
  
  sigma <- sample(sigma_sample[sigma_sample > 0 & sigma_sample < 2], 1)
  
  model <- sample(generating_model, 1)
  
  if(model == "Beverton-Holt"){
    
    data <- bh_function(alpha = alpha, 
                        n0 = n0,
                        sigma = sigma,
                        K = K, 
                        years = years)
    
    data$generating_model <- "Beverton-Holt"
    
  } else{
    
    data <- ric_function(alpha = alpha,
                         K = K, 
                         n0 = n0,
                         sigma = sigma,
                         years = years)
    
    data$generating_model <- "Ricker"
    
  }
  
  data <- data %>% 
    filter(!is.nan(ln_RS), !is.infinite(ln_RS))
  
  # #if data has <2 rows, then go to next simulation
  # if(nrow(data) < 2){
  #   next
  # }
  
  
  true_values <- data %>% 
    # group_by(sigma, alpha, K, Smsy) %>% 
    summarize(sigma = mean(sigma), 
              # forestry_effect = mean(forestry_effect), 
              alpha = mean(alpha), 
              # Smax = mean(Smax),
              # Rk = mean(Rk),
              K = mean(K),
              Smsy = mean(Smsy),
              msy = mean(msy),
              Umsy = mean(Umsy),
              generating_model = first(generating_model)) %>% 
    pivot_longer(cols = c(sigma, alpha, K, Smsy, msy, Umsy), names_to = "parameter", values_to = "true_value")
  
  
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
    
    
    set.seed(12345+i)
    
    
    
    if(fit_model == "Beverton-Holt"){
      
      model_sampling <- sim_bh_model$sample(data = data_list,
                                            iter_sampling  = 1000,
                                            chains = 4,
                                            iter_warmup = 500)
      
      
      
      
      
    } else if(fit_model == "Ricker"){
      
      model_sampling <- sim_ric_model$sample(data = data_list,
                                             iter_sampling  = 1000,
                                             chains = 4,
                                             iter_warmup = 500)
      
      
      
      
      
      
    }
    
    Rhat_values <- data.frame(Rhat = round(model_sampling$summary()$rhat,3)) %>% 
      mutate(parameter = model_sampling$summary()$variable)
    
    
    model_results <- data.frame(model_sampling$draws(variables=c("alpha", "sigma", "K", "Smsy", "msy", "Umsy"),format='draws_matrix')) %>%
      mutate(fitting_model = fit_model, simulation = i) %>% 
      select(fitting_model, alpha, sigma, K, Smsy, msy, Umsy, simulation) %>% 
      pivot_longer(cols = c(alpha, K, sigma, Smsy, msy, Umsy), names_to = "parameter", values_to = "value") %>%
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
    
    
    
    
    
    model_results_df <- model_results_df %>%
      bind_rows(model_results)
  }
  
  
}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

#make a column of simulation number from 1 to i and repeat the each number 8 times

# simulation <- NULL
# for(z in 1:(i)){
#   
#   zz <- rep(z,8)
#   simulation[((z-1)*8 +1):(z*8)] <- zz
#   
# }
# model_results_df$simulation <- simulation

model_results_df_new <- model_results_df %>% 
  group_by(simulation, generating_model, fitting_model) %>%
  mutate(alpha = true_value[parameter == "alpha"],
         sigma = true_value[parameter == "sigma"],
         # b_for = true_value[parameter == "b_for"],
         Smsy = true_value[parameter == "Smsy"],
         msy = true_value[parameter == "msy"],
         Umsy = true_value[parameter == "Umsy"],
         K = true_value[parameter == "K"]
  ) %>%
  mutate(alpha_estimate = estimate_median[parameter == "alpha"],
         sigma_estimate = estimate_median[parameter == "sigma"],
         # b_for_estimate = estimate_median[parameter == "b_for"],
         Smsy_estimate = estimate_median[parameter == "Smsy"],,
         msy_estimate = estimate_median[parameter == "msy"],
         Umsy_estimate = estimate_median[parameter == "Umsy"],
         K_estimate = estimate_median[parameter == "K"]
  ) %>%
  ungroup() 


write_csv(model_results_df_new, here("simulation",
                                 "stan_models",
                                 "output",
                                 paste0("simulation_fitting_results_wo_age_w_msy_",nsims,".csv")))

pal <- PNWColors::pnw_palette("Starfish", 5)

ggplot(model_results_df_new %>% 
         filter(parameter == "alpha"), aes(x = alpha, y = estimate_median)) +
  # geom_point(aes(color = alpha),alpha = 0.5, size = 2) +
  geom_pointrange(aes(ymin= estimate_lower, ymax = estimate_upper, color = sigma), size = 0.5, alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ paste("Data model: ",generating_model) + paste("Fitting model: ",fitting_model)) +
  labs(x = "True alpha", y = "Estimated alpha") +
  # scale_color_gradient2(name = 'alpha',
  #                       low = pal[2], mid = 'gray', high = pal[4], midpoint = 5) +
  scale_color_gradientn(name = 'sigma',
                        colors = pal)+
  theme_classic()

ggplot(model_results_df_new %>% 
         filter(parameter == "K"), aes(x = K, y = estimate_median)) +
  geom_pointrange(aes(ymin= estimate_lower, ymax = estimate_upper, color = sigma), size = 0.5, alpha = 0.5)+
  # geom_point(alpha = 0.5, size = 2) +
  # geom_pointrange(aes(ymin= estimate_lower, ymax = estimate_upper, color = Smsy), size = 0.5, alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_y_continuous(labels = scales::comma, breaks = c(1000, 10000, 100000)) +
  facet_wrap(~ paste("Data model: ",generating_model) + paste("Fitting model: ",fitting_model)) +
  labs(x = "True K", y = "Estimated K") +
  coord_trans(y="log10") +
  # scale_color_gradient2(name = 'alpha',
  #                       low = pal[2], mid = 'gray', high = pal[4], midpoint = 5) +
  scale_color_gradientn(name = 'Estimated S_msy',
                        colors = pal)+
  theme_classic()

ggplot(model_results_df_new %>% 
         filter(parameter == "Smsy"), aes(x = Smsy, y = estimate_median)) +
  # geom_point(aes(color = alpha),alpha = 0.5, size = 2) +
  geom_pointrange(aes(ymin= estimate_lower, ymax = estimate_upper, color = sigma), size = 0.5, alpha = 0.5)+
  coord_trans(y="log10") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_y_continuous(labels = scales::comma, breaks = c(1000, 10000, 100000)) +
  facet_wrap(~ paste("Data model: ",generating_model) + paste("Fitting model: ",fitting_model), scale = "free_y") +
  labs(x = "True Smsy", y = "Estimated Smsy") +
  scale_color_gradientn(name = 'sigma',
                        colors = pal)+
  theme_classic()

ggplot(model_results_df_new %>% 
         filter(parameter == "msy"), aes(x = msy, y = estimate_median)) +
  # geom_point(aes(color = alpha),alpha = 0.5, size = 2) +
  geom_pointrange(aes(ymin= estimate_lower, ymax = estimate_upper, color = sigma), size = 0.5, alpha = 0.5)+
  coord_trans(y="log10") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_y_continuous(labels = scales::comma, breaks = c(1000, 10000, 100000)) +
  facet_wrap(~ paste("Data model: ",generating_model) + paste("Fitting model: ",fitting_model), scale = "free_y") +
  labs(x = "True MSY", y = "Estimated MSY") +
  scale_color_gradientn(name = 'sigma',
                        colors = pal)+
  theme_classic()

ggplot(model_results_df_new %>% 
         filter(parameter == "Umsy"), aes(x = Umsy, y = estimate_median)) +
  # geom_point(aes(color = alpha),alpha = 0.5, size = 2) +
  geom_pointrange(aes(ymin= estimate_lower, ymax = estimate_upper, color = sigma), size = 0.5, alpha = 0.5)+
  # coord_trans(y="log10") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_y_continuous(labels = scales::comma, breaks = c(1000, 10000, 100000)) +
  facet_wrap(~ paste("Data model: ",generating_model) + paste("Fitting model: ",fitting_model), scale = "free_y") +
  labs(x = "True U_MSY", y = "Estimated U_MSY") +
  scale_color_gradientn(name = 'sigma',
                        colors = pal)+
  theme_classic()




