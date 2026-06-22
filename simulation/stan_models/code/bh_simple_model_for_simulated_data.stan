data {
  int<lower=0> N;// number of observations
  array[N] int<lower=0, upper=2025> year; //brood year
  vector[N] spawners; //spawners
  vector[N] ln_RS; //log recruits per spawner, productivity
  // vector[N] forestry;
  real Rk_mean;
  real Rk_sigma;
  real prior_alpha;
  
}

transformed data {
  real log_Rk_pr_sigma;
  real log_Rk_pr_mean;
  log_Rk_pr_sigma = sqrt(log(1+((Rk_sigma)^2)/((Rk_mean)^2)));
  log_Rk_pr_mean = log(Rk_mean) - 0.5*log_Rk_pr_sigma^2;
}

parameters {
  real<lower = 0, upper = 10> alpha;
  real Rk;
  real sigma;
  // real b_for;
}

// transformed parameters {
//   real b;
//   b = 1/Rk;
// }

model {
  vector[N] mu;
  vector[N] e_t;
  // mu = alpha - log(1 + (exp(alpha)/Rk)*spawners) + b_for*forestry;
  mu = alpha - log(1 + (exp(alpha)/Rk)*spawners);
  e_t = ln_RS - mu;
  ln_RS ~ normal(mu, sigma);
  alpha ~ normal(prior_alpha,10);
  Rk  ~ lognormal(log_Rk_pr_mean, log_Rk_pr_sigma);
  // b_for ~ normal(0,1);
  
}

generated quantities {
  // vector[N] yrep;
  // 
  // for(i in 1:N){
  //   yrep[i] = normal_rng(alpha - log(1 + (exp(alpha)/Rk)*spawners[i]) + b_for*forestry[i], sigma);
  // }
  // 
  
  real Smsy;
  real K;

  K = Rk*(exp(alpha) - 1)/exp(alpha);

  Smsy = K/(exp(alpha/2)+1);



}
