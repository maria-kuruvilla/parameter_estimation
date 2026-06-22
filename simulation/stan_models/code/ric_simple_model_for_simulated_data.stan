data {
  int<lower=0> N;// number of observations
  array[N] int<lower=0, upper=2025> year; //brood year
  vector[N] spawners; //spawners
  vector[N] ln_RS; //log recruits per spawner, productivity
  // vector[N] forestry;
  real Smax_mean;
  real Smax_sigma;
  real prior_alpha;
  
}

transformed data {
  real log_Smax_pr_sigma;
  real log_Smax_pr_mean;
  log_Smax_pr_sigma = sqrt(log(1+((Smax_sigma)^2)/((Smax_mean)^2)));
  log_Smax_pr_mean = log(Smax_mean) - 0.5*log_Smax_pr_sigma^2;
}

parameters {
  real<lower = 0, upper = 10> alpha;
  real Smax;
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
  // mu = alpha - spawners/Smax + b_for*forestry;
  mu = alpha - spawners/Smax;
  e_t = ln_RS - mu;
  ln_RS ~ normal(mu, sigma);
  alpha ~ normal(prior_alpha,10);
  Smax  ~ lognormal(log_Smax_pr_mean, log_Smax_pr_sigma);
  // b_for ~ normal(0,1);
  
}

generated quantities {
  // vector[N] yrep;
  // 
  // for(i in 1:N){
  //   yrep[i] = normal_rng(alpha - spawners[i]/Smax + b_for*forestry[i], sigma);
  // }
  real Smsy;
  real z;
  real w;
  real K;
  
  K = Smax*alpha;

  z = exp(1.0-alpha);

  w = lambert_w0(exp(1 - alpha));

  Smsy = (1-w)*Smax;

  
}
