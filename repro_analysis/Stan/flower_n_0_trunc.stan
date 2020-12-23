
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  vector[N] size_l;
  vector[N] clim_val;
  vector[N] native;
  
}

parameters {
  
  real b_intercept;
  real b_size;
  real b_clim;
  real b_native;
  real b_clim_nat;
  
  real<lower=0> shape;  // shape parameter
}


transformed parameters {
 
 real pred_p[N];
 
 for(it in 1:N) {
   
   pred_p[it] = b_intercept    + 
   b_size     * size_l[it]   +
   b_clim     * clim_val[it] +
   b_native   * native[it]   +
   b_clim_nat * native[it] * clim_val[it];
   
 }
 
}


model {
  // priors including all constants
  
  b_intercept ~ normal(0, 100);
  b_size ~ normal(0, 100);
  b_clim ~ normal(0, 100);
  b_native ~ normal(0, 100);
  b_clim_nat ~ normal(0, 100);
  shape ~ gamma(0.01, 0.01);
  
  // 
  for(j in 1:N) {
    
    Y[j] ~ neg_binomial_2_log(pred_p[j], shape);
    
    target += -log1m(neg_binomial_2_log_lpmf(0 | pred_p[j], shape));
  }
}
generated quantities {
}
