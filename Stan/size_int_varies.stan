// generated with brms 2.7.0
functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 
  int Y[N];  // response variable 
  int<lower=1> K;  // number of population-level effects 
  matrix[N, K] X;  // population-level design matrix 
  // data for group-level effects of ID 1
  int<lower=1> J_1[N];
  int<lower=1> N_1;
  int<lower=1> M_1;
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> J_2[N];
  int<lower=1> N_2;
  int<lower=1> M_2;
  vector[N] Z_2_1;
  // data for group-level effects of ID 3
  int<lower=1> J_3[N];
  int<lower=1> N_3;
  int<lower=1> M_3;
  vector[N] Z_3_1;
  // data for group-level effects of ID 4
  int<lower=1> J_4[N];
  int<lower=1> N_4;
  int<lower=1> M_4;
  vector[N] Z_4_1;
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
  int Kc = K - 1; 
  matrix[N, K - 1] Xc;  // centered version of X 
  vector[K - 1] means_X;  // column means of X before centering 
  for (i in 2:K) { 
    means_X[i - 1] = mean(X[, i]); 
    Xc[, i - 1] = X[, i] - means_X[i - 1]; 
  } 
} 
parameters { 
  vector[Kc] b;  // population-level effects 
  real temp_Intercept;  // temporary intercept 
  real<lower=0> shape;  // shape parameter 
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // unscaled group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // unscaled group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  vector[N_3] z_3[M_3];  // unscaled group-level effects
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  vector[N_4] z_4[M_4];  // unscaled group-level effects
} 
transformed parameters { 
  // group-level effects 
  vector[N_1] r_1_1 = sd_1[1] * (z_1[1]);
  // group-level effects 
  vector[N_2] r_2_1 = sd_2[1] * (z_2[1]);
  // group-level effects 
  vector[N_3] r_3_1 = sd_3[1] * (z_3[1]);
  // group-level effects 
  vector[N_4] r_4_1 = sd_4[1] * (z_4[1]);
} 
model { 
  vector[N] mu = temp_Intercept + Xc * b;
  for (n in 1:N) { 
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_3_1[J_3[n]] * Z_3_1[n] + r_4_1[J_4[n]] * Z_4_1[n];
  } 
  // priors including all constants 
  target += student_t_lpdf(temp_Intercept | 3, 1, 10); 
  target += gamma_lpdf(shape | 0.01, 0.01); 
  target += student_t_lpdf(sd_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += normal_lpdf(z_1[1] | 0, 1);
  target += student_t_lpdf(sd_2 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += normal_lpdf(z_2[1] | 0, 1);
  target += student_t_lpdf(sd_3 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += normal_lpdf(z_3[1] | 0, 1);
  target += student_t_lpdf(sd_4 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10); 
  target += normal_lpdf(z_4[1] | 0, 1);
  // likelihood including all constants 
  if (!prior_only) { 
    target += neg_binomial_2_log_lpmf(Y | mu, shape);
  } 
} 
generated quantities { 
  // actual population-level intercept 
  real b_Intercept = temp_Intercept - dot_product(means_X, b); 
} 
