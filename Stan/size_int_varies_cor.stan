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
  vector[N] Z_1_2;
  int<lower=1> NC_1;
  // data for group-level effects of ID 2
  int<lower=1> J_2[N];
  int<lower=1> N_2;
  int<lower=1> M_2;
  vector[N] Z_2_1;
  vector[N] Z_2_2;
  int<lower=1> NC_2;
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
  matrix[M_1, N_1] z_1;  // unscaled group-level effects
  // cholesky factor of correlation matrix
  cholesky_factor_corr[M_1] L_1;
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  matrix[M_2, N_2] z_2;  // unscaled group-level effects
  // cholesky factor of correlation matrix
  cholesky_factor_corr[M_2] L_2;
} 
transformed parameters { 
  // group-level effects 
  matrix[N_1, M_1] r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)';
  vector[N_1] r_1_1 = r_1[, 1];
  vector[N_1] r_1_2 = r_1[, 2];
  // group-level effects 
  matrix[N_2, M_2] r_2 = (diag_pre_multiply(sd_2, L_2) * z_2)';
  vector[N_2] r_2_1 = r_2[, 1];
  vector[N_2] r_2_2 = r_2[, 2];
} 
model { 
  vector[N] mu = temp_Intercept + Xc * b;
  for (n in 1:N) { 
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_2_2[J_2[n]] * Z_2_2[n];
  } 
  // priors including all constants 
  target += student_t_lpdf(temp_Intercept | 3, 1, 10); 
  target += gamma_lpdf(shape | 0.01, 0.01); 
  target += student_t_lpdf(sd_1 | 3, 0, 10)
    - 2 * student_t_lccdf(0 | 3, 0, 10); 
  target += normal_lpdf(to_vector(z_1) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_1 | 1); 
  target += student_t_lpdf(sd_2 | 3, 0, 10)
    - 2 * student_t_lccdf(0 | 3, 0, 10); 
  target += normal_lpdf(to_vector(z_2) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_2 | 1); 
  // likelihood including all constants 
  if (!prior_only) { 
    target += neg_binomial_2_log_lpmf(Y | mu, shape);
  } 
} 
generated quantities { 
  // actual population-level intercept 
  real b_Intercept = temp_Intercept - dot_product(means_X, b); 
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  corr_matrix[M_2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
  vector<lower=-1,upper=1>[NC_2] cor_2;
  // take only relevant parts of correlation matrix
  cor_1[1] = Cor_1[1,2]; 
  // take only relevant parts of correlation matrix
  cor_2[1] = Cor_2[1,2]; 
} 
