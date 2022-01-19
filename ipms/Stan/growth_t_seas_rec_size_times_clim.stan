// generated with brms 2.16.3
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_sigma;  // number of population-level effects
  matrix[N, K_sigma] X_sigma;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  int Kc_sigma = K_sigma - 1;
  matrix[N, Kc_sigma] Xc_sigma;  // centered version of X_sigma without an intercept
  vector[Kc_sigma] means_X_sigma;  // column means of X_sigma before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 2:K_sigma) {
    means_X_sigma[i - 1] = mean(X_sigma[, i]);
    Xc_sigma[, i - 1] = X_sigma[, i] - means_X_sigma[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector[Kc_sigma] b_sigma;  // population-level effects
  real Intercept_sigma;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  r_1_1 = (sd_1[1] * (z_1[1]));
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b;
    // initialize linear predictor term
    vector[N] sigma = Intercept_sigma + Xc_sigma * b_sigma;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    }
    for (n in 1:N) {
      // apply the inverse link function
      sigma[n] = exp(sigma[n]);
    }
    target += normal_lpdf(Y | mu, sigma);
  }
  // priors including constants
  target += student_t_lpdf(Intercept | 3, -2.4, 2.5);
  target += student_t_lpdf(Intercept_sigma | 3, 0, 2.5);
  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_sigma_Intercept = Intercept_sigma - dot_product(means_X_sigma, b_sigma);
}
