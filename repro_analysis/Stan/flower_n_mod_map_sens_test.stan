// generated with brms 2.11.1
functions {
}
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  // temporary intercept for centered predictors
  real Intercept;
  real<lower=0> shape;  // shape parameter
}
transformed parameters {
}
model {
  // priors including all constants
  target += student_t_lpdf(Intercept | 3, 1, 10);
  target += gamma_lpdf(shape | 0.01, 0.01);
  // likelihood including all constants
  if (!prior_only) {
    target += neg_binomial_2_log_glm_lpmf(Y | Xc, Intercept, b, shape);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
