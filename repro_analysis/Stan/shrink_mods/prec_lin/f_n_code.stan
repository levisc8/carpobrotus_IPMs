// generated with brms 2.12.0
functions {

  /* Efficient computation of the horseshoe prior
   * Args:
   *   zb: standardized population-level coefficients
   *   global: global horseshoe parameters
   *   local: local horseshoe parameters
   *   scale_global: global scale of the horseshoe prior
   *   c2: positive real number for regularization
   * Returns:
   *   population-level coefficients following the horseshoe prior
   */
  vector horseshoe(vector zb, vector[] local, real[] global,
                   real scale_global, real c2) {
    int K = rows(zb);
    vector[K] lambda = local[1] .* sqrt(local[2]);
    vector[K] lambda2 = square(lambda);
    real tau = global[1] * sqrt(global[2]) * scale_global;
    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
    return zb .* lambda_tilde * tau;
  }
}
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  int lb[N];  // lower truncation bounds;
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for the horseshoe prior
  real<lower=0> hs_df;  // local degrees of freedom
  real<lower=0> hs_df_global;  // global degrees of freedom
  real<lower=0> hs_df_slab;  // slab degrees of freedom
  real<lower=0> hs_scale_global;  // global prior scale
  real<lower=0> hs_scale_slab;  // slab prior scale
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
  // local parameters for horseshoe prior
  vector[Kc] zb;
  vector<lower=0>[Kc] hs_local[2];
  real Intercept;  // temporary intercept for centered predictors
  // horseshoe shrinkage parameters
  real<lower=0> hs_global[2];  // global shrinkage parameters
  real<lower=0> hs_c2;  // slab regularization parameter
}
transformed parameters {
  vector[Kc] b;  // population-level effects
  // compute actual regression coefficients
  b = horseshoe(zb, hs_local, hs_global, hs_scale_global, hs_scale_slab^2 * hs_c2);
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + Xc * b;
  for (n in 1:N) {
    // apply the inverse link function
    mu[n] = exp(mu[n]);
  }
  // priors including all constants
  target += normal_lpdf(zb | 0, 1);
  target += normal_lpdf(hs_local[1] | 0, 1)
    - 37 * log(0.5);
  target += inv_gamma_lpdf(hs_local[2] | 0.5 * hs_df, 0.5 * hs_df);
  target += student_t_lpdf(Intercept | 3, 1, 10);
  target += normal_lpdf(hs_global[1] | 0, 1)
    - 1 * log(0.5);
  target += inv_gamma_lpdf(hs_global[2] | 0.5 * hs_df_global, 0.5 * hs_df_global);
  target += inv_gamma_lpdf(hs_c2 | 0.5 * hs_df_slab, 0.5 * hs_df_slab);
  // likelihood including all constants
  if (!prior_only) {
    for (n in 1:N) {
      target += poisson_lpmf(Y[n] | mu[n]) -
        poisson_lccdf(lb[n] - 1 | mu[n]);
    }
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
