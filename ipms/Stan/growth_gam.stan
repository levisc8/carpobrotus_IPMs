// generated with brms 2.16.3
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for spline t2(seas_temp_t, log_size, bs = "cs", k = 4)
  int nb_1;  // number of bases
  int knots_1[nb_1];  // number of knots
  // basis function matrices
  matrix[N, knots_1[1]] Zs_1_1;
  // data for spline t2(total_prec_t, log_size, bs = "cs", k = 4)
  int nb_2;  // number of bases
  int knots_2[nb_2];  // number of knots
  // basis function matrices
  matrix[N, knots_2[1]] Zs_2_1;
  // data for spline t2(seas_prec_t, log_size, bs = "cs", k = 4)
  int nb_3;  // number of bases
  int knots_3[nb_3];  // number of knots
  // basis function matrices
  matrix[N, knots_3[1]] Zs_3_1;
  // data for spline s(mean_sw2_t, bs = "cs", k = 4)
  int nb_4;  // number of bases
  int knots_4[nb_4];  // number of knots
  // basis function matrices
  matrix[N, knots_4[1]] Zs_4_1;
  // data for spline s(mean_sw1_t, bs = "cs", k = 4)
  int nb_5;  // number of bases
  int knots_5[nb_5];  // number of knots
  // basis function matrices
  matrix[N, knots_5[1]] Zs_5_1;
  int<lower=1> K_sigma;  // number of population-level effects
  matrix[N, K_sigma] X_sigma;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_sigma_1;
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
  // parameters for spline t2(seas_temp_t, log_size, bs = "cs", k = 4)
  // standarized spline coefficients
  vector[knots_1[1]] zs_1_1;
  real<lower=0> sds_1_1;  // standard deviations of spline coefficients
  // parameters for spline t2(total_prec_t, log_size, bs = "cs", k = 4)
  // standarized spline coefficients
  vector[knots_2[1]] zs_2_1;
  real<lower=0> sds_2_1;  // standard deviations of spline coefficients
  // parameters for spline t2(seas_prec_t, log_size, bs = "cs", k = 4)
  // standarized spline coefficients
  vector[knots_3[1]] zs_3_1;
  real<lower=0> sds_3_1;  // standard deviations of spline coefficients
  // parameters for spline s(mean_sw2_t, bs = "cs", k = 4)
  // standarized spline coefficients
  vector[knots_4[1]] zs_4_1;
  real<lower=0> sds_4_1;  // standard deviations of spline coefficients
  // parameters for spline s(mean_sw1_t, bs = "cs", k = 4)
  // standarized spline coefficients
  vector[knots_5[1]] zs_5_1;
  real<lower=0> sds_5_1;  // standard deviations of spline coefficients
  vector[Kc_sigma] b_sigma;  // population-level effects
  real Intercept_sigma;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
}
transformed parameters {
  // actual spline coefficients
  vector[knots_1[1]] s_1_1;
  // actual spline coefficients
  vector[knots_2[1]] s_2_1;
  // actual spline coefficients
  vector[knots_3[1]] s_3_1;
  // actual spline coefficients
  vector[knots_4[1]] s_4_1;
  // actual spline coefficients
  vector[knots_5[1]] s_5_1;
  vector[N_1] r_1_1;  // actual group-level effects
  vector[N_2] r_2_sigma_1;  // actual group-level effects
  // compute actual spline coefficients
  s_1_1 = sds_1_1 * zs_1_1;
  // compute actual spline coefficients
  s_2_1 = sds_2_1 * zs_2_1;
  // compute actual spline coefficients
  s_3_1 = sds_3_1 * zs_3_1;
  // compute actual spline coefficients
  s_4_1 = sds_4_1 * zs_4_1;
  // compute actual spline coefficients
  s_5_1 = sds_5_1 * zs_5_1;
  r_1_1 = (sd_1[1] * (z_1[1]));
  r_2_sigma_1 = (sd_2[1] * (z_2[1]));
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + Xc * b + Zs_1_1 * s_1_1 + Zs_2_1 * s_2_1 + Zs_3_1 * s_3_1 + Zs_4_1 * s_4_1 + Zs_5_1 * s_5_1;
    // initialize linear predictor term
    vector[N] sigma = Intercept_sigma + Xc_sigma * b_sigma;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      sigma[n] += r_2_sigma_1[J_2[n]] * Z_2_sigma_1[n];
    }
    for (n in 1:N) {
      // apply the inverse link function
      sigma[n] = exp(sigma[n]);
    }
    target += normal_lpdf(Y | mu, sigma);
  }
  // priors including constants
  target += student_t_lpdf(Intercept | 3, -2.3, 2.5);
  target += student_t_lpdf(sds_1_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_1_1);
  target += student_t_lpdf(sds_2_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_2_1);
  target += student_t_lpdf(sds_3_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_3_1);
  target += student_t_lpdf(sds_4_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_4_1);
  target += student_t_lpdf(sds_5_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(zs_5_1);
  target += student_t_lpdf(Intercept_sigma | 3, 0, 2.5);
  target += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_1[1]);
  target += student_t_lpdf(sd_2 | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  target += std_normal_lpdf(z_2[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_sigma_Intercept = Intercept_sigma - dot_product(means_X_sigma, b_sigma);
}
