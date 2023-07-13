data {
  int<lower=1> K;
  int<lower=0> N; 
  array[N] vector[K] y;
  }
parameters {
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] L_sigma;
  vector[K] mu;
}
model {
  matrix[K, K] L_Sigma;

  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);

  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);

  y ~ multi_normal_cholesky(mu, L_Sigma);
}
generated quantities {
  matrix[K, K] mat_Omega;
  mat_Omega = L_Omega * L_Omega';
}

