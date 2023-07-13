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
  
  // Temporary variables for the pairwise loop:
  matrix[2, 2] L_Omega_temp;
  matrix[K, K] L_Sigma_temp;
  vector[2] mu_temp;
  vector[2] sigma_temp;
  
  array[N] vector[2] y_temp;
  
  
  for ( rrow in 2:K){
  for (ccol in 1:(rrow-1)){
    y_temp[:,1] = y[:,rrow];
    y_temp[:,2] = y[:,ccol];
    
    mu_temp[1] = mu[rrow];
    mu_temp[2] = mu[ccol];
    
    sigma_temp[1] = L_sigma[rrow];
    sigma_temp[2] = L_sigma[ccol];
    L_Omega_temp[1,1] = L_Omega[rrow, rrow];
    L_Omega_temp[2,2] = L_Omega[ccol, ccol];
    L_Omega_temp[2,1] = L_Omega[rrow, ccol];
     
    L_Sigma_temp = diag_pre_multiply(sigma_temp, L_Omega_temp);
    y_temp ~ multi_normal_cholesky(mu_temp, L_Sigma_temp);
  }}
  

  // Priors
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);
}
generated quantities {
  matrix[K, K] mat_Omega;
  mat_Omega = L_Omega * L_Omega';
}

