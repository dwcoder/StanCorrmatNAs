data {
  int<lower=1> K; // The number of y_s
  int<lower=0> N_pairs; // The number of pairs = K * (K -1) /2
  int<lower=0> max_N; // The maximum size of the struct we will keep the data in 
  
  array[N_pairs] int<lower=0> num_non_missing;
  
  array[N_pairs, max_N] vector[2] data_struct_y;
  }
parameters {
  corr_matrix[K] Omega; 
  vector<lower=0>[K] L_sigma;
  vector[K] mu;
}
model {
  matrix[K, K] L_Sigma;
  
  // Temporary variables for the pairwise loop:
  matrix[2, 2] Omega_temp;
  matrix[2, 2] Sigma_temp;
  vector[2] mu_temp;
  vector[2] sigma_temp;
  int ii;
  
  
  for ( rrow in 2:K){
  for (ccol in 1:(rrow-1)){
    ii = (rrow-1)*(rrow-2)/2 + ccol; # Magic formula to convert row/col into a pairnumber
    
    array[num_non_missing[ii]] vector[2] y_temp;

    y_temp = data_struct_y[ii, 1:num_non_missing[ii]]; # will have shape y[n, 2]
    
    mu_temp[1] = mu[rrow];
    mu_temp[2] = mu[ccol];
    
    sigma_temp[1] = L_sigma[rrow];
    sigma_temp[2] = L_sigma[ccol];
    Omega_temp[1,1] = 1.0;
    Omega_temp[2,2] = 1.0;
    Omega_temp[2,1] = Omega[rrow, ccol]; // rho
    Omega_temp[1,2] = Omega[ccol, rrow];
     
    Sigma_temp = quad_form_diag(Omega_temp, sigma_temp);
    y_temp ~ multi_normal(mu_temp, Sigma_temp);
  }}
  

  // Priors
  Omega ~ lkj_corr(1);
  L_sigma ~ cauchy(0, 10);
}
