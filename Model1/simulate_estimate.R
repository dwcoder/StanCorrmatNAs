#!/usr/bin/env Rscript

# We are using arrays, and for that we need the most recent version of stan.
# We can't get that from cran, so we have to install from Source:
# https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-Source
# https://discourse.mc-stan.org/t/how-to-update-my-stan-version/31053/2

# remotes::install_github("stan-dev/rstan@experimental", subdir = "StanHeaders")
# remotes::install_github("stan-dev/rstan@experimental", subdir = "rstan/rstan")
# The above might ask you to update other packages as well, definitely choose "yes"
# If you get errors, make sure to install things like r-cran-curl


library(MASS)
library(rstan)

options(mc.cores = 7)
rstan_options(auto_write = TRUE)


corrmat <- matrix(
  c( 1.0, 0.3, 0.6,
     0.3, 1.0, 0.4,
     0.6, 0.4, 1.0), 3,3, byrow = TRUE)

sigma_vec <- c(10, 12, 15)
mu_vec <- c(5, 10, 12)

covmat <- diag(sigma_vec) %*% corrmat %*% diag(sigma_vec)

K <- 3
N <- 5000
y <- mvrnorm(N, mu = mu_vec, Sigma = covmat)


stanlist <- list(
    y = y,
    N = N,
    K = K)

samples <- stan(file='corrmodel1.stan', data=stanlist)

plot(samples, pars = c( 'mat_Omega[1,2]', 
                        'mat_Omega[1,3]', 
                        'mat_Omega[2,3]'
                        ))


plot(samples, pars = c( 'mu' ))
plot(samples, pars = c( 'L_sigma' ))



## Test the second file

samples <- stan(file='corrmodel_pairwise.stan', data=stanlist, chains = 1)
plot(samples, pars = c( 'mat_Omega[1,2]', 
                        'mat_Omega[1,3]', 
                        'mat_Omega[2,3]'
                        ))


plot(samples, pars = c( 'mu' ))
plot(samples, pars = c( 'L_sigma' ))

##########################################################
## Model 3: With missing values
##########################################################
N <- 500
K <- 3
y <- mvrnorm(N, mu = mu_vec, Sigma = covmat)

# Make some missing data
x1mis <- sample.int(N, size = 5)
x2mis <- sample.int(N, size = 15)
x3mis <- sample.int(N, size = 25)

y[x1mis,1] <- NA
y[x2mis,2] <- NA
y[x3mis,3] <- NA



# The easiest datastructure to use is to pass just the 
# Amount we need to stan.
# Thus, we need a way to map the K*(K-1)/2 pairs to numbers, and back


## Test the rho to matrix mapping
##
test <- function(K) {
rhos = 1:((K)*(K-1)/2)
L_Omega <- matrix(0, K,K)
for ( rrow in 2:K){
  for (ccol in 1:(rrow-1)){
    ii <-(rrow-1)*(rrow-2)/2 + ccol
    L_Omega[rrow, ccol] <-  rhos[ ii ];
  }
}
L_Omega
}
test(6)
##


N_pairs <- K*(K-1)/2 # number of rho pairs
max_N <- sort(colSums(!is.na(y)), partial=K )[K-1] # We can actually use the second largest

data_struct <- array(NA, dim=c(N_pairs, max_N, 2))

num_non_missing <- rep(NA, N_pairs)

for ( rrow in 2:K){
  for (ccol in 1:(rrow-1)){
    ii <-(rrow-1)*(rrow-2)/2 + ccol
    y_temp <- y[,c(rrow,ccol)]
    
    c_mask <- complete.cases(y_temp)
    n_temp <- sum(c_mask) 
    y_temp[c_mask,]
    
    data_struct[ii,1:n_temp,] <- y_temp[c_mask]
    num_non_missing[ii] <- n_temp
    
}
}

data_struct[1,,]
data_struct[2,,]
data_struct[3,,]

data_struct[is.na(data_struct)] <- 0

stan_data <- 
  list( K = K,
        N_pairs = N_pairs,
        max_N = max_N,
        num_non_missing= num_non_missing,
        data_struct_y = data_struct)

simple_corrmat <- cor(y, use = 'complete.obs')
L_Omega = chol(simple_corrmat)
simple_sd <- apply(y, 2, function(x) sd(x, na.rm=TRUE))
simple_mu <- apply(y, 2, function(x) mean(x, na.rm=TRUE))

inits <- list(mu= simple_mu, 
              L_sigma = simple_sd,
              Omega = simple_corrmat)
init_f <- function(chainnum) return(inits)

samples <- stan(file='corrmodel_pairwise_w_missing.stan', data=stan_data,init=init_f, chains = 1)

plot(samples, pars = c( 'Omega[1,2]', # 0.3
                        'Omega[1,3]', # 0.6
                        'Omega[2,3]'  # 0.4
                        ))

plot(samples, pars = c('mu')) # 5, 10, 15
plot(samples, pars = c('L_sigma')) #  10, 12, 15


