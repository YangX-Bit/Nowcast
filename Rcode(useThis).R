setwd("/Users/xiaoy0a/Desktop/Github/Dengue/Nowcast")

# import simulation data
source("../Nowcast/functions.R")

library(rstan)
#sims data
set.seed(1)

D <- 15
N_obs <- 58
p1 <- simsP(Type = "GD", gd_alpha = seq(3,3.3, by = 0.01), gd_beta =  rep(30,30), 
            D = D, days = 58)


data <- simsDataGenQ( days = N_obs)

# sims
# avg_lambda_t <- apply(data$lambda, 1, mean)  
# avg_true_cases <- apply(data$case_true, 1, mean) 
# 
# case_reported_acc <- data$case_reported_acc
# case_reported_acc[,16]
data$case_reported
data$case_true
data$lambda

# p_to_q <- function(p){
#   # convert report prob to acc-report prob
#   len <- length(p)
#   q <- numeric(len)
#   for (i in 1:len) {
#     q[i] <- sum(head(p,i))
#   }
#   q
# }
# 
# q <- p_to_q(p)

stan_data <- list(N_obs = N_obs, D = D + 1, Y = round(data$case_reported))

fit <- stan(
  file = "stan_model_code.stan",  
  data = stan_data, 
  iter = 2000, chains = 3, seed = 123
)

print(fit)

# posterior of lambda_t and b
lambda_samples <- rstan::extract(fit)$lambda_t
b_samples <- rstan::extract(fit)$b_t

lambda_samples[,1]
rpois(3000,lambda_samples[,1])

hist(rpois(3000,lambda_samples[,1]))
# distribution for Nt
# N_t_posterior <- function(lambda_t_sample, b_sample, D) {
#   N_t_estimates <- numeric(length(lambda_t_sample))
#   for (i in 1:length(lambda_t_sample)) {
#     q_d <- 1 - exp(-1 * b_sample[i] * D)
#     N_t_estimates[i] <- sum(rpois(1, lambda_t_sample[i] * q_d))  # Nt
#   }
#   return(N_t_estimates)
# }
# 
# N_t_estimates <- N_t_posterior(lambda_samples, b_samples, D)
# hist(N_t_estimates) 

# Function to calculate the Maximum A Posteriori (MAP) estimate of N_t

N_t_posterior_map <- function(lambda_t_samples, b_samples, D) {
  sampling_size <- nrow(b_samples)
  time_length <- ncol(b_samples)
  N_t_estimates <- matrix(0, nrow = sampling_size, ncol = time_length)
  
  q_d_map <- numeric(sampling_size)
  for (t in 1:time_length) {
    # Calculate the MAP estimate of q_d for each time point using the median of b_samples
    q_d_map <- 1 - exp(- b_samples[,t] * D)
    
    # Calculate the MAP estimate of N_t as the product of the median of lambda_t and q_d_map
    N_t_estimates[,t] <- rpois(sampling_size,lambda_t_samples[,t])

  }
  
  return(N_t_estimates)  # Return the MAP estimates of N_t for each time point
}
N_t_map_estimates <- N_t_posterior_map(lambda_samples, b_samples, D)

round((apply(N_t_map_estimates, 2, mean) - data$case_true),0) 


################# try p2
data_2 <- simsDataGen(p = p2, days = N_obs)
stan_data_2 <- list(N_obs = N_obs, D = D + 1, Y = round(data_2$case_reported_acc))

fit_2 <- stan(
  file = "stan_model_code.stan",  
  data = stan_data_2, 
  iter = 2000, chains = 3, seed = 123
)

lambda_samples_2 <- rstan::extract(fit_2)$lambda_t
b_samples_2 <- rstan::extract(fit_2)$b_t

N_t_map_estimates_2 <- N_t_posterior_map(lambda_samples_2, b_samples_2, D)

round(apply(N_t_map_estimates_2, 2, mean) - apply(data_2$case_true, 1, mean),2) 




################# try p3
data_3 <- simsDataGen(p = p3, days = N_obs)
stan_data_3 <- list(N_obs = N_obs, D = D + 1, Y = round(data_3$case_reported_acc))

fit_3 <- stan(
  file = "stan_model_code.stan",  
  data = stan_data_3, 
  iter = 2000, chains = 3, seed = 123
)

lambda_samples_3 <- rstan::extract(fit_3)$lambda_t
b_samples_3 <- rstan::extract(fit_3)$b_t

N_t_map_estimates_3 <- N_t_posterior_map(lambda_samples_3, b_samples_3, D)

round(apply(N_t_map_estimates_3, 2, mean) - apply(data_3$case_true, 1, mean),2) 

