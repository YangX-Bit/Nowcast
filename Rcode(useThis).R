setwd("/Users/xiaoy0a/Desktop/Github/Dengue/Nowcast")

# import simulation data
source("../Nowcast/functions.R")

library(rstan)
#sims data
set.seed(1)

D <- 16
N_obs <- 58
p1 <- simsP(Type = "GD", gd_alpha = 1:20, gd_beta = 20:1, D = D, days = 58)

data <- simsDataGen(p = p1, days = N_obs)

case_reported <- data$case_reported
p <- p1[,c(1:16)]


stan_data <- list(N_obs = N_obs, D = D, Y = case_reported, q = p)

fit <- stan(
  file = "stan_model_code.stan",  
  data = stan_data, 
  iter = 100, chains = 2, seed = 123
)

print(fit)

data$case_true
