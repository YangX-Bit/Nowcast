library(rstan)

setwd("/Users/xiaoy0a/Desktop/Github/Dengue/Nowcast/Project")

#read data
# scenario 1
data_const <- readRDS("../Project/data/simulated/constant.Rdata")
# scenario 2
data_tv_1 <- readRDS("../Project/data/simulated/time_varying_1.Rdata")

# seting
D <- ncol(data_const$case_reported) - 1 # if D=3, we consder 0,1,2,3, so the data length is 4
N_obs <- nrow(data_const$case_reported)



##### fit scenario 1 #####
stan_data_1 <- list(N_obs = N_obs, D = D + 1, Y = round(data_const$case_reported))

fit_1 <- stan(
  file = "../Project/source/stan_model_constantB.stan",  
  data = stan_data_1, 
  iter = 2000, chains = 3, seed = 123
)

print(fit_1)


##### fit scenario 1 #####
# b0 = -2; b1 = 0.5
stan_data_2 <- list(N_obs = N_obs, D = D + 1, Y = round(data_tv_1$case_reported))

fit_2 <- stan(
  file = "../Project/source/stan_model_b_varying.stan",  
  data = stan_data_2, 
  iter = 3000, chains = 3, seed = 123
)

print(fit_2)

