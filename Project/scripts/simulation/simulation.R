setwd("/Users/xiaoy0a/Desktop/Github/Dengue/Nowcast/Project")

# import simulation function
source("../Project/scripts/simulation/simulations_functions.R")

#sims data
set.seed(1)

# setting: delay,  number of days
D <- 15
N_obs <- 58

### scenario 1 ###
# constant b = 3
data_const <- simsDataGenQ(days = N_obs, method = "constant", b = 3)


### scenario 2 ###
# time_varying
data_time_v_1 <- simsDataGenQ(days = N_obs, method = "time_varying",
                              beta0 = -2, beta1 = 0.5, num_sims = 1)


data_time_v_1$case_reported

### write data ###

#scenario 1
saveRDS(data_const, "../Project/data/simulated/constant.Rdata")
#readRDS("../Project/data/simulated/constant.Rdata")

#scenario 2
saveRDS(data_const, "../Project/data/simulated/time_varying_1.Rdata")


beta0 = -2
beta1 = 0.01

seq_delay = 0:15
seq_time = 1:58

b_t = exp(beta0 + beta1 * seq_time)
plot(seq_time, b_t)

q_td = 1 - exp(-b_t[40] * seq_delay)
plot(seq_delay, q_td)
