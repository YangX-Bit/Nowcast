library(dplyr)
library(tidyr)

### functions to generate Generalized Dirichlet distribution
# Parameters:
#  n - number of samples
#  alpha - alpha
#  beta - beta
###
rGeneralizedDirichlet <- function(n = 1, alpha, beta) {
  # n is the number of samples to generate
  # alpha and beta are the parameter vectors for the Generalized Dirichlet distribution
  k <- length(alpha)  # The dimension of the vector
  X <- matrix(0, n, k)  # Initialize a matrix to store the results
  
  # Generate Beta-distributed random numbers sequentially according to the definition
  for (i in 1:n) {
    x <- numeric(k)  # Store the k components of each sample
    x[1] <- rbeta(1, alpha[1], beta[1])  # Generate X_1
    prod_term <- 1 - x[1]  # The remaining product term
    for (j in 2:(k - 1)) {
      x[j] <- rbeta(1, alpha[j], beta[j]) * prod_term  # Generate X_j
      prod_term <- prod_term - x[j]  # Update the remaining product term
    }
    x[k] <- prod_term  # The last component
    X[i, ] <- x  # Store the generated sample
  }
  
  return(X)
}

### functions to generate p(delay probability)
# Parameters:
#  n - number of samples
#  alpha - to control the rbeta(). If alpha is large, than the value is large.
#  beta - to control the rbeta(). If beta is large, than the value is large.
###
generate_controlled_sum <- function(n, alpha, beta, order = "increasing") {
  # n random numbers
  values <- rbeta(n, alpha, beta)
  
  # normalize to  1
  normalized_values <- values / sum(values)
  
  # Arrange order
  if (order == "increasing") {
    normalized_values <- sort(normalized_values)
  } else if (order == "decreasing") {
    normalized_values <- sort(normalized_values, decreasing = TRUE)
  } 
  
  return(normalized_values)
}

### functions to split n into three shares
# Parameters:
#  n - number of samples
###
split_into_three <- function(n) {
  base <- n %/% 3
  remainder <- n %% 3
  return(c(rep(base + 1, remainder), rep(base, 3 - remainder)))
}

### functions to generate p(delay probability)
# Parameters:
#  n - number of samples
#  alpha - alpha
#  beta - beta
###
simsP <- function(Type = "GD", gd_alpha = 1:20, gd_beta = 20:1, D = 15, days = 30) {
  Type <- match.arg(Type, c("GD", "Multi_GD", "basin"))
  
  if (Type == "GD") {
    # fixed prob
    p <- rGeneralizedDirichlet(1, gd_alpha, gd_beta)
    
  } else if (Type == "Multi_GD") {
    # multi
    p <- rGeneralizedDirichlet(days, gd_alpha, gd_beta)
    
  } else {  # Type == "basin"
    # spit into three parts
    n <- split_into_three(days)
    p <- matrix(0, nrow = days, ncol = D + 1)
    
    for (i in 1:days) {
      if (i <= n[1]) {
        p[i, ] <- generate_controlled_sum(D + 1, alpha = 1, beta = 1000, order = "decreasing")
        
      } else if (i <= (n[1] + n[2])) {
        p[i, ] <- generate_controlled_sum(D + 1, alpha = 100, beta = 1)
        
      } else {
        p[i, ] <- generate_controlled_sum(D + 1, alpha = 1, beta = 1000, order = "increasing")
      }
    }
  }
  
  return(p)
}

p1 <- simsP(Type = "GD", gd_alpha = 1:20, gd_beta = 20:1, D = 15, days = 58)
p2 <- simsP(Type = "Multi_GD", gd_alpha = 1:20, gd_beta = 20:1, D = 15, days = 58)
p3 <- simsP(Type = "basin", D = 15, days = 58)

par(mfrow = c(5,6))
for (i in 1:30) {
  plot(p3[i,], type = "l")
}
# 
for (i in 1:30) {
  plot(t(apply(p3, 1, cumsum))[i,], type = "l")
}

### functions to generate simulation data
# Parameters:
#  n - number of samples
#  alpha - alpha
#  beta - beta
###
simsDataGen <- function(alpha =  c(1:10, seq(10, 120, by = 4), seq(120, 3, by = -6) ), beta = 0.5,
                        days = 30, D = 15, seed = 123, p = NULL){
  if(length(alpha) < days){
    stop("Error! The length of alpha cannot be less than days!")
  }
  
  set.seed(seed)
  # Generate the probability of case being reported at a certain day
  # p <- rGeneralizedDirichlet(1, gd_alpha, gd_beta)
  # Cut the probability according to the Maximal delay that we concern
  p_cut <- p[, c(1:(D+1))]
  # Initialize storage for true cases and reported cases
  true_cases <- numeric(days)  # Actual number of cases per day
  reported_cases <- matrix(0, nrow = days, ncol = D + 1)  # Reported cases matrix for delays 0 to D days
  
  # Simulate the true number of cases per day
  for (t in 1:days) {
    # Draw the Poisson intensity parameter lambda_t from a Gamma distribution
    lambda_t <- rgamma(1, shape = alpha[t], rate = beta)
    
    # Draw the actual number of cases N(t, ∞) from a Poisson distribution
    true_cases[t] <- rpois(1, lambda = lambda_t)
  }

  # REPORTED CASES
  N_tT <- matrix(0, ncol = D+1, nrow = days)
  if(nrow(p) == 1){
    for (i in 1:days) {
      N_tT[i,] = rmultinom(1, size = true_cases[i], prob = p_cut)
    }
  }else{
    for (i in 1:days) {
      N_tT[i,] = rmultinom(1, size = true_cases[i], prob = p_cut[i,])
    }
  }
  return(N_tT)
}
simsDataGen(p = p3, days = 58)

### functions to transfer the simulation data to the form of data in the paper
# Parameters:
#  data - Matrix of cases in each day with delays
###
dataTransform <- function(data, 
                          start_date = as.Date("2011-07-04")){
  
  # sequence of the start date
  admission_dates <- start_date - 0:(nrow(data) - 1)
  
  
  df <- as.data.frame(data)
  D <- ncol(data) - 1
  colnames(df) <- paste0("delay", 0:D)
  
  # long data
  long_df <- df %>%
    mutate(admission_date = admission_dates) %>%
    pivot_longer(cols = starts_with("delay"), 
                 names_to = "delay", 
                 names_prefix = "delay",
                 values_to = "reported_cases") %>%
    filter(reported_cases > 0) %>%
    mutate(delay = as.numeric(delay),
           report_date = admission_date + delay) %>%
    uncount(reported_cases)
  
  data_out <- long_df %>% mutate( dHosp= admission_date,
                                     dReport= report_date) %>%
    select(dHosp, dReport) %>%
    as.data.frame()
  
  return(data_out)
}

