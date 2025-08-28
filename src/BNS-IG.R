# This file contains the function that simulate log-returns from BNS-IG model.

# --- Helper Function to Simulate the IG-OU BDLP Increment ---
simulate_ig_ou_bdlp_increment <- function(n, dt, a_ig, b_ig) {
  
  # --- Component 1: An IG-Lévy process z^(1) ---
  A <- (a_ig / 2) * dt
  B <- b_ig
  
  # Convert to the (nu, lambda) parameterization for rinvGauss
  nu_param <- A / B
  lambda_param <- A^2
  
  z1_increments <- SuppDists::rinvGauss(
    n, 
    nu = nu_param, 
    lambda = lambda_param
  )
  
  # --- Component 2: A compound Poisson process z^(2) ---
  jump_intensity_z2 <- 1/(a_ig * b_ig / 2) * dt
  num_jumps_z2 <- rpois(n, jump_intensity_z2)
  
  # Jump sizes are b_ig^-2 * v_n^2, where v_n are standard normal
  jump_sizes_z2 <- sapply(num_jumps_z2, function(nj) {
    if (nj == 0) return(0)
    sum((1 / b_ig^2) * (rnorm(nj)^2))
  })
  
  return(z1_increments + jump_sizes_z2)
}

# --- BNS Simulation Function with IG Stochastic Volatility ---
simulate_bns_ig_sv <- function(params, n_steps, dt) {
  # Unpack parameters
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3] 
  a_ig      <- params[4] 
  b_ig      <- params[5] 
  sigma0_sq <- params[6] 
  
  # Initialize vectors for the simulation path
  log_S <- numeric(n_steps + 1)
  sigma_sq <- numeric(n_steps + 1)
  log_S[1] <- 0
  sigma_sq[1] <- sigma0_sq
  
  # --- Simulate all Background Driving Lévy Process (BDLP) increments ---
  bdlp_increments <- simulate_ig_ou_bdlp_increment(
    n = n_steps, 
    dt = lambda * dt, 
    a_ig = a_ig, 
    b_ig = b_ig
  )
  
  # Generate standard normal random variables for the main Brownian motion
  dW <- rnorm(n_steps, mean = 0, sd = sqrt(dt))
  
  # --- Euler-Maruyama Discretization using the BDLP ---
  for (t in 1:n_steps) {
    sigma_sq_prev <- max(1e-9, sigma_sq[t])
    
    # Update variance process (IG-OU)
    # d(sigma_t^2) = -lambda * sigma_t^2 * dt + d(z_{lambda*t})
    sigma_sq[t+1] <- sigma_sq[t] - lambda * sigma_sq_prev * dt + bdlp_increments[t]
    
    # Update log-return process (BNS)
    # d(log S_t) = (mu - 0.5*sigma_t^2)dt + sigma_t*dW_t + rho*d(z_{lambda*t})
    log_S[t+1] <- log_S[t] + (mu - 0.5 * sigma_sq_prev) * dt + 
      sqrt(sigma_sq_prev) * dW[t] + 
      rho * bdlp_increments[t]
  }
  
  # Return the simulated log-returns
  return(diff(log_S))
}