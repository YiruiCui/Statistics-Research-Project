# This file contains the function that simulate log-returns from BNS-Gamma model.

# --- BNS Simulation Function ---
simulate_bns_gamma_sv <- function(params, n_steps, dt) {
  # Unpack parameters from the input vector
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3]
  a         <- params[4]
  b         <- params[5]
  sigma0_sq <- params[6]
  
  # Initialize vectors for the simulation path
  log_S <- numeric(n_steps + 1)
  sigma_sq <- numeric(n_steps + 1)
  log_S[1] <- 0
  sigma_sq[1] <- sigma0_sq
  
  # --- Simulate all Background Driving Lévy Process (BDLP) increments ---
  num_jumps <- rpois(n_steps, lambda * a * dt)
  bdlp_increments <- sapply(num_jumps, function(nj) {
    if (nj == 0) return(0)
    # The total increment is the sum of the individual exponential jumps
    sum(rgamma(nj, shape = 1, rate = b))
  })
  
  # Generate standard normal random variables for the main Brownian motion
  dW <- rnorm(n_steps, mean = 0, sd = sqrt(dt))
  
  # --- Euler-Maruyama Discretization using the BDLP ---
  for (t in 1:n_steps) {
    sigma_sq_prev <- max(1e-9, sigma_sq[t])
    
    # Update variance process
    # d(sigma_t^2) = -lambda * sigma_t^2 * dt + d(z_{lambda*t})
    sigma_sq[t+1] <- sigma_sq[t] - lambda * sigma_sq_prev * dt + bdlp_increments[t]
    
    # Update log-price process
    # d(log S_t) = (mu - 0.5*sigma_t^2)dt + sigma_t*dW_t + rho*d(z_{lambda*t})
    log_S[t+1] <- log_S[t] + (mu - 0.5 * sigma_sq_prev) * dt + 
      sqrt(sigma_sq_prev) * dW[t] + 
      rho * bdlp_increments[t]
  }
  
  # Return the simulated log-returns
  return(diff(log_S))
}