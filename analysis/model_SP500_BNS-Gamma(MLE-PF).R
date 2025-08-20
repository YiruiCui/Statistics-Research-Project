library(here)
library(tidyverse)
library(ghyp)
library(gridExtra)

set.seed(7914) # For reproducibility

# Identify project location
here::i_am("analysis/model_SP500_BNS-Gamma(MLE-PF).R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))

# --- BNS Simulation Function ---
simulate_bns_gamma_sv <- function(params, n_steps, dt) {
  # Unpack parameters from the input vector
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3]
  a         <- params[4]
  b         <- params[5]
  sigma0_sq <- params[6]
  
  log_S <- numeric(n_steps + 1)
  sigma_sq <- numeric(n_steps + 1)
  log_S[1] <- 0
  sigma_sq[1] <- sigma0_sq
  
  num_jumps <- rpois(n_steps, lambda * a * dt)
  bdlp_increments <- sapply(num_jumps, function(nj) {
    if (nj == 0) return(0)
    # The total increment is the sum of the individual exponential jumps
    sum(rgamma(nj, shape = 1, rate = b))
  })
  
  dW <- rnorm(n_steps, mean = 0, sd = sqrt(dt))
  
  # --- Euler-Maruyama Discretization using the correct BDLP ---
  for (t in 1:n_steps) {
    sigma_sq_prev <- max(1e-9, sigma_sq[t])
    
    # Update variance process (dσ² = -λσ²dt + dz)
    sigma_sq[t+1] <- sigma_sq[t] - lambda * sigma_sq_prev * dt + bdlp_increments[t]
    
    # Update log-price process (dZ = (μ - σ²/2)dt + σdW + ρdz) [cite: 2491]
    log_S[t+1] <- log_S[t] + (mu - 0.5 * sigma_sq_prev) * dt + 
      sqrt(sigma_sq_prev) * dW[t] + 
      rho * bdlp_increments[t]
  }
  
  return(diff(log_S))
}

# --- 2. Particle Filter Log-Likelihood Function ---
# This function approximates the log-likelihood of the data for a given
# set of BNS parameters.

particle_filter_log_likelihood <- function(params, data, n_particles) {
  # Unpack parameters
  mu        <- params[1]
  rho       <- params[2]
  lambda    <- params[3]
  a         <- params[4]
  b         <- params[5]
  sigma0_sq <- params[6]
  
  n_obs <- length(data)
  log_likelihood <- 0
  dt <- 1 # Daily data
  
  # --- Initialization ---
  particles_sigma_sq <- rep(sigma0_sq, n_particles)
  
  # --- Main Filtering Loop ---
  for (t in 1:n_obs) {
    # --- Step 1: Prediction/Propagation ---
    num_jumps <- rpois(n_particles, lambda * a * dt)
    bdlp_increments <- sapply(num_jumps, function(nj) {
      if (nj == 0) return(0)
      # The increment is the sum of Gamma(1,b) i.e. Exponential(b) jumps
      sum(rgamma(nj, shape = 1, rate = b))
    })
    
    # Propagate each particle's state (variance) forward using its unique BDLP increment
    # dσ² = -λσ²dt + dz
    particles_sigma_sq <- pmax(1e-9, particles_sigma_sq * (1 - lambda * dt) + bdlp_increments)
    
    # --- Step 2: Weighting ---
    
    # The weight is the probability of the observed return data[t], conditional
    # on the particle's variance AND the simulated jump (bdlp_increments).
    # We use the observation equation:
    # Return = (μ - σ²/2)dt + σdW + ρdz
    # So, the Gaussian part is: Return - (μ - σ²/2)dt - ρdz = σdW
    
    # This is the mean of the Gaussian component of the return
    particle_means <- (mu - 0.5 * particles_sigma_sq) * dt + rho * bdlp_increments
    particle_sds <- sqrt(particles_sigma_sq * dt)
    
    # The weight is the density of the observed return under this conditional distribution
    log_weights <- dnorm(data[t], mean = particle_means, sd = particle_sds, log = TRUE)
    
    # Avoid numerical underflow by shifting log-weights
    max_log_weight <- max(log_weights)
    weights <- exp(log_weights - max_log_weight)
    
    # --- Step 3: Update Log-Likelihood & Resample ---
    if (sum(weights) == 0 || !is.finite(sum(weights))) {
      # If all weights are zero, the filter has failed. Return a large penalty.
      return(1e9) 
    }
    
    log_likelihood <- log_likelihood + max_log_weight + log(mean(weights))
    
    # Resample particles based on their weights
    normalized_weights <- weights / sum(weights)
    indices <- sample(1:n_particles, size = n_particles, replace = TRUE, prob = normalized_weights)
    particles_sigma_sq <- particles_sigma_sq[indices]
  }
  
  return(-log_likelihood) # Return negative log-likelihood for minimization
}


# --- 3. Define the MLE Objective Function ---
mle_objective_function <- function(scaled_params, data, n_particles) {
  
  # Unscale parameters to their natural domain
  params <- numeric(6)
  params[1] <- scaled_params[1]               # mu
  params[2] <- -exp(scaled_params[2])         # rho <=0
  params[3] <- exp(scaled_params[3])          # lambda > 0
  params[4] <- exp(scaled_params[4])          # a > 0
  params[5] <- exp(scaled_params[5])          # b > 0
  params[6] <- exp(scaled_params[6])          # sigma0_sq > 0
  
  # Calculate the negative log-likelihood
  neg_log_lik <- particle_filter_log_likelihood(params, data, n_particles)
  
  # Print progress
  cat("Testing Params:", round(params, 4), " | -LogLik:", round(neg_log_lik, 4), "\n")
  
  if (!is.finite(neg_log_lik)) {
    return(1e9) # Return a large penalty if likelihood is not finite
  }
  
  return(neg_log_lik)
}

# --- 4. Run the Optimization ---
# Starting values on the transformed scale
initial_scaled_params <- c(
  mu = mean(log_returns),
  rho_trans = log(1.0), 
  log_lambda = log(0.05),
  log_a = log(1.0),
  log_b = log(8000),
  log_sigma0_sq = log(var(log_returns))
)

# NOTE: Particle filtering is computationally intensive.
# For a real estimation, n_particles should be higher (e.g., 5000+) and
# the optimization will take a very long time.
n_particles <- 10000 

cat("\nStarting MLE optimization via Particle Filter (this will be very slow)...\n")
mle_results <- optim(
  par = initial_scaled_params,
  fn = mle_objective_function,
  # Additional arguments for the objective function
  data = log_returns,
  n_particles = n_particles,
  method = "Nelder-Mead",
  control = list(maxit = 200, trace = 1)
)

# --- 5. Display Results ---
cat("\n--- MLE Optimization Finished ---\n")
print("Optimal Scaled Parameters Found:")
print(mle_results$par)

# Transform parameters back to their original scale
estimated_params_mle <- numeric(6)
estimated_params_mle[1] <- mle_results$par[1]
estimated_params_mle[2] <- -exp(mle_results$par[2])
estimated_params_mle[3] <- exp(mle_results$par[3])
estimated_params_mle[4] <- exp(mle_results$par[4])
estimated_params_mle[5] <- exp(mle_results$par[5])
estimated_params_mle[6] <- exp(mle_results$par[6])
names(estimated_params_mle) <- c("mu", "rho", "lambda", "a", "b", "sigma0_sq")

print("Optimal Interpretable Parameters (MLE):")
print(estimated_params_mle)

cat("\nFinal Minimized Negative Log-Likelihood:", mle_results$value, "\n")

cat("\nSimulating final BNS model path with estimated parameters...\n")

bns_returns <- simulate_bns_gamma_sv(
  params = estimated_params_mle,
  n_steps = length(log_returns),
  dt = 1
)
cat("BNS simulation complete.\n")


# --- Comparison Plots ---
cat("Generating comparison plots...\n")

# A. Density Overlay Plot
x_grid <- seq(min(log_returns), max(log_returns), length.out = 500)
gh_density_data <- data.frame(x = x_grid, y = dghyp(x_grid, object = gh_fit))

density_plot_hist_style <- ggplot(data.frame(value = log_returns), aes(x = value)) +
  # Empirical data as a histogram.
  geom_histogram(aes(y = after_stat(density), fill = "S&P 500 Empirical"), 
                 binwidth = 0.001, 
                 alpha = 0.7, 
                 color = "white") +
  # BNS simulated density as a line
  geom_density(data = data.frame(value = bns_returns), 
               aes(color = "BNS (Simulated)"), linewidth = 1.2) +
  # GH fitted density as a line
  geom_line(data = gh_density_data, 
            aes(x = x, y = y, color = "GH (Fitted)"), 
            linewidth = 1.2) +
  # Labels and Titles
  labs(
    title = "Static Comparison: Unconditional Return Distributions", 
    x = "Log Return", 
    y = "Density"
  ) +
  # Colors and Theme
  scale_fill_manual(values = c("S&P 500 Empirical" = "lightblue")) +
  scale_color_manual(values = c("BNS (Simulated)" = "blue", "GH (Fitted)" = "red")) +
  theme_minimal(base_size = 14) +
  coord_cartesian(xlim = quantile(log_returns, c(0.005, 0.995))) +
  guides(fill = guide_legend(title = NULL), color = guide_legend(title = "Model Overlays"))

print(density_plot_hist_style)

# B. Q-Q Plot Comparison
# Generate random sample from the fitted GH model
gh_samples <- rghyp(length(log_returns), object = gh_fit)

# Prepare data for Q-Q plots
qq_data <- data.frame(
  empirical = sort(log_returns), 
  gh = sort(gh_samples), 
  bns = sort(bns_returns))

qq_plot_gh <- ggplot(qq_data, aes(x = gh, y = empirical)) +
  geom_point(alpha = 0.5, color = "red") +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  labs(
    title = "Q-Q Plot: Empirical vs. Fitted GH", 
    x = "Theoretical Quantiles (GH)", 
    y = "Empirical Quantiles (S&P 500)"
  ) +
  theme_minimal(base_size = 14)

qq_plot_bns <- ggplot(qq_data, aes(x = bns, y = empirical)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  labs(
    title = "Q-Q Plot: Empirical vs. Simulated BNS", 
    x = "Theoretical Quantiles (BNS)", 
    y = "Empirical Quantiles (S&P 500)"
  ) +
  theme_minimal(base_size = 14)

# Arrange plots side-by-side
QQplots <- grid.arrange(qq_plot_gh, qq_plot_bns, ncol = 2)
print(QQplots)

# C. Dynamic Comparison (Volatility Clustering)
# Calculate ACF for Absolute returns for all three series
acf_empirical <- acf(abs(log_returns), plot = FALSE, lag.max = 50)
acf_bns <- acf(abs(bns_returns), plot = FALSE, lag.max = 50)
acf_gh <- acf(abs(gh_samples), plot = FALSE, lag.max = 50)

# Create a data frame for the ACF plot
acf_data <- data.frame(
  Lag = acf_empirical$lag,
  ACF = c(acf_empirical$acf, acf_bns$acf, acf_gh$acf),
  Model = factor(rep(c("S&P 500 Empirical", "BNS (Simulated)", "GH (Static)"), 
                     each = length(acf_empirical$lag)))
)

# Confidence interval for ACF plot (for white noise)
ci <- qnorm(0.975) / sqrt(length(log_returns))

acf_plot <- ggplot(acf_data, aes(x = Lag, y = ACF, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = c(ci, -ci), linetype = "dashed", color = "black") +
  labs(
    title = "Dynamic Comparison: Autocorrelation of Absolute Returns", 
    subtitle = "Demonstrating Volatility Clustering", 
    x = "Lag (Days)", 
    y = "Autocorrelation"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("S&P 500 Empirical" = "black", 
                               "BNS (Simulated)" = "blue", 
                               "GH (Static)" = "red")) +
  theme(legend.position = "bottom")

print(acf_plot)

## Save plots
ggsave(
  filename = here("outputs", "BNS-gamma(MLE-PF)&GH_fit.png"),
  plot = density_plot_hist_style,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-gamma(MLE-PF)&GH_QQplot.png"),
  plot = QQplots,
  width = 3000,
  height = 1500,
  units = "px",
  dpi = 300
)

ggsave(
  filename = here("outputs", "BNS-gamma(MLE-PF)&GH_acf.png"),
  plot = acf_plot,
  width = 2000,
  height = 1200,
  units = "px",
  dpi = 300
)

