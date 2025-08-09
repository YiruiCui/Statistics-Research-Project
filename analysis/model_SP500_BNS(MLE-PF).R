library(here)
library(tidyverse)
library(ghyp)

set.seed(7914) # For reproducibility

# Identify project location
here::i_am("analysis/model_SP500_BNS(MLE-PF).R")

# --- Load daily price data ---
data <- readr::read_csv(here("data", "SP500.csv"))
price_vector <- data$Last_Price[18080:nrow(data)]

# --- Compute log-returns ---
log_returns <- diff(log(price_vector))


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
  
  # --- Step 1: Initialization ---
  # Start with a cloud of 'particles' representing possible initial variance values.
  # We draw from the stationary distribution of the Gamma-OU process if possible,
  # or simply start them all at the initial guess.
  particles_sigma_sq <- rep(sigma0_sq, n_particles)
  
  # --- Step 2: Main Loop (Filtering) ---
  for (t in 1:n_obs) {
    # --- Step 2a: Prediction/Propagation ---
    # Move each particle forward one step in time according to the state equation
    # (the Gamma-OU process for variance).
    jump_intensity <- a * lambda * 1 # dt = 1 for daily
    jumps <- rpois(n_particles, jump_intensity)
    jump_sizes <- sapply(jumps, function(nj) {
      if (nj == 0) return(0)
      sum(rgamma(nj, shape = 1, rate = b))
    })
    
    particles_sigma_sq <- pmax(1e-9, particles_sigma_sq * (1 - lambda) + jump_sizes)
    
    # --- Step 2b: Weighting ---
    # Calculate the likelihood of observing the actual return data[t] for each particle.
    # This is the "observation density".
    # Note: This step can be improved by handling the leverage term more explicitly.
    # For simplicity here, we use the standard BNS return equation.
    particle_means <- mu - 0.5 * particles_sigma_sq
    particle_sds <- sqrt(particles_sigma_sq)
    
    weights <- dnorm(data[t], mean = particle_means, sd = particle_sds, log = TRUE)
    
    # Avoid numerical underflow by shifting log-weights
    max_log_weight <- max(weights)
    weights <- exp(weights - max_log_weight)
    
    # --- Step 2c: Update Log-Likelihood ---
    # The likelihood for this time step is the average of the particle weights.
    log_likelihood <- log_likelihood + max_log_weight + log(mean(weights))
    
    # --- Step 2d: Resampling ---
    # Resample the particles based on their weights. Particles with higher
    # weights (i.e., those that better explain the data) are more likely to be chosen.
    # This prevents particle degeneracy.
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
  params[2] <- scaled_params[2]               # rho in (-1, 1)
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
  rho_trans = -1.0, 
  log_lambda = log(0.05),
  log_a = log(1.0),
  log_b = log(8000),
  log_sigma0_sq = log(var(log_returns))
)

# NOTE: Particle filtering is computationally intensive.
# For a real estimation, n_particles should be higher (e.g., 5000+) and
# the optimization will take a very long time.
n_particles <- 1000 

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
estimated_params_mle[2] <- mle_results$par[2]
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

